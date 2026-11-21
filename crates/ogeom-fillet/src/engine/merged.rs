//! Faces a blend leaves in pieces on one carrier, merged after it.
//!
//! A route that applies its strips by a boolean, or closes a corner over
//! a cell, can leave two faces sharing an edge and one carrier facing one
//! way: one surface standing as two faces with a line between them that no
//! feature of the part draws. Mainstream CAD shows no such line after a
//! blend. The pieces the blend made are merged with their neighbours on
//! one carrier ([`ogeom_heal::unify_same_domain_around`]), and the merge's
//! history follows the blend's, so every input element is told of the face
//! that now holds it. Splits the input already held between faces the
//! blend left alone stay as they were, whether the route kept those faces'
//! nodes or rebuilt them where they stood: the trimmed assembly rebuilds
//! every face of the solid, and a face it put back on its carrier within
//! its old boundary is no face the blend made.

use ogeom_algo::Built;
use ogeom_core::{FastSet, OgeomError, OgeomResult, Tolerances};
use ogeom_math::Point;
use ogeom_topo::{Model, Shape, ShapeType, ancestors_of, explore_unique};

/// `built`, a blend of `input`, with the faces it made merged with their
/// neighbours on one carrier where any two such share an edge. `built` as
/// it is where none do, or where the merge does not apply (a cluster whose
/// boundary does not re-chain is left split by the merge itself).
///
/// # Errors
///
/// [`OgeomError::Cancelled`] when cancelled, and as reading the solids.
pub(crate) fn merged(
    model: &mut Model,
    input: &Shape,
    built: Built,
    tol: Tolerances,
) -> OgeomResult<Built> {
    let given = explore_unique(model, input, ShapeType::Face)?;
    let before: FastSet<_> = given.iter().map(Shape::node).collect();
    let mut standing = Vec::with_capacity(given.len());
    for face in &given {
        standing.push(Standing::of(model, face)?);
    }
    let mut made: Vec<Shape> = Vec::new();
    for face in explore_unique(model, &built.shape, ShapeType::Face)? {
        if before.contains(&face.node())
            || rebuilt_as_given(model, &face, &given, &standing, tol)?
        {
            continue;
        }
        made.push(face);
    }
    if made.is_empty() || !any_split(model, &built.shape, &made, tol)? {
        return Ok(built);
    }
    let unified = match ogeom_heal::unify_same_domain_around(model, &built.shape, &made, tol) {
        Ok((_, 0)) => return Ok(built),
        Ok((unified, _)) => unified,
        Err(e @ OgeomError::Cancelled) => return Err(e),
        Err(_) => return Ok(built),
    };
    if !ogeom_algo::check(model, &unified.shape, tol)?.is_valid() {
        return Ok(built);
    }
    Ok(Built::new(
        unified.shape,
        built.history.then(&unified.history),
    ))
}

/// A face read where it stands: the points of its vertices and the number
/// of its edges.
struct Standing {
    points: Vec<Point>,
    edges: usize,
}

impl Standing {
    fn of(model: &Model, face: &Shape) -> OgeomResult<Self> {
        let mut points = Vec::new();
        for v in explore_unique(model, face, ShapeType::Vertex)? {
            let Some(data) = model.node(&v).and_then(|n| n.data().as_vertex()) else {
                continue;
            };
            points.push(v.transform(model.datums())?.apply(data.point));
        }
        Ok(Self {
            points,
            edges: explore_unique(model, face, ShapeType::Edge)?.len(),
        })
    }

    /// The same vertices, each within `near` of one of the other's, and as
    /// many edges.
    fn alike(&self, other: &Self, near: f64) -> bool {
        self.edges == other.edges
            && self.points.len() == other.points.len()
            && self
                .points
                .iter()
                .all(|p| other.points.iter().any(|q| q.distance(*p) <= near))
            && other
                .points
                .iter()
                .all(|p| self.points.iter().any(|q| q.distance(*p) <= near))
    }
}

/// Whether `face`, a face of the blend's result on a node the input does
/// not hold, is a face of the input put back where it stood: on its
/// carrier facing its way ([`ogeom_heal::same_domain`]), through the same
/// vertices with as many edges. A route rebuilding the whole solid (the
/// trimmed assembly) leaves every face on a new node; the faces it left
/// alone are no faces the blend made, and the splits between them are the
/// input's own. A face with no vertex (a whole sphere) is not matched.
fn rebuilt_as_given(
    model: &Model,
    face: &Shape,
    given: &[Shape],
    standing: &[Standing],
    tol: Tolerances,
) -> OgeomResult<bool> {
    let here = Standing::of(model, face)?;
    if here.points.is_empty() {
        return Ok(false);
    }
    for (old, was) in given.iter().zip(standing) {
        if here.alike(was, tol.confusion())
            && ogeom_heal::same_domain(model, face, old, tol)?
        {
            return Ok(true);
        }
    }
    Ok(false)
}

/// Whether a face of `made` shares an edge of `solid` with another face on
/// its carrier facing its way.
fn any_split(model: &Model, solid: &Shape, made: &[Shape], tol: Tolerances) -> OgeomResult<bool> {
    let made_nodes: FastSet<_> = made.iter().map(Shape::node).collect();
    for face in made {
        for edge in explore_unique(model, face, ShapeType::Edge)? {
            let faces = ancestors_of(model, solid, &edge, ShapeType::Face)?;
            let [a, b] = faces.as_slice() else {
                continue;
            };
            if a.node() == b.node() {
                continue;
            }
            if !(made_nodes.contains(&a.node()) || made_nodes.contains(&b.node())) {
                continue;
            }
            if ogeom_heal::same_domain(model, a, b, tol)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
    use super::merged;
    use ogeom_core::Tolerances;
    use ogeom_math::{Direction, Frame, Point};
    use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

    const T: Tolerances = Tolerances::millimetres();

    fn block(m: &mut Model, at: (f64, f64, f64), size: (f64, f64, f64)) -> Shape {
        let frame =
            Frame::new(Point::new(at.0, at.1, at.2), Direction::Z, Direction::X, T).unwrap();
        ogeom_algo::make_box(m, frame, size, T).unwrap().shape
    }

    fn faces(m: &Model, s: &Shape) -> usize {
        explore_unique(m, s, ShapeType::Face).unwrap().len()
    }

    /// Two blocks fused end to end and left split where they meet: their
    /// top, bottom and sides each stand as two faces on one plane. A route
    /// rebuilding every face of the solid (the trimmed assembly) hands back
    /// all of them on new nodes, each where it stood; none is a face the
    /// blend made, and the merge after it leaves the input's splits as they
    /// were. Merged, the four splits fell together: a part rounded far from
    /// its splits came back with fewer faces than it was given, and a
    /// feature after the round that named one piece of a split face read
    /// the whole.
    #[test]
    fn the_input_s_own_splits_stand_when_a_route_rebuilds_every_face() {
        let mut m = Model::new();
        let a = block(&mut m, (0.0, 0.0, 0.0), (2.0, 1.0, 1.0));
        let b = block(&mut m, (2.0, 0.0, 0.0), (2.0, 1.0, 1.0));
        let fused = ogeom_bool::fuse(&mut m, &a, &b, T).unwrap().shape;
        let solid = explore_unique(&m, &fused, ShapeType::Solid).unwrap()[0].clone();
        assert_eq!(faces(&m, &solid), 10, "the fuse leaves its splits");
        let rebuilt = ogeom_algo::copied(&mut m, &solid).unwrap();
        let after = merged(&mut m, &solid, rebuilt, T).unwrap();
        assert_eq!(faces(&m, &after.shape), 10);
        assert!(ogeom_algo::check(&m, &after.shape, T).unwrap().is_valid());
    }
}
