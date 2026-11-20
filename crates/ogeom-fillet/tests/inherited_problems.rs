//! A body that comes in failing its check, broken far from the edges asked:
//! the rounds are judged on what they did, the input's own problems left as
//! they came.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]
use ogeom_core::Tolerances;
use ogeom_geom::{Curve3d as _, Offset2d, PlanarCurve, PlaneSurface};
use ogeom_math::{Direction, Frame, Plane, Point, Vector};
use ogeom_topo::{EdgeRepr, Model, NodeData, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The edge of `solid` whose middle stands at `at`.
fn edge_at(model: &Model, solid: &Shape, at: Point) -> Shape {
    explore_unique(model, solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let (c, (a, b)) = ogeom_algo::placed_edge_curve(model, e, T).unwrap().unwrap();
            c.point_at(0.5 * (a + b), T).unwrap().distance(at) < 1e-9
        })
        .expect("no edge there")
}

/// A Y standing in `x, z`, one long in `y`: a stem whose walls meet the
/// arms' undersides in concave creases of 135 degrees, and a square notch
/// between the arms. The three creases are rounded at one radius, as long
/// as they are; the stem's foot, nine below, has the pcurve of one of its
/// edges moved a hundredth off in its face's chart, which the check reads
/// as four problems the rounds never come near.
#[test]
fn a_round_far_from_a_broken_face_keeps_the_input_s_problems_and_no_more() {
    let mut model = Model::new();
    let profile: Vec<Point> = [
        (-1.0, 0.3),
        (-3.9, 3.2),
        (-3.2, 3.9),
        (0.0, 0.7),
        (3.2, 3.9),
        (3.9, 3.2),
        (1.0, 0.3),
        (1.0, -8.0),
        (0.0, -8.9),
        (-1.0, -8.0),
    ]
    .iter()
    .map(|&(x, z)| Point::new(x, 0.0, z))
    .collect();
    let wire = ogeom_algo::make_polygon(&mut model, &profile, true, T)
        .unwrap()
        .shape;
    let frame = Frame::from_axes(
        Point::ORIGIN,
        Direction::X,
        Direction::Z,
        Direction::new(Vector::new(0.0, -1.0, 0.0), T).unwrap(),
        T,
    )
    .unwrap();
    let face = ogeom_algo::make_face(
        &mut model,
        PlaneSurface::new(Plane::new(frame)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let prism = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 1.0, 0.0), T)
        .unwrap()
        .shape;
    // Every node its own, as an imported body's are: a prism's far end is
    // its profile placed again.
    let solid = ogeom_algo::baked_where_placed(&mut model, &prism, T)
        .unwrap()
        .shape;
    let creases: Vec<Shape> = [(-1.0, 0.3), (0.0, 0.7), (1.0, 0.3)]
        .iter()
        .map(|&(x, z)| edge_at(&model, &solid, Point::new(x, 0.5, z)))
        .collect();

    let foot = edge_at(&model, &solid, Point::new(-0.5, 0.0, -8.45));
    let NodeData::Edge(data) = model.node(&foot).unwrap().data().clone() else {
        panic!("not an edge")
    };
    let original = data
        .representations
        .iter()
        .find_map(|r| match r {
            EdgeRepr::PCurve { curve, .. } => Some(*curve),
            _ => None,
        })
        .unwrap();
    let basis = model.geometry().pcurve(original).unwrap().clone();
    let moved = model
        .geometry_mut()
        .add_pcurve(PlanarCurve::Offset(Box::new(
            Offset2d::new(basis, 0.01).unwrap(),
        )));
    if let Some(NodeData::Edge(data)) = model.node_mut(&foot).map(ogeom_topo::TShape::data_mut) {
        for r in &mut data.representations {
            if let EdgeRepr::PCurve { curve, .. } = r
                && *curve == original
            {
                *curve = moved;
            }
        }
    }
    let before: Vec<String> = ogeom_algo::check(&model, &solid, T)
        .unwrap()
        .problems
        .into_iter()
        .map(|p| p.what)
        .collect();
    assert_eq!(before.len(), 4, "{before:?}");

    let rounded = ogeom_fillet::fillet_edges(&mut model, &solid, &creases, 1.0, T)
        .unwrap()
        .shape;
    let mut after: Vec<String> = ogeom_algo::check(&model, &rounded, T)
        .unwrap()
        .problems
        .into_iter()
        .map(|p| p.what)
        .collect();
    let mut want = before;
    after.sort();
    want.sort();
    assert_eq!(after, want);
    // The three bands are there: a cylinder each.
    let cylinders = explore_unique(&model, &rounded, ShapeType::Face)
        .unwrap()
        .into_iter()
        .filter(|f| {
            let NodeData::Face(d) = model.node(f).unwrap().data() else {
                return false;
            };
            matches!(
                model.geometry().surface(d.surface),
                Some(ogeom_geom::SurfaceGeometry::Cylinder(_))
            )
        })
        .count();
    assert_eq!(cylinders, 3);
}
