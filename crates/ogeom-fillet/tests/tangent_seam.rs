//! A concave band whose far contact crosses a tangent seam: one of its
//! supports is two faces meeting G1, the seam running oblique to the
//! crease, so the face the edge lies on narrows to a corner between the
//! crease and the seam and the ball, near that corner, touches the face
//! across the seam instead.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::{Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The middle of `edge`'s curve.
fn middle(model: &Model, edge: &Shape) -> Option<Point> {
    let (curve, (a, b)) = ogeom_algo::placed_edge_curve(model, edge, T).ok()??;
    curve.point_at(0.5 * (a + b), T).ok()
}

/// A 20 × 20 × 10 block whose top edge along X at the back is rounded at
/// 5 (the top's plane and the round's drum meet tangent along the line
/// y = 15, z = 10), with an upright post of radius 2.5 standing on it at
/// (10, 17), straddling that seam: most of the post's foot lies on the
/// drum, a short arc of it on the top's plane. The crease the post makes
/// with the drum is the band's edge; the drum narrows to a corner at each
/// of its ends, between the crease and the seam, and the ball rolling
/// along the crease touches the plane across the seam near those corners.
fn post_across_a_tangent_seam() -> (Model, Shape, Vec<Shape>) {
    let mut model = Model::new();
    let block = ogeom_algo::make_box(&mut model, Frame::WORLD, (20.0, 20.0, 10.0), T).unwrap();
    let back = explore_unique(&model, &block.shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            middle(&model, e).is_some_and(|m| {
                m.distance(Point::new(10.0, 20.0, 10.0)) < 1e-6
            })
        })
        .expect("the block's back top edge");
    let rounded = ogeom_fillet::fillet_edge(&mut model, &block.shape, &back, 5.0, T).unwrap();
    let foot = Frame::new(
        Point::new(10.0, 17.0, 0.0),
        ogeom_math::Direction::new(ogeom_math::Vector::Z, T).unwrap(),
        ogeom_math::Direction::new(ogeom_math::Vector::X, T).unwrap(),
        T,
    )
    .unwrap();
    let post = ogeom_algo::make_cylinder(&mut model, foot, 2.5, 15.0, T).unwrap();
    let joined = ogeom_bool::fuse(&mut model, &rounded.shape, &post.shape, T).unwrap();
    // The crease on the drum: on the post's wall, behind the seam, above
    // the drum's axis.
    let crease = explore_unique(&model, &joined.shape, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .filter(|e| {
            middle(&model, e).is_some_and(|m| {
                let off_axis = ((m.x - 10.0).powi(2) + (m.y - 17.0).powi(2)).sqrt();
                (off_axis - 2.5).abs() < 1e-3 && m.y > 15.0 + 1e-3 && m.z > 5.0 && m.z < 10.0
            })
        })
        .collect::<Vec<_>>();
    // The post's own seam splits it in two.
    assert_eq!(crease.len(), 2, "the crease on the drum, in two arcs");
    (model, joined.shape, crease)
}

/// The band would have to hand its support over from the drum to the
/// plane at the seam, a transition the march does not make: the round is
/// refused, and the refusal names that seam rather than a fill laid over
/// nothing (the ball there does touch a face of the solid).
#[test]
fn a_concave_band_whose_contact_crosses_a_tangent_seam_is_refused_by_that_name() {
    let (mut model, solid, crease) = post_across_a_tangent_seam();
    let Err(refusal) = ogeom_fillet::fillet_edges(&mut model, &solid, &crease, 1.0, T) else {
        panic!("the band's support transition is not built yet; a round here is news");
    };
    let why = refusal.to_string();
    assert!(why.contains("crosses a tangent seam"), "{why}");
    assert!(!why.contains("over nothing"), "{why}");
}
