//! Records the `.brep` reader takes as the format's specification spells
//! them, each written by hand into an otherwise plain file: a ten by ten
//! square face whose first edge carries the record under test.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_geom::Curve3d as _;
use ogeom_math::Point;
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// The square, with its first edge's chart curve, curve and parameter range
/// (lowest, highest) given.
fn square(chart: &str, curve: &str, range: (&str, &str)) -> String {
    let (lo, hi) = range;
    format!(
        "\
DBRep_DrawableShape

CASCADE Topology V3, (c) Open Cascade
Locations 0
Curve2ds 4
{chart}
1 10 0 0 1
1 10 10 -1 0
1 0 10 0 -1
Curves 4
{curve}
1 10 0 0 0 1 0
1 10 10 0 -1 0 0
1 0 10 0 0 -1 0
Polygon3D 0
PolygonOnTriangulations 0
Surfaces 1
1 0 0 0 0 0 1 1 0 0 0 1 0
Triangulations 0

TShapes 10
Ve
1e-07
0 0 0
0 0

0101101
*
Ve
1e-07
10 0 0
0 0

0101101
*
Ed
 1e-07 1 1 0
1 1 0 {lo} {hi}
2 1 1 0 {lo} {hi}
0

0101000
+10 0 -9 0 *
Ve
1e-07
10 10 0
0 0

0101101
*
Ed
 1e-07 1 1 0
1 2 0 0 10
2 2 1 0 0 10
0

0101000
+9 0 -7 0 *
Ve
1e-07
0 10 0
0 0

0101101
*
Ed
 1e-07 1 1 0
1 3 0 0 10
2 3 1 0 0 10
0

0101000
+7 0 -5 0 *
Ed
 1e-07 1 1 0
1 4 0 0 10
2 4 1 0 0 10
0

0101000
+5 0 -10 0 *
Wi

0101100
+8 0 +6 0 +4 0 +3 0 *
Fa
0 1e-07 1 0

0111000
+2 0 *

+1 0
"
    )
}

/// The face read, checked to be the square: one face, four edges, an area
/// of a hundred.
fn read_square(text: &str) -> (Model, Shape) {
    let (model, face) = ogeom_io::brep::read(text, T).unwrap();
    assert_eq!(model.kind_of(&face).unwrap(), ShapeType::Face);
    let edges = explore_unique(&model, &face, ShapeType::Edge).unwrap();
    assert_eq!(edges.len(), 4);
    let area = ogeom_algo::surface_properties(
        &model,
        &face,
        ogeom_mesh::Deflection::default(),
        T,
    )
    .unwrap()
    .mass;
    assert!((area - 100.0).abs() < 1e-9, "a ten by ten square: {area}");
    (model, face)
}

#[test]
fn bezier_curves_read_as_the_splines_they_are() {
    // A rational quadratic chart curve with equal weights and a linear one
    // in space, both running (0, 0) to (10, 0) over [0, 1]: record 6 in both
    // tables, the weights following each pole.
    let text = square("6 1 2 0 0 2 5 0 2 10 0 2", "6 0 1 0 0 0 10 0 0", ("0", "1"));
    let (model, face) = read_square(&text);
    let mut found = false;
    for edge in explore_unique(&model, &face, ShapeType::Edge).unwrap() {
        let Some((curve, (a, b))) = ogeom_algo::placed_edge_curve(&model, &edge, T).unwrap() else {
            continue;
        };
        let mid = curve.point_at(0.5 * (a + b), T).unwrap();
        if mid.is_equal(Point::new(5.0, 0.0, 0.0), T) {
            assert!((a, b) == (0.0, 1.0), "the Bezier edge runs [{a}, {b}]");
            found = true;
        }
    }
    assert!(found, "the Bezier edge is in the square");

    // The degree is bounded as the format bounds it, so a hostile count
    // cannot size the knot vector.
    let text = square("6 0 26 0 0 10 0", "1 0 0 0 1 0 0", ("0", "10"));
    assert!(ogeom_io::brep::read(&text, T).is_err());
}
