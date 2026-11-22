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

#[test]
fn a_trim_wider_than_the_open_window_reads() {
    // A line has no window of its own; the trim around it is the file's
    // word on how far it runs, here twenty times past the window a bare
    // line is read over, in both tables.
    let text = square(
        "8 -20000000 20000000\n1 0 0 1 0",
        "8 -20000000 20000000\n1 0 0 0 1 0 0",
        ("0", "10"),
    );
    read_square(&text);

    // A curve that has a domain of its own still holds the trim to it: a
    // spline over [0, 1] trimmed to [-1, 2] is refused.
    let text = square(
        "8 -1 2\n7 0 0 1 2 2 0 0 10 0 0 2 1 2",
        "1 0 0 0 1 0 0",
        ("0", "10"),
    );
    let refused = ogeom_io::brep::read(&text, T).unwrap_err().to_string();
    assert!(refused.contains("trim range"), "{refused}");
}

#[test]
fn a_surface_trimmed_one_way_reads_its_other_way_open() {
    // The writing kernel spells a direction left untrimmed with its own
    // infinity: the plane here is trimmed in v only.
    let plane = "Surfaces 1
1 0 0 0 0 0 1 1 0 0 0 1 0";
    let square_on = |surface: &str| {
        square("1 0 0 1 0", "1 0 0 0 1 0 0", ("0", "10")).replace(plane, surface)
    };
    read_square(&square_on(
        "Surfaces 1
10 -2e+100 2e+100 -5 15
1 0 0 0 0 0 1 1 0 0 0 1 0",
    ));

    // A finite bound past the window a plane is read over is a trim, and
    // still refused.
    let refused = ogeom_io::brep::read(
        &square_on("Surfaces 1
10 -2e+100 2e+100 -5 2e9
1 0 0 0 0 0 1 1 0 0 0 1 0"),
        T,
    )
    .unwrap_err()
    .to_string();
    assert!(refused.contains("v range"), "{refused}");
}

/// A quarter of a cone's band, `u` over `[0, pi/2]` and `v` over `[0, 2]`
/// along the generator, on the surface record `surface`.
fn cone_band(surface: &str) -> String {
    format!(
        "\
DBRep_DrawableShape

CASCADE Topology V3, (c) Open Cascade
Locations 0
Curve2ds 4
1 0 0 1 0
1 1.5707963267948966 0 0 1
1 0 2 1 0
1 0 0 0 1
Curves 4
2 0 0 0 0 0 1 1 0 0 0 1 0 10
1 0 10 0 0 0.70710678118654757 0.70710678118654757
2 0 0 1.4142135623730951 0 0 1 1 0 0 0 1 0 11.414213562373096
1 10 0 0 0.70710678118654757 0 0.70710678118654757
Polygon3D 0
PolygonOnTriangulations 0
Surfaces 1
{surface}
Triangulations 0

TShapes 10
Ve
1e-07
10 0 0
0 0

0101101
*
Ve
1e-07
0 10 0
0 0

0101101
*
Ed
 1e-07 1 1 0
1 1 0 0 1.5707963267948966
2 1 1 0 0 1.5707963267948966
0

0101000
+10 0 -9 0 *
Ve
1e-07
0 11.414213562373096 1.4142135623730951
0 0

0101101
*
Ed
 1e-07 1 1 0
1 2 0 0 2
2 2 1 0 0 2
0

0101000
+9 0 -7 0 *
Ve
1e-07
11.414213562373096 0 1.4142135623730951
0 0

0101101
*
Ed
 1e-07 1 1 0
1 3 0 0 1.5707963267948966
2 3 1 0 0 1.5707963267948966
0

0101000
+5 0 -7 0 *
Ed
 1e-07 1 1 0
1 4 0 0 2
2 4 1 0 0 2
0

0101000
+10 0 -5 0 *
Wi

0101100
+8 0 +6 0 -4 0 -3 0 *
Fa
0 1e-07 1 0

0111000
+2 0 *

+1 0
"
    )
}

#[test]
fn a_trimmed_cone_reads_its_window_along_the_generator() {
    // The format measures a cone's `v` along its generator, and so does a
    // trim around one: a band trimmed to [-1, 3] of that `v`, or trimmed in
    // `u` alone, is the same band as on the bare cone.
    let cone = "3 0 0 0 0 0 1 1 0 0 0 1 0 10 0.78539816339744828";
    // (pi / 2) times the mean radius 10 + sin(pi / 4), times the width 2.
    let expected = core::f64::consts::FRAC_PI_2 * (20.0 + 2.0_f64.sqrt());
    for surface in [
        cone.to_string(),
        format!("10 0 6.2831853071795862 -1 3\n{cone}"),
        format!("10 0 6.2831853071795862 -2e+100 2e+100\n{cone}"),
    ] {
        let (model, face) = ogeom_io::brep::read(&cone_band(&surface), T).unwrap();
        let diagnosis = ogeom_algo::check(&model, &face, T).unwrap();
        assert!(diagnosis.is_valid(), "{surface}: {:?}", diagnosis.problems);
        let area = ogeom_algo::surface_properties(
            &model,
            &face,
            ogeom_mesh::Deflection::with_chord(1e-4).unwrap(),
            T,
        )
        .unwrap()
        .mass;
        assert!((area - expected).abs() < 1e-3, "{surface}: {area} against {expected}");
    }
}
