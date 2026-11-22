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

/// The ring of a circle of radius ten as a periodic rational quadratic:
/// four quarter arcs, each opening on a pole on the circle, `w` the weight
/// of the corner poles between them. Each pole is `x y` and its weight;
/// `z` is put after `x y` when given.
fn circle_ring(z: Option<&str>) -> String {
    let w = "0.70710678118654757";
    [
        ("10", "0", "1"),
        ("10", "10", w),
        ("0", "10", "1"),
        ("-10", "10", w),
        ("-10", "0", "1"),
        ("-10", "-10", w),
        ("0", "-10", "1"),
        ("10", "-10", w),
    ]
    .iter()
    .map(|(x, y, weight)| match z {
        Some(z) => format!("{x} {y} {z} {weight}"),
        None => format!("{x} {y} {weight}"),
    })
    .collect::<Vec<_>>()
    .join(" ")
}

#[test]
fn a_periodic_spline_reads_as_the_period_it_states() {
    // A disc bounded by one closed edge, its curve and its chart the same
    // periodic circle: eight poles over knots 0 to 4, each knot doubled,
    // the multiplicities but the last's summing to the poles.
    let knots = "0 2 1 2 2 2 3 2 4 2";
    let text = format!(
        "\
DBRep_DrawableShape

CASCADE Topology V3, (c) Open Cascade
Locations 0
Curve2ds 1
7 1 1 2 8 5 {chart} {knots}
Curves 1
7 1 1 2 8 5 {curve} {knots}
Polygon3D 0
PolygonOnTriangulations 0
Surfaces 1
1 0 0 0 0 0 1 1 0 0 0 1 0
Triangulations 0

TShapes 4
Ve
1e-07
10 0 0
0 0

0101101
*
Ed
 1e-07 1 1 0
1 1 0 0 4
2 1 1 0 0 4
0

0101000
+4 0 -4 0 *
Wi

0101100
+3 0 *
Fa
0 1e-07 1 0

0111000
+2 0 *

+1 0
",
        chart = circle_ring(None),
        curve = circle_ring(Some("0")),
    );
    let (model, face) = ogeom_io::brep::read(&text, T).unwrap();
    let diagnosis = ogeom_algo::check(&model, &face, T).unwrap();
    assert!(diagnosis.is_valid(), "{:?}", diagnosis.problems);
    let area = ogeom_algo::surface_properties(
        &model,
        &face,
        ogeom_mesh::Deflection::with_chord(1e-4).unwrap(),
        T,
    )
    .unwrap()
    .mass;
    let disc = core::f64::consts::PI * 100.0;
    assert!((area - disc).abs() < 1e-2, "a disc of radius ten: {area}");

    // Every point of the curve read is on the circle, the period's ends
    // and its middle on the poles the record puts there.
    let (_, curve) = model.geometry().curves().next().unwrap();
    for i in 0..=40 {
        let u = 0.1 * f64::from(i);
        let p = curve.point_at(u, T).unwrap();
        assert!((p.x.hypot(p.y) - 10.0).abs() < 1e-9, "u = {u}: {p:?}");
    }
    assert!(curve.point_at(0.0, T).unwrap().is_equal(Point::new(10.0, 0.0, 0.0), T));
    assert!(curve.point_at(2.0, T).unwrap().is_equal(Point::new(-10.0, 0.0, 0.0), T));
    assert!(curve.point_at(4.0, T).unwrap().is_equal(Point::new(10.0, 0.0, 0.0), T));
}

#[test]
fn a_surface_periodic_one_way_reads_as_its_period() {
    // A cylinder of radius ten and height three, periodic round its `u`:
    // the same ring as the circle at each of two heights, the poles running
    // `v` fastest.
    let low = circle_ring(Some("0"));
    let high = circle_ring(Some("3"));
    let low: Vec<&str> = low.split(' ').collect();
    let high: Vec<&str> = high.split(' ').collect();
    let mut net = Vec::new();
    for (a, b) in low.chunks(4).zip(high.chunks(4)) {
        net.push(a.join(" "));
        net.push(b.join(" "));
    }
    let text = format!(
        "\
DBRep_DrawableShape

CASCADE Topology V3, (c) Open Cascade
Locations 0
Curve2ds 0
Curves 0
Polygon3D 0
PolygonOnTriangulations 0
Surfaces 1
9 1 1 1 0 2 1 8 2 5 2 {net} 0 2 1 2 2 2 3 2 4 2 0 2 1 2
Triangulations 0

TShapes 1
Ve
1e-07
10 0 0
0 0

0101101
*

+1 0
",
        net = net.join(" "),
    );
    let (model, _) = ogeom_io::brep::read(&text, T).unwrap();
    let (_, surface) = model.geometry().surfaces().next().unwrap();
    use ogeom_geom::Surface as _;
    assert_eq!(surface.domain(), ((0.0, 4.0), (0.0, 1.0)));
    for i in 0..=40 {
        for j in 0..=4 {
            let (u, v) = (0.1 * f64::from(i), 0.25 * f64::from(j));
            let p = surface.point_at(u, v, T).unwrap();
            assert!((p.x.hypot(p.y) - 10.0).abs() < 1e-9, "({u}, {v}): {p:?}");
            assert!((p.z - 3.0 * v).abs() < 1e-9, "({u}, {v}): {p:?}");
        }
    }
    let p = surface.point_at(1.0, 0.0, T).unwrap();
    assert!(p.is_equal(Point::new(0.0, 10.0, 0.0), T), "{p:?}");
}

#[test]
fn a_uniform_periodic_spline_reads_as_its_ring() {
    // Simple knots: the basis reaches a full degree past each end, so the
    // ring's first poles come round again after its last. The kernel's own
    // periodic spline through the same ring is the same curve, its domain
    // opening at the degree's knot rather than at zero.
    let ring = [
        Point::new(0.0, 0.0, 0.0),
        Point::new(4.0, 1.0, 0.0),
        Point::new(5.0, 5.0, 1.0),
        Point::new(1.0, 6.0, 2.0),
        Point::new(-2.0, 3.0, 1.0),
    ];
    let poles: Vec<String> = ring.iter().map(|p| format!("{} {} {}", p.x, p.y, p.z)).collect();
    let text = format!(
        "\
DBRep_DrawableShape

CASCADE Topology V3, (c) Open Cascade
Locations 0
Curve2ds 0
Curves 1
7 0 1 3 5 6 {} 0 1 1 1 2 1 3 1 4 1 5 1
Polygon3D 0
PolygonOnTriangulations 0
Surfaces 0
Triangulations 0

TShapes 1
Ve
1e-07
0 0 0
0 0

0101101
*

+1 0
",
        poles.join(" ")
    );
    let (model, _) = ogeom_io::brep::read(&text, T).unwrap();
    let (_, curve) = model.geometry().curves().next().unwrap();
    let own = ogeom_geom::BSplineCurve::periodic(&ring, 3, T).unwrap();
    for i in 0..=50 {
        let u = 0.1 * f64::from(i);
        let (read, expected) = (curve.point_at(u, T).unwrap(), own.point_at(u + 3.0, T).unwrap());
        assert!(read.distance(expected) < 1e-12, "u = {u}: {read:?} against {expected:?}");
    }
}
