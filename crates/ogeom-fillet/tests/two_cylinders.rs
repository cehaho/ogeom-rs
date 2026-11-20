//! Concave rounds between two walls along their rulings: three
//! cylinders fused side by side, the two concave seams on one side rounded,
//! each strip the exact cylinder whose axis stands a radius off both walls.
#![allow(clippy::unwrap_used, clippy::expect_used, reason = "test code")]

use ogeom_core::Tolerances;
use ogeom_math::{Direction, Frame, Point};
use ogeom_topo::{Model, Shape, ShapeType, explore_unique};

const T: Tolerances = Tolerances::millimetres();

/// Three cylinders of radius `wall`, `height` tall, their axes `pitch` apart
/// along `x`, fused, and the two concave seams on the side `y < 0`.
fn lobes(model: &mut Model, wall: f64, pitch: f64, height: f64) -> (Shape, Vec<Shape>) {
    let mut lobes = Vec::new();
    for i in 0..3 {
        let frame = Frame::new(
            Point::new(f64::from(i) * pitch, 0.0, 0.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        lobes.push(
            ogeom_algo::make_cylinder(model, frame, wall, height, T)
                .unwrap()
                .shape,
        );
    }
    let two = ogeom_bool::fuse(model, &lobes[0], &lobes[1], T)
        .unwrap()
        .shape;
    let solid = ogeom_bool::fuse(model, &two, &lobes[2], T).unwrap().shape;
    let y = -(wall * wall - pitch * pitch / 4.0).sqrt();
    let point = |model: &Model, v: &Shape| {
        v.transform(model.datums())
            .unwrap()
            .apply(model.node(v).unwrap().data().as_vertex().unwrap().point)
    };
    let seams = [0.5 * pitch, 1.5 * pitch]
        .iter()
        .map(|&x| {
            explore_unique(model, &solid, ShapeType::Edge)
                .unwrap()
                .into_iter()
                .find(|e| {
                    ogeom_algo::edge_vertices(model, e)
                        .unwrap()
                        .is_some_and(|(a, b)| {
                            let (a, b) = (point(model, &a), point(model, &b));
                            [a, b]
                                .iter()
                                .all(|p| (p.x - x).abs() < 1e-7 && (p.y - y).abs() < 1e-7)
                                && (a.z - b.z).abs() > 0.5 * height
                        })
                })
                .expect("a concave seam")
        })
        .collect();
    (solid, seams)
}

/// The rounded section's area: three discs, less the two lenses where
/// neighbours overlap, plus each round's fill, the region between the two
/// walls and the ball's arc. A fill is twice its half from the seam's
/// line `x = pitch / 2` out to the contact, under the chord of the wall
/// and above the ball's circle, integrated in closed form.
fn section_area(wall: f64, pitch: f64, radius: f64) -> f64 {
    // The area under `sqrt(r^2 - x^2)` from nought to `x`.
    let under = |r: f64, x: f64| 0.5 * (x * (r * r - x * x).sqrt() + r * r * (x / r).asin());
    let half = 0.5 * pitch;
    let reach = wall + radius;
    // The ball's centre below the seam, and the contact's offset from the
    // wall's axis along `x`.
    let depth = (reach * reach - half * half).sqrt();
    let contact = wall * half / reach;
    let lens =
        2.0 * wall * wall * (half / wall).acos() - pitch * (wall * wall - half * half).sqrt();
    // Between `x = contact` and `x = half` (from the wall's axis): below the
    // wall's lower arc and above the ball's upper arc about `(half,
    // -depth)`, each measured down from `y = 0` as the area over it.
    let below_wall = under(wall, half) - under(wall, contact);
    let below_ball =
        depth * (half - contact) - (under(radius, half - contact) - under(radius, 0.0));
    let fill = 2.0 * (below_ball - below_wall);
    3.0 * core::f64::consts::PI * wall * wall - 2.0 * lens + 2.0 * fill
}

/// Two cases: walls 3.58 apart 5.5 rounded at 8, and 4.85 apart 8.5
/// at 4. Built by the engine with no boolean, valid, the same solid either
/// way round, its volume the section's area times the height within 1e-6,
/// and each round the cylinder about the ball's exact centre line.
#[test]
fn concave_rounds_between_fused_cylinders_are_exact_cylinders() {
    for (wall, pitch, height, radius) in [(3.58_f64, 5.5_f64, 14.0, 8.0), (4.85, 8.5, 12.0, 4.0)] {
        let mut volumes = Vec::new();
        for flip in [false, true] {
            let mut model = Model::new();
            let (solid, mut seams) = lobes(&mut model, wall, pitch, height);
            if flip {
                seams.reverse();
            }
            let (built, work) = ogeom_bool::work::counted(|| {
                ogeom_fillet::fillet_edges(&mut model, &solid, &seams, radius, T).unwrap()
            });
            assert_eq!(work.face_pairs, 0, "built with no boolean");
            let check = ogeom_algo::check(&model, &built.shape, T).unwrap();
            assert!(check.is_valid(), "{:?}", check.problems);
            volumes.push(
                ogeom_algo::volume_properties(
                    &model,
                    &built.shape,
                    ogeom_mesh::Deflection::default(),
                    T,
                )
                .unwrap()
                .mass,
            );
            let depth = ((wall + radius).powi(2) - 0.25 * pitch * pitch).sqrt();
            for x in [0.5 * pitch, 1.5 * pitch] {
                let centre = Point::new(x, -depth, 0.5 * height);
                let found = explore_unique(&model, &built.shape, ShapeType::Face)
                    .unwrap()
                    .into_iter()
                    .filter(|f| {
                        let data = model.node(f).unwrap().data().as_face().unwrap();
                        let Some(ogeom_geom::SurfaceGeometry::Cylinder(c)) =
                            model.geometry().surface(data.surface)
                        else {
                            return false;
                        };
                        let axis = c.cylinder().axis();
                        let placed = f.transform(model.datums()).unwrap();
                        let (o, d) = (
                            placed.apply(axis.location),
                            placed.apply_vector(axis.direction.vector()),
                        );
                        (c.cylinder().radius() - radius).abs() < T.confusion()
                            && (centre - o).cross(d).magnitude() < T.confusion()
                    })
                    .count();
                assert_eq!(found, 1, "the round at x = {x} is its exact cylinder");
            }
        }
        let want = section_area(wall, pitch, radius) * height;
        for v in &volumes {
            assert!((v - want).abs() < 1e-6, "{v} against {want}");
        }
        assert!((volumes[0] - volumes[1]).abs() < 1e-9, "either order");
    }
}

/// A concave seam between two fused walls standing on a plate whose top
/// leans off square to the seam (2 in 100, and 1 across), a slot run under the top past
/// the ball's contacts: the round is local (the slot's ceiling stands
/// behind the plate's top, which the round ends on, though above the plane
/// square to the seam through its foot), its foot cut by the leaning top
/// exactly, and built by the engine with no boolean, its fill the
/// section's times the seam's length within 1e-3. Refused before the
/// leaning end was read: the ceiling stood in the square-ended box, and
/// the strip's end was taken flush only on a plane square to it.
#[test]
fn a_concave_round_ends_on_a_leaning_floor_beside_a_hole() {
    let (wall, pitch, height, radius, lean) = (3.58_f64, 5.5_f64, 10.0, 1.0, 0.02);
    let mut model = Model::new();
    let mut lobes = Vec::new();
    for i in 0..2 {
        let frame = Frame::new(
            Point::new(f64::from(i) * pitch, 0.0, -1.0),
            Direction::Z,
            Direction::X,
            T,
        )
        .unwrap();
        lobes.push(
            ogeom_algo::make_cylinder(&mut model, frame, wall, height, T)
                .unwrap()
                .shape,
        );
    }
    let two = ogeom_bool::fuse(&mut model, &lobes[0], &lobes[1], T)
        .unwrap()
        .shape;
    // The plate's top `z = -lean (y + (x - pitch / 2) / 2)`, through the
    // seam's line at `y = 0`: leaning across the seam too, so that the
    // round's two corners on it stand at two heights.
    let up = Direction::new(ogeom_math::Vector::new(0.5 * lean, lean, 1.0), T).unwrap();
    let frame = Frame::new(Point::ORIGIN, up, Direction::X, T).unwrap();
    let (x, y, z) = (frame.x().vector(), frame.y().vector(), up.vector());
    let corner = Point::new(0.5 * pitch, 0.0, 0.0) - z * 3.0 - x * 10.0 - y * 10.0;
    let plate = ogeom_algo::make_box(
        &mut model,
        Frame::new(corner, up, Direction::X, T).unwrap(),
        (20.0, 20.0, 3.0),
        T,
    )
    .unwrap()
    .shape;
    let stood = ogeom_bool::fuse(&mut model, &two, &plate, T).unwrap().shape;
    let seam_y = -(wall * wall - pitch * pitch / 4.0).sqrt();
    // A slot run in from the plate's side under its top, a skin of 0.01
    // over it, from past the ball's contacts on out: its ceiling within the
    // box the round can take up, behind the top.
    let skin = 0.01;
    let mouth = Point::new(0.5 * pitch, 0.0, 0.0) - x * 1.0 - y * 11.0 - z * (skin + 0.5);
    let slot = ogeom_algo::make_box(
        &mut model,
        Frame::new(mouth, up, Direction::X, T).unwrap(),
        (2.0, 11.0 + seam_y - 0.8 * radius, 0.5),
        T,
    )
    .unwrap()
    .shape;
    let solid = ogeom_bool::cut(&mut model, &stood, &slot, T).unwrap().shape;
    let point = |model: &Model, v: &Shape| {
        v.transform(model.datums())
            .unwrap()
            .apply(model.node(v).unwrap().data().as_vertex().unwrap().point)
    };
    let seam = explore_unique(&model, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            ogeom_algo::edge_vertices(&model, e)
                .unwrap()
                .is_some_and(|(a, b)| {
                    let (a, b) = (point(&model, &a), point(&model, &b));
                    [a, b].iter().all(|p| {
                        (p.x - 0.5 * pitch).abs() < 1e-7 && (p.y - seam_y).abs() < 1e-7
                    }) && (a.z - b.z).abs() > 0.5 * height
                })
        })
        .expect("the concave seam");
    let volume = |model: &Model, s: &Shape| {
        ogeom_algo::volume_properties(model, s, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let before = volume(&model, &solid);
    let (built, work) = ogeom_bool::work::counted(|| {
        ogeom_fillet::fillet_edges(&mut model, &solid, &[seam], radius, T).unwrap()
    });
    assert_eq!(work.face_pairs, 0, "built with no boolean");
    let check = ogeom_algo::check(&model, &built.shape, T).unwrap();
    assert!(check.is_valid(), "{:?}", check.problems);
    // One fill: half the two the section holds past its bare discs.
    let half = 0.5 * pitch;
    let lens =
        2.0 * wall * wall * (half / wall).acos() - pitch * (wall * wall - half * half).sqrt();
    let bare = 3.0 * core::f64::consts::PI * wall * wall - 2.0 * lens;
    let fill = 0.5 * (section_area(wall, pitch, radius) - bare);
    // From the leaning floor under the fill's middle to the lobes' tops.
    let length = height - 1.0 + lean * (seam_y - 0.5 * radius);
    let grown = volume(&model, &built.shape) - before;
    assert!(
        (grown - fill * length).abs() < 1e-3 * fill * length,
        "{grown} against {}",
        fill * length
    );
}

/// A crescent's upper half, ten tall: inside a drum of radius 5 and outside
/// one of radius 3 standing in it tangent at `(5, 0)`, cut along `y = 0`.
/// Its horn runs out to nothing at `(5, 0)`: the two walls meet there
/// facing opposite ways, a cusp of the body as given. The corner where the
/// outer drum meets the cut at `(-5, 0)` rounds at 0.5 through the engine
/// with no boolean, valid, the cusp left as it was: the guard against faces
/// the blend folds flat reads a fold the given body already holds there as
/// none of the blend's (it refused the round before).
#[test]
fn a_round_beside_a_cusp_of_the_given_body_is_not_refused_for_it() {
    use ogeom_geom::{CircleCurve, Curve, LineCurve, PlaneSurface};
    use ogeom_math::{Circle, Plane, Vector};
    let mut model = Model::new();
    let (tip, outer_end, inner_end) = (
        Point::new(5.0, 0.0, 0.0),
        Point::new(-5.0, 0.0, 0.0),
        Point::new(-1.0, 0.0, 0.0),
    );
    let vertices: Vec<Shape> = [tip, outer_end, inner_end]
        .iter()
        .map(|&p| ogeom_algo::make_vertex(&mut model, p).shape)
        .collect();
    let pi = core::f64::consts::PI;
    // The outer drum's arc from the tip over to `(-5, 0)`, the cut, and the
    // inner drum's arc from `(-1, 0)` over back to the tip.
    let outer = Circle::new(Frame::WORLD, 5.0, T).unwrap();
    let inner = Circle::new(
        Frame::new(
            Point::new(2.0, 0.0, 0.0),
            -Direction::Z,
            -Direction::X,
            T,
        )
        .unwrap(),
        3.0,
        T,
    )
    .unwrap();
    let sides = [
        (Curve::Circle(CircleCurve::new(outer)), (0.0, pi), 0, 1),
        (
            Curve::Line(LineCurve::segment(outer_end, inner_end, T).unwrap()),
            (0.0, 4.0),
            1,
            2,
        ),
        (Curve::Circle(CircleCurve::new(inner)), (0.0, pi), 2, 0),
    ];
    let mut profile = Vec::new();
    for (curve, range, from, to) in sides {
        profile.push(
            ogeom_algo::make_edge_between(
                &mut model,
                curve,
                range,
                &vertices[from],
                &vertices[to],
                T,
            )
            .unwrap()
            .shape,
        );
    }
    let wire = ogeom_algo::make_wire(&mut model, &profile, T)
        .unwrap()
        .shape;
    let face = ogeom_algo::make_face(
        &mut model,
        PlaneSurface::new(Plane::new(Frame::WORLD)).into(),
        &[wire],
        T,
    )
    .unwrap()
    .shape;
    let solid = ogeom_algo::make_prism(&mut model, &face, Vector::new(0.0, 0.0, 10.0), T)
        .unwrap()
        .shape;
    let given = ogeom_algo::check(&model, &solid, T).unwrap();
    assert!(given.is_valid(), "{:?}", given.problems);
    let corner = explore_unique(&model, &solid, ShapeType::Edge)
        .unwrap()
        .into_iter()
        .find(|e| {
            let bb = ogeom_algo::shape_bounds(&model, e, T).unwrap();
            let (l, h) = (bb.low().unwrap(), bb.high().unwrap());
            (l.x + 5.0).abs() < 1e-6
                && (h.x + 5.0).abs() < 1e-6
                && l.y.abs() < 1e-6
                && h.y.abs() < 1e-6
                && (h.z - l.z - 10.0).abs() < 1e-6
        })
        .expect("the corner at (-5, 0)");
    let (built, work) = ogeom_bool::work::counted(|| {
        ogeom_fillet::fillet_edges(&mut model, &solid, &[corner], 0.5, T).unwrap()
    });
    assert_eq!(work.face_pairs, 0, "built with no boolean");
    let check = ogeom_algo::check(&model, &built.shape, T).unwrap();
    assert!(check.is_valid(), "{:?}", check.problems);
    // The spandrel taken off: a quarter-turn corner between the cut and a
    // drum of 5, its ball's centre at `(-4.5 cos a, 0.5)` off the axis with
    // `|c| = 4.5`, the foot on the cut `(x0, 0)` and on the drum `c 5 / 4.5`.
    let cx = -(4.5_f64 * 4.5 - 0.25).sqrt();
    let foot = Point::new(cx * 5.0 / 4.5, 0.5 * 5.0 / 4.5, 0.0);
    let (p, l) = (outer_end, Point::new(cx, 0.0, 0.0));
    let triangle = 0.5 * ((l.x - p.x) * (foot.y - p.y) - (foot.x - p.x) * (l.y - p.y)).abs();
    let segment = |radius: f64, cos: f64| {
        let phi = cos.clamp(-1.0, 1.0).acos();
        0.5 * radius * radius * (phi - phi.sin())
    };
    let wall = segment(5.0, (foot.x * p.x + foot.y * p.y) / 25.0);
    let ball = {
        let (u, w) = (l - Point::new(cx, 0.5, 0.0), foot - Point::new(cx, 0.5, 0.0));
        segment(0.5, u.dot(w) / 0.25)
    };
    let spandrel = triangle + wall - ball;
    let volume = |s: &Shape| {
        ogeom_algo::volume_properties(&model, s, ogeom_mesh::Deflection::default(), T)
            .unwrap()
            .mass
    };
    let taken = volume(&solid) - volume(&built.shape);
    assert!(
        (taken - spandrel * 10.0).abs() < 1e-6,
        "{taken} against {}",
        spandrel * 10.0
    );
}
