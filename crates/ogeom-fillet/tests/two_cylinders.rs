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
