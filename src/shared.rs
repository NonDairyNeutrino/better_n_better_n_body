//! This is a module documentation

use std::array::from_fn;

pub const G: f64 = 6.67 * 1e-11;
pub const ORIGIN: [f64; 3] = [0.0_f64, 0.0_f64, 0.0_f64];

/// This is the docstring for Object
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Body {
    pub mass: f64,
    pub pos: [f64; 3],
    pub vel: [f64; 3],
}
impl Body {
    /// Creates a new physical Object
    pub fn new() -> Body {
        return Body {
            pos: [0f64, 0f64, 0f64],
            vel: [0f64, 0f64, 0f64],
            mass: 1f64,
        };
    }
    /// Calculates the Euclidean distance between an object and a position
    pub fn distance(&self, pos: &[f64; 3]) -> f64 {
        self.pos
            .iter()
            .zip(pos)
            .map(|pair| (*pair.0 - pair.1).powi(2))
            .sum::<f64>()
            .sqrt()
    }
    pub fn accelerate(&mut self, acc: &[f64; 3], dt: f64) {
        // evolve via symplectic-euler
        for (i, a) in acc.iter().enumerate() {
            self.vel[i] += a * dt;
            self.pos[i] += self.vel[i] * dt;
        }
    }
}

/// A collection of bodies in a cubical space with a center and half-width
#[derive(Debug)]
pub struct System<'a> {
    pub bodies: Vec<&'a Body>,
    pub center: [f64; 3],
    pub hw: f64,
}
impl<'a> System<'a> {
    pub fn new(center: [f64; 3], half_width: f64) -> System<'a> {
        System {
            bodies: Vec::new(),
            center: center,
            hw: half_width,
        }
    }
    pub fn from_bodies(bodies: &'a [Body]) -> System<'a> {
        let bods: Vec<&'a Body> = bodies.iter().collect();
        let (center, hw): ([f64; 3], f64) = find_bounding_box(bodies);
        System {
            bodies: bods,
            center: center,
            hw: hw,
        }
    }
    pub fn add_body(&mut self, body: &'a Body) {
        self.bodies.push(body);
    }
    pub fn add_bodies(&mut self, bodies: &'a [Body]) {
        // know how much more space is needed, so allocate first, then push each
        self.bodies.reserve(bodies.len());
        for body in bodies {
            self.bodies.push(body);
        }
    }
    pub fn join(&mut self, system: &'a System) {
        self.bodies.extend_from_slice(&system.bodies);
    }
    pub fn find_com(&self) -> [f64; 3] {
        let mut com: [f64; 3] = ORIGIN;
        if self.bodies.is_empty() {
            return com;
        }

        // find the mass-weighted total
        for body in &self.bodies {
            for dim in 0..com.len() {
                com[dim] += body.mass * body.pos[dim];
            }
        }

        // normalize by total mass
        let total_mass: f64 = self.bodies.iter().map(|b| b.mass).sum();
        for dim in 0..com.len() {
            com[dim] /= total_mass;
        }
        return com;
    }
    pub fn partition(self) -> [System<'a>; 8] {
        let mut systems: [System; 8] =
            from_fn(|i| System::new(child_center(self.center, self.hw, i), self.hw / 2.0));

        for body in self.bodies {
            let idx: usize = child_index(self.center, body.pos);
            systems[idx].add_body(body);
        }

        return systems;
    }
}

fn child_index(center: [f64; 3], pos: [f64; 3]) -> usize {
    let mut idx: usize = 0;
    if pos[0] >= center[0] {
        idx += 1;
    }
    if pos[1] >= center[1] {
        idx += 2;
    }
    if pos[2] >= center[2] {
        idx += 4;
    }
    return idx;
}
fn child_center(center: [f64; 3], half_width: f64, index: usize) -> [f64; 3] {
    let mut chi_cen: [f64; 3] = ORIGIN;
    for dim in 0..center.len() {
        let bit: i32 = ((index / 2_usize.pow(dim as u32)) % 2) as i32;
        chi_cen[dim] = center[dim] + (-1.0_f64).powi(bit + 1) * (half_width / 2.0_f64);
    }
    return chi_cen;
}

fn find_center(system: &[Body]) -> [f64; 3] {
    let mut center: [f64; 3] = ORIGIN;
    let n_bodies: f64 = system.len() as f64;
    if system.is_empty() {
        return ORIGIN;
    }

    for obj in system {
        for dim in 0..3 {
            center[dim] += obj.pos[dim];
        }
    }

    // average
    for dim in 0..center.len() {
        center[dim] /= n_bodies;
    }
    return center;
}

fn find_half_width(system: &[Body], center: [f64; 3]) -> f64 {
    let mut half_width: f64 = 0.0_f64;
    for obj in system {
        let d: f64 = obj.distance(&center);
        if d > half_width {
            half_width = d;
        }
    }
    return half_width;
}

/// Find the center and half-width of a cube that contains the system.
pub fn find_bounding_box(system: &[Body]) -> ([f64; 3], f64) {
    let center: [f64; 3] = find_center(system);
    let width: f64 = find_half_width(system, center);
    return (center, width);
}

pub fn gravity(object: &Body, pos: &[f64; 3]) -> [f64; 3] {
    let strength: f64 = G * object.mass / object.distance(pos).powi(3); // 3 is on purpose
    let mut grav: [f64; 3] = [0.0_f64, 0.0_f64, 0.0_f64];

    for dim in 0..grav.len() {
        grav[dim] = strength * (object.pos[dim] - pos[dim]);
    }

    return grav;
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-9;

    // ---------------------------------------------------------------- helpers

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < TOL
    }

    fn close3(a: [f64; 3], b: [f64; 3]) -> bool {
        a.iter().zip(b.iter()).all(|(x, y)| close(*x, *y))
    }

    /// A unit-mass body at rest at the given position.
    fn body_at(x: f64, y: f64, z: f64) -> Body {
        Body {
            pos: [x, y, z],
            ..Body::new()
        }
    }

    fn body_with(mass: f64, x: f64, y: f64, z: f64) -> Body {
        Body {
            mass,
            pos: [x, y, z],
            ..Body::new()
        }
    }

    /// True if `pos` lies inside the cube (`center`, `hw`) on every axis.
    fn inside(center: [f64; 3], hw: f64, pos: [f64; 3]) -> bool {
        (0..3).all(|d| (pos[d] - center[d]).abs() <= hw + TOL)
    }

    /// Deterministic spread of `n` bodies, all strictly inside the cube [-1, 1]^3.
    fn scattered(n: usize) -> Vec<Body> {
        (0..n)
            .map(|i| {
                body_at(
                    ((i * 7) % 13) as f64 / 6.5 - 1.0,
                    ((i * 5) % 11) as f64 / 5.5 - 1.0,
                    ((i * 3) % 17) as f64 / 8.5 - 1.0,
                )
            })
            .collect()
    }

    // ------------------------------------------------------------------- Body

    #[test]
    fn distance_is_euclidean() {
        let body = Body {
            pos: [3.0, 4.0, 0.0],
            ..Body::new()
        };
        let d = body.distance(&[0.0, 0.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-9);
    }

    #[test]
    fn distance_uses_all_three_axes() {
        // (1, 2, 2) is exactly 3 from the origin
        let body = body_at(1.0, 2.0, 2.0);
        assert!(close(body.distance(&ORIGIN), 3.0));
        // and the offset is measured from the given point, not from the origin
        assert!(close(body.distance(&[1.0, 2.0, 0.0]), 2.0));
    }

    #[test]
    fn accelerate_without_acceleration_moves_at_constant_velocity() {
        let mut body = Body {
            vel: [1.0, 2.0, 3.0],
            ..Body::new()
        };
        body.accelerate(&[0.0, 0.0, 0.0], 2.0);
        assert!(close3(body.pos, [2.0, 4.0, 6.0]));
        assert!(close3(body.vel, [1.0, 2.0, 3.0]));
    }

    /// The comment in `accelerate` says symplectic Euler: update the velocity
    /// first, then move by the NEW velocity (x += v_new * dt). dt is deliberately
    /// not 1, so a wrong power of dt shows up, and every axis is exercised.
    #[test]
    fn accelerate_position_follows_symplectic_euler() {
        let mut body = Body::new();
        body.accelerate(&[2.0, 4.0, 6.0], 0.5);
        assert!(close3(body.vel, [1.0, 2.0, 3.0]));
        assert!(close3(body.pos, [0.5, 1.0, 1.5]), "pos was {:?}", body.pos);
    }

    // ---------------------------------------------------------------- gravity

    #[test]
    fn gravity_points_toward_the_source() {
        let source = body_at(2.0, 0.0, 0.0);
        let g = gravity(&source, &ORIGIN);
        assert!(g[0] > 0.0);
        // divide by G: gravity here is ~1e-11, far below `close`'s absolute tolerance,
        // so comparing the raw components to zero would pass for any vector
        assert!(close(g[1] / G, 0.0));
        assert!(close(g[2] / G, 0.0));
    }

    #[test]
    fn gravity_magnitude_is_g_m_over_r_squared() {
        let source = body_with(5.0, 0.0, 3.0, 0.0);
        let g = gravity(&source, &ORIGIN);
        // divide out G so the comparison isn't dwarfed by the 1e-11 scale
        assert!(close(g[1] / G, 5.0 / 9.0), "g/G was {}", g[1] / G);
    }

    #[test]
    fn gravity_has_the_right_components_in_3d() {
        // source at (1, 2, 2): distance 3, so |g| = G/9 along (1/3, 2/3, 2/3)
        let g = gravity(&body_at(1.0, 2.0, 2.0), &ORIGIN);
        assert!(close(g[0] / G, 1.0 / 27.0), "gx/G was {}", g[0] / G);
        assert!(close(g[1] / G, 2.0 / 27.0), "gy/G was {}", g[1] / G);
        assert!(close(g[2] / G, 2.0 / 27.0), "gz/G was {}", g[2] / G);
    }

    // ----------------------------------------------------------------- System

    #[test]
    fn from_bodies_references_every_body_in_master_order() {
        let master = scattered(6);
        let system = System::from_bodies(&master);
        assert_eq!(system.bodies.len(), master.len());
        for (i, body) in system.bodies.iter().enumerate() {
            assert!(std::ptr::eq(*body, &master[i]), "body {i} out of order");
        }
    }

    #[test]
    fn add_body_and_add_bodies_append() {
        let first = body_at(0.0, 0.0, 0.0);
        let rest = [body_at(1.0, 0.0, 0.0), body_at(2.0, 0.0, 0.0)];
        let mut system = System::new(ORIGIN, 1.0);
        system.add_body(&first);
        system.add_bodies(&rest);
        assert_eq!(system.bodies.len(), 3);
        assert!(std::ptr::eq(system.bodies[2], &rest[1]));
    }

    #[test]
    fn add_body_does_not_change_the_region() {
        let body = body_at(100.0, 100.0, 100.0);
        let mut system = System::new([1.0, 2.0, 3.0], 4.0);
        system.add_body(&body);
        assert_eq!(system.center, [1.0, 2.0, 3.0]);
        assert_eq!(system.hw, 4.0);
    }

    #[test]
    fn join_appends_the_other_systems_bodies() {
        let a = [body_at(0.0, 0.0, 0.0)];
        let b = [body_at(1.0, 0.0, 0.0), body_at(2.0, 0.0, 0.0)];
        let other = System::from_bodies(&b);
        let mut system = System::from_bodies(&a);
        system.join(&other);
        assert_eq!(system.bodies.len(), 3);
    }

    // -------------------------------------------------------- center of mass

    #[test]
    fn com_is_mass_weighted() {
        let bodies = [body_with(1.0, 0.0, 0.0, 0.0), body_with(3.0, 4.0, 0.0, 0.0)];
        let system = System::from_bodies(&bodies);
        assert!(close3(system.find_com(), [3.0, 0.0, 0.0]));
    }

    #[test]
    fn com_of_equal_masses_is_the_midpoint() {
        let bodies = [body_at(0.0, 0.0, 0.0), body_at(2.0, 4.0, 6.0)];
        let system = System::from_bodies(&bodies);
        assert!(close3(system.find_com(), [1.0, 2.0, 3.0]));
    }

    #[test]
    fn com_of_a_single_body_is_its_position() {
        let bodies = [body_with(7.0, 1.0, -2.0, 3.0)];
        let system = System::from_bodies(&bodies);
        assert!(close3(system.find_com(), [1.0, -2.0, 3.0]));
    }

    #[test]
    fn com_of_an_empty_system_is_the_origin() {
        let system = System::new([5.0, 5.0, 5.0], 1.0);
        assert_eq!(system.find_com(), ORIGIN);
    }

    // ----------------------------------------------------------- bounding box

    fn assert_box_contains_all(bodies: &[Body]) {
        let (center, hw) = find_bounding_box(bodies);
        for body in bodies {
            assert!(
                inside(center, hw, body.pos),
                "{:?} is outside cube center={center:?} hw={hw}",
                body.pos
            );
        }
    }

    #[test]
    fn bounding_box_contains_an_asymmetric_system() {
        // the case where an average-based center and an origin-based width disagree
        assert_box_contains_all(&[
            body_at(10.0, 0.0, 0.0),
            body_at(-10.0, 0.0, 0.0),
            body_at(-10.0, 0.0, 0.0),
        ]);
    }

    #[test]
    fn bounding_box_of_a_single_body_is_a_point() {
        let (center, hw) = find_bounding_box(&[body_at(1.0, 2.0, 3.0)]);
        assert!(close3(center, [1.0, 2.0, 3.0]));
        assert!(close(hw, 0.0));
    }

    #[test]
    fn bounding_box_of_nothing_is_the_origin() {
        let (center, hw) = find_bounding_box(&[]);
        assert_eq!(center, ORIGIN);
        assert_eq!(hw, 0.0);
    }

    // ------------------------------------------------------------ child_index

    #[test]
    fn child_index_labels_all_eight_octants() {
        let cases: [([f64; 3], usize); 8] = [
            ([-1.0, -1.0, -1.0], 0),
            ([1.0, -1.0, -1.0], 1),
            ([-1.0, 1.0, -1.0], 2),
            ([1.0, 1.0, -1.0], 3),
            ([-1.0, -1.0, 1.0], 4),
            ([1.0, -1.0, 1.0], 5),
            ([-1.0, 1.0, 1.0], 6),
            ([1.0, 1.0, 1.0], 7),
        ];
        for (pos, expected) in cases {
            assert_eq!(child_index(ORIGIN, pos), expected, "pos {pos:?}");
        }
    }

    #[test]
    fn child_index_sends_the_exact_center_to_the_upper_octant() {
        assert_eq!(child_index([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]), 7);
    }

    #[test]
    fn child_index_is_relative_to_the_center() {
        // In each case ONE axis is below the center (10, 10, 10) while still being
        // above the origin, so comparing that axis to the origin gives the wrong octant.
        let c = [10.0; 3];
        assert_eq!(child_index(c, [9.0, 11.0, 11.0]), 6); // x below
        assert_eq!(child_index(c, [11.0, 9.0, 11.0]), 5); // y below
        assert_eq!(child_index(c, [11.0, 11.0, 9.0]), 3); // z below
    }

    // ----------------------------------------------------------- child_center

    #[test]
    fn child_center_places_octants_around_the_origin() {
        assert!(close3(child_center(ORIGIN, 1.0, 0), [-0.5, -0.5, -0.5]));
        assert!(close3(child_center(ORIGIN, 1.0, 5), [0.5, -0.5, 0.5]));
        assert!(close3(child_center(ORIGIN, 1.0, 7), [0.5, 0.5, 0.5]));
    }

    #[test]
    fn child_center_is_relative_to_the_parent() {
        assert!(close3(
            child_center([10.0, 20.0, 30.0], 4.0, 5),
            [12.0, 18.0, 32.0]
        ));
    }

    #[test]
    fn child_center_and_child_index_agree() {
        let center = [10.0, -20.0, 30.0];
        for i in 0..8 {
            let c = child_center(center, 4.0, i);
            assert_eq!(child_index(center, c), i, "octant {i} center {c:?}");
        }
    }

    // -------------------------------------------------------------- partition

    /// Two roots to run the partition tests on: one at the origin, one far from it.
    /// Tests that only use an origin-centered root can't notice a partition that
    /// forgets the parent's center.
    const ROOTS: [([f64; 3], f64); 2] = [(ORIGIN, 1.0), ([10.0, -5.0, 3.0], 4.0)];

    /// `scattered(n)` rescaled to fill the cube (`center`, `hw`).
    fn scattered_in(center: [f64; 3], hw: f64, n: usize) -> Vec<Body> {
        scattered(n)
            .iter()
            .map(|b| {
                body_at(
                    center[0] + hw * b.pos[0],
                    center[1] + hw * b.pos[1],
                    center[2] + hw * b.pos[2],
                )
            })
            .collect()
    }

    /// A body at the middle of octant `i` of the cube (`center`, `2 * hw`);
    /// bit 0 = x, bit 1 = y, bit 2 = z, and a set bit means "above".
    fn octant_body(center: [f64; 3], hw: f64, i: usize) -> Body {
        let d = hw / 2.0;
        let s = |bit: usize| if (i >> bit) & 1 == 1 { d } else { -d };
        body_at(center[0] + s(0), center[1] + s(1), center[2] + s(2))
    }

    #[test]
    fn partition_separates_the_eight_octants() {
        for (center, hw) in ROOTS {
            let bodies: Vec<Body> = (0..8).map(|i| octant_body(center, hw, i)).collect();
            let mut root = System::new(center, hw);
            for b in &bodies {
                root.add_body(b);
            }
            for (i, child) in root.partition().iter().enumerate() {
                assert_eq!(child.bodies.len(), 1, "root {center:?}, octant {i}");
                assert!(std::ptr::eq(child.bodies[0], &bodies[i]), "octant {i}");
                // the body sits at the middle of its octant, so that is the child's center
                assert!(
                    close3(child.center, bodies[i].pos),
                    "root {center:?}, octant {i}: child center {:?}",
                    child.center
                );
            }
        }
    }

    #[test]
    fn partition_keeps_every_body_exactly_once() {
        for (center, hw) in ROOTS {
            let master = scattered_in(center, hw, 40);
            let mut root = System::new(center, hw);
            for b in &master {
                root.add_body(b);
            }
            let children = root.partition();

            let total: usize = children.iter().map(|c| c.bodies.len()).sum();
            assert_eq!(total, master.len());
            for (i, body) in master.iter().enumerate() {
                let appearances = children
                    .iter()
                    .flat_map(|c| c.bodies.iter())
                    .filter(|b| std::ptr::eq(**b, body))
                    .count();
                assert_eq!(appearances, 1, "root {center:?}, body {i}");
            }
        }
    }

    #[test]
    fn partition_children_contain_their_own_bodies() {
        for (center, hw) in ROOTS {
            let master = scattered_in(center, hw, 40);
            let mut root = System::new(center, hw);
            for b in &master {
                root.add_body(b);
            }
            for (i, child) in root.partition().iter().enumerate() {
                for body in &child.bodies {
                    assert!(
                        inside(child.center, child.hw, body.pos),
                        "root {center:?}, octant {i}: {:?} outside center={:?} hw={}",
                        body.pos,
                        child.center,
                        child.hw
                    );
                }
            }
        }
    }

    #[test]
    fn partition_children_have_half_the_parents_half_width() {
        let root = System::new([1.0, 2.0, 3.0], 8.0);
        for child in root.partition().iter() {
            assert!(close(child.hw, 4.0));
        }
    }
}
