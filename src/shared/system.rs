use super::{Body, ORIGIN, Vec3};
use std::array::from_fn;

/// A collection of bodies in a cubical space with a center and half-width.
///
/// A `System` borrows its bodies (it holds `&Body`, never a copy), so it is cheap
/// to split and rejoin. That is how the Barnes-Hut and FMM trees divide space:
/// [`partition`](System::partition) sorts the bodies of a cell into its eight
/// child cells without moving any body.
///
/// The cube is the region `center ± hw` on each axis. The bodies are expected to
/// lie inside it, but nothing checks that.
#[derive(Debug)]
pub struct System<'a> {
    /// The bodies in this region, borrowed from the caller's slice.
    pub bodies: Vec<&'a Body>,
    /// Center of the cube.
    pub center: Vec3,
    /// Half of the cube's side length.
    pub hw: f64,
}
impl<'a> System<'a> {
    /// Creates an empty system occupying the cube around `center` with half-width `half_width`.
    pub fn new(center: Vec3, half_width: f64) -> System<'a> {
        System {
            bodies: Vec::new(),
            center: center,
            hw: half_width,
        }
    }
    /// Creates a system holding every body in `bodies`, in a cube that contains all of them.
    ///
    /// See [`find_bounding_box`] for how the cube is chosen.
    pub fn from_bodies(bodies: &'a [Body]) -> System<'a> {
        let bods: Vec<&'a Body> = bodies.iter().collect();
        let (center, hw): (Vec3, f64) = find_bounding_box(bodies);
        System {
            bodies: bods,
            center: center,
            hw: hw,
        }
    }
    /// Adds one body to the system.
    pub fn add_body(&mut self, body: &'a Body) {
        self.bodies.push(body);
    }
    /// Adds every body in `bodies` to the system.
    pub fn add_bodies(&mut self, bodies: &'a [Body]) {
        // know how much more space is needed, so allocate first, then push each
        self.bodies.reserve(bodies.len());
        for body in bodies {
            self.bodies.push(body);
        }
    }
    /// Adds every body of `system` to this one. The cube of this system is unchanged.
    pub fn join(&mut self, system: &'a System) {
        self.bodies.extend_from_slice(&system.bodies);
    }
    /// The sum of the masses of all bodies in the system (0 if it is empty).
    pub fn find_total_mass(&self) -> f64 {
        self.bodies.iter().map(|b| b.mass).sum()
    }
    /// The center of mass: the average position of the bodies, weighted by mass.
    ///
    /// An empty system returns [`ORIGIN`].
    pub fn find_com(&self) -> Vec3 {
        let mut com: Vec3 = ORIGIN;
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
    /// Splits the system into its eight octants, consuming it.
    ///
    /// Each octant is a cube with half the half-width, centered a quarter of the
    /// side length away from the parent's center. The index of an octant is a
    /// 3-bit number: bit 0 is set for the `+x` side, bit 1 for `+y` and bit 2 for
    /// `+z`. A body exactly on a dividing plane goes to the `+` side.
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

/// Which octant (0..8) of the cube around `center` the point `pos` falls in. See [`System::partition`] for the numbering.
fn child_index(center: Vec3, pos: Vec3) -> usize {
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
/// The center of octant `index` of the cube (`center`, `half_width`), numbered as in [`System::partition`].
fn child_center(center: Vec3, half_width: f64, index: usize) -> Vec3 {
    let mut chi_cen: Vec3 = ORIGIN;
    for dim in 0..center.len() {
        let bit: i32 = ((index / 2_usize.pow(dim as u32)) % 2) as i32;
        chi_cen[dim] = center[dim] + (-1.0_f64).powi(bit + 1) * (half_width / 2.0_f64);
    }
    return chi_cen;
}

/// The average position of the bodies (not weighted by mass). `ORIGIN` for an empty slice.
fn find_center(system: &[Body]) -> Vec3 {
    let mut center: Vec3 = ORIGIN;
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

/// The largest distance from `center` to any body.
fn find_half_width(system: &[Body], center: Vec3) -> f64 {
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
///
/// The center is the average position of the bodies. The half-width is the
/// distance from that center to the farthest body. A cube with that half-width
/// contains the whole sphere of that radius, so it contains every body, with some
/// room to spare in the corners. An empty slice gives `(ORIGIN, 0.0)`.
pub fn find_bounding_box(system: &[Body]) -> (Vec3, f64) {
    let center: Vec3 = find_center(system);
    let width: f64 = find_half_width(system, center);
    return (center, width);
}

#[cfg(test)]
mod tests {
    use super::super::test_util::*;
    use super::*;

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
        let cases: [(Vec3, usize); 8] = [
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
    const ROOTS: [(Vec3, f64); 2] = [(ORIGIN, 1.0), ([10.0, -5.0, 3.0], 4.0)];

    /// `scattered(n)` rescaled to fill the cube (`center`, `hw`).
    fn scattered_in(center: Vec3, hw: f64, n: usize) -> Vec<Body> {
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
    fn octant_body(center: Vec3, hw: f64, i: usize) -> Body {
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
