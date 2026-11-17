use crate::shared::{Body, ORIGIN, System, Vec3, gravity};

const MIN_CELL_SIZE: f64 = 1e-10;
const DIAGONAL: f64 = 3.4641016151377544; // 2 * sqrt(3)

pub struct Tree<'a> {
    system: System<'a>,
    mass: f64,
    com: Vec3,
    children: Vec<Tree<'a>>,
}
impl<'a> Tree<'a> {
    pub fn new(system: System<'a>) -> Tree<'a> {
        let mass: f64 = system.find_total_mass();
        let com: Vec3 = system.find_com();

        if should_split(&system) {
            let sys: System = System::new(system.center, system.hw);

            let parts: [System; 8] = system.partition();
            let mut children: Vec<Tree<'a>> = Vec::with_capacity(parts.len());
            for part in parts {
                children.push(Tree::new(part))
            }

            return Tree {
                system: sys,
                mass,
                com,
                children,
            };
        } else {
            Tree {
                system,
                mass,
                com,
                children: Vec::new(),
            }
        }
    }

    fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }

    pub fn accel_at(&self, pos: &Vec3, resolution: f64) -> Vec3 {
        // if there's no mass, there's no acceleration due to gravity in this region
        if self.mass == 0.0 {
            return ORIGIN;
        }

        // if this node is a leaf, it has one body that causes all the gravitational acceleration
        if self.is_leaf() {
            let only_body: &Body = self.system.bodies[0];
            if only_body.pos != *pos {
                return gravity(only_body, pos);
            } else {
                return ORIGIN;
            }
        }

        // if this region is far enough away, treat it as a point mass
        let width: f64 = 2.0_f64 * self.system.hw;
        let distance: f64 = self
            .com
            .iter()
            .zip(pos)
            .map(|xy| (xy.0 - xy.1).powi(2))
            .sum::<f64>()
            .sqrt();
        if width <= resolution * distance {
            let body: Body = Body {
                mass: self.mass,
                pos: self.com,
                ..Body::new()
            };
            return gravity(&body, pos);
        }

        // if this region is close and big enough to the point, then recurse to its children
        let mut acc: Vec3 = ORIGIN;
        for child in self.children.iter() {
            let local_acc: Vec3 = child.accel_at(pos, resolution);
            for (a, la) in acc.iter_mut().zip(local_acc) {
                *a += la;
            }
        }
        return acc;
    }
}

fn should_split(system: &System) -> bool {
    let is_big_enough: bool = DIAGONAL * system.hw >= MIN_CELL_SIZE;
    let enough_bodies: bool = system.bodies.len() > 1_usize;
    return is_big_enough && enough_bodies;
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::shared::test_util::*;
    use crate::shared::{Body, ORIGIN};

    // ---------------------------------------------------------------- helpers

    /// Every leaf under `tree`, in depth-first order.
    fn leaves<'t, 'a>(tree: &'t Tree<'a>) -> Vec<&'t Tree<'a>> {
        if tree.is_leaf() {
            vec![tree]
        } else {
            tree.children.iter().flat_map(|c| leaves(c)).collect()
        }
    }

    /// Every node (internal and leaf) under and including `tree`.
    fn all_nodes<'t, 'a>(tree: &'t Tree<'a>) -> Vec<&'t Tree<'a>> {
        let mut out = vec![tree];
        for child in &tree.children {
            out.extend(all_nodes(child));
        }
        out
    }

    /// Number of levels below `tree` (a lone leaf has depth 0).
    fn depth(tree: &Tree) -> usize {
        tree.children
            .iter()
            .map(|c| 1 + depth(c))
            .max()
            .unwrap_or(0)
    }

    /// `scattered(n)` with distinct masses, so mass weighting actually matters.
    fn weighted(n: usize) -> Vec<Body> {
        scattered(n)
            .iter()
            .enumerate()
            .map(|(i, b)| body_with(1.0 + i as f64, b.pos[0], b.pos[1], b.pos[2]))
            .collect()
    }

    // ------------------------------------------------------------ leaf cases

    #[test]
    fn an_empty_system_is_an_empty_leaf() {
        let tree = Tree::new(System::new(ORIGIN, 1.0));
        assert!(tree.is_leaf());
        assert_eq!(tree.mass, 0.0);
        assert!(tree.system.bodies.is_empty());
        assert!(
            tree.com.iter().all(|x| x.is_finite()),
            "com was {:?}",
            tree.com
        );
    }

    /// Most octants of a real tree are empty, so empty leaves are everywhere.
    /// A 0/0 in any mass or com calculation would put a NaN into the force pass.
    #[test]
    fn no_node_anywhere_has_a_nan_or_infinite_mass_or_com() {
        let bodies = weighted(40);
        let tree = Tree::new(System::from_bodies(&bodies));
        let empties = leaves(&tree)
            .iter()
            .filter(|l| l.system.bodies.is_empty())
            .count();
        assert!(empties > 0, "test is vacuous: the tree has no empty leaves");
        for node in all_nodes(&tree) {
            assert!(node.mass.is_finite());
            assert!(node.com.iter().all(|x| x.is_finite()), "com {:?}", node.com);
        }
    }

    #[test]
    fn a_single_body_is_a_leaf_at_that_body() {
        let bodies = [body_with(7.0, 1.0, -2.0, 3.0)];
        let tree = Tree::new(System::from_bodies(&bodies));
        assert!(tree.is_leaf());
        assert!(close(tree.mass, 7.0));
        assert!(close3(tree.com, [1.0, -2.0, 3.0]));
        assert_eq!(tree.system.bodies.len(), 1);
        assert!(std::ptr::eq(tree.system.bodies[0], &bodies[0]));
    }

    /// Opposite corners: one body per leaf, in two different leaves, with the
    /// other six octants of the root empty.
    #[test]
    fn two_separated_bodies_land_in_different_single_body_leaves() {
        let bodies = [body_at(-1.0, -1.0, -1.0), body_at(1.0, 1.0, 1.0)];
        let tree = Tree::new(System::from_bodies(&bodies));
        assert!(!tree.is_leaf());
        let occupied: Vec<&Tree> = tree
            .children
            .iter()
            .filter(|c| !c.system.bodies.is_empty())
            .collect();
        assert_eq!(occupied.len(), 2);
        for leaf in occupied {
            assert!(leaf.is_leaf());
            assert_eq!(leaf.system.bodies.len(), 1);
        }
    }

    // ------------------------------------------------------------- the shape

    #[test]
    fn internal_nodes_have_eight_children_and_no_bodies() {
        let bodies = weighted(40);
        let tree = Tree::new(System::from_bodies(&bodies));
        for node in all_nodes(&tree).into_iter().filter(|n| !n.is_leaf()) {
            assert_eq!(node.children.len(), 8);
            assert!(
                node.system.bodies.is_empty(),
                "an internal node kept its bodies"
            );
        }
    }

    #[test]
    fn leaves_above_the_resolution_limit_hold_at_most_one_body() {
        let bodies = weighted(40);
        let tree = Tree::new(System::from_bodies(&bodies));
        for leaf in leaves(&tree) {
            assert!(
                leaf.system.bodies.len() <= 1,
                "leaf with {} bodies",
                leaf.system.bodies.len()
            );
        }
    }

    // ------------------------------------------------------ bodies conserved

    #[test]
    fn every_body_is_in_exactly_one_leaf() {
        let bodies = weighted(40);
        let tree = Tree::new(System::from_bodies(&bodies));
        let found: Vec<&Tree> = leaves(&tree);
        for (i, body) in bodies.iter().enumerate() {
            let appearances = found
                .iter()
                .flat_map(|leaf| leaf.system.bodies.iter())
                .filter(|b| std::ptr::eq(**b, body))
                .count();
            assert_eq!(appearances, 1, "body {i}");
        }
    }

    #[test]
    fn every_body_lies_inside_its_leaf() {
        let bodies = weighted(40);
        let tree = Tree::new(System::from_bodies(&bodies));
        for leaf in leaves(&tree) {
            for body in &leaf.system.bodies {
                assert!(inside(leaf.system.center, leaf.system.hw, body.pos));
            }
        }
    }

    // ---------------------------------------------------------- mass and com

    #[test]
    fn root_mass_is_the_total_mass() {
        let bodies = weighted(40);
        let total: f64 = bodies.iter().map(|b| b.mass).sum();
        let tree = Tree::new(System::from_bodies(&bodies));
        assert!(close(tree.mass, total), "mass was {}", tree.mass);
    }

    /// Worked by hand: total mass 8, x = 3*4/8 = 1.5, y = 4*2/8 = 1, z = 0.
    /// The expected value must not come from `find_com`, or this would pass
    /// for any `find_com`, however wrong.
    #[test]
    fn root_com_matches_a_hand_computed_value() {
        let bodies = [
            body_with(1.0, 0.0, 0.0, 0.0),
            body_with(3.0, 4.0, 0.0, 0.0),
            body_with(4.0, 0.0, 2.0, 0.0),
        ];
        let tree = Tree::new(System::from_bodies(&bodies));
        assert!(!tree.is_leaf());
        assert!(close(tree.mass, 8.0));
        assert!(close3(tree.com, [1.5, 1.0, 0.0]), "com was {:?}", tree.com);
    }

    #[test]
    fn an_internal_nodes_mass_is_the_sum_of_its_children() {
        let bodies = weighted(40);
        let tree = Tree::new(System::from_bodies(&bodies));
        for node in all_nodes(&tree).into_iter().filter(|n| !n.is_leaf()) {
            let sum: f64 = node.children.iter().map(|c| c.mass).sum();
            assert!(close(node.mass, sum), "{} vs {}", node.mass, sum);
        }
    }

    #[test]
    fn an_internal_nodes_com_is_the_mass_weighted_com_of_its_children() {
        let bodies = weighted(40);
        let tree = Tree::new(System::from_bodies(&bodies));
        for node in all_nodes(&tree).into_iter().filter(|n| !n.is_leaf()) {
            let mut weighted_sum = ORIGIN;
            for child in &node.children {
                for d in 0..3 {
                    weighted_sum[d] += child.mass * child.com[d];
                }
            }
            let expected = weighted_sum.map(|x| x / node.mass);
            assert!(
                close3(node.com, expected),
                "{:?} vs {:?}",
                node.com,
                expected
            );
        }
    }

    // ----------------------------------------------- a root away from the origin

    #[test]
    fn a_root_far_from_the_origin_still_holds_every_body() {
        let bodies: Vec<Body> = weighted(40)
            .iter()
            .map(|b| body_with(b.mass, b.pos[0] + 50.0, b.pos[1] - 30.0, b.pos[2] + 10.0))
            .collect();
        let tree = Tree::new(System::from_bodies(&bodies));
        let in_leaves: usize = leaves(&tree).iter().map(|l| l.system.bodies.len()).sum();
        assert_eq!(in_leaves, bodies.len());
        for leaf in leaves(&tree) {
            for body in &leaf.system.bodies {
                assert!(inside(leaf.system.center, leaf.system.hw, body.pos));
            }
        }
    }

    // ------------------------------------------------------- resolution limit

    /// Two bodies closer than MIN_CELL_SIZE plus one far away. The far body
    /// keeps the root large, so the pair has to be chased down through many
    /// levels until the cell drops below the limit.
    ///
    /// The root is built by hand at the origin so the pair sits on the lower
    /// corner of every cell it falls into. A mean-based root (`from_bodies`)
    /// can put a cell boundary between two very close bodies and separate
    /// them early, which is correct behaviour but makes this test flaky.
    #[test]
    fn bodies_closer_than_the_resolution_limit_share_a_leaf() {
        let bodies = [
            body_at(0.0, 0.0, 0.0),
            body_at(1e-11, 0.0, 0.0),
            body_at(0.5, 0.5, 0.5),
        ];
        let mut root = System::new(ORIGIN, 1.0);
        root.add_bodies(&bodies);
        let tree = Tree::new(root);
        let sizes: Vec<usize> = leaves(&tree)
            .iter()
            .map(|l| l.system.bodies.len())
            .collect();
        assert_eq!(sizes.iter().sum::<usize>(), 3);
        assert!(sizes.contains(&2), "leaf sizes were {sizes:?}");
    }

    /// A leaf that merges two bodies acts as one point mass: their total mass,
    /// at their mass-weighted center. Same hand-built geometry as the test above
    /// (distinct bodies 1e-11 apart; exact coincidence is outside the input
    /// contract). The expected com is 3 * 1e-11 / 5 = 6e-12, not the unweighted
    /// midpoint 5e-12, so the check is on the weighting.
    #[test]
    fn a_merged_leaf_acts_as_one_point_mass() {
        let bodies = [
            body_with(2.0, 0.0, 0.0, 0.0),
            body_with(3.0, 1e-11, 0.0, 0.0),
            body_at(0.5, 0.5, 0.5),
        ];
        let mut root = System::new(ORIGIN, 1.0);
        root.add_bodies(&bodies);
        let tree = Tree::new(root);
        let pair = leaves(&tree)
            .into_iter()
            .find(|l| l.system.bodies.len() == 2)
            .expect("the close pair should share a leaf");
        assert!(close(pair.mass, 5.0));
        // `close3` would pass here for any com: everything is below its 1e-9 tolerance
        assert!(
            (pair.com[0] - 6e-12).abs() < 1e-18,
            "com was {:?}",
            pair.com
        );
        assert_eq!(pair.com[1], 0.0);
        assert_eq!(pair.com[2], 0.0);
    }

    /// Same hand-built geometry as the test above, but the gap is 1e-9, ten
    /// times the 1e-10 limit, so some cell boundary must fall between the pair
    /// and separate them. The gap is a literal on purpose: written as a multiple
    /// of MIN_CELL_SIZE it would move with the constant, and a limit of 1e-6
    /// would go unnoticed.
    #[test]
    fn bodies_ten_times_the_limit_apart_are_resolved_individually() {
        let bodies = [
            body_at(0.0, 0.0, 0.0),
            body_at(1e-9, 0.0, 0.0),
            body_at(0.5, 0.5, 0.5),
        ];
        let mut root = System::new(ORIGIN, 1.0);
        root.add_bodies(&bodies);
        let tree = Tree::new(root);
        assert!(
            leaves(&tree).iter().all(|l| l.system.bodies.len() <= 1),
            "bodies at 10x the limit should be resolved individually"
        );
    }

    /// A pair too close to separate is the worst case: only the resolution
    /// limit stops it. In a root of half-width 1 the cell diagonal drops below
    /// 1e-10 after about log2(2*sqrt(3) / 1e-10) = 35 halvings, so the depth is
    /// bounded and lands close to that figure. (The pair is 1e-11 apart, not
    /// coincident: exact coincidence is outside the input contract.)
    #[test]
    fn recursion_depth_is_set_by_the_resolution_limit() {
        let bodies = [
            body_at(0.0, 0.0, 0.0),
            body_at(1e-11, 0.0, 0.0),
            body_at(0.5, 0.5, 0.5),
        ];
        let mut root = System::new(ORIGIN, 1.0);
        root.add_bodies(&bodies);
        let d = depth(&Tree::new(root));
        assert!((33..=37).contains(&d), "depth was {d}");
    }

    // -------------------------------------------------------------- accel_at

    /// Relative difference between two vectors. G is ~1e-11, so an absolute
    /// tolerance like `close3`'s would pass for any vector.
    fn rel_diff(a: Vec3, b: Vec3) -> f64 {
        let diff: f64 = (0..3).map(|k| (a[k] - b[k]).powi(2)).sum::<f64>().sqrt();
        let size: f64 = (0..3).map(|k| b[k].powi(2)).sum::<f64>().sqrt();
        diff / size
    }

    /// A two-body cluster seen from far away. The root cell has width 2 and the
    /// target sits about 9.75 from the cluster's center of mass, so the
    /// width/distance ratio is about 0.205. θ = 0.25 is just above that (accept
    /// the node as one point mass) and θ = 0.15 is below it (open it and sum the
    /// bodies). The expected answers are built from `gravity` directly, so this
    /// pins the exact switch-over point, not just the trend.
    #[test]
    fn a_node_is_accepted_exactly_when_width_over_distance_is_below_theta() {
        let bodies = [
            body_with(1.0, -0.5, -0.5, -0.5),
            body_with(3.0, 0.5, 0.5, 0.5),
        ];
        let mut root = System::new(ORIGIN, 1.0);
        root.add_bodies(&bodies);
        let tree = Tree::new(root);
        assert!(!tree.is_leaf());
        let target = [10.0, 0.1, 0.2];

        let exact = {
            let a = gravity(&bodies[0], &target);
            let b = gravity(&bodies[1], &target);
            [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
        };
        let as_point = gravity(&body_with(4.0, 0.25, 0.25, 0.25), &target);
        // the two answers must differ, or the test cannot tell them apart
        assert!(rel_diff(as_point, exact) > 1e-4);

        let accepted = tree.accel_at(&target, 0.25);
        assert!(
            rel_diff(accepted, as_point) < 1e-9,
            "should be the point mass"
        );

        let opened = tree.accel_at(&target, 0.15);
        assert!(rel_diff(opened, exact) < 1e-9, "should be the exact sum");
    }

    /// A target that is itself in the tree must not feel itself.
    #[test]
    fn a_body_in_the_tree_feels_only_the_others() {
        let bodies = [
            body_with(2.0, -1.0, 0.0, 0.0),
            body_with(3.0, 1.0, 0.5, 0.25),
        ];
        let tree = Tree::new(System::from_bodies(&bodies));
        let felt = tree.accel_at(&bodies[0].pos, 0.0);
        assert!(rel_diff(felt, gravity(&bodies[1], &bodies[0].pos)) < 1e-9);
    }

    // ----------------------------------------------------------- should_split

    #[test]
    fn should_split_needs_more_than_one_body_and_a_big_enough_cell() {
        let two = [body_at(0.0, 0.0, 0.0), body_at(1.0, 0.0, 0.0)];
        let one = [body_at(0.0, 0.0, 0.0)];

        let mut big_two = System::new(ORIGIN, 1.0);
        big_two.add_bodies(&two);
        let mut tiny_two = System::new(ORIGIN, MIN_CELL_SIZE / 10.0);
        tiny_two.add_bodies(&two);
        let mut big_one = System::new(ORIGIN, 1.0);
        big_one.add_bodies(&one);

        assert!(should_split(&big_two));
        assert!(!should_split(&tiny_two), "below the limit must stop");
        assert!(!should_split(&big_one), "a single body is a leaf");
        assert!(!should_split(&System::new(ORIGIN, 1.0)), "empty is a leaf");
    }

    /// DIAGONAL is a typed-in literal for 2*sqrt(3); check it against the real thing.
    #[test]
    fn diagonal_is_twice_root_three() {
        assert!((DIAGONAL - 2.0 * 3.0_f64.sqrt()).abs() < 1e-12);
    }

    /// The cutoff is "cell diagonal < MIN_CELL_SIZE", so it sits at
    /// hw = MIN_CELL_SIZE / (2*sqrt(3)). Just above splits, just below stops.
    #[test]
    fn the_cutoff_is_on_the_cell_diagonal() {
        let cutoff = MIN_CELL_SIZE / (2.0 * 3.0_f64.sqrt());
        let two = [body_at(0.0, 0.0, 0.0), body_at(1.0, 0.0, 0.0)];

        let mut above = System::new(ORIGIN, cutoff * 1.01);
        above.add_bodies(&two);
        let mut below = System::new(ORIGIN, cutoff * 0.99);
        below.add_bodies(&two);

        assert!(should_split(&above));
        assert!(!should_split(&below));
    }
}
