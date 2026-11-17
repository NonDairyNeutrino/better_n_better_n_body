use crate::shared::{System, Vec3};

const MIN_CELL_SIZE: f64 = 1e-10;
const DIAGONAL: f64 = 3.4641016151377544; // 2 * sqrt(3)

pub struct Tree<'a> {
    system: System<'a>,
    mass: f64,
    com: Vec3,
    children: Vec<Tree<'a>>,
}
impl<'a> Tree<'a> {
    pub fn new(system: System<'a>) -> Tree {
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
        tree.children.iter().map(|c| 1 + depth(c)).max().unwrap_or(0)
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
        assert!(tree.com.iter().all(|x| x.is_finite()), "com was {:?}", tree.com);
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
        let sizes: Vec<usize> = leaves(&tree).iter().map(|l| l.system.bodies.len()).collect();
        assert_eq!(sizes.iter().sum::<usize>(), 3);
        assert!(sizes.contains(&2), "leaf sizes were {sizes:?}");
    }

    #[test]
    fn coincident_bodies_terminate_in_a_shared_leaf() {
        let bodies = [
            body_with(2.0, 0.0, 0.0, 0.0),
            body_with(3.0, 0.0, 0.0, 0.0),
            body_at(1.0, 1.0, 1.0),
        ];
        let tree = Tree::new(System::from_bodies(&bodies));
        let pair = leaves(&tree)
            .into_iter()
            .find(|l| l.system.bodies.len() == 2)
            .expect("the coincident pair should share a leaf");
        // the shared leaf acts as one point mass
        assert!(close(pair.mass, 5.0));
        assert!(close3(pair.com, [0.0, 0.0, 0.0]));
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

    /// A coincident pair is the worst case: it can only be stopped by the
    /// resolution limit. In a root of half-width 1 the cell diagonal drops
    /// below 1e-10 after about log2(2*sqrt(3) / 1e-10) = 35 halvings, so the
    /// depth is bounded and lands close to that figure.
    #[test]
    fn recursion_depth_is_set_by_the_resolution_limit() {
        let bodies = [
            body_at(0.0, 0.0, 0.0),
            body_at(0.0, 0.0, 0.0),
            body_at(0.5, 0.5, 0.5),
        ];
        let mut root = System::new(ORIGIN, 1.0);
        root.add_bodies(&bodies);
        let d = depth(&Tree::new(root));
        assert!((33..=37).contains(&d), "depth was {d}");
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
