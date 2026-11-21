use std::array::from_fn;

use crate::fmm::harmonic::SolidHarmonic;
use crate::fmm::operators::{multipole_to_multipole, particle_to_multipole};
use crate::shared::{System, Vec3};

const MAX_LEAF_BODIES: usize = 2;
pub struct FmmTree<'a> {
    system: System<'a>,
    mass_distribution: SolidHarmonic, // the mass distribution in the cell as moments
    external_potential: SolidHarmonic, // the potential field felt in the cell
    children: Vec<FmmTree<'a>>,
}
impl<'a> FmmTree<'a> {
    pub fn new(system: System<'a>, p: usize) -> FmmTree<'a> {
        let mass_distribution: SolidHarmonic = SolidHarmonic::irregular(p);
        let external_potential: SolidHarmonic = SolidHarmonic::regular(p);

        if MAX_LEAF_BODIES < system.bodies.len() {
            // split
            let new_system: System<'_> = System::new(system.center, system.hw);
            let mut children: Vec<FmmTree<'a>> = Vec::with_capacity(8);
            let systems: [System; 8] = system.partition();
            for child in systems {
                let child_tree: FmmTree = FmmTree::new(child, p);
                children.push(child_tree);
            }

            return FmmTree {
                system: new_system,
                children,
                mass_distribution,
                external_potential,
            };
        } else {
            // don't split
            let children: Vec<FmmTree<'a>> = Vec::new();

            return FmmTree {
                system,
                children,
                mass_distribution,
                external_potential,
            };
        }
    }
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
    pub fn accumulate_mass_distribution(&mut self) {
        if self.is_leaf() {
            self.mass_distribution = particle_to_multipole(
                self.system.bodies.iter().copied(),
                self.system.center,
                self.mass_distribution.order(),
            );
        } else {
            self.mass_distribution = SolidHarmonic::irregular(self.mass_distribution.order());
            for child in &mut self.children {
                child.accumulate_mass_distribution();
                let shift: Vec3 = from_fn(|i| child.system.center[i] - self.system.center[i]);
                self.mass_distribution += &multipole_to_multipole(&child.mass_distribution, shift);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shared::test_util::*;
    use crate::shared::{Body, ORIGIN, RNG, Vec3};

    // --------------------------------------------------------------- fixtures

    /// `n` bodies uniform in the cube [-1, 1]^3. Unlike a lattice, nothing here
    /// lands on a cell face, so the input contract genuinely holds (checked by
    /// `fixtures_keep_every_body_off_the_cell_faces`).
    fn random_bodies(seed: u64, n: usize) -> Vec<Body> {
        let mut rng = RNG::new(seed);
        (0..n).map(|_| Body::rand(&mut rng, 1.0)).collect()
    }

    /// A tight clump (width 2e-3) plus one body far away. A uniform tree is
    /// shallow; this one must recurse deep and unevenly to separate the clump.
    fn clustered(clump: usize) -> Vec<Body> {
        let mut rng = RNG::new(7);
        let mut bodies: Vec<Body> = (0..clump)
            .map(|_| {
                let b = Body::rand(&mut rng, 1e-3);
                body_at(0.4 + b.pos[0], -0.3 + b.pos[1], 0.2 + b.pos[2])
            })
            .collect();
        bodies.push(body_at(-0.9, 0.8, -0.7));
        bodies
    }

    /// Every fixture the structural tests run over. There is no single-body
    /// fixture on purpose: `System::from_bodies` gives one body a box of
    /// half-width 0, which puts it on its own box's face by definition, so no
    /// geometric claim can hold. `a_single_body_is_one_leaf_holding_it` covers it.
    fn fixtures() -> Vec<(String, Vec<Body>)> {
        let mut out: Vec<(String, Vec<Body>)> = [2, 3, 20, 300, 1000]
            .iter()
            .map(|&n| (format!("random({n})"), random_bodies(100 + n as u64, n)))
            .collect();
        out.push(("clustered(30)".to_string(), clustered(30)));
        out
    }

    fn tree_of(bodies: &[Body], p: usize) -> FmmTree<'_> {
        FmmTree::new(System::from_bodies(bodies), p)
    }

    // ---------------------------------------------------------------- helpers

    /// Every leaf under `tree`, in depth-first order.
    fn leaves<'t, 'a>(tree: &'t FmmTree<'a>) -> Vec<&'t FmmTree<'a>> {
        if tree.children.is_empty() {
            vec![tree]
        } else {
            tree.children.iter().flat_map(|c| leaves(c)).collect()
        }
    }

    /// Every node (internal and leaf) under and including `tree`.
    fn all_nodes<'t, 'a>(tree: &'t FmmTree<'a>) -> Vec<&'t FmmTree<'a>> {
        let mut out = vec![tree];
        for child in &tree.children {
            out.extend(all_nodes(child));
        }
        out
    }

    /// Number of levels below `tree` (a lone leaf has depth 0).
    fn depth(tree: &FmmTree) -> usize {
        tree.children
            .iter()
            .map(|c| 1 + depth(c))
            .max()
            .unwrap_or(0)
    }

    /// The tree's shape alone: (number of children, number of bodies) per node,
    /// depth-first.
    fn shape(tree: &FmmTree) -> Vec<(usize, usize)> {
        all_nodes(tree)
            .iter()
            .map(|n| (n.children.len(), n.system.bodies.len()))
            .collect()
    }

    fn addresses<'b>(bodies: impl Iterator<Item = &'b Body>) -> Vec<*const Body> {
        let mut out: Vec<*const Body> = bodies.map(|b| b as *const Body).collect();
        out.sort();
        out
    }

    /// Strictly inside: no tolerance, so a body on a face does not count.
    fn strictly_inside(center: Vec3, hw: f64, pos: Vec3) -> bool {
        (0..3).all(|d| (pos[d] - center[d]).abs() < hw)
    }

    /// Smallest distance from any leaf body to any face of its own leaf.
    fn smallest_gap_to_a_face(tree: &FmmTree) -> f64 {
        let mut gap = f64::MAX;
        for leaf in leaves(tree) {
            let (c, hw) = (leaf.system.center, leaf.system.hw);
            for b in &leaf.system.bodies {
                for d in 0..3 {
                    gap = gap.min((b.pos[d] - (c[d] - hw)).abs());
                    gap = gap.min((b.pos[d] - (c[d] + hw)).abs());
                }
            }
        }
        gap
    }

    // ------------------------------------------------------------ leaf cases

    #[test]
    fn an_empty_system_is_an_empty_leaf() {
        let tree = FmmTree::new(System::new(ORIGIN, 1.0), 2);
        assert!(tree.children.is_empty());
        assert!(tree.system.bodies.is_empty());
    }

    #[test]
    fn a_single_body_is_one_leaf_holding_it() {
        let bodies = random_bodies(1, 1);
        let tree = tree_of(&bodies, 2);
        assert!(tree.children.is_empty());
        assert_eq!(tree.system.bodies.len(), 1);
    }

    #[test]
    fn a_system_at_the_limit_is_one_leaf_holding_every_body() {
        let bodies = random_bodies(2, MAX_LEAF_BODIES);
        let tree = tree_of(&bodies, 2);
        assert!(tree.children.is_empty(), "exactly the limit must not split");
        assert_eq!(tree.system.bodies.len(), MAX_LEAF_BODIES);
    }

    #[test]
    fn one_body_over_the_limit_splits() {
        let bodies = random_bodies(3, MAX_LEAF_BODIES + 1);
        let tree = tree_of(&bodies, 2);
        assert_eq!(tree.children.len(), 8);
    }

    // ------------------------------------------------------- tree structure

    /// The structural tests below are only about the contract if it holds: this
    /// is what stops them quietly passing on a boundary case through tolerance.
    /// The bar separates "on a face" from "merely near one": the old lattice
    /// fixture had a body 1e-16 from a face (rounding noise), while random
    /// bodies come within about 1e-7 of one now and then, which is harmless.
    #[test]
    fn fixtures_keep_every_body_off_the_cell_faces() {
        for (name, bodies) in fixtures() {
            let gap = smallest_gap_to_a_face(&tree_of(&bodies, 2));
            assert!(gap > 1e-9, "{name}: a body is {gap:e} from a cell face");
        }
    }

    #[test]
    fn every_internal_node_has_eight_children_and_no_bodies() {
        let mut internal_seen = 0;
        for (name, bodies) in fixtures() {
            let tree = tree_of(&bodies, 2);
            for node in all_nodes(&tree)
                .into_iter()
                .filter(|n| !n.children.is_empty())
            {
                internal_seen += 1;
                assert_eq!(node.children.len(), 8, "{name}");
                assert!(
                    node.system.bodies.is_empty(),
                    "{name}: internal nodes are bodyless"
                );
            }
        }
        assert!(
            internal_seen > 100,
            "vacuous: only {internal_seen} internal nodes"
        );
    }

    #[test]
    fn no_leaf_holds_more_than_the_limit_and_some_hold_exactly_the_limit() {
        let mut biggest = 0;
        for (name, bodies) in fixtures() {
            for leaf in leaves(&tree_of(&bodies, 2)) {
                let n = leaf.system.bodies.len();
                assert!(n <= MAX_LEAF_BODIES, "{name}: a leaf holds {n}");
                biggest = biggest.max(n);
            }
        }
        assert_eq!(
            biggest, MAX_LEAF_BODIES,
            "the limit is never reached: a vacuous check"
        );
    }

    #[test]
    fn the_leaves_hold_every_body_exactly_once() {
        for (name, bodies) in fixtures() {
            let tree = tree_of(&bodies, 2);
            let got = addresses(
                leaves(&tree)
                    .iter()
                    .flat_map(|l| l.system.bodies.iter().copied()),
            );
            let want = addresses(bodies.iter());
            assert_eq!(got, want, "{name}: bodies were lost or duplicated");
        }
    }

    #[test]
    fn every_body_lies_strictly_inside_its_leaf() {
        for (name, bodies) in fixtures() {
            let tree = tree_of(&bodies, 2);
            for leaf in leaves(&tree) {
                for b in &leaf.system.bodies {
                    assert!(
                        strictly_inside(leaf.system.center, leaf.system.hw, b.pos),
                        "{name}"
                    );
                }
            }
        }
    }

    /// Child `i` is the octant whose bit 0 / 1 / 2 is set when x / y / z is on
    /// the high side of the parent's center (the order `partition` documents).
    #[test]
    fn child_i_is_the_octant_with_bit_pattern_i() {
        let bodies = random_bodies(5, 300);
        let tree = tree_of(&bodies, 2);
        for node in all_nodes(&tree)
            .into_iter()
            .filter(|n| !n.children.is_empty())
        {
            for (i, child) in node.children.iter().enumerate() {
                assert!(close(child.system.hw, node.system.hw / 2.0));
                for d in 0..3 {
                    let offset = child.system.center[d] - node.system.center[d];
                    assert!(close(offset.abs(), node.system.hw / 2.0), "axis {d}");
                    let high = (i >> d) & 1 == 1;
                    assert_eq!(offset > 0.0, high, "child {i}, axis {d}");
                }
            }
        }
    }

    #[test]
    fn the_root_keeps_the_region_it_was_built_from() {
        let bodies = random_bodies(6, 50);
        let system = System::from_bodies(&bodies);
        let (center, hw) = (system.center, system.hw);
        let tree = FmmTree::new(system, 2);
        assert!(close3(tree.system.center, center));
        assert!(close(tree.system.hw, hw));
    }

    #[test]
    fn the_shape_of_the_tree_does_not_depend_on_the_order() {
        for (name, bodies) in fixtures() {
            let reference = shape(&tree_of(&bodies, 0));
            for p in [1, 2, 5] {
                assert_eq!(shape(&tree_of(&bodies, p)), reference, "{name}, p={p}");
            }
        }
    }

    #[test]
    fn a_clump_makes_a_deep_uneven_tree() {
        let clump = clustered(30);
        let uniform = random_bodies(8, clump.len());
        let (deep, shallow) = (depth(&tree_of(&clump, 2)), depth(&tree_of(&uniform, 2)));
        // separating bodies 1e-3 apart inside a box of half-width ~1 takes about
        // ten halvings; a uniform cloud of the same size separates in a few
        assert!(deep >= 8, "clump depth was only {deep}");
        assert!(deep > shallow + 3, "clump {deep} vs uniform {shallow}");
    }

    // ---------------------------------------------------------- coefficients

    /// What a `SolidHarmonic` is lives in `harmonic.rs`. The tree's job is to give
    /// every node a fresh, empty one of the right order and kind, so compare with
    /// a freshly built one: that pins the kind, the length and the zeros at once.
    /// mass_distribution (the cell's own mass, valid outside the cell) is irregular;
    /// external_potential (the potential from well-separated cells, valid inside the
    /// cell) is regular. This checks the code against that convention, not the
    /// convention against the physics.
    #[test]
    fn every_node_starts_with_an_empty_irregular_mass_distribution_of_order_p() {
        for p in [0, 1, 3] {
            let bodies = random_bodies(9, 60);
            let tree = tree_of(&bodies, p);
            for node in all_nodes(&tree) {
                assert_eq!(node.mass_distribution, SolidHarmonic::irregular(p), "p={p}");
            }
        }
    }

    #[test]
    fn every_node_starts_with_an_empty_regular_external_potential_of_order_p() {
        for p in [0, 1, 3] {
            let bodies = random_bodies(9, 60);
            let tree = tree_of(&bodies, p);
            for node in all_nodes(&tree) {
                assert_eq!(node.external_potential, SolidHarmonic::regular(p), "p={p}");
            }
        }
    }

    // -------------------------------------------- accumulating the leaf summaries

    /// Every leaf, mutably, so the accumulation can be run on each one by itself.
    fn leaves_mut<'t, 'a>(tree: &'t mut FmmTree<'a>) -> Vec<&'t mut FmmTree<'a>> {
        if tree.children.is_empty() {
            vec![tree]
        } else {
            tree.children
                .iter_mut()
                .flat_map(|child| leaves_mut(child))
                .collect()
        }
    }

    /// Slot 0 of a summary is the total mass, so each leaf must hold exactly the
    /// mass of its own bodies and the leaves together must hold all of it. The
    /// expected values come from summing masses, not from the operator.
    #[test]
    fn each_leaf_accumulates_the_mass_of_its_own_bodies() {
        for (name, bodies) in fixtures() {
            let total: f64 = bodies.iter().map(|b| b.mass).sum();
            let mut tree = tree_of(&bodies, 3);
            let mut sum_of_leaves = 0.0;
            let mut crowded = 0;
            for leaf in leaves_mut(&mut tree) {
                leaf.accumulate_mass_distribution();
                let own: f64 = leaf.system.bodies.iter().map(|b| b.mass).sum();
                let c00 = leaf.mass_distribution.get(0, 0);
                assert!(
                    (c00.re - own).abs() <= 1e-12 * own.max(1.0) && c00.im == 0.0,
                    "{name}: leaf holds {} bodies of total mass {own}, summary says {c00:?}",
                    leaf.system.bodies.len()
                );
                sum_of_leaves += c00.re;
                if leaf.system.bodies.len() >= 2 {
                    crowded += 1;
                }
            }
            assert!((sum_of_leaves - total).abs() <= 1e-9 * total, "{name}");
            assert!(crowded > 0, "{name}: no leaf holds several bodies");
        }
    }

    /// The summary of a leaf, evaluated well outside the leaf's own cell, must
    /// reproduce the direct potential of that leaf's bodies. This is what fails if
    /// the wrong center, bodies or order are passed to the operator.
    #[test]
    fn a_leaf_summary_reproduces_the_potential_of_its_bodies_outside_its_cell() {
        let dir = [0.3, -0.5, 0.8];
        let norm = (0.3f64.powi(2) + 0.5f64.powi(2) + 0.8f64.powi(2)).sqrt();
        for (name, bodies) in fixtures() {
            let mut tree = tree_of(&bodies, 6);
            let mut checked = 0;
            for leaf in leaves_mut(&mut tree) {
                leaf.accumulate_mass_distribution();
                if leaf.system.bodies.is_empty() {
                    continue;
                }
                let (center, hw) = (leaf.system.center, leaf.system.hw);
                // six half-widths away: the farthest body is at most sqrt(3) hw out
                let rel: Vec3 = std::array::from_fn(|i| 6.0 * hw * dir[i] / norm);
                let x: Vec3 = std::array::from_fn(|i| center[i] + rel[i]);
                let direct: f64 = leaf
                    .system
                    .bodies
                    .iter()
                    .map(|b| b.mass / b.distance(&x))
                    .sum();
                let got = leaf.mass_distribution.value(rel);
                assert!(
                    (got.re - direct).abs() <= 1e-3 * direct && got.im.abs() <= 1e-3 * direct,
                    "{name}: summary gives {got:?}, direct sum is {direct}"
                );
                if leaf.system.bodies.len() >= 2 {
                    checked += 1;
                }
            }
            assert!(
                checked > 0,
                "{name}: no leaf with several bodies was checked"
            );
        }
    }

    /// Accumulating must fill only the mass distribution: the external potential
    /// stays empty, the kinds stay as constructed and the order is unchanged.
    #[test]
    fn accumulating_touches_only_the_mass_distribution_and_keeps_its_kind_and_order() {
        let p = 3;
        let bodies = random_bodies(5, 40);
        let mut tree = tree_of(&bodies, p);
        for leaf in leaves_mut(&mut tree) {
            leaf.accumulate_mass_distribution();
            assert!(matches!(
                leaf.mass_distribution,
                SolidHarmonic::Irregular(_)
            ));
            assert_eq!(leaf.mass_distribution.order(), p);
            assert_eq!(leaf.external_potential, SolidHarmonic::regular(p));
        }
    }

    /// A single body gets a cell of half-width 0 centered on itself, so it sits
    /// exactly at its own center, the case with no direction.
    #[test]
    fn a_single_body_root_accumulates_a_finite_summary_holding_its_mass() {
        let bodies = vec![body_with(2.5, 0.3, -0.2, 0.5)];
        let mut tree = tree_of(&bodies, 3);
        assert!(tree.is_leaf());
        tree.accumulate_mass_distribution();
        for (k, c) in tree.mass_distribution.coeffs().iter().enumerate() {
            let want = if k == 0 { 2.5 } else { 0.0 };
            assert!(
                c.re.is_finite() && c.im.is_finite(),
                "slot {k} is not finite"
            );
            assert_eq!((c.re, c.im), (want, 0.0), "slot {k}");
        }
    }

    // ------------------------------------------- accumulating up the whole tree

    /// The bodies in every leaf beneath `node`.
    fn bodies_beneath<'t>(node: &'t FmmTree) -> Vec<&'t Body> {
        leaves(node)
            .iter()
            .flat_map(|leaf| leaf.system.bodies.iter().copied())
            .collect()
    }

    /// Slot 0 is the total mass and shifting never changes it, so after the whole
    /// upward pass every node, at every level, must hold the mass of exactly the
    /// bodies beneath it. A forgotten or doubled child shows up here.
    #[test]
    fn every_node_holds_the_mass_of_all_the_bodies_beneath_it() {
        for (name, bodies) in fixtures() {
            let mut tree = tree_of(&bodies, 3);
            tree.accumulate_mass_distribution();
            let mut deepest_internal = 0;
            for node in all_nodes(&tree) {
                let beneath: f64 = bodies_beneath(node).iter().map(|b| b.mass).sum();
                let c00 = node.mass_distribution.get(0, 0);
                assert!(
                    (c00.re - beneath).abs() <= 1e-12 * beneath.max(1.0) && c00.im == 0.0,
                    "{name}: node over {} bodies holds mass {c00:?}, expected {beneath}",
                    bodies_beneath(node).len()
                );
                if !node.is_leaf() {
                    deepest_internal = deepest_internal.max(depth(node));
                }
            }
            if bodies.len() >= 300 {
                assert!(
                    deepest_internal >= 2,
                    "{name}: tree too shallow to test several levels"
                );
            }
        }
    }

    /// Each node's summary, read well outside that node's own cell, must give the
    /// potential of all the bodies beneath it. Internal nodes got theirs only
    /// through the shifts, so this fails if a shift has the wrong sign or size.
    #[test]
    fn every_node_summary_reproduces_the_potential_of_the_bodies_beneath_it() {
        let dir = [0.3, -0.5, 0.8];
        let norm = (0.3f64.powi(2) + 0.5f64.powi(2) + 0.8f64.powi(2)).sqrt();
        for (name, bodies) in fixtures() {
            let mut tree = tree_of(&bodies, 8);
            tree.accumulate_mass_distribution();
            let mut internal_checked = 0;
            for node in all_nodes(&tree) {
                let beneath = bodies_beneath(node);
                if beneath.is_empty() {
                    continue;
                }
                let (center, hw) = (node.system.center, node.system.hw);
                // eight half-widths out: every body is within sqrt(3) hw of the center
                let rel: Vec3 = from_fn(|i| 8.0 * hw * dir[i] / norm);
                let x: Vec3 = from_fn(|i| center[i] + rel[i]);
                let direct: f64 = beneath.iter().map(|b| b.mass / b.distance(&x)).sum();
                let got = node.mass_distribution.value(rel);
                assert!(
                    (got.re - direct).abs() <= 1e-4 * direct && got.im.abs() <= 1e-4 * direct,
                    "{name}: node over {} bodies gives {got:?}, direct sum is {direct}",
                    beneath.len()
                );
                if !node.is_leaf() {
                    internal_checked += 1;
                }
            }
            if bodies.len() > MAX_LEAF_BODIES {
                assert!(internal_checked > 0, "{name}: no internal node was checked");
            }
        }
    }
}
