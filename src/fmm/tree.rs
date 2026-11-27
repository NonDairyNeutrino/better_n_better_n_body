use std::array::from_fn;
use std::cell::RefCell;

use crate::fmm::harmonic::SolidHarmonic;
use crate::fmm::operators::m2l::M2LScalings;
use crate::fmm::operators::{
    local_to_local, local_to_particle, multipole_to_local, multipole_to_multipole,
    particle_to_multipole, particle_to_particle,
};
use crate::shared::{Body, ORIGIN, System, Vec3};

const MAX_LEAF_BODIES: usize = 64;
pub struct FmmTree<'a> {
    system: System<'a>,
    mass_distribution: SolidHarmonic, // the mass distribution in the cell as moments
    external_potential: RefCell<SolidHarmonic>, // the potential field felt in the cell
    children: Vec<FmmTree<'a>>,
    acceleration: RefCell<Vec<Vec3>>,
}
impl<'a> FmmTree<'a> {
    pub fn new(system: System<'a>, p: usize) -> FmmTree<'a> {
        FmmTree::with_leaf_limit(system, p, MAX_LEAF_BODIES)
    }
    /// Like `new`, but a cell is split only when it holds more than `leaf_limit` bodies.
    pub fn with_leaf_limit(system: System<'a>, p: usize, leaf_limit: usize) -> FmmTree<'a> {
        let mass_distribution: SolidHarmonic = SolidHarmonic::irregular(p);
        let external_potential: RefCell<SolidHarmonic> = RefCell::new(SolidHarmonic::regular(p));

        if leaf_limit < system.bodies.len() {
            // split
            let new_system: System<'_> = System::new(system.center, system.hw);
            let mut children: Vec<FmmTree<'a>> = Vec::with_capacity(8);
            let systems: [System; 8] = system.partition();
            for child in systems {
                let child_tree: FmmTree = FmmTree::with_leaf_limit(child, p, leaf_limit);
                children.push(child_tree);
            }
            let acceleration: RefCell<Vec<[f64; 3]>> = RefCell::new(Vec::new());

            return FmmTree {
                system: new_system,
                children,
                mass_distribution,
                external_potential,
                acceleration,
            };
        } else {
            // don't split
            let children: Vec<FmmTree<'a>> = Vec::new();
            let acceleration: RefCell<Vec<[f64; 3]>> =
                RefCell::new(vec![[0.0; 3]; system.bodies.len()]);

            return FmmTree {
                system,
                children,
                mass_distribution,
                external_potential,
                acceleration,
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

    pub fn assemble_potential(&self) {
        let scale = M2LScalings::new(self.mass_distribution.order());
        self.interact(self, &scale);
    }
    fn interact(&self, source: &FmmTree, scale: &M2LScalings) {
        if well_separated(&self.system, &source.system) {
            let shift = from_fn(|i| self.system.center[i] - source.system.center[i]);
            *self.external_potential.borrow_mut() +=
                &multipole_to_local(&source.mass_distribution, shift, scale);
        } else if self.is_leaf() && source.is_leaf() {
            // if the cells are close enough, just add the bodies directly
            particle_to_particle(
                &self.system.bodies,
                &source.system.bodies,
                &mut self.acceleration.borrow_mut(),
            );
        } else if self.is_leaf() {
            for child in source.children.iter() {
                self.interact(child, scale);
            }
        } else if source.is_leaf() {
            for child in self.children.iter() {
                child.interact(source, scale);
            }
        } else if self.system.hw >= source.system.hw {
            // split the larger one
            for child in self.children.iter() {
                child.interact(source, scale);
            }
        } else {
            for child in source.children.iter() {
                self.interact(child, scale);
            }
        }
    }
    pub fn inherit_potential(&self) {
        if self.is_leaf() {
            local_to_particle(
                &self.external_potential.borrow(),
                self.system.center,
                &self.system.bodies,
                &mut self.acceleration.borrow_mut(),
            );
        } else {
            for child in self.children.iter() {
                let shift: Vec3 = from_fn(|i| child.system.center[i] - self.system.center[i]);
                *child.external_potential.borrow_mut() +=
                    &local_to_local(&self.external_potential.borrow(), shift);
                child.inherit_potential();
            }
        }
    }
    pub fn accelerations(&self, bodies: &[Body]) -> Vec<Vec3> {
        let mut out: Vec<Vec3> = vec![ORIGIN; bodies.len()];
        self.collect_accelerations(bodies, &mut out);
        return out;
    }
    fn collect_accelerations(&self, bodies: &[Body], out: &mut [Vec3]) {
        if self.is_leaf() {
            for (body, acc) in self
                .system
                .bodies
                .iter()
                .zip(self.acceleration.borrow().iter())
            {
                let slot = index_in(*body, bodies);
                out[slot].iter_mut().zip(acc).for_each(|(e, a)| *e += a);
            }
        } else {
            for child in self.children.iter() {
                child.collect_accelerations(bodies, out);
            }
        }
    }
}

const THETA: f64 = 0.5;
fn well_separated(target: &System, source: &System) -> bool {
    let dist2: f64 = target
        .center
        .iter()
        .zip(source.center)
        .map(|(t, s)| (t - s).powi(2))
        .sum::<f64>();

    3.0_f64 * (target.hw + source.hw).powi(2) < THETA.powi(2) * dist2
}

fn index_in(body: &Body, bodies: &[Body]) -> usize {
    let body_addr: usize = std::ptr::from_ref(body).addr();
    let start_addr: usize = bodies.as_ptr().addr();
    let offset: usize = body_addr
        .checked_sub(start_addr)
        .expect("Body comes before slice");
    let index = offset / size_of::<Body>();

    assert!(index < bodies.len());
    return index;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fmm::test_util::{assert_close, pull_from};
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

    /// The leaf limit the structural tests build with. It is small so that the
    /// small fixtures grow trees deep enough to test something; the production
    /// limit is covered by `the_default_constructor_uses_the_production_limit`
    /// and by the tests at the end of this module.
    const TEST_LEAF: usize = 2;

    fn tree_of(bodies: &[Body], p: usize) -> FmmTree<'_> {
        FmmTree::with_leaf_limit(System::from_bodies(bodies), p, TEST_LEAF)
    }

    /// A tree built the way the program builds it.
    fn production_tree_of(bodies: &[Body], p: usize) -> FmmTree<'_> {
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
        let tree = FmmTree::with_leaf_limit(System::new(ORIGIN, 1.0), 2, TEST_LEAF);
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
        let bodies = random_bodies(2, TEST_LEAF);
        let tree = tree_of(&bodies, 2);
        assert!(tree.children.is_empty(), "exactly the limit must not split");
        assert_eq!(tree.system.bodies.len(), TEST_LEAF);
    }

    #[test]
    fn one_body_over_the_limit_splits() {
        let bodies = random_bodies(3, TEST_LEAF + 1);
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
                assert!(n <= TEST_LEAF, "{name}: a leaf holds {n}");
                biggest = biggest.max(n);
            }
        }
        assert_eq!(
            biggest, TEST_LEAF,
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
        let tree = FmmTree::with_leaf_limit(system, 2, TEST_LEAF);
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
                assert_eq!(
                    *node.external_potential.borrow(),
                    SolidHarmonic::regular(p),
                    "p={p}"
                );
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
            assert_eq!(*leaf.external_potential.borrow(), SolidHarmonic::regular(p));
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
            if bodies.len() > TEST_LEAF {
                assert!(internal_checked > 0, "{name}: no internal node was checked");
            }
        }
    }

    /// Accumulating is a pure function of the bodies, so doing it again must
    /// leave every node exactly as it was. Internal nodes sum their children into
    /// whatever they already hold, so without a reset they would double.
    #[test]
    fn accumulating_twice_gives_the_same_result_as_accumulating_once() {
        for (name, bodies) in fixtures() {
            let mut tree = tree_of(&bodies, 4);
            tree.accumulate_mass_distribution();
            let once: Vec<SolidHarmonic> = all_nodes(&tree)
                .iter()
                .map(|n| n.mass_distribution.clone())
                .collect();
            tree.accumulate_mass_distribution();
            let twice: Vec<SolidHarmonic> = all_nodes(&tree)
                .iter()
                .map(|n| n.mass_distribution.clone())
                .collect();
            assert_eq!(once, twice, "{name}");
            if bodies.len() > TEST_LEAF {
                let internal_with_mass = all_nodes(&tree)
                    .iter()
                    .filter(|n| !n.is_leaf() && n.mass_distribution.get(0, 0).re > 0.0)
                    .count();
                assert!(internal_with_mass > 0, "{name}: nothing internal to test");
            }
        }
    }

    // ----------------------------------------------------------- interacting

    /// Fixtures small enough to convert every far pair at a realistic order in
    /// a debug build: a few bodies, a random cloud, and the deep uneven clump.
    fn interaction_fixtures() -> Vec<(String, Vec<Body>)> {
        vec![
            ("random(20)".to_string(), random_bodies(120, 20)),
            ("random(100)".to_string(), random_bodies(300, 100)),
            ("clustered(30)".to_string(), clustered(30)),
        ]
    }

    /// The pairs `interact` must produce, found again from the same rules with a
    /// work list instead of recursion: (far pairs, near pairs) as (target, source).
    /// A well-separated pair is converted (far); two leaves that are not well
    /// separated are summed directly (near); otherwise the larger cell is split,
    /// the target on a tie.
    #[allow(clippy::type_complexity)]
    fn expected_pairs<'t, 'a>(
        root: &'t FmmTree<'a>,
    ) -> (
        Vec<(&'t FmmTree<'a>, &'t FmmTree<'a>)>,
        Vec<(&'t FmmTree<'a>, &'t FmmTree<'a>)>,
    ) {
        let (mut far, mut near) = (vec![], vec![]);
        let mut work = vec![(root, root)];
        while let Some((t, s)) = work.pop() {
            if well_separated(&t.system, &s.system) {
                far.push((t, s));
            } else if t.is_leaf() && s.is_leaf() {
                near.push((t, s));
            } else {
                let split_target = if t.is_leaf() {
                    false
                } else if s.is_leaf() {
                    true
                } else {
                    t.system.hw >= s.system.hw
                };
                if split_target {
                    work.extend(t.children.iter().map(|c| (c, s)));
                } else {
                    work.extend(s.children.iter().map(|c| (t, c)));
                }
            }
        }
        (far, near)
    }

    /// The nodes from the root down to `target`, both included.
    fn chain_to<'t, 'a>(root: &'t FmmTree<'a>, target: &FmmTree) -> Option<Vec<&'t FmmTree<'a>>> {
        if std::ptr::eq(root, target) {
            return Some(vec![root]);
        }
        for child in &root.children {
            if let Some(mut chain) = chain_to(child, target) {
                chain.insert(0, root);
                return Some(chain);
            }
        }
        None
    }

    /// For every leaf and every body, the number of times that body is accounted
    /// for: through a far pair at one of the leaf's ancestors, or through a near
    /// pair at the leaf itself. It must be exactly once, which is what makes the
    /// split of the potential into near and far exact.
    #[test]
    fn every_leaf_accounts_for_every_body_exactly_once() {
        for (name, bodies) in interaction_fixtures() {
            let tree = tree_of(&bodies, 2);
            let (far, near) = expected_pairs(&tree);
            let index = |b: &Body| {
                (b as *const Body as usize - bodies.as_ptr() as usize) / std::mem::size_of::<Body>()
            };
            for leaf in leaves(&tree) {
                let chain = chain_to(&tree, leaf).unwrap();
                let mut seen = vec![0usize; bodies.len()];
                for (t, s) in &far {
                    if chain.iter().any(|a| std::ptr::eq(*a, *t)) {
                        for b in bodies_beneath(s) {
                            seen[index(b)] += 1;
                        }
                    }
                }
                for (t, s) in &near {
                    if std::ptr::eq(*t, leaf) {
                        for b in bodies_beneath(s) {
                            seen[index(b)] += 1;
                        }
                    }
                }
                assert!(
                    seen.iter().all(|&c| c == 1),
                    "{name}: a leaf of {} bodies sees some body {:?} times",
                    leaf.system.bodies.len(),
                    seen.iter().find(|&&c| c != 1)
                );
            }
            assert!(!near.is_empty(), "{name}: no near pairs at all");
            if bodies.len() >= 100 {
                assert!(!far.is_empty(), "{name}: no far pairs, nothing is tested");
            }
        }
    }

    /// Each node's external potential must be exactly the conversions of the far
    /// partners the rules give it, one per partner and none for any other.
    #[test]
    fn interacting_gives_each_node_the_conversions_of_its_far_partners() {
        for (name, bodies) in interaction_fixtures() {
            let p = 4;
            let mut tree = tree_of(&bodies, p);
            tree.accumulate_mass_distribution();
            tree.interact(&tree, &M2LScalings::new(p));
            let (far, _) = expected_pairs(&tree);
            let mut nonzero_nodes = 0;
            for node in all_nodes(&tree) {
                let mut want = SolidHarmonic::regular(p);
                for (t, s) in &far {
                    if std::ptr::eq(*t, node) {
                        let shift: Vec3 = from_fn(|i| t.system.center[i] - s.system.center[i]);
                        want += &multipole_to_local(&s.mass_distribution, shift, &M2LScalings::new(p));
                    }
                }
                let got = node.external_potential.borrow();
                for (k, (g, w)) in got.coeffs().iter().zip(want.coeffs()).enumerate() {
                    assert_close(
                        *g,
                        *w,
                        1e-12 * (1.0 + w.re.abs() + w.im.abs()),
                        &format!(
                            "{name}: node of {} bodies, slot {k}",
                            bodies_beneath(node).len()
                        ),
                    );
                }
                if want.coeffs().iter().any(|c| c.re != 0.0 || c.im != 0.0) {
                    nonzero_nodes += 1;
                }
            }
            if bodies.len() >= 100 {
                assert!(
                    nonzero_nodes > 0,
                    "{name}: every external potential is zero"
                );
            }
        }
    }

    /// The physics of the whole stage. At a point inside a leaf, the external
    /// potentials of the leaf and all its ancestors, plus the bodies of the leaf's
    /// near partners summed directly, must give the potential of all the bodies.
    #[test]
    fn external_potentials_plus_near_bodies_give_the_potential_of_all_bodies() {
        for (name, bodies) in interaction_fixtures() {
            if bodies.len() < 20 {
                continue;
            }
            let p = 6;
            let mut tree = tree_of(&bodies, p);
            tree.accumulate_mass_distribution();
            tree.interact(&tree, &M2LScalings::new(p));
            let (_, near) = expected_pairs(&tree);
            let mut worst: f64 = 0.0;
            let mut leaves_with_far_part = 0;
            for leaf in leaves(&tree) {
                if leaf.system.bodies.is_empty() {
                    continue;
                }
                let chain = chain_to(&tree, leaf).unwrap();
                let x: Vec3 =
                    from_fn(|i| leaf.system.center[i] + [0.3, -0.2, 0.4][i] * leaf.system.hw);
                let far_part: f64 = chain
                    .iter()
                    .map(|a| {
                        let rel: Vec3 = from_fn(|i| x[i] - a.system.center[i]);
                        a.external_potential.borrow().value(rel).re
                    })
                    .sum();
                let near_part: f64 = near
                    .iter()
                    .filter(|(t, _)| std::ptr::eq(*t, leaf))
                    .flat_map(|(_, s)| s.system.bodies.iter())
                    .map(|b| b.mass / b.distance(&x))
                    .sum();
                let all: f64 = bodies.iter().map(|b| b.mass / b.distance(&x)).sum();
                worst = worst.max((far_part + near_part - all).abs() / all);
                if far_part.abs() > 0.01 * all {
                    leaves_with_far_part += 1;
                }
            }
            assert!(
                leaves_with_far_part > 0,
                "{name}: no leaf got a real far part"
            );
        }
    }

    /// After interacting, each leaf holds, for each of its bodies, the pull of
    /// every body in the leaf's near partners (itself left out), with Newton's
    /// law written out. Internal nodes hold no accelerations at all.
    #[test]
    fn interacting_gives_each_leaf_the_pull_of_its_near_partners() {
        for (name, bodies) in interaction_fixtures() {
            let mut tree = tree_of(&bodies, 2);
            tree.accumulate_mass_distribution();
            tree.interact(&tree, &M2LScalings::new(2));
            let (_, near) = expected_pairs(&tree);
            let mut nonzero_bodies = 0;
            for node in all_nodes(&tree) {
                let got = node.acceleration.borrow();
                if !node.is_leaf() {
                    assert!(
                        got.is_empty(),
                        "{name}: an internal node holds accelerations"
                    );
                    continue;
                }
                assert_eq!(got.len(), node.system.bodies.len(), "{name}");
                for (i, body) in node.system.bodies.iter().enumerate() {
                    let mut want = [0.0; 3];
                    for (t, s) in &near {
                        if std::ptr::eq(*t, node) {
                            for other in s.system.bodies.iter() {
                                if !std::ptr::eq(*other, *body) {
                                    let pull = pull_from(body.pos, other);
                                    for d in 0..3 {
                                        want[d] += pull[d];
                                    }
                                }
                            }
                        }
                    }
                    for d in 0..3 {
                        assert!(
                            (got[i][d] - want[d]).abs() <= 1e-12 * (1.0 + want[d].abs()),
                            "{name}: body {i} component {d}: got {}, want {}",
                            got[i][d],
                            want[d]
                        );
                    }
                    if want.iter().any(|c| c.abs() > 0.0) {
                        nonzero_bodies += 1;
                    }
                }
            }
            assert!(nonzero_bodies > 0, "{name}: no body has a near pull");
        }
    }

    #[test]
    fn the_default_constructor_uses_the_production_limit() {
        let at_limit = random_bodies(20, MAX_LEAF_BODIES);
        assert!(production_tree_of(&at_limit, 2).is_leaf());
        let over = random_bodies(21, MAX_LEAF_BODIES + 1);
        assert!(!production_tree_of(&over, 2).is_leaf());
    }

    #[test]
    fn a_larger_limit_gives_a_shallower_tree_and_a_smaller_one_a_deeper_tree() {
        let bodies = random_bodies(22, 300);
        let at = |limit: usize| {
            depth(&FmmTree::with_leaf_limit(
                System::from_bodies(&bodies),
                1,
                limit,
            ))
        };
        assert!(at(2) > at(20));
        assert!(at(20) > at(150));
        assert_eq!(at(300), 0, "everything fits in one leaf");
    }

    // ============================================================ index_in

    #[test]
    fn index_in_finds_every_bodys_own_slot() {
        let bodies = random_bodies(1, 40);
        for k in 0..bodies.len() {
            assert_eq!(index_in(&bodies[k], &bodies), k, "body {k}");
        }
    }

    #[test]
    fn index_in_counts_from_the_start_of_the_slice_it_is_given() {
        let bodies = random_bodies(2, 10);
        assert_eq!(index_in(&bodies[7], &bodies[5..]), 2);
        assert_eq!(index_in(&bodies[5], &bodies[5..]), 0);
        assert_eq!(index_in(&bodies[9], &bodies[5..]), 4);
    }

    #[test]
    #[should_panic(expected = "before slice")]
    fn index_in_refuses_a_body_before_the_slice() {
        let bodies = random_bodies(3, 10);
        index_in(&bodies[2], &bodies[5..]);
    }

    #[test]
    #[should_panic(expected = "index < bodies.len()")]
    fn index_in_refuses_a_body_past_the_end_of_the_slice() {
        let bodies = random_bodies(4, 10);
        index_in(&bodies[9], &bodies[..5]);
    }

    // ======================================================= well_separated

    // Cells of half-width a and b count as far apart when
    //   3 (a + b)^2 < THETA^2 * distance^2,
    // that is, when the distance exceeds sqrt(3) (a + b) / THETA. With THETA = 0.5
    // that is 6.93 for two unit cells and 5.20 for half-widths 1 and 0.5.

    fn cell(x: f64, y: f64, z: f64, hw: f64) -> System<'static> {
        System::new([x, y, z], hw)
    }

    #[test]
    fn two_unit_cells_are_far_just_beyond_the_threshold_and_near_just_inside_it() {
        let a = cell(0.0, 0.0, 0.0, 1.0);
        assert!(!well_separated(&a, &cell(6.9, 0.0, 0.0, 1.0)));
        assert!(well_separated(&a, &cell(7.0, 0.0, 0.0, 1.0)));
    }

    #[test]
    fn the_threshold_scales_with_the_sum_of_the_half_widths_not_their_squares() {
        let a = cell(0.0, 0.0, 0.0, 1.0);
        // threshold sqrt(3) * 1.5 / 0.5 = 5.196
        assert!(!well_separated(&a, &cell(5.1, 0.0, 0.0, 0.5)));
        assert!(well_separated(&a, &cell(5.3, 0.0, 0.0, 0.5)));
    }

    #[test]
    fn only_the_distance_between_centers_matters_not_the_direction() {
        let a = cell(1.0, -2.0, 0.5, 1.0);
        let along_x = cell(1.0 + 7.0, -2.0, 0.5, 1.0);
        let diagonal = cell(
            1.0 + 7.0 / 3.0_f64.sqrt(),
            -2.0 + 7.0 / 3.0_f64.sqrt(),
            0.5 + 7.0 / 3.0_f64.sqrt(),
            1.0,
        );
        assert!(well_separated(&a, &along_x));
        assert!(well_separated(&a, &diagonal));
        let near_diagonal = cell(
            1.0 + 6.9 / 3.0_f64.sqrt(),
            -2.0 + 6.9 / 3.0_f64.sqrt(),
            0.5 + 6.9 / 3.0_f64.sqrt(),
            1.0,
        );
        assert!(!well_separated(&a, &near_diagonal));
    }

    #[test]
    fn well_separated_is_symmetric() {
        let a = cell(0.0, 0.0, 0.0, 1.0);
        for (x, hw) in [(5.1, 0.5), (5.3, 0.5), (6.9, 1.0), (7.0, 1.0), (3.0, 0.1)] {
            let b = cell(x, 0.0, 0.0, hw);
            assert_eq!(well_separated(&a, &b), well_separated(&b, &a), "x = {x}");
        }
    }

    #[test]
    fn a_cell_is_never_well_separated_from_itself() {
        let a = cell(0.3, 0.2, -0.1, 0.5);
        assert!(!well_separated(&a, &cell(0.3, 0.2, -0.1, 0.5)));
    }

    // ===================================================== accelerations

    #[test]
    fn accelerations_puts_each_leaf_slot_at_its_bodys_original_index() {
        // Tag each leaf slot with the body's original index, then read back.
        // No physics is involved, so this only tests the bookkeeping.
        let bodies = random_bodies(5, 600);
        let tree = production_tree_of(&bodies, 1);
        assert!(
            leaves(&tree).len() > 1,
            "need several leaves to mean anything"
        );

        for leaf in leaves(&tree) {
            let mut slots = leaf.acceleration.borrow_mut();
            for (body, slot) in leaf.system.bodies.iter().zip(slots.iter_mut()) {
                let i = index_in(body, &bodies) as f64;
                *slot = [i, 2.0 * i, -i];
            }
        }

        let got = tree.accelerations(&bodies);
        assert_eq!(got.len(), bodies.len());
        for (i, acc) in got.iter().enumerate() {
            let i = i as f64;
            assert_eq!(*acc, [i, 2.0 * i, -i], "body {i}");
        }
    }

    #[test]
    fn accelerations_of_an_untouched_tree_are_all_zero() {
        let bodies = random_bodies(6, 300);
        let tree = production_tree_of(&bodies, 1);
        let got = tree.accelerations(&bodies);
        assert_eq!(got, vec![[0.0; 3]; bodies.len()]);
    }

    #[test]
    fn accelerations_does_not_change_the_tree() {
        let bodies = random_bodies(7, 300);
        let tree = production_tree_of(&bodies, 1);
        for leaf in leaves(&tree) {
            for slot in leaf.acceleration.borrow_mut().iter_mut() {
                *slot = [1.0, 2.0, 3.0];
            }
        }
        let first = tree.accelerations(&bodies);
        let second = tree.accelerations(&bodies);
        assert_eq!(first, second);
    }

    // =================================================== inherit_potential

    use crate::fmm::complex::Complex;
    use crate::fmm::harmonic::to_index;
    use crate::fmm::test_util::distinct;

    fn field(p: usize, offset: f64) -> SolidHarmonic {
        let mut l = SolidHarmonic::regular(p);
        for (k, c) in l.coeffs_mut().iter_mut().enumerate() {
            *c = distinct(k, offset);
        }
        l
    }

    fn max_diff(a: &SolidHarmonic, b: &SolidHarmonic) -> f64 {
        a.coeffs()
            .iter()
            .zip(b.coeffs())
            .map(|(x, y)| (x.re - y.re).abs().max((x.im - y.im).abs()))
            .fold(0.0, f64::max)
    }

    #[test]
    fn every_node_ends_up_with_its_own_field_plus_its_parents_shifted_to_it() {
        let p = 3;
        let bodies = random_bodies(8, 700);
        let tree = production_tree_of(&bodies, p);
        let nodes = all_nodes(&tree);
        assert!(
            depth(&tree) >= 2,
            "need at least two levels to test inheriting twice"
        );

        // Give every node its own distinct field first, then inherit.
        let own: Vec<SolidHarmonic> = (0..nodes.len()).map(|k| field(p, 0.1 * k as f64)).collect();
        for (node, f) in nodes.iter().zip(&own) {
            *node.external_potential.borrow_mut() = f.clone();
        }
        tree.inherit_potential();

        let slot_of = |n: &FmmTree| nodes.iter().position(|m| std::ptr::eq(*m, n)).unwrap();
        for parent in &nodes {
            for child in &parent.children {
                let shift: Vec3 = from_fn(|i| child.system.center[i] - parent.system.center[i]);
                let mut want = own[slot_of(child)].clone();
                want += &local_to_local(&parent.external_potential.borrow(), shift);
                let got = child.external_potential.borrow();
                assert!(
                    max_diff(&got, &want) < 1e-9,
                    "child {} of parent {}",
                    slot_of(child),
                    slot_of(parent)
                );
            }
        }
        // the root has nobody above it, so it is unchanged
        assert!(max_diff(&tree.external_potential.borrow(), &own[0]) == 0.0);
    }

    #[test]
    fn a_children_inherit_from_the_already_updated_parent_not_the_stale_one() {
        // Only the root has a field. If a grandchild read its parent before the
        // parent had inherited, it would get nothing.
        let p = 2;
        let bodies = random_bodies(9, 700);
        let tree = production_tree_of(&bodies, p);
        assert!(depth(&tree) >= 2);
        *tree.external_potential.borrow_mut() = field(p, 0.4);
        tree.inherit_potential();

        for leaf in leaves(&tree) {
            let f = leaf.external_potential.borrow();
            assert!(
                f.coeffs().iter().any(|c| c.re != 0.0 || c.im != 0.0),
                "a leaf at depth below the root got nothing"
            );
        }
    }

    #[test]
    fn leaves_read_their_field_into_the_accelerations_they_already_hold() {
        // A field of c * z about the root has slope c along z everywhere, and
        // shifting its center cannot change that. So every body, in every leaf,
        // must gain exactly (0, 0, c), on top of whatever its slot already held.
        let p = 2;
        let c = 0.35;
        let bodies = random_bodies(10, 700);
        let tree = production_tree_of(&bodies, p);
        let mut linear = SolidHarmonic::regular(p);
        linear.coeffs_mut()[to_index(1, 0)] = Complex { re: c, im: 0.0 };
        *tree.external_potential.borrow_mut() = linear;
        for leaf in leaves(&tree) {
            for slot in leaf.acceleration.borrow_mut().iter_mut() {
                *slot = [1.0, 2.0, 3.0];
            }
        }

        tree.inherit_potential();

        let mut seen = 0;
        for leaf in leaves(&tree) {
            for slot in leaf.acceleration.borrow().iter() {
                assert!((slot[0] - 1.0).abs() < 1e-9, "x: {slot:?}");
                assert!((slot[1] - 2.0).abs() < 1e-9, "y: {slot:?}");
                assert!(
                    (slot[2] - (3.0 + crate::shared::G * c)).abs() < 1e-9,
                    "z: {slot:?}"
                );
                seen += 1;
            }
        }
        assert_eq!(seen, bodies.len());
    }

    #[test]
    fn inheriting_with_no_field_anywhere_leaves_the_accelerations_alone() {
        let bodies = random_bodies(11, 300);
        let tree = production_tree_of(&bodies, 2);
        tree.inherit_potential();
        assert_eq!(tree.accelerations(&bodies), vec![[0.0; 3]; bodies.len()]);
    }

    #[test]
    fn a_quadratic_field_at_the_root_pulls_every_body_in_proportion_to_its_offset_from_the_root() {
        // The field u (z^2 - (x^2 + y^2)/2) about the root has slope
        // u (-dx, -dy, 2 dz) at offset d from the root, however many levels of
        // inheriting lie between the root and the leaf. This follows the field
        // down the whole tree to a closed-form answer.
        let p = 3;
        let u = 0.8;
        let bodies = random_bodies(12, 700);
        let tree = production_tree_of(&bodies, p);
        assert!(depth(&tree) >= 2);
        let mut quadratic = SolidHarmonic::regular(p);
        quadratic.coeffs_mut()[to_index(2, 0)] = Complex { re: u, im: 0.0 };
        *tree.external_potential.borrow_mut() = quadratic;

        tree.inherit_potential();

        let got = tree.accelerations(&bodies);
        let mut biggest: f64 = 0.0;
        for (body, acc) in bodies.iter().zip(&got) {
            let d: Vec3 = from_fn(|i| body.pos[i] - tree.system.center[i]);
            let want = [-u * d[0], -u * d[1], 2.0 * u * d[2]];
            for i in 0..3 {
                assert!(
                    (acc[i] - crate::shared::G * want[i]).abs() < 1e-9,
                    "body at {:?}: got {acc:?}, want {want:?}",
                    body.pos
                );
            }
            biggest = biggest.max(want.iter().fold(0.0, |m: f64, c| m.max(c.abs())));
        }
        assert!(biggest > 0.1, "fixture too weak to mean anything");
    }
}
