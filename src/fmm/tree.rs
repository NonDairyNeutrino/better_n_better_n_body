use crate::shared::System;

#[derive(Debug, Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
}
impl Complex {
    fn zero() -> Complex {
        Complex {
            re: 0.0_f64,
            im: 0.0_f64,
        }
    }
}

#[derive(Debug, Clone)]
enum SolidHarmonic {
    Regular(Vec<Complex>),
    Irregular(Vec<Complex>),
}
impl SolidHarmonic {
    fn regular(p: usize) -> SolidHarmonic {
        let coeffs: Vec<Complex> = vec![Complex::zero(); (p + 1).pow(2)];
        return SolidHarmonic::Regular(coeffs);
    }
    fn irregular(p: usize) -> SolidHarmonic {
        let coeffs: Vec<Complex> = vec![Complex::zero(); (p + 1).pow(2)];
        return SolidHarmonic::Irregular(coeffs);
    }
    fn coeffs(&self) -> &[Complex] {
        match self {
            SolidHarmonic::Regular(v) | SolidHarmonic::Irregular(v) => v,
        }
    }
    fn coeffs_mut(&mut self) -> &mut [Complex] {
        match self {
            SolidHarmonic::Regular(v) | SolidHarmonic::Irregular(v) => v,
        }
    }
    fn get(&self, n: i64, m: i64) -> Complex {
        assert!(0 <= n, "n must be non-negative, got n={n}");
        assert!(-n <= m && m <= n, "need |m| <= n, got n={n}, m={m}");
        let idx: usize = (n.pow(2) + n + m) as usize;
        assert!(
            idx < self.coeffs().len(),
            "(n={n}, m={m}) is beyond the order of this solid harmonic"
        );
        let coeff: Complex = self.coeffs()[idx];
        return coeff;
    }
}

const MAX_LEAF_BODIES: usize = 2;
pub struct FmmTree<'a> {
    system: System<'a>,
    far_field: SolidHarmonic,   // the field far from the observer
    local_field: SolidHarmonic, // the field near the observer
    children: Vec<FmmTree<'a>>,
}
impl<'a> FmmTree<'a> {
    pub fn new(system: System<'a>, p: usize) -> FmmTree<'a> {
        let far_field: SolidHarmonic = SolidHarmonic::irregular(p);
        let local_field: SolidHarmonic = SolidHarmonic::regular(p);

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
                far_field,
                local_field,
            };
        } else {
            // don't split
            let children: Vec<FmmTree<'a>> = Vec::new();

            return FmmTree {
                system,
                children,
                far_field,
                local_field,
            };
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

    /// The number of (n, m) pairs with 0 <= n <= p and |m| <= n, counted directly
    /// rather than with the closed formula the code uses.
    fn pair_count(p: usize) -> usize {
        (0..=p).map(|n| 2 * n + 1).sum()
    }

    // ------------------------------------------------------- SolidHarmonic

    #[test]
    fn a_solid_harmonic_has_one_coefficient_per_pair_n_m() {
        for p in 0..=8 {
            assert_eq!(SolidHarmonic::regular(p).coeffs().len(), pair_count(p), "regular p={p}");
            assert_eq!(SolidHarmonic::irregular(p).coeffs().len(), pair_count(p), "irregular p={p}");
        }
    }

    #[test]
    fn the_constructors_make_the_variant_they_are_named_for() {
        assert!(matches!(SolidHarmonic::regular(2), SolidHarmonic::Regular(_)));
        assert!(matches!(SolidHarmonic::irregular(2), SolidHarmonic::Irregular(_)));
    }

    #[test]
    fn a_new_solid_harmonic_is_all_zeros() {
        for sh in [SolidHarmonic::regular(4), SolidHarmonic::irregular(4)] {
            assert!(!sh.coeffs().is_empty());
            assert!(sh.coeffs().iter().all(|c| c.re == 0.0 && c.im == 0.0));
        }
    }

    /// Fill slot k with the number k, so `get(n, m)` reports which slot it read.
    fn numbered(p: usize) -> SolidHarmonic {
        let mut sh = SolidHarmonic::regular(p);
        for (k, c) in sh.coeffs_mut().iter_mut().enumerate() {
            *c = Complex { re: k as f64, im: -(k as f64) };
        }
        sh
    }

    #[test]
    fn get_reads_the_slots_in_the_documented_order() {
        // (0,0) | (1,-1) (1,0) (1,1) | (2,-2) (2,-1) (2,0) (2,1) (2,2) | ...
        let sh = numbered(3);
        let expected = [
            (0, 0, 0.0),
            (1, -1, 1.0),
            (1, 0, 2.0),
            (1, 1, 3.0),
            (2, -2, 4.0),
            (2, 0, 6.0),
            (2, 2, 8.0),
            (3, -3, 9.0),
            (3, 3, 15.0),
        ];
        for (n, m, slot) in expected {
            assert_eq!(sh.get(n, m).re, slot, "(n={n}, m={m})");
            assert_eq!(sh.get(n, m).im, -slot, "(n={n}, m={m}) imaginary part");
        }
    }

    #[test]
    fn every_pair_n_m_has_its_own_slot_and_every_slot_is_used() {
        let p: i64 = 6;
        let sh = numbered(p as usize);
        let mut seen: Vec<usize> = Vec::new();
        for n in 0..=p {
            for m in -n..=n {
                seen.push(sh.get(n, m).re as usize);
            }
        }
        seen.sort();
        let want: Vec<usize> = (0..pair_count(p as usize)).collect();
        assert_eq!(seen, want);
    }

    /// The formula keeps working far past the orders anyone will use.
    #[test]
    fn the_first_and_last_pairs_of_a_high_order_are_the_first_and_last_slots() {
        let p: i64 = 60;
        let sh = numbered(p as usize);
        assert_eq!(sh.get(0, 0).re as usize, 0);
        assert_eq!(sh.get(p, -p).re as usize, (p * p) as usize);
        assert_eq!(sh.get(p, p).re as usize, pair_count(p as usize) - 1);
    }

    #[test]
    #[should_panic(expected = "|m| <= n")]
    fn get_rejects_m_larger_than_n() {
        numbered(3).get(1, 2);
    }

    #[test]
    #[should_panic(expected = "|m| <= n")]
    fn get_rejects_m_smaller_than_minus_n() {
        numbered(3).get(1, -2);
    }

    #[test]
    #[should_panic(expected = "n must be non-negative")]
    fn get_rejects_a_negative_n() {
        numbered(3).get(-1, 0);
    }

    #[test]
    #[should_panic(expected = "beyond the order")]
    fn get_rejects_an_n_beyond_the_order() {
        numbered(3).get(4, 0);
    }

    #[test]
    fn writes_through_coeffs_mut_are_seen_by_coeffs() {
        let mut sh = SolidHarmonic::irregular(2);
        sh.coeffs_mut()[5] = Complex { re: 1.5, im: -2.5 };
        assert_eq!(sh.coeffs()[5].re, 1.5);
        assert_eq!(sh.coeffs()[5].im, -2.5);
        assert!(matches!(sh, SolidHarmonic::Irregular(_)), "kind must not change");
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
            for node in all_nodes(&tree).into_iter().filter(|n| !n.children.is_empty()) {
                internal_seen += 1;
                assert_eq!(node.children.len(), 8, "{name}");
                assert!(node.system.bodies.is_empty(), "{name}: internal nodes are bodyless");
            }
        }
        assert!(internal_seen > 100, "vacuous: only {internal_seen} internal nodes");
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
        assert_eq!(biggest, MAX_LEAF_BODIES, "the limit is never reached: a vacuous check");
    }

    #[test]
    fn the_leaves_hold_every_body_exactly_once() {
        for (name, bodies) in fixtures() {
            let tree = tree_of(&bodies, 2);
            let got = addresses(leaves(&tree).iter().flat_map(|l| l.system.bodies.iter().copied()));
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
                    assert!(strictly_inside(leaf.system.center, leaf.system.hw, b.pos), "{name}");
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
        for node in all_nodes(&tree).into_iter().filter(|n| !n.children.is_empty()) {
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

    #[test]
    fn every_node_starts_with_zeroed_coefficients_of_the_right_length() {
        for p in [0, 1, 3] {
            let bodies = random_bodies(9, 60);
            let tree = tree_of(&bodies, p);
            for node in all_nodes(&tree) {
                for sh in [&node.far_field, &node.local_field] {
                    assert_eq!(sh.coeffs().len(), pair_count(p), "p={p}");
                    assert!(sh.coeffs().iter().all(|c| c.re == 0.0 && c.im == 0.0));
                }
            }
        }
    }

    /// far_field (valid far from the cell, made by its own mass) is irregular;
    /// local_field (valid near the cell, made by far mass) is regular. This checks
    /// the code against that convention, not the convention against the physics.
    #[test]
    fn every_node_has_the_right_kind_in_each_field() {
        let bodies = random_bodies(9, 60);
        let tree = tree_of(&bodies, 2);
        for node in all_nodes(&tree) {
            assert!(matches!(node.far_field, SolidHarmonic::Irregular(_)));
            assert!(matches!(node.local_field, SolidHarmonic::Regular(_)));
        }
    }
}
