//! The Fast Multipole Method (FMM).
//!
//! Barnes-Hut still makes every body consult the far cells one by one. The FMM
//! does that work once per *cell*: each cell builds the smooth field that all the
//! distant masses create over it, and every body in the cell reads its
//! acceleration from that field. Bodies in neighbouring cells are summed
//! directly. The work per body stops growing with N, so a step costs O(N).
//!
//! # One time step
//!
//! The steps below are [`FmmTree`] methods, called in this order
//! by `simulate_step`. The standard FMM name of each operator is in brackets;
//! see [`operators`] for what each does.
//!
//! 1. **Build** the octree. A cell is split while it holds more than a fixed
//!    number of bodies.
//! 2. **Up** (`accumulate_mass_distribution`): summarize the mass in every cell,
//!    as seen from far away. Leaves summarize their bodies
//!    ([P2M](operators::particle_to_multipole)); parents combine their children's
//!    summaries ([M2M](operators::multipole_to_multipole)).
//! 3. **Across** (`assemble_potential`): one traversal over pairs of cells. Cells
//!    that are far enough apart convert the source's summary into a field over
//!    the target ([M2L](operators::multipole_to_local)). Near leaves sum their
//!    bodies directly ([P2P](operators::particle_to_particle)).
//! 4. **Down** (`inherit_potential`): each cell passes its field to its
//!    children ([L2L](operators::local_to_local)), and each leaf reads an
//!    acceleration for every body from its field
//!    ([L2P](operators::local_to_particle)).
//!
//! # Accuracy
//!
//! The knob is the expansion order `p`. Every cell stores (p + 1)² coefficients
//! per expansion; a higher order is slower and more accurate. See [`harmonic`]
//! for what the coefficients are.
//!
//! For the idea in everyday terms see `src/fmm/README.md`, and for the physics
//! in words, `src/fmm/PHYSICS.md`.

mod complex;
mod harmonic;
mod operators;
#[cfg(test)]
mod test_util;
mod tree;
use crate::shared::{Body, System, Vec3};
use tree::FmmTree;

/// One FMM time step: build the tree, run the up / across / down passes, then
/// move the bodies with the accelerations that come out.
fn simulate_step(bodies: &mut [Body], time_step: f64, p: usize) {
    // calculate accelerations
    let system: System = System::from_bodies(bodies);
    let mut tree: FmmTree = FmmTree::new(system, p);
    // other steps
    tree.accumulate_mass_distribution();
    tree.assemble_potential();
    tree.inherit_potential();
    let accelerations: Vec<Vec3> = tree.accelerations(bodies);

    for (body, acc) in bodies.iter_mut().zip(accelerations) {
        body.accelerate(&acc, time_step);
    }
}

/// Advances `system` by `total_time` in steps of `time_step` using the Fast Multipole Method.
///
/// `p` is the expansion order. Higher is more accurate and slower. At order 0
/// each far field is a constant, which has no gradient, so only the near-field
/// (direct) pulls remain; use order 1 or more for a meaningful simulation.
///
/// The same input assumptions as [`barnes_hut::simulate`](crate::barnes_hut::simulate)
/// apply: positive masses, and no two bodies at the same position.
pub fn simulate(system: &mut [Body], time_step: f64, total_time: f64, p: usize) {
    let nsteps: u64 = (total_time / time_step).floor() as u64;
    for step in 0..nsteps {
        eprint!("simulating step {step}/{nsteps}\r");
        simulate_step(system, time_step, p);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::direct;
    use crate::fmm::test_util::{bodies_around, pull_from};
    use crate::shared::ORIGIN;

    // --------------------------------------------------------------- fixtures

    /// `n` bodies at rest in a cube of half-width 1, so that after one step of
    /// length `dt` each velocity is exactly `acceleration * dt`.
    fn at_rest(n: usize, seed: u64) -> Vec<Body> {
        bodies_around(ORIGIN, 1.0, n, seed)
            .into_iter()
            .map(|b| Body { vel: [0.0; 3], ..b })
            .collect()
    }

    /// The acceleration of every body, read back from one FMM step.
    fn fmm_accelerations(bodies: &[Body], p: usize) -> Vec<Vec3> {
        const DT: f64 = 1.0;
        let mut stepped: Vec<Body> = bodies.to_vec();
        simulate_step(&mut stepped, DT, p);
        stepped.iter().map(|b| b.vel.map(|v| v / DT)).collect()
    }

    /// Newton's law written out, summed over every other body. Independent of
    /// `direct` and of `shared::gravity`.
    fn newton(bodies: &[Body]) -> Vec<Vec3> {
        bodies
            .iter()
            .enumerate()
            .map(|(i, target)| {
                let mut acc = [0.0; 3];
                for (j, source) in bodies.iter().enumerate() {
                    if i != j {
                        let a = pull_from(target.pos, source);
                        for d in 0..3 {
                            acc[d] += a[d];
                        }
                    }
                }
                acc
            })
            .collect()
    }

    /// Root-mean-square error of `got` relative to the root-mean-square size of `want`.
    fn relative_error(got: &[Vec3], want: &[Vec3]) -> f64 {
        let mut err = 0.0;
        let mut size = 0.0;
        for (g, w) in got.iter().zip(want) {
            for d in 0..3 {
                err += (g[d] - w[d]).powi(2);
                size += w[d].powi(2);
            }
        }
        (err / size).sqrt()
    }

    // ------------------------------------------------------------------ exact

    // At or below the leaf limit the tree is a single leaf, so every pair is
    // handled by direct summation and there is nothing left to approximate.

    #[test]
    fn a_system_that_fits_in_one_leaf_matches_newton_exactly() {
        let bodies = at_rest(30, 11);
        let got = fmm_accelerations(&bodies, 3);
        let want = newton(&bodies);
        assert!(relative_error(&got, &want) < 1e-12);
    }

    #[test]
    fn simulate_agrees_with_direct_when_the_tree_is_one_leaf() {
        // 0.35 / 0.1 is three whole steps; the same bodies stepped the same way
        // by `direct` must end up in the same place.
        let start = at_rest(20, 12);
        let mut by_fmm = start.clone();
        let mut by_direct = start.clone();
        simulate(&mut by_fmm, 0.1, 0.35, 2);
        direct::simulate(&mut by_direct, 0.1, 0.35);
        for (a, b) in by_fmm.iter().zip(&by_direct) {
            for d in 0..3 {
                assert!(
                    (a.pos[d] - b.pos[d]).abs() < 1e-9,
                    "{:?} vs {:?}",
                    a.pos,
                    b.pos
                );
                assert!(
                    (a.vel[d] - b.vel[d]).abs() < 1e-9,
                    "{:?} vs {:?}",
                    a.vel,
                    b.vel
                );
            }
        }
    }

    // ------------------------------------------------------------- approximate

    // The bounds below are about twice the worst error seen over several seeds
    // (1.4e-2, 2.8e-3, 5.7e-4, 1.1e-4 for p = 1..4), so they catch a real loss of
    // accuracy without failing on a different random fixture. The reference is
    // Newton's law summed over every pair.

    #[test]
    fn the_error_against_newton_stays_within_a_bound_at_each_order() {
        let bodies = at_rest(1500, 21);
        let want = newton(&bodies);
        for (p, bound) in [(1, 3e-2), (2, 6e-3), (3, 1.2e-3), (4, 2.5e-4)] {
            let err = relative_error(&fmm_accelerations(&bodies, p), &want);
            assert!(err < bound, "p = {p}: error {err:.3e} exceeds {bound:.1e}");
        }
    }

    #[test]
    fn raising_the_order_shrinks_the_error() {
        let bodies = at_rest(1500, 21);
        let want = newton(&bodies);
        let errors: Vec<f64> = (1..=4)
            .map(|p| relative_error(&fmm_accelerations(&bodies, p), &want))
            .collect();
        for pair in errors.windows(2) {
            assert!(
                pair[1] < pair[0] / 2.0,
                "error did not at least halve: {errors:?}"
            );
        }
    }

    #[test]
    fn the_far_field_really_is_approximated_at_the_lowest_order() {
        // Guard against the bounds above passing only because the tree never
        // made a far-field conversion: at p = 1 the answer must visibly differ
        // from the exact one.
        let bodies = at_rest(1500, 21);
        let want = newton(&bodies);
        let err = relative_error(&fmm_accelerations(&bodies, 1), &want);
        assert!(err > 1e-3, "p = 1 was suspiciously exact: {err:.3e}");
    }

    // ------------------------------------------------------------------ order

    #[test]
    fn each_body_gets_its_own_acceleration_whatever_order_the_input_is_in() {
        let bodies = at_rest(400, 31);
        let forward = fmm_accelerations(&bodies, 2);

        let mut reversed_bodies = bodies.clone();
        reversed_bodies.reverse();
        let mut backward = fmm_accelerations(&reversed_bodies, 2);
        backward.reverse();

        // Same bodies, same box, same tree: only the order of the sums differs.
        assert!(relative_error(&backward, &forward) < 1e-9);
    }

    #[test]
    fn a_step_changes_every_velocity_and_every_position() {
        // Nothing should be skipped: a zero acceleration or an unmoved body
        // means a slot was never filled.
        let mut bodies = at_rest(400, 32);
        let before = bodies.clone();
        simulate_step(&mut bodies, 0.01, 2);
        for (a, b) in bodies.iter().zip(&before) {
            assert!(
                a.vel.iter().any(|v| *v != 0.0),
                "velocity unchanged: {:?}",
                a.vel
            );
            assert!(a.pos != b.pos, "position unchanged: {:?}", a.pos);
        }
    }
}
