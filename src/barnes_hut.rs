mod tree;
use crate::shared::{Body, System, Vec3};
use tree::Tree;

fn simulate_step(bodies: &mut [Body], dt: f64, resolution: f64) {
    // build the root
    let universe: System = System::from_bodies(bodies);
    let tree: Tree = Tree::new(universe);

    // calculate accelerations
    let mut accelerations: Vec<Vec3> = Vec::with_capacity(bodies.len());
    for this in bodies.iter() {
        let acc: Vec3 = tree.accel_at(&this.pos, resolution);
        accelerations.push(acc);
    }

    // update the bodies
    for (this, acc) in bodies.iter_mut().zip(accelerations) {
        this.accelerate(&acc, dt);
    }
}

/// Advances `system` by `total_time` in steps of `time_step` using Barnes-Hut.
///
/// `resolution` is the opening angle θ: a cell is treated as one point mass when
/// `width <= resolution * distance`. θ = 0 reproduces direct summation; larger
/// values are faster and less accurate.
///
/// # Input contract
/// These are assumed, not checked:
/// - all masses are positive;
/// - no two bodies share a position (self-skip compares positions);
/// - bodies are at least 1e-10 apart (closer ones merge into one point mass);
/// - no body lies exactly on a cell boundary;
/// - units are nondimensionalized (G = 1).
pub fn simulate(system: &mut [Body], time_step: f64, total_time: f64, resolution: f64) {
    let nsteps: u64 = (total_time / time_step).floor() as u64;
    for _ in 0..nsteps {
        simulate_step(system, time_step, resolution);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::direct;
    use crate::shared::test_util::*;

    const DT: f64 = 1.0;

    /// `scattered(n)` with masses 1..=5, so mass weighting matters. Note that
    /// the input contract holds: positive masses, no two bodies at one position,
    /// and `scattered` keeps bodies off the cell faces.
    fn universe(n: usize) -> Vec<Body> {
        scattered(n)
            .iter()
            .enumerate()
            .map(|(i, b)| body_with(1.0 + (i % 5) as f64, b.pos[0], b.pos[1], b.pos[2]))
            .collect()
    }

    /// Size of the velocity difference between two evolved universes, relative to
    /// the size of the reference velocities, over ALL bodies at once. (A per-body
    /// relative error is useless here: a body near the middle of the cloud has an
    /// almost-zero net force, so any error looks enormous next to it.)
    ///
    /// Velocities carry the force signal. G is ~1e-11, so the positions barely
    /// move and would match for any force law.
    fn velocity_error(approx: &[Body], reference: &[Body]) -> f64 {
        let diff: f64 = approx
            .iter()
            .zip(reference)
            .map(|(a, r)| (0..3).map(|k| (a.vel[k] - r.vel[k]).powi(2)).sum::<f64>())
            .sum();
        let size: f64 = reference
            .iter()
            .map(|r| (0..3).map(|k| r.vel[k].powi(2)).sum::<f64>())
            .sum();
        assert!(size > 0.0, "reference velocities are all zero: vacuous");
        (diff / size).sqrt()
    }

    fn one_step_error(n: usize, theta: f64) -> f64 {
        let mut reference = universe(n);
        direct::simulate_step(&mut reference, DT);
        let mut approx = universe(n);
        simulate_step(&mut approx, DT, theta);
        velocity_error(&approx, &reference)
    }

    // -------------------------------------------------------------- theta = 0

    #[test]
    fn theta_zero_matches_direct() {
        // the tree is never trusted, so every body sums every other body exactly
        for n in [2, 8, 20, 100, 300] {
            let err = one_step_error(n, 0.0);
            assert!(err < 1e-12, "n={n}: error {err:e}");
        }
    }

    #[test]
    fn theta_zero_matches_direct_over_several_steps() {
        // also checks that the tree is rebuilt from the moved bodies every step
        let mut reference = universe(50);
        let mut approx = universe(50);
        direct::simulate(&mut reference, DT, 5.0);
        simulate(&mut approx, DT, 5.0, 0.0);
        let err = velocity_error(&approx, &reference);
        assert!(err < 1e-10, "error {err:e}");
    }

    // --------------------------------------------------------------- theta > 0

    /// Not vacuous: a bound on the error means nothing if θ > 0 is secretly
    /// the same computation as θ = 0.
    #[test]
    fn a_large_theta_really_approximates() {
        let err = one_step_error(300, 1.0);
        assert!(
            err > 1e-3,
            "θ=1 error was only {err:e}; nothing is being approximated"
        );
    }

    #[test]
    fn a_small_theta_stays_close_to_direct() {
        // measured about 0.2% for this fixture; allow 1%
        let err = one_step_error(300, 0.3);
        assert!(err < 0.01, "error {err:e}");
    }

    #[test]
    fn error_grows_as_theta_grows() {
        let errs: Vec<f64> = [0.3, 0.5, 0.7, 1.0]
            .iter()
            .map(|&t| one_step_error(300, t))
            .collect();
        assert!(errs.windows(2).all(|w| w[0] < w[1]), "errors were {errs:?}");
    }

    // ---------------------------------------------------------------- hygiene

    #[test]
    fn no_velocity_is_nan_or_infinite_for_any_theta() {
        // a body that feels itself divides by zero and poisons everything
        for theta in [0.0, 0.5, 1.0, 5.0] {
            let mut bodies = universe(100);
            simulate_step(&mut bodies, DT, theta);
            for (i, b) in bodies.iter().enumerate() {
                assert!(b.vel.iter().all(|v| v.is_finite()), "θ={theta}, body {i}");
            }
        }
    }
}
