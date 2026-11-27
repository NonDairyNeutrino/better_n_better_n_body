//! Direct summation: every body feels every other body.
//!
//! This is the simplest method and the reference the faster ones are checked
//! against. It is exact up to the time integrator, and costs O(N²) per step.

use crate::shared::{Body, ORIGIN, Vec3, gravity};

/// The gravitational acceleration at `pos` due to every body in `system`.
///
/// A body whose position equals `pos` is skipped. That is how a body avoids
/// feeling itself, so it also means two bodies at exactly the same place ignore
/// each other.
pub fn calculate_acceleration_from_system<'b>(
    system: impl IntoIterator<Item = &'b Body>,
    pos: &Vec3,
) -> Vec3 {
    let mut this_acc: Vec3 = ORIGIN;

    for that in system {
        if that.pos != *pos {
            let that_acc_at_this: Vec3 = gravity(&that, &pos);

            for dim in 0..this_acc.len() {
                this_acc[dim] += that_acc_at_this[dim];
            }
        }
    }
    return this_acc;
}

/// Advances every body by one time step `time_step`.
///
/// All accelerations are computed from the positions at the start of the step,
/// and only then are the bodies moved, so the answer does not depend on the
/// order of the bodies.
pub fn simulate_step(system: &mut [Body], time_step: f64) {
    let mut accelerations: Vec<Vec3> = Vec::with_capacity(system.len());

    for this in system.iter() {
        let this_acc: Vec3 = calculate_acceleration_from_system(system.iter(), &this.pos);
        accelerations.push(this_acc);
    }

    for (this, acc) in system.iter_mut().zip(accelerations) {
        this.accelerate(&acc, time_step);
    }
}

/// Advances `system` by `total_time` in steps of `time_step` using direct summation.
///
/// The number of steps is `total_time / time_step`, rounded down. Progress is
/// printed to standard error.
pub fn simulate(system: &mut [Body], time_step: f64, total_time: f64) {
    let nsteps: u64 = (total_time / time_step).floor() as u64;
    for step in 0..nsteps {
        eprint!("simulating step {step}/{nsteps}\r");
        simulate_step(system, time_step);
    }
}
