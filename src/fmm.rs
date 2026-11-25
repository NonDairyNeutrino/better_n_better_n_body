mod complex;
mod harmonic;
mod operators;
#[cfg(test)]
mod test_util;
mod tree;
use crate::shared::{Body, System, Vec3};
use tree::FmmTree;

fn simulate_step(bodies: &mut [Body], time_step: f64, p: usize) {
    // calculate accelerations
    let mut accelerations: Vec<Vec3> = Vec::with_capacity(bodies.len());
    let system: System = System::from_bodies(bodies);
    let mut tree: FmmTree = FmmTree::new(system, p);
    // other steps
    tree.accumulate_mass_distribution();
    tree.assemble_potential();
    tree.inherit_potential();

    for (body, acc) in bodies.iter_mut().zip(accelerations) {
        body.accelerate(&acc, time_step);
    }
}

pub fn simulate(system: &mut [Body], time_step: f64, total_time: f64, p: usize) {
    let nsteps: u64 = (total_time / time_step).floor() as u64;
    for _ in 0..nsteps {
        simulate_step(system, time_step, p);
    }
}
