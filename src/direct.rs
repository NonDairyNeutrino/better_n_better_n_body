use crate::shared::{Body, ORIGIN, Vec3, gravity};

fn calculate_acceleration_from_system(system: &[Body], pos: &Vec3) -> Vec3 {
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

pub fn simulate_step(system: &mut [Body], time_step: f64) {
    let mut accelerations: Vec<Vec3> = Vec::with_capacity(system.len());

    for this in system.iter() {
        let this_acc: Vec3 = calculate_acceleration_from_system(system, &this.pos);
        accelerations.push(this_acc);
    }

    for (this, acc) in system.iter_mut().zip(accelerations) {
        this.accelerate(&acc, time_step);
    }
}

pub fn simulate(system: &mut [Body], time_step: f64, total_time: f64) {
    let nsteps: u64 = (total_time / time_step).floor() as u64;
    for _ in 0..nsteps {
        simulate_step(system, time_step);
    }
}
