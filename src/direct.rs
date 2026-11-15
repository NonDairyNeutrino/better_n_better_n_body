use crate::shared::{Body, gravity};

fn calculate_acceleration_from_system(system: &[Body], pos: &[f64; 3]) -> [f64; 3] {
    let mut this_acc: [f64; 3] = [0.0_f64, 0.0_f64, 0.0_f64];

    for that in system {
        if that.pos != *pos {
            let that_acc_at_this: [f64; 3] = gravity(&that, &pos);

            for dim in 0..this_acc.len() {
                this_acc[dim] += that_acc_at_this[dim];
            }
        }
    }
    return this_acc;
}

pub fn simulate_step(system: &mut [Body], time_step: f64) {
    // need old_system or else later objects will use new positions
    let mut old_system: Vec<Body> = Vec::new();
    old_system.resize(system.len(), Body::new());
    old_system.clone_from_slice(system);

    for this in system {
        let this_acc: [f64; 3] =
            calculate_acceleration_from_system(old_system.as_slice(), &this.pos);

        this.accelerate(&this_acc, time_step);
    }
}

pub fn simulate(system: &mut [Body], time_step: f64, total_time: f64) {
    let nsteps: u64 = (total_time / time_step).floor() as u64;
    println!("Stepping {nsteps} times");
    for _ in 0..nsteps {
        simulate_step(system, time_step);
    }
}
