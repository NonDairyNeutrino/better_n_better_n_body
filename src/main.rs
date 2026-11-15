mod barnes_hut;
mod direct;
mod shared;
use std::time::Instant;

use shared::{Body, System};

fn main() {
    let dt: f64 = 0.01; // seconds
    let total_time: f64 = 1000.0f64;
    let mut system: System = System::from_bodies(
        vec![
            Body {
                mass: 1.0,
                pos: [1.0_f64, 0.0_f64, 0.0_f64],
                vel: [0.0_f64, 0.0_f64, 0.0_f64],
            },
            Body {
                mass: 1.0,
                pos: [0.0_f64, 1.0_f64, 0.0_f64],
                vel: [0.0_f64, 0.0_f64, 0.0_f64],
            },
            Body {
                mass: 1.0,
                pos: [0.0_f64, 0.0_f64, 1.0_f64],
                vel: [0.0_f64, 0.0_f64, 0.0_f64],
            },
        ]
        .as_slice(),
    );
    println!("Starting system: {system:?}");

    let start: Instant = Instant::now();
    direct::simulate(system.bodies.as_mut_slice(), dt, total_time);
    let stop: Instant = Instant::now();
    println!("Final system: {system:?} in {:?}", stop - start);

    let bounding_box: ([f64; 3], f64) = shared::find_bounding_box(system.bodies.as_slice());
    println!("The bounding box for the system is {bounding_box:?}")
}
