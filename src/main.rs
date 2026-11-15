mod barnes_hut;
mod direct;
mod shared;
use std::time::Instant;

use shared::{Body, System};

fn main() {
    let dt: f64 = 0.01; // seconds
    let total_time: f64 = 1000.0f64;
    let mut bodies: Vec<Body> = vec![
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
    ];

    // direct simulation
    println!("Starting universe: {bodies:?}");
    let start: Instant = Instant::now();
    direct::simulate(&mut bodies, dt, total_time);
    let stop: Instant = Instant::now();
    println!("Final system: {bodies:?} in {:?}", stop - start);

    let universe: System = System::from_bodies(bodies.as_slice());
}
