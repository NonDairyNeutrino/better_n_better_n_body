mod barnes_hut;
mod direct;
mod shared;
use shared::{Body, System};
use std::time::Instant;

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
    let start_direct: Instant = Instant::now();
    direct::simulate(&mut bodies, dt, total_time);
    println!(
        "Final system: {bodies:?} in {:?}",
        Instant::now() - start_direct
    );

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
    let resolution: f64 = 0.1_f64;
    println!("Simulating with fresh system using Barnes-Hut");
    let start_bh: Instant = Instant::now();
    barnes_hut::simulate(&mut bodies, dt, total_time, resolution);
    println!(
        "Final Barnes-Hut: {:?}\nTime taken: {:?}",
        bodies,
        Instant::now() - start_bh
    )
}
