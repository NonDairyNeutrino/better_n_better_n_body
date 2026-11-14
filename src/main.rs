mod barnes_hut;
mod direct;
mod shared;
use std::time::Instant;

use shared::Object;

fn main() {
    let dt: f64 = 0.01; // seconds
    let mut system: [Object; 3] = [
        Object {
            mass: 1.0,
            pos: [1.0_f64, 0.0_f64, 0.0_f64],
            vel: [0.0_f64, 0.0_f64, 0.0_f64],
        },
        Object {
            mass: 1.0,
            pos: [0.0_f64, 1.0_f64, 0.0_f64],
            vel: [0.0_f64, 0.0_f64, 0.0_f64],
        },
        Object {
            mass: 1.0,
            pos: [0.0_f64, 0.0_f64, 1.0_f64],
            vel: [0.0_f64, 0.0_f64, 0.0_f64],
        },
    ];
    println!("Starting system: {system:?}");
    let total_time: f64 = 1000.0f64;
    let start: Instant = Instant::now();
    direct::simulate(&mut system, dt, total_time);
    let stop: Instant = Instant::now();
    println!("Final system: {system:?} in {:?}", stop - start);

    let bounding_box: ([f64; 3], f64) = barnes_hut::find_bounding_box(system.as_slice());
    println!("The bounding box for the system is {bounding_box:?}")
}
