mod direct;
mod shared;
use shared::{Object, gravity};

use crate::direct::simulate_step;

fn main() {
    let mut ob: Object = Object::new();
    println!("{ob:?}");
    let g: [f64; 3] = gravity(&ob, &[1f64, 0f64, 0f64]);
    println!("gravity 1m away from 1kg is {g:?}");
    let dt: f64 = 0.01; // seconds
    ob.accelerate(&g, dt);
    println!("{ob:?}");
    let mut system: [Object; 2] = [
        ob,
        Object {
            mass: 1.0,
            pos: [1.0_f64, 0.0_f64, 0.0_f64],
            vel: [0.0_f64, 0.0_f64, 0.0_f64],
        },
    ];
    simulate_step(system.as_mut_slice(), dt);
    println!("{system:?}")
}
