mod direct;
mod shared;
use shared::Object;

use crate::direct::simulate;

fn main() {
    let ob: Object = Object::new();
    let dt: f64 = 0.01; // seconds
    let mut system: [Object; 2] = [
        ob,
        Object {
            mass: 1.0,
            pos: [1.0_f64, 0.0_f64, 0.0_f64],
            vel: [0.0_f64, 0.0_f64, 0.0_f64],
        },
    ];
    let total_time: f64 = 1000.0f64;
    simulate(&mut system, dt, total_time);
    println!("{system:?}")
}
