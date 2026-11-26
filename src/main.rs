mod barnes_hut;
mod direct;
mod fmm;
mod shared;
use shared::Body;
use std::{array::from_fn, time::Instant};

use crate::shared::RNG;

fn main() {
    const TIME_STEP: f64 = 0.01;
    const TOTAL_TIME: f64 = 1.0f64;
    const NBODIES: usize = 1e4 as usize;
    const HALF_WIDTH: f64 = 1.0;
    let mut rng: RNG = RNG::new(0);
    let bodies: Vec<Body> =
        from_fn::<Body, NBODIES, _>(|_| Body::rand(&mut rng, HALF_WIDTH)).to_vec();
    println!("Simulating {NBODIES} bodies");

    // direct simulation
    let mut start: Instant = Instant::now();
    direct::simulate(&mut bodies.clone(), TIME_STEP, TOTAL_TIME);
    println!("direct time: {:?}", Instant::now() - start);

    // Barnes-Hut
    for resolution in [0.4_f64, 0.6, 0.8, 1.0] {
        println!("Simulating using Barnes-Hut with cell limit theta = {resolution}");
        start = Instant::now();
        barnes_hut::simulate(&mut bodies.clone(), TIME_STEP, TOTAL_TIME, resolution);
        println!("Barnes-Hut time: {:?}", Instant::now() - start)
    }

    // Fast-Multipole Method
    println!("Simulating using Fast Multipole Method");
    let order: usize = 4;
    let start: Instant = Instant::now();
    fmm::simulate(&mut bodies.clone(), TIME_STEP, TOTAL_TIME, order);
    println!("Barnes-Hut time: {:?}", Instant::now() - start)
}
