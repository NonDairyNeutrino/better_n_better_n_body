mod barnes_hut;
mod direct;
mod shared;
use shared::Body;
use std::{array::from_fn, time::Instant};

use crate::shared::RNG;

fn main() {
    const TIME_STEP: f64 = 0.01;
    const TOTAL_TIME: f64 = 1.0f64;
    const NBODIES: usize = 10_000;
    const HALF_WIDTH: f64 = 1.0;
    let mut rng: RNG = RNG::new(0);
    let mut bodies: Vec<Body> =
        from_fn::<Body, NBODIES, _>(|_| Body::rand(&mut rng, HALF_WIDTH)).to_vec();

    // direct simulation
    let start_direct: Instant = Instant::now();
    direct::simulate(&mut bodies, TIME_STEP, TOTAL_TIME);
    println!("direct time: {:?}", Instant::now() - start_direct);

    // Barnes-Hut
    let mut rng: RNG = RNG::new(0);
    let mut bodies: Vec<Body> =
        from_fn::<Body, NBODIES, _>(|_| Body::rand(&mut rng, HALF_WIDTH)).to_vec();
    let resolution: f64 = 0.5_f64;
    println!("Simulating with fresh system using Barnes-Hut");
    let start_bh: Instant = Instant::now();
    barnes_hut::simulate(&mut bodies, TIME_STEP, TOTAL_TIME, resolution);
    println!("Barnes-Hut time: {:?}", Instant::now() - start_bh)
}
