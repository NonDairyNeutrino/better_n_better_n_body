#![doc = include_str!("../README.md")]
//!
#![doc = include_str!("../PHYSICS.md")]
//!
//! # Reading the code
//!
//! The modules, simplest first: [`direct`], [`barnes_hut`], [`fmm`]. What all
//! three share (the [`Body`] type, the [`System`](shared::System) container and
//! the force law [`gravity`](shared::gravity)) is in [`shared`].
//!
//! # Building these docs
//!
//! This is a binary crate and its modules are private, so plain `cargo doc`
//! renders almost nothing. Use
//!
//! ```sh
//! cargo doc --document-private-items --open
//! ```

mod barnes_hut;
mod direct;
mod fmm;
mod shared;
use shared::Body;
use std::time::Instant;

use crate::shared::RNG;

fn main() {
    const NBODIES: usize = 1e5 as usize;
    const TIME_STEP: f64 = 0.1;
    const TOTAL_TIME: f64 = 1.0f64;
    const HALF_WIDTH: f64 = 1.0;
    let mut rng: RNG = RNG::new(0);
    let bodies: Vec<Body> = (0..NBODIES)
        .map(|_| Body::rand(&mut rng, HALF_WIDTH))
        .collect();
    println!("Simulating {NBODIES} bodies");

    // direct simulation
    println!("Simulating using direct method");
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
    let order: usize = 2;
    println!("Simulating using Fast Multipole Method with order {order}");
    start = Instant::now();
    fmm::simulate(&mut bodies.clone(), TIME_STEP, TOTAL_TIME, order);
    println!("Fast-Multipole time: {:?}", Instant::now() - start)
}
