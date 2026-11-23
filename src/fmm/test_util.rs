//! Helpers shared by the FMM test modules.

use crate::{
    fmm::complex::Complex,
    shared::{Body, RNG, Vec3},
};
use std::array::from_fn;

/// Asserts that both parts of `got` are within `tol` of `want`.
/// `what` names the case in the failure message.
pub fn assert_close(got: Complex, want: Complex, tol: f64, what: &str) {
    assert!(
        (got.re - want.re).abs() <= tol && (got.im - want.im).abs() <= tol,
        "{what}: got ({}, {}), want ({}, {})",
        got.re,
        got.im,
        want.re,
        want.im
    );
}

/// Euclidean distance between two points.
pub fn dist3(a: Vec3, b: Vec3) -> f64 {
    (0..3).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt()
}

/// `n` random bodies with positions in a cube of half-width `hw` around `center`.
pub fn bodies_around(center: Vec3, hw: f64, n: usize, seed: u64) -> Vec<Body> {
    let mut rng = RNG::new(seed);
    (0..n)
        .map(|_| {
            let b = Body::rand(&mut rng, hw);
            Body {
                pos: from_fn(|i| center[i] + b.pos[i]),
                ..b
            }
        })
        .collect()
}

/// The potential sum_j m_j / |x - x_j| (G left out, as in the expansions).
pub fn direct_potential(bodies: &[Body], x: Vec3) -> f64 {
    bodies.iter().map(|b| b.mass / b.distance(&x)).sum()
}
