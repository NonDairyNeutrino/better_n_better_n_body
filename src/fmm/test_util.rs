//! Helpers shared by the FMM test modules.

use crate::{fmm::complex::Complex, shared::Vec3};

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
