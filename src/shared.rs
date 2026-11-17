//! This is a module documentation
mod body;
mod system;
pub use body::Body;
pub use system::System;

pub const G: f64 = 1.0_f64; // non-dimensionalized, otherwise 6.67 * 1e-11 in SI
pub type Vec3 = [f64; 3];
pub const ORIGIN: Vec3 = [0.0_f64, 0.0_f64, 0.0_f64];

// -------------------------------------------------------------------------------- gravity
pub fn gravity(object: &Body, pos: &Vec3) -> Vec3 {
    let strength: f64 = G * object.mass / object.distance(pos).powi(3); // 3 is on purpose
    let mut grav: Vec3 = [0.0_f64, 0.0_f64, 0.0_f64];

    for dim in 0..grav.len() {
        grav[dim] = strength * (object.pos[dim] - pos[dim]);
    }

    return grav;
}

// -------------------------------------------------------------------------------------RNG
// standard LCG constants
const A: u64 = 6364136223846793005;
const B: u64 = 1442695040888963407;
// C = 2^-53
const C: f64 = 1.1102230246251565e-16;

pub struct RNG {
    state: u64,
}
impl RNG {
    pub fn new(state: u64) -> RNG {
        return RNG { state };
    }
    // Linear-Congruentional Generator
    pub fn rand(&mut self) -> f64 {
        self.state = self.state.wrapping_mul(A).wrapping_add(B);
        return ((self.state >> 11) as f64) * C;
    }
}
// -----------------------------------------------------------------------------------tests
#[cfg(test)]
pub mod test_util {
    use super::*;

    const TOL: f64 = 1e-9;

    // ---------------------------------------------------------------- helpers

    pub fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < TOL
    }

    pub fn close3(a: Vec3, b: Vec3) -> bool {
        a.iter().zip(b.iter()).all(|(x, y)| close(*x, *y))
    }

    /// A unit-mass body at rest at the given position.
    pub fn body_at(x: f64, y: f64, z: f64) -> Body {
        Body {
            pos: [x, y, z],
            ..Body::new()
        }
    }

    pub fn body_with(mass: f64, x: f64, y: f64, z: f64) -> Body {
        Body {
            mass,
            pos: [x, y, z],
            ..Body::new()
        }
    }

    /// True if `pos` lies inside the cube (`center`, `hw`) on every axis.
    pub fn inside(center: Vec3, hw: f64, pos: Vec3) -> bool {
        (0..3).all(|d| (pos[d] - center[d]).abs() <= hw + TOL)
    }

    /// Deterministic spread of `n` bodies, all strictly inside the cube [-1, 1]^3.
    pub fn scattered(n: usize) -> Vec<Body> {
        (0..n)
            .map(|i| {
                body_at(
                    ((i * 7) % 13) as f64 / 6.5 - 1.0,
                    ((i * 5) % 11) as f64 / 5.5 - 1.0,
                    ((i * 3) % 17) as f64 / 8.5 - 1.0,
                )
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::test_util::*;
    use super::*;

    #[test]
    fn gravity_points_toward_the_source() {
        let source = body_at(2.0, 0.0, 0.0);
        let g = gravity(&source, &ORIGIN);
        assert!(g[0] > 0.0);
        // divide by G: gravity here is ~1e-11, far below `close`'s absolute tolerance,
        // so comparing the raw components to zero would pass for any vector
        assert!(close(g[1] / G, 0.0));
        assert!(close(g[2] / G, 0.0));
    }

    #[test]
    fn gravity_magnitude_is_g_m_over_r_squared() {
        let source = body_with(5.0, 0.0, 3.0, 0.0);
        let g = gravity(&source, &ORIGIN);
        // divide out G so the comparison isn't dwarfed by the 1e-11 scale
        assert!(close(g[1] / G, 5.0 / 9.0), "g/G was {}", g[1] / G);
    }

    // -------------------------------------------------------------------- RNG

    fn draws(seed: u64, n: usize) -> Vec<f64> {
        let mut rng = RNG::new(seed);
        (0..n).map(|_| rng.rand()).collect()
    }

    /// Expected values come from an independent implementation (Python integers,
    /// same recurrence), so this pins the multiplier, the increment, the shift and
    /// the scale, and cannot pass by re-running the Rust arithmetic.
    #[test]
    fn rng_matches_an_independent_implementation() {
        let expected = [0.07820865487829387, 0.10169876029679303, 0.6053233226252335];
        let got = draws(0, 3);
        for (g, e) in got.iter().zip(expected) {
            assert!((g - e).abs() < 1e-15, "got {g}, expected {e}");
        }
    }

    #[test]
    fn rng_is_reproducible_from_its_seed() {
        assert_eq!(draws(42, 100), draws(42, 100));
    }

    #[test]
    fn rng_seeds_give_different_sequences() {
        assert_ne!(draws(1, 10), draws(2, 10));
    }

    #[test]
    fn rng_draws_are_in_the_half_open_unit_interval() {
        let v = draws(7, 100_000);
        assert!(v.iter().all(|&x| (0.0..1.0).contains(&x)));
    }

    /// A generator stuck on one value still passes every range check, so this
    /// asks for 1000 distinct draws.
    #[test]
    fn rng_does_not_repeat_itself() {
        let mut v = draws(7, 1000);
        v.sort_by(|a, b| a.partial_cmp(b).unwrap());
        v.dedup();
        assert_eq!(v.len(), 1000);
    }

    /// 100k draws in 10 equal bins: about 10k each, with a standard deviation of
    /// about 95, so 9.5k..10.5k is a wide margin. This is what catches a
    /// generator that is always near 0 (the integer-division bug) or lopsided.
    #[test]
    fn rng_draws_are_spread_evenly() {
        let n = 100_000;
        let mut bins = [0usize; 10];
        for x in draws(7, n) {
            bins[(x * 10.0) as usize] += 1;
        }
        for (i, count) in bins.iter().enumerate() {
            assert!((9_500..=10_500).contains(count), "bin {i} had {count}");
        }
        let mean: f64 = draws(7, n).iter().sum::<f64>() / n as f64;
        assert!((mean - 0.5).abs() < 0.01, "mean was {mean}");
    }

    #[test]
    fn gravity_has_the_right_components_in_3d() {
        // source at (1, 2, 2): distance 3, so |g| = G/9 along (1/3, 2/3, 2/3)
        let g = gravity(&body_at(1.0, 2.0, 2.0), &ORIGIN);
        assert!(close(g[0] / G, 1.0 / 27.0), "gx/G was {}", g[0] / G);
        assert!(close(g[1] / G, 2.0 / 27.0), "gy/G was {}", g[1] / G);
        assert!(close(g[2] / G, 2.0 / 27.0), "gz/G was {}", g[2] / G);
    }
}
