use std::array;

/// This is the docstring for Object
use crate::shared::{RNG, Vec3};
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Body {
    pub mass: f64,
    pub pos: Vec3,
    pub vel: Vec3,
}
impl Body {
    /// Creates a new physical Object
    pub fn new() -> Body {
        return Body {
            pos: [0f64, 0f64, 0f64],
            vel: [0f64, 0f64, 0f64],
            mass: 1f64,
        };
    }
    pub fn rand(rng: &mut RNG, half_width: f64) -> Body {
        let mass: f64 = rand_positive(rng, half_width);
        let pos: Vec3 = rand_vec3(rng, half_width);
        let vel: Vec3 = rand_vec3(rng, half_width);
        return Body { mass, pos, vel };
    }
    /// Calculates the Euclidean distance between an object and a position
    pub fn distance(&self, pos: &Vec3) -> f64 {
        self.pos
            .iter()
            .zip(pos)
            .map(|pair| (*pair.0 - pair.1).powi(2))
            .sum::<f64>()
            .sqrt()
    }
    pub fn accelerate(&mut self, acc: &Vec3, dt: f64) {
        // evolve via symplectic-euler
        for (i, a) in acc.iter().enumerate() {
            self.vel[i] += a * dt;
            self.pos[i] += self.vel[i] * dt;
        }
    }
}

fn rand_positive(rng: &mut RNG, half_width: f64) -> f64 {
    half_width * (1.0_f64 - rng.rand())
}
fn rand_vec3(rng: &mut RNG, half_width: f64) -> Vec3 {
    let mut r: Vec3 = array::from_fn(|_| rng.rand());
    for ri in r.iter_mut() {
        *ri = half_width * (2.0_f64 * *ri - 1.0_f64);
    }
    return r;
}

#[cfg(test)]
mod tests {
    use super::super::test_util::*;
    use super::*;
    use crate::shared::ORIGIN;

    const TOL: f64 = 1e-9;

    #[test]
    fn distance_is_euclidean() {
        let body = Body {
            pos: [3.0, 4.0, 0.0],
            ..Body::new()
        };
        let d = body.distance(&[0.0, 0.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-9);
    }

    #[test]
    fn distance_uses_all_three_axes() {
        // (1, 2, 2) is exactly 3 from the origin
        let body = body_at(1.0, 2.0, 2.0);
        assert!(close(body.distance(&ORIGIN), 3.0));
        // and the offset is measured from the given point, not from the origin
        assert!(close(body.distance(&[1.0, 2.0, 0.0]), 2.0));
    }

    #[test]
    fn accelerate_without_acceleration_moves_at_constant_velocity() {
        let mut body = Body {
            vel: [1.0, 2.0, 3.0],
            ..Body::new()
        };
        body.accelerate(&[0.0, 0.0, 0.0], 2.0);
        assert!(close3(body.pos, [2.0, 4.0, 6.0]));
        assert!(close3(body.vel, [1.0, 2.0, 3.0]));
    }

    // -------------------------------------------------------------- Body::rand

    fn random_bodies(seed: u64, n: usize, size: f64) -> Vec<Body> {
        let mut rng = RNG::new(seed);
        (0..n).map(|_| Body::rand(&mut rng, size)).collect()
    }

    #[test]
    fn rand_bodies_have_positive_mass_and_stay_inside_the_cube() {
        // a size that is not 1, so a range that ignores `size` shows up
        let size = 2.5;
        for b in random_bodies(3, 5_000, size) {
            assert!(b.mass > 0.0 && b.mass <= size, "mass {}", b.mass);
            for d in 0..3 {
                assert!((-size..size).contains(&b.pos[d]), "pos {:?}", b.pos);
                assert!((-size..size).contains(&b.vel[d]), "vel {:?}", b.vel);
            }
        }
    }

    /// Staying inside the range is easy; filling it is the real check. A body
    /// generator that returned tiny values would pass the test above.
    #[test]
    fn rand_bodies_fill_the_whole_range_and_are_centered() {
        let size = 2.5;
        let n = 10_000;
        let bodies = random_bodies(3, n, size);
        let max_pos = bodies.iter().flat_map(|b| b.pos).fold(f64::MIN, f64::max);
        let min_pos = bodies.iter().flat_map(|b| b.pos).fold(f64::MAX, f64::min);
        assert!(max_pos > 0.99 * size && min_pos < -0.99 * size, "{min_pos}..{max_pos}");

        // symmetric about 0 (std of the mean is about 0.58*size/sqrt(n) = 0.015;
        // 0.1 is a wide margin), and mass is uniform on (0, size] so its mean is size/2
        let mean = |f: &dyn Fn(&Body) -> f64| bodies.iter().map(f).sum::<f64>() / n as f64;
        for d in 0..3 {
            assert!(mean(&|b| b.pos[d]).abs() < 0.1, "pos mean, axis {d}");
            assert!(mean(&|b| b.vel[d]).abs() < 0.1, "vel mean, axis {d}");
        }
        assert!((mean(&|b| b.mass) - size / 2.0).abs() < 0.1);
    }

    #[test]
    fn rand_scales_with_the_cube_size() {
        let small = random_bodies(3, 1_000, 0.01);
        assert!(small.iter().all(|b| b.mass <= 0.01 && b.pos.iter().all(|x| x.abs() < 0.01)));
        let big = random_bodies(3, 1_000, 1000.0);
        assert!(big.iter().any(|b| b.pos.iter().any(|x| x.abs() > 100.0)));
    }

    #[test]
    fn rand_is_reproducible_from_the_seed() {
        assert_eq!(random_bodies(11, 50, 1.0), random_bodies(11, 50, 1.0));
        assert_ne!(random_bodies(11, 50, 1.0), random_bodies(12, 50, 1.0));
    }

    /// Each field must come from its own draws. If, say, the position reused one
    /// number for all three axes, or velocity copied position, the cloud would
    /// collapse onto a line or a plane.
    #[test]
    fn rand_fields_are_independent_draws() {
        for b in random_bodies(5, 200, 1.0) {
            assert!(b.pos[0] != b.pos[1] && b.pos[1] != b.pos[2] && b.pos[0] != b.pos[2]);
            assert_ne!(b.pos, b.vel);
        }
    }

    /// The generator has to hand the simulations inputs that honour their contract:
    /// no two bodies at one position and none closer than the 1e-10 resolution
    /// limit. Checked on the same fixture size the benchmark will use.
    #[test]
    fn rand_bodies_respect_the_input_contract() {
        let bodies = random_bodies(9, 2_000, 1.0);
        let mut closest = f64::MAX;
        for i in 0..bodies.len() {
            for j in (i + 1)..bodies.len() {
                closest = closest.min(bodies[i].distance(&bodies[j].pos));
            }
        }
        assert!(closest > 1e-10, "closest pair was {closest:e} apart");
    }

    /// The comment in `accelerate` says symplectic Euler: update the velocity
    /// first, then move by the NEW velocity (x += v_new * dt). dt is deliberately
    /// not 1, so a wrong power of dt shows up, and every axis is exercised.
    #[test]
    fn accelerate_position_follows_symplectic_euler() {
        let mut body = Body::new();
        body.accelerate(&[2.0, 4.0, 6.0], 0.5);
        assert!(close3(body.vel, [1.0, 2.0, 3.0]));
        assert!(close3(body.pos, [0.5, 1.0, 1.5]), "pos was {:?}", body.pos);
    }
}
