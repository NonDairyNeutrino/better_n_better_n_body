/// This is the docstring for Object
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Body {
    pub mass: f64,
    pub pos: [f64; 3],
    pub vel: [f64; 3],
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
    /// Calculates the Euclidean distance between an object and a position
    pub fn distance(&self, pos: &[f64; 3]) -> f64 {
        self.pos
            .iter()
            .zip(pos)
            .map(|pair| (*pair.0 - pair.1).powi(2))
            .sum::<f64>()
            .sqrt()
    }
    pub fn accelerate(&mut self, acc: &[f64; 3], dt: f64) {
        // evolve via symplectic-euler
        for (i, a) in acc.iter().enumerate() {
            self.vel[i] += a * dt;
            self.pos[i] += self.vel[i] * dt;
        }
    }
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
