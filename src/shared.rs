//! This is a module documentation

pub const G: f64 = 6.67 * 1e-11;

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
            self.pos[i] += self.vel[i] * dt + 0.5_f64 * a * dt.powi(2);
        }
    }
}

pub fn gravity(object: &Body, pos: &[f64; 3]) -> [f64; 3] {
    let strength: f64 = G * object.mass / object.distance(pos).powi(3); // 3 is on purpose
    let mut grav: [f64; 3] = [0.0_f64, 0.0_f64, 0.0_f64];

    for dim in 0..grav.len() {
        grav[dim] = strength * (object.pos[dim] - pos[dim]);
    }

    return grav;
}
