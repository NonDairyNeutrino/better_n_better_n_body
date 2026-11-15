//! This is a module documentation

use std::array::from_fn;

pub const G: f64 = 6.67 * 1e-11;
pub const ORIGIN: [f64; 3] = [0.0_f64, 0.0_f64, 0.0_f64];

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

/// A collection of bodies in a cubical space with a center and half-width
#[derive(Debug)]
pub struct System<'a> {
    pub bodies: Vec<&'a Body>,
    pub center: [f64; 3],
    pub hw: f64,
}
impl<'a> System<'a> {
    pub fn new(center: [f64; 3], half_width: f64) -> System<'a> {
        System {
            bodies: Vec::new(),
            center: center,
            hw: half_width,
        }
    }
    pub fn from_bodies(bodies: &'a [Body]) -> System {
        let bods: Vec<&'a Body> = bodies.iter().collect();
        let center: [f64; 3] = find_center(bodies);
        let hw: f64 = find_half_width(bodies, center);
        System {
            bodies: bods,
            center: center,
            hw: hw,
        }
    }
    pub fn add_body(&mut self, body: &'a Body) {
        self.bodies.push(body);
    }
    pub fn add_bodies(&mut self, bodies: &'a [Body]) {
        // know how much more space is needed, so allocate first, then push each
        self.bodies.reserve(bodies.len());
        for body in bodies {
            self.bodies.push(body);
        }
    }
    pub fn join(&mut self, system: &'a System) {
        self.bodies.extend_from_slice(&system.bodies);
    }
    pub fn find_com(&self) -> [f64; 3] {
        let mut com: [f64; 3] = ORIGIN;
        if self.bodies.is_empty() {
            return com;
        }

        // find the mass-weighted total
        for body in &self.bodies {
            for dim in 0..com.len() {
                com[dim] += body.mass * body.pos[dim];
            }
        }

        // normalize by total mass
        let total_mass: f64 = self.bodies.iter().map(|b| b.mass).sum();
        for dim in 0..com.len() {
            com[dim] /= total_mass;
        }
        return com;
    }
    pub fn partition(self) -> [System<'a>; 8] {
        let mut systems: [System; 8] =
            from_fn(|i| System::new(child_center(self.center, self.hw, i), self.hw / 2.0));

        for body in self.bodies {
            let idx: usize = child_index(self.center, body.pos);
            systems[idx].add_body(body);
        }

        return systems;
    }
}

fn child_index(center: [f64; 3], pos: [f64; 3]) -> usize {
    let mut idx: usize = 0;
    if pos[0] >= center[0] {
        idx += 1;
    }
    if pos[1] >= center[1] {
        idx += 2;
    }
    if pos[2] >= center[2] {
        idx += 4;
    }
    return idx;
}
fn child_center(center: [f64; 3], half_width: f64, index: usize) -> [f64; 3] {
    let mut chi_cen: [f64; 3] = ORIGIN;
    for dim in 0..center.len() {
        let bit: i32 = ((index / 2_usize.pow(dim as u32)) % 2) as i32;
        chi_cen[dim] = center[dim] + (-1.0_f64).powi(bit + 1) * (half_width / 2.0_f64);
    }
    return chi_cen;
}

fn find_center(system: &[Body]) -> [f64; 3] {
    let mut center: [f64; 3] = ORIGIN;
    let n_bodies: f64 = system.len() as f64;
    if system.is_empty() {
        return ORIGIN;
    }

    for obj in system {
        for dim in 0..3 {
            center[dim] += obj.pos[dim];
        }
    }

    // average
    for dim in 0..center.len() {
        center[dim] /= n_bodies;
    }
    return center;
}

fn find_half_width(system: &[Body], center: [f64; 3]) -> f64 {
    let mut half_width: f64 = 0.0_f64;
    for obj in system {
        let d: f64 = obj.distance(&center);
        if d > half_width {
            half_width = d;
        }
    }
    return half_width;
}

/// Find the center and half-width of a cube that contains the system.
pub fn find_bounding_box(system: &[Body]) -> ([f64; 3], f64) {
    let center: [f64; 3] = find_center(system);
    let width: f64 = find_half_width(system, center);
    return (center, width);
}

pub fn gravity(object: &Body, pos: &[f64; 3]) -> [f64; 3] {
    let strength: f64 = G * object.mass / object.distance(pos).powi(3); // 3 is on purpose
    let mut grav: [f64; 3] = [0.0_f64, 0.0_f64, 0.0_f64];

    for dim in 0..grav.len() {
        grav[dim] = strength * (object.pos[dim] - pos[dim]);
    }

    return grav;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_is_euclidean() {
        let body = Body {
            pos: [3.0, 4.0, 0.0],
            ..Body::new()
        };
        let d = body.distance(&[0.0, 0.0, 0.0]);
        assert!((d - 5.0).abs() < 1e-9);
    }
}
