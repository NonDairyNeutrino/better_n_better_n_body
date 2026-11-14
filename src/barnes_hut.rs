use crate::shared::Object;

fn find_center(system: &[Object]) -> [f64; 3] {
    let mut center: [f64; 3] = [0.0_f64, 0.0_f64, 0.0_f64];
    let n_bodies: f64 = system.len() as f64;

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

fn find_half_width(system: &[Object]) -> f64 {
    let mut half_width: f64 = 0.0_f64;
    for obj in system {
        let d: f64 = obj.distance(&[0.0_f64, 0.0_f64, 0.0_f64]);
        if d > half_width {
            half_width = d;
        }
    }
    return half_width;
}

/// Find the center and half-width of a cube that contains the system
pub fn find_bounding_box(system: &[Object]) -> ([f64; 3], f64) {
    let width: f64 = find_half_width(system);
    let center: [f64; 3] = find_center(system);
    return (center, width);
}
