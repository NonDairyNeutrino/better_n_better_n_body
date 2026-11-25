use std::array::from_fn;

use crate::{
    fmm::{
        complex::Complex,
        harmonic::{SolidHarmonic, to_index},
        operators::local_to_local,
    },
    shared::{Body, G, Vec3},
};

pub fn local_to_particle(local: &SolidHarmonic, center: Vec3, bodies: &[&Body], out: &mut [Vec3]) {
    assert_eq!(bodies.len(), out.len());

    if let SolidHarmonic::Irregular(_) = local {
        panic!("L2P requires regular harmonic as input");
    }

    for (body, slot) in bodies.iter().zip(out.iter_mut()) {
        let shift: Vec3 = from_fn(|i| body.pos[i] - center[i]);
        let shifted: SolidHarmonic = local_to_local(&local, shift);
        let coeffs: &[Complex] = shifted.coeffs();

        let a: Complex = coeffs[to_index(1, 0)];
        let b: Complex = coeffs[to_index(1, 1)];
        let c: Complex = coeffs[to_index(1, -1)];

        let grad_x = (-(b + c) * 2_f64.powf(-0.5)).re;
        let grad_y = (b.im - c.im/* Re(i(c - b)) */) * 2_f64.powf(-0.5);
        let grad_z = a.re;

        slot[0] += G * grad_x;
        slot[1] += G * grad_y;
        slot[2] += G * grad_z;
    }
}
