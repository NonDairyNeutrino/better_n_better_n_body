use std::array::from_fn;

use crate::{
    fmm::{
        complex::Complex,
        harmonic::{SolidHarmonic, from_index, harmonics},
    },
    shared::{Body, Vec3},
};

pub fn particle_to_multipole<'b>(
    bodies: impl IntoIterator<Item = &'b Body>,
    center: Vec3,
    p: usize,
) -> SolidHarmonic {
    let mut mp: SolidHarmonic = SolidHarmonic::irregular(p);
    let coeffs: &mut [Complex] = mp.coeffs_mut();

    for body in bodies {
        let mass: f64 = body.mass;
        let sep: Vec3 = from_fn(|i| body.pos[i] - center[i]);
        let r: f64 = sep.iter().map(|e| e.powi(2)).sum::<f64>().sqrt();
        if r == 0.0 {
            coeffs[0] += Complex { re: mass, im: 0.0 };
            continue;
        }
        let harms: Vec<Complex> = harmonics(sep, p);
        for (i, (c, h)) in coeffs.iter_mut().zip(harms).enumerate() {
            let (n, _): (i64, i64) = from_index(i);
            *c += h.conj() * mass * r.powi(n as i32);
        }
    }
    return mp;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        fmm::test_util::{assert_close, dist3},
        shared::{RNG, System, test_util::body_with},
    };

    // ---------------------------------------------------------------- helpers

    /// `n` random bodies with positions in a cube of half-width `hw` around `center`.
    fn bodies_around(center: Vec3, hw: f64, n: usize, seed: u64) -> Vec<Body> {
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
    fn direct_potential(bodies: &[Body], x: Vec3) -> f64 {
        bodies.iter().map(|b| b.mass / b.distance(&x)).sum()
    }

    // ------------------------------------------------------------------ shape

    #[test]
    fn an_empty_leaf_summarizes_to_all_zeros_of_the_right_length() {
        for p in [0, 1, 4] {
            let mp = particle_to_multipole(&[], [1.0, 2.0, 3.0], p);
            assert!(
                matches!(mp, SolidHarmonic::Irregular(_)),
                "a summary is irregular"
            );
            assert_eq!(mp.coeffs().len(), (p + 1).pow(2));
            assert!(mp.coeffs().iter().all(|c| c.re == 0.0 && c.im == 0.0));
        }
    }

    // ---------------------------------------------------- hand-computed values

    /// On the z-axis Y_n^0 = (+-1)^n and every m != 0 vanishes, so
    /// c_n0 = mass * (+-d)^n, which can be written down without the code.
    #[test]
    fn a_body_on_the_z_axis_gives_pure_m_zero_coefficients() {
        let center = [1.0, -1.0, 0.5];
        for dz in [2.0, -2.0] {
            let b = body_with(3.0, center[0], center[1], center[2] + dz);
            let mp = particle_to_multipole(&[b], center, 3);
            for idx in 0..mp.coeffs().len() {
                let (n, m) = from_index(idx);
                let want = if m == 0 { 3.0 * dz.powi(n as i32) } else { 0.0 };
                assert_close(
                    mp.coeffs()[idx],
                    Complex { re: want, im: 0.0 },
                    1e-12,
                    &format!("dz={dz} (n={n}, m={m})"),
                );
            }
        }
    }

    /// d = (1, 2, 2): R_1^0 = z, R_1^1 = -(x+iy)/sqrt2, R_1^-1 = conj(R_1^1).
    /// The summary stores conj(R_n^m).
    #[test]
    fn an_off_axis_body_gives_the_conjugated_solid_harmonics() {
        let mp = particle_to_multipole(&[body_with(1.0, 1.0, 2.0, 2.0)], [0.0; 3], 1);
        let s = 1.0 / 2f64.sqrt();
        assert_close(mp.get(0, 0), Complex { re: 1.0, im: 0.0 }, 1e-14, "c00");
        assert_close(mp.get(1, 0), Complex { re: 2.0, im: 0.0 }, 1e-14, "c10");
        assert_close(
            mp.get(1, 1),
            Complex {
                re: -s,
                im: 2.0 * s,
            },
            1e-14,
            "c11 = conj(R11)",
        );
        assert_close(
            mp.get(1, -1),
            Complex {
                re: -s,
                im: -2.0 * s,
            },
            1e-14,
            "c1-1 = conj(R1-1)",
        );
    }

    // ------------------------------------------------------- the physics check

    /// The summary evaluated far away must reproduce the direct sum
    /// sum_j m_j / |x - x_j|, with error below the geometric tail
    /// sum_j m_j rho_j^(p+1) / (R - r_j), and clearly better at high order than at p = 0.
    #[test]
    fn a_leaf_summary_matches_the_direct_potential_far_away_within_the_tail_bound() {
        let center = [1.0, -2.0, 0.5];
        let bodies = bodies_around(center, 0.4, 12, 1);
        let far = 3.0;
        let r_max = bodies
            .iter()
            .map(|b| b.distance(&center))
            .fold(0.0, f64::max);
        assert!(r_max / far < 0.3, "fixture is meant to be well separated");
        let dirs: [Vec3; 5] = [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, -1.0],
            [1.0, 1.0, 1.0],
            [-1.0, 0.5, -0.3],
        ];
        let mut worst: Vec<f64> = Vec::new();
        for p in [0usize, 1, 2, 4, 6, 8, 10] {
            let mp = particle_to_multipole(&bodies, center, p);
            let mut worst_p: f64 = 0.0;
            for d in dirs {
                let norm = dist3(d, [0.0; 3]);
                let x: Vec3 = from_fn(|i| center[i] + far * d[i] / norm);
                let v = mp.value(from_fn(|i| x[i] - center[i]));
                let err = (v.re - direct_potential(&bodies, x)).abs();
                let bound: f64 = bodies
                    .iter()
                    .map(|b| {
                        let rj = b.distance(&center);
                        b.mass * (rj / far).powi(p as i32 + 1) / (far - rj)
                    })
                    .sum();
                assert!(
                    v.im.abs() < 1e-12,
                    "p={p}: potential must be real, im={}",
                    v.im
                );
                assert!(
                    err <= bound * (1.0 + 1e-9) + 1e-13,
                    "p={p}, dir={d:?}: err {err} > tail bound {bound}"
                );
                worst_p = worst_p.max(err);
            }
            worst.push(worst_p);
        }
        assert!(
            worst[0] > 1e-3,
            "p=0 must be visibly wrong, got {}",
            worst[0]
        );
        assert!(
            *worst.last().unwrap() < 1e-5,
            "p=10 worst err {}",
            worst.last().unwrap()
        );
    }

    /// About the center of mass the first moment vanishes, so n = 1 is zero.
    /// Moving the center away must make it nonzero, otherwise this proves nothing.
    #[test]
    fn about_the_center_of_mass_the_dipole_terms_vanish() {
        let bodies = bodies_around([0.0, 0.0, 0.0], 0.5, 10, 7);
        let com = System::from_bodies(&bodies).find_com();
        let mp = particle_to_multipole(&bodies, com, 2);
        for m in -1..=1 {
            assert_close(
                mp.get(1, m),
                Complex { re: 0.0, im: 0.0 },
                1e-12,
                &format!("dipole m={m} at the COM"),
            );
        }
        let shifted = particle_to_multipole(&bodies, [com[0] + 0.3, com[1], com[2]], 2);
        let size: f64 = (-1..=1)
            .map(|m| shifted.get(1, m).re.abs() + shifted.get(1, m).im.abs())
            .sum();
        assert!(size > 1e-2, "dipole about a shifted center is only {size}");
    }

    // ----------------------------------------------------------- symmetries

    #[test]
    fn summaries_add_over_bodies_and_scale_with_mass() {
        let center = [0.0, 0.5, 0.0];
        let a = bodies_around(center, 0.4, 5, 21);
        let b = bodies_around(center, 0.4, 4, 22);
        let all: Vec<Body> = a.iter().chain(&b).copied().collect();
        let (sa, sb, sall) = (
            particle_to_multipole(&a, center, 4),
            particle_to_multipole(&b, center, 4),
            particle_to_multipole(&all, center, 4),
        );
        for idx in 0..sall.coeffs().len() {
            assert_close(
                sall.coeffs()[idx],
                sa.coeffs()[idx] + sb.coeffs()[idx],
                1e-12,
                &format!("sum, slot {idx}"),
            );
        }
        let heavy: Vec<Body> = a
            .iter()
            .map(|x| Body {
                mass: 2.5 * x.mass,
                ..*x
            })
            .collect();
        let sh = particle_to_multipole(&heavy, center, 4);
        for idx in 0..sh.coeffs().len() {
            assert_close(
                sh.coeffs()[idx],
                sa.coeffs()[idx] * 2.5,
                1e-12,
                &format!("mass x2.5, slot {idx}"),
            );
        }
    }

    // ----------------------------------------------------------- edge case

    /// A body exactly at the center has no direction; only the n = 0 term
    /// (its mass) can contribute, and nothing may become NaN.
    #[test]
    fn a_body_exactly_at_the_center_contributes_only_its_mass() {
        let center = [0.7, -0.2, 0.1];
        let mp = particle_to_multipole(
            &[body_with(2.0, center[0], center[1], center[2])],
            center,
            3,
        );
        for idx in 0..mp.coeffs().len() {
            let want = if idx == 0 { 2.0 } else { 0.0 };
            assert_close(
                mp.coeffs()[idx],
                Complex { re: want, im: 0.0 },
                1e-15,
                &format!("slot {idx}"),
            );
        }
    }

    /// The center body adds its mass to slot 0 and nothing above; it must not
    /// erase what the other bodies contributed. Hand values: c00 = 1 + 2 + 4 and
    /// c10 = sum of m z = 1 * 0.4 + 4 * (-0.6) (the center body has z = 0).
    #[test]
    fn a_body_at_the_center_among_others_keeps_their_contributions() {
        let center = [0.0; 3];
        let bodies = [
            body_with(1.0, 0.5, 0.2, 0.4),
            body_with(2.0, 0.0, 0.0, 0.0),
            body_with(4.0, -0.1, 0.3, -0.6),
        ];
        let mp = particle_to_multipole(&bodies, center, 2);
        assert_close(mp.get(0, 0), Complex { re: 7.0, im: 0.0 }, 1e-12, "c00");
        assert_close(mp.get(1, 0), Complex { re: -2.0, im: 0.0 }, 1e-12, "c10");
        let swapped = [bodies[1], bodies[0], bodies[2]];
        let mp2 = particle_to_multipole(&swapped, center, 2);
        assert_close(
            mp2.get(0, 0),
            Complex { re: 7.0, im: 0.0 },
            1e-12,
            "c00, center body first",
        );
        assert_close(
            mp2.get(1, 0),
            Complex { re: -2.0, im: 0.0 },
            1e-12,
            "c10, center body first",
        );
    }
}
