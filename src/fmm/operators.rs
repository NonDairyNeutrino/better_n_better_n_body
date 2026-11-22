use std::array::from_fn;

use crate::{
    fmm::{
        complex::Complex,
        harmonic::{SolidHarmonic, from_index, harmonics, to_index},
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

/// Create a multipole centered at a point from the multipole centered at another point
pub fn multipole_to_multipole(mp: &SolidHarmonic, shift: Vec3) -> SolidHarmonic {
    if let SolidHarmonic::Regular(_) = mp {
        panic!("M2M not defined for Regular harmonics")
    }
    let mut mp_shifted: SolidHarmonic = SolidHarmonic::irregular(mp.order());
    let harmonics_shifted: Vec<Complex> = shift_scale_harmonics(mp.order(), shift);
    mix_harmonics(mp, &mut mp_shifted, harmonics_shifted.as_slice());
    return mp_shifted;
}

fn shift_scale_harmonics(order: usize, shift: Vec3) -> Vec<Complex> {
    let shift_size: f64 = shift.map(|dim| dim.powi(2)).iter().sum::<f64>().sqrt();
    if shift_size == 0.0 {
        let mut h: Vec<Complex> = harmonics(shift, order);
        for m in h.iter_mut() {
            *m = Complex { re: 0.0, im: 0.0 }
        }
        h[0] = Complex { re: 1.0, im: 0.0 };
        h
    } else {
        harmonics(shift, order)
            .iter()
            .enumerate()
            .map(|(i, h)| *h * shift_size.powi(from_index(i).0 as i32))
            .collect()
    }
}

fn mix_harmonics(mp: &SolidHarmonic, mp_shifted: &mut SolidHarmonic, shifted_moments: &[Complex]) {
    for n in 0..=mp_shifted.order() {
        let ni: i64 = n as i64;

        for m in -ni..=ni {
            let idx: usize = to_index(ni, m);
            let slot: &mut Complex = &mut mp_shifted.coeffs_mut()[idx];
            *slot += get_mixed_moment(mp, shifted_moments, ni, m);
        }
    }
}

fn get_mixed_moment(mp: &SolidHarmonic, shifted_moments: &[Complex], n: i64, m: i64) -> Complex {
    let mut acc: Complex = Complex::zero();
    for lambda in 0..=n {
        let lb: i64 = (-(lambda as i64)).max(m - (n - lambda));
        let ub: i64 = (lambda as i64).min(m + (n - lambda));
        for mu in lb..=ub {
            let factor: f64 = get_mix_coefficient_factor(n, m, lambda, mu);
            let h: Complex = get_mix_coefficient_direction(mp, shifted_moments, n, m, lambda, mu);
            acc += h * factor;
        }
    }
    return acc;
}

fn get_mix_coefficient_factor(n: i64, m: i64, lambda: i64, mu: i64) -> f64 {
    let ep: f64 = epsilon(m, mu);
    let bin: f64 = (binomial(n + m, lambda + mu) * binomial(n - m, lambda - mu)).sqrt();
    return ep * bin;
}

fn get_mix_coefficient_direction(
    mp: &SolidHarmonic,
    shifted_moments: &[Complex],
    n: i64,
    m: i64,
    lambda: i64,
    mu: i64,
) -> Complex {
    let shifted: Complex = shifted_moments[to_index(lambda, mu)].conj();
    let unshifted: Complex = mp.get(n - lambda, m - mu);
    return shifted * unshifted;
}

// small helper because the hermonics aren't using standard convention for negative m
fn sigma(k: i64) -> f64 {
    if k < 0 && k.rem_euclid(2) == 1 {
        -1.0
    } else {
        1.0
    }
}

fn epsilon(m: i64, mu: i64) -> f64 {
    sigma(m) * sigma(mu) * sigma(m - mu)
}

fn binomial(n: i64, k: i64) -> f64 {
    if k < 0 || n < k {
        return 0.0;
    } else {
        let mut bin: f64 = 1.0;
        for i in 0..=(k - 1) {
            bin *= ((n - i) as f64) / ((i + 1) as f64);
        }
        return bin;
    }
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

    // ------------------------------------------------------ multipole_to_multipole

    /// Largest difference between corresponding slots, scaled by the slot size.
    fn max_diff(a: &SolidHarmonic, b: &SolidHarmonic) -> f64 {
        a.coeffs()
            .iter()
            .zip(b.coeffs())
            .map(|(x, y)| {
                ((x.re - y.re).abs() + (x.im - y.im).abs()) / (1.0 + y.re.abs() + y.im.abs())
            })
            .fold(0.0, f64::max)
    }

    fn assert_same(a: &SolidHarmonic, b: &SolidHarmonic, tol: f64, what: &str) {
        assert_eq!(a.order(), b.order(), "{what}: order");
        for (k, (x, y)) in a.coeffs().iter().zip(b.coeffs()).enumerate() {
            let (n, m) = from_index(k);
            assert_close(*x, *y, tol, &format!("{what}, slot (n={n}, m={m})"));
        }
    }

    /// Pascal's triangle up to row `max`: c[n][k] = C(n, k), by addition only.
    fn pascal(max: usize) -> Vec<Vec<f64>> {
        let mut c: Vec<Vec<f64>> = vec![vec![1.0]];
        for n in 1..=max {
            let prev = &c[n - 1];
            let row = (0..=n)
                .map(|k| {
                    let left = if k > 0 { prev[k - 1] } else { 0.0 };
                    let right = if k < n { prev[k] } else { 0.0 };
                    left + right
                })
                .collect();
            c.push(row);
        }
        c
    }

    /// The strong check: the moments about the old center, shifted, must equal
    /// the moments of the same bodies computed directly about the new center.
    /// Shifts along +z, -z and x hit the axes, where the direction is special.
    #[test]
    fn shifting_a_summary_equals_summarizing_about_the_new_center_directly() {
        let child = [0.4, -0.3, 0.2];
        let bodies = bodies_around(child, 0.3, 9, 5);
        let shifts: [Vec3; 4] = [
            [0.5, 0.2, -0.7],
            [0.0, 0.0, 0.6],
            [0.0, 0.0, -0.6],
            [0.8, 0.0, 0.0],
        ];
        for t in shifts {
            let parent: Vec3 = from_fn(|i| child[i] - t[i]);
            for p in [0, 1, 2, 4, 8] {
                let old = particle_to_multipole(&bodies, child, p);
                let direct = particle_to_multipole(&bodies, parent, p);
                assert_same(
                    &multipole_to_multipole(&old, t),
                    &direct,
                    1e-12,
                    &format!("shift {t:?}, p={p}"),
                );
                if p >= 2 {
                    assert!(
                        max_diff(&old, &direct) > 1e-2,
                        "shift {t:?}, p={p}: the shift changes nothing, so this proves nothing"
                    );
                }
            }
        }
    }

    #[test]
    fn a_zero_shift_returns_the_input() {
        let center = [0.2, 0.1, 0.0];
        let mp = particle_to_multipole(&bodies_around(center, 0.3, 7, 9), center, 5);
        assert!(mp.coeffs()[1..].iter().any(|c| c.re != 0.0 || c.im != 0.0));
        assert_same(
            &multipole_to_multipole(&mp, [0.0; 3]),
            &mp,
            1e-14,
            "zero shift",
        );
    }

    #[test]
    fn shifting_twice_equals_shifting_by_the_sum() {
        let center = [0.2, 0.1, 0.0];
        let mp = particle_to_multipole(&bodies_around(center, 0.3, 7, 9), center, 6);
        let (a, b): (Vec3, Vec3) = ([0.3, -0.1, 0.2], [-0.2, 0.5, 0.1]);
        let sum: Vec3 = from_fn(|i| a[i] + b[i]);
        let twice = multipole_to_multipole(&multipole_to_multipole(&mp, a), b);
        let once = multipole_to_multipole(&mp, sum);
        assert_same(&twice, &once, 1e-12, "a then b versus a+b");
        assert!(
            max_diff(&once, &mp) > 1e-2,
            "the shift must change the moments"
        );
    }

    /// Re-expressing the same mass about another center leaves the potential far
    /// away unchanged, up to the truncation of the two order-8 series.
    #[test]
    fn a_shifted_summary_gives_the_same_potential_far_away() {
        let center = [0.2, 0.1, 0.0];
        let bodies = bodies_around(center, 0.3, 7, 9);
        let t: Vec3 = [0.5, 0.4, -0.3];
        let new_center: Vec3 = from_fn(|i| center[i] - t[i]);
        let moved = multipole_to_multipole(&particle_to_multipole(&bodies, center, 8), t);
        let x: Vec3 = [4.0, -3.0, 5.0];
        let direct = direct_potential(&bodies, x);
        let got = moved.value(from_fn(|i| x[i] - new_center[i]));
        assert!(
            (got.re - direct).abs() <= 1e-6 * direct,
            "{} vs {direct}",
            got.re
        );
        assert!(got.im.abs() <= 1e-6 * direct);
    }

    #[test]
    fn shifting_keeps_the_mass_the_kind_and_the_order() {
        let center = [0.1, 0.2, 0.3];
        let bodies = bodies_around(center, 0.3, 6, 13);
        let total: f64 = bodies.iter().map(|b| b.mass).sum();
        let moved =
            multipole_to_multipole(&particle_to_multipole(&bodies, center, 4), [0.7, -0.2, 0.4]);
        assert!(matches!(moved, SolidHarmonic::Irregular(_)));
        assert_eq!(moved.order(), 4);
        assert_close(
            moved.get(0, 0),
            Complex { re: total, im: 0.0 },
            1e-12,
            "mass",
        );
    }

    /// Hand values at n = 1 for t = (0.3, -0.4, 0.5), with
    /// R_1^0(t) = t_z and R_1^1(t) = -(t_x + i t_y)/sqrt2, R_1^-1 = conj(R_1^1):
    ///   M'_1^0  = M_1^0  + t_z M_0^0
    ///   M'_1^1  = M_1^1  + conj(R_1^1(t))  M_0^0
    ///   M'_1^-1 = M_1^-1 + conj(R_1^-1(t)) M_0^0
    #[test]
    fn shifting_changes_the_first_order_moments_by_the_mass_times_the_shift() {
        let mut mp = SolidHarmonic::irregular(1);
        mp.coeffs_mut()[to_index(0, 0)] = Complex { re: 2.0, im: 0.0 };
        mp.coeffs_mut()[to_index(1, -1)] = Complex { re: 0.5, im: -0.25 };
        mp.coeffs_mut()[to_index(1, 0)] = Complex { re: -0.3, im: 0.7 };
        mp.coeffs_mut()[to_index(1, 1)] = Complex { re: 0.4, im: 0.1 };
        let moved = multipole_to_multipole(&mp, [0.3, -0.4, 0.5]);
        let s = 1.0 / 2f64.sqrt();
        assert_close(moved.get(0, 0), Complex { re: 2.0, im: 0.0 }, 1e-14, "M00");
        assert_close(moved.get(1, 0), Complex { re: 0.7, im: 0.7 }, 1e-14, "M10");
        assert_close(
            moved.get(1, 1),
            Complex {
                re: 0.4 - 2.0 * 0.3 * s,
                im: 0.1 - 2.0 * 0.4 * s,
            },
            1e-14,
            "M11",
        );
        assert_close(
            moved.get(1, -1),
            Complex {
                re: 0.5 - 2.0 * 0.3 * s,
                im: -0.25 + 2.0 * 0.4 * s,
            },
            1e-14,
            "M1-1",
        );
    }

    /// Bodies on the z-axis have only m = 0 moments, M_n = sum m z^n (for either
    /// sign of z). Shifting along z by t then follows the plain binomial rule
    /// M'_n = sum_k C(n,k) t^(n-k) M_k, which needs no spherical harmonics at all.
    #[test]
    fn shifting_along_the_z_axis_follows_the_one_dimensional_binomial_rule() {
        let p = 5;
        let (masses, zs): ([f64; 4], [f64; 4]) = ([1.0, 2.0, 0.5, 1.5], [0.3, -0.2, 0.6, -0.5]);
        let moment = |k: usize| -> f64 { (0..4).map(|j| masses[j] * zs[j].powi(k as i32)).sum() };
        let mut mp = SolidHarmonic::irregular(p);
        for n in 0..=p {
            mp.coeffs_mut()[to_index(n as i64, 0)] = Complex {
                re: moment(n),
                im: 0.0,
            };
        }
        let c = pascal(p);
        for t in [0.7, -0.45] {
            let moved = multipole_to_multipole(&mp, [0.0, 0.0, t]);
            for n in 0..=p {
                let want: f64 = (0..=n)
                    .map(|k| c[n][k] * t.powi((n - k) as i32) * moment(k))
                    .sum();
                assert_close(
                    moved.get(n as i64, 0),
                    Complex { re: want, im: 0.0 },
                    1e-12,
                    &format!("t={t}, n={n}"),
                );
            }
            for idx in 0..moved.coeffs().len() {
                let (n, m) = from_index(idx);
                if m != 0 {
                    assert_close(
                        moved.coeffs()[idx],
                        Complex { re: 0.0, im: 0.0 },
                        1e-13,
                        &format!("t={t} (n={n}, m={m})"),
                    );
                }
            }
        }
    }

    // ------------------------------------------------------------------ helpers

    #[test]
    fn sigma_is_minus_one_only_for_negative_odd_numbers() {
        for (k, want) in [
            (-4, 1.0),
            (-3, -1.0),
            (-2, 1.0),
            (-1, -1.0),
            (0, 1.0),
            (1, 1.0),
            (2, 1.0),
            (3, 1.0),
        ] {
            assert_eq!(sigma(k), want, "sigma({k})");
        }
    }

    #[test]
    fn epsilon_is_the_product_of_the_three_signs() {
        // (m, mu, expected), worked out by hand from sigma(m) sigma(mu) sigma(m-mu)
        for (m, mu, want) in [
            (1, 1, 1.0),
            (-1, -1, 1.0),
            (-1, 0, 1.0),
            (0, -1, -1.0),
            (0, -3, -1.0),
            (-2, 1, -1.0),
            (2, -1, -1.0),
            (-3, -1, 1.0),
            (0, 0, 1.0),
        ] {
            assert_eq!(epsilon(m, mu), want, "epsilon({m}, {mu})");
        }
    }

    #[test]
    fn binomial_matches_known_values_and_is_zero_outside_the_range() {
        for (n, k, want) in [
            (0, 0, 1.0),
            (5, 0, 1.0),
            (5, 5, 1.0),
            (5, 1, 5.0),
            (5, 2, 10.0),
            (6, 3, 20.0),
            (10, 4, 210.0),
            (10, 6, 210.0),
            (11, 5, 462.0),
            (40, 20, 137846528820.0),
            (3, 4, 0.0),
            (3, -1, 0.0),
        ] {
            let got = binomial(n, k);
            assert!(
                (got - want).abs() <= 1e-12 * want.max(1.0),
                "binomial({n}, {k}) = {got}, want {want}"
            );
        }
    }

    #[test]
    fn binomial_agrees_with_pascals_triangle() {
        let c = pascal(30);
        for n in 0..=30usize {
            for k in 0..=n {
                let got = binomial(n as i64, k as i64);
                assert!(
                    (got - c[n][k]).abs() <= 1e-12 * c[n][k],
                    "binomial({n}, {k})"
                );
            }
        }
    }

    // ------------------------------------------------------------ misuse

    /// A shift is defined for multipole moments only; a local expansion has a
    /// different translation rule, so a Regular input must be refused, whether
    /// or not the shift is zero.
    #[test]
    #[should_panic(expected = "not defined for Regular")]
    fn m2m_refuses_a_regular_input() {
        multipole_to_multipole(&SolidHarmonic::regular(2), [0.1, -0.2, 0.3]);
    }

    #[test]
    #[should_panic(expected = "not defined for Regular")]
    fn m2m_refuses_a_regular_input_even_for_a_zero_shift() {
        multipole_to_multipole(&SolidHarmonic::regular(2), [0.0; 3]);
    }
}
