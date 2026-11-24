use crate::{
    fmm::{
        complex::Complex,
        harmonic::{SolidHarmonic, from_index, harmonics, to_index},
        operators::m2m::binomial,
    },
    shared::Vec3,
};

/// Turn the mass distribution of a source cell into the external potential at a target cell.
pub fn multipole_to_local(mp: &SolidHarmonic, shift: Vec3) -> SolidHarmonic {
    if let SolidHarmonic::Regular(_) = mp {
        panic!("M2L requires irregular harmonic as input")
    }

    let table: Vec<Complex> = irregular_table(mp.order(), shift);
    let mut local_potential: SolidHarmonic = SolidHarmonic::regular(mp.order());
    convert_moments(mp, &mut local_potential, &table);
    return local_potential;
}

fn irregular_table(order: usize, shift: Vec3) -> Vec<Complex> {
    let r: f64 = shift.map(|e| e.powi(2)).iter().sum::<f64>().sqrt();
    if r == 0.0 {
        panic!("Cannot use zero shift")
    }
    let harms: Vec<Complex> = harmonics(shift, 2 * order);
    return harms
        .iter()
        .enumerate()
        .map(|(i, e)| *e * r.powi(-(from_index(i).0 + 1) as i32))
        .collect();
}

fn convert_moments(mp: &SolidHarmonic, local: &mut SolidHarmonic, table: &[Complex]) {
    for j in 0..=local.order() {
        let ji: i64 = j as i64;

        for k in -ji..=ji {
            let idx: usize = to_index(ji, k);
            let slot: &mut Complex = &mut local.coeffs_mut()[idx];
            *slot += get_local_coefficient(mp, table, ji, k);
        }
    }
}

fn get_local_coefficient(mp: &SolidHarmonic, table: &[Complex], j: i64, k: i64) -> Complex {
    let mut acc: Complex = Complex::zero();
    for n in 0..=(mp.order() as i64) {
        for m in -n..=n {
            let factor: f64 = conversion_sign(j, k, m) * get_conversion_factor(n, m, j, k);
            let h: Complex = get_conversion_direction(mp, table, n, m, j, k);
            acc += h * factor;
        }
    }
    return acc;
}

fn get_conversion_direction(
    mp: &SolidHarmonic,
    table: &[Complex],
    n: i64,
    m: i64,
    j: i64,
    k: i64,
) -> Complex {
    let tabulated: Complex = table[to_index(j + n, m - k)];
    let moment: Complex = mp.get(n, m);
    return tabulated * moment;
}

fn get_conversion_factor(n: i64, m: i64, j: i64, k: i64) -> f64 {
    let bin1: f64 = binomial(n + j + m - k, n + m);
    let bin2: f64 = binomial(n + j - m + k, n - m);
    return (bin1 * bin2).sqrt();
}

fn conversion_sign(j: i64, k: i64, m: i64) -> f64 {
    let signj: f64 = (-1.0_f64).powi(j as i32);
    let signkm: f64 = if k.signum() == m.signum() {
        (-1.0_f64).powi(k.abs().min(m.abs()) as i32)
    } else {
        1.0_f64
    };
    return signj * signkm;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fmm::{
        operators::p2m::particle_to_multipole,
        test_util::{
            assert_close, bodies_around, direct_potential, distinct, distinct_moments, table_of,
        },
    };

    // ---------------------------------------------------------- irregular_table

    /// Degree up to 2p, so (2p + 1)^2 entries, whatever the shift.
    #[test]
    fn the_irregular_table_reaches_degree_two_p() {
        for p in 0..=4 {
            let table = irregular_table(p, [0.3, -0.4, 0.5]);
            assert_eq!(table.len(), (2 * p + 1) * (2 * p + 1), "order {p}");
        }
    }

    /// Hand values up to degree 2 for T = (0.3, -0.4, 0.5), with w = x + iy and
    /// R_N^M = r^N Y_N^M the regular harmonics (written out by hand):
    ///   R_0^0 = 1, R_1^0 = z, R_1^1 = -w/sqrt2, R_2^0 = (3z^2 - r^2)/2,
    ///   R_2^1 = -3zw/sqrt6, R_2^2 = 3w^2/sqrt24, and R_N^-M = conj(R_N^M).
    /// Since S_N^M = Y_N^M / r^(N+1) = R_N^M / r^(2N+1), each entry is R over that.
    #[test]
    fn the_irregular_table_matches_hand_values_up_to_degree_two() {
        let (x, y, z): (f64, f64, f64) = (0.3, -0.4, 0.5);
        let r2 = x * x + y * y + z * z;
        let r = r2.sqrt();
        let w = Complex { re: x, im: y };
        let w2 = Complex {
            re: x * x - y * y,
            im: 2.0 * x * y,
        };
        let one = |re: f64| Complex { re, im: 0.0 };
        let r11 = w * (-1.0 / 2f64.sqrt());
        let r21 = w * (-3.0 * z / 6f64.sqrt());
        let r22 = w2 * (3.0 / 24f64.sqrt());
        let want = [
            (0, 0, one(1.0) * r.powi(-1)),
            (1, 0, one(z) * r.powi(-3)),
            (1, 1, r11 * r.powi(-3)),
            (1, -1, r11.conj() * r.powi(-3)),
            (2, 0, one((3.0 * z * z - r2) / 2.0) * r.powi(-5)),
            (2, 1, r21 * r.powi(-5)),
            (2, -1, r21.conj() * r.powi(-5)),
            (2, 2, r22 * r.powi(-5)),
            (2, -2, r22.conj() * r.powi(-5)),
        ];
        let table = irregular_table(2, [x, y, z]);
        for (n, m, w) in want {
            assert_close(table[to_index(n, m)], w, 1e-14, &format!("S_{n}^{m}"));
        }
    }

    /// On the z-axis only M = 0 survives. S_N^0 = P_N(+-1) / |z|^(N+1), and
    /// P_N(-1) = (-1)^N, so the sign of z shows up as (-1)^N.
    #[test]
    fn on_the_z_axis_the_irregular_table_is_the_inverse_powers_of_z() {
        for z in [0.7_f64, -0.7] {
            let table = irregular_table(3, [0.0, 0.0, z]);
            for n in 0..=6i64 {
                for m in -n..=n {
                    let want = if m == 0 {
                        let sign = if z < 0.0 && n % 2 == 1 { -1.0 } else { 1.0 };
                        sign / z.abs().powi(n as i32 + 1)
                    } else {
                        0.0
                    };
                    assert_close(
                        table[to_index(n, m)],
                        Complex { re: want, im: 0.0 },
                        1e-12 * want.abs().max(1.0),
                        &format!("z={z} S_{n}^{m}"),
                    );
                }
            }
        }
    }

    /// S_N^M is homogeneous of degree -(N+1): scaling the shift by s scales
    /// degree N by s^-(N+1).
    #[test]
    fn scaling_the_shift_scales_degree_n_by_s_to_the_minus_n_plus_one() {
        let t = [0.3, -0.4, 0.5];
        let s = 2.5;
        let base = irregular_table(3, t);
        let scaled = irregular_table(3, t.map(|c| s * c));
        let mut compared_above_degree_zero = false;
        for k in 0..base.len() {
            let n = from_index(k).0 as i32;
            compared_above_degree_zero |= n > 0 && base[k].re.abs() + base[k].im.abs() > 1e-3;
            assert_close(
                scaled[k],
                base[k] * s.powi(-(n + 1)),
                1e-12 * (1.0 + base[k].re.abs() + base[k].im.abs()),
                &format!("slot {k}"),
            );
        }
        assert!(
            compared_above_degree_zero,
            "nothing above degree 0 compared"
        );
    }

    /// A lower-order table is the start of a higher-order one, since the
    /// degrees below do not depend on how far the table goes.
    #[test]
    fn a_lower_order_table_is_a_prefix_of_a_higher_order_one() {
        let t = [0.3, -0.4, 0.5];
        let big = irregular_table(3, t);
        for p in 0..3 {
            let small = irregular_table(p, t);
            for (k, c) in small.iter().enumerate() {
                assert_close(*c, big[k], 1e-14, &format!("p={p} slot {k}"));
            }
        }
    }

    /// Every S_N^M is a harmonic function away from the origin: its Laplacian
    /// vanishes. A wrong radial power breaks this, so it checks the exponent
    /// without leaning on the formula for the scale.
    #[test]
    fn every_irregular_entry_satisfies_laplaces_equation() {
        let c = [0.7, -0.4, 0.9];
        let h = 1e-3;
        let at = |d: [f64; 3]| irregular_table(2, [c[0] + d[0], c[1] + d[1], c[2] + d[2]]);
        let centre = at([0.0; 3]);
        let mut neighbours = vec![];
        for axis in 0..3 {
            for sign in [-1.0, 1.0] {
                let mut d = [0.0; 3];
                d[axis] = sign * h;
                neighbours.push(at(d));
            }
        }
        let mut largest = 0.0f64;
        for k in 0..centre.len() {
            let lap_re =
                (neighbours.iter().map(|t| t[k].re).sum::<f64>() - 6.0 * centre[k].re) / (h * h);
            let lap_im =
                (neighbours.iter().map(|t| t[k].im).sum::<f64>() - 6.0 * centre[k].im) / (h * h);
            assert!(
                lap_re.abs() < 1e-4 && lap_im.abs() < 1e-4,
                "slot {k}: Laplacian ({lap_re}, {lap_im})"
            );
            largest = largest.max(centre[k].re.abs()).max(centre[k].im.abs());
        }
        assert!(
            largest > 0.1,
            "entries too small for the check to mean anything"
        );
    }

    #[test]
    #[should_panic(expected = "zero shift")]
    fn the_irregular_table_refuses_a_zero_shift() {
        irregular_table(2, [0.0; 3]);
    }

    // ----------------------------------------------------------- conversion_sign

    /// (-1)^j times a min-based sign that applies only when k and m are both
    /// positive or both negative.
    #[test]
    fn the_sign_matches_hand_values() {
        // (j, k, m) -> sign
        let cases: [((i64, i64, i64), f64); 10] = [
            ((0, 0, 0), 1.0),
            ((1, 0, 0), -1.0), // (-1)^1
            ((2, 0, 3), 1.0),  // k = 0: only (-1)^j
            ((0, 1, 1), -1.0), // same sign, min = 1
            ((0, 1, -1), 1.0), // opposite signs
            ((0, -1, 1), 1.0),
            ((2, 2, 1), -1.0),   // (+1) * (-1)^min(2,1)
            ((1, -2, -3), -1.0), // (-1) * (-1)^2
            ((1, -1, -1), 1.0),  // (-1) * (-1)^1
            ((3, 2, 2), -1.0),   // (-1) * (-1)^2
        ];
        for ((j, k, m), want) in cases {
            assert_eq!(conversion_sign(j, k, m), want, "(j,k,m)=({j},{k},{m})");
        }
    }

    // ------------------------------------------------------ get_conversion_factor

    #[test]
    fn the_factor_matches_hand_values() {
        // (n, m, j, k) -> sqrt(C(n+j+m-k, n+m) * C(n+j-m+k, n-m))
        let cases: [((i64, i64, i64, i64), f64); 6] = [
            ((0, 0, 0, 0), 1.0),
            ((1, 0, 1, 0), 2.0), // sqrt(C(2,1) C(2,1))
            ((1, 1, 0, 0), 1.0), // sqrt(C(2,2) C(0,0))
            ((1, 0, 0, 0), 1.0),
            ((1, 1, 1, 0), 3f64.sqrt()),  // sqrt(C(3,2) C(1,0))
            ((1, -1, 1, 0), 3f64.sqrt()), // sqrt(C(1,0) C(3,2))
        ];
        for ((n, m, j, k), want) in cases {
            let got = get_conversion_factor(n, m, j, k);
            assert!(
                (got - want).abs() < 1e-14,
                "(n,m,j,k)=({n},{m},{j},{k}): got {got}, want {want}"
            );
        }
    }

    /// The same weight written with factorials instead of binomials, checked over
    /// every (n, m, j, k) up to degree 5: sqrt((n+j+m-k)!(n+j-m+k)! /
    /// ((n+m)!(n-m)!(j-k)!(j+k)!)).
    #[test]
    fn the_factor_equals_the_factorial_formula() {
        fn fact(x: i64) -> f64 {
            (1..=x).map(|v| v as f64).product()
        }
        let mut checked = 0;
        for n in 0..=5 {
            for m in -n..=n {
                for j in 0..=5 {
                    for k in -j..=j {
                        let want = (fact(n + j + m - k) * fact(n + j - m + k)
                            / (fact(n + m) * fact(n - m) * fact(j - k) * fact(j + k)))
                        .sqrt();
                        let got = get_conversion_factor(n, m, j, k);
                        assert!(
                            (got - want).abs() <= 1e-12 * want,
                            "(n,m,j,k)=({n},{m},{j},{k}): got {got}, want {want}"
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert_eq!(checked, 36 * 36);
    }

    // ------------------------------------------------- get_conversion_direction

    /// n=1, m=1, j=1, k=0 pairs table slot (j+n, m-k) = (2, 1) = 7 with moment
    /// slot (n, m) = (1, 1) = 3. The slots are written out as numbers, and there
    /// is no conjugate.
    #[test]
    fn the_direction_pairs_the_right_table_entry_with_the_right_moment() {
        let mp = distinct_moments(2, 1.0);
        let table = table_of(4, 4.0);
        let got = get_conversion_direction(&mp, &table, 1, 1, 1, 0);
        assert_close(got, table[7] * mp.coeffs()[3], 1e-14, "direction");
        // k shifts the table's azimuthal index: (j+n, m-k) = (2, 0) = 6
        let got = get_conversion_direction(&mp, &table, 1, 1, 1, 1);
        assert_close(got, table[6] * mp.coeffs()[3], 1e-14, "direction, k = 1");
        // the conjugate would give a different number
        let wrong = table[7].conj() * mp.coeffs()[3];
        let right = table[7] * mp.coeffs()[3];
        assert!(
            (wrong.im - right.im).abs() > 1e-6,
            "conjugate not detectable"
        );
    }

    // -------------------------------------------------------- get_local_coefficient

    /// j = 0, k = 0: sign +1 and every weight is 1, so
    /// L_0^0 = sum over (n, m) of M_n^m * S_n^m(T). For order 1 that is
    /// slots 0..=3 of the moments against slots 0..=3 of the table.
    #[test]
    fn the_coefficient_at_degree_zero_is_the_moments_dotted_with_the_table() {
        let mp = distinct_moments(1, 1.0);
        let table = table_of(2, 4.0);
        let m = mp.coeffs();
        let want = m[0] * table[0] + m[1] * table[1] + m[2] * table[2] + m[3] * table[3];
        assert_close(
            get_local_coefficient(&mp, &table, 0, 0),
            want,
            1e-12,
            "L_0^0",
        );
    }

    /// j = 1, k = 0, order 1: the sign is -1, and the terms are
    ///   n=0:        weight 1,   table (1, 0)  = slot 2
    ///   n=1, m=-1:  sqrt(3),    table (2, -1) = slot 5
    ///   n=1, m=0:   2,          table (2, 0)  = slot 6
    ///   n=1, m=1:   sqrt(3),    table (2, 1)  = slot 7
    /// against moment slots 0, 1, 2, 3.
    #[test]
    fn the_coefficient_at_degree_one_matches_the_hand_sum() {
        let mp = distinct_moments(1, 1.0);
        let t = table_of(2, 4.0);
        let m = mp.coeffs();
        let want = (m[0] * t[2]
            + m[1] * t[5] * 3f64.sqrt()
            + m[2] * t[6] * 2.0
            + m[3] * t[7] * 3f64.sqrt())
            * -1.0;
        assert_close(get_local_coefficient(&mp, &t, 1, 0), want, 1e-12, "L_1^0");
    }

    // ------------------------------------------------------------ convert_moments

    #[test]
    fn converting_fills_every_slot_with_that_slots_coefficient() {
        for p in 0..=3 {
            let mp = distinct_moments(p, 1.0);
            let table = table_of(2 * p, 4.0);
            let mut local = SolidHarmonic::regular(p);
            convert_moments(&mp, &mut local, &table);
            assert!(matches!(local, SolidHarmonic::Regular(_)));
            for (idx, got) in local.coeffs().iter().enumerate() {
                let (j, k) = from_index(idx);
                assert_close(
                    *got,
                    get_local_coefficient(&mp, &table, j, k),
                    1e-12,
                    &format!("p={p} slot {idx}"),
                );
                assert!(got.re.abs() + got.im.abs() > 0.0, "slot {idx} left empty");
            }
        }
    }

    // ----------------------------------------------------- multipole_to_local

    #[test]
    fn the_result_is_regular_and_has_the_order_of_the_input() {
        for p in 0..=3 {
            let local = multipole_to_local(&distinct_moments(p, 1.0), [3.0, -2.0, 2.5]);
            assert!(matches!(local, SolidHarmonic::Regular(_)));
            assert_eq!(local.order(), p);
        }
    }

    /// With only the mass (p = 0) the result is the potential of a point mass at
    /// the target center: mass / |T|.
    #[test]
    fn at_order_zero_the_result_is_the_mass_over_the_distance() {
        let mp = distinct_moments(0, 3.0);
        let local = multipole_to_local(&mp, [3.0, 4.0, 0.0]);
        let want = mp.coeffs()[0] * (1.0 / 5.0);
        assert_close(local.coeffs()[0], want, 1e-14, "L_0^0");
    }

    /// The map from moments to coefficients is linear: converting a sum of
    /// moments gives the sum of what each converts to.
    #[test]
    fn converting_a_sum_of_moments_gives_the_sum_of_the_conversions() {
        let t = [3.1, -2.2, 2.6];
        let a = distinct_moments(4, 1.0);
        let b = distinct_moments(4, -2.0);
        let mut sum = a.clone();
        sum += &b;
        let mut want = multipole_to_local(&a, t);
        want += &multipole_to_local(&b, t);
        let got = multipole_to_local(&sum, t);
        for (k, (g, w)) in got.coeffs().iter().zip(want.coeffs()).enumerate() {
            assert_close(
                *g,
                *w,
                1e-9 * (1.0 + w.re.abs() + w.im.abs()),
                &format!("slot {k}"),
            );
        }
    }

    /// The physics check. Bodies in a cell at the origin, a target center `t`
    /// away, evaluation points inside the target cell: the converted expansion
    /// must give the direct potential, with an error that shrinks as p grows.
    #[test]
    fn the_converted_expansion_gives_the_direct_potential_inside_the_target() {
        let bodies = bodies_around([0.0; 3], 0.5, 6, 3);
        // order -> allowed worst relative error, about five times the measured one
        let allowed = [(2, 1e-2), (4, 1e-4), (6, 3e-6), (8, 1e-7)];
        let points = [
            [0.2, -0.1, 0.3],
            [-0.4, 0.3, 0.1],
            [0.0, 0.0, 0.35],
            [0.1, 0.2, -0.3],
        ];
        for (name, t) in [
            ("generic", [3.1, -2.2, 2.6]),
            ("+z", [0.0, 0.0, 4.0]),
            ("-z", [0.0, 0.0, -4.0]),
            ("+x", [4.0, 0.0, 0.0]),
            ("diagonal", [-3.0, 3.0, 0.0]),
        ] {
            let mut worst_at = vec![];
            for p in [0usize, 2, 4, 6, 8] {
                let mp = particle_to_multipole(bodies.iter(), [0.0; 3], p);
                let local = multipole_to_local(&mp, t);
                let mut worst: f64 = 0.0;
                for y in points {
                    let x = [t[0] + y[0], t[1] + y[1], t[2] + y[2]];
                    let want = direct_potential(&bodies, x);
                    let got = local.value(y);
                    assert!(
                        got.im.abs() < 1e-12 * want,
                        "{name} p={p}: imaginary part {}",
                        got.im
                    );
                    worst = worst.max((got.re - want).abs() / want);
                }
                worst_at.push((p, worst));
            }
            for ((p, bound), (q, got)) in allowed.iter().zip(worst_at.iter().skip(1)) {
                assert_eq!(p, q);
                assert!(got < bound, "{name} p={p}: error {got} above {bound}");
            }
            assert!(
                worst_at[0].1 > 1e-2,
                "{name}: p = 0 should be a visibly rough approximation"
            );
            assert!(
                worst_at[4].1 < worst_at[1].1 * 1e-2,
                "{name}: error is not shrinking with p"
            );
        }
    }

    #[test]
    #[should_panic(expected = "requires irregular")]
    fn m2l_refuses_a_regular_input() {
        multipole_to_local(&SolidHarmonic::regular(2), [3.0, 0.0, 0.0]);
    }

    #[test]
    #[should_panic(expected = "zero shift")]
    fn m2l_refuses_a_zero_shift() {
        multipole_to_local(&SolidHarmonic::irregular(2), [0.0; 3]);
    }
}
