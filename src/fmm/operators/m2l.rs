use crate::{
    fmm::{
        complex::Complex,
        harmonic::{SolidHarmonic, from_index, harmonics},
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
    todo!()
}

fn get_local_coefficient(mp: &SolidHarmonic, table: &[Complex], j: i64, k: i64) -> Complex {
    todo!()
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
    use crate::fmm::{harmonic::to_index, test_util::assert_close};

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
}
