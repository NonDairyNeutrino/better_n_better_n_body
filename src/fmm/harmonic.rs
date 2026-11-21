use std::ops::AddAssign;

use crate::{fmm::complex::Complex, shared::Vec3};

#[derive(Debug, Clone, PartialEq)]
pub enum SolidHarmonic {
    Regular(Vec<Complex>),
    Irregular(Vec<Complex>),
}
impl SolidHarmonic {
    pub fn regular(p: usize) -> SolidHarmonic {
        let coeffs: Vec<Complex> = vec![Complex::zero(); (p + 1).pow(2)];
        return SolidHarmonic::Regular(coeffs);
    }
    pub fn irregular(p: usize) -> SolidHarmonic {
        let coeffs: Vec<Complex> = vec![Complex::zero(); (p + 1).pow(2)];
        return SolidHarmonic::Irregular(coeffs);
    }
    pub fn coeffs(&self) -> &[Complex] {
        match self {
            SolidHarmonic::Regular(v) | SolidHarmonic::Irregular(v) => v,
        }
    }
    pub fn coeffs_mut(&mut self) -> &mut [Complex] {
        match self {
            SolidHarmonic::Regular(v) | SolidHarmonic::Irregular(v) => v,
        }
    }
    pub fn get(&self, n: i64, m: i64) -> Complex {
        assert!(0 <= n, "n must be non-negative, got n={n}");
        assert!(-n <= m && m <= n, "need |m| <= n, got n={n}, m={m}");
        let idx: usize = to_index(n, m);
        assert!(
            idx < self.coeffs().len(),
            "(n={n}, m={m}) is beyond the order of this solid harmonic"
        );
        let coeff: Complex = self.coeffs()[idx];
        return coeff;
    }
    pub fn order(&self) -> usize {
        return self.coeffs().len().isqrt() - 1;
    }
    /// Evaluates the series at `pos`, measured from the expansion center.
    ///
    /// Regular:   sum c_nm * r^n     * Y_nm(direction)
    /// Irregular: sum c_nm * r^-(n+1) * Y_nm(direction)
    pub fn value(&self, pos: Vec3) -> Complex {
        let coeffs: &[Complex] = self.coeffs();
        let r: f64 = pos.iter().map(|d| d * d).sum::<f64>().sqrt();

        if r == 0.0 {
            return match self {
                // Only the n = 0 term survives; its Y_00 = 1.
                SolidHarmonic::Regular(v) => v[0],
                SolidHarmonic::Irregular(_) => {
                    panic!("an irregular solid harmonic diverges at its center")
                }
            };
        }

        let bases: Vec<Complex> = harmonics(pos, self.order());
        let mut val: Complex = Complex::zero();
        for (idx, (c, y)) in coeffs.iter().zip(bases).enumerate() {
            let (n, _) = from_index(idx);
            let scale: f64 = match self {
                SolidHarmonic::Regular(_) => r.powi(n as i32),
                SolidHarmonic::Irregular(_) => r.powi(-(n as i32) - 1),
            };
            val += *c * y * scale;
        }
        return val;
    }
}
impl AddAssign<&SolidHarmonic> for SolidHarmonic {
    fn add_assign(&mut self, rhs: &SolidHarmonic) {
        assert!(
            std::mem::discriminant(self) == std::mem::discriminant(rhs),
            "cannot add a regular and an irregular solid harmonic"
        );
        assert_eq!(
            self.order(),
            rhs.order(),
            "cannot add solid harmonics of different order"
        );
        for (a, b) in self.coeffs_mut().iter_mut().zip(rhs.coeffs()) {
            *a += *b;
        }
    }
}

pub fn to_index(n: i64, m: i64) -> usize {
    return (n.pow(2) + n + m) as usize;
}

pub fn from_index(idx: usize) -> (i64, i64) {
    let n: i64 = idx.isqrt() as i64;
    let m: i64 = (idx as i64) - n.pow(2) - n;
    return (n, m);
}

/// Y_n^m(theta, phi) for every 0 <= n <= p, |m| <= n, in `to_index` order.
/// Depends only on the direction of `pos`, not its length; requires pos != 0.
///
/// Y_n^m = Pbar_n^m(cos theta) e^{i m phi},  Y_n^-m = conj(Y_n^m), where
/// Pbar_n^m = sqrt((n-m)!/(n+m)!) P_n^m uses the Condon-Shortley sign
/// (-1)^m inside P_n^m. No 4pi/(2n+1) factor.
pub fn harmonics(pos: Vec3, p: usize) -> Vec<Complex> {
    let [x, y, z] = pos;
    let r: f64 = (x * x + y * y + z * z).sqrt();
    let rho: f64 = (x * x + y * y).sqrt();
    let cos_t: f64 = z / r;
    let sin_t: f64 = rho / r;
    // e^{i phi}; on the z-axis phi is undefined and sin_t = 0 kills m != 0.
    let e_iphi: Complex = if rho == 0.0 {
        Complex { re: 1.0, im: 0.0 }
    } else {
        Complex {
            re: x / rho,
            im: y / rho,
        }
    };

    let mut out: Vec<Complex> = vec![Complex::zero(); (p + 1).pow(2)];
    let mut pmm: f64 = 1.0; // Pbar_m^m
    let mut phase: Complex = Complex { re: 1.0, im: 0.0 }; // e^{i m phi}

    for m in 0..=p {
        let mf: f64 = m as f64;
        let mut prev2: f64 = 0.0; // Pbar_{n-2}^m
        let mut prev: f64 = 0.0; // Pbar_{n-1}^m
        for n in m..=p {
            let nf: f64 = n as f64;
            let value: f64 = if n == m {
                pmm
            } else if n == m + 1 {
                (2.0 * mf + 1.0).sqrt() * cos_t * pmm
            } else {
                ((2.0 * nf - 1.0) * cos_t * prev - ((nf - 1.0).powi(2) - mf * mf).sqrt() * prev2)
                    / (nf * nf - mf * mf).sqrt()
            };
            prev2 = prev;
            prev = value;

            let y_nm: Complex = phase * value;
            out[to_index(n as i64, m as i64)] = y_nm;
            out[to_index(n as i64, -(m as i64))] = y_nm.conj();
        }
        pmm *= -((2.0 * mf + 1.0) / (2.0 * mf + 2.0)).sqrt() * sin_t;
        phase = phase * e_iphi;
    }
    return out;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fmm::test_util::{assert_close, dist3};

    // ---------------------------------------------------------------- helpers

    /// The number of (n, m) pairs with 0 <= n <= p and |m| <= n, counted directly
    /// rather than with the closed formula the code uses.
    fn pair_count(p: usize) -> usize {
        (0..=p).map(|n| 2 * n + 1).sum()
    }

    // ---------------------------------------------------------- SolidHarmonic

    #[test]
    fn a_solid_harmonic_has_one_coefficient_per_pair_n_m() {
        for p in 0..=8 {
            assert_eq!(
                SolidHarmonic::regular(p).coeffs().len(),
                pair_count(p),
                "regular p={p}"
            );
            assert_eq!(
                SolidHarmonic::irregular(p).coeffs().len(),
                pair_count(p),
                "irregular p={p}"
            );
        }
    }

    #[test]
    fn the_constructors_make_the_variant_they_are_named_for() {
        assert!(matches!(
            SolidHarmonic::regular(2),
            SolidHarmonic::Regular(_)
        ));
        assert!(matches!(
            SolidHarmonic::irregular(2),
            SolidHarmonic::Irregular(_)
        ));
    }

    #[test]
    fn a_new_solid_harmonic_is_all_zeros() {
        for sh in [SolidHarmonic::regular(4), SolidHarmonic::irregular(4)] {
            assert!(!sh.coeffs().is_empty());
            assert!(sh.coeffs().iter().all(|c| c.re == 0.0 && c.im == 0.0));
        }
    }

    /// Fill slot k with the number k, so `get(n, m)` reports which slot it read.
    fn numbered(p: usize) -> SolidHarmonic {
        let mut sh = SolidHarmonic::regular(p);
        for (k, c) in sh.coeffs_mut().iter_mut().enumerate() {
            *c = Complex {
                re: k as f64,
                im: -(k as f64),
            };
        }
        sh
    }

    #[test]
    fn get_reads_the_slots_in_the_documented_order() {
        // (0,0) | (1,-1) (1,0) (1,1) | (2,-2) (2,-1) (2,0) (2,1) (2,2) | ...
        let sh = numbered(3);
        let expected = [
            (0, 0, 0.0),
            (1, -1, 1.0),
            (1, 0, 2.0),
            (1, 1, 3.0),
            (2, -2, 4.0),
            (2, 0, 6.0),
            (2, 2, 8.0),
            (3, -3, 9.0),
            (3, 3, 15.0),
        ];
        for (n, m, slot) in expected {
            assert_eq!(sh.get(n, m).re, slot, "(n={n}, m={m})");
            assert_eq!(sh.get(n, m).im, -slot, "(n={n}, m={m}) imaginary part");
        }
    }

    /// The formula keeps working far past the orders anyone will use.
    #[test]
    fn the_first_and_last_pairs_of_a_high_order_are_the_first_and_last_slots() {
        let p: i64 = 60;
        let sh = numbered(p as usize);
        assert_eq!(sh.get(0, 0).re as usize, 0);
        assert_eq!(sh.get(p, -p).re as usize, (p * p) as usize);
        assert_eq!(sh.get(p, p).re as usize, pair_count(p as usize) - 1);
    }

    #[test]
    #[should_panic(expected = "|m| <= n")]
    fn get_rejects_m_larger_than_n() {
        numbered(3).get(1, 2);
    }

    #[test]
    #[should_panic(expected = "|m| <= n")]
    fn get_rejects_m_smaller_than_minus_n() {
        numbered(3).get(1, -2);
    }

    #[test]
    #[should_panic(expected = "n must be non-negative")]
    fn get_rejects_a_negative_n() {
        numbered(3).get(-1, 0);
    }

    #[test]
    #[should_panic(expected = "beyond the order")]
    fn get_rejects_an_n_beyond_the_order() {
        numbered(3).get(4, 0);
    }

    // ------------------------------------------------- harmonics and value

    /// P_n(x) by Bonnet's recurrence: a different method from the code under test.
    fn legendre(n: usize, x: f64) -> f64 {
        let (mut a, mut b) = (1.0, x);
        if n == 0 {
            return a;
        }
        for k in 1..n {
            let kf = k as f64;
            (a, b) = (b, ((2.0 * kf + 1.0) * x * b - kf * a) / (kf + 1.0));
        }
        b
    }

    /// A solid harmonic of order p with a single coefficient (n, m) equal to 1.
    fn single(regular: bool, p: usize, n: i64, m: i64) -> SolidHarmonic {
        let mut sh = if regular {
            SolidHarmonic::regular(p)
        } else {
            SolidHarmonic::irregular(p)
        };
        sh.coeffs_mut()[to_index(n, m)] = Complex { re: 1.0, im: 0.0 };
        sh
    }

    /// Second-difference Laplacian of a complex-valued function of position.
    fn laplacian(f: impl Fn(Vec3) -> Complex, at: Vec3) -> Complex {
        let h = 1e-3;
        let mut sum = Complex { re: 0.0, im: 0.0 };
        for axis in 0..3 {
            let (mut up, mut down) = (at, at);
            up[axis] += h;
            down[axis] -= h;
            sum += f(up) + f(down) - f(at) * 2.0;
        }
        sum * (1.0 / (h * h))
    }

    // ---- index helpers

    #[test]
    fn to_index_and_from_index_are_inverses() {
        for idx in 0..400 {
            let (n, m) = from_index(idx);
            assert!(m.abs() <= n, "slot {idx} gave (n={n}, m={m})");
            assert_eq!(to_index(n, m), idx);
        }
        for n in 0..20 {
            for m in -n..=n {
                assert_eq!(from_index(to_index(n, m)), (n, m));
            }
        }
    }

    #[test]
    fn order_is_the_p_the_solid_harmonic_was_built_with() {
        for p in 0..10 {
            assert_eq!(SolidHarmonic::regular(p).order(), p);
            assert_eq!(SolidHarmonic::irregular(p).order(), p);
        }
    }

    // ---- harmonics

    /// Textbook closed forms at (1, 2, 2): r=3, cos t=2/3, sin t=sqrt5/3,
    /// e^{i phi}=(1+2i)/sqrt5. These are not computed with the recurrence.
    #[test]
    fn harmonics_match_hand_computed_values_up_to_degree_two() {
        let y = harmonics([1.0, 2.0, 2.0], 2);
        let y11 = Complex {
            re: -1.0 / (3.0 * 2f64.sqrt()),
            im: -2.0 / (3.0 * 2f64.sqrt()),
        };
        let k21 = -(1.5f64).sqrt() * 2.0 / 9.0;
        let y21 = Complex {
            re: k21,
            im: 2.0 * k21,
        };
        let k22 = (3.0f64 / 8.0).sqrt() / 9.0;
        let y22 = Complex {
            re: -3.0 * k22,
            im: 4.0 * k22,
        };
        let want = [
            (0, 0, Complex { re: 1.0, im: 0.0 }),
            (1, -1, y11.conj()),
            (
                1,
                0,
                Complex {
                    re: 2.0 / 3.0,
                    im: 0.0,
                },
            ),
            (1, 1, y11),
            (2, -2, y22.conj()),
            (2, -1, y21.conj()),
            (
                2,
                0,
                Complex {
                    re: 1.0 / 6.0,
                    im: 0.0,
                },
            ),
            (2, 1, y21),
            (2, 2, y22),
        ];
        for (n, m, w) in want {
            assert_close(y[to_index(n, m)], w, 1e-14, &format!("Y(n={n}, m={m})"));
        }
    }

    #[test]
    fn harmonics_do_not_depend_on_the_length_of_pos() {
        let dir = [0.4, -0.7, 0.2];
        let base = harmonics(dir, 6);
        for k in [0.01, 3.0, 1e3] {
            let scaled = harmonics([dir[0] * k, dir[1] * k, dir[2] * k], 6);
            for (idx, (a, b)) in base.iter().zip(&scaled).enumerate() {
                assert_close(*b, *a, 1e-12, &format!("slot {idx}, scale {k}"));
            }
        }
    }

    /// sum_m Y_n^m(a) conj(Y_n^m(b)) = P_n(cos gamma), with gamma the angle
    /// between a and b. Checks every degree up to 12 against an independent P_n.
    #[test]
    fn addition_theorem_holds_through_degree_twelve() {
        let a = [1.0, 2.0, 2.0];
        let b = [-0.5, 0.3, 1.7];
        let dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let cos_gamma = dot / (dist3(a, [0.0; 3]) * dist3(b, [0.0; 3]));
        assert!(
            cos_gamma.abs() < 0.95,
            "angle too degenerate to be a good test"
        );
        let (ya, yb) = (harmonics(a, 12), harmonics(b, 12));
        for n in 0..=12i64 {
            let mut sum = Complex { re: 0.0, im: 0.0 };
            for m in -n..=n {
                sum += ya[to_index(n, m)] * yb[to_index(n, m)].conj();
            }
            assert_close(
                sum,
                Complex {
                    re: legendre(n as usize, cos_gamma),
                    im: 0.0,
                },
                1e-12,
                &format!("n={n}"),
            );
        }
    }

    #[test]
    fn on_the_z_axis_only_m_zero_survives_and_it_is_plus_or_minus_one_to_the_n() {
        for z in [2.0, -3.0] {
            let y = harmonics([0.0, 0.0, z], 8);
            for (idx, v) in y.iter().enumerate() {
                let (n, m) = from_index(idx);
                assert!(
                    v.re.is_finite() && v.im.is_finite(),
                    "(n={n}, m={m}) not finite"
                );
                let want = if m != 0 {
                    0.0
                } else if z < 0.0 && n % 2 == 1 {
                    -1.0
                } else {
                    1.0
                };
                assert_close(
                    *v,
                    Complex { re: want, im: 0.0 },
                    1e-12,
                    &format!("z={z} (n={n}, m={m})"),
                );
            }
        }
    }

    // ---- value

    #[test]
    fn value_of_hand_built_series_at_1_2_2() {
        let at = [1.0, 2.0, 2.0]; // r = 3
        let mut m0 = SolidHarmonic::irregular(2);
        m0.coeffs_mut()[to_index(0, 0)] = Complex { re: 2.0, im: 0.0 };
        assert_close(
            m0.value(at),
            Complex {
                re: 2.0 / 3.0,
                im: 0.0,
            },
            1e-14,
            "irregular c00=2 is 2/r",
        );
        assert_close(
            single(false, 2, 1, 0).value(at),
            Complex {
                re: 2.0 / 27.0,
                im: 0.0,
            },
            1e-14,
            "irregular c10 is z/r^3",
        );
        assert_close(
            single(false, 2, 2, 0).value(at),
            Complex {
                re: 1.0 / 162.0,
                im: 0.0,
            },
            1e-14,
            "irregular c20 is (3z^2-r^2)/(2 r^5)",
        );
        assert_close(
            single(true, 2, 0, 0).value(at),
            Complex { re: 1.0, im: 0.0 },
            1e-14,
            "regular c00 is 1",
        );
        assert_close(
            single(true, 2, 1, 0).value(at),
            Complex { re: 2.0, im: 0.0 },
            1e-14,
            "regular c10 is z",
        );
        assert_close(
            single(true, 2, 2, 0).value(at),
            Complex { re: 1.5, im: 0.0 },
            1e-14,
            "regular c20 is (3z^2-r^2)/2",
        );
        let s = 1.0 / 2f64.sqrt();
        assert_close(
            single(true, 2, 1, 1).value(at),
            Complex {
                re: -s,
                im: -2.0 * s,
            },
            1e-14,
            "regular c11 is -(x+iy)/sqrt2",
        );
        assert_close(
            single(true, 2, 1, -1).value(at),
            Complex {
                re: -s,
                im: 2.0 * s,
            },
            1e-14,
            "regular c1,-1 is -(x-iy)/sqrt2",
        );
    }

    #[test]
    fn regular_terms_grow_as_r_to_the_n_and_irregular_fall_as_r_to_minus_n_plus_1() {
        let dir = [0.3, 0.5, -0.4];
        for n in 0..=4i64 {
            for m in -n..=n {
                let (reg, irr) = (single(true, 4, n, m), single(false, 4, n, m));
                let (v1r, v1i) = (reg.value(dir), irr.value(dir));
                for v in [v1r, v1i] {
                    assert!(
                        v.re.abs() + v.im.abs() > 1e-6,
                        "(n={n}, m={m}) term is zero; test would be vacuous"
                    );
                }
                for k in [0.5, 2.0] {
                    let far = [dir[0] * k, dir[1] * k, dir[2] * k];
                    assert_close(
                        reg.value(far),
                        v1r * k.powi(n as i32),
                        1e-12,
                        &format!("regular (n={n}, m={m}) x{k}"),
                    );
                    assert_close(
                        irr.value(far),
                        v1i * k.powi(-(n as i32) - 1),
                        1e-12,
                        &format!("irregular (n={n}, m={m}) x{k}"),
                    );
                }
            }
        }
    }

    /// Both kinds solve Laplace's equation away from the origin; a wrong power
    /// of r in the scale breaks this for every n except by accident.
    #[test]
    fn both_kinds_satisfy_laplaces_equation() {
        let at = [0.7, -0.4, 0.9];
        for regular in [true, false] {
            for n in 0..=4i64 {
                for m in -n..=n {
                    let sh = single(regular, 4, n, m);
                    let v = sh.value(at);
                    assert!(
                        v.re.abs() + v.im.abs() > 1e-6,
                        "(n={n}, m={m}) term is zero; test would be vacuous"
                    );
                    let lap = laplacian(|q| sh.value(q), at);
                    assert_close(
                        lap,
                        Complex { re: 0.0, im: 0.0 },
                        1e-3,
                        &format!("regular={regular} (n={n}, m={m})"),
                    );
                }
            }
        }
        // Non-vacuity: the same check must fail for a function that is not harmonic.
        let not_harmonic = laplacian(
            |q| Complex {
                re: q[0] * q[0] + q[1] * q[1] + q[2] * q[2],
                im: 0.0,
            },
            at,
        );
        assert!(not_harmonic.re > 1.0);
    }

    #[test]
    fn regular_value_at_the_center_is_the_n_zero_coefficient() {
        let mut sh = SolidHarmonic::regular(3);
        for (k, c) in sh.coeffs_mut().iter_mut().enumerate() {
            *c = Complex {
                re: k as f64 + 3.0,
                im: -1.0,
            };
        }
        assert_close(
            sh.value([0.0; 3]),
            Complex { re: 3.0, im: -1.0 },
            0.0,
            "center",
        );
    }

    #[test]
    #[should_panic(expected = "diverges")]
    fn irregular_value_at_the_center_panics() {
        SolidHarmonic::irregular(2).value([0.0; 3]);
    }

    /// Rebuilds 1/|x - src| as a series about the origin and checks every order
    /// against the geometric tail bound rho^(p+1) / (far (1 - rho)), rho = near/far,
    /// where near and far are the smaller and larger of |src| and |at|.
    /// `regular` selects the series valid near the origin (|at| < |src|); otherwise
    /// the one valid far from it (|src| < |at|).
    fn check_point_source_series(regular: bool, src: Vec3, at: Vec3) {
        let (rs, r) = (dist3(src, [0.0; 3]), dist3(at, [0.0; 3]));
        let (near, far) = if regular { (r, rs) } else { (rs, r) };
        assert!(near < far, "the series only converges for near < far");
        let rho = near / far;
        let want = 1.0 / dist3(at, src);
        let mut errs = Vec::new();
        for p in [0usize, 1, 2, 4, 8, 12] {
            let y = harmonics(src, p);
            let mut sh = if regular {
                SolidHarmonic::regular(p)
            } else {
                SolidHarmonic::irregular(p)
            };
            for (idx, c) in sh.coeffs_mut().iter_mut().enumerate() {
                let (n, _) = from_index(idx);
                let radial = if regular {
                    rs.powi(-(n as i32) - 1)
                } else {
                    rs.powi(n as i32)
                };
                *c = y[idx].conj() * radial;
            }
            let v = sh.value(at);
            let err = (v.re - want).abs();
            let bound = rho.powi(p as i32 + 1) / (far * (1.0 - rho));
            assert!(
                v.im.abs() < 1e-12,
                "p={p}: potential must be real, got im={}",
                v.im
            );
            assert!(
                err <= bound * (1.0 + 1e-9) + 1e-14,
                "p={p}: err {err} exceeds tail bound {bound}"
            );
            errs.push(err);
        }
        assert!(errs[0] > 1e-2, "p=0 must be visibly wrong, err={}", errs[0]);
        assert!(
            *errs.last().unwrap() < 1e-4,
            "p=12 err={}",
            errs.last().unwrap()
        );
    }

    /// A source at x' seen from a farther point x:
    /// 1/|x-x'| = sum_nm conj(R_nm(x')) S_nm(x).
    #[test]
    fn irregular_series_of_a_point_source_converges_to_one_over_distance() {
        check_point_source_series(false, [0.3, -0.2, 0.4], [0.9, 0.6, 0.7]);
    }

    /// The mirror image, a far source seen from a point near the center:
    /// 1/|x-x'| = sum_nm conj(S_nm(x')) R_nm(x).
    #[test]
    fn regular_series_of_a_far_source_converges_to_one_over_distance() {
        check_point_source_series(true, [0.9, 0.6, 0.7], [0.3, -0.2, 0.4]);
    }

    // ---- AddAssign

    /// A solid harmonic of the given kind and order whose slot k holds `f(k)`.
    fn filled(regular: bool, p: usize, f: impl Fn(usize) -> Complex) -> SolidHarmonic {
        let mut sh = if regular {
            SolidHarmonic::regular(p)
        } else {
            SolidHarmonic::irregular(p)
        };
        for (k, c) in sh.coeffs_mut().iter_mut().enumerate() {
            *c = f(k);
        }
        sh
    }

    #[test]
    fn adding_solid_harmonics_adds_slot_by_slot_and_keeps_the_kind() {
        for regular in [true, false] {
            let mut a = filled(regular, 3, |k| Complex {
                re: k as f64,
                im: -(k as f64),
            });
            let b = filled(regular, 3, |k| Complex {
                re: 2.0 * k as f64 + 1.0,
                im: 3.0,
            });
            a += &b;
            for k in 0..16 {
                assert_close(
                    a.coeffs()[k],
                    Complex {
                        re: 3.0 * k as f64 + 1.0,
                        im: 3.0 - k as f64,
                    },
                    0.0,
                    &format!("regular={regular}, slot {k}"),
                );
            }
            assert_eq!(matches!(a, SolidHarmonic::Regular(_)), regular);
        }
    }

    #[test]
    fn adding_an_empty_solid_harmonic_changes_nothing() {
        let mut a = filled(false, 2, |k| Complex {
            re: k as f64 + 1.0,
            im: 0.5,
        });
        let before = a.clone();
        assert_ne!(before, SolidHarmonic::irregular(2), "need a nonzero start");
        a += &SolidHarmonic::irregular(2);
        assert_eq!(a, before);
    }

    #[test]
    #[should_panic(expected = "regular and an irregular")]
    fn adding_an_irregular_to_a_regular_panics() {
        let mut a = SolidHarmonic::regular(2);
        a += &SolidHarmonic::irregular(2);
    }

    #[test]
    #[should_panic(expected = "regular and an irregular")]
    fn adding_a_regular_to_an_irregular_panics() {
        let mut a = SolidHarmonic::irregular(2);
        a += &SolidHarmonic::regular(2);
    }

    #[test]
    #[should_panic(expected = "different order")]
    fn adding_solid_harmonics_of_different_order_panics() {
        let mut a = SolidHarmonic::irregular(2);
        a += &SolidHarmonic::irregular(3);
    }

    #[test]
    #[should_panic(expected = "different order")]
    fn adding_a_lower_order_to_a_higher_order_panics_too() {
        let mut a = SolidHarmonic::regular(3);
        a += &SolidHarmonic::regular(2);
    }
}
