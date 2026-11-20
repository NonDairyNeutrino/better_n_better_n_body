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

fn to_index(n: i64, m: i64) -> usize {
    return (n.pow(2) + n + m) as usize;
}

fn from_index(idx: usize) -> (i64, i64) {
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
fn harmonics(pos: Vec3, p: usize) -> Vec<Complex> {
    let [x, y, z] = pos;
    let r: f64 = (x * x + y * y + z * z).sqrt();
    let rho: f64 = (x * x + y * y).sqrt();
    let cos_t: f64 = z / r;
    let sin_t: f64 = rho / r;
    // e^{i phi}; on the z-axis phi is undefined and sin_t = 0 kills m != 0.
    let e_iphi: Complex = if rho == 0.0 {
        Complex { re: 1.0, im: 0.0 }
    } else {
        Complex { re: x / rho, im: y / rho }
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
                ((2.0 * nf - 1.0) * cos_t * prev
                    - ((nf - 1.0).powi(2) - mf * mf).sqrt() * prev2)
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

    #[test]
    fn every_pair_n_m_has_its_own_slot_and_every_slot_is_used() {
        let p: i64 = 6;
        let sh = numbered(p as usize);
        let mut seen: Vec<usize> = Vec::new();
        for n in 0..=p {
            for m in -n..=n {
                seen.push(sh.get(n, m).re as usize);
            }
        }
        seen.sort();
        let want: Vec<usize> = (0..pair_count(p as usize)).collect();
        assert_eq!(seen, want);
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

    #[test]
    fn writes_through_coeffs_mut_are_seen_by_coeffs() {
        let mut sh = SolidHarmonic::irregular(2);
        sh.coeffs_mut()[5] = Complex { re: 1.5, im: -2.5 };
        assert_eq!(sh.coeffs()[5].re, 1.5);
        assert_eq!(sh.coeffs()[5].im, -2.5);
        assert!(
            matches!(sh, SolidHarmonic::Irregular(_)),
            "kind must not change"
        );
    }
}
