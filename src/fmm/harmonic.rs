use crate::fmm::complex::Complex;

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
    fn coeffs(&self) -> &[Complex] {
        match self {
            SolidHarmonic::Regular(v) | SolidHarmonic::Irregular(v) => v,
        }
    }
    fn coeffs_mut(&mut self) -> &mut [Complex] {
        match self {
            SolidHarmonic::Regular(v) | SolidHarmonic::Irregular(v) => v,
        }
    }
    fn get(&self, n: i64, m: i64) -> Complex {
        assert!(0 <= n, "n must be non-negative, got n={n}");
        assert!(-n <= m && m <= n, "need |m| <= n, got n={n}, m={m}");
        let idx: usize = (n.pow(2) + n + m) as usize;
        assert!(
            idx < self.coeffs().len(),
            "(n={n}, m={m}) is beyond the order of this solid harmonic"
        );
        let coeff: Complex = self.coeffs()[idx];
        return coeff;
    }
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
