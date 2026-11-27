//! A minimal complex number type, with just the arithmetic the FMM needs.

use std::{
    iter::Sum,
    ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

/// A complex number `re + i * im`.
///
/// Supports `+`, `-`, `*` (by another `Complex` or by an `f64`), unary `-`, their
/// assigning forms, and `sum()` over iterators. It is `Copy`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}
impl Complex {
    /// The number 0.
    pub fn zero() -> Complex {
        Complex {
            re: 0.0_f64,
            im: 0.0_f64,
        }
    }
    /// The complex conjugate: the same real part and the opposite imaginary part.
    pub fn conj(self) -> Complex {
        Complex {
            re: self.re,
            im: -self.im,
        }
    }
}

impl Add for Complex {
    type Output = Complex;
    fn add(self, rhs: Self) -> Self::Output {
        Complex {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl Sum for Complex {
    fn sum<I: Iterator<Item = Complex>>(iter: I) -> Complex {
        iter.fold(Complex::zero(), |acc, x| acc + x)
    }
}

impl<'a> Sum<&'a Complex> for Complex {
    fn sum<I: Iterator<Item = &'a Complex>>(iter: I) -> Complex {
        iter.fold(Complex::zero(), |acc, x| acc + *x)
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, rhs: Self) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}

impl Sub for Complex {
    type Output = Complex;
    fn sub(self, rhs: Self) -> Self::Output {
        Complex {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}

impl SubAssign for Complex {
    fn sub_assign(&mut self, rhs: Self) {
        self.re -= rhs.re;
        self.im -= rhs.im;
    }
}

impl Mul<Complex> for Complex {
    type Output = Complex;
    fn mul(self, rhs: Self) -> Self::Output {
        let a = self.re;
        let b = self.im;
        let c = rhs.re;
        let d = rhs.im;
        Complex {
            re: a * c - b * d,
            im: a * d + b * c,
        }
    }
}

impl Mul<f64> for Complex {
    type Output = Complex;
    fn mul(self, rhs: f64) -> Self::Output {
        Complex {
            re: rhs * self.re,
            im: rhs * self.im,
        }
    }
}

impl MulAssign for Complex {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Neg for Complex {
    type Output = Complex;
    fn neg(self) -> Self::Output {
        Complex {
            re: -self.re,
            im: -self.im,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---------------------------------------------------------------- helpers

    fn c(re: f64, im: f64) -> Complex {
        Complex { re, im }
    }

    /// Equal up to rounding, relative to the size of the numbers involved.
    fn close(a: Complex, b: Complex) -> bool {
        let scale = 1.0 + a.re.abs().max(a.im.abs()).max(b.re.abs()).max(b.im.abs());
        (a.re - b.re).abs() <= 1e-12 * scale && (a.im - b.im).abs() <= 1e-12 * scale
    }

    /// Values with mixed signs, a zero part, and non-trivial parts, so that a
    /// formula that swaps or drops a part cannot hide.
    fn samples() -> Vec<Complex> {
        vec![
            c(1.5, -2.0),
            c(0.0, 3.25),
            c(-4.0, 0.5),
            c(2.0, 2.0),
            c(-0.75, -1.25),
            c(3.0, 0.0),
        ]
    }

    // ------------------------------------------------------- known answers

    #[test]
    fn zero_is_zero() {
        assert_eq!(Complex::zero(), c(0.0, 0.0));
    }

    #[test]
    fn addition_adds_each_part() {
        assert_eq!(c(1.0, 2.0) + c(3.0, -4.0), c(4.0, -2.0));
    }

    #[test]
    fn subtraction_subtracts_each_part() {
        assert_eq!(c(1.0, 2.0) - c(3.0, -4.0), c(-2.0, 6.0));
    }

    #[test]
    fn multiplication_follows_the_complex_rule() {
        // (1 + 2i)(3 - i) = 3 - i + 6i - 2i^2 = 5 + 5i
        assert_eq!(c(1.0, 2.0) * c(3.0, -1.0), c(5.0, 5.0));
    }

    #[test]
    fn i_squared_is_minus_one() {
        assert_eq!(c(0.0, 1.0) * c(0.0, 1.0), c(-1.0, 0.0));
    }

    #[test]
    fn i_to_the_fourth_is_one() {
        let i = c(0.0, 1.0);
        assert_eq!(i * i * i * i, c(1.0, 0.0));
    }

    #[test]
    fn scaling_by_a_real_scales_each_part() {
        assert_eq!(c(1.0, 2.0) * 3.0, c(3.0, 6.0));
        assert_eq!(c(1.0, -2.0) * -0.5, c(-0.5, 1.0));
        assert_eq!(c(1.0, 2.0) * 0.0, c(0.0, 0.0));
    }

    #[test]
    fn negation_flips_both_parts() {
        assert_eq!(-c(1.0, -2.0), c(-1.0, 2.0));
    }

    #[test]
    fn the_conjugate_flips_only_the_imaginary_part() {
        assert_eq!(c(1.0, 2.0).conj(), c(1.0, -2.0));
        assert_eq!(c(-3.0, -0.5).conj(), c(-3.0, 0.5));
    }

    // ------------------------------------------------------ algebraic laws

    #[test]
    fn addition_and_multiplication_commute() {
        for a in samples() {
            for b in samples() {
                assert!(close(a + b, b + a), "{a:?} + {b:?}");
                assert!(close(a * b, b * a), "{a:?} * {b:?}");
            }
        }
    }

    #[test]
    fn multiplication_is_associative_and_distributes_over_addition() {
        for a in samples() {
            for b in samples() {
                for d in samples() {
                    assert!(close((a * b) * d, a * (b * d)), "assoc {a:?} {b:?} {d:?}");
                    assert!(
                        close(a * (b + d), a * b + a * d),
                        "distrib {a:?} {b:?} {d:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn zero_and_one_are_the_identities() {
        let one = c(1.0, 0.0);
        for a in samples() {
            assert_eq!(a + Complex::zero(), a);
            assert!(close(a * one, a));
            assert_eq!(a - a, Complex::zero());
            assert_eq!(a + (-a), Complex::zero());
            assert_eq!(a * Complex::zero(), Complex::zero());
        }
    }

    #[test]
    fn subtraction_is_adding_the_negative() {
        for a in samples() {
            for b in samples() {
                assert!(close(a - b, a + (-b)));
            }
        }
    }

    #[test]
    fn scaling_by_a_real_agrees_with_multiplying_by_a_real_complex() {
        for a in samples() {
            for s in [-2.5, 0.0, 1.0, 4.0] {
                assert!(close(a * s, a * c(s, 0.0)), "{a:?} * {s}");
            }
        }
    }

    // --------------------------------------------------------- conjugation

    #[test]
    fn conjugating_twice_is_the_identity() {
        for a in samples() {
            assert_eq!(a.conj().conj(), a);
        }
    }

    #[test]
    fn a_number_times_its_conjugate_is_its_squared_size() {
        for a in samples() {
            let want = c(a.re * a.re + a.im * a.im, 0.0);
            assert!(close(a * a.conj(), want), "{a:?}");
        }
    }

    #[test]
    fn the_conjugate_of_a_product_is_the_product_of_the_conjugates() {
        for a in samples() {
            for b in samples() {
                assert!(close((a * b).conj(), a.conj() * b.conj()));
                assert!(close((a + b).conj(), a.conj() + b.conj()));
            }
        }
    }

    // ------------------------------------------- independent of the formula

    /// |a b|^2 = |a|^2 |b|^2. Neither side uses the multiplication rule's shape.
    #[test]
    fn multiplying_multiplies_the_sizes() {
        let size2 = |z: Complex| z.re * z.re + z.im * z.im;
        for a in samples() {
            for b in samples() {
                let (lhs, rhs) = (size2(a * b), size2(a) * size2(b));
                assert!((lhs - rhs).abs() <= 1e-12 * (1.0 + rhs), "{a:?} {b:?}");
            }
        }
    }

    /// e^{ia} e^{ib} = e^{i(a+b)}, with the unit numbers built from cos and sin
    /// and checked against the angle-addition formulas.
    #[test]
    fn multiplying_unit_numbers_adds_their_angles() {
        for (a, b) in [(0.3, 1.1), (2.0, -0.7), (-1.4, -2.2), (3.0, 3.0)] {
            let product = c(f64::cos(a), f64::sin(a)) * c(f64::cos(b), f64::sin(b));
            assert!(
                close(product, c(f64::cos(a + b), f64::sin(a + b))),
                "{a} {b}"
            );
        }
    }

    /// The spherical harmonics need e^{i m phi} for m = 0..p, built by repeated
    /// multiplication, so check that this stays accurate (De Moivre).
    #[test]
    fn repeated_multiplication_builds_the_angle_multiples() {
        let phi: f64 = 0.7;
        let unit = c(f64::cos(phi), f64::sin(phi));
        let mut power = c(1.0, 0.0);
        for m in 0..=30 {
            let want = c(f64::cos(m as f64 * phi), f64::sin(m as f64 * phi));
            assert!(close(power, want), "m={m}");
            power *= unit;
        }
    }

    // ------------------------------------------------- the assigning forms

    #[test]
    fn the_assigning_operators_agree_with_the_plain_ones() {
        for a in samples() {
            for b in samples() {
                let mut x = a;
                x += b;
                assert!(close(x, a + b), "{a:?} += {b:?}");
                let mut x = a;
                x -= b;
                assert!(close(x, a - b), "{a:?} -= {b:?}");
                let mut x = a;
                x *= b;
                assert!(close(x, a * b), "{a:?} *= {b:?}");
            }
        }
    }

    #[test]
    fn accumulating_with_add_assign_matches_summing() {
        let mut total = Complex::zero();
        for a in samples() {
            total += a;
        }
        let summed = samples().into_iter().fold(Complex::zero(), |s, a| s + a);
        assert!(close(total, summed));
    }

    // ------------------------------------------------------------------ Sum

    #[test]
    fn summing_known_values_gives_the_known_total() {
        let list = [c(1.0, 2.0), c(3.0, -4.0), c(-1.0, 0.0)];
        assert_eq!(list.iter().sum::<Complex>(), c(3.0, -2.0));
        assert_eq!(list.iter().copied().sum::<Complex>(), c(3.0, -2.0));
    }

    #[test]
    fn summing_nothing_gives_zero() {
        let none: [Complex; 0] = [];
        assert_eq!(none.iter().sum::<Complex>(), Complex::zero());
        assert_eq!(none.iter().copied().sum::<Complex>(), Complex::zero());
    }

    #[test]
    fn summing_one_value_gives_that_value() {
        for a in samples() {
            assert_eq!([a].iter().sum::<Complex>(), a);
            assert_eq!([a].iter().copied().sum::<Complex>(), a);
        }
    }

    #[test]
    fn summing_borrowed_and_owned_values_agree_and_match_the_loop() {
        let list = samples();
        let mut looped = Complex::zero();
        for a in &list {
            looped += *a;
        }
        let borrowed: Complex = list.iter().sum();
        let owned: Complex = list.iter().copied().sum();
        assert_eq!(borrowed, owned);
        assert!(close(borrowed, looped));
    }

    /// Both parts are added: a sum that drops or mixes one part cannot match this,
    /// because the real and imaginary totals are different numbers.
    #[test]
    fn summing_keeps_the_real_and_imaginary_totals_separate() {
        let list = samples();
        let want = c(
            list.iter().map(|z| z.re).sum::<f64>(),
            list.iter().map(|z| z.im).sum::<f64>(),
        );
        assert!(close(list.iter().sum(), want));
        assert!(close(list.iter().copied().sum(), want));
    }

    #[test]
    fn the_order_of_the_terms_does_not_matter() {
        let forward: Complex = samples().iter().sum();
        let backward: Complex = samples().iter().rev().sum();
        assert!(close(forward, backward));
    }

    #[test]
    fn a_sum_works_at_the_end_of_an_iterator_chain() {
        // the shape `value` uses: pair things up, multiply, add them all
        let a = samples();
        let b: Vec<Complex> = samples().into_iter().rev().collect();
        let total: Complex = a.iter().zip(b.iter()).map(|(x, y)| *x * *y).sum();
        let mut looped = Complex::zero();
        for i in 0..a.len() {
            looped += a[i] * b[i];
        }
        assert!(close(total, looped));
    }
}
