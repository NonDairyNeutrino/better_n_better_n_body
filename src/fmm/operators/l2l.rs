use crate::fmm::complex::Complex;
use crate::fmm::harmonic::{SolidHarmonic, to_index};
use crate::fmm::operators::m2m::{binomial, epsilon, shift_scale_harmonics};
use crate::shared::Vec3;

pub fn local_to_local(local: &SolidHarmonic, shift: Vec3) -> SolidHarmonic {
    if let SolidHarmonic::Irregular(_) = local {
        panic!("L2L requires regular harmonic as input")
    }

    let p = local.order();
    let table: Vec<Complex> = shift_scale_harmonics(p, shift);
    let mut res: SolidHarmonic = SolidHarmonic::regular(p);
    let coeffs: &mut [Complex] = res.coeffs_mut();

    for j in 0..=p {
        let ji: i64 = j as i64;
        for k in -ji..=ji {
            let idx: usize = to_index(ji, k);
            coeffs[idx] = local_coefficient(local, &table, ji, k);
        }
    }
    return res;
}

fn local_coefficient(local: &SolidHarmonic, table: &[Complex], j: i64, k: i64) -> Complex {
    let mut sum: Complex = Complex::zero();
    let p: i64 = local.order() as i64;
    for lambda in 0..=(p - j) {
        let li = lambda as i64;
        for mu in -li..=li {
            let n = j + li;
            let m = k + mu;
            let coeffs: &[Complex] = local.coeffs();
            sum += table[to_index(li, mu)] * coeffs[to_index(n, m)] * weight(n, m, li, mu);
        }
    }
    return sum;
}

fn weight(n: i64, m: i64, lambda: i64, mu: i64) -> f64 {
    let ep: f64 = epsilon(m, mu);
    let binp: f64 = binomial(n + m, lambda + mu);
    let binm: f64 = binomial(n - m, lambda - mu);
    return ep * (binp * binm).sqrt();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fmm::complex::Complex;
    use crate::fmm::harmonic::to_index;
    use crate::fmm::operators::m2l::M2LScalings;
    use crate::fmm::operators::{multipole_to_local, particle_to_multipole};
    use crate::fmm::test_util::{assert_close, bodies_around, direct_potential, distinct};

    // ---------------------------------------------------------------- helpers

    /// A regular harmonic of order `p` holding `distinct(k, offset)` in slot `k`.
    fn regular_moments(p: usize, offset: f64) -> SolidHarmonic {
        let mut l = SolidHarmonic::regular(p);
        for (k, c) in l.coeffs_mut().iter_mut().enumerate() {
            *c = distinct(k, offset);
        }
        l
    }

    /// A regular harmonic of order `p` that is zero except for one coefficient.
    fn only(p: usize, n: i64, m: i64, c: Complex) -> SolidHarmonic {
        let mut l = SolidHarmonic::regular(p);
        l.coeffs_mut()[to_index(n, m)] = c;
        l
    }

    fn max_diff(a: &SolidHarmonic, b: &SolidHarmonic) -> f64 {
        assert_eq!(a.order(), b.order(), "orders differ");
        a.coeffs()
            .iter()
            .zip(b.coeffs())
            .map(|(x, y)| (x.re - y.re).abs().max((x.im - y.im).abs()))
            .fold(0.0, f64::max)
    }

    fn sub(a: Vec3, b: Vec3) -> Vec3 {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }

    // ------------------------------------------------------------------ shape

    #[test]
    fn the_result_is_regular_and_has_the_order_of_the_input() {
        for p in [0, 1, 3, 5] {
            let out = local_to_local(&regular_moments(p, 0.4), [0.2, -0.3, 0.1]);
            assert!(
                matches!(out, SolidHarmonic::Regular(_)),
                "p = {p}: not Regular"
            );
            assert_eq!(out.order(), p);
            assert_eq!(out.coeffs().len(), (p + 1) * (p + 1));
        }
    }

    #[test]
    #[should_panic(expected = "regular harmonic")]
    fn l2l_refuses_an_irregular_input() {
        let irregular = SolidHarmonic::irregular(3);
        local_to_local(&irregular, [0.1, 0.2, 0.3]);
    }

    // ------------------------------------------------------------- hand cases

    // A field that is a single coefficient is a simple polynomial, so shifting
    // its center can be expanded by hand. With z' = z - t (the new height):
    //   R_2^0 = z^2 - (x^2 + y^2)/2   becomes  R_2^0' + 2 t z' + t^2
    //   R_1^0 = z                     becomes  z' + t
    // so the coefficients of the new field are read straight off.

    #[test]
    fn a_pure_quadratic_shifted_along_z_matches_the_hand_expansion() {
        let t = 0.5;
        let c = Complex { re: 1.3, im: -0.4 };
        let out = local_to_local(&only(2, 2, 0, c), [0.0, 0.0, t]);

        // degree 2 is unchanged, degree 1 gains 2 t c, degree 0 gains t^2 c
        let want = |n: i64, m: i64| match (n, m) {
            (2, 0) => c,
            (1, 0) => Complex {
                re: 2.0 * t * c.re,
                im: 2.0 * t * c.im,
            },
            (0, 0) => Complex {
                re: t * t * c.re,
                im: t * t * c.im,
            },
            _ => Complex { re: 0.0, im: 0.0 },
        };
        for n in 0..=2 {
            for m in -n..=n {
                assert_close(out.get(n, m), want(n, m), 1e-12, &format!("L'_{n}^{m}"));
            }
        }
    }

    #[test]
    fn a_pure_linear_term_picks_up_a_constant_from_a_z_shift() {
        let t = -0.8;
        let c = Complex { re: 0.7, im: 0.2 };
        let out = local_to_local(&only(1, 1, 0, c), [0.0, 0.0, t]);

        assert_close(out.get(1, 0), c, 1e-12, "L'_1^0");
        let constant = Complex {
            re: t * c.re,
            im: t * c.im,
        };
        assert_close(out.get(0, 0), constant, 1e-12, "L'_0^0");
        for m in [-1, 1] {
            let zero = Complex { re: 0.0, im: 0.0 };
            assert_close(out.get(1, m), zero, 1e-12, &format!("L'_1^{m}"));
        }
    }

    #[test]
    fn a_shift_across_the_gradient_leaves_a_linear_field_alone() {
        // The field c z does not depend on x or y, so moving its center
        // sideways cannot change it.
        let c = Complex { re: 0.7, im: 0.2 };
        let field = only(1, 1, 0, c);
        let out = local_to_local(&field, [0.6, -0.9, 0.0]);
        assert!(max_diff(&out, &field) < 1e-12);
    }

    #[test]
    fn order_zero_is_a_constant_and_no_shift_changes_it() {
        let field = only(0, 0, 0, Complex { re: 2.5, im: -1.0 });
        let out = local_to_local(&field, [1.0, -2.0, 3.0]);
        assert!(max_diff(&out, &field) < 1e-12);
    }

    // ------------------------------------------------------------- invariants

    #[test]
    fn a_zero_shift_changes_nothing() {
        let field = regular_moments(5, 0.3);
        let out = local_to_local(&field, [0.0, 0.0, 0.0]);
        assert!(max_diff(&out, &field) < 1e-12);
    }

    #[test]
    fn the_top_degree_is_unchanged_by_any_shift() {
        // The highest-degree terms cannot absorb anything from above them.
        let p = 4;
        let field = regular_moments(p, 0.6);
        let out = local_to_local(&field, [0.3, -0.5, 0.4]);
        for m in -(p as i64)..=(p as i64) {
            assert_close(
                out.get(p as i64, m),
                field.get(p as i64, m),
                1e-12,
                &format!("degree {p}, order {m}"),
            );
        }
    }

    #[test]
    fn a_nonzero_shift_changes_the_lower_degrees() {
        // Guard against an implementation that returns its input.
        let field = regular_moments(4, 0.6);
        let out = local_to_local(&field, [0.3, -0.5, 0.4]);
        assert!(max_diff(&out, &field) > 0.1);
    }

    // Moving the center must not move the field. Evaluating the new expansion
    // at (x - t) must give what the old one gives at x. This is the defining
    // property of L2L, and it fails if the shift has the wrong sign, a
    // coefficient is mispaired, or a term is dropped.

    #[test]
    fn the_shifted_field_agrees_with_the_original_at_every_point() {
        let p = 4;
        let field = regular_moments(p, 0.5);
        let shifts: [Vec3; 5] = [
            [0.4, 0.0, 0.0],
            [0.0, -0.3, 0.0],
            [0.0, 0.0, 0.5],
            [0.3, -0.2, 0.25],
            [-0.4, 0.5, -0.35],
        ];
        let points: [Vec3; 4] = [
            [0.1, 0.2, -0.1],
            [-0.3, 0.25, 0.4],
            [0.5, -0.45, 0.05],
            [0.0, 0.0, 0.0],
        ];
        for t in shifts {
            let out = local_to_local(&field, t);
            for x in points {
                let old = field.value(x);
                let new = out.value(sub(x, t));
                assert_close(new, old, 1e-9, &format!("shift {t:?} at {x:?}"));
            }
        }
    }

    #[test]
    fn shifting_there_and_back_returns_the_input() {
        let field = regular_moments(5, 0.2);
        let t: Vec3 = [0.3, -0.6, 0.2];
        let back = local_to_local(&local_to_local(&field, t), [-t[0], -t[1], -t[2]]);
        assert!(max_diff(&back, &field) < 1e-9);
    }

    #[test]
    fn two_shifts_equal_one_shift_by_their_sum() {
        let field = regular_moments(5, 0.2);
        let a: Vec3 = [0.3, -0.1, 0.2];
        let b: Vec3 = [-0.2, 0.4, 0.15];
        let both = [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
        let two_steps = local_to_local(&local_to_local(&field, a), b);
        let one_step = local_to_local(&field, both);
        assert!(max_diff(&two_steps, &one_step) < 1e-9);
    }

    #[test]
    fn shifting_a_sum_gives_the_sum_of_the_shifts() {
        let p = 4;
        let t: Vec3 = [0.25, 0.1, -0.3];

        let mut sum = regular_moments(p, 0.3);
        sum += &regular_moments(p, -0.8);
        let shifted_sum = local_to_local(&sum, t);

        let mut sum_of_shifted = local_to_local(&regular_moments(p, 0.3), t);
        sum_of_shifted += &local_to_local(&regular_moments(p, -0.8), t);

        assert!(max_diff(&shifted_sum, &sum_of_shifted) < 1e-9);
    }

    // ---------------------------------------------------------------- physics

    // The way the tree uses it: a far source group creates an ambient field
    // about a parent cell (M2L); the child inherits it (L2L). Inside the child
    // the inherited field must be the same field, and must still be the true
    // potential of the sources up to the series truncation.

    #[test]
    fn an_inherited_far_field_is_the_same_field_and_the_direct_potential() {
        let p = 10;
        let source_center: Vec3 = [10.0, 0.0, 0.0];
        let sources = bodies_around(source_center, 1.0, 6, 5);
        let mp = particle_to_multipole(sources.iter(), source_center, p);

        let parent: Vec3 = [0.0, 0.0, 0.0];
        let local = multipole_to_local(&mp, sub(parent, source_center), &M2LScalings::new(p));

        let child: Vec3 = [0.5, 0.5, -0.5];
        let inherited = local_to_local(&local, sub(child, parent));

        let offsets: [Vec3; 4] = [
            [0.2, -0.3, 0.1],
            [-0.4, 0.4, 0.45],
            [0.0, 0.0, 0.0],
            [0.45, 0.1, -0.3],
        ];
        for off in offsets {
            let x: Vec3 = [child[0] + off[0], child[1] + off[1], child[2] + off[2]];
            let from_parent = local.value(sub(x, parent));
            let from_child = inherited.value(off);

            // Same field: shifting loses nothing.
            assert_close(
                from_child,
                from_parent,
                1e-9,
                &format!("same field at {x:?}"),
            );

            // And it is still the true potential, up to truncation.
            let truth = direct_potential(&sources, x);
            assert!(
                (from_child.re - truth).abs() < 1e-4 * truth.abs(),
                "at {x:?}: got {}, direct sum {truth}",
                from_child.re
            );
            assert!(from_child.im.abs() < 1e-9, "imaginary part at {x:?}");
        }
    }
}
