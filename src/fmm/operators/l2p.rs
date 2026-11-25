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

    if local.order() == 0 {
        // p == 0 => constant field => gradient == 0 => don't do anything to out
        return;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fmm::harmonic::to_index;
    use crate::fmm::operators::{multipole_to_local, particle_to_multipole};
    use crate::fmm::test_util::{bodies_around, pull_from};

    // ---------------------------------------------------------------- helpers

    fn body_at(pos: Vec3) -> Body {
        Body {
            mass: 1.0,
            pos,
            vel: [0.0; 3],
        }
    }

    /// A regular harmonic of order `p` that is zero except for one coefficient.
    fn only(p: usize, n: i64, m: i64, c: Complex) -> SolidHarmonic {
        let mut l = SolidHarmonic::regular(p);
        l.coeffs_mut()[to_index(n, m)] = c;
        l
    }

    /// The acceleration of a single body at `pos` from `field` expanded about `center`.
    fn pull_at(field: &SolidHarmonic, center: Vec3, pos: Vec3) -> Vec3 {
        let body = body_at(pos);
        let mut out = [[0.0; 3]];
        local_to_particle(field, center, &[&body], &mut out);
        out[0]
    }

    fn assert_vec(got: Vec3, want: Vec3, tol: f64, what: &str) {
        for i in 0..3 {
            assert!(
                (got[i] - want[i]).abs() <= tol,
                "{what}: component {i}: got {got:?}, want {want:?}"
            );
        }
    }

    const CENTER: Vec3 = [1.0, 2.0, 3.0];

    // ------------------------------------------------------------- hand cases

    // The acceleration is G times the slope of the field (the real part of the
    // expansion). For single-coefficient fields the slope is a polynomial that
    // can be differentiated by hand. With d = body position - center:
    //   R_1^0 = z                      slope (0, 0, 1)
    //   R_1^1 = -(x + iy)/sqrt 2       so with L = u + iv, Re(L R) has
    //                                  slope (-u, +v, 0)/sqrt 2
    //   R_1^-1 = -(x - iy)/sqrt 2      slope (-u, -v, 0)/sqrt 2
    //   R_2^0 = z^2 - (x^2 + y^2)/2    slope (-x, -y, 2z)

    #[test]
    fn a_linear_z_field_pulls_along_z_everywhere() {
        let field = only(1, 1, 0, Complex { re: 0.7, im: 0.2 });
        for pos in [[1.3, 1.6, 3.4], [0.0, 5.0, -2.0], CENTER] {
            assert_vec(
                pull_at(&field, CENTER, pos),
                [0.0, 0.0, G * 0.7],
                1e-12,
                &format!("body at {pos:?}"),
            );
        }
    }

    #[test]
    fn the_degree_one_order_plus_one_term_pulls_in_the_xy_plane() {
        let (u, v) = (0.7, -0.4);
        let field = only(1, 1, 1, Complex { re: u, im: v });
        let s = 2.0_f64.sqrt();
        assert_vec(
            pull_at(&field, CENTER, [1.3, 1.6, 3.4]),
            [G * -u / s, G * v / s, 0.0],
            1e-12,
            "L_1^1",
        );
    }

    #[test]
    fn the_degree_one_order_minus_one_term_pulls_in_the_xy_plane() {
        let (u, v) = (0.7, -0.4);
        let field = only(1, 1, -1, Complex { re: u, im: v });
        let s = 2.0_f64.sqrt();
        assert_vec(
            pull_at(&field, CENTER, [1.3, 1.6, 3.4]),
            [G * -u / s, G * -v / s, 0.0],
            1e-12,
            "L_1^-1",
        );
    }

    #[test]
    fn a_quadratic_field_pulls_in_proportion_to_the_offset_from_the_center() {
        // slope of u (z^2 - (x^2 + y^2)/2) at offset (dx, dy, dz) is u (-dx, -dy, 2 dz)
        let u = 0.9;
        let field = only(2, 2, 0, Complex { re: u, im: 0.0 });
        let d: Vec3 = [0.3, -0.4, 0.5];
        let pos: Vec3 = [CENTER[0] + d[0], CENTER[1] + d[1], CENTER[2] + d[2]];
        assert_vec(
            pull_at(&field, CENTER, pos),
            [G * -u * d[0], G * -u * d[1], G * 2.0 * u * d[2]],
            1e-12,
            "quadratic",
        );
    }

    #[test]
    fn a_body_at_the_center_of_a_quadratic_field_feels_nothing() {
        let field = only(2, 2, 0, Complex { re: 0.9, im: 0.0 });
        assert_vec(
            pull_at(&field, CENTER, CENTER),
            [0.0; 3],
            1e-12,
            "at center",
        );
    }

    #[test]
    fn a_constant_field_pulls_nothing() {
        let field = only(3, 0, 0, Complex { re: 2.5, im: -1.0 });
        assert_vec(
            pull_at(&field, CENTER, [1.4, 1.7, 3.2]),
            [0.0; 3],
            1e-12,
            "constant",
        );
    }

    #[test]
    fn an_order_zero_field_pulls_nothing() {
        // There is no degree-1 part to read, which must mean zero, not a crash.
        let field = only(0, 0, 0, Complex { re: 2.5, im: 0.0 });
        assert_vec(
            pull_at(&field, CENTER, [1.4, 1.7, 3.2]),
            [0.0; 3],
            1e-12,
            "order 0",
        );
    }

    #[test]
    fn the_pull_is_the_sum_of_the_pulls_of_the_terms() {
        let p = 3;
        let mut sum = only(p, 1, 0, Complex { re: 0.7, im: 0.0 });
        sum += &only(p, 2, 0, Complex { re: 0.9, im: 0.0 });
        sum += &only(p, 1, 1, Complex { re: 0.3, im: -0.2 });
        let pos: Vec3 = [1.3, 1.6, 3.4];

        let parts = [
            only(p, 1, 0, Complex { re: 0.7, im: 0.0 }),
            only(p, 2, 0, Complex { re: 0.9, im: 0.0 }),
            only(p, 1, 1, Complex { re: 0.3, im: -0.2 }),
        ];
        let mut want = [0.0; 3];
        for part in &parts {
            let a = pull_at(part, CENTER, pos);
            for i in 0..3 {
                want[i] += a[i];
            }
        }
        assert_vec(pull_at(&sum, CENTER, pos), want, 1e-12, "sum of terms");
    }

    // ------------------------------------------------------------ bookkeeping

    #[test]
    fn each_body_gets_its_own_answer_in_its_own_slot() {
        let field = only(2, 2, 0, Complex { re: 0.9, im: 0.0 });
        let near = body_at([1.1, 2.0, 3.0]);
        let far = body_at([1.0, 2.0, 3.5]);
        let mut out = [[0.0; 3]; 2];
        local_to_particle(&field, CENTER, &[&near, &far], &mut out);
        assert_vec(out[0], pull_at(&field, CENTER, near.pos), 1e-12, "slot 0");
        assert_vec(out[1], pull_at(&field, CENTER, far.pos), 1e-12, "slot 1");
        // non-vacuity: the two bodies really do feel different pulls
        assert!((out[0][0] - out[1][0]).abs() > 1e-3 || (out[0][2] - out[1][2]).abs() > 1e-3);
    }

    #[test]
    fn the_pull_is_added_to_what_the_slot_already_holds() {
        let field = only(1, 1, 0, Complex { re: 0.7, im: 0.0 });
        let body = body_at([1.3, 1.6, 3.4]);
        let mut out = [[1.0, -2.0, 3.0]];
        local_to_particle(&field, CENTER, &[&body], &mut out);
        assert_vec(out[0], [1.0, -2.0, 3.0 + G * 0.7], 1e-12, "accumulated");
    }

    #[test]
    fn no_bodies_is_a_no_op() {
        let field = only(1, 1, 0, Complex { re: 0.7, im: 0.0 });
        local_to_particle(&field, CENTER, &[], &mut []);
    }

    #[test]
    #[should_panic(expected = "L2P requires regular harmonic")]
    fn l2p_refuses_an_irregular_input() {
        let irregular = SolidHarmonic::irregular(2);
        let body = body_at([1.3, 1.6, 3.4]);
        local_to_particle(&irregular, CENTER, &[&body], &mut [[0.0; 3]]);
    }

    #[test]
    #[should_panic(expected = "assertion")]
    fn l2p_refuses_a_slot_list_of_the_wrong_length() {
        let field = only(1, 1, 0, Complex { re: 0.7, im: 0.0 });
        let body = body_at([1.3, 1.6, 3.4]);
        local_to_particle(&field, CENTER, &[&body], &mut [[0.0; 3], [0.0; 3]]);
    }

    // ---------------------------------------------------------------- physics

    // The way the tree uses it: a far group of masses is summarized (P2M),
    // converted to a field about the leaf (M2L), and the leaf's bodies read
    // their pull from it. That must be Newton's pull of the far masses, with
    // the law written out independently in `pull_from`.

    #[test]
    fn the_pull_from_a_far_group_matches_newtons_law() {
        let p = 10;
        let source_center: Vec3 = [10.0, 0.0, 0.0];
        let sources = bodies_around(source_center, 1.0, 6, 5);
        let mp = particle_to_multipole(sources.iter(), source_center, p);

        let leaf_center: Vec3 = [0.0, 0.0, 0.0];
        let field = multipole_to_local(&mp, [-10.0, 0.0, 0.0]);

        let targets = bodies_around(leaf_center, 0.5, 5, 9);
        let refs: Vec<&Body> = targets.iter().collect();
        let mut got = vec![[0.0; 3]; refs.len()];
        local_to_particle(&field, leaf_center, &refs, &mut got);

        for (target, g) in targets.iter().zip(&got) {
            let mut want = [0.0; 3];
            for s in &sources {
                let a = pull_from(target.pos, s);
                for i in 0..3 {
                    want[i] += a[i];
                }
            }
            let size = want.iter().map(|c| c * c).sum::<f64>().sqrt();
            assert!(size > 1e-3, "fixture too weak to mean anything: {size}");
            assert_vec(*g, want, 1e-4 * size, &format!("body at {:?}", target.pos));
        }
    }
}
