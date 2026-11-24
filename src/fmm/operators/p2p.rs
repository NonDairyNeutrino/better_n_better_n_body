use crate::{
    direct::calculate_acceleration_from_system,
    shared::{Body, Vec3},
};

pub fn particle_to_particle(targets: &[&Body], sources: &[&Body], out: &mut [Vec3]) {
    assert_eq!(
        targets.len(),
        out.len(),
        "accumulation list must be same length as target list"
    );
    for (target, slot) in targets.iter().zip(out.iter_mut()) {
        let acc: Vec3 = calculate_acceleration_from_system(sources.iter().copied(), &target.pos);
        for (slot_elem, ae) in slot.iter_mut().zip(acc) {
            *slot_elem += ae;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        fmm::test_util::{bodies_around, pull_from},
        shared::{G, test_util::body_with},
    };

    /// `a` and `b` agree component by component to a relative 1e-12.
    fn same(a: Vec3, b: Vec3, what: &str) {
        for d in 0..3 {
            assert!(
                (a[d] - b[d]).abs() <= 1e-12 * (1.0 + b[d].abs()),
                "{what}: component {d}: got {}, want {}",
                a[d],
                b[d]
            );
        }
    }

    /// A mass 3 at distance 2 on the x-axis pulls with G * 3 / 4 toward itself.
    #[test]
    fn one_source_pulls_with_g_m_over_r_squared_toward_itself() {
        let target = body_with(1.0, 0.0, 0.0, 0.0);
        let source = body_with(3.0, 2.0, 0.0, 0.0);
        let mut out = [[0.0; 3]];
        particle_to_particle(&[&target], &[&source], &mut out);
        same(out[0], [G * 0.75, 0.0, 0.0], "pull on the target");
    }

    /// Mass 2 at (1, 2, 2), distance 3: the pull is G * 2 / 27 times (1, 2, 2).
    #[test]
    fn a_source_off_the_axes_pulls_along_the_line_to_it() {
        let target = body_with(1.0, 0.0, 0.0, 0.0);
        let source = body_with(2.0, 1.0, 2.0, 2.0);
        let mut out = [[0.0; 3]];
        particle_to_particle(&[&target], &[&source], &mut out);
        let k = G * 2.0 / 27.0;
        same(out[0], [k, 2.0 * k, 2.0 * k], "pull on the target");
    }

    /// Newton's third law: the force (mass times acceleration) one body feels
    /// from the other is equal and opposite.
    #[test]
    fn the_forces_between_two_bodies_are_equal_and_opposite() {
        let a = body_with(1.5, 0.3, -0.2, 0.9);
        let b = body_with(4.0, -1.1, 0.8, 0.1);
        let mut on_a = [[0.0; 3]];
        let mut on_b = [[0.0; 3]];
        particle_to_particle(&[&a], &[&b], &mut on_a);
        particle_to_particle(&[&b], &[&a], &mut on_b);
        let (fa, fb): (Vec3, Vec3) = (
            std::array::from_fn(|i| a.mass * on_a[0][i]),
            std::array::from_fn(|i| b.mass * on_b[0][i]),
        );
        assert!(
            fa.iter().any(|c| c.abs() > 1e-3),
            "forces too small to test"
        );
        same(
            fa,
            std::array::from_fn(|i| -fb[i]),
            "force on a vs force on b",
        );
    }

    /// With every body as both target and source, each body gets the sum of the
    /// pulls of all the others, written out with Newton's law.
    #[test]
    fn with_all_bodies_as_targets_and_sources_each_gets_the_pull_of_the_others() {
        let bodies = bodies_around([0.0; 3], 1.0, 12, 5);
        let all: Vec<&Body> = bodies.iter().collect();
        let mut out = vec![[0.0; 3]; all.len()];
        particle_to_particle(&all, &all, &mut out);
        for (i, target) in bodies.iter().enumerate() {
            let mut want = [0.0; 3];
            for (j, source) in bodies.iter().enumerate() {
                if i != j {
                    let pull = pull_from(target.pos, source);
                    for d in 0..3 {
                        want[d] += pull[d];
                    }
                }
            }
            same(out[i], want, &format!("body {i}"));
        }
        // each target gets its own pull, not a shared one
        assert!((out[0][0] - out[1][0]).abs() > 1e-6, "slots look identical");
    }

    /// A body that is both target and source must not pull on itself, which
    /// would be a division by zero.
    #[test]
    fn a_body_does_not_pull_on_itself() {
        let body = body_with(2.0, 0.4, 0.5, 0.6);
        let mut out = [[0.0; 3]];
        particle_to_particle(&[&body], &[&body], &mut out);
        assert_eq!(out[0], [0.0; 3]);
    }

    /// The result is added to what `out` already holds, since a leaf collects
    /// from several source leaves.
    #[test]
    fn the_pull_is_added_to_what_the_output_already_holds() {
        let target = body_with(1.0, 0.0, 0.0, 0.0);
        let source = body_with(3.0, 2.0, 0.0, 0.0);
        let mut out = [[1.0, 2.0, 3.0]];
        particle_to_particle(&[&target], &[&source], &mut out);
        same(out[0], [1.0 + G * 0.75, 2.0, 3.0], "accumulated");
    }

    /// Two calls with two source lists give the same as one call with both lists.
    #[test]
    fn pulls_from_two_source_lists_add_up_to_the_pull_from_both() {
        let bodies = bodies_around([0.0; 3], 1.0, 9, 8);
        let target = &bodies[0];
        let first: Vec<&Body> = bodies[1..5].iter().collect();
        let second: Vec<&Body> = bodies[5..].iter().collect();
        let both: Vec<&Body> = bodies[1..].iter().collect();
        let mut in_two = [[0.0; 3]];
        particle_to_particle(&[target], &first, &mut in_two);
        particle_to_particle(&[target], &second, &mut in_two);
        let mut in_one = [[0.0; 3]];
        particle_to_particle(&[target], &both, &mut in_one);
        same(in_two[0], in_one[0], "two calls vs one");
    }

    #[test]
    fn no_sources_leave_the_output_unchanged() {
        let target = body_with(1.0, 0.0, 0.0, 0.0);
        let mut out = [[0.5, -0.5, 2.0]];
        particle_to_particle(&[&target], &[], &mut out);
        assert_eq!(out[0], [0.5, -0.5, 2.0]);
    }

    #[test]
    #[should_panic(expected = "same length")]
    fn an_output_list_of_the_wrong_length_is_refused() {
        let a = body_with(1.0, 0.0, 0.0, 0.0);
        let b = body_with(1.0, 1.0, 0.0, 0.0);
        let mut out = [[0.0; 3]];
        particle_to_particle(&[&a, &b], &[&b], &mut out);
    }
}
