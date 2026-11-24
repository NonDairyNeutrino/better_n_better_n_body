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
