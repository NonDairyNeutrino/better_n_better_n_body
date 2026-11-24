use crate::{
    fmm::{complex::Complex, harmonic::SolidHarmonic},
    shared::Vec3,
};

/// Turn the mass distribution of a source cell into the external potential at a target cell.
pub fn multipole_to_local(mp: &SolidHarmonic, shift: Vec3) -> SolidHarmonic {
    if let SolidHarmonic::Irregular(_) = mp {
        panic!("Cannot convert Irregular harmonic to Regular")
    }

    let table: Vec<Complex> = irregular_table(mp.order(), shift);
    let mut local_potential: SolidHarmonic = SolidHarmonic::regular(mp.order());
    convert_moments(mp, &mut local_potential, &table);
    return local_potential;
}

fn irregular_table(order: usize, shift: Vec3) -> Vec<Complex> {
    todo!()
}

fn convert_moments(mp: &SolidHarmonic, local: &mut SolidHarmonic, table: &[Complex]) {
    todo!()
}

fn get_local_coefficient(mp: &SolidHarmonic, table: &[Complex], j: i64, k: i64) -> Complex {
    todo!()
}

fn get_conversion_factor(n: i64, m: i64, j: i64, k: i64) -> f64 {
    todo!()
}

fn conversion_sign(j: i64, k: i64, m: i64) -> f64 {
    todo!()
}
