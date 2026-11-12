use crate::lib::{G, Object};

pub fn gravity<T>(object: Object, pos: [f64; 3]) -> f64 {
    return G * object.mass / object.distance(pos).powi(2);
}
