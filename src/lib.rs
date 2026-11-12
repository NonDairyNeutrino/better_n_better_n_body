//! This is a module documentation

pub const G: f64 = 6.67 * 1e-11;

/// This is the docstring for Object
pub struct Object {
    pub mass: f64,
    pub pos: [f64; 3],
    pub vel: [f64; 3],
}
impl Object {
    /// Creates a new physical Object
    pub fn new() -> Object {
        return Object {
            pos: [0f64, 0f64, 0f64],
            vel: [0f64, 0f64, 0f64],
            mass: 1f64,
        };
    }
    pub fn distance(self, pos: [f64; 3]) -> f64 {
        self.pos
            .iter()
            .zip(pos)
            .map(|pair| (*pair.0 - pair.1).powi(2))
            .sum()
    }
}
