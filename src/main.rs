use crate::lib::Object;
mod direct;
mod lib;

fn main() {
    let ob: Object = Object::new();
    let g: f64 = direct::gravity::<f64>(ob, [1f64, 0f64, 0f64]);
    println!("{g}")
}
