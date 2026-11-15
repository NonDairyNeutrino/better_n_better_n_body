use crate::shared::{Body, ORIGIN};

fn find_center(system: &[Body]) -> [f64; 3] {
    let mut center: [f64; 3] = ORIGIN;
    let n_bodies: f64 = system.len() as f64;

    for obj in system {
        for dim in 0..3 {
            center[dim] += obj.pos[dim];
        }
    }

    // average
    for dim in 0..center.len() {
        center[dim] /= n_bodies;
    }
    return center;
}

fn find_half_width(system: &[Body], center: [f64; 3]) -> f64 {
    let mut half_width: f64 = 0.0_f64;
    for obj in system {
        let d: f64 = obj.distance(&center);
        if d > half_width {
            half_width = d;
        }
    }
    return half_width;
}

/// Find the center and half-width of a cube that contains the system.
pub fn find_bounding_box(system: &[Body]) -> ([f64; 3], f64) {
    let center: [f64; 3] = find_center(system);
    let width: f64 = find_half_width(system, center);
    return (center, width);
}

#[derive(Debug, Clone)]
enum Node {
    Body0 {
        center: [f64; 3],
        hw: f64,
    },
    Body1 {
        center: [f64; 3],
        hw: f64,
        obj: Body,
    },
    BodyN {
        center: [f64; 3],
        hw: f64,
        mass: f64,
        com: [f64; 3],
        children: Box<[Node; 8]>,
    },
}
impl Node {
    fn new(center: [f64; 3], half_width: f64) -> Node {
        return Node::Body0 {
            center: center,
            hw: half_width,
        };
    }
    fn insert(&mut self, new_body: Body) {
        match self {
            Node::Body0 { center, hw } => {
                *self = Node::Body1 {
                    center: *center,
                    hw: *hw,
                    obj: new_body,
                }
            }
            Node::Body1 { center, hw, obj } => {
                *self = Node::BodyN {
                    center: *center,
                    hw: *hw,
                    mass: obj.mass + new_body.mass,
                    com: todo!(),
                    children: todo!(),
                }
            }
            Node::BodyN {
                center,
                hw,
                mass,
                com,
                children,
            } => {
                todo!()
            }
        }
    }
}

fn child_index(center: [f64; 3], pos: [f64; 3]) -> usize {
    todo!()
}
fn child_center(center: [f64; 3], half_width: f64, index: usize) -> [f64; 3] {
    todo!()
}
