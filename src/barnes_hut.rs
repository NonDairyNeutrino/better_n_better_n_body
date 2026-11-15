use crate::shared::{Body, ORIGIN};

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
