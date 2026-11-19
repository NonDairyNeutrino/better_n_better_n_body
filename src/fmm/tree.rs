use crate::shared::System;

#[derive(Debug, Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
}
impl Complex {
    fn zero() -> Complex {
        Complex {
            re: 0.0_f64,
            im: 0.0_f64,
        }
    }
}

#[derive(Debug, Clone)]
enum SolidHarmonic {
    Regular(Vec<Complex>),
    Irregular(Vec<Complex>),
}
impl SolidHarmonic {
    fn regular(p: usize) -> SolidHarmonic {
        let coeffs: Vec<Complex> = vec![Complex::zero(); (p + 1).pow(2)];
        return SolidHarmonic::Regular(coeffs);
    }
    fn irregular(p: usize) -> SolidHarmonic {
        let coeffs: Vec<Complex> = vec![Complex::zero(); (p + 1).pow(2)];
        return SolidHarmonic::Irregular(coeffs);
    }
    fn coeffs(&self) -> &[Complex] {
        match self {
            SolidHarmonic::Regular(v) | SolidHarmonic::Irregular(v) => v,
        }
    }
    fn coeffs_mut(&mut self) -> &mut [Complex] {
        match self {
            SolidHarmonic::Regular(v) | SolidHarmonic::Irregular(v) => v,
        }
    }
    fn get(&self, n: i64, m: i64) -> Complex {
        assert!(0 <= n);
        assert!(-n <= m && m <= n);
        let idx: usize = (n.pow(2) + n + m) as usize;
        let coeff: Complex = self.coeffs()[idx];
        return coeff;
    }
}

const MAX_LEAF_BODIES: usize = 2;
pub struct FmmTree<'a> {
    system: System<'a>,
    far_field: SolidHarmonic,   // the field far from the observer
    local_field: SolidHarmonic, // the field near the observer
    children: Vec<FmmTree<'a>>,
}
impl<'a> FmmTree<'a> {
    pub fn new(system: System<'a>, p: usize) -> FmmTree<'a> {
        let far_field: SolidHarmonic = SolidHarmonic::irregular(p);
        let local_field: SolidHarmonic = SolidHarmonic::regular(p);

        if MAX_LEAF_BODIES < system.bodies.len() {
            // split
            let new_system: System<'_> = System::new(system.center, system.hw);
            let mut children: Vec<FmmTree<'a>> = Vec::with_capacity(8);
            let systems: [System; 8] = system.partition();
            for child in systems {
                let child_tree: FmmTree = FmmTree::new(child, p);
                children.push(child_tree);
            }

            return FmmTree {
                system: new_system,
                children,
                far_field,
                local_field,
            };
        } else {
            // don't split
            let children: Vec<FmmTree<'a>> = Vec::new();

            return FmmTree {
                system,
                children,
                far_field,
                local_field,
            };
        }
    }
}
