use std::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}
impl Complex {
    pub fn zero() -> Complex {
        Complex {
            re: 0.0_f64,
            im: 0.0_f64,
        }
    }
    pub fn conj(self) -> Complex {
        Complex {
            re: self.re,
            im: -self.im,
        }
    }
}

impl Add for Complex {
    type Output = Complex;
    fn add(self, rhs: Self) -> Self::Output {
        Complex {
            re: self.re + rhs.re,
            im: self.im + rhs.im,
        }
    }
}

impl Sub for Complex {
    type Output = Complex;
    fn sub(self, rhs: Self) -> Self::Output {
        Complex {
            re: self.re - rhs.re,
            im: self.im - rhs.im,
        }
    }
}
impl Mul<Complex> for Complex {
    type Output = Complex;
    fn mul(self, rhs: Self) -> Self::Output {
        let a = self.re;
        let b = self.im;
        let c = rhs.re;
        let d = rhs.im;
        Complex {
            re: a * c - b * d,
            im: a * d + b * c,
        }
    }
}

impl Mul<f64> for Complex {
    type Output = Complex;
    fn mul(self, rhs: f64) -> Self::Output {
        Complex {
            re: rhs * self.re,
            im: rhs * self.im,
        }
    }
}

impl Neg for Complex {
    type Output = Complex;
    fn neg(self) -> Self::Output {
        Complex {
            re: -self.re,
            im: -self.im,
        }
    }
}

impl AddAssign for Complex {
    fn add_assign(&mut self, rhs: Self) {
        self.re += rhs.re;
        self.im += rhs.im;
    }
}
impl SubAssign for Complex {
    fn sub_assign(&mut self, rhs: Self) {
        self.re -= rhs.re;
        self.im -= rhs.im;
    }
}
impl MulAssign for Complex {
    fn mul_assign(&mut self, rhs: Self) {
        self.re *= rhs.re;
        self.im *= rhs.im;
    }
}
