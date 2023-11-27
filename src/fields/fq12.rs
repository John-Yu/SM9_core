#![allow(dead_code)]

use core::ops::{Add, Mul, MulAssign, Neg, Sub};
use rand::Rng;

use crate::fields::{FieldElement, Fq4};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq12 {
    pub(crate) c0: Fq4,
    pub(crate) c1: Fq4,
    pub(crate) c2: Fq4,
}

impl Fq12 {
    pub fn new(c0: Fq4, c1: Fq4, c2: Fq4) -> Self {
        Fq12 { c0, c1, c2 }
    }

    pub fn mul_by_nonresidue(&self) -> Self {
        Fq12 {
            c0: self.c2.mul_by_nonresidue(),
            c1: self.c0,
            c2: self.c1,
        }
    }

    pub fn scale(&self, by: &Fq4) -> Self {
        Fq12 {
            c0: self.c0 * *by,
            c1: self.c1 * *by,
            c2: self.c2 * *by,
        }
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        match power {
            1 => Fq12 {
                c0: self.c0.frobenius_map(10), // case 1 for c0 , so call  10. and so on
                c1: self.c1.frobenius_map(11),
                c2: self.c2.frobenius_map(12),
            },
            2 => Fq12 {
                c0: self.c0.unitary_inverse(),
                c1: self.c1.frobenius_map(21), // case 2 for c1 , so call  21
                c2: self.c2.frobenius_map(22),
            },

            3 => Fq12 {
                c0: self.c0.frobenius_map(30),
                c1: self.c1.frobenius_map(31),
                c2: self.c2.frobenius_map(32),
            },
            6 => Fq12 {
                c0: self.c0.unitary_inverse(),
                c1: -self.c1.unitary_inverse(),
                c2: self.c2.unitary_inverse(),
            },
            _ => unimplemented!(),
        }
    }
    /// Converts an element into a byte representation in
    /// big-endian byte order.
    pub fn to_slice(&self) -> [u8; 384] {
        let mut res = [0u8; 384];
        let b2 = self.c2.to_slice();
        let b1 = self.c1.to_slice();
        let b0 = self.c0.to_slice();
        res[..128].copy_from_slice(&b2);
        res[128..256].copy_from_slice(&b1);
        res[256..].copy_from_slice(&b0);
        res
    }
}

impl FieldElement for Fq12 {
    fn zero() -> Self {
        Fq12 {
            c0: Fq4::zero(),
            c1: Fq4::zero(),
            c2: Fq4::zero(),
        }
    }

    fn one() -> Self {
        Fq12 {
            c0: Fq4::one(),
            c1: Fq4::zero(),
            c2: Fq4::zero(),
        }
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        Fq12 {
            c0: Fq4::random(rng),
            c1: Fq4::random(rng),
            c2: Fq4::random(rng),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }

    fn squared(&self) -> Self {
        let s0 = self.c0.squared();
        let ab = self.c0 * self.c1;
        let s1 = ab + ab;
        let s2 = (self.c0 - self.c1 + self.c2).squared();
        let bc = self.c1 * self.c2;
        let s3 = bc + bc;
        let s4 = self.c2.squared();

        Fq12 {
            c0: s0 + s3.mul_by_nonresidue(),
            c1: s1 + s4.mul_by_nonresidue(),
            c2: s1 + s2 + s3 - s0 - s4,
        }
    }

    fn inverse(self) -> Option<Self> {
        let c0 = self.c0.squared() - self.c1 * self.c2.mul_by_nonresidue();
        let c1 = self.c2.squared().mul_by_nonresidue() - self.c0 * self.c1;
        let c2 = self.c1.squared() - self.c0 * self.c2;
        match ((self.c2 * c1 + self.c1 * c2).mul_by_nonresidue() + self.c0 * c0).inverse() {
            Some(t) => Some(Fq12 {
                c0: t * c0,
                c1: t * c1,
                c2: t * c2,
            }),
            None => None,
        }
    }
}

impl Mul for Fq12 {
    type Output = Fq12;

    fn mul(self, other: Fq12) -> Fq12 {
        let a_a = self.c0 * other.c0;
        let b_b = self.c1 * other.c1;
        let c_c = self.c2 * other.c2;

        Fq12 {
            c0: ((self.c1 + self.c2) * (other.c1 + other.c2) - b_b - c_c).mul_by_nonresidue() + a_a,
            c1: (self.c0 + self.c1) * (other.c0 + other.c1) - a_a - b_b + c_c.mul_by_nonresidue(),
            c2: (self.c0 + self.c2) * (other.c0 + other.c2) - a_a + b_b - c_c,
        }
    }
}
impl MulAssign for Fq12 {
    fn mul_assign(&mut self, rhs: Fq12) {
        *self = *self * rhs;
    }
}

impl Sub for Fq12 {
    type Output = Fq12;

    fn sub(self, other: Fq12) -> Fq12 {
        Fq12 {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
            c2: self.c2 - other.c2,
        }
    }
}

impl Add for Fq12 {
    type Output = Fq12;

    fn add(self, other: Fq12) -> Fq12 {
        Fq12 {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
            c2: self.c2 + other.c2,
        }
    }
}

impl Neg for Fq12 {
    type Output = Fq12;

    fn neg(self) -> Fq12 {
        Fq12 {
            c0: -self.c0,
            c1: -self.c1,
            c2: -self.c2,
        }
    }
}
