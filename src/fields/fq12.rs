#![allow(dead_code)]

use core::ops::{Add, Mul, Neg, Sub};
use rand::Rng;

use crate::fields::{FieldElement, Fq4};
use crate::u256::U256;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq12 {
    pub c0: Fq4,
    pub c1: Fq4,
    pub c2: Fq4,
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

    pub fn scale(&self, by: Fq4) -> Self {
        Fq12 {
            c0: self.c0 * by,
            c1: self.c1 * by,
            c2: self.c2 * by,
        }
    }

    fn final_exponentiation_first_chunk(&self) -> Option<Fq12> {
        match self.inverse() {
            Some(b) => {
                let a = self.frobenius_map(6);
                let c = a * b;
                let d = c.frobenius_map(2);

                Some(d * c)
            }
            None => None,
        }
    }

    fn final_exponentiation_last_chunk(&self) -> Fq12 {
        let a = self.pow(*SM9_A3);
        let b = a.inverse().unwrap();
        let c = b.frobenius_map(1);
        let d = c * b;

        let e = d * b;
        let f = self.frobenius_map(1);
        let g = *self * f;
        let h = g.pow(*SM9_NINE);

        let i = e * h;
        let j = self.squared();
        let k = j.squared();
        let l = k * i;
        let m = f.squared();
        let n = d * m;
        let o = self.frobenius_map(2);
        let p = o * n;

        let q = p.pow(*SM9_A2);
        let r = q * l;
        let s = self.frobenius_map(3);
        let t = s * r;

        t
    }

    pub fn final_exponentiation(&self) -> Option<Fq12> {
        self.final_exponentiation_first_chunk()
            .map(|a| a.final_exponentiation_last_chunk())
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        match power {
            1 => Fq12 {
                c0: self.c0.frobenius_map(10), // case 1 for c0 , so call  10. and so on
                c1: self.c1.frobenius_map(11),
                c2: self.c2.frobenius_map(12),
            },
            2 => {
                Fq12 {
                    c0: self.c0.unitary_inverse(),
                    c1: self.c1.frobenius_map(21), // case 2 for c1 , so call  21
                    c2: self.c2.frobenius_map(22),
                }
            }
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

lazy_static::lazy_static! {

    // a2 = 0xd8000000019062ed 0000b98b0cb27659
    // a3 = 0x2 400000000215d941
    static ref SM9_A2: U256 = U256::from([
            0x0000b98b0cb27659,
            0xd8000000019062ed,
            0x0,
            0x0
        ]);
    static ref SM9_A3: U256 = U256::from([
        0x400000000215d941,
        0x2,
        0x0,
        0x0
    ]);
    static ref SM9_NINE: U256 = U256::from([
        0x9,
        0x0,
        0x0,
        0x0
    ]);
}
