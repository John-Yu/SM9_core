#![allow(dead_code)]

use core::ops::{Add, Mul, Neg, Sub};
use rand::Rng;

use crate::fields::{FieldElement, Fq};
use crate::u256::U256;
use crate::u512::U512;

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
#[repr(C)]
pub struct Fq2 {
    c0: Fq,
    c1: Fq,
}

impl Fq2 {
    pub fn new(c0: Fq, c1: Fq) -> Self {
        Fq2 { c0, c1 }
    }

    //Algorithm 7
    pub fn scale(&self, by: Fq) -> Self {
        Fq2 {
            c0: self.c0 * by,
            c1: self.c1 * by,
        }
    }

    pub fn unitary_inverse(&self) -> Fq2 {
        Fq2::new(self.c0, -self.c1)
    }

    pub fn mul_by_nonresidue(&self) -> Self {
        //c0 = -2 * c1
        //c1 = c0
        Fq2 {
            c0: -(self.c1 + self.c1),
            c1: self.c0,
        }
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        if power % 2 == 0 {
            *self
        } else {
            Fq2 {
                c0: self.c0,
                c1: -self.c1,
            }
        }
    }

    // c^2 * beta , where beta^2 = -2
    pub fn squared_u(&self) -> Self {
        let ab = self.c0 * self.c1;
        let aa = self.c0 * self.c0;
        let bb = self.c1 * self.c1;

        //r0 = -4 * c0 * c1
        //r1 = c0^2 - 2 * c1^2
        Fq2 {
            c0: -(ab + ab + ab + ab),
            c1: aa - bb - bb,
        }
    }

    pub fn tri(&self) -> Self {
        Fq2 {
            c0: self.c0.tri(),
            c1: self.c1.tri(),
        }
    }

    pub fn div2(&self) -> Self {
        Fq2 {
            c0: self.c0.div2(),
            c1: self.c1.div2(),
        }
    }

    pub fn real(&self) -> &Fq {
        &self.c0
    }

    pub fn imaginary(&self) -> &Fq {
        &self.c1
    }

    pub fn i() -> Self {
        Fq2::new(Fq::zero(), Fq::one())
    }

    // self * other * i
    pub fn mul_u(&self, other: Fq2) -> Self {
        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;
        let ab = (self.c0 + self.c1) * (other.c0 + other.c1) - aa - bb;

        Fq2 {
            // c0 =  -2 * (a0 * b1 + a1 * b0)
            // c1 = aa - 2*bb
            c0: -(ab + ab),
            c1: aa - bb - bb,
        }
    }

    pub fn sqrt(&self) -> Option<Self> {
        let a1 = self.pow::<U256>((*FQ_MINUS3_DIV4).into());
        let a1a = a1 * *self;
        let alpha = a1 * a1a;
        let a0 = alpha.pow(*FQ) * alpha;

        if a0 == Fq2::one().neg() {
            return None;
        }

        if alpha == Fq2::one().neg() {
            Some(Self::i() * a1a)
        } else {
            let b = (alpha + Fq2::one()).pow::<U256>((*FQ_MINUS1_DIV2).into());
            Some(b * a1a)
        }
    }

    pub fn to_u512(&self) -> U512 {
        let c0: U256 = (*self.real()).into();
        let c1: U256 = (*self.imaginary()).into();

        U512::new(&c1, &c0, &FQ)
    }
}

impl FieldElement for Fq2 {
    fn zero() -> Self {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        }
    }

    fn one() -> Self {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        }
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        Fq2 {
            c0: Fq::random(rng),
            c1: Fq::random(rng),
        }
    }

    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    fn squared(&self) -> Self {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Complex squaring)

        let ab = self.c0 * self.c1;
        let aa = self.c0 * self.c0;
        let bb = self.c1 * self.c1;

        //c0: c0^2 - 2 * c1^2
        Fq2 {
            // c0: (self.c1 * fq_non_residue() + self.c0) * (self.c0 + self.c1) - ab - ab * fq_non_residue(),
            c0: aa - bb - bb,
            c1: ab + ab,
        }
    }

    fn inverse(self) -> Option<Self> {
        // "High-Speed Software Implementation of the Optimal Ate Pairing
        // over Barreto–Naehrig Curves"; Algorithm 8

        //  t = (c[0]^2 + 2 * c[1]^2)^-1
        match (self.c0.squared() + self.c1.squared() + self.c1.squared()).inverse() {
            Some(t) => Some(Fq2 {
                c0: self.c0 * t,
                c1: -(self.c1 * t),
            }),
            None => None,
        }
    }
}

impl Mul for Fq2 {
    type Output = Fq2;

    fn mul(self, other: Fq2) -> Fq2 {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Karatsuba)

        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;

        Fq2 {
            // beta = -2 ,so c0 = aa - 2*bb
            // c0: bb * fq_non_residue() + aa,
            c0: aa - bb - bb,
            c1: (self.c0 + self.c1) * (other.c0 + other.c1) - aa - bb,
        }
    }
}

impl Sub for Fq2 {
    type Output = Fq2;

    fn sub(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
        }
    }
}

impl Add for Fq2 {
    type Output = Fq2;

    fn add(self, other: Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }
}

impl Neg for Fq2 {
    type Output = Fq2;

    fn neg(self) -> Fq2 {
        Fq2 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

lazy_static::lazy_static! {
    static ref FQ: U256 = U256::from([
        0xE56F9B27E351457D,
        0x21F2934B1A7AEEDB,
        0xD603AB4FF58EC745,
        0xB640000002A3A6F1
    ]);

    static ref FQ_MINUS3_DIV4: Fq =
        Fq::new(3.into()).expect("3 is a valid field element and static; qed").neg() *
        Fq::new(4.into()).expect("4 is a valid field element and static; qed").inverse()
        .expect("4 has inverse in Fq and is static; qed");

    static ref FQ_MINUS1_DIV2: Fq =
        Fq::new(1.into()).expect("1 is a valid field element and static; qed").neg() *
        Fq::new(2.into()).expect("2 is a valid field element and static; qed").inverse()
            .expect("2 has inverse in Fq and is static; qed");
}