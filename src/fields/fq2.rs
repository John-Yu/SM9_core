#![allow(dead_code)]

use core::ops::{Add, Mul, Neg, Sub};
use rand::Rng;

use crate::fields::{FieldElement, Fq, FQ};
use crate::u256::{Error, U256};
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
    pub fn scale(&self, by: &Fq) -> Self {
        Fq2 {
            c0: self.c0 * *by,
            c1: self.c1 * *by,
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
    /// Computes the square root of this element, if it exists.
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General
    // Annex C  C.1.4.2
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero());
        }
        let b = self.c1;
        let a = self.c0;
        let bb = b.squared();
        let aa = a.squared();
        let u = aa + bb + bb;
        let mut y = Fq::zero();
        u.sqrt().and_then(|w| {
            let v = (a + w).div2();
            let m = v.sqrt().map(|t| {
                y = t;
            });
            if m.is_none() {
                let v = (a - w).div2();
                v.sqrt().map(|t| {
                    y = t;
                })?;
            }
            let y2 = y + y;
            let z1 = if y.is_zero() {
                // i^2 = -2
                w.div2().sqrt()?
            } else {
                b * y2.inverse().unwrap()
            };
            let z0 = y;
            Some(Self::new(z0, z1))
        })
    }
    /// Converts an element of `Fq2` into a U512
    pub fn to_u512(self) -> U512 {
        let c0: U256 = (*self.real()).into();
        let c1: U256 = (*self.imaginary()).into();

        U512::new(&c1, &c0, &FQ)
    }
    /// Converts an element of `Fq2` into a byte representation in
    /// big-endian byte order.
    pub fn to_slice(self) -> [u8; 64] {
        let mut res = [0u8; 64];
        let b1 = self.c1.to_slice();
        let b0 = self.c0.to_slice();
        res[..32].copy_from_slice(&b1);
        res[32..].copy_from_slice(&b0);
        res
    }
    /// Converts a &[u8] to an element of `Fq2`.
    pub fn from_slice(s: &[u8]) -> Result<Fq2, Error> {
        if s.len() != 64 {
            return Err(Error::InvalidLength {
                expected: 64,
                actual: s.len(),
            });
        }

        let c1 = Fq::from_slice(&s[..32]).unwrap();
        let c0 = Fq::from_slice(&s[32..]).unwrap();

        Ok(Fq2 { c0, c1 })
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
    /// Returns true if element is zero.
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }
    /// Squares this element
    fn squared(&self) -> Self {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Complex squaring), which takes 2M + 3A + 2B
        let a0 = self.c0;
        let a1 = self.c1;
        let v0 = a0 * a1;
        Fq2 {
            c0: (a0 + a1) * (a0 - a1 - a1) + v0,
            c1: v0 + v0,
        }
    }

    fn inverse(self) -> Option<Self> {
        // "High-Speed Software Implementation of the Optimal Ate Pairing
        // over Barreto–Naehrig Curves"; Algorithm 8
        if self.is_zero() {
            None
        } else {
            //  t = (c[0]^2 + 2 * c[1]^2)^-1
            (self.c0.squared() + self.c1.squared() + self.c1.squared())
                .inverse()
                .map(|t| Self::new(self.c0 * t, -(self.c1 * t)))
        }
    }
}

impl Mul for Fq2 {
    type Output = Fq2;

    fn mul(self, other: Fq2) -> Fq2 {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Karatsuba), which costs 3M + 2A + 4B
        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;

        Fq2 {
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
