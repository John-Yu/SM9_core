#![allow(dead_code)]

use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
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
    /// Multiply by quadratic nonresidue.
    #[inline]
    pub fn mul_by_nonresidue(&self) -> Self {
        Fq12 {
            c0: self.c2.mul_by_nonresidue(),
            c1: self.c0,
            c2: self.c1,
        }
    }

    pub fn scale(&self, by: &Fq4) -> Self {
        Fq12 {
            c0: self.c0 * by,
            c1: self.c1 * by,
            c2: self.c2 * by,
        }
    }

    #[inline]
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
    pub fn to_slice(self) -> [u8; 384] {
        let mut res = [0u8; 384];
        let b2 = self.c2.to_slice();
        let b1 = self.c1.to_slice();
        let b0 = self.c0.to_slice();
        res[..128].copy_from_slice(&b2);
        res[128..256].copy_from_slice(&b1);
        res[256..].copy_from_slice(&b0);
        res
    }

    #[inline]
    fn neg_inplace(&self) -> Fq12 {
        Fq12 {
            c0: self.c0.neg_inplace(),
            c1: self.c1.neg_inplace(),
            c2: self.c2.neg_inplace(),
        }
    }
    #[inline]
    fn add_inplace(&self, rhs: &Fq12) -> Fq12 {
        Fq12 {
            c0: self.c0.add_inplace(&rhs.c0),
            c1: self.c1.add_inplace(&rhs.c1),
            c2: self.c2.add_inplace(&rhs.c2),
        }
    }
    #[inline]
    fn sub_inplace(&self, rhs: &Fq12) -> Fq12 {
        Fq12 {
            c0: self.c0.sub_inplace(&rhs.c0),
            c1: self.c1.sub_inplace(&rhs.c1),
            c2: self.c2.sub_inplace(&rhs.c2),
        }
    }
    // Multiplication and Squaring on Pairing-Friendly Fields.pdf
    // Section 4 (Karatsuba)
    #[inline]
    fn mul_inplace(&self, other: &Fq12) -> Fq12 {
        let a_a = self.c0.mul_inplace(&other.c0);
        let b_b = self.c1.mul_inplace(&other.c1);
        let c_c = self.c2.mul_inplace(&other.c2);

        Fq12 {
            c0: ((self.c1 + self.c2) * (other.c1 + other.c2) - b_b - c_c).mul_by_nonresidue() + a_a,
            c1: (self.c0 + self.c1) * (other.c0 + other.c1) - a_a - b_b + c_c.mul_by_nonresidue(),
            c2: (self.c0 + self.c2) * (other.c0 + other.c2) - a_a + b_b - c_c,
        }
    }
}

impl_binops_additive!(Fq12, Fq12);
impl_binops_multiplicative!(Fq12, Fq12);
impl_binops_negative!(Fq12);

impl FieldElement for Fq12 {
    #[inline]
    fn zero() -> Self {
        Fq12 {
            c0: Fq4::zero(),
            c1: Fq4::zero(),
            c2: Fq4::zero(),
        }
    }

    #[inline]
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

    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero() && self.c2.is_zero()
    }
    /// double this element
    #[inline(always)]
    fn double(&self) -> Self {
        Fq12 {
            c0: self.c0.double(),
            c1: self.c1.double(),
            c2: self.c2.double(),
        }
    }
    /// triple this element
    #[inline(always)]
    fn triple(&self) -> Self {
        Fq12 {
            c0: self.c0.triple(),
            c1: self.c1.triple(),
            c2: self.c2.triple(),
        }
    }
    //Multiplication and Squaring on Pairing-Friendly Fields.pdf
    // Section 4 (CH-SQR2)
    #[inline]
    fn squared(&self) -> Self {
        let s0 = self.c0.squared();
        let s1 = (self.c0 * self.c1).double();
        let s2 = (self.c0 - self.c1 + self.c2).squared();
        let s3 = (self.c1 * self.c2).double();
        let s4 = self.c2.squared();

        Fq12 {
            c0: s0 + s3.mul_by_nonresidue(),
            c1: s1 + s4.mul_by_nonresidue(),
            c2: s1 + s2 + s3 - s0 - s4,
        }
    }
    /// Computes the multiplicative inverse of this field element
    /// return None in the case that this element is zero
    //"High-Speed Software Implementation of the Optimal Ate AbstractPairing over Barreto-Naehrig Curves"
    // Algorithm 17
    #[inline]
    fn inverse(&self) -> Option<Self> {
        let c0 = self.c0.squared() - self.c1 * self.c2.mul_by_nonresidue();
        let c1 = self.c2.squared().mul_by_nonresidue() - self.c0 * self.c1;
        let c2 = self.c1.squared() - self.c0 * self.c2;

        ((self.c2 * c1 + self.c1 * c2).mul_by_nonresidue() + self.c0 * c0)
            .inverse()
            .map(|t| Self::new(t * c0, t * c1, t * c2))
    }
}
