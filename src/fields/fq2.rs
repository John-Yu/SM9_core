#![allow(dead_code)]

use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::TryRngCore;

use crate::{
    One, Zero,
    fields::{FQ, FieldElement, Fq},
    u256::{Error, U256},
    u512::U512,
};

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
#[repr(C)]
pub struct Fq2 {
    pub(crate) c0: Fq,
    pub(crate) c1: Fq,
}

impl Fq2 {
    pub fn new(c0: Fq, c1: Fq) -> Self {
        Fq2 { c0, c1 }
    }

    //Algorithm 7
    pub fn scale(&self, by: &Fq) -> Self {
        Fq2 {
            c0: self.c0 * by,
            c1: self.c1 * by,
        }
    }

    #[inline(always)]
    pub fn unitary_inverse(&self) -> Fq2 {
        Fq2 {
            c0: self.c0,
            c1: -self.c1,
        }
    }

    #[inline(always)]
    pub fn mul_by_nonresidue(&self) -> Self {
        //c0 = -2 * c1
        //c1 = c0
        Fq2 {
            c0: -self.c1.double(),
            c1: self.c0,
        }
    }

    #[inline(always)]
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
        let u = aa + bb.double();
        let mut y = Fq::zero();
        u.sqrt().and_then(|w| {
            let v = (a + w).div2();
            let m = v.sqrt().map(|t| {
                y = t;
            });
            if m.is_none() {
                let v = (a - w).div2();
                y = v.sqrt()?;
            }
            let y2 = y.double();
            let z1 = if y.is_zero() {
                // i^2 = -2
                w.div2().sqrt()?
            } else {
                b * y2.inverse()?
            };
            let z0 = y;
            let sqrt_cand = Self::new(z0, z1);
            // Check if sqrt_cand is actually the square root
            // if not, there exists no square root.
            if sqrt_cand.squared() == *self {
                Some(sqrt_cand)
            } else {
                None
            }
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
    #[inline]
    pub fn neg_inplace(&self) -> Fq2 {
        Fq2 {
            c0: self.c0.neg_inplace(),
            c1: self.c1.neg_inplace(),
        }
    }
    #[inline]
    pub fn sub_inplace(&self, rhs: &Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0.sub_inplace(&rhs.c0),
            c1: self.c1.sub_inplace(&rhs.c1),
        }
    }

    #[inline]
    pub fn add_inplace(&self, rhs: &Fq2) -> Fq2 {
        Fq2 {
            c0: self.c0.add_inplace(&rhs.c0),
            c1: self.c1.add_inplace(&rhs.c1),
        }
    }
    /// Returns `c = self * b`.
    ///
    /// Implements the full-tower interleaving strategy from
    /// [ePrint 2022-376](https://eprint.iacr.org/2022/367).
    #[inline]
    pub fn mul_inplace(&self, b: &Fq2) -> Fq2 {
        // F_{p^2} x F_{p^2} multiplication implemented with operand scanning (schoolbook)
        // computes the result as:
        //
        //   a·b = (a_0 b_0 + a_1 b_1 β) + (a_0 b_1 + a_1 b_0)i
        //
        // In SM9's F_{p^2}, our β is -2, so the resulting F_{p^2} element is:
        //
        //   c_0 = a_0 b_0 - 2 a_1 b_1
        //   c_1 = a_0 b_1 + a_1 b_0
        //
        // Each of these is a "sum of products", which we can compute efficiently.

        let a = self;
        Fq2 {
            c0: Fq::sum_of_products(&[a.c0, -a.c1.double()], &[b.c0, b.c1]),
            c1: Fq::sum_of_products(&[a.c0, a.c1], &[b.c1, b.c0]),
        }
    }
}

impl_binops_additive!(Fq2, Fq2);
impl_binops_multiplicative!(Fq2, Fq2);
impl_binops_negative!(Fq2);

impl Zero for Fq2 {
    #[inline]
    fn zero() -> Self {
        Fq2 {
            c0: Fq::zero(),
            c1: Fq::zero(),
        }
    }
    /// Returns true if element is zero.
    #[inline(always)]
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }
}
impl One for Fq2 {
    #[inline]
    fn one() -> Self {
        Fq2 {
            c0: Fq::one(),
            c1: Fq::zero(),
        }
    }
}
impl FieldElement for Fq2 {
    /// Returns an element chosen uniformly at random using a user-provided RNG.
    fn try_from_rng<R: TryRngCore + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
        let c0 = Fq::try_from_rng(rng)?;
        let c1 = Fq::try_from_rng(rng)?;
        Ok(Fq2 { c0, c1 })
    }
    /// double this element
    #[inline(always)]
    fn double(&self) -> Self {
        Fq2 {
            c0: self.c0.double(),
            c1: self.c1.double(),
        }
    }
    /// triple this element
    #[inline(always)]
    fn triple(&self) -> Self {
        Fq2 {
            c0: self.c0.triple(),
            c1: self.c1.triple(),
        }
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
            c0: (a0 + a1) * (a0 - a1.double()) + v0,
            c1: v0.double(),
        }
    }
    /// Computes the multiplicative inverse of this field element
    /// return None in the case that this element is zero
    fn inverse(&self) -> Option<Self> {
        // "High-Speed Software Implementation of the Optimal Ate Pairing
        // over Barreto–Naehrig Curves"; Algorithm 8
        //  t = (c[0]^2 + 2 * c[1]^2)^-1
        (self.c0.squared() + self.c1.squared().double())
            .inverse()
            .map(|t| Self::new(self.c0 * t, -(self.c1 * t)))
    }
}
