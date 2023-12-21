use alloc::vec::Vec;
use core::fmt;
use core::ops::{Add, AddAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::Rng;

use crate::fields::{
    FieldElement, FQ, FQ_CUBED, FQ_MINUS1_DIV4, FQ_MINUS5_DIV8, FQ_ONE, FQ_SQUARED, FR, FR_CUBED,
    FR_ONE, FR_SQUARED,
};
use crate::u256::U256;
use crate::u512::U512;

macro_rules! field_impl {
    ($name:ident, $modulus:expr, $rsquared:expr, $rcubed:expr, $one:expr, $inv:expr) => {
        #[derive(Copy, Clone, PartialEq, Eq)]
        #[repr(C)]
        pub struct $name(U256);

        impl From<$name> for U256 {
            #[inline]
            fn from(mut a: $name) -> Self {
                a.0.mul(&U256::one(), &$modulus, $inv);

                a.0
            }
        }

        impl $name {
            pub fn from_str(s: &str) -> Option<Self> {
                let ints: Vec<_> = {
                    let mut acc = Self::zero();
                    (0..11)
                        .map(|_| {
                            let tmp = acc;
                            acc += Self::one();
                            tmp
                        })
                        .collect()
                };

                let mut res = Self::zero();
                for c in s.chars() {
                    match c.to_digit(10) {
                        Some(d) => {
                            res *= ints[10];
                            res += ints[d as usize];
                        }
                        None => {
                            return None;
                        }
                    }
                }

                Some(res)
            }

            /// Converts a U256 to an Fp so long as it's below the modulus.
            pub fn new(mut a: U256) -> Option<Self> {
                if a < *$modulus {
                    a.mul(&$rsquared, &$modulus, $inv);

                    Some($name(a))
                } else {
                    None
                }
            }

            /// Converts a &[u8] to an Fp so long as it's below the modulus.
            /// the length of hex must be 32
            pub fn from_slice(hex: &[u8]) -> Option<Self> {
                match U256::from_slice(hex) {
                    Ok(a) => Self::new(a),
                    Err(_) => None,
                }
            }

            /// Converts an element of `Fp` into a byte representation in big-endian byte order.
            pub fn to_slice(self) -> [u8; 32] {
                let mut res = [0u8; 32];
                U256::from(self).to_big_endian(&mut res[..]).unwrap();
                res
            }
            /// Converts a U256 to a Fp regardless of modulus.
            pub fn new_mul_factor(mut a: U256) -> Self {
                a.mul(&$rsquared, &$modulus, $inv);
                $name(a)
            }
            /// Converts a &[u8; 64] to a Fp regardless of modulus.
            pub fn interpret(buf: &[u8; 64]) -> Self {
                $name::new(U512::interpret(buf).divrem(&$modulus).1).unwrap()
            }

            /// Returns the modulus
            #[inline]
            #[allow(dead_code)]
            pub fn modulus() -> U256 {
                *$modulus
            }

            pub fn raw(&self) -> &U256 {
                &self.0
            }

            pub fn set_bit(&mut self, bit: usize, to: bool) {
                self.0.set_bit(bit, to);
            }

            #[inline]
            pub fn add_inplace(&self, other: &$name) -> $name {
                let mut a = self.0;
                a.add(&other.0, &$modulus);

                $name(a)
            }
            #[inline]
            pub fn sub_inplace(&self, other: &$name) -> $name {
                let mut a = self.0;
                a.sub(&other.0, &$modulus);

                $name(a)
            }
            #[inline]
            pub fn mul_inplace(&self, other: &$name) -> $name {
                let mut a = self.0;
                a.mul(&other.0, &$modulus, $inv);

                $name(a)
            }
            #[inline]
            pub fn neg_inplace(&self) -> $name {
                let mut a = self.0;
                a.neg(&$modulus);

                $name(a)
            }
        }

        impl_binops_additive!($name, $name);
        impl_binops_multiplicative!($name, $name);
        impl_binops_negative!($name);

        impl FieldElement for $name {
            #[inline]
            fn zero() -> Self {
                $name(U256::zero())
            }

            #[inline]
            fn one() -> Self {
                $name(*$one)
            }

            fn random<R: Rng>(rng: &mut R) -> Self {
                $name(U256::random(rng, &$modulus))
            }

            #[inline]
            fn is_zero(&self) -> bool {
                self.0.is_zero()
            }
            /// Computes the inverse of this element,
            /// None if the element is zero.
            fn inverse(&self) -> Option<Self> {
                if self.is_zero() {
                    None
                } else {
                    let mut a = self.0;
                    a.invert(&$modulus);
                    a.mul(&$rcubed, &$modulus, $inv);

                    Some(Self(a))
                }
            }
            /// double this element
            #[inline]
            fn double(&self) -> Self {
                let mut a = self.0;
                a.mul2(&$modulus);
                if &a >= &$modulus {
                    a.sub(&$modulus, &$modulus);
                }
                $name(a)
            }
            /// triple this element
            fn triple(&self) -> Self {
                &self.double() + self
            }
            /// Squares this element.
            fn squared(&self) -> Self {
                self * self
            }
        }

        impl From<$name> for [u8; 32] {
            fn from(value: $name) -> [u8; 32] {
                value.to_slice()
            }
        }

        impl<'a> From<&'a $name> for [u8; 32] {
            fn from(value: &'a $name) -> [u8; 32] {
                value.to_slice()
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "Fp({:?})", U256::from(*self))
            }
        }
    };
}

field_impl!(
    Fr,
    FR,
    FR_SQUARED,
    FR_CUBED,
    FR_ONE,
    0xF590740D939A510D1D02662351974B53
);

field_impl!(
    Fq,
    FQ,
    FQ_SQUARED,
    FQ_CUBED,
    FQ_ONE,
    0x181AE39613C8DBAF892BC42C2F2EE42B
);

impl Fq {
    /// Computes the square root of this element, if it exists.
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General
    // Annex C  C.1.4.1, Algorithm 2
    pub fn sqrt(&self) -> Option<Self> {
        if self.is_zero() {
            return Some(Self::zero());
        }
        // 2u+1 = (q-1)/4
        let a1a = self.pow(*FQ_MINUS1_DIV4);
        let mut res = Fq::zero();
        if a1a.is_one() {
            // g^(u+1) = g^u * g
            res = self.pow(*FQ_MINUS5_DIV8) * self;
        } else if (-a1a).is_one() {
            // (q-5)/8
            // 2g
            let a = self.double();
            // (4g)^u
            let b = a.double().pow(*FQ_MINUS5_DIV8);
            res = a * b;
        }
        if res.is_zero() {
            None
        } else {
            // return positive number
            let r = -res;
            if U256::from(r) < U256::from(res) {
                res = r;
            }
            Some(res)
        }
    }

    pub fn div2(mut self) -> Self {
        self.0.div2(&FQ);
        self
    }
    /// Returns true if element is one.
    pub fn is_one(&self) -> bool {
        *self == Self::one()
    }
}

impl Fr {
    /// for H1() and H2()
    // h = (Ha mod (n-1)) + 1;  h in [1, n-1], n is the curve order, Ha is 40 bytes from hash
    pub fn from_hash(ha: &[u8]) -> Option<Self> {
        if ha.len() > 64 {
            return None;
        }
        let mut v = [0u8; 64];
        let start = 64 - ha.len();
        v[start..].copy_from_slice(ha);
        let u512 = U512::interpret(&v);
        // n-1
        let a = U256::from(-Fr::one());
        // h = (Ha mod (n-1)) + 1
        Fr::new(u512.divrem(&a).1).map(|f| f + Fr::one())
    }
}
