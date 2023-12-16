use alloc::vec::Vec;
use core::fmt;
use core::ops::{Add, Mul, Neg, Sub};
use rand::Rng;

use crate::fields::{FieldElement, FQ, FQ_MINUS1_DIV4, FQ_MINUS5_DIV8};
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
                a.0.mul(&U256::one(), &U256::from($modulus), $inv);

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
                            acc = acc + Self::one();
                            tmp
                        })
                        .collect()
                };

                let mut res = Self::zero();
                for c in s.chars() {
                    match c.to_digit(10) {
                        Some(d) => {
                            res = res * ints[10];
                            res = res + ints[d as usize];
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
                if a < U256::from($modulus) {
                    a.mul(&U256::from($rsquared), &U256::from($modulus), $inv);

                    Some($name(a))
                } else {
                    None
                }
            }

            /// Converts a &[u8] to an Fp so long as it's below the modulus.
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
            /// Converts a U256 to an Fr regardless of modulus.
            pub fn new_mul_factor(mut a: U256) -> Self {
                a.mul(&U256::from($rsquared), &U256::from($modulus), $inv);
                $name(a)
            }

            pub fn interpret(buf: &[u8; 64]) -> Self {
                $name::new(U512::interpret(buf).divrem(&U256::from($modulus)).1).unwrap()
            }

            /// Returns the modulus
            #[inline]
            #[allow(dead_code)]
            pub fn modulus() -> U256 {
                U256::from($modulus)
            }

            pub fn raw(&self) -> &U256 {
                &self.0
            }

            pub fn set_bit(&mut self, bit: usize, to: bool) {
                self.0.set_bit(bit, to);
            }
        }

        impl FieldElement for $name {
            #[inline]
            fn zero() -> Self {
                $name(U256::from([0, 0, 0, 0]))
            }

            #[inline]
            fn one() -> Self {
                $name(U256::from($one))
            }

            fn random<R: Rng>(rng: &mut R) -> Self {
                $name(U256::random(rng, &U256::from($modulus)))
            }

            #[inline]
            fn is_zero(&self) -> bool {
                self.0.is_zero()
            }

            fn inverse(mut self) -> Option<Self> {
                if self.is_zero() {
                    None
                } else {
                    self.0.invert(&U256::from($modulus));
                    self.0
                        .mul(&U256::from($rcubed), &U256::from($modulus), $inv);

                    Some(self)
                }
            }
        }

        impl Add for $name {
            type Output = $name;

            #[inline]
            fn add(mut self, other: $name) -> $name {
                self.0.add(&other.0, &U256::from($modulus));

                self
            }
        }

        impl Sub for $name {
            type Output = $name;

            #[inline]
            fn sub(mut self, other: $name) -> $name {
                self.0.sub(&other.0, &U256::from($modulus));

                self
            }
        }

        impl Mul for $name {
            type Output = $name;

            #[inline]
            fn mul(mut self, other: $name) -> $name {
                self.0.mul(&other.0, &U256::from($modulus), $inv);

                self
            }
        }

        impl Neg for $name {
            type Output = $name;

            #[inline]
            fn neg(mut self) -> $name {
                self.0.neg(&U256::from($modulus));

                self
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
    [
        0xE56EE19CD69ECF25,
        0x49F2934B18EA8BEE,
        0xD603AB4FF58EC744,
        0xB640000002A3A6F1
    ],
    [
        0x7598CD79CD750C35,
        0xE4A08110BB6DAEAB,
        0xBFEE4BAE7D78A1F9,
        0x8894F5D163695D0E
    ],
    [
        0xA8EA85210CE29EF9,
        0x8BD11806993E3A54,
        0x0DB935B5F51A6DA4,
        0x85CB2B73F249E8EC
    ],
    [
        0x1A911E63296130DB,
        0xB60D6CB4E7157411,
        0x29FC54B00A7138BB,
        0x49BFFFFFFD5C590E
    ],
    0xF590740D939A510D1D02662351974B53
);

field_impl!(
    Fq,
    [
        0xE56F9B27E351457D,
        0x21F2934B1A7AEEDB,
        0xD603AB4FF58EC745,
        0xB640000002A3A6F1
    ],
    [
        0x27DEA312B417E2D2,
        0x88F8105FAE1A5D3F,
        0xE479B522D6706E7B,
        0x2EA795A656F62FBD
    ],
    [
        0x130257769DF5827E,
        0x36920FC0837EC76E,
        0xCBEC24519C22A142,
        0x219BE84A7C687090
    ],
    [
        0x1A9064D81CAEBA83,
        0xDE0D6CB4E5851124,
        0x29FC54B00A7138BA,
        0x49BFFFFFFD5C590E
    ],
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
            res = self.pow(*FQ_MINUS5_DIV8) * *self;
        } else if (-a1a).is_one() {
            // (q-5)/8
            // 2g
            let a = *self + *self;
            // 4g
            let aa = a + a;
            // (4g)^u
            let b = aa.pow(*FQ_MINUS5_DIV8);
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

    pub fn tri(&self) -> Self {
        *self + *self + *self
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
        let (_, c0) = u512.divrem(&a);
        Fr::new(c0).map(|f| f + Fr::one())
    }
}
