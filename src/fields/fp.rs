use alloc::vec::Vec;
use core::fmt;
use core::ops::{Add, AddAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign};
use rand::TryRngCore;

use crate::{
    One, Zero,
    arith::{adc, mac},
    fields::{
        FQ, FQ_INV, FQ_MINUS1_DIV4, FQ_MINUS5_DIV8, FQ_ONE, FQ_SQUARED, FR, FR_INV, FR_ONE,
        FR_SQUARED, FieldElement,
    },
    u256::U256,
    u512::U512,
};

macro_rules! field_impl {
    ($name:ident, $modulus:expr, $rsquared:expr, $one:expr, $inv:expr) => {
        #[derive(Copy, Clone, PartialEq, Eq)]
        #[repr(C)]
        pub struct $name(pub(crate) U256);

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
                    if !a.is_zero() {
                        a.mul(&$rsquared, &$modulus, $inv);
                    }
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

        impl Zero for $name {
            #[inline]
            fn zero() -> Self {
                $name(U256::zero())
            }

            #[inline]
            fn is_zero(&self) -> bool {
                self.0.is_zero()
            }
        }
        impl One for $name {
            #[inline]
            fn one() -> Self {
                $name(*$one)
            }
        }
        impl FieldElement for $name {
            /// Returns an element chosen uniformly at random using a user-provided RNG.
            fn try_from_rng<R: TryRngCore + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
                Ok($name(U256::try_from_rng(rng, &$modulus)?))
            }

            /// Computes the inverse of this element,
            /// None if the element is zero.
            fn inverse(&self) -> Option<Self> {
                if self.is_zero() {
                    None
                } else {
                    let mut a = self.0;
                    a.invert(&$modulus, &$rsquared);
                    Some($name(a))
                }
            }
            /// double this element
            #[inline]
            fn double(&self) -> Self {
                let mut a = self.0;
                a.mul2(&$modulus);
                $name(a)
            }
            /// triple this element
            fn triple(&self) -> Self {
                &self.double() + self
            }
            /// Squares this element.
            fn squared(&self) -> Self {
                let mut a = self.0;
                a.square(&$modulus, $inv);
                $name(a)
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

        impl Index<usize> for $name {
            type Output = u64;
            #[inline(always)]
            fn index(&self, index: usize) -> &Self::Output {
                &self.0[index]
            }
        }

        impl fmt::Debug for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "Fp({:?})", U256::from(*self))
            }
        }
    };
}

field_impl!(Fr, FR, FR_SQUARED, FR_ONE, *FR_INV);

field_impl!(Fq, FQ, FQ_SQUARED, FQ_ONE, *FQ_INV);

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
    /// Returns `c = a.zip(b).fold(0, |acc, (a_i, b_i)| acc + a_i * b_i)`.
    ///
    /// Implements Algorithm 2 from Patrick Longa's
    /// [ePrint 2022-367](https://eprint.iacr.org/2022/367) $3.
    #[inline]
    pub(crate) fn sum_of_products<const T: usize>(a: &[Fq; T], b: &[Fq; T]) -> Fq {
        // For a single `a x b` multiplication, operand scanning (schoolbook) takes each
        // limb of `a` in turn, and multiplies it by all of the limbs of `b` to compute
        // the result as a double-width intermediate representation, which is then fully
        // reduced at the end. Here however we have pairs of multiplications (a_i, b_i),
        // the results of which are summed.
        //
        // The intuition for this algorithm is two-fold:
        // - We can interleave the operand scanning for each pair, by processing the jth
        //   limb of each `a_i` together. As these have the same offset within the overall
        //   operand scanning flow, their results can be summed directly.
        // - We can interleave the multiplication and reduction steps, resulting in a
        //   single bitshift by the limb size after each iteration. This means we only
        //   need to store a single extra limb overall, instead of keeping around all the
        //   intermediate results and eventually having twice as many limbs.

        // Algorithm 2, line 2
        let (u0, u1, u2, u3, u4) = (0..4).fold((0, 0, 0, 0, 0), |(u0, u1, u2, u3, u4), j| {
            // Algorithm 2, line 3
            // For each pair in the overall sum of products:
            let (t0, t1, t2, t3, t4, t5) =
                (0..T).fold((u0, u1, u2, u3, u4, 0), |(t0, t1, t2, t3, t4, t5), i| {
                    let d = a[i].0[j];
                    let e = b[i].0.as_ref();
                    // Compute digit_j x row and accumulate into `u`.
                    let (t0, carry) = mac(t0, d, e[0], 0);
                    let (t1, carry) = mac(t1, d, e[1], carry);
                    let (t2, carry) = mac(t2, d, e[2], carry);
                    let (t3, carry) = mac(t3, d, e[3], carry);
                    let (t4, carry) = adc(t4, 0, carry);
                    let (t5, _) = adc(t5, 0, carry);

                    (t0, t1, t2, t3, t4, t5)
                });

            // Algorithm 2, lines 4-5
            // This is a single step of the usual Montgomery reduction process.
            // 4：q ← u · p′ mod 2^64,
            let k = t0.wrapping_mul(*FQ_INV);
            // 5：u ← (u + q · p)/2^64，
            let m = FQ.as_ref();
            let (_, carry) = mac(t0, k, m[0], 0);
            let (r1, carry) = mac(t1, k, m[1], carry);
            let (r2, carry) = mac(t2, k, m[2], carry);
            let (r3, carry) = mac(t3, k, m[3], carry);
            let (r4, carry) = adc(t4, 0, carry);
            let (r5, _) = adc(t5, 0, carry);
            (r1, r2, r3, r4, r5)
        });
        // Because we represent F_p elements in non-redundant form, we need a final
        // conditional subtraction to ensure the output is in range.
        let mut r = U256::from([u0, u1, u2, u3]);
        if u4 != 0 {
            // has carry
            for _ in 0..u4 {
                r.add_carry(&FQ);
            }
        }
        // to ensure the output is in range.
        r.subtract_modulus_with_carry(&FQ, false);
        Fq(r)
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
