#![allow(dead_code)]

use byteorder::{BigEndian, ByteOrder};
use core::cmp::Ordering;
use rand::Rng;

use crate::arith::*;
use crate::u256::{Error, U256};

/// 512-bit, stack allocated biginteger for use in extension
/// field serialization and scalar interpretation.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct U512(pub(crate) [u128; 4]);

impl From<[u64; 8]> for U512 {
    fn from(d: [u64; 8]) -> Self {
        let mut a = [0u128; 4];
        a[0] = (d[1] as u128) << 64 | d[0] as u128;
        a[1] = (d[3] as u128) << 64 | d[2] as u128;
        a[2] = (d[5] as u128) << 64 | d[4] as u128;
        a[3] = (d[7] as u128) << 64 | d[6] as u128;
        U512(a)
    }
}

impl U512 {
    /// Multiplies c1 by modulo, adds c0.
    pub fn new(c1: &U256, c0: &U256, modulo: &U256) -> U512 {
        let mut res = mul(&c1.0, &modulo.0);
        let mut carry = false;
        (res[0], carry) = carrying_add(res[0], c0.0[0], carry);
        (res[1], carry) = carrying_add(res[1], c0.0[1], carry);
        (res[2], carry) = carrying_add(res[2], 0, carry);
        (res[3], carry) = carrying_add(res[3], 0, carry);
        debug_assert!(!carry);

        U512(res)
    }

    #[inline]
    pub fn one() -> U512 {
        U512([1, 0, 0, 0])
    }

    pub fn from_slice(s: &[u8]) -> Result<U512, Error> {
        if s.len() != 64 {
            return Err(Error::InvalidLength {
                expected: 64,
                actual: s.len(),
            });
        }

        let mut n = [0; 4];
        for (l, i) in (0..4).rev().zip((0..4).map(|i| i * 16)) {
            n[l] = BigEndian::read_u128(&s[i..]);
        }

        Ok(U512(n))
    }

    /// Get a random U512
    pub fn random<R: Rng>(rng: &mut R) -> U512 {
        U512(rng.gen())
    }
    /// Returns the number of bits necessary to represent this value. Note that zero
    /// is considered to need 0 bits.
    pub fn bit_length(&self) -> usize {
        let mut ret = 512;
        for i in self.0.iter().rev() {
            let leading = i.leading_zeros();
            ret -= leading;
            if leading != u128::BITS {
                break;
            }
        }

        ret as usize
    }
    /// Returns the `n`-th bit where bit 0 is the least significant one. true if bit is 1.
    /// In other words, the bit with weight `2^n`.
    #[inline]
    pub fn get_bit(&self, n: usize) -> Option<bool> {
        if n >= 512 {
            None
        } else {
            let part = n >> 7;
            let bit = n & 0x7F;

            Some(self.0[part] & (1 << bit) > 0)
        }
    }

    /// Divides self by modulo, returning remainder and quotient, if
    /// possible, a quotient smaller than the modulus.
    // Stupid slow base-2 long division taken from
    // https://en.wikipedia.org/wiki/Division_algorithm
    pub fn divrem(&self, modulo: &U256) -> (Option<U256>, U256) {
        let mut q = Some(U256::zero());
        let mut r = U256::zero();
        let bits = self.bit_length();

        for i in (0..bits).rev() {
            let carry = mul2(&mut r.0);
            assert!(r.set_bit(0, self.get_bit(i).unwrap()));
            if &r >= modulo || carry {
                if carry {
                    r.add_carry(modulo);
                }
                if &r >= modulo {
                    sub_noborrow(&mut r.0, &modulo.0);
                }
                if q.is_some() && !q.as_mut().unwrap().set_bit(i, true) {
                    q = None;
                }
            }
        }
        debug_assert!(q.is_none() || *self == Self::new(&q.unwrap(), &r, modulo));
        if q.is_some() && (q.as_ref().unwrap() >= modulo) {
            (None, r)
        } else {
            (q, r)
        }
    }

    pub fn interpret(buf: &[u8; 64]) -> U512 {
        let mut n = [0; 4];
        for (l, i) in (0..4).rev().zip((0..4).map(|i| i * 16)) {
            n[l] = BigEndian::read_u128(&buf[i..]);
        }

        U512(n)
    }
}

impl Ord for U512 {
    #[inline]
    fn cmp(&self, other: &U512) -> Ordering {
        let lhs = self.0.iter().cloned().rev();
        let rhs = other.0.iter().cloned().rev();
        lhs.cmp(rhs)
    }
}

impl PartialOrd for U512 {
    #[inline]
    fn partial_cmp(&self, other: &U512) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
