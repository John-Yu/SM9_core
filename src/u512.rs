#![allow(dead_code)]

use byteorder::{BigEndian, ByteOrder};
use core::cmp::Ordering;
use crunchy::unroll;
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
        let mut res = [0; 4];

        debug_assert_eq!(c1.0.len(), 2);
        unroll! {
            for i in 0..2 {
                mac_digit(i, &mut res, &modulo.0, c1.0[i]);
            }
        }

        let mut carry = 0;

        debug_assert_eq!(res.len(), 4);
        unroll! {
            for i in 0..2 {
                res[i] = adc(res[i], c0.0[i], &mut carry);
            }
        }

        unroll! {
            for i in 0..2 {
                let (a1, a0) = split_u128(res[i + 2]);
                let (c, r0) = split_u128(a0 + carry);
                let (c, r1) = split_u128(a1 + c);
                carry = c;

                res[i + 2] = combine_u128(r1, r0);
            }
        }

        debug_assert!(0 == carry);

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

    pub fn get_bit(&self, n: usize) -> Option<bool> {
        if n >= 512 {
            None
        } else {
            let part = n / 128;
            let bit = n - (128 * part);

            Some(self.0[part] & (1 << bit) > 0)
        }
    }

    /// Divides self by modulo, returning remainder and quotient, if
    /// possible, a quotient smaller than the modulus.
    pub fn divrem(&self, modulo: &U256) -> (Option<U256>, U256) {
        let mut q = Some(U256::zero());
        let mut r = U256::zero();

        for i in (0..512).rev() {
            r.mul2(modulo);
            assert!(r.set_bit(0, self.get_bit(i).unwrap()));
            if &r >= modulo {
                sub_noborrow(&mut r.0, &modulo.0);
                if q.is_some() && !q.as_mut().unwrap().set_bit(i, true) {
                    q = None
                }
            }
        }

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
        for (a, b) in self.0.iter().zip(other.0.iter()).rev() {
            match a.cmp(b) {
                Ordering::Greater => return Ordering::Greater,
                Ordering::Less => return Ordering::Less,
                Ordering::Equal => continue,
            }
        }

        Ordering::Equal
    }
}

impl PartialOrd for U512 {
    #[inline]
    fn partial_cmp(&self, other: &U512) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
