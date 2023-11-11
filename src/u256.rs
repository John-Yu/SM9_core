use byteorder::{BigEndian, ByteOrder};
use core::cmp::Ordering;
use core::fmt;
use rand::Rng;

use crate::arith::*;
use crate::u512::U512;

/// 256-bit, stack allocated biginteger for use in prime field
/// arithmetic.
#[derive(Copy, Clone, PartialEq, Eq)]
#[repr(C)]
pub struct U256(pub [u128; 2]);

impl From<[u64; 4]> for U256 {
    fn from(d: [u64; 4]) -> Self {
        let mut a = [0u128; 2];
        a[0] = (d[1] as u128) << 64 | d[0] as u128;
        a[1] = (d[3] as u128) << 64 | d[2] as u128;
        U256(a)
    }
}

impl From<u64> for U256 {
    fn from(d: u64) -> Self {
        U256::from([d, 0, 0, 0])
    }
}

impl U256 {
    /// Initialize U256 from slice of bytes (big endian)
    pub fn from_slice(s: &[u8]) -> Result<U256, Error> {
        if s.len() != 32 {
            return Err(Error::InvalidLength {
                expected: 32,
                actual: s.len(),
            });
        }

        let mut n = [0; 2];
        for (l, i) in (0..2).rev().zip((0..2).map(|i| i * 16)) {
            n[l] = BigEndian::read_u128(&s[i..]);
        }

        Ok(U256(n))
    }
    /// Converts a U256 into a slice of bytes (big endian)
    pub fn to_big_endian(&self, s: &mut [u8]) -> Result<(), Error> {
        if s.len() != 32 {
            return Err(Error::InvalidLength {
                expected: 32,
                actual: s.len(),
            });
        }

        for (l, i) in (0..2).rev().zip((0..2).map(|i| i * 16)) {
            BigEndian::write_u128(&mut s[i..], self.0[l]);
        }

        Ok(())
    }

    #[inline]
    pub fn zero() -> U256 {
        U256([0, 0])
    }

    #[inline]
    pub fn one() -> U256 {
        U256([1, 0])
    }

    /// Produce a random number (mod `modulo`)
    pub fn random<R: Rng>(rng: &mut R, modulo: &U256) -> U256 {
        U512::random(rng).divrem(modulo).1
    }

    pub fn is_zero(&self) -> bool {
        self.0[0] == 0 && self.0[1] == 0
    }

    pub fn set_bit(&mut self, n: usize, to: bool) -> bool {
        if n >= 256 {
            false
        } else {
            let part = n / 128;
            let bit = n - (128 * part);

            if to {
                self.0[part] |= 1 << bit;
            } else {
                self.0[part] &= !(1 << bit);
            }

            true
        }
    }

    pub fn get_bit(&self, n: usize) -> Option<bool> {
        if n >= 256 {
            None
        } else {
            let part = n / 128;
            let bit = n - (128 * part);

            Some(self.0[part] & (1 << bit) > 0)
        }
    }

    // self = self + 2^256 (mod `modulo`)
    fn add_carry(&mut self, modulo: &U256) {
        let mut a = U256::from([
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF,
            0xFFFFFFFFFFFFFFFF
        ]);
        // TODO: need  self + 2^256 < modulo * 2 
        sub_noborrow(&mut a.0, &modulo.0);
        add_nocarry(&mut self.0, &U256::one().0); //&[1,0]); 
        add_nocarry(&mut self.0, &a.0);   
    }

    /// Add `other` to `self` (mod `modulo`)
    pub fn add(&mut self, other: &U256, modulo: &U256) {
        if add_carry(&mut self.0, &other.0) {
            // has carry
            self.add_carry(modulo);
        }

        if *self >= *modulo {
            sub_noborrow(&mut self.0, &modulo.0);
        }
    }

    /// Subtract `other` from `self` (mod `modulo`)
    pub fn sub(&mut self, other: &U256, modulo: &U256) {
        if *self < *other {
            //(a + q) - b = a + (q - b)
            let mut a = *modulo;
            sub_noborrow(&mut a.0, &other.0);
            add_nocarry(&mut self.0, &a.0);
        }
        else {
            sub_noborrow(&mut self.0, &other.0);
        }
    }

    // a = a * 2 (mod `modulo`)
    pub fn mul2(&mut self, modulo: &U256) {
        if mul2(&mut self.0) {
            // has carry
            self.add_carry(modulo);
        }
    }

    // a = a / 2 (mod `modulo`)
    pub fn div2(&mut self, modulo: &U256) {
        let mut carry = false;
        if  self.is_even() == false {
            carry = add_carry(&mut self.0, &modulo.0);
        }
        div2(&mut self.0);
        if carry {
            self.set_bit(255, true);
        }
    }

    /// Multiply `self` by `other` (mod `modulo`) via the Montgomery
    /// multiplication method.
    pub fn mul(&mut self, other: &U256, modulo: &U256, inv: u128) {
        if mul_reduce(&mut self.0, &other.0, &modulo.0, inv) {
            // has carry
            self.add_carry(modulo);
        }

        if *self >= *modulo {
            sub_noborrow(&mut self.0, &modulo.0);
        }
    }

    /// Turn `self` into its additive inverse (mod `modulo`)
    pub fn neg(&mut self, modulo: &U256) {
        if *self > Self::zero() {
            let mut tmp = modulo.0;
            sub_noborrow(&mut tmp, &self.0);

            self.0 = tmp;
        }
    }

    #[inline]
    pub fn is_even(&self) -> bool {
        self.0[0] & 1 == 0
    }

    /// Turn `self` into its multiplicative inverse (mod `modulo`)
    pub fn invert(&mut self, modulo: &U256) {
        // Guajardo Kumar Paar Pelzl
        // Efficient Software-Implementation of Finite Fields with Applications to Cryptography
        // Algorithm 16 (BEA for Inversion in Fp)

        let mut u = *self;
        let mut v = *modulo;
        let mut b = U256::one();
        let mut c = U256::zero();

        while u != U256::one() && v != U256::one() {
            while u.is_even() {
                div2(&mut u.0);
                b.div2(modulo);
            }
            while v.is_even() {
                div2(&mut v.0);
                c.div2(modulo);
            }

            if u >= v {
                sub_noborrow(&mut u.0, &v.0);
                b.sub(&c, modulo);
            } else {
                sub_noborrow(&mut v.0, &u.0);
                c.sub(&b, modulo);
            }
        }

        if u == U256::one() {
            self.0 = b.0;
        } else {
            self.0 = c.0;
        }
    }

    /// Return an Iterator<Item=bool> over all bits from
    /// MSB to LSB.
    pub fn bits(&self) -> BitIterator {
        BitIterator { int: &self, n: 256 }
    }
}

pub struct BitIterator<'a> {
    int: &'a U256,
    n: usize,
}

impl<'a> Iterator for BitIterator<'a> {
    type Item = bool;

    fn next(&mut self) -> Option<bool> {
        if self.n == 0 {
            None
        } else {
            self.n -= 1;

            self.int.get_bit(self.n)
        }
    }
}

impl Ord for U256 {
    #[inline]
    fn cmp(&self, other: &U256) -> Ordering {
        for (a, b) in self.0.iter().zip(other.0.iter()).rev() {
            if *a < *b {
                return Ordering::Less;
            } else if *a > *b {
                return Ordering::Greater;
            }
        }

        return Ordering::Equal;
    }
}

impl PartialOrd for U256 {
    #[inline]
    fn partial_cmp(&self, other: &U256) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// 
impl fmt::Debug for U256 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "U256({:032X}{:032X})", self.0[1], self.0[0])
    }
}

#[derive(Debug)]
pub enum Error {
    InvalidLength { expected: usize, actual: usize },
}

