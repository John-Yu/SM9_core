use byteorder::{BigEndian, ByteOrder};
use core::{
    cmp::Ordering,
    fmt,
    ops::{Index, IndexMut},
};
use rand::Rng;

use crate::arith::*;
use crate::u512::U512;

/// 256-bit, stack allocated biginteger for use in prime field
/// arithmetic.
#[derive(Copy, Clone, PartialEq, Eq)]
#[repr(C)]
pub struct U256(pub(crate) [u128; 2]);

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
    pub fn to_big_endian(self, s: &mut [u8]) -> Result<(), Error> {
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
    /// Returns true if element is zero.
    pub fn is_zero(&self) -> bool {
        self[0] == 0 && self[1] == 0
    }
    /// Returns true if element is one.
    pub fn is_one(&self) -> bool {
        self[0] == 1 && self[1] == 0
    }

    pub fn set_bit(&mut self, n: usize, to: bool) -> bool {
        if n >= 256 {
            false
        } else {
            let part = n >> 7;
            let bit = n & 0x7F;

            if to {
                self[part] |= 1 << bit;
            } else {
                self[part] &= !(1 << bit);
            }

            true
        }
    }

    pub fn get_bit(&self, n: usize) -> Option<bool> {
        if n >= 256 {
            None
        } else {
            let part = n >> 7;
            let bit = n & 0x7F;

            Some(self[part] & (1 << bit) > 0)
        }
    }

    // self = self + 2^256 (mod `modulo`)
    pub(crate) fn add_carry(&mut self, modulo: &U256) {
        loop {
            if sub_borrow(&mut self.0, &modulo.0) {
                break;
            }
        }
    }
    /// Add `other` to `self` (mod `modulo`)
    pub fn add(&mut self, other: &U256, modulo: &U256) {
        if add_carry(&mut self.0, &other.0) {
            // has carry
            self.add_carry(modulo);
        } else if *self >= *modulo {
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
        } else {
            sub_noborrow(&mut self.0, &other.0);
        }
    }

    /// a = a * 2 (mod `modulo`)
    pub fn mul2(&mut self, modulo: &U256) {
        if mul2(&mut self.0) {
            // has carry
            self.add_carry(modulo);
        } else if *self >= *modulo {
            sub_noborrow(&mut self.0, &modulo.0);
        }
    }

    /// a = a / 2 (mod `modulo`)
    pub fn div2(&mut self, modulo: &U256) {
        let mut carry = false;
        if !self.is_even() {
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
        } else if *self >= *modulo {
            sub_noborrow(&mut self.0, &modulo.0);
        }
    }

    /// Square `self`  (mod `modulo`) via the Montgomery
    pub fn square(&mut self, modulo: &U256, inv: u128) {
        if square_reduce(&mut self.0, &modulo.0, inv) {
            // has carry
            self.add_carry(modulo);
        } else if *self >= *modulo {
            sub_noborrow(&mut self.0, &modulo.0);
        }
    }

    /// Turn `self` into its additive inverse (mod `modulo`)
    pub fn neg(&mut self, modulo: &U256) {
        if !self.is_zero() {
            let mut tmp = modulo.0;
            sub_noborrow(&mut tmp, &self.0);
            self.0 = tmp;
        }
    }

    #[inline]
    pub fn is_even(&self) -> bool {
        self[0] & 1 == 0
    }

    /// Turn `self` into its multiplicative inverse (mod `modulo`)
    // Guajardo Kumar Paar Pelzl
    // Efficient Software-Implementation of Finite Fields with Applications to Cryptography
    // Algorithm 16 (BEA for Inversion in Fp)
    pub fn invert(&mut self, modulo: &U256, rsquared: &U256) {
        let mut u = *self;
        let mut v = *modulo;
        let mut b = *rsquared; // Avoids unnecessary reduction step. perfect!
        let mut c = U256::zero();

        while !u.is_one() && !v.is_one() {
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

        if u.is_one() {
            self.0 = b.0;
        } else {
            self.0 = c.0;
        }
    }

    /// Return an Iterator<Item=bool> over all bits from
    /// MSB to LSB.
    pub fn bits(&self) -> BitIterator {
        BitIterator { int: self, n: 256 }
    }
    /// Construct an iterator that automatically skips any leading zeros.
    /// That is, it skips all zeros before the most-significant one.
    pub fn bits_without_leading_zeros(&self) -> impl Iterator<Item = bool> + '_ {
        BitIterator { int: self, n: 256 }.skip_while(|b| !b)
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
        match self[1].cmp(&other[1]) {
            Ordering::Greater => Ordering::Greater,
            Ordering::Less => Ordering::Less,
            Ordering::Equal => self[0].cmp(&other[0]),
        }
    }
}

impl PartialOrd for U256 {
    #[inline]
    fn partial_cmp(&self, other: &U256) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl AsMut<[u128]> for U256 {
    #[inline]
    fn as_mut(&mut self) -> &mut [u128] {
        &mut self.0
    }
}

impl AsRef<[u128]> for U256 {
    #[inline]
    fn as_ref(&self) -> &[u128] {
        &self.0
    }
}

impl Index<usize> for U256 {
    type Output = u128;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IndexMut<usize> for U256 {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
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
