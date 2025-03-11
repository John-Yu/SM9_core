use crate::{arith::*, u512::U512};
use ark_ff::{biginteger::BigInteger256 as B256, BigInteger as _};
use byteorder::{BigEndian, ByteOrder};
use core::{
    cmp::Ordering,
    fmt,
    ops::{Index, IndexMut},
};
use rand::Rng;

/// 256-bit stack allocated biginteger for use in prime field arithmetic.
#[derive(Default, Copy, Clone, PartialEq, Eq)]
#[repr(C)]
pub struct U256(pub(crate) B256);

impl From<[u64; 4]> for U256 {
    fn from(d: [u64; 4]) -> Self {
        U256(B256::new(d))
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

        let mut d = [0; 4];
        for (j, i) in (0..4).rev().zip((0..4).map(|i| i * 8)) {
            d[j] = BigEndian::read_u64(&s[i..]);
        }

        Ok(U256::from(d))
    }
    /// Converts a U256 into a slice of bytes (big endian)
    pub fn to_big_endian(self, s: &mut [u8]) -> Result<(), Error> {
        if s.len() != 32 {
            return Err(Error::InvalidLength {
                expected: 32,
                actual: s.len(),
            });
        }
        let be = self.0.to_bytes_be();
        s.copy_from_slice(be.as_ref());

        Ok(())
    }

    #[inline]
    pub fn zero() -> U256 {
        U256(B256::zero())
    }
    /// Returns true if element is zero.
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    #[inline]
    pub fn one() -> U256 {
        U256(B256::one())
    }
    /// Returns true if element is one.
    pub fn is_one(&self) -> bool {
        self.0 == B256::one()
    }
    /// Produce a random number (mod `modulo`)
    pub fn random<R: Rng>(rng: &mut R, modulo: &U256) -> U256 {
        U512::random(rng).divrem(modulo).1
    }

    pub fn set_bit(&mut self, n: usize, to: bool) -> bool {
        if n >= 256 {
            false
        } else {
            let limb = n >> 6;
            let bit = n & 0x3f;

            if to {
                self.0 .0[limb] |= 1 << bit;
            } else {
                self.0 .0[limb] &= !(1 << bit);
            }

            true
        }
    }

    pub fn get_bit(&self, n: usize) -> Option<bool> {
        if n >= 256 {
            None
        } else {
            Some(self.0.get_bit(n))
        }
    }

    #[inline]
    pub(crate) fn subtract_modulus_with_carry(&mut self, modulo: &U256, carry: bool) {
        if carry || self.0 >= modulo.0 {
            self.0.sub_with_borrow(&modulo.0);
        }
    }
    // self = self + 2^256 (mod `modulo`)
    pub(crate) fn add_carry(&mut self, modulo: &U256) {
        while !self.0.sub_with_borrow(&modulo.0) {}
    }
    /// Add `other` to `self` (mod `modulo`)
    pub fn add(&mut self, other: &U256, modulo: &U256) {
        let carry = self.0.add_with_carry(&other.0);
        self.subtract_modulus_with_carry(modulo, carry);
    }

    /// Subtract `other` from `self` (mod `modulo`)
    pub fn sub(&mut self, other: &U256, modulo: &U256) {
        // If `other` is larger than `self`, add the modulus to self first.
        if self.0 < other.0 {
            self.0.add_with_carry(&modulo.0);
        }
        self.0.sub_with_borrow(&other.0);
    }

    /// a = a * 2 (mod `modulo`)
    pub fn mul2(&mut self, modulo: &U256) {
        let c = self.0.mul2();
        self.subtract_modulus_with_carry(modulo, c);
    }

    /// a = a / 2 (mod `modulo`)
    pub fn div2(&mut self, modulo: &U256) {
        let mut carry = false;
        // If is odd, add the modulus to self first.
        if self.0.is_odd() {
            carry = self.0.add_with_carry(&modulo.0);
        }
        self.0.div2();
        if carry {
            self.set_bit(255, true);
            self.subtract_modulus_with_carry(modulo, false);
        }
    }

    /// Multiply `self` by `other` (mod `modulo`) via the Montgomery
    /// multiplication method.
    pub fn mul(&mut self, other: &U256, modulo: &U256, inv: u64) {
        let (carry, mut res) = self.mul_without_cond_subtract(other, modulo, inv);
        res.subtract_modulus_with_carry(modulo, carry);
        self.0 = res.0;
    }

    /// Square `self`  (mod `modulo`) via the Montgomery
    pub fn square(&mut self, modulo: &U256, inv: u64) {
        let mut r = MulBuffer::<4>::zeroed();
        let a = self.as_mut();

        let mut carry = 0;
        for i in 0..3 {
            for j in (i + 1)..4 {
                r[i + j] = mac_with_carry!(r[i + j], a[i], a[j], &mut carry);
            }
            r.b1[i] = carry;
            carry = 0;
        }

        r.b1[3] = r.b1[2] >> 63;
        for i in 2..7 {
            r[8 - i] = (r[8 - i] << 1) | (r[8 - (i + 1)] >> 63);
        }
        r.b0[1] <<= 1;
        for (i, ai) in a.iter().enumerate().take(4) {
            r[2 * i] = mac_with_carry!(r[2 * i], *ai, *ai, &mut carry);
            r[2 * i + 1] = adc!(r[2 * i + 1], 0, &mut carry);
        }

        // Montgomery reduction
        let mut carry2 = 0;
        let m = modulo.as_ref();
        for i in 0..4 {
            let k = r[i].wrapping_mul(inv);
            let mut carry = 0;
            mac_discard(r[i], k, m[0], &mut carry);
            for (j, mj) in m.iter().enumerate().take(4).skip(1) {
                r[i + j] = mac_with_carry!(r[i + j], k, *mj, &mut carry);
            }
            r.b1[i] = adc!(r.b1[i], carry, &mut carry2);
        }
        a.copy_from_slice(&r.b1);

        self.subtract_modulus_with_carry(modulo, carry2 != 0);
    }

    /// Turn `self` into its additive inverse (mod `modulo`)
    #[inline]
    pub fn neg(&mut self, modulo: &U256) {
        if !self.is_zero() {
            let mut tmp = modulo.0;
            tmp.sub_with_borrow(&self.0);
            self.0 = tmp;
        }
    }

    #[inline]
    pub fn is_even(&self) -> bool {
        self.0.is_even()
    }
    #[inline]
    pub fn is_odd(&self) -> bool {
        self.0.is_odd()
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
                u.0.div2();
                b.div2(modulo);
            }
            while v.is_even() {
                v.0.div2();
                c.div2(modulo);
            }

            if u >= v {
                u.0.sub_with_borrow(&v.0);
                b.sub(&c, modulo);
            } else {
                v.0.sub_with_borrow(&u.0);
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
    fn mul_without_cond_subtract(mut self, other: &Self, modulo: &U256, inv: u64) -> (bool, Self) {
        let mut r = MulBuffer::<4>::zeroed();
        let e = other.as_ref();
        let d = self.as_ref();
        for (i, di) in d.iter().enumerate().take(4) {
            let mut carry = 0;
            for (j, ej) in e.iter().enumerate().take(4) {
                let k = i + j;
                r[k] = mac_with_carry!(r[k], *di, *ej, &mut carry);
            }
            r.b1[i] = carry;
        }
        // Montgomery reduction
        let mut carry2 = 0;
        let m = modulo.as_ref();
        for i in 0..4 {
            let tmp = r[i].wrapping_mul(inv);
            let mut carry = 0;
            mac_discard(r[i], tmp, m[0], &mut carry);
            for (j, mj) in m.iter().enumerate().take(4).skip(1) {
                let k = i + j;
                r[k] = mac_with_carry!(r[k], tmp, *mj, &mut carry);
            }
            r.b1[i] = adc!(r.b1[i], carry, &mut carry2);
        }
        self.as_mut().copy_from_slice(&r.b1);

        (carry2 != 0, self)
    }
}

pub struct BitIterator<'a> {
    int: &'a U256,
    n: usize,
}

impl Iterator for BitIterator<'_> {
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
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for U256 {
    #[inline]
    fn partial_cmp(&self, other: &U256) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl AsMut<[u64]> for U256 {
    #[inline]
    fn as_mut(&mut self) -> &mut [u64] {
        self.0.as_mut()
    }
}

impl AsRef<[u64]> for U256 {
    #[inline]
    fn as_ref(&self) -> &[u64] {
        self.0.as_ref()
    }
}

impl Index<usize> for U256 {
    type Output = u64;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0 .0[index]
    }
}

impl IndexMut<usize> for U256 {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0 .0[index]
    }
}

//
impl fmt::Debug for U256 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "U256({:016X}{:016X}{:016X}{:016X})",
            self[3], self[2], self[1], self[0]
        )
    }
}

#[derive(Debug)]
pub enum Error {
    InvalidLength { expected: usize, actual: usize },
}
