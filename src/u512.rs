use ark_ff::{BigInt, BigInteger as _};
use byteorder::{BigEndian, ByteOrder};
use core::{cmp::Ordering, fmt, ops::Index};
use rand::Rng;

use crate::u256::{Error, U256};

type B512 = BigInt<8>;
/// 512-bit, stack allocated biginteger for use in extension
/// field serialization and scalar interpretation.
#[derive(Copy, Clone, PartialEq, Eq)]
#[repr(C)]
pub struct U512(pub(crate) B512);

impl From<[u64; 8]> for U512 {
    fn from(d: [u64; 8]) -> Self {
        U512(B512::new(d))
    }
}

impl U512 {
    /// Multiplies c1 by modulo, adds c0.
    pub fn new(c1: &U256, c0: &U256, modulo: &U256) -> U512 {
        let (mut low, mut high) = c1.0.mul(&modulo.0);
        let mut carry = low.add_with_carry(&c0.0);
        if carry {
            carry = high.add_with_carry(&U256::one().0);
        }
        debug_assert!(!carry);
        U512::from_slice(&[high.to_bytes_be(), low.to_bytes_be()].concat()).unwrap()
    }

    #[inline]
    pub fn one() -> U512 {
        U512(B512::one())
    }

    pub fn from_slice(s: &[u8]) -> Result<U512, Error> {
        if s.len() != 64 {
            return Err(Error::InvalidLength {
                expected: 64,
                actual: s.len(),
            });
        }

        let mut d = [0; 8];
        for (j, i) in (0..8).rev().zip((0..8).map(|i| i * 8)) {
            d[j] = BigEndian::read_u64(&s[i..]);
        }

        Ok(U512::from(d))
    }

    /// Get a random U512
    pub fn random<R: Rng>(rng: &mut R) -> U512 {
        U512(rng.gen())
    }
    /// Returns the number of bits necessary to represent this value. Note that zero
    /// is considered to need 0 bits.
    pub fn bit_length(&self) -> usize {
        self.0.num_bits() as usize
    }
    /// Returns the `n`-th bit where bit 0 is the least significant one. true if bit is 1.
    /// In other words, the bit with weight `2^n`.
    #[inline]
    pub fn get_bit(&self, n: usize) -> Option<bool> {
        if n >= 512 {
            None
        } else {
            Some(self.0.get_bit(n))
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
            let carry = r.0.mul2();
            r.set_bit(0, self.0.get_bit(i));
            if &r >= modulo || carry {
                r.0.sub_with_borrow(&modulo.0);
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
        U512::from_slice(buf).unwrap()
    }
}

impl Ord for U512 {
    #[inline]
    fn cmp(&self, other: &U512) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for U512 {
    #[inline]
    fn partial_cmp(&self, other: &U512) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Index<usize> for U512 {
    type Output = u64;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0 .0[index]
    }
}

impl fmt::Debug for U512 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "U512({:016X}{:016X}{:016X}{:016X}{:016X}{:016X}{:016X}{:016X})",
            self[7], self[6], self[5], self[4], self[3], self[2], self[1], self[0]
        )
    }
}
