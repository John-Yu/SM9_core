#![allow(dead_code)]

use core::ops::{Add, Mul, Neg, Sub};
use rand::Rng;

use crate::fields::{FieldElement, Fq, Fq2};
use crate::u256::U256;

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq4 {
    c0: Fq2,
    c1: Fq2,
}

impl Fq4 {
    pub fn new(c0: Fq2, c1: Fq2) -> Self {
        Fq4 { c0, c1 }
    }

    //Algorithm 7
    pub fn scale(&self, by: &Fq2) -> Self {
        Fq4 {
            c0: self.c0 * *by,
            c1: self.c1 * *by,
        }
    }

    pub fn scale_fq(&self, by: &Fq) -> Self {
        Fq4 {
            c0: self.c0.scale(by),
            c1: self.c1.scale(by),
        }
    }

    pub fn mul_by_nonresidue(&self) -> Self {
        //c0 = c1 * u
        //c1 = c0
        Fq4 {
            c0: self.c1.mul_by_nonresidue(),
            c1: self.c0,
        }
    }

    pub fn frobenius_map(&self, power: usize) -> Self {
        match power {
            10 => Fq4 {
                c0: self.c0.unitary_inverse(),
                c1: self
                    .c1
                    .unitary_inverse()
                    .scale(&Fq::new(*SM9_ALPHA3).unwrap()),
            },
            11 => Fq4 {
                c0: self
                    .c0
                    .unitary_inverse()
                    .scale(&Fq::new(*SM9_ALPHA1).unwrap()),
                c1: self
                    .c1
                    .unitary_inverse()
                    .scale(&Fq::new(*SM9_ALPHA4).unwrap()),
            },
            12 => Fq4 {
                c0: self
                    .c0
                    .unitary_inverse()
                    .scale(&Fq::new(*SM9_ALPHA2).unwrap()),
                c1: self
                    .c1
                    .unitary_inverse()
                    .scale(&Fq::new(*SM9_ALPHA5).unwrap()),
            },
            21 => self
                .unitary_inverse()
                .scale_fq(&Fq::new(*SM9_ALPHA2).unwrap()),
            22 => self
                .unitary_inverse()
                .scale_fq(&Fq::new(*SM9_ALPHA4).unwrap()),
            30 => Fq4 {
                c0: self.c0.unitary_inverse(),
                c1: self
                    .c1
                    .unitary_inverse()
                    .mul(Fq2::new(Fq::new(*SM9_BETA).unwrap(), Fq::zero()))
                    .neg(),
            },
            31 => Fq4 {
                c0: self
                    .c0
                    .unitary_inverse()
                    .mul(Fq2::new(Fq::new(*SM9_BETA).unwrap(), Fq::zero())),
                c1: self.c1.unitary_inverse(),
            },
            32 => Fq4 {
                c0: self.c0.unitary_inverse().neg(),
                c1: self
                    .c1
                    .unitary_inverse()
                    .mul(Fq2::new(Fq::new(*SM9_BETA).unwrap(), Fq::zero())),
            },
            _ => unimplemented!(),
        }
    }

    pub fn unitary_inverse(&self) -> Fq4 {
        Fq4::new(self.c0, -self.c1)
    }
    /// Converts an element into a byte representation in
    /// big-endian byte order.
    pub fn to_slice(self) -> [u8; 128] {
        let mut res = [0u8; 128];
        let b1 = self.c1.to_slice();
        let b0 = self.c0.to_slice();
        res[..64].copy_from_slice(&b1);
        res[64..].copy_from_slice(&b0);
        res
    }
}

impl FieldElement for Fq4 {
    fn zero() -> Self {
        Fq4 {
            c0: Fq2::zero(),
            c1: Fq2::zero(),
        }
    }

    fn one() -> Self {
        Fq4 {
            c0: Fq2::one(),
            c1: Fq2::zero(),
        }
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        Fq4 {
            c0: Fq2::random(rng),
            c1: Fq2::random(rng),
        }
    }
    fn is_zero(&self) -> bool {
        self.c0.is_zero() && self.c1.is_zero()
    }

    /// Squares this element
    fn squared(&self) -> Self {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Complex squaring), which takes 2M + 3A + 2B
        let a0 = self.c0;
        let a1 = self.c1;
        let v0 = a0 * a1;
        Fq4 {
            c0: (a0 + a1) * (a0 + a1.mul_by_nonresidue()) - v0 - v0.mul_by_nonresidue(),
            c1: v0 + v0,
        }
    }
    //"High-Speed Software Implementation of the Optimal Ate AbstractPairing over Barreto-Naehrig Curves"
    // Algorithm 23
    fn inverse(self) -> Option<Self> {
        if self.is_zero() {
            None
        } else {
            (self.c0.squared() - self.c1.squared().mul_by_nonresidue())
                .inverse()
                .map(|t| Self::new(self.c0 * t, -(self.c1 * t)))
        }
    }
}

impl Mul for Fq4 {
    type Output = Fq4;

    fn mul(self, other: Fq4) -> Fq4 {
        // Devegili OhEig Scott Dahab
        //     Multiplication and Squaring on Pairing-Friendly Fields.pdf
        //     Section 3 (Karatsuba), which costs 3M + 3A + 2B
        let aa = self.c0 * other.c0;
        let bb = self.c1 * other.c1;

        Fq4 {
            c0: aa + bb.mul_by_nonresidue(),
            c1: (self.c0 + self.c1) * (other.c0 + other.c1) - aa - bb,
        }
    }
}

impl Sub for Fq4 {
    type Output = Fq4;

    fn sub(self, other: Fq4) -> Fq4 {
        Fq4 {
            c0: self.c0 - other.c0,
            c1: self.c1 - other.c1,
        }
    }
}

impl Add for Fq4 {
    type Output = Fq4;

    fn add(self, other: Fq4) -> Fq4 {
        Fq4 {
            c0: self.c0 + other.c0,
            c1: self.c1 + other.c1,
        }
    }
}

impl Neg for Fq4 {
    type Output = Fq4;

    fn neg(self) -> Fq4 {
        Fq4 {
            c0: -self.c0,
            c1: -self.c1,
        }
    }
}

lazy_static::lazy_static! {

// beta   = 0x6c648de5dc0a3f2c f55acc93ee0baf15 9f9d411806dc5177 f5b21fd3da24d011
// alpha1 = 0x3f23ea58e5720bdb 843c6cfa9c086749 47c5c86e0ddd04ed a91d8354377b698b
// alpha2 = 0xf300000002a3a6f2 780272354f8b78f4 d5fc11967be65334
// alpha3 = 0x6c648de5dc0a3f2c f55acc93ee0baf15 9f9d411806dc5177 f5b21fd3da24d011
// alpha4 = 0xf300000002a3a6f2 780272354f8b78f4 d5fc11967be65333
// alpha5 = 0x2d40a38cf6983351 711e5f99520347cc 57d778a9f8ff4c8a 4c949c7fa2a96686
    static ref SM9_BETA: U256 = U256::from([
        0xf5b21fd3da24d011,
        0x9f9d411806dc5177,
        0xf55acc93ee0baf15,
        0x6c648de5dc0a3f2c
    ]);
    static ref SM9_ALPHA1: U256 = U256::from([
        0xa91d8354377b698b,
        0x47c5c86e0ddd04ed,
        0x843c6cfa9c086749,
        0x3f23ea58e5720bdb
    ]);
    static ref SM9_ALPHA2: U256 = U256::from([
        0xd5fc11967be65334,
        0x780272354f8b78f4,
        0xf300000002a3a6f2,
        0x0
    ]);
    static ref SM9_ALPHA3: U256 = U256::from([
        0xf5b21fd3da24d011,
        0x9f9d411806dc5177,
        0xf55acc93ee0baf15,
        0x6c648de5dc0a3f2c
    ]);
    static ref SM9_ALPHA4: U256 = U256::from([
        0xd5fc11967be65333,
        0x780272354f8b78f4,
        0xf300000002a3a6f2,
        0x0
    ]);
    static ref SM9_ALPHA5: U256 = U256::from([
        0x4c949c7fa2a96686,
        0x57d778a9f8ff4c8a,
        0x711e5f99520347cc,
        0x2d40a38cf6983351
    ]);

}
