#![allow(dead_code)]

#[macro_use]
mod utils;
mod fp;
mod fq12;
mod fq2;
mod fq4;

use alloc::fmt::Debug;
use core::ops::{Add, Mul, MulAssign, Neg, Sub};
use rand::Rng;

use crate::{u256::U256, One, Zero};

pub use self::fp::{Fq, Fr};
pub use self::fq12::Fq12;
pub use self::fq2::Fq2;
pub use self::fq4::Fq4;

pub trait FieldElement:
    Sized
    + Copy
    + Clone
    + Zero
    + One
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Neg<Output = Self>
    + PartialEq
    + Eq
    + Debug
{
    //fn one() -> Self;
    fn random<R: Rng>(_: &mut R) -> Self;
    fn squared(&self) -> Self;
    // double this element
    fn double(&self) -> Self;
    // triple this element
    fn triple(&self) -> Self;
    // Computes the inverse of this element
    fn inverse(&self) -> Option<Self>;
    // Left-to-right binary modular exponentiation with square-and-multiply method.
    fn pow<I: Into<U256>>(&self, by: I) -> Self
    where
        for<'a> Self: MulAssign<&'a Self>,
    {
        let mut res = Self::one();
        for i in by.into().bits_without_leading_zeros() {
            res = res.squared();
            if i {
                res *= self;
            }
        }

        res
    }
}

lazy_static::lazy_static! {

    static ref FR: U256 = U256::from([
        0xE56EE19CD69ECF25,
        0x49F2934B18EA8BEE,
        0xD603AB4FF58EC744,
        0xB640000002A3A6F1
    ]);

    static ref FR_SQUARED: U256 = U256::from([
        0x7598CD79CD750C35,
        0xE4A08110BB6DAEAB,
        0xBFEE4BAE7D78A1F9,
        0x8894F5D163695D0E
    ]);

    static ref FR_ONE: U256 = U256::from([
        0x1A911E63296130DB,
        0xB60D6CB4E7157411,
        0x29FC54B00A7138BB,
        0x49BFFFFFFD5C590E
    ]);

    static ref FQ: U256 = U256::from([
        0xE56F9B27E351457D,
        0x21F2934B1A7AEEDB,
        0xD603AB4FF58EC745,
        0xB640000002A3A6F1
    ]);

    static ref FQ_SQUARED: U256 = U256::from([
        0x27DEA312B417E2D2,
        0x88F8105FAE1A5D3F,
        0xE479B522D6706E7B,
        0x2EA795A656F62FBD
    ]);

    static ref FQ_ONE: U256 = U256::from([
        0x1A9064D81CAEBA83,
        0xDE0D6CB4E5851124,
        0x29FC54B00A7138BA,
        0x49BFFFFFFD5C590E
    ]);

    static ref FQ_INV: u64 = 0x892BC42C2F2EE42B;

    static ref FR_INV: u64 = 0x1D02662351974B53;

    pub static ref FQ_MINUS1_DIV4: Fq =
        Fq::new(1.into()).expect("1 is a valid field element and static; qed").neg() *
        Fq::new(4.into()).expect("4 is a valid field element and static; qed").inverse()
            .expect("4 has inverse in Fq and is static; qed");

    pub static ref FQ_MINUS5_DIV8: Fq =
        Fq::new(5.into()).expect("5 is a valid field element and static; qed").neg() *
        Fq::new(8.into()).expect("8 is a valid field element and static; qed").inverse()
            .expect("8 has inverse in Fq and is static; qed");

}

/************************************************************************************************ */
#[cfg(test)]
mod tests {
    use crate::u512::U512;
    use fp::Fq;
    use hex_literal::hex;

    use super::*;

    #[test]
    fn u256_from_slice() {
        // println!("U256 from_slice test");
        let tst = U256::one();
        let mut s = [0u8; 32];
        s[31] = 1; //BigEndian

        let num = U256::from_slice(&s)
            .expect("U256 should initialize ok from slice in `from_slice` test");
        assert_eq!(num, tst);
    }

    #[test]
    fn u256_to_big_endian() {
        // println!("U256 to_big_endian test");
        let num = U256::one();
        let mut s = [0u8; 32];

        num.to_big_endian(&mut s)
            .expect("U256 should convert to bytes ok in `to_big_endian` test");
        assert_eq!(
            s,
            [
                0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8,
                0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 0u8, 1u8,
            ]
        );
    }
    #[test]
    fn u512_from_slice() {
        let tst = U512::one();
        let mut s = [0u8; 64];
        s[63] = 1; //BigEndian

        let num = U512::from_slice(&s)
            .expect("U512 should initialize ok from slice in `from_slice` test");
        assert_eq!(num, tst);
    }

    #[test]
    fn testing_divrem() {
        let modulo = U256::from([
            0x3c208c16d87cfd47,
            0x97816a916871ca8d,
            0xb85045b68181585d,
            0x30644e72e131a029,
        ]);

        let mut s = [0u8; 32];
        modulo.to_big_endian(&mut s).unwrap();
        // Modulus should become 1*q + 0
        let a = U512::from([
            0x3c208c16d87cfd47,
            0x97816a916871ca8d,
            0xb85045b68181585d,
            0x30644e72e131a029,
            0,
            0,
            0,
            0,
        ]);
        let (c1, c0) = a.divrem(&modulo);
        assert_eq!(c1.unwrap(), U256::one());
        assert_eq!(c0, U256::zero());

        {
            // Modulus squared minus 1 should be (q-1) q + q-1
            let a = U512::from([
                0x3b5458a2275d69b0,
                0xa602072d09eac101,
                0x4a50189c6d96cadc,
                0x04689e957a1242c8,
                0x26edfa5c34c6b38d,
                0xb00b855116375606,
                0x599a6f7c0348d21c,
                0x0925c4b8763cbf9c,
            ]);

            let (c1, c0) = a.divrem(&modulo);
            assert_eq!(
                c1.unwrap(),
                U256::from([
                    0x3c208c16d87cfd46,
                    0x97816a916871ca8d,
                    0xb85045b68181585d,
                    0x30644e72e131a029
                ])
            );
            assert_eq!(
                c0,
                U256::from([
                    0x3c208c16d87cfd46,
                    0x97816a916871ca8d,
                    0xb85045b68181585d,
                    0x30644e72e131a029
                ])
            );
        }

        {
            // Modulus squared minus 2 should be (q-1) q + q-2
            let a = U512::from([
                0x3b5458a2275d69af,
                0xa602072d09eac101,
                0x4a50189c6d96cadc,
                0x04689e957a1242c8,
                0x26edfa5c34c6b38d,
                0xb00b855116375606,
                0x599a6f7c0348d21c,
                0x0925c4b8763cbf9c,
            ]);

            let (c1, c0) = a.divrem(&modulo);

            assert_eq!(
                c1.unwrap(),
                U256::from([
                    0x3c208c16d87cfd46,
                    0x97816a916871ca8d,
                    0xb85045b68181585d,
                    0x30644e72e131a029
                ])
            );
            assert_eq!(
                c0,
                U256::from([
                    0x3c208c16d87cfd45,
                    0x97816a916871ca8d,
                    0xb85045b68181585d,
                    0x30644e72e131a029
                ])
            );
        }
        {
            // Ridiculously large number should fail
            let a = U512::from([
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
            ]);

            let (c1, c0) = a.divrem(&modulo);
            assert!(c1.is_none());
            assert_eq!(
                c0,
                U256::from([
                    0xf32cfc5b538afa88,
                    0xb5e71911d44501fb,
                    0x47ab1eff0a417ff6,
                    0x06d89f71cab8351f
                ])
            );
        }

        {
            // Modulus squared should fail
            let a = U512::from([
                0x3b5458a2275d69b1,
                0xa602072d09eac101,
                0x4a50189c6d96cadc,
                0x04689e957a1242c8,
                0x26edfa5c34c6b38d,
                0xb00b855116375606,
                0x599a6f7c0348d21c,
                0x0925c4b8763cbf9c,
            ]);

            let (c1, c0) = a.divrem(&modulo);
            assert!(c1.is_none());
            assert_eq!(c0, U256::zero());
        }

        {
            // Modulus squared plus one should fail
            let a = U512::from([
                0x3b5458a2275d69b2,
                0xa602072d09eac101,
                0x4a50189c6d96cadc,
                0x04689e957a1242c8,
                0x26edfa5c34c6b38d,
                0xb00b855116375606,
                0x599a6f7c0348d21c,
                0x0925c4b8763cbf9c,
            ]);

            let (c1, c0) = a.divrem(&modulo);
            assert!(c1.is_none());
            assert_eq!(c0, U256::one());
        }

        {
            let modulo = U256::from([
                0x43e1f593f0000001,
                0x2833e84879b97091,
                0xb85045b68181585d,
                0x30644e72e131a029,
            ]);

            // Fr modulus masked off is valid
            let a = U512::from([
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0xffffffffffffffff,
                0x07ffffffffffffff,
            ]);

            let (c1, c0) = a.divrem(&modulo);

            assert!(c1.unwrap() < modulo);
            assert!(c0 < modulo);
        }
    }

    fn can_invert<F: FieldElement>() {
        let mut a = F::one();

        for _ in 0..1000 {
            assert_eq!(a * a.inverse().unwrap(), F::one());

            a = a + F::one();
        }

        a = -F::one();
        for _ in 0..1000 {
            assert_eq!(a * a.inverse().unwrap(), F::one());

            a = a - F::one();
        }

        assert_eq!(F::zero().inverse(), None);
    }

    #[test]
    fn test_fq() {
        // println!("can_invert::<Fq> test");
        can_invert::<Fq>();
        assert_eq!(-Fq::zero(), Fq::zero());
        assert_eq!(-Fq::one() + Fq::one(), Fq::zero());
        assert_eq!(Fq::zero() - Fq::zero(), Fq::zero());
        let a = Fq::new(1.into()).unwrap();
        assert_eq!(a, Fq::one());
        let a = Fq::from_str("1").unwrap();
        assert_eq!(a, Fq::one());
    }
    #[test]
    fn test_sum_of_products() {
        use rand::{rngs::StdRng, SeedableRng};
        let seed = [
            0, 0, 0, 0, 0, 0, 64, 13, // 103245
            0, 0, 0, 0, 0, 0, 176, 2, // 191922
            0, 0, 0, 0, 0, 0, 0, 13, // 1293
            0, 0, 0, 0, 0, 0, 96, 7u8, // 192103
        ];
        let mut rng = StdRng::from_seed(seed);
        for _ in 0..500 {
            let a1 = Fq::random(&mut rng);
            let a2 = Fq::random(&mut rng);
            let a3 = Fq::random(&mut rng);
            let a4 = Fq::random(&mut rng);
            let b1 = Fq::random(&mut rng);
            let b2 = Fq::random(&mut rng);
            let b3 = Fq::random(&mut rng);
            let b4 = Fq::random(&mut rng);
            let c1 = Fq::sum_of_products(
                &[a1, a2, a3, a4, a1, a2, a3, a4],
                &[b1, b2, b3, b4, b1, b2, b3, b4],
            );
            let c2 = a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4;
            assert_eq!(c1, c2.double());
        }
    }

    const X: [u8; 32] =
        hex!("85AEF3D0 78640C98 597B6027 B441A01F F1DD2C19 0F5E93C4 54806C11 D8806141");
    const Y: [u8; 32] =
        hex!("41E00A53 DDA532DA 1A7CE027 B7A46F74 1006E85F 5CDFF073 0E75C05F B4E3216D");
    const HEX_IV: [u8; 32] =
        hex!("123456789abcdef00fedcba987654321123456789abcdef00fedcba987654321");
    const R_ADD: [u8; 32] =
        hex!("114efe24536598809df494ff7657484edff1812d51c3955b7d869149aa123d31");
    const R_SUB: [u8; 32] =
        hex!("43cee97c9abed9be3efe7ffffc9d30abe1d643b9b27ea351460aabb2239d3fd4");
    const R_NSUB: [u8; 32] =
        hex!("7271168367e4cd3397052b4ff8f19699401c4f9167fc4b8a9f64ef75bfb405a9");
    const R_MUL: [u8; 32] =
        hex!("9e4d19bb5d94a47352e6f53f4116b2a71b16a1113dc789b26528ee19f46b72e0");
    const R_POW: [u8; 32] =
        hex!("5679a8f0a46ada5b9d48008cde0b8b7a233f882c08afe8f08a36a20ac845bb1a");
    const R_INV: [u8; 32] =
        hex!("7d404b0027a93e3fa8f8bc7ee367a96814c42a3b69feb1845093406948a34753");
    const R_NEG: [u8; 32] =
        hex!("30910c2f8a3f9a597c884b28414d2725301567320b1c5b1790ef2f160ad0e43c");
    const R_SQRT: [u8; 32] =
        hex!("46dc2a5b8853234b341d9c57f9c4ca5709e95bbfef25356812e884e4f38cd0d6");
    const LOOPS: u64 = 0x3ff;

    #[test]
    fn test_fq_add() {
        // println!("test_fq_add test");
        let mut a = Fq::zero();
        for i in 1..LOOPS {
            a = a + Fq::one();
            assert_eq!(
                a,
                Fq::new(i.into()).expect("u64 is a valid field element and static; qed")
            );
        }
        let x = Fq::new(U256::from_slice(&X).unwrap()).unwrap();
        let y = Fq::new(U256::from_slice(&Y).unwrap()).unwrap();
        let r = x + y;
        assert_eq!(r, Fq::new(U256::from_slice(&R_ADD).unwrap()).unwrap());
    }
    #[test]
    fn test_fq_double() {
        let mut x = Fq::new(U256::from_slice(&X).unwrap()).unwrap();
        for _i in 1..LOOPS {
            assert_eq!(x.double(), x + x);
            x = x + x;
        }
    }

    #[test]
    fn test_fq_sub() {
        // println!("test_fq_sub test");
        let mut a = Fq::new(LOOPS.into()).expect("u64 is a valid field element and static; qed");
        for i in (1..LOOPS).rev() {
            a = a - Fq::one();
            assert_eq!(
                a,
                Fq::new(i.into()).expect("u64 is a valid field element and static; qed")
            );
        }
        let x = Fq::new(U256::from_slice(&X).unwrap()).unwrap();
        let y = Fq::new(U256::from_slice(&Y).unwrap()).unwrap();
        let r = x - y;
        assert_eq!(r, Fq::new(U256::from_slice(&R_SUB).unwrap()).unwrap());
        // println!("x = {:X?}", U256::from(x));
        // println!("y = {:X?}", U256::from(y));
        let r = y - x;
        // println!("r    = {:X?}", U256::from(r));
        let nsub = Fq::new(U256::from_slice(&R_NSUB).unwrap()).unwrap();
        // println!("nsub = {:X?}", U256::from(nsub));
        assert_eq!(r, nsub);
    }

    #[test]
    fn test_fq_mul() {
        // println!("test_fq_mul test");
        let x = Fq::new(U256::from_slice(&X).unwrap()).unwrap();
        let y = Fq::new(U256::from_slice(&Y).unwrap()).unwrap();
        let r = x * y;
        assert_eq!(r, Fq::new(U256::from_slice(&R_MUL).unwrap()).unwrap());
    }

    #[test]
    fn test_fq_pow() {
        // println!("test_fq_pow test");
        let x = Fq::new(U256::from_slice(&X).unwrap()).unwrap();
        let y = Fq::new(U256::from_slice(&Y).unwrap()).unwrap();
        let r = x.pow(y);
        assert_eq!(r, Fq::new(U256::from_slice(&R_POW).unwrap()).unwrap());
    }

    #[test]
    fn test_fq_inv() {
        // println!("test_fq_inv test");
        let x = Fq::new(U256::from_slice(&X).unwrap()).unwrap();
        let r = x.inverse().unwrap();
        assert_eq!(r, Fq::new(U256::from_slice(&R_INV).unwrap()).unwrap());
    }

    #[test]
    fn test_fq_neg() {
        // println!("test_fq_neg test");
        let x = Fq::new(U256::from_slice(&X).unwrap()).unwrap();
        let r = -x;
        assert_eq!(r, Fq::new(U256::from_slice(&R_NEG).unwrap()).unwrap());
    }

    #[test]
    fn test_fq_sqr() {
        // println!("test_fq_sqr test");
        let x = Fq::new(U256::from_slice(&X).unwrap()).unwrap();
        // println!("x   = {:?}", x);
        let r = x * x;
        // println!("x*x = {:?}", r);
        let sqr = Fq::new(U256::from_slice(&R_SQRT).unwrap()).unwrap();
        // println!("sqr = {:?}", sqr);

        assert_eq!(r, sqr);
    }
    const H40: [u8; 40] =
        hex!("8795E170FB592DE8BF5B2DB4181E8BA7EE4A3E50281FE2AE0DFA3AFFE018A2C0617FDDF42FB042A1");
    const HASH_H40: [u8; 32] =
        hex!("9CB1F6288CE0E51043CE72344582FFC301E0A812A7F5F2004B85547A24B82716");
    #[test]
    fn test_fr_from_hash() {
        let x = Fr::new(U256::from_slice(&HASH_H40).unwrap()).unwrap();
        let r = Fr::from_hash(&H40).unwrap();

        assert_eq!(r, x);
    }

    #[test]
    fn test_fq2() {
        // println!("can_invert::<Fq2> test");
        can_invert::<Fq2>();
        assert_eq!(-Fq2::zero(), Fq2::zero());
        assert_eq!(-Fq2::one() + Fq2::one(), Fq2::zero());
        assert_eq!(Fq2::zero() - Fq2::zero(), Fq2::zero());
    }

    const X1: [u8; 32] =
        hex!("17509B09 2E845C12 66BA0D26 2CBEE6ED 0736A96F A347C8BD 856DC76B 84EBEB96");
    const X0: [u8; 32] =
        hex!("A7CF28D5 19BE3DA6 5F317015 3D278FF2 47EFBA98 A71A0811 6215BBA5 C999A7C7");
    const Y1: [u8; 32] =
        hex!("9F64080B 3084F733 E48AFF4B 41B56501 1CE0711C 5E392CFB 0AB1B679 1B94C408");
    const Y0: [u8; 32] =
        hex!("29DBA116 152D1F78 6CE843ED 24A3B573 414D2177 386A92DD 8F14D656 96EA5E32");
    const R_ADD1: [u8; 32] =
        hex!("0074a3145c65ac547541612178e584a902248740e70606dcaaafe2bcbd2f6a21");
    const R_ADD0: [u8; 32] =
        hex!("1b6ac9eb2c47b62cf61608b26c3c7e20674a48c4c509ac130bbaf6d47d32c07c");
    const R_MUL0: [u8; 32] =
        hex!("27fe3a559abcc3e1b1fc3f1eb35b4bd5e465f0ef2bcb9997b36e3548637456b6");
    const R_MUL1: [u8; 32] =
        hex!("192eb5c3350a03e4baf23dd035b8804af8d5189c710adda53edd9cc0633f2d67");
    const R_SQR0: [u8; 32] =
        hex!("16bd622a907d7a92e475ed336e8ebca2cc1e38dd2ae69aaf2a96208eba0ee06e");
    const R_SQR1: [u8; 32] =
        hex!("8896d4306fb19d0e4a0e09899240e35cafed70bebb3ad56cf7b07964fefdfb93");
    const R_INV0: [u8; 32] =
        hex!("6face8b958e2bdc0771fd9d700f2703f881ef0d13509f16937f0a0c344647175");
    const R_INV1: [u8; 32] =
        hex!("93ceda7dddd537eb9307a06313598e650a568d931d16ab98ca0a7483c3b502e2");

    #[test]
    fn test_fq2_add() {
        // println!("test_fq2_add test");
        let mut a = Fq2::zero();
        for i in 1..LOOPS {
            a = a + Fq2::one();
            assert_eq!(
                *a.real(),
                Fq::new(i.into()).expect("u64 is a valid field element and static; qed")
            );
        }
        let x = Fq2::new(
            Fq::new(U256::from_slice(&X0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X1).unwrap()).unwrap(),
        );
        let y = Fq2::new(
            Fq::new(U256::from_slice(&Y0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&Y1).unwrap()).unwrap(),
        );
        let r = Fq2::new(
            Fq::new(U256::from_slice(&R_ADD0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_ADD1).unwrap()).unwrap(),
        );
        let re = x + y;
        assert_eq!(r, re);
    }

    #[test]
    fn test_fq2_mul() {
        // println!("test_fq2_mul test");
        let x = Fq2::new(
            Fq::new(U256::from_slice(&X0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X1).unwrap()).unwrap(),
        );
        let y = Fq2::new(
            Fq::new(U256::from_slice(&Y0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&Y1).unwrap()).unwrap(),
        );
        let r = Fq2::new(
            Fq::new(U256::from_slice(&R_MUL0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_MUL1).unwrap()).unwrap(),
        );
        let re = x * y;
        assert_eq!(r, re);
    }

    #[test]
    fn test_fq2_squared() {
        // println!("test_fq2_squared test");
        let x = Fq2::new(
            Fq::new(U256::from_slice(&X0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X1).unwrap()).unwrap(),
        );
        let r = Fq2::new(
            Fq::new(U256::from_slice(&R_SQR0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_SQR1).unwrap()).unwrap(),
        );
        let re = x.squared();
        assert_eq!(r, re);
        // i ^ 2  = -2
        assert_eq!(Fq2::i().squared(), (Fq2::one() + Fq2::one()).neg());
    }

    #[test]
    fn test_fq2_inverse() {
        // println!("test_fq2_inverse test");
        let x = Fq2::new(
            Fq::new(U256::from_slice(&X0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X1).unwrap()).unwrap(),
        );
        let r = Fq2::new(
            Fq::new(U256::from_slice(&R_INV0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_INV1).unwrap()).unwrap(),
        );
        let re = x.inverse().unwrap();
        assert_eq!(r, re);
    }

    #[test]
    fn test_fq2_sub() {
        // println!("test_fq2_sub test");
        let mut a = Fq2::new(
            Fq::new(LOOPS.into()).expect("u64 is a valid field element and static; qed"),
            Fq::zero(),
        );
        for i in (1..LOOPS).rev() {
            a = a - Fq2::one();
            assert_eq!(
                *a.real(),
                Fq::new(i.into()).expect("u64 is a valid field element and static; qed")
            );
        }
    }

    #[test]
    fn test_fq2_sqrt() {
        let a2 = Fq2::one() + Fq2::one();
        let a4 = a2 + a2;
        let a8 = a4 + a4;
        // i is sqrt(-2)
        assert_eq!(a2.neg().sqrt().unwrap(), Fq2::i());
        // 2i is sqrt(-8)
        assert_eq!(a8.neg().sqrt().unwrap(), Fq2::i() + Fq2::i());
        // 2 is sqrt(4)
        assert_eq!(a4.sqrt().unwrap(), a2);
        // 1 is sqrt(1)
        assert_eq!(Fq2::one().sqrt().unwrap(), Fq2::one());
    }

    const X00: [u8; 32] = hex!("3c2f6327ef1c5aa1d06e8cebc4100f0758c04476f40e8a0facb0a0bf09a9dd42");
    const X01: [u8; 32] = hex!("546e5945201b73c6ae44053114761efe351d5884c737301cfc7d2376d349a616");
    const X10: [u8; 32] = hex!("16bd622a907d7a92e475ed336e8ebca2cc1e38dd2ae69aaf2a96208eba0ee06e");
    const X11: [u8; 32] = hex!("8896d4306fb19d0e4a0e09899240e35cafed70bebb3ad56cf7b07964fefdfb93");

    const Y00: [u8; 32] = hex!("1b6ac9eb2c47b62cf61608b26c3c7e20674a48c4c509ac130bbaf6d47d32c07c");
    const Y01: [u8; 32] = hex!("0074a3145c65ac547541612178e584a902248740e70606dcaaafe2bcbd2f6a21");
    const Y10: [u8; 32] = hex!("8aed7a7f47f36b0f718cf99fcc59214c93ea0933c0583a7c5b61fca1962a6c5b");
    const Y11: [u8; 32] = hex!("45f1d11b8b8d1437342e2772863cb4c715a3fc4ee9d75a38904956428ec3c2c2");

    const R_MUL11: [u8; 32] =
        hex!("11d8f3dc2c4a7cd3ff4d557d86871210cff65187190711430b2d898affd61cda");
    const R_MUL10: [u8; 32] =
        hex!("960ee85c0aaacd6cc805053293a4955245ba973c9972b6767d0c68450a905ee7");
    const R_MUL01: [u8; 32] =
        hex!("ac9891b21d82827f6ccc2cd8524179b833239019c0b66cad89d7d8735ee03782");
    const R_MUL00: [u8; 32] =
        hex!("8f456b1cee442d189d01fc42fff7fd8481173dae8dc547d85c01a843005a063e");

    const R_MUL_FP11: [u8; 32] =
        hex!("413b76fe8748ab9130dc2907a55c15da925b496395c2cd82d6311863a4d9cfa8");
    const R_MUL_FP10: [u8; 32] =
        hex!("5cc754d5318f3ed489db7e53f94f3878a527053693983f4d4a61b30f6ea74984");
    const R_MUL_FP01: [u8; 32] =
        hex!("6769891769934201aa8d6de63cc012ec2b722d7b0ad9c9039246a3eea6f3d479");
    const R_MUL_FP00: [u8; 32] =
        hex!("408d33e58a4d3bfaf1d84a7ddad4e4026ca41f2aaa179611d9894584baed89d0");

    const R_MUL_FP211: [u8; 32] =
        hex!("242956015bdff53db568b970d64a7de56a0506309e1309b283317134dd52d53e");
    const R_MUL_FP210: [u8; 32] =
        hex!("5333c472d44677df131eeb1180badb3e1e9f88ba58190d16a92d95f939efb2c3");
    const R_MUL_FP201: [u8; 32] =
        hex!("0ccdaa76a6876ff69de6792161b614ca720bfcee2d5521533fbb28179ec0e31e");
    const R_MUL_FP200: [u8; 32] =
        hex!("2a2d6b832e919c313920f2e13e822795e2ceda8c0d8f4abe78220e4e00aeb6fd");

    const R_MUL_V11: [u8; 32] =
        hex!("ac9891b21d82827f6ccc2cd8524179b833239019c0b66cad89d7d8735ee03782");
    const R_MUL_V10: [u8; 32] =
        hex!("8f456b1cee442d189d01fc42fff7fd8481173dae8dc547d85c01a843005a063e");
    const R_MUL_V01: [u8; 32] =
        hex!("960ee85c0aaacd6cc805053293a4955245ba973c9972b6767d0c68450a905ee7");
    const R_MUL_V00: [u8; 32] =
        hex!("928e1847aa0ead49d7690054e880a3238205f03ce86ccc55cf148811e3a50bc9");

    const R_SQR11: [u8; 32] =
        hex!("8d3bc7848d4ad61017a7cb4efc280103bfe558e240c46c5765f1a4e2ec2e8c54");
    const R_SQR10: [u8; 32] =
        hex!("2f0f2ef9dd3979c7018b67837ba6e73938ba88ae66a101aaa0cf27ee449835ec");
    const R_SQR01: [u8; 32] =
        hex!("93838cbf9e5be34562c5bc031e27357d206f783837a6a921cbf4829292b69441");
    const R_SQR00: [u8; 32] =
        hex!("3681ecc58b68ffc15af31c5b1f1e10e1f3c60bdabb329c0dc7ffb2cc3925f005");

    const R_SQR_V11: [u8; 32] =
        hex!("93838cbf9e5be34562c5bc031e27357d206f783837a6a921cbf4829292b69441");
    const R_SQR_V10: [u8; 32] =
        hex!("3681ecc58b68ffc15af31c5b1f1e10e1f3c60bdabb329c0dc7ffb2cc3925f005");
    const R_SQR_V01: [u8; 32] =
        hex!("2f0f2ef9dd3979c7018b67837ba6e73938ba88ae66a101aaa0cf27ee449835ec");
    const R_SQR_V00: [u8; 32] =
        hex!("520870f6eab1a1c37cb7c001f2cd8c82c41a74d1b36d0508fefbec89ee457252");

    const R_INV11: [u8; 32] =
        hex!("1ec69309f84c5ad450750826fc804b72fb89fb48474222ba05be08bb1765f1d6");
    const R_INV10: [u8; 32] =
        hex!("3f16de331f77f510a3ec06e79319e3be5b3777471f79cd53404652b485133e99");
    const R_INV01: [u8; 32] =
        hex!("1cbf7f3bb04e2389184eade12de2752711cbff452363d2dfaf2bfef40618cebc");
    const R_INV00: [u8; 32] =
        hex!("3a70e829b83dc311970bc8d3e3e652f88a1ecd49b4672aa18c1c613c9a97d86f");

    #[test]
    fn test_fq4_mul() {
        // println!("test_fq4_mul test");
        let x0 = Fq2::new(
            Fq::new(U256::from_slice(&X00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X01).unwrap()).unwrap(),
        );
        let x1 = Fq2::new(
            Fq::new(U256::from_slice(&X10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X11).unwrap()).unwrap(),
        );
        let y0 = Fq2::new(
            Fq::new(U256::from_slice(&Y00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&Y01).unwrap()).unwrap(),
        );
        let y1 = Fq2::new(
            Fq::new(U256::from_slice(&Y10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&Y11).unwrap()).unwrap(),
        );
        let x = Fq4::new(x0, x1);
        let y = Fq4::new(y0, y1);
        let r0 = Fq2::new(
            Fq::new(U256::from_slice(&R_MUL00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_MUL01).unwrap()).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::new(U256::from_slice(&R_MUL10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_MUL11).unwrap()).unwrap(),
        );
        let r = Fq4::new(r0, r1);
        assert_eq!(r, x * y);
    }

    #[test]
    fn test_fq4_mul_fq() {
        // println!("test_fq4_mul_fq test");
        let x0 = Fq2::new(
            Fq::new(U256::from_slice(&X00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X01).unwrap()).unwrap(),
        );
        let x1 = Fq2::new(
            Fq::new(U256::from_slice(&X10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X11).unwrap()).unwrap(),
        );
        let k = Fq::from_slice(&HEX_IV).unwrap();
        let x = Fq4::new(x0, x1);
        let r0 = Fq2::new(
            Fq::new(U256::from_slice(&R_MUL_FP00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_MUL_FP01).unwrap()).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::new(U256::from_slice(&R_MUL_FP10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_MUL_FP11).unwrap()).unwrap(),
        );
        let r = Fq4::new(r0, r1);
        assert_eq!(r, x.scale_fq(&k));
    }

    #[test]
    fn test_fq4_mul_fq2() {
        // println!("test_fq4_mul_fq2 test");
        let x0 = Fq2::new(
            Fq::new(U256::from_slice(&X00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X01).unwrap()).unwrap(),
        );
        let x1 = Fq2::new(
            Fq::new(U256::from_slice(&X10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X11).unwrap()).unwrap(),
        );
        let y = Fq2::new(
            Fq::new(U256::from_slice(&Y0).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&Y1).unwrap()).unwrap(),
        );
        let x = Fq4::new(x0, x1);
        let r0 = Fq2::new(
            Fq::new(U256::from_slice(&R_MUL_FP200).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_MUL_FP201).unwrap()).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::new(U256::from_slice(&R_MUL_FP210).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_MUL_FP211).unwrap()).unwrap(),
        );
        let r = Fq4::new(r0, r1);
        assert_eq!(r, x.scale(&y));
    }

    #[test]
    fn test_fq4_inv() {
        // println!("test_fq4_inv test");
        can_invert::<Fq4>();
        let x0 = Fq2::new(
            Fq::new(U256::from_slice(&X00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X01).unwrap()).unwrap(),
        );
        let x1 = Fq2::new(
            Fq::new(U256::from_slice(&X10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X11).unwrap()).unwrap(),
        );
        let x = Fq4::new(x0, x1);
        let r0 = Fq2::new(
            Fq::new(U256::from_slice(&R_INV00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_INV01).unwrap()).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::new(U256::from_slice(&R_INV10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_INV11).unwrap()).unwrap(),
        );
        let r = Fq4::new(r0, r1);
        assert_eq!(r, x.inverse().unwrap());
    }

    #[test]
    fn test_fq4_sqr() {
        // println!("test_fq4_sqr test");
        let x0 = Fq2::new(
            Fq::new(U256::from_slice(&X00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X01).unwrap()).unwrap(),
        );
        let x1 = Fq2::new(
            Fq::new(U256::from_slice(&X10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&X11).unwrap()).unwrap(),
        );
        let x = Fq4::new(x0, x1);
        let r0 = Fq2::new(
            Fq::new(U256::from_slice(&R_SQR00).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_SQR01).unwrap()).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::new(U256::from_slice(&R_SQR10).unwrap()).unwrap(),
            Fq::new(U256::from_slice(&R_SQR11).unwrap()).unwrap(),
        );
        let r = Fq4::new(r0, r1);
        assert_eq!(r, x.squared());
    }

    const FP12_MULC11: [u8; 32] =
        hex!("058d43459faee14ba2b6a69ff2d8c3ad933a1253e1764dedf5419b144a2ab82b");
    const FP12_MULC10: [u8; 32] =
        hex!("20ef84805ba02ef92a48fb2ae8086e566a644ab0639249f175268f18d8091ad4");
    const FP12_MULC01: [u8; 32] =
        hex!("83cc3be54a699ae24d8f920c87baa395befb424a6dcad1dcdfc2a006765ef8d5");
    const FP12_MULC00: [u8; 32] =
        hex!("1d705169165d9c2386c3bc673df3fa84975afa955a7be27f1b362000a96b8c2c");

    const FP12_MULB11: [u8; 32] =
        hex!("22b910d826f02961ff0fed439beb1e91f45193f87c2cdd9562da539290846ace");
    const FP12_MULB10: [u8; 32] =
        hex!("2c618991ae82d35063cfed629ff7d930b8070ba07d0652ba092f046e133e3491");
    const FP12_MULB01: [u8; 32] =
        hex!("137bc78a9aa182330bd71fb8859314422dd36f5e3c1f6fd36d6c9685fc39419f");
    const FP12_MULB00: [u8; 32] =
        hex!("8d83e7380abe10a2f3677864c2dbbcdad7ae5434e92043a2da3b71f3f9cedd8c");

    const FP12_MULA11: [u8; 32] =
        hex!("850c0562ac08996c05d22ea466cf4b1fa7a7064d4653b5fa725d623254bf7125");
    const FP12_MULA10: [u8; 32] =
        hex!("6dc41016b3ab9b44a4841aa8037e3b4d331cc7c8313abee0c5111a9be5915e90");
    const FP12_MULA01: [u8; 32] =
        hex!("6d1a15e5b765c4b139bf5c6c4a87214c269b26fb709ff5de885c053f405cf626");
    const FP12_MULA00: [u8; 32] =
        hex!("8d4d853489a4a5d809fa77e35627a5351651b926f001e1ee46e95808f9001d24");
    //

    const FP12_SQRC11: [u8; 32] =
        hex!("3592cba3482fb39756b2ed1d3d756685caa005bd5e8288bc92841d29276aa321");
    const FP12_SQRC10: [u8; 32] =
        hex!("8e3a49919e6de83b1ab1a5bb9eb993c3bbd68e8d305aed5c0b88cef0ef41c47f");
    const FP12_SQRC01: [u8; 32] =
        hex!("3d3d9cc8e07619efd21745f6938a26f7cb0a83ad4aa3a9d066e18ad99833e3ac");
    const FP12_SQRC00: [u8; 32] =
        hex!("25195ec7af551c42d7d37a0b120607d4adba6b9377299688b92a8393f3b8c20f");

    const FP12_SQRB11: [u8; 32] =
        hex!("76f676d5d2cb8d1a2cc237fc78c8d544bef1cd560e654236f502aed0d8c9148c");
    const FP12_SQRB10: [u8; 32] =
        hex!("6cde174a5e9d117175a4a163f041b65f868dffa05b5f3474f729b87f92493f2c");
    const FP12_SQRB01: [u8; 32] =
        hex!("667a86d73e8f88a81306f7f0cd28789a55bf7e9cbe155fc6abb300ad027d8801");
    const FP12_SQRB00: [u8; 32] =
        hex!("a49a66d48ec2ef72a9929413a40e316a8aee1d6236a1db8c56496524f1c23f11");

    const FP12_SQRA11: [u8; 32] =
        hex!("1684bc9679aaba4afe35ec8c0852e438f41e15ab37620d9661018f90fe7415f1");
    const FP12_SQRA10: [u8; 32] =
        hex!("8d37fb8b7edf942885b3009cf7e295bea89444d34091fc57380c778395b7c4e4");
    const FP12_SQRA01: [u8; 32] =
        hex!("278b9d9ea61b6b2758e758ed9a64034576b520e65a9d276a0c82f079501a226e");
    const FP12_SQRA00: [u8; 32] =
        hex!("01a333fa4177601de7cd8ed49ea4906f30e23988dcb7cde173da48499fce3ee5");
    //

    const FP12_INVC11: [u8; 32] =
        hex!("47ae900b90945e31afde7fe09f0b69640c468a1648ee52070584a5d13af22bb9");
    const FP12_INVC10: [u8; 32] =
        hex!("8f273655182c3a9f184dc30421161ecdd50655c36a9266c7df1016e410f34102");
    const FP12_INVC01: [u8; 32] =
        hex!("a26e789013203804b5f8f1c5a51dd3fb50176d41108b235d6e66712721060252");
    const FP12_INVC00: [u8; 32] =
        hex!("090aaed5cb83068a0376c6eaca210007744d00c8b4ce53279a67cc069cc519e7");

    const FP12_INVB11: [u8; 32] =
        hex!("80ab89aa446df59ffe2f29cdb917b760d740ceb634c731b93bf1661aa5868b54");
    const FP12_INVB10: [u8; 32] =
        hex!("1e13ab51b3198619cc0016599562ed4d266d1481d0d273d3f97cffe5f8e0dd21");
    const FP12_INVB01: [u8; 32] =
        hex!("5aeb8ed89aafc971a857b8d02f3e3c37ef15ba0e3220e3a7c13c9da8af0c393b");
    const FP12_INVB00: [u8; 32] =
        hex!("518c338b1430e3129c2555650e5d5634d89513f694ba3a5f2aeb444c540f125a");

    const FP12_INVA11: [u8; 32] =
        hex!("aba8c5682695f3feee64772d0e49b432c96470e7d663098e9c271a91d4fc991a");
    const FP12_INVA10: [u8; 32] =
        hex!("0ed800dabe29af5fb41a41cc49fd4084deb02442e8e66f88186607f46395e533");
    const FP12_INVA01: [u8; 32] =
        hex!("a31b642cd5453c7bb16c82bc67bd3b66fa4db58b8e9aa45f9b579860f18d402c");
    const FP12_INVA00: [u8; 32] =
        hex!("798b84002e95753e3b07027a8d68b0a7ab2ac40328fc7ca3ea40780b3428dbc1");

    #[test]
    fn test_fq12_sqr() {
        // println!("test_fq12_sqr test");
        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL00).unwrap(),
            Fq::from_slice(&R_MUL01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL10).unwrap(),
            Fq::from_slice(&R_MUL11).unwrap(),
        );
        let a = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL_FP00).unwrap(),
            Fq::from_slice(&R_MUL_FP01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL_FP10).unwrap(),
            Fq::from_slice(&R_MUL_FP11).unwrap(),
        );
        let b = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL_FP200).unwrap(),
            Fq::from_slice(&R_MUL_FP201).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL_FP210).unwrap(),
            Fq::from_slice(&R_MUL_FP211).unwrap(),
        );
        let c = Fq4::new(r0, r1);
        let x = Fq12::new(a, b, c);

        let r0 = Fq2::new(
            Fq::from_slice(&FP12_SQRA00).unwrap(),
            Fq::from_slice(&FP12_SQRA01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&FP12_SQRA10).unwrap(),
            Fq::from_slice(&FP12_SQRA11).unwrap(),
        );
        let a = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&FP12_SQRB00).unwrap(),
            Fq::from_slice(&FP12_SQRB01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&FP12_SQRB10).unwrap(),
            Fq::from_slice(&FP12_SQRB11).unwrap(),
        );
        let b = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&FP12_SQRC00).unwrap(),
            Fq::from_slice(&FP12_SQRC01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&FP12_SQRC10).unwrap(),
            Fq::from_slice(&FP12_SQRC11).unwrap(),
        );
        let c = Fq4::new(r0, r1);
        let r = Fq12::new(a, b, c);

        assert_eq!(r, x.squared());
    }

    #[test]
    fn test_fq12_mul() {
        // println!("test_fq12_mul test");
        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL00).unwrap(),
            Fq::from_slice(&R_MUL01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL10).unwrap(),
            Fq::from_slice(&R_MUL11).unwrap(),
        );
        let a = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL_FP00).unwrap(),
            Fq::from_slice(&R_MUL_FP01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL_FP10).unwrap(),
            Fq::from_slice(&R_MUL_FP11).unwrap(),
        );
        let b = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL_FP200).unwrap(),
            Fq::from_slice(&R_MUL_FP201).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL_FP210).unwrap(),
            Fq::from_slice(&R_MUL_FP211).unwrap(),
        );
        let c = Fq4::new(r0, r1);
        let x = Fq12::new(a, b, c);

        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL_V00).unwrap(),
            Fq::from_slice(&R_MUL_V01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL_V10).unwrap(),
            Fq::from_slice(&R_MUL_V11).unwrap(),
        );
        let a = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&R_SQR00).unwrap(),
            Fq::from_slice(&R_SQR01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_SQR10).unwrap(),
            Fq::from_slice(&R_SQR11).unwrap(),
        );
        let b = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&R_INV00).unwrap(),
            Fq::from_slice(&R_INV01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_INV10).unwrap(),
            Fq::from_slice(&R_INV11).unwrap(),
        );
        let c = Fq4::new(r0, r1);
        let y = Fq12::new(a, b, c);

        let r0 = Fq2::new(
            Fq::from_slice(&FP12_MULA00).unwrap(),
            Fq::from_slice(&FP12_MULA01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&FP12_MULA10).unwrap(),
            Fq::from_slice(&FP12_MULA11).unwrap(),
        );
        let a = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&FP12_MULB00).unwrap(),
            Fq::from_slice(&FP12_MULB01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&FP12_MULB10).unwrap(),
            Fq::from_slice(&FP12_MULB11).unwrap(),
        );
        let b = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&FP12_MULC00).unwrap(),
            Fq::from_slice(&FP12_MULC01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&FP12_MULC10).unwrap(),
            Fq::from_slice(&FP12_MULC11).unwrap(),
        );
        let c = Fq4::new(r0, r1);
        let r = Fq12::new(a, b, c);

        assert_eq!(r, x * y);
    }

    #[test]
    fn test_fq12_inv() {
        // println!("test_fq12_inv test");
        can_invert::<Fq12>();
        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL00).unwrap(),
            Fq::from_slice(&R_MUL01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL10).unwrap(),
            Fq::from_slice(&R_MUL11).unwrap(),
        );
        let a = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL_FP00).unwrap(),
            Fq::from_slice(&R_MUL_FP01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL_FP10).unwrap(),
            Fq::from_slice(&R_MUL_FP11).unwrap(),
        );
        let b = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&R_MUL_FP200).unwrap(),
            Fq::from_slice(&R_MUL_FP201).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&R_MUL_FP210).unwrap(),
            Fq::from_slice(&R_MUL_FP211).unwrap(),
        );
        let c = Fq4::new(r0, r1);
        let x = Fq12::new(a, b, c);

        let r0 = Fq2::new(
            Fq::from_slice(&FP12_INVA00).unwrap(),
            Fq::from_slice(&FP12_INVA01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&FP12_INVA10).unwrap(),
            Fq::from_slice(&FP12_INVA11).unwrap(),
        );
        let a = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&FP12_INVB00).unwrap(),
            Fq::from_slice(&FP12_INVB01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&FP12_INVB10).unwrap(),
            Fq::from_slice(&FP12_INVB11).unwrap(),
        );
        let b = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&FP12_INVC00).unwrap(),
            Fq::from_slice(&FP12_INVC01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&FP12_INVC10).unwrap(),
            Fq::from_slice(&FP12_INVC11).unwrap(),
        );
        let c = Fq4::new(r0, r1);
        let r = Fq12::new(a, b, c);

        assert_eq!(r, x.inverse().unwrap());
    }
}
