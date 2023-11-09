#![allow(dead_code)]

use core::{
    fmt,
    ops::{Add, Mul, Neg, Sub},
};
use hex_literal::hex;

use crate::fields::{FieldElement, Fq, Fq12, Fq2, Fq4, Fr};
use crate::u256::U256;

// abits = "00100000000000000000000000000000000000010000101100020200101000020";
const SM9_LOOP_COUNT: [u8; 65] = [
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 2, 0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2,
    0,
];

const SM9_P1X: [u8; 32] =
    hex!("93DE051D 62BF718F F5ED0704 487D01D6 E1E40869 09DC3280 E8C4E481 7C66DDDD");
const SM9_P1Y: [u8; 32] =
    hex!("21FE8DDA 4F21E607 63106512 5C395BBC 1C1C00CB FA602435 0C464CD7 0A3EA616");
const SM9_P2X1: [u8; 32] =
    hex!("85AEF3D0 78640C98 597B6027 B441A01F F1DD2C19 0F5E93C4 54806C11 D8806141");
const SM9_P2X0: [u8; 32] =
    hex!("37227552 92130B08 D2AAB97F D34EC120 EE265948 D19C17AB F9B7213B AF82D65B");
const SM9_P2Y1: [u8; 32] =
    hex!("17509B09 2E845C12 66BA0D26 2CBEE6ED 0736A96F A347C8BD 856DC76B 84EBEB96");
const SM9_P2Y0: [u8; 32] =
    hex!("A7CF28D5 19BE3DA6 5F317015 3D278FF2 47EFBA98 A71A0811 6215BBA5 C999A7C7");

lazy_static::lazy_static! {

    // PI1 = 0x3f23ea58e5720bdb 843c6cfa9c086749 47c5c86e0ddd04ed a91d8354377b698b
    // PI2 = 0xf300000002a3a6f2 780272354f8b78f4 d5fc11967be65334
    static ref SM9_PI1: U256 = U256::from([
        0xa91d8354377b698b,
        0x47c5c86e0ddd04ed,
        0x843c6cfa9c086749,
        0x3f23ea58e5720bdb
        ]);
    static ref SM9_PI2: U256 = U256::from([
        0xd5fc11967be65334,
        0x780272354f8b78f4,
        0xf300000002a3a6f2,
        0x0
    ]);
}

pub trait GroupElement:
    Sized
    + Copy
    + Clone
    + PartialEq
    + Eq
    + fmt::Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Mul<Fr, Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn double(&self) -> Self;
}

pub trait GroupParams: Sized + fmt::Debug {
    type Base: FieldElement;

    fn name() -> &'static str;
    fn one() -> G<Self>;
    fn coeff_b() -> Self::Base;
    fn check_order() -> bool {
        false
    }
}

#[repr(C)]
pub struct G<P: GroupParams> {
    x: P::Base,
    y: P::Base,
    z: P::Base,
}

impl<P: GroupParams> G<P> {
    pub fn new(x: P::Base, y: P::Base, z: P::Base) -> Self {
        G { x: x, y: y, z: z }
    }

    pub fn x(&self) -> &P::Base {
        &self.x
    }

    pub fn x_mut(&mut self) -> &mut P::Base {
        &mut self.x
    }

    pub fn y(&self) -> &P::Base {
        &self.y
    }

    pub fn y_mut(&mut self) -> &mut P::Base {
        &mut self.y
    }

    pub fn z(&self) -> &P::Base {
        &self.z
    }

    pub fn z_mut(&mut self) -> &mut P::Base {
        &mut self.z
    }
}

#[derive(Debug)]
pub struct AffineG<P: GroupParams> {
    x: P::Base,
    y: P::Base,
}

#[derive(Debug)]
pub enum Error {
    NotOnCurve,
    NotInSubgroup,
}

impl<P: GroupParams> AffineG<P> {
    pub fn new(x: P::Base, y: P::Base) -> Result<Self, Error> {
        if y.squared() == (x.squared() * x) + P::coeff_b() {
            if P::check_order() {
                let p: G<P> = G {
                    x: x,
                    y: y,
                    z: P::Base::one(),
                };

                if (p * (-Fr::one())) + p != G::zero() {
                    return Err(Error::NotInSubgroup);
                }
            }

            Ok(AffineG { x: x, y: y })
        } else {
            Err(Error::NotOnCurve)
        }
    }

    pub fn x(&self) -> &P::Base {
        &self.x
    }

    pub fn x_mut(&mut self) -> &mut P::Base {
        &mut self.x
    }

    pub fn y(&self) -> &P::Base {
        &self.y
    }

    pub fn y_mut(&mut self) -> &mut P::Base {
        &mut self.y
    }
}

impl<P: GroupParams> PartialEq for AffineG<P> {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x && self.y == other.y
    }
}

impl<P: GroupParams> Eq for AffineG<P> {}

impl<P: GroupParams> fmt::Debug for G<P> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}({:?}, {:?}, {:?})", P::name(), self.x, self.y, self.z)
    }
}

impl<P: GroupParams> Clone for G<P> {
    fn clone(&self) -> Self {
        G {
            x: self.x,
            y: self.y,
            z: self.z,
        }
    }
}

impl<P: GroupParams> Copy for G<P> {}

impl<P: GroupParams> Clone for AffineG<P> {
    fn clone(&self) -> Self {
        AffineG {
            x: self.x,
            y: self.y,
        }
    }
}

impl<P: GroupParams> Copy for AffineG<P> {}

impl<P: GroupParams> PartialEq for G<P> {
    fn eq(&self, other: &Self) -> bool {
        if self.is_zero() {
            return other.is_zero();
        }

        if other.is_zero() {
            return false;
        }

        let z1_squared = self.z.squared();
        let z2_squared = other.z.squared();

        if self.x * z2_squared != other.x * z1_squared {
            return false;
        }

        let z1_cubed = self.z * z1_squared;
        let z2_cubed = other.z * z2_squared;

        if self.y * z2_cubed != other.y * z1_cubed {
            return false;
        }

        return true;
    }
}
impl<P: GroupParams> Eq for G<P> {}

impl<P: GroupParams> G<P> {
    pub fn to_affine(&self) -> Option<AffineG<P>> {
        if self.z.is_zero() {
            None
        } else if self.z == P::Base::one() {
            Some(AffineG {
                x: self.x,
                y: self.y,
            })
        } else {
            let zinv = self.z.inverse().unwrap();
            let zinv_squared = zinv.squared();

            Some(AffineG {
                x: self.x * zinv_squared,
                y: self.y * (zinv_squared * zinv),
            })
        }
    }
}

impl<P: GroupParams> AffineG<P> {
    pub fn to_jacobian(&self) -> G<P> {
        G {
            x: self.x,
            y: self.y,
            z: P::Base::one(),
        }
    }
}

impl<P: GroupParams> GroupElement for G<P> {
    fn zero() -> Self {
        G {
            x: P::Base::zero(),
            y: P::Base::one(),
            z: P::Base::zero(),
        }
    }

    fn one() -> Self {
        P::one()
    }

    fn is_zero(&self) -> bool {
        self.z.is_zero()
    }

    fn double(&self) -> Self {
        let a = self.x.squared();
        let b = self.y.squared();
        let c = b.squared();
        let mut d = (self.x + b).squared() - a - c;
        d = d + d;
        let e = a + a + a;
        let f = e.squared();
        let x3 = f - (d + d);
        let mut eight_c = c + c;
        eight_c = eight_c + eight_c;
        eight_c = eight_c + eight_c;
        let y1z1 = self.y * self.z;

        G {
            x: x3,
            y: e * (d - x3) - eight_c,
            z: y1z1 + y1z1,
        }
    }
}

impl<P: GroupParams> Mul<Fr> for G<P> {
    type Output = G<P>;

    fn mul(self, other: Fr) -> G<P> {
        let mut res = G::zero();
        let mut found_one = false;

        for i in U256::from(other).bits() {
            if found_one {
                res = res.double();
            }

            if i {
                found_one = true;
                res = res + self;
            }
        }

        res
    }
}

impl<P: GroupParams> Add<G<P>> for G<P> {
    type Output = G<P>;

    fn add(self, other: G<P>) -> G<P> {
        if self.is_zero() {
            return other;
        }

        if other.is_zero() {
            return self;
        }

        let z1_squared = self.z.squared();
        let z2_squared = other.z.squared();
        let u1 = self.x * z2_squared;
        let u2 = other.x * z1_squared;
        let z1_cubed = self.z * z1_squared;
        let z2_cubed = other.z * z2_squared;
        let s1 = self.y * z2_cubed;
        let s2 = other.y * z1_cubed;

        if u1 == u2 && s1 == s2 {
            self.double()
        } else {
            let h = u2 - u1;
            let s2_minus_s1 = s2 - s1;
            let i = (h + h).squared();
            let j = h * i;
            let r = s2_minus_s1 + s2_minus_s1;
            let v = u1 * i;
            let s1_j = s1 * j;
            let x3 = r.squared() - j - (v + v);

            G {
                x: x3,
                y: r * (v - x3) - (s1_j + s1_j),
                z: ((self.z + other.z).squared() - z1_squared - z2_squared) * h,
            }
        }
    }
}

impl<P: GroupParams> Neg for G<P> {
    type Output = G<P>;

    fn neg(self) -> G<P> {
        if self.is_zero() {
            self
        } else {
            G {
                x: self.x,
                y: -self.y,
                z: self.z,
            }
        }
    }
}

impl<P: GroupParams> Neg for AffineG<P> {
    type Output = AffineG<P>;

    fn neg(self) -> AffineG<P> {
        AffineG {
            x: self.x,
            y: -self.y,
        }
    }
}

impl<P: GroupParams> Sub<G<P>> for G<P> {
    type Output = G<P>;

    fn sub(self, other: G<P>) -> G<P> {
        self + (-other)
    }
}

#[derive(Debug)]
pub struct G1Params;

impl GroupParams for G1Params {
    type Base = Fq;

    fn name() -> &'static str {
        "G1"
    }
    //generator
    fn one() -> G<Self> {
        G {
            x: Fq::from_hex(&SM9_P1X).unwrap(),
            y: Fq::from_hex(&SM9_P1Y).unwrap(),
            z: Fq::one(),
        }
    }

    fn coeff_b() -> Fq {
        Fq::from_str("5").unwrap()
    }
}

pub type G1 = G<G1Params>;

pub type AffineG1 = AffineG<G1Params>;

#[derive(Debug)]
pub struct G2Params;

impl GroupParams for G2Params {
    type Base = Fq2;

    fn name() -> &'static str {
        "G2"
    }
    // generator
    fn one() -> G<Self> {
        G {
            x: Fq2::new(
                Fq::from_hex(&SM9_P2X0).unwrap(),
                Fq::from_hex(&SM9_P2X1).unwrap(),
            ),
            y: Fq2::new(
                Fq::from_hex(&SM9_P2Y0).unwrap(),
                Fq::from_hex(&SM9_P2Y1).unwrap(),
            ),
            z: Fq2::one(),
        }
    }

    fn coeff_b() -> Fq2 {
        Fq2::one().scale(Fq::from_str("5").unwrap())
    }

    fn check_order() -> bool {
        true
    }
}

pub type G2 = G<G2Params>;

pub type AffineG2 = AffineG<G2Params>;

impl G2 {
    fn add_full(&self, q: &G2) -> Self {
        if self.is_zero() {
            return *q;
        }

        if q.is_zero() {
            return *self;
        }

        let t1 = self.z.squared();
        let t2 = q.z.squared();
        let t3 = q.x * t1;
        let t4 = self.x * t2;
        let t5 = t3 + t4;
        let t3 = t3 - t4;
        let t1 = t1 * self.z;
        let t1 = t1 * q.y;
        let t2 = t2 * q.z;
        let t2 = t2 * self.y;
        let t6 = t1 + t2;
        let t1 = t1 - t2;

        if t1.is_zero() {
            if t3.is_zero() {
                return self.double();
            }
            if t6.is_zero() {
                return Self::zero();
            }
        }

        let t6 = t1.squared();
        let t7 = t3 * self.z;
        let t7 = t7 * q.z;
        let t8 = t3.squared();
        let t5 = t5 * t8;
        let t3 = t3 * t8;
        let t4 = t4 * t8;
        let t6 = t6 - t5;
        let t4 = t4 - t6;
        let t1 = t1 * t4;
        let t2 = t2 * t3;
        let t1 = t1 - t2;

        G2 {
            x: t6,
            y: t1,
            z: t7,
        }
    }

    fn point_pi1(&self) -> Self {
        G2 {
            x: self.x.frobenius_map(1),
            y: self.y.frobenius_map(1),
            z: self.z.frobenius_map(1).scale(Fq::new(*SM9_PI1).unwrap()),
        }
    }

    fn point_pi2(&self) -> Self {
        G2 {
            x: self.x,
            y: self.y,
            z: self.z.scale(Fq::new(*SM9_PI2).unwrap()),
        }
    }

    fn eval_g_tangent(&self, q: &G1) -> (Fq12, Fq12) {
        let mut num = Fq12::zero();
        let mut den = Fq12::zero();

        let t0 = self.z.squared();
        let t1 = t0 * self.z;
        let b1 = t1 * self.y;

        let t2 = b1.scale(q.y);
        let a1 = -t2;

        let t1 = self.x.squared();
        let t0 = t0 * t1;
        let t0 = t0.scale(q.x);
        let t0 = t0.tri();
        let a4 = t0.div2();

        let t1 = t1 * self.x;
        let t1 = t1.tri();
        let t1 = t1.div2();
        let t0 = self.y.squared();
        let a0 = t0 - t1;

        num.c0 = Fq4::new(a0, a1);
        num.c2 = Fq4::new(a4, Fq2::zero());
        den.c0 = Fq4::new(Fq2::zero(), b1);

        (num, den)
    }

    fn eval_g_line(&self, p: &G2, q: &G1) -> (Fq12, Fq12) {
        let mut num = Fq12::zero();
        let mut den = Fq12::zero();

        let t0 = p.z.squared();
        let t1 = t0 * self.x;
        let t0 = t0 * p.z;

        let t2 = self.z.squared();
        let t3 = t2 * p.x;
        let t2 = t2 * self.z;

        let t2 = t2 * p.y;
        let t1 = t1 - t3;
        let t1 = t1 * self.z;

        let t1 = t1 * p.z;
        let t4 = t1 * t0;
        let b1 = t4;

        let t1 = t1 * p.y;
        let t3 = t0 * self.y;
        let t3 = t3 - t2;
        let t0 = t0 * t3;
        let t0 = t0.scale(q.x);
        let a4 = t0;

        let t3 = t3 * p.x;
        let t3 = t3 * p.z;
        let t1 = t1 - t3;
        let a0 = t1;

        let t2 = t4.scale(q.y);
        let a1 = -t2;

        num.c0 = Fq4::new(a0, a1);
        num.c2 = Fq4::new(a4, Fq2::zero());
        den.c0 = Fq4::new(Fq2::zero(), b1);

        (num, den)
    }

    pub fn miller_loop(&self, p: &G1) -> Fq12 {
        let mut t = *self;
        let mut f_num = Fq12::one();
        let mut f_den = Fq12::one();
        let mut g_num;
        let mut g_den;
        let mut q1 = self.neg();

        for i in SM9_LOOP_COUNT.iter() {
            f_num = f_num.squared();
            f_den = f_den.squared();
            (g_num, g_den) = t.eval_g_tangent(p);
            f_num = f_num * g_num;
            f_den = f_den * g_den;

            t = t.double();

            if *i == 1 {
                (g_num, g_den) = t.eval_g_line(self, p);
                f_num = f_num * g_num;
                f_den = f_den * g_den;
                t = t.add_full(self);
            } else if *i == 2 {
                (g_num, g_den) = t.eval_g_line(&q1, p);
                f_num = f_num * g_num;
                f_den = f_den * g_den;
                t = t.add_full(&q1);
            }
        }
        q1 = self.point_pi1();
        let q2 = -self.point_pi2();

        (g_num, g_den) = t.eval_g_line(&q1, p);
        f_num = f_num * g_num;
        f_den = f_den * g_den;
        t = t.add_full(&q1);

        (g_num, g_den) = t.eval_g_line(&q2, p);
        f_num = f_num * g_num;
        f_den = f_den * g_den;
        //    t = t + q2;

        f_den = f_den.inverse().unwrap();
        f_num * f_den
    }
}

//
// R-ate Pairing G2 x G1 -> GT
//
// P is a point of order r in G1. Q(x,y) is a point of order r in G2.
// Note that Q is a point on the sextic twist of the curve over Fp^2, P(x,y) is a point on the
// curve over the base field Fp
// the curve is y^2 = x^3 + 5
pub fn pairing(p: &G1, q: &G2) -> Fq12 {
    q.miller_loop(p)
        .final_exponentiation()
        .expect("miller loop cannot produce zero")
}

/* ************************************************************************************************ */
#[cfg(test)]
mod tests {
    use super::*;

    pub fn group_trials<G: GroupElement>() {
        assert!(G::zero().is_zero());
        assert!((G::one() - G::one()).is_zero());
        assert_eq!(G::one() + G::one(), G::one() * Fr::from_str("2").unwrap());
        assert!(G::zero().double().is_zero());

        assert!((G::one() * (-Fr::one()) + G::one()).is_zero());
    }

    #[test]
    fn test_g1() {
        // println!("test_G1 test");
        group_trials::<G1>();
    }

    const P1X: [u8; 32] = hex!("917be49d159184fba140f4dfc5d653464e94f718fe195b226b3f715829e6e768");
    const P1Y: [u8; 32] = hex!("288578d9505d462867a50acee40ee143b896e72505be10e8ce4c6b0c945b642b");
    const P2X: [u8; 32] = hex!("593417680f252445fd0522383e23c77a54b11fe222de4a886eabc26e16bffa3c");
    const P2Y: [u8; 32] = hex!("38e8fc9a8b60f5ba0c6c411f721c117044435a833757d8fee65828511b8b245d");

    const P_DBL_X: [u8; 32] =
        hex!("268def7968f1e8c51635e277425403df88355fb2ecf16f7920f112eb2a7e50c9");
    const P_DBL_Y: [u8; 32] =
        hex!("5c596b534bbaa85c1d3aecf436e61ff1bfd9f70856f0309c2a63d8248205d84e");
    const P_ADD_X: [u8; 32] =
        hex!("056610cb69f8d5659ea94e4a67bbf3b93fb0bd449672d7ca2525ec3b68c894d1");
    const P_ADD_Y: [u8; 32] =
        hex!("88f3f99ce78ed3ffe6ca1cface5242570cb5d053f16a8e0baae10414babd86a7");
    const P_NEG_X: [u8; 32] =
        hex!("917be49d159184fba140f4dfc5d653464e94f718fe195b226b3f715829e6e768");
    const P_NEG_Y: [u8; 32] =
        hex!("8dba8726b24660c96e5ea081117fe601695bac2614bcddf31723301b4ef5e152");
    const P_SUB_X: [u8; 32] =
        hex!("29e4a54cad98da9939b95f677784bff3b1dd9334c83d93e351e0f8f7c4ce2dc5");
    const P_SUB_Y: [u8; 32] =
        hex!("4473eba3b8ff990b8456c41ec0727b76cb2b0f960495b144949f70bf95643b82");
    const P_MUL_X: [u8; 32] =
        hex!("997fcff625adbae62566f684f9e89181713f972c5a9cd9ce6764636761ba87d1");
    const P_MUL_Y: [u8; 32] =
        hex!("8142a28d1bd109501452a649e2d68f012e265460e0c7d3da743fb036eb23b03b");

    const HEX_IV: [u8; 32] =
        hex!("123456789abcdef00fedcba987654321123456789abcdef00fedcba987654321");

    #[test]
    fn test_g1_dbl() {
        // println!("test_G1_dbl test");
        let p1 = G1::new(
            Fq::from_hex(&P1X).unwrap(),
            Fq::from_hex(&P1Y).unwrap(),
            Fq::one(),
        );
        let r = G1::new(
            Fq::from_hex(&P_DBL_X).unwrap(),
            Fq::from_hex(&P_DBL_Y).unwrap(),
            Fq::one(),
        );
        assert_eq!(r, p1.double());
    }

    #[test]
    fn test_g1_add() {
        // println!("test_G1_add test");
        let p1 = G1::new(
            Fq::from_hex(&P1X).unwrap(),
            Fq::from_hex(&P1Y).unwrap(),
            Fq::one(),
        );
        let p2 = G1::new(
            Fq::from_hex(&P2X).unwrap(),
            Fq::from_hex(&P2Y).unwrap(),
            Fq::one(),
        );
        let r = G1::new(
            Fq::from_hex(&P_ADD_X).unwrap(),
            Fq::from_hex(&P_ADD_Y).unwrap(),
            Fq::one(),
        );
        assert_eq!(r, p1 + p2);
    }

    #[test]
    fn test_g1_sub() {
        // println!("test_G1_sub test");
        let p1 = G1::new(
            Fq::from_hex(&P1X).unwrap(),
            Fq::from_hex(&P1Y).unwrap(),
            Fq::one(),
        );
        let p2 = G1::new(
            Fq::from_hex(&P2X).unwrap(),
            Fq::from_hex(&P2Y).unwrap(),
            Fq::one(),
        );
        let r = G1::new(
            Fq::from_hex(&P_NEG_X).unwrap(),
            Fq::from_hex(&P_NEG_Y).unwrap(),
            Fq::one(),
        );
        assert_eq!(r, -p1);
        let r = G1::new(
            Fq::from_hex(&P_SUB_X).unwrap(),
            Fq::from_hex(&P_SUB_Y).unwrap(),
            Fq::one(),
        );
        assert_eq!(r, p1 - p2);
    }

    #[test]
    fn test_g1_mul() {
        // println!("test_G1_mul test");
        let p1 = G1::new(
            Fq::from_hex(&P1X).unwrap(),
            Fq::from_hex(&P1Y).unwrap(),
            Fq::one(),
        );
        let r = G1::new(
            Fq::from_hex(&P_MUL_X).unwrap(),
            Fq::from_hex(&P_MUL_Y).unwrap(),
            Fq::one(),
        );
        let k = Fr::from_hex(&HEX_IV).unwrap();
        assert_eq!(r, p1 * k);
    }

    const G2_P1XY: [u8; 32] =
        hex!("83f6a65d85d51ec72eacf19bc38384e0369eb22a134a725a0191faa6e4f192ef");
    const G2_P1XX: [u8; 32] =
        hex!("9a79bfd491ef1cb32d9b57f7d0590ccff6b1cfe63dd15c0823d692fafbe96dbc");
    const G2_P1YY: [u8; 32] =
        hex!("9ed11c499291db0454d738555af0ce8a1df960056ee7425a6bf296eae60a5037");
    const G2_P1YX: [u8; 32] =
        hex!("849d4434eb7113fc9fb3809b51d54064fa2f20503423d256bc044905b1eba3fb");

    const G2_P2XY: [u8; 32] =
        hex!("a36232a9713f69157b7cdceef54aa0237b3ba0642a80dbb597af8935aea2c130");
    const G2_P2XX: [u8; 32] =
        hex!("624b19114e49f00281e2aee1f1b9d4f0a081a135868f8bbdb7b7a7b7da5fd6bc");
    const G2_P2YY: [u8; 32] =
        hex!("77966917ec1c5a294dd836c34691ab5e891f8c9f017443902c0a73ec54d449d8");
    const G2_P2YX: [u8; 32] =
        hex!("1be45454b6fa085a53744b22fd398238e400c3e031c8796e59e1bd6222048af0");

    const G2_DBLXY: [u8; 32] =
        hex!("73cbced58a8e76ef5235b480050a74e906e4d27185bd85d7ebdcd43ad24475fd");
    const G2_DBLXX: [u8; 32] =
        hex!("58400f0eb23000d814f5b5d0706749a72909795b7b04f26d6d58b2cf478ad9c9");
    const G2_DBLYY: [u8; 32] =
        hex!("19b460e09ac9ddbb380d6441e078a47bfcaa7d4c3d60b3a6c0d05f896472dc3c");
    const G2_DBLYX: [u8; 32] =
        hex!("1d69f785f47d6f25cb901b131612c37edc5e89ee9ba2dac8c401ced40e340a39");

    const G2_ADDXY: [u8; 32] =
        hex!("5f443752a19e368f404b89abae20a386d2b534c424b93ededdbfd04d4c569e6b");
    const G2_ADDXX: [u8; 32] =
        hex!("a411bbd84ee92a6ee53e5ca9cb81bacc192c6ba406f6fdcb2b04d0ab9c42ae44");
    const G2_ADDYY: [u8; 32] =
        hex!("6a3dadfcaac134e8353dd3abf37d487b206ca28dfab1e0a9376649df748f1605");
    const G2_ADDYX: [u8; 32] =
        hex!("4fa25e5e6100a023d4923df385dd236749c6a7f8e68db55e0bd1e2263fc04d28");

    const G2_MULXY: [u8; 32] =
        hex!("5d704de3261290dbba39dbd14e6bc416025240fd1ed65ec982efed685ae41e8b");
    const G2_MULXX: [u8; 32] =
        hex!("705c9ca4b5ef465c4e5db80ca4880627a6d9d6bcefd4756496baba9d5eaa3304");
    const G2_MULYY: [u8; 32] =
        hex!("4e96eb3543aabf1e9a65cae24177b9d13b0f7fae9472145ba7ae2b14bb447aef");
    const G2_MULYX: [u8; 32] =
        hex!("5d7ba50d7eac49a00b18fee2069afd3cc9719993fa78271e66b7a3efed46ac8b");

    const G2_MULGXY: [u8; 32] =
        hex!("920ef6fb3a2acff52aa0c004c18feca149dfd33d98086f8f402ea9e0de303c49");
    const G2_MULGXX: [u8; 32] =
        hex!("1f97dd359f2b065d63e0987f5bea2f3dc865c2cc112d7d161b46b83451716fd8");
    const G2_MULGYY: [u8; 32] =
        hex!("614881d4d05fef3173a4990465876c5200f58c5015e13354b23ae401c20c4aef");
    const G2_MULGYX: [u8; 32] =
        hex!("18a22e02b7d395a49f0646a79438e79cd37c32f163fe8923c13d56bab668e8a7");

    #[test]
    fn test_g2() {
        // println!("test_G2 test");
        group_trials::<G2>();
    }

    #[test]
    fn test_g2_dbl() {
        // println!("test_G2_dbl test");
        let x = Fq2::new(
            Fq::from_hex(&G2_P1XX).unwrap(),
            Fq::from_hex(&G2_P1XY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_hex(&G2_P1YX).unwrap(),
            Fq::from_hex(&G2_P1YY).unwrap(),
        );
        let p1 = G2::new(x, y, Fq2::one());
        let x = Fq2::new(
            Fq::from_hex(&G2_DBLXX).unwrap(),
            Fq::from_hex(&G2_DBLXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_hex(&G2_DBLYX).unwrap(),
            Fq::from_hex(&G2_DBLYY).unwrap(),
        );
        let r = G2::new(x, y, Fq2::one());
        assert_eq!(r, p1.double());
    }

    #[test]
    fn test_g2_add() {
        // println!("test_G2_add test");
        let x = Fq2::new(
            Fq::from_hex(&G2_P1XX).unwrap(),
            Fq::from_hex(&G2_P1XY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_hex(&G2_P1YX).unwrap(),
            Fq::from_hex(&G2_P1YY).unwrap(),
        );
        let p1 = G2::new(x, y, Fq2::one());
        let x = Fq2::new(
            Fq::from_hex(&G2_P2XX).unwrap(),
            Fq::from_hex(&G2_P2XY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_hex(&G2_P2YX).unwrap(),
            Fq::from_hex(&G2_P2YY).unwrap(),
        );
        let p2 = G2::new(x, y, Fq2::one());
        let x = Fq2::new(
            Fq::from_hex(&G2_ADDXX).unwrap(),
            Fq::from_hex(&G2_ADDXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_hex(&G2_ADDYX).unwrap(),
            Fq::from_hex(&G2_ADDYY).unwrap(),
        );
        let r = G2::new(x, y, Fq2::one());
        assert_eq!(r, p1 + p2);
    }

    #[test]
    fn test_g2_mul() {
        // println!("test_G2_mul test");
        let x = Fq2::new(
            Fq::from_hex(&G2_P1XX).unwrap(),
            Fq::from_hex(&G2_P1XY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_hex(&G2_P1YX).unwrap(),
            Fq::from_hex(&G2_P1YY).unwrap(),
        );
        let p1 = G2::new(x, y, Fq2::one());
        let x = Fq2::new(
            Fq::from_hex(&G2_MULXX).unwrap(),
            Fq::from_hex(&G2_MULXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_hex(&G2_MULYX).unwrap(),
            Fq::from_hex(&G2_MULYY).unwrap(),
        );
        let r = G2::new(x, y, Fq2::one());
        let k = Fr::from_hex(&HEX_IV).unwrap();
        assert_eq!(r, p1 * k);
    }

    #[test]
    fn test_g2_mulg() {
        // println!("test_G2_mulg test");
        let x = Fq2::new(
            Fq::from_hex(&G2_MULGXX).unwrap(),
            Fq::from_hex(&G2_MULGXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_hex(&G2_MULGYX).unwrap(),
            Fq::from_hex(&G2_MULGYY).unwrap(),
        );
        let r = G2::new(x, y, Fq2::one());
        let k = Fr::from_hex(&HEX_IV).unwrap();
        assert_eq!(r, G2::one() * k);
    }

    #[test]
    fn affine_fail() {
        let res = AffineG1::new(Fq::one(), Fq::one());
        assert!(
            res.is_err(),
            "Affine initialization should fail because the point is not on curve"
        );
    }

    #[test]
    fn affine_ok() {
        let res = AffineG1::new(Fq::one(), G1Params::coeff_b());
        assert!(
            res.is_err(),
            "Affine initialization should be ok because the point is on the curve"
        );
        let res = AffineG1::new(
            Fq::from_hex(&SM9_P1X).unwrap(),
            Fq::from_hex(&SM9_P1Y).unwrap(),
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_y_at_point_at_infinity() {
        assert!(G1::zero().y == Fq::one());
        assert!((-G1::zero()).y == Fq::one());

        assert!(G2::zero().y == Fq2::one());
        assert!((-G2::zero()).y == Fq2::one());
    }

    const G1_RAX: [u8; 32] =
        hex!("7CBA5B19069EE66AA79D490413D11846B9BA76DD22567F809CF23B6D964BB265");
    const G1_RAY: [u8; 32] =
        hex!("A9760C99CB6F706343FED05637085864958D6C90902ABA7D405FBEDF7B781599");

    const G2_DEBXY: [u8; 32] =
        hex!("74CCC3AC9C383C60AF083972B96D05C75F12C8907D128A17ADAFBAB8C5A4ACF7");
    const G2_DEBXX: [u8; 32] =
        hex!("01092FF4DE89362670C21711B6DBE52DCD5F8E40C6654B3DECE573C2AB3D29B2");
    const G2_DEBYY: [u8; 32] =
        hex!("44B0294AA04290E1524FF3E3DA8CFD432BB64DE3A8040B5B88D1B5FC86A4EBC1");
    const G2_DEBYX: [u8; 32] =
        hex!("8CFC48FB4FF37F1E27727464F3C34E2153861AD08E972D1625FC1A7BD18D5539");

    const PAIR2C11: [u8; 32] =
        hex!("28542FB6954C84BE6A5F2988A31CB6817BA0781966FA83D9673A9577D3C0C134");
    const PAIR2C10: [u8; 32] =
        hex!("5E27C19FC02ED9AE37F5BB7BE9C03C2B87DE027539CCF03E6B7D36DE4AB45CD1");
    const PAIR2C01: [u8; 32] =
        hex!("A1ABFCD30C57DB0F1A838E3A8F2BF823479C978BD137230506EA6249C891049E");
    const PAIR2C00: [u8; 32] =
        hex!("3497477913AB89F5E2960F382B1B5C8EE09DE0FA498BA95C4409D630D343DA40");

    const PAIR2B11: [u8; 32] =
        hex!("4FEC93472DA33A4DB6599095C0CF895E3A7B993EE5E4EBE3B9AB7D7D5FF2A3D1");
    const PAIR2B10: [u8; 32] =
        hex!("647BA154C3E8E185DFC33657C1F128D480F3F7E3F16801208029E19434C733BB");
    const PAIR2B01: [u8; 32] =
        hex!("73F21693C66FC23724DB26380C526223C705DAF6BA18B763A68623C86A632B05");
    const PAIR2B00: [u8; 32] =
        hex!("0F63A071A6D62EA45B59A1942DFF5335D1A232C9C5664FAD5D6AF54C11418B0D");

    const PAIR2A11: [u8; 32] =
        hex!("8C8E9D8D905780D50E779067F2C4B1C8F83A8B59D735BB52AF35F56730BDE5AC");
    const PAIR2A10: [u8; 32] =
        hex!("861CCD9978617267CE4AD9789F77739E62F2E57B48C2FF26D2E90A79A1D86B93");
    const PAIR2A01: [u8; 32] =
        hex!("9B1CA08F64712E33AEDA3F44BD6CB633E0F722211E344D73EC9BBEBC92142765");
    const PAIR2A00: [u8; 32] =
        hex!("6BA584CE742A2A3AB41C15D3EF94EDEB8EF74A2BDCDAAECC09ABA567981F6437");

    #[test]
    fn test_pairing() {
        let p1 = G1::new(
            Fq::from_hex(&G1_RAX).unwrap(),
            Fq::from_hex(&G1_RAY).unwrap(),
            Fq::one(),
        );

        let x = Fq2::new(
            Fq::from_hex(&G2_DEBXX).unwrap(),
            Fq::from_hex(&G2_DEBXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_hex(&G2_DEBYX).unwrap(),
            Fq::from_hex(&G2_DEBYY).unwrap(),
        );
        let p2 = G2::new(x, y, Fq2::one());

        let r0 = Fq2::new(
            Fq::from_hex(&PAIR2A00).unwrap(),
            Fq::from_hex(&PAIR2A01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_hex(&PAIR2A10).unwrap(),
            Fq::from_hex(&PAIR2A11).unwrap(),
        );
        let a = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_hex(&PAIR2B00).unwrap(),
            Fq::from_hex(&PAIR2B01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_hex(&PAIR2B10).unwrap(),
            Fq::from_hex(&PAIR2B11).unwrap(),
        );
        let b = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_hex(&PAIR2C00).unwrap(),
            Fq::from_hex(&PAIR2C01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_hex(&PAIR2C10).unwrap(),
            Fq::from_hex(&PAIR2C11).unwrap(),
        );
        let c = Fq4::new(r0, r1);
        let r = Fq12::new(a, b, c);

        let alice = pairing(&p1, &p2);

        assert_eq!(alice, r);
    }
}
