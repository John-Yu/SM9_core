#![allow(dead_code)]

use core::{
    fmt,
    ops::{Add, AddAssign, Mul, Neg, Sub},
};
use hex_literal::hex;
use rand::Rng;

use crate::{
    One, Zero,
    fields::{FieldElement, Fq, Fq2, Fr},
    u256::U256,
};

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

pub trait GroupElement:
    Sized
    + Copy
    + Clone
    + Zero
    + PartialEq
    + Eq
    + fmt::Debug
    + Add<Output = Self>
    + Sub<Output = Self>
    + Neg<Output = Self>
    + Mul<Fr, Output = Self>
{
    fn one() -> Self;
    fn random<R: Rng>(rng: &mut R) -> Self;
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
    pub(crate) x: P::Base,
    pub(crate) y: P::Base,
    pub(crate) z: P::Base,
}

impl<P: GroupParams> G<P> {
    pub fn new(x: P::Base, y: P::Base, z: P::Base) -> Self {
        G { x, y, z }
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
                    x,
                    y,
                    z: P::Base::one(),
                };

                if (p * (-Fr::one())) + p != G::zero() {
                    return Err(Error::NotInSubgroup);
                }
            }

            Ok(AffineG { x, y })
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
        *self
    }
}

impl<P: GroupParams> Copy for G<P> {}

impl<P: GroupParams> Clone for AffineG<P> {
    fn clone(&self) -> Self {
        *self
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

        true
    }
}
impl<P: GroupParams> Eq for G<P> {}

impl<P: GroupParams> G<P> {
    pub fn to_affine(self) -> Option<AffineG<P>> {
        if self.z.is_zero() {
            None
        } else if self.z == P::Base::one() {
            Some(AffineG {
                x: self.x,
                y: self.y,
            })
        } else {
            let zinv = self.z.inverse()?;
            let zinv_squared = zinv.squared();

            Some(AffineG {
                x: self.x * zinv_squared,
                y: self.y * (zinv_squared * zinv),
            })
        }
    }
}

impl<P: GroupParams> AffineG<P> {
    pub fn to_jacobian(self) -> G<P> {
        G {
            x: self.x,
            y: self.y,
            z: P::Base::one(),
        }
    }
}

impl<P: GroupParams> Zero for G<P> {
    fn zero() -> Self {
        G {
            x: P::Base::zero(),
            y: P::Base::one(),
            z: P::Base::zero(),
        }
    }

    fn is_zero(&self) -> bool {
        self.z.is_zero()
    }
}

impl<P: GroupParams> GroupElement for G<P> {
    fn one() -> Self {
        P::one()
    }

    fn random<R: Rng>(rng: &mut R) -> Self {
        P::one() * Fr::random(rng)
    }

    // SM9 dentity-based cryptographic algorithms
    // Part 1: General
    // Annex A  A.1.3.3.2
    // 𝜆1=3𝑥1^2, 𝜆2=4𝑥1𝑦1^2, 𝜆3=8𝑦1^4, 𝑥3=𝜆1^2−2𝜆2, 𝑦3=𝜆1(𝜆2−𝑥3)−𝜆3, 𝑧3=2𝑦1𝑧1
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
    fn double(&self) -> Self {
        let a = self.x.squared();
        let b = self.y.squared();
        let c = b.squared();
        let d = ((self.x + b).squared() - a - c).double();
        let e = a.triple();
        let f = e.squared();
        let x3 = f - d.double();
        let eight_c = c.double().double().double();
        let y1z1 = self.y * self.z;

        G {
            x: x3,
            y: e * (d - x3) - eight_c,
            z: y1z1.double(),
        }
    }
}
impl<P: GroupParams> Mul<Fr> for G<P> {
    type Output = G<P>;
    #[allow(clippy::suspicious_arithmetic_impl)]
    /// Standard double-and-add method for multiplication by a scalar.
    fn mul(self, other: Fr) -> G<P> {
        let mut res = G::zero();
        for i in U256::from(other).bits_without_leading_zeros() {
            res = res.double();
            if i {
                res += self;
            }
        }

        res
    }
}
impl<P: GroupParams> AddAssign<G<P>> for G<P> {
    fn add_assign(&mut self, other: G<P>) {
        *self = *self + other;
    }
}
impl<P: GroupParams> AddAssign<&G<P>> for G<P> {
    fn add_assign(&mut self, rhs: &G<P>) {
        *self = *self + rhs;
    }
}
impl<P: GroupParams> Add<G<P>> for G<P> {
    type Output = G<P>;
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General
    // Annex A  A.1.3.3.2
    // 𝜆1=𝑥1𝑧2^2, 𝜆2=𝑥2𝑧1^2, 𝜆3=𝜆1−𝜆2, 𝜆4=𝑦1𝑧2^3, 𝜆5=𝑦2𝑧1^3, 𝜆6=𝜆4−𝜆5, 𝜆7=𝜆1+𝜆2, 𝜆8=𝜆4+𝜆5, 𝜆9=𝜆7𝜆3^2,
    // 𝑥3=𝜆6^2−𝜆9, 𝜆10=𝜆9^2−2𝑥3, 𝑦3=(𝜆10𝜆6−𝜆8𝜆3^3)/2, 𝑧3=𝑧1𝑧2𝜆3
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-1998-cmo-2
    fn add(self, other: G<P>) -> G<P> {
        if self.is_zero() {
            return other;
        }
        if other.is_zero() {
            return self;
        }
        match (self.z == P::Base::one(), other.z == P::Base::one()) {
            (true, true) => {
                let h = other.x - self.x;
                let r = other.y - self.y;
                if r.is_zero() && h.is_zero() {
                    return self.double();
                }
                let hh = h.squared();
                let hhh = h * hh;
                let v = self.x * hh;
                let x = r.squared() - hhh - v.double();
                let y = r * (v - x) - self.y * hhh;
                let z = h;
                G { x, y, z }
            }
            (false, true) => {
                let z1_squared = self.z.squared();
                let u2 = other.x * z1_squared;
                let z1_cubed = self.z * z1_squared;
                let s2 = other.y * z1_cubed;
                let h = u2 - self.x;
                let r = s2 - self.y;
                if r.is_zero() && h.is_zero() {
                    return self.double();
                }
                let hh = h.squared();
                let hhh = h * hh;
                let v = self.x * hh;
                let x = r.squared() - hhh - v.double();
                let y = r * (v - x) - self.y * hhh;
                let z = self.z * h;
                G { x, y, z }
            }
            (true, false) => other + self,
            (false, false) => {
                let z1_squared = self.z.squared();
                let z2_squared = other.z.squared();
                let u1 = self.x * z2_squared;
                let u2 = other.x * z1_squared;
                let z1_cubed = self.z * z1_squared;
                let z2_cubed = other.z * z2_squared;
                let s1 = self.y * z2_cubed;
                let s2 = other.y * z1_cubed;

                let r = s2 - s1;
                let h = u2 - u1;
                let t6 = s1 + s2;
                if r.is_zero() {
                    if h.is_zero() {
                        return self.double();
                    }
                    if t6.is_zero() {
                        return Self::zero();
                    }
                }
                let hh = h.squared();
                let hhh = h * hh;
                let v = u1 * hh;
                let x = r.squared() - hhh - v.double();
                let y = r * (v - x) - s1 * hhh;
                let z = self.z * other.z * h;
                G { x, y, z }
            }
        }
    }
}
impl<P: GroupParams> Add<&G<P>> for G<P> {
    type Output = G<P>;

    fn add(self, rhs: &G<P>) -> Self::Output {
        self + *rhs
    }
}
impl<P: GroupParams> Add<G<P>> for &G<P> {
    type Output = G<P>;

    fn add(self, rhs: G<P>) -> Self::Output {
        *self + rhs
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
            x: Fq::from_slice(&SM9_P1X).unwrap(),
            y: Fq::from_slice(&SM9_P1Y).unwrap(),
            z: Fq::one(),
        }
    }

    fn coeff_b() -> Fq {
        Fq::from_str("5").expect("5 is a valid field element")
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
                Fq::from_slice(&SM9_P2X0).unwrap(),
                Fq::from_slice(&SM9_P2X1).unwrap(),
            ),
            y: Fq2::new(
                Fq::from_slice(&SM9_P2Y0).unwrap(),
                Fq::from_slice(&SM9_P2Y1).unwrap(),
            ),
            z: Fq2::one(),
        }
    }

    fn coeff_b() -> Fq2 {
        Fq2::i().scale(&Fq::from_str("5").expect("5 is a valid field element"))
    }

    fn check_order() -> bool {
        true
    }
}

pub type G2 = G<G2Params>;

pub type AffineG2 = AffineG<G2Params>;

/* ************************************************************************************************ */
#[cfg(test)]
mod tests {
    use super::*;

    fn random_test_addition<G: GroupElement, R: Rng>(rng: &mut R) {
        for _ in 0..50 {
            let r1 = G::random(rng);
            let r2 = G::random(rng);
            let r3 = G::random(rng);

            assert_eq!((r1 + r2) + r3, r1 + (r2 + r3));
            assert!(((r1 + r2 + r3) - r2 - r3 - r1).is_zero());
        }
    }
    fn random_test_doubling<G: GroupElement, R: Rng>(rng: &mut R) {
        for _ in 0..50 {
            let r1 = G::random(rng);
            let r2 = G::random(rng);
            let ti = Fr::from_str("2").unwrap().inverse().unwrap();

            assert_eq!((r1 + r2) + r1, r1.double() + r2);
            assert_eq!(r1, r1.double() * ti);
        }
    }
    pub fn group_trials<G: GroupElement>() {
        assert!(G::zero().is_zero());
        assert!((G::one() - G::one()).is_zero());
        assert_eq!(G::one() + G::one(), G::one() * Fr::from_str("2").unwrap());
        assert!(G::zero().double().is_zero());

        assert!((G::one() * (-Fr::one()) + G::one()).is_zero());
        use rand::{SeedableRng, rngs::StdRng};
        let seed = [
            0, 0, 0, 0, 0, 0, 64, 13, // 103245
            0, 0, 0, 0, 0, 0, 176, 2, // 191922
            0, 0, 0, 0, 0, 0, 0, 13, // 1293
            0, 0, 0, 0, 0, 0, 96, 7u8, // 192103
        ];
        let mut rng = StdRng::from_seed(seed);

        random_test_addition::<G, _>(&mut rng);
        random_test_doubling::<G, _>(&mut rng);
    }

    fn test_projective_negation_and_subtraction<G: GroupElement>() {
        let a = G::one().double();
        assert_eq!(a + (-a), G::zero());
        assert_eq!(a + (-a), a - a);
    }

    fn test_projective_double_and_add<G: GroupElement>() {
        let a = G::one().double().double(); // 4P
        let b = G::one().double(); // 2P
        let c = a + b; //6P

        let mut d = G::one();
        for _ in 0..5 {
            d = d + G::one();
        }
        assert_eq!(c, d);
    }

    #[test]
    fn test_g1_negation_and_subtraction() {
        test_projective_negation_and_subtraction::<G1>();
    }
    #[test]
    fn test_g2_negation_and_subtraction() {
        test_projective_negation_and_subtraction::<G2>();
    }

    #[test]
    fn test_g1() {
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
        let p1 = G1::new(
            Fq::from_slice(&P1X).unwrap(),
            Fq::from_slice(&P1Y).unwrap(),
            Fq::one(),
        );
        let r = G1::new(
            Fq::from_slice(&P_DBL_X).unwrap(),
            Fq::from_slice(&P_DBL_Y).unwrap(),
            Fq::one(),
        );
        assert_eq!(r, p1.double());
    }

    #[test]
    fn test_g1_add() {
        let p1 = G1::new(
            Fq::from_slice(&P1X).unwrap(),
            Fq::from_slice(&P1Y).unwrap(),
            Fq::one(),
        );
        let p2 = G1::new(
            Fq::from_slice(&P2X).unwrap(),
            Fq::from_slice(&P2Y).unwrap(),
            Fq::one(),
        );
        let r = G1::new(
            Fq::from_slice(&P_ADD_X).unwrap(),
            Fq::from_slice(&P_ADD_Y).unwrap(),
            Fq::one(),
        );
        assert_eq!(r, p1 + p2);
        test_projective_double_and_add::<G1>();
    }

    #[test]
    fn test_g1_sub() {
        let p1 = G1::new(
            Fq::from_slice(&P1X).unwrap(),
            Fq::from_slice(&P1Y).unwrap(),
            Fq::one(),
        );
        let p2 = G1::new(
            Fq::from_slice(&P2X).unwrap(),
            Fq::from_slice(&P2Y).unwrap(),
            Fq::one(),
        );
        let r = G1::new(
            Fq::from_slice(&P_NEG_X).unwrap(),
            Fq::from_slice(&P_NEG_Y).unwrap(),
            Fq::one(),
        );
        assert_eq!(r, -p1);
        let r = G1::new(
            Fq::from_slice(&P_SUB_X).unwrap(),
            Fq::from_slice(&P_SUB_Y).unwrap(),
            Fq::one(),
        );
        assert_eq!(r, p1 - p2);
    }

    #[test]
    fn test_g1_mul() {
        let p1 = G1::new(
            Fq::from_slice(&P1X).unwrap(),
            Fq::from_slice(&P1Y).unwrap(),
            Fq::one(),
        );
        let r = G1::new(
            Fq::from_slice(&P_MUL_X).unwrap(),
            Fq::from_slice(&P_MUL_Y).unwrap(),
            Fq::one(),
        );
        let k = Fr::from_slice(&HEX_IV).unwrap();
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
        group_trials::<G2>();
    }

    #[test]
    fn test_g2_dbl() {
        let x = Fq2::new(
            Fq::from_slice(&G2_P1XX).unwrap(),
            Fq::from_slice(&G2_P1XY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&G2_P1YX).unwrap(),
            Fq::from_slice(&G2_P1YY).unwrap(),
        );
        let p1 = G2::new(x, y, Fq2::one());
        let x = Fq2::new(
            Fq::from_slice(&G2_DBLXX).unwrap(),
            Fq::from_slice(&G2_DBLXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&G2_DBLYX).unwrap(),
            Fq::from_slice(&G2_DBLYY).unwrap(),
        );
        let r = G2::new(x, y, Fq2::one());
        assert_eq!(r, p1.double());
    }

    #[test]
    fn test_g2_add() {
        let x = Fq2::new(
            Fq::from_slice(&G2_P1XX).unwrap(),
            Fq::from_slice(&G2_P1XY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&G2_P1YX).unwrap(),
            Fq::from_slice(&G2_P1YY).unwrap(),
        );
        let p1 = G2::new(x, y, Fq2::one());
        let x = Fq2::new(
            Fq::from_slice(&G2_P2XX).unwrap(),
            Fq::from_slice(&G2_P2XY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&G2_P2YX).unwrap(),
            Fq::from_slice(&G2_P2YY).unwrap(),
        );
        let p2 = G2::new(x, y, Fq2::one());
        let x = Fq2::new(
            Fq::from_slice(&G2_ADDXX).unwrap(),
            Fq::from_slice(&G2_ADDXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&G2_ADDYX).unwrap(),
            Fq::from_slice(&G2_ADDYY).unwrap(),
        );
        let r = G2::new(x, y, Fq2::one());
        assert_eq!(r, p1 + p2);
        test_projective_double_and_add::<G2>();
    }

    #[test]
    fn test_g2_mul() {
        let x = Fq2::new(
            Fq::from_slice(&G2_P1XX).unwrap(),
            Fq::from_slice(&G2_P1XY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&G2_P1YX).unwrap(),
            Fq::from_slice(&G2_P1YY).unwrap(),
        );
        let p1 = G2::new(x, y, Fq2::one());
        let x = Fq2::new(
            Fq::from_slice(&G2_MULXX).unwrap(),
            Fq::from_slice(&G2_MULXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&G2_MULYX).unwrap(),
            Fq::from_slice(&G2_MULYY).unwrap(),
        );
        let r = G2::new(x, y, Fq2::one());
        let k = Fr::from_slice(&HEX_IV).unwrap();
        assert_eq!(r, p1 * k);
    }

    #[test]
    fn test_g2_mulg() {
        let x = Fq2::new(
            Fq::from_slice(&G2_MULGXX).unwrap(),
            Fq::from_slice(&G2_MULGXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&G2_MULGYX).unwrap(),
            Fq::from_slice(&G2_MULGYY).unwrap(),
        );
        let r = G2::new(x, y, Fq2::one());
        let k = Fr::from_slice(&HEX_IV).unwrap();
        assert_eq!(r, G2::one() * k);
    }

    #[test]
    fn affine_fail() {
        let res = AffineG1::new(Fq::one(), Fq::one());
        assert!(
            res.is_err(),
            "Affine initialization should fail because the point is not on curve"
        );
        let res = AffineG1::new(Fq::one(), G1Params::coeff_b());
        assert!(
            res.is_err(),
            "Affine initialization should fail because the point is not on curve"
        );
    }

    #[test]
    fn affine_ok() {
        let res = AffineG1::new(
            Fq::from_slice(&SM9_P1X).unwrap(),
            Fq::from_slice(&SM9_P1Y).unwrap(),
        );
        assert!(
            res.is_ok(),
            "Affine initialization should be ok because the point is on the curve"
        );
    }
    #[test]
    fn test_projective_point_equality() {
        let a = G1::one();
        let b = G1::zero();

        assert!(a == a);
        assert!(b == b);
        assert!(a != b);
        assert!(b != a);
    }
    #[test]
    fn test_y_at_point_at_infinity() {
        assert!(G1::zero().y == Fq::one());
        assert!((-G1::zero()).y == Fq::one());

        assert!(G2::zero().y == Fq2::one());
        assert!((-G2::zero()).y == Fq2::one());
    }
}
