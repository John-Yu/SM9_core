//! # sm9_core
//!
//! A pairing cryptography library written in pure Rust.
//! It makes use of the Barreto-Naehrig (BN) curve construction from "SM9 identity-based cryptographic algorithms" to provide two cyclic groups G<sub>1</sub> and G<sub>2</sub> , with a R-ate pairing.
//! - no_std
//! - no unsafe{}
//!
#![no_std]
#![forbid(unsafe_code)]
#![doc = include_str!("../README.md")]

#[cfg(test)]
#[macro_use]
extern crate std;

extern crate alloc;
extern crate hex_literal;
extern crate rand;

mod arith;
mod fields;
mod groups;
mod pairings;
mod u256;
mod u512;

use alloc::fmt::Debug;
use core::ops::{Add, Mul, Neg, Sub};
use rand::Rng;

use crate::fields::FieldElement;
use crate::groups::{G1Params, G2Params, GroupElement, GroupParams};
pub use crate::pairings::G2Prepared;
use crate::u256::U256;

/// Represents an element of the finite field F<sub>r</sub>
// where r = 0xB640000002A3A6F1D603AB4FF58EC74449F2934B18EA8BEEE56EE19CD69ECF25
// Elements of Fr are always in Montgomery form; i.e., Fr(a) = aR mod r, with R = 2^256.

/// # Examples
///
/// ```rust
/// use sm9_core::*;
/// use hex_literal::hex;
///
/// let ks = Fr::from_slice(&hex!("000130E78459D78545CB54C587E02CF480CE0B66340F319F348A1D5B1F2DC5F4")).unwrap();
/// let r = Fr::from_slice(&hex!("00033C86 16B06704 813203DF D0096502 2ED15975 C662337A ED648835 DC4B1CBE")).unwrap();
/// let pub_s = G2::one() * ks;
/// let g = pairing(G1::one(), pub_s).pow(r);
/// let r1 = g.to_slice();
/// println!("{:#?}", g);
/// let r0 = hex!(
/// "81377B8F DBC2839B 4FA2D0E0 F8AA6853 BBBE9E9C 4099608F 8612C607 8ACD7563"
/// "815AEBA2 17AD502D A0F48704 CC73CABB 3C06209B D87142E1 4CBD99E8 BCA1680F"
/// "30DADC5C D9E207AE E32209F6 C3CA3EC0 D800A1A4 2D33C731 53DED47C 70A39D2E"
/// "8EAF5D17 9A1836B3 59A9D1D9 BFC19F2E FCDB8293 28620962 BD3FDF15 F2567F58"
/// "A543D256 09AE9439 20679194 ED30328B B33FD156 60BDE485 C6B79A7B 32B01398"
/// "3F012DB0 4BA59FE8 8DB88932 1CC2373D 4C0C35E8 4F7AB1FF 33679BCA 575D6765"
/// "4F8624EB 435B838C CA77B2D0 347E65D5 E4696441 2A096F41 50D8C5ED E5440DDF"
/// "0656FCB6 63D24731 E8029218 8A2471B8 B68AA993 89926849 9D23C897 55A1A897"
/// "44643CEA D40F0965 F28E1CD2 895C3D11 8E4F65C9 A0E3E741 B6DD52C0 EE2D25F5"
/// "898D6084 8026B7EF B8FCC1B2 442ECF07 95F8A81C EE99A624 8F294C82 C90D26BD"
/// "6A814AAF 475F128A EF43A128 E37F8015 4AE6CB92 CAD7D150 1BAE30F7 50B3A9BD"
/// "1F96B08E 97997363 91131470 5BFB9A9D BB97F755 53EC90FB B2DDAE53 C8F68E42"
/// );
/// assert_eq!(r0, r1);
/// // test  fast_pairing
/// let g1 = fast_pairing(G1::one(), pub_s).pow(r);
/// let r1 = g1.to_slice();
/// assert_eq!(r0, r1);
/// ```

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fr(pub(crate) fields::Fr);

impl Fr {
    /// Returns zero, the additive identity.
    pub fn zero() -> Self {
        Fr(fields::Fr::zero())
    }
    /// Returns one, the multiplicative identity.
    pub fn one() -> Self {
        Fr(fields::Fr::one())
    }
    /// Exponentiates `self` by `exp`, where `exp` is a
    /// Fr element exponent.
    pub fn pow(&self, exp: Fr) -> Self {
        Fr(self.0.pow(exp.0))
    }
    /// Attempts to convert a string base 10 to a element of `Fr`
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Option<Self> {
        fields::Fr::from_str(s).map(Fr)
    }
    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub fn inverse(&self) -> Option<Self> {
        self.0.inverse().map(Fr)
    }
    /// Get a random element
    pub fn random<R: Rng>(rng: &mut R) -> Self {
        Fr(fields::Fr::random(rng))
    }
    /// Returns true if element is the additive identity.
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    pub fn interpret(buf: &[u8; 64]) -> Fr {
        Fr(fields::Fr::interpret(buf))
    }
    /// Attempts to convert a big-endian byte representation of
    /// a field element into an element of `Fr`, failing if the input
    /// is not canonical (is not smaller than r or Invalid Slice Length).
    pub fn from_slice(slice: &[u8]) -> Result<Self, FieldError> {
        U256::from_slice(slice)
            .map_err(|_| FieldError::InvalidSliceLength) // todo: maybe more sensful error handling
            .map(Fr::new_mul_factor)
    }
    /// for H1() and H2(), Attempts to convert a HASH result (40 bytes) to a Fr element
    pub fn from_hash(hex: &[u8]) -> Option<Self> {
        fields::Fr::from_hash(hex).map(Fr)
    }
    /// Converts an element of `Fr` into a byte representation in
    /// big-endian byte order.
    pub fn to_slice(self) -> [u8; 32] {
        self.0.to_slice()
    }
    pub fn new(val: U256) -> Option<Self> {
        fields::Fr::new(val).map(Fr)
    }
    pub fn new_mul_factor(val: U256) -> Self {
        Fr(fields::Fr::new_mul_factor(val))
    }
    /// Converts an element into a 256-bit big-endian integer(U256)
    /// Turn into canonical form by computing: (a.R) / R = a
    pub fn into_u256(self) -> U256 {
        (self.0).into()
    }
    pub fn set_bit(&mut self, bit: usize, to: bool) {
        self.0.set_bit(bit, to);
    }
}

impl Add<Fr> for Fr {
    type Output = Fr;
    /// Adds this element to another element.
    fn add(self, other: Fr) -> Fr {
        Fr(self.0 + other.0)
    }
}

impl Sub<Fr> for Fr {
    type Output = Fr;
    /// Subtracts another element from this element.
    fn sub(self, other: Fr) -> Fr {
        Fr(self.0 - other.0)
    }
}

impl Neg for Fr {
    type Output = Fr;
    /// Negates this element.
    fn neg(self) -> Fr {
        Fr(-self.0)
    }
}

impl Mul for Fr {
    type Output = Fr;
    /// Multiplies this element by another element
    fn mul(self, other: Fr) -> Fr {
        Fr(self.0 * other.0)
    }
}

#[derive(Debug)]
pub enum FieldError {
    InvalidSliceLength,
    InvalidU512Encoding,
    NotMember,
}

#[derive(Debug)]
pub enum CurveError {
    InvalidEncoding,
    NotMember,
    Field(FieldError),
    ToAffineConversion,
}

impl From<FieldError> for CurveError {
    fn from(fe: FieldError) -> Self {
        CurveError::Field(fe)
    }
}

pub use crate::groups::Error as GroupError;

/// Represents an element of the finite field F<sub>q</sub>
// where q = 0xB640000002A3A6F1D603AB4FF58EC74521F2934B1A7AEEDBE56F9B27E351457D
// Elements of Fq are always in Montgomery form; i.e., Fq(a) = aR mod q, with R = 2^256.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq(fields::Fq);

impl Fq {
    /// Returns zero, the additive identity.
    pub fn zero() -> Self {
        Fq(fields::Fq::zero())
    }
    /// Returns one, the multiplicative identity.
    pub fn one() -> Self {
        Fq(fields::Fq::one())
    }
    /// Exponentiates `self` by `exp`, where `exp` is a
    /// Fq element exponent.
    pub fn pow(&self, exp: Fq) -> Self {
        Fq(self.0.pow(exp.0))
    }
    /// Attempts to convert a string base 10 to a element of `Fq`
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(s: &str) -> Option<Self> {
        fields::Fq::from_str(s).map(Fq)
    }
    /// Attempts to convert a big-endian byte representation of
    /// a field element into an element of `Fq`,
    pub fn from_slice(hex: &[u8]) -> Option<Self> {
        fields::Fq::from_slice(hex).map(Fq)
    }
    /// Converts an element of `Fq` into a byte representation in
    /// big-endian byte order.
    pub fn to_slice(self) -> [u8; 32] {
        self.0.to_slice()
    }
    /// Computes the multiplicative inverse of this element,
    /// failing if the element is zero.
    pub fn inverse(&self) -> Option<Self> {
        self.0.inverse().map(Fq)
    }
    /// Returns true if element is the additive identity.
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    pub fn is_even(&self) -> bool {
        self.into_u256().is_even()
    }
    pub fn interpret(buf: &[u8; 64]) -> Fq {
        Fq(fields::Fq::interpret(buf))
    }
    /// Converts an element into a slice of bytes in
    /// big-endian byte order.
    pub fn to_big_endian(self, slice: &mut [u8]) -> Result<(), FieldError> {
        self.into_u256()
            .to_big_endian(slice)
            .map_err(|_| FieldError::InvalidSliceLength)
    }
    pub fn from_u256(u256: U256) -> Result<Self, FieldError> {
        Ok(Fq(fields::Fq::new(u256).ok_or(FieldError::NotMember)?))
    }
    /// Converts an element into a 256-bit big-endian integer(U256)
    /// Turn into canonical form by computing: (a.R) / R = a
    pub fn into_u256(self) -> U256 {
        (self.0).into()
    }
    pub fn modulus() -> U256 {
        fields::Fq::modulus()
    }
    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> Option<Self> {
        self.0.sqrt().map(Fq)
    }
}

impl Add<Fq> for Fq {
    type Output = Fq;

    fn add(self, other: Fq) -> Fq {
        Fq(self.0 + other.0)
    }
}

impl Sub<Fq> for Fq {
    type Output = Fq;

    fn sub(self, other: Fq) -> Fq {
        Fq(self.0 - other.0)
    }
}

impl Neg for Fq {
    type Output = Fq;

    fn neg(self) -> Fq {
        Fq(-self.0)
    }
}

impl Mul for Fq {
    type Output = Fq;

    fn mul(self, other: Fq) -> Fq {
        Fq(self.0 * other.0)
    }
}
/// Represents an element of the finite field F<sub>q</sub>2
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Fq2(fields::Fq2);

impl Fq2 {
    pub fn one() -> Fq2 {
        Fq2(fields::Fq2::one())
    }
    pub fn zero() -> Fq2 {
        Fq2(fields::Fq2::zero())
    }
    /// Initalizes new Fq2(a + bi, a is real coeff, b is imaginary)
    pub fn new(a: Fq, b: Fq) -> Fq2 {
        Fq2(fields::Fq2::new(a.0, b.0))
    }
    /// Returns true if element is the additive identity.
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    pub fn is_even(&self) -> bool {
        self.real().into_u256().is_even()
    }
    /// Exponentiates `self` by `exp`, where `exp` is a
    /// U256 element exponent.
    pub fn pow(&self, exp: U256) -> Self {
        Fq2(self.0.pow(exp))
    }
    pub fn real(&self) -> Fq {
        Fq(*self.0.real())
    }
    pub fn imaginary(&self) -> Fq {
        Fq(*self.0.imaginary())
    }
    /// Computes the square root of this element, if it exists.
    pub fn sqrt(&self) -> Option<Self> {
        self.0.sqrt().map(Fq2)
    }
    /// Attempts to convert a big-endian byte representation of
    /// a field element into an element of `Fq2`,
    pub fn from_slice(hex: &[u8]) -> Option<Self> {
        fields::Fq2::from_slice(hex)
            .map(Fq2)
            .map_err(|_| FieldError::InvalidSliceLength)
            .ok()
    }
    /// Converts an element into a slice of bytes in
    /// big-endian byte order.
    pub fn to_slice(self) -> [u8; 64] {
        self.0.to_slice()
    }
}

impl Add<Fq2> for Fq2 {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Fq2(self.0 + other.0)
    }
}

impl Sub<Fq2> for Fq2 {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Fq2(self.0 - other.0)
    }
}

impl Neg for Fq2 {
    type Output = Self;

    fn neg(self) -> Self {
        Fq2(-self.0)
    }
}

impl Mul for Fq2 {
    type Output = Self;

    fn mul(self, other: Self) -> Self {
        Fq2(self.0 * other.0)
    }
}

pub trait Group:
    Send
    + Sync
    + Copy
    + Clone
    + PartialEq
    + Eq
    + Sized
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Neg<Output = Self>
    + Mul<Fr, Output = Self>
{
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    fn normalize(&mut self);
}

/// an additive cyclic group of prime order r

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct G1(pub groups::G1);

impl G1 {
    pub fn new(x: Fq, y: Fq, z: Fq) -> Self {
        G1(groups::G1::new(x.0, y.0, z.0))
    }

    pub fn x(&self) -> Fq {
        Fq(*self.0.x())
    }

    pub fn set_x(&mut self, x: Fq) {
        *self.0.x_mut() = x.0
    }

    pub fn y(&self) -> Fq {
        Fq(*self.0.y())
    }

    pub fn set_y(&mut self, y: Fq) {
        *self.0.y_mut() = y.0
    }

    pub fn z(&self) -> Fq {
        Fq(*self.0.z())
    }

    pub fn set_z(&mut self, z: Fq) {
        *self.0.z_mut() = z.0
    }

    pub fn b() -> Fq {
        Fq(G1Params::coeff_b())
    }
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General
    // Annex A  A.4.2
    pub fn from_compressed(bytes: &[u8]) -> Result<Self, CurveError> {
        if bytes.len() != 33 {
            return Err(CurveError::InvalidEncoding);
        }

        let sign = bytes[0];
        debug_assert!(sign == 2 || sign == 3);
        let x = Fq::from_slice(&bytes[1..]).ok_or(CurveError::InvalidEncoding)?;
        let y_squared = (x * x * x) + Self::b();
        let mut y = y_squared.sqrt().ok_or(CurveError::NotMember)?;
        let is_even = sign & 1 == 0;

        if is_even != y.is_even() {
            y = -y;
        }

        AffineG1::new(x, y)
            .map_err(|_| CurveError::NotMember)
            .map(Into::into)
    }
    /// Serializes this element into compressed form.
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General 6.2.8
    pub fn to_compressed(self) -> [u8; 33] {
        let mut res = [0u8; 33];
        let ag1 = AffineG1::from_jacobian(self).unwrap();
        res[0] = 2; // compressed
        if !ag1.y().is_even() {
            res[0] |= 1;
        }
        res[1..].copy_from_slice(&ag1.x().to_slice());

        res
    }
    /// Serializes this element into uncompressed form.
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General, 6.2.8 , case 1
    pub fn to_slice(self) -> [u8; 64] {
        let mut res = [0u8; 64];
        let ag1 = AffineG1::from_jacobian(self).unwrap();
        res[..32].copy_from_slice(&ag1.x().to_slice());
        res[32..].copy_from_slice(&ag1.y().to_slice());

        res
    }
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General, 6.2.8 , case 1
    pub fn from_slice(bytes: &[u8]) -> Result<Self, CurveError> {
        if bytes.len() != 64 {
            return Err(CurveError::InvalidEncoding);
        }

        let x = Fq::from_slice(&bytes[..32]).ok_or(CurveError::InvalidEncoding)?;
        let y = Fq::from_slice(&bytes[32..]).ok_or(CurveError::InvalidEncoding)?;

        AffineG1::new(x, y)
            .map_err(|_| CurveError::NotMember)
            .map(Into::into)
    }
}

impl Group for G1 {
    fn zero() -> Self {
        G1(groups::G1::zero())
    }
    fn one() -> Self {
        G1(groups::G1::one())
    }
    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn normalize(&mut self) {
        let new = match self.0.to_affine() {
            Some(a) => a,
            None => return,
        };

        self.0 = new.to_jacobian();
    }
}

impl Add<G1> for G1 {
    type Output = G1;

    fn add(self, other: G1) -> G1 {
        G1(self.0 + other.0)
    }
}

impl Sub<G1> for G1 {
    type Output = G1;

    fn sub(self, other: G1) -> G1 {
        G1(self.0 - other.0)
    }
}

impl Neg for G1 {
    type Output = G1;

    fn neg(self) -> G1 {
        G1(-self.0)
    }
}

impl Mul<Fr> for G1 {
    type Output = G1;

    fn mul(self, other: Fr) -> G1 {
        G1(self.0 * other.0)
    }
}

/// an additive cyclic group of prime order r
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct G2(groups::G2);

impl G2 {
    pub fn new(x: Fq2, y: Fq2, z: Fq2) -> Self {
        G2(groups::G2::new(x.0, y.0, z.0))
    }

    pub fn x(&self) -> Fq2 {
        Fq2(*self.0.x())
    }

    pub fn set_x(&mut self, x: Fq2) {
        *self.0.x_mut() = x.0
    }

    pub fn y(&self) -> Fq2 {
        Fq2(*self.0.y())
    }

    pub fn set_y(&mut self, y: Fq2) {
        *self.0.y_mut() = y.0
    }

    pub fn z(&self) -> Fq2 {
        Fq2(*self.0.z())
    }

    pub fn set_z(&mut self, z: Fq2) {
        *self.0.z_mut() = z.0
    }

    pub fn b() -> Fq2 {
        Fq2(G2Params::coeff_b())
    }
    /// get the element from a compressed form.
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General
    // Annex A  A.4.3
    pub fn from_compressed(bytes: &[u8]) -> Result<Self, CurveError> {
        if bytes.len() != 65 {
            return Err(CurveError::InvalidEncoding);
        }
        let sign = bytes[0];
        debug_assert!(sign == 2 || sign == 3);
        let x = Fq2::from_slice(&bytes[1..]).ok_or(CurveError::InvalidEncoding)?;
        let y_squared = (x * x * x) + Self::b();
        let mut y = y_squared.sqrt().ok_or(CurveError::NotMember)?;
        let is_even = sign & 1 == 0;

        if is_even != y.is_even() {
            y = -y;
        }

        AffineG2::new(x, y)
            .map_err(|_| CurveError::NotMember)
            .map(Into::into)
    }
    /// Serializes this element into compressed form.
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General 6.2.8
    pub fn to_compressed(self) -> [u8; 65] {
        let mut res = [0u8; 65];
        let ag1 = AffineG2::from_jacobian(self).unwrap();
        res[0] = 2; // compressed
        if !ag1.y().is_even() {
            res[0] |= 1;
        }
        res[1..].copy_from_slice(&ag1.x().to_slice());

        res
    }
    /// Serializes this element into uncompressed form.
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General, 6.2.8 , case 1
    pub fn to_slice(self) -> [u8; 128] {
        let mut res = [0u8; 128];
        let ag = AffineG2::from_jacobian(self).unwrap();
        res[..64].copy_from_slice(&ag.x().to_slice());
        res[64..].copy_from_slice(&ag.y().to_slice());

        res
    }
    // SM9 dentity-based cryptographic algorithms
    // Part 1: General, 6.2.8 , case 1
    pub fn from_slice(bytes: &[u8]) -> Result<Self, CurveError> {
        if bytes.len() != 128 {
            return Err(CurveError::InvalidEncoding);
        }

        let x = Fq2::from_slice(&bytes[..64]).ok_or(CurveError::InvalidEncoding)?;
        let y = Fq2::from_slice(&bytes[64..]).ok_or(CurveError::InvalidEncoding)?;

        AffineG2::new(x, y)
            .map_err(|_| CurveError::NotMember)
            .map(Into::into)
    }
}

impl Group for G2 {
    fn zero() -> Self {
        G2(groups::G2::zero())
    }
    fn one() -> Self {
        G2(groups::G2::one())
    }

    fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn normalize(&mut self) {
        let new = match self.0.to_affine() {
            Some(a) => a,
            None => return,
        };

        self.0 = new.to_jacobian();
    }
}

impl Add<G2> for G2 {
    type Output = G2;

    fn add(self, other: G2) -> G2 {
        G2(self.0 + other.0)
    }
}

impl Sub<G2> for G2 {
    type Output = G2;

    fn sub(self, other: G2) -> G2 {
        G2(self.0 - other.0)
    }
}

impl Neg for G2 {
    type Output = G2;

    fn neg(self) -> G2 {
        G2(-self.0)
    }
}

impl Mul<Fr> for G2 {
    type Output = G2;

    fn mul(self, other: Fr) -> G2 {
        G2(self.0 * other.0)
    }
}

/// a multiplicative cyclic group of prime order r
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct Gt(fields::Fq12);

impl Gt {
    pub fn one() -> Self {
        Gt(fields::Fq12::one())
    }
    pub fn pow(&self, exp: Fr) -> Self {
        Gt(self.0.pow(exp.0))
    }
    pub fn inverse(&self) -> Option<Self> {
        self.0.inverse().map(Gt)
    }
    /// Converts an element into a byte representation in
    /// big-endian byte order.
    pub fn to_slice(self) -> [u8; 384] {
        self.0.to_slice()
    }
}

impl Mul<Gt> for Gt {
    type Output = Gt;

    fn mul(self, other: Gt) -> Gt {
        Gt(self.0 * other.0)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
pub struct AffineG1(groups::AffineG1);

impl AffineG1 {
    pub fn new(x: Fq, y: Fq) -> Result<Self, GroupError> {
        Ok(AffineG1(groups::AffineG1::new(x.0, y.0)?))
    }

    pub fn x(&self) -> Fq {
        Fq(*self.0.x())
    }

    pub fn set_x(&mut self, x: Fq) {
        *self.0.x_mut() = x.0
    }

    pub fn y(&self) -> Fq {
        Fq(*self.0.y())
    }

    pub fn set_y(&mut self, y: Fq) {
        *self.0.y_mut() = y.0
    }

    pub fn from_jacobian(g1: G1) -> Option<Self> {
        g1.0.to_affine().map(AffineG1)
    }
}

impl From<AffineG1> for G1 {
    fn from(affine: AffineG1) -> Self {
        G1(affine.0.to_jacobian())
    }
}

#[derive(Copy, Clone, PartialEq, Eq)]
#[repr(C)]
pub struct AffineG2(groups::AffineG2);

impl AffineG2 {
    pub fn new(x: Fq2, y: Fq2) -> Result<Self, GroupError> {
        Ok(AffineG2(groups::AffineG2::new(x.0, y.0)?))
    }

    pub fn x(&self) -> Fq2 {
        Fq2(*self.0.x())
    }

    pub fn set_x(&mut self, x: Fq2) {
        *self.0.x_mut() = x.0
    }

    pub fn y(&self) -> Fq2 {
        Fq2(*self.0.y())
    }

    pub fn set_y(&mut self, y: Fq2) {
        *self.0.y_mut() = y.0
    }

    pub fn from_jacobian(g2: G2) -> Option<Self> {
        g2.0.to_affine().map(AffineG2)
    }
}

impl From<AffineG2> for G2 {
    fn from(affine: AffineG2) -> Self {
        G2(affine.0.to_jacobian())
    }
}
impl From<G2> for G2Prepared {
    fn from(g2: G2) -> G2Prepared {
        let mut g = g2;
        g.normalize();
        G2Prepared::from(g.0)
    }
}
impl G2Prepared {
    pub fn pairing(&self, g1: &G1) -> Gt {
        let mut g = *g1;
        g.normalize();
        Gt(self
            .miller_loop(&g.0)
            .final_exp()
            .expect("miller loop cannot produce zero"))
    }
}
/// compute R-ate Pairing G2 x G1 -> GT
pub fn pairing(p: G1, q: G2) -> Gt {
    Gt(pairings::pairing(&p.0, &q.0))
}
/// compute R-ate Pairing G2 x G1 -> GT
pub fn fast_pairing(mut p: G1, mut q: G2) -> Gt {
    p.normalize();
    q.normalize();
    Gt(pairings::fast_pairing(&p.0, &q.0))
}

/************************************************************************************************ */
#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields;
    use hex_literal::hex;

    #[test]
    fn test_fq_from_to_slice() {
        let f1 = Fq::from_slice(&hex!(
            "A5702F05 CF131530 5E2D6EB6 4B0DEB92 3DB1A0BC F0CAFF90 523AC875 4AA69820"
        ))
        .unwrap();
        let s = f1.to_slice();
        let f2 = Fq::from_slice(&s).unwrap();
        assert_eq!(f1, f2);
    }
    #[test]
    fn test_fq2_from_to_slice() {
        let f1 = Fq2::new(
            Fq::from_slice(&hex!(
                "29DBA116 152D1F78 6CE843ED 24A3B573 414D2177 386A92DD 8F14D656 96EA5E32"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "9F64080B 3084F733 E48AFF4B 41B56501 1CE0711C 5E392CFB 0AB1B679 1B94C408"
            ))
            .unwrap(),
        );
        let s = f1.to_slice();
        let f2 = Fq2::from_slice(&s).unwrap();
        assert_eq!(f1, f2);
    }
    #[test]
    fn test_fq_to_big_endian() {
        let a = Fq::from_str("1").unwrap();
        let mut b = [0u8; 32];
        let mut c = [0u8; 32];
        c[31] = 1;
        a.to_big_endian(&mut b[..]).unwrap();
        assert_eq!(b, c);
    }

    #[test]
    fn lib_g1_mul() {
        let t2 = Fr::from_slice(&hex!(
            "291FE3CA C8F58AD2 DC462C8D 4D578A94 DAFD5624 DDC28E32 8D293668 8A86CF1A"
        ))
        .unwrap();
        let mut ds = G1::one() * t2;
        ds.normalize();

        let r = G1::new(
            Fq::from_slice(&hex!(
                "A5702F05 CF131530 5E2D6EB6 4B0DEB92 3DB1A0BC F0CAFF90 523AC875 4AA69820"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "78559A84 4411F982 5C109F5E E3F52D72 0DD01785 392A727B B1556952 B2B013D3"
            ))
            .unwrap(),
            Fq::one(),
        );
        assert_eq!(r, ds);
    }

    #[test]
    fn lib_g2_mul() {
        let ks = Fr::from_slice(&hex!(
            "000130E78459D78545CB54C587E02CF480CE0B66340F319F348A1D5B1F2DC5F4"
        ))
        .unwrap();
        let mut pub_s = G2::one() * ks;
        pub_s.normalize();

        let x = Fq2::new(
            Fq::from_slice(&hex!(
                "29DBA116 152D1F78 6CE843ED 24A3B573 414D2177 386A92DD 8F14D656 96EA5E32"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "9F64080B 3084F733 E48AFF4B 41B56501 1CE0711C 5E392CFB 0AB1B679 1B94C408"
            ))
            .unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&hex!(
                "41E00A53 DDA532DA 1A7CE027 B7A46F74 1006E85F 5CDFF073 0E75C05F B4E3216D"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "69850938 ABEA0112 B57329F4 47E3A0CB AD3E2FDB 1A77F335 E89E1408 D0EF1C25"
            ))
            .unwrap(),
        );
        let r = G2::new(x, y, Fq2::one());

        println!("pub_s {:?}", pub_s);
        assert_eq!(r, pub_s);
    }

    #[test]
    fn lib_pairing_test() {
        let ks = Fr::from_slice(&hex!(
            "000130E78459D78545CB54C587E02CF480CE0B66340F319F348A1D5B1F2DC5F4"
        ))
        .unwrap();
        let pub_s = G2::one() * ks;
        let g = pairing(G1::one(), pub_s);
        // println!(" {:#?}", g);
        let r0 = Fq2::new(
            Fq::from_slice(&hex!(
                "AAB9F06A 4EEBA432 3A7833DB 202E4E35 639D93FA 3305AF73 F0F071D7 D284FCFB"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "84B87422 330D7936 EABA1109 FA5A7A71 81EE16F2 438B0AEB 2F38FD5F 7554E57A"
            ))
            .unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&hex!(
                "4C744E69 C4A2E1C8 ED72F796 D151A17C E2325B94 3260FC46 0B9F73CB 57C9014B"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "B3129A75 D31D1719 4675A1BC 56947920 898FBF39 0A5BF5D9 31CE6CBB 3340F66D"
            ))
            .unwrap(),
        );
        let a = fields::Fq4::new(r0.0, r1.0);
        let r0 = Fq2::new(
            Fq::from_slice(&hex!(
                "93634F44 FA13AF76 169F3CC8 FBEA880A DAFF8475 D5FD28A7 5DEB83C4 4362B439"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "1604A3FC FA9783E6 67CE9FCB 1062C2A5 C6685C31 6DDA62DE 0548BAA6 BA30038B"
            ))
            .unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&hex!(
                "5A1AE172 102EFD95 DF7338DB C577C66D 8D6C15E0 A0158C75 07228EFB 078F42A6"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "67E0E0C2 EED7A699 3DCE28FE 9AA2EF56 83430786 0839677F 96685F2B 44D0911F"
            ))
            .unwrap(),
        );
        let b = fields::Fq4::new(r0.0, r1.0);
        let r0 = Fq2::new(
            Fq::from_slice(&hex!(
                "A01F2C8B EE817696 09462C69 C96AA923 FD863E20 9D3CE26D D889B55E 2E3873DB"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "38BFFE40 A22D529A 0C66124B 2C308DAC 92299126 56F62B4F ACFCED40 8E02380F"
            ))
            .unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&hex!(
                "28B3404A 61908F5D 6198815C 99AF1990 C8AF3865 5930058C 28C21BB5 39CE0000"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "4E378FB5 561CD066 8F906B73 1AC58FEE 25738EDF 09CADC7A 29C0ABC0 177AEA6D"
            ))
            .unwrap(),
        );
        let c = fields::Fq4::new(r0.0, r1.0);
        let r = fields::Fq12::new(a, b, c);
        assert_eq!(r, g.0);
        let g1 = fast_pairing(G1::one(), pub_s);
        assert_eq!(r, g1.0);
    }

    #[test]
    fn lib_gt_pow_test() {
        let ks = Fr::from_slice(&hex!(
            "000130E78459D78545CB54C587E02CF480CE0B66340F319F348A1D5B1F2DC5F4"
        ))
        .unwrap();
        let r = Fr::from_slice(&hex!(
            "00033C86 16B06704 813203DF D0096502 2ED15975 C662337A ED648835 DC4B1CBE"
        ))
        .unwrap();
        let pub_s = G2::one() * ks;
        let g = pairing(G1::one(), pub_s).pow(r);
        // println!(" {:#?}", g);
        let r0 = Fq2::new(
            Fq::from_slice(&hex!(
                "1F96B08E 97997363 91131470 5BFB9A9D BB97F755 53EC90FB B2DDAE53 C8F68E42"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "6A814AAF 475F128A EF43A128 E37F8015 4AE6CB92 CAD7D150 1BAE30F7 50B3A9BD"
            ))
            .unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&hex!(
                "898D6084 8026B7EF B8FCC1B2 442ECF07 95F8A81C EE99A624 8F294C82 C90D26BD"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "44643CEA D40F0965 F28E1CD2 895C3D11 8E4F65C9 A0E3E741 B6DD52C0 EE2D25F5"
            ))
            .unwrap(),
        );
        let a = fields::Fq4::new(r0.0, r1.0);
        let r0 = Fq2::new(
            Fq::from_slice(&hex!(
                "0656FCB6 63D24731 E8029218 8A2471B8 B68AA993 89926849 9D23C897 55A1A897"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "4F8624EB 435B838C CA77B2D0 347E65D5 E4696441 2A096F41 50D8C5ED E5440DDF"
            ))
            .unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&hex!(
                "3F012DB0 4BA59FE8 8DB88932 1CC2373D 4C0C35E8 4F7AB1FF 33679BCA 575D6765"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "A543D256 09AE9439 20679194 ED30328B B33FD156 60BDE485 C6B79A7B 32B01398"
            ))
            .unwrap(),
        );
        let b = fields::Fq4::new(r0.0, r1.0);
        let r0 = Fq2::new(
            Fq::from_slice(&hex!(
                "8EAF5D17 9A1836B3 59A9D1D9 BFC19F2E FCDB8293 28620962 BD3FDF15 F2567F58"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "30DADC5C D9E207AE E32209F6 C3CA3EC0 D800A1A4 2D33C731 53DED47C 70A39D2E"
            ))
            .unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&hex!(
                "815AEBA2 17AD502D A0F48704 CC73CABB 3C06209B D87142E1 4CBD99E8 BCA1680F"
            ))
            .unwrap(),
            Fq::from_slice(&hex!(
                "81377B8F DBC2839B 4FA2D0E0 F8AA6853 BBBE9E9C 4099608F 8612C607 8ACD7563"
            ))
            .unwrap(),
        );
        let c = fields::Fq4::new(r0.0, r1.0);
        let w = fields::Fq12::new(a, b, c);
        assert_eq!(w, g.0)
    }
}
