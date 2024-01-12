#![allow(dead_code)]

use alloc::vec::Vec;
use core::ops::Neg;

use crate::fields::{FieldElement, Fq, Fq12, Fq2, Fq4};
use crate::groups::{GroupElement, G1, G2};
use crate::u256::U256;

// abits = "00100000000000000000000000000000000000010000101100020200101000020";
const SM9_LOOP_COUNT: [u8; 65] = [
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 2, 0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2,
    0,
];
const SM9_S: u128 = 0x600000000058F98A;
// 6s+2, where s=0x600000000058F98A
const SM9_LOOP_N: u128 = 0x2400000000215D93E;
// a2 = 0xd8000000019062ed 0000b98b0cb27659
const SM9_A2: u128 = 0xd8000000019062ed0000b98b0cb27659;
// a3 = 0x2 400000000215d941
const SM9_A3: u128 = 0x2400000000215d941;
const SM9_NINE: u128 = 0x9;

impl Fq12 {
    /// Raises a value to the power of exp, using exponentiation by squaring.
    #[inline]
    fn pow(&self, mut exp: u128) -> Fq12 {
        if exp == 0 {
            return Fq12::one();
        }
        let mut base = *self;

        while exp & 1 == 0 {
            base = base.squared();
            exp >>= 1;
        }

        if exp == 1 {
            return base;
        }

        let mut acc = base;
        while exp > 1 {
            exp >>= 1;
            base = base.squared();
            if exp & 1 == 1 {
                acc *= &base;
            }
        }
        acc
    }
    fn final_exponentiation_first_chunk(&self) -> Option<Fq12> {
        let b = self.inverse()?;
        let a = self.frobenius_map(6);
        let c = a * b;
        let d = c.frobenius_map(2);

        Some(d * c)
    }

    fn final_exponentiation_last_chunk(&self) -> Fq12 {
        let a = self.pow(SM9_A3);
        let b = a.inverse().unwrap();
        let c = b.frobenius_map(1);
        let d = c * b;

        let e = d * b;
        let f = self.frobenius_map(1);
        let g = self * f;
        let h = g.pow(SM9_NINE);

        let i = e * h;
        let j = self.squared();
        let k = j.squared();
        let l = k * i;
        let m = f.squared();
        let n = d * m;
        let o = self.frobenius_map(2);
        let p = o * n;

        let q = p.pow(SM9_A2);
        let r = q * l;
        let s = self.frobenius_map(3);

        s * r
    }
    /// This performs a "final exponentiation" routine to convert the result
    /// of a Miller loop into an element of `Fq12`
    pub fn final_exponentiation(&self) -> Option<Fq12> {
        self.final_exponentiation_first_chunk()
            .map(|a| a.final_exponentiation_last_chunk())
    }

    //final_exp
    fn final_exp_last_chunk(&self) -> Fq12 {
        let mut t1 = self.pow(SM9_S).inverse().unwrap();
        let mut t0 = self.frobenius_map(1);
        let mut x0 = self.frobenius_map(2);
        let x1 = self.frobenius_map(6);
        let x3 = t1.frobenius_map(1);
        let mut x4 = t1;

        x0 *= self * t0;
        x0 = x0.frobenius_map(1);
        let x5 = t1.pow(SM9_S);
        t1 = x5.inverse().unwrap();
        x4 *= t1.frobenius_map(1).inverse().unwrap();
        let x2 = t1.frobenius_map(2);

        t0 = t1.pow(SM9_S).inverse().unwrap();
        t1 = t0.frobenius_map(1);
        t0 *= t1;
        t0 = t0.squared();
        t0 *= x4 * x5;
        t1 = x3 * x5;
        t1 *= t0;
        t0 *= x2;
        t1 = t1.squared();
        t1 *= t0;
        t1 = t1.squared();
        t0 = t1 * x1;
        t1 *= x0;
        t0 = t0.squared();

        t0 * t1
    }
    pub fn final_exp(&self) -> Option<Fq12> {
        self.final_exponentiation_first_chunk()
            .map(|a| a.final_exp_last_chunk())
    }
}

impl G2 {
    fn point_pi1(&self) -> Self {
        G2::new(
            self.x().unitary_inverse(),
            self.y().unitary_inverse(),
            self.z()
                .unitary_inverse()
                .scale(&Fq::new(*SM9_PI1).unwrap()),
        )
    }

    fn point_pi2(&self) -> Self {
        G2::new(self.x, self.y, self.z().scale(&Fq::new(*SM9_PI2).unwrap()))
    }

    fn eval_g_tangent(&self, q: &G1) -> (Fq12, Fq12) {
        let mut num = Fq12::zero();
        let mut den = Fq12::zero();

        let t0 = self.z().squared();
        let t1 = t0 * self.z();
        let b1 = t1 * self.y();

        let t2 = b1.scale(q.y());
        let a1 = -t2;

        let t1 = self.x().squared();
        let t0 = t0 * t1;
        let t0 = t0.scale(q.x());
        let t0 = t0.triple();
        let a4 = t0.div2();

        let t1 = t1 * self.x();
        let t1 = t1.triple();
        let t1 = t1.div2();
        let t0 = self.y().squared();
        let a0 = t0 - t1;

        num.c0 = Fq4::new(a0, a1);
        num.c2 = Fq4::new(a4, Fq2::zero());
        den.c0 = Fq4::new(Fq2::zero(), b1);

        (num, den)
    }

    fn eval_g_line(&self, p: &G2, q: &G1) -> (Fq12, Fq12) {
        let mut num = Fq12::zero();
        let mut den = Fq12::zero();

        let t0 = p.z().squared();
        let t1 = t0 * self.x();
        let t0 = t0 * p.z();

        let t2 = self.z().squared();
        let t3 = t2 * p.x();
        let t2 = t2 * self.z();

        let t2 = t2 * p.y();
        let t1 = t1 - t3;
        let t1 = t1 * self.z();

        let t1 = t1 * p.z();
        let t4 = t1 * t0;
        let b1 = t4;

        let t1 = t1 * p.y();
        let t3 = t0 * self.y();
        let t3 = t3 - t2;
        let t0 = t0 * t3;
        let t0 = t0.scale(q.x());
        let a4 = t0;

        let t3 = t3 * p.x();
        let t3 = t3 * p.z();
        let t1 = t1 - t3;
        let a0 = t1;

        let t2 = t4.scale(q.y());
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
            f_num *= g_num;
            f_den *= g_den;

            t = t.double();

            if *i == 1 {
                (g_num, g_den) = t.eval_g_line(self, p);
                f_num *= g_num;
                f_den *= g_den;
                t = t + *self;
            } else if *i == 2 {
                (g_num, g_den) = t.eval_g_line(&q1, p);
                f_num *= g_num;
                f_den *= g_den;
                t = t + q1;
            }
        }
        q1 = self.point_pi1();
        let q2 = -self.point_pi2();

        (g_num, g_den) = t.eval_g_line(&q1, p);
        f_num *= g_num;
        f_den *= g_den;
        t = t + q1;

        (g_num, g_den) = t.eval_g_line(&q2, p);
        f_num *= g_num;
        f_den *= g_den;

        f_den = f_den.inverse().unwrap();
        f_num * f_den
    }
    // Fast multiplication of A by q (for Trace-Zero group members only)
    // Calculate q*P. P(X,Y) -> P(X^p,Y^p))
    fn q_power_frobenius(&self, f: &Fq2) -> Option<Self> {
        let r = f.inverse()?;
        let w = r.squared();

        Some(G2::new(
            self.x().unitary_inverse() * w,
            self.y().unitary_inverse() * w * r,
            self.z().unitary_inverse(),
        ))
    }
    fn g_line(&mut self, g2: &G2) -> (Fq2, Fq2, Fq2) {
        let lam = self.z().squared();
        let c2 = self.y() - (lam * self.z()) * g2.y();
        *self = *self + *g2;
        let c0 = self.z;
        let c1 = -c2 * g2.x() - c0 * g2.y();

        (c0, c1, c2)
    }
    fn g_tangent(&mut self) -> (Fq2, Fq2, Fq2) {
        let lam = self.x().squared().triple();
        let extra = self.y().squared().double();
        let zz = self.z().squared();
        let c1 = lam * self.x() - extra;
        let c2 = -(zz * lam);
        *self = self.double();
        let c0 = self.z() * zz;

        (c0, c1, c2)
    }
}
//
/// R-ate Pairing G2 x G1 -> GT
//
// P is a point of order r in G1. Q(x,y) is a point of order r in G2.
// Note that Q is a point on the sextic twist of the curve over Fp^2, P(x,y) is a point on the
// curve over the base field Fp
// the curve is y^2 = x^3 + 5
pub fn pairing(p: &G1, q: &G2) -> Fq12 {
    match (p.to_affine(), q.to_affine()) {
        (None, _) | (_, None) => Fq12::one(),
        (Some(p), Some(q)) => q
            .to_jacobian()
            .miller_loop(&p.to_jacobian())
            .final_exponentiation()
            .expect("miller loop cannot produce zero"),
    }
}

#[derive(Clone, Debug)]
/// This structure contains cached computations pertaining to a G2
/// element as part of the pairing function (specifically, the Miller loop)
// and so should be computed whenever a G2 element is being used in
// multiple pairings or is otherwise known in advance.
pub struct G2Prepared {
    coeffs: Vec<(Fq2, Fq2, Fq2)>,
}

impl From<G2> for G2Prepared {
    fn from(g2: G2) -> G2Prepared {
        let mut coeffs: Vec<(Fq2, Fq2, Fq2)> = Vec::new();
        let mut p = g2;

        for i in (0..65).rev() {
            let coeff = p.g_tangent();
            coeffs.push(coeff);
            if bit(SM9_LOOP_N, i) {
                let coeff = p.g_line(&g2);
                coeffs.push(coeff);
            }
        }
        let frob = Fq2::new(Fq::new(*SM9_PI1).unwrap(), Fq::zero());
        let mut ka = g2.q_power_frobenius(&frob).unwrap();
        let coeff = p.g_line(&ka);
        coeffs.push(coeff);

        ka = -ka.q_power_frobenius(&frob).unwrap();
        let coeff = p.g_line(&ka);
        coeffs.push(coeff);

        G2Prepared { coeffs }
    }
}
impl G2Prepared {
    fn get_fq12(&self, c: &(Fq2, Fq2, Fq2), t1: &Fq2, x: &Fq) -> Fq12 {
        Fq12 {
            c0: Fq4::new(c.0 * t1, c.1),
            c1: Fq4::zero(),
            c2: Fq4::new(Fq2::zero(), c.2.scale(x)),
        }
    }
    pub fn miller_loop(&self, g1: &G1) -> Fq12 {
        let mut f = Fq12::one();
        let t1 = Fq2::new(g1.y, Fq::zero()).mul_by_nonresidue();
        let mut idx = 0;

        for i in (0..65).rev() {
            let c = &self.coeffs[idx];
            idx += 1;
            f = f.squared();
            f = f.mul_015(&self.get_fq12(c, &t1, g1.x()));

            if bit(SM9_LOOP_N, i) {
                let c = &self.coeffs[idx];
                idx += 1;
                f = f.mul_015(&self.get_fq12(c, &t1, g1.x()));
            }
        }

        for _ in 0..2 {
            let c = &self.coeffs[idx];
            idx += 1;
            f = f.mul_015(&self.get_fq12(c, &t1, g1.x()));
        }

        f
    }
}

pub fn fast_pairing(g1: &G1, g2: &G2) -> Fq12 {
    let g2p = G2Prepared::from(*g2);
    g2p.miller_loop(g1)
        .final_exp()
        .expect("miller loop cannot produce zero")
}

#[inline]
fn bit(n: u128, pos: usize) -> bool {
    n & (1 << pos) != 0
}

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

/* ************************************************************************************************ */
#[cfg(test)]
mod tests {
    use super::*;
    use hex_literal::hex;

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
            Fq::from_slice(&G1_RAX).unwrap(),
            Fq::from_slice(&G1_RAY).unwrap(),
            Fq::one(),
        );

        let x = Fq2::new(
            Fq::from_slice(&G2_DEBXX).unwrap(),
            Fq::from_slice(&G2_DEBXY).unwrap(),
        );
        let y = Fq2::new(
            Fq::from_slice(&G2_DEBYX).unwrap(),
            Fq::from_slice(&G2_DEBYY).unwrap(),
        );
        let p2 = G2::new(x, y, Fq2::one());

        let r0 = Fq2::new(
            Fq::from_slice(&PAIR2A00).unwrap(),
            Fq::from_slice(&PAIR2A01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&PAIR2A10).unwrap(),
            Fq::from_slice(&PAIR2A11).unwrap(),
        );
        let a = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&PAIR2B00).unwrap(),
            Fq::from_slice(&PAIR2B01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&PAIR2B10).unwrap(),
            Fq::from_slice(&PAIR2B11).unwrap(),
        );
        let b = Fq4::new(r0, r1);
        let r0 = Fq2::new(
            Fq::from_slice(&PAIR2C00).unwrap(),
            Fq::from_slice(&PAIR2C01).unwrap(),
        );
        let r1 = Fq2::new(
            Fq::from_slice(&PAIR2C10).unwrap(),
            Fq::from_slice(&PAIR2C11).unwrap(),
        );
        let c = Fq4::new(r0, r1);
        let r = Fq12::new(a, b, c);

        let alice = pairing(&p1, &p2);
        let bob = fast_pairing(&p1, &p2);

        assert_eq!(alice, r);
        assert_eq!(bob, r);
    }
    #[test]
    fn test_fast_pairing() {
        let g1 = G1::one();
        let g2 = G2::one();
        fast_pairing(&g1, &g2);
    }
}
