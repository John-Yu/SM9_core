// import commonly used items from the prelude:
//use hex_literal::hex;
use rand::prelude::*;
use sm9_core::*;

//A single round protocol is possible through the use of a bilinear pairing: given Alice's public key aP1 and Bob's public key bP2,
//Carol can compute the shared secret with her private key c by e(aP1, bP2)^c.
// Joux's key agreement protocol
#[test]
fn test_joux() {
    // If we want to be a bit more explicit (and a little more efficient) we can
    // make a handle to the thread-local generator:
    let rng = &mut thread_rng();

    // Generate private keys
    let alice_sk = Fr::random(rng);
    let bob_sk = Fr::random(rng);
    let carol_sk = Fr::random(rng);

    // Generate public keys in G1 and G2
    let (alice_pk1, alice_pk2) = (G1::one() * alice_sk, G2::one() * alice_sk);
    let (bob_pk1, bob_pk2) = (G1::one() * bob_sk, G2::one() * bob_sk);
    let (carol_pk1, carol_pk2) = (G1::one() * carol_sk, G2::one() * carol_sk);

    // Each party computes the shared secret
    let alice_ss = pairing(bob_pk1, carol_pk2).pow(alice_sk);
    let bob_ss = pairing(carol_pk1, alice_pk2).pow(bob_sk);
    let carol_ss = pairing(alice_pk1, bob_pk2).pow(carol_sk);

    // they are equal
    assert!(alice_ss == bob_ss && bob_ss == carol_ss);
}

// This is an example of three-party Diffie-Hellman key exchange
// Requires two rounds
#[test]
fn test_dh() {
    let rng = &mut thread_rng();

    // Generate private keys
    let alice_sk = Fr::random(rng);
    let bob_sk = Fr::random(rng);
    let carol_sk = Fr::random(rng);
    // Construct public keys
    let alice_pk = G1::one() * alice_sk;
    let bob_pk = G1::one() * bob_sk;
    let carol_pk = G1::one() * carol_sk;

    // Round one:
    let alice_dh_1 = bob_pk * carol_sk;
    let bob_dh_1 = carol_pk * alice_sk;
    let carol_dh_1 = alice_pk * bob_sk;

    // Round two:
    let alice_dh_2 = alice_dh_1 * alice_sk;
    let bob_dh_2 = bob_dh_1 * bob_sk;
    let carol_dh_2 = carol_dh_1 * carol_sk;

    // All parties should arrive to the same shared secret
    assert!(alice_dh_2 == bob_dh_2 && bob_dh_2 == carol_dh_2);
}
#[test]
fn test_bilinearity() {
    let rng = &mut thread_rng();
    let a = Fr::random(rng);
    let b = Fr::random(rng);
    let c = a * b;

    let g = G1::one() * a;
    let h = G2::one() * b;
    let p = pairing(g, h);

    assert_eq!(p, pairing(G1::one(), G2::one()).pow(c));
    assert_eq!(p, fast_pairing(G1::one(), G2::one()).pow(c));

    let g2_precomputed = G2Prepared::from(h);
    assert_eq!(p, g2_precomputed.pairing(&G1::one()).pow(a));
    assert_eq!(p, g2_precomputed.pairing(&g));
}

#[test]
fn test_unitary() {
    let g = G1::one();
    let h = G2::one();
    let q = pairing(g, -h);
    let r = pairing(-g, h);

    assert_eq!(q, r);
}

#[test]
fn test_gt_to_slice() {
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
    let r1 = g.to_slice();
    //// println!(" {:#?}", g);
    let r0 = hex!(
        "81377B8F DBC2839B 4FA2D0E0 F8AA6853 BBBE9E9C 4099608F 8612C607 8ACD7563"
        "815AEBA2 17AD502D A0F48704 CC73CABB 3C06209B D87142E1 4CBD99E8 BCA1680F"
        "30DADC5C D9E207AE E32209F6 C3CA3EC0 D800A1A4 2D33C731 53DED47C 70A39D2E"
        "8EAF5D17 9A1836B3 59A9D1D9 BFC19F2E FCDB8293 28620962 BD3FDF15 F2567F58"
        "A543D256 09AE9439 20679194 ED30328B B33FD156 60BDE485 C6B79A7B 32B01398"
        "3F012DB0 4BA59FE8 8DB88932 1CC2373D 4C0C35E8 4F7AB1FF 33679BCA 575D6765"
        "4F8624EB 435B838C CA77B2D0 347E65D5 E4696441 2A096F41 50D8C5ED E5440DDF"
        "0656FCB6 63D24731 E8029218 8A2471B8 B68AA993 89926849 9D23C897 55A1A897"
        "44643CEA D40F0965 F28E1CD2 895C3D11 8E4F65C9 A0E3E741 B6DD52C0 EE2D25F5"
        "898D6084 8026B7EF B8FCC1B2 442ECF07 95F8A81C EE99A624 8F294C82 C90D26BD"
        "6A814AAF 475F128A EF43A128 E37F8015 4AE6CB92 CAD7D150 1BAE30F7 50B3A9BD"
        "1F96B08E 97997363 91131470 5BFB9A9D BB97F755 53EC90FB B2DDAE53 C8F68E42"
    );
    assert_eq!(r0, r1);
    let g1 = fast_pairing(G1::one(), pub_s).pow(r);
    let r1 = g1.to_slice();
    assert_eq!(r0, r1);
    // to_slice, from_slice
    let sx = pub_s.x().to_slice();
    let sy = pub_s.y().to_slice();
    let sz = pub_s.z().to_slice();
    let x = Fq2::from_slice(&sx).unwrap();
    let y = Fq2::from_slice(&sy).unwrap();
    let z = Fq2::from_slice(&sz).unwrap();
    let mut g2 = G2::zero();
    g2.set_x(x);
    g2.set_y(y);
    g2.set_z(z);
    assert_eq!(pub_s, g2);
}

#[test]
fn test_g1_compress() {
    let rng = &mut thread_rng();
    let a = Fr::random(rng);
    let mut g = G1::one() * a;
    g.normalize();
    let c0 = g.to_compressed();
    let g1 = G1::from_compressed(&c0).unwrap();

    assert_eq!(g, g1);
}

#[test]
fn test_g1_uncompress() {
    let rng = &mut thread_rng();
    let a = Fr::random(rng);
    let mut g = G1::one() * a;
    g.normalize();
    let c0 = g.to_slice();
    let g1 = G1::from_slice(&c0).unwrap();

    assert_eq!(g, g1);
}

#[test]
fn test_g2_compress() {
    let rng = &mut thread_rng();
    let a = Fr::random(rng);
    let mut g = G2::one() * a;
    g.normalize();

    let c0 = g.to_compressed();
    let g2 = G2::from_compressed(&c0).unwrap();

    assert_eq!(g, g2);
}

#[test]
fn test_g2_uncompress() {
    let rng = &mut thread_rng();
    let a = Fr::random(rng);
    let mut g = G2::one() * a;
    g.normalize();
    let c0 = g.to_slice();
    let g1 = G2::from_slice(&c0).unwrap();

    assert_eq!(g, g1);
}
#[test]
fn test_fr_from_str() {
    let a = Fr::from_str("1024").unwrap();
    let b = Fr::try_from(hex!("04 00").as_ref()).unwrap();
    assert_eq!(a, b);
    assert_eq!(a + b, b + a);
}
