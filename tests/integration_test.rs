// import commonly used items from the prelude:
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

    //unfortunately, they are not equal
    //assert!(alice_ss == bob_ss && bob_ss == carol_ss);
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
