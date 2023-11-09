# sm9_core[![Build status](https://api.travis-ci.org/zcash/bn.svg)](https://travis-ci.org/zcash/bn)

This is a [pairing cryptography](https://en.wikipedia.org/wiki/Pairing-based_cryptography) library written in pure Rust. It makes use of the Barreto-Naehrig (BN) curve construction from "SM9 identity-based  cryptographic algorithms" to provide two cyclic groups **G `<sub>`1 `</sub>`** and **G `<sub>`2 `</sub>`**, with an R-ate pairing:

*e: G `<sub>`1 `</sub>` × G `<sub>`2 `</sub>` → G `<sub>`T `</sub>`*

## Security warnings

This library, like other pairing cryptography libraries implementing this construction, is not resistant to side-channel attacks.

## API

* `Fr` is an element of F `<sub>`r `</sub>`
* `G1` is a point on the BN curve E/Fq : y `<sup>`2 `</sup>` = x `<sup>`3 `</sup>` + b
* `G2` is a point on the twisted BN curve E'/Fq2 : y `<sup>`2 `</sup>` = x `<sup>`3 `</sup>` + b/xi
* `Gt` is a group element (written multiplicatively) obtained with the `pairing` function over `G1` and `G2`.

### Examples

(See `lib.rs` for the full example.)

```rust
        let ks = Fr::from_slice(&hex!("000130E78459D78545CB54C587E02CF480CE0B66340F319F348A1D5B1F2DC5F4")).unwrap();
        let r = Fr::from_slice(&hex!("00033C86 16B06704 813203DF D0096502 2ED15975 C662337A ED648835 DC4B1CBE")).unwrap();
        let pub_s = G2::one() * ks;
        let g = pairing(G1::one(), pub_s).pow(r);
        println!(" {:#?}", g);

```

## License

Licensed under either of

* MIT license, (http://opensource.org/licenses/MIT)
* Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)

at your option.

Copyright 2023 [John-Yu](https://github.com/John-Yu). 

### Thanks

The fields and groups algorithms come from [zcash - bn](https://github.com/zcash-hackworks/bn), and pairing algorithms come from [GmSSl](https://github.com/guanzhi/GmSSL). 

Thanks for them.

### Authors

* [John-Yu](https://github.com/John-Yu)

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
