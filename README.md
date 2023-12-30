# sm9_core

This is a [pairing cryptography](https://en.wikipedia.org/wiki/Pairing-based_cryptography) library written in pure Rust. 
It makes use of the Barreto-Naehrig (BN) curve construction from "[SM9](https://en.wikipedia.org/wiki/SM9_(cryptography_standard)) identity-based  cryptographic algorithms" as well as [ISO/IEC 11770](https://www.iso.org/standard/82709.html) to provide two cyclic groups **G <sub>1 </sub>** and **G <sub>2 </sub>**, with an R-ate pairing:

*e: G <sub>1 </sub> × G <sub>2 </sub> → G <sub>T </sub>*

## Security warnings

This library, like other pairing cryptography libraries implementing this construction, is not resistant to side-channel attacks.

## Usage

Add the `sm9_core` crate to your dependencies in `Cargo.toml`

```toml
[dependencies]
sm9_core = "0.3.5"
```

## API

* `Fr` is an element of F <sub>r </sub>
* `G1` is a point on the BN curve E/Fq : y <sup>2 </sup> = x <sup>3 </sup> + b
* `G2` is a point on the twisted BN curve E'/Fq2 : y <sup>2 </sup> = x <sup>3 </sup> + b x i
* `Gt` is a group element (written multiplicatively)
* `pairing()` is a  API to compute R-ate Pairing G2 x G1 -> GT
* `fast_pairing()` is another  API to compute R-ate Pairing G2 x G1 -> GT

### Examples

(See `integration_test.rs` for the full example.)

```rust
use hex_literal::hex;
use sm9_core::*;

let ks = Fr::from_slice(&hex!("000130E78459D78545CB54C587E02CF480CE0B66340F319F348A1D5B1F2DC5F4")).unwrap();
let r = Fr::from_slice(&hex!("00033C86 16B06704 813203DF D0096502 2ED15975 C662337A ED648835 DC4B1CBE")).unwrap();
let pub_s = G2::one() * ks;
let g = pairing(G1::one(), pub_s).pow(r);
println!(" {:#?}", g);
let r1 = g.to_slice();
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
// test  fast_pairing
let g1 = fast_pairing(G1::one(), pub_s).pow(r);
let r1 = g1.to_slice();
assert_eq!(r0, r1);
```

## License

Licensed under either of

* [MIT license](http://opensource.org/licenses/MIT)
* [Apache License, Version 2.0](http://www.apache.org/licenses/LICENSE-2.0)

at your option.

Copyright 2023 [John-Yu](https://github.com/John-Yu).

### Thanks

The fields and groups algorithms come from [zcash-bn](https://github.com/zcash-hackworks/bn) and [bls12_381](https://github.com/zkcrypto/bls12_381), and pairing algorithms come from [GmSSL](https://github.com/guanzhi/GmSSL).
The fast_pairing algorithms come from [MIRACL](https://github.com/miracl/MIRACL),  it is 23% faster than pairing().

Thanks to them.

### Benchmark

(OS: windows11, CPU:  i7-8700K 3.70GHz, See `my_benchmark.rs` for the details)

| function | times |
|:-:|:-:|
| pairing | time:   [1.7508 ms 1.7713 ms 1.7957 ms] |
| fast_pairing  |time:   [1.3401 ms 1.3478 ms 1.3562 ms] |
| precomputed_pairing |time:   [1.1537 ms 1.1599 ms 1.1680 ms] |

### Authors

* [John-Yu](https://github.com/John-Yu)

### Contribution

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the Apache-2.0
license, shall be dual licensed as above, without any additional terms or
conditions.
