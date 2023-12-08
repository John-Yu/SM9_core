use criterion::{black_box, criterion_group, criterion_main, Criterion};
use rand::prelude::*;
use sm9_core::*;

fn criterion_benchmark(c: &mut Criterion) {
    let rng = &mut thread_rng();
    let a = Fr::random(rng);
    let b = Fr::random(rng);
    let g = G1::one() * a;
    let h = G2::one() * b;
    let g2_precomputed = G2Prepared::from(h);

    c.bench_function("pairing", |b| {
        b.iter(|| pairing(black_box(g), black_box(h)))
    });
    c.bench_function("fast_pairing", |b| {
        b.iter(|| fast_pairing(black_box(g), black_box(h)))
    });
    c.bench_function("precomputed_pairing", |b| {
        b.iter(|| g2_precomputed.pairing(black_box(&g)))
    });
}

criterion_group! {
    name = benches;
    // This can be any expression that returns a `Criterion` object.
    config = Criterion::default().significance_level(0.1).sample_size(50);
    targets = criterion_benchmark
}
criterion_main!(benches);
