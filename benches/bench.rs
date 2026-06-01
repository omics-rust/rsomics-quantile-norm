use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

use criterion::{Criterion, criterion_group, criterion_main};

fn bench_quantile_norm(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-quantile-norm");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let counts = manifest.join("tests/golden/counts.tsv");
    c.bench_function("rsomics-quantile-norm golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([counts.to_str().unwrap(), "-o", "/dev/null"])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_quantile_norm);
criterion_main!(benches);
