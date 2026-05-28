use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;
use tempfile;

fn bench_featurecounts(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-featurecounts");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bam = manifest.join("tests/golden/small.bam");
    let gtf = manifest.join("tests/golden/adv.gtf");
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("counts.txt");
    c.bench_function("rsomics-featurecounts golden", |b| {
        b.iter(|| {
            let status = Command::new(black_box(bin))
                .args([
                    "-a",
                    gtf.to_str().unwrap(),
                    "-o",
                    out.to_str().unwrap(),
                    bam.to_str().unwrap(),
                ])
                .status()
                .unwrap();
            assert!(status.success());
        });
    });
}

criterion_group!(benches, bench_featurecounts);
criterion_main!(benches);
