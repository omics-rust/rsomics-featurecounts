use std::path::Path;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_featurecounts::{CountOpts, count_reads, load_exons};

const GOLDEN_DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden");

fn bench_count_reads_small(c: &mut Criterion) {
    let bam = Path::new(GOLDEN_DIR).join("small.bam");
    let gff = Path::new(GOLDEN_DIR).join("small.gff");
    let opts = CountOpts::default();
    // Load exons once outside the timed loop — index build is the hot path.
    c.bench_function("count_reads/small_gff", |b| {
        b.iter(|| {
            let exons = load_exons(&gff, &opts).unwrap();
            let _ = count_reads(&bam, &exons, &opts).unwrap();
        });
    });
}

fn bench_count_reads_adv(c: &mut Criterion) {
    let bam = Path::new(GOLDEN_DIR).join("adv.bam");
    let gtf = Path::new(GOLDEN_DIR).join("adv.gtf");
    let opts = CountOpts::default();
    c.bench_function("count_reads/adv_gtf", |b| {
        b.iter(|| {
            let exons = load_exons(&gtf, &opts).unwrap();
            let _ = count_reads(&bam, &exons, &opts).unwrap();
        });
    });
}

// Tier-3 bench: run only when FEATURECOUNTS_BENCH_BAM / FEATURECOUNTS_BENCH_GTF are set.
// Use: FEATURECOUNTS_BENCH_BAM=/path/to/large.bam FEATURECOUNTS_BENCH_GTF=/path/to/annot.gtf
//      cargo bench --bench featurecounts -- tier3
fn bench_count_reads_tier3(c: &mut Criterion) {
    let Some(bam_env) = std::env::var("FEATURECOUNTS_BENCH_BAM").ok() else {
        return;
    };
    let Some(gtf_env) = std::env::var("FEATURECOUNTS_BENCH_GTF").ok() else {
        return;
    };
    let bam = Path::new(&bam_env);
    let gtf = Path::new(&gtf_env);
    if !bam.exists() || !gtf.exists() {
        return;
    }
    let opts = CountOpts::default();
    let exons = load_exons(gtf, &opts).unwrap();
    c.bench_function("count_reads/tier3", |b| {
        b.iter(|| {
            let _ = count_reads(bam, &exons, &opts).unwrap();
        });
    });
}

criterion_group!(
    benches,
    bench_count_reads_small,
    bench_count_reads_adv,
    bench_count_reads_tier3
);
criterion_main!(benches);
