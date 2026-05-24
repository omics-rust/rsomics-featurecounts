//! Compatibility tests: compare rsomics-featurecounts output against featureCounts (subread).
//!
//! Skipped gracefully if featureCounts is not found. Tests both the existing simple
//! golden fixture and the adversarial fixture (CIGAR-N, ambiguous, multimapper, low-mapq).

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden");

fn golden(n: &str) -> PathBuf {
    Path::new(GOLDEN).join(n)
}

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-featurecounts"))
}

fn oracle_bin() -> Option<PathBuf> {
    // Try conda rs-up environment first (the declared oracle location).
    let conda_fc = std::process::Command::new("conda")
        .args([
            "run",
            "-n",
            "rs-up",
            "--no-capture-output",
            "which",
            "featureCounts",
        ])
        .output()
        .ok()
        .filter(|o| o.status.success())
        .and_then(|o| {
            let p = String::from_utf8_lossy(&o.stdout).trim().to_string();
            if p.is_empty() {
                None
            } else {
                Some(PathBuf::from(p))
            }
        });
    if let Some(p) = conda_fc
        && p.exists()
    {
        return Some(p);
    }
    // Fallback: $PATH lookup.
    if Command::new("featureCounts")
        .arg("-v")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok()
    {
        return Some(PathBuf::from("featureCounts"));
    }
    None
}

/// Parse (gene_id → count) from a featureCounts-format counts table, skipping header lines.
fn parse_counts(table: &str) -> Vec<(String, u64)> {
    let mut v: Vec<(String, u64)> = table
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("Geneid") && !l.is_empty())
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            if f.len() < 2 {
                return None;
            }
            let count: u64 = f[f.len() - 1].trim().parse().ok()?;
            Some((f[0].to_string(), count))
        })
        .collect();
    v.sort();
    v
}

/// Parse summary category rows (Status column is skipped, key→value).
fn parse_summary(text: &str) -> Vec<(String, u64)> {
    let mut v: Vec<(String, u64)> = text
        .lines()
        .filter(|l| !l.starts_with("Status") && !l.is_empty())
        .filter_map(|l| {
            let mut cols = l.splitn(2, '\t');
            let key = cols.next()?.trim().to_string();
            let val: u64 = cols.next()?.trim().parse().ok()?;
            Some((key, val))
        })
        .collect();
    v.sort();
    v
}

fn run_compat(gtf: PathBuf, bam: PathBuf) {
    let Some(oracle) = oracle_bin() else {
        eprintln!("SKIP: featureCounts (subread) not found");
        return;
    };

    let dir = tempfile::tempdir().expect("tempdir");

    // --- Oracle ---
    let oracle_counts = dir.path().join("oracle.txt");
    let oracle_status = Command::new(&oracle)
        .args(["-a", gtf.to_str().unwrap(), "-o"])
        .arg(&oracle_counts)
        .arg(&bam)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("failed to run featureCounts");
    assert!(oracle_status.success(), "oracle featureCounts failed");

    // --- Ours ---
    let ours_counts = dir.path().join("ours.txt");
    let ours_status = ours()
        .args(["-a", gtf.to_str().unwrap(), "-o"])
        .arg(&ours_counts)
        .arg(&bam)
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .expect("failed to run rsomics-featurecounts");
    assert!(ours_status.success(), "rsomics-featurecounts failed");

    // --- Compare counts ---
    let oracle_text = std::fs::read_to_string(&oracle_counts).expect("read oracle counts");
    let ours_text = std::fs::read_to_string(&ours_counts).expect("read ours counts");
    let oracle_c = parse_counts(&oracle_text);
    let ours_c = parse_counts(&ours_text);
    assert!(
        !oracle_c.is_empty(),
        "oracle counts table is empty — fixture problem"
    );
    assert_eq!(
        ours_c, oracle_c,
        "per-gene counts mismatch:\n=== ours ===\n{ours_c:#?}\n=== oracle ===\n{oracle_c:#?}"
    );

    // --- Compare summary ---
    let oracle_sum_path = PathBuf::from(format!("{}.summary", oracle_counts.display()));
    let ours_sum_path = PathBuf::from(format!("{}.summary", ours_counts.display()));
    let oracle_sum = std::fs::read_to_string(&oracle_sum_path).expect("read oracle summary");
    let ours_sum = std::fs::read_to_string(&ours_sum_path).expect("read ours summary");
    let oracle_sv = parse_summary(&oracle_sum);
    let ours_sv = parse_summary(&ours_sum);
    assert!(
        !oracle_sv.is_empty(),
        "oracle summary is empty — fixture problem"
    );
    assert_eq!(
        ours_sv, oracle_sv,
        "summary mismatch:\n=== ours ===\n{ours_sum}\n=== oracle ===\n{oracle_sum}"
    );
}

#[test]
fn smoke_runs() {
    let dir = tempfile::tempdir().expect("tempdir");
    let out = ours()
        .args(["-a", golden("small.gff").to_str().unwrap(), "-o"])
        .arg(dir.path().join("counts.txt"))
        .arg(golden("small.bam"))
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "rsomics-featurecounts crashed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn matches_oracle_simple_fixture() {
    run_compat(golden("small.gff"), golden("small.bam"));
}

#[test]
fn matches_oracle_adversarial_fixture() {
    run_compat(golden("adv.gtf"), golden("adv.bam"));
}
