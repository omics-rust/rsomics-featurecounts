use std::process::{Command, Stdio};

fn ours() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-featurecounts"))
}
fn golden(n: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), n)
}

fn featurecounts_available() -> bool {
    Command::new("featureCounts")
        .arg("-v")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .status()
        .is_ok()
}

/// (gene_id, count) pairs from a featureCounts-style table (count is the last column).
fn counts(table: &str) -> Vec<(String, String)> {
    let mut v: Vec<(String, String)> = table
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("Geneid") && !l.is_empty())
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f.len() >= 2).then(|| (f[0].to_string(), f[f.len() - 1].to_string()))
        })
        .collect();
    v.sort();
    v
}

#[test]
fn runs_without_crash() {
    let out = ours()
        .arg(golden("small.bam"))
        .args(["-a", &golden("small.gff")])
        .output()
        .unwrap();
    assert!(out.status.success());
}

// Per-feature read counts must match `featureCounts` (subread, the named upstream).
#[test]
fn matches_subread_featurecounts() {
    if !featurecounts_available() {
        eprintln!("skipping: featureCounts (subread) not found");
        return;
    }
    let dir = std::env::temp_dir().join("rsomics-featurecounts-compat");
    let _ = std::fs::create_dir_all(&dir);
    let ours_out = dir.join("ours.txt");
    assert!(
        ours()
            .args(["-a", &golden("small.gff")])
            .arg(golden("small.bam"))
            .args(["-o"])
            .arg(&ours_out)
            .status()
            .unwrap()
            .success()
    );
    let sub_out = dir.join("sub.txt");
    assert!(
        Command::new("featureCounts")
            .args(["-a", &golden("small.gff"), "-o"])
            .arg(&sub_out)
            .arg(golden("small.bam"))
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .unwrap()
            .success()
    );
    let ours_c = counts(&std::fs::read_to_string(&ours_out).unwrap());
    let sub_c = counts(&std::fs::read_to_string(&sub_out).unwrap());
    assert!(!ours_c.is_empty());
    assert_eq!(ours_c, sub_c, "per-feature counts must match featureCounts");
}
