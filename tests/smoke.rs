use std::path::Path;
use std::process::Command;

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-featurecounts"))
}

fn golden(n: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), n)
}

#[test]
fn count_features_simple() {
    let dir = tempfile::tempdir().expect("tempdir");
    let out_path = dir.path().join("counts.txt");
    let out = bin()
        .args(["-a", &golden("small.gff"), "-o"])
        .arg(&out_path)
        .arg(golden("small.bam"))
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "unexpected failure: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = std::fs::read_to_string(&out_path).expect("counts file not written");
    assert!(
        text.contains("Geneid"),
        "output missing Geneid header: {text}"
    );
    assert!(
        Path::new(&format!("{}.summary", out_path.display())).exists(),
        "summary file not written"
    );
}

#[test]
fn count_features_adversarial() {
    let dir = tempfile::tempdir().expect("tempdir");
    let out_path = dir.path().join("counts.txt");
    let out = bin()
        .args(["-a", &golden("adv.gtf"), "-o"])
        .arg(&out_path)
        .arg(golden("adv.bam"))
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "unexpected failure: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let text = std::fs::read_to_string(&out_path).expect("counts file not written");
    // geneA should have count=3 (r_exon1, r_exon2, r_spliced via CIGAR-N block1).
    // geneB should have count=1 (r_lowmapq0, mapq=0 passes default -Q 0).
    // geneC should have count=0.
    let gene_a_line = text
        .lines()
        .find(|l| l.starts_with("geneA\t"))
        .expect("geneA not in output");
    let gene_a_count: u64 = gene_a_line
        .split('\t')
        .next_back()
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert_eq!(gene_a_count, 3, "geneA count wrong: {gene_a_line}");

    let gene_b_line = text
        .lines()
        .find(|l| l.starts_with("geneB\t"))
        .expect("geneB not in output");
    let gene_b_count: u64 = gene_b_line
        .split('\t')
        .next_back()
        .unwrap()
        .trim()
        .parse()
        .unwrap();
    assert_eq!(gene_b_count, 1, "geneB count wrong: {gene_b_line}");
}
