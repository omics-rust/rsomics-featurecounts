//! featureCounts-compatible read counting over genomic features.
//!
//! SE default: unstranded, one count per read, NH>1 → MultiMapping,
//! multi-gene overlap → Ambiguity.
//!
//! Reference: Liao Y, Smyth GK, Shi W. Bioinformatics 2014;30(7):923-930.
//! doi:10.1093/bioinformatics/btt656

pub mod cigar;
pub mod count;
pub mod gtf;
pub mod index;
pub mod output;

use serde::Serialize;

pub use count::count_reads;
pub use gtf::{build_genes, load_exons};
pub use output::write_output;

#[derive(Debug, Clone)]
pub struct Exon {
    pub gene_id: String,
    pub chrom: String,
    /// 1-based inclusive (GTF coordinate).
    pub start: i32,
    /// 1-based inclusive (GTF coordinate).
    pub end: i32,
    pub strand: char,
}

#[derive(Debug, Clone)]
pub struct Gene {
    pub gene_id: String,
    pub exons: Vec<Exon>,
}

impl Gene {
    /// Sum of exon lengths; does not deduplicate overlapping exons.
    pub fn length(&self) -> u64 {
        self.exons
            .iter()
            .map(|e| (e.end - e.start + 1) as u64)
            .sum()
    }
}

#[derive(Debug, Clone)]
pub struct CountOpts {
    /// GFF/GTF feature type (column 3).
    pub feature_type: String,
    /// GTF attribute key for meta-feature grouping.
    pub attribute: String,
    /// 0 = unstranded, 1 = sense, 2 = antisense.
    pub strand_specific: u8,
    pub min_mapq: u8,
}

impl Default for CountOpts {
    fn default() -> Self {
        Self {
            feature_type: "exon".to_string(),
            attribute: "gene_id".to_string(),
            strand_specific: 0,
            min_mapq: 0,
        }
    }
}

/// Matches featureCounts .summary column order.
#[derive(Debug, Default, Clone, Serialize)]
pub struct CountSummary {
    pub assigned: u64,
    pub unassigned_unmapped: u64,
    pub unassigned_read_type: u64,
    pub unassigned_singleton: u64,
    pub unassigned_mapping_quality: u64,
    pub unassigned_chimera: u64,
    pub unassigned_fragment_length: u64,
    pub unassigned_duplicate: u64,
    pub unassigned_multi_mapping: u64,
    pub unassigned_secondary: u64,
    pub unassigned_non_split: u64,
    pub unassigned_no_features: u64,
    pub unassigned_overlapping_length: u64,
    pub unassigned_ambiguity: u64,
}
