use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::{CountSummary, Gene};

/// Write counts + summary in featureCounts format.
///
/// Counts file: `# Program:…` header, `Geneid\tChr\t…\t<bam_path>` columns, one row per
/// gene with Chr/Start/End/Strand semicolon-joined across exons.
/// Summary file: `Status\t<bam_path>` header + fixed category rows matching .summary order.
pub fn write_output(
    counts: &HashMap<String, u64>,
    genes: &[Gene],
    summary: &CountSummary,
    bam_path: &Path,
    counts_path: &Path,
    summary_path: &Path,
    program_line: &str,
) -> Result<()> {
    let bam_str = bam_path.display().to_string();

    let cf = File::create(counts_path).map_err(RsomicsError::Io)?;
    let mut cw = BufWriter::with_capacity(256 * 1024, cf);
    writeln!(cw, "{program_line}").map_err(RsomicsError::Io)?;
    writeln!(cw, "Geneid\tChr\tStart\tEnd\tStrand\tLength\t{bam_str}").map_err(RsomicsError::Io)?;
    for gene in genes {
        let chrs: Vec<&str> = gene.exons.iter().map(|e| e.chrom.as_str()).collect();
        let starts: Vec<String> = gene.exons.iter().map(|e| e.start.to_string()).collect();
        let ends: Vec<String> = gene.exons.iter().map(|e| e.end.to_string()).collect();
        let strands: Vec<String> = gene.exons.iter().map(|e| e.strand.to_string()).collect();
        let length = gene.length();
        let count = counts.get(&gene.gene_id).copied().unwrap_or(0);
        writeln!(
            cw,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            gene.gene_id,
            chrs.join(";"),
            starts.join(";"),
            ends.join(";"),
            strands.join(";"),
            length,
            count
        )
        .map_err(RsomicsError::Io)?;
    }
    cw.flush().map_err(RsomicsError::Io)?;

    let sf = File::create(summary_path).map_err(RsomicsError::Io)?;
    let mut sw = BufWriter::with_capacity(8 * 1024, sf);
    writeln!(sw, "Status\t{bam_str}").map_err(RsomicsError::Io)?;
    writeln!(sw, "Assigned\t{}", summary.assigned).map_err(RsomicsError::Io)?;
    writeln!(sw, "Unassigned_Unmapped\t{}", summary.unassigned_unmapped)
        .map_err(RsomicsError::Io)?;
    writeln!(sw, "Unassigned_Read_Type\t{}", summary.unassigned_read_type)
        .map_err(RsomicsError::Io)?;
    writeln!(sw, "Unassigned_Singleton\t{}", summary.unassigned_singleton)
        .map_err(RsomicsError::Io)?;
    writeln!(
        sw,
        "Unassigned_MappingQuality\t{}",
        summary.unassigned_mapping_quality
    )
    .map_err(RsomicsError::Io)?;
    writeln!(sw, "Unassigned_Chimera\t{}", summary.unassigned_chimera).map_err(RsomicsError::Io)?;
    writeln!(
        sw,
        "Unassigned_FragmentLength\t{}",
        summary.unassigned_fragment_length
    )
    .map_err(RsomicsError::Io)?;
    writeln!(sw, "Unassigned_Duplicate\t{}", summary.unassigned_duplicate)
        .map_err(RsomicsError::Io)?;
    writeln!(
        sw,
        "Unassigned_MultiMapping\t{}",
        summary.unassigned_multi_mapping
    )
    .map_err(RsomicsError::Io)?;
    writeln!(sw, "Unassigned_Secondary\t{}", summary.unassigned_secondary)
        .map_err(RsomicsError::Io)?;
    writeln!(sw, "Unassigned_NonSplit\t{}", summary.unassigned_non_split)
        .map_err(RsomicsError::Io)?;
    writeln!(
        sw,
        "Unassigned_NoFeatures\t{}",
        summary.unassigned_no_features
    )
    .map_err(RsomicsError::Io)?;
    writeln!(
        sw,
        "Unassigned_Overlapping_Length\t{}",
        summary.unassigned_overlapping_length
    )
    .map_err(RsomicsError::Io)?;
    writeln!(sw, "Unassigned_Ambiguity\t{}", summary.unassigned_ambiguity)
        .map_err(RsomicsError::Io)?;
    sw.flush().map_err(RsomicsError::Io)?;

    Ok(())
}
