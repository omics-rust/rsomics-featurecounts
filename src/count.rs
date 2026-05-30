use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use noodles::bam;
use rsomics_common::{Result, RsomicsError};

use crate::cigar::cigar_ref_blocks_from;
use crate::index::{ExonIndex, GeneHit};
use crate::{CountOpts, CountSummary, Exon};

/// Count reads from one BAM file against the exon index.
///
/// Returns (gene_id → count, summary).
pub fn count_reads(
    bam_path: &Path,
    exons: &[Exon],
    opts: &CountOpts,
) -> Result<(HashMap<String, u64>, CountSummary)> {
    let index = ExonIndex::build(exons.to_vec());

    let file = File::open(bam_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", bam_path.display())))?;
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().map_err(RsomicsError::Io)?;

    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(ToString::to_string)
        .collect();

    let mut counts: HashMap<String, u64> = HashMap::new();
    let mut summary = CountSummary::default();

    for result in reader.records() {
        let record = result.map_err(RsomicsError::Io)?;
        let flags = record.flags();

        if flags.is_unmapped() {
            summary.unassigned_unmapped += 1;
            continue;
        }

        // NH > 1 → multi-mapping read.
        let nh = nh_tag(&record);
        if nh > 1 {
            summary.unassigned_multi_mapping += 1;
            continue;
        }

        let mq = record.mapping_quality().map_or(0, |q| q.get());
        if mq < opts.min_mapq {
            summary.unassigned_mapping_quality += 1;
            continue;
        }

        let Some(tid) = record
            .reference_sequence_id()
            .transpose()
            .map_err(RsomicsError::Io)?
        else {
            summary.unassigned_no_features += 1;
            continue;
        };
        let Some(chrom) = ref_names.get(tid) else {
            summary.unassigned_no_features += 1;
            continue;
        };
        let Some(start) = record
            .alignment_start()
            .transpose()
            .map_err(RsomicsError::Io)?
        else {
            summary.unassigned_no_features += 1;
            continue;
        };
        let read_start_1based = start.get() as i32;

        let blocks = cigar_ref_blocks_from(record.cigar().iter(), read_start_1based);
        if blocks.is_empty() {
            summary.unassigned_no_features += 1;
            continue;
        }

        match index.gene_hits_for_blocks(chrom.as_str(), &blocks) {
            GeneHit::None => summary.unassigned_no_features += 1,
            GeneHit::One(g) => {
                *counts.entry(index.gene_ids[g].clone()).or_insert(0) += 1;
                summary.assigned += 1;
            }
            GeneHit::Multiple => summary.unassigned_ambiguity += 1,
        }
    }

    Ok((counts, summary))
}

/// Extract NH tag; returns 1 if absent (treat as uniquely mapped).
fn nh_tag(record: &bam::Record) -> u32 {
    use noodles::sam::alignment::record::data::field::Tag;
    let data = record.data();
    if let Some(Ok(field)) = data.get(&Tag::ALIGNMENT_HIT_COUNT) {
        use noodles::sam::alignment::record::data::field::Value;
        match field {
            Value::UInt8(v) => return v as u32,
            Value::UInt16(v) => return v as u32,
            Value::UInt32(v) => return v,
            Value::Int8(v) => return v.max(0) as u32,
            Value::Int16(v) => return v.max(0) as u32,
            Value::Int32(v) => return v.max(0) as u32,
            _ => {}
        }
    }
    1
}
