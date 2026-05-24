//! featureCounts-compatible read counting over genomic features.
//!
//! Implements SE default mode: unstranded, one count per read, NH>1 → MultiMapping,
//! ambiguous reads (multi-gene overlap) → Ambiguity. Like featureCounts' default,
//! secondary/supplementary/duplicate alignments are counted (they are only excluded
//! under --primary/--ignoreDup, not yet implemented).
//!
//! Reference: Liao Y, Smyth GK, Shi W. Bioinformatics 2014;30(7):923-930.
//! doi:10.1093/bioinformatics/btt656

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use coitrees::{COITree, Interval as CoiInterval, IntervalTree as _};
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

/// One GTF exon feature with its owning gene (meta-feature).
#[derive(Debug, Clone)]
pub struct Exon {
    pub gene_id: String,
    pub chrom: String,
    /// 1-based inclusive start (GTF coordinate).
    pub start: i32,
    /// 1-based inclusive end (GTF coordinate).
    pub end: i32,
    pub strand: char,
}

/// Per-gene annotation aggregated from its exons, in GTF order.
#[derive(Debug, Clone)]
pub struct Gene {
    pub gene_id: String,
    /// Exons in the order they appear in the GTF.
    pub exons: Vec<Exon>,
}

impl Gene {
    /// Total non-overlapping exonic length (sum of exon lengths, each = end−start+1).
    pub fn length(&self) -> u64 {
        self.exons
            .iter()
            .map(|e| (e.end - e.start + 1) as u64)
            .sum()
    }
}

/// Options controlling read assignment.
#[derive(Debug, Clone)]
pub struct CountOpts {
    /// GFF/GTF feature type column 3 to collect (default "exon").
    pub feature_type: String,
    /// GFF/GTF attribute key for meta-feature grouping (default "gene_id").
    pub attribute: String,
    /// Strandedness: 0=unstranded, 1=sense, 2=antisense.
    pub strand_specific: u8,
    /// Minimum mapping quality (default 0 = all reads pass).
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

/// Assignment category counts — matches featureCounts .summary column order.
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

/// Per-chromosome interval tree index over exons.
/// Metadata index is an exon index; `gene_idx[exon_idx]` maps to the gene-level index.
struct ExonIndex {
    trees: HashMap<String, COITree<usize, u32>>,
    /// Parallel array: `gene_idx[i]` gives the gene-level index for exon `i`.
    gene_idx: Vec<usize>,
    /// Ordered gene list (GTF order), indexed by gene-level index.
    gene_ids: Vec<String>,
}

/// Result of a multi-block gene query — allocation-free for the common single-gene case.
enum GeneHit {
    None,
    /// All overlapping exons belong to the same gene. Carries the gene-list index.
    One(usize),
    Multiple,
}

impl ExonIndex {
    fn build(exons: Vec<Exon>) -> Self {
        // Assign a per-gene index in GTF order.
        let mut gene_order: Vec<&str> = Vec::new();
        let mut gene_map: HashMap<&str, usize> = HashMap::new();
        for e in &exons {
            let id = e.gene_id.as_str();
            if !gene_map.contains_key(id) {
                gene_map.insert(id, gene_order.len());
                gene_order.push(id);
            }
        }
        let gene_idx: Vec<usize> = exons
            .iter()
            .map(|e| *gene_map.get(e.gene_id.as_str()).unwrap())
            .collect();
        let gene_ids: Vec<String> = gene_order.iter().map(|s| (*s).to_string()).collect();

        let mut raw: HashMap<String, Vec<CoiInterval<usize>>> = HashMap::new();
        for (idx, e) in exons.iter().enumerate() {
            // coitrees uses end-inclusive [start, end] with i32.
            raw.entry(e.chrom.clone())
                .or_default()
                .push(CoiInterval::new(e.start, e.end, idx));
        }
        let trees = raw
            .into_iter()
            .map(|(chrom, ivs)| (chrom, COITree::new(&ivs)))
            .collect();
        ExonIndex {
            trees,
            gene_idx,
            gene_ids,
        }
    }

    /// Returns which genes overlap any CIGAR-derived mapped block of the read.
    ///
    /// For the common case (0 or 1 gene), no heap allocation occurs.
    /// Uses 1-based inclusive exon coordinates (coitrees convention).
    fn gene_hits_for_blocks(&self, chrom: &str, blocks: &[(i32, i32)]) -> GeneHit {
        let Some(tree) = self.trees.get(chrom) else {
            return GeneHit::None;
        };
        let mut first_gene: Option<usize> = None;
        let mut ambiguous = false;
        'outer: for &(block_start, block_end) in blocks {
            // Convert half-open [block_start, block_end) to coitrees end-inclusive.
            tree.query(block_start, block_end - 1, |node| {
                if ambiguous {
                    return;
                }
                // coitrees metadata may be `usize` or `&usize` depending on toolchain.
                let exon_idx: usize = *std::borrow::Borrow::<usize>::borrow(&node.metadata);
                let g = self.gene_idx[exon_idx];
                match first_gene {
                    None => first_gene = Some(g),
                    Some(prev) if prev == g => {}
                    Some(_) => ambiguous = true,
                }
            });
            if ambiguous {
                break 'outer;
            }
        }
        if ambiguous {
            GeneHit::Multiple
        } else {
            match first_gene {
                None => GeneHit::None,
                Some(g) => GeneHit::One(g),
            }
        }
    }
}

/// Parse GTF/GFF, collecting features of the requested type, returning flat exon list.
pub fn load_exons(gff_path: &Path, opts: &CountOpts) -> Result<Vec<Exon>> {
    let file = File::open(gff_path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", gff_path.display())))?;
    let reader = BufReader::new(file);
    let mut exons = Vec::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 9 || fields[2] != opts.feature_type {
            continue;
        }
        let chrom = fields[0].to_string();
        let start: i32 = fields[3]
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("line {}: bad start", lineno + 1)))?;
        let end: i32 = fields[4]
            .parse()
            .map_err(|_| RsomicsError::InvalidInput(format!("line {}: bad end", lineno + 1)))?;
        let strand = fields[6].chars().next().unwrap_or('.');
        let gene_id = extract_attr(fields[8], &opts.attribute).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "line {}: attribute '{}' not found",
                lineno + 1,
                opts.attribute
            ))
        })?;
        exons.push(Exon {
            gene_id,
            chrom,
            start,
            end,
            strand,
        });
    }

    Ok(exons)
}

/// Build ordered gene list from flat exon list (GTF order preserved).
pub fn build_genes(exons: &[Exon]) -> Vec<Gene> {
    let mut order: Vec<String> = Vec::new();
    let mut map: HashMap<String, Gene> = HashMap::new();
    for e in exons {
        if !map.contains_key(&e.gene_id) {
            order.push(e.gene_id.clone());
            map.insert(
                e.gene_id.clone(),
                Gene {
                    gene_id: e.gene_id.clone(),
                    exons: Vec::new(),
                },
            );
        }
        map.get_mut(&e.gene_id).unwrap().exons.push(e.clone());
    }
    order
        .into_iter()
        .map(|id| map.remove(&id).unwrap())
        .collect()
}

fn extract_attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        let Some(rest) = part.strip_prefix(key) else {
            continue;
        };
        if rest.starts_with('=') || rest.starts_with(' ') {
            let val = rest[1..].trim().trim_matches('"');
            if !val.is_empty() {
                return Some(val.to_string());
            }
        }
    }
    None
}

/// Compute mapped reference blocks from a CIGAR string, starting at `read_start_1based`.
///
/// M/=/X/D all advance the reference cursor; N (intron skip) splits blocks; I/S/H/P
/// do not advance the cursor. Returns 1-based half-open [block_start, block_end) spans.
/// D is included in the block rather than splitting it, matching featureCounts behaviour
/// (a deletion within an exon does not disqualify the read from assignment).
fn cigar_ref_blocks_from<C>(cigar: C, read_start_1based: i32) -> Vec<(i32, i32)>
where
    C: IntoIterator<Item = std::io::Result<noodles::sam::alignment::record::cigar::Op>>,
{
    let mut blocks: Vec<(i32, i32)> = Vec::new();
    let mut ref_pos = read_start_1based;
    let mut in_block = false;
    let mut block_start = ref_pos;

    for op_result in cigar {
        let Ok(op) = op_result else { continue };
        let len = op.len() as i32;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                if !in_block {
                    block_start = ref_pos;
                    in_block = true;
                }
                ref_pos += len;
            }
            Kind::Skip => {
                if in_block {
                    blocks.push((block_start, ref_pos));
                    in_block = false;
                }
                ref_pos += len;
            }
            Kind::Insertion | Kind::SoftClip | Kind::HardClip | Kind::Pad => {}
        }
    }
    if in_block {
        blocks.push((block_start, ref_pos));
    }
    blocks
}

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

        // NH tag > 1 → multi-mapping read.
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

/// Extract NH (number of hits) tag from a BAM record. Returns 1 if absent.
fn nh_tag(record: &bam::Record) -> u32 {
    use noodles::sam::alignment::record::data::field::Tag;
    let data = record.data();
    // NH tag: number of reported alignments for this read.
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

/// Write counts + summary in featureCounts format.
///
/// The counts file has:
///   - `# Program:featureCounts ...` header line
///   - `Geneid\tChr\tStart\tEnd\tStrand\tLength\t<bam_path>` column header
///   - One row per gene (meta-feature), Chr/Start/End/Strand semicolon-joined per exon
///
/// The summary file has: `Status\t<bam_path>` header + fixed category rows.
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

    // Counts file.
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

    // Summary file.
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
