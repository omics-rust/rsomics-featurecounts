use noodles::sam::alignment::record::cigar::op::Kind;

/// Mapped reference blocks from a CIGAR string, starting at `read_start_1based`.
///
/// M/=/X/D advance the reference cursor; N (intron skip) splits blocks; I/S/H/P do not
/// advance. D is included in the block rather than splitting it — a deletion within an
/// exon does not disqualify the read (featureCounts behaviour). Returns 1-based half-open
/// [start, end) spans.
pub(crate) fn cigar_ref_blocks_from<C>(cigar: C, read_start_1based: i32) -> Vec<(i32, i32)>
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
