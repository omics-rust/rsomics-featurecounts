use std::collections::HashMap;

use coitrees::{COITree, Interval as CoiInterval, IntervalTree as _};

use crate::Exon;

/// Per-chromosome interval tree over exons.
///
/// Metadata stored per interval is a flat exon index; `gene_idx[exon_idx]` maps to the
/// gene-level index. This avoids per-hit String allocation.
pub(crate) struct ExonIndex {
    pub(crate) trees: HashMap<String, COITree<usize, u32>>,
    /// `gene_idx[i]` → gene-level index for exon `i`.
    pub(crate) gene_idx: Vec<usize>,
    /// Gene list in GTF order, indexed by gene-level index.
    pub(crate) gene_ids: Vec<String>,
}

/// Result of a multi-block gene query — no heap allocation for the common single-gene case.
pub(crate) enum GeneHit {
    None,
    /// All overlapping exons belong to one gene. Carries the gene-list index.
    One(usize),
    Multiple,
}

impl ExonIndex {
    pub(crate) fn build(exons: Vec<Exon>) -> Self {
        // Assign per-gene indices in GTF order.
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

    /// Which genes overlap any CIGAR-derived mapped block?
    ///
    /// Allocation-free for the 0-or-1-gene case. Coordinates are 1-based half-open; coitrees
    /// expects end-inclusive, so each block end is decremented by 1 before querying.
    pub(crate) fn gene_hits_for_blocks(&self, chrom: &str, blocks: &[(i32, i32)]) -> GeneHit {
        let Some(tree) = self.trees.get(chrom) else {
            return GeneHit::None;
        };
        let mut first_gene: Option<usize> = None;
        let mut ambiguous = false;
        'outer: for &(block_start, block_end) in blocks {
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
