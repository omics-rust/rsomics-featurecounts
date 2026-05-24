# rsomics-featurecounts

Count reads overlapping genomic features from BAM + GTF/GFF annotation —
a Rust port of featureCounts (Subread).

## Usage

```bash
rsomics-featurecounts -a genes.gtf input.bam -o counts.txt
```

Output: `counts.txt` (tab-separated counts table) and `counts.txt.summary`
(assignment category breakdown), matching featureCounts format.

## Covered modes (byte-exact vs featureCounts v2.1.1)

| Mode | Status |
|---|---|
| SE default (`-t exon -g gene_id`, unstranded, NH>1 skipped) | byte-exact |
| Minimum MAPQ (`-Q N`) | byte-exact |
| Strandedness (`-s 0/1/2`) flag wired | NOT YET (parsed, not enforced in overlap logic) |
| Paired-end fragment counting (`-p --countReadPairs`) | NOT YET |
| Multi-overlap (`-O`) | NOT YET |
| Fractional counting | NOT YET |
| Primary-only / dedup (`--primary`, `--ignoreDup`) | NOT YET (default counts secondary/supplementary/duplicate, matching featureCounts) |

## Algorithm

1. Parse GTF/GFF, extract features matching `--feature-type` (default: `exon`).
2. Build a per-chromosome coitrees interval tree (O(N log N) build, O(log N + k) query).
3. For each mapped, NH=1 BAM record (secondary/supplementary/duplicate included, per featureCounts default):
   - Compute CIGAR-derived mapped blocks (M/=/X/D extend, N splits blocks, I/S/H/P ignored).
   - Query the interval tree for each block; collect distinct gene hits.
   - NH > 1 → `Unassigned_MultiMapping`; 0 genes → `Unassigned_NoFeatures`;
     1 gene → `Assigned`; ≥ 2 genes → `Unassigned_Ambiguity`.
4. Write featureCounts-format counts table and `.summary` file.

The CIGAR-N block splitting ensures spliced RNA reads are evaluated by their
mapped segments, not their full genomic span.

## Origin

This crate is an independent Rust reimplementation of `featureCounts` based on:

- The published method: Liao Y, Smyth GK, Shi W. "featureCounts: an efficient
  general purpose program for assigning sequence reads to genomic features."
  Bioinformatics 2014;30(7):923-930. doi:10.1093/bioinformatics/btt656
- The public GTF/GFF and SAM/BAM format specifications
- Black-box behavior testing against the upstream binary (featureCounts v2.1.1)

No source code from the GPL-3.0 upstream (Subread/featureCounts) was used as
reference during implementation. Test fixtures are independently generated
synthetic data.

License: MIT OR Apache-2.0.
Upstream credit: [featureCounts (Subread)](http://subread.sourceforge.net/) (GPL-3.0).
