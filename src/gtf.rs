use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

use crate::{CountOpts, Exon, Gene};

/// Parse GTF/GFF, collecting features of the requested type.
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

/// GFF/GTF column-9 attribute parser; handles both `key=val` (GFF3) and `key "val"` (GTF).
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
