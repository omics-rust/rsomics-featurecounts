use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_featurecounts::{CountOpts, build_genes, count_reads, load_exons, write_output};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(
    name = "rsomics-featurecounts",
    version,
    about,
    long_about = None,
    disable_help_flag = true
)]
pub struct Cli {
    /// Gene annotation file (GTF/GFF).
    #[arg(short = 'a', long = "annotation")]
    annotation: PathBuf,

    /// Input BAM file(s).
    pub inputs: Vec<PathBuf>,

    /// Output counts file (required; a <output>.summary file is also written).
    #[arg(short = 'o', long = "output")]
    output: PathBuf,

    /// Feature type to count (GTF column 3). Long-only to avoid conflict with -t/--threads.
    #[arg(long = "feature-type", default_value = "exon")]
    feature_type: String,

    /// Attribute to group features by (GTF column 9 key).
    #[arg(long = "attribute", default_value = "gene_id")]
    attribute: String,

    /// Strandedness: 0=unstranded, 1=sense, 2=antisense.
    #[arg(short = 's', long = "strandness", default_value_t = 0)]
    strand_specific: u8,

    /// Minimum mapping quality.
    #[arg(short = 'Q', long = "min-mapq", default_value_t = 0)]
    min_mapq: u8,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn execute(self) -> Result<()> {
        let opts = CountOpts {
            feature_type: self.feature_type,
            attribute: self.attribute,
            strand_specific: self.strand_specific,
            min_mapq: self.min_mapq,
        };

        let exons = load_exons(&self.annotation, &opts)?;
        let genes = build_genes(&exons);

        let program_line = format!(
            "# Program:featureCounts v{}; Command:\"rsomics-featurecounts\"",
            env!("CARGO_PKG_VERSION")
        );

        for bam_path in &self.inputs {
            let (counts, summary) = count_reads(bam_path, &exons, &opts)?;

            let summary_path = PathBuf::from(format!("{}.summary", self.output.display()));

            if self.common.json {
                let j = serde_json::json!({
                    "file": bam_path.display().to_string(),
                    "summary": summary,
                    "counts": counts,
                });
                let stdout = std::io::stdout();
                serde_json::to_writer_pretty(stdout.lock(), &j)
                    .map_err(|e| RsomicsError::InvalidInput(format!("JSON: {e}")))?;
                use std::io::Write;
                writeln!(std::io::stdout().lock()).map_err(RsomicsError::Io)?;
            } else {
                write_output(
                    &counts,
                    &genes,
                    &summary,
                    bam_path,
                    &self.output,
                    &summary_path,
                    &program_line,
                )?;
                if !self.common.quiet {
                    eprintln!(
                        "{}: {} assigned, {} no_feature, {} ambiguous, {} multi_mapping / total",
                        bam_path.display(),
                        summary.assigned,
                        summary.unassigned_no_features,
                        summary.unassigned_ambiguity,
                        summary.unassigned_multi_mapping,
                    );
                }
            }
        }

        Ok(())
    }
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }

    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.execute()
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Count reads over genomic features from BAM + GTF/GFF (featureCounts-compatible).",
    origin: Some(Origin {
        upstream: "featureCounts (Subread)",
        upstream_license: "GPL-3.0",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btt656"),
    }),
    usage_lines: &["-a annotation.gtf <input.bam> -o counts.txt"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('a'),
                long: "annotation",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Gene annotation file (GTF/GFF).",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Output counts file. A <output>.summary file is also written.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "feature-type",
                aliases: &[],
                value: Some("<type>"),
                type_hint: Some("String"),
                required: false,
                default: Some("exon"),
                description: "GFF feature type to count.",
                why_default: Some("featureCounts default"),
            },
            FlagSpec {
                short: None,
                long: "attribute",
                aliases: &[],
                value: Some("<key>"),
                type_hint: Some("String"),
                required: false,
                default: Some("gene_id"),
                description: "GFF attribute to group features by.",
                why_default: Some("featureCounts default"),
            },
            FlagSpec {
                short: Some('s'),
                long: "strandness",
                aliases: &[],
                value: Some("<0|1|2>"),
                type_hint: Some("u8"),
                required: false,
                default: Some("0"),
                description: "Strandedness: 0=unstranded, 1=sense, 2=antisense.",
                why_default: Some("featureCounts default"),
            },
            FlagSpec {
                short: Some('Q'),
                long: "min-mapq",
                aliases: &[],
                value: Some("<int>"),
                type_hint: Some("u8"),
                required: false,
                default: Some("0"),
                description: "Minimum mapping quality threshold.",
                why_default: Some("featureCounts default"),
            },
        ],
    }],
    examples: &[Example {
        description: "Count reads per gene",
        command: "rsomics-featurecounts -a genes.gtf input.bam -o counts.txt",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
