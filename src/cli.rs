use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Section};

use rsomics_quantile_norm::quantile_normalize;

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-quantile-norm", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    pub counts: PathBuf,
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        let n = quantile_normalize(&self.counts, &mut out)?;
        if !self.common.quiet {
            eprintln!("{n} genes quantile-normalized");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Quantile-normalize a gene x sample count matrix (limma normalizeQuantiles).",
    origin: None,
    usage_lines: &["<counts.tsv> [-o out.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[FlagSpec {
            short: Some('o'),
            long: "output",
            aliases: &[],
            value: Some("<path>"),
            type_hint: Some("String"),
            required: false,
            default: Some("-"),
            description: "Output matrix TSV ('-' for stdout).",
            why_default: Some("stdout composes with shell pipelines."),
        }],
    }],
    examples: &[Example {
        description: "Quantile-normalize a count matrix",
        command: "rsomics-quantile-norm counts.tsv -o normalized.tsv",
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
