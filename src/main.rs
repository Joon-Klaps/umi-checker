use anyhow::Result;
use clap::Parser;
use std::path::{Path, PathBuf};

use umi_checker::processing::{process_bam, process_fastq};

#[derive(Parser, Debug)]
#[command(
    author,
    version,
    about = "UMI presence validator - checks if UMI from header exists in read"
)]
struct Args {
    /// Input file (FASTQ, FASTQ.gz, BAM, or SAM)
    #[arg(short, long)]
    input: PathBuf,

    /// Maximum number of mismatches allowed when finding UMI in read (<=3)
    #[arg(short, long, default_value_t = 0, value_parser = clap::value_parser!(u32).range(0..=3))]
    mismatches: u32,

    /// UMI length in base pairs
    #[arg(short = 'l', long, default_value_t = 12)]
    umi_length: usize,

    /// Optional output file prefix (suffix will be derived from the input).
    /// If not provided, no output files will be written.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Number of threads for parallel processing
    #[arg(short, long, default_value_t = 4)]
    threads: usize,

    /// Verbose output (show elapsed time)
    #[arg(short, long, default_value_t = false)]
    verbose: bool,
}

#[derive(Debug, PartialEq, Eq)]
enum FileType {
    Fastq,
    FastqGz,
    Bam,
    Sam,
}

impl FileType {
    /// Determine the input `FileType` from the filename suffix.
    ///
    /// Supports `.fq`, `.fastq`, `.fq.gz`, `.fastq.gz`, `.bam`, and `.sam`.
    fn from_path(path: &Path) -> anyhow::Result<Self> {
        let fname = path
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or_else(|| anyhow::anyhow!("Invalid file name"))?
            .to_lowercase();

        if fname.ends_with(".fq.gz") || fname.ends_with(".fastq.gz") {
            return Ok(FileType::FastqGz);
        }

        if fname.ends_with(".fq") || fname.ends_with(".fastq") {
            return Ok(FileType::Fastq);
        }

        if fname.ends_with(".bam") {
            return Ok(FileType::Bam);
        }

        if fname.ends_with(".sam") {
            return Ok(FileType::Sam);
        }

        anyhow::bail!("Unsupported file type: {}", fname)
    }

    /// Returns the canonical suffix and acceptable suffix variants for this file type.
    fn suffix_info(&self) -> (&'static str, &'static [&'static str]) {
        match self {
            FileType::Fastq => ("fq", &[".fq", ".fastq"]),
            FileType::FastqGz => ("fq.gz", &[".fq.gz", ".fastq.gz"]),
            FileType::Bam => ("bam", &[".bam"]),
            FileType::Sam => ("sam", &[".sam"]),
        }
    }

    /// Build output file paths for the matched and removed sets based on the
    /// provided `out_prefix` and this file type's suffix. The returned pair is
    /// `(matched_path, removed_path)`.
    fn build_output_paths(&self, out_prefix: &Path) -> (PathBuf, PathBuf) {
        let (suffix, candidates) = self.suffix_info();
        let prefix_str = out_prefix.to_string_lossy();

        // If the prefix ends with any of the acceptable variants, trim that variant.
        let base = candidates
            .iter()
            .find(|s| prefix_str.ends_with(*s))
            .map(|s| prefix_str.trim_end_matches(*s).to_string())
            .unwrap_or_else(|| prefix_str.to_string());

        let matched = PathBuf::from(format!("{}.{}", base, suffix));
        let removed = PathBuf::from(format!("{}.removed.{}", base, suffix));

        (matched, removed)
    }
}

/// Extracted business logic - now testable!
/// Returns formatted summary string instead of printing directly.
fn run(args: Args) -> Result<String> {
    // Validate mismatches
    if args.mismatches > 3 {
        anyhow::bail!("Maximum allowed mismatches is 3");
    }

    // Determine file type and process
    let file_type = FileType::from_path(&args.input)?;

    // Build output file paths (matched + removed) based on input suffix and provided prefix.
    // If --output is not provided we won't write output files (use None).
    let (clean_output, removed_output) = if let Some(ref out) = args.output {
        let (c, r) = file_type.build_output_paths(out);
        (Some(c), Some(r))
    } else {
        (None, None)
    };

    // Start timer
    let start = std::time::Instant::now();

    let (total, with_umi, without_umi) = match file_type {
        FileType::Fastq | FileType::FastqGz => process_fastq(
            &args.input,
            clean_output.as_deref(),
            removed_output.as_deref(),
            args.mismatches,
            args.umi_length,
        )?,
        FileType::Bam | FileType::Sam => process_bam(
            &args.input,
            clean_output.as_deref(),
            removed_output.as_deref(),
            args.mismatches,
            args.umi_length,
        )?,
    };

    let elapsed = start.elapsed();

    // Output concise tab-separated summary
    let perc_with = if total > 0 {
        (with_umi as f64 / total as f64) * 100.0
    } else {
        0.0
    };
    let perc_without = if total > 0 {
        (without_umi as f64 / total as f64) * 100.0
    } else {
        0.0
    };

    // Include input filename as first column for easier aggregation in shell loops
    let fname = args
        .input
        .file_name()
        .and_then(|s| s.to_str())
        .map(|s| s.to_string())
        .unwrap_or_else(|| args.input.to_string_lossy().to_string());

    let mut output = format!(
        "{}\t{}\t{}\t{:.2}\t{}\t{:.2}",
        fname, total, with_umi, perc_with, without_umi, perc_without
    );

    if args.verbose {
        output.push_str(&format!("\nElapsed: {:.3}s", elapsed.as_secs_f64()));
    }

    Ok(output)
}

/// CLI entry point: parse args, configure threading, and delegate to run().
fn main() -> Result<()> {
    let args = Args::parse();

    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    let output = run(args)?;
    println!("{}", output);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_file_type_from_path() {
        assert_eq!(
            FileType::from_path(Path::new("test.fastq")).unwrap(),
            FileType::Fastq
        );
        assert_eq!(
            FileType::from_path(Path::new("test.fq")).unwrap(),
            FileType::Fastq
        );
        assert_eq!(
            FileType::from_path(Path::new("test.fastq.gz")).unwrap(),
            FileType::FastqGz
        );
        assert_eq!(
            FileType::from_path(Path::new("test.fq.gz")).unwrap(),
            FileType::FastqGz
        );
        assert_eq!(
            FileType::from_path(Path::new("test.bam")).unwrap(),
            FileType::Bam
        );
        assert_eq!(
            FileType::from_path(Path::new("test.sam")).unwrap(),
            FileType::Sam
        );
        assert!(FileType::from_path(Path::new("test.txt")).is_err());
    }

    #[test]
    fn test_build_output_paths_fastq() {
        let ft = FileType::Fastq;
        let (matched, removed) = ft.build_output_paths(Path::new("output"));
        assert_eq!(matched, PathBuf::from("output.fq"));
        assert_eq!(removed, PathBuf::from("output.removed.fq"));
    }

    #[test]
    fn test_build_output_paths_with_suffix() {
        let ft = FileType::Fastq;
        let (matched, removed) = ft.build_output_paths(Path::new("output.fastq"));
        assert_eq!(matched, PathBuf::from("output.fq"));
        assert_eq!(removed, PathBuf::from("output.removed.fq"));
    }

    #[test]
    fn test_build_output_paths_bam() {
        let ft = FileType::Bam;
        let (matched, removed) = ft.build_output_paths(Path::new("output"));
        assert_eq!(matched, PathBuf::from("output.bam"));
        assert_eq!(removed, PathBuf::from("output.removed.bam"));
    }

    #[test]
    fn test_run_validates_mismatches() {
        let args = Args {
            input: PathBuf::from("test.fastq"),
            mismatches: 4,
            umi_length: 12,
            output: None,
            threads: 1,
            verbose: false,
        };

        let result = run(args);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Maximum allowed mismatches is 3"));
    }

    #[test]
    fn test_run_invalid_file_type() {
        let args = Args {
            input: PathBuf::from("test.txt"),
            mismatches: 1,
            umi_length: 12,
            output: None,
            threads: 1,
            verbose: false,
        };

        let result = run(args);
        assert!(result.is_err());
        assert!(result
            .unwrap_err()
            .to_string()
            .contains("Unsupported file type"));
    }

    #[test]
    fn test_run_with_real_data() {
        use tempfile::NamedTempFile;

        let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.fastq");

        // Skip if test data doesn't exist
        if !data_path.exists() {
            eprintln!("Skipping test - test data not found");
            return;
        }

        let matched_tmp = NamedTempFile::new().expect("create temp file");
        let out_prefix = matched_tmp.path().parent().unwrap().join("test_output");

        let args = Args {
            input: data_path,
            mismatches: 1,
            umi_length: 12,
            output: Some(out_prefix),
            threads: 1,
            verbose: true,
        };

        let result = run(args);
        assert!(result.is_ok());

        let output = result.unwrap();
        assert!(output.contains("example.fastq"));
        assert!(output.contains("\t3\t")); // total reads
        assert!(output.contains("Elapsed:")); // verbose output
    }
}
