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
}

/// Build output file paths for the matched and removed sets based on the
/// provided `out_prefix` and detected input suffix. The returned pair is
/// `(matched_path, removed_path)`.
fn build_output_paths(
    input: &Path,
    out_prefix: &Path,
) -> anyhow::Result<(std::path::PathBuf, std::path::PathBuf)> {
    let fname = input
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow::anyhow!("Invalid input file name"))?
        .to_lowercase();

    let suffix = if fname.ends_with(".fq.gz") {
        "fq.gz"
    } else if fname.ends_with(".fastq.gz") {
        "fastq.gz"
    } else if fname.ends_with(".fq") {
        "fq"
    } else if fname.ends_with(".fastq") {
        "fastq"
    } else if fname.ends_with(".bam") {
        "bam"
    } else if fname.ends_with(".sam") {
        "sam"
    } else {
        anyhow::bail!("Unsupported input file type: {}", fname)
    };

    let prefix_str = out_prefix.to_string_lossy();
    let base = if prefix_str.ends_with(&format!(".{}", suffix)) {
        prefix_str
            .trim_end_matches(&format!(".{}", suffix))
            .to_string()
    } else {
        prefix_str.to_string()
    };

    let matched = std::path::PathBuf::from(format!("{}.{}", base, suffix));
    let removed = std::path::PathBuf::from(format!("{}.removed.{}", base, suffix));

    Ok((matched, removed))
}

/// CLI entry point: parse args, configure threading, and run the
/// appropriate processor for the input file type. Prints a concise
/// tab-separated summary: total, with_umi, percent_with, without_umi, percent_without.
fn main() -> Result<()> {
    let args = Args::parse();

    // Validate mismatches
    if args.mismatches > 3 {
        anyhow::bail!("Maximum allowed mismatches is 3");
    }

    // Set up thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()?;

    // Determine file type and process
    let file_type = FileType::from_path(&args.input)?;

    // Build output file paths (matched + removed) based on input suffix and provided prefix.
    // If --output is not provided we won't write output files (use None).
    let (clean_output, removed_output) = if let Some(ref out) = args.output {
        let (c, r) = build_output_paths(&args.input, out)?;
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

    println!(
        "{}\t{}\t{}\t{:.2}\t{}\t{:.2}",
        fname, total, with_umi, perc_with, without_umi, perc_without
    );

    if args.verbose {
        println!("Elapsed: {:.3}s", elapsed.as_secs_f64());
    }

    Ok(())
}


#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use std::path::Path;

    #[test]
    fn test_filetype_from_path_variants() {
        assert_eq!(FileType::from_path(Path::new("reads.fq")).unwrap(), FileType::Fastq);
        assert_eq!(FileType::from_path(Path::new("reads.fastq")).unwrap(), FileType::Fastq);
        assert_eq!(FileType::from_path(Path::new("reads.fq.gz")).unwrap(), FileType::FastqGz);
        assert_eq!(FileType::from_path(Path::new("reads.fastq.gz")).unwrap(), FileType::FastqGz);
        assert_eq!(FileType::from_path(Path::new("reads.bam")).unwrap(), FileType::Bam);
        assert_eq!(FileType::from_path(Path::new("reads.sam")).unwrap(), FileType::Sam);

        // Unknown suffix should be an error
        assert!(FileType::from_path(Path::new("reads.unknown")).is_err());
    }

    #[test]
    fn test_build_output_paths_suffix_handling() {
        let input = Path::new("sample.fastq");
        let out_prefix = Path::new("outprefix");
        let (matched, removed) = build_output_paths(input, out_prefix).unwrap();
        assert_eq!(matched, Path::new("outprefix.fastq"));
        assert_eq!(removed, Path::new("outprefix.removed.fastq"));

        // If prefix already contains the suffix it should not duplicate it
        let out_prefix2 = Path::new("outprefix.fastq");
        let (m2, r2) = build_output_paths(input, out_prefix2).unwrap();
        assert_eq!(m2, Path::new("outprefix.fastq"));
        assert_eq!(r2, Path::new("outprefix.removed.fastq"));

        // gz suffix handling
        let input_gz = Path::new("reads.fq.gz");
        let (mg, rg) = build_output_paths(input_gz, Path::new("outprefix.fq.gz")).unwrap();
        assert_eq!(mg, Path::new("outprefix.fq.gz"));
        assert_eq!(rg, Path::new("outprefix.removed.fq.gz"));
    }

    #[test]
    fn test_args_parsing_and_validation() {
        // Minimal good parse
        let args = Args::try_parse_from(&["prog", "-i", "reads.fastq"]).unwrap();
        assert_eq!(args.mismatches, 0);
        assert_eq!(args.umi_length, 12);
        assert_eq!(args.threads, 4);
        assert_eq!(args.output, None);

        // Invalid mismatches (>3) should be a parsing error
        let bad = Args::try_parse_from(&["prog", "-i", "reads.fastq", "-m", "5"]);
        assert!(bad.is_err());
    }
}
