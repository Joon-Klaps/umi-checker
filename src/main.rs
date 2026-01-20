use anyhow::Result;
use clap::Parser;
use std::path::{Path, PathBuf};

use umi_checker::processing::{process_bam, process_fastq};

#[derive(Parser, Debug)]
#[command(author, version, about = "UMI presence validator - checks if UMI from header exists in read")]
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

    /// Output file prefix (suffix will be derived from the input).
    /// Example: --output outprefix -> creates outprefix.fastq and outprefix.removed.fastq
    #[arg(short, long)]
    output: PathBuf,

    /// Number of threads for parallel processing
    #[arg(short, long, default_value_t = 4)]
    threads: usize,

    /// Verbose output (show elapsed time)
    #[arg(short, long, default_value_t = false)]
    verbose: bool,
}

#[derive(Debug)]
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
fn build_output_paths(input: &Path, out_prefix: &Path) -> anyhow::Result<(std::path::PathBuf, std::path::PathBuf)> {
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
        prefix_str.trim_end_matches(&format!(".{}", suffix)).to_string()
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

    // Build output file paths (matched + removed) based on input suffix and provided prefix
    let (clean_output, removed_output) = build_output_paths(&args.input, &args.output)?;

    // Start timer
    let start = std::time::Instant::now();

    let (total, with_umi, without_umi) = match file_type {
        FileType::Fastq | FileType::FastqGz => {
            process_fastq(&args.input, &clean_output, &removed_output, args.mismatches, args.umi_length)?
        }
        FileType::Bam | FileType::Sam => {
            process_bam(&args.input, &clean_output, &removed_output, args.mismatches, args.umi_length)?
        }
    };

    let elapsed = start.elapsed();

    // Output concise tab-separated summary
    let perc_with = if total > 0 { (with_umi as f64 / total as f64) * 100.0 } else { 0.0 };
    let perc_without = if total > 0 { (without_umi as f64 / total as f64) * 100.0 } else { 0.0 };

    println!("{}\t{}\t{:.2}\t{}\t{:.2}", total, with_umi, perc_with, without_umi, perc_without);

    if args.verbose {
        println!("Elapsed: {:.3}s", elapsed.as_secs_f64());
    }

    Ok(())
}