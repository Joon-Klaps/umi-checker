use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use rust_htslib::{bam, bam::Read};
use std::path::Path;
use rayon::prelude::*;

use crate::io::{GenericWriter, BioRecord, FastqRecord, BamRecord, create_fastq_writer, create_bam_writer};
use crate::matcher::is_umi_in_read;

const BATCH_SIZE: usize = 10_000;

/// Process a batch of records: perform parallel matching then serial writes.
///
/// The function runs the expensive UMI matching in parallel (with Rayon) and
/// then performs outputs serially to avoid interleaved writes. Returns a tuple
/// `(removed_count, kept_count)` describing how many reads were routed to each
/// output writer.
fn process_batch<R: BioRecord>(
    batch: Vec<R>,
    kept_writer: &mut GenericWriter,
    removed_writer: &mut GenericWriter,
    max_mismatches: u32,
    umi_len: usize,
) -> Result<(usize, usize)> {
    // 1. Parallel compute (Heavy lifting)
    let results: Vec<bool> = batch.par_iter().map(|rec| {
        if let Some(umi) = crate::extract_umi_from_header(rec.header(), umi_len) {
            // This is our SIMD-optimized, zero-allocation matcher
            is_umi_in_read(&umi, rec.seq(), max_mismatches)
        } else {
            false
        }
    }).collect();

    // 2. Serial write (I/O)
    let mut removed = 0;
    let mut kept = 0;
    for (rec, matched) in batch.into_iter().zip(results) {
        if matched {
            removed += 1;
            rec.write_to(removed_writer)?;
        } else {
            kept += 1;
            rec.write_to(kept_writer)?;
        }
    }
    Ok((removed, kept))
}

/// Process an input FASTQ (or gzipped FASTQ) file, separating reads
/// into two outputs: reads containing the UMI (kept) and reads where the UMI
/// was found inside the sequence (removed). Returns `(total, removed, kept)`.
///
/// `max_m` controls allowed mismatches and `umi_len` is the expected UMI length
/// used when extracting the UMI from the read header.
pub fn process_fastq(
    input: &Path, kept_out: &Path, rem_out: &Path, max_m: u32, umi_len: usize,
) -> Result<(usize, usize, usize)> {
    let mut reader = parse_fastx_file(input).context("Failed to parse FASTX file")?;
    let mut kept_w = GenericWriter::Fastq(create_fastq_writer(kept_out)?);
    let mut rem_w = GenericWriter::Fastq(create_fastq_writer(rem_out)?);

    let mut stats = (0, 0, 0); // total, removed, kept
    let mut batch = Vec::with_capacity(BATCH_SIZE);

    while let Some(record) = reader.next() {
        let r = record?;
        stats.0 += 1;
        batch.push(FastqRecord {
            head: r.id().to_vec(),
            seq: r.seq().to_vec(),
            qual: r.qual().map(|q| q.to_vec()),
        });

        if batch.len() >= BATCH_SIZE {
            let (r_inc, k_inc) = process_batch(batch, &mut kept_w, &mut rem_w, max_m, umi_len)?;
            stats.1 += r_inc; stats.2 += k_inc;
            batch = Vec::with_capacity(BATCH_SIZE);
        }
    }
    // Final flush...
    let (r_inc, k_inc) = process_batch(batch, &mut kept_w, &mut rem_w, max_m, umi_len)?;
    stats.1 += r_inc; stats.2 += k_inc;

    Ok(stats)
}

// --- BAM PROCESSOR ---

/// Process an input BAM (or SAM) file, separating reads into `kept_out` and
/// `rem_out` files similarly to `process_fastq`. Uses the BAM header from the
/// input when creating output BAM writers.
pub fn process_bam(
    input: &Path, kept_out: &Path, rem_out: &Path, max_m: u32, umi_len: usize,
) -> Result<(usize, usize, usize)> {
    let mut reader = bam::Reader::from_path(input).context("Failed to open BAM file")?;
    let header = bam::Header::from_template(reader.header());

    // Note: header is used to initialize writers


    let mut kept_w = GenericWriter::Bam(create_bam_writer(kept_out, &header)?);
    let mut rem_w = GenericWriter::Bam(create_bam_writer(rem_out, &header)?);

    let mut stats = (0, 0, 0); // total, removed, kept
    let mut batch = Vec::with_capacity(BATCH_SIZE);

    for result in reader.records() {
        let r = result?;
        stats.0 += 1;
        let seq = r.seq().as_bytes();
        batch.push(BamRecord { rec: r, seq });

        if batch.len() >= BATCH_SIZE {
            let (r_inc, k_inc) = process_batch(batch, &mut kept_w, &mut rem_w, max_m, umi_len)?;
            stats.1 += r_inc; stats.2 += k_inc;
            batch = Vec::with_capacity(BATCH_SIZE);
        }
    }
    // Final flush...
    let (r_inc, k_inc) = process_batch(batch, &mut kept_w, &mut rem_w, max_m, umi_len)?;
    stats.1 += r_inc; stats.2 += k_inc;

    Ok(stats)
}
