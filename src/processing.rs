use anyhow::{Context, Result};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rust_htslib::{bam, bam::Read};
use std::fs;
use std::path::Path;

use crate::io::{
    create_bam_writer, create_fastq_writer, BamRecord, BioRecord, FastqRecord, GenericWriter,
};
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
    if batch.is_empty() {
        return Ok((0, 0));
    }

    // 1. Parallel compute
    let results: Vec<bool> = batch
        .par_iter()
        .map(|rec| {
            if let Some(umi) = crate::extract_umi_from_header(rec.header(), umi_len) {
                is_umi_in_read(&umi, rec.seq(), max_mismatches)
            } else {
                false
            }
        })
        .collect();

    // 2. Serial write
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
    input: &Path,
    kept_out: Option<&Path>,
    rem_out: Option<&Path>,
    max_m: u32,
    umi_len: usize,
) -> Result<(usize, usize, usize)> {
    // Check for 0-byte file BEFORE parsing to avoid parser errors/panics
    if fs::metadata(input)?.len() == 0 {
        // Create empty output if requested, then return
        if let Some(p) = kept_out {
            let _ = create_fastq_writer(p)?;
        }
        return Ok((0, 0, 0));
    }

    let mut reader = match parse_fastx_file(input) {
        Ok(r) => r,
        // If the file is empty the parser returns ParseErrorKind::EmptyFile
        Err(e) if e.kind == needletail::errors::ParseErrorKind::EmptyFile => {
            return Ok((0, 0, 0));
        }
        Err(e) => {
            // Any other parse error is fatal
            return Err(e).context("Failed to parse FASTX file");
        }
    };

    // Initialize writers immediately
    let mut kept_w = match kept_out {
        Some(p) => GenericWriter::Fastq(create_fastq_writer(p)?),
        None => GenericWriter::Sink,
    };
    let mut rem_w = match rem_out {
        Some(p) => GenericWriter::Fastq(create_fastq_writer(p)?),
        None => GenericWriter::Sink,
    };

    let mut stats = (0, 0, 0); // total, removed, kept
    let mut batch = Vec::with_capacity(BATCH_SIZE);

    // Standard loop: no need to peek at the first record manually
    while let Some(record) = reader.next() {
        let r = record?;
        stats.0 += 1;

        // Own the data
        batch.push(FastqRecord {
            head: r.id().to_vec(),
            seq: r.seq().to_vec(),
            qual: r.qual().map(|q| q.to_vec()),
        });

        if batch.len() >= BATCH_SIZE {
            let (r_inc, k_inc) = process_batch(batch, &mut kept_w, &mut rem_w, max_m, umi_len)?;
            stats.1 += r_inc;
            stats.2 += k_inc;
            batch = Vec::with_capacity(BATCH_SIZE);
        }
    }

    // Final flush
    let (r_inc, k_inc) = process_batch(batch, &mut kept_w, &mut rem_w, max_m, umi_len)?;
    stats.1 += r_inc;
    stats.2 += k_inc;

    Ok(stats)
}

// --- BAM PROCESSOR ---

/// Process an input BAM (or SAM) file, separating reads into `kept_out` and
/// `rem_out` files similarly to `process_fastq`. Uses the BAM header from the
/// input when creating output BAM writers.
pub fn process_bam(
    input: &Path,
    kept_out: Option<&Path>,
    rem_out: Option<&Path>,
    max_m: u32,
    umi_len: usize,
) -> Result<(usize, usize, usize)> {
    let mut reader = bam::Reader::from_path(input).context("Failed to open BAM file")?;

    // Read header immediately to setup output writers
    let header = bam::Header::from_template(reader.header());

    // Note: header is used to initialize writers (if provided)
    let mut kept_w = match kept_out {
        Some(p) => GenericWriter::Bam(create_bam_writer(p, &header)?),
        None => GenericWriter::Sink,
    };
    let mut rem_w = match rem_out {
        Some(p) => GenericWriter::Bam(create_bam_writer(p, &header)?),
        None => GenericWriter::Sink,
    };

    let mut stats = (0, 0, 0); // total, removed, kept
    let mut batch = Vec::with_capacity(BATCH_SIZE);

    // Iterate directly. If file is empty (has header but no records),
    // this loop simply won't run, and we flow to the empty final flush.
    for result in reader.records() {
        let r = result?;
        stats.0 += 1;
        let seq = r.seq().as_bytes();
        batch.push(BamRecord { rec: r, seq });

        if batch.len() >= BATCH_SIZE {
            let (r_inc, k_inc) = process_batch(batch, &mut kept_w, &mut rem_w, max_m, umi_len)?;
            stats.1 += r_inc;
            stats.2 += k_inc;
            batch = Vec::with_capacity(BATCH_SIZE);
        }
    }

    // Final flush
    let (r_inc, k_inc) = process_batch(batch, &mut kept_w, &mut rem_w, max_m, umi_len)?;
    stats.1 += r_inc;
    stats.2 += k_inc;

    Ok(stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::FastqRecord;
    use std::io::{Result as IoResult, Write};
    use std::sync::{Arc, Mutex};

    /// Small writer that appends into an Arc<Mutex<Vec<u8>>> so tests can
    /// inspect written bytes after the function under test owns the writer.
    struct SharedWriter(Arc<Mutex<Vec<u8>>>);
    impl Write for SharedWriter {
        fn write(&mut self, buf: &[u8]) -> IoResult<usize> {
            let mut m = self.0.lock().unwrap();
            m.extend_from_slice(buf);
            Ok(buf.len())
        }
        fn flush(&mut self) -> IoResult<()> {
            Ok(())
        }
    }

    #[test]
    fn test_process_batch_fastq_routing() {
        let batch = vec![
            FastqRecord {
                head: b"r1:ACGT".to_vec(),
                seq: b"XXXXACGTYYYY".to_vec(),
                qual: None,
            },
            FastqRecord {
                head: b"r2:TTTT".to_vec(),
                seq: b"AAAAAAAA".to_vec(),
                qual: None,
            },
        ];

        let kept_buf = Arc::new(Mutex::new(Vec::new()));
        let rem_buf = Arc::new(Mutex::new(Vec::new()));
        let mut kept_writer = GenericWriter::Fastq(Box::new(SharedWriter(kept_buf.clone())));
        let mut rem_writer = GenericWriter::Fastq(Box::new(SharedWriter(rem_buf.clone())));

        let (removed, kept) =
            process_batch(batch, &mut kept_writer, &mut rem_writer, 0, 4).unwrap();
        assert_eq!(removed, 1);
        assert_eq!(kept, 1);

        let k = kept_buf.lock().unwrap();
        let r = rem_buf.lock().unwrap();
        assert!(!k.is_empty());
        assert!(!r.is_empty());
        // Check the removed writer contains the expected FASTQ header
        assert!(String::from_utf8_lossy(&r).contains("@r1:ACGT"));
    }
}
