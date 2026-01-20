use anyhow::{Context, Result};
use flate2::write::GzEncoder;
use flate2::Compression;
use rust_htslib::bam;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Generic writer abstraction that can be either a FASTQ writer or a BAM writer.
///
/// This type encapsulates format-specific write logic so higher-level code can
/// work with a single writer type for both FASTQ and BAM outputs.
pub enum GenericWriter {
    Fastq(Box<dyn Write>),
    Bam(bam::Writer),
}

impl GenericWriter {
    /// Write a BAM record to the underlying BAM writer.
    ///
    /// No-op when the `GenericWriter` is not a BAM writer.
    pub fn write_bam(&mut self, rec: &bam::Record) -> Result<()> {
        if let Self::Bam(ref mut w) = self {
            w.write(rec).context("Failed to write BAM record")?;
        }
        Ok(())
    }

    /// Write a FASTQ-formatted record to the underlying writer.
    ///
    /// This writes a single `@<header>\n<seq>\n+\n<qual>` entry; if `qual` is
    /// `None`, a placeholder `+` line is still emitted.
    pub fn write_fastq(&mut self, head: &[u8], seq: &[u8], qual: Option<&[u8]>) -> Result<()> {
        if let Self::Fastq(ref mut w) = self {
            w.write_all(b"@")?;
            w.write_all(head)?;
            w.write_all(b"\n")?;
            w.write_all(seq)?;
            w.write_all(b"\n+\n")?;
            if let Some(q) = qual { w.write_all(q)?; }
            w.write_all(b"\n")?;
        }
        Ok(())
    }
}

/// The common interface for any sequence record.
pub trait BioRecord: Send + Sync {
    fn seq(&self) -> &[u8];
    fn header(&self) -> &[u8];
    fn write_to(self, writer: &mut GenericWriter) -> Result<()>;
}

/// A FASTQ-style in-memory record used for batching and processing.
pub struct FastqRecord {
    /// The header / id field from the FASTQ record (bytes only)
    pub head: Vec<u8>,
    /// Sequence bytes
    pub seq: Vec<u8>,
    /// Optional quality string as bytes
    pub qual: Option<Vec<u8>>,
}

impl BioRecord for FastqRecord {
    fn seq(&self) -> &[u8] { &self.seq }
    fn header(&self) -> &[u8] { &self.head }
    fn write_to(self, writer: &mut GenericWriter) -> Result<()> {
        writer.write_fastq(&self.head, &self.seq, self.qual.as_deref())
    }
}

/// A small wrapper for a BAM record that also stores a copy of the sequence
/// bytes so it can implement `BioRecord` without lifetime issues.
pub struct BamRecord {
    pub rec: bam::Record,
    #[allow(dead_code)] // The seq is read via the trait
    pub seq: Vec<u8>
}

impl BioRecord for BamRecord {
    fn seq(&self) -> &[u8] { &self.seq }
    fn header(&self) -> &[u8] { self.rec.qname() }
    fn write_to(self, writer: &mut GenericWriter) -> Result<()> {
        writer.write_bam(&self.rec)
    }
}

/// Create a writer for FASTQ output. If `path` ends with `.gz`, returns a
/// gzip-wrapped writer.
pub fn create_fastq_writer(path: &Path) -> Result<Box<dyn Write>> {
    let file = File::create(path).with_context(|| format!("Failed to create {}", path.display()))?;
    let writer = BufWriter::new(file);
    if path.extension().map_or(false, |e| e == "gz") {
        Ok(Box::new(GzEncoder::new(writer, Compression::default())))
    } else {
        Ok(Box::new(writer))
    }
}

/// Create a BAM writer from `path` using `header` as a template.
pub fn create_bam_writer(path: &Path, header: &bam::Header) -> Result<bam::Writer> {
    bam::Writer::from_path(path, header, bam::Format::Bam).context("Failed to create BAM writer")
}
