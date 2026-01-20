use std::path::Path;
use tempfile::NamedTempFile;

#[test]
fn test_process_fastq_integration() {
    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.fastq");
    let matched_tmp = NamedTempFile::new().expect("create temp file");
    let removed_tmp = NamedTempFile::new().expect("create temp file");

    // Call processing function
    let (total, with_umi, without_umi) = umi_checker::processing::process_fastq(
        &data_path,
        matched_tmp.path(),
        removed_tmp.path(),
        1, // allow 1 mismatch
        12,
    )
    .expect("processing failed");

    // From our small FASTQ: read1 and read2 contain the UMI in the sequence (read3 does not)
    assert_eq!(total, 3);
    assert_eq!(with_umi, 2);
    assert_eq!(without_umi, 1);
}

#[test]
fn test_process_bam_integration() {
    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.bam");
    let matched_tmp = NamedTempFile::new().expect("create temp file");
    let removed_tmp = NamedTempFile::new().expect("create temp file");

    // Call processing function
    let (total, with_umi, without_umi) = umi_checker::processing::process_bam(
        &data_path,
        matched_tmp.path(),
        removed_tmp.path(),
        2, // allow 1 mismatch
        12,
    )
    .expect("processing failed");

    // From our small FASTQ: read1 and read2 contain the UMI in the sequence (read3 does not)
    assert_eq!(total, 17619);
    assert_eq!(with_umi, 76);
    assert_eq!(without_umi, 17543);
}
