use std::path::Path;
use tempfile::tempdir;
use tempfile::NamedTempFile;

#[test]
fn test_process_fastq_integration() {
    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.fastq");
    let matched_tmp = NamedTempFile::new().expect("create temp file");
    let removed_tmp = NamedTempFile::new().expect("create temp file");

    // Call processing function
    let (total, with_umi, without_umi) = umi_checker::processing::process_fastq(
        &data_path,
        Some(matched_tmp.path()),
        Some(removed_tmp.path()),
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
        Some(matched_tmp.path()),
        Some(removed_tmp.path()),
        2, // allow 1 mismatch
        12,
    )
    .expect("processing failed");

    // From our small FASTQ: read1 and read2 contain the UMI in the sequence (read3 does not)
    assert_eq!(total, 17619);
    assert_eq!(with_umi, 76);
    assert_eq!(without_umi, 17543);
}

// New CLI integration test using a separate process (avoids rayon global build issues).
#[test]
fn test_main_cli_writes_outputs_and_prints_summary() -> Result<(), Box<dyn std::error::Error>> {
    use assert_cmd::assert::OutputAssertExt;
    use assert_cmd::cargo;
    use predicates::prelude::*;
    use std::process::Command;

    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.fastq");
    let tmp = tempdir()?;
    let out_prefix = tmp.path().join("outprefix");

    // Run the compiled binary and assert it prints the expected summary.
    let mut cmd = Command::new(cargo::cargo_bin!(env!("CARGO_PKG_NAME"))); // Use the macro, not cargo::cargo_bin
    cmd.arg("-i")
        .arg(&data_path)
        .arg("-o")
        .arg(&out_prefix)
        .arg("-m")
        .arg("1");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("example.fastq\t3\t2"));

    // Check the output files were written (fastq input -> .fastq outputs)
    let matched = tmp.path().join("outprefix.fq");
    let removed = tmp.path().join("outprefix.removed.fq");
    assert!(matched.exists());
    assert!(removed.exists());

    Ok(())
}
