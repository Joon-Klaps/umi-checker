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
        2, // allow 2 mismatches
        12,
    )
    .expect("processing failed");

    // From our small BAM file
    assert_eq!(total, 17619);
    assert_eq!(with_umi, 76);
    assert_eq!(without_umi, 17543);
}

// CLI integration test using a separate process (avoids rayon global build issues).
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
    let mut cmd = Command::new(cargo::cargo_bin!(env!("CARGO_PKG_NAME")));
    cmd.arg("-i")
        .arg(&data_path)
        .arg("-o")
        .arg(&out_prefix)
        .arg("-m")
        .arg("1");
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("example.fastq\t3\t2"));

    // Check the output files were written (fastq input -> .fq outputs)
    let matched = tmp.path().join("outprefix.fq");
    let removed = tmp.path().join("outprefix.removed.fq");
    assert!(matched.exists());
    assert!(removed.exists());

    Ok(())
}

#[test]
fn test_main_cli_verbose_flag() -> Result<(), Box<dyn std::error::Error>> {
    use assert_cmd::assert::OutputAssertExt;
    use assert_cmd::cargo;
    use predicates::prelude::*;
    use std::process::Command;

    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.fastq");

    let mut cmd = Command::new(cargo::cargo_bin!(env!("CARGO_PKG_NAME")));
    cmd.arg("-i")
        .arg(&data_path)
        .arg("-m")
        .arg("1")
        .arg("--verbose");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Elapsed:"));

    Ok(())
}

#[test]
fn test_main_cli_no_output_files() -> Result<(), Box<dyn std::error::Error>> {
    use assert_cmd::assert::OutputAssertExt;
    use assert_cmd::cargo;
    use predicates::prelude::*;
    use std::process::Command;

    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.fastq");

    // Run without --output flag
    let mut cmd = Command::new(cargo::cargo_bin!(env!("CARGO_PKG_NAME")));
    cmd.arg("-i").arg(&data_path).arg("-m").arg("0");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("example.fastq\t3\t"));

    Ok(())
}

#[test]
fn test_main_cli_invalid_mismatch() -> Result<(), Box<dyn std::error::Error>> {
    use assert_cmd::assert::OutputAssertExt;
    use assert_cmd::cargo;
    use std::process::Command;

    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.fastq");

    let mut cmd = Command::new(cargo::cargo_bin!(env!("CARGO_PKG_NAME")));
    cmd.arg("-i").arg(&data_path).arg("-m").arg("5"); // Too high

    cmd.assert().failure();

    Ok(())
}

#[test]
fn test_main_cli_unsupported_file_type() -> Result<(), Box<dyn std::error::Error>> {
    use assert_cmd::assert::OutputAssertExt;
    use assert_cmd::cargo;
    use predicates::prelude::*;
    use std::process::Command;

    let tmp_file = NamedTempFile::with_suffix(".txt")?;

    let mut cmd = Command::new(cargo::cargo_bin!(env!("CARGO_PKG_NAME")));
    cmd.arg("-i").arg(tmp_file.path()).arg("-m").arg("1");

    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Unsupported file type"));

    Ok(())
}

#[test]
fn test_main_cli_bam_processing() -> Result<(), Box<dyn std::error::Error>> {
    use assert_cmd::assert::OutputAssertExt;
    use assert_cmd::cargo;
    use predicates::prelude::*;
    use std::process::Command;

    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.bam");
    let tmp = tempdir()?;
    let out_prefix = tmp.path().join("outprefix");

    let mut cmd = Command::new(cargo::cargo_bin!(env!("CARGO_PKG_NAME")));
    cmd.arg("-i")
        .arg(&data_path)
        .arg("-o")
        .arg(&out_prefix)
        .arg("-m")
        .arg("2");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("example.bam\t17619\t76"));

    // Check BAM output files were created
    let matched = tmp.path().join("outprefix.bam");
    let removed = tmp.path().join("outprefix.removed.bam");
    assert!(matched.exists());
    assert!(removed.exists());

    Ok(())
}

#[test]
fn test_process_fastq_empty_input_creates_empty_kept() -> Result<(), Box<dyn std::error::Error>> {
    let input = NamedTempFile::new().expect("create temp file");
    let tmp = tempdir().unwrap();
    let matched = tmp.path().join("matched.fq");
    let removed = tmp.path().join("removed.fq");

    let (total, with_umi, without_umi) =
        umi_checker::processing::process_fastq(input.path(), Some(&matched), Some(&removed), 1, 12)
            .expect("processing failed");

    assert_eq!(total, 0);
    assert_eq!(with_umi, 0);
    assert_eq!(without_umi, 0);
    assert!(matched.exists());
    assert_eq!(std::fs::metadata(&matched).unwrap().len(), 0);
    assert!(!removed.exists());

    Ok(())
}

#[test]
fn test_process_bam_empty_input_creates_kept() -> Result<(), Box<dyn std::error::Error>> {
    use rust_htslib::bam::Read;
    let tmp = tempdir().unwrap();
    let input_path = tmp.path().join("empty.sam");
    std::fs::write(&input_path, b"@HD\tVN:1.0\n").expect("write sam header");

    let matched = tmp.path().join("matched.bam");
    let removed = tmp.path().join("removed.bam");

    let (total, with_umi, without_umi) =
        umi_checker::processing::process_bam(&input_path, Some(&matched), Some(&removed), 1, 12)
            .expect("processing failed");

    assert_eq!(total, 0);
    assert_eq!(with_umi, 0);
    assert_eq!(without_umi, 0);
    assert!(matched.exists());
    assert!(removed.exists());

    // Ensure matched BAM has zero records
    let mut reader = rust_htslib::bam::Reader::from_path(&matched)?;
    assert!(reader.records().next().is_none());

    Ok(())
}

#[test]
fn test_main_cli_custom_threads() -> Result<(), Box<dyn std::error::Error>> {
    use assert_cmd::assert::OutputAssertExt;
    use assert_cmd::cargo;
    use predicates::prelude::*;
    use std::process::Command;

    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.fastq");

    let mut cmd = Command::new(cargo::cargo_bin!(env!("CARGO_PKG_NAME")));
    cmd.arg("-i")
        .arg(&data_path)
        .arg("-t")
        .arg("2")
        .arg("-m")
        .arg("1");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("example.fastq\t3\t2"));

    Ok(())
}

#[test]
fn test_main_cli_custom_umi_length() -> Result<(), Box<dyn std::error::Error>> {
    use assert_cmd::assert::OutputAssertExt;
    use assert_cmd::cargo;
    use predicates::prelude::*;
    use std::process::Command;

    let data_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/example.umi10.fastq");

    let mut cmd = Command::new(cargo::cargo_bin!(env!("CARGO_PKG_NAME")));
    cmd.arg("-i")
        .arg(&data_path)
        .arg("-l")
        .arg("10")
        .arg("-m")
        .arg("1");

    cmd.assert()
        .success()
        .stdout(predicate::str::contains("example.umi10.fastq\t3\t"));

    Ok(())
}
