pub mod io;
pub mod matcher;
pub mod processing;

/// Extract the UMI from a read header.
///
/// The function expects headers like `READ_ID:UMI` or `READ_ID_UMI` and returns
/// the UMI as an uppercase `Vec<u8>` when the extracted UMI length matches
/// `expected_length`. Returns `None` for malformed UTF-8 or if no token is
/// found. **Note:** the function will panic if a UMI is found but its length
/// does not equal `expected_length` to enforce caller invariants.
pub fn extract_umi_from_header(header: &[u8], expected_length: usize) -> Option<Vec<u8>> {
    let header_str = std::str::from_utf8(header).ok()?;

    // Try to find UMI after last ':' or '_' but before any whitespace
    let umi_str = header_str
        .split_whitespace()
        .next()?
        .rsplit([':', '_'])
        .next()?;

    if umi_str.len() != expected_length {
        // Throw an exception if UMI length does not match expected length
        panic!(
            "UMI length does not match expected length: expected {}, found {}",
            expected_length,
            umi_str.len()
        );
    }

    Some(umi_str.as_bytes().to_ascii_uppercase())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extract_umi_from_header() {
        let header = b"READ_12345:ACGTACGTACGT";
        let umi = extract_umi_from_header(header, 12);
        assert_eq!(umi.unwrap(), b"ACGTACGTACGT");
    }
}
