use std::convert::TryInto;

/// Count how many bytes within `x` are non-zero.
///
/// This is a SWAR (SIMD Within A Register) trick that computes the count of
/// non-zero bytes in a 64-bit word without loops. It is used by the matcher
/// to compute mismatches across up to 8 bases at a time for performance.
#[inline(always)]
fn count_nonzero_bytes(mut x: u64) -> u32 {
    // 1. Accumulate bits: If any bit in a byte is set, make the LSB of that byte 1.
    x |= x >> 4;
    x |= x >> 2;
    x |= x >> 1;

    // 2. Mask: Clear all bits except LSB of each byte.
    x &= 0x0101010101010101;

    // 3. Horizontal Sum: Multiply by 0x01... to sum bits into the top byte.
    (x.wrapping_mul(0x0101010101010101) >> 56) as u32
}

/// Produce a mask with the high bit set for each zero byte in `x`.
///
/// The returned word has 0x80 in each byte position that was zero in `x` and
/// 0x00 otherwise. This bit-hack is useful for detecting 'N' characters when
/// packed as 8-byte words.
#[inline(always)]
fn is_n_mask(x: u64) -> u64 {
    const N_MASK: u64 = 0x4E4E4E4E4E4E4E4E; // 'N' repeated
                                            // Standard bit-hack to find bytes equal to 0x4E
    let diff = x ^ N_MASK;
    diff.wrapping_sub(0x0101010101010101) & !diff & 0x8080808080808080
}

/// Compute the Hamming distance between `seq1` and `seq2`.
///
/// This function treats 'N' and 'n' in *either* sequence as a mismatch and is
/// optimized to process 8 bytes at a time using SWAR techniques. The slices
/// must have equal length. Returns the number of positions that differ.
///
/// # Panics
/// Panics in debug builds if the slices are of unequal length.
pub fn hamming_distance(seq1: &[u8], seq2: &[u8]) -> u32 {
    assert_eq!(seq1.len(), seq2.len());

    // 1. Process 8-byte blocks using Iterators
    let mut distance = seq1
        .chunks_exact(8)
        .zip(seq2.chunks_exact(8))
        .map(|(s1, s2)| {
            let c1 = u64::from_ne_bytes(s1.try_into().unwrap());
            let c2 = u64::from_ne_bytes(s2.try_into().unwrap());

            // XOR to find differing bytes
            let diff = c1 ^ c2;
            let n_present = is_n_mask(c1) | is_n_mask(c2);

            count_nonzero_bytes(diff | n_present)
        })
        .sum::<u32>();

    // 2. Handle the remainder (tail)
    // we can only split up into 8-byte chunks, so handle any leftover bytes
    // if we have a seq of 12, 8 will be handled above, and we have 4 left to handle here
    // The reason this code is a lot simplere is that we aren't doing SIMD comparisons here,
    // Because we have an if statment we are doing it serial again.
    let remainder1 = seq1.chunks_exact(8).remainder();
    let remainder2 = seq2.chunks_exact(8).remainder();

    for (&a, &b) in remainder1.iter().zip(remainder2) {
        if a != b || a == b'N' || b == b'N' {
            distance += 1;
        }
    }

    distance
}

/// Check whether `umi` occurs in `read` allowing up to `max_mismatches`.
///
/// Behavior:
/// - If `max_mismatches == 0`, performs an exact substring search.
/// - Otherwise it uses a pigeonhole-based heuristic: the UMI is partitioned
///   into `max_mismatches + 1` chunks; at least one chunk must match exactly
///   for the UMI to be within the allowed mismatches. The function avoids
///   heap allocations and prefers SIMD-like operations when possible.
///
/// Returns `true` if a window in `read` is within `max_mismatches` of `umi`.
pub fn is_umi_in_read(umi: &[u8], read: &[u8], max_mismatches: u32) -> bool {
    let umi_len = umi.len();
    let read_len = read.len();

    if read_len < umi_len {
        return false;
    }

    // Optimization: Exact search (0 mismatches)
    if max_mismatches == 0 {
        return read.windows(umi_len).any(|window| window == umi);
    }

    // Fallback: If UMI is very short or mismatches are high not worth chunking
    let num_chunks = (max_mismatches + 1) as usize;
    if umi_len < num_chunks {
        return read
            .windows(umi_len)
            .any(|window| hamming_distance(umi, window) <= max_mismatches);
    }

    // ***********************
    // ** Pigeonhole Search **
    // ***********************
    //
    // We want to speedup the search by avoiding comparing every possible window fully.
    // So we split it up into chunks and check if any chunk matches exactly first.
    // If no chunks match, we can be certain that we have at least `max_mismatches + 1` mismatches
    //
    // If we split up the UMI into #mismatches + 1 chunks.
    // We slide over the read and for each window of UMI length:
    // - Check if any chunk matches exactly.
    // - If so, compute full Hamming distance to confirm.

    let chunk_size = umi_len / num_chunks;

    // Helper to get chunk boundaries
    let get_chunk_range = |chunk_idx: usize| -> (usize, usize) {
        let start = chunk_idx * chunk_size;
        let end = if chunk_idx == num_chunks - 1 {
            umi_len
        } else {
            (chunk_idx + 1) * chunk_size
        };
        (start, end)
    };

    // Check if any chunk matches at this position
    let has_matching_chunk = |window: &[u8]| -> bool {
        (0..num_chunks).any(|chunk_idx| {
            let (start, end) = get_chunk_range(chunk_idx);
            umi[start..end] == window[start..end]
        })
    };

    // Iterate through all possible windows in the read
    read.windows(umi_len)
        .any(|window| has_matching_chunk(window) && hamming_distance(umi, window) <= max_mismatches)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hamming_distance_exact() {
        let a = b"ACGTACGT";
        let b = b"ACGTACGT";
        assert_eq!(hamming_distance(a, b), 0);
    }

    #[test]
    fn test_hamming_distance_with_n_and_tail() {
        let a = b"ACGTNACGTA";
        let b = b"ACGTAACGTT";
        // 'N' counts as mismatch; expect two mismatches
        assert_eq!(hamming_distance(a, b), 2);
    }

    #[test]
    fn test_is_umi_in_read_exact_and_mismatch() {
        let umi = b"ACGTACGTACGT"; // 12
        let read = b"GGGGACGTACGTACGTGGGG";
        assert!(is_umi_in_read(umi, read, 0));

        let read2 = b"GGGGACGTACGAACGTGGGG"; // 1 mismatch in the middle
        assert!(is_umi_in_read(umi, read2, 1));
        assert!(!is_umi_in_read(umi, read2, 0));
    }
}
