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
fn zero_byte_mask(x: u64) -> u64 {
    x.wrapping_sub(0x0101010101010101) & !x & 0x8080808080808080
}

/// Compute the Hamming distance between `seq1` and `seq2`.
///
/// This function treats 'N' and 'n' in *either* sequence as a mismatch and is
/// optimized to process 8 bytes at a time using SWAR techniques. The slices
/// must have equal length. Returns the number of positions that differ.
///
/// # Panics
/// Panics in debug builds if the slices are of unequal length.
pub fn hamming_distance_with_n(seq1: &[u8], seq2: &[u8]) -> u32 {
    debug_assert_eq!(seq1.len(), seq2.len());

    let len = seq1.len();
    let mut distance = 0u32;
    let mut i = 0usize;

    const N_UPPER: u64 = 0x4E4E4E4E4E4E4E4E; // 'N' repeated
    const N_LOWER: u64 = 0x6E6E6E6E6E6E6E6E; // 'n' repeated

    // Process 8 bytes at a time
    while i + 8 <= len {
        // Unsafe is acceptable here for performance if you trust bounds,
        // but try_into().unwrap() is safe and usually optimized out by LLVM.
        let chunk1 = u64::from_ne_bytes(seq1[i..i + 8].try_into().unwrap());
        let chunk2 = u64::from_ne_bytes(seq2[i..i + 8].try_into().unwrap());

        let diff = chunk1 ^ chunk2;

        // Create masks where 0x80 is set if the byte is 'N' or 'n'
        let n_mask = zero_byte_mask(chunk1 ^ N_UPPER)
            | zero_byte_mask(chunk1 ^ N_LOWER)
            | zero_byte_mask(chunk2 ^ N_UPPER)
            | zero_byte_mask(chunk2 ^ N_LOWER);

        // A mismatch is a standard difference OR the presence of an N
        distance += count_nonzero_bytes(diff | n_mask);
        i += 8;
    }

    // Handle tail bytes
    while i < len {
        let x = seq1[i];
        let y = seq2[i];
        if x != y || x == b'N' || x == b'n' || y == b'N' || y == b'n' {
            distance += 1;
        }
        i += 1;
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

    // Fallback: If UMI is very short or mismatches are high, the chunking overhead
    // might outweigh the benefit. Just scan everything.
    let num_chunks = (max_mismatches + 1) as usize;
    if umi_len < num_chunks {
        return scan_all_positions(umi, read, max_mismatches);
    }

    // Pigeonhole Search
    // We do NOT allocate a Vec for chunks. We calculate slices on the fly.
    let chunk_size = umi_len / num_chunks;

    // We iterate over every possible start position in the read
    // But we only perform the expensive Hamming check if a "trigger" chunk matches.
    // To implement this efficiently without scanning the whole read multiple times,
    // we iterate positions `i` in the read, and for each position, we check if
    // any of the `k+1` chunks align with the read at that position.

    // Note: The previous implementation scanned the read for *every* chunk independently.
    // That is O(read_len * num_chunks).
    // Below is the simplified scan which might be slightly slower than a tailored
    // multi-pattern search, but it uses zero allocations.

    for i in 0..=(read_len - umi_len) {
        let window = &read[i..i + umi_len];

        // Optimistic check: Does any chunk match exactly?
        let mut potential_match = false;
        for chunk_idx in 0..num_chunks {
            let start = chunk_idx * chunk_size;
            // The last chunk extends to the end of the UMI
            let end = if chunk_idx == num_chunks - 1 {
                umi_len
            } else {
                (chunk_idx + 1) * chunk_size
            };

            if umi[start..end] == window[start..end] {
                potential_match = true;
                break;
            }
        }

        // Only calculate full hamming distance if a pigeonhole chunk matched
        if potential_match && hamming_distance_with_n(umi, window) <= max_mismatches {
            return true;
        }
    }

    false
}

/// Fallback: scan every possible alignment of `umi` in `read` and check
/// Hamming distance at each position.
///
/// This is used when the UMI is short or chunking would be inefficient.
fn scan_all_positions(umi: &[u8], read: &[u8], max_mismatches: u32) -> bool {
    for i in 0..=(read.len() - umi.len()) {
        if hamming_distance_with_n(umi, &read[i..i + umi.len()]) <= max_mismatches {
            return true;
        }
    }
    false
}
