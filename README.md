# umi-checker

A small, fast CLI tool for checking the presence of UMIs (Unique Molecular Identifiers) in fastq files or BAM files.

## Quick Start ✅

You can install the latest release with a single command:

```bash
curl -fsSL https://raw.githubusercontent.com/Joon-Klaps/umi-checker/main/install.sh | bash
```

```bash
umi-checker --help

UMI presence validator - checks if UMI from header exists in read

Usage: umi-checker [OPTIONS] --input <INPUT> --output <OUTPUT>

Options:
  -i, --input <INPUT>            Input file (FASTQ, FASTQ.gz, BAM, or SAM)
  -m, --mismatches <MISMATCHES>  Maximum number of mismatches allowed when finding UMI in read (<=3) [default: 0]
  -l, --umi-length <UMI_LENGTH>  UMI length in base pairs [default: 12]
  -o, --output <OUTPUT>          Output file prefix (suffix will be derived from the input). Example: --output outprefix -> creates outprefix.fastq and outprefix.removed.fastq
  -t, --threads <THREADS>        Number of threads for parallel processing [default: 4]
  -v, --verbose                  Verbose output (show elapsed time)
  -h, --help                     Print help
  -V, --version                  Print version
```

The output printed to sdout will contain the following tab-separated columns:

- read: Input read file name
- total reads: Total number of reads processed
- reads with umi: Number of reads where UMI was found in the sequence
- % with umi: Percentage of reads with UMI
- reads without umi: Number of reads where UMI was not found in the sequence
- %perc without umi: Percentage of reads without UMI

Example usage:

```bash
echo -e "read \ttotal reads\treads with umi\t% with umi\treads without umi\t%perc without umi" > abundance.tsv
for read in `*.fastq.gz`; do
    umi-checker -i $read -o ${read%.fastq.gz}.umi_checked >> abundance.tsv
done
```

If you only want a tab-separated summary on stdout (for aggregating across many files) and don't want output files created, omit `--output`. The tool will print a single line with the input filename as the first column and will not write any output files.

## Contributing 🧑‍💻

Contributions are welcome! Please open an issue or a pull request with a clear description of the change and tests when applicable.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
