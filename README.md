# umi-checker

A small, fast CLI tool for working with UMIs (Unique Molecular Identifiers).

## Quick Start ✅

- Build locally (debug):

```bash
cargo build
```

- Run tests:

```bash
cargo test --all
```

- Build release binary:

```bash
cargo build --release
```

- Run the CLI (after building):

```bash
./target/debug/umi-checker --help
```

## Development & Tooling 🔧

- Format the code:

```bash
cargo fmt --all
```

- Lint with Clippy (errors on warnings):

```bash
cargo clippy --all-targets --all-features -- -D warnings
```

- Use pre-commit hooks (recommended):

```bash
pip install pre-commit
pre-commit install
pre-commit run --all-files
```

The repository includes a `.pre-commit-config.yaml` that runs basic hygiene hooks and verifies formatting, Clippy, and compilation.

## CI / Releases 🚀

This project includes a GitHub Actions workflow (`.github/workflows/ci.yml`) which:

- Runs format checks, Clippy, and tests on push and PRs
- Builds release artifacts for Linux (x86_64 & aarch64), macOS, and Windows
- Attaches artifacts to a GitHub Release when pushing a tag like `v1.2.3`

Notes about cross-building for `aarch64-unknown-linux-gnu`:

- The workflow uses `cross` (https://github.com/cross-rs/cross). If you want to reproduce locally, install `cross` with `cargo install --locked cross` and run `cross build --target aarch64-unknown-linux-gnu --release`.

## Packaging & Artifacts

Release artifacts are named `umi-checker-<version>-<platform>.zip` where `<version>` is taken from `Cargo.toml`.

## Contributing

Contributions are welcome! Please open an issue or a pull request with a clear description of the change and tests when applicable.

## License

See `LICENSE` (if present) — otherwise contact the maintainers for licensing details.
