# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-04-02

### Added
- **Parallel State-Sum Computation:** Automatic `multiprocessing` support for braids with 12+ crossings.
- **LaTeX Support:** Added `to_latex()` method to `LaurentPolynomial` and `normalized_jones_latex` to `JonesResult`.
- **Expanded Rolfsen Table:** Standard knot notations (3_1, 4_1, 5_1, 5_2, 6_1, 6_2, 6_3, 7_1) added to presets.
- **Improved Performance:** Optimized `LaurentPolynomial` power operations using binary exponentiation.
- Robust handling in `table` command to skip knots that exceed physical metric constraints.

### Changed
- CLI output now includes LaTeX strings in JSON mode.
- `JonesResult` is now a frozen dataclass for better immutability.

## [0.2.0] - 2026-04-02

### Added
- `--table` subcommand to the CLI for comparing all preset knots.
- `--max-crossings` performance guard flag (defaults to 18).
- `CrossingLimitError` for clean failure on large braids.
- `py.typed` marker for Mypy support.
- Comprehensive metadata in `pyproject.toml`.
- Proper library re-exports in `__init__.py`.
- Expanded knot presets: Whitehead link, Borromean rings, Stevedore knot.

### Fixed
- Fixed `--json` flag position (now a sub-argument for `knot` and `braid`).
- Standardized package naming across all files.

## [0.1.1] - 2026-03-25

### Added
- Initial standalone script consolidation.
- Basic CLI structure.
- State-sum implementation of Jones polynomial.
- LCL-833 metrics calibration.

## [0.1.0] - 2026-03-15

### Added
- Initial internal release of the core algorithm.
