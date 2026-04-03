# lcl833-master

[![PyPI version](https://img.shields.io/pypi/v/lcl833-master.svg)](https://pypi.org/project/lcl833-master/)
[![Python versions](https://img.shields.io/pypi/pyversions/lcl833-master.svg)](https://img.shields.io/pypi/pyversions/lcl833-master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**lcl833-master** computes exact Jones polynomials from braid words using the Kauffman bracket state sum, evaluates them at the fifth root of unity $q = \exp(2\pi i/5)$, and derives calibrated LCL-833 / SATI-CODEX protection metrics ($\delta_{eff}, \alpha_{op}, G_{gap}, \epsilon_{eff}, \omega, T_{min}$) directly in your terminal or Python environment.

## Features

- **Exact Jones Polynomials:** Uses `fractions.Fraction` for rational exponents to maintain exactness during intermediate steps.
- **Parallel Execution:** Multiprocessing support for complex braids (>= 12 crossings) to significantly reduce computation time.
- **LaTeX Export:** Generate publication-ready TeX strings for any computed polynomial.
- **Expanded Presets:** Full support for standard knots from the Rolfsen table (3_1 through 7_1).
- **LCL-833 Metrics:** Calibrated against the trefoil knot anchor ($|J(3_1; q_5)| \approx 1.543$).
- **CLI & API:** Use it as a standalone tool or a Python library.

## Installation

```bash
pip install lcl833-master
```

## CLI Usage

### Basic Commands

```bash
# Run a preset knot (3_1, 4_1, figure8, etc.)
lcl833 knot 3_1

# Run a custom braid word (Artin generators σ_i)
lcl833 braid 1 1 1
lcl833 braid 1 -2 1 -2

# Emit JSON output (includes LaTeX strings)
lcl833 knot 3_1 --json
```

### Advanced Options

- `--table`: Show a comparison table of all preset knots.
- `--genus`: Logical genus (default: 5).
- `--gamma`: Khovanov correction factor (default: 0.05).
- `--max-crossings`: Maximum allowed crossings (default: 18).
- `--no-lcl`: Skip calibrated metrics and only show Jones results.

```bash
lcl833 table
```

## Python API Usage

```python
from lcl833_master import compute_jones_result, compute_lcl833_metrics

# Compute Jones polynomial for the trefoil (σ₁³)
res = compute_jones_result((1, 1, 1))
print(f"Jones Polynomial: {res.normalized_jones_pretty}")
print(f"LaTeX: {res.normalized_jones_latex}")
print(f"Magnitude at q5: {res.magnitude_at_q5:.4f}")

# Compute LCL-833 metrics
metrics = compute_lcl833_metrics(v_j=res.magnitude_at_q5, genus=5)
print(f"Alpha_op: {metrics.alpha_op:.4f}")
print(f"T_min: {metrics.t_min}")
```

## Performance Note

Kauffman bracket computation is $O(2^n)$. For braids with $n \ge 12$, the library automatically utilizes all available CPU cores via `multiprocessing`. Use `--max-crossings` to prevent accidental long-running processes on extremely large diagrams.

## License

MIT License. See `LICENSE` for details.
