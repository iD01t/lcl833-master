# lcl833-master

[![PyPI version](https://img.shields.io/pypi/v/lcl833-master.svg)](https://pypi.org/project/lcl833-master/)
[![Python versions](https://img.shields.io/pypi/pyversions/lcl833-master.svg)](https://img.shields.io/pypi/pyversions/lcl833-master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/iD01t/lcl833-master/actions/workflows/ci.yml/badge.svg)](https://github.com/iD01t/lcl833-master/actions)

**lcl833-master** is a high-performance Python engine for computing exact Jones polynomials and calibrated protection metrics ($LCL-833$, $SATI-CODEX$). It leverages parallel state-sum algorithms and exact algebraic representation to provide a world-class research environment for knot theory and topological invariants.

## 🚀 Key Features

- **Exact Algebraic Core:** Uses `fractions.Fraction` for rational exponents and exact coefficient bookkeeping.
- **High-Performance Parallelism:** Automatic `multiprocessing` support for complex braids ($n \ge 12$).
- **Publication-Ready LaTeX:** Instant generation of TeX strings for Laurent polynomials.
- **Rolfsen Table Integration:** Native presets for standard knots (3_1, 4_1, 5_1, 5_2, 6_1, 6_2, 6_3, 7_1, etc.).
- **Calibrated Protection Models:** Implementation of LCL-833/SATI-CODEX metrics calibrated against the trefoil knot anchor ($|J(3_1; q_5)| \approx 1.543$).
- **Universal CLI & API:** Seamless integration into terminal workflows or larger Python research pipelines.

---

## 🛠 Installation

```bash
pip install lcl833-master
```

Requires **Python 3.10+**.

---

## 💻 CLI Usage

The package provides the `lcl833` command-line entry point.

### Basic Commands
```bash
# Analyze a preset knot (e.g., 4_1 or figure-8)
lcl833 knot 4_1

# Analyze a custom braid word (Artin generators σ_i)
lcl833 braid 1 1 1        # Trefoil
lcl833 braid 1 -2 1 -2    # Figure-8

# Export results as machine-readable JSON (includes LaTeX)
lcl833 knot 3_1 --json
```

### Advanced Research Tools
```bash
# Comparison Table: See metrics for all presets in one view
lcl833 table

# Control the execution engine
lcl833 braid 1 2 1 2 1 2 --no-parallel  # Force single-threaded
lcl833 braid 1 2 1 2 1 2 --max-crossings 24 # Increase performance guard
```

---

## 🐍 Python API Usage

```python
from lcl833_master import compute_jones_result, compute_lcl833_metrics

# 1. Compute Jones Polynomial results
res = compute_jones_result((1, 1, 1))
print(f"Jones Polynomial: {res.normalized_jones_pretty}")
print(f"LaTeX: {res.normalized_jones_latex}")
print(f"Magnitude at q5: {res.magnitude_at_q5:.12f}")

# 2. Derive calibrated LCL-833 metrics
metrics = compute_lcl833_metrics(v_j=res.magnitude_at_q5, genus=5)
print(f"Alpha_op (operational alpha): {metrics.alpha_op:.4f}")
print(f"T_min (minimum time): {metrics.t_min}")
```

---

## 🔬 Mathematical Theory

### The Kauffman Bracket State-Sum
The core algorithm implements the Kauffman bracket state-sum:
$$ \langle L \rangle = \sum_{s} A^{\text{ind}(s)} (-A^2 - A^{-2})^{|s|-1} $$
Where $\text{ind}(s) = \text{zeros}(s) - \text{ones}(s)$. The engine then normalizes the bracket into the Jones polynomial $V_L(q)$ via the substitution $q = A^{-4}$.

### LCL-833 Calibration
The LCL-833 model utilizes the magnitude of the Jones polynomial evaluated at the fifth root of unity ($q = e^{2\pi i / 5}$):
- **Effective Delta:** $\delta_{eff}(K) = |J(K; q_5)| + \gamma R_{Kh}(K)$
- **Operational Alpha:** $\alpha_{op} = \alpha_{ref} \frac{\delta_{eff}(K)}{\delta_{eff}(3_1)}$

---

## ⚡ Performance Note

The state-sum algorithm has a time complexity of $O(2^n)$ where $n$ is the number of crossings. 
- **n < 12:** Single-threaded execution (fast).
- **n >= 12:** Automatic multiprocessing (utilizes all CPU cores).
- **n > 18:** May require significant time. Use `--max-crossings` to override guards if necessary.

---

## 🤝 Contributing

We welcome contributions to the **lcl833-master** engine! 
- **Bugs:** Report issues on the [Bug Tracker](https://github.com/iD01t/lcl833-master/issues).
- **Code:** Pull requests are welcome for algorithm optimizations or new preset knots.

---

## 📜 License

Distributed under the **MIT License**. See `LICENSE` for more information.

---

## 🔗 Links

- **PyPI:** [pypi.org/project/lcl833-master/](https://pypi.org/project/lcl833-master/)
- **Source:** [github.com/iD01t/lcl833-master](https://github.com/iD01t/lcl833-master)
- **Author:** [Guillaume Lessard](https://github.com/iD01t)
