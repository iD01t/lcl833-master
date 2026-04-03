#!/usr/bin/env python3
"""
LCL-833 Master Core
===================
Core library for:
- exact Jones polynomial computation from a braid word via the Kauffman bracket state sum,
- high-performance parallel state-sum computation for complex braids,
- publication-ready LaTeX export for Laurent polynomials,
- evaluation at the fifth root of unity q = exp(2πi/5),
- calibrated protection metrics (δ_eff, α_op, G_gap, ε_eff, ω, T_min),
- validation and self-tests.
"""

from __future__ import annotations

import argparse
import cmath
import json
import logging
import math
import multiprocessing
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass
from fractions import Fraction
from typing import DefaultDict, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

LOGGER = logging.getLogger("lcl833")

TREFOIL_MAGNITUDE = 1.543361918426817
TREFOIL_ALPHA_OP = 0.8783
DEFAULT_GAMMA = 0.05
MACHINE_EPSILON_DOUBLE = 2.0 ** -52
FIFTH_ROOT_OF_UNITY = cmath.exp(2j * math.pi / 5)
MAX_CROSSINGS_DEFAULT = 18


class LCL833Error(Exception):
    """Base exception for this module."""


class InputValidationError(LCL833Error):
    """Raised for malformed braid words or invalid model parameters."""


class CrossingLimitError(LCL833Error):
    """Raised when braid word exceeds the specified crossing limit."""


class LaurentPolynomial:
    """Sparse Laurent polynomial with rational exponents and exact coefficient bookkeeping.

    Exponents are represented with ``fractions.Fraction`` so the Jones polynomial can be kept
    exact under the A -> q substitution, even when an intermediate expression contains quarter
    powers before collapsing to integer powers in the standard knot normalization.
    """

    def __init__(self, terms: Mapping[Fraction | int, complex | int | float] | None = None) -> None:
        self.terms: DefaultDict[Fraction, complex] = defaultdict(complex)
        if terms is not None:
            for power, coefficient in terms.items():
                frac_power = power if isinstance(power, Fraction) else Fraction(power)
                complex_coefficient = complex(coefficient)
                if complex_coefficient != 0:
                    self.terms[frac_power] += complex_coefficient
        self._clean()

    def _clean(self) -> None:
        for power in list(self.terms.keys()):
            if abs(self.terms[power]) < 1e-14:
                del self.terms[power]

    def copy(self) -> "LaurentPolynomial":
        return LaurentPolynomial(dict(self.terms))

    def __bool__(self) -> bool:
        return bool(self.terms)

    def __add__(self, other: "LaurentPolynomial | int | float | complex") -> "LaurentPolynomial":
        if not isinstance(other, LaurentPolynomial):
            other = LaurentPolynomial({0: other})
        out = self.copy()
        for power, coefficient in other.terms.items():
            out.terms[power] += coefficient
        out._clean()
        return out

    __radd__ = __add__

    def __neg__(self) -> "LaurentPolynomial":
        return LaurentPolynomial({p: -c for p, c in self.terms.items()})

    def __sub__(self, other: "LaurentPolynomial | int | float | complex") -> "LaurentPolynomial":
        if not isinstance(other, LaurentPolynomial):
            other = LaurentPolynomial({0: other})
        return self + (-other)

    def __rsub__(self, other: "LaurentPolynomial | int | float | complex") -> "LaurentPolynomial":
        if not isinstance(other, LaurentPolynomial):
            other = LaurentPolynomial({0: other})
        return other + (-self)

    def __mul__(self, other: "LaurentPolynomial | int | float | complex") -> "LaurentPolynomial":
        if not isinstance(other, LaurentPolynomial):
            other = LaurentPolynomial({0: other})
        out: DefaultDict[Fraction, complex] = defaultdict(complex)
        for p1, c1 in self.terms.items():
            for p2, c2 in other.terms.items():
                out[p1 + p2] += c1 * c2
        return LaurentPolynomial(out)

    __rmul__ = __mul__

    def __truediv__(self, other: int | float | complex) -> "LaurentPolynomial":
        if isinstance(other, (int, float, complex)):
            if other == 0:
                raise ZeroDivisionError("Cannot divide LaurentPolynomial by zero.")
            return LaurentPolynomial({p: c / other for p, c in self.terms.items()})
        return NotImplemented

    def __pow__(self, power: int) -> "LaurentPolynomial":
        if power < 0:
            raise InputValidationError("Polynomial exponent must be a non-negative integer.")
        result = LaurentPolynomial({0: 1})
        base = self
        while power > 0:
            if power % 2 == 1:
                result = result * base
            base = base * base
            power //= 2
        return result

    def monomial_mul(self, power: Fraction | int, coefficient: complex | int | float = 1) -> "LaurentPolynomial":
        frac_power = power if isinstance(power, Fraction) else Fraction(power)
        coeff = complex(coefficient)
        return LaurentPolynomial({p + frac_power: c * coeff for p, c in self.terms.items()})

    def evaluate(self, value: complex) -> complex:
        total = 0j
        for power, coefficient in self.terms.items():
            total += coefficient * (value ** float(power))
        return total

    def as_serializable_dict(self) -> Dict[str, float | int]:
        """Serialize to JSON-friendly mapping with normalized exponent labels."""
        payload: Dict[str, float | int] = {}
        for power in sorted(self.terms.keys(), reverse=True):
            coeff = self.terms[power]
            if abs(coeff.imag) > 1e-14:
                raise LCL833Error("Polynomial has non-real coefficients; JSON export is disabled for this object.")
            if power.denominator == 1:
                power_key = str(power.numerator)
            else:
                power_key = f"{power.numerator}/{power.denominator}"
            real_coeff = coeff.real
            if abs(real_coeff - round(real_coeff)) < 1e-14:
                payload[power_key] = int(round(real_coeff))
            else:
                payload[power_key] = real_coeff
        return payload

    def to_pretty_string(self, variable: str = "q") -> str:
        if not self.terms:
            return "0"
        parts: List[str] = []
        for power in sorted(self.terms.keys(), reverse=True):
            coeff = self.terms[power]
            if abs(coeff.imag) > 1e-14:
                raise LCL833Error("Pretty formatting expects real coefficients.")
            real_coeff = coeff.real
            coeff_int = int(round(real_coeff)) if abs(real_coeff - round(real_coeff)) < 1e-14 else real_coeff

            if power == 0:
                term = f"{coeff_int}"
            else:
                if coeff_int == 1:
                    coeff_part = ""
                elif coeff_int == -1:
                    coeff_part = "-"
                else:
                    coeff_part = f"{coeff_int}"

                if power == 1:
                    power_part = variable
                else:
                    if power.denominator == 1:
                        power_label = str(power.numerator)
                    else:
                        power_label = f"({power.numerator}/{power.denominator})"
                    power_part = f"{variable}^{power_label}"
                term = f"{coeff_part}{power_part}"
            parts.append(term)
        return " + ".join(parts).replace("+ -", "- ")

    def to_latex(self, variable: str = "q") -> str:
        """Export the polynomial to TeX format for publication."""
        if not self.terms:
            return "0"
        parts: List[str] = []
        for power in sorted(self.terms.keys(), reverse=True):
            coeff = self.terms[power]
            if abs(coeff.imag) > 1e-14:
                raise LCL833Error("LaTeX formatting expects real coefficients.")
            real_coeff = coeff.real
            coeff_val = int(round(real_coeff)) if abs(real_coeff - round(real_coeff)) < 1e-14 else round(real_coeff, 4)

            if power == 0:
                term = f"{coeff_val}"
            else:
                if coeff_val == 1:
                    coeff_part = ""
                elif coeff_val == -1:
                    coeff_part = "-"
                else:
                    coeff_part = f"{coeff_val}"

                if power == 1:
                    power_part = variable
                else:
                    if power.denominator == 1:
                        power_label = str(power.numerator)
                    else:
                        power_label = f"{power.numerator}/{power.denominator}"
                    power_part = f"{variable}^{{{power_label}}}"
                term = f"{coeff_part}{power_part}"
            
            if term.startswith("-"):
                parts.append(term)
            else:
                parts.append("+" + term if parts else term)

        return " ".join(parts).replace("+", "+ ").replace("-", "- ").strip()


class UnionFind:
    """Disjoint-set structure for loop counting in braid closures."""

    def __init__(self, size: int) -> None:
        if size < 0:
            raise InputValidationError("UnionFind size must be non-negative.")
        self.parent = list(range(size))

    def find(self, item: int) -> int:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        # Path compression
        curr = item
        while self.parent[curr] != root:
            next_node = self.parent[curr]
            self.parent[curr] = root
            curr = next_node
        return root

    def union(self, left: int, right: int) -> None:
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left != root_right:
            self.parent[root_left] = root_right


@dataclass(frozen=True)
class JonesResult:
    braid_word: Tuple[int, ...]
    strand_count: int
    writhe: int
    kaufman_bracket_A: Dict[str, float | int]
    normalized_jones_q: Dict[str, float | int]
    normalized_jones_pretty: str
    normalized_jones_latex: str
    value_at_q5_real: float
    value_at_q5_imag: float
    magnitude_at_q5: float


@dataclass(frozen=True)
class LCL833Metrics:
    genus: int
    gamma: float
    v_j: float
    r_kh: float
    delta_eff: float
    alpha_ref: float
    alpha_op: float
    g_gap: float
    eps_eff: float
    omega: float
    t_min: int


PRESET_KNOTS: Dict[str, Tuple[int, ...]] = {
    "unknot": tuple(),
    "trefoil": (1, 1, 1),
    "3_1": (1, 1, 1),
    "figure8": (1, -2, 1, -2),
    "figure-8": (1, -2, 1, -2),
    "4_1": (1, -2, 1, -2),
    "cinquefoil": (1, 1, 1, 1, 1),
    "5_1": (1, 1, 1, 1, 1),
    "5_2": (1, 1, 1, 2, -1, 2),
    "6_1": (1, 1, 1, -2, 1, -2),
    "6_2": (1, 1, 2, -1, 2, -1, 2),
    "6_3": (1, 2, -1, 2, -1, -2, 1, -2),
    "7_1": (1, 1, 1, 1, 1, 1, 1),
    "hopf-link": (1, 1),
    "whitehead-link": (1, 1, -2, 1, 1, -2),
    "borromean-rings": (1, -2, 1, -2, 1, -2),
    "stevedore": (1, 1, 1, 1, -2, 1, -2),
}

# Reduced Khovanov total-rank ratios where explicitly available to the caller.
# By default, the module uses R_Kh = 0 unless the user supplies a value.
PRESET_RKH: Dict[str, float] = {
    "unknot": 0.0,
    "trefoil": 0.0,
    "3_1": 0.0,
    "figure8": 0.0,
    "figure-8": 0.0,
    "4_1": 0.0,
    "cinquefoil": 0.0,
    "5_1": 0.0,
    "hopf-link": 0.0,
}


def validate_braid_word(braid_word: Sequence[int], max_crossings: int = 0) -> Tuple[int, ...]:
    """Validate a braid word.

    Conditions:
    - entries must be non-zero integers,
    - generator indices must form a compatible strand set,
    - the implied strand count must be at least 1.
    """
    if any(not isinstance(entry, int) for entry in braid_word):
        raise InputValidationError("Braid entries must all be integers.")
    if any(entry == 0 for entry in braid_word):
        raise InputValidationError("Braid entries must be non-zero integers (Artin generators).")
    
    crossing_count = len(braid_word)
    if max_crossings > 0 and crossing_count > max_crossings:
        raise CrossingLimitError(f"Braid word length {crossing_count} exceeds limit {max_crossings}.")

    if not braid_word:
        return tuple()

    max_index = max(abs(entry) for entry in braid_word)
    if max_index < 1:
        raise InputValidationError("Braid generator indices must start at 1.")

    # This implementation assumes an m-strand braid with generators σ_i acting on strands i and i+1.
    strand_count = max_index + 1
    if strand_count < 2:
        raise InputValidationError("Non-empty braid words require at least 2 strands.")
    return tuple(braid_word)


def count_loops_in_closure(state: int, braid_word: Sequence[int], strand_count: int) -> int:
    """Count loops in a smoothing state of a closed braid diagram.

    Each crossing gets four ports: TL, TR, BL, BR. The smoothing choice is encoded by the bit of
    ``state`` corresponding to that crossing. The orientation-dependent smoothing convention here
    matches the original calculator logic and reproduces the standard trefoil polynomial under the
    standard braid word σ₁³.
    """
    crossing_count = len(braid_word)
    if crossing_count == 0:
        return strand_count

    uf = UnionFind(4 * crossing_count)
    current_port: List[int | None] = [None] * strand_count
    top_port: List[int | None] = [None] * strand_count

    for crossing_index, crossing in enumerate(braid_word):
        strand = abs(crossing) - 1
        tl, tr, bl, br = 4 * crossing_index, 4 * crossing_index + 1, 4 * crossing_index + 2, 4 * crossing_index + 3

        if current_port[strand] is None:
            top_port[strand] = tl
        else:
            uf.union(current_port[strand], tl)

        if current_port[strand + 1] is None:
            top_port[strand + 1] = tr
        else:
            uf.union(current_port[strand + 1], tr)

        current_port[strand] = bl
        current_port[strand + 1] = br

        bit = (state >> crossing_index) & 1
        is_positive = crossing > 0
        apply_vertical = is_positive != bool(bit)

        if apply_vertical:
            uf.union(tl, bl)
            uf.union(tr, br)
        else:
            uf.union(tl, tr)
            uf.union(bl, br)

    untouched_strands = 0
    for strand in range(strand_count):
        if current_port[strand] is not None:
            assert top_port[strand] is not None
            uf.union(current_port[strand], top_port[strand])
        else:
            untouched_strands += 1

    component_count = len({uf.find(index) for index in range(4 * crossing_count)})
    return component_count + untouched_strands


def _state_worker(args: Tuple[int, int, Sequence[int], int]) -> Dict[Tuple[int, int], int]:
    """Worker function for parallel state-sum computation.
    Returns a mapping of (zeros - ones, loop_count) -> occurrence_count.
    """
    start_state, end_state, braid_word, strand_count = args
    crossing_count = len(braid_word)
    counts: DefaultDict[Tuple[int, int], int] = defaultdict(int)
    for state in range(start_state, end_state):
        loop_count = count_loops_in_closure(state, braid_word, strand_count)
        ones = state.bit_count() if hasattr(int, "bit_count") else bin(state).count("1")
        zeros = crossing_count - ones
        counts[(zeros - ones, loop_count)] += 1
    return dict(counts)


def compute_kauffman_bracket_A(
    braid_word: Sequence[int], max_crossings: int = 0, parallel: bool = True
) -> Tuple[LaurentPolynomial, int, int]:
    """Compute the normalized Kauffman bracket in the A variable.

    Returns
    -------
    normalized_A_polynomial, strand_count, writhe
    """
    braid_word = validate_braid_word(braid_word, max_crossings)
    if not braid_word:
        return LaurentPolynomial({0: 1}), 1, 0

    crossing_count = len(braid_word)
    strand_count = max(abs(crossing) for crossing in braid_word) + 1
    num_states = 1 << crossing_count
    
    # Use multiprocessing for complex braids (>= 12 crossings)
    num_procs = multiprocessing.cpu_count() if parallel and crossing_count >= 12 else 1
    full_counts: DefaultDict[Tuple[int, int], int] = defaultdict(int)

    if num_procs > 1:
        chunk_size = num_states // num_procs
        worker_args = []
        for i in range(num_procs):
            start = i * chunk_size
            end = num_states if i == num_procs - 1 else (i + 1) * chunk_size
            worker_args.append((start, end, braid_word, strand_count))
        
        with multiprocessing.Pool(num_procs) as pool:
            results = pool.map(_state_worker, worker_args)
            for r in results:
                for key, val in r.items():
                    full_counts[key] += val
    else:
        res = _state_worker((0, num_states, braid_word, strand_count))
        for key, val in res.items():
            full_counts[key] += val

    d = LaurentPolynomial({2: -1, -2: -1})
    d_powers = [LaurentPolynomial({0: 1})]
    # Pre-calculate powers of d up to strand_count + crossing_count
    for _ in range(crossing_count + strand_count + 1):
        d_powers.append(d_powers[-1] * d)

    bracket = LaurentPolynomial()
    for (weight_pow, loop_count), count in full_counts.items():
        weight = LaurentPolynomial({weight_pow: count})
        bracket = bracket + (weight * d_powers[loop_count - 1])

    writhe = sum(1 if crossing > 0 else -1 for crossing in braid_word)
    normalized = bracket.monomial_mul(-3 * writhe, (-1) ** writhe)
    return normalized, strand_count, writhe


def convert_A_polynomial_to_q(normalized_A: LaurentPolynomial) -> LaurentPolynomial:
    """Convert A-exponents to q-exponents using q = A^(-4).

    If a monomial is A^p, then its q exponent is -p/4.
    """
    q_poly = LaurentPolynomial()
    for power_A, coefficient in normalized_A.terms.items():
        q_power = -power_A / 4
        q_poly.terms[q_power] += coefficient
    q_poly._clean()
    return q_poly


def compute_jones_result(braid_word: Sequence[int], max_crossings: int = 0, parallel: bool = True) -> JonesResult:
    normalized_A, strand_count, writhe = compute_kauffman_bracket_A(braid_word, max_crossings, parallel=parallel)
    normalized_q = convert_A_polynomial_to_q(normalized_A)
    q_value = normalized_q.evaluate(FIFTH_ROOT_OF_UNITY)

    return JonesResult(
        braid_word=tuple(braid_word),
        strand_count=strand_count,
        writhe=writhe,
        kaufman_bracket_A=normalized_A.as_serializable_dict(),
        normalized_jones_q=normalized_q.as_serializable_dict(),
        normalized_jones_pretty=normalized_q.to_pretty_string("q"),
        normalized_jones_latex=normalized_q.to_latex("q"),
        value_at_q5_real=float(q_value.real),
        value_at_q5_imag=float(q_value.imag),
        magnitude_at_q5=float(abs(q_value)),
    )


def compute_t_min(g_gap: float) -> int:
    if not (0 < g_gap < 1):
        raise InputValidationError("G_gap must satisfy 0 < G_gap < 1 to compute T_min.")
    return math.ceil((53.0 * math.log(2.0)) / abs(math.log(g_gap)))


def compute_lcl833_metrics(
    *,
    v_j: float,
    r_kh: float = 0.0,
    genus: int = 5,
    gamma: float = DEFAULT_GAMMA,
    alpha_ref: float = TREFOIL_ALPHA_OP,
    delta_ref: float = TREFOIL_MAGNITUDE,
) -> LCL833Metrics:
    """Compute the calibrated LCL-833 model metrics.

    The manuscript defines:
        δ_eff(K) = v_J(K) + γ R_Kh(K)
        α_op(g, K) = α_ref(g) * δ_eff(K) / δ_eff(3_1)
        G_gap(g, K) = 1 - α_op(g, K)
        ε_eff = (1 - α_op)/2  under the Pauli-diagonal logical model
        ω(g, K) = (g - 1) α_op(g, K)
        T_min = ceil(53 ln 2 / |ln G_gap|)
    with the trefoil anchor α_op(5,3_1)=0.8783 and v_J(3_1)=1.543361918426817.
    """
    if genus < 1:
        raise InputValidationError("Genus must be a positive integer.")
    if gamma < 0:
        raise InputValidationError("Gamma must be non-negative.")
    if delta_ref <= 0:
        raise InputValidationError("delta_ref must be strictly positive.")
    if alpha_ref <= 0:
        raise InputValidationError("alpha_ref must be strictly positive.")

    delta_eff = float(v_j + gamma * r_kh)
    alpha_op = float(alpha_ref * (delta_eff / delta_ref))
    g_gap = float(1.0 - alpha_op)
    eps_eff = float((1.0 - alpha_op) / 2.0)
    omega = float((genus - 1) * alpha_op)

    if not (0 < alpha_op < 1):
        raise InputValidationError(
            "Calibrated α_op is outside (0,1). Adjust v_J, R_Kh, gamma, or reference calibration."
        )
    if not (0 < g_gap < 1):
        raise InputValidationError(
            "Derived G_gap is outside (0,1), so the contraction-based T_min law is not valid here."
        )

    t_min = compute_t_min(g_gap)

    return LCL833Metrics(
        genus=genus,
        gamma=float(gamma),
        v_j=float(v_j),
        r_kh=float(r_kh),
        delta_eff=delta_eff,
        alpha_ref=float(alpha_ref),
        alpha_op=alpha_op,
        g_gap=g_gap,
        eps_eff=eps_eff,
        omega=omega,
        t_min=t_min,
    )


def build_report(jones: JonesResult, metrics: LCL833Metrics | None = None) -> Dict[str, object]:
    payload: Dict[str, object] = {"jones": asdict(jones)}
    if metrics is not None:
        payload["lcl833"] = asdict(metrics)
    return payload


def print_human_report(jones: JonesResult, metrics: LCL833Metrics | None = None) -> None:
    print("=" * 72)
    print("LCL-833 MASTER")
    print("=" * 72)
    print(f"Braid word            : {list(jones.braid_word)}")
    print(f"Strands               : {jones.strand_count}")
    print(f"Writhe                : {jones.writhe}")
    print(f"Jones polynomial      : {jones.normalized_jones_pretty}")
    print(f"q = exp(2πi/5)        : {FIFTH_ROOT_OF_UNITY.real:.15f} + {FIFTH_ROOT_OF_UNITY.imag:.15f}i")
    print(f"J(q)                  : {jones.value_at_q5_real:.15f} + {jones.value_at_q5_imag:.15f}i")
    print(f"|J(q)|                : {jones.magnitude_at_q5:.15f}")

    if metrics is not None:
        print("-" * 72)
        print("LCL-833 calibrated metrics")
        print(f"genus                 : {metrics.genus}")
        print(f"γ                     : {metrics.gamma:.6f}")
        print(f"R_Kh                  : {metrics.r_kh:.6f}")
        print(f"δ_eff                 : {metrics.delta_eff:.15f}")
        print(f"α_ref                 : {metrics.alpha_ref:.6f}")
        print(f"α_op                  : {metrics.alpha_op:.15f}")
        print(f"G_gap                 : {metrics.g_gap:.15f}")
        print(f"ε_eff                 : {metrics.eps_eff:.15f}")
        print(f"ω                     : {metrics.omega:.15f}")
        print(f"T_min                 : {metrics.t_min}")
    print("=" * 72)


def run_selftest() -> None:
    """Minimal regression suite for production deployment."""
    LOGGER.info("Running self-tests.")

    trefoil = compute_jones_result((1, 1, 1))
    expected_terms = {"4": -1, "3": 1, "1": 1}
    if trefoil.normalized_jones_q != expected_terms:
        raise AssertionError(f"Trefoil polynomial mismatch: {trefoil.normalized_jones_q} != {expected_terms}")
    if abs(trefoil.magnitude_at_q5 - TREFOIL_MAGNITUDE) > 1e-12:
        raise AssertionError(
            f"Trefoil magnitude mismatch: {trefoil.magnitude_at_q5} != {TREFOIL_MAGNITUDE}"
        )

    metrics = compute_lcl833_metrics(v_j=trefoil.magnitude_at_q5, r_kh=0.0, genus=5)
    if abs(metrics.alpha_op - TREFOIL_ALPHA_OP) > 1e-12:
        raise AssertionError(f"Trefoil alpha mismatch: {metrics.alpha_op} != {TREFOIL_ALPHA_OP}")
    if abs(metrics.g_gap - 0.1217) > 1e-12:
        raise AssertionError(f"Trefoil G_gap mismatch: {metrics.g_gap} != 0.1217")
    if metrics.t_min != 18:
        raise AssertionError(f"Trefoil T_min mismatch: {metrics.t_min} != 18")

    unknot = compute_jones_result(())
    if unknot.normalized_jones_q != {"0": 1}:
        raise AssertionError(f"Unknot mismatch: {unknot.normalized_jones_q}")

    figure8 = compute_jones_result((1, -2, 1, -2))
    if figure8.normalized_jones_q != {"2": 1, "1": -1, "0": 1, "-1": -1, "-2": 1}:
        raise AssertionError(f"Figure-8 polynomial mismatch: {figure8.normalized_jones_q}")

    LOGGER.info("All self-tests passed.")
    print("SELFTEST: PASS")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="LCL-833 Jones / calibrated metrics master tool")
    parser.add_argument("--log-level", default="WARNING", help="Logging level (DEBUG, INFO, WARNING, ERROR).")

    subparsers = parser.add_subparsers(dest="command", required=True)

    # Common arguments for knot/braid
    def add_common_args(p: argparse.ArgumentParser):
        p.add_argument("--json", action="store_true", help="Emit JSON instead of a human-readable report.")
        p.add_argument("--genus", type=int, default=5)
        p.add_argument("--gamma", type=float, default=DEFAULT_GAMMA)
        p.add_argument("--rkh", type=float, default=None, help="Override reduced Khovanov rank ratio R_Kh.")
        p.add_argument("--no-lcl", action="store_true", help="Skip LCL-833 calibrated metrics.")
        p.add_argument("--max-crossings", type=int, default=MAX_CROSSINGS_DEFAULT)
        p.add_argument("--parallel", action="store_true", default=True, help="Use multiprocessing for large braids (default: True).")
        p.add_argument("--no-parallel", action="store_false", dest="parallel", help="Disable multiprocessing.")

    preset_parser = subparsers.add_parser("knot", help="Run a preset knot by name.")
    preset_parser.add_argument("name", choices=sorted(PRESET_KNOTS.keys()))
    add_common_args(preset_parser)

    braid_parser = subparsers.add_parser("braid", help="Run a custom braid word.")
    braid_parser.add_argument("entries", nargs="*", type=int, help="Artin generator list, e.g. 1 1 1 or 1 -2 1 -2")
    add_common_args(braid_parser)

    table_parser = subparsers.add_parser("table", help="Output comparison table for all preset knots.")
    table_parser.add_argument("--json", action="store_true", help="Emit JSON comparison.")
    table_parser.add_argument("--genus", type=int, default=5)
    table_parser.add_argument("--gamma", type=float, default=DEFAULT_GAMMA)

    subparsers.add_parser("selftest", help="Run the regression test suite.")
    return parser.parse_args(argv)


def configure_logging(log_level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, log_level.upper(), logging.WARNING),
        format="%(levelname)s: %(message)s",
    )


def execute_from_args(args: argparse.Namespace) -> int:
    if args.command == "selftest":
        run_selftest()
        return 0

    if args.command == "table":
        results = {}
        for name, word in PRESET_KNOTS.items():
            try:
                jones = compute_jones_result(word)
                metrics = compute_lcl833_metrics(
                    v_j=jones.magnitude_at_q5, 
                    r_kh=PRESET_RKH.get(name, 0.0), 
                    genus=args.genus, 
                    gamma=args.gamma
                )
                results[name] = build_report(jones, metrics)
            except InputValidationError as exc:
                LOGGER.warning(f"Skipping {name} in table: {exc}")
                continue
        
        if args.json:
            print(json.dumps(results, indent=2, sort_keys=True))
        else:
            print(f"{'Knot':<20} | {'|J(q5)|':<15} | {'α_op':<10} | {'T_min':<5}")
            print("-" * 60)
            for name in sorted(results.keys()):
                r = results[name]
                vj = r["jones"]["magnitude_at_q5"]
                alpha = r["lcl833"]["alpha_op"]
                tmin = r["lcl833"]["t_min"]
                print(f"{name:<20} | {vj:<15.6f} | {alpha:<10.4f} | {tmin:<5}")
        return 0

    if args.command == "knot":
        braid_word = PRESET_KNOTS[args.name]
        r_kh = PRESET_RKH.get(args.name, 0.0) if args.rkh is None else args.rkh
    elif args.command == "braid":
        braid_word = tuple(args.entries)
        r_kh = args.rkh or 0.0
    else:
        raise InputValidationError(f"Unsupported command: {args.command}")

    jones = compute_jones_result(braid_word, max_crossings=args.max_crossings, parallel=args.parallel)
    metrics = None if args.no_lcl else compute_lcl833_metrics(
        v_j=jones.magnitude_at_q5, r_kh=r_kh, genus=args.genus, gamma=args.gamma
    )

    if args.json:
        print(json.dumps(build_report(jones, metrics), indent=2, sort_keys=True))
    else:
        print_human_report(jones, metrics)
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parse_args(argv)
    configure_logging(args.log_level)
    try:
        return execute_from_args(args)
    except (InputValidationError, CrossingLimitError) as exc:
        LOGGER.error(str(exc))
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    except LCL833Error as exc:
        LOGGER.error(str(exc))
        print(f"ERROR: {exc}", file=sys.stderr)
        return 3


if __name__ == "__main__":
    raise SystemExit(main())
