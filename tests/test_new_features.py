import pytest
from lcl833_master.core import LaurentPolynomial, compute_jones_result, compute_kauffman_bracket_A

def test_laurent_to_latex():
    lp = LaurentPolynomial({2: 1, 0: -3, 0.5: 1})
    latex = lp.to_latex("q")
    # sorted descending: 2, 1/2, 0
    # Current output: 'q^{2} + q^{1/2} - 3'
    assert "q^{2}" in latex
    assert "q^{1/2}" in latex
    assert "- 3" in latex
    assert " + " in latex

def test_jones_result_includes_latex():
    res = compute_jones_result((1, 1, 1))
    assert hasattr(res, "normalized_jones_latex")
    # Trefoil latex: '- q^{4} + q^{3} + q'
    assert "q^{4}" in res.normalized_jones_latex
    assert "q^{3}" in res.normalized_jones_latex
    assert " q" in res.normalized_jones_latex

def test_parallel_computation_smoke():
    # Use a slightly larger braid to trigger parallel code path if enabled, 
    # or just ensure it doesn't crash.
    # We can force parallel=True.
    poly, strands, writhe = compute_kauffman_bracket_A((1, 1, 1), parallel=True)
    assert poly is not None
    assert strands == 2
    assert writhe == 3

def test_laurent_pow_efficiency():
    lp = LaurentPolynomial({1: 1, -1: 1})
    res = lp ** 10
    # (q + q^-1)^10 should have 11 terms
    assert len(res.terms) == 11
    # Check middle term: C(10, 5) = 252
    assert res.terms[0] == 252
