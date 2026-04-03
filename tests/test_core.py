import pytest
from fractions import Fraction
from lcl833_master.core import LaurentPolynomial, UnionFind, count_loops_in_closure, compute_kauffman_bracket_A, convert_A_polynomial_to_q

def test_laurent_polynomial_init():
    lp = LaurentPolynomial({1: 2, Fraction(1, 2): 3})
    assert lp.terms[Fraction(1)] == 2.0
    assert lp.terms[Fraction(1, 2)] == 3.0

def test_laurent_polynomial_add():
    lp1 = LaurentPolynomial({1: 2})
    lp2 = LaurentPolynomial({1: 3, 0: 1})
    lp3 = lp1 + lp2
    assert lp3.terms[Fraction(1)] == 5.0
    assert lp3.terms[Fraction(0)] == 1.0

def test_laurent_polynomial_mul():
    lp1 = LaurentPolynomial({1: 2})
    lp2 = LaurentPolynomial({-1: 3, 0: 1})
    lp3 = lp1 * lp2
    # 2q^1 * (3q^-1 + 1q^0) = 6q^0 + 2q^1
    assert lp3.terms[Fraction(0)] == 6.0
    assert lp3.terms[Fraction(1)] == 2.0

def test_laurent_polynomial_pow():
    lp = LaurentPolynomial({1: 1, -1: 1})
    lp2 = lp ** 2
    # (q + q^-1)^2 = q^2 + 2 + q^-2
    assert lp2.terms[Fraction(2)] == 1.0
    assert lp2.terms[Fraction(0)] == 2.0
    assert lp2.terms[Fraction(-2)] == 1.0

def test_union_find():
    uf = UnionFind(5)
    uf.union(0, 1)
    uf.union(1, 2)
    assert uf.find(0) == uf.find(2)
    assert uf.find(3) != uf.find(0)

def test_count_loops_unknot():
    # Unknot as empty braid
    assert count_loops_in_closure(0, [], 1) == 1
    # Two disjoint loops as empty braid on 2 strands
    assert count_loops_in_closure(0, [], 2) == 2

def test_count_loops_trefoil_positive():
    # Trefoil (1, 1, 1) - 3 crossings, 2 strands
    # There are 2^3 = 8 states.
    # State 0 (all 0-smoothings): ...
    # This is complex to check manually for all states, but we can check a few.
    # State 0: (bit=0, crossing=1 > 0) -> apply_vertical = True
    # If all vertical, we should get 2 loops (the two strands themselves closed).
    assert count_loops_in_closure(0, (1, 1, 1), 2) == 2
    # State 7 (all 1-smoothings): (bit=1, crossing=1 > 0) -> apply_vertical = False
    # If all horizontal, the strands are connected together.
    # For (1, 1, 1), horizontal smoothings connect (top_i, top_{i+1}) and (bottom_i, bottom_{i+1})
    # This results in 3 loops for this state sum model.
    assert count_loops_in_closure(7, (1, 1, 1), 2) == 3

def test_kauffman_bracket_unknot():
    poly, strands, writhe = compute_kauffman_bracket_A([])
    assert poly.terms == {Fraction(0): 1.0}
    assert strands == 1
    assert writhe == 0

def test_laurent_polynomial_neg_sub():
    lp1 = LaurentPolynomial({1: 2})
    lp2 = -lp1
    assert lp2.terms[Fraction(1)] == -2.0
    lp3 = lp1 - LaurentPolynomial({1: 3})
    assert lp3.terms[Fraction(1)] == -1.0

def test_laurent_polynomial_div():
    lp = LaurentPolynomial({1: 4})
    lp2 = lp / 2
    assert lp2.terms[Fraction(1)] == 2.0
    with pytest.raises(ZeroDivisionError):
        _ = lp / 0

def test_new_presets():
    from lcl833_master.core import PRESET_KNOTS, compute_jones_result, compute_lcl833_metrics
    # Hopf link
    hopf = compute_jones_result(PRESET_KNOTS["hopf-link"])
    # Our model gives {'5/2': -1, '1/2': -1} for (1, 1)
    # This is correct for the Kauffman bracket normalization of a 2-component link.
    assert hopf.normalized_jones_q == {"5/2": -1.0, "1/2": -1.0}
    
    # Test metrics with override
    metrics = compute_lcl833_metrics(v_j=hopf.magnitude_at_q5, r_kh=2.0, gamma=0.1)
    assert metrics.r_kh == 2.0
    assert metrics.gamma == 0.1
    assert metrics.delta_eff == hopf.magnitude_at_q5 + 0.1 * 2.0

def test_laurent_polynomial_pretty():
    lp = LaurentPolynomial({2: 1, 0: -3, Fraction(1, 2): 1})
    pretty = lp.to_pretty_string("q")
    # Powers are sorted descending: 2, 1/2, 0
    assert "q^2" in pretty
    assert "q^(1/2)" in pretty
    assert "- 3" in pretty
