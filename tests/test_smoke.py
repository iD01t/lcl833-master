import lcl833_master

def test_version():
    assert lcl833_master.__version__ == "0.3.0"

def test_imports():
    from lcl833_master import compute_jones_result, compute_lcl833_metrics
    assert compute_jones_result is not None
    assert compute_lcl833_metrics is not None
