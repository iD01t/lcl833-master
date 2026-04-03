import json
import pytest
from lcl833_master.core import main

def test_cli_selftest(capsys):
    ret = main(["selftest"])
    assert ret == 0
    captured = capsys.readouterr()
    assert "SELFTEST: PASS" in captured.out

def test_cli_knot_json(capsys):
    # Fixed position: --json is now after the subparser
    ret = main(["knot", "trefoil", "--json"])
    assert ret == 0
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert "jones" in data
    assert "lcl833" in data
    assert data["jones"]["braid_word"] == [1, 1, 1]

def test_cli_braid_rkh(capsys):
    # Fixed position: --json and --rkh after 'braid'
    ret = main(["braid", "1", "1", "1", "--json", "--rkh", "1.0"])
    assert ret == 0
    captured = capsys.readouterr()
    data = json.loads(captured.out)
    assert data["lcl833"]["r_kh"] == 1.0
    assert data["lcl833"]["delta_eff"] > 1.59

def test_cli_table(capsys):
    ret = main(["table"])
    assert ret == 0
    captured = capsys.readouterr()
    assert "trefoil" in captured.out.lower()
    assert "figure8" in captured.out.lower()

def test_cli_max_crossings(capsys):
    # Cinquefoil has 5 crossings. Test with limit 2.
    ret = main(["knot", "cinquefoil", "--max-crossings", "2"])
    assert ret == 2
    captured = capsys.readouterr()
    assert "exceeds limit 2" in captured.err

def test_cli_invalid_knot(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["knot", "invalid-knot-name"])
    assert excinfo.value.code == 2
