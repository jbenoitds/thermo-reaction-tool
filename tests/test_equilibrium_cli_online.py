import json
import urllib.error

import pytest

from thermoreact import Database, Reaction, solve_equilibrium
from thermoreact import online
from thermoreact.cli import main
from thermoreact.database import SpeciesNotFound


# ----------------------------------------------------------------------------- equilibrium
def test_equilibrium_satisfies_K():
    rxn = Reaction.from_string("N2 + 3 H2 -> 2 NH3")
    res = solve_equilibrium(rxn, {"N2": 1, "H2": 3}, T=700, P_bar=200)
    x = res.mole_fractions()
    Q = x["NH3(g)"] ** 2 / (x["N2(g)"] * x["H2(g)"] ** 3) * 200.0 ** -2
    assert Q == pytest.approx(res.K, rel=1e-9)
    assert sum(x.values()) == pytest.approx(1)
    assert res.amounts["N2(g)"] * 3 == pytest.approx(res.amounts["H2(g)"])  # stoichiometric feed stays so


def test_le_chatelier_pressure_and_inert():
    rxn = Reaction.from_string("N2 + 3 H2 -> 2 NH3")
    lo = solve_equilibrium(rxn, {"N2": 1, "H2": 3}, T=700, P_bar=50).extent
    hi = solve_equilibrium(rxn, {"N2": 1, "H2": 3}, T=700, P_bar=300).extent
    diluted = solve_equilibrium(rxn, {"N2": 1, "H2": 3}, T=700, P_bar=50, inert_gas=4).extent
    assert hi > lo > diluted


def test_equilibrium_with_condensed_phase():
    rxn = Reaction.from_string("C(graphite) + CO2 -> 2 CO")
    res = solve_equilibrium(rxn, {"C(graphite)": 5, "CO2": 1}, T=1000, P_bar=1)
    x = res.mole_fractions()
    assert x["CO(g)"] ** 2 / x["CO2(g)"] == pytest.approx(res.K, rel=1e-9)
    assert "C(graphite)" not in x


def test_equilibrium_condensed_species_runs_out():
    rxn = Reaction.from_string("C(graphite) + CO2 -> 2 CO")
    res = solve_equilibrium(rxn, {"C(graphite)": 0.01, "CO2": 1}, T=1300, P_bar=1)
    assert res.amounts["C(graphite)"] == pytest.approx(0, abs=1e-12)
    assert any("used up" in n for n in res.notes)


def test_equilibrium_extreme_K():
    rxn = Reaction.from_string("2 H2 + O2 -> 2 H2O(g)")
    res = solve_equilibrium(rxn, {"H2": 2, "O2": 1}, T=298.15)
    assert res.amounts["H2O(g)"] == pytest.approx(2, rel=1e-9)


def test_decomposition_degenerate_case():
    rxn = Reaction.from_string("CaCO3 -> CaO + CO2")
    assert solve_equilibrium(rxn, {"CaCO3": 1}, T=1200).amounts["CaO(s)"] == pytest.approx(1)
    assert solve_equilibrium(rxn, {"CaCO3": 1}, T=900).amounts["CaO(s)"] == pytest.approx(0)


# ----------------------------------------------------------------------------- CLI
def run(capsys, *argv):
    code = main(list(argv))
    out = capsys.readouterr()
    return code, out.out, out.err


def test_cli_reaction_table(capsys):
    code, out, _ = run(capsys, "reaction", "N2 + 3 H2 -> 2 NH3", "-T", "298.15", "700", "--crossover")
    assert code == 0
    assert "-91.800" in out and "455." in out


def test_cli_reaction_json_and_csv(capsys, tmp_path):
    csv_path = tmp_path / "t.csv"
    code, out, _ = run(capsys, "rxn", "CO + H2O(g) -> CO2 + H2", "--trange", "300", "1100", "200",
                       "--json", "--csv", str(csv_path))
    data = json.loads(out)
    assert code == 0 and len(data["results"]) == 5
    assert data["results"][0]["T_K"] == 300
    assert csv_path.read_text().count("\n") == 6


def test_cli_celsius(capsys):
    code, out, _ = run(capsys, "rxn", "N2 + 3 H2 -> 2 NH3", "-T", "25", "--celsius", "--json")
    assert json.loads(out)["results"][0]["T_K"] == pytest.approx(298.15)


def test_cli_plot(capsys, tmp_path):
    pytest.importorskip("matplotlib")
    png = tmp_path / "p.png"
    code, _, _ = run(capsys, "rxn", "CaCO3 -> CaO + CO2", "--trange", "800", "1300", "50", "--plot", str(png))
    assert code == 0 and png.stat().st_size > 10000


def test_cli_errors_are_friendly(capsys):
    code, _, err = run(capsys, "rxn", "H2 + O2 -> H2O")
    assert code == 2 and "not balanced" in err
    code, _, err = run(capsys, "rxn", "Unobtainium -> H2")
    assert code == 2 and "No species" in err


def test_cli_missing_data_file(capsys):
    code, _, err = run(capsys, "rxn", "N2 + 3 H2 -> 2 NH3", "--data", "does_not_exist.csv")
    assert code == 2 and "does_not_exist.csv" in err


def test_cli_output_redirected_with_legacy_encoding(tmp_path):
    """Windows: 'thermoreact ... > out.txt' uses cp1252, which has no Δ or →."""
    import os
    import subprocess
    import sys

    env = dict(os.environ, PYTHONIOENCODING="cp1252")
    out = tmp_path / "out.txt"
    with open(out, "wb") as fh:
        proc = subprocess.run([sys.executable, "-m", "thermoreact", "rxn", "N2 + 3 H2 -> 2 NH3"],
                              stdout=fh, stderr=subprocess.PIPE, env=env)
    assert proc.returncode == 0, proc.stderr
    assert "ΔrH°" in out.read_text(encoding="utf-8")


def test_cli_other_commands(capsys):
    assert run(capsys, "species", "NH3", "-T", "500", "1000")[1].count("NH3") >= 1
    assert "C2H5OH(l)" in run(capsys, "list", "ethanol")[1]
    assert "-OH (alcohol)" in run(capsys, "estimate", "CCO")[1]
    out = run(capsys, "eq", "N2O4 -> 2 NO2", "--init", "N2O4=1", "--extent", "0.201")[1]
    assert "0.16841" in out
    out = run(capsys, "eq", "A + B -> C + 2 D", "--abstract", "--init", "A=2,B=1,D=3", "--K", "9.58")[1]
    assert "C(g)" in out
    assert "2.2298" in run(capsys, "vanthoff", "--K1", "0.414", "--T1", "1600", "--T2", "2000", "--dH", "112")[1]
    assert "-38.7" in run(capsys, "vanthoff", "--K1", "1", "--T1", "310", "--K2", "0.5", "--T2", "325")[1]
    assert "1.688e+16" in run(capsys, "electrochem", "--E0", "0.48", "--n", "2")[1]


# ----------------------------------------------------------------------------- online (mocked)
@pytest.fixture
def fake_pubchem(monkeypatch):
    calls = []

    def fetch(url):
        calls.append(url)
        if "isopentane" in url:
            return json.dumps({"PropertyTable": {"Properties": [
                {"CID": 6556, "MolecularFormula": "C5H12", "IsomericSMILES": "CCC(C)C",
                 "IUPACName": "2-methylbutane"}]}}).encode()
        if "newkeys" in url:  # newer PubChem responses use 'SMILES'
            return json.dumps({"PropertyTable": {"Properties": [{"CID": 1, "SMILES": "CCCC=O"}]}}).encode()
        raise urllib.error.HTTPError(url, 404, "Not Found", {}, None)

    monkeypatch.setattr(online, "_fetch", fetch)
    online._cache.clear()
    return calls


def test_online_name_resolution(fake_pubchem):
    db = Database.default().copy()
    sp = db.get("isopentane", online=True)
    assert sp.formula == "C5H12" and sp.estimated
    db.get("isopentane", online=True)
    assert len(fake_pubchem) == 1                         # cached
    assert online.resolve_name("newkeys")["smiles"] == "CCCC=O"
    with pytest.raises(SpeciesNotFound, match="PubChem does not know"):
        db.get("notarealmolecule", online=True)


def test_online_network_failure(monkeypatch):
    def boom(url):
        raise urllib.error.URLError("offline")

    monkeypatch.setattr(online, "_fetch", boom)
    online._cache.clear()
    with pytest.raises(SpeciesNotFound, match="Could not reach PubChem"):
        Database.default().copy().get("something", online=True)
