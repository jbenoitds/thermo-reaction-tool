import pytest

from thermoreact.formula import FormulaError, hill_formula, molar_mass, parse_formula
from thermoreact.smiles import SmilesError, parse_smiles


@pytest.mark.parametrize("text, elements, charge", [
    ("H2O", {"H": 2, "O": 1}, 0),
    ("Ca(OH)2", {"Ca": 1, "O": 2, "H": 2}, 0),
    ("Al2(SO4)3", {"Al": 2, "S": 3, "O": 12}, 0),
    ("K4[Fe(CN)6]", {"K": 4, "Fe": 1, "C": 6, "N": 6}, 0),
    ("CuSO4·5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}, 0),
    ("CuSO4*5H2O", {"Cu": 1, "S": 1, "O": 9, "H": 10}, 0),
    ("C₆H₁₂O₆", {"C": 6, "H": 12, "O": 6}, 0),
    ("SO4^2-", {"S": 1, "O": 4}, -2),
    ("Fe^3+", {"Fe": 1}, 3),
    ("Cu++", {"Cu": 1}, 2),
    ("Cu+2", {"Cu": 1}, 2),
    ("Cu[2+]", {"Cu": 1}, 2),
    ("NO3-", {"N": 1, "O": 3}, -1),
    ("NH4+", {"N": 1, "H": 4}, 1),
])
def test_parse_formula(text, elements, charge):
    assert parse_formula(text) == (elements, charge)


@pytest.mark.parametrize("bad", ["", "H2Xx", "Ca(OH", "h2o", "Ca)OH("])
def test_bad_formulas(bad):
    with pytest.raises(FormulaError):
        parse_formula(bad)


def test_hill_and_mass():
    assert hill_formula({"O": 1, "H": 6, "C": 2}) == "C2H6O"
    assert hill_formula({"O": 4, "S": 1}, -2) == "O4S^2-"
    assert hill_formula({"Cl": 1, "Na": 1}) == "ClNa"
    assert molar_mass(parse_formula("H2O")[0]) == pytest.approx(18.015, abs=1e-3)
    assert molar_mass(parse_formula("CuSO4·5H2O")[0]) == pytest.approx(249.68, abs=0.01)


@pytest.mark.parametrize("smiles, formula", [
    ("C", "CH4"), ("CCO", "C2H6O"), ("OCC", "C2H6O"), ("C=C", "C2H4"), ("C#C", "C2H2"),
    ("c1ccccc1", "C6H6"), ("C1=CC=CC=C1", "C6H6"), ("Oc1ccccc1", "C6H6O"),
    ("CC(=O)O", "C2H4O2"), ("c1ccncc1", "C5H5N"), ("c1cc[nH]c1", "C4H5N"), ("c1ccoc1", "C4H4O"),
    ("C[N+](=O)[O-]", "CH3NO2"), ("[NH4+]", "H4N^+"), ("[13CH4]", "CH4"), ("ClC(Cl)(Cl)Cl", "CCl4"),
    ("N[C@@H](C)C(=O)O", "C3H7NO2"), ("C1CC2CCC1CC2", "C8H14"), ("C%10CCCCC%10", "C6H12"),
    ("[Na+].[Cl-]", "ClNa"), ("O=C=O", "CO2"), ("FC(F)=C(F)F", "C2F4"), ("OS(=O)(=O)O", "H2O4S"),
])
def test_smiles_formula(smiles, formula):
    assert parse_smiles(smiles).formula == formula


def test_smiles_rings():
    m = parse_smiles("CC1CCCC1")
    assert not m.in_ring(0) and all(m.in_ring(i) for i in range(1, 6))
    m = parse_smiles("c1ccccc1-c1ccccc1")  # biphenyl: inter-ring bond is not a ring bond
    assert not m.bond_in_ring(5, 6)


@pytest.mark.parametrize("bad", ["C1CC", "C(C", "C)C", "[Xx]", "C^C", ""])
def test_bad_smiles(bad):
    with pytest.raises(SmilesError):
        parse_smiles(bad)


def test_signature_distinguishes_isomers():
    assert parse_smiles("CCO").signature() == parse_smiles("OCC").signature()
    assert parse_smiles("CCO").signature() != parse_smiles("COC").signature()
