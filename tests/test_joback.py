import pytest

from thermoreact import Database, Reaction
from thermoreact.joback import JobackError, assign_groups, estimate
from thermoreact.smiles import parse_smiles


@pytest.mark.parametrize("smiles, groups", [
    ("CCO", {"-CH3": 1, "-CH2-": 1, "-OH (alcohol)": 1}),
    ("CC(C)=O", {"-CH3": 2, ">C=O (nonring)": 1}),
    ("CC(=O)O", {"-CH3": 1, "-COOH (acid)": 1}),
    ("CCOC(C)=O", {"-CH3": 2, "-CH2-": 1, "-COO- (ester)": 1}),
    ("CC=O", {"-CH3": 1, "O=CH- (aldehyde)": 1}),
    ("Oc1ccccc1", {"-OH (phenol)": 1, "=C< (ring)": 1, "=CH- (ring)": 5}),
    ("C1=CC=C(C=C1)O", {"-OH (phenol)": 1, "=C< (ring)": 1, "=CH- (ring)": 5}),  # Kekulé form
    ("C1CCCCC1", {"-CH2- (ring)": 6}),
    ("C1CCOC1", {"-CH2- (ring)": 4, "-O- (ring)": 1}),
    ("O=C1CCCCC1", {">C=O (ring)": 1, "-CH2- (ring)": 5}),
    ("CC#N", {"-CH3": 1, "-CN": 1}),
    ("C[N+](=O)[O-]", {"-CH3": 1, "-NO2": 1}),
    ("CN(C)C", {"-CH3": 3, ">N- (nonring)": 1}),
    ("c1ccncc1", {"=CH- (ring)": 5, "-N= (ring)": 1}),
    ("c1cc[nH]c1", {"=CH- (ring)": 4, ">NH (ring)": 1}),
    ("CSC", {"-CH3": 2, "-S- (nonring)": 1}),
    ("C=C=C", {"=CH2": 2, "=C=": 1}),
    ("CC#C", {"-CH3": 1, "≡C-": 1, "≡CH": 1}),
    ("CC(C)(C)C", {"-CH3": 4, ">C<": 1}),
    ("FC(F)(F)Cl", {"-F": 3, "-Cl": 1, ">C<": 1}),
])
def test_group_assignment(smiles, groups):
    assert assign_groups(parse_smiles(smiles)) == groups


@pytest.mark.parametrize("smiles, dHf_exp, tol", [
    # experimental gas-phase ΔfH°(298 K), kJ/mol (NIST WebBook); Joback is typically within a few kJ/mol
    ("CCC", -104.7, 3), ("CCCCCC", -166.9, 3), ("CCO", -234.8, 4), ("CC(C)=O", -217.1, 3),
    ("CC(=O)O", -432.8, 4), ("c1ccccc1", 82.9, 4), ("Cc1ccccc1", 50.1, 3), ("Oc1ccccc1", -96.4, 3),
    ("CCOC(C)=O", -443.6, 3), ("ClCCl", -95.4, 3), ("C[N+](=O)[O-]", -74.3, 3), ("c1cc[nH]c1", 108.3, 3),
])
def test_joback_enthalpy_vs_experiment(smiles, dHf_exp, tol):
    assert estimate(smiles).dHf_kJ == pytest.approx(dHf_exp, abs=tol)


def test_joback_reference_numbers():
    """Exact arithmetic of the method for ethanol (hand calculation)."""
    r = estimate("CCO")
    assert r.dHf_kJ == pytest.approx(68.29 - 76.45 - 20.64 - 208.04)
    assert r.dGf_kJ == pytest.approx(53.88 - 43.96 + 8.42 - 189.20)
    a, b, c, d = r.cp_coeffs
    assert a == pytest.approx(19.5 - 0.909 + 25.7 - 37.93)
    assert d == pytest.approx(-9.67e-8 + 1.19e-8 - 9.88e-8 + 2.06e-7)
    assert r.S_J == pytest.approx(281.6, abs=5)            # experimental S°(g) = 281.6
    assert r.cp(298.15) == pytest.approx(65.6, abs=2.5)     # experimental Cp(g) = 65.6


@pytest.mark.parametrize("smiles", ["C", "O", "N", "[Na+]", "C=O", "CS(C)=O", "[H][H]", "P"])
def test_joback_refuses_unsupported(smiles):
    with pytest.raises(JobackError):
        estimate(smiles)


def test_smiles_species_use_database_when_known():
    db = Database.default()
    assert db.get("smiles:CCO").id == "C2H5OH(l)"          # standard state at 298 K
    assert db.get("smiles:OCC(g)").id == "C2H5OH(g)"
    assert db.get("smiles:COC").id == "CH3OCH3(g)"          # isomer distinguished
    assert db.get("smiles:O=C=O").id == "CO2(g)"


def test_smiles_species_estimated_when_unknown():
    db = Database.default().copy()
    sp = db.get("smiles:CC(C)CC")                          # isopentane is not tabulated
    assert sp.estimated and sp.phase == "g" and sp.formula == "C5H12"
    assert sp.dHf / 1000 == pytest.approx(-153.7, abs=5)   # exp. -153.6
    rxn = Reaction.from_string("smiles:CC(C)CC + 8 O2 -> 5 CO2 + 6 H2O(l)", db)
    # experimental ΔcH° of gaseous isopentane ≈ -3529 kJ/mol (liquid: -3504; ΔvapH ≈ 25)
    assert rxn.properties(298.15).dH / 1000 == pytest.approx(-3529, abs=6)
