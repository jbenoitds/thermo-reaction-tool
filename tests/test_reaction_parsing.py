import pytest

from thermoreact import AmbiguousSpecies, Database, Reaction, ReactionError, SpeciesNotFound

DB = Database.default()


@pytest.mark.parametrize("text", [
    "CH4 + 2 O2 -> CO2 + 2 H2O(l)",
    "CH4 + 2O2 → CO2 + 2H2O(l)",
    "CH4+2O2=>CO2+2H2O(l)",
    "methane + 2 oxygen <=> carbon dioxide + 2 water",
    "CH4(g) + 2 O2(g) ⇌ CO2(g) + 2 H2O(l)",
    "1/2 CH4 + O2 -> 1/2 CO2 + H2O(l)",
    "0.5 CH4 + O2 = 0.5 CO2 + H2O(l)",
])
def test_equivalent_spellings(text):
    rxn = Reaction.from_string(text)
    ratio = float(rxn.nu("CO2(g)"))
    assert ratio > 0
    assert rxn.properties().dH / 1000 / ratio == pytest.approx(-890.31, abs=0.01)


@pytest.mark.parametrize("text, expected", [
    ("CH4 + O2 -> CO2 + H2O(g)", "CH4(g) + 2 O2(g) → CO2(g) + 2 H2O(g)"),
    ("C3H8 + O2 -> CO2 + H2O(g)", "C3H8(g) + 5 O2(g) → 3 CO2(g) + 4 H2O(g)"),
    ("Fe2O3 + CO -> Fe + CO2", "Fe2O3(s) + 3 CO(g) → 2 Fe(s) + 3 CO2(g)"),
    ("NH3 + O2 -> NO + H2O(g)", "4 NH3(g) + 5 O2(g) → 4 NO(g) + 6 H2O(g)"),
    ("C6H12O6 + O2 -> CO2 + H2O", "C6H12O6(s) + 6 O2(g) → 6 CO2(g) + 6 H2O(l)"),
    ("Cu + Ag+ -> Cu2+ + Ag", "Cu(s) + 2 Ag+(aq) → Cu2+(aq) + 2 Ag(s)"),
    ("Al + O2 -> Al2O3", "4 Al(s) + 3 O2(g) → 2 Al2O3(s)"),
])
def test_auto_balance(text, expected):
    assert Reaction.from_string(text, balance=True).equation() == expected


def test_unbalanced_is_reported():
    with pytest.raises(ReactionError, match="not balanced"):
        Reaction.from_string("H2 + O2 -> H2O")
    with pytest.raises(ReactionError, match="charge"):
        Reaction.from_string("Cu + Ag+ -> Cu2+ + Ag")


def test_unbalanceable():
    with pytest.raises(ReactionError):
        Reaction.from_string("H2 -> O2", balance=True)
    with pytest.raises(ReactionError, match="not unique"):
        Reaction.from_string("CO + O2 + H2 -> CO2 + H2O(g)", balance=True)


def test_species_lookup():
    assert DB.get("H2O").id == "H2O(l)"                    # standard state at 298 K
    assert DB.get("H2O(g)").id == "H2O(g)"
    assert DB.get("steam").id == "H2O(g)"
    assert DB.get("C").id == "C(graphite)"
    assert DB.get("C(diamond)").id == "C(diamond)"
    assert DB.get("C2H5OH").id == "C2H5OH(l)"
    assert DB.get("CuSO4.5H2O").id == "CuSO4·5H2O(s)"
    assert DB.get("SO4 2-").id == "SO4^2-(aq)"
    assert DB.get("Cu2+").id == "Cu2+(aq)"
    assert DB.get("Calcite").id == "CaCO3(s)"
    with pytest.raises(AmbiguousSpecies):
        DB.get("C2H6O")                                     # ethanol or dimethyl ether?
    with pytest.raises(SpeciesNotFound, match="Did you mean"):
        DB.get("glucoze")
    with pytest.raises(SpeciesNotFound, match="phase"):
        DB.get("NaCl(g)")


def test_same_species_on_both_sides_cancels():
    rxn = Reaction.from_string("CO + H2O(g) + H2O(g) -> CO2 + H2 + H2O(g)")
    assert float(rxn.nu("H2O(g)")) == -1
    with pytest.raises(ReactionError, match="no net change"):
        Reaction.from_string("H2O(l) + CO2 -> H2O(l) + CO2")


def test_user_data_file(tmp_path):
    f = tmp_path / "extra.csv"
    f.write_text("# my data\nformula,phase,name,dHf_kJ,S_J,Cp_J\nXeF4,s,xenon tetrafluoride,-267,146,118\n"
                 "Xe,g,xenon,0,169.68,20.79\n", encoding="utf-8")
    db = DB.copy()
    added = db.load_csv(str(f))
    assert [s.id for s in added] == ["XeF4(s)", "Xe(g)"]
    p = Reaction.from_string("Xe + 2 F2 -> XeF4", db).properties()
    assert p.dH / 1000 == pytest.approx(-267)
    assert "XeF4(s)" not in DB                              # default database untouched


def test_user_data_file_from_excel_with_bom(tmp_path):
    f = tmp_path / "excel.csv"
    f.write_bytes(b"\xef\xbb\xbfformula,phase,dHf_kJ,S_J,Cp_J\r\nPbO,s,-218.99,66.5,45.77\r\n")
    db = DB.copy()
    assert [s.id for s in db.load_csv(str(f))] == ["PbO(s)"]


def test_user_data_file_missing_column(tmp_path):
    f = tmp_path / "bad.csv"
    f.write_text("formula,phase,S_J\nPbO,s,66.5\n", encoding="utf-8")
    with pytest.raises(ValueError, match="dHf_kJ"):
        DB.copy().load_csv(str(f))
