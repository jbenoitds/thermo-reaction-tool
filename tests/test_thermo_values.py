"""Results compared with independent reference values (NIST-JANAF tables, textbooks)."""
import pytest

from thermoreact import Database, Reaction

DB = Database.default()


def props(eq, T=298.15, mode="full", **kw):
    return Reaction.from_string(eq, **kw).properties(T, mode)


@pytest.mark.parametrize("species, dGf_kJ", [
    # standard Gibbs energies of formation at 298.15 K (NIST-JANAF / CODATA / NBS)
    ("CO2(g)", -394.36), ("H2O(l)", -237.13), ("H2O(g)", -228.58), ("NH3(g)", -16.37),
    ("CH4(g)", -50.77), ("CO(g)", -137.16), ("NO(g)", 86.60), ("NO2(g)", 51.26),
    ("C2H5OH(l)", -174.78), ("CaCO3(s)", -1128.8), ("Cl-(aq)", -131.23), ("Na+(aq)", -261.91),
    ("OH-(aq)", -157.24), ("Cu2+(aq)", 65.49), ("SO4^2-(aq)", -744.53),
])
def test_gibbs_of_formation_298(species, dGf_kJ):
    sp = DB.get(species)
    assert DB.formation(sp)["dGf"] / 1000 == pytest.approx(dGf_kJ, abs=0.35)


@pytest.mark.parametrize("species, T, dGf_kJ", [
    # NIST-JANAF ΔfG° at high temperature
    ("NH3(g)", 500, 4.800), ("NH3(g)", 700, 27.190), ("NH3(g)", 1000, 61.910),
    ("H2O(g)", 1000, -192.590), ("H2O(g)", 1500, -164.376),
    ("NO(g)", 1000, 77.772), ("HCl(g)", 1000, -100.799),
])
def test_gibbs_of_formation_high_T(species, T, dGf_kJ):
    sp = DB.get(species)
    assert DB.formation(sp, T)["dGf"] / 1000 == pytest.approx(dGf_kJ, abs=0.15)


def test_methane_combustion():
    p = props("CH4 + 2 O2 -> CO2 + 2 H2O(l)")
    assert p.dH / 1000 == pytest.approx(-890.3, abs=0.2)
    assert p.dG / 1000 == pytest.approx(-817.9, abs=0.3)


def test_haber_bosch():
    p = props("N2 + 3 H2 -> 2 NH3")
    assert p.dH / 1000 == pytest.approx(-91.8, abs=0.1)
    assert p.dS == pytest.approx(-198.1, abs=0.2)
    assert p.K == pytest.approx(5.4e5, rel=0.03)
    # K(700 K) from JANAF: ΔfG°(NH3, 700 K) = 27.190 kJ/mol
    assert props("N2 + 3 H2 -> 2 NH3", 700).dG / 1000 == pytest.approx(2 * 27.190, abs=0.1)


def test_water_gas_shift_1000K():
    # JANAF at 1000 K: ΔrG° = -395.886 + 200.275 + 192.590 = -3.021 kJ/mol, K ≈ 1.44
    p = props("CO + H2O(g) -> CO2 + H2", 1000)
    assert p.dG / 1000 == pytest.approx(-3.02, abs=0.1)
    assert p.K == pytest.approx(1.44, rel=0.02)


def test_boudouard_reaction():
    # C + CO2 = 2 CO; JANAF ΔrG°(1000 K) = -4.66 kJ/mol (graphite uses the Atkins a+bT+c/T² Cp)
    p = props("C(graphite) + CO2 -> 2 CO", 1000)
    assert p.dG / 1000 == pytest.approx(-4.66, abs=1.0)


def test_silver_chloride_solubility_product():
    p = props("AgCl -> Ag+ + Cl-")
    assert p.K == pytest.approx(1.77e-10, rel=0.05)


def test_cp_modes_agree_at_298_and_differ_elsewhere():
    rxn = Reaction.from_string("N2 + 3 H2 -> 2 NH3")
    a, b, c = (rxn.properties(298.15, m) for m in ("full", "constant", "none"))
    assert a.dG == pytest.approx(b.dG) == pytest.approx(c.dG)
    hi = [rxn.properties(800, m).dH for m in ("full", "constant", "none")]
    assert hi[2] == pytest.approx(a.dH)                      # 'none' keeps ΔH(298)
    assert hi[1] == pytest.approx(b.dH + b.dCp * (800 - 298.15))  # Kirchhoff with constant ΔrCp
    assert abs(hi[0] - hi[1]) > 1000                          # T-dependent Cp matters


def test_kirchhoff_consistency():
    """d(ΔrH)/dT = ΔrCp and d(ΔrG/T)/d(1/T) = ΔrH (Gibbs-Helmholtz), numerically."""
    rxn = Reaction.from_string("CO + H2O(g) -> CO2 + H2")
    T, h = 900.0, 0.01
    p, p1, p2 = rxn.properties(T), rxn.properties(T - h), rxn.properties(T + h)
    assert (p2.dH - p1.dH) / (2 * h) == pytest.approx(p.dCp, rel=1e-4)
    gh = ((p2.dG / p2.T) - (p1.dG / p1.T)) / (1 / p2.T - 1 / p1.T)
    assert gh == pytest.approx(p.dH, rel=1e-5)


def test_kc_and_kp_atm():
    p = props("N2O4 -> 2 NO2")
    assert p.Kc == pytest.approx(p.K * (1e5 / (1000 * 8.314462618 * 298.15)), rel=1e-9)
    assert p.Kp_atm == pytest.approx(p.K * (1e5 / 101325), rel=1e-9)


def test_crossover_temperature():
    rxn = Reaction.from_string("CaCO3 -> CaO + CO2")
    roots = rxn.crossover_temperature()
    assert len(roots) == 1 and 1050 < roots[0] < 1120
    assert rxn.properties(roots[0]).dG == pytest.approx(0, abs=1.0)


def test_warnings_outside_validity():
    p = props("H2O(l) -> H2O(g)", 450)
    assert any("373" in w for w in p.warnings)
    p = props("C2H6 + 7/2 O2 -> 2 CO2 + 3 H2O(g)", 1200)
    assert any("constant Cp" in w for w in p.warnings)


def test_vaporisation_of_water_normal_boiling_point():
    # ΔvapG° = 0 near 373 K (1 bar standard state -> ~372.8 K)
    rxn = Reaction.from_string("H2O(l) -> H2O(g)")
    roots = rxn.crossover_temperature(t_lo=300, t_hi=450, step=2)
    assert roots[0] == pytest.approx(372.8, abs=1.0)
