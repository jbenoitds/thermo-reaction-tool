"""Worked problems from CHM1402 (Chimie physique 2, Université de Montréal), Série d'exercices #5
« Équilibres chimiques » (inspired by Atkins' Physical Chemistry), checked against the answer key.
"""
import math

import pytest

from thermoreact import Reaction, k_from_extent, solve_equilibrium
from thermoreact.relations import (K_dissociation, K_from_dG, K_from_E0, crossover_T, dG_from_E,
                                   dG_from_K, vant_hoff_dH, vant_hoff_K2, vant_hoff_T2)


def test_ex3_N2O4_degree_of_dissociation():
    """α = 0.201 at 298 K, P = 1.00 bar  ->  K = 0.169 (answer key)."""
    rxn = Reaction.from_string("N2O4 -> 2 NO2")
    K = k_from_extent(rxn, {"N2O4": 1.0}, 0.201, P_bar=1.0)
    assert K == pytest.approx(0.169, rel=0.005)
    assert K == pytest.approx(K_dissociation(0.201, 1.0), rel=1e-12)
    # and the reverse problem: K -> α
    eq = solve_equilibrium(rxn, {"N2O4": 1.0}, K=K, P_bar=1.0)
    assert eq.extent == pytest.approx(0.201, rel=1e-9)


def test_ex3_compared_with_tabulated_data():
    """Tabulated data give K(298 K) ≈ 0.15 for N2O4 = 2 NO2, i.e. α ≈ 0.19 at 1 bar."""
    rxn = Reaction.from_string("N2O4 -> 2 NO2")
    eq = solve_equilibrium(rxn, {"N2O4": 1.0}, T=298.15, P_bar=1.0)
    assert eq.K == pytest.approx(0.148, abs=0.01)
    assert eq.extent == pytest.approx(0.189, abs=0.01)


def test_ex5_Kp_and_dG_from_composition():
    """A + B = C + 2 D: n0 = 2, 1, 0, 3 mol; 0.79 mol C at equilibrium, P = 1 bar, 298 K.
    Answer key: Kp = 9.58, ΔrG° = -5.60 kJ/mol."""
    rxn = Reaction.from_string("A + B -> C + 2 D", abstract=True)
    K = k_from_extent(rxn, {"A": 2.0, "B": 1.0, "D": 3.0}, 0.79, P_bar=1.0)
    assert K == pytest.approx(9.58, rel=0.005)
    assert dG_from_K(K, 298) / 1000 == pytest.approx(-5.60, abs=0.01)


def test_ex8_Br2_dissociation_and_vant_hoff():
    """α = 0.24 at 1600 K and 1.67 bar; ΔrH° = +112 kJ/mol.  Key: K = 0.414, then 2.23 at 2000 K."""
    K1600 = K_dissociation(0.24, 1.67)
    # The key's 0.414 corresponds to P/P° = 1.67 x 1.01325 (pressure taken as 1.67 atm with P° = 1 bar);
    # with P = 1.67 bar exactly, K = 0.408.
    assert K1600 == pytest.approx(0.408, abs=0.001)
    assert K_dissociation(0.24, 1.67 * 1.01325) == pytest.approx(0.414, abs=0.001)
    assert vant_hoff_K2(0.414, 1600, 2000, 112e3) == pytest.approx(2.23, abs=0.005)


def test_ex9_temperature_where_K_exceeds_1():
    """ΔrG°(1120 K) = +22 kJ/mol, ΔrH° = +125 kJ/mol  ->  K > 1 above 1359 K."""
    K1 = K_from_dG(22e3, 1120)
    assert vant_hoff_T2(K1, 1120, 1.0, 125e3) == pytest.approx(1359, abs=0.5)


def test_ex10_dH_from_K_halving_and_T_for_tripling():
    """K halves from 310 K to 325 K  ->  ΔrH° = -38.7 kJ/mol;  K x 3 at T = 289 K."""
    dH = vant_hoff_dH(1.0, 310, 0.5, 325)
    assert dH / 1000 == pytest.approx(-38.7, abs=0.05)
    assert vant_hoff_T2(1.0, 310, 3.0, dH) == pytest.approx(289, abs=0.5)


def test_ex11_CuSO4_pentahydrate_dehydration_temperature():
    """CuSO4·5H2O(s) -> CuSO4(s) + 5 H2O(g) is significant when K = 1: T = 398 K (key)."""
    rxn = Reaction.from_string("CuSO4·5H2O -> CuSO4 + 5 H2O(g)")
    p = rxn.properties(298.15)
    assert crossover_T(p.dH, p.dS) == pytest.approx(398, abs=1.0)        # textbook: ΔH°/ΔS°
    assert rxn.crossover_temperature("none")[0] == pytest.approx(398, abs=1.0)
    assert rxn.crossover_temperature("full")[0] == pytest.approx(398, abs=1.5)  # with Cp(T): still ~398 K


def test_ex12_cell_gibbs_energy():
    """E = +0.4658 V, n = 2  ->  ΔrG = -89.89 kJ/mol."""
    assert dG_from_E(0.4658, 2) / 1000 == pytest.approx(-89.89, abs=0.01)


def test_ex13_K_from_standard_potentials():
    """(a) Sn + Cu2+ -> Cu + Sn2+: E° = 0.34 - (-0.14) V -> K = 1.72e16;
    (b) Cu2+ + Cu -> 2 Cu+: E° = 0.16 - 0.52 V -> K = 8.15e-7 (key, computed at 298 K)."""
    assert K_from_E0(0.34 - (-0.14), 2, 298) == pytest.approx(1.72e16, rel=0.02)
    assert K_from_E0(0.16 - 0.52, 1, 298) == pytest.approx(8.15e-7, rel=0.02)


def test_ex13a_from_tabulated_formation_data():
    """The same equilibrium from ΔfH°/S° tables (no potentials needed)."""
    p = Reaction.from_string("Sn + Cu2+ -> Cu + Sn2+").properties(298.15)
    assert math.log10(p.K) == pytest.approx(math.log10(1.72e16), abs=0.5)


def test_ex14_AgI_solubility_product():
    """Ag|AgI|AgI(aq)|Ag with E° = -0.9509 V for AgI -> Ag+ + I-: Ksp = 8.3e-17, s = 9.1e-9 mol/L."""
    Ksp = K_from_E0(-0.9509, 1, 298.15)
    assert Ksp == pytest.approx(8.3e-17, rel=0.03)
    assert math.sqrt(Ksp) == pytest.approx(9.1e-9, rel=0.02)
