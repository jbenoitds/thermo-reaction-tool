"""Textbook thermodynamic relations used in equilibrium problems.

All energies in J/mol, temperatures in K, potentials in V.
"""
from __future__ import annotations

import math

from .constants import F, P_STD, R


def K_from_dG(dG: float, T: float) -> float:
    """K = exp(−ΔrG°/RT)."""
    return math.exp(-dG / (R * T))


def dG_from_K(K: float, T: float) -> float:
    """ΔrG° = −RT ln K."""
    return -R * T * math.log(K)


def vant_hoff_K2(K1: float, T1: float, T2: float, dH: float) -> float:
    """Integrated van 't Hoff equation with constant ΔrH°: ln(K2/K1) = −ΔrH°/R (1/T2 − 1/T1)."""
    return K1 * math.exp(-dH / R * (1.0 / T2 - 1.0 / T1))


def vant_hoff_dH(K1: float, T1: float, K2: float, T2: float) -> float:
    """ΔrH° from equilibrium constants at two temperatures."""
    return -R * math.log(K2 / K1) / (1.0 / T2 - 1.0 / T1)


def vant_hoff_T2(K1: float, T1: float, K2: float, dH: float) -> float:
    """Temperature at which K reaches K2, given K1 at T1 and constant ΔrH°."""
    inv = 1.0 / T1 - R * math.log(K2 / K1) / dH
    if inv <= 0:
        raise ValueError("No positive temperature reaches K = %g with these values" % K2)
    return 1.0 / inv


def crossover_T(dH: float, dS: float) -> float:
    """T at which ΔrG° = ΔrH° − TΔrS° = 0 (constant ΔrH°, ΔrS°)."""
    if dS == 0:
        raise ValueError("ΔrS° = 0: ΔrG° has the sign of ΔrH° at all T")
    T = dH / dS
    if T <= 0:
        raise ValueError("ΔrH° and ΔrS° have opposite signs: no crossover temperature")
    return T


def dG_from_E(E: float, n: int) -> float:
    """ΔrG = −nFE."""
    return -n * F * E


def K_from_E0(E0: float, n: int, T: float = 298.15) -> float:
    """K = exp(nFE°/RT)."""
    return math.exp(n * F * E0 / (R * T))


def E0_from_K(K: float, n: int, T: float = 298.15) -> float:
    return R * T * math.log(K) / (n * F)


def Kc_from_K(K: float, delta_n_gas: float, T: float, c0: float = 1.0) -> float:
    """Kc (c° in mol/L) from the pressure-based K (P° = 1 bar) for ideal gases."""
    return K * (P_STD / (c0 * 1000.0 * R * T)) ** delta_n_gas


def K_dissociation(alpha: float, P_bar: float, nu: int = 2) -> float:
    """K for A(g) ⇌ ν B(g) from the degree of dissociation α at total pressure P (pure A initially).

    K = (ν α)^ν / ((1 − α)(1 + (ν−1)α)^(ν−1)) · (P/P°)^(ν−1)
    """
    if not 0 <= alpha < 1:
        raise ValueError("alpha must be in [0, 1)")
    ntot = 1 + (nu - 1) * alpha
    return (nu * alpha) ** nu / ((1 - alpha) * ntot ** (nu - 1)) * (P_bar * 1e5 / P_STD) ** (nu - 1)
