"""Equilibrium composition of a single reaction (ideal gases + pure condensed phases).

For a reaction at temperature T and constant total pressure P, the reaction
quotient written with mole fractions of the gas phase is

    Q(ξ) = Π_gas (x_i P / P°)^ν_i ,   x_i = n_i(ξ) / n_gas(ξ),   n_i(ξ) = n_i,0 + ν_i ξ

and pure solids/liquids have activity 1.  The extent of reaction ξ at
equilibrium solves ln Q(ξ) = ln K, which is monotonic in ξ and is found by bisection.
"""
from __future__ import annotations

import math
from typing import Dict, List, Optional

from .constants import P_STD
from .reaction import Reaction, ReactionError


class EquilibriumResult:
    def __init__(self, reaction: Reaction, T: Optional[float], P_bar: float, K: float, extent: float,
                 initial: Dict[str, float], amounts: Dict[str, float], notes: List[str]):
        self.reaction = reaction
        self.T = T
        self.P_bar = P_bar
        self.K = K
        self.extent = extent
        self.initial = initial
        self.amounts = amounts
        self.notes = notes

    @property
    def n_gas(self) -> float:
        return sum(n for sid, n in self.amounts.items() if self._sp(sid).is_gas)

    def _sp(self, sid: str):
        for _, sp in self.reaction.terms:
            if sp.id == sid:
                return sp
        raise KeyError(sid)

    def mole_fractions(self) -> Dict[str, float]:
        ng = self.n_gas
        return {sid: (n / ng if ng > 0 else 0.0) for sid, n in self.amounts.items() if self._sp(sid).is_gas}

    def partial_pressures(self) -> Dict[str, float]:
        return {sid: x * self.P_bar for sid, x in self.mole_fractions().items()}

    def conversion(self) -> Dict[str, float]:
        out = {}
        for nu, sp in self.reaction.terms:
            n0 = self.initial.get(sp.id, 0.0)
            if nu < 0 and n0 > 0:
                out[sp.id] = (n0 - self.amounts[sp.id]) / n0
        return out

    def as_dict(self) -> Dict:
        return {"T_K": self.T, "P_bar": self.P_bar, "K": self.K, "extent_mol": self.extent,
                "amounts_mol": self.amounts, "mole_fractions": self.mole_fractions(),
                "partial_pressures_bar": self.partial_pressures(), "conversion": self.conversion(),
                "notes": self.notes}


def _check_phases(rxn: Reaction) -> None:
    for _, sp in rxn.terms:
        if sp.phase == "aq":
            raise ReactionError("Equilibrium solver handles gases and pure condensed phases only "
                                "(%s is aqueous)" % sp.id)


def _amounts(rxn: Reaction, initial: Dict[str, float], xi: float) -> Dict[str, float]:
    return {sp.id: initial.get(sp.id, 0.0) + float(nu) * xi for nu, sp in rxn.terms}


def ln_q(rxn: Reaction, amounts: Dict[str, float], P_bar: float, inert_gas: float = 0.0) -> float:
    """ln Q for given amounts (mol); gases ideal, condensed phases activity 1."""
    ng = sum(amounts[sp.id] for _, sp in rxn.terms if sp.is_gas) + inert_gas
    total = 0.0
    for nu, sp in rxn.terms:
        if not sp.is_gas:
            continue
        n = amounts[sp.id]
        if n <= 0:
            return -math.inf if nu > 0 else math.inf
        total += float(nu) * math.log(n / ng * P_bar * 1e5 / P_STD)
    return total


def k_from_extent(rxn: Reaction, initial: Dict[str, float], extent: float, P_bar: float = 1.0,
                  inert_gas: float = 0.0) -> float:
    """Equilibrium constant implied by a measured extent of reaction (e.g. degree of dissociation)."""
    _check_phases(rxn)
    initial = _resolve_keys(rxn, initial)
    amounts = _amounts(rxn, initial, extent)
    bad = [k for k, v in amounts.items() if v < -1e-12]
    if bad:
        raise ReactionError("Extent %.6g gives negative amounts for %s" % (extent, ", ".join(bad)))
    return math.exp(ln_q(rxn, amounts, P_bar, inert_gas))


def _resolve_keys(rxn: Reaction, initial: Dict[str, float]) -> Dict[str, float]:
    """Accept initial amounts keyed by id ('N2O4(g)'), formula ('N2O4') or name."""
    out: Dict[str, float] = {}
    for key, val in initial.items():
        match = None
        for _, sp in rxn.terms:
            if key in (sp.id, sp.formula, sp.name) or key in sp.aliases or key.lower() == sp.name.lower():
                match = sp.id
                break
            body = sp.id.rsplit("(", 1)[0]
            if key == body:
                match = sp.id
                break
        if match is None:
            raise ReactionError("Initial amount given for %r, which is not in the reaction" % key)
        out[match] = out.get(match, 0.0) + float(val)
    return out


def solve_equilibrium(rxn: Reaction, initial: Dict[str, float], T: Optional[float] = None,
                      P_bar: float = 1.0, K: Optional[float] = None, cp_mode: str = "full",
                      inert_gas: float = 0.0) -> EquilibriumResult:
    """Equilibrium amounts for one reaction from initial amounts (mol).

    Either *T* (K is then computed from the data) or *K* must be given.
    *inert_gas* adds moles of a non-reacting gas (e.g. N2 diluent in air) to the gas phase.
    """
    _check_phases(rxn)
    notes: List[str] = []
    if K is None:
        if T is None:
            raise ValueError("Give a temperature or an equilibrium constant")
        props = rxn.properties(T, cp_mode)
        lnK = props.lnK
        K = props.K
        notes.extend(props.warnings)
    else:
        if K <= 0:
            raise ValueError("K must be positive")
        lnK = math.log(K)
    initial = _resolve_keys(rxn, initial)
    if any(v < 0 for v in initial.values()):
        raise ValueError("Initial amounts must be non-negative")

    lo, hi = -math.inf, math.inf
    for nu, sp in rxn.terms:
        n0 = initial.get(sp.id, 0.0)
        if nu > 0:
            lo = max(lo, -n0 / float(nu))
        else:
            hi = min(hi, n0 / -float(nu))
    if not (math.isfinite(lo) and math.isfinite(hi)):
        raise ReactionError("Initial amounts leave the extent unbounded")
    if hi - lo <= 0:
        return EquilibriumResult(rxn, T, P_bar, K, 0.0, initial, _amounts(rxn, initial, 0.0),
                                 notes + ["No reaction possible from these initial amounts"])

    def f(xi: float) -> float:
        return ln_q(rxn, _amounts(rxn, initial, xi), P_bar, inert_gas) - lnK

    span = hi - lo
    eps = span * 1e-15
    f_lo, f_hi = f(lo + eps), f(hi - eps)
    gas_terms = [nu for nu, sp in rxn.terms if sp.is_gas]
    if not gas_terms:
        raise ReactionError("No gaseous species: composition is not determined by K "
                            "(the reaction goes to completion in the direction of decreasing G)")
    if abs(f(lo + 0.25 * span) - f(lo + 0.75 * span)) < 1e-12 and math.isfinite(f_lo) and math.isfinite(f_hi):
        # Q independent of ξ (e.g. CaCO3(s) = CaO(s) + CO2(g) with CO2 the only gas, at fixed P)
        xi = hi if f_lo < 0 else lo
        notes.append("Q does not depend on the extent here; the reaction proceeds fully %s"
                     % ("forward" if f_lo < 0 else "backward"))
        return EquilibriumResult(rxn, T, P_bar, K, xi, initial, _amounts(rxn, initial, xi), notes)
    if f_lo >= 0:
        xi = lo
        notes.append("Equilibrium lies at the lower bound of the extent (a condensed species is used up)")
    elif f_hi <= 0:
        xi = hi
        notes.append("Equilibrium lies at the upper bound of the extent (a condensed species is used up)")
    else:
        a, b = lo, hi
        for _ in range(200):
            m = 0.5 * (a + b)
            if f(m) < 0:
                a = m
            else:
                b = m
            if b - a <= 1e-15 * max(1.0, abs(m)):
                break
        xi = 0.5 * (a + b)
    amounts = {k: (0.0 if abs(v) < 1e-14 * max(1.0, span) else v) for k, v in _amounts(rxn, initial, xi).items()}
    return EquilibriumResult(rxn, T, P_bar, K, xi, initial, amounts, notes)
