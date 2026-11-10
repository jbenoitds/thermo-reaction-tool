"""Chemical reactions: parsing, balancing and standard reaction properties at any temperature.

For a reaction Σ ν_i X_i = 0 (ν < 0 for reactants) the tool evaluates

* ΔrH°(T) = Σ ν_i [ΔfH°_i(298.15) + ∫ Cp,i dT]          (Kirchhoff)
* ΔrS°(T) = Σ ν_i [S°_i(298.15) + ∫ Cp,i/T dT]
* ΔrG°(T) = ΔrH°(T) − T ΔrS°(T)
* K(T)    = exp(−ΔrG°(T) / RT)
"""
from __future__ import annotations

import math
import re
from fractions import Fraction
from functools import reduce
from typing import Dict, List, Optional, Sequence, Tuple

from .constants import P_STD, PHASE_ALIASES, R, T_REF
from .database import Database, split_phase
from .species import CP_MODES, Species


class ReactionError(ValueError):
    """Problems with a reaction string or its stoichiometry."""


_ARROWS = ["<=>", "<->", "⇌", "⇄", "↔", "-->", "->", "→", "=>"]
_COEF_RE = re.compile(r"^\s*(\d+/\d+|\d*\.\d+|\d+)\s*(?=\S)")


def _split_sides(text: str) -> Tuple[str, str]:
    for arrow in _ARROWS:
        if arrow in text:
            left, right = text.split(arrow, 1)
            return left, right
    if " = " in text:
        left, right = text.split(" = ", 1)
        return left, right
    raise ReactionError("No reaction arrow found in %r (use '->', '=>', '<=>', '⇌' or ' = ')" % text)


def _split_terms(side: str) -> List[str]:
    side = side.strip()
    if not side:
        raise ReactionError("Empty side of reaction")
    if re.search(r"\s\+\s", side):
        terms = re.split(r"\s+\+\s+", side)
    else:
        # compact form "CH4+2O2": split on '+' followed by a coefficient/formula start
        terms = re.split(r"(?<=[A-Za-z0-9\)\]])\+(?=\d*\s*[A-Z\[(]|smi)", side)
    return [t.strip() for t in terms if t.strip()]


def _parse_term(term: str) -> Tuple[Fraction, str]:
    m = _COEF_RE.match(term)
    if m and not term.lower().startswith(("smiles:", "smi:")):
        rest = term[m.end():]
        # "2H2O" -> 2, "H2O"; but keep formulas that start with a digit-free element
        if rest and (rest[0].isupper() or rest[0] in "[(" or rest.lower().startswith(("smiles:", "smi:"))
                     or rest[0].islower()):
            coef = Fraction(m.group(1)) if "." not in m.group(1) else Fraction(m.group(1)).limit_denominator(1000)
            return coef, rest.strip()
    return Fraction(1), term.strip()


def _abstract_species(text: str) -> Species:
    body, ph = split_phase(text)
    phase = PHASE_ALIASES.get((ph or "g").lower(), "g")
    return Species(id="%s(%s)" % (body, phase), formula=body, phase=phase, elements={}, charge=0,
                   dHf=None, S=None, name=body, source="abstract (no data)")


class Reaction:
    """A reaction ``Σ ν_i X_i = 0``; stoichiometric numbers are negative for reactants."""

    def __init__(self, terms: Sequence[Tuple[float, Species]], text: str = "", abstract: bool = False):
        net: Dict[str, List] = {}
        order: List[str] = []
        for nu, sp in terms:
            if sp.id not in net:
                net[sp.id] = [Fraction(nu).limit_denominator(10000), sp]
                order.append(sp.id)
            else:
                net[sp.id][0] += Fraction(nu).limit_denominator(10000)
        self.terms: List[Tuple[Fraction, Species]] = [(net[i][0], net[i][1]) for i in order if net[i][0] != 0]
        if not self.terms:
            raise ReactionError("Reaction has no net change")
        self.text = text
        self.abstract = abstract
        self.notes: List[str] = []

    # ------------------------------------------------------------------ construction
    @classmethod
    def from_string(cls, text: str, db: Optional[Database] = None, balance: bool = False,
                    online: bool = False, abstract: bool = False, check: bool = True) -> "Reaction":
        """Parse e.g. ``"CH4 + 2 O2 -> CO2 + 2 H2O(l)"``.

        ``abstract=True`` treats every species as a generic substance without data
        (gas unless a phase is given), useful for composition-only problems such as
        ``A + B -> C + 2 D``.
        """
        left, right = _split_sides(text)
        db = db or Database.default()
        terms: List[Tuple[Fraction, Species]] = []
        for sign, side in ((-1, left), (1, right)):
            for term in _split_terms(side):
                coef, name = _parse_term(term)
                sp = _abstract_species(name) if abstract else db.get(name, online=online)
                terms.append((sign * coef, sp))
        rxn = cls(terms, text=text, abstract=abstract)
        if abstract:
            return rxn
        if balance:
            rxn = rxn.balanced()
            rxn.notes.append("Balanced automatically: %s" % rxn.equation())
        elif check:
            problems = rxn.imbalance()
            if problems:
                raise ReactionError("Reaction is not balanced (products minus reactants: %s). "
                                    "Fix the coefficients or use --balance."
                                    % ", ".join("%s %+g" % (k, float(v)) for k, v in problems.items()))
        return rxn

    # ------------------------------------------------------------------ stoichiometry
    @property
    def species(self) -> List[Species]:
        return [sp for _, sp in self.terms]

    @property
    def reactants(self) -> List[Tuple[Fraction, Species]]:
        return [(-nu, sp) for nu, sp in self.terms if nu < 0]

    @property
    def products(self) -> List[Tuple[Fraction, Species]]:
        return [(nu, sp) for nu, sp in self.terms if nu > 0]

    def nu(self, sid: str) -> Fraction:
        for nu, sp in self.terms:
            if sp.id == sid:
                return nu
        return Fraction(0)

    def imbalance(self) -> Dict[str, Fraction]:
        """Net change of each element (and 'charge') across the reaction; empty if balanced."""
        tot: Dict[str, Fraction] = {}
        for nu, sp in self.terms:
            for el, n in sp.elements.items():
                tot[el] = tot.get(el, Fraction(0)) + nu * n
            if sp.charge:
                tot["charge"] = tot.get("charge", Fraction(0)) + nu * sp.charge
        return {k: v for k, v in tot.items() if v != 0}

    def balanced(self) -> "Reaction":
        """Return the reaction with the smallest whole-number coefficients that balance it."""
        sps = self.species
        signs = [1 if nu > 0 else -1 for nu, _ in self.terms]
        keys = sorted({el for sp in sps for el in sp.elements} | ({"charge"} if any(sp.charge for sp in sps) else set()))
        mat = [[Fraction(sp.charge if k == "charge" else sp.elements.get(k, 0)) * s for sp, s in zip(sps, signs)]
               for k in keys]
        basis = _nullspace(mat, len(sps))
        if len(basis) == 0:
            raise ReactionError("These species cannot be balanced (check formulas and phases)")
        if len(basis) > 1:
            raise ReactionError("The coefficients are not unique (%d independent reactions); "
                                "please give them explicitly" % len(basis))
        vec = basis[0]
        if all(v <= 0 for v in vec):
            vec = [-v for v in vec]
        if any(v <= 0 for v in vec):
            raise ReactionError("No balance with all species on the given sides; "
                                "a species may be on the wrong side")
        lcm = reduce(lambda a, b: a * b // math.gcd(a, b), (v.denominator for v in vec), 1)
        ints = [int(v * lcm) for v in vec]
        g = reduce(math.gcd, ints)
        ints = [i // g for i in ints]
        new = Reaction([(s * c, sp) for s, c, sp in zip(signs, ints, sps)], text=self.text)
        new.notes = list(self.notes)
        return new

    def scaled(self, factor: float) -> "Reaction":
        r = Reaction([(nu * Fraction(factor).limit_denominator(10000), sp) for nu, sp in self.terms], self.text)
        r.notes = list(self.notes)
        return r

    @property
    def delta_n_gas(self) -> Fraction:
        return sum((nu for nu, sp in self.terms if sp.is_gas), Fraction(0))

    # ------------------------------------------------------------------ display
    @staticmethod
    def _coef(c: Fraction) -> str:
        if c == 1:
            return ""
        if c.denominator == 1:
            return "%d " % c.numerator
        return "%d/%d " % (c.numerator, c.denominator)

    def equation(self, arrow: str = "→") -> str:
        lhs = " + ".join(self._coef(c) + sp.id for c, sp in self.reactants)
        rhs = " + ".join(self._coef(c) + sp.id for c, sp in self.products)
        return "%s %s %s" % (lhs, arrow, rhs)

    __str__ = equation

    def __repr__(self) -> str:
        return "Reaction(%r)" % self.equation("->")

    # ------------------------------------------------------------------ thermodynamics
    def _check_data(self) -> None:
        missing = [sp.id for sp in self.species if not sp.has_data]
        if missing:
            raise ReactionError("No thermochemical data for: %s" % ", ".join(missing))

    def properties(self, T: float = T_REF, cp_mode: str = "full") -> "ReactionProperties":
        """Standard reaction properties at temperature *T* (K)."""
        if cp_mode not in CP_MODES:
            raise ValueError("cp_mode must be one of %s" % (CP_MODES,))
        if T <= 0:
            raise ValueError("Temperature must be positive (kelvin)")
        self._check_data()
        dH = dS = dCp = 0.0
        cp_known = True
        warnings: List[str] = []
        for nu, sp in self.terms:
            n = float(nu)
            dH += n * sp.H(T, cp_mode)
            dS += n * sp.S_at(T, cp_mode)
            c = sp.cp(T, cp_mode)
            if c is None:
                cp_known = False
            else:
                dCp += n * c
            warnings.extend(sp.warnings_at(T, cp_mode))
        return ReactionProperties(self, T, dH, dS, dCp if cp_known else None, cp_mode, warnings)

    def table(self, temperatures: Sequence[float], cp_mode: str = "full") -> List["ReactionProperties"]:
        return [self.properties(T, cp_mode) for T in temperatures]

    def crossover_temperature(self, cp_mode: str = "full", t_lo: float = 100.0, t_hi: float = 6000.0,
                              step: float = 5.0) -> List[float]:
        """Temperatures at which ΔrG° = 0 (K = 1) in [t_lo, t_hi], found by scanning + bisection."""
        self._check_data()
        roots: List[float] = []

        def g(T: float) -> float:
            return self.properties(T, cp_mode).dG

        T0, g0 = t_lo, g(t_lo)
        T = t_lo
        while T < t_hi:
            T1 = min(T + step, t_hi)
            g1 = g(T1)
            if g0 == 0:
                roots.append(T0)
            elif g0 * g1 < 0:
                a, b, ga = T0, T1, g0
                for _ in range(60):
                    m = 0.5 * (a + b)
                    gm = g(m)
                    if ga * gm <= 0:
                        b = m
                    else:
                        a, ga = m, gm
                roots.append(0.5 * (a + b))
            T0, g0, T = T1, g1, T1
        return roots


class ReactionProperties:
    """ΔrH°, ΔrS°, ΔrG°, ΔrCp° and K for one reaction at one temperature (SI units)."""

    def __init__(self, reaction: Reaction, T: float, dH: float, dS: float, dCp: Optional[float],
                 cp_mode: str, warnings: List[str]):
        self.reaction = reaction
        self.T = T
        self.dH = dH
        self.dS = dS
        self.dG = dH - T * dS
        self.dCp = dCp
        self.cp_mode = cp_mode
        self.warnings = warnings

    @property
    def lnK(self) -> float:
        return -self.dG / (R * self.T)

    @property
    def log10K(self) -> float:
        return self.lnK / math.log(10)

    @property
    def K(self) -> float:
        """Thermodynamic equilibrium constant (standard state 1 bar / 1 mol kg⁻¹ / pure phase)."""
        x = self.lnK
        if x > 709:
            return float("inf")
        if x < -745:
            return 0.0
        return math.exp(x)

    @property
    def Kc(self) -> Optional[float]:
        """Concentration-based constant (c° = 1 mol/L) for the gas-phase part: K (P°/(c° R T))^Δn."""
        dn = float(self.reaction.delta_n_gas)
        c0 = 1000.0  # mol/m^3
        return math.exp(self.lnK + dn * math.log(P_STD / (c0 * R * self.T)))

    @property
    def Kp_atm(self) -> float:
        """K expressed with pressures in atm (standard state 1 atm instead of 1 bar)."""
        dn = float(self.reaction.delta_n_gas)
        return math.exp(self.lnK + dn * math.log(1e5 / 101325.0))

    @property
    def spontaneous(self) -> bool:
        return self.dG < 0

    def as_dict(self) -> Dict:
        return {
            "T_K": self.T, "dH_kJ_per_mol": self.dH / 1000, "dS_J_per_mol_K": self.dS,
            "dG_kJ_per_mol": self.dG / 1000, "dCp_J_per_mol_K": self.dCp, "lnK": self.lnK,
            "log10K": self.log10K, "K": self.K, "cp_mode": self.cp_mode,
        }


def _nullspace(mat: List[List[Fraction]], ncols: int) -> List[List[Fraction]]:
    """Exact rational null space of a matrix (list of rows) by Gauss-Jordan elimination."""
    rows = [row[:] for row in mat]
    pivots: List[int] = []
    r = 0
    for c in range(ncols):
        piv = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        pv = rows[r][c]
        rows[r] = [v / pv for v in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[r])]
        pivots.append(c)
        r += 1
        if r == len(rows):
            break
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fcol in free:
        vec = [Fraction(0)] * ncols
        vec[fcol] = Fraction(1)
        for i, pc in enumerate(pivots):
            vec[pc] = -rows[i][fcol]
        basis.append(vec)
    return basis
