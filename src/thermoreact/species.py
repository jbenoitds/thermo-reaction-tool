"""The :class:`Species` object: standard-state properties of one substance in one phase."""
from __future__ import annotations

from typing import Dict, List, Optional

from .constants import PHASE_NAMES, T_REF
from .formula import hill_formula, molar_mass
from .heat_capacity import PiecewiseCp

#: How heat capacity is treated when moving away from 298.15 K
CP_MODES = ("full", "constant", "none")
CONSTANT_CP_WARN_K = 200.0  # warn if a constant-Cp species is used further than this from 298.15 K


class Species:
    """A chemical species in a definite phase.

    Energies are stored in SI units: ``dHf`` in J/mol, ``S`` and ``cp298`` in J/(mol K).
    """

    def __init__(self, id: str, formula: str, phase: str, elements: Dict[str, int], charge: int,
                 dHf: Optional[float], S: Optional[float], cp298: Optional[float] = None,
                 cp_model: Optional[PiecewiseCp] = None, name: str = "", aliases: Optional[List[str]] = None,
                 cas: str = "", smiles: str = "", std: bool = True, source: str = "",
                 tmin: Optional[float] = None, tmax: Optional[float] = None, estimated: bool = False,
                 notes: Optional[List[str]] = None):
        self.id = id
        self.formula = formula
        self.phase = phase
        self.elements = dict(elements)
        self.charge = charge
        self.dHf = dHf
        self.S = S
        self.cp298 = cp298
        self.cp_model = cp_model
        self.name = name or formula
        self.aliases = list(aliases or [])
        self.cas = cas
        self.smiles = smiles
        self.std = std
        self.source = source
        self.tmin = tmin
        self.tmax = tmax
        self.estimated = estimated
        self.notes = list(notes or [])

    # ------------------------------------------------------------------ basics
    def __repr__(self) -> str:
        return "Species(%r)" % self.id

    @property
    def hill(self) -> str:
        return hill_formula(self.elements, self.charge)

    @property
    def molar_mass(self) -> float:
        return molar_mass(self.elements)

    @property
    def is_gas(self) -> bool:
        return self.phase == "g"

    @property
    def phase_name(self) -> str:
        return PHASE_NAMES.get(self.phase, self.phase)

    @property
    def has_data(self) -> bool:
        return self.dHf is not None and self.S is not None

    # ------------------------------------------------------------------ Cp(T)
    def _cp_value_298(self) -> Optional[float]:
        if self.cp298 is not None:
            return self.cp298
        if self.cp_model is not None:
            return self.cp_model.cp(T_REF)
        return None

    def cp(self, T: float, mode: str = "full") -> Optional[float]:
        """Cp°(T) in J/(mol K); ``None`` if unknown."""
        if mode == "none":
            return 0.0
        if mode == "constant" or self.cp_model is None:
            return self._cp_value_298()
        return self.cp_model.cp(T)

    def int_cp(self, T: float, mode: str = "full") -> float:
        """∫_{298.15}^{T} Cp dT (J/mol); 0 if Cp is unknown or mode is 'none'."""
        if mode == "none" or T == T_REF:
            return 0.0
        if mode == "full" and self.cp_model is not None:
            return self.cp_model.int_cp(T_REF, T)
        c = self._cp_value_298()
        return 0.0 if c is None else c * (T - T_REF)

    def int_cp_over_t(self, T: float, mode: str = "full") -> float:
        """∫_{298.15}^{T} Cp/T dT (J/(mol K))."""
        import math

        if mode == "none" or T == T_REF:
            return 0.0
        if mode == "full" and self.cp_model is not None:
            return self.cp_model.int_cp_over_t(T_REF, T)
        c = self._cp_value_298()
        return 0.0 if c is None else c * math.log(T / T_REF)

    # ------------------------------------------------------------------ H, S, G
    def _require_data(self) -> None:
        if not self.has_data:
            raise ValueError("No thermochemical data for %s" % self.id)

    def H(self, T: float = T_REF, mode: str = "full") -> float:
        """Standard enthalpy on the 'ΔfH°(298 K)' scale: ΔfH°(298.15) + ∫Cp dT  (J/mol)."""
        self._require_data()
        return self.dHf + self.int_cp(T, mode)  # type: ignore[operator]

    def S_at(self, T: float = T_REF, mode: str = "full") -> float:
        """Standard molar entropy at T (J/(mol K))."""
        self._require_data()
        return self.S + self.int_cp_over_t(T, mode)  # type: ignore[operator]

    def G(self, T: float = T_REF, mode: str = "full") -> float:
        """Standard Gibbs energy on the same scale, G = H - T S (J/mol)."""
        return self.H(T, mode) - T * self.S_at(T, mode)

    # ------------------------------------------------------------------ diagnostics
    def warnings_at(self, T: float, mode: str = "full") -> List[str]:
        out: List[str] = []
        if mode == "none" or abs(T - T_REF) < 1e-9:
            pass
        elif self._cp_value_298() is None:
            out.append("%s: no heat-capacity data; treated as Cp = 0 (ΔH, S taken at 298.15 K)" % self.id)
        elif mode == "full" and self.cp_model is not None and not self.cp_model.label.startswith("constant"):
            if self.cp_model.out_of_range(T):
                out.append("%s: T = %.1f K is outside the Cp fit range %.0f-%.0f K (extrapolated)"
                           % (self.id, T, self.cp_model.tmin, self.cp_model.tmax))
        elif abs(T - T_REF) > CONSTANT_CP_WARN_K:
            out.append("%s: constant Cp(298 K) assumed %.0f K away from 298 K (approximation)"
                       % (self.id, abs(T - T_REF)))
        if self.tmax is not None and T > self.tmax + 1e-6:
            out.append("%s: T = %.1f K is above %.1f K, where this phase is normally no longer stable"
                       % (self.id, T, self.tmax))
        if self.tmin is not None and T < self.tmin - 1e-6:
            out.append("%s: T = %.1f K is below %.1f K, where this phase is normally no longer stable"
                       % (self.id, T, self.tmin))
        return out

    def describe_cp(self) -> str:
        if self.cp_model is not None and not self.cp_model.label.startswith("constant"):
            return "%s, %.0f-%.0f K" % (self.cp_model.label, self.cp_model.tmin, self.cp_model.tmax)
        if self.cp298 is not None:
            return "constant Cp(298 K)"
        return "unknown"
