"""Heat-capacity models and their exact integrals.

Every model provides, for temperatures in kelvin,

* ``cp(T)``                     Cp° in J/(mol K)
* ``int_cp(T1, T2)``            ∫ Cp dT        in J/mol      (enthalpy change)
* ``int_cp_over_t(T1, T2)``     ∫ Cp/T dT      in J/(mol K)  (entropy change)

so that  H°(T) = ΔfH°(298.15) + ∫ Cp dT  and  S°(T) = S°(298.15) + ∫ Cp/T dT
(Kirchhoff's law and its entropy analogue).
"""
from __future__ import annotations

import math
from typing import List, Optional, Sequence, Tuple


class CpModel:
    """Base class; subclasses implement the three methods on a single smooth segment."""

    tmin: float = 0.0
    tmax: float = float("inf")
    label: str = "Cp model"

    def cp(self, T: float) -> float:  # pragma: no cover - interface
        raise NotImplementedError

    def int_cp(self, T1: float, T2: float) -> float:  # pragma: no cover - interface
        raise NotImplementedError

    def int_cp_over_t(self, T1: float, T2: float) -> float:  # pragma: no cover - interface
        raise NotImplementedError


class ConstantCp(CpModel):
    label = "constant Cp"

    def __init__(self, cp: float, tmin: float = 0.0, tmax: float = float("inf")):
        self.value = float(cp)
        self.tmin, self.tmax = tmin, tmax

    def cp(self, T: float) -> float:
        return self.value

    def int_cp(self, T1: float, T2: float) -> float:
        return self.value * (T2 - T1)

    def int_cp_over_t(self, T1: float, T2: float) -> float:
        return self.value * math.log(T2 / T1)


class Shomate(CpModel):
    """NIST Shomate form: Cp = A + B t + C t² + D t³ + E/t², with t = T/1000."""

    label = "Shomate"

    def __init__(self, A: float, B: float, C: float, D: float, E: float,
                 tmin: float = 298.0, tmax: float = 6000.0, **extra):
        self.A, self.B, self.C, self.D, self.E = A, B, C, D, E
        self.tmin, self.tmax = float(tmin), float(tmax)
        self.extra = extra  # F, G, H, fit_tmin ... (kept for reference)

    def cp(self, T: float) -> float:
        t = T / 1000.0
        return self.A + self.B * t + self.C * t * t + self.D * t ** 3 + self.E / (t * t)

    def _h(self, T: float) -> float:  # kJ/mol, without F - H
        t = T / 1000.0
        return self.A * t + self.B * t * t / 2 + self.C * t ** 3 / 3 + self.D * t ** 4 / 4 - self.E / t

    def _s(self, T: float) -> float:  # J/(mol K), without G
        t = T / 1000.0
        return (self.A * math.log(t) + self.B * t + self.C * t * t / 2 + self.D * t ** 3 / 3
                - self.E / (2 * t * t))

    def int_cp(self, T1: float, T2: float) -> float:
        return (self._h(T2) - self._h(T1)) * 1000.0

    def int_cp_over_t(self, T1: float, T2: float) -> float:
        return self._s(T2) - self._s(T1)


class Polynomial(CpModel):
    """Cp = a + bT + cT² + dT³ (+ e/T²); covers Joback and textbook a + bT + c/T² forms."""

    label = "polynomial"

    def __init__(self, a: float, b: float = 0.0, c: float = 0.0, d: float = 0.0, e: float = 0.0,
                 tmin: float = 0.0, tmax: float = float("inf")):
        self.a, self.b, self.c, self.d, self.e = a, b, c, d, e
        self.tmin, self.tmax = tmin, tmax

    def cp(self, T: float) -> float:
        return self.a + self.b * T + self.c * T ** 2 + self.d * T ** 3 + self.e / T ** 2

    def int_cp(self, T1: float, T2: float) -> float:
        def F(T: float) -> float:
            return self.a * T + self.b * T ** 2 / 2 + self.c * T ** 3 / 3 + self.d * T ** 4 / 4 - self.e / T
        return F(T2) - F(T1)

    def int_cp_over_t(self, T1: float, T2: float) -> float:
        def F(T: float) -> float:
            return (self.a * math.log(T) + self.b * T + self.c * T ** 2 / 2 + self.d * T ** 3 / 3
                    - self.e / (2 * T ** 2))
        return F(T2) - F(T1)


class PiecewiseCp:
    """A sequence of Cp segments covering a temperature range.

    Below the first segment / above the last one, the nearest segment is
    extrapolated and the temperature is reported by :meth:`out_of_range`.
    """

    def __init__(self, segments: Sequence[CpModel], label: str = ""):
        if not segments:
            raise ValueError("PiecewiseCp needs at least one segment")
        self.segments: List[CpModel] = sorted(segments, key=lambda s: s.tmin)
        self.label = label or self.segments[0].label

    @property
    def tmin(self) -> float:
        return self.segments[0].tmin

    @property
    def tmax(self) -> float:
        return self.segments[-1].tmax

    def out_of_range(self, T: float) -> bool:
        return T < self.tmin - 1e-6 or T > self.tmax + 1e-6

    def _segment(self, T: float) -> CpModel:
        for seg in self.segments:
            if T <= seg.tmax:
                return seg
        return self.segments[-1]

    def cp(self, T: float) -> float:
        return self._segment(T).cp(T)

    def _pieces(self, T1: float, T2: float) -> List[Tuple[CpModel, float, float]]:
        lo, hi = (T1, T2) if T1 <= T2 else (T2, T1)
        cuts = [lo] + [s.tmax for s in self.segments[:-1] if lo < s.tmax < hi] + [hi]
        return [(self._segment(0.5 * (a + b)), a, b) for a, b in zip(cuts, cuts[1:])]

    def int_cp(self, T1: float, T2: float) -> float:
        if T1 == T2:
            return 0.0
        total = sum(seg.int_cp(a, b) for seg, a, b in self._pieces(T1, T2))
        return total if T2 >= T1 else -total

    def int_cp_over_t(self, T1: float, T2: float) -> float:
        if T1 == T2:
            return 0.0
        total = sum(seg.int_cp_over_t(a, b) for seg, a, b in self._pieces(T1, T2))
        return total if T2 >= T1 else -total


def make_cp(cp298: Optional[float], segments: Optional[List[dict]] = None,
            source: str = "") -> Optional[PiecewiseCp]:
    """Build the best available Cp model: Shomate segments if given, else constant Cp(298)."""
    if segments:
        segs = [Shomate(**s) for s in segments]
        return PiecewiseCp(segs, label="Shomate")
    if cp298 is not None:
        return PiecewiseCp([ConstantCp(cp298)], label="constant Cp(298 K)")
    return None
