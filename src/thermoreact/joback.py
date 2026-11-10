"""Joback group-contribution estimates for ideal-gas thermochemistry.

Reference: K. G. Joback and R. C. Reid, "Estimation of pure-component properties
from group-contributions", Chem. Eng. Commun. 57 (1987) 233-243.  Also tabulated in
Poling, Prausnitz & O'Connell, *The Properties of Gases and Liquids*, 5th ed. (2001).

For a molecule split into groups *k* with counts *N_k*:

* ΔfH°(298 K, gas) = 68.29 + Σ N_k ΔH_k                     (kJ/mol)
* ΔfG°(298 K, gas) = 53.88 + Σ N_k ΔG_k                     (kJ/mol)
* Cp°(T) = (Σa − 37.93) + (Σb + 0.210) T + (Σc − 3.91e-4) T² + (Σd + 2.06e-7) T³   (J/(mol K))

The standard entropy follows from ΔfG = ΔfH − T ΔfS and the entropies of the elements:
S°(298) = (ΔfH − ΔfG)/298.15 + Σ_elements n S°(element).

Typical errors are a few kJ/mol for ΔfH and ΔfG, larger for small or highly
functionalised molecules; the Cp polynomial is intended for roughly 273-1000 K.
"""
from __future__ import annotations

from collections import Counter
from typing import Dict, List, Optional, Tuple

from .smiles import Molecule, parse_smiles


class JobackError(ValueError):
    """Raised when a molecule cannot be decomposed into Joback groups."""


# group: (dHf kJ/mol, dGf kJ/mol, a, b, c, d) ; None = no published value
# fmt: off
GROUPS: Dict[str, Tuple[Optional[float], ...]] = {
    "-CH3":                  (-76.45, -43.96, 1.95e1, -8.08e-3, 1.53e-4, -9.67e-8),
    "-CH2-":                 (-20.64, 8.42, -9.09e-1, 9.50e-2, -5.44e-5, 1.19e-8),
    ">CH-":                  (29.89, 58.36, -2.30e1, 2.04e-1, -2.65e-4, 1.20e-7),
    ">C<":                   (82.23, 116.02, -6.62e1, 4.27e-1, -6.41e-4, 3.01e-7),
    "=CH2":                  (-9.63, 3.77, 2.36e1, -3.81e-2, 1.72e-4, -1.03e-7),
    "=CH-":                  (37.97, 48.53, -8.00, 1.05e-1, -9.63e-5, 3.56e-8),
    "=C<":                   (83.99, 92.36, -2.81e1, 2.08e-1, -3.06e-4, 1.46e-7),
    "=C=":                   (142.14, 136.70, 2.74e1, -5.57e-2, 1.01e-4, -5.02e-8),
    "≡CH":                   (79.30, 77.71, 2.45e1, -2.71e-2, 1.11e-4, -6.78e-8),
    "≡C-":                   (115.51, 109.82, 7.87, 2.01e-2, -8.33e-6, 1.39e-9),
    "-CH2- (ring)":          (-26.80, -3.68, -6.03, 8.54e-2, -8.00e-6, -1.80e-8),
    ">CH- (ring)":           (8.67, 40.99, -2.05e1, 1.62e-1, -1.60e-4, 6.24e-8),
    ">C< (ring)":            (79.72, 87.88, -9.09e1, 5.57e-1, -9.00e-4, 4.69e-7),
    "=CH- (ring)":           (2.09, 11.30, -2.14, 5.74e-2, -1.64e-6, -1.59e-8),
    "=C< (ring)":            (46.43, 54.05, -8.25, 1.01e-1, -1.42e-4, 6.78e-8),
    "-F":                    (-251.92, -247.19, 2.65e1, -9.13e-2, 1.91e-4, -1.03e-7),
    "-Cl":                   (-71.55, -64.31, 3.33e1, -9.63e-2, 1.87e-4, -9.96e-8),
    "-Br":                   (-29.48, -38.06, 2.86e1, -6.49e-2, 1.36e-4, -7.45e-8),
    "-I":                    (21.06, 5.74, 3.21e1, -6.41e-2, 1.26e-4, -6.87e-8),
    "-OH (alcohol)":         (-208.04, -189.20, 2.57e1, -6.91e-2, 1.77e-4, -9.88e-8),
    "-OH (phenol)":          (-221.65, -197.37, -2.81, 1.11e-1, -1.16e-4, 4.94e-8),
    "-O- (nonring)":         (-132.22, -105.00, 2.55e1, -6.32e-2, 1.11e-4, -5.48e-8),
    "-O- (ring)":            (-138.16, -98.22, 1.22e1, -1.26e-2, 6.03e-5, -3.86e-8),
    ">C=O (nonring)":        (-133.22, -120.50, 6.45, 6.70e-2, -3.57e-5, 2.86e-9),
    ">C=O (ring)":           (-164.50, -126.27, 3.04e1, -8.29e-2, 2.36e-4, -1.31e-7),
    "O=CH- (aldehyde)":      (-162.03, -143.48, 3.09e1, -3.36e-2, 1.60e-4, -9.88e-8),
    "-COOH (acid)":          (-426.72, -387.87, 2.41e1, 4.27e-2, 8.04e-5, -6.87e-8),
    "-COO- (ester)":         (-337.92, -301.95, 2.45e1, 4.02e-2, 4.02e-5, -4.52e-8),
    "=O (other than above)": (-247.61, -250.83, 6.82, 1.96e-2, 1.27e-5, -1.78e-8),
    "-NH2":                  (-22.02, 14.07, 2.69e1, -4.12e-2, 1.64e-4, -9.76e-8),
    ">NH (nonring)":         (53.47, 89.39, -1.21, 7.62e-2, -4.86e-5, 1.05e-8),
    ">NH (ring)":            (31.65, 75.61, 1.18e1, -2.30e-2, 1.07e-4, -6.28e-8),
    ">N- (nonring)":         (123.34, 163.16, -3.11e1, 2.27e-1, -3.20e-4, 1.46e-7),
    "-N= (nonring)":         (23.61, None, None, None, None, None),
    "-N= (ring)":            (55.52, 79.93, 8.83, -3.84e-3, 4.35e-5, -2.60e-8),
    "=NH":                   (93.70, 119.66, 5.69, -4.12e-3, 1.28e-4, -8.88e-8),
    "-CN":                   (88.43, 89.22, 3.65e1, -7.33e-2, 1.84e-4, -1.03e-7),
    "-NO2":                  (-66.57, -16.83, 2.59e1, -3.74e-3, 1.29e-4, -8.88e-8),
    "-SH":                   (-17.33, -22.99, 3.53e1, -7.58e-2, 1.85e-4, -1.03e-7),
    "-S- (nonring)":         (41.87, 33.12, 1.96e1, -5.61e-3, 4.02e-5, -2.76e-8),
    "-S- (ring)":            (39.10, 27.76, 1.67e1, 4.81e-3, 2.77e-5, -2.11e-8),
}
# fmt: on

#: Standard molar entropy of each element in its reference state, per atom (J/(mol K)),
#: consistent with the bundled database (e.g. H = S°(H2, g)/2).
ELEMENT_ENTROPY_PER_ATOM = {
    "C": 5.740, "H": 130.68 / 2, "O": 205.15 / 2, "N": 191.61 / 2, "F": 202.78 / 2,
    "Cl": 223.08 / 2, "Br": 152.23 / 2, "I": 116.14 / 2, "S": 31.80,
}

CP_VALID_RANGE = (273.0, 1000.0)


class JobackResult:
    """Outcome of a Joback estimate (all values for the ideal gas)."""

    def __init__(self, smiles: str, formula: str, groups: Dict[str, int], dHf_kJ: float,
                 dGf_kJ: Optional[float], S_J: Optional[float],
                 cp_coeffs: Optional[Tuple[float, float, float, float]],
                 element_counts: Dict[str, int]):
        self.smiles = smiles
        self.formula = formula
        self.groups = groups
        self.dHf_kJ = dHf_kJ
        self.dGf_kJ = dGf_kJ
        self.S_J = S_J
        self.cp_coeffs = cp_coeffs
        self.element_counts = element_counts

    def cp(self, T: float) -> float:
        if self.cp_coeffs is None:
            raise JobackError("No Joback heat-capacity data for group(s) in %s" % self.smiles)
        a, b, c, d = self.cp_coeffs
        return a + b * T + c * T ** 2 + d * T ** 3

    def as_dict(self) -> Dict:
        return {
            "smiles": self.smiles, "formula": self.formula, "groups": dict(self.groups),
            "dHf_kJ_per_mol": self.dHf_kJ, "dGf_kJ_per_mol": self.dGf_kJ,
            "S_J_per_mol_K": self.S_J, "Cp_coefficients_abcd": self.cp_coeffs,
        }


def _has_ring_double_bond(mol: Molecule, i: int) -> bool:
    for j, o in mol.neighbors[i]:
        if o in (1.5, 2.0) and mol.bond_in_ring(i, j):
            return True
    return False


def assign_groups(mol: Molecule) -> Dict[str, int]:
    """Decompose a molecule into Joback groups; raise :class:`JobackError` if impossible."""
    if any(a.charge for a in mol.atoms) and not _is_nitro_charge_pattern(mol):
        raise JobackError("Charged species are outside the scope of the Joback method")
    atoms = mol.atoms
    groups: Counter = Counter()
    used = set()

    def nb(i: int) -> List[Tuple[int, float]]:
        return mol.neighbors[i]

    # --- multi-atom groups first -----------------------------------------------------
    for i, a in enumerate(atoms):
        if a.element == "N":
            oxy = [(j, o) for j, o in nb(i) if atoms[j].element == "O" and mol.heavy_degree(j) == 1]
            if len(oxy) == 2 and mol.heavy_degree(i) == 3:
                groups["-NO2"] += 1
                used.update([i, oxy[0][0], oxy[1][0]])
    for i, a in enumerate(atoms):
        if a.element != "C" or i in used:
            continue
        carbonyl_o = [j for j, o in nb(i) if o == 2.0 and atoms[j].element == "O"
                      and mol.heavy_degree(j) == 1 and j not in used]
        nitrile_n = [j for j, o in nb(i) if o == 3.0 and atoms[j].element == "N" and mol.heavy_degree(j) == 1]
        if nitrile_n:
            groups["-CN"] += 1
            used.update([i, nitrile_n[0]])
            continue
        if not carbonyl_o:
            continue
        o_dbl = carbonyl_o[0]
        single_o = [j for j, o in nb(i) if o == 1.0 and atoms[j].element == "O" and j not in used]
        hydroxyl = [j for j in single_o if (atoms[j].hcount or 0) == 1 and mol.heavy_degree(j) == 1]
        ether = [j for j in single_o if (atoms[j].hcount or 0) == 0 and mol.heavy_degree(j) == 2]
        if hydroxyl:
            groups["-COOH (acid)"] += 1
            used.update([i, o_dbl, hydroxyl[0]])
        elif ether:
            groups["-COO- (ester)"] += 1
            used.update([i, o_dbl, ether[0]])
        elif (a.hcount or 0) == 1:
            groups["O=CH- (aldehyde)"] += 1
            used.update([i, o_dbl])
        elif (a.hcount or 0) == 0:
            groups[">C=O (ring)" if mol.in_ring(i) else ">C=O (nonring)"] += 1
            used.update([i, o_dbl])
        else:
            raise JobackError("Formaldehyde-type carbon (H2C=O) has no Joback group")

    # --- single-atom groups ------------------------------------------------------------
    for i, a in enumerate(atoms):
        if i in used:
            continue
        el, h, ring = a.element, a.hcount or 0, mol.in_ring(i)
        orders = [o for _, o in nb(i)]
        n_double = sum(1 for o in orders if o == 2.0)
        n_triple = sum(1 for o in orders if o == 3.0)
        aromatic = a.aromatic or any(o == 1.5 for o in orders)
        if el == "C":
            if any(o == 2.0 and atoms[j].element not in ("C", "N") for j, o in nb(i)):
                raise JobackError("C=%s double bonds are not covered by Joback groups" %
                                  next(atoms[j].element for j, o in nb(i) if o == 2.0 and atoms[j].element not in ("C", "N")))
            if n_triple:
                key = {1: "≡CH", 0: "≡C-"}.get(h)
            elif n_double >= 2:
                key = "=C=" if h == 0 else None
            elif n_double == 1 or aromatic:
                key = {2: "=CH2", 1: "=CH- (ring)" if ring else "=CH-", 0: "=C< (ring)" if ring else "=C<"}.get(h)
            else:
                key = {3: "-CH3", 2: "-CH2- (ring)" if ring else "-CH2-",
                       1: ">CH- (ring)" if ring else ">CH-", 0: ">C< (ring)" if ring else ">C<"}.get(h)
            if key is None:
                raise JobackError("Carbon with %d H (e.g. methane) has no Joback group" % h)
        elif el == "O":
            deg = mol.heavy_degree(i)
            if h == 1 and deg == 1:
                j = nb(i)[0][0]
                phenolic = atoms[j].element == "C" and (atoms[j].aromatic or _has_ring_double_bond(mol, j))
                key = "-OH (phenol)" if phenolic else "-OH (alcohol)"
            elif h == 0 and deg == 2:
                key = "-O- (ring)" if ring else "-O- (nonring)"
            elif h == 0 and deg == 1 and n_double == 1:
                key = "=O (other than above)"
            else:
                raise JobackError("Oxygen environment not covered (water or charged O?)")
        elif el == "N":
            if a.aromatic:
                key = ">NH (ring)" if h == 1 else "-N= (ring)"
            elif n_double or aromatic:
                key = "=NH" if h == 1 else ("-N= (ring)" if ring else "-N= (nonring)")
            else:
                key = {2: "-NH2", 1: ">NH (ring)" if ring else ">NH (nonring)", 0: ">N- (nonring)"}.get(h)
                if key is None:
                    raise JobackError("Ammonia-type nitrogen has no Joback group")
        elif el == "S":
            if h == 1 and mol.heavy_degree(i) == 1:
                key = "-SH"
            elif h == 0 and mol.heavy_degree(i) == 2 and n_double == 0:
                key = "-S- (ring)" if (ring or a.aromatic) else "-S- (nonring)"
            else:
                raise JobackError("Sulfur environment not covered (e.g. S=O)")
        elif el in ("F", "Cl", "Br", "I"):
            key = "-" + el
        else:
            raise JobackError("Element %s is not covered by the Joback method" % el)
        groups[key] += 1
    if not groups:
        raise JobackError("No Joback groups found")
    return dict(groups)


def _is_nitro_charge_pattern(mol: Molecule) -> bool:
    """Charge-separated nitro groups ([N+](=O)[O-]) are fine: net charge zero, N+ and O- only."""
    if mol.charge != 0:
        return False
    for a in mol.atoms:
        if a.charge and not ((a.element == "N" and a.charge == 1) or (a.element == "O" and a.charge == -1)):
            return False
    return True


def estimate(smiles_or_mol, element_entropy: Optional[Dict[str, float]] = None) -> JobackResult:
    """Joback estimate of ΔfH°, ΔfG°, S° (298.15 K) and Cp°(T) for an ideal gas.

    >>> r = estimate("CCO")       # ethanol (exp. ΔfH° = -234.8 kJ/mol)
    >>> round(r.dHf_kJ, 2)
    -236.84
    """
    mol = parse_smiles(smiles_or_mol) if isinstance(smiles_or_mol, str) else smiles_or_mol
    groups = assign_groups(mol)
    sums = [0.0] * 6
    missing = [False] * 6
    for g, n in groups.items():
        vals = GROUPS[g]
        for k in range(6):
            if vals[k] is None:
                missing[k] = True
            else:
                sums[k] += n * vals[k]  # type: ignore[operator]
    dHf = 68.29 + sums[0]
    dGf = None if missing[1] else 53.88 + sums[1]
    cp = None if any(missing[2:]) else (sums[2] - 37.93, sums[3] + 0.210, sums[4] - 3.91e-4, sums[5] + 2.06e-7)
    counts = mol.element_counts()
    S = None
    if dGf is not None:
        ent = element_entropy or ELEMENT_ENTROPY_PER_ATOM
        try:
            S = (dHf - dGf) * 1000.0 / 298.15 + sum(n * ent[el] for el, n in counts.items())
        except KeyError as exc:  # pragma: no cover - all Joback elements are listed
            raise JobackError("No reference entropy for element %s" % exc)
    return JobackResult(mol.smiles, mol.formula, groups, dHf, dGf, S, cp, counts)
