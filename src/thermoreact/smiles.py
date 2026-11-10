"""A small, dependency-free SMILES parser.

It understands the subset of SMILES needed for thermochemistry estimates:

* organic-subset atoms ``B C N O P S F Cl Br I`` and aromatic ``b c n o p s``
* bracket atoms with isotope, hydrogen count and charge: ``[NH4+]``, ``[13CH4]``, ``[nH]``, ``[O-]``
* bonds ``- = # :`` (``/`` and ``\\`` are read as single bonds; stereo is ignored)
* branches ``( )``, ring closures ``1``-``9`` and ``%10``-``%99``, and ``.`` for disconnected parts

Chirality marks (``@``, ``@@``) are accepted and ignored.  Implicit hydrogens are
added following the usual SMILES valence rules.  If RDKit is installed you do not
need it here: this module is intentionally self-contained so the tool runs anywhere.
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Dict, List, Optional, Tuple

from .constants import ATOMIC_MASSES
from .formula import hill_formula


class SmilesError(ValueError):
    """Raised for SMILES strings that cannot be parsed."""


_ORGANIC = {"B": (3,), "C": (4,), "N": (3, 5), "O": (2,), "P": (3, 5), "S": (2, 4, 6),
            "F": (1,), "Cl": (1,), "Br": (1,), "I": (1,)}
_AROMATIC = {"b": "B", "c": "C", "n": "N", "o": "O", "p": "P", "s": "S", "se": "Se", "as": "As"}
_BOND_ORDERS = {"-": 1.0, "=": 2.0, "#": 3.0, "$": 4.0, ":": 1.5, "/": 1.0, "\\": 1.0}
_BRACKET_RE = re.compile(
    r"^\[(?P<iso>\d+)?(?P<el>[A-Z][a-z]?|se|as|[bcnops]|\*)(?P<chiral>@{1,2}(?:TH[12]|AL[12]|SP[123]|TB\d+|OH\d+)?)?"
    r"(?P<h>H\d*)?(?P<charge>[+-]\d*|\+\++|--+)?(?::\d+)?\]$"
)


class Atom:
    __slots__ = ("element", "aromatic", "charge", "hcount", "isotope", "bracket", "index")

    def __init__(self, element: str, aromatic: bool = False, charge: int = 0,
                 hcount: Optional[int] = None, isotope: Optional[int] = None,
                 bracket: bool = False, index: int = 0):
        self.element = element
        self.aromatic = aromatic
        self.charge = charge
        self.hcount = hcount  # explicit H count (bracket atoms) or computed implicit H
        self.isotope = isotope
        self.bracket = bracket
        self.index = index

    def __repr__(self) -> str:  # pragma: no cover - debugging helper
        return "Atom(%s%s, H%s, q=%d)" % (self.element, "*" if self.aromatic else "", self.hcount, self.charge)


class Molecule:
    """Molecular graph: atoms plus bonds ``(i, j, order)``; aromatic bonds have order 1.5."""

    def __init__(self, atoms: List[Atom], bonds: List[Tuple[int, int, float]], smiles: str = ""):
        self.atoms = atoms
        self.bonds = bonds
        self.smiles = smiles
        self.neighbors: List[List[Tuple[int, float]]] = [[] for _ in atoms]
        for i, j, o in bonds:
            self.neighbors[i].append((j, o))
            self.neighbors[j].append((i, o))
        self._ring_atoms: Optional[set] = None
        self._ring_bonds: Optional[set] = None

    # ------------------------------------------------------------------ composition
    def element_counts(self) -> Dict[str, int]:
        c: Counter = Counter()
        for a in self.atoms:
            c[a.element] += 1
            if a.hcount:
                c["H"] += a.hcount
        return dict(c)

    @property
    def charge(self) -> int:
        return sum(a.charge for a in self.atoms)

    @property
    def formula(self) -> str:
        return hill_formula(self.element_counts(), self.charge)

    def heavy_degree(self, i: int) -> int:
        return len(self.neighbors[i])

    # ------------------------------------------------------------------ rings
    def _find_ring_bonds(self) -> None:
        """A bond is in a ring iff it is not a bridge (Tarjan's bridge-finding algorithm)."""
        n = len(self.atoms)
        disc = [-1] * n
        low = [0] * n
        bridges = set()
        counter = [0]
        for root in range(n):
            if disc[root] != -1:
                continue
            # iterative DFS: stack of (node, parent, iterator index)
            stack = [(root, -1, 0)]
            disc[root] = low[root] = counter[0]
            counter[0] += 1
            while stack:
                node, parent, k = stack[-1]
                if k < len(self.neighbors[node]):
                    stack[-1] = (node, parent, k + 1)
                    nb = self.neighbors[node][k][0]
                    if disc[nb] == -1:
                        disc[nb] = low[nb] = counter[0]
                        counter[0] += 1
                        stack.append((nb, node, 0))
                    elif nb != parent:
                        low[node] = min(low[node], disc[nb])
                else:
                    stack.pop()
                    if parent != -1:
                        low[parent] = min(low[parent], low[node])
                        if low[node] > disc[parent]:
                            bridges.add(frozenset((parent, node)))
        ring_bonds = set()
        ring_atoms = set()
        for i, j, _ in self.bonds:
            key = frozenset((i, j))
            if key not in bridges:
                ring_bonds.add(key)
                ring_atoms.update((i, j))
        self._ring_bonds = ring_bonds
        self._ring_atoms = ring_atoms

    def in_ring(self, i: int) -> bool:
        if self._ring_atoms is None:
            self._find_ring_bonds()
        return i in self._ring_atoms  # type: ignore[operator]

    def bond_in_ring(self, i: int, j: int) -> bool:
        if self._ring_bonds is None:
            self._find_ring_bonds()
        return frozenset((i, j)) in self._ring_bonds  # type: ignore[operator]

    # ------------------------------------------------------------------ comparison
    def signature(self) -> Tuple:
        """Order-independent fingerprint (atom environments) used to match database entries."""
        env = []
        for idx, a in enumerate(self.atoms):
            nbrs = tuple(sorted((self.atoms[j].element, o) for j, o in self.neighbors[idx]))
            env.append((a.element, a.aromatic, a.hcount or 0, a.charge, nbrs))
        return tuple(sorted(env))


def _implicit_h(mol_atoms: List[Atom], bond_sum: List[float], aromatic_bonds: List[int]) -> None:
    for i, a in enumerate(mol_atoms):
        if a.bracket:
            if a.hcount is None:
                a.hcount = 0
            continue
        valences = _ORGANIC.get(a.element)
        if valences is None:
            a.hcount = 0
            continue
        used = bond_sum[i]
        if a.aromatic:
            # Aromatic bonds were counted as 1 each.  Aromatic C/B contribute one extra
            # bond's worth to the pi system; aromatic n, o, s, p carry no implicit H
            # (pyrrole-type N must be written [nH], as the SMILES spec requires).
            if a.element in ("C", "B"):
                used += 1
            else:
                a.hcount = 0
                continue
        target = None
        for v in valences:
            if v >= used - 1e-9:
                target = v
                break
        if target is None:
            a.hcount = 0
        else:
            a.hcount = max(0, int(round(target - used)))


def parse_smiles(smiles: str) -> Molecule:
    """Parse a SMILES string into a :class:`Molecule`.

    >>> parse_smiles("CCO").formula
    'C2H6O'
    >>> parse_smiles("c1ccccc1").formula
    'C6H6'
    """
    s = smiles.strip()
    if not s:
        raise SmilesError("Empty SMILES")
    atoms: List[Atom] = []
    bonds: List[Tuple[int, int, float]] = []
    branch_stack: List[int] = []
    rings: Dict[int, Tuple[int, Optional[float]]] = {}
    prev: Optional[int] = None
    pending_bond: Optional[float] = None
    i = 0

    def add_atom(atom: Atom) -> None:
        nonlocal prev, pending_bond
        atom.index = len(atoms)
        atoms.append(atom)
        if prev is not None:
            order = pending_bond
            if order is None:
                order = 1.5 if (atoms[prev].aromatic and atom.aromatic) else 1.0
            bonds.append((prev, atom.index, order))
        prev = atom.index
        pending_bond = None

    while i < len(s):
        ch = s[i]
        if ch == "[":
            j = s.find("]", i)
            if j < 0:
                raise SmilesError("Unclosed bracket atom in %r" % smiles)
            m = _BRACKET_RE.match(s[i: j + 1])
            if not m:
                raise SmilesError("Cannot parse bracket atom %r" % s[i: j + 1])
            el = m.group("el")
            aromatic = el in _AROMATIC
            element = _AROMATIC.get(el, el)
            if element != "*" and element not in ATOMIC_MASSES:
                raise SmilesError("Unknown element %r in %r" % (element, s[i: j + 1]))
            if element == "*":
                raise SmilesError("Wildcard atoms are not supported")
            h = m.group("h")
            hcount = 0 if not h else (int(h[1:]) if len(h) > 1 else 1)
            q = m.group("charge")
            charge = 0
            if q:
                if q in ("+", "-"):
                    charge = 1 if q == "+" else -1
                elif set(q) <= {"+"} or set(q) <= {"-"}:
                    charge = len(q) * (1 if q[0] == "+" else -1)
                else:
                    charge = int(q[1:]) * (1 if q[0] == "+" else -1)
            iso = int(m.group("iso")) if m.group("iso") else None
            add_atom(Atom(element, aromatic, charge, hcount, iso, bracket=True))
            i = j + 1
            continue
        if ch in "BCNOPSFI" or ch in "bcnops":
            if s.startswith("Cl", i) or s.startswith("Br", i):
                add_atom(Atom(s[i: i + 2]))
                i += 2
                continue
            if ch.islower():
                add_atom(Atom(_AROMATIC[ch], aromatic=True))
            else:
                add_atom(Atom(ch))
            i += 1
            continue
        if ch in _BOND_ORDERS:
            pending_bond = _BOND_ORDERS[ch]
            i += 1
            continue
        if ch == "(":
            if prev is None:
                raise SmilesError("Branch before any atom in %r" % smiles)
            branch_stack.append(prev)
            i += 1
            continue
        if ch == ")":
            if not branch_stack:
                raise SmilesError("Unmatched ')' in %r" % smiles)
            prev = branch_stack.pop()
            i += 1
            continue
        if ch.isdigit() or ch == "%":
            if ch == "%":
                digits = s[i + 1: i + 3]
                if len(digits) != 2 or not digits.isdigit():
                    raise SmilesError("Bad ring-closure label in %r" % smiles)
                label = int(digits)
                i += 3
            else:
                label = int(ch)
                i += 1
            if prev is None:
                raise SmilesError("Ring closure before any atom in %r" % smiles)
            if label in rings:
                other, order0 = rings.pop(label)
                order = pending_bond if pending_bond is not None else order0
                if order is None:
                    order = 1.5 if (atoms[other].aromatic and atoms[prev].aromatic) else 1.0
                bonds.append((other, prev, order))
            else:
                rings[label] = (prev, pending_bond)
            pending_bond = None
            continue
        if ch == ".":
            prev = None
            pending_bond = None
            i += 1
            continue
        if ch in "@":
            i += 1
            continue
        raise SmilesError("Unexpected character %r in SMILES %r" % (ch, smiles))

    if rings:
        raise SmilesError("Unclosed ring(s) %s in %r" % (sorted(rings), smiles))
    if branch_stack:
        raise SmilesError("Unclosed branch in %r" % smiles)

    bond_sum = [0.0] * len(atoms)
    aromatic_bonds = [0] * len(atoms)
    for a, b, o in bonds:
        val = 1.0 if o == 1.5 else o
        bond_sum[a] += val
        bond_sum[b] += val
        if o == 1.5:
            aromatic_bonds[a] += 1
            aromatic_bonds[b] += 1
    _implicit_h(atoms, bond_sum, aromatic_bonds)
    return Molecule(atoms, bonds, smiles)


def smiles_formula(smiles: str) -> str:
    """Hill formula of a SMILES string."""
    return parse_smiles(smiles).formula
