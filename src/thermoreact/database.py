"""Species database: bundled data files, lookup by formula/name/SMILES, user extensions."""
from __future__ import annotations

import csv
import io
import json
import os
import re
from typing import Dict, Iterable, List, Optional, Tuple

from .constants import ELEMENT_REFERENCE, PHASE_ALIASES, T_REF
from .formula import FormulaError, hill_formula, parse_formula
from .heat_capacity import PiecewiseCp, Polynomial, make_cp
from .joback import CP_VALID_RANGE, JobackError, estimate
from .smiles import SmilesError, parse_smiles
from .species import Species

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

_PHASE_SUFFIX_RE = re.compile(r"^(?P<body>.+?)\s*\((?P<phase>[A-Za-z]+)\)$")
SMILES_PREFIXES = ("smiles:", "smi:")


class SpeciesNotFound(LookupError):
    """No species matches the given text."""


class AmbiguousSpecies(LookupError):
    """Several species match; the message lists them."""


def split_phase(text: str) -> Tuple[str, Optional[str]]:
    """``'H2O(g)' -> ('H2O', 'g')``; ``'C(graphite)' -> ('C', 'graphite')``; no phase -> ``(text, None)``."""
    text = text.strip()
    m = _PHASE_SUFFIX_RE.match(text)
    if m:
        label = m.group("phase")
        # do not mistake a formula group such as Ca(OH) for a phase
        if label.lower() in PHASE_ALIASES or label.islower():
            return m.group("body").strip(), label
    return text, None


def _float(x: str) -> Optional[float]:
    x = (x or "").strip()
    return float(x) if x else None


class Database:
    """A collection of :class:`Species` with flexible lookup.

    >>> db = Database.default()
    >>> db.get("H2O(g)").dHf
    -241830.0
    """

    _default: Optional["Database"] = None

    def __init__(self) -> None:
        self.species: Dict[str, Species] = {}
        self._by_key: Dict[str, List[str]] = {}
        self._by_name: Dict[str, List[str]] = {}
        self._by_hill: Dict[str, List[str]] = {}

    # ------------------------------------------------------------------ loading
    @classmethod
    def default(cls) -> "Database":
        """The bundled database (cached)."""
        if cls._default is None:
            db = cls()
            db.load_csv(os.path.join(DATA_DIR, "species_298.csv"))
            db.load_shomate(os.path.join(DATA_DIR, "shomate.json"))
            cls._default = db
        return cls._default

    def copy(self) -> "Database":
        new = Database()
        for sp in self.species.values():
            new.add(sp)
        return new

    def load_csv(self, path_or_file, override: bool = True) -> List[Species]:
        """Load species from a CSV file with the columns of ``data/species_298.csv``.

        Required columns: ``formula``, ``phase``, ``dHf_kJ``, ``S_J``.  Optional: ``id``,
        ``name``, ``aliases`` (``;``-separated), ``Cp_J``, ``std``, ``cas``, ``smiles``,
        ``tmin``, ``tmax``, ``source``.  Lines starting with ``#`` are comments.
        """
        if hasattr(path_or_file, "read"):
            text = path_or_file.read()
        else:
            # utf-8-sig also accepts files saved by Excel ("CSV UTF-8", with a byte-order mark)
            with open(path_or_file, encoding="utf-8-sig") as fh:
                text = fh.read()
        text = text.lstrip("\ufeff")
        lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
        added = []
        reader = csv.DictReader(io.StringIO("\n".join(lines)))
        missing = [c for c in ("formula", "phase", "dHf_kJ", "S_J")
                   if c not in [(f or "").strip() for f in (reader.fieldnames or [])]]
        if missing:
            raise ValueError("Species CSV is missing required column(s): %s (found: %s)"
                             % (", ".join(missing), ", ".join(reader.fieldnames or []) or "none"))
        for row in reader:
            row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
            formula = row["formula"]
            phase = PHASE_ALIASES.get(row["phase"].lower(), row["phase"].lower())
            elements, charge = parse_formula(formula)
            sid = row.get("id") or "%s(%s)" % (formula, phase)
            std_txt = row.get("std", "1")
            sp = Species(
                id=sid, formula=formula, phase=phase, elements=elements, charge=charge,
                dHf=_float(row["dHf_kJ"]) * 1000.0 if _float(row["dHf_kJ"]) is not None else None,
                S=_float(row["S_J"]), cp298=_float(row.get("Cp_J", "")),
                name=row.get("name", ""),
                aliases=[a.strip() for a in row.get("aliases", "").split(";") if a.strip()],
                cas=row.get("cas", ""), smiles=row.get("smiles", ""),
                std=std_txt not in ("0", "false", "False", "no"),
                source=row.get("source", ""), tmin=_float(row.get("tmin", "")), tmax=_float(row.get("tmax", "")),
            )
            if sp.cp298 is not None:
                sp.cp_model = make_cp(sp.cp298)
            if sid in self.species and not override:
                continue
            self.add(sp)
            added.append(sp)
        return added

    def load_shomate(self, path: str) -> None:
        with open(path, encoding="utf-8") as fh:
            data = json.load(fh)
        for sid, entry in data["species"].items():
            if sid not in self.species:
                raise KeyError("Shomate data for unknown species %s" % sid)
            sp = self.species[sid]
            sp.cp_model = make_cp(sp.cp298, entry["segments"], source=entry.get("source", ""))
            sp.notes.append("Cp(T): " + entry.get("source", "Shomate"))

    # ------------------------------------------------------------------ indexing
    def add(self, sp: Species) -> Species:
        if sp.id in self.species:
            self._unindex(self.species[sp.id])
        self.species[sp.id] = sp
        base, _ = split_phase(sp.id)
        keys = {sp.id, sp.formula, base} | set(a for a in sp.aliases)
        for k in keys:
            self._by_key.setdefault(k, []).append(sp.id)
        for n in {sp.name.lower()} | {a.lower() for a in sp.aliases}:
            if n:
                self._by_name.setdefault(n, []).append(sp.id)
        self._by_hill.setdefault(sp.hill, []).append(sp.id)
        return sp

    def _unindex(self, sp: Species) -> None:
        for index in (self._by_key, self._by_name, self._by_hill):
            for k in list(index):
                if sp.id in index[k]:
                    index[k] = [x for x in index[k] if x != sp.id]
                    if not index[k]:
                        del index[k]

    def __contains__(self, sid: str) -> bool:
        return sid in self.species

    def __len__(self) -> int:
        return len(self.species)

    def __iter__(self):
        return iter(self.species.values())

    # ------------------------------------------------------------------ lookup
    def get(self, text: str, phase: Optional[str] = None, allow_estimate: bool = True,
            online: bool = False) -> Species:
        """Find a species from user text.

        Accepted forms: ``H2O``, ``H2O(g)``, ``water``, ``C(graphite)``, ``Cu2+(aq)``,
        ``CuSO4·5H2O``, ``smiles:CCO`` (optionally ``smiles:CCO(l)``).  When the phase
        is omitted and several phases exist, the stable phase at 298.15 K is chosen.
        """
        text = text.strip()
        low = text.lower()
        for prefix in SMILES_PREFIXES:
            if low.startswith(prefix):
                smi, ph = split_phase(text[len(prefix):])
                if ph is not None and ph.lower() not in PHASE_ALIASES:
                    smi, ph = text[len(prefix):], None
                return self.from_smiles(smi, phase or ph, allow_estimate=allow_estimate)
        if text in self.species and phase is None:
            return self.species[text]
        body, ph = split_phase(text)
        if ph is not None and ph.lower() not in PHASE_ALIASES:
            # allotrope / named state, e.g. C(graphite)
            if text in self.species:
                return self.species[text]
            raise SpeciesNotFound("Unknown species %r" % text)
        want = phase or ph
        want = PHASE_ALIASES.get(want.lower(), want.lower()) if want else None
        cands = self._candidates(body)
        if not cands and online:
            from .online import resolve_name  # imported lazily: needs network

            info = resolve_name(body)
            return self.from_smiles(info["smiles"], want, allow_estimate=allow_estimate, name=body)
        if not cands:
            raise SpeciesNotFound(self._not_found_message(body))
        if want:
            cands = [c for c in cands if c.phase == want]
            if not cands:
                phases = sorted({c.phase for c in self._candidates(body)})
                raise SpeciesNotFound("%s is in the database only as phase(s) %s, not (%s)"
                                      % (body, ", ".join(phases), want))
        return self._choose(body, cands)

    def _candidates(self, body: str) -> List[Species]:
        ids: List[str] = []
        ids += self._by_key.get(body, [])
        ids += self._by_name.get(body.lower(), [])
        if not ids:
            try:
                elements, charge = parse_formula(body)
                ids += self._by_hill.get(hill_formula(elements, charge), [])
            except FormulaError:
                pass
        seen, out = set(), []
        for i in ids:
            if i not in seen:
                seen.add(i)
                out.append(self.species[i])
        return out

    def _choose(self, body: str, cands: List[Species]) -> Species:
        if len(cands) == 1:
            return cands[0]
        std = [c for c in cands if c.std]
        if len(std) == 1:
            return std[0]
        pool = std or cands
        # several isomers or phases: still ambiguous
        raise AmbiguousSpecies("%r is ambiguous; use one of: %s" % (body, ", ".join(
            "%s [%s]" % (c.id, c.name) for c in pool)))

    def _not_found_message(self, body: str) -> str:
        msg = "No species %r in the database." % body
        close = self.search(body, limit=5)
        if not close:
            import difflib

            names = {}
            for sp in self.species.values():
                for n in [sp.name, sp.formula, sp.id] + sp.aliases:
                    names.setdefault(n.lower(), sp)
            hits = difflib.get_close_matches(body.lower(), list(names), n=5, cutoff=0.75)
            close = list({names[h].id: names[h] for h in hits}.values())
        if close:
            msg += " Did you mean: %s?" % ", ".join(s.id for s in close)
        msg += " You can also give a SMILES string (smiles:...), use --online, or add data with --data FILE.csv."
        return msg

    def search(self, term: str, limit: int = 20) -> List[Species]:
        term_l = term.lower()
        hits = []
        for sp in self.species.values():
            hay = " ".join([sp.id, sp.formula, sp.name, sp.cas] + sp.aliases).lower()
            if term_l in hay:
                hits.append(sp)
        return hits[:limit]

    # ------------------------------------------------------------------ SMILES
    def match_smiles(self, smiles: str, phase: Optional[str] = None) -> List[Species]:
        """Database species whose stored SMILES describes the same molecule."""
        mol = parse_smiles(smiles)
        sig = mol.signature()
        out = []
        for sp in self.species.values():
            if not sp.smiles or sp.hill != mol.formula:
                continue
            try:
                if parse_smiles(sp.smiles).signature() == sig and (phase is None or sp.phase == phase):
                    out.append(sp)
            except SmilesError:
                continue
        return out

    def from_smiles(self, smiles: str, phase: Optional[str] = None, allow_estimate: bool = True,
                    name: str = "") -> Species:
        """Species for a SMILES string: tabulated data if the molecule is known, else a Joback estimate."""
        phase = PHASE_ALIASES.get(phase.lower(), phase.lower()) if phase else None
        try:
            mol = parse_smiles(smiles)
        except SmilesError as exc:
            raise SpeciesNotFound(str(exc))
        matches = self.match_smiles(smiles, phase)
        if matches:
            return self._choose(smiles, matches)
        if not allow_estimate:
            raise SpeciesNotFound("SMILES %s not in database (estimation disabled)" % smiles)
        if phase not in (None, "g"):
            raise SpeciesNotFound("SMILES %s (%s) is not tabulated; Joback estimates exist only for the "
                                  "ideal gas, use phase (g)" % (smiles, phase))
        sid = "smiles:%s(g)" % smiles
        if sid in self.species:
            return self.species[sid]
        try:
            est = estimate(mol)
        except JobackError as exc:
            raise SpeciesNotFound("Cannot estimate %s: %s" % (smiles, exc))
        if est.S_J is None:
            raise SpeciesNotFound("Cannot estimate S° for %s (missing Joback ΔG contribution)" % smiles)
        cp_model = None
        cp298 = None
        if est.cp_coeffs is not None:
            a, b, c, d = est.cp_coeffs
            cp_model = PiecewiseCp([Polynomial(a, b, c, d, tmin=CP_VALID_RANGE[0], tmax=CP_VALID_RANGE[1])],
                                   label="Joback polynomial")
            cp298 = est.cp(T_REF)
        groups = ", ".join("%d×[%s]" % (n, g) for g, n in sorted(est.groups.items()))
        sp = Species(id=sid, formula=mol.formula, phase="g", elements=mol.element_counts(),
                     charge=mol.charge, dHf=est.dHf_kJ * 1000.0, S=est.S_J, cp298=cp298, cp_model=cp_model,
                     name=name or smiles, smiles=smiles, source="Joback estimate", estimated=True,
                     notes=["Joback groups: " + groups])
        self.add(sp)
        return sp

    # ------------------------------------------------------------------ formation properties
    def formation(self, sp: Species, T: float = T_REF, mode: str = "full") -> Dict[str, Optional[float]]:
        """ΔfH°(T) and ΔfG°(T) (J/mol) from the species and its elements' reference states.

        Ions follow the usual convention ΔfH°, ΔfG°, S° of H+(aq) = 0: the formation
        reaction of an ion of charge z is balanced with z H+(aq) and z/2 H2(g).
        """
        dH = sp.H(T, mode)
        dG = sp.G(T, mode)
        if sp.charge:
            if "H+(aq)" not in self.species or "H2(g)" not in self.species:
                return {"dHf": None, "dGf": None}
            h2, hp = self.species["H2(g)"], self.species["H+(aq)"]
            z = sp.charge
            dH += z / 2.0 * h2.H(T, mode) - z * hp.H(T, mode)
            dG += z / 2.0 * h2.G(T, mode) - z * hp.G(T, mode)
        for el, n in sp.elements.items():
            ref = ELEMENT_REFERENCE.get(el)
            if ref is None or ref[0] not in self.species:
                return {"dHf": None, "dGf": None}
            ref_sp = self.species[ref[0]]
            dH -= n / ref[1] * ref_sp.H(T, mode)
            dG -= n / ref[1] * ref_sp.G(T, mode)
        return {"dHf": dH, "dGf": dG}

    def ids(self) -> Iterable[str]:
        return self.species.keys()
