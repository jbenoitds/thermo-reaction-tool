"""Parsing of chemical formulas such as ``H2O``, ``Ca(OH)2``, ``CuSO4·5H2O`` or ``SO4^2-``.

The parser returns element counts and the net charge.  Charges may be written

* with a caret: ``SO4^2-``, ``Cu^2+``, ``Fe^3+``
* as a single trailing sign (charge of 1): ``Cl-``, ``NH4+``, ``NO3-``
* as repeated signs: ``Cu++``, ``SO4--``
* with the sign before the number: ``Cu+2``, ``SO4-2``
* in brackets: ``Cu[2+]``

A bare ``Cu2+`` is ambiguous (Cu2 with 1+ or Cu with 2+); the parser reads it
literally as Cu2 with charge 1+, but the species database lists ``Cu2+`` etc. as
aliases of the ions, so in a reaction they resolve to the expected ions.
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Dict, Tuple

from .constants import ATOMIC_MASSES


class FormulaError(ValueError):
    """Raised when a chemical formula cannot be parsed."""


_ELEMENT_RE = re.compile(r"[A-Z][a-z]?")
_HYDRATE_SEPARATORS = ("·", "•", "∙", "*", ".")
_SUBSCRIPTS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
_SUPERSCRIPTS = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻", "0123456789+-")


def _split_charge(text: str) -> Tuple[str, int]:
    """Separate a trailing charge specification from the formula body."""
    text = text.strip()
    # Cu[2+] / SO4[2-]
    m = re.search(r"\[(\d*)([+-])\]$", text)
    if m:
        n = int(m.group(1) or 1)
        return text[: m.start()], n if m.group(2) == "+" else -n
    # SO4^2- / Cu^2+ / Fe^+3 / Cl^-
    m = re.search(r"\^(\d*)([+-])(\d*)$", text)
    if m:
        if m.group(1) and m.group(3):
            raise FormulaError("Bad charge in %r" % text)
        n = int(m.group(1) or m.group(3) or 1)
        return text[: m.start()], n if m.group(2) == "+" else -n
    # Cu++ / SO4--
    m = re.search(r"(\+{2,}|-{2,})$", text)
    if m:
        n = len(m.group(1))
        return text[: m.start()], n if m.group(1)[0] == "+" else -n
    # Cu+2 / SO4-2
    m = re.search(r"([+-])(\d+)$", text)
    if m:
        n = int(m.group(2))
        return text[: m.start()], n if m.group(1) == "+" else -n
    # Cl- / NH4+
    if text.endswith("+") or text.endswith("-"):
        return text[:-1], 1 if text[-1] == "+" else -1
    return text, 0


def _parse_group(text: str, pos: int, closing: str = "") -> Tuple[Counter, int]:
    counts: Counter = Counter()
    while pos < len(text):
        ch = text[pos]
        if ch in "([{":
            close = {"(": ")", "[": "]", "{": "}"}[ch]
            inner, pos = _parse_group(text, pos + 1, close)
            mult, pos = _read_number(text, pos)
            for el, n in inner.items():
                counts[el] += n * mult
        elif closing and ch == closing:
            return counts, pos + 1
        elif ch in ")]}":
            raise FormulaError("Unbalanced bracket in %r" % text)
        else:
            m = _ELEMENT_RE.match(text, pos)
            if not m:
                raise FormulaError("Unexpected character %r in formula %r" % (ch, text))
            el = m.group(0)
            if el not in ATOMIC_MASSES:
                raise FormulaError("Unknown element %r in formula %r" % (el, text))
            pos = m.end()
            mult, pos = _read_number(text, pos)
            counts[el] += mult
    if closing:
        raise FormulaError("Missing %r in formula %r" % (closing, text))
    return counts, pos


def _read_number(text: str, pos: int) -> Tuple[int, int]:
    m = re.compile(r"\d+").match(text, pos)
    if m:
        return int(m.group(0)), m.end()
    return 1, pos


def parse_formula(formula: str) -> Tuple[Dict[str, int], int]:
    """Return ``(element_counts, charge)`` for a chemical formula.

    >>> parse_formula("Ca(OH)2")
    ({'Ca': 1, 'O': 2, 'H': 2}, 0)
    >>> parse_formula("CuSO4·5H2O")[0]["O"]
    9
    >>> parse_formula("SO4^2-")
    ({'S': 1, 'O': 4}, -2)
    """
    if not formula or not formula.strip():
        raise FormulaError("Empty formula")
    text = formula.strip().translate(_SUBSCRIPTS).translate(_SUPERSCRIPTS).replace(" ", "")
    if text in ("e-", "e^-", "e"):
        return {}, -1
    body, charge = _split_charge(text)
    total: Counter = Counter()
    parts = [body]
    for sep in _HYDRATE_SEPARATORS:
        if sep in body:
            parts = body.split(sep)
            break
    for part in parts:
        if not part:
            raise FormulaError("Empty part in formula %r" % formula)
        m = re.match(r"^(\d+)(?=[A-Z(\[])", part)
        mult = 1
        if m:
            mult = int(m.group(1))
            part = part[m.end():]
        counts, _ = _parse_group(part, 0)
        for el, n in counts.items():
            total[el] += n * mult
    if not total:
        raise FormulaError("No elements found in %r" % formula)
    return dict(total), charge


def hill_formula(counts: Dict[str, int], charge: int = 0) -> str:
    """Hill-system formula (C first, then H, then alphabetical; alphabetical if no C)."""
    def part(el: str) -> str:
        n = counts[el]
        return el + (str(n) if n != 1 else "")

    els = [e for e, n in counts.items() if n]
    if "C" in els:
        order = ["C"] + (["H"] if "H" in els else []) + sorted(e for e in els if e not in ("C", "H"))
    else:
        order = sorted(els)
    s = "".join(part(e) for e in order)
    if charge:
        s += "^" + (str(abs(charge)) if abs(charge) != 1 else "") + ("+" if charge > 0 else "-")
    return s


def molar_mass(counts: Dict[str, int]) -> float:
    """Molar mass in g/mol."""
    return sum(ATOMIC_MASSES[el] * n for el, n in counts.items())
