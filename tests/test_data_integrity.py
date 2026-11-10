"""Consistency checks of the bundled data (catches typos in coefficients and tables)."""
import json
import math
import os

import pytest

from thermoreact import Database
from thermoreact.constants import ELEMENT_REFERENCE, T_REF
from thermoreact.database import DATA_DIR

DB = Database.default()
with open(os.path.join(DATA_DIR, "shomate.json"), encoding="utf-8") as fh:
    SHOMATE = json.load(fh)["species"]


def _H(s, T):  # NIST H - H298 (kJ/mol)
    t = T / 1000
    return s["A"] * t + s["B"] * t * t / 2 + s["C"] * t ** 3 / 3 + s["D"] * t ** 4 / 4 - s["E"] / t + s["F"] - s["H"]


def _S(s, T):
    t = T / 1000
    return s["A"] * math.log(t) + s["B"] * t + s["C"] * t * t / 2 + s["D"] * t ** 3 / 3 - s["E"] / (2 * t * t) + s["G"]


NIST_SPECIES = [k for k, v in SHOMATE.items() if "F" in v["segments"][0]]


@pytest.mark.parametrize("sid", NIST_SPECIES)
def test_shomate_reproduces_298_values(sid):
    seg = SHOMATE[sid]["segments"][0]
    sp = DB.species[sid]
    assert abs(_H(seg, T_REF)) < 0.015          # H - H298 = 0 at 298.15 K (kJ/mol)
    assert _S(seg, T_REF) == pytest.approx(sp.S, abs=0.1)
    assert sp.cp_model.cp(T_REF) == pytest.approx(sp.cp298, rel=0.01)
    assert seg["H"] == pytest.approx(sp.dHf / 1000, abs=0.03)


@pytest.mark.parametrize("sid", [k for k in NIST_SPECIES if len(SHOMATE[k]["segments"]) > 1])
def test_shomate_segments_join(sid):
    segs = SHOMATE[sid]["segments"]
    for a, b in zip(segs, segs[1:]):
        T = a["tmax"]
        assert abs(_H(a, T) - _H(b, T)) < 0.03
        assert abs(_S(a, T) - _S(b, T)) < 0.03


@pytest.mark.parametrize("sid", NIST_SPECIES)
def test_integration_matches_nist_tables(sid):
    """Our piecewise integration of Cp reproduces NIST's tabulated H(T) - H298 and S(T)."""
    sp = DB.species[sid]
    for T in (400.0, 900.0, 1500.0, 2500.0):
        seg = next((s for s in SHOMATE[sid]["segments"] if s["tmin"] <= T <= s["tmax"]), None)
        if seg is None or "F" not in seg:
            continue
        assert sp.int_cp(T) / 1000 == pytest.approx(_H(seg, T), abs=0.05)
        # 0.1 J/(mol K): the tabulated S° (CODATA) may differ slightly from the Shomate-fit value (Chase)
        assert sp.S_at(T) == pytest.approx(_S(seg, T), abs=0.1)


def test_elements_in_reference_state_have_zero_formation_values():
    for el, (sid, n) in ELEMENT_REFERENCE.items():
        sp = DB.species[sid]
        assert sp.dHf == 0, sid
        f = DB.formation(sp, 1000.0)
        assert f["dHf"] == pytest.approx(0, abs=1e-6) and f["dGf"] == pytest.approx(0, abs=1e-6)


def test_every_species_is_complete():
    for sp in DB:
        assert sp.dHf is not None and sp.S is not None, sp.id
        assert sp.phase in ("g", "l", "s", "aq"), sp.id
        if sp.phase != "aq":
            assert sp.S > 0, sp.id
            assert sp.cp298 is not None and sp.cp298 > 0, sp.id
        assert sp.source, sp.id
