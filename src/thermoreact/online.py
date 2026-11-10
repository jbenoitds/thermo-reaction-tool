"""Optional online lookups (PubChem PUG-REST), using only the standard library.

Used when a species name is not in the local database and ``--online`` is given:
the name is resolved to a SMILES string, which is then matched against the
database or estimated with the Joback method.  Results are cached per session.
"""
from __future__ import annotations

import json
import urllib.error
import urllib.parse
import urllib.request
from typing import Callable, Dict, Optional

from .database import SpeciesNotFound

PUBCHEM_URL = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/"
               "IsomericSMILES,CanonicalSMILES,MolecularFormula,IUPACName/JSON")
TIMEOUT_S = 10.0
_cache: Dict[str, Dict[str, str]] = {}

#: Replaceable for testing: a function url -> response bytes
_fetch: Optional[Callable[[str], bytes]] = None


def _http_get(url: str) -> bytes:
    if _fetch is not None:
        return _fetch(url)
    req = urllib.request.Request(url, headers={"User-Agent": "thermoreact/0.1 (educational tool)"})
    with urllib.request.urlopen(req, timeout=TIMEOUT_S) as resp:  # nosec - fixed https host
        return resp.read()


def resolve_name(name: str) -> Dict[str, str]:
    """Return ``{'smiles', 'formula', 'iupac_name', 'cid'}`` for a compound name via PubChem."""
    key = name.strip().lower()
    if key in _cache:
        return _cache[key]
    url = PUBCHEM_URL.format(name=urllib.parse.quote(name.strip(), safe=""))
    try:
        raw = _http_get(url)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            raise SpeciesNotFound("PubChem does not know %r" % name)
        raise SpeciesNotFound("PubChem request failed (%s)" % exc)
    except (urllib.error.URLError, OSError) as exc:
        raise SpeciesNotFound("Could not reach PubChem (%s); check your connection or work offline" % exc)
    try:
        props = json.loads(raw.decode("utf-8"))["PropertyTable"]["Properties"][0]
    except (ValueError, KeyError, IndexError):
        raise SpeciesNotFound("Unexpected PubChem response for %r" % name)
    # PubChem has used several key names for SMILES over the years
    smiles = (props.get("IsomericSMILES") or props.get("SMILES") or props.get("CanonicalSMILES")
              or props.get("ConnectivitySMILES"))
    if not smiles:
        raise SpeciesNotFound("PubChem returned no SMILES for %r" % name)
    info = {"smiles": smiles, "formula": props.get("MolecularFormula", ""),
            "iupac_name": props.get("IUPACName", ""), "cid": str(props.get("CID", ""))}
    _cache[key] = info
    return info
