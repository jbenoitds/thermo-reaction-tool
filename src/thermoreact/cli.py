"""Command-line interface: ``thermoreact <command> ...`` (also ``python -m thermoreact``)."""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from typing import Dict, List, Optional, Sequence

from . import __version__
from .constants import R, T_REF
from .database import AmbiguousSpecies, Database, SpeciesNotFound
from .equilibrium import k_from_extent, solve_equilibrium
from .formula import FormulaError
from .joback import JobackError, estimate
from .reaction import Reaction, ReactionError
from .relations import (E0_from_K, K_from_dG, K_from_E0, dG_from_K, vant_hoff_dH, vant_hoff_K2,
                        vant_hoff_T2)
from .smiles import SmilesError, parse_smiles

ERRORS = (SpeciesNotFound, AmbiguousSpecies, ReactionError, FormulaError, SmilesError, JobackError, ValueError,
          OSError)


def _utf8_output() -> None:
    """Make sure symbols such as Δ, ° and → can be printed.

    A Windows console already accepts Unicode, but output redirected to a file
    (``thermoreact ... > results.txt``) would otherwise use the legacy code page
    (e.g. cp1252), which has no Δ or →.
    """
    for stream in (sys.stdout, sys.stderr):
        enc = (getattr(stream, "encoding", None) or "").lower().replace("-", "").replace("_", "")
        if enc != "utf8" and hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8")
            except (ValueError, OSError):  # pragma: no cover - exotic streams
                pass


# ---------------------------------------------------------------------------- helpers
def fmt_num(x: Optional[float], digits: int = 4) -> str:
    if x is None:
        return "n/a"
    if isinstance(x, float) and (math.isinf(x) or math.isnan(x)):
        return str(x)
    if x == 0:
        return "0"
    a = abs(x)
    if a >= 1e5 or a < 1e-3:
        return "%.*e" % (digits - 1, x)
    return "%.*g" % (digits, x)


def fmt_K(lnK: float) -> str:
    log10 = lnK / math.log(10)
    if -300 < log10 < 300:
        return fmt_num(math.exp(lnK))
    exp = math.floor(log10)
    return "%.3fe%+d" % (10 ** (log10 - exp), exp)


def table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    widths = [max(len(str(h)), *(len(str(r[i])) for r in rows)) if rows else len(str(h))
              for i, h in enumerate(headers)]
    line = "  ".join(str(h).rjust(w) for h, w in zip(headers, widths))
    sep = "  ".join("-" * w for w in widths)
    body = ["  ".join(str(c).rjust(w) for c, w in zip(r, widths)) for r in rows]
    return "\n".join([line, sep] + body)


def to_kelvin(values: Sequence[float], celsius: bool) -> List[float]:
    return [v + 273.15 for v in values] if celsius else list(values)


def temperatures(args) -> List[float]:
    Ts: List[float] = []
    if getattr(args, "T", None):
        Ts += to_kelvin(args.T, args.celsius)
    if getattr(args, "trange", None):
        start, stop, step = args.trange
        if step <= 0:
            raise ValueError("--trange step must be positive")
        n = int(math.floor((stop - start) / step + 1e-9)) + 1
        Ts += to_kelvin([start + i * step for i in range(n)], args.celsius)
    if not Ts:
        Ts = [T_REF]
    for T in Ts:
        if T <= 0:
            raise ValueError("Temperatures must be above 0 K (got %g K)" % T)
    return Ts


def load_db(args) -> Database:
    db = Database.default()
    if getattr(args, "data", None):
        db = db.copy()
        for path in args.data:
            added = db.load_csv(path)
            print("Loaded %d species from %s" % (len(added), path), file=sys.stderr)
    return db


def parse_amounts(items: Optional[Sequence[str]]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for item in items or []:
        for part in item.split(","):
            part = part.strip()
            if not part:
                continue
            if "=" not in part:
                raise ValueError("Initial amounts look like N2=1,H2=3 (got %r)" % part)
            k, v = part.rsplit("=", 1)
            out[k.strip()] = out.get(k.strip(), 0.0) + float(v)
    return out


def dedupe(items: Sequence[str]) -> List[str]:
    seen, out = set(), []
    for it in items:
        if it not in seen:
            seen.add(it)
            out.append(it)
    return out


# ---------------------------------------------------------------------------- commands
def cmd_reaction(args) -> int:
    db = load_db(args)
    rxn = Reaction.from_string(args.equation, db, balance=args.balance, online=args.online)
    Ts = temperatures(args)
    results = rxn.table(Ts, args.cp)
    if args.json:
        out = {"reaction": rxn.equation("->"), "delta_n_gas": float(rxn.delta_n_gas), "cp_mode": args.cp,
               "species": [species_dict(sp, nu) for nu, sp in rxn.terms],
               "results": [r.as_dict() for r in results],
               "warnings": dedupe([w for r in results for w in r.warnings]), "notes": rxn.notes}
        if args.crossover:
            out["crossover_T_K"] = rxn.crossover_temperature(args.cp)
        print(json.dumps(out, indent=2, ensure_ascii=False))
    else:
        print("Reaction: %s" % rxn.equation())
        for n in rxn.notes:
            print("  note: %s" % n)
        print("Δn(gas) = %s   Cp treatment: %s   standard state: 1 bar" % (rxn.delta_n_gas, args.cp))
        print()
        print(species_table(rxn))
        print()
        headers = ["T/K", "ΔrH°/kJ", "ΔrS°/(J/K)", "ΔrG°/kJ", "ΔrCp°/(J/K)", "K", "log10 K"]
        if args.kc:
            headers.append("Kc")
        rows = []
        for r in results:
            row = ["%.2f" % r.T, "%.3f" % (r.dH / 1000), "%.3f" % r.dS, "%.3f" % (r.dG / 1000),
                   "n/a" if r.dCp is None else "%.3f" % r.dCp, fmt_K(r.lnK), "%.3f" % r.log10K]
            if args.kc:
                row.append(fmt_K(math.log(r.Kc)) if r.Kc and r.Kc > 0 else "n/a")
            rows.append(row)
        print(table(headers, rows))
        if args.crossover:
            roots = rxn.crossover_temperature(args.cp)
            print()
            if roots:
                print("ΔrG° = 0 (K = 1) at: " + ", ".join("%.1f K (%.1f °C)" % (t, t - 273.15) for t in roots))
            else:
                print("ΔrG° does not change sign between 100 and 6000 K")
        warns = dedupe([w for r in results for w in r.warnings])
        if warns:
            print()
            for w in warns:
                print("warning: " + w)
    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["T_K", "dH_kJ_per_mol", "dS_J_per_mol_K", "dG_kJ_per_mol", "dCp_J_per_mol_K",
                        "lnK", "log10K", "K"])
            for r in results:
                w.writerow([r.T, r.dH / 1000, r.dS, r.dG / 1000, r.dCp, r.lnK, r.log10K, r.K])
        print("Saved table to %s" % args.csv, file=sys.stderr)
    if args.plot:
        from .plotting import plot_reaction

        Tplot = Ts if len(Ts) > 2 else [T_REF + i * (1200 - T_REF) / 60 for i in range(61)]
        plot_reaction(rxn, Tplot, args.plot, args.cp)
        print("Saved plot to %s" % args.plot, file=sys.stderr)
    return 0


def species_dict(sp, nu=None) -> Dict:
    d = {"id": sp.id, "name": sp.name, "formula": sp.formula, "phase": sp.phase,
         "dHf_kJ_per_mol": None if sp.dHf is None else sp.dHf / 1000, "S_J_per_mol_K": sp.S,
         "Cp298_J_per_mol_K": sp.cp(T_REF), "cp_model": sp.describe_cp(), "source": sp.source,
         "estimated": sp.estimated, "notes": sp.notes}
    if nu is not None:
        d["nu"] = float(nu)
    return d


def species_table(rxn: Reaction) -> str:
    rows = []
    for nu, sp in rxn.terms:
        rows.append([("%g" % float(nu)), sp.id, sp.name[:28], "%.2f" % (sp.dHf / 1000), "%.2f" % sp.S,
                     fmt_num(sp.cp(T_REF)), sp.describe_cp()[:34], sp.source])
    return table(["ν", "species", "name", "ΔfH°/kJ", "S°/(J/K)", "Cp°/(J/K)", "Cp(T) model", "source"], rows)


def cmd_species(args) -> int:
    db = load_db(args)
    Ts = temperatures(args)
    out = []
    for name in args.names:
        sp = db.get(name, online=args.online)
        if args.json:
            d = species_dict(sp)
            d["at_T"] = []
            for T in Ts:
                f = db.formation(sp, T, args.cp)
                d["at_T"].append({"T_K": T, "Cp": sp.cp(T, args.cp), "S": sp.S_at(T, args.cp),
                                  "H_minus_H298_kJ": sp.int_cp(T, args.cp) / 1000,
                                  "dHf_kJ": None if f["dHf"] is None else f["dHf"] / 1000,
                                  "dGf_kJ": None if f["dGf"] is None else f["dGf"] / 1000})
            out.append(d)
            continue
        print("%s — %s (%s)" % (sp.id, sp.name, sp.phase_name))
        print("  formula %s   Hill %s   M = %.3f g/mol%s" % (sp.formula, sp.hill, sp.molar_mass,
                                                              "   CAS %s" % sp.cas if sp.cas else ""))
        if sp.smiles:
            print("  SMILES  %s" % sp.smiles)
        print("  source: %s%s" % (sp.source, "  (ESTIMATED)" if sp.estimated else ""))
        print("  Cp(T):  %s" % sp.describe_cp())
        for n in sp.notes:
            print("  note:   %s" % n)
        rows = []
        warns = []
        for T in Ts:
            f = db.formation(sp, T, args.cp)
            rows.append(["%.2f" % T, fmt_num(sp.cp(T, args.cp)), "%.3f" % sp.S_at(T, args.cp),
                         "%.3f" % (sp.int_cp(T, args.cp) / 1000),
                         "n/a" if f["dHf"] is None else "%.3f" % (f["dHf"] / 1000),
                         "n/a" if f["dGf"] is None else "%.3f" % (f["dGf"] / 1000)])
            warns += sp.warnings_at(T, args.cp)
        print(table(["T/K", "Cp°/(J/K)", "S°/(J/K)", "H°-H°298/kJ", "ΔfH°/kJ", "ΔfG°/kJ"], rows))
        for w in dedupe(warns):
            print("  warning: " + w)
        print()
    if args.json:
        print(json.dumps(out, indent=2, ensure_ascii=False))
    return 0


def cmd_list(args) -> int:
    db = load_db(args)
    items = db.search(args.search, limit=10 ** 6) if args.search else list(db)
    if args.phase:
        items = [s for s in items if s.phase == args.phase]
    rows = [[s.id, s.name[:34], s.phase, "%.2f" % (s.dHf / 1000), "%.2f" % s.S, fmt_num(s.cp298),
             "T-dependent" if (s.cp_model and len(s.cp_model.segments) and s.cp_model.label != "constant Cp(298 K)")
             else "constant", s.source] for s in items]
    print(table(["id", "name", "phase", "ΔfH°/kJ", "S°/(J/K)", "Cp°/(J/K)", "Cp", "source"], rows))
    print("\n%d species" % len(rows))
    return 0


def cmd_estimate(args) -> int:
    db = load_db(args)
    mol = parse_smiles(args.smiles)
    est = estimate(mol)
    if args.json:
        print(json.dumps(est.as_dict(), indent=2, ensure_ascii=False))
        return 0
    print("Joback estimate for %s  (%s, ideal gas)" % (args.smiles, mol.formula))
    print(table(["group", "count"], [[g, str(n)] for g, n in sorted(est.groups.items())]))
    print()
    print("  ΔfH°(298 K) = %.2f kJ/mol" % est.dHf_kJ)
    if est.dGf_kJ is not None:
        print("  ΔfG°(298 K) = %.2f kJ/mol" % est.dGf_kJ)
        print("  S°(298 K)   = %.2f J/(mol K)" % est.S_J)
    if est.cp_coeffs:
        a, b, c, d = est.cp_coeffs
        print("  Cp°(T) = %.4g + %.4g T + %.4g T^2 + %.4g T^3  J/(mol K)   (≈273-1000 K)" % (a, b, c, d))
        print("  Cp°(298 K) = %.2f J/(mol K)" % est.cp(T_REF))
    known = db.match_smiles(args.smiles, "g")
    if known:
        sp = known[0]
        print()
        print("Tabulated value for comparison (%s, %s): ΔfH° = %.2f kJ/mol, S° = %.2f J/(mol K)"
              % (sp.id, sp.source, sp.dHf / 1000, sp.S))
        print("  error: ΔfH° %+.2f kJ/mol, S° %+.2f J/(mol K)" % (est.dHf_kJ - sp.dHf / 1000,
                                                                 (est.S_J or 0) - sp.S))
    return 0


def cmd_equilibrium(args) -> int:
    db = load_db(args)
    rxn = Reaction.from_string(args.equation, db, balance=args.balance, online=args.online,
                               abstract=args.abstract)
    initial = parse_amounts(args.init)
    if not initial:
        raise ValueError("Give initial amounts, e.g. --init N2=1,H2=3")
    T = to_kelvin([args.T], args.celsius)[0] if args.T is not None else None
    if args.extent is not None:
        K = k_from_extent(rxn, initial, args.extent, args.P, args.inert)
        res = {"reaction": rxn.equation("->"), "extent_mol": args.extent, "P_bar": args.P, "K": K}
        print_lines = ["Reaction: %s" % rxn.equation(),
                       "Extent ξ = %g mol at P = %g bar  →  K = %s" % (args.extent, args.P, fmt_num(K, 5))]
        if T is not None:
            dG = dG_from_K(K, T)
            res.update({"T_K": T, "dG_kJ_per_mol": dG / 1000})
            print_lines.append("ΔrG° = −RT ln K = %.3f kJ/mol at %.2f K" % (dG / 1000, T))
        if args.json:
            print(json.dumps(res, indent=2, ensure_ascii=False))
        else:
            print("\n".join(print_lines))
        return 0
    if T is None and args.K is None:
        raise ValueError("Give --T (K from data) or --K")
    if args.abstract and args.K is None:
        raise ValueError("Abstract species have no data: give --K")
    res = solve_equilibrium(rxn, initial, T=T, P_bar=args.P, K=args.K, cp_mode=args.cp, inert_gas=args.inert)
    if args.json:
        d = res.as_dict()
        d["reaction"] = rxn.equation("->")
        print(json.dumps(d, indent=2, ensure_ascii=False))
        return 0
    print("Reaction: %s" % rxn.equation())
    print("T = %s   P = %g bar   K = %s   ξ(eq) = %.6g mol" % (
        "n/a" if T is None else "%.2f K" % T, args.P, fmt_num(res.K, 5), res.extent))
    x = res.mole_fractions()
    rows = []
    for nu, sp in rxn.terms:
        rows.append([sp.id, fmt_num(res.initial.get(sp.id, 0.0)), fmt_num(res.amounts[sp.id], 5),
                     fmt_num(x.get(sp.id)) if sp.is_gas else "(pure %s)" % sp.phase,
                     fmt_num(x.get(sp.id, 0) * args.P) if sp.is_gas else "-"])
    print(table(["species", "n0/mol", "n_eq/mol", "x(gas)", "p/bar"], rows))
    if args.inert:
        print("  plus %g mol inert gas" % args.inert)
    for sid, c in res.conversion().items():
        print("Conversion of %s: %.2f %%" % (sid, 100 * c))
    for n in dedupe(res.notes):
        print("note: " + n)
    return 0


def cmd_vanthoff(args) -> int:
    T1 = to_kelvin([args.T1], args.celsius)[0]
    if args.K1 is None and args.dG1 is None:
        raise ValueError("Give --K1 or --dG1")
    K1 = args.K1 if args.K1 is not None else K_from_dG(args.dG1 * 1000, T1)
    T2 = to_kelvin([args.T2], args.celsius)[0] if args.T2 is not None else None
    lines = ["K1 = %s at T1 = %.2f K" % (fmt_num(K1, 5), T1)]
    result: Dict[str, float] = {"K1": K1, "T1_K": T1}
    if args.dH is not None and T2 is not None:
        K2 = vant_hoff_K2(K1, T1, T2, args.dH * 1000)
        lines.append("K2 = %s at T2 = %.2f K   (ΔrH° = %g kJ/mol assumed constant)" % (fmt_num(K2, 5), T2, args.dH))
        result.update({"K2": K2, "T2_K": T2, "dH_kJ_per_mol": args.dH})
    elif args.K2 is not None and T2 is not None:
        dH = vant_hoff_dH(K1, T1, args.K2, T2)
        lines.append("ΔrH° = %.3f kJ/mol   (from K2 = %s at %.2f K)" % (dH / 1000, fmt_num(args.K2, 5), T2))
        result.update({"K2": args.K2, "T2_K": T2, "dH_kJ_per_mol": dH / 1000})
    elif args.K2 is not None and args.dH is not None:
        T2 = vant_hoff_T2(K1, T1, args.K2, args.dH * 1000)
        lines.append("K reaches %s at T2 = %.2f K (%.2f °C)" % (fmt_num(args.K2, 5), T2, T2 - 273.15))
        result.update({"K2": args.K2, "T2_K": T2, "dH_kJ_per_mol": args.dH})
    else:
        raise ValueError("Give two of --T2, --K2, --dH")
    print(json.dumps(result, indent=2) if args.json else "\n".join(lines))
    return 0


def cmd_electrochem(args) -> int:
    T = to_kelvin([args.T], args.celsius)[0]
    if args.E0 is not None:
        K = K_from_E0(args.E0, args.n, T)
        dG = -args.n * 96485.33212 * args.E0
        res = {"E0_V": args.E0, "n": args.n, "T_K": T, "K": K, "dG_kJ_per_mol": dG / 1000}
        text = "E° = %g V, n = %d, T = %.2f K  →  ΔrG° = %.3f kJ/mol, K = %s" % (
            args.E0, args.n, T, dG / 1000, fmt_K(args.n * 96485.33212 * args.E0 / (R * T)))
    elif args.K is not None:
        E0 = E0_from_K(args.K, args.n, T)
        res = {"K": args.K, "n": args.n, "T_K": T, "E0_V": E0}
        text = "K = %s, n = %d, T = %.2f K  →  E° = %.4f V" % (fmt_num(args.K, 5), args.n, T, E0)
    else:
        raise ValueError("Give --E0 or --K")
    print(json.dumps(res, indent=2) if args.json else text)
    return 0


# ---------------------------------------------------------------------------- parser
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="thermoreact",
        description="Reaction thermochemistry & equilibrium calculator (ΔH°, ΔS°, ΔG°, K vs T).",
        epilog='Example: thermoreact reaction "N2 + 3 H2 -> 2 NH3" -T 298.15 500 700 --crossover')
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    sub = p.add_subparsers(dest="command")

    def common(sp, temps=True):
        if temps:
            sp.add_argument("-T", "--T", type=float, nargs="+", metavar="T", help="temperature(s), K (default 298.15)")
            sp.add_argument("--trange", type=float, nargs=3, metavar=("START", "STOP", "STEP"),
                            help="temperature range, e.g. --trange 300 1500 100")
        sp.add_argument("--celsius", action="store_true", help="temperatures given in °C instead of K")
        sp.add_argument("--cp", choices=["full", "constant", "none"], default="full",
                        help="heat capacity: T-dependent (full, default), Cp(298 K) constant, "
                             "or ignored (ΔH°, ΔS° constant)")
        sp.add_argument("--data", action="append", metavar="FILE.csv", help="extra species data (CSV), repeatable")
        sp.add_argument("--online", action="store_true", help="resolve unknown names through PubChem")
        sp.add_argument("--json", action="store_true", help="machine-readable JSON output")

    r = sub.add_parser("reaction", aliases=["rxn"], help="ΔrH°, ΔrS°, ΔrG°, K at one or more temperatures")
    r.add_argument("equation", help='e.g. "CH4 + 2 O2 -> CO2 + 2 H2O(l)" or "smiles:CCO -> smiles:C=C + H2O(g)"')
    common(r)
    r.add_argument("--balance", action="store_true", help="balance the equation automatically")
    r.add_argument("--crossover", action="store_true", help="find T where ΔrG° = 0 (K = 1)")
    r.add_argument("--kc", action="store_true", help="also print Kc (c° = 1 mol/L)")
    r.add_argument("--csv", metavar="FILE", help="write the results table to CSV")
    r.add_argument("--plot", metavar="FILE.png", help="save ΔG°(T) and van 't Hoff plots (needs matplotlib)")
    r.set_defaults(func=cmd_reaction)

    s = sub.add_parser("species", aliases=["sp"], help="properties of species: S°, Cp°, ΔfH°, ΔfG° vs T")
    s.add_argument("names", nargs="+", help="formulas, names or smiles:... strings")
    common(s)
    s.set_defaults(func=cmd_species)

    li = sub.add_parser("list", help="list or search the species database")
    li.add_argument("search", nargs="?", help="text to search in ids, names, formulas, CAS")
    li.add_argument("--phase", choices=["g", "l", "s", "aq"])
    li.add_argument("--data", action="append", metavar="FILE.csv")
    li.set_defaults(func=cmd_list)

    e = sub.add_parser("estimate", help="Joback group-contribution estimate from SMILES")
    e.add_argument("smiles")
    e.add_argument("--json", action="store_true")
    e.add_argument("--data", action="append", metavar="FILE.csv")
    e.set_defaults(func=cmd_estimate)

    q = sub.add_parser("equilibrium", aliases=["eq"], help="equilibrium composition, or K from a measured extent")
    q.add_argument("equation")
    q.add_argument("--init", action="append", metavar="SPECIES=MOL", help="initial amounts, e.g. N2=1,H2=3")
    q.add_argument("-T", "--T", type=float, help="temperature (K); K is computed from the data")
    q.add_argument("-P", "--P", type=float, default=1.0, help="total pressure, bar (default 1)")
    q.add_argument("--K", type=float, help="use this equilibrium constant instead of the data")
    q.add_argument("--extent", type=float, help="measured extent ξ (mol) → compute K (and ΔrG° if -T given)")
    q.add_argument("--inert", type=float, default=0.0, help="moles of inert gas present")
    q.add_argument("--abstract", action="store_true", help="species are generic gases without data (e.g. A + B -> C)")
    q.add_argument("--balance", action="store_true")
    common(q, temps=False)
    q.set_defaults(func=cmd_equilibrium)

    v = sub.add_parser("vanthoff", help="van 't Hoff equation: K2, ΔrH° or T2 (constant ΔrH°)")
    v.add_argument("--K1", type=float)
    v.add_argument("--dG1", type=float, help="ΔrG° at T1 in kJ/mol (instead of --K1)")
    v.add_argument("--T1", type=float, required=True)
    v.add_argument("--T2", type=float)
    v.add_argument("--K2", type=float)
    v.add_argument("--dH", type=float, help="ΔrH° in kJ/mol")
    v.add_argument("--celsius", action="store_true")
    v.add_argument("--json", action="store_true")
    v.set_defaults(func=cmd_vanthoff)

    c = sub.add_parser("electrochem", help="K and ΔrG° from a standard cell potential (or E° from K)")
    c.add_argument("--E0", type=float, help="standard cell potential, V")
    c.add_argument("--K", type=float)
    c.add_argument("--n", type=int, required=True, help="electrons transferred")
    c.add_argument("-T", "--T", type=float, default=T_REF)
    c.add_argument("--celsius", action="store_true")
    c.add_argument("--json", action="store_true")
    c.set_defaults(func=cmd_electrochem)
    return p


def main(argv: Optional[Sequence[str]] = None) -> int:
    _utf8_output()
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "func", None):
        parser.print_help()
        return 1
    try:
        return args.func(args)
    except ERRORS as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 2
    except RuntimeError as exc:
        print("error: %s" % exc, file=sys.stderr)
        return 3


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
