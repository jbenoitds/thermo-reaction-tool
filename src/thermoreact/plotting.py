"""Plots of reaction properties versus temperature (requires matplotlib)."""
from __future__ import annotations

from typing import Sequence


def plot_reaction(rxn, temperatures: Sequence[float], path: str, cp_mode: str = "full") -> str:
    """Save a two-panel figure: ΔrH°, ΔrG°, TΔrS° vs T, and the van 't Hoff plot ln K vs 1/T."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise RuntimeError("Plotting needs matplotlib: pip install 'thermoreact[plot]'") from exc

    Ts = sorted(temperatures)
    props = [rxn.properties(T, cp_mode) for T in Ts]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))
    ax1.plot(Ts, [p.dH / 1000 for p in props], label=r"$\Delta_rH^\circ$")
    ax1.plot(Ts, [p.dG / 1000 for p in props], label=r"$\Delta_rG^\circ$")
    ax1.plot(Ts, [p.T * p.dS / 1000 for p in props], "--", label=r"$T\Delta_rS^\circ$")
    ax1.axhline(0, color="0.6", lw=0.8)
    ax1.set_xlabel("T / K")
    ax1.set_ylabel("kJ mol$^{-1}$")
    ax1.legend(frameon=False)
    ax1.set_title("Standard reaction properties", fontsize=10)
    ax2.plot([1000.0 / T for T in Ts], [p.lnK for p in props], "o-", ms=3)
    ax2.axhline(0, color="0.6", lw=0.8)
    ax2.set_xlabel("1000 K / T")
    ax2.set_ylabel(r"$\ln K$")
    ax2.set_title("van 't Hoff plot (slope = $-\\Delta_rH^\\circ$/1000R)", fontsize=10)
    fig.suptitle(rxn.equation(), fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path
