"""Design study: how particle size, loading and filler choice shape a PVDF PENG.

Run from the repository root:

    python examples/design_study.py

Writes ``examples/figures/design_study.png`` and prints a summary table.

The beta-nucleation law used here, F_polar(phi) = 0.50 + 0.35 (1 - exp(-phi / 0.03)),
is an illustrative assumption with the saturating shape commonly reported
for oxide-filled PVDF films. Replace it with measured FTIR data (see
pvdf_nano.phases) for a real sample.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pvdf_nano import dielectric, harvester, interface, piezo  # noqa: E402
from pvdf_nano.materials import BATIO3, PVDF_BETA, ZNO  # noqa: E402

# Reference categorical palette, slots 1-3 (validated for all-pairs use).
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

X_C = 0.50  # matrix crystallinity, held constant for clarity


def f_polar(phi):
    """Illustrative beta-nucleation law (see module docstring)."""
    return 0.50 + 0.35 * (1 - np.exp(-np.asarray(phi) / 0.03))


def style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", fontsize=10, color=INK, fontweight="bold")
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)


def device_output(d33, eps_r, area=4e-4, thickness=60e-6, force=50.0, freq=5.0):
    """Charge, V_oc and matched-load power for a compression-mode film."""
    cap = harvester.capacitance(eps_r, area, thickness)
    q0 = harvester.charge_33(abs(d33), force)
    return {
        "C_nF": cap * 1e9,
        "Voc_V": harvester.open_circuit_voltage(q0, cap),
        "R_opt_MOhm": harvester.optimal_load(cap, freq) / 1e6,
        "P_max_uW": harvester.max_power(q0, cap, freq) * 1e6,
        "q0": q0,
        "cap": cap,
    }


def main():
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    m = PVDF_BETA

    # (a) Interphase volume vs particle size --------------------------------
    ax = axes[0, 0]
    diam = np.logspace(np.log10(5e-9), np.log10(2e-6), 200)
    for t, c in ((2e-9, C1), (5e-9, C2), (10e-9, C3)):
        frac = [interface.interphase_fraction(d / 2, t, 0.05) for d in diam]
        ax.plot(diam * 1e9, frac, color=c, lw=2, label=f"t = {t * 1e9:.0f} nm")
    ax.set_xscale("log")
    style(ax, "a  Interphase volume at 5 vol% filler", "Particle diameter (nm)",
          "Interphase volume fraction")
    ax.legend(frameon=False, fontsize=8)

    # (b) Permittivity models, BaTiO3 in PVDF -------------------------------
    ax = axes[0, 1]
    phi = np.linspace(0, 0.30, 151)
    ax.plot(phi * 100, dielectric.maxwell_garnett(m.eps_r, BATIO3.eps_r, phi),
            color=C1, lw=2, label="Maxwell-Garnett")
    ax.plot(phi * 100, dielectric.bruggeman(m.eps_r, BATIO3.eps_r, phi),
            color=C2, lw=2, label="Bruggeman")
    ax.plot(phi * 100, dielectric.interphase_maxwell_garnett(
        m.eps_r, BATIO3.eps_r, 30.0, phi, 50e-9, 10e-9),
        color=C3, lw=2, label="MG + 10 nm interphase (eps_i = 30)")
    style(ax, "b  Permittivity of PVDF/BaTiO3 (100 nm)", "BaTiO3 loading (vol%)",
          "Relative permittivity")
    ax.legend(frameon=False, fontsize=8)

    # (c) d33 vs loading for BaTiO3 under different poling ------------------
    ax = axes[1, 0]
    phi = np.linspace(0, 0.30, 151)
    common = dict(eps_m=m.eps_r, eps_f=BATIO3.eps_r, c_m=m.youngs_modulus,
                  c_f=BATIO3.youngs_modulus, d33_f=BATIO3.d33)
    for mode, c, label in (("antiparallel", C1, "filler poled antiparallel"),
                           ("none", C2, "filler unpoled (nucleation only)"),
                           ("parallel", C3, "filler poled parallel")):
        d = piezo.composite_d33(phi, X_C, f_polar(phi), filler_poling=mode, **common)
        ax.plot(phi * 100, np.abs(d) * 1e12, color=c, lw=2, label=label)
    style(ax, "c  |d33| of PVDF/BaTiO3", "BaTiO3 loading (vol%)", "|d33| (pC/N)")
    ax.legend(frameon=False, fontsize=8)

    # (d) Power vs load resistance for three films --------------------------
    ax = axes[1, 1]
    phi_f = 0.10
    films = {
        "Neat PVDF (F_polar 0.50)": (piezo.matrix_d33(X_C, f_polar(0.0)), m.eps_r),
        "PVDF + 10 vol% ZnO": (
            piezo.composite_d33(phi_f, X_C, f_polar(phi_f), m.eps_r, ZNO.eps_r,
                                m.youngs_modulus, ZNO.youngs_modulus, ZNO.d33,
                                filler_poling="none"),
            dielectric.maxwell_garnett(m.eps_r, ZNO.eps_r, phi_f)),
        "PVDF + 10 vol% BaTiO3": (
            piezo.composite_d33(phi_f, X_C, f_polar(phi_f), filler_poling="none",
                                **common),
            dielectric.maxwell_garnett(m.eps_r, BATIO3.eps_r, phi_f)),
    }
    r = np.logspace(5, 10, 400)
    rows = []
    for (name, (d33, eps)), c in zip(films.items(), (C1, C2, C3)):
        out = device_output(d33, eps)
        ax.plot(r, harvester.load_power(out["q0"], out["cap"], 5.0, r) * 1e6,
                color=c, lw=2, label=name)
        rows.append((name, d33, eps, piezo.harvesting_fom(d33, eps), out))
    ax.set_xscale("log")
    style(ax, "d  Power into a load (50 N, 5 Hz, 2 x 2 cm, 60 um)",
          "Load resistance (ohm)", "Average power (uW)")
    ax.legend(frameon=False, fontsize=8)

    fig.tight_layout()
    out_dir = os.path.join(os.path.dirname(__file__), "figures")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "design_study.png")
    fig.savefig(path, dpi=160, facecolor="white")
    print(f"saved {path}\n")

    header = (f"{'Film':<26}{'d33 pC/N':>10}{'eps_r':>8}{'FoM 1e-12/Pa':>14}"
              f"{'Voc V':>8}{'R_opt MOhm':>12}{'Pmax uW':>9}")
    print(header)
    print("-" * len(header))
    for name, d33, eps, fom, out in rows:
        print(f"{name:<26}{d33 * 1e12:>10.1f}{eps:>8.1f}{fom * 1e12:>14.2f}"
              f"{out['Voc_V']:>8.2f}{out['R_opt_MOhm']:>12.1f}{out['P_max_uW']:>9.3f}")


if __name__ == "__main__":
    main()
