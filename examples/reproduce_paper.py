"""Reproduce every figure and table number in the manuscript.

Run from the repository root:

    python examples/reproduce_paper.py

Outputs
-------
paper/figures/fig*.pdf, fig*.png   manuscript figures
paper/results/*.csv               data behind each figure and table
paper/results/summary.json        every number quoted in the text

All results are model predictions under the stated assumptions; none are
measurements.
"""

import csv
import json
import os
import sys

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pvdf_nano import (dielectric, harvester, interface, nucleation, piezo,  # noqa: E402
                       sensitivity)
from pvdf_nano.materials import BATIO3, PVDF_BETA, ZNO  # noqa: E402

FIG_DIR = os.path.join(ROOT, "paper", "figures")
RES_DIR = os.path.join(ROOT, "paper", "results")

# Palette: one accent, a neutral grey, shades of the accent for ordered series.
ACCENT, ORANGE, GREY = "#2a78d6", "#eb6834", "#8a8984"
INK, INK2, GRIDC, SHADE = "#0b0b0b", "#52514e", "#e4e3df", "#f2f1ee"

M = PVDF_BETA
X_C = 0.50
EM_LIMIT = 0.15  # loading above which dilute effective-medium results are extrapolations

# Reference device: 2 cm x 2 cm, 60 um film, 50 N sinusoidal force.
AREA, THICK, FORCE, FREQ = 4e-4, 60e-6, 50.0, 5.0

summary = {}


# --- helpers -------------------------------------------------------------------

def style(ax, title, xlabel, ylabel):
    ax.set_title(title, loc="left", fontsize=9.5, color=INK, fontweight="bold")
    ax.set_xlabel(xlabel, color=INK2, fontsize=8.5)
    ax.set_ylabel(ylabel, color=INK2, fontsize=8.5)
    ax.grid(True, color=GRIDC, linewidth=0.7)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=7.5)


def shade_extrapolation(ax, xmax_pct):
    ax.axvspan(EM_LIMIT * 100, xmax_pct, color=SHADE, zorder=0)
    ax.text(EM_LIMIT * 100 + 0.5, ax.get_ylim()[1] * 0.97, "dilute-model\nextrapolation",
            fontsize=6.5, color=INK2, va="top")


def save(fig, name):
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=200, facecolor="white")
    plt.close(fig)


def write_csv(name, header, rows):
    with open(os.path.join(RES_DIR, f"{name}.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def filler_args(f):
    return dict(eps_m=M.eps_r, eps_f=f.eps_r, c_m=M.youngs_modulus,
                c_f=f.youngs_modulus, d33_f=f.d33)


def film(d33, eps_r, freq=FREQ, tan_delta=0.0, r_leak=None):
    cap = harvester.capacitance(eps_r, AREA, THICK)
    q0 = harvester.charge_33(abs(d33), FORCE)
    return {
        "d33_pCN": d33 * 1e12,
        "eps_r": float(eps_r),
        "fom_e12": piezo.harvesting_fom(d33, eps_r) * 1e12,
        "voc_V": harvester.open_circuit_voltage(q0, cap),
        "r_opt_MOhm": float(harvester.lossy_optimal_load(cap, freq, tan_delta, r_leak)) / 1e6,
        "p_max_nW": float(harvester.lossy_max_power(q0, cap, freq, tan_delta, r_leak)) * 1e9,
    }


# --- Figure 1: interphase volume vs particle size -------------------------------

def fig_interphase():
    diam = np.geomspace(5e-9, 2e-6, 300)
    fig, ax = plt.subplots(figsize=(6.0, 3.3))
    rows = []
    for t, alpha in ((2e-9, 0.45), (5e-9, 1.0), (10e-9, 0.7)):
        frac = np.array([interface.interphase_fraction(d / 2, t, 0.05) for d in diam])
        ls = "-" if t == 5e-9 else "--"
        ax.plot(diam * 1e9, frac * 100, color=ACCENT, alpha=alpha, lw=2, ls=ls,
                label=f"t = {t * 1e9:.0f} nm")
        rows += [(d * 1e9, t * 1e9, f * 100) for d, f in zip(diam, frac)]
    ax.set_xscale("log")
    style(ax, "Interphase volume at 5 vol% filler", "Particle diameter (nm)",
          "Interphase (% of composite)")
    ax.legend(frameon=False, fontsize=7.5)
    save(fig, "fig1_interphase")
    write_csv("fig1_interphase", ["diameter_nm", "shell_nm", "interphase_pct"], rows)
    summary["interphase_pct_t5"] = {
        str(d): round(interface.interphase_fraction(d / 2 * 1e-9, 5e-9, 0.05) * 100, 2)
        for d in (5, 10, 20, 50, 100, 1000)}
    summary["spacing_20nm_5vol_nm"] = round(interface.interparticle_distance(20, 0.05), 1)


# --- Figure 2 / Table 4: permittivity models -------------------------------------

def fig_permittivity():
    phi = np.linspace(0, 0.30, 151)
    models = {
        "Maxwell-Garnett": dielectric.maxwell_garnett(M.eps_r, BATIO3.eps_r, phi),
        "Bruggeman": dielectric.bruggeman(M.eps_r, BATIO3.eps_r, phi),
        "Lichtenecker": dielectric.lichtenecker(M.eps_r, BATIO3.eps_r, phi),
        "MG + 10 nm interphase": dielectric.interphase_maxwell_garnett(
            M.eps_r, BATIO3.eps_r, 30.0, phi, 50e-9, 10e-9),
    }
    fig, ax = plt.subplots(figsize=(6.0, 3.3))
    colors = {"Maxwell-Garnett": ACCENT, "Bruggeman": GREY, "Lichtenecker": GREY,
              "MG + 10 nm interphase": ORANGE}
    styles = {"Bruggeman": "-", "Lichtenecker": ":"}
    for name, eps in models.items():
        ax.plot(phi * 100, eps, color=colors[name], ls=styles.get(name, "-"), lw=2, label=name)
    ax.set_xlim(0, 30)
    style(ax, "Permittivity of PVDF/BaTiO$_3$ (100 nm particles)", "BaTiO$_3$ loading (vol%)",
          "Relative permittivity")
    shade_extrapolation(ax, 30)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    save(fig, "fig2_permittivity")
    table = []
    for p in (0.05, 0.10, 0.20, 0.30):
        row = [p * 100] + [round(float(np.interp(p, phi, e)), 1) for e in models.values()]
        table.append(row)
    write_csv("table4_permittivity", ["vol_pct"] + list(models), table)
    summary["table4"] = table
    eps10_bt = dielectric.maxwell_garnett(M.eps_r, BATIO3.eps_r, 0.10)
    eps10_zn = dielectric.maxwell_garnett(M.eps_r, ZNO.eps_r, 0.10)
    summary["LE_dilute_BaTiO3"] = round(piezo.local_field_coefficient(M.eps_r, BATIO3.eps_r), 4)
    summary["LE_10vol_BaTiO3"] = round(piezo.local_field_coefficient(eps10_bt, BATIO3.eps_r), 4)
    summary["LE_10vol_ZnO"] = round(piezo.local_field_coefficient(eps10_zn, ZNO.eps_r), 3)
    summary["LT_BaTiO3"] = round(piezo.local_stress_coefficient(M.youngs_modulus, BATIO3.youngs_modulus), 3)
    summary["LT_ZnO"] = round(piezo.local_stress_coefficient(M.youngs_modulus, ZNO.youngs_modulus), 3)


# --- Figure 3: d33 vs loading and poling direction ----------------------------------

def fig_poling():
    phi = np.linspace(0, 0.30, 61)
    f_pol = nucleation.saturating_law(phi, *nucleation.SCENARIOS["nominal"])
    fig, ax = plt.subplots(figsize=(6.0, 3.3))
    rows = []
    neat = abs(piezo.matrix_d33(X_C, f_pol[0])) * 1e12
    ax.axhline(neat, color=GREY, ls="--", lw=1)
    ax.text(0.5, neat - 1.6, f"neat film {neat:.1f} pC/N", ha="left", fontsize=7, color=INK2)
    for mode, color, alpha, label in (("antiparallel", GREY, 1.0, "antiparallel"),
                                      ("none", GREY, 0.5, "filler unpoled"),
                                      ("parallel", ACCENT, 1.0, "co-poled (parallel)")):
        d = piezo.composite_d33(phi, X_C, f_pol, filler_poling=mode, **filler_args(BATIO3))
        ax.plot(phi * 100, np.abs(d) * 1e12, color=color, alpha=alpha, lw=2)
        ax.text(30.4, abs(d[-1]) * 1e12, f"{label}: {abs(d[-1]) * 1e12:.1f}",
                fontsize=7, color=INK, va="center")
        rows += [(p * 100, mode, abs(v) * 1e12) for p, v in zip(phi, d)]
        summary[f"d33_BaTiO3_30vol_{mode}"] = round(abs(d[-1]) * 1e12, 1)
        summary[f"d33_BaTiO3_10vol_{mode}"] = round(abs(float(np.interp(0.10, phi, d))) * 1e12, 1)
    ax.set_xlim(0, 30)
    ax.set_ylim(0, 30)
    style(ax, "|d$_{33}$| of PVDF/BaTiO$_3$ by filler poling direction",
          "BaTiO$_3$ loading (vol%)", "|d$_{33}$| (pC/N)")
    shade_extrapolation(ax, 30)
    fig.subplots_adjust(right=0.78)
    save(fig, "fig3_poling")
    write_csv("fig3_poling", ["vol_pct", "filler_poling", "abs_d33_pCN"], rows)
    summary["d33_neat_pCN"] = round(neat, 1)


# --- Figure 4: nucleation-law sensitivity -------------------------------------------

def fig_nucleation_sensitivity():
    phi = np.linspace(0, 0.20, 81)
    alphas = {"none": 0.3, "weak": 0.5, "nominal": 1.0, "strong": 0.75}
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 5.0), sharex=True)
    rows = []
    for col, f in enumerate((ZNO, BATIO3)):
        eps = dielectric.maxwell_garnett(M.eps_r, f.eps_r, phi)
        for name, (f0, df, s) in nucleation.SCENARIOS.items():
            fp = nucleation.saturating_law(phi, f0, df, s)
            d = piezo.composite_d33(phi, X_C, fp, filler_poling="none", **filler_args(f))
            fom = piezo.harvesting_fom(d, eps)
            ls = "--" if name == "strong" else "-"
            axes[0, col].plot(phi * 100, np.abs(d) * 1e12, color=ACCENT, alpha=alphas[name],
                              ls=ls, lw=1.8, label=name)
            axes[1, col].plot(phi * 100, fom * 1e12, color=ACCENT, alpha=alphas[name],
                              ls=ls, lw=1.8, label=name)
            rows += [(f.name, name, p * 100, abs(a) * 1e12, b * 1e12) for p, a, b in zip(phi, d, fom)]
            i10 = int(np.argmin(np.abs(phi - 0.10)))
            summary[f"sens_{f.name}_{name}_10vol"] = {
                "d33": round(abs(d[i10]) * 1e12, 1), "fom": round(fom[i10] * 1e12, 2)}
        style(axes[0, col], f"PVDF/{f.name}", "", "|d$_{33}$| (pC/N)")
        style(axes[1, col], "", "Filler loading (vol%)", "d$_{33}$g$_{33}$ (10$^{-12}$ m$^2$/N)")
        for r in range(2):
            axes[r, col].axvspan(EM_LIMIT * 100, 20, color=SHADE, zorder=0)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, frameon=False, fontsize=7.5,
               title="Polar-phase gain scenario", title_fontsize=7.5)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"fig4_nucleation_sensitivity.{ext}"), dpi=200,
                    facecolor="white")
    plt.close(fig)
    write_csv("fig4_nucleation_sensitivity",
              ["filler", "scenario", "vol_pct", "abs_d33_pCN", "fom_e12"], rows)


# --- Figure 5: particle-size study --------------------------------------------------

SIZE_BASE = dict(phi=0.05, shell=5e-9, f0=0.50, delta_f=0.35, phi_i_sat=0.05)


def size_point(f, diameter, eps_i=None, shell=None, delta_f=None, phi_i_sat=None,
               eps_m=None, eta=1.0, d_ref=piezo.D33_PVDF_REF, tan_delta=0.0, phi=None,
               eps_f=None, gamma_share=0.0, alpha_gamma=1.0, filler_poling="none",
               d_f=None):
    """One model evaluation for the size study, tornado and global sensitivity.

    The nucleated polar phase is split into beta and gamma (``gamma_share``);
    gamma is weighted by ``alpha_gamma``. The reference film is all-beta.
    """
    p = SIZE_BASE
    phi = p["phi"] if phi is None else phi
    shell = p["shell"] if shell is None else shell
    delta_f = p["delta_f"] if delta_f is None else delta_f
    phi_i_sat = p["phi_i_sat"] if phi_i_sat is None else phi_i_sat
    eps_m = M.eps_r if eps_m is None else eps_m
    eps_i = eps_m if eps_i is None else eps_i
    eps_f = f.eps_r if eps_f is None else eps_f
    d_f = f.d33 if d_f is None else d_f
    r = diameter / 2
    phi_i = interface.interphase_fraction(r, shell, phi)
    fp = float(nucleation.interphase_law(phi_i, p["f0"], delta_f, phi_i_sat))
    f_eff = piezo.effective_polar_fraction(*piezo.split_polar(fp, gamma_share), alpha_gamma)
    eps_c = float(dielectric.interphase_maxwell_garnett(eps_m, eps_f, eps_i, phi, r, shell))
    l_e = piezo.local_field_coefficient(eps_c, eps_f)
    l_t = piezo.local_stress_coefficient(M.youngs_modulus, f.youngs_modulus)
    d_m = piezo.matrix_d33(X_C, f_eff, d33_ref=d_ref, poling_efficiency=eta)
    sign = {"parallel": 1.0, "antiparallel": -1.0, "none": 0.0}[filler_poling]
    d = (1 - phi) * d_m + sign * phi * l_e * l_t * abs(d_f)
    out = film(d, eps_c, tan_delta=tan_delta)
    out.update({"diameter_nm": diameter * 1e9, "phi_i_pct": phi_i * 100, "f_polar": fp,
                "f_eff": f_eff, "L_E": l_e, "L_T": l_t,
                "valid": bool(dielectric.coated_model_valid(phi, r, shell))})
    return out


def fig_size():
    diam = np.geomspace(10e-9, 1e-6, 120)
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 5.0), sharex=True)
    keys = (("d33_pCN", "|d$_{33}$| (pC/N)"), ("eps_r", "Relative permittivity"),
            ("fom_e12", "d$_{33}$g$_{33}$ (10$^{-12}$ m$^2$/N)"), ("p_max_nW", "P$_{max}$ at 5 Hz (nW)"))
    rows = []
    for f, color, eps_i, ls in ((BATIO3, GREY, 30.0, "--"), (ZNO, ACCENT, None, "-")):
        pts = [size_point(f, d, eps_i=eps_i) for d in diam]
        for ax, (k, lab) in zip(axes.flat, keys):
            vals = np.array([abs(p[k]) for p in pts])
            ax.plot(diam * 1e9, vals, color=color, lw=2, ls=ls, label=f"PVDF/{f.name}")
            ax.set_xscale("log")
            style(ax, "", "Particle diameter (nm)" if k in ("fom_e12", "p_max_nW") else "", lab)
        rows += [(f.name, p["diameter_nm"], p["phi_i_pct"], p["f_polar"], abs(p["d33_pCN"]),
                  p["eps_r"], p["fom_e12"], p["p_max_nW"]) for p in pts]
        for dn in (10, 20, 50, 100, 1000):
            p = size_point(f, dn * 1e-9, eps_i=eps_i)
            summary[f"size_{f.name}_{dn}nm"] = {k: round(abs(float(v)), 3) for k, v in p.items()}
    axes[0, 0].legend(frameon=False, fontsize=7, title="|d$_{33}$| identical (filler unpoled)",
                      title_fontsize=6.5)
    fig.suptitle("5 vol% filler, 5 nm interphase, interphase-controlled nucleation",
                 fontsize=8.5, color=INK2, x=0.02, ha="left")
    save(fig, "fig5_particle_size")
    write_csv("fig5_particle_size", ["filler", "diameter_nm", "interphase_pct", "f_polar",
                                     "abs_d33_pCN", "eps_r", "fom_e12", "p_max_nW"], rows)


# --- Figure 6: frequency and dielectric loss ------------------------------------------

def nominal_films(phi_f=0.10):
    fp = nucleation.saturating_law(phi_f, *nucleation.SCENARIOS["nominal"])
    fp0 = nucleation.saturating_law(0.0, *nucleation.SCENARIOS["nominal"])
    return {
        "PVDF/ZnO 10 vol%": (piezo.composite_d33(phi_f, X_C, fp, filler_poling="none",
                                                  **filler_args(ZNO)),
                             dielectric.maxwell_garnett(M.eps_r, ZNO.eps_r, phi_f)),
        "PVDF/BaTiO$_3$ 10 vol%": (piezo.composite_d33(phi_f, X_C, fp, filler_poling="none",
                                                        **filler_args(BATIO3)),
                                   dielectric.maxwell_garnett(M.eps_r, BATIO3.eps_r, phi_f)),
        "Neat PVDF": (piezo.matrix_d33(X_C, fp0), M.eps_r),
    }


def fig_frequency():
    freq = np.geomspace(1, 100, 200)
    films = nominal_films()
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.0))
    rows = []
    for (name, (d33, eps)), color in zip(films.items(), (ACCENT, ORANGE, GREY)):
        cap = harvester.capacitance(eps, AREA, THICK)
        q0 = harvester.charge_33(abs(d33), FORCE)
        p0 = harvester.lossy_max_power(q0, cap, freq)
        tand = harvester.conduction_tan_delta(1e-10, eps, freq) + 0.02
        p1 = harvester.lossy_max_power(q0, cap, freq, tan_delta=tand)
        axes[0].loglog(freq, p0 * 1e9, color=color, lw=2, label=name)
        axes[0].loglog(freq, p1 * 1e9, color=color, lw=1.5, ls="--")
        rows += [(name, f, a * 1e9, b * 1e9) for f, a, b in zip(freq, p0, p1)]
        summary[f"freq_{name}"] = {
            "p1Hz_lossless_nW": round(float(p0[0] * 1e9), 2),
            "p100Hz_lossless_nW": round(float(p0[-1] * 1e9), 1),
            "p1Hz_lossy_nW": round(float(p1[0] * 1e9), 2),
            "retention_1Hz": round(float(p1[0] / p0[0]), 3),
            "retention_100Hz": round(float(p1[-1] / p0[-1]), 3)}
    style(axes[0], "Matched-load power (solid: lossless)", "Frequency (Hz)", "P$_{max}$ (nW)")
    axes[0].legend(frameon=False, fontsize=6.5)
    d33, eps = films["PVDF/ZnO 10 vol%"]
    cap = harvester.capacitance(eps, AREA, THICK)
    q0 = harvester.charge_33(abs(d33), FORCE)
    p0 = harvester.lossy_max_power(q0, cap, freq)
    for sigma, alpha in ((1e-12, 0.35), (1e-11, 0.65), (1e-10, 1.0)):
        tand = harvester.conduction_tan_delta(sigma, eps, freq)
        ratio = harvester.lossy_max_power(q0, cap, freq, tan_delta=tand) / p0
        axes[1].semilogx(freq, ratio * 100, color=ACCENT, alpha=alpha, lw=2,
                         label=f"σ = {sigma:.0e} S/m")
        summary[f"retention_ZnO_sigma{sigma:.0e}_1Hz"] = round(float(ratio[0]), 3)
    axes[1].set_ylim(0, 105)
    style(axes[1], "PVDF/ZnO: power retained", "Frequency (Hz)",
          "P$_{max}$ / lossless (%)")
    axes[1].legend(frameon=False, fontsize=6.5, loc="lower right")
    save(fig, "fig6_frequency_loss")
    write_csv("fig6_frequency", ["film", "freq_Hz", "pmax_lossless_nW",
                                 "pmax_sigma1e-10_tand0.02_nW"], rows)


# --- Table 5: device comparison ------------------------------------------------------

def table_devices():
    rows = []
    for name, (d33, eps) in nominal_films().items():
        r = film(d33, eps)
        rl = film(d33, eps, tan_delta=0.05)
        rows.append([name.replace("$_3$", "3"), round(r["d33_pCN"], 1), round(r["eps_r"], 1),
                     round(r["fom_e12"], 2), round(r["voc_V"], 2), round(r["r_opt_MOhm"], 1),
                     round(r["p_max_nW"], 1), round(rl["p_max_nW"], 1)])
    write_csv("table5_devices", ["film", "d33_pCN", "eps_r", "fom_e12", "voc_V",
                                 "r_opt_MOhm", "p_max_nW", "p_max_tand0.05_nW"], rows)
    summary["table5"] = rows


# --- Figure 7: one-at-a-time sensitivity (tornado) ------------------------------------

def fig_tornado():
    base = size_point(ZNO, 50e-9)["p_max_nW"]
    params = [
        ("Interphase thickness t: 2-10 nm", dict(shell=2e-9), dict(shell=10e-9)),
        ("Particle diameter: 100-20 nm", dict(diameter=100e-9), dict(diameter=20e-9)),
        ("Max. polar gain ΔF: 0.10-0.45", dict(delta_f=0.10), dict(delta_f=0.45)),
        ("Nucleation saturation φ_i,sat: 0.10-0.02", dict(phi_i_sat=0.10), dict(phi_i_sat=0.02)),
        ("Matrix permittivity: 14-10", dict(eps_m=14.0), dict(eps_m=10.0)),
        ("Interphase permittivity: 30-12", dict(eps_i=30.0), dict(eps_i=12.0)),
        ("Poling efficiency η: 0.7-1.0", dict(eta=0.7), dict(eta=1.0)),
        ("Reference d33: -20 to -35 pC/N", dict(d_ref=-20e-12), dict(d_ref=-35e-12)),
        ("Loss tangent: 0.10-0", dict(tan_delta=0.10), dict(tan_delta=0.0)),
        ("ZnO loading: 2-10 vol%", dict(phi=0.02), dict(phi=0.10)),
    ]
    out = []
    for label, lo, hi in params:
        kw_lo = {"diameter": 50e-9, **lo}
        kw_hi = {"diameter": 50e-9, **hi}
        p_lo = size_point(ZNO, kw_lo.pop("diameter"), **kw_lo)["p_max_nW"]
        p_hi = size_point(ZNO, kw_hi.pop("diameter"), **kw_hi)["p_max_nW"]
        out.append((label, (p_lo / base - 1) * 100, (p_hi / base - 1) * 100))
    out.sort(key=lambda r: abs(r[2] - r[1]))
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    y = np.arange(len(out))
    for i, (label, lo, hi) in enumerate(out):
        ax.barh(i, lo, color=GREY, height=0.6, label="first value of range" if i == 0 else None)
        ax.barh(i, hi, color=ACCENT, height=0.6, label="second value of range" if i == 0 else None)
    ax.set_yticks(y, [r[0] for r in out], fontsize=7)
    ax.axvline(0, color=INK2, lw=0.8)
    style(ax, f"Change in P$_{{max}}$ from base case ({base:.1f} nW)", "Change (%)", "")
    ax.tick_params(axis="y", labelsize=7)
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    save(fig, "fig7_tornado")
    write_csv("fig7_tornado", ["parameter", "low_pct", "high_pct"], out)
    summary["tornado_base_nW"] = round(base, 2)
    summary["tornado"] = [(l, round(a, 1), round(b, 1)) for l, a, b in out[::-1]]


# --- Table: beta/gamma weighting ---------------------------------------------------

def table_beta_gamma():
    """Effect of gamma share and gamma activity on the PVDF/ZnO 10 vol% film."""
    phi = 0.10
    fp = float(nucleation.saturating_law(phi, *nucleation.SCENARIOS["nominal"]))
    eps = dielectric.maxwell_garnett(M.eps_r, ZNO.eps_r, phi)
    rows = []
    for share in (0.0, 0.25, 0.5):
        for alpha in (1.0, 0.5, 0.25):
            f_eff = piezo.effective_polar_fraction(*piezo.split_polar(fp, share), alpha)
            d = piezo.composite_d33(phi, X_C, f_eff, filler_poling="none", **filler_args(ZNO))
            rows.append([share, alpha, round(f_eff, 3), round(abs(d) * 1e12, 1),
                         round(piezo.harvesting_fom(d, eps) * 1e12, 2)])
    write_csv("table_beta_gamma", ["gamma_share", "alpha_gamma", "F_eff", "abs_d33_pCN",
                                   "fom_e12"], rows)
    summary["table_beta_gamma"] = rows


# --- Table: BaTiO3 filler permittivity and the local-field argument ------------------

def table_batio3_eps():
    """Filler term at 10 vol% for bulk and nanocrystal BaTiO3 permittivities."""
    phi = 0.10
    fp = float(nucleation.saturating_law(phi, *nucleation.SCENARIOS["nominal"]))
    d_m = abs(piezo.matrix_d33(X_C, fp)) * 1e12
    l_t = piezo.local_stress_coefficient(M.youngs_modulus, BATIO3.youngs_modulus)
    rows = []
    for eps_f in (80.0, 200.0, 1700.0):
        eps_c = dielectric.maxwell_garnett(M.eps_r, eps_f, phi)
        l_e = piezo.local_field_coefficient(eps_c, eps_f)
        for label, d_f in (("fixed 190 pC/N", BATIO3.d33),
                           ("scaled with eps_f", piezo.electrostrictive_d33(eps_f))):
            term = phi * l_e * l_t * d_f * 1e12
            anti = (1 - phi) * d_m + term
            par = (1 - phi) * d_m - term
            rows.append([eps_f, label, round(d_f * 1e12, 1), round(eps_c, 1), round(l_e, 3),
                         round(term, 2), round(anti, 1), round(par, 1),
                         round(piezo.harvesting_fom(anti * 1e-12, eps_c) * 1e12, 2)])
    write_csv("table_batio3_eps", ["eps_f", "d_f_rule", "d_f_pCN", "eps_c", "L_E",
                                   "filler_term_pCN", "abs_d33_anti_pCN", "abs_d33_par_pCN",
                                   "fom_anti_e12"], rows)
    summary["table_batio3_eps"] = rows
    summary["matrix_term_10vol_pCN"] = round((1 - phi) * d_m, 1)
    summary["LE_dilute_eps80"] = round(piezo.local_field_coefficient(M.eps_r, 80.0), 3)
    summary["LE_dilute_eps200"] = round(piezo.local_field_coefficient(M.eps_r, 200.0), 3)


# --- Figure 8: Morris global screening -------------------------------------------------

GLOBAL_BOUNDS = [
    ("diameter", 20e-9, 100e-9, "log"),
    ("shell", 2e-9, 10e-9, "lin"),
    ("phi", 0.02, 0.06, "lin"),
    ("delta_f", 0.10, 0.45, "lin"),
    ("phi_i_sat", 0.02, 0.10, "log"),
    ("eps_m", 10.0, 14.0, "lin"),
    ("eps_i_ratio", 1.0, 2.5, "lin"),
    ("eta", 0.7, 1.0, "lin"),
    ("d_ref", 20e-12, 35e-12, "lin"),
    ("tan_delta", 0.0, 0.10, "lin"),
    ("gamma_share", 0.0, 0.5, "lin"),
    ("alpha_gamma", 0.25, 1.0, "lin"),
    ("eps_f", None, None, "log"),
]
EPS_F_RANGE = {"ZnO": (8.0, 11.0), "BaTiO3": (80.0, 1700.0)}
LABELS = {"diameter": "particle diameter", "shell": "interphase thickness t",
          "phi": "filler loading", "delta_f": "max. polar gain ΔF",
          "phi_i_sat": "nucleation saturation", "eps_m": "matrix permittivity",
          "eps_i_ratio": "interphase permittivity ε$_i$/ε$_m$", "eta": "poling efficiency η",
          "d_ref": "reference |d33|", "tan_delta": "loss tangent",
          "gamma_share": "γ share of polar phase", "alpha_gamma": "γ activity α$_γ$",
          "eps_f": "filler permittivity"}


def bounds_for(f):
    lo, hi = EPS_F_RANGE[f.name]
    return [b if b[0] != "eps_f" else ("eps_f", lo, hi, "log") for b in GLOBAL_BOUNDS]


def global_model(f, filler_poling="none", d_f_rule="scaled"):
    def model(p):
        d_f = None
        if f.name == "BaTiO3" and d_f_rule == "scaled":
            d_f = piezo.electrostrictive_d33(p["eps_f"])
        return size_point(f, p["diameter"], eps_i=p["eps_i_ratio"] * p["eps_m"],
                          shell=p["shell"], delta_f=p["delta_f"], phi_i_sat=p["phi_i_sat"],
                          eps_m=p["eps_m"], eta=p["eta"], d_ref=-p["d_ref"],
                          tan_delta=p["tan_delta"], phi=p["phi"], eps_f=p["eps_f"],
                          gamma_share=p["gamma_share"], alpha_gamma=p["alpha_gamma"],
                          filler_poling=filler_poling, d_f=d_f)["p_max_nW"]
    return model


def fig_morris():
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.8), sharey=True)
    rows = []
    results = {f.name: sensitivity.morris(global_model(f), bounds_for(f), r=60, seed=7)
               for f in (ZNO, BATIO3)}
    order = sorted(results["ZnO"], key=lambda n: results["ZnO"][n]["mu_star"])
    y = np.arange(len(order))
    for ax, f in zip(axes, (ZNO, BATIO3)):
        res = results[f.name]
        mu = [res[n]["mu_star"] for n in order]
        sig = [res[n]["sigma"] for n in order]
        ax.barh(y, mu, color=ACCENT, height=0.6, label="μ* (mean |effect|)")
        ax.plot(sig, y, "o", color=INK2, ms=4, label="σ (nonlinearity/interaction)")
        ax.set_yticks(y, [LABELS[n] for n in order])
        style(ax, f"PVDF/{f.name}", "Effect on P$_{max}$ (nW)", "")
        ax.tick_params(axis="y", labelsize=7)
        for n, v in res.items():
            rows.append((f.name, n, round(v["mu_star"], 3), round(v["mu"], 3), round(v["sigma"], 3)))
        summary[f"morris_{f.name}"] = [(n, round(res[n]["mu_star"], 2), round(res[n]["sigma"], 2))
                                       for n in sorted(res, key=lambda k: -res[k]["mu_star"])]
    axes[1].legend(frameon=False, fontsize=7, loc="lower right")
    save(fig, "fig8_morris")
    write_csv("fig8_morris", ["filler", "parameter", "mu_star_nW", "mu_nW", "sigma_nW"], rows)


# --- Paired Latin-hypercube comparison of ZnO and BaTiO3 -------------------------------

def lhs_comparison(n=2000):
    shared = [b for b in GLOBAL_BOUNDS if b[0] != "eps_f"]
    x = sensitivity.latin_hypercube(shared + [("u_eps_zn", 0.0, 1.0, "lin"),
                                              ("u_eps_bt", 0.0, 1.0, "lin")], n, seed=11)
    names = [b[0] for b in shared]
    out = {}
    for case, poling, rule in (("unpoled", "none", "scaled"),
                               ("antiparallel, d_f scaled", "antiparallel", "scaled"),
                               ("antiparallel, d_f fixed", "antiparallel", "fixed")):
        mz, mb = global_model(ZNO, poling), global_model(BATIO3, poling, rule)
        pz, pb = [], []
        for row in x:
            p = dict(zip(names, row[:len(names)]))
            ez = EPS_F_RANGE["ZnO"][0] + row[-2] * np.diff(EPS_F_RANGE["ZnO"])[0]
            lo, hi = EPS_F_RANGE["BaTiO3"]
            eb = float(np.exp(np.log(lo) + row[-1] * (np.log(hi) - np.log(lo))))
            pz.append(mz({**p, "eps_f": ez}))
            pb.append(mb({**p, "eps_f": eb}))
        pz, pb = np.array(pz), np.array(pb)
        ratio = pz / pb
        out[case] = {"p_zno_gt_bto": round(float(np.mean(pz > pb)), 3),
                     "ratio_p5": round(float(np.percentile(ratio, 5)), 2),
                     "ratio_median": round(float(np.median(ratio)), 2),
                     "ratio_p95": round(float(np.percentile(ratio, 95)), 2),
                     "zno_pmax_p5_p50_p95": [round(float(np.percentile(pz, q)), 1)
                                             for q in (5, 50, 95)],
                     "bto_pmax_p5_p50_p95": [round(float(np.percentile(pb, q)), 1)
                                             for q in (5, 50, 95)]}
    summary["lhs_comparison"] = out
    write_csv("lhs_comparison", ["case", "p_zno_gt_bto", "ratio_p5", "ratio_median",
                                 "ratio_p95"],
              [[k, v["p_zno_gt_bto"], v["ratio_p5"], v["ratio_median"], v["ratio_p95"]]
               for k, v in out.items()])


# --- Figure 9: literature benchmark of the nucleation law ----------------------------

# Bagla et al., arXiv:2502.16547 (2025), Table 1: electrospun PVDF with carbon-coated
# ZnO nanoparticles (NP) and nanorods (NR). wt% -> (X_c %, F_beta %) as reported.
BAGLA = {
    "nanoparticles": [(0.0, 41.9, 83.0), (0.5, 68.1, 91.0), (1.0, 64.9, 92.0), (2.0, 65.0, 95.0)],
    "nanorods": [(0.0, 41.9, 83.0), (0.5, 68.8, 93.0), (1.0, 72.1, 92.0), (2.0, 76.1, 92.0)],
}


def wt_to_vol(w, rho_f=ZNO.density, rho_m=M.density):
    return (w / rho_f) / (w / rho_f + (1 - w) / rho_m)


def fig_benchmark():
    fig, ax = plt.subplots(figsize=(6.0, 3.3))
    fits = {}
    rows = []
    for (name, data), color, marker in zip(BAGLA.items(), (ACCENT, GREY), ("o", "s")):
        phi = np.array([wt_to_vol(w / 100) for w, _, _ in data])
        polar = np.array([xc * fb / 1e4 for _, xc, fb in data])
        fit = nucleation.fit_saturating_law(phi, polar, np.geomspace(1e-5, 0.05, 600))
        fit["phi_sat_resolved"] = fit["phi_sat"] > 2e-5
        fits[name] = {k: (round(v, 5) if not isinstance(v, bool) else v) for k, v in fit.items()}
        grid = np.linspace(0, phi.max() * 1.1, 200)
        ax.plot(phi * 100, polar, marker, color=color, ms=6, label=f"measured, {name}")
        ax.plot(grid * 100, nucleation.saturating_law(grid, fit["f0"], fit["delta_f"],
                                                      fit["phi_sat"]),
                color=color, lw=1.5, label=f"fit, {name}")
        rows += [(name, p * 100, v) for p, v in zip(phi, polar)]
    style(ax, "Polar crystalline content X$_c$F vs ZnO loading (Bagla et al. 2025)",
          "ZnO loading (vol%)", "X$_c$·F$_{EA}$")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    save(fig, "fig9_benchmark")
    write_csv("fig9_benchmark", ["morphology", "vol_pct", "xc_times_f"], rows)
    summary["benchmark_fits"] = fits
    summary["benchmark_vol_pct_of_2wt"] = round(wt_to_vol(0.02) * 100, 3)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    os.makedirs(RES_DIR, exist_ok=True)
    fig_interphase()
    fig_permittivity()
    fig_poling()
    fig_nucleation_sensitivity()
    fig_size()
    fig_frequency()
    table_devices()
    fig_tornado()
    table_beta_gamma()
    table_batio3_eps()
    fig_morris()
    lhs_comparison()
    fig_benchmark()
    with open(os.path.join(RES_DIR, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
