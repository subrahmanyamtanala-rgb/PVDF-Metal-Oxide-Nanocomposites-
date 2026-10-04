"""Polar-phase nucleation laws and fitting to measured phase fractions.

The composite d33 model in :mod:`pvdf_nano.piezo` needs the polar fraction
F_polar of the matrix as an input. Two parameterised laws are provided:

* :func:`saturating_law` gives F_polar as a function of filler volume
  fraction phi. It is the shape most often reported for oxide-filled PVDF.
* :func:`interphase_law` gives F_polar as a function of interphase volume
  fraction phi_i. It ties nucleation to particle size through
  :func:`pvdf_nano.interface.interphase_fraction`.

Both are phenomenological. :func:`fit_saturating_law` fits the first to
measured (phi, F_polar) pairs, for example from FTIR via
:func:`pvdf_nano.phases.electroactive_fraction`, so the model can be driven
by data instead of assumed parameters.
"""

import numpy as np

# Named scenarios used for sensitivity analysis: (f0, delta_f, phi_sat).
SCENARIOS = {
    "none": (0.50, 0.00, 0.03),
    "weak": (0.50, 0.10, 0.03),
    "nominal": (0.50, 0.35, 0.03),
    "strong": (0.50, 0.45, 0.015),
}


def saturating_law(phi, f0=0.50, delta_f=0.35, phi_sat=0.03):
    """F_polar(phi) = f0 + delta_f * (1 - exp(-phi / phi_sat)).

    f0 is the polar fraction of the unfilled matrix, delta_f the maximum
    gain and phi_sat the loading at which 63 % of the gain is reached.
    """
    phi = np.asarray(phi, dtype=float)
    return f0 + delta_f * (1 - np.exp(-phi / phi_sat))


def interphase_law(phi_i, f0=0.50, delta_f=0.35, phi_i_sat=0.05):
    """F_polar(phi_i) = f0 + delta_f * (1 - exp(-phi_i / phi_i_sat)).

    Assumes the polar-phase gain is controlled by the volume of polymer
    within the interphase rather than by filler volume as such.
    """
    phi_i = np.asarray(phi_i, dtype=float)
    return f0 + delta_f * (1 - np.exp(-phi_i / phi_i_sat))


def fit_saturating_law(phi, f_polar, phi_sat_grid=None):
    """Least-squares fit of :func:`saturating_law` to measured data.

    For fixed phi_sat the model is linear in (f0, delta_f), so those are
    solved exactly and phi_sat is found by a 1-D grid search. No SciPy needed.

    Returns a dict with f0, delta_f, phi_sat and rmse.
    """
    phi = np.asarray(phi, dtype=float)
    f = np.asarray(f_polar, dtype=float)
    if phi.shape != f.shape or phi.size < 3:
        raise ValueError("need at least three (phi, f_polar) pairs of equal length")
    if phi_sat_grid is None:
        phi_sat_grid = np.geomspace(1e-3, 0.5, 400)

    best = None
    for s in phi_sat_grid:
        basis = np.column_stack([np.ones_like(phi), 1 - np.exp(-phi / s)])
        coef, *_ = np.linalg.lstsq(basis, f, rcond=None)
        rmse = float(np.sqrt(np.mean((basis @ coef - f) ** 2)))
        if best is None or rmse < best["rmse"]:
            best = {"f0": float(coef[0]), "delta_f": float(coef[1]),
                    "phi_sat": float(s), "rmse": rmse}
    return best
