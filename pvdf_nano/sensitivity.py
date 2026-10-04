"""Global sensitivity tools: Latin hypercube sampling and Morris screening.

Parameters are described by a list of ``(name, low, high, scale)`` tuples,
where ``scale`` is ``"lin"`` or ``"log"``. Models take a dict of parameter
values and return a scalar. Only NumPy is required.
"""

import numpy as np


def _to_physical(u, bounds):
    """Map unit-cube coordinates u (n x k) to physical parameter values."""
    out = np.empty_like(u, dtype=float)
    for j, (_, lo, hi, scale) in enumerate(bounds):
        if scale == "log":
            out[:, j] = np.exp(np.log(lo) + u[:, j] * (np.log(hi) - np.log(lo)))
        else:
            out[:, j] = lo + u[:, j] * (hi - lo)
    return out


def latin_hypercube(bounds, n, seed=0):
    """n x k Latin hypercube sample in physical units (one stratum per row)."""
    rng = np.random.default_rng(seed)
    k = len(bounds)
    u = (rng.permuted(np.tile(np.arange(n), (k, 1)), axis=1).T
         + rng.random((n, k))) / n
    return _to_physical(u, bounds)


def morris(model, bounds, r=40, levels=4, seed=0):
    """Morris elementary-effects screening.

    Builds ``r`` random one-at-a-time trajectories on a ``levels``-level grid
    in the unit cube (step Delta = levels / (2 (levels - 1))). Elementary
    effects are computed in unit-cube coordinates, so they are comparable
    across parameters with different units.

    Returns a dict name -> {"mu_star": mean |EE|, "mu": mean EE,
    "sigma": std EE}. A large sigma relative to mu_star indicates
    nonlinearity or interaction with other parameters.
    """
    rng = np.random.default_rng(seed)
    k = len(bounds)
    delta = levels / (2 * (levels - 1))
    grid = np.arange(levels // 2) / (levels - 1)  # start points that allow +delta
    effects = np.zeros((r, k))
    names = [b[0] for b in bounds]
    for t in range(r):
        x = rng.choice(grid, size=k)
        order = rng.permutation(k)
        point = x.copy()
        y0 = model(dict(zip(names, _to_physical(point[None, :], bounds)[0])))
        for j in order:
            point[j] += delta
            y1 = model(dict(zip(names, _to_physical(point[None, :], bounds)[0])))
            effects[t, j] = (y1 - y0) / delta
            y0 = y1
    return {n: {"mu_star": float(np.mean(np.abs(effects[:, j]))),
                "mu": float(np.mean(effects[:, j])),
                "sigma": float(np.std(effects[:, j], ddof=1))}
            for j, n in enumerate(names)}
