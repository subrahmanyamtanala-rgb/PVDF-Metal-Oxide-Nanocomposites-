"""Dielectric models for 0-3 polymer / oxide composites.

Notation: ``eps_m`` is the matrix permittivity, ``eps_f`` the filler
permittivity and ``phi`` the filler volume fraction. Functions accept real or
complex permittivities and NumPy arrays wherever the formula allows. Passing
complex permittivities from :func:`complex_permittivity` gives the
Maxwell-Wagner-Sillars (MWS) interfacial relaxation of the composite directly.
"""

import numpy as np

from .constants import EPS0


def complex_permittivity(eps_r, conductivity, frequency):
    """Complex relative permittivity eps' - j sigma / (omega eps0) of a lossy dielectric."""
    omega = 2 * np.pi * np.asarray(frequency, dtype=float)
    return eps_r - 1j * conductivity / (omega * EPS0)


def maxwell_garnett(eps_m, eps_f, phi):
    """Maxwell-Garnett effective permittivity for dilute, isolated spheres.

        eps = eps_m * (eps_f + 2 eps_m + 2 phi (eps_f - eps_m))
                    / (eps_f + 2 eps_m -   phi (eps_f - eps_m))

    Accurate for phi below roughly 0.1-0.2; it never predicts percolation.
    """
    num = eps_f + 2 * eps_m + 2 * phi * (eps_f - eps_m)
    den = eps_f + 2 * eps_m - phi * (eps_f - eps_m)
    return eps_m * num / den


def bruggeman(eps_m, eps_f, phi):
    """Symmetric Bruggeman effective-medium permittivity for spheres (real inputs).

    Solves phi (eps_f - eps)/(eps_f + 2 eps) + (1 - phi)(eps_m - eps)/(eps_m + 2 eps) = 0.
    """
    b = (3 * phi - 1) * eps_f + (2 - 3 * phi) * eps_m
    return (b + np.sqrt(b ** 2 + 8 * eps_m * eps_f)) / 4


def lichtenecker(eps_m, eps_f, phi):
    """Lichtenecker logarithmic mixing rule, log eps = phi log eps_f + (1 - phi) log eps_m."""
    return np.exp(phi * np.log(eps_f) + (1 - phi) * np.log(eps_m))


def yamada(eps_m, eps_f, phi, n=3.0):
    """Yamada model for ellipsoidal inclusions with shape parameter ``n``.

        eps = eps_m * (1 + n phi (eps_f - eps_m) / (n eps_m + (eps_f - eps_m)(1 - phi)))

    n = 3 reproduces Maxwell-Garnett for spheres; larger n describes
    particles elongated along the field.
    """
    d = eps_f - eps_m
    return eps_m * (1 + n * phi * d / (n * eps_m + d * (1 - phi)))


def coated_sphere(eps_core, eps_shell, radius, thickness):
    """Effective permittivity of a sphere coated with a concentric shell.

    Used to represent a particle plus its interphase layer as one equivalent
    inclusion. ``q = (r / (r + t))^3`` is the core fraction of the coated
    particle.
    """
    q = (radius / (radius + thickness)) ** 3
    d = eps_core - eps_shell
    s = eps_core + 2 * eps_shell
    return eps_shell * (s + 2 * q * d) / (s - q * d)


COATED_FRACTION_LIMIT = 0.5


def coated_volume_fraction(phi, radius, thickness):
    """Volume fraction of particles plus shells, phi * (1 + t/r)^3."""
    return phi * (1 + thickness / radius) ** 3


def coated_model_valid(phi, radius, thickness, limit=COATED_FRACTION_LIMIT):
    """True when the coated spheres fill at most ``limit`` of the composite.

    Above roughly 0.5 the shells overlap strongly and Maxwell-Garnett mixing of
    isolated coated spheres is no longer meaningful; calculations are stopped
    there rather than extrapolated.
    """
    return coated_volume_fraction(phi, radius, thickness) <= limit


def interphase_maxwell_garnett(eps_m, eps_f, eps_i, phi, radius, thickness):
    """Three-phase (matrix / interphase / filler) Maxwell-Garnett permittivity.

    Each particle is replaced by an equivalent coated sphere whose volume
    fraction is phi * (1 + t/r)^3, capped at 0.74 (close packing of the
    coated spheres). This captures the strong permittivity enhancement that
    an interphase with eps_i > eps_m produces at nanometre particle sizes.
    """
    eps_eq = coated_sphere(eps_f, eps_i, radius, thickness)
    phi_eq = np.minimum(phi * (1 + thickness / radius) ** 3, 0.74)
    return maxwell_garnett(eps_m, eps_eq, phi_eq)


def percolation_permittivity(eps_m, phi, phi_c, s=1.0):
    """Power-law permittivity near a conductive-filler percolation threshold.

        eps = eps_m * |phi_c - phi|^(-s),   phi < phi_c

    Relevant for conducting or semiconducting fillers (Fe3O4, reduced ZnO, CNT
    co-fillers). Diverges at phi_c; only meaningful below it.
    """
    phi = np.asarray(phi, dtype=float)
    if np.any(phi >= phi_c):
        raise ValueError("model only valid below the percolation threshold")
    return eps_m * np.abs(phi_c - phi) ** (-s)


def mws_bilayer(eps1, sigma1, d1, eps2, sigma2, d2, frequency):
    """Complex permittivity of a two-layer series capacitor (MWS bilayer).

    Returns (eps_eff, tau), with eps_eff the complex relative permittivity of
    the stack at ``frequency`` and tau the interfacial relaxation time

        tau = eps0 (eps1 d2 + eps2 d1) / (sigma1 d2 + sigma2 d1).
    """
    e1 = complex_permittivity(eps1, sigma1, frequency)
    e2 = complex_permittivity(eps2, sigma2, frequency)
    eps_eff = (d1 + d2) / (d1 / e1 + d2 / e2)
    tau = EPS0 * (eps1 * d2 + eps2 * d1) / (sigma1 * d2 + sigma2 * d1)
    return eps_eff, tau


def mws_sphere_relaxation_time(eps_m, sigma_m, eps_f, sigma_f, phi):
    """MWS relaxation time (s) for dilute spheres in a matrix.

        tau = eps0 (eps_f + 2 eps_m - phi (eps_f - eps_m))
                  / (sigma_f + 2 sigma_m - phi (sigma_f - sigma_m))

    Obtained from the pole of the complex Maxwell-Garnett expression.
    """
    num = eps_f + 2 * eps_m - phi * (eps_f - eps_m)
    den = sigma_f + 2 * sigma_m - phi * (sigma_f - sigma_m)
    return EPS0 * num / den


def loss_tangent(eps_complex):
    """tan(delta) = eps'' / eps' for eps = eps' - j eps''."""
    eps_complex = np.asarray(eps_complex)
    return -eps_complex.imag / eps_complex.real
