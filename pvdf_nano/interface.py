"""Interfacial physics of polymer / oxide nanocomposites.

Nanocomposite behaviour comes mostly from the *interphase*: a polymer shell
around each particle whose chain conformation, crystal phase, charge
distribution and mobility differ from the bulk. This module provides
geometric and electrostatic estimates for that region.

References
----------
* T. Tanaka et al., IEEE Trans. Dielectr. Electr. Insul. 12, 669 (2005):
  multi-core model (bonded, bound and loose layers plus a diffuse
  Gouy-Chapman double layer).
* T. J. Lewis, IEEE Trans. Dielectr. Electr. Insul. 11, 739 (2004):
  interfaces as electrical double layers.
"""

import math

from .constants import EPS0, KB, E_CHARGE, DEBYE

# Dipole moment of a -CH2-CF2- repeat unit in the all-trans conformation,
# about 2.1 D (7.0e-30 C*m).
MU_CH2CF2 = 2.1 * DEBYE


def specific_interface_area(radius: float, phi: float) -> float:
    """Interfacial area per unit composite volume (m^2/m^3) for spheres, 3 phi / r."""
    if radius <= 0:
        raise ValueError("radius must be positive")
    return 3 * phi / radius


def interphase_fraction(radius: float, thickness: float, phi: float) -> float:
    """Volume fraction of the composite occupied by interphase polymer.

    Each sphere of radius ``r`` is wrapped in a shell of thickness ``t``. In
    the dilute limit the shells occupy phi * ((1 + t/r)^3 - 1). Shell overlap
    at higher loading is handled with a Poisson (random-overlap)
    approximation:

        phi_i = (1 - phi) * (1 - exp(-phi * ((1 + t/r)^3 - 1) / (1 - phi)))

    which tends to the dilute expression for small phi and never exceeds the
    available matrix fraction (1 - phi).
    """
    if not 0 <= phi < 1:
        raise ValueError("phi must be in [0, 1)")
    if radius <= 0 or thickness < 0:
        raise ValueError("radius must be positive and thickness non-negative")
    shell = (1 + thickness / radius) ** 3 - 1
    return (1 - phi) * (1 - math.exp(-phi * shell / (1 - phi)))


def interparticle_distance(diameter: float, phi: float) -> float:
    """Mean surface-to-surface spacing (m) for spheres on a simple cubic lattice.

        s = d * ((pi / (6 phi))^(1/3) - 1)

    A negative value means the particles would touch (phi > pi/6).
    """
    if phi <= 0:
        raise ValueError("phi must be positive")
    return diameter * ((math.pi / (6 * phi)) ** (1 / 3) - 1)


def debye_length(eps_r: float, ion_density: float, temperature: float = 298.15,
                 valence: int = 1) -> float:
    """Debye screening length (m) of the diffuse double layer.

        lambda_D = sqrt(eps_r eps0 k_B T / (2 n z^2 e^2))

    ``ion_density`` is the number density (1/m^3) of each ion species in a
    symmetric z:z electrolyte, e.g. residual ionic impurities in the polymer.
    """
    if ion_density <= 0:
        raise ValueError("ion density must be positive")
    return math.sqrt(eps_r * EPS0 * KB * temperature
                     / (2 * ion_density * valence ** 2 * E_CHARGE ** 2))


def surface_field(surface_charge_density: float, eps_r: float) -> float:
    """Electric field (V/m) just outside a charged surface, sigma / (eps_r eps0)."""
    return surface_charge_density / (eps_r * EPS0)


def dipole_alignment_energy(surface_charge_density: float, eps_r: float,
                            dipole: float = MU_CH2CF2) -> float:
    """Energy (J) to align one -CH2-CF2- dipole in the field of a charged surface, mu E.

    A negatively charged oxide surface (e.g. -OH or -O^- terminations) attracts
    the CH2 side of the chain (H^delta+) and a positive surface attracts CF2
    (F^delta-). In both cases the dipoles are pinned perpendicular to the
    surface, favouring all-trans (beta) segments.
    """
    return dipole * surface_field(abs(surface_charge_density), eps_r)


def alignment_ratio(surface_charge_density: float, eps_r: float,
                    temperature: float = 298.15,
                    dipole: float = MU_CH2CF2) -> float:
    """Ratio of the dipole alignment energy to thermal energy, mu E / k_B T.

    Values above about 1 indicate that surface fields can bias conformational
    statistics near the particle against thermal disorder.
    """
    return dipole_alignment_energy(surface_charge_density, eps_r, dipole) / (KB * temperature)


def langevin_orientation(x: float) -> float:
    """Mean dipole orientation <cos theta> = coth(x) - 1/x for alignment ratio x."""
    if abs(x) < 1e-6:
        return x / 3
    return 1 / math.tanh(x) - 1 / x


def tanaka_layers(radius: float, bonded: float = 1e-9, bound: float = 5e-9,
                  loose: float = 10e-9):
    """Outer radii (m) of the bonded, bound and loose layers of the multi-core model.

    Default thicknesses are typical orders of magnitude: ~1 nm (chemically
    bonded), 2-9 nm (bound, ordered chains) and tens of nm (loose, altered
    mobility). Returns a dict of outer radii.
    """
    r1 = radius + bonded
    r2 = r1 + bound
    r3 = r2 + loose
    return {"core": radius, "bonded": r1, "bound": r2, "loose": r3}
