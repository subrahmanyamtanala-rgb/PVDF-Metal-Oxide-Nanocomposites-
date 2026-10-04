"""Piezoelectric response of PVDF / metal-oxide 0-3 nanocomposites.

Sign convention
---------------
Coefficients here are *signed*. Poled PVDF has d33 < 0 (dimensional effect
of the polymer), while ceramic and wurtzite oxides poled in the same
direction have d33 > 0. If the filler and matrix are poled in the same
direction their contributions therefore partially cancel. This is why
``filler_poling`` is an explicit argument below.

Model
-----
The effective coefficient is the sum of a matrix term and a filler term:

    d33 = (1 - phi) * d_m(X_c, F_polar) + phi * L_E * L_T * d_f

* d_m scales linearly with the polar crystalline content X_c * F_polar,
  anchored to a reference poled film. Particles act on the matrix term
  indirectly, by changing X_c and F_polar (interfacial nucleation).
* L_E and L_T are the local field and local stress coefficients of the
  Furukawa dilute-sphere model (T. Furukawa, IEEE Trans. Electr. Insul. 24,
  375 (1989)):

      L_E = 3 eps_c / (2 eps_c + eps_f)
      L_T = 5 c_f  / (3 c_m  + 2 c_f)

  with eps_c the composite permittivity. A high-permittivity filler in a
  low-permittivity polymer sees only a small fraction of the applied field
  (L_E << 1), so it is hard to pole and contributes little directly.

This is a first-order, phenomenological model: it reproduces trends and
orders of magnitude, not the absolute response of a particular sample.
"""

import numpy as np

from .constants import EPS0
from .dielectric import maxwell_garnett

# Reference poled beta-PVDF film used to anchor the matrix term.
D33_PVDF_REF = -28e-12  # C/N
XC_REF = 0.50
F_POLAR_REF = 0.85


def local_field_coefficient(eps_c, eps_f):
    """Furukawa local field coefficient L_E = 3 eps_c / (2 eps_c + eps_f)."""
    return 3 * eps_c / (2 * eps_c + eps_f)


def local_stress_coefficient(c_m, c_f):
    """Furukawa local stress coefficient L_T = 5 c_f / (3 c_m + 2 c_f)."""
    return 5 * c_f / (3 * c_m + 2 * c_f)


def matrix_d33(x_c, f_polar, d33_ref=D33_PVDF_REF, x_c_ref=XC_REF,
               f_polar_ref=F_POLAR_REF, poling_efficiency=1.0):
    """Matrix d33 scaled by polar crystalline content relative to a reference film.

    ``poling_efficiency`` (0-1) is the fraction of polar crystallites whose
    dipoles are switched into the poling direction.
    """
    return d33_ref * poling_efficiency * (x_c * f_polar) / (x_c_ref * f_polar_ref)


def composite_d33(phi, x_c, f_polar, eps_m, eps_f, c_m, c_f, d33_f,
                  filler_poling="parallel", poling_efficiency=1.0,
                  filler_poling_efficiency=1.0):
    """Effective signed d33 (C/N) of a PVDF / oxide 0-3 composite.

    Parameters
    ----------
    phi : filler volume fraction.
    x_c, f_polar : matrix crystallinity and polar-phase fraction (may be
        functions of phi measured experimentally).
    eps_m, eps_f : relative permittivities of matrix and filler.
    c_m, c_f : elastic moduli (Pa) of matrix and filler.
    d33_f : magnitude of the filler d33 (C/N); 0 for non-piezoelectric oxides.
    filler_poling : "parallel" (filler and matrix polarisation aligned,
        contributions oppose because the signs differ), "antiparallel"
        (contributions add) or "none" (filler not poled).
    """
    eps_c = maxwell_garnett(eps_m, eps_f, phi)
    l_e = local_field_coefficient(eps_c, eps_f)
    l_t = local_stress_coefficient(c_m, c_f)
    d_m = matrix_d33(x_c, f_polar, poling_efficiency=poling_efficiency)

    sign = {"parallel": +1.0, "antiparallel": -1.0, "none": 0.0}[filler_poling]
    d_f = sign * abs(d33_f) * filler_poling_efficiency
    return (1 - phi) * d_m + phi * l_e * l_t * d_f


def g33(d33, eps_r):
    """Piezoelectric voltage coefficient g33 = d33 / (eps0 eps_r), V*m/N."""
    return d33 / (EPS0 * eps_r)


def harvesting_fom(d33, eps_r):
    """Transducer figure of merit FoM = d33 * g33 = d33^2 / (eps0 eps_r), m^2/N."""
    return d33 ** 2 / (EPS0 * eps_r)


def coupling_factor(d33, eps_r, youngs_modulus):
    """Electromechanical coupling k33 = |d33| sqrt(Y / (eps0 eps_r)) (thin-rod approximation)."""
    return np.abs(d33) * np.sqrt(youngs_modulus / (EPS0 * eps_r))


def direct_effect(stress, field, d33, eps_r):
    """Electric displacement D = d33 T + eps0 eps_r E (C/m^2), strain-charge form."""
    return d33 * stress + EPS0 * eps_r * field


def converse_strain(stress, field, d33, youngs_modulus):
    """Strain S = T / Y + d33 E, strain-charge form."""
    return stress / youngs_modulus + d33 * field
