"""Crystalline-phase quantification of PVDF and PVDF nanocomposites.

PVDF crystallises in at least five polymorphs. Three matter for harvesting:

* alpha (form II): TGTG' chains, antiparallel dipoles, non-polar.
* beta  (form I):  all-trans TTTT chains, parallel dipoles, highest
  spontaneous polarisation; the phase that gives PVDF its piezoelectricity.
* gamma (form III): T3GT3G' chains, polar but weaker than beta.

Characteristic signatures used below:

======== ========================== ======================================
Phase    FTIR bands (cm^-1)         XRD 2-theta, Cu K-alpha (deg)
======== ========================== ======================================
alpha    614, 763, 795, 976         17.7 (100), 18.4 (020), 19.9 (110),
                                    26.6 (021)
beta     840, 1275 (beta only)      20.6-20.8 (110)/(200)
gamma    833-840, 1234 (gamma only) 18.5 (020), 19.2 (002), 20.0 (110)
======== ========================== ======================================
"""

import math

from .constants import DH_100_PVDF

# Absorption coefficients at 763 cm^-1 (alpha) and 840 cm^-1 (electroactive),
# cm^2/mol (Gregorio & Cestari, J. Polym. Sci. B 32, 859 (1994)).
K_ALPHA_763 = 6.1e4
K_EA_840 = 7.7e4


def electroactive_fraction(a_763: float, a_840: float) -> float:
    """Fraction of electroactive (beta + gamma) phase in the crystalline part.

    Lambert-Beer based relation
        F_EA = A_840 / ((K_840 / K_763) * A_763 + A_840)

    Parameters are baseline-corrected absorbances at 763 and 840 cm^-1. The
    840 cm^-1 band is shared by beta and gamma, so the result is F(beta+gamma);
    use :func:`split_beta_gamma` to separate them.
    """
    if a_763 < 0 or a_840 < 0:
        raise ValueError("absorbances must be non-negative")
    denom = (K_EA_840 / K_ALPHA_763) * a_763 + a_840
    if denom == 0:
        raise ValueError("both absorbances are zero")
    return a_840 / denom


def electroactive_fraction_uncertainty(a_763: float, a_840: float,
                                       sigma_763: float, sigma_840: float) -> float:
    """Standard uncertainty of F_EA from independent absorbance uncertainties.

    First-order propagation through Eq. F_EA = A_840 / (k A_763 + A_840),
    k = K_840 / K_763:

        dF/dA_840 =  k A_763 / D^2,   dF/dA_763 = -k A_840 / D^2,
        D = k A_763 + A_840.

    sigma_763 and sigma_840 are typically the standard deviations of
    replicate, baseline-corrected spectra.
    """
    k = K_EA_840 / K_ALPHA_763
    d = k * a_763 + a_840
    if d <= 0:
        raise ValueError("absorbances must not both be zero")
    return math.hypot(k * a_763 / d ** 2 * sigma_840, k * a_840 / d ** 2 * sigma_763)


def split_beta_gamma(f_ea: float, dh_1275: float, dh_1234: float):
    """Split the electroactive fraction into beta and gamma parts.

    Uses the peak-to-valley heights of the beta-only 1275 cm^-1 band and the
    gamma-only 1234 cm^-1 band (Cai et al., RSC Adv. 7, 15382 (2017)):

        F(beta)  = F_EA * dH_1275 / (dH_1275 + dH_1234)
        F(gamma) = F_EA * dH_1234 / (dH_1275 + dH_1234)

    Returns ``(f_beta, f_gamma)``.
    """
    if not 0 <= f_ea <= 1:
        raise ValueError("f_ea must be in [0, 1]")
    total = dh_1275 + dh_1234
    if total <= 0:
        raise ValueError("peak heights must sum to a positive value")
    return f_ea * dh_1275 / total, f_ea * dh_1234 / total


def dsc_crystallinity(dh_melt: float, filler_wt_fraction: float = 0.0,
                      dh_cold: float = 0.0, dh_100: float = DH_100_PVDF) -> float:
    """Degree of crystallinity of the PVDF matrix from DSC.

        X_c = (dH_m - dH_cc) / ((1 - w_f) * dH_100)

    ``dh_melt`` and ``dh_cold`` (cold crystallisation) are in J per gram of
    composite; dividing by (1 - w_f) normalises to the polymer mass.
    """
    if not 0 <= filler_wt_fraction < 1:
        raise ValueError("filler weight fraction must be in [0, 1)")
    return (dh_melt - dh_cold) / ((1 - filler_wt_fraction) * dh_100)


def total_polar_content(x_c: float, f_polar: float) -> float:
    """Polar crystalline content of the whole matrix, X_c * F(polar)."""
    return x_c * f_polar


def bragg_d_spacing(two_theta_deg: float, wavelength: float = 1.5406e-10) -> float:
    """Interplanar spacing (m) from Bragg's law, first order."""
    theta = math.radians(two_theta_deg / 2)
    return wavelength / (2 * math.sin(theta))


def scherrer_size(fwhm_deg: float, two_theta_deg: float,
                  wavelength: float = 1.5406e-10, k: float = 0.9) -> float:
    """Crystallite size (m) from the Scherrer equation, D = K lambda / (beta cos theta).

    ``fwhm_deg`` is the instrument-corrected peak width in degrees 2-theta.
    """
    beta = math.radians(fwhm_deg)
    theta = math.radians(two_theta_deg / 2)
    return k * wavelength / (beta * math.cos(theta))


def xrd_crystallinity(crystalline_areas, amorphous_area: float) -> float:
    """Crystallinity from deconvoluted XRD peak areas, sum(A_c) / (sum(A_c) + A_a)."""
    a_c = sum(crystalline_areas)
    total = a_c + amorphous_area
    if total <= 0:
        raise ValueError("total area must be positive")
    return a_c / total
