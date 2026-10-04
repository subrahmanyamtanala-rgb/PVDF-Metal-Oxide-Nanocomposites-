import math

import numpy as np
import pytest

from pvdf_nano import dielectric, harvester, interface, materials, phases, piezo
from pvdf_nano.constants import EPS0, KB


# --- phases ------------------------------------------------------------------

def test_electroactive_fraction_limits():
    assert phases.electroactive_fraction(0.0, 1.0) == 1.0
    assert phases.electroactive_fraction(1.0, 0.0) == 0.0


def test_electroactive_fraction_equal_absorbance():
    # Equal absorbances -> 1 / (1 + K840/K763)
    expected = 1 / (1 + 7.7e4 / 6.1e4)
    assert phases.electroactive_fraction(0.5, 0.5) == pytest.approx(expected)


def test_split_beta_gamma_conserves_total():
    fb, fg = phases.split_beta_gamma(0.8, 3.0, 1.0)
    assert fb + fg == pytest.approx(0.8)
    assert fb == pytest.approx(0.6)


def test_dsc_crystallinity_normalises_filler():
    neat = phases.dsc_crystallinity(52.35)
    filled = phases.dsc_crystallinity(52.35 * 0.9, filler_wt_fraction=0.1)
    assert neat == pytest.approx(0.5)
    assert filled == pytest.approx(neat)


def test_bragg_beta_peak():
    # beta (110)/(200) at 20.6 deg -> d ~ 0.43 nm
    assert phases.bragg_d_spacing(20.6) == pytest.approx(4.31e-10, rel=0.01)


def test_scherrer_positive():
    assert 1e-9 < phases.scherrer_size(0.5, 20.6) < 1e-7


# --- interface -----------------------------------------------------------------

def test_interphase_dilute_limit():
    r, t, phi = 10e-9, 5e-9, 1e-4
    dilute = phi * ((1 + t / r) ** 3 - 1)
    assert interface.interphase_fraction(r, t, phi) == pytest.approx(dilute, rel=1e-3)


def test_interphase_bounded_by_matrix():
    assert interface.interphase_fraction(5e-9, 50e-9, 0.3) <= 0.7


def test_nano_has_more_interphase_than_micro():
    nano = interface.interphase_fraction(10e-9, 5e-9, 0.05)
    micro = interface.interphase_fraction(1e-6, 5e-9, 0.05)
    assert nano > 20 * micro


def test_interparticle_distance_touching():
    assert interface.interparticle_distance(1.0, math.pi / 6) == pytest.approx(0.0, abs=1e-12)


def test_debye_length_scaling():
    l1 = interface.debye_length(10, 1e22)
    l4 = interface.debye_length(10, 4e22)
    assert l1 / l4 == pytest.approx(2.0)


def test_langevin_limits():
    assert interface.langevin_orientation(1e-8) == pytest.approx(1e-8 / 3)
    assert interface.langevin_orientation(1e3) == pytest.approx(1 - 1e-3, rel=1e-6)


def test_alignment_ratio_definition():
    sigma, eps = 0.01, 10
    expected = interface.MU_CH2CF2 * sigma / (eps * EPS0) / (KB * 298.15)
    assert interface.alignment_ratio(sigma, eps) == pytest.approx(expected)


# --- dielectric ----------------------------------------------------------------

@pytest.mark.parametrize("model", [dielectric.maxwell_garnett, dielectric.bruggeman,
                                   dielectric.lichtenecker, dielectric.yamada])
def test_mixing_rule_endpoints(model):
    assert model(10.0, 1000.0, 0.0) == pytest.approx(10.0)
    assert model(10.0, 1000.0, 1.0) == pytest.approx(1000.0)


def test_mixing_rules_within_wiener_bounds():
    em, ef, phi = 10.0, 1700.0, 0.2
    lower = 1 / ((1 - phi) / em + phi / ef)
    upper = (1 - phi) * em + phi * ef
    for model in (dielectric.maxwell_garnett, dielectric.bruggeman,
                  dielectric.lichtenecker, dielectric.yamada):
        assert lower <= model(em, ef, phi) <= upper


def test_yamada_n3_equals_maxwell_garnett():
    assert dielectric.yamada(10, 100, 0.1, n=3) == pytest.approx(
        dielectric.maxwell_garnett(10, 100, 0.1))


def test_coated_sphere_limits():
    assert dielectric.coated_sphere(100, 20, 1.0, 0.0) == pytest.approx(100)
    assert dielectric.coated_sphere(20, 20, 1.0, 0.5) == pytest.approx(20)


def test_interphase_raises_permittivity():
    base = dielectric.maxwell_garnett(10, 1700, 0.05)
    with_i = dielectric.interphase_maxwell_garnett(10, 1700, 30, 0.05, 25e-9, 10e-9)
    assert with_i > base


def test_mws_bilayer_relaxation_peak():
    f = np.logspace(-4, 6, 2001)
    eps, tau = dielectric.mws_bilayer(10, 1e-12, 1e-6, 10, 1e-8, 1e-6, f)
    # Interfacial polarisation shows up as a step in eps' at f ~ 1 / (2 pi tau).
    eps_lo, eps_hi = eps.real[0], eps.real[-1]
    assert eps_lo > eps_hi
    f_half = f[np.argmin(np.abs(eps.real - (eps_lo + eps_hi) / 2))]
    assert f_half == pytest.approx(1 / (2 * np.pi * tau), rel=0.3)


def test_mws_sphere_matches_complex_maxwell_garnett():
    em, sm, ef, sf, phi = 10, 1e-12, 9, 1e-6, 0.1
    tau = dielectric.mws_sphere_relaxation_time(em, sm, ef, sf, phi)
    f = np.logspace(-2, 8, 4001)
    eps = dielectric.maxwell_garnett(dielectric.complex_permittivity(em, sm, f),
                                     dielectric.complex_permittivity(ef, sf, f), phi)
    # Remove the DC-conduction tail, then locate the relaxation loss peak.
    sigma_dc = (eps * 1j * 2 * np.pi * f * EPS0).real[0]
    loss = -eps.imag - sigma_dc / (2 * np.pi * f * EPS0)
    f_peak = f[np.argmax(loss)]
    assert f_peak == pytest.approx(1 / (2 * np.pi * tau), rel=0.1)


# --- piezo ---------------------------------------------------------------------

def test_matrix_d33_reference():
    assert piezo.matrix_d33(piezo.XC_REF, piezo.F_POLAR_REF) == pytest.approx(
        piezo.D33_PVDF_REF)


def test_local_field_small_for_high_k_filler():
    assert piezo.local_field_coefficient(12, 1700) < 0.03
    assert piezo.local_field_coefficient(12, 12) == pytest.approx(1.0)


def test_local_stress_limits():
    assert piezo.local_stress_coefficient(1, 1) == pytest.approx(1.0)
    assert piezo.local_stress_coefficient(1, 1e9) == pytest.approx(2.5, rel=1e-6)


def test_antiparallel_poling_adds():
    m, f = materials.PVDF_BETA, materials.BATIO3
    args = (0.2, 0.5, 0.85, m.eps_r, f.eps_r, m.youngs_modulus, f.youngs_modulus, f.d33)
    par = piezo.composite_d33(*args, filler_poling="parallel")
    anti = piezo.composite_d33(*args, filler_poling="antiparallel")
    none = piezo.composite_d33(*args, filler_poling="none")
    assert abs(anti) > abs(none) > abs(par)


def test_fom_and_g33():
    d, er = -28e-12, 12
    assert piezo.harvesting_fom(d, er) == pytest.approx(d * piezo.g33(d, er))
    assert piezo.g33(d, er) == pytest.approx(-0.2635, rel=0.01)


# --- harvester -----------------------------------------------------------------

def test_matched_load_maximises_power():
    q0, cap, f = 1e-8, 1e-9, 10.0
    r = np.logspace(3, 11, 20001)
    p = harvester.load_power(q0, cap, f, r)
    assert r[np.argmax(p)] == pytest.approx(harvester.optimal_load(cap, f), rel=0.01)
    assert p.max() == pytest.approx(harvester.max_power(q0, cap, f), rel=1e-4)


def test_open_circuit_voltage_consistency():
    d33, er, area, t, force = 28e-12, 12, 1e-4, 50e-6, 10.0
    q = harvester.charge_33(d33, force)
    c = harvester.capacitance(er, area, t)
    assert harvester.open_circuit_voltage(q, c) == pytest.approx(
        harvester.open_circuit_voltage_33(d33, force / area, t, er))


def test_bridge_charging_saturates_at_voc():
    q0, cp, cs = 1e-8, 1e-9, 1e-7
    v = harvester.bridge_charging(q0, cp, cs, 5000, diode_drop=0.3)
    assert v[-1] == pytest.approx(q0 / cp - 0.6, rel=1e-3)
    assert np.all(np.diff(v) >= 0)


def test_power_density_requires_one_normaliser():
    with pytest.raises(ValueError):
        harvester.power_density(1.0)
    assert harvester.power_density(2.0, area=4.0) == 0.5
