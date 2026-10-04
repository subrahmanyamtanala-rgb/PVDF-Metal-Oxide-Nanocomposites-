"""pvdf_nano: physics models for PVDF / metal-oxide piezoelectric nanocomposites.

Modules
-------
constants   Physical constants.
materials   Typical property values for PVDF phases and metal-oxide fillers.
phases      Crystalline-phase quantification (FTIR, DSC, XRD).
interface   Interphase volume, Debye length, surface-charge nucleation metrics.
nucleation  Polar-phase nucleation laws and fitting to measured data.
dielectric  Effective-medium permittivity and Maxwell-Wagner-Sillars relaxation.
piezo       Piezoelectric coupling of 0-3 composites and figures of merit.
harvester   Electrical output of PVDF-based harvesters, with dielectric loss.
sensitivity Latin hypercube sampling and Morris global sensitivity screening.

All functions use SI units unless a docstring states otherwise.
"""

from . import (constants, materials, phases, interface, nucleation, dielectric, piezo,
               harvester, sensitivity)

__all__ = [
    "constants",
    "materials",
    "phases",
    "interface",
    "nucleation",
    "dielectric",
    "piezo",
    "harvester",
    "sensitivity",
]

__version__ = "0.3.0"
