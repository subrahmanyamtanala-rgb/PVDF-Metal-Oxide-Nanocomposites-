"""pvdf_nano: physics models for PVDF / metal-oxide piezoelectric nanocomposites.

Modules
-------
constants   Physical constants.
materials   Typical property values for PVDF phases and metal-oxide fillers.
phases      Crystalline-phase quantification (FTIR, DSC, XRD).
interface   Interphase volume, Debye length, surface-charge nucleation metrics.
dielectric  Effective-medium permittivity and Maxwell-Wagner-Sillars relaxation.
piezo       Piezoelectric coupling of 0-3 composites and figures of merit.
harvester   Electromechanical output of PVDF-based energy harvesters.

All functions use SI units unless a docstring states otherwise.
"""

from . import constants, materials, phases, interface, dielectric, piezo, harvester

__all__ = [
    "constants",
    "materials",
    "phases",
    "interface",
    "dielectric",
    "piezo",
    "harvester",
]

__version__ = "0.1.0"
