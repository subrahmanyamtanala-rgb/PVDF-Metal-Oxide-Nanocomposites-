"""Typical material properties.

Values are representative room-temperature, low-frequency numbers collected
from standard handbooks and reviews. Real values depend strongly on processing
(poling, stretching, crystallinity, particle size), so treat them as starting
points for modelling, not as specifications. Piezoelectric coefficients are
magnitudes; PVDF has a negative d33 (it contracts along the poling direction
under a field parallel to the polarisation).
"""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Material:
    name: str
    eps_r: float  # relative permittivity
    d33: float  # |d33|, C/N
    youngs_modulus: float  # Pa
    density: float  # kg/m^3
    conductivity: float = 1e-14  # S/m
    crystal_structure: str = ""
    piezoelectric: bool = True
    band_gap_ev: Optional[float] = None


PC_PER_N = 1e-12  # 1 pC/N in C/N

# --- Polymer matrices ------------------------------------------------------

PVDF_ALPHA = Material(
    "PVDF (alpha, unpoled)", eps_r=10.0, d33=0.0, youngs_modulus=2.5e9,
    density=1780, conductivity=1e-13, crystal_structure="monoclinic TGTG'",
    piezoelectric=False,
)
PVDF_BETA = Material(
    "PVDF (beta, poled film)", eps_r=12.0, d33=28 * PC_PER_N,
    youngs_modulus=2.5e9, density=1780, conductivity=1e-13,
    crystal_structure="orthorhombic all-trans TTTT",
)
P_VDF_TRFE = Material(
    "P(VDF-TrFE) 70/30 (poled)", eps_r=10.0, d33=30 * PC_PER_N,
    youngs_modulus=1.5e9, density=1880, conductivity=1e-13,
    crystal_structure="orthorhombic all-trans",
)

# --- Metal-oxide fillers ---------------------------------------------------

ZNO = Material(
    "ZnO", eps_r=8.8, d33=12 * PC_PER_N, youngs_modulus=140e9, density=5610,
    conductivity=1e-6, crystal_structure="wurtzite (P6_3mc)", band_gap_ev=3.37,
)
BATIO3 = Material(
    "BaTiO3", eps_r=1700.0, d33=190 * PC_PER_N, youngs_modulus=110e9,
    density=6020, conductivity=1e-10, crystal_structure="tetragonal perovskite",
    band_gap_ev=3.2,
)
PZT = Material(
    "PZT-5H", eps_r=3400.0, d33=590 * PC_PER_N, youngs_modulus=60e9,
    density=7500, conductivity=1e-10, crystal_structure="perovskite (MPB)",
)
KNN = Material(
    "(K,Na)NbO3", eps_r=400.0, d33=110 * PC_PER_N, youngs_modulus=100e9,
    density=4500, conductivity=1e-10, crystal_structure="orthorhombic perovskite",
)
TIO2 = Material(
    "TiO2 (anatase)", eps_r=40.0, d33=0.0, youngs_modulus=230e9, density=3900,
    conductivity=1e-10, crystal_structure="tetragonal anatase",
    piezoelectric=False, band_gap_ev=3.2,
)
FE3O4 = Material(
    "Fe3O4", eps_r=20.0, d33=0.0, youngs_modulus=175e9, density=5170,
    conductivity=1e4, crystal_structure="inverse spinel", piezoelectric=False,
)
AL2O3 = Material(
    "Al2O3", eps_r=9.8, d33=0.0, youngs_modulus=370e9, density=3950,
    conductivity=1e-14, crystal_structure="corundum", piezoelectric=False,
    band_gap_ev=8.8,
)
SIO2 = Material(
    "SiO2 (amorphous)", eps_r=3.9, d33=0.0, youngs_modulus=70e9, density=2200,
    conductivity=1e-15, crystal_structure="amorphous", piezoelectric=False,
    band_gap_ev=9.0,
)

FILLERS = {m.name: m for m in (ZNO, BATIO3, PZT, KNN, TIO2, FE3O4, AL2O3, SIO2)}
MATRICES = {m.name: m for m in (PVDF_ALPHA, PVDF_BETA, P_VDF_TRFE)}
