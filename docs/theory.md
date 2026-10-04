# Interfacial and Piezoelectric Physics of PVDF–Metal Oxide Nanocomposites for Energy Harvesting

This document reviews the physics behind PVDF / metal-oxide nanocomposite
energy harvesters. It covers the polymer's polymorphism, how oxide
interfaces nucleate polar phases, interfacial polarisation, how
piezoelectric coupling is transferred across the particle–polymer boundary,
and how material properties turn into device output. Each section names the
`pvdf_nano` functions that implement its equations.

---

## 1. Why PVDF, and why add metal oxides?

Poly(vinylidene fluoride), –(CH2–CF2)n–, combines a large repeat-unit
dipole (about 2.1 D, from the electronegativity difference between H and F)
with flexibility, chemical stability, low acoustic impedance and easy
processing. Its piezoelectric coefficients are modest (|d33| ≈ 20–35 pC/N),
but its low permittivity (ε_r ≈ 10–12) gives a high voltage coefficient
(|g33| ≈ 0.2–0.33 V·m/N), roughly ten times that of PZT. That makes it well
suited to low-frequency, large-strain mechanical sources such as human
motion, wind, flow and vibration.

The polar β-phase does not form spontaneously from the melt or from most
solutions, and a film's response is limited by (i) how much of it is
crystalline (X_c), (ii) what fraction of the crystals are polar (F(β)),
and (iii) how well those crystals are poled. Metal-oxide nanoparticles
address all three:

| Role | Mechanism | Typical fillers |
|---|---|---|
| Polar-phase nucleation | Surface charge / hydroxyl groups pin chain conformation | ZnO, TiO2, Fe3O4, NiO, Al2O3, SiO2, CuO |
| Intrinsic piezo contribution | Filler is itself piezoelectric | BaTiO3, ZnO, KNN, PZT, BiFeO3 |
| Permittivity / field tailoring | High-ε or semiconducting filler modifies the local field | BaTiO3, TiO2, doped ZnO |
| Stress concentration | Stiff particles concentrate strain in the soft matrix | All oxides |
| Multifunctionality | Magnetoelectric or photo-coupling | CoFe2O4, NiFe2O4, Fe3O4, ZnO |

---

## 2. Polymorphism of PVDF

| Phase | Chain conformation | Unit cell | Polarity | Main FTIR bands (cm⁻¹) | Main XRD peaks, 2θ Cu Kα (°) |
|---|---|---|---|---|---|
| α (II) | TGTG′ | Monoclinic, antiparallel chains | Non-polar | 614, 763, 795, 976 | 17.7, 18.4, 19.9, 26.6 |
| β (I) | All-trans TTTT | Orthorhombic, parallel dipoles | Strongly polar (~0.13 C/m² in the crystal) | 840, 1275 | 20.6–20.8 |
| γ (III) | T3GT3G′ | Monoclinic | Moderately polar | 833–840, 1234 | 18.5, 19.2, 20.0 |
| δ (IV) | TGTG′ (polar α) | Monoclinic | Polar | ≈ α | ≈ α |

The α-phase is the kinetically favoured product of melt crystallisation. β
is thermodynamically stable only under stress, high fields, specific
solvent and temperature windows, or in the presence of nucleating
interfaces. γ forms by slow high-temperature crystallisation and from some
filler surfaces.

### 2.1 Quantifying the phases (`pvdf_nano.phases`)

**FTIR (Lambert–Beer).** The electroactive fraction of the crystalline part is

$$F_{EA} = \frac{A_{840}}{(K_{840}/K_{763})\,A_{763} + A_{840}}, \qquad K_{763}=6.1\times10^4,\ K_{840}=7.7\times10^4\ \mathrm{cm^2\,mol^{-1}}$$

(`electroactive_fraction`). The 840 cm⁻¹ band belongs to both β and γ, so
reporting F_EA as "β-fraction" is a common error. Split it with the
peak-to-valley heights ΔH of the β-only 1275 cm⁻¹ and γ-only 1234 cm⁻¹ bands:

$$F(\beta) = F_{EA}\,\frac{\Delta H_{1275}}{\Delta H_{1275}+\Delta H_{1234}},\qquad F(\gamma) = F_{EA}\,\frac{\Delta H_{1234}}{\Delta H_{1275}+\Delta H_{1234}}$$

(`split_beta_gamma`).

**DSC.** Matrix crystallinity, normalised to polymer mass:

$$X_c = \frac{\Delta H_m - \Delta H_{cc}}{(1-w_f)\,\Delta H_{100}},\qquad \Delta H_{100}\approx 104.7\ \mathrm{J\,g^{-1}}$$

(`dsc_crystallinity`). The useful single number for harvesting is the polar
crystalline content X_c·F(β) (`total_polar_content`): a filler that raises
F(β) but disrupts crystallisation can leave it unchanged.

**XRD.** Bragg spacing (`bragg_d_spacing`), Scherrer crystallite size
D = Kλ/(β cos θ) (`scherrer_size`), and crystallinity from deconvoluted
peak areas (`xrd_crystallinity`).

---

## 3. Interfacial physics

### 3.1 Why "nano" matters: interphase volume

A shell of thickness t around a sphere of radius r occupies, in the dilute limit,

$$\phi_i \approx \phi\left[\left(1+\tfrac{t}{r}\right)^3 - 1\right]$$

At 5 vol% loading with a 5 nm interphase, 1 µm particles give about
0.15 % of the composite as interphase. 20 nm particles give about 11 %,
10 nm particles about 29 %, and 5 nm particles about 70 %, at which point
most of the matrix is interphase (Figure a in `examples/design_study.py`). `interface.interphase_fraction` adds a
random-overlap correction so φ_i never exceeds the available matrix 1 − φ.
The specific interfacial area 3φ/r (`specific_interface_area`) and the
surface-to-surface spacing d[(π/6φ)^{1/3} − 1] (`interparticle_distance`)
show the same scaling: at a few vol% of 20 nm particles, neighbouring
interphases overlap and the whole matrix becomes interfacial.

### 3.2 Structure of the interphase: the multi-core model

Tanaka's multi-core model describes the polymer around a particle as
concentric layers (`interface.tanaka_layers`):

1. **Bonded layer (~1 nm):** chains tied to the surface by covalent,
   hydrogen or ionic bonds (e.g. surface –OH···F–C, coupling agents).
2. **Bound layer (~2–9 nm):** chains strongly ordered and immobilised by the
   bonded layer. This is where conformational templating (Section 3.3) acts.
3. **Loose layer (tens of nm):** chains with altered mobility, free volume
   and crystallisation kinetics.
4. **Diffuse electrical double layer:** a Gouy–Chapman space-charge cloud
   superimposed on the above, with Debye length

$$\lambda_D = \sqrt{\frac{\varepsilon_r\varepsilon_0 k_B T}{2 n z^2 e^2}}$$

(`debye_length`). In a polymer with ppm-level ionic impurities λ_D ranges
from tens of nm to µm, so double layers of neighbouring particles overlap
at modest loadings. Lewis's picture of interfaces as electrical double
layers explains why nanofillers change charge trapping, conduction and
breakdown so strongly.

### 3.3 Mechanisms of polar-phase nucleation at oxide surfaces

The β- and γ-phases are stabilised near oxide surfaces by several
cooperating mechanisms:

**(a) Ion–dipole / surface-charge interaction.** Oxide surfaces in contact
with the processing solvent (DMF, DMAc, NMP) carry net charge set by their
point of zero charge and surface hydroxylation. A negatively charged
surface (–O⁻, –OH with δ⁻ oxygen) attracts the electropositive –CH2–
groups. A positive surface (e.g. under-coordinated Zn²⁺, Fe³⁺, Ti⁴⁺)
attracts the –CF2– groups. Either way the C–F/C–H dipoles near the surface
are pinned normal to it. Adjacent monomers then line up their dipoles,
which only the all-trans (β) or T3G (γ) sequences allow. The surface field
is E = σ/(ε_rε_0) (`surface_field`), and the ratio of dipole alignment
energy to thermal energy

$$x = \frac{\mu E}{k_B T}$$

(`alignment_ratio`) tells you whether this biasing beats thermal disorder.
For σ ≈ 0.01 C/m² and ε_r ≈ 10, x ≈ 0.2 per monomer. Because chain
segments respond cooperatively over several monomers, this is enough to
tip the conformational balance. The mean orientation follows the Langevin
function ⟨cos θ⟩ = coth x − 1/x (`langevin_orientation`).

**(b) Hydrogen bonding.** Surface –OH groups form C–F···H–O bonds,
anchoring CF2 groups. FTIR shifts of the 3400 cm⁻¹ O–H band and the
1170–1180 cm⁻¹ CF2 stretch are the usual evidence. This is why hydrated or
hydroxylated oxides (TiO2, Fe(OH)3-coated Fe3O4, hydrated ZnO) are
especially effective nucleators.

**(c) Epitaxial / lattice matching.** When a filler surface has a lattice
spacing close to the β-PVDF (110)/(200) spacing (~0.43 nm) or the chain
repeat (~0.256 nm), heterogeneous nucleation of β is favoured. Examples
include ZnO and some clays.

**(d) Confinement and stress.** Narrow inter-particle gaps (Section 3.1)
restrict lamellar growth and impose local strain, both of which favour
the extended all-trans conformation, much as mechanical stretching does.

**(e) Altered crystallisation kinetics.** Particles act as heterogeneous
nuclei, raising the crystallisation temperature and the nucleation density.
This produces smaller spherulites with more interfacial area per crystal.

Two consequences matter for design. Nucleation efficiency saturates (and
then falls) once particles agglomerate, because agglomeration removes
interfacial area and introduces voids. Surface functionalisation
(silanes, dopamine, PVP, fluorinated coupling agents) both disperses the
particles and sets the sign and density of surface charge, so it is often
the most important processing parameter.

### 3.4 Interfacial (Maxwell–Wagner–Sillars) polarisation (`pvdf_nano.dielectric`)

When two phases differ in ε/σ, charge accumulates at their boundary under a
field. For a two-layer capacitor this gives a Debye-like relaxation with

$$\tau_{MWS} = \varepsilon_0\,\frac{\varepsilon_1 d_2 + \varepsilon_2 d_1}{\sigma_1 d_2 + \sigma_2 d_1}$$

(`mws_bilayer`). For dilute spheres (from the pole of the complex
Maxwell–Garnett expression)

$$\tau_{MWS} = \varepsilon_0\,\frac{\varepsilon_f + 2\varepsilon_m - \phi(\varepsilon_f-\varepsilon_m)}{\sigma_f + 2\sigma_m - \phi(\sigma_f-\sigma_m)}$$

(`mws_sphere_relaxation_time`). Passing complex permittivities from
`complex_permittivity` to any mixing rule gives the full spectrum.

For harvesting, MWS polarisation cuts both ways:

* **Helpful during poling:** trapped interfacial charge and space charge can
  stabilise oriented dipoles and add a quasi-permanent electret-like
  polarisation, which raises the measured d33. Semiconducting fillers (ZnO,
  Fe3O4) also redistribute the poling field (Section 4.2).
* **Harmful in operation:** a large low-frequency ε' and loss raise
  capacitance and leakage. That lowers g33 and V_oc and lets generated
  charge leak away at the low frequencies of human motion (1–10 Hz). Free
  charge can also screen the piezoelectric polarisation. Check tan δ in the
  0.1–100 Hz range, not just at 1 kHz.

### 3.5 Effective permittivity

| Model | Function | Notes |
|---|---|---|
| Maxwell–Garnett | `maxwell_garnett` | Dilute, isolated spheres; lower bound for well-dispersed composites |
| Bruggeman | `bruggeman` | Symmetric; treats phases equally; rises faster at high φ |
| Lichtenecker | `lichtenecker` | Empirical log mixing; often fits PVDF/BaTiO3 data well |
| Yamada | `yamada` | Shape factor n (n = 3: spheres) |
| Three-phase interphase | `interphase_maxwell_garnett` | Coated-sphere inclusion; captures interphase-driven enhancement |
| Percolation | `percolation_permittivity` | Conducting fillers below φ_c |

The permittivities of BaTiO3 (≈1700) and PVDF (≈12) are very mismatched,
so the composite permittivity at 10–20 vol% rises only to ~15–30. The
field concentrates in the polymer and the filler sees very little of it.
This mismatch is behind most of the piezoelectric coupling problems in the
next section.

---

## 4. Piezoelectric physics of the composite (`pvdf_nano.piezo`)

### 4.1 Constitutive relations and sign

In strain–charge form (one-dimensional, poling axis 3):

$$D_3 = d_{33}\,T_3 + \varepsilon_0\varepsilon_{33}^T E_3, \qquad S_3 = s_{33}^E T_3 + d_{33} E_3$$

(`direct_effect`, `converse_strain`). PVDF has a **negative** d33: under
compression along the polarisation, the polymer's dimensional
(Poisson/thickness) effect dominates, because the soft amorphous regions
change thickness more than the stiff crystals change dipole moment. Oxide
ceramics and ZnO have positive d33. This sign difference is easy to forget
and matters a great deal (Section 4.3).

### 4.2 Local field and stress coupling

For a dilute sphere (Furukawa's 0–3 model):

$$L_E = \frac{3\varepsilon_c}{2\varepsilon_c + \varepsilon_f},\qquad L_T = \frac{5c_f}{3c_m + 2c_f}$$

(`local_field_coefficient`, `local_stress_coefficient`).

* **Field coupling is poor for high-ε fillers.** For BaTiO3 in PVDF,
  L_E ≈ 3·12/(24 + 1700) ≈ 0.02. A 50 MV/m poling field puts only ~1 MV/m
  across the particle, too little to switch it fully. This is why
  BaTiO3-filled PVDF often needs high-temperature, high-field or corona
  poling, and why low-ε piezoelectric fillers (ZnO, ε ≈ 9) couple almost
  perfectly (L_E ≈ 1).
* **Stress coupling is favourable.** Stiff oxides in a soft polymer
  concentrate stress, with L_T → 2.5 for rigid spheres.
* **Semiconducting fillers redistribute the field.** If σ_f ≫ σ_m, the
  particle is screened at DC, but it raises the field in the polymer gaps
  between particles. That helps pole the matrix and can trigger local
  breakdown if particles agglomerate.

### 4.3 Effective d33

`composite_d33` combines the two contributions:

$$d_{33}^{c} = (1-\phi)\,d_m(X_c, F_{polar}) + \phi\,L_E L_T\,d_f$$

with d_m scaled linearly by the polar crystalline content relative to a
reference poled film. Three regimes follow (Figure c):

1. **Low loading (≲ 5–10 vol%), any oxide:** the gain comes from β/γ
   nucleation, i.e. from raising X_c·F(polar). The intrinsic filler term is
   small.
2. **Higher loading:** dilution of the piezoelectric matrix, agglomeration
   and loss of crystallinity take over. d33 peaks and then falls, which
   matches the optimum loading of a few wt% reported in most experimental
   studies.
3. **Co-poled piezoelectric fillers partly cancel the matrix.** Poled in
   the same direction, positive-d33 ceramic and negative-d33 PVDF oppose
   each other. Poling the two phases antiparallel (e.g. pole the ceramic
   above the polymer's Curie/melting point, then pole the polymer in the
   reverse direction at lower temperature/field) makes them add. This was
   demonstrated for PZT/P(VDF-TrFE) composites.

### 4.4 Figures of merit

| Quantity | Expression | Function |
|---|---|---|
| Voltage coefficient | g33 = d33 / (ε0 ε_r) | `g33` |
| Transduction FoM (off-resonance harvesting) | d33·g33 = d33² / (ε0 ε_r) | `harvesting_fom` |
| Coupling factor | k33 = \|d33\| √(Y / ε0 ε_r) | `coupling_factor` |

d33·g33 is the energy density per unit stress² for a quasi-static
harvester. It rewards d33 but penalises permittivity, so a filler that
raises d33 by 30 % while doubling ε_r *lowers* the harvested energy. This is
why ZnO and nucleation-only fillers often beat BaTiO3 in device tests at
equal loading (Figure d and the summary table printed by the example).

---

## 5. From material to device (`pvdf_nano.harvester`)

Below its mechanical resonance a PVDF film behaves as a charge source Q(t)
in parallel with its capacitance C_p = ε0ε_rA/t (`capacitance`).

* **Compression (33) mode:** Q = d33F (`charge_33`),
  V_oc = g33σt (`open_circuit_voltage_33`).
* **Bending / stretching (31) mode:** Q = d31σ1A (`charge_31`). Most
  flexible harvesters work this way; σ1 comes from beam bending
  (σ1 = Y·z/ρ for a layer at distance z from the neutral axis with
  curvature 1/ρ).
* **Energy per half-cycle:** Q²/2C_p (`stored_energy`).
* **Resistive load, sinusoidal force:**

$$P(R) = \frac{(\omega Q_0)^2 R}{2\left[1+(\omega R C_p)^2\right]},\qquad R_{opt} = \frac{1}{\omega C_p},\qquad P_{max} = \frac{\omega Q_0^2}{4C_p}$$

  (`load_power`, `optimal_load`, `max_power`). At 1–10 Hz and nF
  capacitance, R_opt is tens of MΩ. That is why PVDF harvesters are
  usually characterised into 1–100 MΩ loads, and why oscilloscope probe
  impedance (1–10 MΩ) distorts V_oc measurements.
* **Rectification and storage:** with an ideal full-wave bridge each
  half-cycle adds ΔV = [2Q0 − 2C_p(V + 2V_d)]/(C_p + C_s) to the storage
  capacitor, saturating at V_oc − 2V_d (`bridge_charging`). Fillers that
  raise C_p lower this saturation voltage.

### 5.1 Reporting checklist

Device numbers in the literature vary over orders of magnitude, mostly
because test conditions differ. A credible report states:

* force or pressure amplitude, impact area, frequency, waveform and contact
  mode;
* film thickness, electrode area and total active volume;
* load resistance for each power value, and P normalised per area *and*
  volume (`power_density`);
* the polarity-reversal test (the signal inverts when the connections are
  reversed), to rule out triboelectric and electrostatic artefacts;
* d33 measured with a quasi-static (Berlincourt) meter, with clamping
  noted, alongside FTIR-derived F(β), F(γ) and DSC X_c.

---

## 6. Filler-specific notes

| Filler | ε_r | \|d33\| (pC/N) | Dominant effect in PVDF | Remarks |
|---|---|---|---|---|
| ZnO (wurtzite) | ~9 | ~12 | Nucleation + good field coupling | Semiconducting; native-defect conductivity can screen polarisation; doping (Li, Ag, Mg) used to compensate |
| BaTiO3 | ~1700 (bulk) | ~190 | Nucleation + ε increase | Poor L_E; tetragonality and d33 fall below ~50–100 nm (size effect) |
| KNN, NaNbO3 | 300–500 | 80–150 | Lead-free piezo contribution | Better ε match than BaTiO3 |
| PZT | ~3400 | ~590 | Strong intrinsic contribution if poled antiparallel | Lead-containing |
| TiO2 | 40–100 | 0 | Nucleation via –OH groups | Non-piezoelectric; photo-responsive |
| Fe3O4, NiFe2O4, CoFe2O4 | — | 0 | Strong surface-charge nucleation; magnetoelectric coupling | Conductive (Fe3O4): watch percolation and leakage |
| Al2O3, SiO2 | 4–10 | 0 | Nucleation, breakdown-strength gain | Enable higher poling fields |

---

## 7. Open problems

* **Separating mechanisms:** quantitatively attributing a d33 gain to
  nucleation, intrinsic filler response, space-charge electret effects or
  triboelectric artefacts.
* **Interphase metrology:** directly measuring interphase thickness,
  permittivity and polarisation (AFM-IR, PFM, nano-DMA, Kelvin probe).
* **Size effects in ferroelectric fillers:** the loss of tetragonality in
  sub-50 nm BaTiO3 offsets interfacial benefits.
* **Long-term stability:** depolarisation and relaxation of space charge
  under humidity, temperature and cyclic loading (10⁶–10⁸ cycles).
* **Multiscale modelling:** linking DFT/MD of surface–chain interactions to
  phase-field crystallisation and continuum electromechanics.

---

## Key references

1. A. J. Lovinger, "Ferroelectric polymers," *Science* 220, 1115 (1983).
2. R. Gregorio Jr., M. Cestari, *J. Polym. Sci. B: Polym. Phys.* 32, 859 (1994): FTIR quantification of β-PVDF.
3. P. Martins, A. C. Lopes, S. Lanceros-Mendez, "Electroactive phases of poly(vinylidene fluoride): determination, processing and applications," *Prog. Polym. Sci.* 39, 683 (2014).
4. X. Cai, T. Lei, D. Sun, L. Lin, "A critical analysis of the α, β and γ phases in poly(vinylidene fluoride) using FTIR," *RSC Adv.* 7, 15382 (2017).
5. T. Tanaka, M. Kozako, N. Fuse, Y. Ohki, "Proposal of a multi-core model for polymer nanocomposite dielectrics," *IEEE Trans. Dielectr. Electr. Insul.* 12, 669 (2005).
6. T. J. Lewis, "Interfaces are the dominant feature of dielectrics at nanometric dimensions," *IEEE Trans. Dielectr. Electr. Insul.* 11, 739 (2004).
7. T. Furukawa, "Piezoelectricity and pyroelectricity in polymers," *IEEE Trans. Electr. Insul.* 24, 375 (1989).
8. B. Ploss, B. Ploss, F. G. Shin, H. L. W. Chan, C. L. Choy, "Pyroelectric or piezoelectric compensated ferroelectric composites," *Appl. Phys. Lett.* 76, 2776 (2000).
9. S. Yamada, T. Ueda, K. Kitayama, "Piezoelectricity of a high-content lead zirconate titanate/polymer composite," *J. Appl. Phys.* 53, 4328 (1982).
10. I. Katsouras et al., "The negative piezoelectric effect of the ferroelectric polymer poly(vinylidene fluoride)," *Nat. Mater.* 15, 78 (2016).
11. Z. L. Wang, J. Song, "Piezoelectric nanogenerators based on zinc oxide nanowire arrays," *Science* 312, 242 (2006).

Verify bibliographic details against the original sources before citing
them in a manuscript.
