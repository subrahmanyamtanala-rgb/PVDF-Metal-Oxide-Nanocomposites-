# Response to Reviewer — Round 3

We thank the reviewer for the third report and the recommendation of minor revision. Every requested change has been made except one item that needs the authors' decision (R22, below). Section and figure numbers refer to the revised `paper/manuscript.pdf`.

## Essential revisions (reviewer §27)

**1. "Four elements" → "five elements".** Corrected in the Introduction.

**2. "Confirms" → "supports".**
- The benchmark is now described as supporting a saturating form *for the specific electrospun PVDF/ZnO system examined*. This wording appears in the abstract, the Section 8 text, Table 9 and the conclusions.
- No "confirms" remains in the manuscript. The Morris comparison now says "agrees with" the one-at-a-time ranking.

**3. Electrostrictive scaling as an assumption.** Section 5.3 now reads: "the phenomenological relation d33 = 2Qε₀ε₃₃P_s holds. We adopt it as a *scaling assumption*, not as an established law for nanoscale BaTiO₃". The abstract says "under the adopted electrostrictive scaling assumption", and the conclusions say "as the adopted electrostrictive scaling assumes".

**4. Variable-X_c case.** Added as new Section 8.1, Fig. 10.
- *Method:* X_c(φ) and F(φ) are each fitted to the published PVDF/ZnO data with the saturating law. X_c rises from 0.42 to 0.66–0.75 and F from 0.83 to 0.92–0.95.
- *Comparison:* case B (fitted X_c) against case A (fixed X_c = 0.50, nominal law), both as ratios to the unfilled film:

| Quantity | Case A (fixed X_c) | Case B (fitted X_c) |
|---|---|---|
| \|d33\| at 1 vol% | 1.19× | 1.78–1.97× |
| d33·g33 at 1 vol% | 1.41× | 3.2–3.9× |
| \|d33\| at 10 vol% | 1.51× | 1.62–1.79× |

- *Interpretation:* The fixed-X_c results are therefore **conservative** for this system. Crystallinity changes contribute at least as much as polar-fraction changes, and they appear at much lower loading.
- *Caveat:* Applied to the measured saturated values, the linear scaling gives a matrix |d33| of 41–46 pC/N, above typical poled-PVDF values. We therefore report ratios, call case B an upper estimate of the relative gain, and identify this as a limit of the linear matrix scaling.
- The abstract and conclusions summarise the result.

**5. "Model-predicted" qualifiers.**
- Section 7.7 now opens "Under the assumed mechanical and electrical boundary conditions…" before quoting device numbers.
- The paired comparison reads "the model predicts a higher P_max for PVDF/ZnO".
- The optimum size, d33, figure of merit, power and the ZnO–BaTiO₃ ranking all carry "model-predicted" or equivalent qualifiers.

**6. Author placeholders.**
- Reference [11], Furukawa (1989): completed as IEEE Trans. Electr. Insul. 24 (1989) 375–394, and the red note is removed.
- License, release tag and archive DOI: **still open** — see "Outstanding items".

**7. Preprint replaced.** Reference [16] now cites the peer-reviewed article: Bagla et al., *J. Phys. Chem. C* 129 (2025) 5808–5820, doi:10.1021/acs.jpcc.4c07913. The values used were read from the open preprint version, and the text says so.

## Other comments

**R17 (positioning).** The abstract now states that "the framework is thus checked at the level of individual assumptions, with full experimental validation proposed". We avoided calling it "validated", in line with the reviewer's caution against artificial claims of validation.

**R19 (FTIR reproducibility).** Step 1 of the validation protocol now asks authors to:
- report spectral resolution (e.g. 4 cm⁻¹), the number of co-added scans and the baseline algorithm;
- use replicate spectra from at least three films per loading;
- propagate the replicate scatter into F_EA;
- fit X_c(φ) as well as F_polar(φ).

A new function, `phases.electroactive_fraction_uncertainty`, performs first-order uncertainty propagation and is tested against finite differences.

**R28 (title).** Adopted the reviewer's preferred title: *Physics-Based Modeling of Interfacial Polar-Phase Nucleation and Piezoelectric Coupling in PVDF–Metal Oxide Nanocomposites for Energy Harvesting*.

**R25 (full tensor homogenization).** As the reviewer suggests, this is left for future work and is acknowledged in Section 5.3.

## Code

pvdf_nano is now at version 0.4.0, with 48 unit tests. `examples/reproduce_paper.py` regenerates Fig. 10 and its data (`paper/results/fig10_variable_xc.csv`) along with all other results.

## Outstanding items for the authors

1. **License, release tag/commit and archive DOI.** This is the only placeholder left in the manuscript (Code and data availability). Choosing a license and minting a DOI (e.g. Zenodo) are author decisions.
2. **Author names and affiliations, and target-journal formatting.**
