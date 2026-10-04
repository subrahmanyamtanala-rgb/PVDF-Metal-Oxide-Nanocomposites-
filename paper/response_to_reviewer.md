# Response to Reviewer

**Manuscript:** Interfacial Polar-Phase Nucleation and Piezoelectric Coupling in PVDF–Metal Oxide Nanocomposites for Energy Harvesting: A Physics-Based Modeling Framework (revised title)

We thank the reviewer for a detailed and constructive report. We agree with the main assessment: the work is a modeling framework, not an experimentally validated materials study. The revision says so throughout. Section numbers below refer to the revised manuscript (`paper/manuscript.pdf`).

## Summary of changes

- The paper is reframed as a first-order modeling framework. Every result is labelled as a model prediction, and the central claim is qualified.
- The nucleation law is no longer a single assumption. Four scenarios are reported, and a least-squares routine fits the law to measured FTIR data.
- Three new analyses are added: particle size (Fig. 5), frequency and dielectric loss (Fig. 6), and one-at-a-time parameter sensitivity (Fig. 7).
- The device model now includes dielectric loss and leakage.
- Every reference was rechecked, two were corrected, and three were added. The note that references were "cited from memory" is removed.
- Equations are typeset in LaTeX.
- One script regenerates every figure and quoted number.

## Point-by-point responses

**R2. Interphase thickness presented as physical.** Agreed. Section 3.1 now states that *t* is a model parameter, not a measured thickness. Section 3.2 describes Tanaka's model as a conceptual model of interaction zones, not a fixed-thickness shell. The Results text now says these are geometric consequences of the assumed shell, and Fig. 7 varies *t* from 2 to 10 nm.

**R3 / R25. Central claim too strong.** Agreed. The abstract and conclusions now use the reviewer's suggested formulation: *"within the parameter range considered, the model indicates that nanoparticle-induced changes in the PVDF polar-phase content can contribute more…"*. The abstract also states that the conclusion depends on the nucleation law.

**R4. Circularity of the nucleation law (major).** Agreed, and this was the most important point. The revision makes three changes:
- Section 5.4 states the circularity explicitly.
- Section 7.4 and Fig. 4 report four scenarios (none, weak, nominal, strong). With no polar-phase gain, the model predicts that fillers *lower* |d33| (16.5 → 14.8 pC/N at 10 vol%). With a weak gain, PVDF/BaTiO3 still falls below the unfilled film in d33·g33 (2.22 vs 2.55 × 10⁻¹² m²/N).
- `nucleation.fit_saturating_law` fits the law to measured (φ, F_polar) data, and Section 8 gives a validation protocol built on it.

**R5. "ZnO is superior" too broad.** Agreed. Section 7.7 calls it a controlled model comparison under assumed properties. It also lists the real differences in surface chemistry, morphology, agglomeration, conductivity and nucleation efficiency that the model holds fixed.

**R6. Which field.** New Section 5.2 distinguishes the applied, matrix-average, particle-internal, interfacial, poling and operating fields. It states that *L_E* refers only to the average field inside a dilute spherical particle.

**R7. Maxwell–Garnett applicability.** The reviewer's sentence is added to Section 4.1, together with the effects that become significant at 20–30 vol%. All loading plots shade the region above 15 vol% as "dilute-model extrapolation", and Table 4 footnotes the 20 and 30 vol% rows.

**R8. Stress transfer.** Section 5.2 now states that *L_T* assumes perfect bonding and linear elasticity and is an upper estimate. It lists debonding, voids, agglomeration, viscoelasticity and frequency dependence.

**R9. d33 model justification.** Section 5.3 is retitled "First-order phenomenological mixture model". It lists six explicit assumptions, which answer each of the reviewer's questions: volume weighting, linearity in φ, multiplicative coupling, orientation, interphase piezoelectricity and particle shape. It also states that d31, d15 and full tensor homogenization are outside the scope.

**R10. Ploss system.** Section 5.3 now states that Ploss et al. studied lead titanate/P(VDF-TrFE) at 27 vol%. It supports the sign argument but does not directly validate PVDF/BaTiO3 or PVDF/ZnO.

**R11. Device model generality.** Section 6 lists what the linear, off-resonance model omits: resonance, contact mechanics, hysteresis, parasitic capacitance and rectifier losses.

**R12. Predicted power.** Table 6 is titled "Model-predicted device performance under assumed parameters", and Section 7 opens by stating that all results are model predictions.

**R13. Experimental validation (major).** *Partly addressed; this point remains open.* We could not perform measurements. We also did not digitise published datasets: values reproduced from figures without the original data would add uncertainty we could not document. The revision provides:
- an explicit five-step validation protocol (Section 8), built on the fitting routine;
- a statement in the abstract that the framework is not yet validated.

**This is the main remaining gap before submission** (see "Outstanding items").

**R14. References.** Every reference was rechecked against bibliographic records, and the "from memory" sentence is removed. Corrections:
- The Yamada et al. (1982) first author is T. (Takeshi) Yamada, not S. Yamada.
- Ploss et al. used lead titanate in P(VDF-TrFE); the text now says so.
- Furukawa (1989) is cited by start page only, because the end page could not be confirmed.

**R15. Missing citations.**
- *Surface-charge mechanism:* now cited to Martins et al., CrystEngComm 14 (2012) 2807, for the negative-surface/CH₂ interaction. Sebastian et al., RSC Adv. 6 (2016) 113007, is cited for positive-surface/CF₂ nucleation and the role of particle shape. The text presents both as findings of these studies.
- *Lattice matching:* softened. Epitaxy is now described as proposed mainly for layered silicates and is not used quantitatively, and the numerical spacing claim is removed.
- *BaTiO3 tetragonality:* now cited to Hoshina et al., Appl. Phys. Lett. 93 (2008) 192914.

**R16. Debye length.** Section 3.2 now distinguishes mobile ions, trapped space charge, injected carriers, Maxwell–Wagner charge and fixed surface charge. It uses λ_D only as an order-of-magnitude indicator.

**R17. Dielectric loss.** Sections 4.2 and 6 add complex permittivity, tan δ, conduction loss tan δ_σ = σ/(ωε₀ε′), and a lossy matched-load model:

P_max = (ωQ₀)² / {4[√(G² + (ωC)²) + G]}

These are implemented and tested in `harvester.lossy_*`. Table 6 adds a tan δ = 0.05 column.

**R18. Frequency dependence.** New Fig. 6 shows P_max from 1 to 100 Hz with and without loss, and the power retained for σ = 10⁻¹², 10⁻¹¹ and 10⁻¹⁰ S/m. With σ = 10⁻¹⁰ S/m, 84% of the lossless power is retained at 1 Hz and 98% at 100 Hz.

**R19. Particle size.** New Fig. 5 couples nucleation to interphase volume and plots d33, εr, d33·g33 and P_max from 10 nm to 1 µm. A non-trivial result emerges:
- PVDF/ZnO improves monotonically as particles shrink.
- PVDF/BaTiO3 peaks near 17 nm, because a high-permittivity interphase raises the composite permittivity.

**R20. Reproducibility.**
- `examples/reproduce_paper.py` regenerates every figure and writes every quoted number to `paper/results/` (CSV files plus `summary.json`).
- The package is versioned 0.2.0, and the unit tests now number 40.
- `requirements.txt` and a README are included.
- A license and an archive DOI are still to be added (see below).

**R21. Raw LaTeX equations.** Fixed. The manuscript is now authored in LaTeX.

**R22. Title.** Adopted the reviewer's recommended title.

**R3/R24. Minor corrections.**
- Chemical formulas and symbols are used consistently.
- "Particle size" is defined as diameter.
- Symbols are defined at each equation.
- The text explains the choice of X_c = 0.50, the origin of d_ref and the role of η. Parameter ranges are explored in Fig. 7.
- "ZnO film" now reads "PVDF/ZnO composite".
- Predictive statements are prefixed by "the model predicts".

## Outstanding items for the authors

1. **Experimental or literature validation (R13).** Measured FTIR/DSC phase fractions, dielectric spectra and quasi-static d33 for at least one PVDF/oxide series are needed to run the Section 8 protocol. Once the data are available, fitting and re-running the figures takes one command.
2. **License and archive DOI (R20).** Choose a license (for example MIT or BSD-3), tag a release, and archive it on Zenodo to obtain a DOI.
3. **Author names and affiliations,** and formatting for the target journal.
