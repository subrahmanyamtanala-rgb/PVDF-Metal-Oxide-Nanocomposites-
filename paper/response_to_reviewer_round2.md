# Response to Reviewer — Round 2

We thank the reviewer for the second report and for recognising the improvements. Every point has been addressed. Section, table and figure numbers refer to the revised `paper/manuscript.pdf`.

Two changes go beyond wording, and both made the conclusions more conditional:

1. **Benchmarking against published data revealed that the BaTiO₃ local-field argument rested on bulk permittivity.** Li et al. (Nanoscale 2024) report a permittivity of about 80 for BaTiO₃ nanocrystals of 100 nm or less, and 150–300 for those of 200 nm or more. With these values, *L_E* rises from 0.021 to 0.16–0.35. The central claim now holds only if the nanocrystal d33 falls in proportion to its permittivity, as electrostriction implies (d33 = 2Qε₀ε₃₃P_s). The abstract, Sections 5.3 and 7.2, Table 5 and the conclusions now state this condition explicitly.
2. **The β/γ weighting can remove the predicted gain entirely.** With half the polar phase as γ and α_γ = 0.25, the PVDF/ZnO figure of merit falls below that of the unfilled film (Table 7).

## Point-by-point responses

**R9. "Composite d33" wording.** Section 5.3 now introduces Eq. (14) as the "model-predicted effective coefficient", and the text uses "model-predicted" throughout the Results.

**R11. Frequency scaling.** We adopted the reviewer's sentence verbatim in Section 7.6. We also list the real-device effects it excludes: mechanical impedance, resonance, viscoelasticity, rate-dependent hysteresis and contact dynamics.

**R13 / R23. The 17 nm BaTiO₃ optimum.** The abstract, Section 7.5 and the conclusions now read: "For the assumed interphase permittivity ε_i = 30, the model predicts a figure-of-merit maximum near 17 nm … not an experimentally optimized size". Section 7.5 adds that with ε_i = ε_m the BaTiO₃ curve would rise monotonically, as ZnO's does.

**R14. β/γ weighting (essential).**
- *Model change:* Implemented F_eff = F(β) + α_γF(γ) (Eq. 4, `piezo.effective_polar_fraction`). We found no literature value for α_γ that we could verify, so α_γ is treated as uncertain.
- *Sensitivity:* Table 7 covers γ share 0–0.5 × α_γ 0.25–1, and both parameters are included in the Morris screening, where they rank 4th and 7th of 13 for ZnO.
- *Baseline:* The baseline assumes that the nucleated phase is β, as reported for ferrite- and ZnO-nucleated PVDF.

**R15. FTIR preprocessing.** Section 2 now specifies:
- baseline-corrected absorbance heights, not areas, for Eq. (1);
- peak-to-valley heights for the β/γ split;
- consistent baselines and normalization across samples;
- deconvolution of overlapping bands;
- that the preprocessing should be reported together with the fractions.

**R16. Repository verification and environment.**
- Version 0.3.0 has 47 unit tests (`python -m pytest`).
- The new tests cover β/γ weighting, the electrostrictive scaling, the coated-sphere validity criterion, the Latin-hypercube strata, and Morris screening on functions with known sensitivities (linear, inert, interacting).
- `environment-lock.txt` pins Python 3.11.15, NumPy 2.4.6, Matplotlib 3.11.2 and pytest 9.1.1, and `requirements.txt` records the tested versions.
- The full analysis runs in about 5 s.

License, release tag and DOI are still open; see "Outstanding items".

**R17. Validation section title.** Renamed to "Proposed experimental validation protocol" (Section 9).

**R18. Literature validation (strongly recommended).** New Section 8 adds a benchmark table (Table 9) and Fig. 9. All values were read from the full text of the original articles. They test four assumptions:
1. *BaTiO₃ nanocrystal permittivity* (Li et al. 2024): about 80–300, against the bulk 1700. This changed the conclusions; see change 1 above.
2. *Mixing rules at 21–41 vol%* (Li et al. 2024): Maxwell–Garnett gives unphysical (negative) inferred particle permittivity, and Bruggeman is preferred. This supports our restriction of Maxwell–Garnett to about 15 vol% or less.
3. *Dense nanocomposites* (Pylypchuk et al. 2024): Lichtenecker mixing with ceramic permittivity overpredicts 28 vol% BaTiO₃ composites.
4. *Nucleation law* (Bagla et al. 2025, Table 1): we fitted our saturating law to the published X_c·F data for PVDF/ZnO with `nucleation.fit_saturating_law`, at r.m.s. error 0.009–0.011.
   - The saturating form is confirmed.
   - Saturation occurs below about 0.2 vol%, roughly 35 times faster than our nominal φ_sat of 3 vol%.
   - Most of the measured gain comes from crystallinity, which the model holds fixed.
   - Caveats are stated in the text: electrospun fibres, carbon-coated filler, F_EA rather than F(β), and a preprint source.

We did not include published d33 or power-density values. They come from heterogeneous geometries and test conditions that the model cannot be matched to without fabricating comparability.

**R21. Global sensitivity.** New Section 7.9 and Fig. 8 add two analyses.
- *Morris screening* (13 parameters, 60 trajectories) for both fillers:
  - It confirms the top three parameters: reference |d33| (μ* = 7.3 nW), poling efficiency (5.0 nW) and maximum polar gain (4.1 nW).
  - It shows σ ≈ μ* for most parameters, indicating strong interactions; the text states that one-at-a-time ranges understate joint uncertainty.
- *Paired Latin-hypercube comparison* of ZnO and BaTiO₃ (2000 samples, shared parameters), probability that ZnO gives the higher P_max:

| Case | P(ZnO > BaTiO₃) |
|---|---|
| Fillers unpoled | 100% |
| Both poled antiparallel, BaTiO₃ d33 scaled with permittivity | 100% |
| BaTiO₃ nanocrystals assumed to keep bulk d33 | 57% |

**R22. Coated-sphere validity criterion.** This is now explicit (Section 4.1, `dielectric.coated_model_valid`). Calculations stop when φ(1 + t/r)³ > 0.5, which at 5 vol% and t = 5 nm means d ≥ 8.7 nm. The particle-size study therefore starts at 10 nm.

**R24. BaTiO₃ permittivity uncertainty.** Table 4 now gives both the bulk value and the nanocrystal range of 80–1700, citing Li et al. Table 5 shows how the filler term depends on ε_f, under fixed and under permittivity-scaled d_f. ε_f is also included in the global analysis.

**R25. Intrinsic vs effective measured d33.** Section 5.1 defines both terms and states that the model predicts an intrinsic effective coefficient. Step 4 of the validation protocol adds a short-circuited thermal-annealing control, given as an example, to estimate the space-charge contribution.

**R26. Reference audit.**
- DOIs were added for 12 of 17 references; each DOI was taken from the publisher, a repository or the article itself.
- Three references were added: Li 2024 (Nanoscale), Pylypchuk 2024 (arXiv) and Bagla 2025 (arXiv preprint).
- Furukawa (1989), Tanaka (2005) and Hoshina (2008): volume and start page are confirmed, but not every DOI could be retrieved because the publisher sites were not reachable from our environment. **Furukawa's end page and DOI remain unconfirmed** and are flagged in red in the manuscript for the authors to complete.

**R28. Suggested abstract.** We adopted the reviewer's abstract as the base. We added the BaTiO₃ nanocrystal condition, the β/γ result and the benchmark finding.

## Outstanding items for the authors

1. **License, release tag/commit and archive DOI.** These are marked in red in the Code and data availability statement. Choosing a license and minting a DOI (e.g. Zenodo) are author decisions.
2. **Furukawa (1989) page range and DOI.** Marked in red in the reference list.
3. **Measured nucleation data and d33 for the authors' own system.** The benchmark shows that the nominal nucleation law is too slow for at least one published PVDF/ZnO system, so the authors' own FTIR/DSC data should replace it before submission if available.
4. **Author names and affiliations, and journal formatting.**
