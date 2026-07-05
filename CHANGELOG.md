# Changelog — THPE fit code

## v1.7.1 — 2026-07-05

### Performance (critical for usability)

- **~10⁴× speedup.** h_ent(z) does not depend on fit parameters: it is
  now tabulated once and interpolated in log-log space (capturing the
  ~z⁻³ divergence at z→0). The comoving distance is computed as a
  single cumulative-trapezoid integral over a fine grid instead of one
  independent integration per supernova (1588 per likelihood call).
  A full 32×5000 MCMC now takes ~8 minutes on a modest desktop CPU
  (previously estimated in weeks). Relative distance error < 2e-3
  within the prior-allowed region, negligible vs observational errors.

### Fixed

- **Prior physicality range.** Φ(z) > 0 is now enforced from z = 0.011
  (the nearest SN after the mask), not from z = 0.1. With the current
  h_ent normalisation, significantly negative γ makes E²(z) < 0 exactly
  where SNe exist; the old prior admitted these silently.
- **Model comparison reference point.** AIC/BIC are evaluated at the
  maximum-probability sample of the chain instead of the posterior
  median; in degenerate posteriors (Φ₀/α/β) the median can sit far
  from the peak and bias the comparison against the larger model.
- **Windows console encoding.** stdout/stderr forced to UTF-8 so Greek
  characters (χ, Λ, α) no longer crash redirected output on cp1252.

### First full-run results (2026-07-05, DESI DR1 + Pantheon+ full cov.)

- Sanity check: χ²(ΛCDM | BAO) = 20.59 / 12 points (χ²/n = 1.72).
- R̂ ≤ 1.013 on all parameters (32 walkers × 5000 steps).
- Φ₀ = 0.608 ± 0.05, α = 0.068, β = 0.039, γ pinned at 0 by the
  physicality prior + data (horizon-entropy term is not viable under
  the current normalisation — renormalise h_ent or fix γ = 0).
- ΔAIC = +10.5, ΔBIC = +26.7 in favour of ΛCDM: the current data show
  no statistical preference for THPE. Strong Φ₀–α–β degeneracy
  confirmed (Φ₀·(1+α+β) ≈ Ω_Λ).

## v1.7 — 2026-07-01

### Fixed (critical)

- **BAO data vector.** The v1.6 fallback table labelled as D_V/r_d values
  that in DESI 2024 (Table 1, arXiv:2404.03002) are D_M/r_d. Against pure
  ΛCDM (Planck 2018) this produced χ² ≈ 186 over 7 points, with a spurious
  +8.4σ pull at z = 2.33. Any fit run on that vector is invalid: the MCMC
  would absorb the artifact by inflating α, β, γ, producing false support
  for THPE. v1.7 embeds the real DESI DR1 data vector — 12 observables:
  five (D_M/r_d, D_H/r_d) pairs with per-tracer correlation coefficients
  (LRG1, LRG2, LRG3+ELG1, ELG2, Lya) plus D_V/r_d for BGS and QSO.
  Verified: χ²(ΛCDM) = 20.6 / 12 points (χ²/n = 1.72).

- **SNIa likelihood.** v1.6 compared `MU_SH0ES` (SH0ES-calibrated,
  H0 ≈ 73) against a model normalised with H0 = 67.4 (Planck), injecting
  ~0.17 mag of calibration bias into the χ². v1.7 marginalises the constant
  magnitude offset analytically:
  χ²_marg = ΔᵀC⁻¹Δ − (ΔᵀC⁻¹1)² / (1ᵀC⁻¹1).
  Verified: the likelihood is invariant under constant offsets to machine
  precision (Δ log L ≈ 1e-13 for a +0.17 mag shift), and the estimator
  recovers an injected offset exactly.

### Added

- Startup sanity check: the script aborts if χ²(ΛCDM)/n > 5 against the
  BAO vector, so a data/bookkeeping error can never again pass silently.
- Support for the full Pantheon+ STAT+SYS covariance matrix
  (`Pantheon+SH0ES_STAT+SYS.cov`), sliced to the z > 0.01 mask.
- Best-offset estimator `best_offset_SNIa()` for visualisation.
- Separate D_M/r_d and D_H/r_d panels in the results figure; the Hubble
  diagram displays the marginalised ΔM.

### Changed

- Removed the unverified DESI download URL; DESI DR1 values are embedded
  with an explicit instruction to verify them against arXiv:2404.03002.
- Requirements: `pandas` and `requests` documented (both were already used
  by v1.6 but missing from the README).

### Notes

- v1.6 is kept in the repository as a historical record. No result
  produced by v1.6 should be used or cited.
- The harmonic-mean evidence estimator remains approximate; for
  publication use nested sampling (MultiNest, PolyChord, dynesty).

## v1.6 — 2026-06-25

Initial public release. Superseded by v1.7 (see above).
