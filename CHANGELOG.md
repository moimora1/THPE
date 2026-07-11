# Changelog — THPE fit code

## v1.8 — 2026-07-10 — Campaña de confrontación con datos (DR1→DR2→CMB)

Esta versión cierra la campaña empírica completa de la THPE contra los
mejores datos geométricos públicos. Resumen del veredicto: con
BAO (DESI DR2) + SNe Ia (Pantheon+) + ancla de distancia del CMB
(Planck 2018), la evidencia bayesiana favorece **fuertemente** a ΛCDM
sobre todas las variantes de la THPE, y el término de entropía del
horizonte (γ) queda excluido. Las contribuciones informacionales
restantes se acotan por debajo del ~1% de la densidad de energía
oscura. La THPE no queda falsada (es compatible con los datos) pero sí
resulta innecesaria en su parametrización actual.

### Añadido

- **`THPE_fit_v18.py`** — soporte para DESI DR2 (13 observables con
  correlaciones) como dataset por defecto, seleccionable con
  `--dataset dr1|dr2`. Valores verificados contra la fuente primaria.
- **`THPE_dynesty_v1.py`** — evidencia bayesiana por muestreo anidado
  (dynesty). Sustituye a la media armónica (no publicable). Calcula
  ln B para las variantes THPE-4p, THPE-3p (γ=0) y las de un solo
  trazador (solo-α, solo-β) frente a ΛCDM, con control de anidamiento.
- **`THPE_dynesty_v2_cmb.py`** — añade el ancla de distancia del CMB
  (parámetros comprimidos R y ℓ_A de Planck 2018), autocalibrada al
  fondo del pipeline. Incluye un centinela que aborta si el χ² del
  control ΛCDM es incoherente (evita resultados inválidos por errores
  de transcripción o de modelado del fondo).

### Resultados de la campaña

| Datos | Variante | ln B vs ΛCDM | Veredicto |
|---|---|---|---|
| DR1 (BAO+SNe) | THPE-4p | ΔAIC=+10.5 | pro-ΛCDM |
| DR2 (BAO+SNe) | solo-α (f_SFR) | −4.94 ± 0.4 | moderado pro-ΛCDM |
| DR2 (BAO+SNe) | THPE-3p | −11.06 ± 0.3 | fuerte pro-ΛCDM |
| **DR2+CMB** | **solo-α** | **−7.43 ± 0.19** | **fuerte pro-ΛCDM** |
| **DR2+CMB** | **THPE-3p** | **−14.79 ± 0.22** | **fuerte pro-ΛCDM** |

Cotas finales (68%): α < 0.008, β < 0.07, |γ| < ~10⁻³. El mejor ajuste
THPE coincide exactamente con ΛCDM (Δχ² = 0; α=0 en el máximo).

### Nota metodológica (control de calidad)

Durante la campaña, los controles de coherencia detectaron y
descartaron automáticamente dos corridas inválidas: (i) un ln B = +290
espurio a favor de la THPE causado por un truncamiento de la integral
del horizonte del sonido (r_s = 139 en vez de 144.4 Mpc); (ii) un
desajuste común de modelado del fondo (~0.1%) que sesgaba el ancla.
Ambos fueron señalados por el centinela `χ²_CMB(ΛCDM)` antes de
producir resultados. Documentado como advertencia: en comparación
bayesiana de modelos, un control de coherencia sobre el modelo de
referencia es imprescindible.

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
  for THPE. v1.7 embeds the real DESI DR1 data vector — 12 observables.
  Verified: χ²(ΛCDM) = 20.6 / 12 points (χ²/n = 1.72).

- **SNIa likelihood.** v1.6 compared `MU_SH0ES` (SH0ES-calibrated,
  H0 ≈ 73) against a model normalised with H0 = 67.4 (Planck), injecting
  ~0.17 mag of calibration bias into the χ². v1.7 marginalises the
  constant magnitude offset analytically.

### Added

- Startup sanity check: aborts if χ²(ΛCDM)/n > 5 against the BAO vector.
- Support for the full Pantheon+ STAT+SYS covariance matrix.

### Notes

- v1.6 is kept in the repository as a historical record. No result
  produced by v1.6 should be used or cited.

## v1.6 — 2026-06-25

Initial public release. Superseded (see above).
