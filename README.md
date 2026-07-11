# THPE — Holographic Theory of State Persistence
### Cosmological fit code · Código de ajuste cosmológico

**Authors / Autores:** Moisés Mora García, Jorge Ordóñez Mora (aeronautical engineer)
**Location:** Palma de Mallorca, Spain · 2026

This repository contains the statistical fitting code accompanying the paper
*"The Holographic Theory of State Persistence: An Informational Extension of
Holographic Dark Energy"* (Mora García & Ordóñez Mora, 2026).

> **Result of the empirical campaign (v1.8, July 2026).** Confronted with
> DESI DR2 (BAO) + Pantheon+ (SNe Ia) + CMB distance priors (Planck 2018),
> the Bayesian evidence favours ΛCDM **strongly** over every THPE variant
> (ln B = −7.4 for the best single-tracer model, −14.8 for the 3-parameter
> model). The horizon-entropy term (γ) is excluded; residual informational
> couplings are bounded below ~1% of the dark-energy density. The THPE is
> **not falsified** (it is compatible with the data) but is **unnecessary**
> in its present parametrisation. Full numbers and the correction history
> are in `CHANGELOG.md`. This is an honest negative result, published in
> full with its coherence controls and reproducible pipeline.
>
> **Resultado de la campaña empírica (v1.8, julio 2026).** Frente a
> DESI DR2 + Pantheon+ + priors de distancia del CMB, la evidencia bayesiana
> favorece **fuertemente** a ΛCDM sobre toda variante de la THPE
> (ln B = −7.4 el mejor modelo de un trazador; −14.8 el de 3 parámetros).
> El término de entropía del horizonte queda excluido; los acoplamientos
> informacionales restantes se acotan por debajo del ~1% de la densidad de
> energía oscura. La THPE **no queda falsada** (es compatible con los datos)
> pero sí resulta **innecesaria** en su parametrización actual. Cifras
> completas e historial de correcciones en `CHANGELOG.md`.
>
> **Scripts:** `THPE_fit_v18.py` (MCMC fit, DR1/DR2) · `THPE_dynesty_v1.py`
> (Bayesian evidence, BAO+SNe) · `THPE_dynesty_v2_cmb.py` (evidence with the
> CMB distance anchor). Earlier versions (`v16`, `v17`) are kept as a
> historical record; **no result from v1.6 should be used**.

---

## English

### What this code does

The THPE promotes the dark-energy density to a redshift-dependent
holographic-pressure function

```
Φ(z) = Φ₀ · [ 1 + α·f_SFR(z) + β·g_struct(z) + γ·h_ent(z) ]
```

recovering ΛCDM exactly when α = β = γ = 0. `THPE_fit_v18.py` fits it to
DESI baryon acoustic oscillations (DR1 or DR2) and the Pantheon+ Type Ia
supernova sample via MCMC; `THPE_dynesty_v1.py` and `THPE_dynesty_v2_cmb.py`
compute the Bayesian evidence with nested sampling (the latter adding the
Planck CMB distance anchor).

### Requirements

```bash
pip install numpy scipy matplotlib pandas emcee requests corner dynesty
```

Python 3.10+ with NumPy 2.x (the code uses `np.trapezoid`; on NumPy 1.x
replace with `np.trapz` or upgrade).

### Data

**DESI** — DR1 (12 observables) and DR2 (13 observables) BAO vectors, with
per-tracer D_M/r_d–D_H/r_d correlations, are embedded in `THPE_fit_v18.py`.
Select with `--dataset dr1|dr2` (default: dr2). Values transcribed from the
official tables; verify before a definitive fit.

**Pantheon+** — the scripts download the public files from the
`PantheonPlusSH0ES/DataRelease` repository (distance moduli + full STAT+SYS
covariance). If the download fails, place them in the script folder and use
`--no-download`. The full covariance is required for a publishable result.

### How to run

```bash
python THPE_fit_v18.py --quick               # quick MCMC test
python THPE_fit_v18.py                        # full MCMC fit (DR2 by default)
python THPE_dynesty_v1.py                     # Bayesian evidence, BAO+SNe
python THPE_dynesty_v2_cmb.py                 # Bayesian evidence + CMB anchor
```

On startup each script runs a consistency check: χ²(ΛCDM) against the data
vector must be reasonable. If it is anomalously large the script stops — that
signals a data or bookkeeping error, not physics. The CMB script additionally
aborts if the ΛCDM control χ² is incoherent (protecting against transcription
or background-modelling errors).

### Output

- Best-fit parameters (Φ₀, α, β, γ) with credible intervals
- Results and corner figures (`THPE_v18_*.png`)
- Model comparison (AIC/BIC for the MCMC; ln B for the nested-sampling runs)
- Convergence and nesting-control diagnostics

### Methodological notes

The parameters **α, β and γ are strongly degenerate** — the three functions
have similar shapes over the observed redshift range, so different
combinations produce nearly identical expansion histories
(Φ₀·(1+α+β) ≈ Ω_Λ). Ensemble samplers (emcee) do not converge on this
posterior; nested sampling (dynesty) does, and also yields the evidence.

A comparison of models must include a **coherence control on the reference
model** (ΛCDM): during this campaign such controls caught two invalid runs
before they produced results — see `CHANGELOG.md`.

### Version history

**v1.8 (current)** — DESI DR2 support (`--dataset`), Bayesian evidence via
nested sampling, and the CMB distance anchor. Result of the full campaign
(BAO+SNe+CMB): ln B = −7.4 (best single-tracer) to −14.8 (3-parameter) in
favour of ΛCDM; γ excluded; couplings bounded below ~1%. Details in
`CHANGELOG.md`.

**v1.7.1** — ~10⁴× speedup (h_ent tabulated + cumulative distance integral),
physicality prior extended to z = 0.011, AIC/BIC at the chain maximum, UTF-8
console output. First full run (DESI DR1): ΔAIC = +10.5 in favour of ΛCDM.

**v1.7** — corrects two critical v1.6 errors: BAO data vector (D_M/r_d
mislabelled as D_V/r_d, producing a spurious +8.4σ pull) and SNIa likelihood
(SH0ES/Planck calibration bias, fixed by analytic marginalisation of the
magnitude offset). Verified: χ²(ΛCDM) = 20.6 / 12 (χ²/n = 1.72).

**v1.6** — initial release, kept as historical record. **No result produced
by v1.6 should be used** (its BAO data vector was invalid).

---

## Español

### Qué hace este código

La THPE promueve la densidad de energía oscura a una función de presión
holográfica dependiente del desplazamiento al rojo

```
Φ(z) = Φ₀ · [ 1 + α·f_SFR(z) + β·g_struct(z) + γ·h_ent(z) ]
```

que recupera ΛCDM exactamente cuando α = β = γ = 0. `THPE_fit_v18.py` la
ajusta a las oscilaciones acústicas de bariones de DESI (DR1 o DR2) y a las
supernovas de tipo Ia de Pantheon+ mediante MCMC; `THPE_dynesty_v1.py` y
`THPE_dynesty_v2_cmb.py` calculan la evidencia bayesiana por muestreo anidado
(el segundo añade el ancla de distancia del CMB de Planck).

### Requisitos

```bash
pip install numpy scipy matplotlib pandas emcee requests corner dynesty
```

Python 3.10+ con NumPy 2.x (el código usa `np.trapezoid`; con NumPy 1.x
sustituir por `np.trapz` o actualizar).

### Datos

**DESI** — los vectores BAO de DR1 (12 observables) y DR2 (13 observables),
con las correlaciones D_M/r_d–D_H/r_d por trazador, van incrustados en
`THPE_fit_v18.py`. Se seleccionan con `--dataset dr1|dr2` (por defecto: dr2).
Valores transcritos de las tablas oficiales; verificar antes del ajuste
definitivo.

**Pantheon+** — los scripts descargan los archivos públicos del repositorio
`PantheonPlusSH0ES/DataRelease` (módulos de distancia + covarianza STAT+SYS
completa). Si la descarga falla, colocarlos en la carpeta del script y usar
`--no-download`. La covarianza completa es imprescindible para un resultado
publicable.

### Cómo ejecutarlo

```bash
python THPE_fit_v18.py --quick               # prueba rápida MCMC
python THPE_fit_v18.py                        # ajuste MCMC completo (DR2)
python THPE_dynesty_v1.py                     # evidencia bayesiana, BAO+SNe
python THPE_dynesty_v2_cmb.py                 # evidencia bayesiana + ancla CMB
```

Al arrancar, cada script ejecuta una verificación de coherencia: el χ²(ΛCDM)
contra el vector de datos debe ser razonable. Si es anómalamente grande, el
script se detiene — eso indica un error de datos, no de física. El script del
CMB además aborta si el χ² de control de ΛCDM es incoherente (protege contra
errores de transcripción o de modelado del fondo).

### Resultados

- Parámetros de mejor ajuste (Φ₀, α, β, γ) con intervalos de credibilidad
- Figuras de resultados y de esquina (`THPE_v18_*.png`)
- Comparación de modelos (AIC/BIC en el MCMC; ln B en las corridas anidadas)
- Diagnósticos de convergencia y de control de anidamiento

### Notas metodológicas

Los parámetros **α, β y γ son fuertemente degenerados** — las tres funciones
tienen formas similares en el rango observado, de modo que distintas
combinaciones producen historias de expansión casi idénticas
(Φ₀·(1+α+β) ≈ Ω_Λ). Los muestreadores de conjunto (emcee) no convergen en
este posterior; el muestreo anidado (dynesty) sí, y además da la evidencia.

Una comparación de modelos debe incluir un **control de coherencia sobre el
modelo de referencia** (ΛCDM): durante esta campaña, esos controles
detectaron dos corridas inválidas antes de que produjeran resultados — ver
`CHANGELOG.md`.

### Histórico de versiones

**v1.8 (actual)** — soporte para DESI DR2 (`--dataset`), evidencia bayesiana
por muestreo anidado y ancla de distancia del CMB. Resultado de la campaña
completa (BAO+SNe+CMB): ln B = −7.4 (mejor trazador) a −14.8 (3 parámetros)
a favor de ΛCDM; γ excluido; acoplamientos acotados por debajo del ~1%.
Detalles en `CHANGELOG.md`.

**v1.7.1** — aceleración ~10⁴× (h_ent tabulada + integral acumulada de
distancias), prior de fisicalidad extendido hasta z = 0.011, AIC/BIC en el
máximo de la cadena, salida UTF-8. Primera corrida completa (DESI DR1):
ΔAIC = +10.5 a favor de ΛCDM.

**v1.7** — corrige dos errores críticos de la v1.6: el vector de datos BAO
(D_M/r_d etiquetado como D_V/r_d, con un pull espurio de +8.4σ) y la
verosimilitud SNIa (sesgo de calibración SH0ES/Planck, resuelto por
marginalización analítica del offset de magnitud). Verificado:
χ²(ΛCDM) = 20.6 / 12 (χ²/n = 1.72).

**v1.6** — versión inicial, conservada como registro histórico. **Ningún
resultado producido por la v1.6 debe utilizarse** (su vector de datos BAO
era inválido).

---

## Citation / Cita

If you use this code, please cite the accompanying paper:
> Mora García, M. & Ordóñez Mora, J. (2026). *The Holographic Theory of State
> Persistence: An Informational Extension of Holographic Dark Energy.*

## License / Licencia

MIT. The authors intend this work to be shared openly for scientific
scrutiny. / MIT. Los autores desean compartir este trabajo abiertamente para
el escrutinio científico.

---

*Code development, auditing and technical writing assisted by AI (Claude,
Anthropic) as a tool, under the authors' responsibility. / Desarrollo,
auditoría y redacción técnica del código asistidos por IA (Claude, Anthropic)
como herramienta, bajo responsabilidad de los autores.*
