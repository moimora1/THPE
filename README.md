# THPE — Holographic Theory of State Persistence
### Cosmological fit code · Código de ajuste cosmológico

**Authors / Autores:** Moisés Mora García, Jorge Ordóñez Mora (aeronautical engineer)
**Location:** Palma de Mallorca, Spain · 2026

This repository contains the statistical fitting code accompanying the paper
*"The Holographic Theory of State Persistence: An Informational Extension of
Holographic Dark Energy"* (Mora García & Ordóñez Mora, 2026).

> **Current version: `THPE_fit_v17.py`** — v1.7 corrects two errors found in
> v1.6 (see *Version history* below). `THPE_fit_v16.py` is kept for the
> historical record; **no result produced by v1.6 should be used**.
>
> **Versión actual: `THPE_fit_v17.py`** — la v1.7 corrige dos errores
> detectados en la v1.6 (véase el *Histórico de versiones*).
> `THPE_fit_v16.py` se conserva como registro histórico; **ningún resultado
> producido por la v1.6 debe utilizarse**.

---

## English

### What this code does

`THPE_fit_v17.py` fits the THPE holographic-pressure function

```
Φ(z) = Φ₀ · [ 1 + α·f_SFR(z) + β·g_struct(z) + γ·h_ent(z) ]
```

to cosmological data — DESI DR1 baryon acoustic oscillations (BAO) and the
Pantheon+ Type Ia supernova sample — using Markov Chain Monte Carlo (MCMC)
sampling. It compares the THPE model against ΛCDM via the Akaike (AIC) and
Bayesian (BIC) information criteria, and reports Gelman–Rubin convergence
diagnostics.

### Requirements

- Python 3.10 or newer, with NumPy 2.x
  (the code uses `np.trapezoid`; on NumPy 1.x replace it with `np.trapz`
  or upgrade: `pip install -U numpy`)
- `numpy`, `scipy`, `matplotlib`, `pandas`
- `emcee` (MCMC sampler)
- `requests` (automatic download of Pantheon+)
- `corner` (posterior corner plots, optional)

Install everything with:

```bash
pip install numpy scipy matplotlib pandas emcee requests corner
```

### Data

**DESI DR1 BAO** — the 12 measurements (D_M/r_d and D_H/r_d pairs with their
per-tracer correlation coefficients, plus D_V/r_d for BGS and QSO) are
embedded in the code (`DESI_DR1` list), transcribed from Table 1 of
arXiv:2404.03002. **Verify them against the paper before a definitive fit.**

**Pantheon+** — the script attempts to download the public files from the
`PantheonPlusSH0ES/DataRelease` GitHub repository:

- `Pantheon+SH0ES.dat` (distance moduli, ~1700 SNe)
- `Pantheon+SH0ES_STAT+SYS.cov` (full covariance matrix)

If the download fails, place both files manually in the same folder as the
script and run with `--no-download`. Without the `.cov` file the fit runs
with diagonal errors, but **the full covariance is required for a publishable
result**.

If no Pantheon+ file is found, the code falls back to simulated data so the
pipeline can be tested — it says so explicitly on screen.

### How to run

Quick test (≈ 500 steps):

```bash
python THPE_fit_v17.py --quick
```

Full fit (32 walkers × 5000 steps; minutes to hours depending on machine):

```bash
python THPE_fit_v17.py
```

On startup the script runs a consistency check: χ²(ΛCDM) against the BAO
vector must be reasonable (expected: χ²/n ≈ 1.7). If χ²/n > 5 the script
stops — that indicates a data or bookkeeping error, not physics.

### Output

- Best-fit parameters (Φ₀, α, β, γ) with credible intervals
- Six-panel results figure (`THPE_v17_resultados.png`)
- Corner plot of the posterior (`THPE_v17_corner.png`)
- Model comparison table (THPE vs ΛCDM: AIC, BIC, ΔAIC, ΔBIC)
- Convergence diagnostics (Gelman–Rubin R̂; aim for R̂ < 1.01)

### Important methodological note

Preliminary tests show that, with limited data, the parameters **α, β and γ
are partially degenerate** — the three functions have similar shapes over the
observed redshift range, so different combinations can produce nearly
identical expansion histories. Two consequences:

- Use the **complete** datasets, not subsamples, and include the **covariance
matrices**.
- Consider **fixing one parameter** (e.g. γ) to sharpen the prediction. Fewer
free parameters means a more predictive, more defensible model.

A good fit confirms that dark energy is dynamical and compatible with
informational complexity, but does **not** by itself demonstrate the formal
independence of THPE from the Holographic Space-Time framework of Banks. That
is a conceptual question, not a statistical one.

### Version history

**v1.7.1 (current)** — ~10⁴× speedup (h_ent tabulated + cumulative
distance integral: a full MCMC now runs in minutes), physicality prior
extended down to z = 0.011, AIC/BIC evaluated at the chain's
maximum-probability point, UTF-8 console output. First full run
(DESI DR1 + Pantheon+ full covariance): ΔAIC = +10.5 in favour of
ΛCDM, γ pinned at zero by physicality — see `CHANGELOG.md`.

**v1.7** corrects two errors in v1.6:

1. **BAO data vector.** The v1.6 fallback table labelled as D_V/r_d values
   that in DESI 2024 (Table 1, arXiv:2404.03002) are D_M/r_d. Against pure
   ΛCDM this produced χ² ≈ 186 over 7 points, with an artificial +8.4σ pull
   at z = 2.33 that the MCMC would have absorbed by inflating α, β, γ.
   v1.7 uses the real DESI DR1 data vector: (D_M/r_d, D_H/r_d) pairs with
   per-tracer correlations, plus D_V/r_d where DESI reports only that
   observable. Verified: χ²(ΛCDM) = 20.6 over 12 points (χ²/n = 1.72), with
   the largest residual pulls at exactly the points where DESI itself
   reported mild tension with ΛCDM.
2. **SNIa likelihood.** v1.6 compared `MU_SH0ES` (calibrated to the SH0ES
   distance scale, H0 ≈ 73) against a model with H0 = 67.4 fixed (Planck),
   injecting ~0.17 mag of bias directly into the χ². v1.7 marginalises the
   constant magnitude offset analytically
   (χ²_marg = ΔᵀC⁻¹Δ − (ΔᵀC⁻¹1)²/(1ᵀC⁻¹1)), making the fit exactly
   invariant under absolute-calibration offsets. Verified numerically.

Full details in the header docstring of `THPE_fit_v17.py` and in
`CHANGELOG.md`.

**v1.6** — kept as historical record. Do not use its results.

---

## Español

### Qué hace este código

`THPE_fit_v17.py` ajusta la función de presión holográfica de la THPE

```
Φ(z) = Φ₀ · [ 1 + α·f_SFR(z) + β·g_struct(z) + γ·h_ent(z) ]
```

a datos cosmológicos — oscilaciones acústicas de bariones (BAO) de DESI DR1 y
la muestra de supernovas de tipo Ia de Pantheon+ — mediante muestreo por
cadenas de Markov Monte Carlo (MCMC). Compara el modelo THPE con ΛCDM
mediante los criterios de información de Akaike (AIC) y Bayesiano (BIC), y
reporta el diagnóstico de convergencia de Gelman–Rubin.

### Requisitos

- Python 3.10 o posterior, con NumPy 2.x
  (el código usa `np.trapezoid`; con NumPy 1.x sustituir por `np.trapz`
  o actualizar: `pip install -U numpy`)
- `numpy`, `scipy`, `matplotlib`, `pandas`
- `emcee` (muestreador MCMC)
- `requests` (descarga automática de Pantheon+)
- `corner` (gráficos de esquina del posterior, opcional)

Instálalo todo con:

```bash
pip install numpy scipy matplotlib pandas emcee requests corner
```

### Datos

**DESI DR1 BAO** — las 12 medidas (pares D_M/r_d y D_H/r_d con sus
coeficientes de correlación por trazador, más D_V/r_d para BGS y QSO) van
incrustadas en el código (lista `DESI_DR1`), transcritas de la Tabla 1 de
arXiv:2404.03002. **Verificarlas contra el artículo antes del ajuste
definitivo.**

**Pantheon+** — el script intenta descargar los archivos públicos del
repositorio GitHub `PantheonPlusSH0ES/DataRelease`:

- `Pantheon+SH0ES.dat` (módulos de distancia, ~1700 SNe)
- `Pantheon+SH0ES_STAT+SYS.cov` (matriz de covarianza completa)

Si la descarga falla, colocar ambos archivos a mano en la misma carpeta que
el script y ejecutar con `--no-download`. Sin el archivo `.cov` el ajuste
funciona con errores diagonales, pero **la covarianza completa es
imprescindible para un resultado publicable**.

Si no encuentra ningún archivo Pantheon+, el código recurre a datos simulados
para poder probar la maquinaria — y lo dice explícitamente por pantalla.

### Cómo ejecutarlo

Prueba rápida (≈ 500 pasos):

```bash
python THPE_fit_v17.py --quick
```

Ajuste completo (32 walkers × 5000 pasos; de minutos a horas según la
máquina):

```bash
python THPE_fit_v17.py
```

Al arrancar, el script ejecuta una verificación de coherencia: el χ²(ΛCDM)
contra el vector BAO debe ser razonable (esperado: χ²/n ≈ 1.7). Si χ²/n > 5
el script se detiene — eso indica un error de datos, no de física.

### Resultados

- Parámetros de mejor ajuste (Φ₀, α, β, γ) con intervalos de credibilidad
- Figura de resultados de seis paneles (`THPE_v17_resultados.png`)
- Gráfico de esquina del posterior (`THPE_v17_corner.png`)
- Tabla de comparación de modelos (THPE vs ΛCDM: AIC, BIC, ΔAIC, ΔBIC)
- Diagnóstico de convergencia (R̂ de Gelman–Rubin; objetivo R̂ < 1.01)

### Nota metodológica importante

Las pruebas preliminares muestran que, con datos limitados, los parámetros
**α, β y γ son parcialmente degenerados** — las tres funciones tienen formas
similares en el rango de desplazamiento al rojo observado, de modo que
distintas combinaciones pueden producir historias de expansión casi
idénticas. Dos consecuencias:

- Usa los conjuntos de datos **completos**, no submuestras, e incluye las
**matrices de covarianza**.
- Considera **fijar un parámetro** (por ejemplo γ) para afinar la predicción.
Menos parámetros libres significa un modelo más predictivo y más defendible.

Un buen ajuste confirma que la energía oscura es dinámica y compatible con la
complejidad informacional, pero **no** demuestra por sí solo la independencia
formal de la THPE respecto al marco Holographic Space-Time de Banks. Esa es
una cuestión conceptual, no estadística.

### Histórico de versiones

**v1.7.1 (actual)** — aceleración ~10⁴× (h_ent tabulada + integral
acumulada de distancias: el MCMC completo corre en minutos), prior de
fisicalidad extendido hasta z = 0.011, AIC/BIC evaluado en el punto de
máxima probabilidad de la cadena, salida UTF-8. Primera corrida
completa (DESI DR1 + covarianza completa de Pantheon+): ΔAIC = +10.5
a favor de ΛCDM, γ clavado en cero por fisicalidad — ver `CHANGELOG.md`.

**v1.7** corrige dos errores de la v1.6:

1. **Vector de datos BAO.** La tabla de respaldo de v1.6 etiquetaba como
   D_V/r_d valores que en DESI 2024 (Tabla 1, arXiv:2404.03002) son D_M/r_d.
   Contra ΛCDM puro esto producía χ² ≈ 186 con 7 puntos, con un pull
   artificial de +8.4σ en z = 2.33 que el MCMC habría absorbido inflando
   α, β, γ. La v1.7 usa el vector real de DESI DR1: pares (D_M/r_d, D_H/r_d)
   con correlaciones por trazador, más D_V/r_d donde DESI solo publica ese
   observable. Verificado: χ²(ΛCDM) = 20.6 con 12 puntos (χ²/n = 1.72), con
   los residuos mayores exactamente en los puntos donde el propio DESI
   reportó tensión suave con ΛCDM.
2. **Verosimilitud SNIa.** La v1.6 comparaba `MU_SH0ES` (calibrado a la
   escala de distancias SH0ES, H0 ≈ 73) contra un modelo con H0 = 67.4 fijo
   (Planck), inyectando ~0.17 mag de sesgo directo al χ². La v1.7 marginaliza
   analíticamente el offset constante de magnitud
   (χ²_marg = ΔᵀC⁻¹Δ − (ΔᵀC⁻¹1)²/(1ᵀC⁻¹1)), haciendo el ajuste exactamente
   invariante ante offsets de calibración absoluta. Verificado numéricamente.

Detalles completos en el docstring de cabecera de `THPE_fit_v17.py` y en
`CHANGELOG.md`.

**v1.6** — se conserva como registro histórico. No usar sus resultados.

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

*Bibliographic review and code development assisted by Claude (Anthropic) as a
tool. Responsibility for the content lies with the human authors. / Revisión
bibliográfica y desarrollo del código asistidos por Claude (Anthropic) como
herramienta. La responsabilidad sobre el contenido es de los autores humanos.*
