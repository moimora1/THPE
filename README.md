# THPE — Holographic Theory of State Persistence
### Cosmological fit code · Código de ajuste cosmológico

**Authors / Autores:** Moisés Mora García, Jorge Ordóñez Mora (aeronautical engineer)
**Location:** Palma de Mallorca, Spain · 2026

This repository contains the statistical fitting code accompanying the paper
*"The Holographic Theory of State Persistence: An Informational Extension of
Holographic Dark Energy"* (Mora García & Ordóñez Mora, 2026).

---

## English

### What this code does

`THPE_fit_v16.py` fits the THPE holographic-pressure function

```
Φ(z) = Φ₀ · [ 1 + α·f_SFR(z) + β·g_struct(z) + γ·h_ent(z) ]
```

to cosmological data — DESI DR1 baryon acoustic oscillations (BAO) and the
Pantheon+ Type Ia supernova sample — using Markov Chain Monte Carlo (MCMC)
sampling. It compares the THPE model against ΛCDM and standard holographic
dark energy via the Akaike (AIC) and Bayesian (BIC) information criteria, and
reports Gelman–Rubin convergence diagnostics.

### Requirements

- Python 3.9 or newer
- `numpy`, `scipy`, `matplotlib`
- `emcee` (MCMC sampler)
- `corner` (posterior corner plots)

Install everything with:

```bash
pip install numpy scipy matplotlib emcee corner
```

### Data

The code expects two public datasets placed in the same folder:

1. **DESI DR1 BAO** — distance measurements and covariance matrix.
2. **Pantheon+** — supernova distance moduli and covariance matrix.

If the data files are not found, the code falls back to representative values
from the published papers so the pipeline can be tested. **For a real fit, the
complete datasets with their covariance matrices are required** (see the note
on degeneracies below).

### How to run

First, a quick test to confirm everything works (≈ 500 steps):

```bash
python THPE_fit_v16.py --quick
```

Then the full fit (this can take from minutes to hours depending on your
machine):

```bash
python THPE_fit_v16.py
```

### Output

- Best-fit parameters (Φ₀, α, β, γ) with credible intervals
- Corner plot of the posterior (`THPE_v16_corner.png`)
- Model comparison table (THPE vs ΛCDM vs standard HDE)
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

---

## Español

### Qué hace este código

`THPE_fit_v16.py` ajusta la función de presión holográfica de la THPE

```
Φ(z) = Φ₀ · [ 1 + α·f_SFR(z) + β·g_struct(z) + γ·h_ent(z) ]
```

a datos cosmológicos — oscilaciones acústicas de bariones (BAO) de DESI DR1 y
la muestra de supernovas de tipo Ia de Pantheon+ — mediante muestreo por
cadenas de Markov Monte Carlo (MCMC). Compara el modelo THPE con ΛCDM y con la
energía oscura holográfica estándar mediante los criterios de información de
Akaike (AIC) y Bayesiano (BIC), y reporta el diagnóstico de convergencia de
Gelman–Rubin.

### Requisitos

- Python 3.9 o posterior
- `numpy`, `scipy`, `matplotlib`
- `emcee` (muestreador MCMC)
- `corner` (gráficos de esquina del posterior)

Instálalo todo con:

```bash
pip install numpy scipy matplotlib emcee corner
```

### Datos

El código espera dos conjuntos de datos públicos en la misma carpeta:

1. **DESI DR1 BAO** — medidas de distancia y matriz de covarianza.
2. **Pantheon+** — módulos de distancia de supernovas y matriz de covarianza.

Si no encuentra los archivos de datos, el código recurre a valores
representativos de los artículos publicados para poder probar el
funcionamiento. **Para un ajuste real se requieren los conjuntos de datos
completos con sus matrices de covarianza** (véase la nota sobre degeneraciones
más abajo).

### Cómo ejecutarlo

Primero, una prueba rápida para confirmar que todo funciona (≈ 500 pasos):

```bash
python THPE_fit_v16.py --quick
```

Luego el ajuste completo (puede tardar de minutos a horas según la máquina):

```bash
python THPE_fit_v16.py
```

### Resultados

- Parámetros de mejor ajuste (Φ₀, α, β, γ) con intervalos de credibilidad
- Gráfico de esquina del posterior (`THPE_v16_corner.png`)
- Tabla de comparación de modelos (THPE vs ΛCDM vs HDE estándar)
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

---

## Citation / Cita

If you use this code, please cite the accompanying paper:

> Mora García, M. & Ordóñez Mora, J. (2026). *The Holographic Theory of State
> Persistence: An Informational Extension of Holographic Dark Energy.*

## License / Licencia

The authors intend this work to be shared openly for scientific scrutiny.
Los autores desean compartir este trabajo abiertamente para el escrutinio
científico.

---

*Bibliographic review and code development assisted by Claude (Anthropic) as a
tool. Responsibility for the content lies with the human authors. / Revisión
bibliográfica y desarrollo del código asistidos por Claude (Anthropic) como
herramienta. La responsabilidad sobre el contenido es de los autores humanos.*
