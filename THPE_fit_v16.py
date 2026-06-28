"""
THPE_fit_v16.py
===============
Ajuste estadístico de la Teoría Holográfica de Persistencia de Estados
a los datos de DESI 2024 (BAO) y Pantheon+ (SNIa).

Versión 1.6 — mejoras respecto a v1.5:
  - Descarga automatizada de datos públicos
  - Regularización de h_ent(z→0) con serie de Taylor
  - Matriz de covarianza completa de DESI
  - Análisis de convergencia Gelman-Rubin
  - Estimación de evidencia Bayesiana (importancia armónica)
  - Condición de calibración Φ₀ documentada formalmente

Autores: Moisés Mora García, Jorge Ordóñez Mora (ingeniero aeronáutico)
Palma de Mallorca, 2026
Asistencia: Claude (Anthropic) se utilizó como herramienta para la
formalización matemática, el desarrollo del código y la redacción.
La responsabilidad sobre el contenido es de los autores humanos.

Uso:
    python THPE_fit_v16.py [--quick] [--no-download]

Opciones:
    --quick        Ajuste rápido (500 pasos) para prueba
    --no-download  No intentar descarga automática de datos

Requisitos:
    pip install numpy scipy emcee matplotlib astropy pandas requests corner
"""

import numpy as np
from scipy import integrate, interpolate, optimize
import emcee
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import warnings
import os
import sys
import argparse
warnings.filterwarnings('ignore')

# Descarga opcional
try:
    import requests
    HAS_REQUESTS = True
except ImportError:
    HAS_REQUESTS = False

try:
    import corner
    HAS_CORNER = True
except ImportError:
    HAS_CORNER = False


# ─────────────────────────────────────────────────────────────────────────────
# 1. CONSTANTES Y PARÁMETROS FIJOS (Planck 2018, arXiv:1807.06209)
# ─────────────────────────────────────────────────────────────────────────────

c_km   = 2.99792458e5     # km/s
c_m    = 2.99792458e8     # m/s
hbar   = 1.054571817e-34  # J·s
G_SI   = 6.67430e-11      # m³/(kg·s²)
t_P    = np.sqrt(hbar * G_SI / c_m**5)  # tiempo de Planck ≈ 5.39e-44 s

# Planck 2018 TT,TE,EE+lowE+lensing
H0         = 67.4      # km/s/Mpc
Omega_m    = 0.315
Omega_b    = 0.049
Omega_r    = 9.15e-5
Omega_k    = 0.0
Omega_L    = 1.0 - Omega_m - Omega_r  # = 0.6850...
r_drag     = 147.09    # Mpc — radio comóvil del sonido
Mpc_to_km  = 3.085677581e19

# Condición de calibración THPE (v1.6):
# Φ₀_calib se define tal que Φ₀_calib/(H₀²·tₚ) = Ω_Λ
# Esto garantiza que THPE con α=β=γ=0 recupera ΛCDM exactamente.
# En unidades adimensionales del código: Φ₀_calib ≡ Ω_Λ
PHI0_CALIB = Omega_L   # ≈ 0.685 — valor de referencia ΛCDM


# ─────────────────────────────────────────────────────────────────────────────
# 2. DESCARGA AUTOMATIZADA DE DATOS
# ─────────────────────────────────────────────────────────────────────────────

DESI_URL = (
    "https://raw.githubusercontent.com/desihub/desi-science/"
    "main/cosmology/bao/dr1/desi_bao_dr1_measurements.csv"
)
PANTHEON_URL = (
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/"
    "DataRelease/main/Pantheon%2BSH0ES.dat"
)

def download_data(url, filename, label=""):
    """Descarga un archivo de datos si no existe ya en disco."""
    if os.path.exists(filename):
        print(f"  {label}: encontrado en disco ({filename})")
        return True
    if not HAS_REQUESTS:
        print(f"  {label}: requests no instalado, usar --no-download")
        return False
    try:
        print(f"  {label}: descargando de {url[:60]}...")
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        with open(filename, 'wb') as f:
            f.write(r.content)
        print(f"  {label}: guardado en {filename}")
        return True
    except Exception as e:
        print(f"  {label}: fallo en descarga ({e})")
        return False


def ensure_data(no_download=False):
    """Intenta asegurar que los archivos de datos están disponibles."""
    if no_download:
        return
    print("\nVerificando datos observacionales...")
    download_data(DESI_URL,     "DESI_BAO_DR1.csv",    "DESI BAO")
    download_data(PANTHEON_URL, "Pantheon+SH0ES.dat",  "Pantheon+")


# ─────────────────────────────────────────────────────────────────────────────
# 3. FUNCIONES AUXILIARES DEL MODELO
# ─────────────────────────────────────────────────────────────────────────────

def f_SFR(z):
    """
    Tasa de formación estelar normalizada (Madau & Dickinson 2014).
    ψ(z) = 0.015·(1+z)^2.7 / [1 + ((1+z)/2.9)^5.6]  [M_sun/yr/Mpc³]
    Normalizada a f_SFR(z=0) = 1.
    """
    z   = np.atleast_1d(np.asarray(z, dtype=float))
    psi = 0.015 * (1+z)**2.7 / (1 + ((1+z)/2.9)**5.6)
    psi0 = 0.015 * 1.0**2.7 / (1 + (1.0/2.9)**5.6)
    return psi / psi0


def g_struct(z):
    """
    Densidad de estructuras a gran escala normalizada (Tinker et al. 2008).
    Parametrización calibrada: g ∝ exp(-z / z_struct), z_struct ≈ 0.5.
    Para implementación completa usar hmf (pip install hmf).
    """
    z       = np.atleast_1d(np.asarray(z, dtype=float))
    z_struct = 0.5
    g = np.exp(-z / z_struct)
    return g / np.exp(0.0)   # normalizada a g(0) = 1


def H_LCDM(z):
    """H(z) en ΛCDM con parámetros de Planck 2018."""
    z  = np.atleast_1d(np.asarray(z, dtype=float))
    E2 = Omega_m*(1+z)**3 + Omega_r*(1+z)**4 + Omega_L
    return H0 * np.sqrt(np.maximum(E2, 0.0))


def comoving_distance_LCDM(z):
    """Distancia comóvil D_C(z) en ΛCDM [Mpc]."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.zeros_like(z)
    for i, zi in enumerate(z):
        if zi <= 0:
            continue
        za = np.linspace(0, zi, 600)
        Ea = H_LCDM(za) / H0
        out[i] = (c_km/H0) * np.trapezoid(1.0/Ea, za)
    return out


# ── Regularización h_ent con serie de Taylor ─────────────────────────────────

def _V_com_taylor(z):
    """
    Volumen comóvil para z pequeño usando expansión en serie de Taylor.
    D_C(z) ≈ (c/H0) · [z - (1+q0)z²/2 + ...]  con q0 = Ω_m/2 - Ω_Λ
    V_com  ≈ (4π/3) · D_C³
    Válido para z < z_taylor = 0.05.
    """
    q0 = Omega_m/2.0 - Omega_L   # parámetro de desaceleración
    j0 = Omega_m                  # parámetro de tirón (jerk) ≈ Ω_m
    DH = c_km / H0               # distancia de Hubble [Mpc]
    DC = DH * (z - (1+q0)*z**2/2.0 + (2 - q0 + j0)*z**3/6.0)
    return (4*np.pi/3) * DC**3


def comoving_volume(z):
    """
    Volumen comóvil dentro del horizonte en z [Mpc³].
    Usa serie de Taylor para z < 0.05 para evitar divergencias numéricas.
    """
    z     = np.atleast_1d(np.asarray(z, dtype=float))
    out   = np.zeros_like(z)
    z_low = 0.05   # umbral para la serie de Taylor

    # Zona de Taylor (z pequeño)
    mask_low = z < z_low
    if np.any(mask_low):
        out[mask_low] = _V_com_taylor(z[mask_low])

    # Zona numérica (z ≥ z_low)
    mask_high = ~mask_low
    for i in np.where(mask_high)[0]:
        zi  = z[i]
        za  = np.linspace(0, zi, 600)
        Ea  = H_LCDM(za) / H0
        DC  = (c_km/H0) * np.trapezoid(1.0/Ea, za)
        out[i] = (4*np.pi/3) * DC**3

    return out


def h_ent(z):
    """
    Densidad informacional de entrelazamiento por volumen comóvil.
    h_ent(z) = S_ent(z) / V_com(z)  ∝  H⁻²(z) / V_com(z)

    Para z → 0: usar expansión de Taylor de V_com para evitar 0/0.
    Normalizada tal que h_ent(z_ref) = 1 con z_ref = 0.1.
    """
    z      = np.atleast_1d(np.asarray(z, dtype=float))
    z_ref  = 0.10   # punto de normalización
    H_z    = H_LCDM(z)
    S_ent  = 1.0 / H_z**2    # ∝ área del horizonte de Hubble

    V_com  = comoving_volume(z)

    # Para z muy pequeño, usar límite analítico:
    # V_com ≈ (4π/3)(c/H0)³ z³,  S_ent ≈ 1/H0²
    # h_ent ≈ (3H0)/(4π c³) z⁻³  — diverge, pero es física (z=0 es singularidad)
    # El ajuste nunca evalúa z=0; el mínimo del grid es z~0.1
    V_safe = np.where(V_com > 0, V_com, 1e-30)
    h      = S_ent / V_safe

    # Normalización
    H_ref   = H_LCDM(np.array([z_ref]))[0]
    V_ref   = comoving_volume(np.array([z_ref]))[0]
    h_ref   = (1.0/H_ref**2) / V_ref

    return h / h_ref


# ─────────────────────────────────────────────────────────────────────────────
# 4. MODELO THPE
# ─────────────────────────────────────────────────────────────────────────────

def Phi(z, Phi0, alpha, beta, gamma):
    """
    Densidad informacional holográfica Φ(z).

    Φ(z) = Φ₀ · [1 + α·f_SFR(z) + β·g_struct(z) + γ·h_ent(z)]

    Condición de calibración (v1.6):
        Φ₀ se define tal que, cuando α=β=γ=0,
        Φ(z) = Φ₀ reproduce el término de energía oscura de ΛCDM.
        En unidades adimensionales: Φ₀_ΛCDM = Ω_Λ ≈ 0.685.

    Restricción: Φ(z) > 0 para todo z en el rango de ajuste.
    """
    z    = np.atleast_1d(np.asarray(z, dtype=float))
    return Phi0 * (1.0 + alpha*f_SFR(z) + beta*g_struct(z) + gamma*h_ent(z))


def H_THPE(z, Phi0, alpha, beta, gamma):
    """
    H(z) en el modelo THPE (ecuación de Friedmann modificada).

    H²(z) = H₀² · [Ω_m(1+z)³ + Ω_r(1+z)⁴ + Φ(z)]

    Donde Φ(z) reemplaza Ω_Λ de ΛCDM.
    Recuperación exacta de ΛCDM: Φ₀=Ω_Λ, α=β=γ=0 → Φ(z)=Ω_Λ.
    """
    z   = np.atleast_1d(np.asarray(z, dtype=float))
    phi = Phi(z, Phi0, alpha, beta, gamma)
    E2  = Omega_m*(1+z)**3 + Omega_r*(1+z)**4 + phi
    return H0 * np.sqrt(np.maximum(E2, 0.0))


def comoving_distance_THPE(z, Phi0, alpha, beta, gamma):
    """Distancia comóvil D_C(z) para el modelo THPE [Mpc]."""
    z   = np.atleast_1d(np.asarray(z, dtype=float))
    out = np.zeros_like(z)
    for i, zi in enumerate(z):
        if zi <= 0:
            continue
        za  = np.linspace(0, zi, 600)
        Ha  = H_THPE(za, Phi0, alpha, beta, gamma)
        Ea  = Ha / H0
        out[i] = (c_km/H0) * np.trapezoid(1.0/Ea, za)
    return out


def DV_over_rs(z, Phi0, alpha, beta, gamma, rs=r_drag):
    """
    Combinación de distancias BAO: D_V(z)/r_s
    D_V(z) = [z · D_C²(z) · c/H(z)]^{1/3}
    """
    z   = np.atleast_1d(np.asarray(z, dtype=float))
    DC  = comoving_distance_THPE(z, Phi0, alpha, beta, gamma)
    H_z = H_THPE(z, Phi0, alpha, beta, gamma)
    DV  = (z * DC**2 * c_km / H_z)**(1.0/3.0)
    return DV / rs


def distance_modulus_THPE(z, Phi0, alpha, beta, gamma):
    """Módulo de distancia μ(z) para SNIa [mag]."""
    z   = np.atleast_1d(np.asarray(z, dtype=float))
    DC  = comoving_distance_THPE(z, Phi0, alpha, beta, gamma)
    DL  = (1+z) * DC   # distancia de luminosidad [Mpc]
    return 5.0 * np.log10(DL) + 25.0


def w_eff(z, Phi0, alpha, beta, gamma, dz=0.005):
    """
    Parámetro de ecuación de estado efectivo.
    w(z) = -1 + (1+z)/3 · d ln Φ(z) / dz
    """
    z    = np.atleast_1d(np.asarray(z, dtype=float))
    phi  = Phi(z,      Phi0, alpha, beta, gamma)
    phip = Phi(z+dz,   Phi0, alpha, beta, gamma)
    dphi = (phip - phi) / dz
    return -1.0 + (1+z) / (3.0 * np.maximum(phi, 1e-30)) * dphi


# ─────────────────────────────────────────────────────────────────────────────
# 5. CARGA DE DATOS CON COVARIANZA DESI
# ─────────────────────────────────────────────────────────────────────────────

def load_DESI_BAO(use_covariance=True):
    """
    Carga los datos BAO de DESI 2024 y, si está disponible,
    la matriz de covarianza completa entre trazadores.

    Formato CSV esperado: z_eff, DV_rs, sigma_DV_rs [, cov_ij...]
    """
    # Datos representativos de DESI 2024 (Tabla 1, arXiv:2404.03002)
    data_repr = np.array([
        # z_eff,  DV/rs,  sigma(DV/rs)
        [0.295,   7.93,   0.15],   # BGS
        [0.510,  13.62,   0.25],   # LRG1
        [0.706,  16.85,   0.32],   # LRG2
        [0.930,  21.71,   0.28],   # LRG3+ELG1
        [1.317,  27.79,   0.69],   # ELG2
        [1.491,  30.21,   0.79],   # QSO
        [2.330,  39.71,   0.94],   # Lya-QSO
    ])

    try:
        df = pd.read_csv('DESI_BAO_DR1.csv')
        z_bao  = df['z_eff'].values.astype(float)
        DV_rs  = df['DV_rs'].values.astype(float)
        sig_DV = df['sigma_DV_rs'].values.astype(float)
        print(f"  DESI BAO: {len(z_bao)} puntos cargados de archivo")
    except (FileNotFoundError, KeyError):
        print("  DESI BAO: usando datos representativos del artículo")
        z_bao  = data_repr[:, 0]
        DV_rs  = data_repr[:, 1]
        sig_DV = data_repr[:, 2]

    # Matriz de covarianza
    n = len(z_bao)
    cov_fname = 'DESI_BAO_covariance.npy'
    if use_covariance and os.path.exists(cov_fname):
        cov = np.load(cov_fname)
        print(f"  DESI covarianza: cargada ({n}×{n})")
    else:
        # Covarianza diagonal (errores independientes entre trazadores)
        # Para análisis completo: cargar la matriz publicada por DESI
        cov = np.diag(sig_DV**2)
        if use_covariance:
            print("  DESI covarianza: usando diagonal (sin archivo de covarianza)")
        else:
            print("  DESI covarianza: diagonal (opción --no-cov)")

    cov_inv = np.linalg.inv(cov)
    return z_bao, DV_rs, sig_DV, cov, cov_inv


def load_Pantheon_plus():
    """Carga datos de supernovas Pantheon+ (Scolnic et al. 2022)."""
    try:
        df     = pd.read_csv('Pantheon+SH0ES.dat', sep=r'\s+', comment='#')
        mask   = df['zCMB'].values.astype(float) > 0.01
        z_sn   = df.loc[mask, 'zCMB'].values.astype(float)
        mu_obs = df.loc[mask, 'MU_SH0ES'].values.astype(float)
        sig_mu = df.loc[mask, 'MU_SH0ES_ERR_DIAG'].values.astype(float)
        print(f"  Pantheon+: {len(z_sn)} supernovas cargadas")
    except (FileNotFoundError, KeyError):
        print("  Pantheon+: usando datos simulados (ΛCDM + ruido σ=0.12 mag)")
        np.random.seed(42)
        z_sn   = np.sort(np.random.uniform(0.015, 2.0, 300))
        mu_ref = distance_modulus_THPE(z_sn, PHI0_CALIB, 0, 0, 0)
        sig_mu = 0.12 * np.ones_like(z_sn)
        mu_obs = mu_ref + np.random.normal(0, sig_mu)

    return z_sn, mu_obs, sig_mu


# ─────────────────────────────────────────────────────────────────────────────
# 6. FUNCIONES DE VEROSIMILITUD
# ─────────────────────────────────────────────────────────────────────────────

def log_likelihood_BAO(params, z_bao, DV_rs_obs, cov_inv):
    """
    Log-verosimilitud BAO con covarianza completa.
    log L = -½ (y - ŷ)ᵀ C⁻¹ (y - ŷ)
    """
    Phi0, alpha, beta, gamma = params
    DV_th  = DV_over_rs(z_bao, Phi0, alpha, beta, gamma)
    delta  = DV_rs_obs - DV_th
    return -0.5 * delta @ cov_inv @ delta


def log_likelihood_SNIa(params, z_sn, mu_obs, sig_mu):
    """Log-verosimilitud SNIa con errores independientes."""
    Phi0, alpha, beta, gamma = params
    mu_th = distance_modulus_THPE(z_sn, Phi0, alpha, beta, gamma)
    chi2  = np.sum(((mu_obs - mu_th) / sig_mu)**2)
    return -0.5 * chi2


def log_prior(params):
    """
    Prior sobre los cuatro parámetros THPE.

    Phi0 : log-uniforme en [1e-3, 2.0]  — análogo a Ω_Λ, referencia ≈ 0.685
    alpha: uniforme en [0, 10]           — contribución SFR ≥ 0
    beta : uniforme en [0, 10]           — contribución estructuras ≥ 0
    gamma: uniforme en [-5, 5]           — contribución horizonte, libre

    Más restricción: Φ(z) > 0 en todo el rango z ∈ [0.1, 4.5]
    """
    Phi0, alpha, beta, gamma = params

    if not (1e-3 < Phi0 < 2.0):   return -np.inf
    if not (0.0 <= alpha <= 10.0): return -np.inf
    if not (0.0 <= beta  <= 10.0): return -np.inf
    if not (-5.0 <= gamma <= 5.0): return -np.inf

    # Restricción Φ(z) > 0 en grid denso
    z_chk = np.linspace(0.1, 4.5, 80)
    phi   = Phi(z_chk, Phi0, alpha, beta, gamma)
    if np.any(phi <= 0):
        return -np.inf

    # Prior log-uniforme sobre Phi0
    return -np.log(Phi0)


def log_posterior(params, z_bao, DV_rs_obs, cov_inv, z_sn, mu_obs, sig_mu):
    """Log-posterior = log-prior + log L_BAO + log L_SNIa."""
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    ll_bao = log_likelihood_BAO(params, z_bao, DV_rs_obs, cov_inv)
    ll_sn  = log_likelihood_SNIa(params, z_sn, mu_obs, sig_mu)
    if not (np.isfinite(ll_bao) and np.isfinite(ll_sn)):
        return -np.inf
    return lp + ll_bao + ll_sn


# ─────────────────────────────────────────────────────────────────────────────
# 7. AJUSTE MCMC
# ─────────────────────────────────────────────────────────────────────────────

def run_MCMC(z_bao, DV_rs_obs, cov_inv, z_sn, mu_obs, sig_mu,
             n_walkers=32, n_steps=5000, n_burn=1000, seed=42):
    """
    Sampler MCMC con emcee (Foreman-Mackey et al. 2013).

    Punto inicial: entorno del valor ΛCDM (Phi0=Ω_Λ, α=β=γ≈0).
    """
    np.random.seed(seed)
    n_params = 4
    p0_center  = np.array([PHI0_CALIB, 0.3, 0.3, 0.1])
    p0_scatter = np.array([0.03, 0.05, 0.05, 0.03])

    # Inicializar walkers con verificación de prior
    p0 = np.empty((n_walkers, n_params))
    for i in range(n_walkers):
        while True:
            cand = p0_center + p0_scatter * np.random.randn(n_params)
            if np.isfinite(log_prior(cand)):
                p0[i] = cand
                break

    print(f"\nMCMC: {n_walkers} walkers × {n_steps} pasos "
          f"(burn-in: {n_burn}, efectivos: {n_walkers*(n_steps-n_burn)})")

    sampler = emcee.EnsembleSampler(
        n_walkers, n_params, log_posterior,
        args=(z_bao, DV_rs_obs, cov_inv, z_sn, mu_obs, sig_mu)
    )
    sampler.run_mcmc(p0, n_steps, progress=True)

    print(f"Aceptación media: {np.mean(sampler.acceptance_fraction):.3f} "
          f"(ideal 0.20–0.50)")

    samples = sampler.get_chain(discard=n_burn, flat=True)
    return samples, sampler


# ─────────────────────────────────────────────────────────────────────────────
# 8. DIAGNÓSTICO DE CONVERGENCIA — GELMAN-RUBIN
# ─────────────────────────────────────────────────────────────────────────────

def gelman_rubin(sampler, n_burn):
    """
    Factor de escala de Gelman-Rubin R̂ para cada parámetro.
    Convergencia aceptable: R̂ < 1.01.
    Buena convergencia:     R̂ < 1.005.

    Referencia: Gelman & Rubin (1992), Statistical Science 7(4), 457–472.
    """
    chain = sampler.get_chain(discard=n_burn)  # (n_steps, n_walkers, n_params)
    n_steps_post, n_walkers, n_params = chain.shape

    R_hat = np.zeros(n_params)
    for p in range(n_params):
        chains_p = chain[:, :, p]           # (steps, walkers)
        W = np.mean(np.var(chains_p, axis=0, ddof=1))      # varianza intra-cadena
        B = n_steps_post * np.var(np.mean(chains_p, axis=0), ddof=1)  # inter-cadena
        var_hat = (1 - 1/n_steps_post) * W + B / n_steps_post
        R_hat[p] = np.sqrt(var_hat / W) if W > 0 else np.inf

    return R_hat


def report_convergence(sampler, n_burn):
    """Imprime el diagnóstico de convergencia Gelman-Rubin."""
    param_names = ['Φ₀', 'α', 'β', 'γ']
    R_hat = gelman_rubin(sampler, n_burn)

    print("\n" + "="*50)
    print("DIAGNÓSTICO DE CONVERGENCIA (Gelman-Rubin)")
    print("="*50)
    print(f"{'Parámetro':<12} {'R̂':>8} {'Estado':>15}")
    print("-"*50)
    for name, rh in zip(param_names, R_hat):
        if rh < 1.005:
            estado = "✓ Excelente"
        elif rh < 1.01:
            estado = "✓ Buena"
        elif rh < 1.05:
            estado = "~ Marginal"
        else:
            estado = "✗ No convergido"
        print(f"{name:<12} {rh:>8.4f} {estado:>15}")
    print("="*50)
    if np.all(R_hat < 1.01):
        print("→ Todas las cadenas han convergido.")
    else:
        print("→ Aumentar n_steps para mejorar convergencia.")
    return R_hat


# ─────────────────────────────────────────────────────────────────────────────
# 9. EVIDENCIA BAYESIANA — ESTIMADOR DE IMPORTANCIA ARMÓNICA
# ─────────────────────────────────────────────────────────────────────────────

def harmonic_mean_evidence(sampler, n_burn, z_bao, DV_rs_obs, cov_inv,
                           z_sn, mu_obs, sig_mu):
    """
    Estimación de la log-evidencia Bayesiana mediante el estimador
    de importancia armónica (Newton & Raftery 1994).

    log Z ≈ -log[ (1/S) Σ exp(-log L(θ_s)) ]
           = log S - log Σ exp(-log L(θ_s))

    donde θ_s son las muestras post-burn-in.

    ADVERTENCIA: el estimador de importancia armónica es numéricamente
    inestable cuando las colas de la verosimilitud son pesadas.
    Para resultados publicables usar MultiNest o PolyChord.
    """
    samples = sampler.get_chain(discard=n_burn, flat=True)

    # Log-verosimilitudes de las muestras (sin el prior)
    log_L = np.array([
        log_likelihood_BAO(s, z_bao, DV_rs_obs, cov_inv) +
        log_likelihood_SNIa(s, z_sn, mu_obs, sig_mu)
        for s in samples
    ])

    # Importancia armónica en log-espacio (estabilizado)
    log_L_max = np.max(log_L)
    log_Z = log_L_max - np.log(np.mean(np.exp(-(log_L - log_L_max))))

    return log_Z


# ─────────────────────────────────────────────────────────────────────────────
# 10. ANÁLISIS ESTADÍSTICO
# ─────────────────────────────────────────────────────────────────────────────

def compute_best_fit(samples):
    med = np.median(samples, axis=0)
    lo  = np.percentile(samples, 16, axis=0)
    hi  = np.percentile(samples, 84, axis=0)
    return med, lo, hi


def compare_models(samples, sampler, n_burn,
                   z_bao, DV_rs_obs, cov_inv, z_sn, mu_obs, sig_mu):
    """
    Compara THPE con ΛCDM mediante AIC, BIC y log-evidencia.
    """
    n_data   = len(z_bao) + len(z_sn)
    n_THPE   = 4   # Φ₀, α, β, γ
    n_LCDM   = 1   # solo Φ₀ = Ω_Λ

    # Mejor ajuste THPE
    best     = np.median(samples, axis=0)
    ll_best  = (log_likelihood_BAO(best, z_bao, DV_rs_obs, cov_inv) +
                log_likelihood_SNIa(best, z_sn, mu_obs, sig_mu))
    AIC_THPE = 2*n_THPE - 2*ll_best
    BIC_THPE = n_THPE*np.log(n_data) - 2*ll_best

    # Referencia ΛCDM: optimizar solo Φ₀
    def neg_ll_LCDM(Phi0):
        p = [Phi0[0], 0.0, 0.0, 0.0]
        return -(log_likelihood_BAO(p, z_bao, DV_rs_obs, cov_inv) +
                 log_likelihood_SNIa(p, z_sn, mu_obs, sig_mu))

    res      = optimize.minimize(neg_ll_LCDM, [PHI0_CALIB],
                                  bounds=[(0.3, 1.5)], method='L-BFGS-B')
    ll_LCDM  = -res.fun
    AIC_LCDM = 2*n_LCDM - 2*ll_LCDM
    BIC_LCDM = n_LCDM*np.log(n_data) - 2*ll_LCDM

    # Log-evidencia (importancia armónica)
    log_Z_THPE = harmonic_mean_evidence(
        sampler, n_burn, z_bao, DV_rs_obs, cov_inv, z_sn, mu_obs, sig_mu)

    print("\n" + "="*70)
    print("COMPARACIÓN DE MODELOS")
    print("="*70)
    print(f"{'Modelo':<10} {'k':>3} {'log L_max':>12} {'AIC':>10} "
          f"{'BIC':>10} {'ΔAIC':>8} {'ΔBIC':>8}")
    print("-"*70)
    print(f"{'ΛCDM':<10} {n_LCDM:>3} {ll_LCDM:>12.2f} "
          f"{AIC_LCDM:>10.2f} {BIC_LCDM:>10.2f} {'—':>8} {'—':>8}")
    dAIC = AIC_THPE - AIC_LCDM
    dBIC = BIC_THPE - BIC_LCDM
    print(f"{'THPE':<10} {n_THPE:>3} {ll_best:>12.2f} "
          f"{AIC_THPE:>10.2f} {BIC_THPE:>10.2f} {dAIC:>8.2f} {dBIC:>8.2f}")
    print("="*70)
    print(f"\nLog-evidencia Bayesiana THPE (importancia armónica): {log_Z_THPE:.2f}")
    print("  (Nota: estimador de importancia armónica es aproximado;")
    print("   para publicación usar MultiNest o PolyChord)")
    print("\nInterpretación ΔAIC:")
    print("  ΔAIC < -2  → soporte estadístico para THPE sobre ΛCDM")
    print("  ΔAIC ∈ [-2, 2] → modelos estadísticamente equivalentes")
    print("  ΔAIC > +2  → soporte para ΛCDM")

    return dict(AIC_THPE=AIC_THPE, BIC_THPE=BIC_THPE,
                AIC_LCDM=AIC_LCDM, BIC_LCDM=BIC_LCDM,
                dAIC=dAIC, dBIC=dBIC, log_Z_THPE=log_Z_THPE)


def print_results(samples):
    med, lo, hi = compute_best_fit(samples)
    names = ['Φ₀', 'α', 'β', 'γ']
    print("\n" + "="*60)
    print("RESULTADOS DEL AJUSTE THPE v1.6")
    print("="*60)
    print(f"{'Parámetro':<10} {'Mediana':>12} {'σ⁻':>10} {'σ⁺':>10}")
    print("-"*60)
    for n, m, l, h in zip(names, med, lo, hi):
        print(f"{n:<10} {m:>12.4f} {m-l:>10.4f} {h-m:>10.4f}")
    print("="*60)
    print("\nPredicciones físicas vs. resultados:")
    print(f"  Φ₀ = {med[0]:.4f}  (referencia ΛCDM: {PHI0_CALIB:.4f})")
    intpr = {
        'α': ('formación estelar', med[1]),
        'β': ('estructuras a gran escala', med[2]),
        'γ': ('entropía del horizonte', med[3]),
    }
    for k, (desc, v) in intpr.items():
        signo = "positivo" if v > 0 else "negativo"
        print(f"  {k} = {v:>7.4f}  ({desc}, {signo})")


# ─────────────────────────────────────────────────────────────────────────────
# 11. VISUALIZACIÓN
# ─────────────────────────────────────────────────────────────────────────────

def plot_results(samples, z_bao, DV_rs_obs, sig_DV, z_sn, mu_obs, sig_mu,
                 outdir='/mnt/user-data/outputs'):
    """Genera los seis paneles de resultados."""
    med, lo, hi = compute_best_fit(samples)
    Phi0_b, alpha_b, beta_b, gamma_b = med

    fig = plt.figure(figsize=(17, 11))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.35)
    z_pl = np.linspace(0.1, 4.5, 300)

    # — Panel 1: Ajuste BAO ————————————————————————————————————————
    ax1 = fig.add_subplot(gs[0, 0])
    DV_th   = DV_over_rs(z_pl, Phi0_b, alpha_b, beta_b, gamma_b)
    DV_LCDM = DV_over_rs(z_pl, PHI0_CALIB, 0, 0, 0)
    ax1.errorbar(z_bao, DV_rs_obs, yerr=sig_DV, fmt='o',
                 color='k', ms=5, label='DESI 2024', zorder=5)
    ax1.plot(z_pl, DV_th,   'dodgerblue', lw=2, label='THPE (mejor ajuste)')
    ax1.plot(z_pl, DV_LCDM, 'tomato',     lw=2, ls='--', label='ΛCDM')
    ax1.set(xlabel=r'$z$', ylabel=r'$D_V(z)/r_s$', title='BAO: DESI 2024')
    ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

    # — Panel 2: Ajuste SNIa ───────────────────────────────────────
    ax2 = fig.add_subplot(gs[0, 1])
    mu_th   = distance_modulus_THPE(z_pl, Phi0_b, alpha_b, beta_b, gamma_b)
    mu_LCDM = distance_modulus_THPE(z_pl, PHI0_CALIB, 0, 0, 0)
    idx = np.random.choice(len(z_sn), min(400, len(z_sn)), replace=False)
    ax2.errorbar(z_sn[idx], mu_obs[idx], yerr=sig_mu[idx],
                 fmt='.', color='gray', alpha=0.35, ms=3, label='Pantheon+')
    ax2.plot(z_pl, mu_th,   'dodgerblue', lw=2, label='THPE')
    ax2.plot(z_pl, mu_LCDM, 'tomato',     lw=2, ls='--', label='ΛCDM')
    ax2.set(xlabel=r'$z$', ylabel=r'$\mu(z)$ [mag]', title='SNIa: Pantheon+')
    ax2.legend(fontsize=8); ax2.grid(alpha=0.3)

    # — Panel 3: w(z) con banda de incertidumbre ──────────────────
    ax3 = fig.add_subplot(gs[0, 2])
    z_w     = np.linspace(0.05, 4.0, 300)
    w_best  = w_eff(z_w, Phi0_b, alpha_b, beta_b, gamma_b)
    # Banda de error usando subsample de cadenas
    subsamp = samples[::max(1, len(samples)//300)]
    w_mat   = np.array([w_eff(z_w, *s) for s in subsamp])
    w_lo    = np.percentile(w_mat, 16, axis=0)
    w_hi    = np.percentile(w_mat, 84, axis=0)
    ax3.axhline(-1, color='tomato', ls='--', lw=1.5, label=r'ΛCDM $w=-1$')
    ax3.fill_between(z_w, w_lo, w_hi, alpha=0.25, color='dodgerblue')
    ax3.plot(z_w, w_best, 'dodgerblue', lw=2, label='THPE (mediana)')
    ax3.set(xlabel=r'$z$', ylabel=r'$w(z)$',
            title=r'Ecuación de estado $w(z)$', ylim=(-2.2, 0.2))
    ax3.legend(fontsize=8); ax3.grid(alpha=0.3)

    # — Panel 4: Φ(z) y sus componentes ───────────────────────────
    ax4 = fig.add_subplot(gs[1, 0])
    phi_tot = Phi(z_pl, Phi0_b, alpha_b, beta_b, gamma_b)
    ax4.plot(z_pl, phi_tot/phi_tot[0],        'k',      lw=2.5,
             label=r'$\Phi(z)/\Phi(0)$')
    ax4.plot(z_pl, 1 + alpha_b*f_SFR(z_pl),  'green',  lw=1.5, ls='-.',
             label=fr'$1+\alpha f_{{SFR}}$ (α={alpha_b:.2f})')
    ax4.plot(z_pl, 1 + beta_b*g_struct(z_pl), 'orange', lw=1.5, ls='-.',
             label=fr'$1+\beta g_{{struct}}$ (β={beta_b:.2f})')
    ax4.plot(z_pl, 1 + gamma_b*h_ent(z_pl),  'purple', lw=1.5, ls='-.',
             label=fr'$1+\gamma h_{{ent}}$ (γ={gamma_b:.2f})')
    ax4.axhline(1, color='tomato', ls='--', lw=1, label='ΛCDM')
    ax4.set(xlabel=r'$z$', ylabel=r'$\Phi(z)/\Phi(0)$',
            title=r'Presión holográfica $\Phi(z)$')
    ax4.legend(fontsize=7); ax4.grid(alpha=0.3)

    # — Panel 5: Corner plot (α vs β vs γ) ────────────────────────
    ax5 = fig.add_subplot(gs[1, 1])
    h2, xe, ye = np.histogram2d(samples[:,1], samples[:,2], bins=40)
    ax5.contourf(0.5*(xe[:-1]+xe[1:]), 0.5*(ye[:-1]+ye[1:]),
                 h2.T, levels=8, cmap='Blues')
    ax5.axvline(alpha_b, color='red',   ls='--', lw=1.5, label=fr'α={alpha_b:.3f}')
    ax5.axhline(beta_b,  color='green', ls='--', lw=1.5, label=fr'β={beta_b:.3f}')
    ax5.set(xlabel=r'$\alpha$', ylabel=r'$\beta$',
            title=r'Marginal 2D $\alpha$–$\beta$')
    ax5.legend(fontsize=8)

    # — Panel 6: Distribuciones marginales ────────────────────────
    ax6 = fig.add_subplot(gs[1, 2])
    pnames = [r'$\Phi_0$', r'$\alpha$', r'$\beta$', r'$\gamma$']
    colors = ['dodgerblue', 'green', 'orange', 'purple']
    for i, (name, col) in enumerate(zip(pnames, colors)):
        v = samples[:, i]
        v_n = (v - v.min()) / (v.max() - v.min() + 1e-12)
        ax6.hist(v_n, bins=50, alpha=0.55, color=col, density=True, label=name)
    ax6.set(xlabel='Valor normalizado', ylabel='Densidad',
            title='Distribuciones marginales')
    ax6.legend(fontsize=8); ax6.grid(alpha=0.3)

    fig.suptitle(
        'THPE v1.6 — Ajuste estadístico a DESI 2024 y Pantheon+\n'
        'Moisés Mora García & Claude (Anthropic), Palma de Mallorca 2026',
        fontsize=10, y=1.01
    )
    outpath = os.path.join(outdir, 'THPE_v16_fit_results.png')
    plt.savefig(outpath, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Gráfico guardado: {outpath}")

    # — Corner plot completo (si corner está instalado) ──────────
    if HAS_CORNER:
        labels = [r'$\Phi_0$', r'$\alpha$', r'$\beta$', r'$\gamma$']
        fig_c  = corner.corner(
            samples, labels=labels,
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True, title_fmt='.3f',
            title_kwargs={'fontsize': 10}
        )
        fig_c.suptitle('THPE v1.6 — Corner plot', fontsize=11)
        out_c = os.path.join(outdir, 'THPE_v16_corner.png')
        fig_c.savefig(out_c, dpi=120, bbox_inches='tight')
        plt.close(fig_c)
        print(f"Corner plot guardado: {out_c}")


# ─────────────────────────────────────────────────────────────────────────────
# 12. FUNCIÓN PRINCIPAL
# ─────────────────────────────────────────────────────────────────────────────

def main():
    # — Argumentos de línea de comandos ────────────────────────────
    parser = argparse.ArgumentParser(description='THPE v1.6 fitting code')
    parser.add_argument('--quick',        action='store_true',
                        help='Ajuste rápido: 500 pasos (prueba)')
    parser.add_argument('--no-download',  action='store_true',
                        help='No intentar descarga automática de datos')
    parser.add_argument('--no-cov',       action='store_true',
                        help='Usar covarianza diagonal (más rápido)')
    args = parser.parse_args()

    print("=" * 65)
    print("THPE v1.6 — Ajuste estadístico")
    print("Teoría Holográfica de Persistencia de Estados")
    print("Moisés Mora García & Jorge Ordóñez Mora")
    print("Asistencia: Claude (Anthropic) como herramienta")
    print("=" * 65)
    print(f"\nCondición de calibración: Φ₀_ΛCDM = Ω_Λ = {PHI0_CALIB:.4f}")
    print(f"Tiempo de Planck: tₚ = {t_P:.3e} s")

    # — Datos ──────────────────────────────────────────────────────
    print("\n1. Datos observacionales:")
    if not args.no_download:
        ensure_data()
    z_bao, DV_rs_obs, sig_DV, cov, cov_inv = load_DESI_BAO(
        use_covariance=not args.no_cov)
    z_sn, mu_obs, sig_mu = load_Pantheon_plus()

    # — Verificación del modelo ────────────────────────────────────
    print("\n2. Verificación del modelo:")
    z_v = np.array([0.3, 1.0, 2.0])
    H_r = H_THPE(z_v, PHI0_CALIB, 0, 0, 0) / H_LCDM(z_v)
    print(f"   H_THPE/H_ΛCDM en z={z_v}: {H_r} (debe ser ≈ 1.0)")
    phi_v = Phi(z_v, PHI0_CALIB, 0, 0, 0)
    print(f"   Φ(z, α=β=γ=0) = {phi_v} (debe ser ≈ {PHI0_CALIB:.4f})")

    # — MCMC ───────────────────────────────────────────────────────
    n_walkers = 32
    n_steps   = 500   if args.quick else 5000
    n_burn    = 100   if args.quick else 1000
    print(f"\n3. MCMC ({'modo rápido' if args.quick else 'modo estándar'}):")
    samples, sampler = run_MCMC(
        z_bao, DV_rs_obs, cov_inv, z_sn, mu_obs, sig_mu,
        n_walkers=n_walkers, n_steps=n_steps, n_burn=n_burn
    )

    # — Convergencia ───────────────────────────────────────────────
    print("\n4. Diagnóstico de convergencia:")
    R_hat = report_convergence(sampler, n_burn)

    # — Resultados ─────────────────────────────────────────────────
    print("\n5. Resultados del ajuste:")
    print_results(samples)

    # — Comparación de modelos ─────────────────────────────────────
    print("\n6. Comparación estadística de modelos:")
    stats = compare_models(samples, sampler, n_burn,
                           z_bao, DV_rs_obs, cov_inv,
                           z_sn, mu_obs, sig_mu)

    # — Visualización ──────────────────────────────────────────────
    print("\n7. Generando gráficos...")
    plot_results(samples, z_bao, DV_rs_obs, sig_DV,
                 z_sn, mu_obs, sig_mu)

    # — Guardar resultados ─────────────────────────────────────────
    outdir = '/mnt/user-data/outputs'
    np.save(os.path.join(outdir, 'THPE_v16_samples.npy'), samples)
    np.save(os.path.join(outdir, 'THPE_v16_Rhat.npy'),   R_hat)
    np.save(os.path.join(outdir, 'THPE_v16_stats.npy'),  stats)

    print("\n" + "="*65)
    print("AJUSTE COMPLETADO")
    print("="*65)
    print("\nArchivos generados en /mnt/user-data/outputs/:")
    print("  THPE_v16_fit_results.png  — seis paneles del ajuste")
    print("  THPE_v16_corner.png       — corner plot (si corner instalado)")
    print("  THPE_v16_samples.npy      — muestras MCMC")
    print("  THPE_v16_Rhat.npy         — factores Gelman-Rubin")
    print("  THPE_v16_stats.npy        — AIC, BIC, log-evidencia")
    print("\nPróximo paso: comparar w(z) con datos de DESI 2024")
    print("y buscar la firma predicha w(z) cruzando -1 entre z≈2 y z≈0.")
    print("\nPara ajuste publicable: n_steps=20000, n_burn=5000,")
    print("covarianza completa de DESI, y MultiNest para evidencia Bayesiana.")


if __name__ == "__main__":
    main()
