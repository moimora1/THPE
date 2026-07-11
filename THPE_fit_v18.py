"""
THPE_fit_v18.py
===============
Ajuste estadístico de la Teoría Holográfica de Persistencia de Estados
a los datos de DESI DR1 (BAO) y Pantheon+ (SNIa).

Versión 1.8 (05/07/2026) — novedad sobre v1.7.1: vector BAO de
DESI DR2 (arXiv:2503.14738, Tabla IV, 13 observables verificados
contra la fuente) seleccionable con --release dr2 (defecto) o dr1.
Hereda de v1.7/v1.7.1:

  [CRÍTICO] Vector de datos BAO corregido. La tabla de respaldo de v1.6
  etiquetaba como D_V/r_d valores que en DESI 2024 (Tabla 1,
  arXiv:2404.03002) son D_M/r_d. Con parámetros ΛCDM puros esto
  producía χ² ≈ 186 (7 puntos), con un pull artificial de +8.4σ en
  z = 2.33 que el MCMC habría absorbido inflando α, β, γ.
  v1.7 usa el vector de datos real de DESI DR1: pares (D_M/r_d, D_H/r_d)
  con su coeficiente de correlación por trazador, y D_V/r_d solo donde
  DESI lo publica como único observable (BGS y QSO). 12 puntos en total.

  [CRÍTICO] Marginalización analítica sobre el offset de magnitud
  absoluta M_B en la verosimilitud de SNIa. v1.6 comparaba MU_SH0ES
  (calibrado con la escala SH0ES, H0 ≈ 73) contra un modelo con
  H0 = 67.4 fijo (Planck), introduciendo ~0.17 mag de sesgo directo
  al χ². v1.7 marginaliza el offset constante de forma exacta:
      χ²_marg = Δᵀ C⁻¹ Δ − (Δᵀ C⁻¹ 1)² / (1ᵀ C⁻¹ 1)
  lo que hace el ajuste insensible a la calibración absoluta y a la
  elección de H0 dentro del término de normalización.

  [Mejora] Soporte para la matriz de covarianza completa STAT+SYS de
  Pantheon+ (archivo Pantheon+SH0ES_STAT+SYS.cov) si está en disco;
  en su ausencia, diagonal con advertencia explícita.

  [Mejora] Verificación automática al arranque: χ² de ΛCDM contra el
  vector BAO. Si χ²/n > 5, el programa se detiene: indica un error de
  datos u observables, no de física.

  [Nota] Las URLs de descarga automática de v1.6 no estaban
  verificadas. v1.7 mantiene la descarga de Pantheon+ (repositorio
  GitHub PantheonPlusSH0ES/DataRelease, público) y elimina la de DESI:
  los valores DR1 van incrustados, transcritos de la Tabla 1 del
  artículo. VERIFICAR contra arXiv:2404.03002 antes de publicar.

Autores: Moisés Mora García, Jorge Ordóñez Mora (ingeniero aeronáutico)
Palma de Mallorca, 2026
Asistencia: Claude (Anthropic) se utilizó como herramienta para la
formalización matemática, el desarrollo del código y la redacción.
La responsabilidad sobre el contenido es de los autores humanos.

Uso:
    python THPE_fit_v17.py [--quick] [--no-download] [--no-cov]

Opciones:
    --quick        Ajuste rápido (500 pasos) para prueba
    --no-download  No intentar descarga automática de Pantheon+
    --no-cov       Ignorar covarianzas completas (solo diagonales)

Requisitos:
    pip install numpy scipy emcee matplotlib pandas requests corner
"""

import numpy as np
from scipy import optimize
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

# Consolas Windows (cp1252) no soportan caracteres griegos (χ, Λ, α...).
# Forzamos UTF-8 en stdout/stderr; si el terminal no puede, sustituye
# el carácter en vez de abortar.
for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, 'reconfigure'):
        try:
            _stream.reconfigure(encoding='utf-8', errors='replace')
        except Exception:
            pass

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

# Planck 2018 TT,TE,EE+lowE+lensing
H0         = 67.4      # km/s/Mpc
Omega_m    = 0.315
Omega_r    = 9.15e-5
Omega_L    = 1.0 - Omega_m - Omega_r
r_drag     = 147.09    # Mpc — radio comóvil del sonido en el arrastre

# Condición de calibración THPE (heredada de v1.6):
# Φ₀_calib ≡ Ω_Λ garantiza que THPE con α=β=γ=0 recupera ΛCDM exactamente.
PHI0_CALIB = Omega_L


# ─────────────────────────────────────────────────────────────────────────────
# 2. DATOS BAO DE DESI DR1 — VECTOR CORREGIDO
# ─────────────────────────────────────────────────────────────────────────────
#
# Transcrito de DESI Collaboration 2024, arXiv:2404.03002, Tabla 1.
# Observables por trazador:
#   BGS  y QSO : solo D_V/r_d (señal insuficiente para separar DM y DH)
#   resto      : par (D_M/r_d, D_H/r_d) con coeficiente de correlación r
#
# ⚠ VERIFICAR estos números contra la Tabla 1 del artículo antes del
#   ajuste definitivo. Un error de transcripción aquí invalida todo.

DESI_DR1 = [
    # (nombre,        z_eff,  tipo,  valor,  sigma)
    dict(name='BGS',        z=0.295, kind='DV', val=7.93,  sig=0.15),
    dict(name='LRG1',       z=0.510, kind='DM', val=13.62, sig=0.25),
    dict(name='LRG1',       z=0.510, kind='DH', val=20.98, sig=0.61),
    dict(name='LRG2',       z=0.706, kind='DM', val=16.85, sig=0.32),
    dict(name='LRG2',       z=0.706, kind='DH', val=20.08, sig=0.60),
    dict(name='LRG3+ELG1',  z=0.930, kind='DM', val=21.71, sig=0.28),
    dict(name='LRG3+ELG1',  z=0.930, kind='DH', val=17.88, sig=0.35),
    dict(name='ELG2',       z=1.317, kind='DM', val=27.79, sig=0.69),
    dict(name='ELG2',       z=1.317, kind='DH', val=13.82, sig=0.42),
    dict(name='QSO',        z=1.491, kind='DV', val=26.07, sig=0.67),
    dict(name='Lya',        z=2.330, kind='DM', val=39.71, sig=0.94),
    dict(name='Lya',        z=2.330, kind='DH', val=8.52,  sig=0.17),
]

# Correlación DM–DH dentro de cada trazador (Tabla 1 de DESI 2024)
DESI_CORR = {
    'LRG1':      -0.445,
    'LRG2':      -0.420,
    'LRG3+ELG1': -0.389,
    'ELG2':      -0.444,
    'Lya':       -0.477,
}



# ─── DESI DR2 (v1.8) ───
# Transcrito de DESI Collaboration 2025, arXiv:2503.14738 (v3), Tabla IV.
# Baseline: 7 trazadores, 13 observables. En DR2 el QSO ya tiene par
# DM/DH (en DR1 solo DV). Verificado contra la fuente el 05/07/2026.
DESI_DR2 = [
    dict(name='BGS',        z=0.295, kind='DV', val=7.942,  sig=0.075),
    dict(name='LRG1',       z=0.510, kind='DM', val=13.588, sig=0.167),
    dict(name='LRG1',       z=0.510, kind='DH', val=21.863, sig=0.425),
    dict(name='LRG2',       z=0.706, kind='DM', val=17.351, sig=0.177),
    dict(name='LRG2',       z=0.706, kind='DH', val=19.455, sig=0.330),
    dict(name='LRG3+ELG1',  z=0.934, kind='DM', val=21.576, sig=0.152),
    dict(name='LRG3+ELG1',  z=0.934, kind='DH', val=17.641, sig=0.193),
    dict(name='ELG2',       z=1.321, kind='DM', val=27.601, sig=0.318),
    dict(name='ELG2',       z=1.321, kind='DH', val=14.176, sig=0.221),
    dict(name='QSO',        z=1.484, kind='DM', val=30.512, sig=0.760),
    dict(name='QSO',        z=1.484, kind='DH', val=12.817, sig=0.516),
    dict(name='Lya',        z=2.330, kind='DM', val=38.988, sig=0.531),
    dict(name='Lya',        z=2.330, kind='DH', val=8.632,  sig=0.101),
]

DESI_CORR_DR2 = {
    'LRG1':      -0.459,
    'LRG2':      -0.404,
    'LRG3+ELG1': -0.416,
    'ELG2':      -0.434,
    'QSO':       -0.500,
    'Lya':       -0.431,
}

def build_DESI_vector(use_covariance=True, release='dr2'):
    """
    Construye el vector de datos BAO y su covarianza bloque-diagonal.

    Devuelve:
        z_bao   : array de z_eff (uno por observable)
        kinds   : lista de tipos ('DV' | 'DM' | 'DH')
        y_obs   : valores observados
        cov     : matriz de covarianza (bloques 2×2 por trazador DM/DH,
                  1×1 para DV; correlaciones inter-trazador = 0, como
                  aproxima el propio análisis oficial de DESI DR1)
        cov_inv : inversa
    """
    data  = DESI_DR2 if release == 'dr2' else DESI_DR1
    corrs = DESI_CORR_DR2 if release == 'dr2' else DESI_CORR
    z_bao = np.array([d['z']   for d in data])
    kinds = [d['kind']         for d in data]
    y_obs = np.array([d['val'] for d in data])
    sig   = np.array([d['sig'] for d in data])

    n   = len(y_obs)
    cov = np.diag(sig**2)

    if use_covariance:
        for tracer, r in corrs.items():
            idx = [i for i, d in enumerate(data) if d['name'] == tracer]
            if len(idx) == 2:
                i, j = idx
                cov[i, j] = cov[j, i] = r * sig[i] * sig[j]
        print(f"  DESI {release.upper()}: {n} observables "
              f"({sum(k=='DV' for k in kinds)} DV, "
              f"{sum(k=='DM' for k in kinds)} DM, "
              f"{sum(k=='DH' for k in kinds)} DH), "
              f"covarianza bloque-diagonal con correlaciones DM–DH")
    else:
        print(f"  DESI {release.upper()}: {n} observables, covarianza diagonal (--no-cov)")

    cov_inv = np.linalg.inv(cov)
    return z_bao, kinds, y_obs, sig, cov, cov_inv


# ─────────────────────────────────────────────────────────────────────────────
# 3. DESCARGA DE PANTHEON+
# ─────────────────────────────────────────────────────────────────────────────

PANTHEON_URL = (
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/"
    "DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/"
    "Pantheon%2BSH0ES.dat"
)
PANTHEON_COV_URL = (
    "https://raw.githubusercontent.com/PantheonPlusSH0ES/"
    "DataRelease/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/"
    "Pantheon%2BSH0ES_STAT%2BSYS.cov"
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
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        with open(filename, 'wb') as f:
            f.write(r.content)
        print(f"  {label}: guardado en {filename}")
        return True
    except Exception as e:
        print(f"  {label}: fallo en descarga ({e})")
        return False


def ensure_data(no_download=False):
    if no_download:
        return
    print("\nVerificando datos observacionales...")
    download_data(PANTHEON_URL,     "Pantheon+SH0ES.dat",
                  "Pantheon+ distancias")
    download_data(PANTHEON_COV_URL, "Pantheon+SH0ES_STAT+SYS.cov",
                  "Pantheon+ covarianza")


# ─────────────────────────────────────────────────────────────────────────────
# 4. FUNCIONES AUXILIARES DEL MODELO (sin cambios físicos desde v1.6)
# ─────────────────────────────────────────────────────────────────────────────

def f_SFR(z):
    """Tasa de formación estelar normalizada (Madau & Dickinson 2014)."""
    z    = np.atleast_1d(np.asarray(z, dtype=float))
    psi  = 0.015 * (1+z)**2.7 / (1 + ((1+z)/2.9)**5.6)
    psi0 = 0.015 * 1.0**2.7 / (1 + (1.0/2.9)**5.6)
    return psi / psi0


def g_struct(z):
    """Densidad de estructuras normalizada. Parametrización exp(-z/0.5)."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    return np.exp(-z / 0.5)


def H_LCDM(z):
    z  = np.atleast_1d(np.asarray(z, dtype=float))
    E2 = Omega_m*(1+z)**3 + Omega_r*(1+z)**4 + Omega_L
    return H0 * np.sqrt(np.maximum(E2, 0.0))


def _V_com_taylor(z):
    """Volumen comóvil para z pequeño (serie de Taylor, v1.6)."""
    q0 = Omega_m/2.0 - Omega_L
    j0 = Omega_m
    DH = c_km / H0
    DC = DH * (z - (1+q0)*z**2/2.0 + (2 - q0 + j0)*z**3/6.0)
    return (4*np.pi/3) * DC**3


def comoving_volume(z):
    z     = np.atleast_1d(np.asarray(z, dtype=float))
    out   = np.zeros_like(z)
    z_low = 0.05
    mask_low = z < z_low
    if np.any(mask_low):
        out[mask_low] = _V_com_taylor(z[mask_low])
    for i in np.where(~mask_low)[0]:
        za = np.linspace(0, z[i], 600)
        Ea = H_LCDM(za) / H0
        DC = (c_km/H0) * np.trapezoid(1.0/Ea, za)
        out[i] = (4*np.pi/3) * DC**3
    return out


def h_ent(z):
    """Entropía del horizonte por volumen comóvil, normalizada en z=0.1."""
    z     = np.atleast_1d(np.asarray(z, dtype=float))
    z_ref = 0.10
    S_ent = 1.0 / H_LCDM(z)**2
    V_com = comoving_volume(z)
    V_safe = np.where(V_com > 0, V_com, 1e-30)
    h = S_ent / V_safe
    H_ref = H_LCDM(np.array([z_ref]))[0]
    V_ref = comoving_volume(np.array([z_ref]))[0]
    return h / ((1.0/H_ref**2) / V_ref)


# ─────────────────────────────────────────────────────────────────────────────
# 5. MODELO THPE
# ─────────────────────────────────────────────────────────────────────────────

# ─── Aceleración v1.7.1 ───
# h_ent(z) NO depende de los parámetros del ajuste: se tabula UNA sola
# vez y luego se interpola en espacio log-log (capta la divergencia
# ~z^-3 en z→0). Sin esta tabla, cada evaluación del posterior con las
# 1588 SNe tardaba ~35 s (≈6 días el modo --quick); con ella, ~ms.
# Error relativo en distancias < 2e-3 en el rango permitido por el prior.
_HENT_ZGRID = np.concatenate([np.geomspace(1e-6, 0.1, 500),
                              np.linspace(0.1005, 5.0, 1200)])
_HENT_CACHE = {}

def h_ent_fast(z):
    if 'tabla' not in _HENT_CACHE:
        vals = h_ent(_HENT_ZGRID)
        _HENT_CACHE['tabla'] = (np.log(_HENT_ZGRID), np.log(vals))
    LZ, LH = _HENT_CACHE['tabla']
    z = np.atleast_1d(np.asarray(z, dtype=float))
    return np.exp(np.interp(np.log(np.clip(z, 1e-6, None)), LZ, LH))


def Phi(z, Phi0, alpha, beta, gamma):
    """Φ(z) = Φ₀ · [1 + α·f_SFR + β·g_struct + γ·h_ent]. Φ₀_ΛCDM = Ω_Λ."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    return Phi0 * (1.0 + alpha*f_SFR(z) + beta*g_struct(z)
                   + gamma*h_ent_fast(z))


def H_THPE(z, Phi0, alpha, beta, gamma):
    """H²(z) = H₀²·[Ω_m(1+z)³ + Ω_r(1+z)⁴ + Φ(z)]."""
    z   = np.atleast_1d(np.asarray(z, dtype=float))
    phi = Phi(z, Phi0, alpha, beta, gamma)
    E2  = Omega_m*(1+z)**3 + Omega_r*(1+z)**4 + phi
    return H0 * np.sqrt(np.maximum(E2, 1e-30))


def comoving_distance_THPE(z, Phi0, alpha, beta, gamma):
    """D_C(z) [Mpc].
    v1.7.1: integral acumulada sobre una rejilla fina única en vez de
    una integración independiente por cada z (1588 integrales por
    evaluación con Pantheon+). Mismo resultado, ~10⁴ veces más rápido."""
    z = np.atleast_1d(np.asarray(z, dtype=float))
    z_max = max(float(z.max()), 0.01)
    zg  = np.linspace(0.0, z_max, 6000)
    E   = H_THPE(zg, Phi0, alpha, beta, gamma) / H0
    inv = 1.0 / E
    dz  = np.diff(zg)
    cum = np.concatenate([[0.0],
                          np.cumsum(0.5*(inv[1:] + inv[:-1])*dz)])
    return np.interp(z, zg, (c_km/H0) * cum)


def bao_observables(z, kinds, Phi0, alpha, beta, gamma, rs=r_drag):
    """
    Vector teórico BAO en los mismos observables que el vector de datos.
      DM: D_M/r_d = D_C/r_d          (universo plano)
      DH: D_H/r_d = c/(H(z)·r_d)
      DV: D_V/r_d = [z·D_M²·D_H]^{1/3}/r_d
    """
    z   = np.atleast_1d(np.asarray(z, dtype=float))
    DC  = comoving_distance_THPE(z, Phi0, alpha, beta, gamma)
    H_z = H_THPE(z, Phi0, alpha, beta, gamma)
    DH  = c_km / H_z
    out = np.zeros(len(z))
    for i, k in enumerate(kinds):
        if   k == 'DM':
            out[i] = DC[i] / rs
        elif k == 'DH':
            out[i] = DH[i] / rs
        elif k == 'DV':
            out[i] = (z[i] * DC[i]**2 * DH[i])**(1.0/3.0) / rs
        else:
            raise ValueError(f"Observable BAO desconocido: {k}")
    return out


def distance_modulus_THPE(z, Phi0, alpha, beta, gamma):
    """μ(z) [mag] con la normalización del modelo (H0 de Planck).
    El offset absoluto se marginaliza en la verosimilitud."""
    z  = np.atleast_1d(np.asarray(z, dtype=float))
    DC = comoving_distance_THPE(z, Phi0, alpha, beta, gamma)
    DL = (1+z) * DC
    return 5.0 * np.log10(np.maximum(DL, 1e-10)) + 25.0


def w_eff(z, Phi0, alpha, beta, gamma, dz=0.005):
    """w(z) = −1 + (1+z)/(3Φ) · dΦ/dz."""
    z    = np.atleast_1d(np.asarray(z, dtype=float))
    phi  = Phi(z,    Phi0, alpha, beta, gamma)
    phip = Phi(z+dz, Phi0, alpha, beta, gamma)
    return -1.0 + (1+z) / (3.0*np.maximum(phi, 1e-30)) * (phip-phi)/dz


# ─────────────────────────────────────────────────────────────────────────────
# 6. CARGA DE PANTHEON+
# ─────────────────────────────────────────────────────────────────────────────

def load_Pantheon_plus(use_covariance=True, z_min=0.01):
    """
    Carga Pantheon+ (Scolnic et al. 2022).

    Usa la columna MU_SH0ES pero el offset de calibración absoluta se
    marginaliza analíticamente en la verosimilitud, de modo que el
    resultado no depende de la escala SH0ES ni del H0 del modelo.

    Si existe Pantheon+SH0ES_STAT+SYS.cov (formato: primera línea N,
    luego N² valores), se carga y se recorta a la máscara z > z_min.
    """
    try:
        df     = pd.read_csv('Pantheon+SH0ES.dat', sep=r'\s+', comment='#')
        zcmb   = df['zCMB'].values.astype(float)
        mask   = zcmb > z_min
        z_sn   = zcmb[mask]
        mu_obs = df['MU_SH0ES'].values.astype(float)[mask]
        sig_mu = df['MU_SH0ES_ERR_DIAG'].values.astype(float)[mask]
        print(f"  Pantheon+: {len(z_sn)} supernovas cargadas (z > {z_min})")

        cov = None
        if use_covariance and os.path.exists('Pantheon+SH0ES_STAT+SYS.cov'):
            raw = np.loadtxt('Pantheon+SH0ES_STAT+SYS.cov')
            N   = int(raw[0])
            C   = raw[1:].reshape(N, N)
            if N == len(zcmb):
                idx = np.where(mask)[0]
                cov = C[np.ix_(idx, idx)]
                print(f"  Pantheon+ covarianza STAT+SYS: {cov.shape[0]}² cargada")
            else:
                print(f"  Pantheon+ covarianza: dimensión inesperada "
                      f"({N} vs {len(zcmb)}), usando diagonal")
        if cov is None:
            cov = np.diag(sig_mu**2)
            if use_covariance:
                print("  Pantheon+ covarianza: DIAGONAL — para publicación"
                      " es imprescindible la matriz STAT+SYS completa")
    except (FileNotFoundError, KeyError):
        print("  Pantheon+: usando datos simulados (ΛCDM + ruido σ=0.12 mag)")
        np.random.seed(42)
        z_sn   = np.sort(np.random.uniform(0.015, 2.0, 300))
        mu_ref = distance_modulus_THPE(z_sn, PHI0_CALIB, 0, 0, 0)
        sig_mu = 0.12 * np.ones_like(z_sn)
        mu_obs = mu_ref + np.random.normal(0, sig_mu)
        cov    = np.diag(sig_mu**2)

    cov_inv = np.linalg.inv(cov)
    ones    = np.ones(len(z_sn))
    Cinv_1  = cov_inv @ ones
    S_11    = ones @ Cinv_1          # 1ᵀC⁻¹1, escalar > 0
    return z_sn, mu_obs, sig_mu, cov_inv, Cinv_1, S_11


# ─────────────────────────────────────────────────────────────────────────────
# 7. FUNCIONES DE VEROSIMILITUD
# ─────────────────────────────────────────────────────────────────────────────

def log_likelihood_BAO(params, z_bao, kinds, y_obs, cov_inv):
    """log L = −½ (y−ŷ)ᵀ C⁻¹ (y−ŷ) con el vector mixto DM/DH/DV."""
    Phi0, alpha, beta, gamma = params
    y_th  = bao_observables(z_bao, kinds, Phi0, alpha, beta, gamma)
    delta = y_obs - y_th
    return -0.5 * float(delta @ cov_inv @ delta)


def log_likelihood_SNIa(params, z_sn, mu_obs, cov_inv, Cinv_1, S_11):
    """
    Log-verosimilitud SNIa con marginalización analítica sobre un
    offset constante ΔM (calibración absoluta M_B y/o elección de H0):

        χ²_marg = ΔᵀC⁻¹Δ − (ΔᵀC⁻¹1)² / (1ᵀC⁻¹1)

    Equivale a integrar el offset con prior plano (Goliath et al. 2001;
    práctica estándar en análisis de SNIa).
    """
    Phi0, alpha, beta, gamma = params
    mu_th = distance_modulus_THPE(z_sn, Phi0, alpha, beta, gamma)
    d     = mu_obs - mu_th
    A     = float(d @ cov_inv @ d)
    B     = float(d @ Cinv_1)
    return -0.5 * (A - B*B/S_11)


def best_offset_SNIa(params, z_sn, mu_obs, Cinv_1, S_11):
    """Offset ΔM de máxima verosimilitud (solo para visualización)."""
    mu_th = distance_modulus_THPE(z_sn, *params)
    return float((mu_obs - mu_th) @ Cinv_1) / S_11


def log_prior(params):
    """
    Phi0 : log-uniforme en [1e-3, 2.0]
    alpha: uniforme en [0, 10]
    beta : uniforme en [0, 10]
    gamma: uniforme en [-5, 5]
    Más restricción Φ(z) > 0 en z ∈ [0.1, 4.5].
    """
    Phi0, alpha, beta, gamma = params
    if not (1e-3 < Phi0 < 2.0):    return -np.inf
    if not (0.0 <= alpha <= 10.0): return -np.inf
    if not (0.0 <= beta  <= 10.0): return -np.inf
    if not (-5.0 <= gamma <= 5.0): return -np.inf
    # v1.7.1: la positividad de Φ se exige desde z=0.011 (la SN más
    # cercana tras la máscara), no desde z=0.1. Con la normalización de
    # h_ent (divergente en z→0), esto excluye γ significativamente
    # negativos, que harían E²<0 justo donde hay datos. El prior antiguo
    # los dejaba pasar en silencio y envenenaba las distancias SNIa.
    z_chk = np.linspace(0.011, 4.5, 120)
    if np.any(Phi(z_chk, Phi0, alpha, beta, gamma) <= 0):
        return -np.inf
    return -np.log(Phi0)


def log_posterior(params, z_bao, kinds, y_obs, cov_inv_bao,
                  z_sn, mu_obs, cov_inv_sn, Cinv_1, S_11):
    lp = log_prior(params)
    if not np.isfinite(lp):
        return -np.inf
    ll_bao = log_likelihood_BAO(params, z_bao, kinds, y_obs, cov_inv_bao)
    ll_sn  = log_likelihood_SNIa(params, z_sn, mu_obs, cov_inv_sn,
                                 Cinv_1, S_11)
    if not (np.isfinite(ll_bao) and np.isfinite(ll_sn)):
        return -np.inf
    return lp + ll_bao + ll_sn


# ─────────────────────────────────────────────────────────────────────────────
# 8. VERIFICACIÓN DE COHERENCIA AL ARRANQUE
# ─────────────────────────────────────────────────────────────────────────────

def sanity_check(z_bao, kinds, y_obs, cov_inv):
    """
    Comprueba que ΛCDM (Planck 2018) da un χ² razonable contra el
    vector BAO. Un χ²/n grande indica error de datos u observables
    (exactamente el fallo de la v1.6), no de física.
    """
    p_lcdm = [PHI0_CALIB, 0.0, 0.0, 0.0]
    chi2   = -2.0 * log_likelihood_BAO(p_lcdm, z_bao, kinds, y_obs, cov_inv)
    n      = len(y_obs)
    print(f"\nVerificación: χ²(ΛCDM | BAO) = {chi2:.2f} con n = {n} "
          f"→ χ²/n = {chi2/n:.2f}")
    if chi2 / n > 5.0:
        print("  ✗ χ²/n > 5: el vector de datos no es coherente con el")
        print("    pipeline de observables. Revisar DESI_DR1 antes de ajustar.")
        sys.exit(1)
    print("  ✓ Coherencia datos–observables verificada.")
    return chi2


# ─────────────────────────────────────────────────────────────────────────────
# 9. AJUSTE MCMC
# ─────────────────────────────────────────────────────────────────────────────

def run_MCMC(args_like, n_walkers=64, n_steps=20000, n_burn=6000, seed=42):
    """Sampler MCMC con emcee (Foreman-Mackey et al. 2013)."""
    np.random.seed(seed)
    n_params   = 4
    p0_center  = np.array([PHI0_CALIB, 0.3, 0.3, 0.1])
    # v1.8.1: dispersion inicial ampliada. Con DR2 el posterior es
    # multimodal (R-hat 2-5 en la corrida de 5000 pasos del 05/07/2026)
    # y una nube inicial estrecha deja walkers atrapados en modos
    # distintos. Mas walkers + mas pasos + nube ancha = mezcla real.
    p0_scatter = np.array([0.08, 0.15, 0.15, 0.08])

    p0 = np.empty((n_walkers, n_params))
    for i in range(n_walkers):
        while True:
            cand = p0_center + p0_scatter * np.random.randn(n_params)
            if np.isfinite(log_prior(cand)):
                p0[i] = cand
                break

    print(f"\nMCMC: {n_walkers} walkers × {n_steps} pasos "
          f"(burn-in: {n_burn}, efectivos: {n_walkers*(n_steps-n_burn)})")

    sampler = emcee.EnsembleSampler(n_walkers, n_params, log_posterior,
                                    args=args_like)
    sampler.run_mcmc(p0, n_steps, progress=True)
    print(f"Aceptación media: {np.mean(sampler.acceptance_fraction):.3f} "
          f"(ideal 0.20–0.50)")
    samples = sampler.get_chain(discard=n_burn, flat=True)
    return samples, sampler


# ─────────────────────────────────────────────────────────────────────────────
# 10. DIAGNÓSTICO DE CONVERGENCIA — GELMAN-RUBIN (sin cambios desde v1.6)
# ─────────────────────────────────────────────────────────────────────────────

def gelman_rubin(sampler, n_burn):
    chain = sampler.get_chain(discard=n_burn)
    n_steps_post, n_walkers, n_params = chain.shape
    R_hat = np.zeros(n_params)
    for p in range(n_params):
        cp = chain[:, :, p]
        W  = np.mean(np.var(cp, axis=0, ddof=1))
        B  = n_steps_post * np.var(np.mean(cp, axis=0), ddof=1)
        var_hat  = (1 - 1/n_steps_post) * W + B / n_steps_post
        R_hat[p] = np.sqrt(var_hat / W) if W > 0 else np.inf
    return R_hat


def report_convergence(sampler, n_burn):
    param_names = ['Φ₀', 'α', 'β', 'γ']
    R_hat = gelman_rubin(sampler, n_burn)
    print("\n" + "="*50)
    print("DIAGNÓSTICO DE CONVERGENCIA (Gelman-Rubin)")
    print("="*50)
    print(f"{'Parámetro':<12} {'R̂':>8} {'Estado':>15}")
    print("-"*50)
    for name, rh in zip(param_names, R_hat):
        if   rh < 1.005: estado = "✓ Excelente"
        elif rh < 1.01:  estado = "✓ Buena"
        elif rh < 1.05:  estado = "~ Marginal"
        else:            estado = "✗ No convergido"
        print(f"{name:<12} {rh:>8.4f} {estado:>15}")
    print("="*50)
    if np.all(R_hat < 1.01):
        print("→ Todas las cadenas han convergido.")
    else:
        print("→ Aumentar n_steps para mejorar convergencia.")
    return R_hat


# ─────────────────────────────────────────────────────────────────────────────
# 11. EVIDENCIA BAYESIANA — IMPORTANCIA ARMÓNICA (aproximada)
# ─────────────────────────────────────────────────────────────────────────────

def harmonic_mean_evidence(sampler, n_burn, like_args):
    """
    ADVERTENCIA (heredada de v1.6): el estimador de importancia armónica
    es numéricamente inestable con colas pesadas. Para resultados
    publicables usar MultiNest, PolyChord o dynesty.
    """
    samples = sampler.get_chain(discard=n_burn, flat=True)
    (z_bao, kinds, y_obs, cov_inv_bao,
     z_sn, mu_obs, cov_inv_sn, Cinv_1, S_11) = like_args
    log_L = np.array([
        log_likelihood_BAO(s, z_bao, kinds, y_obs, cov_inv_bao) +
        log_likelihood_SNIa(s, z_sn, mu_obs, cov_inv_sn, Cinv_1, S_11)
        for s in samples[::max(1, len(samples)//4000)]
    ])
    log_L_max = np.max(log_L)
    return log_L_max - np.log(np.mean(np.exp(-(log_L - log_L_max))))


# ─────────────────────────────────────────────────────────────────────────────
# 12. ANÁLISIS ESTADÍSTICO
# ─────────────────────────────────────────────────────────────────────────────

def compute_best_fit(samples):
    med = np.median(samples, axis=0)
    lo  = np.percentile(samples, 16, axis=0)
    hi  = np.percentile(samples, 84, axis=0)
    return med, lo, hi


def compare_models(samples, sampler, n_burn, like_args):
    """AIC, BIC y log-evidencia. ΛCDM de referencia: solo Φ₀ libre."""
    (z_bao, kinds, y_obs, cov_inv_bao,
     z_sn, mu_obs, cov_inv_sn, Cinv_1, S_11) = like_args

    n_data = len(z_bao) + len(z_sn)
    n_THPE, n_LCDM = 4, 1

    # v1.7.1: el punto de referencia para AIC/BIC es el de máxima
    # probabilidad de la cadena, no la mediana. En posteriors
    # degenerados (Φ₀/α/β) la mediana puede caer lejos del pico y
    # sesgar la comparación en contra del modelo con más parámetros.
    logp_chain = sampler.get_log_prob(discard=n_burn, flat=True)
    best    = samples[int(np.argmax(logp_chain))]
    ll_best = (log_likelihood_BAO(best, z_bao, kinds, y_obs, cov_inv_bao) +
               log_likelihood_SNIa(best, z_sn, mu_obs, cov_inv_sn,
                                   Cinv_1, S_11))
    AIC_T = 2*n_THPE - 2*ll_best
    BIC_T = n_THPE*np.log(n_data) - 2*ll_best

    def neg_ll_LCDM(x):
        p = [x[0], 0.0, 0.0, 0.0]
        return -(log_likelihood_BAO(p, z_bao, kinds, y_obs, cov_inv_bao) +
                 log_likelihood_SNIa(p, z_sn, mu_obs, cov_inv_sn,
                                     Cinv_1, S_11))
    res    = optimize.minimize(neg_ll_LCDM, [PHI0_CALIB],
                               bounds=[(0.3, 1.5)], method='L-BFGS-B')
    ll_L   = -res.fun
    AIC_L  = 2*n_LCDM - 2*ll_L
    BIC_L  = n_LCDM*np.log(n_data) - 2*ll_L

    log_Z = harmonic_mean_evidence(sampler, n_burn, like_args)

    print("\n" + "="*70)
    print("COMPARACIÓN DE MODELOS")
    print("="*70)
    print(f"{'Modelo':<10} {'k':>3} {'log L_max':>12} {'AIC':>10} "
          f"{'BIC':>10} {'ΔAIC':>8} {'ΔBIC':>8}")
    print("-"*70)
    print(f"{'ΛCDM':<10} {n_LCDM:>3} {ll_L:>12.2f} "
          f"{AIC_L:>10.2f} {BIC_L:>10.2f} {'—':>8} {'—':>8}")
    dAIC, dBIC = AIC_T - AIC_L, BIC_T - BIC_L
    print(f"{'THPE':<10} {n_THPE:>3} {ll_best:>12.2f} "
          f"{AIC_T:>10.2f} {BIC_T:>10.2f} {dAIC:>8.2f} {dBIC:>8.2f}")
    print("="*70)
    print(f"\nLog-evidencia THPE (importancia armónica): {log_Z:.2f}")
    print("  (aproximado; para publicación usar MultiNest/PolyChord/dynesty)")
    print("\nInterpretación ΔAIC:")
    print("  ΔAIC < -2  → soporte estadístico para THPE sobre ΛCDM")
    print("  ΔAIC ∈ [-2, 2] → modelos estadísticamente equivalentes")
    print("  ΔAIC > +2  → soporte para ΛCDM")
    return dict(AIC_THPE=AIC_T, BIC_THPE=BIC_T, AIC_LCDM=AIC_L,
                BIC_LCDM=BIC_L, dAIC=dAIC, dBIC=dBIC, log_Z_THPE=log_Z)


def print_results(samples):
    med, lo, hi = compute_best_fit(samples)
    names = ['Φ₀', 'α', 'β', 'γ']
    print("\n" + "="*60)
    print("RESULTADOS DEL AJUSTE THPE v1.8")
    print("="*60)
    print(f"{'Parámetro':<10} {'Mediana':>12} {'σ⁻':>10} {'σ⁺':>10}")
    print("-"*60)
    for n, m, l, h in zip(names, med, lo, hi):
        print(f"{n:<10} {m:>12.4f} {m-l:>10.4f} {h-m:>10.4f}")
    print("="*60)
    print(f"\n  Φ₀ = {med[0]:.4f}  (referencia ΛCDM: {PHI0_CALIB:.4f})")
    for k, desc, v in [('α', 'formación estelar', med[1]),
                       ('β', 'estructuras a gran escala', med[2]),
                       ('γ', 'entropía del horizonte', med[3])]:
        print(f"  {k} = {v:>7.4f}  ({desc})")


# ─────────────────────────────────────────────────────────────────────────────
# 13. VISUALIZACIÓN
# ─────────────────────────────────────────────────────────────────────────────

def plot_results(samples, z_bao, kinds, y_obs, sig, z_sn, mu_obs, sig_mu,
                 Cinv_1, S_11, outdir='.', tag='v18'):
    med, _, _ = compute_best_fit(samples)
    Phi0_b, alpha_b, beta_b, gamma_b = med
    p_lcdm = [PHI0_CALIB, 0.0, 0.0, 0.0]

    fig  = plt.figure(figsize=(17, 11))
    gs   = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.35)
    z_pl = np.linspace(0.1, 4.5, 300)

    # Panel 1: BAO D_M/r_d
    ax = fig.add_subplot(gs[0, 0])
    for kind, color, label in [('DM', 'k', r'$D_M/r_d$'),
                               ('DV', 'purple', r'$D_V/r_d$')]:
        m = [i for i, k in enumerate(kinds) if k == kind]
        if m:
            ax.errorbar(z_bao[m], y_obs[m], yerr=sig[m], fmt='o',
                        color=color, ms=5, label=f'DESI {label}', zorder=5)
    kd = ['DM']*len(z_pl)
    ax.plot(z_pl, bao_observables(z_pl, kd, *med), 'dodgerblue', lw=2,
            label='THPE $D_M/r_d$')
    ax.plot(z_pl, bao_observables(z_pl, kd, *p_lcdm), 'tomato', lw=2,
            ls='--', label='ΛCDM')
    ax.set(xlabel=r'$z$', ylabel=r'$D_M/r_d$', title='BAO: DESI DR1 (transversal)')
    ax.legend(fontsize=7); ax.grid(alpha=0.3)

    # Panel 2: BAO D_H/r_d
    ax = fig.add_subplot(gs[0, 1])
    m = [i for i, k in enumerate(kinds) if k == 'DH']
    ax.errorbar(z_bao[m], y_obs[m], yerr=sig[m], fmt='s', color='k',
                ms=5, label=r'DESI $D_H/r_d$', zorder=5)
    kd = ['DH']*len(z_pl)
    ax.plot(z_pl, bao_observables(z_pl, kd, *med), 'dodgerblue', lw=2,
            label='THPE')
    ax.plot(z_pl, bao_observables(z_pl, kd, *p_lcdm), 'tomato', lw=2,
            ls='--', label='ΛCDM')
    ax.set(xlabel=r'$z$', ylabel=r'$D_H/r_d$', title='BAO: DESI DR1 (radial)')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # Panel 3: SNIa con offset marginalizado aplicado
    ax = fig.add_subplot(gs[0, 2])
    dM = best_offset_SNIa(med, z_sn, mu_obs, Cinv_1, S_11)
    idx = np.random.choice(len(z_sn), min(400, len(z_sn)), replace=False)
    ax.errorbar(z_sn[idx], mu_obs[idx] - dM, yerr=sig_mu[idx], fmt='.',
                color='gray', alpha=0.35, ms=3, label='Pantheon+ (−ΔM)')
    ax.plot(z_pl, distance_modulus_THPE(z_pl, *med), 'dodgerblue', lw=2,
            label='THPE')
    ax.plot(z_pl, distance_modulus_THPE(z_pl, *p_lcdm), 'tomato', lw=2,
            ls='--', label='ΛCDM')
    ax.set(xlabel=r'$z$', ylabel=r'$\mu(z)$ [mag]',
           title=f'SNIa: Pantheon+ (ΔM={dM:+.3f} mag marginalizado)')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # Panel 4: w(z)
    ax = fig.add_subplot(gs[1, 0])
    z_w    = np.linspace(0.05, 4.0, 300)
    subsmp = samples[::max(1, len(samples)//300)]
    w_mat  = np.array([w_eff(z_w, *s) for s in subsmp])
    ax.axhline(-1, color='tomato', ls='--', lw=1.5, label=r'ΛCDM $w=-1$')
    ax.fill_between(z_w, np.percentile(w_mat, 16, axis=0),
                    np.percentile(w_mat, 84, axis=0),
                    alpha=0.25, color='dodgerblue')
    ax.plot(z_w, w_eff(z_w, *med), 'dodgerblue', lw=2, label='THPE (mediana)')
    ax.set(xlabel=r'$z$', ylabel=r'$w(z)$', ylim=(-2.2, 0.2),
           title=r'Ecuación de estado $w(z)$')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # Panel 5: Φ(z) y componentes
    ax = fig.add_subplot(gs[1, 1])
    phi_tot = Phi(z_pl, *med)
    ax.plot(z_pl, phi_tot/phi_tot[0], 'k', lw=2.5, label=r'$\Phi(z)/\Phi(0)$')
    ax.plot(z_pl, 1 + alpha_b*f_SFR(z_pl), 'green', lw=1.5, ls='-.',
            label=fr'$1+\alpha f_{{SFR}}$ (α={alpha_b:.2f})')
    ax.plot(z_pl, 1 + beta_b*g_struct(z_pl), 'orange', lw=1.5, ls=':',
            label=fr'$1+\beta g_{{struct}}$ (β={beta_b:.2f})')
    ax.plot(z_pl, 1 + gamma_b*h_ent(z_pl), 'purple', lw=1.5, ls='--',
            label=fr'$1+\gamma h_{{ent}}$ (γ={gamma_b:.2f})')
    ax.set(xlabel=r'$z$', ylabel='contribución relativa',
           title=r'Componentes de $\Phi(z)$')
    ax.legend(fontsize=7); ax.grid(alpha=0.3)

    # Panel 6: posteriores marginales
    ax = fig.add_subplot(gs[1, 2])
    for i, (name, color) in enumerate(zip(['Φ₀', 'α', 'β', 'γ'],
                                          ['k', 'green', 'orange', 'purple'])):
        s = samples[:, i]
        h, edges = np.histogram(s, bins=60, density=True)
        ax.plot(0.5*(edges[1:]+edges[:-1]), h/h.max() + i*1.2,
                color=color, lw=1.5)
        ax.text(edges[0], i*1.2 + 0.5, name, fontsize=10, color=color)
    ax.set(title='Posteriores marginales (normalizadas, apiladas)',
           yticks=[])
    ax.grid(alpha=0.2)

    fname = os.path.join(outdir, f'THPE_{tag}_resultados.png')
    fig.savefig(fname, dpi=150, bbox_inches='tight')
    print(f"\nFigura guardada: {fname}")

    if HAS_CORNER:
        fig2 = corner.corner(samples, labels=['Φ₀', 'α', 'β', 'γ'],
                             quantiles=[0.16, 0.5, 0.84], show_titles=True)
        f2 = os.path.join(outdir, f'THPE_{tag}_corner.png')
        fig2.savefig(f2, dpi=150, bbox_inches='tight')
        print(f"Figura corner guardada: {f2}")


# ─────────────────────────────────────────────────────────────────────────────
# 14. PROGRAMA PRINCIPAL
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Ajuste THPE v1.7')
    parser.add_argument('--quick',       action='store_true')
    parser.add_argument('--no-download', action='store_true')
    parser.add_argument('--no-cov',      action='store_true')
    parser.add_argument('--release', choices=['dr1', 'dr2'], default='dr2',
                        help='Vector BAO de DESI: dr2 (2025, 13 obs, defecto) o dr1 (2024, 12 obs)')
    args = parser.parse_args()

    print("="*70)
    print(f"THPE v1.8 — Ajuste a DESI {args.release.upper()} (BAO) + Pantheon+ (SNIa)")
    print("Mora García & Ordóñez Mora, 2026")
    print("="*70)

    ensure_data(no_download=args.no_download)

    use_cov = not args.no_cov
    z_bao, kinds, y_obs, sig, cov_bao, cov_inv_bao = \
        build_DESI_vector(use_covariance=use_cov, release=args.release)
    z_sn, mu_obs, sig_mu, cov_inv_sn, Cinv_1, S_11 = \
        load_Pantheon_plus(use_covariance=use_cov)

    sanity_check(z_bao, kinds, y_obs, cov_inv_bao)

    like_args = (z_bao, kinds, y_obs, cov_inv_bao,
                 z_sn, mu_obs, cov_inv_sn, Cinv_1, S_11)

    n_steps = 500 if args.quick else 20000
    n_burn  = 100 if args.quick else 6000
    samples, sampler = run_MCMC(like_args, n_steps=n_steps, n_burn=n_burn)

    report_convergence(sampler, n_burn)
    print_results(samples)
    compare_models(samples, sampler, n_burn, like_args)
    plot_results(samples, z_bao, kinds, y_obs, sig,
                 z_sn, mu_obs, sig_mu, Cinv_1, S_11,
                 tag=f'v18_{args.release}')

    print("\nAjuste v1.8 completado.")


if __name__ == '__main__':
    main()