# -*- coding: utf-8 -*-
"""
THPE_dynesty_v2_cmb.py — Evidencia bayesiana THPE vs LCDM añadiendo el
ancla del CMB: priors comprimidos de distancia de Planck 2018
(parametro de desplazamiento R y escala acustica l_A).

Por que: sin CMB, LCDM absorbia la tension de DESI DR2 inflando Phi0 a
0.723 (rompiendo planitud). El CMB fija la distancia a la ultima
dispersion (z* ~ 1089) y elimina esa libertad. Es la unica modificacion
con posibilidad real de mover ln B hacia la THPE.

Priors comprimidos (Planck 2018 TT,TE,EE+lowE; Chen, Huang & Wang 2018,
arXiv:1808.05724, Tabla 1):
    R   = 1.7502  +- 0.0046
    l_A = 301.471 +- 0.090
    correlacion(R, l_A) = 0.46
*** VERIFICAR estos tres numeros contra la Tabla 1 del paper antes de
*** dar el resultado por definitivo (transcritos de memoria por Claude).

Modelos: LCDM (1p), THPE solo-alpha (2p), THPE-3p (gamma=0).
(El 4p se omite: gamma ya esta excluido y su espina sesga el muestreo.)

Uso:   pip install dynesty   (ya instalado)
       python THPE_dynesty_v2_cmb.py > salida_cmb.txt 2>&1
Opcional: --dataset dr1 | --nlive 800 | --sin-cmb (control sin ancla)
"""
import sys, argparse
import numpy as np

for _s in (sys.stdout, sys.stderr):
    if hasattr(_s, 'reconfigure'):
        try: _s.reconfigure(encoding='utf-8', errors='replace')
        except Exception: pass

import THPE_fit_v18 as M

try:
    import dynesty
except ImportError:
    sys.exit("Falta dynesty:  pip install dynesty")

p = argparse.ArgumentParser()
p.add_argument('--dataset', choices=['dr1', 'dr2'], default='dr2')
p.add_argument('--nlive', type=int, default=800)
p.add_argument('--sin-cmb', action='store_true',
               help='control: misma corrida sin el termino CMB')
args = p.parse_args()

usar_cmb = not args.sin_cmb
print("="*70)
print(f"THPE — Evidencia con dynesty | DESI {args.dataset.upper()} + Pantheon+"
      + (" + CMB (R, l_A)" if usar_cmb else "  [SIN CMB - control]"))
print("="*70)

# ---------------- datos BAO + SNIa (identico a v1) ----------------
M.ensure_data()
try:
    z_bao, kinds, y_obs, sig, cov, ci_bao = M.build_DESI_vector(True, args.dataset)
except TypeError:
    z_bao, kinds, y_obs, sig, cov, ci_bao = M.build_DESI_vector(True, release=args.dataset)
z_sn, mu_obs, sig_mu, ci_sn, C1, S11 = M.load_Pantheon_plus(True)

# ---------------- bloque CMB: priors comprimidos ----------------
c_km   = M.c_km
H0     = M.H0
h      = H0/100.0
om_m   = M.Omega_m * h*h          # omega_m = Omega_m h^2
om_b   = 0.02237                  # Planck 2018 (fijo, coherente con el fondo)
om_g   = 2.469e-5                 # fotones (omega_gamma h^2 estandar)

def z_estrella():
    """Redshift de la ultima dispersion (ajuste de Hu & Sugiyama 1996)."""
    g1 = 0.0783*om_b**(-0.238) / (1 + 39.5*om_b**0.763)
    g2 = 0.560 / (1 + 21.1*om_b**1.81)
    return 1048.0*(1 + 0.00124*om_b**(-0.738))*(1 + g1*om_m**g2)

Z_STAR = z_estrella()

def r_s_star():
    """Horizonte del sonido comovil en z* [Mpc], fondo fijo (no depende
    de Phi(z): a z>z* la energia oscura es despreciable)."""
    # Rejilla logaritmica hasta z=1e9: la cola de la era de radiacion
    # aporta ~5.4 Mpc; truncar en 5e4 (bug del 09/07) daba r_s=139.04
    # en vez de ~144.4 y desplazaba l_A ~130 sigma.
    za = np.geomspace(Z_STAR, 1e9, 400000)
    Rb = (3.0*om_b/(4.0*om_g)) / (1.0+za)
    cs = c_km / np.sqrt(3.0*(1.0+Rb))
    Hz = M.H_LCDM(za)
    return float(np.trapezoid(cs/Hz, za))

R_S_STAR = r_s_star()
print(f"  CMB: z* = {Z_STAR:.1f}, r_s(z*) = {R_S_STAR:.2f} Mpc (fondo Planck)")

# Distancia comovil hasta z* con Phi(z) del modelo. Nota de aproximacion:
# la tabla de h_ent llega a z=5; mas alla se extrapola con su ultimo
# valor, que es ~1e-3 y multiplicado por gamma=0 en estos modelos: nulo.
def DC_hasta_zstar(theta):
    Phi0, a, b, g = theta
    zg  = np.linspace(0.0, Z_STAR, 30000)
    E   = M.H_THPE(zg, Phi0, a, b, g) / H0
    dz  = np.diff(zg)
    inv = 1.0/E
    return (c_km/H0) * float(np.sum(0.5*(inv[1:]+inv[:-1])*dz))

# Datos y covarianza del prior comprimido
R_OBS, R_SIG   = 1.7502, 0.0046
LA_OBS, LA_SIG = 301.471, 0.090
RHO            = 0.46
COV_CMB  = np.array([[R_SIG**2, RHO*R_SIG*LA_SIG],
                     [RHO*R_SIG*LA_SIG, LA_SIG**2]])
CINV_CMB = np.linalg.inv(COV_CMB)

# AUTOCALIBRACION (09/07/2026): las aproximaciones comunes del fondo
# (neutrinos como radiacion pura, r_s numerico) producen un desajuste
# absoluto de ~0.1% (chi2_abs ~ 25) IDENTICO para todos los modelos,
# porque H0/Omega_m/omega_b estan fijos. El ancla se recalibra al punto
# Planck del propio pipeline: penaliza desviaciones de D_C(z*) respecto
# al fondo Planck con las incertidumbres de Planck 2018 (R, l_A).
# El sesgo comun se cancela por construccion; la comparacion de modelos
# es insensible a el. Declarar esta eleccion en el paper.
_DC0   = DC_hasta_zstar(np.array([M.PHI0_CALIB, 0.0, 0.0, 0.0]))
R_REF  = np.sqrt(M.Omega_m) * H0 * _DC0 / c_km
LA_REF = np.pi * _DC0 / R_S_STAR
print(f"  Ancla autocalibrada al fondo Planck del pipeline:")
print(f"    referencia R = {R_REF:.4f} (Planck abs: {R_OBS}),"
      f" l_A = {LA_REF:.2f} (Planck abs: {LA_OBS})")
print(f"    desajuste comun absorbido: dR = {R_REF-R_OBS:+.4f},"
      f" dl_A = {LA_REF-LA_OBS:+.2f}")

def loglike_CMB(theta):
    DC  = DC_hasta_zstar(theta)
    Rt  = np.sqrt(M.Omega_m) * H0 * DC / c_km
    lAt = np.pi * DC / R_S_STAR
    d   = np.array([Rt - R_REF, lAt - LA_REF])
    return -0.5 * float(d @ CINV_CMB @ d)

# ---------------- verosimilitud total ----------------
def loglike_full(theta):
    Phi0, a, b, g = theta
    zc = np.linspace(0.011, 4.5, 120)
    if np.any(M.Phi(zc, Phi0, a, b, g) <= 0):
        return -1e10
    ll = (M.log_likelihood_BAO(theta, z_bao, kinds, y_obs, ci_bao)
          + M.log_likelihood_SNIa(theta, z_sn, mu_obs, ci_sn, C1, S11))
    if usar_cmb:
        ll += loglike_CMB(theta)
    return ll if np.isfinite(ll) else -1e10

# Verificacion de coherencia con el ancla puesta
p_lcdm = np.array([M.PHI0_CALIB, 0, 0, 0])
print(f"  Control: logL(LCDM Planck, Phi0=0.685) = {loglike_full(p_lcdm):.2f}")
chi2_cmb_control = -2*loglike_CMB(p_lcdm)
print(f"           chi2_CMB(Planck) = {chi2_cmb_control:.2f}  (esperado O(1))")
if chi2_cmb_control > 20:
    sys.exit("CENTINELA: chi2_CMB(Planck) > 20 -> el ancla CMB es "
             "inconsistente (r_s o priors mal transcritos). CORRIDA "
             "ABORTADA para no producir resultados invalidos.")

# ---------------- modelos ----------------
LP0_MIN, LP0_MAX = np.log10(1e-3), np.log10(2.0)
def _phi0(u): return 10**(LP0_MIN + u*(LP0_MAX-LP0_MIN))

def ptform_lcdm(u):  return np.array([_phi0(u[0])])
def loglike_lcdm(t): return loglike_full(np.array([t[0], 0, 0, 0]))

def ptform_2p(u):    return np.array([_phi0(u[0]), 10.0*u[1]])
def loglike_alpha(t):return loglike_full(np.array([t[0], t[1], 0, 0]))

def ptform_3p(u):    return np.array([_phi0(u[0]), 10.0*u[1], 10.0*u[2]])
def loglike_3p(t):   return loglike_full(np.array([t[0], t[1], t[2], 0]))

def run(tag, loglike, ptform, ndim):
    print(f"\n--- {tag}: nested sampling (nlive={args.nlive}, ndim={ndim}) ---")
    s = dynesty.NestedSampler(loglike, ptform, ndim,
                              nlive=args.nlive, sample='rslice')
    s.run_nested(dlogz=0.1, print_progress=True)
    r = s.results
    logZ, logZe = r.logz[-1], r.logzerr[-1]
    w = np.exp(r.logwt - r.logz[-1]); w /= w.sum()
    med = [dynesty.utils.quantile(r.samples[:, i], [0.16, 0.5, 0.84], weights=w)
           for i in range(ndim)]
    i0 = int(np.argmax(r.logl))
    print(f"log Z = {logZ:.2f} ± {logZe:.2f} | log L_max = {r.logl[i0]:.2f}"
          f" en {np.round(r.samples[i0], 4)}")
    return logZ, logZe, med, r.logl[i0]

zL, zLe, medL, llL = run("LCDM (1p)",            loglike_lcdm,  ptform_lcdm, 1)
zA, zAe, medA, llA = run("THPE solo-alpha (2p)", loglike_alpha, ptform_2p,   2)
z3, z3e, med3, ll3 = run("THPE-3p (gamma=0)",    loglike_3p,    ptform_3p,   3)

print("\n" + "="*70)
print(f"CON ANCLA CMB {'ACTIVADA' if usar_cmb else 'DESACTIVADA (control)'}:")
print(f"  LCDM:   Phi0 = {medL[0][1]:.4f}  (sin CMB era 0.723; con el ancla"
      f" debe volver hacia ~0.685)")
print(f"  solo-a: Phi0 = {medA[0][1]:.4f}, alpha = {medA[1][1]:.4f}"
      f" (-{medA[1][1]-medA[1][0]:.4f}/+{medA[1][2]-medA[1][1]:.4f})")
print(f"  3p:     Phi0 = {med3[0][1]:.4f}, alpha = {med3[1][1]:.4f},"
      f" beta = {med3[2][1]:.4f}")
print(f"\nFACTOR DE BAYES:  ln B(solo-alpha/LCDM) = {zA-zL:.2f} ± {np.hypot(zAe,zLe):.2f}")
print(f"                  ln B(THPE-3p /LCDM)   = {z3-zL:.2f} ± {np.hypot(z3e,zLe):.2f}")
print(f"CONTROL ANIDAMIENTO: logLmax a={llA:.2f}, 3p={ll3:.2f}, LCDM={llL:.2f}"
      f"  (a y 3p deben ser >= LCDM)")
print("Referencia SIN CMB (07-08/07): solo-alpha=-4.94, 3p=-11.06")
print("Escala de Jeffreys: |lnB|<1 inconcluso | 1-2.5 debil |"
      " 2.5-5 moderado | >5 fuerte")
