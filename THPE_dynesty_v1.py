# -*- coding: utf-8 -*-
"""
THPE_dynesty_v1.py — Evidencia bayesiana y posterior de la THPE con
nested sampling (dynesty), sobre DESI DR2 + Pantheon+ (cov completa).

Por qué: emcee no converge en este posterior (cresta de degeneración
Phi0-alpha-beta con paredes en 0; R-hat 3-8 incluso con 64x20000).
El nested sampling muestrea este tipo de geometrias de forma nativa
y ademas calcula log Z (evidencia) de forma rigurosa, sustituyendo
a la media armonica.

Requiere:  pip install dynesty
Uso:       python THPE_dynesty_v1.py > salida_dynesty.txt 2>&1
           (opcional: --dataset dr1 | --nlive 800)

Calcula log Z para THPE (4 params) y para LCDM (1 param, Phi0 libre)
con los MISMOS datos y verosimilitud, e imprime el factor de Bayes
ln B = ln Z_THPE - ln Z_LCDM con su interpretacion (escala Jeffreys).
"""
import sys, argparse
import numpy as np

for _s in (sys.stdout, sys.stderr):
    if hasattr(_s, 'reconfigure'):
        try: _s.reconfigure(encoding='utf-8', errors='replace')
        except Exception: pass

import THPE_fit_v18 as M   # reutiliza datos, modelo y verosimilitudes v1.8

try:
    import dynesty
except ImportError:
    sys.exit("Falta dynesty:  pip install dynesty")

# ---------------- configuracion ----------------
p = argparse.ArgumentParser()
p.add_argument('--dataset', choices=['dr1', 'dr2'], default='dr2')
p.add_argument('--nlive', type=int, default=800)
args = p.parse_args()

print("="*70)
print(f"THPE — Evidencia bayesiana con dynesty | DESI {args.dataset.upper()} + Pantheon+")
print("="*70)

M.ensure_data()
z_bao, kinds, y_obs, sig, cov, ci_bao = M.build_DESI_vector(True, args.dataset) \
    if 'dataset' in M.build_DESI_vector.__code__.co_varnames else \
    M.build_DESI_vector(True, release=args.dataset)
z_sn, mu_obs, sig_mu, ci_sn, C1, S11 = M.load_Pantheon_plus(True)

def loglike_full(theta):
    """log-verosimilitud comun (BAO + SNIa marginalizada en M_B)."""
    Phi0, a, b, g = theta
    # constricion fisica: Phi>0 en todo el rango con datos
    zc = np.linspace(0.011, 4.5, 120)
    if np.any(M.Phi(zc, Phi0, a, b, g) <= 0):
        return -1e10
    ll = (M.log_likelihood_BAO(theta, z_bao, kinds, y_obs, ci_bao)
          + M.log_likelihood_SNIa(theta, z_sn, mu_obs, ci_sn, C1, S11))
    return ll if np.isfinite(ll) else -1e10

# ---- THPE: 4 parametros ----
LP0_MIN, LP0_MAX = np.log10(1e-3), np.log10(2.0)
def ptform_thpe(u):
    return np.array([10**(LP0_MIN + u[0]*(LP0_MAX-LP0_MIN)),  # Phi0 log-unif
                     10.0*u[1],                                # alpha [0,10]
                     10.0*u[2],                                # beta  [0,10]
                     -0.0 + 5.0*u[3]])                     # gamma [0,5]: gamma<0 es
                     # no-fisico (Phi<0 en z~0.011) y creaba una meseta de -1e10 en
                     # medio prior que sesgaba logZ. Se elimina del prior directamente.
def loglike_thpe(u_theta):
    return loglike_full(u_theta)

# ---- LCDM: 1 parametro (Phi0), alpha=beta=gamma=0 ----
def ptform_lcdm(u):
    return np.array([10**(LP0_MIN + u[0]*(LP0_MAX-LP0_MIN))])
def loglike_lcdm(t):
    return loglike_full(np.array([t[0], 0.0, 0.0, 0.0]))

def run(tag, loglike, ptform, ndim):
    print(f"\n--- {tag}: nested sampling (nlive={args.nlive}, ndim={ndim}) ---")
    s = dynesty.NestedSampler(loglike, ptform, ndim,
                              nlive=args.nlive, sample='rslice')
    s.run_nested(dlogz=0.1, print_progress=True)
    r = s.results
    logZ, logZerr = r.logz[-1], r.logzerr[-1]
    # posterior ponderado
    w = np.exp(r.logwt - r.logz[-1]); w /= w.sum()
    med = np.array([dynesty.utils.quantile(r.samples[:, i], [0.16, 0.5, 0.84], weights=w)
                    for i in range(ndim)])
    imax = int(np.argmax(r.logl))
    print(f"log Z = {logZ:.2f} ± {logZerr:.2f} | log L_max = {r.logl[imax]:.2f}")
    return logZ, logZerr, med, r.samples[imax], r.logl[imax]

zT, zTe, medT, bestT, llT = run("THPE (4p)", loglike_thpe, ptform_thpe, 4)

# ---- THPE-3p: variante gamma=0 (la que el propio ajuste pide) ----
def ptform_thpe3(u):
    return np.array([10**(LP0_MIN + u[0]*(LP0_MAX-LP0_MIN)),
                     10.0*u[1], 10.0*u[2]])
def loglike_thpe3(t):
    return loglike_full(np.array([t[0], t[1], t[2], 0.0]))

z3, z3e, med3, best3, ll3 = run("THPE-3p (gamma=0)", loglike_thpe3, ptform_thpe3, 3)

# ---- THPE-2p: desagregacion por trazador (un solo acoplamiento) ----
def ptform_2p(u):
    return np.array([10**(LP0_MIN + u[0]*(LP0_MAX-LP0_MIN)), 10.0*u[1]])
def loglike_alpha(t):
    return loglike_full(np.array([t[0], t[1], 0.0, 0.0]))
def loglike_beta(t):
    return loglike_full(np.array([t[0], 0.0, t[1], 0.0]))

zA, zAe, medA, bestA, llA = run("THPE-2p (solo alpha: f_SFR)",   loglike_alpha, ptform_2p, 2)
zB, zBe, medB, bestB, llB = run("THPE-2p (solo beta: g_struct)", loglike_beta,  ptform_2p, 2)
zL, zLe, medL, bestL, llL = run("LCDM (1p)", loglike_lcdm, ptform_lcdm, 1)

names = ['Phi0', 'alpha', 'beta', 'gamma']
print("\n" + "="*70)
print("POSTERIOR THPE (percentiles 16 / 50 / 84, ponderados)")
for n, q in zip(names, medT):
    print(f"  {n:<6} = {q[1]:.4f}  (-{q[1]-q[0]:.4f}/+{q[2]-q[1]:.4f})")
print(f"  Punto de maxima verosimilitud: {np.round(bestT,4)}  (logL={llT:.2f})")
print(f"\nLCDM: Phi0 = {medL[0][1]:.4f}  (logL_max={llL:.2f})")

lnB = zT - zL
err = np.hypot(zTe, zLe)
lnB3 = z3 - zL
err3 = np.hypot(z3e, zLe)
print("\n" + "="*70)
print(f"FACTOR DE BAYES:  ln B(THPE-4p/LCDM) = {lnB:.2f} ± {err:.2f}")
print(f"                  ln B(THPE-3p/LCDM) = {lnB3:.2f} ± {err3:.2f}   (gamma=0)")
print(f"                  ln B(solo-alpha/LCDM) = {zA-zL:.2f} ± {np.hypot(zAe,zLe):.2f}   [f_SFR]")
print(f"                  ln B(solo-beta /LCDM) = {zB-zL:.2f} ± {np.hypot(zBe,zLe):.2f}   [g_struct]")
print(f"POSTERIOR solo-alpha: Phi0={medA[0][1]:.4f}, alpha={medA[1][1]:.4f} (-{medA[1][1]-medA[1][0]:.4f}/+{medA[1][2]-medA[1][1]:.4f})")
print(f"POSTERIOR solo-beta : Phi0={medB[0][1]:.4f}, beta ={medB[1][1]:.4f} (-{medB[1][1]-medB[1][0]:.4f}/+{medB[1][2]-medB[1][1]:.4f})")
print(f"CONTROL ANIDAMIENTO: logLmax 4p={llT:.2f}, 3p={ll3:.2f}, a={llA:.2f}, b={llB:.2f}, LCDM={llL:.2f}")
print("  (los tres maximos deben coincidir ~; si THPE < LCDM, muestreo sesgado)")
print("Escala de Jeffreys: |lnB|<1 inconcluso | 1-2.5 debil | 2.5-5 moderado | >5 fuerte")
if lnB > 1:    v = "los datos favorecen THPE"
elif lnB < -1: v = "los datos favorecen LCDM"
else:          v = "inconcluso: los datos no distinguen los modelos"
print(f"Veredicto: {v}.")
print("\nNota de control: logL_max(THPE) debe ser >= logL_max(LCDM)")
print("(modelos anidados). Si dynesty lo cumple, el muestreo es sano —")
print("exactamente lo que emcee no lograba.")
