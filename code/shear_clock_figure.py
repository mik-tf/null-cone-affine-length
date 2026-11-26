"""Figure 1 of the paper: D_max (sharp), elementary bound L_el and the Appendix D bound (D.52a) over tau, eps = 1e-3.
Run: python3 -I shear_clock_figure.py [output.pdf]   (about 30 s)"""
import sys
import mpmath as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
mp.mp.dps = 30
eps = mp.mpf('1e-3'); Lam = mp.log(1/eps)

def dmax(tau):
    lam0 = tau/eps
    beta_of = lambda al: (lambda r: r/(1-r))(lam0*(1-al)/al)
    f = lambda al: (lambda b: (b*mp.log(1/(1-al)) - al*mp.log(1+b))/(al+b) - Lam)(beta_of(al))
    lo = lam0/(1+lam0)*(1+mp.mpf('1e-25')); hi = 1-mp.mpf('1e-25')
    for _ in range(300):
        mid = (lo+hi)/2
        if f(mid) > 0: lo = mid
        else: hi = mid
    al = (lo+hi)/2; b = beta_of(al)
    return 2*mp.sqrt(al*b)/(al+b)*mp.log((1+b)/(1-al))
def l_el(tau): return 2*mp.sqrt(tau*Lam/(1-eps))
def l_hold(tau):
    best = min(tau**(mp.mpf(k)/1000)*(1-eps)**(1-mp.mpf(k)/1000)/(mp.mpf(k)/1000) for k in range(1, 1001))
    return 2*mp.sqrt(best/eps)
taus = [mp.mpf(10)**(mp.mpf(k)/8 - mp.mpf('0.0')) for k in range(1, 25)]   # 1.33 .. 1000
taus = [(1-eps) + mp.mpf(10)**(mp.mpf(k)/8 - 3) for k in range(0, 49)]
xs = [float(t) for t in taus]
fig, ax = plt.subplots(figsize=(5.2, 3.3))
ax.loglog([x-float(1-eps) for x in xs], [float(dmax(t)) for t in taus], 'k-', lw=2, label=r'sharp $D_{\max}$ (Theorem 1)')
ax.loglog([x-float(1-eps) for x in xs], [float(l_el(t)) for t in taus], 'b--', label=r'elementary $L_{\rm el}$')
ax.loglog([x-float(1-eps) for x in xs], [float(l_hold(t)) for t in taus], 'r:', label=r'Appendix D, (D.52a)')
ax.set_xlabel(r'clock excess $\tau-(1-\varepsilon)$'); ax.set_ylabel(r'largest hyperbolic length')
ax.set_title(r'$\varepsilon=10^{-3}$', fontsize=9); ax.legend(fontsize=8); ax.grid(alpha=.3)
fig.tight_layout(); fig.savefig(sys.argv[1] if len(sys.argv) > 1 else 'fig1_dmax.pdf')
print("figure written; D_max at tau=2, 100:", mp.nstr(dmax(mp.mpf(2)),5), mp.nstr(dmax(mp.mpf(100)),5))
