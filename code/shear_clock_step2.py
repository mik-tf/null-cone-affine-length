"""Step 2 of the shear-clock line: the sharp 1D problem behind Section 5.1 of arXiv 2611.34924.
Problem (t = ln(v/eps) in [0, Lam], Lam = ln(1/eps), X = Xi, a = hyperbolic speed, X' = a^2/4):
   minimise  u*/rho_S = eps * int_0^Lam e^{t+X(t)} dt   subject to  int_0^Lam a dt = D.
Closed form (parametric in alpha in (0,1), beta > 0), see shear-clock-step2.md:
   Lam = (beta ln(1/(1-alpha)) - alpha ln(1+beta)) / (alpha+beta)
   D   = 2 sqrt(alpha beta)/(alpha+beta) * ln((1+beta)/(1-alpha))
   u_min/rho_S = eps * alpha beta / ((1-alpha)(1+beta))
Run: python3 -I shear_clock_step2.py   (a few seconds, mpmath 30 digits)
"""
import mpmath as mp
mp.mp.dps = 30

def Lam_of(al, be): return (be*mp.log(1/(1-al)) - al*mp.log(1+be))/(al+be)
def D_of(al, be):   return 2*mp.sqrt(al*be)/(al+be)*mp.log((1+be)/(1-al))
def lam0_of(al, be): return al*be/((1-al)*(1+be))

def solve(eps, D):
    """(alpha, beta) with Lam=ln(1/eps), D given; variables x=logit(alpha), y=ln(beta)."""
    Lam = mp.log(1/eps)
    f = lambda x, y: [Lam_of(1/(1+mp.e**-x), mp.e**y) - Lam, D_of(1/(1+mp.e**-x), mp.e**y) - D]
    best = None
    for y0 in [mp.mpf(k) for k in (8, 6, 4, 2, 0, -2)]:
        try:
            x, y = mp.findroot(f, (mp.log(1/eps), y0))
            al, be = 1/(1+mp.e**-x), mp.e**y
            if abs(Lam_of(al, be)-Lam) < 1e-20 and abs(D_of(al, be)-D) < 1e-20:
                best = (al, be); break
        except (ValueError, ZeroDivisionError):
            pass
    return best

def u_geo(eps, D):  # constant-speed geodesic: X(t) = D^2 t/(4 Lam^2)
    Lam = mp.log(1/eps); k = 1 + D**2/(4*Lam**2)
    return eps*(mp.e**(k*Lam) - 1)/k

def check_quadratures(al, be):
    """Recompute Lam and D from the extremal F(lam) = c^2/(4 lam) - lam + K by quadrature."""
    l0 = lam0_of(al, be); p, m = l0/al, l0/be; c = 2*mp.sqrt(p*m); K = p - m
    F = lambda l: c**2/(4*l) - l + K
    Lq = mp.quad(lambda l: 1/F(l), [0, l0]); Dq = mp.quad(lambda l: c/(l*F(l)), [0, l0])
    return Lq, Dq, F(l0)

eps = mp.mpf('1e-3')
print("eps = 1e-3.  D | u_min (exact) | u_geo | SLSQP grid optimum (step 1) | quadrature check")
grid = {1: '1.028316', 3: '1.227536', 6: '2.454592'}
for D in (mp.mpf(1), mp.mpf(3), mp.mpf(6), mp.mpf(10)):
    al, be = solve(eps, D)
    um = eps*lam0_of(al, be)
    Lq, Dq, F0 = check_quadratures(al, be)
    ok = abs(Lq - mp.log(1/eps)) < 1e-15 and abs(Dq - D) < 1e-15 and abs(F0 - 1) < 1e-20
    print(f"  D={int(D):2d}  u_min={mp.nstr(um, 15)}  u_geo={mp.nstr(u_geo(eps, D), 10)}  "
          f"grid={grid.get(int(D), '-')}  alpha={mp.nstr(al, 8)} beta={mp.nstr(be, 8)}  quad ok: {ok}")

# Inverse: largest hyperbolic excursion compatible with a clock bound u*/rho_S <= tau
print("Inverse (eps=1e-3): tau -> D_max(tau), versus the geodesic-only value")
for tau in (mp.mpf('1.01'), mp.mpf(2), mp.mpf(10), mp.mpf(100)):
    Dmax = mp.findroot(lambda D: eps*lam0_of(*solve(eps, D)) - tau, mp.mpf(2)*mp.sqrt(mp.log(1/eps)*mp.log(tau)) + 0.3)
    Dgeo = mp.findroot(lambda D: u_geo(eps, D) - tau, Dmax*0.9)
    print(f"  tau={mp.nstr(tau, 4):5s}  D_max={mp.nstr(Dmax, 12)}  geodesic-speed D={mp.nstr(Dgeo, 8)}  "
          f"sqrt(4 Lam ln tau)={mp.nstr(mp.sqrt(4*mp.log(1/eps)*mp.log(tau)), 8)}")

# The extremal has Xi(1) = infinity with a finite clock: e^{-Xi} = -sigma' -> G(0)=0.
al, be = solve(eps, mp.mpf(3)); l0 = lam0_of(al, be); p, m = l0/al, l0/be; c = 2*mp.sqrt(p*m); K = p-m
for lt in ('1e-2', '1e-4', '1e-8'):
    l = mp.mpf(lt); F = c**2/(4*l) - l + K
    t = mp.quad(lambda s: 1/(c**2/(4*s) - s + K), [l, l0])  # elapsed t
    print(f"  D=3 extremal: remaining clock {lt}: Xi = ln F - t = {mp.nstr(mp.log(F) - (mp.log(1/eps) - (mp.log(1/eps)-t)) , 8)}")
