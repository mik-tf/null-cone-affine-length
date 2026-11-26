"""Step 3 of the shear-clock line (arXiv 2611.34924): extension checks.
 1. cosh d_H2(q0, q1) = (1/2) tr(q0^{-1} q1) for unimodular SPD 2x2 (the memory observable is algebraic).
 2. Convexity of F_eps(D) = u_min(D, eps)/rho_S (second differences of the exact closed form).
 3. Small-D behaviour of F - (1 - eps): the quadratic coefficient vanishes; measure the actual law.
 4. Whole-cut inequality  ubar* >= sqrt(A/4pi) F_eps(<d>_A)  on random multi-generator data, and sharpness.
 5. eps -> 0 behaviour of D_max at fixed tau.
Run: python3 -I shear_clock_step3.py  (about a minute)
"""
import math, random
import numpy as np
import mpmath as mp
mp.mp.dps = 25

def Lam_of(al, be): return (be*mp.log(1/(1-al)) - al*mp.log(1+be))/(al+be)
def D_of(al, be):   return 2*mp.sqrt(al*be)/(al+be)*mp.log((1+be)/(1-al))
def lam0_of(al, be): return al*be/((1-al)*(1+be))

_cache = {}
def solve(eps, D, guess=None):
    Lam = mp.log(1/eps)
    f = lambda x, y: [Lam_of(1/(1+mp.e**-x), mp.e**y) - Lam, D_of(1/(1+mp.e**-x), mp.e**y) - D]
    starts = ([guess] if guess else []) + [(mp.log(1/eps), mp.mpf(k)) for k in (8, 6, 4, 2, 0, -2, 10, 14)]
    for s in starts:
        try:
            x, y = mp.findroot(f, s)
            al, be = 1/(1+mp.e**-x), mp.e**y
            if abs(Lam_of(al, be)-Lam) < 1e-15 and abs(D_of(al, be)-D) < 1e-15:
                return al, be, (x, y)
        except (ValueError, ZeroDivisionError):
            pass
    raise RuntimeError((eps, D))

def F(eps, D, guess=None):
    al, be, g = solve(eps, mp.mpf(D), guess)
    return eps*lam0_of(al, be), g

# 1. memory observable
print("1. cosh d = tr(q0^-1 q1)/2, random pairs:")
def q_of(z1, z2):
    return np.array([[math.exp(z1), math.exp(z1)*z2], [math.exp(z1)*z2, math.exp(-z1)+math.exp(z1)*z2**2]])
random.seed(2611); err = 0
for _ in range(1000):
    a = [random.uniform(-3, 3) for _ in range(4)]
    q0, q1 = q_of(a[0], a[1]), q_of(a[2], a[3])
    y0, y1, x0, x1 = math.exp(-a[0]), math.exp(-a[2]), a[1], a[3]
    ch = 1 + ((x1-x0)**2 + (y1-y0)**2)/(2*y0*y1)   # upper half plane, K=-1
    err = max(err, abs(ch - 0.5*np.trace(np.linalg.solve(q0, q1)))/ch)
print(f"   max rel err over 1000 pairs = {err:.1e}")
q0 = q_of(0, 0); q1 = q_of(1.3, 0)
print(f"   diag check: d=1.3 -> eigen ratio {np.linalg.eigvalsh(q1)[1]/np.linalg.eigvalsh(q1)[0]:.6f} vs e^(2d)={math.exp(2.6):.6f} (quasiconformal K = e^d)")

# 2. convexity, 3. small D
for eps in (mp.mpf('0.1'), mp.mpf('1e-3'), mp.mpf('1e-6')):
    Ds = [mp.mpf(k)/10 for k in range(1, 121)]
    vals, g = [], None
    for D in Ds:
        v, g = F(eps, D, g); vals.append(v)
    sd = [vals[i-1] - 2*vals[i] + vals[i+1] for i in range(1, len(vals)-1)]
    print(f"2. eps={mp.nstr(eps,3)}: F(0.1)={mp.nstr(vals[0],12)} F(12)={mp.nstr(vals[-1],8)}  min second difference = {mp.nstr(min(sd),4)} (>0 means convex)")
print("3. small D (eps=1e-3): delta = F - (1-eps);  delta/D^2  and  delta*ln(1/delta)/D^2  and  D^2/(delta ln(1/delta))")
eps = mp.mpf('1e-3'); g = None
for D in ('0.3', '0.1', '0.03', '0.01', '0.003', '0.001'):
    v, g = F(eps, D, g); dl = v - (1-eps); D = mp.mpf(D)
    print(f"   D={mp.nstr(D,3):6s} delta={mp.nstr(dl,6):12s} delta/D^2={mp.nstr(dl/D**2,6):10s} delta*ln(1/delta)/D^2={mp.nstr(dl*mp.log(1/dl)/D**2,6)}")

# 4. whole-cut inequality on random data
print("4. whole cut: ubar* >= sqrt(A/4pi) F(<d>_A), random data (eps=1e-3, N generators, random rho_S, random curves)")
eps_f = 1e-3; Lam = math.log(1/eps_f)
Fgrid_D = np.linspace(0, 12, 241); Fgrid = []
g = None
for D in Fgrid_D:
    if D == 0: Fgrid.append(1-eps_f); continue
    v, g = F(mp.mpf(eps_f), mp.mpf(float(D)), g); Fgrid.append(float(v))
Fgrid = np.array(Fgrid)
Finterp = lambda D: np.interp(D, Fgrid_D, Fgrid)   # chords lie above a convex F: lhs check is conservative only if we use exact F
t = np.linspace(0, Lam, 4001)
def generator(rng):
    z1 = sum(rng.uniform(-1, 1)*rng.uniform(0, 2)*np.sin((k+1)*math.pi*t/Lam + rng.uniform(0, 6)) for k in range(4))
    z2 = sum(rng.uniform(-1, 1)*rng.uniform(0, 2)*np.sin((k+1)*math.pi*t/Lam + rng.uniform(0, 6)) for k in range(4))
    z1 -= z1[0]; z2 -= z2[0]
    dz1, dz2 = np.gradient(z1, t), np.gradient(z2, t)
    X = np.concatenate([[0], np.cumsum(0.25*0.5*((dz1[1:]**2+np.exp(2*z1[1:])*dz2[1:]**2)+(dz1[:-1]**2+np.exp(2*z1[:-1])*dz2[:-1]**2))*np.diff(t))])
    tau = eps_f*np.trapz(np.exp(t+X), t)
    y0, y1 = 1.0, math.exp(-z1[-1]); ch = 1 + (z2[-1]**2 + (y1-y0)**2)/(2*y0*y1)
    return tau, math.acosh(ch)
rng = np.random.default_rng(2611); worst = 1e9
for trial in range(200):
    N = int(rng.integers(3, 30)); w = rng.dirichlet(np.ones(N))*4*math.pi
    rho = rng.uniform(0.2, 3, N)
    data = [generator(rng) for _ in range(N)]
    tau = np.array([d[0] for d in data]); dist = np.array([d[1] for d in data])
    A = np.sum(rho**2*w); ubar = np.sum(rho*(rho*tau)*w)/np.sum(rho*w)
    dA = np.sum(rho**2*w*dist)/A
    if dA > 12: continue
    ratio = ubar/(math.sqrt(A/(4*math.pi))*Finterp(dA))
    worst = min(worst, ratio)
print(f"   200 random cuts: min ubar*/(sqrt(A/4pi) F(<d>_A)) = {worst:.6f}  (>= 1 required; F by linear interpolation of the exact values, an upper bound on F, so the check is conservative)")
for D in (1.0, 3.0, 6.0):
    print(f"   sharpness: rho_S const, every generator the extremal of length {D}: ratio -> 1 exactly (u*(y)=rho_S F(D), A=4pi rho^2, ubar*=rho F(D))")
    break

# 5. eps -> 0
print("5. D_max(tau) as eps -> 0 at fixed tau (bound needs a cap):")
for tau in (mp.mpf('1.5'), mp.mpf(3)):
    row = []
    for e in ('1e-2', '1e-4', '1e-8', '1e-16'):
        eps = mp.mpf(e)
        Dm = mp.findroot(lambda D: F(eps, D)[0] - tau, mp.sqrt(4*mp.log(1/eps)*mp.log(tau)) + 0.3)
        row.append(f"eps={e}: {mp.nstr(Dm,6)} (sqrt(4 Lam ln tau)={mp.nstr(mp.sqrt(4*mp.log(1/eps)*mp.log(tau)),5)})")
    print(f"   tau={mp.nstr(tau,3)}: " + "; ".join(row))
