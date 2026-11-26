"""Kill check for N1 (new-angles-2611.md): shear term Xi of arXiv 2611.34924
(eq 3.26c) as hyperbolic-plane energy, and the Cauchy-Schwarz clock bound.
Paper conventions: Xi_y(v) = (1/8) int_eps^v s (d_s qb)^AB (d_s qb)_AB ds,
indices raised with qb (unimodular), so (dqb)^AB (dqb)_AB = tr((qb^-1 qb')^2).
qb parametrised as in (4.21a); g = diag(1, e^{2 z1}) = (1/2) tr(qb^-1 dqb qb^-1 dqb).
Run: python3 -I shear_clock_check.py   (about 10 s)
"""
import sympy as sp
import numpy as np
trap = getattr(np, 'trapezoid', None) or np.trapz
from scipy.optimize import minimize

s = sp.symbols('s', positive=True)
z1 = sp.Function('z1')(s); z2 = sp.Function('z2')(s)

# ---- 1. exact identity -------------------------------------------------
q = sp.Matrix([[sp.exp(z1), sp.exp(z1)*z2],
               [sp.exp(z1)*z2, sp.exp(-z1) + sp.exp(z1)*z2**2]])
assert sp.simplify(q.det() - 1) == 0, "not unimodular"
qi = q.inv(); dq = q.diff(s)
shear = sp.simplify((qi*dq*qi*dq).trace())      # (dq)^AB (dq)_AB
target = 2*(sp.diff(z1, s)**2 + sp.exp(2*z1)*sp.diff(z2, s)**2)  # 2 g(gamma',gamma'), g curvature -1
assert sp.simplify(shear - target) == 0
print("[1a] tr((q^-1 q')^2) = 2 (z1'^2 + e^{2 z1} z2'^2): exact")

# metric g = (1/2) tr(...) is the hyperbolic plane: y = e^{-z1}, x = z2 -> (dx^2+dy^2)/y^2
x, y = sp.symbols('x y', positive=True)
yy = sp.exp(-z1); xx = z2
uhp = (sp.diff(xx, s)**2 + sp.diff(yy, s)**2)/yy**2
assert sp.simplify(uhp - (sp.diff(z1, s)**2 + sp.exp(2*z1)*sp.diff(z2, s)**2)) == 0
print("[1b] g = diag(1,e^{2 z1}) equals upper half plane (dx^2+dy^2)/y^2 (K=-1): exact")

# hyperboloid form: qb=[[a,b],[b,c]], ac-b^2=1, X=(a+c)/2, Y=(a-c)/2, Z=b
a_, b_, c_ = [sp.Function(n)(s) for n in 'abc']
Xf, Yf, Zf = (a_+c_)/2, (a_-c_)/2, b_
Q = sp.Matrix([[a_, b_], [b_, c_]])
det1 = {c_: (1 + b_**2)/a_}
Qs = Q.subs(det1).doit()
lhs = sp.simplify((Qs.inv()*Qs.diff(s)*Qs.inv()*Qs.diff(s)).trace()/2)
rhs = (sp.diff(Yf, s)**2 + sp.diff(Zf, s)**2 - sp.diff(Xf, s)**2).subs(det1).doit()
assert sp.simplify(lhs - rhs) == 0
print("[1c] (1/2)tr((q^-1 q')^2) = Y'^2+Z'^2-X'^2 on ac-b^2=1: exact (hyperboloid, K=-1)")

# [1e] Section 2 setup: q = rho^2 qb, rho = (det q)^(1/4); Raychaudhuri <=> rho'' = -(1/8) tr((qb^-1 qb')^2) rho
rho = sp.Function('rho', positive=True)(s)
qf = rho**2*q
assert sp.simplify(sp.sqrt(qf.det()) - rho**2) == 0 or sp.simplify(qf.det() - rho**4) == 0
th = sp.simplify((qf.inv()*qf.diff(s)).trace()/2)            # expansion = (1/2) tr(q^-1 q')
assert sp.simplify(th - 2*sp.diff(rho, s)/rho) == 0
Mx = qf.inv()*qf.diff(s); sig = Mx/2 - th*sp.eye(2)/2        # mixed shear tensor
sigsq = sp.simplify((sig*sig).trace())
assert sp.simplify(sigsq - (qi*dq*qi*dq).trace()/4) == 0
ray = sp.simplify(sp.diff(th, s) + th**2/2 + sigsq)
assert sp.simplify(ray - (2*sp.diff(rho, s, 2)/rho + (qi*dq*qi*dq).trace()/4)) == 0
print("[1e] rho=(det q)^(1/4), theta=2 rho'/rho, sigma^2=(1/4)tr((qb^-1 qb')^2); Raychaudhuri <=> rho''=-(1/8)tr(.) rho: exact")

# measure: s ds = ds in w=ln s with weight: s*|gamma_s|^2 ds = |gamma_w|^2 dw
w = sp.symbols('w', real=True)
G = sp.Function('G')
print("[1d] measure: int s |g_s|^2 ds = int |g_w|^2 dw with w = ln s (chain rule):",
      sp.simplify(sp.diff(G(sp.log(s)), s)**2*s**2 - sp.Subs(sp.diff(G(w), w)**2, w, sp.log(s)).doit()) == 0)
print("     Hence Xi(v) = (1/8) int tr(..) dw = (1/4) int |gamma_w|_g^2 dw = E[gamma]/4, g of curvature -1.")
print("     The paper's own (4.22) states the same with factor 1/4.  Candidate's g=(1/8)tr (K=-4) gives Xi = int |g'|^2 dw: same.")

# ---- 2. bound: Xi(v) >= d(v)^2 / (4 ln(v/eps)),  equality iff constant-speed geodesic in w
def hdist(p, r):  # upper half plane points (x,y)
    return np.arccosh(1 + ((p[0]-r[0])**2 + (p[1]-r[1])**2)/(2*p[1]*r[1]))

def curve_stats(wg, z1v, z2v):
    """Xi(w) cumulative, d(w) from start, on grid wg."""
    xs, ys = z2v, np.exp(-z1v)
    dz1 = np.gradient(z1v, wg); dz2 = np.gradient(z2v, wg)
    speed2 = dz1**2 + np.exp(2*z1v)*dz2**2
    Xi = np.concatenate([[0], np.cumsum(0.5*(speed2[1:]+speed2[:-1])*np.diff(wg))])/4
    d = np.array([hdist((xs[0], ys[0]), (xs[i], ys[i])) for i in range(len(wg))])
    return Xi, d

eps = 1e-3
wg = np.linspace(np.log(eps), 0, 4001)
L = wg - wg[0]; v = np.exp(wg)

def report(name, z1v, z2v):
    Xi, d = curve_stats(wg, z1v, z2v)
    m = L > 1e-9
    bound = d[m]**2/(4*L[m])
    slack = np.min(Xi[m] - bound)
    ratio = np.max(bound/np.maximum(Xi[m], 1e-300)) 
    # clock integral int exp(Xi) dv vs int exp(bound) dv  (rho_S = 1)
    U = trap(np.exp(Xi)*v, wg)
    Ub = trap(np.exp(np.where(m, d**2/(4*np.maximum(L,1e-12)), 0))*v, wg)
    print(f"  {name:10s} min(Xi-bound)={slack:+.2e}  max(bound/Xi)={ratio:.6f}  u*={U:.6f} >= bound-u*={Ub:.6f}: {U >= Ub*(1-1e-6)}")
    return slack >= -1e-5 and U >= Ub*(1-1e-6)

rng = np.random.default_rng(2611)
ok = True
print("[2] numerics on 5 random smooth curves (eps=1e-3, 4001 grid pts) and one geodesic:")
for k in range(5):
    th = rng.uniform(0, 2*np.pi, size=(2, 4)); amp = rng.uniform(0.1, 0.5, size=(2, 4))
    f = lambda t, i: sum(amp[i, j]*np.sin((j+1)*t*0.4 + th[i, j]) for j in range(4))
    t = wg - wg[0]
    ok &= report(f"random{k+1}", f(t, 0), f(t, 1))
# geodesic: vertical line x=0, y = e^{-z1}, z1 = D*L/Ltot (hyperbolic arclength = z1), constant speed in w
D = 2.0
ok &= report("geodesic", D*L/L[-1], np.zeros_like(L))
# a non-vertical geodesic (semicircle), constant speed in w: x=tanh-type param
sig = 2.0*L/L[-1] - 1.0     # arclength parameter from -1 to 1
xs, ys = np.tanh(sig), 1/np.cosh(sig)   # unit semicircle, hyperbolic arclength = sig
ok &= report("semicircle", -np.log(ys), xs)
Xi, d = curve_stats(wg, -np.log(ys), xs)
print(f"     geodesic equality: max|Xi - d^2/(4L)| = {np.max(np.abs(Xi[1:]-d[1:]**2/(4*L[1:]))):.2e}")
print("     all inequality checks pass:", ok)

# ---- 3. the 1D problem: min u* = int e^{Xi(w)} e^w dw given total distance D (speed a(w), Xi'=a^2/4)
def ustar(a, dw, ev):
    Xi = np.concatenate([[0], np.cumsum(0.25*0.5*(a[1:]**2 + a[:-1]**2)*dw)])
    return trap(np.exp(Xi)*ev, dx=dw)   # trapezoid rule (the rectangle sum was 0.9% high at D=1)
n = 400; wg2 = np.linspace(np.log(eps), 0, n); dw = wg2[1]-wg2[0]; ev = np.exp(wg2)
Ltot = wg2[-1]-wg2[0]
def solve(D):
    cons = {'type': 'eq', 'fun': lambda a: np.sum(a)*dw - D}
    a0 = np.full(n, D/Ltot)
    r = minimize(lambda a: ustar(a, dw, ev), a0, constraints=[cons], bounds=[(0, None)]*n, method='SLSQP', options={'maxiter': 300})
    return r, a0
print("[3] 1D problem min u*[a] s.t. int a dw = D  (u* uses rho_S=1):")
for D in (1.0, 3.0, 6.0):
    r, a0 = solve(D)
    Uopt = r.fun; Ugeo = ustar(a0, dw, ev)
    Xg = np.concatenate([[0], np.cumsum(0.25*0.5*(r.x[1:]**2 + r.x[:-1]**2)*dw)])
    Lc = wg2 - wg2[0]; dist = np.concatenate([[0], np.cumsum(0.5*(r.x[1:]+r.x[:-1])*dw)])
    cs_ok = np.all(Xg[1:] >= dist[1:]**2/(4*Lc[1:]) - 1e-9)
    print(f"  D={D}: constant-speed geodesic u*={Ugeo:.6f}, optimum over speed profiles u*={Uopt:.6f}, "
          f"optimum a(first)/a(last)={r.x[0]/max(r.x[-1],1e-12):.2f}, CS holds along optimum: {cs_ok}")
