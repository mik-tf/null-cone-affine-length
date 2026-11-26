"""Shear-clock deeper (2026-11-26): direction (1) downstream of arXiv 2611.34924 (D.50)-(D.55).
Sharp patch bound  int_P D^2 dA <= 2 A_P Phi(T/A_P),  Phi(tau) = D_max(tau,eps)^2  (D = sqrt2 * hyperbolic length),
versus the paper's (D.52b)  8 e (1-eps) A_P/eps * [1 + log(A*/A_P)].  Also: clock excess does not control Q[h] (D.40).
Run: python3 -I shear_clock_deeper.py   (well under 1 minute)."""
import mpmath as mp
mp.mp.dps = 60

def dmax(tau, eps):
    """Closed form of Theorem 1 (article D): solve for alpha, beta; return D_max (hyperbolic length)."""
    tau, eps = mp.mpf(tau), mp.mpf(eps)
    lam0, Lam = tau/eps, mp.log(1/eps)
    def beta_of(al):
        r = lam0*(1-al)/al
        return r/(1-r)
    def f(al):
        b = beta_of(al)
        return (b*mp.log(1/(1-al)) - al*mp.log(1+b))/(al+b) - Lam
    lo, hi = lam0/(1+lam0)*(1+mp.mpf(10)**-35), 1-mp.mpf(10)**-35
    al = mp.findroot(f, (lo, hi), solver='anderson')
    b = beta_of(al)
    return 2*mp.sqrt(al*b)/(al+b)*mp.log((1+b)/(1-al))

def Phi(tau, eps):
    return dmax(tau, eps)**2

print("[0] closed form check vs article table (eps=1e-3): tau=1.01,2,10,100 ->",
      [mp.nstr(dmax(t, 1e-3), 5) for t in (1.01, 2, 10, 100)], "(article: 0.8119 5.4142 9.2691 12.653)")

print("\n[1] concavity of Phi = D_max^2 in tau (min second difference on a log grid; >0 would break Jensen)")
for eps in (1e-1, 1e-3, 1e-6):
    taus = [(1-eps) + mp.mpf(10)**(k/4.0) for k in range(-24, 25)]
    vals = [Phi(t, eps) for t in taus]
    worst = mp.mpf(-1e9)
    for i in range(1, len(taus)-1):
        t0, t1, t2 = taus[i-1], taus[i], taus[i+1]
        # divided second difference
        s = ((vals[i+1]-vals[i])/(t2-t1) - (vals[i]-vals[i-1])/(t1-t0))/((t2-t0)/2)
        worst = max(worst, s)
    print(f"  eps={eps:g}: max second divided difference = {mp.nstr(worst, 4)} (negative everywhere = concave)")

print("\n[2] large-tau law: Phi(tau)/(4 ln(1/eps) ln tau) and Phi - 4 ln(1/eps) ln tau")
for eps in (1e-3,):
    L = mp.log(1/mp.mpf(eps))
    for t in (1e2, 1e4, 1e8, 1e16):
        p = Phi(t, eps)
        print(f"  eps={eps:g} tau={t:g}: ratio={mp.nstr(p/(4*L*mp.log(t)),6)}  Phi={mp.nstr(p,6)}")

print("\n[2b] proved lower half: Phi(tau) >= 4 ln(1/eps) ln tau for tau > 1 (constant-speed curve); min ratio on a grid")
for eps in (1e-1, 1e-3, 1e-6):
    L = mp.log(1/mp.mpf(eps))
    r = min(Phi(t, eps)/(4*L*mp.log(t)) for t in [mp.mpf(10)**(k/4.0) for k in range(1, 65)])
    print(f"  eps={eps:g}: min Phi/(4 ln(1/eps) ln tau) over tau=10^0.25..10^16 = {mp.nstr(r,5)} (must be >= 1)")

print("\n[3] sharp patch bound vs paper (D.52b), per unit A_P, x = A*/A_P (clock budget T=(1-eps)A*)")
print("    sharp: 2*Phi((1-eps)x);  paper: 8e(1-eps)/eps*(1+ln x);  ratio paper/sharp")
for eps in (1e-1, 1e-2, 1e-3):
    for x in (1e1, 1e3, 1e6, 1e12):
        s = 2*Phi((1-eps)*x, eps)
        pp = 8*mp.e*(1-eps)/eps*(1+mp.log(x))
        print(f"  eps={eps:g} x={x:g}: sharp={mp.nstr(s,5):>8}  paper={mp.nstr(pp,5):>9}  ratio={mp.nstr(pp/s,4)}")

print("\n[4] minimal cell from (D.55) with the sharp constant: |loop area| <= A_P Phi((1-eps)A*/A_P)/(2 pi G)")
print("    smallest A_P/(G hbar) that can enclose hbar, and cells N = A*/A_P, for s* = A*/(4 G hbar)")
for eps in (1e-3,):
    for sstar in (1e2, 1e6, 1e12):
        Astar = 4*mp.mpf(sstar)
        g = lambda y: mp.e**y*Phi((1-eps)*Astar/mp.e**y, eps)/(2*mp.pi) - 1   # y = ln(A_P/(G hbar))
        lo, hi = mp.log(Astar)-80, mp.log(Astar)-mp.mpf('1e-3')   # g(lo)<0, g(hi)>0 assumed; plain bisection
        for _ in range(60):
            mid = (lo+hi)/2
            lo, hi = (mid, hi) if g(mid) < 0 else (lo, mid)
        y = (lo+hi)/2
        AP = mp.e**y
        print(f"  eps={eps:g} s*={sstar:g}: A_P,min/(G hbar)={mp.nstr(AP,4)}  N=A*/A_P={mp.nstr(Astar/AP,4)}"
              f"  N/(s* ln s*)={mp.nstr(Astar/AP/(sstar*mp.log(sstar)),4)}")

print("\n[5] the clock does not control Q[h] (D.40) shear term: |gamma_v|^2 = A2 on [1-eta,1], K = A2*eta, eta = exp(-K/4)/K")
print("    exact nonlinear clock excess tau-(1-eps) vs Q_shear = int v^2 |gamma_v|^2/8 (normalisation of D.40 up to a fixed factor)")
for K in (1, 10, 40, 100, 200):
    K = mp.mpf(K); eta = mp.e**(-K/4)/K; A2 = K/eta
    Xi = lambda v: A2*((v**2-(1-eta)**2)/8)            # Xi = 1/4 int s |gamma_s|^2 ds on [1-eta, v]
    exc = mp.quad(lambda v: mp.e**Xi(v)-1, [1-eta, 1-eta/2, 1-eta/10, 1])
    Q = A2*(1-(1-eta)**3)/24
    print(f"  K={int(K):>3}: clock excess={mp.nstr(exc,4):>9}  Q_shear={mp.nstr(Q,5):>8}  hyperbolic length={mp.nstr(mp.sqrt(A2)*eta,4)}")
print("  -> clock excess -> 0 while Q_shear ~ K/8 grows without bound.")
