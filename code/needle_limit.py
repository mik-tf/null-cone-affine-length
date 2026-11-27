"""Thread 3 read: the needle limit of Theorem 1 of the null-cone affine-length paper.

Question. At fixed cap ratio eps (Lambda = ln(1/eps)), what minimum normalised affine
length tau = F_eps(d) does an expanding generator need to reach hyperbolic excursion d,
as d -> infinity (shape aspect ratio e^(-d) -> 0, a needle)?

Claims checked here (asserts fixed in advance, measured numbers printed):
  (A) series: with alpha = 1 (error exponentially small in u), the equations of Theorem 1
      give beta = Lambda/u + O(u^-2) and D^2 = 4 Lambda u + 4 Lambda (1 - Lambda) + O(1/u),
      u = ln(1/(1-alpha)).
  (B) large-clock law with its constant:
      D_max(tau,eps)^2 = 4 Lambda (ln tau + ln(u/Lambda) + 1) + o(1)
                       = 4 Lambda (ln tau + ln ln tau + 1 - ln Lambda) + o(1)  (tau -> inf).
  (C) inverse (the needle statement):
      F_eps(d) = (4 Lambda^2 / (e d^2)) exp(d^2/(4 Lambda)) (1 + o(1))   (d -> inf).
  (D) constant-speed curve at the same d needs e times more clock in the limit.
  (E) the extremal profile a = c/lambda, rescaled by Lambda/d, tends to 1 uniformly on
      [0, Lambda - delta]; its non-constant part is a layer of width ~ beta ~ Lambda/u at the cut.
  (F) two exact identities (v1.1): x = u - Lambda = T + ln((1+beta)/(alpha beta)) and
      D^2 = 4 (Lambda + ln(1+beta)) x = 4 (u - Lambda)(Lambda + ln(1+beta)); to 1e-40 at 50 digits.
  (G) explicit two-sided bound Psi + r_- < D^2 <= Psi + r_+ for ln tau >= 3 + Lambda (proved in the
      text; this checks the formulas on an 8 x 8 grid, T up to 1e6) and the remainder decomposition.
  (H) R = D^2 - Psi > 0 for 1 < tau < tau_0 = e^3/eps, on the 400-point grid T = k T0/400 in T = ln tau (computed, not proved).
  (I) needle form F = (4 Lambda^2/(e d^2)) e^y Q, Q = e^-rho (1 + (b+rho)/T_F), its bounds
      q(rho_+) <= Q <= q(rho_-) and the expansion Q = 1 + (1 - Lambda/2)/T_F + O(ln^2 T_F/T_F^2).
Run: nice -n 19 python3 -I needle_limit.py   (light; well under a minute)."""
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
E = mp.e


def solve(tau, eps):
    """Theorem 1 closed form, parametrised by beta (alpha explicit from the first equation)."""
    tau, eps = mp.mpf(tau), mp.mpf(eps)
    lam0, Lam, L = tau / eps, mp.log(1 / eps), mp.log(tau / eps)

    def alpha(b):  # 1 - alpha = alpha*beta/(lam0 (1+beta))
        return 1 / (1 + b / (lam0 * (1 + b)))

    def u_of(b):  # ln(1/(1-alpha)) without cancellation
        return L + mp.log(1 + b) - mp.log(alpha(b)) - mp.log(b)

    def f(b):
        a = alpha(b)
        return (b * u_of(b) - a * mp.log(1 + b)) / (a + b) - Lam

    # f -> -Lam as b -> 0 and f -> ln(1+lam0) - Lam > 0 as b -> inf: bracket by bisection on log b
    lo, hi = mp.mpf(10) ** -40, mp.mpf(10) ** 40
    for _ in range(90):
        mid = mp.sqrt(lo * hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    b = mp.findroot(f, mp.sqrt(lo * hi))
    a, u = alpha(b), u_of(b)
    D = 2 * mp.sqrt(a * b) / (a + b) * (u + mp.log(1 + b))
    return dict(alpha=a, beta=b, u=u, D=D, lam0=lam0, Lam=Lam)


def F_of_D(D, eps):
    """Inverse: minimum clock for excursion D (Lemma 2: D_max is strictly increasing in tau)."""
    D, eps = mp.mpf(D), mp.mpf(eps)
    Lam = mp.log(1 / eps)
    # start from the claimed asymptotic form (C); secant on ln tau
    x0 = D ** 2 / (4 * Lam) + mp.log(4 * Lam ** 2 / (mp.e * D ** 2))
    return mp.e ** mp.findroot(lambda x: solve(mp.e ** x, eps)["D"] - D, (x0, x0 + mp.mpf("0.1")), solver="secant")


def tau_const(D, eps):
    """Constant-speed curve a = 2 sqrt(k-1): length 2 Lam sqrt(k-1), clock eps (e^{k Lam}-1)/k."""
    Lam = mp.log(1 / mp.mpf(eps))
    k = 1 + mp.mpf(D) ** 2 / (4 * Lam ** 2)
    return eps * (mp.e ** (k * Lam) - 1) / k


print("(A) sympy series of Theorem 1's equations at alpha = 1, in w = 1/u")
w, Lm, b2, b3 = sp.symbols("w Lambda b2 b3", positive=True)
beta = Lm * w + b2 * w ** 2 + b3 * w ** 3
eq2 = sp.series(Lm * (1 + beta) * w - beta + w * sp.log(1 + beta), w, 0, 4).removeO()
sol = sp.solve([sp.expand(eq2).coeff(w, 2), sp.expand(eq2).coeff(w, 3)], [b2, b3], dict=True)[0]
beta_s = sp.simplify(beta.subs(sol))
D2 = 4 * beta_s * (1 / w + sp.log(1 + beta_s)) ** 2 / (1 + beta_s) ** 2
D2s = sp.expand(sp.series(D2, w, 0, 2).removeO())
print("  beta =", sp.factor(beta_s))
print("  D^2  =", sp.collect(D2s, w))
lead = sp.simplify(D2s.coeff(w, -1)); const = sp.simplify(D2s.coeff(w, 0))
assert sp.simplify(lead - 4 * Lm) == 0 and sp.simplify(const - 4 * Lm * (1 - Lm)) == 0
print("  PASS: D^2 = 4 Lambda u + 4 Lambda (1 - Lambda) + O(1/u)")

print("\n(0) reproduce the paper's large-clock ratios (eps = 1e-3; paper: 1.258 1.180 1.120 1.076)")
r = [solve(mp.mpf(10) ** k, 1e-3)["D"] ** 2 / (4 * mp.log(1000) * mp.log(mp.mpf(10) ** k)) for k in (2, 4, 8, 16)]
print("  ", [mp.nstr(x, 4) for x in r])
assert [round(float(x), 3) for x in r] == [1.258, 1.180, 1.120, 1.076]

print("\n(B) residuals R1 = D^2 - 4L(ln tau + ln(u/L) + 1),  R2 = D^2 - 4L(ln tau + ln ln tau + 1 - ln L)")
print("    eps     | log10 tau | beta        | D          | R1          | R2          | R2*ln(tau)/ln ln tau")
for eps in (mp.mpf("0.1"), mp.mpf("1e-3"), mp.mpf("1e-6")):
    Lam = mp.log(1 / eps)
    last = None
    for k in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096, 16384):
        tau = mp.mpf(10) ** k
        s = solve(tau, eps)
        lt = mp.log(tau)
        R1 = s["D"] ** 2 - 4 * Lam * (lt + mp.log(s["u"] / Lam) + 1)
        R2 = s["D"] ** 2 - 4 * Lam * (lt + mp.log(lt) + 1 - mp.log(Lam))
        print(f"    {mp.nstr(eps,2):7s} | {k:9d} | {mp.nstr(s['beta'],4):11s} | {mp.nstr(s['D'],8):10s} |"
              f" {mp.nstr(R1,4):11s} | {mp.nstr(R2,4):11s} | {mp.nstr(R2*lt/mp.log(lt),4)}")
        last = (R1, R2)
    assert abs(last[0]) < mp.mpf("0.05") * (1 + Lam) ** 2 / 10 and abs(last[1]) < mp.mpf("0.05") * (1 + Lam) ** 2
print("  PASS: both residuals tend to 0 (R1 like (1+L)^2/ln tau, R2 like ln ln tau / ln tau)")

print("\n(C) needle form: F_eps(d) / [ (4 L^2/(e d^2)) exp(d^2/(4L)) ]   and   (D) tau_const(d)/F_eps(d) -> e")
print("    eps     | d       | F_eps(d)      | ratio (C) | tau_const/F | /e")
for eps in (mp.mpf("0.1"), mp.mpf("1e-3")):
    Lam = mp.log(1 / eps)
    rows = []
    for mult in (2, 4, 8, 16, 32, 64):
        d = mult * 2 * Lam
        Fd = F_of_D(d, eps)
        pred = 4 * Lam ** 2 / (E * d ** 2) * mp.e ** (d ** 2 / (4 * Lam))
        rc = tau_const(d, eps) / Fd
        rows.append((Fd / pred, rc / E))
        print(f"    {mp.nstr(eps,2):7s} | {mp.nstr(d,5):7s} | {mp.nstr(Fd,6):13s} | {mp.nstr(Fd/pred,6):9s} |"
              f" {mp.nstr(rc,6):11s} | {mp.nstr(rc/E,6)}")
    assert abs(rows[-1][0] - 1) < mp.mpf("0.02") and abs(rows[-1][1] - 1) < mp.mpf("0.02")
    assert abs(rows[-1][0] - 1) < abs(rows[0][0] - 1) and abs(rows[-1][1] - 1) < abs(rows[0][1] - 1)
print("  PASS: ratio (C) -> 1, tau_const/F -> e")

print("\n(E) extremal profile a(t) = c/lambda(t); max |a Lambda/d - 1| on t in [0, Lambda - delta]")
print("    eps 1e-3; parametric t(lambda) of Theorem 1; layer width beta")
eps = mp.mpf("1e-3"); Lam = mp.log(1 / eps)
prev = None
for k in (4, 16, 64, 256, 1024):
    s = solve(mp.mpf(10) ** k, eps)
    a_, b_, lam0, d = s["alpha"], s["beta"], s["lam0"], s["D"]
    rp, rm = lam0 / a_, lam0 / b_
    c = 2 * mp.sqrt(rp * rm)

    def t_of_y(y):  # y = r_+ - lambda > 0, avoids cancellation near lambda = lam0
        return Lam - (rp * mp.log(rp / y) - rm * mp.log1p((rp - y) / rm)) / (rp + rm)

    out = {}
    for delta in (mp.mpf("0.1") * Lam, mp.mpf("0.01") * Lam):
        # t increases as y increases (lambda decreases); bisect ln y for t = Lam - delta
        lo, hi = mp.log(b_ / (1 + b_)), mp.log(rp)  # r_+ - lam0 = beta/(1+beta) exactly
        for _ in range(200):
            mid = (lo + hi) / 2
            if t_of_y(mp.e ** mid) > Lam - delta:
                hi = mid
            else:
                lo = mid
        lam_d = rp - mp.e ** ((lo + hi) / 2)
        # a = c/lambda increases along t, so its extremes on [0, Lam - delta] are at the two ends
        out[delta] = max(abs(c / lam0 * Lam / d - 1), abs(c / lam_d * Lam / d - 1))
    print(f"    log10 tau={k:5d}  beta={mp.nstr(b_,4):9s}  delta=0.1L: {mp.nstr(out[mp.mpf('0.1')*Lam],4):10s}"
          f"  delta=0.01L: {mp.nstr(out[mp.mpf('0.01')*Lam],4)}")
    prev = out
assert prev[mp.mpf("0.01") * Lam] < mp.mpf("0.01")
print("  PASS: rescaled extremal speed -> 1 uniformly away from the cut; no concentration")


# ---------------------------------------------------------------------------------------------
# v1.1 additions (F) to (I).  Notation: lam0 = tau/eps, Lam = ln(1/eps), T = ln tau, u = ln(1/(1-alpha)),
# l = ln(1+beta), x = u - Lam, m = Lam + l, Psi(T) = 4 Lam (T + ln T + 1 - ln Lam), R = D^2 - Psi(T).
# ---------------------------------------------------------------------------------------------
DPS = 50
EPS_GRID = [mp.mpf(s_) for s_ in ("0.99", "0.9", "0.5", "0.1", "1e-3", "1e-6", "1e-12", "1e-30")]
TOL = mp.mpf("1e-40")


def solve_T(T, eps):
    """Root beta of equation 2 given T = ln tau (call inside mp.workdps(DPS)).  Equation 1 gives alpha and
    1-alpha = beta/(lam0 (1+beta)+beta) exactly, so u = ln(lam0 (1+beta)+beta) - ln beta needs no cancellation."""
    T, eps = mp.mpf(T), mp.mpf(eps)
    Lam = mp.log(1 / eps)
    lam0 = mp.e ** T / eps

    def alpha(b):
        return lam0 * (1 + b) / (lam0 * (1 + b) + b)

    def u_of(b):
        return mp.log(lam0 * (1 + b) + b) - mp.log(b)

    def f(b):  # equation 2 in the form (b u - a l)/(a + b) - Lam
        a = alpha(b)
        return (b * u_of(b) - a * mp.log1p(b)) / (a + b) - Lam

    lo, hi = mp.mpf(10) ** -60, mp.mpf(10) ** 60
    for _ in range(260):  # bisection on ln b: 2 ln(1e60) / 2^260 is far below 1e-50
        mid = mp.sqrt(lo * hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    b = mp.findroot(f, mp.sqrt(lo * hi))
    a, u = alpha(b), u_of(b)
    l = mp.log1p(b)
    D = 2 * mp.sqrt(a * b) / (a + b) * (u + l)  # Theorem 1
    nla = mp.log1p(b / (lam0 * (1 + b)))  # -ln alpha, no cancellation
    oma = b / (lam0 * (1 + b) + b)  # 1 - alpha
    return dict(alpha=a, beta=b, u=u, l=l, D=D, lam0=lam0, Lam=Lam, T=T, eps=eps, nla=nla, oma=oma)


def Psi(T, Lam):
    return 4 * Lam * (T + mp.log(T) + 1 - mp.log(Lam))


def bounds(T, eps):
    """r_-(T), r_+(T), sbar as in the Proposition."""
    T, eps = mp.mpf(T), mp.mpf(eps)
    Lam = mp.log(1 / eps)
    de = eps / mp.e ** T  # eps/tau
    sbar = (Lam + T) * (mp.log(1 + T / Lam) + 2 * de) / (Lam + T - 1)
    first = 4 * Lam * mp.log(1 + mp.log(1 + (T - 1) / Lam) / T)
    second = max(mp.mpf(0), 2 * Lam ** 2 / ((1 + de) * (T + sbar)) - (4 + 2 * Lam) * Lam ** 2 / (T - 1) ** 2)
    rm = first + second
    rp = 4 * Lam * sbar / T + 2 * Lam ** 2 / (T - 1) + (4 * Lam ** 3 / 3 + 2 * Lam) / (T - 1) ** 2 + 8 * Lam * de
    return rm, rp, sbar


def T_grid(T0, n=8):
    """n points, geometric from T0 to 1e6 inclusive."""
    return [T0 * (mp.mpf(10) ** 6 / T0) ** (mp.mpf(k) / (n - 1)) for k in range(n)]


print("\n(F) exact identities at 50 digits: eps in {0.99 ... 1e-30}, T from T0 = 3 + Lambda to 1e6")
print("    I1: x - [T + ln((1+b)/(a b))]   I2: D^2 - 4 m x   I3: D - 2 sqrt(a/b) m   I4: eq.2  b x - a m   (residuals)")
print("    eps     | T        | beta        | D^2          | |I1|     | |I2|     | |I3|     | |I4|")
worst = mp.mpf(0)
with mp.workdps(DPS):
    for eps in EPS_GRID:
        Lam = mp.log(1 / eps)
        for T in T_grid(3 + Lam):
            s = solve_T(T, eps)
            a, b, l, u, D, Lam_ = s["alpha"], s["beta"], s["l"], s["u"], s["D"], s["Lam"]
            x, m = u - Lam_, Lam_ + l
            r1 = abs(x - (T + mp.log((1 + b) / (a * b))))
            r2 = abs(D ** 2 - 4 * m * x)
            r3 = abs(D - 2 * mp.sqrt(a / b) * m)
            r4 = abs(b * x - a * m)
            worst = max(worst, r1, r2, r3, r4)
            print(f"    {mp.nstr(eps,2):7s} | {mp.nstr(T,6):8s} | {mp.nstr(b,5):11s} | {mp.nstr(D**2,10):12s} |"
                  f" {mp.nstr(r1,2):8s} | {mp.nstr(r2,2):8s} | {mp.nstr(r3,2):8s} | {mp.nstr(r4,2)}")
            assert r1 < TOL and r2 < TOL and r3 < TOL and r4 < TOL
            slack = 1 + mp.mpf("1e-45")  # relative rounding slack; the strict inequalities are proved in the text
            assert 0 < s["nla"] <= mp.log1p(1 / s["lam0"]) * slack and mp.log1p(1 / s["lam0"]) <= eps / mp.e ** T * slack  # 0 < -ln alpha <= ln(1+1/lam0) < eps/tau
            assert s["oma"] <= b / s["lam0"] * slack
    print("    direct check (no use of the closed form for 1-alpha): alpha from eq.1 at T = T0, u = -ln(1-alpha)")
    for eps in EPS_GRID:
        Lam = mp.log(1 / eps)
        T = 3 + Lam
        with mp.workdps(DPS + 10 + int((T + Lam) / 2.302)):  # enough digits to survive the cancellation in 1 - alpha
            s = solve_T(T, eps)
            a, b = s["alpha"], s["beta"]
            lam0 = s["lam0"]
            u_dir = -mp.log(1 - a)
            res1 = abs(lam0 - a * b / ((1 - a) * (1 + b)))  # eq. 1 of Theorem 1
            m = Lam + mp.log1p(b)
            res2 = abs(b * (u_dir - Lam) - a * m)
            res3 = abs(s["D"] ** 2 - 4 * (u_dir - Lam) * m)
            assert res1 < TOL * lam0 and res2 < TOL and res3 < TOL
        print(f"      eps={mp.nstr(eps,2):7s}  eq.1 rel {mp.nstr(res1/lam0,2)}  eq.2 {mp.nstr(res2,2)}  D^2 {mp.nstr(res3,2)}")
print(f"  worst residual {mp.nstr(worst,3)}")
print("  PASS: both identities hold to 1e-40 at 50 digits; 0 < -ln alpha <= ln(1+1/lam0) < eps/tau and 1-alpha <= beta/lam0")

print("\n(G) r_- < R <= r_+ on the 64-point grid (T >= T0 = 3 + Lambda, i.e. tau >= e^3/eps); R = D^2 - Psi(T)")
print("    eps     | T        | R            | r_-          | r_+          | R/r_+   | decomposition residual")
minratio = mp.mpf(10)
with mp.workdps(DPS):
    for eps in EPS_GRID:
        Lam = mp.log(1 / eps)
        for T in T_grid(3 + Lam):
            s = solve_T(T, eps)
            a, b, l = s["alpha"], s["beta"], s["l"]
            R = s["D"] ** 2 - Psi(T, Lam)
            rm, rp, sbar = bounds(T, eps)
            # exact remainder decomposition
            g = (Lam + l) * l / b - Lam + Lam * l - Lam * mp.log(1 + l / Lam)
            x = s["u"] - Lam
            oma, nla = s["oma"], s["nla"]
            dec = 4 * Lam * mp.log(x / T) + 4 * g + 8 * Lam * nla - 4 * oma * (Lam + l) * l / b
            # intermediate claims
            assert x > T and b < Lam / (T - 1) and b < Lam / (2 + Lam) and g > 0
            assert Lam * b / 2 - b ** 2 * (1 + Lam / 2) <= g <= Lam * b / 2 + b ** 2 * (Lam / 3 + 1 / (2 * Lam))
            a_terms = 8 * Lam * nla - 4 * oma * (Lam + l) * l / b
            assert 4 * oma * (Lam - b) <= a_terms <= 8 * Lam * eps / mp.e ** T
            sx = x - T
            assert mp.log(1 + (T - 1) / Lam) < sx <= sbar and b >= a * Lam / (T + sbar) and a >= 1 / (1 + eps / mp.e ** T)
            assert abs(R - dec) < TOL and rm > 0 and rm < R <= rp
            minratio = min(minratio, R / rp)
            print(f"    {mp.nstr(eps,2):7s} | {mp.nstr(T,6):8s} | {mp.nstr(R,8):12s} | {mp.nstr(rm,8):12s} |"
                  f" {mp.nstr(rp,8):12s} | {mp.nstr(R/rp,5):7s} | {mp.nstr(abs(R-dec),2)}")
print(f"  smallest R/r_+ on the grid: {mp.nstr(minratio,3)} (about 0.38 to 0.58 at T0 depending on eps; 0.99994 to 0.99999 at T = 1e6; ratio never exceeds 1)")
print("  PASS: r_- < R <= r_+ at all 64 points, and every intermediate claim of the proof holds at each point")

print("\n(H) R = D^2 - Psi(T) > 0 for 1 < tau < tau_0 (R is defined only for tau > 1, i.e. T = ln tau > 0).")
print("    Grid in T = ln tau: T_k = k T0/400, k = 1..400, T0 = 3 + Lambda (so tau from e^(T0/400) to tau_0 = e^T0). Computed, not proved.")
print("    R > 0 is trivial where Psi(T) <= 0 (small T, from the ln T term); 'nontriv' counts grid points with Psi > 0.")
print("    eps     | T0       | min R on grid | at T      | nontriv")
with mp.workdps(DPS):
    for eps in [mp.mpf(s_) for s_ in ("0.9", "0.1", "1e-3", "1e-6", "1e-12", "1e-30")]:
        Lam = mp.log(1 / eps)
        T0 = 3 + Lam
        pts = [T0 * mp.mpf(k) / 400 for k in range(1, 401)]
        best = None
        nontriv = 0
        for T in pts:
            s = solve_T(T, eps)
            R = s["D"] ** 2 - Psi(T, Lam)
            assert R > 0
            if Psi(T, Lam) > 0:
                nontriv += 1
            if best is None or R < best[0]:
                best = (R, T)
        print(f"    {mp.nstr(eps,2):7s} | {mp.nstr(T0,6):8s} | {mp.nstr(best[0],8):13s} | {mp.nstr(best[1],6):9s} | {nontriv}")
print("  PASS (computed): R > 0 at every grid point T_k in (0, T0] (tau in (1, tau_0]); a numerical observation below tau_0, not a theorem")

print("\n(I) needle form: F = (4 Lam^2/(e d^2)) e^y Q,  Q = e^-rho (1 + (b+rho)/T_F),  q(rho_+) <= Q <= q(rho_-)")
print("    eps     | d        | T_F       | Q            | q(rho_+)     | q(rho_-)     | 1+(1-L/2)/T_F | (Q-that) T_F^2/ln^2 T_F")
scoping = {}
with mp.workdps(DPS):
    for eps, ds in ((mp.mpf("0.1"), (295, 600, 2000)), (mp.mpf("1e-3"), (884, 1500, 4000)),
                    (mp.mpf("1e-12"), (700, 2000, 6000))):
        Lam = mp.log(1 / eps)
        d0 = solve_T(3 + Lam, eps)["D"]
        for d in ds:
            d = mp.mpf(d)
            assert d >= d0
            y = d ** 2 / (4 * Lam)
            x0 = y + mp.log(Lam / (mp.e * y))
            TF = mp.findroot(lambda t: solve_T(t, eps)["D"] - d, (x0, x0 + mp.mpf("0.1")), solver="secant", tol=mp.mpf("1e-90"))
            F = mp.e ** TF
            R = d ** 2 - Psi(TF, Lam)
            rho = R / (4 * Lam)
            bb = mp.log(TF) + 1 - mp.log(Lam)
            q = lambda r: mp.e ** (-r) * (1 + (bb + r) / TF)
            Q = q(rho)
            Q_direct = F / (4 * Lam ** 2 / (mp.e * d ** 2) * mp.e ** y)
            rm, rp, _ = bounds(TF, eps)
            qlo, qhi = q(rp / (4 * Lam)), q(rm / (4 * Lam))
            lead = 1 + (1 - Lam / 2) / TF
            ex = (Q - lead) * TF ** 2 / mp.log(TF) ** 2
            # identity, monotonicity, bounds, expansion, constant-speed ratio
            assert abs(Q - Q_direct) < TOL * Q
            assert qlo <= Q <= qhi
            assert abs(ex) < 2 * (1 + Lam) ** 3
            k = 1 + d ** 2 / (4 * Lam ** 2)
            tau_c = eps * (mp.e ** (k * Lam) - 1) / k
            tau_c2 = (4 * Lam ** 2 / (4 * Lam ** 2 + d ** 2)) * (mp.e ** y - eps)
            ratio = mp.e * (1 - eps * mp.e ** (-y)) / ((1 + 4 * Lam ** 2 / d ** 2) * Q)
            assert abs(tau_c - tau_c2) < TOL * tau_c and abs(tau_c / F - ratio) < TOL * ratio
            scoping[(float(eps), int(d))] = Q
            print(f"    {mp.nstr(eps,2):7s} | {mp.nstr(d,5):8s} | {mp.nstr(TF,8):9s} | {mp.nstr(Q,10):12s} | {mp.nstr(qlo,10):12s} |"
                  f" {mp.nstr(qhi,10):12s} | {mp.nstr(lead,10):13s} | {mp.nstr(ex,4)}")
    # Q decreasing in rho (derivative sign) at the sampled points
    for TF_ in (mp.mpf(10), mp.mpf(100), mp.mpf(1e4)):
        for rho_ in (mp.mpf(0), mp.mpf("0.5"), mp.mpf(3)):
            bb = mp.log(TF_) + 1 - mp.log(Lam)
            dq = mp.e ** (-rho_) * (-(TF_ + mp.log(TF_) - mp.log(Lam) + rho_) / TF_)
            assert dq < 0
assert round(float(scoping[(0.1, 295)]), 6) == 0.999984 and round(float(scoping[(1e-3, 884)]), 5) == 0.99991
print("  scoping check (C): eps=0.1, d=295 gives", mp.nstr(scoping[(0.1, 295)], 6), "; eps=1e-3, d=884 gives", mp.nstr(scoping[(1e-3, 884)], 5))
print("  PASS: exact needle form, Q monotone in rho, q(rho_+) <= Q <= q(rho_-), expansion Q = 1 + (1 - Lam/2)/T_F + O(ln^2 T_F/T_F^2),")
print("        tau_c/F = e (1 - eps e^-y)/((1 + 4 Lam^2/d^2) Q) to 1e-40")
print("\n(J) ratio limit of (D.52b) to the sharp patch bound: 4 e (1-eps)(1 + ln(A*/A_P)) / (eps Phi(x)), x = (1-eps) A*/A_P, -> e (1-eps)/(eps Lam)")
print("    Here Phi(tau) = D^2 (the grid check uses Phi in place of Phi-hat; they are asymptotically equal by the majorant lemma, Phi-hat/Phi -> 1).")
print("    eps     | limit      | T=ln x = 1e2  | 1e4        | 1e6        | 1e9        | rel.err at 1e9")
with mp.workdps(DPS):
    for eps in [mp.mpf(s_) for s_ in ("0.1", "0.01", "1e-3")]:
        Lam = mp.log(1 / eps)
        lim = mp.e * (1 - eps) / (eps * Lam)
        vals = []
        for Tx in (mp.mpf(10) ** 2, mp.mpf(10) ** 4, mp.mpf(10) ** 6, mp.mpf(10) ** 9):
            s = solve_T(Tx, eps)
            Lg = Tx - mp.log(1 - eps)  # ln(A*/A_P) = ln x - ln(1-eps)
            vals.append(4 * mp.e * (1 - eps) * (1 + Lg) / (eps * s["D"] ** 2))
        err = abs(vals[-1] / lim - 1)
        assert err < mp.mpf("1e-6") and abs(vals[2] / lim - 1) < mp.mpf("1e-3")
        assert all(abs(vals[i + 1] / lim - 1) < abs(vals[i] / lim - 1) for i in range(3))  # converging
        print(f"    {mp.nstr(eps,2):7s} | {mp.nstr(lim,7):10s} | " + " | ".join(f"{mp.nstr(v,7):10s}" for v in vals) + f" | {mp.nstr(err,2)}")
print("  majorant lemma check at eps = 1e-3 (Psi~ + C >= Phi on a grid, Psi~ = Psi for tau >= tau_0 and the tangent line at tau_0 below):")
with mp.workdps(DPS):
    eps = mp.mpf("1e-3")
    Lam = mp.log(1 / eps)
    T0 = 3 + Lam
    tau0 = mp.e ** T0
    slope = 4 * Lam * (1 + 1 / T0) / tau0                       # d Psi / d tau at tau_0
    L = lambda tau: Psi(T0, Lam) + slope * (tau - tau0)           # tangent line
    Ps = lambda tau: Psi(mp.log(tau), Lam) if tau >= tau0 else L(tau)
    Ts_hi = [T0 * (mp.mpf(10) ** 6 / T0) ** (mp.mpf(k) / 200) for k in range(201)]
    Cp = max(bounds(T, eps)[1] for T in Ts_hi)
    Phi0 = solve_T(T0, eps)["D"] ** 2
    C = max(Cp, Phi0 - L(1 - eps))
    # concavity of Psi~ in tau: second differences on a geometric tau grid, tau from 1-eps to e^50
    taus = [1 - eps] + [mp.e ** (mp.mpf(k) / 4) for k in range(1, 201)]
    taus = sorted(set(taus + [tau0]))
    ok_conc = all(((Ps(taus[i + 1]) - Ps(taus[i])) / (taus[i + 1] - taus[i]) >= (Ps(taus[i + 2]) - Ps(taus[i + 1])) / (taus[i + 2] - taus[i + 1]) - mp.mpf("1e-30"))
                  for i in range(len(taus) - 2))
    worst = None
    for T in [mp.mpf(k) * T0 / 100 for k in range(1, 100)] + Ts_hi[:60]:
        tau = mp.e ** T
        gap = Ps(tau) + C - solve_T(T, eps)["D"] ** 2
        assert gap >= 0
        worst = gap if worst is None else min(worst, gap)
    assert ok_conc
    print(f"    C = {mp.nstr(C,6)} (sup r_+ on grid {mp.nstr(Cp,5)}), Psi~ concave on the tau grid, min (Psi~ + C - Phi) on the grid = {mp.nstr(worst,5)} >= 0")
print("  PASS: the ratio tends to e(1-eps)/(eps Lam) (relative error below 1e-6 at ln x = 1e9, decreasing along the sequence), and the concave majorant Psi~ + C of Phi exists")

print("\nALL PASS")
