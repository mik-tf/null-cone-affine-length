"""Step 4 of the shear-clock line (referee round 1): numbers behind the revised note.
 A. Elementary affine-parameter bound  L_el = 2 sqrt(tau ln(1/eps)/(1-eps))  versus D_max, the constant-speed
    geodesic bound, and the Hoelder bound from (D.52a) of arXiv 2611.34924 (v1) at its best t in (0,1].
 B. Small-excess law: delta = F - (1-eps) against the parametric form and the leading form D^2/(4 ln(1/(eps D^2))).
 C. Existence curve of Theorem 1: Lambda(alpha) runs from ln(1+lam0) down to 0 (intermediate value theorem).
Run: python3 -I shear_clock_step4.py   (about a minute, mpmath 30 digits)
"""
import mpmath as mp
mp.mp.dps = 30

def Lam_of(al, be): return (be*mp.log(1/(1-al)) - al*mp.log(1+be))/(al+be)
def D_of(al, be):   return 2*mp.sqrt(al*be)/(al+be)*mp.log((1+be)/(1-al))
def lam0_of(al, be): return al*be/((1-al)*(1+be))

def solve_tau(eps, tau):
    """(alpha, beta) for clock tau via the 1-parameter curve of the existence proof (bisection on alpha)."""
    lam0 = tau/eps; Lam = mp.log(1/eps)
    beta_of = lambda al: (lambda r: r/(1-r))(lam0*(1-al)/al)
    f = lambda al: Lam_of(al, beta_of(al)) - Lam
    lo = lam0/(1+lam0)*(1+mp.mpf('1e-25')); hi = 1 - mp.mpf('1e-25')
    # f(lo) > 0, f(hi) < 0
    for _ in range(300):
        mid = (lo+hi)/2
        if f(mid) > 0: lo = mid
        else: hi = mid
    al = (lo+hi)/2
    return al, beta_of(al)

def Dmax(eps, tau):
    al, be = solve_tau(eps, tau); return D_of(al, be)

def L_el(eps, tau): return 2*mp.sqrt(tau*mp.log(1/eps)/(1-eps))

def L_geo(eps, tau):  # constant-speed geodesic of clock tau: tau = eps (e^{k Lam}-1)/k, k = 1 + D^2/(4 Lam^2)
    Lam = mp.log(1/eps)
    k = mp.findroot(lambda k: eps*(mp.e**(k*Lam)-1)/k - tau, 1 + mp.mpf('1e-6') + mp.log(tau)/Lam)
    return 2*Lam*mp.sqrt(k-1)

def L_hoelder(eps, tau):  # ell^2 <= (4/eps) min_t tau^t (1-eps)^(1-t) / t, from (D.52a) with 1-exp(-t Xi(1)) <= 1, D = sqrt2 ell
    f = lambda t: tau**t*(1-eps)**(1-t)/t
    best = min(f(mp.mpf(k)/1000) for k in range(1, 1001))
    return 2*mp.sqrt(best/eps)

eps = mp.mpf('1e-3')
print("A. eps = 1e-3.  tau | D_max | elementary L_el | ratio L_el/D_max | geodesic bound | Hoelder (D.52a)")
for tau in ('1.001', '1.01', '1.1', '1.5', '2', '5', '10', '100', '1000'):
    tau = mp.mpf(tau); d = Dmax(eps, tau)
    print(f"  tau={mp.nstr(tau,5):6s} D_max={mp.nstr(d,6):9s} L_el={mp.nstr(L_el(eps,tau),6):9s} "
          f"L_el/D_max={mp.nstr(L_el(eps,tau)/d,5):8s} geo={mp.nstr(L_geo(eps,tau),6):9s} Hold={mp.nstr(L_hoelder(eps,tau),6)}")
print("   gain over elementary bound at tau=2, eps=1e-3: 1 - D_max/L_el =", mp.nstr(1-Dmax(eps,mp.mpf(2))/L_el(eps,mp.mpf(2)), 4))

print("B. small excess: eps | D | beta | delta | delta/param | delta/(D^2/(4 ln(1/(eps D^2)))) | delta/(D^2/(8 ln(1/(eps D)))) | delta/implicit-B form")
for e in ('1e-1', '1e-3', '1e-6', '1e-12'):
    eps_ = mp.mpf(e); Lam = mp.log(1/eps_)
    for be in ('1e4', '1e7', '1e10'):
        be = mp.mpf(be)
        # alpha from the first equation with lam0 unknown: use Lambda equation to find lam0 for given beta
        g = lambda lam0: Lam_of((1+be)*lam0/((1+be)*lam0+be), be) - Lam
        lam0 = mp.findroot(g, 1/eps_ - 1 + 1/be)
        al = (1+be)*lam0/((1+be)*lam0+be)
        D = D_of(al, be); delta = eps_*lam0 - (1-eps_)
        B = Lam + mp.log(be)
        param = D**2*(B-1)/(4*B**2)
        lead = D**2/(4*mp.log(1/(eps_*D**2)))
        old = D**2/(8*mp.log(1/(eps_*D)))
        Bi = mp.log(1/(eps_*D**2))
        for _ in range(60): Bi = mp.log(1/(eps_*D**2)) + 2*mp.log(Bi) + mp.log(4*(1-eps_))
        impl = D**2*(Bi-1)/(4*Bi**2)
        print(f"  eps={e:6s} D={mp.nstr(D,4):9s} beta={mp.nstr(be,2):6s} delta={mp.nstr(delta,5):11s} "
              f"{mp.nstr(delta/param,6):9s} {mp.nstr(delta/lead,5):8s} {mp.nstr(delta/old,5):8s} {mp.nstr(delta/impl,5)}")

print("C. existence curve, eps=1e-3, tau=2: Lambda(alpha) from ln(1+lam0)=", mp.nstr(mp.log(1+2/eps),6), "to 0; target ln(1/eps)=", mp.nstr(mp.log(1/eps),6))
lam0 = 2/eps
for x in ('1e-12', '1e-6', '1e-3', '0.1', '0.5', '0.9'):
    al = lam0/(1+lam0) + (1-lam0/(1+lam0))*mp.mpf(x) if x not in ('1e-12',) else lam0/(1+lam0)*(1+mp.mpf(x))
    r = lam0*(1-al)/al; be = r/(1-r)
    print(f"   alpha={mp.nstr(al,10)} beta={mp.nstr(be,6)} Lambda={mp.nstr(Lam_of(al,be),8)}")
