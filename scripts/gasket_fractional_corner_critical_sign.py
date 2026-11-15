#!/usr/bin/env python3
"""Theorem EF: the sign at alpha = 4 gamma/3.  EE.3 (iii) reduced ED's open sign of Y_n - Y_oo at alpha = 4gamma/3 (rho^3 = 3) to the
sign of S_N = -nu^3 (1-delta_N) - varrho_N and left it undecided ("no certified nu, delta_N, D'_oo at 4gamma/3").  Here all three
become integrals of objects that do not depend on alpha, and S_N < 0 is certified for EVERY N >= 1:

  EF.1 (DERIVED; every alpha in (gamma,1), N >= 1; dominated convergence as in EE.1, monotone convergence as in ED.1).
       D'_oo = lim rho^-k J_k = int_0^oo x^-a G(x) dx,   h_oo = (1-delta_N) D'_oo = int_0^oo x^-a G(x) eta_N(Phi(5^-N x)) dx,
       nu_a = C_a/A_a = C'/D'_oo = -rho Mt_oo / ((rho-1) D'_oo),   Mt_oo = int_0^oo t^-a r_oo(t) (omega_oo(q(t)) - 1) dt,
       with G, Phi, A, eta_N the level-oo profile objects of EE.1 (alpha-free) and r_oo the corner resolvent (AQ (4) node enclosures).
  EF.2 (THEOREM; exact algebra on EE.3).  S_N = (1-delta_N)(|nu|^3 - e*/D'_oo) - (3/2) kappa_N / D'_oo, so |nu|^3 < e*/D'_oo
       implies S_N < 0 for every N >= 1 (1-delta_N > 0, kappa_N > 0).
  EF.3 (arb, 128 bits, alpha = 4gamma/3 enclosed as a ball; monotone cell bounds + explicit head/tail).  |nu|^3 < e*/D'_oo, so
       S_N < 0 for every N >= 1: at alpha = 4gamma/3, rho^n (Y_n - Y_oo) -> S_N < 0 and Y_n < Y_oo for all large n.  Brackets of
       D'_oo, Mt_oo, C', nu, e*, and (N = 1, 2, 3, 6, 10) 1-delta_N, kappa_N, varrho_N, -nu^3(1-delta_N), S_N.
  EF.4 (consistency, arb at alpha = 0.8).  EF.1's formulas, certified independently, overlap AZ (3)'s A_0.8, C_0.8, nu_0.8 and
       EA.4's 1 - delta_N = C_delta/(-nu).
  EF.5 (OBSERVED, floats).  Trapezoid values of the same integrals.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

FAILS: list[str] = []
RESULT: dict = {}
NS = (1, 2, 3, 6, 10)
GAMMA = math.log(3) / math.log(5)
BI4 = "0.6031"            # BI.4: r_oo(s) <= 0.6031 s^(gamma-1), so G(x) <= 0.6031 x^(gamma-1) (EE.1)


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def _ee():
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import gasket_fractional_corner_remainder_limit as ee
    return ee


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
def rinf_float(t):
    """Floats: r_oo(t) from the affine recursion u_(k+1)(s) = beta(s) u_k(ps) - e(s), u = 1/r, started at the first orbit point
    t_K > 1e30 with u in [t_K, t_K + 3] (AQ (4)'s node enclosures T^K(1/s) >= r_oo >= T^K r_0)."""
    t = np.asarray(t, float)
    K = int(max(0.0, math.log(1 / t.min()) / math.log(5))) + 12
    ts = [t]
    with np.errstate(over="ignore"):
        for _ in range(K):
            ts.append(np.minimum(ts[-1] * (ts[-1] + 5), 1e300))
    v = ts[K] + 1.5
    with np.errstate(over="ignore", invalid="ignore"):
        for j in range(K - 1, -1, -1):
            s = ts[j]
            v = np.where(s > 1e30, s + 1.5, (s + 6) / ((s + 2) * (s + 5)) * v - 2 * s / (2 + s))
    return 1 / v


def mt_float(al: float, n_per: int = 100, lo: int = -40, hi: int = 30):
    """Floats: Mt_oo = int t^-a r_oo (omega_oo(q(t)) - 1) dt (trapezoid in log t)."""
    t = 5.0 ** (np.arange(lo * n_per, hi * n_per + 1) / n_per)
    r = rinf_float(t)
    q = 2 * t / (5 + np.sqrt(25 + 4 * t))
    b = 2 * t / (q + 6)
    L = (1 + al) * np.log1p(q / 5) + np.log1p(q / 2) - np.log1p(2 * q / 5) - np.log1p(q / 6) - np.log1p(-b * r)
    return float(np.trapezoid(t ** (1 - al) * r * np.expm1(L), np.log(t)))


def critical_floats(al: float, Ns=NS, n_per: int = 60, lo: int = -25, hi: int = 8):
    """Floats (OBSERVED): D'_oo, h_N, Mt_oo, C', nu, e*, kappa_N, varrho_N, S_N by EF.1's formulas."""
    ee = _ee()
    rho = 5**al / 3
    x = 5.0 ** (np.arange(lo * n_per, hi * n_per + 1) / n_per)
    u = np.log(x)
    G, A, _, zs = ee.profile_float(x)
    band = x >= 5.0 ** (hi - 1)
    Ib = np.trapezoid((x ** (1 - al) * G)[band], u[band])
    Dp = np.trapezoid(x ** (1 - al) * G, u) + x[0] ** (1 - al) / (1 - al) / 3 + Ib / (rho - 1)
    Mt = mt_float(al)
    Cp = -Mt * rho / (rho - 1)
    nu = Cp / Dp
    lf = ee.limits_float(al, Ns=Ns)
    out = {"Dp": Dp, "Mt": Mt, "Cp": Cp, "nu": nu, "d": lf["d"], "e": lf["e"], "om": {}, "kappa": lf["kappa"], "varrho": {}, "S": {}}
    for N in Ns:
        xx = 5.0 ** (np.arange(lo * n_per, (N + 6) * n_per + 1) / n_per)
        with np.errstate(all="ignore"):
            GG, _, _, zz = ee.profile_float(np.minimum(xx, 5.0**7))
            kk = np.maximum(0.0, np.round(np.log(xx / np.minimum(xx, 5.0**7)) / np.log(5)))
            GG = ee.profile_float(xx / 5.0**kk)[0] * 0.6**kk
            _, _, _, zz = ee.profile_float(xx)
            et = np.nan_to_num(ee.eta_float(N, zz))
        eta0 = 1 - 1.5 * 0.6**N
        h = np.trapezoid(xx ** (1 - al) * GG * et, np.log(xx)) + xx[0] ** (1 - al) / (1 - al) * eta0 / 3
        om = h / Dp
        vr = (om * lf["e"] + 1.5 * lf["kappa"][N]) / Dp
        out["om"][N], out["varrho"][N], out["S"][N] = om, vr, -nu**3 * om - vr
    return out


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    import sympy as sp
    s, t, a, x, k, q, r, b, w = sp.symbols("s t alpha x k q r b w", positive=True)
    rho = 5**a / 3
    check(sp.simplify(rho ** (-k) * 5 ** (k * a) * 5 ** (-k) - (sp.Rational(3, 5)) ** k) == 0,
          "s = 5^-k x: rho^-k int s^-a r_k(s) ds = int x^-a (3/5)^k r_k(5^-k x) dx = int x^-a R_k(x) dx (D'_k), same for H_k")
    check(sp.simplify((3 * 5 ** (-a)) ** k - rho ** (-k)) == 0,
          "x = 5^j y: (5^j)^(1-a) (3/5)^j = (3 5^-a)^j = rho^-j (scaling tail of int x^-a G)")
    a43 = 4 * sp.log(3) / (3 * sp.log(5))
    check(sp.simplify(sp.expand_log((5**a43 / 3) ** 3, force=True) - 3) == 0, "alpha = 4gamma/3: rho = 3^(1/3), rho^3 = 3")
    check(sp.expand((1 + q / 5) * (1 + q / 2) - (1 + 2 * q / 5) * (1 + q / 6) - (2 * q / 15 + q**2 / 30)) == 0,
          "omega^0(q) = (1+q/5)^(1+a)(1+q/2)/((1+2q/5)(1+q/6)) >= 1: (1+q/5)(1+q/2) - (1+2q/5)(1+q/6) = 2q/15 + q^2/30 >= 0")
    beta = (q + 6) / ((q + 2) * (q + 5))
    om0 = (1 + q / 5) ** a / ((1 + 2 * q / 5) * sp.Rational(5, 3) * beta)
    check(sp.simplify(om0 - (1 + q / 5) ** (1 + a) * (1 + q / 2) / ((1 + 2 * q / 5) * (1 + q / 6))) == 0,
          "omega_k(q) = omega^0(q)/(1 - b(q) r_k(t)), omega^0 = (1+q/5)^a/((1+2q/5)(5/3)beta(q)) in the product form above")
    F = r * (w / (1 - b * r) - 1)
    check(sp.simplify(sp.diff(F, r) - (w / (1 - b * r) ** 2 - 1)) == 0,
          "dF/dr = omega^0/(1-br)^2 - 1 >= 0 (omega^0 >= 1): F(t, r) = r(omega_oo - 1) is increasing in r (upper cell bound from r_oo(left node))")
    qq = 2 * t / (5 + sp.sqrt(25 + 4 * t))
    check(sp.simplify(sp.radsimp(qq * (qq + 5)) - t) == 0, "q(t) = 2t/(5 + sqrt(25+4t)) (cancellation-free root of q(q+5) = t), q <= t/5")
    check(sp.expand(6 * (1 + q / 6) - 5 * (1 + q / 5)) == 1 and sp.expand(5 * (1 + 2 * q / 5) - 4 * (1 + q / 2)) == 1,
          "tail: (1+q/5)/(1+q/6) <= 6/5, (1+q/2)/(1+2q/5) <= 5/4, so omega^0 <= (3/2)(1+q/5)^a; b r_oo <= 2/(q+6) <= 1/3 (r_oo <= 1/t)")
    check(sp.simplify(sp.integrate(sp.Rational(9, 4) * t ** (-1 - a / 2), (t, x, sp.oo), conds="none") - sp.Rational(9, 2) * x ** (-a / 2) / a) == 0,
          "tail t > X: t^-a F <= (9/4) t^(-1-a/2) (1 + sqrt(t)/5 <= sqrt(t) for t >= 25/16), integral (9/(2a)) X^(-a/2)")
    eta = sp.Rational(-1, 2)
    ok = True
    for N in range(1, 12):
        eta = sp.Rational(2, 5) + sp.Rational(3, 5) * eta
        ok = ok and sp.simplify(eta - (1 - sp.Rational(3, 2) * sp.Rational(3, 5) ** N)) == 0
    check(ok, "eta_N(0) = alpha_1(0) + beta(0) eta_(N-1)(0) = 2/5 + (3/5) eta_(N-1)(0) = 1 - (3/2)(3/5)^N (head of h_oo)")
    nu, om, es, Dp, kap = sp.symbols("nu om e D kappa")
    vr = (om * es + sp.Rational(3, 2) * kap) / Dp
    check(sp.simplify((-nu**3 * om - vr) - (om * (-nu**3 - es / Dp) - sp.Rational(3, 2) * kap / Dp)) == 0,
          "EF.2: S_N = -nu^3(1-delta_N) - varrho_N = (1-delta_N)(|nu|^3 - e*/D'_oo) - (3/2) kappa_N/D'_oo  (nu < 0)")
    n = sp.Symbol("n", positive=True)
    R = sp.Symbol("R")
    check(sp.simplify((rho**3 / 3) ** n * 3**n * R - rho ** (3 * n) * R) == 0,
          "EE.3 (iii): rho^n (rho^(2n) R_n) = (rho^3/3)^n 3^n R_n -> varrho_N when rho^3 = 3")


# ----------------------------------------------------------------- (B) arb
def _arb():
    from flint import arb, ctx
    return arb, ctx


def alpha_ball():
    arb, ctx = _arb()
    ctx.prec = 128
    return 4 * arb(3).log() / (3 * arb(5).log())


def hull(lo, hi):
    arb, _ = _arb()
    return arb(arb(lo).lower()).union(arb(arb(hi).upper()))


def rinf_arb(t):
    """Enclosure [lo, hi] of r_oo(t) at a point t: AQ (4)'s node enclosures T^K r_0 <= r_oo <= T^K(1/s) along the orbit
    t_j = p^j(t), i.e. u = 1/r with u_K in [t_K, t_K + 3] and u_j = beta(t_j) u_(j+1) - e(t_j) (monotone in u_(j+1))."""
    arb, _ = _arb()
    ts = [t]
    big = arb(2) ** 200
    while not bool(ts[-1] > big):
        z = ts[-1]
        ts.append(z * (z + 5))
    vl, vh = ts[-1], ts[-1] + 3
    for z in reversed(ts[:-1]):
        bz = (z + 6) / ((z + 2) * (z + 5))
        ez = 2 * z / (2 + z)
        vl, vh = bz * vl - ez, bz * vh - ez
    return 1 / vh, 1 / vl


def F_arb(tb, r, al):
    """F(t, r) = r (omega^0(q)/(1 - b r) - 1), t a ball (interval extension), q = q(t), b = 2t/(q+6)."""
    q = 2 * tb / (5 + (25 + 4 * tb).sqrt())
    b = 2 * tb / (q + 6)
    L = (1 + al) * (q / 5).log1p() + (q / 2).log1p() - (2 * q / 5).log1p() - (q / 6).log1p() - (-(b * r)).log1p()
    return r * L.expm1()


def certify_mt(al, n: int = 400, lo: int = -40, hi: int = 30):
    """Mt_oo = int t^-a F(t, r_oo(t)) dt with F >= 0 increasing in r and r_oo nonincreasing (AQ (3)): on a cell [a, b] of a log grid
    (exact float nodes) F in [inf F(t-ball, r^-(b)), sup F(t-ball, r^+(a))], times W = int_cell t^-a.  Head t < 5^lo:
    t^-a F <= 0.6031 c1 t^(2gamma-1-a) (r_oo <= 0.6031 t^(gamma-1), b <= t/3, omega^0 <= e^(0.9q) <= e^(0.18 t)); tail t > 5^hi:
    t^-a F <= (9/4) t^(-1-a/2)."""
    arb, _ = _arb()
    g = arb(3).log() / arb(5).log()
    ts = [arb(float(5.0 ** (lo + i / n))) for i in range((hi - lo) * n + 1)]
    rr = [rinf_arb(tv) for tv in ts]
    Lo, Up = arb(0), arb(0)
    ok = True
    for i in range(len(ts) - 1):
        a, b = ts[i], ts[i + 1]
        tb = a.union(b)
        W = (b ** (1 - al) - a ** (1 - al)) / (1 - al)
        fu = F_arb(tb, rr[i][1].union(rr[i][1]), al)
        fl = F_arb(tb, rr[i + 1][0], al)
        ok = ok and fu.is_finite() and fl.is_finite()
        Up += W * arb(fu.upper())
        if bool(fl.lower() > 0):
            Lo += W * arb(fl.lower())
    eps = ts[0]
    Gc = arb(BI4)
    c1 = (arb("0.18") * (arb("0.18") * eps).exp() * eps ** (1 - g) + Gc / 3) / (1 - Gc / 3 * eps**g)
    head = Gc * c1 * eps ** (2 * g - al) / (2 * g - al)
    X = ts[-1]
    tail = arb(9) / (2 * al) * X ** (-al / 2)
    return hull(Lo, Up + head + tail), {"head": head, "tail": tail, "gap": Up - Lo, "ok": ok}


def certify_profile_integrals(al, Ns=NS, n: int = 400, lo: int = -25, hiD: int = 8):
    """D'_oo = int x^-a G and h_N = int x^-a G eta_N(Phi(5^-N x)), both integrands nonincreasing (EE.1, BL.1): monotone cell
    bounds on a log grid.  D'_oo head x < 5^lo: G in [G(5^lo), 1/3]; tail x > 5^hiD by scaling from the band [5^(hiD-1), 5^hiD]:
    G(5y) = (3/5)(G + D)(y), 0 <= D <= (6/5)(5/9) A G, A nonincreasing, so the tail is in [I_band/(rho-1), I_band (1+tau)/(rho-1-tau)],
    tau = (2/3) A(5^(hiD-1)).  h_N: head [G eta_N(z_N)(5^lo), eta_N(0)/3] X0^(1-a)/(1-a); tail x > 5^(N+6): G <= 0.6031 x^(gamma-1),
    eta_N(z) <= 2/z."""
    arb, _ = _arb()
    ee = _ee()
    g = arb(3).log() / arb(5).log()
    rho = arb(5) ** al / 3
    hi = max(hiD, max(Ns) + 6)
    xs = [arb(float(5.0 ** (lo + i / n))) for i in range((hi - lo) * n + 1)]
    Gs, As, Es, PhiN = [], [], {N: [] for N in Ns}, {N: [] for N in Ns}
    for xv in xs:
        G, A, zs = ee.prof_arb(xv)
        Gs.append(G)
        As.append(A)
        for N in Ns:
            Es[N].append(G * ee.eta_arb(N, zs))
            PhiN[N].append(zs[N])
    W = [(xs[i + 1] ** (1 - al) - xs[i] ** (1 - al)) / (1 - al) for i in range(len(xs) - 1)]
    iD = hiD * n - lo * n       # index of 5^hiD
    iB = (hiD - 1) * n - lo * n
    Lo = sum((W[i] * arb(Gs[i + 1].lower()) for i in range(iD)), arb(0))
    Up = sum((W[i] * arb(Gs[i].upper()) for i in range(iD)), arb(0))
    Ibl = sum((W[i] * arb(Gs[i + 1].lower()) for i in range(iB, iD)), arb(0))
    Ibu = sum((W[i] * arb(Gs[i].upper()) for i in range(iB, iD)), arb(0))
    tau = arb(2) / 3 * arb(As[iB].upper())
    X0w = xs[0] ** (1 - al) / (1 - al)
    Dl = Lo + X0w * arb(Gs[0].lower()) + Ibl / (rho - 1)
    Du = Up + X0w / 3 + Ibu * (1 + tau) / (rho - 1 - tau)
    out = {"Dp": hull(Dl, Du), "tau": tau, "Dgap": Du - Dl, "h": {}, "ok_tail": {}}
    Gc = arb(BI4)
    for N in Ns:
        iN = (N + 6) * n - lo * n
        hl = sum((W[i] * arb(Es[N][i + 1].lower()) for i in range(iN)), arb(0)) + X0w * arb(Es[N][0].lower())
        eta0 = 1 - arb(3) / 2 * (arb(3) / 5) ** N
        Xn = xs[iN]
        P = arb(PhiN[N][iN].lower())
        tailN = 2 / P * Gc * Xn ** (g - al) / (al - g)
        hu = sum((W[i] * arb(Es[N][i].upper()) for i in range(iN)), arb(0)) + X0w * eta0 / 3 + tailN
        out["h"][N] = hull(hl, hu)
        out["ok_tail"][N] = bool(P > arb(10) ** 20)
    return out


def certify_alpha(al, label: str, n_p: int, n_m: int, n_d: int, n_k: int, Ns=NS):
    arb, _ = _arb()
    ee = _ee()
    t0 = time.time()
    rho = arb(5) ** al / 3
    kap_s = 2 * (arb.pi() * al).sin() / arb.pi()
    pi = certify_profile_integrals(al, Ns, n=n_p)
    Dp = pi["Dp"]
    t1 = time.time()
    Mt, mi = certify_mt(al, n=n_m)
    t2 = time.time()
    Cp = -Mt * rho / (rho - 1)
    nu = Cp / Dp
    d, _ = ee.certify_dstar(al, n=n_d)
    es = 3 * d / (2 * (1 - rho / 3))
    t3 = time.time()
    f = lambda v, k=7: f"[{float(v.lower()):.{k}g}, {float(v.upper()):.{k}g}]"
    print(f"   alpha = {label}: D'_oo = int x^-a G in {f(Dp)} (A = kappa D'_oo in {f(kap_s * Dp)}; tau {float(pi['tau'].upper()):.1e}) ({t1 - t0:.1f}s)")
    print(f"     Mt_oo in {f(Mt)} (cell gap {float(mi['gap'].upper()):.1e}, head {float(mi['head'].upper()):.1e}, tail {float(mi['tail'].upper()):.1e}) "
          f"({t2 - t1:.1f}s); C' = -rho Mt_oo/(rho-1) in {f(Cp)}, C = kappa C' in {f(kap_s * Cp)}, nu = C'/D'_oo in {f(nu)}")
    print(f"     d* in {f(d)}, e* in {f(es)} ({t3 - t2:.1f}s)")
    check(mi["ok"], f"alpha = {label}: every Mt cell enclosure finite (1 - b r^+ > 0 on every cell)")
    out = {"Dp": Dp, "Mt": Mt, "Cp": Cp, "nu": nu, "d": d, "e": es, "rho": rho, "kap_s": kap_s, "om": {}, "kappa": {}, "varrho": {},
           "nu3om": {}, "S": {}}
    for N in Ns:
        t4 = time.time()
        om = pi["h"][N] / Dp
        kN, ok_t = ee.certify_kappa(N, al, n=n_k)
        vr = (om * es + 3 * kN / 2) / Dp
        n3 = -nu**3 * om
        S = om * (-nu**3 - es / Dp) - 3 * kN / (2 * Dp)
        out["om"][N], out["kappa"][N], out["varrho"][N], out["nu3om"][N], out["S"][N] = om, kN, vr, n3, S
        print(f"     N = {N:2d}: 1-delta_N in {f(om)}, kappa_N in {f(kN, 6)}, varrho_N in {f(vr, 6)}, -nu^3(1-delta_N) in {f(n3, 6)}, "
              f"S_N in {f(S, 6)} ({time.time() - t4:.1f}s)")
        check(pi["ok_tail"][N] and ok_t and bool(om > 0) and bool(kN > 0), f"alpha = {label}, N = {N}: 1-delta_N > 0, kappa_N > 0 (tails Phi > 1e20)")
    return out


def main(full: bool = False) -> int:
    t0 = time.time()
    print("Theorem EF: the sign of Y_n - Y_oo at alpha = 4gamma/3 (S_N < 0 for every N >= 1). Not thrust.")
    symbolic_checks()
    arb, ctx = _arb()
    ctx.prec = 128
    ee = _ee()
    import gasket_fractional_corner_limit_mobius as eb
    n_p, n_m, n_d, n_k = (4000, 20000, 6000, 2000) if full else (400, 1000, 1000, 400)
    # ---- EF.3 at alpha = 4gamma/3 (ball)
    al = alpha_ball()
    rho = arb(5) ** al / 3
    check(bool(abs(rho**3 - 3) < arb(10) ** -30) and bool(al > arb("0.9101415926")) and bool(al < arb("0.9101415927")),
          "alpha = 4gamma/3 enclosed in a ball of radius < 1e-37 (0.91014159264798...), rho^3 = 3 to 1e-30")
    c = certify_alpha(al, "4gamma/3", n_p, n_m, n_d, n_k)
    RESULT["crit"] = c
    U = c["e"] / c["Dp"] + c["nu"] ** 3
    RESULT["U"] = U
    print(f"   EF.2/EF.3: e*/D'_oo - |nu|^3 in [{float(U.lower()):.5f}, {float(U.upper()):.5f}]  (|nu|^3 in "
          f"[{float((-c['nu'] ** 3).lower()):.5f}, {float((-c['nu'] ** 3).upper()):.5f}], e*/D'_oo in "
          f"[{float((c['e'] / c['Dp']).lower()):.5f}, {float((c['e'] / c['Dp']).upper()):.5f}])")
    check(bool(c["nu"] < 0), "EF.3: nu_(4gamma/3) < 0 certified (DZ: C_a < 0)")
    check(bool(U > 0), "EF.3: |nu|^3 < e*/D'_oo at alpha = 4gamma/3, hence S_N < 0 for EVERY N >= 1 (EF.2)")
    check(all(bool(c["S"][N] < 0) for N in NS), "EF.3: S_N < 0 certified directly for N = 1, 2, 3, 6, 10")
    check(all(bool(c["varrho"][N] > c["nu3om"][N]) for N in NS), "EF.3: varrho_N > -nu^3(1-delta_N) (disjoint brackets) for N = 1, 2, 3, 6, 10")
    # ---- EF.4 consistency at 0.8
    c8 = certify_alpha(arb("0.8"), "0.8", n_p, n_m, n_d // 2, n_k // 2, Ns=(1, 2, 3, 10))
    RESULT["c08"] = c8
    iv = lambda lo, hi: hull(arb(lo), arb(hi))
    A8, C8, nu8 = iv(*eb.A08), iv(*eb.C08), iv(*eb.NU08)
    check(bool((c8["kap_s"] * c8["Dp"]).overlaps(A8)), "EF.4: kappa int x^-a G overlaps AZ (3)'s A_0.8 = 2.39640970397922...")
    check(bool((c8["kap_s"] * c8["Cp"]).overlaps(C8)), "EF.4: -kappa rho Mt_oo/(rho-1) overlaps AZ (3)'s C_0.8 = -1.74325717839...")
    check(bool(c8["nu"].overlaps(nu8)), "EF.4: EF.1's nu overlaps AZ (3)'s nu_0.8 = -0.72744538444...")
    for N in (1, 2, 3, 10):
        check(bool(c8["om"][N].overlaps(iv(*eb.CD_EA[N]) / (-nu8))), f"EF.4: N = {N}: h_oo/D'_oo overlaps EA.4's C_delta(0.8,{N})/(-nu)")
    check(all(bool(c8["varrho"][N].overlaps(hull(*v))) for N, v in
              {1: (arb("0.0287924"), arb("0.0288038")), 2: (arb("0.239231"), arb("0.239288")), 3: (arb("0.705919"), arb("0.706085")),
               10: (arb("446.091"), arb("446.198"))}.items()),
          "EF.4: varrho_N(0.8) through EF.1's D'_oo, 1-delta_N overlaps EE.4's brackets (N = 1, 2, 3, 10)")
    # ---- EF.4b: nu_alpha at other alpha (DZ's open item "a certified nu_alpha at other alpha")
    RESULT["nu_other"] = {}
    for a_s in ("0.7", "0.9", "0.95"):
        a = arb(a_s)
        rho_a = arb(5) ** a / 3
        kap_a = 2 * (arb.pi() * a).sin() / arb.pi()
        Dp = certify_profile_integrals(a, Ns=(1,), n=n_p)["Dp"]
        Mt, _ = certify_mt(a, n=n_m)
        nu = -rho_a * Mt / ((rho_a - 1) * Dp)
        RESULT["nu_other"][a_s] = (kap_a * Dp, kap_a * (-rho_a * Mt / (rho_a - 1)), nu)
        print(f"   alpha = {a_s}: A in [{float((kap_a * Dp).lower()):.6f}, {float((kap_a * Dp).upper()):.6f}], C in "
              f"[{float((kap_a * -rho_a * Mt / (rho_a - 1)).lower()):.6f}, {float((kap_a * -rho_a * Mt / (rho_a - 1)).upper()):.6f}], "
              f"nu in [{float(nu.lower()):.6f}, {float(nu.upper()):.6f}]")
        check(bool(nu < 0), f"EF.4: nu_{a_s} < 0 certified (DZ's C_a < 0)")
    # ---- EF.5 floats
    fl = critical_floats(4 * GAMMA / 3)
    print("   floats alpha = 4gamma/3: D'_oo = %.8f, Mt_oo = %.8f, C' = %.8f, nu = %.8f, e* = %.8f; S_N: %s" % (
        fl["Dp"], fl["Mt"], fl["Cp"], fl["nu"], fl["e"], ", ".join(f"{N}: {fl['S'][N]:.6g}" for N in NS)))
    check(all(bool(c["S"][N].lower() < fl["S"][N] < c["S"][N].upper()) for N in NS) and bool(c["nu"].lower() < fl["nu"] < c["nu"].upper()),
          "EF.5: floats of nu and S_N lie inside the arb brackets")
    RESULT["floats"] = fl
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EF checks passed")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    sys.exit(main(full="--full" in sys.argv))
