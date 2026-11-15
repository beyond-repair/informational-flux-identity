#!/usr/bin/env python3
"""Theorem EE: the 3^-n coefficient has a limit.  ED's kept limit "no limit of (3/rho)^n eps_n or 3^n R_n is proved" is
resolved: for every alpha in (gamma,1),
    (3/rho)^k (Mt_(k+1) - Mt_k) -> d*,    (3/rho)^n eps_n -> e* = 3 d*/(2(1 - rho/3)),
    3^n eta_n -> (3/2) kappa_N,          3^n R_n -> varrho_N = [(1-delta_N) e* + (3/2) kappa_N] / D'_oo > 0,
with d* = (1/3) int_0^oo x^(1-alpha) D(x) (G(x) + Gh(x)) dx, kappa_N = (1/3) int_0^oo x^(1-alpha) G(x)^2 eta_N(Phi(5^-N x)) dx,
where G(x) = lim (3/5)^k r_k(5^-k x) (BK's scaled profile R_k at level oo), Gh(x) = (5/3) G(5x), Phi(x) = lim p^k(5^-k x)
(the Poincare function of p, Phi(5x) = p(Phi(x))), A(x) = prod_(m>=1) (5/3) beta(Phi(5^-m x)) and D = Gh - G = G Gh A (Phi+6)/(Phi+5).

Setting (GASKET.md, Theorems AQ, AZ, BK, BL, DZ, EB, EC, ED): alpha in (gamma,1), rho = 5^alpha/3, rho' = rho/3, p(s) = s(s+5),
q = p^-1, beta(s) = (s+6)/((s+2)(s+5)), e(s) = 2s/(2+s), r_(k+1)(s) = r_k(ps)/(beta(s) - e(s) r_k(ps)), r_0 = 1/(3+s);
Mt_k (ED.1), eps_n (AZ (1)), eta_n (DZ.1/DZ.2), R_n (EB.1), w_k (EC.5), eta_N (BL: eta_N = alpha_1 + beta eta_(N-1) o p, eta_0 = -1/2).

Theorem EE (checked here; exit 1 on any failure):
  EE.1 (THEOREM; every alpha in (gamma,1)).  The limits G, Phi, A exist pointwise (monotone/dominated limits of the exact
       affine recursion u_(k+1)(s) = beta(s) u_k(ps) - e(s) for u_k = 1/r_k), and (3/rho)^k (Mt_(k+1) - Mt_k) -> d* > 0 (dominated
       convergence in x = 5^k t; the region t >= 1 contributes o(1)).
  EE.2 (THEOREM, from EE.1 and ED.1).  (3/rho)^n eps_n -> e* = 3 d*/(2(1-rho')),  rho'^-k (Mt_oo - Mt_k) -> d*/(1-rho').
  EE.3 (THEOREM; DZ.1, DZ.2, EB.1).  3^n eta_n -> (3/2) kappa_N and 3^n R_n -> varrho_N > 0; so R_(n+1)/R_n -> 1/3 (EB.4's and
       EC.5's ratio statements become theorems) and 3^k w_k -> 2 varrho_N (1 - 3/rho)/(9 (1-delta_N)^2) < 0.
  EE.4 (arb at alpha = 0.8, outward; monotone cell bounds on a log grid + explicit head/tail).  d*, e*, kappa_N, varrho_N brackets.
  EE.5 (OBSERVED, floats).  Direct quadrature of (3/rho)^k (Mt_(k+1)-Mt_k) and of kappa_k approaches the limits; EB.1 + EE.3 with
       the leading rho^-k corrections reproduce EC.5's observed w_(k+1)/w_k to ~1e-4.

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


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
def profile_float(x, M: int = 50):
    """Floats: (G, A, Phi, orbit) at x (array, x <= 5^7) from the orbit z_m = Phi(5^-m x), z_M ~ 5^-M x (OBSERVED level)."""
    x = np.asarray(x, float)
    zs = [x * 5.0 ** (-M)]
    for _ in range(M):
        zs.append(np.minimum(zs[-1] * (zs[-1] + 5), 1e150))   # clipped: only used where eta/A are negligible
    zs = zs[::-1]                                  # zs[m] = Phi(5^-m x)
    V = 3 + zs[0]
    A = np.ones_like(x)
    for m in range(1, M + 1):
        z = zs[m]
        fb = (5 / 3) * (z + 6) / ((z + 2) * (z + 5))
        V = fb * V - 2 * z / (2 + z) * (5 / 3) ** m
        A = A * fb
    return 1 / V, A, zs[0], zs


def eta_float(N: int, zs):
    """eta_N(z_N) from the orbit zs[m] = z_m = Phi(5^-m x): BL's recursion eta_i(z) = alpha_1(z) + beta(z) eta_(i-1)(p z), eta_0 = -1/2,
    evaluated inside-out at z_1, z_2, ..., z_N (p(z_m) = z_(m-1))."""
    v = -0.5 * np.ones_like(zs[0])
    for i in range(1, N + 1):
        w = zs[i]
        v = (w + 4) / ((w + 2) * (w + 5)) + (w + 6) / ((w + 2) * (w + 5)) * v
    return v


def limits_float(al: float, Ns=(1, 2, 3), n_per: int = 60, lo: int = -25, hi: int = 7):
    """d*, e*, kappa_N by trapezoid in log x (floats).  For kappa_N the grid stops at 5^(hi): G is not needed beyond (eta_N tail)
    only when N + 6 <= hi, so N <= 1 here is exact to float precision and N = 2, 3 use the scaling G(5x) ~ (3/5) G(x) above 5^7."""
    j = np.arange(lo * n_per, hi * n_per + 1)
    x = 5.0 ** (j / n_per)
    u = np.log(x)
    G, A, Phi, _ = profile_float(x)
    G5 = profile_float(5 * x)[0]
    Gh = 5 / 3 * G5
    d = np.trapezoid(x ** (2 - al) * G * Gh * A * (Phi + 6) / (Phi + 5) * (G + Gh), u) / 3
    rho = 5**al / 3
    out = {"d": d, "e": 3 * d / (2 * (1 - rho / 3)), "kappa": {}}
    for N in Ns:
        jj = np.arange(lo * n_per, (N + 6) * n_per + 1)
        xx = 5.0 ** (jj / n_per)
        xs = np.minimum(xx, 5.0**6)
        k = np.maximum(0.0, np.round(np.log(xx / xs) / np.log(5)))
        Gs, _, _, _ = profile_float(xx / 5.0**k)
        Gx = Gs * 0.6**k                               # G(5^k y) ~ (3/5)^k G(y) for y >= 5^6 (relative error ~1e-28): floats only
        with np.errstate(all="ignore"):
            _, _, _, zs = profile_float(xx)
            et = eta_float(N, zs)
        out["kappa"][N] = np.trapezoid(xx ** (2 - al) * Gx**2 * et, np.log(xx)) / 3
    return out


def direct_increment_float(al: float, k: int, h: float = 2e-3):
    """OBSERVED: (3/rho)^k (Mt_(k+1) - Mt_k) by direct quadrature of the positive form
    (9/5)^k int x^-a (r'-r)[omega^0/((1-br)(1-br')) - 1] dx, r = r_k, r' = r_(k+1) at t = 5^-k x (no J-differences)."""
    y = np.arange(-40.0, math.log(5.0**k) + 6, h)
    x = np.exp(y)
    t = x / 5.0**k
    ss = [t]
    for _ in range(k + 1):
        ss.append(np.minimum(ss[-1] * (ss[-1] + 5), 1e150))
    beta = lambda s: (s + 6) / ((s + 2) * (s + 5))
    e = lambda s: 2 * s / (2 + s)

    def u(n):
        v = 3 + ss[n]
        for j in range(n - 1, -1, -1):
            v = beta(ss[j]) * v - e(ss[j])
        return v
    uk, uk1 = u(k), u(k + 1)
    Pi = np.ones_like(t)
    for j in range(k):
        Pi = Pi * beta(ss[j])
    r, r1 = 1 / uk, 1 / uk1
    dr = r * r1 * Pi * (ss[k] + 6) / (ss[k] + 5)       # r_(k+1) - r_k = r_k r_(k+1) prod beta(s_j) (s_k+6)/(s_k+5)  (exact)
    q = (-5 + np.sqrt(25 + 4 * t)) / 2
    q = np.where(t < 1e-8, t / 5 - t * t / 125, q)
    b = 2 * t / (q + 6)
    om0m1 = np.expm1((1 + al) * np.log1p(q / 5) + np.log1p(q / 2) - np.log1p(2 * q / 5) - np.log1p(q / 6))
    brk = om0m1 + (1 + om0m1) * (b * (r + r1) - b * b * r * r1) / ((1 - b * r) * (1 - b * r1))
    f = x ** (1 - al) * (9 / 5) ** k * dr * brk
    return float(np.trapezoid(f, y))


def ec5_model(al, nu, Dinf, om, es, kap1, ks):
    """EB.1 with eps_n ~ e* rho'^n, eta_n ~ 1.5 kappa 3^-n and the exact rho^-n terms (D'_n = D'_oo + rho^-n (C' + eps_n),
    Moebius part, 1/(1-delta^(n))): predicted w_(k+1)/w_k and 3^k w_k (OBSERVED-level comparison with EC.5)."""
    rho = 5**al / 3
    Cp = nu * Dinf

    def R(k):
        x = rho**-k
        eps = es * (rho / 3) ** k
        Dk = Dinf + x * (Cp + eps)
        return x / Dk * (om * eps / (1 + nu * x) + 1.5 * kap1 * 3.0**-k / x)

    def u(k):
        x = rho**-k
        M = (1 - om + nu * x) / (1 + nu * x)
        return 1 / (1 - M - R(k))
    w = lambda k: (u(k + 1) - u(k + 2)) - (u(k) - u(k + 1)) / rho
    return [(k, w(k + 1) / w(k), 3**k * w(k)) for k in ks]


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    import sympy as sp
    s, t, a, x, y, k = sp.symbols("s t alpha x y k", positive=True)
    beta = (s + 6) / ((s + 2) * (s + 5))
    e = 2 * s / (2 + s)
    p = s * (s + 5)
    u1 = beta * (3 + p) - e
    check(sp.simplify((3 + s) - u1 - (s + 6) / (s + 5)) == 0,
          "u_0 - u_1 = (s+6)/(s+5) in (1, 6/5]: r_(k+1)(t) - r_k(t) = r_k r_(k+1) prod_(j<k) beta(p^j t) (p^k t + 6)/(p^k t + 5)")
    check(sp.simplify(e / beta - 2 * p / (s + 6)) == 0, "b(q) := e(q)/beta(q) = 2t/(q+6) <= t/3 (t = p(q))")
    om0 = (1 + s / 5) ** a / ((1 + 2 * s / 5) * sp.Rational(5, 3) * beta)
    jac = sp.simplify(p ** (-a) * om0 * beta * sp.diff(p, s) - 3 * 5 ** (-a) * s ** (-a))
    check(sp.simplify(jac.subs(a, sp.Rational(4, 5))) == 0 and sp.simplify(jac.subs(a, sp.Rational(7, 10))) == 0,
          "t^-alpha omega^0(q) beta(q) dt = q^-alpha dq / rho (so int_(t>=1) t^-a omega^0 beta delta_(k+1)(q) dt = rho^-1 int_(q>=q(1)) q^-a delta_(k+1))")
    r, r1, b, w = sp.symbols("r r1 b w", positive=True)
    F = lambda rr: rr * (w / (1 - b * rr) - 1)
    check(sp.simplify(F(r1) - F(r) - (r1 - r) * (w / ((1 - b * r) * (1 - b * r1)) - 1)) == 0,
          "F(r') - F(r) = (r'-r)[omega^0/((1-br)(1-br')) - 1] >= 0 (F = r(omega_k(q) - 1), omega_k = omega^0/(1 - b r_k))")
    fb = sp.Rational(5, 3) * beta
    check(sp.simplify(fb - (1 + s / 6) / ((1 + s / 2) * (1 + s / 5))) == 0,
          "(5/3) beta(s) = (1+s/6)/((1+s/2)(1+s/5)): -log((5/3) beta) <= s/2 + s/5 = 0.7 s, and (5/3) beta is decreasing, <= 1")
    check(sp.expand(35 * (y + 2) * (y + 5) - 30 * y * (y + 6) - (5 * y**2 + 65 * y + 350)) == 0,
          "(5/3) beta(y) <= (35/18)/y for every y > 0 (35(y+2)(y+5) - 30 y(y+6) = 5y^2 + 65y + 350 > 0)")
    rho = 5**a / 3
    check(sp.simplify((3 / rho) ** k * 5 ** (-k * (1 - a)) - (sp.Rational(9, 5)) ** k) == 0
          and sp.simplify((sp.Rational(9, 5)) ** k * (sp.Rational(5, 3)) ** k - 3**k) == 0,
          "scaling t = 5^-k x: (3/rho)^k 5^(-k(1-alpha)) = (9/5)^k, (9/5)^k (5/3)^k = 3^k, 3^k (5/3)^k = 5^k")
    qq = (-5 + sp.sqrt(25 + 4 * x / 5**k)) / 2
    lim = sp.limit((5**k * 2 * (x / 5**k) / (qq + 6)).subs(k, sp.Symbol("K", positive=True)), sp.Symbol("K", positive=True), sp.oo)
    check(sp.simplify(lim - x / 3) == 0, "5^k b(q(5^-k x)) -> x/3 and 5^k b <= x (so the pointwise limit of the scaled integrand is x^(1-a)(Gh^2 - G^2)/3)")
    rp = sp.Symbol("rp", positive=True)
    i = sp.Symbol("i", integer=True, nonnegative=True)
    # eps_n = rho^n sum_(k>=n) rho^-k sum_(j>=k) rp^j c_j with c_j -> d*: (3/rho)^n eps_n -> (3/2) d*/(1-rp)
    check(sp.simplify(sp.Rational(3, 2) * sp.Symbol("d") / (1 - rp) - sp.Symbol("d") * sp.summation(sp.Rational(1, 3) ** i, (i, 0, sp.oo)) / (1 - rp)) == 0,
          "sum_(k>=n) 3^-(k-n) = 3/2: (3/rho)^n eps_n -> (3/2) d*/(1-rho'), 3^n eta_n -> (3/2) kappa_N")
    ww = sp.Rational(1, 3) - sp.Rational(1, 9) - (1 - sp.Rational(1, 3)) / sp.Symbol("rho")
    check(sp.simplify(ww - sp.Rational(2, 9) * (1 - 3 / sp.Symbol("rho"))) == 0,
          "w_k on Q_k ~ varrho 3^-k/(1-delta)^2: 3^k w_k -> (varrho/(1-delta)^2)(1/3 - 1/9 - (2/3)/rho) = 2 varrho (1 - 3/rho)/(9(1-delta)^2)")
    g = sp.log(3) / sp.log(5)
    check(bool(sp.N(2 * g - 1) < sp.N(g)) and bool(sp.N(2 * g) > 1),
          "exponents: 2 gamma - alpha > 0 and gamma - alpha - 1 < -1 on (gamma,1) (head/tail integrability of the dominators)")


# ----------------------------------------------------------------- (B) arb: profile enclosures and the limit integrals
def _arb():
    from flint import arb, ctx
    return arb, ctx


def prof_arb(x, extra: int = 40):
    """Rigorous enclosures of (G(x), A(x), orbit z_m = Phi(5^-m x), m = 0..M) at a point x > 0 (arb ball).
    Phi(y) in [y, y exp(y e^0.1/20)] for y <= 1.8 (ED.2 (a)); z_(m-1) = p(z_m) exactly; truncation at m = M with
    A-tail factor in [exp(-0.7 sum_(m>M) z_m), 1] and S-tail in [0, sum_(m>M) z_m (5/3)^m]."""
    arb, _ = _arb()
    xu = x.upper()
    M = int(math.ceil(math.log(max(float(xu), 1e-300)) / math.log(5) + extra))
    five = arb(5)
    c = (arb(1) / 10).exp() / 20
    yv = x / five**M
    z = arb(yv.lower()).union(yv.upper() * (yv.upper() * c).exp())
    zs = [z]
    for _ in range(M):
        z = z * (z + 5)
        zs.append(z)
    zs.reverse()
    V = 3 + zs[0]
    A = arb(1)
    f53 = five / 3
    pw = arb(1)
    for m in range(1, M + 1):
        zm = zs[m]
        fb = f53 * (zm + 6) / ((zm + 2) * (zm + 5))
        pw = pw * f53
        V = fb * V - 2 * zm / (2 + zm) * pw
        A = A * fb
    yM = arb(xu) / five ** (M + 1)
    sfac = (yM * c).exp()
    tau = arb("0.7") * sfac * arb(xu) / five**M / 4
    T = sfac * arb(xu) / arb(3) ** M / 2
    P = arb((-tau).exp().lower()).union(arb(1))
    invG = P * V - arb(0).union(T)
    return 1 / invG, A * P, zs


def eta_arb(N, zs):
    arb, _ = _arb()
    v = arb(-1) / 2
    for i in range(1, N + 1):
        w = zs[i]
        v = (w + 4) / ((w + 2) * (w + 5)) + (w + 6) / ((w + 2) * (w + 5)) * v
    return v


def certify_dstar(al_s: str = "0.8", n: int = 2000, lo: int = -25, hi: int = 8):
    """d* = (1/3) int x^(1-a) Dec(x) dx, Dec = G Gh A (Phi+6)/(Phi+5)(G+Gh) >= 0 nonincreasing; cells [x_i, x_(i+1)] of a log grid
    with exact (float) nodes: int_cell x^(1-a) Dec in [Dec(b) W, Dec(a) W], W = int_cell x^(1-a) (exact).  Head (x < 5^lo):
    Dec <= (1/3)(5/9)(6/5)(1/3+5/9).  Tail (x > 5^hi): Dec <= (6/5)(1/3)(5/9)(8/9) prod_(m<=hi) (35/18) 5^m / x."""
    arb, ctx = _arb()
    ctx.prec = 128
    al = arb(al_s)
    five = arb(5)
    xs = [arb(float(5.0 ** (lo + i / n))) for i in range((hi - lo) * n + 1)]
    dec = []
    for xv in xs:
        G, A, zs = prof_arb(xv)
        Gh = five / 3 * prof_arb(5 * xv)[0]
        Phi = zs[0]
        dec.append(G * Gh * A * (Phi + 6) / (Phi + 5) * (G + Gh))
    Lo, Up = arb(0), arb(0)
    for i in range(len(xs) - 1):
        a, b = xs[i], xs[i + 1]
        W = (b ** (2 - al) - a ** (2 - al)) / (2 - al)
        Lo += W * arb(dec[i + 1].lower())
        Up += W * arb(dec[i].upper())
    c0 = arb(1) / 3 * arb(5) / 9 * arb(6) / 5 * (arb(1) / 3 + arb(5) / 9)
    head = xs[0] ** (2 - al) / (2 - al) * c0
    J = hi
    X = xs[-1]
    cT = arb(6) / 5 * arb(1) / 3 * arb(5) / 9 * arb(8) / 9 * (arb(35) / 18) ** J * five ** (J * (J + 1) // 2)
    tail = cT * X ** (2 - al - J) / (J - 2 + al)
    lo_v = arb(Lo.lower()) / 3
    hi_v = (Up + head + tail) / 3
    return arb(lo_v.lower()).union(arb(hi_v.upper())), {"head": head / 3, "tail": tail / 3, "gap": (Up - Lo) / 3}


def certify_kappa(N: int, al_s: str = "0.8", n: int = 800, lo: int = -25, hi_off: int = 6):
    """kappa_N = (1/3) int x^(1-a) G^2 eta_N(Phi(5^-N x)) dx; G^2 and eta_N o Phi(5^-N .) nonincreasing (AQ, BL.1).  Head: G <= 1/3,
    eta_N <= 1; tail x > X = 5^(N+hi_off): G <= 1/3, eta_N(z) <= 2/z, Phi(5^j y) >= Phi(y)^(2^j)."""
    arb, ctx = _arb()
    ctx.prec = 128
    al = arb(al_s)
    hi = N + hi_off
    xs = [arb(float(5.0 ** (lo + i / n))) for i in range((hi - lo) * n + 1)]
    dec = []
    for xv in xs:
        G, _, zs = prof_arb(xv)
        dec.append(G * G * eta_arb(N, zs))          # eta_N(z_N) = alpha_1(z_N) + beta(z_N) eta_(N-1)(z_(N-1)), z_(m-1) = p(z_m)
    Lo, Up = arb(0), arb(0)
    for i in range(len(xs) - 1):
        a, b = xs[i], xs[i + 1]
        W = (b ** (2 - al) - a ** (2 - al)) / (2 - al)
        Lo += W * arb(dec[i + 1].lower())
        Up += W * arb(dec[i].upper())
    head = xs[0] ** (2 - al) / (2 - al) / 9
    X = xs[-1]
    PhiX = prof_arb(X)[2][N]
    ok_tail = bool(PhiX.lower() > arb(10) ** 20)
    tail = 2 * (5 * X) ** (2 - al) / (2 - al) / 9 * 2 / arb(PhiX.lower())
    lo_v = arb(Lo.lower()) / 3
    hi_v = (Up + head + tail) / 3
    return arb(lo_v.lower()).union(arb(hi_v.upper())), ok_tail


def certify_08(n_d: int = 2000, n_k: int = 800):
    import contextlib
    import io
    arb, ctx = _arb()
    ctx.prec = 128
    import gasket_fractional_corner_limit_mobius as eb

    def iv(lo, hi):
        lo, hi = arb(lo), arb(hi)
        return arb(lo.lower()).union(arb(hi.upper()))
    t0 = time.time()
    d, info = certify_dstar("0.8", n=n_d)
    al = arb("0.8")
    rho = arb(5) ** al / 3
    es = 3 * d / (2 * (1 - rho / 3))
    print(f"   alpha = 0.8: d* in [{float(d.lower()):.7f}, {float(d.upper()):.7f}] (cell gap {float(info['gap'].upper()):.1e}, head {float(info['head'].upper()):.1e}, "
          f"tail {float(info['tail'].upper()):.1e}); e* = 3d*/(2(1-rho/3)) in [{float(es.lower()):.7f}, {float(es.upper()):.7f}]  ({time.time() - t0:.1f}s)")
    check(bool(d.lower() > arb("0.29500")) and bool(d.upper() < arb("0.29516")) and bool(es.lower() > arb("0.74077")) and bool(es.upper() < arb("0.74118")),
          "EE.4: d*(0.8) in [0.29500, 0.29516], e*(0.8) in [0.74077, 0.74118]")
    nu, A = iv(*eb.NU08), iv(*eb.A08)
    kap_s = 2 * (arb.pi() * al).sin() / arb.pi()
    Dinf = A / kap_s
    with contextlib.redirect_stdout(io.StringIO()):
        import gasket_fractional_corner_remainder_lower as ed
        ed.FAILS.clear()
        rows = ed.certify_constants()
    L8, B8 = rows["0.8"]["L"], rows["0.8"]["B"]
    check(bool(es.lower() > L8.upper()) and bool(es.upper() < B8.lower()), "EE.4 inside ED.2's bracket: L_0.8 < e* < B_0.8 (ED.2 re-derived)")
    out = {"d": d, "e": es, "kappa": {}, "varrho": {}, "w": {}}
    for N in NS:
        t1 = time.time()
        kN, ok_tail = certify_kappa(N, "0.8", n=n_k)
        Cd = iv(*eb.CD_EA[N])
        om = Cd / (-nu)                                      # 1 - delta_N = C_delta/(-nu) (EA.4, DZ.2)
        vr = (om * es + 3 * kN / 2) / Dinf
        wl = 2 * vr * (1 - 3 / rho) / (9 * om * om)
        out["kappa"][N], out["varrho"][N], out["w"][N] = kN, vr, wl
        print(f"   N = {N:2d}: kappa_N in [{float(kN.lower()):.6g}, {float(kN.upper()):.6g}], lim 3^n R_n in [{float(vr.lower()):.6g}, {float(vr.upper()):.6g}], "
              f"lim 3^k w_k in [{float(wl.lower()):.5g}, {float(wl.upper()):.5g}]  ({time.time() - t1:.1f}s)")
        check(ok_tail and bool(kN.lower() > 0) and bool(vr.lower() > 0), f"EE.4 N = {N}: kappa_N > 0 and varrho_N > 0 certified (tail Phi(5^-N X) > 1e20)")
    ed4 = {1: (3.716e-5, 13.97), 2: (1.887e-4, 68.44), 3: (3.066e-4, 113.2), 6: (5.480e-4, 280.6), 10: (7.155e-4, 3.624e3)}   # ED.4 (GASKET.md)
    check(all(ed4[N][0] < float(out["varrho"][N].lower()) and float(out["varrho"][N].upper()) < ed4[N][1] for N in NS),
          "EE.4 consistent with ED.4: lim 3^n R_n lies inside ED.4's every-n >= 60 bracket for N = 1, 2, 3, 6, 10")
    lw = out["w"][1]
    check(bool(lw.upper() < 0), "EE.3 at 0.8, N = 1: lim 3^k w_k < 0 (EC.5's sign), so w_(k+1)/w_k -> 1/3 is a theorem")
    return out


def certify_other_alphas(n: int = 400):
    arb, _ = _arb()
    out = {}
    for a in ("0.7", "0.9", "0.95"):
        d, _ = certify_dstar(a, n=n)
        rho = arb(5) ** arb(a) / 3
        es = 3 * d / (2 * (1 - rho / 3))
        out[a] = (d, es)
        print(f"   alpha = {a}: d* in [{float(d.lower()):.5f}, {float(d.upper()):.5f}], e* in [{float(es.lower()):.5f}, {float(es.upper()):.5f}]")
        check(bool(d.lower() > 0), f"EE.4 alpha = {a}: d* > 0 certified (coarse grid)")
    return out


# ----------------------------------------------------------------- main
def main(full: bool = False) -> int:
    t0 = time.time()
    print("Theorem EE: (3/rho)^n eps_n -> e*, 3^n R_n -> varrho_N > 0 (limits exist; arb brackets at 0.8). Not thrust.")
    symbolic_checks()
    n_d, n_k = (40000, 8000) if full else (6000, 2000)
    c8 = certify_08(n_d, n_k)
    RESULT["c8"] = c8
    RESULT["other"] = certify_other_alphas(400 if not full else 2000)
    # ---- (C) floats: OBSERVED
    lf = limits_float(0.8)
    obs = [(k, direct_increment_float(0.8, k)) for k in (8, 12, 16)]
    print("   floats alpha = 0.8: limit formula d* = %.10f, e* = %.10f, kappa_1..3 = %s; direct (3/rho)^k (Mt_(k+1)-Mt_k): %s" % (
        lf["d"], lf["e"], ", ".join(f"{lf['kappa'][N]:.8f}" for N in (1, 2, 3)), ", ".join(f"k={k}: {v:.8f}" for k, v in obs)))
    d = c8["d"]
    check(bool(d.lower() < lf["d"]) and bool(lf["d"] < d.upper()), "EE.5: the float limit formula lies inside the arb bracket of d*")
    check(abs(obs[-1][1] - lf["d"]) < abs(obs[0][1] - lf["d"]) and abs(obs[-1][1] - lf["d"]) < 1e-4,
          "EE.5 (OBSERVED): direct quadrature of (3/rho)^k (Mt_(k+1)-Mt_k) approaches the limit formula (k = 8, 12, 16)")
    for N in (1, 2, 3):
        kN = c8["kappa"][N]
        check(bool(kN.lower() < lf["kappa"][N]) and bool(lf["kappa"][N] < kN.upper()), f"EE.5: float kappa_{N} inside its arb bracket")
    al = 0.8
    nu = -0.7274453844421553
    Dinf = 2.3964097039792205 / (2 * math.sin(math.pi * al) / math.pi)
    om = 0.03131 / -nu
    mod = ec5_model(al, nu, Dinf, om, lf["e"], lf["kappa"][1], (14, 16, 18, 20))
    ec5 = (0.33637, 0.33541, 0.33475, 0.33426)
    print("   EC.5 model (EB.1 + EE.3, leading rho^-k terms), N = 1: " + ", ".join(f"k={k}: ratio {r:.5f}, 3^k w_k {v:.4f}" for k, r, v in mod)
          + "; EC.5 observed ratios " + ", ".join(f"{v:.5f}" for v in ec5) + ", 3^k w_k ~ -5.04")
    check(max(abs(m[1] - o) for m, o in zip(mod, ec5)) < 3e-4 and abs(mod[-1][2] + 5.04) < 0.05,
          "EE.5 (OBSERVED): EB.1 + EE.3's limits reproduce EC.5's observed w_(k+1)/w_k to < 3e-4 and 3^k w_k ~ -5.04")
    RESULT["floats"] = {"lf": lf, "direct": obs, "ec5": mod}
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EE checks passed")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    sys.exit(main(full="--full" in sys.argv))
