#!/usr/bin/env python3
"""Theorem BQ: the rate of the corner deficit below log3/log5.  For alpha < gamma, 1 - delta_N^(n) is of exact order rho^n,
rho = 5^alpha/3 < 1, at EVERY level n (two-sided, not only asymptotically), and rho^-n (1 - delta_N^(n)) -> C_N(alpha) > 0;
certified brackets for C_1 at alpha = 0.01, 0.1, 0.3, 0.45.  Kept Failure BQ.1: rho^-n (1 - delta_1^(n)) is neither monotone
in n nor on one side of its limit (certified witnesses at alpha = 0.45 and 0.3).

Setting (GASKET.md, Theorems AQ, BG, BJ, BK, BL, BP): r_n(s) corner resolvent, p(y) = y(y+5), rho = 5^a/3, gamma = log3/log5,
D_n = int s^-a r_n ds, D'_n = rho^-n D_n, R_n(tau) = (3/5)^n r_n(5^-n tau), psi_k(t) = p^k(5^-k t), eta_N the neighbour profile
(eta_1 = 1/(2(z+5)), eta_N = alpha_1 + beta * eta_(N-1)(p z)), omega_n BK.1's one-level weight.  For N >= 1, n >= N:

    1 - delta_N^(n) = D_n^-1 int s^-a r_n(s) eta_N(p^(n-N) s) ds = rho^n calN^(n) / D_n,
    calN^(n) = calN_N^(n) = int tau^-a R_n(tau) eta_N(psi_(n-N)(5^-N tau)) dtau.

Theorem BQ (checked here; exit 1 on any failure):
  (1) [monotone numerator; every a in (0,1)] BK.1's change of variables tau = 5^(n+1) y, tau' = 5^n p(y) carries the integrand
      of calN^(n+1) to omega_n times that of calN^(n) (F_(n+1)(tau) = F_n(tau')), so
      calN^(n+1) - calN^(n) = rho^-n int (omega_n - 1) eta_N(p^(n+1-N) y) p^-a r_n(p) p' dy > 0
      (omega_n >= (1+y/5)^a > 1 for y > 0, BK.2; eta_N > 0, BL.1).  calN^(n) increases to P_N = int tau^-a R eta_N(psi(5^-N tau))
      (BG.3's dominated convergence), finite for every a.
  (2) [exact order rho^n; a < gamma] V^(n) := rho^-n (1 - delta_N^(n)) = calN^(n)/D_n with calN^(n) increasing to P_N and D_n
      increasing to D_inf < inf (AQ (3)); so for n > m >= N:  calN^(m)/D_inf < V^(n) < P_N/D_m, and V^(n) -> C_N = P_N/D_inf
      in (0, inf) (the limit is BG.3 specialised to the neighbour; the all-n bracket is new); hence
      (1 - delta^(n+1))/(1 - delta^(n)) -> rho.  At a = gamma, n(1 - delta^(n)) -> (2 sin(pi gamma)/pi) P_N / c_R(gamma)
      (derived from AY as stated).
  (3) [certified, N = 1] python-flint arb, outward: enclosures of calN^(n) (n = 1..6, 30), D_n (n = 1..6), D_30 (floor),
      D_inf, and E_30 >= P_1 - calN^(30) (BK.4's omega_hi, sum over k >= 30 of rho^-k 2/(2 + 5^k y)); C_1(a) and the
      all-n bracket for n >= 30 at a = 0.01, 0.1, 0.3, 0.45.
  Kept Failure BQ.1: at a = 0.45, V^(1) < V^(2) > V^(3) and V^(2) > C_1 (above the limit); at a = 0.3, V^(3) < C_1 < V^(4)
      (crosses the limit); at a = 0.01, V^(1) < C_1.  So V^(n) is not monotone in n and not one-sided.
Assumptions used as stated: AQ (3) (r_n increasing to r_inf), BI.4 (G = 0.6031), BJ.1 (1/r_n - 1/r_inf =
prod beta (1+2 eta(p^n s))), BJ.3 (g = s^(1-gamma) r_inf in [0.600945014, 0.600963487] on (0, p(1e-12)]), BJ.4/BH.2 (eta_N <= eta
<= 2/(2+z)), BK.1-BK.4 (transfer, omega >= 1, omega_hi and its head lemma), BL.1 (eta_N positive, strictly decreasing),
BG.3, AY.  These certificates are not re-run here.
Floats (consistency only): V^(n) by recursion quadrature.  Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import sympy as sp
from flint import arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))

FAILS: list[str] = []
G_BI4 = "0.6031"
G_LO, G_HI = "0.600945014", "0.600963487"
ALPHAS = ("0.01", "0.1", "0.3", "0.45")
M_TAIL = 30
LEVELS = (1, 2, 3, 4, 5, 6)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    z, y, a, n, tau = sp.symbols("z y a n tau", positive=True)
    al1 = (z + 4) / ((z + 2) * (z + 5))
    beta = (z + 6) / ((z + 2) * (z + 5))
    eta1 = sp.simplify(al1 - beta / 2)
    check(sp.simplify(eta1 - 1 / (2 * (z + 5))) == 0, "eta_1 = alpha_1 - beta/2 = 1/(2(z+5)) > 0")
    check(sp.factor(sp.simplify(al1 - 1 / (2 * (z + 5)))) == sp.factor((z + 6) / (2 * (z + 2) * (z + 5))),
          "alpha_1 - 1/(2(z+5)) = (z+6)/(2(z+2)(z+5)) > 0, so eta_N = alpha_1 + beta eta_(N-1)(pz) >= 1/(2(z+5)) for every N >= 1")
    rho = 5**a / 3
    s = 5**(-n) * tau
    jac = sp.simplify(rho**(-n) * s**(-a) * sp.diff(s, tau) / ((sp.Rational(3, 5))**n * tau**(-a)))
    check(sp.simplify(sp.powsimp(sp.expand_power_base(jac, force=True), force=True)) == 1,
          "rho^-n s^-a ds = (3/5)^n tau^-a dtau at s = 5^-n tau, so rho^-n int s^-a r_n eta_N(p^(n-N) s) ds = calN^(n)")
    t1, t2 = 5**(n + 1) * y, 5**n * y * (y + 5)
    check(sp.simplify(t2 / t1 - (1 + y / 5)) == 0 and sp.simplify(sp.diff(t2, y) / sp.diff(t1, y) - (1 + 2 * y / 5)) == 0,
          "BK.1's change of variables: tau'/tau = 1 + y/5, dtau'/dtau = 1 + 2y/5; F_(n+1)(tau) = eta_N(p^(n+1-N) y) = F_n(tau')")
    # geometric pieces of E_M
    r_ = sp.symbols("r", positive=True)
    okg = all(sp.simplify(sum(r_**(-k) for k in range(M, k1)) - (r_**(-k1) - r_**(-M)) / (1 / r_ - 1)) == 0
              for M in (1, 2, 30) for k1 in (M + 1, M + 4, M + 9))
    check(okg,
          "sum_(k=M)^(k1-1) rho^-k = (rho^-k1 - rho^-M)/(rho^-1 - 1); the k >= k1 part is 2 5^(N-1) y^-1 sum (5 rho)^-k (2/(2+x) <= 2/x)")
    check(True, "min(1, x) <= x^theta for theta in [0,1]: with theta = (1+gamma-a)/2 > gamma - a, q = rho^-1 5^-theta < 1 and "
                "S(y) <= (2 5^(N-1)/y)^theta q^M/(1-q) (head of E_M; q < 1 and 2gamma - a - theta > 0 checked in arb below)")
    # 1/r_n lower route for r_inf: BJ.1 with eta <= 2/(2+z): 1 + 2 eta <= 1 + 4/(2+z)
    check(sp.simplify((1 + 2 * (2 / (2 + z))) - (1 + 4 / (2 + z))) == 0,
          "BJ.1 with eta <= 2/(2+z): 1/r_inf >= 1/r_K - prod_(i<K) beta(p^i s) (1 + 4/(2+p^K s)) for every K")


# ----------------------------------------------------------------- (B) exact rationals
def _p(x):
    return x * (x + 5)


def r_rec(k, s):
    if k == 0:
        return 1 / (3 + s)
    x = r_rec(k - 1, _p(s))
    return (2 + s) * (5 + s) * x / ((6 + s) - 2 * s * (5 + s) * x)


def eta_rec(N, z):
    if N == 0:
        return Fr(-1, 2)
    D = (z + 2) * (z + 5)
    return (z + 4) / D + (z + 6) / D * eta_rec(N - 1, _p(z))


def exact_checks():
    zs = (Fr(1, 1000), Fr(1, 10), Fr(1), Fr(7, 2), Fr(30))
    ok = True
    for N in range(1, 6):
        vals = [eta_rec(N, zz) for zz in zs]
        ok = ok and all(v >= 1 / (2 * (zz + 5)) and v <= 2 / (2 + zz) for v, zz in zip(vals, zs))
        ok = ok and all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))
    check(ok, "exact rationals N = 1..5, z = 1/1000, 1/10, 1, 7/2, 30: 1/(2(z+5)) <= eta_N(z) <= 2/(2+z), strictly decreasing")
    ok2 = True
    for nn in range(0, 4):
        for yv in (Fr(1, 100), Fr(1, 3), Fr(2), Fr(13)):
            lhs = Fr(3, 5) * r_rec(nn + 1, yv) / r_rec(nn, _p(yv))
            B = (yv + 6) / ((yv + 2) * (yv + 5))
            e = 2 * yv / (2 + yv)
            rhs = 1 / (Fr(5, 3) * (B - e * r_rec(nn, _p(yv))))
            ok2 = ok2 and lhs == rhs and lhs * 1 >= (1 + 2 * yv / 5)   # omega_n (1+2y/5)/(1+y/5)^a = lhs, so omega_n >= (1+y/5)^a
    check(ok2, "exact rationals n = 0..3, y = 1/100, 1/3, 2, 13: (3/5) r_(n+1)(y)/r_n(p) = 1/((5/3)(beta - e r_n(p))) >= 1 + 2y/5 "
               "(BK.1-BK.2: omega_n >= (1+y/5)^a > 1, the increment sign of BQ.1)")
    ok3 = all(r_rec(k + 1, sv) > r_rec(k, sv) for k in range(4) for sv in zs)
    check(ok3, "exact rationals: r_(k+1)(s) > r_k(s), k = 0..3 (consistency with AQ (3): D_n increases to D_inf)")


# ----------------------------------------------------------------- (C) arb certificates


def _gam():
    return arb(3).log() / arb(5).log()


def _orbit(s, n, big=None):
    o = [s]
    while len(o) <= n and (big is None or o[-1] < big):
        o.append(o[-1] * (o[-1] + 5))
    return o


def _r(o, n):
    r = 1 / (3 + o[n])
    for k in range(n - 1, -1, -1):
        z = o[k]
        r = (2 + z) * (5 + z) * r / ((6 + z) - 2 * z * (5 + z) * r)
    return r


def _eta(z, N):
    zs = _orbit(z, N)
    e = arb(-1) / 2
    for k in range(N - 1, -1, -1):
        zz = zs[k]
        D = (zz + 2) * (zz + 5)
        e = (zz + 4) / D + (zz + 6) / D * e
    return e


def _grid(lo, hi, q):
    pts = [arb(lo)]
    q, hi = arb(q), arb(hi)
    while pts[-1] < hi:
        pts.append(pts[-1] * q)
    return pts


def _min_up(u, v):
    uu, vv = arb(u.upper()), arb(v.upper())
    return uu if uu < vv else vv


def numer(al, N, m, G, q="1.0005", tlo="1e-16", thi="1e24"):
    """[lo, hi] for calN_N^(m) = int tau^-a R_m eta_N(psi_(m-N)(5^-N tau)) dtau; integrand decreasing (r_m Stieltjes, BL.1)."""
    a, gam = arb(al), _gam()
    pts = _grid(tlo, thi, q)
    c, f5 = (arb(3) / 5) ** m, arb(5) ** (-m)
    vals = []
    for t in pts:
        o = _orbit(t * f5, m)
        vals.append(t ** (-a) * c * _r(o, m) * _eta(o[m - N], N))
    lo = hi = arb(0)
    for i in range(len(pts) - 1):
        w = pts[i + 1] - pts[i]
        lo += vals[i + 1] * w
        hi += vals[i] * w
    head = pts[0] ** (1 - a) / (3 * (1 - a))                                   # R_m <= 1/3, eta_N <= 1
    tail = 2 * G * arb(5) ** N * pts[-1] ** (gam - 1 - a) / (1 + a - gam)     # R_m <= G tau^(g-1), eta_N(z) <= 2/z, z >= 5^-N tau
    return arb(lo.lower()), arb((hi + head + tail).upper())


def Dlevel(al, m, q="1.0005", slo="1e-32", shi="1e16"):
    """[lo, hi] for D_m = int s^-a r_m ds (r_m decreasing; head r_m(eps) <= r_m <= r_m(0) = (5/3)^m/3; tail 1/(3+s) <= r_m <= 1/s)."""
    a = arb(al)
    pts = _grid(slo, shi, q)
    vals = [t ** (-a) * _r(_orbit(t, m), m) for t in pts]
    lo = hi = arb(0)
    for i in range(len(pts) - 1):
        w = pts[i + 1] - pts[i]
        lo += vals[i + 1] * w
        hi += vals[i] * w
    e0, T = pts[0], pts[-1]
    hlo = _r(_orbit(e0, m), m) * e0 ** (1 - a) / (1 - a)
    hhi = (arb(5) / 3) ** m / 3 * e0 ** (1 - a) / (1 - a)
    return arb((lo + hlo + T ** (-a) / (a * (1 + 3 / T))).lower()), arb((hi + hhi + T ** (-a) / a).upper())


def Dinf(al, K, G, glo, ghi, q="1.0005", slo="1e-32", shi="1e16"):
    """[lo, hi] for D_inf.  Lower: r_inf >= r_k (AQ (3), k = min(K, orbit cap)), D_k(lower) is also a floor for D_K.
    Upper: r_inf <= min(G s^(g-1) [g_hi below p(1e-12), BJ.3], 1/s, 1/L_k) with BJ.1's L_k.  Returns (lo, hi, DK_lo)."""
    a, gam = arb(al), _gam()
    pts = _grid(slo, shi, q)
    cut = arb("1e-12")
    cut = cut * (cut + 5)
    big = arb("1e100")
    up, low = [], []
    for t in pts:
        o = _orbit(t, K, big)
        k = len(o) - 1
        rk = _r(o, k)
        prod = arb(1)
        for i in range(k):
            zz = o[i]
            prod *= (zz + 6) / ((zz + 2) * (zz + 5))
        L = 1 / rk - prod * (1 + 4 / (2 + o[k]))
        u = _min_up((ghi if t <= cut else G) * t ** (gam - 1), 1 / t)
        if L.lower() > 0:
            u = _min_up(u, 1 / L)
        up.append(t ** (-a) * u)
        low.append(t ** (-a) * rk)
    lo = hi = arb(0)
    for i in range(len(pts) - 1):
        w = pts[i + 1] - pts[i]
        lo += low[i + 1] * w
        hi += up[i] * w
    e0, T = pts[0], pts[-1]
    tlo = T ** (-a) / (a * (1 + 3 / T))
    lo_inf = lo + glo * e0 ** (gam - a) / (gam - a) + tlo
    hi_inf = hi + ghi * e0 ** (gam - a) / (gam - a) + T ** (-a) / a
    lo_K = lo + tlo        # r_K >= r_k pointwise (k <= K), head dropped
    return arb(lo_inf.lower()), arb(hi_inf.upper()), arb(lo_K.lower())


def Sbound(y, al, N, M):
    a = arb(al)
    rho = arb(5) ** a / 3
    k1 = M
    while (arb(5) ** (k1 + 1 - N) * y).upper() < 1:
        k1 += 1
    s1 = (rho ** (-k1) - rho ** (-M)) / (1 / rho - 1) if k1 > M else arb(0)
    x = 1 / (5 * rho)
    return s1 + 2 * arb(5) ** (N - 1) / y * x ** k1 / (1 - x)


def Ebound(al, N, M, G, ylo="1e-60", yhi="1e6", ratio="1.005"):
    """Upper bound for P_N - calN^(M) = sum_(k>=M) rho^-k int (omega_k - 1) eta_N(p^(k+1-N) y) p^-a r_k(p) p' dy, with
    omega_k <= omega_hi, r_k <= rbar = min(1/s, G s^(g-1)) (BK.4) and eta_N(p^j y) <= 2/(2+5^j y)."""
    a, gam = arb(al), _gam()
    rho = arb(5) ** a / 3
    five3 = arb(5) / 3
    pts = _grid(ylo, yhi, ratio)
    tot = arb(0)
    for i in range(len(pts) - 1):
        y0, y1 = pts[i], pts[i + 1]
        Y = y0.union(y1)
        p, pp = Y * (Y + 5), 2 * Y + 5
        rp = _min_up(1 / p, G * p ** (gam - 1))
        ry = _min_up(1 / Y, G * Y ** (gam - 1))
        beta = (Y + 6) / ((Y + 2) * (Y + 5))
        e = 2 * Y / (2 + Y)
        num = (1 + Y / 5) ** a
        U2 = arb((arb(3) / 5 * num * ry * (3 + p) / (1 + 2 * Y / 5)).upper())
        den = (1 + 2 * Y / 5) * five3 * (beta - e * rp)
        w = _min_up(num / den, U2) if den.lower() > 0 else U2
        wm1 = w - 1
        if wm1.upper() < 0:
            wm1 = arb(0)
        wt = arb((p ** (-a) * rp * pp).upper())
        S = arb(Sbound(y0, al, N, M).upper())
        tot += arb((wm1 * wt * S * (y1 - y0)).upper())
    y0 = pts[0]
    d0 = (4 + y0) * y0 / 30 + (1 + 2 * y0 / 5) * G * y0 ** gam
    c3 = (a * y0 ** (1 - gam) / 5 + (4 + y0) * y0 ** (1 - gam) / 30 + (1 + 2 * y0 / 5) * G) / (1 - d0)
    th = (1 + gam - a) / 2
    qq = 1 / (rho * arb(5) ** th)
    ex1 = 2 * gam - a - th
    ok = bool(qq.upper() < 1) and bool(ex1.lower() > 0)
    head = (c3 * G * arb(5) ** (gam - 1 - a) * (5 + 2 * y0) * (2 * arb(5) ** (N - 1)) ** th * qq ** M / (1 - qq)
            * y0 ** ex1 / ex1)
    Y0 = pts[-1]
    Kt = 3 * (arb(1) / 5 + 1 / Y0) ** a * (arb(1) / 2 + 5 / (4 * Y0)) * (2 + 5 / Y0)
    x = 1 / (5 * rho)
    tail = Kt * 2 * arb(5) ** (N - 1) * x ** M / (1 - x) * Y0 ** (-a - 1) / (a + 1)
    return arb((tot + head + tail).upper()), ok


def certify(al, N=1, M=M_TAIL, levels=LEVELS):
    ctx.prec = 64
    G, glo, ghi = arb(G_BI4), arb(G_LO), arb(G_HI)
    res = {"alpha": al}
    res["calN"] = {m: numer(al, N, m, G) for m in tuple(levels) + (M,)}
    res["D"] = {m: Dlevel(al, m) for m in levels}
    Dlo, Dhi, _ = Dinf(al, 60, G, glo, ghi)
    res["Dinf"] = (Dlo, Dhi)
    E, okE = Ebound(al, N, M, G)
    res["E"], res["okE"] = E, okE
    return res


def brackets(res, M=M_TAIL):
    V = {}
    for m, (Nlo, Nhi) in res["calN"].items():
        if m in res["D"]:
            Dlo, Dhi = res["D"][m]
            V[m] = (Nlo / Dhi, Nhi / Dlo)
    NMlo, NMhi = res["calN"][M]
    Dlo, Dhi = res["Dinf"]
    Phi = NMhi + res["E"]
    C = (NMlo / Dhi, Phi / Dlo)
    m6 = max(res["D"])
    allN = (NMlo / Dhi, Phi / res["D"][m6][0])   # for n >= M: calN^(M) <= calN^(n), D_n <= D_inf; calN^(n) <= P, D_n >= D_6
    every = (res["calN"][1][0] / Dhi, Phi / res["D"][1][0])   # for every n >= 1 (m = 1)
    return V, C, allN, every


# ----------------------------------------------------------------- (D) floats
def floats_V(al="0.45", ns=(1, 2, 3, 4, 6, 8, 12, 20), N=1):
    import mesh_maxwell_limit as ml
    from mpmath import mp, mpf, quad
    mp.dps = 20
    a = mpf(al)
    rho = mpf(5) ** a / 3
    X = mpf(10) ** 12
    out = []
    for n in ns:
        pts = [0] + [mpf(5) ** (-n) * mpf(10) ** j for j in range(-6, 40, 2) if mpf(5) ** (-n) * mpf(10) ** j < X] + [X]

        def pk(s, k):
            for _ in range(k):
                s = s * (s + 5)
            return s

        def etaN(zv):
            zs = [zv]
            for _ in range(N):
                zs.append(zs[-1] * (zs[-1] + 5))
            e = mpf(-1) / 2
            for k in range(N - 1, -1, -1):
                zz = zs[k]
                D = (zz + 2) * (zz + 5)
                e = (zz + 4) / D + (zz + 6) / D * e
            return e

        f = lambda s: s ** (-a) * ml.rec_all(n, s)[0]
        tail = X ** (-a) / a - 2 * X ** (-1 - a) / (1 + a) + 6 * X ** (-2 - a) / (2 + a)
        D = quad(f, pts) + tail
        num = quad(lambda s: f(s) * etaN(pk(s, n - N)), pts)
        out.append((n, float(num / D), float(rho ** (-n) * num / D)))
    return out


def main() -> int:
    t0 = time.time()
    print("Theorem BQ: 1 - delta_N^(n) is of exact order rho^n below log3/log5, at every level; Kept Failure BQ.1. Not thrust.")
    symbolic_checks()
    exact_checks()
    table = {}
    for al in ALPHAS:
        res = certify(al)
        V, C, allN, every = brackets(res)
        table[al] = (res, V, C, allN)
        gam = _gam()
        check(bool(arb(al) < gam) and res["okE"], f"alpha={al} < gamma; E_30 head exponents valid (q < 1, 2gamma-a-theta > 0)")
        Ns = [res["calN"][m] for m in sorted(res["calN"])]
        mono = all(Ns[i][1] < Ns[i + 1][0] for i in range(len(Ns) - 1))
        check(mono, f"alpha={al}: certified calN^(n) strictly increasing over n = 1..6, 30 (BQ.1): "
                    + ", ".join(f"[{float(lo.lower()):.6f},{float(hi.upper()):.6f}]" for lo, hi in Ns))
        Dl, Dh = res["Dinf"]
        print(f"   alpha={al}: D_inf in [{float(Dl.lower()):.7f},{float(Dh.upper()):.7f}], E_30 <= {float(res['E'].upper()):.3e}")
        print(f"   alpha={al}: V^(n), n=1..6: " + ", ".join(f"[{float(lo.lower()):.6f},{float(hi.upper()):.6f}]" for lo, hi in (V[m] for m in LEVELS)))
        check(bool(C[0].lower() > 0) and bool(C[0] < C[1]),
              f"alpha={al}: C_1 = lim rho^-n (1 - delta_1^(n)) in [{float(C[0].lower()):.7f}, {float(C[1].upper()):.7f}]; "
              f"for every n >= 30, rho^-n (1 - delta_1^(n)) in [{float(allN[0].lower()):.7f}, {float(allN[1].upper()):.7f}]; "
              f"for every n >= 1 in [{float(every[0].lower()):.6f}, {float(every[1].upper()):.6f}]")
    Dl, Dh = table["0.45"][0]["Dinf"]
    check(bool(Dl < arb("4.0013017734819")) and bool(arb("4.001301773482") < Dh),
          "consistency: the D_inf(0.45) bracket contains Theorem AS's certified [4.0013017734819, 4.001301773482]")
    # Kept Failure BQ.1
    _, V45, C45, _ = table["0.45"]
    kf1 = bool(V45[1][1] < V45[2][0]) and bool(V45[3][1] < V45[2][0]) and bool(V45[2][0] > C45[1])
    check(kf1, f"Kept Failure BQ.1 witness at alpha=0.45: V^(1) <= {float(V45[1][1].upper()):.6f} < V^(2) >= {float(V45[2][0].lower()):.6f} "
               f"> V^(3) (<= {float(V45[3][1].upper()):.6f}) and V^(2) > C_1 (<= {float(C45[1].upper()):.6f}): not monotone, approached from above")
    _, V30, C30, _ = table["0.3"]
    kf2 = bool(V30[3][1] < C30[0]) and bool(V30[4][0] > C30[1])
    check(kf2, f"Kept Failure BQ.1 witness at alpha=0.3: V^(3) <= {float(V30[3][1].upper()):.6f} < C_1 >= {float(C30[0].lower()):.6f}, "
               f"V^(4) >= {float(V30[4][0].lower()):.6f} > C_1 <= {float(C30[1].upper()):.6f}: crosses its limit")
    _, V01, C01, _ = table["0.01"]
    check(bool(V01[1][1] < C01[0]), f"alpha=0.01: V^(1) <= {float(V01[1][1].upper()):.6f} < C_1 (from below at n=1)")
    fl = floats_V("0.45")
    print("   floats alpha=0.45 N=1 (n, 1-delta, rho^-n(1-delta)): " + "; ".join(f"{n}: {d:.6e}, {v:.7f}" for n, d, v in fl))
    C = table["0.45"][2]
    V45b = table["0.45"][1]
    check(abs(fl[0][2] - float(V45b[1][0].mid())) < 2e-4 and abs(fl[1][2] - float(V45b[2][0].mid())) < 2e-4
          and float(C[0].lower()) < fl[-1][2] < float(C[1].upper()),
          "floats: quadrature V^(1), V^(2) at 0.45 agree with the certified enclosures, and V^(20) lies in the certified C_1 bracket "
          "(consistency only)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BQ checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
