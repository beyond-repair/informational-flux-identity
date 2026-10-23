#!/usr/bin/env python3
"""Theorem BS: the rate of rho^N C_N -> c_inf below log3/log5, and an exact characterisation of R, P_N and C_N.

Setting (GASKET.md, Theorems BG, BH, BJ, BQ, BR): p(z) = z(z+5), alpha_1 = (z+4)/((z+2)(z+5)), beta = (z+6)/((z+2)(z+5)),
eta_0 = -1/2, eta_N(z) = alpha_1(z) + beta(z) eta_(N-1)(p z), eta = lim eta_N, rho = 5^a/3, gamma = log3/log5, r_n the corner
resolvent, r_inf = lim r_n, R(u) = lim_M (3/5)^M r_M(5^-M u), psi(v) = lim p^k(5^-k v), D_inf = lim D_n, and for alpha < gamma,
N >= 1:  P_N = int u^-a R(u) eta_N(psi(5^-N u)) du,  C_N = P_N/D_inf (BQ.2),  c_inf = lim rho^N C_N (BR.4).

Theorem BS (checked here; exit 1 on any failure):
  (1) [scaled BJ.1; alpha-free] For every u > 0 the limit Q_inf(u) = lim_M (3/5)^M r_inf(5^-M u) exists (AT.1 is NOT used), and
          1/R(u) - 1/Q_inf(u) = Phi(u) := (1 + 2 eta(psi(u))) prod_(k>=1) (5/3) beta(psi(5^-k u))  in (0, 3].
      So R = Q_inf/(1 + Q_inf Phi) exactly.  Q_inf is nonincreasing, Q_inf(5u) = (3/5) Q_inf(u) (u^(1-gamma) Q_inf is
      multiplicatively 5-periodic), and g_lo u^(gamma-1) <= Q_inf(u) <= g_hi u^(gamma-1) for EVERY u > 0 (BJ.3).
  (2) [exact correction; every alpha in (0, gamma), N >= 1] c_inf - rho^N C_N = rho^N chi_N with chi_N = (A_N + B_N)/D_inf,
          A_N = int u^-a h(u) eta(psi(5^-N u)) du,  h = Q_inf - R = Q_inf^2 Phi/(1 + Q_inf Phi) > 0,
          B_N = int u^-a R(u) prod_(k=1..N) beta(psi(5^-k u)) (1/2 + eta(psi(u))) du.
      A_N strictly increases to J = int u^-a h du < inf (monotone convergence, eta(psi(5^-N u)) increases to eta(0) = 1), and
      0 < B_N <= (3/5)^(N-1) B_1.  Hence chi_N -> chi_inf = J/D_inf, and for every N >= 1
          A_1/D_inf <= chi_N <= (J + (3/5)^(N-1) B_1)/D_inf,      C_N = c_inf rho^-N - chi_inf + o(1):
      rho^N C_N increases to c_inf at the EXACT rate rho^N (not faster: chi_N >= A_1/D_inf > 0; not slower: chi_N bounded).
  (3) [certified, python-flint arb, outward] J, A_1, B_1, chi_inf and the every-N chi_N bracket at alpha = 0.01, 0.1, 0.3, 0.45;
      at alpha = 0.45: P_N, C_N, rho^N C_N and chi_N for N = 1..6 individually, c_inf D_inf, and the consistency checks against
      BQ.3's C_1 and BR.5's c_inf.
Assumptions used as stated (not re-run): AQ (3) (r_n increasing to r_inf; r_n, r_inf decreasing Stieltjes), BD.4 (psi_k
increasing to psi, psi(5v) = p(psi(v)), psi(v) >= v), BJ.6's psi(v) <= v e^(v/10), BG.3 (R = lim R_M, R <= 1/3), BH.2 (eta
recursion, eta = eta_N + prod beta (1/2 + eta o p^N), 0 < eta <= 2/(2+z), decreasing), BJ.1 (1/r_k - 1/r_inf =
prod_(i<k) beta(p^i s)(1 + 2 eta(p^k s))), BJ.3 (g = s^(1-gamma) r_inf in [g_lo, g_hi] on (0, p(1e-12)]), BL.1 (eta_N decreasing),
BQ.2 (C_N = P_N/D_inf, D_inf routine), BR.4 (rho^N P_N = int t^-a Q_N eta_N(psi t), c_inf D_inf = lim).
Floats: consistency only.  Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import sympy as sp
from flint import arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasket_fractional_deficit_depth as br  # noqa: E402
import gasket_fractional_deficit_rate_small_alpha as bq  # noqa: E402

FAILS: list[str] = []
ALPHAS = ("0.01", "0.1", "0.3", "0.45")
G_LO, G_HI = "0.600945014", "0.600963487"
NMAX = 6
EPS = "1e-20"
# BQ.3's certified C_1 and BR.5's certified c_inf (consistency only)
C1_BQ = br.C1_BQ
CINF_BR = {"0.01": ("0.010439", "0.010444"), "0.1": ("0.105440", "0.105526"),
           "0.3": ("0.327320", "0.327592"), "0.45": ("0.520387", "0.520824")}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    x, u, t, Q, F, a, N, M, m = sp.symbols("x u t Q Phi a N M m", positive=True)
    beta = (x + 6) / ((x + 2) * (x + 5))
    one_minus = sp.factor(sp.simplify(1 - sp.Rational(5, 3) * beta))
    check(sp.simplify(one_minus - x * (3 * x + 16) / (3 * x**2 + 21 * x + 30)) == 0,
          f"1 - (5/3) beta(x) = x(3x+16)/(3x^2+21x+30) in [0, x(3x+16)/30]: (5/3) beta in (0, 1], -log((5/3)beta(x)) = O(x), "
          "so prod_(k>=1) (5/3) beta(psi(5^-k u)) converges (psi(v) <= v e^(v/10))")
    check(sp.factor(sp.simplify(1 / x - beta)) == sp.factor((x + 10) / (x * (x + 2) * (x + 5))),
          "1/x - beta(x) = (x+10)/(x(x+2)(x+5)) > 0: beta(x) <= 1/x (tails)")
    # scaled BJ.1: with s = 5^-M u, (5/3)^M / r(s) = 1/((3/5)^M r(s)) and the product index k = M - i
    rho = 5**a / 3
    check(sp.simplify(sp.powsimp(sp.expand_power_base((sp.Rational(3, 5))**M * (5**(-M) * u)**(sp.log(3) / sp.log(5) - 1),
                                                      force=True), force=True) - u**(sp.log(3) / sp.log(5) - 1)) == 0,
          "(3/5)^M (5^-M u)^(gamma-1) = u^(gamma-1): BJ.3's bracket for r_inf near 0 becomes g_lo u^(gamma-1) <= Q_inf(u) <= "
          "g_hi u^(gamma-1) for every u > 0")
    check(sp.simplify((sp.Rational(3, 5))**(m + 1) - sp.Rational(3, 5) * (sp.Rational(3, 5))**m) == 0,
          "(3/5)^(m+1) r_inf(5^-m u) = (3/5) (3/5)^m r_inf(5^-m u): Q_inf(5u) = (3/5) Q_inf(u)")
    # h = Q - R with R = Q/(1+Q Phi)
    h = Q - Q / (1 + Q * F)
    check(sp.simplify(h - Q**2 * F / (1 + Q * F)) == 0, "Q - Q/(1+Q Phi) = Q^2 Phi/(1+Q Phi)")
    dQ, dF = sp.factor(sp.diff(Q**2 * F / (1 + Q * F), Q)), sp.factor(sp.diff(Q**2 * F / (1 + Q * F), F))
    check(sp.simplify(dQ - Q * F * (Q * F + 2) / (1 + Q * F)**2) == 0 and sp.simplify(dF - Q**2 / (1 + Q * F)**2) == 0,
          "d/dQ = Q Phi (Q Phi + 2)/(1+Q Phi)^2 > 0, d/dPhi = Q^2/(1+Q Phi)^2 > 0: h increasing in Q and Phi (bracket ends); "
          "R = 1/(1/Q + Phi) increasing in Q, decreasing in Phi")
    # scaling u = 5^N t
    tt = 5**(-N) * u
    fac = sp.simplify(sp.powsimp(sp.expand_power_base(rho**(-N) * tt**(-a) * sp.diff(tt, u) * (sp.Rational(5, 3))**N,
                                                      force=True), force=True))
    check(sp.simplify(fac - u**(-a)) == 0,
          "rho^-N t^-a (5/3)^N dt = u^-a du at t = 5^-N u: rho^-N int t^-a Q_inf(t) eta(psi t) dt = int u^-a Q_inf(u) "
          "eta(psi(5^-N u)) du (Q_inf(5^-N u) = (5/3)^N Q_inf(u)) and rho^-N int t^-a Q_N(t) eta_N(psi t) dt = P_N "
          "(Q_N(5^-N u) = (5/3)^N R(u))")
    # decomposition
    Rr, e, eN = sp.symbols("R e e_N")
    check(sp.expand(Q * e - Rr * eN - ((Q - Rr) * e + Rr * (e - eN))) == 0,
          "Q_inf eta - R eta_N = (Q_inf - R) eta + R (eta - eta_N): rho^-N (c_inf D_inf - rho^N P_N) = A_N + B_N")
    # tail-sum geometric constants of the Phi lower bound
    j = sp.symbols("j", integer=True, positive=True)
    check(sp.summation(sp.Rational(1, 5)**j, (j, 1, sp.oo)) == sp.Rational(1, 4) and
          sp.summation(sp.Rational(1, 25)**j, (j, 1, sp.oo)) == sp.Rational(1, 24),
          "sum_(j>=1) 5^-j = 1/4, sum 25^-j = 1/24 (tail k > K of the Phi lower bound, prod(1-c_k) >= 1 - sum c_k)")
    # head of J: int_0^eps u^(c-1) = eps^c/c
    c = sp.Rational(1, 3)
    check(sp.integrate(u**(c - 1), (u, 0, sp.Rational(1, 100))) == sp.Rational(1, 100)**c / c,
          "int_0^eps u^(c-1) du = eps^c/c (heads; c = gamma - a or 1 - a)")
    # eta_N(0) = 1 - (9/10)(3/5)^(N-1) -> 1 = eta(0)
    ok = all(bq.eta_rec(n, Fr(0)) == 1 - Fr(9, 10) * Fr(3, 5)**(n - 1) for n in range(1, 9))
    check(ok, "eta_N(0) = 1 - (9/10)(3/5)^(N-1) (N = 1..8, exact), so eta(0) = 1 and eta(psi(5^-N u)) increases to 1")


# ----------------------------------------------------------------- (B) exact rationals
def _p(x):
    return x * (x + 5)


def exact_checks():
    """The product structure behind (1): for any two levels, 1/r_M(s) - 1/r_(M+L)(s) = prod_(i<M) beta(p^i s) *
    (1/r_0(p^M s) - 1/r_L(p^M s)) (AQ's affine reciprocal step, the same mechanism as BJ.1)."""
    ok = True
    for M in range(0, 4):
        for L in range(1, 3):
            for s in (Fr(1, 1000), Fr(1, 7), Fr(2), Fr(19, 3)):
                lhs = 1 / bq.r_rec(M, s) - 1 / bq.r_rec(M + L, s)
                prod, z = Fr(1), s
                for _ in range(M):
                    prod *= (z + 6) / ((z + 2) * (z + 5))
                    z = _p(z)
                ok = ok and lhs == prod * (1 / bq.r_rec(0, z) - 1 / bq.r_rec(L, z))
    check(ok, "exact rationals M = 0..3, L = 1..2, s in {1/1000, 1/7, 2, 19/3}: 1/r_M - 1/r_(M+L) = prod_(i<M) beta(p^i s) "
              "(1/r_0 - 1/r_L)(p^M s)")
    ok2 = all(bq.r_rec(k, Fr(0)) == Fr(5, 3)**k / 3 for k in range(0, 6))
    check(ok2, "exact rationals: r_k(0) = (5/3)^k/3, k = 0..5, so R_M(0) = 1/3 and Phi(0) = 3 (R(0+) = 1/3)")
    xs = (Fr(0), Fr(1, 100), Fr(1), Fr(5), Fr(200))
    check(all(Fr(5, 3) * (x + 6) / ((x + 2) * (x + 5)) <= 1 for x in xs) and
          all((x + 6) / ((x + 2) * (x + 5)) <= 1 / x for x in xs if x > 0),
          "exact rationals: (5/3) beta(x) <= 1 and beta(x) <= 1/x")


# ----------------------------------------------------------------- (C) arb
def _eta_bracket(z):
    """eta_K(z) <= eta(z) <= eta_K(z) + prod_(i<K) beta(p^i z)(1/2 + 2/(2+p^K z)) (BH.2 tail form, eta <= 2/(2+z)), with K
    adaptive: at least 14 and 8 steps past p^k z > 1000 (so small z are as tight as large ones)."""
    o = [z]
    while len(o) < 15 or len([q for q in o if q > 1000]) < 8:
        o.append(_p(o[-1]))
    K = len(o) - 1
    e = arb(-1) / 2
    prod = arb(1)
    for k in range(K - 1, -1, -1):
        zz = o[k]
        D = (zz + 2) * (zz + 5)
        e = (zz + 4) / D + (zz + 6) / D * e
    for k in range(K):
        zz = o[k]
        prod *= (zz + 6) / ((zz + 2) * (zz + 5))
    return e, e + prod * (arb(1) / 2 + 2 / (2 + o[K]))


def _point(u, glo, ghi, gam, nmax=NMAX):
    """Brackets at one grid point u (an exact arb point).  Returns dict of (lo, hi) for the integrand ingredients."""
    K = nmax
    tiny = arb("1e-12")
    while (u / arb(5) ** K).upper() > tiny.upper():
        K += 1
    v = u / arb(5) ** K
    w = v * (v / 10).exp()
    zlo, zhi = [None] * (K + 1), [None] * (K + 1)
    zlo[K], zhi[K] = v, w                       # psi(v) in [v, v e^(v/10)]
    for k in range(K - 1, -1, -1):              # psi(5 x) = p(psi(x))
        zlo[k], zhi[k] = _p(zlo[k + 1]), _p(zhi[k + 1])
    # eta(psi(5^-k u)) by the exact recursion eta = alpha_1 + beta eta o p, started at k = 0
    e_up = [arb(_eta_bracket(zlo[0])[1].upper())]
    e_dn = [arb(_eta_bracket(zhi[0])[0].lower())]
    en_up, en_dn = [arb(-1) / 2], [arb(-1) / 2]    # eta_k(psi(5^-k u)), eta_0 = -1/2
    b_up, b_dn = [arb(1)], [arb(1)]                # prod_(j=1..k) beta(psi(5^-j u))
    for k in range(1, nmax + 1):
        for (zz, eu, nu, bu, side) in ((zlo[k], e_up, en_up, b_up, "u"), (zhi[k], e_dn, en_dn, b_dn, "d")):
            D = (zz + 2) * (zz + 5)
            a1, be = (zz + 4) / D, (zz + 6) / D
            ne, nn, nb = a1 + be * eu[-1], a1 + be * nu[-1], bu[-1] * be
            f = (lambda q: arb(q.upper())) if side == "u" else (lambda q: arb(q.lower()))
            eu.append(f(ne))
            nu.append(f(nn))
            bu.append(f(nb))
    five3 = arb(5) / 3
    Pu, Pd = arb(1), arb(1)
    for k in range(1, K + 1):
        Pu *= five3 * (zlo[k] + 6) / ((zlo[k] + 2) * (zlo[k] + 5))
        Pd *= five3 * (zhi[k] + 6) / ((zhi[k] + 2) * (zhi[k] + 5))
    S = (16 * w / 4 + 3 * w * w / 24) / 30       # sum_(k>K) x_k(3x_k+16)/30
    Phi_hi = arb(((1 + 2 * e_up[0]) * Pu).upper())
    Phi_lo = arb(((1 + 2 * e_dn[0]) * Pd * (1 - S)).lower())
    x = u ** (gam - 1)
    Qlo, Qhi = arb((glo * x).lower()), arb((ghi * x).upper())
    h_hi = arb((Qhi * Qhi * Phi_hi / (1 + Qhi * Phi_hi)).upper())
    h_lo = arb((Qlo * Qlo * Phi_lo / (1 + Qlo * Phi_lo)).lower())
    R_hi = arb((1 / (1 / Qhi + Phi_lo)).upper())
    R_lo = arb((1 / (1 / Qlo + Phi_hi)).lower())
    return dict(Q=(Qlo, Qhi), Phi=(Phi_lo, Phi_hi), h=(h_lo, h_hi), R=(R_lo, R_hi), zlo=zlo,
                eta=(e_dn, e_up), etaN=(en_dn, en_up), bprod=(b_dn, b_up))


def _grid_pieces(T):
    pts = bq._grid(EPS, "1e-6", "1.003") + bq._grid("1e-6", "1e8", "1.0002")[1:] + bq._grid("1e8", T, "1.003")[1:]
    return pts
    return pts


def certify(al, T="1e10"):
    """Lower/upper sums (all true integrands are positive and nonincreasing in u) plus closed-form heads and tails."""
    ctx.prec = 64
    a, gam = arb(al), bq._gam()
    glo, ghi = arb(G_LO), arb(G_HI)
    pts = _grid_pieces(T)
    Ns = range(1, NMAX + 1)
    keys = ["J", "M"] + [f"A{N}" for N in Ns] + [f"B{N}" for N in Ns] + [f"P{N}" for N in Ns]
    lo = {k: arb(0) for k in keys}
    hi = {k: arb(0) for k in keys}
    prev = None
    first = None
    for i, u in enumerate(pts):
        d = _point(u, glo, ghi, gam)
        ua = u ** (-a)
        e_dn, e_up = d["eta"]
        n_dn, n_up = d["etaN"]
        b_dn, b_up = d["bprod"]
        vals = {"J": (d["h"][0], d["h"][1]), "M": (d["Q"][0] * e_dn[0], d["Q"][1] * e_up[0])}
        for N in Ns:
            vals[f"A{N}"] = (d["h"][0] * e_dn[N], d["h"][1] * e_up[N])
            vals[f"B{N}"] = (d["R"][0] * b_dn[N] * (arb(1) / 2 + e_dn[0]), d["R"][1] * b_up[N] * (arb(1) / 2 + e_up[0]))
            vals[f"P{N}"] = (d["R"][0] * n_dn[N], d["R"][1] * n_up[N])
        cur = {k: (arb((ua * l).lower()), arb((ua * h).upper())) for k, (l, h) in vals.items()}
        if prev is not None:
            wdt = u - pts[i - 1]
            for k in keys:
                lo[k] += arb((cur[k][0] * wdt).lower())     # right endpoint: lower
                hi[k] += arb((prev[k][1] * wdt).upper())    # left endpoint: upper
        else:
            first = d
        prev = cur
        last = d
    e0, Tt = pts[0], pts[-1]
    c1, c2 = gam - a, 1 - a
    # heads on (0, eps]
    head_Q_hi = ghi * e0 ** c1 / c1
    head_h_lo = glo * e0 ** c1 / c1 - e0 ** c2 / (3 * c2)        # h = Q - R >= Q - 1/3
    hi["J"] += head_Q_hi
    lo["J"] += head_h_lo
    hi["M"] += head_Q_hi
    lo["M"] += first["eta"][0][0] * glo * e0 ** c1 / c1
    for N in Ns:
        hi[f"A{N}"] += head_Q_hi
        lo[f"A{N}"] += first["eta"][0][N] * head_h_lo
        hi[f"B{N}"] += (arb(3) / 5) ** N / 2 * e0 ** c2 / c2      # R <= 1/3, prod beta <= (3/5)^N, 1/2 + eta <= 3/2
        hi[f"P{N}"] += e0 ** c2 / (3 * c2)
    # tails on [T, inf)
    tailJ = 3 * ghi ** 2 * (arb(25) / 9) * 125 * Tt ** (2 * gam - 3 - a) / (3 + a - 2 * gam)
    hi["J"] += tailJ
    # psi(v)/v nondecreasing (psi_k(v)/v = prod_(j<k)(p^j(5^-k v)+5)/5), so for u >= T: psi(5^-k u) >= (u/T) psi(5^-k T),
    # and eta(z), eta_N(z) <= 2/(2+z) <= 2/z, beta(z) <= 1/z; with R, Q_inf <= g_hi u^(gamma-1) each tail is
    # <= g_hi (2 or 3/2) T/psi_lo(5^-k T) * T^(gamma-1-a)/(1+a-gam)
    zT = last["zlo"]
    tailpow = ghi * Tt * Tt ** (gam - 1 - a) / (1 + a - gam)
    hi["M"] += 2 * tailpow / zT[0]
    for N in Ns:
        hi[f"A{N}"] += tailJ
        hi[f"B{N}"] += tailpow * (arb(3) / 2) / zT[1]
        hi[f"P{N}"] += 2 * tailpow / zT[N]
    G = arb(bq.G_BI4)
    Dlo, Dhi, _ = bq.Dinf(al, 60, G, glo, ghi)
    out = {k: (arb(lo[k].lower()), arb(hi[k].upper())) for k in keys}
    out["D"] = (Dlo, Dhi)
    out["npts"] = len(pts)
    return out


def brackets(res):
    """chi_inf, the every-N bracket (from N0 on) and the individual enclosures for N = 1..NMAX."""
    Jl, Jh = res["J"]
    Dl, Dh = res["D"]
    chi_inf = (Jl / Dh, Jh / Dl)

    def every(N0):   # for every N >= N0: A_N >= A_N0 (A increasing), B_N <= B_N0 (B decreasing)
        return (res[f"A{N0}"][0] / Dh, (Jh + res[f"B{N0}"][1]) / Dl)

    indiv = {}
    for N in range(1, NMAX + 1):
        Pl, Ph = res[f"P{N}"]
        Al, Ah = res[f"A{N}"]
        Bl, Bh = res[f"B{N}"]
        indiv[N] = dict(C=(Pl / Dh, Ph / Dl), chi=((Al + Bl) / Dh, (Ah + Bh) / Dl), AB=(Al + Bl, Ah + Bh), P=(Pl, Ph))
    return chi_inf, every, indiv


def _f(x, lo=True, d=6):
    return f"{float(x.lower() if lo else x.upper()):.{d}f}"


# ----------------------------------------------------------------- (D) floats
def floats_chi(al="0.45", Ns=(1, 2, 3)):
    """Quadrature of chi_N and C_N with Q_inf ~ <Theta> u^(gamma-1) (AU's mean, floats only)."""
    from mpmath import inf, log, mp, mpf, quad
    mp.dps = 15
    gam = log(3) / log(5)
    g = mpf("0.6009542511")
    a = mpf(al)

    def psi(v, K=40):
        y = v / mpf(5) ** K
        for _ in range(K):
            y = y * (y + 5)
        return y

    def eta(z, N=60, inf_tail=True):
        zs = [z]
        for _ in range(N):
            if inf_tail and zs[-1] > 1e30:
                break
            zs.append(zs[-1] * (zs[-1] + 5))
        e = mpf(0) if (inf_tail and len(zs) <= N) else mpf(-1) / 2
        for k in range(len(zs) - 2, -1, -1):
            zz = zs[k]
            D = (zz + 2) * (zz + 5)
            e = (zz + 4) / D + (zz + 6) / D * e
        return e

    beta = lambda x: (x + 6) / ((x + 2) * (x + 5))

    def Phi(u):
        pr = 1 + 2 * eta(psi(u))
        for k in range(1, 50):
            pr *= mpf(5) / 3 * beta(psi(u / mpf(5) ** k))
        return pr

    R = lambda u: 1 / (u ** (1 - gam) / g + Phi(u))
    pts = [0] + [mpf(10) ** j for j in range(-8, 26, 3)] + [inf]
    Dinf = mpf("4.0013017734819") if al == "0.45" else None
    rho = mpf(5) ** a / 3
    cinfD = quad(lambda u: u ** (-a) * g * u ** (gam - 1) * eta(psi(u)), pts)
    J = quad(lambda u: u ** (-a) * (g * u ** (gam - 1) - R(u)), pts)
    out = []
    for N in Ns:
        P = quad(lambda u: u ** (-a) * R(u) * eta(psi(u / mpf(5) ** N), N, inf_tail=False), pts)
        out.append((N, float(P / Dinf), float(rho ** N * P / Dinf), float((cinfD - rho ** N * P) / (rho ** N * Dinf))))
    return float(cinfD / Dinf), float(J / Dinf), out


def main() -> int:
    t0 = time.time()
    print("Theorem BS: below gamma, c_inf - rho^N C_N = rho^N chi_N with chi_N -> chi_inf in (0, inf): rho^N C_N -> c_inf at the "
          "exact rate rho^N; 1/R = 1/Q_inf + Phi. Not thrust.")
    symbolic_checks()
    exact_checks()
    table = {}
    for al in ALPHAS:
        res = certify(al)
        rho = arb(5) ** arb(al) / 3
        chi_inf, every, indiv = brackets(res)
        table[al] = (res, chi_inf, every, indiv)
        Jl, Jh = res["J"]
        Dl, Dh = res["D"]
        Ml, Mh = res["M"]
        print(f"   alpha={al}: {res['npts']} grid points; J in [{_f(Jl)},{_f(Jh, False)}], A_1 >= {_f(res['A1'][0])}, "
              f"B_1 <= {_f(res['B1'][1], False)}, c_inf D_inf in [{_f(Ml)},{_f(Mh, False)}], D_inf in [{_f(Dl)},{_f(Dh, False)}]")
        e1, e6 = every(1), every(NMAX)
        check(bool(res["A1"][0] > 0) and bool(chi_inf[0] > 0) and bool(chi_inf[0] < chi_inf[1]) and bool(e6[0] < chi_inf[0])
              and bool(chi_inf[1] < e6[1]),
              f"alpha={al}: chi_inf = J/D_inf in [{_f(chi_inf[0])}, {_f(chi_inf[1], False)}]; every N >= 1: chi_N in "
              f"[{_f(e1[0])}, {_f(e1[1], False)}]; every N >= {NMAX}: [{_f(e6[0])}, {_f(e6[1], False)}] (contains chi_inf)")
        cl, ch = (arb(v) for v in CINF_BR[al])
        check(bool(Ml / Dh < ch) and bool(Mh / Dl > cl),
              f"alpha={al}: consistency, c_inf = D_inf^-1 int u^-a Q_inf eta(psi u) in [{_f(Ml / Dh)},{_f(Mh / Dl, False)}] "
              f"meets BR.5's [{CINF_BR[al][0]},{CINF_BR[al][1]}]")
        c1l, c1h = (arb(v) for v in C1_BQ[al])
        C1 = indiv[1]["C"]
        check(bool(C1[0] < c1h) and bool(C1[1] > c1l),
              f"alpha={al}: C_1 = P_1/D_inf in [{_f(C1[0], d=7)},{_f(C1[1], False, 7)}] meets BQ.3's independent "
              f"[{C1_BQ[al][0]},{C1_BQ[al][1]}] (level-30 numerator route)")
        ok_id = True
        for N in range(1, NMAX + 1):
            Pl, Ph = indiv[N]["P"]
            ABl, ABh = indiv[N]["AB"]
            ok_id = ok_id and bool(Ml < rho**N * (Ph + ABh)) and bool(rho**N * (Pl + ABl) < Mh)
        check(ok_id, f"alpha={al}: identity c_inf D_inf = rho^N (P_N + A_N + B_N), N = 1..{NMAX}: independent enclosures overlap")
        check(all(bool(indiv[N]["P"][1] * rho**N < indiv[N + 1]["P"][0] * rho**(N + 1)) for N in range(1, NMAX)),
              f"alpha={al}: certified rho^N P_N strictly increasing over N = 1..{NMAX} (consistency with BR.4)")
        print(f"   alpha={al}: (N, C_N, rho^N C_N, chi_N) certified: " + "; ".join(
            f"{N}: [{_f(indiv[N]['C'][0])},{_f(indiv[N]['C'][1], False)}], [{_f(rho**N * indiv[N]['C'][0])},"
            f"{_f(rho**N * indiv[N]['C'][1], False)}], [{_f(indiv[N]['chi'][0])},{_f(indiv[N]['chi'][1], False)}]"
            for N in range(1, NMAX + 1)))
    # certified decrease of chi_N at 0.45 for N = 1..5 (common D_inf, so compare A_N + B_N)
    ind45 = table["0.45"][3]
    dec = [N for N in range(1, NMAX) if bool(ind45[N]["AB"][0] > ind45[N + 1]["AB"][1])]
    check(dec == [1, 2, 3, 4],
          "alpha=0.45: A_N + B_N certified strictly decreasing for N = 1..5 (so chi_N decreases there); N = 5 vs 6 not "
          "separated at this resolution (Conjecture BS.4: chi_N decreases for every N)")
    fl = floats_chi()
    print(f"   floats alpha=0.45 (<Theta> for Q_inf): c_inf={fl[0]:.7f}, chi_inf={fl[1]:.7f}; (N, C_N, rho^N C_N, chi_N): "
          + "; ".join(f"{N}: {c:.7f}, {rc:.7f}, {x:.7f}" for N, c, rc, x in fl[2]))
    _, chi45, _, _ = table["0.45"]
    check(float(chi45[0].lower()) < fl[1] < float(chi45[1].upper()) and abs(fl[2][0][1] - 0.0719251) < 1e-5,
          "floats: chi_inf(0.45) inside the certified bracket and C_1(0.45) = 0.0719251 (consistency only)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BS checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
