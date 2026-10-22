#!/usr/bin/env python3
"""Theorem BF: the 3^-n constant c(alpha) of Theorems BD/BE is strictly convex on (0,1) (Conjecture BE.2 proved).

Setting: Assumption BA / BC (README), notation of Theorems BB-BE.  T, Y = scaling limits of BD.4,
c(a) = int tau^-a T Y dtau / int tau^-a Y dtau.  Put u = log tau, eps = 1 - a, P = (T-2) Y, and fix a shift s0.
With v = u - s0:
    c(a) - 2 = F(eps) := eps N(eps) / (1/3 + eps G(eps)),
    N(eps) = int e^(eps v) P dv,   G(eps) = int e^(eps v) (Y - (1/3) 1[v<0]) dv
(e^(eps s0) cancels; the 1/(3 eps) of int_{v<0} e^(eps v)/3 dv is pulled out).  N, G are smooth on a neighbourhood
of [0,1] (P <= tau/18, |Y - 1/3| <= C tau near 0, Y <= tau/(3 psi_6) at infinity), so F is smooth on [0,1] and
c''(a) = F''(1-a), with F'' = (Phi'' - 2F'h' - F h'')/h, Phi = eps N, h = 1/3 + eps G.

Theorem BF (checked here; exit 1 on any failure):
  (1) the representation above, for every shift s0 (sympy for the quotient and shift identities).
  (2) F''(eps) > 0 for every eps in [0,1]: certified on 40 closed eps-intervals of width 1/40 (mpmath.iv, eps as
      an interval, outward cell hulls on u in [-20, log 3125], explicit head and tail bounds); smallest certified
      lower bound 0.52046.  Hence c is strictly convex on (0,1) and extends C^2 to [0,1] with c'' > 0 there
      (Conjecture BE.2 proved).  Behind it: c'(a) = -Cov(T, u), c''(a) = kappa(T, u, u) (third joint cumulant)
      under nu_a ~ tau^-a Y dtau; no general principle fixes the sign of kappa, so (2) is a certificate, not a
      structural argument.
  (3) corollaries (convexity + endpoint data): for every a in (0,1)
        max(2 + 3S(1-a), c_0 + c'(0+) a)  <  c(a)  <  2 + (c_0 - 2)(1-a),
      with 3S = F'(0) (BE.3), c_0 = 2 + F(1) (BE.4), certified c'(0+) = -F'(1) in [-4.341285, -4.227515], so
        max(2 + 2.70632481 (1-a), 5.48553067 - 4.341285 a) < c(a) < 2 + 3.50249013 (1-a).
  Conjecture BF.2 (floats only): c'' is not monotone; it rises from c''(1-) ~ 1.531 to ~ 1.593 near a ~ 0.46
      and falls to c''(0+) ~ 1.557.  No closed form for T, Y, c, S, c_0 or c''.

Not a field solver.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import sympy as sp
from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_noncentrality_limit as bd  # noqa: E402

FAILS: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


CQ = "(e^0.1 + 0.3)/9"   # |Y - 1/3| <= CQ tau for tau <= 1


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    e, s0, d = sp.symbols("e s0 d", positive=True)
    Nf, Gf = sp.Function("N")(e), sp.Function("G")(e)
    F = e * Nf / (sp.Rational(1, 3) + e * Gf)
    Phi, h = e * Nf, sp.Rational(1, 3) + e * Gf
    F1 = (sp.diff(Phi, e) - F * sp.diff(h, e)) / h
    F2 = (sp.diff(Phi, e, 2) - 2 * F1 * sp.diff(h, e) - F * sp.diff(h, e, 2)) / h
    check(sp.simplify(sp.diff(F, e) - F1) == 0 and sp.simplify(sp.diff(F, e, 2) - F2) == 0,
          "F = eps N/(1/3 + eps G): F' = (Phi' - F h')/h, F'' = (Phi'' - 2F'h' - F h'')/h, Phi = eps N, h = 1/3 + eps G")
    # shift: eps int e^(eps u) Y du = e^(eps s0) (1/3 + eps G_s0), int_{-inf}^0 e^(eps v)/3 dv = 1/(3 eps)
    v = sp.Symbol("v", real=True)
    check(sp.simplify(e * sp.integrate(sp.exp(e * v) / 3, (v, -sp.oo, 0)) - sp.Rational(1, 3)) == 0,
          "eps int_{v<0} e^(eps v)/3 dv = 1/3, so eps D(eps) = e^(eps s0)(1/3 + eps G(eps)) for any shift s0")
    # c'' = third joint cumulant under nu ~ e^(eps u) m(du)
    # moments M_j = int e^(eps u) u^j m, MT_j = int e^(eps u) u^j T m; c = MT0/M0
    M = [sp.Function(f"M{j}")(e) for j in range(3)]
    MT = [sp.Function(f"MT{j}")(e) for j in range(3)]
    subs = {sp.Derivative(M[0], e): M[1], sp.Derivative(M[1], e): M[2],
            sp.Derivative(MT[0], e): MT[1], sp.Derivative(MT[1], e): MT[2]}
    c = MT[0] / M[0]
    dc = sp.diff(c, e).subs(subs)
    d2c = sp.diff(dc, e).subs(subs).subs(subs)
    ET, Eu, Eu2, ETu, ETu2 = MT[0] / M[0], M[1] / M[0], M[2] / M[0], MT[1] / M[0], MT[2] / M[0]
    cov = ETu - ET * Eu
    kap = ETu2 - ET * Eu2 - 2 * Eu * ETu + 2 * ET * Eu**2
    check(sp.simplify(dc - cov) == 0 and sp.simplify(d2c - kap) == 0,
          "d<T>/d eps = Cov(T,u), d^2<T>/d eps^2 = kappa(T,u,u) (eps = 1 - a: c'(a) = -Cov, c''(a) = kappa)")
    # head integral: int_{-inf}^{m} |v|^k e^v dv = e^m sum_j k!/j! |m|^j, m < 0
    mm = sp.Symbol("mm", positive=True)
    ok = True
    for kk in range(4):
        I1 = sp.integrate((-v) ** kk * sp.exp(v), (v, -sp.oo, -mm))
        I2 = sp.exp(-mm) * sum(sp.factorial(kk) / sp.factorial(j) * mm**j for j in range(kk + 1))
        ok = ok and sp.simplify(I1 - I2) == 0
    check(ok, "head: int_{-inf}^{-m} |v|^k e^v dv = e^-m sum_{j<=k} k!/j! m^j  (k = 0..3)")
    # |Y - 1/3| bound near 0: Y >= e^-x/(3+tau e^x), x = tau/10 (BE.3); 1/3 - that = (3 + tau e^x - 3e^-x)/(3(3+tau e^x))
    t, x = sp.symbols("t x", positive=True)
    expr = sp.Rational(1, 3) - sp.exp(-x) / (3 + t * sp.exp(x))
    check(sp.simplify(expr - (3 + t * sp.exp(x) - 3 * sp.exp(-x)) / (3 * (3 + t * sp.exp(x)))) == 0,
          "1/3 - e^-x/(3+tau e^x) = (tau e^x + 3(1-e^-x))/(3(3+tau e^x)) <= (tau e^0.1 + 3x)/9 = tau (e^0.1+0.3)/9  (tau <= 1)")


# ----------------------------------------------------------------- (B) cells
UMIN = -20
DU_DEN = 500          # cell width 1/500 in u
N_LEVEL = 60
TAU1 = 3125


def build_cells(den=DU_DEN, umin=UMIN, N=N_LEVEL):
    """Nodes u_i = umin + i/den up to the first node >= log 3125.  Per cell: hulls of P = (T-2)Y and Y."""
    iv.dps = 20
    imax = math.ceil((math.log(TAU1) - umin) * den) + 1
    nodes = []
    for i in range(imax + 1):
        u = iv.mpf(umin) + iv.mpf(i) / den
        tau = iv.exp(u)
        T, Y = bd.TY(N, tau)
        eT = iv.exp(bd.epsT(N, tau.b))
        eY = iv.exp(bd.epsY(N, tau.b))
        nodes.append((u, T, Y, eT, eY))
    cells = []
    for i in range(imax):
        ua, Ta, Ya, _, _ = nodes[i]
        ub, Tb, Yb, eTb, eYb = nodes[i + 1]
        # T increasing with T_N <= T <= T_N e^epsT;  Y decreasing with |log Y - log Y_N| <= epsY (eps increasing)
        Tlo, Thi = Ta.a, (Tb * eTb).b
        Ylo, Yhi = (Yb / eYb).a, (Ya * eYb).b
        if Thi < 2:
            Thi = iv.mpf(2).b
        P = iv.mpf([max((iv.mpf(Tlo) - 2).a, 0), (iv.mpf(Thi) - 2).b]) * iv.mpf([Ylo, Yhi])
        cells.append((i, iv.mpf([ua.a, ub.b]), P, iv.mpf([Ylo, Yhi])))
    return cells, imax


def tail_consts(M=6, K=8):
    """Y <= tau^(1-K) TAU1^K/(3 P_K) on tau >= TAU1 (BD.4); return pre = TAU1^K/(3 P_K)."""
    xx = sp.Symbol("x")
    q = xx
    for _ in range(M):
        q = sp.expand(q * (q + 5))
    q = sp.Poly(q.subs(xx, xx / sp.Integer(5) ** M), xx)
    Pk = sp.Rational(q.eval(sp.Integer(TAU1))) - sum(q.coeff_monomial(xx**j) * sp.Integer(TAU1) ** j for j in range(1, K))
    Pki = iv.mpf(int(Pk.p)) / iv.mpf(int(Pk.q))
    return iv.mpf(TAU1) ** K / (3 * Pki), K


def moments(cells, E, s0_num, s0_den, kmax, den=DU_DEN, umin=UMIN, tail=None):
    """Enclosures of N_k(E) = int e^(eps v) v^k P dv and G_k(E) = int e^(eps v) v^k (Y - 1/3 [v<0]) dv, k < kmax,
    for eps in the interval E (subset of [0,1]), v = u - s0, s0 = s0_num/s0_den a node."""
    s0 = iv.mpf(s0_num) / s0_den
    i0 = (s0_num * den) // s0_den - umin * den      # index of the node u = s0
    assert (s0_num * den) % s0_den == 0
    du = iv.mpf(1) / den
    third = iv.mpf(1) / 3
    N = [iv.mpf(0)] * kmax
    G = [iv.mpf(0)] * kmax
    for i, U, P, Y in cells:
        V = U - s0
        w = iv.exp(E * V) * du
        Q = (Y - third) if i < i0 else Y
        for k in range(kmax):
            vk = iv.mpf(1) if k == 0 else (V if k == 1 else V**k)   # V**2 is the exact even-power hull
            N[k] += w * vk * P
            G[k] += w * vk * Q
    # head u < UMIN: 0 <= P <= tau/18, |Q| <= CQ tau, e^(eps v) <= 1, |v|^k e^u = e^s0 |v|^k e^v
    m = s0 - umin
    CQ_i = (iv.exp(iv.mpf(1) / 10) + iv.mpf(3) / 10) / 9
    for k in range(kmax):
        H = iv.exp(s0) * iv.exp(-m) * sum(iv.mpf(math.factorial(k)) / math.factorial(j) * m**j for j in range(k + 1))
        N[k] += iv.mpf([-1, 1]) * H / 18
        G[k] += iv.mpf([-1, 1]) * H * CQ_i
    # tail tau >= TAU1 (covers the cells' end): e^(eps v) <= tau e^-s0, 0 < v^k <= (log tau)^k <= k! tau,
    # (T-2) <= tau/6, Y <= pre tau^(1-K);  du = dtau/tau
    pre, K = tail
    t1 = iv.mpf(TAU1)
    for k in range(kmax):
        fk = iv.mpf(math.factorial(k)) * iv.exp(-s0) * pre
        N[k] += iv.mpf([0, 1]) * fk / 6 * t1 ** (4 - K) / (K - 4)
        G[k] += iv.mpf([0, 1]) * fk * t1 ** (3 - K) / (K - 3)
    return N, G


def F_derivs(E, N, G, order=2):
    h = [1 / iv.mpf(3) + E * G[0], G[0] + E * G[1], 2 * G[1] + E * G[2]]
    Ph = [E * N[0], N[0] + E * N[1], 2 * N[1] + E * N[2]]
    if not h[0].a > 0:
        return None
    F0 = Ph[0] / h[0]
    F1 = (Ph[1] - F0 * h[1]) / h[0]
    F2 = (Ph[2] - 2 * F1 * h[1] - F0 * h[2]) / h[0]
    return F0, F1, F2


def shift_for(e_mid):
    """s0 = 2.25 + 0.25 floor((5/4) e_mid / 0.25 + 1/2) (a node: multiple of 1/4); chosen once from floats."""
    q = int(math.floor(1.25 * e_mid / 0.25 + 0.5))
    return 9 + q, 4          # (2.25 + q/4) = (9 + q)/4


def convexity_certificate(cells, tail, nint=40, verbose=True):
    lows = []
    for j in range(nint):
        E = iv.mpf([mpf(j) / nint, mpf(j + 1) / nint])
        sn, sd = shift_for((j + 0.5) / nint)
        N, G = moments(cells, E, sn, sd, 3, tail=tail)
        r = F_derivs(E, N, G)
        lo = r[2].a if r is not None else -mp.inf
        lows.append((j, lo, r[2].b if r is not None else mp.inf, sn / sd))
        if verbose:
            print(f"   eps in [{j}/{nint},{j+1}/{nint}] (a in [{1-(j+1)/nint:.3f},{1-j/nint:.3f}]), s0={sn/sd}: "
                  f"c'' in [{float(lows[-1][1]):.5f}, {float(lows[-1][2]):.5f}]")
    return lows


def point(cells, tail, eps, s0=(9, 4)):
    E = iv.mpf(eps)
    N, G = moments(cells, E, s0[0], s0[1], 3, tail=tail)
    return F_derivs(E, N, G)


# ----------------------------------------------------------------- (C) floats (display only)
def float_cpp():
    """c''(a) by double-precision Riemann sums of the cumulant formula (Conjecture BF.2 display, not a proof)."""
    import numpy as np
    u = np.linspace(-45, math.log(2e4), 200001)
    tau = np.exp(u)
    Nl = 40
    s = tau / 5.0**Nl
    o = [s]
    with np.errstate(over="ignore", invalid="ignore"):
        for _ in range(Nl):
            o.append(o[-1] * (o[-1] + 5))
        r = 1 / (3 + o[Nl]); t = np.ones_like(tau); yr = np.ones_like(tau)
        for k in range(Nl - 1, -1, -1):
            z = o[k]
            r = np.where(np.isfinite(z), (2 + z) * (5 + z) * r / ((6 + z) - 2 * z * (5 + z) * r), 0.0)
            t = (2 + z) * t / (6 + z - 2 * t)
            yr = yr / (5 + z)
    T, Y = 3.0**Nl * t, 3.0**Nl * np.nan_to_num(yr * r)
    P, Q = (T - 2) * Y, np.where(u < 0, Y - 1 / 3, Y)
    out = []
    for e in (0.0, 0.25, 0.54, 0.75, 1.0):
        w = np.exp(e * u)
        Nk = [np.sum(w * u**k * P) for k in range(3)]
        Gk = [np.sum(w * u**k * Q) for k in range(3)]
        du = u[1] - u[0]
        Nk = [x * du for x in Nk]; Gk = [x * du for x in Gk]
        h = [1 / 3 + e * Gk[0], Gk[0] + e * Gk[1], 2 * Gk[1] + e * Gk[2]]
        Ph = [e * Nk[0], Nk[0] + e * Nk[1], 2 * Nk[1] + e * Nk[2]]
        F0 = Ph[0] / h[0]; F1 = (Ph[1] - F0 * h[1]) / h[0]; F2 = (Ph[2] - 2 * F1 * h[1] - F0 * h[2]) / h[0]
        out.append((1 - e, F2))
    return out


# ----------------------------------------------------------------- main
def main(nint=40) -> int:
    t0 = time.time()
    print("Theorem BF: c(alpha) is strictly convex on (0,1) (Conjecture BE.2). Not thrust.")
    symbolic_checks()
    t1 = time.time()
    cells, nc = build_cells()
    tail = tail_consts()
    print(f"   {nc} cells on u in [-20, {UMIN + nc / DU_DEN:.3f}] built in {time.time()-t1:.1f}s")

    t1 = time.time()
    lows = convexity_certificate(cells, tail, nint=nint)
    mn = min(lo for _, lo, _, _ in lows)
    print(f"   certificate time {time.time()-t1:.1f}s; min certified lower bound for c'' on (0,1): {float(mn):.5f}")
    check(mn > 0, f"Theorem BF.2: F''(eps) = c''(1-eps) > 0 on all {nint} eps-intervals covering [0,1]: c strictly convex")

    # endpoint and consistency data
    r0 = point(cells, tail, 0)
    r1 = point(cells, tail, 1, s0=(14, 4))
    r55 = point(cells, tail, mpf("0.55"), s0=(12, 4))
    F1_0, F0_1, F1_1, F0_55 = r0[1], r1[0], r1[1], r55[0]
    print(f"   F'(0) = 3S in [{float(F1_0.a):.6f}, {float(F1_0.b):.6f}]; c_0 = 2 + F(1) in [{float((2+F0_1).a):.6f}, "
          f"{float((2+F0_1).b):.6f}]; F'(1) = -c'(0+) in [{float(F1_1.a):.6f}, {float(F1_1.b):.6f}]; "
          f"c(0.45) in [{float((2+F0_55).a):.6f}, {float((2+F0_55).b):.6f}]")
    check(F1_0.a < 2.70903118 and F1_0.b > 2.70632481, "consistency: F'(0) meets BE.3's certified 3S")
    check((2 + F0_1).a < 5.50249013 and (2 + F0_1).b > 5.48553067, "consistency: 2 + F(1) meets BE.4's certified c_0")
    check((2 + F0_55).a < 3.72914181 and (2 + F0_55).b > 3.72198897, "consistency: 2 + F(0.55) meets BD's certified c(0.45)")
    # corollary check at the BD points (both sides must be consistent with the certified c values)
    S3lo, c0hi = mpf("2.70632481"), mpf("5.50249013")
    ok = True
    for a, clo, chi in (("0.45", "3.72198897", "3.72914181"), ("0.8", "2.57126273", "2.57352573")):
        a = mpf(a)
        ok = ok and 2 + S3lo * (1 - a) < mpf(chi) and mpf(clo) < 2 + (c0hi - 2) * (1 - a)
    check(ok, "corollary BF.3 (tangent at a=1 below, chord above) is consistent with BD's certified c(0.45), c(0.8)")
    print(f"   BF.3 certified: for all a in (0,1), max(2 + 2.70632481 (1-a), 5.48553067 - {float(F1_1.b):.6f} a) < c(a) "
          f"< 2 + 3.50249013 (1-a);  c'(0+) in [{float(-F1_1.b):.6f}, {float(-F1_1.a):.6f}], c'(1-) = -3S")
    check(F1_1.a > 4 and F1_1.b < 4.35, "c'(0+) = -F'(1) certified in a window of width < 0.12")
    fl = float_cpp()
    print("   float c''(a): " + ", ".join(f"a={a:.2f}: {v:.4f}" for a, v in fl) + "  (Conjecture BF.2: not monotone)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BF checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
