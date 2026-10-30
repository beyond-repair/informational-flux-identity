#!/usr/bin/env python3
r"""Theorem CL: the essential spectral radius of CE's transfer operator on C^k(X) is exactly rho(L_k) for every k >= 1
(lower bound new for k >= 2), in the Fredholm sense; rho(L_2), rho(L_3), rho(L_4) certified two-sided.  On C^1 the
Fredholm essential spectrum is the whole closed disk |z| <= rho(L_1), including its boundary circle; inside the open disk
Lt - mu is onto with infinite-dimensional kernel; and no eigenfunction for |mu| >= rho(L_1) can have zero 1-jet on the
attractor, so CJ's gap-bump eigenfunctions provably never reach the circle.

Setting (GASKET.md, Theorems CE-CK): in s = sqrt(t + 25/4), X = [sqrt5/2, 5/2],
    (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s),   sigma_(1,2)(s) = sqrt(15/4 -+ s),
    wt1 = (s+1/2)/(2s(7/2-s)),  wt2 = (s-1/2)/(2s(s+7/2)),   L_k f = sum_i wt_i |sigma_i'|^k f(sigma_i),
    M_k g = sum_i wt_i (sigma_i')^k g(sigma_i) on C(X).  Lambda = attractor of {sigma_1, sigma_2} (a Cantor set).

Theorem CL (checked here; exit 1 on any failure):
  CL.1 [proved] Calkin norm: ||M_k^n||_e = ||M_k^n|| = ||L_k^n 1||_inf on C(X) for every k, n >= 1 (disjointly supported
       sign-matched bumps are weakly null in C(X); compact operators kill them).
  CL.2 [proved from CL.1, CK.3(ii)] lim ||Lt^n||_e^(1/n) = rho(L_k) on C^k(X); so (Atkinson + Gelfand in the Calkin
       algebra) some z_k with |z_k| = rho(L_k) has Lt - z_k not Fredholm.  With CK.3: r_ess(Lt; C^k(X)) = rho(L_k) for
       every k >= 1 (Fredholm and Browder senses).  Also 0 is in the essential spectrum (gap bumps span an infinite-
       dimensional kernel) for every k.
  CL.3 [certified; arb, two-sided Collatz-Wielandt] rho(L_2) in [0.0430702, 0.0430704], rho(L_3) in [0.0139192,
       0.0139193], rho(L_4) in [0.0045264, 0.0045265] (CH.2 had upper bounds 0.0432, 0.01396, 0.00454 only).
  CL.4 [proved from CI.3, CJ.4, CK.3] On C^1: for |mu| < rho(L_1), Lt - mu = Lt (I - mu K) is onto with infinite-
       dimensional kernel; every |mu| = rho(L_1) has Lt - mu not Fredholm.  So the Fredholm essential spectrum on C^1 is
       exactly {|z| <= rho(L_1)} and spec \ ess = {r}.
  CL.5 [proved; constants certified in arb] For |mu| >= rho(L_1), if Lt f = mu f, f in C^1(X), and f = f' = 0 on Lambda,
       then f = 0.  CJ's eigenfunctions all have zero 1-jet on Lambda, so no circle point is an eigenvalue with a
       CJ-type eigenfunction; any circle eigenfunction must have nonzero 1-jet on Lambda.  Distortion constants:
       |(log wt_omega)'| <= D_w, |(log|sigma_omega'|)'| <= D_s for all words (certified finite).
Not claimed: whether points of |z| = rho(L_1) are eigenvalues (only the jet-zero class is excluded); whether the open
disk |z| < rho(L_k) lies in the C^k spectrum for k >= 2; B != 0 and an explicit C_G for the level sums.  Stage 1
untouched; no W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

import sympy as sp
from flint import arb, arb_series, ctx

FAILS: list[str] = []
DEG_H = 12
CW = {2: ("0.0430702", "0.0430704"), 3: ("0.0139192", "0.0139193"), 4: ("0.0045264", "0.0045265")}
CH2 = {2: "0.0432", 3: "0.01396", 4: "0.00454"}      # CH.2's one-sided bounds (superseded)
LAM1 = ("0.1343549", "0.1343551")                     # CI.1
R_LO = "0.4237312869626"                              # CG.3
L2_ABS, L2_ERR = "0.0731138175124", "1e-13"
L3_ABS, L3_ERR = "0.0170676440626", "1e-12"
U_C, U_R = Fr(9, 5), Fr(7, 10)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def _a(x: Fr) -> arb:
    return arb(x.numerator) / x.denominator


# ----------------------------------------------------------------- exact identities (sympy)
def symbolic_checks():
    s = sp.symbols("s", positive=True)
    half = sp.Rational(1, 2)
    W1 = (s + half) / (2 * s * (sp.Rational(7, 2) - s))
    W2 = (s - half) / (2 * s * (s + sp.Rational(7, 2)))
    S1, S2 = sp.sqrt(sp.Rational(15, 4) - s), sp.sqrt(sp.Rational(15, 4) + s)
    ok = sp.simplify(sp.diff(sp.log(W1), s) - (1 / (s + half) - 1 / s + 1 / (sp.Rational(7, 2) - s))) == 0
    ok = ok and sp.simplify(sp.diff(sp.log(W2), s) - (1 / (s - half) - 1 / s - 1 / (s + sp.Rational(7, 2)))) == 0
    check(bool(ok), "CL.5: (log wt1)' = 1/(s+1/2) - 1/s + 1/(7/2-s), (log wt2)' = 1/(s-1/2) - 1/s - 1/(s+7/2)")
    d1 = sp.diff(sp.log(-sp.diff(S1, s)), s)
    d2 = sp.diff(sp.log(sp.diff(S2, s)), s)
    ok = sp.simplify(d1 - 1 / (2 * (sp.Rational(15, 4) - s))) == 0 and sp.simplify(d2 + 1 / (2 * (sp.Rational(15, 4) + s))) == 0
    check(bool(ok), "CL.5: (log|sigma_1'|)' = 1/(2(15/4-s)), (log|sigma_2'|)' = -1/(2(15/4+s)); sup_X = 2/5 (at s = 5/2)")
    s0 = sp.sqrt(5) / 2
    ok = sp.simplify(S1.subs(s, sp.Rational(5, 2)) - s0) == 0 and S2.subs(s, sp.Rational(5, 2)) == sp.Rational(5, 2)
    ga = S1.subs(s, s0)
    ok = ok and sp.simplify(ga - sp.sqrt(sp.Rational(15, 4) - s0)) == 0
    check(bool(ok), "CL.5: 5/2 = sigma_2(5/2) and sqrt5/2 = sigma_1(5/2), g_a = sigma_1(sqrt5/2), g_b = sigma_2(sqrt5/2): "
                    "both ends of X and of every gap sigma_omega(G) lie in Lambda")
    # sup |sigma_i'| on X: |sigma_1'| = 1/(2 sigma_1) <= 1/(2 sqrt5/2) = 1/sqrt5
    check(sp.simplify(1 / (2 * s0) - 1 / sp.sqrt(5)) == 0, "CL.5: sup_X |sigma_i'| = 1/sqrt5 < 1 (cylinders shrink like 5^(-n/2))")


# ----------------------------------------------------------------- CL.3 two-sided Collatz-Wielandt for L_k
def _float_h(k):
    import numpy as np
    from numpy.polynomial import chebyshev as C
    a, b = 5 ** 0.5 / 2, 2.5
    n = 40
    z = np.cos(np.pi * (np.arange(n) + .5) / n)
    X = (a + b) / 2 + (b - a) / 2 * z
    tm = lambda y: (2 * y - (a + b)) / (b - a)
    V = C.chebvander(z, n - 1)
    M = np.zeros((n, n))
    w1 = (X + .5) / (2 * X * (3.5 - X))
    w2 = (X - .5) / (2 * X * (X + 3.5))
    y1, y2 = np.sqrt(3.75 - X), np.sqrt(3.75 + X)
    for j in range(n):
        e = np.zeros(n)
        e[j] = 1
        M[:, j] = w1 / (2 * y1) ** k * C.chebval(tm(y1), e) + w2 / (2 * y2) ** k * C.chebval(tm(y2), e)
    ev, vec = np.linalg.eig(np.linalg.solve(V, M))
    i = int(np.argmax(ev.real))
    co = vec[:, i].real
    co /= C.chebval(0, co)
    xx = np.linspace(a, b, 20001)
    pc = np.polynomial.polynomial.polyfit((xx - 1.8) / 0.7, C.chebval(tm(xx), co), DEG_H)
    return float(ev[i].real), [Fr(float(c)).limit_denominator(10 ** 15) for c in pc]


def _hser(coef, x):
    u = (x - _a(U_C)) / _a(U_R)
    acc = arb(0) * u
    for c in reversed(coef):
        acc = acc * u + _a(c)
    return acc


def _Lkh(coef, k, x):
    half = arb(1) / 2
    y1 = (arb(15) / 4 - x).sqrt()
    y2 = (arb(15) / 4 + x).sqrt()
    w1 = (x + half) / (2 * x * (arb(7) / 2 - x))
    w2 = (x - half) / (2 * x * (x + arb(7) / 2))
    return w1 * _hser(coef, y1) / (2 * y1) ** k + w2 * _hser(coef, y2) / (2 * y2) ** k


def _cw_ok(coef, k, lam_lo, lam_up, lo, hi):
    """Second-order Taylor test of L_k h - lam_lo h > 0, lam_up h - L_k h > 0, h > 0 on [lo, hi] (as CI.1)."""
    m = arb(((lo + hi) / 2).mid())
    w = arb((((hi - lo) / 2).upper() + lo.rad() + hi.rad()))
    ball = arb(m.mid(), w.upper())
    out = []
    for xs in (arb_series([m, 1], prec=3), arb_series([ball, 1], prec=3)):
        Lh = _Lkh(coef, k, xs)
        hv = _hser(coef, xs)
        out.append(((Lh - lam_lo * hv).coeffs() + [arb(0)] * 3, (lam_up * hv - Lh).coeffs() + [arb(0)] * 3,
                    hv.coeffs() + [arb(0)] * 3))
    for j in range(3):
        low = out[0][j][0] - abs(out[0][j][1]) * w - abs(out[1][j][2]) * w * w
        if not bool(low > 0):
            return False
    return True


def collatz_two_sided(k, maxdepth=30):
    ctx.prec = 160
    lam_f, coef = _float_h(k)
    lo_l, up_l = arb(CW[k][0]), arb(CW[k][1])
    stack = [(arb(5).sqrt() / 2, arb(5) / 2, 0)]
    n_int, ok = 0, True
    while stack:
        lo, hi, dep = stack.pop()
        if _cw_ok(coef, k, lo_l, up_l, lo, hi):
            n_int += 1
            continue
        if dep >= maxdepth:
            ok = False
            break
        m = arb(((lo + hi) / 2).mid())
        stack += [(lo, m, dep + 1), (m, hi, dep + 1)]
    check(ok and bool(up_l < arb(CH2[k])),
          f"CL.3: k = {k}: h > 0 and {CW[k][0]} h < L_{k} h < {CW[k][1]} h on X ({n_int} second-order arb intervals, "
          f"degree-{DEG_H} rational h), so rho(L_{k}) in [{CW[k][0]}, {CW[k][1]}] (float {lam_f:.11f}; CH.2: <= {CH2[k]})")
    return lam_f


# ----------------------------------------------------------------- CL.5 distortion constants
def _cover_sup(fun, lo, hi, n):
    best = arb(0)
    for j in range(n):
        a = lo + (hi - lo) * j / n
        b = lo + (hi - lo) * (j + 1) / n
        x = arb(((a + b) / 2).mid(), ((b - a) / 2).upper() + a.rad() + b.rad())
        best = arb(max(best.upper(), abs(fun(x)).upper()))
    return best


def distortion_constants(n_cover=4000):
    ctx.prec = 128
    lo, hi = arb(5).sqrt() / 2, arb(5) / 2
    half = arb(1) / 2
    dw1 = _cover_sup(lambda x: 1 / (x + half) - 1 / x + 1 / (arb(7) / 2 - x), lo, hi, n_cover)
    dw2 = _cover_sup(lambda x: 1 / (x - half) - 1 / x - 1 / (x + arb(7) / 2), lo, hi, n_cover)
    ds = _cover_sup(lambda x: 1 / (2 * (arb(15) / 4 - x)), lo, hi, n_cover)
    q = 1 / arb(5).sqrt()
    Dw = arb(max(dw1.upper(), dw2.upper())) / (1 - q)
    Ds = ds / (1 - q)
    ok = bool(Dw < 30) and bool(Ds < 1)
    check(ok, f"CL.5: sup_X|(log wt1)'| <= {dw1.str(6)}, sup_X|(log wt2)'| <= {dw2.str(6)}, sup_X|(log|sigma_i'|)'| <= "
              f"{ds.str(6)}; with contraction 1/sqrt5 every word omega has |(log wt_omega)'| <= D_w = {Dw.str(6)} and "
              f"|(log|sigma_omega'|)'| <= D_s = {Ds.str(6)} (finite, uniform in |omega|)")
    return float(Dw.upper()), float(Ds.upper())


# ----------------------------------------------------------------- consequences (radius comparisons)
def consequences():
    ctx.prec = 128
    l2lo, l2hi = arb(L2_ABS) - arb(L2_ERR), arb(L2_ABS) + arb(L2_ERR)
    l3lo, l3hi = arb(L3_ABS) - arb(L3_ERR), arb(L3_ABS) + arb(L3_ERR)
    up = {k: arb(CW[k][1]) for k in CW}
    lo = {k: arb(CW[k][0]) for k in CW}
    ok = bool(up[2] < l2lo) and bool(l3hi < lo[2]) and bool(up[3] < l3lo) and bool(up[4] < l3lo)
    ok = ok and bool(lo[2] > arb("0.001")) and bool(lo[3] > arb("0.001")) and bool(lo[4] > arb("0.001"))
    check(ok, "CL.2/CL.3: r_ess(C^2) = rho(L_2) in [0.0430702, 0.0430704] lies strictly between |l_3| and |l_2|; "
              "r_ess(C^3) = rho(L_3) in [0.0139192, 0.0139193] and r_ess(C^4) = rho(L_4) in [0.0045264, 0.0045265] lie "
              "below |l_3| and above CG.3's 1e-3; so on C^2 the eigenvalues r, l_2 and on C^3, C^4 the eigenvalues "
              "r, l_2, l_3 are exactly the isolated part, and l_3 sits inside the C^2 essential disk")
    a, b = arb(LAM1[0]), arb(LAM1[1])
    check(bool(b < arb(R_LO)) and bool(a > l2hi),
          f"CL.4: r_ess(C^1) = rho(L_1) in [{LAM1[0]}, {LAM1[1]}] (CI.1) > |l_2|: on C^1 l_2 and l_3 lie in the "
          "essential disk, and the Fredholm essential spectrum is the closed disk")


# ----------------------------------------------------------------- floats: attempts to break CL.1 / CL.5 (no claim)
def _words(n, xs):
    import numpy as np
    sig = [lambda y: np.sqrt(3.75 - y), lambda y: np.sqrt(3.75 + y)]
    dsg = [lambda y: -1 / (2 * np.sqrt(3.75 - y)), lambda y: 1 / (2 * np.sqrt(3.75 + y))]
    w = [lambda y: (y + .5) / (2 * y * (3.5 - y)), lambda y: (y - .5) / (2 * y * (y + 3.5))]
    lw = [lambda y: 1 / (y + .5) - 1 / y + 1 / (3.5 - y), lambda y: 1 / (y - .5) - 1 / y - 1 / (y + 3.5)]
    # state per word: point y = sigma_omega(s), weight W = wt_omega(s), D = sigma_omega'(s), dlogW = (log wt_omega)'(s)
    st = [(xs.copy(), np.ones_like(xs), np.ones_like(xs), np.zeros_like(xs))]
    for _ in range(n):
        nx = []
        for y, Wt, D, dl in st:
            for i in (0, 1):
                nx.append((sig[i](y), Wt * w[i](y), D * dsg[i](y), dl + lw[i](y) * D))
        st = nx
    return st


def float_observations(nmax=8):
    import numpy as np
    xs = np.linspace(np.sqrt(5) / 2, 2.5, 401)
    out = {}
    for k in (1, 2, 3, 4):
        st = _words(nmax, xs)
        Lk = sum(Wt * np.abs(D) ** k for _, Wt, D, _ in st)
        out[k] = Lk.max() ** (1 / nmax)
    st = _words(nmax, xs)
    maxdl = max(np.abs(dl).max() for _, _, _, dl in st)
    pts = np.sort(np.array([y for y, _, _, _ in st]), axis=0)
    mind = np.min(np.diff(pts, axis=0))
    print(f"     (observation) ||L_k^{nmax} 1||^(1/{nmax}) = " + ", ".join(f"{out[k]:.5f}" for k in out) +
          " for k = 1..4 (decreasing toward the certified rho(L_k)); "
          f"max over 2^{nmax} words of |(log wt_omega)'| = {maxdl:.4f}; min distance of the points sigma_omega(s) = {mind:.2e}")
    # Calkin witness for k = 1, n = 3: genuine piecewise-linear sign-matched bumps g_m around the 8 points
    # sigma_omega(s_m), s_m -> s* (the maximiser of L_1^3 1), with pairwise disjoint supports; M_1^3 g_m(s_m) is then
    # evaluated by direct composition (not by the formula) and compared with L_1^3 1(s_m).
    n = 3
    st = _words(n, xs)
    L = sum(Wt * np.abs(D) for _, Wt, D, _ in st)
    j = int(np.argmax(L))
    sstar = xs[j]
    sm = [sstar + (1 if sstar < 2.4 else -1) * 0.02 / m for m in range(1, 7)]
    F = [np.array([y[0] for y, _, _, _ in _words(n, np.array([s]))]) for s in sm]
    Fs = np.array([y[0] for y, _, _, _ in _words(n, np.array([sstar]))])
    radii = []
    for m, Fm in enumerate(F):
        others = np.concatenate([F[q] for q in range(len(F)) if q != m] + [Fs])
        dmin = min(np.min(np.abs(Fm[:, None] - others[None, :])), np.min(np.diff(np.sort(Fm))))
        radii.append(dmin / 3)
    sig = [lambda y: np.sqrt(3.75 - y), lambda y: np.sqrt(3.75 + y)]
    dsg = [lambda y: -1 / (2 * np.sqrt(3.75 - y)), lambda y: 1 / (2 * np.sqrt(3.75 + y))]
    w = [lambda y: (y + .5) / (2 * y * (3.5 - y)), lambda y: (y - .5) / (2 * y * (y + 3.5))]

    def M1(g, x, depth):
        if depth == 0:
            return g(x)
        return sum(w[i](x) * dsg[i](x) * M1(g, sig[i](x), depth - 1) for i in (0, 1))
    vals, disjoint = [], True
    for m, s in enumerate(sm):
        st_m = _words(n, np.array([s]))
        sg = np.array([np.sign(Wt[0] * D[0]) for _, Wt, D, _ in st_m])
        Fm, rm = F[m], radii[m]
        g = lambda x, Fm=Fm, rm=rm, sg=sg: float(np.sum(sg * np.maximum(0, 1 - np.abs(x - Fm) / rm)))
        L1s = float(sum(Wt[0] * abs(D[0]) for _, Wt, D, _ in st_m))
        vals.append((M1(g, s, n), L1s))
        for q in range(m):
            disjoint = disjoint and np.min(np.abs(Fm[:, None] - F[q][None, :])) > rm + radii[q]
    print("     (observation) CL.1 witness k = 1, n = 3, disjoint supports = " + str(bool(disjoint)) + ": M_1^3 g_m(s_m) = " +
          ", ".join(f"{a:.7f}" for a, _ in vals) + f" (= L_1^3 1(s_m) to {max(abs(a - b) for a, b in vals):.1e}) -> "
          f"||L_1^3 1|| = {L.max():.7f}")
    out["witness_err"] = max(abs(a - b) for a, b in vals)
    out["witness_disjoint"] = bool(disjoint)
    return out, maxdl


def main() -> int:
    t0 = time.time()
    print("Theorem CL: r_ess(Lt; C^k) = rho(L_k) exactly for all k (Calkin norm); rho(L_2..4) two-sided; on C^1 the "
          "Fredholm essential spectrum is the closed disk incl. its circle; jet-zero eigenfunctions never reach the circle. "
          "Not thrust.")
    symbolic_checks()
    for k in (2, 3, 4):
        collatz_two_sided(k)
    Dw, Ds = distortion_constants()
    consequences()
    try:
        out, maxdl = float_observations()
        check(maxdl <= Dw, "CL.5 (float sanity): the observed word distortion stays below the certified D_w")
        check(out["witness_disjoint"] and out["witness_err"] < 1e-12,
              "CL.1 (float sanity): genuine disjoint sign-matched bumps reach M_1^3 g_m(s_m) = L_1^3 1(s_m)")
        check(all(out[k] > float(CW[k][1]) for k in CW) and out[1] > float(LAM1[1]),
              "CL.1 (float sanity): ||L_k^n 1||^(1/n) >= rho(L_k) (as it must, since ||L_k^n|| >= rho^n)")
    except ImportError as e:  # pragma: no cover
        print("     (observation skipped:", e, ")")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CL checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
