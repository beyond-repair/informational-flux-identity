#!/usr/bin/env python3
r"""Theorem CI: the C^1 essential spectral radius of CE's transfer operator is rho(L_1), certified to 2e-7 (CH's open
annulus 0.1089 < |mu| <= 0.1345 closed down to 0.1343547 < |mu| <= 0.1343551).

Setting (GASKET.md, Theorems CE-CH): in s = sqrt(t + 25/4), X = [sqrt5/2, 5/2],
    (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s),   sigma_(1,2)(s) = sqrt(15/4 -+ s),
    wt1 = (s+1/2)/(2s(7/2-s)),  wt2 = (s-1/2)/(2s(s+7/2)),   |sigma_i'| = 1/(2 sigma_i),
    L_1 f = wt1 |sigma_1'| f(sigma_1) + wt2 |sigma_2'| f(sigma_2)   (CH's L_1).
sigma_1 X = Y_1 = [sqrt5/2, g_a], sigma_2 X = Y_2 = [g_b, 5/2], gap G = (g_a, g_b) (Lt never reads G).

Theorem CI (checked here; exit 1 on any failure):
  CI.1 [certified; arb, two-sided Collatz-Wielandt] An explicit polynomial h > 0 on X (degree 12, rational coefficients in
       u = (s - 9/5)/(7/10)) satisfies lam_lo h < L_1 h < lam_up h on X, lam_lo = 0.1343549, lam_up = 0.1343551
       (second-order arb Taylor cover).  Hence rho(L_1) in [lam_lo, lam_up] (CH.2 had only rho(L_1) <= 0.1345).
  CI.2 [exact; sympy] With c_i = |sigma_i'| h(sigma_i) / (L_1 h): wt1 c_1 + wt2 c_2 = 1, and 2 sigma_i c_i / h(sigma_i) =
       1/(L_1 h).  Gap profile: p(u) = (1-u)(1-3u) has p(0) = 1, p(1) = 0, int_0^1 p = 0, |p| <= 1,
       max |int_0^u p| = 4/27; psi(t) = 6t(1-t) has psi(0) = psi(1) = 0, int psi = 1, max psi = 3/2, Psi = int psi in [0,1].
  CI.3 [proved; arb] A right inverse with sharp derivative growth.  For f in C^1(X) put Kf(x) = c_i(T_i x) f(T_i x) on Y_i
       (T_i = sigma_i^-1), and on G (tau = (x - g_a)/l, l = g_b - g_a)
           (Kf)'(x) = m_a p(tau/delta) 1[tau<delta] + m_b p((1-tau)/delta) 1[tau>1-delta] + (V_b - V_a) psi(tau)/l,
       Kf(g_a) = V_a, where V_a, V_b, m_a, m_b are the one-sided values and slopes from Y_1, Y_2 (delta = 10^-6).  Then
       Kf in C^1(X), Lt K = I, and with N(f) = (sup|f|, sup|f'|/h) one has N(Kf) <= A N(f) entrywise with an explicit
       2x2 matrix A, A_11 <= kappa / lam_lo (kappa = 1 + O(delta)), and rho(A) < 1/0.1343547.
  CI.4 [proved from CI.3] Every mu with |mu| <= 0.1343547 is an eigenvalue of Lt on C^1(X): v = sum_n mu^n K^n g with
       CH's gap bump g converges in C^1, Lt v = mu v, and v != 0 (on the middle of G, v is a monic quartic plus a cubic).
       This is Kept Failure CH.4 extended from 0.1089 to 0.1343547.
  CI.5 [proved from CI.1, CI.4, CH.3 / Assumption CH] The essential spectral radius of Lt on C^1(X) satisfies
           0.1343547 <= r_ess(C^1) <= rho(L_1) <= 0.1343551,
       every C^1 eigenvalue mu != r has |mu| <= 0.1343551, and sup{|mu| : mu != r a C^1 eigenvalue} lies in
       [0.1343547, 0.1343551] (CH: [0.1089, 0.1345]); ratio to r in [0.3170752, 0.3170762].
Assumption CH (Hennion; literature) is used only for the upper bound r_ess <= rho(L_1), exactly as in CH.3.  CG.3's
r (to 1e-14) and CH.3's identification are imported as stated.  The level sums are untouched (CG.4 concerns the analytic
Lt^4 G; its rate 0.17255 and its non-explicit C stand).  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW enters.
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

import sympy as sp
from flint import arb, arb_series, ctx

FAILS: list[str] = []
DEG_H = 12
LAM_LO = "0.1343549"      # certified: L_1 h > LAM_LO h on X
LAM_UP = "0.1343551"      # certified: L_1 h < LAM_UP h on X
DELTA = Fr(1, 10 ** 6)    # width (in tau) of the two slope windows on G
MU_MAX = "0.1343547"      # every |mu| <= MU_MAX is a C^1 eigenvalue
LAM1_CH = "0.1345"        # CH.2's upper bound (superseded)
MU_CH = "0.1089"          # CH.4's disk (superseded)
R_LO, R_HI = "0.4237312869626", "0.4237312869627"   # r from CG.3 / CE.4
L2_ABS = "0.0731138175124"
U_C, U_R = Fr(9, 5), Fr(7, 10)    # u = (s - 9/5)/(7/10)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def _a(x: Fr) -> arb:
    return arb(x.numerator) / x.denominator


# ----------------------------------------------------------------- CI.2 exact identities
def symbolic_checks():
    s, u, t = sp.symbols("s u t", positive=True)
    H = sp.Function("h")
    half = sp.Rational(1, 2)
    W1 = (s + half) / (2 * s * (sp.Rational(7, 2) - s))
    W2 = (s - half) / (2 * s * (s + sp.Rational(7, 2)))
    S1, S2 = sp.sqrt(sp.Rational(15, 4) - s), sp.sqrt(sp.Rational(15, 4) + s)
    a1, a2 = -sp.diff(S1, s), sp.diff(S2, s)
    check(sp.simplify(a1 - 1 / (2 * S1)) == 0 and sp.simplify(a2 - 1 / (2 * S2)) == 0,
          "CI.2: |sigma_i'| = 1/(2 sigma_i) (sigma_1' < 0 < sigma_2')")
    L1h = W1 * a1 * H(S1) + W2 * a2 * H(S2)
    c1, c2 = a1 * H(S1) / L1h, a2 * H(S2) / L1h
    check(sp.simplify(W1 * c1 + W2 * c2 - 1) == 0, "CI.2: wt1 c_1 + wt2 c_2 = 1 for every positive h (so Lt K = I off G)")
    check(sp.simplify(2 * S1 * c1 / H(S1) - 1 / L1h) == 0 and sp.simplify(2 * S2 * c2 / H(S2) - 1 / L1h) == 0,
          "CI.2: 2 sigma_i c_i / h(sigma_i) = 1/(L_1 h): the f'-coefficient of (Kf)'/h is h/(L_1 h) <= 1/lam_lo")
    x = sp.symbols("x", positive=True)
    s0 = sp.sqrt(5) / 2
    ga, gb = sp.sqrt(sp.Rational(15, 4) - s0), sp.sqrt(sp.Rational(15, 4) + s0)
    T1, T2 = sp.Rational(15, 4) - x ** 2, x ** 2 - sp.Rational(15, 4)
    check(sp.simplify(T1.subs(x, s0) - sp.Rational(5, 2)) == 0 and sp.simplify(T1.subs(x, ga) - s0) == 0
          and sp.simplify(T2.subs(x, gb) - s0) == 0 and T2.subs(x, sp.Rational(5, 2)) == sp.Rational(5, 2),
          "CI.3: T_1 = 15/4 - x^2 maps Y_1 = [sqrt5/2, g_a] onto X, T_2 = x^2 - 15/4 maps Y_2 = [g_b, 5/2] onto X; "
          "both gap endpoints map to sqrt5/2")
    p = (1 - u) * (1 - 3 * u)
    P = sp.integrate(p, (u, 0, u))
    ok = p.subs(u, 0) == 1 and p.subs(u, 1) == 0 and sp.integrate(p, (u, 0, 1)) == 0
    ok = ok and sp.expand(P - u * (1 - u) ** 2) == 0
    crit = [r for r in sp.solve(sp.diff(P, u), u) if 0 <= r <= 1]
    ok = ok and max(abs(P.subs(u, r)) for r in crit + [0, 1]) == sp.Rational(4, 27)
    pc = [r for r in sp.solve(sp.diff(p, u), u) if 0 <= r <= 1]
    ok = ok and max(abs(p.subs(u, r)) for r in pc + [0, 1]) == 1
    psi = 6 * t * (1 - t)
    Psi = sp.integrate(psi, (t, 0, t))
    ok = ok and psi.subs(t, 0) == 0 and psi.subs(t, 1) == 0 and sp.integrate(psi, (t, 0, 1)) == 1
    ok = ok and psi.subs(t, half) == sp.Rational(3, 2) and sp.expand(Psi - (3 * t ** 2 - 2 * t ** 3)) == 0
    ok = ok and sp.expand(sp.diff(Psi, t) - psi) == 0 and sp.factor(psi) == -6 * t * (t - 1) and psi.subs(t, half) > 0   # psi >= 0 on [0,1]: roots 0, 1 only
    check(bool(ok), "CI.2: p = (1-u)(1-3u): p(0) = 1, p(1) = 0, int p = 0, |p| <= 1, max|int_0^u p| = max u(1-u)^2 = 4/27; "
                    "psi = 6t(1-t): psi(0) = psi(1) = 0, int psi = 1, max psi = 3/2, Psi = 3t^2 - 2t^3 increasing from 0 to 1")
    # the gap bump: as in CH.4
    g = ((x - sp.Symbol("A")) * (sp.Symbol("B") - x)) ** 2
    check(sp.Poly(g, x).degree() == 4 and sp.Poly(g, x).LC() == 1,
          "CI.4: CH's bump g = ((x - g_a)(g_b - x))^2 is a monic quartic on G (C^1 on X with g = 0 off G, CH.4)")


# ----------------------------------------------------------------- CI.1 two-sided Collatz-Wielandt
def _float_h():
    """Float proposal only: Chebyshev collocation eigenfunction of L_1, refit to degree DEG_H in u."""
    import numpy as np
    from numpy.polynomial import chebyshev as C
    a, b = 5 ** 0.5 / 2, 2.5

    def L1(f, x):
        w1 = (x + .5) / (2 * x * (3.5 - x))
        w2 = (x - .5) / (2 * x * (x + 3.5))
        y1, y2 = np.sqrt(3.75 - x), np.sqrt(3.75 + x)
        return w1 / (2 * y1) * f(y1) + w2 / (2 * y2) * f(y2)

    n = 40
    z = np.cos(np.pi * (np.arange(n) + .5) / n)
    X = (a + b) / 2 + (b - a) / 2 * z
    tm = lambda y: (2 * y - (a + b)) / (b - a)
    V = C.chebvander(z, n - 1)
    M = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = 1
        M[:, j] = L1(lambda y: C.chebval(tm(y), e), X)
    ev, vec = np.linalg.eig(np.linalg.solve(V, M))
    i = int(np.argmax(abs(ev)))
    co = vec[:, i].real
    co /= C.chebval(0, co)
    xx = np.linspace(a, b, 20001)
    uu = (xx - 1.8) / 0.7
    pc = np.polynomial.polynomial.polyfit(uu, C.chebval(tm(xx), co), DEG_H)
    return float(ev[i].real), [Fr(float(c)).limit_denominator(10 ** 15) for c in pc]


def _horner(coef, u):
    acc = arb(0) * u
    for c in reversed(coef):
        acc = acc * u + _a(c)
    return acc


def _hser(coef, x):
    return _horner(coef, (x - _a(U_C)) / _a(U_R))


def _L1h(coef, x):
    half = arb(1) / 2
    y1 = (arb(15) / 4 - x).sqrt()
    y2 = (arb(15) / 4 + x).sqrt()
    w1 = (x + half) / (2 * x * (arb(7) / 2 - x))
    w2 = (x - half) / (2 * x * (x + arb(7) / 2))
    return w1 * _hser(coef, y1) / (2 * y1) + w2 * _hser(coef, y2) / (2 * y2), y1, y2


def _cw_ok(coef, lam_lo, lam_up, lo, hi):
    """Second-order Taylor test of L_1 h - lam_lo h > 0, lam_up h - L_1 h > 0 and h > 0 on [lo, hi]."""
    m = arb(((lo + hi) / 2).mid())
    w = arb((((hi - lo) / 2).upper() + lo.rad() + hi.rad()))
    ball = arb(m.mid(), w.upper())
    out = []
    for xs, n in ((arb_series([m, 1], prec=3), 3), (arb_series([ball, 1], prec=3), 3)):
        Lh, _, _ = _L1h(coef, xs)
        hv = _hser(coef, xs)
        out.append(((Lh - lam_lo * hv).coeffs() + [arb(0)] * 3, (lam_up * hv - Lh).coeffs() + [arb(0)] * 3,
                    hv.coeffs() + [arb(0)] * 3))
    ok = True
    for k in range(3):
        e0, e1 = out[0][k][0], out[0][k][1]
        e2 = out[1][k][2]
        low = e0 - abs(e1) * w - abs(e2) * w * w
        ok = ok and bool(low > 0)
    return ok


def collatz_two_sided(maxdepth=30):
    ctx.prec = 160
    lam_f, coef = _float_h()
    lo_l, up_l = arb(LAM_LO), arb(LAM_UP)
    stack = [(arb(5).sqrt() / 2, arb(5) / 2, 0)]
    n_int, ok = 0, True
    while stack:
        lo, hi, dep = stack.pop()
        if _cw_ok(coef, lo_l, up_l, lo, hi):
            n_int += 1
            continue
        if dep >= maxdepth:
            ok = False
            break
        m = arb(((lo + hi) / 2).mid())
        stack.append((lo, m, dep + 1))
        stack.append((m, hi, dep + 1))
    check(ok, f"CI.1: h > 0 and {LAM_LO} h < L_1 h < {LAM_UP} h on X ({n_int} second-order arb intervals; degree-{DEG_H} "
              f"rational h), so rho(L_1) in [{LAM_LO}, {LAM_UP}] (float {lam_f:.10f}; CH.2: <= {LAM1_CH})")
    return lam_f, coef


# ----------------------------------------------------------------- CI.3 the majorant of K
def _cover_sup(fun, lo, hi, n):
    """sup over [lo, hi] of an arb-valued fun on balls, n uniform pieces; returns an arb upper bound."""
    best = arb(0)
    for j in range(n):
        a = lo + (hi - lo) * j / n
        b = lo + (hi - lo) * (j + 1) / n
        x = arb(((a + b) / 2).mid(), ((b - a) / 2).upper() + a.rad() + b.rad())
        v = fun(x)
        best = arb(max(best.upper(), v.upper()))
    return best


def _cover_inf(fun, lo, hi, n):
    best = None
    for j in range(n):
        a = lo + (hi - lo) * j / n
        b = lo + (hi - lo) * (j + 1) / n
        x = arb(((a + b) / 2).mid(), ((b - a) / 2).upper() + a.rad() + b.rad())
        v = fun(x).lower()
        best = v if best is None else min(best, v)
    return arb(best)


def _ci(coef, x, i):
    """(c_i(x), c_i'(x)) enclosures on the ball x, via a first-order arb series."""
    xs = arb_series([x, 1], prec=2)
    Lh, y1, y2 = _L1h(coef, xs)
    y = y1 if i == 1 else y2
    c = _hser(coef, y) / (2 * y) / Lh
    cc = c.coeffs() + [arb(0)] * 2
    yy = y.coeffs()[0]
    return cc[0], cc[1], yy


def k_majorant(coef, n_cover=4000):
    ctx.prec = 128
    s0, b = arb(5).sqrt() / 2, arb(5) / 2
    ga, gb = (arb(15) / 4 - s0).sqrt(), (arb(15) / 4 + s0).sqrt()
    l = gb - ga
    dl = _a(DELTA)
    h = lambda x: _hser(coef, x)
    A00Y = arb(0)
    A10Y = arb(0)
    for i in (1, 2):
        A00Y = arb(max(A00Y.upper(), _cover_sup(lambda x: _ci(coef, x, i)[0], s0, b, n_cover).upper()))

        def a10(x, i=i):
            c, dc, y = _ci(coef, x, i)
            return 2 * y * abs(dc) / h(y)
        A10Y = arb(max(A10Y.upper(), _cover_sup(a10, s0, b, n_cover).upper()))
    A11Y = 1 / arb(LAM_LO)                                  # sup h/(L_1 h) < 1/lam_lo by CI.1
    c1s0, _, _ = _ci(coef, arb(s0), 1)
    c2s0, _, _ = _ci(coef, arb(s0), 2)
    Jc = abs(c2s0 - c1s0)
    Vmax = arb(max(c1s0.upper(), c2s0.upper()))
    Hmax = arb(max(h(ga).upper(), h(gb).upper()))
    minG = _cover_inf(h, ga, gb, 2000)
    wa = _cover_inf(h, ga, ga + dl * l, 1)
    wb = _cover_inf(h, gb - dl * l, gb, 1)
    kap = arb(max((h(ga) / wa).upper(), (h(gb) / wb).upper()))
    A00 = arb(max(A00Y.upper(), (Vmax + arb(4) / 27 * dl * l * Hmax * A10Y).upper()))
    A01 = arb(4) / 27 * dl * l * Hmax * A11Y
    A10 = arb(max(A10Y.upper(), (kap * A10Y + arb(3) / 2 * Jc / (l * minG)).upper()))
    A11 = kap * A11Y
    tr, det = A00 + A11, A00 * A11 - A01 * A10
    rho = (tr + (tr * tr - 4 * det).sqrt()) / 2
    check(bool(minG > 0) and bool(tr * tr - 4 * det > 0) and bool(A00 < A11),
          f"CI.3: majorant A = [[{A00.str(6)}, {A01.str(4)}], [{A10.str(6)}, {A11.str(10)}]] (sup c_i <= {A00Y.str(6)}, "
          f"sup 2 sigma_i|c_i'|/h(sigma_i) <= {A10Y.str(6)}, kappa <= {kap.str(10)}, |c_2 - c_1|(s0) <= {Jc.str(6)}, "
          f"min_G h >= {minG.str(6)}, delta = 1e-6)")
    check(bool(arb(MU_MAX) * rho < 1),
          f"CI.3/CI.4: rho(A) <= {rho.str(12)}, 1/rho(A) >= {(1 / rho).str(10)}; {MU_MAX} * rho(A) = "
          f"{(arb(MU_MAX) * rho).str(10)} < 1, so sum mu^n K^n g converges in C^1(X) for |mu| <= {MU_MAX}")
    return rho


def consequences():
    ctx.prec = 128
    lo, up, mu = arb(LAM_LO), arb(LAM_UP), arb(MU_MAX)
    r_lo, r_hi = arb(R_LO), arb(R_HI)
    ok = bool(mu < lo) and bool(up < arb(LAM1_CH)) and bool(mu > arb(MU_CH)) and bool(up < r_lo)
    ok = ok and bool(mu > arb(L2_ABS) + arb("1e-13"))
    check(ok, f"CI.5: {MU_CH} < {MU_MAX} <= r_ess(C^1) <= rho(L_1) <= {LAM_UP} < {LAM1_CH} < r: CH's open annulus "
              f"(0.1089, 0.1345] shrinks to ({MU_MAX}, {LAM_UP}] (width {(up - mu).str(3)})")
    rl, rh = mu / r_hi, up / r_lo
    check(bool(rl > arb("0.3170752")) and bool(rh < arb("0.3170762")),
          f"CI.5: sup{{|mu|/r : mu != r a C^1 eigenvalue}} in [{rl.str(8)}, {rh.str(8)}] (CH: [0.2569, 0.3175])")


def main() -> int:
    t0 = time.time()
    print("Theorem CI: the C^1 essential spectral radius of CE's transfer operator is rho(L_1) in [0.1343547, 0.1343551]; "
          "every |mu| <= 0.1343547 is a C^1 eigenvalue. Not thrust.")
    symbolic_checks()
    _, coef = collatz_two_sided()
    k_majorant(coef)
    consequences()
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CI checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
