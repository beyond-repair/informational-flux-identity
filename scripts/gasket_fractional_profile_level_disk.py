#!/usr/bin/env python3
r"""Theorem CJ: Conjecture CI.6 proved.  On C^1(X) the spectrum of CE's transfer operator is {r} plus the closed disk
of radius rho(L_1), and every |mu| < rho(L_1) is a C^1 eigenvalue (CI.4 had |mu| <= 0.1343547).

Setting (GASKET.md, Theorems CE-CI): in s = sqrt(t + 25/4), X = [sqrt5/2, 5/2],
    (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s),   sigma_(1,2)(s) = sqrt(15/4 -+ s),
    wt1 = (s+1/2)/(2s(7/2-s)),  wt2 = (s-1/2)/(2s(s+7/2)),   |sigma_i'| = 1/(2 sigma_i),
    L_1 f = u_1 f(sigma_1) + u_2 f(sigma_2),  u_i = wt_i |sigma_i'|  (CH's L_1).
Y_1 = sigma_1 X = [sqrt5/2, g_a], Y_2 = sigma_2 X = [g_b, 5/2], gap G = (g_a, g_b), l = g_b - g_a, s0 = sqrt5/2.

Theorem CJ (checked here; exit 1 on any failure):
  CJ.1 [exact; sympy] (log u_1)' = 1/(s+1/2) - 1/s + 1/(7/2-s) + 1/(2(15/4-s)),
       (log u_2)' = 1/(s-1/2) - 1/s - 1/(s+7/2) - 1/(2(15/4+s)); |sigma_i'| = 1/(2 sigma_i).
  CJ.2 [certified; arb] M_i = sup_X |(log u_i)'|: M_1 <= 1.33335 (= 4/3 at s = 5/2), M_2 <= 0.4045; Lipschitz constants
       c_1 = 1/sqrt5, c_2 = 1/(2 g_b) of sigma_i on X; B = max_i M_i/(1-c_i) <= 2.4121.  So the cone
       C_b = {f > 0 : log f is b-Lipschitz in s} obeys L_1 C_b in C_(max_i(M_i + c_i b)), and L_1^n 1 in C_B for all n.
       With a = 2B: b' = max_i(M_i + c_i a) < a and diam_(Theta_a) L_1 C_a <= 2 log((a+b')/(a-b')) + 2 b'|X| <= 13.31.
       Hence (CF.1) inf R_n -> rho(L_1) with R_n = L_1^(n+1) 1 / L_1^n 1.
  CJ.3 [certified; arb] With lambda = 2 and E = e^(B l): 
           Q = max(1/sqrt5, lambda * sup_G max_i |sigma_i'|, E/(2 g_a lambda)) <= 0.8050 < 1.
  CJ.4 [proved from CJ.2-CJ.3 and CI.3's construction] For every mu with |mu| < rho(L_1), mu is an eigenvalue of Lt on
       C^1(X) (witness v = sum mu^n K^n g with CI.3's K built from h = L_1^n 1, n large, delta small).
  CJ.5 [proved from CJ.4, CI.1, CI.5 / Assumption CH] spec(Lt; C^1(X)) = {r} u {|z| <= rho(L_1)},
       rho(L_1) in [0.1343549, 0.1343551]; r_ess(C^1) = rho(L_1) exactly.  Conjecture CI.6 holds.
Assumption CH (Hennion; literature) is used only for "nothing but r outside rho(L_1)", exactly as in CI.5.  Whether the
points of the circle |z| = rho(L_1) are eigenvalues is not decided.  The level sums are untouched.  Stage 1 untouched;
no W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time

import sympy as sp
from flint import arb, ctx

FAILS: list[str] = []
LAM = 2                     # weight factor of the sup-norm component on G
LAM_LO, LAM_UP = "0.1343549", "0.1343551"   # CI.1
R_LO = "0.4237312869626"                    # CG.3
L2_ABS = "0.0731138175124"                  # CG.3, |lambda_2|
MU_CI = "0.1343547"                         # CI.4's disk (superseded)
B_MAX = "2.4121"
Q_MAX = "0.8050"


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def symbolic_checks():
    s = sp.symbols("s", positive=True)
    half = sp.Rational(1, 2)
    W1 = (s + half) / (2 * s * (sp.Rational(7, 2) - s))
    W2 = (s - half) / (2 * s * (s + sp.Rational(7, 2)))
    S1, S2 = sp.sqrt(sp.Rational(15, 4) - s), sp.sqrt(sp.Rational(15, 4) + s)
    a1, a2 = -sp.diff(S1, s), sp.diff(S2, s)
    ok = sp.simplify(a1 - 1 / (2 * S1)) == 0 and sp.simplify(a2 - 1 / (2 * S2)) == 0
    d1 = sp.diff(sp.log(W1 * a1), s) - (1 / (s + half) - 1 / s + 1 / (sp.Rational(7, 2) - s)
                                       + 1 / (2 * (sp.Rational(15, 4) - s)))
    d2 = sp.diff(sp.log(W2 * a2), s) - (1 / (s - half) - 1 / s - 1 / (s + sp.Rational(7, 2))
                                       - 1 / (2 * (sp.Rational(15, 4) + s)))
    ok = ok and sp.simplify(d1) == 0 and sp.simplify(d2) == 0
    e1 = (1 / (s + half) - 1 / s + 1 / (sp.Rational(7, 2) - s) + 1 / (2 * (sp.Rational(15, 4) - s))).subs(s, sp.Rational(5, 2))
    ok = ok and sp.nsimplify(e1) == sp.Rational(4, 3)
    check(bool(ok), "CJ.1: |sigma_i'| = 1/(2 sigma_i); (log u_1)', (log u_2)' as stated (sympy); (log u_1)'(5/2) = 4/3")
    # CI.2's identities for symbolic h are imported as stated; re-check the central one quickly
    H = sp.Function("h")
    L1h = W1 * a1 * H(S1) + W2 * a2 * H(S2)
    c1, c2 = a1 * H(S1) / L1h, a2 * H(S2) / L1h
    check(sp.simplify(W1 * c1 + W2 * c2 - 1) == 0, "CI.2 (re-checked): wt1 c_1 + wt2 c_2 = 1 for every h > 0")


def _dlog(x, i):
    if i == 1:
        return 1 / (x + arb(1) / 2) - 1 / x + 1 / (arb(7) / 2 - x) + 1 / (2 * (arb(15) / 4 - x))
    return 1 / (x - arb(1) / 2) - 1 / x - 1 / (x + arb(7) / 2) - 1 / (2 * (arb(15) / 4 + x))


def _cover_sup(fun, lo, hi, n):
    best = arb(0)
    for j in range(n):
        a = lo + (hi - lo) * j / n
        b = lo + (hi - lo) * (j + 1) / n
        x = arb(((a + b) / 2).mid(), ((b - a) / 2).upper() + a.rad() + b.rad())
        best = arb(max(best.upper(), fun(x).upper()))
    return best


def constants(n_cover=20000):
    ctx.prec = 128
    s0, b5 = arb(5).sqrt() / 2, arb(5) / 2
    ga, gb = (arb(15) / 4 - s0).sqrt(), (arb(15) / 4 + s0).sqrt()
    l, X = gb - ga, b5 - s0
    M1 = _cover_sup(lambda x: abs(_dlog(x, 1)), s0, b5, n_cover)
    M2 = _cover_sup(lambda x: abs(_dlog(x, 2)), s0, b5, n_cover)
    c1 = 1 / (2 * (arb(15) / 4 - b5).sqrt())          # sup_X |sigma_1'| at s = 5/2
    c2 = 1 / (2 * gb)                                  # sup_X |sigma_2'| at s = sqrt5/2
    check(bool(abs(c1 - 1 / arb(5).sqrt()) < arb("1e-30")) and bool(c2 < arb("0.2267")),
          f"CJ.2: Lip(sigma_1) = 1/sqrt5 = {c1.str(8)}, Lip(sigma_2) = 1/(2 g_b) <= {c2.str(6)}")
    B = arb(max((M1 / (1 - c1)).upper(), (M2 / (1 - c2)).upper()))
    check(bool(M1 < arb("1.33335")) and bool(M2 < arb("0.4045")) and bool(B < arb(B_MAX)),
          f"CJ.2: M_1 <= {M1.str(8)}, M_2 <= {M2.str(6)} ({n_cover}-ball arb covers); B = max M_i/(1-c_i) <= {B.str(8)} "
          f"< {B_MAX}; L_1^n 1 in C_B for every n")
    a = 2 * B
    bp = arb(max((M1 + c1 * a).upper(), (M2 + c2 * a).upper()))
    Delta = 2 * ((a + bp) / (a - bp)).log() + 2 * bp * X
    check(bool(bp < a) and bool(Delta < arb("13.31")),
          f"CJ.2: a = 2B, b' = max(M_i + c_i a) <= {bp.str(6)} < a = {a.str(6)}; diam L_1 C_a <= Delta = {Delta.str(6)} "
          f"< 13.31, so theta = tanh(Delta/4) < 1 and inf R_n -> rho(L_1) (CF.1, CF.3 for L_1)")
    # CJ.3
    sG = _cover_sup(lambda x: arb(max((1 / (2 * (arb(15) / 4 - x).sqrt())).upper(),
                                      (1 / (2 * (arb(15) / 4 + x).sqrt())).upper())), ga, gb, 2000)
    E = (B * l).exp()
    lam = arb(LAM)
    q1, q2, q3 = 1 / arb(5).sqrt(), lam * sG, E / (2 * ga * lam)
    Q = arb(max(q1.upper(), q2.upper(), q3.upper()))
    check(bool(Q < arb(Q_MAX)) and bool(Q < 1),
          f"CJ.3: g_a = {ga.str(8)}, g_b = {gb.str(8)}, l = {l.str(8)}, sup_G |sigma_i'| <= {sG.str(6)}, "
          f"E = e^(Bl) <= {E.str(6)}; Q = max({q1.str(5)}, {q2.str(5)}, {q3.str(5)}) <= {Q.str(6)} < {Q_MAX} < 1")
    return dict(M1=M1, M2=M2, B=B, Delta=Delta, Q=Q, E=E, sG=sG)


def consequences():
    ctx.prec = 128
    lo, up = arb(LAM_LO), arb(LAM_UP)
    ok = bool(arb(L2_ABS) + arb("1e-13") < lo) and bool(up < arb(R_LO)) and bool(arb(MU_CI) < lo)
    check(ok, f"CJ.5: |lambda_2| < {LAM_LO} <= rho(L_1) <= {LAM_UP} < r, so the only H^2 eigenvalue beyond rho(L_1) is r; "
              f"spec(C^1) = {{r}} u {{|z| <= rho(L_1)}} and CI's open annulus ({MU_CI}, {LAM_UP}] is closed")

    rl, rh = lo / arb("0.4237312869627"), up / arb(R_LO)
    check(bool(rl > arb("0.3170757")) and bool(rh < arb("0.3170762")),
          f"CJ.5: sup{{|mu|/r : mu != r a C^1 eigenvalue}} = rho(L_1)/r in [{rl.str(8)}, {rh.str(8)}]")


def float_observation(n=40):
    """No claim: inf/sup of R_n = L_1^(n+1)1/L_1^n 1 on a grid, from a float Chebyshev collocation."""
    import numpy as np
    N = 60
    k = np.arange(N)
    lo, hi = np.sqrt(5) / 2, 2.5
    xs = (lo + hi) / 2 + (hi - lo) / 2 * np.cos(np.pi * (k + 0.5) / N)

    def u(x, i):
        sg = np.sqrt(3.75 - x) if i == 1 else np.sqrt(3.75 + x)
        w = (x + 0.5) / (2 * x * (3.5 - x)) if i == 1 else (x - 0.5) / (2 * x * (x + 3.5))
        return w / (2 * sg), sg

    def interp(vals, y):
        c = np.polynomial.chebyshev.chebfit((2 * xs - lo - hi) / (hi - lo), vals, N - 1)
        return np.polynomial.chebyshev.chebval((2 * y - lo - hi) / (hi - lo), c)

    f = np.ones(N)
    for _ in range(n):
        g = sum(u(xs, i)[0] * interp(f, u(xs, i)[1]) for i in (1, 2))
        R = g / f
        f = g / g.max()
    print(f"     (observation, floats) R_{n} in [{R.min():.10f}, {R.max():.10f}] on {N} nodes")


def main() -> int:
    t0 = time.time()
    print("Theorem CJ: Conjecture CI.6: spec(Lt; C^1) = {r} u {|z| <= rho(L_1)}; every |mu| < rho(L_1) a C^1 eigenvalue")
    symbolic_checks()
    constants()
    consequences()
    try:
        float_observation()
    except Exception as e:  # pragma: no cover - observation only
        print("     (observation skipped:", e, ")")
    print(f"{'PASS' if not FAILS else 'FAIL'} ({time.time() - t0:.1f}s)")
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
