#!/usr/bin/env python3
r"""Theorem CH: in which function space is CG's spectrum the spectrum?  (CG's kept limit: holomorphic eigenfunctions only.)
On C^k(X), k >= 2, it is: outside an explicit disk the eigenvalues of CE's transfer operator on C^k(X) are exactly CG's,
with the same multiplicities, so l_2 is the second eigenvalue on C^2(X) and l_3 the third on C^3(X).  On C^1(X) it is
not: every mu with |mu| <= 0.1089 is a C^1 eigenvalue (Kept Failure CH.4, explicit series witness), so l_2 (|l_2| =
0.0731) is not the second C^1 eigenvalue; but every C^1 eigenvalue mu != r has |mu| < 0.1345 = 0.3175 r (CF.5: 0.8855 r).

Setting (GASKET.md, Theorems CE-CG): in s = sqrt(t + 25/4) (a real-analytic diffeomorphism [-5,0] -> X, so C^k is the
same space in t and in s)
    (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s),   sigma_(1,2)(s) = sqrt(15/4 -+ s),   X = [sqrt5/2, 5/2],
    wt1 = (s+1/2)/(2s(7/2-s)),  wt2 = (s-1/2)/(2s(s+7/2)).
sigma_1 X = [sqrt5/2, g_a], sigma_2 X = [g_b, 5/2], g_a = sqrt(15/4 - sqrt5/2), g_b = sqrt(15/4 + sqrt5/2); the gap
G = (g_a, g_b) is not seen by Lt.  L_k is the positive operator with weights wt_i |sigma_i'|^k (L_0 = Lt).

Theorem CH (checked here; exit 1 on any failure):
  CH.1 [proved; exact sympy] |sigma_i'| sigma_i = 1/2 and wt1 + wt2 = 16/(49 - 4 s^2), so h(s) = s gives
       L_1 h = 8/(49 - 4 s^2) exactly and c(s) = (49 - 4 s^2)/16 satisfies wt1 c + wt2 c = 1.  On X:
       3/2 <= c <= 11/4, c' = -s/2, and d(s) = s/(L_1 h)(s) = s(49 - 4 s^2)/8 <= 343 sqrt3/72 = 8.2513 (at s = 7 sqrt3/6).
  CH.2 [certified; arb, Collatz-Wielandt] Polynomials h_k > 0 on X with L_k h_k <= lam_k h_k on X (adaptive arb
       interval cover), lam = (0.1345, 0.0432, 0.01396, 0.00454) for k = 1..4.  Hence rho(L_k on C(X)) <= lam_k
       (floats: 0.1343550, 0.0430703, 0.0139193, 0.0045265).
  CH.3 [proved from CH.2, CG, Assumption CH] For k >= 1 and |z| > lam_k, the spectrum of Lt on C^k(X) consists of
       isolated eigenvalues of finite multiplicity, and they are exactly the eigenvalues of Lt on H^2(D) (CG's disk)
       of modulus > lam_k, with equal algebraic multiplicities.  With CG.3:
         C^1(X): only r has |mu| > 0.1345;  C^2(X): only r, l_2 have |mu| > 0.0432;
         C^3(X) and C^4(X): only r, l_2, l_3 have |mu| > 0.01396 (resp. 0.00454); all simple.
       So CF.5's bound for C^1 eigenfunctions improves from 0.8855 r to 0.1345 < 0.3175 r, and on C^2(X) l_2 is the
       second eigenvalue (CG's kept limit removed for k >= 2).
  CH.4 [Kept Failure; proved, arb] On C^1(X), every mu with |mu| <= 0.1089 is an eigenvalue.  Witness:
       g = ((s - g_a)(g_b - s))^2 on G, 0 elsewhere (C^1, Lt g = 0); K f = c(T s) f(T s) on sigma_i X (T = sigma_i^-1),
       cubic Hermite interpolation on G; then Lt K = I on C^1(X), and with N(f) = (sup|f|, sup|f'|/s),
       N(K f) <= A N(f) entrywise for an explicit 2x2 matrix A with rho(A) < 9.18.  So v = sum_n mu^n K^n g converges
       in C^1(X), Lt v = mu v, and v != 0 (on G it is a quartic plus a cubic).  In particular mu = +-1/10, +-1/12 are C^1
       eigenvalues with |mu| > |l_2| = 0.0731: CG's l_2 is NOT the second eigenvalue of Lt on C^1(X), and
           0.1089 <= sup{|mu| : mu != r a C^1 eigenvalue} <= 0.1345   (ratio to r in [0.2569, 0.3175]).
Assumption CH (literature; used as stated, not re-proved): Hennion's theorem (Proc. AMS 118 (1993)): if B is compactly
embedded in B_w and ||P^n f||_B <= A_n ||f||_B + B_n |f|_(B_w), then the essential spectral radius of P on B is
<= liminf A_n^(1/n); outside it the spectrum is isolated eigenvalues of finite algebraic multiplicity (Riesz projections).
Uses CG.1-CG.3 as stated (imported constants only; CG's certificate is not re-run here).  Level sums are untouched:
CG.4 concerns the analytic Lt^4 G, so its rate 0.17255 stands.  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW enters.
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

import sympy as sp
from flint import arb, ctx

FAILS: list[str] = []
LAM = {1: "0.1345", 2: "0.0432", 3: "0.01396", 4: "0.00454"}     # certified upper bounds for rho(L_k)
DEG_H = 10                # degree of the Collatz-Wielandt test polynomials h_k
RHO_A_MAX = "9.18"        # certified bound for the spectral radius of the 2x2 majorant of K
MU_MAX = "0.1089"         # every |mu| <= MU_MAX is a C^1 eigenvalue
L2_ABS = "0.0731138175124"    # |l_2| from CG (to 1e-13)
L3_ABS = "0.0170676440626"    # |l_3| from CG (to 1e-12)
R_LO = "0.4237312869626"      # r from CE.4 (lower, to 1e-13)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- CH.1 exact identities
def symbolic_checks():
    s, t = sp.symbols("s t", positive=True)
    h = sp.Rational(1, 2)
    W1 = (s + h) / (2 * s * (sp.Rational(7, 2) - s))
    W2 = (s - h) / (2 * s * (s + sp.Rational(7, 2)))
    S1 = sp.sqrt(sp.Rational(15, 4) - s)
    S2 = sp.sqrt(sp.Rational(15, 4) + s)
    d1, d2 = -sp.diff(S1, s), sp.diff(S2, s)
    check(sp.simplify(d1 * S1 - h) == 0 and sp.simplify(d2 * S2 - h) == 0,
          "CH.1: |sigma_i'| sigma_i = 1/2 for both branches (sigma_1' < 0 < sigma_2')")
    check(sp.simplify(W1 + W2 - 16 / (49 - 4 * s ** 2)) == 0, "CH.1: wt1 + wt2 = 16/(49 - 4 s^2) (= 4/(6 - t), CE.2)")
    L1h = sp.simplify(W1 * d1 * S1 + W2 * d2 * S2)
    check(sp.simplify(L1h - 8 / (49 - 4 * s ** 2)) == 0, "CH.1: L_1 h = 8/(49 - 4 s^2) exactly for h(s) = s")
    c = (49 - 4 * s ** 2) / 16
    check(sp.simplify(W1 * c + W2 * c - 1) == 0 and sp.diff(c, s) == -s / 2,
          "CH.1: c = (49 - 4 s^2)/16 gives wt1 c + wt2 c = 1 (so Lt K = I), c' = -s/2")
    a, b = sp.sqrt(5) / 2, sp.Rational(5, 2)
    check(c.subs(s, b) == sp.Rational(3, 2) and sp.simplify(c.subs(s, a) - sp.Rational(11, 4)) == 0,
          "CH.1: c decreasing on X from c(sqrt5/2) = 11/4 to c(5/2) = 3/2")
    d = s * (49 - 4 * s ** 2) / 8
    check(sp.simplify(d - s / (8 / (49 - 4 * s ** 2))) == 0, "CH.1: d = s/(L_1 h) = s(49 - 4 s^2)/8")
    crit = sp.solve(sp.diff(d, s), s)
    crit = [x for x in crit if x.is_positive]
    dmax = sp.nsimplify(sp.simplify(d.subs(s, crit[0])))
    check(len(crit) == 1 and sp.simplify(crit[0] - 7 * sp.sqrt(3) / 6) == 0 and bool(a < crit[0] < b)
          and sp.simplify(dmax - 343 * sp.sqrt(3) / 72) == 0
          and bool(d.subs(s, a) < dmax) and bool(d.subs(s, b) < dmax),
          f"CH.1: sup_X d = 343 sqrt3/72 = {float(343 * sp.sqrt(3) / 72):.6f} at s = 7 sqrt3/6 in X")
    # branch images: sigma_1 X = [sqrt5/2, g_a], sigma_2 X = [g_b, 5/2]
    check(sp.simplify(S1.subs(s, b) - a) == 0 and sp.simplify(S2.subs(s, b) - b) == 0
          and bool(S1.subs(s, a) < S2.subs(s, a)),
          "CH.1: sigma_1(5/2) = sqrt5/2, sigma_2(5/2) = 5/2, g_a = sigma_1(sqrt5/2) < g_b = sigma_2(sqrt5/2): X = sigma_1X u G u sigma_2X")
    # Hermite basis on [0,1]
    t = sp.symbols("t", real=True)
    h00, h01 = 2 * t ** 3 - 3 * t ** 2 + 1, -2 * t ** 3 + 3 * t ** 2
    h10, h11 = t ** 3 - 2 * t ** 2 + t, t ** 3 - t ** 2
    ok = sp.expand(h00 + h01 - 1) == 0 and sp.expand(h00 - (t - 1) ** 2 * (2 * t + 1)) == 0 and sp.expand(h01 - t ** 2 * (3 - 2 * t)) == 0
    ok = ok and sp.expand(h10 - t * (1 - t) ** 2) == 0 and sp.expand(h11 + t ** 2 * (1 - t)) == 0
    # |h10| + |h11| = t(1-t) <= 1/4 ; |h10'| + |h11'| <= 1 piecewise
    p1, p2 = sp.diff(h10, t), sp.diff(h11, t)
    ok = ok and sp.expand(p1 - (t - 1) * (3 * t - 1)) == 0 and sp.expand(p2 - t * (3 * t - 2)) == 0
    pieces = [((0, sp.Rational(1, 3)), p1 - p2), ((sp.Rational(1, 3), sp.Rational(2, 3)), -p1 - p2),
              ((sp.Rational(2, 3), 1), -p1 + p2)]
    for (lo, hi), e in pieces:
        mx = max([e.subs(t, lo), e.subs(t, hi)] + [e.subs(t, x) for x in sp.solve(sp.diff(e, t), t) if lo <= x <= hi])
        ok = ok and mx <= 1
    check(ok, "CH.4: Hermite basis: h00, h01 >= 0 with sum 1; |h10| + |h11| = t(1-t) <= 1/4; |h10'| + |h11'| <= 1 on [0,1]")
    # the bump g
    ga, gb = sp.symbols("g_a g_b", positive=True)
    g = ((s - ga) * (gb - s)) ** 2
    check(all(sp.simplify(e) == 0 for e in (g.subs(s, ga), g.subs(s, gb), sp.diff(g, s).subs(s, ga), sp.diff(g, s).subs(s, gb)))
          and sp.Poly(g, s).degree() == 4 and sp.Poly(g, s).LC() == 1,
          "CH.4: g = ((s-g_a)(g_b-s))^2 on G, 0 off G, is C^1 on X (g = g' = 0 at g_a, g_b) and a monic quartic on G")


# ----------------------------------------------------------------- CH.2 Collatz-Wielandt bounds for rho(L_k)
def _float_eigfun(k, n=4001, it=200):
    import numpy as np
    a, b = 5 ** 0.5 / 2, 2.5
    xs = np.linspace(a, b, n)
    w1 = (xs + .5) / (2 * xs * (3.5 - xs))
    w2 = (xs - .5) / (2 * xs * (xs + 3.5))
    y1, y2 = np.sqrt(3.75 - xs), np.sqrt(3.75 + xs)
    k1, k2 = w1 / (2 * y1) ** k, w2 / (2 * y2) ** k
    f = np.ones(n)
    lam = 0.0
    for _ in range(it):
        g = k1 * np.interp(y1, xs, f) + k2 * np.interp(y2, xs, f)
        lam = g.max() / f.max()
        f = g / g.max()
    u = (xs - 1.7) / 0.7
    coef = np.polynomial.polynomial.polyfit(u, f, DEG_H)
    return lam, [Fr(float(c)).limit_denominator(10 ** 12) for c in coef]


def _h_eval(coef, x):
    u = (x - arb("1.7")) / arb("0.7")
    acc = arb(0)
    for c in reversed(coef):
        acc = acc * u + arb(c.numerator) / c.denominator
    return acc


def _Lk_minus(coef, k, lam, x):
    h = arb(1) / 2
    y1 = (arb(15) / 4 - x).sqrt()
    y2 = (arb(15) / 4 + x).sqrt()
    w1 = (x + h) / (2 * x * (arb(7) / 2 - x))
    w2 = (x - h) / (2 * x * (x + arb(7) / 2))
    Lh = w1 * _h_eval(coef, y1) / (2 * y1) ** k + w2 * _h_eval(coef, y2) / (2 * y2) ** k
    return Lh - lam * _h_eval(coef, x), _h_eval(coef, x)


def collatz_wielandt(k, maxdepth=22):
    ctx.prec = 128
    lam_f, coef = _float_eigfun(k)
    lam = arb(LAM[k])
    lo0, hi0 = arb(5).sqrt() / 2, arb(5) / 2
    stack = [(lo0, hi0, 0)]
    n_int = 0
    ok = True
    while stack:
        lo, hi, dep = stack.pop()
        mid = (lo + hi) / 2                 # ball [mid - R, mid + R] covering [lo, hi] (lo, hi themselves balls)
        x = arb(mid.mid(), ((hi - lo) / 2).upper() + mid.rad() + lo.rad() + hi.rad())
        e, hv = _Lk_minus(coef, k, lam, x)
        if bool(e < 0) and bool(hv > 0):
            n_int += 1
            continue
        if dep >= maxdepth:
            ok = False
            break
        m = arb(((lo + hi) / 2).mid())
        stack.append((lo, m, dep + 1))
        stack.append((m, hi, dep + 1))
    check(ok, f"CH.2: k = {k}: h_{k} > 0 and L_{k} h_{k} < {LAM[k]} h_{k} on X ({n_int} arb intervals), so rho(L_{k}) <= {LAM[k]} "
              f"(float {lam_f:.7f})")
    return lam_f


# ----------------------------------------------------------------- CH.4 the C^1 majorant of K
def k_majorant():
    ctx.prec = 128
    s0 = arb(5).sqrt() / 2
    ga = (arb(15) / 4 - s0).sqrt()
    gb = (arb(15) / 4 + s0).sqrt()
    hG = gb - ga
    c0 = arb(11) / 4                       # c(s0)
    dY = arb(343) * arb(3).sqrt() / 72     # sup_X s/(L_1 h)
    mb = gb * s0                           # |m_j| <= g_j s0 (|f| + (11/2)|f|_h), max at g_b
    A00 = c0 + hG / 4 * mb                 # Y: c <= 11/4; G: |H| <= c0|f(s0)| + (hG/4) max|m_j|
    A01 = hG / 4 * mb * arb(11) / 2
    A10 = arb(5) / 2                       # Y: sup s = 5/2 ; G: max|m_j|/g_a = mb/ga < 5/2 (checked)
    A11 = mb * arb(11) / 2 / ga            # G term; Y term dY (checked smaller)
    ok = bool(mb / ga < A10) and bool(dY < A11)
    tr, det = A00 + A11, A00 * A11 - A01 * A10
    rho = (tr + (tr * tr - 4 * det).sqrt()) / 2
    check(ok and bool(tr * tr - 4 * det > 0), f"CH.4: majorant A = [[{A00.str(6)}, {A01.str(6)}], [{A10.str(6)}, {A11.str(6)}]] "
          f"(Y-entries 11/4, 5/2, {dY.str(6)} dominated where stated); g_a = {ga.str(8)}, g_b = {gb.str(8)}")
    check(bool(rho < arb(RHO_A_MAX)) and bool(arb(MU_MAX) * rho < 1),
          f"CH.4: rho(A) = {rho.str(8)} < {RHO_A_MAX}, so sum mu^n K^n g converges in C^1(X) for |mu| <= {MU_MAX} "
          f"(|mu| rho(A) <= {(arb(MU_MAX) * rho).str(6)} < 1)")
    l2 = arb(L2_ABS) + arb("1e-13")
    check(bool(arb(MU_MAX) > l2) and bool(arb(1) / 10 > l2) and bool(arb(1) / 12 > l2),
          f"CH.4: |l_2| < 1/12 < 1/10 < {MU_MAX}: C^1 eigenvalues of modulus > |l_2| exist (+-1/10, +-1/12); l_2 is not "
          "the second eigenvalue of Lt on C^1(X)")
    return rho


def consequences(lams):
    ctx.prec = 128
    r_lo = arb(R_LO)
    l2, l3 = arb(L2_ABS), arb(L3_ABS)
    ok = bool(arb(LAM[1]) < r_lo) and bool(arb(LAM[1]) > l2 + arb("1e-13"))
    ok = ok and bool(arb(LAM[2]) < l2 - arb("1e-13")) and bool(arb(LAM[2]) > l3 + arb("1e-12"))
    ok = ok and bool(arb(LAM[3]) < l3 - arb("1e-12")) and bool(arb(LAM[4]) < l3 - arb("1e-12")) and bool(arb(LAM[4]) > arb("1e-3"))
    check(ok, "CH.3: with CG.3 (r, l_2, l_3 simple; rest < 1e-3 on H^2): C^1: only r beyond 0.1345; C^2: r, l_2 beyond 0.0432; "
              "C^3, C^4: r, l_2, l_3 beyond 0.01396, 0.00454")
    ratio_hi = arb(LAM[1]) / r_lo
    ratio_lo = arb(MU_MAX) / (arb(R_LO) + arb("1e-13"))
    check(bool(ratio_hi < arb("0.3175")) and bool(ratio_lo > arb("0.2569")),
          f"CH.3/CH.4: sup{{|mu|/r : mu != r a C^1 eigenvalue}} in [{ratio_lo.str(5)}, {ratio_hi.str(5)}] (CF.5: <= 0.8855)")
    check(bool(arb(MU_MAX) <= arb(LAM[1])), "CH.3/CH.4 consistency: the witness disk radius 0.1089 is <= the C^1 bound 0.1345")


def k_growth_observation(nmax=8):
    """Floats, no claim: sup|(K^n g)'|/s grows by about 7.5 per step (1/rho(L_1) = 7.44; certified majorant 9.18)."""
    import numpy as np
    s0 = 5 ** 0.5 / 2
    ga, gb = (3.75 - s0) ** 0.5, (3.75 + s0) ** 0.5
    hG = gb - ga
    xs = np.linspace(s0, 2.5, 200001)

    def c(s):
        return (49 - 4 * s * s) / 16

    def K(f, df):
        F, D = np.empty_like(xs), np.empty_like(xs)
        m1, m2 = xs <= ga, xs >= gb
        mg = ~(m1 | m2)
        for m, sgn in ((m1, -1), (m2, 1)):
            sv = sgn * (xs[m] ** 2 - 3.75)
            fs, dfs = np.interp(sv, xs, f), np.interp(sv, xs, df)
            F[m] = c(sv) * fs
            D[m] = sgn * 2 * xs[m] * (-sv / 2 * fs + c(sv) * dfs)
        y0 = c(s0) * f[0]
        q = -s0 / 2 * f[0] + c(s0) * df[0]
        m0, m1_ = -2 * ga * q, 2 * gb * q
        tt = (xs[mg] - ga) / hG
        F[mg] = y0 + hG * (m0 * (tt ** 3 - 2 * tt ** 2 + tt) + m1_ * (tt ** 3 - tt ** 2))
        D[mg] = m0 * (3 * tt ** 2 - 4 * tt + 1) + m1_ * (3 * tt ** 2 - 2 * tt)
        return F, D

    inG = (xs > ga) & (xs < gb)
    f = np.where(inG, ((xs - ga) * (gb - xs)) ** 2, 0.0)
    df = np.where(inG, 2 * (xs - ga) * (gb - xs) * (gb + ga - 2 * xs), 0.0)
    n0 = np.abs(df / xs).max()
    for _ in range(nmax):
        f, df = K(f, df)
    rate = (np.abs(df / xs).max() / n0) ** (1 / nmax)
    print(f"   observation (no claim): (sup|(K^n g)'|/s / sup|g'|/s)^(1/n) = {rate:.3f} at n = {nmax}; 1/rho(L_1) = 7.443")
    return rate


def main() -> int:
    t0 = time.time()
    print("Theorem CH: CG's spectrum is the C^k spectrum for k >= 2 beyond explicit radii, not on C^1 (Kept Failure CH.4: "
          "every |mu| <= 0.1089 is a C^1 eigenvalue). Not thrust.")
    symbolic_checks()
    lams = {k: collatz_wielandt(k) for k in (1, 2, 3, 4)}
    k_majorant()
    consequences(lams)
    k_growth_observation()
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CH checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
