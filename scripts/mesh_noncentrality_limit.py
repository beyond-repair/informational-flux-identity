#!/usr/bin/env python3
"""Theorem BD: the non-centrality ratio E_n/Y_n of Theorem BC tends to 2 (Conjecture BC.2 proved).

Setting: Assumption BA / BC (README), notation of Theorems BB, BC.  For s > 0, level-n gasket:
  y_n(s) = [(L_n+s)^-1]_{p'c'} - [(L_n+s)^-1]_{p'c''},  e_n(s) = s g_n(s),  t_n(s) = s [(L_n+s)^-1 1_C]_c,
  Y_n = (sin pi a/pi) int s^-a y_n ds,  E_n = (sin pi a/pi) int s^-a e_n ds,  p(s) = s(s+5).

Theorem BD (checked here; exit 1 on any failure):
  (1) exact identity, every n, every s > 0:  e_n(s) = (2 - t_n(s)) y_n(s).
      Hence 2 - E_n/Y_n = <t_n>_n, the mean of t_n under the probability measure prop. to s^-a y_n(s) ds.
  (2) t_n(s) = (1/3) sum_lambda ||P_lambda 1_C||^2 s/(lambda+s) (S3 symmetry), so t_n is increasing and concave,
      t_n(0+) = a_n = 2/(3^n+1), t_n'(0) = b_n = 3^n(5^n-1)/(6(3^n+1)^2) <= (1/6)(5/3)^n, t_n < 1.
      Consequences for every n >= 1, a in (0,1):  1 < E_n/Y_n < 2*3^n/(3^n+1), pointwise e_n/y_n <= 2*3^n/(3^n+1);
      kappa_n^cross = 2^n (4/R_n^2) Y_n (1 + <t_n>) > 0 at EVERY level (Kept Failure BB.1 for all n, all a);
      tan theta_n = -1/(3 sqrt3) - 4 d_n/(3 sqrt3 (3 - d_n)),  d_n = 2 - E_n/Y_n, approached from below.
  (3) two-sided rate, n >= 2:  2/(3^n+1) < 2 - E_n/Y_n <= 2/(3^n+1) + B_a 3^-n,
      B_a = (1-a) J_a e^0.1 (3+e^0.1)/18,  J_a = int_0^inf tau^(1-a) dtau/((1+tau/25)(1+tau/125)) (closed form).
      So E_n/Y_n -> 2 and tan theta_n -> -1/(3 sqrt3) for every a in (0,1): Conjecture BC.2 holds.
  (4) the 3^-n constant exists: T_n(tau) = 3^n t_n(5^-n tau) increases to T(tau) <= 2 + tau/6, and
      Y_n(tau) = 3^n y_n(5^-n tau) -> Y(tau) > 0 with explicit summable log-increments; dominated convergence gives
          c(a) = lim 3^n (2 - E_n/Y_n) = int tau^-a T Y / int tau^-a Y,   2 < c(a) <= 2 + B_a.
      Certified (mpmath.iv, outward Riemann bounds + explicit tails): c(0.45) in [3.72198, 3.72915] (contains the
      float 3.7255634), c(0.8) in [2.57126, 2.57353], so c depends on a (Kept Failure BD.1: c != 2 = lim 3^n t_n(0+)).

Not a field solver.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_maxwell_limit as bb  # noqa: E402
import mesh_centred_stencil as bc  # noqa: E402

FAILS: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    s, t, x, X, Z, y, r = sp.symbols("s t x X Z y r", positive=True)
    ps = s * (s + 5)
    # induction step of (1): Q_{n+1}(s) = (6+s)/(6+s-2t_n(ps)) Q_n(ps) (BC.5) and t_{n+1}(s) = (2+s)t/(6+s-2t) (BB)
    lhs = (6 + s) / (6 + s - 2 * t) * (2 - t)
    rhs = 2 - (2 + s) * t / (6 + s - 2 * t)
    check(sp.simplify(lhs - rhs) == 0, "induction step: (6+s)(2-t)/(6+s-2t) = 2 - (2+s)t/(6+s-2t); base e_0/y_0 = 1 = 2 - t_0")
    # y recursion in the form used for the scaling limit
    c = (6 + s) / ((5 + s) * (2 + s))
    e = 2 * s / (2 + s)
    check(sp.simplify(1 / ((5 + s) * (c - e * r / ps)) - (2 + s) / (6 + s - 2 * r)) == 0,
          "y_{n+1}(s) = (2+s) y_n(ps)/(6+s-2 p(s) r_n(ps))  (BB.2 with AQ's den)")
    # derivative recursion for t'(0): differentiate t_{n+1}(s) = (2+s) T(ps)/(6+s-2T(ps))
    a, b = sp.symbols("a b", positive=True)
    T = a + b * ps
    tn1 = (2 + s) * T / (6 + s - 2 * T)
    d0 = sp.simplify(sp.diff(tn1, s).subs(s, 0))
    check(sp.simplify(d0 - (60 * b + 4 * a - 2 * a**2) / (6 - 2 * a) ** 2) == 0,
          "t_{n+1}'(0) = (60 b_n + 4 a_n - 2 a_n^2)/(6 - 2 a_n)^2,  a_{n+1} = 2a_n/(6-2a_n)")
    # closed form b_n = X(Z-1)/(6(X+1)^2), X = 3^n, Z = 5^n, a_n = 2/(X+1)
    an = 2 / (X + 1)
    bn = X * (Z - 1) / (6 * (X + 1) ** 2)
    bn1 = (3 * X) * (5 * Z - 1) / (6 * (3 * X + 1) ** 2)
    check(sp.simplify((60 * bn + 4 * an - 2 * an**2) / (6 - 2 * an) ** 2 - bn1) == 0
          and sp.simplify(2 * an / (6 - 2 * an) - 2 / (3 * X + 1)) == 0,
          "closed form b_n = 3^n(5^n-1)/(6(3^n+1)^2) satisfies the recursion (b_0 = 0); so b_n <= (1/6)(5/3)^n")
    # tan theta in terms of d = 2 - E/Y
    d = sp.Symbol("d")
    f = lambda q: (q - 3) / (sp.sqrt(3) * (q + 1))
    check(sp.simplify(f(2 - d) - (-1 / (3 * sp.sqrt(3)) - 4 * d / (3 * sp.sqrt(3) * (3 - d)))) == 0,
          "tan theta_n = -1/(3 sqrt3) - 4 d_n/(3 sqrt3 (3 - d_n)),  d_n = 2 - E_n/Y_n")
    # 3Y - E = Y(1 + d): cross part sign
    check(sp.simplify(3 * y - (2 - d) * y - y * (1 + d)) == 0, "3Y_n - E_n = Y_n (1 + d_n) > 0")
    # J_a closed form: 1/((1+t/25)(1+t/125)) = (3125/100)(1/(25+t) - 1/(125+t))
    tt = sp.Symbol("tt", positive=True)
    check(sp.simplify(1 / ((1 + tt / 25) * (1 + tt / 125)) - sp.Rational(3125, 100) * (1 / (25 + tt) - 1 / (125 + tt))) == 0,
          "partial fractions for J_a: (125/4)(1/(25+t) - 1/(125+t))")


def J_closed(al, ctx=mp):
    """J_a = (125/4) (pi / sin(pi(2-a))) (25^(1-a) - 125^(1-a)), by continuation of int t^(mu-1)/(c+t) = c^(mu-1) pi/sin(pi mu)."""
    mu = 2 - al
    return ctx.mpf(125) / 4 * ctx.pi / ctx.sin(ctx.pi * mu) * (ctx.mpf(25) ** (mu - 1) - ctx.mpf(125) ** (mu - 1))


def B_const(al, ctx=mp):
    e1 = ctx.exp(ctx.mpf(1) / 10)
    return (1 - al) * J_closed(al, ctx) * e1 * (3 + e1) / 18


# ----------------------------------------------------------------- (B) exact rationals
def exact_checks():
    ok = True
    for n in range(0, 9):
        for s in (Fraction(1, 7), Fraction(3), Fraction(1, 1000), Fraction(22, 5)):
            r, t, g, y = bb.rec_all(n, s)
            ok = ok and s * g == (2 - t) * y and 0 < t <= 1 and (n == 0 or t < 1)
    check(ok, "exact rationals, n = 0..8, four s: e_n(s) = (2 - t_n(s)) y_n(s); 0 < t_n < 1 for n >= 1 (t_0 = 1)")
    ok = True
    for n in (1, 2, 3):
        P, adj, C, pb, pc = bb.corner_frame(n)
        a, b_, c_ = C
        for s in (Fraction(1, 2), Fraction(7, 5)):
            Wb, Wc, W1 = bb.exact_resolvent_cols(adj, s, [{b_}, {c_}, set(C)])
            y = Wb[pb] - Wc[pb]
            e = Wb[pb] + Wc[pb] - Wb[a] - Wc[a]
            t = s * W1[a]
            ok = ok and e == (2 - t) * y
    check(ok, "exact resolvent solves on the meshes, levels 1-3: e_n = (2 - t_n) y_n (independent of the recursions)")
    a, b = Fraction(1), Fraction(0)
    ok = True
    for n in range(30):
        ok = ok and a == Fraction(2, 3**n + 1) and b == Fraction(3**n * (5**n - 1), 6 * (3**n + 1) ** 2)
        a, b = 2 * a / (6 - 2 * a), (60 * b + 4 * a - 2 * a * a) / (6 - 2 * a) ** 2
    check(ok, "a_n = 2/(3^n+1), b_n = 3^n(5^n-1)/(6(3^n+1)^2), n = 0..29 (exact)")


# ----------------------------------------------------------------- (C) dense floats on the meshes
def dense_checks(levels=(1, 2, 3, 4)):
    ok_w = ok_b = True
    for n in levels:
        P, adj, C, pb, pc = bb.corner_frame(n)
        L = bb.laplacian(adj)
        w, V = np.linalg.eigh(L)
        one = np.zeros(len(adj))
        one[C] = 1
        coef = V.T @ one
        # group eigenvalues; weight of eigenspace at corner c: [P_lam 1_C]_c = (1/3)||P_lam 1_C||^2 >= 0
        lam = np.round(w, 9)
        for lv in np.unique(lam):
            idx = lam == lv
            Pl = V[:, idx] @ coef[idx]
            ok_w = ok_w and all(abs(Pl[c] - (Pl @ Pl) / 3) < 1e-10 for c in C) and Pl[C[0]] > -1e-12
        wv = one - 3 / len(adj)
        Lp = np.linalg.pinv(L)
        bn = 3**n * (5**n - 1) / (6 * (3**n + 1) ** 2)
        ok_b = ok_b and abs((wv @ Lp @ wv) / 3 - bn) < 1e-9 * bn
    check(ok_w, f"levels {levels[0]}-{levels[-1]}: [P_lam 1_C]_c = ||P_lam 1_C||^2/3 >= 0 for every eigenspace (t_n Stieltjes, increasing, concave)")
    check(ok_b, "dense: t_n'(0) = <w, L^+ w>/3 = 3^n(5^n-1)/(6(3^n+1)^2), w = 1_C - 3/N")


def quad_mean_t(n, al="0.45", dps=30):
    """<t_n> = int s^-a t_n y_n / int s^-a y_n by quadrature of the recursions (consistency, not interval)."""
    mp.dps = dps
    a = mpf(al)
    pts = [0] + [mpf(10) ** k for k in range(-30, 31, 2)]
    num = mp.quad(lambda s: s ** (-a) * rec(n, s, "ty"), pts)
    den = mp.quad(lambda s: s ** (-a) * rec(n, s, "y"), pts)
    return num / den


def rec(n, s, what):
    r, t, g, y = bb.rec_all(n, s)
    return t * y if what == "ty" else y


# ----------------------------------------------------------------- (D) interval tools for the scaling limit
def p(x):
    return x * (x + 5)


def TY(N, tau):
    """T_N(tau) = 3^N t_N(5^-N tau), Y_N(tau) = 3^N y_N(5^-N tau), interval arithmetic along the orbit."""
    s = tau / iv.mpf(5) ** N
    o = [s]
    for _ in range(N):
        o.append(p(o[-1]))
    x = o[N]
    r, t, yr = 1 / (3 + x), iv.mpf(1), iv.mpf(1)
    for k in range(N - 1, -1, -1):
        z = o[k]
        r = (2 + z) * (5 + z) * r / ((6 + z) - 2 * z * (5 + z) * r)
        t = (2 + z) * t / (6 + z - 2 * t)
        yr = yr / (5 + z)
    return iv.mpf(3) ** N * t, iv.mpf(3) ** N * yr * r


def psibar(tau):
    """sup_n psi_n(tau), psi_n(tau) = p^n(5^-n tau): psi_m(u) <= u e^(u/10) for u <= 1, psi_n = p^k(psi_{n-k}(5^-k tau))."""
    k = 0
    while tau > 5**k:
        k += 1
    u = iv.mpf(tau) / iv.mpf(5) ** k
    v = u * iv.exp(u / 10)
    for _ in range(k):
        v = p(v)
    return v


def epsT(N, tau):
    """sum_{n>=N} sup log(T_{n+1}/T_n) <= sum [delta_n + s_n/2 + 3^-n (1 + tau(1+delta_n)/12)]."""
    tau = iv.mpf(tau)
    return (tau * iv.mpf(5) ** (-N) * (iv.mpf(1) / 25 + iv.mpf(1) / 10) * iv.mpf(5) / 4
            + iv.mpf(3) / 2 * iv.mpf(3) ** (-N) * (1 + tau / 12) + tau**2 / 300 * iv.mpf(15) ** (-N) * iv.mpf(15) / 14)


def epsY(N, tau):
    """sum_{n>=N} |log Y_{n+1} - log Y_n| <= sum [s_n/2 + tau 3^-n(1+delta_n)/6 + delta_n + psibar eta_n e^eta_n/45]."""
    tau = iv.mpf(tau)
    S5 = iv.mpf(5) ** (-N) * iv.mpf(5) / 4
    S3 = iv.mpf(3) ** (-N) * iv.mpf(3) / 2
    S15 = iv.mpf(15) ** (-N) * iv.mpf(15) / 14
    etaN = iv.mpf(2) ** N * iv.mpf(5) ** (-N - 2) * tau
    Seta = tau / 25 * (iv.mpf(2) / 5) ** N * iv.mpf(5) / 3
    return tau / 5 * S5 / 2 + tau / 6 * (S3 + tau / 25 * S15) + tau / 25 * S5 + psibar(tau) * Seta * iv.exp(etaN) / 45


def c_bracket(al="0.45", N=60, tau0=1e-20, tau1=3125, r1=1.05, r2=1.001, split=0.01, M=6, k=8):
    """Certified enclosure of c(a) = int tau^-a T Y / int tau^-a Y from monotone cell bounds and explicit tails."""
    iv.dps = 20
    a = iv.mpf(al)
    grid = [tau0]
    while grid[-1] < split:
        grid.append(grid[-1] * r1)
    while grid[-1] < tau1:
        grid.append(min(grid[-1] * r2, tau1))
    vals = [TY(N, iv.mpf(g)) for g in grid]
    eT = [iv.exp(epsT(N, g)) for g in grid]
    eY = [iv.exp(epsY(N, g)) for g in grid]
    numL = numU = denL = denU = iv.mpf(0)
    for i in range(len(grid) - 1):
        u, v = iv.mpf(grid[i]), iv.mpf(grid[i + 1])
        w = (v ** (1 - a) - u ** (1 - a)) / (1 - a)
        Tu, Yu = vals[i]
        Tv, Yv = vals[i + 1]
        # T increasing, Y decreasing, eps increasing in tau
        numU += Tv * eT[i + 1] * Yu * eY[i + 1] * w
        numL += Tu * Yv / eY[i + 1] * w
        denU += Yu * eY[i + 1] * w
        denL += Yv / eY[i + 1] * w
    # head [0, tau0]: 0 <= T <= 2 + tau/6, c0 <= Y <= 1/3
    h = iv.mpf(tau0) ** (1 - a) / (1 - a)
    e1 = iv.exp(iv.mpf(1) / 10)
    c0 = 1 / (e1 * (3 + e1))
    numU += (2 + iv.mpf(tau0) / 6) * h / 3
    denU += h / 3
    denL += c0 * h
    # tail [tau1, inf): Y <= tau/(3 psi_M(tau)) and, for tau >= tau1, psi_M(tau) >= (tau/tau1)^k P_k with
    # P_k = sum_{j>=k} c_j tau1^j (psi_M has nonnegative coefficients c_j); T <= 2 + tau/6
    xx = sp.Symbol("x")
    q = xx
    for _ in range(M):
        q = sp.expand(q * (q + 5))
    q = sp.Poly(q.subs(xx, xx / sp.Integer(5) ** M), xx)
    Pk = sp.Rational(q.eval(sp.Integer(tau1))) - sum(q.coeff_monomial(xx**j) * sp.Integer(tau1) ** j for j in range(1, k))
    Pki = iv.mpf(int(Pk.p)) / iv.mpf(int(Pk.q))
    t1 = iv.mpf(tau1)
    pre = t1**k / (3 * Pki)
    numU += pre * (2 * t1 ** (2 - a - k) / (k - 2 + a) + t1 ** (3 - a - k) / (6 * (k - 3 + a)))
    denU += pre * t1 ** (2 - a - k) / (k - 2 + a)
    lo = numL / denU
    hi = numU / denL
    return iv.mpf([lo.a, hi.b]), len(grid)


# ----------------------------------------------------------------- main
def main(cert_levels=(1, 2, 3), r2=1.001) -> int:
    t0 = time.time()
    print("Theorem BD: E_n/Y_n -> 2 on the 0.45 gasket mesh (Conjecture BC.2). Not thrust.")
    symbolic_checks()
    exact_checks()
    dense_checks()

    # B_a, closed form vs quadrature
    mp.dps = 30
    for al in ("0.25", "0.45", "0.8"):
        a = mpf(al)
        Jq = mp.quad(lambda t: t ** (1 - a) / ((1 + t / 25) * (1 + t / 125)), [0, 1, 25, 125, mp.inf])
        check(abs(Jq - J_closed(a)) < 1e-6 * Jq, f"a={al}: J_a closed form = quadrature ({mp.nstr(Jq, 12)})")
    iv.dps = 30
    B45 = B_const(iv.mpf("0.45"), iv)
    print(f"   B_0.45 in {bc.ba.s(B45, 12)}")

    # intervals: certified E_n/Y_n (BC pipeline) inside the bracket of (2)-(3)
    for n in cert_levels:
        rr = bc.certify_level(n)
        ratio = rr["E"] / rr["Y"]
        d = 2 - ratio
        lo = iv.mpf(2) / (iv.mpf(3) ** n + 1)
        ok = d.a > lo.b and ratio.a > 1
        if n >= 2:
            ok = ok and d.b < (lo + B45 / iv.mpf(3) ** n).a
        r3 = iv.sqrt(iv.mpf(3))
        tan = (ratio - 3) / (r3 * (ratio + 1))
        tform = -1 / (3 * r3) - 4 * d / (3 * r3 * (3 - d))
        ok = ok and not bc.ba.disjoint(tan, tform) and tan.b < (-1 / (3 * r3)).a
        check(ok, f"level {n}: certified 2 - E_n/Y_n in {bc.ba.s(d, 16)} lies in (2/(3^n+1), 2/(3^n+1) + B 3^-n]; "
                  f"tan theta_n < -1/(3 sqrt3)")

    # float consistency of (1) after integration: E_n = 2Y_n - <t_n> Y_n
    for n in (1, 2, 4):
        qd = bb.quad_level(n, "0.45")
        mt = quad_mean_t(n)
        dd = bb.dense_level(n)
        check(abs(float(2 - qd["E"] / qd["Y"] - mt)) < 1e-12 and abs(dd["E"] / dd["Y"] - float(2 - mt)) < 1e-8,
              f"level {n}: 2 - E_n/Y_n = <t_n> = {mp.nstr(mt, 14)} (quadrature and dense mesh)")

    # recursion table: 3^n d_n increasing towards c(a), below 2 + B
    mp.dps = 30
    tab = {}
    for al, levels in (("0.45", (5, 10, 20, 30)), ("0.8", (10, 20, 30)), ("0.25", (10, 16, 20))):
        for n in levels:
            q = bb.quad_level(n, al)
            tab[(al, n)] = (2 - q["E"] / q["Y"]) * mpf(3) ** n
            print(f"   a={al} n={n:2d}: 3^n (2 - E_n/Y_n) = {mp.nstr(tab[(al, n)], 12)}")
    check(2 < tab[("0.45", 10)] < tab[("0.45", 20)] < float(B45.a) + 2 and abs(tab[("0.45", 30)] - tab[("0.45", 20)]) < 1e-7,
          "a=0.45: 3^n d_n increases and settles (float 3.7255634) inside (2, 2 + B_0.45]")

    # pointwise scaling limits, monotone T_n (proved), Y_n Cauchy (proved): float display
    iv.dps = 20
    mono = True
    for tau in ("0.5", "5", "50"):
        prev = None
        for n in (5, 10, 20, 30):
            T, Y = TY(n, iv.mpf(tau))
            if prev is not None:
                mono = mono and T.a >= prev.b * (1 - 1e-15)
            prev = T
    check(mono, "T_n(tau) = 3^n t_n(5^-n tau) increases with n at tau = 0.5, 5, 50 (intervals)")

    # certified c(0.45), c(0.8)
    t1 = time.time()
    c45, ncell = c_bracket("0.45", r2=r2)
    print(f"   c(0.45) in [{float(c45.a):.8f}, {float(c45.b):.8f}]  ({ncell} cells, {time.time()-t1:.1f}s)")
    check(c45.a > 2 and c45.b < (2 + B45).a and c45.a < tab[("0.45", 30)] < c45.b,
          "c(0.45) certified: 2 < c(0.45) <= 2 + B and the recursion value 3.7255634 lies in the enclosure")
    c80, _ = c_bracket("0.8", r2=r2)
    print(f"   c(0.8)  in [{float(c80.a):.8f}, {float(c80.b):.8f}]")
    check(c80.b < c45.a and c80.a > 2 and c80.a < tab[("0.8", 30)] < c80.b,
          "Kept Failure BD.1: c(0.8) < c(0.45), both > 2 = lim 3^n t_n(0+): the 3^-n constant depends on a")
    print(f"   limit tan theta = -1/(3 sqrt3) = {-1/(3*math.sqrt(3)):.12f}; 3^n(tan theta_n + 1/(3 sqrt3)) -> -4c/(9 sqrt3) "
          f"= {-4*float(tab[('0.45', 30)])/(9*math.sqrt(3)):.8f} at a=0.45 (float)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BD checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
