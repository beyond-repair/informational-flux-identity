#!/usr/bin/env python3
r"""Theorem CM: for every k >= 1 the open disk |z| < rho(L_k) is filled by C^k eigenvalues of CE's transfer operator,
each of infinite multiplicity (the C^k analogue of CJ.4, proved for all k at once without CJ's cone certificate).  With
CL, the Fredholm essential spectrum of Lt on C^k(X) is exactly the closed disk |z| <= rho(L_k), so
spec(Lt; C^2) = {r, l_2} u {|z| <= rho(L_2)} and spec(Lt; C^3), spec(Lt; C^4) = {r, l_2, l_3} u {|z| <= rho(L_k)}.

Setting (GASKET.md, Theorems CE-CL): in s = sqrt(t + 25/4), X = [sqrt5/2, 5/2],
    (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s),   sigma_(1,2)(s) = sqrt(15/4 -+ s),
    L_k f = sum_i wt_i |sigma_i'|^k f(sigma_i),   Y_i = sigma_i X,  G = X \ (Y_1 u Y_2) = (g_a, g_b),  s0 = sqrt5/2.
    T_1(y) = 15/4 - y^2 on Y_1 and T_2(y) = y^2 - 15/4 on Y_2 invert sigma_1, sigma_2.

Theorem CM (checked here; exit 1 on any failure):
  CM.1 [proved; Harnack constants certified in arb] For every k >= 1 and every a < rho(L_k) there is a real-analytic
       h > 0 on a neighbourhood of X with L_k h >= a h on X: h = sum_{j<J} L_k^j 1 / t^j, t > rho(L_k) close to it,
       J large.  Uses only ||L_k^j|| >= rho(L_k)^j and the Harnack bound L_k^j 1(s) <= C_k L_k^j 1(s'),
       C_k = exp((D_w + k D_s)|X|) (CL.5's distortion constants), giving L_k h / h >= t - C_k (t - rho) - o_J(1).
  CM.2 [proved; sympy identities] With c_i = |sigma_i'|^k h(sigma_i) / L_k h (so sum_i wt_i c_i = 1), the operator
       K f = c_i(T_i .) f(T_i .) on Y_i, and on G the degree-(2k+1) Hermite interpolant of the one-sided k-jets at
       g_a, g_b (both read from the k-jet of f at s0 = T_1 g_a = T_2 g_b), is bounded on C^k(X) and Lt K = I.  Its top
       derivative is (K f)^(k) = N f^(k) + (compact) + (finite rank), with N g(y) = c_i(T_i y) T_i'(y)^k g(T_i y) on Y_i,
       0 on G, acting on bounded Borel functions B(X), and |c_i(T_i y)| |T_i'(y)|^k = h(y) / L_k h(T_i y) exactly.
       Telescoping: ||N^n|| <= (max h / min h) a^(-n), so r(N) <= 1/a.
  CM.3 [proved from CM.1-CM.2 and the analytic Fredholm theorem; textbook] For every 0 < |mu| < rho(L_k),
       dim ker(Lt - mu; C^k(X)) = infinity.  (On F = C^k (+) B(X), K lifts to (0 (+) N) + compact; (I - mu K)^(-1) is
       meromorphic on |mu| < a with finite-rank principal parts and maps the image of C^k into itself; for a gap bump
       g, Lt (I - mu K)^(-1) g = mu (I - mu K)^(-1) g, and at the discrete poles the bumps in a finite-codimension
       subspace still give regular values.)  Also mu = 0 (gap bumps).
  CM.4 [proved from CM.3, CL.2, CK.3, CL.3, CG.3] sigma_e(Lt; C^k(X)) = {|z| <= rho(L_k)} (Fredholm) for every k >= 1,
       and spec(Lt; C^k) = {|z| <= rho(L_k)} u {isolated eigenvalues of modulus > rho(L_k)}: {r, l_2} for k = 2,
       {r, l_2, l_3} for k = 3, 4.  For k = 1 this re-proves CJ.4 without CJ's cone / gap-weight certificates.
Not claimed: whether points of the circle |z| = rho(L_k) are eigenvalues (for any k); B != 0 and an explicit C_G for the
level sums; a rate for inf R_n.  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not
a flux.  Not thrust.
"""
from __future__ import annotations

import math
import sys
import time

import sympy as sp
from flint import arb, ctx

import gasket_fractional_profile_level_calkin as cl

FAILS: list[str] = []
CW = cl.CW                                     # CL.3 two-sided Collatz-Wielandt brackets for k = 2, 3, 4
LAM1 = cl.LAM1                                 # CI.1 bracket for rho(L_1)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- CM.2 exact identities (sympy)
def symbolic_checks(kmax: int = 4):
    s, y = sp.symbols("s y", positive=True)
    half = sp.Rational(1, 2)
    W = [(s + half) / (2 * s * (sp.Rational(7, 2) - s)), (s - half) / (2 * s * (s + sp.Rational(7, 2)))]
    S = [sp.sqrt(sp.Rational(15, 4) - s), sp.sqrt(sp.Rational(15, 4) + s)]
    T = [sp.Rational(15, 4) - y ** 2, y ** 2 - sp.Rational(15, 4)]
    s0 = sp.sqrt(5) / 2
    ok = all(sp.simplify(T[i].subs(y, S[i]) - s) == 0 for i in (0, 1))
    ga, gb = S[0].subs(s, s0), S[1].subs(s, s0)
    ok = ok and sp.simplify(T[0].subs(y, ga) - s0) == 0 and sp.simplify(T[1].subs(y, gb) - s0) == 0
    ok = ok and sp.simplify(T[0].subs(y, s0) - sp.Rational(5, 2)) == 0 and bool(sp.N(ga) < sp.N(gb))
    check(bool(ok), "CM.2: T_1 = 15/4 - y^2, T_2 = y^2 - 15/4 invert sigma_1, sigma_2 (polynomials, analytic); "
                    "T_1(g_a) = T_2(g_b) = s0 = sqrt5/2 and T_1(s0) = 5/2: the gap ends' one-sided jets of K f are read "
                    "from the k-jet of f at s0")
    h, H = sp.Function("h"), sp.Function("H")   # H stands for L_k h
    ok = True
    for k in range(1, kmax + 1):
        Lkh = sum(W[i] * (sp.Abs(sp.diff(S[i], s))) ** k * h(S[i]) for i in (0, 1))
        c = [(sp.Abs(sp.diff(S[i], s))) ** k * h(S[i]) / Lkh for i in (0, 1)]
        ok = ok and sp.simplify(sum(W[i] * c[i] for i in (0, 1)) - 1) == 0
        for i in (0, 1):
            ci = (1 / (2 * S[i])) ** k * h(S[i]) / H(s)          # |sigma_i'| = 1/(2 sigma_i)
            n_i = ci.subs(s, T[i]) * sp.Abs(sp.diff(T[i], y)) ** k
            ok = ok and sp.simplify(n_i - h(y) / H(T[i])) == 0
    check(bool(ok), f"CM.2: sum_i wt_i c_i = 1 for c_i = |sigma_i'|^k h(sigma_i)/L_k h and generic h, and "
                    f"|c_i(T_i y)| |T_i'(y)|^k = h(y) / L_k h(T_i y) exactly (k = 1..{kmax}); so N^n telescopes")
    cc, lam = sp.Function("c"), sp.Symbol("lam")
    ok = True
    for k in range(1, kmax + 1):
        for i in (0, 1):
            # f = exp(lam x): the lam^j coefficient of (c(T) f(T))^(k) e^(-lam T) is the coefficient of f^(j)(T)
            expr = sp.expand(sp.simplify(sp.diff(cc(T[i]) * sp.exp(lam * T[i]), y, k) * sp.exp(-lam * T[i])))
            poly = sp.Poly(expr, lam)
            ok = ok and poly.degree() == k
            ok = ok and sp.simplify(poly.coeff_monomial(lam ** k) - cc(T[i]) * sp.diff(T[i], y) ** k) == 0
    check(bool(ok), f"CM.2: Leibniz/Faa di Bruno: (c(T_i y) f(T_i y))^(k) = c(T_i y) T_i'(y)^k f^(k)(T_i y) + "
                    f"sum_(j<k) b_ij(y) f^(j)(T_i y) with b_ij analytic on Y_i (k = 1..{kmax})")
    ok = True
    for k in range(1, kmax + 1):
        x = sp.Symbol("x")
        basis = [x ** p for p in range(2 * k + 2)]
        rows = [[sp.diff(b, x, d).subs(x, pt) for b in basis] for pt in (ga, gb) for d in range(k + 1)]
        det = sp.Matrix(rows).det()
        ok = ok and abs(float(sp.N(det))) > 0
    check(bool(ok), f"CM.2: the degree-(2k+1) Hermite interpolation of k-jets at g_a != g_b is uniquely solvable "
                    f"(confluent Vandermonde det != 0, k = 1..{kmax}), so K f is C^k across both gap ends")


# ----------------------------------------------------------------- CM.1 Harnack constants (arb)
def harnack_constants(n_cover: int = 2000):
    Dw, Ds = cl.distortion_constants(n_cover=n_cover)
    ctx.prec = 64
    Xlen = arb(5) / 2 - arb(5).sqrt() / 2
    Ck = {}
    for k in range(1, 5):
        Ck[k] = (Xlen * (arb(Dw) + k * arb(Ds))).exp()
    ok = all(bool(Ck[k] < 1e3) for k in Ck)
    check(ok, "CM.1: Harnack: L_k^j 1(s) <= C_k L_k^j 1(s') for all j and s, s' in X, C_k = exp((D_w + k D_s)|X|) <= " +
          ", ".join(f"{float(Ck[k].upper()):.3f}" for k in Ck) + " (k = 1..4; uniform in j)")
    return {k: float(Ck[k].upper()) for k in Ck}


# ----------------------------------------------------------------- consequences (constants of CL.3, CG.3, CI.1)
def consequences():
    ctx.prec = 128
    l2lo = arb(cl.L2_ABS) - arb(cl.L2_ERR)
    l3lo = arb(cl.L3_ABS) - arb(cl.L3_ERR)
    l3hi = arb(cl.L3_ABS) + arb(cl.L3_ERR)
    up = {k: arb(CW[k][1]) for k in CW}
    lo = {k: arb(CW[k][0]) for k in CW}
    ok = bool(up[2] < l2lo) and bool(l3hi < lo[2]) and bool(up[3] < l3lo) and bool(up[4] < l3lo)
    ok = ok and bool(arb(LAM1[0]) > l2lo + 2 * arb(cl.L2_ERR))
    check(ok, "CM.4: spec(C^1) = {r} u {|z| <= rho(L_1)}; spec(C^2) = {r, l_2} u {|z| <= rho(L_2)} with "
              "|l_3| < 0.0430702 <= rho(L_2) (l_3 is one of the infinitely many C^2 eigenvalues of the disk); "
              "spec(C^3), spec(C^4) = {r, l_2, l_3} u {|z| <= rho(L_k)}, rho(L_3) <= 0.0139193, rho(L_4) <= 0.0045265")


# ----------------------------------------------------------------- floats (no claim): attempts to break CM
class _Jets:
    """Truncated Taylor jets (order k) of functions of y, vectorised over points."""

    def __init__(self, k, coef):
        self.O = k + 1
        self.k = k
        self.coef = coef

    def const(self, v, n):
        import numpy as np
        a = np.zeros((self.O, n), dtype=complex if isinstance(v, complex) else float)
        a[0] = v
        return a

    def mul(self, a, b):
        import numpy as np
        c = np.zeros(np.broadcast_shapes(a.shape, b.shape), dtype=np.result_type(a, b))
        for i in range(self.O):
            for j in range(self.O - i):
                c[i + j] += a[i] * b[j]
        return c

    def inv(self, a):
        import numpy as np
        b = np.zeros_like(a)
        b[0] = 1 / a[0]
        for n in range(1, self.O):
            b[n] = -sum(a[j] * b[n - j] for j in range(1, n + 1)) / a[0]
        return b

    def comp(self, a, derivs):
        import numpy as np
        d = a.copy()
        d[0] = 0
        out = np.zeros_like(a)
        p = self.const(1.0, a.shape[1])
        fact = 1
        for j in range(self.O):
            out = out + derivs[j] / fact * p
            p = self.mul(p, d)
            fact *= j + 1
        return out

    def sqrt(self, a):
        x = a[0]
        ds, c = [x ** 0.5], 0.5
        for j in range(1, self.O):
            ds.append(c * x ** (0.5 - j))
            c *= 0.5 - j
        return self.comp(a, ds)

    def pw(self, a, k):
        r = self.const(1.0, a.shape[1])
        for _ in range(k):
            r = self.mul(r, a)
        return r

    def h(self, a):
        n = a.shape[1]
        u = (a - self.const(1.8, n)) / 0.7
        acc = self.const(0.0, n)
        for c in reversed(self.coef):
            acc = self.mul(acc, u) + self.const(c, n)
        return acc

    def Lkh(self, a):
        n = a.shape[1]
        w1 = self.mul(a + self.const(.5, n), self.inv(self.mul(2 * a, self.const(3.5, n) - a)))
        w2 = self.mul(a - self.const(.5, n), self.inv(self.mul(2 * a, a + self.const(3.5, n))))
        tot = self.const(0.0, n)
        for wi, sg in ((w1, -1), (w2, 1)):
            yy = self.sqrt(self.const(3.75, n) + sg * a)
            tot = tot + self.mul(wi, self.mul(self.h(yy), self.inv(self.pw(2 * yy, self.k))))
        return tot

    def bump(self, a, A=1.7, B=2.1):
        import numpy as np
        n = a.shape[1]
        e = -self.inv(self.mul(a - self.const(A, n), self.const(B, n) - a))
        ex = np.exp(e[0])
        return self.comp(e, [ex] * self.O)


def _term_norms(J, mu, m, nx=40, A=1.7, B=2.1):
    """sup_X |(mu^m K^m g)^(j)|, j = 0..k, for the gap bump g on [A, B] and K built from h (exact pointwise formula:
    K^m g(y) = prod_(l<m) c(T^l y) g(T^m y) where the T-orbit of y first enters G at time m)."""
    import numpy as np
    xs = np.linspace(A, B, nx + 2)[1:-1]
    pts = [xs]
    for _ in range(m):
        pts = [np.sqrt(3.75 - p) for p in pts] + [np.sqrt(3.75 + p) for p in pts]
    yv = np.concatenate(pts)
    n = len(yv)
    cur = J.const(0.0, n)
    cur[0], cur[1] = yv, 1
    prod = J.const(1.0, n)
    for _ in range(m):
        sgn = np.where(cur[0] < 2.0, 1.0, -1.0)
        T = sgn * (J.const(3.75, n) - J.mul(cur, cur))
        c = J.mul(J.mul(J.h(cur), J.inv(J.Lkh(T))), J.inv(J.pw(2 * cur, J.k)))
        prod = J.mul(prod, c)
        cur = T
    v = J.mul(prod, J.bump(cur, A, B)) * mu ** m
    return [float(np.abs(v[j]).max()) * math.factorial(j) for j in range(J.O)]


def _N_norm(J, a_lo, nmax=10, npts=4001):
    """sup_y prod_(l<n) h(T^l y) / L_k h(T^(l+1) y) over orbits staying in Y_1 u Y_2 (= ||N^n|| on B(X))."""
    import numpy as np
    xs = np.linspace(5 ** 0.5 / 2, 2.5, npts)
    pts = xs.copy()
    hv = lambda x: J.h(np.vstack([x] + [np.zeros_like(x)] * J.k))[0]
    Lv = lambda x: J.Lkh(np.vstack([x] + [np.zeros_like(x)] * J.k))[0]
    worst = 0.0
    hmax, hmin = hv(xs).max(), hv(xs).min()
    for n in range(1, nmax + 1):
        # points y = sigma_omega(x) for |omega| = n give every orbit that stays n steps in Y
        P = [xs]
        for _ in range(n):
            P = [np.sqrt(3.75 - p) for p in P] + [np.sqrt(3.75 + p) for p in P]
            if sum(len(p) for p in P) > 400000:
                P = [p[::4] for p in P]
        y = np.concatenate(P)
        prod = np.ones_like(y)
        cur = y
        for _ in range(n):
            T = np.where(cur < 2.0, 3.75 - cur ** 2, cur ** 2 - 3.75)
            prod *= hv(cur) / Lv(T)
            cur = T
        worst = max(worst, prod.max() * a_lo ** n / (hmax / hmin))
    return worst


def float_observations(ks=(2, 3)):
    import numpy as np
    out = {}
    for k in ks:
        lam, coef = cl._float_h(k)
        J = _Jets(k, [float(c) for c in coef])
        res = {}
        rho_k = float(CW[k][0])
        for mu in (0.7 * rho_k, -0.93 * rho_k, 0.9j * rho_k, 1.16 * rho_k):
            rows = [_term_norms(J, mu, m) for m in (6, 10)]
            ratio = [(rows[1][j] / rows[0][j]) ** 0.25 for j in range(k + 1)]
            res[mu] = ratio
            print(f"     (observation) k = {k}, mu = {mu:.5f}: growth per level of ||(mu^m K^m g)^(j)||, j = 0..{k}: " +
                  ", ".join(f"{x:.4f}" for x in ratio) + f"  (|mu|/rho(L_k) = {abs(mu) / lam:.4f})")
        out[k] = (lam, res)
        out[("N", k)] = _N_norm(J, float(CW[k][0]))
        print(f"     (observation) k = {k}: max_n<=10 ||N^n|| a_k^n / (max h/min h) = {out[('N', k)]:.4f} (<= 1 by CM.2)")
    # Lt v = mu v for the truncated explicit eigenfunction (k = 2, mu = 0.03), at random points
    lam, coef = cl._float_h(2)
    J = _Jets(2, [float(c) for c in coef])
    rng = np.random.default_rng(1)
    GA, GB = math.sqrt(3.75 - 5 ** 0.5 / 2), math.sqrt(3.75 + 5 ** 0.5 / 2)

    def v_val(y, mu, M=40):
        a = J.const(0.0, 1)
        a[0] = y
        prod = 1.0
        for m in range(M + 1):
            if 1.7 < a[0, 0] < 2.1:
                return prod * mu ** m * float(J.bump(a)[0, 0])
            yv = a[0, 0]
            if 5 ** 0.5 / 2 <= yv <= GA:          # Y_1
                T = 3.75 - yv ** 2
            elif GB <= yv <= 2.5:                 # Y_2
                T = yv ** 2 - 3.75
            else:                                 # in G but outside the bump support: v = 0 there
                return 0.0
            c = float(J.h(a)[0, 0] / J.Lkh(J.const(T, 1))[0, 0]) / (2 * a[0, 0]) ** 2
            prod *= c
            a = J.const(T, 1)
        return 0.0
    err = 0.0
    for s in rng.uniform(5 ** 0.5 / 2, 2.5, 200):
        w1 = (s + .5) / (2 * s * (3.5 - s))
        w2 = (s - .5) / (2 * s * (s + 3.5))
        lhs = w1 * v_val(np.sqrt(3.75 - s), 0.03) + w2 * v_val(np.sqrt(3.75 + s), 0.03)
        err = max(err, abs(lhs - 0.03 * v_val(s, 0.03)))
    out["eig_err"] = err
    print(f"     (observation) explicit v = sum mu^m K^m g (k = 2, mu = 0.03): max |Lt v - mu v| at 200 random s = {err:.1e}")
    # Harnack sanity: max/min of L_k^j 1 for j <= 10
    xs = np.linspace(5 ** .5 / 2, 2.5, 301)
    hr = {}
    for k in (1, 2, 3, 4):
        st = cl._words(10, xs)
        Lk = sum(Wt * np.abs(D) ** k for _, Wt, D, _ in st)
        hr[k] = Lk.max() / Lk.min()
    out["harnack"] = hr
    print("     (observation) max/min of L_k^10 1 = " + ", ".join(f"{hr[k]:.3f}" for k in hr) + " (k = 1..4)")
    # CM.1 construction in floats: h = sum_{j<J} L_k^j 1/t^j, t = rho (1 + d)
    cm1 = {}
    for k in (2,):
        lam_f = out[k][0]
        for d in (0.05, 0.01):
            t = lam_f * (1 + d)
            st_all = [cl._words(j, xs) for j in range(0, 13)]
            hj = sum(sum(Wt * np.abs(D) ** k for _, Wt, D, _ in st_all[j]) / t ** j for j in range(12))
            st1 = st_all[1]
            # L_k h via one more step of words: L_k h(s) = sum_{j<J} L_k^{j+1} 1 / t^j
            Lh = sum(sum(Wt * np.abs(D) ** k for _, Wt, D, _ in st_all[j + 1]) / t ** j for j in range(12))
            cm1[d] = float((Lh / hj).min() / lam_f)
        print(f"     (observation) CM.1 in floats, k = 2, J = 12: min L_2 h / h / rho(L_2) = " +
              ", ".join(f"{cm1[d]:.4f} (t = rho(1+{d}))" for d in cm1))
    out["cm1"] = cm1
    return out


def main() -> int:
    t0 = time.time()
    print("Theorem CM: every |mu| < rho(L_k) is a C^k eigenvalue of infinite multiplicity (all k >= 1); the C^k "
          "Fredholm essential spectrum is the closed disk |z| <= rho(L_k). Not thrust.")
    symbolic_checks()
    Ck = harnack_constants()
    for k in (2, 3, 4):
        cl.collatz_two_sided(k)
    consequences()
    try:
        out = float_observations()
        for k in (2, 3):
            lam, res = out[k]
            ok = True
            for mu, ratio in res.items():
                q = abs(mu) / lam
                ok = ok and abs(ratio[k] - q) < 0.02 * q
            check(ok, f"CM.3 (float sanity, k = {k}): the top-derivative terms of v = sum mu^m K^m g grow per level by "
                      f"|mu|/rho(L_{k}) (to 2%), so v is C^{k} for |mu| < rho(L_{k}) and the series is not C^{k} beyond")
            check(out[("N", k)] <= 1 + 1e-9, f"CM.2 (float sanity, k = {k}): ||N^n|| <= (max h/min h) a_k^(-n)")
        check(out["eig_err"] < 1e-12, "CM.3 (float sanity): the explicit C^2 eigenfunction satisfies Lt v = mu v")
        check(all(out["harnack"][k] <= Ck[k] for k in (1, 2, 3, 4)),
              "CM.1 (float sanity): max/min of L_k^10 1 stays below the certified Harnack constant C_k")
    except ImportError as e:  # pragma: no cover
        print("     (observation skipped:", e, ")")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CM checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
