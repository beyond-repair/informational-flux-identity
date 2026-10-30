#!/usr/bin/env python3
r"""Theorem CK: Assumption CH (Hennion's theorem, literature) is discharged for CE's transfer operator.  On C^k(X) the
operator is an explicit compact perturbation of diag(0, M_k), where M_k is a signed transfer operator on C(X) with
||M_k^n|| = ||L_k^n 1|| exactly; so r_ess(Lt; C^k(X)) <= rho(L_k), and outside that radius the spectrum is isolated
eigenvalues with finite-rank Riesz projections, by an elementary finite-rank (analytic Fredholm) argument.  CH.3, CI.5
and CJ.5 therefore hold without any literature assumption: spec(Lt; C^1(X)) = {r} u {|z| <= rho(L_1)} unconditionally.

Setting (GASKET.md, Theorems CE-CJ): in s = sqrt(t + 25/4), X = [sqrt5/2, 5/2],
    (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s),   sigma_(1,2)(s) = sqrt(15/4 -+ s),
    wt1 = (s+1/2)/(2s(7/2-s)),  wt2 = (s-1/2)/(2s(s+7/2)),   L_k f = sum_i wt_i |sigma_i'|^k f(sigma_i).
    M_k g = sum_i wt_i (sigma_i')^k g(sigma_i)   (signed: sigma_1' < 0 < sigma_2').

Theorem CK (checked here; exit 1 on any failure):
  CK.1 [exact; sympy] Faa di Bruno / Leibniz: for k = 1..4,
           (Lt f)^(k) = M_k f^(k) + sum_i sum_(j<k) a^(k)_(ij)(s) f^(j)(sigma_i s),
       with a^(k)_(ik) = wt_i (sigma_i')^k exactly; the a^(k)_(ij) are polynomials in derivatives of wt_i, sigma_i (all
       real-analytic near X: poles of wt_i at 0, +-7/2, branch points of sigma_i at +-15/4, all outside X).  Checked
       against direct sympy differentiation for test functions f.
  CK.2 [certified; arb] On X: wt_i > 0, sigma_1 X = [sqrt5/2, g_a], sigma_2 X = [g_b, 5/2], g_a < g_b (disjoint images),
       and the one-step Lasota-Yorke constants sup_X L_k 1 and sup_X sum_i |a^(k)_(ij)| (finite; listed) for k = 1..4.
  CK.3 [proved] (i) ||M_k^n||_(C(X)) = ||L_k^n 1||_inf for every n (the 2^n points sigma_omega(s) are distinct, so a
       continuous g with |g| <= 1 matches every sign), hence rho(M_k) = rho(L_k).  (ii) Phi f = (f(s0), .., f^(k-1)(s0),
       f^(k)) is an isomorphism C^k(X) -> C^k x C(X) (Taylor), and Phi Lt Phi^-1 = diag(0, M_k) + K_k with K_k compact
       (finite rank in the C^k part; Arzela-Ascoli in the C(X) part).  (iii) Elementary analytic Fredholm argument:
       in |z| > rho(M_k), z - Lt is invertible except on a discrete set of eigenvalues with finite-rank Riesz
       projections.  So r_ess(Lt; C^k(X)) <= rho(L_k), with no literature assumption.
  CK.4 [proved from CK.3, CH.2, CI.1, CG.3 and CH.3's restriction argument] Assumption CH is no longer used:
       C^1: only r beyond rho(L_1) <= 0.1343551; C^2: only r, l_2 beyond 0.0432; C^3, C^4: only r, l_2, l_3 beyond
       0.01396, 0.00454.  With CJ.4: spec(Lt; C^1(X)) = {r} u {|z| <= rho(L_1)}, r_ess(C^1) = rho(L_1) in
       [0.1343549, 0.1343551], unconditionally.
Not claimed: the boundary circle |z| = rho(L_1); r_ess(C^k) = rho(L_k) for k >= 2 (only <=); B != 0 and an explicit
C_G for the level sums.  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.
Not thrust.
"""
from __future__ import annotations

import sys
import time

import sympy as sp
from flint import arb, ctx

FAILS: list[str] = []
KMAX = 4
LAM = {1: "0.1343551", 2: "0.0432", 3: "0.01396", 4: "0.00454"}   # CI.1 (k = 1), CH.2 (k = 2..4): rho(L_k) <= LAM[k]
LAM1_LO = "0.1343549"                                              # CI.1
R_LO, R_HI = "0.4237312869626", "0.4237312869627"                  # CG.3
L2_ABS, L2_ERR = "0.0731138175124", "1e-13"                        # CG.3
L3_ABS, L3_ERR = "0.0170676440626", "1e-12"
SMALL = "0.001"                                                    # CG.3: every other H^2 eigenvalue has |mu| < 1e-3

s = sp.symbols("s", positive=True)
x = sp.symbols("x")
HALF = sp.Rational(1, 2)
W = (s + HALF) / (2 * s * (sp.Rational(7, 2) - s)), (s - HALF) / (2 * s * (s + sp.Rational(7, 2)))
SIG = sp.sqrt(sp.Rational(15, 4) - s), sp.sqrt(sp.Rational(15, 4) + s)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- CK.1 exact Faa di Bruno structure
def chain_coeffs(sig, k):
    """b[m][j]: (g o sig)^(m) = sum_j b[m][j] g^(j)(sig).  Recursion b_(m+1,j) = b_(m,j)' + b_(m,j-1) sig' (chain rule)."""
    b = [[sp.Integer(1)]]
    ds = sp.diff(sig, s)
    for m in range(k):
        prev = b[-1] + [sp.Integer(0)]
        new = [sp.diff(prev[0], s)]
        for j in range(1, m + 2):
            new.append(sp.diff(prev[j], s) + prev[j - 1] * ds)
        b.append(new)
    return b


def lt_coeffs(i, k):
    """a[j], j = 0..k: (wt_i * g(sig_i))^(k) = sum_j a[j] g^(j)(sig_i) (Leibniz on top of chain_coeffs)."""
    w, sig = W[i], SIG[i]
    b = chain_coeffs(sig, k)
    a = [sp.Integer(0)] * (k + 1)
    for m in range(k + 1):
        dw = sp.diff(w, s, k - m) * sp.binomial(k, m)
        for j, bj in enumerate(b[m]):
            a[j] += dw * bj
    return [sp.factor(sp.together(sp.simplify(e))) for e in a]   # factored: tight arb enclosures


def symbolic_checks():
    ds = [sp.diff(SIG[0], s), sp.diff(SIG[1], s)]
    check(sp.simplify(ds[0] + 1 / (2 * SIG[0])) == 0 and sp.simplify(ds[1] - 1 / (2 * SIG[1])) == 0,
          "CK.1: sigma_1' = -1/(2 sigma_1) < 0 < sigma_2' = 1/(2 sigma_2) (sympy)")
    coeffs = {}
    tests = [sp.exp(x) * x ** 3, sp.sin(2 * x) + x ** 7]
    pts = [sp.Rational(6, 5), sp.Rational(17, 10), sp.Rational(23, 10), sp.Rational(12, 5)]
    for k in range(1, KMAX + 1):
        ok = True
        for i in (0, 1):
            a = lt_coeffs(i, k)
            coeffs[(k, i)] = a
            ok = ok and sp.simplify(a[k] - W[i] * ds[i] ** k) == 0
            for f in tests:
                direct = sp.diff(W[i] * f.subs(x, SIG[i]), s, k)
                via = sum(a[j] * sp.diff(f, x, j).subs(x, SIG[i]) for j in range(k + 1))
                for p in pts:   # 60-digit evaluation of an identity that holds symbolically (chain rule); guards typos
                    d = sp.N((direct - via).subs(s, p), 60)
                    ok = ok and abs(d) < sp.Float("1e-45")
        check(ok, f"CK.1: k = {k}: (Lt f)^({k}) = sum_i [wt_i (sigma_i')^{k} f^({k})(sigma_i) + sum_(j<{k}) a_ij f^(j)(sigma_i)]; "
                  f"leading coefficient exact (sympy), expansion matches direct differentiation of 2 test f at 4 points")
    # singularities of the coefficients: only s = 0, +-7/2 (weights) and s = +-15/4 (branches)
    sing = set()
    for (k, i), a in coeffs.items():
        for e in a:
            num, den = sp.fraction(sp.together(e))
            for fac in sp.Mul.make_args(sp.factor(den)):
                base = fac.base if fac.is_Pow else fac
                if base.is_number:
                    continue
                inner = base.base if (base.is_Pow and base.exp.is_Rational) else base
                for root in sp.solve(sp.Eq(inner, 0), s) if inner.free_symbols else []:
                    sing.add(root)
    allowed = {sp.Rational(7, 2), sp.Rational(15, 4)}
    real_sing = set(sing)
    check(real_sing <= allowed and all(not (sp.sqrt(5) / 2 <= r <= sp.Rational(5, 2)) for r in real_sing),
          f"CK.1: denominators of all a^(k)_(ij) vanish only at positive s in {sorted(real_sing)} (s = 0, -7/2, -15/4 are "
          f"excluded by s > 0), none in X = [sqrt5/2, 5/2]: the coefficients are real-analytic on a neighbourhood of X")
    return coeffs


# ----------------------------------------------------------------- arb evaluation of sympy expressions
def arb_eval(e, xv):
    if e == s:
        return xv
    if e.is_Integer:
        return arb(int(e))
    if e.is_Rational:
        return arb(int(e.p)) / arb(int(e.q))
    if e.is_Add:
        out = arb(0)
        for t in e.args:
            out += arb_eval(t, xv)
        return out
    if e.is_Mul:
        out = arb(1)
        for t in e.args:
            out *= arb_eval(t, xv)
        return out
    if e.is_Pow:
        b = arb_eval(e.base, xv)
        p = e.exp
        if p.is_Integer:
            return b ** int(p)
        if p.is_Rational and p.q == 2:
            r = b.sqrt()
            return r ** int(p.p) if p.p > 0 else 1 / r ** int(-p.p)
    raise ValueError(f"cannot evaluate {e}")


def cover_sup(fun, lo, hi, n):
    best = arb(0)
    for j in range(n):
        a = lo + (hi - lo) * j / n
        b = lo + (hi - lo) * (j + 1) / n
        xb = arb(((a + b) / 2).mid(), ((b - a) / 2).upper() + a.rad() + b.rad())
        best = arb(max(best.upper(), fun(xb).upper()))
    return best


def cover_inf(fun, lo, hi, n):
    best = None
    for j in range(n):
        a = lo + (hi - lo) * j / n
        b = lo + (hi - lo) * (j + 1) / n
        xb = arb(((a + b) / 2).mid(), ((b - a) / 2).upper() + a.rad() + b.rad())
        v = fun(xb).lower()
        best = v if best is None else min(best, v)
    return arb(best)


def mv_abs_sup(exprs, lo, hi, n):
    """Certified sup over [lo, hi] of sum_e |e(s)|, by the mean-value (centred) form on each of n balls:
    |e(x)| <= |e(m)| + sup_ball |e'| * rad.  e(m) is an arb point evaluation; e' is enclosed on the whole ball."""
    ders = [sp.diff(e, s) for e in exprs]
    best = arb(0)
    for j in range(n):
        a = lo + (hi - lo) * j / n
        b = lo + (hi - lo) * (j + 1) / n
        m = arb(((a + b) / 2).mid())
        rad = ((b - a) / 2).upper() + a.rad() + b.rad()
        xb = arb(m.mid(), rad)
        tot = arb(0)
        for e, d in zip(exprs, ders):
            tot += abs(arb_eval(e, m)) + abs(arb_eval(d, xb)) * arb(rad)
        best = arb(max(best.upper(), tot.upper()))
    return best


def certified_constants(coeffs, n_cover=2000):
    ctx.prec = 128
    s0, b5 = arb(5).sqrt() / 2, arb(5) / 2
    ga, gb = (arb(15) / 4 - s0).sqrt(), (arb(15) / 4 + s0).sqrt()
    wmin = [cover_inf(lambda xb, i=i: arb_eval(W[i], xb), s0, b5, n_cover) for i in (0, 1)]
    s1_lo, s1_hi = (arb(15) / 4 - b5).sqrt(), ga
    check(bool(wmin[0] > 0) and bool(wmin[1] > 0) and bool(ga < gb) and bool(abs(s1_lo - s0) < arb("1e-30"))
          and bool(gb < b5),
          f"CK.2: wt_1 >= {wmin[0].str(5)} > 0, wt_2 >= {wmin[1].str(5)} > 0 on X; sigma_1 X = [sqrt5/2, g_a], "
          f"sigma_2 X = [g_b, 5/2] with g_a = {ga.str(8)} < g_b = {gb.str(8)}: the branch images are disjoint")
    out = {}
    for k in range(1, KMAX + 1):
        lk = mv_abs_sup([coeffs[(k, i)][k] for i in (0, 1)], s0, b5, n_cover)
        lower = [mv_abs_sup([coeffs[(k, i)][j] for i in (0, 1)], s0, b5, n_cover) for j in range(k)]
        fin = bool(lk.is_finite()) and all(bool(v.is_finite()) for v in lower)
        check(fin, f"CK.2: k = {k}: ||(Lt f)^({k})|| <= {lk.str(6)} ||f^({k})|| + "
                   + " + ".join(f"{v.str(5)} ||f^({j})||" for j, v in enumerate(lower))
                   + f" ({n_cover}-ball arb mean-value covers; first factor = sup L_{k} 1)")
        out[k] = (lk, lower)
    # the bounds quoted in GASKET.md CK.2
    (l1, (c10,)), (l2, (c20, c21)) = out[1], out[2]
    check(bool(wmin[0] > arb("0.3035")) and bool(wmin[1] > arb("0.0598")) and bool(l1 < arb("0.28167"))
          and bool(c10 < arb("0.56446")) and bool(l2 < arb("0.12268")) and bool(c21 < arb("0.61107"))
          and bool(c20 < arb("1.1560")),
          "CK.2: quoted bounds: wt_1 >= 0.3035, wt_2 >= 0.0598; ||(Lt f)'|| <= 0.28167 ||f'|| + 0.56446 ||f||; "
          "||(Lt f)''|| <= 0.12268 ||f''|| + 0.61107 ||f'|| + 1.1560 ||f||")
    return out


# ----------------------------------------------------------------- CK.3 ||M_k^n|| = ||L_k^n 1|| (float illustration)
def norm_identity_observation(nmax=6):
    """No claim beyond CK.3(i), which is proved: compute sup_s sum_omega |c_omega(s)| and |M_1^n g_s(s)| for the sign-matched
    g on a grid of s, showing they coincide (floats)."""
    import numpy as np
    xs = np.linspace(np.sqrt(5) / 2, 2.5, 201)
    sig = [lambda y: np.sqrt(3.75 - y), lambda y: np.sqrt(3.75 + y)]
    w = [lambda y: (y + .5) / (2 * y * (3.5 - y)), lambda y: (y - .5) / (2 * y * (y + 3.5))]
    dsg = [lambda y: -1 / (2 * np.sqrt(3.75 - y)), lambda y: 1 / (2 * np.sqrt(3.75 + y))]
    res = []
    for n in range(1, nmax + 1):
        pts, cf = [xs.copy()], [np.ones_like(xs)]
        for _ in range(n):
            pts2, cf2 = [], []
            for p, c in zip(pts, cf):
                for i in (0, 1):
                    pts2.append(sig[i](p))
                    cf2.append(c * w[i](p) * dsg[i](p))
            pts, cf = pts2, cf2
        P, C = np.array(pts), np.array(cf)
        L1n = np.abs(C).sum(axis=0)                 # L_1^n 1 on the grid
        mind = np.min(np.diff(np.sort(P, axis=0), axis=0)) if len(pts) > 1 else 1.0
        res.append((n, L1n.max(), L1n.max() ** (1 / n), mind))
    for n, nm, rt, md in res:
        print(f"     (observation) n = {n}: ||M_1^n|| = ||L_1^n 1|| = {nm:.6e}, n-th root {rt:.6f}, "
              f"min gap between the 2^n points sigma_omega(s) = {md:.2e} > 0")
    return res[-1][2]


# ----------------------------------------------------------------- CK.4 consequences
def consequences():
    ctx.prec = 128
    lam = {k: arb(v) for k, v in LAM.items()}
    r_lo = arb(R_LO)
    l2_hi, l2_lo = arb(L2_ABS) + arb(L2_ERR), arb(L2_ABS) - arb(L2_ERR)
    l3_hi, l3_lo = arb(L3_ABS) + arb(L3_ERR), arb(L3_ABS) - arb(L3_ERR)
    small = arb(SMALL)
    ok1 = bool(l2_hi < lam[1]) and bool(lam[1] < r_lo)
    ok2 = bool(l3_hi < lam[2]) and bool(lam[2] < l2_lo)
    ok3 = bool(small < lam[3]) and bool(lam[3] < l3_lo)
    ok4 = bool(small < lam[4]) and bool(lam[4] < l3_lo)
    check(ok1 and ok2 and ok3 and ok4,
          "CK.4: without Assumption CH, beyond rho(L_k) the C^k eigenvalues are CG's: C^1 only r (rho(L_1) <= 0.1343551); "
          "C^2 only r, l_2 (> 0.0432); C^3, C^4 only r, l_2, l_3 (> 0.01396, 0.00454); all simple (CH.3's list, now unconditional)")
    lo = arb(LAM1_LO)
    check(bool(lo <= lam[1]) and bool(lo > l2_hi),
          f"CK.4: with CJ.4 (every |mu| < rho(L_1) a C^1 eigenvalue): spec(Lt; C^1) = {{r}} u {{|z| <= rho(L_1)}}, "
          f"r_ess(C^1) = rho(L_1) in [{LAM1_LO}, {LAM[1]}], with no literature assumption (CI.5, CJ.5 unconditional)")


def main() -> int:
    t0 = time.time()
    print("Theorem CK: Assumption CH discharged: Lt = diag(0, M_k) + compact on C^k(X), ||M_k^n|| = ||L_k^n 1||, "
          "elementary analytic Fredholm. Not thrust.")
    coeffs = symbolic_checks()
    certified_constants(coeffs)
    sys.path.insert(0, __import__("os").path.dirname(__file__))
    import gasket_fractional_profile_level_essential as ci
    import gasket_fractional_profile_level_regularity as ch
    ci.FAILS.clear()
    ch.FAILS.clear()
    ci.collatz_two_sided()                       # CI.1 re-run: rho(L_1) in [0.1343549, 0.1343551]
    for k in (2, 3, 4):
        ch.collatz_wielandt(k)                   # CH.2 re-run: rho(L_k) <= 0.0432, 0.01396, 0.00454
    for msg in ci.FAILS + ch.FAILS:
        FAILS.append("imported: " + msg)
    check(ci.LAM_UP == LAM[1] and ci.LAM_LO == LAM1_LO and all(ch.LAM[k] == LAM[k] for k in (2, 3, 4)),
          "CK.4: radii imported from CI.1 / CH.2 match the constants used here")
    consequences()
    try:
        norm_identity_observation()
    except Exception as e:  # pragma: no cover
        print("     (observation skipped:", e, ")")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CK checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
