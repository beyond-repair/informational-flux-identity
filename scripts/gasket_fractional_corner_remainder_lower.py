#!/usr/bin/env python3
"""Theorem ED: EC's 3^-n is attained.  AZ (1)'s defect eps_n is POSITIVE and >= L_alpha (rho/3)^n for every n >= 0 and every
alpha in (gamma,1), so with EC.2 (|eps_n| <= B_alpha (rho/3)^n) rho^-n eps_n = Theta(3^-n), EB.1's remainder R_n = Theta(3^-n), and
EC.3's thresholds are sharp: the J-term expansion of X_n holds IFF rho^J < 3; above 4 gamma/3 the third term fails and Y_n
approaches Y_oo from BELOW.

Setting (GASKET.md, Theorems AQ, AZ, BK, DZ, EB, EC): alpha in (gamma,1), rho = 5^alpha/3, rho' = rho/3, p(s) = s(s+5), q = p^-1,
beta(s) = (s+6)/((s+2)(s+5)), e(s) = 2s/(2+s), r_(k+1)(s) = r_k(ps)/(beta(s) - e(s) r_k(ps)), r_0 = 1/(3+s); J_k = int s^-alpha r_k,
D'_k = rho^-k J_k, omega_k = BK.1's weight, eps_k = J_k - c rho^k/(rho-1) - C' (AZ (1)), X_k, Y_k, Y_oo, R_k as in EB.

Theorem ED (checked here; exit 1 on any failure):
  ED.1 (THEOREM, exact; every alpha in (gamma,1)).  Mt_k := int t^-alpha r_k(t) (omega_k(q(t)) - 1) dt >= 0 is nondecreasing in k,
       D'_(k+1) - D'_k = rho^-k Mt_k, Mt_k -> Mt_oo < oo, C' = -Mt_oo rho/(rho-1) (< 0: C_alpha < 0 again, without BN.4), and
       eps_n = rho^n sum_(k>=n) rho^-k (Mt_oo - Mt_k) >= 0 for every n >= 0 (AZ's E_n^- may be replaced by 0).
  ED.2 (THEOREM, explicit; every alpha in (gamma,1), every k >= 0).  Mt_(k+1) - Mt_k >= ell_alpha (rho/3)^k with
       ell_alpha = (1/3) int_0^x1 x^(1-alpha) (2 - U(x)) / (3 (3+U(x)) (3+U(x)/5) (1+x/25) (1+x/10)) dx,  U(x) = 5x exp(5x e^(1/10)/20),
       x1 = (2/5) e^(-1/10) (U(x1) = 2); hence eps_n >= L_alpha (rho/3)^n, L_alpha = 3 ell_alpha / (2 (1 - rho/3)) for every n >= 0.
       Certified: ell_gamma >= 0.0016681 (uniform: L_alpha >= 0.0025021 on (gamma,1)), L_0.8 >= 0.0055320, L_0.95 >= 0.0098674;
       with EC.2: (3/rho)^n eps_n in [L_alpha, B_alpha] for every n >= 6 (0.8: [0.005532, 1.94e3]).
  ED.3 (THEOREM; every alpha in (gamma,1), N >= 1, n >= N with 1 + nu rho^-n > 0).  R_n >= (1-delta_N) L_alpha 3^-n / D'_oo > 0, so
       R_n = Theta(3^-n), X_n < C_delta/(1 + nu rho^-n) strictly, and for every J >= 1:
       X_n = C_delta sum_(j<J) (-nu rho^-n)^j + o(rho^-(J-1)n)  IFF  rho^J < 3  (alpha < (J+1) gamma/J);
       for rho^J > 3, rho^((J-1)n) (X_n - sum) = -Theta((rho^J/3)^n) -> -oo.  J = 3: EB/EC's third term holds exactly on (gamma, 4gamma/3);
       for alpha in (4 gamma/3, 1), Y_n < Y_oo for all large n (from below; EC.3: from above below 4 gamma/3).
  ED.4 (certified at alpha = 0.8; arb).  For every n >= 60 and N = 1, 2, 3, 6, 10: 3^n R_n in [(1-delta_N) L_0.8/D'_oo, Lam1_N + Lam2_N (3/(5 rho))^60].
  ED.5 (OBSERVED, floats).  From J-quadrature differences (bias-free): Mt-increments (3/rho)^k (Mt_(k+1) - Mt_k) -> 0.2951 and
       (3/rho)^n eps_n ~ 0.741 at 0.8 (ED.2's L_0.8 is ~134x below): the bracket is honest but loose; no limit of 3^n R_n is claimed.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import contextlib
import io
import math
import sys
import time

import numpy as np

FAILS: list[str] = []
RESULT: dict = {}
NS = (1, 2, 3, 6, 10)


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- floats (numpy only)
def ell_float(al: float, M: int = 4000) -> float:
    """Lower Riemann sum of ED.2's ell_alpha (integrand = x^(1-alpha) * decreasing f), floats."""
    x1 = 0.4 * math.exp(-0.1)
    x = np.linspace(0.0, x1, M + 1)
    U = 5 * x * np.exp(5 * x * math.exp(0.1) / 20)
    f = (2 - U) / (3 * (3 + U) * (3 + U / 5) * (1 + x / 25) * (1 + x / 10))
    return float(np.sum(x[:-1] ** (1 - al) * f[1:] * np.diff(x)) / 3)


def float_increments(al: float, ks, h: float = 1e-3):
    """OBSERVED: Mt_k = (J_(k+1) - rho J_k)/rho from DZ.5-type quadrature (trapezoid in log s); returns
    (3/rho)^k (Mt_(k+1) - Mt_k) and (3/rho)^n eps_n with eps_n = rho^n sum_(k>=n) rho^-k (Mt_oo - Mt_k), Mt_oo ~ Mt_K."""
    import gasket_fractional_corner_limit_mobius as eb
    rho = 5**al / 3
    K = max(ks) + 8

    def J(n):
        with np.errstate(all="ignore"):
            u = np.arange(np.log(5.0 ** (-n - 40)), np.log(1e12), h)
            s = np.exp(u)
            r, _ = eb._rn(s, n)
            return np.trapezoid(s ** (1 - al) * r, u)

    Js = {k: J(k) for k in range(min(ks), K + 2)}
    Mt = {k: (Js[k + 1] - rho * Js[k]) / rho for k in range(min(ks), K + 1)}
    out = []
    for n in ks:
        d = (Mt[n + 1] - Mt[n]) * (3 / rho) ** n
        e = sum(rho ** (n - k) * (Mt[K] - Mt[k]) for k in range(n, K + 1)) * (3 / rho) ** n
        out.append((n, d, e))
    return out, -Mt[K] * rho / (rho - 1)


# ----------------------------------------------------------------- (A) sympy
def symbolic_checks():
    import sympy as sp
    s, m, a, x, t = sp.symbols("s m alpha x t", positive=True)
    beta = (s + 6) / ((s + 2) * (s + 5))
    rho = 5**a / 3
    check(sp.expand(3 * (s + 2) * (s + 5) - (5 + 2 * s) * (s + 6) - (s**2 + 4 * s)) == 0,
          "omega^0 >= 1: (1+2s/5)(5/3)beta(s) = (5+2s)(s+6)/(3(s+2)(s+5)) <= 1, since 3(s+2)(s+5) - (5+2s)(s+6) = s^2 + 4s")
    dlog = sp.together(sp.diff(sp.log(beta), s))
    num = sp.expand(sp.numer(dlog) * sp.sign(sp.denom(dlog).subs(s, 1)))
    check(all(c < 0 for c in sp.Poly(num, s).coeffs()) and sp.simplify(beta.subs(s, 0) - sp.Rational(3, 5)) == 0,
          "beta(0) = 3/5 and beta decreasing on s >= 0 (so 1/beta >= 5/3, r_(k+1)(s) >= (5/3) r_k(ps))")
    check(sp.simplify(sp.Rational(5, 3) / (3 + m) - sp.Rational(1, 3) - (2 - m) / (3 * (3 + m))) == 0, "(5/3)/(3+m) - 1/3 = (2-m)/(3(3+m))")
    w = rho * (1 + s / 5) ** a / ((1 + 2 * s / 5) * sp.Rational(5, 3))
    check(sp.simplify(w - (5 + s) ** a / (5 + 2 * s)) == 0, "transfer: (p(s)/s)^alpha / p'(s) = rho (1+s/5)^alpha / ((1+2s/5)(5/3))  (BK.1)")
    check(sp.simplify(5 ** (-(1 - a)) * sp.Rational(25, 9) / 5 - rho / 3) == 0, "scaling t = 5^-k x: 5^-k(1-alpha) (25/9)^k 5^-k = (rho/3)^k")
    x1 = sp.Rational(2, 5) * sp.exp(-sp.Rational(1, 10))
    U1 = 5 * x1 * sp.exp(5 * x1 * sp.exp(sp.Rational(1, 10)) / 20)
    check(sp.simplify(U1 - 2) == 0, "U(x1) = 2 at x1 = (2/5) e^(-1/10)")
    check(sp.Rational(1, 5) * sum(sp.Rational(1, 5) ** i for i in range(1, 60)) < sp.Rational(1, 20), "sum_(j<n) 5^(j-n-1) < 1/20")
    gm = sp.log(3) / sp.log(5)
    l5 = lambda v: sp.log(v) / sp.log(5)
    for J in (2, 3, 4, 5):
        thr = (J + 1) * gm / J
        check(sp.simplify(thr - l5(3 ** (J + 1)) / J) == 0, f"J = {J}: rho^J < 3 iff alpha < (J+1) gamma/J = {float(sp.N(thr)):.6f}")
    check(bool(sp.N(3 * gm / 2) > 1) and bool(sp.N(4 * gm / 3) < 1), "J = 2 threshold 3gamma/2 > 1 (always), J = 3 threshold 4gamma/3 < 1")
    RESULT["alpha3"] = float(sp.N(4 * gm / 3, 12))


# ----------------------------------------------------------------- (B) arb: ell_alpha, L_alpha
def ell_arb(al, M: int = 4000):
    from flint import arb
    x1 = arb(2) / 5 * (arb(-1) / 10).exp()
    e1 = (arb(1) / 10).exp()
    U = lambda x: 5 * x * (5 * x * e1 / 20).exp()
    f = lambda x: (2 - U(x)) / (3 * (3 + U(x)) * (3 + U(x) / 5) * (1 + x / 25) * (1 + x / 10))
    S = arb(0)
    for i in range(1, M):                      # cell [a,b]: x^(1-alpha) increasing, f decreasing and >= 0 on [0, x1]
        a, b = x1 * i / M, x1 * (i + 1) / M
        S += a ** (1 - al) * f(b) * (b - a)
    return arb(S.lower()) / 3


def certify_constants():
    from flint import arb, ctx
    ctx.prec = 128
    import gasket_fractional_corner_remainder_sharp as ec
    gm = arb(3).log() / arb(5).log()
    lg = ell_arb(gm)
    check(bool(lg.lower() > arb("0.0016681")), f"ED.2: ell_gamma >= 0.0016681 ({float(lg.lower()):.7f}); ell_alpha >= ell_gamma on (gamma,1) (x^(1-alpha) >= x^(1-gamma), x1 < 1)")
    rows = {"gamma": {"ell": lg, "L": arb((3 * lg / 2).lower())}}
    for a in ("0.7", "0.8", "0.9", "0.95"):
        al = arb(a)
        l = ell_arb(al)
        rho = arb(5) ** al / 3
        L = arb((3 * l / (2 * (1 - rho / 3))).lower())
        with contextlib.redirect_stdout(io.StringIO()):
            B = ec.constants(a)["B"]
        rows[a] = {"ell": l, "L": L, "B": B}
        print(f"   alpha = {a}: ell >= {float(l.lower()):.7f}, L = 3 ell/(2(1-rho/3)) >= {float(L.lower()):.7f}; EC.2's B <= {float(B.upper()):.4e}: "
              f"(3/rho)^n eps_n in [{float(L.lower()):.6f}, {float(B.upper()):.4e}] for every n >= 6")
        check(bool(L.lower() > 0) and bool(L.upper() < B.lower()), f"ED.2 at {a}: 0 < L_alpha < B_alpha (two-sided Theta((rho/3)^n) for eps_n)")
    check(bool(rows["0.8"]["L"].lower() > arb("0.0055320")) and bool(rows["0.95"]["L"].lower() > arb("0.0098674")),
          "ED.2: L_0.8 >= 0.0055320, L_0.95 >= 0.0098674 (0.95 > 4 gamma/3)")
    check(bool((3 * lg / 2).lower() > arb("0.0025021")), "ED.2: uniform L_alpha >= (3/2) ell_gamma >= 0.0025021 on (gamma,1)")
    RESULT["const"] = rows
    return rows


# ----------------------------------------------------------------- (C) arb: the pointwise lemmas, tested on exact r_k (attempts to break)
def lemma_spot_checks(al_str="0.8", ks=(0, 1, 3, 6, 10, 20), npts=24):
    from flint import arb, ctx
    ctx.prec = 256
    al = arb(al_str)
    rho = arb(5) ** al / 3
    x1 = arb(2) / 5 * (arb(-1) / 10).exp()
    e1 = (arb(1) / 10).exp()
    U = lambda x: 5 * x * (5 * x * e1 / 20).exp()
    f = lambda x: (2 - U(x)) / (3 * (3 + U(x)) * (3 + U(x) / 5) * (1 + x / 25) * (1 + x / 10))
    beta = lambda s: (s + 6) / ((s + 2) * (s + 5))
    e = lambda s: 2 * s / (2 + s)

    def r(k, s):
        ms = [s]
        for _ in range(k):
            ms.append(ms[-1] * (ms[-1] + 5))
        y = 1 / (3 + ms[k])
        for j in range(k - 1, -1, -1):
            y = y / (beta(ms[j]) - e(ms[j]) * y)
        return y, ms

    ok = True
    worst = None
    for k in ks:
        for i in range(1, npts):
            x = x1 * i / npts
            t = x / arb(5) ** k
            rk, ms = r(k, t)
            rk1, ms1 = r(k + 1, t)
            q = (-5 + (25 + 4 * t).sqrt()) / 2
            om = (1 + q / 5) ** al / ((1 + 2 * q / 5) * arb(5) / 3 * (beta(q) - e(q) * rk))
            c53 = (arb(5) / 3) ** k
            ok &= bool(ms1[k + 1].upper() <= U(x).lower() + 0) and bool(ms1[k + 1].upper() < 2)
            ok &= bool((rk1 - rk).lower() >= (c53 * (2 - U(x)) / (3 * (3 + U(x)))).upper())
            ok &= bool(rk.lower() >= (c53 / (3 + U(x) / 5)).upper()) and bool(rk.upper() <= (c53 / 3).lower())
            ok &= bool((om - 1).lower() >= (arb(5) / 3 * e(q) * rk).upper())
            ok &= bool(e(q).lower() >= (t / (5 * (1 + x / 25) * (1 + x / 10))).upper())
            g_true = t ** (-al) * (rk1 - rk) * (om - 1) / (arb(5) ** k * (rho / 3) ** k)    # integrand in x, divided by (rho/3)^k
            g_low = x ** (1 - al) * f(x) / 3
            ok &= bool(g_true.lower() >= g_low.upper())
            rat = float((g_true / g_low).mid())
            worst = rat if worst is None else min(worst, rat)
    print(f"   exact arb r_k, r_(k+1), omega_k at alpha = {al_str}, k in {ks}, {npts - 1} points x in (0, x1): every pointwise lemma holds; "
          f"min true/bound integrand ratio {worst:.3f}")
    check(ok, "ED.2 lemmas on exact r_k: p^(k+1)(t) <= U(x) < 2, r_(k+1)-r_k >= (5/3)^k (2-U)/(3(3+U)), (5/3)^k/(3+U/5) <= r_k <= (5/3)^k/3, "
              "omega_k(q)-1 >= (5/3) e(q) r_k(t), e(q) >= t/(5(1+x/25)(1+x/10)), integrand >= (rho/3)^k x^(1-alpha) f(x)/3")
    RESULT["spot_worst"] = worst


# ----------------------------------------------------------------- (D) arb at alpha = 0.8: 3^n R_n bracket
def certify_08(L08):
    from flint import arb, ctx
    ctx.prec = 256
    import gasket_fractional_corner_limit_mobius as eb
    import gasket_fractional_corner_remainder_sharp as ec

    def iv(lo, hi):
        lo, hi = arb(lo), arb(hi)
        mm = ((lo + hi) / 2).mid()
        return arb(mm, max((hi.upper() - mm).upper(), (mm - lo.lower()).upper()))

    with contextlib.redirect_stdout(io.StringIO()):
        c8 = ec.constants("0.8")
        n_fail = len(ec.FAILS)
        tab = ec.certify_08(c8)
        ec_ok = len(ec.FAILS) == n_fail
    check(ec_ok, "EC.4's Lam1_N, Lam2_N re-derived (EC's own checks pass)")
    al = arb("0.8")
    rho = arb(5) ** al / 3
    kap = 2 * (arb.pi() * al).sin() / arb.pi()
    A, nu = iv(*eb.A08), iv(*eb.NU08)
    Dinf = A / kap
    check(bool((1 + nu / rho).lower() > 0) and bool(nu.upper() < 0), "alpha = 0.8: 0 < 1 + nu rho^-n <= 1 for every n >= 1")
    out = {}
    for N in NS:
        Cd = iv(*eb.CD_EA[N])
        om = arb((Cd / (-nu)).lower())                    # 1 - delta_N >= C_delta_lo / (-nu)
        lo = arb((om * L08 / arb(Dinf.upper())).lower())
        hi = arb((tab[N]["Lam"] + tab[N]["Lam2"] * (3 / (5 * rho)) ** 60).upper())
        out[N] = (lo, hi)
        print(f"   N = {N:2d}: 3^n R_n in [{float(lo.lower()):.4e}, {float(hi.upper()):.4e}] for every n >= 60 (lower: every n >= N)")
    check(all(bool(out[N][0].lower() > 0) and bool(out[N][0].upper() < out[N][1].lower()) for N in NS),
          "ED.4: 0 < (1-delta_N) L_0.8 / D'_oo <= 3^n R_n <= Lam1_N + Lam2_N (3/(5 rho))^60 for every n >= 60, N = 1,2,3,6,10")
    RESULT["ed4"] = out
    return out


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem ED: eps_n >= L (rho/3)^n > 0, R_n = Theta(3^-n); EC.3's thresholds rho^J < 3 are sharp; Y_n < Y_oo eventually above 4 gamma/3. Not thrust.")
    symbolic_checks()
    rows = certify_constants()
    lemma_spot_checks("0.8")
    lemma_spot_checks("0.95", ks=(0, 2, 8, 16), npts=12)
    certify_08(rows["0.8"]["L"])
    obs, Cp_f = float_increments(0.8, (4, 6, 8, 10, 12))
    print("   floats alpha = 0.8: (3/rho)^k (Mt_(k+1)-Mt_k), (3/rho)^n eps_n = " + ", ".join(f"{n}: {d:.5f}, {e:.5f}" for n, d, e in obs) +
          f"; C' from Mt_oo = {Cp_f:.9f} (AZ (3): -4.658677573)")
    L8 = float(rows["0.8"]["L"].lower())
    l8 = float(rows["0.8"]["ell"].lower())
    check(all(e > L8 and d > l8 for _, d, e in obs) and abs(Cp_f + 4.658677573) < 1e-7,
          "ED.5 (OBSERVED): floats above ED.2's bounds (Mt-increment ~0.295 >= ell_0.8, eps ~0.741 >= L_0.8); C' = -Mt_oo rho/(rho-1) matches AZ (3)")
    RESULT["floats"] = obs
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all ED checks passed")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    sys.exit(main())
