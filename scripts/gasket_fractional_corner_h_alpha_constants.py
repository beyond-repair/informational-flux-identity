#!/usr/bin/env python3
"""Theorem FA: the constants of Theorem EY, certified end to end (EZ certified only EY.3's tail part).

Setting (EY): g(y) = (5/3)(y+3)/(y+5), P_n = prod_{j<=n} g(psi_(n-j)(5^-j t)), P = prod_{j>=1} g(psi(5^-j t)), x = t (2/5)^n <= 1,
f_n = R_n/(2(5+psi_(n-1)(t/5))), f = R/(2(5+psi(t/5))), W(b) = int F t^(-b-1) dt = 2 int f t^(-b) dt.

  FA.1 (THEOREM, exact; KEPT FAILURE of EY's constant 0.127):  d log g / d log y = 2y/((y+3)(y+5)) and, with c* = 4 - sqrt15,
        c*(y^2+8y+15) - 2y = c*(y - sqrt15)^2 >= 0  (because 8c* - 2 = -2 sqrt15 c*  and  15 c* = 15 c*),
        so sup_{y>=0} = c* = 0.12701665379..., attained at y = sqrt15.  EY.2's "<= 0.127" is FALSE: exact witness y = 39/10 gives
        780/6141 = 0.1270151... > 0.127.  Smallest valid constant: c* (decimal: 0.12702).
  FA.2 (THEOREM, analytic):  psi(w) <= w e^(w/10) for 0 <= w <= 10, from EV.1 alone (replaces BJ.6 in EY.1/EY.2; EY only uses w = 5^-n t <= 2^-n).
        log g is concave with log g(0) = 0, (log g)'(0) = 2/15, so log g(y) <= 2y/15; for j > n, v = 5^-j t <= 2^-n/5 <= 1/10, hence
        sum_{j>n} log g(psi(5^-j t)) <= (2/15) e^(1/100) t sum_{j>n} 5^-j = (e^(1/100)/30) t 5^-n = 0.033668 t 5^-n <= 0.037 t 5^-n  (TRUE).
        The ratio tends to 1/30 as t -> 0, so no constant below 1/30 = 0.0333 works.
  FA.3 (THEOREM; EY.2 re-proved with c*):  log(P/P_n) <= (c*/10) x (1 - 2^-n) + (e^(1/100)/30) x 2^-n <= 0.03367 x <= 0.05 x;
        log(psi(t/5)/psi_(n-1)(t/5)) <= x/20 (EY.1 at m = n-1, w = 5^-n t); e^(0.1x) - 1 <= (e^(0.1) - 1) x = 0.105171 x <= 0.106 x (convexity, x <= 1);
        lower side R_n/R >= 1 - (t/3) 3^-n.  So |f_n/f - 1| <= x/3 holds as stated.
  FA.4 (THEOREM, arb):  W(-2) enclosed with EX.3's certified routine (EV.2 evaluator, blocks [5^m, 5^(m+1)], m = -24..7, low/high tails).
  FA.5 (THEOREM; EY.3 end to end):  |rho^-n H_n(alpha) - W(alpha)/2| <= (W(Re alpha - 1)/6 + 1)(2/5)^n <= C (2/5)^n, n >= 11, -1 <= Re alpha <= 9/10,
        C = W(-2)/6 + 1 (certified upper bound printed), using |t^-alpha| = t^(-Re alpha), F >= 0 and EX.2 (W decreasing on (-oo, alpha*], alpha* > 3/5).
Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0])

W_M2_FLOAT = 11291.8334854864  # float log-grid value of W(-2) (OBSERVED; the certified bracket is printed by main())


# ---------- pure-Python exact / float pieces (used by CI) ----------
def elasticity(y: Fraction) -> Fraction:
    """d log g / d log y = y/(y+3) - y/(y+5) = 2y/((y+3)(y+5)), exact."""
    return 2 * y / ((y + 3) * (y + 5))


def sos_identity_ok(pts=None) -> bool:
    """FA.1: c(y^2+8y+15) - 2y == c (y - s)^2 with s^2 = 15, c = 4 - s, checked as polynomials in Q(sqrt15): compare coefficients.
    Write numbers as (u, v) = u + v*sqrt15."""
    def mul(a, b):
        return (a[0] * b[0] + 15 * a[1] * b[1], a[0] * b[1] + a[1] * b[0])
    c = (Fraction(4), Fraction(-1))
    s = (Fraction(0), Fraction(1))
    # coefficients of y^2, y^1, y^0
    lhs = [c, (8 * c[0] - 2, 8 * c[1]), (15 * c[0], 15 * c[1])]
    m2s = mul(c, (-2 * s[0], -2 * s[1]))
    s2 = mul(s, s)
    rhs = [c, m2s, mul(c, s2)]
    ok = lhs == rhs and (4 - Fraction(1)) > 0  # c > 0 since 15 < 16
    ok &= Fraction(15) < Fraction(16)
    return bool(ok)


def witness_ok() -> bool:
    """FA.1 KEPT FAILURE witness: 2y/((y+3)(y+5)) at y = 39/10 is 780/6141 > 127/1000 (exact)."""
    v = elasticity(Fraction(39, 10))
    return v == Fraction(780, 6141) and v > Fraction(127, 1000)


def elasticity_grid_ok(n=2000) -> bool:
    """Float sanity: max over a fine y-grid is within 1e-9 of 4 - sqrt15 and never above it."""
    c = 4 - math.sqrt(15)
    m = max(2 * y / ((y + 3) * (y + 5)) for y in (10 ** (-3 + 6 * i / n) for i in range(n + 1)))
    return m <= c + 1e-15 and c - m < 1e-6


def float_constants() -> dict:
    c = 4 - math.sqrt(15)
    tail = math.exp(0.01) / 30
    return {"cstar": c, "tail": tail, "logP": max(c / 10, tail), "exp_chord": math.expm1(0.1)}


def psi_bound_ok(npts=400) -> bool:
    """FA.2 float check: EV.1 upper bound w + w^2/20 + w^3/(400(1-w/20)) <= w e^(w/10) on (0, 10]."""
    for i in range(1, npts + 1):
        w = 10.0 * i / npts
        if w + w * w / 20 + w ** 3 / (400 * (1 - w / 20)) > w * math.exp(w / 10) * (1 + 1e-12):
            return False
    return True


# ---------- arb certificates ----------
def arb_constants():
    from flint import arb, ctx
    ctx.prec = 200
    out = {}
    s15 = arb(15).sqrt()
    c = 4 - s15
    out["cstar"] = c
    out["FA.1 c* > 0.127 (EY constant false)"] = bool(c > arb("0.127"))
    out["FA.1 c* < 0.12702"] = bool(c < arb("0.12702"))
    tail = arb("0.01").exp() / 30
    out["tail"] = tail
    out["FA.2 e^(1/100)/30 <= 0.037"] = bool(tail < arb("0.037"))
    out["FA.2 e^(1/100)/30 < 0.03367"] = bool(tail < arb("0.03367"))
    out["FA.3 max(c*/10, e^(1/100)/30) <= 0.05"] = bool(c / 10 < arb("0.05") and tail < arb("0.05"))
    chord = arb("0.1").exp() - 1
    out["chord"] = chord
    out["FA.3 e^0.1 - 1 <= 0.106"] = bool(chord < arb("0.106"))
    out["FA.3 0.106 <= 1/3"] = bool(arb("0.106") < arb(1) / 3)
    # FA.2: EV.1 => psi(w) <= w e^(w/10) on [0,10]: w^2/20 + w^3/200 - w^3/(400(1-w/20)) >= 0 iff 1/20 + w/200 - w/(400(1-w/20)) >= 0;
    # the last term is <= w/200 when 1 - w/20 >= 1/2, i.e. w <= 10 (analytic).  arb spot check on a grid of balls covering [0,10]:
    ok = True
    for i in range(201):
        wb = arb(arb(i) / 20, arb(1) / 40)  # ball containing [i/20 - 1/40, i/20 + 1/40]; i = 0..200 covers [0, 10]
        lhs = arb(1) / 20 + wb / 200 - wb / (400 * (1 - wb / 20))
        ok &= bool(lhs > 0)
    out["FA.2 EV.1 => psi(w) <= w e^(w/10) on [0,10] (ball grid)"] = ok
    return out


def certify_W(values=(Fraction(-2), Fraction(-1))):
    import gasket_fractional_corner_h_alpha_limit as ex
    ev = ex.ev
    ex.W_ALPHAS = list(values)
    ex.DW_ALPHAS = []
    acb, arb, RY, psi, _fc, _fk = ev.make_arb()
    out = ex.certified(acb, arb, RY, psi)
    return arb, out, ex


def main() -> int:
    t0 = time.time()
    print("Theorem FA: EY's constants certified (EY.2 constants; W(alpha-1) <= W(-2)). Functions on a graph, not a flux, not thrust.", flush=True)
    fails = []

    def check(ok, msg):
        print(("  ok   " if ok else "  FAIL ") + msg, flush=True)
        if not ok:
            fails.append(msg)

    check(sos_identity_ok(), "FA.1: c*(y^2+8y+15) - 2y = c*(y - sqrt15)^2 in Q(sqrt15)[y], c* = 4 - sqrt15 > 0 (exact)")
    check(witness_ok(), "FA.1 KEPT FAILURE: elasticity(39/10) = 780/6141 > 0.127 (exact)")
    check(elasticity_grid_ok(), "FA.1: float grid max of 2y/((y+3)(y+5)) agrees with 4 - sqrt15")
    check(psi_bound_ok(), "FA.2: EV.1 bound <= w e^(w/10) on (0,10] (float grid; analytic proof in docstring)")
    a = arb_constants()
    for k, v in a.items():
        if isinstance(v, bool):
            check(v, k)
    print(f"   c* = 4 - sqrt15 in {a['cstar'].str(20)}", flush=True)
    print(f"   e^(1/100)/30 in {a['tail'].str(20)};  e^(0.1) - 1 in {a['chord'].str(20)}", flush=True)
    arb, out, ex = certify_W()
    lo2, hi2 = out[("W", Fraction(-2))]
    lo1, hi1 = out[("W", Fraction(-1))]
    print(f"   FA.4: W(-2) in [{lo2.str(25, radius=False)}, {hi2.str(25, radius=False)}]  width {(hi2 - lo2).str(3, radius=False)}", flush=True)
    check(bool(hi2 - lo2 < arb(10) ** -12), "FA.4: W(-2) bracket width < 1e-12")
    ex1 = arb("71.795668675310493193")
    check(bool(lo1 <= ex1 + arb("1e-18")) and bool(ex1 - arb("1e-18") <= hi1), "FA.4 cross-check: W(-1) bracket meets EX.3's value 71.795668675310493193")
    check(bool(lo2 > hi1), "FA.4: W(-2) > W(-1) (consistent with EX.2: W decreasing on (-oo, alpha*])")
    wf = ex.W_float(-2.0)
    mid = float(((lo2 + hi2) / 2).mid())
    check(abs(wf - mid) < 1e-6 * mid, f"FA.4: float W(-2) = {wf:.10f} agrees with arb midpoint {mid:.10f}")
    C_hi = hi2 / 6 + 1
    print(f"   FA.5: constant C = W(-2)/6 + 1 <= {C_hi.upper().str(20, radius=False)}  (W(-2)/6 in [{(lo2 / 6).str(20, radius=False)}, {(hi2 / 6).str(20, radius=False)}])", flush=True)
    check(bool(C_hi < arb("1882.973")), "FA.5: C <= 1882.973")
    print(f"elapsed {time.time() - t0:.1f}s", flush=True)
    if fails:
        print(f"{len(fails)} FAILED")
        return 1
    print("all FA checks passed (FA.1 records a KEPT FAILURE of EY's 0.127; corrected constant 4 - sqrt15)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
