#!/usr/bin/env python3
"""Theorem BE: the 3^-n constant c(alpha) of Theorem BD is strictly decreasing, c(1-) = 2 with a certified slope,
and c(0+) is FINITE (Conjecture BD.2 half proved, half refuted).

Setting: Assumption BA / BC (README), notation of Theorems BB, BC, BD.  d_n(a) = 2 - E_n/Y_n = <t_n> under
s^-a y_n(s) ds;  T, Y = scaling limits of 3^n t_n(5^-n tau), 3^n y_n(5^-n tau) (BD.4);
c(a) = int tau^-a T Y / int tau^-a Y.

Theorem BE (checked here; exit 1 on any failure):
  (1) monotone in alpha, every level: for 0 < a < b < 1 and every n >= 1, d_n(a) > d_n(b), so E_n/Y_n is strictly
      increasing in alpha and tan theta_n strictly increasing (|theta_n| decreasing).  Proof: for any positive
      measure m and nondecreasing f,
        Z_a Z_b (<f>_a - <f>_b) = 1/2 iint (f(x)-f(y)) (x^-a y^-b - x^-b y^-a) m(dx) m(dy) >= 0,
      the kernel has the sign of x - y; strict when f is not m-a.e. constant.  With m = y_n ds, f = t_n
      (strictly increasing Stieltjes, BD.2) this is (1); with m = Y dtau, f = T it gives
  (2) c(alpha) is strictly decreasing on (0,1)  (T nondecreasing, T(0+) = 2, T >= T_1 > 2 for tau > 10).
  (3) c(1-) = 2, with slope:  (c(a) - 2)/(1 - a) -> 3 S,  S = int_0^inf (T(tau) - 2) Y(tau) dtau / tau,
      certified 3S in [2.70632481, 2.70903118] (12903 cells).  Uses 1/3 - O(tau) <= Y <= 1/3 near 0 and T - 2 <= tau/6.
      Also directly c(a) - 2 <= B_a = O(1 - a).
  (4) c(0+) = c_0 := int T Y / int Y < inf (dominated convergence; Y <= tau/(3 psi_6(tau)) decays faster than
      any power), certified c_0 in [5.48553067, 5.50249013].  So sup_{(0,1)} c = c_0.
  Kept Failure BE.1: Conjecture BD.2's c(0+) = inf is FALSE; the 1/a blow-up of B_a is an artifact of the
      polynomial domination (1+tau/25)^-1(1+tau/125)^-1 used in BD.3.
  Conjecture BE.2 (floats only): c is convex on (0,1).  No closed form for T, Y, c, S or c_0.

Not a field solver.  Nothing here is thrust.
"""
from __future__ import annotations

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


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    x, y, a, b = sp.symbols("x y a b", positive=True)
    ker = x ** (-a) * y ** (-b) - x ** (-b) * y ** (-a)
    # kernel = (xy)^-a (y^(a-b) - x^(a-b)); for a < b and x > y it is > 0
    check(sp.simplify(ker - (x * y) ** (-a) * (y ** (a - b) - x ** (a - b))) == 0,
          "kernel x^-a y^-b - x^-b y^-a = (xy)^-a (y^(a-b) - x^(a-b)): sign of x - y when a < b")
    f = sp.Function("f")
    lhs = f(x) * ker
    sym = sp.Rational(1, 2) * (f(x) - f(y)) * ker
    # symmetrised integrand: lhs(x,y) + lhs(y,x) = 2 sym(x,y)
    check(sp.simplify(lhs + lhs.subs({x: y, y: x}, simultaneous=True) - 2 * sym) == 0,
          "Z_a Z_b(<f>_a - <f>_b) = iint f(x) ker = (1/2) iint (f(x)-f(y)) ker  (symmetrisation)")
    al = sp.Symbol("al", positive=True)
    mu = 2 - al
    J = sp.Rational(125, 4) * sp.pi / sp.sin(sp.pi * mu) * (25 ** (mu - 1) - 125 ** (mu - 1))
    J1 = sp.limit(J, al, 1)
    check(sp.simplify(J1 - sp.Rational(125, 4) * sp.log(5)) == 0 and sp.limit((1 - al) * J, al, 1) == 0,
          f"J_1 = (125/4) log 5 = {float(J1):.6f} finite, so B_a = (1-a) J_a e^0.1(3+e^0.1)/18 -> 0 as a -> 1-")
    check(sp.limit(al * J, al, 0, "+") == 3125 / sp.pi * sp.pi,
          "a J_a -> 3125 as a -> 0+: BD.3's bound B_a blows up like 1/a (the bound, not c)")
    t = sp.Symbol("t", positive=True)
    T1 = 3 * (10 + t) / (20 + t)
    check(sp.simplify(T1 - 2 - (t - 10) / (20 + t)) == 0,
          "T >= T_1 = 3(10+tau)/(20+tau) = 2 + (tau-10)/(20+tau) > 2 for tau > 10, while T(0+) = 2: T is not constant")


# ----------------------------------------------------------------- (B) finite levels (floats, consistency)
ALPHAS = ("0.1", "0.25", "0.45", "0.65", "0.8", "0.95")


def level_monotone(levels=(1, 2, 4)):
    out = {}
    for n in levels:
        vals = [bd.quad_mean_t(n, al) for al in ALPHAS]
        out[n] = vals
        ok = all(vals[i] > vals[i + 1] for i in range(len(vals) - 1))
        ok = ok and all(2 / mpf(3**n + 1) < v for v in vals)
        check(ok, f"level {n}: d_n(a) = 2 - E_n/Y_n strictly decreasing over a = {', '.join(ALPHAS)}: "
                  + ", ".join(mp.nstr(v, 8) for v in vals))
    return out


def scaled_table(n=20):
    mp.dps = 30
    vals = []
    for al in ALPHAS:
        q = bd.bb.quad_level(n, al)
        vals.append((2 - q["E"] / q["Y"]) * mpf(3) ** n)
    check(all(vals[i] > vals[i + 1] for i in range(len(vals) - 1)),
          f"n={n}: 3^n d_n(a) decreasing: " + ", ".join(mp.nstr(v, 9) for v in vals))
    sl = [(vals[i + 1] - vals[i]) / (mpf(ALPHAS[i + 1]) - mpf(ALPHAS[i])) for i in range(len(vals) - 1)]
    print("   float secant slopes of 3^20 d_20:", ", ".join(mp.nstr(v, 6) for v in sl),
          "(increasing: Conjecture BE.2, convexity, floats only)")
    return vals, sl


# ----------------------------------------------------------------- (C) certified slope at alpha = 1
def S_bracket(N=60, tau0=1e-8, tau1=3125, r1=1.05, r2=1.001, split=0.01, M=6, k=8):
    """Certified S = int_0^inf (T-2) Y dtau/tau.  (T-2)/tau is nonincreasing (T concave, T(0+) = 2) and Y is
    nonincreasing, so on [u,v] the integrand lies between its bounds at v and at u."""
    iv.dps = 25
    grid = [tau0]
    while grid[-1] < split:
        grid.append(grid[-1] * r1)
    while grid[-1] < tau1:
        grid.append(min(grid[-1] * r2, tau1))
    lo = hi = iv.mpf(0)
    prev = None
    for g in grid:
        G = iv.mpf(g)
        T, Y = bd.TY(N, G)
        eT = iv.exp(bd.epsT(N, G))
        eY = iv.exp(bd.epsY(N, G))
        up = (T * eT - 2) / G * Y * eY          # T_N <= T <= T_N e^epsT, |log Y - log Y_N| <= epsY
        low = (T - 2) / G * Y / eY
        if prev is not None:
            w = G - prev[0]
            hi += prev[1] * w
            lo += low * w
        prev = (G, up)
    hi += iv.mpf(tau0) / 18                     # head: (T-2)/tau <= 1/6, Y <= 1/3
    xx = sp.Symbol("x")
    q = xx
    for _ in range(M):
        q = sp.expand(q * (q + 5))
    q = sp.Poly(q.subs(xx, xx / sp.Integer(5) ** M), xx)
    Pk = sp.Rational(q.eval(sp.Integer(tau1))) - sum(q.coeff_monomial(xx**j) * sp.Integer(tau1) ** j for j in range(1, k))
    Pki = iv.mpf(int(Pk.p)) / iv.mpf(int(Pk.q))
    t1 = iv.mpf(tau1)
    hi += t1**k / (18 * Pki) * t1 ** (2 - k) / (k - 2)   # tail: Y <= tau^(1-k) tau1^k/(3 P_k), (T-2)/tau <= 1/6
    return iv.mpf([lo.a, hi.b]), len(grid)


# ----------------------------------------------------------------- main
def main(r2=1.001) -> int:
    t0 = time.time()
    print("Theorem BE: c(alpha) strictly decreasing, c(1-) = 2, c(0+) finite (Conjecture BD.2). Not thrust.")
    symbolic_checks()
    mp.dps = 30
    level_monotone()
    vals, sl = scaled_table()

    t1 = time.time()
    S, ncell = S_bracket(r2=r2)
    S3 = 3 * S
    print(f"   S in [{float(S.a):.10f}, {float(S.b):.10f}], 3S in [{float(S3.a):.8f}, {float(S3.b):.8f}] "
          f"({ncell} cells, {time.time()-t1:.1f}s)")
    check(S3.a > 0 and S3.b - S3.a < 1e-2, "slope at alpha = 1: (c(a)-2)/(1-a) -> 3S, certified enclosure")
    # consistency with c(0.8) of Theorem BD: convexity is NOT proved, so only print the comparison
    print(f"   compare: BD's c(0.8) secant (c(0.8)-2)/0.2 in [{(2.57126273-2)/0.2:.5f}, {(2.57352573-2)/0.2:.5f}]")

    t1 = time.time()
    c0, ncell = bd.c_bracket("0", r2=r2)
    print(f"   c(0+) = int T Y / int Y in [{float(c0.a):.8f}, {float(c0.b):.8f}] ({ncell} cells, {time.time()-t1:.1f}s)")
    check(c0.b < 6 and c0.a > 3.72914181 and c0.a > vals[0],
          "Kept Failure BE.1: c(0+) is finite (certified), above c(0.45) and above the float 3^20 d_20(0.1)")
    cq, _ = bd.c_bracket("0.25", N=40, tau0=1e-12, tau1=625, r2=1.01, M=6, k=8)
    print(f"   coarse c(0.25) in [{float(cq.a):.6f}, {float(cq.b):.6f}]")
    check(cq.b < c0.a and cq.a > 3.72914181, "certified ordering c(0.45) < c(0.25) < c(0+), as (2) requires")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BE checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
