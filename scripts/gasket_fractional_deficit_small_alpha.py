#!/usr/bin/env python3
"""Theorem BP: the corner deficit delta_N^(n) increases in n for EVERY alpha in (0,1), including the range alpha < a* where
Kept Failure BO.1 showed that BK's one-level weight omega_k is not monotone (so BN.3's Chebyshev route is closed there).
The replacement for Chebyshev is a single-crossing covariance argument: omega_k need not be monotone, it only has to cross
its own mean once, and the mean is pinned from below by AQ's monotone resolvent.

Setting (GASKET.md, Theorems AQ, BK, BL, BN, BO): r_k(s) = (1/2) v^T (L_k + s)^-1 v, p(y) = y(y+5), e(y) = 2y/(2+y),
B(y) = 1 - y(y+4)/(3(y+2)(y+5)), C(y) = (5/3)(1+2y/5) e(y), d_k = B - C r_k(p) > 0, omega_k = (1+y/5)^a / d_k, rho = 5^a/3,
D_n = int s^-a r_n ds, D'_n = rho^-n D_n, mu_k = BK.1's level-k probability measure (positive density), G_k = 1 - eta_N(p^(k-N+1)(y))
(strictly increasing in y, BL.1), delta_N^(k) = <G_k>_(mu_k), delta_N^(k+1) - delta_N^(k) = Cov_(mu_k)(omega_k, G_k)/<omega_k> (BK.1).

Theorem BP (checked here; exit 1 on any failure):
  (1) [mean floor] <omega_k>_(mu_k) = D'_(k+1)/D'_k = rho^-1 D_(k+1)/D_k >= rho^-1 = 3 * 5^-a, because r_(k+1) >= r_k pointwise
      (AQ (3), used as stated; exact rationals re-check r_(k+1) > r_k at k = 0..3 on sample points).
  (2) [tail monotone, every k] BN.2's majorant holds for every a >= 0: d_k' - a d_k/(5+y) <= -P_a/(3(y+2)^2(y+5)^2) with
      P_a = y^2 (2a y + 17a - 1) + 46a y + 40a + 10 (sympy).  So with Y_a = (1-17a)/(2a), P_a > 0 on [Y_a, inf) and omega_k is
      strictly increasing on [Y_a, inf) at every level k, for every a in (0, 1/17].
  (3) [head below the mean, every k] r_k(p) <= 1/p (BN.1) gives d_k >= B - C/p = (y+4)(2y+5)/(3(y+2)(y+5)), and
      (8/5)(y+4)(2y+5) - 3(y+2)(y+5) = (y^2 - y + 10)/5 > 0, so omega_k < (8/5)(1+y/5)^a.  On [0, Y_a]:
      (1+y/5)^a <= 5^-a (5+Y_a)^a, 5 + Y_a = 1/(2a) - 7/2 < 1/(2a), and (2a)^-a <= e^(1/(2e)) < 15/8 (calculus; arb), so
      omega_k < (8/5)(15/8) 5^-a = 3 * 5^-a = rho^-1 <= <omega_k> on [0, Y_a].
  (4) [single crossing] Hence {omega_k >= <omega_k>} = [y_c, inf) for one y_c > Y_a, and with m = <omega_k>,
      Cov(omega_k, G_k) = int (omega_k - m)(G_k - G_k(y_c)) dmu_k > 0 (integrand >= 0, > 0 off y_c).  So
      delta_N^(k+1) > delta_N^(k) for every a in (0, 1/17], N >= 1, k >= N; with BN.3 (a >= 1/17): for EVERY a in (0,1).
      This answers BO's open question: the failure of monotone omega_k below a* (BO.1) does not propagate to the deficit.
  (5) [derived, BG as stated] for a in (0, gamma], u_n -> dbar at every non-corner vertex (BG), so delta_N^(n) increases
      strictly to 1.
Floats (consistency only): omega_k at a = 0.001, 0.005, 0.01, 0.0167 (k = 3, 6, 10): non-monotone (BO.1 dip) yet max on
[0, Y_a] below rho^-1; delta_N^(n) at a = 0.01, 0.001 (N = 1, 2; n = N..N+5) increasing to 1.
Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import sympy as sp
from flint import arb
from mpmath import mp, mpf, quad

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_maxwell_limit as ml  # noqa: E402
import gasket_fractional_monotone_weight as bn  # noqa: E402

FAILS: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    y, a = sp.symbols("y a", positive=True)
    p = y * (y + 5)
    B = 1 - y * (y + 4) / (3 * (y + 2) * (y + 5))
    e = 2 * y / (2 + y)
    C = sp.Rational(5, 3) * (1 + 2 * y / 5) * e
    Pa = 2 * a * y**3 + (17 * a - 1) * y**2 + 46 * a * y + 40 * a + 10
    check(sp.expand(Pa - (y**2 * (2 * a * y + 17 * a - 1) + 46 * a * y + 40 * a + 10)) == 0,
          "P_a = y^2 (2a y + 17a - 1) + 46a y + 40a + 10: for y >= Y_a = (1-17a)/(2a) every term is >= 0 and 40a+10 > 0, so P_a > 0")
    BC = sp.factor(sp.together(B - C / p))
    check(sp.simplify(BC - (y + 4) * (2 * y + 5) / (3 * (y + 2) * (y + 5))) == 0,
          "d_k >= B - C/p = (y+4)(2y+5)/(3(y+2)(y+5)) > 0 (r_k(p) <= 1/p, C >= 0)")
    q = sp.expand(sp.Rational(8, 5) * (y + 4) * (2 * y + 5) - 3 * (y + 2) * (y + 5))
    check(sp.expand(q - (y**2 - y + 10) / 5) == 0 and sp.discriminant(y**2 - y + 10, y) < 0,
          "(8/5)(y+4)(2y+5) - 3(y+2)(y+5) = (y^2-y+10)/5 > 0 (discriminant -39): 1/(B - C/p) < 8/5, so omega_k < (8/5)(1+y/5)^a")
    Ya = (1 - 17 * a) / (2 * a)
    check(sp.simplify(5 + Ya - (1 / (2 * a) - sp.Rational(7, 2))) == 0 and Ya.subs(a, sp.Rational(1, 17)) == 0,
          "5 + Y_a = 1/(2a) - 7/2 < 1/(2a), and Y_a >= 0 exactly for a <= 1/17 ((1+y/5)^a <= 5^-a (5+Y_a)^a on [0,Y_a])")
    f = -a * sp.log(2 * a)
    crit = sp.solve(sp.diff(f, a), a)
    check(crit == [1 / (2 * sp.E)] and sp.simplify(sp.diff(f, a, 2) + 1 / a) == 0
          and sp.simplify(f.subs(a, 1 / (2 * sp.E)) - 1 / (2 * sp.E)) == 0,
          "-a log(2a) is concave (f'' = -1/a) with unique critical point a = 1/(2e), value 1/(2e): (2a)^-a <= e^(1/(2e)) for a > 0")
    ub = (arb(1) / (2 * arb.const_e())).exp()
    check(ub.upper() < arb(15) / 8, f"arb: e^(1/(2e)) <= {float(ub.upper()):.6f} < 15/8, so (8/5)(5+Y_a)^a 5^-a < 3 * 5^-a = rho^-1")
    # single-crossing lemma on finite atoms (the identity used; sign is pointwise)
    m1, m2, m3, w1, w2, w3, g1, g2, g3, g0 = sp.symbols("m1 m2 m3 w1 w2 w3 g1 g2 g3 g0")
    M = m1 + m2 + m3
    mw = (m1 * w1 + m2 * w2 + m3 * w3) / M
    cov = (m1 * w1 * g1 + m2 * w2 * g2 + m3 * w3 * g3) / M - mw * (m1 * g1 + m2 * g2 + m3 * g3) / M
    alt = (m1 * (w1 - mw) * (g1 - g0) + m2 * (w2 - mw) * (g2 - g0) + m3 * (w3 - mw) * (g3 - g0)) / M
    check(sp.simplify(cov - alt) == 0,
          "Cov(omega,G) = <(omega - <omega>)(G - g0)> for every constant g0 (atoms); with g0 = G(y_c) the integrand is >= 0 "
          "when omega - <omega> changes sign once (- to +) at y_c and G is increasing")


# ----------------------------------------------------------------- (B) exact rationals
def exact_checks():
    from fractions import Fraction as Fr
    ok = True
    for k in range(0, 4):
        for sv in (Fr(1, 1000), Fr(1, 10), Fr(1), Fr(7, 2), Fr(30)):
            ok = ok and bn.r_rec(k + 1, sv) > bn.r_rec(k, sv)
    check(ok, "exact rationals: r_(k+1)(s) > r_k(s) at s = 1/1000, 1/10, 1, 7/2, 30 for k = 0..3 (consistency with AQ (3), "
              "r_n increasing in n, which gives D_(k+1) > D_k and <omega_k> > rho^-1)")
    # BK.1's weight at exact rationals: omega_k >= 1 and d_k >= B - C/p
    ok2 = True
    for k in range(0, 4):
        for yv in (Fr(1, 100), Fr(1), Fr(13), Fr(14), Fr(200)):
            pv = yv * (yv + 5)
            Bv = 1 - yv * (yv + 4) / (3 * (yv + 2) * (yv + 5))
            Cv = Fr(5, 3) * (1 + 2 * yv / 5) * 2 * yv / (2 + yv)
            dk = Bv - Cv * bn.r_rec(k, pv)
            ok2 = ok2 and Bv - Cv / pv <= dk <= Bv <= 1 and dk * Fr(8, 5) > 1
    check(ok2, "exact rationals k = 0..3, y = 1/100, 1, 13, 14, 200: B - C/p <= d_k <= B <= 1 and d_k > 5/8 (omega_k < (8/5)(1+y/5)^a)")


# ----------------------------------------------------------------- (C) floats
def float_weight(alphas=("0.001", "0.005", "0.01", "0.0167"), ks=(3, 6, 10)):
    mp.dps = 30
    ok = True
    rows = []
    for al in alphas:
        a = mpf(al)
        Ya = (1 - 17 * a) / (2 * a)
        rinv = 3 / mpf(5) ** a
        for k in ks:
            def om(yv):
                pv = yv * (yv + 5)
                Bv = 1 - yv * (yv + 4) / (3 * (yv + 2) * (yv + 5))
                Cv = mpf(5) / 3 * (1 + 2 * yv / 5) * 2 * yv / (2 + yv)
                return (1 + yv / 5) ** a / (Bv - Cv * ml.rec_all(k, pv)[0])
            head = [om(Ya * j / 400) for j in range(1, 401)]
            tail = [om(Ya * mpf(10) ** (mpf(j) / 50)) for j in range(0, 301)]
            nonmono = any(head[i] > head[i + 1] for i in range(len(head) - 1))
            ok = ok and max(head) < rinv and all(tail[i] < tail[i + 1] for i in range(len(tail) - 1)) and nonmono
            rows.append((al, k, float(max(head)), float(rinv), float(Ya)))
    return ok, rows


def _etaN(z, N):
    zs = [z]
    for _ in range(N):
        zs.append(zs[-1] * (zs[-1] + 5))
    e = mpf(-1) / 2
    for k in range(N - 1, -1, -1):
        zz = zs[k]
        D = (zz + 2) * (zz + 5)
        e = (zz + 4) / D + (zz + 6) / D * e
    return e


def deficit(al, n, N, X=None):
    """Floats: delta_N^(n) = D_n^-1 int s^-a r_n (1 - eta_N(p^(n-N) s)) ds (BI's formula), with the s > X = 1e12 tail done
    analytically (r_n = 1/s - 2/s^2 + 6/s^3 + ..., 1 - eta_N = 1 + O(1/s)); the heavy s^(-1-a) tail is what stalls BI's quadrature."""
    mp.dps = 20
    a = mpf(al)
    X = mpf(10) ** 12 if X is None else X
    pts = [0] + [mpf(5) ** (-n) * mpf(10) ** j for j in range(-6, 40, 2) if mpf(5) ** (-n) * mpf(10) ** j < X] + [X]

    def pk(s, k):
        for _ in range(k):
            s = s * (s + 5)
        return s

    f = lambda s: s ** (-a) * ml.rec_all(n, s)[0]
    tail = X ** (-a) / a - 2 * X ** (-1 - a) / (1 + a) + 6 * X ** (-2 - a) / (2 + a)
    D = quad(f, pts) + tail
    num = quad(lambda s: f(s) * (1 - _etaN(pk(s, n - N), N)), pts) + tail
    return num / D


def main() -> int:
    t0 = time.time()
    print("Theorem BP: delta_N^(n) increases in n for every alpha in (0,1), also below a* (single-crossing covariance). Not thrust.")
    symbolic_checks()
    exact_checks()
    bn.FAILS.clear()
    bn.symbolic_checks()
    bn.exact_checks()
    check(not bn.FAILS, "Theorem BN's sympy and exact checks re-run (BN.1 Stieltjes r_k <= 1/s, 0 <= -r_k' <= r_k/s; BN.2 identity and "
                        "majorant -P_a/(3(y+2)^2(y+5)^2), valid for every a >= 0)")
    okw, rows = float_weight()
    for al, k, mx, ri, Ya in rows:
        print(f"   floats alpha={al} k={k}: Y_a={Ya:.2f}, max omega_k on [0,Y_a] = {mx:.6f} < rho^-1 = {ri:.6f}")
    check(okw, "floats: omega_k is NOT monotone on [0,Y_a] (BO.1 dip) at alpha = 0.001, 0.005, 0.01, 0.0167, k = 3, 6, 10, yet stays "
               "below rho^-1 there and increases on [Y_a, 1e6 Y_a] (consistency with BP.2-BP.3)")
    okd = True
    for al in ("0.01", "0.001"):
        for N in (1, 2):
            vals = [deficit(al, n, N) for n in range(N, N + 6)]
            okd = okd and all(vals[i] < vals[i + 1] for i in range(5)) and vals[-1] < 1
            print(f"   floats alpha={al} N={N}: delta_N^(n), n={N}..{N + 5} = " + ", ".join(f"{float(v):.10f}" for v in vals))
    ref = deficit("0.8", 8, 1)
    print(f"   floats control alpha=0.8: delta_1^(8) = {float(ref):.4f} (BN's run: 0.9487)")
    check(okd and abs(float(ref) - 0.9487) < 1e-3,
          "floats: delta_N^(n) strictly increasing in n and < 1 at alpha = 0.01, 0.001 (N = 1, 2, six levels each), and the "
          "tail-split quadrature reproduces BN's delta_1^(8)(0.8) = 0.9487 (consistency only)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BP checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
