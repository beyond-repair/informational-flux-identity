#!/usr/bin/env python3
"""Theorem ET: the refinement limit of H_n(0) (ES.2's extension of H_n to alpha = 0).  ES kept "floor only, no upper bound/limit/rate"
and computed H_n(0) = int_0^oo h_n = -sum_j h_j ln w_j for n <= 8 only, with no statement about n -> oo.

Setting (as in EQ, ER, ES): p(s) = s(s+5), q_k = p^k + 5, r_n = B_n/A_n = (1/2) v^T (L_n + s)^-1 v, h_n = Q_n/(2 A_n) = r_n/(2 q_(n-1))
(B_n = Q_n q_(n-1), EQ.1), H_n(a) = int_0^oo s^-a h_n.  gamma = log3/log5, rho = 5^a/3.  BG.3's rescaled resolvent
R_n(t) = (3/5)^n r_n(5^-n t) and psi_k(t) = p^k(5^-k t); R = lim R_n, psi = lim psi_k (BG.3, BD.4).
Assumptions used as stated: BG.3 (pointwise R_n -> R; psi_k -> psi), BI.4 (R_n <= min(1/3, G t^(gamma-1)), G = 0.6031, and
R_n >= 1/(3+t), n >= 1), BJ.6 (t <= psi(t) <= t e^(t/10)), EQ.1/ES.2 (h_n = r_n/(2 q_(n-1)), H_n(0) = int_0^oo h_n).

  ET.1 (THEOREM, every n >= 1; exact scaling).  s = 5^-n t gives
          3^n H_n(0) = int_0^oo R_n(t) / (2 (5 + psi_(n-1)(t/5))) dt,
        and for a in (0,1), rho^-n H_n(a) = int_0^oo t^-a R_n(t)/(2(5 + psi_(n-1)(t/5))) dt (BG.3's m = 1 numerator, recalled).
  ET.2 (THEOREM).  3^n H_n(0) -> K := W(0)/2 = int_0^oo R(t)/(2(5 + psi(t/5))) dt, BG.3's W at alpha = 0.  BG.3's dominator
        t^-a (1/3)/(4 + t/5) is NOT integrable at a = 0 (log divergence); BI.4's ceiling G t^(gamma-1) supplies the dominator
        min(1/3, G t^(gamma-1)) * 5/(2(25 + t)) = O(t^(gamma-2)).  Hence H_(n+1)(0)/H_n(0) -> 1/3 and H_n(0) = (K + o(1)) 3^-n.
  ET.3 (THEOREM, uniform in n >= 1).  L <= 3^n H_n(0) <= U and L <= K <= U, with
          U = (5/6) ln((25+t0)/25) + (5G/2) int_t0^oo t^(gamma-1)/(25+t) dt,  t0 = (3G)^(1/(1-gamma)),
          L = int_0^T dt / (2 (3+t) (5 + (t/5) e^(t/50))),  T = 4000   (psi_(n-1) <= psi, since p(y) >= 5y).
        Certified with python-flint arb.
  ET.4 (checks / attempts to break).  (a) the pointwise recursion r_(k+1)(s) = Phi_s(r_k(p(s))) equals EQ.1's B_n/A_n (exact
        rationals, n <= 6); (b) mpmath quadrature of ET.1's integral reproduces ES.4's certified H_n(0), n = 1..8, to 1e-9;
        (c) 3^n H_n(0) lies in [L, U] for n = 1..NOBS; (d) OBSERVED only: the values for n <= NOBS and their increments;
        (e) the imported BI.4 ceiling/floor and BJ.6 psi bracket hold on a grid, n = 1..20.
  ET.5 (CONJECTURE, OBSERVED).  K = W(0)/2 = ln(5)/2 = 0.80471895621705...: 3^n H_n(0) = K - c 3^-n + ..., c ~ 2.67182; Richardson
        (3 v_n - v_(n-1))/2 agrees with ln5/2 to ~1e-22 at n = 31 (separate probe) and to 1e-15 at n = NOBS here.  Not proved.

Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_pole_structure as eq  # noqa: E402

FAILS: list[str] = []
G = Fraction(6031, 10000)  # BI.4's certified sup of s^(1-gamma) r_oo(s)
T_LOW = 4000
NOBS = 24
# ES.4 (b) certified values of H_n(0), n = 1..8 (10 digits as printed in GASKET.md)
ES_H0 = ["0.1579836359", "0.0685090601", "0.0266745120", "0.0095409589", "0.0032661609", "0.0010987844", "0.0003673918",
         "0.0001225894"]


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def phi(s, x):
    """AQ's Moebius step: r_(k+1)(s) = Phi_s(r_k(p(s)))."""
    return (s + 2) * (s + 5) * x / ((s + 6) - 2 * s * (s + 5) * x)


def r_rec(n: int, s):
    """r_n(s) by the decimation recursion (works for Fraction, float or mpf)."""
    t = [s]
    for _ in range(n):
        t.append(t[-1] * (t[-1] + 5))
    x = 1 / (3 + t[n])
    for k in range(n - 1, -1, -1):
        x = phi(t[k], x)
    return x


def h_rec(n: int, s):
    """h_n(s) = r_n(s) / (2 (p^(n-1)(s) + 5))."""
    y = s
    for _ in range(n - 1):
        y = y * (y + 5)
    return r_rec(n, s) / (2 * (y + 5))


def rescaled_integrand(n: int, t):
    """ET.1's integrand R_n(t)/(2(5 + psi_(n-1)(t/5))) = 3^n 5^-n h_n(5^-n t)."""
    import mpmath as mp
    s = t / mp.mpf(5) ** n
    return mp.mpf(3) ** n * h_rec(n, s) / mp.mpf(5) ** n


def exact_check(n: int) -> None:
    A, B = eq.eq1_pair(n)
    ok = True
    for s in (Fraction(0), Fraction(1, 7), Fraction(2), Fraction(13, 3)):
        lhs = r_rec(n, s)
        rhs = Fraction(sum(c * s ** i for i, c in enumerate(B)), sum(c * s ** i for i, c in enumerate(A)))
        ok = ok and lhs == rhs
        # ET.1's scaling identity, exactly: 3^n 5^-n h_n(5^-n t) = R_n(t) / (2 (5 + psi_(n-1)(t/5)))
        t = s
        Rn = Fraction(3, 5) ** n * r_rec(n, t / Fraction(5) ** n)
        psi_n1 = t / Fraction(5) ** n
        for _ in range(n - 1):
            psi_n1 = psi_n1 * (psi_n1 + 5)
        ok = ok and Fraction(3, 5) ** n * h_rec(n, t / Fraction(5) ** n) == Rn / (2 * (5 + psi_n1))
    Q = eq.b_product(n - 1) if n >= 2 else [1]
    s = Fraction(3, 11)
    ok = ok and h_rec(n, s) == Fraction(sum(c * s ** i for i, c in enumerate(Q)), 2 * sum(c * s ** i for i, c in enumerate(A)))
    check(ok, f"ET.4(a) n = {n}: recursion r_n = B_n/A_n and h_n = Q_n/(2A_n) at exact rational points")


def scaled_H0(n: int):
    import mpmath as mp
    pts = [0] + [mp.mpf(5) ** k for k in range(-3, 9)] + [mp.inf]
    return mp.quad(lambda t: rescaled_integrand(n, t), pts)


def certified_bounds():
    from flint import arb, acb, ctx
    ctx.prec = 200
    g = arb(3).log() / arb(5).log()
    Ga = arb(G.numerator) / G.denominator
    t0 = (3 * Ga) ** (1 / (1 - g))
    # int_t0^oo t^(g-1)/(25+t) = 25^(g-1) pi/sin(pi g) - sum_k (-1)^k t0^(g+k)/((g+k) 25^(k+1)), alternating, t0 < 25
    full = arb(25) ** (g - 1) * arb.pi() / (arb.pi() * g).sin()
    head = arb(0)
    K = 400
    for k in range(K):
        head += (-1) ** k * t0 ** (g + k) / ((g + k) * arb(25) ** (k + 1))
    tailterm = t0 ** (g + K) / ((g + K) * arb(25) ** (K + 1))
    head = head + arb(0, tailterm.upper())
    U = arb(5) / 6 * ((25 + t0) / 25).log() + 5 * Ga / 2 * (full - head)
    f = lambda x, _: 1 / (2 * (3 + x) * (5 + (x / 5) * (x / 50).exp()))  # noqa: E731
    L = acb.integral(f, 0, T_LOW).real
    return L, U, t0


def main(nobs: int = NOBS) -> int:
    import mpmath as mp
    t_start = time.time()
    mp.mp.dps = 30
    print("Theorem ET: 3^n H_n(0) -> W(0)/2; certified uniform bracket. Not thrust.", flush=True)
    for n in range(1, 7):
        exact_check(n)
    L, U, t0 = certified_bounds()
    print(f"   t0 = {t0.str(12)}; L = {L.str(15)}; U = {U.str(15)}", flush=True)
    check(L.rad() < 1e-20 and U.rad() < 1e-20 and L.upper() < U.lower() and L.lower() > 0,
          "ET.3 certified (arb): 0 < L < U, radii < 1e-20")
    Lf, Uf = float(L.upper()), float(U.lower())
    # ET.4(e) attempt to break the imported bounds (BI.4 ceiling/floor, BJ.6 psi bracket) on a grid, n = 1..20 (mpmath, 30 digits)
    gam = mp.log(3) / mp.log(5)
    worst = mp.mpf(0)
    okb = True
    for n in range(1, 21):
        for k in range(-60, 41):
            t = mp.mpf(5) ** (mp.mpf(k) / 4)
            Rn = (mp.mpf(3) / 5) ** n * r_rec(n, t / mp.mpf(5) ** n)
            ceil = min(mp.mpf(1) / 3, mp.mpf(G.numerator) / G.denominator * t ** (gam - 1))
            worst = max(worst, Rn / ceil)
            okb = okb and Rn <= ceil and Rn >= 1 / (3 + t)
            if t <= 200:
                y = t / 5 / mp.mpf(5) ** (n - 1)
                for _ in range(n - 1):
                    y = y * (y + 5)
                okb = okb and t / 5 <= y <= (t / 5) * mp.e ** (t / 50)
    check(okb, f"ET.4(e) grid t = 5^(k/4), k = -60..40, n = 1..20: 1/(3+t) <= R_n <= min(1/3, G t^(gamma-1)) (max ratio "
               f"{mp.nstr(worst, 8)}), t/5 <= psi_(n-1)(t/5) <= (t/5) e^(t/50)")
    vals = []
    for n in range(1, nobs + 1):
        v = scaled_H0(n)
        vals.append(v)
        if n <= len(ES_H0):
            d = abs(v / mp.mpf(3) ** n - mp.mpf(ES_H0[n - 1]))
            check(d < mp.mpf("1e-9"), f"ET.4(b) n = {n}: quadrature H_n(0) = {mp.nstr(v / mp.mpf(3) ** n, 12)} vs ES.4 {ES_H0[n - 1]}")
        check(Lf <= float(v) <= Uf, f"ET.4(c) n = {n}: 3^n H_n(0) = {mp.nstr(v, 12)} in [L, U]")
    print("   OBSERVED (not claimed): 3^n H_n(0) =", ", ".join(mp.nstr(v, 10) for v in vals), flush=True)
    print("   OBSERVED increments:", ", ".join(mp.nstr(vals[i + 1] - vals[i], 3) for i in range(len(vals) - 1)), flush=True)
    print("   OBSERVED ratios H_(n+1)(0)/H_n(0):", ", ".join(mp.nstr(vals[i + 1] / vals[i] / 3, 8) for i in range(len(vals) - 1)),
          flush=True)
    # ET.5 (CONJECTURE, not proved): K = W(0)/2 = ln(5)/2.  Richardson with the observed 3^-n increments.
    K5 = mp.log(5) / 2
    rich = (3 * vals[-1] - vals[-2]) / 2
    print(f"   CONJECTURE ET.5 probe: ln5/2 - 3^n H_n(0) at n = {nobs}: {mp.nstr(K5 - vals[-1], 6)}; times 3^n: "
          f"{mp.nstr((K5 - vals[-1]) * mp.mpf(3) ** nobs, 8)}; Richardson (3 v_n - v_(n-1))/2 - ln5/2 = {mp.nstr(rich - K5, 4)}",
          flush=True)
    check(abs(rich - K5) < mp.mpf("1e-15") and Lf <= float(K5) <= Uf,
          "ET.5 OBSERVED (conjecture not contradicted): Richardson limit equals ln5/2 to 1e-15 and ln5/2 lies in [L, U]")
    print(f"elapsed {time.time() - t_start:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all ET checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else NOBS))
