#!/usr/bin/env python3
"""Theorem EW: W(alpha) is analytic on Re alpha < 1 and vanishes at alpha = 2 pi i k / ln5, k != 0 (W(0) = ln5, EU).

Setting (as in EU/EV): p(s) = s(s+5); psi = lim p^k(5^-k t), psi(5v) = p(psi(v)) (BD.4);
R = lim R_n (BG.3), 1/R(t) = 3 prod_(j>=1) (5/3)(psi(5^-j t)+3)/(psi(5^-j t)+5) (EU.4);
F(t) = t R(t)/(5 + psi(t/5)) >= 0, sum_(j in Z) F(5^j t) = 1 (EU.4 partition of unity), W(alpha) = int_0^oo t^-alpha R/(5+psi(t/5)) dt
= int_0^oo F(t) t^(-alpha-1) dt.

  EW.1 (THEOREM).  F(t) <= t/15 and F(t) <= G t^gamma / psi(t/5) (BI.4), psi(t) >= exp(b t^(log_5 2)) with b = ln(psi(1))/2 > 0
        (EV.1: psi(5v) >= psi(v)^2, psi increasing, psi(1) > 1), so F decays faster than any power and W(alpha) = int F t^(-alpha-1) is
        absolutely convergent and holomorphic for Re alpha < 1.
  EW.2 (THEOREM).  With omega_k = 2 pi k / ln5 and k a nonzero integer:  W(i omega_k) = 0.   (W(0) = ln5, EU.)
        Proof: Tonelli (F >= 0, int F dt/t = ln5) and the partition of unity on [1,5]:
        W(i omega_k) = sum_j int_(5^j)^(5^(j+1)) F t^(-i omega_k - 1) dt = sum_j int_1^5 F(5^j t) t^(-i omega_k - 1) dt (5^(-i j omega_k) = 1)
        = int_1^5 t^(-i omega_k - 1) dt = (1 - 5^(-i omega_k)) / (i omega_k) = 0.
  EW.3 (THEOREM).  3^n H_n(i omega_k) -> 0 (n -> oo) for every k != 0, where H_n(alpha) = int_0^oo s^-alpha h_n ds converges absolutely
        (h_n rational, h_n(0) finite, O(s^-2)) and, because 5^(i n omega_k) = 1, ET.1 holds verbatim:
        3^n H_n(i omega_k) = int_0^oo t^(-i omega_k) R_n(t) / (2 (5 + psi_(n-1)(t/5))) dt.  |integrand| = the alpha = 0 integrand, so ET.2's
        dominator (BI.4's ceiling) and BG.3's pointwise limits give the limit (1/2) W(i omega_k) = 0.   (no rate.)
  EW.4 (THEOREM, certified arb).  |W(i omega_k)| < 1e-27 for k = 1, 2, and W(0) contains ln5 (cross-check), by acb_calc_integrate with EV's
        evaluator on [5^-M, 5^8] (M = LOW) plus tails: int_0^(5^-M) <= 5^-M/15, int_(5^8)^oo as in EV.4.
  EW.5 (OBSERVED, not proved).  3^n H_n(i omega_1) for n = 1..6 (float log-grid quadrature of ET.1's integral), and their size.

Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import math
import sys
import time

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_h_zero_rate_constant as ev  # noqa: E402

FAILS = ev.FAILS
check = ev.check
LOW = 40


def block_edges(arb):
    return [arb(5) ** m for m in range(-LOW, 9)]


def F_float(t: float) -> float:
    return t * ev.R_float(t) / (5.0 + ev.psi_float(t / 5.0))


def partition_defect(t: float, J: int = 25) -> float:
    """|sum_j F(5^j t) - 1| in floats (EU.4's partition of unity), j = -J..8."""
    return abs(sum(F_float(t * 5.0 ** j) for j in range(-J, 9)) - 1.0)


def mellin_float(k: int, lo: float = 1e-14, hi_exp: int = 8, n: int = 6001) -> complex:
    """trapezoid in u = ln t of int F(e^u) e^(-i omega_k u) du (float; CI)."""
    om = 2 * math.pi * k / math.log(5)
    u0, u1 = math.log(lo), hi_exp * math.log(5)
    h = (u1 - u0) / (n - 1)
    s = 0j
    for i in range(n):
        u = u0 + i * h
        w = 0.5 if i in (0, n - 1) else 1.0
        s += w * F_float(math.exp(u)) * complex(math.cos(om * u), -math.sin(om * u))
    return s * h


def main() -> int:
    t0 = time.time()
    print("Theorem EW: W(alpha) vanishes at alpha = 2 pi i k / ln5, k != 0. Not thrust.", flush=True)
    acb, arb, RY, psi, f_c, f_k = ev.make_arb()
    # EW.1: psi(1) > 1 and the superpolynomial decay constant
    p1 = psi(1).real
    check(bool(p1 > 1), f"EW.1: psi(1) = {p1.str(10)} > 1, so psi(5^m) >= psi(1)^(2^m) and F decays like exp(-b t^0.43)")
    tb, ratio_ok, R8, p7 = ev.tail_bound(arb, RY, psi)
    k_tail = (2 * 4 * arb(5) ** 8 * R8 / (2 * (5 + p7))).upper()  # EV's tail for int R/(2(5+psi)); W has twice
    w_tail = 2 * k_tail
    low_tail = arb(5) ** (-LOW) / 15
    edges = block_edges(arb)
    out = {}
    for k in (0, 1, 2):
        om = arb(2) * arb.pi() * k / arb(5).log()
        def f(z, analytic, om=om):
            R, y1 = RY(z)
            return R / (5 + y1) * (-(acb(0, 1) * om) * z.log()).exp()
        tot = acb(0)
        for a, b in zip(edges, edges[1:]):
            tot += acb.integral(f, a, b, rel_tol=arb(2) ** -100, abs_tol=arb(2) ** -110)
        slack = arb(0, (w_tail + low_tail).upper())
        tot = tot + acb(slack, slack)
        out[k] = tot
        print(f"   k={k}: W(i omega_k) ball = {tot.real.str(12)} + i {tot.imag.str(12)}", flush=True)
    check(out[0].real.contains(arb(5).log()) and out[0].imag.contains(0), "EW.4 cross-check: W(0) ball contains ln5 (EU)")
    for k in (1, 2):
        a = abs(out[k])
        check(bool(a.upper() < arb("1e-27")), f"EW.4 THEOREM (certified arb): |W(i omega_{k})| <= {a.upper().str(5)} < 1e-27 (k = {k}); exact value 0 by EW.2")
    # EW.5: observed finite-n values via EU's exact closed form
    import mpmath as mp
    import gasket_fractional_corner_h_zero_value as eu
    mp.mp.dps = 20
    om = 2 * math.pi / math.log(5)
    vals = []
    for n in range(1, 7):
        s = 0j
        u0, u1, N = math.log(1e-12), math.log(1e12 if n <= 3 else 1e6), 4001
        h = (u1 - u0) / (N - 1)
        for i in range(N):
            u = u0 + i * h
            x = math.exp(u)
            sx = x / 5.0 ** n
            Rn = (3.0 / 5.0) ** n * float(eu.r_closed(n, sx))
            q = float(eu.pn(n - 1, sx)) + 5.0
            w = 0.5 if i in (0, N - 1) else 1.0
            s += w * Rn / (2 * q) * x * complex(math.cos(om * u), -math.sin(om * u))
        vals.append(s * h)
        print(f"   EW.5 (OBSERVED): n={n}: 3^n H_n(i omega_1) ~ {s*h:.3e}   |.| = {abs(s*h):.3e}", flush=True)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EW checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
