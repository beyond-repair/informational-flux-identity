#!/usr/bin/env python3
"""Theorem EZ: EY.3's tail inequality, proved for ALL n >= 11 and all real a in [-1, 9/10] (EY had a float grid check for n < 40 only).

Tail (EY.3): T(n,a) = (1/3) 5^((J+1)(1-a)) s0^(1-a-2^J) / (2^J+a-1),  J = floor(n/3),  s0 = (5/2)^n / 5^(J+1).
Claim EZ: for n >= 11 and -1 <= a <= 9/10:  s0 > 1,  T(n,a) <= T(n,-1) <= (2/5)^n;  and n = 10 FAILS (T(10,-1) > (2/5)^10), so 11 is the exact threshold of this tail formula.
  EZ.1 (monotone in a):  d ln T / da = -(J+1) ln5 - ln s0 - 1/(2^J+a-1) < 0 (s0 > 1, 2^J+a-1 >= 2^J-2 > 0 for J >= 3).
  EZ.2 (a = -1):  L(n) = ln T(n,-1) + n ln(5/2) = (c+2) B - (c-1) A - ln3 - ln c,  c = 2^J-2, A = n ln(5/2), B = (J+1) ln5.
        d L/d n = -(c-1) ln(5/2) < 0, so for fixed J the worst n is the smallest allowed: J = 3: n = 11 (computed, arb); J >= 4: n = 3J, where
        L <= c (B - 3J ln2.5) + 2B + 3J ln2.5 (dropping -ln3 - ln c < 0)  with  B - 3J ln2.5 <= 1.6095 (J+1) - 2.7486 J = 1.6095 - 1.1391 J <= -0.7 J (J >= 4),
        so L <= -0.7 J (2^J-2) + 3.219 (J+1) + 2.749 J <= -J (0.7 (2^J-2) - 5.97) + 3.22 < 0  since 0.7 (2^J - 2) >= 9.8 > 9.7 for J >= 4 and J*3.8 > 3.22.
  EZ.3 (arb check):  L(n) < 0 for 11 <= n <= 600 with arb balls, and L(10) > 0.
Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import sys
import time

import math

try:  # python-flint is needed for the arb certificates only; the float check below is pure Python
    from flint import arb, ctx

    ctx.prec = 200
except ImportError:  # pragma: no cover
    arb = None


def L_float(n, a=-1.0):
    """Pure-Python float version of L (margins >= 1.7 at n = 10, 11, so float rounding is irrelevant)."""
    J = n // 3
    A, B = n * math.log(2.5), (J + 1) * math.log(5)
    if A - B <= 0 or 2.0**J + a - 1 <= 0:
        return None
    return (J + 1) * (1 - a) * math.log(5) + (1 - a - 2.0**J) * (A - B) - math.log(2.0**J + a - 1) - math.log(3) + A


def L(n, a=-1):
    """ln( T(n,a) (5/2)^n ) as an arb ball; None if s0 <= 1 or 2^J+a-1 <= 0."""
    J = n // 3
    a = arb(a)
    A = n * arb(2.5).log()
    B = (J + 1) * arb(5).log()
    ln_s0 = A - B
    den = arb(2) ** J + a - 1
    if J < 1 or not (ln_s0 > 0) or not (den > 0):
        return None
    return (J + 1) * (1 - a) * arb(5).log() + (1 - a - arb(2) ** J) * ln_s0 - den.log() - arb(3).log() + A


def check_closed_form(ns=range(11, 80)):
    """EZ.2's formula for L(n) at a = -1 equals the direct definition."""
    for n in ns:
        J = n // 3
        c = arb(2) ** J - 2
        A = n * arb(2.5).log()
        B = (J + 1) * arb(5).log()
        F = (c + 2) * B - (c - 1) * A - arb(3).log() - c.log()
        d = F - L(n)
        if not abs(d) < arb(10) ** -40:
            return False
    return True


def check_analytic_constants():
    ln5, ln25 = arb(5).log(), arb(2.5).log()
    ok = ln5 < arb("1.6095") and 3 * ln25 > arb("2.7486") and ln5 > arb("1.6094")
    # B - 3J ln2.5 <= 1.6095 - 1.1391 J <= -0.7 J for J >= 4
    ok &= (arb("1.6095") - arb("1.1391") * 4) < -arb("0.7") * 4
    # 0.7 (2^J - 2) >= 9.8 at J = 4 (and increasing)
    ok &= arb("0.7") * (arb(2) ** 4 - 2) > arb("9.7")
    ok &= arb("3.8") * 4 > arb("3.22")
    return bool(ok)


def check_all(nmax=600):
    for n in range(11, nmax + 1):
        v = L(n)
        if v is None or not (v < 0):
            return False
    return True


def check_n10_fails():
    v = L(10)
    return v is not None and bool(v > 0)


def check_monotone_a(ns=range(11, 60)):
    for n in ns:
        prev = None
        for i in range(0, 20):
            a = -1 + arb(19) / 10 * i / 19
            v = L(n, a)
            if v is None or (prev is not None and not (v < prev)) or not (v < 0):
                return False
            prev = v
    return True


def main():
    t0 = time.time()
    fails = []
    for name, fn in [("EZ.1 monotone in a (arb grid)", check_monotone_a), ("EZ.2 closed form", check_closed_form),
                     ("EZ.2 analytic constants", check_analytic_constants), ("EZ.3 L(n)<0, 11<=n<=600", check_all),
                     ("EZ.3 n=10 fails", check_n10_fails)]:
        ok = fn()
        print(name, ok)
        if not ok:
            fails.append(name)
    for n in (10, 11, 12, 13, 14, 20):
        print("n", n, "ln(T(n,-1)(5/2)^n) =", L(n).str(12))
    print("FAILS:", fails, "time %.1fs" % (time.time() - t0))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
