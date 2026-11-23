#!/usr/bin/env python3
"""Theorem EY: an explicit geometric rate (2/5)^n for EX.1's limit rho^-n H_n(alpha) -> W(alpha)/2 on the strip -1 <= Re alpha <= 9/10.

Notation (EX.1): f_n(t) = R_n(t) / (2 (5 + psi_(n-1)(t/5))), f(t) = R(t) / (2 (5 + psi(t/5))), G_n(alpha) = int t^-alpha f_n = rho^-n H_n(alpha),
W(alpha) = 2 int t^-alpha f.  Closed form (EU.1 with s = 5^-n t, p^k(5^-n t) = psi_k(5^(k-n) t)):
    1/R_n(t) = 3 P_n(t) + (t/2) 3^-n (1 + (3/5)^n P_n(t)),   P_n(t) = prod_(j=1..n) g(psi_(n-j)(5^-j t)),   g(y) = (5/3)(y+3)/(y+5),
    1/R(t)   = 3 P(t),                                         P(t)   = prod_(j>=1)   g(psi(5^-j t)).
  EY.1 (THEOREM; Koenigs-type psi_m vs psi).  0 <= log psi(t) - log psi_m(t) <= (t/10)(2/5)^m.
        psi(t) = p^m(psi(w)), psi_m(t) = p^m(w), w = 5^-m t; BJ.6: log(psi(w)/w) <= w/10; d log p / d log y = 1 + y/(y+5) <= 2.
  EY.2 (THEOREM; pointwise).  For 0 < t <= (5/2)^n:  |f_n(t)/f(t) - 1| <= (t/3)(2/5)^n.
        Upper: R_n/R <= P/P_n <= exp(0.0127 x + 0.037 t 5^-n) <= exp(0.05 x), x = t (2/5)^n <= 1 (d log g/d log y <= 0.127; sum of 2^-j; tail
        of P by log g <= 2y/15, psi(u) <= 1.106 u for u <= 1); (5+psi)/(5+psi_(n-1)) <= psi/psi_(n-1) <= exp(0.05 x); e^(0.1 x) - 1 <= 0.106 x.
        Lower: R_n/R >= 1/(1 + extra/(3 P_n)) >= 1 - (t/3) 3^-n, extra/(3P_n) <= (t/6) 3^-n (1 + 3/5); (5+psi)/(5+psi_(n-1)) >= 1.
  EY.3 (THEOREM; rate).  For every complex alpha with -1 <= Re alpha <= 9/10 and every n >= 11:
        |rho^-n H_n(alpha) - W(alpha)/2| <= (W(alpha-1)/6 + 1) (2/5)^n.
        Proof: on t <= T_n = (5/2)^n EY.2 gives (2/5)^n/3 * int t^(1-a) f = (2/5)^n W(alpha-1)/6 (finite: Re(alpha-1) < 1).  Beyond T_n,
        f_n + f <= 1/(3 (5 + psi_(n-1)(t/5))) (R_n, R <= 1/3, psi >= psi_(n-1)), psi_(n-1)(t/5) >= (t/5^(J+1))^(2^J), J = floor(n/3) <= n-1, so
        the tail <= (1/3) 5^((J+1)(1-a)) s0^(1-a-2^J)/(2^J+a-1), s0 = T_n/5^(J+1) > 1; decreasing in a, so a = -1 is the worst case; it is
        <= (2/5)^n for all n >= 11 (ratio 0.168 at n = 11, 1.6e-8 at n = 12, then super-geometric).
  EY.4 (OBSERVED, not proved).  The true rate is 3^-n: 3^n (W/2 - G_n) -> 304.7 (alpha=-1), 24.89 (-1/2), 0.3976 (1/2), 0.1152 (9/10).
Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import math
import sys
import time


def p(y):
    return y * (y + 5.0)


def psi_m(m, t):
    y = t / 5.0**m
    for _ in range(m):
        y = p(y)
    return y


def psi(t):
    m = max(0, int(math.log(t + 1) / math.log(5)) + 6)
    y = psi_m(30, t / 5.0**m)
    for _ in range(m):
        y = p(y)
        if y > 1e300:
            return math.inf
    return y


def g(y):
    return (5.0 / 3.0) * (y + 3.0) / (y + 5.0)


def inv_R_n(n, t):
    s = t / 5.0**n
    prod, y = 1.0, s
    for _ in range(n):
        prod *= (y + 3.0) / (y + 5.0)
        y = p(y)
        if y > 1e150:
            return math.inf
    return (5.0 / 3.0) ** n * (s / 2 + (s + 6) / 2 * prod)


def inv_R(t):
    P = 1.0
    for j in range(1, 45):
        P *= g(psi(t / 5.0**j))
    return 3.0 * P


def f_n(n, t):
    v = inv_R_n(n, t) * 2 * (5 + psi_m(n - 1, t / 5))
    return 0.0 if not math.isfinite(v) else 1.0 / v


def f_inf(t):
    v = inv_R(t) * 2 * (5 + psi(t / 5))
    return 0.0 if not math.isfinite(v) else 1.0 / v


def check_eps_lemma(ms=range(0, 9), ts=(0.01, 0.1, 1, 3, 10, 30, 100, 300, 1000)):
    for m in ms:
        for t in ts:
            lhs = math.log(psi(t) / psi_m(m, t))
            if not (-1e-12 <= lhs <= t / 10 * 0.4**m * (1 + 1e-9)):
                return False
    return True


def check_pointwise(ns=range(1, 9), pts=21):
    worst = 0.0
    for n in ns:
        T = 2.5**n
        for k in range(pts):
            t = 1e-3 * (T * 1e3) ** (k / (pts - 1))
            r = abs(f_n(n, t) / f_inf(t) - 1) / (t * 0.4**n / 3)
            worst = max(worst, r)
    return worst


def tail_bound(n, a):
    J = n // 3
    s0 = 2.5**n / 5.0 ** (J + 1)
    if J < 1 or s0 <= 1 or 2**J + a - 1 <= 0:
        return None
    return 5.0 ** ((J + 1) * (1 - a)) * s0 ** (1 - a - 2**J) / (2**J + a - 1) / 3


def tail_ok(nmax=40):
    for n in range(11, nmax):
        for a in [-1 + 1.9 * i / 19 for i in range(20)]:
            tb = tail_bound(n, a)
            if tb is None or tb > 0.4**n:
                return False
    return True


def grid(h=0.01):
    lo, hi = math.log(1e-12), math.log(2e8)
    N = int((hi - lo) / h) // 2 * 2
    h = (hi - lo) / N
    X = [lo + i * h for i in range(N + 1)]
    return h, X


def make_integrals(nmax=12):
    h, X = grid()
    T = [math.exp(x) for x in X]
    finf = [f_inf(t) for t in T]
    fn = {n: [f_n(n, t) for t in T] for n in range(1, nmax + 1)}

    def I(F, a):
        s = 0.0
        for i, (t, v) in enumerate(zip(T, F)):
            w = 1 if i in (0, len(T) - 1) else (4 if i % 2 else 2)
            s += w * t ** (1 - a) * v
        return 1e-12 ** (1 - a) / (30 * (1 - a)) + s * h / 3

    return I, finf, fn


def main():
    t0 = time.time()
    fails = []
    if not check_eps_lemma():
        fails.append("EY.1")
    w = check_pointwise()
    print("EY.2 max |f_n/f-1| / bound =", w)
    if w > 1:
        fails.append("EY.2")
    if not tail_ok():
        fails.append("EY.3 tail")
    I, finf, fn = make_integrals()
    for a in (-1.0, -0.5, 0.5, 0.9):
        W = 2 * I(finf, a)
        C = 2 * I(finf, a - 1) / 6
        print("alpha", a, "W", W, "C=W(a-1)/6", C)
        for n in range(1, 13):
            gap = abs(I(fn[n], a) - W / 2)
            q = 0.4**n
            tb = tail_bound(n, a)
            flag = ""
            if tb is not None and gap > C * q + tb:
                flag = "VIOL-explicit"
            if n >= 11 and gap > (C + 1) * q:
                flag = "VIOL-theorem"
            if flag:
                fails.append((a, n, flag))
            print(n, "gap %.4e  gap/q %.4f  (C+1)q %.4e  3^n gap %.4f" % (gap, gap / q, (C + 1) * q, gap * 3**n), flag)
    print("FAILS:", fails, "time %.1fs" % (time.time() - t0))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
