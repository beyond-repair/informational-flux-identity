#!/usr/bin/env python3
"""Theorem EU: the value of ET's limit, K = lim 3^n H_n(0) = ln(5)/2 (Conjecture ET.5 proved), with the rate
3^n H_n(0) = ln5/2 - c 3^-n + o(3^-n), c = (1/4) int_0^oo t R(t)^2 / (5 + psi(t/5)) dt.

Setting (as in EQ-ET): p(s) = s(s+5), q_k = p^k + 5, r_n = B_n/A_n = (1/2) v^T (L_n + s)^-1 v (AQ's recursion
r_(k+1)(s) = Phi_s(r_k(p(s))), r_0 = 1/(3+s)), h_n = r_n/(2 q_(n-1)), H_n(0) = int_0^oo h_n (ES.2).  R_n(t) = (3/5)^n r_n(5^-n t),
psi_k(t) = p^k(5^-k t).

  EU.1 (THEOREM, exact, every n >= 0).  With Theta_n(s) = prod_(k<n) (p^k(s)+3)/(p^k(s)+5),
          1/r_n(s) = s/2 + (s+6) Theta_n(s)/2,   equivalently   2 A_n = s B_n + (s+6) prod_(k<n) (p^k+3),   B_n = prod_(k<n)(p^k+5).
        Proof: AQ's step is 1/r_(k+1)(s) = c(s)/r_k(p(s)) - e(s), c = (s+6)/((s+2)(s+5)), e = 2s/(s+2); and c p/2 - e = s/2, so
        u_n = 1/r_n - s/2 satisfies u_n(s) = c(s) u_(n-1)(p(s)), u_0 = (s+6)/2, and c(s)(p(s)+6) = (s+6)(s+3)/(s+5).
  EU.2 (THEOREM, exact, every n >= 1).  r~_n = 2/((s+6) Theta_n), g_n(s) = 3^(n+1) s r~_n(s)/p^n(s) (g_n(0+) = 1),
        m(s) = (s+6)(2s+5)/(3(s+2)(s+5)).  Then (i) 3^n s r~_n/q_(n-1) = g_(n-1) - g_n; (ii) g_(n-1)(p(s)) p'(s)/p(s) = m(s) g_n(s)/s;
        (iii) Frullani along p (p(e)/e -> 5): int_0^oo (g_(n-1) - m g_n) ds/s = ln 5.  Hence, exactly,
          2 * 3^n H_n(0) = ln5 - E_n - D_n,  E_n = int_0^oo (s+4) g_n(s) / (3 (s+2)(s+5)) ds >= 0,
          D_n = (3^n/2) int_0^oo s r_n r~_n / q_(n-1) ds = (3^-n/2) int_0^oo t R_n R~_n / (5 + psi_(n-1)(t/5)) dt >= 0.
  EU.3 (THEOREM).  0 < g_n <= 1, R_n <= R~_n <= 1/3, psi_n >= psi_2 >= t^4/5^8, so E_n <= 5^(5/3-n) (n >= 2) and
        D_n <= (pi/72) 5^(11/2) 3^-n (n >= 3).  Hence 3^n H_n(0) <= ln5/2 for every n >= 1 and 3^n H_n(0) -> ln5/2:
        K = W(0)/2 = ln5/2, W(0) = ln 5 (Conjecture ET.5).  Self-contained: no BG.3, BI.4, BJ.6 needed for EU.1-EU.3.
  EU.4 (THEOREM; uses BD.4 psi_k increasing to psi and BJ.6 psi(v) <= v e^(v/10)).  R_n -> R pointwise with
          1/R(t) = 3 prod_(j>=1) (5/3)(psi(5^-j t)+3)/(psi(5^-j t)+5),  R(t) (3 + psi(t/5)) = (3/5)(5 + psi(t/5)) R(t/5),
        and 3^n (ln5/2 - 3^n H_n(0)) -> c = (1/4) int_0^oo t R^2/(5 + psi(t/5)) dt (dominated convergence; E_n = O(5^-n)).
        c = 2.6718228119692690733... by quadrature (OBSERVED digits, not an arb enclosure).
  EU.5 (checks).  Exact integer polynomials 2A_n = s B_n + (s+6) prod(p^k+3), n = 0..8; EU.2 (i), (ii), 1/r_n - 1/r~_n = s/2 at
        exact rational points, n = 1..6; quadrature of the exact identity n = 1..NOBS; EU.3 bounds; ES.4's certified H_n(0), n <= 8;
        the product formula vs AQ's recursion; the partition of unity sum_j F(5^j t) = 1, F(t) = t R(t)/(5 + psi(t/5)) (DERIVED).

Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_pole_structure as eq  # noqa: E402
import gasket_fractional_corner_h_zero_limit as et  # noqa: E402

FAILS: list[str] = []
NOBS = 24
C_OBS = "2.6718228119692690733"  # EU.4 quadrature value of c (30 digits working precision), OBSERVED


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def poly_identity(n: int) -> bool:
    """EU.1 as integer polynomials: 2 A_n == s B_n + (s+6) prod_(k<n) (p^k + 3)."""
    A, B = eq.eq1_pair(n)
    prod3, f = [1], [0, 1]
    for _ in range(n):
        prod3 = eq.pmul(prod3, eq.padd(f, [3]))
        f = eq.pcomp(f, eq.P_S)
    rhs = eq.padd(eq.pmul([0, 1], B), eq.pmul([6, 1], prod3))
    return [2 * a for a in A] == rhs and B == eq.b_product(n)


def theta(n: int, s):
    out, y = 1, s
    for _ in range(n):
        out = out * (y + 3) / (y + 5)
        y = y * (y + 5)
    return out


def pn(n: int, s):
    y = s
    for _ in range(n):
        y = y * (y + 5)
    return y


def r_closed(n: int, s):
    """EU.1: r_n(s) = 2 / (s + (s+6) Theta_n(s))."""
    return 2 / (s + (s + 6) * theta(n, s))


def rt(n: int, s):
    return 2 / ((s + 6) * theta(n, s))


def g(n: int, s):
    return 3 ** (n + 1) * s * rt(n, s) / pn(n, s)


def m_fun(s):
    return (s + 6) * (2 * s + 5) / (3 * (s + 2) * (s + 5))


def exact_checks(n: int) -> None:
    ok = True
    for s in (Fraction(1, 7), Fraction(2), Fraction(13, 3), Fraction(5, 11)):
        ok = ok and r_closed(n, s) == et.r_rec(n, s)
        ok = ok and 1 / r_closed(n, s) - 1 / rt(n, s) == s / 2
        q = pn(n - 1, s) + 5
        ok = ok and 3 ** n * s * rt(n, s) / q == g(n - 1, s) - g(n, s)                       # EU.2 (i)
        ps = s * (s + 5)
        ok = ok and g(n - 1, ps) * (2 * s + 5) / ps == m_fun(s) * g(n, s) / s                 # EU.2 (ii)
        ok = ok and 0 < g(n, s) <= 1 and rt(n, s) >= r_closed(n, s)
    check(ok, f"EU.1/EU.2 n = {n}: r_n = 2/(s + (s+6)Theta_n) = AQ's recursion; 1/r_n - 1/r~_n = s/2; g_(n-1) - g_n = 3^n s r~_n/q_(n-1); "
              f"g_(n-1)(p) p'/p = m g_n/s; 0 < g_n <= 1 (exact rationals)")


def scaled_parts(n: int):
    """mpmath: (2 * 3^n H_n(0), E_n, 3^n D_n), all as integrals in t = 5^n s."""
    import mpmath as mp
    five_n = mp.mpf(5) ** n
    pts = [0] + [mp.mpf(5) ** k for k in range(-3, 13)] + [mp.inf]

    def parts(t):
        s = t / five_n
        th, y, ypre = mp.mpf(1), s, None
        for _ in range(n):
            th = th * (y + 3) / (y + 5)
            ypre = y
            y = y * (y + 5)
        psi_n, psi_n1 = y, ypre                          # p^n(s), p^(n-1)(s)
        Rt = (mp.mpf(3) / 5) ** n * 2 / ((s + 6) * th)    # R~_n(t)
        R = (mp.mpf(3) / 5) ** n * 2 / (s + (s + 6) * th)  # R_n(t)
        return R, Rt, psi_n, psi_n1, s

    def f_h(t):
        R, Rt, psi_n, psi_n1, s = parts(t)
        return R / (5 + psi_n1)                           # = 2 * 3^n 5^-n h_n(5^-n t)

    def f_e(t):
        R, Rt, psi_n, psi_n1, s = parts(t)
        k = 3 * t * Rt / psi_n if t > 0 else mp.mpf(1)    # g_n(5^-n t)
        return (s + 4) / (3 * (s + 2) * (s + 5)) * k / five_n

    def f_d(t):
        R, Rt, psi_n, psi_n1, s = parts(t)
        return t * R * Rt / (2 * (5 + psi_n1))

    return mp.quad(f_h, pts), mp.quad(f_e, pts), mp.quad(f_d, pts)


def continuum(nn: int = 70):
    """R(t) by EU.4's product formula, psi(t) ~ psi_nn(t); F(t) = t R(t)/(5 + psi(t/5))."""
    import mpmath as mp

    def psi(t):
        y = t / mp.mpf(5) ** nn
        for _ in range(nn):
            y = y * (y + 5)
        return y

    def R(t):
        prod, j = mp.mpf(1), 1
        while True:
            u = t / mp.mpf(5) ** j
            y = psi(u) if u > mp.mpf(10) ** -20 else u * (1 + u / 5)   # psi(u) = u + u^2/5 + O(u^3) below 1e-20 (BD.4)
            prod *= mp.mpf(5) / 3 * (y + 3) / (y + 5)
            if u < mp.mpf(10) ** -40:
                break
            j += 1
        return 1 / (3 * prod)

    return psi, R


def main(nobs: int = NOBS) -> int:
    import mpmath as mp
    t_start = time.time()
    mp.mp.dps = 30
    print("Theorem EU: K = lim 3^n H_n(0) = ln(5)/2 (Conjecture ET.5 proved); rate c 3^-n. Not thrust.", flush=True)
    for n in range(0, 9):
        check(poly_identity(n), f"EU.1 n = {n}: 2 A_n = s B_n + (s+6) prod_(k<n)(p^k+3), B_n = prod_(k<n)(p^k+5) (exact integers)")
    for n in range(1, 7):
        exact_checks(n)
    ln5 = mp.log(5)
    E_bound = lambda n: mp.mpf(5) ** (mp.mpf(5) / 3 - n)  # noqa: E731
    D3_bound = mp.pi / 72 * mp.mpf(5) ** mp.mpf("5.5")
    vals, d3s = [], []
    for n in range(1, nobs + 1):
        twoV, E, D3 = scaled_parts(n)
        D = D3 / mp.mpf(3) ** n
        v = twoV / 2
        vals.append(v)
        d3s.append(D3)
        check(abs(twoV - (ln5 - E - D)) < mp.mpf("1e-20"),
              f"EU.2 n = {n}: 2*3^n H_n(0) = {mp.nstr(twoV, 15)}; ln5 - E_n - D_n = {mp.nstr(ln5 - E - D, 15)} "
              f"(E_n = {mp.nstr(E, 6)}, 3^n D_n = {mp.nstr(D3, 10)})")
        if n <= len(et.ES_H0):
            check(abs(v / mp.mpf(3) ** n - mp.mpf(et.ES_H0[n - 1])) < mp.mpf("1e-9"), f"EU.5 n = {n}: H_n(0) matches ES.4's certified value")
        okb = v <= ln5 / 2 and E >= 0 and D >= 0
        if n >= 2:
            okb = okb and E <= E_bound(n)
        if n >= 3:
            okb = okb and D3 <= D3_bound
        check(okb, f"EU.3 n = {n}: 3^n H_n(0) <= ln5/2; E_n <= 5^(5/3-n) (n>=2); 3^n D_n <= (pi/72) 5^(11/2) = {mp.nstr(D3_bound, 6)} (n>=3)")
    # EU.4: continuum product formula, functional equation, partition of unity, rate constant c
    psi, R = continuum()
    okc = True
    for t in (mp.mpf("0.3"), mp.mpf(7), mp.mpf(40)):
        Rn = (mp.mpf(3) / 5) ** 60 * et.r_rec(60, t / mp.mpf(5) ** 60)
        y = psi(t / 5)
        okc = okc and abs(R(t) - Rn) < mp.mpf("1e-25") and abs(R(t) * (3 + y) - mp.mpf(3) / 5 * (5 + y) * R(t / 5)) < mp.mpf("1e-25")
    check(okc, "EU.4: product formula 1/R = 3 prod (5/3)(psi+3)/(psi+5) = AQ's R_60 (1e-25); R(t)(3+psi(t/5)) = (3/5)(5+psi(t/5))R(t/5)")
    F = lambda t: t * R(t) / (5 + psi(t / 5))  # noqa: E731
    pu = [mp.fsum(F(mp.mpf(5) ** j * t) for j in range(-40, 12)) for t in (mp.mpf(1), mp.mpf("2.9"))]
    check(all(abs(x - 1) < mp.mpf("1e-15") for x in pu), f"EU.4 DERIVED partition of unity sum_j F(5^j t) = 1: {[mp.nstr(x, 20) for x in pu]}")
    pts = [0] + [mp.mpf(5) ** k for k in range(-6, 13)] + [mp.inf]
    c = mp.quad(lambda t: t * R(t) ** 2 / (4 * (5 + psi(t / 5))), pts)
    Kq = mp.quad(lambda t: R(t) / (2 * (5 + psi(t / 5))), pts)
    print(f"   K by continuum quadrature: {mp.nstr(Kq, 25)}; ln5/2 = {mp.nstr(ln5 / 2, 25)}", flush=True)
    check(abs(Kq - ln5 / 2) < mp.mpf("1e-20"), "EU.3 cross-check: (1/2) int R/(5+psi(t/5)) = ln5/2 by quadrature (1e-20)")
    cn = [(ln5 / 2 - vals[i]) * mp.mpf(3) ** (i + 1) for i in range(len(vals))]
    print(f"   c = (1/4) int t R^2/(5+psi(t/5)) dt = {mp.nstr(c, 22)}", flush=True)
    print("   3^n (ln5/2 - 3^n H_n(0)), n = 1..:", ", ".join(mp.nstr(x, 12) for x in cn), flush=True)
    print("   3^n D_n / 2:", ", ".join(mp.nstr(x / 2, 12) for x in d3s[-4:]), flush=True)
    gaps = [abs(x - c) for x in cn]
    check(abs(c - mp.mpf(C_OBS)) < mp.mpf("1e-18") and gaps[-1] < mp.mpf("1e-4") and gaps[-1] < gaps[-5],
          f"EU.4 rate: 3^n(ln5/2 - 3^n H_n(0)) -> c = {mp.nstr(c, 20)}; gap at n = {nobs}: {mp.nstr(gaps[-1], 4)} < 1e-4 and < gap at n - 4")
    print(f"elapsed {time.time() - t_start:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EU checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else NOBS))
