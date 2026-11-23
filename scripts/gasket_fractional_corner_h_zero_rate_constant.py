#!/usr/bin/env python3
"""Theorem EV: a certified enclosure of EU's rate constant
    c = lim 3^n (ln5/2 - 3^n H_n(0)) = (1/4) int_0^oo t R(t)^2 / (5 + psi(t/5)) dt     (EU.4),
which EU recorded only as quadrature digits (OBSERVED, "not an arb enclosure").

Setting (as in EU): p(s) = s(s+5); psi = lim p^k(5^-k t) (BD.4: psi(5v) = p(psi(v)), psi_k increasing to psi);
1/R(t) = 3 prod_(j>=1) (5/3)(psi(5^-j t)+3)/(psi(5^-j t)+5) (EU.4).

  EV.1 (THEOREM, exact; Taylor coefficients of psi).  psi(w) = sum_(m>=1) a_m w^m with a_1 = 1 and
          a_m (5^m - 5) = sum_(i+j=m, i,j>=1) a_i a_j   (m >= 2),   so a_2 = 1/20, and  0 <= a_m <= 20^(1-m) for every m >= 2
        (induction: sum_(i+j=m) a_i a_j <= (m-1) 20^(2-m) and 20 (m-1) <= 5^m - 5 for m >= 2).  Hence, for |w| < 20,
          |psi(w) - w - w^2/20| <= |w|^3 / (400 (1 - |w|/20)),
        psi >= 0 and increasing on [0, oo), and psi(5v) = psi(v)(psi(v)+5) >= psi(v)^2 there.
        (ASSUMPTION, BD.4 + Koenigs' linearization at the repelling fixed point 0 of p, p'(0) = 5: this series, radius >= 20,
        is BD.4's psi.  Checked: the arb psi agrees with EU's psi_70 iteration, EV.5.)
  EV.2 (THEOREM; truncation of EU.4's product).  Each factor is f_j = 1 + 2 y_j/(3(y_j+5)), y_j = psi(5^-j z).  If |5^-J z| = w <= 2^-130,
        then |y_j| <= 1.01 |5^-j z| for j >= J, sum_(j>J) |5^-j z| <= w/4, so |prod_(j>J) f_j - 1| <= exp(delta) - 1,
        delta = 2 (1.01) (w/4) / (3 (5 - 1.01 w)).  y_J is enclosed by EV.1 and y_(j-1) = p(y_j) (BD.4); y_1 = psi(z/5).
        So z |-> z R(z)^2 / (4(5 + psi(z/5))) has a rigorous complex-ball evaluator (meromorphic: any ball meeting a pole or a zero
        of a denominator returns non-finite, which arb's integrator treats as "not holomorphic here" and subdivides).
  EV.3 (THEOREM; tail).  R is decreasing on [0, oo) (each factor increases with psi >= 0) and R <= 1/3.  For m >= 8 the integrand on
        [5^m, 5^(m+1)] is <= 5^(m+1) R(5^8)^2 / (4 (5 + psi(5^(m-1)))), and psi(5^(m-1)) >= psi(5^7)^(2^(m-8)), so
          int_(5^8)^oo <= 2 * 5^17 R(5^8)^2 / (5 + psi(5^7))   (ratio of consecutive blocks <= 25/psi(5^7) <= 1/2).
  EV.4 (THEOREM, certified arb).  arb's rigorous integrator (acb_calc_integrate via python-flint) on [0,1], [1,5], ..., [5^7, 5^8]
        plus EV.3's tail gives
          c in [C_LO, C_HI]  (width < 1e-30),   c = 2.67182281196926907329576690025...
        Cross-check of the evaluator: the same integrator gives K = (1/2) int R/(5 + psi(t/5)) dt in a ball that contains ln5/2
        (EU.3's theorem, independent of this evaluator).
  EV.5 (checks).  EV.1's bound exactly for m <= 60; the ball R(t) contains EU's mpmath product value and AQ's R_60 at t = 0.3, 7, 40, 300;
        the enclosure contains EU's OBSERVED C_OBS = 2.6718228119692690733.

Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0])

FAILS: list[str] = []
PREC = 224
TINY_EXP = 130
EDGES = [0] + [5 ** k for k in range(0, 9)]
C_LO = "2.6718228119692690732957669002"   # certified (EV.4), rounded outward from the arb ball
C_HI = "2.6718228119692690732957669003"
C_OBS = "2.6718228119692690733"           # EU.4's quadrature value (OBSERVED)


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def psi_coeffs(M: int) -> list[Fraction]:
    """EV.1: a_1..a_M exactly (index 0 unused)."""
    a = [Fraction(0), Fraction(1)]
    for m in range(2, M + 1):
        conv = sum(a[i] * a[m - i] for i in range(1, m))
        a.append(conv / (5 ** m - 5))
    return a


def coeff_bound_ok(M: int = 60) -> bool:
    a = psi_coeffs(M)
    ok = a[2] == Fraction(1, 20)
    for m in range(2, M + 1):
        ok = ok and 0 <= a[m] <= Fraction(1, 20 ** (m - 1)) and 20 * (m - 1) <= 5 ** m - 5
    return ok


def psi_float(t: float) -> float:
    """float psi via EV.1's series at 5^-k t and BD.4's p-iteration (CI helper)."""
    k = 0
    w = t
    while abs(w) > 1e-8:
        w /= 5.0
        k += 1
    y = w + w * w / 20.0
    for _ in range(k):
        y = y * (y + 5.0)
    return y


def R_float(t: float) -> float:
    prod, j = 1.0, 1
    while True:
        u = t / 5.0 ** j
        y = psi_float(u)
        prod *= 5.0 / 3.0 * (y + 3.0) / (y + 5.0)
        if u < 1e-17:
            break
        j += 1
    return 1.0 / (3.0 * prod)


def make_arb():
    from flint import acb, arb, ctx
    ctx.prec = PREC
    tiny = arb(2) ** (-TINY_EXP)

    def psi_small(w):
        a = abs(w)
        assert a.upper() < 1
        err = (a ** 3 / 400 / (1 - a / 20)).upper()
        return w + w * w / 20 + acb(arb(0, err), arb(0, err))

    def RY(z):
        """(R(z), psi(z/5)) as complex balls (EV.1, EV.2)."""
        z = acb(z)
        a = abs(z).upper()
        J = 1
        while a / arb(5) ** J > tiny:
            J += 1
        wJ = z / arb(5) ** J
        y = psi_small(wJ)
        aw = abs(wJ).upper()
        delta = arb(2) * arb(101) / 100 * aw / 4 / (3 * (5 - arb(101) / 100 * aw))
        e = (delta.exp() - 1).upper()
        prod = acb(1) + acb(arb(0, e), arb(0, e))
        for j in range(J, 0, -1):
            prod = prod * (arb(5) / 3) * (y + 3) / (y + 5)
            if j > 1:
                y = y * (y + 5)
        return 1 / (3 * prod), y

    def psi(z):
        return RY(acb(z) * 5)[1]

    def f_c(z, analytic):
        R, y1 = RY(z)
        return z * R * R / (4 * (5 + y1))

    def f_k(z, analytic):
        R, y1 = RY(z)
        return R / (2 * (5 + y1))

    return acb, arb, RY, psi, f_c, f_k


def tail_bound(arb, RY, psi):
    """EV.3: int_(5^8)^oo of the c-integrand, upper bound (and the ratio condition)."""
    R8 = RY(5 ** 8)[0].real
    p7 = psi(5 ** 7).real
    ratio_ok = bool(arb(25) / p7 <= arb(1) / 2)
    return (2 * arb(5) ** 17 * R8 * R8 / (5 + p7)).upper(), ratio_ok, R8, p7


def integrate(acb, arb, f, verbose=True):
    tot = acb(0)
    for a, b in zip(EDGES, EDGES[1:]):
        v = acb.integral(f, a, b, rel_tol=arb(2) ** -110, abs_tol=arb(2) ** -120)
        if verbose:
            print(f"     [{a}, {b}]: {v.real}", flush=True)
        tot += v
    return tot


def main() -> int:
    t0 = time.time()
    print("Theorem EV: certified enclosure of EU's rate constant c = (1/4) int t R^2/(5+psi(t/5)). Not thrust.", flush=True)
    check(coeff_bound_ok(60), "EV.1: a_2 = 1/20 and 0 <= a_m <= 20^(1-m), 20(m-1) <= 5^m - 5, m = 2..60 (exact rationals)")
    acb, arb, RY, psi, f_c, f_k = make_arb()
    # EV.5: evaluator vs EU's mpmath product and AQ's R_60
    import mpmath as mp
    mp.mp.dps = 40
    import gasket_fractional_corner_h_zero_value as eu
    import gasket_fractional_corner_h_zero_limit as et
    mpsi, mR = eu.continuum()
    ok = True
    for t in ("0.3", "7", "40", "300"):
        Rb, yb = RY(acb(t))
        Rm = mR(mp.mpf(t))
        R60 = (mp.mpf(3) / 5) ** 60 * et.r_rec(60, mp.mpf(t) / mp.mpf(5) ** 60)
        ok = ok and abs(mp.mpf(Rb.real.mid().str(40, radius=False)) - Rm) < mp.mpf("1e-25") and abs(Rm - R60) < mp.mpf("1e-25")
        ok = ok and float(Rb.real.rad()) < 1e-40 and Rb.imag.contains(0)
        ok = ok and abs(mp.mpf(yb.real.mid().str(40, radius=False)) - mpsi(mp.mpf(t) / 5)) < mp.mpf("1e-25")
    check(ok, "EV.5: arb R(t), psi(t/5) (radius < 1e-40) agree with EU's mpmath product and AQ's R_60 at t = 0.3, 7, 40, 300 (1e-25)")
    tb, ratio_ok, R8, p7 = tail_bound(arb, RY, psi)
    check(ratio_ok and tb < arb("1e-50"), f"EV.3: R(5^8) = {R8.str(8)}, psi(5^7) = {p7.str(8)}, 25/psi(5^7) <= 1/2; tail beyond 5^8 <= {tb.str(5)}")
    print("   EV.4: integrating c (arb, rigorous):", flush=True)
    cI = integrate(acb, arb, f_c)
    c_ball = cI.real + arb(0, tb) + arb(tb) / 2   # tail in [0, tb]
    check(cI.imag.contains(0), "EV.4: imaginary part of the c-integral contains 0")
    c_lo, c_hi = c_ball.lower(), c_ball.upper()
    print(f"   c in [{c_lo.str(35, radius=False)}, {c_hi.str(35, radius=False)}], ball {c_ball.str(35)}", flush=True)
    check(bool(arb(C_LO) <= c_lo) and bool(c_hi <= arb(C_HI)) and bool((c_hi - c_lo) < arb("1e-30")),
          f"EV.4 THEOREM: c in [{C_LO}, {C_HI}] (certified, width < 1e-30)")
    check(abs(float(arb(C_OBS) - c_ball.mid())) < 1e-19,
          f"EV.5: EU's OBSERVED C_OBS = {C_OBS} agrees with the certified enclosure to its 20 digits")
    print("   EV.4 cross-check: K = (1/2) int R/(5+psi(t/5)) (arb):", flush=True)
    kI = integrate(acb, arb, f_k, verbose=False)
    k_tail = (2 * 4 * arb(5) ** 8 * R8 / (2 * (5 + p7))).upper()   # EV.3's block argument without the t factor
    k_ball = kI.real + arb(0, k_tail) + arb(k_tail) / 2
    check(k_ball.contains(arb(5).log() / 2) and bool(k_ball.rad() < arb("1e-30")),
          f"EV.4 cross-check: K ball {k_ball.str(30)} contains ln5/2 = {(arb(5).log() / 2).str(30)} (EU.3)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EV checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
