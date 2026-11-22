#!/usr/bin/env python3
"""Theorem EO: level 1 in closed form, and the full sign pattern at alpha = 0.88 and 0.89 (N = 1).  EN kept "level 1 is not decided
(EJ.2 (b): +5.057e-3; EN.1 at n = 1 needs U(1), not computed)" and left open "level 1 at 0.88 (and 0.89)".  Level 1 is a finite object:
the tail route (EN.1 / EJ.2 (b)) is not needed there.

  EO.1 (THEOREM, exact; every alpha in (0, 1) for the integrals, every alpha in (gamma, 1) for Y; AQ (1), DZ.1, BL, EB as stated).
       AQ (1) at k = 0 (u_0 = 3 + s, u_1(s) = beta(s) u_0(p(s)) - e(s)) gives (sympy) u_1 = (s^2 + 7s + 9)/(s + 5), r_1 = (s + 5)/(s^2 + 7s + 9),
       and BL's recursion with eta_0 = -1/2 (delta_0 = 3/2) gives eta_1(z) = 1 - z/(z+2) - (3/2) beta(z) = 1/(2(z + 5)).  With
       a_(1,2) = (7 -+ sqrt 13)/2 > 0 (s^2 + 7s + 9 = (s + a_1)(s + a_2)) and int_0^oo s^-a/(s + a) ds = pi a^-a / sin(pi a):
         J_1 = pi/(sqrt13 sin(pi a)) [(5 - a_1) a_1^-a + (a_2 - 5) a_2^-a],   H_1 = pi/(sqrt13 sin(pi a)) (a_1^-a - a_2^-a)/2,
         1 - delta_1^(1) = H_1/J_1 = (a_1^-a - a_2^-a) / (2 [(5 - a_1) a_1^-a + (a_2 - 5) a_2^-a]),   D'_1 = J_1/rho.
       With om = 1 - delta_1, om_1 = 1 - delta_1^(1), nu < 0 (DZ.3): X_1 = rho (om_1 - om), C_delta = -nu om, Y_1 = rho (X_1 - C_delta),
       Y_oo = nu^2 om (EB), so (sympy)
         Y_1 - Y_oo = rho^2 (om_1 - om) - om |nu| (rho + |nu|),
       decreasing in om and in |nu| (om, |nu| > 0).  Only EF.1's om and |nu| enter; no eps_n, eta_n, L, U.
  EO.2 (THEOREM, certified; arb, 128 bits; EF.1 point cells 4000/32000 per unit, EN's certify_point_en): alpha = 0.88, N = 1:
       Y_1 - Y_oo < 0.  With EJ.4 (n = 2..8), EN.2 (n = 9) and EJ.4 / EN.2 (every n >= 10): Y_n < Y_oo for n = 1..9 and Y_n > Y_oo
       for every n >= 10, i.e. the full sign pattern for all n >= 1, exactly one sign change, between 9 and 10.
       alpha = 0.89, N = 1: Y_1 - Y_oo < 0; with EJ.4 / EM.3 (n = 2..14, every n >= 15): exactly one sign change for all n >= 1, 14 -> 15.
  EO.3 (THEOREM, certified; --full; point cells 1000/8000 per unit): the sign of Y_1 - Y_oo (N = 1) at EJ.4's other exponents:
       > 0 at 0.7 and 0.8, < 0 at 0.85 (with EJ.4's onset 6: at 0.85 Y_n - Y_oo changes sign at least once, between 1 and 6), < 0 at
       0.9, 0.905, 0.91 (reproduces EJ.3 / EJ.4's level 1 by a route independent of EJ.1 / EG.1).
  EO.4 (OBSERVED, floats): om_1 by quadrature of EJ's r_1 recursion and by mpmath agrees with EO.1's closed form; direct level-n floats
       of Y_n - Y_oo (trapezoid of x^-a R_n and x^-a R_n eta_1, EN's bounded recursion) for n = 1..12 at all eight exponents: level 1
       equals the closed form, its sign matches EO.2 / EO.3, and at 0.88 / 0.89 the float signs match the certified patterns.

Attempts to break (0.88): EN.1 at n = 1 (U(1) by EL.1, L(1) by EM.1, 1600 cells/unit) must contain the exact D'_1 (it does:
[3.573152, 6.742495] against 3.69222949).  EJ.2 (b) at n = 1 with the exact D'_1 in place of D'_oo was tried as a second, tail-based
route: it lowers the bound from +0.005057 to +8.635e-05 but stays > 0 (undecided; L(1), kappa_lo(1, 1) are loose at k0 = 1), so the
tail route does not reach level 1 and only the closed form decides it (OBSERVED, not a failure: (b) is a one-sided criterion).

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_two_sided_rate as ej  # noqa: E402
import gasket_fractional_corner_exact_0905_n3 as el  # noqa: E402
import gasket_fractional_corner_exact_091 as em  # noqa: E402
import gasket_fractional_corner_exact_088 as en  # noqa: E402

FAILS: list[str] = []
MAIN = {"0.88": (4000, 32000), "0.89": (4000, 32000)}                 # EO.2 (EN's / EM.3's point grids)
SIDE = ("0.7", "0.8", "0.85", "0.9", "0.905", "0.91")                 # EO.3 (--full)
SIDE_GRID = (1000, 8000)
SIGN1 = {"0.7": 1, "0.8": 1, "0.85": -1, "0.88": -1, "0.89": -1, "0.9": -1, "0.905": -1, "0.91": -1}   # EO.2 / EO.3: sign(Y_1 - Y_oo)
EN2 = {"Dp": (5.4676986, 5.4679626), "nu": (0.49916739, 0.4993144), "om": (0.068040179, 0.06804675)}   # EN.2's 0.88 brackets
EM3_NU = (0.4663815, 0.4665156)    # EM.3's |nu| at 0.89
EN_LEVEL1 = 0.005057               # EJ.2 (b) at n = 1, 0.88, with D'_n <= D'_oo (EN.2, undecided)
PATTERN = {"0.88": (9, 10), "0.89": (14, 15)}                         # last negative level, onset (EN.2, EM.3)
EJ4_ONSET85 = 6                    # EJ.4: Y_n > Y_oo for every n >= 6 at 0.85


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
A1, A2 = (7 - np.sqrt(13)) / 2, (7 + np.sqrt(13)) / 2


def om1_float(al: float) -> float:
    """EO.1 in floats: 1 - delta_1^(1) = (a_1^-a - a_2^-a) / (2 [(5 - a_1) a_1^-a + (a_2 - 5) a_2^-a])."""
    p1, p2 = A1 ** (-al), A2 ** (-al)
    return (p1 - p2) / (2 * ((5 - A1) * p1 + (A2 - 5) * p2))


def j1_float(al: float) -> float:
    """EO.1 in floats: J_1 = pi/(sqrt13 sin(pi a)) [(5 - a_1) a_1^-a + (a_2 - 5) a_2^-a]."""
    return np.pi / (np.sqrt(13) * np.sin(np.pi * al)) * ((5 - A1) * A1 ** (-al) + (A2 - 5) * A2 ** (-al))


def level1_gap_float(al: float, nu_abs: float, om: float) -> float:
    """EO.1: Y_1 - Y_oo = rho^2 (om_1 - om) - om |nu| (rho + |nu|)."""
    rho = 5**al / 3
    return rho**2 * (om1_float(al) - om) - om * nu_abs * (rho + nu_abs)


def level_quad_float(al: float, n: int, lo: int = -200, hi: int = 200, m: int = 100):
    """Floats (independent route): D'_n = int x^-a R_n and 1 - delta_1^(n) = int x^-a R_n eta_1(p^(n-1)(5^-n x)) / D'_n, trapezoid in
    log x (EN's bounded recursion for R_n; eta_1(z) = 1/(2(z + 5)))."""
    u = np.linspace(lo, hi, (hi - lo) * m + 1) * np.log(5.0)
    x = np.exp(u)
    R = en.rk_scaled_float(x, n)
    z = x * 5.0 ** (-n)
    with np.errstate(over="ignore"):
        for _ in range(n - 1):
            z = z * (z + 5)
        eta = 0.5 / (z + 5)
    f = x ** (1 - al) * R
    h = u[1] - u[0]
    D = float(np.sum((f[1:] + f[:-1]) / 2) * h)
    g = f * eta
    H = float(np.sum((g[1:] + g[:-1]) / 2) * h)
    return D, H / D


def gap_n_float(al: float, n: int, nu_abs: float, om: float) -> float:
    """Floats: Y_n - Y_oo = rho^(2n) (om_n - om) - rho^n |nu| om - nu^2 om from the direct level-n quadrature."""
    rho = 5**al / 3
    om_n = level_quad_float(al, n)[1]
    return rho ** (2 * n) * (om_n - om) - rho**n * nu_abs * om - nu_abs**2 * om


# ----------------------------------------------------------------- sympy
def symbolic_checks() -> None:
    import sympy as sp
    s, z, r, nu, om, o1 = sp.symbols("s z rho nu om om1", real=True)
    be = (s + 6) / ((s + 2) * (s + 5))
    u1 = sp.simplify(be * (3 + s * (s + 5)) - 2 * s / (s + 2))
    check(sp.simplify(u1 - (s**2 + 7 * s + 9) / (s + 5)) == 0, "EO.1 (sympy): u_1 = beta(s)(3 + p(s)) - e(s) = (s^2 + 7s + 9)/(s + 5)")
    bz = be.subs(s, z)
    check(sp.simplify(1 - z / (z + 2) - sp.Rational(3, 2) * bz - 1 / (2 * (z + 5))) == 0,
          "EO.1 (sympy): eta_1(z) = 1 - z/(z+2) - (3/2) beta(z) = 1/(2(z + 5)) (BL, eta_0 = -1/2)")
    a1, a2 = (7 - sp.sqrt(13)) / 2, (7 + sp.sqrt(13)) / 2
    check(sp.expand((s + a1) * (s + a2) - (s**2 + 7 * s + 9)) == 0 and bool(a1 > 0), "EO.1 (sympy): s^2 + 7s + 9 = (s + a_1)(s + a_2), a_1 = (7 - sqrt13)/2 > 0")
    pf_j = (5 - a1) / (sp.sqrt(13) * (s + a1)) + (a2 - 5) / (sp.sqrt(13) * (s + a2))
    pf_h = (1 / (s + a1) - 1 / (s + a2)) / (2 * sp.sqrt(13))
    check(sp.simplify(pf_j - (s + 5) / (s**2 + 7 * s + 9)) == 0 and sp.simplify(pf_h - 1 / (2 * (s**2 + 7 * s + 9))) == 0
          and bool(5 - a1 > 0) and bool(a2 - 5 > 0),
          "EO.1 (sympy): r_1 and r_1 eta_1 split into 1/(s + a_i) with positive coefficients (5 - a_1)/sqrt13, (a_2 - 5)/sqrt13, +-1/(2 sqrt13)")
    Y1 = r * (r * (o1 - om) - (-nu * om))
    check(sp.simplify(Y1 - nu**2 * om - (r**2 * (o1 - om) + om * (r * nu - nu**2))) == 0,
          "EO.1 (sympy): Y_1 - Y_oo = rho^2 (om_1 - om) + om (rho nu - nu^2) = rho^2 (om_1 - om) - om |nu| (rho + |nu|) for nu < 0")


# ----------------------------------------------------------------- arb
def level1_arb(a_s):
    """EO.1 in arb: om_1 = 1 - delta_1^(1), J_1, D'_1 = J_1/rho (exact closed forms; ball only from rounding)."""
    arb, _ = ej._arb()
    al = arb(a_s)
    s13 = arb(13).sqrt()
    a1, a2 = (7 - s13) / 2, (7 + s13) / 2
    p1, p2 = a1 ** (-al), a2 ** (-al)
    c = arb.pi() / (s13 * (arb.pi() * al).sin())
    J1 = c * ((5 - a1) * p1 + (a2 - 5) * p2)
    H1 = c * (p1 - p2) / 2
    rho = arb(5) ** al / 3
    return {"om1": H1 / J1, "J1": J1, "D1": J1 / rho, "rho": rho}


def level1_gap(P, Q):
    """EO.1: Y_1 - Y_oo as a ball from EF.1's om, |nu| balls (P) and the exact om_1 (Q)."""
    om, nab, rho = P["om"][1], -P["nu"], P["rho"]
    return rho**2 * (Q["om1"] - om) - om * nab * (rho + nab)


def run_point(a_s, n_p, n_m):
    t1 = time.time()
    P = en.certify_point_en(a_s, n_p, n_m)
    Q = level1_arb(a_s)
    gap = level1_gap(P, Q)
    check(P["ok"], f"EF.1 routines valid at alpha = {a_s} ({n_p}/{n_m})")
    print(f"   alpha = {a_s} ({n_p}/{n_m}): D'_oo in {ej.fmt(P['Dp'], 8)}, |nu| in {ej.fmt(-P['nu'], 8)}, 1-delta_1 in "
          f"{ej.fmt(P['om'][1], 8)}; exact 1-delta_1^(1) = {float(Q['om1'].mid()):.10f}, D'_1 = {float(Q['D1'].mid()):.8f}; "
          f"Y_1 - Y_oo in [{float(gap.lower()):+.6e}, {float(gap.upper()):+.6e}] ({time.time() - t1:.1f}s)", flush=True)
    sgn = 1 if bool(gap > 0) else (-1 if bool(gap < 0) else 0)
    check(sgn == SIGN1[a_s], f"alpha = {a_s}, N = 1: Y_1 {'>' if SIGN1[a_s] > 0 else '<'} Y_oo (certified sign of the level-1 gap)")
    return P, Q, gap


def main(full: bool = False) -> int:
    t0 = time.time()
    print(f"Theorem EO: level 1 in closed form; full sign pattern at 0.88 and 0.89 (N = 1). Not thrust. "
          f"({'--full' if full else 'default'}, {em.WORKERS} workers)", flush=True)
    arb, ctx = ej._arb()
    ctx.prec = 128
    symbolic_checks()
    # ---- EO.4 part 1: closed form vs two independent float routes
    import mpmath as mp
    mp.mp.dps = 30
    for al in (0.7, 0.88, 0.91):
        J = mp.quad(lambda u: mp.exp(u * (1 - al)) * (mp.exp(u) + 5) / (mp.exp(2 * u) + 7 * mp.exp(u) + 9), [-mp.inf, 0, mp.inf])
        H = mp.quad(lambda u: mp.exp(u * (1 - al)) * mp.mpf(1) / 2 / (mp.exp(2 * u) + 7 * mp.exp(u) + 9), [-mp.inf, 0, mp.inf])
        uq = np.linspace(-400, 170, 114001)            # s = e^u; p(5s) stays below EJ's 1e150 clip
        sq = np.exp(uq)
        r1 = (5 / 3) * ej.rk_float(5 * sq, 1)          # r_1(s) = (5/3) R_1(5s), EJ's recursion
        Jt = float(np.trapezoid(sq ** (1 - al) * r1, uq))
        Ht = float(np.trapezoid(sq ** (1 - al) * r1 * 0.5 / (sq + 5), uq))
        Q = level1_arb(str(al))
        o1 = float(Q["om1"].mid())
        check(abs(float(H / J) - o1) < 1e-12 and abs(Ht / Jt - o1) < 1e-9 and abs(float(J) - float(Q["J1"].mid())) < 1e-10,
              f"EO.4: alpha = {al}: closed-form 1-delta_1^(1) = {o1:.12f} = mpmath {float(H / J):.12f} = trapezoid of EJ's r_1 {Ht / Jt:.12f}")
    # ---- EO.2
    res = {}
    for a_s, (n_p, n_m) in MAIN.items():
        res[a_s] = run_point(a_s, n_p, n_m)
    P88, Q88, g88 = res["0.88"]
    check(all(abs(float(x.lower()) - lo) < 1e-7 and abs(float(x.upper()) - hi) < 1e-7
              for x, (lo, hi) in ((P88["Dp"], EN2["Dp"]), (-P88["nu"], EN2["nu"]), (P88["om"][1], EN2["om"]))),
          "EO.2: 0.88 point brackets reproduce EN.2's to the quoted digits (same grid)")
    P89 = res["0.89"][0]
    check(EM3_NU[0] - 1e-9 <= float((-P89["nu"]).lower()) and float((-P89["nu"]).upper()) <= EM3_NU[1] + 1e-9,
          "EO.2: 0.89 |nu| bracket reproduces EM.3's (same grid)")
    for a_s, (last, onset) in PATTERN.items():
        if bool(res[a_s][2] < 0):
            print(f"     -> alpha = {a_s}, N = 1: Y_n < Y_oo for n = 1..{last}, Y_n > Y_oo for every n >= {onset}: "
                  f"exactly one sign change for all n >= 1, between {last} and {onset}", flush=True)
    # ---- attempts to break (0.88): EN.1 at n = 1 and EJ.2 (b) with the exact D'_1
    t1 = time.time()
    up = el.certify_upper_el("0.88", 1, 1600, Ns=(1,))
    U1 = up[1]["U"]
    check(up[1]["cond"], "EJ.1 alpha = 0.88: (3/rho) beta(p(1)) <= 1 at k0 = 1")
    low = em.certify_lower_em("0.88", (1,), (1,), 1600)
    L1, k1 = low["L"][1], low["kappa"][(1, 1)]
    d1lo, d1hi = en.dn_bounds(P88, 1, L1, U1)
    print(f"   (break) k0 = 1, 1600/unit: L(1) >= {float(L1.lower()):.6f}, U(1) <= {float(U1.upper()):.6f}, kappa_lo(1, 1) >= "
          f"{float(k1.lower()):.6f}; EN.1 at n = 1: D'_1 in [{float(d1lo.lower()):.6f}, {float(d1hi.upper()):.6f}], exact D'_1 = "
          f"{float(Q88['D1'].mid()):.8f} ({time.time() - t1:.1f}s)", flush=True)
    check(bool(L1 <= U1) and bool(d1lo < Q88["D1"]) and bool(Q88["D1"] < d1hi),
          "attempt to break: EN.1's certified D'_1 bracket (k0 = n = 1) contains EO.1's exact D'_1")
    b_old = en.lower_side(P88, L1, k1, 1, arb(P88["Dp"].upper()))
    b_ex = en.lower_side(P88, L1, k1, 1, arb(Q88["D1"].upper()))
    print(f"     EJ.2 (b) at n = 1, k0 = 1: {b_old:+.4g} with D'_1 <= D'_oo, {b_ex:+.4g} with the exact D'_1 "
          f"(EN.2 with k0 <= 1: +{EN_LEVEL1})", flush=True)
    check(b_ex < b_old, "attempt to break: the exact D'_1 (< D'_oo) lowers EJ.2 (b)'s braces bound at n = 1 (consistency)")
    print(f"     (OBSERVED) EJ.2 (b) at n = 1 stays undecided even with the exact D'_1 ({b_ex:+.4g} > 0: L(1), kappa_lo(1, 1) are too "
          f"loose at k0 = 1); the tail route does not reach level 1, EO.1's closed form does", flush=True)
    # ---- EO.3
    if full:
        for a_s in SIDE:
            res[a_s] = run_point(a_s, *SIDE_GRID)
        if bool(res["0.85"][2] < 0):
            print(f"     -> alpha = 0.85: Y_1 < Y_oo and (EJ.4) Y_n > Y_oo for every n >= {EJ4_ONSET85}: at least one sign change in [1, 6]",
                  flush=True)
    # ---- EO.4 part 2 (OBSERVED): direct level-n floats
    ef = ej._mods()[0]
    for a_s in ("0.7", "0.8", "0.85", "0.88", "0.89", "0.9", "0.905", "0.91"):
        al = float(a_s)
        c = ef.critical_floats(al, Ns=(1,))
        vals = [gap_n_float(al, n, -c["nu"], c["om"][1]) for n in range(1, 13)]
        print(f"   (OBSERVED) alpha = {a_s}: float Y_n - Y_oo, n = 1..12: " + " ".join(f"{v:+.2e}" for v in vals), flush=True)
        cf1 = level1_gap_float(al, -c["nu"], c["om"][1])
        check(abs(cf1 - vals[0]) < 1e-6, f"EO.4: alpha = {a_s}: closed-form float Y_1 - Y_oo {cf1:+.6e} = direct level-1 quadrature {vals[0]:+.6e}")
        check(np.sign(vals[0]) == SIGN1[a_s], f"EO.4: alpha = {a_s}: float sign of Y_1 - Y_oo matches EO.2 / EO.3")
        if a_s in PATTERN:
            last = PATTERN[a_s][0]
            check(all((v < 0) == (n <= last) for n, v in enumerate(vals, 1)),
                  f"EO.4: alpha = {a_s}: float signs for n = 1..12 match the certified pattern (negative exactly for n <= {last})")
        if a_s in res:
            g = res[a_s][2]
            check(float(g.lower()) - 1e-6 < vals[0] < float(g.upper()) + 1e-6, f"EO.4: direct float Y_1 - Y_oo at {a_s} inside the certified ball")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EO checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(full="--full" in sys.argv))
