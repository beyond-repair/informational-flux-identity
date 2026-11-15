#!/usr/bin/env python3
"""Theorem EI: an explicit onset level below 4gamma/3.  EH (and EG before it) kept the limit "one-sided only: no every-n upper
bound, so no explicit level at which Y_n > Y_oo starts below 4gamma/3 (EC.4's n >= 60 at 0.8 stays the only explicit one), and no
exact position of the sign change in EH.4".  The upper side does not need a two-sided rate for c_k: EC.2's crude every-k bound
|eps_k| <= B_alpha (rho/3)^k (k >= 6, every alpha in (gamma, 1)) and DZ.1's eta_k <= etabar_k suffice, because below 4gamma/3 the
term |nu|^3 (3/rho^3)^n grows geometrically.

  EI.1 (THEOREM; every alpha in (gamma, 4gamma/3), N >= 1; EB.1 exact, EC.2, DZ.1, ED.1 as stated).  For every k >= N,
         Y_k - Y_oo = (1-delta_N)/(1+nu rho^-k) [ |nu|^3 rho^-k - rho^k eps_k / D'_k ] - rho^(2k) eta_k / D'_k     (EB.1, sympy)
       = (rho^2/3)^k { (1-delta_N)/(1+nu rho^-k) [ |nu|^3 g^k - (3/rho)^k eps_k / D'_k ] - 3^k eta_k / D'_k },  g = 3/rho^3 > 1.
       With D'_k = D'_oo + rho^-k (C' + eps_k) >= D'_oo + C' rho^-k (ED.1: eps_k >= 0), (3/rho)^k |eps_k| <= B_alpha (EC.2, k >= 6),
       0 <= 3^k eta_k <= H_k = 3^k etabar_k = 1.5 a_N + b_N (3/(5 rho))^k / (1 - 1/(5 rho)) (DZ.1), 0 < 1 + nu rho^-k <= 1:
         Phi_N(k) = om_lo (nu_lo^3 g^k - B/D_k) - H_k/D_k > 0,  D_k = D_lo - rho M_hi rho^-k/(rho-1) > 0
       implies Y_k > Y_oo.  Phi_N is strictly increasing in k (g > 1, D_k increasing, H_k decreasing as 3/(5 rho) < 1), so
       Phi_N(n1) > 0 at one n1 >= max(6, N) gives Y_n > Y_oo for EVERY n >= n1.
  EI.2 (arb, 128 bits; EF.1's cell routines for D'_oo, Mt_oo, h_N at fixed alpha; EC.2's B_alpha; DZ.1's etabar): the smallest
       such n1(alpha, N) at alpha = 0.7, 0.8, 0.85, 0.88, 0.89, 0.9, 0.905, 0.91 and N = 1, 2, 3, 6, 10.
  EI.3 (THEOREM, with EH.4): at 0.9, 0.905, 0.91, for N <= 8, 21, 926, the LAST level with Y_n < Y_oo exists and lies in
       [8, n1 - 1], [21, n1 - 1], [926, n1 - 1]; at 0.8 the explicit onset n >= 60 of EC.4 drops to n1 (13..18).
  EI.4 (OBSERVED, floats): the limit-substituted crossing level (EE's e*, kappa_N in place of the finite-n coefficients) lies
       inside every certified bracket; the bracket is loose because B_alpha ~ 10^3 e*.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

FAILS: list[str] = []
RESULT: dict = {}
GAMMA = math.log(3) / math.log(5)
NS = (1, 2, 3, 6, 10)
ALPHAS = ("0.7", "0.8", "0.85", "0.88", "0.89", "0.9", "0.905", "0.91")
EH4_WINDOW_END = {"0.9": 8, "0.905": 21, "0.91": 926}          # EH.4 (--full): Y_n < Y_oo for n in [max(N, 2 or 3), end]
G_BI4 = 0.6031


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def _mods():
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import gasket_fractional_corner_critical_sign as ef
    return ef


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
def B_float(al: float) -> float:
    """EC.2's B_alpha in floats (same formula as EC's `constants`): B = K I_6/(rho'^6 (1-rho')) + 2 T_6/rho'^6."""
    gm = GAMMA
    G, y0 = G_BI4, 0.2
    ce = 2 * G * (5 + y0) ** gm / 6
    me = 1 - y0**gm * ce
    K1 = y0 ** (1 - gm) * (16 + 3 * y0) / 30 + 2 * ce / me
    K = math.exp(al / 20 + K1 / 2) - 1 + 0.1
    rp = 5**al / 9
    ap = al - gm
    p6 = 1
    for _ in range(6):
        p6 = p6 * (p6 + 5)
    I6 = math.exp(ap / 20 + K1 / 2) * G / (2 * gm - al) + 3e10 * 2 * 1.2**6 / (math.sqrt(p6) * (2**5 - 5 + ap))
    T6 = 3e10 * 2 * 1.2**6 * math.exp(-(2**4) * math.log(6)) / (2**5 - 5 + al)
    return K * I6 / (rp**6 * (1 - rp)) + 2 * T6 / rp**6


def H_float(al: float, N: int, k: int) -> float:
    """3^k etabar_k (DZ.1) in floats: 1.5 a_N + b_N (3/(5 rho))^k / (1 - 1/(5 rho))."""
    gm, G = GAMMA, G_BI4
    m0 = 2 / 3 - 1 / (21 + 6 * math.sqrt(10))
    K1 = ((0.2 + 0.14) * 0.2 ** (1 - gm) + 1.8 * G * 5 ** (gm - 1)) / m0
    rho = 5**al / 3
    Lam = 1 / (2 * gm - al) + 1 / (1 + al - 2 * gm)
    Q = (1 / (al + 1) + 1 / (5 * (al + 0.5))) / m0
    aN = K1 * G * 5 ** (-gm) * Lam * (2 * 5**N) ** (2 * gm - al)
    bN = 2 * 5**N * Q
    return 1.5 * aN + bN * (3 / (5 * rho)) ** k / (1 - 1 / (5 * rho))


def onset_float(al: float, Ns=NS):
    """Floats (OBSERVED): EI.1's criterion with EF's float D'_oo, nu, 1-delta_N: the smallest n1 >= max(6, N) with Phi_N(n1) > 0;
    and the limit-substituted crossing n_as (EE's e*, kappa_N): first n with om (|nu|^3 g^n - e*/D') - 1.5 kappa_N/D' > 0."""
    ef = _mods()
    c = ef.critical_floats(al, Ns=Ns)
    rho = 5**al / 3
    g = 3 / rho**3
    B = B_float(al)
    nu, Dp, Mt = -c["nu"], c["Dp"], c["Mt"]
    out = {"alpha": al, "g": g, "B": B, "nu": nu, "Dp": Dp, "e": c["e"], "n1": {}, "nas": {}, "om": {}}
    for N in Ns:
        om = c["om"][N]
        out["om"][N] = om

        def phi(k):
            Dk = Dp - rho * Mt / (rho - 1) * rho ** (-k)
            return om * (nu**3 * g**k - B / Dk) - H_float(al, N, k) / Dk if Dk > 0 else -1.0

        lo, hi = max(6, N), max(6, N)
        while phi(hi) <= 0:
            lo, hi = hi, 2 * hi
        while lo < hi:
            m = (lo + hi) // 2
            lo, hi = (lo, m) if phi(m) > 0 else (m + 1, hi)
        out["n1"][N] = hi
        n = 1
        while om * (nu**3 * g**n - c["e"] / Dp) - 1.5 * c["kappa"][N] / Dp <= 0 and n < 10**6:
            n += 1
        out["nas"][N] = n
    return out


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    import sympy as sp
    k = sp.symbols("k", positive=True, integer=True)
    nu, om, eps, eta, D, rho, x, a = sp.symbols("nu omega epsilon eta D rho x alpha", real=True)
    Cd = -nu * om
    X = Cd / (1 + nu * x) - (om * eps / (1 + nu * x) + eta / x) / D        # EB.1, x = rho^-k, rho^k eta_k = eta/x
    Y = (X - Cd) / x
    Yinf = nu**2 * om
    rhs = om / (1 + nu * x) * (-nu**3 * x - eps / (x * D)) - eta / (x**2 * D)
    check(sp.simplify(Y - Yinf - rhs) == 0, "sympy: EB.1 => Y_k - Y_oo = (1-delta_N)/(1+nu x)(|nu|^3 x - rho^k eps_k/D'_k) - rho^2k eta_k/D'_k, x = rho^-k (nu < 0)")
    r, T, e = sp.symbols("r T e", positive=True)
    okg = all(sp.simplify((r**2 / 3) ** j * ((3 / r**3) ** j * T - (3 / r) ** j * e) - (T * r ** (-j) - r**j * e)) == 0 for j in range(1, 9))
    check(okg, "sympy: (rho^2/3)^k [(3/rho^3)^k T - (3/rho)^k e] = T rho^-k - rho^k e (the g^k form; k = 1..8 symbolic rho)")
    al = sp.symbols("alpha", positive=True)
    gm = sp.log(3) / sp.log(5)
    check(sp.simplify(sp.log(sp.Rational(81, 1)) / (3 * sp.log(5)) - 4 * gm / 3) == 0, "sympy: 3/rho^3 = 81/5^(3 alpha) > 1 iff alpha < log_5(81)/3 = 4gamma/3")
    check(bool(sp.N(9 - 5 ** (1 + gm)) < 0), "sympy: 3/(5 rho) = 9/5^(1+alpha) < 9/15 < 1 for alpha > gamma (H_k decreasing)")
    check(sp.simplify(sp.diff(rho ** (-k), k) + sp.log(rho) * rho ** (-k)) == 0, "sympy: d/dk rho^-k = -log(rho) rho^-k < 0 (D_k increasing)")


# ----------------------------------------------------------------- (B) arb
def _arb():
    from flint import arb, ctx
    return arb, ctx


def ebar3_arb(al, N: int, k: int, cst, G):
    """3^k etabar_k (DZ.1) as an arb upper bound (DZ's own `ebar`, times 3^k)."""
    import gasket_fractional_corner_limit as dz
    arb, _ = _arb()
    return dz.ebar(al, N, k, cst, G) * arb(3) ** k


def certify_alpha(a: str, n_p: int, n_m: int, cst, G):
    arb, ctx = _arb()
    ef = _mods()
    import gasket_fractional_corner_remainder_sharp as ec
    c = ec.constants(a)                                      # EC.2 at 256 bits
    ctx.prec = 128
    B = arb(c["B"].upper())
    al = arb(a)
    rho = arb(5) ** al / 3
    g = 3 / rho**3
    pi = ef.certify_profile_integrals(al, NS, n=n_p)
    Mt, mi = ef.certify_mt(al, n=n_m)
    Dp = pi["Dp"]
    nu = -rho * Mt / ((rho - 1) * Dp)
    Dlo, Mhi = arb(Dp.lower()), arb(Mt.upper())
    nlo, nhi = arb((-nu).lower()), arb((-nu).upper())
    ok = c["ok"] and mi["ok"] and bool(g > 1) and bool(nlo > 0) and bool(nhi < rho**6)
    row = {"alpha": a, "B": B, "Dp": Dp, "Mt": Mt, "nu": nu, "g": g, "n1": {}, "om": {}, "ok": ok, "Phi": {}}
    for N in NS:
        om = arb((pi["h"][N] / Dp).lower())
        ok_o = pi["ok_tail"][N] and bool(om > 0)

        def phi(k):
            Dk = Dlo - rho * Mhi / (rho - 1) * rho ** (-k)
            if not bool(Dk > 0):
                return None, Dk
            return om * (nlo**3 * g**k - B / Dk) - ebar3_arb(al, N, k, cst, G) / Dk, Dk

        k0 = max(6, N)
        lo, hi = k0, k0
        while not bool((phi(hi)[0] or arb(-1)) > 0):
            lo, hi = hi, 2 * hi
        while lo < hi:
            m = (lo + hi) // 2
            v = phi(m)[0]
            lo, hi = (lo, m) if (v is not None and bool(v > 0)) else (m + 1, hi)
        v1, D1 = phi(hi)
        v0 = phi(hi - 1)[0] if hi > k0 else None
        minimal = hi == k0 or v0 is None or not bool(v0 > 0)
        row["n1"][N], row["om"][N], row["Phi"][N] = hi, om, v1
        row["ok"] = row["ok"] and ok_o and v1 is not None and bool(v1 > 0) and bool(D1 > 0) and minimal
    return row


def fmt(v, k=6):
    return f"[{float(v.lower()):.{k}g}, {float(v.upper()):.{k}g}]"


def main(full: bool = False) -> int:
    t0 = time.time()
    print("Theorem EI: an explicit onset level of Y_n > Y_oo below 4gamma/3 (EB.1 + EC.2 + DZ.1 + ED.1). Not thrust.")
    symbolic_checks()
    arb, ctx = _arb()
    import gasket_fractional_corner_limit as dz
    G = arb("0.6031")
    cst = dz.dz1_constants(G)                                 # DZ.1's m_0, K_1 (re-checked)
    ctx.prec = 128
    n_p, n_m = (1000, 4000) if full else (200, 500)
    rows = []
    for a in ALPHAS:
        t1 = time.time()
        r = certify_alpha(a, n_p, n_m, cst, G)
        rows.append(r)
        print(f"   alpha = {a}: B <= {float(r['B'].upper()):.4e}, g = 3/rho^3 in {fmt(r['g'], 7)}, D'_oo in {fmt(r['Dp'])}, "
              f"Mt_oo in {fmt(r['Mt'])}, |nu| in {fmt(-r['nu'])} ({time.time() - t1:.1f}s)")
        print("     " + ", ".join(f"N={N}: 1-delta_N >= {float(r['om'][N].lower()):.5f}, n1 = {r['n1'][N]}" for N in NS))
        check(r["ok"], f"EI.2 alpha = {a}: Phi_N(n1) > 0 with D_n1 > 0, n1 minimal for the route, for N = 1,2,3,6,10 "
              f"(n1 = {[r['n1'][N] for N in NS]}): Y_n > Y_oo for EVERY n >= n1")
    RESULT["rows"] = rows
    by = {r["alpha"]: r for r in rows}
    # ---- consistency with earlier theorems
    check(bool(by["0.8"]["B"] < arb("1.94e3")) and bool(by["0.9"]["B"] < arb("1.03e3")), "consistency: B_0.8 <= 1.94e3, B_0.9 <= 1.03e3 (EC.2's printed values)")
    nu8 = -by["0.8"]["nu"]
    check(bool(nu8.lower() <= 0.7274454) and bool(nu8.upper() >= 0.7274454), "consistency: |nu_0.8| bracket contains AZ (3)'s 0.7274454")
    nu9 = -by["0.9"]["nu"]
    check(bool(nu9.lower() <= 0.432350) and bool(nu9.upper() >= 0.432534), "consistency: |nu_0.9| bracket contains EF.4's [0.432350, 0.432534]")
    check(all(by["0.8"]["n1"][N] <= 60 for N in NS), f"EI.3: at 0.8 the onset drops below EC.4's 60: n1 = {[by['0.8']['n1'][N] for N in NS]}")
    ok3 = True
    for a, w in EH4_WINDOW_END.items():
        for N in NS:
            if N <= w:
                ok3 = ok3 and by[a]["n1"][N] > w
                print(f"   EI.3 alpha = {a}, N = {N}: last level with Y_n < Y_oo in [{w}, {by[a]['n1'][N] - 1}]")
    check(ok3, "EI.3: every n1 exceeds EH.4's window end (no contradiction; the last negative level is bracketed)")
    check(all(rows[i]["n1"][N] <= rows[i + 1]["n1"][N] for i in range(len(rows) - 1) for N in NS),
          "EI.2: n1 nondecreasing in alpha on the grid (g -> 1 as alpha -> 4gamma/3)")
    # ---- EI.4 floats
    for a in ALPHAS:
        s = onset_float(float(a))
        r = by[a]
        print(f"   floats alpha = {a}: n1 (floats) = {[s['n1'][N] for N in NS]}, limit-substituted crossing n_as = {[s['nas'][N] for N in NS]}")
        check(all(s["n1"][N] <= r["n1"][N] for N in NS), f"EI.4 alpha = {a}: float n1 <= certified n1 (certified route only loses)")
        check(all(s["nas"][N] < r["n1"][N] for N in NS), f"EI.4 (OBSERVED) alpha = {a}: limit-substituted crossing below the certified onset")
        if a in EH4_WINDOW_END:
            w = EH4_WINDOW_END[a]
            check(all(s["nas"][N] > w for N in NS if N <= w), f"EI.4 (OBSERVED) alpha = {a}: limit-substituted crossing after EH.4's window end {w}")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EI checks passed")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    sys.exit(main(full="--full" in sys.argv))
