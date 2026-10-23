#!/usr/bin/env python3
"""Theorem BM: the sharpness of Theorem BL holds at every alpha in (gamma, 1) for all large N (qualitative), with certified
thresholds N_0(alpha) at alpha = 0.7, 0.75, 0.85, 0.9, 0.95.

Setting and notation exactly as in Theorem BL (scripts/gasket_fractional_rate_sharpness.py): p(y) = y(y+5),
rho = 5^a/3, gamma = log3/log5, delta_N^(k) = <G_k>_(mu_k) the data-free corner deficit at the level-N corner neighbours,
delta_N = lim_k delta_N^(k); L(y0), J(y0), I_N, Q, I_0, D'_k as in BL.2-BL.4 (all stated there for every alpha in (gamma,1)).

Theorem BM (checked here; exit 1 on any failure):
  (1) qualitative, every alpha in (gamma, 1): fix any y0 > 0.  Then L(y0) > 0 (omega_lo > 1 for y > 0), J(y0) <= J_alpha
      < inf (BK.4: needs 2 gamma > alpha at 0 and alpha > 0 at inf), Q < inf (head K_S t^(1-gamma) t^(gamma-1-a) = K_S t^-a,
      alpha < 1; tail t^(gamma-1-a), alpha > gamma), I_0 < inf, D'_N >= D'_0 > 0, and g_N = 1 - eta_N(p(y0)) >= p(y0)/(p(y0)+2)
      > 0 for every N >= 1 (first term of BL.1's positive series).  So BL.4(ii)'s ceiling
        Delta_N <= (G Q rho^-N + (3/2)(3/5)^N I_0)/D'_0 -> 0   (rho > 1),
      and BL.2's bracket (g_N - Delta_N) L(y0) - Delta_N J(y0) is > 0 for all N >= some N_0(alpha).  Hence for every
      N >= N_0(alpha) and n >= N: delta_N^(n) strictly increases in n and c_N rho^-n <= delta_N - delta_N^(n) <= C rho^-n with
      c_N > 0, liminf_N c_N >= (p(y0)/(p(y0)+2)) L(y0)/(D'_inf (1 - rho^-1)) > 0.  So BK.4's order rho^-n is exact at
      EVERY alpha in (gamma, 1) (the sharpness half of Conjecture BK.6 / BL.6 for all alpha, at large N).
  (2) certified (python-flint arb, outward): explicit N_0 at alpha = 0.7 (y0 = 1/10), 0.75 (1/10), 0.85 (1/50), 0.9 (1/10),
      0.95 (1/10): direct ceiling BL.4(i) for N_0 <= N < N_u and the N-uniform ceiling BL.4(ii) at N_u (monotone in N, so
      it covers every N >= N_u).
Floats (consistency only): delta_6^(n) at alpha = 0.9 increases in n, and rho^n times the n-error sits in [c_6, ceiling].
Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import sympy as sp
from flint import arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasket_fractional_uniform_limit as bi  # noqa: E402
import gasket_fractional_corner_sharpness as bj  # noqa: E402
import gasket_fractional_uniform_rate as bk  # noqa: E402
import gasket_fractional_rate_sharpness as bl  # noqa: E402

FAILS: list[str] = []

# alpha -> (y0, N_0 (first direct N), N_u (uniform from here on))
CASES = {
    "0.7": ("0.1", 61, 62),
    "0.75": ("0.1", 16, 18),
    "0.85": ("0.02", 7, 10),
    "0.9": ("0.1", 6, 8),
    "0.95": ("0.1", 5, 8),
}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def symbolic_checks():
    a, y, t = sp.symbols("a y t", positive=True)
    gam = sp.log(3) / sp.log(5)
    check(bool(sp.N(2 * gam, 30) > 1) and bool(sp.N(gam, 30) < 1),
          "2 gamma = 2 log3/log5 > 1 > gamma: J_alpha's head exponents gamma+1-a, 2 gamma-a are > 0 for every alpha < 1 (BK.4)")
    Y = sp.symbols("Y", positive=True)
    B = 1 - Y * (Y + 4) / (3 * (Y + 2) * (Y + 5))
    check(sp.simplify(1 - B - Y * (Y + 4) / (3 * (Y + 2) * (Y + 5))) == 0 and sp.limit(B, Y, sp.oo) == sp.Rational(2, 3),
          "omega_lo = (1+y/5)^a / B, 0 < B < 1 for y > 0 (B -> 2/3 at inf): omega_lo > 1 on (0, inf), so L(y0) > 0 for every y0, alpha")
    z = sp.symbols("z", positive=True)
    check(sp.simplify(z / (z + 2) - (1 - 2 / (z + 2))) == 0,
          "j = 0 term of BL.1's positive series: 1 - eta_N(z) >= z/(z+2) for every N >= 1, so g_N >= p(y0)/(p(y0)+2) > 0")
    # Q: head K_S t^(1-gamma) * t^(gamma-1-a) = K_S t^-a; tail int_1^inf t^(gamma-1-a) = 1/(a-gamma)
    g = sp.symbols("g", positive=True)
    check(sp.simplify(t ** (1 - g) * t ** (g - 1 - a) - t ** (-a)) == 0,
          "Q's head integrand is K_S t^-alpha (integrable at 0 for alpha < 1); tail t^(gamma-1-alpha) integrable at inf for alpha > gamma")


def certify(al: str, y0: str, N0: int, Nu: int, G: arb):
    ctx.prec = 96
    a = arb(al)
    gam = bl.gam_arb()
    rho = arb(5) ** a / 3
    check(bool(a > gam) and bool(a < 1), f"alpha = {al} lies in (gamma, 1)")
    Jtot, _, _, _ = bk.J_bound(al, G)
    _, D0, ec, _ = bk.constants(al, Jtot)
    Dup = arb((D0 + ec).upper())
    L0, _ = bl.L_lower(al, "1e-10")
    Ly, _ = bl.L_lower(al, y0)
    Jy, _ = bl.J_upper(al, G, y0)
    ctx.prec = 96
    Dlow = lambda N: arb((D0 + L0 * (1 - rho ** (-N)) / (1 - 1 / rho)).lower())  # noqa: E731
    p0 = arb(y0) * (arb(y0) + 5)
    rows = []
    ok = bool(Ly.lower() > 0)
    for N in range(N0, Nu):
        IN, _ = bl.I_direct(al, G, N)
        ctx.prec = 96
        Dl = arb((IN / Dlow(N)).upper())
        tot, _ = bl.one_minus_eta_series(p0, N)
        g = arb(tot.lower())
        br = (g - Dl) * Ly - Dl * Jy
        cN = arb((br / (Dup * (1 - 1 / rho))).lower())
        ok = ok and bool(Dl < g) and br.lower() > 0
        rows.append((N, float(Dl.upper()), float(g.lower()), float(cN.lower())))
    Q = bl.Q_upper(al)
    ctx.prec = 96
    I0 = bj.I_closed(al, arb(1), G)
    Du = arb(((G * Q * rho ** (-Nu) + arb(3) / 2 * (arb(3) / 5) ** Nu * I0) / Dlow(Nu)).upper())
    _, trunc = bl.one_minus_eta_series(p0, min(Nu, 40))  # truncated positive series: a floor for every N >= Nu
    gu = arb(trunc.lower())
    bru = (gu - Du) * Ly - Du * Jy
    cu = arb((bru / (Dup * (1 - 1 / rho))).lower())
    ok = ok and bool(Du < gu) and bru.lower() > 0
    _, tr = bl.one_minus_eta_series(p0, 40)  # floor for g_N, N >= 40
    casym = arb((arb(tr.lower()) * Ly / (Dup * (1 - 1 / rho))).lower())
    up = arb((arb(3) / 2 * ec / D0).upper())
    print(f"   alpha={al}: rho={float(rho.mid()):.6f}, J_a<={float(Jtot.upper()):.6f}, D'_0={float(D0.mid()):.6f}, D'_inf<={float(Dup.upper()):.6f}, "
          f"L(0)>={float(L0.lower()):.6f}, y0={y0}: L>={float(Ly.lower()):.6f}, J<={float(Jy.upper()):.6f}, Q<={float(Q.upper()):.6f}, I_0<={float(I0.upper()):.6f}")
    for N, Dl, g, cN in rows:
        print(f"      N={N} (direct): Delta_N<={Dl:.6f} < g>={g:.6f}, c_N>={cN:.4e}")
    print(f"      N>={Nu} (uniform): Delta_N<={float(Du.upper()):.6f} < g>={float(gu.lower()):.6f}, c_N>={float(cu.lower()):.4e}; "
          f"liminf_N c_N>={float(casym.lower()):.4f}; ceiling (3/2)(D'_inf-D'_0)/D'_0 coefficient <= {float(up.upper()):.4f}")
    check(ok, f"certified alpha={al}, y0={y0}: bracket > 0 for N={N0}..{Nu - 1} (direct) and every N>={Nu} (uniform): "
              f"delta_N^(n) strictly increases in n>=N and delta_N-delta_N^(n) >= c_N rho^-n for every N>={N0}")
    return {"N0": N0, "Nu": Nu, "rows": rows, "cu": float(cu.lower()), "casym": float(casym.lower()), "up": float(up.upper()),
            "c_first": rows[0][3] if rows else float(cu.lower())}


def float_deficits(al="0.9", N=6, ns=(6, 7, 8, 10, 12, 16, 20), nref=60):
    ref = bi.corner_deficit(al=al, n=nref, Ns=(N,), dps=20)[0][1]
    vals = [(n, bi.corner_deficit(al=al, n=n, Ns=(N,), dps=20)[0][1]) for n in ns]
    return vals, ref


def main() -> int:
    t0 = time.time()
    print("Theorem BM: BL's sharpness at every alpha in (gamma,1) for large N; certified N_0 at 0.7, 0.75, 0.85, 0.9, 0.95. Not thrust.")
    bl.FAILS.clear()
    bl.symbolic_checks()
    bl.exact_checks()
    check(not bl.FAILS, "Theorem BL's sympy and exact-rational checks (BL.1, alpha-free) re-run")
    symbolic_checks()
    G0, _, _ = bi.sup_g()
    G = arb("0.6031")
    okb, _ = bi.bootstrap(G0, G)
    check(okb and G0.upper() < G.lower(), "BI.4 re-run: sup s^(1-gamma) r_inf <= G = 0.6031")
    res = {}
    for al, (y0, N0, Nu) in CASES.items():
        res[al] = certify(al, y0, N0, Nu, G)
    summary = ", ".join(f"N_0({al})={r['N0']}" for al, r in res.items())
    print("   thresholds: " + summary)
    vals, ref = float_deficits()
    rf = 5**0.9 / 3
    inc = all(vals[i][1] < vals[i + 1][1] for i in range(len(vals) - 1)) and vals[-1][1] < ref
    errs = [(n, float(ref - v) * rf**n) for n, v in vals]
    print("   floats alpha=0.9, N=6: rho^n (delta^(60) - delta^(n)) = " + ", ".join(f"{e:.4f} (n={n})" for n, e in errs))
    r9 = res["0.9"]
    check(inc and all(r9["c_first"] <= e <= r9["up"] for n, e in errs),
          "floats: delta_6^(n)(0.9) increases in n and rho^n times the n-error lies in [c_6, BK's ceiling] (consistency)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BM checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
