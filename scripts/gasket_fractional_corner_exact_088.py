#!/usr/bin/env python3
"""Theorem EN: the exact sign-change level at alpha = 0.88 (N = 1).  EM kept (EM.4) "alpha = 0.88: EL.1 upper at k0 = 9 gives onset 10
(Phi(9) <= -0.001656), EM.1 lower-side bound at 9 is +0.0007145: level 9 undecided by both sides, bracket [8, 9] unchanged".  The
lower side (EJ.2 (b)) loses at small n in its denominator, not in its integrals: it replaces D'_n by D'_oo (ED.1: D'_n increases to
D'_oo), and at n = 9, alpha = 0.88, D'_9 is about 2.9 % below D'_oo (rho^-9 |C'| ~ 0.157 against D'_oo ~ 5.47).

  EN.1 (THEOREM; every alpha in (gamma, 1), k0 >= 1, n >= k0; ED.1 / AZ (1) exact, EJ.1 and EM.1 (EG.1) as stated).
         D'_oo + rho^-n C' + 3^-n L(k0)  <=  D'_n  <=  D'_oo + rho^-n C' + 3^-n U(k0),   C' = -rho Mt_oo/(rho - 1).
       Proof: ED.1 gives D'_(k+1) - D'_k = rho^-k Mt_k, so D'_oo - D'_n = sum_(k>=n) rho^-k Mt_k = rho^-n rho Mt_oo/(rho-1)
       - sum_(k>=n) rho^-k (Mt_oo - Mt_k) = -rho^-n C' - rho^-n eps_n (ED.1's eps_n), i.e. D'_n = D'_oo + rho^-n (C' + eps_n) (AZ (1));
       and L(k0) <= (3/rho)^n eps_n <= U(k0) for every n >= k0 (EJ.1 upper, EG.1/EM.1 lower), so rho^-n eps_n in [3^-n L, 3^-n U].
       Hence EJ.2 (b)'s lower-side criterion stays valid with D'_oo's upper end replaced by
         D_n^hi = D'_oo^hi + rho^-n C'^hi + 3^-n U(k0),   C'^hi = -rho Mt_oo^lo/(rho - 1):
       with E_n = (3/rho)^n eps_n >= L(k0) > 0 and 0 < D'_n <= D_n^hi, -E_n/D'_n <= -L/D_n^hi and -3^n eta_n/D'_n <= -1.5 kappa_lo/D_n^hi,
       and the sign split on X = |nu|_hi^3 g^n - L/D_n^hi (om_lo if X < 0, om_hi/(1 - |nu|_hi rho^-n) otherwise) is EJ.2's.  Each n is
       checked on its own (no monotone-in-n claim).  The upper side (EJ.2 (a), Phi increasing in n) is NOT changed.
  EN.2 (THEOREM, certified; arb, 128 bits; EF.1 point cells 4000/32000 (EM's certify_mt_par), EL.1 upper and EM.1 lower at k0 = 9 on
       1600 cells/unit, EM.1 lower k0 in {1, 2, 3, 4, 6, 8} as in EM): at alpha = 0.88, N = 1, EN.1's lower side at n = 9 is < 0, so
       Y_9 < Y_oo; with EJ.2 (b) for n = 2..8 (reproduces EJ.4's window) and EL.1's onset (Y_n > Y_oo for every n >= 10, EJ.4 / EM.4),
       Y_n - Y_oo changes sign exactly once for n >= 2, between levels 9 and 10 (EM.4: [8, 9]).  Level 1 is not decided (kept limit).
  EN.3 (OBSERVED, floats): D'_9 by trapezoid of x^-a R_9 (exact affine recursion, numpy) lies inside EN.1's certified [D_9^lo, D_9^hi];
       EE's float limits sit inside [L(9), U(9)], [kappa_lo, kappa_hi]; the float crossing is at 10.
  EN.4 (OBSERVED; EM's docstring corrected): at 0.905, k0 = 12, EM.1's L(12) on 200, 800, 3200 cells/unit (logged here; EM's docstring
       had said L(14) and quoted 800/3200 values that were never logged).

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

FAILS: list[str] = []
A88 = "0.88"
K0 = 9
K0_LO = (1, 2, 3, 4, 6, 8, 9)
CFG = {"n_p": 4000, "n_m": 32000, "n": 1600}
EJ4_WINDOW = (2, 8)        # EJ.4: Y_n < Y_oo for n = 2..8 at 0.88 (N = 1)
ONSET = 10                 # EJ.4 / EM.4: Y_n > Y_oo for every n >= 10 (EL.1 upper, k0 = 9)
EXPECT = (9, 10)           # EN.2: last negative level, onset
EM4 = {"lower9": 0.0007145, "phi9": -0.001656}   # EM.4's undecided bounds at n = 9 (consistency)


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
def rk_scaled_float(x, k: int):
    """Floats: R_k(x) = (3/5)^k r_k(5^-k x) via w_j = u_j/s_j (u_j = 1/r_j on the orbit s_j = p^j(5^-k x)):
    w_k = 1 + 3/s_k, w_j = w_(j+1) (s_j+6)/(s_j+2) - 2/(s_j+2) (the affine recursion u_j = beta(s_j) u_(j+1) - e(s_j) divided by
    s_j); bounded for large s (s_j may overflow to inf harmlessly), unlike u_j itself."""
    t = np.asarray(x, float) * 5.0 ** (-k)
    with np.errstate(over="ignore", invalid="ignore"):
        ss = [t]
        for _ in range(k):
            ss.append(ss[-1] * (ss[-1] + 5))
        w = 1 + 3 / ss[k]
        for j in range(k - 1, -1, -1):
            w = w * (1 + 4 / (ss[j] + 2)) - 2 / (ss[j] + 2)
    return 0.6**k / (w * t)


def dprime_float(al: float, k: int, lo: int = -150, hi: int = 150, m: int = 100) -> float:
    """Floats: D'_k = rho^-k J_k = int_0^oo x^-a R_k(x) dx (x = 5^k s), trapezoid in log x on [5^lo, 5^hi]."""
    u = np.linspace(lo, hi, (hi - lo) * m + 1) * np.log(5.0)
    x = np.exp(u)
    f = x ** (1 - al) * rk_scaled_float(x, k)
    return float(np.sum((f[1:] + f[:-1]) / 2) * (u[1] - u[0]))


def dn_bounds_float(Dlo, Dhi, Mlo, Mhi, L, U, al: float, n: int):
    """Floats: EN.1's [D_n^lo, D_n^hi] from brackets of D'_oo, Mt_oo and L(k0), U(k0)."""
    rho = 5**al / 3
    return (Dlo - rho * Mhi / (rho - 1) * rho ** (-n) + 3.0 ** (-n) * L,
            Dhi - rho * Mlo / (rho - 1) * rho ** (-n) + 3.0 ** (-n) * U)


# ----------------------------------------------------------------- arb
def certify_point_en(a_s, n_p, n_m):
    """EM's certify_point_em with Mt_oo kept (C' = -rho Mt_oo/(rho - 1) is needed separately by EN.1)."""
    arb, _ = ej._arb()
    ef, _ = ej._mods()
    al = arb(a_s)
    rho = arb(5) ** al / 3
    pi = ef.certify_profile_integrals(al, (1,), n=n_p)
    Mt, mi = em.certify_mt_par(a_s, n_m)
    Dp = pi["Dp"]
    return {"al": al, "rho": rho, "g": 3 / rho**3, "Dp": Dp, "Mt": Mt, "nu": -rho * Mt / ((rho - 1) * Dp),
            "om": {1: pi["h"][1] / Dp}, "ok": mi["ok"] and pi["ok_tail"][1]}


def dn_bounds(P, n, L, U):
    """EN.1 (arb): [D_n^lo, D_n^hi]."""
    arb, _ = ej._arb()
    rho = P["rho"]
    Dlo, Dhi = arb(P["Dp"].lower()), arb(P["Dp"].upper())
    Mlo, Mhi = arb(P["Mt"].lower()), arb(P["Mt"].upper())
    rn = rho ** (-n)
    return Dlo - rho * Mhi / (rho - 1) * rn + L / arb(3) ** n, Dhi - rho * Mlo / (rho - 1) * rn + U / arb(3) ** n


def lower_side(P, L, klo, n, Dh, N=1):
    """EJ.2 (b)'s braces bound at level n with the given upper bound Dh on D'_n (EN.1: D_n^hi; EJ.2: D'_oo^hi); < 0 certifies."""
    arb, _ = ej._arb()
    rho, g = P["rho"], P["g"]
    nhi = arb((-P["nu"]).upper())
    oml, omh = arb(P["om"][N].lower()), arb(P["om"][N].upper())
    den = 1 - nhi * rho ** (-n)
    if not bool(den > 0):
        return None
    X = nhi**3 * g**n - L / Dh
    val = (oml * X if bool(X < 0) else omh / den * X) - 3 * klo / (2 * Dh)
    return float(val.upper())


def l12_0905(ns=(200, 800, 3200)):
    """EN.4: EM.1's certified L(12) at 0.905 (right-node cell bound) on the given grids."""
    return {n: float(em.certify_lower_em("0.905", (12,), (1,), n)["L"][12].lower()) for n in ns}


def main() -> int:
    t0 = time.time()
    print(f"Theorem EN: the exact sign-change level at alpha = 0.88 (N = 1), with EN.1's two-sided D'_n. Not thrust. "
          f"({em.WORKERS} workers)", flush=True)
    arb, ctx = ej._arb()
    ctx.prec = 128
    t1 = time.time()
    P = certify_point_en(A88, CFG["n_p"], CFG["n_m"])
    check(P["ok"], f"EF.1 routines valid at alpha = {A88}")
    print(f"   point 0.88 ({CFG['n_p']}/{CFG['n_m']}): D'_oo in {ej.fmt(P['Dp'], 8)}, Mt_oo in {ej.fmt(P['Mt'], 8)}, |nu| in "
          f"{ej.fmt(-P['nu'], 8)}, 1-delta_1 in {ej.fmt(P['om'][1], 8)} ({time.time() - t1:.1f}s)", flush=True)
    t1 = time.time()
    up = el.certify_upper_el(A88, K0, CFG["n"], Ns=(1,))
    r = up[K0]
    U, khi = r["U"], r["kappa"][1]["hi"]
    check(r["cond"], f"EJ.1 alpha = {A88}: (3/rho) beta(p^k0(1)) <= 1 at k0 = {K0}")
    on = ej.onset_upper(P, up, 1)
    print(f"   upper (EL.1, k0 = {K0}, n = {CFG['n']}): U({K0}) <= {float(U.upper()):.8f}, kappa_hi(1, {K0}) <= {float(khi.upper()):.8f}; "
          f"Phi({on[0]}) >= {float(on[2].lower()):.4g} > 0, Phi({on[0] - 1}) <= {float(on[3].upper()):.4g} ({time.time() - t1:.1f}s)",
          flush=True)
    check(on[0] == ONSET, f"EJ.2 (a) (unchanged): Y_n > Y_oo for every n >= {ONSET} at alpha = {A88} (EJ.4 / EM.4)")
    t1 = time.time()
    low = em.certify_lower_em(A88, K0_LO, (1,), CFG["n"])
    L9, k9 = low["L"][K0], low["kappa"][(1, K0)]
    print(f"   lower (EM.1, n = {CFG['n']}): L(9) >= {float(L9.lower()):.8f}, kappa_lo(1, 9) >= {float(k9.lower()):.8f} "
          f"({time.time() - t1:.1f}s)", flush=True)
    check(bool(L9 <= U) and bool(k9 <= khi), "two-sided: L(9) <= U(9) and kappa_lo(1, 9) <= kappa_hi(1, 9)")
    # ---- EJ.2 (b) as before (D'_n <= D'_oo), n = 1..9, best k0 <= n
    Dhi_old = arb(P["Dp"].upper())
    old = {}
    for n in range(1, K0 + 1):
        vs = [lower_side(P, low["L"][k], low["kappa"][(1, k)], n, Dhi_old) for k in K0_LO if k <= n]
        old[n] = min(v for v in vs if v is not None)
        print(f"     EJ.2 (b) (D'_n <= D'_oo) at n = {n}: {old[n]:+.4g}" + ("  (Y_n < Y_oo)" if old[n] < 0 else "  (undecided)"), flush=True)
    win = [n for n in range(1, K0 + 1) if old[n] < 0]
    check(win == list(range(EJ4_WINDOW[0], EJ4_WINDOW[1] + 1)),
          f"EJ.2 (b) reproduces EJ.4's window: Y_n < Y_oo for n = {EJ4_WINDOW[0]}..{EJ4_WINDOW[1]}; n = 1 and n = 9 undecided by it")
    check(0 < old[9] <= EM4["lower9"] + 1e-7, f"EJ.2 (b) at n = 9: {old[9]:+.4g} (> 0; k0 = 9 added here; EM.4: +{EM4['lower9']} with k0 <= 8)")
    # ---- EN.1 at n = 9 (k0 = 9)
    D9lo, D9hi = dn_bounds(P, K0, L9, U)
    print(f"   EN.1 at n = 9: D'_9 in [{float(D9lo.lower()):.8f}, {float(D9hi.upper()):.8f}] (EJ.2 used D'_9 <= {float(Dhi_old.upper()):.8f})",
          flush=True)
    check(bool(D9hi < Dhi_old) and bool(D9lo < D9hi), "EN.1: D_9^lo < D_9^hi < D'_oo^hi")
    b9 = lower_side(P, L9, k9, K0, D9hi)
    print(f"     EN.1 lower side at n = 9: {b9:+.4g}" + ("  (Y_9 < Y_oo)" if b9 < 0 else "  (undecided)"), flush=True)
    check(b9 < 0, f"EN.2 alpha = {A88}, N = 1: Y_9 < Y_oo (EN.1 lower side at n = 9 < 0; EM.4: +{EM4['lower9']})")
    if b9 < 0 and win == list(range(2, 9)) and on[0] == ONSET:
        print(f"     -> exactly one sign change at alpha = {A88} (N = 1, n >= 2), between levels {EXPECT[0]} and {EXPECT[1]} (EM.4: [8, 9])")
    # ---- EN.3 (OBSERVED, floats)
    d9 = dprime_float(0.88, K0)
    cf = ej._mods()[0].critical_floats(0.88, Ns=(1,))
    print(f"   (OBSERVED) float D'_9 = {d9:.8f} (trapezoid, 100/unit on [5^-150, 5^150]); float D'_oo = {cf['Dp']:.8f}", flush=True)
    check(float(D9lo.lower()) < d9 < float(D9hi.upper()), "EN.3: float D'_9 lies inside EN.1's certified [D_9^lo, D_9^hi]")
    check(float(L9.lower()) < cf["e"] < float(U.upper()) and float(k9.lower()) < cf["kappa"][1] < float(khi.upper()),
          f"EN.3: EE's float limits e* = {cf['e']:.7f}, kappa_1 = {cf['kappa'][1]:.7f} inside [L(9), U(9)], [kappa_lo, kappa_hi]")
    check(ej.crossing_float(0.88, cf["e"], cf["kappa"][1], N=1)[0] == EXPECT[1], "EN.3: float limit crossing at 10")
    # ---- EN.4 (OBSERVED): EM's docstring numbers
    t1 = time.time()
    l12 = l12_0905()
    print("   (OBSERVED, EN.4) 0.905, k0 = 12: EM.1's L(12) >= " + ", ".join(f"{v:.6f}" for v in l12.values())
          + f" on {', '.join(str(k) for k in l12)} cells/unit ({time.time() - t1:.1f}s)", flush=True)
    vals = list(l12.values())
    check(abs(vals[0] - 0.735421) < 2e-6 and vals[0] < vals[1] < vals[2] < 0.738670,
          "EN.4: L(12) at 200/unit reproduces EM's logged 0.735421; increasing in the grid, below EE's float e*(0.905) = 0.738670")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EN checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
