#!/usr/bin/env python3
"""Theorem EK: the exact sign-change level at alpha = 0.905 for N = 2, and a narrower bracket for N = 3.  EJ kept "bracketed at
0.905 (N = 2, 3)": the last level with Y_n < Y_oo was only located in [82, 83] (N = 2) and [101, 103] (N = 3).  No new lemma
here: EJ.1-EJ.2's certificates (scripts/gasket_fractional_corner_two_sided_rate.py, imported unchanged) are run on finer
grids, and the lower side additionally at k0 = 14.

  EK.1 (arb, 128 bits; --full; EJ.1-EJ.2 as stated): at alpha = 0.905 (EF.1 cells n_p = 2000, n_m = 8000; EJ.1 upper k0 = 14 on
       n = 8000; EJ.1 lower k0 <= 12 on n = 2000, and k0 in {12, 14} on n = 4000, hi = 14):
         N = 1: Y_n < Y_oo for n = 1..67, Y_n > Y_oo for every n >= 68 (reproduces EJ.3);
         N = 2: Y_n < Y_oo for n = 2..83, Y_n > Y_oo for every n >= 84: exactly one sign change, between 83 and 84 (EJ: [82, 83]);
         N = 3: Y_n < Y_oo for n = 3..102, Y_n > Y_oo for every n >= 104: last negative level in [102, 103] (EJ: [101, 103]);
                n = 103 is not decided by this method (KEPT LIMIT).
  Default (no --full): only the lower-side certificates (n_p = 2000, n_m = 8000 point, lower k0 in {12, 14} on n = 4000) at the
       three levels N = 2, n = 82, 83 and N = 3, n = 102, plus the (OBSERVED) sign of the method's bound at N = 3, n = 103.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_two_sided_rate as ej  # noqa: E402

FAILS: list[str] = []
A_S = "0.905"
N_P, N_M = 2000, 8000            # EF.1 point cells (EJ --full: 1000, 4000)
N_UX, K0_UP = 8000, 14           # EJ.1 upper side (EJ --full: 4000, k0 = 14)
N_LA, HI_A = 2000, 12            # EJ.1 lower side, k0 in EJ.K0_LO (EJ --full: 1000, hi 12)
N_LB, HI_B, K0_B = 4000, 14, (12, 14)
EJ_BRACKET = {1: (67, 67), 2: (82, 83), 3: (101, 103)}   # EJ.3 at 0.905 (last negative level)
EJ_DP, EJ_NU = (5.86394, 5.86483), (0.414526, 0.415382)   # EJ.3 at 0.905 (--full grid)
EJ_EPS = (0.729805, 0.741544)                              # EJ.3 at 0.905: [L(12), U(14)]
EXPECT = {1: (67, 68), 2: (83, 84), 3: (102, 104)}        # EK.1 (--full): last certified negative level, onset n1


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def lower(n, hi, k0s):
    old = ej.K0_LO
    ej.K0_LO = tuple(k0s)
    try:
        return ej.certify_lower(A_S, n, hi)
    finally:
        ej.K0_LO = old


def bound_at(P, low, N, n):
    """EJ.2's lower-side upper bound on the braces of EB.1 at level n (as in EJ.below_upper), smallest over the k0 in low; < 0 certifies."""
    arb, _ = ej._arb()
    rho, g = P["rho"], P["g"]
    Dhi, nhi = arb(P["Dp"].upper()), arb((-P["nu"]).upper())
    oml, omh = arb(P["om"][N].lower()), arb(P["om"][N].upper())
    best = None
    for k0 in sorted(low["L"]):
        if k0 > n:
            continue
        X = nhi**3 * g**n - low["L"][k0] / Dhi
        den = 1 - nhi * rho ** (-n)
        val = (oml * X if bool(X < 0) else omh / den * X) - 3 * low["kappa"][(N, k0)] / (2 * Dhi)
        if best is None or float(val.upper()) < float(best[1].upper()):
            best = (k0, val)
    return best


def main(full: bool = False) -> int:
    t0 = time.time()
    print(f"Theorem EK: exact sign-change level at alpha = {A_S} (EJ.1-EJ.2 on finer grids). Not thrust. ({'--full' if full else 'default'})")
    arb, ctx = ej._arb()
    ctx.prec = 128
    P = ej.certify_point(A_S, N_P, N_M)
    print(f"   point (n_p = {N_P}, n_m = {N_M}): D'_oo in {ej.fmt(P['Dp'])}, |nu| in {ej.fmt(-P['nu'])}, 1-delta_N in "
          + ", ".join(f"N={N}: {ej.fmt(P['om'][N])}" for N in ej.NS) + f" ({time.time() - t0:.1f}s)", flush=True)
    check(P["ok"], f"EF.1 routines valid at alpha = {A_S}")
    check(EJ_DP[0] <= float(P["Dp"].lower()) and float(P["Dp"].upper()) <= EJ_DP[1] and EJ_NU[0] <= float((-P["nu"]).lower())
          and float((-P["nu"]).upper()) <= EJ_NU[1], f"consistency: D'_oo and |nu| brackets nested inside EJ.3's {EJ_DP}, {EJ_NU}")
    t1 = time.time()
    lowB = lower(N_LB, HI_B, K0_B)
    print(f"   lower (n = {N_LB}, hi = {HI_B}, k0 in {K0_B}): " + ", ".join(f"L({k}) >= {float(v.lower()):.6f}" for k, v in lowB["L"].items())
          + "; " + ", ".join(f"kappa_lo({N},{k}) >= {float(v.lower()):.6f}" for (N, k), v in lowB["kappa"].items() if k == 14)
          + f" ({time.time() - t1:.1f}s)", flush=True)
    for N, n in ((2, 82), (2, 83), (3, 102)):
        r = ej.below_upper(P, lowB, N, n)
        print(f"     N = {N}, n = {n}: " + (f"braces <= {float(r[1].upper()):.4g} < 0 (k0 = {r[0]}) -> Y_n < Y_oo" if r else "not certified"))
        check(r is not None, f"EK.1 (EJ.2 lower side) alpha = {A_S}, N = {N}: Y_{n} < Y_oo")
    b = bound_at(P, lowB, 3, 103)
    print(f"   (OBSERVED) N = 3, n = 103: smallest EJ.2 lower-side bound {float(b[1].upper()):+.4g} (k0 = {b[0]}) > 0 -> not decided by this method")
    check(ej.below_upper(P, lowB, 3, 103) is None and bool(b[1] > 0), "KEPT LIMIT: N = 3, n = 103 not certified negative (bound > 0)")
    if full:
        t1 = time.time()
        up = ej.certify_upper(P["al"], (K0_UP,), n=N_UX)
        row = up[K0_UP]
        print(f"   upper (k0 = {K0_UP}, n = {N_UX}): c_hi <= {float(row['c_hi'].upper()):.6f}, U({K0_UP}) <= {float(row['U'].upper()):.6f}, "
              + ", ".join(f"kappa_hi({N}) <= {float(row['kappa'][N]['hi'].upper()):.6f}" for N in ej.NS) + f" ({time.time() - t1:.1f}s)", flush=True)
        check(row["cond"], f"EJ.1 alpha = {A_S}: (3/rho) beta(p^k0(1)) <= 1")
        t1 = time.time()
        lowA = lower(N_LA, HI_A, ej.K0_LO)
        print(f"   lower (n = {N_LA}, hi = {HI_A}, k0 in {ej.K0_LO}): L(12) >= {float(lowA['L'][12].lower()):.6f} ({time.time() - t1:.1f}s)", flush=True)
        check(bool(lowA["L"][12] < row["U"]) and bool(lowB["L"][14] < row["U"])
              and all(bool(lowB["kappa"][(N, 14)] < row["kappa"][N]["hi"]) for N in ej.NS)
              and EJ_EPS[0] <= float(lowA["L"][12].lower()) and float(row["U"].upper()) <= EJ_EPS[1],
              "two-sided: L(12), L(14) < U(14), [L(12), U(14)] nested inside EJ.3's [0.729805, 0.741544] and kappa_lo(N, 14) < kappa_hi(N, 14) (no contradiction)")
        for N in ej.NS:
            on = ej.onset_upper(P, up, N)
            n1 = on[0] if on else None
            win = [n for n in range(N, n1 or 200) if ej.below_upper(P, lowA, N, n)]
            extra = [n for n in range(N, n1 or 200) if n not in win and ej.below_upper(P, lowB, N, n)]
            neg = sorted(win + extra)
            contig = bool(neg) and neg == list(range(N, neg[-1] + 1))
            last = neg[-1] if neg else None
            exact = n1 is not None and last == n1 - 1
            print(f"     N = {N}: Y_n > Y_oo for EVERY n >= {n1}; Y_n < Y_oo for n in [{neg[0] if neg else '-'}, {last}] (k0 <= 12 grid: up to "
                  f"{win[-1] if win else '-'}, k0 in {{12, 14}} adds {extra})"
                  + (f" -> exactly one sign change, between {last} and {n1}" if exact else f" -> last negative level in [{last}, {n1 - 1 if n1 else '?'}]"))
            lo_, hi_ = EJ_BRACKET[N]
            check(n1 is not None and contig and (last, n1) == EXPECT[N],
                  f"EK.1 alpha = {A_S}, N = {N}: negative for n = {N}..{EXPECT[N][0]}, positive for every n >= {EXPECT[N][1]}")
            check(lo_ <= last <= hi_ and n1 - 1 <= hi_,
                  f"consistency with EJ.3: bracket [{last}, {n1 - 1}] inside EJ's [{lo_}, {hi_}]")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EK checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(full="--full" in sys.argv))
