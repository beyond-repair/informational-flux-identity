#!/usr/bin/env python3
r"""Theorem DN: DM's kept limit "DK.4's level-sum bound is not carried to C_2, because Theta_2(f_0, Lt f_0) is not computed".

Setting (GASKET.md, Theorems CE, CF, DC-DM): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
C_a = {f > 0 : log f a-Lipschitz in s}, Theta_a its Hilbert metric (CF's definition), now a = 2 (DL/DM's cone);
f_0 = F0 = Lt^4 G (DC.2), u = F0, v = Lt F0, R = v/u.  DI's machinery (DI.1's exact quotient identity, DI.2's arb cell
tables, DI.4's cell-pair pass and floor) is re-run unchanged with a = 2 instead of a = 100.

Theorem DN (checked here; exit 1 on any failure):
  DN.1 [certified; DI.2 re-run, arb] (log F0)' in (0.035291, 0.77334), (log Lt F0)' in (0.17682, 0.79038), so
       f_0 in C_0.77334 and Lt f_0 in C_0.79038, both inside C_2 with slope < a: CF.3's hypothesis f in C_b, b < a holds
       for a = 2, and every Phi_u(x, y) = e^(2|x-y|) u(y) - u(x) > 0 (DI.1 applies).
  DN.2 [certified; DI.4 with a = 2: arb cell tables, outward float cell pairs; floor: arb point limits]
       0.14444 < Theta_2(f_0, Lt f_0) <= Theta_N' < 0.14473 (16384 cells; DN_CELLS=65536: < 0.14452).
       DI.3's signed factor lemma gives only 0.9448 at a = 2 (cone term (a+beta)/(a-alpha) is large for small a).
  DN.3 [proved from DN.2 with CF.3's proof (DH.5/DI.5 form), DL.5's theta_L < 0.0108954 (n = 5), DM.5's theta_7 <
       0.0011028 (n = 7) and theta_9 < 0.00011123 (n = 9) used as stated; DL.1: Lt C_2 subset C_2]
       for every k >= 4, |log(s_(k+1)/(r s_k))| <= Theta_N' min(theta_L^floor((k-4)/5), theta_7^floor((k-4)/7),
       theta_9^floor((k-4)/9)); combined with DK.4 (Theta_I' min(theta_G^.., theta_K^..)) by taking the minimum.
       Consistency: CC's s_4..s_13 below the combined bound.
  DN.4 [certified floor -> kept failure] every CF.3-form bound with the C_2 constant Theta_2(f_0, Lt f_0) is > 0.14444
       = 2.22 x DI.4's Theta_I' < 0.064982, so for k = 4..8 (before the first C_2 block step) the C_2 route is worse than
       DI.5/DK.4 and the minimum is needed.
Kept limits: see GASKET.md Theorem DN.  No Psi_info, Stage 2 still blocked.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time

from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_f0_coupled as di          # noqa: E402

FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A_NUM = 2
NC = int(os.environ.get("DN_CELLS", "16384"))
NFP = int(os.environ.get("DN_FLOOR_PTS", "2049"))
# contraction pairs on C_2 used as stated (DL.5, DM.5; DM's scripts are not re-run here)
PAIRS = ((5, "0.0108954"), (7, "0.0011028"), (9, "0.00011123"))
# DK.4's level-sum statement on C_100 (used as stated): Theta_I' (DI.4), theta_G (DG.3, n = 8), theta_K (DK.4, n = 9)
TI_DI, THETA_G, THETA_K = "0.064982", "0.014183", "0.0035000"
TN_MAX = {16384: "0.14473", 65536: "0.14452"}.get(NC, "0.14473")
TFLOOR_MIN = "0.14444"
TS_MIN = "0.9"                     # DI.3 at a = 2 (comparison only; checked to be worse than DN.2)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def main_f0(T):
    ctx.prec = PREC
    a = arb(A_NUM)
    di.NC = NC
    t0 = time.time()
    out = di.cell_tables(T)
    if out is None:
        check(False, "DN.1: DI.2's cell tables failed")
        return False
    tab, h = out
    print(f"   DN.1: {NC} cells ({time.time() - t0:.1f}s)", flush=True)
    c0, d0 = float(tab[:, 2].min()), float(tab[:, 3].max())
    c1, d1 = float(tab[:, 4].min()), float(tab[:, 5].max())
    Rlo, Rhi = float(tab[:, 0].min()), float(tab[:, 1].max())
    logR = arb((arb(Rhi) / arb(Rlo)).log().upper())
    check(bool(arb(c0) > arb("0.035291")) and bool(arb(d0) < arb("0.77334"))
          and bool(arb(c1) > arb("0.17682")) and bool(arb(d1) < arb("0.79038")) and d0 < A_NUM and d1 < A_NUM,
          f"DN.1: (log F0)' in [{c0:.8f}, {d0:.8f}], (log Lt F0)' in [{c1:.8f}, {d1:.8f}] (DI.2 re-run): f_0 in C_0.77334, "
          f"Lt f_0 in C_0.79038, both inside C_2 with slope < a = 2 (CF.3's hypothesis; every Phi_u > 0, DI.1 applies)")
    KK = di.signed_factor(a, arb(c0), arb(d0), arb(c1), arb(d1))
    TS = arb((KK[0].log() + KK[1].log() + logR).upper()) if KK is not None else None
    t0 = time.time()
    beta, alpha, jb, ja, upm, dnm = di.coupled(tab, h, a)
    TN = arb((arb(beta) / arb(alpha)).log().upper())
    print(f"   DN.2: cell pairs ({time.time() - t0:.1f}s); beta <= {beta:.12f} (cell {jb}), alpha >= {alpha:.12f} (cell {ja}), "
          f"sup R = {Rhi:.12f}, inf R = {Rlo:.12f}, excess at the extremal cells: up {upm:.3e}, down {dnm:.3e}", flush=True)
    t0 = time.time()
    fl = di.floor_points(T, a, m=NFP)
    check(bool(TN < arb(TN_MAX)) and bool(fl > arb(TFLOOR_MIN)) and bool(fl <= TN)
          and (TS is None or bool(TS > arb(TS_MIN))),
          f"DN.2: Theta_2(f_0, Lt f_0) <= Theta_N' = log(beta/alpha) = {TN.str(8)} < {TN_MAX} (DI.1/DI.4 with a = 2, {NC} "
          f"cells); certified floor (arb, {NFP} points, limits x -> y+-) {fl.str(8)} > {TFLOOR_MIN} "
          f"({time.time() - t0:.1f}s); Theta_N' - floor <= {(TN - fl).str(4)}; log(sup R/inf R) <= {logR.str(8)}; "
          f"DI.3's signed factor lemma at a = 2: {TS.str(6) if TS is not None else 'hypotheses fail'}")
    RESULT.update(c0=c0, d0=d0, c1=c1, d1=d1, logR=logR, TS=TS, TN=TN, floor=fl, beta=beta, alpha=alpha)
    return True


def bounds(k, TN):
    """(C_2 bound, DK.4 bound) for |log(s_(k+1)/(r s_k))|, k >= 4, as arb upper values."""
    m = k - 4
    c2 = min((arb(th) ** (m // n)).upper() for n, th in PAIRS)
    c100 = min((arb(THETA_G) ** (m // 8)).upper(), (arb(THETA_K) ** (m // 9)).upper())
    return arb((TN * arb(c2)).upper()), arb((arb(TI_DI) * arb(c100)).upper())


def consequences(T):
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    import gasket_fractional_profile_pure_point as cc
    bk = di.bk
    TN, fl = RESULT["TN"], RESULT["floor"]
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    rl = arb(r_lo.numerator) / r_lo.denominator
    rh = arb(r_hi.numerator) / r_hi.denominator
    found, sums, counts, okb = cc.period_atoms(bk.KMAX, verbose=False)
    ctx.prec = PREC
    r = bk._hull(rl, rh) if hasattr(bk, "_hull") else (rl + rh) / 2
    okd, worst = okb, arb(0)
    for k in range(4, bk.KMAX):
        # |log(s_(k+1)/(r s_k))| with r anywhere in [r_lo, r_hi]: take the worse end
        dev = arb(max(abs((sums[k + 1] / (rl * sums[k])).log()).upper(), abs((sums[k + 1] / (rh * sums[k])).log()).upper()))
        b2, b100 = bounds(k, TN)
        bmin = arb(min(b2.upper(), b100.upper()))
        okd = okd and bool(dev < bmin)
        worst = arb(max(worst.upper(), (dev / bmin).upper()))
    facs = {}
    for k in (9, 13, 22, 49, 103):
        b2, b100 = bounds(k, TN)
        facs[k] = b100 / b2
    rate9 = arb(PAIRS[2][1]) ** (arb(1) / 9)
    check(okd and all(bool(facs[k] > 1) for k in (9, 13, 22, 49, 103)) and bool(rate9 < arb("0.363657")),
          "DN.3: for every k >= 4, |log(s_(k+1)/(r s_k))| <= Theta_N' min(theta_L^floor((k-4)/5), theta_7^floor((k-4)/7), "
          "theta_9^floor((k-4)/9)) (CF.3 on (C_2, Lt^n), n = 5, 7, 9; DL.5/DM.5's thetas as stated); DK.4's bound divided by "
          + ", ".join(f"{facs[k].str(4)} (k = {k})" for k in facs)
          + f"; asymptotic per-level rate theta_9^(1/9) <= {rate9.str(7)} < 0.363657 (DK.4: 0.53349). Consistency: on CC's "
          f"s_4..s_{bk.KMAX} every deviation is below min(DN.3, DK.4) (max ratio {worst.str(4)})")
    check(bool(fl > arb(TI_DI)) and bool(fl / arb(TI_DI) > arb("2.22")),
          f"DN.4 (kept failure): every CF.3-form constant on C_2 for f_0 is Theta_2(f_0, Lt f_0) > {fl.str(7)} = "
          f"{(fl / arb(TI_DI)).str(4)} x DI.4's Theta_I' (< {TI_DI}); so for k = 4..8 (no C_2 block step yet) DN.3 is worse "
          f"than DI.5/DK.4 and the minimum is needed")
    RESULT.update(facs=facs, worst=worst, rate9=rate9)


def main() -> int:
    t0 = time.time()
    print("Theorem DN: Theta_2(f_0, Lt f_0) on the small cone C_2 and the level-sum bound on C_2. Not thrust.", flush=True)
    T = di.taylor()
    if not main_f0(T):
        print("FAILED")
        return 1
    consequences(T)
    print(f"elapsed {time.time() - t0:.1f}s")
    bk = di.bk
    nf = len(FAILS) + len(di.FAILS) + len(bk.FAILS) + len(bk.co.FAILS) + len(bk.cn.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DN checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
