#!/usr/bin/env python3
r"""Theorem DW: the cone-free direct bound carried to k = 12..16 (DV's rejected strict claim DV.3' "B_V(k) < min(B_DO(k),
B_U(k)) for every k >= 4", false exactly at k = 13, 16, and DV's kept limit "direct term only for k <= 11").

Setting (GASKET.md, Theorems DO, DT, DU, DV): s in X = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
F_0 = f_0 = Lt^4 G (DC.2; CO.1/CO.2's Taylor data re-run unchanged), F_m = Lt^m F_0, R_m = F_(m+1)/F_m,
osc log R_m = log(sup R_m / inf R_m), s_(4+m) = c0 F_m(x0), so s_(k+1)/s_k = R_(k-4)(x0) (DI.5's consistency).

Theorem DW (checked here; exit 1 on any failure):
  (Default: rows m = 6..9 only, DW.2/DW.3 then use the stated uppers TARGET_UP for m = 10..12, certified by DW_FULL=1.)
  DW.1 [certified; arb jets of order 4 through m + 1 <= 13 layers of Lt (DT/DV's Fjet, unchanged), a FOURTH-order centred
       form of R_m = F_(m+1)/F_m on adaptively bisected cells: R_m(x) in sum_(i<4) R^(i)(c)/i! (x-c)^i + [R^(4)(J)/4!](x-c)^4,
       c the cell midpoint (point jets), J the cell hull, intersected with the direct quotient F_(m+1)(J)/F_m(J); floors:
       arb point values at DW_FLOOR_PTS points, endpoints included] for m = 8..12: F_m > 0 on X and two-sided bounds for
       osc log R_m.  Rows m = 6, 7 are recomputed the same way as a cross-check against DV.1 (the intervals must meet).
  DW.2 [proved from DW.1 and the inclusion r in [inf R_m, sup R_m] of CF.3/DO.3's proof, as in DV.2; no cone membership
       is needed] |log(s_(k+1)/(r s_k))| <= osc log R_(k-4) for k = 12..16; by DO.4 (floor of every CF.3-form bound at
       k = 4 + m is osc log R_m) it is within the factor upper/lower of the best bound of that form.
  DW.3 [exact rationals, outward-rounded inputs; DV's constants as stated] B_W(k) = min(B_V(k), osc log R_(k-4) [12<=k<=16])
       satisfies B_W(k) < B_V(k) for k = 12..16 and B_W(k) < min(B_DO(k), B_U(k)) for EVERY k >= 4 (checked directly for
       k < 390, beyond by the exact 8-periodicity of G_U as in DV.3): DV.3' (the strict form) holds for B_W.
Kept limits: see GASKET.md Theorem DW.  No Psi_info, Stage 2 still blocked.  Stage 1 untouched; no W, 0.08, 0.23 or
30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction as Fr

from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_direct_constants as dv   # noqa: E402  (DV: jets, cells, constants as stated)

dt = dv.dt
di, bk = dv.di, dv.bk
FAILS: list[str] = []
RESULT: dict = {}
PREC = int(os.environ.get("DW_PREC", "128"))
FULL = os.environ.get("DW_FULL", "") not in ("", "0")
# new direct rows m = M0..M1 (k = 12..16); the default computes only m = 8, 9 and uses the stated (DW_FULL-certified)
# upper bounds TARGET_UP for m = 10..12 in DW.2/DW.3, as DV.3 uses DV's stated constants
M0, M1 = 8, int(os.environ.get("DW_MMAX", "12" if FULL else "9"))
MN = 12                                                  # last direct row entering B_W
MX = (6, 7)                                              # cross-check rows (DV.1)
NC0 = int(os.environ.get("DW_CELLS0", "64"))
DEPTH = int(os.environ.get("DW_DEPTH", "8"))
TOL = float(os.environ.get("DW_TOL", "0.002"))
NFP = int(os.environ.get("DW_FLOOR_PTS", "33"))
WORKERS = int(os.environ.get("DW_WORKERS", "4"))
# certified targets (checked, exit 1 if not met), m = 8..12
TARGET_UP = {8: '8.80e-8', 9: '1.519e-8', 10: '2.621e-9', 11: '4.523e-10', 12: '7.805e-11'}
TARGET_LO = {8: '8.70e-8', 9: '1.500e-8', 10: '2.588e-9', 11: '4.465e-10', 12: '7.703e-11'}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


_TJ = None


def cell4(TJ, m, a, b):
    """(lo, hi) of R_m on the cell [a, b]: fourth-order centred form intersected with the direct quotient; None if
    positivity of F_m on the cell is not certified."""
    J = bk._hull(a, b)
    mp = arb(J.mid())
    dx = J - mp
    uJ, vJ = dv._jet(TJ, m, J), dv._jet(TJ, m + 1, J)
    if not (bool(uJ[0] > 0) and bool(vJ[0] > 0)):
        return None
    uM, vM = dv._jet(TJ, m, mp), dv._jet(TJ, m + 1, mp)
    rJ = dt.jmul(vJ, dt.jrec(uJ))
    rM = dt.jmul(vM, dt.jrec(uM))
    d2 = dx * dx                                   # products, not powers (arb's ** on a ball around 0 is not used)
    rc = rM[0] + rM[1] * dx + rM[2] * d2 + rM[3] * d2 * dx + rJ[4] * d2 * d2
    rd = vJ[0] / uJ[0]
    return max(rc.lower(), rd.lower()), min(rc.upper(), rd.upper())


def _job(args):
    m, p, q = args
    ctx.prec = PREC
    r = cell4(_TJ, m, dv._x(p), dv._x(q))
    return (p, q, None if r is None else (dv._fr(r[0]), dv._fr(r[1])))


def floors(TJ, m):
    """arb point values of R_m at NFP points (endpoints included): certified (max of lowers, min of uppers)."""
    ctx.prec = PREC
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    hi = lo = None
    for k in range(NFP):
        s = s0 if k == 0 else (s1 if k == NFP - 1 else arb((s0 + (s1 - s0) * k / (NFP - 1)).mid()))
        r = dv._jet(TJ, m + 1, s)[0] / dv._jet(TJ, m, s)[0]
        hi = r.lower() if hi is None else max(hi, r.lower())
        lo = r.upper() if lo is None else min(lo, r.upper())
    return hi, lo


def row(TJ, m, pool):
    ctx.prec = PREC
    t0 = time.time()
    bf, af = floors(TJ, m)
    low = arb((arb(bf) / arb(af)).log().lower())
    eps = Fr(float(low.lower()) * TOL / 2)
    capb, capa = dv._fr(bf) * (1 + eps), dv._fr(af) * (1 - eps)
    step = 2 ** dv.DM // NC0
    wave = [(m, i * step, (i + 1) * step) for i in range(NC0)]
    acc, nsplit, d = [], 0, 0
    while wave:
        res = pool.map(_job, wave, chunksize=1) if pool else list(map(_job, wave))
        nxt = []
        for p, q, r in res:
            if r is not None and ((r[1] <= capb and r[0] >= capa) or d >= DEPTH):
                acc.append((p, q, r))
            elif d < DEPTH + 4 and q - p >= 2:
                c = (p + q) // 2
                nxt += [(m, p, c), (m, c, q)]
                nsplit += 1
            else:
                check(False, f"DW.1: positivity of F_m not certified on cell ({p}, {q}) (m = {m})")
                return None
        wave, d = nxt, d + 1
    acc.sort()
    tiles = acc[0][0] == 0 and acc[-1][1] == 2 ** dv.DM and all(acc[i][1] == acc[i + 1][0] for i in range(len(acc) - 1))
    if not tiles:
        check(False, f"DW.1: accepted cells do not tile X (m = {m})")
        return None
    b = max(r[1] for _, _, r in acc)
    a = min(r[0] for _, _, r in acc)
    up = arb((dv._arb(b) / dv._arb(a)).log().upper())
    maxd = max(step.bit_length() - (q - p).bit_length() for p, q, _ in acc)
    dt_ = time.time() - t0
    print(f"   m = {m}: osc log R_m in [{low.str(8, radius=False)}, {up.str(8, radius=False)}]; {len(acc)} cells "
          f"(max depth {maxd}, {nsplit} splits), {NFP} floor points ({dt_:.0f}s)", flush=True)
    return dict(up=up, low=low, ncell=len(acc), maxd=maxd, secs=dt_)


def tables(T):
    global _TJ
    TJ = _TJ = dv.TaylorJet4P(T)
    import multiprocessing as mp
    pool = mp.get_context("fork").Pool(WORKERS) if WORKERS > 1 else None
    tab = {}
    try:
        for m in list(MX) + list(range(M0, M1 + 1)):
            r = row(TJ, m, pool)
            if r is None:
                return False
            tab[m] = r
    finally:
        if pool:
            pool.close()
            pool.join()
    RESULT["tab"] = tab
    return True


def consequences():
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    import gasket_fractional_profile_pure_point as cc
    tab = RESULT["tab"]
    _up, knapF, f = dv._up, dv.knapF, (lambda x: arb(x.numerator) / arb(x.denominator))
    # ---- DW.1: rows, targets, cross-check with DV.1 (m = 6, 7)
    xok = all(bool(tab[m]["up"] >= arb(dv.TARGET_OSC_LO[m])) and bool(tab[m]["low"] <= arb(dv.TARGET_OSC_UP[m]))
              and bool(tab[m]["low"] <= tab[m]["up"]) for m in MX)
    check(xok, "DW.1 cross-check (fourth-order form vs DV.1's second-order rows): " + "; ".join(
        f"m={m}: [{tab[m]['low'].str(6, radius=False)}, {tab[m]['up'].str(6, radius=False)}] meets DV.1's "
        f"({dv.TARGET_OSC_LO[m]}, {dv.TARGET_OSC_UP[m]})" for m in MX))
    new = list(range(M0, M1 + 1))
    tg = all(bool(tab[m]["low"] <= tab[m]["up"]) and bool(tab[m]["up"] < arb(TARGET_UP[m]))
             and bool(tab[m]["low"] > arb(TARGET_LO[m])) for m in new if m in TARGET_UP)
    check(tg, "DW.1: F_m > 0 on X and osc log R_m two-sided (arb): " + "; ".join(
        f"m={m}: [{tab[m]['low'].str(6, radius=False)}, {tab[m]['up'].str(6, radius=False)}]" for m in new)
        + " (targets met)")
    ms = sorted(tab)
    rat = {m: (tab[m + 1]["low"] / tab[m]["up"], tab[m + 1]["up"] / tab[m]["low"]) for m in ms if m + 1 in tab}
    RESULT["rat"] = rat
    print("     consecutive ratios osc_(m+1)/osc_m in: " + "; ".join(
        f"m={m}: [{a.str(5, radius=False)}, {b.str(5, radius=False)}]" for m, (a, b) in rat.items()), flush=True)
    # ---- DW.2 / DW.3: exact rationals with DV's stated constants (outward)
    oscV = {m: _up(dv.TARGET_OSC_UP[m]) for m in range(8)}
    thK = {m: _up(dv.TARGET_K_UP[m]) for m in range(8)}
    for m in range(6):
        thK[m] = min(thK[m], _up(dv.DT_UPPER[m]))
    oscW = {m: (_up(tab[m]["up"]) if m in tab else _up(TARGET_UP[m])) for m in range(M0, MN + 1)}
    upW = {m: (tab[m]["up"] if m in tab else arb(TARGET_UP[m])) for m in range(M0, MN + 1)}
    th = {n: _up((arb(D) / 4).tanh()) for n, D in dv.DBAR_U.items()}
    KM = dv.KM
    G = knapF(sorted(th.items()), KM)
    GD = {a: knapF([(n, _up(t)) for n, t in dv.DO_PAIRS[a]], 20) for a in (2, 100)}
    thK_dt = {m: _up(dv.DT_UPPER[m]) for m in range(6)}

    def B_DO(k):
        return min(_up(dv.DO_UPPER[a][m]) * GD[a][k - 4 - m] for m in range(0, min(5, k - 4) + 1) for a in (2, 100))

    def B_U(k):
        return min(thK_dt[m] * G[k - 4 - m] for m in range(0, min(5, k - 4) + 1))

    def B_V(k):
        c = [thK[m] * G[k - 4 - m] for m in range(0, min(7, k - 4) + 1)]
        if k - 4 <= 7:
            c.append(oscV[k - 4])
        return min(c)

    def B_W(k):
        return min(B_V(k), oscW[k - 4]) if k - 4 in oscW else B_V(k)

    prev = {k: (min(B_DO(k), B_U(k)) if k <= 24 else B_U(k)) for k in range(4, KM - 10)}
    okdu = all(B_U(k) < B_DO(k) for k in range(11, 25))
    gainV = {k: B_V(k) / B_W(k) for k in range(4, KM - 10)}
    fac = {k: prev[k] / B_W(k) for k in range(4, KM - 10)}
    kk = [k for k in range(4, KM - 10) if gainV[k] > 1]
    used = [k for k in range(4, KM - 10) if (k - 4) in oscW and B_W(k) == oscW[k - 4]]
    floor = {k: tab[k - 4]["up"] / tab[k - 4]["low"] for k in used if k - 4 in tab}
    check(used == [m + 4 for m in range(M0, MN + 1)] and kk == used,
          "DW.2: for k = " + ", ".join(map(str, used)) + " the direct term osc log R_(k-4) <= "
          + ", ".join(upW[k - 4].str(5, radius=False) + ("" if k - 4 in tab else " (stated)") for k in used) + " is the minimum of B_W and is strictly below "
          "DV's B_V (by the factors " + ", ".join(f"{f(gainV[k]).str(5, radius=False)} (k={k})" for k in used)
          + "); B_W = B_V at every other k; by DO.4 it is within " + ", ".join(f"{floor[k].str(6, radius=False)} (k={k})" for k in floor)
          + " of the floor of every CF.3-form bound")
    RESULT.update(used=used, gainV={k: f(gainV[k]) for k in used}, floor=floor)
    J = None
    for j0 in range(8, KM - 16):
        if all(G[j + 8] == th[8] * G[j] for j in range(j0, j0 + 8)):
            J = j0
            break
    tail_ok = J is not None and J + 19 < KM - 10 and J + 11 > MN + 4   # periodic window beyond the new direct rows
    strict = all(fac[k] > 1 for k in fac)
    per = [fac[k] for k in range(J + 11, J + 19)] if J is not None else [Fr(0)]
    show = (4, 10, 11, 12, 13, 14, 15, 16, 17, 22, 49, 103)
    check(okdu and strict and tail_ok and min(per) > 1,
          "DW.3: B_W(k) < min(B_DO(k), B_U(k)) for EVERY k >= 4 (exact rationals; directly k < " + str(KM - 10)
          + f", and G_U(j + 8) = theta_8 G_U(j) exactly for j >= {J}, so for k >= {(J or 0) + 11} the factor is 8-periodic, in "
          f"[{f(min(per)).str(5, radius=False)}, {f(max(per)).str(5, radius=False)}]); min(B_DO, B_U)/B_W: "
          + ", ".join(f"{f(fac[k]).str(5, radius=False)} (k={k})" for k in show)
          + f"; smallest {f(min(fac.values())).str(5, radius=False)} (k = {min(fac, key=fac.get)}). DV.3' (the strict "
          "form, false for B_V at k = 13, 16) holds for B_W")
    RESULT.update(J=J, fac={k: f(fac[k]) for k in show}, facmin=f(min(fac.values())), permin=f(min(per)),
                  BW={k: f(B_W(k)) for k in range(4, 21)})
    # consistency with CC's level sums (k = 4..KMAX-1; the new rows enter at k = 12)
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    rl, rh = arb(r_lo.numerator) / r_lo.denominator, arb(r_hi.numerator) / r_hi.denominator
    found, sums, counts, okb = cc.period_atoms(bk.KMAX, verbose=False)
    ctx.prec = PREC
    okd, worst, wk = okb, arb(0), None
    for k in range(4, bk.KMAX):
        dev = arb(max(abs((sums[k + 1] / (rl * sums[k])).log()).upper(), abs((sums[k + 1] / (rh * sums[k])).log()).upper()))
        b = f(B_W(k))
        okd = okd and bool(dev < b)
        if bool((dev / b).upper() > worst.upper()):
            worst, wk = arb((dev / b).upper()), k
    check(okd, f"DW.3 consistency: on CC's s_4..s_{bk.KMAX} every |log(s_(k+1)/(r s_k))| is below B_W(k) "
               f"(max ratio {worst.str(4, radius=False)} at k = {wk}; k = 12 uses the new row m = 8)")
    RESULT.update(worst=worst, wk=wk)


def main() -> int:
    t0 = time.time()
    print(f"Theorem DW: direct rows m = {M0}..{M1} (k = {M0 + 4}..{M1 + 4}), fourth-order centred forms ({NC0} initial "
          f"cells, depth <= {DEPTH}, tol {TOL}{', FULL' if FULL else ''}). Not thrust.", flush=True)
    T = di.taylor()
    if not tables(T):
        print("FAILED")
        return 1
    consequences()
    RESULT["elapsed"] = time.time() - t0
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(dv.FAILS) + len(dt.FAILS) + len(di.FAILS) + len(bk.FAILS) + len(bk.co.FAILS) + len(bk.cn.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DW checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
