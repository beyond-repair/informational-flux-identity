#!/usr/bin/env python3
r"""Theorem DV: the level-sum bound for k = 4..10 (DT's kept failure DT.4 "for k = 4..10 DO.3 stays the bound") and
direct constants for m = 6, 7 (DT's and DU's open item "direct constants for m >= 6").

Setting (GASKET.md, Theorems DO, DT, DU): s in X = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
K = K_(2,6) = {f > 0 : |f'| <= 2f, |f''| <= 6f} (DQ.1: Lt K in K), cut out by
ell in {delta_s, 2 delta_s +- delta_s', 6 delta_s +- delta_s''}, Theta_K(u, v) = log(sup_ell rho_ell / inf_ell rho_ell),
rho_ell = ell(v)/ell(u).  f_0 = F_0 = Lt^4 G (DC.2; CO.1/CO.2's Taylor data re-run unchanged), F_m = Lt^m F_0,
R_m = F_(m+1)/F_m, s_(4+m) = c0 F_m(x0), so s_(k+1)/s_k = R_(k-4)(x0) (DI.5's consistency).

Theorem DV (checked here; exit 1 on any failure):
  DV.1 [certified; arb jets of order 4 through m + 1 <= 8 layers of Lt, DT.2's second-order centred forms on adaptively
       bisected cells (a cell is split while its enclosure is far from the floors; the accepted cells cover X outward);
       floors: arb point values at DV_FLOOR_PTS points, endpoints included] for m = 0..7: F_m in int K_(2,6) and two-sided
       bounds for Theta_K(F_m, F_(m+1)) and for osc log R_m = log(sup R_m / inf R_m) (the delta-type alone).
  DV.2 [proved from DV.1 and the inclusion r in [inf R_m, sup R_m] of CF.3/DO.3's proof] cone-free direct bound:
       |log(s_(k+1)/(r s_k))| <= osc log R_(k-4) for k = 4..11.  It beats DO.3 at every k = 4..10 (k = 4..9 by the
       factor DO.4 allows, <= 1.0252; k = 10 because R_6 is now computed directly).  By DO.4 no CF.3-form bound (any cone)
       at k = 4 + m can be below osc log R_m, so for k = 4..11 this is the best bound of that form up to the enclosure width.
  DV.3 [proved from DV.1 with CF.1, CF.3's proof (DT.3), DQ.1 and DU.3's factors as stated]
       B_V(k) = min( osc log R_(k-4) [k <= 11],  min_(m <= min(7, k-4)) Theta_K(F_m, F_(m+1)) G_U(k-4-m) )
       satisfies B_V(k) <= min(B_DO(k), B_U(k)) for every k >= 4, with equality exactly at the k listed by the script
       (k = 13, 16) and strict inequality elsewhere (exact rational comparison with outward-rounded inputs; for large k by the
       exact 8-periodicity of the knapsack G_U, checked on a window).  The strict form "for every k >= 4" is false at
       k = 13, 16 and is kept as a rejected claim (DV.3').
Kept limits: see GASKET.md Theorem DV.  No Psi_info, Stage 2 still blocked.  Stage 1 untouched; no W, 0.08, 0.23 or
30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction as Fr

from math import comb

from flint import arb, arb_poly, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_second_order_sums as dt   # noqa: E402  (DT: order-4 jets, ell, centred form)

di, do = dt.di, dt.do
bk = di.bk
FAILS: list[str] = []
RESULT: dict = {}
PREC = int(os.environ.get("DV_PREC", "128"))      # cell jobs; CO's Taylor data stay at 200 bits
FULL = os.environ.get("DV_FULL", "") not in ("", "0")
MV = 7                                                   # direct rows m = 0..MV
NC0 = int(os.environ.get("DV_CELLS0", "64"))             # initial uniform cells
DEPTH = int(os.environ.get("DV_DEPTH", "10" if FULL else "8"))
TOL = float(os.environ.get("DV_TOL", "0.0005" if FULL else "0.01"))   # relative target slack (log scale) per side
NFP = int(os.environ.get("DV_FLOOR_PTS", "257" if FULL else "33"))
WORKERS = int(os.environ.get("DV_WORKERS", "4"))
ELLS, ELLN = dt.ELLS, dt.ELLN
# DU.3's factors, used as stated: Dbar_n (n = 3..6 from DR.3/DS.3, n = 7, 8 from DU.2)
DBAR_U = {3: "0.25513308", 4: "0.03838589", 5: "0.0098239159", 6: "0.0018377053", 7: "0.00023483", 8: "0.000041213"}
DT_UPPER = dt.TARGETS[1024]["UPPER"]                     # DT.2's uppers for Theta_K(F_m, F_(m+1)), m <= 5, as stated
DO_UPPER = dt.DO_UPPER
DO_PAIRS = dt.DO_PAIRS
# certified targets (checked, exit 1 if not met): uppers for Theta_K and osc log R, m = 0..7; lowers likewise
TARGET_K_UP = ('0.14742', '0.072562', '0.0069590', '0.0010079', '0.00016637', '2.8382e-5', '4.8865e-6', '8.4251e-7')
TARGET_OSC_UP = ('0.063415', '0.036883', '0.0041590', '0.00060873', '0.00010052', '1.7204e-5', '2.9652e-6', '5.1065e-7')
TARGET_K_LO = ('0.14731', '0.072340', '0.0069152', '0.0010008', '0.00016484', '2.8125e-5', '4.8402e-6', '8.3466e-7')
TARGET_OSC_LO = ('0.063386', '0.036827', '0.0041515', '0.00060778', '0.00010040', '1.7143e-5', '2.9507e-6', '5.0884e-7')


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


class TaylorJet4P(dt.TaylorJet4):
    """DT's TaylorJet4 (CO.2's series with DO.1's Cauchy tail, order 4), the Horner sums evaluated by arb_poly (C)."""

    def __init__(self, T):
        super().__init__(T)
        self.P = [arb_poly(self.d[j]) for j in range(dt.KJ)]

    def jet(self, X):
        KJ = dt.KJ
        U0 = (X[0] - self.c) / self.R
        q = arb(abs(U0).upper()) / self.rho
        assert bool(q < 1)
        M = self.N + 1
        c = []
        for j in range(KJ):
            tail = self.S * self.rho ** (-j) * q ** (M - j) * sum(
                (comb(M, j - l) * q ** l / (1 - q) ** (l + 1) for l in range(j + 1)), arb(0))
            c.append(self.P[j](U0) + arb(0, tail.upper()))
        dU = [arb(0)] + [x / self.R for x in X[1:]]
        out, p = [c[0]] + [arb(0)] * (KJ - 1), [arb(1)] + [arb(0)] * (KJ - 1)
        for j in range(1, KJ):
            p = dt.jmul(p, dU)
            out = dt.jadd(out, [c[j] * x for x in p])
        return out


def _jet(TJ, m, x):
    one, zero = arb(1), arb(0)
    return dt.Fjet(TJ, m, [x] + [one] + [zero] * (dt.KJ - 2))


# ----------------------------------------------------------------- DV.1
def floors(TJ, m):
    """arb point values of rho_ell at NFP points (endpoints included): per type, certified (max lower, min upper)."""
    ctx.prec = PREC
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    hi, lo = [None] * 5, [None] * 5
    for k in range(NFP):
        s = s0 if k == 0 else (s1 if k == NFP - 1 else arb((s0 + (s1 - s0) * k / (NFP - 1)).mid()))
        u, v = _jet(TJ, m, s), _jet(TJ, m + 1, s)
        for e, a in enumerate(ELLS):
            r = dt.ell(a, v)[0] / dt.ell(a, u)[0]
            hi[e] = r.lower() if hi[e] is None else max(hi[e], r.lower())
            lo[e] = r.upper() if lo[e] is None else min(lo[e], r.upper())
    return hi, lo


def cell_eval(TJ, m, a, b):
    """Per type (lo, hi) of rho_ell on the cell [a, b] (DT.2's centred form intersected with the direct quotient), and
    the cell's |F_m'|/F_m upper, F_m''/F_m range, min ell(F_m)/F_m."""
    J = bk._hull(a, b)
    mp = arb(J.mid())
    dx = J - mp
    uJ, vJ, uM, vM = _jet(TJ, m, J), _jet(TJ, m + 1, J), _jet(TJ, m, mp), _jet(TJ, m + 1, mp)
    if not bool(uJ[0] > 0):
        return None
    q1 = abs(uJ[1] / uJ[0]).upper()
    q2 = 2 * uJ[2] / uJ[0]
    out, emin = [], None
    for a_ in ELLS:
        AJ, BJ, AM, BM = dt.ell(a_, vJ), dt.ell(a_, uJ), dt.ell(a_, vM), dt.ell(a_, uM)
        if not (bool(BJ[0] > 0) and bool(AJ[0] > 0)):
            return None
        e = (BJ[0] / uJ[0]).lower()
        emin = e if emin is None else min(emin, e)
        rJ, rM = dt.jdiv3(AJ, BJ), dt.jdiv3(AM, BM)
        rc = rM[0] + rM[1] * dx + rJ[2] * dx * dx
        rd = rJ[0]
        out.append((max(rc.lower(), rd.lower()), min(rc.upper(), rd.upper())))
    return out, q1, q2.lower(), q2.upper(), emin


DM = 30                                                  # cell endpoints x(p) = s0 + (s1 - s0) p / 2^DM (p integer)
_TJ = None


def _x(p):
    s0l, s1u = arb((arb(5).sqrt() / 2).lower()), arb((arb(5) / 2).upper())
    return s0l + (s1u - s0l) * p / 2 ** DM


def _fr(x):
    """exact Fraction of an exact arb (an arb .lower()/.upper())."""
    mn, e = x.man_exp()
    return Fr(int(mn)) * Fr(2) ** int(e)


def _cell_job(args):
    m, p, q = args
    ctx.prec = PREC
    r = cell_eval(_TJ, m, _x(p), _x(q))
    if r is None:
        return (p, q, None)
    out, q1, l2, h2, e = r
    return (p, q, ([(_fr(a), _fr(b)) for a, b in out], _fr(q1), _fr(l2), _fr(h2), _fr(e)))


def _pool():
    import multiprocessing as mp
    return mp.get_context("fork").Pool(WORKERS) if WORKERS > 1 else None


def _arb(x: Fr):
    return arb(x.numerator) / arb(x.denominator)


def row(TJ, m, pool):
    """Two-sided Theta_K(F_m, F_(m+1)) and osc log R_m on adaptively bisected cells (waves; all cells of a wave in
    parallel).  A cell is accepted when its enclosures are within TOL * floor (log scale) of the floor extremes, or at
    depth DEPTH; the accepted cells tile [0, 2^DM] in index space, so their hulls cover X."""
    ctx.prec = PREC
    t0 = time.time()
    fhi, flo = floors(TJ, m)
    bK_f, aK_f = max(fhi), min(flo)
    bD_f, aD_f = fhi[0], flo[0]
    lowK = arb((arb(bK_f) / arb(aK_f)).log().lower())
    lowD = arb((arb(bD_f) / arb(aD_f)).log().lower())
    epsK = Fr(float(lowK.lower()) * TOL / 2)
    epsD = Fr(float(lowD.lower()) * TOL / 2)
    capKb, capKa = _fr(bK_f) * (1 + epsK), _fr(aK_f) * (1 - epsK)
    capDb, capDa = _fr(bD_f) * (1 + epsD), _fr(aD_f) * (1 - epsD)
    step = 2 ** DM // NC0
    wave = [(m, i * step, (i + 1) * step) for i in range(NC0)]
    acc = []
    nsplit, d = 0, 0
    while wave:
        res = pool.map(_cell_job, wave, chunksize=max(1, len(wave) // (8 * WORKERS))) if pool else list(map(_cell_job, wave))
        nxt = []
        for p, q, r in res:
            ok = False
            if r is not None:
                out = r[0]
                hiK, loK = max(x[1] for x in out), min(x[0] for x in out)
                ok = hiK <= capKb and loK >= capKa and out[0][1] <= capDb and out[0][0] >= capDa
            if ok or (r is not None and d >= DEPTH):
                acc.append((p, q, r))
            elif d < DEPTH + 4 and q - p >= 2:
                c = (p + q) // 2
                nxt += [(m, p, c), (m, c, q)]
                nsplit += 1
            else:
                check(False, f"DV.1: positivity of F_m or ell(F_m) not certified on cell ({p}, {q}) (m = {m})")
                return None
        wave, d = nxt, d + 1
    acc.sort()
    tiles = acc[0][0] == 0 and acc[-1][1] == 2 ** DM and all(acc[i][1] == acc[i + 1][0] for i in range(len(acc) - 1))
    if not tiles:
        check(False, f"DV.1: accepted cells do not tile X (m = {m})")
        return None
    bK = max(max(x[1] for x in r[0]) for _, _, r in acc)
    aK = min(min(x[0] for x in r[0]) for _, _, r in acc)
    bD = max(r[0][0][1] for _, _, r in acc)
    aD = min(r[0][0][0] for _, _, r in acc)
    upK = arb((_arb(bK) / _arb(aK)).log().upper())
    upD = arb((_arb(bD) / _arb(aD)).log().upper())
    q1max = max(r[1] for _, _, r in acc)
    q2lo = min(r[2] for _, _, r in acc)
    q2hi = max(r[3] for _, _, r in acc)
    emin = min(r[4] for _, _, r in acc)
    maxd = max(step.bit_length() - (q - p).bit_length() for p, q, _ in acc)
    print(f"   m = {m}: Theta_K in [{lowK.str(8, radius=False)}, {upK.str(8, radius=False)}], osc log R_m in "
          f"[{lowD.str(8, radius=False)}, {upD.str(8, radius=False)}]; {len(acc)} cells (max depth {maxd}, {nsplit} splits), "
          f"{NFP} floor points ({time.time() - t0:.0f}s)", flush=True)
    return dict(upK=upK, lowK=lowK, upD=upD, lowD=lowD, q1max=_arb(q1max), q2lo=_arb(q2lo), q2hi=_arb(q2hi),
                emin=_arb(emin), ncell=len(acc), maxd=maxd)


def main_tables(T):
    global _TJ
    TJ = _TJ = TaylorJet4P(T)
    tab = {}
    pool = _pool()
    try:
        for m in range(MV + 1):
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


# ----------------------------------------------------------------- DV.2 / DV.3 (exact rationals)
KM = 400                                                 # k checked directly; beyond, the 8-periodicity of G_U


def _up(x) -> Fr:
    """exact Fraction >= x (x an arb or a decimal string)."""
    ctx.prec = PREC
    a = x if isinstance(x, arb) else arb(x)
    return _fr(arb(a.upper()))


def knapF(pairs, jmax):
    """G(j), j = 0..jmax, exactly: least product of the factors whose block lengths sum to at most j."""
    G = [Fr(1)] * (jmax + 1)
    for j in range(1, jmax + 1):
        G[j] = min([G[j - 1]] + [G[j - n] * t for n, t in pairs if j >= n])
    return G


def consequences():
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    import gasket_fractional_profile_pure_point as cc
    tab = RESULT["tab"]
    # DV.1 checks: membership, two-sided rows, monotonicity Theta_K >= osc, targets
    q1 = max(tab[m]["q1max"].upper() for m in tab)
    q2l = min(tab[m]["q2lo"].lower() for m in tab)
    q2h = max(tab[m]["q2hi"].upper() for m in tab)
    em = min(tab[m]["emin"].lower() for m in tab)
    check(bool(arb(q1) < 2) and bool(arb(q2l) > -6) and bool(arb(q2h) < 6) and bool(arb(em) > 0),
          f"DV.1: |F_m'|/F_m <= {arb(q1).str(6, radius=False)} < 2, F_m''/F_m in [{arb(q2l).str(6, radius=False)}, "
          f"{arb(q2h).str(6, radius=False)}] inside (-6, 6), min_ell ell(F_m)/F_m >= {arb(em).str(6, radius=False)} > 0 for "
          f"m = 0..{MV}: every F_m in int K_(2,6) (consistent with DT.1 for m <= 5 and with DQ.1)")
    ok = all(bool(tab[m]["lowK"] <= tab[m]["upK"]) and bool(tab[m]["lowD"] <= tab[m]["upD"])
             and bool(tab[m]["upD"] <= tab[m]["upK"]) for m in tab)
    tg = all(bool(tab[m]["upK"] < arb(TARGET_K_UP[m])) and bool(tab[m]["upD"] < arb(TARGET_OSC_UP[m]))
             and bool(tab[m]["lowK"] > arb(TARGET_K_LO[m])) and bool(tab[m]["lowD"] > arb(TARGET_OSC_LO[m]))
             for m in tab)
    check(ok and tg, "DV.1: two-sided rows (arb): " + "; ".join(
        f"m={m}: Theta_K in [{tab[m]['lowK'].str(6, radius=False)}, {tab[m]['upK'].str(6, radius=False)}], osc log R_m in "
        f"[{tab[m]['lowD'].str(6, radius=False)}, {tab[m]['upD'].str(6, radius=False)}]" for m in tab)
        + " (osc <= Theta_K since delta is in Lambda; targets met)")
    # inputs as exact rationals (outward)
    osc = {m: _up(tab[m]["upD"]) for m in tab}
    thK = {m: _up(tab[m]["upK"]) for m in tab}
    for m in range(6):
        thK[m] = min(thK[m], _up(DT_UPPER[m]))                       # DT.2's uppers, as stated
    th = {n: _up((arb(D) / 4).tanh()) for n, D in DBAR_U.items()}
    pairsU = sorted(th.items())
    G = knapF(pairsU, KM)
    GD = {a: knapF([(n, _up(t)) for n, t in DO_PAIRS[a]], 20) for a in (2, 100)}

    def B_DO(k):
        return min(_up(DO_UPPER[a][m]) * GD[a][k - 4 - m] for m in range(0, min(5, k - 4) + 1) for a in (2, 100))

    def B_U(k):
        return min(thK_dt[m] * G[k - 4 - m] for m in range(0, min(5, k - 4) + 1))

    def B_V(k):
        c = [thK[m] * G[k - 4 - m] for m in range(0, min(MV, k - 4) + 1)]
        if k - 4 <= MV:
            c.append(osc[k - 4])
        return min(c)

    thK_dt = {m: _up(DT_UPPER[m]) for m in range(6)}
    # the previous bound: min(B_DO, B_U) wherever B_DO's knapsack table reaches (k <= 24), B_U beyond (DU.3: B_U < B_DO
    # for k >= 11, as stated; recomputed here for k = 11..24 as a consistency check)
    prev = {k: (min(B_DO(k), B_U(k)) if k <= 24 else B_U(k)) for k in range(4, KM - 10)}
    okdu = all(B_U(k) < B_DO(k) for k in range(11, 25))
    fac = {k: prev[k] / B_V(k) for k in range(4, KM - 10)}
    # tail: exact 8-periodicity of G (G(j+8) = theta_8 G(j) on a window of 8 consecutive j >= 8 implies it for all later j)
    J = None
    for j0 in range(8, KM - 16):
        if all(G[j + 8] == th[8] * G[j] for j in range(j0, j0 + 8)):
            J = j0
            break
    tail_ok = J is not None and J + 19 < KM - 10
    weak = all(fac[k] >= 1 for k in fac)                      # B_V <= min(B_DO, B_U), checked k < KM - 10
    eq = [k for k in sorted(fac) if fac[k] == 1]               # exact rational equality
    per_eq = [k for k in range(J + 11, J + 19) if fac[k] == 1] if J is not None else []
    show = (4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 22, 49, 103)
    per = [fac[k] for k in range(J + 11, J + 19)] if J is not None else [Fr(0)]

    def f(x, d=4):
        return arb(x.numerator) / arb(x.denominator) if isinstance(x, Fr) else x

    check(okdu and weak and tail_ok and not per_eq,
          "DV.3: B_V(k) = min(osc log R_(k-4) [k <= " + str(4 + MV) + "], min_(m <= min(" + str(MV) + ", k-4)) Theta_K(F_m, "
          "F_(m+1)) G_U(k-4-m)) <= min(B_DO(k), B_U(k)) for every k >= 4, with equality exactly at k in " + str(eq)
          + " and strict inequality at every other k (exact rationals; DO.3's and DU.3's constants as stated, B_U < B_DO "
          "for k = 11..24 recomputed). min(B_DO, B_U)/B_V: "
          + ", ".join(f"{f(fac[k]).str(5, radius=False)} (k={k})" for k in show)
          + f"; G_U(j + 8) = theta_8 G_U(j) exactly for every j >= {J}, so for k >= {(J or 0) + 11} the factor is 8-periodic, "
          f"in [{f(min(per)).str(5, radius=False)}, {f(max(per)).str(5, radius=False)}] (no equality there); B_V(k): "
          + ", ".join(f"{f(B_V(k)).str(5, radius=False)} (k={k})" for k in show))
    # the earlier, mislabeled strict claim, kept as a rejected claim (it is false exactly at the equality k's)
    check(len(eq) > 0 and all(fac[k] > 1 for k in fac if k not in eq),
          "DV.3' (rejected strict form, kept): 'B_V(k) < min(B_DO(k), B_U(k)) for every k >= 4' is FALSE: "
          f"min(B_DO, B_U)/B_V = 1 exactly at k = {eq} (there no m = 6, 7 term and no direct term improves the DU.3 "
          "minimum); strict at every other k")
    RESULT.update(eq=eq, fac_all={k: fac[k] for k in range(4, 41)})
    RESULT.update(fac={k: f(fac[k]) for k in show}, BV={k: f(B_V(k)) for k in show}, J=J,
                  permin=f(min(per)), permax=f(max(per)))
    # DV.2: k = 4..4+MV, the cone-free direct term; DO.4's floor makes it optimal among CF.3-form bounds (up to enclosures)
    opt = {k: tab[k - 4]["upD"] / tab[k - 4]["lowD"] for k in range(4, 5 + MV)}
    used = [k for k in range(4, 5 + MV) if B_V(k) == osc[k - 4]]
    check(used == list(range(4, 5 + MV)) and all(fac[k] > 1 for k in range(4, 11)),
          "DV.2: for k = 4.." + str(4 + MV) + " the minimum in B_V(k) is the cone-free direct term osc log R_(k-4) <= "
          + ", ".join(f"{tab[k - 4]['upD'].str(6, radius=False)}" for k in range(4, 5 + MV))
          + "; it beats DO.3 at every k = 4..10 (DT.4's range): factors "
          + ", ".join(f"{f(fac[k]).str(5, radius=False)}" for k in range(4, 11))
          + "; by DO.4 every CF.3-form bound at k is >= osc log R_(k-4), and the upper is within the factors "
          + ", ".join(f"{opt[k].str(6, radius=False)}" for k in opt) + " of that floor")
    RESULT.update(used=used, opt=opt)
    # consistency with CC's level sums (r at either end of CE.4's interval)
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    rl = arb(r_lo.numerator) / r_lo.denominator
    rh = arb(r_hi.numerator) / r_hi.denominator
    found, sums, counts, okb = cc.period_atoms(bk.KMAX, verbose=False)
    ctx.prec = PREC
    okd, worst, wk = okb, arb(0), None
    for k in range(4, bk.KMAX):
        dev = arb(max(abs((sums[k + 1] / (rl * sums[k])).log()).upper(), abs((sums[k + 1] / (rh * sums[k])).log()).upper()))
        b = f(B_V(k))
        okd = okd and bool(dev < b)
        if bool((dev / b).upper() > worst.upper()):
            worst, wk = arb((dev / b).upper()), k
    check(okd, f"DV.3 consistency: on CC's s_4..s_{bk.KMAX} every deviation |log(s_(k+1)/(r s_k))| is below B_V(k) "
               f"(max ratio {worst.str(4, radius=False)} at k = {wk})")
    RESULT.update(worst=worst, wk=wk)


def main() -> int:
    t0 = time.time()
    print(f"Theorem DV: direct constants m = 0..{MV} and the level-sum bound for k = 4..10 ({NC0} initial cells, depth "
          f"<= {DEPTH}, tol {TOL}{', FULL' if FULL else ''}). Not thrust.", flush=True)
    T = di.taylor()
    if not main_tables(T):
        print("FAILED")
        return 1
    consequences()
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(dt.FAILS) + len(di.FAILS) + len(bk.FAILS) + len(bk.co.FAILS) + len(bk.cn.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DV checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
