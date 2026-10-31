#!/usr/bin/env python3
r"""Theorem DH: DG's kept limits "(a, n) = (100, 8) not re-optimised with DG.1" and "Theta_G' not redone with DG.1".
Same cone C_100, block lengthened from Lt^8 to Lt^9 (512 words); DG.1's factor lemma made asymmetric for the f_0 term.

Setting (GASKET.md, Theorems CE, CF, DC-DG): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
C_a = {f > 0 : log f a-Lipschitz in s}, Theta_a its Hilbert metric (CF's definition), a = 100.

Theorem DH (checked here; exit 1 on any failure):
  DH.1 [proved] (asymmetric factor lemma) If x > y => -beta_u (x-y) <= log u(x)/u(y) <= alpha_u (x-y) and the same for v
       with (alpha_v, beta_v), all four in [0, a), then
           Theta_a(u, v) <= log(sup(v/u)/inf(v/u)) + log K_up + log K_dn,
           K_up = max((a+beta_v)/(a-alpha_u), (a+alpha_v)/(a-beta_u)),  K_dn = max((a+beta_u)/(a-alpha_v), (a+alpha_u)/(a-beta_v)).
       (DG.1(ii) is the case of one common pair of slope bounds; DG.1's signed cases are used unchanged for Lt^9 C_100.)
  DH.2 [certified; DG.2's arb cell tables and outward DP with N = 9, 4096 cells] signed slopes of log Lt^9 f, f in C_100.
  DH.3 [certified; DF.2's adaptive outward joint programme with N = 9 and the Lipschitz pad 2 D+ H (Lt^9 C_100 in
       C_(D+) by DH.2)] W_H >= log(sup(v/u)/inf(v/u)) on Lt^9 C_100; floor from explicit bang-bang f, g (arb).
  DH.4 [certified; arb] Delta_H = 2 log K + W_H, theta_H = tanh(Delta_H/4), per-level rate theta_H^(1/9).
  DH.5 [proved from DH.1-DH.4, CF.1, DF.4's proofs with (a, n) = (100, 9), DC.2's b_0, B_0, log(sup R_0/inf R_0) as
       stated] |mu| <= theta_H^(1/9) r; Theta_100(f_0, Lt f_0) <= Theta_H'; level sums.
Kept limits: see GASKET.md Theorem DH.  No Psi_info, Stage 2 still blocked.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_hilbert as dg           # noqa: E402

df = dg.df
de = dg.de
dd = dg.dd
cf = dg.cf
FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A = 100.0
A_NUM = 100
N = 9
CELLS = int(os.environ.get("DH_CELLS", "4096"))
W_TARGET = os.environ.get("DH_TARGET", "0.0142")
b0_DC = "0.77333706"            # DC.2: sup |(log f_0)'| <= b_0 <= 0.77333706, used as stated
B0_DC = df.B0                   # DC.2: Lt f_0 in C_(phi(b_0)), phi(b_0) = B_0 <= 1.2791803, used as stated
LOGR0 = df.LOGR0                # DC.2: log(sup R_0/inf R_0) <= 0.063592701, used as stated
THETA_G, RATE_G, CONE_G, TG_G, MU_G = "0.014183", "0.58745", "0.012225", "0.11477", "0.24893"
DP_MAX = "0.74455"
DM_MAX = "-0.16327"
W_MAX = "0.0142"
WLOW_MIN = "0.013764"
CT_MAX = "0.011679"
DELTA_MAX = "0.025879"
THETA_MAX = "0.0064696"
RATE_MAX = "0.57117"
MU_MAX = "0.24203"
TH_MAX = "0.10465"
DPLOW_MIN = "0.74356"
DMLOW_MIN = "-0.16376"
INF, NINF = math.inf, -math.inf


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _set_block():
    """DG's and DF's routines read the block length and cone from module globals; set them to (100, 9)."""
    for mod in (dg, df):
        mod.N = N
    dg.A, dg.A_NUM = A, A_NUM
    df.A, df.A_ARB_NUM = A, A_NUM
    df.W_TARGET = W_TARGET
    df.CELLS, df.LEVELS, df.FAC = 128, int(os.environ.get("DH_LEVELS", "6")), 2
    df.CH = int(os.environ.get("DH_CHUNK", "512"))


# ----------------------------------------------------------------- DH.2 signed slopes for Lt^9
def slopes(m=None):
    m = m or CELLS
    t0 = time.time()
    Dp = Dm = NINF
    argp = argm = None
    for j in range(m):
        WL, WH, NP, NM, IH, IL = dg.cell_coeffs(m, j)
        p, q = dg.sup_ratio(NP, WL, WH, IH, IL), dg.sup_ratio(NM, WL, WH, IH, IL)
        if p > Dp:
            Dp, argp = p, j
        if q > Dm:
            Dm, argm = q, j
    print(f"   DH.2: {m} cells, D+ <= {Dp:.8f} (cell {argp}), D- <= {Dm:.8f} (cell {argm}) ({time.time() - t0:.1f}s)",
          flush=True)
    return Dp, Dm, argp, argm, m


def slope_floor(cell, m, sign):
    """Certified lower bound for sup_(f in C_100, x > y) sign (log Lt^9 f(x) - log Lt^9 f(y))/(x - y) from DF.3's
    bang-bang construction at two points of the cell (x - y = 1e-4 or the cell width)."""
    s0, L = math.sqrt(5) / 2, (5 - math.sqrt(5)) / 2
    c = s0 + (cell + 0.5) * L / m
    best = arb(-10)
    ctx.prec = 120
    a = arb(A_NUM)
    for h in (1e-3, 1e-4, 1e-5):
        x, y = min(c + h / 2, 2.5), max(c - h / 2, s0)
        xa, ya = arb(x), arb(y)
        if sign > 0:          # maximise w(x)/w(y)
            p, sg = df._bangbang(x, y, A)
            val = df._eval_ratio_arb(xa, ya, p, sg, a).log() / (xa - ya)
        else:                 # maximise w(y)/w(x)
            p, sg = df._bangbang(y, x, A)
            val = df._eval_ratio_arb(ya, xa, p, sg, a).log() / (xa - ya)
        if bool(val.lower() > best.lower()):
            best = arb(val.lower())
    return best


def w_floor():
    ctx.prec = 120
    s0, s1 = math.sqrt(5) / 2, 2.5
    a = arb(A_NUM)
    best = arb(0)
    for x, y in [(s0 + 1e-9, s1 - 1e-9), (s0 + 1e-9, 2.4), (s0 + 0.01, s1 - 1e-9)]:
        xa, ya = arb(x), arb(y)
        pf, sf = df._bangbang(y, x, A)
        pg, sg = df._bangbang(x, y, A)
        val = df._eval_ratio_arb(ya, xa, pf, sf, a).log() + df._eval_ratio_arb(xa, ya, pg, sg, a).log()
        if bool(val.lower() > best.lower()):
            best = arb(val.lower())
    return best


def main_bound():
    _set_block()
    df.check_separated()
    Dp_f, Dm_f, argp, argm, m = slopes()
    ctx.prec = PREC
    Dp, Dm = arb(Dp_f), arb(Dm_f)
    a = arb(A_NUM)
    check(bool(Dp < arb(DP_MAX)) and bool(Dm < arb(DM_MAX)),
          f"DH.2: D+ <= {Dp_f:.8f} < {DP_MAX}, D- <= {Dm_f:.8f} < {DM_MAX}: for f in C_100 and x > y, "
          f"{-Dm_f:.5f} (x-y) <= log Lt^9 f(x) - log Lt^9 f(y) <= {Dp_f:.5f} (x-y)")
    bn = dd.b_iterate(a, N)
    check(bool(Dp < bn) and bool(-Dm < Dp) and bool(Dm < 0),
          f"DH.2: Lt^9 C_100 subset C_({DP_MAX}) (CF.2's b_9 = {bn.str(8)}); every log Lt^9 f is increasing")
    hyp = (bool(abs(Dp) < a) and bool(abs(Dm) < a) and bool(Dp + Dm > 0) and bool(Dp <= a / 2))
    check(hyp, "DG.1 for Lt^9: hypotheses of the factor lemma hold (|D+-| < a, D+ + D- > 0, D+ <= a/2)")
    K = arb(max((((a + Dm) / (a - Dp))).upper(), ((a + Dp) / (a - Dm)).upper()))
    cone = arb((2 * K.log()).upper())
    check(bool(cone < arb(CT_MAX)), f"DH.4: cone term 2 log K <= {cone.str(8)} < {CT_MAX} (DG, Lt^8: {CONE_G})")
    # floors for D+- (true sup slopes), DF.3's bang-bang construction at the maximising cells
    ctx.prec = PREC
    fp = slope_floor(argp, m, +1)
    fm = slope_floor(argm, m, -1)
    ctx.prec = PREC
    check(bool(fp > arb(DPLOW_MIN)) and bool(fp <= Dp) and bool(fm > arb(DMLOW_MIN)) and bool(fm <= Dm),
          f"DH.2 (floor): explicit f in C_100 give a slope of log Lt^9 f >= {fp.str(8)} > {DPLOW_MIN} and a slope <= "
          f"{(-fm).str(8)} (so D+_true >= {DPLOW_MIN}, D-_true >= {DMLOW_MIN}); D+/D+_low <= {(Dp / fp).str(5)}, "
          f"gap D- - D-_low <= {(Dm - fm).str(5)}")
    # DH.3: joint term with the Lipschitz pad D+
    t0 = time.time()
    Wf, Smax, stats = df.adaptive_W(float(Dp.upper()))
    Wf = float(dg._up(arb(Wf + 1e-12)))
    ctx.prec = PREC
    W = arb(Wf)
    print(f"   DH.3: {time.time() - t0:.1f}s")
    check(bool(W < arb(W_MAX)), f"DH.3: log(sup(v/u)/inf(v/u)) <= W_H = {W.str(8)} < {W_MAX} for all u, v in Lt^9 C_100 "
                                f"(DF.2, Lt^8: 0.04451)")
    Wl = w_floor()
    ctx.prec = PREC
    check(bool(Wl > arb(WLOW_MIN)) and bool(Wl < W),
          f"DH.3 (floor): explicit bang-bang f, g in C_100 give log(sup(v/u)/inf(v/u)) >= {Wl.str(8)} > {WLOW_MIN}; "
          f"W_H/W_low <= {(W / Wl).str(5)}")
    Delta = arb((cone + W).upper())
    theta = arb((Delta / 4).tanh().upper())
    rate = arb((theta.log() / N).exp().upper())
    check(bool(Delta < arb(DELTA_MAX)) and bool(theta < arb(THETA_MAX)) and bool(rate < arb(RATE_MAX)),
          f"DH.4: Delta_H = 2 log K + W_H <= {Delta.str(8)} < {DELTA_MAX}, theta_H <= {theta.str(8)} < {THETA_MAX}, "
          f"theta_H^(1/9) <= {rate.str(8)} < {RATE_MAX} (DG: 0.056735, {THETA_G}, {RATE_G})")
    check(bool(rate < arb(RATE_G)), f"DH.4: the per-level rate is strictly below DG's {RATE_G}")
    RESULT.update(Dp=Dp, Dm=Dm, K=K, cone=cone, W=W, Wlow=Wl, Dplow=fp, Dmlow=fm, Delta=Delta, theta=theta,
                  rate=rate, bn=bn, stats=stats)


# ----------------------------------------------------------------- DH.1 asymmetric factor lemma for the f_0 term
def theta_prime():
    ctx.prec = PREC
    a = arb(A_NUM)
    bu, bv = arb(b0_DC), arb(B0_DC)        # u = f_0 in C_(b_0), v = Lt f_0 in C_(B_0): symmetric slope bounds
    check(bool(bu < a) and bool(bv < a) and bool(bu <= bv), "DH.1: 0 <= b_0 <= B_0 < a (all four slope bounds in [0, a))")
    Kup = (a + bv) / (a - bu)
    Kdn = (a + bu) / (a - bv)
    TH = arb((Kup.log() + Kdn.log() + arb(LOGR0)).upper())
    TG = arb((2 * ((a + bv) / (a - bv)).log() + arb(LOGR0)).upper())
    check(bool(TH < arb(TH_MAX)) and bool(TH < TG),
          f"DH.5: Theta_100(f_0, Lt f_0) <= Theta_H' = log K_up + log K_dn + log(sup R_0/inf R_0) = {TH.str(8)} < {TH_MAX} "
          f"(DF.4/DG.4 cone lemma: Theta_G' = {TG.str(8)}; cone part {(TH - arb(LOGR0)).str(6)} vs "
          f"{(TG - arb(LOGR0)).str(6)})")
    # numerical sanity (not part of the proof): factor bounds on a d-grid in both orientations
    ok = True
    for al_u, be_v, bound in ((bu, bv, Kup), (bv, bu, Kdn)):
        for e in range(-12, 2):
            for mant in (1, 3):
                d = arb(mant) * arb(10) ** e
                g = ((a * d).exp() - (-be_v * d).exp()) / ((a * d).exp() - (al_u * d).exp())
                ok = ok and bool(g <= bound)
    check(ok, "DH.1 (sanity, arb): the asymmetric factor bounds hold on d = 1e-12..30")
    RESULT.update(TH=TH, TGold=TG)


# ----------------------------------------------------------------- DH.5 consequences
def consequences():
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    theta, rate, TH = RESULT["theta"], RESULT["rate"], RESULT["TH"]
    a = arb(A_NUM)
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    mu = rate * dd._arb(r_hi)
    check(bool(mu < arb(MU_MAX)), f"DH.5: |mu| <= theta_H^(1/9) r <= {mu.str(8)} < {MU_MAX} (DG.4: < {MU_G})")
    ok, Th0, w0 = df._consistency_f1(12, theta, r_lo, r_hi, a)
    check(ok, f"DH.5 (consistency): CF.4's f = 1 check with Theta_0(100) = {Th0.str(6)} and theta_H^floor(m/9) for m < 12 "
              f"(max used fraction {w0:.3e})")
    import gasket_fractional_profile_pure_point as cc
    found, sums, counts, okb = cc.period_atoms(dd.KMAX, verbose=False)
    ctx.prec = PREC
    r = (dd._arb(r_lo) + dd._arb(r_hi)) / 2
    thG, TGold = arb(THETA_G), RESULT["TGold"]
    okd, worst, worstG = okb, arb(0), arb(0)
    for k in range(4, dd.KMAX):
        dev = abs((sums[k + 1] / (r * sums[k])).log())
        bH = theta ** ((k - 4) // N) * TH
        bG8 = thG ** ((k - 4) // 8) * TH
        okd = okd and bool(dev < bH) and bool(dev < bG8)
        worst = arb(max(worst.upper(), (dev / bH).upper()))
        worstG = arb(max(worstG.upper(), (dev / bG8).upper()))
    check(okd, f"DH.5 (consistency): on CC's s_4..s_{dd.KMAX}, |log(s_(k+1)/(r s_k))| is below theta_H^floor((k-4)/9) "
               f"Theta_H' (max ratio {worst.str(4)}) and below theta_G^floor((k-4)/8) Theta_H' (max ratio {worstG.str(4)})")
    # comparison: both forms hold with Theta_H'; DH.5's level-sum bound is their minimum
    ratio = TGold / TH
    check(bool(ratio > 1), f"DH.5: with Theta_H', DG.4's 8-block form theta_G^floor((k-4)/8) Theta_H' is DG.4's bound divided by "
                           f"Theta_G'/Theta_H' = {ratio.str(6)} at every k >= 4")
    nine = []
    for k in range(4, 200):
        if bool(theta ** ((k - 4) // N) < thG ** ((k - 4) // 8)):
            nine.append(k)
    # tail: for k >= 200, theta_H^floor((k-4)/9) / theta_G^floor((k-4)/8) <= theta_H^((k-12)/9) / theta_G^((k-4)/8),
    # decreasing because theta_H^(1/9) < theta_G^(1/8)
    k0 = 200
    tail = (theta.log() * (k0 - 12) / N).exp() / (thG.log() * (k0 - 4) / 8).exp()
    check(bool(tail < 1) and bool(rate < (thG.log() / 8).exp()),
          f"DH.5: the 9-block form is strictly smaller than the 8-block form for {len(nine)} of the k in 4..199 (first: "
          f"{nine[:8]}) and for every k >= 200 (decreasing majorant {tail.str(4)} < 1 at k = 200)")
    cmp = {k: (thG ** ((k - 4) // 8) * TGold) / (TH * min(theta ** ((k - 4) // N), thG ** ((k - 4) // 8)))
           for k in (4, 12, 13, 20, 22, 52, 104)}
    print("   DH.5: DG.4 bound / DH.5 bound at k = " + ", ".join(f"{k}: {v.str(4)}" for k, v in cmp.items()))
    worse = [k for k in range(4, 200) if k not in nine]
    RESULT.update(mu=mu, worst=worst, worstG=worstG, nine=nine, cmp=cmp, ratio=ratio)


def main() -> int:
    t0 = time.time()
    print("Theorem DH: block Lt^9 under DG.1 and an asymmetric factor lemma for the f_0 term. Not thrust.", flush=True)
    main_bound()
    theta_prime()
    consequences()
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(dg.FAILS) + len(df.FAILS) + len(cf.FAILS) + len(dd.FAILS) + len(de.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DH checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
