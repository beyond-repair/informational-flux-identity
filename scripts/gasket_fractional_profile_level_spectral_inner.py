#!/usr/bin/env python3
r"""Theorem DY: DX's two-sided spectral rate pushed from m >= 5 to m >= 4 by a sharper inner remainder (DX's open item
"a sharper inner remainder (CQ/CR's lambda_4, lambda_5 and the 1e-8 circle would push m_0 below 5)").

Setting (GASKET.md, Theorems CG, CN, CO, CQ, DX): s in X = [sqrt5/2, 5/2], (Lt f)(s) = wt1 f(sigma_1 s) + wt2 f(sigma_2 s) on
H^2(D), D = {|s - c| < R}, c = 83/50, R = 8/5; F_0 = Lt^4 G (CN.1), F_m = Lt^m F_0, R_m = F_(m+1)/F_m,
osc log R_m = log(sup_X R_m / inf_X R_m); A_j = Pi_j F_0 (j = r, l_2, l_3, l_4), q = l_2/r, phi = A_2/A_1, psi = A_3/A_1.
DX put everything inside |z| = 1e-3 into one remainder with ||Lt^m Pi_in F_0|| <= 1e-3^(m+1) K_in S (K_in <= 1.0141e9).
Here that remainder is split along CQ.2's circles: Pi_in = Pi_4 + Pi_in', Pi_4 the rank-one Riesz projection at the
simple eigenvalue l_4 (circle G_4 = {|z + 2.164466e-4| = 1/20000}), Pi_in' the projection along |z| = rho' = 1/20000.

Theorem DY (checked here; exit 1 on any failure):
  DY.1 [certified; DX.1's spectral data re-run unchanged (CO.1-CO.4 at N = 80, CN.2's circles); CQ.2's two circles re-run
       at N = 80 (rigorous float resolvent bounds, CQ's FastRes); CQ.3's Krawczyk eigenvectors of M_140 at l_4 and its G_4
       circle at N = 140 re-run; CO.2's Taylor coefficients f_0..f_140 (M = 2048); arb on DY_CELLS cells of X]
       ||Lt^m Pi_in' F_0||_H2 <= rho'^(m+1) K' S, Lt^m Pi_4 F_0 = l_4^m A_4, A_4 = A_4^N + d_4 with A_4^N = (w.f)/(w.v) (Lt v)/mu
       (CO.3/CO.4's exact projection formula for Lt P_140) and ||d_4|| <= rad_4 K K^N eps_140 S; on X, sup|A_4|/inf A_1 <= beta_4.
  DY.2 [proved from DY.1 and DX.1, DX.2's argument unchanged] with e_m = |l_3/r|^m sup|psi| + |l_4/r|^m beta_4
       + (rho'/r)^m beta', beta' = ev_X rho' K' S / inf A_1, and DX.2's formulas for U_m, L_m: L_m > 0 first at m = M0Y = 4
       (DX: 5), U_m nonincreasing, L_m nondecreasing (every rate |l_3|/r, |l_4|/r, rho'/r is below |q|); so for every m >= 4
       L_4 |q|^m <= osc log R_m <= U_4 |q|^m, and for m >= 5 the constants are strictly inside DX's.
  DY.3 [consistency] DY.2's intervals meet DV.1's (m = 4..7) and DW.1's (m = 6..12) stated enclosures and lie inside DX's
       intervals recomputed in the same run (m = 5..40).
  DY.4 [exact rationals; DV/DW constants as stated] B_Y(k) = min(B_X(k), U^Y_(k-4) |q|^(k-4)) vs B_X (DX.4's bound,
       recomputed in this run): factors reported; at k = 8, 9 DV.1's direct rows stay (checked).
  DY.5 [DERIVED from DY.2, DO.4, CS.2 as stated] every CF.3-form bound exceeds CS.2's D_4 |q|^k by >= L^Y_5 |q|^-4 / D_4
       for every k >= 9 (DX.5: 8.0706) and by >= L^Y_4 |q|^-4 / D_4 for every k >= 8.
  DY.6 [kept limit, certified] with this split m_0 = 3 is not reachable: at m = 3 the inner term alone,
       2 (rho'/r)^3 beta' / |q|^3, exceeds (1+|q|) osc phi^+ (so L_3 < 0 for this bound); the lambda_4 term is not the obstruction.

Kept limits: a rate of the orbit F_m, not a cone contraction rate (DU's theta_8^(1/8) <= 0.2380 and DP.3's cap unchanged).
CR's lambda_5 and the 1e-8 circle are not used (they need N = 300). No Psi_info, Stage 2 still blocked.  Stage 1 untouched;
no W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time

from flint import acb, acb_mat, arb, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_spectral_osc as dx     # noqa: E402  (DX.1-DX.4, unchanged)
import gasket_fractional_profile_level_fourth as cq           # noqa: E402  (CQ.2, CQ.3, unchanged)

cn, co, cg, bk = dx.cn, dx.co, dx.cg, dx.bk
FAILS: list[str] = []
RESULT: dict = {}
PREC = dx.PREC
N = dx.N                       # 80
N4 = cq.N2                     # 140
NC = int(os.environ.get("DY_CELLS", "1024"))
M0Y = 4
MTOP = 40
KMX = dx.KMX
_a = cn._a


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- DY.1
def spectral_data():
    """DX.1's spectral data (r, l_2, l_3 at N = 80; K_in on |z| = 1e-3), with the Taylor data carried to f_140."""
    t0 = time.time()
    cn.tail_constants()
    cn.box_checks()
    co.big_disk()
    S_rho = co.circle_sup_rho(K=bk.ARCS)
    f = co.taylor_coeffs(S_rho, M=bk.MPTS, Nc=N4)
    if f is None:
        return None
    S = cn.circle_sup(1024)
    cn.spectral_layout()
    ctx.prec = PREC
    M = cg.build_matrix(N, prec=PREC)
    dc = cg.disk_constants()
    eps, beta = cg.tail_constants(dc, N)
    Mc = acb_mat(M)
    Id = acb_mat(N + 1, N + 1)
    for i in range(N + 1):
        Id[i, i] = 1
    eig_disk = {"r": (cn.R_MID, "1e-14"), "l2": (cn.L2_MID, cn.L2_RAD), "l3": (cn.L3_MID, cn.L3_RAD)}
    parts = {}
    for name in ("r", "l2", "l3"):
        ctx.prec = PREC
        v, muN, _ = co.eigvec(M, co.MIDS[name], N + 1)
        w, muL, _ = co.eigvec(M.transpose(), co.MIDS[name], N + 1)
        if v is None or w is None:
            check(False, f"DY.1 (DX.1 re-run): Krawczyk eigenvector enclosure for {name} failed")
            return None
        m_, d_ = eig_disk[name]
        in_disk = bool(abs(muN - arb(m_)) < arb(d_)) and bool(abs(muL - arb(m_)) < arb(d_))
        wv = sum((w[i] * v[i] for i in range(N + 1)), arb(0))
        wf = sum((w[i] * f[i] for i in range(N + 1)), arb(0))
        coef = wf / wv / muN
        ctx.prec = cg.PREC_RES
        cen, rad = cn.CIRCLES[name]
        okc, KN, ns = cg.circle_bound(Mc, Id, acb(arb(cen)), rad, eps, beta)
        KT = KN / (1 - KN * eps)
        ctx.prec = PREC
        err = _a(rad) * KT * KN * eps * S
        check(in_disk and okc and bool(err.is_finite()),
              f"DY.1 (DX.1 re-run, {name}): mu_N = {muN.str(12)}, circle |z - {cen[:9]}| = {rad}: K <= {KT.str(5)} "
              f"({ns} samples); ||Pi_{name} F_0 - Pi_{name}^N F_0||_H2 <= {err.str(4)}")
        parts[name] = (v, coef, err)
    ctx.prec = cg.PREC_RES
    okc, KN, ns = cg.circle_bound(Mc, Id, acb(0), cn.RHO_IN, eps, beta)
    K_in = KN / (1 - KN * eps)
    check(okc, f"DY.1 (CN.2 re-run, for DX's comparison bound): |z| = 1/1000: K_in <= {K_in.str(6)} ({ns} samples)")
    print(f"   DX spectral data {time.time() - t0:.0f} s", flush=True)
    # CQ.2 at N = 80: G_4 and |z| = 1/20000 (rigorous float resolvent bounds)
    ctx.prec = PREC
    circ, eps80, beta80, dc80 = cq.counting(M)
    if not all(v[0] for v in circ.values()):
        check(False, "DY.1 (CQ.2 re-run): circle certificates failed")
        return None
    Kp = arb(circ["in"][2])
    # CQ.3 at N = 140: Krawczyk eigenvectors at l_4, G_4 circle, l_4 enclosure
    M140 = cg.build_matrix(N4, prec=320)
    loc = cq.localise(M140, M, dc80)
    if loc is None:
        return None
    ctx.prec = PREC
    v4, w4, mu4, wv4 = loc["v"], loc["w"], loc["mu"], loc["wv"]
    wf4 = sum((w4[i] * f[i] for i in range(N4 + 1)), arb(0))
    coef4 = wf4 / wv4 / mu4
    err4 = _a(cq.RAD4) * arb(loc["KL"]) * arb(loc["KT"]) * loc["eps"] * S
    inner = _a(cq.RHO_P) * Kp * S
    check(bool(err4 < arb("1e-6")) and bool(inner.is_finite()),
          f"DY.1 (CQ.2/CQ.3 re-run): ||Pi_4 F_0 - Pi_4^N F_0||_H2 <= rad_4 K K^N eps_{N4} S = {err4.str(4)}; along |z| = "
          f"{cq.RHO_P}: sup ||(z - Lt)^-1|| <= K' = {Kp.str(5)}, so ||Lt^m Pi_in' F_0||_H2 <= rho'^(m+1) K' S, "
          f"rho' K' S = {inner.str(5)}; l_4 = {loc['l4'].str(15)}")
    print(f"   spectral data total {time.time() - t0:.0f} s", flush=True)
    return dict(parts=parts, S=S, K_in=K_in, eps=eps, v4=v4, coef4=coef4, err4=err4, Kp=Kp, inner=inner, l4=loc["l4"])


def cells_X(sd):
    """DX.1 on X (unchanged: A_1 > 0, phi monotone, osc phi, sup|phi|, sup|psi|) plus sup |A_4| / inf A_1."""
    dx.NC = NC
    cx = dx.cells_X(sd)
    if cx is None:
        return None
    ctx.prec = PREC
    c, R = _a(cn.C_CEN), _a(cn.R_RAD)
    a0, b0 = arb(5).sqrt() / 2, arb(5) / 2
    lo, hi = arb(a0.lower()), arb(b0.upper())
    hstep = (hi - lo) / NC
    sup4 = arb(0)
    for i in range(NC):
        J = dx.hull(lo + hstep * i, lo + hstep * (i + 1))
        u = arb(abs((J - c) / R).upper())
        ev = 1 / (1 - u * u).sqrt()
        val, _ = dx.lt2(sd["v4"], J)
        x = abs(sd["coef4"] * val) + ev * sd["err4"]
        if not x.is_finite():
            check(False, f"DY.1: A_4 enclosure not finite on cell {i}")
            return None
        xu = arb(x.upper())
        if bool(xu > sup4):
            sup4 = xu
    beta4 = arb((sup4 / cx["infA1"]).upper())
    check(bool(beta4.is_finite()) and bool(sup4 > 0),
          f"DY.1 (arb, {NC} cells): sup_X |A_4| <= {sup4.str(5)}, beta_4 = sup|A_4| / inf A_1 <= {beta4.str(5)}")
    cx = dict(cx)
    cx.update(sup4=sup4, beta4=beta4)
    return cx


# ----------------------------------------------------------------- DY.2 / DY.3 / DY.6
def tables(sd, cx):
    ctx.prec = PREC
    r = arb(cn.R_MID) + arb(0, arb("1e-35"))
    l2 = arb(cn.L2_MID) + arb(0, arb(cn.L2_RAD))
    l3 = arb(cn.L3_MID) + arb(0, arb(cn.L3_RAD))
    q = abs(l2) / r
    t3, t4, tp, tin = abs(l3) / r, abs(sd["l4"]) / r, _a(cq.RHO_P) / r, _a(cn.RHO_IN) / r
    evX, iA = cx["ev_max"], cx["infA1"]
    beta_in = evX * _a(cn.RHO_IN) * sd["K_in"] * sd["S"] / iA          # DX.2's
    beta_p = evX * sd["inner"] / iA
    beta4 = cx["beta4"]
    Phl, Phh = cx["osc_lo"], cx["osc_hi"]
    sp, ps = cx["sup_phi"], cx["sup_psi"]

    def make(em):
        def UL(m):
            a0, a1 = q ** m * sp, q ** (m + 1) * sp
            e0, e1 = em(m), em(m + 1)
            d0, d1 = 1 - a0 - e0, 1 - a1 - e1
            if not (bool(d0 > 0) and bool(d1 > 0)):
                return None
            pert = 2 * (e0 / d0 + e1 / d1) / q ** m
            U = (1 + q) * Phh / ((1 - q * a0) * (1 - a0)) + pert
            L = (1 + q) * Phl / ((1 + q * a0) * (1 + a0)) - pert
            return arb(U.upper()), arb(L.lower())
        tab = {}
        for m in range(0, MTOP + 1):
            x = UL(m)
            if x is not None:
                tab[m] = x
        return tab

    tabX = make(lambda m: t3 ** m * ps + tin ** m * beta_in)
    tabY = make(lambda m: t3 ** m * ps + t4 ** m * beta4 + tp ** m * beta_p)
    firstX = min((m for m in tabX if bool(tabX[m][1] > 0)), default=None)
    firstY = min((m for m in tabY if bool(tabY[m][1] > 0)), default=None)
    order = bool(t3 < q) and bool(t4 < q) and bool(tp < q)
    mono = all(bool(tabY[m + 1][0] <= tabY[m][0]) and bool(tabY[m + 1][1] >= tabY[m][1]) for m in range(M0Y, MTOP))
    inside = all(bool(tabY[m][0] <= tabX[m][0]) and bool(tabY[m][1] >= tabX[m][1]) for m in range(dx.M0X, MTOP + 1))
    strict = [m for m in range(dx.M0X, MTOP + 1) if bool(tabY[m][0] < tabX[m][0]) and bool(tabY[m][1] > tabX[m][1])]
    strict_ok = strict[:8] == list(range(dx.M0X, dx.M0X + 8))
    check(firstX == dx.M0X, f"DY.3: DX.2 recomputed in this run: L_m > 0 first at m = {firstX} (DX: {dx.M0X}); "
                            f"U_5 = {tabX[5][0].str(7, radius=False)}, L_5 = {tabX[5][1].str(7, radius=False)}")
    check(firstY == M0Y and order and mono,
          f"DY.2: |l_3|/r = {t3.str(6)}, |l_4|/r = {t4.str(6)}, rho'/r = {tp.str(6)} < |q| = {q.str(10)}; beta_4 = "
          f"{beta4.str(4)}, beta' = ev_X rho' K' S / inf A_1 = {beta_p.str(4)} (DX's beta_in = {beta_in.str(4)}); L_m > 0 first "
          f"at m = {firstY} (DX: {firstX}); U_m nonincreasing, L_m nondecreasing for m = {M0Y}..{MTOP} (every ingredient monotone)")
    rows = []
    for m in (4, 5, 6, 7, 8, 9, 10, 12, 16):
        UY, LY = tabY[m]
        qm = q ** m
        sx = f"; DX: [{tabX[m][1].str(7, radius=False)}, {tabX[m][0].str(7, radius=False)}]" if m in tabX and m >= 5 else ""
        rows.append(f"m={m}: osc log R_m in [{(LY * qm).str(8, radius=False)}, {(UY * qm).str(8, radius=False)}] "
                    f"(L = {LY.str(8, radius=False)}, U = {UY.str(8, radius=False)}{sx})")
    check(True, "DY.2 two-sided: " + "; ".join(rows))
    gap = {m: (float(tabX[m][0] - tabY[m][0]), float(tabY[m][1] - tabX[m][1])) for m in (5, 6, 8, 12)}
    check(inside and strict_ok,
          f"DY.3: [L^Y_m, U^Y_m] lies inside DX's [L_m, U_m] for every m = {dx.M0X}..{MTOP}, strictly for m in "
          f"{strict[0]}..{strict[-1] if strict == list(range(strict[0], strict[-1] + 1)) else strict} (U^X - U^Y, L^Y - L^X: "
          + ", ".join(f"m={m}: {a:.3e}, {b:.3e}" for m, (a, b) in gap.items()) + ")")
    # consistency with DV.1 / DW.1 stated
    ok = True
    for D in (dx.DV_OSC, dx.DW_OSC):
        for m, (lo_s, hi_s) in D.items():
            if m < M0Y:
                continue
            U, L = tabY[m]
            ok = ok and bool(L * q ** m <= arb(hi_s)) and bool(U * q ** m >= arb(lo_s))
    check(ok, "DY.3: DY.2's intervals meet DV.1's (m = 4..7) and DW.1's (m = 6..12) stated enclosures")
    rat = {}
    for m in range(M0Y, 16):
        U0, L0 = tabY[m]
        U1, L1 = tabY[m + 1]
        rat[m] = (arb((q * L1 / U0).lower()), arb((q * U1 / L0).upper()))
    check(True, "DY.2: osc_(m+1)/osc_m in " + "; ".join(f"m={m}: [{a.str(7, radius=False)}, {b.str(7, radius=False)}]"
                                                      for m, (a, b) in rat.items() if m in (4, 5, 6, 8)))
    # DY.6: m = 3 is out of reach with rho' = 1/20000
    inner3 = 2 * tp ** 3 * beta_p / q ** 3
    main3 = (1 + q) * Phh
    lam4_3 = 2 * t4 ** 3 * beta4 / q ** 3
    check(bool(inner3 > main3) and bool(lam4_3 < main3 / 100) and (3 not in tabY or bool(tabY[3][1] < 0)),
          f"DY.6 (kept limit): at m = 3, 2 (rho'/r)^3 beta' / |q|^3 >= {inner3.str(4, radius=False)} > (1+|q|) osc phi^+ = "
          f"{main3.str(6, radius=False)}, so this split gives no L_3 > 0; the lambda_4 term is only {lam4_3.str(3, radius=False)}")
    return dict(q=q, tabX=tabX, tabY=tabY, firstX=firstX, firstY=firstY, strict=strict, rat=rat, beta_p=beta_p, beta4=beta4,
                beta_in=beta_in, inner3=inner3, main3=main3)


# ----------------------------------------------------------------- DY.4 / DY.5
def _ranges(ks):
    out, i = [], 0
    while i < len(ks):
        j = i
        while j + 1 < len(ks) and ks[j + 1] == ks[j] + 1:
            j += 1
        out.append(f"{ks[i]}" if i == j else f"{ks[i]}..{ks[j]}")
        i = j + 1
    return "{" + ", ".join(out) + "}"


def level_sums(tb):
    import gasket_fractional_profile_level_direct_constants as dv
    q, tabX, tabY = tb["q"], tb["tabX"], tb["tabY"]
    ctx.prec = 128
    up = dv._up
    qU = up(q)
    UX = {m: up(tabX[m][0]) for m in tabX}
    UY = {m: up(tabY[m][0]) for m in tabY}

    def BXd(k):
        m = k - 4
        return (UX[m] if m in UX else UX[MTOP]) * qU ** m

    def BYd(k):
        m = k - 4
        return (UY[m] if m in UY else UY[MTOP]) * qU ** m

    th = {n: up((arb(D) / 4).tanh()) for n, D in dv.DBAR_U.items()}
    G = dv.knapF(sorted(th.items()), KMX + 20)
    thK_dt = {m: up(dv.DT_UPPER[m]) for m in range(6)}
    thK = {m: up(dv.TARGET_K_UP[m]) for m in range(dv.MV + 1)}
    for m in range(6):
        thK[m] = min(thK[m], thK_dt[m])
    oscV = {m: up(dv.TARGET_OSC_UP[m]) for m in range(dv.MV + 1)}
    oscW = {m + 4: up(hi) for m, (lo, hi) in dx.DW_OSC.items() if m >= 8}

    def B_V(k):
        c = [thK[m] * G[k - 4 - m] for m in range(0, min(dv.MV, k - 4) + 1)]
        if k - 4 <= dv.MV:
            c.append(oscV[k - 4])
        return min(c)

    def B_W(k):
        return min(B_V(k), oscW[k]) if k in oscW else B_V(k)

    def B_X(k):
        return min(B_W(k), BXd(k)) if k >= 4 + dx.M0X else B_W(k)

    def B_Y(k):
        return min(B_X(k), BYd(k))

    ks = range(4 + M0Y, KMX)
    fac = {k: B_X(k) / B_Y(k) for k in ks}
    fac_dir = {k: BXd(k) / BYd(k) for k in range(4 + dx.M0X, KMX)}
    fac8 = BYd(8) / B_W(8)
    fac9 = BYd(9) / B_W(9)
    better = [k for k in ks if fac[k] > 1]
    never_worse = all(fac[k] >= 1 for k in ks)
    show = {k: float(fac[k]) for k in (9, 10, 11, 12, 13, 16, 17, 20, 25, 49, 103, 200)}
    showd = {k: float(fac_dir[k]) for k in (9, 10, 11, 12, 16, 17, 25, 103)}
    ok = never_worse and fac8 > 1 and fac9 < 1 and 9 in better and 10 in better
    check(ok,
          "DY.4 (exact rationals, DV/DW constants as stated): B_X(k)/B_Y(k) = "
          + ", ".join(f"k={k}: {v:.8g}" for k, v in show.items())
          + f"; B_Y <= B_X for every k = 8..{KMX - 1}, strictly exactly for k in {_ranges(better)}; k = 9: U^Y_5 |q|^5 / B_W(9) = "
          f"{float(fac9):.6f} < 1 (DV.1's direct row 1.7204e-5 is beaten: DX's kept limit 'for k <= 9 DV.1's rows stay' "
          f"now only at k <= 8); k = 8: U^Y_4 |q|^4 / B_W(8) = {float(fac8):.5f} > 1 (DV.1's row stays); direct spectral "
          "terms U^X_(k-4)/U^Y_(k-4): " + ", ".join(f"k={k}: {v:.8g}" for k, v in showd.items()))
    CS_D4 = arb(dx.CS_D4)
    qa = tb["q"]
    r5 = tabY[5][1] / (qa ** 4 * CS_D4)
    r4 = tabY[4][1] / (qa ** 4 * CS_D4)
    r5x = tabX[5][1] / (qa ** 4 * CS_D4)
    check(bool(r5 > r5x) and bool(r4 > 1),
          f"DY.5 (from DY.2, DO.4, CS.2 as stated): every CF.3-form bound exceeds CS.2's D_4 |q|^k by >= "
          f"L^Y_5 |q|^-4 / D_4 = {r5.str(6)} for every k >= 9 (DX.5: {r5x.str(6)}) and by >= L^Y_4 |q|^-4 / D_4 = {r4.str(6)} "
          "for every k >= 8 (new)")
    return dict(fac=fac, fac_dir=fac_dir, fac8=fac8, fac9=fac9, better=better, show=show, r4=r4, r5=r5, r5x=r5x)


def main() -> int:
    t0 = time.time()
    print("Theorem DY: DX's two-sided spectral rate from m >= 4 (lambda_4 split off, CQ.2's 1/20000 circle). Not thrust.",
          flush=True)
    sd = spectral_data()
    if sd is None:
        return 1
    cx = cells_X(sd)
    if cx is None:
        return 1
    tb = tables(sd, cx)
    ls = level_sums(tb)
    el = time.time() - t0
    RESULT.update(sd=sd, cx=cx, tb=tb, ls=ls, elapsed=el)
    print(f"elapsed {el:.1f}s")
    nf = len(FAILS) + len(dx.FAILS) + len(cq.FAILS) + len(cn.FAILS) + len(co.FAILS) + len(cg.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DY checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
