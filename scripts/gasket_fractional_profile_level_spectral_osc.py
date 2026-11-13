#!/usr/bin/env python3
r"""Theorem DX: the exact rate of osc log R_m, proved by a non-cone (spectral) route (DW's observation "consecutive ratios
osc_(m+1)/osc_m in [0.17239, 0.17276] for m = 6..11 match spectral 0.1725: observed, not a proved rate").

Setting (GASKET.md, Theorems CE-CG, CN, CO, DO, DT-DW): s in X = [sqrt5/2, 5/2], (Lt f)(s) = wt1 f(sigma_1 s) + wt2 f(sigma_2 s)
on H^2(D), D = {|s - c| < R}, c = 83/50, R = 8/5, e_n = u^n, u = (s - c)/R; F_0 = Lt^4 G (CN.1), F_m = Lt^m F_0,
R_m = F_(m+1)/F_m, osc log R_m = log(sup_X R_m / inf_X R_m).  CG.3/CN.2: Lt has exactly the eigenvalues r, l_2, l_3 in
|z| >= 1e-3, all simple; Pi_j their rank-one Riesz projections, Pi_in = (1/2 pi i) oint_(|z| = 1e-3) (z - Lt)^-1 dz.
A_j = Pi_j F_0 (j = 1 for r), q = l_2/r < 0, phi = A_2/A_1, psi = A_3/A_1.

Theorem DX (checked here; exit 1 on any failure):
  DX.1 [certified; CO.1-CO.4 re-run (Taylor data, Krawczyk eigenvectors of M_80, exact projection formula, second-resolvent
       error), CN.2's four circles re-run; arb on 1024 cells of X, value + derivative] A_j = A_j^N + d_j with
       A_j^N = (w.f)/(w.v) (Lt v)/mu_N explicit and ||d_j||_(H^2) <= rad K K^N eps_80 S, |d_j(s)| <= ||d_j||/sqrt(1-u^2).
       On X: A_1 > 0 (inf A_1 >= 0.0028742), phi^N is strictly decreasing ((phi^N)' <= -0.03496 on every cell) and
       |phi - phi^N| <= 8.9e-9, so osc_X phi in [0.09528337, 0.09528341] (endpoint values), sup_X |phi| <= 0.06847674
       (at 5/2), sup_X |psi| <= 0.40617.
  DX.2 [proved from DX.1, CN.2] For every m >= 0, with a_m = |q|^m sup|phi|, e_m = |l_3/r|^m sup|psi|
       + (1e-3/r)^m ev 1e-3 K_in S / inf A_1, eta_m = e_m/(1 - a_m - e_m) (whenever a_m + e_m < 1):
           L_m |q|^m <= osc log R_m <= U_m |q|^m,
           U_m = (1+|q|) osc phi^+ / ((1 - |q| a_m)(1 - a_m)) + 2 (eta_m + eta_(m+1)) / |q|^m,
           L_m = (1+|q|) osc phi^- / ((1 + |q| a_m)(1 + a_m)) - 2 (eta_m + eta_(m+1)) / |q|^m.
       (F_m = r^m A_1 (1 + q^m phi + eps_m), |eps_m| <= e_m; log R_m - log r = g(q^m phi) + O(eta), g(x) = log(1+qx) - log(1+x),
       g'(x) = (q - 1)/((1+qx)(1+x)), mean value theorem.)  U_m is nonincreasing and L_m nondecreasing in m (every
       ingredient is monotone since |l_3|, 1e-3 < |l_2|), L_m > 0 for m >= M0X = 5, and both tend to (1+|q|) osc phi
       = 0.1117245...  Hence osc log R_m / |q|^m -> Omega := (1+|q|) osc_X phi, the exact rate is |l_2|/r = 0.1725476...,
       and for every m >= M0X:  L_(M0X) |q|^m <= osc log R_m <= U_(M0X) |q|^m,
       osc_(m+1)/osc_m in [|q| L_(m+1)/U_m, |q| U_(m+1)/L_m] (-> |q|).
  DX.3 [consistency, stated values] DX.2's intervals meet DV.1's (m = 5..7) and DW.1's (m = 6..12) enclosures.
  DX.4 [proved from DX.2, CF.3/DO.3/DV.2's inclusion r in [inf R_m, sup R_m] (as stated), DV/DW's constants as stated]
       |log(s_(k+1)/(r s_k))| <= U_(k-4) |q|^(k-4) for every k >= 9; it is strictly below DW's B_W(k) for every k >= 10
       (factors 1.0036 .. 1.0012 for k = 10..16, 21.37 at k = 17, 3.5e13 at k = 103; exact rationals for k < KMX, beyond
       by the 8-periodicity of G_U (as in DV.3) and |q|^8 < theta_8), and not at k = 9 (0.9449: DV.1's row stays); also
       below min(B_DO, B_U) for every k >= 17.  So B_X(k) = min(B_W(k), U_(k-4)|q|^(k-4)) < B_W(k) for every k >= 10.
  DX.5 [DERIVED from DX.2, DO.4 and CS.2 (as stated); kept limit of the route] Every CF.3-form level-sum bound at
       k = 4 + m is >= osc log R_m >= L_m |q|^(k-4), while CS.2's sharp bound is D_4 |q|^k, D_4 = 14.59509785; so for every
       k >= 4 + M0X every CF.3-form bound exceeds CS.2's by the factor >= L_(M0X) |q|^-4 / D_4 = 8.0705 (-> >= 8.6358 as
       k -> infinity): the cone/direct route now has the spectral rate but cannot reach CS's constant.
Kept limits: the rate is for the orbit F_m = Lt^m F_0 (DW's ratios), not a cone contraction rate (DP.3's cap 0.3170757 and
DQ's K_(2,6) limit concern cone diameters; theta_8^(1/8) <= 0.2380 is unchanged as a cone statement); constants rest on
CG.3's eigenvalue count and enclosures as stated.  No Psi_info, Stage 2 still blocked.  Stage 1 untouched; no W, 0.08,
0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction as Fr

from flint import acb, acb_mat, arb, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_coefficient as cn      # noqa: E402  (CN.1, CN.2)
import gasket_fractional_profile_level_third as co            # noqa: E402  (CO.1-CO.4)
import gasket_fractional_profile_level_spectrum as cg         # noqa: E402  (CG.1, matrix, circle bounds)
import gasket_fractional_profile_level_birkhoff as bk         # noqa: E402  (CO.1/CO.2's Taylor data, as DT-DW)

FAILS: list[str] = []
RESULT: dict = {}
PREC = 320
N = cn.N_TRUNC
NC = int(os.environ.get("DX_CELLS", "1024"))
M0X = 5                         # first m with L_m > 0 (checked)
MSHOW = 16
KMX = 400                       # k checked directly in DX.4
# stated enclosures of osc log R_m (DV.1, m <= 7, and DW.1, m = 6..12; DW's recomputed rows m = 6, 7), consistency only
DV_OSC = {m: (cn_lo, cn_hi) for m, (cn_lo, cn_hi) in enumerate(zip(
    ('0.063386', '0.036827', '0.0041515', '0.00060778', '0.00010040', '1.7143e-5', '2.9507e-6', '5.0884e-7'),
    ('0.063415', '0.036883', '0.0041590', '0.00060873', '0.00010052', '1.7204e-5', '2.9652e-6', '5.1065e-7')))}
DW_OSC = {6: ('2.95072e-6', '2.95167e-6'), 7: ('5.08846e-7', '5.09024e-7'), 8: ('8.77884e-8', '8.78161e-8'),
          9: ('1.51472e-8', '1.51524e-8'), 10: ('2.61359e-9', '2.61600e-9'), 11: ('4.50968e-10', '4.51126e-10'),
          12: ('7.78134e-11', '7.79092e-11')}
CS_D4 = "14.59509786"           # CS.2: D_4 = 14.59509785 +- 5.9e-10 (upper, as stated)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


_a = cn._a


# ----------------------------------------------------------------- DX.1: the spectral pieces of F_0
def spectral_data():
    """CO.1/CO.2's Taylor data, CN.1's S, CO.3's eigenvectors and CO.4's error for r, l_2, l_3, CN.2's inner circle."""
    t0 = time.time()
    S_rho, f = bk.taylor_data()                 # CN.1 box checks, CO.1, CO.2 (re-run)
    S = cn.circle_sup(1024)                     # CN.1: ||F_0||_H2 <= S
    cn.spectral_layout()                        # CN.2: the circles separate r, l_2, l_3 and |z| < 1e-3 (CG.3 as stated)
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
            check(False, f"DX.1 (CO.3): Krawczyk eigenvector enclosure for {name} failed")
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
              f"DX.1 (CO.3/CO.4 re-run, {name}): eigenvectors of M_{N} at mu_N = {muN.str(12)} (in CG's disk), circle "
              f"|z - {cen[:9]}| = {rad}: K^N <= {KN.str(5)}, K <= {KT.str(5)} ({ns} samples); ||Pi_{name} F_0 - "
              f"Pi_{name}^N F_0||_H2 <= rad K K^N eps_N S = {err.str(4)}")
        parts[name] = (v, coef, err)
    ctx.prec = cg.PREC_RES
    okc, KN, ns = cg.circle_bound(Mc, Id, acb(0), cn.RHO_IN, eps, beta)
    K_in = KN / (1 - KN * eps)
    check(okc, f"DX.1 (CN.2 re-run): |z| = 1/1000: sup ||(z - Lt)^-1|| <= K_in = {K_in.str(6)} ({ns} samples), so "
               f"||Lt^m Pi_in F_0||_H2 <= 1e-3^(m+1) K_in S for every m >= 0")
    ctx.prec = PREC
    print(f"   spectral data {time.time() - t0:.0f} s", flush=True)
    return dict(parts=parts, S=S, K_in=K_in, eps=eps)


def lt2(v, y):
    """(Lt v)(y) and its derivative in y for the polynomial v(s) = sum v_n ((s - c)/R)^n (arb, real y)."""
    c, R = _a(cn.C_CEN), _a(cn.R_RAD)
    q4, h, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    w1 = (y + h) / (2 * y * (s7 - y))
    w2 = (y - h) / (2 * y * (y + s7))
    w1d = w1 * (1 / (y + h) - 1 / y + 1 / (s7 - y))
    w2d = w2 * (1 / (y - h) - 1 / y - 1 / (y + s7))

    def poly(s):
        u = (s - c) / R
        a = arb(0)
        b = arb(0)
        for t in reversed(v):
            b = b * u + a
            a = a * u + t
        return a, b / R
    s1, s2 = (q4 - y).sqrt(), (q4 + y).sqrt()
    p1, p1d = poly(s1)
    p2, p2d = poly(s2)
    return w1 * p1 + w2 * p2, w1d * p1 - w1 * p1d / (2 * s1) + w2d * p2 + w2 * p2d / (2 * s2)


def hull(lo, hi):
    """Ball containing [lo.lower(), hi.upper()] (outward)."""
    a, b = arb(lo.lower()), arb(hi.upper())
    m = (a + b) / 2
    return arb(m.mid(), ((b - a) / 2).upper() + m.rad())


def cells_X(sd):
    """DX.1 on X: A_1 > 0, phi strictly decreasing, osc phi, sup|phi|, sup|psi| (arb, NC cells, endpoints outward)."""
    ctx.prec = PREC
    parts = sd["parts"]
    c, R = _a(cn.C_CEN), _a(cn.R_RAD)
    a0, b0 = arb(5).sqrt() / 2, arb(5) / 2
    lo, hi = arb(a0.lower()), arb(b0.upper())
    hstep = (hi - lo) / NC

    def ev(J):                                   # ||point evaluation|| on H^2(D) at s in J
        u = arb(abs((J - c) / R).upper())
        return 1 / (1 - u * u).sqrt()

    def A(name, J):
        v, coef, err = parts[name]
        val, der = lt2(v, J)
        return coef * val, coef * der, ev(J) * err

    infA1 = None
    sup_psi = arb(0)
    dphi_max = None
    dphi_bound = arb(0)
    mono = True
    for i in range(NC):
        J = hull(lo + hstep * i, lo + hstep * (i + 1))
        A1, A1d, d1 = A("r", J)
        A2, A2d, d2 = A("l2", J)
        A3, _, d3 = A("l3", J)
        a1lo = A1 - d1
        if not bool(a1lo > 0):
            check(False, f"DX.1: A_1 > 0 not certified on cell {i}")
            return None
        lowA1 = arb(a1lo.lower())
        infA1 = lowA1 if infA1 is None else arb(min(infA1, lowA1).lower())
        # phi^N = A_2^N / A_1^N: derivative enclosure on J
        dphi = (A2d * A1 - A2 * A1d) / (A1 * A1)
        mono = mono and bool(dphi < 0)
        dphi_max = arb(dphi.upper()) if dphi_max is None else arb(max(dphi_max, arb(dphi.upper())).upper())
        # |phi - phi^N| <= (|d_2| + |phi^N| |d_1|) / (A_1^N - |d_1|) on J
        dphi_bound = arb(max(dphi_bound, ((d2 + abs(A2 / A1) * d1) / a1lo)).upper())
        sup_psi = arb(max(sup_psi, (abs(A3) + d3) / a1lo).upper())
    # endpoints: phi^N(a), phi^N(b) at the balls of sqrt5/2 and 5/2
    pa = lt2(parts["l2"][0], a0)[0] * parts["l2"][1] / (lt2(parts["r"][0], a0)[0] * parts["r"][1])
    pb = lt2(parts["l2"][0], b0)[0] * parts["l2"][1] / (lt2(parts["r"][0], b0)[0] * parts["r"][1])
    oscN = pa - pb
    osc_lo = arb((oscN - 2 * dphi_bound).lower())
    osc_hi = arb((oscN + 2 * dphi_bound).upper())
    sup_phi = arb((arb(max(abs(pa), abs(pb)).upper()) + dphi_bound).upper())
    check(mono and bool(infA1 > 0),
          f"DX.1 (arb, {NC} cells, value + derivative): A_1 >= {infA1.str(7)} > 0 on X; (phi^N)' <= {dphi_max.str(5)} < 0 on "
          f"every cell, so phi^N is strictly decreasing; |phi - phi^N| <= {dphi_bound.str(3)} on X")
    check(bool(osc_lo > 0),
          f"DX.1: phi(sqrt5/2) = {pa.str(10)}, phi(5/2) = {pb.str(10)}; osc_X phi in [{osc_lo.str(10, radius=False)}, "
          f"{osc_hi.str(10, radius=False)}], sup_X |phi| <= {sup_phi.str(8)}, sup_X |psi| = sup |A_3/A_1| <= {sup_psi.str(6)}")
    ev_max = ev(hull(lo, hi))
    return dict(infA1=infA1, osc_lo=osc_lo, osc_hi=osc_hi, sup_phi=sup_phi, sup_psi=sup_psi, ev_max=ev_max,
                phi_a=pa, phi_b=pb, dphi=dphi_bound)


# ----------------------------------------------------------------- DX.2: the two-sided rate
def rate_bounds(sd, cx):
    ctx.prec = PREC
    r = arb(cn.R_MID) + arb(0, arb("1e-35"))
    l2 = arb(cn.L2_MID) + arb(0, arb(cn.L2_RAD))
    l3 = arb(cn.L3_MID) + arb(0, arb(cn.L3_RAD))
    rho = _a(cn.RHO_IN)
    q = abs(l2) / r
    t3, tin = abs(l3) / r, rho / r
    beta_in = cx["ev_max"] * rho * sd["K_in"] * sd["S"] / cx["infA1"]
    Phl, Phh = cx["osc_lo"], cx["osc_hi"]
    sp, ps = cx["sup_phi"], cx["sup_psi"]

    def am(m):
        return q ** m * sp

    def em(m):
        return t3 ** m * ps + tin ** m * beta_in

    def eta(m):
        d = 1 - am(m) - em(m)
        return em(m) / d if bool(d > 0) else None

    def UL(m):
        e0, e1 = eta(m), eta(m + 1)
        if e0 is None or e1 is None:
            return None, None
        a = am(m)
        pert = 2 * (e0 + e1) / q ** m
        U = (1 + q) * Phh / ((1 - q * a) * (1 - a)) + pert
        L = (1 + q) * Phl / ((1 + q * a) * (1 + a)) - pert
        return arb(U.upper()), arb(L.lower())

    tab = {}
    for m in range(0, 41):
        U, L = UL(m)
        if U is not None:
            tab[m] = (U, L)
    first = min((m for m in tab if bool(tab[m][1] > 0)), default=None)
    ok_order = bool(t3 < q) and bool(tin < q)
    mono = all(bool(tab[m + 1][0] <= tab[m][0]) and bool(tab[m + 1][1] >= tab[m][1]) for m in range(M0X, 40) if m in tab)
    Om_lo, Om_hi = (1 + q) * Phl, (1 + q) * Phh
    check(first == M0X and ok_order and mono,
          f"DX.2: |l_3|/r = {t3.str(6)}, 1e-3/r = {tin.str(6)} < |q| = |l_2|/r = {q.str(10)}; beta_in = ev 1e-3 K_in S / inf A_1 = "
          f"{beta_in.str(4)}; L_m > 0 first at m = {first} (expected {M0X}); U_m nonincreasing, L_m nondecreasing (checked "
          f"m = {M0X}..40, and every ingredient is monotone); limit Omega = (1+|q|) osc phi in [{Om_lo.str(10, radius=False)}, "
          f"{Om_hi.str(10, radius=False)}]")
    rows = []
    for m in range(M0X, MSHOW + 1):
        U, L = tab[m]
        qm = q ** m
        rows.append(f"m={m}: osc log R_m in [{(L * qm).str(8, radius=False)}, {(U * qm).str(8, radius=False)}] "
                    f"(L_m = {L.str(9, radius=False)}, U_m = {U.str(9, radius=False)})")
    check(True, "DX.2 two-sided: " + "; ".join(rows))
    rat = {}
    for m in range(M0X, MSHOW):
        U0, L0 = tab[m]
        U1, L1 = tab[m + 1]
        rat[m] = (arb((q * L1 / U0).lower()), arb((q * U1 / L0).upper()))
    check(all(bool(rat[m][0] > arb("0.1721")) and bool(rat[m][1] < arb("0.17296")) for m in rat if m >= 6),
          "DX.2: osc_(m+1)/osc_m in " + "; ".join(f"m={m}: [{a.str(7, radius=False)}, {b.str(7, radius=False)}]"
                                                    for m, (a, b) in rat.items() if m in (5, 6, 8, 11, 15)) + " (-> |q|)")
    return dict(q=q, tab=tab, rat=rat, Om=(Om_lo, Om_hi), first=first, beta_in=beta_in)


# ----------------------------------------------------------------- DX.3: consistency with DV.1, DW.1
def consistency(rb):
    q, tab = rb["q"], rb["tab"]
    ok = True
    msgs = []
    for src, D in (("DV.1", DV_OSC), ("DW.1", DW_OSC)):
        for m, (lo_s, hi_s) in D.items():
            if m < M0X:
                continue
            U, L = tab[m]
            lo, hi = L * q ** m, U * q ** m
            meet = bool(lo <= arb(hi_s)) and bool(hi >= arb(lo_s))
            ok = ok and meet
            if src == "DW.1":
                msgs.append(f"m={m}: width DX {((hi - lo) / lo).str(2, radius=False)} vs DW "
                            f"{((arb(hi_s) - arb(lo_s)) / arb(lo_s)).str(2, radius=False)} (relative)")
    check(ok, "DX.3 (consistency): DX.2's intervals meet DV.1's (m = 5..7) and DW.1's (m = 6..12) stated enclosures; "
              + "; ".join(msgs))
    return ok


# ----------------------------------------------------------------- DX.4/DX.5: level sums
def level_sums(rb):
    import gasket_fractional_profile_level_direct_constants as dv
    q, tab = rb["q"], rb["tab"]
    ctx.prec = 128

    up = dv._up
    qU = up(q)
    UX = {m: up(tab[m][0]) for m in tab}
    # U_m q^m for m >= 40 bounded by U_40 q^m (U nonincreasing)

    def BX(k):
        m = k - 4
        return (UX[m] if m in UX else UX[40]) * qU ** m

    th = {n: dv._up((arb(D) / 4).tanh()) for n, D in dv.DBAR_U.items()}
    KM = KMX + 20
    G = dv.knapF(sorted(th.items()), KM)
    GD = {a: dv.knapF([(n, dv._up(t)) for n, t in dv.DO_PAIRS[a]], 20) for a in (2, 100)}
    thK_dt = {m: dv._up(dv.DT_UPPER[m]) for m in range(6)}

    def B_DO(k):
        return min(dv._up(dv.DO_UPPER[a][m]) * GD[a][k - 4 - m] for m in range(0, min(5, k - 4) + 1) for a in (2, 100))

    def B_U(k):
        return min(thK_dt[m] * G[k - 4 - m] for m in range(0, min(5, k - 4) + 1))

    # DV.3's B_V with DV's stated constants, and DW.3's B_W (DW.1's stated uppers for k = 12..16)
    thK = {m: dv._up(dv.TARGET_K_UP[m]) for m in range(dv.MV + 1)}
    for m in range(6):
        thK[m] = min(thK[m], thK_dt[m])
    oscV = {m: dv._up(dv.TARGET_OSC_UP[m]) for m in range(dv.MV + 1)}
    oscW = {m + 4: dv._up(hi) for m, (lo, hi) in DW_OSC.items() if m >= 8}

    def B_V(k):
        c = [thK[m] * G[k - 4 - m] for m in range(0, min(dv.MV, k - 4) + 1)]
        if k - 4 <= dv.MV:
            c.append(oscV[k - 4])
        return min(c)

    def B_W(k):
        return min(B_V(k), oscW[k]) if k in oscW else B_V(k)

    ks = range(4 + M0X, KMX)
    facW = {k: B_W(k) / BX(k) for k in ks}
    prev = {k: (min(B_DO(k), B_U(k)) if k <= 24 else B_U(k)) for k in range(17, KMX)}
    fac = {k: prev[k] / BX(k) for k in prev}
    better = [k for k in ks if facW[k] > 1]
    strict = all(f > 1 for f in fac.values()) and better == list(range(10, KMX))
    # tail: G_U 8-periodic (G(j+8) = theta_8 G(j) for j >= J, as in DV.3) => B_W(k+8) = theta_8 B_W(k) for k >= J + 11,
    # while U_(k+4) |q|^(k+4) <= |q|^8 U_(k-4) |q|^(k-4) (U nonincreasing) and |q|^8 < theta_8: the factor grows by 8 steps
    J = None
    for j0 in range(8, KM - 16):
        if all(G[j + 8] == th[8] * G[j] for j in range(j0, j0 + 8)):
            J = j0
            break
    tail_ok = J is not None and J + 30 < KMX and qU ** 8 < th[8] and all(facW[k] > 1 for k in range(J + 11, J + 19))
    show = {k: float(facW[k]) for k in (9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20, 22, 25, 49, 103, 200) if k in facW}
    check(strict and tail_ok,
          f"DX.4 (exact rationals, DV/DW constants as stated): B_W(k) / (U_(k-4) |q|^(k-4)) = "
          + ", ".join(f"k={k}: {v:.6g}" for k, v in show.items())
          + f"; > 1 exactly for k = 10..{KMX - 1} (k = 9: DV.1's direct term stays); also < min(B_DO, B_U) for every "
          f"k = 17..{KMX - 1}; G_U 8-periodic from j = {J}, |q|^8 < theta_8, so B_X < B_W for every k >= 10")
    # DX.5: every CF.3-form bound >= osc log R_(k-4) >= L_m q^m (DO.4); CS.2's sharp bound D_4 q^k
    L0 = tab[M0X][1]
    ratio = L0 / (q ** 4 * arb(CS_D4))
    rlim = rb["Om"][0] / (q ** 4 * arb(CS_D4))
    check(bool(ratio > arb("8.07")) and bool(rlim > arb("8.63")),
          f"DX.5 (from DX.2, DO.4, CS.2 as stated): for every k >= {4 + M0X} every CF.3-form level-sum bound is >= "
          f"L_{M0X} |q|^(k-4) and exceeds CS.2's D_4 |q|^k by the factor >= L_{M0X} |q|^-4 / D_4 = {ratio.str(6)} (as k -> infinity: >= Omega^- |q|^-4 / D_4 = {rlim.str(6)})")
    return dict(fac=fac, facW=facW, J=J, ratio=ratio, rlim=rlim, show=show, better=better)


def main() -> int:
    t0 = time.time()
    print("Theorem DX: the exact rate |l_2|/r of osc log R_m (DW's observed ratios, proved by a spectral, non-cone route). "
          "Not thrust.", flush=True)
    sd = spectral_data()
    if sd is None:
        return 1
    cx = cells_X(sd)
    if cx is None:
        return 1
    rb = rate_bounds(sd, cx)
    consistency(rb)
    ls = level_sums(rb)
    RESULT.update(sd=sd, cx=cx, rb=rb, ls=ls)
    el = time.time() - t0
    RESULT["elapsed"] = el
    print(f"elapsed {el:.1f}s")
    nf = len(FAILS) + len(cn.FAILS) + len(co.FAILS) + len(cg.FAILS) + len(getattr(bk, "FAILS", []))
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DX checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
