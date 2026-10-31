#!/usr/bin/env python3
r"""Theorem DE: the joint (ratio-of-sums) oscillation for CF.2's cone, so the Birkhoff rate improves from
theta_D^(1/5) < 0.84339 (DD) to theta_E^(1/5) < 0.69306 with the SAME cone C_(47/2) and the SAME block Lt^5 (DD's kept
limits "max over words, not the true ratio of sums" and "osc log u + osc log v"; DD's open item "the true ratio-of-sums
oscillation").

Setting (GASKET.md, Theorems CE, CF, DC, DD): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
C_a = {f > 0 : log f a-Lipschitz in s}, a = 47/2, b_5 = 2.0785933 (CF.2).  CF.2/DD bound
diam_(C_a) Lt^5 C_a <= 2 log((a+b_5)/(a-b_5)) + log(sup(v/u)/inf(v/u)) and then bound the last term by
osc log u + osc log v <= 2 V (V_5 = 1.0426526 in CF.2, V_D = 0.73438 in DD).

Theorem DE (checked here; exit 1 on any failure):
  DE.1 [proved; exact] (separation and the line DP.)  For u = Lt^5 f, v = Lt^5 g (f, g in C_a)
           log(sup(v/u)/inf(v/u)) = sup_(x,y) [log(v(x)/v(y)) + log(u(y)/u(x))] <= W := sup_(x,y) [rho(x,y) + rho(y,x)],
       rho(x,y) = sup_(f in C_a) log(Lt^5 f(x)/Lt^5 f(y)) (the sups over f and g separate for fixed x, y).  For fixed
       points, Lt^5 f(x)/Lt^5 f(y) = sum_k A_k u_k / sum_k B_k u_k over the 64 sorted points s_1 < ... < s_64
       (sigma_w x, sigma_w y), u_k = f(s_k), and on a line the a-Lipschitz condition on log u is equivalent to the
       consecutive conditions u_(k+1)/u_k in [e^(-a d_k), e^(a d_k)].  Hence sup ratio <= t iff V_1(t) <= 0, where
           V_64 = c_64,  V_k = c_k + max(e^(a d_k) V_(k+1), e^(-a d_k) V_(k+1)),  c_k = A_k - t B_k
       (exact dynamic programme; V_k is increasing in c_k, V_(k+1) and d_k, so cell enclosures (sup A, inf B,
       an upper bound for each consecutive distance in any fixed order) give a valid upper bound).
  DE.2 [certified; adaptive cover of X_s x X_s by cell pairs (all pairs on the 128-grid, pairs above 0.2904 split
       8 x 8 up to the 65536-grid), arb tables as in DD.2, float DP with outward rounding after every operation and
       1e-12 pads on exp] W <= W_E < 0.29040 and, with CF.2's b_5,
           Delta_E = 2 log((a+b_5)/(a-b_5)) + W_E < 0.64513,  theta_E = tanh(Delta_E/4) < 0.15990,
           theta_E^(1/5) < 0.69306  (DD: 1.82349, 0.42672, 0.84339; CF.2: 2.44004, 0.54413, 0.88541).
  DE.3 [certified; arb] floor: for the explicit bang-bang pair (f, g) of slopes +-a(1 - 1e-6) produced by the DP at the
       maximising points, log(sup(v/u)/inf(v/u)) >= W_low > 0.29009, so no bound of this form goes below W_low; W_E/W_low
       <= 1.0011 (the separated bound 2 V_D was 5.06 times W_low).  Also sup_f osc log Lt^5 f >= V_low' > 0.66 (DD.4's
       floor was 0.51383, DD's V_D = 0.73438), so DD's oscillation bound is within 1.104 of the truth.
  DE.4 [proved from DE.2 with CF.1, CF.3, CF.5, DC.1-DC.3 as stated, as DD.3 was]: every statement of CF.3, CF.5, DC,
       DD.3 holds with theta_E: every C^1 eigenvalue mu != r has |mu| <= theta_E^(1/5) r < 0.29367; for every k >= 4
           |log(s_(k+1)/(r s_k))| <= theta_E^floor((k-4)/5) Theta_G,  Theta_G = 0.28154093 (DC.2),
       C_E = Theta_G theta_E^(-8/5) < 5.28923; consistency on CF.4's f = 1 check and CC's s_4..s_13.
Kept limits: same cone family, same a, same block; the cone-lemma term 2 log((a+b_5)/(a-b_5)) = 0.3547 is untouched
and is now more than half of Delta_E; the rate 0.69306 is still far from 0.1725; C_E carries theta_E^(-8/5).
Observation (floats, point grid, no claim): with the same DP, larger a and n give better per-level rates (a = 100,
n = 8: about 0.643; a = 60, n = 6: about 0.678), because W decays with n while the cone-lemma term decays like 4 b*/a.
Stage 1 untouched; no Psi_info, W (Stage-1 symbol), 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a
flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_cone as dd            # noqa: E402

cf = dd.cf
FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A = 23.5
CELLS = int(os.environ.get("DE_CELLS", "128"))
LEVELS = int(os.environ.get("DE_LEVELS", "4"))
FAC = 8
W_TARGET = os.environ.get("DE_TARGET", "0.2904")
THETA_G = "0.28154093"
W_MAX, DELTA_MAX, THETA_MAX, RATE_MAX, MU_MAX, WLOW_MIN, CE_MAX = "0.29040", "0.64513", "0.15990", "0.69306", "0.29367", "0.29009", "5.28923"
INF, NINF = math.inf, -math.inf


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _up(x):
    return np.nextafter(x, INF)


def _dn(x):
    return np.nextafter(x, NINF)


# ----------------------------------------------------------------- DE.1/DE.2 the relaxed DP on cells
_ORDER = None


def cell_tables(m: int, idx):
    """Outward arb enclosures for the cells J_p = [e_p, e_(p+1)], e_j = s0 + (s1 - s0) j/m, p in idx, for every word
    (sorted by the position of its image sigma_w(X)): Ahi >= sup_(J_p) wt_w, Blo <= inf_(J_p) wt_w, and the hull
    [IL, IH] of sigma_w(J_p) (sigma_w monotone).  Same arithmetic as DD.2's tables."""
    import itertools
    global _ORDER
    ctx.prec = 100
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    q, h, sv = arb(15) / 4, arb(1) / 2, arb(7) / 2
    words = list(itertools.product((0, 1), repeat=5))
    n = len(idx)
    AH, BL, IH, IL = (np.empty((32, n)) for _ in range(4))
    for c, j in enumerate(idx):
        e0j, e1j = s0 + (s1 - s0) * j / m, s0 + (s1 - s0) * (j + 1) / m
        for k, w in enumerate(words):
            e0, e1 = e0j, e1j
            ball = arb.union(e0, e1)
            lw = arb(0)
            for i in w:
                if i == 0:
                    lw += ((ball + h) / (2 * ball * (sv - ball))).log()
                    ball, e0, e1 = (q - ball).sqrt(), (q - e0).sqrt(), (q - e1).sqrt()
                else:
                    lw += ((ball - h) / (2 * ball * (ball + sv))).log()
                    ball, e0, e1 = (q + ball).sqrt(), (q + e0).sqrt(), (q + e1).sqrt()
            ew = lw.exp()
            AH[k, c] = math.nextafter(float(ew.upper()), INF)
            BL[k, c] = math.nextafter(float(ew.lower()), NINF)
            IH[k, c] = max(dd._fup(e0), dd._fup(e1))
            IL[k, c] = min(dd._flo(e0), dd._flo(e1))
    if _ORDER is None:
        _ORDER = np.argsort(IH[:, 0] + IL[:, 0])
    o = _ORDER
    return AH[o], BL[o], IH[o], IL[o]


def check_separated():
    _, _, IH, IL = cell_tables(1, [0])
    ok = bool(np.all(IH[:-1, 0] < IL[1:, 0]))
    check(ok, "DE.1: the 32 images sigma_w(X_s) are pairwise disjoint and ordered (outward hulls), so the 64 points of a "
              "pair (x, y) are ordered word by word")
    return ok


def _V1(t, c_x, b_y, isx, Ep, Em, outward):
    if outward:
        tb = _dn(t[None, :] * b_y)
        c = np.where(isx, c_x, -tb)
        V = c[-1].copy()
        for k in range(c.shape[0] - 2, -1, -1):
            V = _up(c[k] + np.where(V >= 0, _up(Ep[k] * V), _up(Em[k] * V)))
        return V
    c = np.where(isx, c_x, -t[None, :] * b_y)
    V = c[-1].copy()
    for k in range(c.shape[0] - 2, -1, -1):
        V = c[k] + np.where(V >= 0, Ep[k] * V, Em[k] * V)
    return V


def rho_cols(X, Y, a=A, iters=46):
    """Columns are (x-cell, y-cell) pairs; X = (AH, BL, IH, IL) of the x-cells, Y of the y-cells (each (32, C)).
    Returns R >= sup_(x in J, y in J', f in C_a) log(Lt^5 f(x)/Lt^5 f(y)), certified by the outward DP."""
    xA, _, xh, xl = X
    _, yB, yh, yl = Y
    K, C = xh.shape
    xfirst = (xh + xl) <= (yh + yl)
    H = np.empty((2 * K, C)); L = np.empty((2 * K, C)); isx = np.empty((2 * K, C), bool)
    H[0::2], H[1::2] = np.where(xfirst, xh, yh), np.where(xfirst, yh, xh)
    L[0::2], L[1::2] = np.where(xfirst, xl, yl), np.where(xfirst, yl, xl)
    isx[0::2], isx[1::2] = xfirst, ~xfirst
    cA = np.empty((2 * K, C)); bB = np.empty((2 * K, C))
    cA[0::2], cA[1::2] = xA, xA
    bB[0::2], bB[1::2] = yB, yB
    D = np.maximum(np.maximum(_up(H[1:] - L[:-1]), _up(H[:-1] - L[1:])), 0.0)
    aD = _up(a * D)
    Ep = _up(np.exp(aD) * (1 + 1e-12))
    Em = _dn(np.exp(-aD) * (1 - 1e-12))
    lo, hi = np.full(C, math.exp(-3.0)), np.full(C, math.exp(3.0))
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        pos = _V1(mid, cA, bB, isx, Ep, Em, False) > 0
        lo, hi = np.where(pos, mid, lo), np.where(pos, hi, mid)
    t = _up(hi * (1 + 1e-9))
    for _ in range(5):
        bad = _V1(t, cA, bB, isx, Ep, Em, True) > 0
        if not bad.any():
            break
        t = np.where(bad, _up(t * 1.001), t)
    if np.any(_V1(t, cA, bB, isx, Ep, Em, True) > 0):
        raise RuntimeError("outward DP certificate failed")
    return _up(np.log(t) + 1e-14)          # log rounded up (pad >> libm error for t in [e^-3, e^3])


def _sel(T, cols):
    return tuple(a[:, cols] for a in T)


def pair_sums(m, P, Q):
    """S[k] >= sup over x in J_(P[k]), y in J_(Q[k]) of rho(x,y) + rho(y,x) on the m-grid."""
    cells = sorted(set(P) | set(Q))
    pos = {c: i for i, c in enumerate(cells)}
    T = cell_tables(m, cells)
    out = np.empty(len(P))
    Rmax = -INF
    CH = 8192
    for s in range(0, len(P), CH):
        ip = [pos[c] for c in P[s:s + CH]]; iq = [pos[c] for c in Q[s:s + CH]]
        Rxy = rho_cols(_sel(T, ip), _sel(T, iq))
        Ryx = rho_cols(_sel(T, iq), _sel(T, ip))
        out[s:s + CH] = _up(Rxy + Ryx)
        Rmax = max(Rmax, float(Rxy.max()), float(Ryx.max()))
    return out, Rmax


def adaptive_W(m0=None, T=None, levels=None, fac=None):
    """W <= max over a cover of X_s x X_s by cell pairs: all pairs on the m0-grid, then every pair whose bound exceeds
    the target T is split into fac x fac sub-pairs (exactly nested grids), up to `levels` times."""
    m0 = m0 or CELLS; levels = LEVELS if levels is None else levels; fac = fac or FAC
    T = float(T if T is not None else W_TARGET)
    P = [p for p in range(m0) for q in range(m0)]
    Q = [q for p in range(m0) for q in range(m0)]
    m = m0
    kept = -INF
    stats = []
    Vrel = None
    for lev in range(levels + 1):
        S, Rm = pair_sums(m, P, Q)
        if Vrel is None:
            Vrel = Rm
        hot = S > T
        stats.append((m, len(P), int(hot.sum()), float(S.max())))
        print(f"   DE.2: grid {m}: {len(P)} cell pairs, max bound {S.max():.6f}, {int(hot.sum())} above {T}", flush=True)
        if lev == levels or not hot.any():
            kept = max(kept, float(S.max()))
            break
        if (~hot).any():
            kept = max(kept, float(S[~hot].max()))
        P2, Q2 = [], []
        for p, q in zip(np.array(P)[hot], np.array(Q)[hot]):
            for i in range(fac):
                for j in range(fac):
                    P2.append(int(p) * fac + i); Q2.append(int(q) * fac + j)
        P, Q, m = P2, Q2, m * fac
    return kept, Vrel, stats


def joint_bound():
    t0 = time.time()
    check_separated()
    Wf, Vf, stats = adaptive_W()
    Wf = float(_up(Wf + 1e-12))
    ctx.prec = PREC
    a = arb(47) / 2
    b5 = dd.b_iterate(a)
    W = arb(Wf)
    Delta = 2 * ((a + b5) / (a - b5)).log() + W
    theta = (Delta / 4).tanh()
    rate = (theta.log() / 5).exp()
    print(f"   DE.2: {time.time() - t0:.1f}s; relaxed sup osc log Lt^5 f <= {Vf:.6f} on the {CELLS}-grid (DD: 0.73438)")
    check(bool(W < arb(W_MAX)), f"DE.2: log(sup(v/u)/inf(v/u)) <= W_E = {W.str(8)} < {W_MAX} for all u, v in Lt^5 C_(47/2) "
                                f"(DD: 2 V_D = 1.46876; CF.2: 2 V_5 = 2.0853)")
    check(bool(b5 < arb("2.0785934")) and bool(b5 > arb("2.0785932")), f"DE.2: b_5 = {b5.str(10)} (CF.2's value)")
    check(bool(Delta < arb(DELTA_MAX)) and bool(theta < arb(THETA_MAX)) and bool(rate < arb(RATE_MAX)),
          f"DE.2: Delta_E = {Delta.str(8)} < {DELTA_MAX}, theta_E = {theta.str(8)} < {THETA_MAX}, "
          f"theta_E^(1/5) = {rate.str(8)} < {RATE_MAX} (DD: 0.84339)")
    cone_term = 2 * ((a + b5) / (a - b5)).log()
    RESULT.update(W=W, Vrel=arb(Vf), b5=b5, Delta=Delta, theta=theta, rate=rate, cone_term=cone_term, stats=stats)


# ----------------------------------------------------------------- DE.3 certified floor from explicit bang-bang f
def _point_problem(x: float, y: float, a: float):
    """sorted points, coefficients and consecutive distances for the exact point pair (floats, for the DP only)."""
    import itertools
    pts, A_, B_ = [], [], []
    for w in itertools.product((0, 1), repeat=5):
        for z, isx in ((x, True), (y, False)):
            b, lw = z, 0.0
            for i in w:
                if i == 0:
                    lw += math.log((b + .5) / (2 * b * (3.5 - b))); b = math.sqrt(3.75 - b)
                else:
                    lw += math.log((b - .5) / (2 * b * (b + 3.5))); b = math.sqrt(3.75 + b)
            pts.append((b, w, isx, math.exp(lw)))
    pts.sort(key=lambda r: r[0])
    return pts


def _bangbang(x: float, y: float, a: float):
    """the DP's maximiser for sup log(Lt^5 f(x)/Lt^5 f(y)): slope signs between consecutive sorted points."""
    pts = _point_problem(x, y, a)
    s = np.array([p[0] for p in pts]); d = np.diff(s)
    Av = np.array([p[3] if p[2] else 0.0 for p in pts]); Bv = np.array([0.0 if p[2] else p[3] for p in pts])
    lo, hi = math.exp(-3), math.exp(3)
    for _ in range(80):
        t = math.sqrt(lo * hi)
        c = Av - t * Bv; V = c[-1]
        for k in range(len(c) - 2, -1, -1):
            V = c[k] + (math.exp(a * d[k]) if V >= 0 else math.exp(-a * d[k])) * V
        lo, hi = (t, hi) if V > 0 else (lo, t)
    t = lo
    c = Av - t * Bv; V = c[-1]; signs = [0] * (len(c) - 1)
    for k in range(len(c) - 2, -1, -1):
        signs[k] = 1 if V >= 0 else -1
        V = c[k] + math.exp(signs[k] * a * d[k]) * V
    return pts, signs


def _eval_ratio_arb(x: arb, y: arb, pts, signs, a: arb):
    """certified Lt^5 f(x)/Lt^5 f(y) for f = exp(phi), phi piecewise linear through the sorted points with slopes
    signs * a (1 - 1e-6) (constant outside), which is (1 - 1e-6) a-Lipschitz, so f in C_a.  The points are recomputed
    in arb; their order is checked."""
    import itertools
    slope = a * (1 - arb("1e-6"))
    key = {}
    for w in itertools.product((0, 1), repeat=5):
        for z, isx in ((x, True), (y, False)):
            b, wt = z, arb(1)
            for i in w:
                if i == 0:
                    wt *= (b + arb(1) / 2) / (2 * b * (arb(7) / 2 - b)); b = (arb(15) / 4 - b).sqrt()
                else:
                    wt *= (b - arb(1) / 2) / (2 * b * (b + arb(7) / 2)); b = (arb(15) / 4 + b).sqrt()
            key[(w, isx)] = (b, wt)
    seq = [key[(p[1], p[2])] for p in pts]
    for (b0, _), (b1, _) in zip(seq, seq[1:]):
        if not bool(b0 < b1):
            raise RuntimeError("order of points not certified")
    phi = arb(0); num = arb(0); den = arb(0)
    for k, (p, (b, wt)) in enumerate(zip(pts, seq)):
        if k > 0:
            phi = phi + signs[k - 1] * slope * (b - seq[k - 1][0])
        if p[2]:
            num += wt * phi.exp()
        else:
            den += wt * phi.exp()
    return num / den


def floors():
    ctx.prec = 120
    s0, s1 = math.sqrt(5) / 2, 2.5
    a = arb(47) / 2
    cand = [(s0 + 1e-9, s1 - 1e-9), (s0 + 1e-9, 2.4), (s0 + 0.01, s1 - 1e-9)]
    best = arb(0)
    for x, y in cand:
        xa, ya = arb(x), arb(y)
        pf, sf = _bangbang(y, x, A)      # f maximises u(y)/u(x)
        pg, sg = _bangbang(x, y, A)      # g maximises v(x)/v(y)
        val = _eval_ratio_arb(ya, xa, pf, sf, a).log() + _eval_ratio_arb(xa, ya, pg, sg, a).log()
        if bool(val.lower() > best.lower()):
            best = arb(val.lower())
    # floor for the single-function oscillation at the end points
    x, y = s0 + 1e-9, s1 - 1e-9
    vb = arb(0)
    for (u, v) in ((x, y), (y, x)):
        pg, sg = _bangbang(u, v, A)
        val = _eval_ratio_arb(arb(u), arb(v), pg, sg, a).log()
        if bool(val.lower() > vb.lower()):
            vb = arb(val.lower())
    ctx.prec = PREC
    W = RESULT["W"]
    check(bool(best > arb(WLOW_MIN)) and bool(best < W),
          f"DE.3: explicit f, g in C_a give log(sup(v/u)/inf(v/u)) >= W_low = {best.str(8)} > {WLOW_MIN}; "
          f"W_E/W_low <= {(W / best).str(5)} (DD's 2 V_D / W_low = {(arb('1.46876') / best).str(4)})")
    check(bool(vb > arb("0.66")) and bool(vb < arb("0.73438")),
          f"DE.3: sup_f osc log Lt^5 f >= V_low' = {vb.str(8)} > 0.66 (DD.4: 0.51383), so DD's V_D = 0.73438 is "
          f"within {(arb('0.73438') / vb).str(4)} of the truth")
    RESULT.update(Wlow=best, Vlow2=vb)


# ----------------------------------------------------------------- DE.4 consequences and consistency
def consequences():
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    theta, rate = RESULT["theta"], RESULT["rate"]
    r_lo, r_hi = ce.certificate()
    rhi = dd._arb(r_hi)
    check(bool(rate * rhi < arb(MU_MAX)), f"DE.4: |mu| <= theta_E^(1/5) r <= {(rate * rhi).str(8)} < {MU_MAX} "
                                          f"(DD.3: < 0.3575; CF.5: < 0.3752)")
    ok = cf.pointwise_consistency(12, dict(a=dd._arb(dd.A_CONE), b1=None, theta=theta), r_lo, r_hi)
    check(ok, "DE.4 (consistency): CF.4's f = 1 check at five points for m < 12 holds with theta_E")
    ThetaG = arb(THETA_G)
    CE = arb(((ThetaG / (theta.log() * 8 / 5).exp()) * (1 + arb("1e-40"))).upper())
    import gasket_fractional_profile_pure_point as cc
    found, sums, counts, okb = cc.period_atoms(dd.KMAX, verbose=False)
    ctx.prec = PREC
    r = (dd._arb(r_lo) + dd._arb(r_hi)) / 2
    worst = arb(0)
    okd = okb
    for k in range(4, dd.KMAX):
        dev = abs((sums[k + 1] / (r * sums[k])).log())
        bnd = theta ** ((k - 4) // 5) * ThetaG
        okd = okd and bool(dev < bnd) and bool(bnd <= CE * (theta.log() * k / 5).exp())
        worst = arb(max(worst.upper(), (dev / bnd).upper()))
    check(bool(CE < arb(CE_MAX)), f"DE.4: C_E = Theta_G theta_E^(-8/5) = {CE.str(8)} < {CE_MAX} (block factor; the floor form "
                                  f"theta_E^floor((k-4)/5) Theta_G is never larger than DD.3's)")
    check(okd, f"DE.4 (consistency): on CC's s_4..s_{dd.KMAX}, |log(s_(k+1)/(r s_k))| < theta_E^floor((k-4)/5) Theta_G "
               f"<= C_E theta_E^(k/5), C_E = {CE.str(6)} (max ratio deviation/bound {worst.str(4)})")
    thD = arb("0.42672")
    gains = {k: (thD / theta) ** ((k - 4) // 5) for k in (9, 14, 24, 54)}
    print("   DE.4: DD.3 bound / DE.4 bound at k = " + ", ".join(f"{k}: {g.str(5)}" for k, g in gains.items()))
    RESULT.update(CE=CE, worst=worst, gains=gains)


def _obs_point_W(a: float, n: int, N: int = 33):
    """floats only, point grid, no claim: max over grid pairs of rho(x,y) + rho(y,x) for block n and cone a."""
    import itertools
    s = np.linspace(math.sqrt(5) / 2, 2.5, N)
    P, Lw = [], []
    for w in itertools.product((0, 1), repeat=n):
        b = s.copy(); lw = np.zeros(N)
        for i in w:
            if i == 0:
                lw += np.log((b + .5) / (2 * b * (3.5 - b))); b = np.sqrt(3.75 - b)
            else:
                lw += np.log((b - .5) / (2 * b * (b + 3.5))); b = np.sqrt(3.75 + b)
        P.append(b); Lw.append(np.exp(lw))
    P, Wt = np.array(P), np.array(Lw)
    o = np.argsort(P[:, 0]); P, Wt = P[o], Wt[o]
    ii, jj = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
    ii, jj = ii.ravel(), jj.ravel()
    X = (Wt[:, ii], Wt[:, ii], P[:, ii], P[:, ii]); Y = (Wt[:, jj], Wt[:, jj], P[:, jj], P[:, jj])
    R = rho_cols(X, Y, a=a).reshape(N, N)
    return float((R + R.T).max())


def observation():
    print("   observation (floats, 33-point grid, no claim): per-level rate tanh(Delta/4)^(1/n) with the joint term")
    for a, n in ((23.5, 5), (40.0, 5), (60.0, 6), (100.0, 8)):
        Wn = _obs_point_W(a, n)
        b = float(dd.b_iterate(arb(a), n).mid())
        D = 2 * math.log((a + b) / (a - b)) + Wn
        print(f"   a = {a:g}, n = {n}: W ~ {Wn:.4f}, cone term ~ {D - Wn:.4f}, Delta ~ {D:.4f}, "
              f"rate ~ {math.tanh(D / 4) ** (1 / n):.4f}", flush=True)


def main(obs=False) -> int:
    t0 = time.time()
    print("Theorem DE: the joint (ratio-of-sums) oscillation for CF.2's cone C_(47/2) under Lt^5: Birkhoff rate "
          "0.8434 -> 0.6931 per level. Not thrust.", flush=True)
    joint_bound()
    floors()
    consequences()
    if obs:
        observation()
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(cf.FAILS) + len(dd.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DE checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main(obs=True))
