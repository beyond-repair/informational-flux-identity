#!/usr/bin/env python3
r"""Theorem DF: a larger cone and a longer block for CF.2's contraction.  DE's joint (ratio-of-sums) programme, run on
the cone C_100 under the block Lt^8, certifies the per-level Birkhoff rate theta_F^(1/8) < 0.64398 (DE: theta_E^(1/5) <
0.69306 for C_(47/2), Lt^5), and CF.3/CF.5/DC.3 are carried to this cone (DE's open item "certify a larger (a, n) with the
joint term and carry CF.3/DC to that cone").

Setting (GASKET.md, Theorems CE, CF, DC, DD, DE): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
C_a = {f > 0 : log f a-Lipschitz in s}, phi(b) = max(M_1 + c_1 b, M_2 + c_2 b) (CF.2's cone step Lt C_b subset C_phi(b)).
Here a = 100, n = 8 (256 words, 512 points per pair).

Theorem DF (checked here; exit 1 on any failure):
  DF.1 [proved; exact] (i) b_1 = phi(100) = 45.6547 < 100, so Lt C_100 subset C_100; b_8 = phi^8(100) <= 1.8457144.
       (ii) DE.1 for general n: for u = Lt^8 f, v = Lt^8 g (f, g in C_a), log(sup(v/u)/inf(v/u)) <= W := sup_(x,y)
       [rho(x,y) + rho(y,x)], rho(x,y) = sup_(f in C_a) log(Lt^8 f(x)/Lt^8 f(y)), and for fixed x, y, rho is DE.1's exact
       line programme on the 512 sorted points sigma_w x, sigma_w y (nothing in DE.1 uses n = 5 or a = 47/2; the 256
       images sigma_w(X_s) are pairwise disjoint and ordered, outward hulls).  (iii) (Lipschitz reduction) Lt^8 f is in
       C_(b_8), so |rho(x,y) - rho(x',y')| <= b_8 (|x-x'| + |y-y'|): on a square of side H around (x_c, y_c),
       rho(x,y) + rho(y,x) <= rho(x_c,y_c) + rho(y_c,x_c) + 2 b_8 H.
  DF.2 [certified; arb point tables of wt_w(x), sigma_w(x) (prefix tree), DE.2's float DP with outward rounding after every
       operation and 1e-12 pads on exp; adaptive square cover: 128-grid (i <= j by symmetry), squares above 0.0445 split
       2 x 2 down to the 4096-grid] W <= W_F < 0.04451, and with CE.3's cone lemma (lambda = b_8/a)
           Delta_F = 2 log((a+b_8)/(a-b_8)) + W_F < 0.11835 (cone-lemma term 0.0738370),
           theta_F = tanh(Delta_F/4) < 0.02958,  theta_F^(1/8) < 0.64398  (DE: 0.64513, 0.15990, 0.69306).
  DF.3 [certified; arb] floor: explicit bang-bang f, g in C_100 (slopes +-a(1 - 1e-6)) at the end points give
       log(sup(v/u)/inf(v/u)) >= W_low > 0.04297, so W_F/W_low <= 1.0355.
  DF.4 [proved from DF.2 with CF.1, CE.3's cone lemma and DC.2's numbers as stated; CF.3's and CF.5's proofs with
       (47/2, 5) replaced by (100, 8)]: for f in C_b, b < 100, |log(Lt^(m+1) f/Lt^m f) - log r| <= theta_F^floor(m/8)
       Theta_100(f, Lt f); every C^1 eigenvalue mu != r has |mu| <= theta_F^(1/8) r < 0.27288 (DE.4: 0.29367); with DC.2's
       B_0 <= 1.2791803 and log(sup R_0/inf R_0) <= 0.063592701 (used as stated), Theta_100(f_0, Lt f_0) <= Theta_G' =
       2 log((100+B_0)/(100-B_0)) + 0.063592701 < 0.11477, so for every k >= 4
           |log(s_(k+1)/(r s_k))| <= theta_F^floor((k-4)/8) Theta_G' <= C_F theta_F^(k/8),  C_F = Theta_G' theta_F^(-11/8) < 14.531,
       which is below DE.4's theta_E^floor((k-4)/5) 0.28154093 for every k >= 4 except k = 9, 10, 11, 19 (both hold; the
       minimum is the bound).  Consistency: CF.4's f = 1 check with Theta_0(100) and theta_F^floor(m/8) for m < 12; CC's
       s_4..s_13 against DF.4's bound.
Kept limits: same cone family (log-Lipschitz in s); (a, n) = (100, 8) taken from DE's float observation, not optimised
(floats: n = 9, a = 200 about 0.633); the rate 0.644 is still far from the spectral 0.1725 (observation, floats, no
claim: in this family the cone-lemma term is about 4 b*/a while W grows with a, and the scan plateaus near 0.63); DC.1's
general-f constant K(beta) is not recomputed for a = 100; the DP certificate is outward-rounded floats (as DE.2).
Stage 1 untouched; no Psi_info, W (Stage-1 symbol), 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.
Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_joint as de            # noqa: E402

dd = de.dd
cf = dd.cf
FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A = 100.0
A_ARB_NUM = 100
N = 8
CELLS = int(os.environ.get("DF_CELLS", "128"))
LEVELS = int(os.environ.get("DF_LEVELS", "6"))
FAC = int(os.environ.get("DF_FAC", "2"))
W_TARGET = os.environ.get("DF_TARGET", "0.0445")
CH = int(os.environ.get("DF_CHUNK", "1024"))
THETA_G = "0.28154093"          # DC.2, used as stated
B0 = "1.2791803"                # DC.2: B_0 = max(b_0, phi(b_0)) <= 1.2791803, used as stated
LOGR0 = "0.063592701"           # DC.2: log(sup R_0 / inf R_0) <= 0.063592701, used as stated
THETA_E = "0.15990"             # DE.2
W_MAX, DELTA_MAX, THETA_MAX, RATE_MAX, MU_MAX, WLOW_MIN, TG_MAX, CF_MAX = ("0.04451", "0.11835", "0.02958", "0.64398",
                                                                           "0.27288", "0.04297", "0.11477", "14.531")
INF, NINF = math.inf, -math.inf
_up, _dn = de._up, de._dn


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- DF.2 tables by a prefix tree
_ORDER = None


def cell_tables(m: int, idx):
    """As DE's cell_tables, for words of length N: for each cell J_p = [e_p, e_(p+1)] (p in idx) and each word w
    (applied first letter first; sorted by the position of sigma_w(X_s)), Ahi >= sup_(J_p) wt_w, Blo <= inf_(J_p) wt_w
    and the hull [IL, IH] of sigma_w(J_p).  Prefixes are shared (same arb operations as DE, in the same order)."""
    global _ORDER
    ctx.prec = 100
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    q, h, sv = arb(15) / 4, arb(1) / 2, arb(7) / 2
    n = len(idx)
    K = 2 ** N
    AH, BL, IH, IL = (np.empty((K, n)) for _ in range(4))
    for c, j in enumerate(idx):
        e0j, e1j = s0 + (s1 - s0) * j / m, s0 + (s1 - s0) * (j + 1) / m
        layer = [(arb.union(e0j, e1j), e0j, e1j, arb(0))]
        for _ in range(N):
            nxt = []
            for ball, e0, e1, lw in layer:
                nxt.append((( q - ball).sqrt(), (q - e0).sqrt(), (q - e1).sqrt(),
                            lw + ((ball + h) / (2 * ball * (sv - ball))).log()))
                nxt.append(((q + ball).sqrt(), (q + e0).sqrt(), (q + e1).sqrt(),
                            lw + ((ball - h) / (2 * ball * (ball + sv))).log()))
            layer = nxt
        # layer index = binary word (first letter most significant), same as itertools.product order
        for k, (ball, e0, e1, lw) in enumerate(layer):
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
    global _ORDER
    _ORDER = None
    _, _, IH, IL = cell_tables(1, [0])
    ok = bool(np.all(IH[:-1, 0] < IL[1:, 0]))
    gap = float(np.min(IL[1:, 0] - IH[:-1, 0]))
    check(ok, f"DF.1: the {2 ** N} images sigma_w(X_s) are pairwise disjoint and ordered (outward hulls; smallest gap "
              f"{gap:.3e}), so the {2 ** (N + 1)} points of a pair (x, y) are ordered word by word")
    return ok


def point_tables(xs):
    """cell_tables on degenerate cells {x} (x an exact float): thin-input arb enclosures of wt_w(x) and sigma_w(x)."""
    global _ORDER
    ctx.prec = 100
    q, h, sv = arb(15) / 4, arb(1) / 2, arb(7) / 2
    K = 2 ** N
    n = len(xs)
    AH, BL, IH, IL = (np.empty((K, n)) for _ in range(4))
    for c, x in enumerate(xs):
        layer = [(arb(float(x)), arb(0))]
        for _ in range(N):
            nxt = []
            for b, lw in layer:
                nxt.append(((q - b).sqrt(), lw + ((b + h) / (2 * b * (sv - b))).log()))
                nxt.append(((q + b).sqrt(), lw + ((b - h) / (2 * b * (b + sv))).log()))
            layer = nxt
        for k, (b, lw) in enumerate(layer):
            ew = lw.exp()
            AH[k, c] = math.nextafter(float(ew.upper()), INF)
            BL[k, c] = math.nextafter(float(ew.lower()), NINF)
            IH[k, c] = dd._fup(b)
            IL[k, c] = dd._flo(b)
    o = _ORDER
    return AH[o], BL[o], IH[o], IL[o]


_PT_CACHE: dict = {}


def pair_sums(pts_x, pts_y):
    """S[k] >= rho(x_k, y_k) + rho(y_k, x_k) (outward DP certificate) at exact float points."""
    need = sorted(set(pts_x) | set(pts_y))
    new = [z for z in need if z not in _PT_CACHE]
    if new:
        T = point_tables(new)
        for i, z in enumerate(new):
            _PT_CACHE[z] = tuple(t[:, i] for t in T)
    out = np.empty(len(pts_x))
    for s in range(0, len(pts_x), CH):
        X = tuple(np.stack([_PT_CACHE[z][i] for z in pts_x[s:s + CH]], axis=1) for i in range(4))
        Y = tuple(np.stack([_PT_CACHE[z][i] for z in pts_y[s:s + CH]], axis=1) for i in range(4))
        Rxy = de.rho_cols(X, Y, a=A, iters=42)
        Ryx = de.rho_cols(Y, X, a=A, iters=42)
        out[s:s + CH] = _up(Rxy + Ryx)
    return out


def adaptive_W(b_lip: float):
    """Cover X_s x X_s by squares Q = [xc - H/2, xc + H/2] x [yc - H/2, yc + H/2] (centres xc = s0 + (i + 1/2) H on the
    m-grid, H = |X_s|/m, nested by FAC).  For f in C_a, log Lt^N f is b_N-Lipschitz (CF.2's cone step), so on Q
    rho(x,y) + rho(y,x) <= S(xc, yc) + 2 b_N (H + 2 eps) (eps = 1e-12 covers the float centres and end points).  Squares
    whose bound exceeds the target are split FAC x FAC.  S is symmetric, so the base level uses i <= j only."""
    T = float(W_TARGET)
    s0, L = math.sqrt(5) / 2, (5 - math.sqrt(5)) / 2
    m = CELLS
    I = [i for i in range(m) for j in range(i, m)]
    J = [j for i in range(m) for j in range(i, m)]
    kept = -INF
    stats = []
    for lev in range(LEVELS + 1):
        t0 = time.time()
        H = L / m
        xc = [s0 + (i + 0.5) * H for i in range(m)]
        S = pair_sums([xc[i] for i in I], [xc[j] for j in J])
        U = _up(S + _up(2 * b_lip * _up(H * (1 + 1e-12) + 4e-12)))
        hot = U > T
        stats.append((m, len(I), int(hot.sum()), float(S.max()), float(U.max())))
        print(f"   DF.2: grid {m}: {len(I)} squares, max centre value {S.max():.6f}, max bound {U.max():.6f}, "
              f"{int(hot.sum())} above {T} ({time.time() - t0:.1f}s)", flush=True)
        if lev == LEVELS or not hot.any():
            kept = max(kept, float(U.max()))
            break
        if (~hot).any():
            kept = max(kept, float(U[~hot].max()))
        I2, J2 = [], []
        for i, j in zip(np.array(I)[hot], np.array(J)[hot]):
            for u in range(FAC):
                for v in range(FAC):
                    I2.append(int(i) * FAC + u); J2.append(int(j) * FAC + v)
        I, J, m = I2, J2, m * FAC
    Smax = max(st[3] for st in stats)
    return kept, Smax, stats


def joint_bound():
    t0 = time.time()
    check_separated()
    ctx.prec = PREC
    a = arb(A_ARB_NUM)
    b1 = dd.b_iterate(a, 1)
    bn = dd.b_iterate(a, N)
    Wf, Smax, stats = adaptive_W(dd._fup(bn))
    Wf = float(_up(Wf + 1e-12))
    ctx.prec = PREC
    W = arb(Wf)
    cone_term = 2 * ((a + bn) / (a - bn)).log()
    Delta = cone_term + W
    theta = (Delta / 4).tanh()
    rate = (theta.log() / N).exp()
    print(f"   DF.2: {time.time() - t0:.1f}s")
    check(bool(b1 < a), f"DF.1: b_1 = phi(100) = {b1.str(8)} < 100, so Lt C_100 subset C_100 (CF.2's cone step)")
    check(bool(bn < a), f"DF.1: b_8 = phi^8(100) = {bn.str(10)} (b* = 1.6884)")
    check(bool(W < arb(W_MAX)), f"DF.2: log(sup(v/u)/inf(v/u)) <= W_F = {W.str(8)} < {W_MAX} for all u, v in Lt^8 C_100 "
                                f"(DE, for C_(47/2), Lt^5: 0.29040)")
    check(bool(Delta < arb(DELTA_MAX)) and bool(theta < arb(THETA_MAX)) and bool(rate < arb(RATE_MAX)),
          f"DF.2: cone-lemma term {cone_term.str(8)}, Delta_F = {Delta.str(8)} < {DELTA_MAX}, theta_F = {theta.str(8)} < "
          f"{THETA_MAX}, theta_F^(1/8) = {rate.str(8)} < {RATE_MAX} (DE: 0.69306)")
    check(bool(rate < arb("0.69306")), "DF.2: the per-level rate is strictly below DE's 0.69306")
    RESULT.update(W=W, Smax=Smax, b1=b1, bn=bn, Delta=Delta, theta=theta, rate=rate, cone_term=cone_term, stats=stats)


# ----------------------------------------------------------------- DF.3 floor
def _point_problem(x: float, y: float):
    import itertools
    pts = []
    for w in itertools.product((0, 1), repeat=N):
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
    pts = _point_problem(x, y)
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
    """certified Lt^8 f(x)/Lt^8 f(y) for f = exp(phi), phi piecewise linear through the sorted points with slopes
    signs * a (1 - 1e-6), constant outside: (1 - 1e-6) a-Lipschitz, so f in C_a.  Point order is certified."""
    import itertools
    slope = a * (1 - arb("1e-6"))
    key = {}
    for w in itertools.product((0, 1), repeat=N):
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
    a = arb(A_ARB_NUM)
    best = arb(0)
    for x, y in [(s0 + 1e-9, s1 - 1e-9), (s0 + 1e-9, 2.4), (s0 + 0.01, s1 - 1e-9)]:
        xa, ya = arb(x), arb(y)
        pf, sf = _bangbang(y, x, A)
        pg, sg = _bangbang(x, y, A)
        val = _eval_ratio_arb(ya, xa, pf, sf, a).log() + _eval_ratio_arb(xa, ya, pg, sg, a).log()
        if bool(val.lower() > best.lower()):
            best = arb(val.lower())
    ctx.prec = PREC
    W = RESULT["W"]
    check(bool(best > arb(WLOW_MIN)) and bool(best < W),
          f"DF.3: explicit f, g in C_100 give log(sup(v/u)/inf(v/u)) >= W_low = {best.str(8)} > {WLOW_MIN}; "
          f"W_F/W_low <= {(W / best).str(5)}")
    RESULT.update(Wlow=best)


# ----------------------------------------------------------------- DF.4 consequences and consistency
def _consistency_f1(K, theta, r_lo, r_hi, a):
    """CF.4's check with the cone C_a and block N: f = 1 in C_0, Lt 1 in C_(M_1), Theta_a(1, Lt 1) <= Theta_0(a)."""
    from fractions import Fraction as Fr
    ctx.prec = 200
    B = arb(14) / 15
    Theta0 = 2 * ((a + B) / (a - B)).log() + (arb(11) / 6).log()
    lr_lo, lr_hi = (arb(r_lo.numerator) / r_lo.denominator).log(), (arb(r_hi.numerator) / r_hi.denominator).log()
    ok, worst = True, 0.0
    for m in range(K):
        bound = theta ** (m // N) * Theta0
        for t in [Fr(-5), Fr(-3), Fr(-5, 2), Fr(-1), Fr(0)]:
            sv = (arb(t.numerator) / t.denominator + arb(25) / 4).sqrt()
            lR = (cf._Lt_iter(sv, m + 1) / cf._Lt_iter(sv, m)).log()
            dev = arb(max(abs(lR - lr_lo).upper(), abs(lR - lr_hi).upper()))
            ok = ok and bool(dev < bound)
            worst = max(worst, float(dev) / float(bound))
    return ok, Theta0, worst


def consequences():
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    theta, rate = RESULT["theta"], RESULT["rate"]
    a = arb(A_ARB_NUM)
    r_lo, r_hi = ce.certificate()
    rhi = dd._arb(r_hi)
    mu = rate * rhi
    check(bool(mu < arb(MU_MAX)), f"DF.4: |mu| <= theta_F^(1/8) r <= {mu.str(8)} < {MU_MAX} (DE.4: < 0.29367)")
    ok, Th0, w0 = _consistency_f1(12, theta, r_lo, r_hi, a)
    check(ok, f"DF.4 (consistency): CF.4's f = 1 check with Theta_0(100) = {Th0.str(6)} and theta_F^floor(m/8) at five "
              f"points for m < 12 (max used fraction {w0:.3e})")
    ctx.prec = PREC
    B = arb(B0)
    check(bool(B < a), "DF.4: DC.2's B_0 <= 1.2791803 < 100, so DC.2's cone-lemma bound applies in C_100")
    TG = 2 * ((a + B) / (a - B)).log() + arb(LOGR0)
    TG = arb(TG.upper())
    check(bool(TG < arb(TG_MAX)), f"DF.4: Theta_100(f_0, Lt f_0) <= Theta_G' = {TG.str(8)} < {TG_MAX} (DC.2 in C_(47/2): "
                                  f"0.28154093)")
    CF_ = arb((TG / (theta.log() * arb(N + 3) / N).exp()).upper())
    import gasket_fractional_profile_pure_point as cc
    found, sums, counts, okb = cc.period_atoms(dd.KMAX, verbose=False)
    ctx.prec = PREC
    r = (dd._arb(r_lo) + dd._arb(r_hi)) / 2
    thE, ThG = arb(THETA_E), arb(THETA_G)
    worst = arb(0)
    okd = okb
    for k in range(4, dd.KMAX):
        dev = abs((sums[k + 1] / (r * sums[k])).log())
        bF = theta ** ((k - 4) // N) * TG
        okd = okd and bool(dev < bF) and bool(bF < CF_ * (1 + arb("1e-20")) * (theta.log() * k / N).exp())
        worst = arb(max(worst.upper(), (dev / bF).upper()))
    check(bool(CF_ < arb(CF_MAX)), f"DF.4: C_F = Theta_G' theta_F^(-11/8) = {CF_.str(8)} < {CF_MAX} (block factor; "
                                   f"the floor form is the statement)")
    check(okd, f"DF.4 (consistency): on CC's s_4..s_{dd.KMAX}, |log(s_(k+1)/(r s_k))| < theta_F^floor((k-4)/8) Theta_G' "
               f"<= C_F theta_F^(k/8), C_F = Theta_G' theta_F^(-11/8) = {CF_.str(6)} (max ratio deviation/bound "
               f"{worst.str(4)})")
    # comparison with DE.4's bound theta_E^floor((k-4)/5) Theta_G (both bounds hold; DF.4 states their minimum)
    kF_all = None
    for k0 in range(4, 400):
        # for k >= k0: bF/bE <= (Theta_G'/Theta_G) theta_F^((k-11)/8) / theta_E^((k-4)/5), decreasing in k because
        # theta_F^(1/8) < theta_E^(1/5)
        lhs = TG / ThG * (theta.log() * (k0 - 11) / N).exp() / (thE.log() * (k0 - 4) / 5).exp()
        if bool(lhs < 1):
            kF_all = k0
            break
    mono = bool(rate < (thE.log() / 5).exp())
    worse = []
    for k in range(4, kF_all or 4):
        bE = thE ** ((k - 4) // 5) * ThG
        bF = theta ** ((k - 4) // N) * TG
        if not bool(bF < bE):
            worse.append(k)
    gains = {k: (thE ** ((k - 4) // 5) * ThG) / (theta ** ((k - 4) // N) * TG) for k in (4, 9, 12, 20, 24, 54, 104)}
    print("   DF.4: DE.4 bound / DF.4 bound at k = " + ", ".join(f"{k}: {g.str(4)}" for k, g in gains.items()))
    check(kF_all is not None and mono,
          f"DF.4: theta_F^floor((k-4)/8) Theta_G' < theta_E^floor((k-4)/5) Theta_G for every k >= 4 except k in {worse} "
          f"(checked one by one for k < {kF_all}, by the decreasing majorant for k >= {kF_all})")
    RESULT.update(mu=mu, TG=TG, CF=CF_, worst=worst, kF_all=kF_all, worse=worse, gains=gains)


def main(obs=False) -> int:
    t0 = time.time()
    print("Theorem DF: DE's joint programme on the cone C_100 under Lt^8: Birkhoff rate 0.6931 -> 0.64398 per level. "
          "Not thrust.", flush=True)
    joint_bound()
    floors()
    consequences()
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(cf.FAILS) + len(dd.FAILS) + len(de.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DF checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
