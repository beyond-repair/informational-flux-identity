#!/usr/bin/env python3
r"""Theorem DL: DK's kept limit "DJ.3's cap means (C_100, Lt^9) is exhausted; spectral 0.1725 needs another cone family or
block".  A much smaller cone, C_2, and a shorter block, Lt^5 (32 words).  DK's reduction (DK.1) applies verbatim; DJ.2's
floor construction and DK.3's outward cover are re-run for (a, n) = (2, 5).

Setting (GASKET.md, Theorems CE, CF, DC-DK): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
Lt^5 f(t) = sum_w W_w(t) f(b_w(t)) (32 words, b_w = sigma_w), C_a = {f > 0 : log f a-Lipschitz in s}, Theta_a its Hilbert
metric (CF's definition), here a = 2 (CF.2's absorbing value b* = 1.6884 < 2).

Theorem DL (checked here; exit 1 on any failure):
  DL.1 [proved; arb] (admissibility of (C_2, Lt^5)) CF.2's cone step gives phi(2) = max(14/15 + 2/sqrt5, M2 + 2 c2)
       < 1.8278 < 2 and b* = 1.6884 < 2, so Lt C_b subset C_phi(b) subset C_2 for every b <= 2 (CF.3/CF.5's hypotheses
       on a).  From the arb jets (DL.2): |(log Lt^5 f)'| <= sup_w |lambda_w| + 2 sup_w |b_w'| < 0.929 for f in C_2, so
       Lt^5 C_2 subset C_0.929 and every functional l_P (DJ.1) and m^s_t (DK.1) is positive on Lt^5 C_2.  DK.1's
       reduction holds verbatim for (a, n) = (2, 5).
  DL.2 [certified; arb jets to third order on 128 balls of X_s] sb = sup |b_w'| < 5.224e-3, G2 < 8.83,
       a + s lambda_w +- a |b_w'| > 1.074 (every C^(s,+-)_w > 0), b_w' of constant sign on each ball; float word
       tables agree with 200-bit arb at 5 centres.
  DL.3 [certified; arb] (floor, DJ.2's construction) at the near-diagonal pairs P1 = (sqrt5/2 + d, sqrt5/2),
       P2 = (5/2, 5/2 - d), d = 1e-6 (DK.1's (-,-) corner), explicit bang-bang f, g in C_2 (slopes +-2(1 - 1e-6))
       give Theta_2(Lt^5 f, Lt^5 g) >= 0.043581853 > 0.0435818.
  DL.4 [certified; DK.1's relaxation, adaptive dyadic cover, outward DP]  diam_(Theta_2) Lt^5 C_2 <= 0.043583,
       so 0.0435818 < diam Lt^5 C_2 <= 0.043583 (width 1.2e-6).
  DL.5 [certified; arb] theta_L = tanh(0.043583/4) < 0.0108954, per-level rate theta_L^(1/5) < 0.404994 (DK.4: 0.53349;
       DJ.3's cap 0.53342 for (C_100, Lt^9) does not apply to this pair); every CF.1-form rate on (C_2, Lt^5) is
       > 0.404991 (DL.3), so (C_2, Lt^5) is exhausted at 0.40499; every other C^1 eigenvalue
       |mu| <= theta_L^(1/5) r < 0.17161 (DK.4: 0.22606).
  Observation (floats, no claim, DL_OBS=1): DK.1's relaxed (-,-) corner value for (C_2, Lt^n), n = 1..9.
Kept limits: see GASKET.md Theorem DL.  No Psi_info, Stage 2 still blocked.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import itertools
import math
import os
import sys
import time

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A = 2.0
A_NUM = 2
N = 5
NBALL = 128
S0, S1 = math.sqrt(5) / 2, 2.5
EPS_PT = 1e-13                   # pad on float image points
EPS_LOG = 1e-10                  # pad on float log coefficients
DL_MAX = float(os.environ.get("DL_MAX", "0.043583"))     # acceptance threshold of the cover
DL_MAX_S = "0.043583"
DL_LOW = "0.0435818"
D_NEAR = 1e-6
SLOPE_FAC = 1e-6
TH_MAX, RATE_MAX, MU_MAX, RLOW = "0.0108954", "0.404994", "0.17161", "0.404991"
RATE_K, MU_K, CAP_J = "0.53349", "0.22606", "0.53342"     # DK.4, DJ.3 (as stated)
B_IMG = "0.929"                  # Lt^5 C_2 subset C_B_IMG
SB_MAX, G2_MAX, Q_MIN = 5.224e-3, 8.83, 1.074
MAXDEPTH = int(os.environ.get("DL_DEPTH", "22"))
BATCH = int(os.environ.get("DL_BATCH", "4000"))
OBS = os.environ.get("DL_OBS", "0") == "1"
WD = np.array(list(itertools.product((0, 1), repeat=N)))
WORDS = [tuple(w) for w in WD]
INF, NINF = np.inf, -np.inf


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- DL.1: admissibility of the cone
def admissible():
    ctx.prec = PREC
    M1, c1 = arb(14) / 15, 1 / arb(5).sqrt()          # CF.2 (exact)
    M2, c2 = arb("0.5071"), arb("0.2267")             # CF.2's certified bounds, as stated
    a = arb(A_NUM)
    phi = lambda b: arb(max((M1 + c1 * b).upper(), (M2 + c2 * b).upper()))   # noqa: E731
    bstar = M1 / (1 - c1)
    p2 = phi(a)
    check(bool(p2 < arb("1.8278")) and bool(p2 < a) and bool(bstar < a) and bool(M2 / (1 - c2) < bstar),
          f"DL.1: CF.2's cone step phi(2) = max(14/15 + 2/sqrt5, M2 + 2 c2) <= {p2.str(8)} < 1.8278 < 2, b* = 7(5+sqrt5)/30 = "
          f"{bstar.str(8)} < 2: Lt C_b subset C_phi(b) subset C_2 for every b <= 2 (CF.3/CF.5's hypotheses b* < a, phi(b) < a)")
    RESULT.update(phi2=p2, bstar=bstar)


# ----------------------------------------------------------------- DL.2: arb jets (third order) on balls of X_s
def jets(T: arb):
    """forward-mode jets of b_w, log W_w at the ball T over all 32 words (prefix tree):
    (b, b', b'', b''', log W, lambda, lambda', lambda'')"""
    half, q75, q35 = arb(1) / 2, arb(15) / 4, arb(7) / 2
    level = [(T, arb(1), arb(0), arb(0), arb(0), arb(0), arb(0), arb(0))]
    for _ in range(N):
        nxt = []
        for (b, b1, b2, b3, lw, L1, L2, L3) in level:
            for i in (0, 1):
                if i == 0:
                    w = (b + half) / (2 * b * (q35 - b))
                    g = 1 / (b + half) - 1 / b + 1 / (q35 - b)
                    gp = -1 / (b + half) ** 2 + 1 / b ** 2 + 1 / (q35 - b) ** 2
                    gpp = 2 / (b + half) ** 3 - 2 / b ** 3 + 2 / (q35 - b) ** 3
                    F = (q75 - b).sqrt(); Fp = -1 / (2 * F); Fppp = -3 / (8 * F ** 5)
                else:
                    w = (b - half) / (2 * b * (b + q35))
                    g = 1 / (b - half) - 1 / b - 1 / (b + q35)
                    gp = -1 / (b - half) ** 2 + 1 / b ** 2 + 1 / (b + q35) ** 2
                    gpp = 2 / (b - half) ** 3 - 2 / b ** 3 - 2 / (b + q35) ** 3
                    F = (q75 + b).sqrt(); Fp = 1 / (2 * F); Fppp = 3 / (8 * F ** 5)
                Fpp = -1 / (4 * F ** 3)
                nxt.append((F, Fp * b1, Fpp * b1 ** 2 + Fp * b2, Fppp * b1 ** 3 + 3 * Fpp * b1 * b2 + Fp * b3,
                            lw + w.log(), L1 + g * b1, L2 + gp * b1 ** 2 + g * b2,
                            L3 + gpp * b1 ** 3 + 3 * gp * b1 * b2 + g * b3))
        level = nxt
    return level


def constants():
    ctx.prec = 64
    t0 = time.time()
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    L = s1 - s0
    sb = G2 = lam = 0.0
    qmin = 1e300
    sgn = True
    for j in range(NBALL):
        T = (s0 + L * j / NBALL).union(s0 + L * (j + 1) / NBALL)
        for idx, (b, b1, b2, b3, lw, L1, L2, L3) in enumerate(jets(T)):
            if not (bool(b1 > 0) or bool(b1 < 0)):
                sgn = False
                continue
            sg = 1 if bool(b1 > 0) else -1
            sb = max(sb, float(abs(b1).upper()))
            lam = max(lam, float(abs(L1).upper()))
            for s in (1, -1):
                for sig in (1, -1):
                    q = A_NUM + s * L1 + sig * A_NUM * sg * b1
                    qp = s * L2 + sig * A_NUM * sg * b2
                    qpp = s * L3 + sig * A_NUM * sg * b3
                    qmin = min(qmin, float(q.lower()))
                    G2 = max(G2, float(abs(L2 + qpp / q - (qp / q) ** 2).upper()))
    check(sgn, f"DL.2: on each of {NBALL} balls of X_s every b_w' (32 words) has constant sign (so |b_w'| is smooth there)")
    check(sb < SB_MAX and G2 < G2_MAX and qmin > Q_MIN,
          f"DL.2 (arb jets): sup |b_w'| <= {sb:.6e} < {SB_MAX}, sup |(log C^(s,+-)_w)''| <= {G2:.5f} < {G2_MAX}, "
          f"a + s lambda_w +- a|b_w'| >= {qmin:.5f} > {Q_MIN} (so every C^(s,+-)_w > 0) ({time.time() - t0:.1f}s)")
    ctx.prec = PREC
    bimg = arb(lam) + 2 * arb(sb)
    check(bool(bimg < arb(B_IMG)),
          f"DL.1: |(log Lt^5 f)'| <= sup|lambda_w| + a sup|b_w'| <= {float(bimg.upper()):.6f} < {B_IMG} for f in C_2 "
          f"(weighted mean of lambda_w + b_w' (log f)'(b_w)), so Lt^5 C_2 subset C_{B_IMG} and every l_P, m^s_t is positive there")
    RESULT.update(sb=sb, G2=G2, qmin=qmin, lam=lam)
    return SB_MAX, G2_MAX


# ----------------------------------------------------------------- float word tables (DK's, for n = 5)
def words(t):
    t = np.atleast_1d(np.asarray(t, float))
    b = np.repeat(t[:, None], 2 ** N, 1)
    lw = np.zeros_like(b); L1 = np.zeros_like(b); L2 = np.zeros_like(b)
    b1 = np.ones_like(b); b2 = np.zeros_like(b)
    for k in range(N):
        z = (WD[:, k] == 0)[None, :]
        w = np.where(z, (b + .5) / (2 * b * (3.5 - b)), (b - .5) / (2 * b * (b + 3.5)))
        g = np.where(z, 1 / (b + .5) - 1 / b + 1 / (3.5 - b), 1 / (b - .5) - 1 / b - 1 / (b + 3.5))
        gp = np.where(z, -1 / (b + .5) ** 2 + 1 / b ** 2 + 1 / (3.5 - b) ** 2,
                      -1 / (b - .5) ** 2 + 1 / b ** 2 + 1 / (b + 3.5) ** 2)
        lw = lw + np.log(w); L2 = L2 + gp * b1 ** 2 + g * b2; L1 = L1 + g * b1
        F = np.where(z, np.sqrt(3.75 - b), np.sqrt(3.75 + b))
        Fp = np.where(z, -1 / (2 * F), 1 / (2 * F)); Fpp = -1 / (4 * F ** 3)
        b2 = Fpp * b1 ** 2 + Fp * b2; b1 = Fp * b1; b = F
    return b, lw, L1, L2, b1, b2


_PERM = np.argsort(words([2.0])[0][0])


def coef(t, s):
    """points b_w(t) (float order), log C^(s,+)_w, log C^(s,-)_w and their t-derivatives g"""
    b, lw, L1, L2, b1, b2 = words(t)
    out = [b[:, _PERM]]
    for sig in (1, -1):
        q = A + s * L1 + sig * A * np.abs(b1)
        out.append((lw + np.log(q))[:, _PERM])
        out.append((L1 + (s * L2 + sig * A * np.sign(b1) * b2) / q)[:, _PERM])
    return out


def float_vs_arb():
    ctx.prec = PREC
    eb = el = eg = 0.0
    for tc in (S0, 1.4, 1.8, 2.2, S1):
        J = jets(arb(tc))
        for s in (1, -1):
            p, cp, gp, cm, gm = coef([tc], s)
            for idx in range(2 ** N):
                b, b1, b2, b3, lw, L1, L2, L3 = J[idx]
                k = int(np.nonzero(_PERM == idx)[0][0])
                sg = 1 if bool(b1 > 0) else -1
                eb = max(eb, abs(float((b - p[0, k]).mid())))
                for sig, cc, gg in ((1, cp, gp), (-1, cm, gm)):
                    q = A_NUM + s * L1 + sig * A_NUM * sg * b1
                    el = max(el, abs(float((lw + q.log() - cc[0, k]).mid())))
                    eg = max(eg, abs(float((L1 + (s * L2 + sig * A_NUM * sg * b2) / q - gg[0, k]).mid())))
    check(eb < 1e-14 and el < 1e-12 and eg < 1e-10,
          f"DL.2: float tables vs 200-bit arb at 5 centres: |b| err {eb:.1e}, |log C| err {el:.1e}, |g| err {eg:.1e} "
          f"(pads {EPS_PT:.0e} on points, {EPS_LOG:.0e} on log coefficients)")


# ----------------------------------------------------------------- DK's line programme and box bound (a = 2)
def _merge(p1, c1, p2, c2):
    P = np.concatenate([p1, p2], 1)
    An = np.concatenate([c1, np.zeros_like(c2)], 1)
    Bd = np.concatenate([np.zeros_like(c1), c2], 1)
    o = np.argsort(P, axis=1, kind="stable")
    return (np.take_along_axis(P, o, 1), np.take_along_axis(An, o, 1), np.take_along_axis(Bd, o, 1))


def _rho_up(P, An, Bd):
    """certified upper bound for sup_(h in C_a) log(sum An h(x) / sum Bd h(x)) (An upper, Bd lower coefficients,
    float points P with |x - P| <= EPS_PT); DP with gaps d_k + 2 EPS_PT, outward rounding after every operation."""
    M, K = P.shape
    d = np.nextafter(np.maximum(np.diff(P, axis=1), 0) + 2 * EPS_PT, INF)
    ad = np.nextafter(A * d, INF)
    eu, ed = np.exp(ad), np.exp(-ad)
    EU = np.nextafter(eu * (1 + 1e-15), INF)
    ED = np.nextafter(ed * (1 - 1e-15), NINF)

    def V(t):
        c = An - t[:, None] * Bd
        v = c[:, -1].copy()
        for k in range(K - 2, -1, -1):
            v = c[:, k] + np.where(v >= 0, eu[:, k], ed[:, k]) * v
        return v
    t = An.sum(1) / Bd.sum(1)
    for _ in range(3):
        c = An - t[:, None] * Bd
        v = c[:, -1].copy(); sg = np.empty((M, K - 1), bool)
        for k in range(K - 2, -1, -1):
            pos = v >= 0; sg[:, k] = pos
            v = c[:, k] + np.where(pos, eu[:, k], ed[:, k]) * v
        lg = np.concatenate([np.zeros((M, 1)), np.cumsum(np.where(sg, ad, -ad), 1)], 1)
        lg -= lg.max(1, keepdims=True)
        h = np.exp(lg)
        t = np.maximum(t, (An * h).sum(1) / (Bd * h).sum(1))
    lo = np.log(t); hi = lo + 0.02
    while True:
        bad = V(np.exp(hi)) > 0
        if not bad.any():
            break
        hi = np.where(bad, hi + 0.05, hi)
    for _ in range(30):
        mid = (lo + hi) / 2
        pos = V(np.exp(mid)) > 0
        lo = np.where(pos, mid, lo); hi = np.where(pos, hi, mid)
    tup = np.exp(hi) * (1 + 1e-12)
    todo = np.ones(M, bool)
    for _ in range(40):
        cc = np.nextafter(An - np.nextafter(tup[:, None] * Bd, NINF), INF)
        v = cc[:, -1].copy()
        for k in range(K - 2, -1, -1):
            f = np.where(v >= 0, EU[:, k], ED[:, k])
            v = np.nextafter(cc[:, k] + np.nextafter(f * v, INF), INF)
        todo = v > 0
        if not todo.any():
            break
        tup = np.where(todo, tup * (1 + 1e-9), tup)
    if todo.any():
        raise RuntimeError("outward DP certificate failed")
    return np.nextafter(np.log(tup) + 1e-14, INF)


def box_bound(c1, h1, s1, c2, h2, s2, sb, G2):
    """certified upper bound of sup over t1 in [c1 +- h1], t2 in [c2 +- h2] of rho~(Q1,Q2) + rho~(Q2,Q1)"""
    sides = []
    for c, h, s in ((c1, h1, s1), (c2, h2, s2)):
        p, lp, gp, lm, gm = coef(c, s)
        gall = np.concatenate([gp, gm], 1)
        mu = (gall.max(1) + gall.min(1)) / 2
        base = (A * sb + G2 * h / 2)[:, None] * h[:, None]
        up = lp + h[:, None] * np.abs(gp - mu[:, None]) + base + EPS_LOG
        dn = lm - h[:, None] * np.abs(gm - mu[:, None]) - base - EPS_LOG
        sh = up.max(1, keepdims=True)
        Au = np.nextafter(np.exp(up - sh) * (1 + 1e-15), INF) + 1e-300
        Bl = np.nextafter(np.exp(dn - sh) * (1 - 1e-15), NINF)
        sides.append((p, Au, Bl))
    (p1, A1, B1), (p2, A2, B2) = sides
    r12 = _rho_up(*_merge(p1, A1, p2, B2))
    r21 = _rho_up(*_merge(p2, A2, p1, B1))
    return np.nextafter(r12 + r21, INF)


def cover(sb, G2, log=lambda m: print(m, flush=True)):
    t0 = time.time()
    L = S1 - S0
    stats = {}
    Dmax = 0.0
    nbox = 0
    for s1, s2 in ((1, 1), (-1, -1), (1, -1)):
        sym = s1 == s2
        n0 = 8
        H = L / n0
        q = [(i, j) for i in range(n0) for j in range(n0) if (not sym or i <= j)]
        depth = 0
        acc = 0
        smax = 0.0
        while q:
            if depth > MAXDEPTH:
                raise RuntimeError(f"cover did not close at depth {MAXDEPTH} for signs {(s1, s2)}")
            idx = np.array(q, dtype=np.int64)
            c1 = S0 + (idx[:, 0] + 0.5) * H
            c2 = S0 + (idx[:, 1] + 0.5) * H
            hh = np.full(len(idx), H / 2 * (1 + 1e-12))
            Dv = np.empty(len(idx))
            for k in range(0, len(idx), BATCH):
                sl = slice(k, k + BATCH)
                Dv[sl] = box_bound(c1[sl], hh[sl], s1, c2[sl], hh[sl], s2, sb, G2)
            ok = Dv <= DL_MAX
            acc += int(ok.sum())
            if ok.any():
                smax = max(smax, float(Dv[ok].max()))
            nq = []
            for (i, j) in idx[~ok]:
                for a_ in (0, 1):
                    for b_ in (0, 1):
                        ii, jj = 2 * i + a_, 2 * j + b_
                        if sym and ii > jj:
                            continue
                        nq.append((ii, jj))
            log(f"   signs {(s1, s2)} depth {depth}: {len(idx)} boxes, {int(ok.sum())} accepted, max accepted {smax:.8f}, "
                f"max D_box {Dv.max():.6f} ({time.time() - t0:.0f}s)")
            q = nq
            H /= 2
            depth += 1
            nbox += len(idx)
        stats[(s1, s2)] = (acc, depth, smax)
        Dmax = max(Dmax, smax)
    return Dmax, stats, nbox, time.time() - t0


# ----------------------------------------------------------------- DL.3: DJ.2's floor construction for (a, n) = (2, 5)
def _images(z: float):
    B, Wt = [], []
    for w in WORDS:
        b, wt = z, 1.0
        for i in w:
            if i == 0:
                wt *= (b + .5) / (2 * b * (3.5 - b)); b = math.sqrt(3.75 - b)
            else:
                wt *= (b - .5) / (2 * b * (b + 3.5)); b = math.sqrt(3.75 + b)
        B.append(b); Wt.append(wt)
    return np.array(B), np.array(Wt)


def _functional(x, y):
    """coefficients of l_(x,y)(Lt^5 f) = e^(a|x-y|) Lt^5 f(y) - Lt^5 f(x) in f, normalised"""
    By, Wy = _images(y); Bx, Wx = _images(x)
    c = np.concatenate([math.exp(A * abs(x - y)) * Wy, -Wx])
    return np.concatenate([By, Bx]), c / c.sum()


def rho_float(P1, P2):
    """DJ.1's signed ratio programme (floats; chooses f only): sup_f log(l_P1(Lt^5 f)/l_P2(Lt^5 f)) over
    f = exp(phi), phi a(1 - 1e-6)-Lipschitz, with the bang-bang maximiser's knots and slope signs."""
    p1, c1 = _functional(*P1); p2, c2 = _functional(*P2)
    s = np.concatenate([p1, p2]); Av = np.concatenate([c1, 0 * c2]); Bv = np.concatenate([0 * c1, c2])
    o = np.argsort(s, kind="stable"); s, Av, Bv = s[o], Av[o], Bv[o]
    d = np.diff(s); ap = A * (1 - SLOPE_FAC)
    eu, ed = np.exp(ap * d), np.exp(-ap * d)

    def V(t, signs=None):
        c = Av - t * Bv; v = c[-1]
        for k in range(len(c) - 2, -1, -1):
            sg = 1 if v >= 0 else -1
            v = c[k] + (eu[k] if sg > 0 else ed[k]) * v
            if signs is not None:
                signs[k] = sg
        return v
    lo, hi = math.exp(-1), math.exp(1)
    for _ in range(60):
        t = math.sqrt(lo * hi)
        lo, hi = (t, hi) if V(t) > 0 else (lo, t)
    signs = [0] * (len(s) - 1)
    V(lo, signs)
    return math.log(lo), s, signs


def _make_f(s, signs):
    slope = arb(A_NUM) * (1 - arb(SLOPE_FAC))
    ks = [arb(float(x)) for x in s]       # exact dyadic knots
    for k0, k1 in zip(ks, ks[1:]):
        if not bool(k0 <= k1):
            raise RuntimeError("knots not ordered")
    phi = [arb(0)]
    for k in range(1, len(ks)):
        phi.append(phi[-1] + signs[k - 1] * slope * (ks[k] - ks[k - 1]))
    return s, ks, phi, slope, signs


def _phi_at(b: arb, f):
    s, ks, phi, slope, signs = f
    bm = float(b.mid())
    j = int(np.searchsorted(s, bm))
    j = min(max(j, 0), len(s) - 1)
    if j > 0 and abs(s[j - 1] - bm) < abs(s[j] - bm):
        j -= 1
    out = None
    for k in (j - 1, j):
        if k < 0:
            val = phi[0]
        elif k >= len(ks) - 1:
            val = phi[-1]
        else:
            val = phi[k] + signs[k] * slope * (b - ks[k])
        out = val if out is None else out.union(val)
    if not ((j - 1 < 0 or bool(b >= ks[j - 1])) and (j + 1 > len(ks) - 1 or bool(b <= ks[j + 1]))):
        raise RuntimeError("point not localised")
    return out


def _Lt_arb(z: arb, f):
    tot = arb(0)
    for w in WORDS:
        b, wt = z, arb(1)
        for i in w:
            if i == 0:
                wt *= (b + arb(1) / 2) / (2 * b * (arb(7) / 2 - b)); b = (arb(15) / 4 - b).sqrt()
            else:
                wt *= (b - arb(1) / 2) / (2 * b * (b + arb(7) / 2)); b = (arb(15) / 4 + b).sqrt()
        tot += wt * _phi_at(b, f).exp()
    return tot


def _l_arb(x, y, f):
    return (arb(A_NUM) * abs(x - y)).exp() * _Lt_arb(y, f) - _Lt_arb(x, f)


def floor_certificate():
    ctx.prec = PREC
    t0 = time.time()
    s0a, s1a, dn = arb(5).sqrt() / 2, arb(5) / 2, arb(D_NEAR)
    P1f, P2f = (S0 + D_NEAR, S0), (S1, S1 - D_NEAR)           # x > y: l_(x,y) -> d m^-_y as d -> 0
    P1a, P2a = (s0a + dn, s0a), (s1a, s1a - dn)
    r12, sg_s, sg = rho_float(P1f, P2f)                       # g maximises l_P1/l_P2
    r21, sf_s, sf = rho_float(P2f, P1f)                       # f maximises l_P2/l_P1
    print(f"   DL.3: float programme rho(P1,P2) = {r12:.9f}, rho(P2,P1) = {r21:.9f}, sum {r12 + r21:.9f}", flush=True)
    ctx.prec = PREC
    g, f = _make_f(sg_s, sg), _make_f(sf_s, sf)
    vals = {}
    for nm, h in (("u", f), ("v", g)):
        for k, P in (("1", P1a), ("2", P2a)):
            vals[nm, k] = _l_arb(*P, h)
    check(all(bool(v > 0) for v in vals.values()),
          "DL.3: l_P(u), l_P(v) > 0 at both pairs (u = Lt^5 f, v = Lt^5 g; arb; f, g in C_2: slopes +-2(1 - 1e-6))")
    low = (vals["v", "1"] / vals["u", "1"]).log() - (vals["v", "2"] / vals["u", "2"]).log()
    check(bool(low > arb(DL_LOW)),
          f"DL.3: Theta_2(Lt^5 f, Lt^5 g) >= log(Q_P1/Q_P2) = {low.str(10)} > {DL_LOW} (DJ.1's characterisation; "
          f"P1 = (sqrt5/2 + d, sqrt5/2), P2 = (5/2, 5/2 - d), d = {D_NEAR:g}) ({time.time() - t0:.1f}s)")
    RESULT.update(low=low, r12=r12, r21=r21)


# ----------------------------------------------------------------- observation and consequences
def corner_observation():
    out = {}
    for s1, s2 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        z = np.zeros(1)
        out[(s1, s2)] = float(box_bound(np.array([S0]), z, s1, np.array([S1]), z, s2, 0.0, 0.0)[0])
    return out


def consequences():
    ctx.prec = PREC
    Delta = arb(DL_MAX_S)
    theta = arb((Delta / 4).tanh().upper())
    rate = arb((theta.log() / N).exp().upper())
    low = RESULT["low"]
    rate_low = arb(((arb(low.lower()) / 4).tanh().log() / N).exp().lower())
    check(bool(theta < arb(TH_MAX)) and bool(rate < arb(RATE_MAX)) and bool(rate < arb(CAP_J)) and bool(rate > arb(RLOW)),
          f"DL.5: theta_L = tanh(Delta_L/4) <= {theta.str(8)} < {TH_MAX}, per-level rate theta_L^(1/5) <= {rate.str(8)} < "
          f"{RATE_MAX} (DK.4: {RATE_K}; below DJ.3's cap {CAP_J} for (C_100, Lt^9))")
    check(bool(rate_low > arb(RLOW)),
          f"DL.5 (from DL.3): every CF.1 factor for (C_2, Lt^5) is >= tanh(Delta_low/4), so every such per-level rate is "
          f">= {rate_low.str(8)} > {RLOW}")
    import gasket_fractional_profile_level_ratio as ce
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    rh = arb(r_hi.numerator) / r_hi.denominator
    mu = rate * rh
    check(bool(mu < arb(MU_MAX)), f"DL.5: every other C^1 eigenvalue |mu| <= theta_L^(1/5) r <= {mu.str(8)} < {MU_MAX} "
                                  f"(DK.4: {MU_K}; CF.5's proof on C_2)")
    RESULT.update(theta=theta, rate=rate, rate_low=rate_low, mu=mu)


def observation():
    """DK.1's relaxed corner value ((sqrt5/2, -), (5/2, -)) for (C_2, Lt^n), n = 1..9 (floats, no claim)."""
    global N, WD, _PERM
    keep = (N, WD, _PERM)
    out = []
    for n in range(1, 10):
        N = n; WD = np.array(list(itertools.product((0, 1), repeat=n))); _PERM = np.argsort(words([2.0])[0][0])
        z = np.zeros(1)
        D = float(box_bound(np.array([S0]), z, -1, np.array([S1]), z, -1, 0.0, 0.0)[0])
        out.append((n, D, math.tanh(D / 4) ** (1 / n)))
    N, WD, _PERM = keep
    print("   Observation (floats): relaxed corner value and tanh(D/4)^(1/n) for (C_2, Lt^n): "
          + ", ".join(f"n={n}: {D:.6g}, {r:.4f}" for n, D, r in out), flush=True)
    RESULT.update(obs_n=out)


def main() -> int:
    t0 = time.time()
    print("Theorem DL: the small cone C_2 with block Lt^5 (DK's kept limit: (C_100, Lt^9) exhausted). Not thrust.", flush=True)
    admissible()
    sb, G2 = constants()
    float_vs_arb()
    floor_certificate()
    obs = corner_observation()
    print("   observation (floats): relaxed corner sums at ((sqrt5/2, s1), (5/2, s2)): "
          + ", ".join(f"{k}: {v:.8f}" for k, v in obs.items()), flush=True)
    Dmax, stats, nbox, tc = cover(sb, G2)
    ctx.prec = PREC
    Dk = arb(Dmax)
    check(Dmax <= DL_MAX and bool(Dk <= arb(DL_MAX_S)),
          f"DL.4: {nbox} boxes over the sign pairs (+,+), (-,-) (t1 <= t2), (+,-); every accepted box bound <= {Dmax:.8f}, so "
          f"diam_(Theta_2) Lt^5 C_2 <= {Dmax:.8f} <= {DL_MAX_S}; depths "
          + ", ".join(f"{k}: {v[1]}" for k, v in stats.items()) + f" ({tc:.0f}s)")
    check(bool(RESULT["low"] < Dk),
          f"DL.4: with DL.3's floor, {DL_LOW} < diam Lt^5 C_2 <= {Dmax:.8f} (width {Dmax - float(DL_LOW):.2e})")
    check(obs[(-1, -1)] <= Dmax and obs[(-1, -1)] > float(DL_LOW) and obs[(-1, -1)] == max(obs.values()),
          "DL.4 (consistency): the largest relaxed corner value is at (-,-) and lies between DL.3's floor and the bound")
    consequences()
    if OBS:
        observation()
    RESULT.update(D=Dmax, obs=obs, stats=stats, nbox=nbox)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem DL checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
