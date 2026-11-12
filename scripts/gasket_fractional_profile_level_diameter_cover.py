#!/usr/bin/env python3
r"""Theorem DK: DJ's kept limit "no certified upper bound by DJ.1 (needs an outward cover of (x1, y1, x2, y2) normalised
as d -> 0)".  Same cone C_100, same block Lt^9 (512 words).  The 4-parameter cover is replaced by an exact reduction to
2 parameters (derivative functionals), and the 2-parameter family is covered outward.

Setting (GASKET.md, Theorems CE, CF, DC-DJ): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
Lt^9 f(t) = sum_w W_w(t) f(b_w(t)) (512 words, b_w = sigma_w), C_a = {f > 0 : log f a-Lipschitz in s}, Theta_a its
Hilbert metric (CF's definition), a = 100.

Theorem DK (checked here; exit 1 on any failure):
  DK.1 [proved] (reduction to derivative functionals) For t in X_s and s = +-1 put m^s_t(u) = a u(t) + s u'(t).  For
       x < y, l_(x,y)(u) = int_x^y e^(a(tau-x)) m^+_tau(u) dtau; for x > y, l_(x,y)(u) = int_y^x e^(a(x-tau)) m^-_tau(u) dtau;
       u(y) = (m^+_y(u) + m^-_y(u))/(2a) wherever u' exists.  So every quotient l_P(v)/l_P(u) of DJ.1 is at most
       ess sup_(tau, s) m^s_tau(v)/m^s_tau(u), and for u = Lt^9 f (f in C_a, u Lipschitz)
           m^s_t(Lt^9 f) = sum_w W_w [(a + s lambda_w) f(b_w) + s b_w' f'(b_w)],  lambda_w = (log W_w)',
       a.e. in t, with |f'| <= a f.  Hence M^-_(t,s)(f) <= m^s_t(Lt^9 f) <= M^+_(t,s)(f), where
           M^+-_(t,s)(f) = sum_w C^(s,+-)_w(t) f(b_w(t)),  C^(s,+-)_w = W_w (a + s lambda_w +- a |b_w'|) > 0,
       and   diam_(Theta_a) Lt^9 C_a <= sup_(Q1, Q2) [rho~(Q1, Q2) + rho~(Q2, Q1)],
             rho~(Q, Q') = sup_(h in C_a) log(M^+_Q(h)/M^-_Q'(h)),   Q = (t, s) in X_s x {+-1}.
       Each rho~ is DE.1's positive line programme (512 + 512 points); the DP over the points in float order with
       gaps d_k >= |x_(k+1) - x_k| (true positions) is a relaxation, so no ordering certificate is needed.  The sum is
       unchanged by scaling M^+_Q and M^-_Q by one common factor, so on a box T1 x T2 (half-widths h1, h2, centres c)
       with mu_i any real,  log C_w(t) - mu_i (t - c) <= log C_w(c) + h |g_w(c) - mu_i| + G2 h^2/2  (g_w = (log C_w)',
       |(log C_w)''| <= G2) and |log h(b_w(t)) - log h(b_w(c))| <= a sb h  (|b_w'| <= sb) give a box bound by one DP at
       the centres with inflated (numerator) / deflated (denominator) coefficients.
  DK.2 [certified; arb jets to third order on 128 balls of X_s] sb < 6.46e-5, G2 < 1.89, every C^(s,+-)_w > 0
       (a + s lambda_w +- a |b_w'| > 99), b_w' of constant sign on each ball; float word tables agree with 200-bit arb
       at 5 sample centres (b, log C to < 1e-12, g to < 1e-10; pads 1e-13 on points, 1e-10 on log coefficients).
  DK.3 [certified; adaptive dyadic box cover of X_s^2 for the sign pairs (+,+), (-,-) (t1 <= t2 by symmetry) and (+,-)
       (covers (-,+) by symmetry), DE.2's float DP with outward rounding after every operation]
           diam_(Theta_100) Lt^9 C_100 <= DK_MAX = 0.014000,
       so with DJ.2: 0.013997 < diam Lt^9 C_100 <= 0.014000 (DH.4: Delta_H < 0.025879).
  DK.4 [certified; arb] theta_K = tanh(Delta_K/4) < 0.0035000, per-level rate theta_K^(1/9) < 0.53349 (DH.4: 0.57117;
       DJ.3's floor for every CF.1-form rate on (C_100, Lt^9): 0.53342), every other C^1 eigenvalue |mu| <= theta_K^(1/9) r
       < 0.22606 (DH.5: 0.24203), and DI.5's level-sum bound with theta_K in place of theta_H.
  Observation (floats, no claim): the relaxed sum at the corner ((5/2, +), (sqrt5/2, +)) is 0.0139981, equal to DJ.2's
       floor to 1e-7; the relaxation (|f'| <= a f taken pointwise) costs at most 4 sb relative per coefficient.
Kept limits: see GASKET.md Theorem DK.  No Psi_info, Stage 2 still blocked.  Functions on a graph, not a flux.  Not thrust.
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
A = 100.0
A_NUM = 100
N = 9
NBALL = 128
S0, S1 = math.sqrt(5) / 2, 2.5
EPS_PT = 1e-13                   # pad on float image points
EPS_LOG = 1e-10                  # pad on float log coefficients
DK_MAX = float(os.environ.get("DK_MAX", "0.014"))      # acceptance threshold of the cover
DK_MAX_S = "0.014000"
TH_MAX, RATE_MAX, MU_MAX = "0.0035000", "0.53349", "0.22606"
DLOW = "0.013997"                # DJ.2 floor (used as stated)
DELTA_H, THETA_H, RATE_H, MU_H = "0.025879", "0.0064696", "0.57117", "0.24203"   # DH.4, DH.5
RLOW = "0.53342"                 # DJ.3
THETA_G = "0.014183"             # DG.3
TI_MAX = "0.064982"              # DI.4's Theta_I'
SB_MAX, G2_MAX = 6.46e-5, 1.89
MAXDEPTH = int(os.environ.get("DK_DEPTH", "22"))
BATCH = int(os.environ.get("DK_BATCH", "3000"))
WD = np.array(list(itertools.product((0, 1), repeat=N)))
INF, NINF = np.inf, -np.inf


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- DK.2: arb jets (third order) on balls of X_s
def jets(T: arb):
    """forward-mode jets of b_w, log W_w at the ball T over all 512 words (prefix tree):
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
    sb = G2 = 0.0
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
            for s in (1, -1):
                for sig in (1, -1):
                    q = A_NUM + s * L1 + sig * A_NUM * sg * b1
                    qp = s * L2 + sig * A_NUM * sg * b2
                    qpp = s * L3 + sig * A_NUM * sg * b3
                    qmin = min(qmin, float(q.lower()))
                    G2 = max(G2, float(abs(L2 + qpp / q - (qp / q) ** 2).upper()))
    check(sgn, f"DK.2: on each of {NBALL} balls of X_s every b_w' (512 words) has constant sign (so |b_w'| is smooth there)")
    check(sb < SB_MAX and G2 < G2_MAX and qmin > 99,
          f"DK.2 (arb jets): sup |b_w'| <= {sb:.6e} < {SB_MAX}, sup |(log C^(s,+-)_w)''| <= {G2:.5f} < {G2_MAX}, "
          f"a + s lambda_w +- a|b_w'| >= {qmin:.5f} > 99 (so every C^(s,+-)_w > 0) ({time.time() - t0:.1f}s)")
    RESULT.update(sb=sb, G2=G2, qmin=qmin)
    return SB_MAX, G2_MAX


# ----------------------------------------------------------------- float word tables
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
          f"DK.2: float tables vs 200-bit arb at 5 centres: |b| err {eb:.1e}, |log C| err {el:.1e}, |g| err {eg:.1e} "
          f"(pads {EPS_PT:.0e} on points, {EPS_LOG:.0e} on log coefficients)")


# ----------------------------------------------------------------- line programme
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
    # Dinkelbach lower start (3 steps), then bisection in log t
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
    for _ in range(26):
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
            ok = Dv <= DK_MAX
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


# ----------------------------------------------------------------- observation
def corner_observation():
    out = {}
    for s1, s2 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        z = np.zeros(1)
        out[(s1, s2)] = float(box_bound(np.array([S1]), z, s1, np.array([S0]), z, s2, 0.0, 0.0)[0])
    return out


def main() -> int:
    t0 = time.time()
    print("Theorem DK: certified upper bound for the Hilbert diameter of Lt^9 C_100 (DJ's kept limit). Not thrust.", flush=True)
    sb, G2 = constants()
    float_vs_arb()
    obs = corner_observation()
    print(f"   observation (floats): relaxed corner sums at ((5/2, s1), (sqrt5/2, s2)): "
          + ", ".join(f"{k}: {v:.8f}" for k, v in obs.items()), flush=True)
    Dmax, stats, nbox, tc = cover(sb, G2)
    ctx.prec = PREC
    Dk = arb(Dmax)
    check(Dmax <= DK_MAX and bool(Dk <= arb(DK_MAX_S)),
          f"DK.3: {nbox} boxes over the sign pairs (+,+), (-,-) (t1 <= t2), (+,-); every accepted box bound <= {Dmax:.8f}, so "
          f"diam_(Theta_100) Lt^9 C_100 <= {Dmax:.8f} <= {DK_MAX_S} (DH.4: Delta_H < {DELTA_H}); depths "
          + ", ".join(f"{k}: {v[1]}" for k, v in stats.items()) + f" ({tc:.0f}s)")
    check(bool(arb(DLOW) < Dk), f"DK.3: with DJ.2's floor, {DLOW} < diam Lt^9 C_100 <= {Dmax:.8f} (interval width "
                                f"{Dmax - float(DLOW):.2e}; DJ: 0.011881)")
    check(obs[(1, 1)] <= Dmax and obs[(1, 1)] > float(DLOW), "DK.3 (consistency): the relaxed corner value lies between DJ.2's "
                                                         "floor and the certified bound")
    consequences()
    RESULT.update(D=Dmax, obs=obs, stats=stats, nbox=nbox)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem DK checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


def consequences():
    ctx.prec = PREC
    Delta = arb(DK_MAX_S)
    theta = arb((Delta / 4).tanh().upper())
    rate = arb((theta.log() / N).exp().upper())
    check(bool(theta < arb(TH_MAX)) and bool(rate < arb(RATE_MAX)) and bool(rate < arb(RATE_H)) and bool(rate > arb(RLOW)),
          f"DK.4: theta_K = tanh(Delta_K/4) <= {theta.str(8)} < {TH_MAX}, per-level rate theta_K^(1/9) <= {rate.str(8)} < {RATE_MAX} "
          f"(DH.4: {RATE_H}; DJ.3's floor for every CF.1-form rate: {RLOW})")
    import gasket_fractional_profile_level_ratio as ce
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    rh = arb(r_hi.numerator) / r_hi.denominator
    mu = rate * rh
    check(bool(mu < arb(MU_MAX)), f"DK.4: every other C^1 eigenvalue |mu| <= theta_K^(1/9) r <= {mu.str(8)} < {MU_MAX} (DH.5: {MU_H})")
    import gasket_fractional_profile_pure_point as cc
    import gasket_fractional_profile_level_birkhoff as bk
    found, sums, counts, okb = cc.period_atoms(bk.KMAX, verbose=False)
    ctx.prec = PREC
    rl = arb(r_lo.numerator) / r_lo.denominator
    r = (rl + rh) / 2
    TI, thG, thH = arb(TI_MAX), arb(THETA_G), arb(THETA_H)
    okd, worst, better = okb, arb(0), []
    for k in range(4, bk.KMAX):
        dev = abs((sums[k + 1] / (r * sums[k])).log())
        bK = TI * arb(min((thG ** ((k - 4) // 8)).upper(), (theta ** ((k - 4) // 9)).upper()))
        okd = okd and bool(dev < bK)
        worst = arb(max(worst.upper(), (dev / bK).upper()))
    for k in (13, 22, 49, 103):
        old = min((thG ** ((k - 4) // 8)).upper(), (thH ** ((k - 4) // 9)).upper())
        new = min((thG ** ((k - 4) // 8)).upper(), (theta ** ((k - 4) // 9)).upper())
        better.append((k, float(old) / float(new)))
    check(okd, f"DK.4: |log(s_(k+1)/(r s_k))| <= Theta_I' min(theta_G^floor((k-4)/8), theta_K^floor((k-4)/9)) for k >= 4 "
               f"(DI.5 with theta_K for theta_H); improvement factor over DI.5 "
               + ", ".join(f"{fa:.3g} at k = {k}" for k, fa in better)
               + f"; consistency: CC's s_4..s_{bk.KMAX} below the bound (max ratio {worst.str(4)})")
    RESULT.update(theta=theta, rate=rate, mu=mu, worst=worst, better=better)


if __name__ == "__main__":
    sys.exit(main())
