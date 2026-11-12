#!/usr/bin/env python3
r"""Theorem DM: DL's kept limits "longer blocks only floats (box pad scales with h, not D, so cover cost grows ~D_n^-2)"
and "A box bound whose pad scales with D rather than with h".  Same cone C_2, block Lt^5 -> Lt^7 (128 words).

Setting (GASKET.md, Theorems CE, CF, DC-DL): s in X_s = [sqrt5/2, 5/2], Lt^n f(t) = sum_w W_w(t) f(b_w(t)), C_a, Theta_a
(CF's definition), a = 2, DK.1's functionals m^s_t(u) = a u(t) + s u'(t), Q = (t, s), DK.1's relaxation
M^-_Q(f) <= m^s_t(Lt^n f) <= M^+_Q(f), M^(+-)_Q(f) = sum_w C^(s,+-)_w(t) f(b_w(t)), C^(s,+-)_w = W_w (a + s lambda_w +- a|b_w'|),
and diam_(Theta_a) Lt^n C_a <= sup_(Q1,Q2) [rho~(Q1,Q2) + rho~(Q2,Q1)], rho~(Q1,Q2) = sup_(f in C_a) log(M^+_Q1 f / M^-_Q2 f).

Theorem DM (checked here; exit 1 on any failure):
  DM.1 [proved] (corner box bound)  On a box t_i in [c_i - h, c_i + h] (fixed signs s_1, s_2) write
       log C^(s,+-)_w(t) = l^+-_w + g^+-_w delta + R, |R| <= G2 h^2/2, x^+-_w = g^+-_w - mu (any common mu), and
       K = a s_b h + G2 h^2/2.  By convexity e^(x delta) <= cosh(xh) + (delta/h) sinh(xh) (chord, = e^(+-xh) at
       delta = +-h) and e^(x delta) >= 1 + x delta (tangent, > 0 for |x| h < 1).  Then
       rho~(Q1,Q2) + rho~(Q2,Q1) <= 4K + log(N1(f, d1)/D1(g, d1)) + log(N2(g, d2)/D2(f, d2)) maximised over f, g,
       where N_i (chords) and D_i (tangents) are affine in d_i and positive (the common factors e^(mu d_i) cancel
       between the two terms).  For fixed f, g each quotient is linear-fractional in d_i, hence monotone, so the sup over
       the box is attained at one of the four corners d = (e1 h, e2 h): the box bound is
       4K + max_(e1, e2 = +-1) [rho^(e)(1,2) + rho^(e)(2,1)], each a DE.1 line programme with coefficients
       e^(l^+ + e x^+ h) (numerators) and e^(l^-) (1 + e x^- h) (denominators).  The only first-order pad left is
       4 a s_b h (point motion); DK.3's 2 h max|g - mu| is replaced by the second-order chord-tangent gap
       x^2 h^2/2.
  DM.2 [certified; arb jets to third order on 128 balls of X_s] for n = 7: s_b = sup |b_w'| < 5.80e-4, G2 < 8.80,
       a + s lambda_w +- a|b_w'| > 1.081; b_w' of constant sign on each ball; |(log Lt^7 f)'| < 0.9196 for f in C_2,
       so Lt^7 C_2 subset C_0.9196 and DK.1 holds verbatim for (a, n) = (2, 7) (DL.1's cone step); float word tables
       agree with 200-bit arb at 5 centres.
  DM.3 [certified; arb] (floor, DJ.2's construction at DK.1's (-,-) corner, P1 = (sqrt5/2 + 1e-6, sqrt5/2),
       P2 = (5/2, 5/2 - 1e-6), bang-bang f, g with slopes +-2(1 - 1e-6)): Theta_2(Lt^7 f, Lt^7 g) > 0.0044106.
  DM.4 [certified; DM.1's corner bound, adaptive dyadic cover, outward DP]  diam_(Theta_2) Lt^7 C_2 <= 0.004411,
       so 0.0044106 < diam <= 0.004411.
  DM.5 [certified; arb] theta_M = tanh(0.004411/4) < 0.0011028, per-level rate theta_M^(1/7) < 0.378005
       (DL.5: 0.404994); every CF.1-form rate on (C_2, Lt^7) is > 0.377999 (DM.3); every other C^1 eigenvalue
       |mu| <= theta_M^(1/7) r < 0.16018 (DL.5: 0.17161).
  Observation (floats, no claim): DK.3's pad vs DM.1's at the (-,-) corner box, and DM_N=9 (see GASKET.md).
Kept limits: see GASKET.md Theorem DM.  No Psi_info, Stage 2 still blocked.  Functions on a graph, not a flux.  Not thrust.
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
import gasket_fractional_profile_level_small_cone as dl  # noqa: E402

FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A, A_NUM = 2.0, 2
N = int(os.environ.get("DM_N", "7"))
# per-n constants: (SB_MAX, G2_MAX, Q_MIN, B_IMG, cover threshold, floor, theta max, rate max, rate floor, mu max)
TABLE = {
    7: (5.80e-4, 8.80, 1.081, "0.9196", "0.004411", "0.0044106", "0.0011028", "0.378005", "0.377999", "0.16018"),
    9: (6.45e-5, 8.80, 1.081, "0.9196", "0.0004449", "0.00044485", "0.00011123", "0.363655", "0.363650", "0.15410"),
}
SB_MAX, G2_MAX, Q_MIN, B_IMG, DM_MAX_S, DM_LOW, TH_MAX, RATE_MAX, RLOW, MU_MAX = TABLE[N]
DM_MAX = float(os.environ.get("DM_MAX", DM_MAX_S))
RATE_L, MU_L = "0.404994", "0.17161"         # DL.5 (as stated)
INF, NINF = np.inf, -np.inf


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


_DL_KEEP = ("N", "WD", "WORDS", "_PERM", "DL_MAX", "MAXDEPTH", "BATCH", "box_bound")


def configure(n: int) -> dict:
    """point DL's word tables, DP, floor programme and cover at the block Lt^n (a = 2); returns DL's old settings"""
    keep = {k: getattr(dl, k) for k in _DL_KEEP}
    dl.N = n
    dl.WD = np.array(list(itertools.product((0, 1), repeat=n)))
    dl.WORDS = [tuple(w) for w in dl.WD]
    dl._PERM = np.argsort(dl.words([2.0])[0][0])
    dl.DL_MAX = DM_MAX
    dl.MAXDEPTH = int(os.environ.get("DM_DEPTH", "24"))
    dl.BATCH = int(os.environ.get("DM_BATCH", "2000" if n <= 7 else "500"))
    dl.box_bound = box_bound
    return keep


# ----------------------------------------------------------------- DM.2: constants (DL.2's jets at block n)
def constants():
    ctx.prec = 64
    t0 = time.time()
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    L = s1 - s0
    sb = G2 = lam = 0.0
    qmin = 1e300
    sgn = True
    for j in range(dl.NBALL):
        T = (s0 + L * j / dl.NBALL).union(s0 + L * (j + 1) / dl.NBALL)
        for (b, b1, b2, b3, lw, L1, L2, L3) in dl.jets(T):
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
    check(sgn, f"DM.2: on each of {dl.NBALL} balls of X_s every b_w' ({2 ** N} words) has constant sign")
    check(sb < SB_MAX and G2 < G2_MAX and qmin > Q_MIN,
          f"DM.2 (arb jets, n = {N}): sup |b_w'| <= {sb:.6e} < {SB_MAX}, sup |(log C^(s,+-)_w)''| <= {G2:.5f} < {G2_MAX}, "
          f"a + s lambda_w +- a|b_w'| >= {qmin:.5f} > {Q_MIN} ({time.time() - t0:.1f}s)")
    ctx.prec = PREC
    bimg = arb(lam) + 2 * arb(sb)
    check(bool(bimg < arb(B_IMG)),
          f"DM.2: |(log Lt^{N} f)'| <= sup|lambda_w| + 2 sup|b_w'| <= {float(bimg.upper()):.6f} < {B_IMG} for f in C_2, so "
          f"Lt^{N} C_2 subset C_{B_IMG} (< 2) and every l_P, m^s_t is positive there; with DL.1's cone step DK.1 holds for (2, {N})")
    ctx.prec = 64
    RESULT.update(sb=sb, G2=G2, qmin=qmin, lam=lam)
    return SB_MAX, G2_MAX


def float_vs_arb():
    ctx.prec = PREC
    eb = el = eg = 0.0
    for tc in (dl.S0, 1.4, 1.8, 2.2, dl.S1):
        J = dl.jets(arb(tc))
        for s in (1, -1):
            p, cp, gp, cm, gm = dl.coef([tc], s)
            for idx in range(2 ** N):
                b, b1, b2, b3, lw, L1, L2, L3 = J[idx]
                k = int(np.nonzero(dl._PERM == idx)[0][0])
                sg = 1 if bool(b1 > 0) else -1
                eb = max(eb, abs(float((b - p[0, k]).mid())))
                for sig, cc, gg in ((1, cp, gp), (-1, cm, gm)):
                    q = A_NUM + s * L1 + sig * A_NUM * sg * b1
                    el = max(el, abs(float((lw + q.log() - cc[0, k]).mid())))
                    eg = max(eg, abs(float((L1 + (s * L2 + sig * A_NUM * sg * b2) / q - gg[0, k]).mid())))
    check(eb < 1e-14 and el < 1e-12 and eg < 1e-11,
          f"DM.2: float tables vs 200-bit arb at 5 centres: |b| err {eb:.1e}, |log C| err {el:.1e}, |g| err {eg:.1e} "
          f"(pads {dl.EPS_PT:.0e} on points, {dl.EPS_LOG:.0e} on log coefficients, which absorb h|g err| + |log C err|)")
    RESULT.update(eb=eb, el=el, eg=eg)


# ----------------------------------------------------------------- DM.1: the corner box bound
def _scaled(up, dn):
    sh = up.max(1, keepdims=True)
    return (np.nextafter(np.exp(up - sh) * (1 + 1e-15), INF) + 1e-300,
            np.nextafter(np.exp(dn - sh) * (1 - 1e-15), NINF))


def box_bound(c1, h1, s1, c2, h2, s2, sb, G2):
    """DM.1: certified upper bound of sup over t1 in [c1 +- h1], t2 in [c2 +- h2] of rho~(Q1,Q2) + rho~(Q2,Q1),
    as 4K + the max over the four corners of the chord/tangent line programmes"""
    sides = []
    for c, h, s in ((c1, h1, s1), (c2, h2, s2)):
        p, lp, gp, lm, gm = dl.coef(c, s)
        gall = np.concatenate([gp, gm], 1)
        mu = (gall.max(1) + gall.min(1)) / 2
        H = h[:, None]
        K = ((A * sb + G2 * h / 2) * h)[:, None] + dl.EPS_LOG
        corners = []
        for e in (1, -1):
            xp, xm = (gp - mu[:, None]) * e, (gm - mu[:, None]) * e
            if not (np.abs(xm) * H < 0.5).all():
                raise RuntimeError("tangent not positive: box too large")
            up = lp + xp * H + K                                  # chord at d = e h: exactly e^(e x h)
            dn = lm + np.nextafter(np.log1p(xm * H), NINF) - K    # tangent 1 + e x h
            corners.append(_scaled(up, dn))
        sides.append((p, corners))
    (p1, o1), (p2, o2) = sides
    best = np.full(len(c1), NINF)
    for (A1, B1) in o1:
        for (A2, B2) in o2:
            r = dl._rho_up(*dl._merge(p1, A1, p2, B2)) + dl._rho_up(*dl._merge(p2, A2, p1, B1))
            best = np.maximum(best, r)
    return np.nextafter(best, INF)


_DK_BOX = dl.box_bound          # DK.3's box bound (pad h|g - mu| per coefficient), observation only


def dk_box_bound(c1, h1, s1, c2, h2, s2, sb, G2):
    return _DK_BOX(c1, h1, s1, c2, h2, s2, sb, G2)


# ----------------------------------------------------------------- DM.3: floor (DJ.2's construction, block n)
def floor_certificate():
    ctx.prec = PREC
    t0 = time.time()
    s0a, s1a, dn = arb(5).sqrt() / 2, arb(5) / 2, arb(dl.D_NEAR)
    P1f, P2f = (dl.S0 + dl.D_NEAR, dl.S0), (dl.S1, dl.S1 - dl.D_NEAR)
    P1a, P2a = (s0a + dn, s0a), (s1a, s1a - dn)
    r12, sg_s, sg = dl.rho_float(P1f, P2f)
    r21, sf_s, sf = dl.rho_float(P2f, P1f)
    print(f"   DM.3: float programme rho(P1,P2) = {r12:.10f}, rho(P2,P1) = {r21:.10f}, sum {r12 + r21:.10f}", flush=True)
    ctx.prec = PREC
    g, f = dl._make_f(sg_s, sg), dl._make_f(sf_s, sf)
    vals = {}
    for nm, h in (("u", f), ("v", g)):
        for k, P in (("1", P1a), ("2", P2a)):
            vals[nm, k] = dl._l_arb(*P, h)
    check(all(bool(v > 0) for v in vals.values()),
          f"DM.3: l_P(u), l_P(v) > 0 at both pairs (u = Lt^{N} f, v = Lt^{N} g; arb; f, g in C_2: slopes +-2(1 - 1e-6))")
    low = (vals["v", "1"] / vals["u", "1"]).log() - (vals["v", "2"] / vals["u", "2"]).log()
    check(bool(low > arb(DM_LOW)),
          f"DM.3: Theta_2(Lt^{N} f, Lt^{N} g) >= log(Q_P1/Q_P2) = {low.str(10)} > {DM_LOW} (DJ.1's characterisation; "
          f"P1 = (sqrt5/2 + d, sqrt5/2), P2 = (5/2, 5/2 - d), d = {dl.D_NEAR:g}) ({time.time() - t0:.1f}s)")
    RESULT.update(low=low, r12=r12, r21=r21)


# ----------------------------------------------------------------- observation and consequences
def corner_observation():
    out = {}
    z = np.zeros(1)
    for s1, s2 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        out[(s1, s2)] = float(box_bound(np.array([dl.S0]), z, s1, np.array([dl.S1]), z, s2, 0.0, 0.0)[0])
    return out


def pad_observation(sb, G2):
    """excess of the box bound over the corner value at the (-,-) corner box, DK.3's pad vs DM.1's (floats)"""
    out = []
    z = np.zeros(1)
    D0 = float(box_bound(np.array([dl.S0]), z, -1, np.array([dl.S1]), z, -1, 0.0, 0.0)[0])
    for h in (1e-2, 1e-3, 1e-4):
        hh = np.array([h])
        c1, c2 = np.array([dl.S0 + h]), np.array([dl.S1 - h])
        e_dk = float(dk_box_bound(c1, hh, -1, c2, hh, -1, sb, G2)[0]) - D0
        e_dm = float(box_bound(c1, hh, -1, c2, hh, -1, sb, G2)[0]) - D0
        out.append((h, e_dk, e_dm))
    print("   Observation (floats): box bound minus corner value at the (-,-) corner box of half-width h, DK.3 pad vs DM.1: "
          + ", ".join(f"h={h:g}: {a:.3e} vs {b:.3e}" for h, a, b in out), flush=True)
    RESULT.update(pad_obs=out)
    return out


def consequences():
    ctx.prec = PREC
    Delta = arb(DM_MAX_S)
    theta = arb((Delta / 4).tanh().upper())
    rate = arb((theta.log() / N).exp().upper())
    low = RESULT["low"]
    rate_low = arb(((arb(low.lower()) / 4).tanh().log() / N).exp().lower())
    check(bool(theta < arb(TH_MAX)) and bool(rate < arb(RATE_MAX)) and bool(rate < arb(RATE_L)) and bool(rate > arb(RLOW)),
          f"DM.5: theta_M = tanh(Delta_M/4) <= {theta.str(8)} < {TH_MAX}, per-level rate theta_M^(1/{N}) <= {rate.str(8)} < "
          f"{RATE_MAX} (DL.5: {RATE_L})")
    check(bool(rate_low > arb(RLOW)),
          f"DM.5 (from DM.3): every CF.1 factor for (C_2, Lt^{N}) is >= tanh(Delta_low/4), so every such per-level rate is "
          f">= {rate_low.str(8)} > {RLOW}")
    import gasket_fractional_profile_level_ratio as ce
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    rh = arb(r_hi.numerator) / r_hi.denominator
    mu = rate * rh
    check(bool(mu < arb(MU_MAX)) and bool(mu < arb(MU_L)),
          f"DM.5: every other C^1 eigenvalue |mu| <= theta_M^(1/{N}) r <= {mu.str(8)} < {MU_MAX} (DL.5: {MU_L}; CF.5's proof on C_2)")
    RESULT.update(theta=theta, rate=rate, rate_low=rate_low, mu=mu)


def main() -> int:
    t0 = time.time()
    print(f"Theorem DM: a corner box bound on (C_2, Lt^{N}) (DL's kept limit: pad scales with h, not D). Not thrust.", flush=True)
    keep = configure(N)
    try:
        return _run(t0)
    finally:
        for k, v in keep.items():
            setattr(dl, k, v)


def _run(t0: float) -> int:
    dl.admissible()                       # DL.1's cone step (independent of n)
    sb, G2 = constants()
    float_vs_arb()
    floor_certificate()
    obs = corner_observation()
    print("   observation (floats): relaxed corner sums at ((sqrt5/2, s1), (5/2, s2)): "
          + ", ".join(f"{k}: {v:.9f}" for k, v in obs.items()), flush=True)
    pad_observation(sb, G2)
    Dmax, stats, nbox, tc = dl.cover(sb, G2)
    ctx.prec = PREC
    Dk = arb(Dmax)
    check(Dmax <= DM_MAX and bool(Dk <= arb(DM_MAX_S)),
          f"DM.4: {nbox} boxes over the sign pairs (+,+), (-,-) (t1 <= t2), (+,-); every accepted box bound (DM.1) <= "
          f"{Dmax:.9f}, so diam_(Theta_2) Lt^{N} C_2 <= {Dmax:.9f} <= {DM_MAX_S}; depths "
          + ", ".join(f"{k}: {v[1]}" for k, v in stats.items()) + f" ({tc:.0f}s)")
    check(bool(RESULT["low"] < Dk),
          f"DM.4: with DM.3's floor, {DM_LOW} < diam Lt^{N} C_2 <= {Dmax:.9f} (width {Dmax - float(DM_LOW):.2e})")
    check(obs[(-1, -1)] <= Dmax and obs[(-1, -1)] > float(DM_LOW) and obs[(-1, -1)] == max(obs.values()),
          "DM.4 (consistency): the largest relaxed corner value is at (-,-) and lies between DM.3's floor and the bound")
    consequences()
    RESULT.update(D=Dmax, obs=obs, stats=stats, nbox=nbox, n=N)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem DM checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
