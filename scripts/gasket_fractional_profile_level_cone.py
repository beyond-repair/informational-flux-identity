#!/usr/bin/env python3
r"""Theorem DD: a sharper diameter for CF.2's cone, so the Birkhoff rate improves from theta^(1/5) < 0.8855 to
theta_D^(1/5) < 0.84339 with the SAME cone C_(47/2) and the SAME block Lt^5 (DC's kept limit "a better cone /
a, theta not re-optimised"; DB/DC's open item).

Setting (GASKET.md, Theorems CE, CF, DC): s = sqrt(t + 25/4) in X_s = [sqrt5/2, 5/2], (Lt f) = wt1 f o sigma_1 +
wt2 f o sigma_2, sigma_(1,2)(s) = sqrt(15/4 -+ s), wt1 = (s+1/2)/(2s(7/2-s)), wt2 = (s-1/2)/(2s(s+7/2)),
C_a = {f > 0 : log f is a-Lipschitz in s}, a = 47/2.  CF.2 bounded diam_(C_a) Lt^5 C_a by
Delta = 2 log((a+b_5)/(a-b_5)) + 2 V_5 with V_5 = 1.0426526 from a branch-by-branch maximum recursion
(osc_J log Lt^n f <= max_i (Var_J log wt_i + osc_(sigma_i J) log Lt^(n-1) f)): theta = tanh(Delta/4) = 0.5441333.

Theorem DD (checked here; exit 1 on any failure):
  DD.1 [proved; exact] (pairwise word bound.)  For x, y in X_s and f in C_a,
           Lt^5 f(x) / Lt^5 f(y) <= max_w  wt_w(x) f(sigma_w x) / (wt_w(y) f(sigma_w y))
                                  <= max_w exp(log wt_w(x) - log wt_w(y) + a |sigma_w x - sigma_w y|),
       over the 32 words w of length 5 (wt_w the product weight along w, sigma_w the composed branch); this is the
       mediant inequality sum p_w / sum q_w <= max p_w/q_w.  Hence, for any cover of X_s by cells J_p,
           osc_X log Lt^5 f <= V_D := max_(p,q) max_w [sup_(J_p) log wt_w - inf_(J_q) log wt_w
                                                        + a max(sup sigma_w J_p - inf sigma_w J_q, sup sigma_w J_q - inf sigma_w J_p)].
  DD.2 [certified; arb, 2048 cells, outward rounding] V_D <= 0.73438, so (b_5 = 2.0785933 as in CF.2)
           diam_(C_a) Lt^5 C_a <= Delta_D = 2 log((a+b_5)/(a-b_5)) + 2 V_D <= 1.82349,
           theta_D = tanh(Delta_D/4) <= 0.42672,  theta_D^(1/5) <= 0.84339  (CF.2: 0.5441333, 0.8854032).
  DD.3 [proved from DD.2 with CF.1, CF.3, CF.5, DC.1-DC.3 as stated] every statement of CF.3, CF.5 and DC holds with
       theta replaced by theta_D: for f in C_b, b < a, |log(Lt^(m+1) f/Lt^m f) - log r| <= theta_D^floor(m/5) Theta_a(f, Lt f);
       every eigenvalue mu != r with a C^1 eigenfunction has |mu| <= theta_D^(1/5) r < 0.3575; DC.1's
       |log R_m - log r| <= theta_D^floor((m-j)/5) K(beta_j); and for the level sums, for every k >= 4,
           |log(s_(k+1)/(r s_k))| <= theta_D^floor((k-4)/5) Theta_G <= C_D theta_D^(k/5),
       Theta_G = 0.28154093 (DC.2), C_D = Theta_G theta_D^(-8/5) <= 1.09984 (the block factor grows; the bound itself is
       smaller than DC.3's for every k >= 9 and equal for k = 4..8).  Consistency (certified): CF.4's check with f = 1
       for m < 12 and CC's s_4..s_13 against the new bound.
  DD.4 [certified; arb point values] floor of the oscillation term: f_+- = exp(+-a s) lie in C_a, and
           max_+- osc_X log Lt^5 f_+- >= V_low = 0.513834 (certified at 2049 points),
       so no bound on sup_(C_a) osc_X log Lt^5 f can go below V_low; V_D/V_low <= 1.4293 (and Delta >= 2 V_low here).
Observation (floats, no claim): over a in [3, 40], n in [3, 8] the same pairwise bound gives its best per-level rate
near a = 47/2, n = 5 (about 0.847 at 128 cells); for large n the oscillation term stalls near 0.594, so this cone family cannot
approach the spectral rate |lambda_2|/r = 0.1725 (CG).
Kept limits: same cone family (log-Lipschitz in s); the bound still takes the maximum over words (mediant
inequality), not the true ratio of sums; C_D's block factor theta_D^(-8/5) is larger than DC's; the Birkhoff rate
0.8434 is still far from 0.1725.  Assumptions used as stated (not re-run): CF.1-CF.3, CF.5's argument, DC.1-DC.3's
proofs and DC.2's Theta_G <= 0.28154093 (DC's 15-min Taylor scan is not re-run).  CF's cone constants, CE.4's
certificate for r and CC's level sums are re-run.  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW enters.
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import itertools
import math
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_rate as cf            # noqa: E402

FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A_CONE = Fr(47, 2)
N_STEP = 5
CELLS = int(os.environ.get("DD_CELLS", "2048"))
THETA_G = "0.28154093"     # DC.2 (assumption used as stated)
KMAX = 13
V_MAX, DELTA_MAX, THETA_MAX, RATE_MAX, CD_MAX = "0.73438", "1.82349", "0.42672", "0.84339", "1.09984"


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _fup(x: arb) -> float:
    return math.nextafter(float(x.upper()), math.inf)


def _flo(x: arb) -> float:
    return math.nextafter(float(x.lower()), -math.inf)


def _arb(fr: Fr) -> arb:
    return arb(fr.numerator) / fr.denominator


# ----------------------------------------------------------------- DD.1/DD.2 the pairwise word bound
def word_tables(m=CELLS, n=N_STEP):
    """For every word w (tuple of branch indices, applied first to last) and every cell J_p = [e_p, e_(p+1)]:
    an outward enclosure of log wt_w on J_p (interval evaluation, arb ball containing J_p) and of sigma_w(J_p)
    (sigma_w is monotone, so its image is the hull of the endpoint images)."""
    ctx.prec = 100
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    q, h, sv = arb(15) / 4, arb(1) / 2, arb(7) / 2
    edges = [s0 + (s1 - s0) * j / m for j in range(m + 1)]
    tabs = []
    for w in itertools.product((0, 1), repeat=n):
        LH, LL, IH, IL = (np.empty(m) for _ in range(4))
        for j in range(m):
            e0, e1 = edges[j], edges[j + 1]
            ball = arb.union(e0, e1)
            lw = arb(0)
            for i in w:
                if i == 0:
                    lw += ((ball + h) / (2 * ball * (sv - ball))).log()
                    ball, e0, e1 = (q - ball).sqrt(), (q - e0).sqrt(), (q - e1).sqrt()
                else:
                    lw += ((ball - h) / (2 * ball * (ball + sv))).log()
                    ball, e0, e1 = (q + ball).sqrt(), (q + e0).sqrt(), (q + e1).sqrt()
            LH[j], LL[j] = _fup(lw), _flo(lw)
            IH[j], IL[j] = max(_fup(e0), _fup(e1)), min(_flo(e0), _flo(e1))
        tabs.append((w, LH, LL, IH, IL))
    return tabs


def pair_bound(tabs, a: float, chunk=512):
    """max over cell pairs and words; float sums of <= 4 terms of size < 40 are padded by 1e-12 (>> 4 ulp)."""
    best, arg = -1.0, None
    for w, LH, LL, IH, IL in tabs:
        m = len(LH)
        for p0 in range(0, m, chunk):
            sl = slice(p0, min(m, p0 + chunk))
            d = np.maximum(IH[sl, None] - IL[None, :], IH[None, :] - IL[sl, None])
            D = LH[sl, None] - LL[None, :] + a * d
            v = float(D.max())
            if v > best:
                i = np.unravel_index(int(D.argmax()), D.shape)
                best, arg = v, (w, p0 + int(i[0]), int(i[1]))
    return best + 1e-12, arg


def b_iterate(a: arb, n=N_STEP):
    """CF.2's b_n: n-th iterate of b -> max(M_1 + c_1 b, M_2 + c_2 b) from a (same formulas as CF.2)."""
    s0 = arb(5).sqrt() / 2
    M1, M2 = arb(14) / 15, 1 / (s0 - arb(1) / 2) - 1 / s0 - 1 / (s0 + arb(7) / 2)
    c1, c2 = 1 / arb(5).sqrt(), 1 / (2 * (arb(15) / 4 + s0).sqrt())
    b = a
    for _ in range(n):
        b = arb(max((M1 + c1 * b).upper(), (M2 + c2 * b).upper()))
    return b


def sharper_diameter():
    ctx.prec = PREC
    t0 = time.time()
    tabs = word_tables()
    a = _arb(A_CONE)
    Vf, arg = pair_bound(tabs, float(A_CONE))
    ctx.prec = PREC
    V = arb(Vf)
    b5 = b_iterate(a)
    Delta = 2 * ((a + b5) / (a - b5)).log() + 2 * V
    theta = (Delta / 4).tanh()
    rate = (theta.log() / N_STEP).exp()
    w, p, qq = arg
    print(f"   DD.2: {CELLS} cells, 32 words, {time.time() - t0:.1f}s; maximiser word {''.join(str(i + 1) for i in w)}, "
          f"cells {p}, {qq}")
    check(bool(V < arb(V_MAX)), f"DD.2: osc_X log Lt^5 f <= V_D = {V.str(8)} < {V_MAX} for every f in C_(47/2) (CF.2: 1.0426526)")
    check(bool(b5 < arb("2.0785934")) and bool(b5 > arb("2.0785932")), f"DD.2: b_5 = {b5.str(10)} (CF.2's value)")
    check(bool(Delta < arb(DELTA_MAX)) and bool(theta < arb(THETA_MAX)) and bool(rate < arb(RATE_MAX)),
          f"DD.2: Delta_D = {Delta.str(8)} < {DELTA_MAX}, theta_D = tanh(Delta_D/4) = {theta.str(8)} < {THETA_MAX}, "
          f"theta_D^(1/5) = {rate.str(8)} < {RATE_MAX}")
    RESULT.update(V=V, b5=b5, Delta=Delta, theta=theta, rate=rate, tabs=tabs)
    return tabs


# ----------------------------------------------------------------- DD.4 floor of the oscillation term
def _Lt5_exp(sv: arb, sign: int, a: arb, n=N_STEP) -> arb:
    if n == 0:
        return (sign * a * sv).exp()
    q = arb(15) / 4
    w1 = (sv + arb(1) / 2) / (2 * sv * (arb(7) / 2 - sv))
    w2 = (sv - arb(1) / 2) / (2 * sv * (sv + arb(7) / 2))
    return w1 * _Lt5_exp((q - sv).sqrt(), sign, a, n - 1) + w2 * _Lt5_exp((q + sv).sqrt(), sign, a, n - 1)


def oscillation_floor(npts=None):
    ctx.prec = 100
    npts = npts or (CELLS + 1)
    a = _arb(A_CONE)
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    best = arb(0)
    for sign in (1, -1):
        vals = [_Lt5_exp(s0 + (s1 - s0) * j / (npts - 1), sign, a).log() for j in range(npts)]
        hi = max(vals, key=lambda x: float(x.mid()))
        lo = min(vals, key=lambda x: float(x.mid()))
        osc = hi - lo
        if bool(osc.lower() > best.lower()) or best == 0:
            best = arb(osc.lower())
    ctx.prec = PREC
    ratio = RESULT["V"] / best
    check(bool(best > 0) and bool(best < RESULT["V"]),
          f"DD.4: f = exp(+-a s) in C_a gives osc_X log Lt^5 f >= V_low = {best.str(8)} ({npts} points), so "
          f"V_D/V_low <= {ratio.str(5)}; no bound on osc over Lt^5 C_a goes below V_low")
    RESULT.update(Vlow=best, Vratio=ratio)


# ----------------------------------------------------------------- DD.3 consequences and consistency
def consequences():
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    theta, rate = RESULT["theta"], RESULT["rate"]
    r_lo, r_hi = ce.certificate()
    rhi = _arb(r_hi)
    check(bool(rate * rhi < arb("0.3575")), f"DD.3: |mu| <= theta_D^(1/5) r <= {(rate * rhi).str(8)} < 0.3575 (CF.5: < 0.3752)")
    consts = dict(a=_arb(A_CONE), b1=None, theta=theta)
    ok = cf.pointwise_consistency(12, consts, r_lo, r_hi)
    check(ok, "DD.3 (consistency): CF.4's f = 1 check at five points for m < 12 holds with theta_D")
    ThetaG = arb(THETA_G)
    CD = arb(((ThetaG / (theta.log() * 8 / 5).exp()) * (1 + arb("1e-40"))).upper())   # rounded up
    check(bool(CD < arb(CD_MAX)), f"DD.3: C_D = Theta_G theta_D^(-8/5) = {CD.str(8)} < {CD_MAX} (DC.3's C_G^B = 0.74544057 "
                                  f"with theta = 0.5441333; theta_D^floor((k-4)/5) <= theta^floor((k-4)/5) for every k)")
    import gasket_fractional_profile_pure_point as cc
    found, sums, counts, okb = cc.period_atoms(KMAX, verbose=False)
    ctx.prec = PREC
    r = (_arb(r_lo) + _arb(r_hi)) / 2
    worst = arb(0)
    okd = okb
    for k in range(4, KMAX):
        dev = abs((sums[k + 1] / (r * sums[k])).log())
        bnd = theta ** ((k - 4) // 5) * ThetaG
        okd = okd and bool(dev < bnd) and bool(bnd <= CD * (theta.log() * k / 5).exp())
        worst = arb(max(worst.upper(), (dev / bnd).upper()))
    check(okd, f"DD.3 (consistency): on CC's s_4..s_{KMAX}, |log(s_(k+1)/(r s_k))| < theta_D^floor((k-4)/5) Theta_G "
               f"<= C_D theta_D^(k/5) (max ratio deviation/bound {worst.str(4)})")
    # gain over DC.3 at a few k
    th0 = arb("0.5441333187")
    gains = {k: (th0 / theta) ** ((k - 4) // 5) for k in (9, 14, 24, 54)}
    print("   DD.3: DC.3 bound / DD.3 bound at k = " + ", ".join(f"{k}: {g.str(5)}" for k, g in gains.items()))
    RESULT.update(CD=CD, worst=worst, gains=gains)


def observation():
    """floats only: the same pairwise bound over a grid of (a, n); no claim."""
    print("   observation (floats, coarse 128 cells, no claim): best per-level rate of the pairwise bound by n")
    for n in (4, 5, 6):
        tabs = word_tables(128, n)
        row = []
        for a in (12.0, 16.0, 23.5, 30.0, 40.0):
            V, _ = pair_bound(tabs, a)
            b = float(b_iterate(arb(a), n).mid())
            D = 2 * math.log((a + b) / (a - b)) + 2 * V
            row.append((a, math.tanh(D / 4) ** (1 / n)))
        print(f"   n = {n}: " + ", ".join(f"a={a:g}: {r:.4f}" for a, r in row))


def main(obs=True) -> int:
    t0 = time.time()
    print("Theorem DD: a sharper diameter for CF.2's cone C_(47/2) under Lt^5 (pairwise word bound): Birkhoff rate "
          "0.8855 -> 0.8434 per level. Not thrust.", flush=True)
    cf.FAILS.clear()
    base = cf.cone_constants()                 # CF.2 re-run (theta = 0.5441333, for comparison)
    RESULT["theta_CF"] = base["theta"]
    sharper_diameter()
    oscillation_floor()
    consequences()
    if obs:
        observation()
    RESULT.pop("tabs", None)
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(cf.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DD checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
