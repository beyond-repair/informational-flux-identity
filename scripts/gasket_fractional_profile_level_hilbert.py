#!/usr/bin/env python3
r"""Theorem DG: the cone-lemma term of DF is replaced by a signed local-slope bound.  For u, v in Lt^8 C_100 the Hilbert
distance of C_100 satisfies Theta_100(u, v) <= 2 log K + log(sup(v/u)/inf(v/u)), K = max((a+D-)/(a-D+), (a+D+)/(a-D-)),
where D+ = sup (log Lt^8 f)' and D- = sup -(log Lt^8 f)' over f in C_100 (signed, pointwise in s).  Certified:
D+ < 0.76155, D- < -0.15311 (so every log Lt^8 f is strictly increasing in s), cone-lemma term 0.0738370 -> 0.012225,
Delta_G < 0.056735, theta_G < 0.014183, per-level Birkhoff rate theta_G^(1/8) < 0.58745 (DF: 0.64398).  DF's open item
"the true Hilbert diameter (the cone-lemma term)".

Setting (GASKET.md, Theorems CE, CF, DC-DF): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
sigma_1 s = sqrt(15/4 - s), sigma_2 s = sqrt(15/4 + s), wt1 = (s + 1/2)/(2 s (7/2 - s)), wt2 = (s - 1/2)/(2 s (s + 7/2)),
C_a = {f > 0 : log f a-Lipschitz in s}, Theta_a its Hilbert metric (CF's definition).  a = 100, block Lt^8 (256 words).

Theorem DG (checked here; exit 1 on any failure):
  DG.1 [proved] (i) (factor lemma) For u, v in C_a with log u(x) - log u(y) <= alpha (x - y) and log v(x) - log v(y) >=
       -beta (x - y) for all x > y, where |alpha|, |beta| < a, alpha + beta >= 0 and (if beta < 0) alpha <= a/2:
       Phi_v(x,y)/Phi_u(x,y) <= (v(y)/u(y)) (a + beta)/(a - alpha), Phi_u(x,y) = e^(a d) u(y) - u(x), d = x - y > 0;
       the same with x < y and (alpha, beta) -> (D-, D+).  Hence for u, v in Lt^8 C_100:
           Theta_100(u, v) <= log(sup(v/u)/inf(v/u)) + 2 log K.
       (ii) (signed slope) For f in C_a and a.e. s, (log Lt^8 f)'(s) <= sup over positive f_w satisfying the consecutive
       a-Lipschitz constraints at the 256 points sigma_w s of sum_w wt_w (+-(log wt_w)' + a |sigma_w'|) f_w / sum_w wt_w f_w.
  DG.2 [certified; arb cell enclosures (prefix tree, 4096 cells in s), outward float DP as DE.2/DF.2]
       D+ < 0.76155, D- < -0.15311; so Lt^8 C_100 subset C_(0.76155) (DF.1: C_(1.8457144)) and log Lt^8 f is increasing.
  DG.3 [certified; arb] cone-lemma term 2 log K < 0.012225 (DF: 0.0738370); with DF.2's W_F < 0.04451 (used as stated)
       Delta_G < 0.056735, theta_G = tanh(Delta_G/4) < 0.014183, theta_G^(1/8) < 0.58745 (DF: 0.11835, 0.02958, 0.64398).
  DG.4 [proved from DG.3 with CF.1 and DF.4's proofs, (a, n) = (100, 8)]: every other C^1 eigenvalue |mu| < 0.24893
       (DF.4: 0.27288); |log(s_(k+1)/(r s_k))| <= theta_G^floor((k-4)/8) Theta_G' (Theta_G' < 0.11477, DF.4) for k >= 4.
Kept limits: decoupled bound (the diagonal ratio and the two factors are bounded separately), not the true Hilbert
diameter; no floor for D+-; same cone family and (a, n); Theta_G' (f_0 term) not re-done with DG.1; rate 0.587 still far
from the spectral 0.1725; DP certificate is outward-rounded floats.  No Psi_info, Stage 2 still blocked.  Functions on a
graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_block as df            # noqa: E402

de = df.de
dd = df.dd
cf = dd.cf
FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A = 100.0
A_NUM = 100
N = 8
CELLS = int(os.environ.get("DG_CELLS", "4096"))
W_F = "0.04451"                 # DF.2, used as stated
TG_PRIME = "0.11477"            # DF.4: Theta_100(f_0, Lt f_0) <= Theta_G' < 0.11477 (re-checked below)
THETA_F, RATE_F, CONE_F = "0.02958", "0.64398", "0.0738370"
DP_MAX, DM_MAX, CT_MAX, DELTA_MAX, THETA_MAX, RATE_MAX, MU_MAX = ("0.76155", "-0.15311", "0.012225", "0.056735",
                                                                  "0.014183", "0.58745", "0.24893")
INF, NINF = math.inf, -math.inf


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _up(x: arb) -> float:
    return math.nextafter(float(x.upper()), INF)


def _lo(x: arb) -> float:
    return math.nextafter(float(x.lower()), NINF)


# ----------------------------------------------------------------- DG.2 cell coefficients (arb, prefix tree)
def cell_coeffs(m: int, j: int):
    """For the cell J = [e_j, e_(j+1)] of the m-grid of X_s and every word w (length N, first letter applied first):
    WL <= inf_J wt_w, WH >= sup_J wt_w, NP >= sup_J wt_w ((log wt_w)' + a|sigma_w'|), NM >= sup_J wt_w (a|sigma_w'| -
    (log wt_w)'), and the hull [IL, IH] of sigma_w(J) (sigma_w monotone), sorted by position.  Chain rule along the
    orbit: (log wt_w)' = sum_k l_(i_k)(s_(k-1)) J_(k-1), J_k = prod sigma_(i_j)'(s_(j-1)), with l_1 = 1/(s+1/2) - 1/s +
    1/(7/2-s), l_2 = 1/(s-1/2) - 1/s - 1/(s+7/2), sigma_1' = -1/(2 sigma_1), sigma_2' = 1/(2 sigma_2)."""
    ctx.prec = 100
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    q, h, sv = arb(15) / 4, arb(1) / 2, arb(7) / 2
    e0j, e1j = s0 + (s1 - s0) * j / m, s0 + (s1 - s0) * (j + 1) / m
    layer = [(arb.union(e0j, e1j), e0j, e1j, arb(0), arb(0), arb(1))]
    for _ in range(N):
        nxt = []
        for b, e0, e1, lw, dl, J in layer:
            t = (q - b).sqrt()
            nxt.append((t, (q - e0).sqrt(), (q - e1).sqrt(), lw + ((b + h) / (2 * b * (sv - b))).log(),
                        dl + (1 / (b + h) - 1 / b + 1 / (sv - b)) * J, -J / (2 * t)))
            t = (q + b).sqrt()
            nxt.append((t, (q + e0).sqrt(), (q + e1).sqrt(), lw + ((b - h) / (2 * b * (b + sv))).log(),
                        dl + (1 / (b - h) - 1 / b - 1 / (b + sv)) * J, J / (2 * t)))
        layer = nxt
    K = len(layer)
    a = arb(A_NUM)
    WL, WH, NP, NM, IH, IL = (np.empty(K) for _ in range(6))
    for k, (b, e0, e1, lw, dl, J) in enumerate(layer):
        w = lw.exp()
        aJ = a * abs(J)
        WL[k], WH[k] = _lo(w), _up(w)
        NP[k], NM[k] = _up(w * (dl + aJ)), _up(w * (aJ - dl))
        IH[k], IL[k] = max(_up(e0), _up(e1)), min(_lo(e0), _lo(e1))
    o = np.argsort(IH + IL)
    return WL[o], WH[o], NP[o], NM[o], IH[o], IL[o]


def _V0(t, n, wl, wh, Ep, Em):
    """Outward upper bound for sup over the relaxed cone (f_0 = 1, f_(k+1)/f_k in [Em_k, Ep_k]) of sum_k (n_k - t w_k)
    f_k, w_k in [wl_k, wh_k]: backward recursion V_k = c_k + (Ep_k V if V >= 0 else Em_k V), rounded up."""
    if t >= 0:
        c = np.nextafter(n - np.nextafter(t * wl, NINF), INF)
    else:
        c = np.nextafter(n + np.nextafter((-t) * wh, INF), INF)
    V = float(c[-1])
    for k in range(len(c) - 2, -1, -1):
        V = math.nextafter(c[k] + math.nextafter((Ep[k] if V >= 0 else Em[k]) * V, INF), INF)
    return V


def sup_ratio(n, wl, wh, IH, IL):
    """Certified t with sup_f sum n_k f_k / sum w_k f_k <= t (V0(t) <= 0 checked with outward rounding)."""
    if not bool(np.all(IH[:-1] < IL[1:])):
        raise RuntimeError("images not separated")
    D = np.nextafter(IH[1:] - IL[:-1], INF)          # >= |sigma_w s - sigma_w' s| for consecutive words, s in J
    aD = np.nextafter(A * D, INF)
    Ep = np.nextafter(np.exp(aD) * (1 + 1e-12), INF)
    Em = np.nextafter(np.exp(-aD) * (1 - 1e-12), NINF)
    lo, hi = -5.0, 5.0
    if _V0(hi, n, wl, wh, Ep, Em) > 0:
        raise RuntimeError("DP certificate failed at t = 5")
    for _ in range(45):
        t = (lo + hi) / 2
        if _V0(t, n, wl, wh, Ep, Em) > 0:
            lo = t
        else:
            hi = t
    return hi


def slopes(m=None):
    m = m or CELLS
    t0 = time.time()
    Dp = Dm = NINF
    argp = argm = None
    for j in range(m):
        WL, WH, NP, NM, IH, IL = cell_coeffs(m, j)
        p, q = sup_ratio(NP, WL, WH, IH, IL), sup_ratio(NM, WL, WH, IH, IL)
        if p > Dp:
            Dp, argp = p, j
        if q > Dm:
            Dm, argm = q, j
    print(f"   DG.2: {m} cells, D+ <= {Dp:.8f} (cell {argp}), D- <= {Dm:.8f} (cell {argm}) ({time.time() - t0:.1f}s)",
          flush=True)
    return Dp, Dm


def main_bound():
    df.check_separated()
    Dp_f, Dm_f = slopes()
    ctx.prec = PREC
    Dp, Dm = arb(Dp_f), arb(Dm_f)
    a = arb(A_NUM)
    check(bool(Dp < arb(DP_MAX)) and bool(Dm < arb(DM_MAX)),
          f"DG.2: D+ <= {Dp_f:.8f} < {DP_MAX}, D- <= {Dm_f:.8f} < {DM_MAX}: for f in C_100 and x > y, "
          f"{-Dm_f:.5f} (x-y) <= log Lt^8 f(x) - log Lt^8 f(y) <= {Dp_f:.5f} (x-y)")
    bn = dd.b_iterate(a, N)
    check(bool(Dp < bn) and bool(-Dm < bn),
          f"DG.2: Lt^8 C_100 subset C_({DP_MAX}) (DF.1's b_8 = {bn.str(8)})")
    # DG.1 hypotheses for both orientations: (alpha, beta) = (D+, D-) and (D-, D+)
    hyp = (bool(abs(Dp) < a) and bool(abs(Dm) < a) and bool(Dp + Dm > 0)
           and (bool(Dm >= 0) or bool(Dp <= a / 2)) and (bool(Dp >= 0) or bool(Dm <= a / 2)))
    check(hyp, "DG.1: hypotheses of the factor lemma hold (|D+-| < a, D+ + D- > 0, the positive one <= a/2)")
    K = arb(max((((a + Dm) / (a - Dp))).upper(), ((a + Dp) / (a - Dm)).upper()))
    cone = arb((2 * K.log()).upper())
    cone_F = 2 * ((a + bn) / (a - bn)).log()
    check(bool(cone < arb(CT_MAX)), f"DG.3: cone-lemma term 2 log K <= {cone.str(8)} < {CT_MAX} (DF: {cone_F.str(8)})")
    Delta = arb((cone + arb(W_F)).upper())
    theta = arb((Delta / 4).tanh().upper())
    rate = arb((theta.log() / N).exp().upper())
    check(bool(Delta < arb(DELTA_MAX)) and bool(theta < arb(THETA_MAX)) and bool(rate < arb(RATE_MAX)),
          f"DG.3: Delta_G = 2 log K + W_F <= {Delta.str(8)} < {DELTA_MAX}, theta_G <= {theta.str(8)} < {THETA_MAX}, "
          f"theta_G^(1/8) <= {rate.str(8)} < {RATE_MAX} (DF: 0.11835, {THETA_F}, {RATE_F})")
    check(bool(rate < arb(RATE_F)), "DG.3: the per-level rate is strictly below DF's 0.64398")
    RESULT.update(Dp=Dp, Dm=Dm, K=K, cone=cone, Delta=Delta, theta=theta, rate=rate, bn=bn)


# ----------------------------------------------------------------- DG.1(i) numerical sanity of the factor lemma
def lemma_sanity():
    """Not part of the proof: arb check of g(d) = (e^(ad) - e^(-beta d))/(e^(ad) - e^(alpha d)) <= (a+beta)/(a-alpha)
    on a d-grid for both orientations (the proof is F(0) = 0, F' >= 0 in GASKET.md)."""
    ctx.prec = PREC
    a = arb(A_NUM)
    ok = True
    for al, be in ((RESULT["Dp"], RESULT["Dm"]), (RESULT["Dm"], RESULT["Dp"])):
        bound = (a + be) / (a - al)
        for e in range(-12, 2):
            for mant in (1, 3):
                d = arb(mant) * arb(10) ** e
                g = ((a * d).exp() - (-be * d).exp()) / ((a * d).exp() - (al * d).exp())
                ok = ok and bool(g <= bound)
    check(ok, "DG.1 (sanity, arb): the factor bound holds on d = 1e-12..30 in both orientations")


# ----------------------------------------------------------------- DG.4 consequences
def consequences():
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    theta, rate = RESULT["theta"], RESULT["rate"]
    a = arb(A_NUM)
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    mu = rate * dd._arb(r_hi)
    check(bool(mu < arb(MU_MAX)), f"DG.4: |mu| <= theta_G^(1/8) r <= {mu.str(8)} < {MU_MAX} (DF.4: < 0.27288)")
    ok, Th0, w0 = df._consistency_f1(12, theta, r_lo, r_hi, a)
    check(ok, f"DG.4 (consistency): CF.4's f = 1 check with Theta_0(100) = {Th0.str(6)} and theta_G^floor(m/8) for m < 12 "
              f"(max used fraction {w0:.3e})")
    ctx.prec = PREC
    TG = 2 * ((a + arb(df.B0)) / (a - arb(df.B0))).log() + arb(df.LOGR0)
    TG = arb(TG.upper())
    check(bool(TG < arb(TG_PRIME)), f"DG.4: Theta_G' = {TG.str(8)} < {TG_PRIME} (DF.4, re-evaluated)")
    import gasket_fractional_profile_pure_point as cc
    found, sums, counts, okb = cc.period_atoms(dd.KMAX, verbose=False)
    ctx.prec = PREC
    r = (dd._arb(r_lo) + dd._arb(r_hi)) / 2
    okd, worst = okb, arb(0)
    for k in range(4, dd.KMAX):
        dev = abs((sums[k + 1] / (r * sums[k])).log())
        bG = theta ** ((k - 4) // N) * TG
        okd = okd and bool(dev < bG)
        worst = arb(max(worst.upper(), (dev / bG).upper()))
    check(okd, f"DG.4 (consistency): on CC's s_4..s_{dd.KMAX}, |log(s_(k+1)/(r s_k))| < theta_G^floor((k-4)/8) Theta_G' "
               f"(max ratio {worst.str(4)})")
    thF = arb(THETA_F)
    gains = {k: (thF / theta) ** ((k - 4) // N) for k in (12, 20, 28, 52, 104)}
    print("   DG.4: DF.4 bound / DG.4 bound at k = " + ", ".join(f"{k}: {g.str(4)}" for k, g in gains.items()))
    check(bool(theta < thF), "DG.4: theta_G < theta_F, so DG.4's level-sum bound equals DF.4's for k = 4..11 and is "
                             "smaller by (theta_F/theta_G)^floor((k-4)/8) for k >= 12")
    RESULT.update(mu=mu, TG=TG, worst=worst, gains=gains)


def main() -> int:
    t0 = time.time()
    print("Theorem DG: a signed local-slope bound for DF's cone-lemma term; Birkhoff rate 0.64398 -> 0.58745 per level. "
          "Not thrust.", flush=True)
    main_bound()
    lemma_sanity()
    consequences()
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(df.FAILS) + len(cf.FAILS) + len(dd.FAILS) + len(de.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DG checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
