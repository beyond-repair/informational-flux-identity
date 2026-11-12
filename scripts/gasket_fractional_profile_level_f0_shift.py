#!/usr/bin/env python3
r"""Theorem DO: DN's kept failure DN.4 ("the C_2 constant Theta_2(f_0, Lt f_0) > 2.22 x DI.4's, so for k = 4..8 DK.4/DI.5
stay the bound") and the question "does another cone family beat DI.4 at small k?".

Setting (GASKET.md, Theorems CE, CF, DC-DN): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
C_a = {f > 0 : log f a-Lipschitz in s}, Theta_a its Hilbert metric (CF's definition); f_0 = F0 = Lt^4 G (DC.2,
CO.1/CO.2's Taylor data re-run unchanged), F_m = Lt^m F0, R_m = F_(m+1)/F_m, s_(4+m) = c0 F_m(x0) (DI.5's consistency).

Theorem DO (checked here; exit 1 on any failure):
  DO.1 [certified; arb jets of order 2 (value, first and second derivative) through m + 1 <= 6 branches of Lt on DO_CELLS
       cells] (log F_m)' in [0.03529, 0.79032] for m = 0..6: every F_m lies in C_0.79032, inside C_2 and C_100 with slope < a.
  DO.2 [certified; DI.1's exact identity and DI.4's cell-pair pass (outward floats) on arb jet tables, R and R' by
       second-order centred forms; floors: arb, DI.1's quotient at grid pairs and its limits x -> y+-]
       two-sided bounds for Theta_a(F_m, F_(m+1)), m = 0..5, a = 2 and a = 100 (table in GASKET.md), and the cone-free
       floor osc log R_m (every Hilbert metric of a cone of positive functions dominates it).
  DO.3 [proved from DO.2 with CF.1 and the contraction pairs as stated: DL.5/DM.5 theta_L (n = 5), theta_7, theta_9 on C_2;
       DG.3's theta_G (n = 8), DK.4's theta_K (n = 9) on C_100] for every k >= 4,
         |log(s_(k+1)/(r s_k))| <= B(k) = min over m <= min(5, k-4), a in {2, 100} of Theta_a(F_m, F_(m+1)) G_a(k-4-m),
       G_a(j) = min of products of the block factors whose block lengths sum to at most j (blocks may be mixed).
  DO.4 [kept failure, certified witness] for k = 4 + m (m = 0..5) every CF.3-form bound, in any cone of positive functions
       (any a, a toward b*, eigenfunction-weighted, ...), is >= osc log R_m; at k = 4 that floor is 0.0633864, so no cone
       family beats DI.4's 0.064982 by more than the factor 1.0252, and C_2 is never the best cone for k <= 9.
Kept limits: see GASKET.md Theorem DO.  No Psi_info, Stage 2 still blocked.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time
from math import comb

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_f0_coupled as di          # noqa: E402

bk = di.bk
FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
NC = int(os.environ.get("DO_CELLS", "4096"))
NFP = int(os.environ.get("DO_FLOOR_PTS", "1025"))
NPAIR = int(os.environ.get("DO_PAIR_PTS", "257"))
MMAX = 5
AS = (2, 100)
# contraction pairs used as stated
PAIRS = {2: ((5, "0.0108954"), (7, "0.0011028"), (9, "0.00011123")),     # DL.5, DM.5
         100: ((8, "0.014183"), (9, "0.0035000"))}                       # DG.3, DK.4
TN_DN, TI_DI = "0.14473", "0.064982"                                     # DN.2, DI.4 (DN.3/DK.4's constants)
KJET = 3
# certified targets (set after the first runs; checked, exit 1 if not met).  Default 4096 cells; DO_CELLS=16384 is local.
TARGETS = {
    4096: dict(
        UPPER={2: ("0.144458", "0.069733", "0.0069175", "0.0010024", "0.00016633", "0.000029590"),
               100: ("0.0649756", "0.0375169", "0.0042203", "0.00061746", "0.00010201", "0.000017446")},
        LOWER={2: ("0.144449", "0.069688", "0.0069152", "0.0010008", "0.00016484", "0.000028125"),
               100: ("0.0649755", "0.0375166", "0.0042202", "0.00061741", "0.00010197", "0.000017410")}),
    16384: dict(
        UPPER={2: ("0.1444515", "0.0696993", "0.0069156", "0.0010010", "0.00016495", "0.000028218"),
               100: ("0.06497554", "0.0375168", "0.00422022", "0.00061743", "0.00010198", "0.000017413")},
        LOWER={2: ("0.144449", "0.069688", "0.0069152", "0.0010008", "0.00016484", "0.000028125"),
               100: ("0.0649755", "0.0375166", "0.0042202", "0.00061741", "0.00010197", "0.000017410")}),
}
TG = TARGETS.get(NC, TARGETS[4096])
UPPER, LOWER = TG["UPPER"], TG["LOWER"]
OSC_LOWER = ("0.0633864", "0.0368271", "0.0041515", "0.00060778", "0.00010040", "0.000017143")
SLOPE = ("0.03529", "0.79032")


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- jets (f, f', f''/2) with arb coefficients
def jadd(a, b):
    return [x + y for x, y in zip(a, b)]


def jmul(a, b):
    return [sum((a[i] * b[k - i] for i in range(k + 1)), arb(0)) for k in range(KJET)]


def jrec(a):
    b = [1 / a[0]]
    for k in range(1, KJET):
        b.append(-sum((a[i] * b[k - i] for i in range(1, k + 1)), arb(0)) * b[0])
    return b


def jsqrt(a):
    c0 = a[0].sqrt()
    c1 = a[1] / (2 * c0)
    return [c0, c1, (a[2] - c1 * c1) / (2 * c0)]


class TaylorJet:
    """Jets of F0(s) = sum f_n u^n, u = (s - c)/R, from DC.2/CO.2's coefficients: F0^(j)(u)/j! = sum_n f_n C(n,j) u^(n-j)
    plus the Cauchy tail sum_(n>N) C(n,j) S rho^-n |u|^(n-j) = S rho^-j q^(M-j) sum_(l<=j) C(M,j-l) q^l/(1-q)^(l+1)
    (q = |u|/rho, M = N+1; Vandermonde with sum_i C(i,l) q^i = q^l/(1-q)^(l+1)), composed with the input jet."""

    def __init__(self, T):
        self.f, self.S, self.rho, self.c, self.R = T.f, T.S, T.rho, T.c, T.R
        self.N = len(T.f) - 1
        self.d = [[T.f[n] * comb(n, j) for n in range(j, self.N + 1)] for j in range(KJET)]

    def jet(self, X):
        U0 = (X[0] - self.c) / self.R
        U1, U2 = X[1] / self.R, X[2] / self.R
        q = arb(abs(U0).upper()) / self.rho
        assert bool(q < 1)
        M = self.N + 1
        c = []
        for j in range(KJET):
            acc = arb(0)
            for t in reversed(self.d[j]):
                acc = acc * U0 + t
            tail = self.S * self.rho ** (-j) * q ** (M - j) * sum(
                (comb(M, j - l) * q ** l / (1 - q) ** (l + 1) for l in range(j + 1)), arb(0))
            c.append(acc + arb(0, tail.upper()))
        return [c[0], c[1] * U1, c[1] * U2 + c[2] * U1 * U1]


def Fjet(TJ, m, X):
    """Jet of F_m = Lt^m F0 at the jet X: F_m(X) = wt1(X) F_(m-1)(sigma_1(X)) + wt2(X) F_(m-1)(sigma_2(X)),
    wt1 = (s+1/2)/(2s(7/2-s)), wt2 = (s-1/2)/(2s(s+7/2)), sigma_1 = sqrt(15/4 - s), sigma_2 = sqrt(15/4 + s)."""
    if m == 0:
        return TJ.jet(X)
    q4, half, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    rX = [half * x for x in jrec(X)]
    w1 = jmul(jmul([X[0] + half, X[1], X[2]], rX), jrec([s7 - X[0], -X[1], -X[2]]))
    w2 = jmul(jmul([X[0] - half, X[1], X[2]], rX), jrec([s7 + X[0], X[1], X[2]]))
    s1 = jsqrt([q4 - X[0], -X[1], -X[2]])
    s2 = jsqrt([q4 + X[0], X[1], X[2]])
    return jadd(jmul(w1, Fjet(TJ, m - 1, s1)), jmul(w2, Fjet(TJ, m - 1, s2)))


# ----------------------------------------------------------------- DO.1 / DO.2 tables
def _qcentred(wJ, wM, dx):
    """(log w)' on J: direct w'/w intersected with q(c) + q'(J)(J - c), q' = (2 w2 w0 - w1^2)/w0^2 (w2 = w''/2)."""
    qd = wJ[1] / wJ[0]
    qp = (2 * wJ[2] * wJ[0] - wJ[1] * wJ[1]) / (wJ[0] * wJ[0])
    qc = wM[1] / wM[0] + qp * dx
    return bk._hull(arb(max(qd.lower(), qc.lower())), arb(min(qd.upper(), qc.upper())))


def cell_tables(TJ, m):
    """Per cell J: enclosures of R_m, q_m = (log F_m)', q_(m+1), R_m'.  R(J) in R(c) + R'(c)(J-c) + [R''/2](J)(J-c)^2 and
    R'(J) in R'(c) + R''(J)(J-c) (Taylor with Lagrange remainder, c the midpoint), intersected with the direct jets."""
    ctx.prec = PREC
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    s0l, s1u = arb(s0.lower()), arb(s1.upper())
    h = (s1u - s0l) / NC
    tab = np.zeros((NC, 8))
    one, zero = arb(1), arb(0)
    for i in range(NC):
        J = bk._hull(s0l + h * i, s0l + h * (i + 1))
        mp = arb(J.mid())
        uJ, vJ = Fjet(TJ, m, [J, one, zero]), Fjet(TJ, m + 1, [J, one, zero])
        uM, vM = Fjet(TJ, m, [mp, one, zero]), Fjet(TJ, m + 1, [mp, one, zero])
        if not (bool(uJ[0] > 0) and bool(vJ[0] > 0)):
            check(False, f"DO.1: positivity not certified on cell {i} (m = {m})")
            return None
        RJ, RM = jmul(vJ, jrec(uJ)), jmul(vM, jrec(uM))
        dx = J - mp
        Rc = RM[0] + RM[1] * dx + RJ[2] * dx * dx
        R = bk._hull(arb(max(RJ[0].lower(), Rc.lower())), arb(min(RJ[0].upper(), Rc.upper())))
        Rpc = RM[1] + 2 * RJ[2] * dx
        Rp = bk._hull(arb(max(RJ[1].lower(), Rpc.lower())), arb(min(RJ[1].upper(), Rpc.upper())))
        q0, q1 = _qcentred(uJ, uM, dx), _qcentred(vJ, vM, dx)
        tab[i] = (di._fl(R), di._fu(R), di._fl(q0), di._fu(q0), di._fl(q1), di._fu(q1), di._fl(Rp), di._fu(Rp))
    return tab, h


def floors(TJ, m):
    """Certified lower bounds (arb): for a in AS, Theta_a(F_m, F_(m+1)) >= log(sup/inf) of DI.1's quotient values at
    NPAIR x NPAIR grid pairs, its limits x -> y+- (R(y)(a -+ q1)/(a -+ q0)) at NFP points (endpoints included,
    by continuity of the limit functions), and R itself; and
    osc log R_m >= log(max R/min R) over the NFP points."""
    ctx.prec = PREC
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    one, zero = arb(1), arb(0)
    pts = []
    for k in range(NFP):
        s = s0 if k == 0 else (s1 if k == NFP - 1 else arb((s0 + (s1 - s0) * k / (NFP - 1)).mid()))
        u, v = Fjet(TJ, m, [s, one, zero]), Fjet(TJ, m + 1, [s, one, zero])
        pts.append((s, u[0], v[0] / u[0], u[1] / u[0], v[1] / v[0]))
    out = {}
    Rhi_lo = max((p[2].lower() for p in pts))
    Rlo_hi = min((p[2].upper() for p in pts))
    osc = arb((arb(Rhi_lo) / arb(Rlo_hi)).log().lower())
    step = (NFP - 1) // (NPAIR - 1)
    grid = pts[::step]
    for a_num in AS:
        a = arb(a_num)
        hi, lo = Rhi_lo, Rlo_hi
        for k, (s, u, R, q0, q1) in enumerate(pts):
            # L+(y) = lim_(x->y+) and L-(y) = lim_(x->y-) of DI.1's quotient; both are continuous on X_s (F_m analytic
            # near X_s), so their values at the endpoints are limits of values at interior points, hence in [alpha, beta]
            cs = [R * (a - q1) / (a - q0), R * (a + q1) / (a + q0)]
            for c in cs:
                hi, lo = max(hi, c.lower()), min(lo, c.upper())
        for (x, ux, Rx, _, _) in grid:
            for (y, uy, Ry, _, _) in grid:
                if x is y:
                    continue
                Q = Ry + (Ry - Rx) / ((a * abs(x - y)).exp() * uy / ux - 1)
                hi, lo = max(hi, Q.lower()), min(lo, Q.upper())
        out[a_num] = arb((arb(hi) / arb(lo)).log().lower())
    return out, osc


def coupled_mediant(tab, h, a_num):
    """DO.2's cell-pair pass (outward floats): DI.4 with the mean value theorem replaced by integrals.  For y in cell j and
    x in cell l, write [y, x] (or [x, y]) as a partial piece of cell j, the full cells strictly between (``mid''), and a
    partial piece of cell l.  With nu_i >= the relevant signed part of R' on cell i and kappa_i <= a -+ q_0 on cell i,
    |R(y) - R(x)| <= N = int nu and e^(a d) u(y)/u(x) - 1 >= e^D - 1, D = int kappa, so
        |E| <= N/(e^D - 1) = (N/D) D/(e^D - 1) <= max_corners(N/D) phi(D_mid), phi(t) = t/(e^t - 1) decreasing (D >= D_mid):
    with lengths alpha, beta in [0, h] of the partial pieces, N/D = (N_mid + nu_j alpha + nu_l beta)/(D_mid + kappa_j alpha +
    kappa_l beta) is linear-fractional, hence monotone, in each of alpha, beta, so its sup is at one of the 4 corners (when
    D_mid = 0 the mediant inequality gives max(nu_j/kappa_j, nu_l/kappa_l)); and also
    |E| <= (N_mid + (nu_j + nu_l) h)/(e^(D_mid) - 1) (partial pieces have length <= h); the smaller is used per pair.  For l = j (same cell) and l = j +- 1 the
    middle is empty and the bound is max(nu_j/kappa_j, nu_l/kappa_l).  Returns (beta upper, alpha lower)."""
    Rlo, Rhi, q0lo, q0hi, Rplo, Rphi = tab[:, 0], tab[:, 1], tab[:, 2], tab[:, 3], tab[:, 6], tab[:, 7]
    n = len(Rlo)
    hf = float(h.lower()) * (1 - 1e-15)
    af = float(a_num)
    pos = np.maximum(Rphi, 0.0) * (1 + 1e-14)       # >= (R')_+ on the cell
    neg = np.maximum(-Rplo, 0.0) * (1 + 1e-14)      # >= (-R')_+ on the cell
    kr = (af - q0hi) * (1 - 1e-14)                  # x > y: kappa = a - q0
    kl = (af + q0lo) * (1 - 1e-14)                  # x < y: kappa = a + q0
    assert kr.min() > 0 and kl.min() > 0
    UP, DN = 1 + 1e-10, 1 - 1e-10                   # outward pads for prefix sums (n <= 2^20 terms, relative error < 1e-10)

    def pref(v):
        return np.concatenate(([0.0], np.cumsum(v * hf)))

    PN, PP, PKr, PKl = pref(neg), pref(pos), pref(kr), pref(kl)
    up, dn = np.zeros(n), np.zeros(n)
    idx = np.arange(n)
    for j in range(n):
        # x to the right (l > j): up uses nu = neg, dn uses nu = pos, kappa = kr
        # x to the left  (l < j): up uses nu = pos, dn uses nu = neg, kappa = kl
        cands_up, cands_dn = [max(neg[j] / kr[j], pos[j] / kl[j])], [max(pos[j] / kr[j], neg[j] / kl[j])]
        for side in (1, -1):
            if side == 1:
                ls = idx[j + 1:]
                if len(ls) == 0:
                    continue
                Dm = (PKr[ls] - PKr[j + 1]) * DN
                Nu = (PN[ls] - PN[j + 1]) * UP
                Nd = (PP[ls] - PP[j + 1]) * UP
                kk, nu_u, nu_d = kr, neg, pos
            else:
                ls = idx[:j][::-1]
                if len(ls) == 0:
                    continue
                Dm = (PKl[j] - PKl[ls + 1]) * DN
                Nu = (PP[j] - PP[ls + 1]) * UP
                Nd = (PN[j] - PN[ls + 1]) * UP
                kk, nu_u, nu_d = kl, pos, neg
            phi = di._mvfac(Dm)
            kj, kl_ = kk[j] * hf * DN, kk[ls] * hf * DN
            eu = np.maximum(nu_u[j] / kk[j], nu_u[ls] / kk[ls])        # mediant (used when D_mid = 0)
            ed = np.maximum(nu_d[j] / kk[j], nu_d[ls] / kk[ls])
            with np.errstate(divide="ignore", invalid="ignore"):
                for (ea, eb) in ((0, 0), (1, 0), (0, 1), (1, 1)):
                    Dc = Dm + ea * kj + eb * kl_
                    cu = (Nu + (ea * nu_u[j] + eb * nu_u[ls]) * hf * UP) / Dc
                    cd = (Nd + (ea * nu_d[j] + eb * nu_d[ls]) * hf * UP) / Dc
                    if ea == 0 and eb == 0:
                        cu, cd = np.where(Dm > 0, cu, 0.0), np.where(Dm > 0, cd, 0.0)
                        mu_, md_ = cu, cd
                    else:
                        mu_, md_ = np.maximum(mu_, cu), np.maximum(md_, cd)
                corner_u = np.where(Dm > 0, mu_ * UP, eu)
                corner_d = np.where(Dm > 0, md_ * UP, ed)
            eu, ed = corner_u * phi, corner_d * phi
            # second bound per pair: N <= N_mid + (nu_j + nu_l) h, e^D - 1 >= e^(D_mid) - 1
            den = np.expm1(Dm) * DN
            with np.errstate(divide="ignore", invalid="ignore"):
                fu = np.where(Dm > 0, (Nu + (nu_u[j] + nu_u[ls]) * hf * UP) / den, np.inf)
                fd = np.where(Dm > 0, (Nd + (nu_d[j] + nu_d[ls]) * hf * UP) / den, np.inf)
            eu, ed = np.minimum(eu, fu), np.minimum(ed, fd)
            cands_up.append(float(eu.max()))
            cands_dn.append(float(ed.max()))
        up[j] = max(cands_up) * (1 + 1e-12)
        dn[j] = max(cands_dn) * (1 + 1e-12)
    pad = 1 + 1e-12
    beta = float(np.max((Rhi + up) * pad))
    alpha = float(np.min((Rlo - dn) / pad))
    return beta, alpha


def main_tables(T):
    TJ = TaylorJet(T)
    res = {}
    slo, shi = math.inf, -math.inf
    for m in range(MMAX + 1):
        t0 = time.time()
        out = cell_tables(TJ, m)
        if out is None:
            return False
        tab, h = out
        t1 = time.time() - t0
        slo = min(slo, float(tab[:, 2].min()), float(tab[:, 4].min()))
        shi = max(shi, float(tab[:, 3].max()), float(tab[:, 5].max()))
        ups = {}
        for a_num in AS:
            beta1, alpha1, *_ = di.coupled(tab, h, arb(a_num))
            beta2, alpha2 = coupled_mediant(tab, h, a_num)
            beta, alpha = min(beta1, beta2), max(alpha1, alpha2)
            ups[a_num] = arb((arb(beta) / arb(alpha)).log().upper())
            print(f"      a = {a_num}: DI.4 pass {math.log(beta1 / alpha1):.9e}, mediant pass {math.log(beta2 / alpha2):.9e}", flush=True)
        fl, osc = floors(TJ, m)
        res[m] = dict(up=ups, low=fl, osc=osc)
        print(f"   m = {m}: Theta_2 in [{fl[2].str(7)}, {ups[2].str(7)}], Theta_100 in [{fl[100].str(7)}, {ups[100].str(7)}], "
              f"osc log R_m >= {osc.str(7)} (cells {t1:.0f}s, total {time.time() - t0:.0f}s)", flush=True)
    check(bool(arb(slo) > arb(SLOPE[0])) and bool(arb(shi) < arb(SLOPE[1])),
          f"DO.1: (log F_m)' in [{slo:.6f}, {shi:.6f}] for m = 0..{MMAX + 1} ({NC} cells, arb jets): every F_m lies in "
          f"C_{SLOPE[1]}, inside C_2 and C_100 with slope < a (CF.3's hypothesis; DI.1 applies)")
    ok = True
    for m in range(MMAX + 1):
        r = res[m]
        for a_num in AS:
            ok = ok and bool(r["up"][a_num] < arb(UPPER[a_num][m])) and bool(r["low"][a_num] > arb(LOWER[a_num][m]))
            ok = ok and bool(r["low"][a_num] <= r["up"][a_num]) and bool(r["osc"] <= r["low"][a_num])
        ok = ok and bool(r["osc"] > arb(OSC_LOWER[m]))
    check(ok, f"DO.2 ({NC} cells): for m = 0..5, " + "; ".join(
        f"m={m}: {LOWER[2][m]} < Theta_2 < {UPPER[2][m]}, {LOWER[100][m]} < Theta_100 < {UPPER[100][m]}, osc log R_m > "
        f"{OSC_LOWER[m]}" for m in range(MMAX + 1)) + f" (DI.1/DI.4 on arb jet tables, {NC} cells; floors arb)")
    RESULT.update(tab=res, slo=slo, shi=shi)
    return True


# ----------------------------------------------------------------- DO.3 level sums
def block_min(a_num, jmax):
    """G_a(j) for j = 0..jmax: min over block counts with total length <= j of the product of the factors (arb uppers)."""
    G = [arb(1)] * (jmax + 1)
    for j in range(1, jmax + 1):
        best = G[j - 1]
        for n, th in PAIRS[a_num]:
            if j >= n:
                c = arb((G[j - n] * arb(th)).upper())
                if bool(c.upper() < best.upper()):
                    best = c
        G[j] = best
    return G


def B_new(k, G):
    res = RESULT["tab"]
    best = None
    for m in range(0, min(MMAX, k - 4) + 1):
        for a_num in AS:
            c = arb((res[m]["up"][a_num] * G[a_num][k - 4 - m]).upper())
            if best is None or bool(c.upper() < best.upper()):
                best = c
    return best


def B_old(k):
    """min(DN.3, DK.4): Theta_N' min_pairs theta^floor((k-4)/n) on C_2, Theta_I' min(theta_G^.., theta_K^..) on C_100."""
    m = k - 4
    c2 = min((arb(th) ** (m // n)).upper() for n, th in PAIRS[2])
    c100 = min((arb(th) ** (m // n)).upper() for n, th in PAIRS[100])
    return arb(min((arb(TN_DN) * arb(c2)).upper(), (arb(TI_DI) * arb(c100)).upper()))


def consequences():
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    import gasket_fractional_profile_pure_point as cc
    KM = 2000
    G = {a: block_min(a, KM) for a in AS}
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    rl = arb(r_lo.numerator) / r_lo.denominator
    rh = arb(r_hi.numerator) / r_hi.denominator
    found, sums, counts, okb = cc.period_atoms(bk.KMAX, verbose=False)
    ctx.prec = PREC
    okd, worst = okb, arb(0)
    for k in range(4, bk.KMAX):
        dev = arb(max(abs((sums[k + 1] / (rl * sums[k])).log()).upper(), abs((sums[k + 1] / (rh * sums[k])).log()).upper()))
        b = B_new(k, G)
        okd = okd and bool(dev < b)
        worst = arb(max(worst.upper(), (dev / b).upper()))
    facs = {k: B_old(k) / B_new(k, G) for k in range(4, KM + 4 - 9)}
    # all k: B_new(k+9) <= theta_9 B_new(k) (G(j+9) <= theta_9 G(j) by definition of G); for k >= K0 the theta_9 term is the
    # minimum in B_old at k (sufficient: (x/9 - 8/9) log theta_9 + log Theta_N' <= (x/n) log theta_n + log C_n, x = k - 4,
    # for every competitor (n, theta_n, C_n)), so B_old(k+9) = theta_9 B_old(k) there and facs(k+9) >= facs(k).
    l9, lTN = arb(PAIRS[2][2][1]).log(), arb(TN_DN).log()
    x0 = arb(0)
    for (n, th, C) in ((5, PAIRS[2][0][1], TN_DN), (7, PAIRS[2][1][1], TN_DN), (8, PAIRS[100][0][1], TI_DI),
                       (9, PAIRS[100][1][1], TI_DI)):
        sl = l9 / 9 - arb(th).log() / n                       # < 0
        assert bool(sl < 0)
        rhs = arb(C).log() - lTN + l9 * 8 / 9                 # need x sl <= rhs
        x0 = arb(max(x0.upper(), (rhs / sl).upper()))
    K0 = 4 + int(math.ceil(float(x0.upper())))
    assert K0 + 9 < KM + 4 - 9, K0
    fmin_k = min((k for k in facs if k >= 5), key=lambda k: facs[k].lower())
    fmin_tail = min((facs[k].lower() for k in facs if k >= 14))
    show = (5, 6, 7, 8, 9, 13, 22, 49, 103)
    check(okd and all(bool(facs[k] >= 1) for k in facs) and all(bool(facs[k] > 1) for k in facs if k >= 5),
          "DO.3: for every k >= 4, |log(s_(k+1)/(r s_k))| <= B(k) = min_(m <= min(5, k-4), a in {2, 100}) Theta_a(F_m, F_(m+1)) "
          "G_a(k-4-m) (CF.1 with mixed blocks; thetas as stated). B(k) values: "
          + ", ".join(f"k={k}: {B_new(k, G).str(5)}" for k in show)
          + "; min(DN.3, DK.4)/B(k): " + ", ".join(f"{facs[k].str(4)} (k={k})" for k in show)
          + f"; smallest factor over 5 <= k <= {KM - 5}: {facs[fmin_k].str(4)} at k = {fmin_k}; over 14 <= k <= {KM - 5}: "
          f">= {arb(fmin_tail).str(4)}; for k >= K0 = {K0} the factor is non-decreasing along k -> k + 9, so both minima hold "
          f"for every k (checked range covers K0 + 9). Consistency: on CC's s_4..s_{bk.KMAX} every deviation is below B(k) (max ratio "
          f"{worst.str(4)})")
    RESULT.update(K0=K0, facs=facs, worst=worst, fmin_k=fmin_k, fmin=facs[fmin_k], fmin_tail=arb(fmin_tail),
                  B={k: B_new(k, G) for k in show})
    # DO.4 kept failure
    res = RESULT["tab"]
    r0 = arb(TI_DI) / res[0]["osc"]
    rats = [res[m]["up"][100] / res[m]["osc"] for m in range(MMAX + 1)]
    c2bad = [res[m]["low"][2] / res[m]["up"][100] for m in range(MMAX + 1)]
    check(bool(r0 < arb("1.0252")) and all(bool(x > 1) for x in c2bad),
          f"DO.4 (kept failure): for k = 4 + m, every CF.3-form bound in any cone of positive functions is >= osc log R_m "
          f"(Theta_C >= Hilbert metric of the positive cone); m = 0: osc log R_0 > {res[0]['osc'].str(7)}, so DI.4's "
          f"{TI_DI} is <= {r0.str(5)} x any cone's constant at k = 4; Theta_100/osc <= "
          + ", ".join(x.str(5) for x in rats) + " (m = 0..5); Theta_2/Theta_100 >= "
          + ", ".join(x.str(4) for x in c2bad) + ": C_2 is never the better cone for the direct terms")
    RESULT.update(r0=r0, rats=rats, c2bad=c2bad)


def main(T=None) -> int:
    t0 = time.time()
    print("Theorem DO: direct f_m terms Theta_a(Lt^m f_0, Lt^(m+1) f_0), m = 0..5, a = 2, 100, and level sums. Not thrust.",
          flush=True)
    if T is None:
        T = di.taylor()
    if not main_tables(T):
        print("FAILED")
        return 1
    consequences()
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(di.FAILS) + len(bk.FAILS) + len(bk.co.FAILS) + len(bk.cn.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DO checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
