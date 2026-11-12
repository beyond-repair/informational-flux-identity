#!/usr/bin/env python3
r"""Theorem DT: the level-sum bound of CF.3/DN.3/DO.3 form carried to the second-order cone K_(2,6) (DR's and DS's kept
limit "the level-sum bound of CF.3/DN.3/DO.3 form on K_(2,6) is still not computed (it also needs f_0 in K_(2,6) and
Theta(f_0, Lt f_0))").

Setting (GASKET.md, Theorems CF, DO, DQ-DS): s in X = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
K = K_(2,6) = {f > 0 : |f'| <= 2f, |f''| <= 6f} (DQ.1: Lt K in K), cut out by the functionals
ell in {delta_s, 2 delta_s +- delta_s', 6 delta_s +- delta_s''}, so Theta_K(u, v) = log(sup_ell rho_ell / inf_ell rho_ell),
rho_ell = ell(v)/ell(u) (DQ's setting, as in DP.1).  f_0 = F0 = Lt^4 G (DC.2, CO.1/CO.2's Taylor data re-run unchanged),
F_m = Lt^m F0, R_m = F_(m+1)/F_m, s_(4+m) = c0 F_m(x0) (DI.5's consistency).

Theorem DT (checked here; exit 1 on any failure):
  DT.1 [certified; arb jets of order 4 (f, f', ..., f''''/24) through m + 1 <= 6 layers of Lt on DT_CELLS cells]
       |F_m'| / F_m < 0.80 and 0.20 < F_m''/F_m < 1.69 on X for m = 0..5: every F_m lies in the interior of K_(2,6),
       every ell(F_m) > 0, so CF.3's argument can start from each F_m in K.
  DT.2 [certified; per cell rho_ell(J) in rho_ell(c) + rho_ell'(c)(J - c) + [rho_ell''/2](J)(J - c)^2 intersected with the
       direct quotient, all arb;
       floors: arb point values at DT_FLOOR_PTS points, endpoints included] two-sided bounds for Theta_K(F_m, F_(m+1)),
       m = 0..5 (table in GASKET.md).
  DT.3 [proved from DT.2 with CF.1, CF.3's proof, DQ.1 and the K-contraction factors theta_n = tanh(Dbar_n/4) of DR.3
       (n = 3, 5..8) and DS.3 (n = 4), used as stated] for every k >= 4,
         |log(s_(k+1)/(r s_k))| <= B_K(k) = min over m <= min(5, k-4) of Theta_K(F_m, F_(m+1)) G_K(k-4-m),
       G_K(j) the least product of theta_n (n = 3..8, mixed) with block lengths summing to at most j.  Per-level rate
       theta_7^(1/7) <= 0.2694000 (DO.3: theta_9^(1/9) < 0.363657).  B_K(k) < B_DO(k) (DO.3, constants as stated) for
       every k >= 11, by a factor that grows at least like (0.36365/0.26940)^k.
  DT.4 [kept failure, certified witness] for k = 4..10 DO.3 stays the bound: Theta_K >= Theta_2 >= Theta_100 (K in C_2 in
       C_100) and the K-floors exceed DO's constants there; the second-order constraints are not what helps (for m >= 2,
       Theta_K(F_m, F_(m+1)) agrees with DO's Theta_2 to within the enclosures); the gain is DR/DS's contraction factors.
Kept limits: see GASKET.md Theorem DT.  No Psi_info, Stage 2 still blocked.  Stage 1 untouched; no W, 0.08, 0.23 or
30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time
from math import comb

from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_f0_coupled as di          # noqa: E402
import gasket_fractional_profile_level_f0_shift as do            # noqa: E402

bk = di.bk
FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
FULL = os.environ.get("DT_FULL", "") not in ("", "0")
NC = int(os.environ.get("DT_CELLS", "4096" if FULL else "1024"))
NFP = int(os.environ.get("DT_FLOOR_PTS", "4097" if FULL else "1025"))
MMAX = 5
KJ = 5                                    # jet length: f, f', f''/2, f'''/6, f''''/24
KM = 2000                                 # range of k checked directly (k <= KM)
# K-contraction factors: Dbar_n = certified upper bounds for D_n = diam Lt^n K_(2,6), used as stated
DBAR = {3: "0.25513308", 4: "0.03838589", 5: "0.0098239159", 6: "0.0018377053", 7: "0.00041194852", 8: "0.00017542239"}
# DO.3's constants, used as stated (DO.2's 4096-cell table and DO's contraction pairs)
DO_UPPER = do.TARGETS[4096]["UPPER"]
DO_PAIRS = do.PAIRS
# functionals ell = a0 delta + a1 delta' + a2 delta''
ELLS = ((1, 0, 0), (2, 1, 0), (2, -1, 0), (6, 0, 1), (6, 0, -1))
ELLN = ("delta", "2delta+delta'", "2delta-delta'", "6delta+delta''", "6delta-delta''")
# certified targets (checked, exit 1 if not met)
TARGETS = {
    1024: dict(UPPER=("0.147314", "0.0723405", "0.0069153", "0.0010009", "0.00016489", "0.000028163")),
    4096: dict(UPPER=("0.147314", "0.0723405", "0.0069153", "0.0010009", "0.00016489", "0.000028163")),
}
LOWER = ("0.147313", "0.0723403", "0.0069152", "0.0010008", "0.00016484", "0.000028125")
D1_RANGE = ("0.80", "0.20", "1.69")       # |F'/F| <, F''/F >, F''/F <


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- jets of order 4 with arb coefficients
def jadd(a, b):
    return [x + y for x, y in zip(a, b)]


def jmul(a, b):
    return [sum((a[i] * b[k - i] for i in range(k + 1)), arb(0)) for k in range(KJ)]


def jrec(a):
    b = [1 / a[0]]
    for k in range(1, KJ):
        b.append(-sum((a[i] * b[k - i] for i in range(1, k + 1)), arb(0)) * b[0])
    return b


def jsqrt(a):
    b = [a[0].sqrt()]
    for k in range(1, KJ):
        b.append((a[k] - sum((b[i] * b[k - i] for i in range(1, k)), arb(0))) / (2 * b[0]))
    return b


class TaylorJet4:
    """F0(s) = sum f_n u^n, u = (s - c)/R (CO.2): F0^(j)(u)/j! = sum_n f_n C(n,j) u^(n-j) plus DO.1's Cauchy tail
    S rho^-j q^(M-j) sum_(l<=j) C(M,j-l) q^l/(1-q)^(l+1) (q = |u|/rho, M = N+1), j = 0..4, composed with the input jet."""

    def __init__(self, T):
        self.f, self.S, self.rho, self.c, self.R = T.f, T.S, T.rho, T.c, T.R
        self.N = len(T.f) - 1
        self.d = [[T.f[n] * comb(n, j) for n in range(j, self.N + 1)] for j in range(KJ)]

    def jet(self, X):
        U0 = (X[0] - self.c) / self.R
        q = arb(abs(U0).upper()) / self.rho
        assert bool(q < 1)
        M = self.N + 1
        c = []
        for j in range(KJ):
            acc = arb(0)
            for t in reversed(self.d[j]):
                acc = acc * U0 + t
            tail = self.S * self.rho ** (-j) * q ** (M - j) * sum(
                (comb(M, j - l) * q ** l / (1 - q) ** (l + 1) for l in range(j + 1)), arb(0))
            c.append(acc + arb(0, tail.upper()))
        dU = [arb(0)] + [x / self.R for x in X[1:]]
        out, p = [c[0]] + [arb(0)] * (KJ - 1), [arb(1)] + [arb(0)] * (KJ - 1)
        for j in range(1, KJ):
            p = jmul(p, dU)
            out = jadd(out, [c[j] * x for x in p])
        return out


def Fjet(TJ, m, X):
    """Jet of F_m = Lt^m F0 at the jet X (DO's Fjet at order 4):
    wt1 = (s+1/2)/(2s(7/2-s)), wt2 = (s-1/2)/(2s(s+7/2)), sigma_1 = sqrt(15/4 - s), sigma_2 = sqrt(15/4 + s)."""
    if m == 0:
        return TJ.jet(X)
    q4, half, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    neg = [-x for x in X]

    def sh(a, Y):
        return [a + Y[0]] + Y[1:]
    rX = [half * x for x in jrec(X)]
    w1 = jmul(jmul(sh(half, X), rX), jrec(sh(s7, neg)))
    w2 = jmul(jmul(sh(-half, X), rX), jrec(sh(s7, X)))
    s1, s2 = jsqrt(sh(q4, neg)), jsqrt(sh(q4, X))
    return jadd(jmul(w1, Fjet(TJ, m - 1, s1)), jmul(w2, Fjet(TJ, m - 1, s2)))


def ell(a, j):
    """Taylor jet (e_0, e_1, e_2), e_k = ell(F)^(k)/k!, from the Taylor jet j = (F^(i)/i!)_(i<=4):
    ell(F)^(k) = a0 F^(k) + a1 F^(k+1) + a2 F^(k+2)."""
    d = [j[i] * math.factorial(i) for i in range(KJ)]
    return [(a[0] * d[k] + a[1] * d[k + 1] + a[2] * d[k + 2]) / math.factorial(k) for k in range(3)]


def jdiv3(A, B):
    """Taylor jet of A/B to order 2 (A, B jets of length 3, B[0] > 0)."""
    r0 = A[0] / B[0]
    r1 = (A[1] - r0 * B[1]) / B[0]
    r2 = (A[2] - r1 * B[1] - r0 * B[2]) / B[0]
    return [r0, r1, r2]


# ----------------------------------------------------------------- DT.1 / DT.2
def theta_upper(TJ, m):
    """Upper bound for Theta_K(F_m, F_(m+1)) over NC cells covering X (outward), and DT.1's ranges for F_m."""
    ctx.prec = PREC
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    s0l, s1u = arb(s0.lower()), arb(s1.upper())
    h = (s1u - s0l) / NC
    one, zero = arb(1), arb(0)
    beta, alpha = None, None
    q1max, q2lo, q2hi, ellmin = arb(0), None, None, None
    argb = arga = None
    for i in range(NC):
        J = bk._hull(s0l + h * i, s0l + h * (i + 1))
        mp = arb(J.mid())
        dx = J - mp
        uJ, vJ = Fjet(TJ, m, [J] + [one] + [zero] * (KJ - 2)), Fjet(TJ, m + 1, [J] + [one] + [zero] * (KJ - 2))
        uM, vM = Fjet(TJ, m, [mp] + [one] + [zero] * (KJ - 2)), Fjet(TJ, m + 1, [mp] + [one] + [zero] * (KJ - 2))
        if not bool(uJ[0] > 0):
            check(False, f"DT.1: positivity not certified on cell {i} (m = {m})")
            return None
        q1 = abs(uJ[1] / uJ[0])
        q2 = 2 * uJ[2] / uJ[0]
        q1max = arb(max(q1max.upper(), q1.upper()))
        q2lo = q2.lower() if q2lo is None else min(q2lo, q2.lower())
        q2hi = q2.upper() if q2hi is None else max(q2hi, q2.upper())
        for k, a in enumerate(ELLS):
            AJ, BJ, AM, BM = ell(a, vJ), ell(a, uJ), ell(a, vM), ell(a, uM)
            if not (bool(BJ[0] > 0) and bool(AJ[0] > 0)):
                check(False, f"DT.1: ell(F) > 0 not certified on cell {i}, ell = {ELLN[k]} (m = {m})")
                return None
            ellmin = (BJ[0] / uJ[0]).lower() if ellmin is None else min(ellmin, (BJ[0] / uJ[0]).lower())
            rJ, rM = jdiv3(AJ, BJ), jdiv3(AM, BM)
            rc = rM[0] + rM[1] * dx + rJ[2] * dx * dx       # Taylor with Lagrange remainder about the midpoint
            rd = rJ[0]
            lo, hi = max(rc.lower(), rd.lower()), min(rc.upper(), rd.upper())
            if beta is None or hi > beta:
                beta, argb = hi, (i, ELLN[k])
            if alpha is None or lo < alpha:
                alpha, arga = lo, (i, ELLN[k])
    up = arb((arb(beta) / arb(alpha)).log().upper())
    return dict(up=up, q1max=q1max, q2lo=arb(q2lo), q2hi=arb(q2hi), ellmin=arb(ellmin), argb=argb, arga=arga)


def theta_floor(TJ, m):
    """Certified lower bounds (arb): Theta_K(F_m, F_(m+1)) >= log(max rho_ell / min rho_ell) over NFP points of X
    (endpoints included) and the five ell; osc log R_m >= the same with ell = delta only."""
    ctx.prec = PREC
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    one, zero = arb(1), arb(0)
    hi = lo = None
    rhi = rlo = None
    for k in range(NFP):
        s = s0 if k == 0 else (s1 if k == NFP - 1 else arb((s0 + (s1 - s0) * k / (NFP - 1)).mid()))
        u, v = Fjet(TJ, m, [s] + [one] + [zero] * (KJ - 2)), Fjet(TJ, m + 1, [s] + [one] + [zero] * (KJ - 2))
        for e, a in enumerate(ELLS):
            r = ell(a, v)[0] / ell(a, u)[0]
            hi = r.lower() if hi is None else max(hi, r.lower())
            lo = r.upper() if lo is None else min(lo, r.upper())
            if e == 0:
                rhi = r.lower() if rhi is None else max(rhi, r.lower())
                rlo = r.upper() if rlo is None else min(rlo, r.upper())
    low = arb((arb(hi) / arb(lo)).log().lower())
    osc = arb((arb(rhi) / arb(rlo)).log().lower())
    return low, osc


def main_tables(T):
    TJ = TaylorJet4(T)
    res = {}
    q1max, q2lo, q2hi, emin = arb(0), None, None, None
    for m in range(MMAX + 1):
        t0 = time.time()
        out = theta_upper(TJ, m)
        if out is None:
            return False
        q1max = arb(max(q1max.upper(), out["q1max"].upper()))
        q2lo = out["q2lo"] if q2lo is None else arb(min(q2lo.lower(), out["q2lo"].lower()))
        q2hi = out["q2hi"] if q2hi is None else arb(max(q2hi.upper(), out["q2hi"].upper()))
        emin = out["ellmin"] if emin is None else arb(min(emin.lower(), out["ellmin"].lower()))
        low, osc = theta_floor(TJ, m)
        res[m] = dict(up=out["up"], low=low, osc=osc, argb=out["argb"], arga=out["arga"])
        print(f"   m = {m}: Theta_K in [{low.str(7)}, {out['up'].str(7)}], osc log R_m >= {osc.str(7)}; beta at cell "
              f"{out['argb']}, alpha at cell {out['arga']} ({time.time() - t0:.0f}s)", flush=True)
    mm = MMAX
    check(bool(q1max < arb(D1_RANGE[0])) and bool(q2lo > arb(D1_RANGE[1])) and bool(q2hi < arb(D1_RANGE[2])),
          f"DT.1: |F_m'|/F_m <= {q1max.str(6)} < 2 and F_m''/F_m in [{q2lo.str(6)}, {q2hi.str(6)}] inside (-6, 6) for "
          f"m = 0..{mm} ({NC} cells, arb jets of order 4); min_ell ell(F_m)/F_m >= {emin.str(6)} > 0: every F_m lies in the "
          f"interior of K_(2,6) (CF.3's hypothesis on K)")
    up = TARGETS.get(NC, TARGETS[1024] if NC < 4096 else TARGETS[4096])["UPPER"]
    ok = True
    for m in range(MMAX + 1):
        r = res[m]
        ok = ok and bool(r["up"] < arb(up[m])) and bool(r["low"] > arb(LOWER[m])) and bool(r["low"] <= r["up"])
        ok = ok and bool(r["osc"] <= r["low"])
    check(ok, f"DT.2 ({NC} cells): " + "; ".join(f"m={m}: {LOWER[m]} < Theta_K(F_m, F_(m+1)) < {up[m]}"
                                                for m in range(MMAX + 1))
          + f" (arb centred forms per cell and functional; floors arb at {NFP} points)")
    RESULT.update(tab=res, q1max=q1max, q2lo=q2lo, q2hi=q2hi, ellmin=emin)
    return True


# ----------------------------------------------------------------- DT.3 / DT.4 level sums
def thetas():
    ctx.prec = PREC
    return {n: arb((arb(D) / 4).tanh().upper()) for n, D in DBAR.items()}


def knap(pairs, jmax):
    """G(j), j = 0..jmax: least product of the factors whose block lengths sum to at most j (arb uppers)."""
    G = [arb(1)] * (jmax + 1)
    for j in range(1, jmax + 1):
        best = G[j - 1]
        for n, th in pairs:
            if j >= n:
                c = arb((G[j - n] * th).upper())
                if bool(c.upper() < best.upper()):
                    best = c
        G[j] = best
    return G


def B_K(k, G):
    res = RESULT["tab"]
    return min((arb((res[m]["up"] * G[k - 4 - m]).upper()) for m in range(0, min(MMAX, k - 4) + 1)),
               key=lambda x: x.upper())


def B_DO(k, GD):
    """DO.3's B(k) recomputed from DO's stated constants (DO.2's 4096-cell uppers, DO's contraction pairs)."""
    return min((arb((arb(DO_UPPER[a][m]) * GD[a][k - 4 - m]).upper()) for m in range(0, min(5, k - 4) + 1) for a in (2, 100)),
               key=lambda x: x.upper())


def consequences():
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    import gasket_fractional_profile_pure_point as cc
    th = thetas()
    per = {n: th[n] ** (arb(1) / n) for n in th}
    G = knap([(n, th[n]) for n in sorted(th)], KM)
    GD = {a: knap([(n, arb(t)) for n, t in DO_PAIRS[a]], KM) for a in (2, 100)}
    rate = min(per.values(), key=lambda x: x.upper())
    nbest = min(per, key=lambda n: per[n].upper())
    rho9 = arb(DO_PAIRS[2][2][1]) ** (arb(1) / 9)                     # DO.3's per-level rate (smallest of DO's pairs)
    rmin_do = min((arb(t) ** (arb(1) / n) for a in (2, 100) for n, t in DO_PAIRS[a]), key=lambda x: x.lower())
    # consistency with CC's level sums
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    rl = arb(r_lo.numerator) / r_lo.denominator
    rh = arb(r_hi.numerator) / r_hi.denominator
    found, sums, counts, okb = cc.period_atoms(bk.KMAX, verbose=False)
    ctx.prec = PREC
    okd, worst = okb, arb(0)
    for k in range(4, bk.KMAX):
        dev = arb(max(abs((sums[k + 1] / (rl * sums[k])).log()).upper(), abs((sums[k + 1] / (rh * sums[k])).log()).upper()))
        b = B_K(k, G)
        okd = okd and bool(dev < b)
        worst = arb(max(worst.upper(), (dev / b).upper()))
    facs = {k: B_DO(k, GD) / B_K(k, G) for k in range(4, KM + 4 - 10)}
    # tail: G_a(j) >= rmin_do^j for DO's pairs (every factor theta_n >= rmin_do^n), so B_DO(k) >= c_DO rmin_do^(k-4) with
    # c_DO = min_(m, a) DO_UPPER[a][m] rmin_do^(-m) (lower ends); and B_K(k) <= Theta_K(F_5, F_6) theta_7^floor((k-9)/7)
    # <= Theta_K(F_5, F_6) theta_7^((k-15)/7).  So facs(k) >= exp(A + S (k - 4)), S = log rmin_do - log theta_7 / 7 > 0.
    t7 = th[7]
    lr = arb(rmin_do.lower()).log()
    cDO = min((arb(DO_UPPER[a][m]) * arb(rmin_do.lower()) ** (-m) for m in range(6) for a in (2, 100)),
              key=lambda x: x.lower())
    up5 = RESULT["tab"][5]["up"]
    S = lr - t7.log() / 7
    A = cDO.log() - up5.log() - t7.log() * (-11) / 7               # (k-15)/7 = (k-4)/7 - 11/7
    check(bool(S > 0), f"DT.3 tail slope: log(DO rate) - log(theta_7)/7 = {S.str(6)} > 0")
    # smallest k_t with A + S (k_t - 4) > 0, i.e. facs(k) > 1 for every k >= k_t
    kt = 4 + max(0, int(math.ceil(float((-A / S).upper()))))
    assert kt < KM - 20, kt
    kstar = None
    for k in range(4, KM + 4 - 10):
        if all(bool(facs[kk] > 1) for kk in range(k, KM + 4 - 10)):
            kstar = k
            break
    show = (5, 6, 7, 8, 9, 10, 11, 13, 22, 49, 103)
    ok3 = okd and kstar is not None and kt <= KM
    fmin_tail = min((facs[k].lower() for k in facs if kstar is not None and k >= kstar))
    check(ok3,
          "DT.3: for every k >= 4, |log(s_(k+1)/(r s_k))| <= B_K(k) = min_(m <= min(5, k-4)) Theta_K(F_m, F_(m+1)) G_K(k-4-m) "
          "(CF.1 on K_(2,6) with mixed blocks, theta_n = tanh(Dbar_n/4): " + ", ".join(f"theta_{n} <= {th[n].str(8)}" for n in sorted(th))
          + f"); per-level rate min theta_n^(1/n) = theta_{nbest}^(1/{nbest}) <= {rate.str(8)} (DO.3: theta_9^(1/9) <= "
          f"{rho9.str(8)}). B_K(k): " + ", ".join(f"k={k}: {B_K(k, G).str(5)}" for k in show)
          + "; B_DO(k)/B_K(k): " + ", ".join(f"{facs[k].str(4)} (k={k})" for k in show)
          + f"; B_K(k) < B_DO(k) for every k >= {kstar} (checked for k <= {KM - 7}; for k >= {kt} by the tail bound "
          f"facs(k) >= exp({A.str(4)} + {S.str(5)} (k - 4))); smallest factor over k >= {kstar}: {arb(fmin_tail).str(4)}. "
          f"Consistency: on CC's s_4..s_{bk.KMAX} every deviation is below B_K(k) (max ratio {worst.str(4)})")
    RESULT.update(th=th, per=per, rate=rate, nbest=nbest, rho9=rho9, facs=facs, kstar=kstar, kt=kt, S=S, A=A,
                  worst=worst, B={k: B_K(k, G) for k in show}, BDO={k: B_DO(k, GD) for k in show}, fmin_tail=arb(fmin_tail))
    # DT.4 kept failure: for k < kstar the K route does not beat DO.3; the K-floors against DO's constants
    res = RESULT["tab"]
    DO_LOWER2 = do.TARGETS[4096]["LOWER"][2]
    upK_vs_do = [res[m]["up"] / arb(DO_LOWER2[m]) for m in range(MMAX + 1)]
    lowK_vs_do = [res[m]["low"] / arb(DO_UPPER[2][m]) for m in range(MMAX + 1)]
    lowK_vs_100 = [res[m]["low"] / arb(DO_UPPER[100][m]) for m in range(MMAX + 1)]
    # lower bound for every K-route bound at k (direct terms only up to m <= 5): min_m low_K[m] * Gmin(k-4-m), with
    # Gmin(j) >= min over block combos; compare with B_DO(k)
    lowB = {}
    for k in range(4, kstar if kstar else 4):
        lb = min((arb((res[m]["low"] * G[k - 4 - m]).lower()) for m in range(0, min(MMAX, k - 4) + 1)),
                 key=lambda x: x.lower())
        lowB[k] = lb
    okf = all(bool(lowB[k] > B_DO(k, GD)) for k in lowB)
    check(okf and bool(lowK_vs_do[0] > 1) and bool(lowK_vs_do[1] > 1),
          f"DT.4 (kept failure): for k = 4..{(kstar or 5) - 1} every bound of DT.3's form (direct terms m <= 5, K-factors as "
          f"stated) exceeds B_DO(k): K-route floors / B_DO(k) >= "
          + ", ".join(f"{(lowB[k] / B_DO(k, GD)).str(4)} (k={k})" for k in lowB)
          + "; Theta_K-floor / DO's Theta_2-upper: " + ", ".join(x.str(5) for x in lowK_vs_do[:2])
          + " (m = 0, 1: the second-order constraints bind); Theta_K-upper / DO's Theta_2-floor: "
          + ", ".join(x.str(6) for x in upK_vs_do[2:]) + " (m = 2..5: they do not, within the enclosures)"
          + "; Theta_K-floor / DO's Theta_100-upper: " + ", ".join(x.str(4) for x in lowK_vs_100)
          + ". K in C_2 in C_100, so Theta_K >= Theta_2 >= Theta_100; for m >= 2 Theta_K and Theta_2 agree within the "
          "enclosures: the gain of DT.3 is DR/DS's contraction factors, not a smaller constant")
    RESULT.update(lowB=lowB, lowK_vs_do=lowK_vs_do, lowK_vs_100=lowK_vs_100, upK_vs_do=upK_vs_do)
    # corollary: K_(2,6) in C_2, so Theta_2 <= Theta_K and DT.2's uppers bound Theta_2(F_m, F_(m+1)) as well
    sharper = [m for m in range(MMAX + 1) if bool(res[m]["up"] < arb(DO_UPPER[2][m]))]
    check(sharper == [2, 3, 4, 5],
          "DT.2 corollary (K_(2,6) in C_2, so Theta_2 <= Theta_K): Theta_2(F_m, F_(m+1)) <= "
          + ", ".join(f"{res[m]['up'].str(7)} (m={m}, DO.2: {DO_UPPER[2][m]})" for m in sharper)
          + ": DT.2's uppers sharpen DO.2's 4096-cell Theta_2 uppers for m = 2..5")
    RESULT.update(sharper=sharper)


def main(T=None) -> int:
    t0 = time.time()
    print(f"Theorem DT: the CF.3/DO.3 level-sum bound on K_(2,6) ({NC} cells, m = 0..{MMAX}{', FULL' if FULL else ''}). "
          "Not thrust.", flush=True)
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
    print("all Theorem DT checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
