#!/usr/bin/env python3
r"""Theorem DI: DH's kept limit "Theta_H' uses DC.2's symmetric slope bounds for f_0, Lt f_0".
The f_0 term Theta_100(f_0, Lt f_0) of DH.5 is recomputed: signed slopes of f_0 and Lt f_0 (DI.2), DH.1 with them (DI.3),
and a coupled bound from an exact identity for the Hilbert-metric quotient (DI.1, DI.4), with a certified floor.

Setting (GASKET.md, Theorems CE, CF, DC-DH): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
C_a = {f > 0 : log f a-Lipschitz in s}, Theta_a its Hilbert metric (CF's definition), a = 100; f_0 = F0 = Lt^4 G
(DC.2, CO.1/CO.2's Taylor data re-run unchanged), v = Lt F0, R = v/u = Lt F0/F0 (DC.2's R_0).

Theorem DI (checked here; exit 1 on any failure):
  DI.1 [proved; exact] For u, v > 0 with Phi_u(x, y) = e^(a|x-y|) u(y) - u(x) > 0 and R = v/u,
           Phi_v(x, y)/Phi_u(x, y) = R(y) + (R(y) - R(x)) / (e^(a|x-y|) u(y)/u(x) - 1),
       so Theta_a(u, v) = log(beta/alpha) with beta = sup, alpha = inf of this quotient over x != y together with R.
  DI.2 [certified; arb, interval Horner on DC.2's cells] signed slopes of log F0 and log Lt F0 (d/ds) and an enclosure
       of R' = R ((log Lt F0)' - (log F0)').
  DI.3 [certified; arb] DH.1 / DG.1(i) with the signed bounds of DI.2 (asymmetric factor lemma, signed form).
  DI.4 [certified; cell pairs, outward float arithmetic] DI.1 with (R(y) - R(x)) bounded by the mean value theorem on
       the cells between x and y and by sup/inf of R on the cells, and the denominator by DI.2's slopes of log F0:
       Theta_100(f_0, Lt f_0) <= Theta_I'.  Floor (certified, arb): the limits x -> y+- of DI.1 at points.
  DI.5 [proved from DI.4 with DH.5's proof, DG.3's theta_G and DH.4's theta_H as stated] level sums.
Kept limits: see GASKET.md Theorem DI.  No Psi_info, Stage 2 still blocked.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_birkhoff as bk          # noqa: E402

FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A_NUM = 100
NC = int(os.environ.get("DI_CELLS", "16384"))
THETA_G, THETA_H, RATE_H = "0.014183", "0.0064696", "0.57117"   # DG.3, DH.4 (used as stated)
TH_DH = "0.10465"                                                # DH.5's Theta_H'
B0_DC, b0_DC, LOGR0_DC = "1.2791803", "0.77333706", "0.063592701"
# certified targets (set after the first run; checked, exit 1 if not met)
Q0_RANGE = ("0.035291", "0.77334")     # DI.2: (log F0)' in (0.035291, 0.77334)
Q1_RANGE = ("0.17682", "0.79038")      # DI.2: (log Lt F0)' in (0.17682, 0.79038)
LOGR_MAX = "0.063387"                  # DI.2: log(sup R_0/inf R_0) (DC.2: 0.063592701; DC.4's floor 0.0633864)
TS_MAX = "0.078489"                    # DI.3
TI_MAX = "0.064982"                    # DI.4
TFLOOR_MIN = "0.064975"                # DI.4 floor


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _fu(x: arb) -> float:
    """float >= x.upper()"""
    return math.nextafter(float(x.upper()), math.inf)


def _fl(x: arb) -> float:
    """float <= x.lower()"""
    return math.nextafter(float(x.lower()), -math.inf)


def taylor():
    S_rho, f = bk.taylor_data()
    return bk.Taylor(f, S_rho)


# ----------------------------------------------------------------- DI.2 cell tables
def cell_tables(T):
    """Per cell J of DC.2's cover: enclosures of R = Lt F0/F0, q0 = (log F0)', q1 = (log Lt F0)', R' = R (q1 - q0)."""
    ctx.prec = PREC
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    s0l, s1u = arb(s0.lower()), arb(s1.upper())
    h = (s1u - s0l) / NC
    q4, half, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    tab = np.zeros((NC, 8))
    for i in range(NC):
        J = bk._hull(s0l + h * i, s0l + h * (i + 1))
        F, dF = T.val(J), T.der(J)
        t1, t2 = (q4 - J).sqrt(), (q4 + J).sqrt()
        F1, F2 = T.val(t1), T.val(t2)
        if not (bool(F > 0) and bool(F1 > 0) and bool(F2 > 0)):
            check(False, f"DI.2: positivity not certified on cell {i}")
            return None
        w1, w2 = bk.wts(J)
        l1 = 1 / (J + half) - 1 / J + 1 / (s7 - J)          # (log wt1)'
        l2 = 1 / (J - half) - 1 / J - 1 / (J + s7)          # (log wt2)'
        g1 = l1 - T.der(t1) / F1 / (2 * t1)                 # (log(wt1 F0 o sigma_1))', sigma_1' = -1/(2 sigma_1)
        g2 = l2 + T.der(t2) / F2 / (2 * t2)                 # sigma_2' = 1/(2 sigma_2)
        p1 = 1 / (1 + (w2 * F2) / (w1 * F1))                 # weight of word 1 in Lt F0
        q1 = g2 + p1 * (g1 - g2)
        q0 = dF / F
        R = (w1 * F1 + w2 * F2) / F
        Rp = R * (q1 - q0)
        # centred form: R(J) in R(m) + R'(J)(J - m), m the midpoint (mean value theorem); intersect with the direct enclosure
        mpt = arb(J.mid())
        Rm = bk.LtF(T, mpt) / T.val(mpt)
        Rc = Rm + Rp * (J - mpt)
        lo, hi = max(R.lower(), Rc.lower()), min(R.upper(), Rc.upper())
        R = bk._hull(arb(lo), arb(hi))
        tab[i] = (_fl(R), _fu(R), _fl(q0), _fu(q0), _fl(q1), _fu(q1), _fl(Rp), _fu(Rp))
    return tab, h


# ----------------------------------------------------------------- DI.3 signed factor lemma (DG.1(i) per orientation)
def _factor(a, al, be):
    """DG.1(i): (a+beta)/(a-alpha) under its hypotheses, else None."""
    ok = bool(abs(al) < a) and bool(abs(be) < a) and bool(al + be >= 0)
    if bool(be < 0):
        ok = ok and bool(al <= a / 2)
    if bool(al < 0):
        ok = ok and bool(be >= abs(al))
    return (a + be) / (a - al) if ok else None


def signed_factor(a, cu, du, cv, dv):
    """Theta_a(u, v) <= log(sup R/inf R) + log K_up + log K_dn when c_u (x-y) <= log u(x)/u(y) <= d_u (x-y) (x > y), same
    for v.  K_up bounds Phi_v/Phi_u / (v(y)/u(y)): x > y uses (alpha, beta) = (d_u, -c_v), x < y uses (-c_u, d_v);
    K_dn is the same with u, v exchanged."""
    f = [_factor(a, du, -cv), _factor(a, -cu, dv), _factor(a, dv, -cu), _factor(a, -cv, du)]
    if any(x is None for x in f):
        return None
    return arb(max(f[0].upper(), f[1].upper())), arb(max(f[2].upper(), f[3].upper()))


# ----------------------------------------------------------------- DI.4 coupled bound from DI.1 on cell pairs
def _mvfac(x):
    """Upper bound for x/(e^x - 1) (decreasing in x >= 0; 1 at x = 0), padded outward."""
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(x > 1e-8, x / (np.expm1(x) * (1 - 1e-12)), 1.0)
    return np.minimum(f * (1 + 1e-12), 1.0)


def coupled(tab, h, a):
    """Upper bound for beta = sup(R(y) + E) and lower bound for alpha = inf(R(y) + E), E = (R(y)-R(x))/(e^(ad) u(y)/u(x) - 1),
    over all x != y in X_s; d = |x - y|, u = F0.  With q0 = (log F0)' on the cells spanned by [x, y]:
    x > y: u(x)/u(y) <= e^(d sup q0), so den >= e^(k d) - 1 with k = a - sup_span q0;
    x < y: u(x)/u(y) <= e^(-d inf q0), so den >= e^(k d) - 1 with k = a + inf_span q0.
    Numerator: mean value theorem on the cells spanned, |E| <= sup_span |R'| d/(e^(k d) - 1) <= (sup_span |R'|/k) (k g)/(e^(k g) - 1)
    for d >= g (x/(e^x - 1) decreasing), or sup/inf of R on the two cells when they are separated by a gap g > 0
    (then den >= e^(k g) - 1)."""
    Rlo, Rhi, q0lo, q0hi, Rplo, Rphi = tab[:, 0], tab[:, 1], tab[:, 2], tab[:, 3], tab[:, 6], tab[:, 7]
    af = float(a.lower())
    hf = float(h.lower()) * (1 - 1e-15)
    pos = np.maximum(Rphi, 0.0)            # sup (R')_+ per cell
    neg = np.maximum(-Rplo, 0.0)           # sup (-R')_+ per cell
    n = len(Rlo)
    up = np.zeros(n)
    dn = np.zeros(n)
    pad = 1 + 1e-12
    gaps = np.arange(n) * hf               # gap between cells j and j +- m is (m - 1) h; index m-1
    for j in range(n):
        # same cell: x > y with k = a - q0hi_j, x < y with k = a + q0lo_j
        kr0 = (af - q0hi[j]) * (1 - 1e-15)
        kl0 = (af + q0lo[j]) * (1 - 1e-15)
        e_up = max(neg[j] / kr0, pos[j] / kl0)
        e_dn = max(pos[j] / kr0, neg[j] / kl0)
        # x in cells to the right (x > y): R(y) - R(x) <= d sup(-R')_+ on the span, R(x) - R(y) <= d sup(R')_+
        if j + 1 < n:
            mneg = np.maximum.accumulate(np.maximum(neg[j:], neg[j]))[1:]
            mpos = np.maximum.accumulate(np.maximum(pos[j:], pos[j]))[1:]
            k = (af - np.maximum.accumulate(q0hi[j:])[1:]) * (1 - 1e-15)
            g = gaps[: n - j - 1]
            den = np.where(g > 0, np.expm1(k * g) * (1 - 1e-12), 0.0)
            fac = _mvfac(k * g) / k
            with np.errstate(divide="ignore", invalid="ignore"):
                gu = np.where(den > 0, np.maximum(Rhi[j] - Rlo[j + 1:], 0) / den, np.inf)
                gd = np.where(den > 0, np.maximum(Rhi[j + 1:] - Rlo[j], 0) / den, np.inf)
            e_up = max(e_up, float(np.max(np.minimum(mneg * fac, gu))))
            e_dn = max(e_dn, float(np.max(np.minimum(mpos * fac, gd))))
        # x in cells to the left (x < y): R(y) - R(x) <= d sup(R')_+, R(x) - R(y) <= d sup(-R')_+
        if j > 0:
            mpos = np.maximum.accumulate(np.maximum(pos[j::-1], pos[j]))[1:]
            mneg = np.maximum.accumulate(np.maximum(neg[j::-1], neg[j]))[1:]
            k = (af + np.minimum.accumulate(q0lo[j::-1])[1:]) * (1 - 1e-15)
            g = gaps[:j]
            den = np.where(g > 0, np.expm1(k * g) * (1 - 1e-12), 0.0)
            fac = _mvfac(k * g) / k
            with np.errstate(divide="ignore", invalid="ignore"):
                gu = np.where(den > 0, np.maximum(Rhi[j] - Rlo[j - 1::-1], 0) / den, np.inf)
                gd = np.where(den > 0, np.maximum(Rhi[j - 1::-1] - Rlo[j], 0) / den, np.inf)
            e_up = max(e_up, float(np.max(np.minimum(mpos * fac, gu))))
            e_dn = max(e_dn, float(np.max(np.minimum(mneg * fac, gd))))
        up[j] = e_up * pad
        dn[j] = e_dn * pad
    beta = float(np.max((Rhi + up) * pad))
    alpha = float(np.min((Rlo - dn) / pad))
    jb, ja = int(np.argmax(Rhi + up)), int(np.argmin(Rlo - dn))
    return beta, alpha, jb, ja, float(up[jb]), float(dn[ja])


def floor_points(T, a, m=int(os.environ.get("DI_FLOOR_PTS", "2049"))):
    """Certified lower bound for Theta_a(f_0, Lt f_0): beta >= Phi_v/Phi_u -> R(y)(a - q1)/(a - q0) as x -> y+ and
    R(y)(a + q1)/(a + q0) as x -> y-; alpha <= the same; points s_k, k = 0..m-1, equally spaced in X_s."""
    ctx.prec = PREC
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    q4, half, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    blo, ahi = None, None
    best = None
    for k in range(m):
        if k == 0 or k == m - 1:
            s = s0 if k == 0 else s1          # the endpoints themselves (balls; the candidate enclosures contain the values)
        else:
            s = arb((s0 + (s1 - s0) * k / (m - 1)).mid())
        F, dF = T.val(s), T.der(s)
        t1, t2 = (q4 - s).sqrt(), (q4 + s).sqrt()
        F1, F2 = T.val(t1), T.val(t2)
        w1, w2 = bk.wts(s)
        v = w1 * F1 + w2 * F2
        l1 = 1 / (s + half) - 1 / s + 1 / (s7 - s)
        l2 = 1 / (s - half) - 1 / s - 1 / (s + s7)
        dv = w1 * F1 * (l1 - T.der(t1) / F1 / (2 * t1)) + w2 * F2 * (l2 + T.der(t2) / F2 / (2 * t2))
        q0, q1, R = dF / F, dv / v, v / F
        cands = [R]
        if k < m - 1:
            cands.append(R * (a - q1) / (a - q0))        # x -> y+
        if k > 0:
            cands.append(R * (a + q1) / (a + q0))        # x -> y-
        for c in cands:
            if blo is None or bool(c.lower() > blo.lower()):
                blo = arb(c.lower())
            if ahi is None or bool(c.upper() < ahi.upper()):
                ahi = arb(c.upper())
    return arb((blo / ahi).log().lower())


def main_f0(T):
    ctx.prec = PREC
    a = arb(A_NUM)
    t0 = time.time()
    out = cell_tables(T)
    if out is None:
        return False
    tab, h = out
    print(f"   DI.2: {NC} cells ({time.time() - t0:.1f}s)", flush=True)
    c0, d0 = float(tab[:, 2].min()), float(tab[:, 3].max())
    c1, d1 = float(tab[:, 4].min()), float(tab[:, 5].max())
    Rlo, Rhi = float(tab[:, 0].min()), float(tab[:, 1].max())
    logR = arb((arb(Rhi) / arb(Rlo)).log().upper())
    b0 = max(abs(c0), abs(d0))
    B0d = max(abs(c1), abs(d1))
    check(bool(arb(c0) > arb(Q0_RANGE[0])) and bool(arb(d0) < arb(Q0_RANGE[1]))
          and bool(arb(c1) > arb(Q1_RANGE[0])) and bool(arb(d1) < arb(Q1_RANGE[1])),
          f"DI.2: (log F0)' in [{c0:.8f}, {d0:.8f}], (log Lt F0)' in [{c1:.8f}, {d1:.8f}] on X_s "
          f"(DC.2 symmetric: |(log F0)'| <= {b0_DC}, Lt F0 in C_(B_0), B_0 = phi(b_0) <= {B0_DC}); direct sup |(log Lt F0)'| "
          f"<= {B0d:.8f}")
    check(bool(logR < arb(LOGR_MAX)) and bool(logR >= arb("0.0633864")) and b0 <= float(b0_DC) * (1 + 1e-9),
          f"DI.2: log(sup R_0/inf R_0) <= {logR.str(9)} < {LOGR_MAX} (centred-form enclosure of R_0; DC.2: {LOGR0_DC}, DC.4's "
          f"floor 0.0633864); consistency: sup|(log F0)'| <= {b0:.8f} <= DC.2's b_0 = {b0_DC}")
    # DI.3 signed factor lemma
    KK = signed_factor(a, arb(c0), arb(d0), arb(c1), arb(d1))
    check(KK is not None, "DI.3: DG.1(i)'s hypotheses hold in all four orientations with DI.2's signed slopes")
    Kup, Kdn = KK
    TS = arb((Kup.log() + Kdn.log() + logR).upper())
    # symmetric DH.1 with the direct B0 (for comparison)
    bu, bv = arb(b0), arb(B0d)
    TDsym = arb((((a + bv) / (a - bu)).log() + ((a + bu) / (a - bv)).log() + logR).upper())
    check(bool(TS < arb(TS_MAX)) and bool(TS < arb(TH_DH)),
          f"DI.3: Theta_100(f_0, Lt f_0) <= log K_up + log K_dn + log(sup R/inf R) = {TS.str(8)} < {TS_MAX} (cone part "
          f"{(TS - logR).str(6)}; DH.5: Theta_H' < {TH_DH}, cone part 0.0410521; DH.1 with symmetric b_0 and direct "
          f"sup|(log Lt F0)'|: {TDsym.str(8)})")
    # DI.4 coupled
    t0 = time.time()
    beta, alpha, jb, ja, upm, dnm = coupled(tab, h, a)
    TI = arb((arb(beta) / arb(alpha)).log().upper())
    print(f"   DI.4: cell pairs ({time.time() - t0:.1f}s); beta <= {beta:.12f} (cell {jb}), alpha >= {alpha:.12f} (cell {ja}), "
          f"sup R = {Rhi:.12f}, inf R = {Rlo:.12f}, excess at the extremal cells: up {upm:.3e}, down {dnm:.3e}", flush=True)
    fl = floor_points(T, a)
    check(bool(TI < arb(TI_MAX)) and bool(TI < TS) and bool(fl > arb(TFLOOR_MIN)) and bool(fl <= TI),
          f"DI.4: Theta_100(f_0, Lt f_0) <= Theta_I' = log(beta/alpha) = {TI.str(8)} < {TI_MAX} (DI.3: {TS.str(8)}; DH.5: "
          f"{TH_DH}); certified floor {fl.str(8)} > {TFLOOR_MIN}; Theta_I' - floor <= {(TI - fl).str(4)}; "
          f"Theta_I' - log(sup R/inf R) <= {(TI - logR).str(4)}")
    RESULT.update(c0=c0, d0=d0, c1=c1, d1=d1, logR=logR, TS=TS, TDsym=TDsym, TI=TI, floor=fl, beta=beta, alpha=alpha,
                  Kup=Kup, Kdn=Kdn, B0d=B0d)
    return True


# ----------------------------------------------------------------- DI.5 consequences
def consequences(T):
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    import gasket_fractional_profile_pure_point as cc
    TI, TH = RESULT["TI"], arb(TH_DH)
    thG, thH = arb(THETA_G), arb(THETA_H)
    r_lo, r_hi = ce.certificate()
    ctx.prec = PREC
    rl = arb(r_lo.numerator) / r_lo.denominator
    rh = arb(r_hi.numerator) / r_hi.denominator
    found, sums, counts, okb = cc.period_atoms(bk.KMAX, verbose=False)
    ctx.prec = PREC
    r = (rl + rh) / 2
    # consistency with DC.3: s_4, s_5 from F0 itself
    c0 = 3 * ce.sigma_closed(arb(-3), 0)
    x0 = arb(13).sqrt() / 2
    s4p, s5p = c0 * T.val(x0), c0 * bk.LtF(T, x0)
    ok_s = bool(abs(s4p - sums[4]) < arb("1e-12") * sums[4]) and bool(abs(s5p - sums[5]) < arb("1e-12") * sums[5])
    okd, worst = okb, arb(0)
    log4 = abs((sums[5] / (r * sums[4])).log())
    for k in range(4, bk.KMAX):
        dev = abs((sums[k + 1] / (r * sums[k])).log())
        bI = TI * arb(min((thG ** ((k - 4) // 8)).upper(), (thH ** ((k - 4) // 9)).upper()))
        okd = okd and bool(dev < bI)
        worst = arb(max(worst.upper(), (dev / bI).upper()))
    ratio = TH / TI
    check(ok_s and okd and bool(ratio > 1),
          f"DI.5: for every k >= 4, |log(s_(k+1)/(r s_k))| <= Theta_I' min(theta_G^floor((k-4)/8), theta_H^floor((k-4)/9)), "
          f"Theta_I' < {TI_MAX}: DH.5's bound divided by Theta_H'/Theta_I' >= {ratio.str(6)} at every k (same exponents). "
          f"Consistency: c0 F0(x0), c0 Lt F0(x0) = CC's s_4, s_5; on CC's s_4..s_{bk.KMAX} every deviation is below the bound "
          f"(max ratio {worst.str(4)}); at k = 4 the bound is Theta_I' = {(TI / log4).str(5)} |log(s_5/(r s_4))|")
    RESULT.update(ratio=ratio, worst=worst, log4=log4)


def main() -> int:
    t0 = time.time()
    print("Theorem DI: the f_0 term of DH.5 with signed slopes and a coupled Hilbert-quotient bound. Not thrust.", flush=True)
    T = taylor()
    if not main_f0(T):
        print("FAILED")
        return 1
    consequences(T)
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(bk.FAILS) + len(bk.co.FAILS) + len(bk.cn.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DI checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
