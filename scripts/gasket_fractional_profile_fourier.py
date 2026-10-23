#!/usr/bin/env python3
"""Theorem BU: BT's profile F(u) = u^(1-gamma) Q_inf(u) IS AU's theta (so F_k = theta_k for every Fourier mode k), and a
second-order (convexity) cell certificate pins max F, min F, <F> and the first two Fourier modes of F to ~1e-12.

Setting (GASKET.md, Theorems AQ, AT, AU, BD, BN, BT): p(s) = s(s+5), q = p^-1, gamma = log3/log5, r_inf = lim r_k,
g(s) = s^(1-gamma) r_inf(s), Theta = lim g o q^k (BT.1), F(u) = u^(1-gamma) Q_inf(u) = Theta(psi(u)) (BT.2),
psi = psi_BD = lim p^k(5^-k .) (BD.4), psi_AU = lim 5^k q^k (AU's Koenigs coordinate), theta = Theta o psi_AU^-1 (AU),
theta_k = (1/ln5) int_period theta(y) y^(-i omega k) dy/y = a_k - i b_k, omega = 2 pi/ln5.

Theorem BU (checked here; exit 1 on any failure):
  (1) [identification] psi_AU o psi_BD = id on (0, inf) and psi_AU is strictly increasing, so psi_AU^-1 = psi_BD and
      F = theta identically.  Hence F_k = theta_k for EVERY k in Z (BT's open item), and AU's certified table gives
      F_1, ..., F_4 != 0; AU (2)'s Mellin residues of r_inf are exactly the Fourier coefficients of BS's profile.
  (2) [convex cell sandwich] r_inf is convex and decreasing (limit of the Stieltjes functions r_k, BN.1).  On a cell
      [v_i, v_(i+1)] the chord through the node ceilings is a ceiling of r_inf; the two neighbouring chords through
      (node floor, outer node ceiling), extended into the cell, are floors.  For a line c + d v with c > 0 > d,
      phi(v) = v^(1-gamma)(c + d v) is strictly concave, so its extremes over a cell are at the ends, at the crossing of
      the two floor lines, or at the stationary point v* = -(1-gamma)c/((2-gamma)d); its log-integral and its Mellin
      moments are closed forms.  Errors are O(h^2) instead of BT.4's O(h).
  (3) [certified, python-flint arb; 2^18 log-equispaced cells of one period, M = 30]
      max F, min F, <F>, a_1, b_1, a_2 enclosed to ~1e-12; consistency with AU's 34-digit <Theta> and a_1, b_1, a_2.
  Kept Failure BU.4: "F(t + 1/2) = 2<F> - F(t)" (t = log_5 u; the half-period antisymmetry suggested by BT.4's extremal
      nodes 98534 - 32998 = 2^16, and satisfied by every odd-harmonic profile) is false: (max F + min F)/2 - <F> >= 2.4e-11.
Assumptions used as stated (not re-run): AQ (3), BN.1 (r_k Stieltjes of unit mass), BJ.1, BH.2, BD.4 with BJ.6
(psi(5v) = p(psi(v)), v <= psi(v) <= v e^(v/10)), AU's setting (psi_AU = s/A_inf(s), psi_AU o q = psi_AU/5,
0 <= s (log A_inf)'(s) <= s/20), BT.1-BT.2 (incl. the transfer bracket F(v)/g(v)), BI.4 (G = 0.6031); AU's table for consistency.
Floats: consistency only.  Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import sympy as sp
from flint import acb, arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasket_fractional_log_periodic_profile as bt  # noqa: E402  (BT's r_inf bracket and E-bound, unchanged)

FAILS: list[str] = []
NGRID = 2**18
M_DEPTH = 30
KMAX = 4
# AU's certified values (GASKET Theorem AU (3)), consistency only
AU_MEAN = "0.6009542511687736078237051500999942"
AU_A = {1: "4.042278096762267394377e-8", 2: "-1.347221801190789742277e-11"}
AU_B = {1: "-3.663477703263059920737e-6", 2: "6.374978894928635265275e-13"}
BT_MAX = ("0.6009615785", "0.6009631774")
BT_MIN = ("0.6009453250", "0.6009469238")


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    v, c, a, lam, s = sp.symbols("v c a lambda s", positive=True)
    d = sp.symbols("d", negative=True)
    phi = v**a * (c + d * v)
    dd = sp.simplify(sp.diff(phi, v, 2) - a * v**(a - 2) * ((a - 1) * c + (a + 1) * d * v))
    check(dd == 0, "phi = v^a (c + d v): phi'' = a v^(a-2)((a-1)c + (a+1)d v) < 0 for 0 < a < 1, c > 0 > d: phi strictly concave")
    vs = sp.solve(sp.Eq(sp.diff(phi, v), 0), v)
    check(len(vs) == 1 and sp.simplify(vs[0] + a * c / ((a + 1) * d)) == 0,
          "stationary point of phi: v* = -a c/((a+1) d) (the only critical point; concave, so it is the maximum)")
    z = sp.symbols("z")
    Phi = c * v**z / z + d * v**(z + 1) / (z + 1)
    check(sp.simplify(sp.diff(Phi, v) - v**(z - 1) * (c + d * v)) == 0,
          "int v^(z-1)(c + d v) dv = c v^z/z + d v^(z+1)/(z+1) for every complex z (z = a gives <F>, z = a - i omega k gives F_k)")
    check(sp.simplify(sp.diff(1 / (lam + s), s, 2) - 2 / (lam + s)**3) == 0 and
          sp.simplify(sp.diff(1 / (lam + s), s) + 1 / (lam + s)**2) == 0,
          "d/ds 1/(lam+s) = -1/(lam+s)^2 < 0, d^2/ds^2 = 2/(lam+s)^3 > 0: Stieltjes r_k (BN.1) are decreasing and convex, and so "
          "is their pointwise limit r_inf (AQ (3))")
    x0, x1, x2, r0, r1, r2 = sp.symbols("x0 x1 x2 r0 r1 r2", real=True)
    # convexity: for x >= x1 the line through (x0,r0),(x1,r1) lies below r; floor with interval values
    L = r1 + (r1 - r0) * (v - x1) / (x1 - x0)
    lam_ = (v - x1) / (x1 - x0)
    check(sp.simplify(L - (r1 * (1 + lam_) - r0 * lam_)) == 0,
          "extended chord = r_1(1+lam) - r_0 lam, lam = (v-v_1)/(v_1-v_0) >= 0: monotone up in r_1, down in r_0, so (floor r_1, "
          "ceiling r_0) gives a floor of the true extended chord, which is <= r_inf by convexity")
    j = sp.symbols("j", integer=True, positive=True)
    check(sp.summation(sp.Rational(1, 5)**j / 5, (j, 1, sp.oo)) == sp.Rational(1, 20),
          "sum_(i>=1) 5^-i/5 = 1/20: s (log A_inf)'(s) = s sum (q^i)'/(5+q^i) <= s/20 (AU), so d log psi_AU/d log s >= 1 - s/20")
    k, Mm = sp.symbols("k M", integer=True)
    check(sp.simplify(sp.exp(-2 * sp.pi * sp.I * k * Mm)) == 1,
          "u^(-i omega k) = v^(-i omega k) e^(-2 pi i k M) = v^(-i omega k) for u = 5^M v (omega = 2 pi/ln5): Fourier modes "
          "may be computed at depth v = 5^-M u")
    t = sp.symbols("t", real=True)
    f = sp.Function("f")
    # antisymmetry F(t+1/2) = 2m - F(t) kills even modes: int_0^1 F e^{-2 pi i 2l t} = 0
    ll = sp.symbols("l", integer=True)
    check(sp.simplify(sp.exp(-2 * sp.pi * sp.I * 2 * ll * (t + sp.Rational(1, 2))) - sp.exp(-2 * sp.pi * sp.I * 2 * ll * t)) == 0,
          "even modes are half-period invariant: F(t+1/2) = 2m - F(t) forces F_(2l) = -F_(2l), i.e. F_(2l) = 0 (l != 0), and "
          "max F + min F = 2m")


# ----------------------------------------------------------------- (B) exact rationals
def exact_checks():
    ok = True
    pts = [Fr(1, 100) * k for k in range(1, 12)] + [Fr(2), Fr(5), Fr(17, 2), Fr(40)]
    for kk in range(0, 5):
        vals = [bt.r_rec(kk, x) for x in pts]
        for i in range(len(pts) - 2):
            x0, x1, x2 = pts[i], pts[i + 1], pts[i + 2]
            s01 = (vals[i + 1] - vals[i]) / (x1 - x0)
            s12 = (vals[i + 2] - vals[i + 1]) / (x2 - x1)
            ok = ok and s01 < 0 and s12 < 0 and s12 > s01          # decreasing, convex (slopes increase)
    check(ok, "exact rationals k = 0..4, 15 points in [1/100, 40]: r_k strictly decreasing with increasing chord slopes "
              "(convexity of the Stieltjes r_k, BN.1)")
    # psi_BD o psi_AU structure at exact rationals: p^n(5^-n tau) increases in n (BD.4's psi_n increasing), tau = 1, 3, 7
    ok2 = True
    for tau in (Fr(1), Fr(3), Fr(7)):
        prev = None
        for n in range(0, 6):
            x = tau / 5**n
            for _ in range(n):
                x = x * (x + 5)
            ok2 = ok2 and (prev is None or x > prev) and x >= tau
            prev = x
    check(ok2, "exact rationals: psi_n(tau) = p^n(5^-n tau) >= tau and strictly increasing in n = 0..5 for tau = 1, 3, 7 (BD.4)")


# ----------------------------------------------------------------- (C) arb
def _psi_bd(v, n=60):
    x = v / arb(5) ** n
    for _ in range(n):
        x = x * (x + 5)
    return x                                      # psi_n(v) <= psi(v) <= psi_n(v) e^(5^-n v/10) (BD.4/BJ.6 bootstrap)


def _q(x):
    return 2 * x / (5 + (25 + 4 * x).sqrt())     # = (-5 + sqrt(25+4x))/2 without cancellation


def identification_check():
    """Consistency of BU.1's step (not an independent check): 5^40 q^40(psi_60(v)) = 5^40 psi_20(5^-40 v) -> v."""
    worst = 0.0
    for v in ("0.3", "1", "2.5", "4"):
        v = arb(v)
        y = _psi_bd(v)
        for _ in range(40):
            y = _q(y)
        y = y * arb(5) ** 40
        worst = max(worst, float(abs(y - v).upper()))
    return worst


def certify(N=NGRID):
    ctx.prec = 96
    gam = bt._gam()
    a = 1 - gam
    G = arb(bt.G_BI4)
    ln5 = arb(5).log()
    h = ln5 / N
    scale = arb(5) ** (-M_DEPTH)
    omega = 2 * arb.pi() / ln5
    V, R = {}, {}
    for i in range(-1, N + 2):
        v = scale * (h * i).exp()
        V[i] = v
        R[i] = bt.rinf_bracket(v)

    def phi(v, c, d):
        return v ** a * (c + d * v)

    def I0(x, y, c, d):
        return c * (y ** a - x ** a) / a + d * (y ** (a + 1) - x ** (a + 1)) / (a + 1)

    ez = [acb(a, -omega * k) for k in range(KMAX + 1)]
    ez1 = [acb(a + 1, -omega * k) for k in range(KMAX + 1)]

    def Ik(x, y, c, d, k):
        X, Y = acb(x), acb(y)
        return c * (Y ** ez[k] - X ** ez[k]) / ez[k] + d * (Y ** ez1[k] - X ** ez1[k]) / ez1[k]

    signs_ok = True
    below_G = True
    node_max_lo = node_min_hi = None
    cell_max_hi = cell_min_lo = None
    I_hi, I_lo, gap = arb(0), arb(0), arb(0)
    Ck = [acb(0) for _ in range(KMAX + 1)]
    for i in range(N):
        vi, vj = V[i], V[i + 1]
        # node values of F (BT.2 transfer at v_i)
        Yi = vi * (vi / 10).exp()
        nlo = vi ** a * R[i][0] * (-gam * vi / 10).exp()
        nhi = vi ** a * R[i][1] * ((1 - gam) * vi / 10 + bt.E_bound(Yi, gam, G)).exp()
        if node_max_lo is None or bool(nlo.lower() > node_max_lo.lower()):
            node_max_lo, imax = arb(nlo.lower()), i
        if node_min_hi is None or bool(nhi.upper() < node_min_hi.upper()):
            node_min_hi, imin = arb(nhi.upper()), i
        # ceiling chord and the two extended floor chords
        d = (R[i + 1][1] - R[i][1]) / (vj - vi)
        c = R[i][1] - d * vi
        d1 = (R[i][0] - R[i - 1][1]) / (vi - V[i - 1])
        c1 = R[i][0] - d1 * vi
        d2 = (R[i + 2][1] - R[i + 1][0]) / (V[i + 2] - vj)
        c2 = R[i + 1][0] - d2 * vj
        signs_ok = signs_ok and all(bool(x > 0) for x in (c, c1, c2)) and all(bool(x < 0) for x in (d, d1, d2))
        cands = [phi(vi, c, d), phi(vj, c, d)]
        vs = -a * c / ((a + 1) * d)
        if not (bool(vs < vi) or bool(vs > vj)):
            cands.append(phi(vs, c, d))
        up = max(cands, key=lambda x: float(x.upper()))
        below_G = below_G and bool(up < G)
        vc = arb(((c2 - c1) / (d1 - d2)).mid())
        if not (bool(vc > vi) and bool(vc < vj)):
            vc = arb(((vi + vj) / 2).mid())
        lo = min([phi(vi, c1, d1), phi(vc, c1, d1), phi(vc, c2, d2), phi(vj, c2, d2)], key=lambda x: float(x.lower()))
        Yj = vj * (vj / 10).exp()
        fu = ((1 - gam) * vj / 10 + bt.E_bound(Yj, gam, G)).exp() - 1
        fd = 1 - (-gam * vj / 10).exp()
        Fhi, Flo = up * (1 + fu), lo * (1 - fd)
        if cell_max_hi is None or bool(Fhi.upper() > cell_max_hi.upper()):
            cell_max_hi = arb(Fhi.upper())
        if cell_min_lo is None or bool(Flo.lower() < cell_min_lo.lower()):
            cell_min_lo = arb(Flo.lower())
        Ih = I0(vi, vj, c, d)
        Il = I0(vi, vc, c1, d1) + I0(vc, vj, c2, d2)
        I_hi += Ih * (1 + fu)
        I_lo += Il * (1 - fd)
        gap += (Ih - Il) / 2 + G * (fu + fd) * h
        for k in range(1, KMAX + 1):
            Ck[k] += (Ik(vi, vj, c, d, k) + Ik(vi, vc, c1, d1, k) + Ik(vc, vj, c2, d2, k)) / 2
    err = gap / ln5
    modes = {}
    for k in range(1, KMAX + 1):
        th = Ck[k] / ln5
        ak, bk = th.real, -th.imag
        modes[k] = ((arb((ak - err).lower()), arb((ak + err).upper())), (arb((bk - err).lower()), arb((bk + err).upper())))
    return dict(max=(node_max_lo, cell_max_hi), min=(cell_min_lo, node_min_hi), imax=imax, imin=imin,
                mean=(arb((I_lo / ln5).lower()), arb((I_hi / ln5).upper())), err=err, modes=modes,
                signs_ok=signs_ok, below_G=below_G)


def _f(x, lo=True, d=13):
    return f"{float(x.lower() if lo else x.upper()):.{d}f}"


def _in(lohi, val):
    return bool(lohi[0] <= arb(val)) and bool(arb(val) <= lohi[1])


def main() -> int:
    t0 = time.time()
    print("Theorem BU: F = u^(1-gamma) Q_inf is AU's theta (F_k = theta_k for every k); convexity cells pin max F, min F, <F>, "
          "F_1, F_2 to ~1e-12. Not thrust.")
    symbolic_checks()
    exact_checks()
    ctx.prec = 96
    w = identification_check()
    check(w < 1e-12, f"arb consistency of BU.1's step: |5^40 q^40(psi_60(v)) - v| <= {w:.1e} at v = 0.3, 1, 2.5, 4 (psi_AU depth 40, psi_BD depth 60)")
    res = certify()
    check(res["signs_ok"], f"BU.2 hypotheses on all {NGRID} cells: every ceiling/floor line has c > 0 > d (concavity of phi applies)")
    check(res["below_G"], "every cell ceiling of g is below BI.4's G = 0.6031 (used in the transfer slack)")
    (Mxl, Mxh), (Mnl, Mnh), (ml, mh) = res["max"], res["min"], res["mean"]
    print(f"   max F in [{_f(Mxl)}, {_f(Mxh, False)}] (node 5^({res['imax']}/2^18)); min F in [{_f(Mnl)}, {_f(Mnh, False)}] "
          f"(node 5^({res['imin']}/2^18)); <F> in [{_f(ml, d=14)}, {_f(mh, False, 14)}]")
    check(bool(Mxh - Mxl < arb("2e-12")) and bool(Mnh - Mnl < arb("5e-12")) and bool(mh - ml < arb("3e-12")),
          f"BU.3 widths: max {float((Mxh - Mxl).upper()):.2e}, min {float((Mnh - Mnl).upper()):.2e}, mean {float((mh - ml).upper()):.2e} "
          f"(BT.4: 1.6e-6, 1.6e-6, 2.3e-6)")
    check(bool(Mxl >= arb(BT_MAX[0])) and bool(Mxh <= arb(BT_MAX[1])) and bool(Mnl >= arb(BT_MIN[0])) and bool(Mnh <= arb(BT_MIN[1])),
          "BU.3 brackets lie inside BT.4's max/min brackets")
    check(_in(res["mean"], AU_MEAN), f"consistency: <F> bracket contains AU's 34-digit <Theta> = {AU_MEAN[:16]}...")
    amp = (Mxl - Mnh, Mxh - Mnl)
    print(f"   max F - min F in [{float(amp[0].lower()):.10e}, {float(amp[1].upper()):.10e}]")
    check(bool(amp[0] > arb("1.4654e-5")), "Kept Failure BT.5 sharpened: max F - min F >= 1.4654e-5")
    err = res["err"]
    print(f"   Fourier error radius (all k): {float(err.upper()):.3e}")
    for k in range(1, KMAX + 1):
        (al, ah), (bl, bh) = res["modes"][k]
        print(f"   k={k}: a_k in [{float(al.lower()):.6e}, {float(ah.upper()):.6e}], b_k in [{float(bl.lower()):.6e}, {float(bh.upper()):.6e}]")
    (a1, b1), (a2, b2) = res["modes"][1], res["modes"][2]
    check(_in(a1, AU_A[1]) and _in(b1, AU_B[1]) and _in(a2, AU_A[2]) and _in(b2, AU_B[2]),
          "BU.1 numerically: the independently certified a_1, b_1, a_2, b_2 of F contain AU's certified a_1, b_1, a_2, b_2 of theta")
    check(bool(b1[1] < 0) and bool(a1[0] > 0) and bool(a2[1] < 0),
          "BU.3: F_1 != 0 (a_1 > 0 > b_1) and F_2 != 0 (a_2 < 0) certified from F's own nodes, independently of AU's computation")
    half = ((Mxl + Mnl) / 2 - mh, (Mxh + Mnh) / 2 - ml)
    print(f"   (max F + min F)/2 - <F> in [{float(half[0].lower()):.4e}, {float(half[1].upper()):.4e}]")
    check(bool(half[0] > arb("2.4e-11")),
          f"Kept Failure BU.4: (max F + min F)/2 - <F> >= {float(half[0].lower()):.4e} > 0, so F(t+1/2) = 2<F> - F(t) is false")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BU checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
