#!/usr/bin/env python3
"""Theorem BT: the log-periodic profile u^(1-gamma) Q_inf(u) of Theorem BS is AT.1's Theta composed with psi, is continuous,
has AU's mean, is NOT constant, and its maximum is the sharp constant sup_(s>0) s^(1-gamma) r_inf(s) of BI.4.

Setting (GASKET.md, Theorems AQ, AT, AU, BD, BI, BJ, BN, BS): p(s) = s(s+5), q = p^-1, gamma = log3/log5 (5^gamma = 3),
r_inf = lim r_k the corner resolvent, g(s) = s^(1-gamma) r_inf(s), Q_inf(u) = lim_M (3/5)^M r_inf(5^-M u) (BS.1),
F(u) = u^(1-gamma) Q_inf(u), psi(v) = lim p^k(5^-k v) (BD.4), Theta(s) = lim_k g(q^k s) (AT.1's profile).

Theorem BT (checked here; exit 1 on any failure):
  (1) [orbit monotonicity; every s > 0] with X = r_inf(p s) in AQ's fixed-point relation
          g(s)/g(p(s)) = (5+s)^gamma (2+s)/(6+s-2p(s)X)  in  [1, (5+s)^gamma (2+s)/(6+s-2G p(s)^gamma)],
      so g(q^k s) is NONDECREASING in k, Theta(s) = lim_k g(q^k s) = sup_k g(q^k s) exists for every s > 0 (bounded by G),
      Theta o q = Theta and g <= Theta <= g e^(E(s)), E(s) <= s/4 + (G/3) 6^gamma s^gamma for s <= 1e-3
      (AT.1's two-sided [-(1-gamma)y/5, 1.82 y^0.65] becomes one-sided).
  (2) [identification; no AT.1 input] F(u) = Theta(psi(u)) for every u > 0; F(5u) = F(u); for every v <= 1e-3, with
      Y = v e^(v/10):  e^(-gamma v/10) <= F(v)/g(v) <= e^((1-gamma) v/10 + Y/4 + (G/3) 6^gamma Y^gamma);
      F is continuous (uniform limit of g(5^-M .) on [1,5]).  So the hypothesis of BS's derived paragraph holds and
      c_inf = theta_inf^BH (BH.6 as stated).
  (3) [extremes and mean] sup_(s>0) g(s) = max F (BI.4's sharp constant), g(s) in [min F e^(-(1-gamma)s/10 - E), max F e^(gamma s/10)]
      for s <= 1e-3 (sharpens BJ.3), and the log-mean of F over one period equals AU's <Theta>.
  (4) [certified, python-flint arb] F on 2^17 log-equispaced nodes of [1,5) (M = 30), cell bounds from Q_inf nonincreasing and
      u Q_inf nondecreasing: max F, min F, <F>; consistency with AT.1's eight values, BJ.3's bracket, AU's <Theta>, BI.4's G.
  Kept Failure BT.5: "Q_inf(u) = C u^(gamma-1) for one constant C" (equivalently F constant; BS's floats use C = <Theta>) is
      false: max F - min F >= 1.46e-5, certified.
Assumptions used as stated (not re-run): AQ (3) (r_k increasing to r_inf, fixed-point recursion), BN.1 (r_k Stieltjes of unit mass:
r nonincreasing, s r(s) nondecreasing), BI.4 (sup g <= G = 0.6031), BJ.1 (1/r_k - 1/r_inf = prod_(i<k) beta(p^i s)(1+2 eta(p^k s))),
BH.2 (0 < eta <= 2/(2+z)), BD.4 (psi(5v) = p(psi(v)), psi(v) in [v, v e^(v/10)] with BJ.6), BS.1 (Q_inf exists, Q_inf(5u) =
(3/5)Q_inf(u), nonincreasing), AU (<Theta> = 0.6009542511..., consistency only), AT.1's eight values (consistency only).
Floats: consistency only.  Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

import sympy as sp
from flint import arb, ctx

FAILS: list[str] = []
G_BI4 = "0.6031"
G_LO, G_HI = "0.600945014", "0.600963487"           # BJ.3
THETA_AU = "0.6009542511"                             # AU's <Theta> (10 digits as recorded)
AT1_VALUES = ("0.600954332", "0.600949127", "0.600946924", "0.600949013",
              "0.600954170", "0.600959375", "0.600961578", "0.600959489")   # Theta(5^-40 5^(j/8)), AT.1, 9 digits
M_DEPTH = 30
NGRID = 2**17


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    s, X, gm, G, x = sp.symbols("s X gamma G x", positive=True)
    p = s * (s + 5)
    r_s = (2 + s) * (5 + s) * X / ((6 + s) - 2 * p * X)          # AQ's step: r(s) from X = r(p s)
    ratio = s**(1 - gm) * r_s / (p**(1 - gm) * X)
    target = (5 + s)**gm * (2 + s) / (6 + s - 2 * p * X)
    check(sp.simplify(sp.powsimp(sp.expand_power_base(ratio / target, force=True), force=True)) == 1,
          "g(s)/g(p s) = (5+s)^gamma (2+s)/(6+s-2p(s)X), X = r_inf(p s) (AQ fixed-point relation; s/p(s) = 1/(5+s))")
    check(sp.expand(3 * (2 + s) - (6 + s)) == 2 * s,
          "3(2+s) - (6+s) = 2s >= 0, and (5+s)^gamma >= 5^gamma = 3, 6+s-2pX = (2+s)(5+s)X/r_inf(s) > 0: g(s)/g(p s) >= 1")
    check(sp.simplify(sp.expand(((6 + s - 2 * p * X) * r_s) - (2 + s) * (5 + s) * X)) == 0,
          "denominator identity (6+s-2pX) r_inf(s) = (2+s)(5+s) X")
    # upper: (5+s)^gamma(2+s)/(6+s-2Gp^gamma) = (1+s/5)^gamma (1+s/2)/(1+s/6-Gp^gamma/3)
    check(sp.simplify(3 * (2 + s) / (6 + s - 2 * x) - (1 + s / 2) / (1 + s / 6 - x / 3)) == 0,
          "3(2+s)/(6+s-2x) = (1+s/2)/(1+s/6-x/3) (x = G p^gamma >= p X by BI.4): log ratio <= gamma s/5 + s/2 + 2x/3 "
          "for x/3 <= 1/2 (-log(1-y) <= 2y on [0,1/2])")
    y = sp.symbols("y", positive=True)
    check(sp.simplify(sp.diff(-sp.log(1 - y) - 2 * y, y) - (1 / (1 - y) - 2)) == 0 and
          sp.simplify((1 / (1 - y) - 2).subs(y, sp.Rational(1, 2))) == 0,
          "d/dy(-log(1-y) - 2y) = 1/(1-y) - 2 <= 0 on [0,1/2]: -log(1-y) <= 2y there")
    check(sp.expand(6 * s - p - s * (1 - s)) == 0, "6s - p(s) = s(1-s) >= 0: p(s) <= 6s for s <= 1")
    j = sp.symbols("j", integer=True, positive=True)
    check(sp.summation(sp.Rational(1, 5)**j, (j, 1, sp.oo)) == sp.Rational(1, 4) and
          sp.summation(sp.Rational(1, 3)**j, (j, 1, sp.oo)) == sp.Rational(1, 2),
          "sum 5^-j = 1/4, sum (5^-j)^gamma = sum 3^-j = 1/2: E(s) <= (gamma/5+1/2) s/4 + (2G/3) 6^gamma s^gamma/2 "
          "<= s/4 + (G/3) 6^gamma s^gamma (q^k s <= 5^-k s since p(y) >= 5y)")
    M, u = sp.symbols("M u", positive=True)
    lhs = sp.powsimp(sp.expand_power_base((5**(-M) * u)**(1 - gm), force=True), force=True)
    check(sp.simplify(lhs.subs(gm, sp.log(3) / sp.log(5)) - (sp.Rational(3, 5))**M * u**(1 - sp.log(3) / sp.log(5))) == 0,
          "(5^-M u)^(1-gamma) = (3/5)^M u^(1-gamma): g(5^-M u) = u^(1-gamma) (3/5)^M r_inf(5^-M u) -> F(u) (BS.1)")
    check(sp.simplify(sp.Integer(5)**(1 - sp.log(3) / sp.log(5)) * sp.Rational(3, 5) - 1) == 0,
          "5^(1-gamma) (3/5) = 1: F(5u) = F(u) from Q_inf(5u) = (3/5) Q_inf(u)")
    # cell bounds: max_t min(la+(1-g)t, lb+g(h-t)) <= max(la,lb) + g(1-g)h
    la, lb, h, t = sp.symbols("la lb h t", real=True)
    tc = sp.solve(sp.Eq(la + (1 - gm) * t, lb + gm * (h - t)), t)[0]
    val = sp.simplify(la + (1 - gm) * tc)
    check(sp.simplify(val - (gm * la + (1 - gm) * lb + gm * (1 - gm) * h)) == 0,
          "crossing of the two cell ceilings: value = gamma la + (1-gamma) lb + gamma(1-gamma) h <= max(la,lb) + gamma(1-gamma)h "
          "(and symmetrically for floors)")


# ----------------------------------------------------------------- (B) exact rationals
def _p(x):
    return x * (x + 5)


def r_rec(k, s):
    if k == 0:
        return 1 / (3 + s)
    x = r_rec(k - 1, _p(s))
    return (2 + s) * (5 + s) * x / ((6 + s) - 2 * s * (5 + s) * x)


def exact_checks():
    ok = True
    for k in range(1, 5):
        for s in (Fr(1, 1000), Fr(1, 7), Fr(1), Fr(9, 2), Fr(40)):
            X = r_rec(k - 1, _p(s))
            ok = ok and r_rec(k, s) * (6 + s - 2 * _p(s) * X) == (2 + s) * (5 + s) * X
            ok = ok and 6 + s - 2 * _p(s) * X > 0 and 0 < _p(s) * X <= 1
            ok = ok and 3 * (2 + s) >= 6 + s
    check(ok, "exact rationals k = 1..4, s in {1/1000, 1/7, 1, 9/2, 40}: r_k(s)(6+s-2p r_(k-1)(ps)) = (2+s)(5+s) r_(k-1)(ps), "
              "denominator > 0, 0 < p r_(k-1)(p) <= 1 (BN.1 unit mass), 3(2+s) >= 6+s")


# ----------------------------------------------------------------- (C) arb
def _gam():
    return arb(3).log() / arb(5).log()


def rinf_bracket(v, extra=6):
    """r_k(v) <= r_inf(v) <= 1/(1/r_k(v) - prod_(i<k) beta(p^i v)(1 + 4/(2+p^k v))) (AQ (3), BJ.1, BH.2's eta <= 2/(2+z))."""
    o = [v]
    while len([z for z in o if z > 1]) < extra:
        o.append(_p(o[-1]))
    k = len(o) - 1
    r = 1 / (3 + o[k])
    prod = arb(1)
    for i in range(k - 1, -1, -1):
        z = o[i]
        D = (2 + z) * (5 + z)
        r = D * r / ((6 + z) - 2 * z * (5 + z) * r)
        prod *= (z + 6) / D
    hi = 1 / (1 / r - prod * (1 + 4 / (2 + o[k])))
    return arb(r.lower()), arb(hi.upper())


def E_bound(Y, gam, G):
    return Y / 4 + G / 3 * arb(6) ** gam * Y ** gam


def F_bracket(u, gam, G, M=M_DEPTH):
    """BT.2 at v = 5^-M u: F(u) = F(v) in [g(v) e^(-gamma v/10), g(v) e^((1-gamma)v/10 + E(Y))]."""
    v = u / arb(5) ** M
    lo, hi = rinf_bracket(v)
    Y = v * (v / 10).exp()
    vg = v ** (1 - gam)
    L = vg * lo * (-gam * v / 10).exp()
    H = vg * hi * ((1 - gam) * v / 10 + E_bound(Y, gam, G)).exp()
    return arb(L.lower()), arb(H.upper())


def ratio_check(gam, G):
    """Consistency: enclosures of g(s)/g(p s) from r_inf brackets lie above 1 - width and below the BT.1 ceiling."""
    ok = True
    worst = 0.0
    for s in ("1e-20", "1e-8", "0.001", "0.1", "1", "7", "100"):
        s = arb(s)
        a_lo, a_hi = rinf_bracket(s, extra=8)
        b_lo, b_hi = rinf_bracket(_p(s), extra=8)
        ps = _p(s)
        rat_hi = (s / ps) ** (1 - gam) * a_hi / b_lo
        rat_lo = (s / ps) ** (1 - gam) * a_lo / b_hi
        ceil = (5 + s) ** gam * (2 + s) / (6 + s - 2 * G * ps ** gam) if bool(6 + s - 2 * G * ps ** gam > 0) else None
        ok = ok and bool(rat_hi >= 1)
        if ceil is not None:
            ok = ok and bool(rat_lo <= ceil)
        worst = max(worst, float((rat_hi - 1).upper()))
    return ok


def certify():
    ctx.prec = 96
    gam, G = _gam(), arb(G_BI4)
    h = arb(5).log() / NGRID
    vals = []
    for i in range(NGRID + 1):
        u = (h * i).exp()
        vals.append(F_bracket(u, gam, G))
    widths = max(float((b - a).upper()) for a, b in vals)
    # extremes over all of (0, inf): F is 5-periodic, so [1, 5] is a period; cells [u_i, u_(i+1)] of log-width h
    eh_up = (gam * (1 - gam) * h).exp()
    Fmax_lo = max(vals, key=lambda ab: float(ab[0].mid()))[0]
    imax = max(range(len(vals)), key=lambda i: float(vals[i][0].mid()))
    for (a, b) in vals:
        if bool(a > Fmax_lo):
            Fmax_lo = a
    cell_up = arb(0)
    cell_dn = None
    mean_hi, mean_lo = arb(0), arb(0)
    one_m = 1 - gam
    for i in range(NGRID):
        (al, ah), (bl, bh) = vals[i], vals[i + 1]
        up = (ah if bool(ah > bh) else bh) * eh_up
        if bool(up > cell_up):
            cell_up = up
        dn = (al if bool(al < bl) else bl) / eh_up
        if cell_dn is None or bool(dn < cell_dn):
            cell_dn = dn
        # mean: F(u_i e^t) <= F_i e^((1-gamma)t), <= F_(i+1) e^(gamma(h-t)); >= F_i e^(-gamma t), >= F_(i+1) e^(-(1-gamma)(h-t))
        u1 = ah * ((one_m * h).exp() - 1) / one_m
        u2 = bh * ((gam * h).exp() - 1) / gam
        l1 = al * (1 - (-gam * h).exp()) / gam
        l2 = bl * (1 - (-one_m * h).exp()) / one_m
        mean_hi += u1 if bool(u1 < u2) else u2
        mean_lo += l1 if bool(l1 > l2) else l2
    Fmin_hi = min(vals, key=lambda ab: float(ab[1].mid()))[1]
    imin = min(range(len(vals)), key=lambda i: float(vals[i][1].mid()))
    for (a, b) in vals:
        if bool(b < Fmin_hi):
            Fmin_hi = b
    ln5 = arb(5).log()
    return dict(vals=vals, width=widths, max=(Fmax_lo, arb(cell_up.upper())), min=(arb(cell_dn.lower()), Fmin_hi),
                mean=(arb((mean_lo / ln5).lower()), arb((mean_hi / ln5).upper())), imax=imax, imin=imin, h=h,
                ratio_ok=ratio_check(gam, G))


def _f(x, lo=True, d=10):
    return f"{float(x.lower() if lo else x.upper()):.{d}f}"


def floats_fourier(vals, kmax=3):
    """Floats only: first Fourier moduli of F in log_5 u from the certified midpoints (trapezoid)."""
    import cmath
    import math
    n = len(vals) - 1
    mids = [float(((a + b) / 2).mid()) for a, b in vals[:-1]]
    out = []
    for k in range(1, kmax + 1):
        c = sum(m * cmath.exp(-2j * math.pi * k * i / n) for i, m in enumerate(mids)) / n
        out.append(abs(c))
    return sum(mids) / n, out


def main() -> int:
    t0 = time.time()
    print("Theorem BT: u^(1-gamma) Q_inf(u) = Theta(psi(u)) is continuous, log-periodic with AU's mean, not constant, and its "
          "maximum is sup s^(1-gamma) r_inf. Not thrust.")
    symbolic_checks()
    exact_checks()
    res = certify()
    check(res["ratio_ok"], "arb consistency: enclosed g(s)/g(p s) >= 1 and <= the BT.1 ceiling at s = 1e-20, 1e-8, 1e-3, 0.1, 1, 7, "
                           "100 (r_inf brackets from BJ.1)")
    check(res["width"] < 1e-12, f"arb: F(u) enclosed at {NGRID + 1} nodes u = 5^(i/2^17), M = {M_DEPTH}, max width "
                                f"{res['width']:.2e}")
    gam = _gam()
    # AT.1 eight values: Theta(5^-40 5^(j/8)) = F(5^(j/8)) up to the psi shift (here computed at 5^-40 5^(j/8) directly)
    ok_at = True
    G = arb(G_BI4)
    at_vals = []
    for jj, rec in enumerate(AT1_VALUES):
        s = arb(5) ** (arb(jj) / 8 - 40)
        lo, hi = rinf_bracket(s)
        g_lo, g_hi = s ** (1 - gam) * lo, s ** (1 - gam) * hi
        th_hi = g_hi * E_bound(s * (s / 10).exp(), gam, G).exp()     # Theta in [g, g e^E] (BT.1)
        at_vals.append((g_lo, th_hi))
        ok_at = ok_at and bool(g_lo < arb(rec) + arb("5.1e-10")) and bool(th_hi > arb(rec) - arb("5.1e-10"))
    check(ok_at, "consistency: Theta(5^-40 5^(j/8)) in [g, g e^E] (BT.1) agrees with AT.1's eight recorded 9-digit values, "
                 "j = 0..7: " + ", ".join(_f(a, d=9) for a, _ in at_vals))
    (Mxl, Mxh), (Mnl, Mnh), (ml, mh) = res["max"], res["min"], res["mean"]
    print(f"   max F in [{_f(Mxl, d=12)}, {_f(Mxh, False, 12)}] (node u = 5^({res['imax']}/2^17)); min F in "
          f"[{_f(Mnl, d=12)}, {_f(Mnh, False, 12)}] (node u = 5^({res['imin']}/2^17)); <F> in [{_f(ml, d=12)}, {_f(mh, False, 12)}]")
    amp = Mxl - Mnh
    check(bool(amp > arb("1.46e-5")),
          f"Kept Failure BT.5: max F - min F >= {float(amp.lower()):.6e} > 0, so F is not constant: Q_inf is not C u^(gamma-1)")
    check(bool(Mnl > arb(G_LO)) and bool(Mxh < arb(G_HI)),
          f"consistency: [min F, max F] inside BJ.3's [{G_LO}, {G_HI}] (BT.3 sharpens it to width {float((Mxh - Mnl).upper()):.3e})")
    s0 = arb("1e-12")
    s0 = s0 * (s0 + 5)
    Y0 = s0 * (s0 / 10).exp()
    gb_lo = Mnl * (-(1 - gam) * s0 / 10 - E_bound(Y0, gam, G)).exp()
    gb_hi = Mxh * (gam * s0 / 10).exp()
    check(bool(gb_lo > arb("0.600945317")) and bool(gb_hi < arb("0.600963178")) and bool(gb_lo > arb(G_LO))
          and bool(gb_hi < arb(G_HI)),
          f"BT.3: for every s <= p(1e-12), g(s) in [{_f(gb_lo, d=10)}, {_f(gb_hi, False, 10)}] subset [0.600945317, 0.600963178], "
          f"inside BJ.3's bracket (factors monotone in s, worst case at s = p(1e-12))")
    check(bool(Mxh < arb(G_BI4)), f"BT.3: sup_(s>0) s^(1-gamma) r_inf(s) = max F <= {_f(Mxh, False, 9)} < BI.4's G = 0.6031")
    check(bool(ml < arb(THETA_AU) + arb("1e-10")) and bool(mh > arb(THETA_AU) - arb("1e-10")),
          f"BT.3 consistency: certified <F> in [{_f(ml, d=10)}, {_f(mh, False, 10)}] contains AU's <Theta> = {THETA_AU}")
    mean_f, four = floats_fourier(res["vals"])
    print(f"   floats: trapezoid <F> = {mean_f:.13f}; |F_k| (log_5 u Fourier moduli) k=1,2,3: "
          + ", ".join(f"{c:.3e}" for c in four))
    check(abs(mean_f - float(THETA_AU)) < 2e-10 and four[0] > 1e-6 and four[1] < 1e-9,
          "floats: trapezoid mean = AU's <Theta> to 2e-10; first harmonic ~3.7e-6 dominates (consistency only)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BT checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
