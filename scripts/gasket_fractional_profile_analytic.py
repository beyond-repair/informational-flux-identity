#!/usr/bin/env python3
r"""Theorem BV: the profile F = theta of Theorems BT/BU is holomorphic and bounded on the strip |Im t| < pi/ln5
(t = log_5 u), so its Fourier modes decay at least like e^(-2 pi^2 |k|/ln5); theta_5..theta_11 are certified nonzero;
and F has exactly one maximum and one minimum per period, located to ~1e-55 and evaluated to ~1e-60 (Kept Failure BV.5: the
extremes are NOT exactly half a period apart).

Setting (GASKET.md, Theorems AQ, AU, BN, BS, BT, BU): r_k(s) = (1/2) v^T (L_k + s)^-1 v, r_inf = lim r_k (AQ (3)),
Q_inf(u) = lim_M (3/5)^M r_inf(5^-M u) (BS.1), gamma = log3/log5, F(u) = u^(1-gamma) Q_inf(u) = theta(u) (BT.2, BU.1),
F(5u) = F(u); theta_k = int_0^1 F(5^t) e^(-2 pi i k t) dt = a_k - i b_k (AU's normalisation), so
F = <F> + 2 sum_(k>=1) (a_k cos 2 pi k t + b_k sin 2 pi k t).  G* = max F <= 0.6009615785983 (BU.3).

Theorem BV (checked here; exit 1 on any failure):
  (1) [sector bound] every Stieltjes f(z) = int dmu(lam)/(lam+z), mu >= 0 on [0,inf), satisfies
      |f(z)| <= f(|z|)/cos(phi/2) for z = |z| e^(i phi), |phi| < pi, since
      |lam+z|^2 - (lam+|z|)^2 cos^2(phi/2) = (1-cos phi)(lam-|z|)^2/2 >= 0 (sympy).  BN.1: each r_k is Stieltjes.
  (2) [continuation] By Vitali (local bound from (1) and monotonicity in |z|), r_inf and Q_inf extend holomorphically
      to C \ (-inf,0] with the same bound; F(t) = 5^((1-gamma)t) Q_inf(5^t) is holomorphic and 1-periodic on
      S = {|Im t| < pi/ln5}, and |F(x+iy)| <= F(x)/cos(y ln5/2) <= G*/cos(y ln5/2).
  (3) [decay] shifting the contour to Im t = -y: |theta_k| <= G* e^(-2 pi |k| y)/cos(y ln5/2) for every y in (0, pi/ln5);
      optimal tan(y ln5/2) = T_k = 4 pi |k|/ln5 gives
         |theta_k| <= B_k = G* sqrt(1+T_k^2) exp(-T_k arctan T_k) <= e G* (1+T_k) e^(-2 pi^2 |k|/ln5),
      so limsup |theta_k|^(1/|k|) <= e^(-2 pi^2/ln5) = 4.70...e-6 (BU's open "decay rate", upper half).
  (4) [certified, AU's own routine re-run: mpmath.iv outward, 170 digits, 384 cells, sandwich b = u = 20/21, depth
      5^-140] theta_1..theta_11 to +-1.7e-62: theta_5, ..., theta_11 != 0 (AU: k <= 4), and |theta_k| <= B_k for every
      k <= 11 (a check of (3) against an independent certificate).
  (5) [derived from (3)-(4), python-flint arb] |F - S_11| <= e_0 with S_11 the certified 11-mode sum and e_0 from
      (3)'s tail; same for F', F''.  On a 4096-cell partition of one period every cell has F' != 0 or F'' != 0, F' changes
      sign exactly twice, so F has exactly one maximum t_max and one minimum t_min per period (both nondegenerate);
      bisection pins both locations to ~1e-55 and max F = G*, min F to ~1e-60.
  Kept Failure BV.5: "t_min = t_max + 1/2 (mod 1)" (BU.4's extremal nodes 197068 - 65996 = 2^17 exactly) is false:
      t_min - t_max - 1/2 is certified nonzero (value below).
Assumptions used as stated (not re-run): AQ (3) (r_k increasing to r_inf on (0,inf)), BN.1 (r_k Stieltjes, unit mass),
BS.1 (Q_inf exists), BT.2 (F(5u) = F(u)), BU.1 (F = theta), BU.3 (G* <= 0.6009615785983); AU (3)'s certified routine
(AQ (4) node enclosures, AS (2) sandwich, AT.1 orbit lemma, AU's psi corrections), imported unchanged.
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

FAILS: list[str] = []
KMAX = 11              # modes certified by AU's routine
K_PER, ORD, DPS, CAP_EXP = 384, 20, 170, 130
NCELL = 4096           # partition of one period for the critical-point count
G_STAR_HI = "0.6009615785983"   # BU.3 ceiling of max F
BU_MAX = ("0.6009615785971", "0.6009615785983")
BU_MIN = ("0.6009469237907", "0.6009469237943")
BU_NODES = (197068, 65996)      # BU.3's extremal nodes u = 5^(i/2^18)
AU_TABLE_A = {1: "4.042278096762267394377e-8", 2: "-1.347221801190789742277e-11",
              3: "-1.664604061435855694403e-17", 4: "3.253673674396170403376e-22"}
AU_TABLE_B = {1: "-3.663477703263059920737e-6", 2: "6.374978894928635265275e-13",
              3: "8.666201527432484696846e-17", 4: "2.762525729016318114163e-22"}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    lam, R = sp.symbols("lambda R", nonnegative=True)
    ph = sp.symbols("phi", real=True)
    lhs = (lam + R * sp.cos(ph))**2 + (R * sp.sin(ph))**2 - (lam + R)**2 * (1 + sp.cos(ph)) / 2
    check(sp.simplify(sp.expand(lhs - (1 - sp.cos(ph)) * (lam - R)**2 / 2)) == 0,
          "BV.1: |lam + R e^(i phi)|^2 - (lam+R)^2 cos^2(phi/2) = (1 - cos phi)(lam - R)^2/2 >= 0, so |1/(lam+z)| <= "
          "1/((lam+|z|) cos(phi/2)) and |f(z)| <= f(|z|)/cos(phi/2) for every Stieltjes f (a constant term obeys it too)")
    x, y, g, T = sp.symbols("x y gamma T", real=True)
    check(sp.simplify(sp.Abs(sp.exp((1 - g) * (x + sp.I * y) * sp.log(5))) - sp.exp((1 - g) * x * sp.log(5))) == 0,
          "BV.2: |5^((1-gamma)(x+iy))| = 5^((1-gamma)x) and arg 5^(x+iy) = y ln5, so |F(x+iy)| <= F(x)/cos(y ln5/2) on |y| < pi/ln5")
    k = sp.symbols("k", positive=True)
    check(sp.simplify(sp.Abs(sp.exp(-2 * sp.pi * sp.I * k * (x - sp.I * y))) - sp.exp(-2 * sp.pi * k * y)) == 0,
          "BV.3: |e^(-2 pi i k (x - i y))| = e^(-2 pi k y); the vertical sides of the rectangle cancel by 1-periodicity")
    f_ = sp.symbols("f", real=True)          # f = y ln5/2 in (0, pi/2); exponent -2 pi k y = -T f with T = 4 pi k/ln5
    h = -T * f_ - sp.log(sp.cos(f_))
    check(sp.simplify(sp.diff(h, f_) - (-T + sp.tan(f_))) == 0 and sp.simplify(sp.diff(h, f_, 2) - 1 / sp.cos(f_)**2) == 0,
          "BV.3: log of e^(-T f)/cos f is convex in f with stationary point tan f = T, so the optimum is "
          "B_k = G* sqrt(1+T^2) e^(-T arctan T)")
    Tp = sp.symbols("T", positive=True)
    check(sp.simplify(sp.diff(sp.atan(Tp) + sp.atan(1 / Tp), Tp)) == 0 and sp.limit(sp.atan(Tp) + sp.atan(1 / Tp), Tp, 1) == sp.pi / 2
          and sp.simplify(sp.diff(Tp - sp.atan(Tp), Tp) - Tp**2 / (1 + Tp**2)) == 0,
          "arctan T = pi/2 - arctan(1/T) (T > 0) and arctan(1/T) <= 1/T (x - arctan x has derivative x^2/(1+x^2) >= 0) give "
          "T arctan T >= (pi/2)T - 1 = 2 pi^2 k/ln5 - 1; with sqrt(1+T^2) <= 1+T: B_k <= e G* (1+T) e^(-2 pi^2 k/ln5)")
    xx, aa, KK = sp.symbols("x a K", positive=True)
    j = sp.symbols("j", integer=True, nonnegative=True)
    closed = lambda K_: xx**K_ * ((1 + aa * K_) / (1 - xx) + aa * xx / (1 - xx)**2)
    check(sp.simplify(closed(KK) - closed(KK + 1) - (1 + aa * KK) * xx**KK) == 0,
          "tail sum (telescoping; closed form -> 0 as K -> inf for 0 < x < 1): sum_(k>=K) (1 + a k) x^k = "
          "x^K((1+aK)/(1-x) + a x/(1-x)^2)")


# ----------------------------------------------------------------- (B) exact rationals (Gaussian rationals)
class Q2:
    """Exact Gaussian rationals (re, im) for the BN.1 sector check."""

    def __init__(self, re, im=Fr(0)):
        self.re, self.im = Fr(re), Fr(im)

    def __add__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(self.re + o.re, self.im + o.im)

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return Q2(o) - self

    def __mul__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = o if isinstance(o, Q2) else Q2(o)
        n = o.re * o.re + o.im * o.im
        return Q2((self.re * o.re + self.im * o.im) / n, (self.im * o.re - self.re * o.im) / n)

    def __rtruediv__(self, o):
        return Q2(o) / self

    def abs2(self):
        return self.re * self.re + self.im * self.im


def r_rec(n, s):
    """AQ's decimation recursion for r_n (as in AU's r_mp), exact for rationals or Gaussian rationals."""
    orb = [s]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    x = 1 / (3 + orb[-1])
    for kk in range(n - 1, -1, -1):
        o = orb[kk]
        x = (2 + o) * (5 + o) * x / ((6 + o) - 2 * o * (5 + o) * x)
    return x


def exact_checks():
    ok, worst = True, Fr(0)
    for n in range(0, 4):
        for Rv in (Fr(1, 50), Fr(1, 3), Fr(2), Fr(7), Fr(30)):
            for m in (Fr(1, 4), Fr(1), Fr(3), Fr(10), Fr(40)):
                # cos phi = (1-m^2)/(1+m^2), sin phi = 2m/(1+m^2), cos^2(phi/2) = 1/(1+m^2)
                z = Q2(Rv * (1 - m * m) / (1 + m * m), Rv * 2 * m / (1 + m * m))
                lhs = r_rec(n, z).abs2() / (1 + m * m)
                rhs = r_rec(n, Rv) ** 2
                ok = ok and lhs <= rhs
                worst = max(worst, lhs / rhs)
    check(ok, f"exact Gaussian rationals, r_0..r_3 at 125 points z = R e^(i phi) (|phi| up to 2 arctan 40): "
              f"|r_n(z)|^2 cos^2(phi/2) <= r_n(|z|)^2 (BV.1 for the actual resolvents; max ratio {float(worst):.4f})")


# ----------------------------------------------------------------- (C) AU's routine, re-run with more cells
def _mpf_to_arb(x):
    s, man, e, _ = x._mpf_
    if man == 0:
        return arb(0)
    v = arb(int(man)) * arb(2) ** int(e)          # exact value enclosed (arb is rigorous)
    return -v if s else v


def _iv_to_arb(z):
    from mpmath import mpf
    lo, hi = _mpf_to_arb(mpf(z.a)), _mpf_to_arb(mpf(z.b))
    return lo.union(hi)


def au_modes():
    import gasket_fractional_log_mean as au
    from mpmath import iv, mp, mpf
    iv.dps, mp.dps = DPS, DPS
    au.K_PER, au.ORD, au.CAP = K_PER, ORD, mpf(10) ** CAP_EXP
    P065 = iv.mpf(13) / 20
    K065 = 1 / (1 - iv.mpf(5) ** (-P065))
    gam = iv.log(3) / iv.log(5)
    data = au.period_integrals(au.E_DEPTH, iv.mpf(0), gam, KMAX)
    nodes_ok = all(v.a > 0 and mpf(v.b) - mpf(v.a) < mpf(10) ** -100 * mpf(v.a) for v in data["vals"].values())
    mean, coeffs, info = au.corrections(au.E_DEPTH, iv.mpf(0), gam, P065, K065, KMAX, data)
    out = {"mean": _iv_to_arb(mean), "a": {}, "b": {}, "nodes_ok": nodes_ok}
    for kk in range(1, KMAX + 1):
        out["a"][kk] = _iv_to_arb(coeffs[kk][0])
        out["b"][kk] = _iv_to_arb(coeffs[kk][1])
    return out


# ----------------------------------------------------------------- (D) arb: decay bound and extremes
def B_bound(k, G):
    T = 4 * arb.pi() * k / arb(5).log()
    return G * (1 + T * T).sqrt() * (-T * T.atan()).exp()


def tail_bound(K, j, G):
    """2 sum_(k>=K) (2 pi k)^j B_k: explicit terms K..K+299, then a geometric remainder."""
    tot = arb(0)
    for k in range(K, K + 300):
        tot += (2 * arb.pi() * k) ** j * B_bound(k, G)
    K2 = K + 300
    x = (-2 * arb.pi() ** 2 / arb(5).log()).exp()
    a = 4 * arb.pi() / arb(5).log()
    ratio = (arb(K2 + 1) / K2) ** (j + 1) * x          # >= ratio of consecutive k^j (1+ak) x^k for k >= K2
    rem = arb(1).exp() * G * (2 * arb.pi()) ** j * arb(K2) ** j * (1 + a * K2) * x ** K2 / (1 - ratio)
    return 2 * (tot + rem)


def S_der(t, j, coef):
    """j-th derivative of the certified 11-mode sum at the arb ball t (the mean enters only for j = 0)."""
    tot = acb(0)
    for kk in range(1, KMAX + 1):
        th = acb(coef["a"][kk], -coef["b"][kk])
        w = acb(0, 2 * arb.pi() * kk)
        tot += th * w ** j * (w * t).exp()
    val = 2 * tot.real
    return val + coef["mean"] if j == 0 else val


def _sgn(x):
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def extremes(coef, e1, e2):
    def Fp(t):
        return S_der(t, 1, coef) + arb(0, e1)

    def Fpp(t):
        return S_der(t, 2, coef) + arb(0, e2)

    h = arb(1) / NCELL
    zeros, all_ok = [], True
    for i in range(NCELL):
        I = arb(i) * h + h / 2 + arb(0, h / 2)
        d1 = Fp(I)
        if d1 > 0 or d1 < 0:
            continue
        d2 = Fpp(I)
        if not (d2 > 0 or d2 < 0):
            all_ok = False
            continue
        sl, sr = _sgn(Fp(arb(i) * h)), _sgn(Fp(arb(i + 1) * h))
        if sl == 0 or sr == 0:
            all_ok = False
        elif sl != sr:
            zeros.append((i, _sgn(d2)))
    locs = []
    for i, s2 in zeros:
        lo, hi = arb(i) * h, arb(i + 1) * h
        slo = _sgn(Fp(lo))
        for _ in range(400):
            mid = (lo + hi) / 2
            sm = _sgn(Fp(mid))
            if sm == 0:
                break
            if sm == slo:
                lo = mid
            else:
                hi = mid
        locs.append((lo.union(hi), s2, i))
    return all_ok, zeros, locs


def certify():
    ctx.prec = 400
    G = arb(G_STAR_HI)
    t0 = time.time()
    coef = au_modes()
    print(f"   AU's routine re-run ({K_PER} cells, b = u = {ORD}/{ORD + 1}, {DPS} digits, depth 5^-140): {time.time() - t0:.0f}s")
    return coef, G


def main() -> int:
    t0 = time.time()
    print("Theorem BV: F = theta is holomorphic on |Im t| < pi/ln5, so |theta_k| <= B_k ~ e^(-2 pi^2 k/ln5); theta_5..11 != 0; "
          "F has one max and one min per period, located to ~1e-55. Not thrust.")
    symbolic_checks()
    exact_checks()
    coef, G = certify()
    check(coef["nodes_ok"], "AQ (4) node enclosures of r_inf positive with relative width < 1e-100")
    m = coef["mean"]
    check(bool(m.rad() < arb("1e-60")), f"<F> = <Theta> enclosed to radius {float(m.rad()):.1e}: {m.str(40, radius=False)}")
    # (4) modes, nonvanishing, AU's table inside, decay bound check
    rad = max(float(coef["a"][k].rad()) for k in coef["a"])
    print(f"   common radius of a_k, b_k: {rad:.2e}")
    nz, inside, below, mods = True, True, True, {}
    x = (-2 * arb.pi() ** 2 / arb(5).log()).exp()
    for k in range(1, KMAX + 1):
        a, b = coef["a"][k], coef["b"][k]
        mlo = max(arb(0), abs(a).lower() if not a.contains(0) else arb(0))
        nlo = (abs(a).lower() ** 2 + abs(b).lower() ** 2).sqrt() if (not a.contains(0) or not b.contains(0)) else arb(0)
        nhi = (abs(a).upper() ** 2 + abs(b).upper() ** 2).sqrt()
        mods[k] = (nlo, nhi)
        Bk = B_bound(k, G)
        nz = nz and bool(nlo > 0)
        below = below and bool(nhi < Bk)
        if k <= 4:
            for c_, ref in ((a, AU_TABLE_A[k]), (b, AU_TABLE_B[k])):      # AU prints 22 significant digits
                tol = abs(arb(ref)) * arb("1e-21") + arb("2e-35")
                inside = inside and bool(abs(c_ - arb(ref)) < tol)
        print(f"   k={k:2d}: a_k = {a.mid().str(12, radius=False):>20s}  b_k = {b.mid().str(12, radius=False):>20s}  "
              f"|theta_k| >= {float(nlo):.6e}   B_k = {float(Bk.upper()):.4e}   |theta_k|/B_k <= {float((nhi / Bk).upper()):.4f}"
              f"   |theta_k| e^(2pi^2 k/ln5) ~ {float((nhi / x ** k).mid()):.4f}")
    check(nz, f"BV.4: theta_k != 0 certified for every k = 1..{KMAX} (AU certified k <= 4): M has genuine poles at 1-gamma +- i omega k")
    check(inside, "BV.4 consistency: AU's published 22-digit a_k, b_k (k <= 4) agree with the new enclosures to their printed precision")
    check(below, f"BV.3 vs BV.4: |theta_k| <= B_k for every k = 1..{KMAX} (the decay bound is consistent with an independent certificate)")
    ratios = [float((mods[k + 1][1] / mods[k][0]).upper()) for k in range(1, KMAX)]
    print("   floats: |theta_(k+1)/theta_k| = " + ", ".join(f"{r:.3e}" for r in ratios) +
          f"; e^(-2 pi^2/ln5) = {float(x.mid()):.4e}")
    # (5) extremes
    e0, e1, e2 = (tail_bound(KMAX + 1, j, G) for j in (0, 1, 2))
    print(f"   tails (k >= {KMAX + 1}): |F - S| <= {float(e0.upper()):.2e}, |F' - S'| <= {float(e1.upper()):.2e}, "
          f"|F'' - S''| <= {float(e2.upper()):.2e}")
    all_ok, zeros, locs = extremes(coef, e1.upper(), e2.upper())
    check(all_ok, f"BV.5: on all {NCELL} cells of one period, F' != 0 or F'' != 0, and every sign of F' at cell ends is certified")
    check(len(zeros) == 2 and sorted(s for _, s in zeros) == [-1, 1],
          "BV.5: F' changes sign exactly twice per period (one with F'' < 0, one with F'' > 0): exactly one maximum and one "
          "minimum of F per period, both nondegenerate")
    if len(locs) != 2:
        print("   extremes not located; stopping")
        return 1
    (tmax, _, imax), = [l for l in locs if l[1] < 0]
    (tmin, _, imin), = [l for l in locs if l[1] > 0]
    Fmax = S_der(tmax, 0, coef) + arb(0, e0.upper())
    Fmin = S_der(tmin, 0, coef) + arb(0, e0.upper())
    print(f"   t_max in {tmax.str(62)}\n   t_min in {tmin.str(62)}")
    print(f"   max F = G* in {Fmax.str(60)}\n   min F in {Fmin.str(60)}")
    check(bool(tmax.rad() < arb("1e-54")) and bool(tmin.rad() < arb("1e-54")) and bool(Fmax.rad() < arb("1e-59"))
          and bool(Fmin.rad() < arb("1e-59")), "BV.5: locations enclosed to radius < 1e-54, extreme values to < 1e-59")
    check(bool(Fmax > arb(BU_MAX[0])) and bool(Fmax < arb(BU_MAX[1])) and bool(Fmin > arb(BU_MIN[0])) and bool(Fmin < arb(BU_MIN[1])),
          "consistency: max F and min F lie inside BU.3's independent convexity brackets (F's own nodes at depth 5^-30)")
    hN = arb(1) / 2**18
    nmax, nmin = arb(BU_NODES[0]) * hN, arb(BU_NODES[1]) * hN
    check(bool(abs(tmax - nmax) < 2 * hN) and bool(abs(tmin - nmin) < 2 * hN),
          f"consistency: t_max, t_min are within two cells of BU.3's extremal nodes {BU_NODES[0]}/2^18, {BU_NODES[1]}/2^18 "
          f"(offsets {float(((tmax - nmax) / hN).mid()):+.3f}, {float(((tmin - nmin) / hN).mid()):+.3f} cells)")
    gap = tmax - tmin - arb(1) / 2
    print(f"   t_max - t_min - 1/2 in {gap.str(20)}")
    check(bool(gap > 0) or bool(gap < 0),
          f"Kept Failure BV.5: t_max - t_min - 1/2 = {float(gap.mid()):.6e} != 0 certified, so the extremes are not exactly half "
          f"a period apart (BU's nodes differ by exactly 2^17 only at resolution 2^-18)")
    half = (Fmax + Fmin) / 2 - m
    print(f"   (max F + min F)/2 - <F> in {half.str(30)}")
    check(bool(half > arb("2.434e-11")) and bool(half < arb("2.9084e-11")),
          "consistency: (max F + min F)/2 - <F> lies inside Kept Failure BU.4's bracket [2.4340, 2.9084]e-11")
    amp = Fmax - Fmin
    print(f"   max F - min F in {amp.str(30)}")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BV checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
