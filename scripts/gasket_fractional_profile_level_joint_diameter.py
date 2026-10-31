#!/usr/bin/env python3
r"""Theorem DJ: DI's kept limit "DI.1 not applied to the diameter of Lt^9 C_100" and DH's kept limit "no floor for the
cone term".  Same cone C_100, same block Lt^9 (512 words).  No new upper bound: DH.4's Delta_H < 0.025879 stays.

Setting (GASKET.md, Theorems CE, CF, DC-DI): s in X_s = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
C_a = {f > 0 : log f a-Lipschitz in s}, Theta_a its Hilbert metric (CF's definition), a = 100.

Theorem DJ (checked here; exit 1 on any failure):
  DJ.1 [proved, exact] (separation of the Hilbert diameter) With the evaluation functionals l_y(u) = u(y) and
       l_(x,y)(u) = e^(a|x-y|) u(y) - u(x) (x != y), indexed by P in E = X_s u {(x, y) : x != y},
           diam_(Theta_a) Lt^n C_a = sup_(P1, P2 in E) [rho(P1, P2) + rho(P2, P1)],
           rho(P, P') = sup_(f in C_a) log(l_P(Lt^n f)/l_P'(Lt^n f)),
       and each rho(P, P') is a ratio of two signed linear functionals over at most 4 * 2^n points, maximised by a
       bang-bang log-Lipschitz f (DE.1's line programme with signed coefficients).
  DJ.2 [certified; arb] explicit f, g in C_100 (bang-bang, slopes +-a (1 - 1e-6), point order certified) with
       Theta_100(Lt^9 f, Lt^9 g) >= Delta_low, from DI.1's quotients at two near-diagonal pairs.  So
       Delta_low <= diam Lt^9 C_100 <= Delta_H (DH.4) and every CF.1 contraction factor for (a, n) = (100, 9) is
       >= tanh(Delta_low/4).
  Observation (floats, no claim): a 57-functional grid search of DJ.1's formula.
Kept limits: see GASKET.md Theorem DJ.  No Psi_info, Stage 2 still blocked.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import itertools
import math
import os
import sys
import time

import numpy as np
from flint import arb, ctx

FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
A = 100.0
A_NUM = 100
N = 9
SLOPE_FAC = 1e-6                  # f = exp(phi), |phi'| = a (1 - 1e-6): f in C_100
D_NEAR = float(os.environ.get("DJ_D", "1e-6"))
DELTA_H = "0.025879"              # DH.4, used as stated
THETA_H = "0.0064696"             # DH.4
RATE_H = "0.57117"                # DH.4
W_H, W_LOW_H = "0.0142", "0.013764"   # DH.3 upper bound and floor
DLOW_MIN = "0.013997"
TLOW_MIN = "0.0034992"
RLOW_MIN = "0.53342"
OBS = os.environ.get("DJ_OBS", "0") == "1"
WORDS = list(itertools.product((0, 1), repeat=N))
S0, S1 = math.sqrt(5) / 2, 2.5


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- float search (DJ.1's programme)
_IMG: dict = {}


def _images(z: float):
    if z not in _IMG:
        B, Wt = [], []
        for w in WORDS:
            b, wt = z, 1.0
            for i in w:
                if i == 0:
                    wt *= (b + .5) / (2 * b * (3.5 - b)); b = math.sqrt(3.75 - b)
                else:
                    wt *= (b - .5) / (2 * b * (b + 3.5)); b = math.sqrt(3.75 + b)
            B.append(b); Wt.append(wt)
        _IMG[z] = (np.array(B), np.array(Wt))
    return _IMG[z]


def _functional(P):
    """P = ('R', y) or ('Q', x, y); returns the image points and the coefficients of l_P(Lt^9 f) in f, normalised."""
    if P[0] == "R":
        B, Wt = _images(P[1])
        return B, Wt / Wt.sum()
    _, x, y = P
    By, Wy = _images(y); Bx, Wx = _images(x)
    c = np.concatenate([math.exp(A * abs(x - y)) * Wy, -Wx])
    return np.concatenate([By, Bx]), c / c.sum()


def rho(P1, P2, ret=False):
    """sup_f log(l_P1(Lt^9 f)/l_P2(Lt^9 f)) over f = exp(phi), phi a(1-1e-6)-Lipschitz (float bisection + DE.1's DP with
    signed coefficients: sup_f (l1 - t l2)(f) > 0 iff the backward recursion V = c_k + e^(+-a d_k) V ends positive)."""
    p1, c1 = _functional(P1); p2, c2 = _functional(P2)
    s = np.concatenate([p1, p2]); Av = np.concatenate([c1, 0 * c2]); Bv = np.concatenate([0 * c1, c2])
    o = np.argsort(s, kind="stable"); s, Av, Bv = s[o], Av[o], Bv[o]
    d = np.diff(s); ap = A * (1 - SLOPE_FAC)
    eu, ed = np.exp(ap * d), np.exp(-ap * d)

    def V(t, signs=None):
        c = Av - t * Bv; v = c[-1]
        for k in range(len(c) - 2, -1, -1):
            sg = 1 if v >= 0 else -1
            v = c[k] + (eu[k] if sg > 0 else ed[k]) * v
            if signs is not None:
                signs[k] = sg
        return v
    lo, hi = math.exp(-1), math.exp(1)
    for _ in range(60):
        t = math.sqrt(lo * hi)
        lo, hi = (t, hi) if V(t) > 0 else (lo, t)
    if not ret:
        return math.log(lo)
    signs = [0] * (len(s) - 1)
    V(lo, signs)
    return math.log(lo), s, signs


# ----------------------------------------------------------------- arb evaluation of Lt^9 f for an explicit f
def _knots_arb(s, signs):
    a = arb(A_NUM)
    slope = a * (1 - arb("1e-6"))
    ks = [arb(float(x)) for x in s]       # exact dyadic knots
    for k0, k1 in zip(ks, ks[1:]):
        if not bool(k0 <= k1):
            raise RuntimeError("knots not ordered")
    phi = [arb(0)]
    for k in range(1, len(ks)):
        phi.append(phi[-1] + signs[k - 1] * slope * (ks[k] - ks[k - 1]))
    return ks, phi, slope


def _phi_at(b: arb, s, ks, phi, slope, signs):
    """phi(b) for phi piecewise linear through (ks, phi), constant outside: union of the two adjacent pieces (valid for
    every b in the ball, since phi is continuous and equals one of the pieces on each side of the nearest knot)."""
    bm = float(b.mid())
    j = int(np.searchsorted(s, bm))
    j = min(max(j, 0), len(s) - 1)
    if j > 0 and abs(s[j - 1] - bm) < abs(s[j] - bm):
        j -= 1                                 # nearest knot
    out = None
    for k in (j - 1, j):                       # pieces [ks[k], ks[k+1]]
        if k < 0:
            val = phi[0]
        elif k >= len(ks) - 1:
            val = phi[-1]
        else:
            val = phi[k] + signs[k] * slope * (b - ks[k])
        out = val if out is None else out.union(val)
    # the ball must lie within [ks[j-1], ks[j+1]] (or beyond the end knot) so that only these pieces matter
    if not ((j - 1 < 0 or bool(b >= ks[j - 1])) and (j + 1 > len(ks) - 1 or bool(b <= ks[j + 1]))):
        raise RuntimeError("point not localised")
    return out


def _Lt_arb(z: arb, f):
    s, ks, phi, slope, signs = f
    tot = arb(0)
    for w in WORDS:
        b, wt = z, arb(1)
        for i in w:
            if i == 0:
                wt *= (b + arb(1) / 2) / (2 * b * (arb(7) / 2 - b)); b = (arb(15) / 4 - b).sqrt()
            else:
                wt *= (b - arb(1) / 2) / (2 * b * (b + arb(7) / 2)); b = (arb(15) / 4 + b).sqrt()
        tot += wt * _phi_at(b, s, ks, phi, slope, signs).exp()
    return tot


def _l_arb(P, f, pts):
    if P[0] == "R":
        return _Lt_arb(pts[P[1]], f)
    _, x, y = P
    d = abs(pts[x] - pts[y])
    return (arb(A_NUM) * d).exp() * _Lt_arb(pts[y], f) - _Lt_arb(pts[x], f)


def floor_certificate():
    ctx.prec = PREC
    t0 = time.time()
    s0a, s1a = arb(5).sqrt() / 2, arb(5) / 2
    dn = arb(D_NEAR)
    # P1: x < y at the right end, P2: x < y at the left end (best pair of the float search)
    pts_f = {"s0": S0, "s0d": S0 + D_NEAR, "s1": S1, "s1d": S1 - D_NEAR}
    pts = {"s0": s0a, "s0d": s0a + dn, "s1": s1a, "s1d": s1a - dn}
    P1, P2 = ("Q", "s1d", "s1"), ("Q", "s0", "s0d")
    fl = lambda P: (P[0], pts_f[P[1]], pts_f[P[2]])   # noqa: E731
    r12, sg_s, sg = rho(fl(P1), fl(P2), ret=True)     # g: maximises l_P1/l_P2
    r21, sf_s, sf = rho(fl(P2), fl(P1), ret=True)     # f: maximises l_P2/l_P1
    print(f"   DJ.2: float programme rho(P1,P2) = {r12:.9f}, rho(P2,P1) = {r21:.9f}, sum {r12 + r21:.9f} "
          f"({time.time() - t0:.1f}s)", flush=True)
    g = (sg_s,) + _knots_arb(sg_s, sg)[:3] + (sg,)
    f = (sf_s,) + _knots_arb(sf_s, sf)[:3] + (sf,)
    ctx.prec = PREC
    check(bool(pts["s1d"] < pts["s1"]) and bool(pts["s0"] < pts["s0d"]) and bool(pts["s0d"] < pts["s1d"]),
          f"DJ.2: the four points sqrt5/2, sqrt5/2 + d, 5/2 - d, 5/2 (d = {D_NEAR:g}) are ordered in X_s")
    vals = {}
    for nm, h in (("u", f), ("v", g)):
        for P in (P1, P2):
            vals[nm, P] = _l_arb(P, h, pts)
    pos = all(bool(vals[k] > 0) for k in vals)
    check(pos, "DJ.2: l_P(u), l_P(v) > 0 at both pairs (u = Lt^9 f, v = Lt^9 g; arb)")
    Q1 = vals["v", P1] / vals["u", P1]
    Q2 = vals["v", P2] / vals["u", P2]
    low = Q1.log() - Q2.log()
    # DI.1 consistency at the same pairs: Phi_v/Phi_u = R(y) + (R(y) - R(x))/(e^(ad) u(y)/u(x) - 1)
    okid = True
    for P, Qp in ((P1, Q1), (P2, Q2)):
        x, y = pts[P[1]], pts[P[2]]
        ux, uy, vx, vy = _Lt_arb(x, f), _Lt_arb(y, f), _Lt_arb(x, g), _Lt_arb(y, g)
        Rx, Ry = vx / ux, vy / uy
        rhs = Ry + (Ry - Rx) / ((arb(A_NUM) * abs(x - y)).exp() * uy / ux - 1)
        okid = okid and rhs.overlaps(Qp) and bool(rhs.rad() < arb("1e-6") * abs(Qp))
        RESULT.setdefault("R", []).append((Rx, Ry))
    check(okid, "DJ.2 (consistency): DI.1's identity reproduces both quotients (overlapping arb balls, relative radius < 1e-6)")
    dH = arb(DELTA_H)
    check(bool(low > arb(DLOW_MIN)) and bool(low < dH),
          f"DJ.2: Theta_100(Lt^9 f, Lt^9 g) >= log(Q_P1/Q_P2) = {low.str(10)} > {DLOW_MIN}, so "
          f"{DLOW_MIN} < diam Lt^9 C_100 <= Delta_H < {DELTA_H} (DH.4)")
    # the pair's diagonal part alone (DH.3's quantity) and the cone part realised by the pair
    Rs = [r for pair in RESULT["R"] for r in pair]
    Rlog = [r.log() for r in Rs]
    wpart = arb(max(r.upper() for r in Rlog)) - arb(min(r.lower() for r in Rlog))
    gap = arb((arb(DELTA_H) - low).upper())
    check(bool(low > arb(W_LOW_H)) and bool(gap < arb("0.011882")),
          f"DJ.2: Delta_low exceeds DH.3's diagonal floor W_low = {W_LOW_H} by {(low - arb(W_LOW_H)).str(5)}; the true "
          f"diameter lies in an interval of width Delta_H - Delta_low < {float(gap.upper()):.7f} (DH's cone term 2 log K < 0.011679)")
    theta_low = (low / 4).tanh()
    rate_low = (theta_low.log() / N).exp()
    check(bool(theta_low > arb(TLOW_MIN)) and bool(rate_low > arb(RLOW_MIN)),
          f"DJ.2: every CF.1 factor tanh(Delta/4) with Delta >= diam Lt^9 C_100 is >= tanh(Delta_low/4) = "
          f"{theta_low.str(8)} > {TLOW_MIN}; per-level rate >= {rate_low.str(8)} > {RLOW_MIN} (DH.4: theta_H < {THETA_H}, "
          f"rate < {RATE_H}; spectral 0.1725)")
    print(f"   DJ.2: certified in {time.time() - t0:.1f}s")
    RESULT.update(gap=gap, low=low, theta_low=theta_low, rate_low=rate_low, r12=r12, r21=r21, wpart=wpart)


def observation():
    """float grid search of DJ.1's formula (no claim)."""
    ys = [S0 + (S1 - S0) * k / 8 for k in range(9)]
    Ps = []
    for y in ys:
        Ps.append(("R", y))
        for d in (1e-6, 1e-3, 1e-2):
            if y + d <= S1:
                Ps.append(("Q", y + d, y))
            if y - d >= S0:
                Ps.append(("Q", y - d, y))
    best = (0, None, None)
    for i, P in enumerate(Ps):
        for j, Q in enumerate(Ps):
            if i < j:
                v = rho(P, Q) + rho(Q, P)
                if v > best[0]:
                    best = (v, P, Q)
    print(f"   Observation (floats): {len(Ps)} functionals, max rho(P1,P2) + rho(P2,P1) = {best[0]:.6f} at {best[1]}, {best[2]}")
    RESULT.update(obs=best)


def main() -> int:
    t0 = time.time()
    print("Theorem DJ: separation of the Hilbert diameter of Lt^9 C_100 and a certified floor. Not thrust.", flush=True)
    floor_certificate()
    if OBS:
        observation()
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem DJ checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
