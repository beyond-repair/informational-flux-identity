#!/usr/bin/env python3
r"""Theorem DS: DR's open case n = 4 settled.  The Hilbert diameter D_4 = diam Lt^4 K_(2,6) of the second-order cone
K_(2,6) = {f > 0: |f'| <= 2 f, |f''| <= 6 f} (DQ.1) satisfies 0.0360461 <= D_4 <= 0.0383859 < 0.0404321 = 4 artanh(q^4),
q = 0.3170757 (DP's cap), so tanh(D_4/4)^(1/4) in [0.3081036, 0.3129858] is below the cap and DQ.4's threshold is n0 = 4
exactly.  (DR's open item: "n = 4: certified upper bound 0.3433 per level, floats on both sides of the cap".)

Setting (GASKET.md, Theorems CF, CI-CM, DP, DQ, DR): s in X = [sqrt5/2, 5/2], Lt^n f(s) = sum_w W_w(s) f(y_w(s)) over the
2^n words; K = K_(2,6) with Hilbert metric read through Lam = {2 d_s +- d'_s, 6 d_s +- d''_s : s in X}, Phi_l(f) = l(Lt^n f).
For f in K write g = log f, v = g' = f'/f: then |v| <= 2 and v' = f''/f - v^2 in [-6 - v^2, 6 - v^2] (a.e.).

Theorem DS (checked here; exit 1 on any failure):
  DS.1 [proved; sympy identities] Comparison lemma.  U_a(t) = sqrt6 tanh(sqrt6 t + artanh(a/sqrt6)) solves U' = 6 - U^2,
       U(0) = a; W_b(s) = sqrt6 tan(sqrt6 s + arctan(b/sqrt6)) solves W' = 6 + W^2, W(0) = b; int U = log cosh(..)/cosh(c),
       int W = -log cos(..)/cos(c).  For f in K, y < z in X, d = z - y, a = v(y), b = v(z): v(y + t) <= min(2, U_a(t)) and
       v(z - s) <= min(2, W_b(s)) (ODE comparison), so for every t0 in [0, d]
           log f(z) - log f(y) <= B(a, b, d) := int_0^t0 min(2, U_a) + int_0^(d-t0) min(2, W_b),
       B is nondecreasing in a and b; the reversal v(t) -> -v(d - t) maps the constraint set to itself, so
       log f(z) - log f(y) >= -B(-b, -a, d); and v(z) <= min(2, U_a(d)), v(z) >= -min(2, W_(-a)(d)).
  DS.2 [proved] Value-function certificate.  Let Psi = sum_k al_k d_(y_k) + be_k d'_(y_k) + ga_k d''_(y_k), y_0 < ... < y_N,
       and cut [-2, 2] into cells I_i.  Put L_N(i) <= inf_(I_i) (al_N + be_N v) - 6|ga_N| and
           L_k(i) <= inf_(I_i)(al_k + be_k v) - 6|ga_k| + min_(j reachable) { e^(-B(-v_j, -v_i)) L_(k+1)(j) if L_(k+1)(j) >= 0,
                                                                           e^(B(v_(i+1), v_(j+1))) L_(k+1)(j) otherwise }.
       If min_i L_0(i) >= 0 then Psi >= 0 on K (backward induction on S_k = sum_(j>=k) e^(g(y_j) - g(y_k)) (al_j + be_j v_j +
       ga_j f''/f), S_k = c_k + e^(g(y_(k+1)) - g(y_k)) S_(k+1), using DS.1 on each gap).  B is evaluated from arb tables of
       exp(int min(2, U)), exp(int min(2, W)) on a t-grid of step h = 2^-15 (one slack term 2|s - s'| <= 4h, the integrands
       are bounded by 2), and the recursion runs in IEEE doubles with every operation rounded outward (nextafter): a
       rigorous lower bound.  This replaces DR.1's chain (which loses a factor ~1.5 at n = 4).
  DS.3 [certified; DR.2's cell lemma unchanged, 200 cells x 4 types, reference l0 = d_(21/10), 256 v-cells]
       D_4 <= 2 max_cells log(t1 t2) <= 0.0383859, tanh(D_4/4)^(1/4) <= 0.3129858 < 0.3170757 <= rho(L_1)/r.  With DQ.5
       (n = 1, 2, 3 above the cap) and DR.4 (every n >= 5 below): n0 = 4 exactly, and tanh(D_n/4)^(1/n) < 0.3170757 for
       every n >= 4.  (DS_FULL=1 runs all 200 cells; the default runs 8 cells, including the two worst ones.)
  DS.4 [certified lower bound, arb] Explicit f, g in K (piecewise f'' = +-6 f or constant log-slope, C^1, |f'/f| <= 2 checked
       at every breakpoint, monotone Riccati v' = u - v^2 in between, positivity of the cos pieces checked): with
       l = 2 d_s - d'_s at s = sqrt5/2 and l' = 2 d_s + d'_s at s = 5/2,
           D_4 >= Theta(Lt^4 f, Lt^4 g) >= log[(l(Lt^4 f) l'(Lt^4 g))/(l'(Lt^4 f) l(Lt^4 g))] >= 0.0360461,
       so tanh(D_4/4)^(1/4) >= 0.3081036: the margin below the cap is real but < 2.9% per level.  (f, g are float
       proposals from DS.2's recursion with argmax backtracking; only their arb evaluation is claimed.  They lie in
       W^(2,oo); mollifying their constant-slope extensions gives C^oo functions of K with the same ratio in the limit.)
Kept limits: D_4 is enclosed to within 6.5% (0.0360461 .. 0.0383859), not computed; the float recursion (M = 256) at the worst
pair gives 0.03749, refining with M (0.03849 at M = 128), consistent with the lower witness.  The l0 reference is kept;
pairwise cells were not needed.  No Psi_info, Stage 2 still blocked.  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW
enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np
import sympy as sp
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_essential as ci                 # noqa: E402  (rho(L_1), r as in CI.1/CG.3)
import gasket_fractional_profile_level_effective_threshold as dr       # noqa: E402  (DR.2's word jets and cell lemma)

FAILS: list[str] = []
RESULT: dict = {}
N = 4
N_CELLS = int(os.environ.get("DS_CELLS", "200"))
MV = int(os.environ.get("DS_M", "256"))             # v-cells on [-2, 2] (power of 2: nodes exact in binary)
HEXP = int(os.environ.get("DS_H", "15"))            # t-grid step h = 2^-HEXP
S0 = "21/10"                                       # reference functional l0 = d_(21/10)
CAP = "0.3170757"                                   # DP.3 / DQ.3: rho(L_1)/r >= 0.317075713
FULL = os.environ.get("DS_FULL") == "1"
DEFAULT_ONLY = "0,1,2,99,100,197,198,199"
ONLY = None if FULL else [int(x) for x in os.environ.get("DS_ONLY", DEFAULT_ONLY).split(",")]
DN = lambda x: np.nextafter(x, -np.inf)
UP = lambda x: np.nextafter(x, np.inf)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- DS.1 identities (sympy)
def symbolic_checks():
    t, a, b, c = sp.symbols("t a b c", real=True)
    r6 = sp.sqrt(6)
    U = r6 * sp.tanh(r6 * t + sp.atanh(a / r6))
    W = r6 * sp.tan(r6 * t + sp.atan(b / r6))
    ok1 = (sp.simplify(sp.diff(U, t) - (6 - U ** 2)) == 0 and sp.simplify(U.subs(t, 0) - a) == 0
           and sp.simplify(sp.diff(W, t) - (6 + W ** 2)) == 0 and sp.simplify(W.subs(t, 0) - b) == 0)
    ok2 = (sp.simplify(sp.diff(sp.log(sp.cosh(r6 * t + c)), t) - r6 * sp.tanh(r6 * t + c)) == 0
           and sp.simplify(sp.diff(-sp.log(sp.cos(r6 * t + c)), t) - r6 * sp.tan(r6 * t + c)) == 0)
    check(ok1 and ok2, "DS.1 (sympy): U_a = sqrt6 tanh(sqrt6 t + artanh(a/sqrt6)) solves U' = 6 - U^2, U(0) = a; W_b = sqrt6 "
                       "tan(sqrt6 t + arctan(b/sqrt6)) solves W' = 6 + W^2, W(0) = b; (log cosh(sqrt6 t + c))' = U, "
                       "(-log cos(sqrt6 t + c))' = W")
    # the reversal w(s) = -v(d - s) and the log-derivative identity behind v' = f''/f - v^2
    s, d = sp.symbols("s d", real=True)
    v = sp.Function("v")
    f = sp.Function("f")
    x = sp.symbols("x", real=True)
    ok3 = sp.simplify(sp.diff(sp.diff(f(x), x) / f(x), x) - (sp.diff(f(x), x, 2) / f(x) - (sp.diff(f(x), x) / f(x)) ** 2)) == 0
    w = -v(d - s)
    ok4 = sp.simplify(sp.diff(w, s) - sp.Subs(sp.diff(v(x), x), x, d - s).doit()) == 0
    # pieces of the DS.4 witnesses: f'' = u f closed forms
    f0, fp0, om = sp.symbols("f0 fp0 omega", positive=True)
    ch = f0 * sp.cosh(om * t) + fp0 / om * sp.sinh(om * t)
    co = f0 * sp.cos(om * t) + fp0 / om * sp.sin(om * t)
    ok5 = (sp.simplify(sp.diff(ch, t, 2) - om ** 2 * ch) == 0 and sp.simplify(sp.diff(co, t, 2) + om ** 2 * co) == 0
           and ch.subs(t, 0) == f0 and sp.simplify(sp.diff(ch, t).subs(t, 0) - fp0) == 0
           and co.subs(t, 0) == f0 and sp.simplify(sp.diff(co, t).subs(t, 0) - fp0) == 0)
    check(ok3 and ok4 and ok5, "DS.1/DS.4 (sympy): (f'/f)' = f''/f - (f'/f)^2; w(s) = -v(d-s) has w'(s) = v'(d-s) (so the "
                               "set v' in [-6 - v^2, 6 - v^2] is reversal invariant: w' in [-6 - w^2, 6 - w^2]); "
                               "f0 cosh(wt) + (f0'/w) sinh(wt), f0 cos(wt) + (f0'/w) sin(wt) solve f'' = +-w^2 f")


# ----------------------------------------------------------------- DS.2: arb tables
class Tables:
    """EU[p,m] >= exp(int_0^(mh) min(2, U_(v_p))), EW[q,m] >= exp(int_0^(mh) min(2, W_(v_q))),
    RU[p,m] >= min(2, U_(v_p)(mh)), RW[p,m] >= min(2, W_(v_p)(mh)), v_p = -2 + 4p/M."""

    def __init__(self, M, hexp, Tmax):
        ctx.prec = 80
        self.M, self.h = M, 2.0 ** -hexp
        self.v = np.array([-2 + 4 * p / M for p in range(M + 1)])
        self.NT = int(math.ceil(Tmax / self.h)) + 2
        r6, two = arb(6).sqrt(), arb(2)
        EU, EW, RU, RW = (np.empty((M + 1, self.NT)) for _ in range(4))
        for p in range(M + 1):
            a = arb(self.v[p])
            ca = (a / r6).atanh()
            tU = ((two / r6).atanh() - ca) / r6               # U_a reaches 2
            cb = (a / r6).atan()
            sW = ((two / r6).atan() - cb) / r6                # W_a reaches 2 (before its pole)
            lchc, lcc = ca.cosh().log(), cb.cos().log()
            eU2 = (r6 * tU + ca).cosh().log() - lchc
            eW2 = lcc - (r6 * sW + cb).cos().log()
            for m in range(self.NT):
                t = arb(m) * arb(self.h)
                # int min(2,U) <= int U (any t); = int_0^tU U + 2(t - tU) once t > tU.  Same for W.
                lu = eU2 + two * (t - tU) if t > tU else (r6 * t + ca).cosh().log() - lchc
                lw = eW2 + two * (t - sW) if t > sW else lcc - (r6 * t + cb).cos().log()
                EU[p, m] = UP(float(lu.exp().upper()))
                EW[p, m] = UP(float(lw.exp().upper()))
                RU[p, m] = min(2.0, UP(float((r6 * (r6 * t + ca).tanh()).upper())))
                RW[p, m] = min(2.0, UP(float((r6 * (r6 * t + cb).tan()).upper()))) if t < sW else 2.0
        self.EU, self.EW, self.RU, self.RW = EU, EW, RU, RW
        self.C = UP(float((arb(4) * arb(self.h)).exp().upper()))     # e^(2 slack), slack <= 2h


R6 = math.sqrt(6.0)


def _U(a, t):
    return np.minimum(R6 * np.tanh(R6 * t + np.arctanh(a / R6)), 2.0)


def _W(b, s):
    th = R6 * s + np.arctan(b / R6)
    return np.where(th < np.pi / 2, np.minimum(R6 * np.tan(np.minimum(th, 1.5)), 2.0), 2.0)


def _cross(a, b, d):
    """float: split point t0 (any t0 in [0, d] gives a valid bound; this one makes it tight)."""
    a, b = np.broadcast_arrays(a, b)
    lo, hi = np.zeros(a.shape), np.full(a.shape, d)
    for _ in range(40):
        m = (lo + hi) / 2
        go = _U(a, m) - _W(b, d - m) < 0
        lo, hi = np.where(go, m, lo), np.where(go, hi, m)
    return (lo + hi) / 2


class Gap:
    """one gap d in [dlo, dhi]: Eh[i,j] >= e^B(v_(i+1), v_(j+1), d), Ec[i,j] <= e^(-B(-v_j, -v_i, d)); ok[i,j] False only
    if cell j is certainly unreachable from cell i (DS.1)."""

    def __init__(self, T, dlo, dhi):
        assert 0 < dlo <= dhi and dhi - dlo < T.h
        M, h, v = T.M, T.h, T.v
        t0 = _cross(v[:, None], v[None, :], 0.5 * (dlo + dhi))
        m = np.clip(np.floor(t0 / h).astype(int), 0, int(math.floor(dlo / h)))       # t0 := m h <= dlo
        mp = np.ceil((dhi - m * h) / h).astype(int)
        mp = np.where(mp * h < dhi - m * h, mp + 1, mp)                                # m'h >= d - mh
        assert np.all(UP((m + mp) * h - dlo) <= 2 * h)                                 # slack <= 2h
        idx = np.arange(M + 1)
        B = UP(UP(T.EU[idx[:, None], m] * T.EW[idx[None, :], mp]) * T.C)              # >= e^B(v_p, v_q, d)
        self.Eh = B[1:, 1:]
        self.Ec = DN(1.0 / B[::-1, ::-1].T[:-1, :-1])                                  # -v_i = v_(M-i)
        mr = int(math.ceil(dhi / h))
        up, dn = T.RU[:, mr], T.RW[:, mr]
        self.ok = ~((v[None, :-1] > up[1:, None]) | (v[None, 1:] < -dn[::-1][:-1, None]))


def rig_minL(T, gaps, al, be):
    """DS.2 in outward-rounded doubles; al, be: (N, 2) arrays of float intervals.  Lower bound for min_i L_0(i)."""
    v0, v1 = T.v[:-1], T.v[1:]

    def c(k):
        bl, bu = be[k]
        pr = np.minimum.reduce([DN(bl * v0), DN(bl * v1), DN(bu * v0), DN(bu * v1)])
        return DN(al[k][0] + pr)

    L = c(len(al) - 1)
    for k in range(len(al) - 2, -1, -1):
        G, Ln = gaps[k], L[None, :]
        ph = np.where(Ln >= 0, DN(G.Ec * Ln), DN(G.Eh * Ln))
        L = DN(c(k) + np.where(G.ok, ph, np.inf).min(axis=1))
    return float(L.min())


def float_minL(T, gaps, al, be):
    v0, v1 = T.v[:-1], T.v[1:]
    c = lambda k: al[k] + np.minimum(be[k] * v0, be[k] * v1)
    L = c(len(al) - 1)
    for k in range(len(al) - 2, -1, -1):
        G, Ln = gaps[k], L[None, :]
        L = c(k) + np.where(G.ok, np.where(Ln >= 0, G.Ec * Ln, G.Eh * Ln), np.inf).min(axis=1)
    return L.min()


def _iv(x):
    return (DN(float(x.lower())), UP(float(x.upper())))


# ----------------------------------------------------------------- DS.3: cells (DR.2) + DS.2 certificates
def certify4(T, ncell=N_CELLS, only=None):
    ctx.prec = 96
    t0 = time.time()
    AX, BX = arb(5).sqrt() / 2, arb(5) / 2
    Z = dr.word_jets(arb(21) / 10, N)
    y0 = [w[1][0] for w in Z]
    w0 = [w[0][0] for w in Z]
    sQ = arb(1 / float(sum(w0).mid()))
    w0 = [w * sQ for w in w0]
    edges = np.linspace(float(AX.mid()), float(BX.mid()), ncell + 1)
    fl = lambda L: np.array([float(v.mid()) for v in L])
    worst, wcell, bumps, ncert, cover_ok, prev_hi = arb(0), None, 0, 0, True, None
    for ci_ in range(ncell):
        if only is not None and ci_ not in only:
            continue
        sc = (edges[ci_] + edges[ci_ + 1]) / 2
        dl = (edges[ci_ + 1] - edges[ci_]) / 2 * (1 + 1e-9)
        cball = arb(sc, dl)
        if only is None:
            if ci_ == 0:
                cover_ok &= bool(cball.lower() <= AX.lower())
            if ci_ == ncell - 1:
                cover_ok &= bool(cball.upper() >= BX.upper())
            if prev_hi is not None:
                cover_ok &= bool(cball.lower() <= prev_hi)
            prev_hi = cball.upper()
        cen, ball = dr.word_jets(arb(sc), N), dr.word_jets(cball, N)
        ys = [w[1][0] for w in cen]
        yall = ys + y0
        o = np.argsort(fl(yall), kind="stable")
        ya = [yall[k] for k in o]
        gaps = []
        for k in range(len(ya) - 1):
            lo_, hi_ = _iv(ya[k + 1] - ya[k])
            if not lo_ > 0:
                check(False, f"DS.3: points not separated at cell {ci_}")
                return None
            gaps.append(Gap(T, lo_, hi_))
        nP, z0 = len(ys), [arb(0)] * len(y0)
        for typ in ("1+", "1-", "2+", "2-"):
            _, a0, a1, b0, b1, E = dr.cell_functional(typ, cen, ball, dl)
            S = arb(1 / float(sum(a0).mid()))
            a0, a1, b0, b1, E = ([v * S for v in L] for L in (a0, a1, b0, b1, E))
            tv = []
            for direction in (1, 2):
                tmax = 0.0
                for tau in (-dl, dl):
                    Tt = arb(tau)
                    G = [p + Tt * q for p, q in zip(a0, a1)]
                    Gb = [p + Tt * q for p, q in zip(b0, b1)]
                    if direction == 1:       # t Phi_0 - (G + E) >= 0  =>  R(l_s, l0) <= t
                        aP, bP = [-(g + e) for g, e in zip(G, E)] + z0, [-g for g in Gb] + z0
                        aQ, bQ = [arb(0)] * nP + w0, [arb(0)] * (nP + len(y0))
                    else:                    # t (G - E) - Phi_0 >= 0  =>  R(l0, l_s) <= t
                        aP, bP = [arb(0)] * nP + [-w for w in w0], [arb(0)] * (nP + len(y0))
                        aQ, bQ = [g - e for g, e in zip(G, E)] + z0, Gb + z0
                    aP, bP, aQ, bQ = ([L[k] for k in o] for L in (aP, bP, aQ, bQ))
                    faP, fbP, faQ, fbQ = fl(aP), fl(bP), fl(aQ), fl(bQ)
                    feas = lambda t: float_minL(T, gaps, t * faQ + faP, t * fbQ + fbP) >= 0
                    lo, hi = 0.9, 1.0
                    while not feas(hi):
                        hi *= 1.05
                    while hi - lo > 1e-10:
                        mid = (lo + hi) / 2
                        lo, hi = (lo, mid) if feas(mid) else (mid, hi)
                    t = hi * (1 + 1e-9)
                    for bump in range(8):
                        At = arb(t)
                        al = np.array([_iv(At * q + p) for q, p in zip(aQ, aP)])
                        be = np.array([_iv(At * q + p) for q, p in zip(bQ, bP)])
                        if rig_minL(T, gaps, al, be) >= 0:
                            break
                        t, bumps = t * (1 + 1e-8 * 4 ** bump), bumps + 1
                    else:
                        check(False, f"DS.3: certificate failed at cell {ci_} {typ} dir {direction}")
                        return None
                    ncert += 1
                    tmax = max(tmax, t)
                tv.append(tmax)
            val = (arb(tv[0]) * arb(tv[1])).log()
            if val.upper() > worst.upper():
                worst, wcell = arb(val.upper()), (ci_, typ, round(float(sc), 5))
        if only is None and ci_ % 20 == 0:
            print(f"     cell {ci_}: 2 max log(t1 t2) so far {float((2 * worst).upper()):.7f} at {wcell}, "
                  f"{time.time() - t0:.0f}s", flush=True)
    if only is None:
        check(cover_ok, f"DR.2 cover: the {ncell} cell balls cover X = [sqrt5/2, 5/2] (consecutive overlap, ends outside)")
    Dn = arb((2 * worst).upper())
    th = (Dn / 4).tanh()
    print(f"     {ncert} certificates (DS.2, outward-rounded), {bumps} t-bumps, worst cell {wcell}, {time.time() - t0:.0f}s",
          flush=True)
    return dict(D=Dn, theta=th, per=(th.log() / N).exp(), worst=wcell, ncert=ncert)


# ----------------------------------------------------------------- DS.4: explicit lower-bound witnesses
def _step(f, fp, u, t):
    w = abs(u).sqrt()
    if u > 0:
        return f * (w * t).cosh() + fp / w * (w * t).sinh(), f * w * (w * t).sinh() + fp * (w * t).cosh()
    return f * (w * t).cos() + fp / w * (w * t).sin(), -f * w * (w * t).sin() + fp * (w * t).cos()


def build_witness(v0, segs, margin=1e-12):
    """segs [(x_k, u_k)], u in {6, -6, 4 (hold at the cap)}, last entry marks the end.  f(x_0) = 1, f'(x_0) = v0.  A hold is
    realised with constant log-slope (f'' = v^2 f, v = f'/f at its start); an ascent/descent that would leave |v| <= 2 is
    cut (float time, checked in arb) and continued with constant log-slope; outside [x_0, x_end] constant log-slope.
    Between breakpoints v solves v' = u - v^2 with u constant, so v is monotone: endpoint checks give |v| <= 2 throughout.
    Returns (pieces, failed checks)."""
    xs, us = [x for x, _ in segs], [u for _, u in segs]
    bad = []
    f, fp = arb(1), arb(v0)
    if not abs(fp) <= 2:
        bad.append("v0")
    pieces = [(-np.inf, "exp", arb(v0), arb(xs[0]), f)]

    def cospos(f, fp, u, t, k):
        w = (-u).sqrt()     # f0 cos(wt) + (f0'/w) sin(wt) > 0 on [0, t] iff w t < pi/2 + atan(v0/w)
        if not (w * t < arb.pi() / 2 + (fp / f / w).atan()):
            bad.append(f"positivity {k}")

    for k in range(len(xs) - 1):
        x, xe = arb(xs[k]), arb(xs[k + 1])
        if us[k] == 4:
            v = fp / f
            pieces.append((xs[k], "exp", v, x, f))
            e = (v * (xe - x)).exp()
            f, fp = f * e, fp * e
            continue
        u = arb(us[k])
        f2, fp2 = _step(f, fp, u, xe - x)
        if abs(fp2 / f2) <= 2:
            pieces.append((xs[k], "ode", u, x, (f, fp)))
            if us[k] < 0:
                cospos(f, fp, u, xe - x, k)
            f, fp = f2, fp2
            if not (f > 0 and abs(fp / f) <= 2):
                bad.append(f"slope {k}")
            continue
        lo, hi = 0.0, xs[k + 1] - xs[k]
        for _ in range(80):
            mid = (lo + hi) / 2
            g2, gp2 = _step(f, fp, u, arb(mid))
            lo, hi = (mid, hi) if abs(gp2 / g2) <= 2 - margin else (lo, mid)
        xc = xs[k] + lo
        pieces.append((xs[k], "ode", u, x, (f, fp)))
        if us[k] < 0:
            cospos(f, fp, u, arb(xc) - x, k)
        f, fp = _step(f, fp, u, arb(xc) - x)
        if not (f > 0 and abs(fp / f) <= 2):
            bad.append(f"cut {k}")
        v = fp / f
        pieces.append((xc, "exp", v, arb(xc), f))
        e = (v * (xe - arb(xc))).exp()
        f, fp = f * e, fp * e
    pieces.append((xs[-1], "exp", fp / f, arb(xs[-1]), f))
    return pieces, bad


def _feval(pieces, y):
    """f(y), f'(y) for a ball y: union over the pieces whose interval meets it (f is C^1, each piece formula is exact
    on its own interval)."""
    lo, hi = float(y.lower()), float(y.upper())
    starts = [p[0] for p in pieces] + [np.inf]
    res = None
    for i, (_, kind, u, x0, st) in enumerate(pieces):
        if starts[i + 1] >= lo and starts[i] <= hi:
            if kind == "exp":
                fy = st * (u * (y - x0)).exp()
                val = (fy, u * fy)
            else:
                val = _step(st[0], st[1], u, y - x0)
            res = val if res is None else (res[0].union(val[0]), res[1].union(val[1]))
    return res


def _lval(s, typ, pieces):
    tot = arb(0)
    for W, y in dr.word_jets(s, N):
        A, B, _ = dr._AB(typ, W, y)
        fy, fpy = _feval(pieces, y[0])
        tot += A[0] * fy + B[0] * fpy
    return tot


def lower_witness():
    ctx.prec = 96
    sL, sH = arb(5).sqrt() / 2, arb(5) / 2
    vals = {}
    for name, (v0, segs) in (("f", WIT_F), ("g", WIT_G)):
        pieces, bad = build_witness(v0, segs)
        check(not bad, f"DS.4 witness {name}: {len(segs)} breakpoints, f > 0, |f'/f| <= 2 at every breakpoint (monotone "
                       f"Riccati in between), cos pieces positive, |f''| <= 6 f piecewise{'' if not bad else ' ' + str(bad)}")
        vals[name] = (_lval(sL, "1-", pieces), _lval(sH, "1+", pieces))
    (Pf, Qf), (Pg, Qg) = vals["f"], vals["g"]
    pos = all(bool(x > 0) for x in (Pf, Qf, Pg, Qg))
    D = ((Pf / Qf) * (Qg / Pg)).log()
    return pos, D


def main() -> int:
    t0 = time.time()
    print("Theorem DS: D_4 = diam Lt^4 K_(2,6) certified below DP's cap (DR's open case n = 4), so n0 = 4; explicit "
          "lower witness. Not thrust.", flush=True)
    symbolic_checks()
    ctx.prec = 96
    cap = arb(CAP)
    check(bool(arb(ci.LAM_LO) / arb(ci.R_HI) > cap), f"cap: rho(L_1)/r >= {ci.LAM_LO}/{ci.R_HI} > {CAP} (CI.1, CG.3 as "
                                                      "stated in DQ.3)")
    Dcap = 4 * (cap ** 4).atanh()
    tb = time.time()
    T = Tables(MV, HEXP, 0.62)
    print(f"     DS.2 tables: {MV} v-cells, h = 2^-{HEXP}, {T.NT} t-nodes, arb, {time.time() - tb:.0f}s", flush=True)
    r = certify4(T, only=ONLY)
    if r is None:
        print(f"{len(FAILS)} FAILED")
        return 1
    RESULT.update(r)
    below = bool(r["per"] < cap)
    scope = "all 200 cells" if FULL else f"cells {DEFAULT_ONLY if ONLY is None else ','.join(map(str, ONLY))} (DS_FULL=1: all)"
    check(below, f"DS.3 ({scope}): D_4 <= {r['D'].upper().str(7, radius=False)} < 4 artanh(q^4) = "
                 f"{Dcap.lower().str(7, radius=False)}, tanh(D_4/4)^(1/4) <= {r['per'].upper().str(7, radius=False)} < {CAP}")
    if FULL and below:
        check(True, "DS.3: with DQ.5 (tanh(D_n/4)^(1/n) >= tau(Lt^n)^(1/n) > 0.3170762 for n = 1, 2, 3) and DR.4 (< 0.3135158 "
                    "for every n >= 5): tanh(D_n/4)^(1/n) < 0.3170757 for every n >= 4, so n0 = 4 exactly")
    pos, Dlow = lower_witness()
    perlow = ((Dlow / 4).tanh().log() / N).exp()
    RESULT.update(Dlow=Dlow, perlow=perlow)
    check(pos and bool(Dlow > 0), f"DS.4: D_4 >= Theta(Lt^4 f, Lt^4 g) >= {Dlow.lower().str(7, radius=False)}, "
                                  f"tanh(D_4/4)^(1/4) >= {perlow.lower().str(7, radius=False)} (l = 2d - d' at sqrt5/2, "
                                  f"l' = 2d + d' at 5/2)")
    check(bool(Dlow < r["D"]), f"DS.3/DS.4 consistent: {Dlow.lower().str(7, radius=False)} <= D_4 <= "
                               f"{r['D'].upper().str(7, radius=False)}; per level in [{perlow.lower().str(7, radius=False)}, "
                               f"{r['per'].upper().str(7, radius=False)}], cap/upper > "
                               f"{(cap / r['per']).lower().str(5, radius=False)}")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem DS checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


# DS.4 witnesses: float proposals (DS.2's recursion with argmax backtracking, bang-bang u in {6, -6}, 4 = hold at the cap);
# (v0, [(x_k, u_k), ...]); only their arb evaluation above is claimed.
WIT_F = (-0.1399999999999999, [
    (1.118033988749895, 6), (1.1182025361668166, -6), (1.1233499918526528, -6), (1.134218348509089, 6),
    (1.1343295524038288, -6), (1.1346971912358035, 6), (1.1443050506249453, -6), (1.1446466740863193, 6),
    (1.196733269327442, 6), (1.1970481714113521, -6), (1.2115084808056589, 6), (1.2117265777602466, -6),
    (1.2311071719479239, -6), (1.2421757979564427, 6), (1.2424322268284247, 6), (1.2425070253013797, -6),
    (1.458652286180043, 6), (1.4587910768820167, -6), (1.471893283580752, 6), (1.4719891811363988, -6),
    (1.4954197314717812, -6), (1.4999604120353027, 4), (1.5112322337492317, 6), (1.5137198267248655, 6),
    (1.5833634316892131, -6), (1.5835301617498718, 6), (1.5975758718601678, -6), (1.5978944679397817, 6),
    (1.5979310022762243, -6), (1.6142165125456545, -6), (1.618757193109176, 4), (1.6223335080217338, 4),
    (1.9421323325589928, 6), (2.206362161738162, -6), (2.206631061331541, 6), (2.21230763019634, -6),
    (2.2124869837961363, 6), (2.2241252818417045, -6), (2.224185179845221, 6), (2.234375131178385, -6),
    (2.2343898319129853, 6), (2.2822472009359642, 6), (2.2941961583564523, -6), (2.2942798056742917, 6),
    (2.2946201062320015, -6), (2.3094436909675613, -6), (2.3177637004816867, 6), (2.317829482084852, 6),
    (2.318003901411435, -6), (2.4405659511142415, 6), (2.440819215880976, -6), (2.446298250659225, -6),
    (2.455815794247279, 6), (2.4560633544222683, -6), (2.4632083934155267, 6), (2.4632964665433295, 6),
    (2.487824857374012, -6), (2.4880847958046446, 6), (2.4926237381219223, -6), (2.492648484352202, 6),
    (2.4929573632928537, -6), (2.497615822300268, -6), (2.4996904584162896, 6), (2.5, 0),
])
WIT_G = (1.3450000000000002, [
    (1.118033988749895, -6), (1.118213161361215, 6), (1.1233499918526528, 6), (1.1339344425341358, -6),
    (1.1343295524038288, 6), (1.1346258060018437, -6), (1.1443050506249453, 6), (1.144561104985208, -6),
    (1.196733269327442, -6), (1.1970970973845356, 6), (1.2115084808056589, -6), (1.2115280246602382, 6),
    (1.2311071719479239, 6), (1.242166686301819, -6), (1.2424322268284247, 6), (1.4584780854844814, -6),
    (1.458652286180043, -6), (1.4586841135288593, 6), (1.471893283580752, 6), (1.4953976664944926, -6),
    (1.4954197314717812, 6), (1.5003704691902253, -6), (1.5137198267248655, -6), (1.5834788211896553, 6),
    (1.5835301617498718, -6), (1.597856566537614, 6), (1.5978944679397817, -6), (1.5979942893805366, 6),
    (1.6142165125456545, 6), (1.6219753953096594, -6), (1.6223335080217338, -6), (2.093631280581361, 6),
    (2.206362161738162, -6), (2.206529512544118, 6), (2.21230763019634, 6), (2.223715102595927, -6),
    (2.2241252818417045, 6), (2.2245124275799517, -6), (2.234375131178385, 6), (2.2344266172151372, -6),
    (2.2822472009359642, -6), (2.284255216871167, 4), (2.2942798056742917, 4), (2.2948744022446395, 6),
    (2.3094436909675613, 6), (2.317574089218916, -6), (2.317829482084852, -6), (2.3180104435223408, 6),
    (2.4405659511142415, -6), (2.4409367618601654, 6), (2.446298250659225, -6), (2.4466184275076466, 6),
    (2.4560633544222683, 6), (2.46323425172977, -6), (2.4632964665433295, -6), (2.4879278711077046, 6),
    (2.4880847958046446, -6), (2.492486980912437, 6), (2.492648484352202, -6), (2.492955030842477, 6),
    (2.497615822300268, 6), (2.49986240250181, -6), (2.5, 0),
])


if __name__ == "__main__":
    sys.exit(main())
