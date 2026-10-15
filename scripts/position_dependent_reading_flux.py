#!/usr/bin/env python3
"""Theorem Y witness: position-dependent first-derivative isotropic readings.

Family (Assumption Y):
    T^{ij} = f(x, Psi, s) delta^{ij} + g(x, Psi, s) d^i Psi d^j Psi,  s = |grad Psi|^2,
with f, g in C^1 on an open set of R^3 x R x [0, infinity). Theorem V is the
x-independent slice. Same C6 / static massless / isolated device / vacuum
surfaces as Assumption Q/U/V.

On-shell (Lap Psi = 0 in vacuum):
    d_j T^{ij} = df/dx_i|_explicit + (grad_x g · E) E^i
               + (f_Psi + s g_Psi) E^i + (f_s + g/2) d_i s
               + g_s E^i (E · grad s) + g E^i Lap Psi.
Conservation for every harmonic Psi forces the V conditions pointwise in x,
then g = g(x), f = -(g(x)/2) s + h(x), then
    grad h - (1/2)(grad g) |E|^2 + (grad g · E) E ≡ 0 for all E
forces grad h = 0 and grad g = 0, so
    T = kappa * delta + c * Q
with constant kappa, c — same conserved class as Theorems U/V. Position
dependence does not enlarge it. Conserved readings give G = 0 on every vacuum
surface around an isolated device; otherwise G(S) = -int_outside div T when
the surface term at infinity vanishes, and in every case G moves with S.

Special cases checked explicitly: a(x) Q, a(x) delta, a(x) Q + b(x) delta,
and non-polynomial a(x) = 1/(1+|x|^2), exp(-|x|^2). For linear a(x)=x the
exterior-at-infinity integral of div(a delta)=e_x diverges; the witness uses
the exact volume identity G = Vol(V) e_x instead.

Part 0 (exact, Fractions): divergence identity for polynomial a(x) times
Maxwell / delta; conservation of constant Maxwell+kappa; non-conservation
when grad a ≠ 0.
Part 1 (numeric FD): non-polynomial a(x) and mixed soft(x)*Psi delta.
Part 2 (quadrature): Maxwell G ~ 0; soft/gauss/linear position-dependent
readings give different G on three vacuum spheres; exterior match where the
integrand decays; volume identity for x delta; dilation: soft a Q and gauss
decay, soft a delta has a nonzero surface-dependent limit (R.4 sibling),
x Q decays slowly. Maxwell removal leaves the position-dependent remainder.

Exits nonzero if any identity or recorded figure moves. Nothing here is thrust.
Relates to Theorem E only as a sibling: E is the C5 prefactor W(x); Y is the
reading's own coefficients.
"""
from __future__ import annotations

import math
import random
import sys
from fractions import Fraction

import numpy as np

FAIL = []


def check(ok, msg):
    print(("ok   " if ok else "FAIL ") + msg)
    if not ok:
        FAIL.append(msg)


# ------------------------------------------------- exact polynomial helpers --
def padd(p, q, s=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + s * v
    return {k: v for k, v in r.items() if v != 0}


def pmul(p, q):
    r = {}
    for k1, v1 in p.items():
        for k2, v2 in q.items():
            k = tuple(a + b for a, b in zip(k1, k2))
            r[k] = r.get(k, 0) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}


def pscale(p, s):
    return {k: v * s for k, v in p.items() if v * s != 0}


def pdiff(p, i):
    r = {}
    for k, v in p.items():
        if k[i]:
            kk = list(k)
            kk[i] -= 1
            r[tuple(kk)] = r.get(tuple(kk), 0) + v * k[i]
    return {k: v for k, v in r.items() if v != 0}


def peval(p, xyz):
    t = Fraction(0)
    for k, v in p.items():
        m = v
        for a, e in enumerate(k):
            m *= Fraction(xyz[a]) ** e
        t += m
    return t


def lap(p):
    r = {}
    for a in range(3):
        r = padd(r, pdiff(pdiff(p, a), a))
    return r


def grad2(psi):
    r = {}
    for a in range(3):
        g = pdiff(psi, a)
        r = padd(r, pmul(g, g))
    return r


X = {(1, 0, 0): Fraction(1)}
Y = {(0, 1, 0): Fraction(1)}
ONE = {(0, 0, 0): Fraction(1)}


def T_scaled_Q(psi, a, i, j):
    """a(x) * Q^{ij} with Q = E^i E^j - (1/2) s delta^{ij}."""
    s = grad2(psi)
    gi, gj = pdiff(psi, i), pdiff(psi, j)
    t = pmul(a, pmul(gi, gj))
    if i == j:
        t = padd(t, pmul(a, pscale(s, Fraction(-1, 2))))
    return t


def T_scaled_delta(a, i, j):
    if i == j:
        return dict(a)
    return {}


def T_axQ_bx_delta(psi, a, b, i, j):
    return padd(T_scaled_Q(psi, a, i, j), T_scaled_delta(b, i, j))


def div_direct(T_fn, i):
    d = {}
    for j in range(3):
        d = padd(d, pdiff(T_fn(i, j), j))
    return d


def div_aQ_formula(psi, a, i):
    """(E·grad a) E^i - (s/2) d_i a + a E^i Lap."""
    s = grad2(psi)
    E_i = pdiff(psi, i)
    da = pdiff(a, i)
    Edota = {}
    for j in range(3):
        Edota = padd(Edota, pmul(pdiff(psi, j), pdiff(a, j)))
    r = pmul(Edota, E_i)
    r = padd(r, pscale(pmul(s, da), Fraction(-1, 2)))
    r = padd(r, pmul(a, pmul(E_i, lap(psi))))
    return r


def div_a_delta_formula(a, i):
    return pdiff(a, i)


print("=== Part 0: exact divergence identity (polynomial a(x) times Q / delta) ===")
rng = random.Random(20261009)

PSI_H = {(3, 0, 0): 1, (1, 2, 0): -3, (0, 1, 1): 2, (1, 1, 1): 1,
         (4, 0, 0): 1, (2, 2, 0): -6, (0, 4, 0): 1}
check(lap(PSI_H) == {}, "Psi_h of Theorem Q is harmonic")

monos3 = [(x, y, z) for x in range(4) for y in range(4) for z in range(4) if x + y + z <= 3]
monos2 = [(x, y, z) for x in range(3) for y in range(3) for z in range(3) if x + y + z <= 2]

TRIALS = 40
n_ok_Q = n_ok_d = n_ok_mix = 0
for _ in range(TRIALS):
    psi = {m: Fraction(rng.randint(-3, 3)) for m in rng.sample(monos3, 5)}
    psi = {k: v for k, v in psi.items() if v}
    a = {m: Fraction(rng.randint(-3, 3)) for m in rng.sample(monos2, 4)}
    a = {k: v for k, v in a.items() if v} or dict(ONE)
    b = {m: Fraction(rng.randint(-2, 2)) for m in rng.sample(monos2, 3)}
    b = {k: v for k, v in b.items() if v} or dict(ONE)
    n_ok_Q += all(
        div_direct(lambda i, j, _a=a, _p=psi: T_scaled_Q(_p, _a, i, j), i)
        == div_aQ_formula(psi, a, i) for i in range(3))
    n_ok_d += all(
        div_direct(lambda i, j, _a=a: T_scaled_delta(_a, i, j), i)
        == div_a_delta_formula(a, i) for i in range(3))
    n_ok_mix += all(
        div_direct(lambda i, j, _a=a, _b=b, _p=psi: T_axQ_bx_delta(_p, _a, _b, i, j), i)
        == padd(div_aQ_formula(psi, a, i), div_a_delta_formula(b, i))
        for i in range(3))
check(n_ok_Q == TRIALS, f"div(a Q) = (E·∇a)E - (s/2)∇a + a E Lap: {n_ok_Q}/{TRIALS}")
check(n_ok_d == TRIALS, f"div(a δ) = ∇a: {n_ok_d}/{TRIALS}")
check(n_ok_mix == TRIALS, f"div(a Q + b δ) = div(aQ)+∇b: {n_ok_mix}/{TRIALS}")

cons_ok = 0
for k in range(-2, 3):
    for c in range(-2, 3):
        a = pscale(ONE, Fraction(c)) if c else {}
        b = pscale(ONE, Fraction(k)) if k else {}
        ok = all(
            div_direct(lambda i, j, _a=a, _b=b: T_axQ_bx_delta(PSI_H, _a, _b, i, j), i) == {}
            for i in range(3))
        cons_ok += ok
check(cons_ok == 25, f"all 25 affine constant Maxwell+kappa readings conserved on Psi_h: {cons_ok}/25")

a_x = dict(X)
div_aQ = [div_direct(lambda i, j: T_scaled_Q(PSI_H, a_x, i, j), i) for i in range(3)]
check(any(d != {} for d in div_aQ), "a(x)=x times Q: div does not vanish on Psi_h")
check(div_aQ == [div_aQ_formula(PSI_H, a_x, i) for i in range(3)],
      "a(x)=x times Q: on-shell div = (E·e_x)E - (s/2)e_x on Psi_h")

div_ad = [div_direct(lambda i, j: T_scaled_delta(a_x, i, j), i) for i in range(3)]
check(div_ad == [ONE if i == 0 else {} for i in range(3)],
      "a(x)=x times δ: on-shell div = e_x exactly")

a_x, b_y = dict(X), dict(Y)
div_mix = [div_direct(lambda i, j: T_axQ_bx_delta(PSI_H, a_x, b_y, i, j), i) for i in range(3)]
expected_mix = [padd(div_aQ_formula(PSI_H, a_x, i), div_a_delta_formula(b_y, i)) for i in range(3)]
check(div_mix == expected_mix, "a=x Q + b=y δ: div formula exact on Psi_h")
check(any(d != {} for d in div_mix), "a=x Q + b=y δ: not on-shell conserved")

PT = (1, Fraction(1, 2), -1)
spot = tuple(peval(div_aQ[i], PT) for i in range(3))
print(f"     div(x Q) on Psi_h at (1,1/2,-1) = {tuple(str(v) for v in spot)}")
check(spot != (0, 0, 0), "a(x)=x Q is not on-shell conserved (spot eval)")

# Algebraic collapse of the all-E identity for T = h(x)δ + g(x)Q:
# grad h - (1/2)(grad g)|E|^2 + (grad g · E) E = 0 for all E.
# E=0 ⇒ grad h = 0; then E=grad g ⇒ (1/2)|grad g|^2 grad g = 0 ⇒ grad g = 0.
check(True, "conservation criterion collapses: T=h(x)δ+g(x)Q conserved for all E "
      "iff grad h=0 and grad g=0 (constants only)")


# -------------------------------- Part 1: FD check for non-polynomial a(x) --
print("\n=== Part 1: finite-difference divergence for non-polynomial a(x) ===")


def num_grad_hess_lap(fn, x, h=1e-5):
    g = np.zeros(3)
    H = np.zeros((3, 3))
    for i in range(3):
        ep = np.zeros(3)
        ep[i] = h
        g[i] = (fn(x + ep) - fn(x - ep)) / (2 * h)
    for i in range(3):
        for j in range(3):
            ei, ej = np.zeros(3), np.zeros(3)
            ei[i] = h
            ej[j] = h
            H[i, j] = (fn(x + ei + ej) - fn(x + ei - ej)
                       - fn(x - ei + ej) + fn(x - ei - ej)) / (4 * h * h)
    return g, H, np.trace(H)


def T_val(x, Psi, E, f, g):
    s = float(np.dot(E, E))
    return f(x, Psi, s) * np.eye(3) + g(x, Psi, s) * np.outer(E, E)


def div_T_fd(fn, x, f, g, h=1e-5):
    out = np.zeros(3)
    for j in range(3):
        ep = np.zeros(3)
        ep[j] = h
        Ep, _, _ = num_grad_hess_lap(fn, x + ep, h=h * 0.5)
        Em, _, _ = num_grad_hess_lap(fn, x - ep, h=h * 0.5)
        Pp, Pm = fn(x + ep), fn(x - ep)
        Tp = T_val(x + ep, Pp, Ep, f, g)
        Tm = T_val(x - ep, Pm, Em, f, g)
        out += (Tp[:, j] - Tm[:, j]) / (2 * h)
    return out


def div_T_formula_num(x, Psi, E, H, Lap, f, g, f_x, f_psi, f_s, g_x, g_psi, g_s):
    s = float(np.dot(E, E))
    ds = 2 * H @ E
    return (f_x(x, Psi, s)
            + np.dot(g_x(x, Psi, s), E) * E
            + (f_psi(x, Psi, s) + s * g_psi(x, Psi, s)) * E
            + (f_s(x, Psi, s) + g(x, Psi, s) / 2) * ds
            + g_s(x, Psi, s) * np.dot(E, ds) * E
            + g(x, Psi, s) * E * Lap)


def harm_fn(x):
    return float(x[0] ** 3 - 3 * x[0] * x[1] ** 2 + 2 * x[1] * x[2])


def nonharm_fn(x):
    return float(x[0] ** 2 + 0.5 * x[1] * x[2])


def soft_a(x):
    return 1.0 / (1.0 + float(np.dot(x, x)))


def soft_grad(x):
    r2 = float(np.dot(x, x))
    return -2.0 * x / (1.0 + r2) ** 2


def gauss_a(x):
    return math.exp(-float(np.dot(x, x)))


def gauss_grad(x):
    return -2.0 * x * math.exp(-float(np.dot(x, x)))


NP_PAIRS = [
    ("soft a(x) Q",
     lambda x, P, s: -0.5 * soft_a(x) * s, lambda x, P, s: soft_a(x),
     lambda x, P, s: -0.5 * soft_grad(x) * s, lambda x, P, s: 0.0,
     lambda x, P, s: -0.5 * soft_a(x), lambda x, P, s: soft_grad(x),
     lambda x, P, s: 0.0, lambda x, P, s: 0.0),
    ("soft a(x) δ",
     lambda x, P, s: soft_a(x), lambda x, P, s: 0.0,
     lambda x, P, s: soft_grad(x), lambda x, P, s: 0.0,
     lambda x, P, s: 0.0, lambda x, P, s: np.zeros(3),
     lambda x, P, s: 0.0, lambda x, P, s: 0.0),
    ("gauss a(x) Q",
     lambda x, P, s: -0.5 * gauss_a(x) * s, lambda x, P, s: gauss_a(x),
     lambda x, P, s: -0.5 * gauss_grad(x) * s, lambda x, P, s: 0.0,
     lambda x, P, s: -0.5 * gauss_a(x), lambda x, P, s: gauss_grad(x),
     lambda x, P, s: 0.0, lambda x, P, s: 0.0),
    ("gauss a(x) δ",
     lambda x, P, s: gauss_a(x), lambda x, P, s: 0.0,
     lambda x, P, s: gauss_grad(x), lambda x, P, s: 0.0,
     lambda x, P, s: 0.0, lambda x, P, s: np.zeros(3),
     lambda x, P, s: 0.0, lambda x, P, s: 0.0),
    ("soft(x)·Psi δ",
     lambda x, P, s: soft_a(x) * P, lambda x, P, s: 0.0,
     lambda x, P, s: soft_grad(x) * P, lambda x, P, s: soft_a(x),
     lambda x, P, s: 0.0, lambda x, P, s: np.zeros(3),
     lambda x, P, s: 0.0, lambda x, P, s: 0.0),
]

pts = [np.array([0.7, -0.3, 0.4]), np.array([1.2, 0.5, -0.8]), np.array([-0.4, 0.9, 0.2])]
for name, f, g, f_x, f_psi, f_s, g_x, g_psi, g_s in NP_PAIRS:
    max_rel = 0.0
    for fn in (harm_fn, nonharm_fn):
        for x in pts:
            E, H, Lap = num_grad_hess_lap(fn, x)
            Psi = fn(x)
            fd = div_T_fd(fn, x, f, g)
            form = div_T_formula_num(x, Psi, E, H, Lap, f, g,
                                     f_x, f_psi, f_s, g_x, g_psi, g_s)
            scale = max(1e-8, np.linalg.norm(form), np.linalg.norm(fd))
            max_rel = max(max_rel, np.linalg.norm(fd - form) / scale)
    check(max_rel < 5e-4,
          f"{name}: FD div matches formula (harm+nonharm; max rel {max_rel:.2e})")

x0 = pts[0]
E, H, Lap = num_grad_hess_lap(harm_fn, x0)
Psi = harm_fn(x0)
nd = float(np.linalg.norm(div_T_formula_num(x0, Psi, E, H, Lap, *NP_PAIRS[0][1:])))
check(nd > 1e-3, "soft a(x) Q not on-shell conserved (|div|=%.4f at test point)" % nd)
ndg = float(np.linalg.norm(div_T_formula_num(x0, Psi, E, H, Lap, *NP_PAIRS[2][1:])))
check(ndg > 1e-3, "gauss a(x) Q not on-shell conserved (|div|=%.4f)" % ndg)


# ----------------------------------------------------- Part 2: quadrature ----
print("\n=== Part 2: quadrature on Coulomb cluster (static massless C6) ===")
NT, NP = 128, 256
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
S = np.sqrt(1 - U ** 2)
NRM = np.stack([S * np.cos(PH), S * np.sin(PH), U], axis=-1)
Wgrid = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]

CLUSTER = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]


def fields(pts):
    P = np.zeros(pts.shape[:-1])
    g = np.zeros(pts.shape)
    H = np.zeros(pts.shape + (3,))
    for q, c in CLUSTER:
        dd = pts - np.asarray(c)
        r = np.linalg.norm(dd, axis=-1)
        P += q / r
        g += (-q / r ** 3)[..., None] * dd
        H += q * (3 * dd[..., :, None] * dd[..., None, :] / r[..., None, None] ** 5
                  - np.eye(3) / r[..., None, None] ** 3)
    return P, g, H


def flux_and_abs(center, R, f, g):
    pts = np.asarray(center) + R * NRM
    Psi, E, H = fields(pts)
    s = np.sum(E * E, axis=-1)
    En = np.sum(E * NRM, axis=-1)
    fv = f(pts, Psi, s)
    gv = g(pts, Psi, s)
    Tn = fv[..., None] * NRM + gv[..., None] * En[..., None] * E
    dA = R ** 2 * Wgrid
    G = np.einsum("ij,ijk->k", dA, Tn)
    abs_x = np.einsum("ij,ij->", dA, np.abs(Tn[..., 0]))
    return G, abs_x, Psi, E, H, s, pts


def onshell_div(pts, Psi, E, H, f, g, f_x, f_psi, f_s, g_x, g_psi, g_s):
    s = np.sum(E * E, axis=-1)
    ds = 2 * np.einsum("...ik,...k->...i", H, E)
    fx = f_x(pts, Psi, s)
    gx = g_x(pts, Psi, s)
    return (fx
            + np.sum(gx * E, axis=-1)[..., None] * E
            + (f_psi(pts, Psi, s) + s * g_psi(pts, Psi, s))[..., None] * E
            + (f_s(pts, Psi, s) + 0.5 * g(pts, Psi, s))[..., None] * ds
            + (g_s(pts, Psi, s) * np.sum(E * ds, axis=-1))[..., None] * E)


NR = 96
tr, wr = np.polynomial.legendre.leggauss(NR)
tr, wr = 0.5 * (tr + 1), 0.5 * wr


def exterior_int(center, R, f, g, f_x, f_psi, f_s, g_x, g_psi, g_s):
    out = np.zeros(3)
    for t, w in zip(tr, wr):
        r = R / t
        pts = np.asarray(center) + r * NRM
        Psi, E, H = fields(pts)
        d = onshell_div(pts, Psi, E, H, f, g, f_x, f_psi, f_s, g_x, g_psi, g_s)
        jac = r ** 2 * (R / t ** 2) * w
        out -= jac * np.einsum("ij,ijk->k", Wgrid, d)
    return out


def soft_a_arr(pts):
    return 1.0 / (1.0 + np.sum(pts * pts, axis=-1))


def soft_grad_arr(pts):
    r2 = np.sum(pts * pts, axis=-1)
    return -2.0 * pts / (1.0 + r2)[..., None] ** 2


def gauss_a_arr(pts):
    return np.exp(-np.sum(pts * pts, axis=-1))


def gauss_grad_arr(pts):
    return -2.0 * pts * np.exp(-np.sum(pts * pts, axis=-1))[..., None]


def zeros_s(pts, Psi, s):
    return np.zeros_like(s)


def zeros_vec(pts, Psi, s):
    return np.zeros(pts.shape)


def f_max(pts, Psi, s):
    return -0.5 * s


def g_max(pts, Psi, s):
    return np.ones_like(s)


def fx0(pts, Psi, s):
    return np.zeros(pts.shape)


def fs_max(pts, Psi, s):
    return -0.5 * np.ones_like(s)


def f_softQ(pts, Psi, s):
    return -0.5 * soft_a_arr(pts) * s


def g_softQ(pts, Psi, s):
    return soft_a_arr(pts)


def fx_softQ(pts, Psi, s):
    return -0.5 * soft_grad_arr(pts) * s[..., None]


def fs_softQ(pts, Psi, s):
    return -0.5 * soft_a_arr(pts)


def gx_softQ(pts, Psi, s):
    return soft_grad_arr(pts)


def f_softD(pts, Psi, s):
    return soft_a_arr(pts)


def g0(pts, Psi, s):
    return np.zeros_like(s)


def fx_softD(pts, Psi, s):
    return soft_grad_arr(pts)


def f_gaussQ(pts, Psi, s):
    return -0.5 * gauss_a_arr(pts) * s


def g_gaussQ(pts, Psi, s):
    return gauss_a_arr(pts)


def fx_gaussQ(pts, Psi, s):
    return -0.5 * gauss_grad_arr(pts) * s[..., None]


def fs_gaussQ(pts, Psi, s):
    return -0.5 * gauss_a_arr(pts)


def gx_gaussQ(pts, Psi, s):
    return gauss_grad_arr(pts)


def f_gaussD(pts, Psi, s):
    return gauss_a_arr(pts)


def fx_gaussD(pts, Psi, s):
    return gauss_grad_arr(pts)


def f_xQ(pts, Psi, s):
    return -0.5 * pts[..., 0] * s


def g_xQ(pts, Psi, s):
    return pts[..., 0]


def fx_xQ(pts, Psi, s):
    out = np.zeros(pts.shape)
    out[..., 0] = -0.5 * s
    return out


def fs_xQ(pts, Psi, s):
    return -0.5 * pts[..., 0]


def gx_xQ(pts, Psi, s):
    out = np.zeros(pts.shape)
    out[..., 0] = 1.0
    return out


def f_xD(pts, Psi, s):
    return pts[..., 0]


def fx_xD(pts, Psi, s):
    out = np.zeros(pts.shape)
    out[..., 0] = 1.0
    return out


# (name, f, g, f_x, f_psi, f_s, g_x, g_psi, g_s, conserved, exterior_ok)
# exterior_ok=False for x δ: div=e_x does not decay, exterior integral diverges.
READINGS = [
    ("Maxwell Q", f_max, g_max, fx0, zeros_s, fs_max, zeros_vec, zeros_s, zeros_s, True, True),
    ("soft a(x) Q", f_softQ, g_softQ, fx_softQ, zeros_s, fs_softQ, gx_softQ, zeros_s, zeros_s, False, True),
    ("soft a(x) δ", f_softD, g0, fx_softD, zeros_s, zeros_s, zeros_vec, zeros_s, zeros_s, False, True),
    ("gauss a(x) Q", f_gaussQ, g_gaussQ, fx_gaussQ, zeros_s, fs_gaussQ, gx_gaussQ, zeros_s, zeros_s, False, True),
    ("gauss a(x) δ", f_gaussD, g0, fx_gaussD, zeros_s, zeros_s, zeros_vec, zeros_s, zeros_s, False, True),
    ("x Q", f_xQ, g_xQ, fx_xQ, zeros_s, fs_xQ, gx_xQ, zeros_s, zeros_s, False, True),
    ("x δ", f_xD, g0, fx_xD, zeros_s, zeros_s, zeros_vec, zeros_s, zeros_s, False, False),
]

SURF = [("c=0 R=2", (0.0, 0.0, 0.0), 2.0),
        ("c=(0.6,0.2,-0.3) R=2.4", (0.6, 0.2, -0.3), 2.4),
        ("c=(-0.9,0,0.5) R=3", (-0.9, 0.0, 0.5), 3.0)]


def match_exterior(G, ext):
    """Absolute when both tiny (symmetry zeros); else relative."""
    abs_err = float(np.max(np.abs(G - ext)))
    scale = max(np.max(np.abs(G)), np.max(np.abs(ext)))
    if scale < 1e-8:
        return abs_err < 1e-8, abs_err, 0.0
    rel = abs_err / scale
    return rel < 1e-8, abs_err, rel


REC = {}
for rn, f, g, f_x, f_psi, f_s, g_x, g_psi, g_s, conserved, exterior_ok in READINGS:
    vals = []
    for name, cen, R in SURF:
        G, abs_x, Psi, E, H, s, pts = flux_and_abs(cen, R, f, g)
        print(f"     {rn} sphere {name}: G = {np.array2string(G, precision=10)}  "
              f"(abs row-x {abs_x:.6f})", end="")
        if conserved:
            print()
            check(np.max(np.abs(G)) < 1e-11 and abs_x > 1,
                  f"{rn} sphere {name}: |G| = {np.max(np.abs(G)):.1e} vs abs row-x {abs_x:.6f}")
        elif not exterior_ok:
            print("  [no exterior: div does not decay]")
            vals.append(tuple(float(x) for x in G))
        else:
            ext = exterior_int(cen, R, f, g, f_x, f_psi, f_s, g_x, g_psi, g_s)
            ok, abs_err, rel = match_exterior(G, ext)
            print(f"  (vs exterior abs {abs_err:.1e} rel {rel:.1e})")
            check(ok, f"{rn} sphere {name}: G = -int_outside div "
                  f"(abs {abs_err:.1e}, rel {rel:.1e})")
            vals.append(tuple(float(x) for x in G))
    if not conserved:
        REC[rn] = vals
        gap = min(np.max(np.abs(np.array(vals[a]) - np.array(vals[b])))
                  for a in range(3) for b in range(a + 1, 3))
        check(gap > 0.01, f"{rn}: same device, three spheres, three values "
              f"(smallest pairwise gap {gap:.4f})")

# x δ: exact volume identity G = Vol(V) e_x = (4π/3) R^3 e_x for any sphere
print("\n=== x δ volume identity (exterior diverges; use G = ∫_V ∇a) ===")
for s_idx, (name, cen, R) in enumerate(SURF):
    G = np.array(REC["x δ"][s_idx])
    expected = np.array([(4.0 * math.pi / 3.0) * R ** 3, 0.0, 0.0])
    err = np.max(np.abs(G - expected))
    check(err < 1e-8, f"x δ sphere {name}: G = (4π/3)R³ e_x (err {err:.1e})")
    print(f"     {name}: G_x = {G[0]:.10f}, (4π/3)R³ = {expected[0]:.10f}")

print("\n=== Maxwell removal leaves the position-dependent remainder ===")
for name, cen, R in SURF:
    Gq, _, _, _, _, _, _ = flux_and_abs(cen, R, f_max, g_max)
    Gs, _, _, _, _, _, _ = flux_and_abs(cen, R, f_softD, g0)
    Gsum, _, _, _, _, _, _ = flux_and_abs(
        cen, R,
        lambda pts, P, s: f_max(pts, P, s) + f_softD(pts, P, s),
        g_max)
    check(np.max(np.abs(Gsum - (Gq + Gs))) < 1e-12, f"linearity of flux on {name}")
    check(np.max(np.abs(Gsum - Gs)) < 1e-11,
          f"after Maxwell removal on {name}: residual = soft a(x) δ flux "
          f"(max |G_Q|={np.max(np.abs(Gq)):.1e})")

print("\n=== Dilation: soft/gauss Q decay; soft δ -> nonzero surface-dependent limit ===")
C1 = np.array((0.3, -0.2, 0.1))
hat = C1 / np.linalg.norm(C1)
mags_softQ, mags_gaussQ, mags_gaussD, mags_xQ = [], [], [], []
softD_vecs = []
for R in (4.0, 8.0, 16.0, 32.0, 64.0):
    Gsq, _, _, _, _, _, _ = flux_and_abs(R * C1, R, f_softQ, g_softQ)
    Gsd, _, _, _, _, _, _ = flux_and_abs(R * C1, R, f_softD, g0)
    Ggq, _, _, _, _, _, _ = flux_and_abs(R * C1, R, f_gaussQ, g_gaussQ)
    Ggd, _, _, _, _, _, _ = flux_and_abs(R * C1, R, f_gaussD, g0)
    Gxq, _, _, _, _, _, _ = flux_and_abs(R * C1, R, f_xQ, g_xQ)
    mags_softQ.append(float(np.linalg.norm(Gsq)))
    softD_vecs.append(Gsd.copy())
    mags_gaussQ.append(float(np.linalg.norm(Ggq)))
    mags_gaussD.append(float(np.linalg.norm(Ggd)))
    mags_xQ.append(float(np.linalg.norm(Gxq)))
    print(f"     R={R:>4}: |G_softQ|={mags_softQ[-1]:.6e}, "
          f"G_softD={np.array2string(Gsd, precision=6)}, "
          f"|G_gaussQ|={mags_gaussQ[-1]:.6e}, |G_gaussD|={mags_gaussD[-1]:.6e}, "
          f"|G_xQ|={mags_xQ[-1]:.6e}")

check(all(mags_softQ[k + 1] < 0.7 * mags_softQ[k] for k in range(1, 4)),
      "dilation: |G_soft a Q| decreases for R>=8")
check(mags_softQ[-1] < 0.05,
      f"dilation: |G_soft a Q|(R=64) = {mags_softQ[-1]:.3e} -> 0")
check(mags_gaussQ[-1] < 1e-6 and mags_gaussD[-1] < 1e-6,
      f"dilation: gauss a readings at R=64 ~0 "
      f"(|G_Q|={mags_gaussQ[-1]:.2e}, |G_δ|={mags_gaussD[-1]:.2e})")
check(mags_xQ[-1] < mags_xQ[0],
      f"dilation: |G_x Q|(R=64)={mags_xQ[-1]:.3e} < |G|(R=4)={mags_xQ[0]:.3e}")

# soft a δ: nonzero surface-dependent limit, direction -Ĉ1 (R.4 sibling)
lim = softD_vecs[-1]
check(np.linalg.norm(lim) > 1.0,
      f"dilation: soft a δ has nonzero limit (|G|(R=64)={np.linalg.norm(lim):.4f})")
align = float(np.dot(lim, -hat) / np.linalg.norm(lim))
check(align > 0.999,
      f"dilation: soft a δ limit aligns with -Ĉ1 (cos={align:.6f})")
# Stabilizes: last two magnitudes within 1%
m64, m32 = np.linalg.norm(softD_vecs[-1]), np.linalg.norm(softD_vecs[-2])
check(abs(m64 - m32) / m64 < 0.01,
      f"dilation: soft a δ |G| stabilizes (R=32->{m32:.4f}, R=64->{m64:.4f})")
# Different offset direction ⇒ different limit direction (surface-dependent)
C2 = np.array((-0.2, 0.3, 0.1))
G_other, _, _, _, _, _, _ = flux_and_abs(64.0 * C2, 64.0, f_softD, g0)
check(np.linalg.norm(G_other - lim) > 1.0,
      f"dilation: soft a δ limit depends on offset direction "
      f"(|G_C1 - G_C2|={np.linalg.norm(G_other - lim):.3f})")

LOCKED = {
    "soft a(x) Q": [
        (-0.0000000000, -0.0000000000, 0.0000000000),
        (-0.6364258354, -0.2102975036, 0.3359147973),
        (0.2414780083, 0.0154592946, -0.1309581046),
    ],
    "soft a(x) δ": [
        (-0.0000000000, -0.0000000000, -0.0000000000),
        (-1.5237956654, -0.5079318885, 0.7618978327),
        (2.0678503489, 0.0000000000, -1.1488057494),
    ],
    "gauss a(x) Q": [
        (-0.0000000000, -0.0000000000, -0.0000000000),
        (-0.2716189600, -0.0922300380, 0.1508885882),
        (0.0196246838, 0.0046628728, -0.0104428429),
    ],
    "gauss a(x) δ": [
        (-0.0000000000, 0.0000000000, 0.0000000000),
        (-0.3611976657, -0.1203992219, 0.1805988329),
        (0.1381385175, 0.0000000000, -0.0767436208),
    ],
    "x Q": [
        (17.5779417342, -0.1917692645, 0.1652142724),
        (16.3914344267, -0.6150475666, 0.7051558669),
        (11.4162693907, 0.0429072747, 0.0671943269),
    ],
    "x δ": [
        (33.5103216383, 0.0000000000, -0.0000000000),
        (57.9058357910, -0.0000000000, -0.0000000000),
        (113.0973355292, -0.0000000000, -0.0000000000),
    ],
}

print("\n=== Recorded non-conserved fluxes (lock) ===")
for rn, vals in REC.items():
    for s_idx, (got, exp) in enumerate(zip(vals, LOCKED[rn])):
        err = np.max(np.abs(np.array(got) - np.array(exp)))
        check(err < 1e-8, f"{rn} sphere[{s_idx}]: recorded value (err {err:.1e})")
        print(f"     {rn}[{s_idx}] = {np.array2string(np.array(got), precision=10)}")

if FAIL:
    print(f"\n{len(FAIL)} check(s) failed")
    sys.exit(1)
print("\nall Theorem Y checks passed; nothing here is thrust")
