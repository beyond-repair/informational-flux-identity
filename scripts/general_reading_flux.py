#!/usr/bin/env python3
"""Theorem O witness: every two-derivative reading of (grad Psi)^{ij}.

T^{ij} = a d_i d_j Psi + b delta_ij Lap Psi + c d_i Psi d_j Psi + d delta_ij |grad Psi|^2
d_j T^{ij} = (a+b) d_i Lap Psi + c d_i Psi Lap Psi + (c/2 + d) d_i |grad Psi|^2.

Part 0 (linear algebra): O(3)-invariant rank-4 tensors form a 3-dimensional
space, so the family above is every shift-invariant, O(3)-covariant, local,
constant-coefficient polynomial tensor with exactly two derivatives.

Part 1 (exact, Fractions): the divergence identity holds as a polynomial
identity, and the closed-box flux equals the volume integral of the right side,
for seeded random integer (a, b, c, d) and cubic integer polynomials Psi.

Part 2 (quadrature, NumPy): static massless C6, Psi = sum q_k/|x - c_k|
(4 pi eps0 = 1) around an isolated asymmetric cluster. On every vacuum surface
the Hessian part and the Maxwell part give 0 and
G = (c/2 + d) oint |grad Psi|^2 n. That trace flux is nonzero, changes when
the surface is moved, and shrinks like R^-2 under dilation, with
R^2 P -> Q_tot^2 oint_{S1} n/|x|^4.

Exits nonzero if any identity or recorded figure moves. Nothing here is thrust.
"""
import itertools
import random
import sys
from fractions import Fraction

import numpy as np

FAIL = []


def check(ok, msg):
    print(("ok   " if ok else "FAIL ") + msg)
    if not ok:
        FAIL.append(msg)


# ---------------------------------------------------------------- Part 0 ----
def rot_matrix(axis, t):
    a = np.asarray(axis, float)
    a /= np.linalg.norm(a)
    K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K


gens = [rot_matrix((1, 2, 3), 0.7), rot_matrix((-2, 1, 0.5), 1.3), -np.eye(3),
        np.diag([1.0, 1.0, -1.0])]
rows = []
for g in gens:
    G4 = np.einsum("ia,jb,kc,ld->ijklabcd", g, g, g, g).reshape(81, 81)
    rows.append(G4 - np.eye(81))
sv = np.linalg.svd(np.vstack(rows), compute_uv=False)
null_dim = int(np.sum(sv < 1e-10))
check(null_dim == 3, f"O(3)-invariant rank-4 tensors: dimension {null_dim} "
      "(delta_ij delta_kl, delta_ik delta_jl, delta_il delta_jk)")


# ---------------------------------------------------------------- Part 1 ----
def padd(p, q, s=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + s * v
    return {k: v for k, v in r.items() if v != 0}


def pmul(p, q):
    r = {}
    for k1, v1 in p.items():
        for k2, v2 in q.items():
            k = (k1[0] + k2[0], k1[1] + k2[1], k1[2] + k2[2])
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
    return r


def pint(p, box, axes):
    r = {}
    for k, v in p.items():
        kk, vv = list(k), Fraction(v)
        for a in axes:
            lo, hi = box[a]
            vv *= Fraction(hi ** (k[a] + 1) - lo ** (k[a] + 1), k[a] + 1)
            kk[a] = 0
        r[tuple(kk)] = r.get(tuple(kk), 0) + vv
    return r


def peval_axis(p, a, t):
    r = {}
    for k, v in p.items():
        kk = list(k)
        kk[a] = 0
        r[tuple(kk)] = r.get(tuple(kk), 0) + v * Fraction(t) ** k[a]
    return {k: v for k, v in r.items() if v != 0}


def const(p):
    assert all(k == (0, 0, 0) for k in p), p
    return p.get((0, 0, 0), Fraction(0))


def lap(psi):
    r = {}
    for a in range(3):
        r = padd(r, pdiff(pdiff(psi, a), a))
    return r


def grad2(psi):
    r = {}
    for a in range(3):
        g = pdiff(psi, a)
        r = padd(r, pmul(g, g))
    return r


def T(psi, co, i, j):
    a, b, c, d = co
    gi, gj = pdiff(psi, i), pdiff(psi, j)
    t = padd(pscale(pdiff(gi, j), a), pscale(pmul(gi, gj), c))
    if i == j:
        t = padd(t, padd(pscale(lap(psi), b), pscale(grad2(psi), d)))
    return t


def div_formula(psi, co, i):
    a, b, c, d = co
    L = lap(psi)
    r = pscale(pdiff(L, i), a + b)
    r = padd(r, pscale(pmul(pdiff(psi, i), L), c))
    return padd(r, pscale(pdiff(grad2(psi), i), Fraction(c, 2) + d))


BOX = [(0, 2), (-1, 1), (0, 3)]
rng = random.Random(20261008 + 15)
monos = [(x, y, z) for x in range(4) for y in range(4) for z in range(4) if x + y + z <= 3]
n_pt = n_box = 0
TRIALS = 50
for _ in range(TRIALS):
    co = tuple(Fraction(rng.randint(-4, 4)) for _ in range(4))
    psi = {m: Fraction(rng.randint(-5, 5)) for m in rng.sample(monos, 7)}
    psi = {k: v for k, v in psi.items() if v}
    ok_pt = ok_box = True
    for i in range(3):
        div = {}
        for j in range(3):
            div = padd(div, pdiff(T(psi, co, i, j), j))
        f = div_formula(psi, co, i)
        ok_pt &= (div == f)
        flux = Fraction(0)
        for j in range(3):
            others = [a for a in range(3) if a != j]
            lo, hi = BOX[j]
            tij = T(psi, co, i, j)
            flux += const(pint(peval_axis(tij, j, hi), BOX, others))
            flux -= const(pint(peval_axis(tij, j, lo), BOX, others))
        ok_box &= (flux == const(pint(f, BOX, [0, 1, 2])))
    n_pt += ok_pt
    n_box += ok_box
check(n_pt == TRIALS, f"d_j T^ij = (a+b) d_i Lap + c d_i Psi Lap + (c/2+d) d_i|grad|^2 as polynomials: {n_pt}/{TRIALS}")
check(n_box == TRIALS, f"oint T n = int of that on box {BOX}: {n_box}/{TRIALS}")

# ---------------------------------------------------------------- Part 2 ----
NT, NP = 128, 256
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
S = np.sqrt(1 - U ** 2)
NRM = np.stack([S * np.cos(PH), S * np.sin(PH), U], axis=-1)
W = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]

cluster = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]
QTOT = sum(q for q, _ in cluster)


def fields(pts):
    g = np.zeros(pts.shape)
    H = np.zeros(pts.shape + (3,))
    for q, c in cluster:
        dd = pts - np.asarray(c)
        r = np.linalg.norm(dd, axis=-1)
        g += (-q / r ** 3)[..., None] * dd
        H += q * (3 * dd[..., :, None] * dd[..., None, :] / r[..., None, None] ** 5
                  - np.eye(3) / r[..., None, None] ** 3)
    return g, H


def surface(center, R):
    """Return the Hessian flux, Maxwell flux, trace flux P, and abs Maxwell row-x flux."""
    pts = np.asarray(center) + R * NRM
    g, H = fields(pts)
    gn = np.sum(g * NRM, axis=-1)
    g2 = np.sum(g * g, axis=-1)
    dA = R ** 2 * W
    hess = np.einsum("ij,ijk->k", dA, np.einsum("ijkl,ijl->ijk", H, NRM))
    lapl = np.einsum("ij,ij->", dA, np.abs(np.trace(H, axis1=-2, axis2=-1)))
    Mn = g * gn[..., None] - 0.5 * g2[..., None] * NRM
    M = np.einsum("ij,ijk->k", dA, Mn)
    P = np.einsum("ij,ijk->k", dA, g2[..., None] * NRM)
    absP = np.einsum("ij,ij->", dA, g2)
    return hess, lapl, M, P, absP


def G_of(co, hess, M, P):
    a, b, c, d = co
    # a*Hessian flux; b*oint Lap n (= 0 in vacuum); c*(Maxwell + P/2); d*P
    return a * hess + c * (M + 0.5 * P) + d * P


READINGS = {
    "Hessian (1,0,0,0)": (1, 0, 0, 0),
    "Hessian - delta Lap (1,-1,0,0)": (1, -1, 0, 0),
    "quadratic Q (0,0,1,-1/2)": (0, 0, 1, -0.5),
    "grad x grad (0,0,1,0)": (0, 0, 1, 0),
    "trace-free grad x grad (0,0,1,-1/3)": (0, 0, 1, -1 / 3),
    "delta |grad|^2 (0,0,0,1)": (0, 0, 0, 1),
}

SURF = [("sphere c=0 R=2", (0, 0, 0), 2.0),
        ("sphere c=(0.6,0.2,-0.3) R=2.4", (0.6, 0.2, -0.3), 2.4),
        ("sphere c=(-0.9,0,0.5) R=3", (-0.9, 0.0, 0.5), 3.0)]
res = {}
all_readings_ok = True
for name, cen, R in SURF:
    hess, lapl, M, P, absP = surface(cen, R)
    res[name] = P
    check(np.max(np.abs(hess)) < 1e-11 * absP and lapl < 1e-9 * absP,
          f"{name}: Hessian flux {np.max(np.abs(hess)):.1e}, oint|Lap Psi| {lapl:.1e} (vacuum)")
    check(np.max(np.abs(M)) < 1e-11 * absP,
          f"{name}: Maxwell flux {np.max(np.abs(M)):.1e} vs oint|grad Psi|^2 = {absP:.6f} (Coulomb self-force 0)")
    for rn, co in READINGS.items():
        Gv = G_of(co, hess, M, P)
        k = co[2] / 2 + co[3]
        ok = np.max(np.abs(Gv - k * P)) < 1e-10 * absP
        if not ok:
            all_readings_ok = False
            check(False, f"{name}, {rn}: G != (c/2+d) P")
    print(f"     {name}: P = oint |grad Psi|^2 n = {np.round(P, 9).tolist()}")

check(all_readings_ok, f"all {len(READINGS)} listed readings: G = (c/2+d) oint |grad Psi|^2 n on all three surfaces")

REC_P = {"sphere c=0 R=2": (-18.295274824, -6.781652669, 6.231252988),
         "sphere c=(0.6,0.2,-0.3) R=2.4": (-35.04687982, -11.166762011, 14.456027559),
         "sphere c=(-0.9,0,0.5) R=3": (4.23976535, -2.825807999, -3.116738754)}
check(all(np.max(np.abs(res[n] - np.array(v))) < 1e-8 for n, v in REC_P.items()),
      "recorded trace fluxes P on the three surfaces unchanged (1e-8)")
P0, P1, P2 = (res[n] for n, _, _ in SURF)
spread = max(np.linalg.norm(P0 - P1), np.linalg.norm(P0 - P2), np.linalg.norm(P1 - P2))
check(spread > 1e-2 * np.linalg.norm(P0),
      f"surface dependence: |P(S_a) - P(S_b)| up to {spread:.6f} with |P(S_0)| = {np.linalg.norm(P0):.6f}")

# Dilation S_R = R * S_1, S_1 = sphere center (0.3,-0.2,0.1), radius 1.
c1, r1 = np.array([0.3, -0.2, 0.1]), 1.0
pts1 = c1 + r1 * NRM
lim = QTOT ** 2 * np.einsum("ij,ijk->k", r1 ** 2 * W, NRM / np.linalg.norm(pts1, axis=-1)[..., None] ** 4)
scaled = []
for R in (4.0, 8.0, 16.0, 32.0, 64.0):
    *_, P, _ = surface(R * c1, R * r1)
    scaled.append((R, P, R ** 2 * P))
for R, P, s in scaled:
    print(f"     R={R:>4}: |P| = {np.linalg.norm(P):.6e}, R^2 P = {np.round(s, 6).tolist()}")
errs = [np.linalg.norm(s - lim) for _, _, s in scaled]
rich = 2 * scaled[-1][2] - scaled[-2][2]
rel = np.linalg.norm(rich - lim) / np.linalg.norm(lim)
check(all(e2 < 0.55 * e1 for e1, e2 in zip(errs, errs[1:])) and rel < 5e-3,
      f"dilation: R^2 P -> Q_tot^2 oint_S1 n/|x|^4 = {np.round(lim, 6).tolist()}; errors "
      + ", ".join(f"{e:.2e}" for e in errs) + f" (halving, O(1/R)); Richardson 2s(64)-s(32) off by {rel:.1e}")
check(np.linalg.norm(lim) > 0, "limit coefficient nonzero: G(R S_1) = (c/2+d) R^-2 Q_tot^2 oint_S1 n/|x|^4 + O(R^-3), nonzero at every R, -> 0")

if FAIL:
    print(f"{len(FAIL)} check(s) failed")
    sys.exit(1)
print("all Theorem O checks passed; none of this is thrust")
