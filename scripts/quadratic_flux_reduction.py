#!/usr/bin/env python3
"""Theorem N witness: the quadratic reading of (grad Psi)^{ij}.

Q^{ij} = d_i Psi d_j Psi - (1/2) delta_ij |grad Psi|^2,  d_j Q^{ij} = d_i Psi * Lap Psi.

Part 1 (exact, Fractions): on a 3D box, the closed-surface flux of Q equals the
volume integral of d_i Psi Lap Psi for seeded random integer polynomials, and a
harmonic polynomial gives G = 0 exactly while individual faces carry nonzero flux.

Part 2 (quadrature, NumPy): Psi = sum q_k / |x - c_k| (static massless C6 with
4*pi*eps0 = 1, so Lap Psi = -4 pi rho). On a sphere enclosing an asymmetric
cluster, G = 0 (Coulomb self-force vanishes). With one charge outside, G equals
4 pi times the ordinary Coulomb force on the enclosed charges. Yukawa
(screened) potentials leave G = (m^2/2) oint Psi^2 n, which shrinks as the
sphere grows.

Exits nonzero if any identity or recorded figure moves. Nothing here is thrust.
"""
import random
import sys
from fractions import Fraction

import numpy as np

FAIL = []


def check(ok, msg):
    print(("ok   " if ok else "FAIL ") + msg)
    if not ok:
        FAIL.append(msg)


# ---------------------------------------------------------------- Part 1 ----
# Polynomials in (x, y, z): dict {(a, b, c): Fraction}.
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
    """Integrate over the listed axes of the box; other axes are fixed later."""
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


def Q(psi, i, j):
    g = [pdiff(psi, a) for a in range(3)]
    q = pmul(g[i], g[j])
    if i == j:
        s = {}
        for a in range(3):
            s = padd(s, pmul(g[a], g[a]))
        q = padd(q, pscale(s, Fraction(1, 2)), -1)
    return q


def lap(psi):
    r = {}
    for a in range(3):
        r = padd(r, pdiff(pdiff(psi, a), a))
    return r


def face_fluxes(psi, box, i):
    """Signed outward flux of row i of Q through the six faces."""
    out = {}
    for j in range(3):
        others = [a for a in range(3) if a != j]
        qij = Q(psi, i, j)
        lo, hi = box[j]
        out[(j, "+")] = const(pint(peval_axis(qij, j, hi), box, others))
        out[(j, "-")] = -const(pint(peval_axis(qij, j, lo), box, others))
    return out


def volume(psi, box, i):
    return const(pint(pmul(pdiff(psi, i), lap(psi)), box, [0, 1, 2]))


BOX = [(0, 2), (-1, 1), (0, 3)]
rng = random.Random(20261008)
monos = [(a, b, c) for a in range(4) for b in range(4) for c in range(4) if a + b + c <= 3]
n_ok = 0
for _ in range(60):
    psi = {}
    for m in rng.sample(monos, 7):
        psi[m] = Fraction(rng.randint(-5, 5))
    psi = {k: v for k, v in psi.items() if v}
    if all(sum(face_fluxes(psi, BOX, i).values()) == volume(psi, BOX, i) for i in range(3)):
        n_ok += 1
check(n_ok == 60, f"oint Q n = int d_i Psi Lap Psi on box {BOX} for {n_ok}/60 random cubic polynomials")

# Harmonic, lopsided: x^3 - 3 x y^2 + 2 y z.
H = {(3, 0, 0): Fraction(1), (1, 2, 0): Fraction(-3), (0, 1, 1): Fraction(2)}
check(lap(H) == {}, "Psi_H = x^3 - 3xy^2 + 2yz is harmonic")
GH = [sum(face_fluxes(H, BOX, i).values()) for i in range(3)]
check(GH == [0, 0, 0], f"harmonic Psi_H: G = {tuple(str(g) for g in GH)}")
fx = face_fluxes(H, BOX, 0)
abs_total = sum(abs(v) for v in fx.values())
check(fx[(0, "+")] == Fraction(907, 5) and fx[(0, "-")] == Fraction(173, 5)
      and fx[(1, "+")] == -90 and fx[(1, "-")] == -126 and abs_total == 432,
      f"harmonic Psi_H row x faces: x=2 face {fx[(0, '+')]}, total |face| {abs_total}, signed net 0")

# Non-harmonic: Psi = x^3 gives G_x = int 18 x^3 = 18 * 4 * 2 * 3 = 432.
X3 = {(3, 0, 0): Fraction(1)}
gx = sum(face_fluxes(X3, BOX, 0).values())
check(gx == 432 == volume(X3, BOX, 0), f"Psi = x^3: G_x = {gx} = int 3x^2 * 6x dV")

# ---------------------------------------------------------------- Part 2 ----
NT, NP = 96, 192
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
S = np.sqrt(1 - U ** 2)
NRM = np.stack([S * np.cos(PH), S * np.sin(PH), U], axis=-1)
W = (wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :])


def field(charges, pts, m=0.0):
    psi = np.zeros(pts.shape[:-1])
    grad = np.zeros(pts.shape)
    for q, c in charges:
        d = pts - np.asarray(c)
        r = np.linalg.norm(d, axis=-1)
        e = np.exp(-m * r)
        psi += q * e / r
        grad += (-q * e * (1 + m * r) / r ** 3)[..., None] * d
    return psi, grad


def sphere_G(charges, R, m=0.0):
    pts = R * NRM
    psi, g = field(charges, pts, m)
    gn = np.sum(g * NRM, axis=-1)
    g2 = np.sum(g * g, axis=-1)
    Qn = g * gn[..., None] - 0.5 * g2[..., None] * NRM
    G = R ** 2 * np.einsum("ij,ijk->k", W, Qn)
    P = R ** 2 * np.einsum("ij,ijk->k", W, (psi ** 2)[..., None] * NRM)
    absx = R ** 2 * np.einsum("ij,ij->", W, np.abs(Qn[..., 0]))
    aft = R ** 2 * np.einsum("ij,ij->", W, np.where(NRM[..., 0] > 0, np.abs(Qn[..., 0]), 0.0))
    return G, P, absx, aft


cluster = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]
G0, _, absx, aft = sphere_G(cluster, 2.0)
scale = absx
check(np.max(np.abs(G0)) < 1e-10 * scale and abs(absx - 15.131871) < 1e-5,
      f"isolated asymmetric cluster, R=2: |G| = {np.max(np.abs(G0)):.2e} vs oint |Q_x n| = {absx:.6f}, "
      f"x>0 share {aft / absx:.4f}")

ext = (1.5, (4.0, 1.0, -0.5))
F = np.zeros(3)
for q, c in cluster:
    d = np.asarray(c) - np.asarray(ext[1])
    F += q * ext[0] * d / np.linalg.norm(d) ** 3
G1, *_ = sphere_G(cluster + [ext], 2.0)
check(np.max(np.abs(G1 - 4 * np.pi * F)) < 1e-10 * np.max(np.abs(F))
      and abs(G1[0] + 2.9694050402) < 1e-9,
      f"external charge outside R=2: G = {np.round(G1, 10).tolist()} = 4 pi F_Coulomb, "
      f"F = {np.round(F, 10).tolist()}")

m = 0.7
rows = []
for R in (2.0, 3.0, 5.0):
    Gm, Pm, *_ = sphere_G(cluster, R, m)
    rows.append((R, Gm, 0.5 * m ** 2 * Pm))
ok = all(np.max(np.abs(Gm - bm)) < 1e-10 * max(1e-300, np.max(np.abs(bm))) for _, Gm, bm in rows)
check(ok, "Yukawa m=0.7, isolated cluster: G = (m^2/2) oint Psi^2 n at R = 2, 3, 5")
mags = [np.linalg.norm(Gm) for _, Gm, _ in rows]
REC = (1.524401, 0.3233750, 0.01713515)
check(all(abs(a - b) < 1e-6 * b for a, b in zip(mags, REC)) and mags[0] > mags[1] > mags[2],
      "Yukawa |G| at R = 2, 3, 5: " + ", ".join(f"{v:.6e}" for v in mags) + " (nonzero, shrinking)")

if FAIL:
    print(f"{len(FAIL)} check(s) failed")
    sys.exit(1)
print("all Theorem N checks passed; none of this is thrust")
