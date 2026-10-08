#!/usr/bin/env python3
"""Theorem P witness: every linear reading of (grad Psi)^{ij}, of any order.

T^{ij} = p(Lap) d_i d_j Psi + delta_ij q(Lap) Psi,   p, q polynomials.
d_j T^{ij} = d_i r(Lap) Psi  and  oint_S T n = oint_S r(Lap) Psi n,
with r(s) = s p(s) + q(s).

Part 0 (linear algebra): O(3)-invariant tensors A^{ij}_{k1..km}, symmetric in
the k's, have dimension 1, 0, 2, 0, 2 for m = 0..4. So every O(3)-covariant,
constant-coefficient, local, linear reading with up to 6 derivatives is of the
form above (delta_ij |k|^m and k_i k_j |k|^(m-2) in Fourier).

Part 1 (exact, Fractions): the divergence identity holds as a polynomial
identity and the closed-box flux equals oint r(Lap) Psi n, for seeded random
integer p, q of degree <= 2 and random integer Psi of degree <= 7. On a
harmonic Psi, every shift-invariant reading (q(0) = 0) gives G = 0 exactly.

Part 2 (quadrature, NumPy): static massless C6, Psi = sum q_k/|x - c_k|
(4 pi eps0 = 1). The one surviving non-shift-invariant term delta_ij e Psi
gives G = e oint Psi n. On every sphere enclosing the device that equals
e (4 pi / 3)(p - Q x0) (dipole moment about the centre x0); it moves with the
centre, grows like R under dilation of an off-centre sphere, and for a neutral
device two ellipsoids around the same charges give values different from the
sphere value.

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


GENS = [rot_matrix((1, 2, 3), 0.7), rot_matrix((-2, 1, 0.5), 1.3), -np.eye(3),
        np.diag([1.0, 1.0, -1.0])]


def kron_power(g, n):
    out = np.ones((1, 1))
    for _ in range(n):
        out = np.kron(out, g)
    return out


def perm_matrix(n, perm):
    dim = 3 ** n
    P = np.zeros((dim, dim))
    for idx in itertools.product(range(3), repeat=n):
        src = int(np.ravel_multi_index(idx, (3,) * n))
        dst_idx = tuple(idx[perm[a]] for a in range(n))
        P[int(np.ravel_multi_index(dst_idx, (3,) * n)), src] = 1
    return P


EXPECTED_DIM = {0: 1, 1: 0, 2: 2, 3: 0, 4: 2}
dims = {}
for m, want in EXPECTED_DIM.items():
    n = 2 + m
    dim = 3 ** n
    rows = [kron_power(g, n) - np.eye(dim) for g in GENS]
    for a in range(2, n - 1):  # adjacent transpositions of the k slots
        perm = list(range(n))
        perm[a], perm[a + 1] = perm[a + 1], perm[a]
        rows.append(perm_matrix(n, perm) - np.eye(dim))
    sv = np.linalg.svd(np.vstack(rows), compute_uv=False)
    dims[m] = int(np.sum(sv < 1e-9)) + max(0, dim - len(sv))
check(dims == EXPECTED_DIM,
      f"O(3)-invariant A^ij_(k1..km), symmetric in k: dimensions {dims} for m = 0..4 "
      "(m even: delta_ij |k|^m and k_i k_j |k|^(m-2); m odd: none)")


# ---------------------------------------------------------------- Part 1 ----
def padd(p, q, s=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + s * v
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


def poly_of_lap(coeffs, psi):
    """sum_k coeffs[k] Lap^k psi."""
    out, cur = {}, psi
    for c in coeffs:
        out = padd(out, pscale(cur, c))
        cur = lap(cur)
    return out


def T(psi, pc, qc, i, j):
    t = pdiff(pdiff(poly_of_lap(pc, psi), i), j)
    if i == j:
        t = padd(t, poly_of_lap(qc, psi))
    return t


def r_coeffs(pc, qc):
    n = max(len(pc) + 1, len(qc))
    r = [Fraction(0)] * n
    for k, c in enumerate(pc):
        r[k + 1] += c
    for k, c in enumerate(qc):
        r[k] += c
    return r


BOX = [(0, 2), (-1, 1), (0, 3)]


def box_flux_row(psi, pc, qc, i):
    flux = Fraction(0)
    for j in range(3):
        others = [a for a in range(3) if a != j]
        lo, hi = BOX[j]
        tij = T(psi, pc, qc, i, j)
        flux += const(pint(peval_axis(tij, j, hi), BOX, others))
        flux -= const(pint(peval_axis(tij, j, lo), BOX, others))
    return flux


def box_scalar_normal(f, i):
    others = [a for a in range(3) if a != i]
    lo, hi = BOX[i]
    return (const(pint(peval_axis(f, i, hi), BOX, others))
            - const(pint(peval_axis(f, i, lo), BOX, others)))


rng = random.Random(20261008 + 16)
monos = [(x, y, z) for x in range(8) for y in range(8) for z in range(8) if x + y + z <= 7]
TRIALS = 40
n_pt = n_box = 0
for _ in range(TRIALS):
    pc = [Fraction(rng.randint(-4, 4)) for _ in range(3)]
    qc = [Fraction(rng.randint(-4, 4)) for _ in range(3)]
    psi = {m: Fraction(rng.randint(-5, 5)) for m in rng.sample(monos, 9)}
    psi = {k: v for k, v in psi.items() if v}
    rc = r_coeffs(pc, qc)
    rpsi = poly_of_lap(rc, psi)
    ok_pt = ok_box = True
    for i in range(3):
        div = {}
        for j in range(3):
            div = padd(div, pdiff(T(psi, pc, qc, i, j), j))
        ok_pt &= (div == pdiff(rpsi, i))
        ok_box &= (box_flux_row(psi, pc, qc, i) == box_scalar_normal(rpsi, i))
    n_pt += ok_pt
    n_box += ok_box
check(n_pt == TRIALS, f"d_j T^ij = d_i r(Lap) Psi as polynomials, orders up to 6 derivatives: {n_pt}/{TRIALS}")
check(n_box == TRIALS, f"oint T n = oint r(Lap) Psi n on box {BOX}: {n_box}/{TRIALS}")

# Harmonic Psi: every shift-invariant reading gives G = 0; q(0) = e gives e oint Psi n.
harm = {(3, 0, 0): Fraction(1), (1, 2, 0): Fraction(-3), (0, 1, 1): Fraction(2),
        (5, 0, 0): Fraction(1), (3, 2, 0): Fraction(-10), (1, 4, 0): Fraction(5)}
assert lap(harm) == {}
abs_row = Fraction(0)
zero_ok = True
for _ in range(20):
    pc = [Fraction(rng.randint(-4, 4)) for _ in range(3)]
    qc = [Fraction(0)] + [Fraction(rng.randint(-4, 4)) for _ in range(2)]
    G = [box_flux_row(harm, pc, qc, i) for i in range(3)]
    zero_ok &= all(g == 0 for g in G)
pc, qc = [Fraction(1), Fraction(2), Fraction(-3)], [Fraction(0), Fraction(1), Fraction(4)]
for j in range(3):
    others = [a for a in range(3) if a != j]
    for t in BOX[j]:
        abs_row += abs(const(pint(peval_axis(T(harm, pc, qc, 0, j), j, t), BOX, others)))
check(zero_ok, "harmonic Psi = x^3-3xy^2+2yz + x^5-10x^3y^2+5xy^4: G = (0,0,0) exactly for 20 random "
      "shift-invariant readings (q(0) = 0)")
check(abs_row == Fraction(1584), f"  ...while row x of the reading p=(1,2,-3), q=(0,1,4) carries face total {abs_row} in absolute value")
Ge = [box_flux_row(harm, [Fraction(0)], [Fraction(1)], i) for i in range(3)]
Gpsi = [box_scalar_normal(harm, i) for i in range(3)]
check(Ge == Gpsi and Ge == [Fraction(80), Fraction(36), Fraction(0)],
      f"undifferentiated reading delta_ij Psi on the same harmonic Psi: G = oint Psi n = {[str(g) for g in Ge]}")


# ---------------------------------------------------------------- Part 2 ----
NT, NP = 128, 256
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
S = np.sqrt(1 - U ** 2)
W = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]

POS = [(-0.6, 0.1, 0.0), (0.5, 0.35, -0.2), (0.1, -0.45, 0.3)]
CHARGED = [3.0, -1.0, 2.0]   # the cluster of Theorems N and O, Q_tot = 4
NEUTRAL = [3.0, -1.0, -2.0]  # same positions, Q_tot = 0


def psi_at(pts, qs):
    out = np.zeros(pts.shape[:-1])
    for q, c in zip(qs, POS):
        out += q / np.linalg.norm(pts - np.asarray(c), axis=-1)
    return out


def ellipsoid_flux(qs, center, axes):
    a, b, c = axes
    pts = np.asarray(center) + np.stack([a * S * np.cos(PH), b * S * np.sin(PH), c * U], axis=-1)
    ndA = np.stack([b * c * S * np.cos(PH), a * c * S * np.sin(PH), a * b * U], axis=-1) * W[..., None]
    f = psi_at(pts, qs)
    return np.einsum("ij,ijk->k", f, ndA), float(np.sum(np.abs(f)[..., None] * np.abs(ndA))), ndA, pts


def dipole(qs, x0):
    return sum(q * (np.asarray(c) - np.asarray(x0)) for q, c in zip(qs, POS))


# quadrature sanity on an ellipsoid: oint n dA = 0, oint x.n dA = 4 pi abc
_, _, ndA, pts = ellipsoid_flux(NEUTRAL, (0.2, -0.1, 0.3), (2.0, 3.0, 4.0))
check(np.max(np.abs(ndA.sum(axis=(0, 1)))) < 1e-12
      and abs(np.sum(pts * ndA) - 3 * 4 * np.pi * 24 / 3) < 1e-9,
      "ellipsoid quadrature: oint n dA = 0 and oint x.n dA = 3 vol")

SPH = [("sphere c=0 R=2", (0, 0, 0), 2.0),
       ("sphere c=(0.6,0.2,-0.3) R=2.4", (0.6, 0.2, -0.3), 2.4),
       ("sphere c=(-0.9,0,0.5) R=3", (-0.9, 0.0, 0.5), 3.0)]
rec_sph = {}
for name, cen, R in SPH:
    F, scale, _, _ = ellipsoid_flux(CHARGED, cen, (R, R, R))
    want = 4 * np.pi / 3 * dipole(CHARGED, cen)
    rec_sph[name] = F
    check(np.max(np.abs(F - want)) < 1e-11 * scale,
          f"charged device, {name}: oint Psi n = (4pi/3)(p - Q x0) = {np.round(F, 9).tolist()}")
REC_SPH = {"sphere c=0 R=2": (-8.79645943, -3.979350695, 3.351032164),
           "sphere c=(0.6,0.2,-0.3) R=2.4": (-18.849555922, -7.330382858, 8.37758041),
           "sphere c=(-0.9,0,0.5) R=3": (6.283185307, -3.979350695, -5.026548246)}
check(all(np.max(np.abs(rec_sph[n] - np.array(v))) < 1e-8 for n, v in REC_SPH.items()),
      "recorded sphere values unchanged (1e-8): same device, three surfaces, three values")

c1 = np.array([0.3, -0.2, 0.1])
dil_ok, norms = True, []
for R in (4.0, 8.0, 16.0, 32.0):
    F, scale, _, _ = ellipsoid_flux(CHARGED, R * c1, (R, R, R))
    dil_ok &= np.max(np.abs(F - 4 * np.pi / 3 * dipole(CHARGED, R * c1))) < 1e-11 * scale
    norms.append(np.linalg.norm(F))
slope = -4 * np.pi / 3 * sum(CHARGED) * c1
check(dil_ok and all(n2 > 1.7 * n1 for n1, n2 in zip(norms, norms[1:])),
      "dilation S_R = R * (unit sphere at c1=(0.3,-0.2,0.1)): oint Psi n = (4pi/3)(p - Q R c1), |.| = "
      + ", ".join(f"{n:.6f}" for n in norms) + f" at R = 4, 8, 16, 32; slope -(4pi/3) Q c1 = {np.round(slope, 6).tolist()}")

pN = dipole(NEUTRAL, (0, 0, 0))
sph_vals = []
for name, cen, R in SPH:
    F, scale, _, _ = ellipsoid_flux(NEUTRAL, cen, (R, R, R))
    sph_vals.append(F)
want = 4 * np.pi / 3 * pN
check(all(np.max(np.abs(F - want)) < 1e-11 * 10 for F in sph_vals),
      f"neutral device: every enclosing sphere gives (4pi/3) p = {np.round(want, 9).tolist()}")
E1, _, _, _ = ellipsoid_flux(NEUTRAL, (0, 0, 0), (2.0, 3.0, 4.0))
E2, _, _, _ = ellipsoid_flux(NEUTRAL, (0.1, 0.0, -0.2), (3.0, 1.6, 1.8))
print(f"     neutral device, ellipsoid (2,3,4) at 0: oint Psi n = {np.round(E1, 9).tolist()}")
print(f"     neutral device, ellipsoid (3,1.6,1.8) at (0.1,0,-0.2): oint Psi n = {np.round(E2, 9).tolist()}")
REC_E = [(-15.196767609, 3.257898425, -1.061936758), (-6.204064424, 4.570118967, -1.883253718)]
check(np.max(np.abs(E1 - np.array(REC_E[0]))) < 1e-8 and np.max(np.abs(E2 - np.array(REC_E[1]))) < 1e-8,
      "recorded ellipsoid values unchanged (1e-8)")
# Independent route for the (2,3,4) ellipsoid: oint_ell Psi n = oint_sphere(R=5) Psi n - int_shell grad Psi dV,
# with the shell integral done radially from the origin (no surface normal used).
NR = 200
om = np.stack([S * np.cos(PH), S * np.sin(PH), U], axis=-1)
rho = 1 / np.sqrt((om[..., 0] / 2.0) ** 2 + (om[..., 1] / 3.0) ** 2 + (om[..., 2] / 4.0) ** 2)
xr, wr = np.polynomial.legendre.leggauss(NR)
shell = np.zeros(3)
for xi, wi in zip(xr, wr):
    rr = rho + (5.0 - rho) * (xi + 1) / 2
    pts = rr[..., None] * om
    g = np.zeros(pts.shape)
    for q, c in zip(NEUTRAL, POS):
        dd = pts - np.asarray(c)
        g += -q * dd / np.linalg.norm(dd, axis=-1)[..., None] ** 3
    shell += np.einsum("ij,ijk->k", W * wi * (5.0 - rho) / 2 * rr ** 2, g)
check(np.max(np.abs((want - shell) - E1)) < 1e-9,
      "ellipsoid (2,3,4): surface value matches sphere value minus shell volume integral of grad Psi (1e-9)")
check(min(np.linalg.norm(E1 - want), np.linalg.norm(E2 - want), np.linalg.norm(E1 - E2)) > 1e-2 * np.linalg.norm(want),
      "neutral device: the two ellipsoids and the spheres give three different values (surface dependence)")

if FAIL:
    print(f"{len(FAIL)} check(s) failed")
    sys.exit(1)
print("all Theorem P checks passed; none of this is thrust")
