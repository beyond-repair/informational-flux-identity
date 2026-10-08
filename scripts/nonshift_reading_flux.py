#!/usr/bin/env python3
"""Theorem R witness: polynomial readings of (grad Psi)^{ij} that are NOT shift-invariant.

Class: T^{ij} is a finite sum of O(3)-covariant, constant-coefficient monomials,
each a product of m >= 0 factors d^k Psi with k >= 0 (undifferentiated Psi allowed).
Static massless C6, isolated device, S a closed vacuum surface around it.
    Lemma R.0: on-shell conserved  <=>  every (m, D) graded piece is.
    R.2: on-shell conserved  =>  no e delta Psi, no f delta Psi^2 part, and G(S) = 0.
    R.3: otherwise G(S) moves with S (proof of Q.3 unchanged).
    R.4: f delta Psi^2 alone gives a flux that does NOT decay under dilation of an
         off-centre sphere; the limit is Q^2 2 pi (1/a - (1+a^2)/(2a^2) ln((1+a)/(1-a))) C/a.

Readings used here:
    U1 = Psi d_i d_j Psi - 1/2 delta_ij |grad Psi|^2   (m=2, D=2; shift adds d_i d_j Psi)
         d_j U1 = Psi d_i Lap Psi                       -> zero on harmonic Psi
    U2 = delta_ij Psi^2                                (m=2, D=0; screened-mass shape)
         d_j U2 = 2 Psi d_i Psi
    U3 = Psi^2 d_i d_j Psi                             (m=3, D=2)
         d_j U3 = Psi d_i |grad Psi|^2 + Psi^2 d_i Lap Psi

Part 0 (exact, Fractions): divergence identities, non-shift-invariance of U1,
on-shell values on a harmonic polynomial, and the grading step for e, f.
Part 1 (quadrature, NumPy, 4 pi eps0 = 1): point-charge cluster of Theorems N, O, Q.
U1 gives G = 0 on three spheres; U2, U3 give three different values, each equal to
minus the exterior integral of the on-shell divergence; under dilation of an
off-centre sphere R^3 G_U3 tends to its monopole coefficient while G_U2 itself
tends to the closed-form nonzero constant of R.4.

Exits nonzero if any identity or recorded figure moves. Nothing here is thrust.
"""
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


def lap(p):
    r = {}
    for a in range(3):
        r = padd(r, pdiff(pdiff(p, a), a))
    return r


def peval(p, pt):
    s = Fraction(0)
    for k, v in p.items():
        t = Fraction(v)
        for a in range(3):
            t *= Fraction(pt[a]) ** k[a]
        s += t
    return s


def readings(psi):
    g = [pdiff(psi, a) for a in range(3)]
    H = [[pdiff(g[a], b) for b in range(3)] for a in range(3)]
    g2 = {}
    for a in range(3):
        g2 = padd(g2, pmul(g[a], g[a]))
    P2 = pmul(psi, psi)
    T = {1: {}, 2: {}, 3: {}}
    for i in range(3):
        for j in range(3):
            u1 = pmul(psi, H[i][j])
            if i == j:
                u1 = padd(u1, pscale(g2, Fraction(-1, 2)))
            T[1][i, j] = u1
            T[2][i, j] = P2 if i == j else {}
            T[3][i, j] = pmul(P2, H[i][j])
    return g, H, g2, P2, T


def div_direct(T, i):
    r = {}
    for j in range(3):
        r = padd(r, pdiff(T[i, j], j))
    return r


def div_formula(n, psi, g, g2, P2, i):
    dL = pdiff(lap(psi), i)
    if n == 1:
        return pmul(psi, dL)
    if n == 2:
        return pscale(pmul(psi, g[i]), 2)
    return padd(pmul(psi, pdiff(g2, i)), pmul(P2, dL))


rng = random.Random(20261008)
TRIALS = 25
ok_id = {n: 0 for n in (1, 2, 3)}
shift_ok = 0
for _ in range(TRIALS):
    psi = {}
    for _ in range(9):
        k = tuple(rng.randint(0, 3) for _ in range(3))
        if sum(k) <= 5:
            psi[k] = psi.get(k, 0) + rng.randint(-4, 4)
    psi = {k: v for k, v in psi.items() if v}
    g, H, g2, P2, T = readings(psi)
    for n in (1, 2, 3):
        ok_id[n] += all(padd(div_direct(T[n], i), div_formula(n, psi, g, g2, P2, i), -1) == {}
                        for i in range(3))
    # U1[psi + 1] - U1[psi] = d_i d_j psi exactly
    *_, T1s = readings(padd(psi, {(0, 0, 0): 1}))
    shift_ok += all(padd(padd(T1s[1][i, j], T[1][i, j], -1), H[i][j], -1) == {}
                    for i in range(3) for j in range(3))
for n in (1, 2, 3):
    check(ok_id[n] == TRIALS, f"div U{n} formula holds as a polynomial identity: {ok_id[n]}/{TRIALS}")
check(shift_ok == TRIALS, f"U1[Psi+1] - U1[Psi] = d_i d_j Psi (U1 is not shift-invariant): {shift_ok}/{TRIALS}")

PSI_H = {(3, 0, 0): 1, (1, 2, 0): -3, (0, 1, 1): 2, (1, 1, 1): 1,
         (4, 0, 0): 1, (2, 2, 0): -6, (0, 4, 0): 1}
check(lap(PSI_H) == {}, "Psi_h = x^3 - 3xy^2 + 2yz + xyz + x^4 - 6x^2y^2 + y^4 is harmonic")
g, H, g2, P2, T = readings(PSI_H)
check(all(div_direct(T[1], i) == {} for i in range(3)),
      "div U1 vanishes identically on Psi_h (on-shell conserved, not shift-invariant)")
PT = (1, Fraction(1, 2), -1)
EXPECT = {2: (Fraction(-297, 32), Fraction(621, 16), Fraction(-81, 16)),
          3: (Fraction(-3051, 4), Fraction(-3591, 8), Fraction(7155, 64))}
VALS = {}
for n in (2, 3):
    VALS[n] = tuple(peval(div_direct(T[n], i), PT) for i in range(3))
    print(f"     div U{n}(Psi_h) at (1, 1/2, -1) = {tuple(str(v) for v in VALS[n])}")
    check(VALS[n] == EXPECT[n], f"div U{n} is not zero on Psi_h: exact value at (1,1/2,-1) recorded")


# Lemma R.0 grading step: on h = x, div(e delta Psi + f delta Psi^2) = (e + 2 f x, 0, 0)
X = {(1, 0, 0): 1}
zero_pairs = []
for e in range(-3, 4):
    for f in range(-3, 4):
        low = padd(pscale(X, e), pscale(pmul(X, X), f))
        d = pdiff(low, 0)
        if d == {}:
            zero_pairs.append((e, f))
check(zero_pairs == [(0, 0)], "on h = x, e delta Psi + f delta Psi^2 is conserved only for (e, f) = (0, 0) (49 pairs)")

# ----------------------------------------------------- Part 1: quadrature ----
NT, NP = 128, 256
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
S = np.sqrt(1 - U ** 2)
NRM = np.stack([S * np.cos(PH), S * np.sin(PH), U], axis=-1)
W = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]

CLUSTER = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]


def fields(pts, cl):
    P = np.zeros(pts.shape[:-1])
    g = np.zeros(pts.shape)
    H = np.zeros(pts.shape + (3,))
    for q, c in cl:
        dd = pts - np.asarray(c)
        r = np.linalg.norm(dd, axis=-1)
        P += q / r
        g += (-q / r ** 3)[..., None] * dd
        H += q * (3 * dd[..., :, None] * dd[..., None, :] / r[..., None, None] ** 5
                  - np.eye(3) / r[..., None, None] ** 3)
    return P, g, H


def tn(P, g, H):
    Hn = np.einsum("...ij,...j->...i", H, NRM)
    g2 = np.sum(g * g, axis=-1)
    return {1: P[..., None] * Hn - 0.5 * g2[..., None] * NRM,
            2: (P ** 2)[..., None] * NRM,
            3: (P ** 2)[..., None] * Hn}


def divs(P, g, H):
    """On-shell (Lap Psi = 0) divergences of U1, U2, U3."""
    return {1: np.zeros(g.shape),
            2: 2 * P[..., None] * g,
            3: 2 * P[..., None] * np.einsum("...ij,...j->...i", H, g)}


def flux(center, R, cl=CLUSTER):
    P, g, H = fields(np.asarray(center) + R * NRM, cl)
    t = tn(P, g, H)
    dA = R ** 2 * W
    return ({n: np.einsum("ij,ijk->k", dA, t[n]) for n in t},
            {n: np.einsum("ij,ij->", dA, np.abs(t[n][..., 0])) for n in t})


NR = 96
tr, wr = np.polynomial.legendre.leggauss(NR)
tr, wr = 0.5 * (tr + 1), 0.5 * wr


def exterior(center, R):
    """-int_{|x-c|>R} div U dV via r = R/t, t in (0,1], angular first."""
    out = {n: np.zeros(3) for n in (1, 2, 3)}
    for t, w in zip(tr, wr):
        r = R / t
        P, g, H = fields(np.asarray(center) + r * NRM, CLUSTER)
        d = divs(P, g, H)
        jac = r ** 2 * (R / t ** 2) * w
        for n in out:
            out[n] -= jac * np.einsum("ij,ijk->k", W, d[n])
    return out


SURF = [("c=0 R=2", (0, 0, 0), 2.0),
        ("c=(0.6,0.2,-0.3) R=2.4", (0.6, 0.2, -0.3), 2.4),
        ("c=(-0.9,0,0.5) R=3", (-0.9, 0.0, 0.5), 3.0)]
REC = {
    2: [(-35.4343595363, -15.3746788139, 13.1778307095),
        (-69.8470404027, -25.997156695, 30.5781414053),
        (17.0626744284, -11.0562983204, -13.460682377)],
    3: [(-39.1942920352, -13.866363922, 12.9663781189),
        (-77.2995955976, -23.6858511431, 31.324181094),
        (5.9282613979, -4.0210589129, -4.2759207583)],
}
for s, (name, cen, R) in enumerate(SURF):
    G, A = flux(cen, R)
    ext = exterior(cen, R)
    check(np.max(np.abs(G[1])) < 1e-11 and A[1] > 1,
          f"U1 (conserved, not shift-invariant) sphere {name}: |G| = {np.max(np.abs(G[1])):.1e} vs abs row-x flux {A[1]:.6f}")
    for n in (2, 3):
        rel = np.max(np.abs(G[n] - ext[n])) / np.max(np.abs(G[n]))
        print(f"     U{n} sphere {name}: G = {np.array2string(G[n], precision=10)}  (abs row-x {A[n]:.6f})")
        check(rel < 1e-9, f"U{n} sphere {name}: G = -int_outside div U to {rel:.1e} relative")
        if REC is not None:
            check(np.max(np.abs(G[n] - np.array(REC[n][s]))) < 1e-8, f"U{n} sphere {name}: recorded value")
if REC is not None:
    for n in (2, 3):
        vals = [np.array(v) for v in REC[n]]
        gap = min(np.max(np.abs(vals[a] - vals[b])) for a in range(3) for b in range(a + 1, 3))
        check(gap > 0.05, f"U{n}: same device, three spheres, three values (smallest pairwise gap {gap:.3f})")

# dilation S_R = R S_1 of the off-centre unit sphere S_1 (centre C1, radius 1)
C1 = np.array((0.3, -0.2, 0.1))
QTOT = sum(q for q, _ in CLUSTER)
LEAD, _ = flux(C1, 1.0, [(QTOT, (0.0, 0.0, 0.0))])
a = float(np.linalg.norm(C1))
closed = QTOT ** 2 * 2 * math.pi * (1 / a - (1 + a * a) / (2 * a * a) * math.log((1 + a) / (1 - a))) * C1 / a
print(f"     R.4 closed form for U2: {np.array2string(closed, precision=9)}; monopole quadrature "
      f"{np.array2string(LEAD[2], precision=9)}")
check(np.max(np.abs(closed - LEAD[2])) < 1e-9 and np.linalg.norm(closed) > 50,
      "R.4: monopole flux of delta Psi^2 on the off-centre unit sphere matches the closed form")
for n, p in ((2, 0), (3, 3)):
    sc = []
    for R in (4.0, 8.0, 16.0, 32.0, 64.0):
        G, _ = flux(R * C1, R)
        sc.append(R ** p * G[n])
    err = [np.linalg.norm(v - LEAD[n]) for v in sc]
    rat = [err[k + 1] / err[k] for k in range(4)]
    rich = 2 * sc[-1] - sc[-2]
    rel = np.linalg.norm(rich - LEAD[n]) / np.linalg.norm(LEAD[n])
    print(f"     U{n}: R^{p} G on R S_1, R=4..64: " + ", ".join(np.array2string(v, precision=5) for v in sc))
    print(f"     U{n}: leading coefficient {np.array2string(LEAD[n], precision=6)}; "
          f"Richardson 2s(64)-s(32) within {rel:.1e} relative")
    check(np.linalg.norm(LEAD[n]) > 50 and all(x < 0.6 for x in rat) and rat[-1] > 0.4 and rel < 5e-3,
          f"U{n}: R^{p} G on the dilated sphere tends to the nonzero leading coefficient "
          "(error ratios " + ", ".join(f"{x:.3f}" for x in rat) + ")")

if FAIL:
    print(f"\n{len(FAIL)} check(s) failed")
    sys.exit(1)
print("\nall checks passed; nothing here is thrust")
