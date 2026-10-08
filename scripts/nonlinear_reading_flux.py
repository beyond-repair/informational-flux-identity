#!/usr/bin/env python3
"""Theorem Q witness: nonlinear readings of (grad Psi)^{ij} with any number of derivatives.

Class: T^{ij} is a finite sum of O(3)-covariant, constant-coefficient monomials,
each a product of m >= 1 factors d^k Psi with every k >= 1 (shift-invariant).
Static massless C6, isolated device, S a closed vacuum surface around it:
    G(S) = - int_{outside S} d_j T^{ij} dV                              (Q.1)
    G(S) = 0 for every such S if d_j T^{ij} vanishes on harmonic Psi    (Q.2)
    otherwise G(S) moves with S                                         (Q.3)

Readings used here:
    T1 = d_i d_k Psi d_j d_k Psi - 1/2 delta_ij |dd Psi|^2   (m=2, D=4)
         d_j T1 = d_i d_k Psi d_k Lap Psi                    -> zero on harmonic Psi
    T2 = d_i Psi d_j Lap Psi                    (m=2, D=4; named in Theorem P)
         d_j T2 = d_i d_j Psi d_j Lap Psi + d_i Psi Lap^2 Psi -> zero on harmonic Psi
    T3 = |grad Psi|^2 d_i Psi d_j Psi           (m=4, D=4)
         d_j T3 = d_i Psi d_j Psi d_j|grad|^2 + |grad|^2 (1/2 d_i|grad|^2 + d_i Psi Lap Psi)
    T4 = |grad Psi|^2 d_i d_j Psi               (m=3, D=4)
         d_j T4 = d_j|grad|^2 d_i d_j Psi + |grad|^2 d_i Lap Psi

Part 0 (exact, Fractions): the four divergence formulas hold as polynomial
identities; on a harmonic polynomial d_j T1 and d_j T2 vanish identically and
d_j T3, d_j T4 do not.
Part 1 (exact, Fractions): Lemma Q.0. The k-jets at y of 1/|. - x|, x in a ball
away from y, span every harmonic k-jet (rank 16 for k=3, 25 for k=4).
Part 2 (quadrature, NumPy, 4 pi eps0 = 1): point-charge cluster of Theorems N, O.
T1 gives G = 0 on three spheres; T3, T4 give three different values, each equal to
minus the exterior integral of the divergence, and decaying like R^-6, R^-5
under dilation of an off-centre sphere.

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
    H2 = {}
    for a in range(3):
        for b in range(3):
            H2 = padd(H2, pmul(H[a][b], H[a][b]))
    L = lap(psi)
    T = {1: {}, 2: {}, 3: {}, 4: {}}
    for i in range(3):
        for j in range(3):
            t1 = {}
            for k in range(3):
                t1 = padd(t1, pmul(H[i][k], H[j][k]))
            if i == j:
                t1 = padd(t1, pscale(H2, Fraction(-1, 2)))
            T[1][i, j] = t1
            T[2][i, j] = pmul(g[i], pdiff(L, j))
            T[3][i, j] = pmul(g2, pmul(g[i], g[j]))
            T[4][i, j] = pmul(g2, H[i][j])
    return g, H, g2, L, T


def div_direct(T, i):
    r = {}
    for j in range(3):
        r = padd(r, pdiff(T[i, j], j))
    return r


def div_formula(n, g, H, g2, L, i):
    dL = [pdiff(L, a) for a in range(3)]
    dg2 = [pdiff(g2, a) for a in range(3)]
    if n == 1:
        r = {}
        for k in range(3):
            r = padd(r, pmul(H[i][k], dL[k]))
        return r
    if n == 2:
        r = pmul(g[i], lap(L))
        for j in range(3):
            r = padd(r, pmul(H[i][j], dL[j]))
        return r
    if n == 3:
        r = {}
        for j in range(3):
            r = padd(r, pmul(pmul(g[i], g[j]), dg2[j]))
        r = padd(r, pmul(g2, padd(pscale(dg2[i], Fraction(1, 2)), pmul(g[i], L))))
        return r
    r = pmul(g2, dL[i])
    for j in range(3):
        r = padd(r, pmul(dg2[j], H[i][j]))
    return r


rng = random.Random(20261008)
TRIALS = 25
ok_id = {n: 0 for n in (1, 2, 3, 4)}
for _ in range(TRIALS):
    psi = {}
    for _ in range(9):
        k = tuple(rng.randint(0, 3) for _ in range(3))
        if sum(k) <= 5:
            psi[k] = psi.get(k, 0) + rng.randint(-4, 4)
    psi = {k: v for k, v in psi.items() if v}
    g, H, g2, L, T = readings(psi)
    for n in (1, 2, 3, 4):
        ok_id[n] += all(padd(div_direct(T[n], i), div_formula(n, g, H, g2, L, i), -1) == {}
                        for i in range(3))
for n in (1, 2, 3, 4):
    check(ok_id[n] == TRIALS, f"div T{n} formula holds as a polynomial identity: {ok_id[n]}/{TRIALS}")

# harmonic polynomial, Lap = 0 exactly
PSI_H = {(3, 0, 0): 1, (1, 2, 0): -3, (0, 1, 1): 2, (1, 1, 1): 1,
         (4, 0, 0): 1, (2, 2, 0): -6, (0, 4, 0): 1}
check(lap(PSI_H) == {}, "Psi_h = x^3 - 3xy^2 + 2yz + xyz + x^4 - 6x^2y^2 + y^4 is harmonic")
g, H, g2, L, T = readings(PSI_H)
for n in (1, 2):
    check(all(div_direct(T[n], i) == {} for i in range(3)),
          f"div T{n} vanishes identically on Psi_h (on-shell conserved)")
PT = (1, Fraction(1, 2), -1)
EXPECT = {3: (Fraction(858843, 32), Fraction(327369, 8), Fraction(-970097, 128)),
          4: (Fraction(19927, 8), Fraction(-45683, 4), Fraction(1024))}
for n in (3, 4):
    val = tuple(peval(div_direct(T[n], i), PT) for i in range(3))
    print(f"     div T{n}(Psi_h) at (1, 1/2, -1) = {tuple(str(v) for v in val)}")
    check(val == EXPECT[n], f"div T{n} is not zero on Psi_h: exact value at (1,1/2,-1) recorded")


# --------------------------------------------------------- Part 1: Lemma Q.0
def monos(k):
    return [(a, b, c) for a in range(k + 1) for b in range(k + 1 - a) for c in range(k + 1 - a - b)]


def jet_inv_r(z, k):
    """Taylor coefficients in h of |z+h|^{-1}, divided by |z|^{-1} (exact rational)."""
    a = sum(Fraction(t) ** 2 for t in z)
    p = {(1, 0, 0): 2 * Fraction(z[0]), (0, 1, 0): 2 * Fraction(z[1]), (0, 0, 1): 2 * Fraction(z[2]),
         (2, 0, 0): Fraction(1), (0, 2, 0): Fraction(1), (0, 0, 2): Fraction(1)}
    trunc = lambda q: {m: v for m, v in q.items() if sum(m) <= k}
    out, term, binom = {(0, 0, 0): Fraction(1)}, {(0, 0, 0): Fraction(1)}, Fraction(1)
    for n in range(1, k + 1):
        binom *= Fraction(-1, 2) - (n - 1)
        binom /= n
        term = trunc(pmul(term, pscale(p, 1 / a)))
        out = padd(out, pscale(term, binom))
    return trunc(out)


def rank_q(rows):
    M = [list(r) for r in rows]
    rk, col, ncol = 0, 0, len(M[0])
    while rk < len(M) and col < ncol:
        piv = next((r for r in range(rk, len(M)) if M[r][col] != 0), None)
        if piv is None:
            col += 1
            continue
        M[rk], M[piv] = M[piv], M[rk]
        for r in range(len(M)):
            if r != rk and M[r][col] != 0:
                f = M[r][col] / M[rk][col]
                M[r] = [x - f * y for x, y in zip(M[r], M[rk])]
        rk += 1
        col += 1
    return rk


Y = (3, Fraction(1, 2), -1)
for k, dimH in ((3, 16), (4, 25)):
    ms = monos(k)
    rows, harm = [], True
    for _ in range(dimH + 8):
        x = tuple(Fraction(rng.randint(-5, 5), 10) for _ in range(3))   # inside the unit ball
        jet = jet_inv_r(tuple(Y[a] - x[a] for a in range(3)), k)
        harm &= {m: v for m, v in lap(jet).items() if sum(m) <= k - 2} == {}
        rows.append([jet.get(m, Fraction(0)) for m in ms])
    check(harm, f"k={k}: every jet of 1/|y-x| is a harmonic Taylor polynomial")
    rk = rank_q(rows)
    check(rk == dimH, f"k={k}: jets of 1/|y-x| at y=(3,1/2,-1), x in |x|<1, have rank {rk} = dim of harmonic {k}-jets {dimH}")

# ----------------------------------------------------- Part 2: quadrature ----
NT, NP = 128, 256
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
S = np.sqrt(1 - U ** 2)
NRM = np.stack([S * np.cos(PH), S * np.sin(PH), U], axis=-1)
W = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]

cluster = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]


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


def tn(g, H):
    """T^{ij} n_j for T1, T3, T4 on a sphere (T2 vanishes pointwise in vacuum)."""
    g2 = np.sum(g * g, axis=-1)
    gn = np.sum(g * NRM, axis=-1)
    Hn = np.einsum("...ij,...j->...i", H, NRM)
    t1 = np.einsum("...ik,...k->...i", H, Hn) - 0.5 * np.sum(H * H, axis=(-2, -1))[..., None] * NRM
    t3 = (g2 * gn)[..., None] * g
    t4 = g2[..., None] * Hn
    return {1: t1, 3: t3, 4: t4}


def divs(g, H):
    """On-shell (Lap Psi = 0) divergences of T1, T3, T4."""
    g2 = np.sum(g * g, axis=-1)
    dg2 = 2 * np.einsum("...jk,...k->...j", H, g)
    d1 = np.zeros(g.shape)
    d3 = np.sum(g * dg2, axis=-1)[..., None] * g + 0.5 * g2[..., None] * dg2
    d4 = np.einsum("...j,...ij->...i", dg2, H)
    return {1: d1, 3: d3, 4: d4}


def flux(center, R):
    pts = np.asarray(center) + R * NRM
    g, H = fields(pts)
    t = tn(g, H)
    dA = R ** 2 * W
    return ({n: np.einsum("ij,ijk->k", dA, t[n]) for n in t},
            {n: np.einsum("ij,ij->", dA, np.abs(t[n][..., 0])) for n in t})


NR = 96
tr, wr = np.polynomial.legendre.leggauss(NR)
tr, wr = 0.5 * (tr + 1), 0.5 * wr


def exterior(center, R):
    """-int_{|x-c|>R} div T dV via r = R/t, t in (0,1]."""
    out = {1: np.zeros(3), 3: np.zeros(3), 4: np.zeros(3)}
    for t, w in zip(tr, wr):
        r = R / t
        pts = np.asarray(center) + r * NRM
        g, H = fields(pts)
        d = divs(g, H)
        jac = r ** 2 * (R / t ** 2) * w
        for n in out:
            out[n] -= jac * np.einsum("ij,ijk->k", W, d[n])
    return out


SURF = [("c=0 R=2", (0, 0, 0), 2.0),
        ("c=(0.6,0.2,-0.3) R=2.4", (0.6, 0.2, -0.3), 2.4),
        ("c=(-0.9,0,0.5) R=3", (-0.9, 0.0, 0.5), 3.0)]
REC = {
    3: [(-37.0759432675, -9.2196731911, 10.0472161178),
        (-95.5534175465, -21.3721686011, 33.8903975982),
        (1.7271019689, -1.2018086361, -1.0826483846)],
    4: [(-21.7317583295, -5.4809207986, 5.9790206965),
        (-50.8676131315, -11.4164465132, 18.1281743757),
        (1.6845711842, -1.1109562809, -1.0556329878)],
}
for s, (name, cen, R) in enumerate(SURF):
    G, A = flux(cen, R)
    ext = exterior(cen, R)
    check(np.max(np.abs(G[1])) < 1e-11 and A[1] > 1,
          f"T1 (on-shell conserved) sphere {name}: |G| = {np.max(np.abs(G[1])):.1e} vs abs row-x flux {A[1]:.6f}")
    for n in (3, 4):
        rel = np.max(np.abs(G[n] - ext[n])) / np.max(np.abs(G[n]))
        print(f"     T{n} sphere {name}: G = {np.array2string(G[n], precision=10)}  (abs row-x {A[n]:.6f})")
        check(rel < 1e-9, f"T{n} sphere {name}: G = -int_outside div T to {rel:.1e} relative")
        check(np.max(np.abs(G[n] - np.array(REC[n][s]))) < 1e-8, f"T{n} sphere {name}: recorded value")
for n in (3, 4):
    vals = [np.array(v) for v in REC[n]]
    gap = min(np.max(np.abs(vals[a] - vals[b])) for a in range(3) for b in range(a + 1, 3))
    check(gap > 0.05, f"T{n}: same device, three spheres, three values (smallest pairwise gap {gap:.3f})")

# decay under dilation S_R = R S_1 of the off-centre unit sphere S_1 (centre C1, radius 1).
# Leading term: the monopole Q_tot/|x| alone on S_1 (T3 is homogeneous of degree -8, T4 of -7).
C1 = np.array((0.3, -0.2, 0.1))
QTOT = sum(q for q, _ in cluster)
saved = cluster
cluster = [(QTOT, (0.0, 0.0, 0.0))]
LEAD, _ = flux(C1, 1.0)
cluster = saved
for n, p in ((3, 6), (4, 5)):
    sc = []
    for R in (4.0, 8.0, 16.0, 32.0, 64.0):
        G, _ = flux(R * C1, R)
        sc.append(R ** p * G[n])
    err = [np.linalg.norm(v - LEAD[n]) for v in sc]
    rat = [err[a + 1] / err[a] for a in range(4)]
    rich = 2 * sc[-1] - sc[-2]
    rel = np.linalg.norm(rich - LEAD[n]) / np.linalg.norm(LEAD[n])
    print(f"     T{n}: R^{p} G on R S_1, R=4..64: " + ", ".join(np.array2string(v, precision=5) for v in sc))
    print(f"     T{n}: monopole leading coefficient {np.array2string(LEAD[n], precision=6)}; "
          f"Richardson 2s(64)-s(32) within {rel:.1e} relative")
    check(np.linalg.norm(LEAD[n]) > 100 and all(x < 0.6 for x in rat) and rat[-1] > 0.4 and rel < 5e-3,
          f"T{n}: R^{p} G on the dilated sphere tends to the nonzero monopole coefficient "
          "(error ratios " + ", ".join(f"{x:.3f}" for x in rat) + ")")

if FAIL:
    print(f"\n{len(FAIL)} check(s) failed")
    sys.exit(1)
print("\nall checks passed; nothing here is thrust")
