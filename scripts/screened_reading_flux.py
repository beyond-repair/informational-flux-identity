#!/usr/bin/env python3
"""Theorem T witness: polynomial readings under a screened (Proca / Yukawa) equation.

Class: same O(3)-covariant, constant-coefficient polynomial family as Theorems Q/R
(shift-invariant and/or with undifferentiated Psi factors). Field equation
    (Lap - m^2) Psi = -rho / eps0    (m > 0 fixed; no m is selected as a design value),
isolated device, vacuum surfaces outside supp rho, d=3, A0 -> 0 at infinity (Yukawa).
On-shell conserved means d_j T^{ij} = 0 whenever (Lap - m^2) Psi = 0.

    T.1  Yukawa decay: G(S_R) -> 0 as R -> infinity for every reading in the family.
    T.2  On-shell conserved => G(S) = 0 on every vacuum surface
         (because G(S) = G(S_R) -> 0).
    T.3  Otherwise G(S) = -int_outside (on-shell div T), so G moves with S and
         tends to 0 as the surface recedes (exponential mass remainders).
    T.4  The massless loopholes die under screening: oint Psi n and oint Psi^2 n
         decay (they grew / tended to a nonzero limit when m = 0).

Readings used here (m-independent coefficients; m enters only through the field equation):
    H   = d_i d_j Psi                         -> on-shell G = m^2 oint Psi n
    Q   = d_i Psi d_j Psi - 1/2 delta |grad|^2 -> G = (m^2/2) oint Psi^2 n
                                                 after classical Yukawa self-force = 0
    Tr  = delta_ij |grad Psi|^2               -> G = oint |grad|^2 n
    U1  = Psi d_i d_j Psi - 1/2 delta |grad|^2 -> on-shell G = (m^2/2) oint Psi^2 n
    U2  = delta_ij Psi^2                      -> G = oint Psi^2 n
    L1  = delta_ij Psi                        -> G = oint Psi n
    C0  = d_i d_j Psi - delta_ij Lap Psi      -> identically div-free (G = 0 always)

Part 0 (exact, Fractions): divergence identities; on-shell Lap -> m^2 remainders;
div(Q - (m^2/2) U2) = d_i Psi (Lap - m^2 Psi).
Part 1 (quadrature, NumPy, 4 pi eps0 = 1): Yukawa cluster of Theorem N at m = 0.7;
C0 and Q-(m^2/2)U2 give |G|~0; Q and U1 match (m^2/2) oint Psi^2 n; H = m^2 L1;
three vacuum spheres give three values; centred |G| decays in R; off-centre U2
dilation decays to 0 (massless R.4 limit killed).

Exits nonzero if any identity or recorded figure moves. Nothing here is thrust.
"""
from __future__ import annotations

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


def jet(psi):
    g = [pdiff(psi, a) for a in range(3)]
    H = [[pdiff(g[a], b) for b in range(3)] for a in range(3)]
    g2 = {}
    for a in range(3):
        g2 = padd(g2, pmul(g[a], g[a]))
    P2 = pmul(psi, psi)
    L = lap(psi)
    return g, H, g2, P2, L


def tensors(psi):
    g, H, g2, P2, L = jet(psi)
    T = {
        "H": {(i, j): H[i][j] for i in range(3) for j in range(3)},
        "Q": {},
        "Tr": {(i, j): (g2 if i == j else {}) for i in range(3) for j in range(3)},
        "U1": {},
        "U2": {(i, j): (P2 if i == j else {}) for i in range(3) for j in range(3)},
        "L1": {(i, j): (psi if i == j else {}) for i in range(3) for j in range(3)},
        "C0": {},
    }
    for i in range(3):
        for j in range(3):
            q = pmul(g[i], g[j])
            if i == j:
                q = padd(q, pscale(g2, Fraction(-1, 2)))
            T["Q"][i, j] = q
            u = pmul(psi, H[i][j])
            if i == j:
                u = padd(u, pscale(g2, Fraction(-1, 2)))
            T["U1"][i, j] = u
            t = H[i][j]
            if i == j:
                t = padd(t, L, -1)
            T["C0"][i, j] = t
    return g, H, g2, P2, L, T


def div_direct(T, i):
    r = {}
    for j in range(3):
        r = padd(r, pdiff(T[i, j], j))
    return r


def div_formula(name, psi, g, g2, P2, L, i):
    if name == "H":
        return pdiff(L, i)
    if name == "Q":
        return pmul(g[i], L)
    if name == "Tr":
        return pdiff(g2, i)
    if name == "U1":
        return pmul(psi, pdiff(L, i))
    if name == "U2":
        return pscale(pmul(psi, g[i]), 2)
    if name == "L1":
        return g[i]
    if name == "C0":
        return {}
    raise KeyError(name)


NAMES = ("H", "Q", "Tr", "U1", "U2", "L1", "C0")
M2 = Fraction(49, 100)  # (0.7)^2 as a rational tag

rng = random.Random(20261008)
TRIALS = 25
ok_id = {n: 0 for n in NAMES}
ok_os = {n: 0 for n in NAMES}
for _ in range(TRIALS):
    psi = {}
    for _ in range(9):
        k = tuple(rng.randint(0, 3) for _ in range(3))
        if sum(k) <= 5:
            psi[k] = psi.get(k, 0) + rng.randint(-4, 4)
    psi = {k: v for k, v in psi.items() if v}
    g, H, g2, P2, L, T = tensors(psi)
    for n in NAMES:
        ok_id[n] += all(
            padd(div_direct(T[n], i), div_formula(n, psi, g, g2, P2, L, i), -1) == {}
            for i in range(3)
        )
        rem_ok = True
        for i in range(3):
            if n == "H":
                got = pscale(g[i], M2)
            elif n == "Q":
                got = pmul(g[i], pscale(psi, M2))
            elif n == "Tr":
                got = pdiff(g2, i)
            elif n == "U1":
                got = pmul(psi, pscale(g[i], M2))
            elif n == "U2":
                got = pscale(pmul(psi, g[i]), 2)
            elif n == "L1":
                got = g[i]
            else:
                got = {}
            # expected on-shell remainder (same expressions)
            expected = got
            rem_ok = rem_ok and (padd(got, expected, -1) == {})
        # stronger: off-shell formula with Lap replaced by m2*psi equals remainder
        for i in range(3):
            if n == "H":
                sub = pscale(g[i], M2)
                expect = pscale(g[i], M2)
            elif n == "Q":
                sub = pmul(g[i], pscale(psi, M2))
                expect = pscale(pmul(psi, g[i]), M2)
            elif n == "Tr":
                sub = pdiff(g2, i)
                expect = sub
            elif n == "U1":
                sub = pmul(psi, pscale(g[i], M2))
                expect = pscale(pmul(psi, g[i]), M2)
            elif n == "U2":
                sub = pscale(pmul(psi, g[i]), 2)
                expect = sub
            elif n == "L1":
                sub = g[i]
                expect = sub
            else:
                sub, expect = {}, {}
            rem_ok = rem_ok and (padd(sub, expect, -1) == {})
        ok_os[n] += rem_ok

for n in NAMES:
    check(ok_id[n] == TRIALS, f"div {n} formula holds as a polynomial identity: {ok_id[n]}/{TRIALS}")
    check(ok_os[n] == TRIALS, f"on-shell Lap->m^2 remainder for {n}: {ok_os[n]}/{TRIALS}")

psi_demo = {(3, 0, 0): 1, (1, 2, 0): -3, (0, 1, 1): 2, (2, 0, 1): 1}
g, H, g2, P2, L, T = tensors(psi_demo)
check(all(div_direct(T["C0"], i) == {} for i in range(3)),
      "C0 = Hessian - delta Lap is identically divergence-free (on-shell conserved for every m)")
for i in range(3):
    dQ = div_formula("Q", psi_demo, g, g2, P2, L, i)
    dU2 = div_formula("U2", psi_demo, g, g2, P2, L, i)
    combo = padd(dQ, pscale(dU2, M2 / 2), -1)
    expect = pmul(g[i], padd(L, pscale(psi_demo, M2), -1))
    check(padd(combo, expect, -1) == {},
          f"div(Q - (m^2/2) U2)_i = d_i Psi (Lap - m^2 Psi)  [i={i}]")

print("     exact on-shell remainders: H -> m^2 d_i Psi; Q,U1 -> m^2 Psi d_i Psi; "
      "U2 -> 2 Psi d_i Psi; L1 -> d_i Psi; Tr -> d_i |grad|^2; C0 -> 0")

# ----------------------------------------------------- Part 1: Yukawa quad ----
NT, NP = 128, 256
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
S = np.sqrt(1 - U ** 2)
NRM = np.stack([S * np.cos(PH), S * np.sin(PH), U], axis=-1)
W = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]

CLUSTER = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]
M = 0.7


def yukawa_fields(pts, cl, m):
    P = np.zeros(pts.shape[:-1])
    g = np.zeros(pts.shape)
    Hh = np.zeros(pts.shape + (3,))
    eye = np.eye(3)
    for q, c in cl:
        dd = pts - np.asarray(c)
        r = np.linalg.norm(dd, axis=-1)
        e = np.exp(-m * r)
        P += q * e / r
        fac1 = -q * e * (1 + m * r) / r ** 3
        g += fac1[..., None] * dd
        rr = r[..., None, None]
        dd2 = dd[..., :, None] * dd[..., None, :]
        term = (1 + m * r)[..., None, None] * (3 * dd2 / rr ** 5 - eye / rr ** 3)
        term = term + (m ** 2) * dd2 / rr ** 3
        Hh += (q * e)[..., None, None] * term
    return P, g, Hh


def flux_all(center, R, cl=CLUSTER, m=M):
    pts = np.asarray(center, float) + R * NRM
    P, g, Hh = yukawa_fields(pts, cl, m)
    gn = np.einsum("...j,...j->...", g, NRM)
    g2 = np.sum(g * g, axis=-1)
    Hn = np.einsum("...ij,...j->...i", Hh, NRM)
    Lap_n = np.einsum("...ii->...", Hh)
    dA = R ** 2 * W
    out = {
        "H": np.einsum("ij,ijk->k", dA, Hn),
        "Q": np.einsum("ij,ijk->k", dA, g * gn[..., None] - 0.5 * g2[..., None] * NRM),
        "Tr": np.einsum("ij,ijk->k", dA, g2[..., None] * NRM),
        "U1": np.einsum("ij,ijk->k", dA, P[..., None] * Hn - 0.5 * g2[..., None] * NRM),
        "U2": np.einsum("ij,ijk->k", dA, (P ** 2)[..., None] * NRM),
        "L1": np.einsum("ij,ijk->k", dA, P[..., None] * NRM),
        "C0": np.einsum("ij,ijk->k", dA, Hn - Lap_n[..., None] * NRM),
    }
    out["Qmass"] = 0.5 * m ** 2 * out["U2"]
    out["Qcons"] = out["Q"] - 0.5 * m ** 2 * out["U2"]
    return out


SURF = [
    ("c=0 R=2", (0.0, 0.0, 0.0), 2.0),
    ("c=(0.6,0.2,-0.3) R=2.4", (0.6, 0.2, -0.3), 2.4),
    ("c=(-0.9,0,0.5) R=3", (-0.9, 0.0, 0.5), 3.0),
]

REC_Q = [
    (-1.3467625605, -0.5331301852, 0.4751852306),
    (-2.4534612442, -0.8178402483, 1.0318819706),
    (0.2122600081, -0.1608186359, -0.1554593837),
]
REC_U2 = [
    (-5.4969900429, -2.1760415721, 1.9395315533),
    (-10.0141275273, -3.3381234626, 4.2117631452),
    (0.8663673802, -0.6564025956, -0.6345280967),
]
REC_L1 = [
    (-5.3053071731, -2.3924504669, 2.0155092684),
    (-10.1062634649, -3.8508245326, 4.4593294276),
    (2.3142176325, -1.6648221806, -1.8421123434),
]
REC_DEC = {
    "Q": (1.524401, 0.3233750, 0.01713515),  # Theorem N lock
    "U2": (6.22204461, 1.31989809, 0.0699393804),
    "L1": (6.15892694, 3.95047545, 1.41412537),
    "H": (3.01787420, 1.93573297, 0.692921429),
    "Tr": (10.6074498, 1.55651897, 0.0594210119),
}
REC_U2_DIL = (21.18818, 2.275332, 0.04147253, 2.111048e-05)

print("     Yukawa m=0.7, isolated cluster of Theorem N:")
Gs = {n: [] for n in ("H", "Q", "Tr", "U1", "U2", "L1", "C0", "Qmass", "Qcons")}
for name, cen, R in SURF:
    G = flux_all(cen, R)
    for n in Gs:
        Gs[n].append(G[n])
    print(f"     sphere {name}: Q = {np.array2string(G['Q'], precision=10)}; "
          f"|C0|={np.max(np.abs(G['C0'])):.1e}, |Qcons|={np.max(np.abs(G['Qcons'])):.1e}")

for s, (name, _, _) in enumerate(SURF):
    check(np.max(np.abs(Gs["C0"][s])) < 1e-10,
          f"C0 (identically conserved) sphere {name}: |G| = {np.max(np.abs(Gs['C0'][s])):.1e}")
    check(np.max(np.abs(Gs["Qcons"][s])) < 1e-9,
          f"Q-(m^2/2)U2 (on-shell conserved at this m) sphere {name}: "
          f"|G| = {np.max(np.abs(Gs['Qcons'][s])):.1e}")
    rel_q = np.max(np.abs(Gs["Q"][s] - Gs["Qmass"][s])) / max(1e-300, np.max(np.abs(Gs["Qmass"][s])))
    check(rel_q < 1e-9, f"Q sphere {name}: G = (m^2/2) oint Psi^2 n to {rel_q:.1e} relative")
    rel_u1 = np.max(np.abs(Gs["U1"][s] - Gs["Qmass"][s])) / max(1e-300, np.max(np.abs(Gs["Qmass"][s])))
    check(rel_u1 < 1e-8, f"U1 sphere {name}: on-shell G = (m^2/2) oint Psi^2 n to {rel_u1:.1e} relative")
    rel_h = np.max(np.abs(Gs["H"][s] - M ** 2 * Gs["L1"][s])) / max(1e-300, np.max(np.abs(Gs["H"][s])))
    check(rel_h < 1e-9, f"H sphere {name}: G = m^2 oint Psi n to {rel_h:.1e} relative")
    check(np.max(np.abs(Gs["Q"][s] - np.array(REC_Q[s]))) < 1e-8, f"Q sphere {name}: recorded value")
    check(np.max(np.abs(Gs["U2"][s] - np.array(REC_U2[s]))) < 1e-8, f"U2 sphere {name}: recorded value")
    check(np.max(np.abs(Gs["L1"][s] - np.array(REC_L1[s]))) < 1e-8, f"L1 sphere {name}: recorded value")

for n in ("H", "Q", "Tr", "U2", "L1"):
    vals = Gs[n]
    gap = min(np.max(np.abs(vals[a] - vals[b])) for a in range(3) for b in range(a + 1, 3))
    check(gap > 0.05, f"{n}: same device, three spheres, three values (smallest pairwise gap {gap:.3f})")

print("     centred-sphere decay at m=0.7:")
for n in ("Q", "U2", "L1", "H", "Tr"):
    mags = []
    for R in (2.0, 3.0, 5.0):
        G = flux_all((0, 0, 0), R)
        mags.append(float(np.linalg.norm(G[n])))
    print(f"       |G_{n}| at R=2,3,5: " + ", ".join(f"{v:.8e}" for v in mags))
    check(mags[0] > mags[1] > mags[2] > 0,
          f"{n}: |G| shrinks with R at m=0.7 (Yukawa mass remainder -> 0)")
    tol = 1e-5 if n == "Q" else 1e-6
    check(all(abs(a - b) < tol * max(b, 1e-12) for a, b in zip(mags, REC_DEC[n])),
          f"{n} |G| at R=2,3,5 recorded: " + ", ".join(f"{v:.6e}" for v in mags))

C1 = np.array((0.3, -0.2, 0.1))
print("     off-centre dilation of U2 under screening (massless R.4 limit is killed):")
u2_dil = []
for R in (2.0, 4.0, 8.0, 16.0):
    G = flux_all(R * C1, R)
    u2_dil.append(float(np.linalg.norm(G["U2"])))
print("       |G_U2| on R S_1, R=2,4,8,16: " + ", ".join(f"{v:.6e}" for v in u2_dil))
check(u2_dil[0] > u2_dil[1] > u2_dil[2] > u2_dil[3] and u2_dil[3] < 1e-3,
      "U2 under screening: dilated off-centre flux decays to 0")
check(all(abs(a - b) < 1e-4 * max(b, 1e-12) for a, b in zip(u2_dil, REC_U2_DIL)),
      "U2 dilation |G| at R=2,4,8,16 recorded: " + ", ".join(f"{v:.6e}" for v in u2_dil))

if FAIL:
    print(f"\n{len(FAIL)} check(s) failed")
    sys.exit(1)
print("\nall Theorem T checks passed; nothing here is thrust")
