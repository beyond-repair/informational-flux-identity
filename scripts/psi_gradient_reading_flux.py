#!/usr/bin/env python3
"""Theorem V witness: first-derivative O(3)-covariant readings depending on Psi and grad Psi.

Family (Assumption V):
    T^{ij} = f(Psi, s) delta^{ij} + g(Psi, s) d^i Psi d^j Psi,   s = |grad Psi|^2,
with f, g in C^1(R x [0, infinity)). This is the general isotropic tensor function
of one scalar and one vector. Theorem U is the Psi-independent slice.

On-shell (static massless C6, Lap Psi = 0 in vacuum):
    d_j T^{ij} = (f_Psi + s g_Psi) E^i + (f_s + g/2) d_i s + g_s E^i (E · grad s)
               + g E^i Lap Psi.
Conserved for every harmonic field iff g_s = 0, f_s + g/2 = 0, and
f_Psi + s g_Psi = 0, which forces
    T = kappa * delta + c * Q
(exactly Maxwell plus constant delta — same conserved class as Theorem U).
Psi-dependence does not enlarge the conserved class. Then G = 0 on every vacuum
surface around an isolated device. Otherwise G(S) = -int_outside div T, moves
with S. Remainders already treated in U decay like R^{-2}; the exp(-Psi^2) delta
Taylor piece recovers the nonzero but surface-dependent R.4 limit of Psi^2 delta.

Part 0 (exact, Fractions): divergence identity for low-degree polynomial f, g
of (Psi, s); conservation of Maxwell+kappa on Psi_h; pseudo-Maxwell
(g=Psi, f=-(Psi/2)s) has on-shell div = (s/2) E; pure f=Psi (g=0) has div = E.
Part 1 (numeric FD): non-polynomial examples f=exp(-Psi^2), g=0;
f=0, g=1/(1+Psi^2); f=-(s/2)exp(-Psi^2), g=exp(-Psi^2) (satisfies f_s+g/2=0
and g_s=0 but not f_Psi+s g_Psi=0).
Part 2 (quadrature): Maxwell G ~ 0; Psi-dependent non-conserved readings give
different G on three vacuum spheres; exterior-integral match; dilation:
pseudo-Maxwell and soft E⊗E decay; exp(-Psi^2) delta -> -R.4 limit.

Exits nonzero if any identity or recorded figure moves. Nothing here is thrust.
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


def ppow(p, n):
    if n == 0:
        return {(0, 0, 0): Fraction(1)}
    r = {(0, 0, 0): Fraction(1)}
    for _ in range(n):
        r = pmul(r, p)
    return r


def poly_fg(coeffs, psi, s):
    """f = sum coeffs[(p,q)] * Psi^p * s^q."""
    r = {}
    for (p, q), c in coeffs.items():
        if c:
            r = padd(r, pscale(pmul(ppow(psi, p), ppow(s, q)), c))
    return r


def T_poly(psi, fcoe, gcoe, i, j):
    s = grad2(psi)
    f = poly_fg(fcoe, psi, s)
    g = poly_fg(gcoe, psi, s)
    gi, gj = pdiff(psi, i), pdiff(psi, j)
    t = pmul(g, pmul(gi, gj))
    if i == j:
        t = padd(t, f)
    return t


def div_poly_direct(psi, fcoe, gcoe, i):
    d = {}
    for j in range(3):
        d = padd(d, pdiff(T_poly(psi, fcoe, gcoe, i, j), j))
    return d


def d_partials(coeffs, var):
    """Partial of poly in (Psi, s): var='Psi' or 's'."""
    r = {}
    for (p, q), c in coeffs.items():
        if var == "Psi" and p:
            r[(p - 1, q)] = r.get((p - 1, q), 0) + c * p
        elif var == "s" and q:
            r[(p, q - 1)] = r.get((p, q - 1), 0) + c * q
    return {k: v for k, v in r.items() if v != 0}


def div_poly_formula(psi, fcoe, gcoe, i):
    """(f_Psi + s g_Psi) E^i + (f_s + g/2) d_i s + g_s E^i (E·grad s) + g E^i Lap."""
    s = grad2(psi)
    f = poly_fg(fcoe, psi, s)
    g = poly_fg(gcoe, psi, s)
    f_psi = poly_fg(d_partials(fcoe, "Psi"), psi, s)
    g_psi = poly_fg(d_partials(gcoe, "Psi"), psi, s)
    f_s = poly_fg(d_partials(fcoe, "s"), psi, s)
    g_s = poly_fg(d_partials(gcoe, "s"), psi, s)
    E_i = pdiff(psi, i)
    ds_i = pdiff(s, i)
    L = lap(psi)
    r = pmul(padd(f_psi, pmul(s, g_psi)), E_i)
    r = padd(r, pmul(padd(f_s, pscale(g, Fraction(1, 2))), ds_i))
    Edots = {}
    for j in range(3):
        Edots = padd(Edots, pmul(pdiff(psi, j), pdiff(s, j)))
    r = padd(r, pmul(g_s, pmul(E_i, Edots)))
    r = padd(r, pmul(g, pmul(E_i, L)))
    return r


def rand_fg_coeffs(rng, max_p=2, max_q=2):
    r = {}
    for p in range(max_p + 1):
        for q in range(max_q + 1):
            if rng.random() < 0.45:
                c = Fraction(rng.randint(-3, 3))
                if c:
                    r[(p, q)] = c
    if not r:
        r[(0, 0)] = Fraction(rng.randint(-2, 2) or 1)
    return r


print("=== Part 0: exact divergence identity (polynomial f(Psi,s), g(Psi,s)) ===")
rng = random.Random(20261009)
monos = [(x, y, z) for x in range(4) for y in range(4) for z in range(4) if x + y + z <= 3]
TRIALS = 40
n_ok = 0
for _ in range(TRIALS):
    fcoe = rand_fg_coeffs(rng, 2, 2)
    gcoe = rand_fg_coeffs(rng, 1, 1)
    psi = {m: Fraction(rng.randint(-4, 4)) for m in rng.sample(monos, 6)}
    psi = {k: v for k, v in psi.items() if v}
    ok = True
    for i in range(3):
        if div_poly_direct(psi, fcoe, gcoe, i) != div_poly_formula(psi, fcoe, gcoe, i):
            ok = False
            break
    n_ok += ok
check(n_ok == TRIALS,
      f"d_j T^ij = (f_Psi+s g_Psi)E^i + (f_s+g/2)d_i s + g_s E^i(E·∇s) + g E^i Lap: "
      f"{n_ok}/{TRIALS}")

PSI_H = {(3, 0, 0): 1, (1, 2, 0): -3, (0, 1, 1): 2, (1, 1, 1): 1,
         (4, 0, 0): 1, (2, 2, 0): -6, (0, 4, 0): 1}
check(lap(PSI_H) == {}, "Psi_h of Theorem Q is harmonic")

fQ = {(0, 1): Fraction(-1, 2)}
gQ = {(0, 0): Fraction(1)}
check(all(div_poly_direct(PSI_H, fQ, gQ, i) == {} for i in range(3)),
      "Maxwell Q (f=-s/2, g=1): div vanishes on Psi_h")
check(all(div_poly_direct(PSI_H, {(0, 0): Fraction(1)}, {}, i) == {} for i in range(3)),
      "kappa delta (f=1, g=0): div vanishes on Psi_h")

cons_ok = 0
for k in range(-2, 3):
    for c in range(-2, 3):
        fcoe = {(0, 0): Fraction(k), (0, 1): Fraction(-c, 2)}
        gcoe = {(0, 0): Fraction(c)} if c else {}
        fcoe = {kk: vv for kk, vv in fcoe.items() if vv}
        if all(div_poly_direct(PSI_H, fcoe, gcoe, i) == {} for i in range(3)):
            cons_ok += 1
check(cons_ok == 25, f"all 25 affine Maxwell+kappa readings conserved on Psi_h: {cons_ok}/25")

# Pseudo-Maxwell: g=Psi, f=-(Psi/2)s — on-shell div = (s/2) E
f_pseudo_p = {(1, 1): Fraction(-1, 2)}
g_pseudo_p = {(1, 0): Fraction(1)}
div_pseudo = [div_poly_direct(PSI_H, f_pseudo_p, g_pseudo_p, i) for i in range(3)]
check(any(d != {} for d in div_pseudo),
      "pseudo-Maxwell (g=Psi, f=-(Psi/2)s): div does not vanish on Psi_h")
s_h = grad2(PSI_H)
expected_pseudo = [pscale(pmul(s_h, pdiff(PSI_H, i)), Fraction(1, 2)) for i in range(3)]
check(div_pseudo == expected_pseudo,
      "pseudo-Maxwell on-shell div = (s/2) E exactly on Psi_h")

f_psi_only = {(1, 0): Fraction(1)}
div_fpsi = [div_poly_direct(PSI_H, f_psi_only, {}, i) for i in range(3)]
expected_fpsi = [pdiff(PSI_H, i) for i in range(3)]
check(div_fpsi == expected_fpsi, "pure f=Psi (g=0): on-shell div = E on Psi_h")

PT = (1, Fraction(1, 2), -1)
div_pseudo_pt = tuple(peval(div_pseudo[i], PT) for i in range(3))
print(f"     div(pseudo-Maxwell) on Psi_h at (1,1/2,-1) = {tuple(str(v) for v in div_pseudo_pt)}")
check(div_pseudo_pt != (0, 0, 0), "pseudo-Maxwell is not on-shell conserved (spot eval)")


# -------------------------------- Part 1: FD check for non-polynomial f, g --
print("\n=== Part 1: finite-difference divergence for non-polynomial f(Psi,s), g ===")


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


def T_val(Psi, E, f, g):
    s = float(np.dot(E, E))
    return f(Psi, s) * np.eye(3) + g(Psi, s) * np.outer(E, E)


def div_T_fd(fn, x, f, g, h=1e-5):
    out = np.zeros(3)
    for j in range(3):
        ep = np.zeros(3)
        ep[j] = h
        Ep, _, _ = num_grad_hess_lap(fn, x + ep, h=h * 0.5)
        Em, _, _ = num_grad_hess_lap(fn, x - ep, h=h * 0.5)
        Pp, Pm = fn(x + ep), fn(x - ep)
        Tp = T_val(Pp, Ep, f, g)
        Tm = T_val(Pm, Em, f, g)
        out += (Tp[:, j] - Tm[:, j]) / (2 * h)
    return out


def div_T_formula_num(Psi, E, H, Lap, f, g, f_psi, f_s, g_psi, g_s):
    s = float(np.dot(E, E))
    ds = 2 * H @ E
    return ((f_psi(Psi, s) + s * g_psi(Psi, s)) * E
            + (f_s(Psi, s) + g(Psi, s) / 2) * ds
            + g_s(Psi, s) * np.dot(E, ds) * E
            + g(Psi, s) * E * Lap)


def harm_fn(x):
    return float(x[0] ** 3 - 3 * x[0] * x[1] ** 2 + 2 * x[1] * x[2])


def nonharm_fn(x):
    return float(x[0] ** 2 + 0.5 * x[1] * x[2])  # Lap = 2


NP_PAIRS = [
    ("exp(-Psi^2) delta",
     lambda P, s: math.exp(-P * P), lambda P, s: 0.0,
     lambda P, s: -2 * P * math.exp(-P * P), lambda P, s: 0.0,
     lambda P, s: 0.0, lambda P, s: 0.0),
    ("soft-Psi E⊗E",
     lambda P, s: 0.0, lambda P, s: 1.0 / (1.0 + P * P),
     lambda P, s: 0.0, lambda P, s: 0.0,
     lambda P, s: -2 * P / (1.0 + P * P) ** 2, lambda P, s: 0.0),
    ("pseudo-exp-Maxwell",
     lambda P, s: -0.5 * s * math.exp(-P * P), lambda P, s: math.exp(-P * P),
     lambda P, s: P * s * math.exp(-P * P), lambda P, s: -0.5 * math.exp(-P * P),
     lambda P, s: -2 * P * math.exp(-P * P), lambda P, s: 0.0),
]

pts = [np.array([0.7, -0.3, 0.4]), np.array([1.2, 0.5, -0.8]), np.array([-0.4, 0.9, 0.2])]
for name, f, g, f_psi, f_s, g_psi, g_s in NP_PAIRS:
    max_rel = 0.0
    for fn in (harm_fn, nonharm_fn):
        for x in pts:
            E, H, Lap = num_grad_hess_lap(fn, x)
            Psi = fn(x)
            fd = div_T_fd(fn, x, f, g)
            form = div_T_formula_num(Psi, E, H, Lap, f, g, f_psi, f_s, g_psi, g_s)
            scale = max(1e-8, np.linalg.norm(form), np.linalg.norm(fd))
            max_rel = max(max_rel, np.linalg.norm(fd - form) / scale)
    check(max_rel < 5e-4,
          f"{name}: FD div matches formula (harm+nonharm; max rel {max_rel:.2e})")

x0 = pts[0]
E, H, Lap = num_grad_hess_lap(harm_fn, x0)
Psi = harm_fn(x0)
form = div_T_formula_num(Psi, E, H, Lap, *NP_PAIRS[2][1:])
nd = float(np.linalg.norm(form))
check(nd > 1e-3,
      "pseudo-exp-Maxwell (f=-(s/2)e^{-Psi^2}, g=e^{-Psi^2}) not on-shell conserved "
      "(|div|=%.4f at test point)" % nd)

s0 = float(np.dot(E, E))
fs_check = NP_PAIRS[2][4](Psi, s0) + 0.5 * NP_PAIRS[2][2](Psi, s0)
gs_check = NP_PAIRS[2][6](Psi, s0)
check(abs(fs_check) < 1e-15 and abs(gs_check) < 1e-15,
      "pseudo-exp-Maxwell satisfies f_s+g/2=0 and g_s=0 (fails only f_Psi+s g_Psi)")


# ----------------------------------------------------- Part 2: quadrature ----
print("\n=== Part 2: quadrature on Coulomb cluster (static massless C6) ===")
NT, NP = 128, 256
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
S = np.sqrt(1 - U ** 2)
NRM = np.stack([S * np.cos(PH), S * np.sin(PH), U], axis=-1)
W = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]

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
    Tn = f(Psi, s)[..., None] * NRM + g(Psi, s)[..., None] * En[..., None] * E
    dA = R ** 2 * W
    G = np.einsum("ij,ijk->k", dA, Tn)
    abs_x = np.einsum("ij,ij->", dA, np.abs(Tn[..., 0]))
    return G, abs_x, Psi, E, H, s


def onshell_div(Psi, E, H, f, g, f_psi, f_s, g_psi, g_s):
    s = np.sum(E * E, axis=-1)
    ds = 2 * np.einsum("...ik,...k->...i", H, E)
    return ((f_psi(Psi, s) + s * g_psi(Psi, s))[..., None] * E
            + (f_s(Psi, s) + 0.5 * g(Psi, s))[..., None] * ds
            + (g_s(Psi, s) * np.sum(E * ds, axis=-1))[..., None] * E)


NR = 96
tr, wr = np.polynomial.legendre.leggauss(NR)
tr, wr = 0.5 * (tr + 1), 0.5 * wr


def exterior_int(center, R, f, g, f_psi, f_s, g_psi, g_s):
    out = np.zeros(3)
    for t, w in zip(tr, wr):
        r = R / t
        Psi, E, H = fields(np.asarray(center) + r * NRM)
        d = onshell_div(Psi, E, H, f, g, f_psi, f_s, g_psi, g_s)
        jac = r ** 2 * (R / t ** 2) * w
        out -= jac * np.einsum("ij,ijk->k", W, d)
    return out


def f_max(P, s):
    return -0.5 * s


def g_max(P, s):
    return np.ones_like(s)


def fp_max(P, s):
    return np.zeros_like(s)


def fs_max(P, s):
    return -0.5 * np.ones_like(s)


def gp0(P, s):
    return np.zeros_like(s)


def gs0(P, s):
    return np.zeros_like(s)


def f_expP(P, s):
    return np.exp(-P * P)


def f_psi_expP(P, s):
    return -2 * P * np.exp(-P * P)


def fs0(P, s):
    return np.zeros_like(s)


def f0(P, s):
    return np.zeros_like(s)


def g_softP(P, s):
    return 1.0 / (1.0 + P * P)


def g_psi_softP(P, s):
    return -2 * P / (1.0 + P * P) ** 2


def f_pseudo(P, s):
    return -0.5 * P * s


def g_pseudo(P, s):
    return P


def f_psi_pseudo(P, s):
    return -0.5 * s


def fs_pseudo(P, s):
    return -0.5 * P


def g_psi_pseudo(P, s):
    return np.ones_like(s)


READINGS = [
    ("Maxwell Q", f_max, g_max, fp_max, fs_max, gp0, gs0, True),
    ("exp(-Psi^2) delta", f_expP, f0, f_psi_expP, fs0, gp0, gs0, False),
    ("E⊗E/(1+Psi^2)", f0, g_softP, fp_max, fs0, g_psi_softP, gs0, False),
    ("pseudo-Maxwell Psi", f_pseudo, g_pseudo, f_psi_pseudo, fs_pseudo, g_psi_pseudo, gs0, False),
]

SURF = [("c=0 R=2", (0.0, 0.0, 0.0), 2.0),
        ("c=(0.6,0.2,-0.3) R=2.4", (0.6, 0.2, -0.3), 2.4),
        ("c=(-0.9,0,0.5) R=3", (-0.9, 0.0, 0.5), 3.0)]

REC = {}
for rn, f, g, f_psi, f_s, g_psi, g_s, conserved in READINGS:
    vals = []
    for name, cen, R in SURF:
        G, abs_x, Psi, E, H, s = flux_and_abs(cen, R, f, g)
        ext = exterior_int(cen, R, f, g, f_psi, f_s, g_psi, g_s)
        rel = np.max(np.abs(G - ext)) / max(1e-30, np.max(np.abs(G)), np.max(np.abs(ext)))
        print(f"     {rn} sphere {name}: G = {np.array2string(G, precision=10)}  "
              f"(abs row-x {abs_x:.6f}; vs exterior rel {rel:.1e})")
        if conserved:
            check(np.max(np.abs(G)) < 1e-11 and abs_x > 1,
                  f"{rn} sphere {name}: |G| = {np.max(np.abs(G)):.1e} vs abs row-x {abs_x:.6f}")
        else:
            check(rel < 1e-9, f"{rn} sphere {name}: G = -int_outside div to {rel:.1e} relative")
            vals.append(tuple(float(x) for x in G))
    if not conserved:
        REC[rn] = vals
        gap = min(np.max(np.abs(np.array(vals[a]) - np.array(vals[b])))
                  for a in range(3) for b in range(a + 1, 3))
        check(gap > 0.01, f"{rn}: same device, three spheres, three values "
              f"(smallest pairwise gap {gap:.4f})")

print("\n=== Maxwell removal leaves the Psi-dependent remainder ===")
for s_idx, (name, cen, R) in enumerate(SURF):
    Gq, _, _, _, _, _ = flux_and_abs(cen, R, f_max, g_max)
    Ge, _, _, _, _, _ = flux_and_abs(cen, R, f_expP, f0)
    Gsum, _, _, _, _, _ = flux_and_abs(cen, R,
                                      lambda P, s: f_max(P, s) + f_expP(P, s),
                                      g_max)
    check(np.max(np.abs(Gsum - (Gq + Ge))) < 1e-12,
          f"linearity of flux on {name}")
    check(np.max(np.abs(Gsum - Ge)) < 1e-11,
          f"after Maxwell removal on {name}: residual = exp(-Psi^2) delta flux "
          f"(max |G_Q|={np.max(np.abs(Gq)):.1e})")

# Dilation: exp(-Psi^2) ~ 1 - Psi^2 -> flux tends to -R.4 limit of Psi^2 delta
# (nonzero, surface-dependent). Pseudo-Maxwell and soft E⊗E decay like R^{-2}.
print("\n=== Dilation: exp(-Psi^2) -> -R.4 limit; decaying remainders -> 0 ===")
C1 = np.array((0.3, -0.2, 0.1))
QTOT = sum(q for q, _ in CLUSTER)
a = float(np.linalg.norm(C1))
r4 = (QTOT ** 2 * 2 * math.pi
      * (1 / a - (1 + a * a) / (2 * a * a) * math.log((1 + a) / (1 - a)))
      * C1 / a)
lim_exp = -r4  # exp(-Psi^2) = 1 - Psi^2 + O(Psi^4); oint 1 n = 0
print(f"     R.4 closed form for Psi^2 delta: {np.array2string(r4, precision=9)}")
print(f"     expected lim exp(-Psi^2) delta = -R.4: {np.array2string(lim_exp, precision=9)}")

mags_pseudo = []
mags_soft = []
errs_exp = []
for R in (4.0, 8.0, 16.0, 32.0, 64.0):
    Ge, _, _, _, _, _ = flux_and_abs(R * C1, R, f_expP, f0)
    Gp, _, _, _, _, _ = flux_and_abs(R * C1, R, f_pseudo, g_pseudo)
    Gs, _, _, _, _, _ = flux_and_abs(R * C1, R, f0, g_softP)
    err = float(np.linalg.norm(Ge - lim_exp))
    errs_exp.append(err)
    mags_pseudo.append(float(np.linalg.norm(Gp)))
    mags_soft.append(float(np.linalg.norm(Gs)))
    print(f"     R={R:>4}: |G_exp+R.4|={err:.6e}, |G_pseudo|={np.linalg.norm(Gp):.6e}, "
          f"|G_soft|={np.linalg.norm(Gs):.6e}")

check(np.linalg.norm(lim_exp) > 50 and errs_exp[-1] / np.linalg.norm(lim_exp) < 5e-2,
      f"dilation: G_exp(-Psi^2) -> -R.4 limit "
      f"(rel err at R=64 = {errs_exp[-1] / np.linalg.norm(lim_exp):.2e})")
check(all(errs_exp[k + 1] < 0.6 * errs_exp[k] for k in range(4)),
      "dilation: error of G_exp to -R.4 halves roughly each doubling")
check(all(mags_pseudo[k + 1] < 0.6 * mags_pseudo[k] for k in range(4)),
      "dilation: |G_pseudo-Maxwell| decreases ~O(R^{-2})")
check(mags_pseudo[-1] < 0.1,
      f"dilation: |G_pseudo|(R=64) = {mags_pseudo[-1]:.3e} -> 0")
check(all(mags_soft[k + 1] < 0.6 * mags_soft[k] for k in range(4)),
      "dilation: |G_soft E⊗E/(1+Psi^2)| decreases ~O(R^{-2})")
check(mags_soft[-1] < 0.05,
      f"dilation: |G_soft|(R=64) = {mags_soft[-1]:.3e} -> 0")

LOCKED = {
    "exp(-Psi^2) delta": [
        (0.8677036002, 0.4738897793, -0.3661840665),
        (3.8753254205, 1.7138521206, -1.7879281576),
        (-2.7289157145, 1.6232717555, 2.2585882375),
    ],
    "E⊗E/(1+Psi^2)": [
        (-0.3438719407, -0.1605731743, 0.1357384023),
        (-0.7373169450, -0.2964568398, 0.3352190949),
        (0.2564455823, -0.1551684989, -0.1982562974),
    ],
    "pseudo-Maxwell Psi": [
        (-4.8241970412, -1.6431402777, 1.5606256234),
        (-9.6387393650, -2.7964529504, 3.8240961681),
        (0.7309160876, -0.5067393424, -0.5251162281),
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
print("\nall Theorem V checks passed; nothing here is thrust")
