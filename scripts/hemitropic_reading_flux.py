#!/usr/bin/env python3
"""Theorem Z witness: hemitropic (pseudotensor) first-derivative readings.

Family (Assumption Z):
    T^{ij} = f(x,Psi,s) delta^{ij} + g(x,Psi,s) d^iPsi d^jPsi + h(x,Psi,s) eps^{ijk} d_kPsi,
s = |grad Psi|^2, f, g, h in C^1. This is the general SO(3)-covariant (hemitropic)
second-order tensor function of one scalar and one vector; Theorem Y is the h = 0
slice. Same C6 / static massless / isolated device / vacuum surfaces as Q/U/V/Y.

Divergence of the new piece (off-shell, no Lap Psi term):
    d_j(h eps^{ijk} E_k) = ((grad_x h + h_Psi E + h_s grad s) x E)^i
                         = (grad_x h x E)^i + 2 h_s (H E x E)^i.
On-shell conservation forces Y's conditions plus h_s = 0 and grad_x h = 0, so
    T = kappa delta + c Q + h(Psi) eps.E = kappa delta + c Q + eps.grad K(Psi),  K' = h.
The new conserved piece has G_i = oint K'(Psi) (n x grad Psi)_i = 0 on EVERY closed
surface near which Psi is C^2 (Stokes), whatever is inside. Non-conserved readings
move with the surface. Exits nonzero if any identity or recorded figure moves.
Nothing here is thrust.
"""
import sys
import random
import numpy as np
import sympy as sp

FAIL = []


def check(ok, msg):
    print(("  ok   " if ok else "  FAIL ") + msg)
    if not ok:
        FAIL.append(msg)


x, y, z, P, S = sp.symbols("x y z Psi s")
XS = (x, y, z)


def grad(e):
    return sp.Matrix([sp.diff(e, v) for v in XS])


def div_hemi(psi, hexpr):
    """Direct divergence d_j (h eps^{ijk} d_k psi) with h(x,Psi,s) substituted."""
    E = grad(psi)
    s = (E.T * E)[0]
    hv = hexpr.subs({P: psi, S: s})
    out = []
    for i in range(3):
        acc = 0
        for j in range(3):
            for k in range(3):
                e = sp.LeviCivita(i, j, k)
                if e:
                    acc += e * sp.diff(hv * E[k], XS[j])
        out.append(sp.expand(acc))
    return sp.Matrix(out)


def div_hemi_formula(psi, hexpr):
    E = grad(psi)
    s = (E.T * E)[0]
    sub = {P: psi, S: s}
    hx = sp.Matrix([sp.diff(hexpr, v).subs(sub) for v in XS])
    hP = sp.diff(hexpr, P).subs(sub)
    hS = sp.diff(hexpr, S).subs(sub)
    gs = grad(s)
    return (hx + hP * E + hS * gs).cross(E)


# --------------------------------------------------- Part 0: exact identities
print("=== Part 0: exact divergence identity and conservation criterion ===")
rng = random.Random(20261015)
mons3 = [x**a * y**b * z**c for a in range(4) for b in range(4) for c in range(4) if a + b + c <= 3]
monsh = [m * P**p * S**q for m in [1, x, y, z] for p in range(3) for q in range(3) if p + q <= 2]
TR = 30
good = 0
for _ in range(TR):
    psi = sum(rng.randint(-3, 3) * m for m in rng.sample(mons3, 6))
    hh = sum(rng.randint(-3, 3) * m for m in rng.sample(monsh, 4))
    d = div_hemi(psi, hh) - div_hemi_formula(psi, hh)
    good += all(sp.expand(c) == 0 for c in d)
check(good == TR, f"Z.1 divergence formula exact for {good}/{TR} random polynomial Psi (deg<=3), h(x,Psi,s)")

# non-polynomial h(Psi) on a NON-harmonic field: identically divergence-free
psi_nh = x**2 * y + sp.Rational(1, 2) * z**3 - x * z + y
for name, hh in [("Psi^2", P**2), ("sin Psi", sp.sin(P)), ("exp(-Psi)", sp.exp(-P)), ("1/(1+Psi^2)", 1 / (1 + P**2))]:
    d = div_hemi(psi_nh, hh)
    check(all(sp.simplify(c) == 0 for c in d), f"h(Psi)={name}: eps.E piece divergence-free off-shell (non-harmonic Psi)")

# harmonic polynomial; non-conserved choices
psi_h = x**3 - 3 * x * y**2 + 2 * y * z + x * y * z + z**2 - x**2
check(sp.expand(sum(sp.diff(psi_h, v, 2) for v in XS)) == 0, "test field Psi_h is harmonic")
pt = {x: 1, y: sp.Rational(1, 2), z: -1}
E0 = grad(psi_h).subs(pt)
H0 = sp.Matrix(3, 3, lambda i, j: sp.diff(psi_h, XS[i], XS[j])).subs(pt)
d_s = div_hemi(psi_h, S).subs(pt)
check(d_s == 2 * (H0 * E0).cross(E0) and d_s != sp.zeros(3, 1),
      f"h=s: on-shell div = 2 HE x E = {list(d_s)} != 0 (exact)")
d_x = div_hemi(psi_h, x).subs(pt)
check(d_x == sp.Matrix([1, 0, 0]).cross(E0) and d_x != sp.zeros(3, 1),
      f"h=x: on-shell div = e_x x E = {list(d_x)} != 0 (exact)")
d_ps = div_hemi(psi_h, P * S).subs(pt)
check(d_ps == 2 * psi_h.subs(pt) * (H0 * E0).cross(E0) and d_ps != sp.zeros(3, 1),
      "h=Psi*s: on-shell div = 2 Psi HE x E != 0 (exact; h_Psi part drops)")

# Lemma: H -> HE onto R^3 for traceless symmetric H, E != 0; E, HE, HE x E independent
ranks, dets = 0, 0
for _ in range(20):
    Ev = sp.Matrix([rng.randint(-5, 5) for _ in range(3)])
    if Ev == sp.zeros(3, 1):
        Ev[0] = 1
    basis = []
    for (i, j) in [(0, 1), (0, 2), (1, 2)]:
        M = sp.zeros(3, 3); M[i, j] = M[j, i] = 1; basis.append(M)
    for i in (0, 1):
        M = sp.zeros(3, 3); M[i, i] = 1; M[2, 2] = -1; basis.append(M)
    A = sp.Matrix.hstack(*[B * Ev for B in basis])
    ranks += A.rank() == 3
    Hr = sum((rng.randint(-4, 4) * B for B in basis), sp.zeros(3, 3))
    v = Hr * Ev
    dets += sp.Matrix.hstack(Ev, v, v.cross(Ev)).det() != 0 or v.cross(Ev) == sp.zeros(3, 1)
check(ranks == 20, f"traceless H -> HE has rank 3 for {ranks}/20 random integer E != 0")
Ev = sp.Matrix([2, -1, 3]); w = sp.Matrix([1, 2, 0]); s0 = (Ev.T * Ev)[0]
check((Ev * w.T + w * Ev.T).trace() == 0 and (Ev * w.T + w * Ev.T) * Ev == s0 * w,
      "explicit H = E w^T + w E^T (w perp E) is traceless with HE = s w")
H2 = Ev * Ev.T - s0 / 3 * sp.eye(3)
check(H2.trace() == 0 and H2 * Ev == sp.Rational(2, 3) * s0 * Ev, "explicit H = E E^T - (s/3) I gives HE = (2s/3) E")

# ------------------------------------------------------ Part 2: quadrature ---
print("\n=== Part 2: quadrature on Coulomb cluster (static massless C6, 4 pi eps0 = 1) ===")
NT, NP = 128, 256
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
SN = np.sqrt(1 - U ** 2)
NRM = np.stack([SN * np.cos(PH), SN * np.sin(PH), U], axis=-1)
Wgrid = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]
CLUSTER = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]


def fields(pts):
    Pv = np.zeros(pts.shape[:-1]); g = np.zeros(pts.shape); H = np.zeros(pts.shape + (3,))
    for q, c in CLUSTER:
        dd = pts - np.asarray(c); r = np.linalg.norm(dd, axis=-1)
        Pv += q / r
        g += (-q / r ** 3)[..., None] * dd
        H += q * (3 * dd[..., :, None] * dd[..., None, :] / r[..., None, None] ** 5
                  - np.eye(3) / r[..., None, None] ** 3)
    return Pv, g, H


def flux(center, R, reading):
    pts = np.asarray(center) + R * NRM
    Pv, E, H = fields(pts)
    s = np.sum(E * E, axis=-1)
    Tn = reading(pts, Pv, E, s, NRM)
    dA = R ** 2 * Wgrid
    return np.einsum("ij,ijk->k", dA, Tn), np.einsum("ij,ij->", dA, np.abs(Tn[..., 0]))


def maxwell(pts, Pv, E, s, n):
    En = np.sum(E * n, axis=-1)
    return En[..., None] * E - 0.5 * s[..., None] * n


def hemi(hfun):
    def T(pts, Pv, E, s, n):  # T^{ij} n_j = h eps^{ijk} n_j E_k = h (n x E)^i
        return hfun(pts, Pv, s)[..., None] * np.cross(n, E)
    return T


R_PSI = hemi(lambda p, Pv, s: Pv)
R_EXP = hemi(lambda p, Pv, s: np.exp(-Pv))
R_S = hemi(lambda p, Pv, s: s)
R_X = hemi(lambda p, Pv, s: p[..., 0])
MIX = lambda p, Pv, E, s, n: maxwell(p, Pv, E, s, n) + R_PSI(p, Pv, E, s, n)
MIXS = lambda p, Pv, E, s, n: maxwell(p, Pv, E, s, n) + R_S(p, Pv, E, s, n)

SURF = [((0.0, 0.0, 0.0), 2.0), ((0.6, 0.2, -0.3), 2.4), ((-0.9, 0.0, 0.5), 3.0)]
for nm, rd in [("Maxwell Q", maxwell), ("Psi eps.E", R_PSI), ("exp(-Psi) eps.E", R_EXP), ("Q + Psi eps.E", MIX)]:
    for c, R in SURF:
        G, ab = flux(c, R, rd)
        check(np.max(np.abs(G)) < 1e-11, f"{nm}: |G|={np.max(np.abs(G)):.1e} on c={c} R={R} (abs row-x {ab:.4f})")

# split surface: sphere around the q=3 charge only (not a vacuum surface for the device)
cs, Rs = (-0.6, 0.1, 0.0), 0.3
Gm, _ = flux(cs, Rs, maxwell)
Gp, abp = flux(cs, Rs, R_PSI)
check(np.max(np.abs(Gp)) < 1e-9 and np.max(np.abs(Gm)) > 1.0,
      f"split sphere around q=3 only: Psi eps.E gives |G|={np.max(np.abs(Gp)):.1e} (abs {abp:.3f}); Maxwell gives {np.array2string(Gm, precision=6)} (force on enclosed charge, not zero)")


def onshell_div(pts, kind):
    Pv, E, H = fields(pts)
    HE = np.einsum("...ik,...k->...i", H, E)
    if kind == "s":
        return 2 * np.cross(HE, E)
    if kind == "x":
        return np.cross(np.broadcast_to([1.0, 0, 0], E.shape), E)


NR = 96
tr, wr = np.polynomial.legendre.leggauss(NR)
tr, wr = 0.5 * (tr + 1), 0.5 * wr


def exterior_int(center, R, kind):
    out = np.zeros(3)
    for t, w in zip(tr, wr):
        r = R / t
        pts = np.asarray(center) + r * NRM
        d = onshell_div(pts, kind)
        out -= r ** 2 * (R / t ** 2) * w * np.einsum("ij,ijk->k", Wgrid, d)
    return out


REC = {}
for nm, rd, kind in [("s eps.E", R_S, "s"), ("x eps.E", R_X, None)]:
    REC[nm] = []
    for c, R in SURF:
        G, _ = flux(c, R, rd)
        REC[nm].append(G)
        print(f"     {nm} on c={c} R={R}: G = {np.array2string(G, precision=10)}")
        if kind:
            ext = exterior_int(c, R, kind)
            err = np.max(np.abs(G - ext))
            check(err < 1e-10, f"{nm}: G = -int_outside div (abs err {err:.1e})")
    spread = max(np.max(np.abs(REC[nm][a] - REC[nm][b])) for a in range(3) for b in range(3))
    check(spread > 1e-3, f"{nm}: same device, three vacuum spheres, values differ (max spread {spread:.4f})")

for c, R in SURF:
    a, _ = flux(c, R, MIXS); b, _ = flux(c, R, maxwell); d, _ = flux(c, R, R_S)
    check(np.max(np.abs(a - b - d)) < 1e-12, f"Maxwell removal from Q + s eps.E leaves exactly the s eps.E flux on c={c}")

# x eps.E reduces exactly to the P.2 reading: Stokes on K = x Psi gives
#   oint x (n x grad Psi) = e_x x oint Psi n,  and on any sphere enclosing the charges
#   oint Psi n = (4 pi/3) p_c,  p_c = sum q (c_q - c)  (depends on the centre since Qtot = 4).
print("\n=== x eps.E: exact reduction to oint Psi n (P.2 sibling) ===")
def psi_n(pts, Pv, E, s, n):
    return Pv[..., None] * n
for c, R in SURF + [((0.0, 0.0, 0.0), 50.0), ((0.6, 0.2, -0.3), 200.0)]:
    G, _ = flux(c, R, R_X)
    Pn, _ = flux(c, R, psi_n)
    pc = sum(q * (np.asarray(cq) - np.asarray(c)) for q, cq in CLUSTER)
    lim = 4 * np.pi / 3 * np.cross([1.0, 0, 0], pc)
    e1 = np.max(np.abs(G - np.cross([1.0, 0, 0], Pn)))
    e2 = np.max(np.abs(G - lim))
    check(e1 < 1e-11 and e2 < 1e-11,
          f"x eps.E on c={c} R={R}: G = e_x x oint Psi n = (4pi/3) e_x x p_c = {np.array2string(lim, precision=6)} (err {e1:.1e}, {e2:.1e})")

print("\n=== Dilation ===")
for c in [(0.0, 0.0, 0.0), (0.6, 0.2, -0.3)]:
    Gs = [np.max(np.abs(flux(c, R, R_S)[0])) for R in (5.0, 10.0, 20.0)]
    check(Gs[2] < Gs[1] < Gs[0], f"s eps.E on c={c}: |G| decays under dilation ({Gs[0]:.2e}, {Gs[1]:.2e}, {Gs[2]:.2e})")

LOCKED = {
    "s eps.E": [
        (0.0396736069, 0.1242382122, 0.2473772397),
        (0.2264880458, 0.5105596407, 0.9389671858),
        (0.0058596383, 0.0104726951, -0.0085588517),
    ],
    "x eps.E": [
        (0.0, -3.3510321638, -3.9793506945),
        (0.0, -8.3775804096, -7.3303828584),
        (0.0, 5.0265482457, -3.9793506945),
    ],
}
if LOCKED is not None:
    print("\n=== Recorded non-conserved fluxes (lock) ===")
    for nm, vals in REC.items():
        for i, (got, exp) in enumerate(zip(vals, LOCKED[nm])):
            err = np.max(np.abs(np.array(got) - np.array(exp)))
            check(err < 1e-8, f"{nm} sphere[{i}] recorded value (err {err:.1e})")

if FAIL:
    print(f"\n{len(FAIL)} check(s) failed"); sys.exit(1)
print("\nall Theorem Z checks passed; nothing here is thrust")
