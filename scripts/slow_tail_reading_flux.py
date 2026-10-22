#!/usr/bin/env python3
"""Theorem AW witness: slow-tail smoothing readings Phi = (-Lap)^{-tau} Psi of a static isolated source.

Family (Assumption AW). Psi = A_0 is the static massless potential of a compactly supported charge
(point charges and dipoles allowed); units 4 pi eps0 = 1, a unit point charge has Psi = 1/r, and
rho~ := -Lap Psi.  The reading is Phi = K_tau * rho~ where K_tau is the radial fundamental solution of
(-Lap)^{1+tau}: K_tau = C_tau r^{2tau-1}/(4 pi) with C_tau = Gamma(1/2-tau)/(4^tau sqrt(pi) Gamma(1+tau))
(a log kernel at tau = 1/2), i.e. hat Phi = |k|^{-2tau} hat Psi on k != 0.  The literal multiplier
|k|^{-2-2tau} hat rho~ is locally integrable only for tau < (1+l)/2 (l = order of the first nonzero
multipole), so beyond that the literal (-Lap)^{-tau} Psi does not exist (REJECTED as a definition) and
K_tau is its finite-part continuation.  Effective source sigma = -Lap Phi, hat sigma = |k|^{-2tau} hat rho~.
Maxwell form T[Phi,Phi], div T = grad Phi Lap Phi = -sigma grad Phi.  G(S) = oint_S T n.

AW.1  one radial source, sphere radius R, centre offset b (no decay hypothesis needed):
      G = pi e_b int_{R-b}^{R+b} m m'(1-u0^2)/r^2 dr = D_tau q^2 J_tau(R,b) e_b,
      D_tau = 2 pi tau (2tau-1)^2 C_tau^2, J_tau = int r^{4tau-3}(1-u0^2) dr ~ (4b/3) R^{4tau-3}.
AW.2  for every tau > 0: G is the same on every vacuum surface only if div T = 0 on the exterior,
      hence only if G == 0; for tau not an integer the self reading is surface-dependent for EVERY
      nonzero device (Fourier: |k|^{2tau} is not real-analytic).
AW.3  marginal tau = 3/4: one charge, G = (16/9) q^2 (c - a) EXACTLY on every off-centre sphere (R-free);
      a charged cluster tends to (16/9) Q^2 (c - x_Q).  Neutral dipole: marginal tau = 5/4,
      G = (16/75)(|p|^2 b - 8 (p.b) p) exactly, R-free.  Centred spheres give 0.
Exits nonzero if any identity or recorded figure moves. Nothing here is thrust.
"""
import sys
import numpy as np
import sympy as sp
import mpmath as mp
from math import pi

FAIL = []


def check(ok, msg):
    print(("  ok   " if ok else "  FAIL ") + msg)
    if not ok:
        FAIL.append(msg)


mp.mp.dps = 40
ts = sp.Symbol("tau", positive=True)
rr = sp.Symbol("r", positive=True)
Rs, bs = sp.symbols("R b", positive=True)
CT = sp.gamma(sp.Rational(1, 2) - ts) / (4**ts * sp.sqrt(sp.pi) * sp.gamma(1 + ts))


def Cnum(t):
    return mp.gamma(mp.mpf(1) / 2 - t) / (4**t * mp.sqrt(mp.pi) * mp.gamma(1 + t))


def riesz(alpha):  # (-Lap)^{-alpha/2} kernel constant in R^3
    return mp.gamma((3 - alpha) / 2) / (2**alpha * mp.pi**1.5 * mp.gamma(alpha / 2))


def gfix(e):
    """Exact: shift Gamma(x), x<0 rational, up by the recurrence Gamma(x)=Gamma(x+1)/x, then simplify."""
    def up(g):
        a = g.args[0]; out = sp.Integer(1)
        while a < 0:
            out /= a; a += 1
        return out * sp.gamma(a)
    return sp.simplify(e.replace(lambda u: isinstance(u, sp.gamma) and u.args[0].is_Rational and u.args[0] < 0, up))


def lap_rad(f):
    return sp.diff(f, rr, 2) + 2 * sp.diff(f, rr) / rr


print("=== Part 0: exact kernels and existence ===")
# sigma for a unit charge must be (-Lap)^{-tau}(4 pi delta) = 4 pi c_{2tau} r^{2tau-3}
for t in [sp.Rational(1, 4), sp.Rational(3, 4), sp.Rational(1, 1), sp.Rational(5, 4), sp.Rational(7, 5)]:
    sig = sp.simplify(-lap_rad(CT.subs(ts, t) * rr**(2 * t - 1)) / rr**(2 * t - 3))
    ref = 4 * mp.pi * riesz(2 * mp.mpf(t.p) / t.q)
    check(abs(mp.mpf(sp.N(sig, 45)) - ref) < mp.mpf(10)**-30,
          f"tau={t}: -Lap(C_tau r^(2tau-1)) = 4 pi c_(2tau) r^(2tau-3), coefficient {mp.nstr(ref, 15)}")
sig_half = sp.simplify(-lap_rad(-2 / sp.pi * sp.log(rr)))
check(sp.simplify(sig_half - 2 / (sp.pi * rr**2)) == 0 and abs(4 * mp.pi * riesz(1) - 2 / mp.pi) < mp.mpf(10)**-35,
      "tau=1/2 (log kernel): Phi=-(2/pi) ln r, sigma = 2/(pi r^2) = (-Lap)^(-1/2)(4 pi delta)")
check(sp.simplify(CT.subs(ts, 1)) == sp.Rational(-1, 2) and sp.simplify(-lap_rad(-rr / 2) - 1 / rr) == 0,
      "tau=1: Phi = -r/2, sigma = 1/r = Psi exactly (sigma is the Coulomb potential itself)")
# literal existence of |k|^{-2-2tau} hat rho near k=0: radial integrand k^2 * k^{-2-2tau+l}
kk = sp.Symbol("k", positive=True)
for t, l, conv in [(sp.Rational(1, 4), 0, True), (sp.Rational(1, 2), 0, False), (sp.Rational(3, 4), 0, False),
                   (sp.Rational(3, 4), 1, True), (sp.Rational(1, 1), 1, False), (sp.Rational(5, 4), 1, False)]:
    I = sp.integrate(kk**(-2 * t + l), (kk, 0, 1), conds="none")
    ok = (I.is_finite is True) == conv and (t < sp.Rational(1 + l, 2)) == conv
    check(ok, f"literal (-Lap)^(-tau)Psi, multipole order l={l}, tau={t}: " + ("exists" if conv else "diverges at k=0 (REJECTED; finite-part kernel used)"))

print("\n=== Part 0b: exact one-charge off-centre formula (AW.1) ===")
u0 = (rr**2 + bs**2 - Rs**2) / (2 * rr * bs)
DEXACT = {sp.Rational(1, 4): 1, sp.Rational(3, 4): sp.Rational(4, 3), sp.Rational(1, 1): sp.pi / 2, sp.Rational(5, 4): sp.Rational(16, 5)}
for t, d in DEXACT.items():
    D = sp.simplify(sp.gammasimp(2 * sp.pi * t * (2 * t - 1)**2 * CT.subs(ts, t)**2))
    check(gfix(D - d) == 0, f"D_tau at tau={t} = {d}")
Dhalf = sp.limit(2 * sp.pi * ts * (2 * ts - 1)**2 * CT**2, ts, sp.Rational(1, 2))
check(sp.simplify(Dhalf - 4 / sp.pi) == 0, "D_tau is continuous at tau=1/2: limit 4/pi = log-kernel value pi*(2/pi)^2")
JEX = {}
for t in [sp.Rational(1, 4), sp.Rational(1, 2), sp.Rational(3, 4), sp.Rational(1, 1), sp.Rational(5, 4)]:
    JEX[t] = sp.simplify(sp.integrate(sp.expand(rr**(4 * t - 3) * (1 - u0**2)), (rr, Rs - bs, Rs + bs)))
check(sp.simplify(JEX[sp.Rational(1, 4)] - 4 * bs / (3 * (Rs**2 - bs**2))) == 0, "J_(1/4) = 4b/(3(R^2-b^2)): G = 4 q^2 b/(3(R^2-b^2)) e_b")
check(sp.simplify(JEX[sp.Rational(3, 4)] - 4 * bs / 3) == 0, "J_(3/4) = 4b/3 for every R>b: G = (16/9) q^2 b e_b, independent of the radius")
check(sp.simplify(JEX[sp.Rational(5, 4)] - 4 * bs * (5 * Rs**2 - bs**2) / 15) == 0, "J_(5/4) = 4b(5R^2-b^2)/15: G = (64/75) q^2 b(5R^2-b^2) e_b, grows like R^2")
eps = sp.Symbol("eps", positive=True)
for t in [sp.Rational(1, 4), sp.Rational(1, 1), sp.Rational(5, 4)]:
    lead = sp.limit(JEX[t].subs(Rs, 1 / eps) * eps**(4 * t - 3), eps, 0)
    check(sp.simplify(lead - 4 * bs / 3) == 0, f"tau={t}: J_tau ~ (4b/3) R^(4tau-3) under dilation")

print("\n=== Part 0c: exact neutral dipole at the marginal tau=5/4 ===")
X, Y, Z = sp.symbols("x y z", real=True)
uu, ph = sp.symbols("u phi", real=True)
RAD = sp.sqrt(X**2 + Y**2 + Z**2)


def dipole_G_exact(t, p):
    """G on the sphere centre (0,0,b), radius R, for Phi = p.grad K_tau (sign irrelevant), via
    G = G(centred) - int_{between} div T with G(centred) = 0 by parity."""
    K = CT.subs(ts, t) * RAD**(2 * t - 1)
    Phi = sum(pi_ * sp.diff(K, v) for pi_, v in zip(p, (X, Y, Z)))
    L = sum(sp.diff(Phi, v, 2) for v in (X, Y, Z))
    sub = {X: rr * sp.sqrt(1 - uu**2) * sp.cos(ph), Y: rr * sp.sqrt(1 - uu**2) * sp.sin(ph), Z: rr * uu}
    out = []
    for v in (X, Y, Z):
        h = sp.simplify((sp.diff(Phi, v) * L).subs(sub) * rr**2)
        h = sp.integrate(sp.expand(h), (ph, 0, 2 * sp.pi))
        h = sp.integrate(sp.expand(sp.simplify(h)), (uu, -1, u0))
        h = sp.integrate(sp.expand(sp.simplify(h)), (rr, Rs - bs, Rs + bs))
        out.append(sp.simplify(sp.gammasimp(-h)))
    return out


t54 = sp.Rational(5, 4)
for p in [(0, 0, 1), (1, 0, 0), (1, 0, 1)]:
    g = dipole_G_exact(t54, p)
    pv = sp.Matrix(p); bv = sp.Matrix([0, 0, bs])
    pred = sp.Rational(16, 75) * ((pv.dot(pv)) * bv - 8 * pv.dot(bv) * pv)
    check(all(gfix(g[i] - pred[i]) == 0 for i in range(3)),
          f"tau=5/4, p={p}, offset b e_z: G = {[gfix(gi) for gi in g]} = (16/75)(|p|^2 b - 8(p.b)p), R-free")
g34 = dipole_G_exact(sp.Rational(3, 4), (0, 0, 1))
check(gfix(g34[2] + sp.Rational(4, 15) * bs / (Rs**2 - bs**2)) == 0, "tau=3/4 dipole along offset: G = -(4/15) b/(R^2-b^2) e_z -> 0 like R^-2 = R^(4tau-5)")

# ------------------------------------------------------------- numerics
NU, NP = 160, 320
gu, wu = np.polynomial.legendre.leggauss(NU)
PH = 2 * pi * np.arange(NP) / NP
U = gu[:, None] * np.ones(NP); SN = np.sqrt(1 - U**2)
NRM = np.stack([SN * np.cos(PH), SN * np.sin(PH), U], -1)
WG = wu[:, None] * np.full(NP, 2 * pi / NP)


def Cf(t):
    return float(Cnum(mp.mpf(t)))


def grad_unit(t):
    if abs(t - 0.5) < 1e-12:
        return lambda d, r: (-(2 / pi) / r**2)[..., None] * d
    C, e = Cf(t), 2 * t - 1
    return lambda d, r: (C * e * r**(e - 2))[..., None] * d


def mono(t, charges):
    g = grad_unit(t)

    def F(p):
        E = np.zeros(p.shape)
        for q, c in charges:
            d = p - np.asarray(c); E += q * g(d, np.linalg.norm(d, axis=-1))
        return E
    return F


def dip(t, pv):
    C, e, pv = Cf(t), 2 * t - 1, np.asarray(pv, float)

    def F(x):
        r = np.linalg.norm(x, axis=-1); px = x @ pv
        return C * e * ((r**(e - 2))[..., None] * pv + ((e - 2) * r**(e - 4) * px)[..., None] * x)
    return F


def surf(F1, c, R, F2=None):
    pts = np.asarray(c, float) + R * NRM
    E1 = F1(pts); E2 = E1 if F2 is None else F2(pts)
    Tn = 0.5 * (np.sum(E2 * NRM, -1)[..., None] * E1 + np.sum(E1 * NRM, -1)[..., None] * E2) - 0.5 * np.sum(E1 * E2, -1)[..., None] * NRM
    return np.einsum("ij,ijk->k", R**2 * WG, Tn), np.sum(R**2 * WG * np.linalg.norm(Tn, axis=-1))


def Dnum(t):
    return 4 / pi if abs(t - 0.5) < 1e-12 else 2 * pi * t * (2 * t - 1)**2 * Cf(t)**2


def Jnum(t, R, b):
    if b == 0:
        return 0.0
    x, w = np.polynomial.legendre.leggauss(200)
    r = R + b * x
    return b * np.sum(w * r**(4 * t - 3) * (1 - ((r**2 + b**2 - R**2) / (2 * r * b))**2))


print("\n=== Part 1: one unit charge, surface quadrature of T vs D_tau J_tau (AW.1) ===")
REC1 = {}
for t in [0.25, 0.5, 0.6, 0.75, 1.0, 1.25]:
    F = mono(t, [(1.0, (0, 0, 0))])
    for R, b in [(2.0, 0.0), (2.0, 0.5), (3.0, 1.0), (10.0, 1.0)]:
        G, _ = surf(F, (0, 0, b), R)
        pred = Dnum(t) * Jnum(t, R, b)
        REC1[(t, R, b)] = G[2]
        check(abs(G[2] - pred) < 1e-9 * max(1, abs(pred)) and abs(G[0]) + abs(G[1]) < 1e-9 * max(1, abs(pred)),
              f"tau={t} (R,b)=({R},{b}): G.e_b = {G[2]:.12g} (formula {pred:.12g})")
check(all(abs(REC1[(0.75, R, b)] - 16 / 9 * b) < 1e-11 for R, b in [(2.0, 0.5), (3.0, 1.0), (10.0, 1.0)]),
      "tau=3/4: G = (16/9) b on spheres of radius 2, 3, 10 (radius drops out)")

print("\n=== Part 2: (H) fails yet G -> 0 for 1/2 <= tau < 3/4 ===")
F = mono(0.6, [(1.0, (0, 0, 0))])
A, B = [], []
for R in [10.0, 20.0, 40.0, 80.0]:
    G, aT = surf(F, (0, 0, 1.0), R); A.append(aT); B.append(G[2])
check(A[0] < A[1] < A[2] < A[3] and abs(A[3] / A[2] - 2**0.4) < 0.01,
      f"tau=0.6: oint|T n| = {', '.join(f'{a:.4f}' for a in A)} grows like R^0.4 (decay part of (H) fails)")
check(B[0] > B[1] > B[2] > B[3] and abs(B[3] / B[2] - 2**-0.6) < 0.01,
      f"tau=0.6: signed G = {', '.join(f'{v:.6f}' for v in B)} decays like R^-0.6 = (4/3)D q^2 b R^(4tau-3)")

print("\n=== Part 3: Coulomb cluster q=(3,-1,2) of Theorems N-AV ===")
CL = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]
SURF = [((0.0, 0.0, 0.0), 2.0), ((0.6, 0.2, -0.3), 2.4), ((-0.9, 0.0, 0.5), 3.0)]
Q = sum(q for q, _ in CL); xQ = sum(q * np.array(c) for q, c in CL) / Q
REC3 = {}
for t in [0.6, 0.75, 1.0, 1.25]:
    F = mono(t, CL)
    REC3[t] = [surf(F, c, R)[0] for c, R in SURF]
    spread = max(np.max(np.abs(REC3[t][i] - REC3[t][j])) for i in range(3) for j in range(3))
    check(spread > 1.0, f"tau={t}: three vacuum spheres, three values " + "; ".join(np.array2string(v, precision=6) for v in REC3[t]))
F = mono(0.75, CL); c0 = np.array((0.6, 0.2, -0.3)); lim = 16 / 9 * Q**2 * (c0 - xQ)
errs, cen = [], []
for R in [10.0, 20.0, 40.0, 80.0]:
    errs.append(np.max(np.abs(surf(F, c0, R)[0] - lim))); cen.append(np.max(np.abs(surf(F, xQ, R)[0])))
check(errs[3] < errs[2] < errs[1] < errs[0] < 0.05 and errs[0] * 10 < 4.0,
      f"tau=3/4, centre {c0.tolist()}: G -> (16/9)Q^2(c-x_Q) = {np.array2string(lim, precision=6)}; error {', '.join(f'{e:.2e}' for e in errs)} at R=10..80")
check(cen[3] < cen[2] < cen[1] < cen[0] < 1e-2, f"tau=3/4, spheres centred at x_Q: |G| = {', '.join(f'{e:.2e}' for e in cen)} -> 0")
F1 = mono(1.0, CL); g1 = [np.max(np.abs(surf(F1, c0, R)[0])) for R in (20.0, 40.0, 80.0)]
check(abs(g1[2] / g1[1] - 2) < 0.01 and abs(g1[1] / g1[0] - 2) < 0.01, f"tau=1: |G| = {', '.join(f'{v:.2f}' for v in g1)} grows like R")
# cross reading T[Psi, Phi]
def COUL(p):  # grad Psi for the cluster
    return sum(q * (-(p - np.asarray(c)) / np.linalg.norm(p - np.asarray(c), axis=-1)[..., None]**3) for q, c in CL)
XR = [surf(COUL, c, R, F)[0] for c, R in SURF]
check(max(np.max(np.abs(XR[i] - XR[j])) for i in range(3) for j in range(3)) > 0.1,
      "cross T[Psi, (-Lap)^(-3/4)Psi]: " + "; ".join(np.array2string(v, precision=6) for v in XR))
GM = max(np.max(np.abs(surf(COUL, c, R)[0])) for c, R in SURF)
check(GM < 1e-11, f"Maxwell piece T[Psi,Psi] on the same spheres: {GM:.1e}; removing it leaves the values above")

print("\n=== Part 4: neutral dipole ===")
for t in [0.75, 1.0, 1.25, 1.4]:
    F = dip(t, (0, 0, 1.0))
    check(np.max(np.abs(surf(F, (0, 0, 0), 3.0)[0])) < 1e-12, f"tau={t}: centred sphere G = 0 (parity)")
F = dip(1.25, (0.3, -0.5, 0.8)); pv = np.array((0.3, -0.5, 0.8))
for c, R in [((0.4, 0.7, -0.2), 2.5), ((0.8, 1.4, -0.4), 7.0), ((0.4, 0.7, -0.2), 30.0)]:
    G = surf(F, c, R)[0]; c = np.array(c); pred = 16 / 75 * (pv @ pv * c - 8 * (pv @ c) * pv)
    check(np.max(np.abs(G - pred)) < 1e-11, f"tau=5/4, p={pv.tolist()}, centre {c.tolist()}, R={R}: G = {np.array2string(G, precision=8)}")
for t, ex in [(1.0, -1.0), (1.4, 0.6)]:
    F = dip(t, (0, 0, 1.0)); v = [surf(F, (0, 0, 1.0), R)[0][2] for R in (20.0, 40.0)]
    check(abs(v[1] / v[0] - 2**ex) < 0.02, f"tau={t}: dipole G along p scales like R^(4tau-5): ratio {v[1] / v[0]:.4f} vs {2**ex:.4f}")

LOCK3 = {0.75: [(14.65468537, 6.4214493, -5.47896618), (31.39891094, 11.73411026, -13.76696834), (-10.61325382, 7.03900622, 8.41002305)],
         1.0: [(34.30191405, 15.3130513, -12.96761438), (84.49466522, 32.40620401, -37.38564363), (-37.52238943, 24.00462997, 29.9091709)]}
print("\n=== Recorded values (lock) ===")
for t, vals in LOCK3.items():
    for i, v in enumerate(vals):
        check(np.max(np.abs(REC3[t][i] - np.array(v))) < 1e-7, f"cluster tau={t} sphere[{i}] recorded value")
for k, v in {(0.5, 2.0, 0.5): 0.4354705905236808, (1.0, 3.0, 1.0): 6.1412561794947615, (1.25, 10.0, 1.0): 425.81333333337}.items():
    check(abs(REC1[k] - v) < 1e-9 * max(1, v), f"single charge {k}: {REC1[k]:.12g}")

if FAIL:
    print(f"\n{len(FAIL)} check(s) failed"); sys.exit(1)
print("\nall Theorem AW checks passed; nothing here is thrust")
