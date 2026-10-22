#!/usr/bin/env python3
"""Theorem AV witness: nonlocal (convolution / fractional) readings of a static isolated source.

Family (Assumption AV). Psi = A_0 is the static massless potential of a compactly
supported charge (point charges allowed). A nonlocal reading replaces Psi by
Phi_a = K_a * Psi, K_a translation-invariant (a finite radial Borel measure, e.g.
the heat kernel e^{t Lap} or the Yukawa smearing mu^2/(mu^2 - Lap), or the
fractional multiplier (-Lap)^{s/2}, 0<s<2), and uses the symmetric Maxwell form
    T[Phi1,Phi2]^{ij} = 1/2 (d^iPhi1 d^jPhi2 + d^iPhi2 d^jPhi1) - 1/2 delta^{ij} grad Phi1 . grad Phi2,
    d_j T^{ij} = 1/2 (d^iPhi1 Lap Phi2 + d^iPhi2 Lap Phi1) = -1/2 (sig2 d^iPhi1 + sig1 d^iPhi2),
with effective sources sig_a := -Lap Phi_a = K_a * rho (NOT compactly supported unless K is).
G(S) = oint_S T n.

AV.1  G(S) = -int_{outside S} div T  (decay at infinity), and when sig in L1 n Linf n C1,
      int_{R^3} div T = 0 (antisymmetry of grad G: Newton's third law / the k=0 mode),
      so also G(S) = int_{inside S} div T.
AV.2  G(S_k) -> 0 along vacuum surfaces running to infinity; so a value that is the same
      on every vacuum surface is 0 (holds for ANY reading with o(R^-2) decay).
AV.3  For a finite radial measure K: div T = 0 on the vacuum region for every source
      iff G = 0 on every vacuum surface for every source iff K = c delta (Maxwell).
AV.4  One radial source, sphere of radius R whose centre is b from it (exact):
      G = pi e_b int_{R-b}^{R+b} m m' (1 - u0^2) / r^2 dr,  m = -r^2 Phi',  u0 = (r^2+b^2-R^2)/(2rb).
      (-Lap)^{1/2}: G = -(64/(3 pi)) q^2 R^3 b / (R^2-b^2)^4 e_b, exact.
Exits nonzero if any identity or recorded figure moves. Nothing here is thrust.
"""
import sys
import random
import numpy as np
import sympy as sp
import mpmath as mp
from math import erf, sqrt, pi, exp, gamma, sin

FAIL = []


def check(ok, msg):
    print(("  ok   " if ok else "  FAIL ") + msg)
    if not ok:
        FAIL.append(msg)


x, y, z = sp.symbols("x y z")
XS = (x, y, z)
rs, Rs, bs, ts, mus, qs = sp.symbols("r R b t mu q", positive=True)


def lap(e):
    return sum(sp.diff(e, v, 2) for v in XS)


# ------------------------------------------------------------- Part 0: exact
print("=== Part 0: exact identities ===")
rng = random.Random(20261022)
mons = [x**a * y**b * z**c for a in range(4) for b in range(4) for c in range(4) if a + b + c <= 3]
TR, good = 30, 0
for _ in range(TR):
    p1 = sum(rng.randint(-3, 3) * m for m in rng.sample(mons, 6))
    p2 = sum(rng.randint(-3, 3) * m for m in rng.sample(mons, 6))
    g1 = [sp.diff(p1, v) for v in XS]; g2 = [sp.diff(p2, v) for v in XS]
    dot = sum(a * c for a, c in zip(g1, g2))
    ok = True
    for i in range(3):
        div = sum(sp.diff(sp.Rational(1, 2) * (g1[i] * g2[j] + g2[i] * g1[j]) - (sp.Rational(1, 2) * dot if i == j else 0), XS[j]) for j in range(3))
        ok &= sp.expand(div - sp.Rational(1, 2) * (g1[i] * lap(p2) + g2[i] * lap(p1))) == 0
    good += ok
check(good == TR, f"div T[Phi1,Phi2] = 1/2(grad Phi1 Lap Phi2 + grad Phi2 Lap Phi1) for {good}/{TR} random cubics (off-shell)")

# radial kernels acting on q/r: Phi, m = -r^2 Phi', sig = -Lap Phi = m'/r^2
lapr = lambda f: sp.diff(rs**2 * sp.diff(f, rs), rs) / rs**2
PhiG = qs * sp.erf(rs / (2 * sp.sqrt(ts))) / rs                # heat kernel e^{t Lap}
PhiY = qs * (1 - sp.exp(-mus * rs)) / rs                        # Yukawa mu^2/(mu^2 - Lap)
check(sp.simplify(sp.diff(PhiG, ts) - lapr(PhiG)) == 0, "heat-kernel reading: d_t Phi = Lap Phi (Phi = e^{t Lap} q/r)")
check(sp.limit(PhiG * rs, ts, 0) == qs, "heat-kernel reading -> q/r as t -> 0")
check(sp.simplify(mus**2 * PhiY - lapr(PhiY) - mus**2 * qs / rs) == 0, "Yukawa reading: (mu^2 - Lap) Phi = mu^2 q/r off the charge")
mG = sp.simplify(-rs**2 * sp.diff(PhiG, rs))
mY = sp.simplify(-rs**2 * sp.diff(PhiY, rs))
check(sp.simplify(mG - qs * (sp.erf(rs / (2 * sp.sqrt(ts))) - rs / sp.sqrt(sp.pi * ts) * sp.exp(-rs**2 / (4 * ts)))) == 0,
      "heat kernel: m(r) = q(erf(r/2sqrt t) - r e^{-r^2/4t}/sqrt(pi t))")
check(sp.simplify(mY - qs * (1 - (1 + mus * rs) * sp.exp(-mus * rs))) == 0, "Yukawa: m(r) = q(1 - (1+mu r) e^{-mu r})")
check(sp.simplify(sp.diff(mG, rs) - qs * rs**2 * sp.exp(-rs**2 / (4 * ts)) / (2 * sp.sqrt(sp.pi) * ts**sp.Rational(3, 2))) == 0
      and sp.simplify(sp.diff(mY, rs) - qs * mus**2 * rs * sp.exp(-mus * rs)) == 0,
      "m' > 0 for both positive kernels (sig = m'/r^2 > 0 everywhere: smeared charge leaks into the vacuum)")
for f, nm in [(PhiG, "heat"), (PhiY, "Yukawa")]:
    mm = -rs**2 * sp.diff(f, rs)
    check(sp.simplify(-lapr(f) - sp.diff(mm, rs) / rs**2) == 0, f"{nm}: sig = -Lap Phi = m'/r^2")

# fractional (-Lap)^{s/2} q/r = q C_s r^{-1-s},  C_s = (2/pi) Gamma(s) sin(pi s/2)
#   (radial inverse FT of 4 pi q k^{s-2}; Dirichlet/Abel integral int_0^inf k^{s-1} sin(kr) dk = Gamma(s) sin(pi s/2) r^{-s})
Cs = lambda s: 2 / mp.pi * mp.gamma(s) * mp.sin(mp.pi * s / 2)
riesz = lambda s: 2**s * mp.gamma((1 + s) / 2) / (mp.sqrt(mp.pi) * mp.gamma(1 - s / 2))
mp.mp.dps = 30
okC = all(abs(Cs(mp.mpf(s)) - riesz(mp.mpf(s))) < mp.mpf(10)**-25 for s in ("0.25", "0.5", "1", "1.5", "1.75"))
check(okC, "C_s = (2/pi)Gamma(s)sin(pi s/2) equals the Riesz constant 2^s Gamma((1+s)/2)/(sqrt(pi)Gamma(1-s/2)) at s=1/4,1/2,1,3/2,7/4")
dirich = 2 * mp.quadosc(lambda u: mp.sin(u**2), [0, mp.inf], zeros=lambda n: mp.sqrt(mp.pi * n))  # k = u^2
check(abs(dirich - mp.gamma(0.5) * mp.sin(mp.pi / 4)) < 1e-20, f"Dirichlet integral at s=1/2: {mp.nstr(dirich, 15)} = Gamma(1/2) sin(pi/4)")
check(sp.simplify(sp.Rational(2, 1) / sp.pi * sp.gamma(1) * sp.sin(sp.pi / 2) - 2 / sp.pi) == 0, "C_1 = 2/pi: (-Lap)^{1/2} q/r = 2q/(pi r^2)")
s_ = sp.Symbol("s", positive=True)
PhiF = qs * sp.Symbol("C") * rs**(-1 - s_)
mF = sp.simplify(-rs**2 * sp.diff(PhiF, rs)); sigF = sp.simplify(-lapr(PhiF))
check(sp.simplify(mF - (1 + s_) * sp.Symbol("C") * qs * rs**(-s_)) == 0 and sp.simplify(sigF + s_ * (1 + s_) * sp.Symbol("C") * qs * rs**(-3 - s_)) == 0,
      "fractional: m = (1+s)C_s q r^{-s}, sig = -s(1+s)C_s q r^{-3-s} (nonzero at EVERY vacuum point, sign opposite to q)")

# AV.4 exact for s = 1
u0 = (rs**2 + bs**2 - Rs**2) / (2 * rs * bs)
I5 = sp.simplify(sp.integrate(sp.expand(rs**-5 * (1 - u0**2)), (rs, Rs - bs, Rs + bs)))
check(sp.simplify(I5 - sp.Rational(4, 3) * Rs**3 * bs / (Rs**2 - bs**2)**4) == 0, "int_{R-b}^{R+b} r^-5 (1-u0^2) dr = 4R^3 b / (3 (R^2-b^2)^4) exactly")
pref = sp.pi * 1 * (1 + 1)**2 * (2 / sp.pi)**2          # pi s (1+s)^2 C_s^2 at s = 1
Gexact = sp.simplify(-pref * I5)
check(sp.simplify(Gexact + sp.Rational(64, 3) / sp.pi * Rs**3 * bs / (Rs**2 - bs**2)**4) == 0,
      "AV.4 at s=1, unit charge: G = -(64/(3 pi)) R^3 b/(R^2-b^2)^4 e_b")
EXACT = {}
for RR, bb in [(2, 0), (2, sp.Rational(1, 2)), (2, 1), (3, 1)]:
    v = sp.nsimplify(Gexact.subs({Rs: RR, bs: bb}) * sp.pi)
    EXACT[(RR, bb)] = float(Gexact.subs({Rs: RR, bs: bb}))
    print(f"     s=1, R={RR}, b={bb}: G.e_b = ({v})/pi = {EXACT[(RR, bb)]:.15f}")
check(sp.nsimplify(Gexact.subs({Rs: 2, bs: sp.Rational(1, 2)}) * sp.pi) == sp.Rational(-65536, 151875)
      and sp.nsimplify(Gexact.subs({Rs: 2, bs: 1}) * sp.pi) == sp.Rational(-512, 243)
      and sp.nsimplify(Gexact.subs({Rs: 3, bs: 1}) * sp.pi) == sp.Rational(-9, 64),
      "exact rationals: G.e_b = -65536/(151875 pi), -512/(243 pi), -9/(64 pi) at (R,b)=(2,1/2),(2,1),(3,1)")
# AV.4 exact for s = 1/2: pi s (1+s)^2 C_s^2 = 9/4 (C_{1/2}^2 = 2/pi), rational in R, b
I4 = sp.simplify(sp.integrate(sp.expand(rs**-4 * (1 - u0**2)), (rs, Rs - bs, Rs + bs)))
Ghalf = sp.simplify(-sp.Rational(9, 4) * I4)
check(sp.simplify((2 / sp.pi * sp.gamma(sp.Rational(1, 2)) * sp.sin(sp.pi / 4))**2 - 2 / sp.pi) == 0
      and sp.simplify(Ghalf + 3 * bs * (5 * Rs**2 - bs**2) / (5 * (Rs**2 - bs**2)**3)) == 0,
      "AV.4 at s=1/2, unit charge: G = -3b(5R^2-b^2)/(5(R^2-b^2)^3) e_b (rational)")
EXACT_HALF = {(2, 0.5): sp.Rational(-632, 5625), (2, 1.0): sp.Rational(-19, 45), (3, 1.0): sp.Rational(-33, 640)}
check(all(Ghalf.subs({Rs: k[0], bs: sp.nsimplify(k[1])}) == v for k, v in EXACT_HALF.items()),
      "exact: G.e_b = -632/5625, -19/45, -33/640 at (R,b)=(2,1/2),(2,1),(3,1)")
e_ = sp.Symbol("e", positive=True)
ser = sp.series(Gexact.subs(Rs, 1 / e_), e_, 0, 7).removeO()
check(sp.simplify(ser + sp.Rational(64, 3) / sp.pi * bs * e_**5) == 0, "dilation: G.e_b = -(64 b/(3 pi)) R^-5 + O(R^-7) -> 0, set by the centre offset b")

# --------------------------------------------- Part 1: single charge, 1D vs 2D
print("\n=== Part 1: one unit charge at the origin; sphere centre b e_x, radius R (1D formula vs surface quadrature) ===")
NT, NP = 128, 256
uu, wu = np.polynomial.legendre.leggauss(NT)
phis = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(uu, phis, indexing="ij")
SN = np.sqrt(1 - U**2)
NRM = np.stack([SN * np.cos(PH), SN * np.sin(PH), U], axis=-1)
Wgrid = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]
T_HEAT, MU, S_HALF = 0.25, 1.0, 0.5
derf = np.vectorize(erf)

KERNELS = {  # unit-charge m(r) (vectorised), m'(r)
    "Maxwell (K=delta)": (lambda r: np.ones_like(r), lambda r: np.zeros_like(r)),
    "heat t=1/4": (lambda r: derf(r / (2 * sqrt(T_HEAT))) - r / np.sqrt(pi * T_HEAT) * np.exp(-r**2 / (4 * T_HEAT)),
                   lambda r: r**2 * np.exp(-r**2 / (4 * T_HEAT)) / (2 * sqrt(pi) * T_HEAT**1.5)),
    "Yukawa mu=1": (lambda r: 1 - (1 + MU * r) * np.exp(-MU * r), lambda r: MU**2 * r * np.exp(-MU * r)),
    "(-Lap)^{1/2}": (lambda r: 2 * (2 / pi) / r, lambda r: -2 * (2 / pi) / r**2),
    "(-Lap)^{1/4}": (lambda r: 1.5 * float(Cs(0.5)) * r**-0.5, lambda r: -0.75 * float(Cs(0.5)) * r**-1.5),
}
MPK = {  # high-precision versions for the 1D formula
    "Maxwell (K=delta)": (lambda r: 1, lambda r: 0),
    "heat t=1/4": (lambda r: mp.erf(r / (2 * mp.sqrt(T_HEAT))) - r / mp.sqrt(mp.pi * T_HEAT) * mp.exp(-r**2 / (4 * T_HEAT)),
                   lambda r: r**2 * mp.exp(-r**2 / (4 * T_HEAT)) / (2 * mp.sqrt(mp.pi) * mp.mpf(T_HEAT)**1.5)),
    "Yukawa mu=1": (lambda r: 1 - (1 + r) * mp.exp(-r), lambda r: r * mp.exp(-r)),
    "(-Lap)^{1/2}": (lambda r: 4 / mp.pi / r, lambda r: -4 / mp.pi / r**2),
    "(-Lap)^{1/4}": (lambda r: mp.mpf(1.5) * Cs(mp.mpf(0.5)) * r**-0.5, lambda r: -mp.mpf(0.75) * Cs(mp.mpf(0.5)) * r**-1.5),
}


def G1d(name, R, b):
    if b == 0:
        return mp.mpf(0)
    m, mpr = MPK[name]
    f = lambda r: m(r) * mpr(r) * (1 - ((r**2 + b**2 - R**2) / (2 * r * b))**2) / r**2
    return mp.pi * mp.quad(f, [R - b, R, R + b])


def surface_G(fieldfun, center, R):
    pts = np.asarray(center) + R * NRM
    E = fieldfun(pts)
    En = np.sum(E * NRM, axis=-1)
    Tn = En[..., None] * E - 0.5 * np.sum(E * E, axis=-1)[..., None] * NRM
    return np.einsum("ij,ijk->k", R**2 * Wgrid, Tn)


def radial_field(mfun, charges):
    def F(pts):
        E = np.zeros(pts.shape)
        for q, c in charges:
            d = pts - np.asarray(c); r = np.linalg.norm(d, axis=-1)
            E += (-q * mfun(r) / r**3)[..., None] * d
        return E
    return F


REC1 = {}
for name, (mf, _) in KERNELS.items():
    for R, b in [(2.0, 0.0), (2.0, 0.5), (2.0, 1.0), (3.0, 1.0)]:
        g1 = float(G1d(name, mp.mpf(R), mp.mpf(b)))
        g2 = surface_G(radial_field(mf, [(1.0, (0, 0, 0))]), (b, 0, 0), R)
        REC1[(name, R, b)] = g1
        err = max(abs(g2[0] - g1), abs(g2[1]), abs(g2[2]))
        check(err < 1e-11, f"{name:18s} R={R} b={b}: 1D G.e_b = {g1:+.12e}, surface quadrature agrees ({err:.1e})")
        if name == "(-Lap)^{1/2}":
            check(abs(g1 - EXACT[(int(R), sp.Rational(b).limit_denominator(4))]) < 1e-13, "   matches the exact AV.4 closed form")
        if name == "(-Lap)^{1/4}" and b > 0:
            check(abs(g1 - float(EXACT_HALF[(int(R), b)])) < 1e-13, f"   matches the exact rational {EXACT_HALF[(int(R), b)]}")
for name in ["heat t=1/4", "Yukawa mu=1"]:
    check(all(REC1[(name, R, b)] > 0 for R, b in [(2.0, 0.5), (2.0, 1.0), (3.0, 1.0)]) and REC1[(name, 2.0, 0.0)] == 0,
          f"{name}: centred sphere 0, every off-centre sphere G.e_b > 0 (positive kernel): surface-dependent")
for name in ["(-Lap)^{1/2}", "(-Lap)^{1/4}"]:
    check(all(REC1[(name, R, b)] < 0 for R, b in [(2.0, 0.5), (2.0, 1.0), (3.0, 1.0)]), f"{name}: off-centre spheres give G.e_b < 0: surface-dependent")
check(all(abs(REC1[("Maxwell (K=delta)", R, b)]) == 0 for R, b in [(2.0, 0.5), (2.0, 1.0), (3.0, 1.0)]), "Maxwell (K=delta): m' = 0, G = 0 on every sphere")
# positive-kernel ceiling 0 < G.e_b <= pi q (q - m(R-b)) / (R-b)^2
for name in ["heat t=1/4", "Yukawa mu=1"]:
    for R, b in [(2.0, 1.0), (3.0, 1.0), (6.0, 1.0)]:
        g = float(G1d(name, mp.mpf(R), mp.mpf(b)))
        ceil = float(mp.pi * (1 - MPK[name][0](mp.mpf(R - b))) / (R - b)**2)
        check(0 < g <= ceil, f"{name} R={R} b=1: 0 < {g:.4e} <= pi(q-m(R-b))/(R-b)^2 = {ceil:.4e}")
# dilation: fractional power law vs exponential
gq = [float(G1d("(-Lap)^{1/4}", mp.mpf(R), mp.mpf(1))) for R in (10, 20, 40)]
lead = lambda R: -float(mp.pi * 0.5 * 1.5**2 * Cs(mp.mpf(0.5))**2) * (4.0 / 3.0) * R**-4.0
check(all(abs(g / lead(R) - 1) < 0.05 for g, R in zip(gq, (10, 20, 40))) and abs(gq[2] / lead(40) - 1) < abs(gq[0] / lead(10) - 1),
      f"(-Lap)^{{1/4}}, b=1: G ~ -pi s(1+s)^2 C_s^2 (4b/3) R^(-3-2s) (ratios {[round(g / lead(R), 5) for g, R in zip(gq, (10, 20, 40))]})")

# ------------------------------------------- Part 2: the Coulomb cluster
print("\n=== Part 2: cluster q=(3,-1,2) of Theorems N-Z; three vacuum spheres ===")
CLUSTER = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]
SURF = [((0.0, 0.0, 0.0), 2.0), ((0.6, 0.2, -0.3), 2.4), ((-0.9, 0.0, 0.5), 3.0)]
COUL = radial_field(lambda r: np.ones_like(r), CLUSTER)


def sig_field(mpfun):
    def S(pts):
        out = np.zeros(pts.shape[:-1])
        for q, c in CLUSTER:
            r = np.linalg.norm(pts - np.asarray(c), axis=-1)
            out += q * mpfun(r) / r**2
        return out
    return S


# coarser sphere grid for the exterior volume integrals (integrand smooth there)
uE, wuE = np.polynomial.legendre.leggauss(64)
phE = 2 * np.pi * np.arange(128) / 128
UE, PHE = np.meshgrid(uE, phE, indexing="ij")
SNE = np.sqrt(1 - UE**2)
NRME = np.stack([SNE * np.cos(PHE), SNE * np.sin(PHE), UE], axis=-1)
WgridE = wuE[:, None] * np.full(128, 2 * np.pi / 128)[None, :]


def exterior(center, R, integrand, power_law, length=16.0):
    """int_{|x-c|>R} integrand dV via radial GL panels (exponential) or r = R/tau (power law)."""
    out = np.zeros(3)
    if power_law:
        tr, wr = np.polynomial.legendre.leggauss(160); tr, wr = 0.5 * (tr + 1), 0.5 * wr
        nodes = [(R / t, R / t**2 * w) for t, w in zip(tr, wr)]
    else:
        g, w0 = np.polynomial.legendre.leggauss(32); nodes = []
        for a in np.arange(R, R + length, 1.0):
            nodes += [(a + 0.5 * (gg + 1), 0.5 * ww) for gg, ww in zip(g, w0)]
    for r, w in nodes:
        pts = np.asarray(center) + r * NRME
        out += r**2 * w * np.einsum("ij,ijk->k", WgridE, integrand(pts))
    return out


REC2 = {}
for name, pl, L in [("heat t=1/4", False, 12.0), ("Yukawa mu=1", False, 44.0), ("(-Lap)^{1/2}", True, 0.0)]:
    mf, mpf_ = KERNELS[name]
    Ef, Sf = radial_field(mf, CLUSTER), sig_field(mpf_)
    REC2[name] = []
    for c, R in SURF:
        G = surface_G(Ef, c, R)
        ext = exterior(c, R, lambda p: Sf(p)[..., None] * Ef(p), pl, L)     # AV.1: G = int_ext sig grad Phi
        REC2[name].append(G)
        check(np.max(np.abs(G - ext)) < 1e-10, f"{name}: c={c} R={R}: G = {np.array2string(G, precision=10)} = int_outside sig grad Phi ({np.max(np.abs(G - ext)):.1e})")
        Gm = surface_G(COUL, c, R)
        check(np.max(np.abs(Gm)) < 1e-11, f"   Maxwell piece on the same sphere: |G_Q| = {np.max(np.abs(Gm)):.1e}; removal leaves G_K")
    spread = max(np.max(np.abs(REC2[name][i] - REC2[name][j])) for i in range(3) for j in range(3))
    check(spread > 1e-4, f"{name}: same device, three vacuum spheres, three values (spread {spread:.4e})")

# AV.1 interior form (Newton's third law / k=0 mode) for the smooth heat-kernel reading
mf, mpf_ = KERNELS["heat t=1/4"]
Ef, Sf = radial_field(mf, CLUSTER), sig_field(mpf_)
for c, R in SURF[:2]:
    gr, wr = np.polynomial.legendre.leggauss(96)
    vol = np.zeros(3)
    for gg, ww in zip(gr, wr):
        r = 0.5 * R * (gg + 1); pts = np.asarray(c) + r * NRM
        vol -= r**2 * 0.5 * R * ww * np.einsum("ij,ijk->k", Wgrid, Sf(pts)[..., None] * Ef(pts))
    G = surface_G(Ef, c, R)
    check(np.max(np.abs(G - vol)) < 1e-9, f"heat: G(S) = int_inside div T = -int_V sig grad Phi on c={c} R={R} ({np.max(np.abs(G - vol)):.1e}); so int_R3 div T = 0")

# cross reading T[Psi, Phi]: on vacuum div = -1/2 sig grad Psi
for c, R in SURF:
    pts = np.asarray(c) + R * NRM
    E1, E2 = COUL(pts), Ef(pts)
    Tn = 0.5 * (np.sum(E2 * NRM, -1)[..., None] * E1 + np.sum(E1 * NRM, -1)[..., None] * E2) - 0.5 * np.sum(E1 * E2, -1)[..., None] * NRM
    G = np.einsum("ij,ijk->k", R**2 * Wgrid, Tn)
    ext = exterior(c, R, lambda p: 0.5 * Sf(p)[..., None] * COUL(p), False, 12.0)
    check(np.max(np.abs(G - ext)) < 1e-10 and np.max(np.abs(G)) > 1e-4,
          f"cross T[Psi, e^(Lap/4)Psi] c={c}: G = {np.array2string(G, precision=8)} = 1/2 int_out sig grad Psi")

print("\n=== Dilation (AV.2): G -> 0, so no surface-independent nonzero value ===")
for name, Rlist in [("heat t=1/4", (3.0, 4.0, 5.0)), ("Yukawa mu=1", (5.0, 10.0, 20.0)), ("(-Lap)^{1/2}", (5.0, 10.0, 20.0))]:
    Ef = radial_field(KERNELS[name][0], CLUSTER)
    Gs = [np.max(np.abs(surface_G(Ef, (0.6, 0.2, -0.3), R))) for R in Rlist]
    check(Gs[2] < Gs[1] < Gs[0], f"{name}: |G| on dilated spheres R={Rlist}: {Gs[0]:.2e}, {Gs[1]:.2e}, {Gs[2]:.2e}")

LOCKED = {
    "heat t=1/4": [(1.9216953999, 0.6370531893, -0.6085148589), (3.539183942, 0.8309375122, -1.314652441), (0.0831547402, 0.1220272808, -0.0600820331)],
    "Yukawa mu=1": [(1.4344502109, 0.5915081429, -0.5173454202), (2.5199205479, 0.8488218105, -1.0655074819), (-0.3265550523, 0.2917522815, 0.2468140837)],
    "(-Lap)^{1/2}": [(-2.2374241597, -0.7720884765, 0.7284680158), (-4.7064016253, -1.2940127131, 1.8225234138), (0.1814520415, -0.1797628144, -0.1251353099)],
}
LOCK1 = {("heat t=1/4", 2.0, 0.5): 0.1163628425531, ("heat t=1/4", 2.0, 1.0): 0.4313643009608, ("heat t=1/4", 3.0, 1.0): 0.01228902003482,
         ("Yukawa mu=1", 2.0, 0.5): 0.08832298677601, ("Yukawa mu=1", 2.0, 1.0): 0.2021386130048, ("Yukawa mu=1", 3.0, 1.0): 0.06782121850442}
if LOCK1 is not None:
    print("\n=== Recorded values (lock) ===")
    for k, v in LOCK1.items():
        check(abs(REC1[k] - v) < 1e-10, f"single charge {k}: {REC1[k]:.12e}")
    for nm, vals in LOCKED.items():
        for i, (got, exp_) in enumerate(zip(REC2[nm], vals)):
            check(np.max(np.abs(np.array(got) - np.array(exp_))) < 1e-8, f"{nm} sphere[{i}] recorded value")

if FAIL:
    print(f"\n{len(FAIL)} check(s) failed"); sys.exit(1)
print("\nall Theorem AV checks passed; nothing here is thrust")
