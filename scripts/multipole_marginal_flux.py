#!/usr/bin/env python3
"""Theorem AX witness: marginal slow-tail readings of a static 2^l-pole (proves Conjecture AW.4).

Setting (Assumption AW). Psi = A_0, units 4 pi eps0 = 1.  A point source of symbol P (a real homogeneous
polynomial of degree l) has Psi = P(d)(1/r) and slow-tail reading Phi = P(d) K_tau, K_tau = C_tau r^(2tau-1)
(finite-part kernel of Theorem AW).  A cluster sum_i q_i delta_(a_i) has Phi = sum_n P_n(d) K_tau with
P_n(xi) = sum_i q_i (-a_i.xi)^n / n!.  Maxwell form T[Phi,Phi]; div T = f := grad Phi Lap Phi.
S_(c,R) = sphere of centre c and radius R > |c| (the source sits at 0).

AX.1  every tau (pure point source): G(S_(c,R)) = int_(S^2) f(w) (rho(w)^k - R^k)/k dw, k = 4tau-2-2l,
      rho(w) = w.c + sqrt(R^2-|c|^2+(w.c)^2); hence G = R^(k-1) M_tau c + O(|c|^3 R^(k-3)), M_tau = oint_(S^2) f (x) n.
AX.2  marginal tau_l = (3+2l)/4 (k = 1): G = M c EXACTLY for every R > |c|.  tau_l - (1+l)/2 = 1/4, so the
      literal (-Lap)^(-tau)Psi never exists at the marginal value; tau_l is never an integer or half-integer.
AX.3  P = sum_j |xi|^(2j) H_(l-2j) (H_m harmonic):  M = sum_j M_(l-2j)[H_(l-2j)],
      M_m[H] = 8/(pi (2m+3)^2) ( (N/2) I - 4/(4m^2-1) B ),  N = oint H^2,  B = oint grad H (x) grad H,
      tr M_m[H] = -4N/(pi (2m-1)(2m+3)).  Quadrupole P = xi.A xi:
      M = (16/81)(tr A)^2 I + (32/735) tr(A0^2) I - (512/2205) A0^2,  A0 = A - (tr A/3) I.
AX.4  M != 0 (l>=1) if P has no |xi|^l part (all traces < 0), and for every l = 2 symbol; M = 0 for the l = 4 symbol
      P* = H4 + (4 sqrt 11/385)|xi|^4, H4 = sum xi_i^4 - (3/5)|xi|^4: G = 0 on every enclosing sphere, not on cubes.
AX.5  P_l + P_(l+1) at tau_l: G = M c + g0 + (1/2) ln(1-|c|^2/R^2) v + M' c/(R^2-|c|^2) exactly.
Kept Failure AX.1: P* + xi_x^5 gives G = g0 = (64 sqrt11/8085 - 5888/266805) e_x on every centred sphere and
      G -> g0 on spheres of every centre; dilated cubes grow like lambda.  Nothing here is thrust.
Exits nonzero if any identity or recorded figure moves.
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


x, y, z, rho = sp.symbols("x y z rho", real=True)
X = (x, y, z)
R2 = x**2 + y**2 + z**2
lsym = sp.Symbol("l", nonnegative=True, integer=True)


def Ct(t):
    return sp.gamma(sp.Rational(1, 2) - t) / (4**t * sp.sqrt(sp.pi) * sp.gamma(1 + t))


def gfix(e):
    """Shift Gamma(x), x<0 rational, up by Gamma(x)=Gamma(x+1)/x, then simplify (exact)."""
    def up(g):
        a = g.args[0]; out = sp.Integer(1)
        while a < 0:
            out /= a; a += 1
        return out * sp.gamma(a)
    e = sp.sympify(e)
    return sp.simplify(e.replace(lambda u: isinstance(u, sp.gamma) and u.args[0].is_Rational and u.args[0] < 0, up))


def D(F, i):
    """d/dx_i of an expression in x, y, z and rho = |x|."""
    return sp.diff(F, X[i]) + sp.diff(F, rho) * X[i] / rho


def Dpoly(F, P):
    """P(d) F for a polynomial P in x, y, z."""
    out = 0
    for mon, cf in sp.Poly(P, *X).terms():
        G = F
        for i, e in enumerate(mon):
            for _ in range(e):
                G = D(G, i)
        out += cf * G
    return sp.expand(out)


def sph(expr):
    """Exact oint_(S^2) of a polynomial: oint x^a y^b z^c = 4 pi (a-1)!!(b-1)!!(c-1)!!/(a+b+c+1)!!."""
    p = sp.Poly(sp.expand(expr), *X); s = 0
    for (a, b, c), cf in p.terms():
        if a % 2 or b % 2 or c % 2:
            continue
        s += cf * 4 * sp.factorial2(a - 1) * sp.factorial2(b - 1) * sp.factorial2(c - 1) / sp.factorial2(a + b + c + 1)
    return sp.expand(s) * sp.pi


def unit(e):
    return sp.expand(e.subs(rho, 1))


def fields(P, t):
    Phi = Ct(t) * Dpoly(rho**(2 * t - 1), sp.expand(P))
    g = [sp.expand(D(Phi, i)) for i in range(3)]
    L = sp.expand(sum(D(g[i], i) for i in range(3)))
    return g, L


def tau_l(l):
    return sp.Rational(3 + 2 * l, 4)


def Mmat(P, l, t=None):
    """M = oint_(S^2) f (x) n with f = grad Phi Lap Phi, exactly."""
    g, L = fields(P, tau_l(l) if t is None else t)
    return sp.Matrix(3, 3, lambda i, j: gfix(sph(unit(g[i] * L * X[j]))))


def closed(H, m):
    B = sp.Matrix(3, 3, lambda i, j: sph(sp.diff(H, X[i]) * sp.diff(H, X[j])))
    N = sph(H**2)
    return sp.simplify(8 / (sp.pi * (2 * m + 3)**2) * (N / 2 * sp.eye(3) - sp.Rational(4, 4 * m * m - 1) * B))


def harm_rand(l, seed):
    """Random real harmonic homogeneous polynomial of degree l: sums of Re((u+iv).x)^l, u.u=v.v, u.v=0."""
    rnd = random.Random(seed); P = 0
    if l == 0:
        return sp.Integer(rnd.randint(1, 3))
    while P == 0:
        for _ in range(3):
            u = sp.Matrix([rnd.randint(-3, 3) for _ in range(3)]); w = sp.Matrix([rnd.randint(-3, 3) for _ in range(3)])
            v = u.cross(w)
            if v.norm() == 0 or u.norm() == 0:
                continue
            v = v * u.norm() / v.norm()
            P += rnd.randint(1, 3) * sp.re(sp.expand(sum((u[i] + sp.I * v[i]) * X[i] for i in range(3))**l))
        P = sp.expand(P)
    return P


def lap(P):
    return sp.expand(sum(sp.diff(P, v, 2) for v in X))


print("=== Part 0: marginal exponent, existence, kernel constants ===")
check(sp.simplify((3 + 2 * lsym) / sp.Integer(4) - (1 + lsym) / sp.Integer(2)) == sp.Rational(1, 4),
      "tau_l - (1+l)/2 = 1/4 for every l: the literal (-Lap)^(-tau)Psi diverges at every marginal value (finite part used)")
check(all((tau_l(l) * 2).q != 1 for l in range(40)),
      "tau_l = (3+2l)/4 is never an integer or half-integer (l<40; 3+2l is odd): pure power kernel, no r^(2j) ln r")
check(all(4 * ((1 + l) / sp.Integer(2)) - 3 - 2 * l == -1 for l in range(40)),
      "where the literal reading exists (tau<(1+l)/2) the AX.1 exponent 4tau-3-2l is < -1, so G -> 0 under dilation")
rr_ = sp.Symbol("r", positive=True)
for t in [sp.Rational(7, 4), sp.Rational(9, 4), sp.Rational(11, 4), sp.Rational(13, 4)]:
    K = Ct(t) * rr_**(2 * t - 1)
    lhs = -(sp.diff(K, rr_, 2) + 2 / rr_ * sp.diff(K, rr_))
    check(gfix(lhs - Ct(t - 1) * rr_**(2 * t - 3)) == 0, f"-Lap K_tau = K_(tau-1) at tau={t} (so |xi|^2j H(d) K_tau_l = (-1)^j H(d) K_tau_(l-2j))")
for l in range(0, 9):
    cl = sp.prod([sp.Rational(2 * l + 1, 2) - 2 * k for k in range(l)])
    lhs = gfix(Ct(tau_l(l)) * cl)
    rhs = (-1)**l * 4 * sp.sqrt(2) / (sp.sqrt(sp.pi) * (2 * l - 1) * (2 * l + 3))
    kap = gfix(Ct(tau_l(l))**2 * cl**2 * (sp.Rational(1, 2) - l) * (l + sp.Rational(3, 2)))
    check(sp.simplify(lhs - rhs) == 0 and sp.simplify(kap + 8 / (sp.pi * (2 * l - 1) * (2 * l + 3))) == 0,
          f"l={l}: C_tau_l c_l = (-1)^l 4 sqrt2/(sqrt(pi)(2l-1)(2l+3)), kappa_l = -8/(pi(2l-1)(2l+3)) = {sp.nsimplify(kap * sp.pi)}/pi")

print("\n=== Part 1: exact M for harmonic symbols (AX.3), l = 0..5 ===")
pp = sp.Matrix([1, 2, -1])
check(Mmat(sp.Integer(1), 0) == sp.Rational(16, 9) * sp.eye(3), "l=0: M = (16/9) I  (Theorem AW, one charge)")
M1 = Mmat(-(x + 2 * y - z), 1)
check(sp.simplify(M1 - sp.Rational(16, 75) * ((pp.dot(pp)) * sp.eye(3) - 8 * pp * pp.T)) == sp.zeros(3),
      "l=1, p=(1,2,-1): M = (16/75)(|p|^2 I - 8 p p^T)  (Theorem AW, dipole)")
for l in range(1, 6):
    H = harm_rand(l, 7 + l)
    M = Mmat(H, l); Mc = closed(H, l)
    trp = -4 * sph(H**2) / (sp.pi * (2 * l - 1) * (2 * l + 3))
    check(lap(H) == 0 and sp.simplify(M - Mc) == sp.zeros(3) and sp.simplify(M.trace() - trp) == 0 and sp.simplify(M - M.T) == sp.zeros(3),
          f"l={l}, random harmonic H ({len(sp.Poly(H, *X).terms())} terms): sphere integral = closed form, tr M = -4N/(pi(2l-1)(2l+3)) = {sp.nsimplify(M.trace())} < 0")
Qz = 2 * z**2 - x**2 - y**2
MQ = Mmat(Qz, 2)
check(MQ == sp.diag(sp.Rational(64, 2205), sp.Rational(64, 2205), sp.Rational(-1472, 2205)),
      "axial quadrupole 2z^2-x^2-y^2: M = diag(64, 64, -1472)/2205, so G = (64 c_x, 64 c_y, -1472 c_z)/2205 on every radius")
check(Mmat(x * y, 2) == sp.diag(sp.Rational(-16, 441), sp.Rational(-16, 441), sp.Rational(16, 735)),
      "quadrupole xy: M = diag(-16/441, -16/441, 16/735)")

print("\n=== Part 1b: non-harmonic symbols; harmonic degrees decouple ===")
q = sp.symbols("q0:6")
Asym = sp.Matrix([[q[0], q[3], q[5]], [q[3], q[1], q[4]], [q[5], q[4], q[2]]])
Pg = sp.expand((sp.Matrix(X).T * Asym * sp.Matrix(X))[0])
Mg = Mmat(Pg, 2)
A0 = Asym - Asym.trace() / 3 * sp.eye(3)
pred = sp.Rational(16, 81) * Asym.trace()**2 * sp.eye(3) + sp.Rational(32, 735) * (A0 * A0).trace() * sp.eye(3) - sp.Rational(512, 2205) * A0 * A0
check(sp.simplify(Mg - pred) == sp.zeros(3),
      "general quadrupole xi.A xi: M = (16/81)(trA)^2 I + (32/735)tr(A0^2) I - (512/2205)A0^2 (6 free coefficients, exact)")
check(Mmat(R2, 2) == sp.Rational(16, 9) * sp.eye(3),
      "P=|xi|^2 (Psi = -4 pi delta: zero exterior Coulomb field) gives M = 16/9 I, the unit-charge value")
a1, a2, be = sp.symbols("a1 a2 beta", real=True)
Pd = sp.expand(a1 * x**2 + a2 * y**2 - (a1 + a2) * z**2 + be * R2)
Md = Mmat(Pd, 2)
gb = sp.groebner([sp.numer(sp.together(Md[i, i])) for i in range(3)], a1, a2, be, order="lex")
gbe = [sp.expand(g_) for g_ in gb.exprs]
print("       Groebner basis:", gbe)
check(Md.is_diagonal() and be**4 in gbe and any(sp.expand(g_ - (a1**2 - 35 * be**2 / 2) * sp.Poly(g_, a1, a2, be).LC()) == 0 for g_ in gbe),
      "l=2, diagonalised A: Groebner basis of M=0 contains beta^4 and a1^2-35beta^2/2, so M=0 only for P=0")
for l, seed in [(3, 1), (4, 2), (5, 3)]:
    Ha, Hb = harm_rand(l, 30 + seed), harm_rand(l - 2, 40 + seed)
    Hc = harm_rand(l - 4, 50 + seed) if l >= 4 else sp.Integer(0)
    Mt = Mmat(sp.expand(Ha + R2 * Hb + R2**2 * Hc), l)
    Mpre = closed(Ha, l) + closed(Hb, l - 2) + (closed(Hc, l - 4) if l >= 4 else sp.zeros(3))
    check(sp.simplify(Mt - Mpre) == sp.zeros(3), f"l={l}: M[H_l + |xi|^2 H_(l-2)" + (" + |xi|^4 H_(l-4)]" if l >= 4 else "]") + " = sum of M_m[H_m] (no cross terms)")

Pn4 = sp.expand(harm_rand(4, 61) + R2 * harm_rand(2, 62))
trn = Mmat(Pn4, 4).trace()
check(trn < 0, f"l=4 symbol H_4 + |xi|^2 H_2 (no |xi|^4 part): tr M = {sp.nsimplify(trn)} < 0, so M != 0")

print("\n=== Part 1c: M can be singular (l=2) and can vanish (l=4) ===")
tt = (-13 + sp.sqrt(69)) / 10
Psg = sp.expand(x**2 + tt * y**2 - (1 + tt) * z**2)
Msg = sp.simplify(Mmat(Psg, 2))
check(sp.simplify(Msg[2, 2]) == 0 and Msg.is_diagonal() and sp.simplify(Msg[0, 0]) != 0,
      "harmonic quadrupole x^2+t y^2-(1+t)z^2, t=(-13+sqrt69)/10: M_zz = 0, so G = 0 on every sphere centred on the z-axis")
H4 = sp.expand(x**4 + y**4 + z**4 - sp.Rational(3, 5) * R2**2)
check(lap(H4) == 0 and Mmat(H4, 4) == sp.Rational(-256, 121275) * sp.eye(3), "cubic hexadecapole H4: M = -(256/121275) I")
bstar = 4 * sp.sqrt(11) / 385
Pstar = sp.expand(H4 + bstar * R2**2)
check(sp.simplify(Mmat(Pstar, 4)) == sp.zeros(3),
      "P* = H4 + (4 sqrt11/385)|xi|^4: M = 0, so G = 0 on EVERY sphere enclosing the point (all centres, all radii)")

print("\n=== Part 1d: exact P_l + P_(l+1) pieces (AX.5) and Kept Failure AX.1 witness ===")


def pair_data(Pa, Pb, t):
    ga, La = fields(Pa, t); gb_, Lb = fields(Pb, t)
    g0 = [gfix(sph(unit(sum(ga[k] * X[k] for k in range(3)) * gb_[i] + sum(gb_[k] * X[k] for k in range(3)) * ga[i]
                        - sum(ga[k] * gb_[k] for k in range(3)) * X[i]))) for i in range(3)]
    v = [gfix(sph(unit(ga[i] * Lb + gb_[i] * La))) for i in range(3)]
    Mb = sp.Matrix(3, 3, lambda i, j: gfix(sph(unit(gb_[i] * Lb * X[j]))))
    return sp.Matrix(g0), sp.Matrix(v), Mb


t4 = tau_l(4)
g0x, v_x, M5x = pair_data(Pstar, x**5, t4)
G0EX = 64 * sp.sqrt(11) / 8085 - sp.Rational(5888, 266805)
check(sp.simplify(g0x[0] - G0EX) == 0 and g0x[1] == 0 and g0x[2] == 0 and v_x == sp.zeros(3, 1),
      f"P* + xi_x^5 at tau=11/4: g0 = (64 sqrt11/8085 - 5888/266805) e_x = {float(G0EX):.12f} e_x, v = 0")
gq, vq, Mq = pair_data(Qz, x**3 - 3 * x * y**2 + z**3, tau_l(2))
check(any(c_ != 0 for c_ in gq), f"quadrupole 2z^2-x^2-y^2 + octupole x^3-3xy^2+z^3 at tau=7/4: g0 = {[float(c_) for c_ in gq]} (shifts the zero of G off the device)")

# ------------------------------------------------------------- numerics
NU, NP = 160, 320
gu, wu = np.polynomial.legendre.leggauss(NU)
PH = 2 * np.pi * np.arange(NP) / NP
U = gu[:, None] * np.ones(NP); SN = np.sqrt(1 - U**2)
NRM = np.stack([SN * np.cos(PH), SN * np.sin(PH), U], -1)
WG = wu[:, None] * np.full(NP, 2 * np.pi / NP)


def gradfun(P, t):
    g, _ = fields(P, t)
    f = sp.lambdify((x, y, z, rho), g, "numpy")

    def F(p):
        r = np.linalg.norm(p, axis=-1)
        return np.stack([np.broadcast_to(np.asarray(o, float), r.shape) for o in f(p[..., 0], p[..., 1], p[..., 2], r)], -1)
    return F


def Tn(E, n):
    return np.sum(E * n, -1)[..., None] * E - 0.5 * np.sum(E * E, -1)[..., None] * n


def sphere(F, c, R):
    pts = np.asarray(c, float) + R * NRM
    return np.einsum("ij,ijk->k", R**2 * WG, Tn(F(pts), NRM))


def cube(F, c, h, n=96):
    g, w = np.polynomial.legendre.leggauss(n); G = np.zeros(3)
    for ax in range(3):
        o = [a for a in range(3) if a != ax]
        for sgn in (1, -1):
            pts = np.zeros((n, n, 3)); pts[..., ax] = sgn * h; pts[..., o[0]] = h * g[:, None]; pts[..., o[1]] = h * g[None, :]
            nv = np.zeros(3); nv[ax] = sgn
            G += np.einsum("ij,ijk->k", h * h * w[:, None] * w[None, :], Tn(F(pts + np.asarray(c, float)), np.broadcast_to(nv, pts.shape)))
    return G


def ax1(P, t, c, R):
    """AX.1: int_(S^2) f(w)(rho^k - R^k)/k dw on the quadrature grid."""
    l = sp.Poly(P, *X).total_degree()
    g, L = fields(P, t)
    f = sp.lambdify((x, y, z), [unit(gi * L) for gi in g], "numpy")
    W = np.stack([np.broadcast_to(np.asarray(o, float), U.shape) for o in f(NRM[..., 0], NRM[..., 1], NRM[..., 2])], -1)
    c = np.asarray(c, float); s = NRM @ c
    rh = s + np.sqrt(R * R - c @ c + s * s); k = float(4 * t - 2 - 2 * l)
    return np.einsum("ij,ijk->k", WG * (rh**k - R**k) / k, W)


print("\n=== Part 2: surface quadrature on off-centre spheres, pure point sources ===")
CASES = [((0.4, -0.3, 0.5), 1.5), ((0.4, -0.3, 0.5), 6.0), ((0.8, -0.6, 1.0), 3.0), ((0.4, -0.3, 0.5), 40.0)]
REC = {}
for lab, P, l in [("2z^2-x^2-y^2", Qz, 2), ("x^2 (non-harmonic)", x**2, 2), ("random harmonic l=2", harm_rand(2, 9), 2),
                  ("random harmonic l=3", harm_rand(3, 10), 3), ("l=3 non-harmonic x^2 y + y z^2 - x^3", x**2 * y + y * z**2 - x**3, 3)]:
    F = gradfun(P, tau_l(l)); Mn = np.array(Mmat(P, l).evalf(), float)
    errs = []
    for c, R in CASES:
        G = sphere(F, c, R); REC[(lab, c, R)] = G
        errs.append(np.max(np.abs(G - Mn @ np.array(c))) / np.max(np.abs(Mn @ np.array(c))))
    check(max(errs) < 1e-10, f"{lab} at tau={tau_l(l)}: G = M c on R = 1.5, 6, 40 (and 2c on R=3), rel. err {max(errs):.1e}; "
          f"G(c=(0.4,-0.3,0.5)) = {np.array2string(REC[(lab, CASES[0][0], CASES[0][1])], precision=9)}")
F = gradfun(Qz, tau_l(2))
check(np.max(np.abs(sphere(F, (0, 0, 0), 2.0))) < 1e-14, "centred sphere: G = 0 (parity)")
for t in [sp.Rational(8, 5), sp.Rational(19, 10)]:
    F = gradfun(Qz, t); c = (0.3, 0.2, -0.5)
    e1 = max(np.max(np.abs(sphere(F, c, R) - ax1(Qz, t, c, R))) / np.max(np.abs(sphere(F, c, R))) for R in (1.2, 4.0))
    v = [sphere(F, (0, 0, 1.0), R)[2] for R in (20.0, 40.0)]
    check(e1 < 1e-10 and abs(v[1] / v[0] - 2**(4 * float(t) - 7)) < 0.01,
          f"quadrupole off the marginal value, tau={t}: AX.1 integral matches quadrature (rel {e1:.1e}); R->2R ratio {v[1] / v[0]:.4f} vs 2^(4tau-7) = {2**(4 * float(t) - 7):.4f}")
Fsg = gradfun(Psg, tau_l(2))
check(max(np.max(np.abs(sphere(Fsg, (0, 0, b_), R))) for b_, R in [(0.7, 1.5), (2.0, 5.0)]) < 1e-13,
      "singular quadrupole: G = 0 on spheres centred at (0,0,0.7) and (0,0,2)")

print("\n=== Part 3: P* (M=0): zero on every sphere, nonzero and growing on cubes ===")
Fs = gradfun(Pstar, t4)
zs = max(np.max(np.abs(sphere(Fs, c, R))) for c, R in [((0.3, -0.2, 0.5), 2.0), ((0.3, -0.2, 0.5), 9.0), ((1.5, 0.4, -1.0), 3.0)])
check(zs < 1e-14, f"P* on three off-centre spheres: |G| <= {zs:.1e}")
CUBE = [cube(Fs, (0.3 * s_, -0.2 * s_, 0.5 * s_), 2.0 * s_) for s_ in (1, 2, 4)]
check(np.max(np.abs(CUBE[1] - 2 * CUBE[0])) < 1e-12 and np.max(np.abs(CUBE[2] - 4 * CUBE[0])) < 1e-12 and np.max(np.abs(CUBE[0])) > 1e-4,
      f"P* on cubes of half-side 2, 4, 8 (centre scaled with them): G = {np.array2string(CUBE[0], precision=8)} x 1, 2, 4 (grows like lambda)")

print("\n=== Part 4: Kept Failure AX.1, P* + xi_x^5 ===")
F5 = gradfun(x**5, t4)


def Ft(p):
    return Fs(p) + F5(p)


g0n = np.array([float(c_) for c_ in g0x]); M5n = np.array(M5x.evalf(), float)
cen = [sphere(Ft, (0, 0, 0), R) for R in (1.0, 3.0, 30.0)]
check(max(np.max(np.abs(cv - g0n)) for cv in cen) < 1e-14, f"centred spheres R = 1, 3, 30: G = {np.array2string(cen[1], precision=12)} = g0 (radius-free, NOT zero)")
ex = []
for c, R in [((0.5, 1.0, -0.7), 3.0), ((0.5, 1.0, -0.7), 12.0), ((-3.0, 2.0, 1.0), 6.0), ((-3.0, 2.0, 1.0), 600.0)]:
    c = np.array(c); G = sphere(Ft, c, R); pr = g0n + M5n @ c / (R * R - c @ c)
    ex.append((np.max(np.abs(G - pr)), G))
check(max(e for e, _ in ex) < 1e-14, "off-centre spheres: G = g0 + M' c/(R^2-|c|^2) exactly (AX.5 with M=0, v=0)")
check(np.max(np.abs(ex[3][1] - g0n)) < 3e-6 and np.max(np.abs(ex[2][1] - g0n)) > 1e-2,
      f"centre (-3,2,1), R=600: G = {np.array2string(ex[3][1], precision=8)} -> g0 for every centre")
CT5 = [cube(Ft, (0, 0, 0), h) for h in (2.0, 4.0, 8.0)]
check(max(np.max(np.abs(v_ - CT5[0])) for v_ in CT5) < 1e-12 and abs(CT5[0][0] - g0n[0]) > 1e-2,
      "same device on centred cubes h = 2, 4, 8: G = " + np.array2string(CT5[0], precision=10) + " on every size: a second radius-free value, not g0")
CO5 = [cube(Ft, (0.3 * s_, -0.2 * s_, 0.5 * s_), 2.0 * s_) for s_ in (1, 2, 4)]
DL = [CO5[i] - (1, 2, 4)[i] * CUBE[0] for i in range(3)]
check(np.max(np.abs((DL[0] - DL[1]) - 2 * (DL[1] - DL[2]))) < 1e-12 and np.max(np.abs(CUBE[0])) > 1e-4,
      "off-centre cubes (centre and side scaled by lambda = 1, 2, 4): G_x = " + ", ".join(f"{v_[0]:.8f}" for v_ in CO5)
      + " = lambda G[P*] + B + C/lambda exactly (degrees 1, 0, -1): grows without bound")

print("\n=== Part 5: extended devices (point-charge clusters) ===")
Kc, ex_ = float(Ct(tau_l(2))), float(2 * tau_l(2) - 1)


def clusterF(cl):
    def F(p):
        E = np.zeros(p.shape)
        for qq, a in cl:
            dd = p - np.asarray(a, float); r = np.linalg.norm(dd, axis=-1)
            E += qq * (Kc * ex_ * r**(ex_ - 2))[..., None] * dd
        return E
    return F


def moments(cl, n):
    return sp.expand(sum(sp.nsimplify(qq) * (-sum(sp.nsimplify(a[i]) * X[i] for i in range(3)))**n for qq, a in cl) / sp.factorial(n))


d = sp.Rational(1, 2)
CL = [(1.0, (0, 0, -0.5)), (-2.0, (0, 0, 0)), (1.0, (0, 0, 0.5))]
Mlin = np.array(Mmat(moments(CL, 2), 2).evalf(), float)
c0 = np.array((0.3, -0.4, 0.6)); er = [np.max(np.abs(sphere(clusterF(CL), c0, R) - Mlin @ c0)) for R in (5.0, 10.0, 20.0)]
check(moments(CL, 2) == d**2 * z**2 and er[0] > er[1] > er[2] and abs(er[1] / er[2] - 4) < 0.6,
      f"linear quadrupole (1,-2,1), P_2 = z^2/4 (non-harmonic), tau=7/4: G -> M c, error {', '.join(f'{e:.2e}' for e in er)} at R=5,10,20 (~R^-2)")
CL2 = [(1.0, (0.5, 0.1, 0.0)), (1.0, (-0.4, 0.3, 0.2)), (-1.0, (0.2, -0.5, 0.1)), (-1.0, (-0.3, 0.1, -0.3))]
CL3 = CL2 + [(-qq, (a[0], a[1], a[2] + 0.6)) for qq, a in CL2]
P3 = [moments(CL3, n) for n in range(4)]
check(P3[0] == 0 and P3[1] == 0 and P3[2] != 0 and P3[3] != 0, "8-charge device: charge 0, dipole 0, quadrupole and octupole nonzero")
M3 = np.array(Mmat(P3[2], 2).evalf(), float)
g03, v3, _ = pair_data(P3[2], P3[3], tau_l(2))
g03n = np.array([float(c_) for c_ in g03]); xs = -np.linalg.solve(M3, g03n)
F3 = clusterF(CL3)
e3 = [np.max(np.abs(sphere(F3, c0, R) - (M3 @ c0 + g03n))) for R in (6.0, 12.0, 24.0)]
z3 = [np.max(np.abs(sphere(F3, xs, R))) for R in (6.0, 12.0, 24.0)]
check(e3[0] > e3[1] > e3[2] and e3[2] < 2e-3 * np.max(np.abs(M3 @ c0 + g03n)),
      f"device at tau=7/4: G -> M c + g0, error {', '.join(f'{e:.2e}' for e in e3)} at R=6,12,24; g0 = {np.array2string(g03n, precision=8)}")
check(z3[0] > z3[1] > z3[2], f"spheres centred at x* = -M^-1 g0 = {np.array2string(xs, precision=6)}: |G| = {', '.join(f'{e:.2e}' for e in z3)} -> 0")
REC["cluster_centred_origin"] = sphere(F3, (0, 0, 0), 24.0)
print(f"       spheres centred at the origin, R=24: G = {np.array2string(REC['cluster_centred_origin'], precision=8)} (not 0)")

print("\n=== Recorded values (lock) ===")
print("RECORD", np.array2string(REC[("2z^2-x^2-y^2", CASES[0][0], CASES[0][1])], precision=12), np.array2string(CUBE[0], precision=12), np.array2string(g03n, precision=12))
check(abs(float(G0EX) - 0.004185497112987405) < 1e-15, f"g0 of the Kept Failure witness = {float(G0EX):.15f}")

if FAIL:
    print(f"\n{len(FAIL)} check(s) failed"); sys.exit(1)
print("\nall Theorem AX checks passed; nothing here is thrust")
