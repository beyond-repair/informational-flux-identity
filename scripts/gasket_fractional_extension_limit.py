#!/usr/bin/env python3
"""Theorem BG: the refinement limit of the L^alpha Dirichlet extension itself (not only of a scalar audit norm).

Setting (GASKET.md, Theorem BG): level-n combinatorial gasket V_n, corner data d on C = {a, b, c}, dbar = mean(d),
x = d - dbar, u_n the L_n^alpha-Dirichlet extension (Assumption 4), alpha in (0,1), rho = 5^alpha/3,
gamma = log3/log5.  Notation of Theorems AQ-AZ and BB-BD: r_n(s) = 1/2 v^T (L_n+s)^-1 v, p(s) = s(s+5),
psi_n(tau) = p^n(5^-n tau) increasing to psi, Y_n(tau) = 3^n y_n(5^-n tau) -> Y (BD.4),
R_alpha(n) = (2 sin pi a/pi) D_n, D_n = int s^-a r_n ds.  For x in V_m (m <= n):
  h_m(z; x) = [(L_m+z)^-1 x_C]_x / r_m(z)   (x_C = x on the corners, 0 elsewhere).

Theorem BG (checked here; exit 1 on any failure):
  (1) exact profile representation, every n >= m, every alpha, every x in V_m:
        u_n(x) = dbar + int s^-a r_n(s) h_m(p^(n-m)(s); x) ds / D_n,
      from [(L_n+s)^-1 x_C]|_{V_m} = (r_n(s)/r_m(p^(n-m)s)) [(L_m+p^(n-m)s)^-1 x_C] (AQ decimation + Woodbury
      on the zero-sum corner mode); h_{m+1}(z; .) = h_m(p(z); .) on V_m; h_m(0; .) + dbar is the harmonic
      extension, h_m(z; c) = x_c, |h_m(z; x)| <= 4 ||x||_inf/(4+z) off the corners (maximum principle);
      at level 1, h_1(z; m_ab) = -x_c/(5+z) (sympy).
  (2) flat pointwise limit for alpha <= gamma, with an explicit bound valid for every alpha, m, n >= m:
        |u_n(x) - dbar| <= (8/3) 4^-a 5^(m(1-a)) ||x||_inf rho^n / R_alpha(n)   for x in V_m \\ V_0,
      so u_n -> dbar on every vertex of V_* \\ V_0: geometrically (rho^n) for alpha < gamma, like 1/n at gamma.
  (3) the renormalised profile exists, every alpha: rho^-n int s^-a r_n h_m(p^(n-m)s; x) ds ->
        P_m(x) = int_0^inf tau^-a Rr(tau) h_m(psi(5^-m tau); x) dtau,  Rr = Y psi/tau = lim (3/5)^n r_n(5^-n tau),
      so for alpha < gamma  rho^-n (u_n(x) - dbar) -> P_m(x)/D_inf, and for alpha > gamma
      u_n(x) -> dbar + P_m(x)/(lim rho^-n D_n) (Theorem AZ: R_alpha(n) = A_alpha rho^n + C_alpha + o(1)): NOT flat.
      At the level-1 midpoint m_ab: P_1 = -x_c W(alpha), W = int tau^-a Rr(tau)/(5+psi(tau/5)) dtau.
      Certified (mpmath.iv, monotone cells, BD.4 tail factors, explicit psi tail):
        W(0.45), U(0.45) = W/D_inf (D_inf from Theorem AS, re-run live), W(0.8), kappa_0.8 = (2 sin(0.8 pi)/pi) W/A_0.8.
  (4) renormalised energy (DERIVED): rho^n Q_alpha(u_n) = 2 rho^n ||x||^2/R_alpha(n) -> 0 for alpha <= gamma,
      -> 2||x||^2/A_alpha > 0 for alpha > gamma.
Kept Failure BG.1: the flat limit is not uniform: the mean of the two neighbours of corner a equals
  1 - 5||F(n)||/(2 sqrt21) for data (1,-1/2,0) and tends to 1 - 5 L(0.45)/(2 sqrt21) (certified), so
  sup_{V_n \\ V_0} |u_n - dbar| stays above 0.3853 in the limit: a corner boundary layer, not a continuous limit.

Not a field solver.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_maxwell_limit as bb  # noqa: E402
import mesh_maxwell_residual as ba  # noqa: E402
import mesh_noncentrality_limit as bd  # noqa: E402
import gasket_fractional_limit_digits as asd  # noqa: E402

FAILS: list[str] = []
GAMMA = math.log(3) / math.log(5)
# Theorem AZ (3), certified A_0.8 = c_R(0.8)/(rho-1) (width < 1e-20), used as stated, padded outward
AZ_A08 = ("2.39640970397922052703", "2.39640970397922052705")


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def level1_graph():
    V = ["a", "b", "c", "ab", "bc", "ca"]
    E = [("a", "ab"), ("a", "ca"), ("b", "ab"), ("b", "bc"), ("c", "bc"), ("c", "ca"), ("ab", "bc"), ("bc", "ca"), ("ca", "ab")]
    idx = {v: i for i, v in enumerate(V)}
    L = sp.zeros(6)
    for i, j in E:
        I, J = idx[i], idx[j]
        L[I, I] += 1
        L[J, J] += 1
        L[I, J] -= 1
        L[J, I] -= 1
    return L


def symbolic_checks():
    z, s = sp.symbols("z s", positive=True)
    xa, xb = sp.symbols("x_a x_b")
    xc = -xa - xb
    L = level1_graph()
    M = L + z * sp.eye(6)
    w = M.LUsolve(sp.Matrix([xa, xb, xc, 0, 0, 0]))
    v = sp.Matrix([1, -1, 0, 0, 0, 0])
    r1 = sp.simplify((v.T * M.inv() * v)[0] / 2)
    rr = (2 + z) * (5 + z) * (1 / (3 + z * (z + 5))) / ((6 + z) - 2 * z * (5 + z) / (3 + z * (z + 5)))
    check(sp.simplify(r1 - rr) == 0, "level 1: r_1(z) = Phi_z(r_0(p(z))) = (z+5)/(z^2+7z+9)")
    check(all(sp.simplify(w[k] - r1 * [xa, xb, xc][k]) == 0 for k in range(3)),
          "level 1: corner block of (L_1+z)^-1 on the zero-sum plane is r_1(z) I, so h_1(z; c) = x_c")
    h = [sp.simplify(w[3 + k] / r1) for k in range(3)]
    check(sp.simplify(h[0] + xc / (5 + z)) == 0 and sp.simplify(h[1] + xa / (5 + z)) == 0 and sp.simplify(h[2] + xb / (5 + z)) == 0,
          "level 1: h_1(z; m_ab) = -x_c/(5+z), h_1(z; m_bc) = -x_a/(5+z), h_1(z; m_ca) = -x_b/(5+z)")
    check(sp.simplify(h[0].subs(z, 0) - (2 * xa + 2 * xb + xc) / 5) == 0,
          "z = 0: h_1(0; m_ab) + dbar = (2d_a+2d_b+d_c)/5, the harmonic (alpha = 1) value")
    # Woodbury step: S = c Z^-1 - e P_C, Z q = r q on the zero-sum corner mode  =>  S^-1 q = Z q/(c - e r)
    c = (6 + s) / ((5 + s) * (2 + s))
    e = 2 * s / (2 + s)
    check(sp.simplify(c * s * (s + 5) - (2 * e + s)) == 0,
          "c(s) p(s) = 2e(s) + s: AQ's Schur complement c L_n + (2e+s)I - eP_C = c(L_n+p(s)) - eP_C")
    r, q = sp.symbols("r q")
    check(sp.simplify((c * q - e * r * q) / (c - e * r) - q) == 0,
          "(c Z^-1 - e P_C)(Z q/(c - e r)) = (c q - e r q)/(c - e r) = q when P_C Z q = r q")
    P = sp.Symbol("P", positive=True)
    check(sp.simplify(1 / (5 + P) - P / (P * (P + 5))) == 0,
          "1/(5+psi(tau/5)) = psi(tau/5)/psi(tau) (psi(5t) = p(psi(t))), so W = int tau^(-1-a) Y(tau) psi(tau/5) dtau")
    a = sp.Rational(9, 20)
    K = sp.Integer(25)
    t = sp.Symbol("t", positive=True)
    val = sp.integrate(t ** (-a) / (4 + t / K), (t, 0, sp.oo))
    check(abs(float(val) - float(K * (4 * K) ** (-a) * sp.pi / sp.sin(sp.pi * a))) < 1e-10,
          "int_0^inf t^-a/(4+t/K) dt = K (4K)^-a pi/sin(pi a) (a = 0.45, K = 25)")
    xs = [sp.Rational(5, 6), sp.Rational(-2, 3), sp.Rational(-1, 6)]
    Rs = sp.Symbol("R", positive=True)
    check(sp.simplify((2 / Rs) ** 2 * (Rs / 2) * sum(x_ ** 2 for x_ in xs) - sp.Rational(7, 3) / Rs) == 0,
          "Q_a(u_n) = q^T G q = 2|x|^2/R_a(n) = 7/(3R_a(n)) for data (1,-1/2,0) (Theorem AD's pairing Q R = 7/3)")


# ----------------------------------------------------------------- (B) exact rationals on the meshes
def exact_solve(adj, s: Fraction, rhs):
    N = len(adj)
    A = [[Fraction(0)] * N + [rhs[i]] for i in range(N)]
    for i in range(N):
        A[i][i] = Fraction(len(adj[i])) + s
        for j in adj[i]:
            A[i][j] -= 1
    for k in range(N):
        piv = next(i for i in range(k, N) if A[i][k] != 0)
        A[k], A[piv] = A[piv], A[k]
        for i in range(N):
            if i != k and A[i][k] != 0:
                f = A[i][k] / A[k][k]
                A[i] = [x - f * y for x, y in zip(A[i], A[k])]
    return [A[i][N] / A[i][i] for i in range(N)]


def mesh(n):
    P, adj, C, lat = ba.gasket_setup(n)
    return adj, C, {lat[i]: i for i in range(len(adj))}


def corner_solve(m, z, xC):
    adj, C, pos = mesh(m)
    rhs = [Fraction(0)] * len(adj)
    for k, cc in enumerate(C):
        rhs[cc] = xC[k]
    return exact_solve(adj, z, rhs), pos, C


def exact_checks(levels=(1, 2, 3)):
    xC = (Fraction(3), Fraction(-1), Fraction(-2))  # a generic zero-sum corner vector, ||x||_inf = 3
    ok_rep = ok_bnd = ok_corner = True
    for n in levels:
        for s in (Fraction(1, 2), Fraction(7, 5)):
            wn, posn, Cn = corner_solve(n, s, xC)
            rn = bb.rec_all(n, s)[0]
            ok_corner = ok_corner and all(wn[cc] == rn * xC[k] for k, cc in enumerate(Cn))
            ok_bnd = ok_bnd and all(abs(wn[i] / rn) <= 4 * 3 / (4 + s) for i in range(len(wn)) if i not in Cn)
            for m in range(1, n):
                z = bb.orbit(s, n - m)[n - m]
                wm, pos_m, _ = corner_solve(m, z, xC)
                rm = bb.rec_all(m, z)[0]
                ok_rep = ok_rep and all(wn[posn[key]] == rn / rm * wm[i] for key, i in pos_m.items())
    check(ok_corner, f"exact, levels {levels[0]}-{levels[-1]}, s in {{1/2, 7/5}}: [(L_n+s)^-1 x_C]_c = r_n(s) x_c (zero-sum corner mode)")
    check(ok_bnd, "exact: |h_n(s; x)| <= 4 ||x||_inf/(4+s) at every non-corner vertex (maximum principle)")
    check(ok_rep, "exact: [(L_n+s)^-1 x_C]|_{V_m} = (r_n(s)/r_m(p^(n-m)s)) [(L_m+p^(n-m)s)^-1 x_C] for every m < n (BG.1)")


# ----------------------------------------------------------------- (C) dense floats: u_n itself
def dense_u(n, al, d=(1.0, -0.5, 0.0)):
    P, adj, C, lat = ba.gasket_setup(n)
    L = bb.laplacian(adj)
    w, V = np.linalg.eigh(L)
    keep = np.arange(len(w)) > 0
    w[0] = 1.0
    A = (V * np.where(keep, w ** al, 0.0)) @ V.T
    G = (V * np.where(keep, w ** -al, 0.0)) @ V.T
    R = G[C[0], C[0]] + G[C[1], C[1]] - 2 * G[C[0], C[1]]
    I = [i for i in range(len(adj)) if i not in C]
    u = np.zeros(len(adj))
    u[C] = d
    u[I] = np.linalg.solve(A[np.ix_(I, I)], -A[np.ix_(I, C)] @ np.array(d))
    I_a = 2 * u[C[0]] - sum(u[j] for j in adj[C[0]])
    return dict(u=u, pos={lat[i]: i for i in range(len(adj))}, C=C, adj=adj, I_a=I_a, R=R)


def mid_ratio(n, al, dps=30):
    """W_n, D_n with W_n = int s^-a r_n/(5+p^(n-1)(s)), D_n = int s^-a r_n (recursion quadrature, consistency only)."""
    mp.dps = dps
    a = mpf(al)
    pts = [0] + [mpf(10) ** k for k in range(-30, 31, 2)]
    X = pts[-1]

    def fW(s):
        return s ** (-a) * bb.rec_all(n, s)[0] / (5 + bb.orbit(s, n - 1)[n - 1])

    W = mp.quad(fW, pts)
    D = mp.quad(lambda s: s ** (-a) * bb.rec_all(n, s)[0], pts) + X ** (-a) / a
    return W, D


def dense_checks(levels=(1, 2, 3, 4, 5, 6), alphas=(0.45, 0.8), quad_levels=(1, 2, 4)):
    xs = (5 / 6, -2 / 3, -1 / 6)
    dbar = 1 / 6
    out = {}
    for al in alphas:
        rho = 5**al / 3
        ok_mid = ok_bnd = True
        for n in levels:
            r = dense_u(n, al)
            u, pos, C = r["u"], r["pos"], r["C"]
            mab = pos[(Fraction(1, 2), Fraction(0))]
            r["dev"] = u[mab] - dbar
            out[(al, n)] = r
            if n in quad_levels:
                W, D = mid_ratio(n, str(al))
                ok_mid = ok_mid and abs(r["dev"] - (-xs[2]) * float(W / D)) < 1e-10
            for m in range(1, min(n, 3) + 1):
                bound = 8 / 3 * 4 ** (-al) * 5 ** (m * (1 - al)) * max(abs(x) for x in xs) * rho**n / r["R"]
                for key in mesh(m)[2]:
                    i = pos[key]
                    if i not in C:
                        ok_bnd = ok_bnd and abs(u[i] - dbar) <= bound
        check(ok_mid, f"alpha={al}: dense u_n(m_ab) - 1/6 = (1/6) W_n/D_n (BG.1 with h_1 = -x_c/(5+z)) at levels {quad_levels}, |err| < 1e-10")
        check(ok_bnd, f"alpha={al}: |u_n(x) - dbar| <= (8/3)4^-a 5^(m(1-a)) ||x|| rho^n/R_a(n) at every x in V_m minus corners, "
                      f"m <= 3, levels {levels[0]}-{levels[-1]} (BG.2)")
    return out


# ----------------------------------------------------------------- (D) interval certificate for W(alpha)
def R_psi(N, tau):
    """Rr_N(tau) = (3/5)^N r_N(5^-N tau) and psi_{N-1}(tau/5) = p^(N-1)(5^-N tau), interval arithmetic along the orbit."""
    s = tau / iv.mpf(5) ** N
    o = [s]
    for _ in range(N):
        o.append(o[-1] * (o[-1] + 5))
    r = 1 / (3 + o[N])
    for k in range(N - 1, -1, -1):
        z = o[k]
        r = (2 + z) * (5 + z) * r / ((6 + z) - 2 * z * (5 + z) * r)
    return (iv.mpf(3) / 5) ** N * r, o[N - 1]


def eps_psi(N, tau):
    """psi_N <= psi <= psi_N e^(2^N 5^-N tau/10): 5^N p^k(5^(-N-k) tau) <= tau e^(5^-N tau/10), deg psi_N = 2^N."""
    return iv.mpf(2) ** N * iv.mpf(5) ** (-N) * iv.mpf(tau) / 10


def build_cells(N=60, tau0=1e-20, tau1=3125, r1=1.05, r2=1.001, split=0.01):
    """Ceilings of Rr(g)/(5+psi(g/5)) and floors at every grid point (alpha-free, reused for every alpha)."""
    iv.dps = 20
    grid = [tau0]
    while grid[-1] < split:
        grid.append(grid[-1] * r1)
    while grid[-1] < tau1:
        grid.append(min(grid[-1] * r2, tau1))
    up, lo = [], []
    for g in grid:
        t = iv.mpf(g)
        R, ps = R_psi(N, t)
        eY = bd.epsY(N, t)
        # Rr = Y psi/tau: Y within e^(+-epsY) of Y_N (BD.4), psi_N <= psi <= psi_N e^(eps_psi)
        up.append(R * iv.exp(eY + eps_psi(N, t)) / (5 + ps))
        lo.append(R * iv.exp(-eY) / (5 + ps * iv.exp(eps_psi(N - 1, t / 5))))
    return dict(grid=grid, up=up, lo=lo, N=N, tau0=tau0, tau1=tau1)


def tail_P(t1, M=6, k=8):
    """psi_M(t) >= (t/t1)^k P_k for t >= t1, P_k = sum_{j>=k} c_j t1^j (psi_M has nonnegative coefficients)."""
    xx = sp.Symbol("x")
    q = xx
    for _ in range(M):
        q = sp.expand(q * (q + 5))
    q = sp.Poly(q.subs(xx, xx / sp.Integer(5) ** M), xx)
    Pk = sp.Rational(q.eval(sp.Integer(t1))) - sum(q.coeff_monomial(xx**j) * sp.Integer(t1) ** j for j in range(1, k))
    return iv.mpf(int(Pk.p)) / iv.mpf(int(Pk.q))


def W_bracket(cells, al, M=6, k=8):
    """Certified W(a) = int_0^inf tau^-a Rr(tau)/(5+psi(tau/5)) dtau: tau^-a times a nonincreasing positive function."""
    iv.dps = 20
    a = iv.mpf(al)
    grid, up, lo = cells["grid"], cells["up"], cells["lo"]
    U = Lw = iv.mpf(0)
    for i in range(len(grid) - 1):
        u, v = iv.mpf(grid[i]), iv.mpf(grid[i + 1])
        w = (v ** (1 - a) - u ** (1 - a)) / (1 - a)
        U += up[i] * w
        Lw += lo[i + 1] * w
    h = iv.mpf(cells["tau0"]) ** (1 - a) / (1 - a)
    U += h / 15                      # head: Rr <= 1/3, 1/(5+psi) <= 1/5
    Lw += lo[0] * h                  # head: monotone floor at tau0
    t1 = cells["tau1"]
    Pk = tail_P(t1 // 5, M, k)       # tail: Rr <= 1/3, psi(tau/5) >= psi_M(tau/5) >= (tau/tau1)^k P_k(tau1/5)
    U += iv.mpf(t1) ** (1 - a) / (3 * Pk * (k - 1 + a))
    return iv.mpf([Lw.a, U.b]), len(grid)


# ----------------------------------------------------------------- main
def main(r2=1.001) -> int:
    t0 = time.time()
    print("Theorem BG: refinement limit of the L^alpha Dirichlet extension on the gasket. Not thrust.")
    symbolic_checks()
    exact_checks()
    dd = dense_checks()

    # Kept Failure BG.1 on the meshes: corner-neighbour mean = 1 - I_a/2 = 1 - 5||F||/(2 sqrt21), increasing in n
    means = [sum(dd[(0.45, n)]["u"][j] for j in dd[(0.45, n)]["adj"][dd[(0.45, n)]["C"][0]]) / 2 for n in range(1, 7)]
    check(all(abs(means[n - 1] - (1 - dd[(0.45, n)]["I_a"] / 2)) < 1e-12 for n in range(1, 7))
          and all(means[i] < means[i + 1] for i in range(5)) and means[0] - 1 / 6 > 0.04,
          f"alpha=0.45, levels 1-6: corner-neighbour mean = {', '.join(f'{m:.6f}' for m in means)}: increasing, > dbar + 0.04")

    # recursion quadrature: renormalised midpoint deviation
    tab = {}
    for al in ("0.45", "0.8", repr(GAMMA)):
        rho = mpf(5) ** mpf(al) / 3
        for n in (10, 20, 30):
            W, D = mid_ratio(n, al)
            tab[(al, n)] = (W / rho**n, W / D)
            print(f"   a={al[:6]} n={n}: rho^-n W_n = {mp.nstr(W / rho**n, 14)}, W_n/D_n = {mp.nstr(W / D, 12)}, "
                  f"rho^-n W_n/D_n = {mp.nstr(W / D / rho**n, 12)}, n W_n/D_n = {mp.nstr(n * W / D, 8)}")
    check(abs(tab[("0.45", 30)][0] - tab[("0.45", 20)][0]) < 1e-9 and abs(tab[("0.8", 30)][0] - tab[("0.8", 20)][0]) < 1e-9,
          "rho^-n W_n settles (BG.3 dominated convergence) at 0.45 and 0.8")
    g = repr(GAMMA)
    check(tab[(g, 10)][1] > tab[(g, 20)][1] > tab[(g, 30)][1] > 0 and 30 * tab[(g, 30)][1] < 1.5 * 10 * tab[(g, 10)][1],
          "alpha = gamma: W_n/D_n -> 0 with n W_n/D_n bounded (flat limit at rate 1/n)")

    # certified W(0.45), W(0.8); D_inf(0.45) from Theorem AS (re-run live)
    t1 = time.time()
    cells = build_cells(r2=r2)
    W45, ncell = W_bracket(cells, "0.45")
    W80, _ = W_bracket(cells, "0.8")
    print(f"   W(0.45) in [{float(W45.a):.10f}, {float(W45.b):.10f}], W(0.8) in [{float(W80.a):.10f}, {float(W80.b):.10f}] "
          f"({ncell} cells, {time.time() - t1:.1f}s)")
    check(W45.a < tab[("0.45", 30)][0] < W45.b and W80.a < tab[("0.8", 30)][0] < W80.b and W45.a > 0 and W80.a > 0,
          "W(0.45), W(0.8) certified > 0; the level-30 recursion values lie inside")
    iv.dps = 50
    mp.dps = 50
    B, okB, _ = asd.power_ceiling()
    check(okB, "Theorem AS power ceiling re-certified")
    qn, s0, vals, N = asd.build_nodes()
    cAS = asd.certify_limit("0.45", qn, s0, vals, N, B)
    iv.dps = 20
    Dinf = iv.mpf([cAS["Jr"].a, cAS["Jr"].b])
    Lim = iv.mpf([cAS["L"].a, cAS["L"].b])
    U45 = W45 / Dinf
    rho45 = mpf(5) ** mpf("0.45") / 3
    print(f"   D_inf(0.45) = int s^-a r_inf in {asd.fmt(Dinf, 14)};  U(0.45) = W/D_inf in [{float(U45.a):.9f}, {float(U45.b):.9f}]")
    check(U45.a < tab[("0.45", 30)][1] / rho45**30 < U45.b,
          "U(0.45) = lim rho^-n (u_n(m_ab) - dbar)/(-x_c) certified; the level-30 value lies inside")
    dev6 = dd[(0.45, 6)]["dev"] / (5**0.45 / 3) ** 6
    print(f"   dense level 6: rho^-6 (u_6(m_ab) - 1/6) = {dev6:.6f}; limit U/6 in [{float(U45.a)/6:.6f}, {float(U45.b)/6:.6f}]")
    A08 = iv.mpf([AZ_A08[0], AZ_A08[1]])
    a8 = iv.mpf("0.8")
    kap = 2 * iv.sin(iv.pi * a8) / iv.pi * W80 / A08
    uinf = iv.mpf(1) / 6 + kap / 6
    print(f"   kappa_0.8 in [{float(kap.a):.9f}, {float(kap.b):.9f}];  lim u_n(m_ab) at 0.8, data (1,-1/2,0): "
          f"[{float(uinf.a):.9f}, {float(uinf.b):.9f}] (harmonic: 0.2)")
    check(kap.a > 0.08 and kap.b < 0.2 and float(kap.a) < float(tab[("0.8", 30)][1]),
          "alpha = 0.8 > gamma: lim u_n(m_ab) = dbar - kappa x_c with certified kappa_0.8 > 0.08 (NOT flat), < 1/5 (harmonic); "
          "level 30 above the limit (C_0.8 < 0)")
    gap = iv.mpf(5) / 6 - 5 * Lim / (2 * iv.sqrt(iv.mpf(21)))
    print(f"   Kept Failure BG.1: lim (corner-neighbour mean - dbar) = 5/6 - 5L/(2 sqrt21) in [{float(gap.a):.10f}, {float(gap.b):.10f}]")
    check(gap.a > 0.3853 and means[-1] - 1 / 6 < gap.a,
          "Kept Failure BG.1: sup over V_n minus corners of |u_n - dbar| does not tend to 0 (certified liminf > 0.3853); "
          "levels 1-6 approach from below")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BG checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
