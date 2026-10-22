#!/usr/bin/env python3
"""Theorem BB: exact decimation of the mesh residual constant kappa_n(alpha) of Theorem BA.

Setting: Assumption BA (README).  Level-n gasket, A = L^alpha, G = A^+, Psi = u the A-Dirichlet
extension of corner data d, q = A u on the corners, one-sided corner-cell gradient g_c, f_c = -q_c g_c,
D = sum_c f_c = kappa_n(alpha) Q(d).  Notation of Theorems AQ-AS (GASKET.md):
  r_n(s) = 1/2 v^T (L_n+s)^-1 v (v = e_a - e_b),  R_n = v^T G v,  sigma_n = 1/2 v^T L^(1-alpha) v,
  p(s) = s(s+5),  r_{n+1}(s) = Phi_s(r_n(p(s))),  r_0 = 1/(3+s).
For a corner c with neighbours p' (towards c') and p'' (towards c''):
  Y_n = G_{p'c'} - G_{p'c''},   E_n = G_{p'c'} + G_{p'c''} - G_{cc'} - G_{cc''}.

Theorem BB (checked here; exit 1 on any failure):
  (1) kappa_n = 2^n (4/R_n^2)(3 Y_n - sigma_n);  self part kappa_n^self = -2^n (4/R_n^2)(sigma_n - E_n);
      cross part kappa_n^cross = kappa_n - kappa_n^self = 2^n (4/R_n^2)(3 Y_n - E_n).
  (2) y_n(s) := [(L_n+s)^-1]_{p'c'} - [(L_n+s)^-1]_{p'c''} = r_n(s) prod_{j<n} 1/(5+p^j(s)),
      so 0 < Y_n <= 5^-n R_n / 2.
  (3) e_n(s) := resolvent version of E_n = s g_n(s), g_n(s) = [(L_n+s)^-1]_{cc'};
      g_{n+1}(s) = F_n(s) g_n(p(s)), F_n = (2+s)(5+s)(6+s)/((6+s-2t_n(ps))(6+s-2p(s) r_n(ps))),
      t_n(s) = s [(L_n+s)^-1 1_C]_c, t_{n+1}(s) = (2+s) t_n(ps)/(6+s-2 t_n(ps)), t_0 = 1, g_0 = 1/(s(3+s));
      F_n <= (3/20)(5+s)^2, hence g_n(s) <= (3/20)^n / s^2 and
      E_n <= (sin pi a/pi) (3/20)^(n(1-a)) (1/(3(1-a)) + 1/a).
  (4) Conjecture BA.2 holds: for alpha < log3/log5, 2^-n kappa_n -> K(alpha) = -4 sigma_inf / R_inf^2 < 0;
      certified K(0.45) from Theorem AS's enclosures (re-certified live).  For alpha >= log3/log5 the
      limit is 0 (R_n -> inf).
  (5) for alpha < 1 - log2/log(20/3) (~0.6346), kappa_n^cross -> 0: the whole 2^n growth and the whole
      limit are the self-force term.  Kept Failure BB.1: kappa_n^cross(0.45) > 0 certified at n = 1..4.

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
import mesh_maxwell_residual as ba  # noqa: E402
import gasket_fractional_limit_digits as asd  # noqa: E402

FAILS: list[str] = []
GAMMA = math.log(3) / math.log(5)
ALPHA_STAR = 1 - math.log(2) / math.log(20 / 3)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- graph helpers
def corner_frame(n):
    """Positions, adjacency, corners (a,b,c), and for corner a its neighbours toward b and toward c."""
    P, adj, C, lat = ba.gasket_setup(n)
    a, b, c = C
    nb = adj[a]
    pb = min(nb, key=lambda j: np.linalg.norm(P[j] - P[b]))
    pc = [j for j in nb if j != pb][0]
    return P, adj, C, pb, pc


def laplacian(adj):
    N = len(adj)
    L = np.zeros((N, N))
    for i in range(N):
        L[i, i] = len(adj[i])
        for j in adj[i]:
            L[i, j] -= 1
    return L


# ----------------------------------------------------------------- recursions (any number type)
def orbit(s, n):
    o = [s]
    for _ in range(n):
        o.append(o[-1] * (o[-1] + 5))
    return o


def rec_all(n, s):
    """(r_n, t_n, g_n, y_n) at s from the level-0 data along the orbit s, p(s), ..., p^n(s)."""
    o = orbit(s, n)
    x = o[n]
    r, t, g, yr = 1 / (3 + x), 1, 1 / (x * (3 + x)), 1
    for k in range(n - 1, -1, -1):
        z, pz = o[k], o[k + 1]
        F = (2 + z) * (5 + z) * (6 + z) / ((6 + z - 2 * t) * (6 + z - 2 * pz * r))
        g = F * g
        r = (2 + z) * (5 + z) * r / ((6 + z) - 2 * z * (5 + z) * r)
        t = (2 + z) * t / (6 + z - 2 * t)
        yr = yr / (5 + z)
    return r, t, g, yr * r


# ----------------------------------------------------------------- (A) symbolic steps
def symbolic_checks():
    s, t, x, X = sp.symbols("s t x X", positive=True)
    c = (6 + s) / ((5 + s) * (2 + s))
    e = 2 * s / (2 + s)
    ps = s * (s + 5)
    # Woodbury mode factor on a corner mode with corner-block eigenvalue lam at p(s): lam -> lam/(c - e lam)
    # zero-sum mode: lam = r = x/p(s);  1_C mode: lam = u = t/p(s)
    den_r = sp.factor(sp.simplify(c - e * x / ps))
    den_u = sp.factor(sp.simplify(c - e * t / ps))
    check(sp.simplify(den_r - (6 + s - 2 * x) / ((2 + s) * (5 + s))) == 0,
          "c - e r_n(ps) = (6+s-2 p(s) r_n(ps))/((2+s)(5+s))  (AQ's den)")
    check(sp.simplify(den_u - (6 + s - 2 * t) / ((2 + s) * (5 + s))) == 0,
          "c - e u_n(ps) = (6+s-2 t_n(ps))/((2+s)(5+s)),  t = s u")
    # t recursion: t_{n+1}(s) = s * u/(c - e u), u = t/p(s)
    tn1 = sp.simplify(s * (t / ps) / (c - e * t / ps))
    check(sp.simplify(tn1 - (2 + s) * t / (6 + s - 2 * t)) == 0, "t_{n+1}(s) = (2+s) t_n(ps)/(6+s-2t_n(ps))")
    # g = (u - r)/3:  g_{n+1} = c (u'-r')/(3 (c-eu')(c-er'))
    U, Rr = sp.symbols("U Rr", positive=True)
    lhs = (U / (c - e * U) - Rr / (c - e * Rr)) / 3
    rhs = c * (U - Rr) / (3 * (c - e * U) * (c - e * Rr))
    check(sp.simplify(lhs - rhs) == 0, "g_{n+1}(s) = c g_n(ps)/((c-eu')(c-er'))  (corner block, both modes)")
    F = sp.simplify(c / (den_u * den_r))
    check(sp.simplify(F - (2 + s) * (5 + s) * (6 + s) / ((6 + s - 2 * t) * (6 + s - 2 * x))) == 0,
          "F_n(s) = (2+s)(5+s)(6+s)/((6+s-2t)(6+s-2 p r))")
    # universal ceiling: t <= 1 and p r <= 1 give F <= M(s) = (2+s)(5+s)(6+s)/(4+s)^2 <= (3/20)(5+s)^2
    poly = sp.expand(3 * (5 + s) * (4 + s) ** 2 - 20 * (6 + s) * (2 + s))
    check(poly == 3 * s**3 + 19 * s**2 + 8 * s, "3(5+s)(4+s)^2 - 20(6+s)(2+s) = 3s^3+19s^2+8s >= 0, so M(s) <= (3/20)(5+s)^2")
    # midpoint formula: M^-1 = (5+s)^-1 (I + J/(2+s)); p' value minus p'' value for an antisymmetric column
    Mi = ((sp.eye(3) + sp.ones(3, 3) / (2 + s)) / (5 + s))
    # midpoints m0=(c,b'), m1=(c,e'), m2=(b',e'); old values (Wc, Wb, We)
    Wc, Wb, We = sp.symbols("Wc Wb We")
    Bw = sp.Matrix([Wc + Wb, Wc + We, Wb + We])
    mid = Mi * Bw
    anti = sp.simplify(mid[0].subs({Wc: 0, We: -Wb}))
    check(sp.simplify(anti - Wb / (5 + s)) == 0,
          "antisymmetric column e_c'-e_c'' (W_c=0, W_e'=-W_b'): W_p' = W_b'/(5+s), so y_{n+1} = y_n(ps)/((5+s)(c-er'))")
    # base cases on the triangle G_0: (L_0+s)^-1 = (3I-J+sI)^-1
    T = (3 * sp.eye(3) - sp.ones(3, 3) + s * sp.eye(3)).inv()
    check(sp.simplify(T[0, 1] - 1 / (s * (3 + s))) == 0 and sp.simplify(T[0, 0] - T[0, 1] - 1 / (3 + s)) == 0,
          "level 0: g_0 = 1/(s(3+s)), r_0 = 1/(3+s), y_0 = r_0, t_0 = s(g0+2g1) = 1")
    # cell-gradient geometry: sum_c x_c (x_c' - x_c'')(P_c' - P_c'') = -3 Q(x) on the zero-sum plane
    a_, b_ = sp.symbols("a b", real=True)
    xv = [a_, b_, -a_ - b_]
    Pc = [sp.Matrix([0, 0]), sp.Matrix([1, 0]), sp.Matrix([sp.Rational(1, 2), sp.sqrt(3) / 2])]
    Pb = (Pc[0] + Pc[1] + Pc[2]) / 3
    Q = sum((xv[k] ** 2 * (Pc[k] - Pb) for k in range(3)), sp.zeros(2, 1))
    S = sum((xv[k] * (xv[(k + 1) % 3] - xv[(k + 2) % 3]) * (Pc[(k + 1) % 3] - Pc[(k + 2) % 3]) for k in range(3)), sp.zeros(2, 1))
    check(sp.simplify(S + 3 * Q) == sp.zeros(2, 1), "sum_c x_c(x_c'-x_c'')(P_c'-P_c'') = -3 Q(x) (cyclic c -> c' -> c'')")
    # alpha = 1 consistency: R_n(1) = (2/3)(5/3)^n, sigma = 1, Y_n(1) = 5^-n R_n/2 -> BA.4
    n = sp.Symbol("n", positive=True, integer=True)
    R1 = sp.Rational(2, 3) * sp.Rational(5, 3) ** n
    Y1 = sp.Rational(1, 5) ** n * R1 / 2
    kap = 2**n * 4 / R1**2 * (3 * Y1 - 1)
    check(sp.simplify(sp.powsimp(sp.expand(kap + 9 * (sp.Rational(18, 25) ** n - sp.Rational(6, 25) ** n)), force=True)) == 0
          or all(sp.simplify((kap + 9 * (sp.Rational(18, 25) ** n - sp.Rational(6, 25) ** n)).subs(n, k)) == 0 for k in range(1, 15)),
          "alpha=1: 2^n (4/R^2)(3Y - sigma) with R=(2/3)(5/3)^n, Y=5^-n R/2, sigma=1 equals BA.4's -9((18/25)^n-(6/25)^n)")


# ----------------------------------------------------------------- (B) exact rational checks
def exact_resolvent_cols(adj, s: Fraction, cols):
    N = len(adj)
    M = [[Fraction(0)] * N for _ in range(N)]
    for i in range(N):
        M[i][i] = Fraction(len(adj[i])) + s
        for j in adj[i]:
            M[i][j] -= 1
    out = []
    for col in cols:
        A = [row[:] + [Fraction(1 if i in col else 0)] for i, row in enumerate(M)]
        for k in range(N):
            piv = next(i for i in range(k, N) if A[i][k] != 0)
            A[k], A[piv] = A[piv], A[k]
            for i in range(N):
                if i != k and A[i][k] != 0:
                    f = A[i][k] / A[k][k]
                    A[i] = [x - f * y for x, y in zip(A[i], A[k])]
        out.append([A[i][N] / A[i][i] for i in range(N)])
    return out


def exact_checks():
    ok = True
    for n in (1, 2, 3):
        P, adj, C, pb, pc = corner_frame(n)
        a, b, c = C
        for s in (Fraction(1, 2), Fraction(3), Fraction(7, 5)):
            Wa, Wb, Wc, W1 = exact_resolvent_cols(adj, s, [{a}, {b}, {c}, set(C)])
            r = Wa[a] - Wa[b]
            y = Wb[pb] - Wc[pb]
            g = Wa[b]
            e = Wb[pb] + Wc[pb] - Wb[a] - Wc[a]
            t = s * W1[a]
            rr, tt, gg, yy = rec_all(n, s)
            ok = ok and r == rr and y == yy and g == gg and t == tt and e == s * g
    check(ok, "exact rationals, levels 1-3, s in {1/2, 3, 7/5}: r_n, y_n = r_n prod 1/(5+p^j), g_n (F-product), t_n (psi), e_n = s g_n")


# ----------------------------------------------------------------- (C) float identities on the actual meshes
def dense_level(n, al=0.45):
    P, adj, C, pb, pc = corner_frame(n)
    a, b, c = C
    L = laplacian(adj)
    w, V = np.linalg.eigh(L)
    keep = np.arange(len(w)) > 0
    w[0] = 1.0
    G = (V * np.where(keep, w ** -al, 0.0)) @ V.T
    Lb = (V * np.where(keep, w ** (1 - al), 0.0)) @ V.T
    R = G[a, a] + G[b, b] - 2 * G[a, b]
    sig = Lb[a, a] - Lb[a, b]
    Y = G[pb, b] - G[pb, c]
    E = G[pb, b] + G[pb, c] - G[a, b] - G[a, c]
    # self part directly: sum_c -q_c^2 g_c[G e_c] (x component) / Q_x for d = (1,-1/2,-1/2)
    A = (V * np.where(keep, w ** al, 0.0)) @ V.T
    I = [i for i in range(len(adj)) if i not in C]
    d = np.array([1, -0.5, -0.5])
    u = np.zeros(len(adj))
    u[C] = d
    u[I] = np.linalg.solve(A[np.ix_(I, I)], -A[np.ix_(I, C)] @ d)
    q = (A @ u)[C]
    Dself = np.zeros(2)
    for k, cc in enumerate(C):
        Em = np.array([P[j] - P[cc] for j in adj[cc]])
        col = G[:, cc]
        Dself -= q[k] ** 2 * np.linalg.solve(Em, np.array([col[j] - col[cc] for j in adj[cc]]))
    Wq = P[C] - P[C].mean(0)
    Qv = ((d - d.mean())[:, None] ** 2 * Wq).sum(0)
    kself = float(Dself @ Qv / (Qv @ Qv))
    return dict(R=R, sig=sig, Y=Y, E=E, kself=kself, q0=q[0])


def quad_level(n, al: str, dps=30):
    """R_n, sigma_n, Y_n, E_n by mpmath quadrature of the recursions (consistency, not interval)."""
    mp.dps = dps
    a = mpf(al)
    pref = mp.sin(mp.pi * a) / mp.pi
    pts = [0] + [mpf(10) ** k for k in range(-30, 31, 2)]
    X = pts[-1]
    Jr = mp.quad(lambda s: s ** (-a) * rec_all(n, s)[0], pts) + X ** (-a) / a
    Jw = mp.quad(lambda s: s ** (-a) * (1 - s * rec_all(n, s)[0]), pts) + 2 * X ** (-a) / a
    Jy = mp.quad(lambda s: s ** (-a) * rec_all(n, s)[3], pts)
    Je = mp.quad(lambda s: s ** (1 - a) * rec_all(n, s)[2], pts)
    return dict(R=2 * pref * Jr, sig=pref * Jw, Y=pref * Jy, E=pref * Je)


def kap(d):
    R, sig, Y, E = d["R"], d["sig"], d["Y"], d["E"]
    return 4 / R**2 * (3 * Y - sig), -4 / R**2 * (sig - E), 4 / R**2 * (3 * Y - E)  # 2^-n kappa, 2^-n self, 2^-n cross


def ebound(n, al):
    return math.sin(math.pi * al) / math.pi * (3 / 20) ** (n * (1 - al)) * (1 / (3 * (1 - al)) + 1 / al)


# ----------------------------------------------------------------- main
def main(cert_levels=(1, 2, 3, 4), dense_levels=6) -> int:
    t0 = time.time()
    print("Theorem BB: decimation of the mesh residual kappa_n(alpha) of Theorem BA; Conjecture BA.2. Not thrust.")
    symbolic_checks()
    exact_checks()

    # (C) dense float identities, alpha = 0.45 and two other exponents
    for al in (0.45, 0.25, 0.6):
        for n in range(1, (dense_levels if al == 0.45 else 4) + 1):
            dd = dense_level(n, al)
            k_ba = ba.float_kappa(n, al)[0]
            k, ks, kc = kap(dd)
            ok = abs(2**n * k - k_ba) < 1e-9 * max(1, abs(k_ba)) and abs(2**n * ks - dd["kself"]) < 1e-9 * abs(k_ba)
            ok = ok and abs(dd["q0"] - 2 / dd["R"]) < 1e-10
            ok = ok and 0 < dd["Y"] <= 5.0**-n * dd["R"] / 2 and 0 < dd["E"] <= ebound(n, al)
            if al == 0.45:
                print(f"  level {n}: 2^-n kappa = {k:.12f}, 2^-n kappa_self = {ks:.12f}, kappa_cross = {2**n*kc:.3e}, "
                      f"Y = {dd['Y']:.3e} (<= {5.0**-n*dd['R']/2:.3e}), E = {dd['E']:.3e} (<= {ebound(n, al):.3e})")
            check(ok, f"alpha={al} level {n}: kappa_n = 2^n(4/R^2)(3Y-sigma) and self = -2^n(4/R^2)(sigma-E) match BA's dense values; "
                      f"q_c0 = 2/R; 0 < Y <= 5^-n R/2; 0 < E <= (3/20)^(n(1-a)) bound")

    # recursion quadrature reproduces dense values and continues to high levels
    print("  recursion (mpmath quadrature, consistency only):")
    prev = None
    rows = {}
    for n in (1, 4, 6, 10, 20, 40):
        qd = quad_level(n, "0.45")
        k, ks, kc = kap(qd)
        rows[n] = (k, ks, kc)
        print(f"    n={n:2d}: 2^-n kappa_n = {mp.nstr(k, 15)}, 2^-n kappa_self = {mp.nstr(ks, 15)}, kappa_cross = {mp.nstr(2**n * kc, 6)}")
        if n in (1, 4, 6):
            dd = dense_level(n)
            check(abs(float(k) - kap(dd)[0]) < 1e-9, f"recursion quadrature reproduces the dense mesh at level {n}")
    mp.dps = 50

    # (D) interval: cross part certified > 0 at levels 1..4 (Kept Failure BB.1), identity inside enclosures
    cert = {}
    for n in cert_levels:
        r = ba.certify_gasket(n)
        kc = r["kappa"] - r["kappa_self"]
        cert[n] = r
        print(f"  level {n}: certified kappa_cross = kappa - kappa_self in {ba.s(kc, 20)}")
        check(kc.a > 0, f"level {n}: kappa_cross(0.45) > 0 certified (self-force removal does not give 0 at a finite level)")
        dd = dense_level(n)
        check(abs(float(kc.mid) - 2**n * kap(dd)[2]) < 1e-9, f"level {n}: certified cross part = 2^n(4/R^2)(3Y-E)")

    # (E) the limit, certified from Theorem AS's sandwich (re-run live at alpha = 0.45)
    iv.dps = 50
    mp.dps = 50
    B, okB, _ = asd.power_ceiling()
    check(okB, "AQ (4) power ceiling re-certified (r_n <= B s^-0.35 on (0, 0.01])")
    qn, s0, vals, N = asd.build_nodes()
    c = asd.certify_limit("0.45", qn, s0, vals, N, B)
    Rinf, Sinf = c["R"], c["S"]
    K = -4 * Sinf / Rinf**2
    print(f"  R_inf(0.45) in {asd.fmt(Rinf, 15)}, sigma_inf(0.45) in {asd.fmt(Sinf, 15)}")
    print(f"  K(0.45) = lim 2^-n kappa_n = -4 sigma_inf/R_inf^2 in {asd.fmt(K, 15)}")
    check(K.b < 0 and float(K.b - K.a) < 1e-11, "Conjecture BA.2 PROVED at 0.45: lim 2^-n kappa_n = K < 0, certified width < 1e-11")
    check(abs(float(K.mid) - float(rows[40][0])) < 1e-6 and abs(float(rows[40][0]) - float(rows[40][1])) < 1e-20,
          "recursion at n=40: 2^-n kappa_n within 1e-6 (~rho^40) of K; self part equals it to 1e-20")
    check(all(float(rows[n][0]) < float(K.a) for n in rows), "2^-n kappa_n < K on every computed level (approach from below)")
    # derived rate: 2^-n kappa_n - K ~ -(8 sigma c_R/(R^3 (1-rho))) rho^n, c_R(0.45) from Theorem AT
    rho = 5**0.45 / 3
    cR = 0.49780265076356035
    ck = 8 * float(Sinf.mid) * cR / (float(Rinf.mid) ** 3 * (1 - rho))
    obs = float((rows[40][0] - K.mid) / rho**40)
    print(f"  rate: (2^-40 kappa_40 - K)/rho^40 = {obs:.6f}, predicted -8 sigma c_R/(R^3(1-rho)) = {-ck:.6f}")
    check(abs(obs + ck) < 1e-3 * ck, "derived rate constant matches the recursion at n=40 (AR (4a) + AT's c_R)")
    check(float(rows[40][2]) * 2**40 < 1e-12 and all(float(rows[n][2]) > 0 for n in rows),
          "kappa_cross > 0 and -> 0 (2^40 cross part < 1e-12)")
    # alpha >= gamma: R_n -> inf so 2^-n kappa -> 0 (recursion witness at alpha = 0.75)
    q75 = [kap(quad_level(n, "0.75"))[0] for n in (10, 20, 30)]
    print(f"  alpha=0.75: 2^-n kappa_n at n=10,20,30 = {[mp.nstr(x, 6) for x in q75]}")
    check(abs(q75[2]) < abs(q75[1]) < abs(q75[0]) and abs(q75[2]) < 0.03,
          "alpha=0.75 > log3/log5: 2^-n kappa_n -> 0 (BA.2 is specific to alpha < log3/log5)")
    check(ALPHA_STAR > 0.45 and ALPHA_STAR < GAMMA, f"alpha_* = 1 - log2/log(20/3) = {ALPHA_STAR:.6f} covers 0.45")
    mp.dps = 15
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BB checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    lv = tuple(range(1, int(sys.argv[1]) + 1)) if len(sys.argv) > 1 else (1, 2, 3, 4)
    sys.exit(main(cert_levels=lv))
