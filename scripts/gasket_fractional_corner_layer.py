#!/usr/bin/env python3
"""Theorem BH: the corner boundary layer of Kept Failure BG.1 has a limit profile.

Setting (GASKET.md, Theorem BH): level-n gasket V_n with corner a at the origin (skew lattice coordinates),
corner data d on C = {a, b, c}, x = d - dbar, u_n the L_n^alpha-Dirichlet extension, r_n, p(s) = s(s+5), D_n,
h_m as in Theorems AQ and BG.  Layer coordinates: xi in the infinite one-sided gasket graph Gamma_a, the vertex
2^-n xi of V_n (for n large).  Level-1 s-harmonic weights (sympy):  from m_ab: a, b -> (s+4)/Dd, c -> 2/Dd,
Dd = (s+2)(s+5);  alpha1(s) = (s+4)/Dd,  beta(s) = (s+6)/Dd <= 3/5.

Theorem BH (checked here; exit 1 on any failure):
  (1) corner-cell decomposition, exact, every n >= k, y in V_k:  [(L_n+s)^-1 x_C](2^-(n-k) y) is the s-harmonic
      extension on V_k of the three cell-corner values (r_n x_a, r_n h_{n-k}(p^k s; nu_b), r_n h_{n-k}(p^k s; nu_c)).
  (2) corner-neighbour profile: with nu_b, nu_c the two neighbours of a in V_N,
        h_N(z; nu_b) - h_N(z; nu_c) = (x_b - x_c) prod_{i<N} 1/(5 + p^i z)          (exact),
        |h_N(z; nu) - x_a eta(z)| <= 2 ||x||_inf prod_{i<N} beta(p^i z) <= 2 ||x||_inf (3/5)^N,
      eta(z) = sum_j alpha1(p^j z) prod_{i<j} beta(p^i z), continuous, strictly decreasing, eta(0) = 1,
      0 < eta(z) <= 2/(2+z); and r_inf(s) = lim r_n(s) = 1/(2 + s - 2 eta(s)) (corner equation).
      Corollary: u_n(nu_b) - u_n(nu_c) = (x_b - x_c) I_n/D_n exactly, I_n = int s^-a y_n ds (BB's y_n), 0 < I_n/D_n <= 5^-n.
  (3) layer limit for alpha < gamma: u_n(2^-(n-k) y) -> dbar + x_a Theta_k(y),
        Theta_k(y) = int s^-a r_inf(s) phi_k(s; y) ds / J_r,  phi_k = H_k(s; a, y) + eta(p^k s)(H_k(s; b, y) + H_k(s; c, y)),
      with |error| <= 2||x||_inf (3/5)^(n-k) + 2|x_a| (J_r - D_n)/D_n; universal (x_b, x_c drop out).
  (4) no layer for alpha >= gamma: u_n(2^-n xi) -> d_a for every fixed xi.
  (5) decay into the bulk: theta_k = Theta at xi = 2^k e_b, = int s^-a r_inf eta(p^k s)/J_r, strictly decreasing to 0;
      theta_0 = 1 - J_w/(2 J_r) = 1 - 3L/sqrt21, so x_a theta_0 is BG.1's gap, now at each neighbour separately.
      Certified at 0.45 (python-flint arb balls on monotone cells, J_r and J_w from Theorem AS re-run live):
      theta_0, theta_1 (xi = (2,0)), theta_2 (xi = (4,0)), theta_Y (xi = (1,1)).
  (6) far field (DERIVED from Theorem AT's small-y lemma and BD.4): rho^-k theta_k -> theta_inf > 0 (floats ~0.52061).

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
from flint import arb, ctx
from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_maxwell_limit as bb  # noqa: E402
import mesh_maxwell_residual as ba  # noqa: E402
import gasket_fractional_limit_digits as asd  # noqa: E402

FAILS: list[str] = []
GAMMA = math.log(3) / math.log(5)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    z, s, t, r, R, e0 = sp.symbols("z s t r R eta0", positive=True)
    A, B, C, X, Y, Z = sp.symbols("A B C X Y Z")
    sol = sp.solve([(4 + z) * X - (A + B + Y + Z), (4 + z) * Y - (B + C + X + Z), (4 + z) * Z - (C + A + X + Y)], [X, Y, Z])
    Dd = (z + 2) * (z + 5)
    al1, be = (z + 4) / Dd, (z + 6) / Dd
    w = [sp.simplify(sp.diff(sol[X], v)) for v in (A, B, C)]
    check(sp.simplify(w[0] - al1) == 0 and sp.simplify(w[1] - al1) == 0 and sp.simplify(w[2] - 2 / Dd) == 0,
          "level-1 s-harmonic measure from m_ab: weights (s+4)/((s+2)(s+5)) on a, b and 2/((s+2)(s+5)) on c")
    check(sp.simplify(al1 + be - 2 / (z + 2)) == 0 and sp.simplify(al1 - 2 / Dd - 1 / (z + 5)) == 0,
          "alpha1 + beta = 2/(2+z) (so eta <= 1 => eta <= 2/(2+z)); antisymmetric factor alpha1 - 2/Dd = 1/(z+5)")
    dbe = sp.factor(sp.together(sp.diff(be, z)))
    check(sp.simplify(be.subs(z, 0) - sp.Rational(3, 5)) == 0 and sp.simplify(dbe + (z**2 + 12 * z + 32) / Dd**2) == 0,
          "beta(0) = 3/5 and beta' = -(z^2+12z+32)/Dd^2 < 0, so beta <= 3/5 on [0, inf) (the contraction factor)")
    j = sp.Symbol("j", integer=True, nonnegative=True)
    check(sp.summation(sp.Rational(2, 5) * sp.Rational(3, 5) ** j, (j, 0, sp.oo)) == 1, "eta(0) = sum (2/5)(3/5)^j = 1")
    # corner equation r(2+s-2 eta)=1 and the fixed point eta(s) = alpha1(s) + beta(s) eta(p s) reproduce AQ's recursion
    ps = s * (s + 5)
    eta_s = lambda rr, ss: (2 + ss - 1 / rr) / 2  # noqa: E731
    lhs = 1 / (2 + s - 2 * (al1.subs(z, s) + be.subs(z, s) * eta_s(R, ps)))
    cq = (6 + s) / ((5 + s) * (2 + s))
    eq = 2 * s / (2 + s)
    check(sp.simplify(lhs - R / (cq - eq * R)) == 0,
          "r = 1/(2+s-2eta) with eta(s) = alpha1 + beta eta(ps) is AQ's fixed point r(s) = r(ps)/(c - e r(ps)), c = beta")
    check(sp.simplify(be.subs(z, s) - eq / ps - (s + 4) / ((s + 2) * (s + 5))) == 0,
          "c - e/p = (s+4)/((s+2)(s+5)): with r(ps) <= 1/(ps), r(s)/r(ps) <= (s+2)(s+5)/(s+4) (head lemma)")
    phiY = 2 * (t + (al1 + be * t).subs(z, s)) / (4 + s)
    check(sp.simplify(phiY - (2 + 2 * (s + 4) * t) / ((s + 2) * (s + 5))) == 0
          and sp.simplify(1 - phiY - (s * (s + 5) + 2 * (s + 4) * (1 - t)) / ((s + 2) * (s + 5))) == 0,
          "cell-midpoint xi=(1,1): phi_Y = (2+2(s+4)eta(ps))/((s+2)(s+5)), 1-phi_Y = (s(s+5)+2(s+4)(1-eta(ps)))/((s+2)(s+5))")
    check(sp.simplify(r * eta_s(r, s) - (r * (2 + s) - 1) / 2) == 0 and sp.simplify(1 - eta_s(r, s) - (1 - s * r) / (2 * r)) == 0,
          "r eta = r - (1-sr)/2, 1 - eta = (1-sr)/(2r): theta_0 = 1 - J_w/(2J_r)")
    Ls = sp.Symbol("L", positive=True)
    Jw = 6 * Ls * R / sp.sqrt(21)  # AS: L = sqrt21 J_w/(6 J_r)
    th0 = 1 - Jw / (2 * R)
    check(sp.simplify(sp.Rational(5, 6) * th0 - (sp.Rational(5, 6) - 5 * Ls / (2 * sp.sqrt(21)))) == 0,
          "with AS's L = sqrt21 J_w/(6J_r): x_a theta_0 = 5/6 - 5L/(2 sqrt21) for (1,-1/2,0), i.e. Kept Failure BG.1's gap")


# ----------------------------------------------------------------- arb helpers
def eta_arb(z: arb) -> arb:
    """eta(z) = sum_j alpha1(p^j z) prod_{i<j} beta(p^i z), with tail in [0, P * 2/(2+p^J z)] (eta <= 2/(2+.))."""
    S, P = arb(0), arb(1)
    while True:
        D = (z + 2) * (z + 5)
        S += P * (z + 4) / D
        P *= (z + 6) / D
        z = z * (z + 5)
        if P.upper() < arb("1e-30"):
            break
    lo, hi = S.lower(), (S + P * 2 / (2 + z)).upper()
    return (lo + hi) / 2 + arb(0, (hi - lo) / 2 * (1 + arb("1e-20")))


def beta_prod(z: arb, N: int) -> arb:
    P = arb(1)
    for _ in range(N):
        P *= (z + 6) / ((z + 2) * (z + 5))
        z = z * (z + 5)
    return P


def F(q: Fraction) -> arb:
    return arb(q.numerator) / arb(q.denominator)


# ----------------------------------------------------------------- (B) exact rationals
def exact_solve(adj, s: Fraction, rhs, fixed=None):
    """(L+s) f = rhs on free vertices; fixed = {i: value} imposes Dirichlet values."""
    N = len(adj)
    fixed = fixed or {}
    A = [[Fraction(0)] * N + [Fraction(0)] for _ in range(N)]
    for i in range(N):
        if i in fixed:
            A[i][i] = Fraction(1)
            A[i][N] = fixed[i]
            continue
        A[i][i] = Fraction(len(adj[i])) + s
        for jj in adj[i]:
            A[i][jj] -= 1
        A[i][N] = rhs[i]
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


def exact_checks(levels=(2, 3, 4)):
    ctx.prec = 100
    xC = (Fraction(3), Fraction(-1), Fraction(-2))  # zero-sum, ||x|| = 3, x_b - x_c = 1
    xn = 3
    ok_cell = ok_asym = ok_eta = ok_r = True
    for n in levels:
        adj, C, pos = mesh(n)
        for s in (Fraction(1, 2), Fraction(7, 5)):
            rhs = [Fraction(0)] * len(adj)
            for kk, cc in enumerate(C):
                rhs[cc] = xC[kk]
            w = exact_solve(adj, s, rhs)
            rn = bb.rec_all(n, s)[0]
            # (1) corner cells of every size k < n
            for k in range(1, n):
                adjk, Ck, posk = mesh(k)
                sc = Fraction(1, 2 ** (n - k))
                fixed = {Ck[0]: w[pos[(Fraction(0), Fraction(0))]], Ck[1]: w[pos[(sc, Fraction(0))]],
                         Ck[2]: w[pos[(Fraction(0), sc)]]}
                f = exact_solve(adjk, s, [Fraction(0)] * len(adjk), fixed)
                ok_cell = ok_cell and all(w[pos[(key[0] * sc, key[1] * sc)]] == f[i] for key, i in posk.items())
            # (2) neighbours of a at level n: asymmetry identity and contraction to x_a eta
            hb = w[pos[(Fraction(1, 2**n), Fraction(0))]] / rn
            hc = w[pos[(Fraction(0), Fraction(1, 2**n))]] / rn
            prod, zz = Fraction(1), s
            for _ in range(n):
                prod /= 5 + zz
                zz = zz * (zz + 5)
            ok_asym = ok_asym and hb - hc == (xC[1] - xC[2]) * prod
            et = eta_arb(F(s))
            bound = 2 * xn * beta_prod(F(s), n)
            for hv in (hb, hc):
                ok_eta = ok_eta and bool(abs(F(hv) - F(xC[0]) * et) < bound)
            # corner equation: r_n = x_a/((2+s)x_a - h_b - h_c)
            ok_r = ok_r and rn == xC[0] / ((2 + s) * xC[0] - hb - hc)
    check(ok_cell, f"exact, levels {levels[0]}-{levels[-1]}, s in {{1/2, 7/5}}: on every corner cell 2^-(n-k) V_k (k<n) the corner-data "
                   "resolvent is the s-harmonic extension of its three cell-corner values (BH.1)")
    check(ok_asym, "exact: h_n(s; nu_b) - h_n(s; nu_c) = (x_b - x_c) prod_{i<n} 1/(5+p^i s) (BH.2 asymmetry)")
    check(ok_eta, "arb: |h_n(s; nu) - x_a eta(s)| < 2||x|| prod_{i<n} beta(p^i s) at both neighbours (BH.2 contraction)")
    check(ok_r, "exact: corner equation r_n(s) = x_a/((2+s)x_a - h_n(s; nu_b) - h_n(s; nu_c))")
    ok = True
    for sv in ("0.001", "1", "10"):
        rinf = 1 / (2 + arb(sv) - 2 * eta_arb(arb(sv)))
        r60 = bb.rec_all(60, mpf(sv))[0]
        ok = ok and abs(float(rinf.mid()) - float(r60)) < 1e-12 * float(r60)
    check(ok, "r_inf = 1/(2+s-2eta(s)) agrees with AQ's level-60 recursion at s = 0.001, 1, 10")


# ----------------------------------------------------------------- (C) dense floats
def dense_u(n, al, d=(1.0, -0.5, 0.0)):
    P, adj, C, lat = ba.gasket_setup(n)
    L = bb.laplacian(adj)
    w, V = np.linalg.eigh(L)
    keep = np.arange(len(w)) > 0
    w[0] = 1.0
    A = (V * np.where(keep, w**al, 0.0)) @ V.T
    I = [i for i in range(len(adj)) if i not in C]
    u = np.zeros(len(adj))
    u[C] = d
    u[I] = np.linalg.solve(A[np.ix_(I, I)], -A[np.ix_(I, C)] @ np.array(d))
    return u, {lat[i]: i for i in range(len(adj))}


def DY(n, al, dps=25):
    mp.dps = dps
    a = mpf(al)
    pts = [0] + [mpf(10) ** k for k in range(-30, 31, 2)]
    X = pts[-1]
    D = mp.quad(lambda s: s ** (-a) * bb.rec_all(n, s)[0], pts) + X ** (-a) / a
    Y = mp.quad(lambda s: s ** (-a) * bb.rec_all(n, s)[3], pts)
    return D, Y


XI = {"nu_b (1,0)": (1, 0, 0), "nu_c (0,1)": (0, 1, 0), "(2,0)": (2, 0, 1), "(1,1)": (1, 1, 1), "(4,0)": (4, 0, 2)}


def dense_checks(levels=(1, 2, 3, 4, 5, 6), quad_levels=(1, 2, 3, 4), Jr=4.0013017734819):
    xa, xb, xc, dbar = 5 / 6, -2 / 3, -1 / 6, 1 / 6
    out = {}
    for al in (0.45, 0.8):
        for n in levels:
            u, pos = dense_u(n, al)
            vals = {}
            for name, (i, j, k) in XI.items():
                key = (Fraction(i, 2**n), Fraction(j, 2**n))
                if key in pos and i + j <= 2**n:
                    vals[name] = (u[pos[key]] - dbar) / xa
            out[(al, n)] = (u, pos, vals)
    ok = True
    for n in quad_levels:
        D, Y = DY(n, "0.45")
        u, pos, _ = out[(0.45, n)]
        diff = u[pos[(Fraction(1, 2**n), Fraction(0))]] - u[pos[(Fraction(0), Fraction(1, 2**n))]]
        ok = ok and abs(diff - (xb - xc) * float(Y / D)) < 1e-10 and 0 < float(Y / D) <= 5.0**-n
        out[("D", n)] = float(D)
    check(ok, f"alpha=0.45, dense levels {quad_levels}: u_n(nu_b) - u_n(nu_c) = (x_b-x_c) I_n/D_n to 1e-10, 0 < I_n/D_n <= 5^-n")
    return out


# ----------------------------------------------------------------- (D) certificate
def grid(rf=1.00002, rm=1.0005, rc=1.01):
    pts, s = [], 1e-16
    while s < 1e24:
        pts.append(s)
        if 1e-3 <= s < 1e3:
            s *= rf
        elif 1e-8 <= s < 1e8:
            s *= rm
        else:
            s *= rc
    pts.append(1e24)
    return pts


def M_bracket(al: str, pts=None):
    """Enclosures of M_phi = int_0^inf s^-a r_inf(s)(1 - phi(s)) ds for phi in {eta, eta o p, eta o p^2, phi_Y}.
    r_inf decreasing, 1-phi increasing (eta decreasing): cell value in [r(s1)(1-phi(s0)), r(s0)(1-phi(s1))] int s^-a.
    Head (0, s_lo]: r(1-phi) <= 4 (BH head lemma); tail [s_hi, inf): r(1-phi) <= r <= 1/s."""
    ctx.prec = 80
    a = arb(al)
    pts = pts or grid()
    vals = []
    for sf in pts:
        s = arb(sf)
        e0 = eta_arb(s)
        ps = s * (s + 5)
        e1 = eta_arb(ps)
        e2 = eta_arb(ps * (ps + 5))
        r = 1 / (2 + s - 2 * e0)
        cY = (s * (s + 5) + 2 * (s + 4) * (1 - e1)) / ((s + 2) * (s + 5))
        vals.append((s ** (1 - a), r, (1 - e0, 1 - e1, 1 - e2, cY)))
    lo = [arb(0)] * 4
    hi = [arb(0)] * 4
    for (w0, r0, c0), (w1, r1, c1) in zip(vals, vals[1:]):
        wt = (w1 - w0) / (1 - a)
        for j in range(4):
            lo[j] += r1 * c0[j] * wt
            hi[j] += r0 * c1[j] * wt
    head = 4 * arb(pts[0]) ** (1 - a) / (1 - a)
    tail = arb(pts[-1]) ** (-a) / a
    out = []
    for j in range(4):
        L_, H_ = lo[j].lower(), (hi[j] + head + tail).upper()
        # float() rounds to nearest (relative error < 2^-52); widen by 1e-15 relative (all values positive)
        out.append(iv.mpf([float(L_) * (1 - 1e-15), float(H_) * (1 + 1e-15)]))
    return out, len(pts)


def float_far_field(ks=(0, 2, 4, 6, 8), al="0.45", Jr="4.0013017734819922628"):
    mp.dps = 20
    a = mpf(al)
    rho = mpf(5) ** a / 3

    def eta(z):
        tot, pr = mpf(0), mpf(1)
        while pr > mpf(10) ** -22:
            D = (z + 2) * (z + 5)
            tot += pr * (z + 4) / D
            pr *= (z + 6) / D
            z = z * (z + 5)
        return tot

    def pk(s, k):
        for _ in range(k):
            s = s * (s + 5)
        return s

    res = []
    for k in ks:
        pts = [0] + [mpf(5) ** (-k) * mpf(10) ** j for j in range(-24, 8)] + [mp.inf]
        th = mp.quad(lambda s: s ** -a / (2 + s - 2 * eta(s)) * eta(pk(s, k)), pts) / mpf(Jr)
        res.append((k, th, th / rho**k))
    return res


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem BH: limit profile of the corner boundary layer of the fractional Dirichlet extension. Not thrust.")
    symbolic_checks()
    exact_checks()
    dd = dense_checks(levels=(1, 2, 3, 4, 5, 6, 7))

    t1 = time.time()
    M, npts = M_bracket("0.45")
    print(f"   M_phi enclosures ({npts} nodes, {time.time() - t1:.1f}s): " + ", ".join(f"[{float(m.a):.9f}, {float(m.b):.9f}]" for m in M))
    iv.dps = 50
    mp.dps = 50
    B, okB, _ = asd.power_ceiling()
    check(okB, "Theorem AS power ceiling re-certified")
    qn, s0, vals, N = asd.build_nodes()
    cAS = asd.certify_limit("0.45", qn, s0, vals, N, B)
    iv.dps = 20
    Jr = iv.mpf([cAS["Jr"].a, cAS["Jr"].b])
    Jw = iv.mpf([cAS["Jw"].a, cAS["Jw"].b])
    Lim = iv.mpf([cAS["L"].a, cAS["L"].b])
    print(f"   AS re-run: J_r in {asd.fmt(Jr, 14)}, J_w in {asd.fmt(Jw, 14)}")
    half = Jw / 2
    check(M[0].a <= half.b and half.a <= M[0].b,
          "cross-check: the eta-series certificate of M_eta = int s^-a r_inf (1-eta) contains AS's J_w/2 (independent routes)")
    th0 = 1 - Jw / (2 * Jr)
    th = [1 - m / Jr for m in M]
    names = ["theta_0 (xi=(1,0))", "theta_1 (xi=(2,0))", "theta_2 (xi=(4,0))", "theta_Y (xi=(1,1))"]
    th[0] = iv.mpf([max(th[0].a, th0.a), min(th[0].b, th0.b)])
    for nm, x in zip(names, th):
        print(f"   {nm} in [{mp.nstr(mpf(x.a), 13)}, {mp.nstr(mpf(x.b), 13)}]")
    quoted = [("0.4624385835869", "0.462438583587"), ("0.340381375", "0.340407195"), ("0.241594653", "0.241622693"),
              ("0.382431398", "0.382456436")]
    check(all(mpf(lo) <= x.a and x.b <= mpf(hi) for x, (lo, hi) in zip(th, quoted)),
          "the enclosures quoted in GASKET.md (BH.5) contain the certified ones (outward rounding)")
    gap = iv.mpf(5) / 6 - 5 * Lim / (2 * iv.sqrt(iv.mpf(21)))
    check(abs(float((iv.mpf(5) / 6 * th0).a) - float(gap.a)) < 1e-12,
          f"x_a theta_0 = BG.1 gap: (5/6) theta_0 in [{float((iv.mpf(5)/6*th0).a):.10f}, {float((iv.mpf(5)/6*th0).b):.10f}]")
    check(th[0].a > th[3].b and th[3].a > th[1].b and th[1].a > th[2].b and th[2].a > 0.24,
          "certified ordering theta_0 > theta_Y > theta_1 > theta_2 > 0.24: a genuine layer, decaying into the bulk")

    # dense levels vs certified limits at 0.45 (increasing, below the limit, inside the BH.3 error bound)
    key = {"nu_b (1,0)": th[0], "nu_c (0,1)": th[0], "(2,0)": th[1], "(4,0)": th[2], "(1,1)": th[3]}
    okinc = okbd = True
    for nm, lim in key.items():
        k = XI[nm][2]
        seq = [dd[(0.45, n)][2][nm] for n in range(1, 8) if nm in dd[(0.45, n)][2]]
        print(f"   dense alpha=0.45 {nm}: " + ", ".join(f"{v:.6f}" for v in seq) + f"  -> [{float(lim.a):.6f}, {float(lim.b):.6f}]")
        okinc = okinc and all(seq[i] < seq[i + 1] for i in range(len(seq) - 1)) and seq[-1] < float(lim.a)
        for n in range(1, 5):
            if nm in dd[(0.45, n)][2]:
                bnd = 2 * (5 / 6) / (5 / 6) * 0.6 ** (n - k) + 2 * (float(Jr.b) - dd[("D", n)]) / dd[("D", n)]
                okbd = okbd and abs(dd[(0.45, n)][2][nm] - float(lim.a)) <= bnd
    check(okinc, "alpha=0.45, dense levels to 7: each layer value increases with n and stays below its certified limit")
    check(okbd, "alpha=0.45, levels 1-4: |(u_n - dbar)/x_a - Theta| <= 2(||x||/x_a)(3/5)^(n-k) + 2(J_r - D_n)/D_n (BH.3 rate)")
    ok8 = True
    for nm in ("nu_b (1,0)", "(2,0)", "(1,1)"):
        seq = [dd[(0.8, n)][2][nm] for n in range(1, 8) if nm in dd[(0.8, n)][2]]
        print(f"   dense alpha=0.8 {nm}: " + ", ".join(f"{v:.6f}" for v in seq) + "  -> 1 (BH.4)")
        ok8 = ok8 and all(seq[i] < seq[i + 1] for i in range(len(seq) - 1)) and seq[-1] > float(key[nm].b) + 0.3
    check(ok8, "alpha=0.8 > gamma: normalised layer values increase with n and exceed the 0.45 layer limits by > 0.3 at level 7 (consistent with BH.4: -> 1, no layer; slow, no rate proved)")

    ff = float_far_field()
    for k, v, w in ff:
        print(f"   floats: theta_{k} = {mp.nstr(v, 12)}, rho^-k theta_{k} = {mp.nstr(w, 12)}")
    check(abs(float(ff[0][1]) - float(th[0].mid.a)) < 1e-9, "floats: theta_0 by 20-digit quadrature agrees with the certified value to 1e-9")
    dif = [float(ff[i + 1][2] - ff[i][2]) for i in range(len(ff) - 1)]
    check(all(x > 0 for x in dif) and all(dif[i + 1] < dif[i] / 5 for i in range(len(dif) - 1)) and abs(float(ff[-1][2]) - 0.52061) < 2e-4,
          "floats: rho^-k theta_k increases with increments shrinking ~9x per 2 levels, settling near 0.52061 (BH.6 evidence, not a certificate)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BH checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
