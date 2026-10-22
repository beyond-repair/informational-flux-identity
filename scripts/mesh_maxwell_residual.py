#!/usr/bin/env python3
"""Theorem BA: the Maxwell-form reading on the frozen 0.45 mesh.

Mesh: level-n combinatorial gasket, A = L^alpha (spectral calculus, 0^alpha = 0), alpha = 0.45.
Psi = u, the A-Dirichlet extension of corner data d.  Discrete charge q = A u (supported on the
three corners, Theorem H).  Discrete field E_c = -g_c, g_c = gradient of the affine interpolant of
u on the unique level-n cell at corner c.  Discrete divergence of the quadratic (Maxwell-form)
reading at vertex i: f_i = q_i E_i (continuum: d_j T^{ij} = d^i Psi Lap Psi = rho E).

Pair split: u = G q + const with the symmetric Green kernel G = A^+, so f_c = sum_k phi_ck with
phi_ck = -q_c q_k H_ck, H_ck = g_c[G e_k].  Maxwell (action-reaction) piece = antisymmetric part of
phi, residual = symmetric part (diagonal self terms included).

Checks (exit 1 on any failure):
  * exact minimal polynomial of L_n (Krylov over Z, verified m(L)=0 on every basis vector),
    rational root isolation, L^alpha and A^+ as interval interpolation polynomials (mpmath.iv);
  * Krawczyk-verified Dirichlet solve;
  * neutrality sum q = 0, sum of Maxwell piece = 0, residual total = D (u = Gq + const);
  * D = kappa_n(alpha) Q(d), Q(d) = sum_c (d_c - dbar)^2 (P_c - Pbar), |Q| = |d - dbar|^2/(2 sqrt 3);
  * certified kappa_n(0.45) < 0 for n = 1..4 (so D != 0 for every nonconstant d);
  * alpha = 1 closed form kappa_n(1) = -9 (0.72^n - 0.24^n) (symbolic in n, exact at levels 1-3);
  * torus control: centred stencil gives D = 0 exactly (odd stencil, even translation-invariant
    kernel), one-sided cell stencil gives a certified nonzero D on the same translation-invariant
    torus (Kept Failure BA.1);
  * float table kappa_n(0.45), n = 1..6 (2^-n kappa_n decreasing in modulus, Conjecture BA.2).

Not a field solver.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv, mp, nstr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_currents import build_gasket  # noqa: E402

FAILS: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- exact graph data
class Graph:
    """Integer Laplacian L = D - A_adj given by adjacency lists (multi-edges not allowed)."""

    def __init__(self, adj: list[list[int]]):
        self.adj = adj
        self.N = len(adj)
        self.deg = [len(a) for a in adj]

    def mul(self, v: list[int]) -> list[int]:
        return [self.deg[i] * v[i] - sum(v[j] for j in self.adj[i]) for i in range(self.N)]


def krylov_minpoly(g: Graph, seed: int = 7) -> list[int]:
    """Monic integer polynomial m (coeffs low->high) with m(L) v = 0 for a random v; then verified
    to annihilate every basis vector exactly, so spec(L) is contained in the roots of m."""
    rnd = random.Random(seed)
    v = [rnd.randint(-50, 50) for _ in range(g.N)]
    vecs = [v]
    basis: list[tuple[int, list[Fraction], list[Fraction]]] = []  # (pivot, reduced row, combination)
    k = 0
    while True:
        w = [Fraction(x) for x in vecs[k]]
        comb = [Fraction(0)] * (k + 1)
        comb[k] = Fraction(1)
        for piv, row, cmb in basis:
            if w[piv] != 0:
                f = w[piv] / row[piv]
                w = [a - f * b for a, b in zip(w, row)]
                for t in range(len(cmb)):
                    comb[t] -= f * cmb[t]
        piv = next((i for i, x in enumerate(w) if x != 0), None)
        if piv is None:
            # sum comb_t L^t v = 0 with comb_k = 1
            den = math.lcm(*[c.denominator for c in comb])
            coeffs = [int(c * den) for c in comb]
            assert coeffs[-1] == den
            m = [c // den for c in coeffs] if all(c % den == 0 for c in coeffs) else None
            assert m is not None, "minimal polynomial not integral"
            return m
        basis.append((piv, w, comb))
        vecs.append(g.mul(vecs[k]))
        k += 1


def krylov_columns(g: Graph, m: int) -> list[list[list[int]]]:
    """K[j][t] = L^t e_j (exact integers), t = 0..m."""
    out = []
    for j in range(g.N):
        e = [0] * g.N
        e[j] = 1
        seq = [e]
        for _ in range(m):
            seq.append(g.mul(seq[-1]))
        out.append(seq)
    return out


def verify_annihilates(K, mpoly) -> bool:
    deg = len(mpoly) - 1
    for seq in K:
        N = len(seq[0])
        for i in range(N):
            if sum(mpoly[t] * seq[t][i] for t in range(deg + 1)) != 0:
                return False
    return True


def isolate_roots(mpoly, eps=Fraction(1, 10**70)):
    x = sp.Symbol("x")
    P = sp.Poly(list(reversed(mpoly)), x)
    assert sp.gcd(P, P.diff(x)).degree() == 0, "minimal polynomial not squarefree"
    ivs = P.intervals(eps=sp.Rational(eps.numerator, eps.denominator))
    assert len(ivs) == P.degree(), "minimal polynomial has non-real roots"
    return [(Fraction(int(a.p), int(a.q)), Fraction(int(b.p), int(b.q))) for (a, b), _ in ivs]


def ivq(f: Fraction):
    return iv.mpf(f.numerator) / f.denominator


def interp_coeffs(roots, alpha, power_sign):
    """Monomial coefficients (intervals) of the interpolant of x^(power_sign*alpha) (0 at x=0)."""
    lam, vals = [], []
    for a, b in roots:
        if a == b:
            l = ivq(a)
        else:
            lo, hi = ivq(a), ivq(b)
            l = iv.mpf([lo.a, hi.b])
        lam.append(l)
        if a == b == 0:
            vals.append(iv.mpf(0))
        else:
            assert a > 0, "isolating interval must lie in (0, inf) for a nonzero eigenvalue"
            vals.append(iv.exp(power_sign * alpha * iv.log(l)))
    m = len(lam)
    dd = list(vals)
    for j in range(1, m):
        for i in range(m - 1, j - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) / (lam[i] - lam[i - j])
    coef = [iv.mpf(0)] * m
    poly = [iv.mpf(1)]
    for k in range(m):
        for i, c in enumerate(poly):
            coef[i] += dd[k] * c
        if k < m - 1:
            new = [iv.mpf(0)] * (len(poly) + 1)
            for i, c in enumerate(poly):
                new[i + 1] += c
                new[i] -= lam[k] * c
            poly = new
    return coef


def matfun(K, coef):
    N = len(K)
    A = [[iv.mpf(0)] * N for _ in range(N)]
    for j in range(N):
        seq = K[j]
        for t, c in enumerate(coef):
            col = seq[t]
            for i in range(N):
                if col[i]:
                    A[i][j] += c * col[i]
    return A


def verified_solve(A, b):
    """Enclosure of the solution of A x = b for every A, b inside the interval data
    (midpoint inverse R, ||I - R A||_inf < 1, Krawczyk/Rump a-posteriori bound)."""
    n = len(b)
    mp.prec = iv.prec
    Am = mp.matrix([[A[i][j].mid for j in range(n)] for i in range(n)])
    bm = mp.matrix([b[i].mid for i in range(n)])
    R = mp.inverse(Am)
    xt = R * bm
    Riv = [[iv.mpf(R[i, j]) for j in range(n)] for i in range(n)]
    xiv = [iv.mpf(xt[i]) for i in range(n)]
    r = [b[i] - sum((A[i][j] * xiv[j] for j in range(n)), iv.mpf(0)) for i in range(n)]
    Z = [sum((Riv[i][k] * r[k] for k in range(n)), iv.mpf(0)) for i in range(n)]
    beta = 0
    for i in range(n):
        row = iv.mpf(0)
        for j in range(n):
            s = sum((Riv[i][k] * A[k][j] for k in range(n)), iv.mpf(0))
            row += abs((iv.mpf(1) if i == j else iv.mpf(0)) - s)
        beta = max(beta, row.b)
    beta = iv.mpf(beta)
    if not beta.b < 1:
        raise RuntimeError("Krawczyk contraction failed")
    delta = iv.mpf(max(abs(z).b for z in Z)) / (1 - beta)
    box = iv.mpf([-delta.b, delta.b])
    return [xiv[i] + box for i in range(n)]


# ----------------------------------------------------------------- reading on a graph
def cell_gradient(evec, du):
    """Gradient of the affine function with increments du along the two edge vectors evec."""
    (ax, ay), (bx, by) = evec
    det = ax * by - ay * bx
    return ((du[0] * by - du[1] * ay) / det, (ax * du[1] - bx * du[0]) / det)


def maxwell_reading(K, roots, alpha, dev, data, grad):
    """dev: device vertices (Dirichlet set); grad(u, c) -> vector gradient at device vertex c.
    Returns q, f, D, Maxwell and residual totals, self term."""
    N = len(K)
    A = matfun(K, interp_coeffs(roots, alpha, +1))
    G = matfun(K, interp_coeffs(roots, alpha, -1))
    I = [i for i in range(N) if i not in dev]
    rhs = [-sum((A[i][dev[k]] * data[k] for k in range(len(dev))), iv.mpf(0)) for i in I]
    uI = verified_solve([[A[i][j] for j in I] for i in I], rhs)
    u = [iv.mpf(0)] * N
    for k, c in enumerate(dev):
        u[c] = data[k]
    for t, i in enumerate(I):
        u[i] = uI[t]
    q = [sum((A[c][j] * u[j] for j in range(N)), iv.mpf(0)) for c in dev]
    qint = max(abs(sum((A[i][j] * u[j] for j in range(N)), iv.mpf(0))).b for i in I)
    gs = [grad(u, c) for c in dev]
    f = [(-q[k] * gs[k][0], -q[k] * gs[k][1]) for k in range(len(dev))]
    D = tuple(sum((x[s] for x in f), iv.mpf(0)) for s in (0, 1))
    nd = len(dev)
    H = [[grad([G[j][dev[k]] for j in range(N)], dev[c]) for k in range(nd)] for c in range(nd)]
    phi = [[tuple(-q[c] * q[k] * H[c][k][s] for s in (0, 1)) for k in range(nd)] for c in range(nd)]
    Mtot = tuple(sum(((phi[c][k][s] - phi[k][c][s]) / 2 for c in range(nd) for k in range(nd)), iv.mpf(0)) for s in (0, 1))
    Rtot = tuple(sum(((phi[c][k][s] + phi[k][c][s]) / 2 for c in range(nd) for k in range(nd)), iv.mpf(0)) for s in (0, 1))
    self_ = tuple(sum((phi[c][c][s] for c in range(nd)), iv.mpf(0)) for s in (0, 1))
    return dict(q=q, qint=qint, f=f, D=D, M=Mtot, R=Rtot, self=self_, u=u)


def contains0(x) -> bool:
    return x.a <= 0 <= x.b


def excludes0(x) -> bool:
    return x.b < 0 or x.a > 0


def disjoint(x, y) -> bool:
    return x.b < y.a or y.b < x.a


# ----------------------------------------------------------------- the 0.45 mesh
def gasket_setup(n):
    P, E, C = build_gasket(n)
    N = len(P)
    adj = [[] for _ in range(N)]
    for i, j in E:
        adj[i].append(j)
        adj[j].append(i)
    s = 2**n
    lat = []
    for x, y in P:  # oblique coordinates: P = a(1,0) + b(1/2, sqrt3/2), a, b in 2^-n Z
        b = round(2 * y / math.sqrt(3) * s)
        a = round((x - y / math.sqrt(3)) * s)
        lat.append((Fraction(a, s), Fraction(b, s)))
    return P, adj, [int(c) for c in C], lat


def euclid(ab):
    a, b = ab
    r3 = iv.sqrt(iv.mpf(3))
    return (ivq(a) + ivq(b) / 2, ivq(b) * r3 / 2)


def certify_gasket(n, alpha=Fraction(45, 100), data=(Fraction(1), Fraction(-1, 2), Fraction(-1, 2)), dps=120):
    iv.dps = dps
    P, adj, C, lat = gasket_setup(n)
    g = Graph(adj)
    mpoly = krylov_minpoly(g)
    K = krylov_columns(g, len(mpoly) - 1)
    ann = verify_annihilates(K, mpoly)
    roots = isolate_roots(mpoly)
    pos = [euclid(t) for t in lat]
    for c in C:
        assert len(adj[c]) == 2

    def grad(u, c):
        j1, j2 = adj[c]
        ev = [(pos[j][0] - pos[c][0], pos[j][1] - pos[c][1]) for j in (j1, j2)]
        return cell_gradient(ev, [u[j1] - u[c], u[j2] - u[c]])

    al = ivq(alpha)
    dv = [ivq(x) for x in data]
    out = maxwell_reading(K, roots, al, C, dv, grad)
    dbar = sum(dv, iv.mpf(0)) / 3
    cen = [sum((pos[c][s] for c in C), iv.mpf(0)) / 3 for s in (0, 1)]
    Q = tuple(sum(((dv[k] - dbar) ** 2 * (pos[C[k]][s] - cen[s]) for k in range(3)), iv.mpf(0)) for s in (0, 1))
    out.update(Q=Q, kappa=out["D"][0] / Q[0], cross=out["D"][0] * Q[1] - out["D"][1] * Q[0],
               kappa_self=out["self"][0] / Q[0], deg=len(mpoly) - 1, ann=ann, N=g.N)
    return out


# ----------------------------------------------------------------- alpha = 1, symbolic in n
def alpha_one_symbolic():
    """Harmonic extension: corner-cell values (u_c, u_toward c', u_toward c'') = M^n (d_c, d_c', d_c'')."""
    nn = sp.Symbol("n", positive=True, integer=True)
    M = sp.Matrix([[1, 0, 0], [sp.Rational(2, 5), sp.Rational(2, 5), sp.Rational(1, 5)],
                   [sp.Rational(2, 5), sp.Rational(1, 5), sp.Rational(2, 5)]])
    Pm, J = M.diagonalize()
    Mn = Pm * sp.diag(*[J[i, i] ** nn for i in range(3)]) * Pm.inv()
    r3 = sp.sqrt(3)
    Pc = [sp.Matrix([0, 0]), sp.Matrix([1, 0]), sp.Matrix([sp.Rational(1, 2), r3 / 2])]
    d = [sp.Integer(1), -sp.Rational(1, 2), -sp.Rational(1, 2)]
    h = sp.Rational(1, 2) ** nn
    D = sp.zeros(2, 1)
    for c in range(3):
        o = [(c + 1) % 3, (c + 2) % 3]
        vals = Mn * sp.Matrix([d[c], d[o[0]], d[o[1]]])
        q = 2 * vals[0] - vals[1] - vals[2]
        Em = sp.Matrix([[h * (Pc[o[0]] - Pc[c])[0], h * (Pc[o[0]] - Pc[c])[1]],
                        [h * (Pc[o[1]] - Pc[c])[0], h * (Pc[o[1]] - Pc[c])[1]]])
        gvec = Em.inv() * sp.Matrix([vals[1] - vals[0], vals[2] - vals[0]])
        D += -q * gvec
    Qx = -sp.Rational(3, 8)
    Qy = -r3 / 8
    kap = sp.simplify(D[0] / Qx)
    claimed = -9 * (sp.Rational(18, 25) ** nn - sp.Rational(6, 25) ** nn)
    # symbolic in n (powsimp merges (3/5)^n 2^n ... into single bases), plus n = 1..11 by substitution
    ok = sp.simplify(sp.expand(sp.powsimp(sp.expand(kap - claimed), force=True))) == 0
    ok = ok and all(sp.nsimplify(sp.simplify((kap - claimed).subs(nn, k))) == 0 for k in range(1, 12))
    ok = ok and all(sp.simplify((D[1] - kap * Qy).subs(nn, k)) == 0 for k in range(1, 12))
    return kap, claimed, ok


def alpha_one_exact_level(n):
    """Exact rational Dirichlet solve of L at level n, kappa_n(1) from D_x / Q_x."""
    P, adj, C, lat = gasket_setup(n)
    N = len(adj)
    I = [i for i in range(N) if i not in C]
    d = {C[0]: Fraction(1), C[1]: Fraction(-1, 2), C[2]: Fraction(-1, 2)}
    Lm = sp.zeros(len(I), len(I))
    rhs = sp.zeros(len(I), 1)
    ix = {v: t for t, v in enumerate(I)}
    for v in I:
        Lm[ix[v], ix[v]] = len(adj[v])
        for w in adj[v]:
            if w in ix:
                Lm[ix[v], ix[w]] -= 1
            else:
                rhs[ix[v]] += sp.Rational(d[w].numerator, d[w].denominator)
    sol = Lm.LUsolve(rhs)
    u = {v: sp.Rational(d[v].numerator, d[v].denominator) for v in C}
    for v in I:
        u[v] = sol[ix[v]]
    # x-component of the cell gradient: positions x = a + b/2, y = b sqrt3/2
    Dx = 0
    for c in C:
        j1, j2 = adj[c]
        ev = []
        for j in (j1, j2):
            da, db = lat[j][0] - lat[c][0], lat[j][1] - lat[c][1]
            ev.append((sp.Rational(da.numerator, da.denominator) + sp.Rational(db.numerator, db.denominator) / 2,
                       sp.Rational(db.numerator, db.denominator) * sp.sqrt(3) / 2))
        du = [u[j1] - u[c], u[j2] - u[c]]
        gx = sp.simplify((du[0] * ev[1][1] - du[1] * ev[0][1]) / (ev[0][0] * ev[1][1] - ev[0][1] * ev[1][0]))
        q = sum(u[c] - u[w] for w in adj[c])
        Dx += -q * gx
    return sp.nsimplify(sp.simplify(Dx / (-sp.Rational(3, 8))))


# ----------------------------------------------------------------- torus control
def certify_torus(m=6, alpha=Fraction(45, 100), rule="centred", data=(Fraction(1), Fraction(-1, 2), Fraction(-1, 2)),
                  device=((0, 0), (2, 0), (0, 2)), dps=80):
    iv.dps = dps
    nb = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, 1), (1, -1)]
    idx = lambda i, j: (i % m) * m + (j % m)
    adj = [[] for _ in range(m * m)]
    for i in range(m):
        for j in range(m):
            for a, b in nb:
                adj[idx(i, j)].append(idx(i + a, j + b))
    g = Graph(adj)
    mpoly = krylov_minpoly(g)
    K = krylov_columns(g, len(mpoly) - 1)
    ann = verify_annihilates(K, mpoly)
    roots = isolate_roots(mpoly)
    r3 = iv.sqrt(iv.mpf(3))
    vec = lambda a, b: (iv.mpf(a) + iv.mpf(b) / 2, iv.mpf(b) * r3 / 2)
    coord = {idx(i, j): (i, j) for i in range(m) for j in range(m)}

    def grad(u, c):
        i, j = coord[c]
        if rule == "centred":  # least squares over the six neighbours = (1/3) sum_e e (u(x+e)-u(x))
            gx = sum(((u[idx(i + a, j + b)] - u[c]) * vec(a, b)[0] for a, b in nb), iv.mpf(0)) / 3
            gy = sum(((u[idx(i + a, j + b)] - u[c]) * vec(a, b)[1] for a, b in nb), iv.mpf(0)) / 3
            return (gx, gy)
        return cell_gradient([vec(1, 0), vec(0, 1)], [u[idx(i + 1, j)] - u[c], u[idx(i, j + 1)] - u[c]])

    dev = [idx(*p) for p in device]
    out = maxwell_reading(K, roots, ivq(alpha), dev, [ivq(x) for x in data], grad)
    out.update(ann=ann, deg=len(mpoly) - 1)
    return out


# ----------------------------------------------------------------- float table
def float_kappa(n, alpha=0.45):
    P, adj, C, lat = gasket_setup(n)
    N = len(adj)
    L = np.zeros((N, N))
    for i in range(N):
        L[i, i] = len(adj[i])
        for j in adj[i]:
            L[i, j] -= 1
    w, V = np.linalg.eigh(L)
    w[0] = 0.0
    keep = w > 1e-9
    A = (V * np.where(keep, np.abs(w) ** alpha, 0.0)) @ V.T
    I = [i for i in range(N) if i not in C]
    res = []
    for d in (np.array([1, -0.5, -0.5]), np.array([1, -0.5, 0.0])):
        u = np.zeros(N)
        u[C] = d
        u[I] = np.linalg.solve(A[np.ix_(I, I)], -A[np.ix_(I, C)] @ d)
        q = (A @ u)[C]
        D = np.zeros(2)
        for k, c in enumerate(C):
            Em = np.array([P[j] - P[c] for j in adj[c]])
            D -= q[k] * np.linalg.solve(Em, np.array([u[j] - u[c] for j in adj[c]]))
        x = d - d.mean()
        W = P[C] - P[C].mean(0)
        Q = (x[:, None] ** 2 * W).sum(0)
        res.append((D, Q))
    (D1, Q1), (D2, Q2) = res
    kap = float(D1 @ Q1 / (Q1 @ Q1))
    return kap, float(np.linalg.norm(D2 - kap * Q2) / np.linalg.norm(D2))


def s(x, k=24):
    lo, hi = x._mpi_
    mp.prec = max(mp.prec, iv.prec)
    return f"[{nstr(mp.make_mpf(lo), k)}, {nstr(mp.make_mpf(hi), k)}]"


def main(levels=(1, 2, 3, 4), float_levels=6) -> int:
    t0 = time.time()
    print("Theorem BA: Maxwell-form reading on the 0.45 mesh (L^0.45 gasket). Not thrust.")
    # symbolic |Q| identity
    a, b = sp.symbols("a b", real=True)
    x = [a, b, -a - b]
    W = [(-sp.Rational(1, 2), -sp.sqrt(3) / 6), (sp.Rational(1, 2), -sp.sqrt(3) / 6), (0, sp.sqrt(3) / 3)]
    Qx = sum(x[c] ** 2 * W[c][0] for c in range(3))
    Qy = sum(x[c] ** 2 * W[c][1] for c in range(3))
    check(sp.simplify(Qx**2 + Qy**2 - (sum(t**2 for t in x)) ** 2 / 12) == 0,
          "|Q(d)| = |d-dbar|^2/(2 sqrt3): Q vanishes only for constant data")

    kappas = {}
    for n in levels:
        t = time.time()
        r = certify_gasket(n)
        kap = r["kappa"]
        kappas[n] = kap
        print(f"level {n}: N={r['N']}, deg minpoly={r['deg']}, kappa_n(0.45) in {s(kap)} (width {float(kap.delta):.1e}, {time.time()-t:.1f}s)")
        check(r["ann"], f"level {n}: Krylov polynomial annihilates L_n exactly (spec(L_n) inside its roots)")
        check(r["qint"] < 1e-30, f"level {n}: (A u)_i = 0 on the interior (Krawczyk enclosure, max {float(r['qint']):.1e})")
        check(contains0(sum(r["q"], iv.mpf(0))), f"level {n}: corner charges neutral (Theorem H)")
        check(contains0(r["M"][0]) and contains0(r["M"][1]), f"level {n}: Maxwell (antisymmetric pair) piece sums to 0")
        check(not disjoint(r["R"][0], r["D"][0]) and not disjoint(r["R"][1], r["D"][1]),
              f"level {n}: residual total = D (u = A^+ q + const)")
        check(contains0(r["cross"]), f"level {n}: D parallel to Q(d) (S3 Schur)")
        check(kap.b < 0, f"level {n}: kappa_n(0.45) < 0 certified, so D != 0")
        print(f"   D = ({s(r['D'][0], 20)}, {s(r['D'][1], 20)}) for d=(1,-1/2,-1/2)")
        print(f"   self (diagonal) part kappa_self in {s(r['kappa_self'], 16)}; corner charges q_0 in {s(r['q'][0], 16)}")
        fx = [r["f"][k][0] for k in range(3)]
        check(disjoint(fx[0], r["D"][0]) and disjoint(fx[0] + fx[1], r["D"][0]),
              f"level {n}: sets enclosing one or two corners give other values than D (surface-dependent below C)")
    check(all(kappas[n + 1].b < kappas[n].a for n in levels if n + 1 in kappas),
          "|kappa_n(0.45)| strictly increases on the certified levels")

    # alpha = 1 control
    kap1, claimed, ok = alpha_one_symbolic()
    check(ok, "alpha=1: kappa_n(1) = -9((18/25)^n - (6/25)^n) symbolically in n (harmonic extension matrix M^n)")
    for n in (1, 2, 3):
        ex = alpha_one_exact_level(n)
        check(ex == -9 * (Fraction(18, 25) ** n - Fraction(6, 25) ** n), f"alpha=1 level {n}: exact solve gives kappa={ex}")
    r1 = certify_gasket(2, alpha=Fraction(1), dps=60)
    check(contains0(r1["kappa"] - ivq(-9 * (Fraction(18, 25) ** 2 - Fraction(6, 25) ** 2))),
          "alpha=1 level 2: interval pipeline reproduces -9(0.72^2-0.24^2)")

    # torus control
    for rule in ("centred", "onesided"):
        for data in ((Fraction(1), Fraction(-1, 2), Fraction(-1, 2)), (Fraction(1), Fraction(-1, 2), Fraction(0))):
            rt = certify_torus(rule=rule, data=data)
            Dx, Dy = rt["D"]
            print(f"torus 6x6 {rule:8s} d={tuple(str(v) for v in data)}: D = ({s(Dx, 14)}, {s(Dy, 14)})")
            check(rt["ann"], "torus: Krylov polynomial annihilates L exactly")
            check(contains0(rt["M"][0]) and contains0(rt["M"][1]), f"torus {rule}: Maxwell piece sums to 0")
            if rule == "centred":
                check(contains0(Dx) and contains0(Dy) and max(Dx.delta, Dy.delta) < 1e-40,
                      f"torus centred: D = 0 (odd stencil, even translation-invariant kernel)")
            else:
                check(excludes0(Dx) and excludes0(Dy), "torus one-sided cell stencil: D != 0 certified (Kept Failure BA.1)")
    rt3 = certify_torus(rule="centred", data=(Fraction(1), Fraction(-1, 2), Fraction(0)), device=((0, 0), (3, 0), (1, 2)))
    check(contains0(rt3["D"][0]) and contains0(rt3["D"][1]), "torus centred, scalene device (0,0),(3,0),(1,2): D = 0")

    # float table
    prev = None
    for n in range(1, float_levels + 1):
        k, rel = float_kappa(n)
        print(f"float level {n}: kappa_n(0.45) = {k:.12f}, 2^-n kappa_n = {k / 2**n:.9f}, Schur defect (second datum) {rel:.1e}")
        check(rel < 1e-9, f"float level {n}: second datum (1,-1/2,0) obeys D = kappa Q(d)")
        if n in kappas:
            check(abs(k - float(kappas[n].mid)) < 1e-9, f"float level {n} matches the certificate")
        if prev is not None:
            check(abs(k) / 2**n < abs(prev) / 2 ** (n - 1), f"float level {n}: |2^-n kappa_n| decreases")
            check(abs(k) > abs(prev), f"float level {n}: |kappa_n| increases")
        prev = k
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BA checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    lv = tuple(range(1, int(sys.argv[1]) + 1)) if len(sys.argv) > 1 else (1, 2, 3, 4)
    sys.exit(main(levels=lv))
