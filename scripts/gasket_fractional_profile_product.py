#!/usr/bin/env python3
r"""Theorem CB: an exact infinite-product formula for the corner resolvent and for 1/Q_inf, and the exact atom set of
BX's measure nu off -Z.  With p(z) = z(z+5), H := 1/r - z/2 obeys the purely multiplicative law H(z) = c(z) H(p(z));
since c(w) p(w)/w = (w+6)/(w+2), this telescopes to

    1/r_n(s)   = (s/2) [1 + (1 + 6/p^n(s)) prod_(j<n) (p^j s + 6)/(p^j s + 2)]      (every level n >= 0, exact),
    1/r_inf(z) = (z/2) [1 + P(z)],   P(z) = prod_(j>=0) (p^j z + 6)/(p^j z + 2)        (z in C \ K_p),
    1/Q_inf(u) = (psi(u)/2) prod_(j>=0) (psi(5^j u)+6)/(psi(5^j u)+2) * prod_(j>=1) (5/3) c(psi(5^-j u))   (u in C \ Z),

i.e. BZ.3's S equals Pi psi/2 identically.  Consequences: the atoms of nu in (0,inf) \ (-Z) are exactly Lambda_O
(Conjecture CA.6 holds off -Z; what remains of CA.6 is nu(-Z) = 0); no lam with psi(-5^k lam) in the -2 branch is an
atom (CA.5 for the whole branch); and the atom masses have a closed product form in the orbit psi(5^i u) alone.

Setting (GASKET.md, Theorems AQ, BD, BJ, BS, BX, BZ, CA): r_n(s) = (1/2) v^T (L_n + s)^-1 v (corner zero-sum v,
|v|^2 = 2), r_(n+1)(s) = Phi_s(r_n(p(s))), r_0 = 1/(3+s) (AQ); c(z) = (z+6)/((z+2)(z+5)), e(z) = 2z/(2+z),
D = 1/r_inf, D = c D(p) - e on C \ K_p (BZ.2); psi (BD.4); f = 1/Q_inf = Pi D(psi) - S (BZ.3);
K_p (filled Julia set of p, a null Cantor set in [-5,0], BZ.1), Z = psi^-1(K_p); O = U_k p^-k(-3);
Lambda_O = {lam > 0 : psi(-5^k lam) = -3 for some k in Z} (CA).

Theorem CB (checked here; exit 1 on any failure):
  CB.1 [proved; sympy + exact rationals] H_n := 1/r_n - s/2 satisfies H_(n+1)(s) = c(s) H_n(p(s)), H_0 = 3 + s/2,
       hence the finite-level product.  Checks: sympy against AQ's recursion for n <= 4; exact rational resolvents
       of the actual gasket graphs G_1, G_2, G_3 (built here from vertex coordinates; no AQ input) at s = 1/2, 3.
  CB.2 [proved from CB.1, AQ (3), BZ.1-BZ.2] 1/r_inf = (z/2)(1 + P(z)) on C \ K_p; P is meromorphic there (p^j z -> inf
       locally uniformly, BZ.1), zeros exactly at U_j p^-j(-6), poles at U_j p^-j(-2) (cancelling on the -2 branch).
       arb: agrees with BZ's D_of (an independent enclosure using the far-field bound) to 1e-50 at real and complex points.
  CB.3 [proved from CB.2, BZ.3] S(u) = Pi(u) psi(u)/2 identically, so f = Pi (psi/2) P(psi) on C \ Z.
       arb: agrees with BZ's inv_Q to 1e-40; reproduces BW's G* = F(t_max), min F = F(t_min) to 1e-38 without D_of.
  CB.4 [proved from CB.3, BZ.4] for u < 0, u not in Z: f(u) = 0 iff psi(5^k u) = -3 for some k in Z.  Hence the atoms
       of nu in (0,inf) \ (-Z) are exactly Lambda_O, and no lam with psi(-5^k lam) in U_j p^-j(-2) is an atom.
       Certified instances: CA.5's witness (psi(u) = -2) re-derived from the product: f(u) in [2.1979467595912,
       2.1979467595913]; a second -2-branch point psi(u) = (-5+sqrt17)/2 has f(u) != 0.
  CB.5 [proved from CB.3; certified] if psi(5^k u) = -3 (k >= 0), the atom lam = -u has
         nu({lam}) = -4 / [ Pi(u) (x_0/2) prod_(i>=0, i != k+1) (x_i+6)/(x_i+2) * prod_(i=0..k) (2x_i+5) * psi'(u) ],
       x_i = psi(5^i u), psi'(u) = prod_(i>=1) (1 + 2 psi(5^-i u)/5).  arb: reproduces BZ.5's masses of the three
       largest atoms (0.479640692929..., 0.0359221043246..., 0.0136278303937...) and BZ's own certify_atom.
Conjecture CA.6 is reduced, not proved: what remains open is nu(-Z) = 0 (no atom on -Z, no singular continuous part).
Assumptions used as stated (not re-run): AQ (1)-(3), BD.4, BJ.6, BS.1, BX.1-BX.2, BZ.1-BZ.4, CA.1-CA.2 (the -3 branch).
BW's digits of G*, min F, t_max, t_min and BZ.5's printed masses for consistency only.  No literature assumption.
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction

import sympy as sp
from flint import acb, arb, ctx

import gasket_fractional_profile_atoms as bz
import gasket_fractional_profile_atom_law as ca

FAILS: list[str] = []
J = 200
W_ESC = 60
BZ_MASSES = {0: ("3.61525671988770", "0.479640692929"), 2: ("1.171910250998837", "0.0359221043246"),
             3: ("1.016316417446789", "0.0136278303937")}
F_TWO = ca.F_TWO


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- symbolic
def r_closed_sym(n, s):
    pj, prod = s, sp.Integer(1)
    for _ in range(n):
        prod *= (pj + 6) / (pj + 2)
        pj = pj * (pj + 5)
    return 1 / (s / 2 * (1 + (1 + 6 / pj) * prod))


def symbolic_checks():
    z, w, Hn, x = sp.symbols("z w Hn x")
    c = (z + 6) / ((z + 2) * (z + 5))
    e = 2 * z / (2 + z)
    pz = z * (z + 5)
    # H = D - z/2: D' = c D(p) - e  ->  H' = c H(p)  because c p/2 - e = z/2
    check(sp.simplify(c * (Hn + pz / 2) - e - z / 2 - c * Hn) == 0, "CB.1: D = cD(p) - e  <=>  H = cH(p) for H = D - z/2 (c p/2 - e = z/2)")
    check(sp.simplify(c * pz / z - (z + 6) / (z + 2)) == 0, "CB.1: c(w) p(w)/w = (w+6)/(w+2) (telescoping factor)")
    check(sp.simplify((3 + z) - z / 2 - z / 2 * (1 + 6 / z)) == 0, "CB.1: H_0 = 1/r_0 - s/2 = 3 + s/2 = (s/2)(1 + 6/s)")
    s = sp.Symbol("s", positive=True)
    Phi = lambda X: (2 + s) * (5 + s) * X / ((6 + s) - 2 * s * (5 + s) * X)
    r = 1 / (3 + s)
    ok = True
    for n in range(0, 5):
        ok = ok and sp.cancel(r - r_closed_sym(n, s)) == 0
        r = sp.cancel(Phi(r.subs(s, s * (s + 5))))
    check(ok, "CB.1: 1/r_n = (s/2)[1 + (1+6/p^n s) prod_(j<n)(p^j s+6)/(p^j s+2)] equals AQ's recursion r_(n+1) = Phi_s(r_n(p s)) for n = 0..4 (sympy)")
    check(sp.expand(pz + 6 - (z + 2) * (z + 3)) == 0, "CB.2/CB.4: p(z)+6 = (z+2)(z+3), so on the -2 branch the zero of (p z+6) cancels the pole of 1/(z+2)")
    check(sp.simplify((z + 6) / (z + 2) * (pz + 6) / (pz + 2) - (z + 6) * (z + 3) / (pz + 2)) == 0,
          "CB.2: (z+6)(pz+6)/((z+2)(pz+2)) = (z+6)(z+3)/(pz+2) (regular pair at z = -2, value -1)")
    roots = sp.solve(sp.Eq(pz, -2), z)
    check(set(roots) == {sp.Rational(-5, 2) - sp.sqrt(17) / 2, sp.Rational(-5, 2) + sp.sqrt(17) / 2} and all(rt not in (-2, -3, -5, -6) for rt in roots),
          "CB.4: p^-1(-2) = (-5 +- sqrt17)/2, neither -2, -3, -5 nor -6")
    check(pz.subs(z, -6) == 6 and all(v > 0 for v in (6, 66)), "CB.4: the forward orbit of -6 is 6, 66, ... > 0 (never -6, -2 or -5 again)")
    check(sp.diff((z + 6) / (z + 2), z).subs(z, -6) == sp.Rational(-1, 4) and sp.diff(sp.Rational(5, 3) * c, z).subs(z, -6) == sp.Rational(5, 12),
          "CB.5: d/dx (x+6)/(x+2) = -1/4 and d/dx (5/3)c(x) = 5/12 at x = -6")
    psi = sp.Function("psi")
    check(sp.simplify(sp.diff(psi(w) * (psi(w) + 5), w) - (2 * psi(w) + 5) * sp.diff(psi(w), w)) == 0,
          "CB.5: psi(5w) = p(psi(w)) gives 5 psi'(5w) = (2 psi(w)+5) psi'(w), so psi'(u) = psi'(5^-J u) prod_(i<=J)(1 + 2 psi(5^-i u)/5)")
    # tail constants
    check(sp.Rational(4, 8) <= sp.Rational(1, 2) and 10 * (10 - 5) >= 2 * 10,
          "CB.2: for |w| >= 10: |4/(w+2)| <= 1/2, |log(1+t)| <= 2|t|, |p(w)| >= 2|w|; tail of P bounded by exp(16/(|w|-2))")
    t = sp.Symbol("t", positive=True)
    check(sp.limit((2 * sp.exp(t / 5) - 1 - 1) / t, t, 0) == sp.Rational(2, 5),
          "CB.5: psi'(t) <= (psi(2t)-psi(t))/t <= 2e^(t/5) - 1 (convexity, BJ.6, psi >= t), so |psi'(v) - 1| <= 2(e^(|v|/5) - 1)")


# ----------------------------------------------------------------- the actual graphs (independent of AQ)
def gasket_graph(n):
    tri = [((0, 0), (1, 0)), ((1, 0), (0, 1)), ((0, 1), (0, 0))]
    edges = set(tri)
    size = 1
    for _ in range(n):
        new = set()
        for off in ((0, 0), (size, 0), (0, size)):
            for (a, b) in edges:
                new.add(((a[0] + off[0], a[1] + off[1]), (b[0] + off[0], b[1] + off[1])))
        edges, size = new, 2 * size
    verts = sorted({v for ed in edges for v in ed})
    idx = {v: i for i, v in enumerate(verts)}
    und = {tuple(sorted((idx[a], idx[b]))) for a, b in edges}
    return len(verts), und, idx[(0, 0)], idx[(size, 0)], idx[(0, size)]


def resolvent_exact(n, s):
    N, E, a, b, c3 = gasket_graph(n)
    A = [[Fraction(0)] * N for _ in range(N)]
    for i, j in E:
        A[i][j] -= 1
        A[j][i] -= 1
        A[i][i] += 1
        A[j][j] += 1
    for i in range(N):
        A[i][i] += s
    rhs = [Fraction(0)] * N
    rhs[a], rhs[b] = Fraction(1), Fraction(-1)
    # Gaussian elimination (exact)
    M = [row[:] + [rhs[i]] for i, row in enumerate(A)]
    for col in range(N):
        piv = next(r for r in range(col, N) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        for r in range(col + 1, N):
            if M[r][col] != 0:
                fct = M[r][col] / pv
                M[r] = [M[r][k] - fct * M[col][k] for k in range(N + 1)]
    xs = [Fraction(0)] * N
    for r in range(N - 1, -1, -1):
        xs[r] = (M[r][N] - sum(M[r][k] * xs[k] for k in range(r + 1, N))) / M[r][r]
    degs = sorted({sum(1 for ed in E if i in ed) for i in range(N)})
    return (xs[a] - xs[b]) / 2, N, degs


def graph_checks():
    ok, info = True, []
    for n in (1, 2, 3):
        for s in (Fraction(1, 2), Fraction(3)):
            rg, N, degs = resolvent_exact(n, s)
            rc = Fraction(str(sp.Rational(r_closed_sym(n, sp.Rational(s.numerator, s.denominator)))))
            ok = ok and rg == rc
        info.append(f"G_{n}: {N} vertices, degrees {degs}")
    print("   " + "; ".join(info))
    check(ok, "CB.1: exact resolvents (1/2) v^T (L_n + s)^-1 v of the gasket graphs G_1, G_2, G_3 (v = e_a - e_b) equal the product formula at s = 1/2, 3")


# ----------------------------------------------------------------- interval product
def _isz(w):
    return w.contains(0)


def P_enc(z, cplx=False):
    """Enclosure of P(z) = prod_(j>=0)(p^j z+6)/(p^j z+2) for z off K_p; a factor whose (w+2) encloses 0 is paired
    with the next one in the regular form (w+6)(w+3)/(p(w)+2)."""
    Wb = arb(10) ** W_ESC
    big = (lambda y: abs(y).lower() > Wb) if cplx else (lambda y: abs(y).lower() > Wb)
    w, prod, steps = z, (acb(1) if cplx else arb(1)), 0
    while not big(w):
        if _isz(w + 2):
            pw = bz._p(w)
            bz._nz(pw + 2)
            prod *= (w + 6) * (w + 3) / (pw + 2)
            w = bz._p(pw)
        else:
            prod *= (w + 6) / (w + 2)
            w = bz._p(w)
        steps += 1
        if steps > 400:
            raise bz.Fail("orbit does not escape")
    tl = (16 / (abs(w).lower() - 2))
    tl = arb(tl.upper())
    tail = (acb(arb(0, tl), arb(0, tl)) if cplx else arb(0, tl)).exp()
    return prod * tail


def D_prod(z, cplx=False):
    return z / 2 * (1 + P_enc(z, cplx))


def orbit(u, cplx=False):
    """x[j] = psi(5^-j u) for j = 0..J (BZ's construction) and the Pi tail ball exponent eps."""
    five = arb(5)
    ball = (lambda r: acb(arb(0, r), arb(0, r))) if cplx else (lambda r: arb(0, r))
    v = u / five ** J
    av = arb(abs(v).upper())
    dl = (av * ((av / 10).exp() - 1)).upper()
    x = [None] * (J + 1)
    x[J] = v + ball(dl)
    for j in range(J, 0, -1):
        x[j - 1] = bz._p(x[j])
    eps = (abs(u).upper() * five ** (-J)).upper()
    return x, eps, ball, av


def Pi_enc(x, eps, ball, cplx=False):
    f3 = arb(5) / 3
    Pi = acb(1) if cplx else arb(1)
    for j in range(1, J + 1):
        bz._nz(x[j] + 2)
        bz._nz(x[j] + 5)
        Pi *= f3 * bz._c(x[j])
    return Pi * (acb(arb(0, eps)) if cplx else arb(0, eps)).exp()


def f_prod(u, cplx=False):
    """1/Q_inf(u) = Pi(u) (psi(u)/2) P(psi(u))  (CB.3)."""
    x, eps, ball, _ = orbit(u, cplx)
    return Pi_enc(x, eps, ball, cplx) * x[0] / 2 * P_enc(x[0], cplx)


def dpsi_enc(u):
    """psi'(u) = psi'(5^-J u) prod_(i=1..J)(1 + 2x_i/5), |psi'(v) - 1| <= 2(e^(|v|/5) - 1) (CB.5)."""
    x, _, _, av = orbit(u)
    d = 1 + arb(0, (2 * ((av / 5).exp() - 1)).upper())
    for i in range(1, J + 1):
        d *= 1 + 2 * x[i] / 5
    return d


def mass_closed(u, k):
    """nu({-u}) by CB.5, for u (ball) with psi(5^k u) = -3, k >= 0: the vanishing factor is index k+1 of P(psi(u))."""
    x, eps, ball, _ = orbit(u)
    Pi = Pi_enc(x, eps, ball)
    w = x[0]
    orb = [w]
    for _ in range(k + 1):
        orb.append(bz._p(orb[-1]))
    rest = arb(1)
    for i in range(k + 1):
        rest *= (orb[i] + 6) / (orb[i] + 2)
    rest /= (orb[k + 1] + 2)                         # d/dx (x+6)/(x+2) at x = -6 is 1/(x+2) = -1/4
    rest *= P_enc(bz._p(orb[k + 1]))                  # factors i >= k+2
    dx = dpsi_enc(u)
    for i in range(k + 1):
        dx *= 2 * orb[i] + 5                          # d/du p^(k+1)(psi(u))
    fp = Pi * w / 2 * rest * dx
    return 1 / fp, orb[k]


# ----------------------------------------------------------------- checks
def consistency():
    ctx.prec = 300
    ok = True
    for zs in ("6", "1", "0.01", "-2.3", "-3.5", "-5.5", "-1.5", "-0.3"):
        z = arb(zs)
        ok = ok and bool(abs(D_prod(z) - bz.D_of(z, False)) < arb("1e-50"))
    for zc in (acb(-1, "0.5"), acb("-4.2", "0.01"), acb(3, 7)):
        ok = ok and bool(abs(D_prod(zc, True) - bz.D_of(zc, True)) < arb("1e-50"))
    print(f"   1/r_inf(6) = {D_prod(arb(6)).str(25)},  1/r_inf(1) = {D_prod(arb(1)).str(25)}")
    check(ok, "CB.2 (consistency): (z/2)(1+P(z)) agrees with BZ's D_of to 1e-50 at 8 real and 3 complex points")
    ok = True
    for us in ("1.7", "0.3", "-0.37", "-2.3", "-4.1", "12.5"):
        u = arb(us)
        ok = ok and bool(abs(f_prod(u) - bz.inv_Q(u)) < arb("1e-40"))
    for uc in (acb("-1.3", "0.2"), acb("2", "-3")):
        ok = ok and bool(abs(f_prod(uc, True) - bz.inv_Q(uc, cplx=True)) < arb("1e-40"))
    check(ok, "CB.3 (consistency): Pi (psi/2) P(psi) agrees with BZ's inv_Q = Pi D(psi) - S to 1e-40 (6 real, 2 complex points), i.e. S = Pi psi/2")
    ctx.prec = 400
    g = arb(3).log() / arb(5).log()
    ok = True
    for t, ref in ((bz.T_MAX, bz.G_STAR), (bz.T_MIN, bz.MIN_F)):
        u = arb(5) ** arb(t)
        Fv = u ** (1 - g) / f_prod(u)
        print(f"   F(5^{t[:12]}...) = {Fv.str(40)}  (product formula)")
        ok = ok and bool(abs(Fv - arb(ref)) < arb("1e-38"))
    check(ok, "CB.3 (consistency): the product formula alone gives BW's G* = F(t_max) and min F = F(t_min) to 1e-38")


def branch_checks():
    ctx.prec = 300
    a, b = ca.bracket_root(ca.psi_enc, -2.25, -2.23, arb(-2))
    ball = arb((a + b) / 2, ((b - a) / 2).upper())
    fv = f_prod(ball)
    print(f"   psi(u) = -2 at u in {ball.str(22)}:  f(u) = {fv.str(16)} (product)")
    check(bool(fv > arb(F_TWO[0])) and bool(fv < arb(F_TWO[1])),
          "CB.4: CA.5's witness from the product alone: f(u) in [2.1979467595912, 2.1979467595913] where psi(u) = -2")
    tgt = (-arb(5) + arb(17).sqrt()) / 2
    # locate u < 0 with psi(u) = (-5+sqrt17)/2 = -0.4384... by a float scan, then certify
    import math
    grid = [-0.2 - 0.0005 * i for i in range(4000)]
    vals = [float(ca.psi_enc(arb(g_), J=60).mid()) - float(tgt.mid()) for g_ in grid]
    i0 = next(i for i in range(len(grid) - 1) if vals[i] * vals[i + 1] < 0)
    a, b = ca.bracket_root(ca.psi_enc, grid[i0 + 1], grid[i0], tgt)
    ball = arb((a + b) / 2, ((b - a) / 2).upper())
    fv = f_prod(ball)
    print(f"   psi(u) = (-5+sqrt17)/2 at u in {ball.str(22)}:  f(u) = {fv.str(16)}")
    check(not fv.contains(0) and math.isfinite(float(fv.mid())),
          "CB.4: a second -2-branch point (psi(u) in p^-1(-2)) has f(u) != 0: no atom there either")


def mass_checks():
    ctx.prec = 300
    ok, msgs = True, []
    brackets = {0: (-3.62, -3.61), 2: (-1.1769, -1.1669), 3: (-1.0213, -1.0113)}
    for k, (lo, hi) in brackets.items():
        tgt = arb(-3)
        fun = (lambda kk: (lambda uu: _pk(ca.psi_enc(uu), kk)))(k)
        a, b = ca.bracket_root(fun, lo, hi, tgt)
        u = arb((a + b) / 2, ((b - a) / 2).upper())
        m, xk = mass_closed(u, k)
        lam_ref, m_ref = BZ_MASSES[k]
        ok = ok and bool(abs(-u - arb(lam_ref)) < arb("1e-13")) and bool(abs(m - arb(m_ref)) < arb("2e-12")) and bool(m > 0)
        msgs.append(f"k={k}: lam = {(-u).str(18)}, nu = {m.str(20)}")
    print("   " + "\n   ".join(msgs))
    check(ok, "CB.5: the closed product form of nu({lam}) reproduces BZ.5's masses of the three largest atoms (k = 0, 2, 3) to 2e-12, all > 0")
    five = arb(5)
    lo, hi = (-five ** (arb(3195) / 4000)).mid(), (-five ** (arb(3194) / 4000)).mid()
    lam, wbz = bz.certify_atom(lo, hi)
    a, b = ca.bracket_root(ca.psi_enc, -3.62, -3.61, arb(-3))
    m, _ = mass_closed(arb((a + b) / 2, ((b - a) / 2).upper()), 0)
    print(f"   largest atom: BZ certify_atom nu = {wbz.str(20)};  CB.5 nu = {m.str(30)}")
    check(bool(abs(m - wbz) < 10 * wbz.rad() + arb("1e-25")), "CB.5: agrees with BZ's Cauchy-estimate certificate for the largest atom")


def _pk(x, k):
    for _ in range(k):
        x = bz._p(x)
    return x


def main() -> int:
    t0 = time.time()
    print("Theorem CB: 1/r = (z/2)(1 + prod (p^j z+6)/(p^j z+2)) exactly; 1/Q_inf = Pi (psi/2) P(psi); atoms of nu off -Z are "
          "exactly Lambda_O; closed product form of the atom masses. Not thrust.")
    symbolic_checks()
    graph_checks()
    consistency()
    branch_checks()
    mass_checks()
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CB checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
