#!/usr/bin/env python3
r"""Theorem CA: an exact location law for atoms of BX's measure nu.  On the backward orbit O = U_k p^-k(-3) of -3
(the -3 branch of p^-1(-6)), 1/r_inf = D is exactly z/2; consequently every lam > 0 with psi(-5^k lam) = -3 for some
integer k is an atom of nu.  These atoms accumulate at every point -x of -Z at which psi - psi(x) changes sign; in
particular nu has infinitely many atoms in one scale period, accumulating at lam* = 1.0546444769654709... (psi(-lam*)
= -1, p(-1) = -4 = fixed point).  All 311 atoms certified by Theorem BZ obey the law.  Kept Failure CA.5: the other
branch -2 of p^-1(-6) carries no atom (lam = 2.2421331852464538... with psi(-lam) = -2 has 1/Q_inf(-lam) =
2.19794675959126... != 0).

Setting (GASKET.md, Theorems AQ, BD, BS, BX, BZ): p(z) = z(z+5), c(z) = (z+6)/((z+2)(z+5)), e(z) = 2z/(2+z),
D = 1/r_inf with D(z) = c(z) D(p(z)) - e(z) on C \ K_p (BZ.2), psi (BD.4), Q_inf = int dnu/(lam+u) (BX.1),
f = 1/Q_inf = Pi D(psi) - S (BZ.3), Z = psi^-1(K_p), b(x) = (5/3)c(x), x_j = psi(5^-j u).

Theorem CA (checked here; exit 1 on any failure):
  CA.1 [proved; sympy] O = U_(k>=0) p^-k(-3) is a subset of (-5,0) \ K_p avoiding -2, -5, -6, and
       D(z) = z/2 for every z in O u {-6}  (c(-6) = 0 gives D(-6) = -e(-6) = -3; c(z) p(z)/2 - e(z) = z/2).
  CA.2 [proved from CA.1, BZ.3, BZ.4] if u < 0 and psi(u) in O, then u is not in Z and f(u) = 0: with X_M the
       unrolled recursion, f(u) = Pi(5^-M u) (5/3)^M x_M/2 - sum_(k>M)(5/3)^k e(x_k) Pi(5^-k u) -> 0 since
       |x_M| <= 2 5^-M |u|.  So lam = -u is an atom of nu, and by nu(5A) = 3 nu(A) so is every
       lam in Lambda_O = {lam > 0 : psi(-5^k lam) = -3 for some k in Z}.
  CA.3 [proved + certified] each level-n piece of K_p (BZ.1) contains a point of O, so for every x in Z at which
       psi - psi(x) changes sign, -x is a limit of points of Lambda_O.  Certified instance: psi(-1.055) < -1 <
       psi(-1.054) (arb), -1 in K_p, and h_n = g_+(g_-^n(-3)) in O tends to -1; so nu has infinitely many atoms in
       (1.054, 1.055), accumulating at lam* (psi(-lam*) = -1).  The first two are BZ's lam = 1.1719102509988...
       and 1.0163164174467...; lam_1 = 3.6152567198877... has psi(-lam_1) = -3.
  CA.4 [certified, arb] re-running BZ.5's scan and certificate: each of the 311 certified atoms lam_j in (1,5)
       lies in Lambda_O: an interval sign change of psi(-5^k_j lam) + 3 inside the box where f is holomorphic and
       strictly increasing (so has exactly one zero) identifies lam_j with a point of Lambda_O.
  Kept Failure CA.5 (certified witness): "every lam with psi(-5^k lam) in p^-1(-6) = {-2,-3} is an atom" is false:
       at u with psi(u) = -2 (u in [-2.24213318524645386322..]), f(u) in [2.1979467595912, 2.1979467595913] != 0
       (D at -2 by the removable form D(z) = (z+6)(z+3)D(p^2 z)/((z+5)(pz+2)(pz+5)) - 2z(z+4)/(pz+2), sympy).
Conjecture CA.6 (not claimed): the atoms of nu are exactly Lambda_O, and nu is carried by them (pure point).
Assumptions used as stated (not re-run): AQ (1)-(3), BD.4, BJ.6, BS.1, BX.1-BX.2, BZ.1-BZ.4 (in particular the
formula f = Pi D(psi) - S, its enclosure routine inv_Q, imported unchanged, and BZ.4's description of nu off -Z).
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time

import sympy as sp
from flint import arb, ctx

import gasket_fractional_profile_atoms as bz

FAILS: list[str] = []
N_ATOMS_MIN = 311
K_MAX = 14
LAM_STAR = ("1.05464447696547", "1.05464447696548")
U_TWO = "-2.24213318524645386322626210854"
F_TWO = ("2.1979467595912", "2.1979467595913")
BZ_LAMS = ("1.171910250998837", "1.016316417446789")


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- symbolic
def symbolic_checks():
    z, w, Dpp = sp.symbols("z w Dpp")
    c = (z + 6) / ((z + 2) * (z + 5))
    e = 2 * z / (2 + z)
    pz = z * (z + 5)
    check(sp.simplify(c * pz / 2 - e - z / 2) == 0, "CA.1: c(z) p(z)/2 - e(z) = z/2 (D(w) = w/2 at w = p(z) propagates to D(z) = z/2)")
    check(c.subs(z, -6) == 0 and -e.subs(z, -6) == -3, "CA.1: c(-6) = 0, so D(-6) = -e(-6) = -3 = (-6)/2 (D(6) = 1/r_inf(6) finite)")
    check(set(sp.solve(sp.Eq(pz, -6), z)) == {-2, -3} and set(sp.solve(sp.Eq(pz, -3), z)) == {sp.Rational(-5, 2) - sp.sqrt(13) / 2, sp.Rational(-5, 2) + sp.sqrt(13) / 2},
          "CA.1: p^-1(-6) = {-2, -3}; p^-1(-3) = -5/2 +- sqrt(13)/2 (real, in (-5,0))")
    # p^-1([-5,0]) in [-5,0]: roots -5/2 +- sqrt(25/4 + w), 0 <= 25/4 + w <= 25/4
    check(sp.Rational(25, 4) - 5 > 0 and sp.sqrt(sp.Rational(25, 4)) == sp.Rational(5, 2),
          "CA.1: for w in [-5,0], p^-1(w) = -5/2 +- sqrt(25/4+w) is real and lies in [-5,0]")
    check(pz.subs(z, -2) == -6 and pz.subs(z, -5) == 0 and pz.subs(z, -6) == 6 and pz.subs(z, -3) == -6,
          "CA.1: p(-2) = -6, p(-5) = 0, p(-6) = 6 are not in O (O lies in [-5,0], 0 is fixed), so O avoids -2, -5, -6")
    check(pz.subs(z, -3) < -5, "CA.1: p(-3) = -6 < -5, so -3 escapes and O misses K_p (p^-1(K_p) = K_p)")
    # unrolled recursion with D(x_k) = x_k/2 (generic b_j, e_j, M = 3) -- the identity used in CA.2
    b = sp.symbols("b1:4")
    ee = sp.symbols("e1:4")
    X0 = sp.Symbol("X0")
    X = X0
    for k in range(1, 4):
        X = b[k - 1] * X - sp.Rational(5, 3) ** k * ee[k - 1]
    closed = b[0] * b[1] * b[2] * X0 - sum(sp.Rational(5, 3) ** k * ee[k - 1] * sp.prod(b[k:3]) for k in range(1, 4))
    check(sp.expand(X - closed) == 0, "CA.2: (5/3)^M D(x_M) = prod b_j D(x_0) - sum (5/3)^k e_k prod_(j>k) b_j (BZ.3's unrolling, M = 3)")
    # removable form of D at -2 (two steps of BZ.2)
    Dp = c.subs(z, pz) * Dpp - e.subs(z, pz)          # D(p z) = c(pz) D(p^2 z) - e(pz)
    lhs = c * Dp - e
    rhs = (z + 6) * (z + 3) * Dpp / ((z + 5) * (pz + 2) * (pz + 5)) - 2 * z * (z + 4) / (pz + 2)
    check(sp.simplify(lhs - rhs) == 0, "CA.5: D(z) = (z+6)(z+3)D(p^2 z)/((z+5)(pz+2)(pz+5)) - 2z(z+4)/(pz+2), regular at z = -2")
    check(sp.simplify(rhs.subs(z, -2) - (Dpp / 3 - 2)) == 0, "CA.5: D(-2) = D(6)/3 - 2")
    gp = sp.Rational(-5, 2) + sp.sqrt(w + sp.Rational(25, 4))
    gm = sp.Rational(-5, 2) - sp.sqrt(w + sp.Rational(25, 4))
    check(gm.subs(w, -4) == -4 and gp.subs(w, -4) == -1 and pz.subs(z, -1) == -4 and sp.diff(gm, w).subs(w, -4) == sp.Rational(-1, 3),
          "CA.3: g_-(-4) = -4 (fixed, g_-'(-4) = -1/3), g_+(-4) = -1, p(-1) = -4: -1 in K_p")


# ----------------------------------------------------------------- interval tools
def psi_enc(u, J=200):
    """Outward enclosure of psi(u) for a real ball u (BD.4, BJ.6), as in BZ's inv_Q."""
    v = u / arb(5) ** J
    av = arb(abs(v).upper())
    dl = (av * ((av / 10).exp() - 1)).upper()
    x = v + arb(0, dl)
    for _ in range(J):
        # p(x): x(x+5) while |x| < 1/2 (no cancellation), else (x+5/2)^2 - 25/4 (one occurrence of x, no dependency)
        x = x * (x + 5) if bool(abs(x) < arb(1) / 2) else (x + arb(5) / 2) ** 2 - arb(25) / 4
    return x


def inv_Q_reg(u, J=bz.J_TAIL):
    """BZ's enclosure of f = 1/Q_inf (copied from bz.inv_Q), with D(x_0) by CA.5's removable form near -2."""
    five = arb(5)
    f3 = arb(5) / 3
    v = u / five ** J
    av = arb(abs(v).upper())
    dl = (av * ((av / 10).exp() - 1)).upper()
    x = [None] * (J + 1)
    x[J] = v + arb(0, dl)
    for j in range(J, 0, -1):
        x[j - 1] = bz._p(x[j])
    for j in range(1, J + 1):
        bz._nz(x[j] + 2)
        bz._nz(x[j] + 5)
    bb = [None] + [f3 * bz._c(x[j]) for j in range(1, J + 1)]
    eps = (abs(u).upper() * five ** (-J)).upper()
    tail = arb(0, eps).exp()
    Pi = arb(1)
    for j in range(1, J + 1):
        Pi *= bb[j]
    S, tl = arb(0), tail
    for k in range(J, 0, -1):
        S += f3 ** k * bz._e(x[k]) * tl
        tl *= bb[k]
    S += arb(0, (arb(42) / 10 * abs(u).upper() * arb(3) ** (-J)).upper())
    z0 = x[0]
    pz = bz._p(z0)
    for d in (z0 + 5, pz + 2, pz + 5):
        bz._nz(d)
    Dz = (z0 + 6) * (z0 + 3) * bz.D_of(bz._p(pz), False) / ((z0 + 5) * (pz + 2) * (pz + 5)) - 2 * z0 * (z0 + 4) / (pz + 2)
    return Pi * tail * Dz - S


def bracket_root(fun, lo, hi, target, it=280):
    """Certified a < b with (fun(a)-target)(fun(b)-target) < 0, by bisection on midpoints (lo, hi floats/arb)."""
    a, b = arb(lo), arb(hi)
    sa = fun(a) - target
    sb = fun(b) - target
    if not ((sa < 0 and sb > 0) or (sa > 0 and sb < 0)):
        raise bz.Fail("no sign change")
    neg_left = bool(sa < 0)
    for _ in range(it):
        m = ((a + b) / 2).mid()
        sm = fun(m) - target
        if sm.contains(0):
            break
        if bool(sm < 0) == neg_left:
            a = m
        else:
            b = m
    sa, sb = fun(a) - target, fun(b) - target
    if not ((sa < 0 and sb > 0) or (sa > 0 and sb < 0)):
        raise bz.Fail("lost sign change")
    return a, b


# ----------------------------------------------------------------- CA.3
def accumulation():
    ctx.prec = 300
    a0, b0 = arb("-1.055"), arb("-1.054")
    pa, pb = psi_enc(a0), psi_enc(b0)
    print(f"   psi(-1.055) = {pa.str(12)},  psi(-1.054) = {pb.str(12)}")
    check(bool(pa < -1) and bool(pb > -1), "CA.3: psi(-1.055) < -1 < psi(-1.054) (arb), so (ψ(a0), ψ(b0)) is a neighbourhood of -1 in K_p")
    a, b = bracket_root(psi_enc, a0, b0, arb(-1))
    lam_star = -(a + b) / 2
    print(f"   lam* in [{(-b).str(20)}, {(-a).str(20)}]  (psi(-lam*) = -1)")
    check(bool(-b > arb(LAM_STAR[0])) and bool(-a < arb(LAM_STAR[1])), "CA.3: lam* = 1.0546444769654709... certified by a sign change of psi(-lam)+1")
    gp = lambda t: -arb(5) / 2 + (t + arb(25) / 4).sqrt()
    gm = lambda t: -arb(5) / 2 - (t + arb(25) / 4).sqrt()
    # |g_+-'| <= 1/sqrt5 on [-5,0] (BZ.1) and g_-(-4) = -4, g_+(-4) = -1: |h_n + 1| <= 5^(-1/2) |g_-^n(-3) + 4| <= 5^(-(n+1)/2)
    gap = min(-1 - pa, pb + 1)
    n0 = next(n for n in range(100) if bool(arb(5) ** (-arb(n + 1) / 2) < gap))
    w, ok_in, ok_rate, seen, hs = arb(-3), True, True, [], []
    for n in range(41):
        h = gp(w)
        ok_rate = ok_rate and bool(abs(h + 1) <= arb(5) ** (-arb(n + 1) / 2))
        if n >= n0:
            ok_in = ok_in and bool(h > pa) and bool(h < pb)
        seen.append(h)
        hs.append(h)
        w = gm(w)
    distinct = all(bool(seen[i] < seen[j]) or bool(seen[i] > seen[j]) for i in range(len(seen)) for j in range(i))
    check(ok_rate and ok_in and distinct,
          f"CA.3: h_n = g_+(g_-^n(-3)) in O are distinct, |h_n + 1| <= 5^(-(n+1)/2) (n <= 40), so h_n in (psi(a0), psi(b0)) for every n >= {n0}")
    # consistency: the first points of Lambda_O near lam* are zeros of f; n = 1, 2 are BZ's second and third largest atoms
    ok_f, ok_bz, cnt = True, True, 0
    for n in range(0, 8):
        t = float(hs[n].mid())
        u0 = -1.0546 if n >= 3 else {0: -0.7230513, 1: -1.1719102, 2: -1.0163164}[n]
        a, b = bracket_root(psi_enc, u0 - 5e-3 if n < 3 else -1.08, u0 + 5e-3 if n < 3 else -1.0, hs[n])
        mid = ((a + b) / 2).mid()
        fv = bz.inv_Q(mid)
        ok_f = ok_f and bool(abs(fv) < arb("1e-60"))
        if n in (1, 2):
            ok_bz = ok_bz and bool(abs(-mid - arb(BZ_LAMS[n - 1])) < arb("1e-14"))
        cnt += 1
        print(f"     n = {n}: h_n = {t:+.12f}, lam = {(-mid).str(22)},  f(-lam) = {fv.str(5)}")
    check(ok_f and ok_bz and cnt == 8, "CA.3 (consistency): f vanishes to 1e-60 at the points of Lambda_O for n = 0..7; n = 1, 2 are BZ's lam = 1.17191025..., 1.01631641...")
    a, b = bracket_root(psi_enc, -3.62, -3.61, arb(-3))
    check(bool(abs((a + b) / 2 + arb(bz.LAM_BIG[0])) < arb("1e-13")), "CA.2 (consistency): BZ's largest atom lam_1 = 3.61525671988770... has psi(-lam_1) = -3")


# ----------------------------------------------------------------- CA.5
def minus_two_witness():
    ctx.prec = 300
    a, b = bracket_root(psi_enc, -2.25, -2.23, arb(-2))
    ball = arb((a + b) / 2, ((b - a) / 2).upper())
    fv = inv_Q_reg(ball)
    print(f"   u with psi(u) = -2 in {ball.str(25)};  f(u) in {fv.str(16)}")
    check(bool(abs(ball - arb(U_TWO)) < arb("1e-28")) and bool(fv > arb(F_TWO[0])) and bool(fv < arb(F_TWO[1])),
          "Kept Failure CA.5 witness: psi(u) = -2 at u = -2.2421331852464538..., f(u) in [2.1979467595912, 2.1979467595913], not 0: no atom")
    # cross-check of the removable form against BZ's routine away from -2
    u = arb("-2.3")
    check(bool(abs(inv_Q_reg(u) - bz.inv_Q(u)) < arb("1e-40")), "CA.5 (consistency): removable form agrees with BZ's inv_Q at u = -2.3")


# ----------------------------------------------------------------- CA.4
def certify_with_box(lo, hi):
    """BZ's certify_atom (same steps), also returning the midpoint c, radius du and the holomorphy half-side rho."""
    a, b = arb(lo), arb(hi)
    for _ in range(110):
        m_ = ((a + b) / 2).mid()
        v = bz.inv_Q(m_)
        if v > 0:
            b = m_
        elif v < 0:
            a = m_
        else:
            break
    if not (bz.inv_Q(a) < 0 and bz.inv_Q(b) > 0):
        raise bz.Fail("no certified sign change")
    c = ((a + b) / 2).mid()
    du = ((b - a) / 2).upper()
    rho = (min(c - arb(lo), arb(hi) - c) / 4).mid()
    from flint import acb
    while True:
        if not rho > 4 * du:
            raise bz.Fail("box too small")
        try:
            M = abs(bz.inv_Q(acb(arb(c, rho.upper()), arb(0, rho.upper())), cplx=True)).upper()
            break
        except bz.Fail:
            rho = (rho / 8).mid()
    h = arb("1e-30")
    fd = (bz.inv_Q(c + h) - bz.inv_Q(c - h)) / (2 * h)
    err = h ** 2 * M / (rho - h) ** 3 + du * 2 * M / (rho - 2 * du) ** 2
    fp = fd + arb(0, err.upper())
    if not fp > 0:
        raise bz.Fail("derivative not certified positive")
    return c, du, rho, 1 / fp


def identify(c, du, rho):
    """Find k <= K_MAX with a sign change of psi(5^k u) + 3 on [c - s, c + s], s in {du, rho/2}."""
    five = arb(5)
    for k in range(K_MAX + 1):
        y = psi_enc(five ** k * c)
        if not bool(abs(y + 3) < arb("1e-6")):
            continue
        for s in (arb(du), (rho / 2).mid()):
            ya, yb = psi_enc(five ** k * (c - s)) + 3, psi_enc(five ** k * (c + s)) + 3
            if (ya < 0 and yb > 0) or (ya > 0 and yb < 0):
                return k
    return None


def law_for_bz_atoms():
    ctx.prec = 300
    g = arb(3).log() / arb(5).log()
    five = arb(5)
    t0 = time.time()
    us = [(-five ** (arb(i) / bz.N_SCAN)).mid() for i in range(bz.N_SCAN + 1)]
    vals = []
    for u in us:
        try:
            vals.append(bz.inv_Q(u))
        except bz.Fail:
            vals.append(None)
    cand = [i for i in range(bz.N_SCAN) if vals[i] is not None and vals[i + 1] is not None and vals[i] > 0 and vals[i + 1] < 0]
    print(f"   BZ.5 scan re-run: {len(cand)} candidate cells ({time.time() - t0:.0f}s)")
    found, ks, unid, tot = [], {}, 0, arb(0)
    for i in cand:
        try:
            c, du, rho, w = certify_with_box(us[i + 1], us[i])
        except bz.Fail:
            continue
        k = identify(c, du, rho)
        if k is None:
            unid += 1
        else:
            ks[k] = ks.get(k, 0) + 1
        lam = -(arb(c) + arb(0, du))
        found.append((lam, w, k))
        tot += w * lam ** (-g)
    mT = (arb(bz.C_R) + arb(0, arb("1e-59"))) / 2
    print(f"   certified atoms: {len(found)}; identified with Lambda_O: {len(found) - unid}; by k: {dict(sorted(ks.items()))}  ({time.time() - t0:.0f}s)")
    print(f"   sum lam^-gamma nu = {tot.str(14)}, m(T) - sum <= {float((mT - tot).upper()):.3e}")
    check(len(found) >= N_ATOMS_MIN and unid == 0,
          f"CA.4: all {len(found)} atoms of BZ.5 (>= {N_ATOMS_MIN}) satisfy psi(-5^k lam) = -3 for some 0 <= k <= {K_MAX} (unique zero of f in each holomorphy box)")


def main() -> int:
    t0 = time.time()
    print("Theorem CA: D = z/2 on the backward orbit of -3; every lam with psi(-5^k lam) = -3 is an atom of nu; infinitely many "
          "atoms per period, accumulating at lam* = 1.05464447696547...; all 311 BZ atoms obey the law; Kept Failure CA.5 (-2 branch). Not thrust.")
    symbolic_checks()
    accumulation()
    minus_two_witness()
    law_for_bz_atoms()
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CA checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
