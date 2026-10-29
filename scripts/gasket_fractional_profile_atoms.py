#!/usr/bin/env python3
r"""Theorem BZ: Q_inf continues meromorphically across (-inf,0) off a Lebesgue-null closed set Z, by an exact formula
1/Q_inf = Pi * (1/r_inf o psi) - S; so BX's measure nu has no absolutely continuous part, every singularity of Q_inf on
(-inf,0) \ Z is a genuine pole, and 311 atoms of nu are certified in one scale period, carrying all but 4.03e-5 of the
mass c_R(gamma)/2 of BX's periodic measure m; by almost periodicity, limsup rho_k >= 0.99968 and liminf rho_k <= 0.59497,
so Conjecture BV.6's asymptotic part is false (Kept Failure BZ.6) while its sharp-rate part holds (Kept Failure BZ.5:
nu is not absolutely continuous and not carried by <= 310 geometric pole series; extends Kept Failure BY.3 from 19).

Setting (GASKET.md, Theorems AQ, BD, BJ, BS, BX): p(z) = z(z+5), c(z) = (z+6)/((z+2)(z+5)), e(z) = 2z/(2+z),
r_inf = lim r_n, D = 1/r_inf, psi(v) = lim p^k(5^-k v) (BD.4), Q_inf(u) = lim (3/5)^M r_inf(5^-M u) (BS.1),
Q_inf(u) = int dnu(lam)/(lam+u) (BX.1), m = lam^-gamma nu on one period [1,5) in t = log_5 lam (BX.2), m(T) = c_R/2.
K_p = {z : p^n(z) bounded} (filled Julia set of p), Z = psi^-1(K_p).

Theorem BZ (checked here; exit 1 on any failure):
  BZ.1 [proved; sympy] p(z)+5/2 = (z+5/2)^2-15/4; |zeta|>5/2 escapes; the inverse branches g_+-(w) = +-sqrt(w+15/4)
       map the disc |w|<=5/2 into itself with |g'| <= 1/sqrt5; so K_p is a Cantor subset of [-5,0] covered by 2^n
       intervals of total length <= 5 (2/sqrt5)^n -> 0 (Lebesgue null), p^-1(K_p) = K_p.
  BZ.2 [proved; sympy] AQ's fixed point gives D(z) = c(z) D(p(z)) - e(z); D is holomorphic for |z| > 16 and
       |D(w)-w| <= 9 for |w| >= 144; hence r_inf is meromorphic on C \ K_p.
  BZ.3 [proved; sympy + arb] for every u in C \ (-inf,0]:
         1/Q_inf(u) = Pi(u) D(psi(u)) - S(u),  Pi(u) = prod_(j>=1) (5/3) c(psi(5^-j u)),
         S(u) = sum_(k>=1) (5/3)^k e(psi(5^-k u)) Pi(5^-k u);
       Pi, S are meromorphic on C, Z is a closed, Lebesgue-null subset of (-inf,0] with 5Z = Z, and Q_inf continues
       meromorphically to C \ Z.  Checks: the finite-level identity (5/3)^n/r_n(5^-n u) = Pi_n (3+psi_n) - S_n
       (AQ's recursion, exact up to arb radius); the formula reproduces BW's G* = F(t_max) and min F = F(t_min) to 1e-38.
  BZ.4 [proved from BZ.3] nu restricted to (0,inf) \ (-Z) is a sum of point masses with no accumulation point
       there, -Z is Lebesgue null, so nu and m have no absolutely continuous part (m is singular on the circle).
  BZ.5 [certified, arb/acb outward] scanning t = log_5 lam on a grid of 128000 cells of one period and certifying
       each candidate (finite acb enclosure of 1/Q_inf on a box around the zero = holomorphy, sign change at two
       real points, Cauchy-estimate enclosure of (1/Q_inf)' > 0) gives 311 distinct atoms lam_j in (1,5) of nu with
       nu({lam_j}) = 1/(1/Q_inf)'(-lam_j) enclosed; sum_j lam_j^-gamma nu({lam_j}) >= 0.25854 against
       m(T) = c_R(gamma)/2 = 0.25858291731525..., so the part of m not in these atoms is <= 4.1e-5.
  BZ.6 [proved from BZ.5; arb] m = m_pp + m_c with m_c(T) <= delta := m(T) - sum <= 4.03e-5; the coefficients a_k of
       m_pp are Bohr almost periodic in k (Dirichlet), |m_k - a_k| <= m_c(T); hence limsup rho_k >= 1 - 2 delta/m(T)
       >= 0.99968 and liminf rho_k <= (min_(k<=100) |a_k^cert| + delta)/m(T) <= 0.59497 (k = 35).
  Kept Failure BZ.6 (certified witness): Conjecture BV.6's "|theta_k| ~ const e^(-2 pi^2 k/ln5)", equivalently (BX.3)
       "rho_k converges to a nonzero limit", is false; its sharp-rate part holds: limsup |theta_k| e^(2pi^2k/ln5)
       >= (pi c_R/ln5)(1 - 2 delta/m(T)), so F continues to no wider strip.
  Kept Failure BZ.5 (certified witness): "nu is absolutely continuous" and "nu is carried by r <= 310 geometric
       series lam_j 5^Z" are false (311 certified atoms in one period; the largest at lam = 3.6152567198877...).
Observation (no claim): the uncaptured mass falls 6.9e-4, 1.8e-4, 4.0e-5 for grids of 4000, 16000, 128000 cells;
consistent with nu purely atomic, which is not proved (nu(-Z) = 0 is open).
Assumptions used as stated (not re-run): AQ (1)-(3) (fixed point r_inf(s) = r_inf(ps)/(c - e r_inf(ps)); r_n Stieltjes
with total mass 1 and spectrum in [0, 8] (Gershgorin, max degree 4); r_n up to r_inf; r_inf nonincreasing),
BD.4 (psi_k up to psi, psi(5v) = p(psi(v)), psi_k with nonnegative coefficients, psi(v) >= v), BJ.6 (psi(v) <= v e^(v/10)),
BS.1 (Q_inf = lim (3/5)^M r_inf(5^-M u)), BX.1-BX.2 (Stieltjes representation, m(T) = c_R(gamma)/2), BW's digits of
G*, min F, t_max, t_min and c_R(gamma) (consistency and mass comparison only).
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time

import sympy as sp
from flint import acb, arb, ctx

FAILS: list[str] = []
N_SCAN = 128000
J_TAIL = 200
W_ESC = 60            # escape radius 10^60
D_FAR = 9             # |D(w) - w| <= 9 for |w| >= 144
N_ATOMS_MIN = 311
MASS_FLOOR = "0.25854"
K_AP = 100
SUP_FLOOR = "0.9996"
INF_CEIL = "0.5950"
BY_RHO = {1: "0.769641861455", 2: "0.600852173802", 3: "0.833706407581", 4: "0.855147103961", 5: "0.721589627521"}
LAM_BIG = ("3.61525671988770", "3.61525671988771")
# BW (Theorem BW.2/BW.3) digits, used for consistency only
C_R = "0.517165834630505243925980485673127200453783829877158701845373"
G_STAR = "0.600961578597101089042056184134616926770040525917258197021856877160250655142022957500717"
MIN_F = "0.600946923794265612104226188172129248411629691787677697008509741029502767706538876230717"
T_MAX = "0.7517559617903175587780293327382082565622019261165521457856"
T_MIN = "0.2517561241922139288493028770006731010657456148245398158230"


class Fail(Exception):
    pass


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- symbolic
def symbolic_checks():
    z, R, w, x, a = sp.symbols("z R w x a")
    c = (z + 6) / ((z + 2) * (z + 5))
    e = 2 * z / (2 + z)
    pz = z * (z + 5)
    # AQ fixed point r = R/(c - e R), R = r(p z)  <=>  1/r = c/R - e
    check(sp.simplify(1 / (R / (c - e * R)) - (c / R - e)) == 0, "BZ.2: AQ's fixed point is D(z) = c(z) D(p(z)) - e(z), D = 1/r_inf")
    check(sp.simplify(c - sp.Rational(3, 5) * (1 + z / 6) / ((1 + z / 2) * (1 + z / 5))) == 0,
          "BZ.3: (5/3)c(x) = (1+x/6)/((1+x/2)(1+x/5)) -> 1 (tail factor)")
    # unrolled identity, M = 3, generic b_j, e_j:  X_k = b_k X_(k-1) - (5/3)^k e_k
    b = sp.symbols("b1:4")
    ee = sp.symbols("e1:4")
    X0 = sp.Symbol("X0")
    X = X0
    for k in range(1, 4):
        X = b[k - 1] * X - sp.Rational(5, 3) ** k * ee[k - 1]
    closed = b[0] * b[1] * b[2] * X0 - sum(sp.Rational(5, 3) ** k * ee[k - 1] * sp.prod(b[k:3]) for k in range(1, 4))
    check(sp.expand(X - closed) == 0, "BZ.3: unrolled recursion (5/3)^M D(x_M) = prod b_j D(x_0) - sum (5/3)^k e_k prod_(j>k) b_j (M=3)")
    # conjugacy, fixed point, escape, branch derivative bound
    check(sp.expand(pz + sp.Rational(5, 2) - ((z + sp.Rational(5, 2)) ** 2 - sp.Rational(15, 4))) == 0,
          "BZ.1: p(z)+5/2 = (z+5/2)^2 - 15/4 (c = -15/4 < -2), fixed point beta = 5/2")
    check(sp.expand((x - sp.Rational(5, 2)) * (x + sp.Rational(3, 2)) - (x ** 2 - sp.Rational(15, 4) - x)) == 0,
          "BZ.1: |f(zeta)| - |zeta| >= |zeta|^2 - 15/4 - |zeta| = (|zeta|-5/2)(|zeta|+3/2) > 0 for |zeta| > 5/2")
    check(sp.Rational(15, 4) - sp.Rational(5, 2) == sp.Rational(5, 4) and sp.Rational(1, 2) / sp.sqrt(sp.Rational(5, 4)) == 1 / sp.sqrt(5),
          "BZ.1: |w + 15/4| >= 5/4 on |w| <= 5/2, so |g_+-'| = 1/(2 sqrt|w+15/4|) <= 1/sqrt5 < 1; 2/sqrt5 < 1")
    check(sp.sqrt(sp.Rational(5, 4)) + 0 <= sp.Rational(5, 2) and sp.sqrt(sp.Rational(25, 4)) == sp.Rational(5, 2),
          "BZ.1: g_+-([-5/2,5/2]) = +-[sqrt5/2, 5/2] lie in [-5/2,5/2] (real branches)")
    # psi(v) = v + v^2/20 + ...: psi(5v) = p(psi(v)) at order 2
    eq = sp.expand((5 * w + a * 25 * w ** 2) - ((w + a * w ** 2) * (w + a * w ** 2 + 5))).coeff(w, 2)
    check(sp.solve(eq, a) == [sp.Rational(1, 20)], "BZ.3: psi(v) = v + v^2/20 + O(v^3) (only used for the float scan)")
    # far field of D: r = (1/w) int dmu/(1+lam/w), lam in [0,8], |w| >= 144: |D - w| <= 8|w|/(|w|-16) <= 9
    W = sp.Integer(144)
    check(8 * W / (W - 16) <= D_FAR, "BZ.2: |D(w) - w| <= 8|w|/(|w|-16) <= 9 for |w| >= 144 (mu probability on [0,8])")
    # tail constants: |log b| <= 1.06 (1/6+1/2+1/5)|x| <= |x| for |x| <= 0.1 (|y| <= 0.05); |e(x)| <= 2|x|/1.9
    check(sp.Rational(106, 100) * (sp.Rational(1, 6) + sp.Rational(1, 2) + sp.Rational(1, 5)) < 1
          and 1 / (1 - sp.Rational(5, 100)) < sp.Rational(106, 100), "BZ.3: tail constants |log (5/3)c(x)| <= |x|, |e(x)| <= 1.06|x| for |x| <= 0.1")


# ----------------------------------------------------------------- the formula in interval arithmetic
F3 = None


def _p(z):
    return z * (z + 5)


def _c(z):
    return (z + 6) / ((z + 2) * (z + 5))


def _e(z):
    return 2 * z / (2 + z)


def _nz(z):
    if z.contains(0):
        raise Fail("denominator encloses 0")


def D_of(x, cplx):
    """Enclosure of D(x) = 1/r_inf(x): iterate p until |p^N x| > 10^60, use |D(w)-w| <= 9 there (BZ.2), recurse back."""
    Wb = arb(10) ** W_ESC
    orb = [x]
    while not (abs(orb[-1]).lower() > Wb if cplx else orb[-1] > Wb):
        orb.append(_p(orb[-1]))
        if len(orb) > 300:
            raise Fail("orbit does not escape")
    y = orb[-1]
    d = y + (acb(arb(0, D_FAR), arb(0, D_FAR)) if cplx else arb(0, D_FAR))
    for zz in reversed(orb[:-1]):
        _nz(zz + 2)
        _nz(zz + 5)
        d = _c(zz) * d - _e(zz)
    return d


def inv_Q(u, cplx=False, J=J_TAIL):
    """Outward enclosure of 1/Q_inf(u) = Pi(u) D(psi(u)) - S(u) (BZ.3) for a real or complex ball u."""
    five = arb(5)
    f3 = arb(5) / 3
    ball = (lambda r: acb(arb(0, r), arb(0, r))) if cplx else (lambda r: arb(0, r))
    v = u / five ** J
    av = arb(abs(v).upper())
    dl = (av * ((av / 10).exp() - 1)).upper()        # |psi(v)-v| <= psi(|v|)-|v| <= |v|(e^(|v|/10)-1)  (BD.4, BJ.6)
    x = [None] * (J + 1)
    x[J] = v + ball(dl)
    for j in range(J, 0, -1):
        x[j - 1] = _p(x[j])                          # x_j = psi(5^-j u), psi(5v) = p(psi(v))
    for j in range(1, J + 1):
        _nz(x[j] + 2)
        _nz(x[j] + 5)
    b = [None] + [f3 * _c(x[j]) for j in range(1, J + 1)]
    eps = (abs(u).upper() * five ** (-J)).upper()     # sum_(j>J) |x_j| <= |u| 5^-J / 2
    tail = (acb(arb(0, eps)) if cplx else arb(0, eps)).exp()
    Pi = acb(1) if cplx else arb(1)
    for j in range(1, J + 1):
        Pi *= b[j]
    S = acb(0) if cplx else arb(0)
    tl = tail
    for k in range(J, 0, -1):
        S += f3 ** k * _e(x[k]) * tl
        tl *= b[k]
    S += ball((arb(42) / 10 * abs(u).upper() * arb(3) ** (-J)).upper())   # sum_(k>J) <= 1.1 |u| 3^-J
    return Pi * tail * D_of(x[0], cplx) - S


def finite_level_check():
    """(5/3)^n / r_n(5^-n u) = Pi_n (3 + psi_n(u)) - S_n with AQ's r_n: the same unrolling with D_0 = 3 + z."""
    ctx.prec = 400
    f3 = arb(5) / 3
    ok = True
    for n, u in ((12, arb("1.7")), (20, arb("0.3")), (9, arb("-0.37"))):
        def r(k, s):
            if k == 0:
                return 1 / (3 + s)
            R = r(k - 1, _p(s))
            return (2 + s) * (5 + s) * R / ((6 + s) - 2 * s * (5 + s) * R)
        xn = u / arb(5) ** n
        lhs = f3 ** n / r(n, xn)
        x = [None] * (n + 1)
        x[n] = xn
        for j in range(n, 0, -1):
            x[j - 1] = _p(x[j])
        bb = [None] + [f3 * _c(x[j]) for j in range(1, n + 1)]
        P = arb(1)
        for j in range(1, n + 1):
            P *= bb[j]
        S, tl = arb(0), arb(1)
        for k in range(n, 0, -1):
            S += f3 ** k * _e(x[k]) * tl
            tl *= bb[k]
        rhs = P * (3 + x[0]) - S
        ok = ok and bool(abs(lhs - rhs) < arb("1e-80") * (1 + abs(lhs)))
    check(ok, "BZ.3: finite-level identity (5/3)^n/r_n(5^-n u) = Pi_n(3+psi_n) - S_n (AQ recursion; n=12, 20, 9; u=1.7, 0.3, -0.37)")


def profile_checks():
    ctx.prec = 400
    g = arb(3).log() / arb(5).log()
    ok = True
    for t, ref in ((T_MAX, G_STAR), (T_MIN, MIN_F)):
        u = arb(5) ** arb(t)
        Fv = u ** (1 - g) / inv_Q(u)
        print(f"   F(5^{t[:12]}...) = {Fv.str(40)}")
        ok = ok and bool(abs(Fv - arb(ref)) < arb("1e-38"))
    check(ok, "BZ.3: the formula gives F(t_max) = G* and F(t_min) = min F of Theorem BW to 1e-38 (consistency)")


# ----------------------------------------------------------------- atoms
def certify_atom(lo, hi):
    """lo < hi < 0 with 1/Q(lo) < 0 < 1/Q(hi) (floats from the scan).  Returns (lam ball, nu({lam}) ball) or raises."""
    a, b = arb(lo), arb(hi)
    for _ in range(110):
        m_ = ((a + b) / 2).mid()
        v = inv_Q(m_)
        if v > 0:
            b = m_
        elif v < 0:
            a = m_
        else:
            break
    if not (inv_Q(a) < 0 and inv_Q(b) > 0):
        raise Fail("no certified sign change")
    c = ((a + b) / 2).mid()
    du = ((b - a) / 2).upper()
    rho = (min(c - arb(lo), arb(hi) - c) / 4).mid()
    while True:
        if not rho > 4 * du:
            raise Fail("box too small")
        try:
            M = abs(inv_Q(acb(arb(c, rho.upper()), arb(0, rho.upper())), cplx=True)).upper()
            break
        except Fail:
            rho = (rho / 8).mid()
    # f = 1/Q holomorphic on the box (finite enclosure) -> real-analytic on [c-rho, c+rho], zero u* in (a, b).
    h = arb("1e-30")
    fd = (inv_Q(c + h) - inv_Q(c - h)) / (2 * h)
    err = h ** 2 * M / (rho - h) ** 3 + du * 2 * M / (rho - 2 * du) ** 2     # Cauchy estimates for f''' and f''
    fp = fd + arb(0, err.upper())
    if not fp > 0:
        raise Fail("derivative not certified positive")
    return -(arb(c) + arb(0, du)), 1 / fp


def atoms():
    ctx.prec = 300
    g = arb(3).log() / arb(5).log()
    five = arb(5)
    t0 = time.time()
    us = [(-five ** (arb(i) / N_SCAN)).mid() for i in range(N_SCAN + 1)]
    vals = []
    for u in us:
        try:
            vals.append(inv_Q(u))
        except Fail:
            vals.append(None)
    cand = [i for i in range(N_SCAN) if vals[i] is not None and vals[i + 1] is not None and vals[i] > 0 and vals[i + 1] < 0]
    print(f"   scan of {N_SCAN} cells of t = log_5 lam in [0,1]: {len(cand)} candidate sign changes ({time.time() - t0:.0f}s)")
    found, tot, fails, big = [], arb(0), 0, None
    for i in cand:
        try:
            lam, w = certify_atom(us[i + 1], us[i])
        except Fail:
            fails += 1
            continue
        found.append((lam, w))
        tot += w * lam ** (-g)
        if big is None or w.mid() > big[1].mid():
            big = (lam, w)
    print(f"   certified atoms: {len(found)} (uncertified candidates {fails}); {time.time() - t0:.0f}s")
    for lam, w in sorted(found, key=lambda q: -float(q[1].mid()))[:8]:
        print(f"     lam = {lam.str(16)}   nu({{lam}}) = {w.str(12)}   lam^-gamma nu = {(w * lam ** (-g)).str(10)}")
    mT = (arb(C_R) + arb(0, arb("1e-59"))) / 2
    gap = mT - tot
    print(f"   sum lam^-gamma nu over certified atoms = {tot.str(14)};  m(T) = c_R/2 = {mT.str(16)};  rest <= {float(gap.upper()):.3e}")
    inside = all(bool(l > 1) and bool(l < 5) for l, _ in found)
    distinct = all(bool(found[k][0] < found[k + 1][0]) or bool(found[k][0] > found[k + 1][0]) for k in range(len(found) - 1))
    check(inside and distinct and len(found) >= N_ATOMS_MIN, f"BZ.5: {len(found)} >= {N_ATOMS_MIN} distinct atoms of nu certified in (1,5) (one period)")
    check(bool(tot.lower() > arb(MASS_FLOOR)) and bool(tot < mT), f"BZ.5: sum of lam^-gamma nu({{lam}}) over them in [{MASS_FLOOR}, m(T)), m(T)-sum <= 4.1e-5")
    check(bool(gap.upper() < arb("4.1e-5")), "BZ.5: part of m outside the certified atoms <= 4.1e-5 (0.016% of m(T))")
    bv6(found, mT, gap)
    check(big is not None and bool(big[0] > arb(LAM_BIG[0])) and bool(big[0] < arb(LAM_BIG[1])) and bool(big[1] > 0),
          f"Kept Failure BZ.5 witness: atom at lam = {big[0].str(15) if big else '?'} with nu({{lam}}) = {big[1].str(10) if big else '?'} > 0")


def bv6(found, mT, gap):
    """BZ.6: almost periodicity of the atomic part's coefficients against Conjecture BV.6."""
    two_pi = 2 * arb.pi()
    l5 = arb(5).log()
    ts = [(lam.log() / l5, w * lam ** (-(arb(3).log() / l5))) for lam, w in found]
    delta = arb(gap.upper())
    cons, amin, kmin = True, None, None
    for k in range(1, K_AP + 1):
        a = acb(0)
        for t, w in ts:
            a += acb(0, -two_pi * k * t).exp() * w
        if k in BY_RHO:
            cons = cons and bool(abs(abs(a) / mT - arb(BY_RHO[k])) < delta / mT + arb("1e-11"))
        if amin is None or abs(a).upper() < amin:
            amin, kmin = abs(a).upper(), k
    sup_lo = 1 - 2 * delta / mT
    inf_hi = (arb(amin) + delta) / mT
    print(f"   |a_k|/m(T) of the certified atoms agree with BY's rho_1..rho_5 within delta/m(T) = {float((delta / mT).upper()):.2e}")
    print(f"   BZ.6: limsup rho_k >= {float(sup_lo.lower()):.8f}  "
          f"liminf rho_k <= {float(inf_hi.upper()):.8f} (from k = {kmin})")
    check(cons, "consistency: |a_k| of the certified atoms is within delta of BY's certified rho_k m(T), k = 1..5 (normalisation of m)")
    check(bool(sup_lo > arb(SUP_FLOOR)) and bool(inf_hi < arb(INF_CEIL)),
          f"BZ.6: limsup rho_k >= {SUP_FLOOR} and liminf rho_k <= {INF_CEIL}, so rho_k has no limit: Conjecture BV.6 is false")


def main() -> int:
    t0 = time.time()
    print("Theorem BZ: 1/Q_inf = Pi (1/r_inf o psi) - S, Q_inf meromorphic off a null set Z; nu has no a.c. part; "
          "311 certified atoms carry all but 4.03e-5 of m; limsup rho_k >= 0.9996 > 0.5950 >= liminf rho_k. Kept Failures BZ.5, BZ.6. Not thrust.")
    symbolic_checks()
    finite_level_check()
    profile_checks()
    atoms()
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BZ checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
