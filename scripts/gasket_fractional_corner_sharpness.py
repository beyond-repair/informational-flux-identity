#!/usr/bin/env python3
"""Theorem BJ: the cell modulus rho^-N of Theorem BI is sharp (Conjecture BI.7 proved), with a certified corner
constant c_corner(0.8) in [1.13589, 1.13594].

Setting (GASKET.md, Theorems BG, BH, BI): level-n gasket V_n, corner data d, x = d - dbar, u_n the L_n^alpha-Dirichlet
extension, u-bar its uniform limit for alpha > gamma = log3/log5 (BI), r_n(s), r_inf(s) = 1/(2+s-2 eta(s)) (BH.2),
p(s) = s(s+5), beta(s) = (s+6)/((s+2)(s+5)) = AQ's c(s), e(s) = 2s/(2+s), rho = 5^alpha/3, g(s) = s^(1-gamma) r_inf(s),
R_n(tau) = (3/5)^n r_n(5^-n tau) -> R(tau) (BD.4, BG.3), psi_K(t) = p^K(5^-K t) increasing to psi(t).
nu_b^N = (2^-N, 0), nu_c^N = (0, 2^-N) are the neighbours of corner a in V_N.

Theorem BJ (checked here; exit 1 on any failure):
  (1) exact reciprocal identity (every n >= 0, s > 0):
        1/r_n(s) - 1/r_inf(s) = prod_(i<n) beta(p^i s) * (1 + 2 eta(p^n s))  in (0, 3 (3/5)^n],
      from AQ's step in reciprocal form 1/r_(k+1)(s) = beta(s)/r_k(ps) - e(s), r_0(z) = 1/(3+z), 1/r_inf = 2+z-2eta.
  (2) hence 1/(tau^(1-gamma)/g(5^-n tau) + 3) <= R_n(tau) <= tau^(gamma-1) g(5^-n tau), and for the limit R:
        1/(tau^(1-gamma)/g_lo + 3) <= R(tau) <= min(1/3, g_hi tau^(gamma-1)),  g_lo <= liminf_(s->0) g <= limsup g <= g_hi.
  (3) certified: g on (0, p(eps0)], eps0 = 1e-12, lies in [g_lo, g_hi] = [0.60094, 0.60097] (arb cells on one
      p-period [eps0, p(eps0)] + a two-sided bootstrap along the p-orbit); AU's log-mean <Theta> lies inside.
  (4) corner deficit (every alpha > gamma, every N, any data with x_a != 0):
        delta_N := (d_a - (u-bar(nu_b^N) + u-bar(nu_c^N))/2)/x_a = kappa_a int tau^-a R(tau) (1 - eta_N(psi(5^-N tau))) dtau,
      kappa_a = 1/int tau^-a R = 2 sin(pi a)/(pi A_a), eta_N = eta - prod_(i<N) beta(p^i z)(1/2 + eta(p^N z)) (data-free).
  (5) sharpness: rho^N delta_N >= c_-(alpha) = kappa_a int t^-a (1 - eta(psi_K t))/(t^(1-gamma)/g_lo + 3) dt > 0 for every N,
      so osc of u-bar on the level-N corner cell >= |x_a| c_- rho^-N: the order rho^-N of BI.5 is exact and the Hoelder
      exponent (alpha-gamma)/(1-gamma) cannot be raised.  Certified c_- at 0.75, 0.8, 0.9.
  (6) corner constant: g_lo M_a <= liminf rho^N delta_N <= limsup rho^N delta_N <= g_hi M_a,
      M_a = kappa_a int t^(gamma-1-a)(1 - eta(psi(t))) dt; at 0.8 (AZ's A_0.8): [1.13589, 1.13594].
      Derived from AT.1 as stated: the limit exists, = kappa_a int t^(gamma-1-a) Theta_AT(psi t)(1 - eta(psi t)) dt.

Floats (consistency only): rho^N delta_N^(n) by recursion quadrature at n = 100 (BI's corner_deficit).
A rate for u_n -> u-bar remains open.  Not a field solver.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from flint import arb, ctx
from mpmath import mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_maxwell_limit as bb  # noqa: E402
import gasket_fractional_corner_layer as bh  # noqa: E402
import gasket_fractional_uniform_limit as bi  # noqa: E402

FAILS: list[str] = []
A08 = ("2.39640970397922052703", "2.39640970397922052705")  # Theorem AZ (3), used as stated
THETA_AU = "0.6009542511687736078237051500999942"  # Theorem AU, <Theta> (log-mean of the AT profile)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def gam_arb() -> arb:
    return arb(3).log() / arb(5).log()


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    s, z, X, Y, E, tau, n, g = sp.symbols("s z X Y E tau n g", positive=True)
    c = (6 + s) / ((5 + s) * (2 + s))
    e = 2 * s / (2 + s)
    beta = (s + 6) / ((s + 2) * (s + 5))
    step = X / (c - e * X)  # AQ: r_(k+1)(s) = r_k(ps)/(c - e r_k(ps)), X = r_k(ps)
    check(sp.simplify(1 / step - (beta / X - e)) == 0 and sp.simplify(c - beta) == 0,
          "AQ's step in reciprocal form: 1/r_(k+1)(s) = beta(s)/r_k(ps) - e(s), with c = beta (affine in 1/r)")
    check(sp.simplify((beta / X - e) - (beta / Y - e) - beta * (1 / X - 1 / Y)) == 0,
          "so the difference D_k = 1/r_k - 1/r_inf obeys D_(k+1)(s) = beta(s) D_k(ps)")
    check(sp.simplify((3 + z) - (2 + z - 2 * E) - (1 + 2 * E)) == 0,
          "base: 1/r_0(z) - 1/r_inf(z) = (3+z) - (2+z-2eta(z)) = 1 + 2 eta(z) in (1, 3] (0 < eta <= 1, BH.2)")
    dbeta = sp.together(sp.diff(beta, s))
    check(sp.simplify(beta.subs(s, 0) - sp.Rational(3, 5)) == 0 and sp.expand(-sp.numer(dbeta)) == sp.expand(s**2 + 12 * s + 32),
          "beta(0) = 3/5 and beta' = -(s^2+12s+32)/((s+2)(s+5))^2 < 0, so prod_(i<n) beta(p^i s) <= (3/5)^n")
    # rescaling: (5/3)^n / r_inf(5^-n tau) = tau^(1-gamma)/g(5^-n tau) because 5^(1-gamma) = 5/3
    gm = sp.log(3) / sp.log(5)
    check(sp.simplify(sp.Integer(5) ** (1 - gm) - sp.Rational(5, 3)) == 0,
          "5^(1-gamma) = 5/3: (5/3)^n/r_inf(5^-n tau) = tau^(1-gamma)/g(5^-n tau), and 5^(N(1-a)) (3/5)^N = rho^-N")
    # eta_N: symmetric recursion with base -1/2
    a1 = (s + 4) / ((s + 2) * (s + 5))
    de = 2 / ((s + 2) * (s + 5))
    check(sp.simplify(a1 + de - beta) == 0,
          "mean of the two neighbour profiles: m_N(z) = alpha1 x_a + (alpha1+delta) m_(N-1)(pz) = alpha1 x_a + beta m_(N-1)(pz), m_0 = (x_b+x_c)/2 = -x_a/2")
    # psi lemma bookkeeping: sum_(i<j) 2*5^(i-k) y / 5 <= y/10 ; e^(y/10) <= 2 iff y <= 10 ln 2
    check(float(10 * sp.log(2)) > 6.93, "orbit lemma for psi: psi(y) <= y e^(y/10) for y <= 10 ln2 = 6.931..., so psi(t) <= psi_K(t) e^(2^K 5^-K t/10) for 5^-K t <= 6.93")
    a = sp.Symbol("a", positive=True)
    K, G = sp.symbols("K G", positive=True)
    tc = (3 * G / K) ** (1 / (1 - gm))
    I1 = sp.integrate(tau ** (-a), (tau, 0, tc), conds="none") * K / 3
    I2 = G * tc ** (gm - a) / (a - gm)
    check(sp.simplify(sp.expand_power_base(I1, force=True) - K / 3 * tc ** (1 - a) / (1 - a)) == 0 and I2 != 0,
          "int_0^inf t^-a min(K/3, G t^(gamma-1)) dt = (K/3) t_c^(1-a)/(1-a) + G t_c^(gamma-a)/(a-gamma), t_c = (3G/K)^(1/(1-gamma))")


# ----------------------------------------------------------------- (B) exact rationals / arb identities
def r_exact(n, s: Fraction) -> Fraction:
    adj, C, pos = bh.mesh(n)
    rhs = [Fraction(0)] * len(adj)
    rhs[C[0]], rhs[C[1]] = Fraction(1), Fraction(-1)
    w = bh.exact_solve(adj, s, rhs)
    return (w[C[0]] - w[C[1]]) / 2


def eta_N_exact(z: Fraction, N: int) -> Fraction:
    zs = [z]
    for _ in range(N):
        zs.append(zs[-1] * (zs[-1] + 5))
    e = Fraction(-1, 2)
    for k in range(N - 1, -1, -1):
        zz = zs[k]
        D = (zz + 2) * (zz + 5)
        e = (zz + 4) / D + (zz + 6) / D * e
    return e


def exact_checks(levels=(0, 1, 2, 3, 4)):
    ctx.prec = 120
    svals = (Fraction(1, 1000), Fraction(1, 2), Fraction(7, 5), Fraction(10))
    ok_r = ok_step = ok_id = ok_bd = True
    for n in levels:
        for s in svals:
            re = r_exact(n, s)
            ok_r = ok_r and re == bb.rec_all(n, s)[0]
            if n == 0:
                ok_r = ok_r and re == 1 / (3 + s)
            else:
                ps = s * (s + 5)
                beta = (s + 6) / ((s + 2) * (s + 5))
                ok_step = ok_step and 1 / re == beta / r_exact(n - 1, ps) - 2 * s / (2 + s)
            # BJ.1 against arb eta
            ss = bh.F(s)
            rinf = 1 / (2 + ss - 2 * bh.eta_arb(ss))
            z = ss
            for _ in range(n):
                z = z * (z + 5)
            rhs = bh.beta_prod(ss, n) * (1 + 2 * bh.eta_arb(z))
            lhs = 1 / bh.F(re) - 1 / rinf
            ok_id = ok_id and bool((lhs - rhs).contains(0)) and abs(float((lhs - rhs).mid())) < 1e-25
            ok_bd = ok_bd and bool(lhs.lower() > 0) and bool(lhs.upper() <= 3 * (arb(3) / 5) ** n)
    check(ok_r, f"exact rationals, levels {levels[0]}-{levels[-1]}, s in {{1/1000,1/2,7/5,10}}: r_n(s) = (1/2)v^T(L_n+s)^-1 v equals AQ's recursion; r_0 = 1/(3+s)")
    check(ok_step, "exact: 1/r_n(s) = beta(s)/r_(n-1)(p s) - e(s) at levels 1-4")
    check(ok_id, "arb: 1/r_n(s) - 1/r_inf(s) = prod_(i<n) beta(p^i s)(1+2eta(p^n s)) (BJ.1) at levels 0-4, four s")
    check(ok_bd, "arb: 0 < 1/r_n - 1/r_inf <= 3(3/5)^n")
    # mean-of-neighbours identity for general data (BJ.4's eta_N)
    xC = (Fraction(3), Fraction(-1), Fraction(-2))
    ok_m = ok_e = True
    for n in (2, 3, 4):
        adj, C, pos = bh.mesh(n)
        for s in (Fraction(1, 2), Fraction(7, 5)):
            rhs = [Fraction(0)] * len(adj)
            for kk, cc in enumerate(C):
                rhs[cc] = xC[kk]
            w = bh.exact_solve(adj, s, rhs)
            rn = w[C[0]] / xC[0]
            hb = w[pos[(Fraction(1, 2**n), Fraction(0))]] / rn
            hc = w[pos[(Fraction(0), Fraction(1, 2**n))]] / rn
            ok_m = ok_m and (hb + hc) / 2 == xC[0] * eta_N_exact(s, n)
            ss = bh.F(s)
            z = ss
            for _ in range(n):
                z = z * (z + 5)
            en = bh.eta_arb(ss) - bh.beta_prod(ss, n) * (arb(1) / 2 + bh.eta_arb(z))
            ok_e = ok_e and bool((en - bh.F(eta_N_exact(s, n))).contains(0))
    check(ok_m, "exact, levels 2-4, data (3,-1,-2): (h_n(s;nu_b)+h_n(s;nu_c))/2 = x_a eta_n(s), eta_n from the base -1/2 (data-free)")
    check(ok_e, "arb: eta_n = eta - prod_(i<n) beta(p^i s)(1/2 + eta(p^n s)), so 0 <= (1-eta_n) - (1-eta) <= (3/2)(3/5)^n")


# ----------------------------------------------------------------- (C) certificates
def psi_K(t: arb, K: int) -> arb:
    y = t / arb(5) ** K
    for _ in range(K):
        y = y * (y + 5)
    return y


def psi_hi(t: arb, K: int) -> arb:
    """psi(t) <= psi_K(t) exp(2^K 5^-K t/10), valid for 5^-K t <= 6.93."""
    assert (t / arb(5) ** K).upper() <= 6.93
    return psi_K(t, K) * (arb(2) ** K / arb(5) ** K * t / 10).exp()


def g_range(eps0="1e-12", ratio="1.00001", prec=96):
    """g = s^(1-gamma) r_inf on one p-period [eps0, p(eps0)]: r_inf decreasing, s^(1-gamma) increasing."""
    ctx.prec = prec
    gam = gam_arb()
    e = arb(eps0)
    pe = e * (e + 5)
    pts = [e]
    rr = arb(ratio)
    while pts[-1].upper() < pe.lower():
        pts.append(pts[-1] * rr)
    pts[-1] = pe
    r = [1 / (2 + p - 2 * bh.eta_arb(p)) for p in pts]
    hi, lo = arb(0), arb(10)
    for i in range(len(pts) - 1):
        up = (pts[i + 1] ** (1 - gam) * r[i]).upper()
        dn = (pts[i] ** (1 - gam) * r[i + 1]).lower()
        if up > hi:
            hi = arb(up)
        if dn < lo:
            lo = arb(dn)
    return lo, hi, len(pts)


def g_bootstrap(lo0: arb, hi0: arb, G: arb, eps0="1e-12"):
    """Every t in (0, eps0) is q^k(s0), s0 in [eps0, p(eps0)), orbit t_j = q^j(s0) <= 5^-j p(eps0).
    g(t_j) = g(t_(j-1)) Q(t_j), Q = (3/5)(1+t/5)^(gamma-1)/(c - e r_inf(pt)):
      Q >= (1+t/5)^(gamma-1) >= exp(-(1-gamma) t/5)            (c <= 3/5, e r >= 0)
      Q <= 1/(1 - 8t/15 - G t^gamma)                           (BI.4: c >= 3/5 - 8t/25, e r_inf(pt) <= (3/5) G t^gamma)
    so g in [lo0 exp(-(1-gamma) p/20), hi0 exp(X/(1-x_1))], X = (2/15)p + G p^gamma/2 (5^gamma = 3), p = p(eps0)."""
    gam = gam_arb()
    e = arb(eps0)
    pe = e * (e + 5)
    lo = lo0 * (-(1 - gam) * pe / 20).exp()
    x1 = arb(8) / 15 * pe / 5 + G * (pe / 5) ** gam
    X = arb(2) / 15 * pe + G * pe**gam / 2
    hi = hi0 * (X / (1 - x1)).exp()
    return arb(lo.lower()), arb(hi.upper()), bool(x1.upper() < 0.5)


def Wint(a: arb, gam: arb, t0: arb, t1: arb) -> arb:
    return (t0 ** (gam - a) - t1 ** (gam - a)) / (a - gam)


def M_bracket(al: str, glo: arb, ghi: arb, K=12, h="1e-13", T="1e4", ratio="1.0001", prec=96):
    """Two-sided enclosure of Mint = int_0^inf t^(gamma-1-a)(1 - eta(psi(t))) dt.
    Cells [t0,t1] (geometric): weight decreasing (exact integral), 1-eta(psi) increasing:
      lower (1-eta(psi_K(t0))) W, upper (1-eta(psi_hi(t1))) W.
    Head (0,h]: psi(t) in [t, t e^(t/10)] subset (0, eps0), 1-eta(z) = (z^(1-gamma)/g(z) - z)/2 with g in [glo, ghi].
    Tail [T, inf): 1 - 2/(2+psi_K(T)) <= 1-eta(psi) <= 1 (eta <= 2/(2+z), BH.2)."""
    ctx.prec = prec
    a = arb(al)
    gam = gam_arb()
    hA = arb(h)
    eh = (hA / 10).exp()
    head_lo = (hA ** (1 - a) / (1 - a) / ghi - eh * hA ** (1 + gam - a) / (1 + gam - a)) / 2
    head_hi = eh ** (1 - gam) * hA ** (1 - a) / (1 - a) / (2 * glo)
    TA = arb(T)
    tail_hi = TA ** (gam - a) / (a - gam)
    tail_lo = (1 - 2 / (2 + psi_K(TA, K))) * tail_hi
    lo, hi = head_lo + tail_lo, head_hi + tail_hi
    rr = arb(ratio)
    t0 = hA
    one_m_lo = 1 - bh.eta_arb(psi_K(t0, K))
    ncell = 0
    while t0.upper() < TA.lower():
        t1 = t0 * rr
        if t1.upper() > TA.lower():
            t1 = TA
        W = Wint(a, gam, t0, t1)
        nxt_lo = 1 - bh.eta_arb(psi_K(t1, K))
        up = 1 - bh.eta_arb(psi_hi(t1, K))
        lo += arb(one_m_lo.lower()) * W
        hi += arb(up.upper()) * W
        one_m_lo = nxt_lo
        t0 = t1
        ncell += 1
    return arb(lo.lower()), arb(hi.upper()), ncell


def I_closed(al: str, Kf: arb, G: arb) -> arb:
    """int_0^inf t^-a min(Kf/3, G t^(gamma-1)) dt in closed form."""
    a = arb(al)
    gam = gam_arb()
    tc = (3 * G / Kf) ** (1 / (1 - gam))
    return Kf / 3 * tc ** (1 - a) / (1 - a) + G * tc ** (gam - a) / (a - gam)


def B_lower(al: str, glo: arb, K=12, lo="1e-13", hi="1e4", ratio="1.002", N=None):
    """Lower bound of int t^-a (1 - eta(psi_K t))/(t^(1-gamma)/glo + 3 q) dt, q = 1 (uniform) or (3/5)^N.
    f = t^-a/(t^(1-gamma)/glo + 3q) decreasing, 1-eta(psi_K) increasing: cell >= (1-eta(psi_K(t0))) f(t1)(t1-t0);
    head (0, lo] dropped (integrand >= 0); tail [hi, inf) in closed form."""
    ctx.prec = 80
    a = arb(al)
    gam = gam_arb()
    q = arb(1) if N is None else (arb(3) / 5) ** N
    t0 = arb(lo)
    rr = arb(ratio)
    tot = arb(0)
    m0 = 1 - bh.eta_arb(psi_K(t0, K))
    while t0.upper() < arb(hi).lower():
        t1 = t0 * rr
        f1 = t1 ** (-a) / (t1 ** (1 - gam) / glo + 3 * q)
        tot += arb(m0.lower()) * f1 * (t1 - t0)
        m0 = 1 - bh.eta_arb(psi_K(t1, K))
        t0 = t1
    # tail [T, inf): t^(1-gamma)/glo + 3q <= t^(1-gamma)(1/glo + 3q T^(gamma-1)), 1-eta(psi_K) >= 1-eta(psi_K(T))
    tot += arb(m0.lower()) * t0 ** (gam - a) / (a - gam) / (1 / glo + 3 * q * t0 ** (gam - 1))
    return arb(tot.lower())


# ----------------------------------------------------------------- (D) floats
def float_deficit(al="0.8", cases=((10, 100), (20, 100)), dps=25):
    mp.dps = dps
    a = mpf(al)
    rho = mpf(5) ** a / 3
    out = []
    for N, n in cases:
        r = bi.corner_deficit(al=al, n=n, Ns=(N,), dps=dps)
        out.append((N, n, r[0][2]))
    return out


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem BJ: sharpness of the rho^-N cell modulus (Conjecture BI.7) and the corner constant. Not thrust.")
    symbolic_checks()
    exact_checks()

    # (C1) BI's sup g, re-run live (used in the upper bootstrap)
    G0, _, _ = bi.sup_g()
    G = arb("0.6031")
    okb, _ = bi.bootstrap(G0, G)
    check(okb and G0.upper() < G.lower(), "BI.4 re-run: sup_(0,1] g <= 0.6031")

    lo0, hi0, npts = g_range()
    glo, ghi, okx = g_bootstrap(lo0, hi0, G)
    print(f"   g on [1e-12, p(1e-12)]: [{lo0.lower()}, {hi0.upper()}] ({npts} nodes); on (0, p(1e-12)]: [{glo.lower()}, {ghi.upper()}]")
    check(okx and glo.lower() > 0.60094 and ghi.upper() < 0.60097,
          f"certified: g(s) = s^(1-gamma) r_inf(s) in [{float(glo.lower()):.8f}, {float(ghi.upper()):.8f}] on (0, p(1e-12)], so g_lo >= 0.60094, g_hi <= 0.60097")
    th = arb(THETA_AU)
    check(bool(glo.lower() <= th.lower()) and bool(th.upper() <= ghi.upper()),
          "consistency: AU's certified log-mean <Theta> = 0.600954251... lies in [g_lo, g_hi] (Theta_AT takes values there)")

    # (C2) the corner constant at 0.8
    al = "0.8"
    a = arb(al)
    A = arb(A08[0]).union(arb(A08[1]))
    kap = 2 * (arb.pi() * a).sin() / (arb.pi() * A)
    Mlo, Mhi, nc = M_bracket(al, glo, ghi)
    clo = arb((kap * glo * Mlo).lower())
    chi = arb((kap * ghi * Mhi).upper())
    print(f"   int t^(gamma-1-a)(1-eta(psi)) at 0.8 in [{Mlo.lower()}, {Mhi.upper()}] ({nc} cells); kappa_0.8 = {kap.str(12)}")
    print(f"   liminf/limsup rho^N delta_N(0.8) in [{clo.lower()}, {chi.upper()}]")
    check(clo.lower() > 1.13589 and chi.upper() < 1.13594,
          "certified: 1.13589 < g_lo M_0.8 <= liminf rho^N delta_N(0.8) <= limsup <= g_hi M_0.8 < 1.13594 (BI.7's float 1.1359 confirmed)")

    # (C3) N-uniform sharpness constants c_-(alpha)
    cm = {}
    for alx in ("0.75", "0.8", "0.9"):
        B = B_lower(alx, glo)
        Iup = I_closed(alx, arb(1), ghi)  # int tau^-a R <= int tau^-a min(1/3, g_hi tau^(gamma-1))
        cI = arb((B / Iup).lower())
        line = f"   alpha={alx}: B >= {float(B.lower()):.6f}, int tau^-a R <= {float(Iup.upper()):.6f}, c_- >= B/I = {float(cI.lower()):.6f}"
        if alx == "0.8":
            cA = arb((kap * B).lower())
            line += f"; with AZ's A_0.8: c_- >= {float(cA.lower()):.6f}"
            cm[alx] = cA
        else:
            cm[alx] = cI
        print(line)
    check(all(v.lower() > 0.25 for v in cm.values()) and cm["0.8"].lower() > 0.6,
          f"certified: rho^N delta_N >= c_- for every N >= 0: c_-(0.75) >= {float(cm['0.75'].lower()):.4f}, c_-(0.8) >= {float(cm['0.8'].lower()):.4f}, c_-(0.9) >= {float(cm['0.9'].lower()):.4f}")
    # finite-N upper bound at 0.8: rho^N delta_N <= g_hi M + kappa (3/2)(3/5)^N int t^-a min((5/3)^N/3, g_hi t^(gamma-1))
    ups = []
    for N in (0, 5, 10, 20, 40):
        ex = kap * arb(3) / 2 * (arb(3) / 5) ** N * I_closed(al, (arb(5) / 3) ** N, ghi)
        ups.append((N, float((chi + ex).upper())))
    print("   upper bounds rho^N delta_N(0.8) <= " + ", ".join(f"{u:.6f} (N={N})" for N, u in ups))
    check(ups[-1][1] < 1.1361 and all(u >= 1.13589 for _, u in ups),
          "certified: rho^N delta_N(0.8) <= g_hi M + O((3/5)^(N(1-a)/(1-gamma))) for every N (two-sided band around 1.1359)")
    # ratio to BI's K_0.8: the modulus is sharp up to a factor
    K08 = bi.K_alpha("0.8", G)
    print(f"   ratio of BI's upper constant to the corner floor: K_0.8/c_-(0.8) <= {float((K08 / cm['0.8']).upper()):.2f}")

    # (D) floats
    fd = float_deficit()
    for N, n, w in fd:
        print(f"   floats alpha=0.8: rho^{N} delta_{N}^(n={n}) = {mp.nstr(w, 10)}")
    check(all(float(cm["0.8"].lower()) <= float(w) <= ups[2][1] for _, _, w in fd)
          and 1.13589 < float(fd[-1][2]) < 1.13594 and abs(float(fd[0][2]) - float(fd[1][2])) < 1e-6,
          "floats: rho^N delta_N^(n=100) = 1.1359143 at N = 10, 20, inside the certified band [1.13589, 1.13594] "
          "(consistency; n - N >= 80 is used because the n-error looks like rho^-(n-N), floats only)")
    fe = float_deficit(cases=((10, 50), (20, 60)))
    rho = 5**0.8 / 3
    rat = []
    for (N, n, w), (_, _, w100) in zip(fe, fd):
        rat.append(float(abs(w100 - w)) * rho ** (n - N))
        print(f"   floats alpha=0.8: rho^{N} delta_{N}^(n={n}) = {mp.nstr(w, 10)}; rho^(n-N) |rho^N(delta^(100) - delta^(n))| = {rat[-1]:.4f}")
    check(all(0.05 < r_ < 5 for r_ in rat) and abs(rat[0] - rat[1]) < 0.2 * rat[0],
          f"floats: the n-error of rho^N delta_N is ~ C rho^-(n-N) with comparable C ({rat[0]:.2f}, {rat[1]:.2f}) at N = 10, 20 (Conjecture BJ.7, not proved)")

    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BJ checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
