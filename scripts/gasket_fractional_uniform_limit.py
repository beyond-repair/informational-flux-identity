#!/usr/bin/env python3
"""Theorem BI: above log3/log5 the fractional Dirichlet extension converges uniformly to a continuous limit
(Conjecture BG.2 proved), with an explicit cell modulus of order rho^-N.

Setting (GASKET.md, Theorems BG, BH): level-n gasket V_n, corner data d on C = {a, b, c}, x = d - dbar,
u_n the L_n^alpha-Dirichlet extension, w_s = (L_n+s)^-1 x_C, r_n(s), p(s) = s(s+5), D_n = int s^-a r_n,
rho = 5^alpha/3, gamma = log3/log5, R_n(tau) = (3/5)^n r_n(5^-n tau), g(s) = s^(1-gamma) r_inf(s).
A level-i cell of V_n is a copy 2^-i V_(n-i) (3^i of them); omega_i(s) = max over level-i cells of the
largest difference of w_s between two corners of the cell.

Theorem BI (checked here; exit 1 on any failure):
  (1) cell recursion (every alpha, n, s > 0): omega_0 <= 2 r_n ||x||, and for i < n
        omega_(i+1) <= (3/5) omega_i + kappa(z_i) r_n ||x||,   z_i = p^(n-i-1)(s),  kappa(z) = z/(z+2),
      because the midpoints of a level-i cell are the level-1 z_i-harmonic extension of its corners
      (weights alpha1, alpha1, delta, total 2/(z+2); m_kl - m_km = (f_l - f_m)/(z+5)) and |w_s| <= r_n ||x||;
      a chain through nested cells gives osc_cell_N(w_s) <= min(2, 2 sum_(i>=N) omega_i/(r_n||x||)) r_n ||x||.
  (2) orbit lemma: 5^k s <= 1  =>  p^k(s) <= e^(1/2) 5^k s; so kappa(z_l) <= min(1, 5^-(l+1) tau), tau = 5^n s.
  (3) hence, for every alpha, n >= max(N,1), every level-N cell Delta:
        osc_Delta u_n <= ||x|| (rho^n/D_n) int tau^-a R_n(tau) F_N(tau) dtau,
        F_N = min(2, 10(3/5)^N + 5 Phi_N),  Phi_N(tau) = sum_(l<N) (3/5)^(N-1-l) m(5^-(l+1) tau) + sum_(l>=N) m(5^-(l+1) tau).
  (4) two-sided bounds on R_n (n >= 1): R_n(tau) >= 1/(3+tau) (Dirichlet principle with the harmonic test function),
      so rho^n/D_n <= sin(pi a)/(3^-a pi); and R_n(tau) <= min(1/3, G tau^(gamma-1)) with G >= sup g, certified
      (arb cells on [1e-5, 1], a bootstrap along the p-orbit below 1e-5, and g <= s^-gamma (2+s)/(4+s) <= 3/5 on [1, inf)).
  (5) for alpha > gamma: osc_Delta u_n <= ||x|| eps_N(alpha) with eps_N <= K_alpha rho^-N -> 0, uniformly in n;
      with BG.3's pointwise limit on V_*, u_n -> u_inf uniformly and u_inf extends continuously to the gasket with
      u_inf = d on C (Conjecture BG.2).  Certified eps_N(0.8) for several N and K_alpha at 0.75, 0.8, 0.9.
  For alpha <= gamma the weight tau^(gamma-1-alpha) is not integrable: no modulus (Kept Failure BG.1 is the reason).

Floats (consistency only): dense L^0.8 cell oscillations at levels 5-7 lie under the bound; the corner deficit
(d_a - u_n(nu_N))/x_a times rho^N levels off (sharpness of rho^-N, Conjecture BI.7, not certified).

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
from mpmath import mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_maxwell_limit as bb  # noqa: E402
import gasket_fractional_corner_layer as bh  # noqa: E402

FAILS: list[str] = []
GAMMA = math.log(3) / math.log(5)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    z, s, t, E, Q, tt = sp.symbols("z s t E Q tt", positive=True)
    g = sp.Symbol("gamma", positive=True)
    Dd = (z + 2) * (z + 5)
    al1, de = (z + 4) / Dd, 2 / Dd
    check(sp.simplify(2 * al1 + de - 2 / (z + 2)) == 0 and sp.simplify(1 - (2 * al1 + de) - z / (z + 2)) == 0,
          "level-1 z-harmonic measure: total mass 2 al1 + delta = 2/(z+2), killed mass kappa(z) = z/(z+2) <= min(1, z/2)")
    check(sp.simplify(al1 - de - 1 / (z + 5)) == 0 and sp.simplify(al1 + de - (z + 6) / Dd) == 0,
          "m_kl - m_km = (al1 - delta)(f_l - f_m) = (f_l - f_m)/(z+5); |m_kl - f_k| <= beta omega + kappa |f_k|, beta = al1+delta <= 3/5 (BH)")
    c = (6 + s) / ((5 + s) * (2 + s))
    check(sp.expand(sp.together((c - (sp.Rational(3, 5) - sp.Rational(8, 25) * s)) * (5 + s) * (2 + s)))
          == sp.expand(sp.Rational(41, 25) * s**2 + sp.Rational(8, 25) * s**3),
          "c(s) - (3/5 - 8s/25) = (41 s^2/25 + 8 s^3/25)/((s+2)(s+5)) >= 0 (bootstrap lemma)")
    ratio = sp.powsimp(sp.expand_power_base(s ** (1 - g) * (s * (s + 5)) ** (g - 1), force=True), force=True)
    check(sp.simplify(ratio - (s + 5) ** (g - 1)) == 0 and abs(float(5 ** (sp.log(3) / sp.log(5) - 1)) - 0.6) < 1e-15,
          "g(s)/g(ps) = (s+5)^(gamma-1)/(c - e r(ps)) = (3/5)(1+s/5)^(gamma-1)/(c - e r(ps)) (5^(gamma-1) = 3/5)")
    # g <= s^-gamma (2+s)/(4+s) on [1, inf), decreasing: d/ds log = -gamma/s + 2/((2+s)(4+s))
    dl = sp.together(-g / s + 2 / ((2 + s) * (4 + s)))
    num = sp.expand(sp.numer(dl) * -1)
    check(sp.expand(num - (g * (s**2 + 6 * s + 8) - 2 * s)) == 0 and all(
        float((sp.log(3) / sp.log(5)) * (v**2 + 6 * v + 8) - 2 * v) > 0 for v in np.linspace(0, 50, 501)),
          "eta <= 2/(2+s) gives r_inf <= (2+s)/(s(s+4)); f = s^-gamma (2+s)/(4+s) has log-derivative of sign -(gamma(s^2+6s+8)-2s) < 0 "
          "(discriminant (6gamma-2)^2 - 32gamma^2 < 0), so g <= f(1) = 3/5 on [1, inf)")
    gm = sp.log(3) / sp.log(5)
    check(float((6 * gm - 2) ** 2 - 32 * gm**2) < 0, "discriminant (6 gamma - 2)^2 - 32 gamma^2 < 0")
    # Dirichlet principle
    check(sp.simplify(sp.diff(tt - tt**2 * (E + s * Q) / 2, tt).subs(tt, 1 / (E + s * Q))) == 0
          and sp.simplify((tt - tt**2 * (E + s * Q) / 2).subs(tt, 1 / (E + s * Q)) - 1 / (2 * (E + s * Q))) == 0,
          "1/2 v^T A^-1 v >= max_t (t v.h - t^2 h^T A h / 2) = 1/(2 h^T A h) with v.h = 1")
    n = sp.Symbol("n", positive=True, integer=True)
    lb = (sp.Rational(3, 5) ** n) / (3 * sp.Rational(3, 5) ** n + 5 ** (-n) * tt * 3 * (3**n + 1) / 4)
    check(sp.simplify(lb - 1 / (3 + tt * sp.Rational(3, 4) * (1 + 3 ** (-n)))) == 0,
          "R_n(tau) >= (3/5)^n/(3(3/5)^n + 5^-n tau 3(3^n+1)/4) = 1/(3 + (3/4)(1+3^-n) tau) >= 1/(3+tau) for n >= 1")
    a = sp.Rational(4, 5)
    val = sp.integrate(tt ** (-a) / (3 + tt), (tt, 0, sp.oo))
    check(abs(float(val) - float(3 ** (-a) * sp.pi / sp.sin(sp.pi * a))) < 1e-12,
          "int_0^inf tau^-a/(3+tau) = 3^-a pi/sin(pi a) (checked at a = 4/5)")
    check(math.exp(0.5) / 4 < 0.5 and math.exp(0.5) / 2 < 1,
          "orbit lemma constants: e^(1/2)/4 < 1/2 (induction closes), e^(1/2)/2 < 1 (kappa(z) <= min(1, 5^k s))")


# ----------------------------------------------------------------- (B) exact rationals
def cells(N):
    out = [(Fraction(0), Fraction(0))]
    for k in range(N):
        hk = Fraction(1, 2 ** (k + 1))
        out = [(x + dx, y + dy) for (x, y) in out for (dx, dy) in ((0, 0), (hk, 0), (0, hk))]
    return out


def F_N(N: int, tau, one=1):
    """F_N(tau) = min(2, 10 (3/5)^N + 5 Phi_N(tau)); exact for Fractions, generic for arb (upper bound)."""
    th = Fraction(3, 5) if isinstance(tau, Fraction) else arb(3) / 5
    five = 5 if isinstance(tau, Fraction) else arb(5)

    def m(v):
        return v if v < 1 else one

    A = 0 * tau
    for l in range(N):
        A += th ** (N - 1 - l) * m(tau / five ** (l + 1))
    B, l = 0 * tau, N
    while tau / five ** (l + 1) >= 1:
        B += 1
        l += 1
    B += tau / five ** (l + 1) * five / 4  # sum_{j>=l} tau 5^-(j+1) = tau 5^-(l+1) * 5/4, all terms < 1
    val = 10 * th**N + 5 * (A + B)
    return val if val < 2 else 2 * one


def exact_checks(levels=(2, 3, 4), svals=(Fraction(1, 100), Fraction(1, 2), Fraction(7, 5))):
    xC = (Fraction(3), Fraction(-1), Fraction(-2))
    xn = 3
    ok_mid = ok_rec = ok_osc = ok_bound = ok_max = True
    for n in levels:
        adj, C, pos = bh.mesh(n)
        for s in svals:
            rhs = [Fraction(0)] * len(adj)
            for kk, cc in enumerate(C):
                rhs[cc] = xC[kk]
            w = bh.exact_solve(adj, s, rhs)
            rn = bb.rec_all(n, s)[0]
            M = rn * xn
            ok_max = ok_max and max(abs(v) for v in w) == M
            orb = [s]
            for _ in range(n):
                orb.append(orb[-1] * (orb[-1] + 5))
            om = []
            for i in range(n + 1):
                h = Fraction(1, 2**i)
                best = Fraction(0)
                for (x0, y0) in cells(i):
                    f = [w[pos[(x0, y0)]], w[pos[(x0 + h, y0)]], w[pos[(x0, y0 + h)]]]
                    best = max(best, max(f) - min(f))
                    if i < n:
                        zz = orb[n - i - 1]
                        a1, dl = (zz + 4) / ((zz + 2) * (zz + 5)), 2 / ((zz + 2) * (zz + 5))
                        hh = h / 2
                        mids = [w[pos[(x0 + hh, y0)]], w[pos[(x0 + hh, y0 + hh)]], w[pos[(x0, y0 + hh)]]]
                        pred = [a1 * (f[0] + f[1]) + dl * f[2], a1 * (f[1] + f[2]) + dl * f[0], a1 * (f[2] + f[0]) + dl * f[1]]
                        ok_mid = ok_mid and mids == pred
                om.append(best)
            ok_rec = ok_rec and om[0] <= 2 * M and all(
                om[i + 1] <= Fraction(3, 5) * om[i] + orb[n - i - 1] / (orb[n - i - 1] + 2) * M for i in range(n))
            tau = s * 5**n
            for N in range(n + 1):
                h = Fraction(1, 2**N)
                S = sum(om[N:])
                FN = F_N(N, tau)
                for (x0, y0) in cells(N):
                    vals = [w[i] for key, i in pos.items() if key[0] >= x0 and key[1] >= y0 and key[0] - x0 + key[1] - y0 <= h]
                    osc = max(vals) - min(vals)
                    ok_osc = ok_osc and osc <= min(2 * M, 2 * S)
                    ok_bound = ok_bound and osc <= M * FN
    check(ok_max, f"exact, levels {levels[0]}-{levels[-1]}, s in {{1/100,1/2,7/5}}: max |w_s| = r_n ||x|| (attained on the corners)")
    check(ok_mid, "exact: on every level-i cell the three midpoints are the level-1 z-harmonic extension of the cell corners, "
                  "z = p^(n-i-1)(s) (decimation inside a cell)")
    check(ok_rec, "exact: omega_0 <= 2 r_n||x||, omega_(i+1) <= (3/5) omega_i + z_i/(z_i+2) r_n ||x|| (BI.1)")
    check(ok_osc, "exact: osc of w_s on every level-N cell <= min(2 r_n||x||, 2 sum_(i>=N) omega_i) (chain through nested cells)")
    check(ok_bound, "exact: osc of w_s on every level-N cell <= r_n ||x|| F_N(5^n s) (BI.1 + orbit lemma)")
    # orbit lemma and Dirichlet lower bound, exact
    ok = True
    for k in range(0, 12):
        for q in (Fraction(1, 5**k), Fraction(1, 3 * 5**k), Fraction(1, 1000 * 5**k)):
            v = q
            for _ in range(k):
                v = v * (v + 5)
            ok = ok and float(v) <= math.exp(0.5) * float(q * 5**k) * (1 + 1e-15)
    check(ok, "orbit lemma p^k(s) <= e^(1/2) 5^k s when 5^k s <= 1 (exact orbits, k <= 11)")
    ok = True
    for n in (1, 2, 3, 4, 5, 6):
        for s in (Fraction(0), Fraction(1, 7), Fraction(3), Fraction(1, 5**n)):
            rn = bb.rec_all(n, s)[0] if s > 0 else Fraction(5**n, 3 * 3**n)
            ok = ok and rn >= 1 / (3 * Fraction(3, 5) ** n + s * Fraction(3 * (3**n + 1), 4))
    check(ok, "exact: r_n(s) >= 1/(3(3/5)^n + s 3(3^n+1)/4) at levels 1-6 (Dirichlet principle; equality at s = 0)")
    ok = True
    for n in (1, 2, 3, 4):
        adj, C, pos = bh.mesh(n)
        hv = bh.exact_solve(adj, Fraction(0), [Fraction(0)] * len(adj), {C[0]: Fraction(1, 2), C[1]: Fraction(-1, 2), C[2]: Fraction(0)})
        En = sum((hv[i] - hv[j]) ** 2 for i in range(len(adj)) for j in adj[i] if i < j)
        ok = ok and En == Fraction(3, 2) * Fraction(3, 5) ** n and len(adj) == 3 * (3**n + 1) // 2 and max(abs(v) for v in hv) == Fraction(1, 2)
    check(ok, "exact: the harmonic test function with corner values (1/2,-1/2,0) has energy (3/2)(3/5)^n, |h| <= 1/2, |V_n| = 3(3^n+1)/2")


# ----------------------------------------------------------------- (C) certificates
def sup_g(eps0="1e-5", ratio="1.002", prec=80):
    """Upper bound G0 for g = s^(1-gamma) r_inf on [eps0, 1] (r_inf decreasing, s^(1-gamma) increasing on cells)."""
    ctx.prec = prec
    gam = arb(3).log() / arb(5).log()
    pts = [arb(eps0)]
    rr = arb(ratio)
    while pts[-1].upper() < 1:
        pts.append(pts[-1] * rr)
    pts[-1] = arb(1)
    rup = [(1 / (2 + p - 2 * bh.eta_arb(p))) for p in pts]
    G0 = arb(0)
    lo = arb(10)
    for i in range(len(pts) - 1):
        up = (pts[i + 1] ** (1 - gam) * rup[i]).upper()
        dn = (pts[i] ** (1 - gam) * rup[i + 1]).lower()
        if up > G0:
            G0 = arb(up)
        if dn < lo:
            lo = arb(dn)
    return G0, lo, len(pts)


def bootstrap(G0: arb, G: arb, eps0="1e-5"):
    gam = arb(3).log() / arb(5).log()
    e = arb(eps0)
    pe = e * (e + 5)
    lhs = G0 * (arb(4) / 15 * pe + G * pe**gam).exp()
    xmax = arb(8) / 15 * pe / 5 + G * (pe / 5) ** gam
    return bool(lhs.upper() <= G.lower()) and bool(xmax.upper() < 0.5) and bool(pe.upper() <= 1), lhs


def eps_bound(al: str, N: int, G: arb, lo="1e-40", hi="1e140", ratio="1.02"):
    """Upper bound for eps_N(alpha) = S_a int_0^inf w(tau) F_N(tau) dtau, w = tau^-a min(1/3, G tau^(gamma-1)),
    S_a = sin(pi a)/(3^-a pi).  F_N nondecreasing: cell [t0,t1] <= F_N(t1) int w (exact antiderivatives);
    head [0, lo] <= F_N(lo) lo^(1-a)/(3(1-a)); tail [hi, inf) <= 2 G hi^(gamma-a)/(a-gamma)."""
    ctx.prec = 80
    a = arb(al)
    gam = arb(3).log() / arb(5).log()
    tc = (3 * G) ** (1 / (1 - gam))  # min switches at G tau^(gamma-1) = 1/3

    def Wint(t0, t1):  # int_t0^t1 w, t0 < t1, both on one side of tc
        if t1.upper() <= tc.lower():
            return (t1 ** (1 - a) - t0 ** (1 - a)) / (3 * (1 - a))
        if t0.lower() >= tc.upper():
            return G * (t0 ** (gam - a) - t1 ** (gam - a)) / (a - gam)
        # straddles: bound w by the larger branch on the cell: tau^-a/3 >= w and G tau^(g-1-a) >= w
        return arb(min((t1 ** (1 - a) - t0 ** (1 - a)) / (3 * (1 - a)), G * (t0 ** (gam - a) - t1 ** (gam - a)) / (a - gam), key=lambda v: v.upper()).upper())

    S = (arb.pi() * a).sin() / (arb(3) ** (-a) * arb.pi())
    t = arb(lo)
    rr = arb(ratio)
    tot = F_N(N, t, arb(1)) * t ** (1 - a) / (3 * (1 - a))
    hiA = arb(hi)
    while t.upper() < hiA.lower():
        t1 = t * rr
        tot += F_N(N, t1, arb(1)) * Wint(t, t1)
        t = t1
    tot += 2 * G * t ** (gam - a) / (a - gam)
    return arb((S * tot).upper())


def K_alpha(al: str, G: arb):
    """Analytic constant: eps_N <= S_a (10 (3/5)^N I_w + G rho^-N J_a) <= K_a rho^-N,
    I_w = int w, J_a = int t^(gamma-1-a) min(2, 21.25 t^(1-gamma)) dt."""
    ctx.prec = 80
    a = arb(al)
    gam = arb(3).log() / arb(5).log()
    tc = (3 * G) ** (1 / (1 - gam))
    Iw = tc ** (1 - a) / (3 * (1 - a)) + G * tc ** (gam - a) / (a - gam)
    ts = (2 / arb("21.25")) ** (1 / (1 - gam))
    Ja = arb("21.25") * ts ** (1 - a) / (1 - a) + 2 * ts ** (gam - a) / (a - gam)
    S = (arb.pi() * a).sin() / (arb(3) ** (-a) * arb.pi())
    return arb((S * (10 * Iw + G * Ja)).upper())


# ----------------------------------------------------------------- (D) floats
def dense_osc(n, al, d=(1.0, -0.5, 0.0)):
    u, pos = bh.dense_u(n, al, d)
    keys = list(pos.keys())
    X = np.array([float(k[0]) for k in keys])
    Y = np.array([float(k[1]) for k in keys])
    U = np.array([u[pos[k]] for k in keys])
    out = []
    for N in range(n + 1):
        h = 2.0**-N
        best = 0.0
        for (x0, y0) in cells(N):
            msk = (X >= float(x0) - 1e-12) & (Y >= float(y0) - 1e-12) & (X - float(x0) + Y - float(y0) <= h + 1e-12)
            v = U[msk]
            best = max(best, v.max() - v.min())
        out.append(best)
    return out, u, pos


def corner_deficit(al="0.8", n=40, Ns=range(0, 16), dps=20):
    """Floats: delta_N = (d_a - u_n(nu_N))/x_a for symmetric data x_b = x_c (h_N(z) = x_a eta_N(z), exact recursion
    eta_N = alpha1 + beta eta_(N-1)(pz), eta_0 = -1/2): delta_N = D_n^-1 int s^-a r_n (1 - eta_N(p^(n-N) s))."""
    mp.dps = dps
    a = mpf(al)

    def etaN(z, N):
        zs = [z]
        for _ in range(N):
            zs.append(zs[-1] * (zs[-1] + 5))
        e = mpf(-1) / 2
        for k in range(N - 1, -1, -1):
            zz = zs[k]
            D = (zz + 2) * (zz + 5)
            e = (zz + 4) / D + (zz + 6) / D * e
        return e

    def pk(s, k):
        for _ in range(k):
            s = s * (s + 5)
        return s

    pts = [0] + [mpf(5) ** (-n) * mpf(10) ** j for j in range(-6, 60, 2)] + [mp.inf]
    D = mp.quad(lambda s: s ** (-a) * bb.rec_all(n, s)[0], pts)
    rho = mpf(5) ** a / 3
    res = []
    for N in Ns:
        v = mp.quad(lambda s: s ** (-a) * bb.rec_all(n, s)[0] * (1 - etaN(pk(s, n - N), N)), pts) / D
        res.append((N, v, v * rho**N))
    return res


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem BI: uniform convergence and continuity of the fractional Dirichlet extension above log3/log5. Not thrust.")
    symbolic_checks()
    exact_checks()

    G0, gmin, npts = sup_g()
    print(f"   sup g on [1e-5, 1] <= {G0.upper()} ({npts} nodes); inf g there >= {gmin.lower()}")
    G = arb("0.6031")
    okb, lhs = bootstrap(G0, G)
    check(okb and G0.upper() < G.lower(), f"bootstrap below 1e-5: G0 exp((4/15)p(eps0) + G p(eps0)^gamma) <= {lhs.upper()} <= G = 0.6031, x_j < 1/2")
    check(bool(arb(3) / 5 <= G), "G = 0.6031 also covers [1, inf) (g <= 3/5 there); so R_n(tau) <= min(1/3, 0.6031 tau^(gamma-1)) for all n, tau")

    rho = 5**0.8 / 3
    epsv = {}
    for N in (0, 5, 10, 20, 30, 40, 60):
        e = eps_bound("0.8", N, G)
        epsv[N] = float(e.upper())
        print(f"   eps_{N}(0.8) <= {float(e.upper()):.6g}   (rho^N eps_N <= {float(e.upper()) * rho**N:.4f})")
    check(all(epsv[a] >= epsv[b] for a, b in zip(list(epsv), list(epsv)[1:])) and epsv[40] < 0.05 and epsv[60] < 3e-3,
          "certified (arb, outward): eps_N(0.8) decreases, eps_40(0.8) < 0.05, eps_60(0.8) < 3e-3")
    Ks = {}
    for al in ("0.75", "0.8", "0.9"):
        Ks[al] = K_alpha(al, G)
        print(f"   K_{al} <= {float(Ks[al].upper()):.6g}: eps_N <= K rho^-N for every N")
    check(all(epsv[N] <= float(Ks["0.8"].upper()) * rho ** (-N) * (1 + 1e-12) for N in epsv),
          "the numerically integrated eps_N(0.8) sit under the analytic K_0.8 rho^-N")

    # dense floats at 0.8: oscillations under the certified bound; uniform Cauchy
    xn = 5 / 6
    okd = True
    prev = None
    for n in (5, 6, 7):
        osc, u, pos = dense_osc(n, 0.8)
        print(f"   dense alpha=0.8 level {n}: osc on level-N cells, N=0..{n}: " + ", ".join(f"{v:.4f}" for v in osc))
        for N in range(n + 1):
            bnd = xn * (epsv[N] if N in epsv else float(Ks["0.8"].upper()) * rho ** (-N))
            okd = okd and osc[N] <= min(bnd, 2 * xn) + 1e-9
        if prev is not None:
            up, pp = prev
            diff = max(abs(u[pos[k]] - up[pp[k]]) for k in pp)
            print(f"   max over V_{n-1} |u_{n} - u_{n-1}| = {diff:.5f}")
        prev = (u, pos)
    check(okd, "floats: dense L^0.8 cell oscillations at levels 5-7 are below ||x|| min(2, eps_N) (consistency)")
    osc45, _, _ = dense_osc(6, 0.45)
    check(osc45[6] > 0.45, f"floats: at alpha=0.45 the finest-cell oscillation stays large ({osc45[6]:.4f} at level 6): the corner layer, no modulus below gamma")

    cd = corner_deficit(n=40, Ns=range(0, 16))
    for N, v, w in cd[::3]:
        print(f"   floats alpha=0.8, n=40: delta_{N} = {mp.nstr(v, 10)}, rho^N delta_N = {mp.nstr(w, 10)}")
    cd60 = corner_deficit(n=60, Ns=(5, 10))
    for N, v, w in cd60:
        print(f"   floats alpha=0.8, n=60: delta_{N} = {mp.nstr(v, 10)}, rho^N delta_N = {mp.nstr(w, 10)}")
    ws = [float(w) for _, _, w in cd]
    check(all(v > 0 for _, v, _ in cd) and abs(ws[-1] - ws[-2]) < 2e-3 * ws[-1] and 0.05 < ws[-1] < 5
          and abs(float(cd60[0][2]) - float(cd60[1][2])) < 2e-4 and abs(float(cd60[1][2]) - 1.1359) < 1e-3,
          "floats: the corner deficit delta_N decays like rho^-N (rho^N delta_N ~ 1.1359 at n=60): the modulus order looks sharp (Conjecture BI.7)")

    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BI checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
