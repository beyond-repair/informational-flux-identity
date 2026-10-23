#!/usr/bin/env python3
"""Theorem BK: the uniform convergence u_n -> u-bar of Theorem BI has rate rho^-n on all of V_n (Conjecture BJ.7 proved,
with exponent c = 1), by an exact one-level transfer identity whose weight is >= 1.

Setting (GASKET.md, Theorems BG-BJ): level-n gasket V_n, corner data d, x = d - dbar, u_n the L_n^alpha-Dirichlet
extension, u-bar its limit, r_n(s) the corner resolvent (AQ), p(y) = y(y+5), beta(y) = (y+6)/((y+2)(y+5)) (= AQ's c),
e(y) = 2y/(2+y), rho = 5^alpha/3, gamma = log3/log5, R_n(tau) = (3/5)^n r_n(5^-n tau), D'_n = int tau^-a R_n = rho^-n D_n.

Theorem BK (checked here; exit 1 on any failure):
  (1) exact transfer (every alpha in (0,1), n >= m, x in V_m): with mu_n = tau^-a R_n(tau) dtau / D'_n and
      F_n(tau) = h_m(p^(n-m)(5^-n tau); x) (BG's profile, u_n(x) - dbar = <F_n>_(mu_n)),
        u_(n+1)(x) - dbar = <omega_n F_n>_(mu_n) / <omega_n>_(mu_n),   <omega_n>_(mu_n) = D'_(n+1)/D'_n,
      omega_n(tau') = (3/5) (1+y/5)^a r_(n+1)(y) / ((1+2y/5) r_n(p(y))), tau' = 5^n p(y)
                    = (1+y/5)^a / ((1+2y/5)(5/3)(beta(y) - e(y) r_n(p(y))))    (AQ's step).
  (2) omega_n >= 1 everywhere, every alpha >= 0: (1+2y/5)(5/3)beta(y) = 1 - y(y+4)/(3(y+2)(y+5)) <= 1 (sympy).
      Hence D'_n = rho^-n D_n is strictly increasing in n for every alpha (derived), and, since
      |Cov(omega,F)| <= <|omega - <omega>|> sup|F - c| <= 2(<omega>-1) sup|F - c|,
        |u_(n+1)(x) - u_n(x)| <= 2 sup|F_n - c| (1 - D'_n/D'_(n+1)) <= 2 sup|F_n - c| log(D'_(n+1)/D'_n).
  (3) for alpha > gamma (D'_n increases to D'_inf = pi A_a/(2 sin pi a) < inf): for every n and x in V_n,
        |u_n(x) - u-bar(x)| <= 2 ||x|| log(D'_inf/D'_n),
      and at the corner neighbours nu^N (N <= n) the data-split form
        |u_n(nu_b^N) - u-bar(nu_b^N)| <= ((3/2)|x_a| + 5^-N |x_b - x_c|/2) log(D'_inf/D'_n)
      (h_N = x_a eta_N + ((x_b-x_c)/2) prod_(i<N)(5+p^i z)^-1, eta_N in [-1/2, 1]); data-free for the deficit:
        |delta_N^(n) - delta_N| <= (3/2) log(D'_inf/D'_n).
  (4) certified rate: D'_(k+1) - D'_k = rho^-k int_0^inf (omega_k(y)-1) p^-a r_k(p) p' dy <= J_a rho^-k, with
      J_a = int (omega_hi - 1) p^-a rbar(p) p' dy, rbar(s) = min(1/s, G s^(gamma-1)) (BI.4's G = 0.6031, r_k <= r_inf,
      1/r_inf = 2+s-2eta >= s), omega_hi = min((1+y/5)^a/((1+2y/5)(5/3)(beta - e rbar(p))), (3/5)(1+y/5)^a rbar(y)(3+p)/(1+2y/5))
      (r_k >= r_0 = 1/(3+.)).  So log(D'_inf/D'_n) <= (D'_inf - D'_n)/D'_n <= J_a rho^-n/((1-1/rho) D'_0),
      D'_0 = 3^-a pi / sin(pi a):  sup_(V_n)|u_n - u-bar| <= C_a ||x|| rho^-n for every n.  Certified (arb) J_a and C_a
      at 0.75, 0.8, 0.9; with AZ's A_0.8 the asymptotic constant 2 J/((1-1/rho) D'_inf).

Derived from AZ (3) as stated (C_0.8, A_0.8): rho^n log(D'_inf/D'_n) -> -C/A = -nu_0.8 = 0.72744..., so
limsup rho^n sup|u_n - u-bar| / ||x|| <= 1.45489.  Floats (consistency only): dense L^0.8 step differences under the
per-step bound; deficit n-errors under the bound and ~ 0.6 rho^-n.  Nothing here is thrust.
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
import gasket_fractional_uniform_limit as bi  # noqa: E402
import gasket_fractional_corner_sharpness as bj  # noqa: E402

FAILS: list[str] = []
A08 = ("2.39640970397922052703", "2.39640970397922052705")  # Theorem AZ (3), used as stated
C08 = "-1.74325717839207572395"  # Theorem AZ intercept C_0.8, used as stated (derived part only)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    y, a, X, n, T = sp.symbols("y a X n T", positive=True)
    beta = (y + 6) / ((y + 2) * (y + 5))
    e = 2 * y / (2 + y)
    p = y * (y + 5)
    B = sp.Rational(5, 3) * beta
    check(sp.simplify((1 + 2 * y / 5) * B - (1 - y * (y + 4) / (3 * (y + 2) * (y + 5)))) == 0,
          "(1+2y/5)(5/3)beta(y) = 1 - y(y+4)/(3(y+2)(y+5)) <= 1, so omega_n >= (1+y/5)^a >= 1 for every alpha >= 0")
    # change of variables tau = 5^(n+1) y, tau' = 5^n p(y)
    tau, taup = 5 ** (n + 1) * y, 5**n * p
    check(sp.simplify(taup / tau - (1 + y / 5)) == 0 and sp.simplify(sp.diff(taup, y) / sp.diff(tau, y) - (1 + 2 * y / 5)) == 0,
          "tau' = 5^n p(5^-(n+1) tau): tau'/tau = 1 + y/5, dtau'/dtau = 1 + 2y/5, and p^(n+1-m)(5^-(n+1) tau) = p^(n-m)(5^-n tau')")
    step = X / (beta - e * X)  # AQ: r_(n+1)(y) = r_n(p y)/(beta - e r_n(p y)), X = r_n(p(y))
    om1 = sp.Rational(3, 5) * step / ((1 + 2 * y / 5) * X)
    om2 = 1 / ((1 + 2 * y / 5) * sp.Rational(5, 3) * (beta - e * X))
    check(sp.simplify(om1 - om2) == 0,
          "omega_n / (1+y/5)^a = (3/5) r_(n+1)(y)/((1+2y/5) r_n(p)) = 1/((1+2y/5)(5/3)(beta - e r_n(p))): increasing in r_n(p)")
    # the measure: tau'^-a R_n(tau') dtau' = rho^-n p^-a r_n(p) p' dy
    k = sp.symbols("k", positive=True)
    lhs = (5**k * p) ** (-a) * sp.Rational(3, 5) ** k * 5**k
    check(sp.simplify(sp.expand_power_base(lhs / (3**k * 5 ** (-a * k) * p ** (-a)), force=True)) == 1,
          "tau'^-a R_k(tau') dtau' = rho^-k p(y)^-a r_k(p(y)) p'(y) dy with rho^-k = 3^k 5^(-a k)")
    # Cov inequality: (1/<w>) |<(w-<w>)(F-c)>| <= 2 (<w>-1)/<w> sup|F-c| when w >= 1; and 1 - 1/t <= log t
    t = sp.symbols("t", positive=True)
    f = sp.log(t) - (1 - 1 / t)
    check(sp.simplify(sp.diff(f, t) - (t - 1) / t**2) == 0 and f.subs(t, 1) == 0,
          "1 - 1/t <= log t for t >= 1 (derivative (t-1)/t^2 >= 0), so the steps telescope to 2 sup|F-c| log(D'_inf/D'_n)")
    # head: (1+2y/5)(5/3) e(y) G (5y)^(gamma-1) <= (1+2y/5) G y^gamma (e <= y, 5^gamma = 3)
    check(sp.simplify(y - e - y**2 / (2 + y)) == 0,
          "e(y) = y - y^2/(2+y) <= y; p >= 5y; with 5^(gamma-1) = 3/5 the head term is <= (1+2y/5) G y^gamma")
    # tail: (y^2+5y+3) <= (2y+5)(y/2+5/4)
    check(sp.expand((2 * y + 5) * (y / 2 + sp.Rational(5, 4)) - (y**2 + 5 * y + 3)) == sp.Rational(13, 4),
          "tail: (3+p)/(1+2y/5) = 5(y^2+5y+3)/(2y+5) <= 5(y/2+5/4)")
    # eta_N in [-1/2, 1]: eta_N = alpha1 + beta eta_(N-1)(pz), eta_0 = -1/2, alpha1, beta > 0, eta_N <= eta <= 1
    z = sp.symbols("z", positive=True)
    al1 = (z + 4) / ((z + 2) * (z + 5))
    bz = (z + 6) / ((z + 2) * (z + 5))
    check(sp.simplify(al1 - bz / 2 - (z + 2) / (2 * (z + 2) * (z + 5))) == 0,
          "eta_N >= alpha1 - beta/2 = 1/(2(z+5)) > 0 for N >= 1 (eta_0 = -1/2) and eta_N <= eta <= 1 (BJ.4): eta_N in [-1/2, 1]")


# ----------------------------------------------------------------- (B) exact rationals
def exact_checks(levels=(0, 1, 2, 3)):
    ys = (Fraction(1, 1000), Fraction(1, 10), Fraction(1), Fraction(7, 2), Fraction(20))
    ok1 = ok2 = ok3 = True
    for n in levels:
        for y in ys:
            p = y * (y + 5)
            r1 = bj.r_exact(n + 1, y)
            r0 = bj.r_exact(n, p)
            beta = (y + 6) / ((y + 2) * (y + 5))
            e = 2 * y / (2 + y)
            ok1 = ok1 and r1 == r0 / (beta - e * r0)
            w = Fraction(3, 5) * r1 / ((1 + 2 * y / 5) * r0)  # omega / (1+y/5)^a
            ok2 = ok2 and w >= 1
            # crude ceiling U2 ingredients: r_(n+1)(y) <= 1/y and r_n(p) >= 1/(3+p)
            ok3 = ok3 and r1 <= 1 / y and r0 >= 1 / (3 + p)
    check(ok1, f"exact rationals, levels {levels[0]}-{levels[-1]}+1: r_(n+1)(y) = r_n(p(y))/(beta(y) - e(y) r_n(p(y))) at 5 rational y (AQ's step)")
    check(ok2, "exact rationals: (3/5) r_(n+1)(y)/((1+2y/5) r_n(p(y))) >= 1, i.e. omega_n >= 1 (BK.2)")
    check(ok3, "exact rationals: r_(n+1)(y) <= 1/y and r_n(p) >= r_0(p) = 1/(3+p) (the ingredients of omega_hi)")


# ----------------------------------------------------------------- (C) certificates
def _min_up(u: arb, v: arb) -> arb:
    uu, vv = arb(u.upper()), arb(v.upper())
    return uu if uu < vv else vv


def J_bound(al: str, G: arb, ylo="1e-10", yhi="1e10", ratio="1.001", prec=96):
    """Upper bound for J_a = int_0^inf (omega_hi - 1) p^-a rbar(p) p' dy: arb cells on [ylo, yhi] (ball evaluation of
    the whole integrand on each cell), closed-form head and tail."""
    ctx.prec = prec
    a = arb(al)
    gam = arb(3).log() / arb(5).log()
    five3 = arb(5) / 3
    pts = [arb(ylo)]
    rr = arb(ratio)
    while pts[-1].upper() < arb(yhi):
        pts.append(pts[-1] * rr)
    tot = arb(0)
    for i in range(len(pts) - 1):
        y0, y1 = pts[i], pts[i + 1]
        Y = y0.union(y1)
        p = Y * (Y + 5)
        pp = 2 * Y + 5
        rp = _min_up(1 / p, G * p ** (gam - 1))
        ry = _min_up(1 / Y, G * Y ** (gam - 1))
        beta = (Y + 6) / ((Y + 2) * (Y + 5))
        e = 2 * Y / (2 + Y)
        num = (1 + Y / 5) ** a
        U2 = arb((arb(3) / 5 * num * ry * (3 + p) / (1 + 2 * Y / 5)).upper())
        den = (1 + 2 * Y / 5) * five3 * (beta - e * rp)
        w = U2
        if den.lower() > 0:
            w = _min_up(num / den, U2)
        wm1 = w - 1
        if wm1.upper() < 0:
            wm1 = arb(0)
        wt = arb((p ** (-a) * rp * pp).upper())
        tot += arb((wm1 * wt * (y1 - y0)).upper())
    # head (0, y0]: omega_hi - 1 <= (a y/5 + delta(y))/(1 - delta(y0)), delta = (4+y0) y/30 + (1+2y0/5) G y^gamma;
    # weight <= G (5y)^(gamma-1-a) (5+2y0)
    y0 = arb(ylo)
    d0 = (4 + y0) * y0 / 30 + (1 + 2 * y0 / 5) * G * y0**gam
    c1 = a / 5 + (4 + y0) / 30
    c2 = (1 + 2 * y0 / 5) * G
    pref = G * arb(5) ** (gam - 1 - a) * (5 + 2 * y0) / (1 - d0)
    head = pref * (c1 * y0 ** (gam + 1 - a) / (gam + 1 - a) + c2 * y0 ** (2 * gam - a) / (2 * gam - a))
    # tail [Y0, inf): integrand <= 3 (1/5+1/Y0)^a (1/2+5/(4Y0)) (2+5/Y0) y^(-a-1)
    Y0 = pts[-1]
    K = 3 * (arb(1) / 5 + 1 / Y0) ** a * (arb(1) / 2 + 5 / (4 * Y0)) * (2 + 5 / Y0)
    tail = K * Y0 ** (-a) / a
    J = arb((tot + head + tail).upper())
    return J, len(pts) - 1, float(head.upper()), float(tail.upper())


def constants(al: str, J: arb):
    a = arb(al)
    rho = arb(5) ** a / 3
    D0 = arb(3) ** (-a) * arb.pi() / (arb.pi() * a).sin()
    eps_coef = J / (1 - 1 / rho)  # D'_inf - D'_n <= eps_coef rho^-n
    C = 2 * eps_coef / D0
    return rho, D0, eps_coef, C


# ----------------------------------------------------------------- (D) floats
def r_float(n, s):
    return bb.rec_all(n, s)[0]


def Dprime(n, al, dps=20):
    mp.dps = dps
    a = mpf(al)
    pts = [0] + [mpf(10) ** j for j in range(-8, 40)] + [mp.inf]
    sc = (mpf(3) / 5) ** n
    return mp.quad(lambda t: t ** (-a) * sc * r_float(n, mpf(5) ** (-n) * t), pts)


def dense_steps(al=0.8, levels=(1, 2, 3, 4, 5, 6), d=(1.0, -0.5, 0.0)):
    """max over V_n of |u_(n+1) - u_n| (dense L^alpha) against 2||x|| log(D'_(n+1)/D'_n) (quadrature)."""
    xn = max(abs(v - sum(d) / 3) for v in d)
    out = []
    prev = None
    Dp = {n: Dprime(n, str(al)) for n in levels}
    for n in levels:
        u, pos = bh.dense_u(n, al, d)
        if prev is not None:
            up, pp, m = prev
            diff = max(abs(u[pos[k]] - up[pp[k]]) for k in pp)
            bound = 2 * xn * float(mp.log(Dp[n] / Dp[m]))
            out.append((m, n, diff, bound))
        prev = (u, pos, n)
    return out


def deficit_errors(al="0.8", N=10, ns=(15, 20, 30, 40, 50), nref=100, dps=25):
    ref = bi.corner_deficit(al=al, n=nref, Ns=(N,), dps=dps)[0][1]
    out = []
    for n in ns:
        v = bi.corner_deficit(al=al, n=n, Ns=(N,), dps=dps)[0][1]
        out.append((n, v, ref))
    return out


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem BK: rate rho^-n for the uniform convergence u_n -> u-bar above log3/log5 (Conjecture BJ.7). Not thrust.")
    symbolic_checks()
    exact_checks()

    # (C1) BI.4 re-run: sup_(s>0) s^(1-gamma) r_inf <= G
    G0, _, _ = bi.sup_g()
    G = arb("0.6031")
    okb, _ = bi.bootstrap(G0, G)
    check(okb and G0.upper() < G.lower() and bool(arb(3) / 5 <= G), "BI.4 re-run: sup_(s>0) g <= G = 0.6031, so rbar(s) = min(1/s, G s^(gamma-1)) >= r_inf(s) >= r_k(s)")

    # (C2) J_a and the constants
    res = {}
    for al in ("0.75", "0.8", "0.9"):
        J, nc, hd, tl = J_bound(al, G)
        rho, D0, ec, C = constants(al, J)
        res[al] = (J, rho, D0, ec, C)
        print(f"   alpha={al}: J <= {float(J.upper()):.6f} ({nc} cells, head {hd:.2e}, tail {tl:.2e}); D'_0 = {D0.str(10)}; "
              f"D'_inf - D'_n <= {float(ec.upper()):.6f} rho^-n; C_a = 2J/((1-1/rho)D'_0) <= {float(C.upper()):.6f}")
    check(all(r[0].upper() < 1.3 for r in res.values()) and res["0.8"][4].upper() < 6.3,
          f"certified: sup_(V_n)|u_n - u-bar| <= C_a ||x|| rho^-n for every n: C_0.75 <= {float(res['0.75'][4].upper()):.4f}, "
          f"C_0.8 <= {float(res['0.8'][4].upper()):.4f}, C_0.9 <= {float(res['0.9'][4].upper()):.4f}")

    # (C3) with AZ's A_0.8: D'_inf and the asymptotic constant; explicit n-bounds
    J, rho, D0, ec, C = res["0.8"]
    a = arb("0.8")
    A = arb(A08[0]).union(arb(A08[1]))
    Dinf = arb.pi() * A / (2 * (arb.pi() * a).sin())
    Casy = 2 * ec / Dinf
    print(f"   alpha=0.8 with AZ's A_0.8: D'_inf = {Dinf.str(12)}; limsup rho^n sup|u_n-u-bar|/||x|| <= 2J/((1-1/rho)D'_inf) <= {float(Casy.upper()):.6f}")
    check(bool(Dinf.lower() > D0.upper()) and bool((Dinf - D0).upper() <= ec.upper()),
          "consistency: D'_0 < D'_inf and D'_inf - D'_0 <= J/(1-1/rho) (the n = 0 case of BK.4)")
    bnds = []
    for nn in (0, 10, 20, 40, 60):
        epsn = ec * rho ** (-nn)
        Dlo = Dinf - epsn
        if Dlo.lower() < D0.lower():
            Dlo = D0
        bnds.append((nn, float((2 * (Dinf / Dlo).log()).upper())))
    print("   certified sup_(V_n)|u_n - u-bar|/||x|| <= 2 log(D'_inf/max(D'_0, D'_inf - eps_n)): " + ", ".join(f"{b:.4g} (n={nn})" for nn, b in bnds))
    check(bnds[2][1] < 0.05 and bnds[3][1] < 1.2e-3 and bnds[4][1] < 2.6e-5, "certified at 0.8: sup|u_n - u-bar|/||x|| < 0.05 (n=20), < 1.2e-3 (n=40), < 2.6e-5 (n=60)")
    # derived from AZ's intercept as stated
    Cz = arb(C08)
    lim = -2 * Cz / A
    gapc = -arb.pi() * Cz / (2 * (arb.pi() * a).sin())
    print(f"   derived from AZ (C_0.8 as stated): rho^n (D'_inf - D'_n) -> {gapc.str(10)}, so limsup rho^n sup|u_n-u-bar|/||x|| <= -2C/A = {lim.str(10)}")
    check(bool(gapc.upper() <= ec.lower()) and abs(float(lim.mid()) - 1.4548906) < 1e-6,
          "consistency: AZ's -pi C/(2 sin pi a) = 4.6587 sits under the certified J/(1-1/rho); derived constant 2|nu_0.8| = 1.45489")

    # (D) floats
    Dp = {nn: Dprime(nn, "0.8") for nn in (0, 10, 20, 30)}
    okm = all(Dp[x] < Dp[y] for x, y in [(0, 10), (10, 20), (20, 30)]) and abs(Dp[0] - float(D0.mid())) < 1e-5
    gaps = [float((Dinf.mid() - arb(str(Dp[nn]))).mid()) * float(rho.mid()) ** nn for nn in (10, 20, 30)]
    print("   floats: rho^n (D'_inf - D'_n) = " + ", ".join(f"{g_:.5f}" for g_ in gaps) + " (n = 10, 20, 30)")
    check(okm and all(abs(g_ - 4.6587) < 2e-3 for g_ in gaps),
          "floats: D'_n increases (BK.2), D'_0 = 3^-a pi/sin(pi a), and rho^n(D'_inf - D'_n) ~ 4.6587 = AZ's -pi C/(2 sin pi a)")
    st = dense_steps()
    for m, nn, diff, bound in st:
        print(f"   dense alpha=0.8: max_(V_{m}) |u_{nn} - u_{m}| = {diff:.6f} <= 2||x|| log(D'_{nn}/D'_{m}) = {bound:.6f}")
    check(all(diff <= bound + 1e-9 for _, _, diff, bound in st),
          "floats: dense L^0.8 one-level differences at levels 1-6 sit under the BK.2 per-step bound (consistency)")
    de = deficit_errors()
    rhof = 5**0.8 / 3
    Dref = Dinf.mid()
    okd = True
    lines = []
    for nn, v, ref in de:
        err = float(ref - v)
        bnd = 1.5 * float(2 * ec.upper()) / 2 * rhof ** (-nn) / float(D0.lower())
        okd = okd and 0 < err <= bnd
        lines.append(f"n={nn}: delta^(100)-delta^(n) = {err:.4e}, x rho^n = {err * rhof**nn:.4f}, bound (3/2)eps_n/D'_0 = {bnd:.3e}")
    for ln in lines:
        print("   floats alpha=0.8, N=10: " + ln)
    check(okd and 0.4 < float(de[-1][2] - de[-1][1]) * rhof ** de[-1][0] < 0.9,
          "floats: the deficit n-error is positive, under the certified BK.3/BK.4 bound, and ~ 0.6 rho^-n (the order rho^-n looks attained; Conjecture BK.6)")

    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BK checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
