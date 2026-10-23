#!/usr/bin/env python3
"""Theorem BL: the rate rho^-n of Theorem BK is attained (sharpness half of Conjecture BK.6): at alpha = 0.8, for every
N >= 9 the data-free corner deficit delta_N^(n) strictly increases in n >= N, and delta_N - delta_N^(n) >= c_N rho^-n.

Setting (GASKET.md, Theorems BG-BK): p(y) = y(y+5), beta(z) = (z+6)/((z+2)(z+5)), alpha1(z) = (z+4)/((z+2)(z+5)),
e(y) = 2y/(2+y), rho = 5^a/3, gamma = log3/log5, R_k(tau) = (3/5)^k r_k(5^-k tau), D'_k = int tau^-a R_k,
mu_k = tau^-a R_k dtau / D'_k.  eta_N: eta_N(z) = alpha1(z) + beta(z) eta_(N-1)(p z), eta_0 = -1/2 (BJ.4).
G_k(tau) = 1 - eta_N(p^(k-N)(5^-k tau)) (k >= N), delta_N^(k) = <G_k>_(mu_k) (BG.3/BJ.4: the data-free corner deficit
(d_a - (u_k(nu_b^N) + u_k(nu_c^N))/2)/x_a), delta_N = lim_k delta_N^(k) (BG.3, BK.3).

Theorem BL (checked here; exit 1 on any failure):
  (1) exact (every N >= 1): 1 - eta_N(z) = z/(z+2) + beta(z)(1 - eta_(N-1)(p z)), i.e.
        1 - eta_N(z) = sum_(j<N) prod_(i<j) beta(p^i z) * p^j z/(p^j z + 2) + (3/2) prod_(i<N) beta(p^i z)   (all terms >= 0),
      eta_N is strictly decreasing on [0, inf) (alpha1, beta > 0 decreasing, eta_1 = 1/(2(z+5)) > 0, induction),
      eta_N(0) = 1 - (9/10)(3/5)^(N-1); so G_k is strictly increasing in tau with values in (0, 1].
  (2) covariance floor (every alpha in (gamma,1), N >= 1, k >= N, y0 > 0): with BK.1's transfer,
        delta^(k+1) - delta^(k) = <(omega_k - 1)(G_k - <G_k>)>/<omega_k> >= rho^-k [(g - Delta) L(y0) - Delta J(y0)] / D'_(k+1)
      whenever Delta >= delta^(k) and g <= 1 - eta_N(p(y0)) satisfy Delta < g; here
        L(y0) = int_(y0)^inf (omega_lo - 1) p^-a (3+p)^-1 p' dy,  omega_lo = (1+y/5)^a / (1 - y(y+4)/(3(y+2)(y+5))) <= omega_k,
        J(y0) = int_0^(y0) (omega_hi - 1) p^-a rbar(p) p' dy   (BK.4's ceiling restricted to y < y0).
  (3) floors for D'_k (AZ-free): D'_(k+1) - D'_k >= rho^-k L(0), so D'_N >= D'_0 + L(0)(1 - rho^-N)/(1 - rho^-1);
      ceiling D'_inf <= D'_0 + J_a/(1 - rho^-1) (BK.4).
  (4) ceilings for delta^(k), k >= N: (i) direct, Delta_N = I_N / D'_N(lower),
        I_N = int tau^-a min(1/3, G tau^(gamma-1)) H_N(tau) dtau, H_N = 1 - eta_N(e^(1/2) 5^-N tau) (tau <= 5^N), 1 beyond
      (orbit lemma p^j(s) <= e^(1/2) 5^j s for 5^j s <= 1, BI.1; R_k <= min(1/3, G tau^(gamma-1)), BI.4);
      (ii) N-uniform, Delta_N = (G Q rho^-N + (3/2)(3/5)^N I_0)/D'_N(lower), Q = int_0^inf t^(gamma-1-a) min(1, S(t)) dt,
        S(t) = sum_(l>=0) (3/5)^l phi(5^l t), phi(x) = e^(1/2) x/2 (x <= 1), 1 (x > 1); I_0 = int tau^-a min(1/3, G tau^(gamma-1)).
  (5) certified at alpha = 0.8 (python-flint arb): with y0 = 1/10, the bracket (g - Delta)L - Delta J is > 0 for N = 9..13
      by (4i) and for every N >= 14 by (4ii) (monotone in N); hence for every N >= 9 and n >= N:
        delta_N^(n+1) > delta_N^(n),   delta_N - delta_N^(n) >= c_N rho^-n,  c_N = [(g-Delta_N)L - Delta_N J]/(D'_inf(upper)(1-rho^-1)),
      and with symmetric data (x_b = x_c, ||x|| = |x_a|): sup_(V_n)|u_n - u-bar| >= c_N ||x|| rho^-n.  So the order rho^-n of
      BK.4 is exact (the sharpness half of Conjecture BK.6).  Existence of lim rho^n(delta_N - delta_N^(n)) stays open.

Floats (consistency only): delta_10^(n) increases in n, and the n-errors sit between c_10 rho^-n and BK's ceiling.
Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import sympy as sp
from flint import arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasket_fractional_uniform_limit as bi  # noqa: E402
import gasket_fractional_corner_sharpness as bj  # noqa: E402
import gasket_fractional_uniform_rate as bk  # noqa: E402

FAILS: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    z, y, X = sp.symbols("z y X", positive=True)
    al1 = (z + 4) / ((z + 2) * (z + 5))
    bz = (z + 6) / ((z + 2) * (z + 5))
    check(sp.simplify(1 - al1 - bz - z / (z + 2)) == 0,
          "1 - alpha1 - beta = z/(z+2), so 1 - eta_N(z) = z/(z+2) + beta(z)(1 - eta_(N-1)(pz)) (exact)")
    d1, d2 = sp.together(sp.diff(al1, z)), sp.together(sp.diff(bz, z))
    n1, n2 = sp.numer(d1), sp.numer(d2)
    check(all(c < 0 for c in sp.Poly(sp.expand(n1), z).coeffs()) and all(c < 0 for c in sp.Poly(sp.expand(n2), z).coeffs())
          and sp.Poly(sp.expand(sp.denom(d1)), z).coeffs()[0] > 0,
          f"alpha1' and beta' have numerators {sp.factor(n1)}, {sp.factor(n2)} with negative coefficients: alpha1, beta > 0 strictly decreasing")
    check(sp.simplify(al1 - bz / 2 - 1 / (2 * (z + 5))) == 0,
          "eta_1 = alpha1 - beta/2 = 1/(2(z+5)) > 0, strictly decreasing; induction: eta_N = alpha1 + beta*(eta_(N-1) o p), a sum of a "
          "decreasing function and a product of two positive decreasing ones, so every eta_N (N >= 1) is positive and strictly decreasing")
    # eta_N(0): eta_N(0) = 2/5 + (3/5) eta_(N-1)(0), eta_1(0) = 1/10
    e0 = sp.Rational(1, 10)
    ok = True
    for N in range(2, 12):
        e0 = sp.Rational(2, 5) + sp.Rational(3, 5) * e0
        ok = ok and e0 == 1 - sp.Rational(9, 10) * sp.Rational(3, 5) ** (N - 1)
    check(ok and al1.subs(z, 0) == sp.Rational(2, 5) and bz.subs(z, 0) == sp.Rational(3, 5),
          "eta_N(0) = 1 - (9/10)(3/5)^(N-1) (N = 1..11), eta_N(inf) = 0: G_k ranges in ((9/10)(3/5)^(N-1), 1)")
    # omega_lo
    beta_y = (y + 6) / ((y + 2) * (y + 5))
    check(sp.simplify((1 + 2 * y / 5) * sp.Rational(5, 3) * beta_y - (1 - y * (y + 4) / (3 * (y + 2) * (y + 5)))) == 0,
          "omega_k = (1+y/5)^a/((1+2y/5)(5/3)(beta - e r_k(p))) >= (1+y/5)^a/(1 - y(y+4)/(3(y+2)(y+5))) = omega_lo >= 1 (e r_k >= 0)")
    # covariance identity
    w, G, mw, mG, mwG = sp.symbols("w G mw mG mwG")
    lhs = mwG / mw - mG
    rhs = (mwG - mw * mG - mG + mG) / mw  # <(w-1)(G-<G>)> = <wG> - <w><G> - <G> + <G>
    check(sp.simplify(lhs - rhs) == 0,
          "<omega G>/<omega> - <G> = <(omega-1)(G-<G>)>/<omega>: the increment is a covariance against omega - 1 >= 0")
    # S(t) head bound: c t (3^(L+1)-1)/2 + (5/2)(3/5)^(L+1) <= (3c/2 + 5/2) t^(1-gamma) when 5^-(L+1) < t <= 5^-L
    c = sp.exp(sp.Rational(1, 2)) / 2
    check(c < 1 and sp.simplify(sum(sp.Rational(3, 5) ** l * 5**l for l in range(4)) - (3**4 - 1) / 2) == 0,
          "phi(x) = e^(1/2)x/2 <= 1 on [0,1]; sum_(l<=L) 3^l = (3^(L+1)-1)/2, so S(t) <= (3c/2 + 5/2) t^(1-gamma) on (0,1] (3^L <= t^-gamma, (3/5)^(L+1) < t^(1-gamma))")


# ----------------------------------------------------------------- (B) exact rationals
def one_minus_eta_series(z, N):
    """1 - eta_N(z) = sum_(j<N) prod_(i<j) beta(p^i z) w_j/(w_j+2) + (3/2) prod_(i<N) beta(p^i z); returns (total, truncated sum)."""
    S, P, w = 0 * z, 1 + 0 * z, z
    for _ in range(N):
        S += P * w / (w + 2)
        P *= (w + 6) / ((w + 2) * (w + 5))
        w = w * (w + 5)
    return S + P * 3 / 2, S


def exact_checks():
    zs = [Fraction(1, 1000), Fraction(1, 10), Fraction(1), Fraction(7, 3), Fraction(20)]
    ok_rec = ok_mono = ok_orbit = True
    for N in range(1, 7):
        vals = []
        for z in zs:
            tot, _ = one_minus_eta_series(z, N)
            ok_rec = ok_rec and tot == 1 - bj.eta_N_exact(z, N)
            vals.append(1 - bj.eta_N_exact(z, N))
        ok_mono = ok_mono and all(vals[i] < vals[i + 1] for i in range(len(vals) - 1)) and 0 < vals[0] and vals[-1] < 1
    for j in range(0, 6):
        for s in (Fraction(1, 5**j), Fraction(1, 3 * 5**j), Fraction(1, 1000 * 5**j)):
            w = s
            for _ in range(j):
                w = w * (w + 5)
            ok_orbit = ok_orbit and float(w) <= math.exp(0.5) * 5**j * float(s) * (1 - 1e-12)
    check(ok_rec, "exact rationals: the series for 1 - eta_N equals 1 - eta_N from BJ.4's recursion (N = 1..6, five rational z)")
    check(ok_mono, "exact rationals: 1 - eta_N strictly increasing in z with values in (0,1) at the sample points (N = 1..6)")
    check(ok_orbit, "orbit lemma (BI.1, re-checked): p^j(s) <= e^(1/2) 5^j s at 5^j s in {1, 1/3, 1/1000}, j = 0..5")


# ----------------------------------------------------------------- (C) certificates
def gam_arb():
    return arb(3).log() / arb(5).log()


def L_lower(al: str, y0: str, yhi="1e8", ratio="1.001", prec=96):
    """Lower bound of L(y0) = int_(y0)^inf (omega_lo - 1) p^-a (3+p)^-1 p' dy (cells on [y0, yhi], rest dropped)."""
    ctx.prec = prec
    a = arb(al)
    pts = [arb(y0)]
    rr = arb(ratio)
    while pts[-1].upper() < arb(yhi):
        pts.append(pts[-1] * rr)
    tot = arb(0)
    for i in range(len(pts) - 1):
        Y = pts[i].union(pts[i + 1])
        p = Y * (Y + 5)
        B = 1 - Y * (Y + 4) / (3 * (Y + 2) * (Y + 5))
        wl = (1 + Y / 5) ** a / B - 1
        f = wl * p ** (-a) / (3 + p) * (2 * Y + 5)
        lo = f.lower()
        if lo > 0:
            tot += arb(lo) * (pts[i + 1] - pts[i])
    return arb(tot.lower()), len(pts) - 1


def J_upper(al: str, G: arb, y0: str, ylo="1e-10", ratio="1.001", prec=96):
    """Upper bound of J(y0) = int_0^(y0) (omega_hi - 1) p^-a rbar(p) p' dy: BK.4's cell integrand on [ylo, y0] and BK's head."""
    ctx.prec = prec
    a = arb(al)
    gam = gam_arb()
    five3 = arb(5) / 3
    pts = [arb(ylo)]
    rr = arb(ratio)
    while pts[-1].upper() < arb(y0).lower():
        pts.append(pts[-1] * rr)
    pts[-1] = arb(y0)
    tot = arb(0)
    for i in range(len(pts) - 1):
        Y = pts[i].union(pts[i + 1])
        p = Y * (Y + 5)
        pp = 2 * Y + 5
        rp = bk._min_up(1 / p, G * p ** (gam - 1))
        ry = bk._min_up(1 / Y, G * Y ** (gam - 1))
        beta = (Y + 6) / ((Y + 2) * (Y + 5))
        e = 2 * Y / (2 + Y)
        num = (1 + Y / 5) ** a
        U2 = arb((arb(3) / 5 * num * ry * (3 + p) / (1 + 2 * Y / 5)).upper())
        den = (1 + 2 * Y / 5) * five3 * (beta - e * rp)
        w = U2
        if den.lower() > 0:
            w = bk._min_up(num / den, U2)
        wm1 = w - 1
        if wm1.upper() < 0:
            wm1 = arb(0)
        wt = arb((p ** (-a) * rp * pp).upper())
        tot += arb((wm1 * wt * (pts[i + 1] - pts[i])).upper())
    yl = arb(ylo)
    d0 = (4 + yl) * yl / 30 + (1 + 2 * yl / 5) * G * yl**gam
    c1 = a / 5 + (4 + yl) / 30
    c2 = (1 + 2 * yl / 5) * G
    pref = G * arb(5) ** (gam - 1 - a) * (5 + 2 * yl) / (1 - d0)
    head = pref * (c1 * yl ** (gam + 1 - a) / (gam + 1 - a) + c2 * yl ** (2 * gam - a) / (2 * gam - a))
    return arb((tot + head).upper()), len(pts) - 1


def H_up(tau: arb, N: int) -> arb:
    """Upper value of H_N(tau) = 1 - eta_N(e^(1/2) 5^-N tau) at a point (series, positive terms)."""
    z = arb(0.5).exp() * tau / arb(5) ** N
    tot, _ = one_minus_eta_series(z, N)
    return arb(tot.upper())


def I_direct(al: str, G: arb, N: int, th="1e-6", ratio="1.002", prec=80):
    """Upper bound of I_N = int tau^-a min(1/3, G tau^(gamma-1)) H_N(tau) dtau: weight decreasing, H_N increasing (BL.1)."""
    ctx.prec = prec
    a = arb(al)
    gam = gam_arb()
    T = arb(5) ** N
    t0 = arb(th)
    head = H_up(t0, N) * t0 ** (1 - a) / (3 * (1 - a))
    tail = G * T ** (gam - a) / (a - gam)
    rr = arb(ratio)
    tot = arb(0)
    nc = 0
    while t0.upper() < T.lower():
        t1 = t0 * rr
        if t1.upper() > T.lower():
            t1 = T
        wt = bk._min_up(t0 ** (-a) / 3, G * t0 ** (gam - 1 - a))
        tot += arb((wt * H_up(t1, N) * (t1 - t0)).upper())
        t0 = t1
        nc += 1
    return arb((head + tot + tail).upper()), nc


def S_inf(t: arb) -> arb:
    """S(t) = sum_l (3/5)^l phi(5^l t) at a point t in (0,1]: c t (3^(L+1)-1)/2 + (5/2)(3/5)^(L+1), 5^L t <= 1 < 5^(L+1) t."""
    c = arb(0.5).exp() / 2
    L = int(math.floor(-math.log(float(t.mid())) / math.log(5)))
    while (arb(5) ** (L + 1) * t).lower() <= 1:
        L += 1
    while L > 0 and (arb(5) ** L * t).upper() > 1:
        L -= 1
    assert (arb(5) ** L * t).upper() <= 1 and (arb(5) ** (L + 1) * t).lower() > 1
    return c * t * (arb(3) ** (L + 1) - 1) / 2 + arb(5) / 2 * (arb(3) / 5) ** (L + 1)


def Q_upper(al: str, th="1e-12", ratio="1.002", prec=80):
    """Upper bound of Q = int_0^inf t^(gamma-1-a) min(1, S(t)) dt: S increasing; head K_S t^(1-gamma); [1, inf) exact."""
    ctx.prec = prec
    a = arb(al)
    gam = gam_arb()
    c = arb(0.5).exp() / 2
    KS = 3 * c / 2 + arb(5) / 2
    t0 = arb(th)
    head = KS * t0 ** (1 - a) / (1 - a)
    tail = 1 / (a - gam)
    rr = arb(ratio)
    tot = arb(0)
    while t0.upper() < 1:
        t1 = t0 * rr
        if t1.upper() > 1:
            t1 = arb(1)
        s = S_inf(t1)
        m = arb(1) if s.lower() >= 1 else arb(s.upper())
        tot += arb((m * bj.Wint(a, gam, t0, t1)).upper())
        t0 = t1
    return arb((head + tot + tail).upper())


# ----------------------------------------------------------------- (D) floats
def float_deficits(N=10, ns=(10, 11, 12, 14, 16, 20, 25, 30, 40), nref=100):
    ref = bi.corner_deficit(al="0.8", n=nref, Ns=(N,), dps=25)[0][1]
    vals = [(n, bi.corner_deficit(al="0.8", n=n, Ns=(N,), dps=25)[0][1]) for n in ns]
    return vals, ref


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem BL: the rate rho^-n of Theorem BK is attained at the corner neighbours (alpha = 0.8, N >= 9). Not thrust.")
    symbolic_checks()
    exact_checks()

    al = "0.8"
    a = arb(al)
    gam = gam_arb()
    rho = arb(5) ** a / 3
    # BI.4 re-run (G) and BK.4 re-run (J_0.8)
    G0, _, _ = bi.sup_g()
    G = arb("0.6031")
    okb, _ = bi.bootstrap(G0, G)
    check(okb and G0.upper() < G.lower(), "BI.4 re-run: sup_(s>0) s^(1-gamma) r_inf <= G = 0.6031, so R_k <= min(1/3, G tau^(gamma-1)) and r_k <= rbar")
    Jtot, _, _, _ = bk.J_bound(al, G)
    _, D0, ec, _ = bk.constants(al, Jtot)
    Dup = arb((D0 + ec).upper())
    print(f"   BK.4 re-run: J_0.8 <= {float(Jtot.upper()):.6f}; D'_0 = {D0.str(10)}; D'_inf <= D'_0 + J/(1-1/rho) <= {float(Dup.upper()):.6f}")
    check(Jtot.upper() < 1.1882, "BK.4 re-run: J_0.8 <= 1.188154")

    # L(0), L(y0), J(y0)
    y0 = "0.1"
    L0, nL0 = L_lower(al, "1e-10")
    Ly, nLy = L_lower(al, y0)
    Jy, nJy = J_upper(al, G, y0)
    print(f"   L(0) >= {float(L0.lower()):.6f} ({nL0} cells), L(1/10) >= {float(Ly.lower()):.6f} ({nLy} cells), J(1/10) <= {float(Jy.upper()):.6f} ({nJy} cells)")
    check(L0.lower() > 0.341 and Ly.lower() > 0.3345 and Jy.upper() < 0.1786, "certified: L(0) > 0.341, L(1/10) > 0.3345, J(1/10) < 0.1786")
    Dlow = lambda N: arb((D0 + L0 * (1 - rho ** (-N)) / (1 - 1 / rho)).lower())  # noqa: E731
    print(f"   D'_N >= D'_0 + L(0)(1-rho^-N)/(1-rho^-1): {float(Dlow(9).lower()):.6f} (N=9), {float(Dlow(14).lower()):.6f} (N=14), "
          f"D'_inf >= {float((D0 + L0 / (1 - 1 / rho)).lower()):.6f} (AZ-free)")
    check(bool(Dlow(9) > D0) and bool((D0 + L0 / (1 - 1 / rho)).lower() < 6.40416) and bool(Dup > 6.40417),
          "consistency: AZ's D'_inf(0.8) = 6.40416 lies between the AZ-free floor and the BK ceiling")

    def bracket(g, Dl):
        return (g - Dl) * Ly - Dl * Jy

    p0 = arb(y0) * (arb(y0) + 5)
    rows = []
    ok_all = True
    # (4i) direct, N = 9..13
    for N in range(9, 14):
        IN, nc = I_direct(al, G, N)
        Dl = arb((IN / Dlow(N)).upper())
        tot, _ = one_minus_eta_series(p0, N)
        g = arb(tot.lower())
        br = bracket(g, Dl)
        cN = arb((br / (Dup * (1 - 1 / rho))).lower())
        ok_all = ok_all and bool(Dl < g) and br.lower() > 0
        rows.append((N, float(IN.upper()) * float(rho.mid()) ** N, float(Dl.upper()), float(g.lower()), float(br.lower()), float(cN.lower())))
    # (4ii) uniform, N >= 14
    Q = Q_upper(al)
    I0 = bj.I_closed(al, arb(1), G)
    N2 = 14
    Du = arb(((G * Q * rho ** (-N2) + arb(3) / 2 * (arb(3) / 5) ** N2 * I0) / Dlow(N2)).upper())
    _, trunc = one_minus_eta_series(p0, N2)
    gu = arb(trunc.lower())
    bru = bracket(gu, Du)
    cu = arb((bru / (Dup * (1 - 1 / rho))).lower())
    print(f"   Q <= {float(Q.upper()):.6f}, G Q <= {float((G * Q).upper()):.6f}, I_0 <= {float(I0.upper()):.6f}")
    for N, rI, Dl, g, br, cN in rows:
        print(f"   N={N} (direct): rho^N I_N <= {rI:.4f}, Delta_N <= {Dl:.6f} < g = 1-eta_N(p(1/10)) >= {g:.6f}; bracket >= {br:.3e}; c_N >= {cN:.3e}")
    print(f"   N>={N2} (uniform): Delta_N <= {float(Du.upper()):.6f} < g >= {float(gu.lower()):.6f}; bracket >= {float(bru.lower()):.3e}; c_N >= {float(cu.lower()):.3e}")
    ok_all = ok_all and bool(Du < gu) and bru.lower() > 0
    # asymptotic constant: Delta_N -> 0 as N -> inf, so c_N -> g_inf L/(D'_inf(1-1/rho)), g_inf = 1 - eta(p(1/10)) >= truncated series
    _, tr40 = one_minus_eta_series(p0, 40)
    casym = arb((tr40 * Ly / (Dup * (1 - 1 / rho))).lower())
    Nbig = 40
    Db = arb(((G * Q * rho ** (-Nbig) + arb(3) / 2 * (arb(3) / 5) ** Nbig * I0) / Dlow(Nbig)).upper())
    cbig = arb((bracket(arb(tr40.lower()), Db) / (Dup * (1 - 1 / rho))).lower())
    print(f"   c_N >= {float(cbig.lower()):.4f} for N >= {Nbig}; liminf over N of the floor >= {float(casym.lower()):.4f}")
    check(ok_all and float(cu.lower()) > 0,
          f"certified (alpha=0.8, y0=1/10): (g - Delta_N) L - Delta_N J > 0 for N = 9..13 (direct) and every N >= 14 (uniform, monotone in N); "
          f"so delta_N^(n) strictly increases in n >= N and delta_N - delta_N^(n) >= c_N rho^-n with c_9 >= {rows[0][5]:.2e}, "
          f"c_14.. >= {float(cu.lower()):.2e}, c_N >= {float(cbig.lower()):.4f} for N >= 40")
    # ratio to BK's ceiling
    up = arb(3) / 2 * ec / D0  # |delta - delta^(n)| <= (3/2) log(D'_inf/D'_n) <= (3/2)(D'_inf-D'_n)/D'_0
    print(f"   two-sided: c_N rho^-n <= delta_N - delta_N^(n) <= {float(up.upper()):.4f} rho^-n (BK.3-BK.4)")

    # (D) floats
    vals, ref = float_deficits()
    rf = 5**0.8 / 3
    inc = all(vals[i][1] < vals[i + 1][1] for i in range(len(vals) - 1)) and vals[-1][1] < ref
    errs = [(n, float(ref - v) * rf**n) for n, v in vals]
    print("   floats alpha=0.8, N=10: rho^n (delta^(100) - delta^(n)) = " + ", ".join(f"{e:.4f} (n={n})" for n, e in errs))
    check(inc and all(rows[1][5] <= e <= float(up.upper()) for _, e in errs),
          "floats: delta_10^(n) increases in n (BL.5) and rho^n times the n-error lies in [c_10, BK's ceiling] (consistency)")

    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BL checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
