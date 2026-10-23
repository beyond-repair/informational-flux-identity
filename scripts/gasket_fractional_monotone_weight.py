#!/usr/bin/env python3
"""Theorem BN: BK's one-level weight omega_k is increasing (for every alpha in [1/17,1) and every level k), so the corner
deficit delta_N^(n) increases in n for every N >= 1, and the rate rho^-n is attained for every N >= 1 at every alpha in
(gamma, 1) (sharpness half of Conjecture BL.6 for all N >= 1; removes BL/BM's thresholds N_0(alpha)).

Setting (GASKET.md, Theorems AQ, BG-BM): r_k(s) = (1/2) v^T (L_k + s)^-1 v (AQ), p(y) = y(y+5), beta(y) = (y+6)/((y+2)(y+5)),
e(y) = 2y/(2+y), B(y) = 1 - y(y+4)/(3(y+2)(y+5)), C(y) = (5/3)(1+2y/5) e(y), rho = 5^a/3, gamma = log3/log5.
BK.1's weight omega_k(y) = (1+y/5)^a / d_k(y), d_k(y) = (1+2y/5)(5/3)(beta(y) - e(y) r_k(p(y))) = B(y) - C(y) r_k(p(y)) > 0.
G_k = 1 - eta_N(p^(k-N+1)(y)) in BK.1's coordinate (tau' = 5^k p(y)), increasing (BL.1); delta_N^(k) = <G_k>_(mu_k).

Theorem BN (checked here; exit 1 on any failure):
  (1) [spectral lemma] r_k is a Stieltjes function of unit mass: r_k(s) = sum_i w_i/(lambda_i + s), w_i >= 0, sum w_i = |v|^2/2 = 1,
      lambda_i >= 0 (L_k symmetric PSD).  Hence r_k(s) <= 1/s and 0 <= -r_k'(s) <= r_k(s)/s for every k, s > 0.
  (2) [monotone weight] d_k' - a d_k/(5+y) <= B' - aB/(5+y) + r K,  K = K_0 + aC/(5+y),  K_0 = -C' + C p'/p
      = 2y(2y^2+10y+15)/(3(y+2)^2(y+5)) >= 0 (sympy), so with r = r_k(p) <= 1/p:
        d_k' - a d_k/(5+y) <= Psi_a(y) = -P_a(y)/(3(y+2)^2(y+5)^2),  P_a(y) = 2a y^3 + (17a-1) y^2 + 46a y + 40a + 10,
      and P_a > 0 on y > 0 for a >= 1/17 (all coefficients >= 0).  So (log omega_k)' = a/(5+y) - d_k'/d_k > 0:
      omega_k is strictly increasing on (0, inf) for every k >= 0 and every alpha in [1/17, 1).  (Was floats-only in BL.)
  (3) [monotone deficit] Chebyshev's covariance inequality: delta_N^(k+1) - delta_N^(k) = Cov_(mu_k)(omega_k, G_k)/<omega_k> > 0
      for every alpha in [1/17,1), N >= 1, k >= N (both functions strictly increasing in y, mu_k has a positive density).
      Monotonicity half of Conjecture BL.6, for every N >= 1 and every alpha in [1/17,1) (above and below gamma).
  (4) [rate floor, every N >= 1] Cov = int int_(y'<y) (omega(y)-omega(y'))(G(y)-G(y')) dmu dmu, all terms >= 0.  Keep
      y' in A = {tau' <= tau0} (tau0 <= p(y0), so y' <= y0, omega_k(y') <= omega_k(y0) <= W0 = BK.4's omega_hi(y0), and
      G_k(y') <= Gbar_N = 1 - eta_N(e^(1/2) 5^-N tau0) by the orbit lemma, 5^-N tau0 <= 1) and y >= y0 (omega_k >= omega_lo,
      G_k >= H_N(y) = 1 - eta_N(p(y)), r_k(p) >= 1/(3+p)); mu_k(A) >= M/D'_k, M = int_0^tau0 t^-a/(3+t) (R_k >= 1/(3+t), BI.3).
      So delta^(k+1) - delta^(k) >= rho^-k M Lambda_N / (D'_k D'_(k+1)),
        Lambda_N = int_(y0)^inf (omega_lo - W0)^+ (H_N - Gbar_N)^+ p^-a p'/(3+p) dy,
      and for alpha in (gamma,1): delta_N - delta_N^(n) >= c_N rho^-n, c_N = M Lambda_N / (D'_inf^2 (1 - rho^-1)) > 0 for EVERY N >= 1
      (Lambda_N > 0: omega_lo -> inf, H_N -> 1 > Gbar_N).  For N >= N1: H_N >= the series truncated at j < N1 and
      Gbar_N <= K_S (5^-N1 tau0)^(1-gamma) + (3/2)(3/5)^N1 (BL.4(ii)'s pointwise bound), one certificate for all N >= N1.
  (5) certified (python-flint arb, outward) at alpha = 0.7, 0.75, 0.8, 0.85, 0.9, 0.95 with y0 = 1/10, tau0 = p(1/10) = 51/100:
      c_N > 0 for N = 1..19 (direct) and every N >= 20 (uniform), so with BL/BM's large-N result the two-sided bound
      c_N rho^-n <= delta_N - delta_N^(n) <= C rho^-n holds for every N >= 1 at these exponents.

Floats (consistency only): omega_k increasing on a grid (k <= 12, three exponents incl. 0.45 < gamma); delta_N^(n) for N = 1..8
increasing in n at 0.8 with gaps above the certified floors.  Existence of lim rho^n (delta_N - delta_N^(n)) stays open.
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
import gasket_fractional_uniform_rate as bk  # noqa: E402
import gasket_fractional_rate_sharpness as bl  # noqa: E402
import mesh_maxwell_limit as ml  # noqa: E402

FAILS: list[str] = []
ALPHAS = ["0.7", "0.75", "0.8", "0.85", "0.9", "0.95"]
Y0 = "0.1"
N1 = 20


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    y, a = sp.symbols("y a", positive=True)
    p = y * (y + 5)
    pp = sp.diff(p, y)
    B = 1 - y * (y + 4) / (3 * (y + 2) * (y + 5))
    beta = (y + 6) / ((y + 2) * (y + 5))
    e = 2 * y / (2 + y)
    C = sp.Rational(5, 3) * (1 + 2 * y / 5) * e
    check(sp.simplify((1 + 2 * y / 5) * sp.Rational(5, 3) * beta - B) == 0,
          "d_k = (1+2y/5)(5/3)(beta - e r_k(p)) = B - C r_k(p), C = (5/3)(1+2y/5) e >= 0")
    # chain rule identity with an abstract r(.)
    r = sp.Function("r")
    d = B - C * r(p)
    R0, R1 = sp.symbols("R0 R1")
    lhs_s = sp.diff(d, y) - a * d / (5 + y)
    lhs_s = lhs_s.replace(lambda ex: isinstance(ex, sp.Subs), lambda ex: R1).subs(r(p), R0)
    rhs = sp.diff(B, y) - a * B / (5 + y) + R0 * (-sp.diff(C, y) + a * C / (5 + y)) - C * pp * R1
    ok_chain = sp.simplify(lhs_s - rhs) == 0
    check(ok_chain, "d_k' - a d_k/(5+y) = B' - aB/(5+y) + r(-C' + aC/(5+y)) + C p' (-r'(p)), r = r_k(p(y)) (sympy)")
    K0 = sp.factor(sp.together(-sp.diff(C, y) + C * pp / p))
    n0, d0 = sp.fraction(K0)
    check(sp.simplify(K0 - 2 * y * (2 * y**2 + 10 * y + 15) / (3 * (y + 2) ** 2 * (y + 5))) == 0,
          f"K_0 = -C' + C p'/p = {K0} >= 0, so K = K_0 + aC/(5+y) >= 0 and (with -r' <= r/p, r <= 1/p) the r-terms are <= K/p")
    T0 = sp.factor(sp.together(sp.diff(B, y) + K0 / p))
    BC = sp.factor(sp.together(B - C / p))
    check(sp.simplify(T0 - (y**2 - 10) / (3 * (y + 2) ** 2 * (y + 5) ** 2)) == 0
          and sp.simplify(BC - (y + 4) * (2 * y + 5) / (3 * (y + 2) * (y + 5))) == 0,
          f"B' + K_0/p = {T0} (positive for y > sqrt(10): the a = 0 part alone does NOT close), B - C/p = {BC} > 0")
    Psi = sp.together(T0 - a * BC / (5 + y))
    Pa = 2 * a * y**3 + (17 * a - 1) * y**2 + 46 * a * y + 40 * a + 10
    check(sp.simplify(Psi + Pa / (3 * (y + 2) ** 2 * (y + 5) ** 2)) == 0,
          "Psi_a = B' - aB/(5+y) + K/p = -P_a/(3(y+2)^2(y+5)^2), P_a = 2a y^3 + (17a-1) y^2 + 46a y + 40a + 10 (sympy)")
    aa = sp.Rational(1, 17)
    check(all(c >= 0 for c in sp.Poly(Pa.subs(a, aa), y).all_coeffs()) and Pa.subs(a, aa).subs(y, 0) > 0,
          "for a >= 1/17 every coefficient of P_a is >= 0 and P_a(0) = 40a+10 > 0: Psi_a < 0 on y > 0, so (log omega_k)' > 0")
    # witness that the bound is genuinely restricted: tiny a fails for this majorant
    a_bad = sp.Rational(1, 100)
    ys = [sp.Integer(k) for k in range(1, 200)]
    neg = [yy for yy in ys if Pa.subs({a: a_bad, y: yy}) < 0]
    check(len(neg) > 0, f"remark: at a = 1/100 the majorant P_a is negative at y = {neg[0]}..{neg[-1]} (integers): the method (not "
                        "necessarily omega_k) needs a bounded away from 0; a >= 1/17 covers every alpha >= gamma")
    t, t0 = sp.symbols("t t0", positive=True)
    check(sp.simplify(sp.integrate(t ** (-a) / (3 + t0), (t, 0, t0), conds="none") - t0 ** (1 - a) / ((1 - a) * (3 + t0))) == 0,
          "M = int_0^tau0 t^-a/(3+t) dt >= int_0^tau0 t^-a/(3+tau0) dt = tau0^(1-a)/((1-a)(3+tau0))")
    x, x2, w1, w2, g1, g2 = sp.symbols("x x2 w1 w2 g1 g2")
    # Chebyshev pair identity for a two-point measure (the general case is the same algebra integrated)
    m1, m2 = sp.symbols("m1 m2", positive=True)
    cov = (m1 * w1 * g1 + m2 * w2 * g2) / (m1 + m2) - (m1 * w1 + m2 * w2) * (m1 * g1 + m2 * g2) / (m1 + m2) ** 2
    check(sp.simplify(cov - m1 * m2 * (w2 - w1) * (g2 - g1) / (m1 + m2) ** 2) == 0,
          "Cov(omega,G) = (1/2) int int (omega(y)-omega(y'))(G(y)-G(y')) dmu dmu (checked on two atoms): >= 0 for co-monotone pairs")


# ----------------------------------------------------------------- (B) exact
def r_rec(k, s):
    if k == 0:
        return 1 / (3 + s)
    x = r_rec(k - 1, s * (s + 5))
    return (2 + s) * (5 + s) * x / ((6 + s) - 2 * s * (5 + s) * x)


def exact_checks():
    s = sp.symbols("s", positive=True)
    ok_graph = True
    for n in range(0, 4):
        P, adj, Cn, _, _ = ml.corner_frame(n) if n > 0 else (None, {0: [1, 2], 1: [0, 2], 2: [0, 1]}, [0, 1, 2], None, None)
        if n == 0:
            adj = [[1, 2], [0, 2], [0, 1]]
        Lm = sp.Matrix(ml.laplacian(adj).astype(int).tolist())
        a_, b_ = Cn[0], Cn[1]
        v = sp.zeros(Lm.shape[0], 1)
        v[a_], v[b_] = 1, -1
        for sv in (sp.Rational(1, 3), sp.Integer(2), sp.Rational(7, 2)):
            val = (v.T * (Lm + sv * sp.eye(Lm.shape[0])).LUsolve(v))[0] / 2
            ok_graph = ok_graph and sp.nsimplify(val) == sp.nsimplify(r_rec(n, sv))
        ok_graph = ok_graph and Lm.is_symmetric()
    check(ok_graph, "exact rationals: (1/2) v^T (L_n + s)^-1 v = AQ's recursion r_n(s) at s = 1/3, 2, 7/2 for n = 0..3, L_n symmetric "
                    "(so r_n is Stieltjes with unit mass (1/2)|v|^2 = 1, spectrum of L_n in [0, inf))")
    ok_st = True
    for k in range(0, 4):
        rk = sp.cancel(r_rec(k, s))
        num, den = sp.fraction(rk)
        dp = sp.Poly(den, s)
        roots = sp.Poly(den, s).real_roots()
        ok_st = ok_st and len(roots) == dp.degree() and all(rt < 0 for rt in roots)
        res = [sp.N(num.subs(s, rt) / sp.diff(den, s).subs(s, rt), 30) for rt in roots]
        ok_st = ok_st and all(rv > 0 for rv in res) and abs(sum(res) - 1) < sp.Float("1e-25")
        drk = sp.diff(rk, s)
        for sv in (sp.Rational(1, 1000), sp.Rational(1, 10), sp.Integer(1), sp.Integer(30)):
            ok_st = ok_st and rk.subs(s, sv) <= 1 / sv and 0 <= -drk.subs(s, sv) <= rk.subs(s, sv) / sv
    check(ok_st, "exact: r_k (k = 0..3) has only simple negative real poles with positive residues summing to 1, and r_k <= 1/s, "
                 "0 <= -r_k' <= r_k/s at s = 1/1000, 1/10, 1, 30")


# ----------------------------------------------------------------- (C) certificates
def series_lo(z: arb, N: int, cap=arb("1e30")) -> arb:
    """Lower bound of 1 - eta_N(z) (BL.1's positive series; stops, dropping positive terms, once p^j z > cap)."""
    S, P, w = arb(0), arb(1), z
    for j in range(N):
        S += P * w / (w + 2)
        P *= (w + 6) / ((w + 2) * (w + 5))
        if w.lower() > cap:
            return arb(S.lower())
        w = w * (w + 5)
    return arb((S + P * 3 / 2).lower())


def series_trunc_lo(z: arb, J: int, cap=arb("1e30")) -> arb:
    """Lower bound of sum_(j<J) prod beta * w_j/(w_j+2): a floor for 1 - eta_N(z) for every N >= J."""
    S, P, w = arb(0), arb(1), z
    for j in range(J):
        S += P * w / (w + 2)
        P *= (w + 6) / ((w + 2) * (w + 5))
        if w.lower() > cap:
            break
        w = w * (w + 5)
    return arb(S.lower())


def omega_hi_point(al: str, G: arb, y: arb) -> arb:
    """BK.4's k-uniform ceiling of omega_k at the point y (two forms, take the smaller)."""
    a = arb(al)
    gam = bl.gam_arb()
    p = y * (y + 5)
    rp = bk._min_up(1 / p, G * p ** (gam - 1))
    ry = bk._min_up(1 / y, G * y ** (gam - 1))
    beta = (y + 6) / ((y + 2) * (y + 5))
    e = 2 * y / (2 + y)
    num = (1 + y / 5) ** a
    U2 = arb((arb(3) / 5 * num * ry * (3 + p) / (1 + 2 * y / 5)).upper())
    den = (1 + 2 * y / 5) * arb(5) / 3 * (beta - e * rp)
    if den.lower() > 0:
        return bk._min_up(num / den, U2)
    return U2


def Lambda_lower(al: str, W0: arb, Gbar: arb, Hfun, y0=Y0, yhi="1e6", ratio="1.01", prec=80):
    """Lower bound of int_(y0)^inf (omega_lo - W0)^+ (H(y) - Gbar)^+ p^-a p'/(3+p) dy; H increasing (left endpoint); beyond yhi dropped."""
    ctx.prec = prec
    a = arb(al)
    t0 = arb(y0)
    rr = arb(ratio)
    tot = arb(0)
    while t0.upper() < arb(yhi):
        t1 = t0 * rr
        Y = t0.union(t1)
        p = Y * (Y + 5)
        Bc = 1 - Y * (Y + 4) / (3 * (Y + 2) * (Y + 5))
        wl = ((1 + Y / 5) ** a / Bc - W0).lower()
        hl = (Hfun(t0 * (t0 + 5)) - Gbar).lower()
        if wl > 0 and hl > 0:
            wt = (p ** (-a) * (2 * Y + 5) / (3 + p)).lower()
            tot += arb(wl) * arb(hl) * arb(wt) * (t1 - t0)
        t0 = t1
    return arb(tot.lower())


def certify(al: str, G: arb):
    ctx.prec = 96
    a = arb(al)
    gam = bl.gam_arb()
    rho = arb(5) ** a / 3
    Jtot, _, _, _ = bk.J_bound(al, G)
    ctx.prec = 96
    _, D0, ec, _ = bk.constants(al, Jtot)
    Dup = arb((D0 + ec).upper())
    y0 = arb(Y0)
    tau0 = y0 * (y0 + 5)
    W0 = omega_hi_point(al, G, y0)
    M = arb((tau0 ** (1 - a) / ((1 - a) * (3 + tau0))).lower())
    den = Dup**2 * (1 - 1 / rho)
    rows = []
    ok = bool(W0.upper() < 2) and M.lower() > 0
    for N in range(1, N1):
        Gbar = arb(bl.one_minus_eta_series(arb(0.5).exp() * arb(5) ** (-N) * tau0, N)[0].upper())
        Lam = Lambda_lower(al, W0, Gbar, lambda z, N=N: series_lo(z, N))
        ctx.prec = 96
        cN = arb((M * Lam / den).lower())
        ok = ok and cN.lower() > 0
        rows.append((N, float(Gbar.upper()), float(Lam.lower()), float(cN.lower())))
    c = arb(0.5).exp() / 2
    KS = 3 * c / 2 + arb(5) / 2
    zt = arb(5) ** (-N1) * tau0
    Gbu = arb((KS * zt ** (1 - gam) + arb(3) / 2 * (arb(3) / 5) ** N1).upper())
    Lu = Lambda_lower(al, W0, Gbu, lambda z: series_trunc_lo(z, N1))
    ctx.prec = 96
    cu = arb((M * Lu / den).lower())
    ok = ok and cu.lower() > 0 and bool(zt < 1)
    print(f"   alpha={al}: rho={float(rho.mid()):.6f}, J_a<={float(Jtot.upper()):.6f}, D'_inf<={float(Dup.upper()):.6f}, "
          f"W0=omega_hi(1/10)<={float(W0.upper()):.6f}, M>={float(M.lower()):.6f}")
    print("      direct c_N >= " + ", ".join(f"{cN:.3e} (N={N})" for N, _, _, cN in rows))
    print(f"      N>={N1} (uniform): Gbar<={float(Gbu.upper()):.3e}, Lambda>={float(Lu.lower()):.6f}, c_N>={float(cu.lower()):.3e}")
    check(ok, f"certified alpha={al}: c_N > 0 for N = 1..{N1 - 1} (direct) and every N >= {N1} (uniform): "
              f"delta_N - delta_N^(n) >= c_N rho^-n for every N >= 1 and n >= N")
    return {"rows": rows, "cu": float(cu.lower()), "rho": float(rho.mid())}


# ----------------------------------------------------------------- (D) floats
def float_omega_monotone():
    """mpmath floats (40 digits; omega_k - 1 is O(y^gamma) at small y, so double precision is not enough)."""
    from mpmath import mp, mpf
    mp.dps = 40
    ok = True
    ys = [mpf(10) ** (mpf(j) / 40 - 5) for j in range(0, 321)]  # 1e-5 .. 1e3
    for al in (mpf("0.45"), mpf("0.8"), mpf("0.95")):
        for k in (0, 1, 2, 4, 8, 12):
            vals = []
            for yv in ys:
                p = yv * (yv + 5)
                rk = ml.rec_all(k, p)[0]
                beta = (yv + 6) / ((yv + 2) * (yv + 5))
                e = 2 * yv / (2 + yv)
                vals.append((1 + yv / 5) ** al / ((1 + 2 * yv / 5) * mpf(5) / 3 * (beta - e * rk)))
            ok = ok and all(vals[i] < vals[i + 1] for i in range(len(vals) - 1)) and vals[0] >= 1
    return ok


def float_deficits(ns=(8, 12, 20)):
    out = {}
    for n in ns:
        out[n] = {N: v for N, v, _ in bi.corner_deficit(al="0.8", n=n, Ns=range(1, 9), dps=20)}
    return out


def main() -> int:
    t0 = time.time()
    print("Theorem BN: omega_k is increasing; delta_N^(n) increases in n and the rate rho^-n is attained for every N >= 1. Not thrust.")
    symbolic_checks()
    exact_checks()
    bl.FAILS.clear()
    bl.symbolic_checks()
    bl.exact_checks()
    check(not bl.FAILS, "Theorem BL's sympy and exact-rational checks (BL.1: 1 - eta_N increasing, positive series) re-run")
    G0, _, _ = bi.sup_g()
    G = arb("0.6031")
    okb, _ = bi.bootstrap(G0, G)
    check(okb and G0.upper() < G.lower(), "BI.4 re-run: sup s^(1-gamma) r_inf <= G = 0.6031 (r_k <= min(1/s, G s^(gamma-1)))")
    res = {al: certify(al, G) for al in ALPHAS}
    ok_om = float_omega_monotone()
    check(ok_om, "floats: omega_k increasing in y on 321 log-grid points of [1e-5,1e3] (mpmath, 40 digits) for k = 0,1,2,4,8,12 at alpha = 0.45, 0.8, 0.95 "
                 "(consistency; 0.45 < gamma is covered by BN.2 since 0.45 >= 1/17)")
    dd = float_deficits()
    rows8 = {N: cN for N, _, _, cN in res["0.8"]["rows"]}
    rf = res["0.8"]["rho"]
    inc = all(dd[8][N] < dd[12][N] < dd[20][N] for N in range(1, 9))
    gap_ok = all(float(dd[20][N] - dd[8][N]) >= rows8[N] * (rf ** -8 - rf ** -20) for N in range(1, 9))
    print("   floats alpha=0.8: delta_N^(8), delta_N^(20) = " + ", ".join(f"N={N}: {float(dd[8][N]):.4f}->{float(dd[20][N]):.4f}" for N in range(1, 9)))
    check(inc and gap_ok, "floats: delta_N^(n)(0.8) increases over n = 8, 12, 20 for N = 1..8 and delta^(20)-delta^(8) >= c_N(rho^-8-rho^-20) "
                          "(consistency with BN.3-BN.4)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BN checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
