#!/usr/bin/env python3
"""Theorem BO: the exact exponent threshold for Theorem BN's monotone weight.  BK's one-level weight omega_k is increasing in y
at every level k iff alpha >= a_k* := sup_(y>0) Phi_k(y), Phi_k = (5+y) d_k'/d_k, and the level-uniform threshold
a* = sup_k a_k* lies in [0.016778880172852, 0.0167788801732] (BN proved alpha >= 1/17 = 0.0588...).  Kept Failure BO.1: for
alpha < 0.016778880172852 the weight omega_k is NOT monotone at any level k >= 3 (explicit witness y* = 13.778505644, exact
rationals), so BN.2 as a statement about all levels fails below a*, and BN's Chebyshev route (BN.3) is unavailable there.

Setting (GASKET.md, Theorems AQ, BK, BN): r_k(s) = (1/2) v^T (L_k + s)^-1 v, v = e_a - e_b, p(y) = y(y+5), e(y) = 2y/(2+y),
B(y) = 1 - y(y+4)/(3(y+2)(y+5)), C(y) = (5/3)(1+2y/5) e(y), d_k(y) = B(y) - C(y) r_k(p(y)) > 0, omega_k = (1+y/5)^a / d_k.

Theorem BO (checked here; exit 1 on any failure):
  (1) [reduction] (log omega_k)' = (a - Phi_k(y))/(5+y), so omega_k is strictly increasing on (0,inf) iff a >= a_k* = sup Phi_k
      (if a < a_k*, omega_k strictly decreases on an interval).  The condition is linear in a, so each witness covers an interval.
  (2) [moment bounds, every k] r_k(s) = sum_i w_i/(lambda_i+s), w_i >= 0, sum w_i = 1 (BN.1), moments m_j = sum w_i lambda_i^j
      = (1/2) v^T L_k^j v = (L_k^j)_(aa) for j < 2^k.  Exact geometric remainders give, for every J >= 1 and s > 0,
        r_k(s) - sum_(j<J) (-1)^j m_j s^-(j+1)  and  -r_k'(s) - sum_(j<J) (j+1)(-1)^j m_j s^-(j+2)  of sign (-1)^J
      (odd J: upper bounds, even J: lower bounds).  Spectrum in [0, 8] (Gershgorin, max degree 4) and m_1 >= 2 give the
      k-uniform chord bounds r_k <= 1/s - 2/(s(s+8)),  -r_k' <= r_k/s - 2/(s(s+8)^2).
  (3) [locality] m_j(k) = m_j(k+1) whenever j < 2^k (walks of length j from corner a stay at distance <= j/2 < 2^k inside the
      corner subcell, whose degrees agree with V_k there, and cannot reach b).  Hence m_0..m_7 = 1,2,6,22,92,418,1994,9780 for
      EVERY k >= 3 (matrix moments at k = 3..6 and the s -> inf expansion of AQ's recursion at k = 3 agree; k = 0,1,2 have
      3^j; 1,2,5,17,74,365,1889,9938; 1,2,6,22,91,408,1928,9419).
  (4) [upper certificate, every k >= 0] with A_hi = 0.0167788801732 (exact rational):
      (i) coarse (all k, chord bounds of (2) inside BN.2's identity d' - a d/(5+y) = B' - aB/(5+y) + r(K_0 + aC/(5+y))
          + C p'((-r') - r/p)): d_k' - A_hi d_k/(5+y) < 0 for y outside (12, 15) (exact polynomial, two positive roots
          12.5408..., 14.7913...);
      (ii) fine (k >= 3, J = 7/8 moment bounds, C' > 0, C p' > 0): (5+y) d'_hi - A_hi d_lo < 0 on [12, 15] (exact
          polynomial of degree 17, no root on [12,15], Sturm);
      (iii) exact (k = 0, 1, 2, AQ's recursion, rational d_k): A d_k - (5+y) d_k' has no positive root at A = a_k,hi.
      So omega_k is strictly increasing for every k >= 0 and every alpha >= A_hi (extends BN.2 from 1/17).
  (5) [lower certificate = Kept Failure BO.1] at y* = 13.778505644 the J = 7/8 bounds give, for every k >= 3, exactly
      Phi_k(y*) >= 0.016778880172852..., so (log omega_k)'(y*) < 0 for every alpha < A_lo = 0.016778880172852 and every k >= 3;
      also a_0* in [0.0155298932937, 0.0155298932938], a_1* in [0.0167694660209, 0.0167694660210],
      a_2* in [0.0167788799215, 0.0167788799216] (exact rationals: sign change of the cleared numerator).
      Hence a* = sup_k a_k* = a_k* (k >= 3) in [A_lo, A_hi], width < 3.5e-13.
Floats (consistency only): mpmath sup Phi_k at k = 0..10 (k >= 3 agree to 16 digits, 0.0167788801728761 at y = 13.7785056).
Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_maxwell_limit as ml  # noqa: E402

FAILS: list[str] = []
y = sp.symbols("y", positive=True)
P = y * (y + 5)
PP = sp.diff(P, y)
B = 1 - y * (y + 4) / (3 * (y + 2) * (y + 5))
E = 2 * y / (2 + y)
C = sp.Rational(5, 3) * (1 + 2 * y / 5) * E
BP, CP = sp.diff(B, y), sp.diff(C, y)
K0 = -CP + C * PP / P
M_STABLE = [1, 2, 6, 22, 92, 418, 1994, 9780]
M_LOW = {0: [3**j for j in range(8)], 1: [1, 2, 5, 17, 74, 365, 1889, 9938], 2: [1, 2, 6, 22, 91, 408, 1928, 9419]}
A_HI = sp.Rational(167788801732, 10**13)
A_LO = sp.Rational(16778880172852, 10**15)
Y_STAR = sp.Rational(137785056440, 10**10)
WIN = (sp.Integer(12), sp.Integer(15))
LEVEL_BRACKETS = {  # k: (lo, hi, rational witness point)
    0: (sp.Rational(155298932937, 10**13), sp.Rational(155298932938, 10**13), sp.Rational(154553084579, 10**10)),
    1: (sp.Rational(167694660209, 10**13), sp.Rational(167694660210, 10**13), sp.Rational(137981744772, 10**10)),
    2: (sp.Rational(167788799215, 10**13), sp.Rational(167788799216, 10**13), sp.Rational(137785065875, 10**10)),
}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def r_rec(k, s):
    """AQ's recursion (r_0 = 1/(3+s)); matches (1/2) v^T (L_k+s)^-1 v (BN's exact check, n = 0..3)."""
    if k == 0:
        return 1 / (3 + s)
    x = r_rec(k - 1, s * (s + 5))
    return (2 + s) * (5 + s) * x / ((6 + s) - 2 * s * (5 + s) * x)


def taylor_r(S, M, J):
    return sum((-1) ** j * sp.Integer(M[j]) / S ** (j + 1) for j in range(J))


def taylor_mrp(S, M, J):
    return sum((j + 1) * (-1) ** j * sp.Integer(M[j]) / S ** (j + 2) for j in range(J))


def cleared(expr):
    n, d = sp.fraction(sp.cancel(sp.together(expr)))
    return sp.Poly(sp.expand(n), y), sp.Poly(sp.expand(d), y)


def no_root_sign(expr, lo, hi, probe):
    """Exact: numerator and denominator of expr have no root on [lo,hi] (hi=None: [lo,inf)); returns sign of expr there."""
    n, d = cleared(expr)
    cn = n.count_roots(lo, hi) if hi is not None else len([r for r in n.real_roots() if r >= lo])
    cd = d.count_roots(lo, hi) if hi is not None else len([r for r in d.real_roots() if r >= lo])
    return cn == 0 and cd == 0, sp.sign(n.eval(probe)) * sp.sign(d.eval(probe)), n.degree()


# ----------------------------------------------------------------- (A) symbolic lemmas
def symbolic_checks():
    s, lam, a = sp.symbols("s lambda a", positive=True)
    ok = True
    for Jv in range(1, 9):
        part = sum((-lam) ** j / s ** (j + 1) for j in range(Jv))
        rem = (-lam) ** Jv / (s**Jv * (s + lam))
        ok = ok and sp.simplify(part + rem - 1 / (s + lam)) == 0
        part2 = sum((j + 1) * (-lam) ** j / s ** (j + 2) for j in range(Jv))
        rem2 = (-lam) ** Jv * (Jv / (s ** (Jv + 1) * (s + lam)) + 1 / (s**Jv * (s + lam) ** 2))
        ok = ok and sp.simplify(part2 + rem2 - 1 / (s + lam) ** 2) == 0
    check(ok, "1/(s+l) = sum_(j<J) (-l)^j s^-(j+1) + (-l)^J/(s^J(s+l)) and 1/(s+l)^2 = sum_(j<J) (j+1)(-l)^j s^-(j+2) + (-l)^J(J/(s^(J+1)(s+l)) "
              "+ 1/(s^J(s+l)^2)) for J = 1..8 (sympy): remainders have sign (-1)^J for l >= 0, so odd J gives upper, even J lower bounds")
    Lm = sp.Integer(8)
    chord = 1 / s - lam / (s * (s + Lm)) - 1 / (s + lam)
    chord2 = lam / (s + lam) ** 2 - lam / (s + Lm) ** 2
    check(sp.factor(sp.together(chord)) == sp.factor(lam * (Lm - lam) / (s * (s + lam) * (s + Lm)))
          and sp.factor(sp.together(chord2)) == sp.factor(lam * (Lm - lam) * (2 * s + lam + Lm) / ((s + lam) ** 2 * (s + Lm) ** 2)),
          "chord lemmas on [0,8]: 1/(s+l) <= 1/s - l/(s(s+8)) and l/(s+l)^2 >= l/(s+8)^2; with 1/(s+l)^2 = 1/(s(s+l)) - l/(s(s+l)^2) "
          "and m_1 >= 2: r_k <= 1/s - 2/(s(s+8)), -r_k' <= r_k/s - 2/(s(s+8)^2) for every k")
    # BN.2's identity re-used
    r = sp.Function("r")
    d = B - C * r(P)
    R0, R1 = sp.symbols("R0 R1")
    lhs = (sp.diff(d, y) - a * d / (5 + y)).replace(lambda ex: isinstance(ex, sp.Subs), lambda ex: R1).subs(r(P), R0)
    rhs = BP - a * B / (5 + y) + R0 * (K0 + a * C / (5 + y)) + C * PP * ((-R1) - R0 / P)
    check(sp.simplify(lhs - rhs) == 0 and sp.simplify(K0 - 2 * y * (2 * y**2 + 10 * y + 15) / (3 * (y + 2) ** 2 * (y + 5))) == 0,
          "BN.2 identity d' - a d/(5+y) = B' - aB/(5+y) + r(K_0 + aC/(5+y)) + C p'((-r') - r/p), K_0 >= 0 (sympy)")
    nC, dC = cleared(CP)
    check(all(c > 0 for c in nC.all_coeffs()) and all(c > 0 for c in dC.all_coeffs()),
          f"C' = {sp.factor(CP)} > 0 and C p' > 0 on y > 0: d' is increasing in -r' and decreasing in r")
    # reduction
    w = (1 + y / 5) ** a
    Dk = sp.Function("D")(y)
    lw = sp.diff(sp.log(w / Dk), y)
    check(sp.simplify(lw - (a - (5 + y) * sp.diff(Dk, y) / Dk) / (5 + y)) == 0,
          "(log omega_k)' = (a - Phi_k)/(5+y), Phi_k = (5+y) d_k'/d_k (sympy): omega_k increasing iff a >= sup Phi_k")


# ----------------------------------------------------------------- (B) moments
def moment_checks():
    ok = True
    for n in range(1, 7):
        _, adj, Cn, _, _ = ml.corner_frame(n)
        L = ml.laplacian(adj).astype(object)
        v = np.zeros(L.shape[0], dtype=object)
        v[Cn[0]], v[Cn[1]] = 1, -1
        w, ms = v.copy(), [1]
        for _ in range(7):
            w = L.dot(w)
            ms.append(int(v.dot(w)) // 2)
        ok = ok and int(max(L[i, i] for i in range(L.shape[0]))) == 4 and ms[1] == 2
        ok = ok and (ms == M_STABLE if n >= 3 else ms == M_LOW[n])
    s, t = sp.symbols("s t", positive=True)
    for k in range(0, 4):
        ser = sp.expand(sp.series(sp.cancel(r_rec(k, s)).subs(s, 1 / t), t, 0, 9).removeO())
        mk = [ser.coeff(t, j + 1) * (-1) ** j for j in range(8)]
        ok = ok and mk == (M_STABLE if k >= 3 else M_LOW[k])
    check(ok, "moments m_j = (1/2) v^T L_k^j v, j <= 7: graph matrices at k = 1..6 (max degree 4, m_1 = 2) and the s -> inf expansion of "
              "AQ's recursion at k = 0..3 agree; m = 1,2,6,22,92,418,1994,9780 at k = 3..6 (locality: constant for all k >= 3 since 7 < 2^3)")


# ----------------------------------------------------------------- (C) upper certificate
def upper_checks():
    rhat = 1 / P - 2 / (P * (P + 8))
    U = BP + rhat * K0 - 2 * C * PP / (P * (P + 8) ** 2)
    V = (B - C * rhat) / (5 + y)
    n, d = cleared(U - A_HI * V)
    roots = [rt for rt in n.real_roots() if rt > 0]
    okw = len(roots) == 2 and WIN[0] < roots[0] and roots[1] < WIN[1] and n.eval(0) != 0 and len([r for r in d.real_roots() if r >= 0]) == 0
    check(okw and sp.sign(n.eval(1)) * sp.sign(d.eval(1)) < 0,
          f"(i) coarse, every k: U - A_hi V < 0 on (0,inf) minus ({float(roots[0]):.4f}, {float(roots[1]):.4f}) inside (12,15); "
          f"B - C r_hat >= B - C/p > 0 (BN), so d_k' - A_hi d_k/(5+y) < 0 there")
    rlo, rhi = taylor_r(P, M_STABLE, 8), taylor_r(P, M_STABLE, 7)
    mhi = taylor_mrp(P, M_STABLE, 7)
    expr = (5 + y) * (BP - CP * rlo + C * PP * mhi) - A_HI * (B - C * rhi)
    ok2, sg, deg = no_root_sign(expr, *WIN, sp.Integer(13))
    check(ok2 and sg < 0, f"(ii) fine, every k >= 3: (5+y) d'_hi - A_hi d_lo < 0 on [12,15] (J = 7/8 moment bounds, cleared numerator "
                          f"degree {deg}, no root, Sturm): Phi_k < A_hi for all y > 0")
    ok3 = True
    for k, (lo, hi, yw) in LEVEL_BRACKETS.items():
        dk = sp.cancel(B - C * r_rec(k, P))
        ddk = sp.cancel(sp.diff(dk, y))
        okh, sgh, _ = no_root_sign(hi * dk - (5 + y) * ddk, sp.Integer(0), None, sp.Integer(1))
        wl = sp.cancel(lo * dk - (5 + y) * ddk).subs(y, yw)
        ok3 = ok3 and okh and sgh > 0 and wl < 0 and hi < A_HI
        print(f"   k={k}: a_k* in [{float(lo):.13f}, {float(hi):.13f}] (witness y = {float(yw):.10f})")
    check(ok3, "(iii) exact, k = 0,1,2: a_k,hi d_k - (5+y) d_k' > 0 on [0,inf) (no root) and a_k,lo d_k - (5+y) d_k' < 0 at the witness; "
               "all a_k,hi < A_hi")


# ----------------------------------------------------------------- (D) lower certificate (Kept Failure BO.1)
def lower_checks():
    S = P.subs(y, Y_STAR)
    rlo, rhi = taylor_r(S, M_STABLE, 8), taylor_r(S, M_STABLE, 7)
    mlo = taylor_mrp(S, M_STABLE, 8)
    sub = lambda ex: ex.subs(y, Y_STAR)
    dp_lo = sub(BP) - sub(CP) * rhi + sub(C * PP) * mlo
    d_hi = sub(B) - sub(C) * rlo
    phi_lo = (5 + Y_STAR) * dp_lo / d_hi
    print(f"   Phi_k(y*) >= {sp.N(phi_lo, 18)} for every k >= 3 (exact rational, y* = {Y_STAR})")
    check(dp_lo > 0 and d_hi > 0 and phi_lo > A_LO and A_LO < A_HI and A_HI - A_LO < sp.Rational(35, 10**14),
          "KEPT FAILURE BO.1: for every k >= 3 and alpha < A_lo = 0.016778880172852, (log omega_k)'(y*) < 0: omega_k is not monotone; "
          "a* in [A_lo, A_hi], width < 3.5e-13")


# ----------------------------------------------------------------- (E) floats
def float_sup():
    import mpmath as mp
    mp.mp.dps = 40
    out = {}
    for k in (0, 1, 2, 3, 5, 10):
        def phi(t, k=k):
            g = lambda u: (1 - u * (u + 4) / (3 * (u + 2) * (u + 5))
                           - mp.mpf(5) / 3 * (1 + 2 * u / 5) * 2 * u / (2 + u) * ml.rec_all(k, u * (u + 5))[0])
            return (5 + t) * mp.diff(g, t) / g(t)
        yb = mp.findroot(lambda t: mp.diff(phi, t), mp.mpf(14))
        out[k] = (phi(yb), yb)
    return out


def main() -> int:
    t0 = time.time()
    print("Theorem BO: the sharp exponent threshold for BN's monotone weight omega_k. Not thrust.")
    symbolic_checks()
    moment_checks()
    upper_checks()
    lower_checks()
    fl = float_sup()
    for k, (v, yb) in fl.items():
        print(f"   floats: sup Phi_{k} = {float(v):.16f} at y = {float(yb):.10f}")
    check(all(LEVEL_BRACKETS[k][0] <= fl[k][0] <= LEVEL_BRACKETS[k][1] for k in (0, 1, 2))
          and all(A_LO <= fl[k][0] <= A_HI for k in (3, 5, 10)),
          "floats (consistency only): mpmath sup Phi_k inside the certified brackets at k = 0,1,2,3,5,10")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BO checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
