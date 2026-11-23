#!/usr/bin/env python3
"""Theorem EX: the alpha-weighted refinement limit rho^-n H_n(alpha) -> W(alpha)/2 on the whole half-plane Re alpha < 1, and certified
values of W at real alpha (W is positive and strictly convex on (-oo, 1), with its unique minimum at some alpha* in (3/5, 7/10)).

Setting (as in ET-EW): p(s) = s(s+5), psi_k(t) = p^k(5^-k t) -> psi, R_n(t) = (3/5)^n r_n(5^-n t) -> R (BG.3/BD.4), h_n = r_n/(2 q_(n-1)),
H_n(alpha) = int_0^oo s^-alpha h_n ds, rho = 5^alpha/3, F(t) = t R(t)/(5 + psi(t/5)) >= 0, W(alpha) = int_0^oo F(t) t^(-alpha-1) dt.

  EX.1 (THEOREM).  For every complex alpha with Re alpha < 1 and every n with 2^(n-1) > -Re alpha, H_n(alpha) converges absolutely and
        rho^-n H_n(alpha) = int_0^oo t^-alpha R_n(t) / (2 (5 + psi_(n-1)(t/5))) dt   (ET.1's substitution, exact);
        rho^-n H_n(alpha) -> W(alpha)/2, locally uniformly on Re alpha < 1 (with all alpha-derivatives).
        Dominator: R_n <= 1/3 (EU.1/EU.3) and psi_k(t) >= psi_(k-J)(t/5^J)^(2^J) >= (t/5^J)^(2^J) (p(y) >= y^2, p(y) >= 5y on y >= 0), so for
        n - 1 >= J the integrand is <= t^(-Re alpha) / (6 (5 + (t/5^(J+1))^(2^J))), integrable when Re alpha < 1 < 2^J + Re alpha.
        Vitali (holomorphic, locally bounded, pointwise convergent) gives local uniformity.
  EX.2 (THEOREM).  On real alpha < 1: W(alpha) > 0, W'' > 0 (strictly convex), log W convex, (1 - alpha) W(alpha) -> 1/15 as alpha -> 1-,
        and H_(n+1)(alpha)/H_n(alpha) -> 5^alpha/3.
  EX.3 (THEOREM, certified arb).  Enclosures of W(alpha) at alpha = -1, -1/2, 0 (contains ln5), 1/2, 9/10 and of W'(alpha) at 3/5, 7/10:
        W'(3/5) < 0 < W'(7/10), so the unique minimiser alpha* of W lies in (3/5, 7/10).
        Pieces: arb integration over [5^m, 5^(m+1)], m = -LOW..7; low tail int_0^eps by L_eps <= F/t <= 1/15 (BI.4 floor, EV.1 psi bound);
        high tail beyond 5^8 by R <= R(5^8) (EV.3, R decreasing) and psi(5^(m-1)) >= psi(5^7)^(2^(m-8)) (EV.1), a ratio-1/2 series.
  EX.4 (OBSERVED, not proved).  Float quadrature of rho^-n H_n(alpha), alpha = -1/2, 1/2, n = 1..8: increasing towards W(alpha)/2, gaps
        shrinking by about 1/3 per level.  No rate is claimed for alpha != 0.

Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_h_zero_rate_constant as ev  # noqa: E402
import gasket_fractional_corner_h_zero_value as eu  # noqa: E402

FAILS = ev.FAILS
check = ev.check
LOW = 24
W_ALPHAS = [Fraction(-1), Fraction(-1, 2), Fraction(0), Fraction(1, 2), Fraction(9, 10)]
DW_ALPHAS = [Fraction(3, 5), Fraction(7, 10)]


def psi_k_exact(k: int, t: Fraction) -> Fraction:
    y = t / 5 ** k
    for _ in range(k):
        y = y * (y + 5)
    return y


def dominator_ok(kmax: int = 7) -> bool:
    """EX.1: psi_k(t) >= psi_(k-J)(t/5^J)^(2^J) >= (t/5^J)^(2^J), exact rationals."""
    pts = [Fraction(1, 7), Fraction(1), Fraction(13, 3), Fraction(40), Fraction(625), Fraction(10 ** 4)]
    for k in range(0, kmax + 1):
        for J in range(0, k + 1):
            for t in pts:
                a = psi_k_exact(k, t)
                b = psi_k_exact(k - J, t / 5 ** J) ** (2 ** J)
                c = (t / 5 ** J) ** (2 ** J)
                if not (a >= b >= c):
                    return False
    return True


def F_float(t: float) -> float:
    return t * ev.R_float(t) / (5.0 + ev.psi_float(t / 5.0))


def W_float(a: float, lo: float = 1e-30, hi: float = 5.0 ** 8, n: int = 8001, deriv: bool = False) -> float:
    """trapezoid in u = ln t of int F(e^u) e^(-a u) (-u)^deriv du, plus the leading low tail (CI)."""
    u0, u1 = math.log(lo), math.log(hi)
    h = (u1 - u0) / (n - 1)
    s = 0.0
    for i in range(n):
        u = u0 + i * h
        w = 0.5 if i in (0, n - 1) else 1.0
        s += w * F_float(math.exp(u)) * math.exp(-a * u) * ((-u) if deriv else 1.0)
    tail = lo ** (1 - a) / (15 * (1 - a))
    if deriv:
        tail = lo ** (1 - a) * (-math.log(lo) / (1 - a) + 1 / (1 - a) ** 2) / 15
    return s * h + tail


def scaled_H_float(n: int, a: float, N: int = 6001) -> float:
    """rho^-n H_n(a) = int t^-a R_n(t)/(2(5+psi_(n-1)(t/5))) dt (EX.1 identity), float log-grid trapezoid."""
    hi = 1e8 if n <= 3 else 1e6
    u0, u1 = math.log(1e-14), math.log(hi)
    h = (u1 - u0) / (N - 1)
    s = 0.0
    for i in range(N):
        u = u0 + i * h
        x = math.exp(u)
        sx = x / 5.0 ** n
        Rn = (3.0 / 5.0) ** n * float(eu.r_closed(n, sx))
        q = float(eu.pn(n - 1, sx)) + 5.0
        w = 0.5 if i in (0, N - 1) else 1.0
        s += w * Rn / (2 * q) * x * math.exp(-a * u)
    return s * h


def identity_exact(n: int, a: Fraction, t: Fraction) -> bool:
    """EX.1 integrand identity at a rational point for integer a: 5^-n s^-a h_n(s) * rho^-n = t^-a R_n(t)/(2(5+psi_(n-1)(t/5))), s = 5^-n t."""
    s = t / 5 ** n
    hn = eu.r_closed(n, s) / (2 * (eu.pn(n - 1, s) + 5))
    lhs = Fraction(3) ** n * Fraction(5) ** (-n * a) * Fraction(1, 5 ** n) * s ** (-a) * hn
    Rn = Fraction(3, 5) ** n * eu.r_closed(n, s)
    rhs = t ** (-a) * Rn / (2 * (5 + psi_k_exact(n - 1, t / 5)))
    return lhs == rhs


def certified(acb, arb, RY, psi):
    eps = arb(5) ** (-LOW)
    edges = [arb(5) ** m for m in range(-LOW, 9)]
    R8 = RY(5 ** 8)[0].real.upper()
    p7 = psi(5 ** 7).real.lower()
    w = eps / 5
    psi_up = w + w * w / 20 + w ** 3 / (400 * (1 - w / 20))
    L_eps = (1 / ((3 + eps) * (5 + psi_up))).lower()
    out = {}

    def run(a: Fraction, deriv: bool):
        al = arb(a.numerator) / a.denominator

        def f(z, analytic):
            R, y1 = RY(z)
            v = R / (5 + y1) * (-al * z.log()).exp()
            return v * (-z.log()) if deriv else v

        tot = arb(0)
        for lo_, hi_ in zip(edges, edges[1:]):
            tot += acb.integral(f, lo_, hi_, rel_tol=arb(2) ** -70, abs_tol=arb(2) ** -80).real
        one = 1 - al
        I = eps ** one / one if not deriv else eps ** one * (-eps.log() / one + 1 / one ** 2)
        low_lo, low_hi = (L_eps * I).lower(), (I / 15).upper()
        # high tail: T_m = R8 * 4 * 5^m * 5^(max(0,-a)(m+1)) * (deriv ? (m+1) ln5 : 1) / psi(5^7)^(2^(m-8)); ratio <= 1/2 => sum <= 2 T_8
        neg = max(Fraction(0), -a)
        negb = arb(neg.numerator) / neg.denominator
        ratio = arb(5) ** (1 + negb) * (arb(10) / 9 if deriv else 1) / p7
        check(bool(ratio <= arb(1) / 2), f"EX.3: high-tail ratio {ratio.upper().str(5)} <= 1/2 (alpha = {a}, deriv = {deriv})")
        T8 = R8 * 4 * arb(5) ** 8 * arb(5) ** (negb * 9) * (9 * arb(5).log() if deriv else 1) / p7
        hi_tail = (2 * T8).upper()
        # W: integrand >= 0 so the high tail is in [0, hi_tail]; W': -ln t < 0 there, so the high tail is in [-hi_tail, 0]
        if deriv:
            lo_b, hi_b = tot.lower() + low_lo - hi_tail, tot.upper() + low_hi
        else:
            lo_b, hi_b = tot.lower() + low_lo, tot.upper() + low_hi + hi_tail
        return lo_b, hi_b

    for a in W_ALPHAS:
        lo_b, hi_b = run(a, False)
        out[("W", a)] = (lo_b, hi_b)
        print(f"   EX.3: W({a}) in [{lo_b.str(22, radius=False)}, {hi_b.str(22, radius=False)}]  width {(hi_b - lo_b).str(3, radius=False)}", flush=True)
    for a in DW_ALPHAS:
        lo_b, hi_b = run(a, True)
        out[("dW", a)] = (lo_b, hi_b)
        print(f"   EX.3: W'({a}) in [{lo_b.str(22, radius=False)}, {hi_b.str(22, radius=False)}]  width {(hi_b - lo_b).str(3, radius=False)}", flush=True)
    return out


def main() -> int:
    t0 = time.time()
    print("Theorem EX: rho^-n H_n(alpha) -> W(alpha)/2 for Re alpha < 1; certified W at real alpha. Not thrust.", flush=True)
    # EX.1 exact pieces
    check(dominator_ok(), "EX.1: psi_k(t) >= psi_(k-J)(t/5^J)^(2^J) >= (t/5^J)^(2^J), exact, k <= 7, 6 points")
    ok = all(identity_exact(n, a, t) for n in range(1, 6) for a in (-2, -1, 0, 1) for t in (Fraction(1, 3), Fraction(7), Fraction(200)))
    check(ok, "EX.1: rho^-n 5^-n s^-alpha h_n(s) = t^-alpha R_n(t)/(2(5+psi_(n-1)(t/5))), s = 5^-n t, exact, n <= 5, alpha in {-2,-1,0,1}")
    okR = all(Fraction(3, 5) ** n * eu.r_closed(n, Fraction(t, 5 ** n)) <= Fraction(1, 3) for n in range(1, 7) for t in (0, 1, 9, 100, 5000))
    check(okR, "EX.1: R_n(t) <= 1/3 (EU.1/EU.3), exact spot check n <= 6")
    # EX.2: (1 - alpha) W(alpha) -> 1/15 (observed trend only as a check; the limit itself is proved from F(t)/t -> 1/15)
    print(f"   EX.2 (check): F(t)/t at t = 1e-8: {F_float(1e-8) / 1e-8:.12f} (1/15 = {1/15:.12f})", flush=True)
    check(abs(F_float(1e-8) / 1e-8 - 1 / 15) < 1e-8, "EX.2: F(t)/t -> 1/15 as t -> 0")
    # EX.3 certified
    acb, arb, RY, psi, _fc, _fk = ev.make_arb()
    out = certified(acb, arb, RY, psi)
    ln5 = arb(5).log()
    w0 = out[("W", Fraction(0))]
    check(bool(w0[0] <= ln5) and bool(ln5 <= w0[1]), "EX.3: W(0) bracket contains ln5 (EU, cross-check)")
    for a in W_ALPHAS:
        check(bool(out[("W", a)][0] > 0), f"EX.3: W({a}) > 0 certified")
    wh = out[("W", Fraction(1, 2))]
    check(bool(wh[1] < out[("W", Fraction(9, 10))][0]) and bool(wh[1] < w0[0]), "EX.3: W(1/2) < W(0), W(1/2) < W(9/10)")
    d1, d2 = out[("dW", Fraction(3, 5))], out[("dW", Fraction(7, 10))]
    check(bool(d1[1] < 0), "EX.3: W'(3/5) < 0 certified")
    check(bool(d2[0] > 0), "EX.3: W'(7/10) > 0 certified")
    for k, (lo_b, hi_b) in out.items():
        check(bool((hi_b - lo_b) < arb(10) ** -12), f"EX.3: bracket width for {k} < 1e-12")
    # float cross-check of the certified values
    for a in W_ALPHAS:
        lo_b, hi_b = out[("W", a)]
        mid = float(((lo_b + hi_b) / 2).mid())
        wf = W_float(float(a))
        check(abs(wf - mid) < 1e-6 * max(1, mid), f"EX.3: float W({a}) = {wf:.10f} agrees with arb midpoint {mid:.10f}")
    # EX.4 observed
    for a in (-0.5, 0.5):
        vals = [2 * scaled_H_float(n, a) for n in range(1, 9)]
        mid = float(((out[("W", Fraction(a).limit_denominator(10))][0] + out[("W", Fraction(a).limit_denominator(10))][1]) / 2).mid())
        gaps = [mid - v for v in vals]
        print(f"   EX.4 (OBSERVED): alpha={a}: 2 rho^-n H_n = " + ", ".join(f"{v:.8f}" for v in vals) + f"; W = {mid:.10f}", flush=True)
        print("        gaps " + ", ".join(f"{g:.3e}" for g in gaps) + "; ratios " + ", ".join(f"{gaps[i+1]/gaps[i]:.3f}" for i in range(len(gaps) - 1)), flush=True)
        check(all(g > 0 for g in gaps) and abs(gaps[-1]) < 0.01 * mid, f"EX.4: alpha={a}: finite-n values below W and within 1% by n = 8 (observation)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EX checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
