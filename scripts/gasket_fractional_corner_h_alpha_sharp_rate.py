#!/usr/bin/env python3
"""Theorem FB: the sharp rate 3^-n for EX.1's limit rho^-n H_n(alpha) -> W(alpha)/2 on -1 <= Re alpha < 1, its exact leading coefficient,
and approach from below for real alpha (EY.4's observed 3^-n and EX.4's observed 'increasing towards W/2', for n >= 11).

Setting (EY/FA): G_n(alpha) = rho^-n H_n(alpha) = int t^-alpha f_n, W(alpha)/2 = int t^-alpha f,  f_n = R_n/(2(5+psi_(n-1)(t/5))),
f = R/(2(5+psi(t/5))), g(y) = (5/3)(y+3)/(y+5), P_n = prod_{j<=n} g(psi_(n-j)(5^-j t)), P = prod_{j>=1} g(psi(5^-j t)) = 1/(3R),
1/R_n = 3 P_n + (t/2) 3^-n (1 + (3/5)^n P_n)  (EU.1 closed form, as in EY).  c* = 4 - sqrt15 (FA.1), gam = log_5 2.

  FB.1 (THEOREM; exact split + pointwise bounds).  f_n/f = e^Delta/(1+E) exactly, with
        E = (t/6) 3^-n (1/P_n + (3/5)^n),   (t/2) 3^-n R(t) <= E <= (t/6) 3^-n (1 + (3/5)^n)     (1 <= P_n <= P since g >= 1, psi_m <= psi),
        Delta = log(P/P_n) + log((5+psi(t/5))/(5+psi_(n-1)(t/5))),
        0 <= Delta <= Lam_n(t) := w [ (c* n + 1) kap M / 10 + e^(w/50)/30 ],  w = 5^-n t <= 10,  M = max(1,t)^gam,  kap = exp(e^(1/10)/4).
     Proof of the Delta bound: log psi(u) - log psi_m(u) = log p^m(psi(v)) - log p^m(v), v = 5^-m u; mean value in log y with
     d log p^m / d log y = prod_k (1 + y_k/(y_k+5)), increasing in y, at most Pi(u) := prod_{i>=1} (1 + psi(5^-i u)/(psi(5^-i u)+5)) (Koenigs
     psi(5v) = p(psi(v))); log psi(v) - log v <= v/10 (FA.2, v <= 10).  Pi(u) <= kap for u <= 5 (log(1+x/(x+5)) <= x/5, psi(s) <= s e^(s/10)),
     Pi(u) <= 2^k kap <= kap u^gam on (5^k, 5^(k+1)].  Hence log(psi(u)/psi_m(u)) <= (5^-m u/10) kap max(1,u)^gam.  With u = 5^-j t, m = n-j:
     5^-m u = w; elasticity of g <= c* (FA.1); j > n part of log P <= e^(w/50) w/30 (FA.2); second factor <= log(psi/psi_(n-1))(t/5) with m = n-1.
     So -E <= f_n/f - 1 <= e^Lam - 1  and  1 - f_n/f >= E/(1+E) - (e^Lam - 1).
  FB.2 (THEOREM; the 3^-n rate).  For n >= 9 and every complex alpha with -1 <= a = Re alpha < 1:
        |rho^-n H_n(alpha) - W(alpha)/2| <= C_n(a) 3^-n,
        C_n(a) = (1 + (3/5)^n) W(a-1)/12 + del_n (W(a-1) + W(a-3/2))/2 + tau_n,   del_n = (3/5)^n chi* beta_n,
        beta_n = (c* n + 1) kap/10 + e^(1/50)/30, chi* = (e^L* - 1)/L*, L* >= sup_{n>=9} (3^(1+gam)/5)^n beta_n,  tau_n = 3^n * EZ-type tail beyond T_n = 3^n.
     On t <= 3^n: |f_n/f - 1| <= E + (e^Lam - 1), e^Lam - 1 <= chi* Lam, M <= max(1,t)^(1/2), int t^(1-a) max(1,t)^(1/2) f <= (W(a-1)+W(a-3/2))/2.
     Beyond 3^n: f_n + f <= 1/(3(5+psi_(n-1)(t/5))), psi_(n-1)(t/5) >= (t/5^(J+1))^(2^J), t^-a <= t (EX.1/EY.3).
     Uniform: W decreasing on (-oo, alpha*] (EX.2) gives C_n(a) <= C_9(-1)-type constants printed below (certified).
  FB.3 (THEOREM; exact leading coefficient).  For -1 <= Re alpha < 1:  3^n (W(alpha)/2 - rho^-n H_n(alpha)) -> V(alpha) :=
        (1/2) int t^(1-alpha) f R dt = (1/4) int_0^oo t^(1-alpha) R(t)^2/(5+psi(t/5)) dt   (dominated convergence with FB.1's bounds).
     V(0) = EV's c = 2.67182281196926907...; V(alpha) > 0 for real alpha.
  FB.4 (THEOREM, arb; approach from below and two-sided Theta(3^-n)).  For every real alpha in [-1, 1) and every n >= 11:
        W(alpha)/2 - rho^-n H_n(alpha) >= m_i 3^-n > 0  on each partition interval [a_i, a_(i+1)] (m_i certified below).
  FB.5 (OBSERVED, not proved).  mpmath (40 digits) values of 3^n (W/2 - G_n) for n <= 40 at alpha in {-1,-1/2,0,1/2,9/10,99/100} approach V(alpha);
        n <= 10 approach from below is observed only.
Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0])

PART = [Fraction(-1), Fraction(-9, 10), Fraction(-4, 5), Fraction(-7, 10), Fraction(-3, 5), Fraction(-1, 2), Fraction(-3, 8),
        Fraction(-1, 4), Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1)]
N0_RATE = 9      # FB.2 threshold
N0_BELOW = 11    # FB.4 threshold
X_EXP = 5        # FB.4 uses t <= 5^5
EV_C = 2.6718228119692690733


# ---------------- pure-Python float pieces (CI) ----------------
def fconst():
    cs = 4 - math.sqrt(15)
    kap = math.exp(math.exp(0.1) / 4)
    gam = math.log(2) / math.log(5)
    q = 3 ** (1 + gam) / 5
    return cs, kap, gam, q


def beta_f(n):
    cs, kap, _, _ = fconst()
    return (cs * n + 1) * kap / 10 + math.exp(0.02) / 30


def lam_bar_f(n):
    return fconst()[3] ** n * beta_f(n)


def best_J(n, a=-1.0):
    best, bJ = None, None
    for J in range(2, n):
        s0 = 3.0 ** n / 5.0 ** (J + 1)
        if s0 <= 1:
            break
        v = n * math.log(3) + math.log(1 / 3) + (J + 1) * (1 - a) * math.log(5) + (1 - a - 2 ** J) * math.log(s0) - math.log(2 ** J + a - 1)
        if best is None or v < best:
            best, bJ = v, J
    return bJ, best  # log(tau_n) at a = -1


def tail_log_analytic_ok(n):
    """n >= 30, J = floor(n/3): log tau_n <= n(ln3 + 2 ln5/3) + 2 ln5 - (2^J - 2)(n(ln3 - ln5/3) - ln5) < 0 (decreasing in n)."""
    J = n // 3
    l3, l5 = math.log(3), math.log(5)
    ls0 = n * l3 - (n / 3 + 1) * l5
    return J >= 10 and ls0 > 0 and n * (l3 + 2 * l5 / 3) + 2 * l5 - (2 ** J - 2) * ls0 < -1000


def split_float(n, t):
    """FB.1 float sanity: returns (E, Delta, Lam) via the closed forms (EY helpers)."""
    import gasket_fractional_corner_h_alpha_rate as ey
    cs, kap, gam, _ = fconst()
    s = t / 5.0 ** n
    Pn, z = 1.0, s
    for k in range(n):
        Pn *= ey.g(z)
        if k < n - 1:
            z = ey.p(z)
    P = ey.inv_R(t) / 3
    E = t / 6 * 3.0 ** -n * (1 / Pn + 0.6 ** n)
    Delta = math.log(P / Pn) + math.log((5 + ey.psi(t / 5)) / (5 + z))
    w = t / 5.0 ** n
    M = max(1.0, t) ** gam
    Lam = w * ((cs * n + 1) * kap * M / 10 + math.exp(w / 50) / 30)
    ratio = ey.f_n(n, t) / ey.f_inf(t)
    return E, Delta, Lam, ratio, P


def V_float(a, lo=1e-12, hi=5.0 ** 8, N=4000):
    import gasket_fractional_corner_h_alpha_rate as ey
    x0, x1 = math.log(lo), math.log(hi)
    h = (x1 - x0) / N
    s = 0.0
    for i in range(N + 1):
        t = math.exp(x0 + i * h)
        R = 1 / ey.inv_R(t)
        s += t ** (2 - a) * R * R / (4 * (5 + ey.psi(t / 5))) * (1 if i in (0, N) else (4 if i % 2 else 2))
    return s * h / 3


# ---------------- arb certificates ----------------
def arb_setup():
    import gasket_fractional_corner_h_zero_rate_constant as ev
    ev.PREC = 128
    ev.TINY_EXP = 90
    acb, arb, RY, psi, _fc, _fk = ev.make_arb()
    return ev, acb, arb, RY, psi


def W_cert(acb, arb, RY, psi, a: Fraction, low=24, tol=40):
    """Enclosure of W(a) for rational a < 1 (EX.3's routine, looser tolerance)."""
    eps = arb(5) ** (-low)
    edges = [arb(5) ** m for m in range(-low, 9)]
    R8 = RY(5 ** 8)[0].real.upper()
    p7 = psi(5 ** 7).real.lower()
    w = eps / 5
    psi_up = w + w * w / 20 + w ** 3 / (400 * (1 - w / 20))
    L_eps = (1 / ((3 + eps) * (5 + psi_up))).lower()
    al = arb(a.numerator) / a.denominator

    def f(z, analytic):
        R, y1 = RY(z)
        return R / (5 + y1) * (-al * z.log()).exp()

    tot = arb(0)
    for lo_, hi_ in zip(edges, edges[1:]):
        tot += acb.integral(f, lo_, hi_, rel_tol=arb(2) ** -tol, abs_tol=arb(2) ** -(tol + 10)).real
    one = 1 - al
    I = eps ** one / one
    neg = max(Fraction(0), -a)
    negb = arb(neg.numerator) / neg.denominator
    ratio = arb(5) ** (1 + negb) / p7
    assert bool(ratio <= arb(1) / 2)
    T8 = R8 * 4 * arb(5) ** 8 * arb(5) ** (negb * 9) / p7
    return tot.lower() + (L_eps * I).lower(), tot.upper() + (I / 15).upper() + (2 * T8).upper()


def V_piece_lower(acb, arb, RY, b: Fraction, mlo: int, mhi: int, tol=30):
    """Lower bound of (1/4) int_{5^mlo}^{5^mhi} t^(1-b) R^2/(5+psi(t/5)) dt (integrand >= 0)."""
    bb = arb(b.numerator) / b.denominator

    def f(z, analytic):
        R, y1 = RY(z)
        return R * R / (4 * (5 + y1)) * ((1 - bb) * z.log()).exp()

    tot = arb(0)
    for m in range(mlo, mhi):
        tot += acb.integral(f, arb(5) ** m, arb(5) ** (m + 1), rel_tol=arb(2) ** -tol, abs_tol=arb(2) ** -(tol + 10)).real
    return tot.lower()


def main() -> int:
    t0 = time.time()
    print("Theorem FB: 3^-n rate, leading coefficient V(alpha), approach from below. Functions on a graph, not a flux, not thrust.", flush=True)
    fails = []

    def check(ok, msg):
        print(("  ok   " if ok else "  FAIL ") + msg, flush=True)
        if not ok:
            fails.append(msg)

    # ---- float sanity of FB.1 (not part of the proof) ----
    worst = 0.0
    okE = True
    for n in range(2, 13):
        for k in range(41):
            t = 1e-3 * (3.0 ** n * 1e3) ** (k / 40)
            E, D, L, ratio, P = split_float(n, t)
            okE &= (t / 2 * 3.0 ** -n / (3 * P) <= E * (1 + 1e-12)) and E <= t / 6 * 3.0 ** -n * (1 + 0.6 ** n) * (1 + 1e-12)
            okE &= abs(ratio - math.exp(D) / (1 + E)) <= 1e-10 * ratio
            okE &= D >= -1e-12
            if L > 0:
                worst = max(worst, D / L)
    check(okE and worst <= 1, f"FB.1 float sanity: f_n/f = e^Delta/(1+E), E bounds, 0 <= Delta <= Lam (max Delta/Lam = {worst:.3f}), n = 2..12, t <= 3^n")

    ev, acb, arb, RY, psi = arb_setup()
    from flint import arb as A
    cs = 4 - A(15).sqrt()
    kap = (A(1) / 10).exp().__truediv__(4).exp()
    gam = A(2).log() / A(5).log()
    q = A(3) ** (1 + gam) / 5
    check(bool(cs < A("0.12702")) and bool(kap < A("1.3183")) and bool(q < A("0.9637")), "constants: c* < 0.12702, kap = exp(e^0.1/4) < 1.3183, q = 3^(1+gam)/5 < 0.9637")
    # FB.1: kap bound pieces: e^(1/10) * 5 * sum_{j>=1} 5^-j / 5 = e^0.1/4 (exact geometric sum); 2^k = 5^(k gam) exact.
    A_ = cs * kap / 10
    B_ = kap / 10 + (A(1) / 50).exp() / 30

    def beta(n):
        return A_ * n + B_

    # sup_{n>=9} q^n beta_n: arb for n = 9..60, and for n >= 27 the ratio q (1 + A/(A n + B)) < 1 (decreasing)
    Lmax = max((q ** n * beta(n)).upper() for n in range(N0_RATE, 61))
    dec = all(bool(q * (1 + A_ / (A_ * n + B_)) < 1) for n in (27,))
    check(dec, "FB.2: q^n beta_n decreasing for n >= 27 (q (1 + A/(A n + B)) < 1 at n = 27, hence for all larger n)")
    Lstar = A("0.24")
    check(bool(Lmax < Lstar), f"FB.2: sup_(n>=9) Lam_bar_n = sup q^n beta_n <= {Lmax.str(6, radius=False)} < L* = 0.24")
    chi = (Lstar.exp() - 1) / Lstar
    print(f"   chi* = (e^0.24 - 1)/0.24 in {chi.str(10)}", flush=True)

    def delta(n):
        return A(3) ** n / A(5) ** n * chi * beta(n)

    check(all(bool(A(3) / 5 * (1 + A_ / (A_ * n + B_)) < 1) for n in (1,)), "FB.2: del_n = (3/5)^n chi* beta_n decreasing in n (ratio (3/5)(1 + A/(A+B)) < 1)")

    # tail tau_n, a = -1 worst (t^-a <= t on t >= 1)
    def tau(n):
        J, _ = best_J(n)
        s0 = A(3) ** n / A(5) ** (J + 1)
        assert bool(s0 > 1) and J <= n - 1 and J >= 2
        return (A(3) ** n * A(5) ** (2 * (J + 1)) * s0 ** (2 - 2 ** J) / (2 ** J - 2) / 3).upper()

    taus = {n: tau(n) for n in range(N0_RATE, 30)}
    tau9 = taus[N0_RATE]
    tau11 = max(taus[n] for n in range(N0_BELOW, 30))
    check(all(bool(taus[n] <= tau9) for n in taus) and all(tail_log_analytic_ok(n) for n in range(30, 400)),
          f"FB.2: tau_n <= tau_9 <= {tau9.str(5, radius=False)} for 9 <= n <= 29 (arb), tau_n < e^-1000 for n >= 30 (analytic, J = floor(n/3))")
    # exact analytic part for n >= 30: coefficient of n negative and value at n = 30 negative
    l3, l5 = A(3).log(), A(5).log()
    slope = l3 + 2 * l5 / 3 - 1022 * (l3 - l5 / 3)
    val30 = 30 * (l3 + 2 * l5 / 3) + 2 * l5 - 1022 * (30 * (l3 - l5 / 3) - l5)
    check(bool(slope < 0) and bool(val30 < -1000) and bool(30 * (l3 - l5 / 3) - l5 > 0),
          "FB.2: n >= 30 tail: log tau_n <= n(ln3 + 2ln5/3) + 2ln5 - 1022(n(ln3 - ln5/3) - ln5), slope < 0, value at 30 < -1000")
    print(f"   tau_11..29 <= {tau11.str(5, radius=False)}", flush=True)

    # ---- W values needed ----
    need = sorted({a - 1 for a in PART[:-1]} | {a - Fraction(3, 2) for a in PART[:-1]})
    Wv = {}
    for b in need:
        lo, hi = W_cert(acb, arb, RY, psi, b)
        Wv[b] = (lo, hi)
        print(f"   W({b}) in [{lo.str(12, radius=False)}, {hi.str(12, radius=False)}]  ({time.time() - t0:.0f}s)", flush=True)
    lo2, hi2 = Wv[Fraction(-2)]
    check(bool(lo2 < A("11291.83348548641137405205")) and bool(hi2 > A("11291.83348548641137405157")), "cross-check: W(-2) bracket meets FA.4")
    check(all(bool(Wv[x][0] > Wv[y][1]) for x, y in zip(need, need[1:])), "W strictly decreasing along the certified arguments (consistent with EX.2)")
    W2u, W52u = hi2, Wv[Fraction(-5, 2)][1]
    for n0 in (N0_RATE, N0_BELOW):
        C = (1 + A(3) ** n0 / A(5) ** n0) * W2u / 12 + delta(n0) * (W2u + W52u) / 2 + (tau9 if n0 == N0_RATE else tau11)
        print(f"   FB.2: for n >= {n0}, -1 <= Re alpha < 1:  |rho^-n H_n - W/2| <= C 3^-n,  C <= {C.upper().str(8, radius=False)}", flush=True)
    C9 = ((1 + A(3) ** 9 / A(5) ** 9) * W2u / 12 + delta(9) * (W2u + W52u) / 2 + tau9).upper()
    C11 = ((1 + A(3) ** 11 / A(5) ** 11) * W2u / 12 + delta(11) * (W2u + W52u) / 2 + tau11).upper()
    check(bool(C9 < 1331) and bool(C11 < 1096), "FB.2: uniform constants C(n>=9) < 1331, C(n>=11) < 1096")

    # ---- FB.4 ----
    Xb = A(5) ** X_EXP
    Eup = Xb / 6 / A(3) ** N0_BELOW * (1 + A(3) ** N0_BELOW / A(5) ** N0_BELOW)
    check(bool(Xb <= A(3) ** N0_BELOW), "FB.4: X = 5^5 <= 3^11 = T_11")
    margins = []
    for a0, a1 in zip(PART, PART[1:]):
        v = V_piece_lower(acb, arb, RY, a0, -24, 0) + V_piece_lower(acb, arb, RY, a1, 0, X_EXP)
        K = (Wv[a0 - 1][1] + Wv[a0 - Fraction(3, 2)][1]) / 2
        m = v / (1 + Eup) - delta(N0_BELOW) * K - tau11
        margins.append((a0, a1, m))
        print(f"   FB.4 [{a0}, {a1}]: V_lo >= {A(v).str(8, radius=False)}, K <= {K.upper().str(8, radius=False)}, m >= {m.lower().str(6, radius=False)}  ({time.time() - t0:.0f}s)", flush=True)
        check(bool(m > 0), f"FB.4: for n >= 11 and real alpha in [{a0}, {a1}]: W/2 - rho^-n H_n >= m 3^-n, m > 0")
    # cross-checks (float)
    v0 = V_float(0.0)
    check(abs(v0 - EV_C) < 1e-8, f"FB.3 cross-check: float V(0) = {v0:.12f} equals EV's c = {EV_C}")
    print(f"   float V(-1) = {V_float(-1.0):.9f}, V(-1/2) = {V_float(-0.5):.9f}, V(1/2) = {V_float(0.5):.9f}, V(9/10) = {V_float(0.9):.9f}", flush=True)
    print(f"elapsed {time.time() - t0:.1f}s", flush=True)
    if fails:
        print(f"{len(fails)} FAILED")
        return 1
    print("all FB checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
