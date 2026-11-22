#!/usr/bin/env python3
"""Theorem ES: the sum and the low moments of H_n's residues at every level n, and their magnitude.  ER kept "signs and counts only:
no magnitude, no sum" for H_n's residues.

Setting (EQ.1, EQ.2, ER.1 as proved): p(s) = s(s+5), q_k = p^k + 5, r_n = B_n/A_n, A_n monic of degree 2^n with 2^n simple real roots
z_j < 0 (w_j = -z_j), c_j = B_n(z_j)/A_n'(z_j) > 0, sum c_j = 1; H_n's integrand h_n = Q_n/(2 A_n), Q_n = prod_(k<=n-2) q_k (monic, degree
2^(n-1) - 1, Q_n(0) = 5^(n-1)), reduced; residues h_j = c_j/(2 q_(n-1)(z_j)); H_n = int_0^oo s^-a h_n = pi/sin(pi a) sum_j h_j w_j^-a (EP.1/EQ.2).

  ES.1 (THEOREM, every n >= 1).  deg(2 A_n) - deg Q_n = 2^(n-1) + 1 =: d_n >= 2, so (residue theorem at infinity / Laurent expansion
         h_n = sum_(k>=1) (sum_j h_j z_j^(k-1)) s^-k):
           sum_j h_j z_j^m = 0   for m = 0, 1, ..., 2^(n-1) - 1,        sum_j h_j z_j^(2^(n-1)) = 1/2,
           sum_j h_j / w_j = h_n(0) = 5^(n-1) / (2 * 3^(n+1)).
         In particular sum_j h_j = 0: the 2^(n-2) negative residues (ER.1; one for n = 1) cancel the positive ones exactly,
         S_n^- := sum_(h_j<0) |h_j| = S_n^+ := sum_(h_j>0) h_j = (1/2) sum_j |h_j|.
  ES.2 (THEOREM, every n >= 1; corollary of ES.1).  F_n(a) := sum_j h_j w_j^-a is entire in a, F_n(-m) = 0 for m = 0..2^(n-1)-1,
         F_n(1) = h_n(0); on (0,1) H_n = pi F_n(a)/sin(pi a); H_n extends continuously to a = 0 with
           H_n(0) = int_0^oo h_n(s) ds = -sum_j h_j ln w_j  (= F_n'(0); dominated convergence and l'Hopital, F_n(0) = 0).
  ES.3 (THEOREM, every n >= 2).  On each interval of {p^(n-1) < -5}, p^(n-1) = p(p^(n-2)) >= min p = -25/4, so -5/4 <= q_(n-1) < 0 there and
         |h_j| >= (2/5) c_j for each negative-residue pole: S_n^- >= (2/5) C_n^-, C_n^- := sum_(h_j<0) c_j (the r_n-mass on those poles).
         n = 1 exactly: S_1^- = S_1^+ = 1/(2 sqrt13) (r_1 = (s+5)/(s^2+7s+9), h_1 = 1/(2(s^2+7s+9))).
  ES.4 (checks / attempts to break).  (a) EXACT (pure-Python integers), n = 1..NEX (= 10): the Laurent coefficients of h_n at infinity
         (= the residue power sums) are 0 for m < 2^(n-1) and 1/2 at m = 2^(n-1); h_n(0) = 5^(n-1)/(2 3^(n+1)).  (b) CERTIFIED (arb, 1024
         bits for n <= 7, 2048 bits for n = 8), n = 1..NCERT (= 8): residues from isolated roots; every power sum for m <= 2^(n-1) contains 0 (resp. 1/2) with radius
         < 1e-100; sum h_j/w_j contains h_n(0); S_n^-, S_n^+ enclosed and overlapping; S_n^- >= (2/5) C_n^-; S_1^- contains 1/(2 sqrt13).
         A first run at 1024 bits for every level failed ES.4(b) at n = 8 on the radius threshold only (all containments, signs and the
         floor held; residue radii ~5.7e-96, m = 128 power-sum radius ~1.5e-6): precision, not a counterexample; hence 2048 bits at n = 8.
         (c) CHECKED (mpmath quadrature, n = 1..4): int_0^oo h_n = -sum h_j ln w_j, and H_n(a) at a = 0.3, 0.85 = residue sum.
         OBSERVED (not claimed): S_n^- = 0.1387, 0.1721, 0.0957, 0.0514, 0.0255, 0.0120, 0.0056, 0.0025, decreasing for n >= 2 with ratio
         S_(n+1)^-/S_n^- near 0.45-0.56 (n <= 8); no limit, rate or decay is proved.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_pole_structure as eq  # noqa: E402
import gasket_fractional_corner_h_residue_signs as er  # noqa: E402

FAILS: list[str] = []
NEX = 10
NCERT = 8
NQUAD = 4
PREC_HI_FROM = 8      # levels n >= 8 are certified at 2048 bits: at 1024 bits (n = 8) every containment, sign and floor check holds
PREC = {}             # but the residue radii are ~6e-96 and the m = 128 power sum's ~1.5e-6, above the 1e-100 threshold (precision only)


def prec_for(n: int) -> int:
    return PREC.get(n, 2048 if n >= PREC_HI_FROM else 1024)
RESULTS: dict = {}


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def q_pair(n: int):
    """(Q_n, A_n) as integer coefficient lists (low first): Q_n = prod_(k<=n-2) (p^k + 5)."""
    A, _ = eq.eq1_pair(n)
    return eq.b_product(n - 1), A


def power_sums_exact(n: int, K: int):
    """sum_j h_j z_j^m for m = 0..K-1, exactly, from the Laurent expansion of h_n = Q_n/(2 A_n) at infinity.
    With u = 1/s, a(u) = u^d A(1/u), b(u) = u^d Q(1/u): Q/A = b/a = sum_(k>=1) P_(k-1) u^k (a(0) = 1, so integer coefficients)."""
    Q, A = q_pair(n)
    d = len(A) - 1
    a = A[::-1]                                  # a[i] = coefficient of u^i
    b = [0] * (d + 1)
    for i, c in enumerate(Q):
        b[d - i] = c
    ser = []
    for k in range(K + 1):
        v = (b[k] if k <= d else 0) - sum(a[i] * ser[k - i] for i in range(1, min(k, d) + 1))
        ser.append(v)                            # a[0] = 1
    return [Fraction(ser[m + 1], 2) for m in range(K)]


def exact_part(n: int) -> None:
    Q, A = q_pair(n)
    d = (len(A) - 1) - (len(Q) - 1)
    half = 2 ** (n - 1)
    ps = power_sums_exact(n, half + 1)
    h0 = Fraction(eq.peval(Q, 0), 2 * eq.peval(A, 0))
    ok = d == half + 1 and Q[-1] == 1 and A[-1] == 1 and all(x == 0 for x in ps[:half]) and ps[half] == Fraction(1, 2) \
        and h0 == Fraction(5 ** (n - 1), 2 * 3 ** (n + 1))
    check(ok, f"ES.4(a) n = {n} (exact): deg gap {d} = 2^{n - 1}+1; sum h_j z_j^m = 0 for m < {half}, = 1/2 at m = {half}; "
              f"h_n(0) = 5^{n - 1}/(2*3^{n + 1}) = {h0}")


def certify(n: int) -> None:
    from flint import arb, ctx, fmpz_poly
    t1 = time.time()
    old = ctx.prec
    bits = prec_for(n)
    ctx.prec = bits
    try:
        Ql, Al = q_pair(n)
        A, Q, q = fmpz_poly(Al), fmpz_poly(Ql), fmpz_poly(er.q_poly(n))
        _, Bl = eq.eq1_pair(n)
        B = fmpz_poly(Bl)
        Ad = A.derivative()
        poles = sorted((r.real for r, _ in A.complex_roots()), key=lambda x: float(x.mid()))
        c = [B(z) / Ad(z) for z in poles]
        h = [Q(z) / (2 * Ad(z)) for z in poles]
        h2 = [cj / (2 * q(z)) for cj, z in zip(c, poles)]
        same = all(bool(abs(x - y) < arb(10) ** -100) for x, y in zip(h, h2))
        half = 2 ** (n - 1)
        tiny = arb(10) ** -100
        sums, ok_m = [], True
        pw = [arb(1)] * len(poles)
        for m in range(half + 1):
            s = sum((hj * pj for hj, pj in zip(h, pw)), arb(0))
            target = arb(0) if m < half else arb(1) / 2
            ok_m = ok_m and bool(abs(s - target) < tiny) and s.contains(target)
            sums.append(s)
            pw = [pj * z for pj, z in zip(pw, poles)]
        hw = sum((hj / (-z) for hj, z in zip(h, poles)), arb(0))
        ok_h0 = hw.contains(arb(5 ** (n - 1)) / (2 * 3 ** (n + 1))) and bool(hw.rad() < tiny)
        decided = all(bool(x > 0) or bool(x < 0) for x in h)
        neg = [i for i, x in enumerate(h) if bool(x < 0)]
        Sm = -sum((h[i] for i in neg), arb(0))
        Sp = sum((x for x in h if bool(x > 0)), arb(0))
        Cm = sum((c[i] for i in neg), arb(0))
        ok_s = decided and len(neg) == er.predicted_neg(n) and bool(abs(Sm - Sp) < tiny)
        ok_lb = True if n == 1 else bool(Sm >= arb(2) / 5 * Cm)
        ok_n1 = True if n != 1 else bool(abs(Sm - 1 / (2 * arb(13).sqrt())) < tiny)
        # ES.2: H_n(0) = -sum h_j ln w_j
        H0 = -sum((hj * (-z).log() for hj, z in zip(h, poles)), arb(0))
        RESULTS[n] = {"S": Sm, "C": Cm, "H0": H0, "h0": hw}
        print(f"     n = {n}: S_n^- = {Sm.str(10, radius=False)}  C_n^- = {Cm.str(10, radius=False)}  S/C = {(Sm / Cm).str(6, radius=False)}  "
              f"H_n(0) = -sum h ln w = {H0.str(10, radius=False)}  (rad S {float(Sm.rad()):.1e})", flush=True)
    finally:
        ctx.prec = old
    check(same and ok_m and ok_h0 and ok_s and ok_lb and ok_n1,
          f"ES.4(b) n = {n} (arb, {bits} bits): Q_n/(2A_n') = c_j/(2q_(n-1)) at every pole; power sums m < {2 ** (n - 1)} contain 0, m = {2 ** (n - 1)} "
          f"contains 1/2 (radius < 1e-100); sum h_j/w_j contains h_n(0); S_n^- = S_n^+; {er.predicted_neg(n)} negative"
          + ("" if n > 1 else "; S_1^- = 1/(2 sqrt13)") + ("; S_n^- >= (2/5) C_n^-" if n > 1 else "") + f" ({time.time() - t1:.1f}s)")


def quad_check(n: int) -> None:
    import mpmath as mp
    mp.mp.dps = 30
    Q, A = q_pair(n)
    f = lambda s: mp.polyval(Q[::-1], s) / (2 * mp.polyval(A[::-1], s))  # noqa: E731
    roots = sorted(mp.polyroots(A[::-1], maxsteps=200, extraprec=200), key=lambda z: mp.re(z))
    Ad = [i * A[i] for i in range(1, len(A))]
    h = [mp.re(mp.polyval(Q[::-1], z) / (2 * mp.polyval(Ad[::-1], z))) for z in roots]
    w = [-mp.re(z) for z in roots]
    I0 = mp.quad(f, [0, 1, 10, mp.inf])
    R0 = -sum(hj * mp.log(wj) for hj, wj in zip(h, w))
    ok = abs(I0 - R0) < mp.mpf(10) ** -20
    for a in (mp.mpf("0.3"), mp.mpf("0.85")):
        Ia = mp.quad(lambda t: f(t ** (1 / (1 - a))) / (1 - a), [0, 1, 10, mp.inf])   # s = t^(1/(1-a)) removes s^-a
        Ra = mp.pi / mp.sin(mp.pi * a) * sum(hj * wj ** (-a) for hj, wj in zip(h, w))
        ok = ok and abs(Ia - Ra) < mp.mpf(10) ** -15
    if n in RESULTS:
        ok = ok and abs(R0 - mp.mpf(RESULTS[n]["H0"].mid().str(30, radius=False))) < mp.mpf(10) ** -20
    check(ok, f"ES.4(c) n = {n} (mpmath quadrature): int_0^oo h_n = {mp.nstr(I0, 12)} = -sum h_j ln w_j; H_n(0.3), H_n(0.85) = residue sums")


def main(ncert: int = NCERT) -> int:
    t0 = time.time()
    print("Theorem ES: sum, low moments and magnitude of H_n's residues at every level n. Not thrust.", flush=True)
    for n in range(1, NEX + 1):
        exact_part(n)
    for n in range(1, ncert + 1):
        certify(n)
    for n in range(1, NQUAD + 1):
        quad_check(n)
    if len(RESULTS) >= 3:
        S = [float(RESULTS[n]["S"].mid()) for n in sorted(RESULTS)]
        print("   OBSERVED (not claimed): S_n^- =", ", ".join(f"{x:.6g}" for x in S), "; ratios",
              ", ".join(f"{S[i + 1] / S[i]:.4f}" for i in range(len(S) - 1)), flush=True)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all ES checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else NCERT))
