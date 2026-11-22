#!/usr/bin/env python3
"""Theorem ER: the signs of H_n's residues at every level n.  EQ kept "no sign for H_n's residues" (EQ.2: 'residues of mixed sign, not claimed').

Setting (EQ.1, EQ.2 as proved): p(s) = s(s+5); r_n = B_n/A_n, A_n has 2^n simple real negative roots z_j with residues c_j = B_n(z_j)/A_n'(z_j) > 0,
B_n = prod_(k<n) q_k, q_k = p^k + 5; H_n's integrand is Q_n/(2 A_n), Q_n = prod_(k<=n-2) q_k = B_n/q_(n-1), reduced, so its residue at z_j is
h_j = c_j / (2 q_(n-1)(z_j)), sign h_j = sign q_(n-1)(z_j) (q_(n-1)(z_j) != 0 since gcd(A_n, B_n) = 1).

  ER.1 (THEOREM, every n >= 1; uses only EQ.1, EQ.2).  Let N_n = number of poles z_j with q_(n-1)(z_j) < 0 (negative residues of H_n).  Then
         N_1 = 1 and N_n = 2^(n-2) for n >= 2;  the other 2^n - N_n residues are positive.  More precisely the negative residues sit
         one in each of the N_n bounded open intervals {z : p^(n-1)(z) < -5} (for n = 1 the ray (-inf, -5)).
       Proof.  (i) r_n = sum c_j/(s - z_j), c_j > 0, so r_n' < 0 on each gap and r_n runs from +inf to -inf: between consecutive poles there is exactly
       one zero of B_n, giving 2^n - 1 = deg B_n distinct real simple zeros; poles and zeros strictly alternate, with a pole at each end.
       (ii) q_(n-1) | B_n so its roots are among these zeros; they are all real and simple.  p(y) = -5 has the two roots (-5 +- sqrt5)/2 in
       (-3.62, -1.38) and p(y) = t has two real roots for every t > -25/4; all iterated preimages of -5 lie in (-5, 0) subset (-25/4, inf),
       so q_(n-1) has 2^(n-1) real roots, and as q_(n-1) > 0 near +-inf (n >= 2, even degree, monic) its negative set is 2^(n-2) disjoint open
       intervals (a_i, b_i), a_i < b_i adjacent roots (for n = 1, q_0 = s + 5 < 0 on (-inf, -5)).
       (iii) A zero w of q_k (k < n-1) satisfies p^(n-1)(w) = p^(n-1-k)(-5) = 0 > -5, so w is not in any (a_i, b_i); hence B_n has no zero inside
       (a_i, b_i), and a_i, b_i are adjacent zeros of B_n, so by (i) exactly one pole lies in (a_i, b_i).  For n = 1 the ray (-inf, -5) has the pole at left end
       (alternation: pole < zero -5 < pole) and contains exactly one.  Poles outside the negative set have q_(n-1) > 0 (not 0 by coprimality): positive
       residues.  QED.  No sum of H_n's residues, magnitude or location beyond this counting is claimed.
  ER.2 (checks / attempts to break; ER.1 is proved for every n, these only certify / check it at small n): (a) CERTIFIED, arb (1024 bits), n = 1..NCERT (= 8) only: isolate A_n's roots, evaluate q_(n-1) at the isolating balls: signs decided
       for every root, counts of negatives 1, 1, 2, 4, ... = ER.1; exactly one negative-residue pole in each interval (a_i, b_i) with a_i, b_i the roots of q_(n-1) from
       isolation; B_n's zeros alternate with the poles; (b) CHECKED, not certified: exact rationals at float Laplacian eigenvalues (n <= 4), independent of EQ's polynomials for the
       poles; (c) the ER.1 argument's inputs: p^m(-5) = 0 for m >= 1 exact, q_(n-1) has 2^(n-1) real roots for n <= 8.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_pole_structure as eq  # noqa: E402

FAILS: list[str] = []
NCERT = 8
NLAP = 4


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def q_poly(n: int):
    """q_(n-1) = p^(n-1) + 5 as an integer coefficient list."""
    f = [0, 1]
    for _ in range(n - 1):
        f = eq.pcomp(f, eq.P_S)
    return eq.padd(f, [5])


def predicted_neg(n: int) -> int:
    return 1 if n == 1 else 2 ** (n - 2)


def certify(n: int) -> None:
    from flint import arb, ctx, fmpz_poly
    t1 = time.time()
    old = ctx.prec
    ctx.prec = 1024
    try:
        A, B = (fmpz_poly(c) for c in eq.eq1_pair(n))
        q = fmpz_poly(q_poly(n))
        poles = sorted((r.real for r, _ in A.complex_roots()), key=lambda x: float(x.mid()))
        sg = [q(z) for z in poles]
        decided = all(bool(v > 0) or bool(v < 0) for v in sg)
        neg_idx = [i for i, v in enumerate(sg) if bool(v < 0)]
        zeros = sorted((r.real for r, m in B.complex_roots() if r.imag == 0), key=lambda x: float(x.mid()))
        rootsq = sorted((r for r, _ in q.complex_roots()), key=lambda x: float(x.real.mid()))
        real_q = all(r.imag == 0 for r in rootsq) and len(rootsq) == q.degree() == 2 ** (n - 1)
        alt = len(zeros) == len(poles) - 1 and all(bool(poles[i] < zeros[i]) and bool(zeros[i] < poles[i + 1]) for i in range(len(zeros)))
        # one negative pole in each interval between adjacent roots of q where q < 0 (n >= 2); for n = 1 the ray left of -5
        if n == 1:
            one_each = len(neg_idx) == 1 and bool(poles[neg_idx[0]] < arb(-5))
        else:
            ivs = [(rootsq[2 * i].real, rootsq[2 * i + 1].real) for i in range(len(rootsq) // 2)]
            one_each = len(ivs) == 2 ** (n - 2) and all(
                sum(1 for z in poles if bool(z > a) and bool(z < b)) == 1 for a, b in ivs) and \
                all(any(bool(poles[i] > a) and bool(poles[i] < b) for a, b in ivs) for i in neg_idx)
        ok = decided and len(neg_idx) == predicted_neg(n) and real_q and alt and one_each
        print(f"     n = {n}: {len(poles)} poles, {len(neg_idx)} with q_(n-1) < 0 (predicted {predicted_neg(n)}), zeros of B_n {len(zeros)}", flush=True)
    finally:
        ctx.prec = old
    check(ok, f"ER.2(a) n = {n} (arb, 1024 bits): signs of q_(n-1)(z_j) decided; {predicted_neg(n)} negative, {2**n - predicted_neg(n)} positive; "
          f"q_(n-1) has 2^{n - 1} real roots; B_n zeros alternate with poles; one negative-residue pole per interval ({time.time() - t1:.1f}s)")


def main() -> int:
    t0 = time.time()
    print("Theorem ER: signs of H_n's residues at every level n. Not thrust.", flush=True)
    # ER.2 (c): p^m(-5) = 0 (exact), and p(y) = t > -25/4 has two real roots: min p = p(-5/2) = -25/4 with p(-5/2) < -5 < 0
    pm = [-5]
    for m in range(1, 8):
        pm.append(pm[-1] * (pm[-1] + 5))
    check(all(v == 0 for v in pm[1:]) and Fraction(-5, 2) * Fraction(5, 2) == Fraction(-25, 4),
          "ER.2(c): p^m(-5) = 0 for m >= 1 (zeros of q_k, k < n-1, are never in {p^(n-1) < -5}); min p = -25/4")
    for n in range(1, NCERT + 1):
        certify(n)
    # (b) independent: float Laplacian eigenvalues as poles, exact rational evaluation of q_(n-1)
    for n in range(1, NLAP + 1):
        nv, groups = eq.laplacian_weights(n)
        w = [g for g in groups if g[1] > 1e-12]
        qc = q_poly(n)
        vals = [eq.peval(qc, Fraction(-float(mu))) for mu, _ in w]
        neg = sum(1 for v in vals if v < 0)
        check(len(w) == 2**n and neg == predicted_neg(n) and all(v != 0 for v in vals),
              f"ER.2(b) n = {n} (float Laplacian eigenvalues, exact rationals): {len(w)} poles, {neg} negative H-residues (predicted {predicted_neg(n)})")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all ER checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
