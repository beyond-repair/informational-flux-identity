#!/usr/bin/env python3
"""Theorem EQ: EP.1's pole structure for EVERY level n.  EP kept "EP.1's structure (degrees 2^n, simple real negative poles) is
certified for n <= 6 only, not proved for all n"; EP's open item "a proof of EP.1's pole structure for every n".

  EQ.1 (THEOREM, proof for every n >= 0; AQ (1) as stated).  Put p(s) = s(s+5), A_0 = s + 3, B_0 = 1 and
         A_(n+1)(s) = [(s + 6) A_n(p(s)) - 2 s (s + 5) B_n(p(s))] / (s + 2),     B_(n+1)(s) = (s + 5) B_n(p(s)).
       Then the division is exact in Z[s], A_n and B_n are monic and coprime, r_n = B_n/A_n (so this is the gcd-reduced form EP
       computes), deg A_n = 2^n, deg B_n = 2^n - 1, B_n = prod_(k=0)^(n-1) (p^k(s) + 5), A_n(0) = 3^(n+1), B_n(0) = 5^n,
       A_n(-6) = -3 B_n(-6) != 0 (r_n(-6) = -1/3).  Proof by induction with D = (s+6) A_n(p) - 2s(s+5) B_n(p):
       D(-2) = 4 (A_n(-6) + 3 B_n(-6)) = 0;  D(-5) = A_n(0);  D(-6) = -12 B_n(6);  D(0) = 6 A_n(0) (p(-2) = -6, p(-5) = p(0) = 0,
       p(-6) = 6); B_n(6) = prod (p^k(6) + 5) > 0.  A common root z of A_(n+1), B_(n+1): z = -5 gives A_(n+1)(-5) = -A_n(0)/3 != 0;
       z = -2 gives B_(n+1)(-2) = 3 B_n(-6) != 0; otherwise B_n(p(z)) = 0, so A_n(p(z)) != 0 and D(z) = (z+6) A_n(p(z)) = 0 forces
       z = -6, B_n(6) = 0: impossible.  Degrees: deg (s+6) A_n(p) = 2^(n+1) + 1 > deg s(s+5) B_n(p) = 2^(n+1).
  EQ.2 (THEOREM, every n >= 0; AQ's setting r_n = 1/2 v^T (L_n + s)^-1 v, L_n the graph Laplacian of the connected G_n,
       v = e_a - e_b).  Grouping an orthonormal eigenbasis by distinct eigenvalues mu, r_n = sum_mu m_mu/(s + mu), m_mu = |P_mu v|^2/2;
       v is orthogonal to ker L_n = constants, so only mu > 0 with m_mu > 0 occur, and this sum is already reduced.  By uniqueness of
       the reduced form and EQ.1: A_n = prod (s + mu) over exactly 2^n distinct eigenvalues mu > 0 of L_n with P_mu v != 0.  So for
       every n: A_n has 2^n simple real roots z_j = -mu_j < 0, residues c_j = B_n(z_j)/A_n'(z_j) = m_mu > 0, sum_j c_j = |v|^2/2 = 1,
       sum_j c_j/w_j = r_n(0) = (5/3)^n/3.  H_n's integrand r_n eta_1(p^(n-1)(s)) = prod_(k=0)^(n-2) (p^k + 5) / (2 A_n) (EQ.1's
       product; eta_1 = 1/(2(z+5)), EO.1) is reduced with denominator 2 A_n: the same 2^n simple real negative poles (residues
       c_j/(2(p^(n-1)(z_j) + 5)), of mixed sign: n = 1, z_2 = -(7 + sqrt13)/2 < -5).  Hence EP.1's residue formulas for J_n, H_n,
       om_n, D'_n hold at every level n >= 1.
  EQ.3 (checks / attempts to break): EQ.1's recursion equals EP's gcd-reduced rn_exact for n = 1..7 and its identities hold exactly
       for n = 0..NEX; h_exact's reduced pair is (prod_(k<=n-2)(p^k+5), 2 A_n) for n = 1..7; certified isolation (arb, 1024 bits) at n = 7, 8
       (degrees 128, 256): 2^n simple real negative roots, residues > 0, sum c_j = 1 and sum c_j/w_j = (5/3)^n/3 to 1e-30; the gasket
       Laplacian (floats, n = 1..5; 6, 15, 42, 123, 366 vertices) has 4, 8, 18, 36, 70 distinct eigenvalues of which exactly 2^n carry
       v-weight (the rest < 1e-20), and these are -roots of A_n.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])

FAILS: list[str] = []
NEX = 10          # EQ.3: exact identities for n = 0..NEX (deg A_10 = 1024)
NEP = 7           # EQ.3: comparison with EP's rn_exact / h_exact, certified isolation at n = NEP
NCERT = (7, 8)    # EQ.3: certified isolation
NLAP = 5          # EQ.3: gasket Laplacian floats


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- pure-Python integer polynomials (CI; coefficient lists, low first)
def padd(a, b):
    n = max(len(a), len(b))
    out = [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def pmul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return out


def pcomp(a, q):
    """a(q(s)) by Horner."""
    out = [a[-1]]
    for c in reversed(a[:-1]):
        out = padd(pmul(out, q), [c])
    return out


def pdiv_linear(a, r):
    """Divide a by (s - r) (synthetic division); returns (quotient, remainder)."""
    hi = list(reversed(a))
    q = [hi[0]]
    for c in hi[1:]:
        q.append(c + r * q[-1])
    rem = q.pop()
    return list(reversed(q)), rem


def peval(a, x):
    v = 0
    for c in reversed(a):
        v = v * x + c
    return v


P_S = [0, 5, 1]   # p(s) = s^2 + 5 s


def eq1_pair(n: int):
    """EQ.1's recursion: returns (A_n, B_n) as integer coefficient lists; asserts the division by s + 2 is exact."""
    A, B = [3, 1], [1]
    for _ in range(n):
        Ap, Bp = pcomp(A, P_S), pcomp(B, P_S)
        D = padd(pmul([6, 1], Ap), pmul([0, -10, -2], Bp))
        A, rem = pdiv_linear(D, -2)
        if rem != 0:
            raise ArithmeticError("EQ.1: s + 2 does not divide D")
        B = pmul([5, 1], Bp)
    return A, B


def b_product(n: int):
    """EQ.1: prod_(k=0)^(n-1) (p^k(s) + 5)."""
    out, f = [1], [0, 1]
    for _ in range(n):
        out = pmul(out, padd(f, [5]))
        f = pcomp(f, P_S)
    return out


def eq1_identities(n: int) -> dict:
    A, B = eq1_pair(n)
    return {"degA": len(A) - 1, "degB": len(B) - 1, "monic": A[-1] == 1 and B[-1] == 1, "A0": A[0], "B0": B[0],
            "m6": peval(A, -6) + 3 * peval(B, -6), "Am6": peval(A, -6), "Bprod": B == b_product(n), "B6pos": peval(B, 6) > 0,
            "poscoef": all(c > 0 for c in A) and all(c > 0 for c in B)}


def laplacian_weights(n: int, tol: float = 1e-7):
    """EQ.3 floats: distinct eigenvalues of L_n with their v-weights m_mu = |P_mu v|^2 / 2."""
    import gasket_fractional_corner_exact_levels as ep
    L, v = ep.gasket_laplacian(n)
    lam, phi = np.linalg.eigh(L)
    c = 0.5 * (phi.T @ v) ** 2
    groups: list = []
    for lv, cv in sorted(zip(lam, c)):
        if groups and abs(lv - groups[-1][0]) < tol:
            groups[-1][1] += cv
        else:
            groups.append([lv, cv])
    return len(v), groups


def certify_level(nc: int) -> None:
    """EQ.3: certified isolation (arb_fmpz_poly_complex_roots, 1024 bits) of A_nc, residues, and the two moment identities."""
    from flint import acb_poly, arb, ctx, fmpz_poly
    t1 = time.time()
    old = ctx.prec
    ctx.prec = 1024
    try:
        A, B = (fmpz_poly(c) for c in eq1_pair(nc))
        rts = A.complex_roots()
        real_neg = all(r.imag == 0 and bool(r.real < 0) and m == 1 for r, m in rts)
        N, dA = acb_poly([int(c) for c in B.coeffs()]), acb_poly([int(c) for c in A.derivative().coeffs()])
        cs = [(N(r) / dA(r)).real for r, _ in rts]
        pos = all(bool(c > 0) for c in cs)
        tot = sum(cs, arb(0))
        tw = sum((c / (-r.real) for c, (r, _) in zip(cs, rts)), arb(0))
        target = arb(5) ** nc / arb(3) ** nc / 3
        ok = len(rts) == 2**nc and real_neg and pos and bool(abs(tot - 1) < arb("1e-30")) and bool(abs(tw - target) < arb("1e-30"))
        print(f"     n = {nc}: {len(rts)} roots in [{min(float(r.real.mid()) for r, _ in rts):.6g}, "
              f"{max(float(r.real.mid()) for r, _ in rts):.6g}], min residue {min(float(c.mid()) for c in cs):.3e}", flush=True)
    finally:
        ctx.prec = old
    check(ok, f"EQ.3 n = {nc} (arb, 1024 bits): 2^{nc} simple real roots < 0, residues > 0, sum c_j = 1, sum c_j/w_j = (5/3)^{nc}/3 "
          f"({time.time() - t1:.1f}s)")


def main() -> int:
    t0 = time.time()
    print("Theorem EQ: EP.1's pole structure at every level n. Not thrust.", flush=True)
    import sympy as sp
    s = sp.symbols("s")
    Af, Bf = sp.Function("A"), sp.Function("B")
    p = lambda x: x * (x + 5)  # noqa: E731
    D = lambda x: (x + 6) * Af(p(x)) - 2 * x * (x + 5) * Bf(p(x))  # noqa: E731
    check(sp.simplify(D(sp.Integer(-2)) - 4 * (Af(-6) + 3 * Bf(-6))) == 0 and sp.simplify(D(sp.Integer(-5)) - Af(0)) == 0
          and sp.simplify(D(sp.Integer(-6)) + 12 * Bf(6)) == 0 and sp.simplify(D(sp.Integer(0)) - 6 * Af(0)) == 0,
          "EQ.1 (sympy, generic A, B): D(-2) = 4 (A(-6) + 3 B(-6)), D(-5) = A(0), D(-6) = -12 B(6), D(0) = 6 A(0)")
    x = sp.symbols("x")
    phi = (2 + s) * (5 + s) * x / ((6 + s) - 2 * s * (5 + s) * x)
    Ab, Bb = sp.symbols("Ab Bb")
    check(sp.simplify(phi.subs(x, Bb / Ab) - (s + 2) * (s + 5) * Bb / ((s + 6) * Ab - 2 * s * (s + 5) * Bb)) == 0,
          "EQ.1 (sympy): Phi_s(B(p)/A(p)) = (s+2)(s+5) B(p) / ((s+6) A(p) - 2s(s+5) B(p)) (AQ (1))")
    # ---- EQ.1 exact, n = 0..NEX (pure Python integers)
    for n in range(NEX + 1):
        t1 = time.time()
        I = eq1_identities(n)
        check(I["degA"] == 2**n and I["degB"] == 2**n - 1 and I["monic"] and I["A0"] == 3 ** (n + 1) and I["B0"] == 5**n
              and I["m6"] == 0 and I["Am6"] != 0 and I["Bprod"] and I["B6pos"] and I["poscoef"],
              f"EQ.1 n = {n}: s + 2 | D exactly; A_n, B_n monic, deg {I['degA']}, {I['degB']}; A_n(0) = 3^{n + 1}, B_n(0) = 5^{n}, "
              f"A_n(-6) = -3 B_n(-6) != 0; B_n = prod (p^k + 5); all coefficients > 0 ({time.time() - t1:.2f}s)")
    # ---- EQ.3 exact (python-flint): coprimality, squarefree, EP's gcd-reduced pair, h_exact
    from flint import fmpz_poly
    import gasket_fractional_corner_exact_levels as ep
    for n in range(NEX + 1):
        A, B = (fmpz_poly(c) for c in eq1_pair(n))
        ok = A.gcd(B).degree() == 0 and A.gcd(A.derivative()).degree() == 0
        if 1 <= n <= NEP:
            Be, Ae = ep.rn_exact(n)
            hn, hd = ep.h_exact(n)
            ok = ok and Ae == A and Be == B and hd == 2 * A and hn == fmpz_poly(b_product(n - 1))
        check(ok, f"EQ.3 n = {n}: gcd(A_n, B_n) = 1, gcd(A_n, A_n') = 1" + (" ; equals EP's rn_exact; h_exact = (prod_(k<=n-2)"
              "(p^k+5), 2 A_n)" if 1 <= n <= NEP else ""))
    # ---- EQ.2: H residues are of mixed sign at n = 1 (so EQ.2 claims no sign for them)
    z2 = -(7 + sp.sqrt(13)) / 2
    check(bool(z2 + 5 < 0) and bool((-(7 - sp.sqrt(13)) / 2) + 5 > 0),
          "EQ.2: n = 1, z_2 = -(7 + sqrt13)/2 < -5 < z_1: H_1's residues c_j/(2(z_j + 5)) have opposite signs")
    # ---- EQ.3 certified isolation at n = NEP
    for nc in NCERT:
        certify_level(nc)

    # ---- EQ.3 floats: the gasket Laplacian
    for n in range(1, NLAP + 1):
        nv, groups = laplacian_weights(n)
        w = [g for g in groups if g[1] > 1e-12]
        zero = max((g[1] for g in groups if g[1] <= 1e-12), default=0.0)
        A, _ = eq1_pair(n)
        fa = fmpz_poly(A)
        rr = sorted(float(r.real.mid()) for r, _ in fa.complex_roots())
        mus = sorted(-g[0] for g in w)
        dev = max(abs(a - b) for a, b in zip(rr, mus)) if len(rr) == len(mus) else float("inf")
        check(nv == (3 ** (n + 1) + 3) // 2 and len(w) == 2**n and zero < 1e-20 and min(g[0] for g in w) > 1e-6 and dev < 1e-7,
              f"EQ.3 (floats) n = {n}: L_n ({nv} vertices) has {len(groups)} distinct eigenvalues, exactly {len(w)} = 2^{n} carry v-weight "
              f"(min {min(g[1] for g in w):.2e}; others <= {zero:.1e}), and they are -roots of A_n (max dev {dev:.1e})")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EQ checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
