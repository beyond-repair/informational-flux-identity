#!/usr/bin/env python3
r"""Theorem CP: explicit n_0 and explicit constants for CO.5, and the exact rate of CO.5's convergence.

Setting (GASKET.md, Theorems CC, CG, CN, CO): CC's level sums s_n, q_j = l_j/r, and CN.2
    s_(k+4) = a_1 r^k + a_2 l_2^k + a_3 l_3^k + e_k,   |e_k| <= eps_0 (1/1000)^k,   eps_0 <= 2.889e4  (CN.2).
With B = (a_2/a_1) q_2^-4 and B_3 = (a_3/a_1) q_3^-4 (CO.4), write s_n = a_1 r^(n-4) (1 + beta_n),
beta_n = B q_2^n + B_3 q_3^n + eps_n, |eps_n| <= (eps_0/a_1)(rho/r)^(n-4), rho = 1/1000, and
    T_n = log(s_(n+1)/(r s_n)) - B q_2^n (q_2 - 1)     (CO.5's second-order residual).

Theorem CP (checked here; exit 1 on any failure):
  CP.1 [certified; arb, CC's s_4..s_13] For 4 <= n <= 12, T_n has sign (-1)^n, and T_n / q_3^n is enclosed (in [1530, 1735] for n <= 11).
  CP.2 [proved; explicit] For every n >= n_0 = 12,
           T_n / q_3^n = B_3 (q_3 - 1) + (B^2/2)(1 - q_2^2) w^n + R_n,   w = q_2^2/q_3,
       |R_n| <= E_2(n) = P_1 th^n + P_2 |q_2|^n + P_3 |q_3|^n + P_4 (|q_2|^3/|q_3|)^n,  th = rho/(r |q_3|),
       with certified P_i; every base is < 1, so E_2 decreases.  E_1(12) = (B^2/2)(1 - q_2^2)|w|^12 + E_2(12) < |B_3 (q_3 - 1)|.
       Hence (with CP.1) sign(T_n) = (-1)^n for EVERY n >= 4, and c |q_3|^n <= |T_n| <= C |q_3|^n for every n >= 4 (explicit c, C).
  CP.3 [proved; explicit] Every base of E_2(n)/|w|^n is < 1 and E_2(15)/|w|^15 < (B^2/2)(1 - q_2^2).  So for every n >= n_1 = 15
       the deviation T_n/q_3^n - B_3(q_3 - 1) has sign (-1)^n, and
           (T_n/q_3^n - B_3 (q_3 - 1)) / w^n  ->  (B^2/2)(1 - q_2^2) = 32.7 != 0:
       CO.5's convergence has the exact rate q_2^2/|q_3| = l_2^2/(r |l_3|) = 0.73915.
Kept limits: the lambda_4-term sits inside eps_n (CN.2's crude 1e-3 circle); a_4, C_G, circle eigenvalues stay open.  The
limit B_3 (q_3 - 1) is enclosed only in [1630.29, 1636.66]; floats cannot see the 32.7 w^n term beyond n ~ 9.
Assumptions used as stated (not re-run): CG.3's enclosures of r, l_2, l_3 (re-evaluated via CN's spectral_layout), CN.2's
expansion with eps_0 <= 2.889e4 (rounded up here to 2.9e4), CO.4's enclosures of a_1, a_2, a_3.  CC's level sums are re-run.
Textbook: |log(1+x) - x + x^2/2| <= |x|^3/(3(1 - |x|)), x - log(1+x) <= x^2/(2(1-|x|)).  Stage 1 untouched; no W, 0.08, 0.23
or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time

from flint import arb, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_coefficient as cn     # noqa: E402

FAILS: list[str] = []
PREC = 200
N0, N1 = 12, 15
# CO.4 (as stated): projection mids, radii rounded up
CO4 = {"r": ("0.00571256871089104", "1.6e-13"), "l2": ("-4.15763912993700e-5", "4.9e-11"),
       "l3": ("-2.36115967884500e-5", "4.6e-8")}
EPS0 = "2.9e4"          # CN.2: eps_0 <= 2.889e4, rounded up
RHO = arb(1) / 1000     # CN.2's inner circle


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def ball(m, r):
    return arb(m) + arb(0, arb(r).upper())


def constants():
    ctx.prec = PREC
    eig = cn.spectral_layout()
    r, l2, l3 = eig["r"], eig["l2"], eig["l3"]
    a1, a2, a3 = (ball(*CO4[k]) for k in ("r", "l2", "l3"))
    q2, q3 = l2 / r, l3 / r
    B, B3 = a2 / a1 / q2 ** 4, a3 / a1 / q3 ** 4
    k = dict(r=r, q2=q2, q3=q3, a1=a1, B=B, B3=B3, A=abs(B), A3=abs(B3), Q2=abs(q2), Q3=abs(q3),
             e4=arb(EPS0) / arb(a1.lower()), th=RHO / (r * abs(q3)), w=q2 * q2 / q3)
    k["lim"] = B3 * (q3 - 1)
    k["sec"] = B * B / 2 * (1 - q2 * q2)
    check(bool(a1 > 0) and bool(q2 < 0) and bool(q3 < 0) and not B3.contains(0),
          f"CP (inputs): a_1 > 0, q_2 = {q2.str(10)} < 0, q_3 = {q3.str(10)} < 0; B = {B.str(8)}, B_3 = {B3.str(6)}; "
          f"B_3 (q_3 - 1) = {k['lim'].str(6)}, (B^2/2)(1 - q_2^2) = {k['sec'].str(6)}, w = q_2^2/q_3 = {k['w'].str(8)}")
    return k


def data_part(k, sums):
    """CP.1: T_n for 4 <= n <= N0 from CC's certified level sums."""
    ctx.prec = PREC
    r, q2, q3, B = k["r"], k["q2"], k["q3"], k["B"]
    ok, lo, hi, vals = True, None, None, []
    for n in range(4, N0 + 1):
        T = (sums[n + 1] / (r * sums[n])).log() - B * q2 ** n * (q2 - 1)
        ok = ok and (bool(T > 0) if n % 2 == 0 else bool(T < 0))
        v = T / q3 ** n
        vals.append(v)
        if n >= N0:
            continue
        lo = v.lower() if lo is None else min(lo, v.lower())
        hi = v.upper() if hi is None else max(hi, v.upper())
    print("   T_n / q_3^n, n = 4..12: " + ", ".join(v.str(6) for v in vals))
    check(ok, f"CP.1 (arb, CC's s_4..s_13, CO.4's B): sign(T_n) = (-1)^n for 4 <= n <= {N0}, and T_n/q_3^n in "
              f"[{arb(lo).str(6)}, {arb(hi).str(6)}] for 4 <= n < {N0}")
    return arb(lo), arb(hi), vals


def tail_coeffs(k, n0=N0):
    """CP.2: E_2(n) = P1 th^n + P2 |q2|^n + P3 |q3|^n + P4 (|q2|^3/|q3|)^n for n >= n0."""
    ctx.prec = PREC
    A, A3, Q2, Q3, e4, r, th = k["A"], k["A3"], k["Q2"], k["Q3"], k["e4"], k["r"], k["th"]
    # e_m = e4 (rho/r)^(m-4) = E (rho/r)^m with E = e4 (r/rho)^4;  e_m / |q3|^m = E th^m, decreasing
    E = e4 * (r / RHO) ** 4
    Gs = A3 + E * th ** n0                 # g_m = A3|q3|^m + e_m <= Gs |q3|^m       (m >= n0)
    Bs = A + Gs * (Q3 / Q2) ** n0          # b_m = A|q2|^m + g_m <= Bs |q2|^m        (m >= n0)
    bbar = Bs * Q2 ** n0                   # b_m <= bbar < 1/2                       (m >= n0)
    bases = {"th": th, "|q3/q2|": Q3 / Q2, "|q2|": Q2, "|q3|": Q3, "|q2|^3/|q3|": Q2 ** 3 / Q3}
    ok = all(bool(v < 1) for v in bases.values()) and bool(bbar < arb(1) / 2)
    P1 = E * (1 + RHO / r)                                 # (e_n + e_(n+1)) / |q3|^n
    P2 = A * Gs * (1 + Q2 * Q3)                            # B q2^m gamma_m, m = n, n+1
    P3 = Gs * Gs / 2 * (1 + Q3 * Q3)                       # gamma_m^2 / 2
    P4 = Bs ** 3 / (3 * (1 - bbar)) * (1 + Q2 ** 3)        # cubic remainders of log
    check(ok, "CP.2 (arb): bases " + ", ".join(f"{n} = {v.str(5)}" for n, v in bases.items()) + " all < 1; for m >= "
              f"{n0}: |B_3 q_3^m + eps_m| <= {Gs.str(6)} |q_3|^m, |beta_m| <= {Bs.str(6)} |q_2|^m <= {bbar.str(4)} < 1/2")
    P = (P1, P2, P3, P4)
    print("   P_1..P_4 = " + ", ".join(p.str(6) for p in P))
    return P


def E2(k, P, n):
    return P[0] * k["th"] ** n + P[1] * k["Q2"] ** n + P[2] * k["Q3"] ** n + P[3] * (k["Q2"] ** 3 / k["Q3"]) ** n


def conclude(k, P, lo, hi):
    ctx.prec = PREC
    W = abs(k["w"])
    lim_abs = abs(k["lim"])
    E1 = k["sec"] * W ** N0 + E2(k, P, N0)
    E1m1 = k["sec"] * W ** (N0 - 1) + E2(k, P, N0 - 1)
    check(bool(E1 < lim_abs) and bool(k["lim"] > 0),
          f"CP.2: for every n >= n_0 = {N0}, |T_n/q_3^n - B_3(q_3 - 1)| <= E_1(n) <= E_1({N0}) = {E1.str(5)} < "
          f"B_3(q_3 - 1) = {k['lim'].str(6)}, so sign(T_n) = sign(q_3^n) = (-1)^n  (E_1({N0 - 1}) = {E1m1.str(4)}: the "
          f"tail bound alone needs n >= {N0})")
    C = arb(max(hi.upper(), (lim_abs + E1).upper()))
    c = arb(min(lo.lower(), (lim_abs - E1).lower()))
    check(bool(c > 0), f"CP.2: with CP.1, sign(T_n) = (-1)^n for EVERY n >= 4, and {c.str(5)} |q_3|^n <= |T_n| <= "
                       f"{C.str(5)} |q_3|^n for every n >= 4 (explicit)")
    # CP.3: rate of convergence
    bases = {"th/|w|": k["th"] / W, "|q2|/|w|": k["Q2"] / W, "|q3|/|w|": k["Q3"] / W, "|q2|^3/(|q3||w|)": k["Q2"] ** 3 / k["Q3"] / W}
    okb = all(bool(v < 1) for v in bases.values())
    rel = E2(k, P, N1) / W ** N1
    relm1 = E2(k, P, N1 - 1) / W ** (N1 - 1)
    check(okb and bool(rel < k["sec"]),
          "CP.3 (arb): bases " + ", ".join(f"{n} = {v.str(4)}" for n, v in bases.items()) + " all < 1, so E_2(n)/|w|^n "
          f"decreases; E_2({N1})/|w|^{N1} = {rel.str(4)} < (B^2/2)(1 - q_2^2) = {k['sec'].str(6)}  (at n = {N1 - 1}: "
          f"{relm1.str(4)}). So for n >= {N1} the deviation T_n/q_3^n - B_3(q_3 - 1) has sign sign(w^n) = (-1)^n, and "
          f"divided by w^n it tends to (B^2/2)(1 - q_2^2) != 0: CO.5's convergence has the exact rate "
          f"q_2^2/|q_3| = {W.str(10)}")
    return dict(E1=E1, C=C, c=c, rel=rel)


def main() -> int:
    t0 = time.time()
    print("Theorem CP: explicit n_0 and constants in CO.5, exact rate of CO.5's convergence. Not thrust.")
    k = constants()
    sums = cn.level_data(13)
    lo, hi, _ = data_part(k, sums)
    P = tail_coeffs(k)
    conclude(k, P, lo, hi)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS or cn.FAILS:
        print(f"{len(FAILS) + len(cn.FAILS)} FAILED")
        return 1
    print("all Theorem CP checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
