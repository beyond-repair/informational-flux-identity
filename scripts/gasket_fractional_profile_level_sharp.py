#!/usr/bin/env python3
r"""Theorem CS: sharp constants for the level-sum rate, and an explicit (sharp) Birkhoff constant C_G in CF.3.

Setting (GASKET.md, Theorems CC, CF, CG, CN, CO, CP): CC's level sums s_n, q_j = l_j/r, B = (a_2/a_1) q_2^-4 (CO.4),
    D_n = log(s_(n+1)/(r s_n)) / q_2^n,     L = B (q_2 - 1),     T_n = log(s_(n+1)/(r s_n)) - B q_2^n (q_2 - 1)  (CP),
so D_n = L + U_n with U_n = T_n / q_2^n.

Theorem CS (checked here; exit 1 on any failure):
  CS.1 [proved from CP.2] sign T_n = (-1)^n = sign q_2^n for every n >= 4, so U_n = |T_n|/|q_2|^n > 0 and, by CP.2's
       c |q_3|^n <= |T_n| <= C |q_3|^n (c = 1360.9, C = 1906.1),  c |q_3/q_2|^n <= D_n - L <= C |q_3/q_2|^n  for every n >= 4.
  CS.2 [proved; arb] U_(n+1)/U_n <= (C/c) |q_3/q_2| < 1, so D_n is STRICTLY DECREASING for n >= 4, D_n -> L, and D_n > L.
       Hence, for every n >= 4,
           L |q_2|^n < |log(s_(n+1)/(r s_n))| <= D_4 |q_2|^n,
       with L = B(q_2 - 1) = 9.6274 +- 3.6e-5 (not attained, sharp as the infimum) and D_4 = 14.59509785 +- 5.9e-10
       (attained at n = 4).
       These are the best constants in CN.4's two-sided bound (CN.4: 14.60 for n >= 4, 9.626 only for n >= 13).
       By CS.1, D_n - L has the exact rate |q_3/q_2| = |l_3/l_2| = 0.23343938, and (D_n - L)/(q_3/q_2)^n -> B_3(q_3 - 1)
       by CO.5.
  CS.3 [proved; arb] CF.3's level-sum inequality |log(s_(k+1)/(r s_k))| <= C_G theta^(k/5) (k >= 4) holds with the
       explicit sharp constant C_G* = |log(s_5/(r s_4))| theta^(-4/5) (equality at k = 4), since
       |log(...)| theta^(-k/5) = D_k (|q_2|/theta^(1/5))^k is a product of positive strictly decreasing factors.
       With CF.2's certified theta = tanh(Delta/4) = 0.5441333187: C_G* = 0.02105127810 +- 2.1e-12; for the rounded
       theta-bar = 0.5442 the same argument gives the sharp constant 0.0210492 (< C_G*, since theta-bar > theta).
Data consistency (CC's s_4..s_13): D_4 > D_5 > ... as far as the enclosures separate.
Assumptions used as stated (not re-run): CP.2's constants c = 1360.9, C = 1906.1 (CP's script is imported and re-run
unchanged for its inputs: CG.3's r, l_2, l_3, CO.4's a_1, a_2, a_3), CO.5 for the rate statement, CF.2's theta (CF's
cone_constants re-run unchanged).  CC's level sums re-run.  Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW
enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time

from flint import arb, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_threshold as cp     # noqa: E402
import gasket_fractional_profile_level_rate as cf          # noqa: E402

FAILS: list[str] = []
PREC = 200
C_LO, C_HI = "1360.9", "1906.1"      # CP.2 (as stated)
NDATA = 12                           # CC's s_4..s_13 give D_4..D_12
RESULT: dict = {}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def sharp_rate(k, sums):
    """CS.1, CS.2: sharp two-sided constants for |log(s_(n+1)/(r s_n))| / |q_2|^n."""
    ctx.prec = PREC
    r, q2, q3, B = k["r"], k["q2"], k["q3"], k["B"]
    L = B * (q2 - 1)
    c, C = arb(C_LO), arb(C_HI)
    rho = abs(q3 / q2)
    contr = C / c * rho
    check(bool(L > 0) and bool(contr < 1),
          f"CS.1/CS.2 (arb): L = B(q_2 - 1) = {L.str(10)} > 0; |q_3/q_2| = {rho.str(8)}; U_(n+1)/U_n <= (C/c)|q_3/q_2| = "
          f"{contr.str(6)} < 1 (CP.2: c = {C_LO}, C = {C_HI}, sign T_n = (-1)^n = sign q_2^n, so U_n = T_n/q_2^n > 0)")
    D = {}
    for n in range(4, NDATA + 1):
        D[n] = (sums[n + 1] / (r * sums[n])).log() / q2 ** n
    D4 = D[4]
    # D_n = L + U_n, U_n in [c rho^n, C rho^n]: check consistency at n = 4 and the data decrease where decided
    U4 = D4 - L
    check(bool(U4 > c * rho ** 4) and bool(U4 < C * rho ** 4),
          f"CS.1 (consistency, arb): D_4 - L = {U4.str(8)} in [c rho^4, C rho^4] = [{(c * rho ** 4).str(6)}, {(C * rho ** 4).str(6)}]")
    decided = [n for n in range(4, NDATA) if bool(D[n] > D[n + 1])]
    undecided = [n for n in range(4, NDATA) if not bool(D[n] > D[n + 1])]
    bad = [n for n in range(4, NDATA) if bool(D[n] < D[n + 1])]
    print("   D_n, n = 4..12: " + ", ".join(v.str(8) for v in D.values()))
    check(not bad, f"CS.2 (data, CC's s_4..s_13): D_n > D_(n+1) certified directly for n in {decided}; not separated by the "
                   f"enclosures (decided by CS.2's proof only) for n in {undecided}; no reversal")
    check(bool(D4 > L) and bool(D4 > arb(14)) and bool(D4 < arb("14.60")),
          f"CS.2: for every n >= 4, L |q_2|^n < |log(s_(n+1)/(r s_n))| <= D_4 |q_2|^n with L = {L.str(10)} (sharp infimum, "
          f"not attained) and D_4 = log(s_5/(r s_4))/q_2^4 = {D4.str(10)} (attained at n = 4); CN.4 had 14.60 and 9.626 (n >= 13)")
    return dict(L=L, D=D, D4=D4, rho=rho, contr=contr)


def sharp_cg(k, sums, D4):
    """CS.3: the sharp constant in CF.3's level-sum bound."""
    ctx.prec = PREC
    cons = cf.cone_constants()
    theta = cons["theta"]
    q2 = k["q2"]
    t5 = (theta.log() / 5).exp()
    base = abs(q2) / t5
    log4 = abs((sums[5] / (k["r"] * sums[4])).log())
    CG = log4 / (theta.log() * 4 / 5).exp()
    tbar = arb(cf.THETA_MAX)
    CGbar = log4 / (tbar.log() * 4 / 5).exp()
    check(bool(base < 1) and bool(CG > 0),
          f"CS.3 (arb): |q_2|/theta^(1/5) = {base.str(8)} < 1, so |log(s_(k+1)/(r s_k))| theta^(-k/5) = D_k (|q_2|/theta^(1/5))^k "
          f"strictly decreases for k >= 4; sharp C_G* = |log(s_5/(r s_4))| theta^(-4/5) = {CG.str(10)} with CF.2's theta = "
          f"{theta.str(10)}; with the rounded theta < {cf.THETA_MAX}: C_G <= {arb(CGbar.upper()).str(6)} (sharp for theta-bar)")
    return dict(theta=theta, CG=CG, CGbar=CGbar, log4=log4, base=base)


def main() -> int:
    t0 = time.time()
    print("Theorem CS: sharp constants for the level-sum rate and an explicit sharp C_G in CF.3. Not thrust.")
    k = cp.constants()
    sums = cp.cn.level_data(13)
    res = sharp_rate(k, sums)
    res.update(sharp_cg(k, sums, res["D4"]))
    RESULT.update(res)
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(cp.FAILS) + len(cp.cn.FAILS) + len(cf.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem CS checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
