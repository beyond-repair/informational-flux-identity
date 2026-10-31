#!/usr/bin/env python3
"""Theorem CV: the exact rate of inf R_n^(k) -> rho(L_k) for k = 5, ..., 12 (CU's open item "k >= 5"), and Kept Failure
CV.6: |nu_2^(k)/rho(L_k)| is NOT monotone in k (CU's observation "|q^(k)| increases with k" fails between k = 8 and 9).

Setting exactly as in Theorem CU (GASKET.md): s = sqrt(t + 25/4), X = [sqrt5/2, 5/2], sigma_(1,2) = sqrt(15/4 -+ s),
(L_k f)(s) = u_1 f(sigma_1) + u_2 f(sigma_2), u_i = wt_i / (2 sigma_i)^k, h_n = L_k^n 1, R_n = h_(n+1)/h_n on X; L_k acts on
CG's H^2(D), D = {|s - 83/50| < 8/5}.  The routine is CU's script, imported and run unchanged per k, with three parameter
changes only: N = 180 (CU: 120), Neumann power m = 20 (CU/CR: 7), and per-k inner radii rho_in below.  CU's script gained
one guard: for k >= 5 there is no CL.3 interval, so its consistency check is skipped (k = 2, 3, 4 behave exactly as before).

Theorem CV (checked here; exit 1 on any failure).  For each k = 5, ..., 12:
  CV.1-CV.3 [CU.1-CU.3's proofs] ||L_k (I - P_180)|| <= eps_180 (arb), Krawczyk eigenpairs of M_180 at rho and nu_2, the
       Neumann bound ||A^20||_F < rho_in^20, and CG.3's homotopy: L_k on H^2(D) has exactly two eigenvalues in |z| > rho_in,
       real and simple, rho_H = rho(L_k) (point evaluation; min_X A > 0) and nu_2^(k) < 0, each to +- 6e-20.
  CV.4-CV.5 [CU.4-CU.5's proofs] h_n = rho^n A + nu_2^n B + E_n, Gamma_+- = max_X(+-B/A) enclosed, and for every n >= n_0(k)
       rho(1-q)|q|^n (Gamma_par - theta_n) <= rho(L_k) - inf_X R_n^(k) <= rho(1-q)|q|^n (Gamma_par + theta_n):
       the exact rate is |q^(k)| = |nu_2^(k)|/rho(L_k).
  CV.6 [Kept Failure, interval witness] With CT.3's q^(1) and CU.3's q^(2..4) (used as stated) and CV.3's q^(5..12):
       |q^(k)| strictly increases for k = 1..8 and strictly decreases for k = 8..12; witness |q^(8)| - |q^(9)| > 7e-4.
       So max_(k <= 12) |q^(k)| = |q^(8)| = 0.358835329184..., and CU's observation does not extend.
  CV.7 [proved, finite range] 3^k rho(L_k) strictly decreases for k = 1..12 and stays > 1/3.
Observation (no claim): s* = 3/2 is the fixed point of sigma_1, with sigma_1'(s*) = -1/3 and wt_1(s*) = 1/3; the
single-branch heuristic suggests 3^k rho(L_k) -> 1/3 and q^(k) -> -1/3 as k -> infinity (|q^(12)| - 1/3 = 0.01889).
Neither limit is proved here.
Kept limits: only k <= 12; n < n_0(k) not covered (n_0 grows with k); nothing uniform in k; no k -> infinity limit.
Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction as Fr

from flint import arb

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_infrate_k as cu    # noqa: E402

FAILS: list[str] = []
ALL: dict = {}
N_CV = 180
M_POW_CV = 20
# k: (float seeds for rho, nu_2 (Krawczyk certifies them), inner radius rho_in ~ sqrt(|nu_2 nu_3|))
SEED = {
    5: ("1.479255923963e-03", "-5.161358132105e-04", Fr(1, 5000)),
    6: ("4.853393939727e-04", "-1.723958039714e-04", Fr(1, 12500)),
    7: ("1.597483711719e-04", "-5.721115355469e-05", Fr(1, 35000)),
    8: ("5.271694947509e-05", "-1.891670391848e-05", Fr(1, 100000)),
    9: ("1.743290380131e-05", "-6.242980983375e-06", Fr(1, 300000)),
    10: ("5.774526639771e-06", "-2.058808531062e-06", Fr(1, 900000)),
    11: ("1.915332062748e-06", "-6.789469804583e-07", Fr(1, 2600000)),
    12: ("6.359658135257e-07", "-2.240022214049e-07", Fr(1, 7800000)),
}
KS = [int(x) for x in os.environ.get("CV_KS", ",".join(str(k) for k in SEED)).split(",")]
# CT.3 and CU.3's certified q = nu_2/rho (midpoint, radius), used as stated for CV.6/CV.7
Q_PRIOR = {1: ("-0.244080070117426", "4e-16"), 2: ("-0.289567063896138", "2.7e-16"),
           3: ("-0.318994978288970", "1.9e-16"), 4: ("-0.337629300012356", "4.2e-16")}
RHO_PRIOR = {1: ("0.1343549829703998131", "2.9e-20"), 2: ("0.0430703112227101271", "1.3e-20"),
             3: ("0.0139192547691468646", "2.1e-20"), 4: ("0.0045264690852562123", "2.8e-20")}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def conf(k):
    r, nu, rin = SEED[k]
    rf, nf = float(r), abs(float(nu))
    return dict(seeds=[r, nu], rho_in=rin, proj=[Fr(1, round(3 / rf)), Fr(1, round(2 / (nf - float(rin))))],
                cl3=None, qcap="0.36")


def main() -> int:
    print("Theorem CV: the exact rate of inf R_n^(k) -> rho(L_k) for k = 5..12; Kept Failure CV.6: |q^(k)| not monotone. "
          "Not thrust.", flush=True)
    t0 = time.time()
    ALL.clear()
    old = (cu.N, cu.cr.M_POW)
    cu.N, cu.cr.M_POW = N_CV, M_POW_CV
    try:
        for k in KS:
            cu.CONF[k] = conf(k)
            nf0 = len(cu.FAILS) + len(cu.cr.FAILS)
            rc = cu.run_k(k)
            ALL[k] = dict(cu.RESULT)
            check(rc == 0 and len(cu.FAILS) + len(cu.cr.FAILS) == nf0 and ALL[k].get("n0") is not None,
                  f"CV.1-CV.5 (k = {k}): every CU check passed at N = {N_CV}, m = {M_POW_CV}, rho_in = {SEED[k][2]}")
    finally:
        cu.N, cu.cr.M_POW = old
    for k in list(cu.CONF):
        if k >= 5:
            del cu.CONF[k]
    if any(ALL[k].get("q") is None for k in KS):
        check(False, "CV: some k not certified; CV.6/CV.7 skipped")
    else:
        rows = []
        for k in KS:
            a = ALL[k]
            rows.append(f"k={k}: rho = {a['rho'].str(17)}, nu_2 = {a['nu'].str(17)}, |q| = {abs(a['q']).str(12)}, "
                        f"Gamma_+- = {a['Gp'].str(7)}, {a['Gm'].str(7)}, c_+- = {a['cp'].str(7)}, {a['cm'].str(7)}, "
                        f"n_0 = {a['n0']}, even [{a['ev'][0].str(5)}, {a['ev'][1].str(5)}] n >= {a['ev'][2]}, "
                        f"odd [{a['od'][0].str(5)}, {a['od'][1].str(5)}] n >= {a['od'][2]}")
        print("   CV summary:\n   " + "\n   ".join(rows), flush=True)
        if KS == list(SEED):
            Q = {k: arb(m) + arb(0, r) for k, (m, r) in Q_PRIOR.items()}
            RH = {k: arb(m) + arb(0, r) for k, (m, r) in RHO_PRIOR.items()}
            for k in KS:
                Q[k], RH[k] = ALL[k]["q"], ALL[k]["rho"]
            aq = {k: abs(Q[k]) for k in Q}
            inc = all(bool(aq[k + 1] > aq[k]) for k in range(1, 8))
            dec = all(bool(aq[k + 1] < aq[k]) for k in range(8, 12))
            gap = aq[8] - aq[9]
            check(inc and dec and bool(gap > arb("7e-4")),
                  "CV.6 (Kept Failure, interval witness): |q^(k)| strictly increases for k = 1..8 and strictly decreases for "
                  f"k = 8..12: " + ", ".join(f"{aq[k].str(10)}" for k in sorted(aq)) + f"; |q^(8)| - |q^(9)| = {gap.str(6)} "
                  "> 7e-4.  CU's observation '|q^(k)| increases with k' fails at k = 9; max_(k<=12) |q^(k)| = |q^(8)|")
            t3 = {k: RH[k] * arb(3) ** k for k in RH}
            mono = all(bool(t3[k + 1] < t3[k]) for k in range(1, 12)) and all(bool(t3[k] > arb(1) / 3) for k in t3)
            check(mono, "CV.7: 3^k rho(L_k) strictly decreases for k = 1..12 and stays > 1/3: "
                        + ", ".join(t3[k].str(8) for k in sorted(t3)))
            ALL["aq"], ALL["t3"], ALL["gap"] = aq, t3, gap
            print("   observation (no claim): sigma_1(3/2) = 3/2, sigma_1'(3/2) = -1/3, wt_1(3/2) = 1/3; the single-branch "
                  "heuristic suggests 3^k rho(L_k) -> 1/3 and q^(k) -> -1/3; |q^(12)| - 1/3 = "
                  f"{(aq[12] - arb(1) / 3).str(6)}.  Not proved.", flush=True)
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(cu.FAILS) + len(cu.cr.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem CV checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
