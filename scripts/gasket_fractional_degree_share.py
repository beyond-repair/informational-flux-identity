#!/usr/bin/env python3
"""Degree-gap factorization of the two-neighbor fractional share.

Witness for Theorem AE and Kept Failure AE.1 in GASKET.md.

For the L^α-Dirichlet minimizer u with corner data (1,-1/2,0), write
  w_n      = -(L^α)_{a p} = -(L^α)_{a q}   (equal by altitude reflection),
  diag     = (L^α)_{a a} = sum_j w_{a j},
  λ_n      = 2 w_n / diag                   (two-neighbor degree share),
  W_far    = diag - 2 w_n,
  I_loc    = w_n I_a,   I_far = I_a^α - I_loc,
  ḡ_nbr    = I_a / 2,   ḡ_far = I_far / W_far,
  γ_n      = ḡ_far / ḡ_nbr                  (far-to-neighbor mean-gap ratio),
  ρ_n      = I_loc / I_a^α                  (two-neighbor current share).

Theorem AE:
  (1) sharp ceiling w_n ≤ 8^α / 2  (AB had ≤ 8^α);
  (2) degree-gap factorization  ρ_n = λ_n / (λ_n + (1-λ_n) γ_n);
  (3) resistance identity from AD: I_a^α = 5/(3 R_α), so
        ρ_n = (3/5) w_n I_a R_α = (√21/7) w_n ||F|| R_α;
  (4) for α < log3/log5, ||F||→0 iff ρ→0 iff ((1-λ)/λ) γ → +∞.

Kept Failure AE.1: the law "λ_{n+1}/λ_n ≤ 0.95 for every n, hence λ→0
geometrically and (with any uniform positive lower bound on γ) ρ→0 and
||F(n,0.45)||→0 by AB+AA" is false (ratio already >0.99 by level 3→4).
Not thrust. Does not decide zero vs positive for ||F||.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_limit import (  # noqa: E402
    BC,
    build_gasket,
    combinatorial_laplacian,
    corner_flux_norm,
)

ALPHA = 0.45
W_MIN = ALPHA * (4.0 ** (ALPHA - 1.0))
EIGHT_A = 8.0 ** ALPHA
EIGHT_A_HALF = EIGHT_A / 2.0
SQRT21_5 = math.sqrt(21.0) / 5.0
SQRT21_7 = math.sqrt(21.0) / 7.0

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def spectral_power(L: np.ndarray, alpha: float) -> np.ndarray:
    w, V = np.linalg.eigh(L)
    w = np.clip(w, 0.0, None)
    w[w < 1e-12] = 0.0
    return (V * (w ** alpha)) @ V.T


def level_row(level: int, alpha: float) -> dict:
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    A = spectral_power(L, alpha)
    a = int(corners[0])
    nbrs = sorted(
        {j for (i, j) in edges if i == a} | {i for (i, j) in edges if j == a}
    )
    assert len(nbrs) == 2
    interior = np.setdiff1d(np.arange(N), corners)
    u = np.zeros(N)
    u[corners] = BC
    u[interior] = np.linalg.solve(
        A[np.ix_(interior, interior)], -A[np.ix_(interior, corners)] @ BC
    )
    wa = np.array([-A[a, j] for j in nbrs], dtype=np.float64)
    gaps = np.array([u[a] - u[j] for j in nbrs], dtype=np.float64)
    Ia = float(gaps.sum())
    Iaf = float((A @ u)[a])
    wn = float(wa[0])
    Iloc = float(np.dot(wa, gaps))
    Ifar = Iaf - Iloc
    diag = float(A[a, a])
    Wfar = diag - 2.0 * wn
    g_nbr = Ia / 2.0
    g_far = Ifar / Wfar
    lam = 2.0 * wn / diag
    gamma = g_far / g_nbr
    rho = Iloc / Iaf
    Q = float(u @ A @ u)
    R = (7.0 / 3.0) / Q
    If = (A @ u)[corners]
    I = (L @ u)[corners]
    return {
        "N": N,
        "wn": wn,
        "wa_equal": abs(wa[0] - wa[1]),
        "diag": diag,
        "Wfar": Wfar,
        "lam": lam,
        "gamma": gamma,
        "g_nbr": g_nbr,
        "g_far": g_far,
        "Ia": Ia,
        "Iaf": Iaf,
        "Iloc": Iloc,
        "Ifar": Ifar,
        "rho": rho,
        "phi": 1.0 - rho,
        "Q": Q,
        "R": R,
        "nF": corner_flux_norm(P, corners, I),
        "nF_frac": corner_flux_norm(P, corners, If),
        "lam_max": float(np.linalg.eigvalsh(L).max()),
        "sumI": float(I.sum()),
        "sumIf": float(If.sum()),
    }


def main(through: int = 6) -> int:
    assert through >= 4
    print("Theorem AE: degree-gap factorization of the two-neighbor share")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0).")
    print("Not thrust. Does not decide zero vs positive for ||F||.\n")
    print(
        f"w_min(0.45)={W_MIN:.8f}  8^a/2={EIGHT_A_HALF:.8f}  "
        f"(AB ceiling was 8^a={EIGHT_A:.8f})"
    )

    rows = []
    for n in range(1, through + 1):
        r = level_row(n, ALPHA)
        rows.append(r)
        print(
            f"level {n} N={r['N']} lam={r['lam']:.8f} gamma={r['gamma']:.8f} "
            f"rho={r['rho']:.8f} w_n={r['wn']:.8f} R={r['R']:.8f}",
            flush=True,
        )

    for n, r in enumerate(rows, start=1):
        check(
            r["wa_equal"] < 1e-12,
            f"level {n}: equal neighbor weights (diff {r['wa_equal']:.1e})",
        )
        check(
            W_MIN - 1e-12 <= r["wn"] <= EIGHT_A_HALF + 1e-12,
            f"level {n}: w_min <= w_n={r['wn']:.8f} <= 8^a/2 "
            f"(sharp ceiling; AB had 8^a)",
        )
        check(
            r["diag"] <= EIGHT_A + 1e-12 and r["lam_max"] <= 8.0 + 1e-9,
            f"level {n}: diag={r['diag']:.8f} <= 8^a; lam_max={r['lam_max']:.6f} <= 8",
        )
        check(
            abs(r["diag"] - (2.0 * r["wn"] + r["Wfar"])) < 1e-12 and r["Wfar"] > 0,
            f"level {n}: diag = 2 w_n + W_far with W_far={r['Wfar']:.8f} > 0",
        )
        check(
            0.0 < r["lam"] <= 1.0 + 1e-12,
            f"level {n}: lambda={r['lam']:.8f} in (0,1]",
        )
        # Degree-gap factorization.
        rho_fac = r["lam"] / (r["lam"] + (1.0 - r["lam"]) * r["gamma"])
        check(
            abs(r["rho"] - rho_fac) < 1e-12,
            f"level {n}: rho = lambda/(lambda+(1-lambda)gamma) "
            f"({r['rho']:.12f} vs {rho_fac:.12f})",
        )
        phi_fac = (1.0 - r["lam"]) * r["gamma"] / (
            r["lam"] + (1.0 - r["lam"]) * r["gamma"]
        )
        check(
            abs(r["phi"] - phi_fac) < 1e-12 and abs(r["rho"] + r["phi"] - 1.0) < 1e-12,
            f"level {n}: phi = 1-rho = (1-lambda)gamma/(lambda+(1-lambda)gamma)",
        )
        # Resistance identity from AD: I_a^α = 5/(3 R), Q R = 7/3.
        Iaf_AD = 5.0 / (3.0 * r["R"])
        Q_from_Iaf = (7.0 / 5.0) * r["Iaf"]
        check(
            abs(r["Iaf"] - Iaf_AD) < 1e-10 and abs(r["Q"] - Q_from_Iaf) < 1e-10,
            f"level {n}: I_a^a = 5/(3R) and Q = (7/5) I_a^a "
            f"(AD pairing live: Q R = {r['Q']*r['R']:.12f})",
        )
        check(
            abs(r["Q"] * r["R"] - 7.0 / 3.0) < 1e-10,
            f"level {n}: Q R = 7/3 exactly (AD)",
        )
        rho_AD = (3.0 / 5.0) * r["wn"] * r["Ia"] * r["R"]
        rho_AD2 = SQRT21_7 * r["wn"] * r["nF"] * r["R"]
        check(
            abs(r["rho"] - rho_AD) < 1e-10 and abs(r["rho"] - rho_AD2) < 1e-10,
            f"level {n}: rho = (3/5) w_n I_a R = (sqrt21/7) w_n ||F|| R",
        )
        check(
            abs(r["nF"] - r["Ia"] * SQRT21_5) < 1e-10
            and abs(r["nF_frac"] - r["Iaf"] * SQRT21_5) < 1e-10,
            f"level {n}: Schur audit norms match Ia, Iaf",
        )
        check(
            abs(r["sumI"]) < 1e-8 and abs(r["sumIf"]) < 1e-8,
            f"level {n}: corner-current neutrality",
        )
        check(
            abs(r["Iloc"] - r["wn"] * r["Ia"]) < 1e-12
            and abs(r["Ifar"] - (r["Iaf"] - r["Iloc"])) < 1e-12,
            f"level {n}: I_loc = w_n I_a; I_far = I_a^a - I_loc",
        )

    # Level-2 lock.
    check(abs(rows[1]["nF"] - 1.165909) < 5e-6, "level 2 ||F|| lock 1.165909")

    # Recorded AE table through level 6.
    lam_rec = (0.646893, 0.553135, 0.533143, 0.529482, 0.528861, 0.528759)
    gam_rec = (1.443835, 1.390438, 1.378096, 1.375485, 1.374963, 1.374864)
    for i, r in enumerate(rows[:6]):
        check(
            abs(r["lam"] - lam_rec[i]) < 5e-7 and abs(r["gamma"] - gam_rec[i]) < 5e-7,
            f"level {i+1}: lambda, gamma match recorded AE table",
        )

    # Kept Failure AE.1: no uniform 0.95 geometric wipeout of lambda.
    lams = [r["lam"] for r in rows]
    lam_ratios = [lams[i] / lams[i - 1] for i in range(1, len(lams))]
    print("lambda_ratios", [f"{x:.8f}" for x in lam_ratios])
    check(lam_ratios[0] < 0.90, "level 1→2 lambda ratio < 0.90 (early drop exists)")
    check(
        lam_ratios[2] > 0.99,  # level 3 → 4
        f"Kept Failure AE.1: lambda ratio level 3→4 = {lam_ratios[2]:.8f} > 0.99 "
        "(rejects lambda_{n+1}/lambda_n ≤ 0.95 for every n)",
    )
    if through >= 5:
        check(
            lam_ratios[3] > 0.998,  # level 4 → 5
            f"level 4→5 lambda ratio = {lam_ratios[3]:.8f} > 0.998",
        )
    if through >= 6:
        check(
            lam_ratios[4] > 0.999,  # level 5 → 6
            f"level 5→6 lambda ratio = {lam_ratios[4]:.8f} > 0.999",
        )

    # Sanity: gamma stays O(1) on this range (observation, not a limit claim).
    check(
        all(1.0 < r["gamma"] < 2.0 for r in rows),
        "gamma_n in (1,2) on this range (observation; not a limit theorem)",
    )

    print()
    print("Theorem AE: rho = lambda/(lambda+(1-lambda)gamma); w_n <= 8^a/2;")
    print("rho = (3/5) w_n I_a R_alpha from AD; ||F||->0 iff rho->0 iff ((1-lam)/lam)gam->inf.")
    print("Kept Failure AE.1: no uniform 0.95 geometric wipeout of neighbor-degree share lambda.")
    print("Zero vs positive for combinatorial ||F|| stays OPEN. Not thrust.")
    if failures:
        print(f"\n{len(failures)} check(s) failed")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    through = 6
    if len(sys.argv) > 1:
        through = int(sys.argv[1])
    raise SystemExit(main(through))
