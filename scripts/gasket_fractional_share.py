#!/usr/bin/env python3
"""Combinatorial-to-fractional corner share on the gasket.

Witness for Theorem AB and Kept Failure AB.1 in GASKET.md.

For the L^α-Dirichlet minimizer u with corner data (1,-1/2,0), write
  I_a      = combinatorial corner current at the value-1 corner,
  I_a^α    = (L^α u)_a  (fractional corner current),
  σ_n      = I_a / I_a^α = ||F|| / ||F_frac||,
  w_n      = -(L^α)_{a p} = -(L^α)_{a q}  (equal by altitude reflection),
  I_loc    = w_n I_a  (two-neighbor piece of the fractional current),
  ρ_n      = I_loc / I_a^α  (two-neighbor share).

Theorem AB: σ_n = ρ_n / w_n, with w_min(α) ≤ w_n ≤ 8^α, hence
  ρ_n / 8^α ≤ σ_n ≤ ρ_n / w_min ≤ 1/w_min.
For α < 1/2, Theorem AA floors ||F_frac|| away from 0, so
  ||F|| → 0  iff  σ_n → 0  iff  ρ_n → 0,
and liminf ||F|| > 0 iff liminf ρ_n > 0. Zero vs positive stays OPEN;
this only names the remaining route after AA.1.

Kept Failure AB.1: the law "ρ_{n+1}/ρ_n ≤ 0.95 for every n, hence ρ→0
geometrically and ||F||→0 by AB+AA" is false (ratio already >0.99 by
level 4→5). Not thrust.
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
SQRT21_5 = math.sqrt(21.0) / 5.0

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


def level_share(level: int, alpha: float) -> dict:
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
    Iloc = float(np.dot(wa, gaps))
    If = (A @ u)[corners]
    I = (L @ u)[corners]
    diag = float(A[a, a])
    return {
        "N": N,
        "wn": float(wa[0]),
        "wa_equal": abs(wa[0] - wa[1]),
        "Ia": Ia,
        "Iaf": Iaf,
        "Iloc": Iloc,
        "rho": Iloc / Iaf,
        "sigma": Ia / Iaf,
        "nF": corner_flux_norm(P, corners, I),
        "nF_frac": corner_flux_norm(P, corners, If),
        "diag": diag,
        "lam_max": float(np.linalg.eigvalsh(L).max()),
        "sumI": float(I.sum()),
        "sumIf": float(If.sum()),
    }


def main(through: int = 6) -> int:
    assert through >= 4
    print("Theorem AB: combinatorial/fractional corner share factorization")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0).")
    print("Not thrust. Does not decide zero vs positive for ||F||.\n")
    print(f"w_min(0.45)={W_MIN:.8f}  8^a={EIGHT_A:.8f}  1/w_min={1.0/W_MIN:.8f}")

    rows = []
    for n in range(1, through + 1):
        r = level_share(n, ALPHA)
        rows.append(r)
        print(
            f"level {n} N={r['N']} w_n={r['wn']:.8f} rho={r['rho']:.8f} "
            f"sigma={r['sigma']:.8f} ||F||={r['nF']:.8f} ||F_frac||={r['nF_frac']:.8f}",
            flush=True,
        )

    # Reflection: equal neighbor weights; factorization; sandwich.
    for n, r in enumerate(rows, start=1):
        check(
            r["wa_equal"] < 1e-12,
            f"level {n}: altitude reflection gives equal neighbor weights "
            f"(diff {r['wa_equal']:.1e})",
        )
        check(
            abs(r["sigma"] - r["rho"] / r["wn"]) < 1e-12,
            f"level {n}: sigma = rho/w_n ({r['sigma']:.10f} vs {r['rho']/r['wn']:.10f})",
        )
        check(
            W_MIN - 1e-12 <= r["wn"] <= EIGHT_A + 1e-12,
            f"level {n}: w_min <= w_n={r['wn']:.8f} <= 8^a",
        )
        lo = r["rho"] / EIGHT_A
        hi = r["rho"] / W_MIN
        check(
            lo - 1e-12 <= r["sigma"] <= hi + 1e-12 and r["sigma"] <= 1.0 / W_MIN + 1e-12,
            f"level {n}: rho/8^a={lo:.6f} <= sigma={r['sigma']:.6f} "
            f"<= rho/w_min={hi:.6f} <= 1/w_min",
        )
        check(
            0.0 < r["rho"] <= 1.0 + 1e-12,
            f"level {n}: rho={r['rho']:.8f} in (0,1]",
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
            r["diag"] + 1e-12 >= 2.0 * r["wn"] and r["lam_max"] <= 8.0 + 1e-9,
            f"level {n}: diag={r['diag']:.6f} >= 2 w_n; lam_max={r['lam_max']:.6f} <= 8",
        )

    # Level-2 lock from the related audit.
    check(abs(rows[1]["nF"] - 1.165909) < 5e-6, "level 2 ||F|| lock 1.165909")

    # Kept Failure AB.1: no uniform 0.95 geometric collapse of rho.
    rhos = [r["rho"] for r in rows]
    rho_ratios = [rhos[i] / rhos[i - 1] for i in range(1, len(rhos))]
    print("rho_ratios", [f"{x:.8f}" for x in rho_ratios])
    check(rho_ratios[0] < 0.95, "level 1→2 rho ratio < 0.95 (early drop exists)")
    check(
        rho_ratios[3] > 0.99,  # level 4 → 5
        f"Kept Failure AB.1: rho ratio level 4→5 = {rho_ratios[3]:.8f} > 0.99 "
        "(rejects rho_{n+1}/rho_n ≤ 0.95 for every n)",
    )
    if through >= 6:
        check(
            rho_ratios[4] > 0.999,  # level 5 → 6
            f"level 5→6 rho ratio = {rho_ratios[4]:.8f} > 0.999",
        )
    # Share stays above 1 (combinatorial current exceeds fractional).
    check(
        all(r["sigma"] > 1.35 for r in rows[2:]),
        "sigma_n > 1.35 for every n≥3 on this range (not a limit claim)",
    )
    # AA floor still sits under live ||F_frac|| (sanity, not a re-proof of AA).
    floor_F = (9.0 / 4.0) / (
        4.0 ** (-ALPHA) * (8.0 + (4.0 / math.pi) * math.gamma(0.5 - ALPHA) * math.gamma(1.5) / math.gamma(2.0 - ALPHA))
    ) * (math.sqrt(21.0) / 7.0)
    # Use the published AA floor constant rather than recomputing Beta from scratch:
    floor_F_pub = 0.084015
    check(
        all(r["nF_frac"] >= floor_F_pub for r in rows),
        f"live ||F_frac|| stays above AA floor {floor_F_pub} (AB uses AA, does not re-prove it)",
    )

    print()
    print("Theorem AB: sigma = rho/w_n with w_min <= w_n <= 8^a;")
    print("for alpha<1/2, AA floors ||F_frac||, so ||F||->0 iff rho->0.")
    print("Kept Failure AB.1: no uniform 0.95 geometric collapse of rho.")
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
