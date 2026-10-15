#!/usr/bin/env python3
"""Uniform floor on the neighbor-degree share; no geometric gap-ratio blowup.

Witness for Theorem AF and Kept Failure AF.1 in GASKET.md.

Notation of Theorem AE: w_n = -(L^α)_{ap} = -(L^α)_{aq}, diag = (L^α)_{aa},
λ_n = 2 w_n / diag, γ_n = ḡ_far / ḡ_nbr, ρ_n = λ_n/(λ_n + (1-λ_n) γ_n).
w_min(α) = α 4^{α-1} (Theorem J; equality w_n = 1 at α = 1).

Theorem AF (every α ∈ (0,1], every level n ≥ 1):
  (1) floor            λ_n ≥ λ_min(α) := 2 w_min / 8^α = 2α 4^{α-1} / 8^α > 0
                       (diag ≤ λ_max(L)^α ≤ 8^α, as in AE);
  (1') Jensen floor    diag ≤ (L_aa)^α = 2^α (concavity of t^α against the
                       unit spectral measure of e_a), so
                       λ_n ≥ λ_J(α) := 2 w_min / 2^α = α 2^{α-1} ≥ λ_min,
                       and w_n ≤ 2^{α-1} (sharpens AE's 8^α/2);
  (2) for α ∈ (0, log3/log5): ρ_n ≥ 1/(1 + C_J γ_n), C_J = 1/λ_J - 1, hence
      ||F|| → 0  ⇒  γ_n → +∞   (one-sided; the converse needs limsup λ_n < 1,
      which is NOT proved here).

Kept Failure AF.1: the law "γ_{n+1}/γ_n ≥ 1.05 for every n, hence γ → ∞
geometrically and ||F(n,0.45)|| → 0 by AF + AE" is false: every live ratio
through level 6 is < 1.  This does not prove a positive limit.

Observation only (not a limit theorem): on levels 1-6 at α = 0.45, γ_n ∈ (1,2)
decreasing toward ≈1.375, λ_n decreasing toward ≈0.529.

Zero vs positive for combinatorial ||F(n,0.45)|| stays OPEN. Conjecture J.1
untouched. Not thrust.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_degree_share import level_row, spectral_power  # noqa: E402
from gasket_fractional_limit import build_gasket, combinatorial_laplacian  # noqa: E402

ALPHA = 0.45
GAMMA_C = math.log(3.0) / math.log(5.0)

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def w_min(alpha: float) -> float:
    return alpha * 4.0 ** (alpha - 1.0)


def lam_min(alpha: float) -> float:
    """Theorem AF (1): 2 w_min / 8^alpha."""
    return 2.0 * w_min(alpha) / 8.0 ** alpha


def lam_jensen(alpha: float) -> float:
    """Theorem AF (1'): 2 w_min / 2^alpha = alpha 2^(alpha-1)."""
    return 2.0 * w_min(alpha) / 2.0 ** alpha


def corner_diag_alpha1(level: int) -> tuple[float, float, float]:
    """At alpha=1: (diag, w_n, lambda) straight from L."""
    P, edges, corners = build_gasket(level)
    L = combinatorial_laplacian(len(P), edges)
    a = int(corners[0])
    nbrs = [j for j in range(len(P)) if j != a and L[a, j] != 0]
    assert len(nbrs) == 2
    diag = float(L[a, a])
    wn = float(-L[a, nbrs[0]])
    return diag, wn, 2.0 * wn / diag


def main(through: int = 6) -> int:
    assert through >= 4
    print("Theorem AF: uniform floor on the neighbor-degree share lambda_n")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0).")
    print("Not thrust. Does not decide zero vs positive for ||F||.\n")

    lm = lam_min(ALPHA)
    lj = lam_jensen(ALPHA)
    lm_closed = 2.0 * ALPHA * 4.0 ** (ALPHA - 1.0) / 8.0 ** ALPHA
    lm_pow2 = ALPHA * 2.0 ** (2.0 * ALPHA - 1.0 - 3.0 * ALPHA)  # = 0.45 * 2^(-1.45)
    print(f"lambda_min(0.45) = 2*0.45*4^(-0.55)/8^0.45 = {lm:.12f}")
    print(f"                 = 0.45 * 2^(-1.45)       = {lm_pow2:.12f}")
    print(f"lambda_J(0.45)   = 0.45 * 2^(-0.55)       = {lj:.12f}")
    C_min = 1.0 / lm - 1.0
    C_J = 1.0 / lj - 1.0
    print(f"C = 1/lambda_min - 1 = {C_min:.12f};  C_J = 1/lambda_J - 1 = {C_J:.12f}\n")
    check(abs(lm - lm_closed) < 1e-15 and abs(lm - lm_pow2) < 1e-15,
          "lambda_min closed forms agree (2a 4^(a-1)/8^a = a 2^(-1-a))")
    check(abs(lj - ALPHA * 2.0 ** (ALPHA - 1.0)) < 1e-15,
          "lambda_J = a 2^(a-1)")
    check(0.0 < lm <= lj < 1.0, "0 < lambda_min <= lambda_J < 1 at a=0.45")

    # Floors are positive and below 1 on a grid of alpha in (0,1].
    grid = np.linspace(0.01, 1.0, 100)
    check(all(0.0 < lam_min(a) <= lam_jensen(a) <= 1.0 + 1e-15 for a in grid),
          "0 < lambda_min(a) <= lambda_J(a) <= 1 on a grid of a in (0,1]")

    # Live table at alpha = 0.45.
    rows = []
    for n in range(1, through + 1):
        r = level_row(n, ALPHA)
        rows.append(r)
        print(
            f"level {n} N={r['N']} lam={r['lam']:.8f} gamma={r['gamma']:.8f} "
            f"rho={r['rho']:.8f} w_n={r['wn']:.8f} diag={r['diag']:.8f}",
            flush=True,
        )

    for n, r in enumerate(rows, start=1):
        check(r["diag"] <= 8.0 ** ALPHA + 1e-12,
              f"level {n}: diag={r['diag']:.8f} <= 8^a (AE)")
        check(r["diag"] <= 2.0 ** ALPHA + 1e-12,
              f"level {n}: diag={r['diag']:.8f} <= 2^a={2.0**ALPHA:.8f} (Jensen, AF 1')")
        check(r["wn"] >= w_min(ALPHA) - 1e-12 and r["wn"] <= 2.0 ** (ALPHA - 1.0) + 1e-12,
              f"level {n}: w_min <= w_n={r['wn']:.8f} <= 2^(a-1)={2.0**(ALPHA-1):.8f}")
        check(r["lam"] >= lj - 1e-12 and lj >= lm,
              f"level {n}: lambda={r['lam']:.8f} >= lambda_J={lj:.8f} >= lambda_min={lm:.8f}")
        rho_fac = r["lam"] / (r["lam"] + (1.0 - r["lam"]) * r["gamma"])
        check(abs(r["rho"] - rho_fac) < 1e-12,
              f"level {n}: AE factorization rho = lam/(lam+(1-lam)gam) (sanity)")
        check(r["Wfar"] > 0 and r["gamma"] > 0,
              f"level {n}: W_far > 0 and gamma > 0")
        check((1.0 - r["lam"]) / r["lam"] <= C_J + 1e-12,
              f"level {n}: (1-lam)/lam={(1-r['lam'])/r['lam']:.8f} <= C_J")
        check(r["rho"] >= 1.0 / (1.0 + C_J * r["gamma"]) - 1e-12,
              f"level {n}: rho={r['rho']:.8f} >= 1/(1+C_J gamma)="
              f"{1.0/(1.0+C_J*r['gamma']):.8f} (AF 2)")

    # Other alpha in (0, 1): floors hold live through level 5.
    for a in (0.10, 0.30, 0.60, 0.68, 0.80, 0.95):
        for n in range(1, 6):
            r = level_row(n, a)
            ok = (r["lam"] >= lam_jensen(a) - 1e-12 and lam_jensen(a) >= lam_min(a)
                  and r["diag"] <= 2.0 ** a + 1e-12)
            check(ok, f"a={a:.2f} level {n}: lambda={r['lam']:.6f} >= "
                      f"lambda_J={lam_jensen(a):.6f} >= lambda_min={lam_min(a):.6f}")
    # alpha = 1: lambda = 1 exactly (both neighbors are the whole degree).
    for n in range(1, 5):
        diag, wn, lam = corner_diag_alpha1(n)
        check(diag == 2.0 and wn == 1.0 and lam == 1.0 and lam >= lam_min(1.0),
              f"a=1 level {n}: diag=2, w_n=1=w_min(1), lambda=1 >= lambda_min(1)={lam_min(1.0):.4f}")
    # Spot check of the Jensen inequality on the full diagonal at alpha=0.45, level 3.
    P, edges, corners = build_gasket(3)
    L = combinatorial_laplacian(len(P), edges)
    A = spectral_power(L, ALPHA)
    check(bool(np.all(np.diag(A) <= np.diag(L) ** ALPHA + 1e-12)),
          "level 3: (L^a)_ii <= (L_ii)^a at every vertex (Jensen)")

    # Recorded AE table (unchanged).
    lam_rec = (0.646893, 0.553135, 0.533143, 0.529482, 0.528861, 0.528759)
    gam_rec = (1.443835, 1.390438, 1.378096, 1.375485, 1.374963, 1.374864)
    for i, r in enumerate(rows[:6]):
        check(abs(r["lam"] - lam_rec[i]) < 5e-7 and abs(r["gamma"] - gam_rec[i]) < 5e-7,
              f"level {i+1}: lambda, gamma match recorded AE table")

    # Kept Failure AF.1: no geometric gamma blowup.
    gams = [r["gamma"] for r in rows]
    gam_ratios = [gams[i] / gams[i - 1] for i in range(1, len(gams))]
    print("gamma_ratios", [f"{x:.8f}" for x in gam_ratios])
    check(gam_ratios[0] < 1.0,
          f"Kept Failure AF.1: gamma ratio level 1->2 = {gam_ratios[0]:.8f} < 1 < 1.05")
    check(gam_ratios[1] < 1.0 and gam_ratios[2] < 1.0,
          "Kept Failure AF.1: gamma ratios level 2->3 and 3->4 < 1")
    check(all(x < 1.0 for x in gam_ratios),
          "every gamma ratio through the witness range is < 1 (observation)")
    check(all(gam_ratios[i] < gam_ratios[i + 1] for i in range(len(gam_ratios) - 1)),
          "gamma ratios increase toward 1 from below (observation)")
    check(all(1.0 < g < 2.0 for g in gams),
          "gamma_n in (1,2) on this range (observation; not a limit theorem)")
    lams = [r["lam"] for r in rows]
    check(all(lams[i] > lams[i + 1] for i in range(len(lams) - 1)),
          "lambda_n decreasing on this range (observation; not a limit theorem)")

    print()
    print("Theorem AF: lambda_n >= lambda_min = 2a 4^(a-1)/8^a; Jensen: lambda_n >= a 2^(a-1);")
    print("for a < log3/log5, ||F||->0 implies gamma_n -> infinity (one-sided).")
    print("Kept Failure AF.1: no geometric gamma blowup (ratios < 1 through level 6).")
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
