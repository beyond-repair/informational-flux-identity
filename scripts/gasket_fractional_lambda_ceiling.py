#!/usr/bin/env python3
"""Uniform ceiling on the neighbor-degree share; two-sided vanishing split.

Witness for Theorem AG and Kept Failure AG.1 in GASKET.md.

Notation of Theorems AE–AF: w_n = -(L^α)_{ap} = -(L^α)_{aq},
diag = (L^α)_{aa}, λ_n = 2 w_n / diag, W_far = diag - 2 w_n,
γ_n = ḡ_far / ḡ_nbr, ρ_n = λ_n / (λ_n + (1-λ_n) γ_n).

Theorem AG (every α ∈ (0,1), every level n ≥ 1):
  (1) two-step heat-kernel floor
        W_far ≥ W_far_min(α) := 2 α (1-α) / 4^{2-α} > 0
      (three distance-2 vertices from the corner; M=4I-L path counts
      (M^2)_{a3}=(M^2)_{a5}=1, (M^2)_{a4}=2, level-independent);
  (2) with AF's Jensen diag ≤ 2^α,
        λ_n ≤ λ_ceil(α) := 1 - α(1-α)/2^{3-α} < 1
      (at α=0.45: λ_ceil ≈ 0.957738);
  (3) for α ∈ (0, log3/log5): combining with AF's λ_n ≥ λ_J(α) > 0,
        ||F|| → 0  ⟺  γ_n → +∞
      (two-sided; upgrades AF's one-sided implication).

Kept Failure AG.1: the law "λ_{n+1}/λ_n ≥ 1.01 for every n, hence λ → 1
and AF's converse stays blocked" is false: every live λ ratio through
level 6 is < 1 (λ decreasing). This does not prove a positive ||F|| limit.

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
from gasket_fractional_lambda_floor import lam_jensen, w_min  # noqa: E402
from gasket_fractional_limit import build_gasket, combinatorial_laplacian  # noqa: E402

ALPHA = 0.45
GAMMA_C = math.log(3.0) / math.log(5.0)

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def w_far_min(alpha: float) -> float:
    """Theorem AG (1): 2 α (1-α) / 4^{2-α}."""
    return 2.0 * alpha * (1.0 - alpha) / (4.0 ** (2.0 - alpha))


def lam_ceil(alpha: float) -> float:
    """Theorem AG (2): 1 - α(1-α)/2^{3-α}."""
    return 1.0 - alpha * (1.0 - alpha) / (2.0 ** (3.0 - alpha))


def single_path_bound(alpha: float) -> float:
    """α(1-α)/(2 * 4^{2-α}) from one length-2 walk."""
    return alpha * (1.0 - alpha) / (2.0 * 4.0 ** (2.0 - alpha))


def main(through: int = 6) -> int:
    assert through >= 4
    print("Theorem AG: uniform ceiling on the neighbor-degree share lambda_n")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0).")
    print("Not thrust. Does not decide zero vs positive for ||F||.\n")

    wf = w_far_min(ALPHA)
    lc = lam_ceil(ALPHA)
    lj = lam_jensen(ALPHA)
    sp = single_path_bound(ALPHA)
    # Closed forms agree: W_far_min / 2^α = α(1-α)/2^{3-α} = 1 - λ_ceil.
    gap = ALPHA * (1.0 - ALPHA) / (2.0 ** (3.0 - ALPHA))
    print(f"W_far_min(0.45) = 2*0.45*0.55/4^1.55 = {wf:.12f}")
    print(f"single-path     = 0.45*0.55/(2*4^1.55) = {sp:.12f}")
    print(f"lambda_ceil     = 1 - 0.45*0.55/2^2.55 = {lc:.12f}")
    print(f"lambda_J(0.45)  = {lj:.12f}")
    c_star = (1.0 - lc) / lc
    C_J = 1.0 / lj - 1.0
    print(f"c_* = (1-lambda_ceil)/lambda_ceil = {c_star:.12f}")
    print(f"C_J = 1/lambda_J - 1               = {C_J:.12f}\n")

    check(abs(wf - 4.0 * sp) < 1e-15, "W_far_min = 4 * single-path (1+2+1 walks)")
    check(abs(wf / (2.0 ** ALPHA) - gap) < 1e-15, "W_far_min/2^a = a(1-a)/2^(3-a)")
    check(abs(lc - (1.0 - gap)) < 1e-15, "lambda_ceil = 1 - a(1-a)/2^(3-a)")
    check(0.0 < lj < lc < 1.0, "0 < lambda_J < lambda_ceil < 1 at a=0.45")
    check(wf > 0.0 and gap > 0.0, "W_far_min > 0 and gap > 0")

    # Ceiling strictly below 1 on a grid of alpha in (0,1).
    grid = np.linspace(0.01, 0.99, 99)
    check(
        all(0.0 < lam_jensen(a) < lam_ceil(a) < 1.0 for a in grid),
        "0 < lambda_J(a) < lambda_ceil(a) < 1 on a grid of a in (0,1)",
    )
    # At alpha -> 1 the ceiling tends to 1 (escape hatch closes only for alpha < 1).
    check(abs(lam_ceil(1.0 - 1e-12) - 1.0) < 1e-4, "lambda_ceil(a)->1 as a->1-")

    # Structural: distance-2 path counts are level-independent.
    for n in (1, 2, 3, 4):
        P, edges, corners = build_gasket(n)
        N = len(P)
        L = combinatorial_laplacian(N, edges)
        M = 4.0 * np.eye(N) - L
        check(bool(np.all(M >= -1e-14)), f"level {n}: M=4I-L entrywise nonnegative")
        a = int(corners[0])
        nbrs = sorted(
            {j for (i, j) in edges if i == a} | {i for (i, j) in edges if j == a}
        )
        check(nbrs == [1, 2] or set(nbrs) == {1, 2}, f"level {n}: corner neighbors are vids 1,2")
        # Far vertices at graph distance 2 among the stable local patch {3,4,5}.
        M2 = M @ M
        check(
            abs(M2[a, 3] - 1.0) < 1e-14
            and abs(M2[a, 5] - 1.0) < 1e-14
            and abs(M2[a, 4] - 2.0) < 1e-14,
            f"level {n}: (M^2)_{{a3}}=1, (M^2)_{{a5}}=1, (M^2)_{{a4}}=2",
        )
        check(
            L[a, 3] == 0.0 and L[a, 4] == 0.0 and L[a, 5] == 0.0,
            f"level {n}: vids 3,4,5 are not combinatorial neighbors of the corner",
        )
        # Corner and its two neighbors have no edges leaving {0..5}.
        for v in (a, *nbrs):
            outs = [j for j in range(N) if j > 5 and L[v, j] != 0.0]
            check(outs == [], f"level {n}: vertex {v} has no edge outside local patch")

    # Live table at alpha = 0.45.
    rows = []
    for n in range(1, through + 1):
        r = level_row(n, ALPHA)
        rows.append(r)
        print(
            f"level {n} N={r['N']} lam={r['lam']:.8f} gamma={r['gamma']:.8f} "
            f"Wfar={r['Wfar']:.8f} diag={r['diag']:.8f} rho={r['rho']:.8f}",
            flush=True,
        )

    for n, r in enumerate(rows, start=1):
        check(
            r["Wfar"] >= wf - 1e-12,
            f"level {n}: W_far={r['Wfar']:.8f} >= W_far_min={wf:.8f}",
        )
        check(
            r["diag"] <= 2.0 ** ALPHA + 1e-12,
            f"level {n}: diag={r['diag']:.8f} <= 2^a (AF Jensen)",
        )
        check(
            r["lam"] <= lc + 1e-12,
            f"level {n}: lambda={r['lam']:.8f} <= lambda_ceil={lc:.8f}",
        )
        check(
            r["lam"] >= lj - 1e-12,
            f"level {n}: lambda={r['lam']:.8f} >= lambda_J={lj:.8f} (AF)",
        )
        # Two-sided sandwich on (1-λ)/λ.
        ratio = (1.0 - r["lam"]) / r["lam"]
        check(
            c_star - 1e-12 <= ratio <= C_J + 1e-12,
            f"level {n}: c_* <= (1-lam)/lam={ratio:.8f} <= C_J",
        )
        # Hence 1/(1+C_J γ) ≤ ρ ≤ 1/(1+c_* γ).
        lo = 1.0 / (1.0 + C_J * r["gamma"])
        hi = 1.0 / (1.0 + c_star * r["gamma"])
        check(
            lo - 1e-12 <= r["rho"] <= hi + 1e-12,
            f"level {n}: {lo:.8f} <= rho={r['rho']:.8f} <= {hi:.8f}",
        )
        rho_fac = r["lam"] / (r["lam"] + (1.0 - r["lam"]) * r["gamma"])
        check(abs(r["rho"] - rho_fac) < 1e-12, f"level {n}: AE factorization (sanity)")

    # Other alpha in (0,1): ceiling holds live through level 5.
    for a in (0.10, 0.30, 0.60, 0.68, 0.80, 0.95):
        for n in range(1, 6):
            r = level_row(n, a)
            ok = (
                r["Wfar"] >= w_far_min(a) - 1e-12
                and r["lam"] <= lam_ceil(a) + 1e-12
                and r["lam"] >= lam_jensen(a) - 1e-12
            )
            check(
                ok,
                f"a={a:.2f} level {n}: Wfar={r['Wfar']:.6f}>={w_far_min(a):.6f}; "
                f"lambda_J={lam_jensen(a):.6f}<=lambda={r['lam']:.6f}<="
                f"lambda_ceil={lam_ceil(a):.6f}",
            )

    # Spot-check off-diagonal lower bound vs heat-kernel integral at level 3.
    P, edges, corners = build_gasket(3)
    L = combinatorial_laplacian(len(P), edges)
    A = spectral_power(L, ALPHA)
    a = int(corners[0])
    for j, m2 in ((3, 1.0), (4, 2.0), (5, 1.0)):
        lower = m2 * sp  # m2 * single_path with unit (M^2)=1 bound uses sp; for m2=2 use 2*sp
        # single_path assumes (M^2)=1; for general m2 the bound is m2 * sp
        check(
            -float(A[a, j]) >= lower - 1e-12,
            f"level 3: -(L^a)_{{a{j}}}={-float(A[a,j]):.8f} >= {lower:.8f} (heat kernel)",
        )

    # Recorded AE/AF table (unchanged).
    lam_rec = (0.646893, 0.553135, 0.533143, 0.529482, 0.528861, 0.528759)
    gam_rec = (1.443835, 1.390438, 1.378096, 1.375485, 1.374963, 1.374864)
    for i, r in enumerate(rows[:6]):
        check(
            abs(r["lam"] - lam_rec[i]) < 5e-7 and abs(r["gamma"] - gam_rec[i]) < 5e-7,
            f"level {i+1}: lambda, gamma match recorded AE/AF table",
        )

    # Kept Failure AG.1: no geometric climb of lambda toward 1.
    lams = [r["lam"] for r in rows]
    lam_ratios = [lams[i] / lams[i - 1] for i in range(1, len(lams))]
    print("lambda_ratios", [f"{x:.8f}" for x in lam_ratios])
    check(
        lam_ratios[0] < 1.0,
        f"Kept Failure AG.1: lambda ratio level 1->2 = {lam_ratios[0]:.8f} < 1 < 1.01",
    )
    check(
        lam_ratios[1] < 1.0 and lam_ratios[2] < 1.0,
        "Kept Failure AG.1: lambda ratios level 2->3 and 3->4 < 1",
    )
    check(
        all(x < 1.0 for x in lam_ratios),
        "every lambda ratio through the witness range is < 1 (observation)",
    )
    check(
        all(lams[i] > lams[i + 1] for i in range(len(lams) - 1)),
        "lambda_n decreasing on this range (observation; not a limit theorem)",
    )
    check(
        all(r["lam"] < 0.70 for r in rows),
        "lambda_n < 0.70 on this range (observation; well below lambda_ceil)",
    )

    print()
    print("Theorem AG: W_far >= 2 a (1-a)/4^(2-a); lambda_n <= 1 - a(1-a)/2^(3-a);")
    print("for a < log3/log5, ||F||->0 iff gamma_n -> infinity (two-sided, with AF).")
    print("Kept Failure AG.1: no geometric lambda climb toward 1 (ratios < 1 through level 6).")
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
