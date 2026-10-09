#!/usr/bin/env python3
"""Multi-scale spectral L^α on the gasket: neutrality, no freeze inheritance.

Witness for Theorem X and Kept Failures X.1–X.2 in GASKET.md.

Vertices and multi-scale edges as in Assumption item 7 / Theorem W
(L_ms = Σ_g λ^g L^(g)). Spectral calculus as in item 4 applied to L_ms.
The L_ms^α-Dirichlet solution is not the combinatorial harmonic for α≠1,
so the λ=1 harmonic freeze at 5√21/4 does not freeze the α=0.45 audit,
and multi-scale fractional does not coincide with finest-only fractional.
Not thrust.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np

# Reuse the multi-scale builder and harmonic solver from Theorem W.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_multiscale_weights import (  # noqa: E402
    BC,
    RESIST,
    SQRT21_2,
    build_multiscale,
    closed_form,
    laplacian,
    solve_harmonic,
)

Edge = Tuple[int, int]
ALPHA = 0.45
FREEZE_L1 = 5.0 * math.sqrt(21.0) / 4.0  # Theorem W limit at λ=1
SCHUR = np.array([1.0, -0.8, -0.2])
# Recorded L_ms-audit norms of L_ms^{0.45} at λ=1 (levels 1..6)
F_MS_L = (
    3.749327948061,
    5.059027321162,
    6.363923868385,
    7.716285305812,
    9.127122921290,
    10.587687763937,
)
# Recorded fractional-current norms (L_ms^α currents)
F_MS_A = (
    1.554895275953,
    1.754213131387,
    1.924314598501,
    2.087851104087,
    2.250027663254,
    2.409734082231,
)
# Finest-edge combinatorial fractional L-audit (Theorem J / GASKET table)
F_FIN = (
    1.439599390,
    1.165908607,
    1.032271686,
    0.956205578,
    0.909692382,
    0.880108995,
)


def corner_flux_norm(P: np.ndarray, corners: np.ndarray, currents: np.ndarray) -> float:
    centroid = P.mean(axis=0)
    F = np.zeros(2, dtype=np.float64)
    for k, c in enumerate(corners):
        F += float(currents[k]) * (P[int(c)] - centroid)
    return float(np.linalg.norm(F))


def schur_err(I: np.ndarray) -> float:
    if abs(I[0]) < 1e-15:
        return float("inf")
    return float(np.max(np.abs(I / I[0] - SCHUR)))


def solve_fractional(L: np.ndarray, corners: np.ndarray, alpha: float):
    n = L.shape[0]
    w, V = np.linalg.eigh(L)
    w = np.clip(w, 0.0, None)
    w[w < 1e-12] = 0.0
    A = (V * (w ** alpha)) @ V.T
    mask = np.ones(n, dtype=bool)
    mask[corners] = False
    interior = np.where(mask)[0]
    u = np.zeros(n, dtype=np.float64)
    u[corners] = BC
    if len(interior):
        u[interior] = np.linalg.solve(
            A[np.ix_(interior, interior)],
            -A[np.ix_(interior, corners)] @ BC,
        )
    I_alpha = (A @ u)[corners]
    I_L = (L @ u)[corners]
    resid = (
        float(np.max(np.abs((A @ u)[interior]))) if len(interior) else 0.0
    )
    Q = float(u @ (A @ u))
    return {
        "u": u,
        "I_alpha": np.asarray(I_alpha, dtype=np.float64),
        "I_L": np.asarray(I_L, dtype=np.float64),
        "Q": Q,
        "resid": resid,
        "A": A,
        "w": w,
    }


def audit_ms(level: int, lam: float, alpha: float):
    P, edge_gen, finest, corners = build_multiscale(level)
    n = len(P)
    ew_ms = {e: lam ** g for e, g in edge_gen.items()}
    ew_fin = {e: 1.0 for e in finest}
    Lms = laplacian(n, ew_ms)
    Lf = laplacian(n, ew_fin)
    ms = solve_fractional(Lms, corners, alpha)
    fin = solve_fractional(Lf, corners, alpha)
    uh, Ih = solve_harmonic(Lms, corners)
    # Residual of L_ms^α on the combinatorial harmonic (nested-harmonicity test)
    A = ms["A"]
    mask = np.ones(n, dtype=bool)
    mask[corners] = False
    nest_resid = float(np.max(np.abs((A @ uh)[mask]))) if mask.any() else 0.0
    return {
        "level": level,
        "N": n,
        "nF_L": corner_flux_norm(P, corners, ms["I_L"]),
        "nF_alpha": corner_flux_norm(P, corners, ms["I_alpha"]),
        "nF_fin": corner_flux_norm(P, corners, fin["I_L"]),
        "nF_h": corner_flux_norm(P, corners, Ih),
        "I_L": ms["I_L"],
        "I_alpha": ms["I_alpha"],
        "sum_L": float(ms["I_L"].sum()),
        "sum_alpha": float(ms["I_alpha"].sum()),
        "Q": ms["Q"],
        "resid": ms["resid"],
        "nest_resid": nest_resid,
        "u_uh": float(np.max(np.abs(ms["u"] - uh))),
        "u_uf": float(np.max(np.abs(ms["u"] - fin["u"]))),
        "pred_W": closed_form(lam, level),
    }


def main() -> int:
    fail = 0

    def check(ok: bool, msg: str) -> None:
        nonlocal fail
        print(("ok   " if ok else "FAIL ") + msg)
        if not ok:
            fail += 1

    print("=== Theorem X: α=1 recovers Theorem W (multi-scale harmonic) ===")
    for lam in (1.0, RESIST, 2.0):
        for level in range(1, 6):
            row = audit_ms(level, lam, 1.0)
            check(
                abs(row["nF_L"] - row["pred_W"]) < 1e-9 * max(1.0, row["pred_W"]),
                f"λ={lam} n={level}: α=1 ||F_L||={row['nF_L']:.12f} = W-closed-form "
                f"(err {abs(row['nF_L']-row['pred_W']):.2e})",
            )
            check(
                abs(row["nF_L"] - row["nF_alpha"]) < 1e-9,
                f"λ={lam} n={level}: α=1 ⇒ L-currents = L^α-currents",
            )
            check(
                row["u_uh"] < 1e-9,
                f"λ={lam} n={level}: α=1 ⇒ u = combinatorial harmonic "
                f"(||u-uh||_∞={row['u_uh']:.2e})",
            )
            check(
                row["nest_resid"] < 1e-9,
                f"λ={lam} n={level}: α=1 nested residual {row['nest_resid']:.2e}",
            )

    print("\n=== Theorem X: neutrality + Schur line for L_ms^{0.45} (λ=1) ===")
    rows = []
    for level in range(1, 7):
        row = audit_ms(level, 1.0, ALPHA)
        rows.append(row)
        check(abs(row["sum_L"]) < 1e-8, f"n={level}: Σ I_L = {row['sum_L']:.2e}")
        check(
            abs(row["sum_alpha"]) < 1e-8,
            f"n={level}: Σ I_α = {row['sum_alpha']:.2e}",
        )
        check(
            schur_err(row["I_L"]) < 1e-8,
            f"n={level}: I_L on Schur line (err {schur_err(row['I_L']):.2e})",
        )
        check(
            schur_err(row["I_alpha"]) < 1e-8,
            f"n={level}: I_α on Schur line (err {schur_err(row['I_alpha']):.2e})",
        )
        check(
            row["resid"] < 1e-10,
            f"n={level}: Dirichlet residual {row['resid']:.2e}",
        )
        check(
            abs(row["nF_L"] - F_MS_L[level - 1]) < 5e-9,
            f"n={level}: ||F_L||={row['nF_L']:.12f} matches record "
            f"{F_MS_L[level-1]:.12f}",
        )
        check(
            abs(row["nF_alpha"] - F_MS_A[level - 1]) < 5e-9,
            f"n={level}: ||F_α||={row['nF_alpha']:.12f} matches record "
            f"{F_MS_A[level-1]:.12f}",
        )
        print(
            f"  n={level} N={row['N']}: ||F_L||={row['nF_L']:.10f} "
            f"||F_α||={row['nF_alpha']:.10f} Q={row['Q']:.8f} "
            f"||u-uh||={row['u_uh']:.3e} nest={row['nest_resid']:.3e}"
        )

    print("\n=== Theorem X: nested harmonicity fails for α=0.45 ===")
    for row in rows:
        check(
            row["nest_resid"] > 0.05,
            f"n={row['level']}: |(L_ms^{{0.45}} uh)_int|_∞ = "
            f"{row['nest_resid']:.6f} ≫ 0 (not harmonic for L_ms^α)",
        )
        check(
            row["u_uh"] > 0.05,
            f"n={row['level']}: ||u_ms,α - uh||_∞ = {row['u_uh']:.6f} ≫ 0",
        )

    print("\n=== Kept Failure X.1: λ=1 freeze does not transfer to α=0.45 ===")
    vals = [r["nF_L"] for r in rows]
    vals_a = [r["nF_alpha"] for r in rows]
    qs = [r["Q"] for r in rows]
    check(
        all(vals[i] < vals[i + 1] for i in range(5)),
        f"||F_L|| strictly increasing through n=6: {[f'{v:.6f}' for v in vals]}",
    )
    check(
        all(vals_a[i] < vals_a[i + 1] for i in range(5)),
        f"||F_α|| strictly increasing through n=6",
    )
    check(
        all(qs[i] < qs[i + 1] for i in range(5)),
        f"Q_0.45 strictly increasing through n=6",
    )
    check(
        vals[2] > FREEZE_L1,
        f"n=3 ||F_L||={vals[2]:.6f} already exceeds λ=1 harmonic freeze "
        f"5√21/4={FREEZE_L1:.6f}",
    )
    check(
        vals[-1] > FREEZE_L1 + 1.0,
        f"n=6 ||F_L||={vals[-1]:.6f} well past freeze (not approaching it)",
    )
    # Harmonic at same levels approaches freeze from below
    for level in (1, 3, 5):
        row_h = audit_ms(level, 1.0, 1.0)
        check(
            row_h["nF_L"] < FREEZE_L1,
            f"contrast: harmonic λ=1 n={level} ||F||={row_h['nF_L']:.6f} "
            f"< freeze {FREEZE_L1:.6f}",
        )

    print("\n=== Kept Failure X.2: multi-scale fractional ≠ finest-only ===")
    for i, row in enumerate(rows[:5]):  # levels 1..5 vs recorded finest
        check(
            abs(row["nF_fin"] - F_FIN[i]) < 5e-6,
            f"n={row['level']}: finest ||F||={row['nF_fin']:.9f} matches "
            f"Theorem J record {F_FIN[i]:.9f}",
        )
    check(
        all(rows[i]["nF_fin"] > rows[i + 1]["nF_fin"] for i in range(4)),
        "finest-edge ||F|| decreases through n=1..5 (Theorem J sequence)",
    )
    check(
        all(rows[i]["nF_L"] < rows[i + 1]["nF_L"] for i in range(4)),
        "multi-scale ||F_L|| increases through n=1..5 (opposite monotonicity)",
    )
    ratios = [rows[i]["nF_L"] / rows[i]["nF_fin"] for i in range(5)]
    check(
        all(ratios[i] < ratios[i + 1] for i in range(4)),
        f"||F_msL||/||F_fin|| strictly increasing: "
        f"{[f'{r:.4f}' for r in ratios]}",
    )
    check(
        ratios[-1] > 5.0,
        f"ratio at n=5 is {ratios[-1]:.4f} ≫ 1 (no level-independent constant)",
    )
    check(
        all(rows[i]["u_uf"] > 0.01 for i in range(5)),
        "pointwise ||u_ms - u_fin||_∞ > 0.01 at every level 1..5",
    )

    print("\n=== Bonus: multi-scale resistance fractional also diverges ===")
    resist_vals = []
    for level in range(1, 6):
        row = audit_ms(level, RESIST, ALPHA)
        resist_vals.append(row["nF_L"])
        check(abs(row["sum_L"]) < 1e-7, f"λ=5/3 n={level}: neutral")
        check(
            schur_err(row["I_L"]) < 1e-7,
            f"λ=5/3 n={level}: Schur",
        )
    check(
        all(resist_vals[i] < resist_vals[i + 1] for i in range(4)),
        f"λ=5/3 ||F_L|| increasing: {[f'{v:.4f}' for v in resist_vals]}",
    )
    check(
        resist_vals[-1] > 30.0,
        f"λ=5/3 n=5 ||F_L||={resist_vals[-1]:.4f} diverges (no freeze)",
    )

    print()
    if fail:
        print(f"{fail} check(s) failed")
        return 1
    print("all Theorem X / Kept Failure X.1–X.2 checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
