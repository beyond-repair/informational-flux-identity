#!/usr/bin/env python3
"""Uniform energy floor for spectral L^α on the gasket when α < 1/2.

Witness for Theorem AA and Kept Failure AA.1 in GASKET.md.

Löwner–Heinz plus a boundary-path subgraph: the side a–b of the level-n
gasket is a path P_m with m = 2^n edges, its Laplacian M satisfies
0 ≤ M ≤ L, so M^α ≤ L^α for α ∈ (0,1). Cauchy–Schwarz against
v = e_a − e_b gives Q_α(u_n) ≥ (3/2)^2 / R_α(m), where
R_α(m) = v^T (P_m^α)^+ v has a closed-form cosine sum. For α < 1/2 that
sum is bounded uniformly in m by C(α) = 4^{-α}[8 + (4/π) B(1/2−α, 3/2)].
Hence Q_{0.45}(u_n) and the fractional-current audit ||F_frac|| = Q √21/7
stay bounded away from 0 at every level: the energy-collapse route of
Theorem L to ||F(n,0.45)|| → 0 is closed (Kept Failure AA.1).

It does NOT decide the combinatorial audit norm ||F(n,0.45)||. Not thrust.
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
SQRT21_7 = math.sqrt(21.0) / 7.0
GAP_AB = 1.5  # |u_a − u_b| for corner data (1, −1/2, 0)
# Recorded minimizer energies (GASKET.md, Kept Failures J.3 and J.4).
Q_RECORDED = {9: 0.948091, 10: 0.941537}

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def beta(x: float, y: float) -> float:
    return math.gamma(x) * math.gamma(y) / math.gamma(x + y)


def path_R_closed(m: int, alpha: float) -> float:
    """R_α(m) = (e_0−e_m)^T (P_m^α)^+ (e_0−e_m) by the cosine eigenbasis."""
    k = np.arange(1, m + 1, 2, dtype=np.float64)
    th = np.pi * k / (2.0 * (m + 1))
    return float(
        np.sum((4.0 * np.sin(th) ** 2) ** (-alpha) * (8.0 / (m + 1)) * np.cos(th) ** 2)
    )


def path_R_dense(m: int, alpha: float) -> float:
    n = m + 1
    P = np.zeros((n, n))
    for j in range(m):
        P[j, j] += 1.0
        P[j + 1, j + 1] += 1.0
        P[j, j + 1] -= 1.0
        P[j + 1, j] -= 1.0
    w, V = np.linalg.eigh(P)
    v = np.zeros(n)
    v[0], v[-1] = 1.0, -1.0
    c = V.T @ v
    keep = w > 1e-12
    return float(np.sum(c[keep] ** 2 * w[keep] ** (-alpha)))


def C_bound(alpha: float) -> float:
    return 4.0 ** (-alpha) * (8.0 + (4.0 / math.pi) * beta(0.5 - alpha, 1.5))


def R_limit(alpha: float) -> float:
    return 4.0 ** (-alpha) * (4.0 / math.pi) * beta(0.5 - alpha, 1.5)


def spectral_power(L: np.ndarray, alpha: float) -> np.ndarray:
    w, V = np.linalg.eigh(L)
    w = np.clip(w, 0.0, None)
    w[w < 1e-12] = 0.0
    return (V * (w ** alpha)) @ V.T


def level_record(level: int, alpha: float):
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    A = spectral_power(L, alpha)
    a, b = int(corners[0]), int(corners[1])
    side = sorted(
        (i for i in range(N) if abs(P[i, 1]) < 1e-12), key=lambda i: P[i, 0]
    )
    side_set = set(side)
    side_edges = [(i, j) for (i, j) in edges if i in side_set and j in side_set]
    consecutive = all(
        (min(side[t], side[t + 1]), max(side[t], side[t + 1])) in set(side_edges)
        for t in range(len(side) - 1)
    )
    M = combinatorial_laplacian(N, side_edges)
    Ma = spectral_power(M, alpha)
    lh_min = float(np.linalg.eigvalsh(A - Ma).min())
    interior = np.setdiff1d(np.arange(N), corners)
    u = np.zeros(N)
    u[corners] = BC
    u[interior] = np.linalg.solve(
        A[np.ix_(interior, interior)], -A[np.ix_(interior, corners)] @ BC
    )
    Q = float(u @ A @ u)
    Q_path = float(u @ Ma @ u)
    If = (A @ u)[corners]
    lam_max = float(np.linalg.eigvalsh(L).max())
    return {
        "N": N,
        "side_len": len(side) - 1,
        "side_edges": len(side_edges),
        "path_ok": consecutive and side[0] == a and side[-1] == b,
        "lh_min": lh_min,
        "Q": Q,
        "Q_path": Q_path,
        "nF_frac": corner_flux_norm(P, corners, If),
        "lam_max": lam_max,
    }


def main() -> int:
    print("Theorem AA: energy floor for L^alpha on the gasket, alpha < 1/2")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0).")
    print("Not thrust. Does not decide the combinatorial audit norm ||F||.\n")

    # 1. Closed-form path resistance vs dense pseudo-inverse.
    worst = 0.0
    for alpha in (0.10, 0.25, 0.45, 0.49, 0.9):
        for n in range(0, 8):
            m = 2 ** n
            worst = max(worst, abs(path_R_closed(m, alpha) - path_R_dense(m, alpha)))
    check(worst < 1e-9, f"closed-form R_alpha(m) matches dense (P^alpha)^+ for m=2^0..2^7 (max err {worst:.1e})")

    # 2. Uniform bound R_alpha(2^n) <= C(alpha) for alpha < 1/2, and approach to the limit.
    for alpha in (0.10, 0.25, 0.45, 0.49):
        C = C_bound(alpha)
        Rs = [path_R_closed(2 ** n, alpha) for n in range(0, 23)]
        check(max(Rs) <= C, f"alpha={alpha}: max_n<=22 R(2^n)={max(Rs):.6f} <= C(alpha)={C:.6f}")
        check(all(Rs[i] < Rs[i + 1] for i in range(len(Rs) - 1)) and Rs[-1] < R_limit(alpha),
              f"alpha={alpha}: R(2^n) increases toward the Riemann-sum limit {R_limit(alpha):.6f} (R(2^22)={Rs[-1]:.6f})")
    C45 = C_bound(ALPHA)
    floor_Q = GAP_AB ** 2 / C45
    floor_F = floor_Q * SQRT21_7
    ceil_Q = (7.0 / 6.0) * 8.0 ** ALPHA
    print(f"     alpha=0.45: C={C45:.6f}, Q floor={floor_Q:.6f}, ||F_frac|| floor={floor_F:.6f}, Q ceiling={ceil_Q:.6f}")

    # 3. Control: at alpha = 1/2 the path route fails (log growth), so the proof needs alpha < 1/2.
    R_half = [path_R_closed(2 ** n, 0.5) for n in (10, 14, 18, 22)]
    incs = [R_half[i + 1] - R_half[i] for i in range(3)]
    check(min(incs) > 3.0, f"control alpha=1/2: R(2^n) grows without bound, increments per 4 levels {[round(x, 4) for x in incs]}")

    # 4. Live gasket levels: Löwner–Heinz, path duality, ceiling, pairing.
    for n in range(1, 7):
        r = level_record(n, ALPHA)
        Rn = path_R_closed(2 ** n, ALPHA)
        check(r["path_ok"] and r["side_len"] == 2 ** n and r["side_edges"] == 2 ** n,
              f"level {n}: side a-b is a path of 2^{n} edges inside the gasket (N={r['N']})")
        check(r["lh_min"] > -1e-10, f"level {n}: lambda_min(L^a - M^a) = {r['lh_min']:.1e} >= 0 (Loewner-Heinz)")
        check(r["Q_path"] >= GAP_AB ** 2 / Rn - 1e-10 and r["Q"] >= r["Q_path"] - 1e-10,
              f"level {n}: Q={r['Q']:.6f} >= u^T M^a u={r['Q_path']:.6f} >= (3/2)^2/R(2^{n})={GAP_AB**2/Rn:.6f} >= floor {floor_Q:.6f}")
        check(r["Q"] <= ceil_Q and r["lam_max"] <= 8.0 + 1e-9,
              f"level {n}: Q <= (7/6)8^a={ceil_Q:.6f}; lambda_max(L)={r['lam_max']:.6f} <= 8")
        check(abs(r["nF_frac"] - r["Q"] * SQRT21_7) < 1e-10 and r["nF_frac"] >= floor_F,
              f"level {n}: ||F_frac||={r['nF_frac']:.6f} = Q*sqrt21/7 >= {floor_F:.6f}")
    for n, q in Q_RECORDED.items():
        check(floor_Q < q < ceil_Q, f"recorded level {n} energy Q={q} lies inside [{floor_Q:.6f}, {ceil_Q:.6f}]")

    # 5. Other alpha below 1/2 at level 4 (live).
    for alpha in (0.10, 0.25, 0.49):
        r = level_record(4, alpha)
        fl = GAP_AB ** 2 / C_bound(alpha)
        check(r["Q"] >= fl and r["lh_min"] > -1e-10,
              f"alpha={alpha}, level 4: Q={r['Q']:.6f} >= floor {fl:.6f}")

    print()
    print("Kept Failure AA.1: energy collapse Q_0.45(u_n) -> 0 (the Theorem L route to")
    print("||F(n,0.45)|| -> 0) is false at every level; ||F_frac|| -> 0 is false.")
    print("The combinatorial audit ||F(n,0.45)|| zero vs positive stays OPEN. Not thrust.")
    if failures:
        print(f"\n{len(failures)} check(s) failed")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
