#!/usr/bin/env python3
"""Fractional gasket corner currents: neutrality and failure of the 3/5 law.

Checks Theorem H and Kept Failure H.1 in GASKET.md. No thrust claim.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np

Edge = Tuple[int, int]


def build_gasket(level: int):
    c0 = np.array([0.0, 0.0])
    c1 = np.array([1.0, 0.0])
    c2 = np.array([0.5, math.sqrt(3.0) / 2.0])
    pos_map: Dict[Tuple[float, float], int] = {}
    positions: List[np.ndarray] = []
    edges: set[Edge] = set()

    def key(p: np.ndarray) -> Tuple[float, float]:
        return (round(float(p[0]), 10), round(float(p[1]), 10))

    def vid(p: np.ndarray) -> int:
        k = key(p)
        if k not in pos_map:
            pos_map[k] = len(positions)
            positions.append(p.copy())
        return pos_map[k]

    def rec(a, b, c, d):
        if d == 0:
            ia, ib, ic = vid(a), vid(b), vid(c)
            for i, j in ((ia, ib), (ib, ic), (ic, ia)):
                if i != j:
                    edges.add((min(i, j), max(i, j)))
            return
        rec(a, 0.5 * (a + b), 0.5 * (c + a), d - 1)
        rec(0.5 * (a + b), b, 0.5 * (b + c), d - 1)
        rec(0.5 * (c + a), 0.5 * (b + c), c, d - 1)

    rec(c0, c1, c2, level)
    P = np.vstack(positions)
    corners = np.array(
        [pos_map[key(c0)], pos_map[key(c1)], pos_map[key(c2)]], dtype=int
    )
    return P, sorted(edges), corners


def combinatorial_laplacian(n: int, edges: List[Edge]) -> np.ndarray:
    L = np.zeros((n, n), dtype=float)
    for i, j in edges:
        L[i, i] += 1.0
        L[j, j] += 1.0
        L[i, j] -= 1.0
        L[j, i] -= 1.0
    return L


def fractional_laplacian(L: np.ndarray, alpha: float) -> np.ndarray:
    w, V = np.linalg.eigh(L)
    w = np.clip(w, 0.0, None)
    return (V * (w**alpha)) @ V.T


def dirichlet_solve(A: np.ndarray, corners: np.ndarray, u_c: np.ndarray) -> np.ndarray:
    n = A.shape[0]
    mask = np.ones(n, dtype=bool)
    mask[corners] = False
    interior = np.where(mask)[0]
    u = np.zeros(n, dtype=float)
    u[corners] = u_c
    rhs = -A[np.ix_(interior, corners)] @ u_c
    u[interior] = np.linalg.solve(A[np.ix_(interior, interior)], rhs)
    return u


def corner_flux_norm(P: np.ndarray, corners: np.ndarray, currents: np.ndarray) -> float:
    centroid = P.mean(axis=0)
    F = np.zeros(2, dtype=float)
    for k, c in enumerate(corners):
        F += float(currents[k]) * (P[int(c)] - centroid)
    return float(np.linalg.norm(F))


def audit_at_level(level: int, alpha: float, bc: np.ndarray):
    P, edges, corners = build_gasket(level)
    L = combinatorial_laplacian(len(P), edges)
    A = fractional_laplacian(L, alpha)
    u = dirichlet_solve(A, corners, bc)
    Au = A @ u
    Lu = L @ u
    I_frac = Au[corners]
    I_comb = Lu[corners]
    mask = np.ones(len(P), dtype=bool)
    mask[corners] = False
    interior = np.where(mask)[0]
    return {
        "level": level,
        "N": len(P),
        "I_frac": I_frac,
        "I_comb": I_comb,
        "sum_frac": float(I_frac.sum()),
        "sum_comb": float(I_comb.sum()),
        "nF_frac": corner_flux_norm(P, corners, I_frac),
        "nF_comb": corner_flux_norm(P, corners, I_comb),
        "res_A": float(np.max(np.abs(Au[interior]))),
        "A_ones": float(np.max(np.abs(A @ np.ones(len(P))))),
    }


def main() -> None:
    bc = np.array([1.0, -0.5, 0.0])
    alpha = 0.45

    # Spectral assumption check: L^α 1 = 0 at a representative level.
    P2, E2, _ = build_gasket(2)
    A2 = fractional_laplacian(combinatorial_laplacian(len(P2), E2), alpha)
    assert np.max(np.abs(A2 @ np.ones(len(P2)))) < 1e-6

    rows = [audit_at_level(level, alpha, bc) for level in (1, 2, 3, 4)]

    # Theorem H: fractional currents neutral (within eigh noise).
    for row in rows:
        assert row["res_A"] < 1e-10, row
        assert abs(row["sum_frac"]) < 1e-4, row  # level-4 eigh noise ~1e-6
        assert abs(row["sum_comb"]) < 1e-4, row

    # Match the related-repo audit lock at level 2.
    assert abs(rows[1]["nF_comb"] - 1.165909) < 5e-3

    # Kept Failure H.1: ratios are not the constant 3/5.
    ratios_comb = [
        rows[i + 1]["nF_comb"] / rows[i]["nF_comb"] for i in range(3)
    ]
    ratios_frac = [
        rows[i + 1]["nF_frac"] / rows[i]["nF_frac"] for i in range(3)
    ]
    for r in ratios_comb + ratios_frac:
        assert abs(r - 0.6) > 0.15, r
    # And not even a single common geometric factor across levels.
    assert max(ratios_comb) - min(ratios_comb) > 0.05
    assert max(ratios_frac) - min(ratios_frac) > 0.05

    # α = 1 control: recovers Theorem G norm law.
    expect = [(3 / 5) ** n * math.sqrt(21) / 2 for n in (1, 2, 3, 4)]
    ctrl = [audit_at_level(level, 1.0, bc) for level in (1, 2, 3, 4)]
    for row, ex in zip(ctrl, expect):
        assert abs(row["nF_comb"] - ex) < 1e-9
        assert abs(row["nF_frac"] - ex) < 1e-9
        assert abs(row["sum_frac"]) < 1e-9
    for i in range(3):
        assert abs(ctrl[i + 1]["nF_comb"] / ctrl[i]["nF_comb"] - 0.6) < 1e-9

    print("alpha", alpha)
    print("bc", bc.tolist())
    for row in rows:
        print(
            f"level {row['level']} N={row['N']} "
            f"sum_frac={row['sum_frac']:.3e} sum_comb={row['sum_comb']:.3e} "
            f"||F_frac||={row['nF_frac']:.8f} ||F_comb||={row['nF_comb']:.8f}"
        )
    print("ratios_comb", [f"{r:.8f}" for r in ratios_comb])
    print("ratios_frac", [f"{r:.8f}" for r in ratios_frac])
    print("not_3/5", True)
    print("alpha1_control_ok", True)
    print("limit_nonzero", "OPEN")


if __name__ == "__main__":
    main()
