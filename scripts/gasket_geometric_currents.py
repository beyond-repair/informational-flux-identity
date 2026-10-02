#!/usr/bin/env python3
"""Geometric 1/d^2 gasket corner currents: exact 12/5 growth, no finite limit.

Checks Theorem I and Kept Failure I.1 in GASKET.md. No thrust claim.
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


def geometric_laplacian(P: np.ndarray, edges: List[Edge]) -> np.ndarray:
    n = len(P)
    L = np.zeros((n, n), dtype=float)
    for i, j in edges:
        d = float(np.linalg.norm(P[i] - P[j]))
        w = 1.0 / max(d * d, 1e-18)
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
    return L


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


def audit_at_level(level: int, bc: np.ndarray):
    P, edges, corners = build_gasket(level)
    Lc = combinatorial_laplacian(len(P), edges)
    Lg = geometric_laplacian(P, edges)
    # All finest edges share length 2^{-level}, so Lg = 4^{level} Lc.
    lengths = {round(float(np.linalg.norm(P[i] - P[j])), 12) for i, j in edges}
    u = dirichlet_solve(Lg, corners, bc)
    # Same u solves the combinatorial Dirichlet problem.
    u_c = dirichlet_solve(Lc, corners, bc)
    Ig = (Lg @ u)[corners]
    Ic = (Lc @ u_c)[corners]
    return {
        "level": level,
        "N": len(P),
        "lengths": lengths,
        "scale_err": float(np.max(np.abs(Lg - (4**level) * Lc))),
        "u_err": float(np.max(np.abs(u - u_c))),
        "I_geom": Ig,
        "I_comb": Ic,
        "sum_geom": float(Ig.sum()),
        "nF_geom": corner_flux_norm(P, corners, Ig),
        "nF_comb": corner_flux_norm(P, corners, Ic),
    }


def main() -> None:
    bc = np.array([1.0, -0.5, 0.0])
    rows = [audit_at_level(level, bc) for level in (1, 2, 3, 4, 5)]

    for row in rows:
        n = row["level"]
        assert row["lengths"] == {round(2.0 ** (-n), 12)}, row["lengths"]
        assert row["scale_err"] < 1e-9, row["scale_err"]
        assert row["u_err"] < 1e-9, row["u_err"]
        assert abs(row["sum_geom"]) < 1e-8, row["sum_geom"]
        # I_geom = 4^n I_comb
        assert np.allclose(row["I_geom"], (4**n) * row["I_comb"], atol=1e-8)
        expect = (12 / 5) ** n * math.sqrt(21) / 2
        assert abs(row["nF_geom"] - expect) < 1e-8, (row["nF_geom"], expect)
        # Combinatorial control still matches Theorem G.
        expect_c = (3 / 5) ** n * math.sqrt(21) / 2
        assert abs(row["nF_comb"] - expect_c) < 1e-8

    ratios = [rows[i + 1]["nF_geom"] / rows[i]["nF_geom"] for i in range(4)]
    for r in ratios:
        assert abs(r - 12 / 5) < 1e-8, r

    # Kept Failure I.1: norms grow; they do not approach a finite limit on 1..5.
    assert rows[-1]["nF_geom"] > 10 * rows[0]["nF_geom"]
    assert all(r > 2.0 for r in ratios)

    print("bc", bc.tolist())
    for row in rows:
        print(
            f"level {row['level']} N={row['N']} "
            f"sum_geom={row['sum_geom']:.3e} "
            f"||F_geom||={row['nF_geom']:.8f} "
            f"expect={(12/5)**row['level'] * math.sqrt(21)/2:.8f}"
        )
    print("ratios_geom", [f"{r:.8f}" for r in ratios])
    print("exact_ratio", "12/5")
    print("finite_nonzero_limit", False)
    print("limit_is_infinity", True)


if __name__ == "__main__":
    main()
