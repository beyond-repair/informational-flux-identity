#!/usr/bin/env python3
"""Uniform multiplicative edge weights on the combinatorial gasket.

Witness for Theorem S and Kept Failure S.1 in GASKET.md.
On the level-n finest-edge build, set w_ij = λ^n. Harmonic audit norms are
exactly (λ·3/5)^n √21/2 for corner data (1, -1/2, 0). Resistance λ = 5/3
gives the constant √21/2; combinatorial λ = 1 and geometric λ = 4 recover
Theorems G and I. Not thrust.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np

Edge = Tuple[int, int]
BC = np.array([1.0, -0.5, 0.0])
SQRT21_2 = math.sqrt(21.0) / 2.0
RESIST = 5.0 / 3.0
ALPHA = 0.45

# Recorded combinatorial fractional audit norms (α=0.45) from GASKET.md /
# gasket_fractional_limit.py, levels 1..7. Used only for Kept Failure S.1.
F_FRAC = (
    1.439599390,
    1.165908607,
    1.032271686,
    0.956205578,
    0.909692382,
    0.880108995,
    0.860817772,
)


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


def weighted_laplacian(n: int, edges: List[Edge], weight: float) -> np.ndarray:
    L = np.zeros((n, n), dtype=np.float64)
    for i, j in edges:
        L[i, i] += weight
        L[j, j] += weight
        L[i, j] -= weight
        L[j, i] -= weight
    return L


def corner_flux_norm(P: np.ndarray, corners: np.ndarray, currents: np.ndarray) -> float:
    centroid = P.mean(axis=0)
    F = np.zeros(2, dtype=np.float64)
    for k, c in enumerate(corners):
        F += float(currents[k]) * (P[int(c)] - centroid)
    return float(np.linalg.norm(F))


def harmonic_audit(level: int, lam: float):
    P, edges, corners = build_gasket(level)
    n = len(P)
    w = lam**level
    L = weighted_laplacian(n, edges, w)
    Lc = weighted_laplacian(n, edges, 1.0)
    mask = np.ones(n, dtype=bool)
    mask[corners] = False
    interior = np.where(mask)[0]
    u = np.zeros(n, dtype=np.float64)
    u[corners] = BC
    u[interior] = np.linalg.solve(
        L[np.ix_(interior, interior)], -L[np.ix_(interior, corners)] @ BC
    )
    u_c = np.zeros(n, dtype=np.float64)
    u_c[corners] = BC
    u_c[interior] = np.linalg.solve(
        Lc[np.ix_(interior, interior)], -Lc[np.ix_(interior, corners)] @ BC
    )
    I = (L @ u)[corners]
    Ic = (Lc @ u_c)[corners]
    return {
        "level": level,
        "N": n,
        "nF": corner_flux_norm(P, corners, I),
        "nF_comb": corner_flux_norm(P, corners, Ic),
        "I": np.array(I, dtype=np.float64),
        "sumI": float(I.sum()),
        "u_err": float(np.max(np.abs(u - u_c))),
        "scale_err": float(np.max(np.abs(L - w * Lc))),
        "expect": (lam * 3.0 / 5.0) ** level * SQRT21_2,
    }


def fractional_resistance_live(level: int):
    """L_res^α Dirichlet; audit from L_res currents equals (5/3)^n ||F_comb||."""
    P, edges, corners = build_gasket(level)
    n = len(P)
    rho = RESIST**level
    Lc = weighted_laplacian(n, edges, 1.0)
    w, V = np.linalg.eigh(Lc)
    w = np.clip(w, 0.0, None)
    w[w < 1e-12] = 0.0
    # L_res = rho Lc ⇒ L_res^α = rho^α Lc^α
    A_comb = (V * (w**ALPHA)) @ V.T
    A_res = (rho**ALPHA) * A_comb
    mask = np.ones(n, dtype=bool)
    mask[corners] = False
    interior = np.where(mask)[0]
    u = np.zeros(n, dtype=np.float64)
    u[corners] = BC
    u[interior] = np.linalg.solve(
        A_res[np.ix_(interior, interior)], -A_res[np.ix_(interior, corners)] @ BC
    )
    u2 = np.zeros(n, dtype=np.float64)
    u2[corners] = BC
    u2[interior] = np.linalg.solve(
        A_comb[np.ix_(interior, interior)], -A_comb[np.ix_(interior, corners)] @ BC
    )
    L_res = rho * Lc
    I_res = (L_res @ u)[corners]
    I_comb = (Lc @ u2)[corners]
    nF_res = corner_flux_norm(P, corners, I_res)
    nF_comb = corner_flux_norm(P, corners, I_comb)
    return {
        "u_err": float(np.max(np.abs(u - u2))),
        "nF_res": nF_res,
        "nF_comb": nF_comb,
        "scaled": rho * nF_comb,
        "res": float(np.max(np.abs((A_res @ u)[interior]))),
        "sumI": float(I_res.sum()),
    }


def main() -> None:
    # --- Theorem S: closed form for several λ, levels 1..5 ---
    lambdas = {
        "combinatorial": 1.0,
        "resistance": RESIST,
        "geometric": 4.0,
        "below": 1.5,  # 3/2 < 5/3 ⇒ decay rate 9/10
        "above": 2.0,  # 2 > 5/3 ⇒ growth rate 6/5
    }
    rows: Dict[str, list] = {name: [] for name in lambdas}
    for name, lam in lambdas.items():
        prev = None
        for level in range(1, 6):
            row = harmonic_audit(level, lam)
            rows[name].append(row)
            assert row["scale_err"] < 1e-9, row["scale_err"]
            assert row["u_err"] < 1e-9, row["u_err"]
            assert abs(row["sumI"]) < 1e-8, row["sumI"]
            assert abs(row["nF"] - row["expect"]) < 1e-8, (row["nF"], row["expect"])
            # Currents are λ^n times combinatorial.
            assert abs(row["nF"] - (lam**level) * row["nF_comb"]) < 1e-8
            if prev is not None:
                ratio = row["nF"] / prev
                assert abs(ratio - lam * 3.0 / 5.0) < 1e-8, ratio
            prev = row["nF"]
            print(
                f"{name} λ={lam} level {level} ||F||={row['nF']:.12f} "
                f"expect={row['expect']:.12f}",
                flush=True,
            )

    # Resistance is exactly constant √21/2.
    resist_vals = [r["nF"] for r in rows["resistance"]]
    for v in resist_vals:
        assert abs(v - SQRT21_2) < 1e-8, v
    assert abs(resist_vals[0] - resist_vals[-1]) < 1e-8

    # Threshold trichotomy on levels 1..5.
    below = [r["nF"] for r in rows["below"]]
    above = [r["nF"] for r in rows["above"]]
    assert all(below[i] > below[i + 1] for i in range(4))
    assert all(above[i] < above[i + 1] for i in range(4))
    geom = [r["nF"] for r in rows["geometric"]]
    assert geom[-1] > 10 * geom[0]
    comb = [r["nF"] for r in rows["combinatorial"]]
    assert comb[-1] < 0.2 * comb[0]

    # Recover Theorems G and I as special cases.
    assert abs(comb[1] - (3.0 / 5.0) ** 2 * SQRT21_2) < 1e-8
    assert abs(geom[0] - (12.0 / 5.0) * SQRT21_2) < 1e-8

    # --- Kept Failure S.1: resistance does not freeze α=0.45 ---
    # Live check through level 4: L_res^α shares the combinatorial fractional
    # Dirichlet solution, so ||F_res|| = (5/3)^n ||F_comb||.
    live_scaled = []
    for level in range(1, 5):
        rec = fractional_resistance_live(level)
        assert rec["u_err"] < 1e-8, rec["u_err"]
        assert rec["res"] < 1e-8, rec["res"]
        assert abs(rec["sumI"]) < 1e-6, rec["sumI"]
        assert abs(rec["nF_res"] - rec["scaled"]) < 1e-8
        assert abs(rec["nF_comb"] - F_FRAC[level - 1]) < 5e-6
        live_scaled.append(rec["nF_res"])
        print(
            f"frac_resist level {level} ||F_res||={rec['nF_res']:.12f} "
            f"(5/3)^n||F||={rec['scaled']:.12f}",
            flush=True,
        )
    # Not constant (contrast with harmonic resistance).
    assert live_scaled[0] != live_scaled[-1]
    assert abs(live_scaled[0] - SQRT21_2) > 0.05
    assert all(live_scaled[i] < live_scaled[i + 1] for i in range(3))

    # Recorded combinatorial norms through level 7: resistance-scaled sequence
    # is strictly increasing, so the harmonic cancellation identity fails.
    scaled = [(RESIST**n) * F_FRAC[n - 1] for n in range(1, 8)]
    for i in range(6):
        assert scaled[i] < scaled[i + 1], (scaled[i], scaled[i + 1])
    print("resistance_scaled_0.45", [f"{s:.8f}" for s in scaled])
    print("harmonic_resistance_constant", f"{SQRT21_2:.12f}")
    print("theorem_S_uniform_weights_closed_form", True)
    print("kept_failure_S1_no_fractional_resistance_freeze", True)
    print("alpha_0.45_combinatorial_limit_still_OPEN", True)
    print("not_thrust", True)


if __name__ == "__main__":
    main()
