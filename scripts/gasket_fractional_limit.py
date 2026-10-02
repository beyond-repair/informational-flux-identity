#!/usr/bin/env python3
"""Audit-norm limit for spectral L^α on the combinatorial gasket.

Witness for Theorems J, K, L and Kept Failures J.1–J.3 in GASKET.md.
The norm is the audit norm: combinatorial corner currents of the L^α-Dirichlet
solution, the same net_flux as sierpinski-geometry-045. Not a thrust claim.

Levels 1–6 are the fast check. Level 7 is included because the kept-failure
ratios are not visible at level 4. Pass --through N to recompute a higher level
(level 8 is several tens of minutes).
"""

from __future__ import annotations

import math
import sys
from typing import Dict, List, Tuple

import numpy as np

Edge = Tuple[int, int]
BC = np.array([1.0, -0.5, 0.0])
ALPHA_C = math.log(3.0) / math.log(5.0)  # log(3)/log(5) ≈ 0.682606


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
    L = np.zeros((n, n), dtype=np.float64)
    for i, j in edges:
        L[i, i] += 1.0
        L[j, j] += 1.0
        L[i, j] -= 1.0
        L[j, i] -= 1.0
    return L


def corner_flux_norm(P: np.ndarray, corners: np.ndarray, currents: np.ndarray) -> float:
    centroid = P.mean(axis=0)
    F = np.zeros(2, dtype=np.float64)
    for k, c in enumerate(corners):
        F += float(currents[k]) * (P[int(c)] - centroid)
    return float(np.linalg.norm(F))


def decay_prefactor(alpha: float) -> Tuple[float, float, float]:
    """Return (C, sqrt(r), w_min) for Theorem K: ||F|| ≤ C * sqrt(r)^n."""
    w_min = alpha * (4.0 ** (alpha - 1.0))
    K = (7.0 / 2.0) ** alpha * (3.0 ** (1.0 - alpha))
    r = 3.0 * (5.0 ** (-alpha))
    C = (math.sqrt(21.0) / 5.0) * math.sqrt(2.0 * K / w_min)
    return C, math.sqrt(r), w_min


def audit_level(level: int, alphas: List[float]):
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    deg = np.diag(L)
    w, V = np.linalg.eigh(L)
    w = np.clip(w, 0.0, None)
    w[w < 1e-12] = 0.0
    mask = np.ones(N, dtype=bool)
    mask[corners] = False
    interior = np.where(mask)[0]
    # harmonic extension, for the comparison energy
    u_h = np.zeros(N, dtype=np.float64)
    u_h[corners] = BC
    u_h[interior] = np.linalg.solve(
        L[np.ix_(interior, interior)], -L[np.ix_(interior, corners)] @ BC
    )
    out = {
        "level": level,
        "N": N,
        "max_deg": float(deg.max()),
        "corner_deg": [float(deg[int(c)]) for c in corners],
        "by_alpha": {},
    }
    adj = [[] for _ in range(N)]
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)
    for alpha in alphas:
        A = (V * (w ** alpha)) @ V.T
        u = np.zeros(N, dtype=np.float64)
        u[corners] = BC
        u[interior] = np.linalg.solve(
            A[np.ix_(interior, interior)], -A[np.ix_(interior, corners)] @ BC
        )
        I = (L @ u)[corners]
        If = (A @ u)[corners]
        c0 = int(corners[0])
        nbr_w = [-float(A[c0, j]) for j in adj[c0]]
        Q = float(u @ (A @ u))
        Q_h = float(u_h @ (A @ u_h))
        Q1 = float(u_h @ (L @ u_h))
        out["by_alpha"][alpha] = {
            "nF": corner_flux_norm(P, corners, I),
            "nF_frac": corner_flux_norm(P, corners, If),
            "I": np.array(I, dtype=np.float64),
            "If": np.array(If, dtype=np.float64),
            "res": float(np.max(np.abs((A @ u)[interior]))),
            "umin": float(u.min()),
            "umax": float(u.max()),
            "nbr_w": nbr_w,
            "Q": Q,
            "Q_h": Q_h,
            "Q1": Q1,
            "ones": float(np.max(np.abs(A @ np.ones(N)))),
        }
    return out


def main(through: int = 7) -> None:
    assert through >= 6
    # α=0.45 is the open audit. α=0.9 sits above log(3)/log(5) so Theorem K applies.
    # α=0.10 is only reported through level 6 (same diagonalizations).
    fast_alphas = [0.10, 0.45, 0.90]
    rows = []
    for level in range(1, through + 1):
        alphas = fast_alphas if level <= 6 else [0.45]
        row = audit_level(level, alphas)
        rows.append(row)
        a45 = row["by_alpha"][0.45]
        print(
            f"level {level} N={row['N']} ||F||={a45['nF']:.10f} "
            f"||F_frac||={a45['nF_frac']:.10f} "
            f"I={np.array2string(a45['I'], precision=8, separator=',')} "
            f"sumI={float(a45['I'].sum()):.3e} res={a45['res']:.3e} "
            f"u=[{a45['umin']:.6f},{a45['umax']:.6f}]",
            flush=True,
        )

    # Degree lemma used by the weight bound.
    for row in rows:
        assert row["max_deg"] <= 4.0 + 1e-9, row["max_deg"]
        assert row["corner_deg"] == [2.0, 2.0, 2.0]

    # --- α = 0.45 audit sequence ---
    seq = [row["by_alpha"][0.45] for row in rows]
    # Level-2 lock from sierpinski-geometry-045 test_gasket_flux_audit.py.
    assert abs(seq[1]["nF"] - 1.165909) < 5e-6, seq[1]["nF"]
    # Neutrality and the max principle, through the computed range.
    cap = 3.0 * math.sqrt(21.0) / 5.0  # Theorem J
    for row, rec in zip(rows, seq):
        assert rec["res"] < 1e-8, rec["res"]
        assert abs(float(rec["I"].sum())) < 1e-6, rec["I"]
        assert abs(float(rec["If"].sum())) < 1e-6, rec["If"]
        assert rec["umin"] >= -0.5 - 1e-8
        assert rec["umax"] <= 1.0 + 1e-8
        assert rec["nF"] <= cap + 1e-9
        # Schur direction (1, -4/5, -1/5), up to roundoff.
        Ia = float(rec["I"][0])
        assert abs(float(rec["I"][1]) + 0.8 * Ia) < 1e-6 * max(1.0, abs(Ia))
        assert abs(float(rec["I"][2]) + 0.2 * Ia) < 1e-6 * max(1.0, abs(Ia))
        # Centroid identity: ||F|| = I_a * sqrt(21) / 5 when the triple is Schur's.
        assert abs(rec["nF"] - Ia * math.sqrt(21.0) / 5.0) < 1e-8
        # Neighbor weights beat the integral lower bound.
        w_min = 0.45 * (4.0 ** (0.45 - 1.0))
        assert min(rec["nbr_w"]) >= w_min - 1e-9, rec["nbr_w"]
        assert rec["ones"] < 1e-6

    ratios = [seq[i]["nF"] / seq[i - 1]["nF"] for i in range(1, len(seq))]
    print("ratios_0.45", [f"{r:.10f}" for r in ratios])
    # Strictly decreasing on the computed range, ratios increasing, all below 1.
    for r in ratios:
        assert r < 1.0
    for i in range(1, len(ratios)):
        assert ratios[i] > ratios[i - 1]

    # Kept Failure J.1: ratios do not stay ≤ 0.95.
    assert ratios[3] > 0.95  # level 4 → 5
    assert ratios[4] > 0.96  # level 5 → 6
    if through >= 7:
        assert ratios[5] > 0.97  # level 6 → 7

    # Kept Failure J.2: decrement contraction is not bounded by 5^{-0.45}.
    gaps = [1.0 - r for r in ratios]
    contractions = [gaps[i] / gaps[i - 1] for i in range(1, len(gaps))]
    print("decrement_ratios", [f"{q:.6f}" for q in contractions])
    five_alpha = 5.0 ** (-0.45)
    assert five_alpha < 0.49
    # By the step into level 6 the contraction has already left 5^{-0.45}.
    assert contractions[3] > 0.65  # gaps[4]/gaps[3], uses ||F|| through level 6
    if through >= 7:
        assert contractions[4] > 0.67

    # Theorem L: energy sandwich at every computed α=0.45 level.
    # ||F_frac|| ≥ w_min ||F|| and ||F|| ≤ κ √Q, with w_min = α 4^{α-1}.
    w_min_45 = 0.45 * (4.0 ** (0.45 - 1.0))
    kappa_45 = (math.sqrt(21.0) / 5.0) * math.sqrt(2.0 / w_min_45)
    for rec in seq:
        assert rec["nF_frac"] + 1e-9 >= w_min_45 * rec["nF"], (
            rec["nF_frac"],
            w_min_45 * rec["nF"],
        )
        assert rec["nF"] <= kappa_45 * math.sqrt(rec["Q"]) + 1e-9, (
            rec["nF"],
            kappa_45 * math.sqrt(rec["Q"]),
        )
        # Pairing Q = (7/5) I_a^α identifies energy with fractional currents.
        Ia_f = float(rec["If"][0])
        assert abs(rec["Q"] - (7.0 / 5.0) * Ia_f) < 1e-7 * max(1.0, abs(rec["Q"]))
        assert abs(rec["nF_frac"] - Ia_f * math.sqrt(21.0) / 5.0) < 1e-8

    # Kept Failure J.3: Q_{0.45} ratios do not stay ≤ 0.95.
    Qs = [rec["Q"] for rec in seq]
    q_ratios = [Qs[i] / Qs[i - 1] for i in range(1, len(Qs))]
    print("Q_0.45", [f"{q:.10f}" for q in Qs])
    print("Q_ratios_0.45", [f"{r:.10f}" for r in q_ratios])
    assert q_ratios[3] > 0.95  # level 4 → 5
    assert q_ratios[4] > 0.96  # level 5 → 6
    if through >= 7:
        assert q_ratios[5] > 0.97  # level 6 → 7 (witnessed when run that far)

    # Harmonic graph energy is the closed form; at α=0.45 the harmonic
    # comparison energy Q(u_H) is larger at level 6 than at level 2
    # (so it is not a vanishing upper bound on this range).
    for row in rows:
        n = row["level"]
        Q1 = row["by_alpha"][0.45]["Q1"]
        expect = (7.0 / 2.0) * (3.0 / 5.0) ** n
        assert abs(Q1 - expect) < 1e-8, (Q1, expect)
    assert rows[5]["by_alpha"][0.45]["Q_h"] > rows[1]["by_alpha"][0.45]["Q_h"]

    # --- Theorem K at α=0.9 > log(3)/log(5) ---
    assert 0.90 > ALPHA_C
    C, sqrt_r, w_min = decay_prefactor(0.90)
    assert sqrt_r < 1.0
    print(f"alpha_c={ALPHA_C:.6f} theorem_K α=0.9 C={C:.6f} sqrt_r={sqrt_r:.6f}")
    prev = None
    for row in rows:
        if 0.90 not in row["by_alpha"]:
            continue
        rec = row["by_alpha"][0.90]
        n = row["level"]
        K = (7.0 / 2.0) ** 0.90 * (3.0 ** (1.0 - 0.90))
        rate = 3.0 * (5.0 ** (-0.90))
        assert rec["Q_h"] <= K * (rate ** n) * (1.0 + 1e-8)
        assert min(rec["nbr_w"]) >= w_min - 1e-9
        bound = C * (sqrt_r ** n)
        assert rec["nF"] <= bound * (1.0 + 1e-8), (rec["nF"], bound)
        assert rec["umin"] >= -0.5 - 1e-8 and rec["umax"] <= 1.0 + 1e-8
        if prev is not None:
            assert rec["nF"] < prev
        prev = rec["nF"]
        print(
            f"  α=0.9 level {n} ||F||={rec['nF']:.8f} bound={bound:.8f} "
            f"Q_h={rec['Q_h']:.6f}",
            flush=True,
        )
    # The explicit bound itself decays, and by level 6 it is below the level-1 value.
    assert C * (sqrt_r ** 6) < rows[0]["by_alpha"][0.90]["nF"]

    # α=0.10, levels 1–6: still far above the α=0 value's neighborhood, recorded
    # only as data for the conjecture. Not a limit theorem.
    if 0.10 in rows[0]["by_alpha"]:
        small = [row["by_alpha"][0.10]["nF"] for row in rows if 0.10 in row["by_alpha"]]
        print("nF_alpha_0.10", [f"{x:.8f}" for x in small])
        for i in range(1, len(small)):
            assert small[i] < small[i - 1]
        # Drops after level 3 contract by less than 1/2 on this range.
        drops = [small[i - 1] - small[i] for i in range(1, len(small))]
        for i in range(2, len(drops)):
            assert drops[i] / drops[i - 1] < 0.5

    print("theorem_J_not_infinity", True)
    print("theorem_K_above_threshold_to_zero", True)
    print("theorem_L_energy_controls_audit_norm", True)
    print("kept_failure_J3_no_uniform_Q_ratio_0.95", True)
    print("alpha_0.45_zero_vs_positive", "OPEN")
    print("not_thrust", True)


if __name__ == "__main__":
    through = 7
    if len(sys.argv) > 1:
        through = int(sys.argv[1])
    main(through)
