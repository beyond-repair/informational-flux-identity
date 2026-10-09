#!/usr/bin/env python3
"""Multi-scale edge weights on the gasket: harmonic audit trichotomy.

Witness for Theorem W and Kept Failure W.1 in GASKET.md.

Vertices are the usual level-n gasket points. Edges are retained at every
generation 0..n (cell edges of every scale), with weight w = λ^{gen(e)}.
The combinatorial harmonic extension remains harmonic for every λ, and the
audit norm is exactly (√21/2) Σ_{k=0}^n (λ·3/5)^k. Resistance λ=5/3 makes
||F|| grow like (n+1), so multi-scale resistance does not freeze a finite
nonzero harmonic audit. Not thrust.
"""

from __future__ import annotations

import math
from typing import Dict, List, Tuple

import numpy as np

Edge = Tuple[int, int]
BC = np.array([1.0, -0.5, 0.0])
SQRT21_2 = math.sqrt(21.0) / 2.0
RESIST = 5.0 / 3.0


def build_multiscale(level: int):
    """Level-n vertices; edges of every generation 0..n with their gen index."""
    c0 = np.array([0.0, 0.0])
    c1 = np.array([1.0, 0.0])
    c2 = np.array([0.5, math.sqrt(3.0) / 2.0])
    pos_map: Dict[Tuple[float, float], int] = {}
    positions: List[np.ndarray] = []
    edge_gen: Dict[Edge, int] = {}
    finest: set[Edge] = set()

    def key(p: np.ndarray) -> Tuple[float, float]:
        return (round(float(p[0]), 10), round(float(p[1]), 10))

    def vid(p: np.ndarray) -> int:
        k = key(p)
        if k not in pos_map:
            pos_map[k] = len(positions)
            positions.append(p.copy())
        return pos_map[k]

    def add_edge(a: np.ndarray, b: np.ndarray, gen: int) -> None:
        ia, ib = vid(a), vid(b)
        if ia == ib:
            return
        e = (min(ia, ib), max(ia, ib))
        if e not in edge_gen or gen < edge_gen[e]:
            edge_gen[e] = gen

    def rec(a, b, c, d, gen):
        add_edge(a, b, gen)
        add_edge(b, c, gen)
        add_edge(c, a, gen)
        if d == 0:
            for x, y in ((a, b), (b, c), (c, a)):
                ia, ib = vid(x), vid(y)
                if ia != ib:
                    finest.add((min(ia, ib), max(ia, ib)))
            return
        ab, bc, ca = 0.5 * (a + b), 0.5 * (b + c), 0.5 * (c + a)
        rec(a, ab, ca, d - 1, gen + 1)
        rec(ab, b, bc, d - 1, gen + 1)
        rec(ca, bc, c, d - 1, gen + 1)

    rec(c0, c1, c2, level, 0)
    P = np.vstack(positions)
    corners = np.array(
        [pos_map[key(c0)], pos_map[key(c1)], pos_map[key(c2)]], dtype=int
    )
    return P, edge_gen, sorted(finest), corners


def laplacian(n: int, weights: Dict[Edge, float]) -> np.ndarray:
    L = np.zeros((n, n), dtype=np.float64)
    for (i, j), w in weights.items():
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
    return L


def solve_harmonic(L: np.ndarray, corners: np.ndarray):
    n = L.shape[0]
    mask = np.ones(n, dtype=bool)
    mask[corners] = False
    interior = np.where(mask)[0]
    u = np.zeros(n, dtype=np.float64)
    u[corners] = BC
    if len(interior):
        u[interior] = np.linalg.solve(
            L[np.ix_(interior, interior)], -L[np.ix_(interior, corners)] @ BC
        )
    I = (L @ u)[corners]
    return u, np.asarray(I, dtype=np.float64)


def corner_flux_norm(P: np.ndarray, corners: np.ndarray, currents: np.ndarray) -> float:
    centroid = P.mean(axis=0)
    F = np.zeros(2, dtype=np.float64)
    for k, c in enumerate(corners):
        F += float(currents[k]) * (P[int(c)] - centroid)
    return float(np.linalg.norm(F))


def closed_form(lam: float, n: int) -> float:
    r = lam * 3.0 / 5.0
    return SQRT21_2 * sum(r ** k for k in range(n + 1))


def limit_if_convergent(lam: float):
    if lam >= RESIST - 1e-15:
        return None
    return SQRT21_2 / (1.0 - lam * 3.0 / 5.0)


def audit_multiscale(level: int, lam: float):
    P, edge_gen, finest, corners = build_multiscale(level)
    n = len(P)
    ew_ms = {e: lam ** g for e, g in edge_gen.items()}
    ew_comb = {e: 1.0 for e in finest}
    Lms = laplacian(n, ew_ms)
    Lc = laplacian(n, ew_comb)
    u_ms, I_ms = solve_harmonic(Lms, corners)
    u_c, I_c = solve_harmonic(Lc, corners)
    return {
        "level": level,
        "N": n,
        "nF": corner_flux_norm(P, corners, I_ms),
        "nF_comb": corner_flux_norm(P, corners, I_c),
        "I": I_ms,
        "I_comb": I_c,
        "sumI": float(I_ms.sum()),
        "u_err": float(np.max(np.abs(u_ms - u_c))),
        "pred": closed_form(lam, level),
        "n_edges": len(edge_gen),
        "n_finest": len(finest),
    }


def main() -> int:
    fail = 0

    def check(ok: bool, msg: str) -> None:
        nonlocal fail
        print(("ok   " if ok else "FAIL ") + msg)
        if not ok:
            fail += 1

    print("=== Theorem W: multi-scale harmonic audit closed form ===")
    lambdas = (1.0, 1.2, 1.5, RESIST, 2.0, 3.0, 4.0)
    for lam in lambdas:
        print(f"\nλ = {lam}")
        for level in range(0, 6):
            row = audit_multiscale(level, lam)
            err = abs(row["nF"] - row["pred"])
            check(
                err < 1e-9 * max(1.0, row["pred"]),
                f"n={level}: ||F||={row['nF']:.12f} = Σ (λ·3/5)^k √21/2 "
                f"(err {err:.2e})",
            )
            check(abs(row["sumI"]) < 1e-9, f"n={level}: neutral Σ I = {row['sumI']:.2e}")
            check(
                row["u_err"] < 1e-9,
                f"n={level}: multi-scale harmonic = combinatorial "
                f"(||u-u_c||_∞={row['u_err']:.2e})",
            )

    print("\n=== Nested harmonicity: L^(g) u_c = 0 on fine interior ===")
    # Spot-check: each generation block annihilates the comb harmonic.
    for level in (3, 4):
        P, edge_gen, finest, corners = build_multiscale(level)
        n = len(P)
        ew_comb = {e: 1.0 for e in finest}
        u_c, _ = solve_harmonic(laplacian(n, ew_comb), corners)
        mask = np.ones(n, dtype=bool)
        mask[corners] = False
        for g in range(0, level + 1):
            ew_g = {e: 1.0 for e, gg in edge_gen.items() if gg == g}
            Lg = laplacian(n, ew_g)
            resid = float(np.max(np.abs((Lg @ u_c)[mask])))
            check(
                resid < 1e-9,
                f"n={level} gen {g}: |(L^(g) u_c)_int|_∞ = {resid:.2e}",
            )

    print("\n=== Trichotomy / Kept Failure W.1 ===")
    # λ < 5/3: monotone increase to finite limit
    for lam in (1.0, 1.2, 1.5):
        vals = [audit_multiscale(n, lam)["nF"] for n in range(0, 6)]
        lim = limit_if_convergent(lam)
        assert lim is not None
        check(
            all(vals[k + 1] > vals[k] for k in range(5)),
            f"λ={lam}: ||F|| strictly increasing through n=5",
        )
        check(
            vals[-1] < lim,
            f"λ={lam}: n=5 value {vals[-1]:.6f} < limit {lim:.6f} "
            f"= √21/2/(1-3λ/5)",
        )
        # remainder after n=5 is r^6/(1-r) * √21/2
        r = lam * 3.0 / 5.0
        rem = SQRT21_2 * (r ** 6) / (1.0 - r)
        check(
            abs((lim - vals[-1]) - rem) < 1e-12,
            f"λ={lam}: remainder lim-||F(5)|| = r^6/(1-r)·√21/2 ({rem:.6e})",
        )

    # λ = 5/3: exact linear growth (n+1) √21/2 — does NOT freeze
    resist_vals = []
    for n in range(0, 7):
        row = audit_multiscale(n, RESIST)
        resist_vals.append(row["nF"])
        check(
            abs(row["nF"] - (n + 1) * SQRT21_2) < 1e-9,
            f"λ=5/3 n={n}: ||F||={(n+1)}·√21/2 = {row['nF']:.12f}",
        )
    check(
        resist_vals[-1] > 7.0 * SQRT21_2 - 1e-9,
        f"Kept Failure W.1: resistance multi-scale ||F(6)||="
        f"{resist_vals[-1]:.6f} grows (not frozen at √21/2="
        f"{SQRT21_2:.6f})",
    )

    # λ > 5/3: diverges; successive increments grow by exactly λ·3/5
    for lam in (2.0, 4.0):
        vals = [audit_multiscale(n, lam)["nF"] for n in range(0, 6)]
        r = lam * 3.0 / 5.0
        incs = [vals[k + 1] - vals[k] for k in range(5)]
        check(
            all(vals[k + 1] > vals[k] for k in range(5)),
            f"λ={lam}: ||F|| strictly increasing",
        )
        check(
            all(abs(incs[k + 1] / incs[k] - r) < 1e-9 for k in range(4)),
            f"λ={lam}: successive increments grow by λ·3/5={r:.6f} "
            f"(inc ratios {[incs[k+1]/incs[k] for k in range(4)]})",
        )
        check(
            vals[-1] > 9.0 * SQRT21_2,
            f"λ={lam}: ||F(5)||={vals[-1]:.6f} already ≫ √21/2 (diverges)",
        )

    # Contrast with finest-only Theorem S at resistance
    from gasket_uniform_weights import harmonic_audit as finest_audit

    fo = finest_audit(5, RESIST)
    check(
        abs(fo["nF"] - SQRT21_2) < 1e-9,
        f"contrast Theorem S: finest-only λ=5/3 stays √21/2 ({fo['nF']:.12f})",
    )
    check(
        resist_vals[5] > 5.0 * SQRT21_2,
        "contrast: multi-scale λ=5/3 at n=5 is already > 5·√21/2",
    )

    print()
    if fail:
        print(f"{fail} check(s) failed")
        return 1
    print("all Theorem W checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
