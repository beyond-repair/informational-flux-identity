#!/usr/bin/env python3
"""Depth-2 collar resistance ceiling via period-2 Bloch spectrum (Theorem AO).

Witness for Theorem AO and Kept Failure AO.1 in GASKET.md.

Theorem AN used the depth-1 boundary collar (unit triangles meeting the outer
boundary) and compared it to the abstract fence (apex on every outer edge).
Theorem AO thickens that strip by one layer: the depth-2 collar is the edges of
all unit triangles that share an edge with a depth-1 collar triangle. Along each
outer side that adds every-other apex–apex edge (the SG holes omit the odd
pairs), so the translation-invariant model is a period-2 fence whose 4×4 Bloch
symbol Ĥ(φ) has

  R^{d2}_α(m) < R^{d2}_∞(α) := (1/π) ∫_0^{2π} G_{00}(φ;α) dφ

for every α ∈ (0,1/2), with G_{00} = (Ĥ(φ)^{-α})_{x0,x0}. At α=0.45,
R^{d2}_∞ ≈ 4.270047 (against AN's R^fen_∞ ≈ 4.502104).

The gasket depth-2 collar is a subgraph of G_n, so M ≤ L and M^α ≤ L^α. It
differs from the abstract period-2 fence by a three-corner defect whose
resistance cost δ₂(m) = R^{col2}(m) − R^{d2}(m) stays below 0.55 for every
checked m (increasing toward ≈0.525). Hence

  R^{col2}_α(m) < R^{d2}_∞(α) + 0.55 ≈ 4.820047,

and the three-corner pairing gives Q ≥ (7/3)/R^{col2}_∞, so at α=0.45
Q ≥ 0.484089, ||F_frac|| ≥ 0.316911 (≈1.038× AN), and with σ_AK
||F(n,0.45)|| ≥ 0.222331.

Kept Failure AO.1: still below Conjecture J.1's 0.70; AN.1/AM.1 already showed
no resistance ceiling rescues the σ_AK chain. The depth-2 gain over the fence
is thin (live R≈2.35 still far below 4.82).

Not thrust.
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_cycle_ceiling import (  # noqa: E402
    cycle_R,
    pinv_power,
)
from gasket_fractional_fence_ceiling import (  # noqa: E402
    R_fen_inf,
    fence_R_fourier,
)
from gasket_fractional_limit import (  # noqa: E402
    BC,
    build_gasket,
    combinatorial_laplacian,
    corner_flux_norm,
)
from gasket_fractional_moment_jensen import sigma_AK  # noqa: E402

ALPHA = 0.45
J1_VALUE = 0.70
SQRT21_7 = math.sqrt(21.0) / 7.0
AN_FRAC_FLOOR_045 = 0.305377
AN_F_FLOOR_045 = 0.214239
DELTA_CAP = 0.55  # δ₂(m) < 0.55 for all checked m (witnessed; lim ≈ 0.525)

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


# -------------------------------------------------------- period-2 Bloch fence
def bloch_H(theta: float) -> np.ndarray:
    """4×4 Bloch symbol, basis (x0, x1, a0, a1); supercell of two boundary steps.

    Edges: boundary x0—x1 and x1—x0'; fence spokes a0—{x0,x1}, a1—{x1,x0'};
    plus the depth-2 apex–apex a1—a0' at every even boundary vertex.
    """
    e = np.exp(1j * theta)
    H = np.zeros((4, 4), dtype=complex)

    def add_edge(i: int, j: int, ph: complex = 1 + 0j) -> None:
        H[i, i] += 1.0
        H[j, j] += 1.0
        H[i, j] -= ph
        H[j, i] -= np.conj(ph)

    add_edge(0, 1)  # x0—x1
    add_edge(1, 0, e)  # x1—x0'
    add_edge(2, 0)  # a0—x0
    add_edge(2, 1)  # a0—x1
    add_edge(3, 1)  # a1—x1
    add_edge(3, 0, e)  # a1—x0'
    add_edge(3, 2, e)  # a1—a0' (every-other apex–apex)
    return H


def G00(theta: float, alpha: float) -> float:
    """(Ĥ(θ)^{-α})_{x0,x0}."""
    H = bloch_H(theta)
    evals, evecs = np.linalg.eigh(H)
    g = 0.0
    for i, lam in enumerate(evals):
        if lam > 1e-14:
            g += float(lam ** (-alpha) * abs(evecs[0, i]) ** 2)
    return g


def R_d2_inf(alpha: float, nodes: int = 200_000) -> float:
    """(1/π)∫_0^{2π} G00 = 2∫_0^1 G00(2πt) dt, singularity removed.

    Goldstone: λ_min ∼ (7/32) θ², so G00 ∼ κ θ^{-2α} with κ from small θ.
    """
    th = 1e-6
    kappa = G00(th, alpha) * (th ** (2.0 * alpha))
    A = kappa * (2.0 * math.pi) ** (-2.0 * alpha)
    us = (np.arange(nodes) + 0.5) / nodes
    Gs = np.array([G00(2.0 * math.pi * float(u), alpha) for u in us])
    H = Gs - A * (us ** (-2.0 * alpha) + (1.0 - us) ** (-2.0 * alpha))
    return float(2.0 * (np.mean(H) + 2.0 * A / (1.0 - 2.0 * alpha)))


def depth2_fence_laplacian(m: int) -> np.ndarray:
    """Dense Laplacian of the abstract period-2 depth-2 fence on C_{3m}."""
    N = 3 * m
    V = 2 * N
    Mat = np.zeros((V, V))

    def add(u: int, v: int) -> None:
        Mat[u, u] += 1.0
        Mat[v, v] += 1.0
        Mat[u, v] -= 1.0
        Mat[v, u] -= 1.0

    for i in range(N):
        j = (i + 1) % N
        add(i, j)
        a = N + i
        add(a, i)
        add(a, j)
    for k in range(0, N, 2):
        add(N + ((k - 1) % N), N + k)
    return Mat


def R_of_laplacian(Mat: np.ndarray, m: int, alpha: float) -> float:
    _, Gp = pinv_power(Mat, alpha)
    d = np.zeros(Mat.shape[0])
    d[0] = 1.0
    d[m] = -1.0
    return float(d @ Gp @ d)


def depth2_R(m: int, alpha: float) -> float:
    return R_of_laplacian(depth2_fence_laplacian(m), m, alpha)


def Q_floor_ao(alpha: float, Rd2: float | None = None) -> float:
    if Rd2 is None:
        Rd2 = R_d2_inf(alpha)
    return (7.0 / 3.0) / (Rd2 + DELTA_CAP)


def Ffrac_floor_ao(alpha: float, Rd2: float | None = None) -> float:
    return Q_floor_ao(alpha, Rd2) * SQRT21_7


def F_floor_ao(alpha: float, Rd2: float | None = None) -> float:
    return sigma_AK(alpha) * Ffrac_floor_ao(alpha, Rd2)


# --------------------------------------------------------------- gasket collar
def barycentric(p, c0, c1, c2):
    A = np.array([c1 - c0, c2 - c0]).T
    xy = np.linalg.solve(A, p - c0)
    b1, b2 = xy
    return np.array([1.0 - b1 - b2, b1, b2])


def depth2_collar_edges(P, edges, corners):
    """Edges of unit triangles at graph-distance ≤2 from the outer boundary.

    Depth 1: triangles with ≥2 verts on the outer boundary (AN collar).
    Depth 2: also triangles that share an edge with a depth-1 triangle.
    """
    a, b, c = (int(z) for z in corners)
    c0, c1, c2 = P[a], P[b], P[c]
    bv = {
        i
        for i, p in enumerate(P)
        if np.min(np.abs(barycentric(p, c0, c1, c2))) < 1e-8
    }
    adj: dict[int, set[int]] = defaultdict(set)
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    el = float(
        np.median([np.linalg.norm(P[i] - P[j]) for i, j in list(edges)[:40]])
    )
    tris: list[tuple[int, int, int]] = []
    seen: set[tuple[int, ...]] = set()
    for i, j in edges:
        for k in adj[i] & adj[j]:
            key = tuple(sorted((i, j, k)))
            if key in seen:
                continue
            ds = [
                np.linalg.norm(P[u] - P[v])
                for u, v in ((i, j), (j, k), (k, i))
            ]
            if all(abs(d - el) < 1e-5 for d in ds):
                seen.add(key)
                tris.append(key)  # type: ignore[arg-type]
    coll = {t for t in tris if sum(v in bv for v in t) >= 2}
    new: set[tuple[int, ...]] = set()
    for t in tris:
        if t in coll:
            continue
        et = {
            (min(t[0], t[1]), max(t[0], t[1])),
            (min(t[1], t[2]), max(t[1], t[2])),
            (min(t[2], t[0]), max(t[2], t[0])),
        }
        for prev in coll:
            ep = {
                (min(prev[0], prev[1]), max(prev[0], prev[1])),
                (min(prev[1], prev[2]), max(prev[1], prev[2])),
                (min(prev[2], prev[0]), max(prev[2], prev[0])),
            }
            if et & ep:
                new.add(t)
                break
    coll |= new
    bedges: set[tuple[int, int]] = set()
    for t in coll:
        for u, v in ((t[0], t[1]), (t[1], t[2]), (t[2], t[0])):
            bedges.add((min(u, v), max(u, v)))
    return sorted(bedges), coll


def level_depth2_record(level: int, alpha: float) -> dict:
    P, edges, corners = build_gasket(level)
    N = len(P)
    m = 2 ** level
    L = combinatorial_laplacian(N, edges)
    bedges, tris = depth2_collar_edges(P, edges, corners)
    Eset = {(min(i, j), max(i, j)) for i, j in edges}
    assert set(bedges) <= Eset
    M = combinatorial_laplacian(N, bedges)
    A, Lp = pinv_power(L, alpha)
    Ma, Mp = pinv_power(M, alpha)
    a, b, c = (int(z) for z in corners)
    v = np.zeros(N)
    v[a], v[b] = 1.0, -1.0
    interior = np.setdiff1d(np.arange(N), corners)
    u = np.zeros(N)
    u[corners] = BC
    u[interior] = np.linalg.solve(
        A[np.ix_(interior, interior)], -A[np.ix_(interior, corners)] @ BC
    )
    R_live = float(v @ Lp @ v)
    R_col = float(v @ Mp @ v)
    R_d2 = depth2_R(m, alpha)
    R_fen = fence_R_fourier(m, alpha)
    Q = float(u @ (A @ u))
    Q_floor_level = (7.0 / 3.0) / R_col
    If = (A @ u)[corners]
    I = (L @ u)[corners]
    nF = corner_flux_norm(P, corners, I)
    nFf = corner_flux_norm(P, corners, If)
    return {
        "level": level,
        "m": m,
        "R_live": R_live,
        "R_col": R_col,
        "R_d2": R_d2,
        "R_fen": R_fen,
        "R_cyc": cycle_R(m, alpha),
        "Q": Q,
        "Q_floor_level": Q_floor_level,
        "nF": nF,
        "nFf": nFf,
        "n_bedges": len(bedges),
        "n_tris": len(tris),
        "L_minus_M_psd": bool(np.min(np.linalg.eigvalsh(L - M)) >= -1e-9),
        "La_minus_Ma_psd": bool(np.min(np.linalg.eigvalsh(A - Ma)) >= -1e-8),
    }


def main() -> int:
    alpha = ALPHA
    print("Computing R_d2_inf (200k nodes)...", flush=True)
    Rd2 = R_d2_inf(alpha, nodes=200_000)
    Rceil = Rd2 + DELTA_CAP
    Qf = Q_floor_ao(alpha, Rd2)
    Fff = Ffrac_floor_ao(alpha, Rd2)
    Ff = F_floor_ao(alpha, Rd2)
    Rfen = R_fen_inf(alpha)
    print(f"alpha={alpha}")
    print(f"R_d2_inf={Rd2:.12f}")
    print(f"R_ceil=R_d2_inf+0.55={Rceil:.12f}")
    print(f"R_fen_inf={Rfen:.12f}")
    print(f"Q_floor={Qf:.12f}")
    print(f"Ffrac_floor={Fff:.12f}")
    print(f"sigma_AK={sigma_AK(alpha):.12f}")
    print(f"F_floor={Ff:.12f}")

    # --- small-θ asymptotics λ_min ∼ (7/32) θ² ---
    th = 1e-5
    lm = float(np.linalg.eigvalsh(bloch_H(th)).min())
    check(
        abs(lm / (th * th) - 7.0 / 32.0) < 1e-3,
        f"λ_min/θ² → 7/32: {lm/(th*th):.6f}",
    )

    # --- singularity-removed integral stable ---
    print("Stability check (50k nodes)...", flush=True)
    Rd2b = R_d2_inf(alpha, nodes=50_000)
    check(abs(Rd2 - Rd2b) < 5e-4, f"R_d2_inf stable: {Rd2:.8f} vs {Rd2b:.8f}")
    check(Rd2 < Rfen, f"R_d2_inf={Rd2:.6f} < R_fen_inf={Rfen:.6f}")

    # --- R_d2(m) < R_d2_inf ---
    for m in [2, 4, 8, 16, 32, 64, 128, 256]:
        Rd = depth2_R(m, alpha)
        check(Rd < Rd2, f"R_d2({m})={Rd:.6f} < R_d2_inf={Rd2:.6f}")
        Rf = fence_R_fourier(m, alpha)
        check(Rd < Rf + 1e-9, f"R_d2({m})={Rd:.6f} ≤ R_fen({m})={Rf:.6f}")

    # --- defect δ₂(m) < 0.55 on live gasket levels ---
    max_delta = 0.0
    for n in range(2, 7):
        rec = level_depth2_record(n, alpha)
        delta = rec["R_col"] - rec["R_d2"]
        max_delta = max(max_delta, delta)
        print(
            f"n={n} R_live={rec['R_live']:.6f} R_col2={rec['R_col']:.6f} "
            f"R_d2={rec['R_d2']:.6f} R_fen={rec['R_fen']:.6f} "
            f"δ2={delta:.6f} Q={rec['Q']:.6f} |E|={rec['n_bedges']} tris={rec['n_tris']}"
        )
        check(delta < DELTA_CAP, f"δ2({n})={delta:.6f} < 0.55")
        check(rec["R_col"] < Rceil, f"R_col2({n})={rec['R_col']:.6f} < R_ceil")
        check(rec["R_d2"] < rec["R_col"] + 1e-9, f"R_d2 < R_col2 at level {n}")
        check(rec["L_minus_M_psd"], f"L−M ⪰ 0 at level {n}")
        check(rec["La_minus_Ma_psd"], f"L^α−M^α ⪰ 0 at level {n}")
        check(rec["R_live"] <= rec["R_col"] + 1e-9, f"live R ≤ R_col2 at level {n}")
        # Abstract fence R_fen can sit below gasket R_col2 (corner apexes vs ears);
        # the theorem compares R_col2 to R_d2_∞+δ, not to R_fen.
        check(rec["R_col"] < rec["R_cyc"] + 1e-6, f"R_col2 < R_cyc at level {n}")
        check(rec["Q"] + 1e-9 >= rec["Q_floor_level"], f"Q ≥ (7/3)/R_col2 at level {n}")
        check(rec["Q"] >= Qf, f"Q ≥ Q_floor at level {n}: {rec['Q']:.6f}")
        check(rec["nFf"] >= Fff, f"||F_frac|| ≥ floor at level {n}: {rec['nFf']:.6f}")
        check(rec["nF"] >= Ff, f"||F|| ≥ floor at level {n}: {rec['nF']:.6f}")
    check(max_delta < 0.50, f"max δ2 on levels 2–6 is {max_delta:.6f} < 0.50")

    # --- floor digits ---
    check(Qf > 0.48408 and Qf < 0.48410, f"Q_floor digits {Qf:.12f}")
    check(Fff > 0.31691 and Fff < 0.31692, f"Ffrac_floor digits {Fff:.12f}")
    check(Ff > 0.22233 and Ff < 0.22234, f"F_floor digits {Ff:.12f}")
    check(
        Fff > AN_FRAC_FLOOR_045 * 1.03,
        f"Ffrac floor beats AN by ≥1.03×: {Fff/AN_FRAC_FLOOR_045:.4f}",
    )
    check(Ff > AN_F_FLOOR_045, f"F floor beats AN: {Ff:.6f} > {AN_F_FLOOR_045}")

    # --- Kept Failure AO.1 ---
    check(Ff < J1_VALUE, f"AO floor {Ff:.6f} < J.1's {J1_VALUE}")
    env = 2.0 ** (1.0 - alpha) * Fff
    check(env < J1_VALUE, f"2^(1-α)·Ffrac={env:.6f} < {J1_VALUE}")
    live_share = 1.3525 * Fff
    check(live_share < J1_VALUE, f"live-σ·Ffrac≈{live_share:.6f} < {J1_VALUE}")
    # even with extrapolated δ₂≈0.53
    Rsharp = Rd2 + 0.53
    Fsharp = sigma_AK(alpha) * ((7.0 / 3.0) / Rsharp) * SQRT21_7
    check(Fsharp < J1_VALUE, f"sharp-d2 σ_AK floor {Fsharp:.6f} < {J1_VALUE}")
    print(
        f"Kept Failure AO.1: AO floor={Ff:.6f}, envelope={env:.6f}, "
        f"sharp={Fsharp:.6f}, all < 0.70"
    )

    if failures:
        print(f"\n{len(failures)} FAILURE(S)")
        for f in failures:
            print(" ", f)
        return 1
    print("\nAll checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
