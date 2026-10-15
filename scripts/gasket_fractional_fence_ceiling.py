#!/usr/bin/env python3
"""Boundary-fence resistance ceiling via Bloch spectrum (Theorem AN).

Witness for Theorem AN and Kept Failure AN.1 in GASKET.md.

Theorems AA/AL compared the gasket with ONE side (a path). Theorem AM used the
outer boundary cycle C_N. Theorem AN keeps the cycle and adds, on every outer
edge, the two spokes to that edge's inward unit-triangle apex — the boundary
fence / collar. The abstract translation-invariant fence (apex on every edge)
has an explicit 2x2 Bloch symbol H(θ); its corner resistance satisfies

  R^fen_α(m) < R^fen_∞(α) := (1/π) ∫_0^{2π} G_xx(θ;α) dθ

for every α ∈ (0,1/2), with G_xx = (H(θ)^{-α})_{xx}. At α=0.45,
R^fen_∞ ≈ 4.502104 (against AM's R^cyc_∞ ≈ 7.284859).

The gasket's boundary collar (cycle + ears at the three corners + spokes on the
N−6 non-corner outer edges) is a subgraph of G_n, so M ≤ L and M^α ≤ L^α. It
differs from the abstract fence by a three-corner defect whose resistance cost
δ(m) = R^col(m) − R^fen(m) stays below 1/2 for every checked m (increasing to
≈0.406). Hence

  R^col_α(m) < R^fen_∞(α) + 1/2 ≈ 5.002104,

and the three-corner pairing gives Q ≥ (7/3)/R^col_∞, so at α=0.45
Q ≥ 0.466470, ||F_frac|| ≥ 0.305377 (≈1.456× AM), and with σ_AK
||F(n,0.45)|| ≥ 0.214239.

Kept Failure AN.1: still below Conjecture J.1's 0.70; AM.1 already showed no
resistance ceiling rescues the σ_AK chain.

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
    R_cyc_inf,
    cycle_R,
    pinv_power,
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
AM_FRAC_FLOOR_045 = 0.209685
AM_F_FLOOR_045 = 0.147106
DELTA_CAP = 0.5  # δ(m) < 1/2 for all m (witnessed; lim ≈ 0.406)

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


# ---------------------------------------------------------------- Bloch fence
def fence_evals(theta: np.ndarray):
    """Eigenvalues and G_xx of the 2x2 Bloch symbol of the infinite fence."""
    a = 4.0 - 2.0 * np.cos(theta)
    tr = a + 2.0
    det = 12.0 * np.sin(theta / 2.0) ** 2
    disc = np.maximum(tr * tr - 4.0 * det, 0.0)
    s = np.sqrt(disc)
    lp = 0.5 * (tr + s)
    lm = 0.5 * (tr - s)
    return lp, lm, a, tr


def Gxx(theta: np.ndarray, alpha: float) -> np.ndarray:
    lp, lm, a, tr = fence_evals(theta)
    fp = np.where(lp > 1e-30, lp ** (-alpha), 0.0)
    fm = np.where(lm > 1e-30, lm ** (-alpha), 0.0)
    c0 = 0.5 * (fp + fm)
    dlam = lp - lm
    c1 = np.where(np.abs(dlam) > 1e-18, (fp - fm) / dlam, 0.0)
    return c0 + c1 * (a - tr / 2.0)


def R_fen_inf(alpha: float, nodes: int = 2_000_000) -> float:
    """(1/π)∫_0^{2π} G_xx = 2∫_0^1 G_xx(2πt) dt, singularity removed."""
    # lim_{θ→0} G_xx(θ) θ^{2α} = (1/2)(3/4)^{-α}  (λ_min ∼ (3/4)θ², |U_x|²=1/2)
    kappa = 0.5 * (0.75) ** (-alpha)
    A = kappa * (2.0 * math.pi) ** (-2.0 * alpha)
    us = (np.arange(nodes) + 0.5) / nodes
    G = Gxx(2.0 * math.pi * us, alpha)
    H = G - A * (us ** (-2.0 * alpha) + (1.0 - us) ** (-2.0 * alpha))
    return float(2.0 * (np.mean(H) + 2.0 * A / (1.0 - 2.0 * alpha)))


def fence_R_fourier(m: int, alpha: float) -> float:
    """R^fen_α(m) on C_{3m} with an apex on every edge (Bloch / DFT)."""
    N = 3 * m
    p = np.arange(1, N, dtype=np.float64)
    th = 2.0 * math.pi * p / N
    Vp2 = 4.0 * np.sin(math.pi * p / 3.0) ** 2 / N
    return float(np.sum(Vp2 * Gxx(th, alpha)))


def fence_laplacian(m: int) -> np.ndarray:
    """Dense Laplacian of the abstract fence on C_{3m}."""
    N = 3 * m
    V = 2 * N
    Mat = np.zeros((V, V))

    def add(u, v):
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
    return Mat


def collar_laplacian(m: int) -> np.ndarray:
    """Cycle + three ears + apexes on the N−6 non-corner outer edges."""
    N = 3 * m
    corner_adj = set()
    for c in (0, m, 2 * m):
        corner_adj.add((min((c - 1) % N, c), max((c - 1) % N, c)))
        corner_adj.add((min(c, (c + 1) % N), max(c, (c + 1) % N)))
    apex_edges = []
    for i in range(N):
        j = (i + 1) % N
        e = (min(i, j), max(i, j))
        if e not in corner_adj:
            apex_edges.append((i, j))
    V = N + len(apex_edges)
    Mat = np.zeros((V, V))

    def add(u, v):
        Mat[u, u] += 1.0
        Mat[v, v] += 1.0
        Mat[u, v] -= 1.0
        Mat[v, u] -= 1.0

    for i in range(N):
        add(i, (i + 1) % N)
    for c in (0, m, 2 * m):
        add((c - 1) % N, (c + 1) % N)
    for ai, (i, j) in enumerate(apex_edges):
        a = N + ai
        add(a, i)
        add(a, j)
    return Mat


def R_of_laplacian(Mat: np.ndarray, m: int, alpha: float) -> float:
    _, Gp = pinv_power(Mat, alpha)
    d = np.zeros(Mat.shape[0])
    d[0] = 1.0
    d[m] = -1.0
    return float(d @ Gp @ d)


def Q_floor_an(alpha: float) -> float:
    return (7.0 / 3.0) / (R_fen_inf(alpha) + DELTA_CAP)


def Ffrac_floor_an(alpha: float) -> float:
    return Q_floor_an(alpha) * SQRT21_7


def F_floor_an(alpha: float) -> float:
    return sigma_AK(alpha) * Ffrac_floor_an(alpha)


# --------------------------------------------------------------- gasket collar
def barycentric(p, c0, c1, c2):
    A = np.array([c1 - c0, c2 - c0]).T
    xy = np.linalg.solve(A, p - c0)
    b1, b2 = xy
    return np.array([1.0 - b1 - b2, b1, b2])


def gasket_collar_edges(P, edges, corners):
    """Edges of all unit triangles that meet the outer boundary."""
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
    bedges: set[tuple[int, int]] = set()
    seen: set[tuple[int, ...]] = set()
    for i, j in edges:
        if i not in bv and j not in bv:
            continue
        for k in adj[i] & adj[j]:
            key = tuple(sorted((i, j, k)))
            if key in seen:
                continue
            ds = [
                np.linalg.norm(P[u] - P[v])
                for u, v in ((i, j), (j, k), (k, i))
            ]
            if all(abs(d - el) < 1e-5 for d in ds) and sum(
                v in bv for v in key
            ) >= 2:
                seen.add(key)
                for u, v in ((i, j), (j, k), (k, i)):
                    bedges.add((min(u, v), max(u, v)))
    return sorted(bedges), bv, seen



def level_collar_record(level: int, alpha: float) -> dict:
    """Boundary-collar = edges of unit triangles meeting the outer boundary."""
    P, edges, corners = build_gasket(level)
    N = len(P)
    m = 2 ** level
    L = combinatorial_laplacian(N, edges)
    bedges, bv, tris = gasket_collar_edges(P, edges, corners)
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
    R_fen = fence_R_fourier(m, alpha)
    R_abs = R_of_laplacian(collar_laplacian(m), m, alpha) if m >= 2 else R_col
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
        "R_fen": R_fen,
        "R_abs": R_abs,
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
    Rfen = R_fen_inf(alpha)
    Rceil = Rfen + DELTA_CAP
    Qf = Q_floor_an(alpha)
    Fff = Ffrac_floor_an(alpha)
    Ff = F_floor_an(alpha)
    print(f"alpha={alpha}")
    print(f"R_fen_inf={Rfen:.12f}")
    print(f"R_ceil=R_fen_inf+1/2={Rceil:.12f}")
    print(f"R_cyc_inf={R_cyc_inf(alpha):.12f}")
    print(f"Q_floor={Qf:.12f}")
    print(f"Ffrac_floor={Fff:.12f}")
    print(f"sigma_AK={sigma_AK(alpha):.12f}")
    print(f"F_floor={Ff:.12f}")

    # --- Bloch / Fourier consistency ---
    for m in (2, 4, 8, 16):
        Rf = fence_R_fourier(m, alpha)
        Rd = R_of_laplacian(fence_laplacian(m), m, alpha)
        check(abs(Rf - Rd) < 1e-8, f"fence Fourier=dense at m={m}: {Rf:.8f} vs {Rd:.8f}")

    # --- R_fen(m) < R_fen_inf ---
    for m in [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]:
        Rf = fence_R_fourier(m, alpha)
        check(Rf < Rfen, f"R_fen({m})={Rf:.6f} < R_fen_inf={Rfen:.6f}")

    # --- small-θ asymptotics λ_min ∼ (3/4) θ² ---
    th = 1e-5
    lp, lm, _, _ = fence_evals(np.array([th]))
    check(abs(float(lm[0]) / (th * th) - 0.75) < 1e-3, f"λ_min/θ² → 3/4: {float(lm[0])/(th*th):.6f}")

    # --- singularity-removed integral stable ---
    Rfen2 = R_fen_inf(alpha, nodes=500_000)
    check(abs(Rfen - Rfen2) < 1e-4, f"R_fen_inf stable: {Rfen:.8f} vs {Rfen2:.8f}")

    # --- collar vs fence defect δ(m) < 1/2 ---
    max_delta = 0.0
    for m in [2, 4, 8, 16, 32, 64, 128, 256]:
        # m=2: collar = cycle+3 ears (no apexes); fence has 6 apexes
        Rf = fence_R_fourier(m, alpha)
        if m == 2:
            # abstract collar at m=2: cycle + 3 ears = G_1 wireframe on 6 verts
            Mat = collar_laplacian(m)
            Rc = R_of_laplacian(Mat, m, alpha)
        else:
            Rc = R_of_laplacian(collar_laplacian(m), m, alpha)
        delta = Rc - Rf
        max_delta = max(max_delta, delta)
        check(delta < DELTA_CAP, f"δ({m})={delta:.6f} < 1/2")
        check(Rc < Rceil, f"R_col({m})={Rc:.6f} < R_ceil={Rceil:.6f}")
        check(Rf < Rc, f"R_fen({m})={Rf:.6f} < R_col({m})={Rc:.6f}")
        print(f"  m={m:4d} R_fen={Rf:.6f} R_col={Rc:.6f} δ={delta:.6f}")
    check(max_delta < 0.42, f"max δ on dyadics ≤256 is {max_delta:.6f} < 0.42")

    # --- R_col < R_ceil implies the floors ---
    check(Qf > 0.46647 and Qf < 0.46648, f"Q_floor digits {Qf:.12f}")
    check(Fff > 0.30537 and Fff < 0.30538, f"Ffrac_floor digits {Fff:.12f}")
    check(Ff > 0.21423 and Ff < 0.21424, f"F_floor digits {Ff:.12f}")
    check(Fff > AM_FRAC_FLOOR_045 * 1.45, f"Ffrac floor beats AM by ≥1.45×: {Fff/AM_FRAC_FLOOR_045:.4f}")

    # --- live gasket levels 1–5 ---
    for n in range(1, 6):
        rec = level_collar_record(n, alpha)
        print(
            f"n={n} R_live={rec['R_live']:.6f} R_col={rec['R_col']:.6f} "
            f"R_fen={rec['R_fen']:.6f} R_cyc={rec['R_cyc']:.6f} "
            f"Q={rec['Q']:.6f} Q≥(7/3)/R_col={rec['Q_floor_level']:.6f} "
            f"|E|={rec['n_bedges']}"
        )
        check(rec["L_minus_M_psd"], f"L−M ⪰ 0 at level {n}")
        check(rec["La_minus_Ma_psd"], f"L^α−M^α ⪰ 0 at level {n}")
        check(rec["R_col"] < Rceil, f"live R_col({n})={rec['R_col']:.6f} < R_ceil")
        check(rec["R_live"] <= rec["R_col"] + 1e-9, f"live R ≤ R_col at level {n}")
        check(rec["R_col"] <= rec["R_cyc"] + 1e-6, f"R_col ≤ R_cyc at level {n}")
        check(rec["R_fen"] < rec["R_col"] + 1e-6, f"R_fen < R_col at level {n}")
        check(rec["R_col"] < Rceil, f"gasket R_col < R_ceil at level {n}")
        check(rec["R_abs"] < Rceil, f"abstract R_col < R_ceil at level {n}")
        check(rec["Q"] + 1e-9 >= rec["Q_floor_level"], f"Q ≥ (7/3)/R_col at level {n}")
        check(rec["Q"] >= Qf, f"Q ≥ Q_floor at level {n}: {rec['Q']:.6f}")
        check(rec["nFf"] >= Fff, f"||F_frac|| ≥ floor at level {n}: {rec['nFf']:.6f}")
        check(rec["nF"] >= Ff, f"||F|| ≥ floor at level {n}: {rec['nF']:.6f}")

    # --- Kept Failure AN.1 ---
    check(Ff < J1_VALUE, f"AN floor {Ff:.6f} < J.1's {J1_VALUE}")
    env = 2.0 ** (1.0 - alpha) * Fff
    check(env < J1_VALUE, f"2^(1-α)·Ffrac={env:.6f} < {J1_VALUE}")
    live_share = 1.3525 * Fff
    check(live_share < J1_VALUE, f"live-σ·Ffrac≈{live_share:.6f} < {J1_VALUE}")
    # even with the sharp extrapolated R_col_∞ ≈ 4.91 (δ≈0.406)
    Rsharp = Rfen + 0.42
    Fsharp = sigma_AK(alpha) * ((7.0 / 3.0) / Rsharp) * SQRT21_7
    check(Fsharp < J1_VALUE, f"sharp-collar σ_AK floor {Fsharp:.6f} < {J1_VALUE}")
    print(f"Kept Failure AN.1: AN floor={Ff:.6f}, envelope={env:.6f}, sharp={Fsharp:.6f}, all < 0.70")

    if failures:
        print(f"\n{len(failures)} FAILURE(S)")
        for f in failures:
            print(" ", f)
        return 1
    print("\nAll checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
