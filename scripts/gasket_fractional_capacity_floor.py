#!/usr/bin/env python3
"""Uniform energy floor for spectral L^α on the gasket for every α < log3/log5.

Witness for Theorem AC and Kept Failure AC.1 in GASKET.md.

Theorem AA floored Q_α(u_n) only for α < 1/2, because it compared against a
boundary path (spectral dimension 1). Theorem AC uses the gasket itself:

  (1) cells: a level-k cell is a copy of G_k with (3^{k+1}+3)/2 vertices and
      effective-resistance diameter at most (8/3)(5/3)^k (corner resistance
      (2/3)(5/3)^k plus a nested-corner chain), and resistance in G_n between
      two points of one cell is at most the resistance inside the cell
      (Rayleigh monotonicity);
  (2) low modes are flat: if E(g) <= s||g||^2 then max g^2 <= K_s ||g||^2 with
      K_s = 2/|cell| + 2 R_k s, so the spectral mass of v = e_a - e_b obeys
      mu((0,s]) <= 4 K_s <= 80 s^γ, γ = log3/log5, for 5^{-(n+1)} <= s < 1;
  (3) gap: lambda_1(G_n) >= 2/(N R_max) >= 5^{-n}/4, so mu((0,s]) = 0 below;
  (4) layer cake: R_α(n) = v^T (L^α)^+ v = α ∫ s^{-α-1} mu((0,s]) ds
      <= C_AC(α) = (2γ/(γ-α)) 40^{α/γ}, uniformly in n, iff α < γ;
  (5) duality: Q_α(u_n) >= (3/2)^2 / R_α(n) >= (9/4)/C_AC(α).

Kept Failure AC.1: a level-uniform energy floor at every α ∈ (0,1) is false;
above γ Theorem K sends Q_α(u_n) to 0 (witnessed live at α = 0.75, 0.9).

It does NOT decide the combinatorial audit norm ||F(n,0.45)||. Not thrust.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_limit import (  # noqa: E402
    ALPHA_C,
    BC,
    build_gasket,
    combinatorial_laplacian,
    corner_flux_norm,
)

GAMMA = ALPHA_C  # log 3 / log 5
GAP_AB = 1.5
SQRT21_7 = math.sqrt(21.0) / 7.0
AA_FLOOR_045 = 0.128335  # Theorem AA at alpha = 0.45 (recorded)
Q_RECORDED = {9: 0.948091, 10: 0.941537}

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def C_AC(alpha: float) -> float:
    return (2.0 * GAMMA / (GAMMA - alpha)) * 40.0 ** (alpha / GAMMA)


def eig(level: int):
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    w, V = np.linalg.eigh(L)
    return P, edges, corners, N, L, w, V


def resistance_matrix(L: np.ndarray) -> np.ndarray:
    G = np.linalg.pinv(L, hermitian=True)
    d = np.diag(G)
    return d[:, None] + d[None, :] - 2.0 * G


def main() -> int:
    print("Theorem AC: energy floor for L^alpha on the gasket for every alpha < log3/log5")
    print(f"gamma = log3/log5 = {GAMMA:.6f}. Assumption: spectral calculus (GASKET item 4).")
    print("Not thrust. Does not decide the combinatorial audit norm ||F||.\n")

    # 1. Cell resistance: corner resistance (2/3)(5/3)^k, diameter <= (8/3)(5/3)^k, size.
    for k in range(0, 6):
        P, edges, corners, N, L, w, V = eig(k)
        R = resistance_matrix(L)
        rc = R[corners[0], corners[1]]
        check(abs(rc - (2.0 / 3.0) * (5.0 / 3.0) ** k) < 1e-9,
              f"G_{k}: corner resistance {rc:.9f} = (2/3)(5/3)^{k}")
        check(R.max() <= (8.0 / 3.0) * (5.0 / 3.0) ** k + 1e-9 and N == (3 ** (k + 1) + 3) // 2,
              f"G_{k}: resistance diameter {R.max():.6f} <= (8/3)(5/3)^{k} = {(8/3)*(5/3)**k:.6f}; N={N}")

    # 2-3. Spectral gap and spectral mass of v = e_a - e_b on live levels.
    for n in range(1, 7):
        P, edges, corners, N, L, w, V = eig(n)
        a, b = int(corners[0]), int(corners[1])
        lam1 = w[1]
        check(lam1 >= 5.0 ** (-n) / 4.0, f"level {n}: lambda_1 = {lam1:.6e} >= 5^-{n}/4 = {5.0**(-n)/4:.6e}")
        c2 = (V[a, :] - V[b, :]) ** 2
        pos = w > 1e-10
        lam, mass = w[pos], c2[pos]
        cum = np.cumsum(mass)
        sel = lam < 1.0
        ratio = float(np.max(cum[sel] / lam[sel] ** GAMMA)) if sel.any() else 0.0
        check(ratio <= 80.0 and abs(cum[-1] - 2.0) < 1e-10,
              f"level {n}: max_s<1 mu((0,s])/s^gamma = {ratio:.4f} <= 80; total mass {cum[-1]:.12f} = 2")

    # 4-5. Uniform R_alpha and the floor, including alpha in [1/2, gamma).
    alphas = (0.45, 0.50, 0.55, 0.60, 0.65, 0.68)
    for alpha in alphas:
        print(f"     alpha={alpha}: C_AC={C_AC(alpha):.6f}, Q floor={GAP_AB**2/C_AC(alpha):.6f}")
    levels = {}
    for n in range(1, 7):
        P, edges, corners, N, L, w, V = eig(n)
        levels[n] = (P, corners, N, w, V)
        a, b = int(corners[0]), int(corners[1])
        c2 = (V[a, :] - V[b, :]) ** 2
        pos = w > 1e-10
        interior = np.setdiff1d(np.arange(N), corners)
        for alpha in alphas:
            Ra = float(np.sum(c2[pos] * w[pos] ** (-alpha)))
            wa = np.where(pos, np.clip(w, 0, None) ** alpha, 0.0)
            A = (V * wa) @ V.T
            u = np.zeros(N)
            u[corners] = BC
            u[interior] = np.linalg.solve(A[np.ix_(interior, interior)], -A[np.ix_(interior, corners)] @ BC)
            Q = float(u @ A @ u)
            If = (A @ u)[corners]
            nFf = corner_flux_norm(P, corners, If)
            fl = GAP_AB ** 2 / C_AC(alpha)
            check(Ra <= C_AC(alpha) and Q >= GAP_AB ** 2 / Ra - 1e-10 and GAP_AB ** 2 / Ra >= fl,
                  f"level {n}, alpha={alpha}: R={Ra:.6f} <= C_AC; Q={Q:.6f} >= (9/4)/R={GAP_AB**2/Ra:.6f} >= {fl:.6f}")
            check(abs(nFf - Q * SQRT21_7) < 1e-9, f"level {n}, alpha={alpha}: ||F_frac||={nFf:.6f} = Q*sqrt21/7")
    for n, q in Q_RECORDED.items():
        check(q > GAP_AB ** 2 / C_AC(0.45), f"recorded level {n} energy Q_0.45={q} clears the AC floor")
    check(GAP_AB ** 2 / C_AC(0.45) < AA_FLOOR_045,
          f"at alpha=0.45 the AC floor {GAP_AB**2/C_AC(0.45):.6f} is weaker than AA's {AA_FLOOR_045}; AC's gain is the range [1/2, gamma)")

    # Kept Failure AC.1: above gamma the floor fails (Theorem K), live at 0.75 and 0.9.
    print()
    for alpha in (0.75, 0.90):
        Qs = []
        for n in range(1, 7):
            P, corners, N, w, V = levels[n]
            pos = w > 1e-10
            wa = np.where(pos, np.clip(w, 0, None) ** alpha, 0.0)
            A = (V * wa) @ V.T
            interior = np.setdiff1d(np.arange(N), corners)
            u = np.zeros(N)
            u[corners] = BC
            u[interior] = np.linalg.solve(A[np.ix_(interior, interior)], -A[np.ix_(interior, corners)] @ BC)
            Qs.append(float(u @ A @ u))
        Kb = [(3.5 ** alpha) * 3.0 ** (1 - alpha) * (3.0 * 5.0 ** (-alpha)) ** n for n in range(1, 7)]
        check(all(Qs[i + 1] < Qs[i] for i in range(5)) and all(q <= kb + 1e-12 for q, kb in zip(Qs, Kb)),
              f"alpha={alpha}: Q(u_n), n=1..6, strictly decreasing {[round(q, 6) for q in Qs]} and under the Theorem K bound -> 0")

    print()
    print("Theorem AC: Q_alpha(u_n) >= (9/4)/C_AC(alpha) > 0 at every level for every alpha < log3/log5.")
    print("Kept Failure AC.1: no level-uniform energy floor for alpha > log3/log5 (Theorem K).")
    print("alpha = log3/log5 itself is settled by Theorem AD (collapse); ||F(n,0.45)|| zero vs positive stays OPEN. Not thrust.")
    if failures:
        print(f"\n{len(failures)} check(s) failed")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
