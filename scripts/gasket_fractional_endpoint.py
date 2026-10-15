#!/usr/bin/env python3
"""The endpoint alpha = log3/log5: the fractional energy collapses like 1/n.

Witness for Theorem AD and Kept Failure AD.1 in GASKET.md.

Theorem AC floored Q_alpha(u_n) for every alpha < gamma = log3/log5 and left
alpha = gamma open. Theorem AD decides it:

  (1) exact symmetric Schur identity: the corner Schur complement of L^alpha
      is kappa*(3I - J) by the D3 symmetry of G_n, so for corner data
      (1,-1/2,0), Q_alpha(u_n) = (7/3) / R_alpha(n), R_alpha = v^T (L^alpha)^+ v,
      v = e_a - e_b (sharper than the duality constant 9/4 of AA/AC);
  (2) Stieltjes formula: lambda^{-a} = (sin(pi a)/pi) int_0^inf s^{-a}/(lambda+s) ds,
      so R_alpha = (sin(pi a)/pi) int s^{-a} R^s ds, R^s = v^T (L+s)^{-1} v;
  (3) cell test function: f = harmonic on the level-k corner cell at a with
      f(a)=1, f=0 at its other corners, f=0 elsewhere, has energy 2(3/5)^k and
      ||f||^2 <= N(k) <= 3^{k+1}, so R^s >= (5s)^{gamma-1}/5 for 5^{-n} <= s < 1;
  (4) hence R_gamma(n) >= c_AD * n, c_AD = (sin(pi gamma)/pi)(3/25) ln 5, and with
      AC's spectral mass bound R_gamma(n) <= 80((n+1) ln 3 + 1);
  (5) so (7/3)/(80((n+1)ln3+1)) <= Q_gamma(u_n) <= (7/3)/(c_AD n) -> 0, and by
      Theorem L ||F(n,gamma)|| <= kappa(gamma) sqrt(Q_gamma) -> 0.

Kept Failure AD.1: AC's floor does not extend to alpha = gamma.

Nothing here touches alpha = 0.45: zero vs positive for ||F(n,0.45)|| stays OPEN.
Not thrust.
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

GAMMA = ALPHA_C
C_AD = (math.sin(math.pi * GAMMA) / math.pi) * (3.0 / 25.0) * math.log(5.0)
SQRT21_7 = math.sqrt(21.0) / 7.0

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def kappa(alpha: float) -> float:
    wmin = alpha * 4.0 ** (alpha - 1.0)
    return (math.sqrt(21.0) / 5.0) * math.sqrt(2.0 / wmin)


def upper_R_gamma(n: int) -> float:
    return 80.0 * ((n + 1) * math.log(3.0) + 1.0)


def eig(level: int):
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    w, V = np.linalg.eigh(L)
    return P, corners, N, L, w, V


def frac(w, V, alpha):
    pos = w > 1e-10
    wa = np.where(pos, np.clip(w, 0, None) ** alpha, 0.0)
    return (V * wa) @ V.T


def dirichlet(A, corners, N):
    interior = np.setdiff1d(np.arange(N), corners)
    u = np.zeros(N)
    u[corners] = BC
    u[interior] = np.linalg.solve(A[np.ix_(interior, interior)], -A[np.ix_(interior, corners)] @ BC)
    S = A[np.ix_(corners, corners)] - A[np.ix_(corners, interior)] @ np.linalg.solve(
        A[np.ix_(interior, interior)], A[np.ix_(interior, corners)])
    return u, S


def main() -> int:
    print("Theorem AD: at alpha = log3/log5 the fractional energy and audit norm tend to 0 like 1/n")
    print(f"gamma = {GAMMA:.6f}, c_AD = {C_AD:.6f}. Assumption: spectral calculus (GASKET item 4).")
    print("Not thrust. Does not touch alpha = 0.45.\n")

    # (1) Symmetric Schur identity Q = (7/3)/R_alpha, several alpha, levels 1..6.
    alphas = (0.45, 0.60, GAMMA, 0.75, 0.90, 1.0)
    for n in range(1, 6):
        P, corners, N, L, w, V = eig(n)
        a, b = int(corners[0]), int(corners[1])
        c2 = (V[a, :] - V[b, :]) ** 2
        pos = w > 1e-10
        worst_q, worst_s = 0.0, 0.0
        for alpha in alphas:
            A = frac(w, V, alpha)
            u, S = dirichlet(A, corners, N)
            kap = -S[0, 1]
            worst_s = max(worst_s, float(np.max(np.abs(S - kap * (3 * np.eye(3) - np.ones((3, 3)))))) / kap)
            Q = float(u @ A @ u)
            R = float(np.sum(c2[pos] * w[pos] ** (-alpha)))
            worst_q = max(worst_q, abs(Q * R - 7.0 / 3.0))
        check(worst_s < 1e-8 and worst_q < 1e-8,
              f"level {n}: corner Schur complement = kappa(3I-J) (rel err {worst_s:.1e}) and Q*R_alpha = 7/3 (err {worst_q:.1e}) for alpha in {{0.45,0.6,gamma,0.75,0.9,1}}")

    # (3a) Cell test function energy 2(3/5)^k.
    for k in range(0, 6):
        P, corners, N, L, w, V = eig(k)
        u = np.zeros(N)
        bc = np.array([1.0, 0.0, 0.0])
        interior = np.setdiff1d(np.arange(N), corners)
        u[corners] = bc
        if len(interior):
            u[interior] = np.linalg.solve(L[np.ix_(interior, interior)], -L[np.ix_(interior, corners)] @ bc)
        E = float(u @ L @ u)
        check(abs(E - 2.0 * 0.6 ** k) < 1e-9 and u.min() >= -1e-12 and u.max() <= 1 + 1e-12,
              f"G_{k}: corner bump energy {E:.9f} = 2(3/5)^{k}, 0<=f<=1, N(k)={N}")

    # (3b) Massive resistance lower bound R^s >= (5s)^{gamma-1}/5 on 5^{-n} <= s < 1.
    data = {}
    for n in range(1, 8):
        P, corners, N, L, w, V = eig(n)
        data[n] = (P, corners, N, L, w, V)
        a, b = int(corners[0]), int(corners[1])
        c2 = (V[a, :] - V[b, :]) ** 2
        pos = w > 1e-10
        ss = np.logspace(-n * math.log10(5.0), -1e-9, 400)
        Rs = np.array([np.sum(c2[pos] / (w[pos] + s)) for s in ss])
        lb = (5.0 * ss) ** (GAMMA - 1.0) / 5.0
        ratio = float(np.min(Rs / lb))
        check(ratio >= 1.0, f"level {n}: min over s in [5^-{n},1) of R^s / ((5s)^(gamma-1)/5) = {ratio:.4f} >= 1")

    # (2) Stieltjes formula check at level 3 by log-grid quadrature.
    P, corners, N, L, w, V = data[3]
    a, b = int(corners[0]), int(corners[1])
    c2 = (V[a, :] - V[b, :]) ** 2
    pos = w > 1e-10
    t = np.linspace(-90.0, 45.0, 270001)
    s = np.exp(t)
    Rs = (c2[pos][None, :] / (w[pos][None, :] + s[:, None])).sum(axis=1)
    integ = np.trapezoid(s ** (1 - GAMMA) * Rs, t) * math.sin(math.pi * GAMMA) / math.pi
    Rg = float(np.sum(c2[pos] * w[pos] ** (-GAMMA)))
    check(abs(integ - Rg) < 1e-6, f"level 3: (sin pi g/pi) int s^-g R^s ds = {integ:.8f} = R_gamma = {Rg:.8f}")

    # (4)-(5) Endpoint sandwich, audit norm, levels 1..7.
    print()
    Rs_g, Qs_g, Fs_g, QH = [], [], [], []
    for n in range(1, 8):
        P, corners, N, L, w, V = data[n]
        a, b = int(corners[0]), int(corners[1])
        c2 = (V[a, :] - V[b, :]) ** 2
        pos = w > 1e-10
        R = float(np.sum(c2[pos] * w[pos] ** (-GAMMA)))
        A = frac(w, V, GAMMA)
        u, _ = dirichlet(A, corners, N)
        Q = float(u @ A @ u)
        I = (L @ u)[corners]
        nF = corner_flux_norm(P, corners, I)
        nFf = corner_flux_norm(P, corners, (A @ u)[corners])
        uH, _ = dirichlet(L, corners, N)
        QH.append(float(uH @ A @ uH))
        Rs_g.append(R); Qs_g.append(Q); Fs_g.append(nF)
        lo, hi = 7.0 / 3.0 / upper_R_gamma(n), 7.0 / 3.0 / (C_AD * n)
        check(C_AD * n <= R <= upper_R_gamma(n) and lo <= Q <= hi + 1e-12,
              f"level {n}: c_AD n = {C_AD*n:.4f} <= R_gamma = {R:.6f} <= {upper_R_gamma(n):.2f}; {lo:.6f} <= Q_gamma = {Q:.6f} <= {hi:.4f}")
        check(nF <= kappa(GAMMA) * math.sqrt(Q) and abs(nFf - Q * SQRT21_7) < 1e-9,
              f"level {n}: ||F(n,gamma)|| = {nF:.6f} <= kappa sqrt(Q) = {kappa(GAMMA)*math.sqrt(Q):.6f}; ||F_frac|| = Q sqrt21/7 = {nFf:.6f}")
    inc = [Rs_g[i + 1] - Rs_g[i] for i in range(6)]
    print(f"     R_gamma(n+1)-R_gamma(n), n=1..6: {[round(x, 6) for x in inc]}")
    check(all(x > 0 for x in inc) and all(Qs_g[i + 1] < Qs_g[i] for i in range(6)),
          f"R_gamma strictly increasing, Q_gamma strictly decreasing through level 7: {[round(q, 6) for q in Qs_g]}")

    # Kept Failure AD.1 and the Theorem K majorant at the endpoint.
    print()
    Kc = (3.5 ** GAMMA) * 3.0 ** (1 - GAMMA)
    check(abs(3.0 * 5.0 ** (-GAMMA) - 1.0) < 1e-12 and all(q <= Kc + 1e-12 for q in QH),
          f"Theorem K majorant at gamma is the constant {Kc:.6f} (3*5^-gamma = 1); live Q_gamma(u_H), n=1..7: {[round(q, 4) for q in QH]}")
    check(QH[-1] > Qs_g[-1] and QH[-1] > 0.5 * QH[0],
          "the harmonic comparison energy does not decay at gamma, so Theorem K cannot decide the endpoint; Theorem AD does")
    print("     AC floor (9/4)/C_AC(alpha) -> 0 as alpha -> gamma:",
          [round(2.25 / ((2 * GAMMA / (GAMMA - x)) * 40 ** (x / GAMMA)), 6) for x in (0.6, 0.65, 0.68, 0.6826)])

    print()
    print("Theorem AD: Q_gamma(u_n) = Theta(1/n) -> 0 and ||F(n,gamma)|| -> 0 at alpha = log3/log5.")
    print("Kept Failure AD.1: the AC energy floor does not extend to the endpoint alpha = gamma.")
    print("||F(n,0.45)|| zero vs positive stays OPEN. Not thrust.")
    if failures:
        print(f"\n{len(failures)} check(s) failed")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
