#!/usr/bin/env python3
"""Two-moment (Hermite) Jensen sharpens the diagonal ceilings of Theorems AF–AJ.

Witness for Theorem AK and Kept Failure AK.1 in GASKET.md.

Notation of Theorems AE–AJ: a is the value-1 corner (deg a = 2), p, q its
two neighbors (deg 4, p ~ q), Far = V \\ {a, p, q}, w_xy = -(L^α)_{xy} > 0,
I_a = v_p + v_q with v = u_a - u ≥ 0, I_a^α = w_n I_a + I_far,
σ_n = I_a / I_a^α = ||F|| / ||F_frac||.

Theorem AK (every α ∈ (0,1), every level n ≥ 1):
  (0) line-graph identity: L = 3 diag(c) - NᵀN with N the (level-n cell
      triangle) × vertex incidence and c_x ∈ {1,2} the number of cells at x,
      so spec(L) ⊂ [0, 6] (attained at the checked levels 2-4);
  (1) two-moment ceiling: if μ is a probability measure on [0, 6] with mean
      m and variance s², then ∫ λ^α dμ ≤ π 6^α + (1-π) t^α with
      t = m - s²/(6-m), π = s²/(s² + (6-m)²)  (Hermite majorant: λ^α has
      positive third derivative). For e_a: m = 2, s² = 2, so
      D_a := (L^α)_aa ≤ D_a⁺(α) := 6^α/9 + (8/9)(3/2)^α < 2^α;
      for e_p: m = 4, s² = 4, so (L^α)_pp ≤ D_p⁺(α) := (6^α + 2^α)/2 < 4^α;
  (2) w_n = (D_a - W_far)/2 ≤ w_n⁺ := (D_a⁺ - W_far_min)/2 < 2^{α-1};
  (3) I_far ≤ ((D_p - w_pq)/2) I_a ≤ H_AK I_a, H_AK := (D_p⁺ - α4^{α-1})/2;
      γ_n ≤ Γ_AK := (D_p⁺ - α 4^{α-1}) / W_far_min;
  (4) σ_n ≥ σ_AK := 1/(w_n⁺ + H_AK), ρ_n ≥ ρ_AK := w_min/(w_min + H_AK);
  (5) for α < log3/log5: ||F(n,α)|| ≥ σ_AK · 9√21/(28 C_AC(α)) > 0;
      at α = 0.45 with AA's 0.084015: ||F(n,0.45)|| ≥ 0.058941.
  Existence of the limit and Conjecture J.1's value L ≥ 0.70 are NOT proved.

Kept Failure AK.1: "the one-step neighbor-Harnack chain (σ_n ≥ 1/(w_n +
(D_p - w_pq)/2)) reaches 0.70 once its constants are sharp" is false: with
every constant replaced by its LIVE value the chain gives σ ≈ 0.95286 < 1
at every computed level, and times the live ||F_frac(n)|| it falls below
0.70 from level 4 (≈0.6201 at level 6). The live σ_n ≈ 1.3525 > 1 is out
of reach of that family, so sharpening its constants cannot close J.1.

Not thrust.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_capacity_floor import C_AC  # noqa: E402
from gasket_fractional_lambda_ceiling import w_far_min  # noqa: E402
from gasket_fractional_lambda_floor import w_min  # noqa: E402
from gasket_fractional_limit import build_gasket, combinatorial_laplacian  # noqa: E402
from gasket_fractional_wpq_harnack import (  # noqa: E402
    harnack_const as harnack_const_AJ,
    level_data,
    sigma_AJ,
)

ALPHA = 0.45
GAMMA_C = math.log(3.0) / math.log(5.0)
AA_FRAC_FLOOR_045 = 0.084015
J1_VALUE = 0.70
LMAX = 6.0

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def hermite_two_point(m: float, var: float, M: float = LMAX) -> tuple[float, float]:
    """Atoms t and M (weight π on M) matching mean m and variance var."""
    t = m - var / (M - m)
    pi = var / (var + (M - m) ** 2)
    return t, pi


def moment_ceiling(alpha: float, m: float, var: float, M: float = LMAX) -> float:
    t, pi = hermite_two_point(m, var, M)
    return pi * M ** alpha + (1.0 - pi) * t ** alpha


def hermite_quadratic(alpha: float, t: float, M: float = LMAX):
    """Quadratic g with g(M)=M^α, g(t)=t^α, g'(t)=α t^{α-1}."""
    f = lambda x: x ** alpha  # noqa: E731
    d1 = alpha * t ** (alpha - 1.0)
    c2 = (f(M) - f(t) - d1 * (M - t)) / (M - t) ** 2
    return lambda x: f(t) + d1 * (x - t) + c2 * (x - t) ** 2


def Da_plus(alpha: float) -> float:
    return LMAX ** alpha / 9.0 + (8.0 / 9.0) * 1.5 ** alpha


def Dp_plus(alpha: float) -> float:
    return (LMAX ** alpha + 2.0 ** alpha) / 2.0


def wn_plus(alpha: float) -> float:
    return (Da_plus(alpha) - w_far_min(alpha)) / 2.0


def harnack_const(alpha: float) -> float:
    return (Dp_plus(alpha) - w_min(alpha)) / 2.0


def gamma_ceil(alpha: float) -> float:
    return (Dp_plus(alpha) - w_min(alpha)) / w_far_min(alpha)


def sigma_AK(alpha: float) -> float:
    return 1.0 / (wn_plus(alpha) + harnack_const(alpha))


def rho_AK(alpha: float) -> float:
    return w_min(alpha) / (w_min(alpha) + harnack_const(alpha))


def nF_floor_AC(alpha: float) -> float:
    return sigma_AK(alpha) * 9.0 * math.sqrt(21.0) / (28.0 * C_AC(alpha))


def cell_incidence(level: int):
    """Level-n cells as vertex triples: maximal triangles in the edge graph."""
    P, edges, corners = build_gasket(level)
    N = len(P)
    adj = [set() for _ in range(N)]
    for i, j in edges:
        adj[int(i)].add(int(j))
        adj[int(j)].add(int(i))
    tris = set()
    for i, j in edges:
        i, j = int(i), int(j)
        for k in adj[i] & adj[j]:
            tri = tuple(sorted((i, j, k)))
            ys = sorted(P[list(tri), 1])
            # Cells point up (two base vertices at the lower height); the
            # down-pointing 3-cycles are removed holes, not cells.
            if abs(ys[0] - ys[1]) < 1e-12 and ys[2] > ys[1] + 1e-12:
                tris.add(tri)
    Nmat = np.zeros((len(tris), N))
    for r, tri in enumerate(sorted(tris)):
        Nmat[r, list(tri)] = 1.0
    return P, edges, corners, Nmat


def structure_checks(level: int) -> None:
    P, edges, corners, Nmat = cell_incidence(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    c = Nmat.sum(axis=0)
    check(Nmat.shape[0] == 3 ** level, f"level {level}: {Nmat.shape[0]} cells = 3^n")
    check(set(c.tolist()) <= {1.0, 2.0} and int((c == 1).sum()) == 3,
          f"level {level}: every vertex in 2 cells except the 3 corners (1 cell)")
    check(np.allclose(np.diag(L), 2.0 * c), f"level {level}: deg = 2 c")
    check(np.array_equal(L, 3.0 * np.diag(c) - Nmat.T @ Nmat),
          f"level {level}: L = 3 diag(c) - NᵀN exactly")
    lam, V = np.linalg.eigh(L)
    check(lam.max() <= LMAX + 1e-9, f"level {level}: λ_max = {lam.max():.12f} ≤ 6")
    if level >= 2:
        check(abs(lam.max() - LMAX) < 1e-9, f"level {level}: λ_max = 6 (ceiling attained at this level)")
    a = int(corners[0])
    p = [j for j in range(N) if j != a and L[a, j] != 0.0][0]
    L2 = L @ L
    check(L[a, a] == 2.0 and L2[a, a] == 6.0, f"level {level}: (L)_aa = 2, (L²)_aa = 6")
    check(L[p, p] == 4.0 and L2[p, p] == 20.0, f"level {level}: (L)_pp = 4, (L²)_pp = 20")
    lam = np.clip(lam, 0.0, None)
    for alpha in (0.10, 0.45, 0.68, 0.90):
        Daa = float((lam ** alpha) @ (V[a] ** 2))
        Dpp = float((lam ** alpha) @ (V[p] ** 2))
        ok = Daa <= Da_plus(alpha) + 1e-12 and Dpp <= Dp_plus(alpha) + 1e-12
        check(ok, f"level {level}, α={alpha}: D_a={Daa:.6f} ≤ {Da_plus(alpha):.6f}, "
                  f"D_p={Dpp:.6f} ≤ {Dp_plus(alpha):.6f}")
        for x0, name in ((a, "a"), (p, "p")):
            m = float(L[x0, x0]); var = float(L2[x0, x0]) - m * m
            t, _ = hermite_two_point(m, var)
            g = hermite_quadratic(alpha, t)
            gap = g(lam) - lam ** alpha
            check(gap.min() >= -1e-12,
                  f"level {level}, α={alpha}: Hermite quadratic ≥ λ^α on spec(L) ({name})")
            Eg = float(g(lam) @ (V[x0] ** 2))
            check(abs(Eg - moment_ceiling(alpha, m, var)) < 1e-10,
                  f"level {level}, α={alpha}: E_{name} g = two-point value (moments only)")


def main(through: int = 6) -> int:
    assert through >= 4
    print("Theorem AK: two-moment Jensen on [0,6] sharpens the diagonal ceilings")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0); Theorems J, AG, AJ.")
    print("Not thrust.\n")

    a = ALPHA
    print(f"D_a+ = {Da_plus(a):.9f} (Jensen 2^a = {2**a:.9f});  D_p+ = {Dp_plus(a):.9f} (4^a = {4**a:.9f})")
    print(f"w_n+ = {wn_plus(a):.9f} (AF 2^(a-1) = {2**(a-1):.9f});  H_AK = {harnack_const(a):.9f} (H_AJ = {harnack_const_AJ(a):.9f})")
    print(f"Gamma_AK = {gamma_ceil(a):.6f};  sigma_AK = {sigma_AK(a):.9f} (AJ {sigma_AJ(a):.9f});  rho_AK = {rho_AK(a):.9f}")
    fl = sigma_AK(a) * AA_FRAC_FLOOR_045
    flAJ = sigma_AJ(a) * AA_FRAC_FLOOR_045
    print(f"||F(n,0.45)|| floor = {fl:.9f} (AJ {flAJ:.9f}; gain x{fl/flAJ:.4f});  AC route {nF_floor_AC(a):.6e}\n")

    # Closed forms and strict improvement on a grid of α.
    check(abs(moment_ceiling(a, 2.0, 2.0) - Da_plus(a)) < 1e-14, "corner two-point: t=3/2, π=1/9")
    check(abs(moment_ceiling(a, 4.0, 4.0) - Dp_plus(a)) < 1e-14, "neighbor two-point: t=2, π=1/2")
    check(hermite_two_point(2.0, 2.0) == (1.5, 1.0 / 9.0) and hermite_two_point(4.0, 4.0) == (2.0, 0.5),
          "two-point atoms exact")
    grid = np.linspace(0.01, 0.99, 99)
    check(all(Da_plus(x) < 2 ** x and Dp_plus(x) < 4 ** x and harnack_const(x) < harnack_const_AJ(x)
              and sigma_AK(x) > sigma_AJ(x) and wn_plus(x) < 2 ** (x - 1) and wn_plus(x) > w_min(x)
              for x in grid),
          "on α∈[0.01,0.99]: D_a+<2^α, D_p+<4^α, w_min<w_n+<2^(α-1), H_AK<H_AJ, σ_AK>σ_AJ")
    check(abs(fl - 0.058941) < 1e-6 and fl > flAJ, "||F(n,0.45)|| floor ≈ 0.058941 > AJ's 0.055599")
    check(0.0 < nF_floor_AC(a) < fl, "AC route floor positive, weaker than AA route at 0.45")

    for level in (1, 2, 3, 4):
        structure_checks(level)

    # Live chain checks.
    live = []
    for level in range(1, through + 1):
        d = level_data(level, a)
        Da = 2.0 * d["wn"] + d["Wfar"]
        live_sigma_chain = 1.0 / (d["wn"] + (d["diag_p"] - d["wpq"]) / 2.0)
        live.append((level, d, live_sigma_chain))
        check(Da <= Da_plus(a) + 1e-12 and d["diag_p"] <= Dp_plus(a) + 1e-12
              and d["diag_q"] <= Dp_plus(a) + 1e-12,
              f"level {level}: D_a={Da:.6f} ≤ D_a+, D_p=D_q={d['diag_p']:.6f} ≤ D_p+")
        check(d["wn"] <= wn_plus(a) + 1e-12 and d["Wfar"] >= w_far_min(a) - 1e-12,
              f"level {level}: w_n={d['wn']:.6f} ≤ w_n+, W_far={d['Wfar']:.6f} ≥ W_far_min")
        check(d["Ifar"] <= (d["diag_p"] - d["wpq"]) / 2.0 * d["Ia"] + 1e-12
              and d["Ifar"] <= harnack_const(a) * d["Ia"] + 1e-12,
              f"level {level}: I_far/I_a={d['Ifar']/d['Ia']:.6f} ≤ (D_p-w_pq)/2 ≤ H_AK")
        check(d["gamma"] <= gamma_ceil(a) + 1e-9 and d["sigma"] >= sigma_AK(a) - 1e-12
              and d["rho"] >= rho_AK(a) - 1e-12,
              f"level {level}: γ={d['gamma']:.4f} ≤ Γ_AK, σ={d['sigma']:.6f} ≥ σ_AK, ρ ≥ ρ_AK")
        check(d["nF"] >= sigma_AK(a) * d["nF_frac"] - 1e-12 and d["nF"] >= fl,
              f"level {level}: ||F||={d['nF']:.6f} ≥ σ_AK ||F_frac|| ≥ {fl:.6f}")
    for alpha in (0.10, 0.30, 0.50, 0.60, 0.68, 0.80, 0.95):
        for level in range(1, 6):
            d = level_data(level, alpha)
            ok = (2.0 * d["wn"] + d["Wfar"] <= Da_plus(alpha) + 1e-12
                  and d["diag_p"] <= Dp_plus(alpha) + 1e-12
                  and d["Ifar"] <= harnack_const(alpha) * d["Ia"] + 1e-12
                  and d["sigma"] >= sigma_AK(alpha) - 1e-12
                  and d["gamma"] <= gamma_ceil(alpha) + 1e-9)
            if alpha < GAMMA_C:
                ok = ok and d["nF"] >= nF_floor_AC(alpha) - 1e-12
            check(ok, f"α={alpha}, level {level}: D_a, D_p, Harnack, γ, σ"
                      + (" and AC-route floor" if alpha < GAMMA_C else "") + " hold")

    # Kept Failure AK.1.
    print("\nKept Failure AK.1: live-constant one-step neighbor-Harnack chain")
    for level, d, s in live:
        print(f"  level {level}: live σ_n={d['sigma']:.6f}, chain σ={s:.6f}, "
              f"chain·||F_frac||={s*d['nF_frac']:.6f}, ||F||={d['nF']:.6f}")
    check(all(s < 1.0 for _, _, s in live), "live-constant chain σ < 1 at every level")
    check(all(d["sigma"] > 1.3 for _, d, _ in live), "live σ_n > 1.3 at every level (out of the chain's reach)")
    check(all(s * d["nF_frac"] < J1_VALUE for lv, d, s in live if lv >= 4),
          "live-constant chain · live ||F_frac|| < 0.70 from level 4 (Kept Failure AK.1)")
    check(fl < J1_VALUE and 2.0 ** (1.0 - a) * AA_FRAC_FLOOR_045 < J1_VALUE,
          "AK floor and the Harnack-0 ceiling stay below 0.70")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        return 1
    print("All Theorem AK / Kept Failure AK.1 checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    through = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    raise SystemExit(main(through))
