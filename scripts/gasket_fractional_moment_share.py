#!/usr/bin/env python3
"""Exact share identity σ_n = (L^{1-α})_aa + ω_n and a corner-moment floor.

Witness for Theorem AP and Kept Failure AP.1 in GASKET.md.

Notation of Theorems AB–AK: corners a, b, c with data (1, -1/2, 0); u is the
L^α-Dirichlet minimizer (L^α u = 0 off the corners); I = (L u) on corners,
I^α = (L^α u) on corners; σ_n = I_a / I_a^α = ||F|| / ||F_frac||; β = 1 - α.

Theorem AP (every α ∈ (0,1), every level n ≥ 1):
  (1) exact share identity: L u = L^{1-α} (L^α u) (spectral calculus), and
      by Theorem AD (i) I^α = I_a^α (1, -4/5, -1/5), so
          σ_n = (L^β)_aa + ω_n,   ω_n := -(L^β)_ab = -(L^β)_ac > 0;
  (2) Radau moment floor: the spectral measure μ_a of e_a for L has moments
      (L^k)_aa = 1, 2, 6, 22, 92 (k = 0..4) at every level n ≥ 1, and
      λ^β has positive fifth derivative on (0,∞), so
          σ_n > (L^β)_aa ≥ S_AP(α) := (1 + 33^{-1/2}) t₁^{-α} + (1 - 33^{-1/2}) t₂^{-α},
      t₁,₂ = (13 ∓ √33)/4 (Gauss–Radau nodes of μ_a with fixed node 0);
  (3) Lobatto refinement (n ≥ 2): (L^5)_aa = 418 and spec(L) ⊂ [0,6]
      (Theorem AK (0)); λ^β has negative sixth derivative, so
          (L^β)_aa ≥ S_Lob(α) := Σ_i ν_i s_i^β + ν_6 6^β,
      s₁,₂ = (61 ∓ √685)/22, ν the {0, s₁, s₂, 6} law matching moments 0..5;
  (4) at α = 0.45: S_AP ≈ 1.310267, S_Lob ≈ 1.319222 (n ≥ 2), against
      Theorem AK's σ_AK ≈ 0.701556 and Kept Failure AK.1's cap ≈ 0.952863;
      with Theorem AM's fully proved ||F_frac|| ≥ 0.209685,
      ||F(n,0.45)|| ≥ S_AP · 0.209685 ≈ 0.274743 at every level, and with
      Theorem AO's floor as stated there (0.316911), ≥ 0.415238.
  Existence of the limit and Conjecture J.1's value L ≥ 0.70 are NOT proved.

Kept Failure AP.1: "a sharp enough share floor, times the proved fractional
floors (AM–AO), yields ||F(n,0.45)|| ≥ 0.70" is false: any level-uniform
share floor is ≤ σ_6 ≈ 1.352483, and 1.352483 · 0.316911 ≈ 0.4286 < 0.70.
The share side is now within ≈3.1% of the live σ_n; the remaining gap is the
fractional floor: with S_AP, J.1 would follow from R_α(n) ≤ 2.859 uniformly
(proved ceilings: ≥ 4.82; live R_α(6) ≈ 2.347). That is not proved.

Not thrust.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_cycle_ceiling import Ffrac_floor_cyc  # noqa: E402
from gasket_fractional_depth2_ceiling import Ffrac_floor_ao  # noqa: E402
from gasket_fractional_limit import (  # noqa: E402
    BC,
    build_gasket,
    combinatorial_laplacian,
    corner_flux_norm,
)
from gasket_fractional_moment_jensen import sigma_AK  # noqa: E402

ALPHA = 0.45
J1_VALUE = 0.70
AK1_CAP = 0.952863           # Kept Failure AK.1: live-constant chain σ at level 6
AO_FRAC_FLOOR_045 = 0.316911  # Theorem AO as stated
AM_FRAC_FLOOR_045 = 0.209685  # Theorem AM (Fourier, fully proved)
FF9, FF10 = 0.620670938843, 0.616380576325  # recorded ||F_frac(9)||, ||F_frac(10)||
SQRT21_3 = math.sqrt(21.0) / 3.0
MOM = (1.0, 2.0, 6.0, 22.0, 92.0)  # (L^k)_aa, k = 0..4, every n ≥ 1
M5 = 418.0                          # (L^5)_aa for n ≥ 2 (414 at n = 1)

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


# ------------------------------------------------------------- moment floors
R33 = math.sqrt(33.0)
T1, T2 = (13.0 - R33) / 4.0, (13.0 + R33) / 4.0
NU1, NU2 = (1.0 + 1.0 / R33) / T1, (1.0 - 1.0 / R33) / T2  # Radau weights at t₁, t₂
NU0 = 1.0 - NU1 - NU2                                       # = 3/17

R685 = math.sqrt(685.0)
S1, S2 = (61.0 - R685) / 22.0, (61.0 + R685) / 22.0


def lobatto_weights() -> np.ndarray:
    X = np.array([[1.0, 1.0, 1.0, 1.0],
                  [0.0, S1, S2, 6.0],
                  [0.0, S1 ** 2, S2 ** 2, 36.0],
                  [0.0, S1 ** 3, S2 ** 3, 216.0]])
    return np.linalg.solve(X, np.array(MOM[:4]))


LOB_W = lobatto_weights()
LOB_X = np.array([0.0, S1, S2, 6.0])


def S_AP(alpha: float) -> float:
    """Radau floor: (1 + 33^{-1/2}) t₁^{-α} + (1 - 33^{-1/2}) t₂^{-α}."""
    return (1.0 + 1.0 / R33) * T1 ** (-alpha) + (1.0 - 1.0 / R33) * T2 ** (-alpha)


def S_Lob(alpha: float) -> float:
    b = 1.0 - alpha
    return float(LOB_W[1] * S1 ** b + LOB_W[2] * S2 ** b + LOB_W[3] * 6.0 ** b)


def hermite_poly(beta: float, nodes: list[float]) -> np.poly1d:
    """Hermite interpolant of x^β: simple node 0 (and 6 if present), double t's."""
    rows, rhs = [], []
    deg = len(nodes) - 1
    seen: dict[float, int] = {}
    for x in nodes:
        k = seen.get(x, 0)
        seen[x] = k + 1
        if k == 0:
            rows.append([x ** j for j in range(deg + 1)])
            rhs.append(x ** beta if x > 0 else 0.0)
        else:
            rows.append([j * x ** (j - 1) if j else 0.0 for j in range(deg + 1)])
            rhs.append(beta * x ** (beta - 1.0))
    c = np.linalg.solve(np.array(rows), np.array(rhs))
    return np.poly1d(c[::-1])


def deriv_sign(beta: float, k: int) -> float:
    return math.prod(beta - i for i in range(k))


# ---------------------------------------------------------------- live data
def local_moments(L: np.ndarray, a: int, kmax: int = 5) -> list[float]:
    x = np.zeros(L.shape[0]); x[a] = 1.0
    out = [1.0]
    y = x.copy()
    for _ in range(kmax):
        y = L @ y
        out.append(float(y[a]))
    return out


def level_record(level: int, alphas) -> list[dict]:
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    lam, V = np.linalg.eigh(L)
    lam = np.clip(lam, 0.0, None)
    lam[lam < 1e-12] = 0.0
    a, b, c = (int(x) for x in corners)
    interior = np.setdiff1d(np.arange(N), corners)
    out = []
    for alpha in alphas:
        beta = 1.0 - alpha
        A = (V * lam ** alpha) @ V.T
        B = (V * lam ** beta) @ V.T
        u = np.zeros(N)
        u[corners] = BC
        u[interior] = np.linalg.solve(A[np.ix_(interior, interior)],
                                      -A[np.ix_(interior, corners)] @ BC)
        Lu, Au = L @ u, A @ u
        I, If = Lu[corners], Au[corners]
        out.append({
            "level": level, "alpha": alpha, "N": N,
            "lam_max": float(lam.max()),
            "mom": local_moments(L, a),
            "harm": float(np.abs(Au[interior]).max()),
            "LuBA": float(np.abs(Lu - B @ Au).max()),
            "ratio_b": float(If[1] / If[0]), "ratio_c": float(If[2] / If[0]),
            "Baa": float(B[a, a]), "omega_b": float(-B[a, b]), "omega_c": float(-B[a, c]),
            "Ia": float(I[0]), "Iaf": float(If[0]), "sigma": float(I[0] / If[0]),
            "nF": corner_flux_norm(P, corners, I),
            "nFf": corner_flux_norm(P, corners, If),
            "R": float((5.0 / 3.0) / If[0]),
            "spec_gap_radau": float(min(np.polyval(hermite_poly(beta, [0.0, T1, T1, T2, T2]), lam) * -1 + lam ** beta)),
            "spec_gap_lob": float(min(np.polyval(hermite_poly(beta, [0.0, S1, S1, S2, S2, 6.0]), lam) * -1 + lam ** beta)),
        })
    return out


def main(through: int = 6) -> int:
    assert through >= 4
    a = ALPHA
    print("Theorem AP: exact share identity σ_n = (L^{1-α})_aa + ω_n and a corner-moment floor")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0); Theorems J, AD (i), AK (0), AM, AO.")
    print("Not thrust.\n")
    fam, fao = Ffrac_floor_cyc(a), Ffrac_floor_ao(a)
    print(f"t1,t2 = {T1:.9f}, {T2:.9f};  Radau weights ν0,ν1,ν2 = {NU0:.9f}, {NU1:.9f}, {NU2:.9f}")
    print(f"s1,s2 = {S1:.9f}, {S2:.9f};  Lobatto weights = {LOB_W}")
    print(f"S_AP(0.45) = {S_AP(a):.9f};  S_Lob(0.45) = {S_Lob(a):.9f};  σ_AK = {sigma_AK(a):.9f};  AK.1 cap {AK1_CAP}")
    print(f"||F|| floors: AM route {S_AP(a)*fam:.9f}, AO route {S_AP(a)*fao:.9f}, AO+Lob {S_Lob(a)*fao:.9f}\n")

    # Closed forms.
    check(abs(NU0 - 3.0 / 17.0) < 1e-14, "Radau weight at 0 is 3/17")
    check(min(NU0, NU1, NU2) > 0.0, "Radau weights positive")
    for k, m in enumerate(MOM):
        val = NU1 * T1 ** k + NU2 * T2 ** k + (NU0 if k == 0 else 0.0)
        check(abs(val - m) < 1e-11, f"Radau law {{0,t1,t2}} has moment {k} = {m:g}")
    check(abs(T1 * T1 - 6.5 * T1 + 8.5) < 1e-12 and abs(T2 * T2 - 6.5 * T2 + 8.5) < 1e-12,
          "t1,t2 are roots of x² - 13x/2 + 17/2 (orthogonal for x dμ_a)")
    check(abs(11 * S1 ** 2 - 61 * S1 + 69) < 1e-11 and abs(11 * S2 ** 2 - 61 * S2 + 69) < 1e-11,
          "s1,s2 are roots of 11x² - 61x + 69 (orthogonal for x(6-x) dμ_a)")
    lob_mom = list(MOM) + [M5]
    for k, m in enumerate(lob_mom):
        val = float(LOB_W @ np.where(LOB_X == 0.0, 1.0 if k == 0 else 0.0, LOB_X ** k))
        check(abs(val - m) < 1e-9, f"Lobatto law {{0,s1,s2,6}} has moment {k} = {m:g}")
    check(min(LOB_W) > 0.0, "Lobatto weights positive")
    check(abs(S_AP(a) - (NU1 * T1 ** (1 - a) + NU2 * T2 ** (1 - a))) < 1e-14, "closed form of S_AP")
    # Derivative signs and majorant on a fine grid.
    for alpha in (0.01, 0.10, 0.30, 0.45, 0.60, 0.68, 0.80, 0.95, 0.99):
        beta = 1.0 - alpha
        check(deriv_sign(beta, 5) > 0.0 and deriv_sign(beta, 6) < 0.0,
              f"α={alpha}: d⁵/dx⁵ x^β > 0, d⁶/dx⁶ x^β < 0 on (0,∞)")
        gR = hermite_poly(beta, [0.0, T1, T1, T2, T2])
        gL = hermite_poly(beta, [0.0, S1, S1, S2, S2, 6.0])
        xs = np.concatenate([np.geomspace(1e-12, 1e-2, 400), np.linspace(1e-2, 60.0, 60001)])
        check(float((xs ** beta - gR(xs)).min()) >= -1e-10, f"α={alpha}: Radau quartic ≤ x^β on (0,60]")
        xs6 = xs[xs <= 6.0]
        check(float((xs6 ** beta - gL(xs6)).min()) >= -1e-10, f"α={alpha}: Lobatto quintic ≤ x^β on (0,6]")
    grid = np.linspace(0.01, 0.99, 99)
    check(all(S_Lob(x) > S_AP(x) > sigma_AK(x) for x in grid),
          "on α∈[0.01,0.99]: S_Lob > S_AP > σ_AK (strict improvement on Theorem AK)")
    check(S_AP(a) > 1.0 and S_AP(a) > AK1_CAP, "S_AP(0.45) > 1 and > Kept Failure AK.1's cap 0.952863")
    check(abs(S_AP(a) - 1.310267) < 1e-6 and abs(S_Lob(a) - 1.319222) < 1e-6, "S_AP ≈ 1.310267, S_Lob ≈ 1.319222")
    check(abs(fam - AM_FRAC_FLOOR_045) < 1e-6 and abs(fao - AO_FRAC_FLOOR_045) < 1e-6,
          "imported AM / AO fractional floors match 0.209685 / 0.316911")
    flAM, flAO = S_AP(a) * AM_FRAC_FLOOR_045, S_AP(a) * AO_FRAC_FLOOR_045
    check(abs(flAM - 0.274743) < 1e-6 and abs(flAO - 0.415238) < 1e-6 and flAM > 0.222331,
          "||F(n,0.45)|| floors 0.274743 (AM route) and 0.415238 (AO route) beat AO's 0.222331")

    # Live levels.
    others = (0.10, 0.30, 0.50, 0.60, 0.68, 0.80, 0.95)
    live045 = []
    for level in range(1, through + 1):
        alphas = (a,) + (others if level <= 5 else ())
        for d in level_record(level, alphas):
            al, lv = d["alpha"], d["level"]
            tag = f"α={al}, level {lv}"
            m = d["mom"]
            ok_m = all(abs(m[k] - MOM[k]) < 1e-9 for k in range(5)) and abs(m[5] - (414.0 if lv == 1 else M5)) < 1e-9
            ok = (ok_m and d["lam_max"] <= 6.0 + 1e-9 and d["harm"] < 1e-9 and d["LuBA"] < 1e-9
                  and abs(d["ratio_b"] + 0.8) < 1e-9 and abs(d["ratio_c"] + 0.2) < 1e-9
                  and abs(d["omega_b"] - d["omega_c"]) < 1e-10 and d["omega_b"] > 0.0
                  and abs(d["sigma"] - (d["Baa"] + d["omega_b"])) < 1e-9
                  and d["spec_gap_radau"] >= -1e-10
                  and d["sigma"] > d["Baa"] >= S_AP(al) - 1e-12
                  and d["nF"] >= S_AP(al) * d["nFf"] - 1e-12)
            if lv >= 2:
                ok = ok and d["spec_gap_lob"] >= -1e-10 and d["Baa"] >= S_Lob(al) - 1e-12
            else:
                ok = ok and d["sigma"] >= S_Lob(al) - 1e-12  # level 1: direct finite computation
            if al == a:
                ok = ok and d["nF"] >= flAO - 1e-12 and d["nFf"] >= AO_FRAC_FLOOR_045 - 1e-9
                live045.append(d)
            check(ok, f"{tag}: moments, Lu=L^βL^αu, I^α∝(1,-4/5,-1/5), σ=(L^β)_aa+ω, "
                      f"(L^β)_aa={d['Baa']:.6f} ≥ S_AP{'/S_Lob' if lv >= 2 else ''}, ω={d['omega_b']:.3e}")
    check(abs(live045[1]["nF"] - 1.165909) < 5e-6, f"level-2 lock ||F|| = {live045[1]['nF']:.6f}")
    d1 = live045[0]
    check(d1["Baa"] < S_Lob(a) < d1["sigma"],
          f"level 1 (m5 = 414): (L^β)_aa = {d1['Baa']:.6f} < S_Lob < σ_1 = {d1['sigma']:.6f} "
          "(Lobatto needs n ≥ 2; level 1 covered directly)")

    print("\nLive at α=0.45:")
    for d in live045:
        print(f"  level {d['level']}: σ={d['sigma']:.6f} = (L^β)_aa {d['Baa']:.6f} + ω {d['omega_b']:.3e};"
              f"  R_α={d['R']:.6f};  ||F_frac||={d['nFf']:.6f};  S_AP·||F_frac||={S_AP(a)*d['nFf']:.6f};  ||F||={d['nF']:.6f}")

    # Kept Failure AP.1.
    print("\nKept Failure AP.1: share floors cannot lift the proved fractional floors to 0.70")
    smin = min(d["sigma"] for d in live045)
    print(f"  min live σ_n (levels 1-{through}) = {smin:.6f};  × AO floor = {smin*AO_FRAC_FLOOR_045:.6f};"
          f"  × AM floor = {smin*AM_FRAC_FLOOR_045:.6f}")
    check(smin * AO_FRAC_FLOOR_045 < J1_VALUE and smin * AM_FRAC_FLOOR_045 < J1_VALUE,
          "any uniform share floor (≤ min live σ_n) × AO/AM floor < 0.70")
    check(S_Lob(a) * AO_FRAC_FLOOR_045 < J1_VALUE and 2.0 ** (1 - a) * AO_FRAC_FLOOR_045 < J1_VALUE,
          "S_Lob·0.316911 and 2^{1-α}·0.316911 < 0.70")
    need_F = J1_VALUE / S_AP(a)
    need_R = SQRT21_3 / need_F
    print(f"  with S_AP, J.1 would follow from ||F_frac|| ≥ {need_F:.6f}, i.e. R_α(n) ≤ {need_R:.6f}")
    check(abs(need_R - 2.859237) < 1e-5, "needed resistance ceiling ≈ 2.859237")
    check(need_R < 4.820047 and all(d["R"] < need_R for d in live045),
          "needed ceiling is below every proved ceiling (AO 4.820047) but above live R_α(n) (gap OPEN)")
    check(S_AP(a) * FF9 > J1_VALUE and S_AP(a) * FF10 > J1_VALUE,
          f"observation only: S_AP × recorded ||F_frac(9)||, ||F_frac(10)|| = {S_AP(a)*FF9:.6f}, {S_AP(a)*FF10:.6f} > 0.70")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        return 1
    print("All Theorem AP / Kept Failure AP.1 checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    through = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    raise SystemExit(main(through))
