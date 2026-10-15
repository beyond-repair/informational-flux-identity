#!/usr/bin/env python3
"""Neighbor-edge weight floors the Harnack constant below AI's 4^α/2.

Witness for Theorem AJ and Kept Failure AJ.1 in GASKET.md.

Notation of Theorems AE–AI: a is the value-1 corner (deg a = 2), p, q its
two neighbors (deg 4, p ~ q), Far = V \\ {a, p, q}, w_xy = -(L^α)_{xy} > 0
for x != y, I_a = (L u)_a = v_p + v_q with v = u_a - u ≥ 0,
I_a^α = (L^α u)_a, I_far = Σ_{j∈Far} w_aj v_j, W_far = Σ_{j∈Far} w_aj,
γ_n = 2 I_far / (W_far I_a), σ_n = I_a / I_a^α = ||F|| / ||F_frac||,
w_pq = -(L^α)_{pq}.

Exact identity from L^α-harmonicity of v at p and q (kept cross terms):
    Σ_Far (w_pj + w_qj) v_j = (L^α)_pp v_p + (L^α)_qq v_q - w_pq (v_p + v_q).

Theorem AJ (every α ∈ (0,1), every level n ≥ 1):
  (1) neighbor-edge weight floor: w_pq ≥ w_min(α) := α 4^{α-1}
      (Theorem J on the edge p~q: H_t(p,q) ≥ t e^{-4t});
  (2) sharpened one-step Harnack: I_far ≤ ((4^α - w_pq)/2) I_a
      ≤ H_AJ(α) I_a, with H_AJ(α) := (4^α - α 4^{α-1})/2
      (AI's sum w_pj+w_qj ≥ 2 w_aj, Jensen (L^α)_pp,(L^α)_qq ≤ 4^α,
      and the exact identity above; strictly sharper than AI's 4^α/2);
  (3) γ_n ≤ (4^α - α 4^{α-1}) / W_far ≤ Γ_AJ(α) := 2(4-α)/(α(1-α));
  (4) σ_n ≥ σ_AJ(α) := 2/(2^α + 4^α - α 4^{α-1}),
      ρ_n ≥ ρ_AJ(α) := 2α/(α+4);
  (5) for α < log3/log5: ||F(n,α)|| ≥ σ_AJ · 9√21/(28 C_AC(α)) > 0 (AC);
      at α = 0.45 with AA's 0.084015: ||F(n,0.45)|| ≥ σ_AJ(0.45)·0.084015
      ≈ 0.055599 (about 1.069× AI's 0.051988).
  Existence of the limit and Conjecture J.1's value L ≥ 0.70 are NOT proved.

Kept Failure AJ.1: "keeping the w_pq cross term and combining with AA's floor
yields ||F(n,0.45)|| ≥ 0.70 (Conjecture J.1's value)" is false: the chain
gives σ_AJ · 0.084015 ≈ 0.0556, and even with the live ||F_frac(n)|| in
place of AA's floor, σ_AJ ||F_frac(n)|| < 0.70 at every computed level; the
Harnack-constant-0 ceiling 2^{1-α}·0.084015 ≈ 0.1230 < 0.70 still applies.

Not thrust.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_capacity_floor import C_AC  # noqa: E402
from gasket_fractional_degree_share import spectral_power  # noqa: E402
from gasket_fractional_harnack_sum import (  # noqa: E402
    gamma_ceil as gamma_ceil_AI,
    harnack_const as harnack_const_AI,
    sigma_AI,
)
from gasket_fractional_lambda_ceiling import w_far_min  # noqa: E402
from gasket_fractional_lambda_floor import w_min  # noqa: E402
from gasket_fractional_limit import (  # noqa: E402
    BC,
    build_gasket,
    combinatorial_laplacian,
    corner_flux_norm,
)

ALPHA = 0.45
GAMMA_C = math.log(3.0) / math.log(5.0)
AA_FRAC_FLOOR_045 = 0.084015  # Theorem AA: ||F_frac(n,0.45)|| >= 0.084015
J1_VALUE = 0.70               # Conjecture J.1 (not a theorem)
SQRT21_5 = math.sqrt(21.0) / 5.0

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def harnack_const(alpha: float) -> float:
    """Theorem AJ (2): I_far ≤ H_AJ I_a with H_AJ = (4^α - α 4^{α-1})/2."""
    return (4.0 ** alpha - alpha * 4.0 ** (alpha - 1.0)) / 2.0


def gamma_ceil(alpha: float) -> float:
    """Theorem AJ (3): Γ_AJ(α) = (4^α - α 4^{α-1}) / W_far_min = 2(4-α)/(α(1-α))."""
    return (4.0 ** alpha - alpha * 4.0 ** (alpha - 1.0)) / w_far_min(alpha)


def sigma_AJ(alpha: float) -> float:
    """Theorem AJ (4): σ_n ≥ 2/(2^α + 4^α - α 4^{α-1})."""
    return 1.0 / (2.0 ** (alpha - 1.0) + harnack_const(alpha))


def rho_AJ(alpha: float) -> float:
    """Theorem AJ (4): ρ_n ≥ 2α/(α+4)."""
    return w_min(alpha) / (w_min(alpha) + harnack_const(alpha))


def nF_floor_AC(alpha: float) -> float:
    """Theorem AJ (5) via AC: σ_AJ · 9√21/(28 C_AC)."""
    return sigma_AJ(alpha) * 9.0 * math.sqrt(21.0) / (28.0 * C_AC(alpha))


def level_data(level: int, alpha: float) -> dict:
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    A = spectral_power(L, alpha)
    a = int(corners[0])
    nbrs = [j for j in range(N) if j != a and L[a, j] != 0.0]
    assert len(nbrs) == 2
    p, q = nbrs
    far = np.array([j for j in range(N) if j not in (a, p, q)], dtype=int)
    interior = np.setdiff1d(np.arange(N), corners)
    u = np.zeros(N)
    u[corners] = BC
    u[interior] = np.linalg.solve(
        A[np.ix_(interior, interior)], -A[np.ix_(interior, corners)] @ BC
    )
    v = u[a] - u
    Au = A @ u
    W = -A
    wn = float(W[a, p])
    wpq = float(W[p, q])
    Ia = float((L @ u)[a])
    Iaf = float(Au[a])
    Ifar = float(np.dot(W[a, far], v[far]))
    Wfar = float(W[a, far].sum())
    LHS = float(np.dot(W[p, far] + W[q, far], v[far]))
    RHS = float(A[p, p] * v[p] + A[q, q] * v[q] - wpq * (v[p] + v[q]))
    sum_ratio = (W[p, far] + W[q, far]) / (2.0 * W[a, far])
    I = (L @ u)[corners]
    If = Au[corners]
    return {
        "N": N,
        "deg_a": float(L[a, a]),
        "deg_p": float(L[p, p]),
        "deg_q": float(L[q, q]),
        "p_adj_q": bool(L[p, q] == -1.0),
        "diag_p": float(A[p, p]),
        "diag_q": float(A[q, q]),
        "vmin": float(v.min()),
        "res_p": abs(float(Au[p])),
        "res_q": abs(float(Au[q])),
        "wn": wn,
        "wpq": wpq,
        "Ia": Ia,
        "vp_plus_vq": float(v[p] + v[q]),
        "Iaf": Iaf,
        "Ifar": Ifar,
        "Wfar": Wfar,
        "gamma": 2.0 * Ifar / (Wfar * Ia),
        "sigma": Ia / Iaf,
        "rho": wn * Ia / Iaf,
        "sum_ratio_min": float(sum_ratio.min()),
        "sum_slack_min": float(np.min(W[p, far] + W[q, far] - 2.0 * W[a, far])),
        "LHS": LHS,
        "RHS": RHS,
        "harn_live": (4.0 ** alpha - wpq) / 2.0,
        "nF": corner_flux_norm(P, corners, I),
        "nF_frac": corner_flux_norm(P, corners, If),
        "Iloc_sum": float(wn * Ia + Ifar),
    }


def heat_kernel_wpq(level: int, alpha: float) -> None:
    """Quadrature check: H_t(p,q) ≥ t e^{-4t} and w_pq ≥ α 4^{α-1}."""
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    lam, V = np.linalg.eigh(L)
    lam = np.clip(lam, 0.0, None)
    A = spectral_power(L, alpha)
    a = int(corners[0])
    nbrs = [j for j in range(N) if j != a and L[a, j] != 0.0]
    p, q = nbrs
    check(L[p, q] == -1.0, f"level {level}: p ~ q (L_pq = -1)")
    M = 4.0 * np.eye(N) - L
    check(M[p, q] == 1.0 and bool(np.all(M >= -1e-15)),
          f"level {level}: M = 4I-L >= 0 entrywise and M_pq = 1")

    # H_t(p,q) >= t e^{-4t} on a grid of t (series for small t, eigh beyond).
    worst_ratio = 1.0
    for t in (0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0):
        Ht = float((np.exp(-t * lam) * V[p]) @ V[q])
        lower = t * math.exp(-4.0 * t)
        worst_ratio = min(worst_ratio, Ht / lower if lower > 0 else 1.0)
        check(Ht + 1e-15 >= lower,
              f"level {level}: H_{t}(p,q)={Ht:.6e} >= t e^{{-4t}}={lower:.6e}")
    check(worst_ratio >= 1.0 - 1e-12,
          f"level {level}: min H_t(p,q)/(t e^{{-4t}}) = {worst_ratio:.6f} >= 1")

    # Analytic integral of the pointwise lower bound recovers alpha 4^(alpha-1):
    # c_a ∫_0^∞ (t e^{-4t}) t^{-1-a} dt = (a/Gamma(1-a)) * Gamma(1-a) / 4^{1-a}.
    c_alpha = alpha / math.gamma(1.0 - alpha)
    bound_analytic = c_alpha * math.gamma(1.0 - alpha) / 4.0 ** (1.0 - alpha)
    bound = alpha * 4.0 ** (alpha - 1.0)
    check(abs(bound_analytic - bound) < 1e-14,
          f"level {level}: ∫ t e^{{-4t}} t^{{-1-a}} recovers alpha 4^(a-1) ({bound:.8f})")
    wpq = float(-A[p, q])
    check(wpq >= bound - 1e-12,
          f"level {level}: w_pq={wpq:.8f} >= alpha 4^(alpha-1)={bound:.8f}")


def main(through: int = 6) -> int:
    assert through >= 4
    print("Theorem AJ: neighbor-edge weight sharpens the Harnack constant below AI's 4^a/2")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0); Theorems J, AI.")
    print("Not thrust.\n")

    hc = harnack_const(ALPHA)
    hc_AI = harnack_const_AI(ALPHA)
    Gc = gamma_ceil(ALPHA)
    sm = sigma_AJ(ALPHA)
    sm_AI = sigma_AI(ALPHA)
    rm = rho_AJ(ALPHA)
    fl = sm * AA_FRAC_FLOOR_045
    fl_AI = sm_AI * AA_FRAC_FLOOR_045
    flAC = nF_floor_AC(ALPHA)
    wm = w_min(ALPHA)
    print(f"w_min = alpha 4^(a-1) = {wm:.9f}")
    print(f"H_AJ = (4^a - alpha 4^(a-1))/2 = {hc:.9f}  (AI: 4^a/2 = {hc_AI:.9f})")
    print(f"Gamma_AJ(0.45) = {Gc:.6f}  (Gamma_AI = {gamma_ceil_AI(ALPHA):.6f})")
    print(f"sigma_AJ = {sm:.9f};  rho_AJ = {rm:.9f}  (sigma_AI = {sm_AI:.9f})")
    print(f"||F|| floor (AA) = sigma_AJ * 0.084015 = {fl:.9f}  (AI: {fl_AI:.9f}; gain x{fl/fl_AI:.4f})")
    print(f"||F|| floor (AC) = {flAC:.9e}\n")

    check(abs(wm - ALPHA * 4.0 ** (ALPHA - 1.0)) < 1e-15, "w_min = alpha 4^(a-1)")
    check(abs(hc - (4.0 ** ALPHA - wm) / 2.0) < 1e-15 and hc < hc_AI,
          "H_AJ = (4^a - w_min)/2 < AI's 4^a/2")
    check(abs(sm - 2.0 / (2.0 ** ALPHA + 4.0 ** ALPHA - ALPHA * 4.0 ** (ALPHA - 1.0))) < 1e-14,
          "sigma_AJ closed form 2/(2^a + 4^a - a 4^(a-1))")
    check(abs(sm - 0.6617753) < 1e-7, "sigma_AJ(0.45) ~ 0.6617753")
    check(abs(fl - 0.055599) < 1e-6 and fl > fl_AI, "||F(n,0.45)|| floor ~ 0.055599 > AI's 0.051988")
    check(abs(Gc - 2.0 * (4.0 - ALPHA) / (ALPHA * (1.0 - ALPHA))) < 1e-12 and Gc < gamma_ceil_AI(ALPHA),
          "Gamma_AJ = 2(4-a)/(a(1-a)) < Gamma_AI")
    check(abs(rm - 2.0 * ALPHA / (ALPHA + 4.0)) < 1e-14, "rho_AJ = 2a/(a+4)")
    check(0.0 < flAC < fl, "AC route floor positive and weaker than AA route at 0.45")

    grid0 = np.linspace(0.01, 0.99, 99)
    check(all(abs(gamma_ceil(a_) - 2.0 * (4.0 - a_) / (a_ * (1.0 - a_))) < 1e-9 * gamma_ceil(a_)
              and abs(rho_AJ(a_) - 2.0 * a_ / (a_ + 4.0)) < 1e-14
              and abs(sigma_AJ(a_) - 2.0 / (2.0 ** a_ + 4.0 ** a_ - a_ * 4.0 ** (a_ - 1.0))) < 1e-12
              and abs(gamma_ceil(a_) - gamma_ceil_AI(a_) * (1.0 - a_ / 4.0)) < 1e-9 * gamma_ceil(a_)
              and abs(harnack_const(a_) - harnack_const_AI(a_) * (1.0 - a_ / 4.0)) < 1e-12
              for a_ in grid0),
          "closed forms: Gamma_AJ = 2(4-a)/(a(1-a)) = Gamma_AI*(1-a/4), "
          "rho_AJ = 2a/(a+4), H_AJ = (4^a/2)*(1-a/4)")
    check(all(sigma_AJ(a_) > sigma_AI(a_) and gamma_ceil(a_) < gamma_ceil_AI(a_)
              and harnack_const(a_) < harnack_const_AI(a_) for a_ in grid0),
          "sigma_AJ > sigma_AI, Gamma_AJ < Gamma_AI, H_AJ < H_AI on a grid of a in (0,1)")
    check(all(nF_floor_AC(a_) > 0.0 for a_ in grid0 if a_ < GAMMA_C),
          "||F|| floor via AC positive on the grid below log3/log5")

    # Heat-kernel lower bound for the edge p~q.
    for n in (1, 2, 3):
        heat_kernel_wpq(n, ALPHA)
    heat_kernel_wpq(2, 0.10)
    heat_kernel_wpq(2, 0.90)

    # Live table at alpha = 0.45.
    rows = []
    for n in range(1, through + 1):
        r = level_data(n, ALPHA)
        rows.append(r)
        print(
            f"level {n} N={r['N']} w_pq={r['wpq']:.8f} H_live={(4**ALPHA - r['wpq'])/2:.8f} "
            f"Ifar/Ia={r['Ifar']/r['Ia']:.8f} gamma={r['gamma']:.8f} "
            f"sigma={r['sigma']:.8f} ||F||={r['nF']:.8f} ||F_frac||={r['nF_frac']:.8f}",
            flush=True,
        )

    for n, r in enumerate(rows, start=1):
        check(r["deg_a"] == 2.0 and r["deg_p"] == 4.0 and r["deg_q"] == 4.0 and r["p_adj_q"],
              f"level {n}: deg a = 2, deg p = deg q = 4, p ~ q")
        check(abs(r["vp_plus_vq"] - r["Ia"]) < 1e-12,
              f"level {n}: v_p + v_q = I_a (combinatorial)")
        check(abs(r["LHS"] - r["RHS"]) < 1e-10,
              f"level {n}: exact identity LHS = (L^a)_pp v_p + (L^a)_qq v_q - w_pq I_a")
        check(r["sum_slack_min"] >= 0.0 and r["sum_ratio_min"] >= 1.0,
              f"level {n}: w_pj + w_qj >= 2 w_aj for every far j (min ratio {r['sum_ratio_min']:.8f})")
        check(r["wpq"] >= wm - 1e-12,
              f"level {n}: w_pq={r['wpq']:.8f} >= w_min={wm:.8f}")
        check(r["diag_p"] <= 4.0 ** ALPHA + 1e-12 and r["diag_q"] <= 4.0 ** ALPHA + 1e-12,
              f"level {n}: (L^a)_pp={r['diag_p']:.8f} <= 4^a (Jensen)")
        check(r["vmin"] >= -1e-12, f"level {n}: v = u_a - u >= 0 (max principle)")
        check(r["res_p"] < 1e-10 and r["res_q"] < 1e-10, f"level {n}: v L^a-harmonic at p, q")
        # Live sharpened bound with actual w_pq, then uniform H_AJ.
        check(r["Ifar"] <= r["harn_live"] * r["Ia"] + 1e-12,
              f"level {n}: I_far <= ((4^a - w_pq)/2) I_a "
              f"({r['Ifar']:.8f} <= {r['harn_live']*r['Ia']:.8f})")
        check(r["Ifar"] <= hc * r["Ia"] + 1e-12,
              f"level {n}: I_far={r['Ifar']:.8f} <= H_AJ I_a={hc*r['Ia']:.8f}")
        check(r["Ifar"] <= hc_AI * r["Ia"] + 1e-12,
              f"level {n}: I_far still <= AI's (4^a/2) I_a (consistency)")
        check(abs(r["Iloc_sum"] - r["Iaf"]) < 1e-10, f"level {n}: I_a^a = w_n I_a + I_far")
        check(r["gamma"] <= (4.0 ** ALPHA - wm) / r["Wfar"] + 1e-9 and r["gamma"] <= Gc,
              f"level {n}: gamma={r['gamma']:.8f} <= (4^a-w_min)/W_far "
              f"={(4**ALPHA-wm)/r['Wfar']:.4f} <= Gamma_AJ={Gc:.4f}")
        check(r["sigma"] >= sm and r["rho"] >= rm,
              f"level {n}: sigma={r['sigma']:.6f} >= sigma_AJ; rho={r['rho']:.6f} >= rho_AJ")
        check(abs(r["nF"] - SQRT21_5 * r["Ia"]) < 1e-10 and abs(r["nF_frac"] - SQRT21_5 * r["Iaf"]) < 1e-10,
              f"level {n}: ||F|| = (sqrt21/5) I_a, ||F_frac|| = (sqrt21/5) I_a^a")
        check(r["nF_frac"] >= AA_FRAC_FLOOR_045 and r["nF"] >= sm * r["nF_frac"] >= fl,
              f"level {n}: ||F||={r['nF']:.8f} >= sigma_AJ ||F_frac|| >= {fl:.6f}")

    gam_rec = (1.443835, 1.390438, 1.378096, 1.375485, 1.374963, 1.374864)
    for i, r in enumerate(rows[:6]):
        check(abs(r["gamma"] - gam_rec[i]) < 5e-7, f"level {i+1}: gamma matches recorded AE table")
    check(abs(rows[1]["nF"] - 1.165909) < 5e-6, "level 2 ||F|| lock 1.165909")

    # Other alpha in (0,1), through level 5.
    for a_ in (0.10, 0.30, 0.50, 0.60, 0.68, 0.80, 0.95):
        for n in range(1, 6):
            r = level_data(n, a_)
            ok = (
                r["wpq"] >= w_min(a_) - 1e-12
                and abs(r["LHS"] - r["RHS"]) < 1e-9
                and abs(r["vp_plus_vq"] - r["Ia"]) < 1e-12
                and r["sum_slack_min"] >= 0.0
                and r["Ifar"] <= r["harn_live"] * r["Ia"] + 1e-12
                and r["Ifar"] <= harnack_const(a_) * r["Ia"] + 1e-12
                and r["gamma"] <= gamma_ceil(a_)
                and r["sigma"] >= sigma_AJ(a_)
                and r["rho"] >= rho_AJ(a_)
                and r["vmin"] >= -1e-12
            )
            if a_ < GAMMA_C:
                ok = ok and r["nF"] >= nF_floor_AC(a_)
            check(ok, f"a={a_:.2f} level {n}: w_pq={r['wpq']:.6f} >= {w_min(a_):.6f}; "
                      f"gamma={r['gamma']:.6f} <= Gamma_AJ={gamma_ceil(a_):.3f}; "
                      f"sigma={r['sigma']:.6f} >= {sigma_AJ(a_):.6f}"
                      + (f"; ||F||={r['nF']:.6f} >= {nF_floor_AC(a_):.3e}" if a_ < GAMMA_C else ""))

    # Kept Failure AJ.1: w_pq cross term plus AA does not reach Conjecture J.1's 0.70.
    print(f"\nKept Failure AJ.1: sigma_AJ * 0.084015 = {fl:.6f} vs 0.70")
    chain_live = [sm * r["nF_frac"] for r in rows]
    print("sigma_AJ * live ||F_frac(n)||:", [f"{x:.6f}" for x in chain_live])
    check(fl < J1_VALUE, "Kept Failure AJ.1: sigma_AJ * 0.084015 < 0.70 (chain does not give J.1)")
    check(all(x < J1_VALUE for x in chain_live),
          "Kept Failure AJ.1: even sigma_AJ * live ||F_frac(n)|| < 0.70 at every computed level")
    check(sm < 1.0 and AA_FRAC_FLOOR_045 < J1_VALUE, "Kept Failure AJ.1: sigma_AJ < 1 and AA floor < 0.70")
    perfect = AA_FRAC_FLOOR_045 / 2.0 ** (ALPHA - 1.0)
    check(perfect < J1_VALUE,
          f"Kept Failure AJ.1: even Harnack constant 0 gives only 2^(1-a) * 0.084015 = {perfect:.6f} < 0.70")

    # Observations (not theorems).
    wpqs = [r["wpq"] for r in rows]
    print("live w_pq:", [f"{x:.6f}" for x in wpqs])
    check(all(wpqs[i + 1] < wpqs[i] for i in range(len(wpqs) - 1)) and wpqs[-1] > wm,
          "observation: w_pq decreases toward a positive limit above w_min (uniform floor does not vanish)")
    check(all(r["Ifar"] / r["Ia"] < 0.55 * hc for r in rows),
          "observation: live I_far/I_a < 0.55 H_AJ (Harnack constant still not sharp)")
    check(all(r["nF"] > 0.70 for r in rows), "observation: ||F|| > 0.70 through the range (J.1 not proved)")
    # Live bound with actual w_pq is sharper than H_AJ but still far from live Ifar/Ia.
    live_h = [(4.0 ** ALPHA - r["wpq"]) / 2.0 for r in rows]
    print("live ((4^a - w_pq)/2):", [f"{x:.6f}" for x in live_h])
    check(all(live_h[i] < hc_AI and live_h[i] > rows[i]["Ifar"] / rows[i]["Ia"] for i in range(len(rows))),
          "observation: live ((4^a-w_pq)/2) sits between Ifar/Ia and AI's 4^a/2")

    print()
    print("Theorem AJ: w_pq >= alpha 4^(a-1); I_far <= ((4^a - w_pq)/2) I_a <= H_AJ I_a;")
    print("gamma_n <= 2(4-a)/(a(1-a)); sigma_n >= 2/(2^a + 4^a - a 4^(a-1));")
    print("for a < log3/log5, ||F(n,a)|| >= sigma_AJ ||F_frac|| > 0.")
    print(f"At a=0.45: ||F(n,0.45)|| >= {fl:.6f} at every level (AI: 0.051988).")
    print("Kept Failure AJ.1: w_pq cross term with AA's floor gives ~0.0556, not Conjecture J.1's 0.70.")
    print("Existence of the limit and Conjecture J.1 (L >= 0.70) stay OPEN. Not thrust.")
    if failures:
        print(f"\n{len(failures)} check(s) failed")
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    through = 6
    if len(sys.argv) > 1:
        through = int(sys.argv[1])
    raise SystemExit(main(through))
