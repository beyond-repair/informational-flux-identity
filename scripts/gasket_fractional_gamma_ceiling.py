#!/usr/bin/env python3
"""Uniform ceiling on the far-to-neighbor gap ratio; positive audit floor.

Witness for Theorem AH and Kept Failure AH.1 in GASKET.md.

Notation of Theorems AE–AG: a is the value-1 corner, p, q its two
neighbors, Far = V \\ {a, p, q} (every vertex at graph distance >= 2 from a,
including the corners b, c), w_xy = -(L^α)_{xy} > 0 for x != y,
I_a = (L u)_a, I_a^α = (L^α u)_a, I_far = Σ_{j∈Far} w_aj (u_a - u_j),
W_far = Σ_{j∈Far} w_aj, γ_n = ḡ_far/ḡ_nbr = 2 I_far / (W_far I_a),
σ_n = I_a / I_a^α = ||F|| / ||F_frac||.

Theorem AH (every α ∈ (0,1), every level n ≥ 1, fixed s = 1/10):
  (1) kernel comparison: for every j ∈ Far and x ∈ {p, q},
        w_xj ≥ κ(α) w_aj,
        κ(α) = s e^{-4s} 2^{-1-α} / (1 + K(α)),
        K(α) = e^{8s} α / ((2-α)(1-2^{-α}));
      (semigroup H_{τ+s}(p,j) ≥ H_s(p,a) H_τ(a,j), H_s(p,a) ≥ s e^{-4s},
      and a small-time bound ∫_0^s ≤ K ∫_s^{2s} from (M^k)_{aj} = 0, k ≤ 1);
  (2) one-step Harnack at the neighbors: v = u_a - u ≥ 0 (max principle) is
      L^α-harmonic at p, q and vanishes at a, so with (L^α)_pp ≤ 4^α (Jensen,
      deg p = 4),  I_far ≤ (4^α / (2κ)) I_a;
  (3) uniform gap-ratio ceiling  γ_n ≤ 4^α/(κ W_far) ≤ Γ_AH(α) := 4^α/(κ W_far_min)
      (W_far_min from AG);
  (4) share floors  σ_n ≥ σ_min = 1/(2^{α-1} + 4^α/(2κ)),
                    ρ_n ≥ ρ_min = w_min/(w_min + 4^α/(2κ));
  (5) for α ∈ (0, log3/log5): ||F(n,α)|| ≥ σ_min ||F_frac|| ≥
      σ_min · 9√21/(28 C_AC(α)) > 0 (AC); at α = 0.45 with AA's floor 0.084015,
      ||F(n,0.45)|| ≥ σ_min(0.45) · 0.084015 ≈ 6.44e-4 at every level.
  So ||F(n,0.45)|| does NOT tend to 0. Existence of the limit and Conjecture
  J.1's value L ≥ 0.70 are NOT proved.

Kept Failure AH.1: the shortcut "Theorem J's gap cap u_a - u_j ≤ 3/2 alone
floors I_a, via I_a ≥ (I_a^α - (3/2) W_far)/w_n > 0" is false: live
I_a^α < (3/2) W_far at levels 3–6, so that lower bound is negative there.
A Harnack-type comparison (Theorem AH (1)–(2)) is needed.

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
S = 0.1
AA_FRAC_FLOOR_045 = 0.084015  # Theorem AA: ||F_frac(n,0.45)|| >= 0.084015
SQRT21_5 = math.sqrt(21.0) / 5.0

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def K_small(alpha: float, s: float = S) -> float:
    """Small-time constant: ∫_0^s H τ^{-1-α} ≤ K ∫_s^{2s} H τ^{-1-α}."""
    return math.exp(8.0 * s) * alpha / ((2.0 - alpha) * (1.0 - 2.0 ** (-alpha)))


def kappa(alpha: float, s: float = S) -> float:
    """Theorem AH (1): w_pj ≥ κ w_aj for every j ∈ Far."""
    return s * math.exp(-4.0 * s) * 2.0 ** (-1.0 - alpha) / (1.0 + K_small(alpha, s))


def harnack_const(alpha: float) -> float:
    """Theorem AH (2): I_far ≤ (4^α/(2κ)) I_a."""
    return 4.0 ** alpha / (2.0 * kappa(alpha))


def gamma_ceil(alpha: float) -> float:
    """Theorem AH (3): Γ_AH(α) = 4^α / (κ W_far_min)."""
    return 4.0 ** alpha / (kappa(alpha) * w_far_min(alpha))


def sigma_min(alpha: float) -> float:
    """Theorem AH (4): σ_n ≥ 1/(2^{α-1} + 4^α/(2κ))."""
    return 1.0 / (2.0 ** (alpha - 1.0) + harnack_const(alpha))


def rho_min(alpha: float) -> float:
    wm = w_min(alpha)
    return wm / (wm + harnack_const(alpha))


def nF_floor_AC(alpha: float) -> float:
    """Theorem AH (5) via AC: σ_min · 9√21/(28 C_AC)."""
    return sigma_min(alpha) * 9.0 * math.sqrt(21.0) / (28.0 * C_AC(alpha))


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
    Ia = float((L @ u)[a])
    Iaf = float(Au[a])
    Ifar = float(np.dot(W[a, far], v[far]))
    Wfar = float(W[a, far].sum())
    ratio_p = float(np.min(W[p, far] / W[a, far]))
    ratio_q = float(np.min(W[q, far] / W[a, far]))
    I = (L @ u)[corners]
    If = Au[corners]
    return {
        "N": N,
        "a": a,
        "p": p,
        "q": q,
        "far_dist2": bool(np.all(L[a, far] == 0.0)),
        "deg_p": float(L[p, p]),
        "deg_q": float(L[q, q]),
        "diag_p": float(A[p, p]),
        "diag_q": float(A[q, q]),
        "vmin": float(v.min()),
        "res_p": abs(float(Au[p])),
        "res_q": abs(float(Au[q])),
        "wn": wn,
        "Ia": Ia,
        "Iaf": Iaf,
        "Ifar": Ifar,
        "Wfar": Wfar,
        "gamma": 2.0 * Ifar / (Wfar * Ia),
        "sigma": Ia / Iaf,
        "rho": wn * Ia / Iaf,
        "ratio_p": ratio_p,
        "ratio_q": ratio_q,
        "nF": corner_flux_norm(P, corners, I),
        "nF_frac": corner_flux_norm(P, corners, If),
        "Iloc_sum": float(wn * Ia + Ifar),
    }


def heat_kernel_checks(level: int, alpha: float) -> None:
    """Quadrature check of the three heat-kernel steps behind κ."""
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    lam, V = np.linalg.eigh(L)
    lam = np.clip(lam, 0.0, None)
    A = spectral_power(L, alpha)
    a = int(corners[0])
    nbrs = [j for j in range(N) if j != a and L[a, j] != 0.0]
    p = nbrs[0]
    far = [j for j in range(N) if j not in (a, *nbrs)]
    c_alpha = alpha / math.gamma(1.0 - alpha)

    # H_t = e^{-4t} e^{tM}, M = 4I - L >= 0 entrywise: for t <= T_SER use the
    # nonnegative series (no cancellation for tiny far entries), else eigh.
    M = 4.0 * np.eye(N) - L
    KMAX = 160
    T_SER = 4.0
    pows = {}
    for x0 in (a, p):
        rows_k = np.zeros((KMAX + 1, N))
        rows_k[0, x0] = 1.0
        for kk in range(1, KMAX + 1):
            rows_k[kk] = M @ rows_k[kk - 1]
        pows[x0] = rows_k
    ks = np.arange(KMAX + 1)
    logfact = np.array([math.lgamma(kk + 1.0) for kk in ks])

    def H(t: np.ndarray, x: int) -> np.ndarray:
        # rows: times, cols: vertices; H_t(x, ·)
        t = np.asarray(t, dtype=np.float64)
        out = np.empty((t.size, N))
        small = t <= T_SER
        if np.any(small):
            ts = t[small]
            coef = np.exp(ks[None, :] * np.log(ts)[:, None] - logfact[None, :] - 4.0 * ts[:, None])
            out[small] = coef @ pows[x]
        if np.any(~small):
            out[~small] = (np.exp(-np.outer(t[~small], lam)) * V[x]) @ V.T
        return out

    y = np.linspace(0.0, 60.0, 30001)
    t_lo = S * np.exp(-y)        # (0, s]
    t_hi = S * np.exp(y)         # [s, ∞)
    Ha_lo, Ha_hi = H(t_lo, a), H(t_hi, a)

    def integ(vals: np.ndarray, t: np.ndarray, pw: float) -> np.ndarray:
        # ∫ f(t) t^{-1-pw} dt with t = s e^{±y}: dt/t = dy
        return np.trapezoid(vals * (t ** (-pw))[:, None], y, axis=0)

    I_lo = integ(Ha_lo, t_lo, alpha)
    # analytic tail beyond T_end, where H_t(a,·) = 1/N up to e^{-λ_1 T_end}
    T_end = float(t_hi[-1])
    tail = (1.0 / N) * T_end ** (-alpha) / alpha
    I_hi = integ(Ha_hi, t_hi, alpha) + tail
    sub_a = c_alpha * (I_lo + I_hi)
    # subordination identity reproduces the off-diagonal weights
    err = max(abs(sub_a[j] + A[a, j]) / abs(A[a, j]) for j in far)
    check(err < 1e-5, f"level {level}: subordination quadrature reproduces w_aj (rel err {err:.1e})")
    # small-time bound ∫_0^s ≤ K ∫_s^{2s}
    mask = t_hi <= 2.0 * S
    I_s2s = np.trapezoid(Ha_hi[mask] * (t_hi[mask] ** (-alpha))[:, None], y[mask], axis=0)
    K = K_small(alpha)
    worst = max(I_lo[j] / I_s2s[j] for j in far)
    check(worst <= K, f"level {level}: max_j ∫_0^s / ∫_s^2s = {worst:.6f} <= K = {K:.6f}")
    # pointwise lemma H_τ ≤ (τ/s)^2 e^{8s} H_τ' for τ ≤ s ≤ τ' ≤ 2s
    taus = S * np.linspace(0.01, 1.0, 40)
    tps = S * np.linspace(1.0, 2.0, 21)
    Ht, Htp = H(taus, a)[:, far], H(tps, a)[:, far]
    ok = bool(np.all(Ht[:, None, :] <= (taus[:, None, None] / S) ** 2
                     * math.exp(8.0 * S) * Htp[None, :, :] * (1 + 1e-12)))
    check(ok, f"level {level}: H_tau(a,j) <= (tau/s)^2 e^(8s) H_tau'(a,j) on a grid")
    # H_s(p,a) ≥ s e^{-4s}
    Hs_pa = float(H(np.array([S]), p)[0, a])
    check(Hs_pa >= S * math.exp(-4.0 * S), f"level {level}: H_s(p,a)={Hs_pa:.6f} >= s e^(-4s)={S*math.exp(-4*S):.6f}")
    # chain: w_pj ≥ c H_s(p,a) ∫_s^∞ H(a,j) (τ+s)^{-1-α} dτ ≥ κ w_aj
    shifted = np.trapezoid(Ha_hi * (t_hi / (t_hi + S) ** (1.0 + alpha))[:, None], y, axis=0)
    shifted = shifted + tail * (T_end / (T_end + S)) ** (1.0 + alpha)
    chain = c_alpha * Hs_pa * shifted
    ok = all(-A[p, j] >= chain[j] - 1e-12 and chain[j] >= kappa(alpha) * (-A[a, j]) - 1e-12 for j in far)
    check(ok, f"level {level}: w_pj >= c H_s(p,a) ∫_s^∞ H_tau(a,j)(tau+s)^(-1-a) >= kappa w_aj for every far j")
    # tail split ∫_0^∞ ≤ (1+K) ∫_s^∞
    ok = all(I_lo[j] + I_hi[j] <= (1.0 + K) * I_hi[j] for j in far)
    check(ok, f"level {level}: ∫_0^∞ <= (1+K) ∫_s^∞ for every far j")


def main(through: int = 6) -> int:
    assert through >= 4
    print("Theorem AH: uniform ceiling on the gap ratio gamma_n; positive audit floor")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0).")
    print("Not thrust.\n")

    k = kappa(ALPHA)
    Kc = K_small(ALPHA)
    hc = harnack_const(ALPHA)
    Gc = gamma_ceil(ALPHA)
    sm = sigma_min(ALPHA)
    rm = rho_min(ALPHA)
    fl = sm * AA_FRAC_FLOOR_045
    flAC = nF_floor_AC(ALPHA)
    print(f"s = {S};  K(0.45) = {Kc:.9f};  kappa(0.45) = {k:.9e}")
    print(f"4^a/(2 kappa) = {hc:.6f};  Gamma(0.45) = {Gc:.4f}")
    print(f"sigma_min = {sm:.9e};  rho_min = {rm:.9e}")
    print(f"||F|| floor (AA) = sigma_min * 0.084015 = {fl:.9e}")
    print(f"||F|| floor (AC) = {flAC:.9e}\n")

    check(abs(Kc - math.exp(0.8) * 0.45 / (1.55 * (1.0 - 2.0 ** -0.45))) < 1e-14, "K closed form")
    check(abs(k - 0.1 * math.exp(-0.4) * 2.0 ** -1.45 / (1.0 + Kc)) < 1e-16, "kappa closed form")
    check(abs(k - 7.1923e-3) < 1e-7, "kappa(0.45) ~ 7.1923e-3")
    check(abs(Gc - 4.0 ** ALPHA / (k * w_far_min(ALPHA))) < 1e-9 and Gc < 4500.0, "Gamma(0.45) < 4500")
    check(abs(fl - 6.442e-4) < 1e-6 and fl > 0.0, "||F(n,0.45)|| floor ~ 6.442e-4 > 0")
    check(0.0 < flAC < fl, "AC route floor positive and weaker than AA route at 0.45")
    grid = np.linspace(0.01, 0.99, 99)
    check(all(kappa(a) > 0.0 and gamma_ceil(a) < math.inf and sigma_min(a) > 0.0 for a in grid),
          "kappa, Gamma, sigma_min finite positive on a grid of a in (0,1)")
    check(all(nF_floor_AC(a) > 0.0 for a in grid if a < GAMMA_C),
          "||F|| floor via AC positive on the grid below log3/log5")

    # Heat-kernel steps, by quadrature.
    for n in (1, 2, 3):
        heat_kernel_checks(n, ALPHA)
    heat_kernel_checks(2, 0.10)
    heat_kernel_checks(2, 0.90)

    # Live table at alpha = 0.45.
    rows = []
    for n in range(1, through + 1):
        r = level_data(n, ALPHA)
        rows.append(r)
        print(
            f"level {n} N={r['N']} gamma={r['gamma']:.8f} sigma={r['sigma']:.8f} "
            f"Ifar/Ia={r['Ifar']/r['Ia']:.8f} min w_pj/w_aj={r['ratio_p']:.8f} "
            f"||F||={r['nF']:.8f}",
            flush=True,
        )

    for n, r in enumerate(rows, start=1):
        check(r["far_dist2"], f"level {n}: Far vertices are non-neighbors of a")
        check(r["deg_p"] == 4.0 and r["deg_q"] == 4.0, f"level {n}: deg p = deg q = 4")
        check(r["diag_p"] <= 4.0 ** ALPHA + 1e-12 and r["diag_q"] <= 4.0 ** ALPHA + 1e-12,
              f"level {n}: (L^a)_pp={r['diag_p']:.8f} <= 4^a (Jensen)")
        check(r["vmin"] >= -1e-12, f"level {n}: v = u_a - u >= 0 (max principle)")
        check(r["res_p"] < 1e-10 and r["res_q"] < 1e-10, f"level {n}: v L^a-harmonic at p, q")
        check(min(r["ratio_p"], r["ratio_q"]) >= k,
              f"level {n}: min_j w_pj/w_aj = {min(r['ratio_p'], r['ratio_q']):.8f} >= kappa")
        check(r["Ifar"] <= hc * r["Ia"] + 1e-12,
              f"level {n}: I_far={r['Ifar']:.8f} <= (4^a/(2kappa)) I_a={hc*r['Ia']:.4f}")
        check(abs(r["Iloc_sum"] - r["Iaf"]) < 1e-10, f"level {n}: I_a^a = w_n I_a + I_far")
        check(r["gamma"] <= 4.0 ** ALPHA / (k * r["Wfar"]) + 1e-9 and r["gamma"] <= Gc,
              f"level {n}: gamma={r['gamma']:.8f} <= 4^a/(kappa W_far)={4**ALPHA/(k*r['Wfar']):.2f} <= Gamma")
        check(r["sigma"] >= sm and r["rho"] >= rm,
              f"level {n}: sigma={r['sigma']:.6f} >= sigma_min; rho={r['rho']:.6f} >= rho_min")
        check(abs(r["nF"] - SQRT21_5 * r["Ia"]) < 1e-10 and abs(r["nF_frac"] - SQRT21_5 * r["Iaf"]) < 1e-10,
              f"level {n}: ||F|| = (sqrt21/5) I_a, ||F_frac|| = (sqrt21/5) I_a^a")
        check(r["nF_frac"] >= AA_FRAC_FLOOR_045 and r["nF"] >= sm * r["nF_frac"] >= fl,
              f"level {n}: ||F||={r['nF']:.8f} >= sigma_min ||F_frac|| >= {fl:.3e}")

    # Recorded AE table (unchanged).
    gam_rec = (1.443835, 1.390438, 1.378096, 1.375485, 1.374963, 1.374864)
    for i, r in enumerate(rows[:6]):
        check(abs(r["gamma"] - gam_rec[i]) < 5e-7, f"level {i+1}: gamma matches recorded AE table")
    check(abs(rows[1]["nF"] - 1.165909) < 5e-6, "level 2 ||F|| lock 1.165909")

    # Other alpha in (0,1), through level 5.
    for a_ in (0.10, 0.30, 0.50, 0.60, 0.68, 0.80, 0.95):
        for n in range(1, 6):
            r = level_data(n, a_)
            ok = (
                min(r["ratio_p"], r["ratio_q"]) >= kappa(a_)
                and r["Ifar"] <= harnack_const(a_) * r["Ia"] + 1e-12
                and r["gamma"] <= gamma_ceil(a_)
                and r["sigma"] >= sigma_min(a_)
                and r["vmin"] >= -1e-12
            )
            if a_ < GAMMA_C:
                ok = ok and r["nF"] >= nF_floor_AC(a_)
            check(ok, f"a={a_:.2f} level {n}: kernel ratio {min(r['ratio_p'], r['ratio_q']):.6f} >= "
                      f"kappa={kappa(a_):.3e}; gamma={r['gamma']:.6f} <= Gamma={gamma_ceil(a_):.1f}; "
                      f"sigma={r['sigma']:.6f} >= {sigma_min(a_):.3e}"
                      + (f"; ||F||={r['nF']:.6f} >= {nF_floor_AC(a_):.3e}" if a_ < GAMMA_C else ""))

    # Kept Failure AH.1: Theorem J's gap cap alone does not floor I_a.
    shortcut = [(r["Iaf"] - 1.5 * r["Wfar"]) / r["wn"] for r in rows]
    print("shortcut lower bounds (I_a^a - 1.5 W_far)/w_n:", [f"{x:.6f}" for x in shortcut])
    check(shortcut[0] > 0.0 and shortcut[1] > 0.0, "shortcut bound positive at levels 1, 2")
    check(all(rows[i]["Iaf"] < 1.5 * rows[i]["Wfar"] for i in range(2, len(rows))),
          "Kept Failure AH.1: I_a^a < (3/2) W_far at levels 3..6 (shortcut bound negative)")
    check(all(x < 0.0 for x in shortcut[2:]), "Kept Failure AH.1: shortcut gives no floor from level 3")

    # Observations (not theorems).
    rp = [min(r["ratio_p"], r["ratio_q"]) for r in rows]
    check(all(x >= 1.0 for x in rp),
          "observation: live min_j w_pj/w_aj >= 1 through the range (proved kappa ~ 7.2e-3 is crude)")
    check(all(r["nF"] > 0.70 for r in rows), "observation: ||F|| > 0.70 through the range (J.1 not proved)")

    print()
    print("Theorem AH: w_pj >= kappa w_aj on Far; I_far <= (4^a/(2kappa)) I_a; gamma_n <= Gamma(a);")
    print("sigma_n >= sigma_min; for a < log3/log5, ||F(n,a)|| >= sigma_min ||F_frac|| > 0 uniformly.")
    print("At a=0.45: ||F(n,0.45)|| >= 6.44e-4 at every level, so ||F|| does not tend to 0.")
    print("Kept Failure AH.1: Theorem J's 3/2 gap cap alone gives no floor from level 3.")
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
