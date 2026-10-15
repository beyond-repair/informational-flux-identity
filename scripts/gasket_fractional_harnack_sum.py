#!/usr/bin/env python3
"""Two-neighbor heat-kernel sum; sharpened neighbor Harnack and audit floor.

Witness for Theorem AI and Kept Failure AI.1 in GASKET.md.

Notation of Theorems AE–AH: a is the value-1 corner (deg a = 2), p, q its
two neighbors (deg 4, p ~ q), Far = V \\ {a, p, q}, w_xy = -(L^α)_{xy} > 0
for x != y, H_t = e^{-tL}, I_a = (L u)_a, I_a^α = (L^α u)_a,
I_far = Σ_{j∈Far} w_aj (u_a - u_j), W_far = Σ_{j∈Far} w_aj,
γ_n = 2 I_far / (W_far I_a), σ_n = I_a / I_a^α = ||F|| / ||F_frac||.

Key identity (backward equation at the degree-2 corner): for j != a,
    ∂_t H_t(a,j) + 2 H_t(a,j) = H_t(p,j) + H_t(q,j),
so H_t(a,j) = ∫_0^t e^{-2(t-s)} (H_s(p,j) + H_s(q,j)) ds.

Theorem AI (every α ∈ (0,1), every level n ≥ 1):
  (1) two-neighbor kernel sum: for every j ∈ Far,
        w_pj + w_qj ≥ 2 w_aj
      (subordination, Fubini, K(s) = ∫_s^∞ e^{-2(t-s)} t^{-1-α} dt ≤ s^{-1-α}/2);
  (2) sharpened one-step Harnack:  I_far ≤ (4^α / 2) I_a
      (AH (2) with the p and q identities added before comparing kernels);
  (3) γ_n ≤ 4^α / W_far ≤ Γ_AI(α) := 4^α / W_far_min;
  (4) σ_n ≥ σ_AI(α) := 1/(2^{α-1} + 4^α/2),
      ρ_n ≥ ρ_AI(α) := w_min/(w_min + 4^α/2);
  (5) for α < log3/log5: ||F(n,α)|| ≥ σ_AI · 9√21/(28 C_AC(α)) > 0 (AC);
      at α = 0.45 with AA's 0.084015: ||F(n,0.45)|| ≥ σ_AI(0.45)·0.084015 ≈ 0.0520.
  Existence of the limit and Conjecture J.1's value L ≥ 0.70 are NOT proved.

Kept Failure AI.1: "replacing κ by 1 in AH and combining with AA's floor
yields ||F(n,0.45)|| ≥ 0.70 (Conjecture J.1's value)" is false: the chain
gives σ_AI · 0.084015 ≈ 0.0520, and even with the live ||F_frac(n)|| in
place of AA's floor, σ_AI ||F_frac(n)|| < 0.70 at every computed level, and
no Harnack constant alone rescues the AA route (constant 0 gives
2^{1-α} · 0.084015 ≈ 0.1230 < 0.70).

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
from gasket_fractional_gamma_ceiling import (  # noqa: E402
    gamma_ceil as gamma_ceil_AH,
    harnack_const as harnack_const_AH,
    kappa,
    sigma_min as sigma_min_AH,
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
    """Theorem AI (2): I_far ≤ (4^α/2) I_a."""
    return 4.0 ** alpha / 2.0


def gamma_ceil(alpha: float) -> float:
    """Theorem AI (3): Γ_AI(α) = 4^α / W_far_min."""
    return 4.0 ** alpha / w_far_min(alpha)


def sigma_AI(alpha: float) -> float:
    """Theorem AI (4): σ_n ≥ 1/(2^{α-1} + 4^α/2)."""
    return 1.0 / (2.0 ** (alpha - 1.0) + harnack_const(alpha))


def rho_AI(alpha: float) -> float:
    wm = w_min(alpha)
    return wm / (wm + harnack_const(alpha))


def nF_floor_AC(alpha: float) -> float:
    """Theorem AI (5) via AC: σ_AI · 9√21/(28 C_AC)."""
    return sigma_AI(alpha) * 9.0 * math.sqrt(21.0) / (28.0 * C_AC(alpha))


_K_CACHE: dict[float, tuple[np.ndarray, np.ndarray]] = {}
Z_GRID = np.linspace(math.log(1e-30), math.log(1e12), 14001)  # s = e^z


def K_weight(alpha: float) -> tuple[np.ndarray, np.ndarray]:
    """K(s) = ∫_s^∞ e^{-2(t-s)} t^{-1-α} dt = ∫_0^∞ e^{-2u} (s+u)^{-1-α} du on s = e^Z_GRID.

    Quadrature in u = e^x, x ∈ [-80, 5] (smooth integrand; the piece u < e^{-80}
    is below e^{-80} s^{-1-α}, negligible against K(s) on the grid).
    """
    if alpha not in _K_CACHE:
        s = np.exp(Z_GRID)
        x = np.linspace(-80.0, 5.0, 8501)
        uu = np.exp(x)
        out = np.empty(s.size)
        for i0 in range(0, s.size, 512):
            ss = s[i0:i0 + 512, None]
            vals = np.exp(-2.0 * uu)[None, :] * (ss + uu[None, :]) ** (-1.0 - alpha) * uu[None, :]
            out[i0:i0 + 512] = np.trapezoid(vals, x, axis=1)
        _K_CACHE[alpha] = (s, out)
    return _K_CACHE[alpha]


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
        "Ia": Ia,
        "Iaf": Iaf,
        "Ifar": Ifar,
        "Wfar": Wfar,
        "gamma": 2.0 * Ifar / (Wfar * Ia),
        "sigma": Ia / Iaf,
        "rho": wn * Ia / Iaf,
        "sum_ratio_min": float(sum_ratio.min()),
        "ratio_p": float(np.min(W[p, far] / W[a, far])),
        "ratio_q": float(np.min(W[q, far] / W[a, far])),
        "sum_slack_min": float(np.min(W[p, far] + W[q, far] - 2.0 * W[a, far])),
        "harn_p": float(A[p, p] * v[p] - np.dot(W[p, far], v[far])),
        "harn_q": float(A[q, q] * v[q] - np.dot(W[q, far], v[far])),
        "nF": corner_flux_norm(P, corners, I),
        "nF_frac": corner_flux_norm(P, corners, If),
        "Iloc_sum": float(wn * Ia + Ifar),
    }


def heat_kernel_checks(level: int, alpha: float) -> None:
    """Quadrature checks of the first-step identity behind Theorem AI (1)."""
    P, edges, corners = build_gasket(level)
    N = len(P)
    L = combinatorial_laplacian(N, edges)
    lam, V = np.linalg.eigh(L)
    lam = np.clip(lam, 0.0, None)
    A = spectral_power(L, alpha)
    a = int(corners[0])
    nbrs = [j for j in range(N) if j != a and L[a, j] != 0.0]
    p, q = nbrs
    far = [j for j in range(N) if j not in (a, p, q)]
    c_alpha = alpha / math.gamma(1.0 - alpha)

    # Row structure of L at the degree-2 corner: L_a = 2 e_a - e_p - e_q.
    row = np.zeros(N)
    row[a], row[p], row[q] = 2.0, -1.0, -1.0
    check(np.array_equal(L[a], row), f"level {level}: L row at a is 2e_a - e_p - e_q (deg a = 2)")

    # H_t = e^{-4t} e^{tM}, M = 4I - L >= 0: nonnegative series for t <= 4, eigh beyond.
    M = 4.0 * np.eye(N) - L
    KMAX, T_SER = 160, 4.0
    pows = {}
    for x0 in (a, p, q):
        rk = np.zeros((KMAX + 1, N))
        rk[0, x0] = 1.0
        for kk in range(1, KMAX + 1):
            rk[kk] = M @ rk[kk - 1]
        pows[x0] = rk
    ks = np.arange(KMAX + 1)
    logfact = np.array([math.lgamma(kk + 1.0) for kk in ks])

    def H(t: np.ndarray, x: int) -> np.ndarray:
        t = np.atleast_1d(np.asarray(t, dtype=np.float64))
        out = np.empty((t.size, N))
        zero = t == 0.0
        small = (t <= T_SER) & ~zero
        if np.any(zero):
            out[zero] = pows[x][0]
        if np.any(small):
            ts = t[small]
            coef = np.exp(ks[None, :] * np.log(ts)[:, None] - logfact[None, :] - 4.0 * ts[:, None])
            out[small] = coef @ pows[x]
        big = t > T_SER
        if np.any(big):
            out[big] = (np.exp(-np.outer(t[big], lam)) * V[x]) @ V.T
        return out

    # Backward equation  ∂_t H(a,j) + 2 H(a,j) = H(p,j) + H(q,j)  (j != a), via eigh.
    worst = 0.0
    for t in (0.03, 0.4, 2.5, 30.0):
        Ht = (np.exp(-t * lam) * V) @ V.T
        dHt = -((lam * np.exp(-t * lam)) * V) @ V.T
        resid = dHt[a] + 2.0 * Ht[a] - Ht[p] - Ht[q]
        resid[a] = 0.0
        worst = max(worst, float(np.max(np.abs(resid))))
    check(worst < 1e-12, f"level {level}: backward equation at a holds (max resid {worst:.1e})")

    # Duhamel  H_t(a,j) = ∫_0^t e^{-2(t-s)} (H_s(p,j) + H_s(q,j)) ds  at a few t.
    worst = 0.0
    for t in (0.2, 1.0, 3.0):
        ss = np.linspace(0.0, t, 4001)
        f = H(ss, p) + H(ss, q)
        duh = np.trapezoid(np.exp(-2.0 * (t - ss))[:, None] * f, ss, axis=0)
        Hat = H(np.array([t]), a)[0]
        worst = max(worst, max(abs(duh[j] - Hat[j]) / Hat[j] for j in far))
    check(bool(np.isfinite(worst)) and worst < 1e-6, f"level {level}: Duhamel formula reproduces H_t(a,j) on Far (rel err {worst:.1e})")

    # Weight bound K(s) <= s^{-1-a}/2 on a log grid; and the exact identity
    # w_aj = c_a ∫ (H_s(p,j) + H_s(q,j)) K(s) ds.
    s, Ks = K_weight(alpha)
    z = Z_GRID
    bound = 0.5 * s ** (-1.0 - alpha)
    check(bool(np.all(Ks <= bound)) and bool(np.all(Ks > 0.0)),
          f"level {level}: 0 < K(s) <= s^(-1-a)/2 on s in [1e-30, 1e12] (max K/bound {np.max(Ks/bound):.9f})")
    f = H(s, p) + H(s, q)
    T_end, s0 = float(s[-1]), float(s[0])
    # tails: f -> 2/N and K ~ s^{-1-a}/2 beyond T_end; below s0 the piece is O(s0^{2-a})
    tail = (2.0 / N) * 0.5 * T_end ** (-alpha) / alpha
    rep = c_alpha * (np.trapezoid(f * (Ks * s)[:, None], z, axis=0) + tail)
    err = max(abs(rep[j] + A[a, j]) / abs(A[a, j]) for j in far)
    check(err < 1e-6, f"level {level}: w_aj = c_a ∫ (H_s(p,j)+H_s(q,j)) K(s) ds on Far (rel err {err:.1e})")

    def sub(x0: int) -> np.ndarray:
        # w_x0,j = c_a ∫ H_s(x0,j) s^{-1-a} ds; below s0, H_s(x0,j) ~ s M_x0j
        body = np.trapezoid(H(s, x0) * s[:, None] ** (-alpha), z, axis=0)
        return c_alpha * (body + (1.0 / N) * T_end ** (-alpha) / alpha
                          + M[x0] * s0 ** (1.0 - alpha) / (1.0 - alpha))

    sub_p, sub_q = sub(p), sub(q)
    errp = max(max(abs(sub_p[j] + A[p, j]) / abs(A[p, j]), abs(sub_q[j] + A[q, j]) / abs(A[q, j])) for j in far)
    ok = all(sub_p[j] + sub_q[j] >= 2.0 * rep[j] for j in far)
    check(ok and errp < 1e-6,
          f"level {level}: quadrature w_pj + w_qj >= 2 c_a ∫ f K for every far j (w_pj, w_qj rel err {errp:.1e})")


def main(through: int = 6) -> int:
    assert through >= 4
    print("Theorem AI: two-neighbor heat-kernel sum; sharpened neighbor Harnack and audit floor")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0).")
    print("Not thrust.\n")

    hc = harnack_const(ALPHA)
    Gc = gamma_ceil(ALPHA)
    sm = sigma_AI(ALPHA)
    rm = rho_AI(ALPHA)
    fl = sm * AA_FRAC_FLOOR_045
    flAC = nF_floor_AC(ALPHA)
    fl_AH = sigma_min_AH(ALPHA) * AA_FRAC_FLOOR_045
    print(f"2^(a-1) = {2.0 ** (ALPHA - 1.0):.9f};  4^a/2 = {hc:.9f}")
    print(f"Gamma_AI(0.45) = {Gc:.6f}  (Gamma_AH = {gamma_ceil_AH(ALPHA):.4f}; ratio 1/kappa = {1.0/kappa(ALPHA):.4f})")
    print(f"sigma_AI = {sm:.9f};  rho_AI = {rm:.9f}  (sigma_AH = {sigma_min_AH(ALPHA):.9e})")
    print(f"||F|| floor (AA) = sigma_AI * 0.084015 = {fl:.9f}  (AH: {fl_AH:.6e}; gain x{fl/fl_AH:.2f})")
    print(f"||F|| floor (AC) = {flAC:.9e}\n")

    check(abs(hc - 0.5 * 4.0 ** 0.45) < 1e-15 and abs(hc - harnack_const_AH(ALPHA) * kappa(ALPHA)) < 1e-12,
          "Harnack constant 4^a/2 is AH's 4^a/(2 kappa) with kappa replaced by 1")
    check(abs(sm - 1.0 / (2.0 ** -0.55 + 0.5 * 4.0 ** 0.45)) < 1e-15, "sigma_AI closed form")
    check(abs(sm - 0.6187915) < 1e-7, "sigma_AI(0.45) ~ 0.6187915")
    check(abs(fl - 0.051988) < 1e-6 and fl > fl_AH, "||F(n,0.45)|| floor ~ 0.051988 > AH's 6.442e-4")
    check(abs(Gc - 4.0 ** ALPHA / w_far_min(ALPHA)) < 1e-12 and Gc < gamma_ceil_AH(ALPHA), "Gamma_AI < Gamma_AH")
    check(0.0 < flAC < fl, "AC route floor positive and weaker than AA route at 0.45")
    grid0 = np.linspace(0.01, 0.99, 99)
    check(all(abs(gamma_ceil(a_) - 8.0 / (a_ * (1.0 - a_))) < 1e-9 * gamma_ceil(a_)
              and abs(rho_AI(a_) - a_ / (a_ + 2.0)) < 1e-14
              and abs(sigma_AI(a_) - 2.0 / (2.0 ** a_ * (1.0 + 2.0 ** a_))) < 1e-14
              and abs(gamma_ceil_AH(a_) * kappa(a_) - gamma_ceil(a_)) < 1e-9 * gamma_ceil(a_)
              for a_ in grid0),
          "closed forms: Gamma_AI = 8/(a(1-a)), rho_AI = a/(a+2), sigma_AI = 2/(2^a(1+2^a)), Gamma_AH = Gamma_AI/kappa")
    grid = np.linspace(0.01, 0.99, 99)
    check(all(sigma_AI(a_) > sigma_min_AH(a_) and gamma_ceil(a_) < gamma_ceil_AH(a_) for a_ in grid),
          "sigma_AI > sigma_AH and Gamma_AI < Gamma_AH on a grid of a in (0,1)")
    check(all(nF_floor_AC(a_) > 0.0 for a_ in grid if a_ < GAMMA_C),
          "||F|| floor via AC positive on the grid below log3/log5")

    # First-step identity and weight bound, by quadrature.
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
            f"level {n} N={r['N']} min(w_pj+w_qj)/(2w_aj)={r['sum_ratio_min']:.8f} "
            f"min w_pj/w_aj={min(r['ratio_p'], r['ratio_q']):.8f} Ifar/Ia={r['Ifar']/r['Ia']:.8f} "
            f"gamma={r['gamma']:.8f} sigma={r['sigma']:.8f} ||F||={r['nF']:.8f} ||F_frac||={r['nF_frac']:.8f}",
            flush=True,
        )

    for n, r in enumerate(rows, start=1):
        check(r["deg_a"] == 2.0 and r["deg_p"] == 4.0 and r["deg_q"] == 4.0 and r["p_adj_q"],
              f"level {n}: deg a = 2, deg p = deg q = 4, p ~ q")
        check(r["sum_slack_min"] >= 0.0 and r["sum_ratio_min"] >= 1.0,
              f"level {n}: w_pj + w_qj >= 2 w_aj for every far j (min ratio {r['sum_ratio_min']:.8f})")
        check(r["diag_p"] <= 4.0 ** ALPHA + 1e-12 and r["diag_q"] <= 4.0 ** ALPHA + 1e-12,
              f"level {n}: (L^a)_pp={r['diag_p']:.8f} <= 4^a (Jensen)")
        check(r["vmin"] >= -1e-12, f"level {n}: v = u_a - u >= 0 (max principle)")
        check(r["res_p"] < 1e-10 and r["res_q"] < 1e-10, f"level {n}: v L^a-harmonic at p, q")
        check(r["harn_p"] >= -1e-12 and r["harn_q"] >= -1e-12,
              f"level {n}: (L^a)_pp v_p >= sum_Far w_pj v_j (same at q)")
        check(r["Ifar"] <= hc * r["Ia"] + 1e-12,
              f"level {n}: I_far={r['Ifar']:.8f} <= (4^a/2) I_a={hc*r['Ia']:.8f}")
        check(abs(r["Iloc_sum"] - r["Iaf"]) < 1e-10, f"level {n}: I_a^a = w_n I_a + I_far")
        check(r["gamma"] <= 4.0 ** ALPHA / r["Wfar"] + 1e-9 and r["gamma"] <= Gc,
              f"level {n}: gamma={r['gamma']:.8f} <= 4^a/W_far={4**ALPHA/r['Wfar']:.4f} <= Gamma_AI={Gc:.4f}")
        check(r["sigma"] >= sm and r["rho"] >= rm,
              f"level {n}: sigma={r['sigma']:.6f} >= sigma_AI; rho={r['rho']:.6f} >= rho_AI")
        check(abs(r["nF"] - SQRT21_5 * r["Ia"]) < 1e-10 and abs(r["nF_frac"] - SQRT21_5 * r["Iaf"]) < 1e-10,
              f"level {n}: ||F|| = (sqrt21/5) I_a, ||F_frac|| = (sqrt21/5) I_a^a")
        check(r["nF_frac"] >= AA_FRAC_FLOOR_045 and r["nF"] >= sm * r["nF_frac"] >= fl,
              f"level {n}: ||F||={r['nF']:.8f} >= sigma_AI ||F_frac|| >= {fl:.6f}")

    gam_rec = (1.443835, 1.390438, 1.378096, 1.375485, 1.374963, 1.374864)
    for i, r in enumerate(rows[:6]):
        check(abs(r["gamma"] - gam_rec[i]) < 5e-7, f"level {i+1}: gamma matches recorded AE table")
    check(abs(rows[1]["nF"] - 1.165909) < 5e-6, "level 2 ||F|| lock 1.165909")

    # Other alpha in (0,1), through level 5.
    for a_ in (0.10, 0.30, 0.50, 0.60, 0.68, 0.80, 0.95):
        for n in range(1, 6):
            r = level_data(n, a_)
            ok = (
                r["sum_slack_min"] >= 0.0
                and r["harn_p"] >= -1e-12 and r["harn_q"] >= -1e-12
                and r["Ifar"] <= harnack_const(a_) * r["Ia"] + 1e-12
                and r["gamma"] <= gamma_ceil(a_)
                and r["sigma"] >= sigma_AI(a_)
                and r["rho"] >= rho_AI(a_)
                and r["vmin"] >= -1e-12
            )
            if a_ < GAMMA_C:
                ok = ok and r["nF"] >= nF_floor_AC(a_)
            check(ok, f"a={a_:.2f} level {n}: min(w_pj+w_qj)/(2w_aj)={r['sum_ratio_min']:.6f} >= 1; "
                      f"gamma={r['gamma']:.6f} <= Gamma_AI={gamma_ceil(a_):.3f}; "
                      f"sigma={r['sigma']:.6f} >= {sigma_AI(a_):.6f}"
                      + (f"; ||F||={r['nF']:.6f} >= {nF_floor_AC(a_):.3e}" if a_ < GAMMA_C else ""))
    # Kept Failure AI.1: kappa -> 1 plus AA does not reach Conjecture J.1's 0.70.
    print(f"\nKept Failure AI.1: sigma_AI * 0.084015 = {fl:.6f} vs 0.70")
    chain_live = [sm * r["nF_frac"] for r in rows]
    print("sigma_AI * live ||F_frac(n)||:", [f"{x:.6f}" for x in chain_live])
    check(fl < J1_VALUE, "Kept Failure AI.1: sigma_AI * 0.084015 < 0.70 (chain does not give J.1)")
    check(all(x < J1_VALUE for x in chain_live),
          "Kept Failure AI.1: even sigma_AI * live ||F_frac(n)|| < 0.70 at every computed level")
    check(sm < 1.0 and AA_FRAC_FLOOR_045 < J1_VALUE, "Kept Failure AI.1: sigma_AI < 1 and AA floor < 0.70")
    perfect = AA_FRAC_FLOOR_045 / 2.0 ** (ALPHA - 1.0)
    check(perfect < J1_VALUE,
          f"Kept Failure AI.1: even Harnack constant 0 gives only 2^(1-a) * 0.084015 = {perfect:.6f} < 0.70")

    # Observations (not theorems).
    rp = [min(r["ratio_p"], r["ratio_q"]) for r in rows]
    rs = [r["sum_ratio_min"] for r in rows]
    print("live min_j w_pj/w_aj:", [f"{x:.6f}" for x in rp])
    print("live min_j (w_pj+w_qj)/(2 w_aj):", [f"{x:.6f}" for x in rs])
    check(all(x >= 1.0 for x in rp),
          "observation: pointwise w_pj >= w_aj (kappa = 1) holds live through the range; not proved")
    check(all(rs[i + 1] < rs[i] for i in range(len(rs) - 1)) and rs[-1] < 1.001,
          "observation: min sum ratio decreases toward 1 (constant 2 nearly sharp on the range)")
    check(all(r["Ifar"] / r["Ia"] < 0.45 * hc for r in rows),
          "observation: live I_far/I_a < 0.45 (4^a/2) (Harnack constant still not sharp)")
    check(all(r["nF"] > 0.70 for r in rows), "observation: ||F|| > 0.70 through the range (J.1 not proved)")

    print()
    print("Theorem AI: w_pj + w_qj >= 2 w_aj on Far; I_far <= (4^a/2) I_a; gamma_n <= 4^a/W_far_min;")
    print("sigma_n >= 1/(2^(a-1) + 4^a/2); for a < log3/log5, ||F(n,a)|| >= sigma_AI ||F_frac|| > 0.")
    print(f"At a=0.45: ||F(n,0.45)|| >= {fl:.6f} at every level (AH: 6.44e-4).")
    print("Kept Failure AI.1: kappa -> 1 with AA's floor gives ~0.052, not Conjecture J.1's 0.70.")
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
