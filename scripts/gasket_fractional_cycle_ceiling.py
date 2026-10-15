#!/usr/bin/env python3
"""Boundary-cycle resistance ceiling via Fourier and the heat kernel (Theorem AM).

Witness for Theorem AM and Kept Failure AM.1 in GASKET.md.

Theorems AA/AL compared the gasket with ONE side a–b (a path P_m, m = 2^n),
whose endpoint resistance R_α(m) is capped by R_∞(α) = 4^{-α}(4/π)B(1/2−α,3/2).
Theorem AM keeps all three outer sides: the boundary cycle C_N, N = 3m, with
the corners a, b, c at mutual cycle distance m. Its Laplacian M (boundary edges
only, interior vertices isolated) satisfies 0 ≤ M ≤ L, so M^α ≤ L^α
(Löwner–Heinz, α ∈ (0,1)). On the cycle,

  R^cyc_α(m) := (e_a−e_b)^T (M^α)^+ (e_a−e_b)
             = (2/N) Σ_{k=1}^{N−1} λ_k^{-α} (1 − cos(2πk/3)),  λ_k = 4 sin²(πk/N),

and the heat-kernel representation λ^{-α} = Γ(α)^{-1} ∫ t^{α−1} e^{−tλ} dt with
the periodized walk p^C_t(x) = Σ_j p^Z_t(x + jN) gives, because p^Z_t(x) is
strictly decreasing in |x|,

  p^C_t(0) − p^C_t(m) = p^Z_t(0) − Σ_{r≥0}[p((3r+1)m) + p((3r+2)m) − 2p((3r+3)m)]
                     < p^Z_t(0),

hence for every α ∈ (0, 1/2) and every m ≥ 1

  R^cyc_α(m) < R^cyc_∞(α) := 2 ∫_0^1 (2 sin πx)^{-2α} dx
                          = (2/π) 4^{-α} B(1/2−α, 1/2) = (1−α) R_∞(α).

The cycle also carries all three corners; rotation by m makes the corner block
of (M^α)^+ circulant, so the best corner-supported test vector gives
  Q_α(u_n) ≥ u^T M^α u ≥ 2|u_C − ū_C|² / R^cyc = (7/3)/R^cyc   (data (1,−1/2,0)).
At α = 0.45: R^cyc_∞ ≈ 7.284859, Q ≥ 0.320299, ||F_frac|| ≥ 0.209685
(≈1.886× AL's 0.111208), and with σ_AK, ||F(n,0.45)|| ≥ 0.147106.

Kept Failure AM.1: still below Conjecture J.1's 0.70; and the σ_AK × ||F_frac||
chain cannot reach 0.70 with ANY resistance ceiling, because even the live
energy gives σ_AK·||F_frac||_live < 0.70 from level 2 on (≈0.457 at level 6).

Not thrust.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_energy_floor import (  # noqa: E402
    R_limit,
    beta,
    path_R_closed,
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
AL_FRAC_FLOOR_045 = 0.111208
AL_F_FLOOR_045 = 0.078019
Q_RECORDED = {9: 0.948091, 10: 0.941537}  # GASKET.md, Kept Failures J.3/J.4
SQ3 = math.sqrt(3.0)

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


# ---------------------------------------------------------------- cycle side
def cycle_R(m: int, alpha: float) -> float:
    """R^cyc_α(m) for corners at distance m on C_{3m} (Fourier; 3∤k terms only)."""
    N = 3 * m
    k = np.arange(1, N, dtype=np.float64)
    k = k[np.arange(1, N) % 3 != 0]
    lam = 4.0 * np.sin(np.pi * k / N) ** 2
    return float((3.0 / N) * np.sum(lam ** (-alpha)))


def cycle_R_fourier_raw(m: int, alpha: float) -> float:
    """Literal (2/N) Σ λ_k^{-α}(1 − cos(2πk m/N)), no simplification."""
    N = 3 * m
    k = np.arange(1, N, dtype=np.float64)
    lam = 2.0 - 2.0 * np.cos(2.0 * np.pi * k / N)
    return float((2.0 / N) * np.sum(lam ** (-alpha) * (1.0 - np.cos(2.0 * np.pi * k * m / N))))


def cycle_laplacian(N: int) -> np.ndarray:
    C = 2.0 * np.eye(N)
    for i in range(N):
        C[i, (i + 1) % N] -= 1.0
        C[(i + 1) % N, i] -= 1.0
    return C


def pinv_power(Mat: np.ndarray, alpha: float):
    w, V = np.linalg.eigh(Mat)
    w = np.clip(w, 0.0, None)
    keep = w > 1e-10
    Ap = (V[:, keep] * (w[keep] ** alpha)) @ V[:, keep].T
    Gp = (V[:, keep] * (w[keep] ** (-alpha))) @ V[:, keep].T
    return Ap, Gp


def R_cyc_inf(alpha: float) -> float:
    return (2.0 / math.pi) * 4.0 ** (-alpha) * beta(0.5 - alpha, 0.5)


def R_cyc_inf_quadrature(alpha: float, nodes: int = 400) -> float:
    """2∫_0^1 (2 sin πx)^{-2α} dx = (4/π)∫_0^{π/2} (2 sin θ)^{-2α} dθ, singularity removed."""
    b = 2.0 * alpha
    x, wts = np.polynomial.legendre.leggauss(nodes)
    u = 0.5 * (x + 1.0)
    wu = 0.5 * wts
    p = 1.0 / (1.0 - b)
    th = (math.pi / 2.0) * u ** p
    jac = (math.pi / 2.0) * p * u ** (p - 1.0)
    return float((4.0 / math.pi) * np.sum(wu * (2.0 * np.sin(th)) ** (-b) * jac))


def pZ(t: float, xs: np.ndarray, K: int = 8192) -> np.ndarray:
    """Continuous-time SRW heat kernel on Z: ∫_0^1 e^{-2t(1-cos 2πθ)} cos(2πxθ) dθ."""
    th = (np.arange(K) + 0.5) / K
    g = np.exp(-2.0 * t * (1.0 - np.cos(2.0 * np.pi * th)))
    return np.array([float(np.mean(g * np.cos(2.0 * np.pi * x * th))) for x in xs])


def log_pZ(t: float, x: int, terms: int = 3000) -> float:
    """log p^Z_t(x) = log(e^{-2t} I_|x|(2t)) by the Bessel series in log space."""
    x = abs(int(x))
    k = np.arange(terms, dtype=np.float64)
    lg = np.array([math.lgamma(kk + 1.0) + math.lgamma(kk + x + 1.0) for kk in k])
    terms_log = (2.0 * k + x) * math.log(t) - lg  # (z/2)^{2k+x} with z = 2t
    mx = float(terms_log.max())
    return -2.0 * t + mx + math.log(float(np.sum(np.exp(terms_log - mx))))


def pC_diff(t: float, m: int) -> float:
    """p^C_t(0) − p^C_t(m) on C_{3m}, directly from the cycle spectrum."""
    N = 3 * m
    k = np.arange(N, dtype=np.float64)
    lam = 4.0 * np.sin(np.pi * k / N) ** 2
    return float(np.sum(np.exp(-t * lam) * (1.0 - np.cos(2.0 * np.pi * k / 3.0))) / N)


def Q_floor_cyc(alpha: float) -> float:
    return (7.0 / 3.0) / R_cyc_inf(alpha)


def Ffrac_floor_cyc(alpha: float) -> float:
    return Q_floor_cyc(alpha) * SQRT21_7


# --------------------------------------------------------------- gasket side
def boundary_cycle(P: np.ndarray, edges, corners):
    """Vertices/edges of the three outer sides; edge kept iff both ends lie on one side."""
    x, y = P[:, 0], P[:, 1]
    on = [np.abs(y) < 1e-9, np.abs(y - SQ3 * x) < 1e-9, np.abs(y - SQ3 * (1.0 - x)) < 1e-9]
    bset = set(np.flatnonzero(on[0] | on[1] | on[2]).tolist())
    bedges = [(i, j) for (i, j) in edges if any(s[i] and s[j] for s in on)]
    return sorted(bset), bedges


def cycle_order(bverts, bedges, start: int):
    adj: dict[int, list[int]] = {v: [] for v in bverts}
    for i, j in bedges:
        adj[i].append(j)
        adj[j].append(i)
    if any(len(adj[v]) != 2 for v in bverts):
        return None
    order = [start]
    prev, cur = None, start
    while True:
        nxt = [w for w in adj[cur] if w != prev]
        nxt = nxt[0] if prev is not None else adj[cur][0]
        if nxt == start:
            break
        order.append(nxt)
        prev, cur = cur, nxt
        if len(order) > len(bverts):
            return None
    return order


def level_cycle_record(level: int, alpha: float) -> dict:
    P, edges, corners = build_gasket(level)
    N = len(P)
    m = 2 ** level
    L = combinatorial_laplacian(N, edges)
    bverts, bedges = boundary_cycle(P, edges, corners)
    order = cycle_order(bverts, bedges, int(corners[0]))
    pos = {v: i for i, v in enumerate(order)} if order else {}
    M = combinatorial_laplacian(N, bedges)
    A, Lp = pinv_power(L, alpha)
    Ma, Mp = pinv_power(M, alpha)
    a, b, c = (int(z) for z in corners)
    v = np.zeros(N)
    v[a], v[b] = 1.0, -1.0
    interior = np.setdiff1d(np.arange(N), corners)
    u = np.zeros(N)
    u[corners] = BC
    u[interior] = np.linalg.solve(A[np.ix_(interior, interior)], -A[np.ix_(interior, corners)] @ BC)
    G = Mp[np.ix_(corners, corners)]
    # Schur complement of M^α on the corners = min boundary-extension energy
    rest = np.setdiff1d(np.arange(N), corners)
    bmask = np.zeros(N, bool)
    bmask[bverts] = True
    rest_b = np.array([i for i in rest if bmask[i]], dtype=int)
    S = Ma[np.ix_(corners, corners)] - Ma[np.ix_(corners, rest_b)] @ np.linalg.solve(
        Ma[np.ix_(rest_b, rest_b)], Ma[np.ix_(rest_b, corners)]
    )
    If = (A @ u)[corners]
    I = (L @ u)[corners]
    return {
        "N": N,
        "m": m,
        "nb": len(bverts),
        "ne": len(bedges),
        "cycle_ok": order is not None and len(order) == 3 * m,
        "corner_pos": sorted((pos.get(b, -1), pos.get(c, -1))) if pos else None,
        "LminusM": float(np.linalg.eigvalsh(L - M).min()),
        "LH": float(np.linalg.eigvalsh(A - Ma).min()),
        "Rcyc_live": float(v @ Mp @ v),
        "Rgasket": float(v @ Lp @ v),
        "G": G,
        "Q": float(u @ A @ u),
        "Q_cyc": float(u @ Ma @ u),
        "schur_energy": float(BC @ S @ BC),
        "nF_frac": corner_flux_norm(P, corners, If),
        "nF": corner_flux_norm(P, corners, I),
    }


def main() -> int:
    print("Theorem AM: boundary-cycle resistance ceiling via Fourier / heat kernel")
    print("Assumption: spectral calculus (GASKET item 4); Löwner–Heinz; corner data")
    print("(1,-1/2,0); Theorems AA, AL, AK. Not thrust.\n")

    # 1. Fourier formula vs dense cycle pseudo-inverse.
    worst = 0.0
    for alpha in (0.05, 0.25, 0.45, 0.49, 0.9):
        for m in (1, 2, 3, 4, 5, 8, 16, 32):
            Nn = 3 * m
            _, Gp = pinv_power(cycle_laplacian(Nn), alpha)
            dense = Gp[0, 0] + Gp[m, m] - 2.0 * Gp[0, m]
            worst = max(worst, abs(dense - cycle_R(m, alpha)), abs(dense - cycle_R_fourier_raw(m, alpha)))
    check(worst < 1e-9, f"Fourier R^cyc(m) (raw and 3∤k forms) match dense (C_3m^a)^+ (max err {worst:.1e})")

    # 2. Closed form of the limit: (2/π)4^{-a}B(1/2-a,1/2) = (1-a) R_inf(a) = 2∫(2 sin πx)^{-2a}.
    for alpha in (0.02, 0.05, 0.10, 0.25, 0.45, 0.49):
        Rc, Rp = R_cyc_inf(alpha), R_limit(alpha)
        quad = R_cyc_inf_quadrature(alpha)
        check(abs(Rc - (1.0 - alpha) * Rp) < 1e-10 * Rp and abs(Rc - quad) < 1e-8 * Rc,
              f"alpha={alpha}: R^cyc_inf={Rc:.10f} = (1-a)R_inf={(1-alpha)*Rp:.10f} = quadrature {quad:.10f}")

    # 3. Heat-kernel step: p^Z_t decreasing in |x| (log-space Bessel series, no
    #    underflow); periodization identity; bracket terms > 0; p^C(0)-p^C(m) <= p^Z(0).
    for t in (0.05, 0.3, 1.0, 3.0, 10.0, 40.0, 150.0):
        lp = [log_pZ(t, x) for x in range(0, 61)]
        check(all(lp[i + 1] < lp[i] for i in range(60)),
              f"t={t}: log p^Z_t(x) strictly decreasing in x=0..60 (I_x(2t) decreasing in x)")
        qz = pZ(t, np.array([0.0]))[0]
        check(abs(qz - math.exp(lp[0])) < 1e-12, f"t={t}: p^Z_t(0) quadrature {qz:.12f} = series")
        for m in (1, 2, 3, 5, 8, 16):
            Nn = 3 * m
            js = np.arange(-40, 41)
            per = float(np.sum(pZ(t, np.abs(js * Nn))) - np.sum(pZ(t, np.abs(js * Nn + m))))
            direct = pC_diff(t, m)
            br_ok = all(
                log_pZ(t, (3 * r + 1) * m) > log_pZ(t, (3 * r + 3) * m)
                and log_pZ(t, (3 * r + 2) * m) > log_pZ(t, (3 * r + 3) * m)
                for r in range(0, 4)
            )
            check(abs(per - direct) < 1e-10 and br_ok and direct <= qz + 1e-14,
                  f"t={t}, m={m}: p^C(0)-p^C(m)={direct:.6e} = periodized {per:.6e}; "
                  f"bracket terms > 0; <= p^Z(0)={qz:.6e}")

    # 4. Uniform ceiling R^cyc(m) < (1-a) R_inf for all tested m; cycle beats path.
    for alpha in (0.05, 0.10, 0.25, 0.45, 0.49):
        lim = R_cyc_inf(alpha)
        ms = list(range(1, 129)) + [2 ** k for k in range(8, 21)]
        Rs = [cycle_R(m, alpha) for m in ms]
        check(max(Rs) < lim, f"alpha={alpha}: max_m R^cyc(m)={max(Rs):.8f} < R^cyc_inf={lim:.8f} (m up to 2^20)")
        Rp2 = [cycle_R(2 ** n, alpha) for n in range(0, 21)]
        check(all(Rp2[i] < Rp2[i + 1] for i in range(20)),
              f"alpha={alpha}: R^cyc(2^n) increasing toward the limit (R^cyc(2^20)={Rp2[-1]:.6f}; observation)")
        check(all(cycle_R(2 ** n, alpha) < path_R_closed(2 ** n, alpha) for n in range(0, 16)),
              f"alpha={alpha}: R^cyc(2^n) < path R(2^n) for n=0..15")
        check(lim < R_limit(alpha), f"alpha={alpha}: R^cyc_inf={lim:.6f} < R_inf={R_limit(alpha):.6f}")
    # Control at alpha = 1/2: log growth, the ceiling needs alpha < 1/2.
    Rh = [cycle_R(2 ** n, 0.5) for n in (8, 12, 16, 20)]
    check(min(Rh[i + 1] - Rh[i] for i in range(3)) > 1.0,
          f"control alpha=1/2: R^cyc(2^n) grows without bound {[round(r, 4) for r in Rh]}")

    # 5. Constants at alpha = 0.45.
    Rc = R_cyc_inf(ALPHA)
    Qf = Q_floor_cyc(ALPHA)
    Q2 = 2.25 / Rc
    Ff = Ffrac_floor_cyc(ALPHA)
    sig = sigma_AK(ALPHA)
    Fn = sig * Ff
    print(f"\n     alpha=0.45: R^cyc_inf={Rc:.6f}, (9/4)/R={Q2:.6f}, Q floor (7/3)/R={Qf:.6f}, "
          f"||F_frac|| floor={Ff:.6f}, sigma_AK={sig:.6f}, ||F|| floor={Fn:.6f}")
    check(abs(Rc - 7.284859) < 5e-7, f"R^cyc_inf(0.45) = {Rc:.9f} (~7.284859)")
    check(abs(Q2 - 0.308860) < 5e-7, f"two-corner floor (9/4)/R^cyc_inf = {Q2:.9f} (~0.308860)")
    check(abs(Qf - 0.320299) < 5e-7, f"three-corner floor (7/3)/R^cyc_inf = {Qf:.9f} (~0.320299)")
    check(abs(Ff - 0.209685) < 5e-7, f"||F_frac|| floor = {Ff:.9f} (~0.209685)")
    check(abs(Fn - 0.147106) < 5e-7, f"||F|| floor = sigma_AK*{Ff:.6f} = {Fn:.9f} (~0.147106)")
    ratio = Ff / AL_FRAC_FLOOR_045
    check(abs(ratio - (28.0 / 27.0) / (1.0 - ALPHA)) < 1e-5 and ratio > 1.88,
          f"||F_frac|| floor / AL = {ratio:.6f} = (28/27)/(1-a) > 1.88")
    check(Fn > AL_F_FLOOR_045, f"||F|| floor {Fn:.6f} > AL's {AL_F_FLOOR_045}")

    # 6. Live gasket levels 1..6.
    print("\nLive gasket checks at alpha=0.45 (levels 1–6):")
    live_chain = {}
    for n in range(1, 7):
        r = level_cycle_record(n, ALPHA)
        m = r["m"]
        Rn = cycle_R(m, ALPHA)
        g = r["G"]
        offd = [g[0, 1], g[0, 2], g[1, 2]]
        check(r["cycle_ok"] and r["nb"] == 3 * m and r["ne"] == 3 * m and r["corner_pos"] == [m, 2 * m],
              f"level {n}: outer boundary is one cycle C_{3*m} (3·2^{n} edges), corners at positions 0,{m},{2*m}")
        check(r["LminusM"] > -1e-10 and r["LH"] > -1e-9,
              f"level {n}: lambda_min(L-M)={r['LminusM']:.1e}, lambda_min(L^a-M^a)={r['LH']:.1e} >= 0")
        check(abs(r["Rcyc_live"] - Rn) < 1e-9 and Rn < Rc,
              f"level {n}: (e_a-e_b)^T(M^a)^+(e_a-e_b)={r['Rcyc_live']:.9f} = Fourier {Rn:.9f} < {Rc:.6f}")
        check(max(abs(g[i, i] - g[0, 0]) for i in range(3)) < 1e-9 and max(offd) - min(offd) < 1e-9,
              f"level {n}: corner block of (M^a)^+ is circulant (g0={g[0,0]:.6f}, g1={offd[0]:.6f})")
        three = (7.0 / 3.0) / Rn
        check(abs(r["schur_energy"] - three) < 1e-8,
              f"level {n}: min boundary-extension M^a energy = {r['schur_energy']:.9f} = (7/3)/R^cyc(m)")
        check(r["Q"] + 1e-10 >= r["Q_cyc"] >= three - 1e-9 and three >= Qf,
              f"level {n}: Q={r['Q']:.6f} >= u^T M^a u={r['Q_cyc']:.6f} >= (7/3)/R^cyc={three:.6f} >= {Qf:.6f}")
        check(three > 2.25 / path_R_closed(m, ALPHA),
              f"level {n}: cycle floor {three:.6f} > path floor {2.25/path_R_closed(m, ALPHA):.6f}")
        check(abs(r["nF_frac"] - r["Q"] * SQRT21_7) < 1e-9 and r["nF_frac"] >= Ff,
              f"level {n}: ||F_frac||={r['nF_frac']:.6f} >= {Ff:.6f}")
        check(r["nF"] + 1e-10 >= Fn, f"level {n}: ||F||={r['nF']:.6f} >= sigma_AK*floor={Fn:.6f}")
        check(r["Rgasket"] < Rn, f"level {n}: live gasket R={r['Rgasket']:.6f} < R^cyc(m)={Rn:.6f}")
        live_chain[n] = (sig * r["nF_frac"], r["Rgasket"])
    for n, q in Q_RECORDED.items():
        check(q > Qf, f"recorded level {n} energy Q={q} > cycle floor {Qf:.6f}")

    # 7. Kept Failure AM.1.
    print("\nKept Failure AM.1: the boundary-cycle ceiling does not reach Conjecture J.1")
    check(Fn < J1_VALUE, f"sigma_AK * Ffrac_floor = {Fn:.6f} < 0.70")
    env = 2.0 ** (1.0 - ALPHA) * Ff
    check(env < J1_VALUE, f"2^(1-a) * Ffrac_floor = {env:.6f} < 0.70 (absolute chain envelope)")
    check(1.3525 * Ff < J1_VALUE, f"live sigma~1.3525 * Ffrac_floor = {1.3525*Ff:.6f} < 0.70")
    for n in range(2, 7):
        check(live_chain[n][0] < J1_VALUE,
              f"level {n}: even live sigma_AK*||F_frac||_live = {live_chain[n][0]:.6f} < 0.70 "
              f"(no resistance ceiling can rescue the sigma_AK chain)")
    for n, q in Q_RECORDED.items():
        v = sig * q * SQRT21_7
        check(v < J1_VALUE, f"recorded level {n}: sigma_AK*Q*sqrt21/7 = {v:.6f} < 0.70")
    check(Rc > 3.0 * live_chain[6][1],
          f"R^cyc_inf={Rc:.4f} > 3 x live gasket R(level 6)={live_chain[6][1]:.4f} (cycle still lossy)")

    print()
    if failures:
        print(f"FAILED {len(failures)} checks:")
        for msg in failures:
            print(" -", msg)
        return 1
    print("All checks passed. Theorem AM / Kept Failure AM.1. Not thrust.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
