#!/usr/bin/env python3
"""Sharp path-resistance ceiling via Hermite–Hadamard (Theorem AL).

Witness for Theorem AL and Kept Failure AL.1 in GASKET.md.

Theorem AA bounded the boundary-path resistance by the crude constant
  C(α) = 4^{-α}(8 + (4/π) B(1/2−α, 3/2)),
separating the j=0 Riemann term with sin x ≥ 2x/π. The same cosine sum is
a midpoint Riemann sum for the integrand
  f(θ) = (2 sin θ)^{-2α} cos² θ
on (0, π/2). For α ≥ α_* := (14−8√2)/34 ≈ 0.079009, f is convex on (0, π/2)
(second-derivative minimum nonnegative), so Hermite–Hadamard gives
  R_α(m) ≤ R_∞(α) := 4^{-α}(4/π) B(1/2−α, 3/2) < C(α)
for every path length m ≥ 1. Löwner–Heinz on the gasket side then sharpens
AA's floors: at α = 0.45,
  Q ≥ 0.169873,  ||F_frac|| ≥ 0.111208
(about 1.324× AA's 0.084015). With Theorem AK's share floor,
  ||F(n,0.45)|| ≥ σ_AK · 0.111208 ≈ 0.078019.

Kept Failure AL.1: replacing C by R_∞ in the AA–AK chain does not reach
Conjecture J.1's 0.70 (≈0.078; even 2^{1−α}·0.111208 ≈ 0.163). The path
comparison cannot: live gasket R_α ≈ 2.35 ≪ R_∞ ≈ 13.25.

Not thrust. This is a different comparison from the one-step neighbor-Harnack
chain closed in Kept Failure AK.1.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_energy_floor import (  # noqa: E402
    C_bound,
    R_limit,
    beta,
    path_R_closed,
    path_R_dense,
    level_record,
)
from gasket_fractional_moment_jensen import sigma_AK  # noqa: E402

ALPHA = 0.45
ALPHA_STAR = (14.0 - 8.0 * math.sqrt(2.0)) / 34.0
AA_FRAC_FLOOR_045 = 0.084015
J1_VALUE = 0.70
SQRT21_7 = math.sqrt(21.0) / 7.0

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def d2_over_g(x: float, b: float) -> float:
    """f''/positive for f = sin^{-b} cos^2, in the coordinate x = sin^2 θ."""
    return b * (b + 1.0) / x + 2.0 / (1.0 - x) + 4.0 * b - 4.0 - b * b


def min_d2_over_g(b: float) -> float:
    """Exact critical value: 5b - 2 + 2√(2b(b+1))."""
    return 5.0 * b - 2.0 + 2.0 * math.sqrt(2.0 * b * (b + 1.0))


def Q_floor(alpha: float) -> float:
    return 2.25 / R_limit(alpha)


def Ffrac_floor(alpha: float) -> float:
    return Q_floor(alpha) * SQRT21_7


def main() -> int:
    print("Theorem AL: Hermite–Hadamard sharpens the path-resistance ceiling of AA")
    print("Assumption: spectral calculus (GASKET item 4); Löwner–Heinz; corner data")
    print("(1,-1/2,0); Theorems AA, AK. Not thrust.\n")

    # --- α_* and convexity threshold ---
    b_star = 2.0 * ALPHA_STAR
    check(
        abs(b_star - (14.0 - 8.0 * math.sqrt(2.0)) / 17.0) < 1e-14,
        f"alpha_* = (14-8√2)/34 = {ALPHA_STAR:.12f} (b_* = 2 alpha_* = {b_star:.12f})",
    )
    check(
        abs(min_d2_over_g(b_star)) < 1e-12,
        f"min (f''/pos) at b_* = {min_d2_over_g(b_star):.3e} = 0 (threshold)",
    )
    check(
        min_d2_over_g(2.0 * ALPHA) > 6.0,
        f"at alpha=0.45: min (f''/pos) = {min_d2_over_g(0.9):.6f} > 0 (convex)",
    )
    # Grid check of analytic critical formula vs brute force
    for b in (0.2, 0.5, 0.9, 0.99):
        s = math.sqrt(b * (b + 1.0))
        x_crit = s / (s + math.sqrt(2.0))
        analytic = min_d2_over_g(b)
        at_crit = d2_over_g(x_crit, b)
        xs = np.linspace(1e-6, 1.0 - 1e-6, 20000)
        brute = float(np.min(d2_over_g(xs, b)))
        check(
            abs(analytic - at_crit) < 1e-12 and brute >= analytic - 1e-4,
            f"b={b}: analytic min={analytic:.8f} matches crit and brute ({brute:.8f})",
        )
    # Below α_*: convexity fails (so AL's Hermite proof does not cover it)
    check(
        min_d2_over_g(0.1) < 0,
        f"Kept note: at b=0.1 (alpha=0.05) min (f''/pos)={min_d2_over_g(0.1):.4f} < 0 "
        f"(convexity fails below alpha_*)",
    )

    # --- R_∞ < C and numerical floors at 0.45 ---
    Rinf = R_limit(ALPHA)
    C = C_bound(ALPHA)
    Qf = Q_floor(ALPHA)
    Ff = Ffrac_floor(ALPHA)
    sig = sigma_AK(ALPHA)
    Fn = sig * Ff
    check(Rinf < C - 1.0, f"R_inf={Rinf:.6f} < C(0.45)={C:.6f} (factor {C/Rinf:.6f})")
    check(
        abs(Rinf - 4.0 ** (-ALPHA) * (4.0 / math.pi) * beta(0.5 - ALPHA, 1.5)) < 1e-12,
        f"R_inf formula matches 4^{{-a}}(4/pi)B(1/2-a,3/2) = {Rinf:.12f}",
    )
    check(abs(Qf - 0.169873) < 5e-7, f"Q floor = (9/4)/R_inf = {Qf:.12f} (~0.169873)")
    check(abs(Ff - 0.111208) < 5e-7, f"||F_frac|| floor = {Ff:.12f} (~0.111208)")
    check(
        Ff / AA_FRAC_FLOOR_045 > 1.32,
        f"||F_frac|| floor / AA = {Ff/AA_FRAC_FLOOR_045:.6f} > 1.32",
    )
    check(abs(Fn - 0.078019) < 5e-7, f"||F|| floor via sigma_AK = {Fn:.12f} (~0.078019)")

    # --- Path resistance: R(m) ≤ R_∞ for alpha ≥ alpha_* ---
    alphas = [0.10, 0.25, 0.45, 0.49]
    for a in alphas:
        check(a >= ALPHA_STAR - 1e-15, f"alpha={a} >= alpha_* (Hermite range)")
        lim = R_limit(a)
        worst_ratio = 0.0
        for m in list(range(1, 65)) + [2**k for k in range(6, 12)]:
            r = path_R_closed(m, a)
            worst_ratio = max(worst_ratio, r / lim)
            check(r <= lim + 1e-9, f"alpha={a}: R({m})={r:.8f} <= R_inf={lim:.8f}")
        check(
            worst_ratio <= 1.0 + 1e-9,
            f"alpha={a}: max_m R(m)/R_inf = {worst_ratio:.10f} <= 1",
        )
        # Monotone in m on a short range (observation used only as check, not claimed)
        Rs = [path_R_closed(m, a) for m in range(1, 33)]
        check(
            all(Rs[i + 1] + 1e-14 >= Rs[i] for i in range(len(Rs) - 1)),
            f"alpha={a}: R(m) increasing for m=1..32 (witness observation)",
        )

    # Closed form vs dense for small m
    for m in (1, 2, 3, 4, 8, 16):
        for a in (0.10, 0.45, 0.49):
            rc, rd = path_R_closed(m, a), path_R_dense(m, a)
            check(
                abs(rc - rd) < 1e-9,
                f"alpha={a}, m={m}: closed R={rc:.10f} matches dense {rd:.10f}",
            )

    # Midpoint / Hermite structure at odd m (exact cover of [0, pi/2])
    a = ALPHA
    for m in (1, 3, 5, 9, 17, 33):
        # m odd => (m+1)/2 intervals of width Delta = pi/(m+1) cover [0, pi/2]
        check(
            m % 2 == 1,
            f"m={m} odd: midpoint sum covers [0, pi/2] exactly",
        )
        N = m + 1
        Delta = math.pi / N
        n_terms = (m + 1) // 2
        check(
            abs(n_terms * Delta - math.pi / 2) < 1e-14,
            f"m={m}: {n_terms} * Delta = pi/2 (exact cover)",
        )

    # --- Live gasket: AA route with new floor ---
    from gasket_fractional_wpq_harnack import level_data  # noqa: E402

    print("\nLive gasket checks at alpha=0.45 (levels 1–6):")
    for n in range(1, 7):
        rec = level_record(n, ALPHA)
        d = level_data(n, ALPHA)
        Q = rec["Q"]
        nF_frac = rec["nF_frac"]
        check(Q + 1e-10 >= Qf, f"level {n}: Q={Q:.6f} >= Q_floor={Qf:.6f}")
        check(
            nF_frac + 1e-10 >= Ff,
            f"level {n}: ||F_frac||={nF_frac:.6f} >= {Ff:.6f}",
        )
        check(
            d["nF"] + 1e-10 >= Fn,
            f"level {n}: ||F||={d['nF']:.6f} >= sigma_AK*Ffrac_floor={Fn:.6f}",
        )
        check(
            d["nF"] + 1e-9 >= sig * nF_frac - 1e-6,
            f"level {n}: ||F||={d['nF']:.6f} >= sigma_AK*||F_frac||={sig*nF_frac:.6f}",
        )
        # Path energy still dominates the floor, and live beats path
        m = 2**n
        Rpath = path_R_closed(m, ALPHA)
        Q_path = 2.25 / Rpath
        check(
            Q + 1e-9 >= Q_path >= Qf - 1e-9,
            f"level {n}: Q >= path energy {Q_path:.6f} >= Q_floor",
        )

    # Other alphas in Hermite range on level 4
    for a in (0.10, 0.25, 0.49):
        rec = level_record(4, a)
        Qf_a = Q_floor(a)
        Ff_a = Ffrac_floor(a)
        check(rec["Q"] >= Qf_a - 1e-9, f"alpha={a} level 4: Q={rec['Q']:.6f} >= {Qf_a:.6f}")
        check(
            rec["nF_frac"] >= Ff_a - 1e-9,
            f"alpha={a} level 4: ||F_frac||={rec['nF_frac']:.6f} >= {Ff_a:.6f}",
        )

    # --- Kept Failure AL.1 ---
    print("\nKept Failure AL.1: sharp path ceiling does not reach Conjecture J.1")
    chain = Fn
    envelope = (2.0 ** (1.0 - ALPHA)) * Ff
    live_sigma_proxy = 1.3525 * Ff
    check(
        chain < J1_VALUE,
        f"sigma_AK * Ffrac_floor = {chain:.6f} < 0.70",
    )
    check(
        envelope < J1_VALUE,
        f"2^{{1-a}} * Ffrac_floor = {envelope:.6f} < 0.70 (absolute chain envelope)",
    )
    check(
        live_sigma_proxy < J1_VALUE,
        f"live sigma~1.3525 * Ffrac_floor = {live_sigma_proxy:.6f} < 0.70",
    )
    # Path cannot reach the live gasket resistance
    check(
        Rinf > 5.0 * 2.5,
        f"R_inf={Rinf:.4f} >> live gasket R_alpha~2.35 (path comparison is lossy)",
    )

    print()
    if failures:
        print(f"FAILED {len(failures)} checks:")
        for msg in failures:
            print(" -", msg)
        return 1
    print("All checks passed. Theorem AL / Kept Failure AL.1. Not thrust.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
