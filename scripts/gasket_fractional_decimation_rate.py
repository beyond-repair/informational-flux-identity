#!/usr/bin/env python3
"""Exact refinement rates of the corner resolvent integrals (Theorem AR).

Witness for Theorem AR and Kept Failure AR.1 in GASKET.md.

Notation as in scripts/gasket_fractional_decimation.py (Theorem AQ):
r_n(s) = ½ vᵀ(L_n + s)⁻¹ v, r_{n+1}(s) = Φ_s(r_n(p(s))), p(s) = s(s+5),
q = p⁻¹ on [0,∞), q(u) = 2u/(5+√(25+4u)), δ_n = r_n - r_{n-1},
ρ = 5^α/3, μ = 5^{α-1}/3, γ = log3/log5.

Theorem AR (every α ∈ (0,1)).
  (1) δ_{n+1}(q^n u) = P_n(u) δ_1(u), P_n(u) = Π_{j=1}^n D_j(q^j u), where
      D_j(s) is the secant slope of Φ_s between r_{j-1}(p s) and r_j(p s):
          D_j(s) = (2+s)(5+s)(6+s) / (den(r_{j-1}(ps)) den(r_j(ps))),
          den(x) = (6+s) - 2 s(5+s) x.
      D_j ≥ (2+s)(5+s)/(6+s) ≥ 5/3, and with the AQ(4) ceiling
      r ≤ B s^{-β'} (β' = 0.35, certified, α-free) (3/5)D_j(s) = 1 + O(s^{0.65}).
  (2) R_α(n+1) - R_α(n) = c_R(α) ρ^n (1 + o(1)),
      σ_n - σ_{n+1}     = c_S(α) μ^n (1 + o(1)),
      c_R = (2 sin πα/π) K_R, c_S = (sin πβ/π) K_S,
      K_R = ∫ Ψ^{-α} Ψ' P_∞ δ_1 du,  K_S = ∫ Ψ^{1-α} Ψ' P_∞ δ_1 du,
      Ψ = lim 5^n q^n (Koenigs function of q at its fixed point 0),
      P_∞ = lim (3/5)^n P_n.  0 < K_R, K_S < ∞.
  (3) σ_∞ ≥ 2·4^{-α}.
  Corollaries: α < γ: ||F(n)|| - L = Θ(ρ^n) with ratio exactly 5^α/3;
  α = γ: R_γ(n) = c_R(γ) n + o(n) (sharpens Theorem AD's Θ(n));
  α > γ: ||F(n,α)|| = Θ((3·5^{-α})^n).

Kept Failure AR.1. Without the AQ(4) power ceiling, the universal bound
den ≥ 4+s gives only (3/5)D_j ≤ (3/5)(60/16) = 2.25 near s = 0, which yields
no rate at all; the rate proof needs the certified ceiling.

Not thrust. Does not touch Stage 1 or Stage 2.
"""
from __future__ import annotations

import math
import sys

import numpy as np
import sympy as sp
from mpmath import mp, mpf

failures: list[str] = []
GAMMA = math.log(3) / math.log(5)


def check(ok: bool, msg: str) -> None:
    print(("  PASS " if ok else "  FAIL ") + msg)
    if not ok:
        failures.append(msg)


def r_rec_np(n: int, s: np.ndarray) -> np.ndarray:
    """Exact recursion in float64; for orbit points > 1e40 use r(o) = 1/(3+o) (relative error O(1/o))."""
    orb = [s]
    for _ in range(n):
        orb.append(np.minimum(orb[-1] * (orb[-1] + 5), 1e90))
    x = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        o = orb[k]
        x = np.where(o > 1e40, 1 / (3 + o), (2 + o) * (5 + o) * x / ((6 + o) - 2 * o * (5 + o) * x))
    return x


def r_rec_mp(n: int, s):
    orb = [s]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    x = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        o = orb[k]
        x = (2 + o) * (5 + o) * x / ((6 + o) - 2 * o * (5 + o) * x)
    return x


def qf(u):
    return 2 * u / (5 + np.sqrt(25 + 4 * u))


def main() -> int:
    print("Theorem AR: exact refinement rates of R_α(n) and σ_n\n")

    # (1) symbolic: secant slope of the Möbius map Φ_s, and its floor
    s, x, y = sp.symbols("s x y", positive=True)
    phi = lambda z: (2 + s) * (5 + s) * z / ((6 + s) - 2 * s * (5 + s) * z)
    den = lambda z: (6 + s) - 2 * s * (5 + s) * z
    sec = sp.simplify((phi(x) - phi(y)) / (x - y) - (2 + s) * (5 + s) * (6 + s) / (den(x) * den(y)))
    check(sec == 0, "secant slope of Φ_s = (2+s)(5+s)(6+s)/(den(x)den(y))")
    floor_expr = sp.simplify((2 + s) * (5 + s) / (6 + s) - sp.Rational(5, 3))
    check(sp.simplify(floor_expr - s * (3 * s + 16) / (3 * (6 + s))) == 0, "(2+s)(5+s)/(6+s) - 5/3 = s(3s+16)/(3(6+s))")

    # (1) numeric: factorization δ_{n+1}(q^n u) = P_n(u) δ_1(u)
    mp.dps = 50
    worst = mpf(0)
    for u in [mpf("0.003"), mpf("0.7"), mpf(4), mpf(60)]:
        for n in [1, 3, 6]:
            qs = [u]
            for _ in range(n):
                qs.append(2 * qs[-1] / (5 + mp.sqrt(25 + 4 * qs[-1])))
            P = mpf(1)
            for j in range(1, n + 1):
                sj = qs[j]
                xa, xb = r_rec_mp(j - 1, qs[j - 1]), r_rec_mp(j, qs[j - 1])
                dd = lambda z: (6 + sj) - 2 * sj * (5 + sj) * z
                P *= (2 + sj) * (5 + sj) * (6 + sj) / (dd(xa) * dd(xb))
            d1 = (u + 6) / ((u + 3) * (u ** 2 + 7 * u + 9))
            lhs = r_rec_mp(n + 1, qs[n]) - r_rec_mp(n, qs[n])
            worst = max(worst, abs(lhs / (P * d1) - 1))
    check(worst < mpf(10) ** -30, f"δ_(n+1)(q^n u) = P_n(u) δ_1(u) (max rel err {float(worst):.1e})")

    # (2) limit constants vs scaled increments
    t = np.linspace(np.log(1e-60), np.log(1e8), 120001)
    sg = np.exp(t)
    h = t[1] - t[0]
    N = 45
    with np.errstate(all="ignore"):
        qs = [sg]
        for _ in range(N):
            qs.append(qf(qs[-1]))
        P = np.ones_like(sg)
        Pmax_dev = 0.0
        for j in range(1, N + 1):
            sj = qs[j]
            xa, xb = r_rec_np(j - 1, qs[j - 1]), r_rec_np(j, qs[j - 1])
            D = (2 + sj) * (5 + sj) * (6 + sj) / (((6 + sj) - 2 * sj * (5 + sj) * xa) * ((6 + sj) - 2 * sj * (5 + sj) * xb))
            if j >= 30:
                Pmax_dev = max(Pmax_dev, float(np.max(np.abs(0.6 * D - 1))))
            check_min = float(np.min(0.6 * D))
            if check_min < 1 - 1e-12:
                failures.append(f"(3/5)D_{j} < 1 somewhere ({check_min})")
            P *= 0.6 * D
        Psi = 5.0 ** N * qs[N]
        dPsi = np.ones_like(sg)
        for j in range(N):
            dPsi *= 5 / np.sqrt(25 + 4 * qs[j])
        d1 = (sg + 6) / ((sg + 3) * (sg ** 2 + 7 * sg + 9))
        check(Pmax_dev < 1e-9, f"(3/5)D_j → 1 along the orbit (max |·-1| for j ≥ 30: {Pmax_dev:.1e})")
        check(not any("D_" in f for f in failures), "(3/5)D_j ≥ 1 on the whole grid (j = 1..45)")

        print("\n  α        ρ=5^α/3    c_R(α)       scaled ΔR (n=16..20)            c_S(α)      scaled Δσ (n=8..10)")
        cR_gamma = None
        for a in [0.25, 0.45, 0.60, GAMMA, 0.75]:
            KR = np.sum(sg * Psi ** (-a) * dPsi * P * d1) * h
            KS = np.sum(sg * Psi ** (1 - a) * dPsi * P * d1) * h
            cR = 2 * math.sin(math.pi * a) / math.pi * KR
            cS = math.sin(math.pi * (1 - a)) / math.pi * KS
            rho, mu = 5 ** a / 3, 5 ** (a - 1) / 3
            # tails beyond X = 1e8: r_n(s) = 1/s + O(s^-2), 1 - s r_n(s) = 2/s + O(s^-2) (n-independent to that order)
            tR = 2 * math.sin(math.pi * a) / math.pi * 1e8 ** (-a) / a
            tS = math.sin(math.pi * (1 - a)) / math.pi * 2 * 1e8 ** (-a) / a
            Rf = lambda n: 2 * math.sin(math.pi * a) / math.pi * np.sum(sg ** (1 - a) * r_rec_np(n, sg)) * h + tR
            Sf = lambda n: math.sin(math.pi * (1 - a)) / math.pi * np.sum(sg ** (1 - a) * (1 - sg * r_rec_np(n, sg))) * h + tS
            R = [Rf(n) for n in range(0, 26)]
            S = [Sf(n) for n in range(0, 26)]
            sR = [(R[n + 1] - R[n]) / rho ** n for n in range(16, 21)]
            sS = [(S[n] - S[n + 1]) / mu ** n for n in range(8, 11)]
            print(f"  {a:.6f} {rho:.6f}  {cR:.8f}  {[round(float(v), 8) for v in sR]}  {cS:.6f}  {[round(float(v), 6) for v in sS]}")
            check(abs(sR[-1] / cR - 1) < 1e-6, f"α={a:.4f}: (R(n+1)-R(n))/ρ^n → c_R (rel {abs(sR[-1]/cR-1):.1e})")
            check(abs(sS[-1] / cS - 1) < 1e-3, f"α={a:.4f}: (σ_n-σ_(n+1))/μ^n → c_S (rel {abs(sS[-1]/cS-1):.1e})")
            check(S[11] >= 2 * 4 ** (-a), f"α={a:.4f}: σ_11 = {S[11]:.6f} ≥ 2·4^-α = {2*4**(-a):.6f}")
            if abs(a - GAMMA) < 1e-12:
                cR_gamma = cR
            if abs(a - 0.45) < 1e-12:
                # ||F|| - L ratio → ρ (σ gap is μ^n = (ρ/5)^n, negligible)
                F = [math.sqrt(21) / 3 * S[n] / R[n] for n in range(0, 26)]
                Lim = math.sqrt(21) / 3 * Sf(60) / Rf(60)
                rat = (F[22] - Lim) / (F[21] - Lim)
                print(f"  α=0.45: ||F(10)|| = {F[10]:.6f} (recorded 0.833640), L ≈ {Lim:.6f}")
                check(abs(F[10] - 0.833640) < 5e-6, "α=0.45: recursion reproduces recorded ||F(10)||")
                check(abs(rat - rho) < 2e-3, f"α=0.45: (||F(22)||-L)/(||F(21)||-L) = {rat:.5f} vs ρ = {rho:.5f}")
                check(0.819255 <= Lim <= 0.822959, "α=0.45: recursion limit sits inside Theorem AQ's bracket")
        print(f"\n  α=γ: per-level increment c_R(γ) = {cR_gamma:.6f} (Kept Failure AQ.1 recorded ≈ 0.5172 on n=8–24);"
              f" Theorem AD floor c_AD = {math.sin(math.pi*GAMMA)/math.pi*3/25*math.log(5):.6f}")
        check(abs(cR_gamma - 0.5172) < 2e-3, "α=γ: R_γ(n) = c_R(γ) n + o(n), c_R(γ) ≈ 0.51717 matches AQ.1's measured slope")

    # Kept Failure AR.1
    naive = 0.6 * (2 * 5 * 6) / 16
    print(f"\n  Kept Failure AR.1: naive (3/5)·sup D at s→0 with den ≥ 4+s only: {naive:.4f} > 1")
    check(naive > 1, "AR.1: without the AQ(4) power ceiling the product bound gives no rate")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        for f in failures:
            print("   ", f)
        return 1
    print("All Theorem AR / Kept Failure AR.1 checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
