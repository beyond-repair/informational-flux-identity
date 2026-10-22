#!/usr/bin/env python3
"""Spectral decimation of the corner resolvent: monotone audit limit at α = 0.45.

Witness for Theorem AQ and Kept Failure AQ.1 in GASKET.md.

Notation: G_n is the level-n combinatorial gasket, L its graph Laplacian,
v = e_a - e_b, r_n(s) = ½ vᵀ(L + s)⁻¹ v (s > 0), R_α(n) = vᵀ(L^α)⁺v,
σ_n = ||F|| / ||F_frac|| (Theorem AB), β = 1 - α.

Theorem AQ.
  (1) Decimation.  Eliminating the three midpoints of every level-1 cell of
      G_{n+1} (Schur complement of L_{n+1} + s onto V_n) gives exactly
          c(s) L_n + (2e(s) + s) I - e(s) P_C,
      c = (6+s)/((5+s)(2+s)),  e = 2s/(2+s),  (2e+s)/c = s(s+5),
      P_C the corner projector.  Because the corner block is S_3-circulant
      the antisymmetric mode decouples and
          r_{n+1}(s) = Φ_s(r_n(s(s+5))),  r_0(s) = 1/(3+s),
          Φ_s(x) = (2+s)(5+s) x / ((6+s) - 2s(5+s) x).
  (2) Stieltjes representations (v ⊥ ker L):
          R_α(n) = (2 sin πα/π) ∫_0^∞ s^{-α} r_n(s) ds,
          σ_n    = ½ vᵀL^β v = (sin πβ/π) ∫_0^∞ s^{β-1} (1 - s r_n(s)) ds.
  (3) Monotonicity.  r_1 - r_0 = (s+6)/((s+3)(s²+7s+9)) > 0 and Φ_s is
      strictly increasing on [0, 1/(s(s+5))] ∋ r_n(s(s+5)), so r_n(s) is
      strictly increasing in n for every s > 0.  Hence R_α(n) strictly
      increases, σ_n strictly decreases, ||F_frac(n)|| = √21/(3R_α(n))
      strictly decreases, and ||F(n,α)|| = σ_n ||F_frac(n)|| strictly
      decreases in n, for every α ∈ (0,1).  All converge.
  (4) Uniform ceiling (α < log3/log5).  1/σ is an exact supersolution
      (Φ_s(1/(s(s+5))) = (2+s)/(s(4+s)) < 1/s) and r_m ≤ 1/σ, so
      r_n ≤ T^k(1/σ) pointwise for all n, k; a power supersolution
      B s^{-β'} on (0, σ0] (β' ∈ (1-log3/log5, 1-α)) covers s → 0.
      Certified by outward-rounded interval arithmetic (mpmath.iv) below.
  (5) At α = 0.45: R_α(n) < R̄ ≈ 2.518 for every n, σ_n > σ_lo ≈ 1.35247,
      so ||F(n,0.45)|| > L_lo ≈ 0.8204 at every level, and
      ||F(n,0.45)|| ↓ L ∈ [L_lo, L_hi] with L_hi ≈ 0.8217.  This proves
      Conjecture J.1 (decrease, existence of the limit, and L ≥ 0.70).

Kept Failure AQ.1: the level-uniform ceiling and a positive monotone limit
do not extend to α ≥ log3/log5: the power-supersolution window
(1 - log3/log5, 1 - α) is empty there, and the same recursion shows R_α(n)
growing linearly at α = log3/log5 and geometrically (ratio → 5^α/3) at
α = 0.75 (consistent with Theorems AD and K: the audit norm tends to 0).

Not thrust.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_limit import ALPHA_C, build_gasket, combinatorial_laplacian  # noqa: E402
from gasket_fractional_moment_share import S_AP, level_record  # noqa: E402

ALPHA = 0.45
J1_VALUE = 0.70
AP1_NEEDED_R = 2.859237      # Kept Failure AP.1: R ceiling needed with S_AP
AO_CEILING = 4.820047        # Theorem AO as stated
NF10 = 0.833640              # recorded ||F(10, 0.45)|| (GASKET Conjecture J.1)
FF9, FF10 = 0.620670938843, 0.616380576325  # recorded ||F_frac(9)||, ||F_frac(10)||
SQRT21_3 = math.sqrt(21.0) / 3.0

failures: list[str] = []


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


# ----------------------------------------------------------------- recursion
def phi(s, x):
    return (2 + s) * (5 + s) * x / ((6 + s) - 2 * s * (5 + s) * x)


def r_rec(n: int, s):
    """Exact r_n(s) by the decimation recursion (any number type)."""
    orb = [s]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    x = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        x = phi(orb[k], x)
    return x


def dense_r(L: np.ndarray, a: int, b: int, s: float) -> float:
    v = np.zeros(L.shape[0])
    v[a], v[b] = 1.0, -1.0
    return 0.5 * float(v @ np.linalg.solve(L + s * np.eye(L.shape[0]), v))


# ------------------------------------------------------------- certificate
CAP = mpf(10) ** 60


def iv_bounds(g: mpf):
    """Interval upper bound hi ≥ r_n(g) for every n, and lower bound lo ≤ r_n(g)
    for every n ≥ k (k = depth used).  hi = T^k(1/σ)(g), lo = T^k(r_0)(g)."""
    orb = [iv.mpf(g)]
    while not orb[-1].a > CAP:
        orb.append(orb[-1] * (orb[-1] + 5))
    k = len(orb) - 1
    hi = 1 / orb[-1]
    lo = 1 / (3 + orb[-1])
    for j in range(k - 1, -1, -1):
        s = orb[j]
        for name in ("hi", "lo"):
            x = hi if name == "hi" else lo
            den = (6 + s) - 2 * s * (5 + s) * x
            if not den.a > 0:
                raise RuntimeError("denominator not certified positive")
            y = (2 + s) * (5 + s) * x / den
            if name == "hi":
                hi = y
            else:
                lo = y
    return hi, lo, k


def certify(alpha: float, beta_p: float, sig0: float, s0_exp: int, g_exp: int, ratio: str):
    """Rigorous (interval) bounds; see Theorem AQ (4)."""
    iv.dps = 30
    mp.dps = 30
    a = iv.mpf(alpha)
    bt = 1 - a                     # β = 1 - α
    bp = iv.mpf(beta_p)            # power exponent β'
    S0 = iv.mpf(sig0)
    out = {"alpha": alpha, "beta_p": beta_p}
    out["window"] = (1 - ALPHA_C < beta_p < 1 - alpha)
    # sstar: s(s+5) = sig0
    sstar = (iv.sqrt(25 + 4 * S0) - 5) / 2
    # X ≥ r_n(sig0) ≥ r_n(s') for s' ≥ sig0 (r_n decreasing in s), every n
    X, _, _ = iv_bounds(mpf(sig0))
    Xb = iv.mpf(X.b)
    # choose B from a float scan of C2, then certify
    M = 400
    lo_e, hi_e = mpf(sstar.a), mpf(sig0)
    need = mpf(0)
    cells = []
    for i in range(M):
        u = lo_e + (hi_e - lo_e) * i / M
        w = lo_e + (hi_e - lo_e) * (i + 1) / M if i < M - 1 else hi_e
        U, W = iv.mpf(u), iv.mpf(w)
        den = (6 + U) - 2 * W * (5 + W) * Xb
        val = (2 + W) * (5 + W) * Xb / den
        cells.append((U, W, den, val))
        need = max(need, mpf((val * W ** bp).b))
    B = iv.mpf(need * mpf("1.0001"))
    ok2 = all(den.a > 0 and val.b <= (B * W ** (-bp)).a for U, W, den, val in cells)
    ss = iv.mpf(sstar.b)
    lhs = (2 + ss) * (5 + ss) ** (1 - bp) + 2 * B * ss ** (1 - bp) * (5 + ss) ** (1 - bp)
    ok1 = lhs.b < 6
    ok3 = (S0 ** bp / 3).b <= B.a
    out.update(C1=ok1, C2=ok2, C3=ok3, B=float(B.b), X=float(X.b), C1_lhs=float(lhs.b))
    # grid
    s0 = mpf(10) ** s0_exp
    G = mpf(10) ** g_exp
    q = mpf(ratio)
    assert s0 <= mpf(sig0)
    P = lambda x: iv.mpf(x) ** (1 - a) / (1 - a)   # ∫ s^{-α}
    Qf = lambda x: iv.mpf(x) ** bt / bt             # ∫ s^{β-1}
    R_up = B * iv.mpf(s0) ** (1 - a - bp) / (1 - a - bp)
    R_lo = iv.mpf(0)
    S_lo = iv.mpf(0)
    S_up = iv.mpf(s0) ** bt / bt
    g = s0
    hi_g, lo_g, kmax = iv_bounds(g)
    npts = 1
    while g < G:
        g2 = g * q
        hi2, lo2, k2 = iv_bounds(g2)
        kmax = max(kmax, k2)
        dP = P(g2) - P(g)
        dQ = Qf(g2) - Qf(g)
        R_up += hi_g * dP
        R_lo += lo2 * dP
        S_lo += (1 - iv.mpf(g2) * hi2) * dQ
        S_up += (1 - iv.mpf(g) * lo_g) * dQ
        g, hi_g, lo_g = g2, hi2, lo2
        npts += 1
    R_up += iv.mpf(g) ** (-a) / a
    S_up += 3 * iv.mpf(g) ** (bt - 1) / (1 - bt)
    cR = 2 * iv.sin(iv.pi * a) / iv.pi
    cS = iv.sin(iv.pi * bt) / iv.pi
    out["R_up"] = float((cR * R_up).b)
    out["R_lo"] = float((cR * R_lo).a)
    out["S_lo"] = float((cS * S_lo).a)
    out["S_up"] = float((cS * S_up).b)
    rt = iv.sqrt(21) / 3
    out["L_lo"] = float((cS * S_lo * rt / (cR * R_up)).a)
    out["L_hi"] = float((cS * S_up * rt / (cR * R_lo)).b)
    out["Ffrac_lo"] = float((rt / (cR * R_up)).a)
    out["kmax"] = kmax
    out["npts"] = npts
    return out


def quad_R_sigma(n: int, alpha: float):
    mp.dps = 45
    a = mpf(alpha)
    b = 1 - a
    X = mpf(10) ** 30
    pts = [0] + [mpf(10) ** k for k in range(-16, 31, 2)]
    # tail s > X: r_n(s) = 1/s + O(s^-2), 1 - s r_n(s) = 2/s + O(s^-2); error O(X^{-1-α})
    R = 2 * mp.sin(mp.pi * a) / mp.pi * (mp.quad(lambda s: s ** (-a) * r_rec(n, s), pts) + X ** (-a) / a)
    S = mp.sin(mp.pi * b) / mp.pi * (mp.quad(lambda s: s ** (b - 1) * (1 - s * r_rec(n, s)), pts) + 2 * X ** (-a) / a)
    return float(R), float(S)


def main(through: int = 6) -> int:
    assert through >= 4
    print("Theorem AQ: spectral decimation of the corner resolvent; monotone audit limit at α = 0.45")
    print("Assumption: spectral calculus (GASKET item 4); corner data (1,-1/2,0); Theorems AA, AB, AD (i), AP (1)–(3).")
    print("Not thrust.\n")

    # (1) symbolic cell Schur complement and recursion algebra
    s = sp.symbols("s", positive=True)
    Lc = sp.zeros(6, 6)
    cell_edges = [(0, 3), (0, 5), (1, 3), (1, 4), (2, 4), (2, 5), (3, 4), (4, 5), (3, 5)]
    for i, j in cell_edges:
        Lc[i, i] += 1; Lc[j, j] += 1; Lc[i, j] -= 1; Lc[j, i] -= 1
    Ms = Lc[3:, 3:] + s * sp.eye(3)
    Sch = sp.simplify(Lc[:3, :3] - Lc[:3, 3:] * Ms.inv() * Lc[3:, :3])
    c_s = (6 + s) / ((5 + s) * (2 + s))
    e_s = 2 * s / (2 + s)
    target = c_s * (3 * sp.eye(3) - sp.ones(3, 3)) + e_s * sp.eye(3)
    check(sp.simplify(Sch - target) == sp.zeros(3, 3),
          "cell Schur complement (one G_1 cell, midpoints eliminated, s on midpoints) = c(3I-J) + eI, "
          "c=(6+s)/((5+s)(2+s)), e=2s/(2+s)")
    check(sp.simplify((2 * e_s + s) / c_s - s * (s + 5)) == 0, "(2e+s)/c = s(s+5) (decimation map)")
    x = sp.symbols("x")
    Phi = sp.simplify(x / (c_s - e_s * x))
    check(sp.simplify(Phi - (2 + s) * (5 + s) * x / ((6 + s) - 2 * s * (5 + s) * x)) == 0,
          "antisymmetric corner mode: μ_{n+1}(s) = μ_n(s')/(c - e μ_n(s')) = Φ_s(μ_n(s'))")
    r0 = 1 / (3 + s)
    r1 = sp.simplify(Phi.subs(x, r0.subs(s, s * (s + 5))))
    check(sp.simplify(r1 - (5 + s) / (s ** 2 + 7 * s + 9)) == 0, "r_1 = (5+s)/(s²+7s+9)")
    check(sp.simplify(r1 - r0 - (s + 6) / ((s + 3) * (s ** 2 + 7 * s + 9))) == 0,
          "r_1 - r_0 = (s+6)/((s+3)(s²+7s+9)) > 0")
    check(sp.simplify(sp.diff(Phi, x) - (2 + s) * (5 + s) * (6 + s) / ((6 + s) - 2 * s * (5 + s) * x) ** 2) == 0,
          "Φ_s'(x) = (2+s)(5+s)(6+s)/den² > 0 where den > 0")
    check(sp.simplify(((6 + s) - 2 * s * (5 + s) / (s * (s + 5))) - (4 + s)) == 0,
          "den(1/s') = 4+s > 0: Φ_s increasing on [0, 1/s'] ∋ r_n(s')")
    check(sp.simplify(Phi.subs(x, 1 / (s * (s + 5))) - (2 + s) / (s * (4 + s))) == 0,
          "Φ_s(1/s') = (2+s)/(s(4+s)) < 1/s: 1/σ is a supersolution")
    check(sp.simplify(r0.subs(s, 0) * 2 - sp.Rational(2, 3)) == 0, "r_0(0) = 1/3 (R_eff = (2/3)(5/3)^n at s=0)")

    # (1') decimation on actual graphs
    for n in range(1, 5):
        P0, E0, C0 = build_gasket(n)
        P1, E1, C1 = build_gasket(n + 1)
        L0 = combinatorial_laplacian(len(P0), E0)
        L1 = combinatorial_laplacian(len(P1), E1)
        key = {(round(p[0], 9), round(p[1], 9)): i for i, p in enumerate(P1)}
        idx = np.array([key[(round(p[0], 9), round(p[1], 9))] for p in P0])
        rest = np.setdiff1d(np.arange(len(P1)), idx)
        PC = np.zeros((len(P0), len(P0)))
        PC[C0, C0] = 1.0
        worst = 0.0
        for sv in (0.03, 0.7, 4.0):
            A = L1 + sv * np.eye(len(P1))
            S = A[np.ix_(idx, idx)] - A[np.ix_(idx, rest)] @ np.linalg.solve(A[np.ix_(rest, rest)], A[np.ix_(rest, idx)])
            cv, ev = (6 + sv) / ((5 + sv) * (2 + sv)), 2 * sv / (2 + sv)
            T = cv * L0 + (2 * ev + sv) * np.eye(len(P0)) - ev * PC
            worst = max(worst, float(np.abs(S - T).max()))
        check(worst < 1e-10, f"level {n}->{n+1}: Schur(L_{n+1}+s → V_n) = cL_n + (2e+s)I - eP_C (max err {worst:.1e})")

    # (1'') recursion vs dense resolvent, levels 0..through
    for n in range(0, through + 1):
        P0, E0, C0 = build_gasket(n)
        L0 = combinatorial_laplacian(len(P0), E0)
        errs = [abs(dense_r(L0, int(C0[0]), int(C0[1]), sv) - r_rec(n, sv)) / r_rec(n, sv)
                for sv in (1e-4, 0.01, 0.3, 2.0, 50.0)]
        check(max(errs) < 1e-9, f"level {n}: r_n(s) from recursion = ½vᵀ(L+s)⁻¹v (rel err {max(errs):.1e})")

    # (2)–(3) Stieltjes formulas vs live, and monotonicity
    others = (0.10, 0.25, 0.49, 0.55)
    live045 = []
    for level in range(1, through + 1):
        alphas = (ALPHA,) + (others if level <= 5 else ())
        for d in level_record(level, alphas):
            Rq, Sq = quad_R_sigma(level, d["alpha"])
            ok = (abs(Rq - d["R"]) < 1e-7 * d["R"] and abs(Sq - d["sigma"]) < 1e-7
                  and abs(d["nFf"] - SQRT21_3 / d["R"]) < 1e-9 and abs(d["nF"] - d["sigma"] * d["nFf"]) < 1e-9)
            check(ok, f"α={d['alpha']}, level {level}: R_α = (2sinπα/π)∫s^-α r_n = {Rq:.6f}, "
                      f"σ_n = (sinπβ/π)∫s^(β-1)(1-s r_n) = {Sq:.6f}, ||F_frac|| = √21/(3R), ||F|| = σ||F_frac||")
            if d["alpha"] == ALPHA:
                live045.append(d)
    Rs = [d["R"] for d in live045]
    Ss = [d["sigma"] for d in live045]
    Fs = [d["nF"] for d in live045]
    check(all(Rs[i] < Rs[i + 1] for i in range(len(Rs) - 1)), "live R_0.45(n) strictly increasing, levels 1-%d" % through)
    check(all(Ss[i] > Ss[i + 1] for i in range(len(Ss) - 1)), "live σ_n strictly decreasing")
    check(all(Fs[i] > Fs[i + 1] for i in range(len(Fs) - 1)), "live ||F(n,0.45)|| strictly decreasing")
    check(abs(Fs[1] - 1.165909) < 5e-6, f"level-2 lock ||F|| = {Fs[1]:.6f}")
    sgrid = [mpf(10) ** k for k in range(-12, 7)]
    mp.dps = 30
    mono = all(r_rec(n + 1, sv) - r_rec(n, sv) >= -mpf(10) ** -25 * r_rec(n, sv) for n in range(0, 40) for sv in sgrid)
    check(mono, "recursion: r_n(s) ≤ r_{n+1}(s) for n < 40 on s ∈ {1e-12..1e6} (to 1e-25 rel.; strictness is the algebra above)")
    R9, S9 = quad_R_sigma(9, ALPHA)
    R10, S10 = quad_R_sigma(10, ALPHA)
    check(abs(SQRT21_3 / R9 - FF9) < 1e-8 and abs(SQRT21_3 / R10 - FF10) < 1e-8,
          f"recursion reproduces recorded matrix-free ||F_frac(9)||, ||F_frac(10)|| = {SQRT21_3/R9:.9f}, {SQRT21_3/R10:.9f}")
    check(abs(S10 * SQRT21_3 / R10 - NF10) < 1e-6,
          f"recursion reproduces recorded ||F(10,0.45)|| = {S10*SQRT21_3/R10:.6f}")

    # (4)–(5) certificate at α = 0.45
    print("\nCertificate (mpmath.iv, outward rounding) at α = 0.45 ...")
    c = certify(ALPHA, 0.35, 0.01, -20, 12, "1.004")
    print(f"  β'=0.35, σ0=0.01, B={c['B']:.6f}, X=T^k(1/σ)(σ0)={c['X']:.6f}, C1 lhs={c['C1_lhs']:.6f} < 6, "
          f"grid points {c['npts']}, depth k ≤ {c['kmax']}")
    print(f"  R̄ = {c['R_up']:.6f} (all n);  R_lo = {c['R_lo']:.6f} (n ≥ {c['kmax']})")
    print(f"  σ_lo = {c['S_lo']:.9f} (all n);  σ_up = {c['S_up']:.9f} (n ≥ {c['kmax']})")
    print(f"  ||F_frac(n)|| ≥ {c['Ffrac_lo']:.6f};  L ∈ [{c['L_lo']:.6f}, {c['L_hi']:.6f}]")
    check(c["window"] and c["C1"] and c["C2"] and c["C3"],
          "power supersolution B s^-β' on (0,σ0]: C1 (s ≤ s*), C2 (s* < s ≤ σ0), C3 (base r_0 ≤ Bs^-β')")
    check(c["R_lo"] <= c["R_up"] and c["S_lo"] <= c["S_up"], "certified brackets are consistent")
    check(c["R_up"] < 2.52 and c["R_up"] < AP1_NEEDED_R and c["R_up"] < AO_CEILING,
          f"R_0.45(n) ≤ R̄ = {c['R_up']:.6f} < 2.52 < 2.859237 (AP.1 target) < 4.820047 (AO) for every n")
    check(c["S_lo"] > S_AP(ALPHA), f"σ_n ≥ σ_lo = {c['S_lo']:.6f} > S_AP = {S_AP(ALPHA):.6f} for every n")
    check(c["L_lo"] > 0.819 and c["L_lo"] > J1_VALUE and c["L_hi"] < 0.823,
          f"||F(n,0.45)|| ↓ L with {c['L_lo']:.6f} ≤ L ≤ {c['L_hi']:.6f}; L ≥ 0.70 (Conjecture J.1)")
    check(S_AP(ALPHA) * c["Ffrac_lo"] > J1_VALUE,
          f"also S_AP·√21/(3R̄) = {S_AP(ALPHA)*c['Ffrac_lo']:.6f} ≥ 0.70 (AP.1's open step closed)")
    check(all(d["R"] <= c["R_up"] and d["sigma"] >= c["S_lo"] and d["nF"] >= c["L_lo"] for d in live045)
          and R10 <= c["R_up"] and NF10 >= c["L_lo"],
          "live levels 1-%d and recorded level 10 respect R ≤ R̄, σ ≥ σ_lo, ||F|| ≥ L_lo" % through)
    check(Fs[-1] >= c["L_hi"] - 1e-12 and NF10 >= c["L_hi"] - 1e-6,
          "||F(6)||, ||F(10)|| lie above the certified limit bracket (consistent with decrease)")

    # other α
    print()
    for al, sexp in ((0.10, -20), (0.25, -20), (0.49, -30), (0.55, -40)):
        co = certify(al, 0.35, 0.01, sexp, 12, "1.02")
        lv = [d for lvl in range(1, 6) for d in level_record(lvl, (al,))]
        ok = (co["window"] and co["C1"] and co["C2"] and co["C3"]
              and all(d["R"] <= co["R_up"] and d["sigma"] >= co["S_lo"] and d["nF"] >= co["L_lo"] for d in lv)
              and co["L_lo"] <= co["L_hi"] and co["L_lo"] > 0)
        check(ok, f"α={al}: certified R̄={co['R_up']:.6f}, σ_lo={co['S_lo']:.6f}, L∈[{co['L_lo']:.6f},{co['L_hi']:.6f}]; "
                  "live levels 1-5 inside")

    # Kept Failure AQ.1
    print("\nKept Failure AQ.1: no uniform ceiling for α ≥ log3/log5")
    for al in (ALPHA_C, 0.75):
        check(not (1 - ALPHA_C < 0.35 < 1 - al) and 1 - ALPHA_C >= 1 - al,
              f"α={al:.6f}: power window (1-log3/log5, 1-α) is empty")
    ns = [8, 12, 16, 20, 24]
    Rg = [quad_R_sigma(n, ALPHA_C)[0] for n in ns]
    inc = [(Rg[i + 1] - Rg[i]) / 4 for i in range(len(ns) - 1)]
    print(f"  α=log3/log5: R(n) at n={ns}: {[round(x, 4) for x in Rg]}; per-level increments {[round(x, 4) for x in inc]}")
    check(min(inc) > 0.5 * max(inc) and min(inc) > 0.05, "α=log3/log5: R_α(n) grows linearly (Theorem AD's Θ(n))")
    R75 = [quad_R_sigma(n, 0.75)[0] for n in (16, 20, 24)]
    rat = (R75[2] / R75[1]) ** 0.25
    print(f"  α=0.75: R(16,20,24) = {[round(x, 3) for x in R75]}; per-level ratio {rat:.5f} vs 5^α/3 = {5**0.75/3:.5f}")
    check(abs(rat - 5 ** 0.75 / 3) < 0.02 and rat > 1.05, "α=0.75: R_α(n) grows geometrically (ratio → 5^α/3 > 1)")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        return 1
    print("All Theorem AQ / Kept Failure AQ.1 checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    through = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    raise SystemExit(main(through))
