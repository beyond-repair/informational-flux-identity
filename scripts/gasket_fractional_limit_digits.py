#!/usr/bin/env python3
"""Certified digits of the audit limit via complete monotonicity (Theorem AS).

Witness for Theorem AS and Kept Failure AS.1 in GASKET.md.

Notation as in scripts/gasket_fractional_decimation.py (Theorem AQ):
r_n(s) = ½ vᵀ(L_n + s)⁻¹ v, r_{n+1}(s) = Φ_s(r_n(s(s+5))), r_0 = 1/(3+s),
w_n(s) = 1 - s r_n(s), r_∞ = lim r_n (increasing), w_∞ = 1 - s r_∞.
By AQ (2) and monotone convergence (sin πβ = sin πα):
    R_∞(α) = lim R_α(n) = (2 sin πα/π) J_r,   σ_∞(α) = lim σ_n = (sin πα/π) J_w,
    J_r = ∫_0^∞ s^{-α} r_∞ ds,  J_w = ∫_0^∞ s^{-α} w_∞ ds,
    L(α) = lim ||F(n,α)|| = (√21/3) σ_∞/R_∞ = (√21/6) J_w/J_r.

Theorem AS.
  (1) r_n and w_n are completely monotone (finite sums c_i/(λ_i+s) and
      c_i λ_i/(λ_i+s), c_i ≥ 0), so every divided difference of order k has
      sign (-1)^k; this passes to the pointwise limits r_∞, w_∞.
  (2) Interpolation sandwich.  On a cell (x_c, x_{c+1}) of a node set with b
      nodes ≤ x_c and u nodes ≥ x_{c+1}, f - p = f[nodes, x]·Π(x - x_i) has sign
      (-1)^b for f ∈ {r_∞, w_∞}: the interpolant is a floor for even b and a
      ceiling for odd b.  Used with b = u = 10 (floor) and b = u = 11 (ceiling).
  (3) Node values: r_∞(x) ∈ [T^k r_0(x), T^k(1/σ)(x)] (AQ (4)); head (0, s0]
      by AQ's power ceiling B s^{-0.35} (α-free, re-certified here) and
      monotonicity; tail [X, ∞) by r ∈ [1/s - 2/s², 1/s], w ∈ [2/s - 6/s², 2/s]
      (½vᵀLv = 2, ½vᵀL²v = 6 for n ≥ 2).
  (4) Outward-rounded mpmath.iv sums give J_r, J_w, R_∞, σ_∞, L to ~1e-12 at
      α ∈ {0.10, 0.25, 0.45, 0.49, 0.55}; at α = 0.45 L ∈ [0.82113862713, 0.82113862714]
      (12 digits; AQ's bracket was [0.819255, 0.822959]).
  (5) Conjecture J.2 (||F(n,0.10)|| ≥ 1.410 for every n) is PROVED:
      ||F(n,0.10)|| decreases (AQ (3)) to L(0.10) > 1.410.

Kept Failure AS.1.  Claim: with the exact gap ratio ρ = 5^α/3 of Theorem AR,
one-step geometric extrapolation of the recorded levels 9, 10 recovers L(0.45)
to 1e-4.  False: it gives 0.820866..., off by 2.7e-4 from the certified L
(Aitken on 8, 9, 10 is off by 3.9e-4); at level 10 the subleading corrections
to AR's leading term are still ~2% of the gap.

Not thrust. Does not touch Stage 1 or Stage 2.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_limit import build_gasket, combinatorial_laplacian  # noqa: E402

import numpy as np  # noqa: E402

failures: list[str] = []
SQRT21_3 = math.sqrt(21.0) / 3.0
NF10, NF9 = 0.833640, 0.839442      # recorded ||F(10,0.45)||, recursion ||F(9,0.45)|| (AQ table / AR witness)
AQ_BRACKET = (0.819255, 0.822959)
J2_LEVELS = (1.50450860, 1.44384308, 1.42364897, 1.41577938, 1.41266163, 1.41143248)  # GASKET Conjecture J.2

Q_NUM, Q_DEN = 17, 16          # node ratio q = 17/16 (exact rational)
EX0 = 250                      # s0 = 2^-250 (head cut)
X_TAIL = mpf(10) ** 25         # tail cut
FLOOR_BU = (10, 10)            # (b, u) for the floor interpolant (b even)
CEIL_BU = (11, 11)             # (b, u) for the ceiling interpolant (b odd)
CAP = mpf(10) ** 30


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def node_bounds(x):
    """Enclosures lo ≤ r_∞(x) ≤ hi (x an interval): hi = T^k(1/σ)(x), lo = T^k r_0(x)."""
    orb = [x]
    while not orb[-1].a > CAP:
        orb.append(orb[-1] * (orb[-1] + 5))
    hi = 1 / orb[-1]
    lo = 1 / (3 + orb[-1])
    for j in range(len(orb) - 2, -1, -1):
        s = orb[j]
        A = (2 + s) * (5 + s)
        C = 2 * s * (5 + s)
        dh, dl = (6 + s) - C * hi, (6 + s) - C * lo
        if not (dh.a > 0 and dl.a > 0):
            raise RuntimeError("denominator not certified positive")
        hi, lo = A * hi / dh, A * lo / dl
    return lo, hi


def lagrange_weights(q, ms, a):
    """W_m = ∫_1^q y^{-α} ℓ_m(y) dy for the Lagrange basis on nodes y = q^m (m ∈ ms)."""
    ys = [q ** m for m in ms]
    W = []
    for i, yi in enumerate(ys):
        coef, den = [iv.mpf(1)], iv.mpf(1)
        for j, yj in enumerate(ys):
            if j == i:
                continue
            new = [iv.mpf(0)] * (len(coef) + 1)
            for d, c in enumerate(coef):
                new[d + 1] += c
                new[d] -= c * yj
            coef, den = new, den * (yi - yj)
        tot = iv.mpf(0)
        for d, c in enumerate(coef):
            e = d + 1 - a
            tot += c * (q ** e - 1) / e
        W.append(tot / den)
    return W


def power_ceiling(beta_p: float = 0.35, sig0: float = 0.01):
    """Re-certify AQ (4): r_n ≤ B s^{-β'} on (0, σ0] for every n (C1–C3).  α-free."""
    bp, S0 = iv.mpf(beta_p), iv.mpf(sig0)
    sstar = (iv.sqrt(25 + 4 * S0) - 5) / 2
    _, X = node_bounds(iv.mpf(sig0))
    Xb = iv.mpf(X.b)
    M = 400
    lo_e, hi_e = mpf(sstar.a), mpf(sig0)
    need, cells = mpf(0), []
    for i in range(M):
        U = iv.mpf(lo_e + (hi_e - lo_e) * i / M)
        W = iv.mpf(lo_e + (hi_e - lo_e) * (i + 1) / M) if i < M - 1 else iv.mpf(hi_e)
        den = (6 + U) - 2 * W * (5 + W) * Xb
        val = (2 + W) * (5 + W) * Xb / den
        cells.append((W, den, val))
        need = max(need, mpf((val * W ** bp).b))
    B = iv.mpf(need * mpf("1.0001"))
    ok2 = all(den.a > 0 and val.b <= (B * W ** (-bp)).a for W, den, val in cells)
    ss = iv.mpf(sstar.b)
    lhs = (2 + ss) * (5 + ss) ** (1 - bp) + 2 * B * ss ** (1 - bp) * (5 + ss) ** (1 - bp)
    ok1 = lhs.b < 6
    ok3 = (S0 ** bp / 3).b <= B.a
    return B, ok1 and ok2 and ok3, float(lhs.b)


def build_nodes():
    q = iv.mpf(Q_NUM) / Q_DEN
    s0 = iv.mpf(2) ** (-EX0)
    M = max(FLOOR_BU + CEIL_BU)
    vals, k = {}, -M
    while True:
        x = s0 * q ** k
        lo, hi = node_bounds(x)
        vals[k] = (x, lo, hi)
        if k >= M and vals[k - M][0].a > X_TAIL:
            break
        k += 1
    return q, s0, vals, k - M  # cells c = 0..N-1 cover [x_0, x_N], x_N > X_TAIL


def certify_limit(alpha: str, q, s0, vals, N, B):
    a = iv.mpf(alpha)
    bl, ul = FLOOR_BU
    bu, uu = CEIL_BU
    Wf = lagrange_weights(q, list(range(1 - bl, 1 + ul)), a)
    Wc = lagrange_weights(q, list(range(1 - bu, 1 + uu)), a)

    def rv(k):
        _, lo, hi = vals[k]
        return iv.mpf([lo.a, hi.b])

    def wv(k):
        x, lo, hi = vals[k]
        return iv.mpf([(1 - x * hi).a, (1 - x * lo).b])

    J = {}
    xN, x0 = vals[N][0], vals[0][0]
    for name, f in (("r", rv), ("w", wv)):
        lo_sum, up_sum = iv.mpf(0), iv.mpf(0)
        for c in range(N):
            sc = vals[c][0] ** (1 - a)
            fl = sum((Wf[i] * f(c + 1 - bl + i) for i in range(bl + ul)), iv.mpf(0))
            cu = sum((Wc[i] * f(c + 1 - bu + i) for i in range(bu + uu)), iv.mpf(0))
            lo_sum += sc * iv.mpf(fl.a)
            up_sum += sc * iv.mpf(cu.b)
        J[name] = [lo_sum, up_sum]
    # head (0, x0]: f decreasing, r ≤ B s^-0.35, w ≤ 1
    bp = iv.mpf("0.35")
    h1 = x0 ** (1 - a) / (1 - a)
    J["r"][0] += rv(0).a * h1
    J["r"][1] += B * x0 ** (1 - a - bp) / (1 - a - bp)
    J["w"][0] += wv(0).a * h1
    J["w"][1] += h1
    # tail [xN, ∞): r ∈ [1/s - 2/s², 1/s], w ∈ [2/s - 6/s², 2/s]
    t1 = xN ** (-a) / a
    t2 = xN ** (-1 - a) / (1 + a)
    J["r"][0] += t1 - 2 * t2
    J["r"][1] += t1
    J["w"][0] += 2 * t1 - 6 * t2
    J["w"][1] += 2 * t1
    Jr = iv.mpf([J["r"][0].a, J["r"][1].b])
    Jw = iv.mpf([J["w"][0].a, J["w"][1].b])
    pref = iv.sin(iv.pi * a) / iv.pi
    R = 2 * pref * Jr
    S = pref * Jw
    Lim = iv.sqrt(21) / 6 * Jw / Jr
    return {"Jr": Jr, "Jw": Jw, "R": R, "S": S, "L": Lim}


def fmt(x, d=13):
    return f"[{mp.nstr(mpf(x.a), d)}, {mp.nstr(mpf(x.b), d)}]"


def r_rec(n: int, s):
    orb = [s]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    x = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        o = orb[k]
        x = (2 + o) * (5 + o) * x / ((6 + o) - 2 * o * (5 + o) * x)
    return x


def float_F(n: int, alpha: float) -> float:
    """||F(n,α)|| from the exact recursion by log-grid quadrature (float; consistency only)."""
    t = np.linspace(np.log(1e-60), np.log(1e8), 200001)
    sg = np.exp(t)
    h = t[1] - t[0]
    orb = [sg]
    for _ in range(n):
        orb.append(np.minimum(orb[-1] * (orb[-1] + 5), 1e90))
    x = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        o = orb[k]
        x = np.where(o > 1e40, 1 / (3 + o), (2 + o) * (5 + o) * x / ((6 + o) - 2 * o * (5 + o) * x))
    a = alpha
    # tails beyond 1e8: r = 1/s + O(s^-2), 1 - s r = 2/s + O(s^-2) (float cancellation forbids going higher)
    Jr = np.sum(sg ** (1 - a) * x) * h + 1e8 ** (-a) / a
    Jw = np.sum(sg ** (1 - a) * (1 - sg * x)) * h + 2 * 1e8 ** (-a) / a
    return math.sqrt(21) / 6 * Jw / Jr


def mp_F(n: int, alpha: str):
    """||F(n,α)|| by mpmath quadrature of the exact recursion (not interval; consistency only)."""
    mp.dps = 30
    a = mpf(alpha)
    X = mpf(10) ** 30
    pts = [0] + [mpf(10) ** k for k in range(-30, 31, 2)]
    Jr = mp.quad(lambda s: s ** (-a) * r_rec(n, s), pts) + X ** (-a) / a
    Jw = mp.quad(lambda s: s ** (-a) * (1 - s * r_rec(n, s)), pts) + 2 * X ** (-a) / a
    out = mp.sqrt(21) / 6 * Jw / Jr
    mp.dps = 50
    return out


def main() -> int:
    iv.dps = 50
    mp.dps = 50
    print("Theorem AS: certified digits of R_∞, σ_∞ and L = lim ||F(n,α)|| (complete-monotone interpolation sandwich)")
    print("Assumption: Theorem AQ (1)–(4) as stated (decimation, Stieltjes formulas, monotonicity, pointwise ceilings).")
    print("Not thrust.\n")

    # (1)–(2) sign rule for the interpolation sandwich, on actual CM functions
    for n in (2, 3):
        P, E, C = build_gasket(n)
        L = combinatorial_laplacian(len(P), E)
        v = np.zeros(len(P)); v[int(C[0])], v[int(C[1])] = 1.0, -1.0
        check(abs(0.5 * v @ L @ v - 2) < 1e-12 and abs(0.5 * v @ L @ L @ v - 6) < 1e-12,
              f"level {n}: ½vᵀLv = 2, ½vᵀL²v = 6 (tail moments; corners at distance ≥ 4 for n ≥ 2)")
        lam, V = np.linalg.eigh(L)
        c = 0.5 * (V.T @ v) ** 2
        check(abs(c.sum() - 1) < 1e-12 and c.min() > -1e-15,
              f"level {n}: r_n = Σ c_i/(λ_i+s), w_n = Σ c_iλ_i/(λ_i+s), c_i ≥ 0, Σc_i = 1 (complete monotonicity)")
    mp.dps = 60
    worst_sign = True
    for n in (3, 8):
        for (b, u) in ((2, 2), (3, 3), (10, 10), (11, 11)):
            xs = [mpf("0.37") * mpf(17) ** m / mpf(16) ** m for m in range(1 - b, 1 + u)]
            for which in ("r", "w"):
                f = (lambda s: r_rec(n, s)) if which == "r" else (lambda s: 1 - s * r_rec(n, s))
                fx = [f(x) for x in xs]
                for tt in (mpf("0.25"), mpf("0.5"), mpf("0.75")):
                    xm = xs[b - 1] + tt * (xs[b] - xs[b - 1])
                    p = mpf(0)
                    for i, xi in enumerate(xs):
                        li = mpf(1)
                        for j, xj in enumerate(xs):
                            if j != i:
                                li *= (xm - xj) / (xi - xj)
                        p += fx[i] * li
                    d = f(xm) - p
                    if (b % 2 == 0 and d < 0) or (b % 2 == 1 and d > 0):
                        worst_sign = False
    check(worst_sign, "sandwich sign rule sign(f - p) = (-1)^b holds for r_n, w_n (n = 3, 8; b = u ∈ {2,3,10,11})")
    mp.dps = 50

    # power ceiling (AQ (4)), α-free
    B, okB, c1 = power_ceiling()
    check(okB, f"AQ (4) power ceiling re-certified: r_n ≤ B s^-0.35 on (0, 0.01], B = {mp.nstr(mpf(B.b), 9)}, C1 lhs {c1:.6f} < 6")

    t0 = time.time()
    q, s0, vals, N = build_nodes()
    print(f"  nodes: x_k = 2^-{EX0}·(17/16)^k, {len(vals)} nodes, cells 0..{N-1}, x_N > 1e25 "
          f"(node enclosures in {time.time()-t0:.1f}s)")
    check(all(vals[k][1].a <= vals[k][2].b and vals[k][1].a > 0 for k in vals), "node enclosures 0 < lo ≤ hi")

    res = {}
    print("\n  α      R_∞ = lim R_α(n)                    σ_∞ = lim σ_n                       L = lim ||F(n,α)||")
    for al in ("0.10", "0.25", "0.45", "0.49", "0.55"):
        c = certify_limit(al, q, s0, vals, N, B)
        res[al] = c
        print(f"  {al}  {fmt(c['R'])}  {fmt(c['S'])}  {fmt(c['L'])}")
        wL = float(mpf(c["L"].b) - mpf(c["L"].a))
        tol = 1e-10 if al != "0.55" else 1e-7   # at 0.55 the head term B s0^(0.65-α) = 2^-25 dominates
        check(wL < tol, f"α={al}: certified width of L is {wL:.1e} < {tol:.0e}")

    c = res["0.45"]
    La, Lb = float(mpf(c["L"].a)), float(mpf(c["L"].b))
    check(AQ_BRACKET[0] <= La and Lb <= AQ_BRACKET[1], "α=0.45: L lies inside AQ's bracket [0.819255, 0.822959]")
    check(float(mpf(c["R"].b)) < 2.518955 and float(mpf(c["S"].a)) > 1.350986, "α=0.45: R_∞, σ_∞ inside AQ's ceilings/floors")
    check(float(mpf(c["S"].a)) > 2 * 4 ** -0.45, "α=0.45: σ_∞ ≥ 2·4^-α (Theorem AR (3))")
    F60 = mp_F(60, "0.45")
    print(f"  recursion (mpmath quad, 30 digits): ||F(60,0.45)|| = {mp.nstr(F60, 15)}; ||F(60)|| - L = {float(F60 - mpf(c['L'].b)):.2e}")
    check(F60 > mpf(c["L"].b) and F60 - mpf(c["L"].a) < 1e-9,
          "α=0.45: L < ||F(60)|| < L + 1e-9 (decrease to L; AR's gap ~ρ^60 ≈ 1.7e-10·const)")
    check(abs(La - 0.821138627) < 1e-9, "α=0.45: L = 0.821138627... (AQ's recursion value 0.821139 at n = 40)")

    # (5) Conjecture J.2
    c10 = res["0.10"]
    L10 = float(mpf(c10["L"].a))
    print(f"\n  Conjecture J.2: L(0.10) ≥ {L10:.12f}; recorded ||F(6,0.10)|| = {J2_LEVELS[-1]}")
    check(L10 > 1.410, "Conjecture J.2 PROVED: ||F(n,0.10)|| > L(0.10) > 1.410 at every level (AQ (3) monotone decrease)")
    check(all(J2_LEVELS[i] > J2_LEVELS[i + 1] for i in range(5)) and J2_LEVELS[-1] > float(mpf(c10["L"].b)),
          "recorded levels 1-6 at α=0.10 decrease and stay above the certified limit")
    for al in ("0.10", "0.25", "0.49", "0.55"):
        Fn = float_F(30, float(al))
        ok = Fn >= float(mpf(res[al]["L"].a)) - 1e-10
        check(ok, f"α={al}: float ||F(30)|| = {Fn:.10f} ≥ certified L (decrease to L)")

    # Kept Failure AS.1
    rho = 5 ** 0.45 / 3
    ext = NF10 - (NF9 - NF10) * rho / (1 - rho)
    F8, F9, F10 = float_F(8, 0.45), float_F(9, 0.45), float_F(10, 0.45)
    ait = F10 - (F10 - F9) ** 2 / ((F10 - F9) - (F9 - F8))
    eff = (F10 - La) / (F9 - La)
    print(f"\n  Kept Failure AS.1: ρ = 5^0.45/3 = {rho:.6f}; ρ-extrapolation from levels 9,10: {ext:.6f}; "
          f"Aitken(8,9,10): {ait:.6f}; certified L = {La:.9f}")
    print(f"  effective ratio (||F(10)||-L)/(||F(9)||-L) = {eff:.5f} vs ρ = {rho:.5f}")
    check(abs(F10 - NF10) < 5e-6 and abs(F9 - NF9) < 5e-6, "float recursion reproduces ||F(9)||, ||F(10)||")
    check(abs(ext - La) > 1e-4 and abs(ait - La) > 1e-4,
          f"AS.1: ρ-extrapolation misses L by {abs(ext-La):.2e} and Aitken by {abs(ait-La):.2e} (> 1e-4): claim falsified")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        for f in failures:
            print("   ", f)
        return 1
    print("All Theorem AS / Kept Failure AS.1 checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
