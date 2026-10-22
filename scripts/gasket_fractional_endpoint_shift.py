#!/usr/bin/env python3
"""Certified intercept C_γ and the level shift of the endpoint audit norm (Theorem AY).

Witness for Theorem AY and Kept Failure AY.1 in GASKET.md.

Notation as in Theorems AQ-AU: r_n(s) = ½ vᵀ(L_n + s)⁻¹ v, r_{n+1}(s) = Φ_s(r_n(p(s))),
p(s) = s(s+5), q = p⁻¹, r_∞ = lim r_n, δ_{j+1} = r_{j+1} - r_j, γ = log3/log5,
κ = 2 sin πγ/π, J_n = ∫_0^∞ s^{-γ} r_n ds (so R_γ(n) = κ J_n, AQ (2)),
c = ln5 <Θ> (so c_R(γ) = κ c, Theorem AU (1)).  At α = γ, ρ = 5^γ/3 = 1.
Bracket weights of Theorem AT (2), all independent of the level:
    B_∞ = Π_{i≥0} 5/√(25+4q^i t) ≤ 1,  A_∞ = Π_{i≥1}(1+q^i t/5) ≥ 1,  C̄ = Π_{i≥1} (3/5) D̄(q^i t) ≥ 1.

Theorem AY.
  (1) Summable defect.  For every j ≥ 0, both Δ_j = J_{j+1} - J_j and c lie in
      [∫ t^{-γ} δ_{j+1} B_∞, ∫ t^{-γ} δ_{j+1} A_∞^γ C̄] (AT (2) at ρ = 1; Δ_j is the k = 0 term), so
        -∫ t^{-γ} δ_{j+1} (1 - B_∞) ≤ c - Δ_j ≤ ∫ t^{-γ} δ_{j+1} (A_∞^γ C̄ - 1).
      Summing over j ≥ n (Tonelli, Σ_{j≥n} δ_{j+1} = r_∞ - r_n): C' = lim_m (J_m - m c) exists and
        J_n - n c - C' ∈ [-E_n^-, E_n^+],  E_n^± = ∫_0^∞ t^{-γ}(r_∞ - r_n) W_±,
        W_+ = A_∞^γ C̄ - 1,  W_- = 1 - B_∞.
      E_m^± ≤ E_n^± for m ≥ n because r_m increases in m.  So with C_γ = κ C',
        R_γ(m) - c_R(γ) m - C_γ ∈ κ[-E_n^-, E_n^+]  for every m ≥ n.
  (2) Certified (outward-rounded mpmath.iv, 80 digits; complete-monotone sandwich of Theorem AS (2)
      for r_n and r_∞ on a geometric grid of ratio 65/64 from 5^-60 2^-110 to 1e25, floor b = u = 14,
      ceiling b = u = 15; W_± bounded on blocks of 16 cells with Theorem AT's weight routine; head by
      AQ's power ceiling, tail beyond 1 by AT's crude tail; level n = 60, cross-checked at n = 40):
        E_60^+ ≤ 8.7e-28, E_60^- ≤ 5.9e-34;  E_40^+ ≤ 1.1e-18,
        C_γ    ∈ [0.70307046242185742992904816, 0.70307046242185742992904817],
        ν_γ = C_γ/c_R(γ) = C'/c ∈ [1.35946811514371162883271469, 1.35946811514371162883271470],
        σ_∞(γ) = (sin πγ/π) ∫ s^{-γ}(1 - s r_∞) ∈ [1.15191861167405417640547243, ...244],
        K_γ = √21 σ_∞(γ)/(3 c_R(γ)) ∈ [3.40236076383265507999048519, ...520].
      Theorem AU's float observation C_γ ≈ 0.70307046242184 was low by 1.7e-14.
  (3) Level-shift law (DERIVED).  ||F(n,γ)|| = (√21/3) σ_n/R_γ(n), σ_n - σ_∞ = O(5^{-n})
      (AR/AT at μ = 1/5) and (1) (rate O(5^{-0.65 n}) by AT (4a)) give
        (n + ν_γ) ||F(n,γ)|| = K_γ + O(5^{-0.65 n}).

Kept Failure AY.1.  Claim: at the endpoint, as below it, the renormalized audit norm converges
exponentially fast without a shift, n ||F(n,γ)|| = K_γ + O(θ^n) for some θ < 1.  False:
(3) gives n ||F(n,γ)|| - K_γ = -ν_γ K_γ/(n+ν_γ) + O(5^{-0.65 n}), and ν_γ K_γ ≥ 4.6254 is
certified, so n(n||F|| - K_γ) → -ν_γ K_γ ≠ 0: the error is Θ(1/n), not exponential.

Not thrust. Does not touch Stage 1 or Stage 2.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_limit_digits import lagrange_weights, power_ceiling  # noqa: E402
import gasket_fractional_rate_constants as AT  # noqa: E402

failures: list[str] = []

N_LEVEL = 60                 # main level n
N_CROSS = 40                 # cross-check level
Q_NUM, Q_DEN = 65, 64        # geometric node ratio
ORD = 14                     # floor b = u = 14 (even), ceiling b = u = 15 (odd)
EX0 = 110                    # x_0 = 5^-60 2^-110
X_TAIL = mpf(10) ** 25       # tail cut
CAP = mpf(10) ** 90          # orbit cap in the AQ (4) node enclosure of r_∞
BLOCK = 16                   # cells per block for the bracket weights W_±
AU_THETA = ("0.6009542511687736078237051500999942", "0.6009542511687736078237051500999943")  # Theorem AU (3)
AT_CRG = ("0.517165834630505243", "0.517165834630505245")       # Theorem AT (3)
AT_CSG = ("0.481399003768729118", "0.481399003768729121")       # Theorem AT (3), c_S(γ)
AU_OBS_C = "0.70307046242184"                                   # Theorem AU witness observation (float)


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def r_exact(n: int, x):
    orb = [x]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    r = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        o = orb[k]
        r = (2 + o) * (5 + o) * r / ((6 + o) - 2 * o * (5 + o) * r)
    return r


def rinf_bounds(x):
    """AQ (4): T^k r_0(x) ≤ r_∞(x) ≤ T^k(1/s)(x), iterated until the orbit passes CAP."""
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
    return iv.mpf([lo.a, hi.b])


def build(levels):
    q = iv.mpf(Q_NUM) / Q_DEN
    x0 = iv.mpf(5) ** (-N_LEVEL) * iv.mpf(2) ** (-EX0)
    M = ORD + 1
    vals, k = {}, -M
    while True:
        x = x0 * q ** k
        vals[k] = {"x": x, "inf": rinf_bounds(x)}
        for n in levels:
            vals[k][n] = r_exact(n, x)
        if k >= M and vals[k - M]["x"].a > X_TAIL:
            break
        k += 1
    N = k - M                                   # cells c = 0..N-1 cover [x_0, x_N], x_N > X_TAIL
    N1 = max(c for c in range(N) if vals[c + 1]["x"].b <= 1)  + 1   # cells c < N1 lie in (x_0, 1]
    return q, x0, vals, N, N1


def sandwich_cells(q, vals, gam, key_lo, key_hi, cells):
    """Per cell c: [floor ∫ s^{-γ} f_lo, ceil ∫ s^{-γ} f_hi] (Theorem AS (2)); f given by node keys."""
    Wf = lagrange_weights(q, list(range(1 - ORD, 1 + ORD)), gam)
    Wc = lagrange_weights(q, list(range(-ORD, 2 + ORD)), gam)
    out = {}
    for c in cells:
        sc = vals[c]["x"] ** (1 - gam)
        fl = sum((Wf[i] * key_lo(c + 1 - ORD + i) for i in range(2 * ORD)), iv.mpf(0))
        cu = sum((Wc[i] * key_hi(c - ORD + i) for i in range(2 * ORD + 2)), iv.mpf(0))
        out[c] = (mpf((sc * fl).a), mpf((sc * cu).b))
    return out


def J_level(n, q, x0, vals, N, gam, cells):
    """Enclosure of J_n = ∫_0^∞ s^{-γ} r_n: sandwich on [x_0, x_N], head by monotonicity, tail by moments."""
    lo = sum((iv.mpf(cells[c][0]) for c in range(N)), iv.mpf(0))
    hi = sum((iv.mpf(cells[c][1]) for c in range(N)), iv.mpf(0))
    h = x0 ** (1 - gam) / (1 - gam)            # head: r_n(x_0) ≤ r_n ≤ r_n(0) = (5/3)^n/3
    lo += iv.mpf(vals[0][n].a) * iv.mpf(h.a)
    hi += (iv.mpf(5) / 3) ** n / 3 * h
    X = vals[N]["x"]                            # tail: 1/s - 2/s² ≤ r_n ≤ 1/s (½vᵀLv = 2, n ≥ 1)
    lo += X ** (-gam) / gam - 2 * X ** (-1 - gam) / (1 + gam)
    hi += X ** (-gam) / gam
    return iv.mpf([lo.a, hi.b])


def weight_blocks(q, vals, N1, gam):
    """Upper bounds of W_+ = A_∞^γ C̄ - 1 and W_- = 1 - B_∞ on blocks of BLOCK cells in (x_0, 1]."""
    out = {}
    for b0 in range(0, N1, BLOCK):
        b1 = min(b0 + BLOCK, N1)
        lam_b, lam_psi, cbar = AT.w_data(vals[b0]["x"], vals[b1]["x"])
        wp = lam_psi ** (-gam) * cbar - 1      # lam_psi ≤ A_∞^{-1}(x_hi) ≤ A_∞^{-1} on the block
        wm = 1 - lam_b
        for c in range(b0, b1):
            out[c] = (mpf(wp.b), mpf(wm.b))
    return out


def defect(n, q, x0, vals, N1, gam, cinf, cn, wts, B, P065, K065):
    """E_n^± = ∫ t^{-γ}(r_∞ - r_n) W_± over (0,∞): head (0,x_0], cells (x_0, x_N1], tail [x_N1, ∞)."""
    Ep, Em = iv.mpf(0), iv.mpf(0)
    for c in range(N1):
        d = max(cinf[c][1] - cn[c][0], mpf(0))   # ceil ∫ s^-γ r_∞ - floor ∫ s^-γ r_n ≥ ∫ s^-γ (r_∞ - r_n) ≥ 0
        Ep += iv.mpf(d) * iv.mpf(wts[c][0])
        Em += iv.mpf(d) * iv.mpf(wts[c][1])
    # head (0, x_0]: r_∞ - r_n ≤ r_∞ ≤ B t^-0.35 (AQ (4)); log A_∞ ≤ t/20 (AU); log C̄ ≤ 1.82 K t^0.65 with
    # K = Σ_{i≥1} 5^{-0.65 i} (AT (3) small-y lemma, q^i t ≤ t 5^-i); e^u - 1 ≤ u e^u; 1 - B_∞ ≤ t/10.
    u0 = gam * x0 / 20 + iv.mpf("1.82") * (K065 - 1) * x0 ** P065
    Kp = (gam * x0 ** (1 - P065) / 20 + iv.mpf("1.82") * (K065 - 1)) * iv.exp(u0)    # W_+ ≤ Kp t^0.65 on (0, x_0]
    e1 = 1 - gam - iv.mpf("0.35") + P065                                           # exponent of t^{-γ} t^{-0.35} t^{0.65}, +1
    Ep += B * Kp * x0 ** e1 / e1
    e2 = 1 - gam - iv.mpf("0.35") + 1
    Em += B / 10 * x0 ** e2 / e2
    # [x_N1, 1]: one partial cell, same bound with the block weight of the last block
    # (covered: x_N1 ≤ 1 < x_{N1+1}; the cell [x_N1, x_{N1+1}] ⊃ [x_N1, 1] uses the sandwich of cell N1)
    d = max(cinf[N1][1] - cn[N1][0], mpf(0))
    lam_b, lam_psi, cbar = AT.w_data(vals[N1]["x"], vals[N1 + 1]["x"])
    Ep += iv.mpf(d) * (lam_psi ** (-gam) * cbar - 1)
    Em += iv.mpf(d) * (1 - lam_b)
    return Ep, Em


def tail_log10(n):
    """log10 of Σ_{j≥n} ∫_{x≥1} t^{-γ} δ_{j+1} A_∞^γ C̄ ≤ Σ_j 6e10·1.2^j 6^{-2^{j-1}}/(2^{j-1}-5) ≤ 2·(j = n term)."""
    return (math.log10(2 * 6e10) + n * math.log10(1.2) - 2 ** (n - 1) * math.log10(6) - math.log10(2 ** (n - 1) - 5))


def r_mp(n, s):
    orb = [s]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    x = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        o = orb[k]
        x = (2 + o) * (5 + o) * x / ((6 + o) - 2 * o * (5 + o) * x)
    return x


def float_level(m, g):
    """(J_m, σ_m) by mpmath quadrature (consistency only, not certified)."""
    pts = [mpf(0)] + [mpf(5) ** (-m) * mpf(10) ** k for k in range(-14, 2)] + [mpf(10) ** k for k in range(-6, 31, 2)]
    pts = sorted(set(pts))
    Jm = mp.quad(lambda s: s ** (-g) * r_mp(m, s), pts) + mpf(10) ** (-30 * g) / g
    Jw = mp.quad(lambda s: s ** (-g) * (1 - s * r_mp(m, s)), pts) + 2 * mpf(10) ** (-30 * g) / g
    return Jm, mp.sin(mp.pi * g) / mp.pi * Jw


def fmt(x, d=26):
    sc = mpf(10) ** d
    return f"[{mp.nstr(mp.floor(mpf(x.a) * sc) / sc, d + 2)}, {mp.nstr(mp.ceil(mpf(x.b) * sc) / sc, d + 2)}]"


def main() -> int:
    iv.dps = 80
    mp.dps = 80
    print("Theorem AY: certified intercept C_γ and the level shift of the endpoint audit norm")
    print("Assumptions: Theorems AQ (1)-(4), AR (1)-(2), AS (2), AT (2)-(3) and AU (1),(3), all as stated. Not thrust.\n")
    P065 = iv.mpf(13) / 20
    K065 = 1 / (1 - iv.mpf(5) ** (-P065))
    gam = iv.log(3) / iv.log(5)
    kap = 2 * iv.sin(iv.pi * gam) / iv.pi

    iv.dps = 60
    B, okB, c1 = power_ceiling()
    iv.dps = 80
    BH = iv.mpf(B.b)
    AT.P065, AT.K065, AT.BH = P065, K065, BH
    check(okB, f"AQ (4) power ceiling re-certified: r_m ≤ B s^-0.35 on (0,0.01], B = {mp.nstr(mpf(B.b), 9)} (C1 lhs {c1:.6f} < 6)")
    y = iv.mpf("0.001")
    eps = 2 * BH * (y * (5 + y)) ** P065 / 6
    lem = (y * (3 * y + 16) / (3 * (6 + y)) + 2 * eps / (1 - eps)) / y ** P065
    check(lem.b < 1.82, f"AT (3) small-y lemma: log((3/5)D̄(y)) ≤ 1.82 y^0.65 for y ≤ 1e-3 (ratio {mp.nstr(mpf(lem.b), 5)})")
    xs = [iv.mpf(t) for t in ("1e-30", "0.001", "1")]
    check(all((-iv.log(5 / iv.sqrt(25 + 4 * x))).b <= (2 * x / 25).a for x in xs),
          "½ log(1+4y/25) ≤ 2y/25 (so 1 - B_∞(t) ≤ Σ_i 2 q^i(t)/25 ≤ t/10), spot check")

    c_th = iv.mpf(list(AU_THETA))
    cc = iv.log(5) * c_th                        # c = ln5 <Θ>
    check(mpf((kap * cc).a) >= mpf(AT_CRG[0]) and mpf((kap * cc).b) <= mpf(AT_CRG[1]),
          "κ ln5 <Θ> (Theorem AU) lies inside Theorem AT's c_R(γ): the two cited inputs agree")

    t0 = time.time()
    q, x0, vals, N, N1 = build((N_LEVEL, N_CROSS))
    print(f"  grid: {len(vals)} nodes (ratio {Q_NUM}/{Q_DEN}), {N} cells up to {mp.nstr(mpf(vals[N]['x'].a), 3)}, "
          f"{N1} cells in (x_0, 1] ({time.time()-t0:.0f}s)")
    check(all(vals[k]["inf"].a > 0 and mpf(vals[k]["inf"].b) - mpf(vals[k]["inf"].a) < mpf(10) ** -60 * mpf(vals[k]["inf"].a)
              for k in vals), "node enclosures of r_∞ positive with relative width < 1e-60")
    check(all(vals[k][N_CROSS].a <= vals[k][N_LEVEL].b and vals[k][N_LEVEL].a <= vals[k]["inf"].b for k in vals),
          f"node order r_{N_CROSS} ≤ r_{N_LEVEL} ≤ r_∞ (AQ (3))")

    allc = range(N + 1)
    cinf = sandwich_cells(q, vals, gam, lambda k: iv.mpf(vals[k]["inf"].a), lambda k: iv.mpf(vals[k]["inf"].b), allc)
    wts = weight_blocks(q, vals, N1, gam)
    print(f"  sandwich and weights ({time.time()-t0:.0f}s); W_+ ≤ {mp.nstr(max(w[0] for w in wts.values()), 4)}, "
          f"W_- ≤ {mp.nstr(max(w[1] for w in wts.values()), 4)} on (0, 1]")
    check(all(w[0] >= 0 and w[1] >= 0 for w in wts.values()), "bracket weights W_± ≥ 0 (B_∞ ≤ 1 ≤ A_∞^γ C̄)")

    res = {}
    for n in (N_LEVEL, N_CROSS):
        cn = sandwich_cells(q, vals, gam, lambda k: iv.mpf(vals[k][n].a), lambda k: iv.mpf(vals[k][n].b), allc)
        Jn = J_level(n, q, x0, vals, N, gam, cn)
        Ep, Em = defect(n, q, x0, vals, N1, gam, cinf, cn, wts, BH, P065, K065)
        tl = tail_log10(n)
        tail = iv.mpf(10) ** iv.mpf(tl)
        Ep, Em = Ep + tail, Em + tail
        base = Jn - n * cc
        Cp = iv.mpf([(base - Ep).a, (base + Em).b])
        res[n] = (Jn, Ep, Em, Cp)
        print(f"\n  n = {n}: J_n ∈ {fmt(Jn, 30)} (width {float(mpf(Jn.b) - mpf(Jn.a)):.1e})")
        print(f"          J_n - n c ∈ {fmt(base, 30)}")
        print(f"          E_n^+ ≤ {mp.nstr(mpf(Ep.b), 4)},  E_n^- ≤ {mp.nstr(mpf(Em.b), 4)},  tail [1,∞) ≤ 10^{tl:.4g}")
        print(f"          C' = lim (J_m - m c) ∈ {fmt(Cp, 30)}")
    Cp60, Cp40 = res[N_LEVEL][3], res[N_CROSS][3]
    check(max(mpf(Cp60.a), mpf(Cp40.a)) <= min(mpf(Cp60.b), mpf(Cp40.b)),
          f"C' enclosures from levels {N_LEVEL} and {N_CROSS} overlap (independent sandwiches of r_n)")
    E60 = max(mpf(res[N_LEVEL][1].b), mpf(res[N_LEVEL][2].b))
    E40 = max(mpf(res[N_CROSS][1].b), mpf(res[N_CROSS][2].b))
    check(E40 < 5 ** (-0.6 * 40) and E60 < 5 ** (-0.6 * 60),
          f"defect bounds E_40 = {float(E40):.2e}, E_60 = {float(E60):.2e} below 5^(-0.6 n); E_m ≤ E_n for m ≥ n (r_m ↑)")
    check(mpf(res[N_LEVEL][2].b) < mpf(res[N_LEVEL][1].b), "W_- side is the smaller one (1 - B_∞ = O(t), W_+ = O(t^0.65))")
    Cp = iv.mpf([max(mpf(Cp60.a), mpf(Cp40.a)), min(mpf(Cp60.b), mpf(Cp40.b))])
    Cg = kap * Cp
    nu = Cp / cc
    wC = mpf(Cg.b) - mpf(Cg.a)
    print(f"\n  C_γ = κ C'               ∈ {fmt(Cg)}   (width {float(wC):.1e})")
    print(f"  ν_γ = C_γ/c_R(γ) = C'/c  ∈ {fmt(nu)}")
    check(wC < 1e-25, "certified width of C_γ < 1e-25")
    obs = mpf(AU_OBS_C)
    dev = mpf(Cg.a) - obs
    check(1e-15 < dev < 1e-13,
          f"Theorem AU's float observation C_γ ≈ {AU_OBS_C} agrees to 13 digits but is low by {float(dev):.2e} "
          f"(its 14th digit is wrong: certified {mp.nstr(mpf(Cg.a), 15)}...)")

    # σ_∞(γ) = (sin πγ/π) ∫ s^-γ (1 - s r_∞): w_∞ = 1 - s r_∞ completely monotone (AS (2))
    wlo = lambda k: iv.mpf((1 - vals[k]["x"] * iv.mpf(vals[k]["inf"].b)).a)
    whi = lambda k: iv.mpf((1 - vals[k]["x"] * iv.mpf(vals[k]["inf"].a)).b)
    cw = sandwich_cells(q, vals, gam, wlo, whi, range(N))
    lo = sum((iv.mpf(cw[c][0]) for c in range(N)), iv.mpf(0))
    hi = sum((iv.mpf(cw[c][1]) for c in range(N)), iv.mpf(0))
    h = x0 ** (1 - gam) / (1 - gam)            # head: w_∞(x_0) ≤ w_∞ ≤ 1
    lo += wlo(0) * iv.mpf(h.a)
    hi += h
    X = vals[N]["x"]                            # tail: 2/s - 6/s² ≤ w_∞ ≤ 2/s (AS (3))
    lo += 2 * X ** (-gam) / gam - 6 * X ** (-1 - gam) / (1 + gam)
    hi += 2 * X ** (-gam) / gam
    sig = iv.sin(iv.pi * gam) / iv.pi * iv.mpf([lo.a, hi.b])
    cR = kap * cc
    Kg = iv.sqrt(21) / 3 * sig / cR
    print(f"  σ_∞(γ)                   ∈ {fmt(sig)}")
    print(f"  K_γ = √21 σ_∞/(3 c_R)    ∈ {fmt(Kg)}")
    check(mpf(sig.b) - mpf(sig.a) < 1e-25 and mpf(sig.a) >= 2 * 4 ** (-float(gam.a)),
          "σ_∞(γ) certified (width < 1e-25) and ≥ 2·4^-γ (Theorem AR (3))")

    # consistency (float, 30-digit quadrature): J_m - m c → C', σ_m → σ_∞, level-shift law
    mp.dps = 30
    g = mp.log(3) / mp.log(5)
    cf, Cpf, nuf, Kf, sigf = (mpf(z.a) for z in (cc, Cp, nu, Kg, sig))
    cRf = mpf(cR.a)
    print("\n   m    J_m - m c - C'      σ_m - σ_∞     (m+ν)||F(m,γ)|| - K_γ    m||F(m,γ)|| - K_γ    m(m||F|| - K_γ)")
    rows = []
    for m in (2, 5, 10, 20, 30):
        Jm, sm = float_level(m, g)
        F = mp.sqrt(21) / 3 * sm / (2 * mp.sin(mp.pi * g) / mp.pi * Jm)
        rows.append((m, Jm - m * cf - Cpf, sm - sigf, (m + nuf) * F - Kf, m * F - Kf))
        print(f"  {m:2d}   {float(rows[-1][1]):+.3e}      {float(rows[-1][2]):+.3e}       {float(rows[-1][3]):+.3e}"
              f"              {float(rows[-1][4]):+.6f}          {float(m * rows[-1][4]):+.6f}")
    check(all(0 < r[1] < 5 * 5 ** (-0.65 * r[0]) for r in rows[1:]),
          "0 < J_m - m c - C' < 5·5^(-0.65 m) at m = 5..30 (Theorem AY (1); AT (4a) rate)")
    check(all(abs(r[3]) < 5 * 5 ** (-0.65 * r[0]) for r in rows[1:]),
          "|(m + ν_γ)||F(m,γ)|| - K_γ| < 5·5^(-0.65 m) at m = 5..30 (Theorem AY (3))")
    # Kept Failure AY.1: no-shift renormalization
    lim = -nuf * Kf
    check(mpf(nu.a) > 1 and mpf(Kg.a) > 0, f"AY.1: ν_γ K_γ ≥ {mp.nstr(mpf(nu.a) * mpf(Kg.a), 8)} > 0 certified")
    check(all(r[4] < 0 for r in rows) and abs(30 * rows[-1][4] - lim * 30 / (30 + nuf)) < 1e-10,
          f"AY.1: m||F(m,γ)|| - K_γ < 0 at every m and m(m||F|| - K_γ) → -ν_γ K_γ = {mp.nstr(lim, 10)} "
          f"(at m = 30: {mp.nstr(30 * rows[-1][4], 10)}): the error is Θ(1/m), not exponential; no-shift claim falsified")
    ratios = [rows[i + 1][1] / rows[i][1] for i in range(1, 3)]
    print(f"  observation: (J_m - mc - C') ratios over 5 levels: {[mp.nstr(r_, 6) for r_ in ratios]} "
          f"(3^-5 = {mp.nstr(mpf(3) ** -5, 6)}); not proved")
    mp.dps = 80

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        for f in failures:
            print("   ", f)
        return 1
    print("All Theorem AY / Kept Failure AY.1 checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
