#!/usr/bin/env python3
"""Supercritical decay law of the audit norm and its Mellin intercept (Theorem AZ).

Witness for Theorem AZ and Kept Failure AZ.1 in GASKET.md.

Notation as in Theorems AQ-AY: r_n(s) = ½ vᵀ(L_n + s)⁻¹ v, r_{n+1}(s) = Φ_s(r_n(p(s))),
p(s) = s(s+5), q = p⁻¹, r_∞ = lim r_n, δ_{j+1} = r_{j+1} - r_j, γ = log3/log5,
α ∈ (γ, 1) (supercritical: ρ = 5^α/3 > 1), κ = 2 sin πα/π, J_n = ∫ s^{-α} r_n
(R_α(n) = κ J_n, AQ (2)), Δ_j = J_{j+1} - J_j, c = c_R(α)/κ = lim ρ^{-j} Δ_j (AR (2)).
M(z) = ∫ s^{z-1} r_∞ ds (1-γ < Re z < 1), continued meromorphically (AU (2)).
Bracket weights of Theorem AT (2) (level-free): W_+ = A_∞^α C̄ - 1 ≥ 0, W_- = 1 - B_∞ ≥ 0.

Theorem AZ.
  (1) Geometric expansion with a summable defect.  AT (2) at level j gives
        -∫ t^{-α} δ_{j+1} W_- ≤ c ρ^j - Δ_j ≤ ∫ t^{-α} δ_{j+1} W_+,
      so (Tonelli) C' = lim_m (J_m - c ρ^m/(ρ-1)) exists and for every n
        J_n - c ρ^n/(ρ-1) - C' ∈ [-E_n^-, E_n^+],  E_n^± = ∫ t^{-α}(r_∞ - r_n) W_±,
      i.e. R_α(m) = A_α ρ^m + C_α + O(E_m) with A_α = c_R/(ρ-1), C_α = κ C'.
  (2) Mellin intercept.  For 0 < Re z < 1, ρ_z = 5^{1-z}/3:
        M(z) = I_n(z) + c(z) ρ_z^n/(1-ρ_z) + Σ_{j≥n} e_j(z),  I_n(z) = ∫ s^{z-1} r_n,
      c(z) = lim ρ_z^{-j} ∫ s^{z-1} δ_{j+1} (holomorphic), e_j = ∫ s^{z-1}δ_{j+1} - c ρ_z^j.
      At z = 1-α this is C_α = κ M(1-α): the intercept is the analytic continuation of the
      limit resistance R_∞(α) = κ M(1-α) (α < γ) across the simple pole at γ.
      Consequences: (α-γ) C_α → -c_R(γ)/ln5 as α↓γ, and C_α + c_R(α)/(ρ-1) → C_γ (Theorem AY)
      as α → γ from either side.
  (3) Certified at α = 0.8 (and 0.75): grid 65/64 from 5^-60 2^-110 to 1e25, CM sandwich
      b = u = 14/15 (Theorem AS (2)), levels 60 (main) and 40 (cross-check):
        c_R(0.8), A_0.8, C_0.8 < 0, σ_∞(0.8), ν_0.8 = C/A < 0, K_0.8 = √21 σ_∞/(3 A).
  (4) Decay law (DERIVED).  ||F(n,α)|| = (√21/3) σ_n/R_α(n), σ_n - σ_∞ = O(μ^n), μ = ρ/5:
        (ρ^n + ν_α) ||F(n,α)|| = K_α + O(5^{-0.65 n}),  so ||F(n,α)|| ~ K_α (3·5^{-α})^n,
      with no log-periodic modulation in n (ρ^n is exact at integer n).

Kept Failure AZ.1.  Claim: the endpoint level shift of Theorem AY carries over to α > γ with its
sign, i.e. R_α(n) ≥ A_α ρ^n (C_α ≥ 0) and ρ^n ||F(n,α)|| increases to K_α from below, with error
O(5^{-0.65 n}) or at least below the leading shift.  False: C_0.8 < 0 and ν_0.8 < 0 are certified,
so ρ^n(ρ^n ||F|| - K) → -ν K > 0: ρ^n ||F(n,0.8)|| decreases to K_0.8 from above, the error of the
unshifted renormalization is Θ((3·5^{-0.8})^n), and the sign flip is forced by the simple Mellin
pole at α = γ with positive residue <Θ>.

Not thrust. Does not touch Stage 1 or Stage 2.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_limit_digits import power_ceiling  # noqa: E402
import gasket_fractional_rate_constants as AT  # noqa: E402
import gasket_fractional_endpoint_shift as AY  # noqa: E402

failures: list[str] = []

N_LEVEL, N_CROSS = 60, 40
AT_CR75 = ("0.507985128123073180", "0.507985128123073181")      # Theorem AT (3), c_R(0.75)
AY_CG = ("0.70307046242185742992904816", "0.70307046242185742992904817")   # Theorem AY (2)
AU_THETA = ("0.6009542511687736078237051500999942", "0.6009542511687736078237051500999943")  # Theorem AU (3)


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def c_bracket(n, x0, vals, N1, a, cn, cn1, wts, B, P065, K065):
    """c ρ^n = lim ρ^{n-m} Δ_m ∈ [∫ t^-α δ_{n+1} B_∞, ∫ t^-α δ_{n+1} A_∞^α C̄] (Theorem AT (2))."""
    lo, hi = iv.mpf(0), iv.mpf(0)
    for c in range(N1):
        dlo = max(cn1[c][0] - cn[c][1], mpf(0))   # floor ∫ r_{n+1} - ceil ∫ r_n ≤ ∫ s^-α δ_{n+1} on the cell
        dhi = max(cn1[c][1] - cn[c][0], mpf(0))
        lo += iv.mpf(dlo) * (1 - iv.mpf(wts[c][1]))   # B_∞ ≥ 1 - W_-^sup on the block
        hi += iv.mpf(dhi) * (1 + iv.mpf(wts[c][0]))   # A_∞^α C̄ ≤ 1 + W_+^sup
    # partial cell [x_N1, x_{N1+1}] ⊃ [x_N1, 1] (upper only)
    _, lam_psi, cbar = AT.w_data(vals[N1]["x"], vals[N1 + 1]["x"])
    hi += iv.mpf(max(cn1[N1][1] - cn[N1][0], mpf(0))) * lam_psi ** (-a) * cbar
    # head (0, x_0]: r_{n+1}(x_0) - r_n(0) ≤ δ_{n+1} ≤ r_{n+1}(0) - r_n(x_0); weights as in AY's head
    h = x0 ** (1 - a) / (1 - a)
    r0n = (iv.mpf(5) / 3) ** n / 3
    r0n1 = (iv.mpf(5) / 3) ** (n + 1) / 3
    dl = iv.mpf(AY.r_exact(n + 1, x0).a) - r0n
    if dl.a > 0:
        lo += dl * (1 - x0 / 10) * iv.mpf(h.a)
    u0 = a * x0 / 20 + iv.mpf("1.82") * (K065 - 1) * x0 ** P065
    hi += (r0n1 - iv.mpf(AY.r_exact(n, x0).a)) * iv.exp(u0) * h
    # tail [1, ∞) (upper only): AT's crude tail, same bound as AY.tail_log10 (t^-α ≤ 1, A_∞^α C̄ ≤ 3e10 t^5)
    hi += iv.mpf(10) ** iv.mpf(AY.tail_log10(n))
    return iv.mpf([lo.a, hi.b])


def sigma_inf(q, x0, vals, N, a):
    """σ_∞(α) = (sin πα/π) ∫ s^{-α}(1 - s r_∞) (AS (1), β - 1 = -α); w_∞ completely monotone."""
    wlo = lambda k: iv.mpf((1 - vals[k]["x"] * iv.mpf(vals[k]["inf"].b)).a)
    whi = lambda k: iv.mpf((1 - vals[k]["x"] * iv.mpf(vals[k]["inf"].a)).b)
    cw = AY.sandwich_cells(q, vals, a, wlo, whi, range(N))
    lo = sum((iv.mpf(cw[c][0]) for c in range(N)), iv.mpf(0))
    hi = sum((iv.mpf(cw[c][1]) for c in range(N)), iv.mpf(0))
    h = x0 ** (1 - a) / (1 - a)
    lo += wlo(0) * iv.mpf(h.a)
    hi += h
    X = vals[N]["x"]                            # tail: 2/s - 6/s² ≤ w_∞ ≤ 2/s (AS (3))
    lo += 2 * X ** (-a) / a - 6 * X ** (-1 - a) / (1 + a)
    hi += 2 * X ** (-a) / a
    return iv.sin(iv.pi * a) / iv.pi * iv.mpf([lo.a, hi.b])


def certify_alpha(name, a, levels, q, x0, vals, N, N1, B, P065, K065):
    rho = iv.mpf(5) ** a / 3
    kap = 2 * iv.sin(iv.pi * a) / iv.pi
    allc = range(N + 1)
    cinf = AY.sandwich_cells(q, vals, a, lambda k: iv.mpf(vals[k]["inf"].a), lambda k: iv.mpf(vals[k]["inf"].b), allc)
    wts = AY.weight_blocks(q, vals, N1, a)
    check(all(w[0] >= 0 and w[1] >= 0 for w in wts.values()), f"α={name}: bracket weights W_± ≥ 0 (B_∞ ≤ 1 ≤ A_∞^α C̄)")
    out = {}
    for n in levels:
        cn = AY.sandwich_cells(q, vals, a, lambda k: iv.mpf(vals[k][n].a), lambda k: iv.mpf(vals[k][n].b), allc)
        cn1 = AY.sandwich_cells(q, vals, a, lambda k: iv.mpf(vals[k][n + 1].a), lambda k: iv.mpf(vals[k][n + 1].b), allc)
        Jn = AY.J_level(n, q, x0, vals, N, a, cn)
        Ep, Em = AY.defect(n, q, x0, vals, N1, a, cinf, cn, wts, B, P065, K065)
        tail = iv.mpf(10) ** iv.mpf(AY.tail_log10(n))
        Ep, Em = Ep + tail, Em + tail
        crn = c_bracket(n, x0, vals, N1, a, cn, cn1, wts, B, P065, K065)   # c ρ^n
        c = crn / rho ** n
        base = Jn - crn / (rho - 1)
        Cp = iv.mpf([(base - Ep).a, (base + Em).b])
        out[n] = dict(J=Jn, Ep=Ep, Em=Em, c=c, Cp=Cp)
        print(f"  α = {name}, n = {n}: J_n ∈ {AY.fmt(Jn, 24)}")
        print(f"      c ρ^n ∈ {AY.fmt(crn, 24)};  E_n^+ ≤ {mp.nstr(mpf(Ep.b), 3)}, E_n^- ≤ {mp.nstr(mpf(Em.b), 3)}")
        print(f"      C' = lim (J_m - c ρ^m/(ρ-1)) ∈ {AY.fmt(Cp, 24)}")
    sig = sigma_inf(q, x0, vals, N, a)
    return rho, kap, out, sig


def float_RS(m, a):
    """(R_α(m), σ_m) by mpmath quadrature in log s, unit breakpoints (consistency only, not certified).

    Head (0, e^lo], lo = -m ln5 - 40: r_m ≈ r_m(0), 1 - s r_m ≈ 1; tail [1e40, ∞): r_m ≈ 1/s, 1 - s r_m ≈ 2/s."""
    lo = -m * mp.log(5) - 40
    hi = 40 * mp.log(10)
    pts = [lo + k for k in range(int(hi - lo) + 1)] + [hi]
    fr = lambda x: mp.exp((1 - a) * x) * AY.r_mp(m, mp.exp(x))
    fw = lambda x: mp.exp((1 - a) * x) * (1 - mp.exp(x) * AY.r_mp(m, mp.exp(x)))
    h = mp.exp(lo) ** (1 - a) / (1 - a)
    J = mp.quad(fr, pts) + AY.r_mp(m, mpf(0)) * h + mp.exp(-a * hi) / a
    W = mp.quad(fw, pts) + h + 2 * mp.exp(-a * hi) / a
    return 2 * mp.sin(mp.pi * a) / mp.pi * J, mp.sin(mp.pi * a) / mp.pi * W


def main() -> int:
    iv.dps = 80
    mp.dps = 80
    print("Theorem AZ: supercritical decay law of the audit norm and its Mellin intercept")
    print("Assumptions: Theorems AQ (1)-(4), AR (1)-(3), AS (1)-(3), AT (1)-(3), AU (2)-(3), AY (2), all as stated. Not thrust.\n")
    P065 = iv.mpf(13) / 20
    K065 = 1 / (1 - iv.mpf(5) ** (-P065))
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

    t0 = time.time()
    levels = (N_LEVEL, N_LEVEL + 1, N_CROSS, N_CROSS + 1)
    q, x0, vals, N, N1 = AY.build(levels)
    print(f"  grid: {len(vals)} nodes, {N} cells, {N1} cells in (x_0, 1] ({time.time()-t0:.0f}s)")
    check(all(vals[k]["inf"].a > 0 and mpf(vals[k]["inf"].b) - mpf(vals[k]["inf"].a) < mpf(10) ** -60 * mpf(vals[k]["inf"].a)
              for k in vals), "node enclosures of r_∞ positive with relative width < 1e-60")
    check(all(vals[k][40].a <= vals[k][41].b and vals[k][41].a <= vals[k][60].b and vals[k][60].a <= vals[k][61].b
              and vals[k][61].a <= vals[k]["inf"].b for k in vals), "node order r_40 ≤ r_41 ≤ r_60 ≤ r_61 ≤ r_∞ (AQ (3))")

    gam = iv.log(3) / iv.log(5)
    res = {}
    for name, a, lv in (("0.8", iv.mpf("0.8"), (N_LEVEL, N_CROSS)), ("0.75", iv.mpf("0.75"), (N_LEVEL,))):
        res[name] = certify_alpha(name, a, lv, q, x0, vals, N, N1, BH, P065, K065)
        print(f"  ({time.time()-t0:.0f}s)")

    # α = 0.8: assemble
    rho, kap, out, sig = res["0.8"]
    c60, c40 = out[60]["c"], out[40]["c"]
    check(max(mpf(c60.a), mpf(c40.a)) <= min(mpf(c60.b), mpf(c40.b)), "α=0.8: c enclosures from levels 60 and 40 overlap")
    Cp60, Cp40 = out[60]["Cp"], out[40]["Cp"]
    check(max(mpf(Cp60.a), mpf(Cp40.a)) <= min(mpf(Cp60.b), mpf(Cp40.b)), "α=0.8: C' enclosures from levels 60 and 40 overlap")
    E60 = max(mpf(out[60]["Ep"].b), mpf(out[60]["Em"].b))
    E40 = max(mpf(out[40]["Ep"].b), mpf(out[40]["Em"].b))
    rf = float(mpf(rho.a))
    check(E40 < rf ** 40 * 5 ** (-0.6 * 40) and E60 < rf ** 60 * 5 ** (-0.6 * 60),
          f"α=0.8: defect bounds E_40 = {float(E40):.2e}, E_60 = {float(E60):.2e} below ρ^n 5^(-0.6 n)")
    c = iv.mpf([max(mpf(c60.a), mpf(c40.a)), min(mpf(c60.b), mpf(c40.b))])
    Cp = iv.mpf([max(mpf(Cp60.a), mpf(Cp40.a)), min(mpf(Cp60.b), mpf(Cp40.b))])
    cR = kap * c
    A = cR / (rho - 1)
    Ca = kap * Cp
    nu = Ca / A
    K = iv.sqrt(21) / 3 * sig / A
    print(f"\n  α = 0.8, ρ = 5^0.8/3 = {mp.nstr(mpf(rho.a), 20)}, decay ratio 3·5^-0.8 = {mp.nstr(mpf((1/rho).a), 20)}")
    for lab, v in (("c_R(0.8)", cR), ("A_0.8 = c_R/(ρ-1)", A), ("C_0.8 = κ M(0.2)", Ca), ("σ_∞(0.8)", sig),
                   ("ν_0.8 = C/A", nu), ("K_0.8 = √21 σ_∞/(3A)", K)):
        print(f"  {lab:24s} ∈ {AY.fmt(v, 22)}  (width {float(mpf(v.b) - mpf(v.a)):.1e})")
    check(all(mpf(v.b) - mpf(v.a) < 1e-20 for v in (cR, Ca, sig, K)), "α=0.8: certified widths of c_R, C, σ_∞, K below 1e-20")
    check(mpf(Ca.b) < 0 and mpf(nu.b) < 0, f"α=0.8: C_0.8 < 0 and ν_0.8 < 0 certified (C ≤ {mp.nstr(mpf(Ca.b), 12)})")
    check(mpf(sig.a) >= 2 * 4 ** (-0.8), "σ_∞(0.8) ≥ 2·4^-0.8 (Theorem AR (3))")

    # α = 0.75: c_R against Theorem AT's independent certificate
    rho75, kap75, out75, sig75 = res["0.75"]
    cR75 = kap75 * out75[60]["c"]
    check(mpf(cR75.a) <= mpf(AT_CR75[1]) and mpf(cR75.b) >= mpf(AT_CR75[0]),
          f"α=0.75: c_R ∈ {AY.fmt(cR75, 20)} meets Theorem AT's independent enclosure {AT_CR75}")
    A75 = cR75 / (rho75 - 1)
    C75 = kap75 * out75[60]["Cp"]
    K75 = iv.sqrt(21) / 3 * sig75 / A75
    print(f"  α = 0.75: A ∈ {AY.fmt(A75, 18)}, C ∈ {AY.fmt(C75, 18)}, σ_∞ ∈ {AY.fmt(sig75, 18)}, K ∈ {AY.fmt(K75, 18)}")
    check(mpf(C75.b) < 0, "α=0.75: C_0.75 < 0 certified")

    # Mellin pole: C_α ≈ -c_R(γ)/(ln5 (α-γ)); continuity of C_α + c_R(α)/(ρ-1) through γ (float consistency)
    mp.dps = 40
    g = mp.log(3) / mp.log(5)
    cRg = 2 * mp.sin(mp.pi * g) / mp.pi * mp.log(5) * mpf(AU_THETA[0])
    print("\n   m    R(m) - Aρ^m - C       (ρ^m+ν)||F|| - K       ρ^m||F|| - K         ρ^m(ρ^m||F|| - K)     [α = 0.8, 40-digit quadrature]")
    af = mpf("0.8")
    rhof = mpf(5) ** af / 3
    Af, Cf, nuf, Kf = (mpf(v.a) for v in (A, Ca, nu, K))
    rows = []
    for m in (5, 10, 20, 30):
        R, s = float_RS(m, af)
        F = mp.sqrt(21) / 3 * s / R
        rows.append((m, R - Af * rhof ** m - Cf, (rhof ** m + nuf) * F - Kf, rhof ** m * F - Kf))
        r = rows[-1]
        print(f"  {m:2d}   {float(r[1]):+.3e}           {float(r[2]):+.3e}            {float(r[3]):+.6e}       {float(rhof ** m * r[3]):+.8f}")
    check(all(abs(r[1]) < 5 * rhof ** r[0] * 5 ** (-0.65 * r[0]) for r in rows),
          "|R(m) - A ρ^m - C| < 5 ρ^m 5^(-0.65 m) at m = 5..30 (Theorem AZ (1))")
    check(all(abs(r[2]) < 5 * 5 ** (-0.65 * r[0]) for r in rows),
          "|(ρ^m + ν)||F(m,0.8)|| - K| < 5·5^(-0.65 m) at m = 5..30 (Theorem AZ (4))")
    lim = -nuf * Kf
    check(mpf(nu.b) < 0 and all(r[3] > 0 for r in rows) and abs(rhof ** 30 * rows[-1][3] - lim * rhof ** 30 / (rhof ** 30 + nuf)) < 1e-8,
          f"AZ.1: ρ^m||F|| - K > 0 at every m (approach from above) and ρ^m(ρ^m||F|| - K) → -ν K = {mp.nstr(lim, 10)} > 0 "
          f"(m = 30: {mp.nstr(rhof ** 30 * rows[-1][3], 10)}): error Θ((3·5^-0.8)^m), sign opposite to AY; claim falsified")
    obs = rows[3][1] / rows[2][1]
    print(f"  observation: (R(m) - Aρ^m - C) ratio from m = 20 to 30: {mp.nstr(obs, 6)} vs (ρ/3)^10 = {mp.nstr((rhof / 3) ** 10, 6)} "
          f"(the 3^-m correction of AT/AY times ρ^m); not proved")
    mono = [rows[i][3] > rows[i + 1][3] for i in range(len(rows) - 1)]
    check(all(mono), "AZ.1 witness: ρ^m||F(m,0.8)|| decreases toward K_0.8 on m = 5, 10, 20, 30")

    # pole and continuity through γ (float, quadrature at m = 29, 30; not part of the certificate)
    def G_of(al):
        rho_ = mpf(5) ** al / 3
        R29, _ = float_RS(29, al)
        R30, _ = float_RS(30, al)
        cR_ = (R30 - R29) / rho_ ** 29
        Cal = R30 - cR_ * rho_ ** 30 / (rho_ - 1)
        return cR_, Cal, Cal + cR_ / (rho_ - 1)
    eps_ = mpf("0.01")
    cp, Cpl, Gp = G_of(g + eps_)
    cm, Cmn, Gm = G_of(g - eps_)
    CG = mpf(AY_CG[0])
    print(f"  α = γ+0.01: C_α = {mp.nstr(Cpl, 10)}, (α-γ) C_α = {mp.nstr(eps_ * Cpl, 8)} vs -c_R(γ)/ln5 = {mp.nstr(-cRg / mp.log(5), 8)}")
    print(f"  α = γ-0.01: R_∞(α) = {mp.nstr(Cmn, 10)};  G(γ±0.01) = {mp.nstr(Gp, 10)}, {mp.nstr(Gm, 10)};  C_γ (AY) = {mp.nstr(CG, 10)}")
    check(Cpl < 0 < Cmn, "the intercept changes sign across γ: C_{γ+0.01} < 0 < R_∞(γ-0.01) (simple pole, positive residue)")
    check(abs(eps_ * Cpl + cRg / mp.log(5)) < 0.05, "(α-γ) C_α ≈ -c_R(γ)/ln5 at α = γ+0.01 (Mellin residue <Θ>, Theorem AU (2))")
    check(abs((Gp + Gm) / 2 - CG) < 1e-3 and abs(Gp - Gm) < 0.05,
          "C_α + c_R(α)/(ρ-1) is continuous through γ with value C_γ of Theorem AY (symmetric mean within 1e-3)")
    mp.dps = 80

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        for f in failures:
            print("   ", f)
        return 1
    print("All Theorem AZ / Kept Failure AZ.1 checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
