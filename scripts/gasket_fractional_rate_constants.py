#!/usr/bin/env python3
"""Certified refinement-rate constants c_R, c_S of the corner resolvent (Theorem AT).

Witness for Theorem AT and Kept Failure AT.1 in GASKET.md.

Notation as in Theorems AQ-AS: r_n(s) = ½ vᵀ(L_n + s)⁻¹ v, r_{n+1}(s) =
Φ_s(r_n(p(s))), p(s) = s(s+5), q = p⁻¹, r_0 = 1/(3+s), δ_{n+1} = r_{n+1} - r_n,
ρ = 5^α/3, μ = ρ/5, D_j(y) the Möbius secant slope of Theorem AR (1),
(3/5) D_j ≥ 1.  c_R, c_S are Theorem AR's constants:
    R_α(m+1) - R_α(m) = c_R ρ^m (1+o(1)),   σ_m - σ_{m+1} = c_S μ^m (1+o(1)).

Theorem AT (every α ∈ (0,1), every n ≥ 0).
  (1) Exact transfer.  Substituting s = q^k(t) (k = m - n) in AQ's Stieltjes
      integrals and iterating δ_{j+1}(s) = D_j(s) δ_j(p s):
        ρ^{-m} ∫ s^{-α} δ_{m+1} ds = ρ^{-n} ∫ t^{-α} δ_{n+1}(t) A_k(t)^{α} B_k(t) C_k(t) dt,
        μ^{-m} ∫ s^{1-α} δ_{m+1} ds = μ^{-n} ∫ t^{1-α} δ_{n+1}(t) A_k(t)^{α-1} B_k(t) C_k(t) dt,
      A_k = Π_{i=1}^k (1 + q^i t/5) ∈ [1, A_∞],  B_k = Π_{i=0}^{k-1} 5/√(25+4q^i t) ∈ [B_∞, 1],
      C_k = Π_{i=1}^k (3/5) D_{n+i}(q^i t) ∈ [1, C̄(t)],  C̄ = Π_{i≥1} (3/5) D̄(q^i t),
      D̄ the AQ-ceiling bound on every D_j.  Letting m → ∞ (Theorem AR (2)):
        c_R ∈ (2 sin πα/π) ρ^{-n} [∫ t^{-α} δ_{n+1} B_∞,  ∫ t^{-α} δ_{n+1} A_∞^{α} C̄],
        c_S ∈ (sin πα/π) μ^{-n} [∫ t^{1-α} δ_{n+1} A_∞^{α-1} B_∞,  ∫ t^{1-α} δ_{n+1} C̄].
      The bracket width is O(5^{-0.65 n}) relative; hence (derived)
      ρ^{-n}(R(n+1)-R(n)) = c_R + O(5^{-0.65 n}), μ^{-n}(σ_n-σ_{n+1}) = c_S + O(5^{-0.65 n}),
      and at γ: R_γ(n) = c_R(γ) n + C_γ + O(5^{-0.65 n}).
  (2) Certified (n = 40; outward-rounded mpmath.iv, 60 digits; r_n, r_{n+1}
      completely monotone, so Theorem AS's interpolation sandwich bounds
      ∫ s^{-α} δ_{n+1} cellwise; δ_{n+1} ≤ Π D̄ δ_1 beyond s = 1.2e4·5^{-40};
      an explicit doubly-exponential tail beyond s = 1):
        c_R(0.45), c_S(0.45), c_R(γ), c_S(γ), and c_R, c_S at 0.25, 0.60, 0.75,
      with widths ≤ 2e-17; and c_F = (√21/3) σ_∞ c_R/(R_∞²(1-ρ)) at 0.45
      (the leading coefficient of ||F(n,0.45)|| - L) from Theorem AS's enclosures.

Kept Failure AT.1.  Claim: r_∞(s) ~ C s^{γ-1} as s → 0 with one constant C,
so that the slope at γ is the pure-power value (2 sin πγ/π)·ln5·C.  False:
g(s) = s^{1-γ} r_∞(s) has q-orbit limits Θ(u) = lim_k g(q^k u) that are
certified to differ at two points (by ≈ 1.4e-5), so lim_{s→0} g(s) does not
exist; s^{1-γ} r_∞ is asymptotically log-periodic.

Conjecture AT.2 (observation, not proved): c_R(γ) = (2 sin πγ/π)·ln5·<Θ>, the
log-mean of Θ over one period; a 32-point trapezoid agrees with the certified
c_R(γ) to 5e-22.

Not thrust. Does not touch Stage 1 or Stage 2.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_limit_digits import lagrange_weights, node_bounds, power_ceiling  # noqa: E402

failures: list[str] = []

N_LEVEL = 40                 # n in Theorem AT (1)
Q_NUM, Q_DEN = 33, 32        # geometric node ratio
ORD = 12                     # floor b = u = 12 (even), ceiling b = u = 13 (odd)
EX0 = 110                    # head cut x_0 = 5^-n 2^-110
T_SPLIT = 12000              # interpolation up to 5^-n T_SPLIT, product bound above, crude tail above 1
TINY = mpf(10) ** -90        # orbit cut for the infinite products
AS_R45 = ("2.51594623787706", "2.51594623787714")      # Theorem AS enclosures at α = 0.45
AS_S45 = ("1.35247562325228", "1.35247562325233")
AS_L45 = ("0.82113862713346", "0.82113862713352")
AR_TABLE = {"0.25": (0.37699178, 1.082999), "0.45": (0.49780265, 0.927091), "0.60": (0.52084941, 0.642994),
            "gamma": (0.51716583, 0.481399), "0.75": (0.50798513, 0.357291)}   # Theorem AR float table


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


P065 = None
K065 = None
BH = None


def qf(x):
    return 2 * x / (5 + iv.sqrt(25 + 4 * x))


def r_exact(n: int, x):
    """r_n(x) by the exact recursion of Theorem AQ (1), in interval arithmetic."""
    orb = [x]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    r = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        o = orb[k]
        r = (2 + o) * (5 + o) * r / ((6 + o) - 2 * o * (5 + o) * r)
    return r


def zr_upper(z_hi) -> mpf:
    """Upper bound of z·r_m(z) for every m and every z ≤ z_hi (s·r is increasing; r_m ≤ r_∞)."""
    z = iv.mpf(z_hi)
    if z_hi <= mpf("0.01"):
        u = BH * z ** P065                       # AQ (4) power ceiling, α-free
    else:
        _, hi = node_bounds(z)                   # AQ (4): r_∞ ≤ T^k(1/σ)
        u = z * hi
    return min(mpf(u.b), mpf(1))


def dbar(yl: mpf, yh: mpf):
    """sup of D_j(y), all j ≥ 1, y ∈ [yl, yh]: (2+y)(5+y)(6+y)/den², den ≥ max(6+y-2·sup z r(z), 4+y), z = p(y)."""
    Yl, Yh = iv.mpf(yl), iv.mpf(yh)
    zh = mpf((Yh * (Yh + 5)).b)
    U = zr_upper(zh)
    den = max(mpf(((6 + Yl) - 2 * iv.mpf(U)).a), mpf((4 + Yl).a))
    return (2 + Yh) * (5 + Yh) * (6 + Yh) / iv.mpf(den) ** 2


def w_data(xlo, xhi):
    """B_∞(xhi) (lower), A_∞(xhi)^{-1} = Π 5/(5+q^i) (lower), C̄ over the cell (upper)."""
    lam_b = 5 / iv.sqrt(25 + 4 * xhi)
    lam_psi = iv.mpf(1)
    cbar = iv.mpf(1)
    ylo, yhi = iv.mpf(xlo), iv.mpf(xhi)
    while True:
        ylo, yhi = qf(ylo), qf(yhi)
        if yhi.b < TINY:
            break
        lam_b *= 5 / iv.sqrt(25 + 4 * yhi)
        lam_psi *= 5 / (5 + yhi)
        cbar *= iv.mpf(3) / 5 * dbar(mpf(ylo.a), mpf(yhi.b))
    yt = iv.mpf(yhi.b)      # tails i ≥ I, y_{i+1} ≤ y_i/5, y_I ≤ 1e-90 ≤ 1e-3
    lam_b *= iv.exp(-yt / 10)
    lam_psi *= iv.exp(-yt / 4)
    cbar *= iv.exp(2 * yt ** P065 * K065)
    return iv.mpf(lam_b.a), iv.mpf(lam_psi.a), iv.mpf(cbar.b)


def delta_sup(n: int, xlo, xhi) -> mpf:
    """sup of δ_{n+1} on [xlo, xhi] via δ_{n+1}(s) = Π_{k<n} D_{n-k}(p^k s) · δ_1(p^n s) (Theorem AR (1))."""
    Yl, Yh = iv.mpf(xlo.a), iv.mpf(xhi.b)
    prod = iv.mpf(1)
    for _ in range(n):
        prod *= dbar(mpf(Yl.a), mpf(Yh.b))
        Zl, Zh = Yl * (Yl + 5), Yh * (Yh + 5)
        Yl, Yh = iv.mpf(Zl.a), iv.mpf(Zh.b)
    u = Yl
    return mpf((prod * (u + 6) / ((u + 3) * (u * u + 7 * u + 9))).b)


def build(n: int):
    q = iv.mpf(Q_NUM) / Q_DEN
    x0 = iv.mpf(5) ** (-n) * iv.mpf(2) ** (-EX0)
    xs = iv.mpf(5) ** (-n) * T_SPLIT
    M = ORD + 1
    vals, k = {}, -M
    while True:
        x = x0 * q ** k
        vals[k] = (x, r_exact(n, x), r_exact(n + 1, x))
        if k >= M and vals[k - M][0].a > xs.b:
            break
        k += 1
    N = k - M
    cells = [w_data(vals[c][0], vals[c + 1][0]) for c in range(N)]
    head_w = w_data(iv.mpf(0), vals[0][0])   # λ's at x_0 (monotone); C̄ over [0, x_0]
    mid = []
    xc = vals[N][0]
    while xc.a < 1:
        xn = xc * q
        mid.append((xc, xn, delta_sup(n, xc, xn), w_data(xc, xn)))
        xc = xn
    S = iv.mpf(xc.a)
    lg = iv.log(S)                               # log p^n(S) ≥ ... (p(y) = y(y+5) ≥ y²)
    for _ in range(n):
        lg = lg + iv.log(iv.exp(lg) + 5) if lg.b < 60 else 2 * lg
    tail_log10 = mpf(((iv.log(iv.mpf("6e10")) + n * iv.log(iv.mpf("1.2")) + (iv.log(S) - lg) / 2) / iv.log(10)).b)
    return q, vals, N, cells, head_w, mid, tail_log10


def certify(n: int, a, kind: str, data):
    """Enclosure of K = ∫ t^{e} δ_{n+1} W dt (e = -α for R, 1-α for S) with the bracket weights of (1)."""
    q, vals, N, cells, head_w, mid, tail_log10 = data
    expo = a if kind == "R" else a - 1           # lagrange_weights integrates y^{-expo}
    Wf = lagrange_weights(q, list(range(1 - ORD, 1 + ORD)), expo)
    Wc = lagrange_weights(q, list(range(-ORD, 2 + ORD)), expo)

    def wlo_whi(wd):
        lb, lp, cb = wd
        if kind == "R":
            return lb, lp ** (-a) * cb
        return lp ** (1 - a) * lb, cb

    lo, hi, gap = iv.mpf(0), iv.mpf(0), iv.mpf(0)
    for c in range(N):
        sc = vals[c][0] ** (1 - expo)
        fl1 = sum((Wf[i] * vals[c + 1 - ORD + i][2] for i in range(2 * ORD)), iv.mpf(0))
        fl0 = sum((Wf[i] * vals[c + 1 - ORD + i][1] for i in range(2 * ORD)), iv.mpf(0))
        cu1 = sum((Wc[i] * vals[c - ORD + i][2] for i in range(2 * ORD + 2)), iv.mpf(0))
        cu0 = sum((Wc[i] * vals[c - ORD + i][1] for i in range(2 * ORD + 2)), iv.mpf(0))
        dlo = max(mpf((sc * (fl1 - cu0)).a), mpf(0))
        dhi = mpf((sc * (cu1 - fl0)).b)
        wl, wh = wlo_whi(cells[c])
        lo += iv.mpf(dlo) * iv.mpf(wl.a)
        hi += iv.mpf(dhi) * iv.mpf(wh.b)
        gap += iv.mpf(dhi) - iv.mpf(dlo)
    # head (0, x_0]: r_m decreasing, r_m(0) = (5/3)^m/3, so r_{n+1}(x_0) - r_n(0) ≤ δ_{n+1} ≤ r_{n+1}(0) - r_n(x_0)
    x0 = vals[0][0]
    h1 = x0 ** (1 - expo) / (1 - expo)
    r0n, r0n1 = (iv.mpf(5) / 3) ** n / 3, (iv.mpf(5) / 3) ** (n + 1) / 3
    dl = max(mpf((vals[0][2] - r0n).a), mpf(0))
    dh = mpf((r0n1 - vals[0][1]).b)
    wl, wh = wlo_whi(head_w)
    lo += iv.mpf(dl) * iv.mpf(wl.a) * iv.mpf(h1.a)
    head = iv.mpf(dh) * wh * h1 - iv.mpf(dl) * iv.mpf(wl.a) * iv.mpf(h1.a)
    hi += iv.mpf(dh) * wh * h1
    msum = iv.mpf(0)
    for xc, xn, ds, wd in mid:
        msum += iv.mpf(ds) * wlo_whi(wd)[1] * (xn ** (1 - expo) - xc ** (1 - expo)) / (1 - expo)
    hi += msum
    hi += iv.mpf(10) ** iv.mpf(tail_log10)
    return iv.mpf([lo.a, hi.b]), {"head": mpf(head.b), "mid": mpf(msum.b), "gap": mpf(gap.b)}


def r_mp(n, s):
    orb = [s]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    x = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        o = orb[k]
        x = (2 + o) * (5 + o) * x / ((6 + o) - 2 * o * (5 + o) * x)
    return x


def scaled_increment(m: int, a: mpf, kind: str):
    """ρ^{-m}(R(m+1)-R(m)) or μ^{-m}(σ_m-σ_{m+1}) by mpmath quadrature of δ_{m+1} (consistency only)."""
    e = -a if kind == "R" else 1 - a
    lam = mpf(5) ** (-m)
    pts = [0] + [lam * mpf(10) ** k for k in range(-12, 7)]
    I = mp.quad(lambda s: s ** e * (r_mp(m + 1, s) - r_mp(m, s)), pts)
    if kind == "R":
        return 2 * mp.sin(mp.pi * a) / mp.pi * I / (mpf(5) ** a / 3) ** m
    return mp.sin(mp.pi * a) / mp.pi * I / (mpf(5) ** (a - 1) / 3) ** m


def transfer_identity_check(a: mpf, n: int, k: int) -> mpf:
    """Relative error of the exact transfer (1) at finite m = n + k (mpmath, not interval)."""
    m = n + k
    lhs = mp.quad(lambda s: s ** (-a) * (r_mp(m + 1, s) - r_mp(m, s)),
                  [0] + [mpf(5) ** (-m) * mpf(10) ** j for j in range(-10, 8)]) / (mpf(5) ** a / 3) ** m

    def W(t):
        A, Bk, C = mpf(1), mpf(1), mpf(1)
        y = t
        for i in range(k):
            Bk *= 5 / mp.sqrt(25 + 4 * y)
            y = 2 * y / (5 + mp.sqrt(25 + 4 * y))      # y = q^{i+1} t
            A *= 1 + y / 5
            xa, xb = r_mp(n + i, y * (y + 5)), r_mp(n + i + 1, y * (y + 5))
            den = lambda z: (6 + y) - 2 * y * (5 + y) * z
            C *= mpf(3) / 5 * (2 + y) * (5 + y) * (6 + y) / (den(xa) * den(xb))
        return A ** a * Bk * C

    rhs = mp.quad(lambda t: t ** (-a) * (r_mp(n + 1, t) - r_mp(n, t)) * W(t),
                  [0] + [mpf(5) ** (-n) * mpf(10) ** j for j in range(-10, 8)]) / (mpf(5) ** a / 3) ** n
    return abs(lhs / rhs - 1)


def theta_enclosure(s):
    """Θ(p^k(s)) = lim_j g(q^j p^k s) ∈ g(s)·exp([-dl, dh]), g = s^{1-γ} r_∞(s) (Kept Failure AT.1)."""
    gam = iv.log(3) / iv.log(5)
    lo, hi = node_bounds(s)
    g = iv.mpf([(s ** (1 - gam) * lo).a, (s ** (1 - gam) * hi).b])
    y = qf(s)       # log g(y)/g(p y) ∈ [-(1-γ) y/5, 1.81 y^0.65] for y ≤ 1e-3; y_{j+1} ≤ y_j/5
    dl = (1 - gam) * y / 4
    dh = 2 * y ** P065 * K065
    return iv.mpf([(g * iv.exp(-dl)).a, (g * iv.exp(dh)).b])


def fmt(x, d=18):
    """Outward-rounded decimal enclosure with d digits after the point."""
    sc = mpf(10) ** d
    lo = mp.floor(mpf(x.a) * sc)
    hi = mp.ceil(mpf(x.b) * sc)
    return f"[{mp.nstr(lo / sc, d + 3)}, {mp.nstr(hi / sc, d + 3)}]"


def main() -> int:
    global P065, K065, BH
    iv.dps = 60
    mp.dps = 40
    P065 = iv.mpf(13) / 20
    K065 = 1 / (1 - iv.mpf(5) ** (-P065))
    print("Theorem AT: certified refinement-rate constants c_R, c_S (exact transfer + CM sandwich)")
    print("Assumptions: Theorems AQ (1)-(4), AR (1)-(2), AS (2)-(3) as stated. Not thrust.\n")

    # small-y lemma constants used for the infinite-product tails (y ≤ 1e-3)
    B, okB, c1 = power_ceiling()
    BH = iv.mpf(B.b)
    check(okB, f"AQ (4) power ceiling re-certified: r_m ≤ B s^-0.35 on (0,0.01], B = {mp.nstr(mpf(B.b), 9)} (C1 lhs {c1:.6f} < 6)")
    y = iv.mpf("0.001")
    eps = 2 * BH * (y * (5 + y)) ** P065 / 6
    lem = (y * (3 * y + 16) / (3 * (6 + y)) + 2 * eps / (1 - eps)) / y ** P065
    check(lem.b < 1.82 and eps.b < 0.01,
          f"tail lemma: log((3/5)D̄(y)) ≤ y + 2ε/(1-ε) ≤ 1.82 y^0.65 for y ≤ 1e-3 (ratio {mp.nstr(mpf(lem.b), 5)}, ε ≤ {mp.nstr(mpf(eps.b), 3)})")
    # crude-tail constant: ln W ≤ 2.668(9 + ln s) + 2 ln s + 0.04 ≤ 24.06 + 4.67 ln s, so W ≤ 3e10 s^5 (s ≥ 1)
    cst = 2.668 * 9 + 0.04
    check(math.log(3.6) + 2 * math.log(2) < 2.668 and 2.1 * 0.001 ** 0.65 / (1 - 5 ** -0.65) < 0.04
          and math.exp(cst) < 3e10, f"crude tail lemma constants: e^{cst:.2f} < 3e10")

    # (1) exact transfer at finite m, mpmath
    for (al, n, k) in (("0.45", 2, 3), ("0.75", 3, 2)):
        err = transfer_identity_check(mpf(al), n, k)
        check(err < mpf(10) ** -15, f"transfer identity (1) at α={al}, n={n}, m={n+k}: rel err {float(err):.1e}")

    t0 = time.time()
    data = build(N_LEVEL)
    q, vals, N, cells, head_w, mid, tail_log10 = data
    print(f"  n = {N_LEVEL}: {len(vals)} nodes (ratio {Q_NUM}/{Q_DEN}), {N} interpolation cells, {len(mid)} product-bound cells, "
          f"crude tail < 10^{mp.nstr(tail_log10, 6)} ({time.time()-t0:.0f}s)")
    check(all(0 < vals[k][1].a <= vals[k][2].a for k in vals), "node enclosures: 0 < r_n ≤ r_{n+1} (AQ (3))")
    check(all(c[0].a > 0 and c[1].a > 0 and c[2].b >= 1 for c in cells), "bracket weights: B_∞ > 0, A_∞^{-1} > 0, C̄ ≥ 1")

    gam = iv.log(3) / iv.log(5)
    alphas = [("0.25", iv.mpf("0.25")), ("0.45", iv.mpf("0.45")), ("0.60", iv.mpf("0.60")), ("gamma", gam), ("0.75", iv.mpf("0.75"))]
    res = {}
    print("\n  α        c_R(α) certified                                   c_S(α) certified")
    for name, a in alphas:
        rho = iv.mpf(5) ** a / 3
        KR, iR = certify(N_LEVEL, a, "R", data)
        KS, iS = certify(N_LEVEL, a, "S", data)
        cR = 2 * iv.sin(iv.pi * a) / iv.pi * rho ** (-N_LEVEL) * KR
        cS = iv.sin(iv.pi * a) / iv.pi * (rho / 5) ** (-N_LEVEL) * KS
        res[name] = (cR, cS)
        print(f"  {name:6s} {fmt(cR)}  {fmt(cS)}")
        wR, wS = mpf(cR.b) - mpf(cR.a), mpf(cS.b) - mpf(cS.a)
        check(wR < 1e-15 and wS < 1e-15, f"α={name}: widths {float(wR):.1e} (c_R), {float(wS):.1e} (c_S) < 1e-15 "
              f"[head width {float(iR['head']):.0e}, product cells {float(iR['mid']):.0e}, sandwich gap {float(iR['gap']):.0e}]")
        fR, fS = AR_TABLE[name]
        check(abs(float(mpf(cR.a)) - fR) < 1e-8, f"α={name}: AR float c_R = {fR} agrees to 1e-8")
        print(f"         AR float c_S = {fS}; certified - AR = {float(mpf(cS.a)) - fS:+.2e}")

    # consistency: scaled increments at finite m approach the certified constants (rate ~ 3^-m, inside the O(5^-0.65 m) bound)
    a = mpf("0.45")
    cR45 = mpf(res["0.45"][0].a)
    cS45 = mpf(res["0.45"][1].a)
    devs = []
    for m in (4, 8, 12, 16):
        dR = scaled_increment(m, a, "R") - cR45
        dS = scaled_increment(m, a, "S") - cS45
        devs.append((m, dR, dS))
        print(f"  m = {m:2d}: ρ^-m ΔR_m - c_R = {float(dR):+.3e},  μ^-m Δσ_m - c_S = {float(dS):+.3e}")
    check(all(abs(d[1]) < 5 * 5 ** (-0.65 * d[0]) and abs(d[2]) < 5 * 5 ** (-0.65 * d[0]) for d in devs),
          "α=0.45: |scaled increment - certified constant| < 5·5^(-0.65 m) at m = 4, 8, 12, 16")
    check(all(abs(devs[i + 1][1]) < abs(devs[i][1]) for i in range(3)), "α=0.45: R deviations decrease in m")

    # corollary: leading coefficient of ||F(n,0.45)|| - L
    cR = res["0.45"][0]
    rho = iv.mpf(5) ** iv.mpf("0.45") / 3
    Rinf = iv.mpf(list(AS_R45))
    Sinf = iv.mpf(list(AS_S45))
    cF = iv.sqrt(21) / 3 * Sinf * cR / (Rinf ** 2 * (1 - rho))
    cRinf = cR / (1 - rho)
    print(f"\n  α=0.45: R_∞ - R(n) ~ c_R/(1-ρ) ρ^n, c_R/(1-ρ) ∈ {fmt(cRinf, 14)}; ||F(n)|| - L ~ c_F ρ^n, c_F ∈ {fmt(cF, 14)}")
    mp.dps = 30
    Fm = lambda m: mp.sqrt(21) / 3 * (mp.sin(mp.pi * a) / mp.pi * mp.quad(lambda s: s ** (-a) * (1 - s * r_mp(m, s)), [0] + [mpf(10) ** k for k in range(-30, 31, 2)]) + 2 * mpf(10) ** (-30 * a) / a * mp.sin(mp.pi * a) / mp.pi) \
        / (2 * mp.sin(mp.pi * a) / mp.pi * (mp.quad(lambda s: s ** (-a) * r_mp(m, s), [0] + [mpf(10) ** k for k in range(-30, 31, 2)]) + mpf(10) ** (-30 * a) / a))
    Lmid = (mpf(AS_L45[0]) + mpf(AS_L45[1])) / 2
    rhof = mpf(5) ** a / 3
    ratios = [((Fm(m) - Lmid) / rhof ** m) for m in (20, 30)]
    print(f"  (||F(m)|| - L)/ρ^m at m = 20, 30: {[mp.nstr(x, 8) for x in ratios]}")
    check(abs(ratios[1] - mpf(cF.a)) < 1e-4 * mpf(cF.a), "α=0.45: (||F(30)|| - L)/ρ^30 matches certified c_F to 1e-4")
    mp.dps = 40

    # Kept Failure AT.1: s^{1-γ} r_∞(s) is not asymptotically constant
    th = []
    for j in range(8):
        s = iv.mpf(5) ** (-40) * iv.mpf(5) ** (iv.mpf(j) / 8)
        th.append(theta_enclosure(s))
    lo_max = max(mpf(t.a) for t in th)
    hi_min = min(mpf(t.b) for t in th)
    width = max(mpf(t.b) - mpf(t.a) for t in th)
    print(f"\n  Kept Failure AT.1: Θ at 8 points of one log5-period: " + ", ".join(mp.nstr(mpf(t.a), 12) for t in th))
    check(lo_max - hi_min > 1e-5 and width < 1e-15,
          f"AT.1: certified Θ_max - Θ_min > {float(lo_max - hi_min):.3e} (enclosure widths < {float(width):.0e}); "
          "lim s^(1-γ) r_∞(s) does not exist: pure-power claim falsified")
    # observation: c_R(γ) = (2 sin πγ/π) ln5 <Θ> (log-mean, 32-point trapezoid on a period; not proved)
    gm = mp.log(3) / mp.log(5)
    pts = []
    for j in range(32):
        s = mpf(5) ** (-40) * mpf(5) ** (mpf(j) / 32)
        pts.append(s ** (1 - gm) * r_mp(160, s))
    slope = 2 * mp.sin(mp.pi * gm) / mp.pi * mp.log(5) * sum(pts) / 32
    cRg = res["gamma"][0]
    print(f"  observation: (2 sin πγ/π) ln5 <Θ> = {mp.nstr(slope, 18)} vs certified c_R(γ) ∈ {fmt(cRg, 18)} "
          f"(difference {float(slope - mpf(cRg.a)):.1e})")
    check(abs(slope - mpf(cRg.a)) < 1e-17, "observation: log-mean slope matches certified c_R(γ) to 1e-17 (not a proof; Conjecture AT.2)")

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        for f in failures:
            print("   ", f)
        return 1
    print("All Theorem AT / Kept Failure AT.1 checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
