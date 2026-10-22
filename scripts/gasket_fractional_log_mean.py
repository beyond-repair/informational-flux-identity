#!/usr/bin/env python3
"""Log-mean slope at log3/log5 and the Mellin poles of the corner resolvent (Theorem AU).

Witness for Theorem AU (proof of Conjecture AT.2) and Kept Failure AU.1 in GASKET.md.

Notation as in Theorems AQ-AT: r_n(s) = ½ vᵀ(L_n + s)⁻¹ v, r_{n+1}(s) = Φ_s(r_n(p(s))),
p(s) = s(s+5), q = p⁻¹, r_∞ = lim r_n, γ = log3/log5, g(s) = s^{1-γ} r_∞(s),
Θ(u) = lim_k g(q^k u) (Kept Failure AT.1: exists, Θ∘q = Θ),
ψ(s) = lim_k 5^k q^k(s) = s/A_∞(s) (Koenigs coordinate of q: ψ∘q = ψ/5),
θ = Θ∘ψ⁻¹ (exactly 5-periodic: θ(x/5) = θ(x)), ω = 2π/ln5,
<Θ> = (1/ln5) ∫_{x/5}^{x} θ dlog y,  θ̂_k = (1/ln5) ∫_{x/5}^{x} θ(y) y^{-iωk} dlog y = a_k - i b_k.

Theorem AU.
  (1) Slope formula (Conjecture AT.2 PROVED).  ∫_0^∞ s^{-γ} r_m ds = m ln5 <Θ> + O(1), so
        c_R(γ) = (2 sin πγ/π) ln5 <Θ>.
      Proof pieces: (0,5^-m] ≤ r_m(0) 5^{-m(1-γ)}/(1-γ) = 1/(3(1-γ)); [1,∞) ≤ 1/γ;
      ∫_{5^-m}^1 s^{-γ}(r_∞-r_m) ≤ 5^{mε} Σ_{j≥m} ∫ s^{-(γ-ε)} δ_{j+1} = O(1) (AR at γ-ε);
      ∫_{5^-m}^1 g dlog s = m ln5 <Θ> + O(1) (|g-Θ| = O(s^0.65), ψ change of variables).
  (2) Mellin poles (derived).  M(z) = ∫_0^∞ s^{z-1} r_∞ ds, holomorphic on 1-γ < Re z < 1,
      continues meromorphically to 1-γ-0.65 < Re z < 1 with at most simple poles, only at
      z_k = 1-γ+iωk, residue θ̂_{-k} = conj(θ̂_k).  R_∞(α) = (2 sin πα/π) M(1-α), so
      lim_{α↑γ} (γ-α) R_∞(α) = (2 sin πγ/π) <Θ> = c_R(γ)/ln5.
  (3) Certified (outward-rounded mpmath.iv): complete-monotone sandwich (Theorem AS (2)) for
      ∫ s^{-γ} r_∞ (and against the nonnegative weights 1 ± cos, 1 ± sin) on one period at
      u = 5^-140, plus explicit Θ/g, ψ/s and endpoint corrections:
        <Θ>    ∈ [0.6009542511687736078237051500999942, ...9943],
        c_R(γ) ∈ [0.5171658346305052439259804856731271, ...1273]
      (inside Theorem AT's independent [0.517165834630505243, ...245]); the Fourier
      coefficients θ̂_1..θ̂_4 to ±2e-35, all certified nonzero: genuine complex poles
      of M at 1-γ ± iωk, k = 1..4.

Kept Failure AU.1.  Claim: the log-periodic part is a single harmonic, Θ = <Θ> + 2 Re(θ̂_1 x^{iω}).
False: |θ̂_2| is certified positive, so the k = 2 pole is present too.

Not thrust. Does not touch Stage 1 or Stage 2.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

from mpmath import iv, mp, mpf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_limit_digits import lagrange_weights, power_ceiling  # noqa: E402

failures: list[str] = []

K_PER = 128                 # cells per log5-period (node ratio 5^(1/K))
ORD = 14                    # floor b = u = 14 (even), ceiling b = u = 15 (odd)
E_DEPTH = 140               # period [5^-E-1, 5^-E]
CAP = mpf(10) ** 90         # orbit cap in the AQ (4) node enclosure
KMAX = 4                    # Fourier modes certified
AT_CRG = ("0.517165834630505243", "0.517165834630505245")   # Theorem AT enclosure of c_R(γ)
AT1_THETA = ("0.600954332", "0.600949127", "0.600946924", "0.600949013",
             "0.600954170", "0.600959375", "0.600961578", "0.600959489")   # Kept Failure AT.1 (5^-40 5^{j/8})


def check(ok: bool, msg: str) -> None:
    print(("PASS " if ok else "FAIL ") + msg)
    if not ok:
        failures.append(msg)


def node_bounds(x):
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


def lagrange_trig_weights(Q, ms, a, tau):
    """(C_m, S_m) = ∫_1^Q y^{-a} ℓ_m(y) (cos, sin)(τ ln y) dy, Lagrange basis on y = Q^m."""
    ys = [Q ** m for m in ms]
    lnQ = iv.log(Q)
    cph, sph = iv.cos(tau * lnQ), iv.sin(tau * lnQ)
    out = []
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
        cr, ci = iv.mpf(0), iv.mpf(0)
        for d, c in enumerate(coef):
            e = d + 1 - a                     # ∫_1^Q y^{e-1+iτ} dy = (Q^e e^{iτ lnQ} - 1)/(e + iτ)
            Qe = Q ** e
            nr, ni = Qe * cph - 1, Qe * sph
            dd = e * e + tau * tau
            cr += c * (nr * e + ni * tau) / dd
            ci += c * (ni * e - nr * tau) / dd
        out.append((cr / den, ci / den))
    return out


def period_integrals(E: int, shift, gam, kmax: int):
    """Enclosures over s ∈ [x0, 5 x0], x0 = 5^{-E-1+shift}, of
       I0 = ∫ s^{-γ} r_∞ ds and, for k ≤ kmax, Ic_k = ∫ s^{-γ} r_∞ cos(ωk ln s) ds, Is_k (sin)."""
    Q = iv.mpf(5) ** (iv.mpf(1) / K_PER)
    x0 = iv.mpf(5) ** (-E - 1 + shift)
    M = ORD + 1
    vals = {k: node_bounds(x0 * Q ** k) for k in range(-M, K_PER + M + 1)}
    msf, msc = list(range(1 - ORD, 1 + ORD)), list(range(-ORD, 2 + ORD))
    Wf, Wc = lagrange_weights(Q, msf, gam), lagrange_weights(Q, msc, gam)
    om = 2 * iv.pi / iv.log(5)
    trig = {k: (lagrange_trig_weights(Q, msf, gam, om * k), lagrange_trig_weights(Q, msc, gam, om * k))
            for k in range(1, kmax + 1)}
    lnx0 = iv.log(x0)

    def sandwich(wf_of, wc_of):
        lo, hi = iv.mpf(0), iv.mpf(0)
        for c in range(K_PER):
            sc = (x0 * Q ** c) ** (1 - gam)
            fl = sum((wf_of(c, i) * vals[c + 1 - ORD + i] for i in range(2 * ORD)), iv.mpf(0))
            cu = sum((wc_of(c, i) * vals[c - ORD + i] for i in range(2 * ORD + 2)), iv.mpf(0))
            lo += sc * iv.mpf(fl.a)
            hi += sc * iv.mpf(cu.b)
        return iv.mpf([lo.a, hi.b])

    I0 = sandwich(lambda c, i: Wf[i], lambda c, i: Wc[i])
    out = {"I0": I0, "vals": vals, "x0": x0}
    for k in range(1, kmax + 1):
        tf, tc = trig[k]
        ph = {c: om * k * (lnx0 + c * iv.log(Q)) for c in range(K_PER)}
        cs = {c: (iv.cos(ph[c]), iv.sin(ph[c])) for c in ph}
        res = {}
        # cos(τ ln s) = cosφ cos(τ ln y) - sinφ sin(τ ln y); sin(τ ln s) = sinφ cos + cosφ sin
        for name, comb in (("c", lambda cc, ss, C, S: cc * C - ss * S), ("s", lambda cc, ss, C, S: ss * C + cc * S)):
            parts = []
            for sign in (1, -1):        # nonnegative weights 1 ± cos (resp. 1 ± sin)
                J = sandwich(lambda c, i: Wf[i] + sign * comb(cs[c][0], cs[c][1], tf[i][0], tf[i][1]),
                             lambda c, i: Wc[i] + sign * comb(cs[c][0], cs[c][1], tc[i][0], tc[i][1]))
                parts.append(sign * (J - I0))
            # two enclosures of the same number; intersect
            lo_, hi_ = max(mpf(parts[0].a), mpf(parts[1].a)), min(mpf(parts[0].b), mpf(parts[1].b))
            res[name] = iv.mpf([lo_, hi_])
        out[k] = res
    return out


def corrections(E: int, shift, gam, P065, K065, kmax: int, data):
    """Turn the s-integrals over [x0, 5x0] into enclosures of <Θ>, a_k, b_k (Θ/g, ψ/s, endpoint)."""
    ln5 = iv.log(5)
    u = data["x0"] * 5
    qu = 2 * u / (5 + iv.sqrt(25 + 4 * u))
    dl = (1 - gam) * qu / 4                       # Θ(s) ∈ g(s) exp([-dl, dh]) for s ≤ u (AT.1 lemma)
    dh = 2 * qu ** P065 * K065
    eps = u / 20                                  # 0 ≤ log A_∞ ≤ s/20, 0 ≤ s (log A_∞)' ≤ s/20
    I0 = data["I0"]
    # endpoint piece ∫_{q u}^{u/5} g dlog s ≤ sup g · log(u/(5 q u)) ≤ (u/5)^{1-γ} r_∞(q u) · u/25
    r_qu = node_bounds(iv.mpf(qu.a))
    end = (u / 5) ** (1 - gam) * iv.mpf(r_qu.b) * u / 25
    mean = iv.mpf([((I0 * iv.exp(-dl) * (1 - eps)) / ln5).a, (((I0 + end) * iv.exp(dh)) / ln5).b])
    om = 2 * iv.pi / ln5
    Ihi = iv.mpf(I0.b)
    eta = iv.mpf(max(mpf((iv.exp(dh) - 1).b), mpf((1 - iv.exp(-dl)).b)))
    coeffs = {}
    for k in range(1, kmax + 1):
        err = Ihi * (eta + iv.exp(dh) * (1 + om * k) * eps) / ln5 + end * iv.exp(dh) / ln5
        e = mpf(err.b)
        ak = data[k]["c"] / ln5
        bk = data[k]["s"] / ln5
        coeffs[k] = (iv.mpf([mpf(ak.a) - e, mpf(ak.b) + e]), iv.mpf([mpf(bk.a) - e, mpf(bk.b) + e]))
    return mean, coeffs, {"dh": mpf(dh.b), "eps": mpf(eps.b), "end": mpf(end.b)}


def r_mp(n, s):
    orb = [s]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    x = 1 / (3 + orb[-1])
    for k in range(n - 1, -1, -1):
        o = orb[k]
        x = (2 + o) * (5 + o) * x / ((6 + o) - 2 * o * (5 + o) * x)
    return x


def fmt(x, d=34):
    sc = mpf(10) ** d
    return f"[{mp.nstr(mp.floor(mpf(x.a) * sc) / sc, d + 2)}, {mp.nstr(mp.ceil(mpf(x.b) * sc) / sc, d + 2)}]"


def main() -> int:
    iv.dps = 120
    mp.dps = 120
    print("Theorem AU: log-mean slope at log3/log5 (Conjecture AT.2) and Mellin poles of r_∞")
    print("Assumptions: Theorems AQ (1)-(4), AR (1)-(2), AS (2)-(3), AT (4a) and the orbit lemma of")
    print("Kept Failure AT.1, all as stated. Not thrust.\n")
    P065 = iv.mpf(13) / 20
    K065 = 1 / (1 - iv.mpf(5) ** (-P065))
    gam = iv.log(3) / iv.log(5)

    # re-certify the ingredients of the Θ enclosure (AQ (4) power ceiling; AT.1 orbit lemma, y ≤ 1e-3)
    iv.dps = 60
    B, okB, c1 = power_ceiling()
    BH = iv.mpf(B.b)
    y = iv.mpf("0.001")
    epsl = 2 * BH * (y * (5 + y)) ** P065 / 6
    lem = (y * (3 * y + 16) / (3 * (6 + y)) + 2 * epsl / (1 - epsl)) / y ** P065
    iv.dps = 120
    check(okB, f"AQ (4) power ceiling re-certified: r_m ≤ B s^-0.35 on (0,0.01], B = {mp.nstr(mpf(B.b), 9)} (C1 lhs {c1:.6f} < 6)")
    check(lem.b < 1.82, f"AT.1 orbit lemma: log(g(y)/g(p y)) ≤ 1.82 y^0.65 for y ≤ 1e-3 (ratio {mp.nstr(mpf(lem.b), 5)})")
    # Koenigs bounds: q(x) ≤ x/5 and q'(x) ≤ 1/5, so 0 ≤ log A_∞(s) ≤ s/20 and 0 ≤ s (log A_∞)'(s) ≤ s/20
    xs = [iv.mpf(t) for t in ("1e-30", "0.001", "1", "100")]
    check(all((2 * x / (5 + iv.sqrt(25 + 4 * x))).b <= (x / 5).a and (1 / iv.sqrt(25 + 4 * x)).b <= mpf("0.2") for x in xs),
          "Koenigs bounds q(x) ≤ x/5, q'(x) ≤ 1/5 (spot check; both are immediate from q(x)(5+q(x)) = x)")

    t0 = time.time()
    data = period_integrals(E_DEPTH, iv.mpf(0), gam, KMAX)
    mean, coeffs, info = corrections(E_DEPTH, iv.mpf(0), gam, P065, K065, KMAX, data)
    print(f"  period [5^-{E_DEPTH+1}, 5^-{E_DEPTH}]: {K_PER} cells, {len(data['vals'])} nodes, sandwich b = {ORD}/{ORD+1} "
          f"({time.time()-t0:.0f}s); corrections: Θ/g ≤ e^{mp.nstr(info['dh'], 3)}, ψ/s ≤ {mp.nstr(info['eps'], 3)}, "
          f"endpoint ≤ {mp.nstr(info['end'], 3)}")
    check(all(v.a > 0 and mpf(v.b) - mpf(v.a) < mpf(10) ** -60 * mpf(v.a) for v in data["vals"].values()),
          "node enclosures of r_∞ positive with relative width < 1e-60")
    wmean = mpf(mean.b) - mpf(mean.a)
    print(f"\n  <Θ> ∈ {fmt(mean)}  (width {float(wmean):.1e})")
    check(wmean < 1e-30, "certified width of <Θ> < 1e-30")

    # independence of the base point: a second, shifted period (5^-E-1 · 5^0.37), shallower depth
    data2 = period_integrals(E_DEPTH - 20, iv.mpf("0.37"), gam, 0)
    mean2, _, _ = corrections(E_DEPTH - 20, iv.mpf("0.37"), gam, P065, K065, 0, data2)
    check(max(mpf(mean.a), mpf(mean2.a)) <= min(mpf(mean.b), mpf(mean2.b)),
          f"<Θ> from the shifted period [5^-{E_DEPTH-19}·5^0.37/5, ...] overlaps (well-defined mean; ψ-periodicity)")

    # (1) slope formula and agreement with Theorem AT
    pref = 2 * iv.sin(iv.pi * gam) / iv.pi
    cR = pref * iv.log(5) * mean
    AT = iv.mpf(list(AT_CRG))
    print(f"\n  (2 sin πγ/π) ln5 <Θ> ∈ {fmt(cR)}")
    print(f"  Theorem AT: c_R(γ)    ∈ [{AT_CRG[0]}, {AT_CRG[1]}]")
    check(mpf(cR.a) >= mpf(AT.a) and mpf(cR.b) <= mpf(AT.b),
          "Theorem AU (1) consistency: the log-mean enclosure lies inside Theorem AT's independent enclosure of c_R(γ)")
    resid = cR / iv.log(5)
    print(f"  (2) residue: lim_(α↑γ) (γ-α) R_∞(α) = c_R(γ)/ln5 ∈ {fmt(resid, 30)}")

    # (3) Fourier coefficients of θ
    print("\n  k   a_k = <θ cos(ωk ln x)>                         b_k = <θ sin(ωk ln x)>              |θ̂_k| ≥")
    mods = {}
    for k in range(1, KMAX + 1):
        ak, bk = coeffs[k]
        lo2 = lambda z: mpf(0) if z.a <= 0 <= z.b else min(abs(mpf(z.a)), abs(mpf(z.b)))
        mlo = mp.sqrt(lo2(ak) ** 2 + lo2(bk) ** 2)
        mods[k] = (mlo, mp.sqrt(max(abs(mpf(ak.a)), abs(mpf(ak.b))) ** 2 + max(abs(mpf(bk.a)), abs(mpf(bk.b))) ** 2))
        print(f"  {k}  {mp.nstr(mpf(ak.a), 22)} ± {float((mpf(ak.b)-mpf(ak.a))/2):.0e}   "
              f"{mp.nstr(mpf(bk.a), 22)} ± {float((mpf(bk.b)-mpf(bk.a))/2):.0e}   {mp.nstr(mlo, 8)}")
    for k in range(1, KMAX + 1):
        check(mods[k][0] > 0, f"θ̂_{k} ≠ 0 certified (|θ̂_{k}| ≥ {mp.nstr(mods[k][0], 6)}): M has a genuine pole at 1-γ ± {k}iω")
    check(all(mods[k][0] > 1e4 * mods[k + 1][1] for k in range(1, KMAX)),
          "observation: |θ̂_k| ≥ 1e4 |θ̂_(k+1)| for k = 1..3 (each mode ~1e-5.3 of the previous)")
    # Kept Failure AU.1: single-harmonic claim
    check(mods[2][0] > 0, f"AU.1: |θ̂_2| ≥ {mp.nstr(mods[2][0], 6)} > 0, so Θ is not <Θ> + one harmonic: single-harmonic claim falsified")

    # reconstruction against Kept Failure AT.1's eight Θ values at 5^-40·5^{j/8} (float check)
    m0 = mpf(mean.a)
    worst = mpf(0)
    for j, tv in enumerate(AT1_THETA):
        ph = 2 * mp.pi * (mpf(j) / 8)          # ω ln(5^-40 5^{j/8}) ≡ 2π j/8 (mod 2π); ψ(s)/s - 1 ~ 1e-29 there
        val = m0 + sum(2 * (mpf(coeffs[k][0].a) * mp.cos(k * ph) + mpf(coeffs[k][1].a) * mp.sin(k * ph)) for k in coeffs)
        worst = max(worst, abs(val - mpf(tv)))
    check(worst < 1e-9, f"Fourier reconstruction reproduces Kept Failure AT.1's eight Θ values to {float(worst):.1e} (9-digit record)")

    # observation: J_m - m ln5 <Θ> stays bounded (piece bounds of the proof) and converges (C_γ)
    mp.dps = 30
    g = mp.log(3) / mp.log(5)
    print("\n  m    J_m - m ln5 <Θ>   (J_m = ∫ s^-γ r_m ds, mpmath quadrature, 30 digits; observation)")
    diffs = []
    for m in (5, 10, 20, 30, 40):
        pts = [mpf(0)] + [mpf(5) ** (-m) * mpf(10) ** k for k in range(-14, 2)] + [mpf(10) ** k for k in range(-6, 31, 2)]
        pts = sorted(set(pts))
        Jm = mp.quad(lambda s: s ** (-g) * r_mp(m, s), pts) + mpf(10) ** (-30 * g) / g
        d = Jm - m * mp.log(5) * mpf(mean.a)
        diffs.append(d)
        print(f"  {m:2d}   {mp.nstr(d, 16)}")
    bound = 1 / (3 * (1 - g)) + 1 / g
    check(all(-50 < d < 50 for d in diffs), "J_m - m ln5 <Θ> bounded on m = 5..40 (Theorem AU (1))")
    check(abs(diffs[-1] - diffs[-2]) < 1e-10, "J_m - m ln5 <Θ> converges (|J_40 - J_30 - 10 ln5 <Θ>| < 1e-10)")
    Cg = 2 * mp.sin(mp.pi * g) / mp.pi * diffs[-1]
    print(f"  ⇒ C_γ = lim (R_γ(m) - c_R(γ) m) ≈ {mp.nstr(Cg, 14)} (float quadrature; not certified; piece bound {mp.nstr(bound, 6)})")
    mp.dps = 120

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED")
        for f in failures:
            print("   ", f)
        return 1
    print("All Theorem AU / Kept Failure AU.1 checks passed. Not thrust.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
