#!/usr/bin/env python3
r"""Theorem BW: AU's certified mode routine, re-run deeper (period at 5^-260) with a finer sandwich (b = u = 32/33), certifies
theta_k != 0 for every k <= KMAX (BV: k <= 11), <F> = <Theta> and c_R(gamma) to ~1e-113, and, with BV.3's tail, the
extremes of F (G* = max F, min F to ~1e-109; t_max, t_min to ~1e-103).

Setting (GASKET.md, Theorems AQ, AU, BN, BS-BV): r_inf = lim r_k, gamma = log3/log5, F(u) = u^(1-gamma) Q_inf(u) = theta(u)
(BT.2, BU.1), F(5u) = F(u), t = log_5 u, theta_k = int_0^1 F(5^t) e^(-2 pi i k t) dt = a_k - i b_k (AU's normalisation),
F = <F> + 2 sum_(k>=1)(a_k cos 2 pi k t + b_k sin 2 pi k t), G* = max F.

Theorem BW (checked here; exit 1 on any failure):
  (1) [certified, AU (3)'s routine unchanged, new parameters: mpmath.iv outward, DPS digits, K_PER cells,
      sandwich b = u = ORD/(ORD+1), depth 5^-E_DEPTH, node cap 10^CAP_EXP] a_k, b_k for k <= KMAX and <F> to the radius
      printed; theta_k != 0 for every k <= KMAX; |theta_k| <= B_k (BV.3) for every such k; the k <= 4 rows reproduce AU's
      22-digit table and the k = 5..11 rows agree with BV.4's 12-digit table within BV.4's radius 8.2e-63; |theta_k| e^(2 pi^2 k/ln5) in [0.60, 0.93].
  (2) [derived] c_R(gamma) = (2 sin(pi gamma)/pi) ln5 <Theta> (AU (1)) to the same relative accuracy, inside AU's
      34-digit enclosure.
  (3) [derived from BV.3 + (1), python-flint arb] BV.5's argument with the KMAX-mode sum: exactly one max and one min per
      period (re-proved), and t_max, t_min, G*, min F enclosed far inside BV.5's balls.
Assumptions used as stated (not re-run): AQ (3), BN.1, BS.1, BT.2, BU.1, BU.3 (G* <= 0.6009615785983), BV.3 (decay bound);
AU (3)'s certified routine (AQ (4) node enclosures, AS (2) sandwich for any b, u, AT.1 orbit lemma, AU's psi and endpoint
corrections), imported unchanged; AU's and BV's tables, BV.5's balls for consistency only.
Floats: none enter a check.  Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

from flint import arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasket_fractional_profile_analytic as bv  # noqa: E402

FAILS: list[str] = []
KMAX = 20
K_PER, ORD, DPS, CAP_EXP, E_DEPTH = 1024, 32, 300, 240, 260
BV_TABLE = {   # BV.4's printed 12-digit a_k, b_k (GASKET.md), k = 5..11
    5: ("1.43682897682e-27", "-9.05475087838e-28"), 6: ("-3.25583867317e-33", "-7.98359194922e-33"),
    7: ("-4.31465378795e-38", "-2.43242062331e-39"), 8: ("-5.50475216207e-44", "1.55340459757e-43"),
    9: ("8.38526711281e-49", "3.79766118274e-49"), 10: ("4.12559211982e-54", "-2.89164189910e-54"),
    11: ("-3.03451509014e-60", "-1.78263934544e-59")}
BV_BALLS = {   # BV.5's printed centres and radii
    "tmax": ("0.75175596179031755877802933273820825656220192611655", "9.1e-56"),
    "tmin": ("0.25175612419221392884930287700067310106574561482", "3.5e-55"),
    "Gmax": ("0.60096157859710108904205618413461692677004052591725819702", "4.5e-61"),
    "Fmin": ("0.60094692379426561210422618817212924841162969178767769700", "6.7e-61")}
BV_RADIUS = "8.21e-63"     # BV.4's common radius of a_k, b_k
AU_MEAN = ("0.6009542511687736078237051500999942", "0.6009542511687736078237051500999943")
AU_CRG = ("0.5171658346305052439259804856731271", "0.5171658346305052439259804856731273")


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def configure():
    import gasket_fractional_log_mean as au
    bv.KMAX, bv.K_PER, bv.ORD, bv.DPS, bv.CAP_EXP = KMAX, K_PER, ORD, DPS, CAP_EXP
    au.E_DEPTH = E_DEPTH


def main() -> int:
    t0 = time.time()
    print(f"Theorem BW: theta_k != 0 for k <= {KMAX}; <F>, c_R(gamma) to ~1e-113, G*, min F to ~1e-109, t_max, t_min to ~1e-103. Not thrust.")
    ctx.prec = 1200
    configure()
    G = arb(bv.G_STAR_HI)
    t1 = time.time()
    coef = bv.au_modes()
    print(f"   AU's routine ({K_PER} cells, b = u = {ORD}/{ORD + 1}, {DPS} digits, depth 5^-{E_DEPTH}, cap 1e{CAP_EXP}): {time.time() - t1:.0f}s")
    check(coef["nodes_ok"], "AQ (4) node enclosures of r_inf positive with relative width < 1e-100")
    m = coef["mean"]
    rad = max(float(coef[c][k].rad()) for c in ("a", "b") for k in coef["a"])
    print(f"   <F> = {m.str(95, radius=False)}  (radius {float(m.rad()):.1e});  common radius of a_k, b_k {rad:.2e}")
    check(bool(m.rad() < arb("1e-112")) and rad < 1e-112, "BW.1: <F> and every a_k, b_k (k <= KMAX) enclosed to radius < 1e-112")
    check(bool(m > arb(AU_MEAN[0])) and bool(m < arb(AU_MEAN[1])), "consistency: <F> inside AU's 34-digit enclosure of <Theta>")
    x = (-2 * arb.pi() ** 2 / arb(5).log()).exp()
    nz, below, tab, norm_ok = True, True, True, True
    rows = []
    for k in range(1, KMAX + 1):
        a, b = coef["a"][k], coef["b"][k]
        nlo = (abs(a).lower() ** 2 + abs(b).lower() ** 2).sqrt() if (not a.contains(0) or not b.contains(0)) else arb(0)
        nhi = (abs(a).upper() ** 2 + abs(b).upper() ** 2).sqrt()
        Bk = bv.B_bound(k, G)
        nz = nz and bool(nlo > 0)
        below = below and bool(nhi < Bk)
        cl, ch = nlo / x ** k, nhi / x ** k
        norm_ok = norm_ok and bool(cl > arb("0.60")) and bool(ch < arb("0.93"))
        if k <= 4:
            for c_, ref in ((a, bv.AU_TABLE_A[k]), (b, bv.AU_TABLE_B[k])):
                tab = tab and bool(abs(c_ - arb(ref)) < abs(arb(ref)) * arb("1e-21") + arb("2e-35"))
        elif k in BV_TABLE:      # BV prints 12 significant digits; its enclosures had radius 8.2e-63
            for c_, ref in zip((a, b), BV_TABLE[k]):
                tab = tab and bool(abs(c_ - arb(ref)) < abs(arb(ref)) * arb("5.01e-12") + arb(BV_RADIUS))
        rows.append((k, a, b, nlo, Bk, cl, ch))
        print(f"   k={k:2d}: a_k = {a.mid().str(12, radius=False):>20s}  b_k = {b.mid().str(12, radius=False):>20s}  "
              f"|theta_k| >= {float(nlo):.6e}  B_k = {float(Bk.upper()):.4e}  |theta_k| e^(2pi^2 k/ln5) in "
              f"[{float(cl.lower()):.4f}, {float(ch.upper()):.4f}]")
    check(nz, f"BW.1: theta_k != 0 certified for every k = 1..{KMAX} (BV.4: k <= 11), so AU (2)'s Mellin transform has genuine "
              f"poles at 1-gamma +- i omega k for every |k| <= {KMAX}, and F is not a trigonometric polynomial of degree < {KMAX}")
    check(below, f"BV.3 vs BW.1: |theta_k| <= B_k for every k = 1..{KMAX}")
    check(tab, "consistency: AU's 22-digit a_k, b_k (k <= 4) and BV.4's 12-digit a_k, b_k (k = 5..11) agree with the new enclosures")
    check(norm_ok, f"BW.1: |theta_k| e^(2 pi^2 k/ln5) in [0.60, 0.93] for every k = 1..{KMAX} (evidence for Conjecture BV.6, not a proof)")
    # (2) c_R(gamma)
    gam = arb(3).log() / arb(5).log()
    crg = 2 * (arb.pi() * gam).sin() / arb.pi() * arb(5).log() * m
    print(f"   c_R(gamma) = {crg.str(95, radius=False)}  (radius {float(crg.rad()):.1e})")
    check(bool(crg > arb(AU_CRG[0])) and bool(crg < arb(AU_CRG[1])) and bool(crg.rad() < arb("1e-112")),
          "BW.2: c_R(gamma) = (2 sin(pi gamma)/pi) ln5 <Theta> enclosed to < 1e-112, inside AU's 34-digit enclosure")
    # (3) extremes with the KMAX-mode sum
    e0, e1, e2 = (bv.tail_bound(KMAX + 1, j, G) for j in (0, 1, 2))
    print(f"   tails (k >= {KMAX + 1}): |F - S| <= {float(e0.upper()):.2e}, |F' - S'| <= {float(e1.upper()):.2e}, "
          f"|F'' - S''| <= {float(e2.upper()):.2e}")
    all_ok, zeros, locs = bv.extremes(coef, e1.upper(), e2.upper())
    check(all_ok and len(zeros) == 2 and sorted(s for _, s in zeros) == [-1, 1],
          f"BW.3: on all {bv.NCELL} cells F' != 0 or F'' != 0 and F' changes sign exactly twice: one max and one min per period")
    if len(locs) != 2:
        print("   extremes not located; stopping")
        return 1
    (tmax, _, _), = [l for l in locs if l[1] < 0]
    (tmin, _, _), = [l for l in locs if l[1] > 0]
    Fmax = bv.S_der(tmax, 0, coef) + arb(0, e0.upper())
    Fmin = bv.S_der(tmin, 0, coef) + arb(0, e0.upper())
    for name, v in (("t_max", tmax), ("t_min", tmin), ("G* = max F", Fmax), ("min F", Fmin)):
        print(f"   {name} = {v.str(95, radius=False)}  (radius {float(v.rad()):.1e})")
    check(bool(tmax.rad() < arb("1e-103")) and bool(tmin.rad() < arb("1e-103")) and bool(Fmax.rad() < arb("1e-108"))
          and bool(Fmin.rad() < arb("1e-108")), "BW.3: t_max, t_min enclosed to < 1e-103 and G*, min F to < 1e-108")
    inside = True
    for v, key in ((tmax, "tmax"), (tmin, "tmin"), (Fmax, "Gmax"), (Fmin, "Fmin")):
        c_, r_ = BV_BALLS[key]      # BV prints truncated centres: allow the radius plus one unit in the last printed digit
        ulp = arb(10) ** (-(len(c_) - 2))
        inside = inside and bool(abs(v - arb(c_)) < arb(r_) + ulp)
    check(inside, "consistency: t_max, t_min, G*, min F agree with BV.5's printed digits (radius plus last-digit truncation)")
    gap = tmax - tmin - arb(1) / 2
    print(f"   t_max - t_min - 1/2 = {gap.str(40, radius=False)};  max F - min F = {(Fmax - Fmin).str(60, radius=False)}")
    check(bool(gap < 0), "consistency: Kept Failure BV.5's witness t_max - t_min - 1/2 < 0 reproduced")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BW checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
