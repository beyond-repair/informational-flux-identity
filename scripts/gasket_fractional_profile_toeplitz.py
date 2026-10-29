#!/usr/bin/env python3
r"""Theorem BY: AU's certified mode routine, re-run with sandwich order 40/41 at depth 5^-340 (380 digits), certifies
theta_21 != 0 (BW: k <= 20) and the Toeplitz sections T_18, T_19, T_20 of BX's spectral measure m positive definite
(BX: T_17), so m has at least 20 support points per scale period and nu is not carried by <= 19 geometric pole series
(Kept Failure BY.3, extends Kept Failure BX.5 from 16 to 19).

Setting (GASKET.md, Theorems AQ, AU, BN, BS-BX): gamma = log3/log5, omega = 2 pi/ln5, F(5^t) = theta(5^t) 1-periodic,
theta_k = a_k - i b_k (AU's normalisation), m the positive 1-periodic measure of BX.2 with
theta_k = (pi/ln5) m_k/sin(pi(gamma + i omega k)), rho_k = |m_k|/m_0.

Theorem BY (checked here; exit 1 on any failure):
  (1) [certified, AU (3)'s routine unchanged, new parameters: mpmath.iv outward, 380 digits, 1024 cells, sandwich
      b = u = 40/41, depth 5^-340, node cap 1e310] <F>, a_k, b_k (k <= 21) enclosed to radius < 1e-132 (BW: 2e-113);
      theta_k != 0 for every 1 <= |k| <= 21; the k <= 4 rows reproduce AU's 22-digit table, k = 5..11 BV.4's table.
  (2) [certified + BX.2, BX.3 as stated] rho_k in [0.60, 0.92] for every k <= 21 (rho_21 = 0.70637...), and
      |theta_k| <= <F> sin(pi gamma)/|sin(pi(gamma + i omega k))| for every k <= 21.
  (3) [certified + BX.5's rank lemma] T_n = (m_(i-j))_(i,j<n) is positive definite for n = 18, 19, 20 (exact-midpoint
      LDL* of T_n - 2||E||_F I in arb, E the enclosure radii; Weyl).
  Kept Failure BY.3 (certified witness): "nu is carried by r <= 19 geometric series lam_j 5^Z", i.e.
      Q_inf(u) = sum_(j<=r) c_j sum_n 3^n/(lam_j 5^n + u), is false: r atoms give rank T_(r+1) <= r, but T_20 > 0.
Observation (no claim): T_21 is not decided at this precision.
Assumptions used as stated (not re-run): AQ (3), BN.1, BS.1, BT.2, BU.1, BV.2, AU (1), Assumption BX (Nevanlinna),
BX.2-BX.3 and BX.5's rank lemma; AU (3)'s certified routine (AQ (4), AS (2) for any b, u, AT.1, AU's psi and endpoint
corrections), imported unchanged with new parameters.  Floats: none enter a check.
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

from flint import acb, arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasket_fractional_profile_spectral as bx  # noqa: E402

bw = bx.bw
bv = bw.bv
FAILS: list[str] = []
KMAX = 21
K_PER, ORD, DPS, CAP_EXP, E_DEPTH = 1024, 40, 380, 310, 340
N_TOEP = (18, 19, 20)
RHO_BAND = ("0.60", "0.92")
RAD_MAX = "1e-132"


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
    print(f"Theorem BY: theta_k != 0 for k <= {KMAX}; T_18, T_19, T_20 of BX's spectral measure positive definite; "
          "Kept Failure BY.3: nu is not carried by <= 19 geometric pole series. Not thrust.")
    bx.FAILS.clear()
    bx.toeplitz_toy_checks()
    check(not bx.FAILS, "certificate sanity (BX): for an exact 3-atom measure T_3 is certified and T_4 is not")
    ctx.prec = 1500
    configure()
    t1 = time.time()
    coef = bv.au_modes()
    print(f"   AU's routine ({K_PER} cells, b = u = {ORD}/{ORD + 1}, {DPS} digits, depth 5^-{E_DEPTH}, cap 1e{CAP_EXP}): "
          f"{time.time() - t1:.0f}s")
    check(coef["nodes_ok"], "AQ (4) node enclosures of r_inf positive with relative width < 1e-100")
    mean = coef["mean"]
    rad = max(float(coef[c][k].rad()) for c in ("a", "b") for k in coef["a"])
    print(f"   <F> = {mean.str(60, radius=False)}  (radius {float(mean.rad()):.1e}); common radius of a_k, b_k {rad:.2e}")
    check(bool(mean.rad() < arb(RAD_MAX)) and rad < float(RAD_MAX),
          f"BY.1: <F> and every a_k, b_k (k <= {KMAX}) enclosed to radius < {RAD_MAX}")
    check(bool(mean > arb(bw.AU_MEAN[0])) and bool(mean < arb(bw.AU_MEAN[1])), "consistency: <F> inside AU's 34-digit enclosure")
    nz, tab = True, True
    for k in range(1, KMAX + 1):
        a, b = coef["a"][k], coef["b"][k]
        nz = nz and ((not a.contains(0)) or (not b.contains(0)))
        if k <= 4:
            for c_, ref in ((a, bv.AU_TABLE_A[k]), (b, bv.AU_TABLE_B[k])):
                tab = tab and bool(abs(c_ - arb(ref)) < abs(arb(ref)) * arb("1e-21") + arb("2e-35"))
        elif k in bw.BV_TABLE:
            for c_, ref in zip((a, b), bw.BV_TABLE[k]):
                tab = tab and bool(abs(c_ - arb(ref)) < abs(arb(ref)) * arb("5.01e-12") + arb(bw.BV_RADIUS))
    a21, b21 = coef["a"][KMAX], coef["b"][KMAX]
    lo21 = (abs(a21).lower() ** 2 + abs(b21).lower() ** 2).sqrt()
    print(f"   k=21: a_21 = {a21.mid().str(15, radius=False)}  b_21 = {b21.mid().str(15, radius=False)}  "
          f"|theta_21| >= {float(lo21):.6e}")
    check(nz, f"BY.1: theta_k != 0 certified for every k = 1..{KMAX} (BW: k <= 20)")
    check(tab, "consistency: AU's (k <= 4) and BV.4's (k = 5..11) published a_k, b_k agree with the new enclosures")
    pi, l5, g, w = bx.consts()
    m = {0: bx.m_from_theta(0, acb(mean))}
    for k in range(1, KMAX + 1):
        m[k] = bx.m_from_theta(k, acb(coef["a"][k], -coef["b"][k]))
    m0 = m[0].real
    band, bd = True, True
    for k in range(1, KMAX + 1):
        rho = abs(m[k]) / m0
        band = band and bool(rho > arb(RHO_BAND[0])) and bool(rho < arb(RHO_BAND[1]))
        bd = bd and bool(abs(acb(coef["a"][k], -coef["b"][k])) < bx.bound_k(k, mean))
        print(f"   k={k:2d}: rho_k = {rho.str(12, radius=False):>16s} (rad {float(rho.rad()):.1e})  "
              f"arg m_k/2pi = {float((m[k].arg() / (2 * pi)).mid()):+.10f}")
    check(band, f"BY.2: rho_k in [{RHO_BAND[0]}, {RHO_BAND[1]}] for every k = 1..{KMAX} (certified)")
    check(bd, f"BY.2: |theta_k| <= BX.3's bound for every k = 1..{KMAX}")
    for n in N_TOEP:
        ok, E, D = bx.toeplitz_certificate(m, n)
        dmin = min(D, key=lambda d: float(d.mid())) if D else arb(0)
        print(f"   T_{n}: ||E||_F <= {float(E.upper()):.2e}; smallest pivot of T - 2||E|| I = {float(dmin.mid()):.3e}")
        check(ok, f"BY.3: T_{n} = (m_(i-j)) positive definite (min eigenvalue > ||E||_F)")
    print(f"   Kept Failure BY.3: T_{N_TOEP[-1]} > 0, so m has >= {N_TOEP[-1]} support points per period and nu is not "
          f"carried by <= {N_TOEP[-1] - 1} geometric series lam_j 5^Z")
    ok21, E21, _ = bx.toeplitz_certificate(m, N_TOEP[-1] + 1)
    print(f"   observation (no claim): T_{N_TOEP[-1] + 1} certified: {ok21} (||E||_F <= {float(E21.upper()):.1e})")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BY checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
