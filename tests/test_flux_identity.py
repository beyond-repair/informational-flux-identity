"""Sweep-168: lock the rectangle identity without rewriting witness.json.

Written as unittest.TestCase so the CI command
`python -m unittest tests/test_flux_identity.py` collects these checks
(plain test functions are invisible to unittest, which then exits 5 with
"NO TESTS RAN"). pytest collects the same cases.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import flux_identity as fi


class RectangleIdentityTests(unittest.TestCase):
    def test_divergence_free_witness_matches_readme(self):
        p = fi.potential()
        fx, fy = fi.fluxes(p)
        ledger = fi.ledger(fx, fy)
        self.assertEqual(ledger["div_max_abs"], 0)
        self.assertEqual(ledger["div_sum"], 0)
        self.assertEqual(ledger["net"], 0)
        self.assertEqual(ledger["signed"], {"left": 0, "right": 11, "bottom": -2, "top": -9})
        self.assertEqual(ledger["absolute"], {"left": 0, "right": 349, "bottom": 2, "top": 15})
        self.assertEqual(ledger["absolute_total"], 366)

    def test_theorem_a_on_random_integer_flux(self):
        rng = np.random.default_rng(168)
        fx = rng.integers(-5, 6, size=(8, 7), dtype=np.int64)
        fy = rng.integers(-5, 6, size=(7, 8), dtype=np.int64)
        div = fi.divergence(fx, fy)
        signed = int(fx[-1, :].sum() - fx[0, :].sum() + fy[:, -1].sum() - fy[:, 0].sum())
        self.assertEqual(int(div.sum()), signed)

    def test_boundary_source_matches_added_edge(self):
        p = fi.potential()
        fx, fy = fi.fluxes(p)
        fx = fx.copy()
        fx[-1, 3] = fx[-1, 3] + 17
        sourced = fi.ledger(fx, fy)
        self.assertEqual(sourced["div_sum"], 17)
        self.assertEqual(sourced["net"], 17)



class CornerMobiusFloatTests(unittest.TestCase):
    def test_theorem_eb_mobius_floats(self):
        """Theorem EB (OBSERVED part, numpy only): u_k = 1/(1 - delta_N^(k)) is linear in rho^-k up to O(5^-0.65k), so
        rho^k (u_k - u_(k+1)) / ((1 - 1/rho) u_oo) -> nu_0.8 (AZ (3): -0.72744538...) for N = 1 and N = 2 alike."""
        import gasket_fractional_corner_limit_mobius as eb

        rows = eb.nu_from_levels(0.8, (1, 2), (16, 22), 120)
        v1, v2 = rows[1][-1][1], rows[2][-1][1]
        self.assertLess(abs(v1 - v2), 1e-8)                      # N-independent
        self.assertLess(abs(v1 - (-0.7274453844421553)), 2e-6)  # AZ's nu (float quadrature bias ~5e-7)
        self.assertLess(abs(rows[1][-1][1] - v1), abs(rows[1][0][1] - v1) + 1e-15)
        self.assertLess(abs(rows[1][0][1] - v1), 2e-6)

class CornerRemainderFloatTests(unittest.TestCase):
    def test_theorem_ec_remainder_floats(self):
        """Theorem EC (OBSERVED part + exponent bookkeeping, numpy only): w_k = (u_(k+1) - u_(k+2)) - (u_k - u_(k+1))/rho kills
        EB.1's Moebius part, and w_(k+1)/w_k -> 1/3 (R_n ~ 3^-n); EC's proved threshold for the third term is 4 gamma/3 = 0.910142."""
        import math
        import gasket_fractional_corner_remainder_sharp as ec

        g = math.log(3) / math.log(5)
        self.assertAlmostEqual(4 * g / 3, 0.9101415926, places=9)
        self.assertLess(5 ** (3 * 0.9) / 27, 3)                   # rho(0.9)^3 < 3: third term now proved at 0.9
        self.assertGreater(5 ** (3 * 0.95) / 27, 3)               # rho(0.95)^3 > 3: not claimed
        ce = 2 * 0.6031 * 5.2**g / 6
        K1 = 0.2 ** (1 - g) * 16.6 / 30 + 2 * ce / (1 - 0.2**g * ce)
        self.assertLess(abs(K1 - 1.893324), 1e-5)                 # EC.1's constant (floats)
        for al in (0.8, 0.95):
            rows = ec.w_ratios(al, 1, (14, 16))
            self.assertLess(abs(rows[-1][1] - 1 / 3), 3e-3)
            self.assertLess(abs(rows[-1][1] - 1 / 3), abs(5**al / 9 - 1 / 3))  # not rho/3: the Moebius part is gone


class CornerRemainderLowerFloatTests(unittest.TestCase):
    def test_theorem_ed_remainder_lower_floats(self):
        """Theorem ED (numpy only): ED.2's constant ell_alpha (lower Riemann sum, floats) is positive and increasing in alpha,
        U(x1) = 2, the J-term thresholds rho^J < 3 <=> alpha < (J+1) gamma/J, and the floats of ED.5 sit above ED.2's bounds
        (Mt-increments (3/rho)^k (Mt_(k+1) - Mt_k) ~ 0.295 >= ell_0.8, (3/rho)^n eps_n ~ 0.74 >= L_0.8 > 0)."""
        import math
        import gasket_fractional_corner_remainder_lower as ed

        g = math.log(3) / math.log(5)
        x1 = 0.4 * math.exp(-0.1)
        self.assertAlmostEqual(5 * x1 * math.exp(5 * x1 * math.exp(0.1) / 20), 2.0, places=12)
        for J, thr in ((3, 4 * g / 3), (4, 5 * g / 4)):
            self.assertLess((5 ** (thr - 1e-9) / 3) ** J, 3)
            self.assertGreater((5 ** (thr + 1e-9) / 3) ** J, 3)
        lg, l8, l95 = ed.ell_float(g), ed.ell_float(0.8), ed.ell_float(0.95)
        self.assertGreater(lg, 0.0016681)
        self.assertLess(lg, l8)
        self.assertLess(l8, l95)
        L8 = 3 * l8 / (2 * (1 - 5**0.8 / 9))
        self.assertGreater(L8, 0.005532)
        rows, Cp = ed.float_increments(0.8, (4, 6))
        self.assertLess(abs(Cp + 4.658677573), 2e-5)            # C' = -Mt_oo rho/(rho-1) (ED.1) vs AZ (3)
        for n, d, e in rows:
            self.assertGreater(d, l8)
            self.assertGreater(e, L8)
            self.assertLess(abs(d - 0.297), 0.01)
            self.assertLess(abs(e - 0.746), 0.01)


class CornerRemainderLimitFloatTests(unittest.TestCase):
    def test_theorem_ee_remainder_limit_floats(self):
        """Theorem EE (numpy only): the limit formulas d* = (1/3) int x^(1-a) D (G + Gh), e* = 3 d*/(2(1-rho/3)) and kappa_N
        (floats, trapezoid in log x) sit inside EE.4's arb brackets at 0.8; the direct quadrature of (3/rho)^k (Mt_(k+1) - Mt_k)
        approaches d*; u_0 - u_1 = (s+6)/(s+5); EB.1 + EE.3 reproduce EC.5's observed w_(k+1)/w_k."""
        import math
        import gasket_fractional_corner_remainder_limit as ee

        s = np.array([0.0, 0.3, 2.0, 50.0])
        beta = (s + 6) / ((s + 2) * (s + 5))
        u1 = beta * (3 + s * (s + 5)) - 2 * s / (2 + s)
        self.assertTrue(np.allclose(3 + s - u1, (s + 6) / (s + 5)))
        lf = ee.limits_float(0.8, Ns=(1, 2, 3), n_per=40)
        self.assertTrue(0.2950346 < lf["d"] < 0.2951296)
        self.assertTrue(0.7408653 < lf["e"] < 0.7411039)
        self.assertAlmostEqual(lf["e"], 3 * lf["d"] / (2 * (1 - 5**0.8 / 9)), places=14)
        for N, (lo, hi) in {1: (0.101641, 0.101739), 2: (0.91313, 0.914013), 3: (2.83752, 2.84026)}.items():
            self.assertTrue(lo < lf["kappa"][N] < hi)
        c8, c12 = ee.direct_increment_float(0.8, 8), ee.direct_increment_float(0.8, 12)
        self.assertLess(abs(c12 - lf["d"]), abs(c8 - lf["d"]))
        self.assertLess(abs(c12 - lf["d"]), 3e-4)
        self.assertGreater(c12, 0.005532 * 2 * (1 - 5**0.8 / 9) / 3)          # above ED.2's ell_0.8 bound
        al, nu = 0.8, -0.7274453844421553
        Dinf = 2.3964097039792205 / (2 * math.sin(math.pi * al) / math.pi)
        mod = ee.ec5_model(al, nu, Dinf, 0.03131 / -nu, lf["e"], lf["kappa"][1], (14, 16, 18, 20))
        for (_, r, _), o in zip(mod, (0.33637, 0.33541, 0.33475, 0.33426)):
            self.assertLess(abs(r - o), 3e-4)
        self.assertLess(abs(mod[-1][2] + 5.04), 0.05)


class CornerCriticalSignFloatTests(unittest.TestCase):
    def test_theorem_ef_critical_sign_floats(self):
        """Theorem EF (numpy only): EF.1's formulas D'_oo = int x^-a G, nu = -rho Mt_oo/((rho-1) D'_oo), 1 - delta_N = h_oo/D'_oo
        reproduce AZ (3)'s A_0.8, C_0.8, nu_0.8 and EA.4's C_delta(0.8, 1) in floats; at alpha = 4gamma/3 the floats of nu and
        S_N = -nu^3(1-delta_N) - varrho_N sit inside EF.3's arb brackets and |nu|^3 < e*/D'_oo (S_N < 0 for every N)."""
        import math
        import gasket_fractional_corner_critical_sign as ef

        q = np.array([0.0, 0.01, 1.0, 30.0, 1e4])
        self.assertTrue(np.all((1 + q / 5) * (1 + q / 2) >= (1 + 2 * q / 5) * (1 + q / 6)))     # omega^0 >= 1
        t = np.array([1e-9, 0.3, 7.0, 1e6])
        qq = 2 * t / (5 + np.sqrt(25 + 4 * t))
        self.assertTrue(np.allclose(qq * (qq + 5), t, rtol=1e-12))
        f8 = ef.critical_floats(0.8, Ns=(1,), n_per=40)
        kap = 2 * math.sin(math.pi * 0.8) / math.pi
        self.assertLess(abs(f8["nu"] + 0.7274453844421553), 1e-7)           # AZ (3)
        self.assertLess(abs(kap * f8["Dp"] - 2.3964097039792205), 1e-7)
        self.assertLess(abs(kap * f8["Cp"] + 1.7432571783920757), 1e-7)
        self.assertTrue(0.031293 < -f8["nu"] * f8["om"][1] < 0.031327)      # EA.4: C_delta = -nu (1 - delta_1)
        fc = ef.critical_floats(4 * ef.GAMMA / 3, n_per=40)
        self.assertTrue(-0.396698 < fc["nu"] < -0.3965332)
        self.assertGreater(fc["e"] / fc["Dp"], -fc["nu"] ** 3)
        for N, (lo, hi) in {1: (-0.0233127, -0.0232863), 2: (-0.174681, -0.174497), 3: (-0.454529, -0.454077),
                            6: (-4.69362, -4.68928), 10: (-89.0757, -88.9943)}.items():
            self.assertTrue(lo < fc["S"][N] < hi)


class CornerExplicitThresholdFloatTests(unittest.TestCase):
    def test_theorem_eg_explicit_threshold_floats(self):
        """Theorem EG (numpy only): the every-n lower bound 1/V_k0 <= R_k (k >= k0) holds against the exact recursion, c_lo(k0)
        sits below the direct quadrature of c_k at 0.8 and increases to EE's d*, and at alpha = 4gamma/3 the float L(2) clears
        |nu|^3 D'_oo (EF.3's upper end 0.06243 * 6.005766) with the margin the arb run certifies."""
        import math
        import gasket_fractional_corner_explicit_threshold as eg
        import gasket_fractional_corner_remainder_limit as ee

        x = np.array([1e-3, 0.3, 1.0, 4.0, 20.0, 150.0])
        for k0 in (1, 2, 4):
            V = eg.vplus_float(x, k0)[0]
            for k in (k0, k0 + 1, k0 + 5, k0 + 15):
                self.assertTrue(np.all(1 / V <= eg.rk_float(x, k) * (1 + 1e-12)))
        al = 4 * eg.GAMMA / 3
        L = {k0: eg.lower_float(al, k0, n=60)["L"] for k0 in (1, 2, 4, 8)}
        self.assertTrue(L[1] < L[2] < L[4] < L[8] < 0.7399065)             # below EF.3's e*(4gamma/3)
        self.assertGreater(L[2], 0.06243 * 6.005766 + 0.08)                  # L(2) > |nu|^3 D'_oo with margin
        self.assertLess(L[1], 0.06235 * 6.005548)                            # k0 = 1 alone does not suffice
        self.assertGreater(L[8], 0.739)
        for k in (2, 4):
            self.assertLess(eg.lower_float(0.8, k, n=60)["c"], ee.direct_increment_float(0.8, k, h=4e-3))


class CornerUniformThresholdFloatTests(unittest.TestCase):
    def test_theorem_eh_uniform_threshold_floats(self):
        """Theorem EH (numpy only): on a grid of [4gamma/3, 1) EF/EG's float routines give (3/rho^3)^2 |nu|^3 D'_oo < L(2) and an
        n = 1 bound < 0 (L(1) alone from about 0.92 on, the eta route below), |nu| decreasing to 0 as alpha -> 1, inside the
        arb cell brackets the script certifies."""
        import gasket_fractional_corner_uniform_threshold as eh

        al = [4 * eh.GAMMA / 3 + 1e-9, 0.92, 0.94, 0.96, 0.98, 0.995]
        rows = eh.survey_float(al)
        for r in rows:
            self.assertLess(r["f"] ** 2 * r["T"], r["L2"] - 0.1)
            self.assertLess(r["S11"], 0)
            self.assertLess(r["nu"], 0)
            self.assertLess(-r["nu"], r["rho"])
        self.assertGreater(rows[0]["f"] * rows[0]["T"], rows[0]["L1"])      # at 4gamma/3, n = 1 needs the eta term
        for r in rows[2:]:
            self.assertLess(r["f"] * r["T"], r["L1"])                       # from 0.94 on L(1) alone suffices
        self.assertTrue(all(rows[i + 1]["nu"] > rows[i]["nu"] for i in range(len(rows) - 1)))
        self.assertTrue(0.36830 < -rows[0]["nu"] < 0.40822)                 # EH.3's first-cell |nu| bracket
        self.assertTrue(0.02631 < -eh.survey_float([0.9925])[0]["nu"] < 0.05411)   # EH.3's cell [0.99, 0.995]
        # EH.5 below 4gamma/3: (3/rho^3)^n |nu|^3 D'_oo < L(8) on the certified window n = 3..8 at 0.9, and it fails for large n
        r9 = eh.survey_float([0.9])[0]
        L8 = eh._mods()[1].lower_float(0.9, 8, n=60)["L"]
        g = 3 / r9["rho"] ** 3
        self.assertGreater(g, 1)
        self.assertTrue(all(g**n * r9["T"] < L8 for n in range(3, 9)))
        self.assertGreater(g**12 * r9["T"], 0.7380)                          # beyond EE.4's e*(0.9) upper end: no window there


class CornerOnsetLevelFloatTests(unittest.TestCase):
    def test_theorem_ei_onset_level_floats(self):
        """Theorem EI (numpy only): EC.2's B_alpha and DZ.1's 3^k etabar_k in floats; EI.1's criterion Phi_N(k) with EF's float
        D'_oo, nu, 1-delta_N is increasing in k and first positive at the levels n1 the script certifies with arb (floats may only be
        smaller); below 4gamma/3 the onset grows as g = 3/rho^3 -> 1, and every n1 exceeds EH.4's window end."""
        import gasket_fractional_corner_onset_level as ei

        self.assertTrue(1937.0 < ei.B_float(0.8) < 1938.0)                  # EC.2: B_0.8 <= 1.94e3 (1.9375e3)
        self.assertTrue(1028.0 < ei.B_float(0.9) < 1029.0)                  # EC.2: B_0.9 <= 1.03e3
        self.assertGreater(ei.H_float(0.9, 1, 6), ei.H_float(0.9, 1, 7))   # H_k decreasing in k
        cert = {0.8: [13, 13, 13, 14, 18], 0.9: [159, 158, 159, 162, 191]}   # EI.2 (arb, --full; default: 0.9 one level later)
        prev = None
        for al in (0.8, 0.9):
            s = ei.onset_float(al)
            self.assertGreater(s["g"], 1.0)
            n1 = [s["n1"][N] for N in ei.NS]
            for a, b in zip(n1, cert[al]):
                self.assertLessEqual(a, b)
                self.assertGreaterEqual(a, b - 1)
            if prev is not None:
                self.assertTrue(all(a > b for a, b in zip(n1, prev)))
            prev = n1
        self.assertGreater(min(n1[:4]), 8)                                   # beyond EH.4's window end 8 at 0.9
        self.assertTrue(all(s["nas"][N] > 8 for N in (1, 2, 3, 6)))         # limit-substituted crossing after the window (OBSERVED)
        self.assertTrue(all(s["nas"][N] < s["n1"][N] for N in ei.NS))



class CornerTwoSidedRateFloatTests(unittest.TestCase):
    def test_theorem_ej_two_sided_rate_floats(self):
        """Theorem EJ (numpy only): EJ.1's every-level profile ceiling R^hi(x) = min(1/3, G x^(gamma-1), 1/V^lo) lies above the exact
        R_k(x) (affine recursion) for k = k0..k0+15; the float upper integrals U(14), kappa_hi(1, 14) sit just above EE's e*, kappa_1
        at 0.9, and EJ.2's criterion with them first holds at n = 33 (the level the script certifies with arb on --full)."""
        import numpy as np
        import gasket_fractional_corner_two_sided_rate as ej

        x = np.array([1e-3, 0.3, 1.0, 4.0, 20.0, 150.0, 3000.0])
        for k0 in (2, 6):
            R = ej.rhi_float(x, k0)[0]
            for k in range(k0, k0 + 16):
                self.assertTrue(np.all(ej.rk_float(x, k) <= R * (1 + 1e-12)))
        uf = ej.upper_float(0.9, 14, Ns=(1,), n=40)
        n, c = ej.crossing_float(0.9, uf["U"], uf["kappa"][1])
        self.assertTrue(c["e"] < uf["U"] < 1.03 * c["e"])
        self.assertTrue(c["kappa"][1] < uf["kappa"][1] < 1.03 * c["kappa"][1])
        self.assertEqual(n, 33)

class CornerExact0905FloatTests(unittest.TestCase):
    def test_theorem_ek_exact_0905_floats(self):
        """Theorem EK (numpy only): at alpha = 0.905 the float version of EJ.2's criterion with EJ.1's float upper integrals (k0 = 14)
        first holds at n = 68, 84, 104 for N = 1, 2, 3 -- the onset levels EK.1 certifies with arb on --full -- and the float U(14),
        kappa_hi(2), kappa_hi(3) lie above EK's certified lower bounds L(12) >= 0.73421, kappa_lo(2, 14) >= 0.618019,
        kappa_lo(3, 14) >= 1.70926 and U below EK's certified U(14) <= 0.740363."""
        import gasket_fractional_corner_exact_0905 as ek
        import gasket_fractional_corner_two_sided_rate as ej

        self.assertEqual(ek.A_S, "0.905")
        self.assertEqual(ek.EXPECT, {1: (67, 68), 2: (83, 84), 3: (102, 104)})
        uf = ej.upper_float(0.905, 14, n=40)
        self.assertTrue(0.73421 < uf["U"] < 0.740363)
        self.assertTrue(uf["kappa"][2] > 0.618019 and uf["kappa"][3] > 1.70926)
        for N in ej.NS:
            n, _ = ej.crossing_float(0.905, uf["U"], uf["kappa"][N], N=N)
            self.assertEqual(n, ek.EXPECT[N][1])
            self.assertTrue(ek.EJ_BRACKET[N][0] <= ek.EXPECT[N][0] <= ek.EJ_BRACKET[N][1])


class CornerExact0905N3FloatTests(unittest.TestCase):
    def test_theorem_el_exact_0905_n3_floats(self):
        """Theorem EL (numpy only): EL.1's left-node cell bound in floats -- on cells [a, a 5^(1/1600)] the exact R_k (affine recursion,
        k = 26..41) never exceeds EJ.1's ceiling R^hi(a) at the left node, and R_k is nonincreasing across the cell; EE's float limits
        e*, kappa_N put the first positive level at 68, 84, 103 for N = 1, 2, 3 at 0.905 (EL.2 certifies 103 with arb) and at 2586
        at 0.91 (EL.3 certifies every n >= 2587); EL's certified U(26) <= 0.739060, kappa_hi(3, 26) <= 1.722784 sit above the float
        limits and still give onset 103 in floats."""
        import numpy as np
        import gasket_fractional_corner_exact_0905_n3 as el
        import gasket_fractional_corner_two_sided_rate as ej

        self.assertEqual(el.EXPECT, {1: 68, 2: 84, 3: 103})
        self.assertEqual((el.EXPECT_091, el.EJ_091), (2587, (2564, 2616)))
        self.assertEqual(el.EK_LAST[3] + 1, el.EXPECT[3])          # with EK.1: one sign change 102 -> 103
        for t in (1.0, 2.5, 4.0, 5.0):
            a = 5.0**t
            x = a * 5.0 ** (np.linspace(0.0, 1.0, 6) / 1600)
            Ra = ej.rhi_float(np.array([a]), 26)[0][0]
            for k in range(26, 42):
                r = ej.rk_float(x, k)
                self.assertTrue(np.all(r <= Ra * (1 + 1e-12)))
                self.assertTrue(np.all(np.diff(r) <= 1e-15 * r[0]))
        ef, _ = ej._mods()
        c = ef.critical_floats(0.905, Ns=(1, 2, 3))
        for N in ej.NS:
            self.assertEqual(ej.crossing_float(0.905, c["e"], c["kappa"][N], N=N)[0], el.EXPECT[N])
        self.assertTrue(c["e"] < 0.739060 and c["kappa"][3] < 1.722784)
        self.assertEqual(ej.crossing_float(0.905, 0.739060, 1.722784, N=3)[0], 103)
        c91 = ef.critical_floats(0.91, Ns=(1,))
        n91 = ej.crossing_float(0.91, c91["e"], c91["kappa"][1], N=1)[0]
        self.assertTrue(el.EJ_091[0] < n91 < el.EXPECT_091)


class CornerExact091FloatTests(unittest.TestCase):
    def test_theorem_em_exact_091_floats(self):
        """Theorem EM (numpy only): EM.1's right-node lower cell bound in floats -- on cells [b 5^(-1/1600), b] the exact R_k (affine
        recursion, k = 14..29) never falls below EG.1's 1/V_14(b) at the right node, and R_k is nonincreasing across the cell; EE's
        float limits put the first positive level at 2586 at 0.91 and 15 at 0.89 (N = 1), as EM.2 / EM.3 certify with arb; EM's
        certified U(26) <= 0.73986384, kappa_hi(1, 26) <= 0.07466495 (0.91) and U(15) <= 0.737190, kappa_hi(1) <= 0.079005 (0.89)
        sit above the float limits and still give onsets 2586, 15 in floats."""
        import numpy as np
        import gasket_fractional_corner_exact_091 as em
        import gasket_fractional_corner_two_sided_rate as ej

        self.assertEqual((em.EXPECT91, em.EXPECT89), ((2585, 2586), (14, 15)))
        self.assertEqual((em.EJ4_LAST91, em.EL3_ONSET91), (2564, 2587))
        for t in (-3.0, 0.0, 2.5, 4.0, 5.0):
            b = 5.0**t
            x = b * 5.0 ** (-np.arange(5, -1, -1) / (5 * 1600))
            Rb = em.vlo_right_float(np.array([b]), 14)[0]
            for k in range(14, 30):
                r = ej.rk_float(x, k)
                self.assertTrue(np.all(r >= Rb * (1 - 1e-12)))
                self.assertTrue(np.all(np.diff(r) <= 1e-15 * r[0]))
        ef, _ = ej._mods()
        c91 = ef.critical_floats(0.91, Ns=(1,))
        self.assertEqual(ej.crossing_float(0.91, c91["e"], c91["kappa"][1], N=1)[0], em.EXPECT91[1])
        self.assertTrue(0.73971486 < c91["e"] < 0.73986384 and 0.07464970 < c91["kappa"][1] < 0.07466495)
        self.assertEqual(ej.crossing_float(0.91, 0.73986384, 0.07466495, N=1)[0], em.EXPECT91[1])
        c89 = ef.critical_floats(0.89, Ns=(1,))
        self.assertEqual(ej.crossing_float(0.89, c89["e"], c89["kappa"][1], N=1)[0], em.EXPECT89[1])
        self.assertTrue(c89["e"] < 0.737190 and c89["kappa"][1] < 0.079005)
        self.assertEqual(ej.crossing_float(0.89, 0.737190, 0.079005, N=1)[0], em.EXPECT89[1])


class CornerExact088FloatTests(unittest.TestCase):
    def test_theorem_en_exact_088_floats(self):
        """Theorem EN (numpy only): EN.1's two-sided D'_n in floats at alpha = 0.88, n = k0 = 9.  The scaled recursion w_j = u_j/s_j
        reproduces EJ's R_k where both are finite; the trapezoid D'_9 = int x^-a R_9 lies inside EN.1's bracket rebuilt from the
        certified D'_oo, Mt_oo, L(9), U(9) of EN.2 (and well below D'_oo); EJ.2 (b)'s braces bound at n = 9 recomputed from the certified
        endpoints is > 0 with D'_9 <= D'_oo (EM.4: undecided) and < 0 with EN.1's D_9^hi (EN.2: Y_9 < Y_oo); the float crossing is 10."""
        import numpy as np
        import gasket_fractional_corner_exact_088 as en
        import gasket_fractional_corner_two_sided_rate as ej

        self.assertEqual((en.EXPECT, en.ONSET, en.EJ4_WINDOW, en.K0), ((9, 10), 10, (2, 8), 9))
        x = 5.0 ** np.linspace(-6.0, 4.0, 41)
        for k in (1, 5, 9, 14):
            a, b = en.rk_scaled_float(x, k), ej.rk_float(x, k)
            self.assertTrue(np.all(np.abs(a - b) <= 1e-9 * np.abs(b)))
        al, n = 0.88, 9
        D = (5.4676986, 5.4679626)          # EN.2: D'_oo (rounded outward)
        M = (0.74287912, 0.74306202)        # EN.2: Mt_oo
        L, U = 0.73446007, 0.74222855       # EN.2: L(9) (lower), U(9) (upper)
        nu_hi, om = 0.4993144, (0.068040179, 0.06804675)
        klo = 0.08095188
        d9lo, d9hi = en.dn_bounds_float(D[0], D[1], M[0], M[1], L, U, al, n)
        self.assertTrue(abs(d9lo - 5.31126720) < 2e-8 and abs(d9hi - 5.31157011) < 2e-8)
        d9 = en.dprime_float(al, n)
        self.assertTrue(d9lo < d9 < d9hi < D[0])
        rho = 5**al / 3
        g = 3 / rho**3

        def braces(Dh):
            X = nu_hi**3 * g**n - L / Dh
            return (om[0] * X if X < 0 else om[1] / (1 - nu_hi * rho ** (-n)) * X) - 1.5 * klo / Dh

        self.assertGreater(braces(D[1]), 5e-4)     # EJ.2 (b) with D'_9 <= D'_oo: +6.965e-4, undecided
        self.assertLess(braces(d9hi), -2e-4)       # EN.1: -2.344e-4 < 0
        ef, _ = ej._mods()
        c = ef.critical_floats(al, Ns=(1,))
        self.assertTrue(L < c["e"] < U and klo < c["kappa"][1] < 0.08322902)
        self.assertEqual(ej.crossing_float(al, c["e"], c["kappa"][1], N=1)[0], en.EXPECT[1])


class CornerLevelOneFloatTests(unittest.TestCase):
    def test_theorem_eo_level_one_floats(self):
        """Theorem EO (numpy only): level 1 in closed form.  r_1 from the level-1 gasket Laplacian (6 vertices) equals AQ (1)'s
        (s + 5)/(s^2 + 7s + 9) and EJ's R_1; EO.1's closed forms of D'_1 and 1 - delta_1^(1) equal the trapezoid of x^-a R_1 and
        x^-a R_1 eta_1; Y_1 - Y_oo = rho^2 (om_1 - om) - om |nu| (rho + |nu|) recomputed from EN.2's / EO.2's certified om, |nu| endpoints
        (worst corner) is < 0 at 0.88 and 0.89; float signs of Y_n - Y_oo for n = 1..12 match the certified patterns."""
        import numpy as np
        import gasket_fractional_corner_level_one as eo
        import gasket_fractional_corner_two_sided_rate as ej

        E = [(0, 3), (3, 1), (1, 4), (4, 2), (2, 5), (5, 0), (3, 4), (4, 5), (5, 3)]
        L = np.zeros((6, 6))
        for i, j in E:
            L[i, i] += 1
            L[j, j] += 1
            L[i, j] -= 1
            L[j, i] -= 1
        v = np.array([1.0, -1.0, 0, 0, 0, 0])
        for s in (1e-3, 0.3, 1.0, 4.0, 50.0):
            r = 0.5 * v @ np.linalg.solve(L + s * np.eye(6), v)
            self.assertAlmostEqual(r, (s + 5) / (s * s + 7 * s + 9), places=13)
            self.assertAlmostEqual(r, (5 / 3) * float(ej.rk_float(np.array([5 * s]), 1)[0]), places=12)
        for al in (0.7, 0.8, 0.85, 0.88, 0.89, 0.9, 0.905, 0.91):
            rho = 5**al / 3
            D, om1 = eo.level_quad_float(al, 1)
            self.assertLess(abs(D - eo.j1_float(al) / rho), 1e-8)
            self.assertLess(abs(om1 - eo.om1_float(al)), 1e-10)
        # EO.2 recomputed from the certified endpoints (outward-rounded brackets; gap decreasing in om and |nu|)
        for al, om, nab, om1, ub in ((0.88, (0.068040179, 0.06804675), (0.49916739, 0.4993144), 0.0927170037, -1.703396e-02),
                                     (0.89, (0.070945361, 0.070951642), (0.46638151, 0.46651559), 0.0933602177, -1.793204e-02)):
            self.assertAlmostEqual(eo.om1_float(al), om1, places=9)
            worst = eo.level1_gap_float(al, nab[0], om[0])
            self.assertLess(abs(worst - ub), 1e-7)
            self.assertLess(worst, -0.017)
        ef, _ = ej._mods()
        for al, last in ((0.88, 9), (0.89, 14)):
            c = ef.critical_floats(al, Ns=(1,))
            vals = [eo.gap_n_float(al, n, -c["nu"], c["om"][1]) for n in range(1, 13)]
            self.assertEqual([n for n, x in enumerate(vals, 1) if x < 0], list(range(1, min(last, 12) + 1)))
            self.assertLess(abs(vals[0] - eo.level1_gap_float(al, -c["nu"], c["om"][1])), 1e-8)


class CornerExactLevelsFloatTests(unittest.TestCase):
    def test_theorem_ep_exact_levels_floats(self):
        """Theorem EP (numpy only): levels n = 1..6 at alpha = 0.85.  The level-n gasket Laplacian (n = 1..3; 6, 15, 42 vertices) gives
        the same r_n as AQ (1)'s recursion, and its eigen-expansion r_n = sum c_i/(s + lambda_i) (c_i = (v.phi_i)^2/2 >= 0) gives
        J_n = pi/sin(pi a) sum c_i lambda_i^-a = rho^n x trapezoid D'_n; the trapezoid 1 - delta_1^(n) equals EP.2's exact om_n
        (n = 1..6); Y_n - Y_oo = rho^(2n) (om_n - om) - om |nu| (rho^n + |nu|) recomputed at the worst corner of EP.2's certified
        om, |nu| brackets (widened by 1e-8) is < 0 for n = 1, 2, 3 and > 0 for n = 4, 5, 6."""
        import numpy as np
        import gasket_fractional_corner_exact_levels as ep
        import gasket_fractional_corner_level_one as eo

        self.assertEqual((ep.EXPECT85, ep.ONSET85, ep.SIGN85), ((3, 4), 5, {1: -1, 2: -1, 3: -1, 4: 1, 5: 1, 6: 1}))
        al = 0.85
        rho = 5**al / 3
        for n in (1, 2, 3):
            L, v = ep.gasket_laplacian(n)
            self.assertEqual(len(v), (3 ** (n + 1) + 3) // 2)
            for s in (1e-3, 0.3, 1.0, 4.0, 50.0):
                self.assertLess(abs(ep.r_graph(n, s) - ep.r_recursion(n, s)), 1e-12)
            lam, phi = np.linalg.eigh(L)
            c = 0.5 * (phi.T @ v) ** 2
            keep = lam > 1e-9
            self.assertLess(float(np.sum(c[~keep])), 1e-20)            # v is orthogonal to the constants
            J = np.pi / np.sin(np.pi * al) * float(np.sum(c[keep] * lam[keep] ** (-al)))
            D, _ = eo.level_quad_float(al, n)
            self.assertLess(abs(J / rho**n - D), 1e-8)
        om, nab = (0.059038919 - 1e-8, 0.059046046 + 1e-8), (0.59110183 - 1e-8, 0.59128836 + 1e-8)   # EP.2's brackets, widened
        for n, om_n in ep.OM_EXACT85.items():
            self.assertLess(abs(eo.level_quad_float(al, n)[1] - om_n), 1e-8)
            if ep.SIGN85[n] < 0:
                self.assertLess(ep.gap_float(al, n, nab[0], om[0], om_n), -1.5e-3)
            else:
                self.assertGreater(ep.gap_float(al, n, nab[1], om[1], om_n), 4e-4)


if __name__ == "__main__":
    unittest.main()
