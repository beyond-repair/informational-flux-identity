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


if __name__ == "__main__":
    unittest.main()
