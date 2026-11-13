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


if __name__ == "__main__":
    unittest.main()
