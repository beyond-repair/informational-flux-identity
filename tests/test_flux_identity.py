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


if __name__ == "__main__":
    unittest.main()
