"""Sweep-168: lock the rectangle identity without rewriting witness.json."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import flux_identity as fi


def test_divergence_free_witness_matches_readme():
    p = fi.potential()
    fx, fy = fi.fluxes(p)
    ledger = fi.ledger(fx, fy)
    assert ledger["div_max_abs"] == 0
    assert ledger["div_sum"] == 0
    assert ledger["net"] == 0
    assert ledger["signed"] == {"left": 0, "right": 11, "bottom": -2, "top": -9}
    assert ledger["absolute"] == {"left": 0, "right": 349, "bottom": 2, "top": 15}
    assert ledger["absolute_total"] == 366


def test_theorem_a_on_random_integer_flux():
    rng = np.random.default_rng(168)
    fx = rng.integers(-5, 6, size=(8, 7), dtype=np.int64)
    fy = rng.integers(-5, 6, size=(7, 8), dtype=np.int64)
    div = fi.divergence(fx, fy)
    signed = int(fx[-1, :].sum() - fx[0, :].sum() + fy[:, -1].sum() - fy[:, 0].sum())
    assert int(div.sum()) == signed


def test_boundary_source_matches_added_edge():
    p = fi.potential()
    fx, fy = fi.fluxes(p)
    fx = fx.copy()
    fx[-1, 3] = fx[-1, 3] + 17
    sourced = fi.ledger(fx, fy)
    assert sourced["div_sum"] == 17
    assert sourced["net"] == 17
