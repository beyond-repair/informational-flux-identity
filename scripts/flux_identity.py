#!/usr/bin/env python3
"""Exact integer witness for the closed-rectangle flux identity.

Prints the numbers recorded in README.md. Exits nonzero if a recorded
figure changes.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


NX = NY = 16
ROW_NM2 = [0, -4, -4, 3, 2, 3, 0, 3, -2, 0, 3, -3, -2, -3, 0, 4, -3]
ROW_NM1 = [-2, 1, 16, 27, -28, -22, 20, 27, -15, -11, 23, -5, -14, 20, -15, -6, 9]


def potential() -> np.ndarray:
    p = np.zeros((NX + 1, NY + 1), dtype=np.int64)
    p[-2, :] = np.array(ROW_NM2, dtype=np.int64)
    p[-1, :] = np.array(ROW_NM1, dtype=np.int64)
    return p


def fluxes(p: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    fx = p[:, 1:] - p[:, :-1]
    fy = p[:-1, :] - p[1:, :]
    return fx, fy


def divergence(fx: np.ndarray, fy: np.ndarray) -> np.ndarray:
    return (fx[1:, :] - fx[:-1, :]) + (fy[:, 1:] - fy[:, :-1])


def ledger(fx: np.ndarray, fy: np.ndarray) -> dict:
    signed = {
        "left": int((-fx[0, :]).sum()),
        "right": int(fx[-1, :].sum()),
        "bottom": int((-fy[:, 0]).sum()),
        "top": int(fy[:, -1].sum()),
    }
    absolute = {
        "left": int(np.abs(fx[0, :]).sum()),
        "right": int(np.abs(fx[-1, :]).sum()),
        "bottom": int(np.abs(fy[:, 0]).sum()),
        "top": int(np.abs(fy[:, -1]).sum()),
    }
    div = divergence(fx, fy)
    net = sum(signed.values())
    abs_total = sum(absolute.values())
    return {
        "div_max_abs": int(np.max(np.abs(div))),
        "div_sum": int(div.sum()),
        "net": net,
        "signed": signed,
        "absolute": absolute,
        "absolute_total": abs_total,
        "aft_absolute_fraction": [absolute["right"], abs_total],
    }


def main() -> None:
    p = potential()
    fx, fy = fluxes(p)
    div_free = ledger(fx, fy)

    fx_src = fx.copy()
    fy_src = fy.copy()
    fx_src[-1, 3] = fx_src[-1, 3] + 17
    sourced = ledger(fx_src, fy_src)

    # 1D product rule, exact integers.
    # sum_{i=0}^{n-1} s[i] * (w[i+1] - w[i]) 
    #   = s[n-1]*w[n] - s[0]*w[0] - sum_{i=1}^{n-1} (s[i] - s[i-1]) * w[i]
    # s has one entry per edge; w has one entry per vertex. len(w) == len(s) + 1.
    # In one dimension a vanishing edge-difference means s is constant, so the
    # only compact divergence-free flux is s = 0. The nontrivial vanishing
    # witness is the two-dimensional example above.
    s_const = np.array([4, 4, 4, 4], dtype=np.int64)
    w_matched = np.array([3, 8, -1, 6, 3], dtype=np.int64)  # ends equal
    w_open = np.array([3, 8, -1, 6, 9], dtype=np.int64)  # ends differ
    s_var = np.array([0, 4, -2, 5], dtype=np.int64)
    w_var = np.array([3, 1, 4, -5, 2], dtype=np.int64)

    def pairing(sv: np.ndarray, wv: np.ndarray) -> dict:
        if len(wv) != len(sv) + 1:
            raise ValueError("w must be one longer than s")
        vol = int(np.sum(sv * (wv[1:] - wv[:-1])))
        expanded = int(sv[-1] * wv[-1] - sv[0] * wv[0] - np.sum((sv[1:] - sv[:-1]) * wv[1:-1]))
        return {"volume": vol, "integration_by_parts": expanded, "match": vol == expanded}

    one_d = {
        "constant_flux_equal_ends": pairing(s_const, w_matched)["volume"],
        "constant_flux_unequal_ends": pairing(s_const, w_open)["volume"],
        "variable_flux_match": pairing(s_var, w_var)["match"],
        "constant_equal_match": pairing(s_const, w_matched)["match"],
        "constant_unequal_match": pairing(s_const, w_open)["match"],
    }

    witness = {"divergence_free": div_free, "boundary_source_17": sourced, "one_dimensional": one_d}
    out = Path(__file__).resolve().parents[1] / "witness.json"
    out.write_text(json.dumps(witness, indent=2) + "\n", encoding="utf-8")

    assert div_free["div_max_abs"] == 0
    assert div_free["div_sum"] == 0
    assert div_free["net"] == 0
    assert div_free["signed"] == {"left": 0, "right": 11, "bottom": -2, "top": -9}
    assert div_free["absolute"] == {"left": 0, "right": 349, "bottom": 2, "top": 15}
    assert div_free["absolute_total"] == 366
    assert sourced["div_sum"] == 17
    assert sourced["net"] == 17
    assert one_d["constant_equal_match"] and one_d["constant_unequal_match"] and one_d["variable_flux_match"]
    assert one_d["constant_flux_equal_ends"] == 0
    assert one_d["constant_flux_unequal_ends"] == 4 * (9 - 3)
    print(json.dumps(witness, indent=2))


if __name__ == "__main__":
    main()
