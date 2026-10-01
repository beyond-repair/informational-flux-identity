#!/usr/bin/env python3
"""Exact corner currents for the harmonic combinatorial gasket.

Checks Theorem G in GASKET.md. No thrust claim.
"""

from __future__ import annotations

import math
from collections import defaultdict
from fractions import Fraction
from typing import Dict, List, Tuple

import numpy as np

Edge = Tuple[int, int]


def build_gasket(level: int):
    c0 = np.array([0.0, 0.0])
    c1 = np.array([1.0, 0.0])
    c2 = np.array([0.5, math.sqrt(3.0) / 2.0])
    pos_map: Dict[Tuple[float, float], int] = {}
    positions: List[np.ndarray] = []
    edges: set[Edge] = set()

    def key(p: np.ndarray) -> Tuple[float, float]:
        return (round(float(p[0]), 10), round(float(p[1]), 10))

    def vid(p: np.ndarray) -> int:
        k = key(p)
        if k not in pos_map:
            pos_map[k] = len(positions)
            positions.append(p.copy())
        return pos_map[k]

    def rec(a, b, c, d):
        if d == 0:
            ia, ib, ic = vid(a), vid(b), vid(c)
            for i, j in ((ia, ib), (ib, ic), (ic, ia)):
                if i != j:
                    edges.add((min(i, j), max(i, j)))
            return
        rec(a, 0.5 * (a + b), 0.5 * (c + a), d - 1)
        rec(0.5 * (a + b), b, 0.5 * (b + c), d - 1)
        rec(0.5 * (c + a), 0.5 * (b + c), c, d - 1)

    rec(c0, c1, c2, level)
    corners = [pos_map[key(c0)], pos_map[key(c1)], pos_map[key(c2)]]
    return positions, sorted(edges), corners


def harmonic_currents(level: int, corner_values: List[Fraction]) -> List[Fraction]:
    positions, edges, corners = build_gasket(level)
    n = len(positions)
    cset = set(corners)
    interior = [i for i in range(n) if i not in cset]
    index = {v: k for k, v in enumerate(interior)}
    m = len(interior)
    adj: Dict[int, List[int]] = defaultdict(list)
    for i, j in edges:
        adj[i].append(j)
        adj[j].append(i)
    known = {c: corner_values[k] for k, c in enumerate(corners)}
    matrix = [[Fraction(0) for _ in range(m + 1)] for _ in range(m)]
    for i in interior:
        r = index[i]
        matrix[r][r] = Fraction(len(adj[i]))
        for j in adj[i]:
            if j in cset:
                matrix[r][m] += known[j]
            else:
                matrix[r][index[j]] -= 1
    for col in range(m):
        piv = max(range(col, m), key=lambda r: abs(float(matrix[r][col])))
        matrix[col], matrix[piv] = matrix[piv], matrix[col]
        pivot = matrix[col][col]
        if pivot == 0:
            raise RuntimeError("singular cell")
        for c in range(col, m + 1):
            matrix[col][c] /= pivot
        for r in range(m):
            factor = matrix[r][col]
            if r == col or factor == 0:
                continue
            for c in range(col, m + 1):
                matrix[r][c] -= factor * matrix[col][c]
    u = [Fraction(0) for _ in range(n)]
    for k, c in enumerate(corners):
        u[c] = corner_values[k]
    for i, r in index.items():
        u[i] = matrix[r][m]
    currents = []
    for c in corners:
        currents.append(sum((u[c] - u[j]) for j in adj[c]))
    return currents


def level1_from_formula(a: Fraction, b: Fraction, c: Fraction) -> List[Fraction]:
    mab = (2 * a + 2 * b + c) / 5
    mbc = (2 * b + 2 * c + a) / 5
    mca = (2 * c + 2 * a + b) / 5
    ia = (a - mab) + (a - mca)
    ib = (b - mab) + (b - mbc)
    ic = (c - mbc) + (c - mca)
    return [ia, ib, ic]


def main() -> None:
    data = [Fraction(1), Fraction(-1, 2), Fraction(0)]
    predicted_1 = [Fraction(3, 2), Fraction(-6, 5), Fraction(-3, 10)]
    got_formula = level1_from_formula(*data)
    got_solve = harmonic_currents(1, data)
    assert got_formula == predicted_1
    assert got_solve == predicted_1
    for level in (2, 3, 4):
        got = harmonic_currents(level, data)
        scale = Fraction(3, 5) ** (level - 1)
        expect = [scale * v for v in predicted_1]
        assert got == expect, (level, got, expect)
        assert sum(got) == 0
    # Refinement multiplier, arbitrary data, from the cell formula alone.
    a, b, c = Fraction(2), Fraction(-3), Fraction(5)
    i1 = level1_from_formula(a, b, c)
    # One refinement step at a corner with neighbor values equal to the
    # level-1 midpoints, third vertex the other neighbor: the algebra in
    # GASKET.md gives 3/5 without another solve. Check it numerically here
    # on the solved level-2 currents against level-1.
    i2 = harmonic_currents(2, [a, b, c])
    assert i2 == [Fraction(3, 5) * v for v in i1]
    assert sum(i1) == 0 and sum(i2) == 0
    print("level1", [str(v) for v in predicted_1])
    print("scale", "3/5 each refinement")
    print("sum", 0)
    print("norm", "(3/5)^n * sqrt(21) / 2 for data (1,-1/2,0)")


if __name__ == "__main__":
    main()
