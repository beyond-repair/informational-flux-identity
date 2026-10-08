#!/usr/bin/env python3
"""Theorem M witness: the Hessian reading of (grad Psi)^{ij} reduces the
closed-surface flux G_i to the boundary values of the Laplacian of Psi.

Integer arithmetic only (stdlib). Exits nonzero if any identity or recorded
integer moves. Nothing here is thrust.

Grid: Nx x Ny cells. P is an integer potential on vertices
x in [-1, Nx+1], y in [-1, Ny+1]. For row i, the cell field g^i is the forward
difference D_i P, and the face fluxes are F^x = g(x,y) - g(x-1,y),
F^y = g(x,y) - g(x,y-1), so row i of the array is H^{ij} = D_j D_i P.
"""
import random
import sys

NX, NY = 8, 8


def lap(P, x, y):
    return P[x + 1][y] + P[x - 1][y] + P[x][y + 1] + P[x][y - 1] - 4 * P[x][y]


class Grid:
    """Potential stored with an offset of 1 so index -1 is valid."""

    def __init__(self, f):
        self.v = {(x, y): f(x, y) for x in range(-1, NX + 2) for y in range(-1, NY + 2)}

    def __getitem__(self, x):
        v = self.v
        return _Col(v, x)


class _Col:
    def __init__(self, v, x):
        self.v, self.x = v, x

    def __getitem__(self, y):
        return self.v[(self.x, y)]


def g_field(P, i):
    if i == "x":
        return lambda x, y: P[x + 1][y] - P[x][y]
    return lambda x, y: P[x][y + 1] - P[x][y]


def row_flux(P, i):
    """Face fluxes of row i, cell divergence sum, signed outward flux, abs pieces."""
    g = g_field(P, i)
    Fx = {(xf, y): g(xf, y) - g(xf - 1, y) for xf in range(0, NX + 1) for y in range(NY)}
    Fy = {(x, yf): g(x, yf) - g(x, yf - 1) for x in range(NX) for yf in range(0, NY + 1)}
    div = {
        (x, y): Fx[(x + 1, y)] - Fx[(x, y)] + Fy[(x, y + 1)] - Fy[(x, y)]
        for x in range(NX)
        for y in range(NY)
    }
    faces = {
        "left": -sum(Fx[(0, y)] for y in range(NY)),
        "right": sum(Fx[(NX, y)] for y in range(NY)),
        "bottom": -sum(Fy[(x, 0)] for x in range(NX)),
        "top": sum(Fy[(x, NY)] for x in range(NX)),
    }
    absf = {
        "left": sum(abs(Fx[(0, y)]) for y in range(NY)),
        "right": sum(abs(Fx[(NX, y)]) for y in range(NY)),
        "bottom": sum(abs(Fy[(x, 0)]) for x in range(NX)),
        "top": sum(abs(Fy[(x, NY)]) for x in range(NX)),
    }
    return div, sum(faces.values()), faces, absf


def boundary_laplacian(P, i):
    if i == "x":
        return sum(lap(P, NX, y) - lap(P, 0, y) for y in range(NY))
    return sum(lap(P, x, NY) - lap(P, x, 0) for x in range(NX))


def check(P, label):
    out = {}
    for i in ("x", "y"):
        div, signed, faces, absf = row_flux(P, i)
        sdiv = sum(div.values())
        # Theorem A on row i, then the commutation D_j D_j D_i = D_i Lap.
        assert signed == sdiv, (label, i, signed, sdiv)
        for (x, y), d in div.items():
            d2 = (lap(P, x + 1, y) - lap(P, x, y)) if i == "x" else (lap(P, x, y + 1) - lap(P, x, y))
            assert d == d2, (label, i, x, y)
        bl = boundary_laplacian(P, i)
        assert signed == bl, (label, i, signed, bl)
        out[i] = {
            "G": signed,
            "boundary_laplacian": bl,
            "sum_abs_cell_div": sum(abs(d) for d in div.values()),
            "abs_faces": absf,
        }
    return out


def main():
    ok = True
    # 1. Arbitrary integer potentials: identity holds every time.
    rng = random.Random(20261008)
    for _ in range(200):
        P = Grid(lambda x, y: rng.randint(-9, 9))
        check(P, "random")

    # 2. Lopsided interior source, Laplacian zero on the boundary layers.
    lump = {(2, 3): 7, (3, 3): 11, (3, 4): -2, (6, 5): 13, (6, 6): 5}
    P2 = Grid(lambda x, y: lump.get((x, y), 0))
    assert all(lap(P2, x, y) == 0 for x in (0, NX) for y in range(NY))
    assert all(lap(P2, x, y) == 0 for y in (0, NY) for x in range(NX))
    r2 = check(P2, "lump")
    print("lopsided interior source:", {k: (v["G"], v["sum_abs_cell_div"]) for k, v in r2.items()})
    if not (r2["x"]["G"] == 0 and r2["y"]["G"] == 0 and r2["x"]["sum_abs_cell_div"] > 0 and r2["y"]["sum_abs_cell_div"] > 0):
        ok = False

    # 3. Laplacian nonzero on the surface: P = x^3, discrete Lap = 6x.
    P3 = Grid(lambda x, y: x ** 3)
    r3 = check(P3, "cubic")
    print("P = x^3:", r3["x"]["G"], r3["y"]["G"], "expected", 6 * NX * NY, 0)
    if not (r3["x"]["G"] == 6 * NX * NY == 384 and r3["y"]["G"] == 0):
        ok = False

    # 4. Discrete screened ("Proca-like") surface: Lap P = m2 * P on the boundary
    #    layers only, with m2 = 1, built on the right column. G_x is then
    #    m2 * sum_y P(NX, y) - m2 * sum_y P(0, y) by the identity.
    #    Use P = 0 except column NX+1 set so that Lap P(NX,y) = P(NX,y) with P(NX,y)=1.
    #    Lap(NX,y) = P(NX+1,y) + P(NX-1,y) + P(NX,y+1) + P(NX,y-1) - 4.
    def f4(x, y):
        if x == NX and 0 <= y <= NY - 1:
            return 1
        if x == NX + 1 and 0 <= y <= NY - 1:
            # neighbours in column NX inside rows 0..NY-1
            nb = (1 if 0 <= y + 1 <= NY - 1 else 0) + (1 if 0 <= y - 1 <= NY - 1 else 0)
            return 5 - nb
        return 0
    P4 = Grid(f4)
    assert all(lap(P4, NX, y) == P4[NX][y] for y in range(NY))
    r4 = check(P4, "screened")
    expected4 = sum(P4[NX][y] for y in range(NY)) - sum(P4[0][y] for y in range(NY))
    print("screened surface, m2=1:", r4["x"]["G"], "expected", expected4)
    if r4["x"]["G"] != expected4 or expected4 != NY:
        ok = False

    print("Theorem M: G_i = boundary sum of Lap P n_i on 200 random grids and 3 named witnesses")
    print("Not thrust.")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
