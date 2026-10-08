#!/usr/bin/env python3
"""Matrix-free L^α Dirichlet audit on the combinatorial gasket.

Method matches scripts/gasket_fractional_level9.py notes:
  exact spectral action on eigenpairs with λ ≤ λ_cut;
  Chebyshev of λ^α on [λ_cut, λ_max] for the orthogonal complement.
Not thrust.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from numpy.polynomial.chebyshev import chebfit, chebval

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gasket_fractional_limit import BC, build_gasket, corner_flux_norm

ALPHA = 0.45
LAM_MAX = 8.0


def sparse_laplacian(n: int, edges):
    deg = np.zeros(n, dtype=np.float64)
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for i, j in edges:
        deg[i] += 1.0
        deg[j] += 1.0
        rows.extend([i, j])
        cols.extend([j, i])
        data.extend([-1.0, -1.0])
    for i in range(n):
        rows.append(i)
        cols.append(i)
        data.append(float(deg[i]))
    return sp.csr_matrix((data, (rows, cols)), shape=(n, n)), deg


def cheb_coeffs(alpha: float, a: float, b: float, deg: int):
    xs = 0.5 * (a + b) + 0.5 * (b - a) * np.cos(np.pi * np.arange(deg + 1) / deg)
    fs = xs**alpha
    xi = (2.0 * xs - (b + a)) / (b - a)
    coef = chebfit(xi, fs, deg)
    xt = np.geomspace(a, b, 40000)
    xit = (2.0 * xt - (b + a)) / (b - a)
    err = float(np.max(np.abs(chebval(xit, coef) - xt**alpha)))
    return coef, err


def low_eigenpairs(L: sp.csr_matrix, lam_cut: float, k_hint: int, tol: float = 1e-12):
    n = L.shape[0]
    k = min(max(k_hint, 32), n - 2)
    sigma = 1e-6
    while True:
        t0 = time.time()
        print(f"  eigsh k={k} sigma={sigma} ...", flush=True)
        w, V = sla.eigsh(L, k=k, sigma=sigma, which="LM", tol=tol)
        order = np.argsort(w)
        w = np.clip(w[order], 0.0, None)
        V = V[:, order]
        n_below = int(np.sum(w <= lam_cut + 1e-14))
        print(
            f"  got k={k} in {time.time()-t0:.1f}s; "
            f"range [{w[0]:.3e},{w[-1]:.6f}]; below_cut={n_below}",
            flush=True,
        )
        if w[-1] > lam_cut and n_below < k:
            break
        if k >= n - 2:
            break
        k = min(n - 2, max(k + max(64, k // 2), n_below + 64))
    mask = w <= lam_cut + 1e-14
    w = w[mask]
    V = V[:, mask]
    R = L @ V - V * w
    res = np.linalg.norm(R, axis=0)
    print(
        f"  kept {w.size} modes ≤ {lam_cut}; max|res|={res.max():.3e}; "
        f"w_max={w[-1]:.9f}",
        flush=True,
    )
    return w, V, float(res.max())


class FractionalOperator:
    def __init__(self, L, w_low, V_low, alpha, lam_cut, lam_max, deg):
        self.L = L
        self.V = V_low
        self.a = lam_cut
        self.b = lam_max
        self.coef, self.scalar_err = cheb_coeffs(alpha, lam_cut, lam_max, deg)
        self.c = 0.5 * (lam_max + lam_cut)
        self.d = 0.5 * (lam_max - lam_cut)
        fw = np.zeros_like(w_low)
        pos = w_low > 0
        fw[pos] = w_low[pos] ** alpha
        self.fw = fw

    def _scaled_matvec(self, x):
        return (self.L @ x - self.c * x) / self.d

    def _cheb_apply(self, x):
        """Clenshaw evaluation of sum coef[k] T_k(S) x."""
        coef = self.coef
        m = len(coef) - 1
        b1 = np.zeros_like(x)
        b2 = np.zeros_like(x)
        for k in range(m, 0, -1):
            b0 = coef[k] * x + 2.0 * self._scaled_matvec(b1) - b2
            b2 = b1
            b1 = b0
        return coef[0] * x + self._scaled_matvec(b1) - b2

    def apply(self, x):
        coeffs = self.V.T @ x
        low = self.V @ (self.fw * coeffs)
        x_orth = x - self.V @ coeffs
        return low + self._cheb_apply(x_orth)


def solve_dirichlet(op: FractionalOperator, corners, bc, tol=1e-12, maxiter=5000):
    n = op.L.shape[0]
    mask = np.ones(n, dtype=bool)
    mask[corners] = False
    interior = np.where(mask)[0]
    u = np.zeros(n, dtype=np.float64)
    u[corners] = bc

    def matvec(x_I):
        full = np.zeros(n, dtype=np.float64)
        full[interior] = x_I
        return op.apply(full)[interior]

    A = sla.LinearOperator((interior.size, interior.size), matvec=matvec)
    uB = np.zeros(n, dtype=np.float64)
    uB[corners] = bc
    rhs = -op.apply(uB)[interior]
    t0 = time.time()
    # callback to count iterations
    n_iter = [0]

    def cb(_xk):
        n_iter[0] += 1

    x, info = sla.cg(A, rhs, rtol=tol, atol=0.0, maxiter=maxiter, callback=cb)
    print(f"  CG info={info} iters={n_iter[0]} time={time.time()-t0:.1f}s", flush=True)
    if info != 0:
        raise RuntimeError(f"CG failed info={info}")
    u[interior] = x
    return u


def audit_at_level(level, lam_cut, deg, k_hint):
    print(f"=== level {level} cut={lam_cut} deg={deg} ===", flush=True)
    t_all = time.time()
    P, edges, corners = build_gasket(level)
    n = len(P)
    L, degv = sparse_laplacian(n, edges)
    assert float(degv.max()) <= 4.0 + 1e-9
    print(f"  N={n} edges={len(edges)}", flush=True)
    w, V, max_res = low_eigenpairs(L, lam_cut, k_hint)
    op = FractionalOperator(L, w, V, ALPHA, lam_cut, LAM_MAX, deg)
    print(f"  cheb scalar_err={op.scalar_err:.3e}", flush=True)
    u = solve_dirichlet(op, corners, BC)
    Au = op.apply(u)
    interior = np.setdiff1d(np.arange(n), corners)
    res = float(np.max(np.abs(Au[interior])))
    I = (L @ u)[corners]
    If = Au[corners]
    nF = corner_flux_norm(P, corners, I)
    nFf = corner_flux_norm(P, corners, If)
    Q = float(u @ Au)
    print(
        f"  ||F||={nF:.12f} ||F_frac||={nFf:.12f} Q={Q:.12f} "
        f"sumI={float(I.sum()):.3e} res={res:.3e} "
        f"u=[{u.min():.6f},{u.max():.6f}] time={time.time()-t_all:.1f}s",
        flush=True,
    )
    print(f"  I={I}", flush=True)
    # Schur / centroid checks
    ia = float(I[0])
    schur_err = max(abs(float(I[1]) + 0.8 * ia), abs(float(I[2]) + 0.2 * ia))
    cent_err = abs(nF - ia * (21.0**0.5) / 5.0)
    print(f"  schur_err={schur_err:.3e} centroid_err={cent_err:.3e}", flush=True)
    return {
        "level": level,
        "N": n,
        "nF": float(nF),
        "nFf": float(nFf),
        "Q": float(Q),
        "I": np.array(I, dtype=np.float64),
        "If": np.array(If, dtype=np.float64),
        "res": res,
        "umin": float(u.min()),
        "umax": float(u.max()),
        "max_eig_res": max_res,
        "scalar_err": op.scalar_err,
        "n_modes": int(w.size),
        "lam_cut": lam_cut,
        "schur_err": schur_err,
        "cent_err": cent_err,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--cut", type=float, default=0.05)
    ap.add_argument("--deg", type=int, default=150)
    ap.add_argument("--k-hint", type=int, default=256)
    args = ap.parse_args()
    rec = audit_at_level(args.level, args.cut, args.deg, args.k_hint)
    print(
        "RECORD",
        f"F={rec['nF']:.12f}",
        f"Q={rec['Q']:.12f}",
        f"Ff={rec['nFf']:.12f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
