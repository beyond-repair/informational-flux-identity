#!/usr/bin/env python3
"""Theorem FH: monotonicity in n for every n >= 11, G_(n+1)(alpha) > G_n(alpha) for every real alpha in [-1, 1)
(FF's kept limit 'for n >= 11, FB.2/FB.4 do not give D_n > D_(n+1) (needs m > C/3), monotonicity for n >= 11 stays OPEN').

Setting (FB): G_n = rho^-n H_n = int t^-alpha f_n, W/2 = int t^-alpha f, D_n = W/2 - G_n = int t^-alpha (f - f_n),
G_(n+1) - G_n = D_n - D_(n+1).  FB.1: f_n/f = e^Delta/(1+E), E = (t/6)3^-n (1/P_n + (3/5)^n), 0 <= log(P/P_n) <= Delta <= Lam_n(t),
Lam_n(t) <= 5^-n t beta_n M (w = 5^-n t <= 1), M = max(1,t)^gam, gam = log_5 2, 1/P = 3R; FB.2: Lam_n <= 0.24 on t <= 3^n (n >= 9),
so e^Lam - 1 <= chi* Lam; tau_n = 3^n int_(3^n)^oo t (f + f_n) <= tau_11max for n >= 11.

  FH.1 (THEOREM; a sharp every-level UPPER bound).  For real alpha in [a0, a1] subset [-1, 1) and every m >= 12:
        3^m D_m(alpha) <= U = Vup + (chi*/2) 5^-12 beta_12 J + (3/5)^12 K0/6 + tau_11max,
     Vup = (1/2) int t^(1-a') f R, J = int t^(2-a') f R M, K0 = int t^(1-a') f, with a' = a1 on t < 1 and a' = a0 on t > 1 (sup over the cell).
     Proof: f_n/f >= 1/(1+E) >= 1 - E, so f - f_n <= f E; 1/P_n = (1/P)(P/P_n) <= 3R e^Lam; on t <= 3^n
     3^n f E <= f[(t/2) R (1 + chi* 5^-n t beta_n M) + (t/6)(3/5)^n]; beyond 3^n, 3^n(f - f_n) <= 3^n f (tau).  5^-m beta_m, (3/5)^m decrease.
     (FB.2 had the uniform C <= 1095.75: there 1/P_n <= 1 and e^Lam - 1 multiplied 3^n instead of E.)
  FH.2 (THEOREM; lower bound, FB.4 on a finer partition with the exact M).  For every n >= 11:
        3^n D_n(alpha) >= L = Vlo/(1+E_up) - del_11 K - tau_11max,  K = int t^(1-a') M f (a' as in FH.1), del_11 = (3/5)^11 chi* beta_11,
     Vlo FB.4's lower piece (t <= 5^5; a0 on t < 1, a1 on t > 1).  (FB.4 used M <= max(1,t)^(1/2) <= 1 + t^(1/2), i.e. (W(a0-1)+W(a0-3/2))/2.)
  FH.3 (THEOREM, certified arb).  On every cell of the partition of [-1, 1] into 40 pieces of width 1/20: L > U/3.  Hence for n >= 11:
        3^n (G_(n+1) - G_n) = 3^n D_n - 3^(n+1) D_(n+1)/3 >= L - U/3 > 0.
     With FF.1 (n = 2..10) and FG.4 (n = 1): G_2 < G_3 < ... (every n >= 2) on [-1, 1), and G_1 < G_2 < ... on (beta_c, 1).
  FH.4 (THEOREM, consequence).  L <= 3^n D_n <= U on each cell for n >= 12 (two-sided, both close to V(alpha)).
Tails (t >= 5^8): EX.3's bound as used in FB's W_cert, int_(5^8)^oo t^nu R/(5+psi(t/5)) <= 2 T8(nu), T8 = R(5^8) 4 5^8 5^(9 nu)/psi(5^7),
valid when 5^(1+nu)/psi(5^7) <= 1/2; weights t^gam <= t^(1/2), R <= 1, f <= R/(5+psi) there.  Small t: f <= 1/30, R <= 1/3.
Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, __file__.rsplit("/", 1)[0])

import gasket_fractional_corner_h_alpha_sharp_rate as fb  # noqa: E402

PART = [Fr(-1) + Fr(k, 20) for k in range(41)]
EPS_EXP, HI_EXP = -24, 8
TOL = 30


def _ctx():
    ev, acb, arb, RY, psi = fb.arb_setup()
    return acb, arb, RY, psi


def cell_upper(args):
    """Upper bounds (Vup, J, K0, K) on the cell [a0, a1] and the lower Vlo (FB.4's routine)."""
    a0, a1 = args
    acb, arb, RY, psi = _ctx()
    gam = arb(2).log() / arb(5).log()
    A0 = arb(a0.numerator) / a0.denominator
    A1 = arb(a1.numerator) / a1.denominator
    R8 = RY(5 ** 8)[0].real.upper()
    p7 = psi(5 ** 7).real.lower()
    eps = arb(5) ** EPS_EXP

    def make(kind, small):
        ex = A1 if small else A0

        def F(z, analytic):
            R, y1 = RY(z)
            f = R / (2 * (5 + y1))
            lz = z.log()
            if kind == "V":
                g = f * R / 2 * ((1 - ex) * lz).exp()
            elif kind == "J":
                g = f * R * ((2 - ex) * lz).exp()
            elif kind == "K0":
                g = f * ((1 - ex) * lz).exp()
            else:  # K
                g = f * ((1 - ex) * lz).exp()
            if kind in ("J", "K") and not small:
                g = g * (gam * lz).exp()
            return g
        return F

    def tail(nu: Fr):
        nub = arb(nu.numerator) / nu.denominator
        assert bool(arb(5) ** (1 + nub) / p7 <= arb(1) / 2)
        return (2 * R8 * 4 * arb(5) ** 8 * arb(5) ** (nub * 9) / p7).upper()

    out = {}
    for kind, e_small, nu, rfac in (("V", 1 - a1, 1 - a0, Fr(1, 6)), ("J", 2 - a1, Fr(5, 2) - a0, Fr(1, 3)),
                                    ("K0", 1 - a1, 1 - a0, Fr(1)), ("K", 1 - a1, Fr(3, 2) - a0, Fr(1))):
        tot = arb(0)
        for m in range(EPS_EXP, HI_EXP):
            F = make(kind, m < 0)
            tot += acb.integral(F, arb(5) ** m, arb(5) ** (m + 1), rel_tol=arb(2) ** -TOL, abs_tol=arb(2) ** -(TOL + 10)).real
        e1 = arb(e_small.numerator) / e_small.denominator + 1
        small = eps ** e1 / e1 / 30 * (arb(1) / 3 if kind in ("V", "J") else 1)
        tl = tail(nu) * (arb(rfac.numerator) / rfac.denominator)
        out[kind] = (tot.upper() + small.upper() + tl).upper()
    out["Vlo"] = fb.V_piece_lower(acb, arb, RY, a0, -24, 0) + fb.V_piece_lower(acb, arb, RY, a1, 0, fb.X_EXP)
    return a0, a1, {k: v.str(30, radius=False) if hasattr(v, "str") else str(v) for k, v in out.items()}


def float_margin(a, N=3000):
    """Float (non-rigorous) pointwise L/(U/3) at alpha = a (CI sanity)."""
    import gasket_fractional_corner_h_alpha_rate as ey
    gam = math.log(2) / math.log(5)
    cs, kap = 4 - math.sqrt(15), math.exp(math.exp(0.1) / 4)
    beta = lambda n: (cs * n + 1) * kap / 10 + math.exp(0.02) / 30  # noqa: E731
    chi = (math.exp(0.24) - 1) / 0.24
    x0, x1 = math.log(1e-14), math.log(5.0 ** 9)
    h = (x1 - x0) / N
    V = Vlo = J = K0 = K = 0.0
    for i in range(N + 1):
        t = math.exp(x0 + i * h)
        w = (1 if i in (0, N) else (4 if i % 2 else 2)) * t
        R = 1 / ey.inv_R(t)
        f = R / (2 * (5 + ey.psi(t / 5)))
        M = max(1.0, t) ** gam
        V += w * t ** (1 - a) * f * R / 2
        Vlo += w * t ** (1 - a) * f * R / 2 if t <= 5 ** 5 else 0
        J += w * t ** (2 - a) * f * R * M
        K0 += w * t ** (1 - a) * f
        K += w * t ** (1 - a) * f * M
    V, Vlo, J, K0, K = (x * h / 3 for x in (V, Vlo, J, K0, K))
    Eup = 5 ** 5 / 6 * 3.0 ** -11 * (1 + 0.6 ** 11)
    L = Vlo / (1 + Eup) - 0.6 ** 11 * chi * beta(11) * K
    U = V + chi / 2 * 5.0 ** -12 * beta(12) * J + 0.6 ** 12 * K0 / 6
    return L, U, V


def main() -> int:
    import multiprocessing as mp
    import os
    t0 = time.time()
    print("Theorem FH: G_(n+1) > G_n for every n >= 11 and real alpha in [-1,1). Functions on a graph, not a flux, not thrust.", flush=True)
    fails = []

    def check(ok, msg):
        print(("  ok   " if ok else "  FAIL ") + msg, flush=True)
        if not ok:
            fails.append(msg)

    from flint import arb as A, ctx
    ctx.prec = 128
    cs = 4 - A(15).sqrt()
    kap = ((A(1) / 10).exp() / 4).exp()
    gam = A(2).log() / A(5).log()
    q = A(3) ** (1 + gam) / 5
    A_ = cs * kap / 10
    B_ = kap / 10 + (A(1) / 50).exp() / 30

    def beta(n):
        return A_ * n + B_

    Lmax = max((q ** n * beta(n)).upper() for n in range(9, 61))
    check(bool(Lmax < A("0.24")) and bool(q * (1 + A_ / (A_ * 27 + B_)) < 1), "FB.2 reused: sup_(n>=9) Lam_n(3^n) <= q^n beta_n < 0.24")
    chi = (A("0.24").exp() - 1) / A("0.24")
    check(bool(A(1) / 5 * (1 + A_ / (A_ * 12 + B_)) < 1), "FH.1: 5^-m beta_m decreasing for m >= 12")
    check(bool(A(3) / 5 * (1 + A_ / (A_ * 11 + B_)) < 1), "FH.2: (3/5)^n beta_n decreasing for n >= 11")

    def tau(n):
        J, _ = fb.best_J(n)
        s0 = A(3) ** n / A(5) ** (J + 1)
        assert bool(s0 > 1) and 2 <= J <= n - 1
        return (A(3) ** n * A(5) ** (2 * (J + 1)) * s0 ** (2 - 2 ** J) / (2 ** J - 2) / 3).upper()

    tmax = max(tau(n) for n in range(11, 30))
    check(all(fb.tail_log_analytic_ok(n) for n in range(30, 400)) and bool(tmax < A("1e-15")),
          f"tau_n <= {tmax.str(5, radius=False)} for 11 <= n <= 29 (arb), < e^-1000 for n >= 30 (FB.2)")
    Eup = A(5) ** fb.X_EXP / 6 / A(3) ** 11 * (1 + A(3) ** 11 / A(5) ** 11)
    d11 = (A(3) / 5) ** 11 * chi * beta(11)
    cu = chi / 2 * A(5) ** -12 * beta(12)
    q12 = (A(3) / 5) ** 12 / 6
    print(f"   del_11 <= {d11.upper().str(8, radius=False)}, (chi/2)5^-12 beta_12 <= {cu.upper().str(6, radius=False)}, E_up <= {Eup.upper().str(6, radius=False)}", flush=True)

    jobs = int(os.environ.get("FH_JOBS", "6"))
    cells = list(zip(PART, PART[1:]))
    with mp.Pool(jobs) as pool:
        res = pool.map(cell_upper, cells, chunksize=1)
    worst = None
    for a0, a1, d in res:
        g = {k: A(v) for k, v in d.items()}
        U = g["V"] + cu * g["J"] + q12 * g["K0"] + tmax
        L = g["Vlo"] / (1 + Eup) - d11 * g["K"] - tmax
        marg = (L - U / 3).lower()
        worst = marg if worst is None or bool(marg < worst) else worst
        print(f"   [{float(a0):+.2f},{float(a1):+.2f}]  L >= {L.lower().str(7, radius=False)}  U <= {U.upper().str(7, radius=False)}  L - U/3 >= {marg.str(6, radius=False)}", flush=True)
        check(bool(marg > 0), f"FH.3 cell [{a0}, {a1}]: 3^n(G_(n+1) - G_n) >= L - U/3 > 0 for every n >= 11")
    print(f"   worst margin min_cells (L - U/3) >= {worst.str(6, radius=False)}  ({time.time() - t0:.0f}s)", flush=True)
    for a in (-1.0, 0.0, 0.9):
        L, U, V = float_margin(a)
        print(f"   float alpha = {a}: L ~ {L:.6g}, V ~ {V:.6g}, U ~ {U:.6g}", flush=True)
    print(f"elapsed {time.time() - t0:.1f}s", flush=True)
    if fails:
        print(f"{len(fails)} FAILED")
        return 1
    print("all FH checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
