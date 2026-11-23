#!/usr/bin/env python3
"""Theorem FG: the n = 1 step of monotonicity in n, i.e. the sign of G_2 - G_1 on (-13/20, -3/8) (FF.3 OPEN).

Setting (FC/FF): G_n(alpha) = rho^-n H_n(alpha) = int_0^oo t^-alpha f_n dt, h_n = f - f_n,
  Delta(alpha) := G_2(alpha) - G_1(alpha) = int_0^oo t^-alpha (h_1 - h_2) dt = int t^-alpha (f_2 - f_1) dt   (real alpha in (-1, 1)).
FE (D_1, D_2 holomorphic on -1 < Re alpha < 1, derivatives under the integral) gives
  Delta'(alpha) = D_1' - D_2' = int_0^oo (-ln t) t^-alpha (h_1 - h_2) dt.

  FG.1 (THEOREM, certified arb).  Delta(B_LO) < 0 < Delta(B_HI), B_LO = -0.5679, B_HI = -0.5677 (enclosures printed).
  FG.2 (THEOREM, certified arb).  Delta'(alpha) >= m > 0 for every alpha in [-13/20, -1/2] (m printed).  Integrand is monotone in alpha pointwise,
        so the min over [a0, a1] is acb.real_min of the endpoint values on blocks [5^k, 5^(k+1)], k = -24..6;
        |int_0^(5^-24)| <= 2 int_0^eps |ln t| t^(1-a1)/30 (FB.1: |h_n| <= t/30);
        beyond X = 5^7: (-ln t) t^-alpha (f_2 - f_1) >= ln t t^-a1 f_1 - ln t t^-a0 f_2 (t > 1, a0 <= alpha <= a1); the f_1 part is an Arb lower bound on
        [5^7, 5^60] (rest dropped, >= 0), the f_2 part an Arb upper bound on [5^7, 5^20] plus analytic tail
        int_T^oo ln t 9*5^4 t^(-a0-3) dt = 9*5^4 T^-p (ln T/p + 1/p^2), p = a0 + 2 (f_2 <= 3^2/(t psi_1(t/5)), psi_1(t/5) >= (t/25)^2, FC.2).
  FG.3 (THEOREM, certified arb).  Delta(alpha) > 0 for every alpha in [-1/2, -3/8] (FF's 'lo' evaluator, adaptive bisection; min printed).
  FG.4 (THEOREM, consequence).  With FF.2 (Delta < 0 on (-1, -13/20], Delta > 0 on [-3/8, 1)):  Delta has exactly one zero beta_c on (-1, 1),
        beta_c in (-0.5679, -0.5677), simple (Delta'(beta_c) >= m > 0);  G_2 < G_1 on (-1, beta_c), G_2 > G_1 on (beta_c, 1).
        FF.3 is settled: the n = 1 step of monotonicity in n holds exactly for alpha > beta_c.
Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, __file__.rsplit("/", 1)[0])

import gasket_fractional_corner_h_alpha_below_small_n as fc  # noqa: E402
import gasket_fractional_corner_h_alpha_monotone_n as ff  # noqa: E402

B_LO, B_HI = Fr(-5679, 10000), Fr(-5677, 10000)
D_LO, D_HI = ff.N1_NEG, Fr(-1, 2)      # Delta' > 0 on [-13/20, -1/2]
P_LO, P_HI = Fr(-1, 2), ff.N1_POS      # Delta > 0 on [-1/2, -3/8]
F2_T = 20                              # f_2 integrated by Arb to 5^20


def Delta_float(a, N=4000, hexp=40):
    """Float (non-rigorous) (Delta(a), Delta'(a)) by Simpson in log t on [1e-14, 5^hexp] plus the 15/t^2 asymptote of f_1 beyond (CI sanity)."""
    import math
    import gasket_fractional_corner_h_alpha_rate as ey
    x0, x1 = math.log(1e-14), hexp * math.log(5)
    h = (x1 - x0) / N
    s = d = 0.0
    for i in range(N + 1):
        t = math.exp(x0 + i * h)
        w = 1 if i in (0, N) else (4 if i % 2 else 2)
        g = w * t ** (1 - a) * (ey.f_n(2, t) - ey.f_n(1, t))
        s += g
        d += g * (-math.log(t))
    T, p = 5.0 ** hexp, 1 + a
    return s * h / 3 - 15 * T ** -p / p, d * h / 3 - 15 * T ** -p * (math.log(T) / p + 1 / p ** 2)


def dcore(a0, a1):
    """int_(5^-24)^(5^7) of min over alpha in [a0,a1] of (-ln t) t^-alpha (h_1 - h_2)."""
    from flint import arb
    A0, A1 = fc.q(a0), fc.q(a1)

    def F(z, analytic):
        g = fc.h_eval(1, z)[0] - fc.h_eval(2, z)[0]
        lz = z.log()
        base = g * (-lz)
        return (base * (-A0 * lz).exp()).real_min(base * (-A1 * lz).exp(), analytic=analytic)

    tot = arb(0)
    for m in range(fc.EPS_EXP, fc.X_EXP):
        tot += fc._integ(F, arb(5) ** m, arb(5) ** (m + 1), fc.TOL)
    return tot


def dtail_lower(a0, a1):
    """Lower bound of int_(5^7)^oo (-ln t) t^-alpha (f_2 - f_1) dt for alpha in [a0, a1]."""
    from flint import arb
    A0, A1 = fc.q(a0), fc.q(a1)

    def F1(z, analytic):
        return fc.fn_eval(1, z) * (-A1 * z.log()).exp() * z.log()

    def F2(z, analytic):
        return fc.fn_eval(2, z) * (-A0 * z.log()).exp() * z.log()

    t1 = arb(0)
    for m in range(fc.X_EXP, fc.T_EXP[1]):
        t1 += fc._integ(F1, arb(5) ** m, arb(5) ** (m + 1), fc.TOL)
    t2 = arb(0)
    for m in range(fc.X_EXP, F2_T):
        t2 += fc._integ(F2, arb(5) ** m, arb(5) ** (m + 1), fc.TOL)
    T = arb(5) ** F2_T
    p = A0 + 2
    an2 = 9 * arb(5) ** 4 * T ** (-p) * (T.log() / p + 1 / (p * p))
    return t1.lower() - (t2 + an2).upper()


def small_part(a1):
    """Upper bound of |int_0^eps (-ln t) t^-alpha (h_1 - h_2)| <= 2 int_0^eps |ln t| t^(1-a1)/30 = 2 eps^p (|ln eps| + 1/p)/(30 p), p = 2 - a1."""
    from flint import arb
    eps = arb(5) ** fc.EPS_EXP
    p = 2 - fc.q(a1)
    return (2 * eps ** p * (-eps.log() + 1 / p) / (30 * p)).upper()


def dtask(args):
    a0, a1 = args
    fc.arb_env()
    t0 = time.time()
    val = (dcore(a0, a1) + dtail_lower(a0, a1) - small_part(a1)).lower()
    return a0, a1, val.str(30, radius=True), time.time() - t0


def main() -> int:
    import multiprocessing as mp
    jobs = 6
    for a in sys.argv[1:]:
        if a.startswith("--jobs="):
            jobs = int(a.split("=")[1])
    t0 = time.time()
    print("Theorem FG: G_2 - G_1 has exactly one zero on (-1,1), in (-0.5679,-0.5677). Functions on a graph, not a flux, not thrust.", flush=True)
    fails = []

    def check(ok, msg):
        print(("  ok   " if ok else "  FAIL ") + msg, flush=True)
        if not ok:
            fails.append(msg)

    check(all(fc.small_t_ok(n) for n in (1, 2)), "FB.1 small-t constants (n = 1, 2): |h_n| <= t/30 at t <= 5^-24")
    check(D_LO < B_LO < B_HI < D_HI == P_LO < P_HI, "bracket inside the Delta' interval; Delta' and Delta intervals abut at -1/2")
    fc.arb_env()
    from flint import arb
    with mp.Pool(jobs) as pool:
        res = pool.map(ff.task, [("pt", 1, B_LO, B_LO), ("pt", 1, B_HI, B_HI)])
        fc.arb_env()
        pts = {r[2]: tuple(arb(x) for x in r[4]) for r in res}
        lo, hi = pts[B_LO]
        print(f"   Delta({float(B_LO)}) in [{lo.str(10, radius=False)}, {hi.str(10, radius=False)}]", flush=True)
        check(bool(hi < 0), f"FG.1: Delta({float(B_LO)}) <= {hi.str(8, radius=False)} < 0")
        lo, hi = pts[B_HI]
        print(f"   Delta({float(B_HI)}) in [{lo.str(10, radius=False)}, {hi.str(10, radius=False)}]", flush=True)
        check(bool(lo > 0), f"FG.1: Delta({float(B_HI)}) >= {lo.str(8, radius=False)} > 0")

        # FG.2 and FG.3, adaptive bisection, run together
        init_d = [D_LO + (D_HI - D_LO) * k / 6 for k in range(7)]
        init_p = [P_LO + (P_HI - P_LO) * k / 4 for k in range(5)]
        pend_d = list(zip(init_d, init_d[1:]))
        pend_p = [("lo", 1, a, b) for a, b in zip(init_p, init_p[1:])]
        done_d, done_p, depth = [], [], 0
        while (pend_d or pend_p) and depth < 7:
            rd = pool.map_async(dtask, pend_d)
            rp = pool.map_async(ff.task, pend_p)
            nd, np_ = [], []
            for a0, a1, v, dt in rd.get():
                fc.arb_env()
                val = arb(v)
                ok = bool(val > 0)
                print(f"   Delta' >= {val.lower().str(8, radius=False)} on [{float(a0):.6f}, {float(a1):.6f}]  ({'ok' if ok else 'split'}, {dt:.0f}s)", flush=True)
                if ok:
                    done_d.append((a0, a1, val.lower()))
                else:
                    m = (a0 + a1) / 2
                    nd += [(a0, m), (m, a1)]
            for kind, n, a0, a1, out, dt in rp.get():
                fc.arb_env()
                val = arb(out[0])
                ok = bool(val > 0)
                print(f"   Delta >= {val.str(8, radius=False)} on [{float(a0):.6f}, {float(a1):.6f}]  ({'ok' if ok else 'split'}, {dt:.0f}s)", flush=True)
                if ok:
                    done_p.append((a0, a1, val))
                else:
                    m = (a0 + a1) / 2
                    np_ += [(kind, n, a0, m), (kind, n, m, a1)]
            pend_d, pend_p, depth = nd, np_, depth + 1

    def covered(iv, lo, hi):
        x = lo
        for a, b, _ in sorted(iv):
            if a <= x:
                x = max(x, b)
        return x >= hi

    okd = not pend_d and covered(done_d, D_LO, D_HI)
    okp = not pend_p and covered(done_p, P_LO, P_HI)
    md = min(d[2] for d in done_d) if okd else None
    mp_ = min(d[2] for d in done_p) if okp else None
    check(okd, f"FG.2: Delta' > 0 on [{D_LO}, {D_HI}] ({len(done_d)} pieces)" + (f"; Delta' >= {md.str(8, radius=False)}" if okd else ""))
    check(okp, f"FG.3: Delta > 0 on [{P_LO}, {P_HI}] ({len(done_p)} pieces)" + (f"; Delta >= {mp_.str(8, radius=False)}" if okp else ""))
    if not fails:
        check(True, f"FG.4: with FF.2, exactly one zero beta_c of G_2 - G_1 on (-1,1), beta_c in ({float(B_LO)}, {float(B_HI)}), simple")
    print(f"elapsed {time.time() - t0:.1f}s", flush=True)
    if fails:
        print(f"{len(fails)} FAILED")
        return 1
    print("all FG checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
