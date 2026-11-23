"""Theorem FD: the n = 1 crossing of D_1(alpha) = W/2 - rho^-1 H_1 is unique on [-1,1) and lies in [-0.7408, -0.7406].

Builds on FC (same evaluator h_eval, same Arb tails). D_1(a) = int_0^oo t^-a h_1 dt, h_1 = f - f_1 (EX.1, a > -1).
FD.1 (arb): D_1(A_LO) < 0 < D_1(A_HI), A_LO=-0.7408, A_HI=-0.7406.
FD.2 (arb): D_1'(a) = int t^-a (-ln t) h_1 dt >= m > 0 for all a in [-3/4, -73/100]
    (differentiation under the integral: |t^-a ln t h_1| integrable locally uniformly, h_1 = O(t) at 0, f_1 ~ 15/t^2... at oo for a > -1).
With FC.3 (D_1 < 0 on (-1,-3/4], D_1 > 0 on [-73/100,1)): exactly one zero of D_1 on (-1,1), in [A_LO, A_HI].
Functions on a graph, not a flux, not thrust.
"""
from __future__ import annotations
import sys, time
from fractions import Fraction as Fr
sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_h_alpha_below_small_n as fc

A_LO, A_HI = Fr(-7408, 10000), Fr(-7406, 10000)
G0, G1 = Fr(-3, 4), Fr(-73, 100)


def dcore(a0, a1):
    """min over alpha in [a0,a1] of int_(5^-24)^(5^7) t^-alpha (-ln t) h_1 dt (integrand monotone in alpha)."""
    from flint import arb
    A0, A1 = fc.q(a0), fc.q(a1)

    def F(z, analytic):
        h, _ = fc.h_eval(1, z)
        lz = z.log()
        base = h * (-lz)
        u = base * (-A0 * lz).exp()
        v = base * (-A1 * lz).exp()
        return u.real_min(v, analytic=analytic)
    tot = arb(0)
    for m in range(fc.EPS_EXP, fc.X_EXP):
        tot += fc._integ(F, arb(5) ** m, arb(5) ** (m + 1), fc.TOL)
    return tot


def dtail_lower(a1):
    """lower bound of int_(5^7)^oo t^-alpha ln t (f_1 - f) dt, alpha<=a1: Arb on [5^7,5^60] of ln t t^-a1 f_1 minus f tail (alpha>=-3/4)."""
    from flint import arb
    A1 = fc.q(a1)

    def F(z, analytic):
        return fc.fn_eval(1, z) * (-A1 * z.log()).exp() * z.log()
    tot = arb(0)
    for m in range(fc.X_EXP, fc.T_EXP[1]):
        tot += fc._integ(F, arb(5) ** m, arb(5) ** (m + 1), fc.TOL)
    X = arb(5) ** fc.X_EXP
    # f <= (5^5/t)^16/6, ln t <= 4 t^0.1, t^-alpha <= t^0.75:  int_X^oo <= 4*5^80/6 * X^-14.15/14.15
    ftail = 4 * arb(5) ** 80 / 6 * X ** (arb(-1415) / 100) / (arb(1415) / 100)
    return (tot.lower() - ftail.upper())


def small_part(a0):
    """|int_0^eps t^-alpha ln t h_1| <= int_0^eps |ln t| t^(1-alpha)/30 <= eps^2.73 (|ln eps| + 1/2.73)/(2.73 * 30) (alpha >= -3/4 >= -0.73-.02)."""
    from flint import arb
    eps = arb(5) ** fc.EPS_EXP
    return ((eps ** (arb(273) / 100) * (-eps.log() + arb(100) / 273) / (arb(273) / 100 * 30))).upper()


def dtask(args):
    a0, a1 = args
    fc.arb_env()
    t0 = time.time()
    val = (dcore(a0, a1) + dtail_lower(a1) - small_part(a0)).lower()
    return a0, a1, val.str(20, radius=False), time.time() - t0


def main() -> int:
    import multiprocessing as mp
    t0 = time.time()
    fails = []

    def check(ok, msg):
        print(("PASS " if ok else "FAIL ") + msg, flush=True)
        if not ok:
            fails.append(msg)
    print("Theorem FD: n = 1 crossing of D_1 is unique on [-1,1), in [-0.7408,-0.7406]. Functions on a graph, not a flux, not thrust.", flush=True)
    jobs = [("lo", A_HI), ("hi", A_LO)]
    with mp.Pool(6) as pool:
        r = pool.map(fc.task, [("pt", 1, A_HI, A_HI), ("pt", 1, A_LO, A_LO)])
        fc.arb_env()
        from flint import arb
        lo_hi = arb(r[0][4][0], ) if False else None
        pts = {}
        for kind, n, a0, a1, out, dt in r:
            pts[a0] = (arb(out[0]) if out[0] else None, arb(out[1]) if out[1] else None)
        pHI, pLO = pts[A_HI], pts[A_LO]
        check(bool(pHI[0] > 0), f"FD.1: D_1({A_HI}) >= {pHI[0].str(8, radius=False)} > 0")
        check(bool(pLO[1] < 0), f"FD.1: D_1({A_LO}) <= {pLO[1].str(8, radius=False)} < 0")
        # FD.2 adaptive
        pending, done, depth = [(G0, G1)], [], 0
        while pending and depth < 8:
            res = pool.map(dtask, pending)
            nxt = []
            for a0, a1, v, dt in res:
                fc.arb_env()
                val = arb(v)
                print(f"   D_1' >= {v} on [{float(a0):.6f},{float(a1):.6f}] ({dt:.0f}s)", flush=True)
                if bool(val > 0):
                    done.append((a0, a1, val))
                else:
                    m = (a0 + a1) / 2
                    nxt += [(a0, m), (m, a1)]
            pending, depth = nxt, depth + 1
        check(not pending, f"FD.2: D_1' > 0 on all of [{G0},{G1}] (adaptive depth {depth}, {len(done)} pieces)")
        if done:
            mn = min(d[2] for d in done)
            check(True, f"FD.2: min over pieces of certified D_1' lower bound = {mn.str(6, radius=False)}")
    print(f"elapsed {time.time() - t0:.1f}s")
    if fails:
        print(f"{len(fails)} FAILED"); return 1
    print("all FD checks passed (with FC.3: unique zero of D_1 on (-1,1), in [A_LO, A_HI])")
    return 0


if __name__ == "__main__":
    sys.exit(main())
