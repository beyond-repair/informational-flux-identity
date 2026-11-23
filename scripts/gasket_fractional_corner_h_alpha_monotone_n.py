#!/usr/bin/env python3
"""Theorem FF: monotonicity in n of G_n(alpha) = rho^-n H_n(alpha) for 2 <= n <= 10 on real alpha in [-1, 1)
(FC.4/FD/FE's kept limit "G_(n+1) > G_n OBSERVED only"), and KEPT FAILURE FF.2 at n = 1.

Setting (EY/FA/FB/FC): G_n(alpha) = int_0^oo t^-alpha f_n dt, W(alpha)/2 = int t^-alpha f dt (EX.1), h_n = f - f_n, D_n = W/2 - G_n.
  G_(n+1) - G_n = D_n - D_(n+1) = int_0^oo t^-alpha (h_n - h_(n+1)) dt = int t^-alpha (f_(n+1) - f_n) dt.

  FF.1 (THEOREM, certified arb).  For every n in {2,...,10} and every real alpha in [-1, 1):  G_(n+1)(alpha) > G_n(alpha),
        with certified 3^n (G_(n+1) - G_n) >= mu_n on [-1, 1] (printed).  Method = FC.2's: integrand h_n - h_(n+1) from FC.1's
        cancellation-free ball evaluator (both terms have relative width ~2^-100 and their difference is ~(2/3) h_n, so no cancellation);
        min over alpha in [a0,a1] of t^-alpha g is min(t^-a0 g, t^-a1 g); Arb integrator with acb.real_min on [5^m, 5^(m+1)], m = -24..6;
        |int_0^(5^-24)| <= 2 * 5^-24/30 (FB.1 for n and n+1); beyond X = 5^7: f_(n+1) >= 0, so int_X^oo t^-alpha (f_(n+1) - f_n) >= -int_X^oo t^-a0 f_n
        (FC.2's Arb + analytic tail, t > 1, alpha >= a0).
  FF.2 (KEPT FAILURE).  n = 1 is NOT monotone: G_2(alpha) < G_1(alpha) for every alpha in (-1, alpha_c], alpha_c in (-0.7408, -0.7406) (FD),
        because G_2 < W/2 (FC.2) and G_1 >= W/2 there (FC.3, FD.3).  Witness alpha = -4/5: G_2 - G_1 = D_1 - D_2 < D_1(-4/5) <= -3.3948 (FC.3).
        Certified here (arb): G_2 - G_1 < 0 on (-1, N1_NEG] (direct, hi mode), and G_2 - G_1 > 0 on [N1_POS, 1).
        The sign of G_2 - G_1 on (-13/20, -3/8) is OPEN (float crossing ~ -0.56, OBSERVED; the alpha-interval min loses too much
        near the crossing because of f_1's slow tail, so [-1/2, -3/8] and [-11/20, -1/2] did not certify at one bisection level).
Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, __file__.rsplit("/", 1)[0])

import gasket_fractional_corner_h_alpha_below_small_n as fc  # noqa: E402

N_POS = range(2, 11)
N1_NEG = Fr(-13, 20)   # G_2 - G_1 < 0 on (-1, -13/20]
N1_POS = Fr(-3, 8)    # G_2 - G_1 > 0 on [-3/8, 1)
WITNESS = Fr(-4, 5)


def dG_float(n, a, N=3000):
    """Float (non-rigorous) G_(n+1) - G_n = D_n - D_(n+1) (CI sanity)."""
    return fc.D_float(n, a, N=N) - fc.D_float(n + 1, a, N=N)


def core(n, a0, a1, mode, tol=fc.TOL):
    """int_(5^-24)^(5^7) of min ('lo') / max ('hi') over alpha in [a0,a1] of t^-alpha (h_n - h_(n+1)), or at alpha = a0 ('pt')."""
    from flint import arb
    A0, A1 = fc.q(a0), fc.q(a1)

    def F(z, analytic):
        g = fc.h_eval(n, z)[0] - fc.h_eval(n + 1, z)[0]
        lz = z.log()
        u = g * (-A0 * lz).exp()
        if mode == "pt":
            return u
        v = g * (-A1 * lz).exp()
        return u.real_min(v, analytic=analytic) if mode == "lo" else u.real_max(v, analytic=analytic)

    tot = arb(0)
    for m in range(fc.EPS_EXP, fc.X_EXP):
        tot += fc._integ(F, arb(5) ** m, arb(5) ** (m + 1), tol)
    return tot


def task(args):
    kind, n, a0, a1 = args
    fc.arb_env()
    t0 = time.time()
    eps2 = 2 * fc.eps_part()
    if kind == "lo":
        val = (core(n, a0, a1, "lo") - eps2 - fc.tail_fn_arb(n, a0)).lower()
        out = (val, None)
    elif kind == "hi":
        # beyond X (t > 1, alpha in [a0,a1]): t^-alpha (f_(n+1) - f_n) <= t^-a0 f_(n+1) - t^-a1 f_n
        val = (core(n, a0, a1, "hi") + eps2 + fc.tail_fn_arb(n + 1, a0) - fc.tail_fn_arb(n, a1, lower=True)).upper()
        out = (None, val)
    else:
        c = core(n, a0, a0, "pt")
        out = ((c - eps2 - fc.tail_fn_arb(n, a0)).lower(),
               (c + eps2 + fc.tail_fn_arb(n + 1, a0) - fc.tail_fn_arb(n, a0, lower=True)).upper())
    out = tuple(None if x is None else x.str(30, radius=True) for x in out)
    return kind, n, a0, a1, out, time.time() - t0


def main() -> int:
    import multiprocessing as mp
    jobs, nlist = 6, list(N_POS)
    do_n1 = True
    for a in sys.argv[1:]:
        if a.startswith("--jobs="):
            jobs = int(a.split("=")[1])
        elif a.startswith("--n="):
            nlist = [int(x) for x in a.split("=")[1].split(",") if x]
        elif a == "--no-n1":
            do_n1 = False
    t0 = time.time()
    print("Theorem FF: G_(n+1) > G_n for n = 2..10 on real alpha in [-1,1); kept failure at n = 1. Functions on a graph, not a flux, not thrust.", flush=True)
    fails = []

    def check(ok, msg):
        print(("  ok   " if ok else "  FAIL ") + msg, flush=True)
        if not ok:
            fails.append(msg)

    check(all(fc.small_t_ok(n) for n in range(1, 12)), "FF.1: small-t part (FB.1, n = 1..11): |h_n|, |h_(n+1)| <= t/30 at t <= 5^-24, so |int_0^eps t^-alpha (h_n - h_(n+1))| <= 2 eps/30")
    fc.arb_env()
    from flint import arb
    import gasket_fractional_corner_h_alpha_rate as ey
    okp = True
    for n in (1, 2, 10):
        for t in (1e-3, 0.5, 3.0, 40.0, 900.0):
            g = fc.h_eval(n, arb(t))[0] - fc.h_eval(n + 1, arb(t))[0]
            okp &= abs(float(g.real.mid()) - (ey.f_n(n + 1, t) - ey.f_n(n, t))) <= 1e-9 * abs(ey.f_inf(t)) + 1e-15
    check(okp, "FC.1 evaluator: h_n - h_(n+1) agrees with float f_(n+1) - f_n (n = 1,2,10; 5 points)")

    init = [Fr(k, 4) for k in range(-4, 5)]
    pending = [("lo", n, a, b) for n in nlist for a, b in zip(init, init[1:])]
    if do_n1:
        n1p = [N1_POS, Fr(-1, 4), Fr(0), Fr(1, 2), Fr(1)]
        pending += [("lo", 1, a, b) for a, b in zip(n1p, n1p[1:])]
        n1n = [Fr(-1), Fr(-9, 10), Fr(-4, 5), Fr(-3, 4), Fr(-7, 10), N1_NEG]
        pending += [("hi", 1, a, b) for a, b in zip(n1n, n1n[1:])]
        pending.append(("pt", 1, WITNESS, WITNESS))
    done, depth = [], 0
    with mp.Pool(jobs) as pool:
        while pending and depth < 9:
            nxt = []
            for kind, n, a0, a1, out, dt in pool.imap_unordered(task, pending):
                lo, hi = (None if x is None else arb(x) for x in out)
                if kind == "lo":
                    ok = bool(lo > 0)
                    print(f"   n={n} [{a0}, {a1}]: 3^n (G_(n+1) - G_n) >= {(lo * 3 ** n).str(6, radius=False)}  ({'ok' if ok else 'split'}, {dt:.0f}s, t={time.time() - t0:.0f}s)", flush=True)
                elif kind == "hi":
                    ok = bool(hi < 0)
                    print(f"   n=1 ({a0}, {a1}]: 3 (G_2 - G_1) <= {(hi * 3).str(6, radius=False)}  ({'ok' if ok else 'split'}, {dt:.0f}s, t={time.time() - t0:.0f}s)", flush=True)
                else:
                    ok = bool(hi < 0)
                    print(f"   witness n=1, alpha={a0}: G_2 - G_1 in [{lo.str(10, radius=False)}, {hi.str(10, radius=False)}]  ({dt:.0f}s)", flush=True)
                if ok:
                    done.append((kind, n, a0, a1, lo, hi))
                elif kind == "pt":
                    fails.append("witness")
                else:
                    m = (a0 + a1) / 2
                    nxt += [(kind, n, a0, m), (kind, n, m, a1)]
            pending, depth = nxt, depth + 1
    check(not pending, f"all intervals certified (adaptive depth {depth})")

    def covered(kind, n, lo, hi):
        iv = sorted((d[2], d[3]) for d in done if d[0] == kind and d[1] == n)
        x = lo
        for a, b in iv:
            if a <= x:
                x = max(x, b)
        return x >= hi

    for n in nlist:
        ok = covered("lo", n, Fr(-1), Fr(1))
        mn = min((d[4] * 3 ** n) for d in done if d[0] == "lo" and d[1] == n) if ok else None
        check(ok, f"FF.1: n = {n}: G_{n + 1} > G_{n} for every real alpha in [-1, 1)" + (f"; 3^n (G_(n+1) - G_n) >= {mn.str(5, radius=False)}" if ok else ""))
    if do_n1:
        check(covered("lo", 1, N1_POS, Fr(1)), f"FF.2: n = 1: G_2 > G_1 on [{N1_POS}, 1)")
        check(covered("hi", 1, Fr(-1), N1_NEG), f"FF.2 KEPT FAILURE: n = 1: G_2 < G_1 on (-1, {N1_NEG}]")
        check(any(d[0] == "pt" for d in done), f"FF.2 witness: G_2({WITNESS}) - G_1({WITNESS}) < 0 certified")
    print(f"elapsed {time.time() - t0:.1f}s", flush=True)
    if fails:
        print(f"{len(fails)} FAILED")
        return 1
    print("all FF checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
