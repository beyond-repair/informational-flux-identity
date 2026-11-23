#!/usr/bin/env python3
"""Theorem FC: approach from below for EX.1's limit at small n (FB's kept limit "n <= 10 from below OBSERVED only"), and KEPT FAILURE FC.3 at n = 1.

Setting (EY/FA/FB): G_n(alpha) = rho^-n H_n(alpha) = int_0^oo t^-alpha f_n dt, W(alpha)/2 = int t^-alpha f dt (EX.1, -2^(n-1) < Re alpha < 1),
f_n = R_n/(2(5+psi_(n-1)(t/5))), f = R/(2(5+psi(t/5))), so D_n(alpha) := W(alpha)/2 - G_n(alpha) = int_0^oo t^-alpha h_n dt, h_n = f - f_n.

  FC.1 (THEOREM; cancellation-free ball evaluator of h_n).  With s = 5^-n t, Z_k = p^k(s), Y_k = psi(5^k s) = p^k(psi(s)), d_k = Y_k - Z_k:
        d_0 = psi(s) - s (EV.1 series, from 5^-L s by e -> psi^2 + 5 e), d_(k+1) = d_k (Y_k + Z_k + 5),
        log(P/P_n) = sum_(k<n) log1p(2 d_k/((Y_k+5)(Z_k+3))) + sum_(k<0) log g(Y_k) + [EV.2 truncation ball, radius 2 delta],
        Delta = log(P/P_n) + log1p(d_(n-1)/(5+Z_(n-1))),  E = (t/6) 3^-n (1/P_n + (3/5)^n),  h_n = f (E - expm1(Delta))/(1+E)  (FB.1 split).
  FC.2 (THEOREM, certified arb).  For every n in {2,...,10} and every real alpha in [-1, 1):  D_n(alpha) = W(alpha)/2 - rho^-n H_n(alpha) > 0,
        with certified 3^n D_n >= m_n on [-1, 1] (printed).  With FB.4 (n >= 11): rho^-n H_n(alpha) < W(alpha)/2 for every n >= 2 and every real
        alpha in [-1, 1).  Method: min over alpha in [a0,a1] of t^-alpha h is min(t^-a0 h, t^-a1 h) (monotone in alpha), integrated by Arb's
        rigorous integrator with acb.real_min on blocks [5^m, 5^(m+1)], m = -24..6; |int_0^(5^-24)| <= 5^-24/30 (FB.1: |h| <= f max(E, e^Lam - 1) <= t/30);
        beyond X = 5^7 drop f >= 0 and subtract int_X^T t^-a0 f_n (Arb) + analytic tail (f_n <= 3^n/(t psi_(n-1)(t/5)), psi_(n-1)(t/5) >= (t/5^(J+1))^(2^J)).
  FC.3 (KEPT FAILURE, certified arb).  n = 1 is NOT from below near alpha = -1: D_1(alpha) < 0 for every alpha in (-1, -3/4], certified witness
        D_1(-4/5) < 0 (enclosure printed); and D_1(alpha) > 0 for every alpha in [-73/100, 1).  (G_1 -> +oo as alpha -> -1^+ since f_1 ~ 15/t^2.)
        Float crossing alpha_c = -0.74070 (OBSERVED, not certified beyond the bracket (-3/4, -73/100)).
Functions on a graph, not a flux.  Nothing here is thrust.  No Psi_info; Stage 2 still blocked.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, __file__.rsplit("/", 1)[0])

EPS_EXP = -24   # small-t cut 5^-24
X_EXP = 7       # h integrated up to 5^7
T_EXP = {1: 60}  # f_n integrated to 5^T (n=1), else 5^20
T_DEFAULT = 20
TOL = 26
N_POS = range(2, 11)
N1_NEG = Fr(-3, 4)
N1_POS = Fr(-73, 100)
WITNESS = Fr(-4, 5)


# ---------------- pure-Python float pieces (CI) ----------------
def D_float(n, a, lo=1e-14, hi_exp=13, N=3000):
    """Simpson in log t of int t^-a (f - f_n) (non-rigorous, CI sanity)."""
    import gasket_fractional_corner_h_alpha_rate as ey
    x0, x1 = math.log(lo), hi_exp * math.log(5)
    h = (x1 - x0) / N
    s = 0.0
    for i in range(N + 1):
        t = math.exp(x0 + i * h)
        fi = ey.f_inf(t) if t < 5.0 ** 8 else 0.0
        w = 1 if i in (0, N) else (4 if i % 2 else 2)
        s += w * t ** (1 - a) * (fi - ey.f_n(n, t))
    return s * h / 3


def small_t_ok(n, eps=5.0 ** EPS_EXP):
    """FB.1 constants at t <= eps: E_up <= t and e^Lam - 1 <= 2 Lam <= t (so |h| <= f t <= t/30)."""
    cs, kap = 4 - math.sqrt(15), math.exp(math.exp(0.1) / 4)
    Eup = eps / 6 * 3.0 ** -n * (1 + 0.6 ** n)
    w = eps / 5.0 ** n
    Lam = w * ((cs * n + 1) * kap / 10 + math.exp(w / 50) / 30)
    return Eup <= eps and Lam <= 0.5 and 2 * Lam <= eps


def tail_fn_analytic(n, a0, Texp):
    """Upper bound of int_(5^Texp)^oo t^-a0 f_n dt (float form, exact arb in tail_fn_arb)."""
    if n == 1:
        return 15 * 5.0 ** (Texp * (-1 - a0)) / (1 + a0)
    J = min(n - 1, 3)
    return 3.0 ** n * 5.0 ** ((J + 1) * 2 ** J) * 5.0 ** (Texp * (-a0 - 2 ** J)) / (2 ** J + a0)


# ---------------- arb ----------------
def arb_env(prec=128):
    from flint import acb, arb, ctx
    ctx.prec = prec
    return acb, arb


def q(x):
    from flint import arb
    x = Fr(x)
    return arb(x.numerator) / x.denominator


def psi_small_e(v):
    from flint import acb, arb
    a = abs(v)
    assert a.upper() < 1
    err = (a ** 3 / 400 / (1 - a / 20)).upper()
    return v * v / 20 + acb(arb(0, err), arb(0, err))


def h_eval(n, t):
    """FC.1: (h_n, f) as complex balls, h_n = f - f_n without cancellation."""
    from flint import acb, arb
    five = arb(5)
    t = acb(t)
    s = t / five ** n
    a = abs(s).upper()
    L = 0
    while a / five ** L > arb(2) ** -90 or L < 26:
        L += 1
    v = s / five ** L
    e = psi_small_e(v)
    y = v + e
    aw = abs(v).upper()
    delta = (arb(2) * arb(101) / 100 * aw / 4 / (3 * (5 - arb(101) / 100 * aw))).upper()
    # EV.2: |f_j - 1| = |eps_j| with sum_(j>J) |eps_j| <= delta; |log(1+x)| <= 2|x| for |x| <= 1/2, so |sum log f_j| <= 2 delta
    assert bool(delta < arb(1) / 4)
    logext = acb(arb(0, 2 * delta), arb(0, 2 * delta))
    for _ in range(L):
        logext += (2 * y / (3 * (y + 5))).log1p()
        e = y * y + 5 * e
        y = y * (y + 5)
    Z, Y, d = s, y, e
    logratio, logP, Pn = logext, logext, acb(1)
    for k in range(n):
        logratio += (2 * d / ((Y + 5) * (Z + 3))).log1p()
        Pn = Pn * (five / 3) * (Z + 3) / (Z + 5)
        logP += (2 * Y / (3 * (Y + 5))).log1p()
        if k < n - 1:
            d = d * (Y + Z + 5)
            Y = Y * (Y + 5)
            Z = Z * (Z + 5)
    Delta = logratio + (d / (5 + Z)).log1p()
    E = t / 6 / arb(3) ** n * (1 / Pn + (arb(3) / 5) ** n)
    f = 1 / (3 * logP.exp()) / (2 * (5 + Y))
    return f * (E - Delta.expm1()) / (1 + E), f


def fn_eval(n, t):
    from flint import acb, arb
    five = arb(5)
    t = acb(t)
    Z, Pn = t / five ** n, acb(1)
    for k in range(n):
        Pn = Pn * (five / 3) * (Z + 3) / (Z + 5)
        if k < n - 1:
            Z = Z * (Z + 5)
    invR = 3 * Pn + t / 2 / arb(3) ** n * (1 + (arb(3) / 5) ** n * Pn)
    return 1 / (invR * 2 * (5 + Z))


def _integ(F, lo, hi, tol):
    from flint import acb, arb
    return acb.integral(F, lo, hi, rel_tol=arb(2) ** -tol, abs_tol=arb(2) ** -(tol + 12)).real


def core(n, a0, a1, mode, tol=TOL):
    """int_(5^-24)^(5^7) of min (mode 'lo'), max ('hi') over alpha in [a0,a1] of t^-alpha h_n, or t^-a0 h_n ('pt', a0 = a1)."""
    from flint import arb
    A0, A1 = q(a0), q(a1)

    def F(z, analytic):
        h, _ = h_eval(n, z)
        lz = z.log()
        u = h * (-A0 * lz).exp()
        if mode == "pt":
            return u
        v = h * (-A1 * lz).exp()
        return u.real_min(v, analytic=analytic) if mode == "lo" else u.real_max(v, analytic=analytic)

    tot = arb(0)
    for m in range(EPS_EXP, X_EXP):
        tot += _integ(F, arb(5) ** m, arb(5) ** (m + 1), tol)
    return tot


def tail_fn_arb(n, a0, tol=TOL, lower=False):
    """Upper bound of int_(5^7)^oo t^-a0 f_n: Arb on [5^7, 5^T] plus FC.2's analytic tail (lower=True: Arb lower bound of int_(5^7)^(5^T), f_n >= 0)."""
    from flint import arb
    A0 = q(a0)
    Texp = T_EXP.get(n, T_DEFAULT)

    def F(z, analytic):
        return fn_eval(n, z) * (-A0 * z.log()).exp()

    tot = arb(0)
    for m in range(X_EXP, Texp):
        tot += _integ(F, arb(5) ** m, arb(5) ** (m + 1), tol)
    if lower:
        return tot.lower()
    T = arb(5) ** Texp
    if n == 1:
        assert a0 > -1
        an = 15 * T ** (-1 - A0) / (1 + A0)
    else:
        J = min(n - 1, 3)
        an = arb(3) ** n * arb(5) ** ((J + 1) * 2 ** J) * T ** (-A0 - 2 ** J) / (2 ** J + A0)
    return (tot + an).upper()


def f_tail_upper():
    """int_(5^7)^oo t f <= 5^80 X^-14/84 (f <= 1/(6 psi_4(t/5)) <= (5^5/t)^16/6)."""
    from flint import arb
    return (arb(5) ** 80 * (arb(5) ** X_EXP) ** -14 / 84).upper()


def eps_part():
    from flint import arb
    return (arb(5) ** EPS_EXP / 30).upper()


def task(args):
    kind, n, a0, a1 = args
    arb_env()
    t0 = time.time()
    if kind == "lo":
        val = (core(n, a0, a1, "lo") - eps_part() - tail_fn_arb(n, a0)).lower()
        out = (val, None)
    elif kind == "hi":
        # beyond X: t^-alpha h <= t^-alpha f - t^-a1 f_n  (t > 1, alpha <= a1), so subtract a lower bound of int_X^(5^T) t^-a1 f_n
        val = (core(n, a0, a1, "hi") + eps_part() + f_tail_upper() - tail_fn_arb(n, a1, lower=True)).upper()
        out = (None, val)
    else:
        c = core(n, a0, a0, "pt")
        out = ((c - eps_part() - tail_fn_arb(n, a0)).lower(), (c + eps_part() + f_tail_upper() - tail_fn_arb(n, a0, lower=True)).upper())
    # arb balls are not picklable: ship rigorous decimal enclosures (arb_printn) and re-parse them (arb_set_str)
    out = tuple(None if x is None else x.str(30, radius=True) for x in out)
    return kind, n, a0, a1, out, time.time() - t0


def main() -> int:
    import multiprocessing as mp
    jobs = 6
    for a in sys.argv[1:]:
        if a.startswith("--jobs="):
            jobs = int(a.split("=")[1])
    t0 = time.time()
    print("Theorem FC: approach from below for n <= 10 (n = 2..10), kept failure at n = 1. Functions on a graph, not a flux, not thrust.", flush=True)
    fails = []

    def check(ok, msg):
        print(("  ok   " if ok else "  FAIL ") + msg, flush=True)
        if not ok:
            fails.append(msg)

    check(all(small_t_ok(n) for n in range(1, 11)), "FC.2: small-t part: E_up <= t, e^Lam - 1 <= 2 Lam <= t at t <= 5^-24 (FB.1), so |int_0^eps t^-alpha h| <= eps/30")
    arb_env()
    from flint import arb
    # pointwise cross-check of FC.1 evaluator vs direct f - f_n
    import gasket_fractional_corner_h_alpha_rate as ey
    okp = True
    for n in (1, 2, 6, 10):
        for t in (1e-3, 0.5, 3.0, 40.0, 900.0):
            h, _ = h_eval(n, arb(t))
            okp &= abs(float(h.real.mid()) - (ey.f_inf(t) - ey.f_n(n, t))) <= 1e-9 * abs(ey.f_inf(t)) + 1e-15
    check(okp, "FC.1: ball evaluator agrees with float f - f_n (n = 1,2,6,10; 5 points)")
    check(bool(f_tail_upper() < arb(10) ** -14) and bool(eps_part() < arb(10) ** -17), "tails: int_X^oo t f <= 5^80 X^-14/84 < 1e-14, eps/30 < 1e-17")

    # task list
    init_pos = [Fr(k, 4) for k in range(-4, 5)]
    pending = [("lo", n, a, b) for n in N_POS for a, b in zip(init_pos, init_pos[1:])]
    n1p = [N1_POS, Fr(-18, 25), Fr(-7, 10), Fr(-13, 20), Fr(-3, 5), Fr(-1, 2), Fr(-1, 4), Fr(0), Fr(1, 2), Fr(1)]
    pending += [("lo", 1, a, b) for a, b in zip(n1p, n1p[1:])]
    n1n = [Fr(-1), Fr(-9, 10), Fr(-17, 20), Fr(-4, 5), Fr(-31, 40), Fr(-19, 25), N1_NEG]
    pending += [("hi", 1, a, b) for a, b in zip(n1n, n1n[1:])]
    pending.append(("pt", 1, WITNESS, WITNESS))
    done = []
    depth = 0
    with mp.Pool(jobs) as pool:
        while pending and depth < 9:
            nxt = []
            for kind, n, a0, a1, out, dt in pool.imap_unordered(task, pending):
                lo, hi = (None if x is None else arb(x) for x in out)
                if kind == "lo":
                    ok = bool(lo > 0)
                    print(f"   n={n} [{a0}, {a1}]: 3^n D_n >= {(lo * 3 ** n).str(6, radius=False)}  ({'ok' if ok else 'split'}, {dt:.0f}s, t={time.time() - t0:.0f}s)", flush=True)
                elif kind == "hi":
                    ok = bool(hi < 0)
                    print(f"   n=1 ({a0}, {a1}]: 3 D_1 <= {(hi * 3).str(6, radius=False)}  ({'ok' if ok else 'split'}, {dt:.0f}s, t={time.time() - t0:.0f}s)", flush=True)
                else:
                    ok = bool(hi < 0)
                    print(f"   witness n=1, alpha={a0}: D_1 in [{lo.str(10, radius=False)}, {hi.str(10, radius=False)}]  ({dt:.0f}s)", flush=True)
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

    for n in N_POS:
        ok = covered("lo", n, Fr(-1), Fr(1))
        mn = min((d[4] * 3 ** n) for d in done if d[0] == "lo" and d[1] == n) if ok else None
        check(ok, f"FC.2: n = {n}: W/2 - rho^-n H_n > 0 for every real alpha in [-1, 1)" + (f"; 3^n D_n >= {mn.str(5, radius=False)}" if ok else ""))
    check(covered("lo", 1, N1_POS, Fr(1)), f"FC.3: n = 1: D_1 > 0 on [{N1_POS}, 1)")
    check(covered("hi", 1, Fr(-1), N1_NEG), f"FC.3 KEPT FAILURE: n = 1: D_1 < 0 on (-1, {N1_NEG}] (G_1 above W/2)")
    check(any(d[0] == "pt" for d in done), f"FC.3 witness: D_1({WITNESS}) < 0 certified")
    # float cross-checks
    check(D_float(1, -0.8) < 0 < D_float(1, -0.73) and min(D_float(n, a) for n in (2, 10) for a in (-1.0, 0.0, 0.999)) > 0,
          "float cross-check: D_1(-0.8) < 0 < D_1(-0.73); D_2, D_10 > 0 at alpha = -1, 0, 0.999")
    print(f"elapsed {time.time() - t0:.1f}s", flush=True)
    if fails:
        print(f"{len(fails)} FAILED")
        return 1
    print("all FC checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
