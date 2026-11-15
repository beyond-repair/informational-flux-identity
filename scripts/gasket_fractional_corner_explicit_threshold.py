#!/usr/bin/env python3
"""Theorem EG: an explicit threshold.  EF's kept limit "no explicit n_0 (EE's limits carry no rate)" is removed by a one-sided,
every-n version of EE.1: for every alpha in (gamma,1) and every k0 >= 1 there is an explicit c_lo(k0) (an integral of level-k0
and level-oo profile objects) with
    c_k = (3/rho)^k (Mt_(k+1) - Mt_k) >= c_lo(k0)   for EVERY k >= k0,       c_lo(k0) -> d*  (k0 -> oo),
hence (3/rho)^n eps_n >= L(k0) = 3 c_lo(k0)/(2(1-rho/3)) and 3^n eta_n >= (3/2) kappa_lo(N,k0) for every n >= max(N,k0).

  EG.1 (THEOREM; every alpha in (gamma,1), k >= k0 >= 1).  Pointwise in x = 5^k t:
       Psi_k >= x^-a 5^k b (Rh_k - R_k)(R_k + Rh_k)       (omega^0 >= 1, 1/((1-u)(1-v)) >= 1+u+v),
       Rh_k - R_k = R_k Rh_k A_k (phi_k+6)/(phi_k+5)       (EE.1's product form),  A_k >= A,  phi_k <= Phi,
       1/R_k <= V_k0(x) = A_k0(x)(3 + Phi(x)) - sum_(m=1..k0) e(phi_(k0-m)(5^-m x)) (5/3)^m P_m(x)   (phi_j nondecreasing in j),
       5^k b >= 2x/(6 + 5^(-k0-1) x).  So c_k >= c_lo(k0) = int x^(1-a) 2/(6+5^(-k0-1)x) A (Phi+6)/(Phi+5) Rl Rhl (Rl+Rhl),
       Rl = 1/V_k0(x), Rhl = (5/3)/V_(k0+1)(5x); and the same for eta: kappa_lo(N,k0) = int x^(1-a) 2/(6+5^(-k0-1)x) Rl^2 eta_N(Phi(5^-N x)).
  EG.2 (THEOREM, exact algebra on EB.1, n >= N): at rho^3 = 3,
       rho^n (Y_n - Y_oo) <= (1-delta_N)(|nu|^3 - L(k0)/D'_oo)/(1 + nu rho^-n) - (3/2) kappa_lo(N,k0)/D'_oo   (n >= max(N,k0)),
       and for rho^3 >= 3, |nu|^3 < L(k0)/D'_oo gives Y_n < Y_oo for every n >= max(N,k0).
  EG.3 (arb at alpha = 4gamma/3 (ball)): L(2) > |nu|^3 D'_oo, and at N = n = 1 the bound with k0 = 1 is < 0: Y_n < Y_oo for EVERY
       N >= 1 and EVERY n >= N (explicit n_0 = N).  Every-n brackets: rho^n (Y_n - Y_oo) <= S_N^+ < 0 for n >= max(N, 6).
  EG.4 (arb at 0.7, 0.8, 0.9, 0.95): L(k0), every-n lower bounds for 3^n R_n at 0.8 (ED.4's lower column sharpened); at 0.95
       Y_n < Y_oo for every n >= N.
  EG.5 (OBSERVED, floats): c_lo(k0) sits below the direct quadrature of c_k (EE.5's routine) and approaches d* like ~3^-k0.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

FAILS: list[str] = []
RESULT: dict = {}
GAMMA = math.log(3) / math.log(5)
NS = (1, 2, 3, 6, 10)


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
def _be(s):
    return (s + 6) / ((s + 2) * (s + 5))


def _e(s):
    return 2 * s / (2 + s)


def vplus_float(x, k0: int, M: int = 50):
    """Floats: V_k0(x) = A_k0(x)(3+Phi(x)) - sum_(m<=k0) e(phi_(k0-m)(5^-m x)) (5/3)^m P_m(x), plus A(x), Phi(x), orbit z_m = Phi(5^-m x)."""
    x = np.asarray(x, float)
    zs = [x * 5.0 ** (-M)]
    for _ in range(M):
        zs.append(np.minimum(zs[-1] * (zs[-1] + 5), 1e150))
    zs = zs[::-1]
    P = [None] * (M + 1)
    P[M] = np.ones_like(x)
    for m in range(M - 1, -1, -1):
        P[m] = P[m + 1] * (5 / 3) * _be(zs[m + 1])
    s = [x * 5.0 ** (-k0)]
    for _ in range(k0):
        s.append(np.minimum(s[-1] * (s[-1] + 5), 1e150))
    Ak = np.ones_like(x)
    for m in range(1, k0 + 1):
        Ak = Ak * (5 / 3) * _be(s[k0 - m])
    V = Ak * (3 + zs[0])
    for m in range(1, k0 + 1):
        V = V - _e(s[k0 - m]) * (5 / 3) ** m * P[m]
    return V, P[0], zs[0], zs


def lower_float(al: float, k0: int, Ns=(), n: int = 200, lo: int = -25, hi: int = 8):
    """Floats (trapezoid in log x) of c_lo(k0), L(k0) = 3 c_lo/(2(1-rho/3)) and kappa_lo(N,k0)."""
    j = np.arange(lo * n, hi * n + 1)
    x = 5.0 ** (j / n)
    V, A, Phi, _ = vplus_float(x, k0)
    V1 = vplus_float(5 * x, k0 + 1)[0]
    R, Rh = 1 / V, (5 / 3) / V1
    wb = 2 / (6 + 5.0 ** (-k0 - 1) * x)
    c = float(np.trapezoid(x ** (2 - al) * wb * A * (Phi + 6) / (Phi + 5) * R * Rh * (R + Rh), np.log(x)))
    rho = 5**al / 3
    out = {"c": c, "L": 1.5 * c / (1 - rho / 3), "kappa": {}}
    for N in Ns:
        jj = np.arange(lo * n, (N + 6) * n + 1)
        xx = 5.0 ** (jj / n)
        with np.errstate(all="ignore"):
            VV, _, _, zs = vplus_float(xx, k0)
            v = -0.5 * np.ones_like(xx)
            for i in range(1, N + 1):
                w = zs[i]
                v = (w + 4) / ((w + 2) * (w + 5)) + (w + 6) / ((w + 2) * (w + 5)) * v
            f = xx ** (2 - al) * 2 / (6 + 5.0 ** (-k0 - 1) * xx) * v / VV**2
        out["kappa"][N] = float(np.trapezoid(np.nan_to_num(f), np.log(xx)))
    return out


def rk_float(x, k: int):
    """Floats: R_k(x) = (3/5)^k r_k(5^-k x) from u_k = 1/r_k (exact affine recursion along the orbit s_j = p^j(t))."""
    t = np.asarray(x, float) * 5.0 ** (-k)
    ss = [t]
    for _ in range(k):
        ss.append(np.minimum(ss[-1] * (ss[-1] + 5), 1e150))
    v = 3 + ss[k]
    for j in range(k - 1, -1, -1):
        v = _be(ss[j]) * v - _e(ss[j])
    return 0.6**k / v


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    import sympy as sp
    u, v, w, x, k, a = sp.symbols("u v w x k alpha", positive=True)
    lhs = 1 - (1 + u + v) * (1 - u) * (1 - v)
    check(sp.expand(lhs - (u**2 * (1 - v) + v**2 * (1 - u) + u * v)) == 0,
          "1/((1-u)(1-v)) - 1 - u - v = [u^2(1-v) + v^2(1-u) + uv]/((1-u)(1-v)) >= 0 on [0,1)^2, so with omega^0 >= 1: "
          "omega^0/((1-br)(1-br')) - 1 >= b(r + r')  (EE.1's F(r') - F(r) = (r'-r)[...])")
    check(sp.simplify((sp.Rational(9, 5)) ** k * (sp.Rational(5, 3)) ** (2 * k) - 5**k) == 0,
          "scaling: (9/5)^k (r_(k+1) - r_k)(r_k + r_(k+1)) = 5^k (Rh_k - R_k)(R_k + Rh_k)  (r_j = (5/3)^k R at t = 5^-k x)")
    rho = 5**a / 3
    check(sp.simplify(rho ** (-k) * 5 ** (-k * (1 - a)) * (sp.Rational(25, 9)) ** k * 5 ** (-k) - 3 ** (-k)) == 0,
          "scaling for eta: rho^-j 5^(-j(1-a)) (5/3)^(2j) 5^-j = 3^-j  (rho^-j int t^-a b r_j^2 eta dt = 3^-j int x^-a (5^j b) R_j^2 eta dx)")
    s = sp.Symbol("s", positive=True)
    check(sp.expand(s * (s + 5) - 5 * s - s**2) == 0, "p(q) = q^2 + 5q >= 5q, so q(t) <= t/5 and 5^k b = 2x/(q+6) >= 2x/(6 + 5^(-k-1) x) for t = 5^-k x")
    i = sp.Symbol("i", integer=True, nonnegative=True)
    rp = sp.Symbol("rp", positive=True)
    check(sp.simplify(sp.summation(sp.Rational(1, 3) ** i, (i, 0, sp.oo)) - sp.Rational(3, 2)) == 0,
          "(3/rho)^n eps_n = 3^n sum_(k>=n) rho^-k sum_(j>=k) rho'^j c_j >= c_lo 3^n sum_(k>=n) 3^-k/(1-rho') = (3/2) c_lo/(1-rho')")
    # EB.1 at rho^3 = 3: rho^n (Y_n - Y_oo) = (1-d)/(1+nu X)(|nu|^3 - E/Dn) - H/Dn, E = (3/rho)^n eps_n, H = 3^n eta_n
    om, nu, X, E, H, Dn, Cd = sp.symbols("om nu X E H Dn Cd", real=True)
    r = sp.Symbol("r", positive=True)
    n = sp.Symbol("n", positive=True)
    Cdel = -nu * om
    Xn = Cdel / (1 + nu * X) - (om * E * (r / 3) ** n / (1 + nu * X) + r**n * H * 3 ** (-n)) / Dn        # X = rho^-n, eps = E rho'^n, eta = H 3^-n
    Yn = r**n * (Xn - Cdel)
    Yoo = nu**2 * om
    expr = (r**n * (Yn - Yoo)).subs(X, r ** (-n))
    target = om / (1 + nu * r ** (-n)) * (-nu**3 - E / Dn) - H / Dn
    check(sp.simplify((expr - target).subs(r, sp.Integer(3) ** sp.Rational(1, 3))) == 0,
          "EB.1 at rho^3 = 3: rho^n (Y_n - Y_oo) = (1-delta_N)(|nu|^3 - (3/rho)^n eps_n/D'_n)/(1 + nu rho^-n) - 3^n eta_n/D'_n  (nu < 0)")
    gen = (Yn - Yoo).subs(X, r ** (-n))
    tgen = (r**2 / 3) ** n * (om / (1 + nu * r ** (-n)) * (-nu**3 * (3 / r**3) ** n - E / Dn) - H / Dn)
    check(sp.simplify((gen - tgen).subs({r: sp.Rational(3, 2), n: 7})) == 0 and sp.simplify((gen - tgen).subs({r: sp.Rational(6, 5), n: 3})) == 0,
          "EB.1, any rho: Y_n - Y_oo = (rho^2/3)^n [(1-delta_N)(|nu|^3 (3/rho^3)^n - (3/rho)^n eps_n/D'_n)/(1+nu rho^-n) - 3^n eta_n/D'_n]")
    # exact formula for 1/R_k (EE.1 proof) at k = 2 against the recursion
    t = sp.Symbol("t", positive=True)
    be = lambda z: (z + 6) / ((z + 2) * (z + 5))
    ee_ = lambda z: 2 * z / (2 + z)
    p = lambda z: z * (z + 5)
    s0, s1, s2 = t, p(t), p(p(t))
    u2_rec = be(s0) * (be(s1) * (3 + s2) - ee_(s1)) - ee_(s0)
    u2_sum = be(s0) * be(s1) * (3 + s2) - ee_(s0) - ee_(s1) * be(s0)
    check(sp.simplify(u2_rec - u2_sum) == 0,
          "u_k = prod_(j<k) beta(s_j)(3+s_k) - sum_(j<k) e(s_j) prod_(i<j) beta(s_i) (checked at k = 2): 1/R_k = A_k(3+phi_k) - sum_m e(s_(k-m))(5/3)^m P_(k,m)")


# ----------------------------------------------------------------- (B) arb
def _arb():
    from flint import arb, ctx
    return arb, ctx


def _ef():
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import gasket_fractional_corner_critical_sign as ef
    return ef


def node_arb(x, k0s, Nmax: int = 0, extra: int = 40):
    """Rigorous data at a point x (arb): Phi(x) (ball), P_m(x) = prod_(m'>m) (5/3) beta(Phi(5^-m' x)) for m = 0..max(k0s)+1 (P_0 = A),
    eta-orbit z_m, and for each k0: A_k0(x) and e(phi_(k0-m)(5^-m x)), m = 1..k0 (exact finite orbits).  Phi(y) in [y, y e^(y e^0.1/20)]
    for y <= 1.8 (ED.2 (a)), z_(m-1) = p(z_m); P-tail factor in [exp(-0.7 sum_(m>M) z_m), 1] (EE's prof_arb)."""
    arb, _ = _arb()
    xu = x.upper()
    M = int(math.ceil(math.log(max(float(xu), 1e-300)) / math.log(5) + extra))
    five = arb(5)
    c = (arb(1) / 10).exp() / 20
    yv = x / five**M
    z = arb(yv.lower()).union(yv.upper() * (yv.upper() * c).exp())
    zs = [z]
    for _ in range(M):
        z = z * (z + 5)
        zs.append(z)
    zs.reverse()
    f53 = five / 3
    yM = arb(xu) / five ** (M + 1)
    tau = arb("0.7") * (yM * c).exp() * arb(xu) / five**M / 4
    Ptail = arb((-tau).exp().lower()).union(arb(1))
    kmax = max(k0s) + 1
    P = [None] * (M + 1)
    acc = Ptail
    P[M] = acc
    for m in range(M - 1, -1, -1):
        zm = zs[m + 1]
        acc = acc * f53 * (zm + 6) / ((zm + 2) * (zm + 5))
        P[m] = acc
    out = {"Phi": zs[0], "P": P[: kmax + 2], "zs": zs[: max(Nmax, 1) + 1], "Ak": {}, "em": {}}
    for k0 in k0s:
        s = [x / five**k0]
        for _ in range(k0):
            s.append(s[-1] * (s[-1] + 5))
        Ak = arb(1)
        em = []
        for m in range(1, k0 + 1):
            sm = s[k0 - m]
            Ak = Ak * f53 * (sm + 6) / ((sm + 2) * (sm + 5))
            em.append(2 * sm / (2 + sm))
        out["Ak"][k0], out["em"][k0] = Ak, em
    return out


def vcell(na, nb, k0):
    """Upper bound of 1/R_k on the cell [a, b] for every k >= k0: A_k0(a)(3 + Phi(b)) - sum e_m(a) (5/3)^m P_m(b)
    (A_k0, P_m nonincreasing, e_m, Phi nondecreasing in x)."""
    arb, _ = _arb()
    f53 = arb(5) / 3
    V = arb(na["Ak"][k0].upper()) * (3 + arb(nb["Phi"].upper()))
    for m in range(1, k0 + 1):
        V -= arb(na["em"][k0][m - 1].lower()) * f53**m * arb(nb["P"][m].lower())
    return arb(V.upper())


def eta_arb(N, zs):
    arb, _ = _arb()
    v = arb(-1) / 2
    for i in range(1, N + 1):
        w = zs[i]
        v = (w + 4) / ((w + 2) * (w + 5)) + (w + 6) / ((w + 2) * (w + 5)) * v
    return v


def certify_lower(al, k0s, Ns=(), n: int = 200, lo: int = -25, hi: int = 8):
    """c_lo(k0) and kappa_lo(N, k0) as certified LOWER bounds: on each cell [a, b] of a log grid every factor is replaced by its
    value at the end that makes it smaller (monotone, see vcell), times W = int_cell x^(1-a); head (x < 5^lo) and tail (x > 5^hi,
    resp. 5^(N+6)) are dropped (integrands >= 0)."""
    arb, _ = _arb()
    hiK = max([hi] + [N + 6 for N in Ns])
    xs = [arb(float(5.0 ** (lo + i / n))) for i in range((hiK - lo) * n + 1)]
    ks = sorted(set(k0s) | {k + 1 for k in k0s})
    Nmax = max(Ns) if Ns else 0
    nd = [node_arb(xv, ks, Nmax) for xv in xs]
    iD = (hi - lo) * n
    nd5 = [node_arb(5 * xv, ks) for xv in xs[: iD + 1]]
    W = [(xs[i + 1] ** (2 - al) - xs[i] ** (2 - al)) / (2 - al) for i in range(len(xs) - 1)]
    out = {"c": {}, "kappa": {}}
    five = arb(5)
    for k0 in k0s:
        tot = arb(0)
        for i in range(iD):
            na, nb = nd[i], nd[i + 1]
            Rl = 1 / vcell(na, nb, k0)
            Rhl = (five / 3) / vcell(nd5[i], nd5[i + 1], k0 + 1)
            g = 2 / (6 + xs[i + 1] / five ** (k0 + 1)) * arb(nb["P"][0].lower()) * (1 + 1 / (arb(nb["Phi"].upper()) + 5)) * Rl * Rhl * (Rl + Rhl)
            tot += W[i] * arb(g.lower())
        out["c"][k0] = arb(tot.lower())
        for N in Ns:
            iN = (N + 6 - lo) * n
            tk = arb(0)
            for i in range(iN):
                na, nb = nd[i], nd[i + 1]
                Rl = 1 / vcell(na, nb, k0)
                et = eta_arb(N, nb["zs"])
                if not bool(et.lower() > 0):
                    continue                      # eta_N > 0 (BL): dropping the cell keeps a lower bound
                g = 2 / (6 + xs[i + 1] / five ** (k0 + 1)) * Rl * Rl * arb(et.lower())
                tk += W[i] * arb(g.lower())
            out["kappa"][(N, k0)] = arb(tk.lower())
    return out


def rk_arb(x, k):
    arb, _ = _arb()
    t = x / arb(5) ** k
    ss = [t]
    for _ in range(k):
        ss.append(ss[-1] * (ss[-1] + 5))
    v = 3 + ss[k]
    for j in range(k - 1, -1, -1):
        z = ss[j]
        v = (z + 6) / ((z + 2) * (z + 5)) * v - 2 * z / (2 + z)
    return (arb(3) / 5) ** k / v


def fmt(v, k=7):
    return f"[{float(v.lower()):.{k}g}, {float(v.upper()):.{k}g}]"


# ----------------------------------------------------------------- main
def main(full: bool = False) -> int:
    t0 = time.time()
    print("Theorem EG: every-n lower bounds for EE's coefficients and an explicit threshold at alpha = 4gamma/3. Not thrust.")
    symbolic_checks()
    arb, ctx = _arb()
    ctx.prec = 128
    ef = _ef()
    import gasket_fractional_corner_limit_mobius as eb
    n_g, n_p, n_m = (1200, 4000, 20000) if full else (120, 400, 1000)
    # ---- (i) pointwise sanity: 1/V_k0 <= R_k for k >= k0 (exact arb r_k)
    okp, worst = True, 0.0
    for k0 in (1, 2, 4):
        for xf in (1e-3, 0.3, 1.0, 4.0, 20.0, 150.0):
            xv = arb(xf)
            ndd = node_arb(xv, [k0])
            Vb = vcell(ndd, ndd, k0)
            for k in (k0, k0 + 1, k0 + 3, k0 + 8, k0 + 20):
                Rk = rk_arb(xv, k)
                okp = okp and bool(1 / Vb <= Rk)
                worst = max(worst, float((1 / Vb / Rk).upper()))
    check(okp, f"EG.1 pointwise: 1/V_k0(x) <= R_k(x) for k = k0, k0+1, k0+3, k0+8, k0+20 (k0 = 1, 2, 4; six points; exact arb r_k); max ratio {worst:.4f}")
    # ---- EG.3 at alpha = 4gamma/3
    al = ef.alpha_ball()
    rho = arb(5) ** al / 3
    check(bool(abs(rho**3 - 3) < arb(10) ** -30), "alpha = 4gamma/3 as the ball 4 log3/(3 log5); rho^3 = 3 to 1e-30")
    t1 = time.time()
    pi = ef.certify_profile_integrals(al, NS, n=n_p)
    Dp = pi["Dp"]
    Mt, mi = ef.certify_mt(al, n=n_m)
    nu = -rho * Mt / ((rho - 1) * Dp)
    om = {N: pi["h"][N] / Dp for N in NS}
    print(f"   alpha = 4gamma/3 (EF.1 routines re-run): D'_oo in {fmt(Dp)}, nu in {fmt(nu)}, |nu|^3 D'_oo in {fmt(-nu**3 * Dp, 6)}  ({time.time() - t1:.1f}s)")
    check(mi["ok"] and bool(nu < 0) and all(pi["ok_tail"][N] and bool(om[N] > 0) for N in NS), "EF.1 data at 4gamma/3 finite; nu < 0; 1 - delta_N > 0")
    k0s = (1, 2, 3, 4, 6, 8)
    t2 = time.time()
    lc = certify_lower(al, k0s, Ns=NS, n=n_g)
    rp = rho / 3
    Ls = {k0: 3 * lc["c"][k0] / (2 * (1 - rp)) for k0 in k0s}
    print(f"   c_lo(k0), L(k0) at 4gamma/3 ({time.time() - t2:.1f}s): " + "; ".join(f"k0={k0}: c >= {float(lc['c'][k0].lower()):.6f}, L >= {float(Ls[k0].lower()):.6f}" for k0 in k0s))
    thr = -nu**3 * Dp
    check(bool(Ls[2] > thr), f"EG.3: L(2) > |nu|^3 D'_oo ({float(Ls[2].lower()):.5f} > {float(thr.upper()):.5f}): Y_n < Y_oo for every N >= 1 and every n >= max(N, 2)")
    check(not bool(Ls[1] > thr), "EG.3: L(1) alone does not exceed |nu|^3 D'_oo (so n = 1 needs the eta term)")
    x1 = 1 / rho
    check(bool(1 + nu * x1 > 0), f"EB.1 applies at every n >= 1: 1 + nu rho^-n >= 1 + nu/rho in {fmt(1 + nu * x1, 6)} > 0 (nu < 0, so also <= 1)")
    S11 = om[1] * (-nu**3 - Ls[1] / Dp) / (1 + nu * x1) - 3 * lc["kappa"][(1, 1)] / (2 * Dp)
    print(f"   N = n = 1, k0 = 1: kappa_lo(1,1) >= {float(lc['kappa'][(1, 1)].lower()):.6f}; rho (Y_1 - Y_oo) <= {float(S11.upper()):.6g}")
    check(bool(S11 < 0), "EG.3: at N = 1, n = 1: (1-delta_1)(|nu|^3 - L(1)/D'_oo)/(1+nu/rho) - 1.5 kappa_lo(1,1)/D'_oo < 0, so Y_1 < Y_oo")
    check(all(bool(Ls[k0] > thr) for k0 in k0s if k0 >= 2), "EG.3: L(k0) > |nu|^3 D'_oo for every certified k0 >= 2 (monotone check)")
    Splus, rlo = {}, {}
    for N in NS:
        Splus[N] = om[N] * (-nu**3 - Ls[6] / Dp) - 3 * lc["kappa"][(N, 6)] / (2 * Dp)
        rlo[N] = (om[N] * Ls[6] + 3 * lc["kappa"][(N, 6)] / 2) / Dp
    print("   every n >= max(N, 6): rho^n (Y_n - Y_oo) <= S_N^+, 3^n R_n >= varrho_lo: " + "; ".join(
        f"N={N}: S^+ <= {float(Splus[N].upper()):.6g}, varrho_lo >= {float(rlo[N].lower()):.6g}" for N in NS))
    check(all(bool(Splus[N] < 0) for N in NS), "EG.3: S_N^+ < 0 for N = 1, 2, 3, 6, 10 (every-n quantitative bound, n >= max(N, 6))")
    # consistency with EF.3's limits
    efS = {1: "-0.0232863", 2: "-0.174497", 3: "-0.454077", 6: "-4.68928", 10: "-88.9943"}
    efE_hi, efR_hi = arb("0.7399065"), {1: "0.0280926", 2: "0.197604", 3: "0.489632", 6: "4.7469"}
    check(all(bool(Ls[k0] < efE_hi) for k0 in k0s), "consistency: every L(k0) < EF.3's upper end of e*(4gamma/3) = 0.7399065")
    check(all(bool(Splus[N].upper() >= arb(efS[N])) for N in NS) and all(bool(rlo[N] < arb(efR_hi[N])) for N in efR_hi),
          "consistency: S_N^+ >= EF.3's S_N and varrho_lo <= EF.3's varrho_N (limits obey the every-n bounds)")
    RESULT["crit"] = {"Dp": Dp, "nu": nu, "om": om, "L": Ls, "c": lc["c"], "kappa": lc["kappa"], "S11": S11, "Splus": Splus, "rlo": rlo, "thr": thr}
    # ---- EG.4 other alphas
    RESULT["other"] = {}
    for a_s in ("0.7", "0.8", "0.9", "0.95"):
        t3 = time.time()
        a = arb(a_s)
        ra = arb(5) ** a / 3
        ks = (1, 2, 4, 6) if a_s != "0.8" else (1, 2, 4, 6, 8)
        Ns_a = NS if a_s == "0.8" else ()
        la = certify_lower(a, ks, Ns=Ns_a, n=n_g)
        La = {k0: 3 * la["c"][k0] / (2 * (1 - ra / 3)) for k0 in ks}
        RESULT["other"][a_s] = {"L": La, "c": la["c"], "kappa": la["kappa"]}
        print(f"   alpha = {a_s} ({time.time() - t3:.1f}s): " + "; ".join(f"k0={k0}: L >= {float(La[k0].lower()):.6f}" for k0 in ks))
        check(bool(La[1] > 0) and all(bool(La[ks[i + 1]] > La[ks[i]]) for i in range(len(ks) - 1)), f"EG.4 alpha = {a_s}: L(k0) > 0 and increasing in k0 (certified)")
    # 0.8: every-n lower bounds for 3^n R_n vs ED.4 and EE.4
    nu8, A8 = ef.hull(arb(eb.NU08[0]), arb(eb.NU08[1])), ef.hull(arb(eb.A08[0]), arb(eb.A08[1]))
    a8 = arb("0.8")
    D8 = A8 / (2 * (arb.pi() * a8).sin() / arb.pi())
    o8 = RESULT["other"]["0.8"]
    ed4lo = {1: 3.716e-5, 2: 1.887e-4, 3: 3.066e-4, 6: 5.480e-4, 10: 7.155e-4}
    ee4hi = {1: 0.0288038, 2: 0.239288, 3: 0.706085, 6: 11.6907, 10: 446.198}
    r8 = {}
    for N in NS:
        omN = ef.hull(arb(eb.CD_EA[N][0]), arb(eb.CD_EA[N][1])) / (-nu8)
        r8[N] = (omN * o8["L"][8] + 3 * o8["kappa"][(N, 8)] / 2) / D8
    print("   alpha = 0.8, every n >= max(N, 8): 3^n R_n >= " + "; ".join(f"N={N}: {float(r8[N].lower()):.6g} (ED.4: {ed4lo[N]:.4g})" for N in NS))
    check(all(float(r8[N].lower()) > ed4lo[N] and float(r8[N].upper()) < ee4hi[N] for N in NS),
          "EG.4 at 0.8: every-n lower bound of 3^n R_n beats ED.4's lower column and stays below EE.4's limit (N = 1, 2, 3, 6, 10)")
    check(bool(o8["L"][8] < arb("0.7410025")) and bool(o8["L"][8] > arb("0.0055320") * 100),
          "EG.4 at 0.8: L(8) below EE.4's e* (0.7410025) and > 100 x ED.2's L_0.8 = 0.005532")
    RESULT["r08"] = r8
    # 0.95: rho^3 > 3, |nu|^3 < L(1)/D'_oo => Y_n < Y_oo for every n >= N
    a95 = arb("0.95")
    r95 = arb(5) ** a95 / 3
    pi95 = ef.certify_profile_integrals(a95, Ns=(1,), n=n_p)
    Mt95, _ = ef.certify_mt(a95, n=n_m)
    nu95 = -r95 * Mt95 / ((r95 - 1) * pi95["Dp"])
    L95 = RESULT["other"]["0.95"]["L"][1]
    print(f"   alpha = 0.95: rho^3 in {fmt(r95 ** 3, 6)}, nu in {fmt(nu95, 6)}, |nu|^3 D'_oo in {fmt(-nu95 ** 3 * pi95['Dp'], 6)} vs L(1) >= {float(L95.lower()):.6f}")
    check(bool(r95**3 > 3) and bool(1 + nu95 / r95 > 0) and bool(L95 > -nu95**3 * pi95["Dp"]), "EG.4 at 0.95: rho^3 > 3, 1 + nu/rho > 0 and L(1) > |nu|^3 D'_oo, so Y_n < Y_oo for every N >= 1 and every n >= N")
    RESULT["a095"] = {"nu": nu95, "Dp": pi95["Dp"], "L1": L95}
    # ---- EG.5 floats
    fl = {k0: lower_float(4 * GAMMA / 3, k0, Ns=(1,) if k0 == 1 else ())["L"] for k0 in (1, 2, 3, 4, 6, 8, 10)}
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import gasket_fractional_corner_remainder_limit as ee
    al8 = 0.8
    d8 = {k: ee.direct_increment_float(al8, k) for k in (2, 4, 6)}
    c8f = {k: lower_float(al8, k)["c"] for k in (2, 4, 6)}
    print("   floats: L(k0) at 4gamma/3: " + ", ".join(f"{k}: {v:.6f}" for k, v in fl.items())
          + "; at 0.8 c_lo(k) vs direct c_k: " + ", ".join(f"k={k}: {c8f[k]:.6f} <= {d8[k]:.6f}" for k in d8))
    check(all(c8f[k] < d8[k] for k in d8), "EG.5 (OBSERVED): c_lo(k) below the direct quadrature of c_k at 0.8 (k = 2, 4, 6)")
    check(all(bool(Ls[k0].lower() <= fl[k0] * (1 + 1e-9)) for k0 in k0s if k0 in fl), "EG.5: float L(k0) at or above the certified lower bound (grid loss only)")
    RESULT["floats"] = {"L": fl, "c8": c8f, "d8": d8}
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EG checks passed")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    sys.exit(main(full="--full" in sys.argv))
