#!/usr/bin/env python3
"""Theorem EJ: a two-sided every-level rate, and the exact sign-change level at alpha = 0.9.  EI kept "upper side only through EC.2's
crude B_alpha (~1.4e3 e*) and DZ.1's crude etabar, so n1 is an upper bound on the onset, not the onset; exact sign-change position
and number of changes undecided"; EG kept the open item "an every-n upper bound for c_k (a two-sided rate)".  EG.1's monotone
route runs backwards: every factor of EE.1's scaled integrand has an explicit every-level UPPER bound as well.

  EJ.1 (THEOREM; every alpha in (gamma, 1), k0 >= 1, 1 <= N <= k0; AQ (1)-(3), BI.4, BK.1-2, BL, DZ.1-2, ED.1, EE.1's exact formulas
       as stated; spectral r_k(s) <= 1/s).  For every k >= k0, x > 0:  R_k(x) <= R^hi(x) = min(1/3, G x^(gamma-1), 1/V^lo(x)),
         V^lo = A(x)(3 + phi_k0(x)) - sum_(m>=1) e(Phi(5^-m x)) (5/3)^m Q_m(x),  Q_m = prod_(m<m'<=k0) (5/3) beta(phi_(k0-m')(5^-m' x)),
       Rh_k - R_k <= R^hi Rh^hi A_k0 (phi_k0+6)/(phi_k0+5); with omega^0 <= e^(0.9q), b r_k <= 2/(q+6) <= 1/3 and, for t >= 1,
       u_k - u_(k+1) = prod_(j<k) beta(p^j t) (p^k t + 6)/(p^k t + 5) (sympy), explicit c_hi(k0), kappa_hi(N, k0) with
         c_k = (3/rho)^k (Mt_(k+1) - Mt_k) <= c_hi(k0),   (3/rho)^j int t^-a F(r_j) eta_N(p^(j-N) t) dt <= kappa_hi(N, k0)   (k, j >= k0),
       hence for every n >= k0:  L(k0) <= (3/rho)^n eps_n <= U(k0) = 3 c_hi/(2(1-rho/3)),  1.5 kappa_lo <= 3^n eta_n <= 1.5 kappa_hi.
  EJ.2 (THEOREM, exact on EB.1): EI.1's criterion with U, 1.5 kappa_hi in place of B_alpha, H_k (monotone in k, so one level n1
       suffices), and the mirror criterion with L, kappa_lo for Y_n < Y_oo (sign-split on |nu|^3 g^n - L/D', no monotone-in-n claim
       needed: each n is checked).
  EJ.3 (arb, 128 bits; --full): at alpha = 0.9: Y_n < Y_oo for n = 1..32 and Y_n > Y_oo for every n >= 33 (N = 1): exactly one
       sign change, the last level with Y_n < Y_oo is n = 32 (EI: [8, 158]); N = 2: 2..40 / >= 41; N = 3: 3..50 / >= 51; at 0.905
       (N = 1): 1..67 / >= 68.
  EJ.4 (arb): onset levels n1 and negative windows at EI's eight exponents (N = 1), e.g. 0.88: [2, 8] and >= 10 (EH.4: no window).
  EJ.5 (OBSERVED, floats): the float bounds and EI.4's limit-substituted crossing.

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
G_BI4 = 0.6031
ALPHAS = ("0.7", "0.8", "0.85", "0.88", "0.89", "0.9", "0.905", "0.91")
EI_N1 = {"0.7": 7, "0.8": 13, "0.85": 25, "0.88": 52, "0.89": 79, "0.9": 159, "0.905": 316, "0.91": 11584}   # EI.2, N = 1 (--full)
EXACT = ("0.9", "0.905")
NS = (1, 2, 3)
K0_LO = (1, 2, 3, 4, 6, 8, 10, 12)


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def _mods():
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import gasket_fractional_corner_critical_sign as ef
    import gasket_fractional_corner_explicit_threshold as eg
    return ef, eg


def _be(s):
    return (s + 6) / ((s + 2) * (s + 5))


def _e(s):
    return 2 * s / (2 + s)


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
def _orbit(x, k0):
    s = [np.asarray(x, float) * 5.0 ** (-k0)]
    for _ in range(k0):
        s.append(np.minimum(s[-1] * (s[-1] + 5), 1e150))
    return s


def rhi_float(x, k0: int, M: int = 60):
    """Floats: R^hi(x) = min(1/3, G x^(gamma-1), 1/V^lo(x)) (EJ.1), V^lo = A (3 + phi_k0) - sum_m e(Phi(5^-m x)) (5/3)^m Q_m."""
    x = np.asarray(x, float)
    zs = [x * 5.0 ** (-M)]
    for _ in range(M):
        zs.append(np.minimum(zs[-1] * (zs[-1] + 5), 1e150))
    zs = zs[::-1]
    A = np.ones_like(x)
    for m in range(1, M + 1):
        A = A * (5 / 3) * _be(zs[m])
    s = _orbit(x, k0)
    Q = [np.ones_like(x) for _ in range(M + 1)]
    acc = np.ones_like(x)
    for m in range(k0 - 1, -1, -1):
        acc = acc * (5 / 3) * _be(s[k0 - m - 1])
        Q[m] = acc
    V = A * (3 + s[k0])
    for m in range(1, M + 1):
        V = V - _e(zs[m]) * (5 / 3) ** m * Q[m]
    with np.errstate(divide="ignore"):
        inv = np.where(V > 0, 1 / np.where(V > 0, V, 1), np.inf)
    return np.minimum(np.minimum(1 / 3, G_BI4 * x ** (GAMMA - 1)), inv), Q[0], s


def upper_float(al: float, k0: int, Ns=NS, n: int = 100, lo: int = -25, hi: int = 32):
    """Floats (trapezoid in log x) of EJ.1's c_hi(k0), U(k0) and kappa_hi(N, k0) (main parts only, T = 1)."""
    j = np.arange(lo * n, hi * n + 1)
    x = 5.0 ** (j / n)
    R, Q0, s = rhi_float(x, k0)
    Rh = (5 / 3) * rhi_float(5 * x, k0 + 1)[0]
    E = np.exp(0.18 * np.minimum(1.0, x * 5.0 ** (-k0)))
    u = np.minimum(1 / 3, 3.0 ** (-k0) * x / 3 * R)
    v = np.minimum(1 / 3, 3.0 ** (-k0) * x / 3 * Rh)
    D = R * Rh * Q0 * (1 + 1 / (s[k0] + 5))
    brk = 0.18 * x * 0.6**k0 * E + E * x / 3 * (R + Rh) / ((1 - u) * (1 - v))
    c = float(np.trapezoid(x ** (1 - al) * D * brk, np.log(x)))
    rho = 5**al / 3
    out = {"c": c, "U": 1.5 * c / (1 - rho / 3), "kappa": {}}
    for N in Ns:
        et = -0.5 * np.ones_like(x)
        for i in range(1, N + 1):
            w = s[k0 - i]
            et = (w + 4) / ((w + 2) * (w + 5)) + (w + 6) / ((w + 2) * (w + 5)) * et
        f = x ** (1 - al) * R * (0.18 * x * 0.6**k0 * E + x / 3 * R) / (1 - u) * np.minimum(1.0, et)
        out["kappa"][N] = float(np.trapezoid(f, np.log(x)))
    return out


def rk_float(x, k: int):
    """Floats: R_k(x) = (3/5)^k r_k(5^-k x) from the exact affine recursion u_(k+1)(s) = beta(s) u_k(ps) - e(s), u_0 = 3 + s."""
    ss = _orbit(x, k)
    v = 3 + ss[k]
    for j in range(k - 1, -1, -1):
        v = _be(ss[j]) * v - _e(ss[j])
    return 0.6**k / v


def crossing_float(al: float, U: float, kap: float, N: int = 1):
    """Floats: EJ.2's two criteria with EF's float D'_oo, nu, 1-delta_N and the given U, kappa (upper) -- first n with
    om (|nu|^3 g^n - U/D') - 1.5 kappa/D' > 0."""
    ef, _ = _mods()
    c = ef.critical_floats(al, Ns=(N,))
    g = 3 / (5**al / 3) ** 3
    nu, D, om = -c["nu"], c["Dp"], c["om"][N]
    n = 1
    while om * (nu**3 * g**n - U / D) - 1.5 * kap / D <= 0 and n < 10**6:
        n += 1
    return n, c


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    import sympy as sp
    s, q, t, u, v, a, x, k = sp.symbols("s q t u v alpha x k", positive=True)
    be = (s + 6) / ((s + 2) * (s + 5))
    e = 2 * s / (2 + s)
    d01 = sp.simplify((3 + s) - (be * (3 + s * (s + 5)) - e))
    check(sp.simplify(d01 - (s + 6) / (s + 5)) == 0, "sympy: u_0 - u_1 = (3+s) - (beta(s)(3+p(s)) - e(s)) = (s+6)/(s+5) in [1, 6/5] for s >= 0")
    check(sp.simplify(3 * (s + 2) * (s + 5) - s * (s + 6) - (2 * s**2 + 15 * s + 30)) == 0,
          "sympy: 3(s+2)(s+5) - s(s+6) = 2s^2 + 15s + 30 > 0, so beta(s) <= 3/s; beta decreasing, beta(0) = 3/5")
    expr = 1 / ((1 - u) * (1 - v)) - 1 - (u + v) / ((1 - u) * (1 - v))
    check(sp.simplify(expr + u * v / ((1 - u) * (1 - v))) == 0, "sympy: 1/((1-u)(1-v)) - 1 = (u+v)/((1-u)(1-v)) - uv/((1-u)(1-v)) <= (u+v)/((1-u)(1-v))")
    lw = (1 + a) * sp.log(1 + q / 5) + sp.log(1 + q / 2) - sp.log(1 + 2 * q / 5) - sp.log(1 + q / 6)
    check(sp.simplify(sp.diff((1 + a) * q / 5 + q / 2 - (1 + a) * sp.log(1 + q / 5) - sp.log(1 + q / 2), q)
                      - ((1 + a) * q / (5 * (q + 5)) + q / (2 * (q + 2)))) == 0,
          "sympy: d/dq[(1+a)q/5 + q/2 - (1+a)log(1+q/5) - log(1+q/2)] >= 0, so log omega^0 <= (1+a)q/5 + q/2 <= 0.9 q (a < 1)")
    check(sp.simplify(5 * (1 + 2 * q / 5) - 4 * (1 + q / 2)) == 1 and sp.simplify(6 * (1 + q / 6) - 5 * (1 + q / 5)) == 1,
          "sympy: (1+q/2)/(1+2q/5) <= 5/4 and (1+q/5)/(1+q/6) <= 6/5, so omega^0 <= (1+q/5)^2 (5/4)/(1+q/6) <= 1.5 (1+q/5)")
    qt = 2 * t / (5 + sp.sqrt(25 + 4 * t))
    check(sp.simplify(qt * (qt + 5) - t) == 0, "sympy: q(t) = 2t/(5+sqrt(25+4t)) inverts p; q <= t/5 (p(q) >= 5q) and q <= sqrt(t) (p(sqrt t) >= t)")
    check(sp.simplify(sp.Rational(9, 5) ** k * sp.Rational(5, 3) ** k - 3**k) == 0 and sp.simplify(3**k * sp.Rational(5, 3) ** k - 5**k) == 0,
          "sympy: (9/5)^k (5/3)^k = 3^k (Psi_k = 3^k x^-a (Rh-R)[..]) and 3^k (5/3)^k = 5^k (3^k u = 5^k b R_k)")
    g = sp.log(3) / sp.log(5)
    check(sp.simplify(sp.Rational(5, 3) * 5 ** (g - 1) - 1) == 0, "sympy: (5/3) 5^(gamma-1) = 1, so Rh_k = (5/3) R_(k+1)(5x) <= G x^(gamma-1) too")
    rho = sp.symbols("rho", positive=True)
    w = 3**k * sum(rho ** (-j) for j in range(0, 4))
    check(sp.simplify(sp.Sum(3 ** (-sp.Symbol("i")), (sp.Symbol("i"), 0, sp.oo)).doit() - sp.Rational(3, 2)) == 0,
          "sympy: sum_(i>=0) 3^-i = 3/2 (eps_n and eta_n series as in EG.1 (e))")
    del w, x


# ----------------------------------------------------------------- (B) arb
def _arb():
    from flint import arb, ctx
    return arb, ctx


def _min(vals):
    arb, _ = _arb()
    best = None
    for v in vals:
        w = arb(v.upper())
        if best is None or bool(w < best):
            best = w
    return best


def node_phi(x, M):
    """Phi-orbit data at x (arb): z_m = Phi(5^-m x) balls (ED.2 (a) enclosure at 5^-M x, z_(m-1) = p(z_m)), A(x) lower bound
    (tail factor e^-tau, EE's prof_arb), e(z_m) and the S-tail sum_(m>M) e(z_m)(5/3)^m <= sfac x 3^-M / 2."""
    arb, _ = _arb()
    five = arb(5)
    c = (arb(1) / 10).exp() / 20
    xu = x.upper()
    yv = x / five**M
    z = arb(yv.lower()).union(yv.upper() * (yv.upper() * c).exp())
    zs = [z]
    for _ in range(M):
        z = z * (z + 5)
        zs.append(z)
    zs.reverse()
    f53 = five / 3
    A = arb(1)
    for m in range(1, M + 1):
        A = A * f53 * _be(zs[m])
    sfac = (arb(xu) / five ** (M + 1) * c).exp()
    tau = arb("0.7") * sfac * arb(xu) / five**M / 4
    return {"Alo": arb((A * (-tau).exp()).lower()), "eu": [None] + [arb(_e(zs[m]).upper()) for m in range(1, M + 1)],
            "Stail": arb((sfac * arb(xu) / arb(3) ** M / 2).upper()), "M": M}


def node_k0(x, k0):
    """Level-k0 orbit s_j = p^j(5^-k0 x), Q_m (m < k0; Q_0 = A_k0) at x."""
    arb, _ = _arb()
    five = arb(5)
    s = [x / five**k0]
    for _ in range(k0):
        s.append(s[-1] * (s[-1] + 5))
    Q = {}
    acc = arb(1)
    for m in range(k0 - 1, -1, -1):
        acc = acc * five / 3 * _be(s[k0 - m - 1])
        Q[m] = arb(acc.upper())
    return {"s": s, "Q": Q}


def vlo_cell(pa, pb, ka, k0):
    """Lower bound of 1/R_k on the cell [a, b] for every k >= k0: A(b)(3 + phi_k0(a)) - sum_m e(z_m(b)) (5/3)^m Q_m(a) - S-tail(b)
    (A, Q_m nonincreasing; phi_k0, e(Phi) nondecreasing in x)."""
    arb, _ = _arb()
    f53 = arb(5) / 3
    V = pb["Alo"] * (3 + arb(ka["s"][k0].lower()))
    pw = arb(1)
    for m in range(1, pb["M"] + 1):
        pw = pw * f53
        V -= pb["eu"][m] * pw * (ka["Q"][m] if m < k0 else 1)
    V -= pb["Stail"]
    return arb(V.lower())


def eta_fin(N, s, k0):
    """eta_N at z = phi_(k0-N)(5^-N x) = s_(k0-N), by BL's recursion along the finite orbit (z_i = s_(k0-i))."""
    arb, _ = _arb()
    v = arb(-1) / 2
    for i in range(1, N + 1):
        w = s[k0 - i]
        v = (w + 4) / ((w + 2) * (w + 5)) + (w + 6) / ((w + 2) * (w + 5)) * v
    return v


def certify_upper(al, k0s, Ns=NS, n: int = 400, lo: int = -25, hi: int = 32):
    """EJ.1's c_hi(k0), U(k0), kappa_hi(N, k0) as certified UPPER bounds (T = 1): cells of a log grid, every factor at the end that
    makes it larger; explicit head (x < 5^lo), tail (x > 5^hi) and t >= 1 pieces."""
    arb, _ = _arb()
    g = arb(3).log() / arb(5).log()
    five, third, G = arb(5), arb(1) / 3, arb("0.6031")
    rho = five**al / 3
    T = arb(1)
    xs = [arb(float(5.0 ** (lo + i / n))) for i in range((hi - lo) * n + 1)]
    M = hi + 42
    ph = [node_phi(xv, M) for xv in xs]
    ph5 = [node_phi(5 * xv, M + 1) for xv in xs]
    W0 = [(xs[i + 1] ** (1 - al) - xs[i] ** (1 - al)) / (1 - al) for i in range(len(xs) - 1)]
    out = {}
    for k0 in k0s:
        kn = [node_k0(xv, k0) for xv in xs]
        kn5 = [node_k0(5 * xv, k0 + 1) for xv in xs]
        f35, i3 = (arb(3) / 5) ** k0, arb(3) ** (-k0)
        c = arb(0)
        kap = {N: arb(0) for N in Ns if N <= k0}
        nneg = 0
        for i in range(len(xs) - 1):
            a, b = xs[i], xs[i + 1]
            V = vlo_cell(ph[i], ph[i + 1], kn[i], k0)
            R = _min([third, G * a ** (g - 1)] + ([1 / V] if bool(V > 0) else []))
            V5 = vlo_cell(ph5[i], ph5[i + 1], kn5[i], k0 + 1)
            R5 = _min([third, G * (5 * a) ** (g - 1)] + ([1 / V5] if bool(V5 > 0) else []))
            nneg += 0 if bool(V > 0) else 1
            Rh = arb((five / 3 * R5).upper())
            E = (arb("0.18") * _min([T, b / five**k0])).exp()
            u = _min([third, i3 * b / 3 * R])
            v = _min([third, i3 * b / 3 * Rh])
            sk = arb(kn[i]["s"][k0].lower())
            Dl = R * Rh * kn[i]["Q"][0] * (1 + 1 / (sk + 5))
            brk = arb("0.18") * b * f35 * E + E * b / 3 * (R + Rh) / ((1 - u) * (1 - v))
            c += W0[i] * Dl * brk
            for N in kap:
                eh = _min([arb(1), eta_fin(N, kn[i]["s"], k0)])
                kap[N] += W0[i] * R * (arb("0.18") * b * f35 * E + b / 3 * R) / (1 - u) * eh
        X0, X1 = xs[0], xs[-1]
        E1, ET = (arb("0.18") * X0).exp(), (arb("0.18") * T).exp()
        head = arb(2) / 9 * E1 * (arb("0.18") + third * (third + arb(5) / 9) * arb(9) / 4) * X0 ** (2 - al) / (2 - al)
        tail = G**2 * 30 * ET * (arb("0.18") * f35 * X1 ** (2 * g - 1 - al) / (1 + al - 2 * g)
                                 + G * arb(3) / 2 * X1 ** (3 * g - 2 - al) / (2 + al - 3 * g))
        pT = [T]
        for _ in range(k0 + 1):
            pT.append(pT[-1] * (pT[-1] + 5))
        Pb = arb(1)
        for jj in range(k0):
            Pb = Pb * _be(pT[jj])
        cond = bool((3 / rho) * _be(pT[k0]) <= 1)
        Iw = arb("1.5") * (T ** (-al - 1) / (al + 1) + T ** (-al - arb(1) / 2) / (5 * (al + arb(1) / 2)))
        ttail = (3 / rho) ** k0 * arb(27) / 10 * Pb * Iw
        chi = arb((c + head + tail + ttail).upper())
        row = {"c_main": c, "head": head, "tail": tail, "ttail": ttail, "c_hi": chi, "U": arb((3 * chi / (2 * (1 - rho / 3))).upper()),
               "cond": cond, "nneg": nneg, "kappa": {}}
        for N in kap:
            hk = third * (arb("0.18") * E1 + arb(1) / 9) * arb(3) / 2 * X0 ** (2 - al) / (2 - al)
            tk = 3 * five**N * G * (arb("0.18") * f35 * ET * X1 ** (g - al) / (al - g) + G / 3 * X1 ** (2 * g - 1 - al) / (1 + al - 2 * g))
            Ft = arb("2.25") * (T ** (-al) / al + T ** (arb(1) / 2 - al) / (5 * (al - arb(1) / 2)))
            ttk = (3 / rho) ** k0 * 2 / pT[k0 - N] * Ft
            row["kappa"][N] = {"main": kap[N], "rest": hk + tk + ttk, "hi": arb((kap[N] + hk + tk + ttk).upper())}
        out[k0] = row
    return out


def certify_point(a_s, n_p, n_m):
    """EF.1's cell routines at alpha: D'_oo, Mt_oo, h_N (balls) -> nu, 1 - delta_N."""
    arb, _ = _arb()
    ef, _ = _mods()
    al = arb(a_s)
    rho = arb(5) ** al / 3
    pi = ef.certify_profile_integrals(al, NS, n=n_p)
    Mt, mi = ef.certify_mt(al, n=n_m)
    Dp = pi["Dp"]
    nu = -rho * Mt / ((rho - 1) * Dp)
    om = {N: pi["h"][N] / Dp for N in NS}
    ok = mi["ok"] and all(pi["ok_tail"][N] for N in NS)
    return {"al": al, "rho": rho, "g": 3 / rho**3, "Dp": Dp, "Mt": Mt, "nu": nu, "om": om, "ok": ok}


def onset_upper(P, up, N):
    """EJ.2 (upper side): the smallest n1 >= max(N, k0) (over the given k0) with
    Phi(n) = om_lo (nu_lo^3 g^n - U/D_n) - 1.5 kappa_hi/D_n > 0, D_n = D_lo - rho M_hi rho^-n/(rho-1) > 0 (Phi increasing in n)."""
    arb, _ = _arb()
    rho, g = P["rho"], P["g"]
    Dlo, Mhi = arb(P["Dp"].lower()), arb(P["Mt"].upper())
    nlo = arb((-P["nu"]).lower())
    om = arb(P["om"][N].lower())
    best = None
    for k0, row in up.items():
        if N not in row["kappa"] or not row["cond"]:
            continue
        U, H = row["U"], 3 * row["kappa"][N]["hi"] / 2

        nhi = arb((-P["nu"]).upper())

        def phi(k):
            Dk = Dlo - rho * Mhi / (rho - 1) * rho ** (-k)
            ok = bool(Dk > 0) and bool(nhi < rho**k)      # 0 < 1 + nu rho^-k <= 1
            return om * (nlo**3 * g**k - U / Dk) - H / Dk if ok else arb(-1)

        k1 = max(N, k0)
        hi_ = k1
        while not bool(phi(hi_) > 0):
            hi_ *= 2
            if hi_ > 10**7:
                break
        lo_ = k1
        while lo_ < hi_:
            m = (lo_ + hi_) // 2
            lo_, hi_ = (lo_, m) if bool(phi(m) > 0) else (m + 1, hi_)
        if bool(phi(hi_) > 0) and (best is None or hi_ < best[0]):
            best = (hi_, k0, phi(hi_), phi(hi_ - 1) if hi_ > k1 else None)
    return best


def below_upper(P, low, N, n):
    """EJ.2 (lower side) at one level n >= N: some k0 <= n with braces <= om_X X - 1.5 kappa_lo/D_hi < 0, X = nu_hi^3 g^n - L/D_hi,
    om_X = om_lo if X < 0 else om_hi/(1 - nu_hi rho^-n)."""
    arb, _ = _arb()
    rho, g = P["rho"], P["g"]
    Dhi = arb(P["Dp"].upper())
    nhi = arb((-P["nu"]).upper())
    oml, omh = arb(P["om"][N].lower()), arb(P["om"][N].upper())
    for k0 in sorted(low["L"]):
        if k0 > n or (N, k0) not in low["kappa"]:
            continue
        X = nhi**3 * g**n - low["L"][k0] / Dhi
        den = 1 - nhi * rho ** (-n)                    # 1 + nu rho^-n >= den > 0 needed in both branches
        if not bool(den > 0):
            continue
        if bool(X < 0):
            val = oml * X - 3 * low["kappa"][(N, k0)] / (2 * Dhi)
        else:
            val = omh / den * X - 3 * low["kappa"][(N, k0)] / (2 * Dhi)
        if bool(val < 0):
            return k0, val
    return None


def certify_lower(a_s, n, hi):
    arb, _ = _arb()
    _, eg = _mods()
    al = arb(a_s)
    rho = arb(5) ** al / 3
    o = eg.certify_lower(al, K0_LO, Ns=NS, n=n, hi=hi)
    return {"L": {k: arb((3 * v / (2 * (1 - rho / 3))).lower()) for k, v in o["c"].items()}, "c": o["c"], "kappa": o["kappa"]}


def fmt(v, k=6):
    return f"[{float(v.lower()):.{k}g}, {float(v.upper()):.{k}g}]"


def rk_arb(x, k):
    arb, _ = _arb()
    t = x / arb(5) ** k
    ss = [t]
    for _ in range(k):
        ss.append(ss[-1] * (ss[-1] + 5))
    v = 3 + ss[k]
    for j in range(k - 1, -1, -1):
        v = _be(ss[j]) * v - _e(ss[j])
    return (arb(3) / 5) ** k / v


# ----------------------------------------------------------------- main
def main(full: bool = False) -> int:
    t0 = time.time()
    print("Theorem EJ: a two-sided every-level rate and the exact sign-change level (EG.1 reversed + EB.1). Not thrust.")
    symbolic_checks()
    arb, ctx = _arb()
    ctx.prec = 128
    n_p, n_m = (1000, 4000) if full else (200, 500)
    n_ux, n_u, n_l = (4000, 600, 1000) if full else (1000, 150, 300)
    # ---- attempts to break EJ.1's profile bound: R^hi >= exact R_k (affine recursion) for k = k0..k0+20
    al9 = arb("0.9")
    ok_sp = True
    for k0 in (2, 6, 14):
        for xv in ("0.001", "0.3", "1", "4", "20", "150", "3000", "1e6"):
            x = arb(xv)
            M = 50
            pa, ka = node_phi(x, M), node_k0(x, k0)
            V = vlo_cell(pa, pa, ka, k0)
            R = _min([arb(1) / 3, arb("0.6031") * x ** (arb(3).log() / arb(5).log() - 1)] + ([1 / V] if bool(V > 0) else []))
            for k in range(k0, k0 + 21):
                ok_sp = ok_sp and bool(rk_arb(x, k) <= R)
    check(ok_sp, "EJ.1 spot check: R^hi(x) >= R_k(x) (exact arb affine recursion) at x = 1e-3..1e6, k0 = 2, 6, 14, k = k0..k0+20")
    # ---- EJ.3: exact sign-change level at 0.9 (and 0.905)
    RESULT["exact"] = {}
    for a_s in EXACT:
        t1 = time.time()
        P = certify_point(a_s, n_p, n_m)
        up = certify_upper(P["al"], (14,), n=n_ux)
        low = certify_lower(a_s, n_l, 12)
        row = up[14]
        print(f"   alpha = {a_s}: D'_oo in {fmt(P['Dp'])}, |nu| in {fmt(-P['nu'])}, 1-delta_N in "
              + ", ".join(f"N={N}: {fmt(P['om'][N])}" for N in NS) + f" ({time.time() - t1:.1f}s)")
        print(f"     k0 = 14: c_hi <= {float(row['c_hi'].upper()):.6f} (head {float(row['head'].upper()):.1e}, tail {float(row['tail'].upper()):.1e}, "
              f"t >= 1 {float(row['ttail'].upper()):.1e}), U <= {float(row['U'].upper()):.6f}; L(12) >= {float(low['L'][12].lower()):.6f}; "
              + ", ".join(f"kappa_{N} in [{float(low['kappa'][(N, 12)].lower()):.6f}, {float(row['kappa'][N]['hi'].upper()):.6f}]" for N in NS))
        check(P["ok"] and row["cond"], f"EJ.1 alpha = {a_s}: EF's routines valid; (3/rho) beta(p^k0(1)) <= 1 (t >= 1 piece uniform in k >= k0)")
        check(bool(low["L"][12] < row["U"]) and all(bool(low["kappa"][(N, 12)] < row["kappa"][N]["hi"]) for N in NS),
              f"EJ.1 alpha = {a_s}: two-sided: L(12) < U(14) and kappa_lo < kappa_hi (no contradiction)")
        res = {}
        for N in NS:
            on = onset_upper(P, up, N)
            n1 = on[0] if on else None
            win = [nn for nn in range(N, (n1 or 200)) if below_upper(P, low, N, nn)]
            res[N] = {"n1": n1, "win": win}
            last = win[-1] if win else None
            exact = n1 is not None and last == n1 - 1
            contig = bool(win) and win == list(range(win[0], win[-1] + 1))
            print(f"     N = {N}: Y_n > Y_oo for EVERY n >= {n1} (EI.2: {EI_N1[a_s] if N == 1 else '-'}); Y_n < Y_oo certified for n in "
                  f"[{win[0] if win else '-'}, {last}]" + (" -> last level with Y_n < Y_oo is EXACTLY n = %d" % last if exact else
                                                            f" -> last negative level in [{last}, {n1 - 1 if n1 else '?'}]"))
            check(n1 is not None and contig, f"EJ.3 alpha = {a_s}, N = {N}: onset n1 = {n1} certified; negative window contiguous")
        RESULT["exact"][a_s] = {"P": P, "up": row, "low": low, "res": res}
    r9 = RESULT["exact"]["0.9"]["res"][1]
    if full:
        check(r9["n1"] == 33 and r9["win"] and r9["win"][-1] == 32 and r9["win"][0] == 1,
              "EJ.3 alpha = 0.9, N = 1: Y_n < Y_oo for n = 1..32 and Y_n > Y_oo for every n >= 33: exactly one sign change, between 32 and 33")
    else:
        check(r9["n1"] is not None and r9["n1"] <= 35 and r9["win"] and r9["win"][-1] >= 31,
              f"EJ.3 (default grid) alpha = 0.9, N = 1: last negative level in [{r9['win'][-1]}, {r9['n1'] - 1}] (--full: exactly 32)")
    # ---- EJ.4: improved onset levels at EI's exponents (N = 1)
    RESULT["n1"] = {}
    for a_s in ALPHAS:
        t1 = time.time()
        if a_s in RESULT["exact"]:
            n1 = RESULT["exact"][a_s]["res"][1]["n1"]
        else:
            P = certify_point(a_s, n_p, n_m)
            up = certify_upper(P["al"], (6, 10, 14), Ns=(1,), n=n_u)
            on = onset_upper(P, up, 1)
            n1 = on[0] if on else None
        RESULT["n1"][a_s] = n1
        if a_s in RESULT["exact"]:
            win = RESULT["exact"][a_s]["res"][1]["win"]
        else:
            low = certify_lower(a_s, n_l, 12)
            win = [nn for nn in range(1, n1) if below_upper(P, low, 1, nn)] if n1 else []
        RESULT.setdefault("win", {})[a_s] = win
        contig = not win or win == list(range(win[0], win[-1] + 1))
        print(f"   alpha = {a_s}: n1 = {n1} (EI.2: {EI_N1[a_s]}); Y_n < Y_oo certified for n in "
              + (f"[{win[0]}, {win[-1]}] -> last negative level in [{win[-1]}, {n1 - 1}]" if win else "(none)") + f" ({time.time() - t1:.1f}s)")
        check(contig, f"EJ.4 alpha = {a_s}: certified negative window contiguous")
        check(n1 is not None and n1 <= EI_N1[a_s] + (0 if full else 2), f"EJ.4 alpha = {a_s}: Y_n > Y_oo for every n >= {n1} (N = 1), not later than EI's {EI_N1[a_s]}")
    # ---- consistency / attempts to break
    e9 = RESULT["exact"]["0.9"]
    check(bool(e9["low"]["L"][12] <= arb("0.73802")) and bool(e9["up"]["U"] >= arb("0.73737")),
          "consistency: L(12) <= EG.4's e*(0.9) upper end 0.73802 and U(14) >= its lower end 0.73737 (the bracket contains e*)")
    nu9 = -e9["P"]["nu"]
    check(bool(nu9.lower() <= 0.432350) and bool(nu9.upper() >= 0.432534), "consistency: |nu_0.9| bracket contains EF.4's [0.432350, 0.432534]")
    # ---- EJ.5 floats
    uf = upper_float(0.9, 14, n=100)
    nf, c = crossing_float(0.9, uf["U"], uf["kappa"][1])
    print(f"   floats alpha = 0.9, k0 = 14: U ~ {uf['U']:.5f}, kappa_hi(1) ~ {uf['kappa'][1]:.5f} (EE: e* ~ {c['e']:.5f}, kappa_1 ~ {c['kappa'][1]:.5f}); float crossing {nf}")
    check(c["e"] < uf["U"] < 1.02 * c["e"] and c["kappa"][1] < uf["kappa"][1] < 1.02 * c["kappa"][1],
          "EJ.5 (OBSERVED): the float upper integrals sit within 2% above EE's limits e*, kappa_1")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EJ checks passed")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    sys.exit(main(full="--full" in sys.argv))
