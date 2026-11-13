#!/usr/bin/env python3
"""Theorem EA: the second-order term of the corner deficit above log3/log5.  rho^N delta_N converges (BJ.6's limit, without
AT.1), strictly from above, and the correction has the exact order (rho/5)^N with an explicit limit constant; so
rho^N delta_N is eventually strictly decreasing for every alpha in (gamma, 1) (BJ's open item), and DZ's limit
C_delta(alpha, N) = -nu_alpha (1 - delta_N) reaches -nu_alpha as  -nu_alpha - C_delta = -nu_alpha (c_corner rho^-N + K 5^-N + o(5^-N)).

Setting (GASKET.md, Theorems BH, BJ, BR, BS, CB, DZ): alpha in (gamma, 1), gamma = log3/log5, rho = 5^alpha/3 > 1,
p(z) = z(z+5), beta(z) = (z+6)/((z+2)(z+5)), eta_N = alpha_1 + beta eta_(N-1) o p, eta_0 = -1/2, eta = lim eta_N,
psi(v) = lim p^k(5^-k v), x_j = psi(5^-j u), R(u) = lim (3/5)^M r_M(5^-M u), Q_oo(u) = lim (3/5)^M r_oo(5^-M u) (BS.1),
Pi(u) = prod_(j>=1) (5/3) beta(x_j), Phi(u) = (1 + 2 eta(x_0)) Pi(u), h = Q_oo - R = Q_oo^2 Phi/(1 + Q_oo Phi) (BS.1-BS.2),
kappa = 2 sin(pi alpha)/(pi A_alpha), delta_N = kappa int tau^-alpha R(tau) (1 - eta_N(psi(5^-N tau))) dtau (BJ.4).
New symbols: Xi(t) = prod_(j>=1) (x_j + 6)/(3 (x_j + 2)) in (0, 1] (x_j = psi(5^-j t)), so Pi(t) psi(t) = t Xi(t);
pi_N(z) = prod_(i<N) beta(p^i z) (1/2 + eta(p^N z)) = eta(z) - eta_N(z) (BH.2's remainder: 1 - eta_N = 1 - eta + pi_N);
b_N(u) = prod_(k=1..N) (5/3) beta(x_k) (1/2 + eta(x_0)), decreasing in N to Phi(u)/2.

Theorem EA (checked here; exit 1 on any failure):
  EA.1 (THEOREM, exact; every alpha in (gamma,1), N >= 0).  With c_corner = kappa int t^-alpha Q_oo(t) (1 - eta(psi(t))) dt,
         rho^N delta_N = c_corner + E_N,  E_N = kappa int_0^oo t^(1-alpha) Q_oo^2 pi_N(psi) Xi / (2 (1 + 2 Q_oo pi_N(psi) Pi)) dt > 0.
       (BJ.4 + BR.4's scaling, BS.1's 1/Q_N = 1/Q_oo + (3/5)^N Phi(5^N t) = 1/Q_oo + 2 pi_N(psi(t)) Pi(t), BH.2, and CB.3:
       1/Q_oo = Pi (2 (1 - eta(psi)) + psi/2).)  So lim rho^N delta_N = c_corner exists without AT.1, and rho^N delta_N >
       c_corner for every N (strictly from above).
  EA.2 (THEOREM; exact rate).  E_N = kappa (rho/5)^N I_N,  I_N = int u^(1-alpha) Q_oo^2 b_N Xi(5^-N u) / (2 (1 + 2 Q_oo b_N
       Pi(5^-N u))) du  ->  K_corner/kappa,  K_corner = (kappa/4) int_0^oo u^(1-alpha) h(u) du in (0, oo);  for every N >= N0:
       K^-_(N0) <= kappa I_N <= K^+_(N0), K^-_(N0) = (kappa/4) int u^(1-alpha) h Xi(5^-N0 u), K^+_(N0) = (kappa/2) int u^(1-alpha)
       Q^2 b_N0 / (1 + 2 Q b_N0 Pi(5^-N0 u)), both -> K_corner.  So delta_N = c_corner rho^-N + K_corner 5^-N + o(5^-N), and since
       K^+_(N0)/K^-_(N0) -> 1 < 5/rho, rho^N delta_N is eventually strictly decreasing (BJ's open item), for every alpha.
  EA.3 (COROLLARY with DZ.2/DZ.3).  -nu_alpha - C_delta(alpha, N) = -nu_alpha delta_N = -nu_alpha rho^-N (c_corner + E_N).
  EA.4 (certified at alpha = 0.8; arb, CB.3's product for Q_oo, cells on [1e-15, 5^8] with closed heads and tails):
       K_corner(0.8), K^+-(N0), I_N; every-N brackets of rho^N delta_N(0.8) and of C_delta(0.8, N) (meet DZ.4's);
       (rho/5) K^+_1 < K^-_1, so rho^N delta_N(0.8) strictly decreases for every N >= 0.
  Floats (consistency only): direct quadrature of delta_N against EA.1 (N = 1..10), E_N (5/rho)^N decreasing toward K.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp
from flint import arb, ctx

FAILS: list[str] = []
RESULT: dict = {}

# stated inputs
A08 = ("2.3964097039792205270396", "2.3964097039792205270397")       # AZ (3)
NU08 = ("-0.7274453844421553283991", "-0.7274453844421553283990")    # AZ (3)
GLO, GHI = "0.600945014", "0.600963487"                              # BJ.3 (Q_oo in [g_lo, g_hi] u^(gamma-1), BS.1)
CC08 = ("1.135890346", "1.135938783")                                # BJ.6: liminf / limsup rho^N delta_N(0.8)
DZ4 = {1: ("0.030733", "0.031869"), 2: ("0.158436", "0.159369"), 3: ("0.257779", "0.258555"),
       6: ("0.461214", "0.461662"), 10: ("0.602405", "0.602621")}   # DZ.4 certified C_delta(0.8, N)
BJ_FLOAT = {10: 1.135914525, 20: 1.135914308}                       # BJ's floats rho^N delta_N^(100)
NS = (1, 2, 3, 4, 6, 8, 10, 12)
N0S = (1, 3, 6)
KORB = 40                    # backward orbit depth: x_KORB enclosed by BD.4/BJ.6's [v, v e^(v/10)]
U_LO, U_HI_EXP = arb("1e-15"), 8   # cells on [1e-15, 5^8]


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def hull(lo: arb, hi: arb) -> arb:
    m = ((lo + hi) / 2).mid()
    r = max((hi.upper() - m).upper(), (m - lo.lower()).upper())
    return arb(m, r)


def iv(lo: str, hi: str) -> arb:
    return hull(arb(lo), arb(hi))


# ----------------------------------------------------------------- sympy
def symbolic_checks():
    Q, pi_, Pi, om, psi, x, c, b, Xi, t = sp.symbols("Q pi Pi omega psi x c b Xi t", positive=True)
    al, N = sp.symbols("alpha N", positive=True)
    # EA.1 algebra: 1/Q_N = 1/Q + 2 pi Pi, 1/Q = Pi (2 omega + psi/2), omega = 1 - eta
    QN = 1 / (1 / Q + 2 * pi_ * Pi)
    lhs = QN * (om + pi_) - Q * om
    rhs = Q**2 * pi_ * Pi * psi / (2 * (1 + 2 * Q * pi_ * Pi))
    ok1 = sp.simplify((lhs - rhs).subs(om, (1 / Q / Pi - psi / 2) / 2)) == 0
    check(ok1, "sympy: with 1/Q = Pi (2(1-eta) + psi/2): Q_N (1-eta+pi_N) - Q (1-eta) = Q^2 pi_N Pi psi / (2 (1 + 2 Q pi_N Pi))")
    # 1 - eta = z (P - 1)/4 from BH.2 (1/r = 2 + z - 2 eta) and CB.2 (1/r = z (1 + P)/2)
    z, P, eta = sp.symbols("z P eta")
    sol = sp.solve(sp.Eq(2 + z - 2 * eta, z * (1 + P) / 2), eta)[0]
    check(sp.simplify(1 - sol - z * (P - 1) / 4) == 0, "sympy: BH.2 + CB.2 give 1 - eta(z) = z (P(z) - 1)/4")
    # CB.3 (1/Q = Pi psi P(psi)/2) in the form Pi (2 (1 - eta(psi)) + psi/2)
    check(sp.simplify(Pi * psi * P / 2 - Pi * (2 * (psi * (P - 1) / 4) + psi / 2)) == 0,
          "sympy: CB.3's Pi psi P(psi)/2 = Pi (2 (1 - eta(psi)) + psi/2)")
    bet = (x + 6) / ((x + 2) * (x + 5))
    check(sp.simplify(sp.Rational(5, 3) * bet * (x + 5) / 5 - (x + 6) / (3 * (x + 2))) == 0,
          "sympy: (5/3) beta(x) (x+5)/5 = (x+6)/(3(x+2)), so Pi(t) psi(t) = t Xi(t) (psi(t)/t = prod (x_j+5)/5)")
    check(sp.simplify(1 - (x + 6) / (3 * x + 6) - 2 * x / (3 * x + 6)) == 0 and
          sp.simplify(1 - sp.Rational(5, 3) * bet - x * (3 * x + 16) / (3 * x**2 + 21 * x + 30)) == 0,
          "sympy: 1 - (x+6)/(3x+6) = 2x/(3x+6) <= x/3; 1 - (5/3) beta = x(3x+16)/(3x^2+21x+30) (BS.1)")
    # u-form scaling: t = 5^-N u, Q(5^-N u) = (5/3)^N Q(u), pi_N(psi(5^-N u)) = (3/5)^N b_N(u)
    rho = 5**al / 3
    fac = 5 ** (-N * (2 - al)) * sp.Rational(5, 3) ** (2 * N) * sp.Rational(3, 5) ** N
    check(sp.simplify(sp.expand_power_base(sp.powsimp(fac / (rho / 5) ** N, force=True), force=True)) == 1
          or abs(float((fac / (rho / 5) ** N).subs({al: 0.8, N: 7})) - 1) < 1e-12,
          "sympy: 5^(-N(2-alpha)) (5/3)^(2N) (3/5)^N = (rho/5)^N")
    # limit integrand: Q^2 (Phi/2) / (2 (1 + Q Phi)) = h/4 with h = Q^2 Phi/(1 + Q Phi)
    Ph = sp.symbols("Phi", positive=True)
    check(sp.simplify(Q**2 * (Ph / 2) / (2 * (1 + 2 * Q * (Ph / 2))) - Q**2 * Ph / (1 + Q * Ph) / 4) == 0,
          "sympy: lim integrand Q^2 b / (2(1 + 2 Q b)) at b = Phi/2 is h/4")
    # monotonicity of G = Q b Xi / (2 (1 + 2 Q b Pi))
    G = Q * b * Xi / (2 * (1 + 2 * Q * b * Pi))
    dQ, db, dX, dP = (sp.simplify(sp.diff(G, v)) for v in (Q, b, Xi, Pi))
    okm = all(sp.simplify(d.subs({Q: 0.7, b: 0.3, Xi: 0.9, Pi: 0.5})) > 0 for d in (dQ, db, dX)) and \
        sp.simplify(dP.subs({Q: 0.7, b: 0.3, Xi: 0.9, Pi: 0.5})) < 0
    okm &= sp.factor(dQ) == sp.factor(b * Xi / (2 * (1 + 2 * Q * b * Pi) ** 2)) and sp.factor(db) == sp.factor(Q * Xi / (2 * (1 + 2 * Q * b * Pi) ** 2))
    check(okm, "sympy: G = Q b Xi/(2(1+2 Q b Pi)) has dG/dQ = b Xi/(2(1+2QbPi)^2) > 0, dG/db > 0, dG/dXi > 0, dG/dPi < 0")
    # integrability exponents for alpha in (gamma, 1): 2gamma - 1 - alpha > -1 and 2gamma - 2 - alpha < -1, gamma - alpha > -1
    gm = math.log(3) / math.log(5)
    okx = all(2 * gm - 1 - a > -1 and 2 * gm - 2 - a < -1 and gm - a > -1 for a in np.linspace(gm + 1e-9, 1, 200))
    check(okx and 2 * gm > 1, "exponents: u^(2gamma-1-alpha) integrable at 0 and u^(2gamma-2-alpha) at oo for every alpha in (gamma,1) (alpha < 2gamma, alpha > 2gamma-1)")


# ----------------------------------------------------------------- arb point enclosures
def _p(x):
    return x * (x + 5)


def _beta(x):
    return (x + 6) / ((x + 2) * (x + 5))


def orbit(u):
    """x_j = psi(5^-j u), j = 0..KORB (BD.4/BJ.6: psi(v) in [v, v e^(v/10)]), and sigma >= sum_(j>KORB) x_j."""
    v = u / arb(5) ** KORB
    vu = arb(v.upper())
    x = [None] * (KORB + 1)
    x[KORB] = hull(v, v * (vu / 10).exp())
    for j in range(KORB, 0, -1):
        x[j - 1] = _p(x[j])
    sig = arb(((vu / 4) * (vu / 10).exp()).upper())
    return x, sig


def one_minus_eta(z):
    """1 - eta(z) = z (P(z) - 1)/4 (CB.2 + BH.2), P - 1 = expm1(sum log1p(4/(w_j + 2))), at least two factors, tail in [0, 8/w_m]."""
    L = arb(0)
    w = z
    big = arb("1e30")
    n = 0
    while n < 2 or not bool(w.lower() > big):
        L += (4 / (w + 2)).log1p()
        w = _p(w)
        n += 1
        if n > 400:
            raise RuntimeError("orbit does not escape")
    L += hull(arb(0), arb((8 / w).upper()))         # w_(j+1) >= w_j^2: sum_(j>=m) log1p(4/(w_j+2)) in [0, 8/w_m]
    return z * L.expm1() / 4


def point(u):
    x, sig = orbit(u)
    f3 = arb(5) / 3
    fb = [None] + [f3 * _beta(x[j]) for j in range(1, KORB + 1)]
    fx = [None] + [(x[j] + 6) / (3 * x[j] + 6) for j in range(1, KORB + 1)]
    tpi = hull(1 - arb("0.54") * sig, arb(1))           # prod_(j>KORB) (5/3) beta(x_j): 1 - factor <= x(3x+16)/30 <= 0.54 x
    txi = hull(1 - sig / 3, arb(1))                     # prod_(j>KORB) (x_j+6)/(3x_j+6): 1 - factor <= x/3
    # suffix products: Pi(5^-N u) = prod_(j>N) fb_j, Xi(5^-N u) = prod_(j>N) fx_j
    sufP = [None] * (KORB + 1)
    sufX = [None] * (KORB + 1)
    sufP[KORB], sufX[KORB] = tpi, txi
    for j in range(KORB, 0, -1):
        sufP[j - 1] = sufP[j] * fb[j]
        sufX[j - 1] = sufX[j] * fx[j]
    om = one_minus_eta(x[0])
    eta0 = 1 - om
    Pi = sufP[0]
    Q = 1 / (Pi * x[0] / 2 * (1 + 4 * om / x[0]))      # CB.3: 1/Q = Pi x0 P(x0)/2, P = 1 + 4 (1-eta)/x0
    Phi = (1 + 2 * eta0) * Pi
    pre = [None] * (KORB + 1)                            # prod_(k=1..N) fb_k
    pre[0] = arb(1)
    for k in range(1, KORB + 1):
        pre[k] = pre[k - 1] * fb[k]
    b = {N: pre[N] * (arb(1) / 2 + eta0) for N in set(NS) | set(N0S)}
    return {"x": x, "Q": Q, "Phi": Phi, "Pi": Pi, "om": om, "eta0": eta0, "b": b,
            "PiN": {N: sufP[N] for N in set(NS) | set(N0S)}, "XiN": {N: sufX[N] for N in set(NS) | set(N0S)},
            "Xi0": sufX[0]}


def grid():
    """geometric grid on [1e-15, 5^8]: ln-step 0.0004 on [1e-3, 5^8] and 0.004 below (where Q Phi >> 1)."""
    pts = []
    a, m, b = math.log(1e-15), math.log(1e-3), U_HI_EXP * math.log(5)
    n1 = int(math.ceil((m - a) / 0.004))
    n2 = int(math.ceil((b - m) / 0.0004))
    for i in range(n1):
        pts.append(a + (m - a) * i / n1)
    for i in range(n2 + 1):
        pts.append(m + (b - m) * i / n2)
    us = [U_LO] + [arb(math.exp(s)) for s in pts[1:-1]] + [arb(5) ** U_HI_EXP]
    return us


def certify(al: arb):
    gm = arb(3).log() / arb(5).log()
    glo, ghi = arb(GLO), arb(GHI)
    us = grid()
    t1 = time.time()
    P = [point(u) for u in us]
    tev = time.time() - t1
    e1 = 1 + gm - al
    w = [(us[i + 1] ** e1 - us[i] ** e1) / e1 for i in range(len(us) - 1)]   # exact int u^(gamma-alpha) over the cell

    def sig(y):
        return y / (1 + y)

    def G(Q, b, Xi, Pi):
        return Q * b * Xi / (2 * (1 + 2 * Q * b * Pi))

    acc = {"K": [arb(0), arb(0)]}
    for N in NS:
        acc[("I", N)] = [arb(0), arb(0)]
    for N0 in N0S:
        acc[("Km", N0)] = [arb(0), arb(0)]
        acc[("Kp", N0)] = [arb(0), arb(0)]
    for i in range(len(us) - 1):
        u0, u1, A, B = us[i], us[i + 1], P[i], P[i + 1]
        # Q(u) in u^(gamma-1) [Flo, Fhi] on the cell (Q nonincreasing, BS.1's every-u bracket)
        Fhi = min(ghi.upper(), (u1 ** (1 - gm) * A["Q"]).upper())
        Flo = max(glo.lower(), (u0 ** (1 - gm) * B["Q"]).lower())
        Fhi, Flo = arb(Fhi), arb(Flo)
        Qhi, Qlo = arb(A["Q"].upper()), arb(B["Q"].lower())
        wl, wh = Flo * w[i], Fhi * w[i]
        # K: u^(1-a) h = u^(1-a) Q sig(Q Phi), sig increasing; Q, Phi nonincreasing in u
        acc["K"][0] += wl * sig(Qlo * arb(B["Phi"].lower()))
        acc["K"][1] += wh * sig(Qhi * arb(A["Phi"].upper()))
        for N in NS:
            lo = G(Qlo, arb(B["b"][N].lower()), arb(B["XiN"][N].lower()), arb(A["PiN"][N].upper()))
            hi = G(Qhi, arb(A["b"][N].upper()), arb(A["XiN"][N].upper()), arb(B["PiN"][N].lower()))
            acc[("I", N)][0] += wl * arb(lo.lower())
            acc[("I", N)][1] += wh * arb(hi.upper())
        for N0 in N0S:
            acc[("Km", N0)][0] += wl * arb((sig(Qlo * arb(B["Phi"].lower())) * B["XiN"][N0]).lower())
            acc[("Km", N0)][1] += wh * arb((sig(Qhi * arb(A["Phi"].upper())) * A["XiN"][N0]).upper())
            lo = G(Qlo, arb(B["b"][N0].lower()), arb(1), arb(A["PiN"][N0].upper()))
            hi = G(Qhi, arb(A["b"][N0].upper()), arb(1), arb(B["PiN"][N0].lower()))
            acc[("Kp", N0)][0] += wl * arb(lo.lower())
            acc[("Kp", N0)][1] += wh * arb(hi.upper())
    # heads (0, 1e-15]: h <= Q <= g_hi u^(gamma-1); G <= Xi/(4 Pi) <= 1/(4 Pi(1e-15)) (Pi nonincreasing)
    eps = us[0]
    headK = ghi * eps ** e1 / e1
    headG = headK / (4 * arb(P[0]["Pi"].lower()))
    # tails [T, oo), T = 5^8: Phi <= 3 Pi <= 5/x_1, b_N <= (5/2)/x_1, x_1(u) >= x_1(T) u/T (psi(v)/v nondecreasing)
    T = us[-1]
    x1T = arb(P[-1]["x"][1].lower())
    e2 = 1 + al - 2 * gm
    tailK = 5 * ghi**2 * T / x1T * T ** (-e2) / e2
    tailG = arb(5) / 4 * ghi**2 * T / x1T * T ** (-e2) / e2
    out = {"K": hull(acc["K"][0], acc["K"][1] + headK + tailK)}
    for N in NS:
        out[("I", N)] = hull(acc[("I", N)][0], acc[("I", N)][1] + headG + tailG)
    for N0 in N0S:
        out[("Km", N0)] = hull(acc[("Km", N0)][0], acc[("Km", N0)][1] + headK + tailK)
        out[("Kp", N0)] = hull(acc[("Kp", N0)][0], acc[("Kp", N0)][1] + headG + tailG)
    out["ncell"] = len(us) - 1
    out["tev"] = tev
    out["heads"] = (headK, headG)
    out["tails"] = (tailK, tailG)
    out["F"] = (min(float((p["Q"] * u ** (1 - gm)).lower()) for p, u in zip(P, us)),
                max(float((p["Q"] * u ** (1 - gm)).upper()) for p, u in zip(P, us)))
    return out


def identity_spot_checks(al: arb):
    """EA.1 pointwise, independently: Q_N(t) = (5/3)^N R(5^N t) from a separate orbit at 5^N t (BS.1: R = 1/(1/Q + Phi)),
    eta_N(psi(t)) by the finite recursion from eta_0 = -1/2; against the closed correction integrand."""
    ok = True
    worst = 0.0
    for ts in ("0.001", "0.37", "2.5", "40"):
        t = arb(ts)
        A = point(t)
        x0 = A["x"][0]
        for N in (1, 3, 6):
            B = point(arb(5) ** N * t)
            R = 1 / (1 / B["Q"] + B["Phi"])
            QN = (arb(5) / 3) ** N * R
            # eta_N(x0) by recursion
            orb = [x0]
            for _ in range(N):
                orb.append(_p(orb[-1]))
            e = arb(-1) / 2
            for i in range(N - 1, -1, -1):
                z = orb[i]
                e = (z + 4) / ((z + 2) * (z + 5)) + _beta(z) * e
            lhs = t ** (-al) * (QN * (1 - e) - A["Q"] * A["om"])
            prodb = arb(1)
            for i in range(N):
                prodb *= _beta(orb[i])
            piN = prodb * (arb(1) / 2 + 1 - one_minus_eta(orb[N]))
            rhs = t ** (1 - al) * A["Q"] ** 2 * piN * A["Xi0"] / (2 * (1 + 2 * A["Q"] * piN * A["Pi"]))
            d = abs(lhs - rhs)
            worst = max(worst, float(d.upper()))
            ok &= bool(d.upper() < 1e-20 + 1e-12 * abs(rhs).upper()) and bool(rhs.lower() > 0)
            # (3/5)^N Phi(5^N t) = 2 pi_N(psi(t)) Pi(t)
            ok &= bool(abs((arb(3) / 5) ** N * B["Phi"] - 2 * piN * A["Pi"]).upper() < 1e-12 * abs(B["Phi"]).upper() * (arb(3) / 5) ** N)
            # Pi psi = t Xi
            ok &= bool(abs(A["Pi"] * x0 - t * A["Xi0"]).upper() < 1e-15 * abs(t * A["Xi0"]).upper())
    check(ok, f"EA.1 pointwise at t = 0.001, 0.37, 2.5, 40, N = 1, 3, 6 (Q_N from an independent orbit at 5^N t, eta_N by recursion): "
              f"LHS = RHS > 0 (|diff| <= {worst:.1e}); (3/5)^N Phi(5^N t) = 2 pi_N Pi; Pi psi = t Xi")
    return worst


# ----------------------------------------------------------------- floats
def floats(al=0.8, Nmax=10):
    gm = math.log(3) / math.log(5)
    rho = 5**al / 3
    kap = 2 * math.sin(math.pi * al) / (math.pi * float(A08[0])) if al == 0.8 else None
    J = 40

    def p(z):
        return z * (z + 5)

    def beta(z):
        return (z + 6) / ((z + 2) * (z + 5))

    def Pm1(z):
        S = np.zeros_like(z)
        w = np.minimum(z, 1e150)
        for _ in range(90):
            m = w < 1e150
            S = np.where(m, S + np.log1p(4 / (w + 2)), S)
            w = np.where(m, p(np.minimum(w, 1e150)), w)
        return np.expm1(S)

    def ometa(z):
        zz = np.minimum(z, 1e100)
        return np.where(z > 1e100, 1.0, zz * Pm1(zz) / 4)

    def base0(u):
        x = [None] * (J + 1)
        x[J] = u / 5.0**J
        for j in range(J, 0, -1):
            x[j - 1] = p(x[j])
        Pi = np.ones_like(u)
        Xi = np.ones_like(u)
        for j in range(1, J + 1):
            Pi *= 5 / 3 * beta(x[j])
            Xi *= (x[j] + 6) / (3 * x[j] + 6)
        z = x[0]
        om = ometa(z)
        iQ = Pi * z / 2 * (1 + 4 * om / z)
        Q = 1 / iQ
        Phi = (3 - 2 * om) * Pi
        return dict(x0=z, Pi=Pi, Xi=Xi, Q=Q, om=om, Phi=Phi, R=1 / (iQ + Phi), h=Q * Q * Phi / (1 + Q * Phi))

    def base(u):
        big = u > 5.0**7
        m = np.where(big, np.ceil(np.log(np.maximum(u, 1) / 5.0**6) / np.log(5)), 0)
        o = base0(u / 5.0**m)
        sc = (3 / 5) ** m
        Q = o["Q"] * sc
        o["Q"] = Q
        o["R"] = np.where(big, Q, o["R"] * sc)
        o["h"] = np.where(big, 0.0, o["h"] * sc)
        o["om"] = np.where(big, 1.0, o["om"])
        o["Pi"] = np.where(big, 0.0, o["Pi"])
        o["Xi"] = np.where(big, 0.0, o["Xi"])
        o["x0"] = np.where(big, 1e300, o["x0"])
        return o

    def piN(N, z):
        pr = np.ones_like(z)
        w = z.copy()
        for _ in range(N):
            pr *= beta(w)
            w = p(np.minimum(w, 1e150))
        return pr * (1.5 - ometa(w))

    gx, gw = np.polynomial.legendre.leggauss(30)

    def LI(f, a, b, step=np.log(5) / 8):
        n = int(np.ceil((b - a) / step))
        e = np.linspace(a, b, n + 1)
        lo, hi = e[:-1][:, None], e[1:][:, None]
        xx = (lo + hi) / 2 + (hi - lo) / 2 * gx[None, :]
        ww = (hi - lo) / 2 * gw[None, :]
        t = np.exp(xx).ravel()
        return float(np.sum(f(t) * t * ww.ravel()))

    def full(f, N):
        T = 5.0 ** (N + 14)
        return LI(f, -150.0, np.log(T)) + LI(lambda u: u ** (-al) * base(u)["Q"], np.log(T), np.log(5 * T)) / (1 - 1 / rho)

    with np.errstate(all="ignore"):
        kR = kap * full(lambda u: u ** (-al) * base(u)["R"], 0)
        c = kap * full(lambda t: t ** (-al) * base(t)["Q"] * base(t)["om"], 0)
        K = kap / 4 * LI(lambda u: u ** (1 - al) * base(u)["h"], -150.0, np.log(5.0**25))
        rows = {}
        for N in range(1, Nmax + 1):
            def fd(u, N=N):
                B = base(u)
                z = base(u / 5.0**N)["x0"]
                return u ** (-al) * B["R"] * (ometa(z) + piN(N, z))
            d = kap * full(fd, N)

            def fe(t, N=N):
                B = base(t)
                pn = piN(N, B["x0"])
                return t ** (1 - al) * B["Q"] ** 2 * pn * B["Xi"] / (2 * (1 + 2 * B["Q"] * pn * B["Pi"]))
            E = kap * LI(fe, -150.0, np.log(5.0**20))
            rows[N] = (rho**N * d, rho**N * d - c, E, E / (rho / 5) ** N)
    return {"kR": kR, "c": c, "K": K, "rows": rows}


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem EA: rho^N delta_N = c_corner + E_N, E_N > 0 exactly, E_N (5/rho)^N -> K_corner; rho^N delta_N eventually "
          "strictly decreasing (BJ's open item); -nu - C_delta(alpha,N) = -nu delta_N. Not thrust.")
    symbolic_checks()
    ctx.prec = 96
    al = arb("0.8")
    gm = arb(3).log() / arb(5).log()
    rho = arb(5) ** al / 3
    kap = 2 * (arb.pi() * al).sin() / (arb.pi() * iv(*A08))
    nu = iv(*NU08)
    worst = identity_spot_checks(al)
    ctx.prec = 64
    cert = certify(al)
    print(f"   arb: {cert['ncell']} cells on [1e-15, 5^8], {len(NS) + 2 * len(N0S) + 1} integrals, point evaluations {cert['tev']:.1f}s; "
          f"heads <= {float(cert['heads'][1].upper()):.1e}, tails <= {float(cert['tails'][0].upper()):.1e}")
    Flo, Fhi = cert["F"]
    check(Flo >= float(GLO) - 1e-9 and Fhi <= float(GHI) + 1e-9,
          f"consistency: CB.3's Q_oo u^(1-gamma) at the grid points lies in BJ.3/BS.1's [g_lo, g_hi] ([{Flo:.9f}, {Fhi:.9f}])")
    K = kap / 4 * cert["K"]
    Km = {N0: kap / 4 * cert[("Km", N0)] for N0 in N0S}
    Kp = {N0: kap * cert[("Kp", N0)] for N0 in N0S}
    I = {N: kap * cert[("I", N)] for N in NS}          # kappa I_N
    print(f"   K_corner(0.8) = (kappa/4) int u^(1-a) h in [{float(K.lower()):.7f}, {float(K.upper()):.7f}]")
    for N0 in N0S:
        print(f"   N0 = {N0}: K^-_N0 in [{float(Km[N0].lower()):.7f}, {float(Km[N0].upper()):.7f}], "
              f"K^+_N0 in [{float(Kp[N0].lower()):.7f}, {float(Kp[N0].upper()):.7f}]")
    okb = all(bool(Km[N0].lower() <= K.upper()) and bool(K.lower() <= Kp[N0].upper()) for N0 in N0S)
    okb &= all(bool(Km[N0].lower() <= I[N].upper()) and bool(I[N].lower() <= Kp[N0].upper()) for N0 in N0S for N in NS if N >= N0)
    check(okb, "EA.2 consistency: K_corner and every certified kappa I_N (N >= N0) lie in [K^-_N0, K^+_N0]")
    check(bool(K.lower() > 0) and bool((K.upper() - K.lower()) < 1e-3 * K.lower()),
          f"EA.4: K_corner(0.8) > 0, relative width < 1e-3")
    mono = bool((rho / 5 * Kp[1]).upper() < Km[1].lower())
    check(mono, f"EA.4: (rho/5) K^+_1 = {float((rho / 5 * Kp[1]).upper()):.5f} < K^-_1 = {float(Km[1].lower()):.5f}: "
                "E_(N+1) < E_N for every N >= 1, so rho^N delta_N(0.8) strictly decreases for N >= 1")
    # every-N brackets
    cc = iv(*CC08)
    tab = {}
    for N in NS:
        E = (rho / 5) ** N * I[N]
        rd = hull(cc.lower() + E.lower(), cc.upper() + E.upper())
        dN = rd / rho**N
        Cd = hull((-nu * (1 - dN.upper())).lower(), (-nu * (1 - dN.lower())).upper())
        tab[N] = (E, rd, dN, Cd)
        print(f"   N = {N:2d}: kappa I_N in [{float(I[N].lower()):.7f}, {float(I[N].upper()):.7f}], E_N in [{float(E.lower()):.4e}, {float(E.upper()):.4e}], "
              f"rho^N delta_N in [{float(rd.lower()):.7f}, {float(rd.upper()):.7f}], C_delta(0.8,N) in [{float(Cd.lower()):.6f}, {float(Cd.upper()):.6f}]")
    check(bool(tab[1][1].upper() < 1.5), "rho^0 delta_0 = 3/2 (BL.0) > rho delta_1: with the line above, rho^N delta_N(0.8) strictly decreases for every N >= 0")
    okdz = all(bool(tab[N][3].lower() <= arb(DZ4[N][1]).upper()) and bool(arb(DZ4[N][0]).lower() <= tab[N][3].upper()) for N in DZ4 if N in tab)
    check(okdz, "EA.4 vs DZ.4 (independent route: level-60 brackets): the C_delta(0.8, N) enclosures meet for N = 1, 2, 3, 6, 10")
    narrower = all(bool((tab[N][3].upper() - tab[N][3].lower()) < (arb(DZ4[N][1]) - arb(DZ4[N][0])).lower()) for N in DZ4 if N in tab)
    check(narrower, "EA.4: the new C_delta(0.8, N) enclosures are narrower than DZ.4's for N = 1, 2, 3, 6, 10")
    lim = -nu * cc
    print(f"   EA.3 at 0.8: rho^N (-nu - C_delta(0.8,N)) -> -nu c_corner in [{float(lim.lower()):.6f}, {float(lim.upper()):.6f}] from above; "
          f"(5/rho)^N [rho^N (-nu - C_delta) + nu c_corner] -> -nu K in [{float((-nu * K).lower()):.6f}, {float((-nu * K).upper()):.6f}]")
    okbj = all(bool(tab[N][1].lower() <= BJ_FLOAT[N] + 1e-7) and bool(BJ_FLOAT[N] - 1e-7 <= tab[N][1].upper()) for N in (10,) if N in tab)
    check(okbj, "BJ's float rho^10 delta_10^(100) = 1.135914525 lies in the certified rho^10 delta_10 (up to 1e-7 for the finite level)")
    check(bool(tab[10][1].upper() < arb("1.197144")) and bool(tab[10][1].lower() > arb("1.135890")),
          "every-N bracket at N = 10 inside BJ.6's ceiling 1.197144 and above BJ.6's liminf floor")
    RESULT["ea"] = {"K": K, "Km": Km, "Kp": Kp, "I": I, "tab": tab, "ncell": cert["ncell"], "mono": mono, "worst": worst, "lim": lim}
    # floats
    fl = floats()
    RESULT["floats"] = fl
    print(f"   floats: kappa int R = {fl['kR']:.12f}, c_corner = {fl['c']:.10f}, K = {fl['K']:.10f}")
    okf = abs(fl["kR"] - 1) < 1e-9 and bool(cc.lower() <= fl["c"] <= cc.upper()) and bool(K.lower() - 1e-9 <= fl["K"] <= K.upper() + 1e-9)
    rows = fl["rows"]
    for N in sorted(rows):
        rd, Ed, Ef, sc = rows[N]
        print(f"     N = {N:2d}: rho^N delta_N = {rd:.12f}, direct rho^N delta_N - c = {Ed:.6e}, EA.1 integral = {Ef:.6e}, E_N (5/rho)^N = {sc:.9f}")
    okf &= all(abs(rows[N][1] - rows[N][2]) < 1e-9 * max(1.0, 0) + 2e-12 + 1e-6 * rows[N][2] for N in rows)
    okf &= all(rows[N + 1][3] < rows[N][3] for N in range(1, max(rows)))
    okf &= all(bool(tab[N][0].lower() <= rows[N][2] <= tab[N][0].upper()) for N in NS if N in rows)
    check(okf, "floats: kappa int R = 1; float c_corner in BJ.6's band and K in EA.4's; direct quadrature of rho^N delta_N - c_corner "
               "= EA.1's integral (N = 1..10); E_N (5/rho)^N decreasing (observed only); float E_N inside the certified brackets")
    dt_ = time.time() - t0
    print(f"Theorem EA: {'all checks passed' if not FAILS else str(len(FAILS)) + ' FAILED'} ({dt_:.1f}s). Functions on a graph, not a flux. Not thrust.")
    return 0 if not FAILS else 1


if __name__ == "__main__":
    sys.exit(main())
