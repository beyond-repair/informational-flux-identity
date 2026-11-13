#!/usr/bin/env python3
"""Theorem EB: the second-order term in n of DZ's renormalised corner deficit.  The level-n deficit is a Moebius function
of rho^-n up to O(5^-0.65 n):  delta_N^(n) = (delta_N + nu_alpha rho^-n)/(1 + nu_alpha rho^-n) + R_n,  so
rho^n (delta_N - delta_N^(n)) = C_delta/(1 + nu rho^-n) + O(rho^n 5^-0.65n)  and the second-order term is
rho^n [rho^n (delta_N - delta_N^(n)) - C_delta(alpha, N)] -> -nu_alpha C_delta(alpha, N) = nu_alpha^2 (1 - delta_N) > 0.

Setting (GASKET.md, Theorems AZ, DZ, EA): alpha in (gamma, 1), gamma = log3/log5, rho = 5^alpha/3 > 1, J_k = int s^-alpha r_k,
D'_k = rho^-k J_k -> D'_oo, AZ (1): D'_k = D'_oo + rho^-k (C' + eps_k), eps_k in [-E_k^-, E_k^+], eps_k = O(rho^k 5^-0.65k),
nu = nu_alpha = C'/D'_oo = C_alpha/A_alpha (< 0 by DZ.3); DZ.1: rho^-k H_k = h_oo - eta_k, 0 <= eta_k <= ebar_k = O(3^-k + (5 rho)^-k),
h_oo = (1 - delta_N) D'_oo; DZ.2: X_k = rho^k (delta_N - delta_N^(k)) = [-(1 - delta_N)(C' + eps_k) - rho^k eta_k] / D'_k,
C_delta = C_delta(alpha, N) = -nu (1 - delta_N).

Theorem EB (checked here; exit 1 on any failure):
  EB.1 (THEOREM, exact; every alpha in (gamma,1), N >= 1, k >= N with 1 + nu rho^-k != 0):
         X_k = C_delta / (1 + nu rho^-k) - [(1 - delta_N) eps_k / (1 + nu rho^-k) + rho^k eta_k] / D'_k,
       equivalently delta_N^(k) = (delta_N + nu rho^-k)/(1 + nu rho^-k) + R_k, R_k = rho^-k [(1-delta_N) eps_k/(1+nu rho^-k)
       + rho^k eta_k] / D'_k, and 1/(1 - delta_N^(k)) = (1 + nu rho^-k)/(1 - delta_N) + O(R_k).
  EB.2 (THEOREM; every alpha in (gamma,1), N >= 1).  R_n = O(5^-0.65n) (AZ (1)'s eps rate, DZ.1's ebar, 5^0.65 < 3), and since
       rho^2 5^-0.65 < 1 on (gamma, 1):  rho^n (X_n - C_delta) -> Y_oo(alpha, N) = -nu C_delta = nu^2 (1 - delta_N) > 0.
       So X_n approaches C_delta strictly from above at the exact rate rho^-n and is eventually strictly decreasing;
       delta_N^(n+1) - delta_N^(n) = C_delta (1 - 1/rho) rho^-n + Y_oo (1 - rho^-2) rho^-2n + o(rho^-2n).
       More: X_n = C_delta sum_(j<J) (-nu rho^-n)^j + o(rho^-(J-1) n) whenever rho^-(J-1) > rho 5^-0.65; in particular the third
       term (rho^n (X_n - C_delta) - Y_oo ~ -nu Y_oo rho^-n > 0, Y_n eventually decreasing) for alpha < alpha* =
       (log_5 27 + 0.65)/3 = 0.89926...
  EB.3 (certified at alpha = 0.8; arb; AZ (3), EA.4, DZ.1 as stated):  for every k >= 60, N in {1,2,3,6,10}:
       -L_N <= X_k - C_delta/(1 + nu rho^-k) <= U_N with L_N <= 6.7e-21, U_N <= 8.7e-35; Y_oo(0.8, N) = -nu C_delta brackets;
       X_k > C_delta and X_k > X_(k+1) for every k in [60, K_N]; Y_k > Y_oo and Y_k > Y_(k+1) for k in [60, K'_N].
  Floats (OBSERVED): nu_k(N) := rho^k (u_k - u_(k+1)) / ((1 - 1/rho) u_oo), u_k = 1/(1 - delta_N^(k)) (u_oo from u_K via EB.1), recovers the
       N-independent nu_alpha (0.7, 0.8, 0.9, 0.95; N = 1, 2, 6, 10), with error ratio ~ rho/3 per level (R_n ~ 3^-n observed).

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time

import numpy as np

FAILS: list[str] = []
RESULT: dict = {}

# stated inputs
A08 = ("2.3964097039792205270396", "2.3964097039792205270397")        # AZ (3)
C08 = ("-1.7432571783920757239548", "-1.7432571783920757239546")      # AZ (3)
NU08 = ("-0.7274453844421553283991", "-0.7274453844421553283990")     # AZ (3)
E60P, E60M = "7.22e-23", "6.7e-34"                                    # AZ (3)
G_BI4 = "0.6031"                                                      # BI.4 (re-run by DZ)
CD_EA = {1: ("0.031293", "0.031327"), 2: ("0.158926", "0.158951"), 3: ("0.258205", "0.258225"),   # EA.4: C_delta(0.8, N)
         6: ("0.461475", "0.461487"), 10: ("0.602532", "0.602538")}
DZ4 = {1: ("0.030733", "0.031869"), 2: ("0.158436", "0.159369"), 3: ("0.257779", "0.258555"),    # DZ.4 (level 60)
       6: ("0.461214", "0.461662"), 10: ("0.602405", "0.602621")}
NS = (1, 2, 3, 6, 10)
KMAX = 400


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- floats (numpy only; same recursion as DZ's script)
def _rn(s, n):
    zs = [s]
    for _ in range(n):
        zs.append(np.minimum(zs[-1] * (zs[-1] + 5), 1e150))
    r = 1 / (3 + zs[n])
    for j in range(n - 1, -1, -1):
        y = zs[j]
        r = r / ((y + 6) / ((y + 2) * (y + 5)) - 2 * y / (2 + y) * r)
        r = np.where(y > 1e60, 1 / (y + 2), r)
    return r, zs


def _eta(zs, n, N):
    e = -0.5 * np.ones_like(zs[0])
    for k in range(n - 1, n - N - 1, -1):
        zz = zs[k]
        D = (zz + 2) * (zz + 5)
        e = (zz + 4) / D + (zz + 6) / D * e
    return e


def float_deltas(al, n, Ns, h=2e-3):
    """delta_N^(n) by trapezoid in log s (DZ.5's quadrature)."""
    with np.errstate(all="ignore"):
        u = np.arange(np.log(5.0 ** (-n - 40)), np.log(1e12), h)
        s = np.exp(u)
        r, zs = _rn(s, n)
        wgt = s ** (1 - al) * r
        J = np.trapezoid(wgt, u)
        return {N: np.trapezoid(wgt * (1 - _eta(zs, n, N)), u) / J for N in Ns if N <= n}


def nu_from_levels(al, Ns, ks, K):
    """OBSERVED: nu_k(N) = rho^k (u_k - u_(k+1)) / ((1 - 1/rho) u_K), u_k = 1/(1 - delta_N^(k)); EB.1 => -> nu_alpha."""
    rho = 5**al / 3
    lv = sorted(set(ks) | {k + 1 for k in ks} | {K})
    d = {k: float_deltas(al, k, Ns) for k in lv}
    out = {}
    for N in Ns:
        uK = 1 / (1 - d[K][N])
        # u_K = u_oo (1 + nu rho^-K) by EB.1, so with c = q_k/u_K: nu ~ c / (1 - c rho^-K)  (matters only at alpha = 0.7)
        rows = []
        for k in ks:
            if k >= N:
                c = rho**k * (1 / (1 - d[k][N]) - 1 / (1 - d[k + 1][N])) / ((1 - 1 / rho) * uK)
                rows.append((k, c / (1 - c * rho ** (-K))))
        out[N] = rows
    return out


# ----------------------------------------------------------------- (A) sympy
def symbolic_checks():
    import sympy as sp
    d, e, eta, x, nu, Dinf, rk, a = sp.symbols("delta eps eta x nu Dinf rk alpha")
    Cp = nu * Dinf
    Dk = Dinf + x * (Cp + e)
    Cd = -nu * (1 - d)
    X = (-(1 - d) * (Cp + e) - rk * eta) / Dk                       # DZ.2's identity with h_oo = (1 - delta) D'_oo
    h = (1 - d) * Dinf
    X_dz = (-h * (Cp + e) - rk * eta * Dinf) / (Dk * Dinf)
    check(sp.simplify(X - X_dz) == 0, "DZ.2's X_k = [-h(C'+eps) - rho^k eta D'_oo]/(D'_k D'_oo) = [-(1-delta)(C'+eps) - rho^k eta]/D'_k")
    eb1 = Cd / (1 + nu * x) - ((1 - d) * e / (1 + nu * x) + rk * eta) / Dk
    check(sp.simplify(X - eb1) == 0, "EB.1: X_k = C_delta/(1+nu x) - [(1-delta) eps/(1+nu x) + rho^k eta]/D'_k exactly (x = rho^-k; D'_k = D'_oo(1+nu x) + x eps)")
    mob = (d + nu * x) / (1 + nu * x)
    check(sp.simplify(mob - (d - x * Cd / (1 + nu * x))) == 0 and sp.simplify(1 / (1 - mob) - (1 + nu * x) / (1 - d)) == 0,
          "EB.1: delta - x C_delta/(1+nu x) = (delta + nu x)/(1 + nu x), and 1/(1 - that) = (1 + nu x)/(1 - delta) (linear in rho^-k)")
    ser = sp.series(Cd / (1 + nu * x), x, 0, 4).removeO()
    check(sp.expand(ser - Cd * (1 - nu * x + nu**2 * x**2 - nu**3 * x**3)) == 0 and sp.simplify(-nu * Cd - nu**2 * (1 - d)) == 0,
          "C_delta/(1+nu x) = C_delta (1 - nu x + nu^2 x^2 - ...); second-order coefficient -nu C_delta = nu^2 (1 - delta)")
    r = sp.symbols("rho", positive=True)
    W = (1 - (1 / r + r**-2) + r**-3) - (1 - 1 / r) * (1 - r**-2)
    check(sp.simplify(W) == 0, "1 - rho^-1 - rho^-2 + rho^-3 = (1 - rho^-1)(1 - rho^-2) (second-difference normalisation)")
    # increments: X_k = C + Y x^k => rho^k (delta^(k+1) - delta^(k)) = C(1 - 1/rho) + Y (1 - rho^-2) rho^-k
    C, Y, k = sp.symbols("C Y k")
    Xk = lambda j: C + Y * r ** (-j)
    inc = r**k * (r ** (-k) * Xk(k) - r ** (-k - 1) * Xk(k + 1))
    check(sp.simplify(inc - (C * (1 - 1 / r) + Y * (1 - r**-2) * r ** (-k))) == 0,
          "delta^(k+1) - delta^(k) = rho^-k X_k - rho^-(k+1) X_(k+1) = C(1 - 1/rho) rho^-k + Y (1 - rho^-2) rho^-2k (+ higher)")
    # exponent conditions
    l5 = lambda v: sp.log(v) / sp.log(5)
    ast = (l5(27) + sp.Rational(13, 20)) / 3
    conds = {
        "5^0.65 < 3 (so 3^-n = o(5^-0.65n))": sp.N(sp.Rational(13, 20) - l5(3)) < 0,
        "rho^2 5^-0.65 < 1 for every alpha < 1.0076 (alpha <= 1: 2 - 0.65 < log_5 9)": sp.N(2 - sp.Rational(13, 20) - l5(9)) < 0,
        "rho^2/3 < 1 for alpha < log_5 27 / 2 = 1.0239": sp.N(l5(27) / 2) > 1,
        "rho/5 < rho^-1 for alpha < log_5 45 / 2 = 1.18": sp.N(l5(45) / 2) > 1,
        "alpha* = (log_5 27 + 0.65)/3 in (gamma, 1)": bool(sp.N(ast) > sp.N(l5(3))) and bool(sp.N(ast) < 1),
        "rho^2/3 < rho^-1 for alpha < log_5 81/3 = 0.9101 > alpha*": sp.N(l5(81) / 3 - ast) > 0,
    }
    for m, c in conds.items():
        check(bool(c), m)
    RESULT["alpha_star"] = float(sp.N(ast, 12))
    print(f"   alpha* = (log_5 27 + 0.65)/3 = {RESULT['alpha_star']:.8f}")


# ----------------------------------------------------------------- (B) arb at alpha = 0.8
def certify_08():
    from flint import arb, ctx
    ctx.prec = 256

    def iv(lo, hi):
        lo, hi = arb(lo), arb(hi)
        m = ((lo + hi) / 2).mid()
        return arb(m, max((hi.upper() - m).upper(), (m - lo.lower()).upper()))

    al = arb("0.8")
    rho = arb(5) ** al / 3
    gm = arb(3).log() / arb(5).log()
    kap = 2 * (arb.pi() * al).sin() / arb.pi()
    A, Cal, nu = iv(*A08), iv(*C08), iv(*NU08)
    Dinf, Cp = A / kap, Cal / kap
    check(bool(abs(Cal / A - nu).upper() < arb("1e-20")) and bool(nu.upper() < 0), "AZ (3) as stated: nu_0.8 = C/A < 0, D'_oo = A/kappa, C' = C/kappa")
    G = arb(G_BI4)
    m0 = arb(2) / 3 - 1 / (21 + 6 * arb(10).sqrt())
    K1 = ((arb(1) / 5 + arb("0.14")) * (arb(1) / 5) ** (1 - gm) + arb("1.8") * G * arb(5) ** (gm - 1)) / m0
    check(bool(K1.upper() < arb("1.333045")) and bool(m0.lower() > arb("0.641650")), f"DZ.1's constants: K_1 <= {float(K1.upper()):.6f}, m_0 >= {float(m0.lower()):.6f}")

    def ebar(N, k):
        Lam = 1 / (2 * gm - al) + 1 / (1 + al - 2 * gm)
        Q = (1 / (al + 1) + 1 / (5 * (al + arb(1) / 2))) / m0
        aN = K1 * G * arb(5) ** (-gm) * Lam * (2 * arb(5) ** N) ** (2 * gm - al)
        bN = 2 * arb(5) ** N * Q
        return arb(3) / 2 * aN * arb(3) ** (-k) + bN * (5 * rho) ** (-k) / (1 - 1 / (5 * rho))

    Ep, Em = arb(E60P), arb(E60M)
    k0 = 60
    x0 = rho ** (-k0)
    Dlo = Dinf + x0 * (Cp - Em)              # D'_k >= this for k >= 60 (C' < 0)
    m60 = 1 + nu * x0                        # 1 + nu rho^-k >= this for k >= 60 (nu < 0)
    check(bool(Dlo.lower() > 0) and bool(m60.lower() > 0) and bool((rho / 3).upper() < 1),
          f"k >= 60: D'_k >= {float(Dlo.lower()):.6f} > 0, 1 + nu rho^-k >= {float(m60.lower()):.12f} > 0; rho^k ebar_k decreasing (rho/3 < 1, 1/5 < 1)")
    tab = {}
    for N in NS:
        Cd = iv(*CD_EA[N])
        om = Cd / (-nu)                       # 1 - delta_N
        rb = rho**k0 * ebar(N, k0)
        L = (om * Ep / m60 + rb) / Dlo        # X_k >= C_delta/(1+nu x) - L
        U = om * Em / (m60 * Dlo)             # X_k <= C_delta/(1+nu x) + U
        L, U = arb(L.upper()), arb(U.upper())
        Yoo = -nu * Cd
        Cdlo = arb(Cd.lower())

        def Xb(k):
            x = rho ** (-k)
            base = Cd / (1 + nu * x)
            return base, x

        # windows (dependency-free on C_delta: everything is C_delta times an exact function of x, minus L or plus U)
        Kpos = Kmono = KY = None
        for k in range(k0, KMAX + 1):
            x = rho ** (-k)
            x1 = rho ** (-k - 1)
            # X_k - C_delta >= C_delta (-nu x)/(1 + nu x) - L
            pos = bool((Cdlo * (-nu) * x / (1 + nu * x) - L).lower() > 0)
            # X_k - X_(k+1) >= C_delta (-nu)(x - x1)/((1+nu x)(1+nu x1)) - L - U
            mono = bool((Cdlo * (-nu) * (x - x1) / ((1 + nu * x) * (1 + nu * x1)) - L - U).lower() > 0)
            # Y_k / Y_oo = 1/(1 + nu x) + rho^k (X_k - C_delta/(1+nu x)) / Y_oo ;  Y_k > Y_oo and Y_k > Y_(k+1)
            Yl = 1 / (1 + nu * x) - rho**k * L / (-nu * Cdlo)
            Y1u = 1 / (1 + nu * x1) + rho ** (k + 1) * U / (-nu * Cdlo)
            ydec = bool(Yl.lower() > 1) and bool((Yl - Y1u).lower() > 0)
            if not pos and Kpos is None:
                Kpos = k - 1
            if not mono and Kmono is None:
                Kmono = k - 1
            if not ydec and KY is None:
                KY = k - 1
        # third-order coefficient check at k = 60: rho^k (Y_k/Y_oo - 1)/(-nu) in [lo, hi] ~ 1
        x = x0
        Yl = 1 / (1 + nu * x) - rho**k0 * L / (-nu * Cdlo)
        Yu = 1 / (1 + nu * x) + rho**k0 * U / (-nu * Cdlo)
        t3 = (arb((rho**k0 * (Yl - 1) / (-nu)).lower()), arb((rho**k0 * (Yu - 1) / (-nu)).upper()))
        X60 = (Cd / m60 - L, Cd / m60 + U)
        tab[N] = {"Cd": Cd, "Yoo": Yoo, "L": L, "U": U, "Kpos": Kpos, "Kmono": Kmono, "KY": KY, "t3": t3, "X60": X60}
        print(f"   N = {N:2d}: Y_oo = -nu C_delta in [{float(Yoo.lower()):.6f}, {float(Yoo.upper()):.6f}]; |X_k - C_delta/(1+nu rho^-k)|: "
              f"-{float(L):.3e} .. +{float(U):.3e} (every k >= 60); X_k > C_delta for k in [60, {Kpos}], X_k decreasing on [60, {Kmono}], "
              f"Y_k > Y_oo and decreasing on [60, {KY}]; rho^60 (Y_60/Y_oo - 1)/(-nu) in [{float(t3[0]):.9f}, {float(t3[1]):.9f}]")
    check(all(bool(tab[N]["L"] < arb("6.7e-21")) and bool(tab[N]["U"] < arb("8.7e-35")) for N in NS),
          "EB.3: for every k >= 60 and N in {1,2,3,6,10}: -6.7e-21 <= X_k - C_delta/(1+nu rho^-k) <= 8.7e-35")
    check(all(tab[N]["Kpos"] is not None and tab[N]["Kpos"] >= 200 and tab[N]["Kmono"] >= 200 for N in NS),
          "EB.3: X_k > C_delta and X_k > X_(k+1) for every k in [60, 200] (all five N)")
    check(all(tab[N]["KY"] is not None and tab[N]["KY"] >= 100 for N in NS),
          "EB.3: Y_k = rho^k (X_k - C_delta) > Y_oo and Y_k > Y_(k+1) for every k in [60, 100] (third-order sign, alpha = 0.8 < alpha*)")
    check(all(bool(tab[N]["t3"][0] > arb("1.000008")) and bool(tab[N]["t3"][1] < arb("1.000010")) for N in NS),
          "EB.3: rho^60 (Y_60/Y_oo - 1)/(-nu) in [1.000008, 1.000010] (= 1/(1 + nu rho^-60): third-order coefficient -nu Y_oo)")
    # consistency with DZ.4's level-60 C_delta enclosures and DZ.2's |X_k - C_delta| <= 5.24e-6
    okdz = True
    for N in NS:
        lo, hi = arb(DZ4[N][0]), arb(DZ4[N][1])
        X60lo, X60hi = tab[N]["X60"]
        okdz &= bool(X60lo.upper() <= hi.upper() + arb("1e-5")) and bool(lo.lower() - arb("1e-5") <= X60hi.lower())
        okdz &= bool((tab[N]["Cd"] * (-nu) * x0 / m60 + tab[N]["U"]).upper() < arb("5.24e-6"))
    check(okdz, "consistency: 0 < X_60 - C_delta <= C_delta (-nu) rho^-60/(1+nu rho^-60) + U < 5.24e-6 (DZ.2's rate) and X_60 meets DZ.4's brackets")
    RESULT["eb3"] = tab
    RESULT["nu"] = nu
    return tab


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem EB: delta_N^(n) = (delta_N + nu rho^-n)/(1 + nu rho^-n) + O(5^-0.65n); rho^n (X_n - C_delta) -> nu^2 (1 - delta_N) > 0. Not thrust.")
    symbolic_checks()
    certify_08()
    # floats (OBSERVED)
    fl = {}
    okf = True
    nu08 = float(RESULT["nu"].mid())
    for al, Ns, K, ks in ((0.8, (1, 2, 10), 150, (12, 16, 20, 24, 32)), (0.9, (1, 6), 100, (12, 16, 20, 24, 30)),
                          (0.7, (1, 2), 300, (12, 16, 20, 24)), (0.95, (1, 2), 90, (12, 16, 20, 24, 30))):
        rows = nu_from_levels(al, Ns, ks, K)
        fl[al] = rows
        last = [rows[N][-1][1] for N in Ns]
        spread = max(last) - min(last)
        print(f"   floats alpha = {al}: nu_k(N) = " + "; ".join(f"N={N}: " + ", ".join(f"{v:.9f}" for _, v in rows[N]) for N in Ns))
        okf &= spread < 1e-8 and all(v < 0 for v in last)
        # error ratio per level ~ rho/3 (R_n ~ 3^-n observed)
        N = Ns[0]
        errs = [abs(v - last[0]) for _, v in rows[N][:3]]
        okf &= errs[1] < errs[0] and errs[2] < errs[1]
        if al == 0.8:
            okf &= abs(last[0] - nu08) < 2e-6
    check(okf, "floats: nu_k(N) is N-independent at the last k (24 at 0.7, 30 at 0.9/0.95, 32 at 0.8; spread < 1e-8) at alpha = 0.7, 0.8, 0.9, 0.95, and at 0.8 within 2e-6 of AZ's nu (DZ.5's quadrature bias ~5e-7)")
    RESULT["floats"] = fl
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EB checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
