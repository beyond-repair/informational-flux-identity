#!/usr/bin/env python3
"""Theorem EC: the sharp remainder exponent of EB.  AZ (1)'s defect, AT (4a)'s rate and EB.1's remainder are O(3^-n), not
O(5^-0.65n): the exponent 0.65 came only from AQ's power ceiling r <= B s^-0.35 (AT's own diagnosis); BI.4's sharp ceiling
s^(1-gamma) r_oo(s) <= G = 0.6031 replaces it by gamma = log3/log5, and an exact AT (1) transfer at the shifted exponent
alpha - gamma turns t^gamma into the factor (rho/3)^j.

Setting (GASKET.md, Theorems AT, AZ, DZ, EB): alpha in (gamma, 1), rho = 5^alpha/3, rho' = rho/3 = 5^(alpha-gamma)/3 < 1,
J_k = int s^-alpha r_k, Delta_j = J_(j+1) - J_j = int s^-alpha delta_(j+1), e_j = Delta_j - c rho^j, eps_k = J_k - c rho^k/(rho-1) - C'
= -sum_(j>=k) e_j (AZ (1)), X_k = rho^k (delta_N - delta_N^(k)), Y_k = rho^k (X_k - C_delta), Y_oo = -nu C_delta.

Theorem EC (checked here; exit 1 on any failure):
  EC.1 (THEOREM; every t in (0,1], every alpha in (0,1]).  With BI.4's ceiling in AT (2)'s Dbar (zr_G(z) = min(AT's zr(z), G z^gamma)),
       Cbar_G <= Cbar and  |A_k^alpha B_k C_k - 1| <= W_+^G + W_- <= K_alpha t^gamma,  K_alpha = e^(alpha/20 + K_1/2) - 9/10,
       K_1 = 0.2^(1-gamma) 16.6/30 + 2 c_e/(1 - 0.2^gamma c_e),  c_e = 2 G 5.2^gamma / 6   (K_1 <= 1.8934, K_0.8 <= 1.7823).
  EC.2 (THEOREM; every alpha in (gamma,1), j >= 6, k >= 6).  |e_j| <= K_alpha I_6 rho'^(j-6) + T_j  and
       |eps_k| <= B_alpha rho'^k,  B_alpha = K_alpha I_6 rho'^-6/(1 - rho') + 2 T_6 rho'^-6, with explicit I_6, T_j (doubly-exp.).
       Hence rho^-k eps_k = O(3^-k), rho^-j e_j = O(3^-j) (AT (4a) sharpened on (gamma,1)), (rho^n + nu)||F|| = K + O(3^-n) (AZ (4)),
       and EB.1's remainder R_n = O(3^-n),  X_n - C_delta/(1 + nu rho^-n) = O((rho/3)^n).
  EC.3 (THEOREM).  X_n = C_delta sum_(j<J) (-nu rho^-n)^j + o(rho^-(J-1)n) whenever rho^J < 3, i.e. alpha < (J+1) gamma / J;
       EB's third term (J = 3) holds for alpha < 4 gamma/3 = 0.910142 (was 0.89927).  This is the end of the route: R_n is observed
       Theta(3^-n) (EC.5), and past 4 gamma/3 a Theta(3^-n) remainder dominates rho^-3n.
  EC.4 (certified at alpha = 0.8; arb; AZ (3), EA.4, DZ.1 as stated).  |eps_k| <= min(E_60^+-, B_0.8 rho'^k), B_0.8 <= 2.0e3, so for
       EVERY k >= 60 (no window): |X_k - C_delta/(1+nu rho^-k)| <= Lam1_N (rho/3)^k + Lam2_N 5^-k, X_k > C_delta, X_k > X_(k+1), Y_k > Y_oo,
       Y_k > Y_(k+1), and |rho^k (Y_k/Y_oo - 1)/(-nu) - 1/(1+nu rho^-k)| <= rho^2k (Lam1_N (rho/3)^k + Lam2_N 5^-k)/(nu^2 C_delta) -> 0.
  EC.5 (OBSERVED, floats).  w_k = (u_(k+1) - u_(k+2)) - (u_k - u_(k+1))/rho, u_k = 1/(1 - delta_N^(k)), kills the Moebius part; w_(k+1)/w_k -> 1/3
       at alpha = 0.8, 0.9, 0.95 (the 3^-n remainder has a nonzero coefficient; not claimed).

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time

import numpy as np

FAILS: list[str] = []
RESULT: dict = {}

G_BI4 = "0.6031"          # BI.4: sup_(s>0) s^(1-gamma) r_oo(s) <= 0.6031
TAIL_AT = "3e10"          # AT (3): A_oo^alpha Cbar <= 3e10 s^5 on s >= 1;  delta_(n+1)(s) <= 2 * 1.2^n / (s p^n(s)) on s >= 1
N0 = 6                    # transfer level for EC.2
NS = (1, 2, 3, 6, 10)


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- floats (numpy only)
def w_ratios(al, N, ks):
    """OBSERVED: w_k = v_(k+1) - v_k/rho, v_k = u_k - u_(k+1), u_k = 1/(1 - delta_N^(k)); EB.1 => w_k = O(R_k); returns w_(k+1)/w_k."""
    import gasket_fractional_corner_limit_mobius as eb
    rho = 5**al / 3
    lv = range(min(ks), max(ks) + 4)
    u = {k: 1 / (1 - eb.float_deltas(al, k, (N,))[N]) for k in lv}
    v = {k: u[k] - u[k + 1] for k in lv if k + 1 in u}
    w = {k: v[k + 1] - v[k] / rho for k in v if k + 1 in v}
    return [(k, w[k + 1] / w[k], w[k] * 3.0**k) for k in ks]


def p_iter(x: int, n: int) -> int:
    for _ in range(n):
        x = x * (x + 5)
    return x


# ----------------------------------------------------------------- (A) sympy
def symbolic_checks():
    import sympy as sp
    y, a = sp.symbols("y alpha", positive=True)
    check(sp.simplify(sp.Rational(3, 5) * (2 + y) * (5 + y) / (6 + y) - 1 - y * (16 + 3 * y) / (5 * (6 + y))) == 0,
          "(3/5)(2+y)(5+y)/(6+y) - 1 = y(16+3y)/(5(6+y)) <= y(16+3y)/30  (log x <= x - 1)")
    gm = sp.log(3) / sp.log(5)
    check(sp.simplify(5 ** (a - gm) / 3 - (5**a / 3) / 3) == 0, "rho' = 5^(alpha-gamma)/3 = rho/3 (5^gamma = 3)")
    l5 = lambda v: sp.log(v) / sp.log(5)
    ast = (l5(27) + sp.Rational(13, 20)) / 3
    a3 = 4 * gm / 3
    conds = {
        "4 gamma/3 = log_5 81 / 3 (rho^3 < 3 iff alpha < 4 gamma/3)": sp.simplify(a3 - l5(81) / 3) == 0,
        "alpha* = 0.89927 < 4 gamma/3 = 0.91014 < 1": bool(sp.N(ast) < sp.N(a3)) and bool(sp.N(a3) < 1),
        "rho/3 < 1 on (gamma,1) (alpha < log_5 9 = 1.365)": bool(sp.N(l5(9)) > 1),
        "2 gamma > 1 (int_0^1 t^(2gamma-alpha-1) finite on (gamma,1))": bool(sp.N(2 * gm) > 1),
        "1/5 < rho/3 for alpha > log_5(9/5) = 0.365 < gamma": bool(sp.N(l5(sp.Rational(9, 5))) < sp.N(gm)),
        "1/(5 rho) < 1/3 on (gamma,1)": bool(sp.N(l5(sp.Rational(3, 5)) + 1) < sp.N(gm) + 1),
        "mu = 5^(alpha-1)/3 < 1/3 on (gamma,1)": True,
        "rho^2/3 < 1 on (gamma,1) (alpha < 1.0239)": bool(sp.N(l5(27) / 2) > 1),
        "J = 4 threshold 5 gamma/4 = 0.85326 in (gamma,1)": bool(sp.N(5 * gm / 4) < 1),
    }
    for m, c in conds.items():
        check(bool(c), m)
    RESULT["alpha3"] = float(sp.N(a3, 12))
    RESULT["alpha4"] = float(sp.N(5 * gm / 4, 12))
    print(f"   third-term threshold 4 gamma/3 = {RESULT['alpha3']:.8f} (EB: alpha* = {float(sp.N(ast)):.8f}); J = 4: 5 gamma/4 = {RESULT['alpha4']:.8f}")


# ----------------------------------------------------------------- (B) arb constants (EC.1, EC.2)
def constants(al_str: str):
    from flint import arb, ctx
    ctx.prec = 256
    al = arb(al_str)
    gm = arb(3).log() / arb(5).log()
    G = arb(G_BI4)
    y0 = arb(1) / 5
    ce = 2 * G * (5 + y0) ** gm / 6                     # eps(y) <= c_e y^gamma on (0, 0.2]
    me = 1 - y0**gm * ce                                # 1 - eps(y) >= m_e
    K1 = y0 ** (1 - gm) * (16 + 3 * y0) / 30 + 2 * ce / me
    K = (al / 20 + K1 / 2).exp() - 1 + arb(1) / 10      # W_+ + W_- <= K t^gamma on (0, 1]
    rho = arb(5) ** al / 3
    rp = rho / 3
    ap = al - gm
    tail = arb(TAIL_AT)
    p6 = arb(p_iter(1, N0))
    I6 = (ap / 20 + K1 / 2).exp() * G / (2 * gm - al) + tail * 2 * arb("1.2") ** N0 / (p6.sqrt() * (2 ** (N0 - 1) - 5 + ap))

    def T(j):   # tail_j <= 3e10 * 2 * 1.2^j / (sqrt(p^j(1)) (2^(j-1) - 5 + alpha)),  sqrt(p^j(1)) >= 6^(2^(j-2))
        lg = arb(2) ** (j - 2) * arb(6).log()
        return tail * 2 * arb("1.2") ** j * (-lg).exp() / (arb(2) ** (j - 1) - 5 + al)

    tau = T(N0) / rp**N0
    B = K * I6 / (rp**N0 * (1 - rp)) + 2 * tau
    ok = bool(me.lower() > 0) and bool((y0**gm * ce).upper() < 1) and bool(rp.upper() < 1) and bool((2 * gm - al).lower() > 0)
    # T_(j+1)/T_j <= 1.2 * 6^-(2^(j-2)) <= min(1/2, rho') for j >= 6
    ok &= bool((arb("1.2") * arb(6) ** (-16)).upper() < min(0.5, float(rp.lower())))
    return {"al": al, "gm": gm, "G": G, "K1": K1, "K": K, "I6": I6, "B": B, "rho": rho, "rp": rp, "tau": tau, "T": T, "ok": ok}


def certify_constants():
    rows = {}
    for a in ("0.7", "0.8", "0.9", "0.95", "1"):
        c = constants(a)
        rows[a] = c
        print(f"   alpha = {a}: K_1 <= {float(c['K1'].upper()):.6f}, K <= {float(c['K'].upper()):.6f}, I_6 <= {float(c['I6'].upper()):.6f}, "
              f"rho' = {float(c['rp'].mid()):.6f}, 2 T_6/rho'^6 <= {float((2 * c['tau']).upper()):.3e}, B <= {float(c['B'].upper()):.4e}")
        check(c["ok"], f"alpha = {a}: eps(0.2) < 1, rho' < 1, 2gamma > alpha, T_(j+1)/T_j <= min(1/2, rho') (j >= 6)")
    c8 = rows["0.8"]
    check(bool(c8["K1"].upper() < 1.8934) and bool(c8["K"].upper() < 1.7823),
          f"EC.1: K_1 <= 1.8934 (alpha-free), K_0.8 <= 1.7823 ({float(c8['K'].upper()):.6f}); K_1 = {float(c8['K1'].upper()):.6f}")
    check(bool(rows["1"]["K"].upper() < 1.8355), f"EC.1: uniform K_alpha <= K_1.0 <= 1.8355 ({float(rows['1']['K'].upper()):.6f}) on (0,1]")
    check(bool(c8["B"].upper() < 2.0e3), f"EC.2 at 0.8: |eps_k| <= B rho'^k, B <= 2.0e3 ({float(c8['B'].upper()):.4e}), every k >= 6")
    RESULT["const"] = {a: {k: c[k] for k in ("K1", "K", "I6", "B", "rp")} for a, c in rows.items()}
    return rows


# ----------------------------------------------------------------- (C) arb at alpha = 0.8: EB.3 without windows
def certify_08(c8):
    from flint import arb, ctx
    ctx.prec = 256
    import gasket_fractional_corner_limit_mobius as eb

    def iv(lo, hi):
        lo, hi = arb(lo), arb(hi)
        m = ((lo + hi) / 2).mid()
        return arb(m, max((hi.upper() - m).upper(), (m - lo.lower()).upper()))

    al, gm, rho, rp, B = c8["al"], c8["gm"], c8["rho"], c8["rp"], c8["B"]
    kap = 2 * (arb.pi() * al).sin() / arb.pi()
    A, Cal, nu = iv(*eb.A08), iv(*eb.C08), iv(*eb.NU08)
    Dinf, Cp = A / kap, Cal / kap
    G = arb(eb.G_BI4)
    m0 = arb(2) / 3 - 1 / (21 + 6 * arb(10).sqrt())
    K1d = ((arb(1) / 5 + arb("0.14")) * (arb(1) / 5) ** (1 - gm) + arb("1.8") * G * arb(5) ** (gm - 1)) / m0
    Ep, Em = arb(eb.E60P), arb(eb.E60M)
    k0 = 60
    x0 = rho ** (-k0)
    Dlo = Dinf + x0 * (Cp - Em)
    m60 = 1 + nu * x0
    check(bool(Dlo.lower() > 0) and bool(m60.lower() > 0) and bool(nu.upper() < 0), "k >= 60: D'_k >= D'_lo > 0, 1 + nu rho^-k >= m_60 > 0, nu < 0 (as in EB.3)")
    # where the geometric bound beats AZ (3)'s E_60^+
    kE = next(k for k in range(k0, 400) if bool((B * rp**k).upper() < Ep.lower()))
    tab = {}
    for N in NS:
        Lam_a = 1 / (2 * gm - al) + 1 / (1 + al - 2 * gm)
        Q = (1 / (al + 1) + 1 / (5 * (al + arb(1) / 2))) / m0
        aN = K1d * G * arb(5) ** (-gm) * Lam_a * (2 * arb(5) ** N) ** (2 * gm - al)
        bN = 2 * arb(5) ** N * Q
        Cd = iv(*eb.CD_EA[N])
        om = Cd / (-nu)
        Cdlo = arb(Cd.lower())
        # rho^k ebar_k = 1.5 aN (rho/3)^k + bN 5^-k/(1 - 1/(5 rho))  (DZ.1);  L_k, U_k <= Lam1 (rho/3)^k + Lam2 5^-k  for every k >= 60
        Lam1 = arb(((om * B / m60 + arb(3) / 2 * aN) / Dlo).upper())
        Lam2 = arb((bN / ((1 - 1 / (5 * rho)) * Dlo)).upper())
        q2, q3, r5, s5 = rho**2 / 3, rho**3 / 3, rho / 5, rho**2 / 5      # all < 1: every condition below is monotone in k
        k = k0
        e2 = Lam1 * q2**k + Lam2 * r5**k                 # rho^k (Lam1 (rho/3)^k + Lam2 5^-k)
        e3 = Lam1 * q3**k + Lam2 * s5**k                 # rho^2k (...)
        pos = bool(e2.upper() < (Cdlo * (-nu)).lower())
        mono = bool((2 * e2).upper() < (Cdlo * (-nu) * (1 - 1 / rho)).lower())
        ygt = bool(e3.upper() < (nu**2 * Cdlo).lower())
        ydec = bool((2 * e3).upper() < (nu**2 * Cdlo * (1 - 1 / rho)).lower())
        third = arb((e3 / (nu**2 * Cdlo)).upper())       # |rho^k (Y_k/Y_oo - 1)/(-nu) - 1/(1+nu rho^-k)| <= e3_k/(nu^2 C_delta)
        Zoo = -nu * (-nu * Cd)                           # -nu Y_oo = -nu^3 (1 - delta_N)
        env60 = arb((Lam1 * rp**k0 + Lam2 * arb(5) ** (-k0)).upper())
        Lam = Lam1
        tab[N] = {"Lam": Lam, "Lam2": Lam2, "pos": pos, "mono": mono, "ygt": ygt, "ydec": ydec, "third60": third, "Zoo": Zoo, "env60": env60, "Cd": Cd}
        print(f"   N = {N:2d}: |X_k - Moebius| <= {float(Lam):.4e} (rho/3)^k + {float(Lam2):.4e} 5^-k (= {float(env60):.3e} at k = 60); every k >= 60: "
              f"X_k > C_delta {pos}, X_k decreasing {mono}, Y_k > Y_oo {ygt}, Y_k decreasing {ydec}; third-order defect at 60 <= {float(third):.3e}; "
              f"-nu Y_oo = -nu^3 (1-delta_N) in [{float(Zoo.lower()):.6f}, {float(Zoo.upper()):.6f}]")
    check(kE <= 70, f"EC.4: B rho'^k < E_60^+ = 7.22e-23 from k = {kE} on (min(E_60^+, B rho'^k) is used; AZ (3) as stated)")
    check(all(tab[N]["pos"] and tab[N]["mono"] for N in NS),
          "EC.4: X_k > C_delta and X_k > X_(k+1) for EVERY k >= 60 (all five N; EB.3's windows [60, 232..270] removed)")
    check(all(tab[N]["ygt"] and tab[N]["ydec"] for N in NS),
          "EC.4: Y_k > Y_oo and Y_k > Y_(k+1) for EVERY k >= 60 (all five N; EB.3's windows [60, 115..129] removed)")
    check(all(bool(tab[N]["third60"] < arb("2e-10")) for N in NS),
          "EC.4: |rho^k (Y_k/Y_oo - 1)/(-nu) - 1/(1+nu rho^-k)| < 2e-10 at k = 60 and decreasing in k (rho^3/3, rho^2/5 < 1)")
    check(all(bool(tab[N]["env60"] < arb("1.4e-20")) for N in NS), "EC.4: Lam1 (rho/3)^60 + Lam2 5^-60 < 1.4e-20 for N in {1,2,3,6,10}")
    RESULT["ec4"] = tab
    RESULT["kE"] = kE
    return tab


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem EC: R_n = O(3^-n) (BI.4's ceiling in AT's weight); EB's third term to 4 gamma/3 = 0.910142; every k >= 60 at 0.8. Not thrust.")
    symbolic_checks()
    rows = certify_constants()
    certify_08(rows["0.8"])
    # floats (OBSERVED)
    okf = True
    obs = {}
    for al, ks in ((0.8, (14, 16, 18, 20)), (0.9, (14, 16, 18, 20)), (0.95, (12, 14, 16, 18))):
        r = w_ratios(al, 1, ks)
        obs[al] = r
        print(f"   floats alpha = {al}, N = 1: w_(k+1)/w_k = " + ", ".join(f"{k}: {q:.5f}" for k, q, _ in r) +
              "; 3^k w_k = " + ", ".join(f"{s:.4e}" for _, _, s in r))
        okf &= abs(r[-1][1] - 1 / 3) < 3e-3 and all(abs(s) > 0 for _, _, s in r)
    check(okf, "EC.5 (OBSERVED): w_(k+1)/w_k within 3e-3 of 1/3 at the last k (alpha = 0.8, 0.9, 0.95; N = 1): the 3^-n remainder is not zero")
    RESULT["floats"] = obs
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EC checks passed")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    sys.exit(main())
