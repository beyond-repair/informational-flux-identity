#!/usr/bin/env python3
"""Theorem DZ: the renormalised corner deficit has a limit, and the limit depends on N.  Conjecture BL.6's limit part
proved with the limit identified; its clause "independent of N" is false (kept failure, certified at alpha = 0.8);
and C_alpha < 0 on all of (gamma, 1) (AZ's open item).

Setting (GASKET.md, Theorems AZ, BK, BL, BN): alpha in (gamma, 1), gamma = log3/log5, rho = 5^alpha/3 > 1,
p(s) = s(s+5), q = p^-1, beta(s) = (s+6)/((s+2)(s+5)), e(s) = 2s/(2+s), alpha_1(z) = (z+4)/((z+2)(z+5)), r_k the corner
resolvent (AQ (1): r_(k+1)(s) = r_k(ps)/(beta(s) - e(s) r_k(ps)), r_0 = 1/(3+s), r_k increasing in k, decreasing in s),
J_k = int_0^oo s^-alpha r_k ds (= BK's D_k, D'_k = rho^-k J_k), eta_N = alpha_1 + beta eta_(N-1) o p, eta_0 = -1/2.
BL's corner deficit: delta_N^(k) = J_k^-1 int s^-alpha r_k(s) (1 - eta_N(p^(k-N) s)) ds (k >= N), delta_N = lim_k.

Theorem DZ (checked here; exit 1 on any failure):
  DZ.1 (exact transfer and a 3^-k gap; every alpha in (gamma,1), N >= 1, k >= N).  H_k := J_k (1 - delta_N^(k)) =
       int s^-alpha r_k(s) eta_N(p^(k-N) s) ds obeys  H_(k+1) = rho int t^-alpha r_k(t) eta_N(p^(k-N) t) omega_k(q t) dt
       (t = p(s); omega_k = BK.1's weight >= 1, BK.2).  So rho^-k H_k increases to h_oo = (1 - delta_N) D'_oo and
         0 <= h_oo - rho^-k H_k <= ebar_k = (3/2) a_N 3^-k + b_N (5 rho)^-k / (1 - 1/(5 rho)),
       a_N = K_1 G 5^-gamma Lambda_alpha (2 5^N)^(2 gamma - alpha), b_N = 2 5^N Q_alpha, from omega_k - 1 <= K_1 s^gamma
       (s <= q(1)), omega_k <= (1 + sqrt(t)/5)/m_0 (t >= 1), eta_N(z) <= min(1, 2/z), r_k(t) <= min(G t^(gamma-1), 1/t).
  DZ.2 (THEOREM, from DZ.1 and AZ (1) as stated).  X_k := rho^k (delta_N - delta_N^(k)) -> C_delta(alpha, N)
       = -nu_alpha (1 - delta_N),  nu_alpha = C_alpha / A_alpha  (AZ), with the exact identity
         X_k = [-h_oo (C' + eps_k) - rho^k eta_k D'_oo] / (D'_k D'_oo),  eps_k in [-E_k^-, E_k^+],  0 <= eta_k <= ebar_k.
  DZ.3 (COROLLARY, with BN.4 as stated).  c_N <= -nu_alpha (1 - delta_N) for every N, so nu_alpha < 0 and C_alpha < 0
       for every alpha in (gamma, 1) (AZ: proved only near gamma, certified at 0.75, 0.8).
  DZ.4 (KEPT FAILURE of BL.6's "independent of N", certified at alpha = 0.8).  arb monotone brackets of delta_N^(60)
       (N = 1, 2, 3, 6, 10; 69052 cells) and DZ.2's own ceiling X_60 <= (1 - delta_N^(60)) (|C'| + E_60^-) / D'_60 give
       disjoint enclosures of C_delta(0.8, N): [0.030733, 0.031869], [0.158436, 0.159369], [0.257779, 0.258555],
       [0.461214, 0.461662], [0.602405, 0.602621]; C_delta(N) -> -nu_alpha as N -> oo (delta_N = O(rho^-N), BJ (6)).
       BL's float 0.6025 (N = 10) lies in C_delta(0.8, 10).  Explicit rate from AZ (3)'s E_60: |X_k - C_delta(0.8, N)|
       <= 5.24e-6 for every k >= 60 and N <= 10.
  Floats (consistency only): X_k for alpha = 0.8, 0.9 against -nu (1 - delta_N); rho^-k H_k increasing; gap <= ebar_k.

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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasket_fractional_uniform_limit as bi  # noqa: E402

FAILS: list[str] = []
RESULT: dict = {}

# Theorem AZ (3), used as stated (alpha = 0.8)
A08 = ("2.3964097039792205270396", "2.3964097039792205270397")
C08 = ("-1.7432571783920757239548", "-1.7432571783920757239546")
NU08 = ("-0.7274453844421553283991", "-0.7274453844421553283990")
E60P, E60M = "7.22e-23", "6.7e-34"
# BK.4 / BN.4 ceiling at alpha = 0.8 as stated: delta_N - delta_N^(n) <= 4.6644 rho^-n
BK08 = "4.6644"
# BN.5's certified floors c_N at alpha = 0.8 (stated)
BN_C08 = {1: "1.758e-3", 2: "8.641e-3", 4: "1.566e-2", 8: "1.914e-2"}
BL_FLOAT_08 = 0.6025     # BL.6's float (N = 10)
BM_FLOAT_09 = 0.3643     # BM's float rho^n (delta^(60) - delta^(n)) at alpha = 0.9, N = 6, n = 20
NS_CERT = (1, 2, 3, 6, 10)
N_LEVEL = 60


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def iv(lo: str, hi: str) -> arb:
    return hull(arb(lo), arb(hi))


def hull(lo: arb, hi: arb) -> arb:
    m = ((lo + hi) / 2).mid()
    r = max((hi.upper() - m).upper(), (m - lo.lower()).upper())
    return arb(m, r)


def gam_arb() -> arb:
    return arb(3).log() / arb(5).log()


# ----------------------------------------------------------------- (A) sympy identities
def symbolic_checks():
    s, z = sp.symbols("s z", positive=True)
    p = s * (s + 5)
    be = (s + 6) / ((s + 2) * (s + 5))
    e = 2 * s / (2 + s)
    B = 1 - s * (s + 4) / (3 * (s + 2) * (s + 5))
    C = sp.Rational(5, 3) * (1 + 2 * s / 5) * e
    check(sp.simplify((1 + 2 * s / 5) * sp.Rational(5, 3) * be - B) == 0, "(1+2s/5)(5/3) beta = B = 1 - s(s+4)/(3(s+2)(s+5)) (BK.2's form)")
    check(sp.simplify(B - C / p - (s + 4) * (2 * s + 5) / (3 * (s + 2) * (s + 5))) == 0, "B - C/p = (s+4)(2s+5)/(3(s+2)(s+5)) > 0 (BN.2's identity)")
    check(sp.simplify(sp.Rational(2, 3) - (s + 4) * (2 * s + 5) / (3 * (s + 2) * (s + 5)) - s / (3 * s**2 + 21 * s + 30)) == 0,
          "2/3 - (B - C/p) = s/(3s^2+21s+30) <= 1/(21 + 6 sqrt 10) (AM-GM: 3s + 30/s >= 6 sqrt 10), so B - C/p >= m_0")
    a1 = (z + 4) / ((z + 2) * (z + 5))
    bz = (z + 6) / ((z + 2) * (z + 5))
    check(sp.expand((z + 2) * (z + 5) - z * (z + 4)) == 3 * z + 10 and sp.expand((z + 2) * (z + 5) - z * (z + 6)) == z + 10,
          "alpha_1(z) <= 1/z and beta(z) <= 1/z, so eta_N(z) <= 2/z (N >= 2; eta_1 = 1/(2(z+5)) <= 1/(2z)); with BL.1, eta_N <= min(1, 2/z)")
    check(sp.simplify(a1 + bz * (-sp.Rational(1, 2)) - 1 / (2 * (z + 5))) == 0, "eta_1 = 1/(2(z+5)) (BL.1)")
    check(sp.simplify((1 - B) / s - (s + 4) / (3 * (s + 2) * (s + 5))) == 0 and sp.simplify(C / s - 2 * (5 + 2 * s) / (3 * (2 + s))) == 0,
          "(1-B)/s = (s+4)/(3(s+2)(s+5)) <= 0.14 and C/s = 2(5+2s)/(3(2+s)) <= 1.8 on (0, 1/5]")
    # monotonicity on (0, 1/5] of the two quotients: (s+4)/((s+2)(s+5)) decreasing, (5+2s)/(2+s) = 2 + 1/(2+s) decreasing
    d1 = sp.diff((s + 4) / ((s + 2) * (s + 5)), s)
    check(sp.factor(sp.numer(sp.together(d1))) == sp.factor(-(s**2 + 8 * s + 18)), "(s+4)/((s+2)(s+5)) decreasing (so its sup on (0,1/5] is at 0: 4/30 <= 0.14)")
    check(sp.simplify((5 + 2 * s) / (2 + s) - 2 - 1 / (2 + s)) == 0, "(5+2s)/(2+s) = 2 + 1/(2+s) decreasing (sup 5/2 at 0: C/s <= 5/3 <= 1.8)")
    check(sp.simplify(q1_exact() - (sp.sqrt(29) - 5) / 2) == 0 and sp.N((sp.sqrt(29) - 5) / 2) < sp.Rational(1, 5),
          "q(1) = (sqrt 29 - 5)/2 < 1/5, and q(t) <= t/5, p(s) >= 5s, p^j(t) >= 5^j t")


def q1_exact():
    x = sp.symbols("x", positive=True)
    return [r for r in sp.solve(x * (x + 5) - 1, x) if r.is_positive][0]


# ----------------------------------------------------------------- (B) exact transfer, numeric spot checks (arb)
def r_point(s: arb, n: int, cap: arb):
    """Enclosure of r_n(s) and the orbit z_j = p^j s (stopped once z_j > cap: r_(n-j)(z_j) in [1/(3+z_j), 1/z_j])."""
    zs = [s]
    while len(zs) <= n and zs[-1].upper() <= cap.upper():
        if len(zs) == n + 1:
            break
        z = zs[-1]
        zs.append(z * (z + 5))
    J = len(zs) - 1
    zJ = zs[J]
    if J == n:
        r = 1 / (3 + zJ)
    else:
        r = hull(1 / (3 + zJ), 1 / zJ)
    for j in range(J - 1, -1, -1):
        y = zs[j]
        be = (y + 6) / ((y + 2) * (y + 5))
        ee = 2 * y / (2 + y)
        r = r / (be - ee * r)
    return r, zs


def eta_on_orbit(zs, n: int, N: int) -> arb:
    """Enclosure of eta_N(z_(n-N)) from the orbit; tails eta_M(z) in [0, 2/z] (M >= 1)."""
    J = len(zs) - 1
    if J >= n:
        e = arb(-1) / 2
        top = n
    elif J <= n - N:
        return hull(arb(0), 2 / zs[J])
    else:
        e = hull(arb(0), 2 / zs[J])
        top = J
    for k in range(top - 1, n - N - 1, -1):
        zz = zs[k]
        D = (zz + 2) * (zz + 5)
        e = (zz + 4) / D + (zz + 6) / D * e
    return e


def transfer_spot_checks():
    """omega_k(s) * rho t^-alpha / (s^-alpha / ((2s+5)(beta - e r_k(t)))) = 1 at sample s, k, alpha (arb, r_k(t) exact recursion)."""
    ctx.prec = 120
    ok = True
    for al in ("0.7", "0.8", "0.95"):
        a = arb(al)
        rho = arb(5) ** a / 3
        for k in (0, 1, 3, 6):
            for sv in ("0.001", "0.3", "2", "17"):
                s = arb(sv)
                t = s * (s + 5)
                rk, _ = r_point(t, k, arb(2) ** 400)
                rk1, _ = r_point(s, k + 1, arb(2) ** 400)
                be = (s + 6) / ((s + 2) * (s + 5))
                ee = 2 * s / (2 + s)
                om = (1 + s / 5) ** a / ((1 + 2 * s / 5) * arb(5) / 3 * (be - ee * rk))
                lhs = s ** (-a) * rk1 / (2 * s + 5)              # integrand of H_(k+1) in t (ds = dt/(2s+5)), eta factor common
                rhs = rho * t ** (-a) * rk * om
                ok &= bool(abs(lhs / rhs - 1).upper() < arb("1e-25")) and bool(om.lower() >= 1)
    check(ok, "exact transfer: s^-a r_(k+1)(s) ds = rho t^-a r_k(t) omega_k(s) dt, t = p(s), and omega_k >= 1 (alpha = 0.7, 0.8, 0.95; k = 0,1,3,6; 4 points)")


# ----------------------------------------------------------------- (C) DZ.1 constants (arb)
def dz1_constants(G: arb):
    ctx.prec = 120
    gm = gam_arb()
    m0 = arb(2) / 3 - 1 / (21 + 6 * arb(10).sqrt())
    K1 = ((arb(1) / 5 + arb("0.14")) * (arb(1) / 5) ** (1 - gm) + arb("1.8") * G * arb(5) ** (gm - 1)) / m0
    check(bool(m0.lower() > arb("0.6416")), f"m_0 = 2/3 - 1/(21+6 sqrt10) >= {float(m0.lower()):.6f}")
    check(bool(K1.upper() < arb("1.34")), f"K_1 = ((1/5+0.14)(1/5)^(1-gamma) + 1.8 G 5^(gamma-1))/m_0 <= {float(K1.upper()):.6f}: omega_k - 1 <= K_1 s^gamma on s <= 1/5")
    # direct arb check of omega_hi - 1 <= K_1 s^gamma on a grid of (0, 1/5] with the majorant d >= B - C min(G p^(gamma-1), 1/p)
    ok = True
    for j in range(1, 401):
        s = arb(j) / 2000
        p = s * (s + 5)
        B = 1 - s * (s + 4) / (3 * (s + 2) * (s + 5))
        Cc = arb(5) / 3 * (1 + 2 * s / 5) * 2 * s / (2 + s)
        rb = G * p ** (gm - 1)
        rb = rb if bool(rb.upper() < (1 / p).lower()) else 1 / p
        d = B - Cc * rb
        om = (1 + s / 5) ** arb(1) / d          # alpha <= 1 worst case
        ok &= bool((om - 1).upper() <= (K1 * s ** gm).lower())
    check(ok, "grid check (400 points of (0, 1/5], alpha = 1 worst case): omega_hi(s) - 1 <= K_1 s^gamma")
    return {"m0": m0, "K1": K1, "gm": gm}


def ebar(al: arb, N: int, k: int, cst, G: arb) -> arb:
    gm, m0, K1 = cst["gm"], cst["m0"], cst["K1"]
    rho = arb(5) ** al / 3
    Lam = 1 / (2 * gm - al) + 1 / (1 + al - 2 * gm)
    Q = (1 / (al + 1) + 1 / (5 * (al + arb(1) / 2))) / m0
    aN = K1 * G * arb(5) ** (-gm) * Lam * (2 * arb(5) ** N) ** (2 * gm - al)
    bN = 2 * arb(5) ** N * Q
    return arb(3) / 2 * aN * arb(3) ** (-k) + bN * (5 * rho) ** (-k) / (1 - 1 / (5 * rho))


# ----------------------------------------------------------------- (D) certified brackets of delta_N^(n) at alpha = 0.8
def delta_brackets(al="0.8", n=N_LEVEL, Ns=NS_CERT, lam=1.002, a0exp=20, b1=1e4, prec=80):
    """Monotone brackets: r_n decreasing, g = 1 - eta_N(p^(n-N) s) increasing (BL.1); cells with exact dyadic endpoints."""
    ctx.prec = prec
    a = arb(al)
    a0 = 5.0 ** (-n - a0exp)
    m = int(math.ceil(math.log(b1 / a0) / math.log(lam)))
    xs = [arb(a0 * lam**i).mid() for i in range(m)] + [arb(b1).mid()]
    check(all(bool(xs[i] < xs[i + 1]) for i in range(m)), f"bracket cells: {m} cells with exact dyadic endpoints, strictly increasing, covering [5^-{n + a0exp}, {b1:g}]")
    cap = arb(2) ** 200
    rs, gs = [], {N: [] for N in Ns}
    for x in xs:
        r, zs = r_point(x, n, cap)
        rs.append(r)
        for N in Ns:
            gs[N].append(1 - eta_on_orbit(zs, n, N))
    one_a = 1 - a
    pw = [x ** one_a for x in xs]
    w = [(pw[i + 1] - pw[i]) / one_a for i in range(m)]
    den_lo = sum((arb(rs[i + 1].lower()) * arb(w[i].lower()) for i in range(m)), arb(0))
    den_hi = sum((arb(rs[i].upper()) * arb(w[i].upper()) for i in range(m)), arb(0))
    r0 = (arb(5) / 3) ** n / 3
    tail0 = r0 * xs[0] ** one_a / one_a
    tailinf = xs[-1] ** (-a) / a
    den_hi = den_hi + tail0 + tailinf
    out = {}
    for N in Ns:
        g = gs[N]
        num_lo = sum((arb(rs[i + 1].lower()) * arb(g[i].lower()) * arb(w[i].lower()) for i in range(m)), arb(0))
        num_hi = sum((arb(rs[i].upper()) * arb(g[i + 1].upper()) * arb(w[i].upper()) for i in range(m)), arb(0))
        num_hi = num_hi + tail0 * arb(g[0].upper()) + tailinf
        out[N] = (arb((num_lo / den_hi).lower()), arb((num_hi / den_lo).upper()))
    return out, (arb(den_lo.lower()), arb(den_hi.upper())), m


# ----------------------------------------------------------------- (E) floats
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
    u = np.arange(np.log(5.0 ** (-n - 40)), np.log(1e12), h)
    s = np.exp(u)
    r, zs = _rn(s, n)
    wgt = s ** (1 - al) * r
    J = np.trapezoid(wgt, u)
    return J, {N: np.trapezoid(wgt * (1 - _eta(zs, n, N)), u) / J for N in Ns if N <= n}


def floats(nu08: float, cst, G: arb):
    out = {}
    with np.errstate(all="ignore"):
        for al, Ns, ns in ((0.8, (1, 2, 10), (10, 20, 30, 40, 50)), (0.9, (1, 6), (10, 20, 30, 40, 50))):
            rho = 5**al / 3
            Js, ds = {}, {}
            for n in sorted(set(ns) | {70, 80}):
                Js[n], ds[n] = float_deltas(al, n, Ns)
            Cp = (Js[80] - rho**10 * Js[70]) / (1 - rho**10)
            Dinf = (Js[70] - Cp) / rho**70
            nu = Cp / Dinf
            out[al] = {"nu": nu, "rows": {}}
            for N in Ns:
                dl = ds[80][N]
                X = [(n, rho**n * (dl - ds[n][N])) for n in ns if n >= N]
                out[al]["rows"][N] = (dl, -nu * (1 - dl), X)
            # BM's float at 0.9, N = 6: rho^20 (delta^(60) - delta^(20))
            if al == 0.9:
                _, d60 = float_deltas(al, 60, (6,))
                out[al]["bm20"] = rho**20 * (d60[6] - ds[20][6])
        # DZ.1 gap (floats): rho^-k H_k increasing for k = N..N+12, gap to k = 70 below ebar_k
        al, N = 0.8, 2
        rho = 5**al / 3
        hk = {}
        for k in list(range(N, N + 13)) + [70]:
            J, d = float_deltas(al, k, (N,))
            hk[k] = J * (1 - d[N]) / rho**k
        inc = all(hk[k + 1] > hk[k] for k in range(N, N + 12))
        gap_ok = all(hk[70] - hk[k] <= float(ebar(arb("0.8"), N, k, cst, G).upper()) for k in range(N, N + 13))
        out["gap"] = (inc, gap_ok, hk)
    return out


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem DZ: rho^n (delta_N - delta_N^(n)) -> -nu_alpha (1 - delta_N) (BL.6's limit part; 'independent of N' false; C_alpha < 0 on (gamma,1)). Not thrust.")
    symbolic_checks()
    transfer_spot_checks()
    ctx.prec = 80
    G0, _, _ = bi.sup_g()
    G = arb("0.6031")
    okb, _ = bi.bootstrap(G0, G)
    check(okb and bool(G0.upper() < G.lower()), "BI.4 re-run: sup s^(1-gamma) r_inf <= G = 0.6031")
    cst = dz1_constants(G)
    ctx.prec = 120
    a08 = arb("0.8")
    rho = arb(5) ** a08 / 3
    gm = cst["gm"]
    eb = {(N, k): ebar(a08, N, k, cst, G) for N in NS_CERT for k in (N, 20, 40, 60)}
    for N in (1, 10):
        print(f"   DZ.1 alpha = 0.8, N = {N}: ebar_N <= {float(eb[(N, N)].upper()):.4e}, ebar_40 <= {float(eb[(N, 40)].upper()):.4e}, "
              f"rho^60 ebar_60 <= {float((rho ** 60 * eb[(N, 60)]).upper()):.4e}")
    # the two ratios that drive DZ.1/DZ.2: rho 5^(2 gamma - alpha) = 3 and rho/3 < 1, 1/(5 rho) < 1/rho
    okr = True
    for al in ("0.69", "0.75", "0.8", "0.9", "0.99"):
        a = arb(al)
        rr = arb(5) ** a / 3
        okr &= bool(abs(rr * arb(5) ** (2 * gm - a) - 3).upper() < arb("1e-30")) and bool((rr / 3).upper() < 1) \
            and bool((2 * gm - a).lower() > 0) and bool((1 + a - 2 * gm).lower() > 0)
    check(okr, "rho 5^(2gamma-alpha) = 3, rho/3 < 1, 2gamma - alpha > 0, 1 + alpha - 2gamma > 0 (alpha = 0.69..0.99): rho^k ebar_k -> 0 for every alpha in (gamma,1)")
    RESULT["dz1"] = {"K1": cst["K1"], "m0": cst["m0"], "ebar": eb}

    # DZ.2 / DZ.4 at alpha = 0.8
    A = iv(*A08)
    Cal = iv(*C08)
    nu = iv(*NU08)
    kap = 2 * (arb.pi() * a08).sin() / arb.pi()
    Dinf = A / kap
    Cp = Cal / kap
    check(bool(abs(Cal / A - nu).upper() < arb("1e-20")), "AZ (3) as stated: nu_0.8 = C_0.8/A_0.8; D'_oo = A/kappa, C' = C/kappa")
    t1 = time.time()
    br, den, ncell = delta_brackets()
    print(f"   delta_N^({N_LEVEL})(0.8) brackets ({ncell} cells, {time.time() - t1:.1f}s); J_{N_LEVEL} in [{float(den[0]):.6f}, {float(den[1]):.6f}]")
    # consistency: rho^-n J_n vs D'_oo + rho^-n C' (AZ (1))
    Dk_pred = Dinf + rho ** (-N_LEVEL) * Cp
    check(bool(den[0].lower() <= (rho ** N_LEVEL * Dk_pred).upper() + arb("1e-12")) and bool((rho ** N_LEVEL * Dk_pred).lower() <= den[1].upper()),
          f"J_{N_LEVEL} bracket contains rho^{N_LEVEL} (D'_oo + rho^-{N_LEVEL} C') = {float((rho ** N_LEVEL * Dk_pred).mid()):.6f} (AZ (1), consistency)")
    # DZ.2's identity gives the ceiling X_n <= (1 - delta_N^(n)) (|C'| + E_n^-) / D'_n (AZ (3): E_60^- <= 6.7e-34)
    nL = N_LEVEL
    Dn_lo = Dinf + rho ** (-nL) * (Cp - arb(E60M))
    bk = arb(BK08) * rho ** (-nL)
    tab = {}
    okbk = True
    for N in NS_CERT:
        lo, hi = br[N]
        up = rho ** (-nL) * (1 - lo) * (abs(Cp) + arb(E60M)) / Dn_lo
        dN = hull(lo, hi + up)
        okbk &= bool((hi + up).upper() <= (hi + bk).upper())
        Cd = hull((-nu * (1 - dN.upper())).lower(), (-nu * (1 - dN.lower())).upper())
        tab[N] = (lo, hi, dN, Cd)
        print(f"   N = {N:2d}: delta_N^({nL}) in [{float(lo):.6f}, {float(hi):.6f}], delta_N in [{float(dN.lower()):.6f}, {float(dN.upper()):.6f}], "
              f"C_delta(0.8,N) in [{float(Cd.lower()):.6f}, {float(Cd.upper()):.6f}]")
    check(okbk, f"DZ.2's ceiling rho^-{nL} (1-delta^({nL}))(|C'|+E^-)/D'_{nL} <= {float(up.upper()):.3e} is inside BK/BN's stated 4.6644 rho^-{nL} = {float(bk.upper()):.3e}")
    dis = all(bool(tab[NS_CERT[i]][3].upper() < tab[NS_CERT[i + 1]][3].lower()) for i in range(len(NS_CERT) - 1))
    check(dis, "DZ.4: the enclosures of C_delta(0.8, N), N = 1, 2, 3, 6, 10, are disjoint and increasing: BL.6's 'independent of N' is false")
    check(bool(tab[10][3].lower() <= BL_FLOAT_08 <= tab[10][3].upper()), "BL.6's float C_delta(0.8) ~ 0.6025 (N = 10) lies in the certified C_delta(0.8, 10)")
    check(all(bool(arb(BN_C08[N]).upper() <= tab[N][3].lower()) for N in (1, 2) if N in tab)
          and bool(tab[10][3].upper() <= arb(BK08).lower()) and bool(-nu.lower() > tab[10][3].upper()),
          "DZ.3 consistency at 0.8: BN.5's c_1, c_2 <= C_delta(N) <= 4.6644; C_delta(N) < -nu_0.8 = lim_N C_delta(N)")
    # explicit rate for k >= 60 (AZ (3): E_60^+ <= 7.22e-23, E_60^- <= 6.7e-34; E_m <= E_60 for m >= 60)
    E = arb(E60P)
    rate = {}
    for N in NS_CERT:
        k = 60
        Dk_lo = Dinf + rho ** (-k) * (Cp - arb(E60M))
        bnd = ((1 - tab[N][2].lower()) * (abs(Cp) * rho ** (-k) * (abs(Cp) + E) / Dinf + E) + rho ** k * eb[(N, 60)]) / Dk_lo
        rate[N] = bnd
    print("   DZ.2 explicit rate at 0.8: |X_k - C_delta(N)| <= " + ", ".join(f"{float(rate[N].upper()):.3e} (N={N})" for N in NS_CERT) + " for every k >= 60")
    check(all(bool(rate[N].upper() < arb("1e-5")) for N in NS_CERT), "DZ.2 at 0.8: |X_k - C_delta(0.8,N)| < 1e-5 for every k >= 60, N <= 10 (bound decreasing in k)")
    RESULT["dz4"] = {"tab": tab, "rate": rate, "nu": nu, "ncell": ncell}

    # floats
    fl = floats(float(nu.mid()), cst, G)
    RESULT["floats"] = fl
    f8 = fl[0.8]
    print(f"   floats alpha = 0.8: nu ~ {f8['nu']:.7f} (AZ: {float(nu.mid()):.7f})")
    okf = abs(f8["nu"] - float(nu.mid())) < 1e-5
    for N, (dl, pred, X) in f8["rows"].items():
        print(f"     N = {N:2d}: delta_N ~ {dl:.6f}, -nu(1-delta_N) ~ {pred:.6f}; X_n = " + ", ".join(f"{x:.5f} (n={n})" for n, x in X))
        okf &= abs(X[-2][1] - pred) < 2e-3 and bool(tab[N][3].lower() - 1e-4 <= pred <= tab[N][3].upper() + 1e-4)
    f9 = fl[0.9]
    print(f"   floats alpha = 0.9: nu ~ {f9['nu']:.7f}; BM's rho^20 (delta_6^(60) - delta_6^(20)) ~ {f9['bm20']:.4f}")
    for N, (dl, pred, X) in f9["rows"].items():
        print(f"     N = {N:2d}: delta_N ~ {dl:.6f}, -nu(1-delta_N) ~ {pred:.6f}; X_n = " + ", ".join(f"{x:.5f} (n={n})" for n, x in X))
        okf &= abs(X[-1][1] - pred) < 2e-4
    okf &= abs(f9["bm20"] - BM_FLOAT_09) < 5e-4 and abs(f9["rows"][6][1] - BM_FLOAT_09) < 2e-3 and f9["nu"] < 0
    check(okf, "floats: X_n -> -nu(1 - delta_N) at 0.8 (N = 1, 2, 10) and 0.9 (N = 1, 6); BL's 0.6025 and BM's 0.3643 reproduced; nu_0.9 < 0")
    inc, gap_ok, _ = fl["gap"]
    check(inc and gap_ok, "floats (alpha = 0.8, N = 2, k = 2..14): rho^-k H_k increases and rho^-70 H_70 - rho^-k H_k <= ebar_k (DZ.1)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all DZ checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
