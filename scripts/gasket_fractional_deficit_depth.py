#!/usr/bin/env python3
"""Theorem BR: the depth dependence of the corner deficit.  At every level n and every alpha in (0,1), delta_N^(n) strictly
DEcreases in the depth N (0 <= N < n); below log3/log5 the constant C_N = lim rho^-n (1 - delta_N^(n)) = P_N/D_inf of
Theorem BQ strictly increases in N and is of exact order rho^-N: certified two-sided brackets for rho^N C_N, uniform in N >= 1,
and for liminf/limsup_N rho^N C_N, at alpha = 0.01, 0.1, 0.3, 0.45.

Setting (GASKET.md, Theorems BG, BH, BJ, BQ): p(z) = z(z+5), alpha_1 = (z+4)/((z+2)(z+5)), beta = (z+6)/((z+2)(z+5)),
eta_0 = -1/2, eta_N(z) = alpha_1(z) + beta(z) eta_(N-1)(p z), eta = lim eta_N, rho = 5^a/3, gamma = log3/log5,
R(tau) = lim (3/5)^n r_n(5^-n tau), psi(t) = lim p^k(5^-k t), and for N >= 0, n >= N (BG.1, BJ.4, as in BQ's setting)
    1 - delta_N^(n) = D_n^-1 int s^-a r_n(s) eta_N(p^(n-N) s) ds,      P_N = int tau^-a R(tau) eta_N(psi(5^-N tau)) dtau.

Theorem BR (checked here; exit 1 on any failure):
  (1) [ceiling, proved by induction, every N >= 0] eta_N(z) <= 2/(2+z) on [0, inf)  (was quoted from BH.2/BJ.4 in BQ).
  (2) [depth step, every N >= 0, w >= 0]  eta_(N+1)(w) - eta_N(p(w)) >= w/(w^2+5w+2)  (sympy: alpha_1 - (1-beta) 2/(2+p) =
      w/(w^2+5w+2)).
  (3) [every level, every alpha in (0,1)] for 0 <= N < n:  delta_N^(n) - delta_(N+1)^(n) >= D_n^-1 int s^-a r_n(s) w/(w^2+5w+2) ds
      > 0 with w = p^(n-N-1)(s).  So delta_N^(n) strictly decreases in N; at N = 0 it is BL.0's 3/2.
  (4) [alpha < gamma] P_(N+1) > P_N (same pointwise step at w = psi(5^-(N+1) tau)), so C_N = P_N/D_inf strictly increases.
      Scaling tau = 5^N t: rho^N P_N = int t^-a Q_N(t) eta_N(psi(t)) dt with Q_N = (5/3)^N R(5^N t) in BJ.2's bracket, so for
      every N >= 1   L_1 <= rho^N P_N <= g_hi M,   M = int t^(gamma-1-a) eta(psi(t)) dt <= 2^c pi/sin(pi c), c = gamma - a,
      L_1 = int t^-a eta_1(psi t)/(t^(1-gamma)/g_lo + 9/5) dt  (eta_N >= eta_1, increasing in N; 3(3/5)^N <= 9/5),
      and g_lo M <= liminf rho^N P_N <= limsup rho^N P_N <= g_hi M (monotone convergence).  C_N is of exact order rho^-N.
      Moreover Q_N(t) = lim_m (3/5)^m r_(m+N)(5^-m t) is nondecreasing in N (AQ (3)) and eta_N(psi t) strictly increases in N,
      so rho^N C_N strictly increases in N and converges: c_inf(alpha) = lim rho^N C_N exists in [g_lo M/D_hi, g_hi M/D_lo].
  (5) [certified, python-flint arb, outward] M, L_1, D_inf at alpha = 0.01, 0.1, 0.3, 0.45; the brackets for rho^N C_N.
Assumptions used as stated (not re-run): BG.1, BG.3, BJ.4 (the deficit formulas), BJ.2 (Q_N bracket), BJ.3 (g_lo, g_hi),
BJ.6's psi(y) <= y e^(y/10), BD.4 (psi_k increasing to psi, psi(5u) = p(psi(u))), BH.2 (eta = lim eta_N, eta = eta_N +
prod beta (1/2 + eta o p^N)), AQ (3), BI.4 (G = 0.6031), BQ.1-BQ.3 (D_inf routine, C_1).  Floats: consistency only.
Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr
from pathlib import Path

import sympy as sp
from flint import arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasket_fractional_deficit_rate_small_alpha as bq  # noqa: E402

FAILS: list[str] = []
ALPHAS = ("0.01", "0.1", "0.3", "0.45")
G_LO, G_HI = "0.600945014", "0.600963487"
# BQ.3's certified C_1 brackets (consistency only)
C1_BQ = {"0.01": ("0.0078312", "0.0078363"), "0.1": ("0.0606923", "0.0607522"),
         "0.3": ("0.0933210", "0.0934143"), "0.45": ("0.0718892", "0.0719611")}
K_ORB = 14


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    z, w, e, t, c = sp.symbols("z w e t c", positive=True)
    D = (z + 2) * (z + 5)
    al1, beta = (z + 4) / D, (z + 6) / D
    check(sp.simplify(al1 + beta - 2 / (2 + z)) == 0, "alpha_1 + beta = 2/(2+z)")
    # (1) induction: eta_(N-1) <= 2/(2+pz), beta > 0  =>  eta_N <= alpha_1 + beta 2/(2+pz) <= 2/(2+z)
    pz = z * (z + 5)
    gap = sp.factor(sp.simplify(2 / (2 + z) - (al1 + beta * 2 / (2 + pz))))
    check(sp.simplify(gap - z * (z + 6) / ((z + 2) * (z**2 + 5 * z + 2))) == 0,
          f"2/(2+z) - (alpha_1 + beta 2/(2+p)) = {gap} >= 0")
    num, den = sp.fraction(sp.together(gap))
    check(all(cf >= 0 for cf in sp.Poly(sp.expand(num), z).all_coeffs()) and
          all(cf >= 0 for cf in sp.Poly(sp.expand(den), z).all_coeffs()),
          "  numerator and denominator have nonnegative coefficients (so the gap is >= 0 on z >= 0); base eta_0 = -1/2 <= 1: "
          "eta_N <= 2/(2+z) for every N >= 0 (BR.1)")
    # (2) depth step
    Dw = (w + 2) * (w + 5)
    a1w, bw = (w + 4) / Dw, (w + 6) / Dw
    pw = w * (w + 5)
    step = sp.simplify(a1w - (1 - bw) * 2 / (2 + pw))
    check(sp.simplify(step - w / (w**2 + 5 * w + 2)) == 0,
          "alpha_1(w) - (1 - beta(w)) 2/(2+p(w)) = w/(w^2+5w+2) (BR.2)")
    check(sp.factor(sp.simplify(1 - bw)) == sp.factor((w**2 + 6 * w + 4) / Dw),
          "1 - beta(w) = (w^2+6w+4)/((w+2)(w+5)) > 0, so eta_(N+1)(w) - eta_N(pw) = alpha_1 - (1-beta) eta_N(pw) >= w/(w^2+5w+2)")
    # monotone in N at fixed z: eta_(N+1)-eta_N = prod beta * (1/2 + eta_1(x)), x = p^N z, via BH.2's tail form
    x = sp.symbols("x", positive=True)
    Dx = (x + 2) * (x + 5)
    tailstep = sp.simplify((sp.Rational(1, 2) + (x + 4) / Dx + (x + 6) / Dx * e) - (x + 6) / Dx * (sp.Rational(1, 2) + e))
    check(sp.simplify(tailstep - (sp.Rational(1, 2) + 1 / (2 * (x + 5)))) == 0,
          "with eta = alpha_1 + beta eta o p:  (1/2 + eta(x)) - beta(x)(1/2 + eta(px)) = 1/2 + 1/(2(x+5)) > 0, so "
          "eta_N = eta - prod_(i<N) beta (1/2 + eta o p^N) increases in N at fixed z (eta_N >= eta_1 = 1/(2(z+5)) for N >= 1)")
    # 1 - eta_0 = 3/2 and the depth identity at N = 0 reproduces BL.0's 3/2
    check(1 - sp.Rational(-1, 2) == sp.Rational(3, 2), "N = 0: eta_0 = -1/2 gives delta_0^(n) = 3/2 (BL.0)")
    # scaling
    a, N = sp.symbols("a N", positive=True)
    rho = 5**a / 3
    tau = 5**N * t
    fac = sp.simplify(sp.powsimp(sp.expand_power_base(rho**N * tau**(-a) * sp.diff(tau, t) * (sp.Rational(3, 5))**N, force=True),
                                 force=True))
    check(sp.simplify(fac - t**(-a)) == 0,
          "rho^N tau^-a dtau = t^-a (3/5)^N dt at tau = 5^N t, so rho^N P_N = int t^-a [(5/3)^N R(5^N t)] eta_N(psi(t)) dt "
          "(psi(5^-N tau) = psi(t))")
    # closed form of the M ceiling: int t^(c-1) 2/(2+t) = 2^c pi/sin(pi c)
    cc = sp.Rational(1, 3)
    val = sp.integrate(t**(cc - 1) * 2 / (2 + t), (t, 0, sp.oo))
    check(sp.simplify(val - 2**cc * sp.pi / sp.sin(sp.pi * cc)) == 0,
          "int_0^inf t^(c-1) 2/(2+t) dt = 2^c pi/sin(pi c) (sympy at c = 1/3; Beta integral, 0 < c < 1)")
    # psi_K <= psi <= p^K(5^-K t e^(5^-K t/10)): p increasing
    check(sp.simplify(sp.diff(z * (z + 5), z) - (2 * z + 5)) == 0,
          "p' = 2z+5 > 0: psi = p^K o psi o 5^-K is sandwiched by p^K(5^-K t) and p^K(5^-K t e^(5^-K t/10)) (BJ.6's psi bound)")


# ----------------------------------------------------------------- (B) exact rationals
def exact_checks():
    zs = (Fr(0), Fr(1, 1000), Fr(1, 10), Fr(1), Fr(7, 2), Fr(30), Fr(10**6))
    ok1 = all(bq.eta_rec(N, zz) <= Fr(2) / (2 + zz) for N in range(0, 8) for zz in zs)
    check(ok1, "exact rationals N = 0..7, z in {0, 1/1000, 1/10, 1, 7/2, 30, 1e6}: eta_N(z) <= 2/(2+z)")
    ok2 = all(bq.eta_rec(N + 1, ww) - bq.eta_rec(N, ww * (ww + 5)) >= ww / (ww * ww + 5 * ww + 2)
              for N in range(0, 7) for ww in zs)
    check(ok2, "exact rationals N = 0..6: eta_(N+1)(w) - eta_N(p w) >= w/(w^2+5w+2) (BR.2)")
    ok3 = all(bq.eta_rec(N + 1, zz) > bq.eta_rec(N, zz) for N in range(0, 7) for zz in zs)
    check(ok3, "exact rationals N = 0..6: eta_(N+1)(z) > eta_N(z) at fixed z (eta_N increases to eta)")


# ----------------------------------------------------------------- (C) arb
def _psi_bracket(t, K=K_ORB):
    f = arb(5) ** (-K)
    u = t * f
    lo, hi = u, u * (u / 10).exp()
    for _ in range(K):
        lo, hi = lo * (lo + 5), hi * (hi + 5)
    return lo, hi


def _eta_hi(z, K=K_ORB):
    """eta(z) <= eta_K(z) + prod_(i<K) beta(p^i z) (1/2 + 2/(2+p^K z))  (BH.2 tail form, eta <= 2/(2+z))."""
    o = bq._orbit(z, K)
    prod = arb(1)
    for i in range(K):
        zz = o[i]
        prod *= (zz + 6) / ((zz + 2) * (zz + 5))
    return bq._eta(z, K) + prod * (arb(1) / 2 + 2 / (2 + o[K]))


def M_bracket(al, q="1.0003", tlo="1e-40", thi="1e4"):
    """[lo, hi] for M = int t^(gamma-1-a) eta(psi(t)) dt; integrand decreasing (eta decreasing, psi increasing, c < 1)."""
    a, gam = arb(al), bq._gam()
    c = gam - a
    pts = bq._grid(tlo, thi, q)
    up, low = [], []
    for t in pts:
        plo, phi = _psi_bracket(t)
        tc = t ** (c - 1)
        up.append(tc * _eta_hi(plo))
        low.append(tc * bq._eta(phi, K_ORB))          # eta >= eta_K, eta decreasing, psi <= phi
    lo = hi = arb(0)
    for i in range(len(pts) - 1):
        wdt = pts[i + 1] - pts[i]
        lo += low[i + 1] * wdt
        hi += up[i] * wdt
    e0, T = pts[0], pts[-1]
    # head: eta <= 1.  tail: psi(t)/t is nondecreasing (psi_k(t)/t = prod_(j<k) (p^j(5^-k t) + 5)/5), so for t >= T
    # eta(psi t) <= 2/psi(t) <= 2 T/(psi_K(T) t) and int_T^inf t^(c-2) = T^(c-1)/(1-c)
    hi += e0 ** c / c + 2 * T ** c / (_psi_bracket(T)[0] * (1 - c))
    return arb(lo.lower()), arb(hi.upper())


def L1_lower(al, glo, q="1.0005", tlo="1e-30", thi="1e4"):
    """lower bound for int t^-a eta_1(psi t)/(t^(1-gamma)/g_lo + 9/5) dt (decreasing integrand, right endpoints)."""
    a, gam = arb(al), bq._gam()
    pts = bq._grid(tlo, thi, q)
    vals = []
    for t in pts:
        _, phi = _psi_bracket(t)
        vals.append(t ** (-a) / (2 * (phi + 5)) / (t ** (1 - gam) / glo + arb(9) / 5))
    lo = arb(0)
    for i in range(len(pts) - 1):
        lo += vals[i + 1] * (pts[i + 1] - pts[i])
    return arb(lo.lower())


def certify(al):
    ctx.prec = 64
    glo, ghi, G = arb(G_LO), arb(G_HI), arb(bq.G_BI4)
    Mlo, Mhi = M_bracket(al)
    L1 = L1_lower(al, glo)
    Dlo, Dhi, _ = bq.Dinf(al, 60, G, glo, ghi)
    c = bq._gam() - arb(al)
    closed = arb(2) ** c * arb.pi() / (arb.pi() * c).sin()
    return dict(M=(Mlo, Mhi), L1=L1, D=(Dlo, Dhi), closed=closed,
                every=(L1 / Dhi, ghi * Mhi / Dlo), lim=(glo * Mlo / Dhi, ghi * Mhi / Dlo))


# ----------------------------------------------------------------- (D) floats
def floats_depth(al="0.45", n=6, Ns=(1, 2, 3, 4, 5)):
    """delta_N^(n) at one level n for several depths N (quadrature; consistency with BR.3 only)."""
    import mesh_maxwell_limit as ml
    from mpmath import mp, mpf, quad
    mp.dps = 20
    a = mpf(al)
    X = mpf(10) ** 12
    pts = [0] + [mpf(5) ** (-n) * mpf(10) ** j for j in range(-6, 40, 2) if mpf(5) ** (-n) * mpf(10) ** j < X] + [X]

    def pk(s, k):
        for _ in range(k):
            s = s * (s + 5)
        return s

    def etaN(zv, N):
        zs = [zv]
        for _ in range(N):
            zs.append(zs[-1] * (zs[-1] + 5))
        e = mpf(-1) / 2
        for k in range(N - 1, -1, -1):
            zz = zs[k]
            Dd = (zz + 2) * (zz + 5)
            e = (zz + 4) / Dd + (zz + 6) / Dd * e
        return e

    f = lambda s: s ** (-a) * ml.rec_all(n, s)[0]
    tail = X ** (-a) / a - 2 * X ** (-1 - a) / (1 + a) + 6 * X ** (-2 - a) / (2 + a)
    D = quad(f, pts) + tail
    return [(N, float(1 - quad(lambda s: f(s) * etaN(pk(s, n - N), N), pts) / D)) for N in Ns]


def main() -> int:
    t0 = time.time()
    print("Theorem BR: delta_N^(n) decreases in the depth N at every level; below gamma rho^N C_N increases in N to a finite "
          "positive limit. Not thrust.")
    symbolic_checks()
    exact_checks()
    for al in ALPHAS:
        r = certify(al)
        Ml, Mh = r["M"]
        Dl, Dh = r["D"]
        ev, lm = r["every"], r["lim"]
        rho = arb(5) ** arb(al) / 3
        print(f"   alpha={al}: M in [{float(Ml.lower()):.7f},{float(Mh.upper()):.7f}] (closed ceiling 2^c pi/sin(pi c) = "
              f"{float(r['closed'].mid()):.6f}), L_1 >= {float(r['L1'].lower()):.6f}, D_inf in [{float(Dl.lower()):.6f},"
              f"{float(Dh.upper()):.6f}]")
        check(bool(Ml.lower() > 0) and bool(Mh < r["closed"]) and bool(r["L1"] > 0),
              f"alpha={al}: 0 < M, M below the closed ceiling 2^c pi/sin(pi c), L_1 > 0")
        check(bool(ev[0] > 0) and bool(ev[0] < ev[1]),
              f"alpha={al}: every N >= 1: rho^N C_N in [{float(ev[0].lower()):.6f}, {float(ev[1].upper()):.6f}]; "
              f"liminf/limsup_N rho^N C_N in [{float(lm[0].lower()):.6f}, {float(lm[1].upper()):.6f}]")
        c1lo, c1hi = (arb(v) for v in C1_BQ[al])
        check(bool(rho * c1lo > ev[0]) and bool(rho * c1hi < ev[1]),
              f"alpha={al}: consistency, BQ.3's rho C_1 in [{float((rho * c1lo).lower()):.6f},{float((rho * c1hi).upper()):.6f}] "
              f"lies in the every-N bracket")
        check(bool(rho * c1hi < lm[0]),
              f"alpha={al}: rho C_1 <= {float((rho * c1hi).upper()):.6f} < c_inf = lim_N rho^N C_N "
              f"(>= {float(lm[0].lower()):.6f}): the increase of rho^N C_N is genuine; every N >= 1: rho^N C_N in "
              f"[{float((rho * c1lo).lower()):.6f}, {float(ev[1].upper()):.6f}] with BQ.3's C_1")
    fl = floats_depth()
    print("   floats alpha=0.45, n=6 (N, delta_N^(6)): " + "; ".join(f"{N}: {d:.7f}" for N, d in fl))
    check(all(fl[i][1] > fl[i + 1][1] for i in range(len(fl) - 1)), "floats: delta_N^(6)(0.45) decreasing in N (consistency only)")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BR checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
