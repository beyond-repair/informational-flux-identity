#!/usr/bin/env python3
r"""Theorem BX: the log-periodic profile F = theta of Theorems BT-BW is a convolution of one explicit kernel with a positive
1-periodic "spectral" measure m of total mass c_R(gamma)/2, so every Fourier mode factorises exactly,
    theta_k = (pi/ln5) m_k / sin(pi (gamma + i omega k)),   omega = 2 pi/ln5,   m_k = int_0^1 e^(-2 pi i k t) dm(t),
which gives the sharp bound |theta_k| <= <F> sin(pi gamma)/|sin(pi(gamma + i omega k))| (no polynomial loss, unlike BV.3)
and reduces Conjecture BV.6 to the behaviour of m_k; Kept Failure BX.5: m is NOT supported on <= 16 points per period, so
Q_inf is not a superposition of <= 16 self-similar pole series (in particular not a single one).

Setting (GASKET.md, Theorems AQ, AU, BN, BS-BW): r_inf = lim r_k, Q_inf(u) = lim_M (3/5)^M r_inf(5^-M u),
gamma = log3/log5, F(u) = u^(1-gamma) Q_inf(u) = theta(u), F(5u) = F(u), t = log_5 u,
theta_k = int_0^1 F(5^t) e^(-2 pi i k t) dt = a_k - i b_k (AU's normalisation).

Theorem BX (checked here; exit 1 on any failure):
  (1) [Stieltjes representation] Q_inf(u) = int_(0,inf) dnu(lam)/(lam+u) with nu >= 0, no atom at 0, no constant term,
      and nu(5A) = 3 nu(A).  Proof: BV.2's locally uniform limit of the Stieltjes functions (3/5)^M r_inf(5^-M z) keeps
      Im z * Im f(z) <= 0 and f >= 0 on (0,inf), so Q_inf is Stieltjes (Assumption BX: the Nevanlinna characterisation of
      Stieltjes functions); BS.1's Q_inf(5u) = (3/5) Q_inf(u) and uniqueness of the representation give a/5 = 3a/5,
      b = 3b/5, nu(5A) = 3 nu(A) (sympy).
  (2) [factorisation] sigma = lam^-gamma nu, written in t = log_5 lam, is 1-periodic (3 * 5^-gamma = 1); its restriction m
      to one period is a positive measure and F(tau) = int K(t - tau) dsigma(t) with K(x) = 5^(gamma x)/(5^x + 1) (sympy).
      Fubini and int_0^inf v^(a-1)/(1+v) dv = pi/sin(pi a) (0 < Re a < 1) give theta_k = (pi/ln5) m_k/sin(pi(gamma+i omega k)).
      k = 0: m(T) = <F> ln5 sin(pi gamma)/pi = c_R(gamma)/2 (AU (1)).
  (3) [sharp bound; derived] |m_k| <= m_0, |sin(pi(gamma+iy))|^2 = sin^2(pi gamma) + sinh^2(pi y) (sympy), so
      |theta_k| <= <F> sin(pi gamma)/sqrt(sin^2 pi gamma + sinh^2(2 pi^2 k/ln5)) <= (pi c_R(gamma)/ln5) e^(-2 pi^2 k/ln5)/(1-e^(-4 pi^2 k/ln5)),
      so limsup |theta_k| e^(2 pi^2 k/ln5) <= pi c_R(gamma)/ln5 = 1.00949802... (BV.3: e G*(1+4 pi k/ln5), unbounded in k).
      BV.6 holds iff |m_k| converges to a nonzero limit.
  (4) [certified, BW's routine re-run unchanged] rho_k = |m_k|/m_0 for k <= 20 enclosed; rho_k in [0.60, 0.92];
      |theta_k| <= the bound of (3) for every k <= 20; the Hermitian Toeplitz matrix T_17 = (m_(i-j))_(i,j<17) is
      certified positive definite (exact-midpoint LDL* in arb at 1200 bits with shift 2||E||_F, E the enclosure radii).
  Kept Failure BX.5 (certified witness): "nu is carried by r <= 16 geometric series lam_j 5^Z", i.e.
      Q_inf(u) = sum_(j<=r) c_j sum_(n in Z) 3^n/(lam_j 5^n + u) (r = 1: one self-similar pole series) is false:
      such an m has <= r atoms, so rank T_(r+1) <= r, contradicting T_17 > 0; for r = 1 directly rho_1 <= 0.7697 < 1.
Assumptions used as stated (not re-run): AQ (3), BN.1, BS.1, BT.2, BU.1, BV.2 (holomorphic continuation of Q_inf by a
locally uniform limit), AU (1); Assumption BX (Nevanlinna characterisation of Stieltjes functions, Schilling-Song-
Vondracek, Bernstein Functions, ch. 7); BW's certified routine (AU (3) unchanged with BW's parameters), imported unchanged.
Floats: consistency only (toy single-series check of the normalisation).  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import sympy as sp
from flint import acb, arb, ctx

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gasket_fractional_profile_modes as bw  # noqa: E402

bv = bw.bv
FAILS: list[str] = []
KMAX = bw.KMAX          # 20
N_TOEP = 17             # certified size of the positive definite Toeplitz section
RHO_BAND = ("0.60", "0.92")
LIMSUP_CONST = "1.0095"    # upper bound for pi c_R(gamma)/ln5 (checked below)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    a, b, c = sp.symbols("a b c", nonnegative=True)
    sol = sp.solve([sp.Eq(a / 5, 3 * a / 5), sp.Eq(b, 3 * b / 5)], [a, b], dict=True)
    check(sol == [{a: 0, b: 0}], "BX.1: Q(5u) = (3/5)Q(u) with Q = a/u + b + int dnu/(lam+u) forces a = b = 0 "
                                 "(and nu(5A) = 3 nu(A) by uniqueness)")
    g = sp.log(3) / sp.log(5)
    check(sp.simplify(3 * sp.Integer(5) ** (-g) - 1) == 0, "BX.2: 3 * 5^-gamma = 1, so sigma = lam^-gamma nu is invariant under lam -> 5 lam")
    t, tau = sp.symbols("t tau", real=True)
    G = sp.symbols("G", positive=True)
    lhs = 5 ** ((1 - G) * tau) * 5 ** (G * t) / (5 ** t + 5 ** tau)
    K = 5 ** (G * (t - tau)) / (5 ** (t - tau) + 1)
    check(sp.simplify(sp.powsimp(sp.expand_power_base(lhs / K, force=True), force=True) - 1) == 0,
          "BX.2: 5^((1-gamma)tau) lam^gamma/(lam + 5^tau) with lam = 5^t equals K(t - tau) = 5^(gamma(t-tau))/(5^(t-tau)+1)")
    v = sp.symbols("v", positive=True)
    s = sp.Rational(3, 10)
    I1 = sp.integrate(v ** (s - 1) / (1 + v), (v, 0, sp.oo))
    check(sp.simplify(I1 - sp.pi / sp.sin(sp.pi * s)) == 0, "BX.2: int_0^inf v^(a-1)/(1+v) dv = pi/sin(pi a) (sympy at a = 3/10; "
                                                              "analytic in 0 < Re a < 1, so the identity holds there)")
    x, y = sp.symbols("x y", real=True)
    mod2 = sp.expand(sp.expand_complex(sp.sin(x + sp.I * y) * sp.sin(x - sp.I * y)))
    check(sp.simplify(mod2 - (sp.sin(x) ** 2 + sp.sinh(y) ** 2)) == 0, "BX.3: |sin(x+iy)|^2 = sin^2 x + sinh^2 y")
    check(sp.simplify((sp.exp(y) / (2 * sp.sinh(y))).rewrite(sp.exp) - 1 / (1 - sp.exp(-2 * y))) == 0,
          "BX.3: e^y/(2 sinh y) = 1/(1 - e^(-2y)), the factor in the normalised bound")
    # Toeplitz rank lemma: r atoms => rank T_(r+1) <= r (symbolic, r = 2)
    w1, w2 = sp.symbols("w1 w2", positive=True)
    p1, p2 = sp.symbols("p1 p2", real=True)
    z1, z2 = sp.exp(-2 * sp.pi * sp.I * p1), sp.exp(-2 * sp.pi * sp.I * p2)
    T = sp.Matrix(3, 3, lambda i, j: w1 * z1 ** (i - j) + w2 * z2 ** (i - j))
    check(sp.simplify(T.det()) == 0, "BX.5: for m = w1 delta_p1 + w2 delta_p2 the 3x3 Toeplitz section is singular (rank <= number of atoms)")


# ----------------------------------------------------------------- (B) toy normalisation check (floats)
def toy_checks():
    """nu = sum_n 3^n delta_(5^n): m = delta_0 (mass 1), so theta_k must equal (pi/ln5)/sin(pi(gamma+i omega k))."""
    from mpmath import mp, mpf, mpc, quad, sin, pi, log, exp
    mp.dps = 30
    g = log(3) / log(5)
    w = 2 * pi / log(5)

    def F(tau):
        u = mpf(5) ** tau
        return u ** (1 - g) * sum(mpf(3) ** n / (mpf(5) ** n + u) for n in range(-200, 120))
    ok = True
    for k in (0, 1, 2):
        th = quad(lambda t: F(t) * exp(-2j * pi * k * t), [0, 0.25, 0.5, 0.75, 1])
        pred = pi / log(5) / sin(pi * mpc(g, w * k))
        ok = ok and abs(th - pred) < mpf("1e-18") * max(1, abs(pred) * 1e6)
        print(f"   toy k={k}: quadrature {mp.nstr(th, 12)}  formula {mp.nstr(pred, 12)}")
    check(ok, "BX.2 (floats, normalisation only): a single geometric series (m = delta_0) reproduces "
              "theta_k = (pi/ln5)/sin(pi(gamma+i omega k)) for k = 0, 1, 2")


# ----------------------------------------------------------------- (C) arb helpers
def consts():
    pi, l5 = arb.pi(), arb(5).log()
    g = arb(3).log() / l5
    return pi, l5, g, 2 * pi / l5


def bound_k(k, mean):
    """BX.3: |theta_k| <= <F> sin(pi gamma)/|sin(pi(gamma + i omega k))|."""
    pi, l5, g, w = consts()
    sg = (pi * g).sin()
    return mean * sg / (sg ** 2 + (pi * w * k).sinh() ** 2).sqrt()


def m_from_theta(k, th):
    pi, l5, g, w = consts()
    return th * l5 * (pi * acb(g, w * k)).sin() / pi


def table_checks():
    """Published AU (k <= 4) and BV.4 (k = 5..11) tables: rho_k < 1, inside the band, and the BX.3 bound beats BV.3's B_k."""
    ctx.prec = 300
    mean = arb(bw.AU_MEAN[0])
    tab = {k: (bv.AU_TABLE_A[k], bv.AU_TABLE_B[k]) for k in range(1, 5)}
    tab.update(bw.BV_TABLE)
    G = arb(bv.G_STAR_HI)
    m0 = m_from_theta(0, acb(mean)).real
    ok_rho, ok_bd = True, True
    for k, (a, b) in sorted(tab.items()):
        th = acb(arb(a), -arb(b))
        rho = abs(m_from_theta(k, th)) / m0
        ok_rho = ok_rho and bool(rho < 1) and bool(rho > arb(RHO_BAND[0])) and bool(rho < arb(RHO_BAND[1]))
        ok_bd = ok_bd and bool(abs(th) < bound_k(k, mean)) and bool(bound_k(k, mean) < bv.B_bound(k, G))
    check(ok_rho, f"BX.4 (published tables, k <= 11): rho_k = |m_k|/m_0 in [{RHO_BAND[0]}, {RHO_BAND[1]}] (< 1 as BX.3 requires)")
    check(ok_bd, "BX.3 vs published tables: |theta_k| <= <F> sin(pi gamma)/|sin(pi(gamma+i omega k))| < B_k (BV.3) for k <= 11")
    pi, l5, g, w = consts()
    lim = pi * arb(bw.AU_CRG[1]) / l5
    check(bool(lim < arb(LIMSUP_CONST)) and bool(pi * arb(bw.AU_CRG[0]) / l5 > arb("1.0094")),
          f"BX.3: pi c_R(gamma)/ln5 = 1.0094980205... in [1.0094, {LIMSUP_CONST}] (AU's c_R enclosure)")
    return ok_rho and ok_bd


# ----------------------------------------------------------------- (D) Toeplitz certificate
def ldl_pd(mid, n, shift):
    """Exact-midpoint Hermitian LDL* of T_n - shift*I in acb; True iff every pivot is certified > 0."""
    def T(i, j):
        d = i - j
        v = mid[d] if d >= 0 else mid[-d].conjugate()
        return v - shift if i == j else v
    L = [[acb(0)] * n for _ in range(n)]
    D = []
    for j in range(n):
        s = T(j, j)
        for p in range(j):
            s -= L[j][p] * L[j][p].conjugate() * D[p]
        if not bool(s.real > 0):
            return False, D
        D.append(s.real)
        for i in range(j + 1, n):
            s = T(i, j)
            for p in range(j):
                s -= L[i][p] * L[j][p].conjugate() * D[p]
            L[i][j] = s / D[j]
    return True, D


def toeplitz_certificate(m, n):
    mid = {k: acb(v.real.mid(), v.imag.mid()) for k, v in m.items()}
    E = arb(0)
    for i in range(n):
        for j in range(n):
            v = m[abs(i - j)]
            E += (v.real.rad() + v.imag.rad()) ** 2      # |m_k - mid_k| <= rad(Re m_k) + rad(Im m_k)
    E = E.sqrt()
    ok, D = ldl_pd(mid, n, 2 * E.upper())
    return ok, E, D


def toeplitz_toy_checks():
    """Sanity of the certificate on exact toy measures: three atoms give T_3 > 0 certified and T_4 not certified."""
    ctx.prec = 300
    pts, wts = (arb(0), arb(1) / 7, arb(3) / 5), (arb(1), arb(2), arb(3))
    m = {k: sum((acb(0, -2 * arb.pi() * k * p).exp() * w for p, w in zip(pts, wts)), acb(0)) for k in range(6)}
    ok3, _, _ = toeplitz_certificate(m, 3)
    ok4, _, _ = toeplitz_certificate(m, 4)
    check(ok3 and not ok4, "BX.5 certificate sanity: for a 3-atom measure T_3 is certified positive definite and T_4 is not")


# ----------------------------------------------------------------- main
def main() -> int:
    t0 = time.time()
    print("Theorem BX: theta_k = (pi/ln5) m_k/sin(pi(gamma+i omega k)), m >= 0 periodic of mass c_R(gamma)/2; sharp mode bound; "
          "Kept Failure BX.5: m has >= 17 support points per period. Not thrust.")
    symbolic_checks()
    toy_checks()
    table_checks()
    toeplitz_toy_checks()
    ctx.prec = 1200
    bw.configure()
    t1 = time.time()
    coef = bv.au_modes()
    print(f"   BW's routine ({bw.K_PER} cells, b = u = {bw.ORD}/{bw.ORD + 1}, {bw.DPS} digits, depth 5^-{bw.E_DEPTH}): {time.time() - t1:.0f}s")
    check(coef["nodes_ok"], "AQ (4) node enclosures of r_inf positive with relative width < 1e-100")
    mean = coef["mean"]
    pi, l5, g, w = consts()
    m = {0: m_from_theta(0, acb(mean))}
    for k in range(1, KMAX + 1):
        m[k] = m_from_theta(k, acb(coef["a"][k], -coef["b"][k]))
    m0 = m[0].real
    crg = 2 * (pi * g).sin() / pi * l5 * mean
    print(f"   m_0 = m(T) = {m0.str(60, radius=False)}  (radius {float(m0.rad()):.1e})")
    check(bool(abs(2 * m0 - crg) < arb("1e-110")) and bool(2 * m0 > arb(bw.AU_CRG[0])) and bool(2 * m0 < arb(bw.AU_CRG[1])),
          "BX.2: m(T) = c_R(gamma)/2, inside AU's 34-digit enclosure of c_R(gamma)")
    rows, band, bd_ok, beat = [], True, True, True
    G = arb(bv.G_STAR_HI)
    x = (-2 * pi ** 2 / l5).exp()
    for k in range(1, KMAX + 1):
        rho = abs(m[k]) / m0
        th_abs = abs(acb(coef["a"][k], -coef["b"][k]))
        bk = bound_k(k, mean)
        band = band and bool(rho > arb(RHO_BAND[0])) and bool(rho < arb(RHO_BAND[1]))
        bd_ok = bd_ok and bool(th_abs < bk)
        beat = beat and bool(bk < bv.B_bound(k, G))
        ph = m[k].arg() / (2 * pi)
        rows.append((k, rho, ph, bk))
        print(f"   k={k:2d}: rho_k = |m_k|/m_0 = {rho.str(12, radius=False):>16s} (rad {float(rho.rad()):.1e})  arg m_k/2pi = "
              f"{float(ph.mid()):+.10f}  BX.3 bound {float(bk.upper()):.4e}  B_k/bound {float((bv.B_bound(k, G) / bk).mid()):7.1f}"
              f"  |theta_k|e^(2pi^2k/ln5) <= {float((bk / x ** k).upper()):.5f}")
    check(band, f"BX.4: rho_k in [{RHO_BAND[0]}, {RHO_BAND[1]}] for every k = 1..{KMAX} (certified)")
    check(bd_ok, f"BX.3 vs BW.1: |theta_k| <= <F> sin(pi gamma)/|sin(pi(gamma + i omega k))| for every k = 1..{KMAX}")
    check(beat, f"BX.3 is sharper than BV.3: bound_k < B_k for every k = 1..{KMAX}")
    rho1 = abs(m[1]) / m0
    check(bool(rho1 < arb("0.7697")), f"Kept Failure BX.5, r = 1: rho_1 <= {float(rho1.upper()):.10f} < 1, so m is not one atom "
                                      "(Q_inf is not c sum_n 3^n/(lam_0 5^n + u))")
    ok, E, D = toeplitz_certificate(m, N_TOEP)
    dmin = min(D, key=lambda d: float(d.mid())) if D else arb(0)
    print(f"   T_{N_TOEP}: ||E||_F <= {float(E.upper()):.2e}; pivots of T - 2||E|| I certified positive: {ok}; "
          f"smallest pivot {float(dmin.mid()):.3e}")
    check(ok, f"BX.4 / Kept Failure BX.5: T_{N_TOEP} = (m_(i-j)) is positive definite (min eigenvalue > ||E||_F), so m has "
              f">= {N_TOEP} support points per period and nu is not carried by <= {N_TOEP - 1} geometric series lam_j 5^Z")
    ok18, E18, _ = toeplitz_certificate(m, N_TOEP + 1)
    print(f"   observation (no claim): T_{N_TOEP + 1} undecided at this precision (||E||_F <= {float(E18.upper()):.1e}); "
          f"certified: {ok18}")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BX checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
