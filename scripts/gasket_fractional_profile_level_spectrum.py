#!/usr/bin/env python3
r"""Theorem CG: the top of the spectrum of CE's transfer operator is certified (CF's open item (a)).  On the Hardy space
of an explicit disk, Lt has exactly three eigenvalues of modulus > 10^-3, all real and simple:
    r   = 0.42373128696260794...   +- 1e-14   (= CE's r; CE.4's 2e-38 enclosure lies inside),
    l_2 = -0.07311381751239958...  +- 1e-13,   l_2/r = -0.1725476021...  (CF.5 only gave |mu| <= 0.8855 r),
    l_3 = -0.01706764406262372...  +- 1e-12,   l_3/r = -0.0402794049...,
and every other eigenvalue has |mu| < 10^-3.  Consequently |log(s_(k+1)/(r s_k))| <= C 0.17255^k (CF.3: 0.8855^k).

Setting (GASKET.md, Theorems CE, CF): in s = sqrt(t + 25/4) the twisted operator is
    (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s),   sigma_(1,2)(s) = sqrt(15/4 -+ s),
    wt1 = (s+1/2)/(2s(7/2-s)),  wt2 = (s-1/2)/(2s(s+7/2)),   X = [sqrt5/2, 5/2] (t in [-5,0]).
Lt acts on H^2(D), D = {|s - c| < R}, c = 83/50, R = 8/5, with orthonormal basis e_n = u^n, u = (s - c)/R
(||f||^2 = sum |a_n|^2 <= sup_D |f|^2).  X lies inside D.

Theorem CG (checked here; exit 1 on any failure):
  CG.1 [proved; exact + arb] The closed disk avoids the poles s = 0, +-7/2 of the weights, and Re(15/4 -+ s) >= 49/100 > 0
       on it, so wt_i and v_i = (sigma_i - c)/R are holomorphic near the closed unit u-disk.  By the maximum principle and
       an arb cover of |u| = 1 by K_ARC arcs: |v_1| <= q_1 < 0.6001, |v_2| <= q_2 < 0.6174, |wt_1| <= W_1 < 2.4048,
       |wt_2| <= W_2 < 1.0440.  With P_N the orthogonal projection onto span(e_0..e_N):
           ||Lt (I - P_N)|| <= eps_N = sum_i W_i q_i^(N+1) / sqrt(1 - q_i^2) = 1.741e-17 (N = 80),
           ||Lt P_N|| <= beta_N = (sum_(n<=N) (W_1 q_1^n + W_2 q_2^n)^2)^(1/2) = 4.3326.
       So Lt is compact on H^2(D) (norm limit of the finite-rank Lt P_N).
  CG.2 [certified; arb power series, Pellet] M_N = P_N Lt P_N in the basis e_n (Taylor coefficients of wt_1 v_1^n +
       wt_2 v_2^n, arb series at 1100 bits).  Its characteristic polynomial chi (ball coefficients) passes Pellet's test
       |a_k| rho^k > sum_(j != k) |a_j| rho^j at rho = 3/100 (k = N-1) and rho = 1/1000 (k = N-2): M_N has exactly 2 and 3
       eigenvalues of modulus > rho.  chi changes sign at m_j -+ delta_j/2 and the disks |z - m_j| <= delta_j lie in the
       right regions, so each disk holds exactly one eigenvalue of M_N.  Lt P_N = [[M_N, 0], [B, 0]] on
       span(e_0..e_N) + its complement, so its nonzero eigenvalues are those of M_N, and for z != 0
           ||(z - Lt P_N)^-1|| <= ||(z - M_N)^-1|| (1 + beta_N/|z|) + 1/|z|.
  CG.3 [certified; arb] On each circle |z| = 3/100, |z| = 1/1000, |z - m_j| = delta_j (delta = 1e-14, 1e-13, 1e-12) the
       bound sup ||(z - Lt P_N)^-1|| eps_N < 1 holds (arb inverses of z - M_N, Frobenius norm, adaptive samples with
       step 1/(2K)).  Then z - Lt P_N - tau (Lt - Lt P_N) is invertible on the circle for tau in [0,1], the Riesz
       projections move continuously, and Lt has as many eigenvalues as M_N inside (or outside) each circle.  Hence Lt
       on H^2(D) has exactly one eigenvalue in each disk |z - m_j| < delta_j, real (Lt commutes with f -> conj f(conj s))
       and algebraically simple, and every other eigenvalue has |mu| < 10^-3.  The one near 0.4237 is CE's r: its
       eigenfunction restricted to X is a nonzero C^1 eigenfunction, so CF.5 forces mu = r; CE.4's [lo, hi] lies in the
       disk (consistency).
  CG.4 [proved from CG.3, CE.3, CF.3] For the level sums: |log(s_(k+1)/(r s_k))| <= C (|l_2|/r)^k, |l_2|/r < 0.17255
       (C not explicit).  If the l_2-component of the level-sum functional is nonzero, the rate is exactly |l_2|/r with
       alternating sign (OPEN; CE.5's certified deviations alternate and shrink by about -0.1725 per level).
Observation (no claim): tr M_N = 0.333333333333333..., matching the holomorphic fixed-point trace
sum_i wt_i(x_i)/(1 - sigma_i'(x_i)) = 1/4 + 1/12 = 1/3 (x_1 = 3/2, x_2 = 5/2); r + l_2 + l_3 - 1/3 = 2.16492e-4, the
negative of the float remainder (next eigenvalues ~ -2.1645e-4, ~ -4.5e-8).
Assumptions used as stated (not re-run): CE.1-CE.3 (s_k = 3 sigma_-3 (Lt^k G)(-3); Lt^4 G real-analytic near X), CF.3,
CF.5.  CE.4's certificate is imported and re-run unchanged.  No closed form of l_2 is claimed.  Not claimed: that l_2 is the
second eigenvalue of Lt on C^1(X) (only that it is one; analytic eigenfunctions only).  Stage 1 untouched; no W, 0.08, 0.23
or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

import sympy as sp
from flint import acb, acb_mat, arb, arb_mat, arb_series, ctx

FAILS: list[str] = []
C_CEN = Fr(83, 50)       # disk centre in s
R_RAD = Fr(8, 5)         # disk radius in s
N_TRUNC = 80             # Taylor truncation (matrix size N+1)
RADII_FAST = [Fr(3, 100)]    # pytest run: the |z| = 1/1000 circle (about 80 s) is left to `reproduce`
PREC_HI = 1100           # bits for the matrix and its characteristic polynomial (Pellet needs ~1e-316)
PREC_RES = 128           # bits for the resolvent inverses
K_ARC = 2048             # arcs covering |u| = 1
RHO = [Fr(3, 100), Fr(1, 1000)]     # Pellet / transfer radii: 2, 3 eigenvalues outside
RATIO_MAX = "0.17255"


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def _a(x: Fr) -> arb:
    return arb(x.numerator) / x.denominator


# ----------------------------------------------------------------- CG.1 the disk (exact + sympy)
def symbolic_checks():
    c, R = C_CEN, R_RAD
    check(c - R > 0 and Fr(7, 2) - c > R and Fr(7, 2) + c > R,
          "CG.1: |0 - c|, |7/2 - c|, |-7/2 - c| > R: the poles of wt_1, wt_2 lie outside the closed disk")
    check(Fr(15, 4) - c - R == Fr(49, 100) and Fr(15, 4) + c - R > 0,
          "CG.1: Re(15/4 - s) >= 15/4 - c - R = 49/100 > 0 and Re(15/4 + s) > 0 on the closed disk: principal sqrt holomorphic")
    s0sq = Fr(5, 4)
    check((c - R) ** 2 < s0sq and Fr(5, 2) < c + R, "CG.1: X = [sqrt5/2, 5/2] lies inside D = (c - R, c + R) on the real line")
    s = sp.symbols("s")
    h = sp.Rational(1, 2)
    W1 = (s + h) / (2 * s * (sp.Rational(7, 2) - s))
    W2 = (s - h) / (2 * s * (s + sp.Rational(7, 2)))
    x1 = sp.solve(sp.Eq(s ** 2, sp.Rational(15, 4) - s), s)
    x2 = sp.solve(sp.Eq(s ** 2, sp.Rational(15, 4) + s), s)
    x1 = [v for v in x1 if v > 0][0]
    x2 = [v for v in x2 if v > 0][0]
    d1 = -1 / (2 * x1)
    d2 = 1 / (2 * x2)
    tr = sp.nsimplify(W1.subs(s, x1) / (1 - d1) + W2.subs(s, x2) / (1 - d2))
    check(x1 == sp.Rational(3, 2) and x2 == sp.Rational(5, 2) and tr == sp.Rational(1, 3),
          "observation: fixed points sigma_1(3/2) = 3/2, sigma_2(5/2) = 5/2; holomorphic trace sum wt_i/(1 - sigma_i') = 1/4 + 1/12 = 1/3")
    return tr


def disk_constants(K=K_ARC):
    """Sup of |v_i| and |wt_i| on the closed unit u-disk (= on |u| = 1, maximum principle), arb arcs."""
    ctx.prec = 128
    c, R = _a(C_CEN), _a(R_RAD)
    q4, h, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    pi = arb.pi()
    q1 = q2 = W1 = W2 = arb(0)
    for j in range(K):
        th = arb(2 * j + 1) * pi / K + arb(0, (pi / K).upper())      # ball covering the arc [2 pi j/K, 2 pi (j+1)/K]
        u = acb(0, th).exp()
        s = c + R * u
        v1 = abs(((q4 - s).sqrt() - c) / R).upper()
        v2 = abs(((q4 + s).sqrt() - c) / R).upper()
        w1 = abs((s + h) / (2 * s * (s7 - s))).upper()
        w2 = abs((s - h) / (2 * s * (s + s7))).upper()
        q1, q2, W1, W2 = (arb(max(q1, v1)), arb(max(q2, v2)), arb(max(W1, w1)), arb(max(W2, w2)))
    check(bool(q1 < arb("0.6001")) and bool(q2 < arb("0.6174")),
          f"CG.1: sup_D |v_1| <= q_1 = {q1.str(6)} < 0.6001, sup_D |v_2| <= q_2 = {q2.str(6)} < 0.6174 (sigma_i(D) in D(c, q_i R))")
    check(bool(W1 < arb("2.4048")) and bool(W2 < arb("1.0440")),
          f"CG.1: sup_D |wt_1| <= W_1 = {W1.str(6)} < 2.4048, sup_D |wt_2| <= W_2 = {W2.str(6)} < 1.0440 ({K} arcs)")
    return dict(q1=q1, q2=q2, W1=W1, W2=W2)


def tail_constants(dc, N):
    ctx.prec = 128
    q1, q2, W1, W2 = dc["q1"], dc["q2"], dc["W1"], dc["W2"]
    eps = W1 * q1 ** (N + 1) / (1 - q1 ** 2).sqrt() + W2 * q2 ** (N + 1) / (1 - q2 ** 2).sqrt()
    beta2 = arb(0)
    for n in range(N + 1):
        beta2 += (W1 * q1 ** n + W2 * q2 ** n) ** 2
    beta = beta2.sqrt()
    eps, beta = arb(eps.upper()), arb(beta.upper())
    check(bool(eps < arb("1e-12")), f"CG.1: ||Lt(I - P_N)|| <= eps_N = {eps.str(4)}, ||Lt P_N|| <= beta_N = {beta.str(6)} (N = {N})")
    return eps, beta


# ----------------------------------------------------------------- CG.2 the matrix and its characteristic polynomial
def build_matrix(N, prec=PREC_HI):
    ctx.prec = prec
    ctx.cap = N + 1
    c, R = _a(C_CEN), _a(R_RAD)
    u = arb_series([0, 1], prec=N + 1)
    s = c + R * u
    q4 = arb(15) / 4
    v1 = ((q4 - s).sqrt() - c) / R
    v2 = ((q4 + s).sqrt() - c) / R
    w1 = (s + arb(1) / 2) / (2 * s * (arb(7) / 2 - s))
    w2 = (s - arb(1) / 2) / (2 * s * (s + arb(7) / 2))
    M = arb_mat(N + 1, N + 1)
    p1 = arb_series([1], prec=N + 1)
    p2 = arb_series([1], prec=N + 1)
    for n in range(N + 1):
        col = (w1 * p1 + w2 * p2).coeffs()
        col = col + [arb(0)] * (N + 1 - len(col))
        for m in range(N + 1):
            M[m, n] = col[m]
        p1 = p1 * v1
        p2 = p2 * v2
    ctx.cap = 10
    return M


def _float_eigs(M):
    import numpy as np
    A = np.array([[float(M[i, j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])
    ev = sorted(np.linalg.eigvals(A).real, key=lambda z: -abs(z))
    return ev[:3]


def _refine(P, x0, steps=60):
    """Newton on the mid polynomial (proposal only; the claim rests on sign changes)."""
    a = [arb(v.mid()) for v in P.coeffs()]
    x = arb(x0)
    for _ in range(steps):
        f = df = arb(0)
        for k in range(len(a) - 1, -1, -1):
            df = df * x + f
            f = f * x + a[k]
        x = arb((x - f / df).mid())
    return x


def charpoly_counts(M, N):
    ctx.prec = PREC_HI
    P = M.charpoly()
    a = P.coeffs()
    check(len(a) == N + 2 and a[N + 1] == 1, f"CG.2: chi_M is monic of degree {N + 1}")
    for i, rho in enumerate(RHO):
        k = N + 1 - (i + 2)
        rr = _a(rho)
        lhs = abs(a[k]).lower() * rr ** k
        rhs = arb(0)
        for j in range(len(a)):
            if j != k:
                rhs += arb(abs(a[j]).upper()) * rr ** j
        check(bool(lhs > rhs), f"CG.2: Pellet at rho = {rho}: |a_{k}| rho^{k} > sum_(j != {k}) |a_j| rho^j, so exactly "
                               f"{i + 2} eigenvalues of M_N have |z| > {rho} (none on |z| = {rho})")
    mids = []
    for x0 in _float_eigs(M):
        mids.append(_refine(P, float(x0)))
    return P, mids


def sign_changes(P, mids, deltas):
    """chi_M changes sign on [m_j - delta_j/2, m_j + delta_j/2]; the disk |z - m_j| <= delta_j lies in its Pellet region."""
    ctx.prec = PREC_HI
    big = arb(10)
    regions = [(_a(RHO[0]), big), (_a(RHO[0]), big), (_a(RHO[1]), _a(RHO[0]))]
    ok = True
    for (m, d), (a_lo, a_hi) in zip(zip(mids, deltas), regions):
        dd = _a(d)
        lo, hi = arb((m - dd / 2).mid()), arb((m + dd / 2).mid())
        sc = bool(P(lo) * P(hi) < 0)
        inside = bool(abs(m) - dd > a_lo) and bool(abs(m) + dd < a_hi)
        ok = ok and sc and inside
    ok = ok and bool(mids[0] - _a(deltas[0]) > 0) and bool(mids[1] + _a(deltas[1]) < 0)
    check(ok, "CG.2: chi_M changes sign on [m_j - delta_j/2, m_j + delta_j/2], j = 1..3, and the disks |z - m_j| <= delta_j lie "
              "in their Pellet regions (j = 1, 2: |z| > 3/100, on opposite sides of 0; j = 3: 1/1000 < |z| < 3/100). So each disk "
              "holds exactly one eigenvalue of M_N, real and simple, and these are all eigenvalues of modulus > 1/1000")
    return ok


# ----------------------------------------------------------------- CG.3 resolvent bounds on circles
def _res_norm(Mc, I, z):
    """||(z - M)^-1||_2 <= Frobenius norm = sqrt(tr(B^* B)), B an arb enclosure of the inverse."""
    B = (z * I - Mc).inv()
    F = (B.conjugate().transpose() * B).trace().real
    return arb(F.sqrt().upper())


def circle_bound(Mc, I, centre, radius, eps, beta):
    """Certify sup_Gamma ||(z - Lt P_N)^-1|| eps < 1 on Gamma = {|z - centre| = radius}.
    Adaptive cover: at the sample z_j with ||(z_j - M)^-1|| <= K_j, every z with |z - z_j| <= 1/(2 K_j) has
    ||(z - M)^-1|| <= 2 K_j; the next sample is at arc distance 1/(2 K_j) (chord <= arc), so the arcs cover Gamma."""
    ctx.prec = PREC_RES
    cen, rad = centre, _a(radius)
    pi = arb.pi()
    absmin = arb((abs(cen) - rad).lower()) if bool(abs(cen) > rad) else rad
    two_pi = 2 * pi
    ang = arb(0)
    worst = arb(0)
    n = 0
    while True:
        z = cen + rad * acb(0, ang).exp()
        K = _res_norm(Mc, I, z)
        KT = 2 * K * (1 + beta / absmin) + 1 / absmin
        worst = arb(max(worst.upper(), KT.upper()))
        n += 1
        step = arb((1 / (2 * K * rad)).lower())      # arc length 1/(2K) in angle
        ang = arb((ang + step).lower())
        if bool(ang >= two_pi):
            break
        if n > 20000:
            return False, worst, n
    return bool(worst * eps < 1), worst, n


def transfer(M, N, eps, beta, mids, deltas, radii=None):
    ctx.prec = PREC_RES
    Mc = acb_mat(M)
    I = acb_mat(N + 1, N + 1)
    for i in range(N + 1):
        I[i, i] = 1
    for rho in (RHO if radii is None else radii):
        ok, w, n = circle_bound(Mc, I, acb(0), rho, eps, beta)
        check(ok, f"CG.3: |z| = {rho}: sup ||(z - Lt P_N)^-1|| <= {w.str(4)}, times eps_N < 1 ({n} samples)")
    names = ["r", "l_2", "l_3"]
    for nm, m, d in zip(names, mids, deltas):
        ok, w, n = circle_bound(Mc, I, acb(m.mid()), d, eps, beta)
        check(ok, f"CG.3: |z - m({nm})| = {float(d):.0e}: sup ||(z - Lt P_N)^-1|| <= {w.str(4)}, times eps_N < 1 ({n} samples)")


def deltas_for(N):
    return [Fr(1, 10 ** 14), Fr(1, 10 ** 13), Fr(1, 10 ** 12)] if N >= 80 else \
           [Fr(1, 10 ** 10), Fr(1, 10 ** 9), Fr(1, 10 ** 8)]


def summary(mids, deltas, r_lo, r_hi, trM):
    ctx.prec = 200
    r_m, l2 = mids[0], mids[1]
    lo, hi = _a(r_lo), _a(r_hi)
    d0 = _a(deltas[0])
    check(bool(abs(lo - r_m) < d0) and bool(abs(hi - r_m) < d0),
          f"CG.3: CE.4's [lo, hi] lies in the certified disk around m(r) = {r_m.str(17)}: the top eigenvalue is CE's r")
    ratio_hi = (abs(l2) + _a(deltas[1])) / lo
    check(bool(ratio_hi < arb(RATIO_MAX)) and bool(abs(l2) + _a(deltas[1]) < arb("0.3752")),
          f"CG.4: |l_2|/r <= {ratio_hi.str(12)} < {RATIO_MAX} (CF.5's bound 0.8855 improved); l_2 < 0 so deviations alternate if "
          "the l_2-coefficient is nonzero")
    for nm, m, d in zip(["r", "l_2", "l_3"], mids, deltas):
        print(f"   {nm:3s} in [{(m - _a(d)).str(20)}, {(m + _a(d)).str(20)}]  (ratio to r: {(m / lo).str(12)})")
    tot = mids[0] + mids[1] + mids[2]
    print(f"   observation (no claim): tr M_N = {trM.str(15)}, holomorphic trace 1/3; r + l_2 + l_3 - 1/3 = "
          f"{(tot - arb(1) / 3).str(6)} (floats: next eigenvalue -2.1645e-4, then ~ -4.5e-8)")


def main(N=N_TRUNC) -> int:
    t0 = time.time()
    print("Theorem CG: certified top spectrum of CE's transfer operator on H^2(D): r, l_2 = -0.0731138175..., l_3 "
          "simple and real, all else |mu| < 1e-3; level-sum rate 0.17255. Not thrust.")
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import gasket_fractional_profile_level_ratio as ce
    symbolic_checks()
    dc = disk_constants()
    eps, beta = tail_constants(dc, N)
    M = build_matrix(N)
    ctx.prec = PREC_HI
    trM = sum((M[i, i] for i in range(N + 1)), arb(0))
    P, mids = charpoly_counts(M, N)
    deltas = deltas_for(N)
    sign_changes(P, mids, deltas)
    transfer(M, N, eps, beta, mids, deltas)
    r_lo, r_hi = ce.certificate()
    summary(mids, deltas, r_lo, r_hi, trM)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS or ce.FAILS:
        print(f"{len(FAILS) + len(ce.FAILS)} FAILED")
        return 1
    print("all Theorem CG checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else N_TRUNC))
