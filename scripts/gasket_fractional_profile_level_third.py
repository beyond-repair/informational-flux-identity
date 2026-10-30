#!/usr/bin/env python3
r"""Theorem CO: a_3 != 0 in CN.2.  The l_3-component of the level sums is certified nonzero, a_3 = -2.36116e-5 +- 1e-9 (approx.),
so the level sums have the exact second-order rate |l_3|/r.

Setting (GASKET.md, Theorems CE-CG, CN): in s = sqrt(t + 25/4), (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s) on
H^2(D), D = {|s - c| < R}, c = 83/50, R = 8/5, orthonormal basis e_n = u^n, u = (s - c)/R; M_N = P_N Lt P_N (N = 80),
||Lt (I - P_N)|| <= eps_N = 1.741e-17 (CG.1).  F0 = Lt^4 G (CN.1), x0 = sqrt13/2, c0 = 3 sigma_-3, and CN.2:
    s_(k+4) = a_1 r^k + a_2 l_2^k + a_3 l_3^k + e_k,   a_j = c0 (Pi_j F0)(x0),   |e_k| <= eps_0 1e-3^k.
CN.3 decided a_1 > 0 and a_2 != 0 from the level sums; a_3 was left open (|Q_7| = 2.46e-19 against a 9e-19 remainder).
Here a_j is computed from the projection itself, not from the level sums.

Theorem CO (checked here; exit 1 on any failure):
  CO.1 [certified; exact + acb] On the larger disk D_rho = {|s - c| < rho R}, rho = 51/50: c - rho R > 0 and
       Re(15/4 - s) > 0 (exact), sup_(|u| = rho) |v_i| < 1 (arcs), so the weights and branches are holomorphic on the closed
       D_rho and sigma_i(D_rho) lies in D; Krawczyk boxes (CN.1's test, CN.1's box V) on every arc of |s - c| = rho R put
       t(Phi s) in psi(V°), so F0 is holomorphic near the closed D_rho with sup_(|u| = rho) |F0| <= S_rho.  Cauchy:
       |f_n| <= S_rho rho^-n for the Taylor coefficients f_n of F0 in u.
  CO.2 [certified; acb] The M-point trapezoid sums on |u| = 1 (certified point values of F0) give f_0..f_N with
       |f~_n - f_n| <= S_rho rho^-(n+M) / (1 - rho^-M) (aliasing), M = 2048.
  CO.3 [certified; arb Krawczyk] For mu in {r, l_2, l_3}: a right and a left eigenvector v, w of M_N for its eigenvalue
       mu_N in CG.2's disk are enclosed (Krawczyk on (A - mu I) v = 0, v_k = 1).  The Riesz projection of Lt P_N at mu_N is
       Pi^N h = (w . P_N h)/(w . v) phi, phi = (Lt v)/mu_N, so a_j^N = c0 (w . f)/(w . v) (Lt v)(x0)/mu_N exactly.
  CO.4 [proved] (z - Lt)^-1 - (z - Lt P_N)^-1 = (z - Lt)^-1 Lt (I - P_N) (z - Lt P_N)^-1, so on CN's circles
       ||Pi_j - Pi_j^N|| <= rad K_j K_j^N eps_N (K_j, K_j^N certified resolvent bounds) and
           |a_j - a_j^N| <= c0 ||ev_x0|| rad K_j K_j^N eps_N S,   ||ev_x0|| = 1/sqrt(1 - u0^2).
       Result: a_3 != 0 (negative); a_1, a_2 agree with CN.3's independent level-sum enclosures (consistency).
  CO.5 [proved from CO.4, CN.2, CG.3] With q_j = l_j/r, B = (a_2/a_1) q_2^-4 (CN.3), B_3 = (a_3/a_1) q_3^-4 != 0 and
       q_2^2 < |q_3| (certified):
           (log(s_(n+1)/(r s_n)) - B q_2^n (q_2 - 1)) / q_3^n  ->  B_3 (q_3 - 1) != 0,
       so the correction to CN.4's leading term has the exact rate |l_3|/r (not (l_2/r)^2, not 1e-3/r).
Kept limits: no explicit constant n_0 from which the second-order sign pattern holds; a_4 (the next eigenvalue, about
-2.16e-4 inside |z| < 1e-3) is not touched; C_G, circle eigenvalues |z| = rho(L_k), closed forms stay open.
Assumptions used as stated (not re-run): CE.1-CE.3, CG.2-CG.3 (enclosures and count of r, l_2, l_3, simple), CN.1's
box checks (re-run in the main run) and CN.2's circles (re-run).  Textbook facts: Cauchy estimates, trapezoid aliasing
for power series, Krawczyk/Brouwer, Riesz projections, the second resolvent identity.  Stage 1 untouched; no W, 0.08, 0.23
or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
from flint import acb, acb_mat, arb, arb_mat, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_coefficient as cn     # noqa: E402
import gasket_fractional_profile_level_spectrum as cg        # noqa: E402

FAILS: list[str] = []
PREC = 200
RHO_BIG = Fr(51, 50)
M_PTS = 2048
N = cn.N_TRUNC
MIDS = {"r": cn.R_MID, "l2": cn.L2_MID, "l3": cn.L3_MID}
# CN.3's level-sum enclosures (consistency only; not used in the proof of CO)
CN3 = {"r": ("0.00571256871089", "1.2e-15"), "l2": ("-4.158e-5", "4.4e-9")}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


_a = cn._a


# ----------------------------------------------------------------- CO.1: the larger disk
def big_disk(rho=RHO_BIG, K=4096):
    c, R = cn.C_CEN, cn.R_RAD
    ok_exact = c - rho * R > 0 and Fr(15, 4) - (c + rho * R) > 0 and c + rho * R < Fr(7, 2) and rho > 1
    ctx.prec = 128
    ca, Ra, rh = _a(c), _a(R), _a(rho)
    q4 = arb(15) / 4
    pi = arb.pi()
    q1 = q2 = arb(0)
    for j in range(K):
        th = arb(2 * j + 1) * pi / K + arb(0, (pi / K).upper())
        s = ca + Ra * rh * acb(0, th).exp()
        q1 = arb(max(q1, abs(((q4 - s).sqrt() - ca) / Ra).upper()))
        q2 = arb(max(q2, abs(((q4 + s).sqrt() - ca) / Ra).upper()))
    check(ok_exact and bool(q1 < 1) and bool(q2 < 1),
          f"CO.1 (exact + arb, {K} arcs): rho = {rho}: c - rho R = {float(c - rho * R):.4f} > 0, 15/4 - c - rho R = "
          f"{float(Fr(15, 4) - c - rho * R):.4f} > 0, c + rho R < 7/2 (weights, branches holomorphic on the closed D_rho); "
          f"sup_(|u|=rho) |v_1| <= {q1.str(5)}, |v_2| <= {q2.str(5)} < 1, so sigma_i(D_rho) lies in D (CG.1 from there on)")
    return ok_exact


def circle_sup_rho(rho=RHO_BIG, K=4096):
    """CN.1's Krawczyk arc test on |s - c| = rho R."""
    ctx.prec = PREC
    mp.mp.dps = 30
    gam = arb(3).log() / arb(5).log()
    c, R = _a(cn.C_CEN), _a(cn.R_RAD) * _a(rho)
    pi = arb.pi()
    guesses = _guesses(acb(c + R))
    S = arb(0)
    ok = True
    for j in range(K):
        th = arb(2 * j + 1) * pi / K + arb(0, (pi / K).upper())
        y = c + R * acb(0, th).exp()
        val, guesses = cn.F0_bound(y, guesses, gam)
        if val is None:
            ok = False
            print(f"   arc {j} failed")
            break
        S = arb(max(S, arb(abs(val).upper())).upper())
    check(ok and bool(S.is_finite()), f"CO.1 (acb, {K} arcs x 5 words): Krawczyk boxes in int V on |s - c| = {rho} R, so F0 is "
                                      f"holomorphic near the closed D_rho and sup_(|u| = rho) |F0| <= S_rho = {S.str(6)}")
    return S


def _guesses(y0):
    out = []
    for w in cn.INNER:
        _, s = cn.word_eval(y0, w)
        t = float((s * s - arb(25) / 4).real.mid())
        u0 = mp.findroot(lambda u: cn._psi_mp(u) - t, (-5, -1), solver="anderson")
        out.append(acb(arb(mp.nstr(u0, 30))))
    return out


# ----------------------------------------------------------------- CO.2: Taylor coefficients of F0
def taylor_coeffs(S_rho, rho=RHO_BIG, M=M_PTS, Nc=N):
    ctx.prec = PREC
    mp.mp.dps = 30
    gam = arb(3).log() / arb(5).log()
    c, R = _a(cn.C_CEN), _a(cn.R_RAD)
    pi = arb.pi()
    roots = [acb(0, 2 * pi * j / M).exp() for j in range(M)]
    vals = []
    guesses = _guesses(acb(c + R))
    for j in range(M):
        val, guesses = cn.F0_bound(c + R * roots[j], guesses, gam, rad=arb("1e-30"))
        if val is None:
            check(False, f"CO.2: point evaluation {j} failed")
            return None
        vals.append(val)
    rh = _a(rho)
    tail = S_rho * rh ** (-M) / (1 - rh ** (-M))
    f = []
    for n in range(Nc + 1):
        acc = acb(0)
        for j in range(M):
            acc += vals[j] * roots[(-j * n) % M]
        fn = (acc / M).real
        f.append(fn + arb(0, (tail * rh ** (-n)).upper()))
    width = max(float(x.rad()) for x in f)
    check(bool(tail < arb("1e-15")) and width < 1e-15,
          f"CO.2 (acb, M = {M}): f_0..f_{Nc} enclosed, aliasing <= S_rho rho^-M/(1 - rho^-M) = {tail.str(3)}, max radius "
          f"{width:.2e}; f_0 = {f[0].str(10)}, f_{Nc} = {f[Nc].str(6)}")
    return f


# ----------------------------------------------------------------- CO.3: eigenvectors of M_N
def eigvec(A, lam_mid, n):
    """Krawczyk enclosure of (v, mu) with (A - mu I) v = 0, v_k = 1, near the float eigenpair at lam_mid."""
    import numpy as np
    ctx.prec = 320
    Af = np.array([[float(A[i, j].mid()) for j in range(n)] for i in range(n)])
    ev, V = np.linalg.eig(Af)
    i = int(np.argmin(abs(ev - float(lam_mid))))
    v0 = V[:, i].real
    k = int(np.argmax(abs(v0)))
    v0 = v0 / v0[k]
    x = [arb(float(t)) for t in v0]
    lam = arb(float(ev[i].real))

    def F(xv, lm):
        out = []
        for r_ in range(n):
            acc = -lm * xv[r_]
            for j in range(n):
                acc += A[r_, j] * xv[j]
            out.append(acc)
        return out

    def J(xv, lm):
        Jm = arb_mat(n, n)
        for r_ in range(n):
            for j in range(n):
                if j == k:
                    Jm[r_, j] = -xv[r_]
                else:
                    Jm[r_, j] = A[r_, j] - (lm if r_ == j else 0)
        return Jm
    # Newton on mid values (unknowns: v_j for j != k, and mu in slot k)
    for _ in range(8):
        Fx = F(x, lam)
        d = J(x, lam).solve(arb_mat([[t] for t in Fx]))
        x = [arb((x[j] - d[j, 0]).mid()) if j != k else arb(1) for j in range(n)]
        lam = arb((lam - d[k, 0]).mid())
    Jm = J(x, lam)
    Jmid = arb_mat([[arb(Jm[i_, j].mid()) for j in range(n)] for i_ in range(n)])
    Y = Jmid.inv()
    Y = arb_mat([[arb(Y[i_, j].mid()) for j in range(n)] for i_ in range(n)])
    Fx = arb_mat([[t] for t in F(x, lam)])
    YF = Y * Fx
    for e in (arb("1e-40"), arb("1e-30"), arb("1e-24"), arb("1e-18")):
        X = [x[j] + arb(0, e) if j != k else arb(1) for j in range(n)]
        Lm = lam + arb(0, e)
        JX = J(X, Lm)
        Dm = arb_mat(n, n)
        for i_ in range(n):
            Dm[i_, i_] = 1
        Dm = Dm - Y * JX
        dx = arb_mat([[arb(0, e)] for _ in range(n)])
        Kx = Dm * dx
        ok = True
        res = []
        for j in range(n):
            kj = -YF[j, 0] + Kx[j, 0]          # K(X) - x~ in slot j
            ok = ok and bool(abs(kj) < e)
            res.append(kj)
        if ok:
            vv = [x[j] + arb(0, e) if j != k else arb(1) for j in range(n)]
            return vv, lam + arb(0, e), e
    return None, None, None


def lt_eval(v, y):
    """(Lt v)(y) for the polynomial v(s) = sum v_n ((s - c)/R)^n."""
    c, R = _a(cn.C_CEN), _a(cn.R_RAD)
    q4, h, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    w1 = (y + h) / (2 * y * (s7 - y))
    w2 = (y - h) / (2 * y * (y + s7))

    def poly(s):
        u = (s - c) / R
        acc = arb(0)
        for t in reversed(v):
            acc = acc * u + t
        return acc
    return w1 * poly((q4 - y).sqrt()) + w2 * poly((q4 + y).sqrt())


def projection_coeffs(f, S):
    ctx.prec = 320
    M = cg.build_matrix(N, prec=320)
    MT = M.transpose()
    import gasket_fractional_profile_level_ratio as ce
    c0 = 3 * ce.sigma_closed(arb(-3), 0)
    x0 = arb(13).sqrt() / 2
    u0 = (x0 - _a(cn.C_CEN)) / _a(cn.R_RAD)
    ev = 1 / (1 - u0 * u0).sqrt()
    dc = cg.disk_constants()
    eps, beta = cg.tail_constants(dc, N)
    Mc = acb_mat(M)
    Id = acb_mat(N + 1, N + 1)
    for i in range(N + 1):
        Id[i, i] = 1
    eig_disk = {"r": (cn.R_MID, "1e-14"), "l2": (cn.L2_MID, cn.L2_RAD), "l3": (cn.L3_MID, cn.L3_RAD)}
    out = {}
    for name in ("r", "l2", "l3"):
        ctx.prec = 320
        v, muN, e1 = eigvec(M, MIDS[name], N + 1)
        w, muL, e2 = eigvec(MT, MIDS[name], N + 1)
        ok_v = v is not None and w is not None
        if not ok_v:
            check(False, f"CO.3: Krawczyk eigenvector enclosure for {name} failed")
            continue
        m, d = eig_disk[name]
        dsk = arb(m) + arb(0, arb(d))
        in_disk = bool(abs(muN - arb(m)) < arb(d)) and bool(abs(muL - arb(m)) < arb(d))
        check(in_disk, f"CO.3 (arb Krawczyk, 320 bits): right/left eigenvectors of M_{N} for mu_N = {muN.str(15)} enclosed "
                       f"(radius {float(e1):.0e}/{float(e2):.0e}); mu_N lies in CG.2's disk {dsk.str(6)} (its only eigenvalue there)")
        wv = sum((w[i] * v[i] for i in range(N + 1)), arb(0))
        wf = sum((w[i] * f[i] for i in range(N + 1)), arb(0))
        phi = lt_eval(v, x0) / muN
        aN = c0 * wf / wv * phi
        ctx.prec = cg.PREC_RES
        cen, rad = cn.CIRCLES[name]
        okc, KN, nsamp = cg.circle_bound(Mc, Id, acb(arb(cen)), rad, eps, beta)
        KT = KN / (1 - KN * eps)
        ctx.prec = 320
        err = c0 * ev * _a(rad) * KT * KN * eps * S
        a = aN + arb(0, err.upper())
        check(okc, f"CO.4: circle |z - {cen[:10]}| = {rad}: K^N <= {KN.str(5)}, K <= {KT.str(5)} ({nsamp} samples); "
                   f"|a_{name} - a_{name}^N| <= c0 ||ev|| rad K K^N eps_N S = {err.str(3)}")
        print(f"   a_{name}^N = {aN.str(15)};  a_{name} in {a.str(10)}")
        out[name] = a
    return out, c0


def compare_and_conclude(a):
    ctx.prec = 320
    for name, (m, r_) in CN3.items():
        enc = arb(m) + arb(0, arb(r_))
        check(bool(a[name].overlaps(enc)), f"CO.4 (consistency): a_{name} = {a[name].str(12)} overlaps CN.3's level-sum "
                                           f"enclosure {enc.str(10)} (independent method)")
    a3 = a["l3"]
    check(not a3.contains(0) and bool(a3 < 0), f"CO.4: a_3 = {a3.str(8)} != 0 (negative): the l_3-component of the level "
                                               "sums is not zero (CN.3 left this open)")
    check(bool(abs(a3) < arb("1.125e-4")), "CO.4 (consistency): |a_3| < CN.3's data bound 1.125e-4")
    r = arb(cn.R_MID) + arb(0, arb("1e-35"))
    l2 = arb(cn.L2_MID) + arb(0, arb(cn.L2_RAD))
    l3 = arb(cn.L3_MID) + arb(0, arb(cn.L3_RAD))
    q2, q3 = l2 / r, l3 / r
    B = a["l2"] / a["r"] / q2 ** 4
    B3 = a3 / a["r"] / q3 ** 4
    lim = B3 * (q3 - 1)
    check(bool(q2 * q2 < abs(q3)) and bool(arb("1e-3") / r < abs(q3)) and not B3.contains(0),
          f"CO.5: q_2^2 = {(q2 * q2).str(8)} < |q_3| = {abs(q3).str(8)}, 1e-3/r = {(arb('1e-3') / r).str(6)} < |q_3|; "
          f"B = {B.str(8)}, B_3 = (a_3/a_1) q_3^-4 = {B3.str(8)} != 0, so (log(s_(n+1)/(r s_n)) - B q_2^n (q_2 - 1)) / q_3^n -> "
          f"B_3 (q_3 - 1) = {lim.str(8)}: exact second-order rate |l_3|/r = {abs(q3).str(10)}")
    return B, B3


def main(arcs=4096, M=M_PTS) -> int:
    t0 = time.time()
    print("Theorem CO: a_3 != 0 in CN.2 (the l_3-component of the level sums), and the exact second-order rate |l_3|/r. "
          "Not thrust.")
    cn.tail_constants()
    cn.box_checks()
    S = cn.circle_sup(1024)
    big_disk()
    S_rho = circle_sup_rho(K=arcs)
    f = taylor_coeffs(S_rho, M=M)
    if f is None:
        return 1
    a, c0 = projection_coeffs(f, S)
    if len(a) == 3:
        compare_and_conclude(a)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS or cn.FAILS or cg.FAILS:
        print(f"{len(FAILS) + len(cn.FAILS) + len(cg.FAILS)} FAILED")
        return 1
    print("all Theorem CO checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
