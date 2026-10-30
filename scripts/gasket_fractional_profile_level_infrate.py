#!/usr/bin/env python3
r"""Theorem CT: the exact rate of inf R_n -> rho(L_1) (CJ.2's open item), with explicit two-sided constants.

Setting (GASKET.md, Theorems CI, CJ): in s = sqrt(t + 25/4), X = [sqrt5/2, 5/2], sigma_(1,2)(s) = sqrt(15/4 -+ s),
    (L_1 f)(s) = u_1(s) f(sigma_1 s) + u_2(s) f(sigma_2 s),   u_i = wt_i |sigma_i'| = wt_i / (2 sigma_i),
    wt_1 = (s+1/2)/(2s(7/2-s)),  wt_2 = (s-1/2)/(2s(s+7/2)),
h_n = L_1^n 1, R_n = h_(n+1)/h_n on X.  CJ.2 proved inf_X R_n <= rho(L_1) <= sup_X R_n and
log(sup R_n / inf R_n) <= tanh(Delta/4)^(n-1) Delta (Delta < 13.31, tanh(Delta/4) ~ 0.9974); no rate for inf R_n was
claimed (CJ's Observation).  Here L_1 acts on H^2(D), D = {|s - c| < R}, c = 83/50, R = 8/5 (CG's disk), e_n = u^n.

Theorem CT (checked here; exit 1 on any failure).  N = 120, 512-bit arb.
  CT.1 [proved; exact + arb] CG.1's disk facts (poles of wt_i outside, Re(15/4 -+ s) >= 49/100) also give |sigma_i| >= 7/10
       on the closed disk, so u_i is holomorphic there; arb arcs on |u| = 1 give sup|v_i| <= q_i (CG.1's), sup|u_i| <= U_i,
       hence ||L_1 (I - P_N)|| <= eps_N, ||L_1 P_N|| <= beta_N (CG.1's formula with U_i in place of W_i): L_1 is compact.
  CT.2 [certified; CR.1/CR.2's method] arb Krawczyk right/left eigenpairs of M_N = P_N L_1 P_N at rho and nu_2; certified
       Neumann bound on A = M Q puts every other eigenvalue of M_N in |z| < 1/100 and bounds the resolvent on circles.
  CT.3 [proved; CG.3's homotopy] K_T eps_N < 1 on |z| = 1/100, |z| = R_out = 6, and on circles of radius 1e-20 about rho and
       nu_2: L_1 on H^2(D) has exactly two eigenvalues in |z| > 1/100, both real and algebraically simple,
           rho_H = 0.1343549829703998131 +- 2.9e-20,     nu_2 = -0.0327933736640408162 +- 4.0e-20,
       q = nu_2/rho_H = -0.244080070117426 +- 4e-16; every other eigenvalue has |mu| < 1/100.  rho_H = rho(L_1) on C(X):
       <= because sup_X |f| <= e_X ||f||_H2 and 1 in H^2(D) (Gelfand: rho(L_1) = lim ||L_1^n 1||^(1/n)), >= because the
       rho_H-eigenfunction is analytic and not identically 0 on X, hence a C(X) eigenfunction.  This sharpens CI.1's
       [0.1343549, 0.1343551] to +- 2.9e-20 (CI.1's interval contains it: consistency).
  CT.4 [proved; second resolvent identity] With the Riesz projections P_rho, P_nu of L_1, A = P_rho 1, B = P_nu 1:
           h_n = rho^n A + nu_2^n B + E_n  (n >= 0),   sup_X |E_n| <= e_X K_in (1/100)^(n+1),  K_in <= 8.298e5,
       e_X = sup (1 - |u|^2)^(-1/2) <= 1.17495.  A(x), B(x) are enclosed pointwise by (L_1 p)(x)/mu_N, p = the M_N-Riesz
       image of e_0, with error e_X rad K K_T eps_N (<= 6.8e-20 and 5.4e-18; projection circles of radius 1/20 and 1/100).
       An arb branch-and-bound cover of X gives min_X A >= 0.833891, sup_X |B| <= 0.984690, and
           Gamma_+ = max_X B/A in [0.2000641183, 0.2000642055],   Gamma_- = max_X (-B/A) in [0.5249039315, 0.5249041376],
       both > 0: B/A changes sign on X.
  CT.5 [proved from CT.3, CT.4] For x in X and n >= 0, with G = B/A, eta_n = E_n/rho^n, |eta_n| <= tau_1 (rho_in/rho)^n,
       tau_1 = e_X K_in/100 <= 9750.0,
           (rho - R_n(x)) / (rho (1-q) |q|^n) = [(-1)^n G(x) + (eta_n - eta_(n+1))/((1-q)|q|^n A)] / [1 + (q^n B + eta_n)/A],
       so it lies within theta_n = (max|G| eps_n + Delta_n)/(1 - eps_n) of (-1)^n G(x), where eps_n = (|q|^n sup|B| + tau_n)/min A
       and Delta_n = (tau_n + tau_(n+1))/((1-q)|q|^n min A); theta_n decreases, theta_10 <= 0.0702049 < Gamma_+.  Hence for
       every n >= n_0 = 10:
           rho (1-q) |q|^n (Gamma_par(n) - theta_n) <= rho(L_1) - inf_X R_n <= rho (1-q) |q|^n (Gamma_par(n) + theta_n)
       (par = + for n even, - for n odd), explicitly 0.0217058 |q|^n <= rho - inf R_n <= 0.0451750 |q|^n (n even) and
       0.0841585 |q|^n <= rho - inf R_n <= 0.0913152 |q|^n (n odd), and the same for sup_X R_n - rho with the opposite
       parity.  So inf R_n -> rho(L_1) at the EXACT rate |nu_2|/rho(L_1) = 0.2440800701 (CJ.2's Birkhoff rate:
       tanh(Delta/4) ~ 0.9974), with (rho - inf R_(2m))/|q|^(2m) -> c_+ = rho(1-q)Gamma_+ = 0.0334404 +- 1.2e-8 and
       (rho - inf R_(2m+1))/|q|^(2m+1) -> c_- = rho(1-q)Gamma_- = 0.0877368 +- 6.4e-8; the oscillation
       sup R_n - inf R_n lies within rho(1-q)|q|^n (Gamma_+ + Gamma_- -+ 2 theta_n), limit constant 0.1211772 +- 6.7e-8 for both
       parities, so log(sup R_n / inf R_n) also has the exact rate |q|.
Observation (no claim): floats of M_N^n e_0 on a 2001-point grid give (rho - min R_n)/(rho(1-q)|q|^n) = 0.2000640, 0.524903
for n = 8, 9, matching Gamma_+, Gamma_-; the extrema of B/A sit at the endpoints s = sqrt5/2 (Gamma_+) and s = 5/2 (Gamma_-).
Kept limits: n < 10 is not covered (theta_n >= Gamma_+ there with these constants); Gamma_+- are known to ~1e-7, not
closed form; nothing is claimed about L_k for k >= 2, the circle |z| = rho(L_k), or the level sums.
Assumptions used as stated (not re-run): CJ.2's definition of R_n, CI.1's interval (consistency only), CG.1's disk and
tail formula (constants re-run for the new weights), CG.3's homotopy argument, CR.1/CR.2's Krawczyk and Neumann routines
(imported unchanged).  Textbook: Riesz projections, Neumann series, second resolvent identity, Gelfand's formula,
point evaluation in H^2.  Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph,
not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction as Fr

from flint import acb, arb, arb_mat, arb_series, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_fifth as cr          # noqa: E402

FAILS: list[str] = []
RESULT: dict = {}
N = int(os.environ.get("CT_N", "120"))
PREC = 512
C_CEN, R_RAD = Fr(83, 50), Fr(8, 5)
RHO_IN = Fr(1, 100)
R_OUT = Fr(6)
SEEDS = ["0.134354982970399880", "-0.032793373664040810"]
NAMES = ["rho", "nu2"]
G_RAD = [Fr(1, 10 ** 20), Fr(1, 10 ** 20)]
PROJ_RAD = [Fr(1, 20), Fr(1, 100)]
CI1 = ("0.1343549", "0.1343551")
K_ARC = 2048
COVER0 = 256          # initial subintervals of X
WMIN = Fr(1, 10 ** 7)  # branch-and-bound stops at this width
NMAX = 60


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _a(x: Fr) -> arb:
    return arb(x.numerator) / x.denominator


# ----------------------------------------------------------------- CT.1 the disk
def disk_constants(K=K_ARC):
    ctx.prec = 128
    c, R = _a(C_CEN), _a(R_RAD)
    check(C_CEN - R_RAD > 0 and Fr(7, 2) - C_CEN > R_RAD and Fr(15, 4) - C_CEN - R_RAD == Fr(49, 100),
          "CT.1: poles 0, 7/2 of wt_i outside the closed disk; Re(15/4 - s) >= 49/100, so |sigma_i| >= 7/10 and u_i holomorphic")
    q4, h, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    pi = arb.pi()
    q1 = q2 = U1 = U2 = arb(0)
    for j in range(K):
        th = arb(2 * j + 1) * pi / K + arb(0, (pi / K).upper())
        s = c + R * acb(0, th).exp()
        g1, g2 = (q4 - s).sqrt(), (q4 + s).sqrt()
        v1 = abs((g1 - c) / R).upper()
        v2 = abs((g2 - c) / R).upper()
        w1 = abs((s + h) / (2 * s * (s7 - s)) / (2 * g1)).upper()
        w2 = abs((s - h) / (2 * s * (s + s7)) / (2 * g2)).upper()
        q1, q2, U1, U2 = (arb(max(q1, v1)), arb(max(q2, v2)), arb(max(U1, w1)), arb(max(U2, w2)))
    check(bool(q1 < 1) and bool(q2 < 1) and U1.is_finite() and U2.is_finite(),
          f"CT.1: sup_D |v_1| <= {q1.str(5)}, sup_D |v_2| <= {q2.str(5)}, sup_D |u_1| <= U_1 = {U1.str(5)}, "
          f"sup_D |u_2| <= U_2 = {U2.str(5)} ({K} arcs)")
    return q1, q2, U1, U2


def tail_constants(dc, n):
    ctx.prec = 128
    q1, q2, U1, U2 = dc
    eps = U1 * q1 ** (n + 1) / (1 - q1 ** 2).sqrt() + U2 * q2 ** (n + 1) / (1 - q2 ** 2).sqrt()
    b2 = arb(0)
    for k in range(n + 1):
        b2 += (U1 * q1 ** k + U2 * q2 ** k) ** 2
    eps, beta = arb(eps.upper()), arb(b2.sqrt().upper())
    check(bool(eps < arb("1e-20")), f"CT.1: ||L_1 (I - P_N)|| <= eps_N = {eps.str(4)}, ||L_1 P_N|| <= beta_N = {beta.str(5)} (N = {n})")
    return eps, beta


def build_matrix(n, prec=PREC):
    ctx.prec = prec
    ctx.cap = n + 1
    c, R = _a(C_CEN), _a(R_RAD)
    u = arb_series([0, 1], prec=n + 1)
    s = c + R * u
    q4 = arb(15) / 4
    g1, g2 = (q4 - s).sqrt(), (q4 + s).sqrt()
    v1, v2 = (g1 - c) / R, (g2 - c) / R
    w1 = (s + arb(1) / 2) / (2 * s * (arb(7) / 2 - s)) / (2 * g1)
    w2 = (s - arb(1) / 2) / (2 * s * (s + arb(7) / 2)) / (2 * g2)
    M = arb_mat(n + 1, n + 1)
    p1 = arb_series([1], prec=n + 1)
    p2 = arb_series([1], prec=n + 1)
    for k in range(n + 1):
        col = (w1 * p1 + w2 * p2).coeffs()
        col = col + [arb(0)] * (n + 1 - len(col))
        for m in range(n + 1):
            M[m, k] = col[m]
        p1 = p1 * v1
        p2 = p2 * v2
    ctx.cap = 10
    return M


def l1_eval(p, y):
    """(L_1 p)(y) for the polynomial p(s) = sum p_k ((s - c)/R)^k; y an arb (point or ball in X)."""
    c, R = _a(C_CEN), _a(R_RAD)
    q4, h, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    g1, g2 = (q4 - y).sqrt(), (q4 + y).sqrt()
    w1 = (y + h) / (2 * y * (s7 - y)) / (2 * g1)
    w2 = (y - h) / (2 * y * (y + s7)) / (2 * g2)

    def poly(t):
        uu = (t - c) / R
        acc = arb(0)
        for a in reversed(p):
            acc = acc * uu + a
        return acc
    return w1 * poly(g1) + w2 * poly(g2)


# ----------------------------------------------------------------- CT.2 / CT.3
def spectral_split(M):
    n = N + 1
    ctx.prec = PREC
    cr.PREC = PREC
    MT = M.transpose()
    pairs = []
    for name, sd in zip(NAMES, SEEDS):
        v, mu, e1 = cr.eigpair(M, sd, n)
        w, muL, e2 = cr.eigpair(MT, sd, n)
        if v is None or w is None or not bool(mu.overlaps(muL)):
            check(False, f"CT.2: Krawczyk enclosure for {name} failed")
            return None
        wv = (w.transpose() * v)[0, 0]
        Pn = arb((cr._fro(v) * cr._fro(w) / abs(wv)).upper())
        pairs.append(dict(name=name, v=v, w=w, mu=mu, wv=wv, P=Pn))
        check(not wv.contains(0), f"CT.2 (arb Krawczyk, {PREC} bits, N = {N}): {name}: mu_N = {mu.str(22)} (radii "
                                  f"{float(e1):.0e}/{float(e2):.0e}), w.v != 0, ||P|| <= {Pn.str(5)}")
    check(not pairs[0]["mu"].overlaps(pairs[1]["mu"]), "CT.2: the two enclosures are real and disjoint")
    Id = arb_mat(n, n)
    for i in range(n):
        Id[i, i] = 1
    Q = Id
    for p in pairs:
        Q = Q - p["v"] * p["w"].transpose() * (1 / p["wv"])
    A = M * Q
    T = []
    X = Q
    for _ in range(cr.M_POW):
        T.append(cr._fro(X))
        X = A * X
    alpha = cr._fro(X)
    check(all(t.is_finite() for t in T) and bool(alpha / _a(RHO_IN) ** cr.M_POW < 1),
          f"CT.2 (arb): ||A^j Q||_F <= " + ", ".join(t.str(3) for t in T) + f", ||A^{cr.M_POW}||_F <= {alpha.str(3)} "
          f"< (1/100)^{cr.M_POW}: every other eigenvalue of M_N has |z| < 1/100")
    return pairs, T, alpha


def certify(pairs, T, alpha, eps, beta):
    out = {}
    circles = [("in", arb(0), _a(RHO_IN)), ("out", arb(0), _a(R_OUT))]
    for p, rad in zip(pairs, G_RAD):
        circles.append((p["name"], arb(p["mu"].mid()), _a(rad)))
    for p, rad in zip(pairs, PROJ_RAD):
        circles.append((p["name"] + "proj", arb(p["mu"].mid()), _a(rad)))
    for name, cen, rad in circles:
        KM, absmin = cr.circle_K(pairs, T, alpha, cen, rad)
        if KM is None:
            check(False, f"CT.3: circle {name} not certified")
            out[name] = None
            continue
        KT = KM * (1 + beta / absmin) + 1 / absmin
        ok = bool(KT * eps < 1)
        KL = arb((KT / (1 - KT * eps)).upper()) if ok else None
        check(ok, f"CT.3: |z - {cen.str(12)}| = {rad.str(3)}: sup ||(z - M)^-1|| <= {KM.str(4)}, K_T = {KT.str(4)}, "
                  f"K_T eps_N = {(KT * eps).str(3)} < 1" + (f", sup ||(z - L_1)^-1|| <= {KL.str(4)}" if ok else ""))
        out[name] = (KM, KT, KL, cen, rad) if ok else None
    rin = _a(RHO_IN)
    lay = bool(beta + eps < _a(R_OUT))
    for p, rad in zip(pairs, G_RAD):
        lay = lay and bool(abs(p["mu"] - arb(p["mu"].mid())) < _a(rad)) and bool(abs(p["mu"]) - _a(rad) > rin)
    lay = lay and bool(abs(pairs[0]["mu"] - pairs[1]["mu"]) > _a(PROJ_RAD[0]) + _a(PROJ_RAD[1]))
    for p, rad in zip(pairs, PROJ_RAD):
        lay = lay and bool(abs(p["mu"]) - _a(rad) > rin)
    allok = lay and all(out[k] is not None for k in out)
    rho = arb(pairs[0]["mu"].mid()) + arb(0, (_a(G_RAD[0]) + pairs[0]["mu"].rad()).upper())
    nu = arb(pairs[1]["mu"].mid()) + arb(0, (_a(G_RAD[1]) + pairs[1]["mu"].rad()).upper())
    check(allok, f"CT.3: beta_N + eps_N < {R_OUT}; the small and projection circles hold mu_i only and stay in |z| > 1/100. "
                 f"CG.3's homotopy: L_1 on H^2(D) has exactly two eigenvalues in |z| > 1/100, real and simple: "
                 f"rho_H = {rho.str(20)}, nu_2 = {nu.str(20)} (+- 1e-20); every other eigenvalue has |mu| < 1/100")
    lo, hi = arb(CI1[0]), arb(CI1[1])
    check(bool(rho > lo) and bool(rho < hi), f"CT.3: rho_H lies in CI.1's [{CI1[0]}, {CI1[1]}] (consistency); rho_H = rho(L_1) "
                                             "on C(X) (point evaluation <= e_X ||.||_H2 and the analytic eigenfunction is not 0 on X)")
    q = nu / rho
    check(bool(q < 0) and bool(abs(q) < arb("0.2441")), f"CT.3: q = nu_2/rho(L_1) = {q.str(15)}, |q| < 0.2441")
    return out, rho, nu, q, allok


# ----------------------------------------------------------------- CT.4 pointwise A, B and the cover
def ev_X():
    """e_X >= sup over the cover of X of (1 - |u|^2)^(-1/2): |u| <= max((c - sqrt5/2), (5/2 - c))/R."""
    ctx.prec = 128
    c, R = _a(C_CEN), _a(R_RAD)
    um = arb(max((c - arb(5).sqrt() / 2).upper(), (arb(5) / 2 - c).upper())) / R + arb("1e-6")
    return arb((1 / (1 - um * um).sqrt()).upper())


def riesz_polys(pairs):
    ctx.prec = PREC
    polys = []
    for p in pairs:
        v = cr._col(p["v"])
        w0 = p["w"][0, 0]
        polys.append(([vi * w0 / p["wv"] for vi in v], p["mu"]))
    return polys


def AB(polys, y, errA, errB):
    (pa, ma), (pb, mb) = polys
    A = l1_eval(pa, y) / ma
    B = l1_eval(pb, y) / mb
    return A + arb(0, errA), B + arb(0, errB)


def cover(polys, errA, errB):
    """min A, max |B| on X; max G and max -G with G = B/A, by arb branch and bound (enclosure over balls)."""
    ctx.prec = 128
    lo_x, hi_x = arb(5).sqrt() / 2, arb(5) / 2
    lo_f, hi_f = arb(lo_x.lower()), arb(hi_x.upper())
    w = (hi_f - lo_f) / COVER0
    boxes = [(lo_f + w * k, lo_f + w * (k + 1)) for k in range(COVER0)]

    def ball(a, b):
        return arb(((a + b) / 2).mid(), ((b - a) / 2 + arb("1e-30")).upper())
    Amin, Bmax = arb("inf"), arb(0)
    for a, b in boxes:
        A, B = AB(polys, ball(a, b), errA, errB)
        Amin = arb(min(Amin.lower(), A.lower()))
        Bmax = arb(max(Bmax.upper(), abs(B).upper()))
    res = {}
    pts_inside = lambda t: bool(t >= lo_x) and bool(t <= hi_x)    # noqa: E731
    for sgn, key in ((1, "Gp"), (-1, "Gm")):
        live = list(boxes)
        best_lo = arb("-inf")
        arg = None
        it = 0
        while live and it < 60:
            it += 1
            nxt = []
            ups = []
            for a, b in live:
                A, B = AB(polys, ball(a, b), errA, errB)
                G = sgn * B / A
                m = arb(((a + b) / 2).mid())
                if not pts_inside(m):
                    m = arb(lo_x.mid()) if not bool(m >= lo_x) else arb(hi_x.mid())
                    if not pts_inside(m):
                        m = arb(((a + b) / 2).mid())
                Am, Bm = AB(polys, m, errA, errB)
                gm = sgn * Bm / Am
                if pts_inside(m) and bool(gm.lower() > best_lo):
                    best_lo, arg = arb(gm.lower()), m
                ups.append((G.upper(), a, b))
            for u, a, b in ups:
                if u >= best_lo:
                    if bool(b - a < _a(WMIN)):
                        nxt.append((a, b, u))
                    else:
                        mid = arb(((a + b) / 2).mid())
                        nxt.append((a, mid, None))
                        nxt.append((mid, b, None))
            done = [t for t in nxt if t[2] is not None]
            todo = [(a, b) for a, b, u in nxt if u is None]
            if not todo:
                up = max(arb(t[2]) for t in done) if done else best_lo
                res[key] = (best_lo, arb(up.upper()), arg)
                break
            live = todo + [(a, b) for a, b, _ in done]
            if all(bool(b - a < _a(WMIN)) for a, b in live):
                up = arb(0)
                ups2 = []
                for a, b in live:
                    A, B = AB(polys, ball(a, b), errA, errB)
                    ups2.append((sgn * B / A).upper())
                up = max(arb(x) for x in ups2)
                res[key] = (best_lo, arb(up.upper()), arg)
                break
        if key not in res:
            res[key] = (best_lo, arb("inf"), arg)
    return Amin, Bmax, res


def theta(n, q, rho, tau1, Amin, Bmax, gmax):
    """theta_n with |Y_n(x) - (-1)^n G(x)| <= theta_n; tau1 = e_X K_in (1/100) so tau_n = tau1 kappa^n."""
    kap = _a(RHO_IN) / rho
    aq = abs(q)
    tn, tn1 = tau1 * kap ** n, tau1 * kap ** (n + 1)
    eps_n = (aq ** n * Bmax + tn) / Amin
    Dn = (tn + tn1) / ((1 - q) * aq ** n * Amin)
    if not bool(eps_n < 1):
        return None, eps_n
    return (gmax * eps_n + Dn) / (1 - eps_n), eps_n


def float_observation(M, rho, q, Gp, Gm):
    import numpy as np
    Mf = np.array([[float(M[i, j].mid()) for j in range(M.ncols())] for i in range(M.nrows())])
    c, R = 83 / 50, 8 / 5
    xs = np.linspace(np.sqrt(5) / 2, 2.5, 2001)
    g1, g2 = np.sqrt(15 / 4 - xs), np.sqrt(15 / 4 + xs)
    w1 = (xs + .5) / (2 * xs * (3.5 - xs)) / (2 * g1)
    w2 = (xs - .5) / (2 * xs * (xs + 3.5)) / (2 * g2)

    def ev(p):
        return w1 * np.polyval(p[::-1], (g1 - c) / R) + w2 * np.polyval(p[::-1], (g2 - c) / R)
    rf, qf = float(rho.mid()), float(q.mid())
    p = np.zeros(M.nrows()); p[0] = 1.0
    hcur = ev(p)            # h_1 = L_1 1  (h_0 = 1)
    hprev = np.ones_like(xs)
    out = []
    for n in range(0, 13):
        R_n = hcur / hprev
        y = (rf - R_n.min()) / (rf * (1 - qf) * abs(qf) ** n)
        out.append((n, y))
        hprev = hcur
        p = Mf @ p
        hcur = ev(p)
    s = ", ".join(f"{n}:{y:.6f}" for n, y in out)
    print(f"   observation (floats, no claim): (rho - min R_n)/(rho(1-q)|q|^n), n = 0..12: {s}; Gamma_+ = {Gp.str(8)}, "
          f"Gamma_- = {Gm.str(8)}")
    return out


def main() -> int:
    t0 = time.time()
    print("Theorem CT: the exact rate of inf R_n -> rho(L_1): |nu_2|/rho(L_1), with explicit constants. Not thrust.", flush=True)
    dc = disk_constants()
    eps, beta = tail_constants(dc, N)
    M = build_matrix(N)
    sp = spectral_split(M)
    if sp is None:
        return 1
    pairs, T, alpha = sp
    circ, rho, nu, q, ok = certify(pairs, T, alpha, eps, beta)
    print(f"   [{time.time() - t0:.0f}s]", flush=True)
    if not ok:
        return 1
    ctx.prec = PREC
    eX = ev_X()
    errs = []
    for p in ("rhoproj", "nu2proj"):
        KM, KT, KL, cen, rad = circ[p]
        errs.append(arb((eX * rad * KL * KT * eps).upper()))
    KL_in = circ["in"][2]
    tau1 = arb((eX * KL_in * _a(RHO_IN)).upper())
    check(bool(errs[0] < arb("1e-12")) and bool(errs[1] < arb("1e-12")),
          f"CT.4: e_X = {eX.str(6)}; pointwise |A - A^N| <= {errs[0].str(3)}, |B - B^N| <= {errs[1].str(3)} on X "
          f"(e_X rad K K_T eps_N); sup_X |E_n| <= tau_1 rho_in^n, tau_1 = e_X K_in/100 = {tau1.str(4)}")
    polys = riesz_polys(pairs)
    Amin, Bmax, ext = cover(polys, errs[0], errs[1])
    Gp_lo, Gp_hi, xp = ext["Gp"]
    Gm_lo, Gm_hi, xm = ext["Gm"]
    Gp = arb(((Gp_lo + Gp_hi) / 2).mid(), ((Gp_hi - Gp_lo) / 2).upper())
    Gm = arb(((Gm_lo + Gm_hi) / 2).mid(), ((Gm_hi - Gm_lo) / 2).upper())
    check(bool(Amin > 0) and Bmax.is_finite(), f"CT.4: min_X A >= {Amin.str(6)} > 0, sup_X |B| <= {Bmax.str(6)} (cover of X)")
    check(bool(Gp_lo > 0) and bool(Gm_lo > 0) and Gp_hi.is_finite() and Gm_hi.is_finite(),
          f"CT.4: Gamma_+ = max_X B/A in [{Gp_lo.str(10)}, {Gp_hi.str(10)}] (near s = {xp.str(8)}), Gamma_- = max_X(-B/A) in "
          f"[{Gm_lo.str(10)}, {Gm_hi.str(10)}] (near s = {xm.str(8)}): both > 0, B/A changes sign on X")
    gmax = arb(max(Gp_hi, Gm_hi).upper())
    c0 = rho * (1 - q)
    n0 = None
    th = {}
    for n in range(NMAX + 1):
        t, en = theta(n, q, rho, tau1, Amin, Bmax, gmax)
        th[n] = t
        if n0 is None and t is not None and bool(t < Gp_lo) and bool(t < Gm_lo):
            n0 = n
    check(n0 is not None, f"CT.5: theta_n < min(Gamma_+, Gamma_-) from n_0 = {n0} on (theta_n decreasing: |q|^n and "
                          f"(rho_in/(rho|q|))^n terms, eps_n < 1)")
    if n0 is not None:
        for n in range(n0, NMAX):
            if not bool(th[n + 1] <= th[n]):
                check(False, f"CT.5: theta not monotone at {n}")
                break
        cp_ = c0 * Gp
        cm_ = c0 * Gm
        lines = []
        for n in sorted({n0, n0 + 1, 20, 21, 40, 41}):
            G = Gp if n % 2 == 0 else Gm
            lo = c0 * (arb(G.lower()) - th[n])
            hi = c0 * (arb(G.upper()) + th[n])
            lines.append(f"n={n}: [{arb(lo.lower()).str(6)}, {arb(hi.upper()).str(6)}]")
        check(bool(cp_ > 0) and bool(cm_ > 0),
              f"CT.5: for every n >= {n0}: rho(1-q)|q|^n (Gamma_par - theta_n) <= rho - inf_X R_n <= rho(1-q)|q|^n (Gamma_par + "
              f"theta_n); limits (rho - inf R_(2m))/|q|^(2m) -> c_+ = rho(1-q)Gamma_+ = {cp_.str(8)}, "
              f"(rho - inf R_(2m+1))/|q|^(2m+1) -> c_- = rho(1-q)Gamma_- = {cm_.str(8)}; exact rate |q| = {abs(q).str(10)} "
              f"(CJ.2's Birkhoff rate ~ 0.9974). Windows for (rho - inf R_n)/|q|^n: " + "; ".join(lines))
        osc = c0 * (Gp + Gm)
        olo = c0 * (arb(Gp.lower()) + arb(Gm.lower()) - 2 * th[n0])
        check(bool(olo > 0), f"CT.5: oscillation: for every n >= {n0}, rho(1-q)|q|^n (Gamma_+ + Gamma_- - 2 theta_n) <= sup_X R_n - "
                             f"inf_X R_n <= rho(1-q)|q|^n (Gamma_+ + Gamma_- + 2 theta_n), lower constant at n_0 >= {olo.str(5)} > 0; "
                             f"(sup R_n - inf R_n)/|q|^n -> rho(1-q)(Gamma_+ + Gamma_-) = {osc.str(8)} (both parities), so "
                             f"log(sup R_n/inf R_n) has the exact rate |q| (CJ.2: tanh(Delta/4)^(n-1) Delta)")
        RESULT.update(cp=cp_, cm=cm_, osc=osc, n0=n0, theta_n0=th[n0])
    RESULT.update(rho=rho, nu=nu, q=q, Gp=Gp, Gm=Gm, Amin=Amin, Bmax=Bmax, eps=eps, errs=errs, tau1=tau1)
    obs = float_observation(M, rho, q, Gp, Gm)
    RESULT.update(obs=obs)
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(cr.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem CT checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
