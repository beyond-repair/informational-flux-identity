#!/usr/bin/env python3
"""Theorem CU: the exact rate of inf R_n^(k) -> rho(L_k) for k = 2, 3, 4 (CT's open item "L_k, k >= 2"), with explicit
two-sided constants.

Setting (GASKET.md, Theorems CH, CK, CL, CT): in s = sqrt(t + 25/4), X = [sqrt5/2, 5/2], sigma_(1,2)(s) = sqrt(15/4 -+ s),
    (L_k f)(s) = u_1(s) f(sigma_1 s) + u_2(s) f(sigma_2 s),   u_i = wt_i |sigma_i'|^k = wt_i / (2 sigma_i)^k,
    wt_1 = (s+1/2)/(2s(7/2-s)),  wt_2 = (s-1/2)/(2s(s+7/2)),
h_n = L_k^n 1 > 0, R_n = R_n^(k) = h_(n+1)/h_n on X.  CL.3 certified rho(L_k) on C(X) to ~1e-7 (two-sided Collatz-Wielandt);
CT proved the exact rate of inf R_n -> rho(L_1) and listed k >= 2 as open.  Here L_k acts on CG's H^2(D), D = {|s - c| < R},
c = 83/50, R = 8/5, e_n = u^n, M_N = P_N L_k P_N, N = 120, 512-bit arb; inner radius rho_in = 1/500, 1/1000, 1/2500 for
k = 2, 3, 4.  The routine is CT's script with the weight exponent k as a parameter (CT's script itself is unchanged).

Theorem CU (checked here for k = 2, 3, 4; exit 1 on any failure).
  CU.1 [proved; exact + arb] |sigma_i| >= 7/10 on the closed disk (CT.1), so u_i is holomorphic there; 2048 arb arcs give
       sup|u_1| <= 1.2283, 0.87833, 0.62835 and sup|u_2| <= 0.068530, 0.017558, 0.0044986, hence
       ||L_k (I - P_120)|| <= 6.08e-27, 2.57e-27, 1.39e-27 (k = 2, 3, 4).
  CU.2 [certified; CR.1/CR.2's method] arb Krawczyk right/left eigenpairs of M_120 at rho and nu_2, real and disjoint; the
       certified Neumann bound ||A^7||_F <= 5.09e-21, 2.23e-23, 1.02e-24 < rho_in^7 puts every other eigenvalue of M_120 in
       |z| < rho_in.
  CU.3 [proved; CG.3's homotopy] K_T eps_N < 1 on |z| = rho_in, |z| = 6 and on circles of radius 1e-20 about rho and nu_2:
       L_k on H^2(D) has exactly two eigenvalues in |z| > rho_in, real and algebraically simple:
         k = 2: rho = 0.0430703112227101271 +- 1.3e-20, nu_2 = -0.0124717435618530418 +- 3.8e-20, q = -0.289567063896138 +- 2.7e-16
         k = 3: rho = 0.0139192547691468646 +- 2.1e-20, nu_2 = -0.0044401723728826436 +- 1.4e-20, q = -0.318994978288970 +- 1.9e-16
         k = 4: rho = 0.0045264690852562123 +- 2.8e-20, nu_2 = -0.0015282685887826225 +- 3.6e-20, q = -0.337629300012356 +- 4.2e-16
       and rho_H = rho(L_k) on C(X) (CT.3's argument: point evaluation and the analytic eigenfunction); each lies inside
       CL.3's interval (consistency), which it sharpens from ~1e-7 to ~3e-20.
  CU.4 [proved; second resolvent identity] h_n = rho^n A + nu_2^n B + E_n, sup_X |E_n| <= tau_1 rho_in^n with
       tau_1 = e_X K_in rho_in <= 8.346e4, 1.633e5, 1.453e6; A, B enclosed pointwise (errors <= 5.8e-17, 1.3e-15, 7.2e-14);
       an arb cover of X gives min_X A >= 0.795593, 0.759466, 0.724643 and, with G = B/A,
         k = 2: Gamma_+ = max G in [0.2546353828, 0.2546354891],  Gamma_- = max(-G) in [0.6796088690, 0.6796091443]
         k = 3: Gamma_+ in [0.3087742424, 0.3087743752],           Gamma_- in [0.8310492070, 0.8310495724]
         k = 4: Gamma_+ in [0.3625750669, 0.3625752290],           Gamma_- in [0.9801674176, 0.9801678993]
  CU.5 [proved from CU.3, CU.4; CT.5's algebra] theta_n (CT.5's formula, bases |q|, rho_in/rho, rho_in/|nu_2| all < 1)
       decreases and theta_n < min Gamma_+- from n_0 = 7, 9, 12.  For every n >= n_0:
           rho (1-q)|q|^n (Gamma_par - theta_n) <= rho(L_k) - inf_X R_n <= rho (1-q)|q|^n (Gamma_par + theta_n),
       explicitly
         k = 2: 0.0120710 |q|^n <= rho - inf R_n <= 0.0162149 |q|^n (even n >= 8), 0.0248320 .. 0.0506618 (odd n >= 7)
         k = 3: 0.00459104 .. 0.00674680 (even n >= 10), 0.0104718 .. 0.0200434 (odd n >= 9)
         k = 4: 0.00117440 .. 0.00321619 (even n >= 12), 0.00566745 .. 0.00620186 (odd n >= 13)
       so inf R_n^(k) -> rho(L_k) at the EXACT rate |nu_2/rho| = 0.2895670639, 0.3189949783, 0.3376293000, with limit
       constants c_+ = 0.01414298, 0.00566892, 0.002195297 (even n) and c_- = 0.03774688, 0.01525759, 0.00593466 (odd n);
       sup R_n - inf R_n has the same rate, limit constants 0.0518899, 0.02092651, 0.00812996.
Observation (no claim): floats of M_N^n e_0 match Gamma_+- from n ~ 6 on; the extrema of B/A sit at the endpoints
s = sqrt5/2 (Gamma_+) and s = 5/2 (Gamma_-), as for k = 1.
Kept limits: n < n_0 not covered; Gamma_+- known to ~1e-7; only k = 2, 3, 4 (no uniform-in-k statement); nu_2 is an
H^2(D) eigenvalue, nothing about the C^k spectrum of L_k is claimed beyond rho_H = rho(L_k).
Assumptions used as stated (not re-run): CJ.2's definition of R_n (here with L_k), CL.3's intervals (consistency only), CG.1's
disk and tail formula (constants re-run for the new weights), CG.3's homotopy argument, CR.1/CR.2's Krawczyk and Neumann
routines (imported unchanged), CT.4/CT.5's algebra.  Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW
enters.  Functions on a graph, not a flux.  Not thrust.
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
ALL: dict = {}
N = int(os.environ.get("CU_N", "120"))
KD = 2                 # derivative order k: L_k has weights wt_i |sigma_i'|^k; set by configure(k)
PREC = 512
C_CEN, R_RAD = Fr(83, 50), Fr(8, 5)
R_OUT = Fr(6)
NAMES = ["rho", "nu2"]
G_RAD = [Fr(1, 10 ** 20), Fr(1, 10 ** 20)]
# per-k data: float seeds (Krawczyk certifies), inner radius, projection radii, CL.3's interval (consistency only), |q| cap
CONF = {
    2: dict(seeds=["0.0430703112", "-0.0124717436"], rho_in=Fr(1, 500), proj=[Fr(1, 100), Fr(1, 200)],
            cl3=("0.0430702", "0.0430704"), qcap="0.29"),
    3: dict(seeds=["0.0139192548", "-0.00444017237"], rho_in=Fr(1, 1000), proj=[Fr(1, 200), Fr(1, 500)],
            cl3=("0.0139192", "0.0139193"), qcap="0.32"),
    4: dict(seeds=["0.00452646909", "-0.00152826859"], rho_in=Fr(1, 2500), proj=[Fr(1, 500), Fr(1, 1250)],
            cl3=("0.0045264", "0.0045265"), qcap="0.34"),
}
KS = [int(x) for x in os.environ.get("CU_KS", "2,3,4").split(",")]
RHO_IN = SEEDS = PROJ_RAD = CI1 = QCAP = None


def configure(k: int) -> None:
    global KD, RHO_IN, SEEDS, PROJ_RAD, CI1, QCAP
    c = CONF[k]
    KD, RHO_IN, SEEDS, PROJ_RAD, CI1, QCAP = k, c["rho_in"], c["seeds"], c["proj"], c["cl3"], c["qcap"]
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


# ----------------------------------------------------------------- CU.1 the disk
def disk_constants(K=K_ARC):
    ctx.prec = 128
    c, R = _a(C_CEN), _a(R_RAD)
    check(C_CEN - R_RAD > 0 and Fr(7, 2) - C_CEN > R_RAD and Fr(15, 4) - C_CEN - R_RAD == Fr(49, 100),
          "CU.1: poles 0, 7/2 of wt_i outside the closed disk; Re(15/4 - s) >= 49/100, so |sigma_i| >= 7/10 and u_i holomorphic")
    q4, h, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    pi = arb.pi()
    q1 = q2 = U1 = U2 = arb(0)
    for j in range(K):
        th = arb(2 * j + 1) * pi / K + arb(0, (pi / K).upper())
        s = c + R * acb(0, th).exp()
        g1, g2 = (q4 - s).sqrt(), (q4 + s).sqrt()
        v1 = abs((g1 - c) / R).upper()
        v2 = abs((g2 - c) / R).upper()
        w1 = abs((s + h) / (2 * s * (s7 - s)) / (2 * g1) ** KD).upper()
        w2 = abs((s - h) / (2 * s * (s + s7)) / (2 * g2) ** KD).upper()
        q1, q2, U1, U2 = (arb(max(q1, v1)), arb(max(q2, v2)), arb(max(U1, w1)), arb(max(U2, w2)))
    check(bool(q1 < 1) and bool(q2 < 1) and U1.is_finite() and U2.is_finite(),
          f"CU.1: sup_D |v_1| <= {q1.str(5)}, sup_D |v_2| <= {q2.str(5)}, sup_D |u_1| <= U_1 = {U1.str(5)}, "
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
    check(bool(eps < arb("1e-20")), f"CU.1: ||L_k (I - P_N)|| <= eps_N = {eps.str(4)}, ||L_k P_N|| <= beta_N = {beta.str(5)} (N = {n})")
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
    w1 = (s + arb(1) / 2) / (2 * s * (arb(7) / 2 - s)) / (2 * g1) ** KD
    w2 = (s - arb(1) / 2) / (2 * s * (s + arb(7) / 2)) / (2 * g2) ** KD
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
    """(L_k p)(y) for the polynomial p(s) = sum p_k ((s - c)/R)^k; y an arb (point or ball in X)."""
    c, R = _a(C_CEN), _a(R_RAD)
    q4, h, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    g1, g2 = (q4 - y).sqrt(), (q4 + y).sqrt()
    w1 = (y + h) / (2 * y * (s7 - y)) / (2 * g1) ** KD
    w2 = (y - h) / (2 * y * (y + s7)) / (2 * g2) ** KD

    def poly(t):
        uu = (t - c) / R
        acc = arb(0)
        for a in reversed(p):
            acc = acc * uu + a
        return acc
    return w1 * poly(g1) + w2 * poly(g2)


# ----------------------------------------------------------------- CU.2 / CU.3
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
            check(False, f"CU.2: Krawczyk enclosure for {name} failed")
            return None
        wv = (w.transpose() * v)[0, 0]
        Pn = arb((cr._fro(v) * cr._fro(w) / abs(wv)).upper())
        pairs.append(dict(name=name, v=v, w=w, mu=mu, wv=wv, P=Pn))
        check(not wv.contains(0), f"CU.2 (arb Krawczyk, {PREC} bits, N = {N}): {name}: mu_N = {mu.str(22)} (radii "
                                  f"{float(e1):.0e}/{float(e2):.0e}), w.v != 0, ||P|| <= {Pn.str(5)}")
    check(not pairs[0]["mu"].overlaps(pairs[1]["mu"]), "CU.2: the two enclosures are real and disjoint")
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
          f"CU.2 (arb): ||A^j Q||_F <= " + ", ".join(t.str(3) for t in T) + f", ||A^{cr.M_POW}||_F <= {alpha.str(3)} "
          f"< (rho_in)^{cr.M_POW}: every other eigenvalue of M_N has |z| < rho_in")
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
            check(False, f"CU.3: circle {name} not certified")
            out[name] = None
            continue
        KT = KM * (1 + beta / absmin) + 1 / absmin
        ok = bool(KT * eps < 1)
        KL = arb((KT / (1 - KT * eps)).upper()) if ok else None
        check(ok, f"CU.3: |z - {cen.str(12)}| = {rad.str(3)}: sup ||(z - M)^-1|| <= {KM.str(4)}, K_T = {KT.str(4)}, "
                  f"K_T eps_N = {(KT * eps).str(3)} < 1" + (f", sup ||(z - L_k)^-1|| <= {KL.str(4)}" if ok else ""))
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
    check(allok, f"CU.3: beta_N + eps_N < {R_OUT}; the small and projection circles hold mu_i only and stay in |z| > rho_in. "
                 f"CG.3's homotopy: L_k on H^2(D) has exactly two eigenvalues in |z| > rho_in, real and simple: "
                 f"rho_H = {rho.str(20)}, nu_2 = {nu.str(20)} (+- 1e-20); every other eigenvalue has |mu| < rho_in")
    lo, hi = arb(CI1[0]), arb(CI1[1])
    check(bool(rho > lo) and bool(rho < hi), f"CU.3: rho_H lies in CL.3's [{CI1[0]}, {CI1[1]}] (consistency); rho_H = rho(L_k) "
                                             "on C(X) (point evaluation <= e_X ||.||_H2 and the analytic eigenfunction is not 0 on X)")
    q = nu / rho
    check(bool(q < 0) and bool(abs(q) < arb(QCAP)), f"CU.3: q = nu_2/rho(L_k) = {q.str(15)}, |q| < {QCAP}")
    return out, rho, nu, q, allok


# ----------------------------------------------------------------- CU.4 pointwise A, B and the cover
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
    """theta_n with |Y_n(x) - (-1)^n G(x)| <= theta_n; tau1 = e_X K_in (rho_in) so tau_n = tau1 kappa^n."""
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
    w1 = (xs + .5) / (2 * xs * (3.5 - xs)) / (2 * g1) ** KD
    w2 = (xs - .5) / (2 * xs * (xs + 3.5)) / (2 * g2) ** KD

    def ev(p):
        return w1 * np.polyval(p[::-1], (g1 - c) / R) + w2 * np.polyval(p[::-1], (g2 - c) / R)
    rf, qf = float(rho.mid()), float(q.mid())
    p = np.zeros(M.nrows()); p[0] = 1.0
    hcur = ev(p)            # h_1 = L_k 1  (h_0 = 1)
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


def run_k(k: int) -> int:
    configure(k)
    RESULT.clear()
    t0 = time.time()
    print(f"--- k = {k}: L_k with weights wt_i |sigma_i'|^{k}, rho_in = {RHO_IN}", flush=True)
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
          f"CU.4: e_X = {eX.str(6)}; pointwise |A - A^N| <= {errs[0].str(3)}, |B - B^N| <= {errs[1].str(3)} on X "
          f"(e_X rad K K_T eps_N); sup_X |E_n| <= tau_1 rho_in^n, tau_1 = e_X K_in rho_in = {tau1.str(4)}")
    polys = riesz_polys(pairs)
    Amin, Bmax, ext = cover(polys, errs[0], errs[1])
    Gp_lo, Gp_hi, xp = ext["Gp"]
    Gm_lo, Gm_hi, xm = ext["Gm"]
    Gp = arb(((Gp_lo + Gp_hi) / 2).mid(), ((Gp_hi - Gp_lo) / 2).upper())
    Gm = arb(((Gm_lo + Gm_hi) / 2).mid(), ((Gm_hi - Gm_lo) / 2).upper())
    check(bool(Amin > 0) and Bmax.is_finite(), f"CU.4: min_X A >= {Amin.str(6)} > 0, sup_X |B| <= {Bmax.str(6)} (cover of X)")
    check(bool(Gp_lo > 0) and bool(Gm_lo > 0) and Gp_hi.is_finite() and Gm_hi.is_finite(),
          f"CU.4: Gamma_+ = max_X B/A in [{Gp_lo.str(10)}, {Gp_hi.str(10)}] (near s = {xp.str(8)}), Gamma_- = max_X(-B/A) in "
          f"[{Gm_lo.str(10)}, {Gm_hi.str(10)}] (near s = {xm.str(8)}): both > 0, B/A changes sign on X")
    gmax = arb(max(Gp_hi, Gm_hi).upper())
    b1, b2 = _a(RHO_IN) / rho, _a(RHO_IN) / abs(nu)
    check(bool(abs(q) < 1) and bool(b1 < 1) and bool(b2 < 1),
          f"CU.5: bases |q| = {abs(q).str(6)}, rho_in/rho = {b1.str(4)}, rho_in/|nu_2| = {b2.str(4)} all < 1, so eps_n, Delta_n and "
          f"theta_n decrease in n")
    c0 = rho * (1 - q)
    n0 = None
    th = {}
    for n in range(NMAX + 1):
        t, en = theta(n, q, rho, tau1, Amin, Bmax, gmax)
        th[n] = t
        if n0 is None and t is not None and bool(t < Gp_lo) and bool(t < Gm_lo):
            n0 = n
    check(n0 is not None, f"CU.5: theta_n < min(Gamma_+, Gamma_-) from n_0 = {n0} on (theta_n decreasing: |q|^n and "
                          f"(rho_in/(rho|q|))^n terms, eps_n < 1)")
    if n0 is not None:
        for n in range(n0, NMAX):
            if not bool(th[n + 1] <= th[n]):
                check(False, f"CU.5: theta not monotone at {n}")
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
              f"CU.5: for every n >= {n0}: rho(1-q)|q|^n (Gamma_par - theta_n) <= rho - inf_X R_n <= rho(1-q)|q|^n (Gamma_par + "
              f"theta_n); limits (rho - inf R_(2m))/|q|^(2m) -> c_+ = rho(1-q)Gamma_+ = {cp_.str(8)}, "
              f"(rho - inf R_(2m+1))/|q|^(2m+1) -> c_- = rho(1-q)Gamma_- = {cm_.str(8)}; exact rate |q| = {abs(q).str(10)} "
              f"Windows for (rho - inf R_n)/|q|^n: " + "; ".join(lines))
        ne = n0 if n0 % 2 == 0 else n0 + 1
        no = n0 if n0 % 2 == 1 else n0 + 1
        ev_lo, ev_hi = arb((c0 * (arb(Gp.lower()) - th[ne])).lower()), arb((c0 * (arb(Gp.upper()) + th[ne])).upper())
        od_lo, od_hi = arb((c0 * (arb(Gm.lower()) - th[no])).lower()), arb((c0 * (arb(Gm.upper()) + th[no])).upper())
        check(bool(ev_lo > 0) and bool(od_lo > 0),
              f"CU.5: uniform constants (theta_n decreasing): {ev_lo.str(6)} |q|^n <= rho - inf R_n <= {ev_hi.str(6)} |q|^n for "
              f"even n >= {ne}; {od_lo.str(6)} |q|^n <= rho - inf R_n <= {od_hi.str(6)} |q|^n for odd n >= {no}")
        RESULT.update(ev=(ev_lo, ev_hi, ne), od=(od_lo, od_hi, no))
        osc = c0 * (Gp + Gm)
        olo = c0 * (arb(Gp.lower()) + arb(Gm.lower()) - 2 * th[n0])
        check(bool(olo > 0), f"CU.5: oscillation: for every n >= {n0}, rho(1-q)|q|^n (Gamma_+ + Gamma_- - 2 theta_n) <= sup_X R_n - "
                             f"inf_X R_n <= rho(1-q)|q|^n (Gamma_+ + Gamma_- + 2 theta_n), lower constant at n_0 >= {olo.str(5)} > 0; "
                             f"(sup R_n - inf R_n)/|q|^n -> rho(1-q)(Gamma_+ + Gamma_-) = {osc.str(8)} (both parities), so "
                             f"log(sup R_n/inf R_n) has the exact rate |q|")
        RESULT.update(cp=cp_, cm=cm_, osc=osc, n0=n0, theta_n0=th[n0])
    RESULT.update(k=k, rho=rho, nu=nu, q=q, Gp=Gp, Gm=Gm, Amin=Amin, Bmax=Bmax, eps=eps, errs=errs, tau1=tau1)
    obs = float_observation(M, rho, q, Gp, Gm)
    RESULT.update(obs=obs)
    print(f"   k = {k}: elapsed {time.time() - t0:.1f}s", flush=True)
    return 0


def main() -> int:
    print("Theorem CU: the exact rate of inf R_n^(k) -> rho(L_k) for k = 2, 3: |nu_2^(k)|/rho(L_k), explicit constants. "
          "Not thrust.", flush=True)
    t0 = time.time()
    ALL.clear()
    for k in KS:
        rc = run_k(k)
        ALL[k] = dict(RESULT)
        if rc:
            check(False, f"CU: k = {k} aborted")
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(cr.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem CU checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
