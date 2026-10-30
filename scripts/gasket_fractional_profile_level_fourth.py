#!/usr/bin/env python3
r"""Theorem CQ: the fourth eigenvalue l_4 of CE's transfer operator, and a_4 != 0 in the level-sum expansion.

Setting (GASKET.md, Theorems CE-CG, CN, CO): (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s) on H^2(D),
D = {|s - c| < R}, c = 83/50, R = 8/5, e_n = u^n, u = (s - c)/R; M_N = P_N Lt P_N, ||Lt (I - P_N)|| <= eps_N (CG.1's
formula, any N).  F0 = Lt^4 G (CN.1), x0 = sqrt13/2, c0 = 3 sigma_-3, a_j = c0 (Pi_j F0)(x0).  CG.3: Lt has exactly three
eigenvalues r, l_2, l_3 of modulus > 1e-3; every other eigenvalue has |mu| < 1e-3.  CN.2 put all of |z| < 1e-3 into one
remainder e_k, |e_k| <= eps_0 1e-3^k.  CO.4 enclosed a_1, a_2, a_3 (N = 80).

Theorem CQ (checked here; exit 1 on any failure):
  CQ.1 [certified; arb charpoly at 1100 bits, Pellet] chi(M_80) passes Pellet's test at rho = 1/1000 (k = N - 2, as CG.2)
       and at rho' = 1/20000 (k = N - 3): M_80 has exactly one eigenvalue mu_4^80 in rho' < |z| < 1/1000, real.
  CQ.2 [certified; rigorous float residual bounds + arb] Resolvent circles for M_80 on G_4 = {|z - m_4| = 1/20000} and on
       |z| = rho': sup ||(z - Lt P_80)^-1|| eps_80 < 1.  So (homotopy z - Lt P_N - tau Lt (I - P_N), as CG.3) Lt has exactly
       one eigenvalue l_4 inside G_4, real and algebraically simple; with CG.3's |z| = 1/1000 (as stated) every eigenvalue
       of Lt other than r, l_2, l_3, l_4 has |mu| < rho' = 5e-5.
  CQ.3 [certified] N = 140: eps_140 from CG.1's formula, Krawczyk right/left eigenvectors of M_140 at mu_4^140 (inside G_4),
       G_4 circle bound for M_140.  Rank-one split (z - M)^-1 = P/(z - mu) + H(z), H holomorphic on the closed G_4 disk:
       by the maximum principle ||(z - M)^-1|| <= ||P||/delta + K_G4 + ||P||/(rad_4 - |mu - m_4|) on |z - mu_4^140| = delta = 1e-18, and
       that bound times eps_140 is < 1.  So l_4 lies within 1e-18 of mu_4^140.
  CQ.4 [proved; certified] CO.3/CO.4's exact projection formula for Lt P_140 and the second-resolvent bound on G_4:
           |a_4 - a_4^N| <= c0 ||ev_x0|| rad_4 K K^N eps_140 S,   so a_4 != 0.
       Taylor coefficients f_0..f_140 as in CO.2 (2048 certified points, aliasing bound).
  CQ.5 [proved from CQ.2, CQ.4, CN.1] s_(k+4) = a_1 r^k + a_2 l_2^k + a_3 l_3^k + a_4 l_4^k + e'_k for every k >= 0,
       |e'_k| <= eps_0' rho'^k, eps_0' = c0 ||ev_x0|| rho' K' S.  Hence CN.2's remainder e_k = a_4 l_4^k + e'_k satisfies
       e_k / l_4^k -> a_4 != 0 and sign(e_k) = (-1)^k sign(a_4) for every k >= k_0 (explicit): its exact rate is |l_4|.
Rigorous float residual bound (CQ.2/CQ.3 circles): for a float approximate inverse Y of C = fl(zI - A), A = float mid of
the arb matrix M with ||M - A||_F <= dA (arb), E = I - Y(zI - M) has ||E||_F <= ||I - fl(YC)||_F (1+u) + g ||fl(|Y||C|)||_F
(1 + g) + ||Y||_F (dA + u ||diag C||), g = sqrt2 (n+2) u / (1 - (n+2) u) (Higham, Thm 3.5 and Sec. 3.6, complex dot
products), and ||(zI - M)^-1|| <= ||Y||_F / (1 - ||E||_F); all float norms are inflated by 1 + 1e-9.  Adaptive cover as CG.3.
Kept limits: no closed form of l_4 or a_4; the next eigenvalues (floats: -1.744e-7 and a complex pair near
-1.16e-8 +- 1.55e-7 i) are only bounded by rho' (Pellet at 1e-5 fails); CP is unchanged (its eps_n already contains the
l_4-term).  C_G, circle eigenvalues |z| = rho(L_k), inf R_n rate stay open.
Assumptions used as stated (not re-run): CE.1-CE.3, CG.3's circle |z| = 1/1000 and its count, CN.1's box argument (box
checks and S re-run), CO.1's continuation (re-run).  Textbook: Pellet, Krawczyk/Brouwer, Riesz projections, second
resolvent identity, maximum principle, Cauchy/trapezoid aliasing, floating-point dot-product error bounds.  Stage 1
untouched; no W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time
from fractions import Fraction as Fr

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")      # small matrices: threading only adds overhead
import numpy as np  # noqa: E402
from flint import arb, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_coefficient as cn     # noqa: E402
import gasket_fractional_profile_level_spectrum as cg        # noqa: E402
import gasket_fractional_profile_level_third as co           # noqa: E402

FAILS: list[str] = []
N1 = 80                    # counting (CG's matrix)
N2 = 140                   # projection and localisation
RHO_P = Fr(1, 20000)       # inner radius rho'
RHO_OUT = Fr(1, 1000)      # CG.3's circle
L4_FLOAT = -2.1644655284e-4
RAD4 = Fr(1, 20000)
DELTA4 = Fr(1, 10 ** 18)
U = 2.0 ** -53
INFL = 1 + 1e-9


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


_a = cn._a


def _up(x: arb) -> float:
    """A float >= x (x a finite arb)."""
    return float(arb(x.upper()).mid()) * INFL + 1e-300


# ----------------------------------------------------------------- CQ.1
def pellet(M):
    ctx.prec = cg.PREC_HI
    a = M.charpoly().coeffs()
    oks = []
    for rho, cnt in ((RHO_OUT, 3), (RHO_P, 4)):
        k = N1 + 1 - cnt
        rr = _a(rho)
        lhs = abs(a[k]).lower() * rr ** k
        rhs = arb(0)
        for j in range(len(a)):
            if j != k:
                rhs += arb(abs(a[j]).upper()) * rr ** j
        ok = bool(lhs > rhs)
        oks.append(ok)
        check(ok, f"CQ.1 (arb, {cg.PREC_HI} bits): Pellet at rho = {rho}, k = {k}: |a_k| rho^k / sum_(j != k) |a_j| rho^j "
                  f"= {(lhs / rhs).str(5)} > 1, so M_{N1} has exactly {cnt} eigenvalues with |z| > {rho}")
    check(all(oks), f"CQ.1: M_{N1} has exactly one eigenvalue in {RHO_P} < |z| < {RHO_OUT}, and it is real (chi real: a "
                    "non-real root would come with its conjugate)")
    return all(oks)


# ----------------------------------------------------------------- rigorous float resolvent bounds
class FastRes:
    def __init__(self, M):
        n = M.nrows()
        ctx.prec = 320
        self.n = n
        self.A = np.array([[float(M[i, j].mid()) for j in range(n)] for i in range(n)])
        d2 = arb(0)
        for i in range(n):
            for j in range(n):
                e = abs(M[i, j] - arb(self.A[i, j]))
                d2 += arb(e.upper()) ** 2
        self.dA = _up(d2.sqrt())
        nn = n + 2
        self.g = math.sqrt(2) * nn * U / (1 - nn * U) * INFL
        self.I = np.eye(n)

    def norm_inv(self, z: complex):
        """Rigorous upper bound for ||(z I - M)^-1||_2 (via Frobenius), or inf."""
        C = z * self.I - self.A
        Y = np.linalg.inv(C)
        R = Y @ C
        D = self.I - R
        G = np.abs(Y) @ np.abs(C)
        nY = np.linalg.norm(Y) * INFL
        e1 = np.linalg.norm(D) * INFL * (1 + U)
        e2 = self.g * np.linalg.norm(G) * INFL * (1 + self.g)
        dT = self.dA + U * np.linalg.norm(np.diag(C)) * INFL
        E = (e1 + e2 + nY * dT) * INFL
        if not E < 0.5:
            return math.inf
        return nY / (1 - E) * INFL


def circle_fast(fr: FastRes, cen: float, radius: Fr, eps: arb, beta: arb, max_samples=2_000_000):
    """CG.3's adaptive cover with rigorous float bounds.  Returns (ok, sup ||(z - Lt P_N)^-1||, sup ||(z - M)^-1||, n)."""
    rad = float(radius)
    absmin = (abs(cen) - rad) * (1 - 1e-12) if abs(cen) > rad else rad * (1 - 1e-12)
    b = _up(beta)
    ang, n, worstM = 0.0, 0, 0.0
    while ang < 2 * math.pi * (1 + 1e-12):
        z = complex(cen + rad * math.cos(ang), rad * math.sin(ang))
        K = fr.norm_inv(z)
        if not math.isfinite(K):
            return False, math.inf, math.inf, n
        worstM = max(worstM, 2 * K)
        n += 1
        if n > max_samples:
            return False, math.inf, worstM, n
        # z is within 1e-15 * |z| of the exact circle point; the ball |z' - z| <= 1/(2K) covers the next arc
        step = (1 / (2 * K) - 4e-16 * (abs(cen) + rad)) / rad * 0.999
        if step <= 0:
            return False, math.inf, worstM, n
        ang += step
    KT = worstM * (1 + b / absmin) * INFL + 1 / absmin * INFL
    return bool(arb(KT) * eps < 1), KT, worstM, n


# ----------------------------------------------------------------- CQ.2 counting with N = 80
def counting(M80):
    dc = cg.disk_constants()
    eps80, beta80 = cg.tail_constants(dc, N1)
    fr = FastRes(M80)
    out = {}
    m4 = L4_FLOAT
    for name, cen, rad in (("G4", m4, RAD4), ("in", 0.0, RHO_P)):
        ok, KT, KM, n = circle_fast(fr, cen, rad, eps80, beta80)
        KL = KT / (1 - KT * float(eps80.upper())) * INFL if ok else math.inf
        check(ok, f"CQ.2 (rigorous floats, N = {N1}, ||M - A||_F <= {fr.dA:.2e}): |z - {cen:.6e}| = {rad}: sup ||(z - M)^-1|| <= "
                  f"{KM:.4e}, sup ||(z - Lt P_N)^-1|| <= {KT:.4e}, times eps_{N1} = {KT * float(eps80.upper()):.2e} < 1, so "
                  f"sup ||(z - Lt)^-1|| <= {KL:.4e} ({n} samples)")
        out[name] = (ok, KT, KL)
    lay = abs(m4) - float(RAD4) > float(RHO_P) and abs(m4) + float(RAD4) < float(RHO_OUT)
    check(lay and all(v[0] for v in out.values()),
          f"CQ.2: G_4 = |z - {m4:.6e}| = {RAD4} lies in {RHO_P} < |z| < {RHO_OUT}; it holds mu_4^{N1} (CQ.1, checked in CQ.3) and "
          f"no other eigenvalue of M_{N1}, so Lt has exactly one eigenvalue l_4 inside G_4: real (Lt commutes with "
          f"f -> conj f(conj s), G_4 symmetric) and algebraically simple.  With CG.3's |z| = 1/1000 circle (as stated) and "
          f"|z| = {RHO_P}: every eigenvalue of Lt other than r, l_2, l_3, l_4 has |mu| < {RHO_P}")
    return out, eps80, beta80, dc


# ----------------------------------------------------------------- CQ.3 localisation with N = 140
def localise(M140, M80, dc):
    eps, beta = cg.tail_constants(dc, N2)
    ctx.prec = 320
    v, mu, e1 = co.eigvec(M140, L4_FLOAT, N2 + 1)
    w, muL, e2 = co.eigvec(M140.transpose(), L4_FLOAT, N2 + 1)
    v80, mu80, _ = co.eigvec(M80, L4_FLOAT, N1 + 1)
    ok = v is not None and w is not None and mu80 is not None
    if not ok:
        check(False, "CQ.3: Krawczyk enclosures failed")
        return None
    inG = bool(abs(mu - arb(L4_FLOAT)) < _a(RAD4)) and bool(abs(mu80 - arb(L4_FLOAT)) < _a(RAD4)) and bool(mu.overlaps(muL))
    check(inG, f"CQ.3 (arb Krawczyk, 320 bits): mu_4^{N2} = {mu.str(22)} (right/left eigenvectors, radius {float(e1):.0e}/"
               f"{float(e2):.0e}) and mu_4^{N1} = {mu80.str(20)} both lie inside G_4")
    fr = FastRes(M140)
    okc, KT, KM, n = circle_fast(fr, L4_FLOAT, RAD4, eps, beta)
    KL = KT / (1 - KT * float(eps.upper())) * INFL
    check(okc, f"CQ.3 (rigorous floats, N = {N2}): G_4: sup ||(z - M_{N2})^-1|| <= {KM:.4e}, sup ||(z - Lt P_N)^-1|| <= {KT:.4e}, "
               f"times eps_{N2} = {eps.str(4)} < 1, so sup ||(z - Lt)^-1|| <= {KL:.4e} ({n} samples); since Lt has one "
               f"eigenvalue inside G_4 (CQ.2), so does M_{N2}: mu_4^{N2}, simple")
    ctx.prec = 320
    wv = sum((w[i] * v[i] for i in range(N2 + 1)), arb(0))
    nv = sum((x * x for x in v), arb(0)).sqrt()
    nw = sum((x * x for x in w), arb(0)).sqrt()
    Pn = nv * nw / abs(wv)
    d = _a(DELTA4)
    rr = _a(RAD4) - abs(mu - arb(L4_FLOAT))       # |z - mu| >= rr on G_4
    Kd = Pn / d + arb(KM) + Pn / rr
    zmin = abs(mu) - d
    KTd = Kd * (1 + beta / zmin) + 1 / zmin
    okd = bool(KTd * eps < 1) and bool(zmin > 0)
    l4 = arb(mu.mid()) + arb(0, (d + mu.rad()).upper())
    check(okd, f"CQ.3 (arb): ||P|| = ||v|| ||w|| / |w.v| = {Pn.str(6)}; on |z - mu_4^{N2}| = {DELTA4}, ||(z - M)^-1|| <= "
               f"||P||/delta + K_G4 + ||P||/(rad_4 - |mu - m_4|) = {Kd.str(4)} (maximum principle for H), so ||(z - Lt P_N)^-1|| eps_{N2} <= "
               f"{(KTd * eps).str(3)} < 1: l_4 = {l4.str(18)} (+- {float(d):.0e}); |l_4|/r = {(abs(l4) / arb(cn.R_MID)).str(10)}; "
               f"|l_4 - mu_4^{N1}| = {abs(l4 - mu80).str(3)}")
    return dict(v=v, w=w, mu=mu, wv=wv, eps=eps, beta=beta, KT=KT, KL=KL, l4=l4, v80=v80, mu80=mu80)


# ----------------------------------------------------------------- CQ.4, CQ.5
def a4_coeff(f, S, loc):
    ctx.prec = 320
    import gasket_fractional_profile_level_ratio as ce
    c0 = 3 * ce.sigma_closed(arb(-3), 0)
    x0 = arb(13).sqrt() / 2
    u0 = (x0 - _a(cn.C_CEN)) / _a(cn.R_RAD)
    ev = 1 / (1 - u0 * u0).sqrt()
    v, w, mu, wv = loc["v"], loc["w"], loc["mu"], loc["wv"]
    wf = sum((w[i] * f[i] for i in range(N2 + 1)), arb(0))
    phi = co.lt_eval(v, x0) / mu
    aN = c0 * wf / wv * phi
    err = c0 * ev * _a(RAD4) * arb(loc["KL"]) * arb(loc["KT"]) * loc["eps"] * S
    a4 = aN + arb(0, err.upper())
    print(f"   a_4^N = {aN.str(15)}  (c0 = {c0.str(10)}, (w.f)/(w.v) = {(wf / wv).str(8)}, phi(x0) = {phi.str(8)})")
    check(not a4.contains(0) and bool(a4 < 0),
          f"CQ.4: |a_4 - a_4^N| <= c0 ||ev|| rad_4 K K^N eps_{N2} S = {err.str(3)}; a_4 = {a4.str(8)} != 0 (negative): the "
          "l_4-component of the level sums is not zero")
    return a4, c0, ev


def rate(a4, c0, ev, S, KLin, l4):
    ctx.prec = 320
    rho = _a(RHO_P)
    eps0p = c0 * ev * rho * arb(KLin) * S
    th = rho / abs(l4)
    k0 = None
    for k in range(200):
        if bool(eps0p * th ** k < abs(a4)):
            k0 = k
            break
    ok = k0 is not None and bool(th < 1)
    check(ok, f"CQ.5: s_(k+4) = a_1 r^k + a_2 l_2^k + a_3 l_3^k + a_4 l_4^k + e'_k, |e'_k| <= eps_0' rho'^k, eps_0' = c0 ||ev|| "
              f"rho' K' S = {eps0p.str(4)}; rho'/|l_4| = {th.str(6)} < 1 and eps_0' (rho'/|l_4|)^k < |a_4| for k >= k_0 = {k0}.  So "
              f"CN.2's remainder e_k = a_4 l_4^k + e'_k has e_k / l_4^k -> a_4 and sign(e_k) = (-1)^(k+1) for every k >= {k0}: "
              f"exact rate |l_4| (not 1e-3)")
    return eps0p, k0


def main(arcs=4096, M_pts=co.M_PTS) -> int:
    t0 = time.time()
    print("Theorem CQ: the fourth eigenvalue l_4 and a_4 != 0 in the level-sum expansion. Not thrust.", flush=True)
    cn.tail_constants()
    cn.box_checks()
    S = cn.circle_sup(1024)
    pellet(cg.build_matrix(N1))
    M80 = cg.build_matrix(N1, prec=320)
    circ, eps80, beta80, dc = counting(M80)
    M140 = cg.build_matrix(N2, prec=320)
    loc = localise(M140, M80, dc)
    if loc is None:
        return 1
    co.big_disk(K=arcs)
    S_rho = co.circle_sup_rho(K=arcs)
    f = co.taylor_coeffs(S_rho, M=M_pts, Nc=N2)
    if f is None:
        return 1
    a4, c0, ev = a4_coeff(f, S, loc)
    rate(a4, c0, ev, S, circ["in"][2], loc["l4"])
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(cn.FAILS) + len(cg.FAILS) + len(co.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem CQ checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
