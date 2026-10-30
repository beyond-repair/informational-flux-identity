#!/usr/bin/env python3
r"""Theorem CR: the fifth eigenvalue l_5 of CE's transfer operator, a_5 != 0, and every further eigenvalue below 1e-8.

Setting (GASKET.md, Theorems CE-CG, CN-CQ): (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s) on H^2(D),
D = {|s - c| < R}, c = 83/50, R = 8/5, e_n = u^n, u = (s - c)/R; M_N = P_N Lt P_N, ||Lt (I - P_N)|| <= eps_N and
||Lt P_N|| <= beta_N (CG.1's formula, any N).  F0 = Lt^4 G (CN.1), x0 = sqrt13/2, c0 = 3 sigma_-3, a_j = c0 (Pi_j F0)(x0).
CQ left every eigenvalue other than r, l_2, l_3, l_4 inside |z| < 5e-5 (Pellet at 1e-5 failed), with float guesses
-1.744e-7 and a complex pair near 1.55e-7 for the next ones.  Those floats were float noise: M_N is extremely non-normal
(||(z - M_N)^-1|| ~ 1e20 on a circle of radius 2e-8 about the true l_5), and double-precision eigenvalues below ~1e-6 are
meaningless.  mpmath at 60/90 digits (M_80/M_140) gives l_5 ~ -4.5501411880e-8 and |l_6| ~ 2.07e-15; that is only the
seed of the certified computation below.

Theorem CR (checked here; exit 1 on any failure).  N = 300 (eps_300 <= 1.141e-63), 512-bit arb.
  CR.1 [certified; arb Krawczyk] Right/left eigenpairs (mu_i, v_i), (mu_i, w_i) of M_N, i = 1..5, enclosed; the five
       enclosures are disjoint, real, and the rank-one Riesz projectors P_i = v_i w_i^T / (w_i . v_i) give the exact
       M-invariant splitting C^(N+1) = span(v_1..v_5) + ran Q, Q = I - sum P_i.
  CR.2 [certified; arb] Spectral remainder: with A = M Q, T_j >= ||A^j Q||_F (j < m) and ||A^m||_F <= alpha_m,
           ||(z - M)^-1 Q|| <= H(|z|) := sum_(j<m) T_j |z|^-(j+1) / (1 - alpha_m |z|^-m)      for |z| >= rho_in = 1e-8,
       so every eigenvalue of M_N other than mu_1..mu_5 has |z| < rho_in, the mu_i are algebraically simple, and on any
       circle G with min |z| >= rho_in, ||(z - M)^-1|| <= sum ||P_i|| / dist(mu_i, G) + H(min_G |z|)   (no cover needed).
  CR.3 [proved; CG.3's homotopy] With K_T = K (1 + beta_N / min|z|) + 1/min|z| >= ||(z - Lt P_N)^-1|| (CG.3's block
       formula): K_T eps_N < 1 on |z| = rho_in, on |z| = R_out > beta_N + eps_N >= ||Lt||, and on five small circles G_i
       about mu_i.  So Lt has exactly five eigenvalues in rho_in < |z|, one inside each G_i, each real and algebraically
       simple: r, l_2, l_3, l_4 (overlapping CG/CQ's enclosures) and l_5; every other eigenvalue has |mu| < 1e-8.
       l_5 is enclosed by its G_5 of radius 1e-21.
  CR.4 [proved; certified] CO.3/CO.4's exact projection formula for Lt P_N and the second resolvent identity on a circle
       G_5' = {|z - mu_5| = 1e-8}:  |a_5 - a_5^N| <= c0 ||ev_x0|| rad K K^N eps_N S <= 2.83e-16,  a_5^N = -2.13235e-11,
       so a_5 = -2.1323e-11 +- 7.9e-16 != 0 (negative).  Taylor coefficients f_0..f_N as in CO.2 (4096 points, aliasing
       <= 1.5e-37).  At N = 200 (eps_200 <= 1.03e-42) the same bound is 2.55e5 and does not decide the sign: K K^N ~ 1e57
       on G_5', so it is eps_N that has to absorb it.
  CR.5 [proved from CR.3, CR.4, CN.1] s_(k+4) = sum_(i<=5) a_i l_i^k + e''_k for every k >= 0, |e''_k| <= eps_0'' rho_in^k,
       eps_0'' = c0 ||ev_x0|| rho_in K_in S.  Hence CQ's remainder e'_k = a_5 l_5^k + e''_k satisfies e'_k / l_5^k -> a_5
       and sign(e'_k) = (-1)^(k+1) for every k >= k_0 = 48 (eps_0'' <= 6.79e20, rho_in/|l_5| = 0.21977): its exact rate
       is |l_5| (not 5e-5).
Kept limits: no closed form of l_5 or a_5; k_0 = 48 is what these bounds give, not optimal; l_6 (mpmath ~ -2.07e-15)
is not certified, only |mu| < 1e-8 for everything past l_5; the rates of CP/CQ are unchanged; C_G, circle eigenvalues |z| = rho(L_k), the inf R_n rate stay open.
Assumptions used as stated (not re-run): CE.1-CE.3, CG.1's tail formula (its constants re-run), CQ.3's l_4 enclosure, CG.3's homotopy argument, CN.1's box argument
(box checks and S re-run), CO.1's continuation (re-run).  Textbook: Krawczyk/Brouwer, Riesz projections, Neumann series,
second resolvent identity, Cauchy/trapezoid aliasing.  Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW
enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction as Fr

from flint import arb, arb_mat, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_coefficient as cn     # noqa: E402
import gasket_fractional_profile_level_spectrum as cg        # noqa: E402
import gasket_fractional_profile_level_third as co           # noqa: E402

FAILS: list[str] = []
RESULT: dict = {}
N = int(os.environ.get("CR_N", "300"))   # truncation; at N = 200 CR.1-CR.3 pass but CR.4 does not (see the docstring)
PREC = 512
M_POW = 7                  # power m in the Neumann tail of CR.2
RHO_IN = Fr(1, 10 ** 8)
R_OUT = Fr(6)
M_TAYLOR = 4096            # trapezoid points for f_0..f_N
# seeds (mid values only; the enclosures below are certified independently)
SEEDS = [cn.R_MID, cn.L2_MID, cn.L3_MID, "-0.00021644655283736052", "-4.5501411880086498e-8"]
NAMES = ["r", "l2", "l3", "l4", "l5"]
# small circles G_i about mu_i (radius) and the earlier certified enclosures for consistency (CG.3 / CQ.3)
G_RAD = [Fr(1, 10 ** 20), Fr(1, 10 ** 20), Fr(1, 10 ** 20), Fr(1, 10 ** 20), Fr(1, 10 ** 21)]
EARLIER = {"r": (cn.R_MID, "1e-14"), "l2": (cn.L2_MID, cn.L2_RAD), "l3": (cn.L3_MID, cn.L3_RAD),
           "l4": ("-0.00021644655283736", "1.6e-18")}
RAD5_PROJ = Fr(1, 10 ** 8)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


_a = cn._a


def _mid(A):
    return arb_mat([[arb(A[i, j].mid()) for j in range(A.ncols())] for i in range(A.nrows())])


def _fro(A):
    s = arb(0)
    for i in range(A.nrows()):
        for j in range(A.ncols()):
            s += A[i, j] ** 2
    return arb(s.sqrt().upper())


def _col(A, j=0):
    return [A[i, j] for i in range(A.nrows())]


# ----------------------------------------------------------------- CR.1 Krawczyk eigenpairs
def eigpair(A, seed: str, n: int):
    """Certified (v, mu) with A v = mu v, v_k = 1 (Krawczyk on the n unknowns v_j (j != k), mu).  Seed by inverse iteration
    in mid arithmetic at PREC bits."""
    ctx.prec = PREC
    Am = _mid(A)
    Id = arb_mat(n, n)
    for i in range(n):
        Id[i, i] = 1
    lam = arb(seed)
    x = arb_mat([[arb(1)] for _ in range(n)])
    k = 0
    for _ in range(8):
        y = _mid((Am - lam * Id).solve(x))
        k = max(range(n), key=lambda i: abs(y[i, 0]).mid())
        y = _mid(y * (1 / y[k, 0]))
        lam = arb(((Am * y)[k, 0]).mid())
        x = y

    def J(xv, lm):
        Jm = A - lm * Id
        for r_ in range(n):
            Jm[r_, k] = -xv[r_, 0]
        return Jm
    # two Newton steps on mids
    for _ in range(2):
        Fx = Am * x - lam * x
        d = _mid(J(x, lam)).solve(Fx)
        x = arb_mat([[arb((x[j, 0] - d[j, 0]).mid()) if j != k else arb(1)] for j in range(n)])
        lam = arb((lam - d[k, 0]).mid())
    Y = _mid(_mid(J(x, lam)).inv())
    YF = Y * (A * x - lam * x)
    for e in (arb("1e-130"), arb("1e-100"), arb("1e-80"), arb("1e-60")):
        X = arb_mat([[x[j, 0] + arb(0, e) if j != k else arb(1)] for j in range(n)])
        Lm = lam + arb(0, e)
        D = Id - Y * J(X, Lm)
        dx = arb_mat([[arb(0, e)] for _ in range(n)])
        Kx = D * dx
        if all(bool(abs(-YF[j, 0] + Kx[j, 0]) < e) for j in range(n)):
            v = arb_mat([[x[j, 0] + arb(0, e) if j != k else arb(1)] for j in range(n)])
            return v, lam + arb(0, e), e
    return None, None, None


def spectral_split(M):
    n = N + 1
    ctx.prec = PREC
    MT = M.transpose()
    pairs = []
    for name, sd in zip(NAMES, SEEDS):
        v, mu, e1 = eigpair(M, sd, n)
        w, muL, e2 = eigpair(MT, sd, n)
        ok = v is not None and w is not None and bool(mu.overlaps(muL))
        if not ok:
            check(False, f"CR.1: Krawczyk enclosure for {name} failed")
            return None
        wv = (w.transpose() * v)[0, 0]
        nv, nw = _fro(v), _fro(w)
        Pn = arb((nv * nw / abs(wv)).upper())
        pairs.append(dict(name=name, v=v, w=w, mu=mu, wv=wv, P=Pn, e=(e1, e2)))
        check(not wv.contains(0), f"CR.1 (arb Krawczyk, {PREC} bits, N = {N}): {name}: mu_N = {mu.str(25)} (right/left "
                                  f"radii {float(e1):.0e}/{float(e2):.0e}), w.v != 0, ||P|| = ||v|| ||w||/|w.v| = {Pn.str(5)}")
    disj = all(not pairs[i]["mu"].overlaps(pairs[j]["mu"]) for i in range(5) for j in range(i + 1, 5))
    check(disj, "CR.1: the five eigenvalue enclosures are disjoint and real (so w_i.v_j = 0 for i != j exactly, the P_i are "
                "complementary rank-one projectors commuting with M)")
    Id = arb_mat(n, n)
    for i in range(n):
        Id[i, i] = 1
    Q = Id
    for p in pairs:
        Q = Q - p["v"] * p["w"].transpose() * (1 / p["wv"])
    A = M * Q
    T = []
    X = Q
    for j in range(M_POW):
        T.append(_fro(X))
        X = A * X
    alpha = _fro(X)                       # X = A^m Q = A^m (A Q = A)
    check(all(t.is_finite() for t in T) and alpha.is_finite(),
          f"CR.2 (arb): ||A^j Q||_F <= " + ", ".join(t.str(3) for t in T) + f" (j = 0..{M_POW - 1}), ||A^{M_POW}||_F <= "
          f"{alpha.str(3)}, A = M Q")
    return pairs, T, alpha


def Hbound(T, alpha, t: arb):
    q = alpha / t ** M_POW
    if not bool(q < 1):
        return None
    s = arb(0)
    for j, tj in enumerate(T):
        s += tj / t ** (j + 1)
    return arb((s / (1 - q)).upper())


def circle_K(pairs, T, alpha, cen: arb, rad: arb):
    """Upper bound for sup_G ||(z - M)^-1||, G = {|z - cen| = rad}, and min_G |z|."""
    absmin = abs(cen) - rad if bool(abs(cen) > rad) else rad - abs(cen)
    absmin = arb(absmin.lower())
    H = Hbound(T, alpha, absmin)          # valid on |z| >= absmin whenever the Neumann tail converges
    if H is None:
        return None, absmin
    K = H
    for p in pairs:
        d = abs(abs(p["mu"] - cen) - rad)
        d = arb(d.lower()) - p["mu"].rad()
        if not bool(d > 0):
            return None, absmin
        K += p["P"] / d
    return arb(K.upper()), absmin


def certify(pairs, T, alpha, eps, beta):
    out = {}
    circles = [("in", arb(0), _a(RHO_IN)), ("out", arb(0), _a(R_OUT))]
    for p, rad in zip(pairs, G_RAD):
        circles.append((p["name"], arb(p["mu"].mid()), _a(rad)))
    circles.append(("l5proj", arb(pairs[4]["mu"].mid()), _a(RAD5_PROJ)))
    for name, cen, rad in circles:
        KM, absmin = circle_K(pairs, T, alpha, cen, rad)
        if KM is None:
            check(False, f"CR.3: circle {name} not certified")
            out[name] = None
            continue
        KT = KM * (1 + beta / absmin) + 1 / absmin
        ok = bool(KT * eps < 1)
        KL = KT / (1 - KT * eps) if ok else None
        check(ok, f"CR.3: |z - {cen.str(12)}| = {rad.str(3)}: sup ||(z - M)^-1|| <= {KM.str(4)}, sup ||(z - Lt P_N)^-1|| <= "
                  f"K_T = {KT.str(4)}, K_T eps_{N} = {(KT * eps).str(3)} < 1" + (f", sup ||(z - Lt)^-1|| <= {KL.str(4)}" if ok else ""))
        out[name] = (KM, KT, KL, cen, rad) if ok else None
    # layout: every mu_i inside its G_i (G_i disjoint from each other and from |z| <= rho_in), all inside |z| < R_OUT
    lay = bool(beta + eps < _a(R_OUT))
    for p, rad in zip(pairs, G_RAD):
        lay = lay and bool(abs(p["mu"] - arb(p["mu"].mid())) < _a(rad)) and bool(abs(p["mu"]) - _a(rad) > _a(RHO_IN))
    for i in range(5):
        for j in range(i + 1, 5):
            lay = lay and bool(abs(pairs[i]["mu"] - pairs[j]["mu"]) > _a(G_RAD[i]) + _a(G_RAD[j]))
    p5 = pairs[4]
    lay = lay and bool(abs(p5["mu"]) - _a(RAD5_PROJ) > _a(RHO_IN)) and bool(abs(pairs[3]["mu"] - p5["mu"]) > _a(RAD5_PROJ))
    allok = lay and all(out[k] is not None for k in out)
    check(allok, f"CR.3: beta_N + eps_N = {(beta + eps).str(5)} < {R_OUT} (so ||Lt|| < R_OUT); G_i hold mu_i and are disjoint and "
                 f"outside |z| <= {RHO_IN}; so is G_5' = |z - mu_5| = {RAD5_PROJ} (only mu_5 inside).  CG.3's homotopy z - Lt P_N - tau Lt (I - P_N): Lt has exactly five eigenvalues in "
                 f"{RHO_IN} < |z|, one inside each G_i, real (G_i symmetric, Lt commutes with f -> conj f(conj s)) and "
                 f"algebraically simple; every other eigenvalue has |mu| < {RHO_IN}")
    # identification with earlier enclosures
    ident = True
    for p, rad in zip(pairs[:4], G_RAD[:4]):
        m, r_ = EARLIER[p["name"]]
        enc = arb(p["mu"].mid()) + arb(0, (_a(rad) + p["mu"].rad()).upper())
        ident = ident and bool(enc.overlaps(arb(m) + arb(0, arb(r_))))
    l5 = arb(pairs[4]["mu"].mid()) + arb(0, (_a(G_RAD[4]) + pairs[4]["mu"].rad()).upper())
    r = arb(cn.R_MID)
    l4 = pairs[3]["mu"]
    check(ident and bool(l5 < 0), f"CR.3: the eigenvalues inside G_1..G_4 overlap CG.3's r, l_2, l_3 and CQ.3's l_4 (so they are those); "
                                  f"l_5 = {l5.str(20)} (+- {float(G_RAD[4]):.0e}), negative; |l_5|/r = {(abs(l5) / r).str(10)}, "
                                  f"|l_5|/|l_4| = {(abs(l5) / abs(l4)).str(10)}")
    return out, l5, allok and ident


# ----------------------------------------------------------------- CR.4, CR.5
def a5_coeff(f, S, pairs, circ, eps):
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    c0 = 3 * ce.sigma_closed(arb(-3), 0)
    x0 = arb(13).sqrt() / 2
    u0 = (x0 - _a(cn.C_CEN)) / _a(cn.R_RAD)
    ev = 1 / (1 - u0 * u0).sqrt()
    p = pairs[4]
    v, w, mu, wv = _col(p["v"]), _col(p["w"]), p["mu"], p["wv"]
    wf = sum((w[i] * f[i] for i in range(N + 1)), arb(0))
    phi = co.lt_eval(v, x0) / mu
    aN = c0 * wf / wv * phi
    KM, KT, KL, cen, rad = circ["l5proj"]
    err = c0 * ev * rad * KL * KT * eps * S
    a5 = aN + arb(0, err.upper())
    print(f"   a_5^N = {aN.str(15)}  (c0 = {c0.str(10)}, (w.f)/(w.v) = {(wf / wv).str(8)}, phi(x0) = {phi.str(8)})")
    check(not a5.contains(0), f"CR.4: on |z - mu_5| = {RAD5_PROJ}: |a_5 - a_5^N| <= c0 ||ev|| rad K K^N eps_{N} S = {err.str(3)}; "
                              f"a_5 = {a5.str(8)} != 0 ({'negative' if bool(a5 < 0) else 'positive'})")
    return a5, c0, ev


def rate(a5, c0, ev, S, circ, l5):
    ctx.prec = PREC
    KM, KT, KL, cen, rad = circ["in"]
    rho = _a(RHO_IN)
    eps0 = c0 * ev * rho * KL * S
    th = rho / abs(l5)
    k0 = None
    for k in range(400):
        if bool(eps0 * th ** k < abs(a5)):
            k0 = k
            break
    sgn = "(-1)^k" if bool(a5 > 0) else "(-1)^(k+1)"
    ok = k0 is not None and bool(th < 1)
    check(ok, f"CR.5: s_(k+4) = sum_(i<=5) a_i l_i^k + e''_k, |e''_k| <= eps_0'' rho_in^k, eps_0'' = c0 ||ev|| rho_in K_in S = "
              f"{eps0.str(4)}; rho_in/|l_5| = {th.str(6)} < 1 and eps_0'' (rho_in/|l_5|)^k < |a_5| for k >= k_0 = {k0}.  So CQ's "
              f"remainder e'_k = a_5 l_5^k + e''_k has e'_k / l_5^k -> a_5 and sign(e'_k) = {sgn} for every k >= {k0}: exact "
              f"rate |l_5| (not 5e-5)")
    return eps0, k0


def main(M_pts=M_TAYLOR) -> int:
    t0 = time.time()
    print("Theorem CR: the fifth eigenvalue l_5, a_5 != 0, every further eigenvalue below 1e-8. Not thrust.", flush=True)
    cn.tail_constants()
    cn.box_checks()
    S = cn.circle_sup(1024)
    dc = cg.disk_constants()
    eps, beta = cg.tail_constants(dc, N)
    M = cg.build_matrix(N, prec=PREC)
    sp = spectral_split(M)
    if sp is None:
        return 1
    pairs, T, alpha = sp
    circ, l5, ok = certify(pairs, T, alpha, eps, beta)
    print(f"   [{time.time() - t0:.0f}s]", flush=True)
    if not ok:
        return 1
    co.big_disk(K=4096)
    S_rho = co.circle_sup_rho(K=4096)
    f = co.taylor_coeffs(S_rho, M=M_pts, Nc=N)
    if f is None:
        return 1
    a5, c0, ev = a5_coeff(f, S, pairs, circ, eps)
    eps0, k0 = rate(a5, c0, ev, S, circ, l5)
    RESULT.update(l5=l5, a5=a5, eps0=eps0, k0=k0, eps=eps)
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(cn.FAILS) + len(cg.FAILS) + len(co.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem CR checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
