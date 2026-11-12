#!/usr/bin/env python3
r"""Theorem DR: DQ.5's threshold made effective.  Certified upper bounds for the Hilbert diameter D_n = diam Lt^n K_(2,6)
of the second-order cone K_(2,6) = {f > 0: |f'| <= 2 f, |f''| <= 6 f} (DQ.1), so tanh(D_n/4)^(1/n) < rho(L_1)/r (DP's cap)
for every n >= 5; with DQ.5 (n = 1, 2, 3 above the cap) the threshold n0 of DQ.4 is 4 or 5.  (DQ's kept failure DQ.5:
"n0 is not effective; no certified D_n on K_(2,6)".)

Setting (GASKET.md, Theorems CF, CI-CM, DP, DQ): s in X = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
wt1 = (s+1/2)/(2s(7/2-s)), wt2 = (s-1/2)/(2s(s+7/2)), sigma_(1,2) = sqrt(15/4 -+ s).  Lt^n f(s) = sum_w W_w(s) f(y_w(s)) over
the 2^n words (W_w products of weights, y_w compositions of the sigma's; the images y_w(X) are disjoint).  K = K_(2,6), its
Hilbert metric Theta (CF's definition) through the functionals l in Lam = {2 d_s +- d'_s, 6 d_s +- d''_s : s in X} (DQ's set;
d_s = ((2d_s + d'_s) + (2d_s - d'_s))/4 is redundant), Phi_l(f) = l(Lt^n f).

Theorem DR (checked here; exit 1 on any failure):
  DR.1 [proved; identities in sympy] Chain lemma.  For f in K and y, z in X, d = z - y:
       m f(y) + mu f'(y) = m f(z) + (mu - m d) f'(z) + int_y^z (m (x - y) - mu) f''(x) dx, and |f''| <= 6 f, f(x) <= f(z) e^(2|z-x|)
       give m f(y) + mu f'(y) >= m f(z) + (mu - m d) f'(z) - 6 (|m| I1(|d|) + |mu| I0(|d|)) f(z), I0 = (e^(2d) - 1)/2,
       I1 = (e^(2d) - 1 - 2d)/4; and m f(y) + mu f'(y) >= c f(y) >= c e^(-+2|d|) f(z) (c = m - 2|mu|, sign of c).  A finite
       point-jet functional Psi = sum_k a_k d_(y_k) + b_k d'_(y_k) is >= 0 on K if the chain (split the carried (m, mu) into a
       Taylor part and a collapsed part at every y_k, move both to y_(k+1), add (a, b) there) ends with m - 2|mu| >= 0.
  DR.2 [proved] Cell lemma.  For s = s_c + tau, |tau| <= Delta, and any kappa with |kappa| Delta < 1, the rescaled functional
       l_s(Lt^n f)/(1 + kappa tau) lies between G_tau -+ E (as functionals of f in K), with G_tau = sum_w (a0 + tau a1) d_(y_c) +
       (b0 + tau b1) d'_(y_c) affine in tau (y_c = y_w(s_c); a0, a1, b0, b1 from arb jets of W_w, y_w at s_c) and E = sum_w
       E_w d_(y_c) (E_w from arb enclosures of the derivatives of W_w, y_w up to order 4 on the cell; it collects the Taylor
       remainders, the f'' terms bounded by 6 f, and the quadratic part of the rescaling).  Since R(l, l') = sup_K Phi_l/Phi_l'
       is convex in affine families, two chains at tau = +-Delta bound the whole cell; with a fixed reference l0 = d_(s0),
       s0 = 9/5: R(l, l') <= R(l, l0) R(l0, l'), so D_n <= 2 max_cells log(t1 t2), t1 >= R(l_s, l0), t2 >= R(l0, l_s).
  DR.3 [certified; arb, chains proposed by a float LP] D_n <= Dbar_n for n = 3, ..., 8 (n = 9 with DR_FULL=1 and DR_NS),
       200 cells x 4 functional types, every chain verified in arb: theta_n = tanh(Dbar_n/4) and per-level theta_n^(1/n) are
       printed; for n = 5, 6, 7, 8 the per-level rate is < 0.3170757 <= rho(L_1)/r (DP.3's cap), for n = 3, 4 above it.
  DR.4 [proved from DR.3, CF.1] For every n >= 5, tanh(D_n/4)^(1/n) <= max_j (Dbar_j/4)^(1/j) < 0.3170757 (n = 5q + j,
       j in {5..9}, D_n <= theta_5^q D_j; D_9 <= theta_3 Dbar_6 unless n = 9 is certified directly), so n0 <= 5 and, with
       DQ.5, n0 in {4, 5}.  Explicit envelope: tanh(D_n/4) <= C rho^n (n >= 5), rho = theta_8^(1/8), C = max_(j=5..12) (Dbar_j/4)/rho^j.
  OPEN: n = 4 (DR.3's chain gives only per-level <= 0.3433; floats: pairwise chain 0.04170 and a float primal 0.03638 at the
       worst pair vs the cap's 4 artanh(q^4) = 0.040432: not certified either way).
Kept limits: the chain is a relaxation (ell_0 reference costs up to a factor 2 in D_n); chains are float LP proposals, only
their arb replay is claimed.  Textbook facts: Taylor's formula with integral remainder, Birkhoff's contraction (CF.1).
No Psi_info, Stage 2 still blocked.  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a
flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time

import numpy as np
import sympy as sp
from flint import arb, arb_series, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_essential as ci          # noqa: E402  (rho(L_1) and r as certified in CI.1/CG.3)

FAILS: list[str] = []
RESULT: dict = {}
ORD = 4
N_CELLS = int(os.environ.get("DR_CELLS", "200"))
S0 = 1.8                                  # reference functional l0 = d_(s0)
CAP = "0.3170757"                         # DP.3 / DQ.3: rho(L_1)/r >= 0.317075713
DEFAULT_NS = "3,4,5,6,7,8" if os.environ.get("DR_FULL") == "1" else "3,5,6"
NS = [int(x) for x in os.environ.get("DR_NS", DEFAULT_NS).split(",")]


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- DR.1 / DR.2 identities (sympy)
def symbolic_checks():
    x, y, z, m, mu, d, tau, kap = sp.symbols("x y z m mu d tau kappa", real=True)
    c = sp.symbols("c0:7", real=True)
    f = sum(c[j] * x ** j for j in range(7))
    fy, fz = f.subs(x, y), f.subs(x, z)
    fpy, fpz = sp.diff(f, x).subs(x, y), sp.diff(f, x).subs(x, z)
    rem = sp.integrate((m * (x - y) - mu) * sp.diff(f, x, 2), (x, y, z))
    ok1 = sp.expand(m * fy + mu * fpy - (m * fz + (mu - m * (z - y)) * fpz + rem)) == 0
    u = sp.symbols("u", positive=True)
    I0 = sp.integrate(sp.exp(2 * (d - u)), (u, 0, d))
    I1 = sp.integrate(u * sp.exp(2 * (d - u)), (u, 0, d))
    ok2 = sp.simplify(I0 - (sp.exp(2 * d) - 1) / 2) == 0 and sp.simplify(I1 - (sp.exp(2 * d) - 1 - 2 * d) / 4) == 0
    check(ok1 and ok2, "DR.1 (sympy): m f(y) + mu f'(y) = m f(z) + (mu - m d) f'(z) + int_y^z (m(x-y) - mu) f'' dx (generic "
                       "degree-6 f); int_0^d e^(2(d-u)) du = (e^(2d)-1)/2, int_0^d u e^(2(d-u)) du = (e^(2d)-1-2d)/4")
    # second derivative of W f(Y) (the type-2 functional) and the derivatives of A, B used by the cell bounds
    s = sp.symbols("s", real=True)
    W, Y, F = sp.Function("W")(s), sp.Function("Y")(s), sp.Function("f")
    t = sp.symbols("t")
    expr = W * F(Y)
    fp = lambda k: sp.Subs(sp.diff(F(t), t, k), t, Y).doit() if k else F(Y)
    d2 = (sp.diff(W, s, 2) * fp(0) + (2 * sp.diff(W, s) * sp.diff(Y, s) + W * sp.diff(Y, s, 2)) * fp(1)
          + W * sp.diff(Y, s) ** 2 * fp(2))
    ok3 = sp.simplify(sp.diff(expr, s, 2).doit() - d2) == 0
    B2 = 2 * sp.diff(W, s) * sp.diff(Y, s) + W * sp.diff(Y, s, 2)
    Wd = lambda k: sp.diff(W, s, k)
    Yd = lambda k: sp.diff(Y, s, k)
    ok4 = (sp.simplify(sp.diff(B2, s) - (2 * Wd(2) * Yd(1) + 3 * Wd(1) * Yd(2) + W * Yd(3))) == 0
           and sp.simplify(sp.diff(B2, s, 2) - (2 * Wd(3) * Yd(1) + 5 * Wd(2) * Yd(2) + 4 * Wd(1) * Yd(3) + W * Yd(4))) == 0
           and sp.simplify(sp.diff(W * Yd(1), s, 2) - (Wd(2) * Yd(1) + 2 * Wd(1) * Yd(2) + W * Yd(3))) == 0)
    G0, G1 = sp.symbols("G0 G1")
    ok5 = sp.simplify((G0 + tau * G1) / (1 + kap * tau) - (G0 + tau * (G1 - kap * G0))
                      + kap * tau ** 2 * (G1 - kap * G0) / (1 + kap * tau)) == 0
    check(ok3 and ok4 and ok5, "DR.2 (sympy): (W f(Y))'' = W'' f + (2W'Y' + W Y'') f' + W Y'^2 f''; derivatives of B = W Y' and "
                               "B = 2W'Y' + W Y'' to order 2; (G0 + tau G1)/(1 + kappa tau) = G0 + tau (G1 - kappa G0) - "
                               "kappa tau^2 (G1 - kappa G0)/(1 + kappa tau)")


# ----------------------------------------------------------------- word jets (arb)
def _wt(i, x):
    h, sv = arb(1) / 2, arb(7) / 2
    return (x + h) / (2 * x * (sv - x)) if i == 0 else (x - h) / (2 * x * (x + sv))


def _sig(i, x):
    q = arb(15) / 4
    return (q - x).sqrt() if i == 0 else (q + x).sqrt()


def _ders(ser):
    c = ser.coeffs() + [arb(0)] * (ORD + 1)
    return [c[j] * math.factorial(j) for j in range(ORD + 1)]


def word_jets(s, n):
    """derivatives in s (orders 0..4) of W_w and y_w for the 2^n words; s an arb (point or cell ball)."""
    out = []

    def rec(W, y, k):
        if k == n:
            out.append((_ders(W), _ders(y)))
            return
        for i in (0, 1):
            rec(W * _wt(i, y), _sig(i, y), k + 1)

    rec(arb_series([1], prec=ORD + 1), arb_series([s, 1], prec=ORD + 1), 0)
    return out


def _AB(typ, W, y):
    sg = 1 if typ[1] == "+" else -1
    if typ[0] == "1":       # 2 d_s +- d'_s
        A = [2 * W[0] + sg * W[1], 2 * W[1] + sg * W[2], 2 * W[2] + sg * W[3]]
        B = [sg * W[0] * y[1], sg * (W[1] * y[1] + W[0] * y[2]), sg * (W[2] * y[1] + 2 * W[1] * y[2] + W[0] * y[3])]
        C = arb(0)
    else:                   # 6 d_s +- d''_s
        A = [6 * W[0] + sg * W[2], 6 * W[1] + sg * W[3], 6 * W[2] + sg * W[4]]
        B = [sg * (2 * W[1] * y[1] + W[0] * y[2]), sg * (2 * W[2] * y[1] + 3 * W[1] * y[2] + W[0] * y[3]),
             sg * (2 * W[3] * y[1] + 5 * W[2] * y[2] + 4 * W[1] * y[3] + W[0] * y[4])]
        C = sg * W[0] * y[1] * y[1]
    return A, B, C


def cell_functional(typ, cen, ball, delta):
    """DR.2: y_c, a0, a1, b0, b1, E for the cell |tau| <= delta (after the 1/(1 + kappa tau) rescaling)."""
    D = arb(delta)
    ys, a0s, a1s, b0s, b1s, Es = [], [], [], [], [], []
    for (Wc, yc), (Wb, yb) in zip(cen, ball):
        A, B, _ = _AB(typ, Wc, yc)
        Ab, Bb, Cb = _AB(typ, Wb, yb)
        ys.append(yc[0])
        a0s.append(A[0])
        a1s.append(A[1])
        b0s.append(B[0])
        b1s.append(B[1] + A[0] * yc[1])
        sA, sA1, sA2 = abs(Ab[0]).upper(), abs(Ab[1]).upper(), abs(Ab[2]).upper()
        sB, sB2, sC = abs(Bb[0]).upper(), abs(Bb[2]).upper(), abs(Cb).upper()
        sy1, sy2 = abs(yb[1]).upper(), abs(yb[2]).upper()
        eb = arb(sy1) * D
        ex = (2 * eb).exp()
        E = (D * D / 2 * sA2 + 2 * D * D * (arb(sA1) * sy1 + abs(A[0]) * sy2 / 2 + arb(sB2) / 2)
             + sA * 3 * eb * eb * ex + sB * 6 * eb * ex + sC * 6 * ex)
        Es.append(arb(E.upper()))
    kap = arb(float((sum(a1s) / sum(a0s)).mid()))
    den = 1 - abs(kap) * D
    assert den > 0
    a1n = [p - kap * q for p, q in zip(a1s, a0s)]
    b1n = [p - kap * q for p, q in zip(b1s, b0s)]
    Es = [arb(((e + abs(kap) * D * D * (abs(p) + 2 * abs(q))) / den).upper()) for e, p, q in zip(Es, a1n, b1n)]
    return ys, a0s, a1n, b0s, b1n, Es


# ----------------------------------------------------------------- DR.1 chain: float LP proposal, arb replay
def lp_split(y, aP, bP, aQ, bQ, eta=1e-8):
    """float proposal only: min t such that the DR.1 chain certifies t Q - P >= 0 (sorted float data)."""
    import scipy.sparse as sps
    from scipy.optimize import linprog
    N = len(y)
    d = np.abs(np.diff(y))
    nv = 9
    idx = lambda k, j: 1 + nv * k + j
    nV = 1 + nv * (N - 1) + 1
    R, Cc, V, rhs, ER, EC, EV, EB = [], [], [], [], [], [], [], []
    r, e = [0], [0]

    def ub(ent, b):
        for c, v in ent:
            R.append(r[0]); Cc.append(c); V.append(v)
        rhs.append(b); r[0] += 1

    def eq(ent, b):
        for c, v in ent:
            ER.append(e[0]); EC.append(c); EV.append(v)
        EB.append(b); e[0] += 1

    Mk = lambda k: ([(0, aQ[k])] + ([(idx(k - 1, 0), 1.0)] if k else []), aP[k])
    Uk = lambda k: ([(0, bQ[k])] + ([(idx(k - 1, 1), 1.0)] if k else []), bP[k])
    for k in range(N - 1):
        dk = d[k]
        I0, I1 = np.expm1(2 * dk) / 2, (np.expm1(2 * dk) - 2 * dk) / 4
        mn, mu, ma, mua, maA, muaA, mubA, p, q = [idx(k, j) for j in range(9)]
        ub([(mn, 1), (ma, -1), (maA, 6 * I1), (muaA, 6 * I0), (p, -np.exp(-2 * dk)), (q, np.exp(2 * dk))], 0)
        Me, Mc = Mk(k)
        ub([(p, 1), (q, -1), (ma, 1), (mubA, 2)] + [(c, -v) for c, v in Me], Mc)
        ub([(ma, 1), (maA, -1)], 0); ub([(ma, -1), (maA, -1)], 0)
        ub([(mua, 1), (muaA, -1)], 0); ub([(mua, -1), (muaA, -1)], 0)
        Ue, Uc = Uk(k)
        ub(Ue + [(mua, -1), (mubA, -1)], -Uc)
        ub([(c, -v) for c, v in Ue] + [(mua, 1), (mubA, -1)], Uc)
        eq([(mu, 1), (mua, -1), (ma, np.sign(y[k + 1] - y[k]) * dk)], 0)
    UA = nV - 1
    Me, Mc = Mk(N - 1)
    Ue, Uc = Uk(N - 1)
    ub([(c, -v) for c, v in Me] + [(UA, 2)], Mc - eta)
    ub(Ue + [(UA, -1)], -Uc)
    ub([(c, -v) for c, v in Ue] + [(UA, -1)], Uc)
    bounds = [(None, None)] * nV
    for k in range(N - 1):
        for j in (4, 5, 6, 7, 8):
            bounds[idx(k, j)] = (0, None)
    bounds[UA] = (0, None)
    cobj = np.zeros(nV)
    cobj[0] = 1
    Aub = sps.csr_matrix((V, (R, Cc)), shape=(r[0], nV))
    Aeq = sps.csr_matrix((EV, (ER, EC)), shape=(e[0], nV))
    for opts in (dict(primal_feasibility_tolerance=1e-10, dual_feasibility_tolerance=1e-10), {},
                 dict(primal_feasibility_tolerance=1e-9, dual_feasibility_tolerance=1e-9, presolve=False)):
        res = linprog(cobj, A_ub=Aub, b_ub=rhs, A_eq=Aeq, b_eq=EB, bounds=bounds, method="highs", options=opts)
        if res.status == 0:
            break
    else:
        return None
    xs = res.x
    return (xs[0], np.array([xs[idx(k, 2)] for k in range(N - 1)]), np.array([xs[idx(k, 3)] for k in range(N - 1)]))


def verify_chain(ya, alpha, beta, ma, mua):
    """DR.1 in arb: lower bound for m - 2|mu| at the end of the chain (Psi >= 0 on K if it is >= 0)."""
    m, mu = arb(0), arb(0)
    six = arb(6)
    for k in range(len(ya) - 1):
        M, U = m + alpha[k], mu + beta[k]
        dd = ya[k + 1] - ya[k]
        da = arb(abs(dd).upper())
        e2 = (2 * da).exp()
        I0, I1 = (e2 - 1) / 2, (e2 - 1 - 2 * da) / 4
        a, u = arb(float(ma[k])), arb(float(mua[k]))
        cost = six * (abs(a) * I1 + abs(u) * I0)
        c = ((M - a) - 2 * abs(U - u)).lower()
        col = arb(c) / e2 if c >= 0 else arb(c) * e2
        m = arb((a - cost + col).lower())
        mu = u - a * dd
    M, U = m + alpha[-1], mu + beta[-1]
    return (M - 2 * abs(U)).lower()


def certify_n(n, ncell=N_CELLS, only=None):
    ctx.prec = 96
    t0 = time.time()
    AX, BX = arb(5).sqrt() / 2, arb(5) / 2
    Z = word_jets(arb(S0), n)
    y0 = [w[1][0] for w in Z]
    w0 = [w[0][0] for w in Z]
    sQ = arb(1 / float(sum(w0).mid()))
    w0 = [w * sQ for w in w0]
    lo, hi = float(AX.mid()), float(BX.mid())
    edges = np.linspace(lo, hi, ncell + 1)
    cover_ok = True
    worst, wcell, nchains, bumps = arb(0), None, 0, 0
    fl = lambda L: np.array([float(v.mid()) for v in L])
    prev_hi = None
    for ci_ in range(ncell):
        if only is not None and ci_ not in only:
            continue
        sc = (edges[ci_] + edges[ci_ + 1]) / 2
        dl = (edges[ci_ + 1] - edges[ci_]) / 2 * (1 + 1e-9)
        cball = arb(sc, dl)
        if ci_ == 0:
            cover_ok &= bool(cball.lower() <= AX.lower())
        if ci_ == ncell - 1:
            cover_ok &= bool(cball.upper() >= BX.upper())
        if prev_hi is not None:
            cover_ok &= bool(cball.lower() <= prev_hi)
        prev_hi = cball.upper()
        cen, ball = word_jets(arb(sc), n), word_jets(cball, n)
        for typ in ("1+", "1-", "2+", "2-"):
            ys, a0, a1, b0, b1, E = cell_functional(typ, cen, ball, dl)
            S = arb(1 / float(sum(a0).mid()))
            a0, a1, b0, b1, E = ([v * S for v in L] for L in (a0, a1, b0, b1, E))
            yall = ys + y0
            o = np.argsort(fl(yall), kind="stable")
            yf = fl(yall)[o]
            ya = [yall[k] for k in o]
            nP, z0 = len(ys), [arb(0)] * len(y0)
            tv = []
            for direction in (1, 2):
                tmax = 0.0
                for tau in (-dl, dl):
                    T = arb(tau)
                    G = [p + T * q for p, q in zip(a0, a1)]
                    Gb = [p + T * q for p, q in zip(b0, b1)]
                    if direction == 1:       # t Phi_0 - (G + E) >= 0  =>  R(l_s, l0) <= t
                        aP, bP = [-(g + e) for g, e in zip(G, E)] + z0, [-g for g in Gb] + z0
                        aQ, bQ = [arb(0)] * nP + w0, [arb(0)] * (nP + len(y0))
                    else:                    # t (G - E) - Phi_0 >= 0  =>  R(l0, l_s) <= t
                        aP, bP = [arb(0)] * nP + [-w for w in w0], [arb(0)] * (nP + len(y0))
                        aQ, bQ = [g - e for g, e in zip(G, E)] + z0, Gb + z0
                    aP, bP, aQ, bQ = ([L[k] for k in o] for L in (aP, bP, aQ, bQ))
                    sol = lp_split(yf, fl(aP), fl(bP), fl(aQ), fl(bQ))
                    if sol is None:
                        check(False, f"DR.3 n={n}: LP proposal failed at cell {ci_} {typ}")
                        return None
                    t, ma, mua = sol
                    t = float(np.nextafter(t, np.inf))
                    for bump in range(7):
                        Tt = arb(t)
                        alpha = [Tt * q + p for q, p in zip(aQ, aP)]
                        beta = [Tt * q + p for q, p in zip(bQ, bP)]
                        if verify_chain(ya, alpha, beta, ma, mua) >= 0:
                            break
                        t, bumps = t * (1 + 1e-8 * 4 ** bump), bumps + 1
                    else:
                        check(False, f"DR.3 n={n}: arb chain replay failed at cell {ci_} {typ} dir {direction}")
                        return None
                    nchains += 1
                    tmax = max(tmax, t)
                tv.append(tmax)
            val = (arb(tv[0]) * arb(tv[1])).log()
            if val.upper() > worst.upper():
                worst, wcell = arb(val.upper()), (ci_, typ, round(float(sc), 5))
    check(cover_ok, f"DR.2 n={n}: the {ncell} cell balls cover X = [sqrt5/2, 5/2] (consecutive overlap, ends outside)")
    Dn = arb((2 * worst).upper())
    th = (Dn / 4).tanh()
    per = (th.log() / n).exp()
    print(f"     n={n}: {nchains} chains replayed in arb ({bumps} t-bumps), worst cell {wcell}, "
          f"{time.time() - t0:.0f}s", flush=True)
    return dict(D=Dn, theta=th, per=per, quarter_root=((Dn / 4).log() / n).exp(), worst=wcell)


def main() -> int:
    t0 = time.time()
    print("Theorem DR: certified diameters of Lt^n K_(2,6); tanh(D_n/4)^(1/n) below DP's cap for every n >= 5 "
          "(DQ.5's n0 made effective). Not thrust.", flush=True)
    symbolic_checks()
    ctx.prec = 96
    cap_lo = arb(ci.LAM_LO) / arb(ci.R_HI)
    check(bool(cap_lo > arb(CAP)), f"cap: rho(L_1)/r >= {ci.LAM_LO}/{ci.R_HI} > {CAP} (CI.1, CG.3 as stated in DQ.3)")
    cap = arb(CAP)
    res = {}
    for n in NS:
        r = certify_n(n)
        if r is None:
            break
        res[n] = r
        below = bool(r["per"] < cap)
        print(f"     DR.3 n={n}: D_n <= {r['D'].upper().str(8, radius=False)}, tanh(D_n/4) <= "
              f"{r['theta'].upper().str(8, radius=False)}, per level <= {r['per'].upper().str(7, radius=False)} "
              f"({'below' if below else 'above'} the cap {CAP})", flush=True)
        if n >= 5:
            check(below, f"DR.3 n={n}: tanh(D_n/4)^(1/n) <= {r['per'].upper().str(7, radius=False)} < {CAP} <= rho(L_1)/r")
    RESULT["D"] = res
    # DR.4: every n >= 5
    if all(j in res for j in (3, 5, 6, 7, 8)):
        Dq = {j: res[j]["D"] / 4 for j in res}
        if 9 not in res:
            Dq[9] = res[3]["theta"] * res[6]["D"] / 4          # D_9 <= tau(Lt^3) D_6 <= theta_3 D_6
        roots = {j: (Dq[j].log() / j).exp() for j in range(5, 10)}
        th5 = res[5]["theta"]
        unif = arb(max(max(v.upper() for v in roots.values()), ((th5.log()) / 5).exp().upper()))
        check(bool(unif < cap), f"DR.4: for every n >= 5 (n = 5q + j, j = 5..9, D_n <= theta_5^q D_j"
              f"{'' if 9 in res else ', D_9 <= theta_3 Dbar_6'}): tanh(D_n/4)^(1/n) <= {unif.upper().str(7, radius=False)} "
              f"< {CAP}; with DQ.5 (n = 1, 2, 3 above the cap) n0 in {{4, 5}}")
        rho = (res[8]["theta"].log() / 8).exp()
        for j in (10, 11, 12):                                 # D_j <= theta_5 D_(j-5)
            Dq[j] = th5 * Dq[j - 5]
        C = arb(max((Dq[j] / rho ** j).upper() for j in range(5, 13)))
        check(bool(rho < cap), f"DR.4 envelope: tanh(D_n/4) <= C rho^n for every n >= 5 (n = 8q + j, j = 5..12, "
              f"D_n <= theta_8^q D_j, D_j <= theta_5 D_(j-5) for j >= 10): rho = theta_8^(1/8) <= "
              f"{rho.upper().str(7, radius=False)}, C <= {C.upper().str(6, radius=False)}")
        RESULT.update(unif=unif, rho=rho, C=C)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem DR checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
