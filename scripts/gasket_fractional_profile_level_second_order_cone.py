#!/usr/bin/env python3
r"""Theorem DQ: the Birkhoff rate of a cone is the spectral radius of its order-unit space.  On every C_a (a > b*) the
cone-diameter rate tends to the C^1 cap rho(L_1)/r exactly (DP.3's cap is attained, not only a floor); on the second-order
cone K_(2,6) = {f > 0 : |f'| <= 2 f, |f''| <= 6 f} (invariant, arb) it tends to CN's spectral rate |l_2|/r = 0.17255, so
the cone-diameter route beats the cap for all long blocks; Kept Failure DQ.5: the threshold is not certified, and blocks
Lt, Lt^2, Lt^3 on K_(2,6) are certified to stay above the cap.  (DP's kept limit "only a C^2-type cone escapes the
obstruction (then the floor is |l_2|/r)": whether such a cone exists, is invariant and attains that floor.)

Setting (GASKET.md, Theorems CE, CF, CI-CM, DL-DP): s in X = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
wt1 = (s+1/2)/(2s(7/2-s)), wt2 = (s-1/2)/(2s(s+7/2)), sigma_(1,2) = sqrt(15/4 -+ s); C_a = {f > 0: log f a-Lipschitz in s}
(CF), equivalently {f in W^(1,oo): f >= 0, a f +- f' >= 0 a.e.} \ {0}; K_(a,b) = {f in W^(2,oo): f >= 0, a f +- f' >= 0,
b f +- f'' >= 0 a.e.} \ {0}.  Theta_K = Hilbert metric (CF's definition), tau(T; K) = sup Theta(Tf,Tg)/Theta(f,g), D_n =
diam_Theta Lt^n K.  E_1 = W^(1,oo)(X), E_2 = W^(2,oo)(X) with ||f|| = sum_(j<=k) sup|f^(j)|.  M_k g = sum wt_i sigma_i'^k g o sigma_i
(so |M_k^n g| <= L_k^n |g|, CH's L_k).

Theorem DQ (checked here; exit 1 on any failure):
  DQ.1 [proved; sympy identities, arb cover] (Lt f)' = M_1 f' + sum wt_i' f o sigma_i and (Lt f)'' = M_2 f'' + sum [wt_i''
       f o sigma_i + (2 wt_i' sigma_i' + wt_i sigma_i'') f' o sigma_i] (sympy).  Hence Lt K_(a,b) subset K_(a',b') with
       a' = max_i sup (|wt_i'/wt_i| + a|sigma_i'|), b' = max_i sup (|wt_i''/wt_i| + (2|wt_i'/wt_i||sigma_i'| + |sigma_i''|) a +
       sigma_i'^2 b); for (a, b) = (2, 6): a' < 1.8280, b' < 5.149, and b_min(2) = max_i sup(...)/(1 - sigma_i'^2) < 4.936, so
       every K_(2,b), b >= 4.936, is invariant.  h in K_(a',b') (closed limit of Lt^n 1/r^n, CF.3), so l(h) >= m > 0 for all
       defining functionals l, and V_h(K_(a,b)) = E_2: K_(a,b) is a C^2-type cone (no C^1 eigenfunction of CI.4 outside E_2
       is h-bounded).
  DQ.2 [proved; CI.1, CL.3 re-run; CJ/CK, CM.4, CF.1-CF.3 as stated; textbook Fredholm/Riesz/Gelfand] On E_k, Lt = S_k + K_k
       with S_k = 0 (+) M_k (f <-> (jet at s0, f^(k))), K_k compact, rho(S_k) <= rho(L_k).  Bootstrap: an E_k eigenfunction
       with |mu| > rho(L_k) has continuous f^(k) (v^(k) = sum mu^(-j-1) M_k^j q), so is a C^k eigenfunction; r is
       algebraically simple on E_k (cone argument).  So spec(Lt; E_1) n {|z| > rho(L_1)} = {r}, spec(Lt; E_2) n {|z| > rho(L_2)}
       = {r, l_2}, and ||Lt^n (1 - Pi)||_(E_k)^(1/n) -> R_k, R_1 = rho(L_1), R_2 = |l_2| (Pi f = nu(f) h).
  DQ.3 [proved from DQ.2, DP.2/DP.3, CF.1] For every a > b*: lim tau(Lt^n; C_a)^(1/n) = lim tanh(D_n/4)^(1/n) = lim D_n^(1/n)
       = rho(L_1)/r in [0.3170757, 0.3170762]: the cap of DP.3 is attained by every C_a (DM's ratio ~0.3175 is this limit).
  DQ.4 [proved from DQ.1, DQ.2, DP.2] For K = K_(2,6): lim tau(Lt^n; K)^(1/n) = lim tanh(D_n/4)^(1/n) = lim D_n^(1/n) =
       |l_2|/r in [0.1725475, 0.1725477] (DP.4's floor, now attained); so there is n0 with tanh(D_n/4)^(1/n) < 0.3170757 for
       every n >= n0, and the cone-diameter route on K_(2,6) reaches CN's exact rate in the limit.
  DQ.5 [kept failure, certified witnesses] n0 is not effective.  On K_(2,6), with an explicit polynomial base p (int K, every
       l(p) > 0.8434) and explicit polynomial directions phi_n (float LP proposals, coefficients below), DP.1 gives
       tau(Lt^n; K) >= N_(Lt^n p)(Lt^n phi_n) / N_p(phi_n) (N_u(phi) = osc_l l(phi)/l(u)) with arb: per-level > 0.5523, > 0.3922,
       > 0.3221 for n = 1, 2, 3, all above the cap 0.3170762: so n0 >= 4 on K_(2,6), and no CF.1-form bound with block <= 3
       on this cone beats DP.3's cap.
Kept limits: DQ.3/DQ.4 are limits (Gelfand), with no explicit constant; n0 and the level-sum constants on K_(2,6) are not
certified; floats (DQ_LP=1, needs scipy, ~6.5 min) put the near-h ratio per level at 0.28560, 0.25942, 0.24259 for n = 4, 5, 6,
below the cap but not certified (grid relaxation).  Textbook facts used: Arzela-Ascoli, compact perturbations and the analytic Fredholm alternative,
Riesz projections, Gelfand's formula.  No Psi_info, Stage 2 still blocked.  Stage 1 untouched; no Psi_info, W, 0.08, 0.23
or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time
from fractions import Fraction as Fr

import sympy as sp
from flint import arb, arb_series, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_essential as ci          # noqa: E402  (CI.1 re-run: rho(L_1), r)
import gasket_fractional_profile_level_calkin as cl             # noqa: E402  (CL.3 re-run: rho(L_2); l_2)

FAILS: list[str] = []
RESULT: dict = {}
A_CONE, B_CONE = 2, 6                 # K_(2,6)
N_CELL_INV = int(os.environ.get("DQ_INV_CELLS", "20000"))
N_CELL_OSC = int(os.environ.get("DQ_OSC_CELLS", "8000"))
TARGET = dict(a1="1.8280", b1="5.149", bmin="4.936", q1="0.3170757", q1_up="0.3170762", q2_lo="0.1725475",
              q2_up="0.1725477", minlp="0.8434", per_level={1: "0.5523", 2: "0.3922", 3: "0.3221"})
# DQ.5 witnesses: Chebyshev coefficients (in t = (2s - (sqrt5/2 + 5/2))/(5/2 - sqrt5/2)) of phi_n, from a float LP
# (lp_proposal below, degree 39, 60 x 5 image functionals, 600 x 5 constraint functionals; h-component removed, max
# coefficient scaled to 1, rounded to 10 digits).  Any phi is a valid witness; only the arb check below is claimed.
PHI = {
    1: [0.0, 1.0, -0.3383807631, -0.02441960787, 0.04028933615, -0.02163907472, 0.002724402369, 0.004039340635, -0.003887943838, 0.002115388299, -0.0001647439279, -0.0009472266851, 0.001301288753, -0.0006305701211, -0.000202516114, 0.0007358055017, -0.0004882658157, -9.785063756e-07, 0.0003825451875, -0.0003006172453, 6.43847151e-05, 0.0001541735482, -0.0001510769271, 7.493461036e-05, 2.843359592e-05, -6.370671329e-05, 7.211276166e-05, -3.406686517e-05, -1.635599069e-05, 5.729396668e-05, -5.168549049e-05, 4.059109214e-06, 3.670434207e-05, -4.261002846e-05, 7.79735759e-06, 1.882022155e-05, -2.691493458e-05, 7.683176617e-06, 2.274245422e-06, -6.849694236e-06],
    2: [0.0, -1.0, 0.6941133071, 0.237055752, -0.1698650747, -0.02494234734, 0.004555256186, -0.03364549477, -0.005016417137, 0.01193387858, -0.001246444521, 0.0009875644126, 0.005864234145, -0.000701115463, -0.002232235875, 0.0006577152362, -0.0008600409701, -0.001413424774, 0.0006259230079, 0.0005037650229, -0.0001978589453, 0.0004433443639, 0.0003258164725, -0.0003106009655, -0.0001085255136, 3.482817126e-05, -0.0001842007502, -4.919239208e-05, 0.0001248646783, 1.969273584e-05, 3.611825802e-06, 6.046874869e-05, -4.653516078e-06, -4.022722638e-05, -3.035600586e-06, -4.341388552e-06, -1.353756735e-05, 5.349684618e-06, 9.032634576e-06, 3.598468621e-07],
    3: [0.0, 1.0, -0.9255389997, -0.2081219255, 0.2015921415, 0.01448283952, 0.008654980129, 0.04130412917, -0.009006437042, -0.01115756309, 0.001354406973, -0.006472896557, -0.005457299249, 0.003843279931, 0.0004229881787, -0.0004490985248, 0.002549690101, -0.0003984025444, -0.001351169691, 7.233602811e-05, -0.0003769111318, -0.0008248672497, 0.00071539232, 0.0003901828336, -4.042496951e-05, 0.000520928116, 8.939434e-05, -0.000241954297, -1.375790276e-05, 4.901190148e-05, -0.0002075849369, 0.0001551306274, 0.0001111487819, -1.632437261e-05, 6.67474305e-05, 2.918328657e-05, -7.549127374e-05, -4.793445872e-05, 1.091707413e-06, -6.330467943e-05],
}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def _q(x) -> arb:
    f = Fr(str(x))
    return arb(f.numerator) / f.denominator


# ----------------------------------------------------------------- DQ.1 identities (sympy)
def symbolic_checks():
    s = sp.symbols("s", positive=True)
    f = sp.Function("f")
    half = sp.Rational(1, 2)
    W = [(s + half) / (2 * s * (sp.Rational(7, 2) - s)), (s - half) / (2 * s * (s + sp.Rational(7, 2)))]
    S = [sp.sqrt(sp.Rational(15, 4) - s), sp.sqrt(sp.Rational(15, 4) + s)]
    Lf = sum(W[i] * f(S[i]) for i in range(2))
    y = sp.symbols("y")
    fp = lambda i, k: sp.diff(f(y), y, k).subs(y, S[i])
    d1 = sum(sp.diff(W[i], s) * f(S[i]) + W[i] * sp.diff(S[i], s) * fp(i, 1) for i in range(2))
    d2 = sum(W[i] * sp.diff(S[i], s) ** 2 * fp(i, 2) + sp.diff(W[i], s, 2) * f(S[i])
             + (2 * sp.diff(W[i], s) * sp.diff(S[i], s) + W[i] * sp.diff(S[i], s, 2)) * fp(i, 1) for i in range(2))
    ok1 = sp.simplify((sp.diff(Lf, s) - d1).doit()) == 0
    ok2 = sp.simplify((sp.diff(Lf, s, 2) - d2).doit()) == 0
    ok3 = all(sp.simplify(sp.diff(S[i], s) ** 2 - 1 / (4 * S[i] ** 2)) == 0 for i in range(2))
    check(ok1 and ok2 and ok3,
          "DQ.1 (sympy): (Lt f)' = M_1 f' + sum wt_i' f o sigma_i, (Lt f)'' = M_2 f'' + sum [wt_i'' f o sigma_i + "
          "(2 wt_i' sigma_i' + wt_i sigma_i'') f' o sigma_i], sigma_i'^2 = 1/(4 sigma_i^2) (M_k g = sum wt_i sigma_i'^k g o sigma_i)")
    # bootstrap identity: if mu v = Lt v then mu v^(k) = M_k v^(k) + q, and v^(k) = sum_(j<J) mu^(-j-1) M_k^j q + mu^(-J) M_k^J v^(k)
    mu, M, q, v = sp.symbols("mu M q v")
    J = 4
    rhs = sum(mu ** (-j - 1) * M ** j * q for j in range(J)) + mu ** (-J) * M ** J * v
    ok4 = sp.simplify(rhs.subs(q, mu * v - M * v) - v) == 0
    check(ok4, "DQ.2 (sympy, commuting symbols): from mu v = M v + q, v = sum_(j<J) mu^(-j-1) M^j q + mu^(-J) M^J v (J = 4); "
               "with ||M_k^J|| <= ||L_k^J 1|| and |mu| > rho(L_k) the tail -> 0, so v^(k) is a uniform limit of continuous functions")


# ----------------------------------------------------------------- DQ.1 invariance (arb)
def _w(i, xs):
    h, sv = arb(1) / 2, arb(7) / 2
    return (xs + h) / (2 * xs * (sv - xs)) if i == 0 else (xs - h) / (2 * xs * (xs + sv))


def _sig(i, xs):
    qv = arb(15) / 4
    return (qv - xs).sqrt() if i == 0 else (qv + xs).sqrt()


def _jets(ser, k):
    c = ser.coeffs() + [arb(0)] * (k + 1)
    return [c[j] * math.factorial(j) for j in range(k + 1)]


A_X, B_X = arb(5).sqrt() / 2, arb(5) / 2


def _cells(n):
    for k in range(n):
        lo = A_X + (B_X - A_X) * k / n
        hi = A_X + (B_X - A_X) * (k + 1) / n
        c = (lo + hi) / 2
        yield arb(c.mid()), ((hi - lo) / 2).upper() + c.rad()


def invariance():
    ctx.prec = 128
    a, b = arb(A_CONE), arb(B_CONE)
    a1 = b1 = bmin = arb(0)
    for c, rad in _cells(N_CELL_INV):
        ball = arb(c.mid(), rad)
        xs = arb_series([ball, 1], prec=3)
        for i in range(2):
            w0, w1, w2 = _jets(_w(i, xs), 2)
            y0, y1, y2 = _jets(_sig(i, xs), 2)
            assert w0 > 0
            l1, l2 = abs(w1 / w0), abs(w2 / w0)
            ea = l1 + a * abs(y1)
            Bi = l2 + (2 * l1 * abs(y1) + abs(y2)) * a
            eb = Bi + y1 * y1 * b
            em = Bi / (1 - y1 * y1)
            a1, b1, bmin = (arb(max(a1.upper(), ea.upper())), arb(max(b1.upper(), eb.upper())),
                            arb(max(bmin.upper(), em.upper())))
    check(bool(a1 < arb(TARGET["a1"])) and bool(b1 < arb(TARGET["b1"])) and bool(bmin < arb(TARGET["bmin"])),
          f"DQ.1 (arb, {N_CELL_INV} second-order balls): Lt K_(2,6) subset K_(a',b'), a' <= {a1.upper().str(6, radius=False)} < "
          f"{TARGET['a1']} < 2, b' <= {b1.upper().str(6, radius=False)} < {TARGET['b1']} < 6; b_min(2) <= "
          f"{bmin.upper().str(6, radius=False)} < {TARGET['bmin']} (every K_(2,b), b >= {TARGET['bmin']}, is invariant)")
    RESULT.update(a1=a1, b1=b1, bmin=bmin)


# ----------------------------------------------------------------- DQ.2-DQ.4 constants (arb)
def spectral_constants():
    ctx.prec = 128
    r_lo, r_hi = arb(ci.R_LO), arb(ci.R_HI)
    q1 = arb(ci.LAM_LO) / r_hi
    q1_up = arb(ci.LAM_UP) / r_lo
    l2 = arb(cl.L2_ABS, float(cl.L2_ERR))
    q2_lo, q2_up = l2 / r_hi, l2 / r_lo
    rl2_up = arb(cl.CW[2][1])
    check(bool(rl2_up < l2) and bool(l2 < arb(ci.LAM_LO)),
          f"DQ.2: rho(L_2) <= {cl.CW[2][1]} < |l_2| = {cl.L2_ABS} < rho(L_1) >= {ci.LAM_LO}: on E_2 the second spectral point is "
          f"l_2 (R_2 = |l_2|); on E_1 l_2 lies inside the essential disk (R_1 = rho(L_1))")
    check(bool(q1 > arb(TARGET["q1"])) and bool(q1_up < arb(TARGET["q1_up"])),
          f"DQ.3: lim tau(Lt^n; C_a)^(1/n) = rho(L_1)/r in [{q1.lower().str(9, radius=False)}, {q1_up.upper().str(9, radius=False)}]")
    check(bool(q2_lo > arb(TARGET["q2_lo"])) and bool(q2_up < arb(TARGET["q2_up"])) and bool(q2_up < q1),
          f"DQ.4: lim tau(Lt^n; K_(2,6))^(1/n) = |l_2|/r in [{q2_lo.lower().str(9, radius=False)}, "
          f"{q2_up.upper().str(9, radius=False)}] < cap; ratio cap/limit > {(q1 / q2_up).lower().str(6, radius=False)}")
    RESULT.update(q1=q1, q1_up=q1_up, q2=(q2_lo, q2_up))


# ----------------------------------------------------------------- DQ.5 witnesses (arb)
def _cheb(coef, xs):
    t = (2 * xs - (A_X + B_X)) / (B_X - A_X)
    b1 = 0 * t
    b2 = 0 * t
    for c in reversed(coef[1:]):
        b1, b2 = 2 * t * b1 - b2 + c, b1
    return t * b1 - b2 + coef[0]


def _Lt_ev(coef, n, xs):
    if n == 0:
        return _cheb(coef, xs)
    return _w(0, xs) * _Lt_ev(coef, n - 1, _sig(0, xs)) + _w(1, xs) * _Lt_ev(coef, n - 1, _sig(1, xs))


def _ells(J):
    a, b = arb(A_CONE), arb(B_CONE)
    return [J[0], a * J[0] + J[1], a * J[0] - J[1], b * J[0] + J[2], b * J[0] - J[2]]


def _ells_d(J):
    a, b = arb(A_CONE), arb(B_CONE)
    return [J[1], a * J[1] + J[2], a * J[1] - J[2], b * J[1] + J[3], b * J[1] - J[3]]


def _ball_jets(coef, c, rad, k):
    """f^(j) on [c - rad, c + rad] (j <= k) from the exact Taylor expansion of the polynomial at c, and at c."""
    D = len(coef)
    a = _cheb(coef, arb_series([c, 1], prec=D)).coeffs() + [arb(0)] * D
    d = arb(0, rad)
    ball = []
    for j in range(k + 1):
        acc = arb(0)
        for m in range(D - 1, j - 1, -1):
            acc = acc * d + a[m] * (math.factorial(m) // math.factorial(m - j))
        ball.append(acc)
    return ball, [a[j] * math.factorial(j) for j in range(k + 1)]


def _float_base(deg=20):
    """Float proposal only: Chebyshev collocation eigenfunction of Lt (base point p ~ h)."""
    import numpy as np
    from numpy.polynomial import chebyshev as C
    a0, b0 = 5 ** 0.5 / 2, 2.5
    tm = lambda y: (2 * y - (a0 + b0)) / (b0 - a0)

    def Lt(fun, x):
        w1 = (x + .5) / (2 * x * (3.5 - x))
        w2 = (x - .5) / (2 * x * (x + 3.5))
        return w1 * fun(np.sqrt(3.75 - x)) + w2 * fun(np.sqrt(3.75 + x))

    n = 50
    z = np.cos(np.pi * (np.arange(n) + .5) / n)
    X = (a0 + b0) / 2 + (b0 - a0) / 2 * z
    V = C.chebvander(z, n - 1)
    M = np.zeros((n, n))
    for j in range(n):
        e = np.zeros(n)
        e[j] = 1
        M[:, j] = Lt(lambda y: C.chebval(tm(y), e), X)
    ev, vec = np.linalg.eig(np.linalg.solve(V, M))
    i = int(np.argmax(abs(ev)))
    co = vec[:, i].real
    co = co / C.chebval(0, co)
    return [Fr(float(c)).limit_denominator(10 ** 15) for c in co[:deg + 1]]


def witnesses():
    ctx.prec = 128
    pco = [arb(c.numerator) / c.denominator for c in _float_base()]
    cap_up = RESULT["q1_up"]
    out = {}
    minlp_all = None
    for n, raw in PHI.items():
        phico = [_q(x) for x in raw]
        sup = inf = minlp = None
        for c, rad in _cells(N_CELL_OSC):
            Jb_p, Jc_p = _ball_jets(pco, c, rad, 3)
            Jb_f, Jc_f = _ball_jets(phico, c, rad, 3)
            lp_b, lf_b, dp_b, df_b = _ells(Jb_p), _ells(Jb_f), _ells_d(Jb_p), _ells_d(Jb_f)
            lp_c, lf_c = _ells(Jc_p), _ells(Jc_f)
            for j in range(5):
                lb = lp_b[j].lower()
                minlp = lb if minlp is None else min(minlp, lb)
                rc = lf_c[j] / lp_c[j]
                dr = (df_b[j] * lp_b[j] - lf_b[j] * dp_b[j]) / (lp_b[j] * lp_b[j])
                hi = (rc + abs(dr) * arb(rad)).upper()
                lo = (rc - abs(dr) * arb(rad)).lower()
                sup = hi if sup is None else max(sup, hi)
                inf = lo if inf is None else min(inf, lo)
        n_up = arb(sup) - arb(inf)
        # image: point functionals at 61 points, the two best
        vals = []
        for k in range(61):
            x = A_X + (B_X - A_X) * k / 60
            Jp = _jets(_Lt_ev(pco, n, arb_series([x, 1], prec=3)), 2)
            Jf = _jets(_Lt_ev(phico, n, arb_series([x, 1], prec=3)), 2)
            lp, lf = _ells(Jp), _ells(Jf)
            ok_pos = all(bool(v > 0) for v in lp)
            if not ok_pos:
                check(False, f"DQ.5: l(Lt^{n} p) > 0 at point {k}")
            for j in range(5):
                vals.append((lf[j] / lp[j], k, j))
        hi = max(vals, key=lambda v: v[0].lower())
        lo = min(vals, key=lambda v: v[0].upper())
        n_lo = hi[0].lower() - lo[0].upper()
        tau_lo = arb(n_lo) / arb(n_up.upper())
        per = arb(tau_lo.lower()) ** (arb(1) / n)
        out[n] = dict(n_up=n_up, n_lo=arb(n_lo), tau=tau_lo, per=per, pair=(hi[1:], lo[1:]))
        minlp_all = minlp if minlp_all is None else min(minlp_all, minlp)
        names = ["d", "a d + d'", "a d - d'", "b d + d''", "b d - d''"]
        check(bool(per.lower() > arb(TARGET["per_level"][n])) and bool(per.lower() > cap_up.upper()),
              f"DQ.5 (arb): n = {n}: N_p(phi_{n}) <= {n_up.upper().str(7, radius=False)} ({N_CELL_OSC} cells), "
              f"N_(Lt^{n} p)(Lt^{n} phi_{n}) >= {arb(n_lo).lower().str(7, radius=False)} "
              f"({names[hi[2]]} at s_{hi[1]} vs {names[lo[2]]} at s_{lo[1]}, s_k = sqrt5/2 + k(5/2 - sqrt5/2)/60), so "
              f"tau(Lt^{n}; K_(2,6)) >= {tau_lo.lower().str(7, radius=False)}, per level >= {per.lower().str(7, radius=False)} "
              f"> {TARGET['per_level'][n]} > cap {TARGET['q1_up']}")
    check(bool(arb(minlp_all) > arb(TARGET["minlp"])),
          f"DQ.5 (arb): base p (degree-20 rational Chebyshev, float eigenfunction proposal) lies in int K_(2,6): every "
          f"l(p) >= {arb(minlp_all).lower().str(6, radius=False)} > {TARGET['minlp']} on X (so p + t phi_n in K for small t, DP.1 applies)")
    RESULT.update(wit=out)


# ----------------------------------------------------------------- optional float LP (proposal of PHI; no claim)
def lp_proposal(nmax=6, D=40, nx=60, nxb=600, step=3):
    import numpy as np
    from numpy.polynomial import chebyshev as C
    from scipy.optimize import linprog
    a0, b0 = 5 ** 0.5 / 2, 2.5
    tm = lambda y: (2 * y - (a0 + b0)) / (b0 - a0)
    sc = 2 / (b0 - a0)
    xs_ = sp.symbols("x")
    Wsym = [(xs_ + sp.Rational(1, 2)) / (2 * xs_ * (sp.Rational(7, 2) - xs_)),
            (xs_ - sp.Rational(1, 2)) / (2 * xs_ * (xs_ + sp.Rational(7, 2)))]
    Wd = [[sp.lambdify(xs_, sp.diff(w, xs_, k), "numpy") for k in range(3)] for w in Wsym]
    hco = [float(c) for c in _float_base(40)]
    r = 0.42373128696260794

    def jet_poly(coef):
        d1 = C.chebder(coef)
        d2 = C.chebder(d1)
        return lambda y: (C.chebval(tm(y), coef), sc * C.chebval(tm(y), d1), sc * sc * C.chebval(tm(y), d2))

    def Ltj(J):
        # jets by the exact chain rule (DQ.1's identities), float proposal only
        def K(x):
            o = [0, 0, 0]
            for i, sg in ((0, -1), (1, 1)):
                y = np.sqrt(3.75 + sg * x)
                y1 = sg / (2 * y)
                y2 = -1 / (4 * y ** 3)
                w0, w1, w2 = (Wd[i][k](x) for k in range(3))
                f0, f1, f2 = J(y)
                o[0] = o[0] + w0 * f0
                o[1] = o[1] + w1 * f0 + w0 * y1 * f1
                o[2] = o[2] + w2 * f0 + 2 * w1 * y1 * f1 + w0 * (y2 * f1 + y1 * y1 * f2)
            return o
        return K

    def ells(J, x):
        f0, f1, f2 = J(x)
        return np.concatenate([f0, A_CONE * f0 + f1, A_CONE * f0 - f1, B_CONE * f0 + f2, B_CONE * f0 - f2])

    for n in range(1, nmax + 1):
        x = np.linspace(a0, b0, nx)
        xb = np.linspace(a0, b0, nxb)
        hj = jet_poly(np.array(hco))
        cols, colsb = [], []
        for j in range(D):
            e = np.zeros(D)
            e[j] = 1
            J = jet_poly(e)
            colsb.append(ells(J, xb))
            for _ in range(n):
                J = Ltj(J)
            cols.append(ells(J, x))
        Am = np.array(cols).T / (r ** n * ells(hj, x))[:, None]
        Bm = np.array(colsb).T / ells(hj, xb)[:, None]
        m = Am.shape[0]
        Aub = np.vstack([np.hstack([-Bm, np.ones((len(Bm), 1))]), np.hstack([Bm, -np.ones((len(Bm), 1))])])
        bub = np.concatenate([np.zeros(len(Bm)), np.ones(len(Bm))])
        best = 0
        for k in range(0, m, step):
            j = (k + m // 2) % m
            seen = set()
            while (k, j) not in seen:
                seen.add((k, j))
                cv = np.zeros(D + 1)
                cv[:D] = -(Am[k] - Am[j])
                res = linprog(cv, A_ub=Aub, b_ub=bub, bounds=[(-1e4, 1e4)] * D + [(None, None)], method="highs")
                if res.status != 0:
                    break
                v = Am @ res.x[:D]
                best = max(best, v.max() - v.min())
                j = int(np.argmin(v))
        print(f"     observation (floats, DQ_LP=1, grid LP, no claim): near-h ratio n = {n}: {best:.6f}, per level {best ** (1 / n):.5f}", flush=True)


def float_M1_observation():
    """Floats only (no claim): top eigenvalues of the signed operator M_1 on smooth functions (Chebyshev collocation)."""
    import numpy as np
    from numpy.polynomial import chebyshev as C
    a0, b0 = 5 ** 0.5 / 2, 2.5
    tm = lambda y: (2 * y - (a0 + b0)) / (b0 - a0)
    n = 50
    z = np.cos(np.pi * (np.arange(n) + .5) / n)
    X = (a0 + b0) / 2 + (b0 - a0) / 2 * z
    V = C.chebvander(z, n - 1)
    out = {}
    for name, sgn in (("M_1", -1.0), ("L_1", 1.0)):
        M = np.zeros((n, n))
        for j in range(n):
            e = np.zeros(n)
            e[j] = 1
            y1, y2 = np.sqrt(3.75 - X), np.sqrt(3.75 + X)
            w1 = (X + .5) / (2 * X * (3.5 - X))
            w2 = (X - .5) / (2 * X * (X + 3.5))
            M[:, j] = sgn * w1 / (2 * y1) * C.chebval(tm(y1), e) + w2 / (2 * y2) * C.chebval(tm(y2), e)
        ev = np.linalg.eigvals(np.linalg.solve(V, M))
        out[name] = ev[np.argsort(-abs(ev))][:2]
    r = float(ci.R_LO)
    m1 = out["M_1"][0].real
    print(f"     observation (floats, no claim): top eigenvalue of M_1 on smooth functions {m1:.5f} (|.|/r = {abs(m1) / r:.4f}); "
          f"L_1: {out['L_1'][0].real:.7f} (= rho(L_1)); a word-by-word absolute bound sees only L_1, hence the cap", flush=True)
    RESULT.update(m1=m1)


def main() -> int:
    t0 = time.time()
    print("Theorem DQ: Birkhoff rates are spectral radii: C_a attains the C^1 cap rho(L_1)/r; the second-order cone K_(2,6) "
          "is invariant and tends to |l_2|/r = 0.17255; n <= 3 blocks on K_(2,6) stay above the cap. Not thrust.", flush=True)
    symbolic_checks()
    ci.collatz_two_sided()          # CI.1 re-run: rho(L_1) in [0.1343549, 0.1343551]
    cl.collatz_two_sided(2)         # CL.3 re-run: rho(L_2) in [0.0430702, 0.0430704]
    spectral_constants()
    invariance()
    witnesses()
    float_M1_observation()
    if os.environ.get("DQ_LP") == "1":
        lp_proposal()
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(ci.FAILS) + len(cl.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DQ checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
