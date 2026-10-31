#!/usr/bin/env python3
"""Theorem CX: the exact geometric rate of 3^(k+1) rho(L_k) - 1 (CW's kept limit "rate only bracketed in
[0.5999, 0.8352]") and the monotonicity of 3^k rho(L_k) for every k (CW's open item "k >= 13").

Setting exactly as in Theorem CW (GASKET.md): X = [sqrt5/2, 5/2], sigma_(1,2)(s) = sqrt(15/4 -+ s),
(L_k f)(s) = u_1 f(sigma_1 s) + u_2 f(sigma_2 s), u_i = wt_i/(2 sigma_i)^k, rho(L_k) = spectral radius on C(X), s* = 3/2.
CW's h_k = exp(k Phi + Psi) conjugates L_k to 3^(-k-1) Lt_k with (Lt_k f)(s) = f(sigma_1 s) + a(s) f(sigma_2 s),
a = G B^k (CW's G, B; their routines are imported from CW's script and run unchanged).  Write
eps_k = 3^(k+1) rho(L_k) - 1 (> 0 by CW.3), Y = sigma_2(X) = [sqrt(15/4 + sqrt5/2), 5/2], x_n(t) = sigma_1^n(t),
B* = B(3/2) = (3/sqrt21) exp(Phi(sqrt21/2)), G* = G(3/2), c_k = eps_k/(G* B*^k).

Renewal test function (exact).  For eps > 0 put u_eps(t) = sum_(n>=0) (1+eps)^(-n-1) a(x_n(t)) (uniformly convergent,
continuous, > 0).  Telescoping gives (1+eps) u_eps = u_eps o sigma_1 + a, hence EXACTLY
    Lt_k u_eps = (1+eps) u_eps + a (u_eps o sigma_2 - 1).
Since a > 0, Collatz-Wielandt for the positive operator L_k on C(X) with test function h_k u_eps gives
    sup_Y u_eps <= 1  =>  eps_k <= eps,          inf_Y u_eps >= 1  =>  eps_k >= eps.                  (R)

Theorem CX (checked here; exit 1 on any failure).
  CX.1 [proved, exact + arb] |sigma_1 x - s*| <= theta0 |x - s*| on X' with theta0 = (3 - sqrt5)/2; on J = [s* - 1/8, s* + 1/8]
       the derivatives of log B, log G (chain rule along the orbit, J_D terms + explicit tails, adaptive arb cover) satisfy
       |(log B)'| <= LAM_B, |(log G)'| <= LAM_G; B*, G* enclosed.
  CX.2 [certified, arb, per integer k = 0..K-1] With a cover of Y, the first N_H orbit terms enclosed per cell and the
       n >= N_H tail bracketed by G* B*^k (1+eps)^(-N_H)/eps * exp(+-(LAM_G + k LAM_B) theta0^N_H), (R) gives an explicit
       bracket c_k in [c_lo(k), c_hi(k)].
  CX.3 [proved from CX.1 and CW.2, every REAL k >= K] With N(k) = least N with (LAM_G + k LAM_B) theta0^N <= TAU,
       sup_Y u <= N(k) 0.2873 0.8352^k + e^TAU/(1+DELTA) <= 1 at eps = (1+DELTA) G* B*^k and
       inf_Y u >= (1 - N(k) eps) e^(-TAU)/(1-DELTA) >= 1 at eps = (1-DELTA) G* B*^k (both monotone in k, checked at K):
       1 - DELTA <= c_k <= 1 + DELTA for every real k >= K.
  CX.4 [proved from CX.1-CX.3] (a) lim_(k->oo) (3^(k+1) rho(L_k) - 1)/(G* B*^k) = 1 (for every delta > 0 the CX.3
       argument with tau = delta/3 closes at some K(delta)); so the exact rate is B* = 0.7773760276... and the constant
       G* = 0.2870652417...  (b) 3^k rho(L_k) > 3^(k+1) rho(L_(k+1)) for every integer k >= 0 and every real k >= K,
       i.e. 3^k rho(L_k) decreases strictly to 1/3 (CV.7 extended from k <= 12 to all k, independently of CV).
  CX.5 [consistency] CE's r and CT/CU/CV's certified rho(L_k), k <= 12 (as stated), lie inside CX.2's brackets.
Kept limits: c_k only bracketed (width ~5e-3 at k = 13, ~1.5e-4 at k = 39; +-DELTA = +-0.02 for real k >= K = 40); nothing about
nu_2^(k), so q^(k) -> -1/3 stays HEURISTIC; no uniform-in-k rate for inf R_n^(k).
Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction as Fr

from flint import arb, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_uniform_k as cw    # noqa: E402

FAILS: list[str] = []
RESULT: dict = {}
PREC = 160
K = 40                         # CX.2 for k = 0..K-1, CX.3 for real k >= K
N_H = 60                       # orbit terms enclosed per cell in CX.2
NCELL = 64                     # cells covering Y
D0 = Fr(1, 8)                  # J = [s* - D0, s* + D0]
J_D = 80                       # orbit terms in the derivative sums
TAU = Fr(8, 1000)
DELTA = Fr(2, 100)
YLO = Fr(2259533, 2 ** 20)     # dyadic <= sigma_2(sqrt5/2) = 2.20663..., >= cw.LO (checked)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


_a = cw._a
ball = cw.ball
sig1, sig2, wt1, wt2 = cw.sig1, cw.sig2, cw.wt1, cw.wt2
G1 = None   # 1/(1 - 5^(-1/2))


def dphipsi(x, dx):
    """Enclosures of Phi'(x0) dx0 and Psi'(x0) dx0 over the ball x (subset of X'), |dx| <= 1 assumed for the tails:
    x_(j+1)' = sigma_1'(x_j) x_j' with |sigma_1'| = 1/(2 sigma_1) <= 1/sqrt5 on X', so |x_j'| <= 5^(-j/2);
    Phi' = -sum x_(j+1)'/x_(j+1) (tail <= (2/sqrt5) 5^(-(J+1)/2) G1), Psi' = sum (log wt_1)'(x_j) x_j' (|.| <= 3: tail <= 3 5^(-J/2) G1)."""
    dp = arb(0)
    dq = arb(0)
    for _ in range(J_D):
        y = sig1(x)
        dy = -dx / (2 * y)
        dp -= dy / y
        dq += (1 / (x + arb(1) / 2) - 1 / x + 1 / (arb(7) / 2 - x)) * dx
        x, dx = y, dy
    s5 = arb(5).sqrt()
    tp = 2 / s5 * arb(5) ** (-arb(J_D + 1) / 2) * G1
    tq = 3 * arb(5) ** (-arb(J_D) / 2) * G1
    return dp + arb(0, 1) * tp, dq + arb(0, 1) * tq


def dlogBG(I):
    y = sig2(I)
    d2 = 1 / (2 * y)                              # sigma_2' in (0, 1/sqrt5]
    p0, q0 = dphipsi(I, arb(1))
    p1, q1 = dphipsi(y, d2)
    db = -d2 / y + p1 - p0
    dg = 1 / (I - arb(1) / 2) - 1 / I - 1 / (I + arb(7) / 2) + q1 - q0
    return db, dg


def _dn(x: Fr) -> str:
    """x rounded DOWN to 6 decimals (outward for a lower end)."""
    return f"{math.floor(x * 10 ** 6) / 10 ** 6:.6f}"


def _up(x: Fr) -> str:
    """x rounded UP to 6 decimals (outward for an upper end)."""
    return f"{math.ceil(x * 10 ** 6) / 10 ** 6:.6f}"


def upper(x: arb) -> arb:
    return arb(x.mid() + x.rad())


def lower(x: arb) -> arb:
    return arb(x.mid() - x.rad())


def cx1() -> bool:
    global G1
    G1 = 1 / (1 - 1 / arb(5).sqrt())
    ok = True
    # theta0: |sigma_1 x - s*| = |x - s*|/(sigma_1 x + 3/2), sigma_1 >= sqrt5/2 on X'  =>  factor <= 2/(3 + sqrt5) = (3-sqrt5)/2
    th0 = (3 - arb(5).sqrt()) / 2
    # Y-cover inside X': cw.LO <= YLO <= sigma_2(sqrt5/2)  <=>  (YLO^2 - 15/4)^2 <= 5/4 with YLO^2 >= 15/4
    ok &= cw.LO <= YLO and YLO * YLO >= Fr(15, 4) and (YLO * YLO - Fr(15, 4)) ** 2 <= Fr(5, 4)
    # J inside X'; sigma_1 maps J into J (contraction towards s* by theta0 < 1); points of X' have |x - s*| <= 1
    ok &= Fr(3, 2) - D0 >= cw.LO
    lam_b = lam_g = arb(0)
    stack = [(Fr(3, 2) - D0, Fr(3, 2) + D0, 0)]
    ncov = 0
    while stack:
        a, b, d = stack.pop()
        db, dg = dlogBG(ball(a, b))
        if (db.rad() > 0.05 or dg.rad() > 0.05) and d < 12:
            m = (a + b) / 2
            stack += [(a, m, d + 1), (m, b, d + 1)]
            continue
        ncov += 1
        lam_b = max(lam_b, upper(abs(db)), key=float)
        lam_g = max(lam_g, upper(abs(dg)), key=float)
    # round up to 3 decimals
    LB = Fr(math.ceil(float(lam_b.mid() + lam_b.rad()) * 1000) + 1, 1000)
    LG = Fr(math.ceil(float(lam_g.mid() + lam_g.rad()) * 1000) + 1, 1000)
    ok &= bool(lam_b < _a(LB)) and bool(lam_g < _a(LG))
    lb, lg = cw.logBG(arb(3) / 2)
    Bs, Gs = lb.exp(), lg.exp()
    # closed form check: B* = (3/sqrt21) exp(Phi(sqrt21/2)), Phi(3/2) = 0
    alt = 3 / arb(21).sqrt() * cw.phipsi(arb(21).sqrt() / 2)[0].exp()
    ok &= bool(abs(alt - Bs) < arb("1e-35"))
    RESULT.update(theta0=th0, LAM_B=LB, LAM_G=LG, Bstar=Bs, Gstar=Gs, dcover=ncov)
    check(ok, f"CX.1: theta0 = (3-sqrt5)/2 = {th0.str(8, radius=False)}; on J = [3/2 -+ {D0}] ({ncov} cells) "
              f"|(log B)'| < {float(LB)}, |(log G)'| < {float(LG)} (enclosure max {lam_b.str(5, radius=False)}, "
              f"{lam_g.str(5, radius=False)}); B* = {Bs.str(20)}, G* = {Gs.str(20)}; "
              f"B* = (3/sqrt21) exp(Phi(sqrt21/2)) checked; Y = sigma_2(X) subset [{float(YLO):.6f}, 5/2] subset X'")
    return ok


def precompute():
    """Per cell of [YLO, 5/2]: enclosures of (log B, log G) at x_n(cell), n < N_H."""
    cells = []
    for i in range(NCELL):
        a = YLO + (Fr(5, 2) - YLO) * i / NCELL
        b = YLO + (Fr(5, 2) - YLO) * (i + 1) / NCELL
        x = ball(a, b)
        terms = []
        for _ in range(N_H):
            terms.append(cw.logBG(x))
            x = sig1(x)
        cells.append(terms)
    return cells


def u_bounds(cells, k, eps):
    """(max over cells of an upper bound, min over cells of a lower bound) for u_eps on Y."""
    Bs, Gs = RESULT["Bstar"], RESULT["Gstar"]
    spread = (_a(RESULT["LAM_G"]) + k * _a(RESULT["LAM_B"])) * RESULT["theta0"] ** N_H
    q = 1 / (1 + eps)
    tail = Gs * Bs ** k * q ** N_H / eps
    t_hi, t_lo = tail * spread.exp(), tail * (-spread).exp()
    hi = lo = None
    for terms in cells:
        s = arb(0)
        w = q
        for lb, lg in terms:
            s += w * (lg + k * lb).exp()
            w *= q
        h, l_ = upper(s + t_hi), lower(s + t_lo)
        hi = h if hi is None or bool(h > hi) else hi
        lo = l_ if lo is None or bool(l_ < lo) else lo
    return hi, lo


def cx2(cells):
    Bs, Gs = RESULT["Bstar"], RESULT["Gstar"]
    br = {}
    ok = True
    for k in range(K + 1):
        scale = Gs * Bs ** k
        # bisection in c on [0.5, 2] for the upper criterion (sup u <= 1) and lower (inf u >= 1)
        def up_ok(c):
            return bool(u_bounds(cells, k, scale * c)[0] <= 1)

        def lo_ok(c):
            return bool(u_bounds(cells, k, scale * c)[1] >= 1)
        a, b = Fr(1, 2), Fr(2)            # want smallest c with up_ok
        ok_k = up_ok(_a(b)) and lo_ok(_a(a))
        for _ in range(24):
            m = (a + b) / 2
            if up_ok(_a(m)):
                b = m
            else:
                a = m
        chi = b
        a, b = Fr(1, 2), Fr(2)            # want largest c with lo_ok
        for _ in range(24):
            m = (a + b) / 2
            if lo_ok(_a(m)):
                a = m
            else:
                b = m
        clo = a
        ok_k &= up_ok(_a(chi)) and lo_ok(_a(clo)) and clo <= chi
        ok &= ok_k
        br[k] = (clo, chi)
    RESULT["bracket"] = br
    rows = [f"k={k}: c_k in [{_dn(br[k][0])}, {_up(br[k][1])}]" for k in (0, 1, 2, 5, 12, 13, 20, 30, K - 1, K)]
    check(ok, f"CX.2: cover of Y by {NCELL} cells, N_H = {N_H}, renewal criterion (R) for every integer k = 0..{K}:\n      "
              + "\n      ".join(rows))
    return ok


def cx3() -> bool:
    """All inequalities in arb (no floats)."""
    th0 = RESULT["theta0"]
    LB, LG = _a(RESULT["LAM_B"]), _a(RESULT["LAM_G"])
    tau, dl = _a(TAU), _a(DELTA)
    bmax, gmax = _a(cw.BOUNDS["B_hi"]), _a(cw.BOUNDS["G_hi"])
    ok = True
    # N(k) = least N with (LG + k LB) theta0^N <= tau satisfies N(k) <= Nb(k) = 1 + log((LG + k LB)/tau)/log(1/theta0)
    # (LG + k LB >= tau); theta0^N(k) <= tau/(LG + K LB) <= D0, so x_n in J for n >= N(k)
    ok &= bool(LG + K * LB >= tau) and bool(tau / (LG + K * LB) <= _a(D0))
    nk = 1 + ((LG + K * LB) / tau).log() / (1 / th0).log()
    F = nk * gmax * bmax ** K
    # F(k) = Nb(k) gmax bmax^k decreases for k >= K: (log Nb)' <= 1/(k log(1/theta0)) < -log bmax
    ok &= bool(1 / (K * (1 / th0).log()) < -bmax.log())
    up = F + tau.exp() / (1 + dl)
    lo = (1 - F) * (-tau).exp() / (1 - dl)
    ok &= bool(up < 1) and bool(lo > 1)
    # monotonicity for real k >= K: eps_(k+1)/eps_k <= B* (1+DELTA)/(1-DELTA) < 1
    mono = RESULT["Bstar"] * (1 + dl) / (1 - dl)
    ok &= bool(mono < 1)
    RESULT.update(NK=nk, F=F, cx3_up=up, cx3_lo=lo, mono=mono)
    check(ok, f"CX.3: for every real k >= K = {K}: N(k) <= {upper(nk).str(5, radius=False)}, N(k) 0.2873 0.8352^k <= "
              f"{upper(F).str(4, radius=False)} (decreasing in k), sup_Y u <= {upper(up).str(7, radius=False)} < 1 at "
              f"eps = (1+{float(DELTA)}) G* B*^k and inf_Y u >= {lower(lo).str(7, radius=False)} > 1 at eps = "
              f"(1-{float(DELTA)}) G* B*^k (TAU = {float(TAU)}); so {float(1 - DELTA)} <= c_k <= {float(1 + DELTA)}, and "
              f"eps_(k+1)/eps_k <= B*(1+DELTA)/(1-DELTA) <= {upper(mono).str(5, radius=False)} < 1")
    return ok


def cx4_cx5() -> bool:
    br = RESULT["bracket"]
    Bs, Gs = RESULT["Bstar"], RESULT["Gstar"]
    ok = True
    # (b) integer k = 0..K-1: eps_(k+1) <= c_hi(k+1) G* B*^(k+1) < c_lo(k) G* B*^k <= eps_k;
    # k = K-1 -> K uses CX.2 at k = K (CX.2 covers k = K as well); k >= K by CX.3 (real k).
    worst = None
    for k in range(K):
        r = _a(br[k + 1][1]) * Bs / _a(br[k][0])
        ok &= bool(r < 1)
        worst = r if worst is None or bool(r > worst) else worst
    rows, cons = [], True
    for k, (m, rr) in cw.RHO_PRIOR.items():
        e = (arb(m) + arb(0, rr)) * arb(3) ** (k + 1) - 1
        c = e / (Gs * Bs ** k)
        inside = bool(c > _a(br[k][0])) and bool(c < _a(br[k][1]))
        cons &= inside
        if k in (0, 1, 6, 12):
            rows.append(f"k={k}: {_dn(br[k][0])} < {c.str(8, radius=False)} < {_up(br[k][1])}")
    RESULT.update(worst_ratio=worst, consistency=cons)
    check(ok, f"CX.4(b): eps_(k+1)/eps_k <= {worst.str(5, radius=False)} < 1 for k = 0..{K - 1} (CX.2), and for real "
              f"k >= {K} (CX.3): 3^k rho(L_k) strictly decreases for every integer k >= 0, to 1/3")
    check(cons, "CX.5 (consistency): CE's r and CT/CU/CV's rho(L_k), k <= 12, give c_k inside CX.2's brackets: "
                + "; ".join(rows))
    for k in (13, 20, 39):
        lo = (1 + _a(br[k][0]) * Gs * Bs ** k) / arb(3) ** (k + 1)
        hi = (1 + _a(br[k][1]) * Gs * Bs ** k) / arb(3) ** (k + 1)
        print(f"   CX.2 at k = {k}: rho(L_{k}) in [{lower(lo).str(12, radius=False)}, {upper(hi).str(12, radius=False)}]",
              flush=True)
    return ok and cons


def main() -> int:
    print("Theorem CX: (3^(k+1) rho(L_k) - 1)/(G* B*^k) -> 1 with B* = B(3/2) = 0.77737602763..., and 3^k rho(L_k) "
          "strictly decreasing for every k. Not thrust.", flush=True)
    t0 = time.time()
    FAILS.clear()
    RESULT.clear()
    old = ctx.prec
    ctx.prec = PREC
    try:
        cw.FAILS.clear()
        ok = cw.cw1()             # CW's tails (imported and run unchanged)
        ok = cw.cw2() and ok      # CW.2's bounds 0.2873, 0.8352 (used in CX.3)
        ok = cx1() and ok
        cells = precompute()
        ok = cx2(cells) and ok
        ok = cx3() and ok
        ok = cx4_cx5() and ok
        check(ok and not cw.FAILS, f"CX.4(a): (3^(k+1) rho(L_k) - 1)/(G* B*^k) -> 1: exact rate B* = "
                                   f"{RESULT['Bstar'].str(15, radius=False)}, constant G* = {RESULT['Gstar'].str(15, radius=False)}")
    finally:
        ctx.prec = old
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CX checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
