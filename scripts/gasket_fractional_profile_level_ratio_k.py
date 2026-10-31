#!/usr/bin/env python3
"""Theorem CY: q^(k) = nu_2^(k)/rho(L_k) -> -1/3 (CV's and CX's HEURISTIC), with an explicit geometric rate, and a
uniform-in-k spectral gap: for every integer k >= K0 the operator L_k has exactly two eigenvalues outside a disk of
radius r_in 3^(-k-1), both real and simple.

Setting exactly as in Theorems CT-CX (GASKET.md): sigma_(1,2)(s) = sqrt(15/4 -+ s), u_i = wt_i/(2 sigma_i)^k,
wt_1 = (s+1/2)/(2s(7/2-s)), wt_2 = (s-1/2)/(2s(s+7/2)), (L_k f) = u_1 f o sigma_1 + u_2 f o sigma_2, s* = 3/2,
D_CG = {|s - 83/50| < 8/5} (CG's disk, on whose H^2 CT-CV certified rho(L_k) and nu_2^(k) for k <= 12).

Construction (exact algebra).  On the disk D = {|s - c| < R} (c = 9/5, R = 4/5 below) put x_j = sigma_1^j(s) and the
convergent products (no logarithms)
    P0 = prod_j 3 wt_1(x_j),  P1 = prod_j 3/(2 sigma_1(x_j)),  kappa = (s - s*) prod_j 3/(sigma_1(x_j) + 3/2),
so P0 o sigma_1 = P0/(3 wt_1), P1 o sigma_1 = P1 (2 sigma_1/3), kappa o sigma_1 = -kappa/3, kappa(s*) = 0, kappa'(s*) = 1.
h_k = P0 P1^k (CW's exp(k Phi + Psi), continued holomorphically) gives EXACTLY
    Lt_k := 3^(k+1) h_k^(-1) L_k h_k = C_1 + M_a C_2,   C_i f = f o sigma_i,   a = G B^k,
    B = (3/(2 sigma_2)) P1(sigma_2 s)/P1(s),  G = 3 wt_2 P0(sigma_2 s)/P0(s)    (CW's B, G on X).
C_1 on H^oo(D) splits along f = f(s*) 1 + f'(s*) kappa + g, g in V_2 = {g(s*) = g'(s*) = 0}: C_1 = 1 (+) (-1/3) (+) C_1|V_2,
and by Schwarz's lemma ||C_1^n g|| <= c_n ||g||, c_n = sup_D |phi(sigma_1^n s)|^2 (phi: D -> unit disk Moebius, phi(s*) = 0).
Hence for |z| >= r_in, |z - 1| >= r_1, |z + 1/3| >= r_2:
    ||(z - C_1)^(-1)|| <= 1/r_1 + beta/r_2 + Q S(r_in),  beta = ||kappa|| |phi'(s*)|,  Q = 2 + beta,  S(r) = sum c_n r^(-n-1).
With ||M_a C_2|| <= eta_k = sup_D|G| (sup_D|B|)^k, a Neumann series plus Kato's homotopy of Riesz projections
(C_1 + t M_a C_2, t in [0, 1]) give the eigenvalue statements below.

Theorem CY (checked here; exit 1 on any failure).
  CY.1 [proved; exact + arb] Disk facts: poles/zeros of wt_i and the branch cuts of sigma_i stay off the closed disk
       D' = {|s - c| <= R + 1/100}; sigma_1(D'), sigma_2(D') are inside D (disk-image lemma); every sigma_1-orbit from D'
       enters {|x - s*| <= 1/4} and then contracts by 7/20 per step, so P0, P1, kappa converge uniformly on D' to
       holomorphic, zero-free (P0, P1) functions; arb over NARC arcs of dD gives sup_D|B| < B_M, sup_D|G| < G_M,
       sup_D|kappa| < K_M, and c_1..c_NC (tails explicit).
  CY.2 [proved from CY.1] the resolvent bound above with r_1 = 1/2 and r_in = RIN; A = Q S(RIN) enclosed.
  CY.3 [proved from CY.1-CY.2, every integer k >= K0] eta_k (2 + A) < 1 and r_2(k) = 2 beta eta_k/(1 - eta_k(2+A)) < 1/3 - RIN.
       Then on H^oo(D), Lt_k has exactly one eigenvalue in |z - 1| < 1/2, one in |z + 1/3| < r_2(k), none else in
       |z| >= RIN; both are real and simple.  The first is 3^(k+1) rho(L_k) (rho on C(X), as in CW/CX).  Call the second
       3^(k+1) nu_2^(k).  So |3^(k+1) nu_2^(k) + 1/3| <= r_2(k) <= C2 B_M^k and, with CW.3 (used as stated),
           |q^(k) + 1/3| <= r_2(k) + (0.2873/3) 0.8352^k  ->  0:   q^(k) -> -1/3 (geometrically, rate <= B_M).
       Every other eigenvalue nu of L_k has |nu| < RIN 3^(-k-1) < RIN rho(L_k): a uniform-in-k gap |nu/nu_2| < RIN/(1/3 - r_2).
  CY.4 [proved; exact + arb] Same nonzero spectrum on CT-CV's space: closure(D) is inside D_CG and sigma_w(closure D_CG)
       is inside D for every word w of length M_WORD, so restriction R: H^2(D_CG) -> H^oo(D) and extension
       E = L_k^M_WORD: H^oo(D) -> H^2(D_CG) intertwine L_k and are injective on generalized eigenspaces of nonzero
       eigenvalues: nonzero eigenvalues and algebraic multiplicities agree.  So for k >= K0 the nu_2^(k) of CY.3 is CT-CV's
       nu_2 (second eigenvalue of L_k on H^2(D_CG)), and CY.3 is CV.3's statement for every k >= K0 with rho_in = RIN 3^(-k-1).
  CY.5 [certified; CU.1-CU.3's proofs, CU's routine imported and run unchanged] for k = 13..23 (= K0 - 1), at N = 180,
       Neumann power 40, rho_in = 1/round(1/(0.7|nu_2|)): L_k on H^2(D_CG) has exactly two eigenvalues in |z| > rho_in,
       real and simple, rho(L_k) and nu_2^(k), each to +- 6e-20 (absolute).
  CY.6 [proved from CY.5 and CV.5 / CX.2 as stated] CX.2's brackets of rho(L_13), rho(L_20) contain CY.5's values;
       |q^(k)| strictly decreases for k = 8..23 (CV.6's decreasing half extended) and |q^(k)| > 1/3 for k = 8..23.
  Together with CG (k = 0), CT (k = 1), CU (2..4), CV (5..12): for EVERY integer k >= 0, nu_2^(k) is a real, simple,
  negative eigenvalue of L_k on H^2(D_CG), isolated from the rest of the spectrum, and q^(k) -> -1/3.
Kept limits: the rate B_M ~ 0.853 in CY.3 is a crude sup over dD, not the exact rate of q^(k) + 1/3 (floats: ratios of
successive q^(k) + 1/3 ~ 0.82 at k = 23); the sign of q^(k) + 1/3 (i.e. |q^(k)| > 1/3) and the monotonicity of |q^(k)| are
certified only for k <= 23; no uniform-in-k rate for inf R_n^(k) is derived here (the gap of CY.3 is uniform, the Riesz
projection constants are not tracked).
Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction as Fr

from flint import acb, arb, ctx

FAILS: list[str] = []
RESULT: dict = {}
PREC = 160
C0 = Fr(9, 5)          # disk centre
R0 = Fr(4, 5)          # disk radius
RP = R0 + Fr(1, 100)   # D' radius (holomorphy on a neighbourhood of the closed disk)
SST = Fr(3, 2)         # s*
J = 60                 # orbit terms kept in P0, P1, kappa
NC = 30                # c_n enclosed for n <= NC
NARC = 4096            # arcs of dD
RIN = Fr(9, 40)        # inner radius r_in
R1 = Fr(1, 2)          # radius around the eigenvalue 1
THETA = Fr(7, 20)      # contraction of sigma_1 on {|x - s*| <= 1/4}
DLT = Fr(1, 4)
D_CG = (Fr(83, 50), Fr(8, 5))
M_WORD = 2
K_MAX_REPORT = 60
CW_EPS = (Fr(2873, 10000), Fr(8352, 10000))   # CW.3 (used as stated): 3^(k+1) rho(L_k) - 1 <= 0.2873 * 0.8352^k
# CY.5: CU's routine (imported unchanged) for k = 13..24 at N = 180, Neumann power 40, rho_in = 1/round(1/(0.7 |nu_2|)).
# Float seeds (rho, nu_2, nu_3) from a 60-node Chebyshev collocation (floats; Krawczyk certifies rho and nu_2; nu_3 unused
# except to place rho_in, whose validity is certified by the Neumann bound).
CU_KS = [int(x) for x in __import__("os").environ.get("CY_KS", ",".join(str(k) for k in range(13, 24))).split(",") if x]
N_CU, M_POW_CU, RIN_FRAC = 180, 40, 0.7
CU_SEEDS = {
    13: ('2.113436488714e-7', '-7.39580416187e-8', '2.477357595204e-8'),
    14: ('7.028021203214e-8', '-2.443989881658e-8', '8.305626083679e-9'),
    15: ('2.338318502613e-8', '-8.083884325709e-9', '2.773427461748e-9'),
    16: ('7.783090016327e-9', '-2.676375244504e-9', '9.234903616131e-10'),
    17: ('2.591430403713e-9', '-8.868751475344e-10', '3.06907269382e-10'),
    18: ('8.630495224557e-10', '-2.941293187223e-10', '1.018676585125e-10'),
    19: ('2.874859786867e-10', '-9.762039934433e-11', '3.378647297635e-11'),
    20: ('9.577754576509e-11', '-3.242150100568e-11', '1.120194379002e-11'),
    21: ('3.191260171016e-11', '-1.077407333413e-11', '3.713727648306e-12'),
    22: ('1.06341011446e-11', '-3.58217438757e-12', '1.231345372578e-12'),
    23: ('3.543810888774e-12', '-1.191521485931e-12', '4.083780197655e-13'),
}
# CV.5's certified |q^(k)| (10 digits, used as stated, radius 1e-10) for the monotonicity check CY.6
QABS_CV = {8: "0.3588353292", 9: "0.3581148072", 10: "0.3565328657", 11: "0.3544800370", 12: "0.3522236835"}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _a(x: Fr) -> arb:
    return arb(x.numerator) / arb(x.denominator)


def upper(x: arb) -> arb:
    return arb(x.mid() + x.rad())


def lower(x: arb) -> arb:
    return arb(x.mid() - x.rad())


def up_str(x: arb, sig: int = 6) -> str:
    """A decimal string >= every point of the ball x (rounded UP to sig significant digits; for printing only)."""
    v = float(upper(x)) * (1 + 1e-12) + 1e-300
    if v <= 0:
        return f"{v:.{sig}g}"
    e = math.floor(math.log10(v)) - sig + 1
    return f"{math.ceil(v / 10 ** e) * 10 ** e:.{sig}g}"


def up_fr(x: arb, den: int = 10 ** 6) -> Fr:
    """A rational >= every point of the ball x (rounded up on a grid of 1/den, plus one step)."""
    return Fr(math.ceil(float(upper(x)) * den) + 1, den)


def sig1(x):
    return (15 * acb(1) / 4 - x).sqrt()


def sig2(x):
    return (15 * acb(1) / 4 + x).sqrt()


def wt1(x):
    return (x + acb(1) / 2) / (2 * x * (acb(7) / 2 - x))


def wt2(x):
    return (x - acb(1) / 2) / (2 * x * (x + acb(7) / 2))


# ----------------------------------------------------------------------------------------------- disk-image lemma
def disk_image(c: Fr, r: Fr, sign: int):
    """sigma(s) = sqrt(15/4 + sign*s) on the closed disk D(c, r), c real, needs 15/4 + sign*c - r > 0 (exact).
    |sigma(s) - sigma(c)| = |s - c|/|sigma(s) + sigma(c)| and Re sigma(s) >= sqrt(Re(15/4 + sign*s)) >= sqrt(15/4+sign*c-r),
    so sigma(D(c, r)) lies in D(sigma(c), r/(sqrt(15/4 + sign c - r) + sigma(c))).  Returns (centre, radius) as arb (upper)."""
    w = Fr(15, 4) + sign * c
    assert w - r > 0
    cen = _a(w).sqrt()
    rad = _a(r) / (_a(w - r).sqrt() + cen)
    return cen, rad


def inside(cen: arb, rad: arb, c: Fr, r: Fr) -> bool:
    """closed D(cen, rad) inside open D(c, r) (cen real)."""
    return bool(abs(cen - _a(c)) + rad < _a(r))


def cy1_disk() -> bool:
    ok = True
    # poles 0, 7/2, -7/2 of wt_i, zero -1/2 of wt_1 off D'; branch cuts: Re(15/4 -+ s) > 0 on D'
    ok &= C0 - RP > 0 and C0 + RP < Fr(7, 2) and Fr(15, 4) - C0 - RP > 0 and Fr(15, 4) + C0 - RP > 0
    ok &= abs(C0 - SST) < R0
    imgs = []
    for sg in (-1, 1):
        cen, rad = disk_image(C0, RP, sg)
        imgs.append((cen, rad))
        ok &= inside(cen, rad, C0, R0)
    # orbit entry: |x - s*| <= d  =>  |sigma_1 x - s*| <= d/(sqrt(9/4 - d) + 3/2)  (needs d < 9/4)
    d = _a(abs(C0 - SST) + RP)
    steps = 0
    while not bool(d <= _a(DLT)) and steps < 20:
        ok &= bool(d < _a(Fr(9, 4)))
        d = d / ((_a(Fr(9, 4)) - d).sqrt() + _a(SST))
        d = upper(d)
        steps += 1
    ok &= bool(d <= _a(DLT))
    # on |x - s*| <= 1/4: sqrt(9/4 - 1/4) + 3/2 = sqrt2 + 3/2 >= 7/5 + 3/2 = 29/10 > 20/7  =>  factor <= 7/20
    ok &= Fr(7, 5) ** 2 <= 2 and Fr(29, 10) * THETA >= 1
    # factor bounds on |x - s*| <= 1/4 (exact rational inequalities):
    #   3 wt_1(x) - 1 = (x - s*)(x - 1/2)/(x (7/2 - x)):        |.| <= (5/4)/((5/4)(7/4)) |x - s*| = (4/7)|x - s*|
    #   3/(2 sigma_1 x) - 1 = 4(x - s*)/(2 sigma_1 (3 + 2 sigma_1)), |sigma_1| >= sqrt2 >= 7/5: |.| <= (4/16.24)|x-s*| <= |x-s*|/4
    #   3/(sigma_1 x + 3/2) - 1 = (x - s*)/(sigma_1 x + 3/2)^2:  |.| <= |x - s*|/(29/10)^2 <= |x - s*|/8
    ok &= Fr(5, 4) / (Fr(5, 4) * Fr(7, 4)) == Fr(4, 7) and 4 / (2 * Fr(7, 5) * (3 + 2 * Fr(7, 5))) <= Fr(1, 4)
    ok &= Fr(29, 10) ** 2 >= 8
    RESULT.update(img1=imgs[0], img2=imgs[1], entry_steps=steps)
    check(ok, f"CY.1(a): D = D({C0}, {R0}), D' = D({C0}, {RP}): wt_i poles/zeros and sigma_i branch cuts off D'; "
              f"sigma_1(D') in D({imgs[0][0].str(5, radius=False)}, {up_str(imgs[0][1], 4)}), "
              f"sigma_2(D') in D({imgs[1][0].str(5, radius=False)}, {up_str(imgs[1][1], 4)}), both inside D; "
              f"sigma_1-orbits from D' reach |x - 3/2| <= 1/4 after {steps} steps and contract by {THETA} there")
    return ok


# ----------------------------------------------------------------------------------------------- products on arcs
FACT = {"P0": Fr(4, 7), "P1": Fr(1, 4), "K": Fr(1, 8)}


def tail_ball(dJ: arb, cfac: Fr) -> acb:
    """prod_(j >= J) f_j with |f_j - 1| <= cfac |x_j - s*|, |x_J - s*| <= dJ <= 1/4, contraction THETA:
    sum_j |f_j - 1| <= cfac dJ/(1 - THETA), each |f_j - 1| <= cfac/4 <= 1/7, |log f_j| <= (7/6)|f_j - 1|;
    so the tail lies in 1 + disk(e^t - 1), t = (7/6) cfac dJ/(1 - THETA)."""
    t = _a(Fr(7, 6) * cfac / (1 - THETA)) * dJ
    e = upper(t.exp() - 1)
    return acb(1) + acb(arb(0, e), arb(0, e))


def orbit_products(s: acb):
    """P0(s), P1(s), kappa(s) enclosures (J terms + tails), and the orbit x_0..x_NC."""
    p0 = acb(1)
    p1 = acb(1)
    kp = s - acb(3) / 2
    x = s
    xs = [s]
    for j in range(J):
        y = sig1(x)
        p0 *= 3 * wt1(x)
        p1 *= 3 / (2 * y)
        kp *= 3 / (y + acb(3) / 2)
        x = y
        if j + 1 <= NC:
            xs.append(x)
    dJ = upper(abs(x - acb(3) / 2))
    if not bool(dJ <= _a(DLT)):
        raise RuntimeError("orbit did not enter |x - s*| <= 1/4")
    return p0 * tail_ball(dJ, FACT["P0"]), p1 * tail_ball(dJ, FACT["P1"]), kp * tail_ball(dJ, FACT["K"]), xs


def arc(i: int, n: int) -> acb:
    """A complex ball containing the i-th of n equal arcs of |s - C0| = R0."""
    pi = arb.pi()
    tm = (2 * i + 1) * pi / n
    rad = _a(R0) * pi / n                      # half arc length >= distance from midpoint to any arc point
    z = acb(_a(C0)) + _a(R0) * acb(tm.cos(), tm.sin())
    return z + acb(arb(0, rad), arb(0, rad))


def phi(w: acb) -> acb:
    """Moebius map of D onto the unit disk with phi(s*) = 0 (s* - c real)."""
    return _a(R0) * (w - acb(3) / 2) / (_a(R0) ** 2 - _a(SST - C0) * (w - _a(C0)))


def cy1_sups() -> bool:
    bm = gm = km = arb(0)
    cmax = [arb(0)] * (NC + 1)
    dNC = arb(0)
    for i in range(NARC):
        s = arc(i, NARC)
        p0, p1, kp, xs = orbit_products(s)
        y = sig2(s)
        q0, q1, _, _ = orbit_products(y)
        B = 3 / (2 * y) * q1 / p1
        G = 3 * wt2(s) * q0 / p0
        bm = max(bm, upper(abs(B)), key=float)
        gm = max(gm, upper(abs(G)), key=float)
        km = max(km, upper(abs(kp)), key=float)
        for n in range(1, NC + 1):
            v = upper(abs(phi(xs[n])))
            cmax[n] = max(cmax[n], upper(v * v), key=float)
        dNC = max(dNC, upper(abs(xs[NC] - acb(3) / 2)), key=float)
    BM, GM, KM = up_fr(bm, 10 ** 4), up_fr(gm, 10 ** 4), up_fr(km, 10 ** 4)
    cn = [arb(1)] + [arb(1) if bool(cmax[n] >= 1) else cmax[n] for n in range(1, NC + 1)]
    e = _a(abs(SST - C0))
    dphi = _a(R0) / (_a(R0) ** 2 - e * e)
    ok = bool(bm < _a(BM)) and bool(gm < _a(GM)) and bool(km < _a(KM)) and BM < 1 and bool(dNC <= _a(DLT))
    ok &= all(bool(cn[n] < 1) for n in range(1, NC + 1))
    RESULT.update(BM=BM, GM=GM, KM=KM, Bsup=bm, Gsup=gm, Ksup=km, cn=cn, dNC=dNC, dphi=dphi)
    check(ok, f"CY.1(b): {NARC} arcs of dD, J = {J} terms + tails: sup|B| <= {up_str(bm, 6)} < {float(BM)}, "
              f"sup|G| <= {up_str(gm, 6)} < {float(GM)}, sup|kappa| <= {up_str(km, 6)} < "
              f"{float(KM)}; |phi'(s*)| = {dphi.str(8, radius=False)}; c_1..c_4 <= "
              + ", ".join(up_str(cn[n], 4) for n in range(1, 5))
              + f"; sup|sigma_1^{NC} - s*| <= {up_str(dNC, 3)}")
    return ok


def cy2() -> bool:
    cn, dNC, dphi = RESULT["cn"], RESULT["dNC"], RESULT["dphi"]
    beta = _a(RESULT["KM"]) * dphi
    Q = 2 + beta
    r = _a(RIN)
    S = arb(0)
    for n in range(NC + 1):
        S += cn[n] / r ** (n + 1)
    # tail n > NC: |x_n - s*| <= dNC THETA^(n-NC), |phi(w)| <= R |w - s*|/(R^2 - e(e + |w - s*|))
    e = _a(abs(SST - C0))
    kq = (_a(R0) * dNC / (_a(R0) ** 2 - e * (e + dNC))) ** 2
    th2 = _a(THETA) ** 2
    ok = bool(th2 < r) and bool(_a(R0) ** 2 - e * (e + dNC) > 0)
    S += kq * r ** (-NC - 1) * (th2 / r) / (1 - th2 / r)
    A = Q * S
    RESULT.update(beta=beta, Q=Q, S=S, A=A)
    check(ok, f"CY.2: on {{|z| >= {RIN}, |z - 1| >= {R1}, |z + 1/3| >= r_2}}: ||(z - C_1)^(-1)|| <= 1/r_1 + beta/r_2 + Q S(r_in), "
              f"beta = ||kappa|| |phi'(s*)| <= {up_str(beta, 6)}, Q <= {up_str(Q, 6)}, "
              f"S({RIN}) <= {up_str(S, 6)}, A = Q S <= {up_str(A, 6)}")
    return ok


def eta(k: int) -> arb:
    return _a(RESULT["GM"]) * _a(RESULT["BM"]) ** k


def r2(k: int):
    et = eta(k)
    den = 1 - et * (2 + RESULT["A"])
    return et, den, 2 * RESULT["beta"] * et / den


def cy3() -> bool:
    K0 = None
    for k in range(0, 400):
        et, den, rr = r2(k)
        if bool(den > 0) and bool(rr < _a(Fr(1, 3) - RIN)):
            K0 = k
            break
    ok = K0 is not None
    if not ok:
        check(False, "CY.3: no threshold K0 found")
        return False
    # monotone: eta_k decreasing (B_M < 1), so den increasing and r_2(k) decreasing in k: the conditions hold for all k >= K0
    et0, den0, rr0 = r2(K0)
    C2 = 2 * RESULT["beta"] * _a(RESULT["GM"]) / den0
    # disks: D(1, 1/2) and D(-1/3, r_2) disjoint and inside |z| > r_in
    ok &= bool(rr0 < _a(Fr(1, 3) - RIN)) and R1 > RIN and Fr(1, 2) + Fr(1, 3) > 0
    # the eigenvalue near 1 is 3^(k+1) rho(L_k) and lies in [1, 1 + 0.2873 * 0.8352^k] (CW.3) -- consistency with |z - 1| < 1/2
    ok &= bool(_a(CW_EPS[0]) < _a(R1))
    rows = []
    qb = {}
    for k in sorted({K0, K0 + 1, 25, 30, 40, 50, K_MAX_REPORT}):
        if k < K0:
            continue
        _, _, rr = r2(k)
        epsk = _a(CW_EPS[0]) * _a(CW_EPS[1]) ** k
        qq = upper(rr + epsk / 3)
        qb[k] = (upper(rr), qq)
        rows.append(f"k={k}: |3^(k+1) nu_2 + 1/3| <= {up_str(rr, 4)}, |q + 1/3| <= {up_str(qq, 4)}")
    gap = _a(RIN) / (_a(Fr(1, 3)) - rr0)
    RESULT.update(K0=K0, C2=C2, qbound=qb, gap=gap, r2K0=rr0)
    check(ok, f"CY.3: threshold K0 = {K0} (eta_K0 (2 + A) = {up_str(et0 * (2 + RESULT['A']), 4)} < 1, "
              f"r_2(K0) <= {up_str(rr0, 5)} < 1/3 - {RIN}); for every integer k >= {K0}: exactly two "
              f"eigenvalues of Lt_k in |z| >= {RIN}, real, simple, 3^(k+1) rho(L_k) in |z - 1| < 1/2 and 3^(k+1) nu_2 in "
              f"|z + 1/3| <= r_2(k) <= {up_str(C2, 5)} * {float(RESULT['BM'])}^k; uniform gap |nu/nu_2| < "
              f"{up_str(gap, 4)}:\n      " + "\n      ".join(rows))
    return ok


def cy4() -> bool:
    c, r = D_CG
    ok = abs(C0 - c) + R0 < r            # closure(D) inside D_CG (restriction bounded, injective)
    disks = [(_a(c), _a(r))]
    worst = arb(-1)
    for _ in range(M_WORD):
        new = []
        for cen, rad in disks:
            for sg in (-1, 1):
                # disk_image with arb centre: use an exact rational enclosure of the (real) centre and radius
                cl, ch = Fr(math.floor(float(lower(cen)) * 10 ** 12), 10 ** 12), Fr(math.ceil(float(upper(cen)) * 10 ** 12), 10 ** 12)
                rr = up_fr(rad, 10 ** 12) + (ch - cl)
                cm = (cl + ch) / 2
                ncen, nrad = disk_image(cm, rr, sg)
                new.append((ncen, upper(nrad)))
        disks = new
    for cen, rad in disks:
        ok &= inside(cen, rad, C0, R0)
        worst = max(worst, upper(abs(cen - _a(C0)) + rad), key=float)
    RESULT.update(word_worst=worst)
    check(ok, f"CY.4: closure(D) inside D_CG = D({c}, {r}); sigma_w(closure D_CG) inside D for all {2 ** M_WORD} words of length "
              f"{M_WORD} (max |centre - {C0}| + radius = {up_str(worst, 5)} < {R0}); so L_k has the same nonzero "
              f"eigenvalues, with multiplicities, on H^2(D_CG) and H^oo(D): CY.3's nu_2^(k) is CT-CV's nu_2")
    return ok


def cy5() -> bool:
    """CU.1-CU.3 (CU's routine, imported and run unchanged) for k in CU_KS: exactly two eigenvalues of L_k on H^2(D_CG)
    in |z| > rho_in, real and simple, rho(L_k) and nu_2^(k)."""
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import gasket_fractional_profile_level_infrate_k as cu    # noqa: E402
    old = (cu.N, cu.cr.M_POW)
    cu.N, cu.cr.M_POW = N_CU, M_POW_CU
    out = {}
    ok = True
    try:
        for k in CU_KS:
            r, nu, _ = CU_SEEDS[k]
            rf, nf = float(r), abs(float(nu))
            rin = Fr(1, round(1 / (RIN_FRAC * nf)))
            cu.CONF[k] = dict(seeds=[r, nu], rho_in=rin, proj=[Fr(1, round(3 / rf)), Fr(1, round(2 / (nf - float(rin))))],
                              cl3=None, qcap="0.36")
            cu.configure(k)
            cu.RESULT.clear()
            nf0 = len(cu.FAILS) + len(cu.cr.FAILS)
            t0 = time.time()
            print(f"--- CY.5, k = {k}: CU.1-CU.3 at N = {N_CU}, m = {M_POW_CU}, rho_in = {rin}", flush=True)
            dc = cu.disk_constants()
            eps, beta = cu.tail_constants(dc, cu.N)
            M = cu.build_matrix(cu.N)
            sp = cu.spectral_split(M)
            good = sp is not None
            if good:
                pairs, T, alpha = sp
                _, rho, nuv, q, good = cu.certify(pairs, T, alpha, eps, beta)
            good = bool(good) and len(cu.FAILS) + len(cu.cr.FAILS) == nf0
            ctx.prec = PREC
            if good:
                out[k] = dict(rho=rho, nu=nuv, q=q, rin=rin)
            check(good, f"CY.5 (k = {k}, {time.time() - t0:.0f}s): L_k on H^2(D_CG) has exactly two eigenvalues in |z| > {rin}, "
                        f"real and simple: rho = {rho.str(13) if good else '?'}, nu_2 = {nuv.str(13) if good else '?'}, "
                        f"q = {q.str(12) if good else '?'}")
            ok &= good
            del cu.CONF[k]
    finally:
        cu.N, cu.cr.M_POW = old
        ctx.prec = PREC
    RESULT["cu"] = out
    return ok


def cy6() -> bool:
    """Consistency with CY.3 / CX.2 and the monotonicity of |q^(k)| beyond CV.6."""
    out = RESULT["cu"]
    K0 = RESULT["K0"]
    ok = True
    rows = []
    for k in sorted(out):
        q = out[k]["q"]
        if k >= K0:
            qb = RESULT["qbound"].get(k)
            if qb is None:
                _, _, rr = r2(k)
                qb = (rr, upper(rr + _a(CW_EPS[0]) * _a(CW_EPS[1]) ** k / 3))
            inside_ = bool(abs(q + arb(1) / 3) < qb[1])
            ok &= inside_
            rows.append(f"k={k}: |q + 1/3| = {abs(q + arb(1) / 3).str(6)} <= CY.3's {up_str(qb[1], 4)}")
    aq = {k: arb(v) + arb(0, 1e-10) for k, v in QABS_CV.items()}
    for k in out:
        aq[k] = abs(out[k]["q"])
    ks = sorted(aq)
    dec = all(bool(aq[ks[i + 1]] < aq[ks[i]]) for i in range(len(ks) - 1)) and ks == list(range(ks[0], ks[-1] + 1))
    above = all(bool(aq[k] > arb(1) / 3) for k in ks)
    ok &= dec and above
    # CX.2's certified brackets (used as stated): rho(L_13) in [2.113399, 2.113510]e-7, rho(L_20) in [9.577739, 9.577786]e-11
    cxb = {13: ("2.113399e-7", "2.113510e-7"), 20: ("9.577739e-11", "9.577786e-11")}
    for k, (lo, hi) in cxb.items():
        if k in out:
            ok &= bool(out[k]["rho"] > arb(lo)) and bool(out[k]["rho"] < arb(hi))
            rows.append(f"rho(L_{k}) = {out[k]['rho'].str(10, radius=False)} in CX.2's [{lo}, {hi}]")
    RESULT.update(aq=aq, dec=dec, above=above)
    check(ok, f"CY.6: consistency ({'; '.join(rows) if rows else 'none'}); with CV.5's |q^(k)| "
              f"(k = 8..12, as stated): |q^(k)| strictly decreases for k = {ks[0]}..{ks[-1]} and stays > 1/3: "
              + ", ".join(f"{k}:{aq[k].str(8, radius=False)}" for k in ks))
    return ok


def main() -> int:
    print("Theorem CY: q^(k) = nu_2^(k)/rho(L_k) -> -1/3 with a geometric rate, and a uniform-in-k spectral gap. "
          "Not thrust.", flush=True)
    t0 = time.time()
    FAILS.clear()
    RESULT.clear()
    old = ctx.prec
    ctx.prec = PREC
    try:
        ok = cy1_disk()
        ok = cy1_sups() and ok
        ok = cy2() and ok
        ok = cy3() and ok
        ok = cy4() and ok
        if CU_KS:
            ok = cy5() and ok
            ok = cy6() and ok
        if ok:
            check(True, f"CY.3 (limit): q^(k) -> -1/3; |q^(k) + 1/3| <= {up_str(RESULT['C2'], 5)} * "
                        f"{float(RESULT['BM'])}^k + {float(CW_EPS[0])}/3 * {float(CW_EPS[1])}^k for every k >= {RESULT['K0']}")
    finally:
        ctx.prec = old
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CY checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
