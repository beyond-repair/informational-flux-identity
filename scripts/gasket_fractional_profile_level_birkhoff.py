#!/usr/bin/env python3
r"""Theorem DC: CF.3's Birkhoff-route constant made explicit, for a general log-Lipschitz f and for the level sums.

Setting (GASKET.md, Theorems CE, CF, CN, CO, CS): in s = sqrt(t + 25/4) in X_s = [sqrt5/2, 5/2],
(Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s), sigma_(1,2) = sqrt(15/4 -+ s), wt1 = (s + 1/2)/(2s(7/2 - s)),
wt2 = (s - 1/2)/(2s(s + 7/2)); C_b = {f > 0 : log f is b-Lipschitz in s}; a = 47/2, theta = tanh(Delta/4) < 0.5442 (CF.2);
M_1 = 14/15, M_2 < 0.5071, c_1 = 1/sqrt5, c_2 < 0.2267, phi(b) = max(M_1 + c_1 b, M_2 + c_2 b), b* = M_1/(1 - c_1)
(CF.2: Lt C_b c C_phi(b)).  CF.3: for f in C_b, b < a, |log(Lt^(m+1) f/Lt^m f)(x) - log r| <= theta^floor(m/5) Theta_a(f, Lt f),
and for the level sums |log(s_(k+1)/(r s_k))| <= C_G theta^(k/5) with C_G "depending on the Lipschitz constant of
log Lt^4 G, not made explicit".  CS made the *sharp* C_G* = 0.02105127810 explicit by the spectral route and kept the
Birkhoff-route constant Theta_a(f_0, Lt f_0), f_0 = Lt^4 G, as not computed.  This script computes it.

Theorem DC (checked here; exit 1 on any failure):
  DC.1 [proved; exact + arb] General f.  (i) phi(b) <= b* for b <= b*, and phi(b) - b* <= c_1 (b - b*) for b >= b*
       (because M_2 + c_2 b* < b* and M_2 < M_1, c_2 < c_1), so Lt^j C_b c C_(beta_j), beta_j = b* + c_1^j (b - b*)_+.
       (ii) For u, v in C_B, B < a, CE.3's cone lemma (lambda = B/a) gives
            Theta_a(u, v) <= 2 log((a + B)/(a - B)) + log(sup(v/u)/inf(v/u)) <= 2 log((a + B)/(a - B)) + osc log u + osc log v.
       Hence, with |X_s| = (5 - sqrt5)/2, for every f > 0 whose log is b-Lipschitz in s, every m >= j0(b) and x in X:
            |log(Lt^(m+1) f/Lt^m f)(x) - log r| <= theta^floor((m - j0)/5) K(beta),
            K(beta) = 2 log((a + B)/(a - B)) + (beta + phi(beta)) |X_s|,  B = max(beta, phi(beta)),
       where j0(b) = 0, beta = b if b < a, and otherwise j0 = min{j : b* + c_1^j (b - b*) < a}, beta = b* + c_1^j0 (b - b*).
       The same holds with any j >= j0 for which beta_j < a.  K is increasing in beta; K(0) <= 1.44879 bounds CF.4's
       Theta_0 = 0.765085 (consistency), K(b*) <= 4.95456, K(23) <= 56.356; b = 100 needs j0 = 2 (K <= 50.068), b = 1e4
       needs j0 = 8 (K <= 40.577).
  DC.2 [certified; acb/arb] Level sums.  CO.1/CO.2's certified data for F0 = Lt^4 G (re-run unchanged): S_rho, and the
       Taylor coefficients f_0..f_80 of F0 in u = (s - c)/R with aliasing bounds; the tail n > 80 by Cauchy, |f_n| <= S_rho rho^-n.
       On a cover of X_s by cells (interval Horner): F0 > 0, b_0 = sup |F0'/F0| (d/ds), osc log F0, osc log Lt F0, and
       sup/inf of R_0 = Lt F0/F0, all certified.  Since b_0 < a, f_0 in C_(b_0), and with B_0 = max(b_0, phi(b_0)),
            Theta_a(f_0, Lt f_0) <= Theta_G := 2 log((a + B_0)/(a - B_0)) + log(sup R_0/inf R_0).
       Values (16384 cells): F0 in [0.0029654816, 0.0048349722], b_0 <= 0.77333706, B_0 = phi(b_0) <= 1.2791803,
       R_0 in [0.4079370173, 0.4347214550], 0.0633864 <= log(sup R_0/inf R_0) <= 0.063592701, Theta_G <= 0.28154093
       (CF.3's osc + osc form: 1.2591208).
  DC.3 [proved from DC.2, CF.1-CF.3] For every k >= 4,
            |log(s_(k+1)/(r s_k))| <= theta^floor((k - 4)/5) Theta_G <= C_G^B theta^(k/5),  C_G^B = Theta_G theta^(-8/5),
       C_G^B <= 0.74544057 (theta = 0.5441333187), so CF.3's C_G is now explicit on the Birkhoff route (no spectral
       input; CS.3's sharp spectral C_G* = 0.02105127810).  Consistency: s_4 = c0 F0(x0) and
       s_5 = c0 (Lt F0)(x0) (c0 = 3 sigma_-3, x0 = sqrt13/2) against CC's certified level sums; DC.3's bound holds on
       CC's s_4..s_13.
  DC.4 [proved; arb] What the Birkhoff route costs: Theta_a(f_0, Lt f_0) >= log(sup R_0/inf R_0) (CF.1's last remark),
       bounded below by certified point values of R_0 (cell ends, midpoints, sqrt5/2, 5/2), so
       no a-cone bound for this f_0 is below that floor; the certified ratios Theta_G/|log(s_5/(r s_4))| and
       C_G^B/C_G* (CS.3) are printed: floor >= 0.0633864 > 4.89 |log(s_5/(r s_4))|, |log(s_5/(r s_4))| = 0.012937272,
       Theta_G/|log(s_5/(r s_4))| = 21.762, C_G^B/C_G* = 35.411.  This is a comparison of two proved upper bounds, not a failure of either; on CC's data the
       largest deviation/bound ratio is below 0.046.
Kept limits: a = 47/2 and theta are CF.2's (not re-optimised); the floor theta^floor((k-4)/5) is CF.3's block form, so
C_G^B carries theta^(-8/5); Theta_G is an upper bound (interval Horner on cells), not the value of Theta_a; the general K(b)
uses only osc log f <= b |X_s|.  Assumptions used as stated (not re-run): CE.1-CE.3 (s_k = 3 sigma_-3 (Lt^k G)(-3)),
CE.3's cone lemma, CF.1-CF.3 (CF.2's cone constants re-run unchanged), CO.1-CO.2 (re-run unchanged), CN.1's box argument
(box checks re-run), CS.3's formula for C_G* (recomputed from CC's s_4, s_5, comparison only).  Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or
30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction as Fr

from flint import arb, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_level_rate as cf            # noqa: E402
import gasket_fractional_profile_level_third as co           # noqa: E402

cn = co.cn
FAILS: list[str] = []
RESULT: dict = {}
PREC = 200
NC = int(os.environ.get("DC_CELLS", "16384"))     # cells covering X_s
ARCS = int(os.environ.get("DC_ARCS", "4096"))    # CO.1's arcs (unchanged default)
MPTS = int(os.environ.get("DC_M", "2048"))       # CO.2's trapezoid points (unchanged default)
KMAX = 13                                        # CC's level sums s_4..s_13 for the consistency check


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _up(x: arb) -> arb:
    return arb(x.upper())


def _hull(lo: arb, hi: arb) -> arb:
    """The real ball [lo.lower(), hi.upper()]."""
    l, h = arb(lo.lower()), arb(hi.upper())
    return (l + h) / 2 + arb(0, ((h - l) / 2).upper())


# ----------------------------------------------------------------- DC.1: the general-f constant
def cone_data():
    ctx.prec = PREC
    cons = cf.cone_constants()
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    M1 = arb(14) / 15
    M2 = 1 / (s0 - arb(1) / 2) - 1 / s0 - 1 / (s0 + arb(7) / 2)        # CF.2's M_2 = (log wt2)'(sqrt5/2) (abs value)
    c1 = 1 / arb(5).sqrt()
    c2 = 1 / (2 * (arb(15) / 4 + s0).sqrt())
    bstar = M1 / (1 - c1)
    X = s1 - s0
    return dict(a=cons["a"], theta=cons["theta"], M1=M1, M2=M2, c1=c1, c2=c2, bstar=bstar, X=X)


def phi(cd, b):
    return _up(arb(max((cd["M1"] + cd["c1"] * b).upper(), (cd["M2"] + cd["c2"] * b).upper())))


def Kconst(cd, beta):
    a = cd["a"]
    pb = phi(cd, beta)
    B = arb(max(arb(beta).upper(), pb.upper()))
    if not bool(B < a):
        return None
    return _up(2 * ((a + B) / (a - B)).log() + (arb(beta) + pb) * cd["X"])


def j0_of(cd, b):
    a, bs, c1 = cd["a"], cd["bstar"], cd["c1"]
    b = arb(b)
    if bool(b < a):
        return 0, b
    j = 0
    while True:
        j += 1
        beta = _up(bs + c1 ** j * (b - bs))
        if bool(beta < a):
            return j, beta


def general_constant(cd):
    ctx.prec = PREC
    a, bs, M1, M2, c1, c2 = cd["a"], cd["bstar"], cd["M1"], cd["M2"], cd["c1"], cd["c2"]
    ok_i = bool(M2 + c2 * bs < bs) and bool(M2 < M1) and bool(c2 < c1) and bool(bs < a)
    check(ok_i, f"DC.1(i) (arb): b* = M_1/(1 - c_1) = {bs.str(10)}; M_2 + c_2 b* = {(M2 + c2 * bs).str(8)} < b*, M_2 = "
                f"{M2.str(8)} < M_1 = 14/15, c_2 = {c2.str(8)} < c_1 = {c1.str(8)}: phi(b) <= b* for b <= b*, "
                f"phi(b) - b* <= c_1 (b - b*) for b >= b*, so Lt^j C_b c C_(b* + c_1^j (b - b*)_+)")
    # K increasing in beta (each term is), sampled table
    tab = {}
    for name, b in (("0", arb(0)), ("M1", M1), ("b*", bs), ("10", arb(10)), ("20", arb(20)), ("23", arb(23)),
                    ("100", arb(100)), ("1e4", arb(10 ** 4))):
        j, beta = j0_of(cd, b)
        K = Kconst(cd, beta)
        tab[name] = dict(b=b, j0=j, beta=beta, K=K)
    ok_tab = all(v["K"] is not None for v in tab.values())
    # CF.4's Theta_0 for f = 1 (B = 14/15, exact osc log Lt1 = log(11/6)) must be <= K(0)
    Theta0_cf4 = 2 * ((a + M1) / (a - M1)).log() + (arb(11) / 6).log()
    ok_cf4 = bool(Theta0_cf4 <= tab["0"]["K"])
    check(ok_tab and ok_cf4,
          "DC.1(ii) (arb): K(beta) = 2 log((a+B)/(a-B)) + (beta + phi(beta))|X_s|, |X_s| = " + cd["X"].str(10) + "; "
          + "; ".join(f"b = {k}: j0 = {v['j0']}, beta = {v['beta'].str(6)}, K <= {v['K'].str(6)}" for k, v in tab.items())
          + f"; CF.4's Theta_0 = {Theta0_cf4.str(8)} <= K(0) (consistency)")
    return tab


# ----------------------------------------------------------------- DC.2: the level-sum function F0 on X_s
def taylor_data():
    cn.tail_constants()
    cn.box_checks()
    co.big_disk()
    S_rho = co.circle_sup_rho(K=ARCS)
    f = co.taylor_coeffs(S_rho, M=MPTS)
    return S_rho, f


class Taylor:
    """F0(s) = sum_n f_n u^n, u = (s - c)/R, with f_0..f_N enclosed and |f_n| <= S_rho rho^-n for n > N (Cauchy)."""

    def __init__(self, f, S_rho, rho=co.RHO_BIG):
        ctx.prec = PREC
        self.f = f
        self.N = len(f) - 1
        self.df = [f[n] * n for n in range(1, self.N + 1)]
        self.S = S_rho
        self.rho = arb(rho.numerator) / rho.denominator
        self.c = arb(cn.C_CEN.numerator) / cn.C_CEN.denominator
        self.R = arb(cn.R_RAD.numerator) / cn.R_RAD.denominator

    def _tails(self, U):
        q = arb(abs(U).upper()) / self.rho        # |u|/rho, must be < 1
        assert bool(q < 1)
        N1 = self.N + 1
        t0 = self.S * q ** N1 / (1 - q)                                       # sum_(n>N) S rho^-n |u|^n
        # sum_(n>N) n S rho^-n |u|^(n-1) = (S/rho) sum_(n>N) n q^(n-1) = (S/rho) q^N (N1 - N q)/(1 - q)^2
        t1 = self.S / self.rho * q ** self.N * (N1 - self.N * q) / (1 - q) ** 2
        return t0, t1

    @staticmethod
    def _horner(coef, U):
        acc = arb(0)
        for t in reversed(coef):
            acc = acc * U + t
        return acc

    def val(self, s):
        U = (s - self.c) / self.R
        t0, _ = self._tails(U)
        return self._horner(self.f, U) + arb(0, t0.upper())

    def der(self, s):
        U = (s - self.c) / self.R
        _, t1 = self._tails(U)
        return (self._horner(self.df, U) + arb(0, t1.upper())) / self.R


def wts(s):
    h, s7 = arb(1) / 2, arb(7) / 2
    return (s + h) / (2 * s * (s7 - s)), (s - h) / (2 * s * (s + s7))


def LtF(T, s):
    q4 = arb(15) / 4
    w1, w2 = wts(s)
    return w1 * T.val((q4 - s).sqrt()) + w2 * T.val((q4 + s).sqrt())


def cell_scan(T, cd):
    ctx.prec = PREC
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    s0l, s1u = arb(s0.lower()), arb(s1.upper())
    h = (s1u - s0l) / NC
    Fmin = Fmax = None
    Lmin = Lmax = None
    Rmin = Rmax = None
    b0 = arb(0)
    Pmax = Pmin = None          # point values of R_0: certified lower bound for log(sup/inf) (DC.4)
    for i in range(NC):
        J = _hull(s0l + h * i, s0l + h * (i + 1))
        F = T.val(J)
        dF = T.der(J)
        LF = LtF(T, J)
        if not (bool(F > 0) and bool(LF > 0)):
            check(False, f"DC.2: positivity not certified on cell {i}")
            return None
        R0 = LF / F
        for xm in (arb(J.mid()), arb(J.lower()), arb(J.upper())):
            if not (bool(xm >= s0) and bool(xm <= s1)):
                continue
            Rp = LtF(T, xm) / T.val(xm)
            Pmax = Rp.lower() if Pmax is None else max(Pmax, Rp.lower())
            Pmin = Rp.upper() if Pmin is None else min(Pmin, Rp.upper())
        b0 = arb(max(b0.upper(), (abs(dF) / F).upper()))
        Fmin = F.lower() if Fmin is None else min(Fmin, F.lower())
        Fmax = F.upper() if Fmax is None else max(Fmax, F.upper())
        Lmin = LF.lower() if Lmin is None else min(Lmin, LF.lower())
        Lmax = LF.upper() if Lmax is None else max(Lmax, LF.upper())
        Rmin = R0.lower() if Rmin is None else min(Rmin, R0.lower())
        Rmax = R0.upper() if Rmax is None else max(Rmax, R0.upper())
    for xm in (s0, s1):
        Rp = LtF(T, xm) / T.val(xm)
        Pmax, Pmin = max(Pmax, Rp.lower()), min(Pmin, Rp.upper())
    Fmin, Fmax, Lmin, Lmax, Rmin, Rmax = (arb(x) for x in (Fmin, Fmax, Lmin, Lmax, Rmin, Rmax))
    logR_lo = arb((arb(Pmax) / arb(Pmin)).log().lower())
    oscF = _up((Fmax / Fmin).log())
    oscL = _up((Lmax / Lmin).log())
    logR = _up((Rmax / Rmin).log())
    a = cd["a"]
    pb = phi(cd, b0)
    B0 = arb(max(b0.upper(), pb.upper()))
    ok = bool(Fmin > 0) and bool(B0 < a)
    ThetaG = _up(2 * ((a + B0) / (a - B0)).log() + logR)
    ThetaCF = _up(2 * ((a + B0) / (a - B0)).log() + oscF + oscL)
    check(ok, f"DC.2 (arb, {NC} cells, Taylor f_0..f_{T.N} + Cauchy tail): F0 in [{Fmin.str(8)}, {Fmax.str(8)}] > 0 on X_s, "
              f"b_0 = sup|F0'/F0| <= {b0.str(8)} < a = {a.str(4)}, phi(b_0) <= {pb.str(8)}, B_0 = {B0.str(8)}; "
              f"osc log F0 <= {oscF.str(8)}, osc log Lt F0 <= {oscL.str(8)}; R_0 = Lt F0/F0 in [{Rmin.str(10)}, {Rmax.str(10)}], "
              f"{logR_lo.str(8)} <= log(sup R_0/inf R_0) <= {logR.str(8)} (lower end from point values at cell ends, midpoints and sqrt5/2, 5/2)")
    check(bool(ThetaG <= ThetaCF),
          f"DC.2: Theta_a(f_0, Lt f_0) <= Theta_G = 2 log((a+B_0)/(a-B_0)) + log(sup R_0/inf R_0) = {ThetaG.str(8)} "
          f"(CF.3's form with osc + osc: {ThetaCF.str(8)})")
    return dict(Fmin=Fmin, Fmax=Fmax, b0=b0, B0=B0, oscF=oscF, oscL=oscL, Rmin=Rmin, Rmax=Rmax, logR=logR, logR_lo=logR_lo,
                ThetaG=ThetaG, ThetaCF=ThetaCF)


# ----------------------------------------------------------------- DC.3 / DC.4
def level_sum_bound(T, cd, sc):
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    theta = cd["theta"]
    ThetaG = sc["ThetaG"]
    eig = cn.spectral_layout()
    r = eig["r"]
    ok_r = bool(r >= sc["Rmin"]) and bool(r <= sc["Rmax"])
    sums = cn.level_data(KMAX)
    c0 = 3 * ce.sigma_closed(arb(-3), 0)
    x0 = arb(13).sqrt() / 2
    s4p, s5p = c0 * T.val(x0), c0 * LtF(T, x0)
    ok_s = bool(abs(s4p - sums[4]) < arb("1e-12") * sums[4]) and bool(abs(s5p - sums[5]) < arb("1e-12") * sums[5])
    check(ok_s and ok_r, f"DC.3 (consistency): c0 F0(x0) = {s4p.str(12)} vs CC's s_4 = {sums[4].str(12)}; c0 (Lt F0)(x0) = "
                         f"{s5p.str(12)} vs s_5 = {sums[5].str(12)}; r = {r.str(12)} in [inf R_0, sup R_0]")
    CB = _up(ThetaG / (theta.log() * 8 / 5).exp())
    worst = arb(0)
    okd = True
    for k in range(4, KMAX):
        dev = abs((sums[k + 1] / (r * sums[k])).log())
        bnd = theta ** ((k - 4) // 5) * ThetaG
        okd = okd and bool(dev < bnd) and bool(bnd <= CB * (theta.log() * k / 5).exp())
        worst = arb(max(worst.upper(), (dev / bnd).upper()))
    check(okd, f"DC.3: for every k >= 4, |log(s_(k+1)/(r s_k))| <= theta^floor((k-4)/5) Theta_G <= C_G^B theta^(k/5), "
               f"C_G^B = Theta_G theta^(-8/5) = {CB.str(8)} (theta = {theta.str(10)}); on CC's s_4..s_{KMAX} every deviation "
               f"is below the bound (max ratio deviation/bound <= {worst.str(4)})")
    # DC.4 comparison with CS.3's sharp spectral constant
    log4 = abs((sums[5] / (r * sums[4])).log())
    CGs = log4 / (theta.log() * 4 / 5).exp()
    floor = sc["logR_lo"]          # a certified LOWER bound of log(sup R_0/inf R_0) (point values)
    check(bool(floor > log4 * arb("4.89")) and bool(ThetaG > floor),
          f"DC.4 (arb): Birkhoff floor log(sup R_0/inf R_0) = {floor.str(8)} > 4.89 |log(s_5/(r s_4))| = 4.89 * {log4.str(8)}; "
          f"Theta_G/|log(s_5/(r s_4))| = {(ThetaG / log4).str(6)}, C_G^B/C_G* = {(CB / CGs).str(6)} (C_G* = {CGs.str(10)}, CS.3)")
    return dict(CB=CB, CGs=CGs, log4=log4, worst=worst, s4p=s4p, s5p=s5p, r=r, ratioT=ThetaG / log4, ratioC=CB / CGs)


def main() -> int:
    t0 = time.time()
    print("Theorem DC: CF.3's Birkhoff-route constant made explicit (general f and the level sums). Not thrust.", flush=True)
    cd = cone_data()
    tab = general_constant(cd)
    RESULT.update(cd=cd, tab=tab)
    S_rho, f = taylor_data()
    if f is None:
        print("1 FAILED")
        return 1
    T = Taylor(f, S_rho)
    sc = cell_scan(T, cd)
    if sc is None:
        print("FAILED")
        return 1
    RESULT.update(sc)
    RESULT.update(level_sum_bound(T, cd, sc))
    RESULT["S_rho"] = S_rho
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(co.FAILS) + len(cn.FAILS) + len(cf.FAILS) + len(co.cg.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DC checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
