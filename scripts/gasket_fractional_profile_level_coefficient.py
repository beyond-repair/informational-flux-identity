#!/usr/bin/env python3
r"""Theorem CN: B != 0 in CG.4.  The level sums of nu have the exact rate |l_2|/r, with alternating signs for every level,
an explicit constant, and a certified coefficient B.

Setting (GASKET.md, Theorems CE-CG): s_k = 3 sigma_-3 (Lt^k G)(-3) (CE.1-CE.2), in s = sqrt(t + 25/4)
    (Lt f)(s) = wt1(s) f(sigma_1 s) + wt2(s) f(sigma_2 s),   sigma_(1,2) = sqrt(15/4 -+ s),   X = [sqrt5/2, 5/2],
G = 1_I g/q, g(t) = lam^-gamma/(Pi(u) psi'(u)), u = psi^-1(t) in (-5,-1), lam = -u, q = -t, gamma = log3/log5.
F0 = Lt^4 G is the sum over CE.3's five inner level-4 words Phi of wt_Phi * (g/q)(t(Phi s)).  H^2(D), D = {|s - c| < R},
c = 83/50, R = 8/5, as in CG; x0 = s(-3) = sqrt13/2.  CG.4 left open whether B != 0 in
    s_(k+1)/s_k = r (1 + B (l_2/r)^k (l_2/r - 1) + O((l_2/r)^(2k) + (l_3/r)^k)),
i.e. whether |l_2|/r is the exact rate of the level sums.  Here B is defined by s_n = kappa r^n (1 + B (l_2/r)^n + ...).

Theorem CN (checked here; exit 1 on any failure):
  CN.1 [certified; arb/acb] On the convex box V = [-26/5, -4/5] x [-3/10, 3/10] (u-plane): Re psi' >= 0.5436 > 0 (so psi
       is injective on V), |psi| >= 0.7647, Pi is finite and nonzero, Re(-u) >= 4/5 (1056 boxes).  For every arc of a
       1024-arc cover of |s - c| = R and each of the five inner words, a complex Krawczyk test puts psi^-1(t(Phi s)) in a
       box U inside V.  Hence F0 continues holomorphically to a neighbourhood of the closed disk (psi^-1 is a holomorphic
       bijection onto V° of the simply connected psi(V°)), and ||F0||_(H^2(D)) <= sup_(|s-c|=R) |F0| <= S = 0.0169081...
       Consistency: 3 sigma_-3 F0(x0) = CC's s_4 to 1e-18.
  CN.2 [proved; resolvent circles certified] sup ||(z - Lt)^-1|| <= 1.0142e9 on |z| = 1e-3, 927.46 on |z - r| = 1/10,
       6.3675e4 on |z - l_2| = 1/50, 2.9324e6 on |z - l_3| = 1/125.  So for every k >= 0
           s_(k+4) = a_1 r^k + a_2 l_2^k + a_3 l_3^k + e_k,   |e_k| <= eps_0 1e-3^k,  eps_0 <= 2.889e4,
       a_j = 3 sigma_-3 (Pi_j F0)(x0), |a_1| <= 2.642, |a_2| <= 36.28, |a_3| <= 668.2 (Riesz projections on H^2(D)).
  CN.3 [certified] Annihilating two of the three exponentials (rational centres) in CC's certified s_4..s_13:
           a_1 = 0.00571256871089 +- 1.2e-15,  a_2 = -4.158e-5 +- 4.4e-9 != 0,  |a_3| <= 1.125e-4,
           B = (a_2/a_1)(r/l_2)^4 = -8.211 +- 4.6e-4 != 0,  kappa = a_1/r^4 = 0.1772022812...
  CN.4 [certified] For every n >= 4:  |log(s_(n+1)/(r s_n))| <= C (|l_2|/r)^n with C = 14.5952; the deviation
       s_(n+1)/s_n - r has the sign of (-1)^n (alternates at every level); for n >= 13 it is >= 9.626 (|l_2|/r)^n and
       log(s_(n+1)/(r s_n)) / (l_2/r)^n -> B (l_2/r - 1) = 9.628...  So |l_2|/r is the exact rate (CG.4's open item).
Kept limits: whether a_3 != 0 is not decided; CF.3's Birkhoff constant C_G is not computed (not needed for the level
sums any more); no closed form of B, kappa, r or l_2.
Assumptions used as stated (not re-run): CE.1-CE.3, CG.1-CG.3 (enclosures of l_2, l_3, the eigenvalue count: exactly
r, l_2, l_3 outside |z| <= 1e-3, all simple), CE.4's enclosure of r, CC.4's certified level sums (re-run unchanged),
CB.5's bound |psi'(v) - 1| <= 2(e^(|v|/5) - 1) and BJ.6's |psi(v) - v| <= |v|(e^(|v|/10) - 1) (nonnegative coefficients,
valid for complex v).  Textbook facts: Noshiro-Warschawski (Re f' > 0 on a convex set => injective), Krawczyk/Brouwer,
maximum principle, Riesz projections.  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW enters.
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

import mpmath as mp
from flint import acb, acb_mat, arb, ctx

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_profile_atoms as bz            # noqa: E402
import gasket_fractional_profile_product as cb          # noqa: E402

FAILS: list[str] = []
J = cb.J                       # 200 halvings 5^-J, as in BZ/CB
PREC = 200
C_CEN, R_RAD = Fr(83, 50), Fr(8, 5)
V_BOX = (Fr(-26, 5), Fr(-4, 5), Fr(3, 10))      # Re u in [-5.2, -0.8], |Im u| <= 0.3
INNER = [(0, 0, 0, 1), (0, 1, 0, 0), (0, 1, 0, 1), (1, 1, 0, 0), (1, 1, 0, 1)]   # innermost map first
# CG.3 enclosures (as stated) and CE.4's r (as stated)
R_MID = "0.42373128696260793985253216132843923"
L2_MID, L2_RAD = "-0.0731138175123995805", "1e-13"
L3_MID, L3_RAD = "-0.0170676440626237159", "1e-12"
# resolvent circles (certified here): (centre, radius)
RHO_IN = Fr(1, 1000)
CIRCLES = {"r": ("0.4237312869626079", Fr(1, 10)), "l2": ("-0.07311381751239958", Fr(1, 50)),
           "l3": ("-0.017067644062623716", Fr(1, 125))}
N_TRUNC = 80
K_DATA = 13


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def _a(x: Fr) -> arb:
    return arb(x.numerator) / x.denominator


def _box(re_lo, re_hi, im_lo, im_hi) -> acb:
    return acb(arb((re_lo + re_hi) / 2, ((re_hi - re_lo) / 2).upper()), arb((im_lo + im_hi) / 2, ((im_hi - im_lo) / 2).upper()))


# ----------------------------------------------------------------- psi, psi', Pi on complex balls
def _orbit(u):
    """x[j] = psi(5^-j u) (BZ/CB construction, complex balls): x[J] = v + B(|v|(e^(|v|/10) - 1)), x[j-1] = p(x[j])."""
    x, eps, ball, av = cb.orbit(u, cplx=True)
    return x, av


def psi_c(u):
    return _orbit(u)[0][0]


def dpsi_c(u, orb=None):
    x, av = orb if orb is not None else _orbit(u)
    dl = arb((2 * ((av / 5).exp() - 1)).upper())
    d = acb(1) + acb(arb(0, dl), arb(0, dl))
    for i in range(1, J + 1):
        d *= 1 + 2 * x[i] / 5
    return d


def Pi_c(u, orb=None):
    """Pi(u) = prod_(j>=1) (5/3) c(psi(5^-j u)), c(x) = (x+6)/((x+2)(x+5)); the tail j > J has |log| <= |u| 5^-J
    (|log((5/3)c(x))| <= 1.072|x| for |x| <= 1/2 and sum_(j>J) |psi(5^-j u)| <= 1.01 |u| 5^-J / 4)."""
    x, av = orb if orb is not None else _orbit(u)
    f3 = arb(5) / 3
    P = acb(1)
    for j in range(1, J + 1):
        bz._nz(x[j] + 2)
        bz._nz(x[j] + 5)
        P *= f3 * bz._c(x[j])
    eps = arb((abs(u).upper() * arb(5) ** (-J)).upper())
    return P * acb(arb(0, eps), arb(0, eps)).exp()


def g_c(u, gam):
    orb = _orbit(u)
    return (-u) ** (-gam) / (Pi_c(u, orb) * dpsi_c(u, orb))


# ----------------------------------------------------------------- CN.1 (a): the box V
def tail_constants():
    """Exact: (5/3)c(x) = (1+x/6)/((1+x/2)(1+x/5)), so for |x| <= 1/2, |log((5/3)c(x))| <= sum |z|/(1-|z|) <= 1.072|x|;
    with |psi(v)| <= |v| e^(|v|/10) <= 1.01|v| for |v| <= 1/10, sum_(j>J) |psi(5^-j u)| <= 1.01 |u| 5^-J / 4, and
    1.072 * 1.01 / 4 < 1: the Pi tail has |log| <= |u| 5^-J (used for complex u)."""
    import sympy as sp
    x = sp.symbols("x")
    ident = sp.simplify(sp.Rational(5, 3) * (x + 6) / ((x + 2) * (x + 5)) - (1 + x / 6) / ((1 + x / 2) * (1 + x / 5))) == 0
    k = Fr(1, 6) / (1 - Fr(1, 12)) + Fr(1, 2) / (1 - Fr(1, 4)) + Fr(1, 5) / (1 - Fr(1, 10))
    check(ident and k <= Fr(1072, 1000) and Fr(1072, 1000) * Fr(101, 100) / 4 < 1,
          f"CN.1 (sympy, exact): (5/3)c(x) = (1+x/6)/((1+x/2)(1+x/5)); |log((5/3)c(x))| <= {float(k):.5f}|x| <= 1.072|x| "
          "(|x| <= 1/2); Pi tail |log| <= 1.072 * 1.01/4 |u| 5^-J < |u| 5^-J, also for complex u")


def box_checks(nx=88, ny=12):
    """Cover V by nx * ny closed boxes; on each: Re psi' > 0, psi != 0, Pi finite and != 0, Re(-u) > 0."""
    ctx.prec = PREC
    lo, hi, h = (_a(v) for v in V_BOX)
    dx, dy = (hi - lo) / nx, 2 * h / ny
    min_re, min_psi = None, None
    ok = True
    for i in range(nx):
        for j in range(ny):
            b = _box(lo + i * dx, lo + (i + 1) * dx, -h + j * dy, -h + (j + 1) * dy)
            try:
                orb = _orbit(b)
                d = dpsi_c(b, orb)
                P = Pi_c(b, orb)
                ps = orb[0][0]
            except bz.Fail:
                ok = False
                continue
            re_d = d.real.lower()
            ap = abs(ps).lower()
            min_re = re_d if min_re is None else min(min_re, re_d)
            min_psi = ap if min_psi is None else min(min_psi, ap)
            ok = ok and bool(d.real > 0) and bool(abs(ps) > 0) and not P.contains(0) and P.is_finite() and bool((-b).real > 0)
    check(ok and bool(min_re > 0), f"CN.1 (acb, {nx * ny} boxes): on V = [-26/5,-4/5] x [-3/10,3/10]: Re psi' >= {arb(min_re).str(4)} > 0 "
          f"(psi injective on the convex V), |psi| >= {arb(min_psi).str(4)}, Pi finite and != 0, Re(-u) >= 4/5")
    return ok


# ----------------------------------------------------------------- CN.1 (b): Krawczyk + sup of F0 on the circle
def _psi_mp(u):
    x = u / mp.mpf(5) ** 60
    for _ in range(60):
        x = x * (x + 5)
    return x


def _acb_mid(z):
    return acb(arb(z.real.mid()), arb(z.imag.mid()))


def krawczyk(T, u_guess, rad):
    """Return a box U = u~ + [-rad, rad]^2 with K(U) inside int U (so psi has a unique zero of psi - t in U for every t in T),
    or None.  K(U) = u~ - Y(psi(u~) - T) + (1 - Y psi'(U))(U - u~), Y = 1/psi'(u~) (mean value form for holomorphic psi;
    the integral of psi' over a segment lies in the convex box enclosure psi'(U))."""
    ut = _acb_mid(u_guess)
    for _ in range(3):                       # polish the proposal (midpoints only)
        ut = _acb_mid(ut - (psi_c(ut) - _acb_mid(T)) / dpsi_c(ut))
    tr = arb(max(T.real.rad(), T.imag.rad()))
    R = arb(max(arb(rad), (3 * tr / abs(dpsi_c(ut)).lower())).upper())
    U = ut + acb(arb(0, R), arb(0, R))
    Y = _acb_mid(1 / dpsi_c(ut))
    K = ut - Y * (psi_c(ut) - T) + (1 - Y * dpsi_c(U)) * (U - ut)
    inside = bool(abs(K.real - ut.real) < R) and bool(abs(K.imag - ut.imag) < R)
    return (U, ut) if inside else (None, ut)


def word_eval(y, w):
    """(weight, s-image) of the composition, innermost map w[0] first."""
    q4, h, s7 = arb(15) / 4, arb(1) / 2, arb(7) / 2
    wt, s = acb(1), y
    for i in w:
        wt *= (s + h) / (2 * s * (s7 - s)) if i == 0 else (s - h) / (2 * s * (s + s7))
        s = (q4 - s).sqrt() if i == 0 else (q4 + s).sqrt()
    return wt, s


def F0_bound(y, guesses, gam, rad=arb("2e-3")):
    """Enclosure of F0 on the ball y (complex), via Krawczyk boxes; also certifies U in int V."""
    lo, hi, hh = (_a(v) for v in V_BOX)
    tot = acb(0)
    new = []
    for w, ug in zip(INNER, guesses):
        wt, s = word_eval(y, w)
        T = s * s - arb(25) / 4
        U, ut = krawczyk(T, ug, rad)
        if U is None:
            return None, guesses
        inV = bool(U.real > lo) and bool(U.real < hi) and bool(abs(U.imag) < hh)
        if not inV:
            return None, guesses
        tot += wt * g_c(U, gam) / (-T)
        new.append(ut)
    return tot, new


def circle_sup(K=256):
    ctx.prec = PREC
    mp.mp.dps = 30
    gam = arb(3).log() / arb(5).log()
    c, R = _a(C_CEN), _a(R_RAD)
    pi = arb.pi()
    # float proposals for theta = 0 (real s = c + R), then continuation along the arcs
    guesses = []
    y0 = acb(c + R)
    for w in INNER:
        _, s = word_eval(y0, w)
        t = float((s * s - arb(25) / 4).real.mid())
        u0 = mp.findroot(lambda u: _psi_mp(u) - t, (-5, -1), solver="anderson")
        guesses.append(acb(arb(mp.nstr(u0, 30))))
    S = arb(0)
    ok = True
    for j in range(K):
        th = arb(2 * j + 1) * pi / K + arb(0, (pi / K).upper())
        y = c + R * acb(0, th).exp()
        val, guesses = F0_bound(y, guesses, gam)
        if val is None:
            ok = False
            break
        S = arb(max(S, arb(abs(val).upper())).upper())
    check(ok, f"CN.1 (acb, {K} arcs x 5 words): Krawczyk boxes U (half-side >= 2e-3, adaptive) with K(U) in int U, U in int V, for every "
              "t = t(Phi s), |s - c| = R; so t(Phi(closed D)) lies in psi(V°) and F0 is holomorphic near the closed disk")
    check(ok and bool(S < arb("0.021")), f"CN.1: sup_(|s-c|=R) |F0| <= S = {S.str(6)} (so ||F0||_H2(D) <= S)")
    return S


def real_consistency(s4):
    """F0(x0) from the same continuation, times 3 sigma_-3, against CC's s_4 (consistency)."""
    ctx.prec = PREC
    gam = arb(3).log() / arb(5).log()
    x0 = acb(arb(13).sqrt() / 2)
    guesses = []
    for w in INNER:
        _, s = word_eval(x0, w)
        t = float((s * s - arb(25) / 4).real.mid())
        u0 = mp.findroot(lambda u: _psi_mp(u) - t, (-5, -1), solver="anderson")
        guesses.append(acb(arb(mp.nstr(u0, 30))))
    val, _ = F0_bound(x0, guesses, gam, rad=arb("1e-20"))
    return val


# ----------------------------------------------------------------- CN.2: resolvent circles on H^2(D)
def resolvent_constants(N=N_TRUNC):
    """Certified sup ||(z - Lt)^-1|| on the four circles, from CG's matrix machinery (re-run here):
    ||(z - Lt)^-1|| <= K/(1 - K eps_N) whenever ||(z - Lt P_N)^-1|| <= K and K eps_N < 1 (Neumann)."""
    import gasket_fractional_profile_level_spectrum as cg
    dc = cg.disk_constants()
    eps, beta = cg.tail_constants(dc, N)
    M = cg.build_matrix(N, prec=320)
    ctx.prec = cg.PREC_RES
    Mc = acb_mat(M)
    Id = acb_mat(N + 1, N + 1)
    for i in range(N + 1):
        Id[i, i] = 1
    out = {}
    circles = [("in", acb(0), RHO_IN)] + [(k, acb(arb(v[0])), v[1]) for k, v in CIRCLES.items()]
    for name, cen, rad in circles:
        ok, w, n = cg.circle_bound(Mc, Id, cen, rad, eps, beta)
        KL = arb((w / (1 - w * eps)).upper())
        check(ok, f"CN.2 (arb): circle |z - {cen.real.str(8) if name != 'in' else 0}| = {rad}: sup ||(z - Lt P_N)^-1|| <= {w.str(5)}, "
                  f"K eps_N < 1, so sup ||(z - Lt)^-1|| <= {KL.str(5)} ({n} samples)")
        out[name] = (KL, _a(rad))
    for nm, f in cg.FAILS:
        pass
    return out, cg


def spectral_layout():
    """The four circles separate CG.3's eigenvalues: the disk around each of r, l_2, l_3 holds exactly that eigenvalue
    (and misses |z| <= 1e-3 and the other two); |z| = 1/1000 encloses exactly the rest (CG.3: all other |mu| < 1e-3)."""
    ctx.prec = PREC
    r = arb(R_MID) + arb(0, arb("1e-35"))
    l2 = arb(L2_MID) + arb(0, arb(L2_RAD))
    l3 = arb(L3_MID) + arb(0, arb(L3_RAD))
    eig = {"r": r, "l2": l2, "l3": l3}
    ok = True
    for k, (c, rad) in CIRCLES.items():
        cc_, rr = arb(c), _a(rad)
        for k2, e in eig.items():
            d = abs(e - cc_)
            ok = ok and (bool(d < rr) if k2 == k else bool(d > rr))
        ok = ok and bool(abs(cc_) - rr > _a(RHO_IN))
    ok = ok and bool(abs(l3) > _a(RHO_IN))
    check(ok, "CN.2: each circle |z - m| = rad (r: 1/10, l_2: 1/50, l_3: 1/125) encloses exactly its eigenvalue and lies in "
              "|z| > 1/1000; |z| = 1/1000 (certified resolvent below, so it misses the spectrum) encloses exactly the spectrum "
              "other than r, l_2, l_3 (CG.3: all other |mu| < 1e-3)")
    return eig


# ----------------------------------------------------------------- CN.2-CN.4
def expansion_constants(S, res):
    ctx.prec = PREC
    import gasket_fractional_profile_level_ratio as ce
    sig3 = ce.sigma_closed(arb(-3), 0)
    c0 = 3 * sig3
    c, R = _a(C_CEN), _a(R_RAD)
    u0 = (arb(13).sqrt() / 2 - c) / R
    ev = 1 / (1 - u0 * u0).sqrt()
    alpha = {k: c0 * ev * res[k][1] * res[k][0] * S for k in ("r", "l2", "l3")}
    KLin, rho = res["in"]
    eps0 = c0 * ev * rho * KLin * S
    check(bool(abs(u0) < 1), f"CN.2: x0 = sqrt13/2 in D, |u0| = {abs(u0).str(6)}; |h(x0)| <= ||h|| / sqrt(1 - u0^2) = {ev.str(6)} ||h|| on H^2(D)")
    print(f"   c0 = 3 sigma_-3 = {c0.str(12)};  alpha_r <= {alpha['r'].str(4)}, alpha_l2 <= {alpha['l2'].str(4)}, "
          f"alpha_l3 <= {alpha['l3'].str(4)}, eps_0 <= {eps0.str(4)}, rho = {rho.str(4)}")
    check(all(bool(v.is_finite()) for v in alpha.values()) and bool(eps0.is_finite()),
          "CN.2: s_(k+4) = a_1 r^k + a_2 l_2^k + a_3 l_3^k + e_k for all k >= 0, |a_j| <= alpha_j, |e_k| <= eps_0 (1/1000)^k "
          "(Riesz projections of F0 in H^2(D); Pi_j rank one with Lt Pi_j = mu_j Pi_j; Lt^k Pi_in = (1/2 pi i) oint z^k (z - Lt)^-1 dz)")
    return dict(c0=c0, ev=ev, alpha=alpha, eps0=eps0, rho=rho)


def level_data(K=K_DATA):
    import gasket_fractional_profile_pure_point as cc
    found, sums, counts, okb = cc.period_atoms(K, verbose=False)
    check(okb and len(sums) == K + 1, f"CN.3: CC.4's certified level sums s_0..s_{K} re-run ({len(found)} atoms, CB.5 masses, arb)")
    return sums


def _annihilate(sums, i, j, k, ec, eig):
    """Enclose a_i from s_(k+4), s_(k+5), s_(k+6), annihilating eigenvalues j1, j2 approximately (rational centres)."""
    names = ["r", "l2", "l3"]
    mids = {"r": arb(R_MID), "l2": arb(L2_MID), "l3": arb(L3_MID)}
    others = [n for n in names if n != i]
    A, Bm = mids[others[0]], mids[others[1]]
    Q = sums[k + 6] - (A + Bm) * sums[k + 5] + A * Bm * sums[k + 4]
    err = arb(0)
    for n in others:
        mu = eig[n]
        err += ec["alpha"][n] * abs(mu) ** k * abs(mu - A) * abs(mu - Bm)
    rho = ec["rho"]
    err += ec["eps0"] * rho ** k * (rho * rho + abs(A + Bm) * rho + abs(A * Bm))
    mu = eig[i]
    den = mu ** k * (mu - A) * (mu - Bm)
    return (Q + arb(0, err.upper())) / den, Q, err


def coefficients(sums, ec, eig, kmax=None):
    ctx.prec = PREC
    kmax = len(sums) - 7 if kmax is None else kmax
    best = {}
    for name in ("r", "l2", "l3"):
        enc = None
        for k in range(0, kmax + 1):
            e, Q, err = _annihilate(sums, name, None, k, ec, eig)
            if enc is None or e.rad() < enc.rad():
                enc, kb, Qb, eb = e, k, Q, err
        best[name] = enc
        print(f"   a_{name}: best k = {kb}: |Q_k| = {abs(Qb).str(4)} vs bound on the other terms {eb.str(4)}; enclosure {enc.str(8)}")
    a1, a2, a3 = best["r"], best["l2"], best["l3"]
    ec["a3"] = arb(abs(a3).upper())
    print(f"   a_3 enclosure (data): {a3.str(6)}" + ("  (excludes 0)" if not a3.contains(0) else "  (contains 0)"))
    check(bool(ec["a3"] < ec["alpha"]["l3"]), f"CN.3: |a_3| <= {ec['a3'].str(4)} (data; the projection bound gave {ec['alpha']['l3'].str(4)})")
    check(bool(a1 > 0), f"CN.3: a_1 = {a1.str(12)} > 0")
    check(not a2.contains(0), f"CN.3: a_2 = {a2.str(10)} != 0  (the l_2-component of the level sums is not zero)")
    r, l2 = eig["r"], eig["l2"]
    B = a2 / a1 * (r / l2) ** 4
    kappa = a1 / r ** 4
    check(not B.contains(0), f"CN.3: B = (a_2/a_1)(r/l_2)^4 = {B.str(8)} != 0 (CG.4's B); kappa = a_1/r^4 = {kappa.str(10)}")
    return a1, a2, B


def exact_rate(sums, ec, eig, a1, a2, B, n1=None):
    """CN.4: explicit C for all n >= 4, alternating signs for all n >= 4, and the limit of the normalised deviation."""
    ctx.prec = PREC
    K = len(sums) - 1
    n1 = K if n1 is None else n1
    r, l2, l3 = eig["r"], eig["l2"], eig["l3"]
    q = abs(l2) / r
    rho = ec["rho"]
    sgnB = 1 if bool(B > 0) else -1
    # data part, 4 <= n < n1
    ok, cmax = True, arb(0)
    for n in range(4, n1):
        lr = (sums[n + 1] / (r * sums[n])).log()
        expect = -sgnB * (-1) ** n                     # sign of B (l2/r)^n (l2/r - 1)
        ok = ok and (bool(lr > 0) if expect > 0 else bool(lr < 0))
        cmax = arb(max(cmax, (abs(lr) / q ** n).upper()).upper())
    check(ok, f"CN.4 (data, arb): for 4 <= n < {n1} the deviation s_(n+1)/s_n - r has sign -sign(B)(-1)^n "
              f"(alternating), max |log(s_(n+1)/(r s_n))| / (|l_2|/r)^n <= {cmax.str(6)}")
    # analytic part, n >= n1
    a1l = arb(abs(a1).lower())
    b2 = abs(a2 / a1)

    def eta(n):
        return ec["a3"] / a1l * (abs(l3) / r) ** (n - 4) + ec["eps0"] / a1l * (rho / r) ** (n - 4)

    def bbar(n):
        return b2 * q ** (n - 4) + eta(n)

    def quad(n):
        x = bbar(n)
        return x * x / (2 * (1 - x))
    rem = eta(n1) + eta(n1 + 1) + quad(n1) + quad(n1 + 1)
    rem_n = rem / q ** n1
    main_ = abs(B) * (1 + q)
    ok2 = bool(bbar(n1) < arb("0.5")) and bool(rem_n < main_)
    Ctail = main_ + rem_n
    check(ok2, f"CN.4: for n >= {n1}: log(s_(n+1)/(r s_n)) = B (l_2/r)^n (l_2/r - 1) + rho_n with |rho_n| / (|l_2|/r)^n <= "
               f"{rem_n.str(4)} (decreasing in n) < |B|(1 + |l_2|/r) = {main_.str(8)}: the sign alternates for every n and "
               "log(s_(n+1)/(r s_n)) / (l_2/r)^n -> B (l_2/r - 1)")
    C = arb(max(cmax, Ctail).upper())
    cl = main_ - rem_n
    check(ok and ok2, f"CN.4: |log(s_(n+1)/(r s_n))| <= C (|l_2|/r)^n for every n >= 4 with C = {C.str(6)} (explicit), and "
                      f">= {cl.str(6)} (|l_2|/r)^n for every n >= {n1}: |l_2|/r = {q.str(12)} is the exact rate")
    return C, rem_n


def main(K=K_DATA, arcs=1024) -> int:  # noqa: C901
    t0 = time.time()
    print("Theorem CN: B != 0 in CG.4; the level sums have the exact rate |l_2|/r with alternating signs and an explicit "
          "constant. Not thrust.")
    tail_constants()
    box_checks()
    S = circle_sup(arcs)
    eig = spectral_layout()
    res, cg = resolvent_constants()
    ec = expansion_constants(S, res)
    sums = level_data(K)
    v = real_consistency(None)
    check(bool(abs(ec["c0"] * v.real - sums[4]) < arb("1e-20")),
          f"CN.1 (consistency): 3 sigma_-3 F0(x0) = {(ec['c0'] * v.real).str(15)} = CC's s_4 = {sums[4].str(15)}")
    a1, a2, B = coefficients(sums, ec, eig)
    exact_rate(sums, ec, eig, a1, a2, B)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS or cg.FAILS:
        print(f"{len(FAILS) + len(cg.FAILS)} FAILED")
        return 1
    print("all Theorem CN checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else K_DATA))
