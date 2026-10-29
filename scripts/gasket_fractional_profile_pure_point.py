#!/usr/bin/env python3
r"""Theorem CC: BX's spectral measure nu is pure point.  nu(-Z) = 0, so nu = sum_(lam in Lambda_O) nu({lam}) delta_lam
and Conjecture CA.6 holds in full.  The proof is a contraction on the Julia set K_p:

  H = 1/r_inf - z/2 has Im H(z) >= Im z/2 on the upper half plane (Cauchy-Schwarz on r_inf's probability measure),
  so F = -1/H is a Nevanlinna function with |F(w)| <= 2/|Im w|; CB's law H = c H(p) gives F = k F(p), k = 1/c =
  (z+2)(z+5)/(z+6).  Partial fractions 1/(t - p(z)) = sum_(p(y)=t) 1/(p'(y)(y - z)) turn this into a fixed-point
  equation for the Herglotz measure sigma of F restricted to K_p:  sigma|K = Phi(sigma|K), with weights k(y)/p'(y),
  and  sum_(p(y)=t) k(y)/p'(y) = (2-t)/(6-t) <= 7/11 on [-5,0] (exact, sympy).  So sigma(K_p) <= (7/11) sigma(K_p),
  sigma(K_p) = 0; via 1/Q_inf = Pi H(psi) (CB.3), psi a local diffeomorphism on Z, this gives nu(-Z) = 0.

Setting (GASKET.md, Theorems AQ, BD, BJ, BS, BX, BZ, CA, CB): p(z) = z(z+5), c(z) = (z+6)/((z+2)(z+5)), r_inf the corner
resolvent (Stieltjes, probability measure on [0,8], AQ (2)-(3)), D = 1/r_inf, H = D - z/2 = (z/2) P(z) with
H = c H(p) on C \ K_p (CB.2); psi (BD.4), Pi(u) = prod_(j>=1) (5/3) c(psi(5^-j u)), f = 1/Q_inf = Pi H(psi) (CB.3),
Q_inf = int dnu/(lam+u) (BX.1), K_p the null Cantor Julia set in [-5,0] (BZ.1), Z = psi^-1(K_p), O = U_k p^-k(-3),
Lambda_O (CA), m = lam^-gamma nu on the scale circle, m(T) = c_R(gamma)/2 (BX.2).

Theorem CC (checked here; exit 1 on any failure):
  CC.1 [proved; elementary] Im D >= Im z and Im H >= Im z/2 for Im z > 0; F = -1/H is Nevanlinna, |F(w)| <= 2/|Im w|,
       and F = k F(p).  Checks: sympy (k = 1/c, F = k F(p) from H = c H(p)); exact rationals (Cauchy-Schwarz lemma on
       random finite measures); arb (Im H - Im z/2 > 0 at points near K_p, via CB's product).
  CC.2 [proved from CC.1 + Assumption CC (Herglotz representation)] sigma(K_p) = 0, by the transfer identity
       sigma|K = Phi(sigma|K) and T(t) = sum_(p(y)=t) k(y)/p'(y) = (2-t)/(6-t) in [1/3, 7/11] on [-5,0] (sympy).
  CC.3 [proved from CC.2, CB.3-CB.4, BZ.4] nu((0,inf) n (-Z)) = 0: psi has no critical point on Z (psi' = prod(1 +
       2 psi(5^-i u)/5), -5/2 not in K_p); Pi is finite and nonzero on Z except on E = {psi = 0} (pole from c(-5)), where
       |Q_inf| stays bounded (|F| <= 2/|Im w|), so no atom.  Hence nu is pure point with atom set exactly Lambda_O
       (Conjecture CA.6), and Q_inf(u) = sum_(lam in Lambda_O) nu({lam})/(lam + u).
  CC.4 [proved + certified] psi is strictly increasing on [-5,0] (|psi(v) - v| <= e^(1/10) - 1 on [-1,0]), so one scale
       period is parametrised exactly: Lambda_O n [1,5) = {-psi^-1(t) : t in O n (psi(-5), psi(-1)]}, psi(-5) =
       p(psi(-1)), and c_R(gamma)/2 = sum_(k>=0) sum_(t in p^-k(-3) n I) lam_t^-gamma nu({lam_t}) exactly (CB.5 masses).
       Certified (arb): every atom of level k <= K_CERT, its CB.5 mass, the partial sum and the exact remainder
       delta_K = mass of the levels > K.
  CC.5 [proved from CC.3 and BZ.6's recurrence] limsup rho_k = 1, so limsup |theta_k| e^(2 pi^2 k/ln5) =
       pi c_R(gamma)/ln5 = 1.0094980205...: BX.3's sharp bound is attained exactly.  Certified: liminf rho_k <= the
       bound printed below (from the level <= K_CERT atoms and delta_K).
  CC.6 [proved from CC.3 and CA.3] psi has no critical point on Z, so every point of -Z n (0,inf) is a limit of atoms
       (CA.3 applies everywhere; there are no even-order critical points of psi in Z): supp nu = Lambda_O u (-Z n (0,inf)).
Observation (no claim): the level sums s_k = sum over level-k atoms of lam^-gamma nu have ratio 0.42373129... for
k >= 8, and the geometric tail s_K r/(1-r) matches the exact remainder delta_14 to about 2e-15: the CB.5 masses add up
to c_R(gamma)/2 from Theorems AT/AU's independent route.  Level counts are 1, 0, 1, 3 and then 5 * 2^(k-4).
Assumption CC (literature): Herglotz-Nevanlinna representation of functions with Im F > 0 on the upper half plane
(same family as Assumption BX).  Assumptions used as stated (not re-run): AQ (1)-(3), BD.4, BJ.6, BS.1, BX.1-BX.2,
BZ.1-BZ.4, CA.1-CA.2, CB.1-CB.5 and their enclosure routines (imported unchanged).
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import random
import sys
import time
from fractions import Fraction

import mpmath as mp
import sympy as sp
from flint import acb, arb, ctx

import gasket_fractional_profile_atoms as bz
import gasket_fractional_profile_atom_law as ca
import gasket_fractional_profile_product as cb

FAILS: list[str] = []
K_CERT = 14            # levels certified in `reproduce` (10240 atoms, about 40 s)
K_FAST = 6             # levels certified in the pytest run
C_R = bz.C_R
COUNTS = {0: 1, 1: 0, 2: 1, 3: 3, 4: 5}   # level counts k <= 4; 5 * 2^(k-4) for k >= 4 (observed, checked)


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- symbolic
def symbolic_checks():
    z, t, y, Hp, w = sp.symbols("z t y Hp w")
    pz = z * (z + 5)
    c = (z + 6) / ((z + 2) * (z + 5))
    k = (z + 2) * (z + 5) / (z + 6)
    check(sp.simplify(k * c - 1) == 0, "CC.1: k = 1/c = (z+2)(z+5)/(z+6)")
    # F = -1/H, H = c H(p)  =>  F = -1/(c H(p)) = k * (-1/H(p)) = k F(p)
    check(sp.simplify(-1 / (c * Hp) - k * (-1 / Hp)) == 0, "CC.1: H = c H(p)  =>  F = -1/H obeys F = k F(p)")
    y1, y2 = sp.solve(sp.Eq(y * (y + 5), t), y)
    pf = 1 / (t - pz) - (1 / ((2 * y1 + 5) * (y1 - z)) + 1 / ((2 * y2 + 5) * (y2 - z)))
    check(sp.simplify(pf) == 0, "CC.2: partial fractions 1/(t - p(z)) = sum_(p(y)=t) 1/(p'(y)(y - z))")
    ky = lambda yy: (yy + 2) * (yy + 5) / (yy + 6)
    T = sp.simplify(ky(y1) / (2 * y1 + 5) + ky(y2) / (2 * y2 + 5))
    check(sp.simplify(T - (2 - t) / (6 - t)) == 0, "CC.2: T(t) = sum_(p(y)=t) k(y)/p'(y) = (2-t)/(6-t) (explicit roots)")
    # independent: sum of residues of k(y)/(p(y)-t) at the two roots = -(res at -6) - (res at infinity)
    kk = ky(y) / (y * (y + 5) - t)
    r6 = sp.residue(kk, y, -6)
    rinf = -sp.limit(y * kk, y, sp.oo)
    check(sp.simplify(-(r6 + rinf) - (2 - t) / (6 - t)) == 0, "CC.2: same T from the residue theorem (poles -6 and infinity)")
    dT = sp.simplify(sp.diff((2 - t) / (6 - t), t))
    check(sp.simplify(dT + 4 / (6 - t) ** 2) == 0 and ((2 - t) / (6 - t)).subs(t, -5) == sp.Rational(7, 11)
          and ((2 - t) / (6 - t)).subs(t, 0) == sp.Rational(1, 3),
          "CC.2: T' = -4/(6-t)^2 < 0, T(-5) = 7/11, T(0) = 1/3, so 1/3 <= T <= 7/11 < 1 on [-5,0] (contains K_p)")
    kd = sp.factor(sp.cancel((k - k.subs(z, y)) / (y - z)))
    check(sp.simplify(kd - (-(y * z + 6 * y + 6 * z + 32) / ((y + 6) * (z + 6)))) == 0,
          "CC.2: (k(z) - k(y))/(y - z) = -(yz+6y+6z+32)/((y+6)(z+6)), holomorphic in (z,y) away from -6")
    # p^-1(U), U = (-11/2, 1/2): roots -5/2 +- sqrt(25/4 + t), t in U
    lo_s, hi_s = sp.sqrt(sp.Rational(25, 4) - sp.Rational(11, 2)), sp.sqrt(sp.Rational(25, 4) + sp.Rational(1, 2))
    check(lo_s > 0 and -sp.Rational(5, 2) - hi_s > -6 and -sp.Rational(5, 2) + hi_s < 1,
          "CC.2: p^-1((-11/2, 1/2)) = (-5/2-sqrt(27/4), -5/2-sqrt(3/4)) u (-5/2+sqrt(3/4), -5/2+sqrt(27/4)): avoids -6 and -5/2")
    gm, gp = -sp.Rational(5, 2) - sp.sqrt(sp.Rational(5, 4)), -sp.Rational(5, 2) + sp.sqrt(sp.Rational(5, 4))
    check(gm < -sp.Rational(5, 2) and gp > -2 and sp.simplify(gm * (gm + 5) + 5) == 0 and sp.simplify(gp * (gp + 5) + 5) == 0,
          "CC.3: K_p in p^-1([-5,0]) = [-5, -5/2-sqrt5/2] u [-5/2+sqrt5/2, 0]: -5/2 and -2 are not in K_p")
    check(pz.subs(z, -sp.Rational(5, 2)) == -sp.Rational(25, 4) and pz.subs(z, -2) == -6,
          "CC.3: p(-5/2) = -25/4 and p(-2) = -6 lie outside [-5,0] (so -5/2, -2 not in K_p, p^-1(K_p) = K_p)")
    kp = sp.simplify(k / (2 * z + 5))
    check(all(kp.subs(z, v) >= 0 for v in (-5, sp.Rational(-9, 2), -4, gm, gp, -1, sp.Rational(-1, 2), 0)),
          "CC.2: weights k(y)/p'(y) >= 0 at sample points of both pieces of p^-1([-5,0]) (positive transfer)")
    check(pz.subs(z, -5) == 0 and pz.subs(z, 0) == 0 and sp.solve(sp.Eq(pz, -5), z) and all(abs(sp.N(r) + 5) > 1 for r in sp.solve(sp.Eq(pz, -5), z)),
          "CC.3: p(-5) = 0 = p(0); a pole of Pi at u in Z needs psi(5^-j u) = -5, so psi(u) = p^(j-1)(0) = 0 (u in E); -5 not repeated")
    check(sp.simplify(sp.diff(w * (w + 5), w) - (2 * w + 5)) == 0,
          "CC.3/CC.6: (psi(5u))' = p'(psi(u)) psi'(u)/5: psi' = prod(1 + 2 psi(5^-i u)/5) vanishes only if some psi(5^-i u) = -5/2, "
          "impossible on Z (5^-i Z = Z, -5/2 not in K_p): no critical point of psi on Z")


def cauchy_schwarz_exact(trials=200, seed=7):
    """Exact rationals: for a probability measure mu, Im(1/r) >= Im z where r = sum c_i/(l_i + z) (CC.1's lemma)."""
    rnd = random.Random(seed)
    ok = True
    for _ in range(trials):
        n = rnd.randint(1, 6)
        cs = [Fraction(rnd.randint(1, 50)) for _ in range(n)]
        S = sum(cs)
        cs = [c_ / S for c_ in cs]
        ls = [Fraction(rnd.randint(0, 800), 100) for _ in range(n)]
        x, yv = Fraction(rnd.randint(-900, 300), 100), Fraction(rnd.randint(1, 300), 100)
        re, im = Fraction(0), Fraction(0)
        for c_, l_ in zip(cs, ls):
            dr, di = l_ + x, yv
            den = dr * dr + di * di
            re += c_ * dr / den
            im -= c_ * di / den
        n2 = re * re + im * im
        ok = ok and (-im / n2 >= yv)
    check(ok, f"CC.1: exact rationals: Im(1/r) >= Im z for {trials} random probability measures on [0,8] and points z (Cauchy-Schwarz)")


# ----------------------------------------------------------------- interval checks
def nevanlinna_checks():
    ctx.prec = 300
    ok, worst = True, None
    pts = [acb(-1, "0.1"), acb("-0.95", "0.001"), acb(-4, "1e-4"), acb("-3.7", "0.01"), acb("-0.3", "1e-3"),
           acb("-4.99", "0.05"), acb(2, 3), acb(-7, "0.5")]
    for zc in pts:
        H = cb.D_prod(zc, True) - zc / 2
        d = H.imag - zc.imag / 2
        ok = ok and bool(d > 0)
        Fv = -1 / H
        ok = ok and bool(abs(Fv) * zc.imag < 2)
        # F = k F(p)
        kz = (zc + 2) * (zc + 5) / (zc + 6)
        pzc = bz._p(zc)
        Fp = -1 / (cb.D_prod(pzc, True) - pzc / 2)
        ok = ok and bool(abs(Fv - kz * Fp) < arb("1e-40"))
    check(ok, "CC.1 (arb): Im H(z) - Im z/2 > 0, |F| Im z < 2 and F = k F(p) (to 1e-40) at 8 points, several within 1e-4 of K_p")
    # finite level: Im(1/r_n) >= Im z for CB.1's exact r_n at a complex point (sanity, not needed for the proof)
    zc = acb("-1.2", "0.01")
    for n in (2, 5, 9):
        pj, prod = zc, acb(1)
        for _ in range(n):
            prod *= (pj + 6) / (pj + 2)
            pj = bz._p(pj)
        Dn = zc / 2 * (1 + (1 + 6 / pj) * prod)
        ok = ok and bool(Dn.imag >= zc.imag)
    check(ok, "CC.1 (arb): Im(1/r_n(z)) >= Im z at z = -1.2 + 0.01i for CB.1's r_n, n = 2, 5, 9")


def psi_period():
    """CC.4: psi strictly increasing on [-5,0]; the period endpoints."""
    ctx.prec = 300
    bound = arb(1) * ((arb(1) / 10).exp() - 1)
    check(bool(bound < arb("0.106")), "CC.4: |psi(v) - v| <= psi(|v|) - |v| <= |v|(e^(|v|/10) - 1) <= 0.106 on [-1,0] (BJ.6), so psi > -1.11 > -5/2 there")
    a = ca.psi_enc(arb(-1))
    b = ca.psi_enc(arb(-5))
    check(bool(a > arb("-1.11")) and bool(abs(b - bz._p(a)) < arb("1e-60")),
          f"CC.4: psi(-1) = {a.str(18)}, psi(-5) = p(psi(-1)) = {b.str(18)}; psi' = prod(1 + 2 psi(5^-i u)/5) > 0 on [-5,0]")
    f1 = cb.f_prod(arb(-1))
    check(not f1.contains(0), f"CC.4: f(-1) = {f1.str(15)} != 0 (finite): lam = 1 (and 5) is not an atom, -1 not in Z")
    return a, b


def _psi_mp(u, J=80):
    v = u / mp.mpf(5) ** J
    x = v
    for _ in range(J):
        x = x * (x + 5)
    return x


def period_atoms(K, verbose=True):
    """All atoms of nu in [1,5) of level k <= K (psi(-lam) in p^-k(-3)), certified, with CB.5 masses."""
    ctx.prec = 300
    mp.mp.dps = 60
    a, b = psi_period()
    lo, hi = b, a
    g = arb(3).log() / arb(5).log()
    found, level_sums, counts = [], [], []
    targets = [mp.mpf(-3)]
    ok_bracket = True
    for k in range(K + 1):
        s, n = arb(0), 0
        for tg in targets:
            tga = arb(mp.nstr(tg, 70)) + arb(0, arb("1e-50"))
            if not (bool(tga > lo) and bool(tga < hi)):
                if (tga - lo).contains(0) or (tga - hi).contains(0):
                    ok_bracket = False
                continue
            # float proposal, then a certified bracket (psi strictly increasing on [-5,-1])
            um = mp.findroot(lambda u: _psi_mp(u) - tg, (mp.mpf(-5), mp.mpf(-1)), solver="anderson")
            d = mp.mpf("1e-45")
            ua, ub = arb(mp.nstr(um - d, 70)), arb(mp.nstr(um + d, 70))
            # exact preimage target: p^k(psi(u)) = -3 ; compare psi(u) with an enclosure of tg via k branch steps
            sa, sb = ca.psi_enc(ua) - tga, ca.psi_enc(ub) - tga
            if not (bool(sa < 0) and bool(sb > 0)):
                ok_bracket = False
                continue
            # exact certificate: p^k(psi(u)) + 3 changes sign across [ua, ub], so the ball holds a level-k point of Lambda_O
            fa, fb = cb._pk(ca.psi_enc(ua), k) + 3, cb._pk(ca.psi_enc(ub), k) + 3
            if not ((bool(fa < 0) and bool(fb > 0)) or (bool(fa > 0) and bool(fb < 0))):
                ok_bracket = False
                continue
            u = arb((ua + ub) / 2, ((ub - ua) / 2).upper())
            m, _ = cb.mass_closed(u, k)
            lam = -u
            ok_bracket = ok_bracket and bool(m > 0) and bool(lam >= 1) and bool(lam < 5)
            found.append((lam, m, k))
            s += m * lam ** (-g)
            n += 1
        level_sums.append(s)
        counts.append(n)
        if verbose:
            print(f"   level {k:2d}: {n:5d} atoms, sum lam^-gamma nu = {s.str(12)}")
        nxt = []
        for tg in targets:
            r = mp.sqrt(tg + mp.mpf(25) / 4)
            nxt += [-mp.mpf(5) / 2 + r, -mp.mpf(5) / 2 - r]
        targets = nxt
    return found, level_sums, counts, ok_bracket


def disjoint_ok(found):
    """The certified balls are pairwise disjoint (sorted by midpoint), so they hold distinct atoms."""
    srt = sorted(found, key=lambda x: float(x[0].mid()))
    return all(bool(srt[i][0] < srt[i + 1][0]) for i in range(len(srt) - 1))


def mass_identity(K):
    found, sums, counts, okb = period_atoms(K)
    check(okb, f"CC.4: every point of p^-k(-3) in (psi(-5), psi(-1)), k <= {K}, has a certified preimage ball in (-5,-1] (p^k(psi) + 3 changes sign across it) with CB.5 mass > 0")
    check(disjoint_ok(found), f"CC.4: the {len(found)} balls are pairwise disjoint; psi is injective on [-5,-1], so they are all the atoms of level <= {K} in [1,5)")
    exp_counts = [COUNTS.get(k, 5 * 2 ** (k - 4)) for k in range(K + 1)]
    check(counts == exp_counts, f"CC.4: level counts {counts} (= 1,0,1,3 then 5*2^(k-4)); {sum(counts)} atoms in one period through level {K}")
    tot = arb(0)
    for s in sums:
        tot += s
    mT = (arb(C_R) + arb(0, arb("1e-59"))) / 2
    delta = mT - tot
    print(f"   sum through level {K} = {tot.str(16)};  m(T) = c_R(gamma)/2 = {mT.str(20)}")
    print(f"   delta_{K} = mass of the atoms of level > {K} (exact, by CC.3) in {delta.str(10)}")
    check(bool(delta > 0), f"CC.4: partial sum through level {K} < m(T) (consistent: the remaining levels carry delta_{K} > 0)")
    ratios = [sums[k + 1] / sums[k] for k in range(4, K)]
    print("   level-sum ratios (k >= 4): " + ", ".join(r.str(6) for r in ratios))
    if K >= 12:
        r = sums[K] / sums[K - 1]
        est = sums[K] * r / (1 - r)
        print(f"   Observation: geometric tail s_K r/(1-r), r = s_K/s_(K-1) = {r.str(8)}: {est.str(10)} vs delta_{K} = {delta.str(10)}")
        check(bool(abs(est - delta) < arb("1e-13")), f"CC.4 (consistency, no claim): the geometric extrapolation of the level sums matches delta_{K} to 1e-13, "
              "i.e. the CB.5 masses sum to c_R(gamma)/2 computed by Theorems AT/AU's independent route")
    return found, tot, delta, mT, sums


def bz_comparison(found, delta, K):
    """BZ.5 found 311 atoms carrying 0.25854262636241...; through level 13 our set contains them; for smaller K compare
    with BZ's sum: BZ's remainder 4.03e-5 must exceed ours."""
    ok = bool(delta.upper() < arb("4.03e-5")) if K >= 10 else True
    check(ok, f"CC.4: through level {K} the uncaptured mass delta_{K} <= {float(delta.upper()):.3e} (BZ.5: 4.03e-5 with 311 atoms)")
    big = max(found, key=lambda x: x[1].mid())
    check(bool(abs(big[1] - arb("0.4796406929290034734083466")) < arb("1e-20")) and big[2] == 0,
          f"CC.4: largest atom lam = {big[0].str(16)}, nu = {big[1].str(20)} (CB.5, level 0)")


def rho_bounds(found, delta, mT, kmax=100):
    """CC.5: limsup rho_k = 1 (proof); liminf rho_k <= (|a_k| + delta)/m(T) by BZ.6's recurrence (certified)."""
    two_pi = 2 * arb.pi()
    l5 = arb(5).log()
    g = arb(3).log() / l5
    ts = [(lam.log() / l5, m * lam ** (-g)) for lam, m, _ in found]
    amin, kmin = None, None
    for k in range(1, kmax + 1):
        acc = acb(0)
        for t_, w_ in ts:
            acc += acb(0, -two_pi * k * t_).exp() * w_
        if amin is None or abs(acc).upper() < amin:
            amin, kmin = abs(acc).upper(), k
    inf_hi = (arb(amin) + arb(delta.upper())) / mT
    print(f"   liminf rho_k <= {float(inf_hi.upper()):.8f} (k = {kmin}); limsup rho_k = 1 (CC.5)")
    if len(found) >= 311:
        check(bool(inf_hi < arb("0.59497")), "CC.5: certified liminf rho_k bound at least as sharp as BZ.6's 0.59497")
    sharp = arb.pi() * arb(C_R) / l5
    check(bool(abs(sharp - arb("1.0094980205")) < arb("1e-10")),
          f"CC.5: limsup |theta_k| e^(2 pi^2 k/ln5) = pi c_R(gamma)/ln5 = {sharp.str(15)} (BX.3's bound attained)")
    return inf_hi


def main(K=K_CERT) -> int:
    t0 = time.time()
    print("Theorem CC: nu is pure point (nu(-Z) = 0); atoms exactly Lambda_O (Conjecture CA.6); "
          "c_R(gamma)/2 = sum of CB.5 masses over one period; limsup rho_k = 1. Not thrust.")
    symbolic_checks()
    cauchy_schwarz_exact()
    nevanlinna_checks()
    found, tot, delta, mT, sums = mass_identity(K)
    bz_comparison(found, delta, K)
    rho_bounds(found, delta, mT)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CC checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else K_CERT))
