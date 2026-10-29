#!/usr/bin/env python3
r"""Theorem CE: the level sums of BX's spectral measure over one scale period have an exact geometric ratio.
s_(k+1)/s_k -> r = 0.42373128696260793985253216132843923...,  the leading eigenvalue of an explicit positive
transfer operator that does not involve psi; r is certified to 2e-38 by an exact rational Collatz-Wielandt sandwich.

Setting (GASKET.md, Theorems CB, CC, CD): p(z) = z(z+5), O_k = p^-k(-3), a = psi(-1), I = (p(a), a]; the atoms of nu in
[1,5) are lam_t = -psi^-1(t), t in O n I (CC.4), with CB.5 masses; s_k = sum_(t in O_k n I) lam_t^-gamma nu({lam_t})
(CC.4's level sums, sum_k s_k = m(T) = c_R(gamma)/2).  X = [-5,0]; the inverse branches of p on X are
phi_(1,2)(t) = -5/2 -+ s, s = sqrt(t + 25/4) in [sqrt5/2, 5/2].

Theorem CE (checked here; exit 1 on any failure):
  CE.1 [proved from CB.5; elementary] For t in O put sigma_t = 1/H'(t) (H = (z/2)P, CB.2), i.e. CB.5's mass with the
       u-dependent factors Pi(u) psi'(u) removed:  nu({lam_t}) = sigma_t / (Pi(u_t) psi'(u_t)), u_t = psi^-1(t).  Then
           sigma_y = w(y) sigma_(p(y))   (y in O_(k+1)),   w(y) = k(y)/p'(y) = (y+2)(y+5)/((y+6)(2y+5)),
       sigma_-3 = (8/9)/prod_(j>=1) F(p^j(6)) > 0 (F(x) = (x+6)/(x+2)), and so s_k = sigma_-3 (L^k G0)(-3) with
       (L h)(t) = sum_(p(y)=t) w(y) h(y) and G0 = 1_I g, g(t) = lam_t^-gamma / (Pi(u_t) psi'(u_t)) > 0 on [p(a), a].
  CE.2 [proved; exact] With q(t) = -t, L(q f) = q Lt f, where Lt f(t) = wt1(s) f(phi_1 t) + wt2(s) f(phi_2 t),
           wt1 = (s+1/2)/(2s(7/2-s)),  wt2 = (s-1/2)/(2s(s+7/2)),  wt1 + wt2 = 4/(6-t),
       both real-analytic and positive on X (removable at t = 0: wt1(0) = 3/5, wt2(0) = 1/15).
       s_k = 3 sigma_-3 (Lt^k G)(-3), G = G0/q.
  CE.3 [proved from CE.2, CD.2 and Assumption CE (Birkhoff)] Exactly 5 of the 16 level-4 branch compositions Phi have
       Phi(X) inside I and the other 11 miss [p(a), a]; so Lt^4 G is positive and C^1 on X.  Lt maps the cone
       C_b = {f > 0 : f(x) <= e^(b|x-y|) f(y)} into C_(M + b/sqrt5) (M = max Lip log wt_i <= 1.23, |phi_i'| <= 1/sqrt5),
       hence C_a into C_(lam a) for lam in (1/sqrt5, 1), a = M/(lam - 1/sqrt5); diam_(C_a) C_(lam a) <=
       2 log((1+lam)/(1-lam)) + 10 lam a (proved here).  Birkhoff's theorem gives a positive eigenfunction ht in C_a,
       Lt ht = r ht, and Lt^(n+1) f / Lt^n f -> r uniformly on X for every f in C_a; therefore
           s_(k+1)/s_k -> r,   and the same limit holds for every positive C^1 weight on [p(a), a] in place of g
       (e.g. the plain nu-masses of the period).  Consequences: delta_K/s_K -> r/(1-r) (delta_K = sum_(k>K) s_k), and the
       mean lam^-gamma-mass of a level-k atom, s_k/(5*2^(k-4)) (CD.3), has ratio -> r/2.
  CE.4 [certified; exact rationals] For an explicit polynomial f > 0 on X (degree 48, rational coefficients):
           lo f <= Lt f <= hi f  on X,   hi - lo = 2e-38,
       checked as positivity of two degree-99 rational polynomials in s by Bernstein coefficients (exact Fractions,
       with subdivision).  Hence (Collatz-Wielandt, proved) r in [lo, hi]:
           r = 0.42373128696260793985253216132843923...  (CC's observed 0.42373129...).
  CE.5 [certified consistency, no claim] CC's certified level sums (arb, CB.5 masses) give s_(k+1)/s_k - r with
       alternating sign and |.| -> 0 (1.6e-8 at k = 11); the successive deviations have ratio near -0.1725 (the next
       eigenvalue of Lt, about -0.0731, over r; floats only).  delta_14/s_14 matches r/(1-r).
  Kept Failure CE.6 [proved] The untwisted operator L is not a strictly positive Ruelle operator on C(X): w(-5) = 0 and
       p^-1(0) = {0,-5}, so (L h)(0) = h(0)/3 for every h; 1/3 is an eigenvalue (left eigenvector delta_0), L^k 1(0) = 3^-k,
       and every r-eigenfunction of L vanishes at 0.  The un-normalised Herglotz masses L^k 1(-3) = sigma(O_k)/sigma_-3
       have ratio 0.4261... at k = 14, still 2.4e-3 from r (the 1/3-mode decays only like (1/(3r))^k = 0.787^k).
       The twist q = -t is what makes the cone argument work.
Assumption CE (literature): Birkhoff's contraction theorem: a positive linear map sending a convex cone C into a subcone
of finite Hilbert-metric diameter Delta contracts the Hilbert metric of C by tanh(Delta/4).
Assumptions used as stated (not re-run): CB.2, CB.5 (closed product masses), CC.4 (one period parametrised by psi on
[-5,-1]; g > 0 analytic there), CD.2-CD.3, and the enclosure routines of BZ, CA, CB, CC (imported unchanged).
No closed form of r is claimed (OPEN).  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW enters.
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr
from math import comb

import mpmath as mp
import sympy as sp
from flint import arb, ctx, fmpq, fmpq_poly

import gasket_fractional_profile_atoms as bz
import gasket_fractional_profile_atom_law as ca
import gasket_fractional_profile_product as cb

FAILS: list[str] = []
DEG = 48                 # degree of the test function f (in x = (t + 5/2)/(5/2))
EPS = Fr(1, 10 ** 38)     # half-width of the certified enclosure of r
K_CERT = 13              # CC level sums compared in `reproduce` (about 20 s)
K_FAST = 8               # CC level sums compared in the pytest run
R_PRINT = "0.42373128696260793985253216132843923"


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- symbolic (CE.1, CE.2, CE.3, CE.6)
def symbolic_checks():
    y, t, s, lam, x = sp.symbols("y t s lam x", real=True)
    p = y * (y + 5)
    k = (y + 2) * (y + 5) / (y + 6)
    w = sp.together(k / sp.diff(p, y))
    check(sp.simplify(w - (y + 2) * (y + 5) / ((y + 6) * (2 * y + 5))) == 0, "CE.1: w = k/p' = (y+2)(y+5)/((y+6)(2y+5)), k = 1/c")
    # sigma law from H = c H(p): at a zero y of H(p(y)), H'(y) = c(y) p'(y) H'(p(y))
    Hf = sp.Function("H")
    c = (y + 6) / ((y + 2) * (y + 5))
    dH = sp.diff(c * Hf(p), y).subs(Hf(p), 0)
    check(sp.simplify(dH - c * sp.diff(p, y) * sp.Subs(sp.Derivative(Hf(x), x), x, p).doit()) == 0,
          "CE.1: H = c H(p) and H(p(y)) = 0 give H'(y) = c(y) p'(y) H'(p(y)), so sigma_y = 1/H'(y) = w(y) sigma_(p(y))")
    # the same law read off CB.5's closed product: sigma_t = -4/((t/2) prod_(i != k+1) F(x_i) prod_(i<=k) (2x_i+5))
    F = lambda z: (z + 6) / (z + 2)
    ratio = (p / 2) / ((y / 2) * F(y) * (2 * y + 5))
    check(sp.simplify(ratio - w) == 0, "CE.1: CB.5's product: sigma_y/sigma_(p(y)) = (p(y)/2)/((y/2) F(y) (2y+5)) = w(y) (orbit shift x'_i = x_(i-1))")
    # twisted weights
    ph1, ph2 = -sp.Rational(5, 2) - s, -sp.Rational(5, 2) + s
    tt = s ** 2 - sp.Rational(25, 4)
    w1 = sp.simplify(w.subs(y, ph1) * ph1 / tt)
    w2 = sp.simplify(w.subs(y, ph2) * ph2 / tt)
    W1 = (s + sp.Rational(1, 2)) / (2 * s * (sp.Rational(7, 2) - s))
    W2 = (s - sp.Rational(1, 2)) / (2 * s * (s + sp.Rational(7, 2)))
    check(sp.simplify(w1 - W1) == 0 and sp.simplify(w2 - W2) == 0,
          "CE.2: w(phi_i) q(phi_i)/q(t) = (s+1/2)/(2s(7/2-s)) and (s-1/2)/(2s(s+7/2)), q = -t, t = s^2 - 25/4 (removable at t = 0)")
    check(sp.simplify(W1 + W2 - 4 / (6 - tt)) == 0, "CE.2: wt1 + wt2 = 4/(6-t) (so Lt 1 in [4/11, 2/3] on X)")
    check(W1.subs(s, sp.Rational(5, 2)) == sp.Rational(3, 5) and W2.subs(s, sp.Rational(5, 2)) == sp.Rational(1, 15),
          "CE.2: wt1(0) = 3/5, wt2(0) = 1/15; on s in [sqrt5/2, 5/2] all factors s +- 1/2, s, 7/2 -+ s are positive")
    check(sp.simplify(sp.diff(ph1, s) / (2 * s) + 1 / (2 * s)) == 0 and sp.simplify(sp.diff(ph2, s) / (2 * s) - 1 / (2 * s)) == 0,
          "CE.3: phi_i'(t) = -+1/(2s), so |phi_i'| <= 1/sqrt5 on X")
    # Kept Failure CE.6
    check(w.subs(y, -5) == 0 and w.subs(y, 0) == sp.Rational(1, 3) and set(sp.solve(p, y)) == {0, -5},
          "CE.6: p^-1(0) = {0,-5}, w(0) = 1/3, w(-5) = 0: (L h)(0) = h(0)/3 for every h (delta_0 L = delta_0/3)")
    # cone-diameter lemma: Fq(x) = (1+lam)(e^x - e^(lam x)) - (1-lam)(e^x - e^(-lam x)) >= 0 for x >= 0
    Fq = (1 + lam) * (sp.exp(x) - sp.exp(lam * x)) - (1 - lam) * (sp.exp(x) - sp.exp(-lam * x))
    dF = sp.expand(sp.diff(Fq, x))
    check(Fq.subs(x, 0) == 0 and sp.simplify(dF - lam * (2 * sp.exp(x) - (1 + lam) * sp.exp(lam * x) - (1 - lam) * sp.exp(-lam * x))) == 0,
          "CE.3: cone lemma: Fq(0) = 0, Fq' = lam(2e^x - (1+lam)e^(lam x) - (1-lam)e^(-lam x)) >= 0 (e^(+-lam x) <= e^x), "
          "so (e^x - e^(lam x))/(e^x - e^(-lam x)) >= (1-lam)/(1+lam)")
    return W1, W2


def cone_constants():
    """CE.3: explicit Lipschitz bound M for log wt_i (in t) and the cone parameters (qualitative use only)."""
    ctx.prec = 200
    s0 = arb(5).sqrt() / 2
    kap = 1 / (2 * s0)
    # d/ds log wt1 = 1/(s+1/2) - 1/s + 1/(7/2-s); d/ds log wt2 = 1/(s-1/2) - 1/s - 1/(s+7/2); dt = 2s ds
    M1 = kap * (1 / (s0 + arb(1) / 2) + 1 / s0 + 1)
    M2 = kap * (1 / (s0 - arb(1) / 2) + 1 / s0 + 1 / (s0 + arb(7) / 2))
    M = max(M1.upper(), M2.upper())
    lamc = (1 + kap) / 2
    a = arb(M) / (lamc - kap)
    Delta = 2 * ((1 + lamc) / (1 - lamc)).log() + 2 * lamc * a * 5
    check(bool(arb(M) < arb("1.23")) and bool(kap < 1),
          f"CE.3: Lip(log wt_i) <= M = {float(arb(M)):.4f} < 1.23 on X, |phi_i'| <= 1/sqrt5 = {kap.str(6)}; Lt C_b c C_(M + b/sqrt5)")
    print(f"   cone: lam = {lamc.str(6)}, a = M/(lam - 1/sqrt5) = {a.str(6)}, diam <= {Delta.str(6)}, "
          f"Birkhoff factor tanh(diam/4) = 1 - {(1 - (Delta / 4).tanh()).str(4)} < 1 (qualitative; CE.4 gives the digits)")
    return M


# ----------------------------------------------------------------- CE.1 sigma law on the preimage tree (arb)
def sigma_closed(t, k):
    """sigma_t = 1/H'(t) for t in O_k from CB's closed product (t-orbit only; no psi)."""
    orb = [t]
    for _ in range(k + 1):
        orb.append(bz._p(orb[-1]))
    rest = arb(1)
    for i in range(k + 1):
        rest *= (orb[i] + 6) / (orb[i] + 2)
    rest /= (orb[k + 1] + 2)
    rest *= cb.P_enc(bz._p(orb[k + 1]))
    d = arb(1)
    for i in range(k + 1):
        d *= 2 * orb[i] + 5
    return 1 / (t / 2 * rest * d)


def _w(y):
    return (y + 2) * (y + 5) / ((y + 6) * (2 * y + 5))


def sigma_law(K):
    ctx.prec = 300
    q = arb(25) / 4
    h = arb(5) / 2
    s3 = sigma_closed(arb(-3), 0)
    prod6 = cb.P_enc(arb(6))
    check(bool(s3 > 0) and bool(abs(s3 - arb(8) / 9 / prod6) < arb("1e-60")),
          f"CE.1: sigma_-3 = (8/9)/prod_(j>=0) F(p^j 6) = {s3.str(20)} > 0")
    lev = [(arb(-3), s3)]
    ok, worst = True, arb(0)
    for k in range(1, K + 1):
        nxt = []
        for tt, st in lev:
            r = (tt + q).sqrt()
            for yy in (-h - r, -h + r):
                sy = sigma_closed(yy, k)
                d = abs(sy / (_w(yy) * st) - 1)
                worst = max(worst, d.upper())
                ok = ok and bool(sy > 0) and bool(d < arb("1e-40"))
                nxt.append((yy, sy))
        lev = nxt
    check(ok, f"CE.1 (arb): sigma_y = w(y) sigma_(p(y)) > 0 for all y in O_k, 1 <= k <= {K} (closed product vs law, rel. err <= {float(arb(worst)):.1e})")
    return s3


def mass_factorisation(found_k):
    """CB.5's mass = sigma_t/(Pi(u) psi'(u)), sigma_t computed from t alone, for the certified atoms supplied."""
    ctx.prec = 300
    ok = True
    for lam, m, k in found_k:
        u = -lam
        x, eps, ball, _ = cb.orbit(u)
        t = x[0]
        Pi = cb.Pi_enc(x, eps, ball)
        sg = sigma_closed(t, k)
        ok = ok and bool(abs(sg / (Pi * cb.dpsi_enc(u)) / m - 1) < arb("1e-25"))
    check(ok, f"CE.1 (arb): nu({{lam_t}}) = sigma_t/(Pi(u_t) psi'(u_t)) for {len(found_k)} certified atoms (sigma_t from t's orbit only)")


# ----------------------------------------------------------------- CE.3 the level-4 branches
def branches4():
    ctx.prec = 300
    a = ca.psi_enc(arb(-1))
    pa = bz._p(a)
    q = arb(25) / 4
    h = arb(5) / 2
    ims = [(arb(-5), arb(0))]
    for _ in range(4):
        nxt = []
        for lo, hi in ims:
            for sg in (-1, 1):
                e1, e2 = -h + sg * (lo + q).sqrt(), -h + sg * (hi + q).sqrt()
                nxt.append((min(e1, e2, key=lambda z: float(z.mid())), max(e1, e2, key=lambda z: float(z.mid()))))
        ims = nxt
    inside = sum(1 for lo, hi in ims if bool(lo > pa) and bool(hi < a))
    outside = sum(1 for lo, hi in ims if bool(hi < pa) or bool(lo > a))
    check(len(ims) == 16 and inside == 5 and outside == 11,
          f"CE.3 (arb): of the 16 level-4 branch images Phi(X), {inside} lie inside I = (p(a), a) and {outside} miss [p(a), a]; "
          "so Lt^4 G = sum over 5 branches of positive analytic terms: positive and C^1 on X")


# ----------------------------------------------------------------- CE.4 the exact Collatz-Wielandt certificate
def eigen_approx(deg=DEG, dps=70):
    """Float proposal only: power iteration for Lt on Chebyshev nodes (mpmath); returns the Chebyshev coefficients."""
    mp.mp.dps = dps
    N = deg + 1
    xs = [mp.cos(mp.pi * (j + mp.mpf(1) / 2) / N) for j in range(N)]

    def coeffs(vals):
        cs = [2 * sum(vals[j] * mp.cos(mp.pi * n * (j + mp.mpf(1) / 2) / N) for j in range(N)) / N for n in range(N)]
        cs[0] /= 2
        return cs

    def clen(cs, x):
        b1 = b2 = mp.mpf(0)
        for a_ in reversed(cs[1:]):
            b1, b2 = 2 * x * b1 - b2 + a_, b1
        return x * b1 - b2 + cs[0]

    vals, rr = [mp.mpf(1)] * N, None
    for _ in range(80):
        cs = coeffs(vals)
        new = []
        for x in xs:
            s = mp.sqrt(-mp.mpf(5) / 2 + mp.mpf(5) / 2 * x + mp.mpf(25) / 4)
            w1 = (s + mp.mpf(1) / 2) / (2 * s * (mp.mpf(7) / 2 - s))
            w2 = (s - mp.mpf(1) / 2) / (2 * s * (s + mp.mpf(7) / 2))
            new.append(w1 * clen(cs, -2 * s / 5) + w2 * clen(cs, 2 * s / 5))
        r_new = new[N // 2] / vals[N // 2]
        vals = [v / new[N // 2] for v in new]
        if rr is not None and abs(r_new - rr) < mp.mpf(10) ** -50:
            rr = r_new
            break
        rr = r_new
    return coeffs(vals), rr


def bernstein_positive(P, a, b, depth=0, maxdepth=14):
    """Exact test: P > 0 on [a, b] if all Bernstein coefficients of P(a + (b-a)u) on [0,1] are > 0 (with subdivision).
    Returns the number of leaf intervals, or None if not certified."""
    Q = P(fmpq_poly([a, b - a]))
    n = Q.degree()
    qc = [Fr(int(c.p), int(c.q)) for c in Q.coeffs()] + [Fr(0)] * (n + 1 - len(Q.coeffs()))
    bc = [sum(Fr(comb(k, i), comb(n, i)) * qc[i] for i in range(k + 1)) for k in range(n + 1)]
    if all(v > 0 for v in bc):
        return 1
    if bc[0] <= 0 or bc[-1] <= 0 or depth >= maxdepth:
        return None
    m = (a + b) / 2
    left = bernstein_positive(P, a, m, depth + 1, maxdepth)
    if left is None:
        return None
    right = bernstein_positive(P, m, b, depth + 1, maxdepth)
    return None if right is None else left + right


def _q(v):
    return fmpq(*Fr(mp.nstr(v, 60)).as_integer_ratio())


def certificate(deg=DEG, eps=EPS):
    cs, rr = eigen_approx(deg)
    X = fmpq_poly([0, 1])
    T = [fmpq_poly([1]), X]
    for _ in range(2, deg + 1):
        T.append(2 * X * T[-1] - T[-2])
    f = fmpq_poly(0)
    for n in range(deg + 1):
        f += _q(cs[n]) * T[n]
    S = fmpq_poly([0, 1])
    half, sev = fmpq(1, 2), fmpq(7, 2)
    # x = (t + 5/2)/(5/2): phi_1 -> -2s/5, phi_2 -> 2s/5, t = s^2 - 25/4 -> (s^2 - 15/4) 2/5
    Ns = (S + half) * (S + sev) * f(-S * fmpq(2, 5)) + (S - half) * (sev - S) * f(S * fmpq(2, 5))
    Ds = 2 * S * (fmpq(49, 4) - S * S) * f((S * S - fmpq(15, 4)) * fmpq(2, 5))
    # Lt f / f = Ns/Ds (common denominator 2s(7/2-s)(7/2+s) > 0)
    mp.mp.dps = 60
    tchk = [mp.mpf(v) for v in ("-5", "-3.7", "-2", "-0.4", "0")]
    ok_id = True
    for tv in tchk:
        s = mp.sqrt(tv + mp.mpf(25) / 4)
        fx = lambda tt: mp.mpf(sum(mp.mpf(int(c.p)) / int(c.q) * ((tt + mp.mpf(5) / 2) / (mp.mpf(5) / 2)) ** i for i, c in enumerate(f.coeffs())))
        w1 = (s + mp.mpf(1) / 2) / (2 * s * (mp.mpf(7) / 2 - s))
        w2 = (s - mp.mpf(1) / 2) / (2 * s * (s + mp.mpf(7) / 2))
        lhs = (w1 * fx(-mp.mpf(5) / 2 - s) + w2 * fx(-mp.mpf(5) / 2 + s)) / fx(tv)
        ev = lambda P: mp.mpf(sum(mp.mpf(int(c.p)) / int(c.q) * s ** i for i, c in enumerate(P.coeffs())))
        ok_id = ok_id and abs(lhs - ev(Ns) / ev(Ds)) < mp.mpf(10) ** -40
    check(ok_id, "CE.4: Lt f/f = Ns(s)/Ds(s) (rational polynomials in s) at 5 points of X (floats, consistency of the algebra)")
    fpos = bernstein_positive(f, fmpq(-1), fmpq(1))
    check(fpos is not None, f"CE.4 (exact): f > 0 on X (Bernstein coefficients of the degree-{f.degree()} f, {fpos} piece(s))")
    rq = Fr(mp.nstr(rr, 50))
    lo = fmpq(*(rq - eps).as_integer_ratio())
    hi = fmpq(*(rq + eps).as_integer_ratio())
    a, b = fmpq(1118, 1000), fmpq(5, 2)          # [1.118, 5/2] contains s(X) = [sqrt5/2, 5/2]
    pl = bernstein_positive(Ns - lo * Ds, a, b)
    ph = bernstein_positive(hi * Ds - Ns, a, b)
    check(pl is not None and ph is not None,
          f"CE.4 (exact): lo f <= Lt f <= hi f on X: Ns - lo Ds > 0 and hi Ds - Ns > 0 on s in [1.118, 5/2] "
          f"(degree {Ns.degree()}/{Ds.degree()}, Bernstein, {pl}+{ph} pieces), hi - lo = {float(Fr(int((hi - lo).p), int((hi - lo).q))):.1e}")
    r_lo, r_hi = Fr(int(lo.p), int(lo.q)), Fr(int(hi.p), int(hi.q))
    check(Fr(R_PRINT) <= r_lo and r_hi < Fr(R_PRINT) + Fr(1, 10 ** 35),
          f"CE.4: r in [{mp.nstr(mp.mpf(r_lo.numerator) / r_lo.denominator, 40)}, {mp.nstr(mp.mpf(r_hi.numerator) / r_hi.denominator, 40)}]")
    check(Fr(4, 11) < r_lo and r_hi < Fr(2, 3), "CE.4: consistent with the crude bound 4/11 <= Lt 1 <= 2/3")
    return r_lo, r_hi


# ----------------------------------------------------------------- CE.5 consistency with CC's certified level sums
def level_sum_consistency(K, r_lo, r_hi):
    import gasket_fractional_profile_pure_point as cc
    found, sums, counts, okb = cc.period_atoms(K, verbose=False)
    check(okb, f"CE.5: CC's certified atoms of level <= {K} recomputed ({len(found)} atoms)")
    ctx.prec = 300
    r = arb(Fr((r_lo + r_hi) / 2).numerator) / Fr((r_lo + r_hi) / 2).denominator
    devs = [sums[k + 1] / sums[k] - r for k in range(4, K)]
    print("   s_(k+1)/s_k - r, k = 4..%d: %s" % (K - 1, ", ".join(d.str(4) for d in devs)))
    alt = all(bool(devs[i] * devs[i + 1] < 0) for i in range(len(devs) - 1))
    dec = all(bool(abs(devs[i + 1]) < abs(devs[i]) / 4) for i in range(len(devs) - 1))
    check(alt and dec, f"CE.5 (consistency, no claim): the deviations alternate in sign and shrink by a factor > 4 per level")
    q = [float((devs[i + 1] / devs[i]).mid()) for i in range(len(devs) - 1)]
    print(f"   successive deviation ratios (floats): {', '.join(f'{v:.4f}' for v in q)}")
    if K >= 13:
        check(bool(abs(devs[-1]) < arb("1e-8")), f"CE.5: |s_{K}/s_{K-1} - r| < 1e-8")
    mT = (arb(bz.C_R) + arb(0, arb("1e-59"))) / 2
    tot = arb(0)
    for sk in sums:
        tot += sk
    delta = mT - tot
    dr = delta / sums[K]
    check(bool(abs(dr - r / (1 - r)) < arb("0.02") * (arb("0.18") ** (K - 4))),
          f"CE.5 (consistency): delta_{K}/s_{K} = {dr.str(10)} vs r/(1-r) = {(r / (1 - r)).str(10)}")
    return found


# ----------------------------------------------------------------- Kept Failure CE.6
def naive_failure(K=14, r_mid=None):
    ctx.prec = 200
    q = arb(25) / 4
    h = arb(5) / 2

    def L1(t, n):
        if n == 0:
            return arb(1)
        s = (t + q).sqrt()
        y1, y2 = -h - s, -h + s
        return _w(y1) * L1(y1, n - 1) + _w(y2) * L1(y2, n - 1)

    ok0 = all(bool(abs(L1(arb(0), n) - arb(3) ** (-n)) < arb("1e-40")) for n in range(1, 11))
    check(ok0, "CE.6: L^n 1(0) = 3^-n for n <= 10 (arb; the -5 branch has weight 0)")
    vK, vK1 = L1(arb(-3), K - 1), L1(arb(-3), K)
    rat = vK1 / vK
    check(bool(abs(rat - arb(r_mid)) > arb("1e-3")),
          f"Kept Failure CE.6: untwisted Herglotz masses: L^{K} 1(-3)/L^{K-1} 1(-3) = {rat.str(8)}, still > 1e-3 away from r")


def main(K=K_CERT) -> int:
    t0 = time.time()
    print("Theorem CE: the level sums of nu over one scale period have the exact ratio r = spectral radius of an explicit "
          "positive transfer operator; r certified to 2e-38. Not thrust.")
    symbolic_checks()
    cone_constants()
    sigma_law(8)
    branches4()
    r_lo, r_hi = certificate()
    found = level_sum_consistency(K, r_lo, r_hi)
    sample = [f_ for f_ in found if f_[2] <= 6][:12]
    mass_factorisation(sample)
    naive_failure(14, float((r_lo + r_hi) / 2))
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CE checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else K_CERT))
