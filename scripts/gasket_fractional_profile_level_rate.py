#!/usr/bin/env python3
r"""Theorem CF: Assumption CE (Birkhoff's contraction theorem) is discharged, with an explicit contraction rate.

Setting (GASKET.md, Theorem CE): X = [-5,0], s = sqrt(t + 25/4) in [sqrt5/2, 5/2], the twisted transfer operator
    (Lt f)(t) = wt1(s) f(phi_1 t) + wt2(s) f(phi_2 t),  phi_(1,2)(t) = -5/2 -+ s,
    wt1 = (s+1/2)/(2s(7/2-s)),  wt2 = (s-1/2)/(2s(s+7/2))   (strictly positive on X, CE.2),
and CE's level sums s_k = 3 sigma_-3 (Lt^k G)(-3) with Lt^4 G positive and C^1 on X (CE.1-CE.3).  CE.3 obtained
s_(k+1)/s_k -> r from Birkhoff's theorem, cited as Assumption CE (literature).  This script removes that assumption.

Theorem CF (checked here; exit 1 on any failure):
  CF.1 [proved; elementary] (Birkhoff contraction, the form CE.3 used.)  Let C be a convex cone of functions whose
       closure is C u {0}, T a positive linear map with T(C) in C, and Theta the Hilbert metric of C
       (Theta(f,g) = log(beta/alpha), alpha = max{t: g - t f in C u {0}}, beta = min{t: t f - g in C u {0}}).
       If diam_Theta T(C) <= Delta < oo then Theta(Tf, Tg) <= tanh(Delta/4) Theta(f, g) for all f, g in C.
       Proof: two linear identities and one inequality (z+m)(z+M) - z(sqrt m + sqrt M)^2 = (z - sqrt(mM))^2 >= 0
       (all checked by sympy).
  CF.2 [proved; arb certified] In the variable s the branch maps are s -> sqrt(15/4 -+ s), with |d/ds| <= c_1 = 1/sqrt5
       and c_2 = 1/(2 sqrt(15/4 + sqrt5/2)) < 0.2267; |d/ds log wt1| <= M_1 = 14/15 (log wt1 increasing) and
       |d/ds log wt2| <= M_2 = d/ds log wt2 (sqrt5/2) < 0.5071 (log wt2 unimodal, maximum at s = 1/2 + sqrt2).  With the
       cone C_b = {f > 0 : f(x) <= e^(b|s(x)-s(y)|) f(y)}: Lt C_b c C_(max_i (M_i + c_i b)), so the fixed point
       b* = M_1/(1 - c_1) = 7(5+sqrt5)/30 = 1.688... is absorbing.  For a = 47/2 and n = 5, the cone-diameter lemma of
       CE.3 plus an exact oscillation recursion over the 32 level-5 branch intervals give
           diam_(C_a) Lt^5(C_a) <= Delta = 2.4401...,  theta = tanh(Delta/4) < 0.5442,  theta^(1/5) < 0.8855.
  CF.3 [proved from CF.1, CF.2] For every f in C_b, b < a, and every m >= 0, x in X:
           |log((Lt^(m+1) f)(x) / (Lt^m f)(x)) - log r| <= theta^floor(m/5) Theta_a(f, Lt f) < oo,
       where r is the spectral radius of Lt on C(X); r^-m Lt^m 1 converges uniformly to a positive eigenfunction
       ht in closure(C_b*) with Lt ht = r ht.  Consequently CE.3 holds with no literature assumption: for the level
       sums, s_(k+1)/s_k -> r with |log(s_(k+1)/(r s_k))| <= C_G theta^(k/5) (C_G depends on Lip log Lt^4 G, not
       explicit), and CE.4's enclosure r in [lo, hi] (width 2e-38) holds without an eigenfunction.
  CF.4 [certified consistency] With f = 1 (Theta_a(1, Lt 1) <= Theta_0 = 2 log((a+14/15)/(a-14/15)) + log(11/6)),
       arb values of Lt^m 1 at five points of X, m <= M_CHECK, obey CF.3's bound with r from CE.4.
  CF.5 [proved from CF.3] Every eigenvalue mu != r of Lt with a C^1 (complex) eigenfunction satisfies
           |mu| <= theta^(1/5) r < 0.8855 r < 0.3752.
       (Floats, no claim: Chebyshev collocation gives the next eigenvalue near -0.0731, ratio 0.1725.)
Assumptions used as stated (not re-run): CE.1, CE.2, CE.4 and the C^1 positivity of Lt^4 G from CE.3 (whose proof
there does not use Assumption CE).  No literature assumption is used.  The rate theta^(1/5) is far from the observed
0.1725 per level: an upper bound, not the value.  Stage 1 untouched; no W, 0.08, 0.23 or 30 microN/kW enters.
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

import sympy as sp
from flint import arb, ctx

FAILS: list[str] = []
A_CONE = Fr(47, 2)       # cone parameter a
N_STEP = 5               # contraction is proved for Lt^N_STEP
M_CHECK = 12             # CF.4 iterates in `reproduce`
M_FAST = 9               # CF.4 iterates in the pytest run
THETA_MAX = "0.5442"
RATE_MAX = "0.8855"


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- CF.1 Birkhoff's contraction, elementary proof
def birkhoff_lemma():
    """CF.1.  F = T f, G = T g; u = T(g - alpha f), v = T(beta f - g) lie in T(C u {0}); if v - m u, M u - v in C u {0}
    then beta(F,G) <= (beta + m alpha)/(1 + m) and alpha(F,G) >= (M alpha + beta)/(1 + M); then calculus."""
    al, be, m, M, F, G, z, y, tau = sp.symbols("alpha beta m M F G z y tau", positive=True)
    u = G - al * F
    v = be * F - G
    check(sp.expand((v - m * u) - ((be + m * al) * F - (1 + m) * G)) == 0,
          "CF.1: v - m u = (beta + m alpha) F - (1 + m) G, so beta(Tf,Tg) <= (beta + m alpha)/(1 + m)")
    check(sp.expand((M * u - v) - ((1 + M) * G - (M * al + be) * F)) == 0,
          "CF.1: M u - v = (1 + M) G - (M alpha + beta) F, so alpha(Tf,Tg) >= (M alpha + beta)/(1 + M)")
    # with x = beta/alpha = e^Theta: Theta(Tf,Tg) <= phi(Theta) - phi(0), phi(y) = log((e^y + m)/(e^y + M))
    phi = sp.log((sp.exp(y) + m) / (sp.exp(y) + M))
    dphi = sp.simplify(sp.diff(phi, y) - sp.exp(y) * (M - m) / ((sp.exp(y) + m) * (sp.exp(y) + M)))
    check(dphi == 0, "CF.1: phi'(y) = e^y (M - m)/((e^y + m)(e^y + M))")
    ident = sp.expand((z + m) * (z + M) - z * (sp.sqrt(m) + sp.sqrt(M)) ** 2 - (z - sp.sqrt(m * M)) ** 2)
    check(ident == 0, "CF.1: (z+m)(z+M) - z(sqrt m + sqrt M)^2 = (z - sqrt(mM))^2 >= 0, so phi' <= (sqrt M - sqrt m)/(sqrt M + sqrt m)")
    th = (sp.sqrt(m * sp.exp(4 * tau)) - sp.sqrt(m)) / (sp.sqrt(m * sp.exp(4 * tau)) + sp.sqrt(m))
    check(sp.simplify((th - sp.tanh(tau)).rewrite(sp.exp)) == 0,
          "CF.1: with M/m = e^(4 tau): (sqrt M - sqrt m)/(sqrt M + sqrt m) = tanh(tau) = tanh(log(M/m)/4) <= tanh(Delta/4)")
    # positive cones: Theta(f,g) >= log(sup(g/f)/inf(g/f)) because C u {0} consists of functions >= 0
    print("   CF.1: integrating phi' from 0 to Theta(f,g) gives Theta(Tf,Tg) <= tanh(Delta/4) Theta(f,g); degenerate cases "
          "u = 0 or v = 0 mean g = alpha f or g = beta f, Theta = 0. C_a u {0} is closed (f >= 0 log-Lipschitz, f(y) = 0 => f = 0).")


# ----------------------------------------------------------------- CF.2 the s-variable cone (sympy)
def symbolic_cone():
    s, t = sp.symbols("s t", positive=True)
    h = sp.Rational(1, 2)
    s0 = sp.sqrt(5) / 2
    sig1, sig2 = sp.sqrt(sp.Rational(15, 4) - s), sp.sqrt(sp.Rational(15, 4) + s)
    ph1, ph2 = -sp.Rational(5, 2) - s, -sp.Rational(5, 2) + s
    check(sp.simplify(ph1 + sp.Rational(25, 4) - sig1 ** 2) == 0 and sp.simplify(ph2 + sp.Rational(25, 4) - sig2 ** 2) == 0,
          "CF.2: in s-coordinates the branches are s -> sqrt(15/4 - s) (phi_1) and s -> sqrt(15/4 + s) (phi_2)")
    check(sp.simplify(sig1.subs(s, sp.Rational(5, 2)) - s0) == 0 and sp.simplify(sig2.subs(s, sp.Rational(5, 2)) - sp.Rational(5, 2)) == 0
          and sp.N(sig1.subs(s, s0)) < 2.5 and sp.N(sig2.subs(s, s0)) > sp.N(s0),
          "CF.2: both branches map [sqrt5/2, 5/2] into itself; sigma_1 >= sqrt5/2, sigma_2 >= sqrt(15/4 + sqrt5/2)")
    check(sp.simplify(sp.diff(sig1, s) + 1 / (2 * sig1)) == 0 and sp.simplify(sp.diff(sig2, s) - 1 / (2 * sig2)) == 0,
          "CF.2: |d sigma_i/ds| = 1/(2 sigma_i) <= c_1 = 1/sqrt5, c_2 = 1/(2 sqrt(15/4 + sqrt5/2))")
    W1 = (s + h) / (2 * s * (sp.Rational(7, 2) - s))
    W2 = (s - h) / (2 * s * (s + sp.Rational(7, 2)))
    d1 = sp.diff(sp.log(W1), s)
    d2 = sp.diff(sp.log(W2), s)
    check(sp.simplify(d1 - (1 / (s + h) - 1 / s + 1 / (sp.Rational(7, 2) - s))) == 0
          and sp.simplify(d2 - (1 / (s - h) - 1 / s - 1 / (s + sp.Rational(7, 2)))) == 0, "CF.2: d/ds log wt_i (sympy)")
    dd1 = sp.diff(d1, s)
    check(sp.simplify(dd1 - (1 / s ** 2 - 1 / (s + h) ** 2 + 1 / (sp.Rational(7, 2) - s) ** 2)) == 0 and d1.subs(s, sp.Rational(5, 2)) == sp.Rational(14, 15)
          and sp.N(d1.subs(s, s0), 30) > 0,
          "CF.2: (log wt1)'' = 1/s^2 - 1/(s+1/2)^2 + 1/(7/2-s)^2 > 0 and (log wt1)'(sqrt5/2) > 0, so log wt1 increases and M_1 = (log wt1)'(5/2) = 14/15")
    n2 = sp.factor(sp.numer(sp.together(d2)))
    check(sp.expand(n2 + 4 * s ** 2 - 4 * s - 7) == 0 and sp.solve(4 * s ** 2 - 4 * s - 7, s) == [h + sp.sqrt(2)],
          "CF.2: (log wt2)' = 0 iff 4s^2 - 4s - 7 = 0, i.e. s = 1/2 + sqrt2 (inside [sqrt5/2, 5/2]); log wt2 is unimodal")
    num, den = sp.fraction(sp.together(sp.diff(d2, s)))
    P = sp.Poly(num, s)
    check(P.count_roots(1, sp.Rational(5, 2)) == 0 and num.subs(s, 2) * den.subs(s, 2) < 0,
          "CF.2: (log wt2)'' has no zero on [1, 5/2] and is negative there, so |(log wt2)'| <= max((log wt2)'(sqrt5/2), 1/15) = (log wt2)'(sqrt5/2)")
    check(d2.subs(s, sp.Rational(5, 2)) == -sp.Rational(1, 15), "CF.2: (log wt2)'(5/2) = -1/15")
    bstar = sp.Rational(14, 15) / (1 - 1 / sp.sqrt(5))
    check(sp.simplify(bstar - 7 * (5 + sp.sqrt(5)) / 30) == 0 and sp.simplify(2 * bstar * (sp.Rational(5, 2) - s0) - sp.Rational(14, 3)) == 0,
          "CF.2: b* = M_1/(1 - c_1) = 7(5+sqrt5)/30; for comparison 2 b* diam_s(X) = 14/3 exactly")
    return W1, W2


# ----------------------------------------------------------------- CF.2 certified cone constants (arb)
def _lw1(s):
    return ((s + arb(1) / 2) / (2 * s * (arb(7) / 2 - s))).log()


def _lw2(s):
    return ((s - arb(1) / 2) / (2 * s * (s + arb(7) / 2))).log()


def _var1(lo, hi):
    return _lw1(hi) - _lw1(lo)          # log wt1 increasing


def _var2(lo, hi, crit):
    if bool(hi < crit) or bool(lo > crit):
        return abs(_lw2(hi) - _lw2(lo))  # monotone piece
    return abs(_lw2(crit) - _lw2(lo)) + abs(_lw2(hi) - _lw2(crit))  # unimodal (also a valid bound if crit is outside)


def _osc(lo, hi, n, a, crit):
    """Upper bound for osc_J log(Lt^n f), J = [lo, hi] in s, f in C_a:
    osc_J log Lt^n f <= max_i (Var_J log wt_i + osc_(sigma_i J) log Lt^(n-1) f), osc_J log f <= a |J|."""
    if n == 0:
        return a * (hi - lo)
    q = arb(15) / 4
    v1 = _var1(lo, hi) + _osc((q - hi).sqrt(), (q - lo).sqrt(), n - 1, a, crit)
    v2 = _var2(lo, hi, crit) + _osc((q + lo).sqrt(), (q + hi).sqrt(), n - 1, a, crit)
    return arb(max(v1.upper(), v2.upper()))


def cone_constants(a=A_CONE, n=N_STEP):
    ctx.prec = 200
    a = arb(a.numerator) / a.denominator
    s0, s1 = arb(5).sqrt() / 2, arb(5) / 2
    crit = arb(1) / 2 + arb(2).sqrt()
    M1 = arb(14) / 15
    M2 = 1 / (s0 - arb(1) / 2) - 1 / s0 - 1 / (s0 + arb(7) / 2)
    c1 = 1 / arb(5).sqrt()
    c2 = 1 / (2 * (arb(15) / 4 + s0).sqrt())
    check(bool(M2 < arb("0.5071")) and bool(c2 < arb("0.2267")) and bool(M2 > arb(1) / 15),
          f"CF.2: M_2 = {M2.str(10)} < 0.5071, c_2 = {c2.str(10)} < 0.2267 (c_1 = 1/sqrt5, M_1 = 14/15)")
    b = a
    bs = []
    for _ in range(n):
        b = arb(max((M1 + c1 * b).upper(), (M2 + c2 * b).upper()))
        bs.append(b)
    check(bool(bs[0] < a), f"CF.2: Lt C_a c C_b1, b1 = {bs[0].str(8)} < a = {a.str(5)} (so Lt maps C_a into itself)")
    V = _osc(s0, s1, n, a, crit)
    Delta = 2 * ((a + b) / (a - b)).log() + 2 * V
    theta = (Delta / 4).tanh()
    rate = (theta.log() / n).exp()
    print(f"   CF.2: b_{n} = {b.str(8)} (b* = {(M1 / (1 - c1)).str(8)}), osc_X log Lt^{n} f <= {V.str(8)} for f in C_a")
    check(bool(theta < arb(THETA_MAX)) and bool(rate < arb(RATE_MAX)),
          f"CF.2: diam_(C_a) Lt^{n} C_a <= Delta = 2 log((a+b_{n})/(a-b_{n})) + 2 osc = {Delta.str(8)}, "
          f"theta = tanh(Delta/4) = {theta.str(8)} < {THETA_MAX}, per level theta^(1/{n}) = {rate.str(8)} < {RATE_MAX}")
    return dict(a=a, b1=bs[0], theta=theta, rate=rate, Delta=Delta)


# ----------------------------------------------------------------- CF.4 certified consistency with f = 1
def _Lt_iter(sv, m):
    """(Lt^m 1)(t(s)) by the full 2^m branch tree (arb)."""
    if m == 0:
        return arb(1)
    q = arb(15) / 4
    w1 = (sv + arb(1) / 2) / (2 * sv * (arb(7) / 2 - sv))
    w2 = (sv - arb(1) / 2) / (2 * sv * (sv + arb(7) / 2))
    return w1 * _Lt_iter((q - sv).sqrt(), m - 1) + w2 * _Lt_iter((q + sv).sqrt(), m - 1)


def pointwise_consistency(K, consts, r_lo, r_hi):
    ctx.prec = 200
    a, b1, theta = consts["a"], consts["b1"], consts["theta"]
    B = arb(14) / 15                     # Lt 1 = 4/(6-t) lies in C_(M_1) (and 1 in C_0)
    Theta0 = 2 * ((a + B) / (a - B)).log() + (arb(11) / 6).log()
    lr_lo, lr_hi = arb(r_lo.numerator).__truediv__(r_lo.denominator).log(), (arb(r_hi.numerator) / r_hi.denominator).log()
    pts = [Fr(-5), Fr(-3), Fr(-5, 2), Fr(-1), Fr(0)]
    worst_ratio = 0.0
    ok = True
    for m in range(K):
        bound = theta ** (m // N_STEP) * Theta0
        devs = []
        for t in pts:
            sv = (arb(t.numerator) / t.denominator + arb(25) / 4).sqrt()
            lR = (_Lt_iter(sv, m + 1) / _Lt_iter(sv, m)).log()
            dev = arb(max(abs(lR - lr_lo).upper(), abs(lR - lr_hi).upper()))
            devs.append(dev)
            ok = ok and bool(dev < bound)
        worst = max(float(d) for d in devs)
        worst_ratio = max(worst_ratio, worst / float(bound))
        if m in (0, 4, 5, K - 1):
            print(f"   m = {m:2d}: max_x |log R_m(x) - log r| = {worst:.3e} <= bound {float(bound):.3e}")
    check(ok, f"CF.4: Theta_0 = {Theta0.str(6)}; for m < {K} at t in {{-5,-3,-5/2,-1,0}}: |log(Lt^(m+1)1/Lt^m 1) - log r| "
              f"<= theta^floor(m/5) Theta_0 (max used fraction {worst_ratio:.3f})")
    return ok


# ----------------------------------------------------------------- floats only: next eigenvalue (no claim)
def collocation(nodes=60):
    import numpy as np
    k = np.arange(nodes)
    x = np.cos(np.pi * (k + 0.5) / nodes)        # Chebyshev nodes on [-1,1], t = 5(x-1)/2

    def cheb_vander(xx):
        return np.polynomial.chebyshev.chebvander(xx, nodes - 1)

    t = 2.5 * (x - 1)
    s = np.sqrt(t + 6.25)
    w1 = (s + .5) / (2 * s * (3.5 - s))
    w2 = (s - .5) / (2 * s * (s + 3.5))
    p1, p2 = -2.5 - s, -2.5 + s
    V = cheb_vander(x)
    Vi = np.linalg.inv(V)
    A = (w1[:, None] * cheb_vander(p1 / 2.5 + 1) + w2[:, None] * cheb_vander(p2 / 2.5 + 1)) @ Vi
    ev = sorted(np.linalg.eigvals(A), key=lambda z: -abs(z))
    return ev[0].real, ev[1]


def eigen_gap(consts, r_hi):
    r0, l2 = collocation()
    rate = float(consts["rate"])
    print(f"   floats (no claim): r ~ {r0:.15f}, next eigenvalue ~ {complex(l2).real:.12f}{complex(l2).imag:+.1e}i, ratio {abs(l2) / r0:.4f}")
    check(abs(l2) / r0 <= rate and abs(r0 - 0.42373128696260794) < 1e-12,
          f"CF.5: bound |mu| <= theta^(1/5) r < {RATE_MAX} * {float(r_hi):.6f} = {rate * float(r_hi):.4f}; collocation ratio "
          f"{abs(l2) / r0:.4f} is consistent with it (the bound is not sharp)")


def main(K=M_CHECK) -> int:
    t0 = time.time()
    print("Theorem CF: Birkhoff's contraction (Assumption CE) proved for CE's transfer operator, with explicit rate "
          "theta^(1/5) < 0.8855 per level. Not thrust.")
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import gasket_fractional_profile_level_ratio as ce
    birkhoff_lemma()
    symbolic_cone()
    consts = cone_constants()
    r_lo, r_hi = ce.certificate()
    pointwise_consistency(K, consts, r_lo, r_hi)
    eigen_gap(consts, r_hi)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CF checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else M_CHECK))
