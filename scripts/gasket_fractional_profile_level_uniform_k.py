#!/usr/bin/env python3
"""Theorem CW: a bound on rho(L_k) that holds for every k (CV's open item "nothing uniform in k"), which proves CV's
observation 3^k rho(L_k) -> 1/3.

Setting exactly as in Theorems CT-CV (GASKET.md): s = sqrt(t + 25/4), X = [sqrt5/2, 5/2], sigma_(1,2)(s) = sqrt(15/4 -+ s),
    (L_k f)(s) = u_1(s) f(sigma_1 s) + u_2(s) f(sigma_2 s),   u_i = wt_i / (2 sigma_i)^k,
    wt_1 = (s+1/2)/(2s(7/2-s)),  wt_2 = (s-1/2)/(2s(s+7/2)),
rho(L_k) = spectral radius on C(X) (= ||L_k^n 1||^(1/n) limit; k = 0 is CE's twisted operator, rho(L_0) = r).
s* = 3/2 is the fixed point of sigma_1, 2 sigma_1(s*) = 3, wt_1(s*) = 1/3 (exact).

Construction (exact algebra).  With x_0 = s, x_(j+1) = sigma_1(x_j),
    Phi(s) = - sum_(j>=0) log(2 sigma_1(x_j)/3),     Psi(s) = sum_(j>=0) log(3 wt_1(x_j)),
both converge uniformly on X (CW.1), and telescoping gives Phi(sigma_1 s) - Phi(s) = log(2 sigma_1(s)/3),
Psi(sigma_1 s) - Psi(s) = -log(3 wt_1(s)).  For h_k = exp(k Phi + Psi) > 0 this gives EXACTLY
    (L_k h_k)/h_k = 3^(-k-1) (1 + G(s) B(s)^k),
    B = (3/(2 sigma_2)) exp(Phi(sigma_2 s) - Phi(s)),   G = 3 wt_2 exp(Psi(sigma_2 s) - Psi(s)),
so the branch sigma_1 contributes exactly 3^(-k-1) at every point.  Collatz-Wielandt for the positive operator L_k on
C(X) (m h <= L h <= M h  =>  m <= rho <= M) then gives for every real k >= 0
    3^(-k-1) (1 + min G (min B)^k) <= rho(L_k) <= 3^(-k-1) (1 + max G (max B)^k).

Theorem CW (checked here; exit 1 on any failure).
  CW.1 [proved; exact] sigma_1, sigma_2 map X' = [lo, 5/2] (lo = dyadic <= sqrt5/2) into X, and sigma_1 is a contraction
       towards s* with |sigma_1 x - sigma_1 y| <= |x - y|/sqrt5 on X'; |log(2 sigma_1(x)/3)| <= (2/5)|x - s*| and
       |log(3 wt_1(x))| <= 3 |x - s*| on X'.  So the j >= J tails of Phi, Psi are at most (2/5, 3) 5^(-J/2)/(1 - 5^(-1/2)).
  CW.2 [certified, arb, adaptive cover of X'] 0.5999 < B < 0.8352 and 0.1999 < G < 0.2873 on X.
  CW.3 [proved from CW.1, CW.2] For EVERY real k >= 0:
           0.1999 * 0.5999^k  <=  3^(k+1) rho(L_k) - 1  <=  0.2873 * 0.8352^k.
       Hence 3^k rho(L_k) > 1/3 for every k (CV.7's lower bound, now for all k), 3^k rho(L_k) -> 1/3 geometrically,
       rho(L_(k+1))/rho(L_k) -> 1/3, and (1/k) log rho(L_k) -> -log 3.
  CW.4 [consistency, finite] CE's r (k = 0) and CT.3/CU.3/CV.3's certified rho(L_k), k = 1..12 (used as stated), satisfy
       CW.3's two-sided bound.
Observation (floats, no claim): min B = 3/5 and min G = 1/5 are attained at s = 5/2 (the fixed point of sigma_2), and
max B ~ 0.835144 at s = sqrt5/2; max G ~ 0.287279 in the interior of X.  CV's data give (3^(k+1) rho - 1) ratios ~ 0.78 at k = 11..12, between
the two bases; the exact geometric rate is not determined.
Kept limits: the rate of 3^(k+1) rho(L_k) - 1 is only bracketed in [0.5999, 0.8352]; nothing about nu_2^(k), so the
heuristic q^(k) -> -1/3 stays OPEN (HEURISTIC); monotonicity of 3^k rho(L_k) beyond k = 12 is not implied.
Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time
from fractions import Fraction as Fr

from flint import arb, ctx

FAILS: list[str] = []
RESULT: dict = {}
PREC = 160
J = 120                       # orbit terms kept; tail <= 3 * 5^(-60) / (1 - 5^(-1/2)) < 1e-41
LO = Fr(1172343, 2 ** 20)     # dyadic, LO <= sqrt5/2 (checked exactly below)
HI = Fr(5, 2)
BOUNDS = dict(B_lo=Fr(5999, 10000), B_hi=Fr(8352, 10000), G_lo=Fr(1999, 10000), G_hi=Fr(2873, 10000))
MAXDEPTH = 40
# certified values used as stated (consistency only, CW.4): k: (mid, radius)
RHO_PRIOR = {0: ("0.42373128696260793985253216132843923", "2e-38"),
             1: ("0.1343549829703998131", "2.9e-20"), 2: ("0.0430703112227101271", "1.3e-20"),
             3: ("0.0139192547691468646", "2.1e-20"), 4: ("0.0045264690852562123", "2.8e-20"),
             5: ("0.0014792559239634748", "6e-20"), 6: ("0.0004853393939726874", "6e-20"),
             7: ("0.0001597483711719488", "6e-20"), 8: ("0.0000527169494750854", "6e-20"),
             9: ("0.0000174329038013139", "6e-20"), 10: ("0.0000057745266397710", "6e-20"),
             11: ("0.0000019153320627478", "6e-20"), 12: ("0.0000006359658135257", "6e-20")}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _a(x: Fr) -> arb:
    return arb(x.numerator) / arb(x.denominator)


def ball(a: Fr, b: Fr) -> arb:
    """An arb enclosing [a, b] (dyadic endpoints: exact)."""
    return (_a(a) + _a(b)) / 2 + (_a(b) - _a(a)) / 2 * arb(0, 1)


def sig1(x):
    return (arb(15) / 4 - x).sqrt()


def sig2(x):
    return (arb(15) / 4 + x).sqrt()


def wt1(x):
    return (x + arb(1) / 2) / (2 * x * (arb(7) / 2 - x))


def wt2(x):
    return (x - arb(1) / 2) / (2 * x * (x + arb(7) / 2))


TAIL_PHI = TAIL_PSI = None


def phipsi(x):
    """Enclosures of Phi(x), Psi(x) for every point of the ball x (subset of X'): J terms plus CW.1's tail."""
    phi = arb(0)
    psi = arb(0)
    for _ in range(J):
        y = sig1(x)
        phi -= (2 * y / 3).log()
        psi += (3 * wt1(x)).log()
        x = y
    return phi + arb(0, 1) * TAIL_PHI, psi + arb(0, 1) * TAIL_PSI


def logBG(I):
    p0, q0 = phipsi(I)
    y = sig2(I)
    p1, q1 = phipsi(y)
    lb = (arb(3) / (2 * y)).log() + p1 - p0
    lg = (3 * wt2(I)).log() + q1 - q0
    return lb, lg


def cw1() -> bool:
    global TAIL_PHI, TAIL_PSI
    ok = True
    # LO <= sqrt5/2  <=>  4 LO^2 <= 5 (exact); LO > 1/2 so wt_2 > 0; |s - s*| <= 1 on X'
    ok &= 4 * LO * LO <= 5 and LO > Fr(1, 2) and Fr(3, 2) - LO <= 1 and HI - Fr(3, 2) <= 1
    # sigma_1(X') subset [sqrt5/2, 5/2] = X: 15/4 - s in [5/4, 15/4 - LO], 15/4 - LO <= 25/4; sigma_2(X') subset X:
    # 15/4 + s in [15/4 + LO, 25/4], 15/4 + LO >= 5/4
    ok &= Fr(15, 4) - LO <= Fr(25, 4) and Fr(15, 4) + LO >= Fr(5, 4)
    # contraction: sigma_1 x - sigma_1 y = (y - x)/(sigma_1 x + sigma_1 y), sigma_1 >= sqrt5/2 on X'  =>  Lip <= 1/sqrt5.
    # |log(2 sigma_1(x)/3)| = |log sigma_1 x - log sigma_1 s*| <= |sigma_1 x - sigma_1 s*| / min sigma_1
    #                       <= |x - s*| / (sqrt5 * sqrt5/2) = (2/5)|x - s*|
    # |log(3 wt_1 x)| = |log wt_1 x - log wt_1 s*| <= sup|(log wt_1)'| |x - s*|, (log wt_1)' = 1/(s+1/2) - 1/s + 1/(7/2-s):
    # on X' each term is at most 1/(LO+1/2), 1/LO, 1/(7/2-5/2) = 1, and the first two have opposite signs.
    lw = max(Fr(1) / (LO + Fr(1, 2)), Fr(1) / LO) + 1
    ok &= lw <= 3
    # x_0 in X' (|x_0 - s*| <= 1), x_j in X for j >= 1, |x_j - s*| <= 5^(-j/2): tail_j>=J <= c 5^(-J/2)/(1 - 5^(-1/2))
    geo = arb(5) ** (-arb(J) / 2) / (1 - 1 / arb(5).sqrt())
    TAIL_PHI, TAIL_PSI = arb(2) / 5 * geo, 3 * geo     # arb(0, 1) * T encloses [-sup T, sup T]
    ok &= bool(TAIL_PSI < arb("1e-40"))
    # exact identity check at s*: Phi(s*) = Psi(s*) = 0 (every term vanishes), and telescoping at a sample point
    p, q = phipsi(arb(3) / 2)
    ok &= bool(abs(p) < arb("1e-39")) and bool(abs(q) < arb("1e-39"))
    s = arb(2)
    p0, q0 = phipsi(s)
    p1, q1 = phipsi(sig1(s))
    ok &= bool(abs(p1 - p0 - (2 * sig1(s) / 3).log()) < arb("1e-38"))
    ok &= bool(abs(q1 - q0 + (3 * wt1(s)).log()) < arb("1e-38"))
    check(ok, f"CW.1: X' = [{float(LO):.9f}, 5/2] contains X, sigma_1(X'), sigma_2(X') subset X, Lip(sigma_1) <= 1/sqrt5, "
              f"|log(2 sigma_1/3)| <= (2/5)|x - s*|, |log 3 wt_1| <= {float(lw):.4f}|x - s*| <= 3|x - s*|; J = {J} terms, "
              f"tails <= {TAIL_PHI.str(3, radius=False)}, {TAIL_PSI.str(3, radius=False)}; Phi(s*) = Psi(s*) = 0 and the "
              "telescoping identities hold at s = 2")
    return ok


def cw2() -> bool:
    lbB, ubB = _a(BOUNDS["B_lo"]).log(), _a(BOUNDS["B_hi"]).log()
    lbG, ubG = _a(BOUNDS["G_lo"]).log(), _a(BOUNDS["G_hi"]).log()
    stack = [(LO, HI, 0)]
    n_ok = 0
    worst = dict(Bmax=None, Bmin=None, Gmax=None, Gmin=None)
    fail = None
    while stack:
        a, b, d = stack.pop()
        lb, lg = logBG(ball(a, b))
        good = bool(lb < ubB) and bool(lb > lbB) and bool(lg < ubG) and bool(lg > lbG)
        if good:
            n_ok += 1
            for key, val, up in (("Bmax", lb, True), ("Bmin", lb, False), ("Gmax", lg, True), ("Gmin", lg, False)):
                e = val.exp()
                v = arb(e.mid() + e.rad()) if up else arb(e.mid() - e.rad())
                if worst[key] is None or (up and bool(v > worst[key])) or (not up and bool(v < worst[key])):
                    worst[key] = v
            continue
        if d >= MAXDEPTH:
            fail = (a, b)
            break
        m = (a + b) / 2
        stack += [(a, m, d + 1), (m, b, d + 1)]
    ok = fail is None
    RESULT.update(cover=n_ok, **worst)
    check(ok, f"CW.2: adaptive arb cover of X' by {n_ok} intervals: {BOUNDS['B_lo']} < B < {BOUNDS['B_hi']}, "
              f"{BOUNDS['G_lo']} < G < {BOUNDS['G_hi']} on X (enclosure extremes: B in [{worst['Bmin'].str(7, radius=False) if worst['Bmin'] else '?'}, "
              f"{worst['Bmax'].str(7, radius=False) if worst['Bmax'] else '?'}], G in [{worst['Gmin'].str(7, radius=False) if worst['Gmin'] else '?'}, "
              f"{worst['Gmax'].str(7, radius=False) if worst['Gmax'] else '?'}])" + ("" if ok else f"; failed on [{float(fail[0])}, {float(fail[1])}]"))
    return ok


def cw3_cw4() -> None:
    bl, bh, gl, gh = (_a(BOUNDS[x]) for x in ("B_lo", "B_hi", "G_lo", "G_hi"))
    rows, ok = [], True
    for k, (m, r) in RHO_PRIOR.items():
        rho = arb(m) + arb(0, r)
        a = rho * arb(3) ** (k + 1) - 1
        lo, hi = gl * bl ** k, gh * bh ** k
        ok &= bool(a > lo) and bool(a < hi)
        rows.append(f"k={k}: {lo.str(4, radius=False)} < {a.str(6, radius=False)} < {hi.str(4, radius=False)}")
    RESULT["consistency"] = ok
    check(ok, "CW.4 (consistency): CE's r and CT/CU/CV's certified rho(L_k) satisfy 0.1999*0.5999^k < 3^(k+1) rho(L_k) - 1 "
              "< 0.2873*0.8352^k:\n      " + "\n      ".join(rows))
    # CW.3 is the Collatz-Wielandt implication of CW.2 (no computation beyond CW.2); print sample consequences
    for k in (13, 20, 50):
        lo = (1 + gl * bl ** k) / arb(3) ** (k + 1)
        hi = (1 + gh * bh ** k) / arb(3) ** (k + 1)
        print(f"   CW.3 at k = {k}: rho(L_{k}) in [{lo.str(10, radius=False)}, {hi.str(10, radius=False)}]; "
              f"3^k rho - 1/3 <= {(gh * bh ** k / 3).str(4, radius=False)}", flush=True)
    ratio = [(arb(RHO_PRIOR[k + 1][0]) * 3 ** (k + 2) - 1) / (arb(RHO_PRIOR[k][0]) * 3 ** (k + 1) - 1) for k in (10, 11)]
    print("   observation (floats, no claim): (3^(k+2) rho_(k+1) - 1)/(3^(k+1) rho_k - 1) = "
          + ", ".join(x.str(5, radius=False) for x in ratio) + " at k = 10, 11 (between 0.5999 and 0.8352)", flush=True)


def float_observation() -> None:
    old = ctx.prec
    ctx.prec = 64
    try:
        pts = {"5/2": arb(5) / 2, "sqrt5/2": arb(5).sqrt() / 2}
        vals = {k: tuple(v.exp().str(7, radius=False) for v in logBG(p)) for k, p in pts.items()}
        grid = [Fr(1118034, 10 ** 6) + Fr(1381966, 10 ** 6) * i / 400 for i in range(401)]
        gv = [(float(logBG(_a(x))[1].exp().mid()), float(x)) for x in grid]
        gmax = max(gv)
        print(f"   observation (no claim): (B, G) at 5/2 = {vals['5/2']} (exactly 3/5, 1/5: fixed point of sigma_2), "
              f"at sqrt5/2 = {vals['sqrt5/2']} (max B); grid max G = {gmax[0]:.6f} at s = {gmax[1]:.4f}", flush=True)
    finally:
        ctx.prec = old


def main() -> int:
    print("Theorem CW: 0.1999*0.5999^k <= 3^(k+1) rho(L_k) - 1 <= 0.2873*0.8352^k for every k >= 0; "
          "3^k rho(L_k) -> 1/3. Not thrust.", flush=True)
    t0 = time.time()
    FAILS.clear()
    RESULT.clear()
    old = ctx.prec
    ctx.prec = PREC
    try:
        ok = cw1()
        ok = cw2() and ok
        check(ok, "CW.3: (L_k h_k)/h_k = 3^(-k-1)(1 + G B^k) exactly with h_k = exp(k Phi + Psi) > 0 continuous, so by "
                  "Collatz-Wielandt on C(X), for every real k >= 0: 0.1999*0.5999^k <= 3^(k+1) rho(L_k) - 1 <= "
                  "0.2873*0.8352^k; hence 3^k rho(L_k) > 1/3 for all k and 3^k rho(L_k) -> 1/3")
        cw3_cw4()
        float_observation()
    finally:
        ctx.prec = old
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CW checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
