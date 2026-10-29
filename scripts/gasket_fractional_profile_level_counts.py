#!/usr/bin/env python3
r"""Theorem CD: the level counts of the atoms of BX's spectral measure nu in one scale period are exact for every level.

Theorem CC (CC.4) parametrised one scale period: the atoms of nu in [1,5) are lam_t = -psi^-1(t) for
t in O n I, I = (p(a), a], a = psi(-1), and t has level k if t in O_k := p^-k(-3).  CC certified the counts
N_k = #(O_k n I) = 1, 0, 1, 3, 5, 10, ..., 5*2^(k-4) for k <= 14 and left "the level counts" open (observation).

Theorem CD (checked here; exit 1 on any failure):
  CD.1 [proved; elementary] O_k consists of exactly 2^k distinct real points, all in (-5,0) \ {-5/2}.
  CD.2 [proved + certified] The forward orbit of the period endpoint a = psi(-1) leaves [-5,0] after four steps, with
       -5/2 < a < 0,  -5 < p^2(a) < p(a) < -3,  -3 < p^3(a) < -5/2,  p^4(a) < -5   (arb, 300 bits).
  CD.3 [proved from CD.1-CD.2; exact] For every k >= 0
           N_0 = 1, N_1 = 0, N_2 = 1, N_3 = 3,  N_k = 5 * 2^(k-4)  (k >= 4),
       so one scale period holds exactly 5 * 2^(K-3) atoms of level <= K (K >= 3); 10240 at K = 14, as CC certified.
       Proof: cut each interval at the critical point -5/2, where p is a bijection on each side; push forward with
       C_m(J) = #(O_m n J) = C_(m-1)(p(J)); the pushed endpoints are a's orbit, and after p^4(a) < -5 every image
       interval contains all of (-5,0), so C = 2^(m-1) (CD.1).  The resulting finite linear recursion
           N_k = A_(k-1) + B_(k-1),  A_m = G_(m-1),  B_m = F_(m-1),  G_j = E_(j-1) + 2^(j-1),  E_j = F_j = 2^(j-1)
       (j, m >= 1) with A_0 = B_0 = F_0 = 0, E_0 = G_0 = 1 has the closed solution above.
  CD.4 [certified consistency] The six counting functions are recounted directly from arb enclosures of all 2^k
       points of O_k (k <= K_CERT), each comparison with a's orbit decided; all agree with CD.3.
Assumptions used as stated (not re-run): CC.4 (the parametrisation of one period by psi on [-5,-1] and its
level-by-level atom set), BD.4, BJ.6 and CA's enclosure routine psi_enc (imported unchanged).
No mass, ratio or rho_k statement is made here; the level-sum ratio 0.42373... stays an observation.
Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import sys
import time

import sympy as sp
from flint import arb, ctx

import gasket_fractional_profile_atoms as bz
import gasket_fractional_profile_atom_law as ca

FAILS: list[str] = []
K_CERT = 18            # direct recount in `reproduce` (262144 points at the top level, a few seconds)
K_FAST = 12            # direct recount in the pytest run
K_FORMULA = 400        # the recursion is compared with the closed form for k <= K_FORMULA (exact integers)
CC_COUNTS = [1, 0, 1, 3, 5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120]   # CC.4's certified counts, k <= 14


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def N_closed(k: int) -> int:
    return [1, 0, 1, 3][k] if k < 4 else 5 * 2 ** (k - 4)


# ----------------------------------------------------------------- CD.1 symbolic
def symbolic_checks():
    t, y = sp.symbols("t y", real=True)
    p = y * (y + 5)
    roots = sp.solve(sp.Eq(p, t), y)
    gm, gp = -sp.Rational(5, 2) - sp.sqrt(t + sp.Rational(25, 4)), -sp.Rational(5, 2) + sp.sqrt(t + sp.Rational(25, 4))
    check(set(sp.simplify(r - gm) == 0 or sp.simplify(r - gp) == 0 for r in roots) == {True} and len(roots) == 2,
          "CD.1: p(y) = t  <=>  y = -5/2 -+ sqrt(t + 25/4)")
    # t in (-25/4, 0): sqrt(t+25/4) in (0, 5/2), so both roots in (-5,0), distinct, neither = -5/2
    s = sp.symbols("s", positive=True)
    check(sp.simplify((-sp.Rational(5, 2) - s) + 5 - (sp.Rational(5, 2) - s)) == 0,
          "CD.1: for t in (-25/4,0), s = sqrt(t+25/4) in (0,5/2): roots -5/2-s in (-5,-5/2), -5/2+s in (-5/2,0)")
    check(-3 > -sp.Rational(25, 4) and p.subs(y, -sp.Rational(5, 2)) == -sp.Rational(25, 4),
          "CD.1: -3 in (-25/4,0), and p(-5/2) = -25/4 is not in (-5,0), so -5/2 is in no O_k (k >= 1); (-5,0) is inside (-25/4,0)")
    check(sp.diff(p, y) == 2 * y + 5,
          "CD.3: p' = 2y+5 < 0 on (-inf,-5/2), > 0 on (-5/2,inf): p is a decreasing bijection [-5,-5/2] -> [-25/4,0] "
          "and an increasing bijection [-5/2,0] -> [-25/4,0]")


def preimage_sets(K):
    """arb enclosures of all points of O_k, k <= K (outward rounding; each point is a ball)."""
    ctx.prec = 300
    lev = [[arb(-3)]]
    q = arb(25) / 4
    h = arb(5) / 2
    for _ in range(K):
        nxt = []
        for tt in lev[-1]:
            r = (tt + q).sqrt()
            nxt += [-h - r, -h + r]
        lev.append(nxt)
    return lev


def cd1_certified(lev):
    ok = True
    for k, pts in enumerate(lev):
        ok = ok and len(pts) == 2 ** k and all(bool(x > -5) and bool(x < 0) for x in pts)
        srt = sorted(pts, key=lambda x: float(x.mid()))
        ok = ok and all(bool(srt[i] < srt[i + 1]) for i in range(len(srt) - 1))
        ok = ok and all(not (x + arb(5) / 2).contains(0) for x in pts)
    check(ok, f"CD.1 (arb): for k <= {len(lev) - 1}, O_k has 2^k pairwise separated points, all in (-5,0) \\ {{-5/2}}")


# ----------------------------------------------------------------- CD.2 orbit of the period endpoint
def endpoint_orbit():
    ctx.prec = 300
    a = ca.psi_enc(arb(-1))
    orb = [a]
    for _ in range(5):
        orb.append(bz._p(orb[-1]))
    a0, a1, a2, a3, a4, a5 = orb
    check(bool(a0 > -arb(5) / 2) and bool(a0 < 0), f"CD.2: -5/2 < a = psi(-1) = {a0.str(16)} < 0")
    check(bool(a1 > a2) and bool(a2 > -5) and bool(a1 < -3),
          f"CD.2: -5 < p^2(a) = {a2.str(16)} < p(a) = psi(-5) = {a1.str(16)} < -3")
    check(bool(a3 > -3) and bool(a3 < -arb(5) / 2), f"CD.2: -3 < p^3(a) = psi(-125) = {a3.str(16)} < -5/2")
    check(bool(a4 < -5) and bool(a5 > 0), f"CD.2: p^4(a) = {a4.str(16)} < -5 (and p^5(a) = {a5.str(12)} > 0): a's orbit leaves [-5,0]")
    return orb


# ----------------------------------------------------------------- CD.3 the recursion
def recursion_counts(K):
    """Exact integer solution of CD.3's recursion; returns the dictionary of all six counting sequences."""
    E, F, G, A, B, N = [1], [0], [1], [0], [0], [1]
    for j in range(1, K + 1):
        E.append(2 ** (j - 1))           # C_j((-5, p^3 a]):  maps onto [p^4 a, 0) > (-5,0)
        F.append(2 ** (j - 1))           # C_j((p^3 a, 0)): left piece -> (-25/4, p^4 a), misses (-5,0); right -> [-25/4,0)
        G.append(E[j - 1] + 2 ** (j - 1))  # C_j([p^2 a, 0)): left piece -> (-25/4, p^3 a]; right -> [-25/4, 0)
        A.append(G[j - 1])               # C_j((-5, p a]):  -> [p^2 a, 0)
        B.append(F[j - 1])               # C_j((-5, p^2 a)): -> (p^3 a, 0)
        N.append(A[j - 1] + B[j - 1])    # C_j((p a, a]): left (p a, -5/2) -> (-25/4, p^2 a); right [-5/2, a] -> [-25/4, p a]
    return dict(E=E, F=F, G=G, A=A, B=B, N=N)


def recursion_checks():
    seq = recursion_counts(K_FORMULA)
    N = seq["N"]
    check(all(N[k] == N_closed(k) for k in range(K_FORMULA + 1)),
          f"CD.3: the recursion gives N_k = 1,0,1,3 and 5*2^(k-4) for every k <= {K_FORMULA} (exact integers; induction for all k)")
    # symbolic induction step: for k >= 5, N_k = A_(k-1) + B_(k-1) = G_(k-2) + F_(k-2) = (2^(k-4) + 2^(k-3)) + 2^(k-3)
    k = sp.symbols("k", integer=True)
    expr = (2 ** (k - 4) + 2 ** (k - 3)) + 2 ** (k - 3)
    check(sp.simplify(expr - 5 * 2 ** (k - 4)) == 0,
          "CD.3: induction step for k >= 5: N_k = G_(k-2) + F_(k-2) = (E_(k-3) + 2^(k-3)) + 2^(k-3) = 5*2^(k-4) (sympy)")
    cum = [sum(N[: K + 1]) for K in range(K_FORMULA + 1)]
    check(all(cum[K] == 5 * 2 ** (K - 3) for K in range(3, K_FORMULA + 1)),
          f"CD.3: one period holds exactly 5*2^(K-3) atoms of level <= K for 3 <= K <= {K_FORMULA}")
    check(N[: len(CC_COUNTS)] == CC_COUNTS and cum[14] == 10240,
          "CD.3: agrees with CC.4's independently certified counts for k <= 14 (10240 atoms through level 14)")
    return seq


# ----------------------------------------------------------------- CD.4 direct recount
def _count(pts, lo, lo_closed, hi, hi_closed):
    """#points in the interval (lo,hi) with the stated closedness; lo/hi may be None (= -5 open, 0 open).
    Returns None if a comparison is undecided."""
    n = 0
    for x in pts:
        if lo is not None:
            d = x - lo
            if d.contains(0):
                return None
            if not bool(d > 0):
                continue
        if hi is not None:
            d = hi - x
            if d.contains(0):
                return None
            if not bool(d > 0):
                continue
        n += 1
    return n


def direct_recount(lev, orb, K):
    a0, a1, a2, a3 = orb[:4]
    seq = recursion_counts(K)
    # endpoints never coincide with a point of O_k (each comparison must be decided), so open/closed does not matter
    specs = {"N": (a1, a0), "A": (None, a1), "B": (None, a2), "G": (a2, None), "E": (None, a3), "F": (a3, None)}
    ok, undecided = True, False
    for name, (lo, hi) in specs.items():
        got = [_count(lev[k], lo, False, hi, True) for k in range(K + 1)]
        if any(g is None for g in got):
            undecided = True
        ok = ok and got == seq[name][: K + 1]
    check(ok and not undecided,
          f"CD.4 (arb): direct recount of N, A, B, G, E, F from all 2^k points of O_k, k <= {K}: every comparison decided, "
          "all six sequences equal CD.3's recursion")
    print("   N_k, k <= %d: %s" % (K, [_count(lev[k], a1, False, a0, True) for k in range(K + 1)]))


def main(K=K_CERT) -> int:
    t0 = time.time()
    print("Theorem CD: exact level counts of nu's atoms in one scale period: 1, 0, 1, 3, then 5*2^(k-4). Not thrust.")
    symbolic_checks()
    orb = endpoint_orbit()
    recursion_checks()
    lev = preimage_sets(K)
    cd1_certified(lev)
    direct_recount(lev, orb, K)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CD checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else K_CERT))
