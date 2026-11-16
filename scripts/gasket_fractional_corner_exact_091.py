#!/usr/bin/env python3
"""Theorem EM: the exact sign-change level at alpha = 0.91 and at alpha = 0.89 (N = 1).  EL kept "at 0.91 the last negative level is
still only bracketed, [2564, 2586] (float braces -+1e-5 at 2585/2586), and the lower side (EG's routine, whose cell bound mixes ends
the same way) is unchanged", "0.88, 0.89 brackets unchanged", and left open "the mirror of EL.1 for EG's lower integrals
(inf_[a,b] R_k >= R^lo(b) at one point)".

  EM.1 (THEOREM; every alpha in (gamma, 1), k0 >= 1, N >= 1; EG.1 and AQ (3) as stated).  On every cell [a, b] (0 < a < b) and every
       k >= k0:  inf_[a,b] R_k >= 1/V_k0(b)  and  inf_[a,b] Rh_k >= (5/3)/V_(k0+1)(5b),  with EG.1's V_k0 evaluated at the single
       right node b.  Proof: R_k(x) = (3/5)^k r_k(5^-k x) is nonincreasing in x (AQ (3)), so inf_[a,b] R_k = R_k(b) >= 1/V_k0(b)
       (EG.1 is pointwise: 1/R_k <= V_k0); Rh_k(x) = (5/3) R_(k+1)(5x) likewise.  Every other factor of EG.1's lower integrands is
       already taken at b by EG (2/(6 + 5^(-k0-1) x), A, (Phi+6)/(Phi+5), eta_N(Phi(5^-N x)) are nonincreasing in x).  So EG.1's cell
       sums for c_lo(k0), kappa_lo(N, k0) stay lower bounds with R, Rh at the right node.  (EG's vcell takes A_k0, e_m at a and Phi,
       P_m at b: the two large terms come from opposite ends.  OBSERVED: at 0.905, k0 = 12, 200 cells/unit, L = 0.696283 with EG's
       cell vs 0.735421 with EM.1; with EM.1, L(12) >= 0.735421, 0.737853, 0.738463 on 200, 800, 3200 cells/unit (the 800/3200
       values logged by Theorem EN's script, EN.4; this line earlier said L(14) and quoted unlogged 0.737856, 0.738466) -> EE's float
       e* = 0.738670, consistent with O(1/n).)
  EM.2 (arb, 128 bits; --full; EF.1 point cells n_p = 16000, n_m = 256000 (EF's certify_mt cell bound, cells split into parallel
       chunks); EM.1 lower k0 in {14, 20, 26} on 12800 cells/unit; EL.1 upper k0 = 26 on 12800 cells/unit): at alpha = 0.91, N = 1:
       EJ.2 (b) gives Y_n < Y_oo for n = 2565..2585 (with EJ.4: for n = 1..2585), and EJ.2 (a) holds at n1 = 2586, so Y_n > Y_oo for
       every n >= 2586.  So Y_n - Y_oo changes sign exactly once, between levels 2585 and 2586 (EL.3: [2564, 2586]).
  EM.3 (arb; default and --full; point 4000/32000, EL.1 upper k0 = 15 on 1600 cells/unit, EM.1 lower on 1600/unit): at alpha = 0.89,
       N = 1: Y_n < Y_oo for n = 2..14 and Y_n > Y_oo for every n >= 15: exactly one sign change, between 14 and 15 (EJ.4: [14, 15]).
  Default (no --full): EM.3, and at 0.91 the lower side on point 4000/32000, EM.1 lower 3200/unit: Y_n < Y_oo up to n = 2584
       (with EL.3's onset 2587, last negative level in [2584, 2586]).

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_two_sided_rate as ej  # noqa: E402
import gasket_fractional_corner_exact_0905_n3 as el  # noqa: E402

FAILS: list[str] = []
WORKERS = int(os.environ.get("EM_WORKERS", "4"))
LO, HI_LO = -25, 12                          # EG's lower grid range (x in [5^-25, 5^12] for c; kappa to 5^(N+6)); head/tail dropped
A91, A89, A88 = "0.91", "0.89", "0.88"
K0_LO91 = (14, 20, 26)
K0_LO_SMALL = (1, 2, 3, 4, 6, 8, 10, 12, 14)
FULL91 = {"n_p": 16000, "n_m": 256000, "n_l": 12800, "n_u": 12800}
DEF91 = {"n_p": 4000, "n_m": 32000, "n_l": 3200}
SMALL = {"n_p": 4000, "n_m": 32000, "n": 1600}
EJ4_LAST91 = 2564          # EJ.4: Y_n < Y_oo for n = 1..2564 at 0.91 (N = 1)
EL3_ONSET91 = 2587         # EL.3: Y_n > Y_oo for every n >= 2587 at 0.91 (N = 1)
EXPECT91 = (2585, 2586)    # EM.2 (--full): last negative level, onset
EXPECT91_DEFAULT = 2584    # default: last certified negative level
EXPECT89 = (14, 15)        # EM.3: last negative level, onset (EJ.4: [14, 15])
K0_UP89 = 15
EL_NU91 = (0.3970713, 0.3971799)   # EL.3's |nu| at 0.91 (consistency)


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _lo(v) -> float:
    return math.nextafter(float(v.lower()), -math.inf)


def _hi(v) -> float:
    return math.nextafter(float(v.upper()), math.inf)


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
def vlo_right_float(x, k0: int):
    """Floats: EM.1's right-node lower bound 1/V_k0(x) for R_k (EG.1's V_k0 at the point x)."""
    _, eg = ej._mods()
    return 1.0 / eg.vplus_float(np.asarray(x, float), k0)[0]


# ----------------------------------------------------------------- EM.1: EG.1's lower cell sums with the right node
def _cells_lo(args):
    """EG's certify_lower cell sums over cells i0 <= i < i1, every R / Rh at the right node b (EM.1); 'mixed' reproduces EG's vcell."""
    a_s, k0s, Ns, n, hi, i0, i1, mode = args
    from flint import arb, ctx
    ctx.prec = 128
    _, eg = ej._mods()
    al, five = arb(a_s), arb(5)
    ks = sorted(set(k0s) | {k + 1 for k in k0s})
    Nmax = max(Ns) if Ns else 0
    ex = max(40, max(k0s) + 2 - LO)          # node depth: P_m needed for m <= k0 + 1 at x >= 5^LO
    iD = (hi - LO) * n
    c = {k0: arb(0) for k0 in k0s}
    kap = {(N, k0): arb(0) for N in Ns for k0 in k0s}
    for i in range(i0, i1):
        xa, xb = arb(float(5.0 ** (LO + i / n))), arb(float(5.0 ** (LO + (i + 1) / n)))
        nb = eg.node_arb(xb, ks, Nmax, extra=ex)
        na = nb if mode == "right" else eg.node_arb(xa, ks, Nmax, extra=ex)
        W = (xb ** (2 - al) - xa ** (2 - al)) / (2 - al)
        if i < iD:
            nb5 = eg.node_arb(5 * xb, ks, extra=ex)
            na5 = nb5 if mode == "right" else eg.node_arb(5 * xa, ks, extra=ex)
        for k0 in k0s:
            Rl = 1 / eg.vcell(na, nb, k0)
            wb = 2 / (6 + xb / five ** (k0 + 1))
            if i < iD:
                Rhl = (five / 3) / eg.vcell(na5, nb5, k0 + 1)
                gg = wb * arb(nb["P"][0].lower()) * (1 + 1 / (arb(nb["Phi"].upper()) + 5)) * Rl * Rhl * (Rl + Rhl)
                c[k0] += W * arb(gg.lower())
            for N in Ns:
                if i < (N + 6 - LO) * n:
                    et = eg.eta_arb(N, nb["zs"])
                    if bool(et.lower() > 0):      # eta_N > 0 (BL): dropping a cell keeps a lower bound
                        kap[(N, k0)] += W * arb((wb * Rl * Rl * arb(et.lower())).lower())
    return {k: _lo(v) for k, v in c.items()}, {k: _lo(v) for k, v in kap.items()}


def certify_lower_em(a_s, k0s, Ns, n, hi=HI_LO, mode="right", workers=WORKERS):
    """c_lo(k0), L(k0) = 3 c_lo/(2(1-rho/3)), kappa_lo(N, k0) with EM.1 (certified LOWER bounds; same keys as EJ's certify_lower)."""
    arb, _ = ej._arb()
    hiK = max([hi] + [N + 6 for N in Ns])
    ncell = (hiK - LO) * n
    step = max(1, ncell // (8 * workers))
    jobs = [(a_s, tuple(k0s), tuple(Ns), n, hi, i, min(i + step, ncell), mode) for i in range(0, ncell, step)]
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_cells_lo, jobs))
    rho = arb(5) ** arb(a_s) / 3
    c = {k0: sum((arb(r[0][k0]) for r in res), arb(0)) for k0 in k0s}
    kap = {k: sum((arb(r[1][k]) for r in res), arb(0)) for k in res[0][1]}
    return {"L": {k: arb((3 * v / (2 * (1 - rho / 3))).lower()) for k, v in c.items()}, "c": c, "kappa": kap}


# ----------------------------------------------------------------- EF.1's point constants, Mt cells in parallel chunks
def _mt_chunk(args):
    """EF's certify_mt cell bound (unchanged) on cells i0 <= i < i1 (nodes recomputed per chunk; memory-light)."""
    a_s, n, lo, i0, i1 = args
    from flint import arb, ctx
    ctx.prec = 128
    ef, _ = ej._mods()
    al = arb(a_s)
    Lo, Up, ok = arb(0), arb(0), True
    a = arb(float(5.0 ** (lo + i0 / n)))
    ra = ef.rinf_arb(a)
    for i in range(i0, i1):
        b = arb(float(5.0 ** (lo + (i + 1) / n)))
        rb = ef.rinf_arb(b)
        tb = a.union(b)
        W = (b ** (1 - al) - a ** (1 - al)) / (1 - al)
        fu = ef.F_arb(tb, ra[1].union(ra[1]), al)
        fl = ef.F_arb(tb, rb[0], al)
        ok = ok and fu.is_finite() and fl.is_finite()
        Up += W * arb(fu.upper())
        if bool(fl.lower() > 0):
            Lo += W * arb(fl.lower())
        a, ra = b, rb
    return _lo(Lo), _hi(Up), ok


def certify_mt_par(a_s, n, lo=-40, hi=30, workers=WORKERS):
    """EF's certify_mt with its cells split into chunks (identical cell bound, head and tail)."""
    arb, _ = ej._arb()
    ef, _ = ej._mods()
    al = arb(a_s)
    g = arb(3).log() / arb(5).log()
    ncell = (hi - lo) * n
    step = max(1, ncell // (16 * workers))
    jobs = [(a_s, n, lo, i, min(i + step, ncell)) for i in range(0, ncell, step)]
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_mt_chunk, jobs))
    Lo = sum((arb(r[0]) for r in res), arb(0))
    Up = sum((arb(r[1]) for r in res), arb(0))
    eps, Gc = arb(float(5.0 ** lo)), arb(ef.BI4)
    c1 = (arb("0.18") * (arb("0.18") * eps).exp() * eps ** (1 - g) + Gc / 3) / (1 - Gc / 3 * eps**g)
    head = Gc * c1 * eps ** (2 * g - al) / (2 * g - al)
    tail = arb(9) / (2 * al) * arb(float(5.0 ** hi)) ** (-al / 2)
    return ef.hull(Lo, Up + head + tail), {"ok": all(r[2] for r in res)}


def certify_point_em(a_s, n_p, n_m, NS=(1,)):
    """EJ's certify_point (EF.1's profile integrals unchanged; Mt via certify_mt_par)."""
    arb, _ = ej._arb()
    ef, _ = ej._mods()
    al = arb(a_s)
    rho = arb(5) ** al / 3
    pi = ef.certify_profile_integrals(al, NS, n=n_p)
    Mt, mi = certify_mt_par(a_s, n_m)
    Dp = pi["Dp"]
    return {"al": al, "rho": rho, "g": 3 / rho**3, "Dp": Dp, "Mt": Mt, "nu": -rho * Mt / ((rho - 1) * Dp),
            "om": {N: pi["h"][N] / Dp for N in NS}, "ok": mi["ok"] and all(pi["ok_tail"][N] for N in NS)}


def lower_bound_at(P, low, N, n):
    """EJ.2 (b)'s upper bound on the braces at level n, smallest over k0 <= n in low (EJ.below_upper's formula); < 0 certifies."""
    arb, _ = ej._arb()
    rho, g = P["rho"], P["g"]
    Dhi, nhi = arb(P["Dp"].upper()), arb((-P["nu"]).upper())
    oml, omh = arb(P["om"][N].lower()), arb(P["om"][N].upper())
    best = None
    for k0 in sorted(low["L"]):
        if k0 > n or (N, k0) not in low["kappa"]:
            continue
        X = nhi**3 * g**n - low["L"][k0] / Dhi
        den = 1 - nhi * rho ** (-n)
        if not bool(den > 0):
            continue
        val = (oml * X if bool(X < 0) else omh / den * X) - 3 * low["kappa"][(N, k0)] / (2 * Dhi)
        v = float(val.upper())
        best = v if best is None else min(best, v)
    return best


def rk_check(k0: int, nodes, ncell: int = 1600) -> bool:
    """arb spot check of EM.1: on cells [b 5^(-1/ncell), b] the exact R_k (affine recursion, k = k0..k0+15) at six points is >= the
    right-node bound 1/V_k0(b) and nonincreasing."""
    arb, _ = ej._arb()
    _, eg = ej._mods()
    ok = True
    for t in nodes:
        b = arb(float(5.0**t))
        nb = eg.node_arb(b, [k0], extra=max(40, k0 + 2 - LO))
        Rlo = 1 / eg.vcell(nb, nb, k0)
        xs = [arb(float(5.0 ** (t - j / (5 * ncell)))) for j in range(5, -1, -1)]
        for k in range(k0, k0 + 16):
            vals = [ej.rk_arb(x, k) for x in xs]
            ok = ok and all(bool(v >= Rlo) for v in vals) and all(bool(vals[j] >= vals[j + 1]) for j in range(5))
    return ok


def main(full: bool = False) -> int:
    t0 = time.time()
    print(f"Theorem EM: exact sign-change levels at alpha = 0.91 and 0.89 (N = 1), with EM.1's right-node lower cell bound. Not thrust. "
          f"({'--full' if full else 'default'}, {WORKERS} workers)", flush=True)
    arb, ctx = ej._arb()
    ctx.prec = 128
    # ---- attempts to break EM.1: exact R_k vs the right-node bound; right node vs EG's two-ended cell
    check(rk_check(14, (-3.0, 0.0, 2.5, 4.0, 5.0)), "EM.1 spot check (arb): exact R_k, k = 14..29, on six points of cells [b 5^(-1/1600), b], "
          "b = 5^-3, 1, 5^2.5, 5^4, 5^5: >= 1/V_14(b) and nonincreasing across the cell")
    t1 = time.time()
    lr = certify_lower_em("0.905", (12,), (1,), 200)
    lm = certify_lower_em("0.905", (12,), (1,), 200, mode="mixed")
    print(f"   0.905, k0 = 12, 200 cells/unit: L(12) >= {float(lr['L'][12].lower()):.6f} (EM.1) vs {float(lm['L'][12].lower()):.6f} "
          f"(EG's two-ended cell) ({time.time() - t1:.1f}s)", flush=True)
    check(bool(lr["L"][12] >= lm["L"][12]) and bool(lr["kappa"][(1, 12)] >= lm["kappa"][(1, 12)]),
          "EM.1 dominates EG's cell bound on the same grid (both certified lower bounds)")
    # ---- serial vs parallel Mt (identical cell bound)
    ef, _ = ej._mods()
    Ms, _ = ef.certify_mt(arb(A91), n=200)
    Mp, _ = certify_mt_par(A91, 200)
    check(bool(Mp.overlaps(Ms)) and abs(float(Mp.lower()) - float(Ms.lower())) < 1e-12 and abs(float(Mp.upper()) - float(Ms.upper())) < 1e-12,
          "certify_mt_par reproduces EF's certify_mt at n_m = 200 (same enclosure to 1e-12)")
    # ---- EM.3: alpha = 0.89
    t1 = time.time()
    P = certify_point_em(A89, SMALL["n_p"], SMALL["n_m"])
    check(P["ok"], f"EF.1 routines valid at alpha = {A89}")
    up = el.certify_upper_el(A89, K0_UP89, SMALL["n"], Ns=(1,))
    check(up[K0_UP89]["cond"], f"EJ.1 alpha = {A89}: (3/rho) beta(p^k0(1)) <= 1 at k0 = {K0_UP89}")
    on = ej.onset_upper(P, up, 1)
    n1 = on[0] if on else None
    low = certify_lower_em(A89, K0_LO_SMALL, (1,), SMALL["n"])
    win = [n for n in range(1, (n1 or 30)) if ej.below_upper(P, low, 1, n)]
    print(f"   EM.3 alpha = {A89}: |nu| in {ej.fmt(-P['nu'], 7)}, U({K0_UP89}) <= {float(up[K0_UP89]['U'].upper()):.6f}, kappa_hi(1) <= "
          f"{float(up[K0_UP89]['kappa'][1]['hi'].upper()):.6f}; Phi({n1}) >= {float(on[2].lower()):.4g} > 0 -> Y_n > Y_oo for EVERY n >= {n1}; "
          f"Y_n < Y_oo certified for n in {win[:1]}..{win[-1:]} ({time.time() - t1:.1f}s)", flush=True)
    check(n1 == EXPECT89[1], f"EM.3 alpha = {A89}, N = 1: Y_n > Y_oo for every n >= {EXPECT89[1]} (EJ.4: >= 16)")
    check(win == list(range(2, EXPECT89[0] + 1)), f"EM.3 alpha = {A89}, N = 1: Y_n < Y_oo for n = 2..{EXPECT89[0]} (reproduces EJ.4's window [2, 14])")
    if n1 == EXPECT89[1] and win and win[-1] == EXPECT89[0]:
        print(f"     -> exactly one sign change at alpha = {A89} (N = 1, n >= 2), between levels {EXPECT89[0]} and {EXPECT89[1]} (EJ.4: [14, 15])")
    # ---- (OBSERVED) alpha = 0.88, level 9 stays undecided
    P8 = certify_point_em(A88, SMALL["n_p"], SMALL["n_m"])
    up8 = el.certify_upper_el(A88, 9, SMALL["n"], Ns=(1,))
    on8 = ej.onset_upper(P8, up8, 1)
    low8 = certify_lower_em(A88, tuple(k for k in K0_LO_SMALL if k <= 9), (1,), SMALL["n"])
    b9 = lower_bound_at(P8, low8, 1, 9)
    print(f"   (OBSERVED) alpha = {A88}: EL.1 upper at k0 = 9 gives onset {on8[0]} (Phi(9) <= {float(on8[3].upper()):.4g}); EM.1 lower-side "
          f"bound at n = 9 is {b9:+.4g} (not < 0): level 9 undecided", flush=True)
    check(on8[0] == 10 and b9 > 0, "KEPT LIMIT: alpha = 0.88, n = 9 undecided by both sides (bracket [8, 9] unchanged)")
    # ---- EM.2: alpha = 0.91
    cfg = FULL91 if full else DEF91
    t1 = time.time()
    P = certify_point_em(A91, cfg["n_p"], cfg["n_m"])
    nu_lo, nu_hi = float((-P["nu"]).lower()), float((-P["nu"]).upper())
    print(f"   point 0.91 ({cfg['n_p']}/{cfg['n_m']}): D'_oo in {ej.fmt(P['Dp'], 9)}, |nu| in {ej.fmt(-P['nu'], 9)}, 1-delta_1 in "
          f"{ej.fmt(P['om'][1], 9)} ({time.time() - t1:.1f}s)", flush=True)
    check(P["ok"], f"EF.1 routines valid at alpha = {A91}")
    check(EL_NU91[0] - 1e-6 <= nu_lo and nu_hi <= EL_NU91[1] + 1e-6, f"consistency: |nu_0.91| bracket nested in EL.3's {list(EL_NU91)} (EL rounded to 7 digits; tol 1e-6)")
    t1 = time.time()
    low = certify_lower_em(A91, K0_LO91, (1,), cfg["n_l"])
    L26, k26 = low["L"][26], low["kappa"][(1, 26)]
    print(f"   lower (EM.1, n = {cfg['n_l']}): L(14, 20, 26) >= " + ", ".join(f"{float(low['L'][k].lower()):.8f}" for k in K0_LO91)
          + f", kappa_lo(1, 26) >= {float(k26.lower()):.8f} ({time.time() - t1:.1f}s)", flush=True)
    neg = [n for n in range(EJ4_LAST91 + 1, EL3_ONSET91) if ej.below_upper(P, low, 1, n)]
    last = neg[-1] if neg else None
    contig = neg == list(range(EJ4_LAST91 + 1, (last or EJ4_LAST91) + 1))
    for n in range(2583, 2587):
        b = lower_bound_at(P, low, 1, n)
        print(f"     lower-side bound at n = {n}: {b:+.4g}" + ("  (Y_n < Y_oo)" if b < 0 else "  (undecided)"), flush=True)
    check(contig and last == (EXPECT91[0] if full else EXPECT91_DEFAULT),
          f"EM.2 alpha = {A91}, N = 1: Y_n < Y_oo for n = {EJ4_LAST91 + 1}..{EXPECT91[0] if full else EXPECT91_DEFAULT} (with EJ.4: for n = 1..)")
    if full:
        t1 = time.time()
        up = el.certify_upper_el(A91, 26, cfg["n_u"], Ns=(1,))
        r = up[26]
        check(r["cond"], f"EJ.1 alpha = {A91}: (3/rho) beta(p^k0(1)) <= 1 at k0 = 26")
        on = ej.onset_upper(P, up, 1)
        n1 = on[0] if on else None
        print(f"   upper (EL.1, k0 = 26, n = {cfg['n_u']}): U(26) <= {float(r['U'].upper()):.8f}, kappa_hi(1, 26) <= "
              f"{float(r['kappa'][1]['hi'].upper()):.8f}; Phi({n1}) >= {float(on[2].lower()):.4g} > 0, Phi({n1 - 1}) <= "
              f"{float(on[3].upper()):.4g} -> Y_n > Y_oo for EVERY n >= {n1} ({time.time() - t1:.1f}s)", flush=True)
        check(n1 == EXPECT91[1], f"EM.2 alpha = {A91}, N = 1: Y_n > Y_oo for every n >= {EXPECT91[1]} (EL.3: >= {EL3_ONSET91})")
        check(bool(L26 <= r["U"]) and bool(k26 <= r["kappa"][1]["hi"]), "two-sided: L(26) <= U(26) and kappa_lo(1, 26) <= kappa_hi(1, 26)")
        c = ef.critical_floats(0.91, Ns=(1,))
        check(float(L26.lower()) < c["e"] < float(r["U"].upper()) and float(k26.lower()) < c["kappa"][1] < float(r["kappa"][1]["hi"].upper()),
              f"EE's float limits e* = {c['e']:.7f}, kappa_1 = {c['kappa'][1]:.7f} lie inside [L(26), U(26)], [kappa_lo, kappa_hi]")
        if n1 == EXPECT91[1] and last == EXPECT91[0]:
            print(f"     -> exactly one sign change at alpha = {A91} (N = 1), between levels {EXPECT91[0]} and {EXPECT91[1]} (EL.3: [2564, 2586])")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EM checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(full="--full" in sys.argv))
