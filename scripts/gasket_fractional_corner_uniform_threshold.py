#!/usr/bin/env python3
"""Theorem EH: the explicit threshold on the whole range [4gamma/3, 1).  EG certified "Y_n < Y_oo for every N >= 1 and every
n >= N" at the two exponents 4gamma/3 and 0.95 only and kept the open item "Y_n < Y_oo for every n >= N on all of [4gamma/3, 1)
(needs nu_alpha, D'_oo as certified functions of alpha)".  Here alpha runs over closed cells [a1, a2] covering [4gamma/3, 1]:

  EH.1 (THEOREM; monotonicity in alpha, every cell [a1, a2] in (gamma, 1]).  Every alpha-dependence of EF.1's and EG.1's integrals
       is a weight x^-alpha (resp. x^(1-alpha)) times an alpha-FREE nonnegative profile factor (G, G eta_1(Phi(x/5)), EG.1's
       c_lo / kappa_lo integrands), except F(t, r) = r(omega^0/(1-br) - 1) in Mt_oo, which is NONDECREASING in alpha
       (omega^0 = (1+q/5)^alpha w1(q), log1p(q/5) >= 0).  On an x-cell inside (0, 1] the weight is nondecreasing in alpha, on an
       x-cell inside [1, oo) nonincreasing (the grid has the exact node x = 1).  Hence, for EVERY alpha in [a1, a2]:
         D'_oo(alpha) >= D_lo,  Mt_oo(alpha) <= M_hi,  c_lo(alpha; k0) >= c_lo,  kappa_lo(alpha; 1, 1) >= k_lo,  h_1(alpha) <= h_hi,
       each a finite sum of cell weights taken at the endpoint that makes it smaller (larger) times the alpha-free cell values,
       plus heads/tails enclosed with alpha as the ball [a1, a2]; rho/(rho-1) is decreasing and 1/(1-rho/3) increasing in rho.
  EH.2 (THEOREM; EG.2 on a cell).  |nu|^3 D'_oo = (rho Mt/(rho-1))^3 / D'^2 <= T_hi, |nu| <= nu_hi.  For alpha >= 4gamma/3,
       3/rho^3 <= f = min(1, 3/rho(a1)^3).  If nu_hi < rho(a1) (so 1 + nu rho^-n > 0, n >= 1), f^2 T_hi < max(L1_lo, L2_lo) and
       [f T_hi < L1_lo  or  (1-delta_1)/(1+nu/rho) (|nu|^3 f - L(1)/D') - 1.5 kappa_lo(1,1)/D' < 0 (upper bound)], then Y_n < Y_oo
       for every N >= 1 and every n >= N, at every alpha in the cell (EG.2: n >= 2 from k0 = 1 or 2, n = 1 = N from k0 = 1).
  EH.3 (arb, 128 bits).  The cells cover [4gamma/3, 1] (first endpoint the ball 4 log3/(3 log5)); every cell passes EH.2.  So
       Y_n < Y_oo for EVERY alpha in [4gamma/3, 1), EVERY N >= 1 and EVERY n >= N.
  EH.4 (arb, below 4gamma/3; points 0.88, 0.89, 0.9, 0.905, 0.91).  3/rho^3 > 1 there, but EG.2 still gives Y_n < Y_oo whenever
       n >= max(N, k0) and (3/rho^3)^n |nu|^3 D'_oo < L(k0): certified initial windows of n (0.9: 3..8, 0.905: 2..21, 0.91: 2..926 on
       --full), while EC.3/ED.3 give Y_n > Y_oo for all large n: Y_n - Y_oo changes sign in n there.
  EH.5 (OBSERVED, floats).  EF/EG's float routines on an alpha grid: the margins of EH.2 and the alpha at which the n = 1
       check stops needing the eta term.

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time

import numpy as np

FAILS: list[str] = []
RESULT: dict = {}
GAMMA = math.log(3) / math.log(5)
K0S = (1, 2, 3, 4, 6, 8)


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


def _mods():
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    import gasket_fractional_corner_critical_sign as ef
    import gasket_fractional_corner_explicit_threshold as eg
    import gasket_fractional_corner_remainder_limit as ee
    return ef, eg, ee


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
def survey_float(alphas):
    """Floats (OBSERVED): |nu|^3 D'_oo, 3/rho^3, L(1), L(2), 1-delta_1, kappa_lo(1,1) and the n = 1 eta-route bound by EF/EG's
    float routines."""
    ef, eg, _ = _mods()
    rows = []
    for al in alphas:
        c = ef.critical_floats(al, Ns=(1,))
        rho = 5**al / 3
        T = -c["nu"] ** 3 * c["Dp"]
        l1 = eg.lower_float(al, 1, Ns=(1,))
        L2 = eg.lower_float(al, 2)["L"]
        f = min(1.0, 3 / rho**3)
        S11 = c["om"][1] * (-c["nu"] ** 3 * f - l1["L"] / c["Dp"]) / (1 + c["nu"] / rho) - 1.5 * l1["kappa"][1] / c["Dp"]
        rows.append({"alpha": al, "rho": rho, "Dp": c["Dp"], "nu": c["nu"], "T": T, "f": f, "L1": l1["L"], "L2": L2,
                     "om1": c["om"][1], "k11": l1["kappa"][1], "S11": S11})
    return rows


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    import sympy as sp
    x, a, X = sp.symbols("x alpha X", positive=True)
    check(sp.simplify(sp.diff(x ** (-a), a) + sp.log(x) * x ** (-a)) == 0,
          "d/dalpha x^-alpha = -log(x) x^-alpha: >= 0 on (0,1], <= 0 on [1,oo) (same for x^(1-alpha)); cell weights monotone in alpha")
    q = sp.Symbol("q", positive=True)
    w1 = (1 + q / 5) * (1 + q / 2) / ((1 + 2 * q / 5) * (1 + q / 6))
    om0 = (1 + q / 5) ** (1 + a) * (1 + q / 2) / ((1 + 2 * q / 5) * (1 + q / 6))
    check(sp.simplify(sp.diff(om0, a) - sp.log(1 + q / 5) * om0) == 0 and sp.simplify(om0 - (1 + q / 5) ** a * w1) == 0,
          "omega^0 = (1+q/5)^alpha w1(q), d omega^0/dalpha = log(1+q/5) omega^0 >= 0: F(t,r) = r(omega^0/(1-br) - 1) nondecreasing in alpha")
    f = X ** (1 - a) / (1 - a)
    check(sp.simplify(sp.diff(sp.log(f), a) - (-sp.log(X) + 1 / (1 - a))) == 0,
          "head weight X0^(1-alpha)/(1-alpha): d/dalpha log = -log X0 + 1/(1-alpha) > 0 for X0 < 1 (nondecreasing in alpha)")
    r = sp.Symbol("r", positive=True)
    check(sp.simplify(sp.diff(r / (r - 1), r) + 1 / (r - 1) ** 2) == 0 and sp.simplify(sp.diff(1 / (1 - r / 3), r) - 1 / (3 * (1 - r / 3) ** 2)) == 0,
          "rho/(rho-1) decreasing, 1/(1-rho/3) increasing in rho; rho = 5^alpha/3 increasing in alpha")
    a43 = 4 * sp.log(3) / (3 * sp.log(5))
    check(sp.simplify(sp.expand_log((5**a43 / 3) ** 3, force=True) - 3) == 0 and bool(sp.N(a43, 30) < 1),
          "rho^3 = 3 at 4gamma/3 < 1, so 3/rho^3 <= 1 on [4gamma/3, 1)")
    nu, D, M, rh = sp.symbols("nu D M rho", positive=True)
    check(sp.simplify((rh * M / ((rh - 1) * D)) ** 3 * D - (rh * M / (rh - 1)) ** 3 / D**2) == 0,
          "|nu|^3 D'_oo = (rho Mt/(rho-1))^3 / D'_oo^2 (EF.1: nu = -rho Mt/((rho-1) D'_oo))")


# ----------------------------------------------------------------- (B) arb
def _arb():
    from flint import arb, ctx
    return arb, ctx


def wint(a, b, e, is_m1=False):
    """int_a^b x^e dx (arb); e = -1 handled by log."""
    arb, _ = _arb()
    if is_m1:
        return (b / a).log()
    return (b ** (e + 1) - a ** (e + 1)) / (e + 1)


def precompute(n_p: int, n_m: int, n_g: int, lo: int = -25, hiD: int = 8, lot: int = -40, hit: int = 30):
    """alpha-free cell data (arb):
    (i)  x-grid 5^(lo + i/n_p): G at nodes (EE's prof_arb), G eta_1(Phi(x/5)), A; cell values G(b) (lower), G(a) (upper) etc.
    (ii) t-grid 5^(lot + i/n_m): r_oo enclosures (EF's rinf_arb); per cell l5 = log1p(q/5) and L0 = l5 + log1p(q/2) - log1p(2q/5)
         - log1p(q/6) - log1p(-b r^+(a)) on the t-cell ball, so F <= r^+(a) expm1(alpha l5 + L0) (upper ends).
    (iii) x-grid 5^(lo + i/n_g): EG.1's c_lo(k0) integrand lower cell values (k0 = 1, 2) and kappa_lo(1, 1)'s."""
    arb, _ = _arb()
    ef, eg, ee = _mods()
    five = arb(5)
    out = {"lo": lo, "hiD": hiD, "lot": lot, "hit": hit}
    # (i)
    hi = max(hiD, 7)
    xs = [arb(float(5.0 ** (lo + i / n_p))) for i in range((hi - lo) * n_p + 1)]
    G, E1, A, Z1 = [], [], [], []
    for xv in xs:
        g, aa, zs = ee.prof_arb(xv)
        G.append(g)
        A.append(aa)
        E1.append(g * ee.eta_arb(1, zs))
        Z1.append(zs[1])
    out["x"] = {"xs": xs, "G": G, "E1": E1, "A": A, "Z1": Z1, "n": n_p, "i1": -lo * n_p, "iD": (hiD - lo) * n_p,
                "iB": (hiD - 1 - lo) * n_p, "iN": (1 + 6 - lo) * n_p}
    check(bool(xs[-lo * n_p] == 1), "x-grid contains the exact node x = 1 (cells lie in (0,1] or [1,oo))")
    # (ii)
    ts = [arb(float(5.0 ** (lot + i / n_m))) for i in range((hit - lot) * n_m + 1)]
    rr = [ef.rinf_arb(tv) for tv in ts]
    l5s, L0s, l5m, L0ms, rms = [], [], [], [], []
    ok = True
    for i in range(len(ts) - 1):
        tb = ts[i].union(ts[i + 1])
        rp = arb(rr[i][1].upper())
        q = 2 * tb / (5 + (25 + 4 * tb).sqrt())
        b = 2 * tb / (q + 6)
        l5 = (q / 5).log1p()
        L0 = l5 + (q / 2).log1p() - (2 * q / 5).log1p() - (q / 6).log1p() - (-(b * rp)).log1p()
        ok = ok and l5.is_finite() and L0.is_finite() and bool(b * rp < 1)
        l5s.append(arb(l5.upper()))
        L0s.append(arb(L0.upper()))
        rm = arb(rr[i + 1][0].lower())                     # r_oo >= r^-(b) on the cell (r_oo nonincreasing)
        L0m = l5 + (q / 2).log1p() - (2 * q / 5).log1p() - (q / 6).log1p() - (-(b * rm)).log1p()
        l5m.append(arb(l5.lower()))
        L0ms.append(arb(L0m.lower()))
        rms.append(rm)
        rr[i] = (rr[i][0], rp)
    out["t"] = {"ts": ts, "rp": [rr[i][1] for i in range(len(ts) - 1)], "l5": l5s, "L0": L0s, "l5m": l5m, "L0m": L0ms, "rm": rms, "n": n_m, "i1": -lot * n_m, "ok": ok}
    check(ok and bool(ts[-lot * n_m] == 1), "t-grid: every cell has 1 - b r^+ > 0, finite enclosures; exact node t = 1")
    # (iii)
    xg = [arb(float(5.0 ** (lo + i / n_g))) for i in range((hiD - lo) * n_g + 1)]
    ks = sorted(set(K0S) | {k + 1 for k in K0S})
    nd = [eg.node_arb(xv, ks, 1) for xv in xg]
    nd5 = [eg.node_arb(5 * xv, ks) for xv in xg]
    gc = {k0: [] for k0 in K0S}
    gk = []
    for i in range(len(xg) - 1):
        na, nb = nd[i], nd[i + 1]
        for k0 in K0S:
            Rl = 1 / eg.vcell(na, nb, k0)
            Rhl = (five / 3) / eg.vcell(nd5[i], nd5[i + 1], k0 + 1)
            g = 2 / (6 + xg[i + 1] / five ** (k0 + 1)) * arb(nb["P"][0].lower()) * (1 + 1 / (arb(nb["Phi"].upper()) + 5)) * Rl * Rhl * (Rl + Rhl)
            gc[k0].append(arb(max(0.0, float(g.lower()))) if not bool(g.lower() > 0) else arb(g.lower()))
        Rl = 1 / eg.vcell(na, nb, 1)
        et = eg.eta_arb(1, nb["zs"])
        if bool(et.lower() > 0):
            gk.append(arb((2 / (6 + xg[i + 1] / five**2) * Rl * Rl * arb(et.lower())).lower()))
        else:
            gk.append(arb(0))
    out["g"] = {"xs": xg, "c": gc, "k": gk, "i1": -lo * n_g}
    check(all(bool(v >= 0) for k0 in K0S for v in gc[k0]) and all(bool(v >= 0) for v in gk),
          "EG.1 cell values (c_lo, kappa_lo integrands) nonnegative (negative lower ends replaced by 0)")
    return out


def cell_bounds(pre, A1, A2, need_hi: bool, all_k0: bool = False):
    """EH.1 bounds valid for every alpha in [A1, A2] (A1, A2 arb balls containing the endpoints; A2 may be exactly 1)."""
    arb, _ = _arb()
    ef, _, _ = _mods()
    one = bool(A2 == 1)
    g = arb(3).log() / arb(5).log()
    Ab = arb(A1.lower()).union(arb(A2.upper()))
    rho1 = arb((arb(5) ** A1 / 3).lower())          # lower bound of rho on the cell
    rho2 = arb((arb(5) ** A2 / 3).upper())          # upper bound
    lo_e = lambda a, b, below1: (A1, False) if below1 else (A2, one)     # endpoint of the MIN weight of x^-alpha
    hi_e = lambda a, b, below1: (A2, one) if below1 else (A1, False)
    # --- D'_oo lower (and upper)
    X = pre["x"]
    xs, G = X["xs"], X["G"]

    def wsum(vals, i0, i1, mode, shift=0, e_off=0):
        s = arb(0)
        for i in range(i0, i1):
            below1 = i < X["i1"]
            Ae, m1 = (lo_e if mode == "min" else hi_e)(None, None, below1)
            w = wint(xs[i], xs[i + 1], -Ae, is_m1=m1)
            w = arb(w.lower()) if mode == "min" else arb(w.upper())
            s += w * vals(i)
        return s

    Dl_cells = wsum(lambda i: arb(G[i + 1].lower()), 0, X["iD"], "min")
    Ib_lo = wsum(lambda i: arb(G[i + 1].lower()), X["iB"], X["iD"], "min")
    X0 = xs[0]
    head_lo = arb((X0 ** (1 - A1) / (1 - A1)).lower()) * arb(G[0].lower())
    Dlo = arb((Dl_cells + head_lo + Ib_lo / (rho2 - 1)).lower())
    out = {"Dlo": Dlo}
    if need_hi:
        assert not one
        Du_cells = wsum(lambda i: arb(G[i].upper()), 0, X["iD"], "max")
        Ib_hi = wsum(lambda i: arb(G[i].upper()), X["iB"], X["iD"], "max")
        tau = arb(2) / 3 * arb(X["A"][X["iB"]].upper())
        Dhi = Du_cells + arb((X0 ** (1 - A2) / (1 - A2)).upper()) / 3 + Ib_hi * (1 + tau) / (rho1 - 1 - tau)
        out["Dhi"] = arb(Dhi.upper())
        E1 = X["E1"]
        h_cells = wsum(lambda i: arb(E1[i].upper()), 0, X["iN"], "max")
        eta0 = 1 - arb(3) / 2 * (arb(3) / 5)
        Xn = xs[X["iN"]]
        P = arb(X["Z1"][X["iN"]].lower())
        tailN = 2 / P * arb(ef.BI4) * Xn ** (g - Ab) / (Ab - g)
        hhi = h_cells + arb((X0 ** (1 - A2) / (1 - A2)).upper()) * eta0 / 3 + tailN
        out["hhi"] = arb(hhi.upper())
        out["ok_tailN"] = bool(P > arb(10) ** 20)
    # --- Mt_oo upper
    T = pre["t"]
    ts = T["ts"]
    Mu, Ml = arb(0), arb(0)
    for i in range(len(ts) - 1):
        below1 = i < T["i1"]
        Ae, m1 = hi_e(None, None, below1)
        w = arb(wint(ts[i], ts[i + 1], -Ae, is_m1=m1).upper())
        F = T["rp"][i] * (A2 * T["l5"][i] + T["L0"][i]).expm1()
        Mu += w * arb(F.upper())
        if need_hi:
            Am, m1m = lo_e(None, None, below1)
            wl = arb(wint(ts[i], ts[i + 1], -Am, is_m1=m1m).lower())
            Fl = T["rm"][i] * (A1 * T["l5m"][i] + T["L0m"][i]).expm1()
            if bool(Fl.lower() > 0):
                Ml += wl * arb(Fl.lower())
    eps = ts[0]
    Gc = arb(ef.BI4)
    c1 = (arb("0.18") * (arb("0.18") * eps).exp() * eps ** (1 - g) + Gc / 3) / (1 - Gc / 3 * eps**g)
    head = Gc * c1 * eps ** (2 * g - Ab) / (2 * g - Ab)
    tail = arb(9) / (2 * Ab) * ts[-1] ** (-Ab / 2)
    Mhi = arb((Mu + head + tail).upper())
    out["Mhi"] = Mhi
    if need_hi:
        out["Mlo"] = arb(Ml.lower())
        rho2l = arb((arb(5) ** A2 / 3).upper())
        out["nu_lo"] = arb((rho2l / (rho2l - 1) * out["Mlo"] / out["Dhi"]).lower())   # |nu| >= (rho/(rho-1) at rho2) Mt_lo / D'_hi
    # --- EG.1 lower bounds
    Gd = pre["g"]
    xg = Gd["xs"]

    def gsum(vals):
        s = arb(0)
        for i in range(len(xg) - 1):
            below1 = i < Gd["i1"]
            if below1:
                w = wint(xg[i], xg[i + 1], 1 - A1)
            else:
                w = (xg[i + 1] - xg[i]) if one else wint(xg[i], xg[i + 1], 1 - A2)
            s += arb(w.lower()) * vals[i]
        return arb(s.lower())

    out["k11"] = gsum(Gd["k"])
    lfac = 3 / (2 * (1 - rho1 / 3))
    out["L"] = {k0: arb((lfac * gsum(Gd["c"][k0])).lower()) for k0 in (K0S if all_k0 else (1, 2))}
    out["L1"], out["L2"] = out["L"][1], out["L"][2]
    # --- EH.2
    Cabs = rho1 / (rho1 - 1) * Mhi
    out["nu_hi"] = arb((Cabs / Dlo).upper())
    out["T_hi"] = arb((Cabs**3 / Dlo**2).upper())
    f = 3 / rho1**3
    out["f"] = arb(1) if bool(f.upper() >= 1) else arb(f.upper())
    out["f_raw"] = arb(f.upper())                   # 3/rho^3 <= f_raw on the cell (used below 4gamma/3, where it is > 1)
    out["rho1"] = rho1
    return out


def decide(cb, need_hi):
    arb, _ = _arb()
    f, T = cb["f"], cb["T_hi"]
    ok_nu = bool(cb["nu_hi"] < cb["rho1"])
    ok2 = bool(f * f * T < cb["L2"]) or bool(f * f * T < cb["L1"])
    ok1 = bool(f * T < cb["L1"])
    route = "L1"
    S11 = None
    if not ok1 and need_hi:
        om_hi = cb["hhi"] / cb["Dlo"]
        par = cb["nu_hi"] ** 3 * f - cb["L1"] / cb["Dhi"]
        par = arb(max(0.0, float(par.upper()))) if bool(par.upper() > 0) else arb(0)
        S11 = om_hi * par / (1 - cb["nu_hi"] / cb["rho1"]) - 3 * cb["k11"] / (2 * cb["Dhi"])
        ok1 = bool(S11 < 0) and cb.get("ok_tailN", False)
        route = "eta"
    return ok_nu and ok2 and ok1, route, S11


def alpha_cells(step: float):
    """Endpoints: the ball 4 log3/(3 log5), then exact binary-friendly decimals k*step, then exactly 1."""
    arb, _ = _arb()
    a43 = 4 * arb(3).log() / (3 * arb(5).log())
    pts = [a43]
    k = math.floor(float(a43.upper()) / step) + 1
    while k * step < 1 - 1e-12:
        pts.append(arb(k) * arb(step))           # step = 1/2^j or decimal; ball contains k*step
        k += 1
    pts.append(arb(1))
    return pts


def main(full: bool = False) -> int:
    t0 = time.time()
    print("Theorem EH: Y_n < Y_oo for every alpha in [4gamma/3, 1), every N >= 1 and every n >= N. Not thrust.")
    symbolic_checks()
    arb, ctx = _arb()
    ctx.prec = 128
    n_p, n_m, n_g, step = (1200, 3000, 300, 0.0025) if full else (200, 300, 60, 0.005)
    pre = precompute(n_p, n_m, n_g)
    print(f"   alpha-free precompute ({time.time() - t0:.1f}s): x-grid n = {n_p}, t-grid n = {n_m}, EG.1 grid n = {n_g}")
    pts = alpha_cells(step)
    rows = []
    allok = True
    eta_hi = None
    queue = [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
    depth_split = 0
    while queue:
        A1, A2 = queue.pop(0)
        need_hi = not bool(A2 == 1)
        cb = cell_bounds(pre, A1, A2, need_hi)
        ok, route, S11 = decide(cb, need_hi)
        if not ok and depth_split < 40 and float((A2 - A1).upper()) > 1e-4:
            mid = arb((float(A1.mid()) + float(A2.mid())) / 2)
            queue[:0] = [(A1, mid), (mid, A2)]
            depth_split += 1
            continue
        allok = allok and ok
        if route == "eta":
            eta_hi = float(A2.upper())
        rows.append((float(A1.lower()), float(A2.upper()), cb, route, S11, ok))
    check(abs(rows[0][0] - 4 * GAMMA / 3) < 1e-12 and rows[-1][1] == 1.0 and all(abs(rows[i][1] - rows[i + 1][0]) < 1e-15 or rows[i][1] >= rows[i + 1][0] for i in range(len(rows) - 1)),
          f"EH.3: {len(rows)} closed alpha-cells cover [4gamma/3, 1] (consecutive endpoints shared; first endpoint the ball 4log3/(3log5))")
    show = rows[:: max(1, len(rows) // 12)]
    for r in show + ([rows[-1]] if rows[-1] is not show[-1] else []):
        a1, a2, cb, route, S11, ok = r
        extra = f", (3/rho^2)(Y_1-Y_oo) <= {float(S11.upper()):.6f}" if S11 is not None else ""
        if "nu_lo" in cb:
            extra += f"; |nu| in [{float(cb['nu_lo'].lower()):.5f}, {float(cb['nu_hi'].upper()):.5f}]"
        print(f"   [{a1:.6f}, {a2:.6f}]: D' >= {float(cb['Dlo'].lower()):.5g}, Mt <= {float(cb['Mhi'].upper()):.6f}, |nu| <= "
              f"{float(cb['nu_hi'].upper()):.5f}, |nu|^3 D' <= {float(cb['T_hi'].upper()):.5f}, f = {float(cb['f'].upper()):.5f}, "
              f"L(1) >= {float(cb['L1'].lower()):.5f}, L(2) >= {float(cb['L2'].lower()):.5f}, n=1 via {route}{extra} {'ok' if ok else 'FAIL'}")
    check(allok, "EH.3: every cell passes EH.2 (nu_hi < rho(a1); f^2 T_hi < L(1) or L(2); n = 1: f T_hi < L(1) or the eta route < 0)")
    worst2 = min(float(r[2]["L2"].lower()) - float((r[2]["f"] ** 2 * r[2]["T_hi"]).upper()) for r in rows)
    print(f"   smallest certified margin L(2) - f^2 |nu|^3 D' over the cells: {worst2:.5f}; eta route needed up to alpha = {eta_hi}")
    m1 = min(float(r[2]["L1"].lower()) - float((r[2]["f"] * r[2]["T_hi"]).upper()) for r in rows if r[3] == "L1")
    se = max(float(r[4].upper()) for r in rows if r[3] == "eta")
    nl = sum(1 for r in rows if r[3] == "eta")
    print(f"   n = 1: {len(rows) - nl} cells by L(1) alone (smallest margin L(1) - f|nu|^3 D' = {m1:.5f}), {nl} cells by the eta route "
          f"(largest bound on (3/rho^2)(Y_1 - Y_oo): {se:.6f})")
    RESULT["m1"], RESULT["se"], RESULT["worst2"] = m1, se, worst2
    RESULT["rows"] = rows
    RESULT["allok"] = allok
    RESULT["eta_hi"] = eta_hi
    RESULT["first"] = rows[0]
    # consistency with EG.3 at 4gamma/3 and EG.4 at 0.95
    c0 = rows[0][2]
    check(bool(c0["T_hi"] >= arb("0.374448")) and bool(c0["L2"] <= arb("0.4864")) and bool(c0["L1"] <= arb("0.2887")),
          "consistency: first cell's |nu|^3 D' upper bound >= EG.3's lower end 0.374448, L(1), L(2) below EG.5's floats at 4gamma/3")
    # ---- EH.4 below 4gamma/3: an initial window of n with Y_n < Y_oo (EG.2 for general rho), while Y_n > Y_oo for all large n (ED.3)
    RESULT["below"] = {}
    for a_s in ("0.88", "0.89", "0.9", "0.905", "0.91"):
        a = arb(a_s)
        cb = cell_bounds(pre, a, a, True, all_k0=True)
        f, T = cb["f_raw"], cb["T_hi"]
        win = []
        for n in range(1, 3001):
            k_ok = [k0 for k0 in K0S if k0 <= n and bool(f ** n * T < cb["L"][k0])]
            if k_ok:
                win.append((n, min(k_ok)))
        ok_nu = bool(cb["nu_hi"] < cb["rho1"])
        nlo = win[0][0] if win else None
        nhi = win[-1][0] if win else None
        contiguous = bool(win) and [w[0] for w in win] == list(range(nlo, nhi + 1))
        RESULT["below"][a_s] = {"win": (nlo, nhi), "cb": cb}
        print(f"   alpha = {a_s}: 3/rho^3 <= {float(f.upper()):.6f}, |nu| in [{float(cb['nu_lo'].lower()):.5f}, {float(cb['nu_hi'].upper()):.5f}], "
              f"|nu|^3 D' <= {float(T.upper()):.5f}, L(k0) >= " + ", ".join(f"{k}: {float(cb['L'][k].lower()):.5f}" for k in K0S)
              + (f"; Y_n < Y_oo certified for every N >= 1 and every n in [max(N, {nlo}), {nhi}]" if win else "; no window certified"))
        check(ok_nu and (not win or contiguous), f"EH.4 alpha = {a_s}: 1 + nu rho^-n > 0; certified n-window contiguous ({nlo}..{nhi})")
    w9 = RESULT["below"]["0.905"]["win"]
    w91 = RESULT["below"]["0.91"]["win"]
    check(w9[0] is not None and w9[1] >= w9[0] + 2 and w91[0] is not None and w91[1] >= 100,
          f"EH.4: windows with Y_n < Y_oo at alpha = 0.905 ({w9[0]}..{w9[1]}) and 0.91 ({w91[0]}..{w91[1]}); with EC.3/ED.3 (Y_n > Y_oo for all "
          "large n on (gamma, 4gamma/3)) Y_n - Y_oo changes sign in n there, for every N <= the window's end")
    # attempts to break it: EF/EG's own fixed-alpha arb routines at points inside cells must obey the cell bounds
    ef, eg, _ = _mods()
    for a_s in ("0.9125", "0.93", "0.9625", "0.9975"):
        a = arb(a_s)
        cells = [r for r in rows if r[0] <= float(a_s) <= r[1]]
        D = ef.certify_profile_integrals(a, Ns=(1,), n=100)["Dp"]
        M, _ = ef.certify_mt(a, n=100)
        ra = arb(5) ** a / 3
        nua = ra * M / ((ra - 1) * D)
        okc = all(bool(r[2]["Dlo"] <= D.upper()) and bool(r[2]["Mhi"] >= M.lower()) and r[2]["L1"].upper() <= eg.lower_float(float(a_s), 1)["L"] * (1 + 1e-6) and
                  (("nu_lo" not in r[2]) or (bool(r[2]["nu_lo"] <= nua.upper()) and bool(r[2]["nu_hi"] >= nua.lower()))) for r in cells)
        print(f"   fixed alpha = {a_s} (EF/EG routines, coarse): D' in [{float(D.lower()):.5g}, {float(D.upper()):.5g}], "
              f"Mt in [{float(M.lower()):.6f}, {float(M.upper()):.6f}], |nu| in [{float(nua.lower()):.5f}, {float(nua.upper()):.5f}]; float L(1) = {eg.lower_float(float(a_s), 1)['L']:.5f}")
        check(okc and len(cells) >= 1, f"attempt to break: the {len(cells)} cell(s) containing alpha = {a_s} bound EF's fixed-alpha enclosures of D', Mt, |nu| and EG's float L(1) consistently")
    nuEF = {"4gamma/3": (0.396533, 0.396698), "0.95": (0.240516, 0.240608)}
    okn = bool(rows[0][2]["nu_lo"] <= arb("0.396533")) and bool(rows[0][2]["nu_hi"] >= arb("0.396698"))
    okn = okn and all(bool(r[2]["nu_lo"] <= arb("0.240516")) and bool(r[2]["nu_hi"] >= arb("0.240608")) for r in rows if r[0] <= 0.95 <= r[1])
    check(okn, "consistency: cell |nu| brackets contain EF.3's nu(4gamma/3) and EG.4's nu(0.95)")
    # ---- EH.5 floats
    al_f = [4 * GAMMA / 3 + 1e-9, 0.915, 0.92, 0.925, 0.93, 0.95, 0.97, 0.99, 0.999]
    sv = survey_float(al_f)
    for s in sv:
        print(f"   floats alpha = {s['alpha']:.5f}: |nu|^3 D' = {s['T']:.5f}, f|nu|^3D' = {s['f'] * s['T']:.5f}, L(1) = {s['L1']:.5f}, "
              f"L(2) = {s['L2']:.5f}, nu = {s['nu']:.5f}, n=1 eta-route bound {s['S11']:.6f}")
    check(all(s["f"] ** 2 * s["T"] < s["L2"] for s in sv) and all(s["S11"] < 0 for s in sv),
          "EH.5 (OBSERVED): float margins positive on the grid; float n = 1 bound < 0")
    inside = all(bool(r[2]["nu_lo"] <= -s_["nu"] <= r[2]["nu_hi"]) for s_ in sv for r in rows if r[0] <= s_["alpha"] <= r[1] and "nu_lo" in r[2])
    inside = inside and all(bool(-s_["nu"] <= r[2]["nu_hi"]) for s_ in sv for r in rows if r[0] <= s_["alpha"] <= r[1])
    check(inside, "EH.5: every float nu lies inside the |nu| bracket of each cell containing its alpha")
    first = next((s["alpha"] for s in sv if s["f"] * s["T"] < s["L1"]), None)
    print(f"   floats: L(1) alone covers n = 1 from about alpha = {first}")
    RESULT["floats"] = sv
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EH checks passed")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0])
    sys.exit(main(full="--full" in sys.argv))
