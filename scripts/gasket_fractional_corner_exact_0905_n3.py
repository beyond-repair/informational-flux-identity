#!/usr/bin/env python3
"""Theorem EL: the exact sign-change level at alpha = 0.905 for N = 3.  EK kept "N = 3 at 0.905 still width 2 (103 undecided by this
method)": the last level with Y_n < Y_oo was located in [102, 103]; EJ.2's lower side at n = 103 gave +0.005534 (not < 0).  Floats
(EE.1's limit integrals) put the true braces at n = 103 at about +0.00125, so the missing certificate is on the UPPER side, and EJ.1's
upper integrals on EK's grid (kappa_hi(3, 14) <= 1.734107) are too wide by about 0.012.

  EL.1 (THEOREM; every alpha in (gamma, 1), k0 >= 1, 1 <= N <= k0; EJ.1 and AQ (3) as stated).  On every cell [a, b] (0 < a < b) and
       every k >= k0:  sup_[a,b] R_k <= R^hi(a) = min(1/3, G a^(gamma-1), 1/V^lo(a))  and  sup_[a,b] Rh_k <= (5/3) R^hi_(k0+1)(5a),
       where V^lo(a) is EJ.1's V^lo evaluated at the single point a (no mixing of the cell's two ends).  Proof: R_k(x) = (3/5)^k
       r_k(5^-k x) is nonincreasing in x (AQ (3): r_k decreasing in s), so sup_[a,b] R_k = R_k(a) <= R^hi(a) (EJ.1 is pointwise);
       likewise Rh_k(x) = (5/3) R_(k+1)(5x).  Hence EJ.1's cell sums for c_hi(k0), kappa_hi(N, k0) stay upper bounds when R, Rh are
       evaluated at the left node.  (EJ's cell routine bounds V^lo on [a, b] by A(b)(3 + phi_k0(a)) - sum e(z_m(b)) ...: the two
       large terms are taken at opposite ends, which leaves an O(1)-in-n excess where Phi is large -- OBSERVED: about 0.75 % of R at
       x = 5^4 even with 3200 cells per unit, and kappa_hi(3, 26) decreasing only like n^(-1/2): 1.73938, 1.73453, 1.73113 for n =
       1000, 2000, 4000.  With EL.1 the excess is O(1/n).)  Nothing else in EJ.1 changes; head, tail and t >= 1 pieces are EJ's.
  EL.2 (arb, 128 bits; --full; EF.1 point cells n_p = 4000, n_m = 32000; EL.1 upper k0 = 26 on n = 3200 cells per unit): at alpha =
       0.905, EJ.2's upper criterion holds at n1 = 103 for N = 3 (Phi(103) > 0, Phi increasing in n), so Y_n > Y_oo for every n >= 103;
       with EK.1 (Y_n < Y_oo for n = 3..102) Y_n - Y_oo changes sign exactly once, between levels 102 and 103.  N = 1, 2: onsets 68, 84
       (reproduce EJ.3 / EK.1).
  EL.3 (arb; --full; same grids): at alpha = 0.91 (N = 1), EJ.2's upper criterion holds at n1 = 2587, so Y_n > Y_oo for every
       n >= 2587 (EJ.4: >= 2616); with EJ.4's window (Y_n < Y_oo up to 2564) the last negative level is in [2564, 2586] (EJ: [2564,
       2615]).  The float limit crossing (EE.1's e*, kappa_1) is at 2586 with braces ~1e-5: not closable on these grids.
  Default (no --full): EL.2 only, on n = 1600 cells per unit with point cells 2000/16000 (it already closes: Phi(103) >= 2.35e-4).

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_two_sided_rate as ej  # noqa: E402

FAILS: list[str] = []
A_S = "0.905"
LO, HI = -25, 32                 # EJ.1's grid range (head x < 5^LO, tail x > 5^HI as in EJ)
K0 = 26
FULL = {"n_p": 4000, "n_m": 32000, "n": 3200}
DEFAULT = {"n_p": 2000, "n_m": 16000, "n": 1600}
EK_LAST = {1: 67, 2: 83, 3: 102}          # EK.1: Y_n < Y_oo for n = N..EK_LAST[N] (certified)
EK_BRACKET3 = (102, 103)                  # EK.1: last negative level at N = 3
EK_POINT = {"Dp": (5.86416, 5.86461), "nu": (0.41474, 0.415168)}   # EK.1's point brackets (consistency only)
EK_UP = {"U": 0.740363, "kappa3": 1.734107}                        # EK.1's U(14), kappa_hi(3, 14)
EK_LO = {"L14": 0.736435, "kappa3": 1.70926}                       # EK.1's L(14), kappa_lo(3, 14) (consistency only)
EXPECT = {1: 68, 2: 84, 3: 103}
A_91, EJ_091, EXPECT_091 = "0.91", (2564, 2616), 2587   # EJ.4 at 0.91 (N = 1): Y_n < Y_oo for n <= 2564 (certified), > for n >= 2616
WORKERS = int(os.environ.get("EL_WORKERS", "4"))


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def _up(v) -> float:
    """A float >= the arb v (round the upper endpoint to nearest, then one ulp up)."""
    return math.nextafter(float(v.upper()), math.inf)


def _cells(args):
    """EJ.1's main cell sums (c, kappa_N) over grid cells i0 <= i < i1, with EL.1's left-node R^hi.  Returns float upper bounds."""
    a_s, k0, n, i0, i1, Ns = args
    arb, ctx = ej._arb()
    ctx.prec = 128
    al = arb(a_s)
    g = arb(3).log() / arb(5).log()
    five, third, G, T = arb(5), arb(1) / 3, arb("0.6031"), arb(1)
    f35, i3 = (arb(3) / 5) ** k0, arb(3) ** (-k0)
    M = HI + 42
    c = arb(0)
    kap = {N: arb(0) for N in Ns if N <= k0}
    nneg = 0
    for i in range(i0, i1):
        a, b = arb(float(5.0 ** (LO + i / n))), arb(float(5.0 ** (LO + (i + 1) / n)))
        pa, pa5 = ej.node_phi(a, M), ej.node_phi(5 * a, M + 1)
        ka, ka5 = ej.node_k0(a, k0), ej.node_k0(5 * a, k0 + 1)
        V = ej.vlo_cell(pa, pa, ka, k0)                    # EL.1: V^lo at the single point a
        R = ej._min([third, G * a ** (g - 1)] + ([1 / V] if bool(V > 0) else []))
        V5 = ej.vlo_cell(pa5, pa5, ka5, k0 + 1)
        R5 = ej._min([third, G * (5 * a) ** (g - 1)] + ([1 / V5] if bool(V5 > 0) else []))
        nneg += 0 if bool(V > 0) else 1
        Rh = arb((five / 3 * R5).upper())
        E = (arb("0.18") * ej._min([T, b / five**k0])).exp()
        u = ej._min([third, i3 * b / 3 * R])
        v = ej._min([third, i3 * b / 3 * Rh])
        sk = arb(ka["s"][k0].lower())
        W0 = (b ** (1 - al) - a ** (1 - al)) / (1 - al)
        Dl = R * Rh * ka["Q"][0] * (1 + 1 / (sk + 5))
        brk = arb("0.18") * b * f35 * E + E * b / 3 * (R + Rh) / ((1 - u) * (1 - v))
        c += W0 * Dl * brk
        for N in kap:
            eh = ej._min([arb(1), ej.eta_fin(N, ka["s"], k0)])
            kap[N] += W0 * R * (arb("0.18") * b * f35 * E + b / 3 * R) / (1 - u) * eh
    return _up(c), {N: _up(v) for N, v in kap.items()}, nneg


def certify_upper_el(a_s, k0, n, Ns=ej.NS, workers=WORKERS):
    """EJ.1's c_hi(k0), U(k0), kappa_hi(N, k0) with EL.1's cell bound (main sums in parallel chunks; head, tail, t >= 1 as in EJ)."""
    arb, _ = ej._arb()
    al = arb(a_s)
    g = arb(3).log() / arb(5).log()
    five, third, G, T = arb(5), arb(1) / 3, arb("0.6031"), arb(1)
    rho = five**al / 3
    f35 = (arb(3) / 5) ** k0
    ncell = (HI - LO) * n
    step = max(1, ncell // (8 * workers))
    jobs = [(a_s, k0, n, i, min(i + step, ncell), tuple(Ns)) for i in range(0, ncell, step)]
    with ProcessPoolExecutor(workers) as ex:
        res = list(ex.map(_cells, jobs))
    c = sum((arb(r[0]) for r in res), arb(0))
    kap = {N: sum((arb(r[1][N]) for r in res), arb(0)) for N in res[0][1]}
    nneg = sum(r[2] for r in res)
    X0, X1 = arb(float(5.0 ** LO)), arb(float(5.0 ** HI))
    E1, ET = (arb("0.18") * X0).exp(), (arb("0.18") * T).exp()
    head = arb(2) / 9 * E1 * (arb("0.18") + third * (third + arb(5) / 9) * arb(9) / 4) * X0 ** (2 - al) / (2 - al)
    tail = G**2 * 30 * ET * (arb("0.18") * f35 * X1 ** (2 * g - 1 - al) / (1 + al - 2 * g)
                             + G * arb(3) / 2 * X1 ** (3 * g - 2 - al) / (2 + al - 3 * g))
    pT = [T]
    for _ in range(k0 + 1):
        pT.append(pT[-1] * (pT[-1] + 5))
    Pb = arb(1)
    for jj in range(k0):
        Pb = Pb * ej._be(pT[jj])
    cond = bool((3 / rho) * ej._be(pT[k0]) <= 1)
    Iw = arb("1.5") * (T ** (-al - 1) / (al + 1) + T ** (-al - arb(1) / 2) / (5 * (al + arb(1) / 2)))
    ttail = (3 / rho) ** k0 * arb(27) / 10 * Pb * Iw
    chi = arb((c + head + tail + ttail).upper())
    row = {"c_main": c, "c_hi": chi, "U": arb((3 * chi / (2 * (1 - rho / 3))).upper()), "cond": cond, "nneg": nneg, "kappa": {}}
    for N in kap:
        hk = third * (arb("0.18") * E1 + arb(1) / 9) * arb(3) / 2 * X0 ** (2 - al) / (2 - al)
        tk = 3 * five**N * G * (arb("0.18") * f35 * ET * X1 ** (g - al) / (al - g) + G / 3 * X1 ** (2 * g - 1 - al) / (1 + al - 2 * g))
        Ft = arb("2.25") * (T ** (-al) / al + T ** (arb(1) / 2 - al) / (5 * (al - arb(1) / 2)))
        ttk = (3 / rho) ** k0 * 2 / pT[k0 - N] * Ft
        row["kappa"][N] = {"main": kap[N], "rest": hk + tk + ttk, "hi": arb((kap[N] + hk + tk + ttk).upper())}
    return {k0: row}


def main(full: bool = False) -> int:
    t0 = time.time()
    cfg = FULL if full else DEFAULT
    print(f"Theorem EL: exact sign-change level at alpha = {A_S}, N = 3 (EJ.1 with EL.1's left-node cell bound). Not thrust. "
          f"({'--full' if full else 'default'}: {cfg}, k0 = {K0}, {WORKERS} workers)", flush=True)
    arb, ctx = ej._arb()
    ctx.prec = 128
    P = ej.certify_point(A_S, cfg["n_p"], cfg["n_m"])
    print(f"   point: D'_oo in {ej.fmt(P['Dp'])}, |nu| in {ej.fmt(-P['nu'], 7)}, 1-delta_N in "
          + ", ".join(f"N={N}: {ej.fmt(P['om'][N])}" for N in ej.NS) + f" ({time.time() - t0:.1f}s)", flush=True)
    check(P["ok"], f"EF.1 routines valid at alpha = {A_S}")
    nu_lo, nu_hi = float((-P["nu"]).lower()), float((-P["nu"]).upper())
    check(EK_POINT["Dp"][0] - 1e-5 <= float(P["Dp"].lower()) and float(P["Dp"].upper()) <= EK_POINT["Dp"][1] + 1e-5
          and EK_POINT["nu"][0] - 1e-5 <= nu_lo and nu_hi <= EK_POINT["nu"][1] + 1e-5,
          "consistency: D'_oo, |nu| brackets overlap-consistent with EK.1's (within 1e-5)")
    t1 = time.time()
    up = certify_upper_el(A_S, K0, cfg["n"])
    row = up[K0]
    print(f"   upper (EL.1, k0 = {K0}, n = {cfg['n']}): c_hi <= {float(row['c_hi'].upper()):.6f}, U({K0}) <= {float(row['U'].upper()):.6f}, "
          + ", ".join(f"kappa_hi({N}) <= {float(row['kappa'][N]['hi'].upper()):.6f}" for N in ej.NS)
          + f" (V^lo <= 0 on {row['nneg']} cells; {time.time() - t1:.1f}s)", flush=True)
    check(row["cond"], f"EJ.1 alpha = {A_S}: (3/rho) beta(p^k0(1)) <= 1")
    check(bool(row["U"] > EK_LO["L14"]) and bool(row["kappa"][3]["hi"] > EK_LO["kappa3"]),
          "two-sided: U(26) > EK's L(14) >= 0.736435 and kappa_hi(3, 26) > EK's kappa_lo(3, 14) >= 1.70926 (no contradiction)")
    print(f"   (comparison) EK.1: U(14) <= {EK_UP['U']}, kappa_hi(3, 14) <= {EK_UP['kappa3']}")
    for N in ej.NS:
        on = ej.onset_upper(P, up, N)
        n1 = on[0] if on else None
        prev = on[3] if on else None
        print(f"     N = {N}: Phi({n1}) >= {float(on[2].lower()):.4g} > 0" + (f", Phi({n1 - 1}) <= {float(prev.upper()):.4g}" if prev is not None else "")
              + f" -> Y_n > Y_oo for EVERY n >= {n1}" if on else f"     N = {N}: no onset found")
        check(n1 == EXPECT[N], f"EL.2 alpha = {A_S}, N = {N}: Y_n > Y_oo for every n >= {EXPECT[N]}")
        check(n1 is not None and n1 - 1 <= EK_LAST[N] + (1 if N == 3 else 0),
              f"consistency with EK.1: onset {n1} - 1 <= EK's last certified negative level {EK_LAST[N]}" + (" + 1" if N == 3 else ""))
        if N == 3 and n1 == 103:
            print("     with EK.1 (Y_n < Y_oo for n = 3..102): exactly one sign change at N = 3, between levels 102 and 103")
    if full:
        t1 = time.time()
        P91 = ej.certify_point(A_91, cfg["n_p"], cfg["n_m"])
        check(P91["ok"], f"EF.1 routines valid at alpha = {A_91}")
        up91 = certify_upper_el(A_91, K0, cfg["n"], Ns=(1,))
        r91 = up91[K0]
        check(r91["cond"], f"EJ.1 alpha = {A_91}: (3/rho) beta(p^k0(1)) <= 1")
        on = ej.onset_upper(P91, up91, 1)
        n1 = on[0] if on else None
        print(f"   EL.3 alpha = {A_91} (N = 1): |nu| in {ej.fmt(-P91['nu'], 7)}, U({K0}) <= {float(r91['U'].upper()):.6f}, kappa_hi(1) <= "
              f"{float(r91['kappa'][1]['hi'].upper()):.6f}; Phi({n1}) >= {float(on[2].lower()):.3g} > 0 -> Y_n > Y_oo for EVERY n >= {n1}; "
              f"with EJ.4 (Y_n < Y_oo up to {EJ_091[0]}) the last negative level is in [{EJ_091[0]}, {n1 - 1}] (EJ: [{EJ_091[0]}, {EJ_091[1] - 1}]) "
              f"({time.time() - t1:.1f}s)", flush=True)
        check(n1 == EXPECT_091, f"EL.3 alpha = {A_91}, N = 1: Y_n > Y_oo for every n >= {EXPECT_091} (EJ.4: >= {EJ_091[1]})")
        check(n1 is not None and EJ_091[0] < n1 <= EJ_091[1], f"consistency with EJ.4: {EJ_091[0]} < onset <= {EJ_091[1]}")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EL checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(full="--full" in sys.argv))
