#!/usr/bin/env python3
r"""Theorem DU: sharper certified diameters D_7, D_8 of the second-order cone, and DT's level-sum bound re-derived with
them (DT's open item "certified D_n on K_(2,6) for n >= 9 (or sharper D_7, D_8) to push the level-sum rate from 0.2694
toward 0.1725").

Setting (GASKET.md, Theorems DQ-DT): s in X = [sqrt5/2, 5/2], Lt^n f(s) = sum_w W_w(s) f(y_w(s)) over the 2^n words,
K = K_(2,6) = {f > 0 : |f'| <= 2f, |f''| <= 6f} (DQ.1: Lt K in K), Hilbert metric read through
Lam = {2 d_s +- d'_s, 6 d_s +- d''_s : s in X}, Phi_l(f) = l(Lt^n f), R(l, l') = sup_K Phi_l / Phi_l',
D_n = diam Lt^n K = sup_(l, l') log R(l, l') R(l', l).

Theorem DU (checked here; exit 1 on any failure):
  DU.1 [proved] Several references.  For any positive functionals r_1, ..., r_K on K (here r_k = d_(s_k) o Lt^n,
       s_k in {9/5, 47/20}) and l, l' in Lam: R(l, l') R(l', l) <= [R(l, r_k) R(r_k, l)] [R(l', r_k) R(r_k, l')] for every k,
       so with p_k(l) := log R(l, r_k) R(r_k, l):  D_n <= max_(l, l') min_k (p_k(l) + p_k(l')).  DR.2 bounds p_k uniformly
       over a cell and a functional type (convexity in tau, two chains per direction); DR used K = 1 (r = d_(9/5)) and
       the cruder max_(l, l') <= 2 max_l.  The bound is invariant under the rescaling of l (it enters R(l, r) R(r, l) only
       through a product).
  DU.2 [certified; DR.1's chain lemma and DR.2's cell lemma unchanged, chains proposed by DR's float LP and replayed in
       arb; n = 7 on 800 cells, n = 8 on 1600 cells, 4 types, 2 references]
         D_7 <= DBAR7,  D_8 <= DBAR8  (DR.3: 0.00041194852, 0.00017542239),
       so theta_8^(1/8) = tanh(D_8/4)^(1/8) <= PER8 (DR.3's best, theta_7^(1/7): 0.2694000).
  DU.3 [proved from DU.2, DR.3, DS.3, CF.1] DR.4's composition and envelope with the new factors; DT.3's level-sum bound
       re-derived: B_U(k) = min_(m <= min(5, k-4)) Theta_K(F_m, F_(m+1)) G_U(k-4-m) with DT.2's uppers as stated and the
       factors theta_n (n = 3..8; n = 7, 8 from DU.2), per-level rate theta_8^(1/8); B_U <= B_K (DT.3) for every k, with
       B_U < B_K exactly where a block 7 or 8 is used.
  DU.4 [certified lower bound, arb] explicit f, g in K (exponentials e^(+-2s)) give D_n >= Theta(Lt^n f, Lt^n g) >= LOWn
       (n = 7, 8): the uppers are not vacuous, and the gap is stated.
DU_FULL=1 runs all cells (about 4.5 h on 6 worker processes on a shared box); the default runs the cells DU_ONLY (the ends and the middle),
which re-checks the replay at the worst cells and the stated bounds on those cells only (not a certificate for D_n).
Kept limits: see GASKET.md Theorem DU.  No Psi_info, Stage 2 still blocked.  Stage 1 untouched; no W, 0.08, 0.23 or
30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

import numpy as np
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_effective_threshold as dr   # noqa: E402  (DR.1 chain, DR.2 cell lemma, LP)

FAILS: list[str] = []
RESULT: dict = {}
FULL = os.environ.get("DU_FULL", "") not in ("", "0")
NCELLS = {7: 800, 8: 1600}
REFS = ((9, 5), (47, 20))                       # reference functionals d_(9/5), d_(47/20) (composed with Lt^n)
TYPES = ("1+", "1-", "2+", "2-")
WORKERS = int(os.environ.get("DU_WORKERS", "6"))
DEFAULT_ONLY = {7: (0, 1, 400, 798, 799), 8: (0, 1, 800, 1598, 1599)}
UP = lambda x: float(np.nextafter(x, np.inf))
CAP = "0.3170757"
# DR.3 / DS.3 factors used as stated (n = 3..8); n = 7, 8 are replaced by DU.2
DBAR_DR = {3: "0.25513308", 4: "0.03838589", 5: "0.0098239159", 6: "0.0018377053", 7: "0.00041194852", 8: "0.00017542239"}
# certified targets (set after the full run; checked, exit 1 if not met)
TARGET = {7: "0.00023483", 8: "0.000041213"}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- DU.2: per cell, per type, per reference
_REFCACHE: dict = {}


def _ref(n, k):
    key = (n, k)
    if key not in _REFCACHE:
        ctx.prec = 96
        a, b = REFS[k]
        Z = dr.word_jets(arb(a) / b, n)
        y0 = [w[1][0] for w in Z]
        w0 = [w[0][0] for w in Z]
        sQ = arb(1 / float(sum(w0).mid()))
        _REFCACHE[key] = (y0, [w * sQ for w in w0])
    return _REFCACHE[key]


def cell_job(args):
    """p[k][typ] >= sup over the cell of log R(l, r_k) R(r_k, l) (DR.2 + DR.1 replay in arb); returns floats rounded up."""
    n, ncell, ci_ = args
    ctx.prec = 96
    t0 = time.time()
    AX, BX = arb(5).sqrt() / 2, arb(5) / 2
    edges = np.linspace(float(AX.mid()), float(BX.mid()), ncell + 1)
    sc = (edges[ci_] + edges[ci_ + 1]) / 2
    dl = (edges[ci_ + 1] - edges[ci_]) / 2 * (1 + 1e-9)
    cball = arb(sc, dl)
    fl = lambda L: np.array([float(v.mid()) for v in L])
    cen, ball = dr.word_jets(arb(sc), n), dr.word_jets(cball, n)
    out, bumps, nch = [], 0, 0
    for k in range(len(REFS)):
        y0, w0 = _ref(n, k)
        row = []
        for typ in TYPES:
            ys, a0, a1, b0, b1, E = dr.cell_functional(typ, cen, ball, dl)
            S = arb(1 / float(sum(a0).mid()))
            a0, a1, b0, b1, E = ([v * S for v in L] for L in (a0, a1, b0, b1, E))
            yall = ys + y0
            o = np.argsort(fl(yall), kind="stable")
            yf = fl(yall)[o]
            ya = [yall[j] for j in o]
            nP, z0 = len(ys), [arb(0)] * len(y0)
            tv = []
            for direction in (1, 2):
                tmax = 0.0
                for tau in (-dl, dl):
                    T = arb(tau)
                    G = [p + T * q for p, q in zip(a0, a1)]
                    Gb = [p + T * q for p, q in zip(b0, b1)]
                    if direction == 1:       # t Phi_r - (G + E) >= 0  =>  R(l_s, r) <= t
                        aP, bP = [-(g + e) for g, e in zip(G, E)] + z0, [-g for g in Gb] + z0
                        aQ, bQ = [arb(0)] * nP + w0, [arb(0)] * (nP + len(y0))
                    else:                    # t (G - E) - Phi_r >= 0  =>  R(r, l_s) <= t
                        aP, bP = [arb(0)] * nP + [-w for w in w0], [arb(0)] * (nP + len(y0))
                        aQ, bQ = [g - e for g, e in zip(G, E)] + z0, Gb + z0
                    aP, bP, aQ, bQ = ([L[j] for j in o] for L in (aP, bP, aQ, bQ))
                    sol = dr.lp_split(yf, fl(aP), fl(bP), fl(aQ), fl(bQ))
                    if sol is None:
                        return dict(ci=ci_, err=f"LP proposal failed ({typ}, ref {k}, dir {direction})")
                    t, ma, mua = sol
                    t = float(np.nextafter(t, np.inf))
                    for bump in range(7):
                        Tt = arb(t)
                        alpha = [Tt * q + p for q, p in zip(aQ, aP)]
                        beta = [Tt * q + p for q, p in zip(bQ, bP)]
                        if dr.verify_chain(ya, alpha, beta, ma, mua) >= 0:
                            break
                        t, bumps = t * (1 + 1e-8 * 4 ** bump), bumps + 1
                    else:
                        return dict(ci=ci_, err=f"arb replay failed ({typ}, ref {k}, dir {direction})")
                    nch += 1
                    tmax = max(tmax, t)
                tv.append(tmax)
            row.append(UP(float((arb(tv[0]) * arb(tv[1])).log().upper())))
        out.append(row)
    return dict(ci=ci_, p=out, bumps=bumps, chains=nch, lo=float(cball.lower()), hi=float(cball.upper()),
                sec=time.time() - t0)


def combine(P):
    """P: (items, K) upper bounds p_k(item) >= 0.  Returns an upper bound for max_(x, y) min_k (P[x,k] + P[y,k]) and the
    maximising pair.  Sums are rounded up by one ulp (they are >= 0, so one ulp covers the float rounding)."""
    best, arg = -1.0, None
    B = 512
    for i0 in range(0, len(P), B):
        A = P[i0:i0 + B]
        S = np.nextafter(A[:, None, :] + P[None, :, :], np.inf).min(axis=2)
        j = int(np.argmax(S))
        v = float(S.flat[j])
        if v > best:
            best, arg = v, (i0 + j // len(P), j % len(P))
    return best, arg


def certify(n, only=None, workers=WORKERS):
    ncell = NCELLS[n]
    cells = list(range(ncell)) if only is None else list(only)
    t0 = time.time()
    jobs = [(n, ncell, c) for c in cells]
    if workers > 1 and len(jobs) > 1:
        import multiprocessing as mp
        with mp.get_context("fork").Pool(workers) as pool:
            res = list(pool.imap(cell_job, jobs, chunksize=1))
    else:
        res = [cell_job(j) for j in jobs]
    errs = [r for r in res if "err" in r]
    if errs:
        check(False, f"DU.2 n={n}: {errs[0]['err']} at cell {errs[0]['ci']} ({len(errs)} failures)")
        return None
    res.sort(key=lambda r: r["ci"])
    AX, BX = arb(5).sqrt() / 2, arb(5) / 2
    if only is None:
        cover = (res[0]["lo"] <= float(AX.lower()) and res[-1]["hi"] >= float(BX.upper())
                 and all(res[i]["lo"] <= res[i - 1]["hi"] for i in range(1, len(res))))
        # float endpoints of arb balls: lo/hi are rounded to nearest; re-check with arb on the ends
        cover = cover and bool(arb(res[0]["lo"]) < AX) and bool(arb(res[-1]["hi"]) > BX)
        check(cover, f"DR.2 n={n}: the {ncell} cell balls cover X = [sqrt5/2, 5/2] (consecutive overlap, ends outside)")
    # items = (cell, type); P[item, k]
    P = np.array([[r["p"][k][t] for k in range(len(REFS))] for r in res for t in range(len(TYPES))])
    names = [(r["ci"], TYPES[t]) for r in res for t in range(len(TYPES))]
    half, harg = combine(P)
    single = UP(2 * P[:, 0].max())                       # DR's bound with r = d_(9/5) only, same cells
    Dn = arb(half).upper()
    Dn = arb(Dn)
    th = (Dn / 4).tanh()
    per = (th.log() / n).exp()
    nb = sum(r["bumps"] for r in res)
    nch = sum(r["chains"] for r in res)
    ws = (names[harg[0]], names[harg[1]])
    print(f"     DU.2 n={n}: {len(res)} of {ncell} cells, {nch} chains replayed in arb ({nb} t-bumps), worst pair {ws}, "
          f"single-reference bound 2 max p_1 = {single:.6e}, {time.time() - t0:.0f}s", flush=True)
    return dict(D=Dn, theta=th, per=per, worst=ws, single=single, P=P, names=names, cells=len(res), chains=nch, bumps=nb,
                pmax=[float(P[:, k].max()) for k in range(len(REFS))])


# ----------------------------------------------------------------- DU.4: certified lower bound (explicit f, g in K)
# f = exp(al s + be s^2): v = f'/f = al + 2 be s (linear), f''/f = 2 be + v^2.  Float proposals (Nelder-Mead on a 57-point
# grid, n = 7); only the arb evaluation below is claimed.
WIT = (("3.8637", "-0.83356"), ("-1.1360", "-0.17278"))


def wit_in_K(al, be):
    """exact interval check: |v| <= 2 and |2 be + v^2| <= 6 on X (v linear, so v^2 is convex: max at the ends, min at an
    end or at v = 0)."""
    A, B = arb(5).sqrt() / 2, arb(5) / 2
    vA, vB = al + 2 * be * A, al + 2 * be * B
    ok = bool(abs(vA) <= 2) and bool(abs(vB) <= 2)
    qs = [2 * be + vA * vA, 2 * be + vB * vB]
    if not (bool(vA > 0) and bool(vB > 0)) and not (bool(vA < 0) and bool(vB < 0)):
        qs.append(2 * be)
    return ok and all(bool(q <= 6) and bool(q >= -6) for q in qs)


def wit_ells(n, s, al, be):
    """(2F + F', 2F - F', 6F + F'', 6F - F'') at s for F = Lt^n f, f = exp(al y + be y^2) (arb)."""
    F = [arb(0)] * 3
    for W, y in dr.word_jets(s, n):
        e = (al * y[0] + be * y[0] * y[0]).exp()
        v = al + 2 * be * y[0]
        fp, fpp = v * e, (2 * be + v * v) * e
        F[0] += W[0] * e
        F[1] += W[1] * e + W[0] * fp * y[1]
        F[2] += W[2] * e + (2 * W[1] * y[1] + W[0] * y[2]) * fp + W[0] * y[1] * y[1] * fpp
    return [2 * F[0] + F[1], 2 * F[0] - F[1], 6 * F[0] + F[2], 6 * F[0] - F[2]]


def lower_witness(n, npts=57):
    ctx.prec = 96
    (a1, b1), (a2, b2) = ((arb(a), arb(b)) for a, b in WIT)
    inK = wit_in_K(a1, b1) and wit_in_K(a2, b2)
    A, B = arb(5).sqrt() / 2, arb(5) / 2
    pts = [A + (B - A) * arb(i) / (npts - 1) for i in range(npts)]
    rho = []
    for s in pts:
        lf, lg = wit_ells(n, s, a1, b1), wit_ells(n, s, a2, b2)
        for t, (x, y) in enumerate(zip(lf, lg)):
            if not (bool(x > 0) and bool(y > 0)):
                return None
            rho.append((x / y, t))
    hi = max(rho, key=lambda z: z[0].mid())
    lo = min(rho, key=lambda z: z[0].mid())
    L = (hi[0] / lo[0]).log()            # = log [l(F) l'(G)] / [l'(F) l(G)] <= Theta(Lt^n f, Lt^n g) <= D_n
    return dict(inK=inK, low=arb(L.lower()))


# ----------------------------------------------------------------- DU.3: composition, envelope, DT.3 re-derived
def consequences(dbar):
    """dbar: {7: str, 8: str} certified (DU.2) upper bounds; DR.3/DS.3's n = 3..6 as stated."""
    import gasket_fractional_profile_level_second_order_sums as dt
    ctx.prec = dt.PREC
    new = dict(DBAR_DR)
    new.update({n: dbar[n] for n in (7, 8)})
    Dq = {n: arb(new[n]) for n in new}
    th = {n: arb((Dq[n] / 4).tanh().upper()) for n in Dq}
    # composition (CF.1): diam Lt^(a+b) K <= tanh(D_a/4) D_b
    for j in range(9, 17):
        Dq[j] = min((arb((th[a] * Dq[j - a]).upper()) for a in range(3, 9) if j - a in Dq), key=lambda x: x.upper())
    rho = (th[8].log() / 8).exp()
    rho = arb(rho.upper())
    C = arb(max(((Dq[j] / 4) / rho ** j).upper() for j in range(5, 13)))
    per = {n: arb(((th[n].log()) / n).exp().upper()) for n in th}
    unif = arb(max(((Dq[j] / 4).log() / j).exp().upper() for j in range(5, 13)))
    check(bool(rho < arb(DBAR_PER7)), f"DU.3: per-level rate theta_8^(1/8) <= {rho.str(8, radius=False)} < DR.3's theta_7^(1/7) "
          f"<= {DBAR_PER7} (theta_7^(1/7) now <= {per[7].str(8, radius=False)}); envelope tanh(D_n/4) <= C rho^n for every "
          f"n >= 5 (n = 8q + j, j = 5..12, D_n <= theta_8^q D_j, D_j for j = 9..12 by CF.1 splits), rho = theta_8^(1/8), "
          f"C <= {C.str(6, radius=False)}; tanh(D_n/4)^(1/n) <= {unif.str(7, radius=False)} for every n >= 5 (DR.4: 0.3135158)")
    RESULT.update(th=th, per=per, rho=rho, C=C, unif=unif, Dq=Dq)
    # DT.3 re-derived: B with the DT.2 uppers as stated, old (DT) vs new (DU) factors
    up = dt.TARGETS[1024]["UPPER"]
    dt.RESULT.clear()
    dt.RESULT["tab"] = {m: {"up": arb(up[m]), "low": arb(dt.LOWER[m])} for m in range(dt.MMAX + 1)}
    KM = 600
    th_old = {n: arb((arb(D) / 4).tanh().upper()) for n, D in DBAR_DR.items()}
    G_old = dt.knap([(n, th_old[n]) for n in sorted(th_old)], KM)
    G_new = dt.knap([(n, th[n]) for n in range(3, 9)], KM)
    ks = list(range(4, KM - 10))
    Bo = {k: dt.B_K(k, G_old) for k in ks}
    Bn = {k: dt.B_K(k, G_new) for k in ks}
    le = all(bool(Bn[k].upper() <= Bo[k].upper()) for k in ks)
    same = [k for k in ks if not bool(Bn[k].upper() < Bo[k].lower())]
    first_strict = min(k for k in ks if bool(Bn[k].upper() < Bo[k].lower()))
    fac = {k: Bo[k] / Bn[k] for k in ks}
    show = (11, 12, 13, 15, 22, 49, 103, 200)
    # tail: B_U(k) <= Theta(F_5, F_6) theta_8^floor((k-9)/8); DT's B_K(k) >= c_K rho_K^(k-4) with rho_K = min DT per level
    eqset = [k for k in same if k >= 11]
    kfrom = max(same) + 1
    check(le and first_strict == 11 and kfrom <= 20,
          f"DU.3: B_U(k) <= B_K(k) for every k = 4..{KM - 11} (G_U <= G_K blockwise); B_U(k) < B_K(k) (arb-separated) for "
          f"k = 11, 12 and every k >= {kfrom}; equal for k <= 10 (no block of length 7 or 8 fits; DT.4 unchanged) and for "
          f"k = {', '.join(map(str, eqset))} (DT's optimum there uses blocks 3..6 only); B_K/B_U: "
          + ", ".join(f"{fac[k].str(4)} (k={k})" for k in show)
          + "; B_U(k): " + ", ".join(f"{Bn[k].str(5)} (k={k})" for k in show))
    RESULT.update(Bo={k: Bo[k] for k in show}, Bn={k: Bn[k] for k in show}, fac={k: fac[k] for k in show},
                  first_strict=first_strict, eqset=eqset, kfrom=kfrom)
    # DT's own checks (B_U < B_DO for k >= kstar, consistency with CC's s_4..s_13, DT.4 floors) re-run with DU's factors
    print("     DU.3: DT.3/DT.4's checks re-run with DU.2's factors (DBAR_7, DBAR_8 replaced):", flush=True)
    dt.FAILS.clear()
    saved = dict(dt.DBAR)
    dt.DBAR.clear()
    dt.DBAR.update(new)
    try:
        dt.consequences()
    finally:
        dt.DBAR.clear()
        dt.DBAR.update(saved)
    check(not dt.FAILS and dt.RESULT.get("kstar") == 11,
          f"DU.3: with DU.2's factors B_U(k) < B_DO(k) for every k >= {dt.RESULT.get('kstar')} (as DT.3: 11), consistency with "
          f"CC's s_4..s_13 (max ratio {dt.RESULT['worst'].str(4) if 'worst' in dt.RESULT else '?'}), DT.4 (k = 4..10) unchanged")
    RESULT.update(kstar=dt.RESULT.get("kstar"), worst_cc=dt.RESULT.get("worst"))


DBAR_PER7 = "0.2694000"


def main() -> int:
    t0 = time.time()
    print(f"Theorem DU: sharper certified D_7, D_8 on K_(2,6) ({'FULL: all cells' if FULL else 'default: cells ' + str(DEFAULT_ONLY)}). "
          "Not thrust.", flush=True)
    ctx.prec = 96
    res = {}
    for n in (7, 8):
        r = certify(n, only=None if FULL else DEFAULT_ONLY[n])
        if r is None:
            break
        res[n] = r
        tgt = arb(TARGET[n])
        check(bool(r["D"] <= tgt),
              f"DU.2 n={n}: {'D_n' if FULL else 'the bound over the default cells'} <= max_(l,l') min_k (p_k(l) + p_k(l')) <= "
              f"{r['D'].str(8, radius=False)} <= {TARGET[n]}{'' if FULL else ' (stated D_n; DU_FULL=1 certifies it)'}; "
              f"per level <= {r['per'].str(7, radius=False)}; worst pair {r['worst']}; max_l p_k(l) = "
              + ", ".join(f"{x:.5e}" for x in r["pmax"]) + f" (k = d_(9/5), d_(47/20)); single reference: {r['single']:.5e}; "
              f"DR.3: {DBAR_DR[n]}")
    RESULT["cells"] = res
    if len(res) == 2:
        for n in (7, 8):
            w = lower_witness(n)
            ok = w is not None and w["inK"] and bool(w["low"] < arb(TARGET[n]))
            check(ok, f"DU.4 n={n}: f = exp(al s + be s^2), (al, be) = {WIT[0]} and {WIT[1]}, both in K_(2,6) (exact interval "
                      f"check); D_{n} >= Theta(Lt^n f, Lt^n g) >= {w['low'].str(6, radius=False) if w else '?'} "
                      f"(per level >= {(((w['low'] / 4).tanh().log()) / n).exp().str(6, radius=False) if w else '?'}); gap to "
                      f"DU.2: factor {(arb(TARGET[n]) / w['low']).str(4, radius=False) if w else '?'}")
            RESULT[f"low{n}"] = w
        consequences(TARGET)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem DU checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
