#!/usr/bin/env python3
"""Theorem DB: the exact rate of inf R_n^(k) -> rho(L_k) for k = 13, ..., 28 (DA's kept limit), so that, with CT, CU, CV
and DA, CT.5's two-sided bound with explicit Gamma_+-^(k) and n_0(k) holds for EVERY integer k >= 1.

Setting exactly as in Theorems CT-CV (GASKET.md): s = sqrt(t + 25/4), X = [sqrt5/2, 5/2], sigma_(1,2) = sqrt(15/4 -+ s),
(L_k f)(s) = u_1 f(sigma_1 s) + u_2 f(sigma_2 s), u_i = wt_i/(2 sigma_i)^k, h_n = L_k^n 1, R_n^(k) = h_(n+1)/h_n on X,
q^(k) = nu_2^(k)/rho(L_k); L_k acts on CG's H^2(D_CG), D_CG = {|s - 83/50| < 8/5}.  CT/CU/CV proved CT.5's bound for k <= 12,
DA for every k >= 29; CY.5/CZ.6 certified only the eigenvalues (CU.1-CU.3) for k = 13..29, not Gamma_+- (CU.4) or the
rate (CU.5).  DA's Neumann bound delta_k exceeds Gamma_+ for k <= 28.

Method.  CU's script is imported and run UNCHANGED per k (CU.1-CU.5, every one of its checks), applied to the rescaled
operator Lh_k = 3^(k+1) L_k.  Rescaling by the exact integer 3^(k+1) changes no eigenvector, Riesz projection, B/A, q,
Gamma_+- or theta_n (every quantity CU.5 compares is homogeneous of degree 0), and rho(Lh_k) - inf_X Rh_n = 3^(k+1)(rho(L_k)
- inf_X R_n).  The rescaling keeps every number O(1); the 'underflow' of CZ's kept limit is in fact interval wrapping: in
X <- A X the entry radii grow by ~35 per step while ||A^j Q|| falls by ~0.11, so at 512 bits (Krawczyk radius 1e-100) the
balls of A^j Q contain 0 after ~35 steps and _fro returns nan.  DB therefore runs CU's routine at 1024 bits (Krawczyk
radius 1e-130), which keeps the Neumann power 40 finite for every k <= 28.  The rescaling is done by wrapping CU's three
weight evaluators (disk sup bounds U_i, the matrix M_N, the pointwise L_k p on X) with the factor 3^(k+1) and CU's outer
circle |z| = 6 by |z| = 6 * 3^(k+1); CU's float observation (which uses unscaled float weights) is replaced by a no-op.  Parameters per k: N(k) = 60 + 10k for k = 14..19,
max(180, 30 + 10k) otherwise (so CU.4's projection error stays < 1e-12, CU's own threshold, which is NOT changed; at
N = 30 + 10k the k = 14..19 runs failed CU.4 with errors 5.7e-14..6.5e-10, because K_T^2 grows faster in k than eps_N
falls; with 30 more modes the largest error over k = 13..28 is 9.87e-13 at k = 20), Neumann power 40, rho_in = 0.6 |nu_2| (rescaled, rational),
eigenvalue circles of radius 1e-15 (rescaled), projection circles as in CY.5, theta_n scanned to n = 400.

Theorem DB (checked here; exit 1 on any failure).
  DB.1-DB.3 [CU.1-CU.3's proofs] for each k = 13..28: Lh_k on H^2(D_CG) has exactly two eigenvalues in |z| > rho_in, real
       and algebraically simple, 3^(k+1) rho(L_k) and 3^(k+1) nu_2^(k); rho_H = rho(L_k) on C(X).  Values agree with
       CY.5/CZ.6 (consistency).
  DB.4 [CU.4's proof] h_n = rho^n A + nu_2^n B + E_n, an arb cover of X encloses Gamma_+ = max_X B/A and Gamma_- = max_X(-B/A)
       (both > 0), min_X A > 0.
  DB.5 [CU.5's proof] theta_n decreases and theta_n < min Gamma_+- from n_0(k); for every n >= n_0(k)
           rho(1-q)|q|^n (Gamma_par - theta_n) <= rho(L_k) - inf_X R_n^(k) <= rho(1-q)|q|^n (Gamma_par + theta_n)
       (+ even n, - odd n): exact rate |q^(k)| for k = 13..28.
  DB.6 [proved; finite interval comparisons, DA.1's kappa_-+ used as stated] for k = 13..28,
       Gamma_+^(k)/((2k+3)/12 kappa_-) and Gamma_-^(k)/((2k+3)/12 kappa_+) are enclosed; both lie in (0.98, 1.02) for every
       k >= K_2pc = 19 (DA.4's leading term is already within 2% where DA's own bound is void).
  DB.7 [proved; with CT.4/CU.4/CV.4's Gamma_+- used as stated] Gamma_+^(k) and Gamma_-^(k) strictly increase for k = 1..28.
Corollary (CT, CU, CV, DB, DA): for every integer k >= 1 there is an explicit n_0(k) such that CT.5's two-sided bound holds
for all n >= n_0(k); inf R_n^(k) -> rho(L_k) at the exact rate |q^(k)| for every k >= 1.
Kept limits: n < n_0(k) not covered (n_0 grows with k, 48..100 here); Gamma_+- for k = 13..28 have no closed form (only
enclosures); DB.6 is a finite comparison, not a rate for Gamma_+- - ((2k+3)/12) kappa_-+; rho_in = 0.6|nu_2| is a
parameter choice, not optimal.
Assumptions used as stated: everything CU uses (CU's script is imported and run unchanged); DA.1's kappa_- (to 1e-14);
CT.4/CU.4/CV.4's Gamma_+- for k <= 12 (DB.7 only); CY.5/CZ.6's q^(k) (consistency only).
Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import multiprocessing as mp
import os
import sys
import time
from fractions import Fraction as Fr

from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

FAILS: list[str] = []
RESULT: dict = {}
KS = [int(x) for x in os.environ.get("DB_KS", ",".join(str(k) for k in range(13, 29))).split(",") if x]
NPROC = int(os.environ.get("DB_PROCS", str(min(8, os.cpu_count() or 1))))
M_POW_DB = 40
RIN_FRAC = 0.6
G_RAD_DB = Fr(1, 10 ** 15)
PREC_DB = int(os.environ.get("DB_PREC", "1024"))
NMAX_DB = 400
# DA.1 (used as stated): kappa_- = -kappa(sqrt5/2), kappa_+ = 3 kappa_-
KAP_M = ("0.370860305428934", "1e-14")
# CT.4, CU.4, CV.4 (used as stated, DB.7 only): Gamma_+, Gamma_- for k = 1..12 (lower, upper)
G_PRIOR = {
    1: (("0.2000641183", "0.2000642055"), ("0.5249039315", "0.5249041376")),
    2: (("0.2546353828", "0.2546354891"), ("0.6796088690", "0.6796091443")),
    3: (("0.3087742424", "0.3087743752"), ("0.8310492070", "0.8310495724")),
    4: (("0.3625750669", "0.3625752290"), ("0.9801674176", "0.9801678993")),
}
# CV.4 gives Gamma_+- to within 6e-6: (mid, radius 6e-6)
for _k, _gp, _gm in ((5, "0.416558", "1.129942"), (6, "0.471227", "1.283207"), (7, "0.526963", "1.442049"),
                     (8, "0.583993", "1.607694"), (9, "0.642394", "1.78057"), (10, "0.702119", "1.96047"),
                     (11, "0.763030", "2.14668"), (12, "0.824925", "2.33819")):
    G_PRIOR[_k] = tuple((str(float(g) - 6e-6), str(float(g) + 6e-6)) for g in (_gp, _gm))
# CY.5 / CZ.6 (used as stated, consistency only): q^(k) to the printed digits
Q_PRIOR = {13: "-0.349942106203", 18: "-0.3408023654", 23: "-0.33622604", 24: "-0.3357008", 25: "-0.3352657",
           26: "-0.3349069", 27: "-0.3346119", 28: "-0.3343701"}
Q_TOL = {13: "2e-12", 18: "2e-10", 23: "1.7e-8", 24: "2e-7", 25: "2e-7", 26: "2e-7", 27: "2e-7", 28: "2e-7"}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


N_OVERRIDE = {int(a): int(b) for a, b in (x.split(":") for x in os.environ.get("DB_N", "").split(",") if x)}


def n_of(k: int) -> int:
    # k = 14..19 need 30 more modes than 30 + 10k: at N = 30 + 10k their CU.4 projection error was 5.7e-14..6.5e-10
    # (> CU's 1e-12 threshold; K_T^2 grows faster in k than eps_N falls).  CU's threshold is unchanged.
    if k in N_OVERRIDE:
        return N_OVERRIDE[k]
    return 60 + 10 * k if 14 <= k <= 19 else max(180, 30 + 10 * k)


def pack(x: arb):
    m, e = x.mid().man_exp()
    r, f = x.rad().man_exp() if x.rad() != 0 else (0, 0)
    return (int(m), int(e), int(r), int(f))


def unpack(t) -> arb:
    m, e, r, f = t
    two = arb(2)
    return arb(m) * two ** e + arb(0, (arb(r) * two ** f).upper() if r else 0)


def run_one(k: int) -> dict:
    """CU.1-CU.5 for Lh_k = 3^(k+1) L_k (CU's script imported unchanged; weights wrapped by the exact factor 3^(k+1))."""
    import numpy as np
    import gasket_fractional_profile_level_infrate_k as cu
    S = 3 ** (k + 1)
    cu.PREC = PREC_DB
    ctx.prec = PREC_DB
    sa = arb(S)
    bm, dc, le = cu.build_matrix, cu.disk_constants, cu.l1_eval

    def build_matrix(n, prec=cu.PREC):
        M = bm(n, prec)
        return M * sa

    def disk_constants(K=cu.K_ARC):
        q1, q2, U1, U2 = dc(K)
        print(f"   (rescaled: U_i multiplied by 3^(k+1) = {S})", flush=True)
        return q1, q2, U1 * sa, U2 * sa

    def l1_eval(p, y):
        return le(p, y) * sa

    th0 = cu.theta

    def theta(*a):
        # CU's theta_n returns None while eps_n >= 1 (n far below n_0); CU.5 then prints windows at n = 20, 21, 40, 41,
        # which for n_0(k) > 41 would be None.  Report +inf there instead (no claim is made below n_0).
        t, e = th0(*a)
        return (arb("inf") if t is None else t), e

    cu.build_matrix, cu.disk_constants, cu.l1_eval, cu.theta = build_matrix, disk_constants, l1_eval, theta
    cu.float_observation = lambda *a, **kw: None
    cu.KD = k
    t0 = time.time()
    M0 = build_matrix(120, prec=128)
    Mf = np.array([[float(M0[i, j].mid()) for j in range(M0.ncols())] for i in range(M0.nrows())])
    ev = sorted(np.linalg.eigvals(Mf), key=lambda z: -abs(z))
    r, nu = repr(float(ev[0].real)), repr(float(ev[1].real))          # float seeds; Krawczyk certifies
    rf, nf = float(r), abs(float(nu))
    rin = Fr(1, round(1 / (RIN_FRAC * nf)))
    cu.N, cu.cr.M_POW, cu.NMAX = n_of(k), M_POW_DB, NMAX_DB
    cu.G_RAD = [G_RAD_DB, G_RAD_DB]
    cu.R_OUT = Fr(6 * S)              # CU's outer circle |z| = 6, rescaled
    cu.CONF[k] = dict(seeds=[r, nu], rho_in=rin, proj=[Fr(1, round(3 / rf)), Fr(1, round(2 / (nf - float(rin))))],
                      cl3=None, qcap="0.36")
    rc = cu.run_k(k)
    R = dict(cu.RESULT)
    out = dict(k=k, rc=rc, fails=list(cu.FAILS) + list(cu.cr.FAILS), N=cu.N, rin=(rin.numerator, rin.denominator),
               secs=time.time() - t0, n0=R.get("n0"))
    if rc == 0 and R.get("n0") is not None:
        ctx.prec = 512
        for key in ("rho", "nu", "q", "Gp", "Gm", "Amin", "Bmax", "tau1", "cp", "cm", "osc", "theta_n0"):
            out[key] = pack(R[key])
        out["ev"] = (pack(R["ev"][0]), pack(R["ev"][1]), R["ev"][2])
        out["od"] = (pack(R["od"][0]), pack(R["od"][1]), R["od"][2])
        out["errs"] = [pack(e) for e in R["errs"]]
    return out


def _worker(k: int) -> dict:
    log = os.environ.get("DB_LOGDIR")
    if log:
        sys.stdout = open(os.path.join(log, f"db_k{k}.log"), "w", buffering=1)
    try:
        return run_one(k)
    except Exception as exc:      # reported as a failure
        import traceback
        traceback.print_exc()
        return dict(k=k, rc=1, fails=[f"exception {exc!r}"], n0=None)


def main() -> int:
    print("Theorem DB: the exact rate of inf R_n^(k) -> rho(L_k) for k = 13..28 (DA's kept limit); with CT, CU, CV, DA "
          "for every k >= 1. Not thrust.", flush=True)
    t0 = time.time()
    FAILS.clear()
    RESULT.clear()
    if NPROC > 1 and len(KS) > 1:
        with mp.get_context("fork").Pool(min(NPROC, len(KS)), maxtasksperchild=1) as pool:
            outs = pool.map(_worker, sorted(KS, key=lambda k: -k), chunksize=1)
    else:
        outs = [run_one(k) for k in KS]
    ctx.prec = 512
    tab = {}
    for o in sorted(outs, key=lambda o: o["k"]):
        k = o["k"]
        good = o["rc"] == 0 and not o["fails"] and o.get("n0") is not None
        check(good, f"DB.1-DB.5 (k = {k}, N = {o.get('N')}, rho_in = 3^-{k + 1} * {o.get('rin')}, {o.get('secs', 0):.0f}s): "
                    f"every CU.1-CU.5 check passed on 3^(k+1) L_k" + ("" if good else f"; failures: {o['fails'][:3]}"))
        if not good:
            continue
        S = arb(3) ** (k + 1)
        b = {key: unpack(o[key]) for key in ("rho", "nu", "q", "Gp", "Gm", "Amin", "Bmax", "tau1", "cp", "cm", "osc",
                                             "theta_n0")}
        b["rho"], b["nu"] = b["rho"] / S, b["nu"] / S
        b["cp"], b["cm"], b["osc"] = b["cp"] / S, b["cm"] / S, b["osc"] / S
        b["ev"] = (unpack(o["ev"][0]) / S, unpack(o["ev"][1]) / S, o["ev"][2])
        b["od"] = (unpack(o["od"][0]) / S, unpack(o["od"][1]) / S, o["od"][2])
        b["errs"] = [unpack(e) for e in o["errs"]]
        b["n0"], b["N"] = o["n0"], o["N"]
        tab[k] = b
    RESULT["tab"] = tab
    if sorted(tab) != sorted(KS):
        check(False, "DB: some k not certified; DB.6/DB.7 skipped")
    else:
        rows = []
        for k in sorted(tab):
            b = tab[k]
            rows.append(f"k={k}: rho = {b['rho'].str(13)}, nu_2 = {b['nu'].str(13)}, q = {b['q'].str(12)}, "
                        f"Gamma_+ = {b['Gp'].str(8)}, Gamma_- = {b['Gm'].str(8)}, n_0 = {b['n0']}, "
                        f"c_+ = {b['cp'].str(7)}, c_- = {b['cm'].str(7)}, even [{b['ev'][0].str(5)}, {b['ev'][1].str(5)}] "
                        f"|q|^n, n >= {b['ev'][2]}; odd [{b['od'][0].str(5)}, {b['od'][1].str(5)}] |q|^n, n >= {b['od'][2]}")
        print("   DB summary:\n   " + "\n   ".join(rows), flush=True)
        # consistency with CY.5 / CZ.6
        cons = all(bool(abs(tab[k]["q"] - arb(Q_PRIOR[k])) <= arb(Q_TOL[k])) for k in Q_PRIOR if k in tab)
        check(cons, "DB.3: q^(k) agrees with CY.5/CZ.6's certified values (consistency) at k = "
                    + ", ".join(str(k) for k in Q_PRIOR if k in tab))
        # DB.6: comparison with DA's leading term
        km = arb(KAP_M[0]) + arb(0, KAP_M[1])
        kp = 3 * km
        ratios = {}
        for k in sorted(tab):
            lead_p, lead_m = (2 * k + 3) * km / 12, (2 * k + 3) * kp / 12
            ratios[k] = (tab[k]["Gp"] / lead_p, tab[k]["Gm"] / lead_m)
        within = [k for k in sorted(tab) if all(bool(abs(r - 1) < arb("0.02")) for r in ratios[k])]
        k2 = None
        for k in sorted(tab):
            if all(j in within for j in tab if j >= k):
                k2 = k
                break
        RESULT.update(ratios=ratios, K_2pc=k2)
        check(k2 is not None and k2 <= 28,
              "DB.6: Gamma_+^(k)/((2k+3)/12 kappa_-), Gamma_-^(k)/((2k+3)/12 kappa_+) = "
              + "; ".join(f"k={k}: {ratios[k][0].str(5)}, {ratios[k][1].str(5)}" for k in sorted(tab))
              + f"; both within 2% of 1 for every k >= {k2} in 13..28 (DA.2's bound is void there)")
        # DB.7: monotonicity of Gamma_+- in k (needs k = 13..28 all present)
        if sorted(tab) == list(range(13, 29)):
            gp = {k: (arb(lo), arb(hi)) for k, ((lo, hi), _) in G_PRIOR.items()}
            gm = {k: (arb(lo), arb(hi)) for k, (_, (lo, hi)) in G_PRIOR.items()}
            for k in tab:
                gp[k] = (arb(tab[k]["Gp"].lower()), arb(tab[k]["Gp"].upper()))
                gm[k] = (arb(tab[k]["Gm"].lower()), arb(tab[k]["Gm"].upper()))
            inc = all(bool(gp[k][1] < gp[k + 1][0]) and bool(gm[k][1] < gm[k + 1][0]) for k in range(1, 28))
            check(inc, "DB.7: Gamma_+^(k) and Gamma_-^(k) strictly increase for k = 1..28 (CT.4/CU.4/CV.4 as stated for "
                       "k <= 12): Gamma_+ = " + ", ".join(f"{float(gp[k][0]):.6f}" for k in range(1, 29)))
            RESULT["inc"] = inc
            n0s = {k: tab[k]["n0"] for k in tab}
            RESULT["n0s"] = n0s
            print(f"   Corollary (CT, CU, CV, DB, DA): CT.5's two-sided bound holds for every k >= 1, n >= n_0(k); here "
                  f"n_0(13..28) = {[n0s[k] for k in sorted(n0s)]}", flush=True)
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem DB checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
