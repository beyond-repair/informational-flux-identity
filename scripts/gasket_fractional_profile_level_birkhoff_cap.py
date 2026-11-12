#!/usr/bin/env python3
r"""Theorem DP: a certified cap for the whole cone-diameter (Birkhoff) route.  The C^1 eigenvalues of CI/CM give every
Birkhoff contraction ratio a floor: on any cone C_a (a > b*), any block Lt^n and any cone whose order-unit space at h
contains C^1(X), the per-level rate tanh(diam/4)^(1/n) is >= rho(L_1)/r > 0.3170757 (DM.5's 0.363655 is within 1.147 of
it; CN's exact level-sum rate 0.17255 is out of reach).  (DO's open item "a certified floor for the cone-diameter route";
DM.5's floors concern one pair (C_2, Lt^n) at a time.)

Setting (GASKET.md, Theorems CE, CF, CI, CM, DL-DO): s in X = [sqrt5/2, 5/2], Lt f = wt1 f o sigma_1 + wt2 f o sigma_2,
C_a = {f > 0 : log f a-Lipschitz in s}, Theta_K the Hilbert metric of a cone K (CF's definition: Theta = log(beta/alpha),
beta = min{l : l u - v in Kbar}, alpha = max{l : v - l u in Kbar}, Kbar = K u {0}); h = CF.3's positive eigenfunction
(Lt h = r h, h in Kbar_(b*)).  For a positive linear T with T K in K, tau(T; K) = sup Theta(Tf, Tg)/Theta(f, g).

Theorem DP (checked here; exit 1 on any failure):
  DP.1 [proved; identities in sympy] Let K be a convex cone of positive functions with Kbar closed along lines, h in K,
       phi real with c h +- phi in Kbar for some c (phi in the order-unit space V_h), phi not a multiple of h.  Put
       s0 = min{s : s h + phi in Kbar}, s1 = min{s : s h - phi in Kbar}.  Then s0 + s1 >= osc(phi/h) > 0 and, exactly,
           Theta_K(h + t phi, h) = log((1 + t s1)/(1 - t s0))   for 0 < t < 1/max(s0, 0),
       (and with phi -> -phi for t < 0), so Theta_K(h + t phi, h) = |t| (s0 + s1) + O(t^2).
  DP.2 [proved from DP.1 and CF.1] If moreover T h = rho h and T phi = mu phi (mu real), then tau(T; K) >= |mu|/rho, so
       tanh(diam_Theta T K / 4) >= |mu|/rho (and diam >= 4 artanh(|mu|/rho)).
  DP.3 [kept failure, certified witness; DP.2 with CI.4/CM.3's real C^1 eigenfunctions] For every a > b* = 1.68842, every
       n >= 1 and every product of mixed blocks of total length n (DO.3's G_a), and for every Lt^n-invariant cone K with
       C^1(X; R) in V_h (e.g. h-weighted cones {log(f/h) a-Lipschitz}), tau(Lt^n; K) >= (rho(L_1)/r)^n:
           every CF.1-form per-level rate is >= rho(L_1)/r >= 0.1343549/r > 0.3170757
       (CI.4's explicit eigenfunction at mu = 0.1343547 alone gives > 0.3170752).  Diameter floors 4 artanh(q^n) for
       n = 5, 7, 9 lie below DL.3/DM.3/DJ.2's floors (consistent).  DM.5's 0.363655 is < 1.1470 times the cap, DK.4's
       0.53349 < 1.6826 times; CN's exact rate |l_2|/r = 0.17255 is below the cap by a factor > 1.8375, so no level-sum bound
       of CF.3/DN.3/DO.3 form on such cones reaches the spectral rate.
  DP.4 [proved from DP.2 and CM.4] For any invariant cone with the (real-analytic, real) l_2-eigenfunction in V_h the floor is
       |l_2|/r > 0.172547 (universal); if V_h lies in C^2(X), DP.2's mechanism gives nothing above |l_2|/r, since the C^2
       spectrum is {r, l_2} u {|z| <= rho(L_2)} and rho(L_2)/r < 0.101646.  So the cap of DP.3 is a C^1 phenomenon.
Kept limits: DP.3 bounds tanh(diam/4), i.e. every argument that goes through CF.1 or a contraction ratio on such a cone;
bounds not of contraction form (CN's spectral expansion) are not covered.  No Psi_info, Stage 2 still blocked.
Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time

import sympy as sp
from flint import arb, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_essential as ci          # noqa: E402  (CI.1 re-run, CI.4's mu, r)
import gasket_fractional_profile_level_calkin as cl             # noqa: E402  (CL.3 re-run for rho(L_2), l_2)

FAILS: list[str] = []
RESULT: dict = {}
BSTAR = "1.688416"           # b* = 7(5+sqrt5)/30 = 1.688415... (CF.2); checked below
# certified floors of single pairs (as stated): DL.3, DM.3 (n = 7, 9) on C_2; DJ.2 on C_100 (n = 9)
PAIR_FLOORS = {("C_2", 5): "0.0435818", ("C_2", 7): "0.0044106", ("C_2", 9): "0.00044485", ("C_100", 9): "0.013997"}
# DM.4's diameter enclosures on C_2 (as stated), for the observation
DIAM_C2 = {5: ("0.0435818", "0.043583"), 7: ("0.0044106", "0.004411"), 9: ("0.00044485", "0.0004449")}
# per-level rates as stated: DL.5, DM.5 (n = 7, 9), DK.4, DH.4, DG
RATES = {"DL.5 (C_2, n=5)": "0.404994", "DM.5 (C_2, n=7)": "0.378005", "DM.5 (C_2, n=9)": "0.363655",
         "DK.4 (C_100, n=9)": "0.53349", "DH.4 (C_100, n=9)": "0.57117", "DG (C_100, n=8)": "0.58745"}
MU_DM5 = "0.15410"           # DM.5's |mu| bound for other C^1 eigenvalues
CN_C, CS_LIM = "14.60", "9.6274"   # CN.4's constant; CS's limit B(q2 - 1) (comparison only)
B_DO = {13: "1.8123e-6", 22: "2.0158e-10", 49: "2.7740e-22"}   # DO.3's B(k) (comparison only)
TARGET = dict(q_w="0.3170752", q1="0.3170757", q1_up="0.3170762", ratio_dm="1.1470", ratio_dk="1.6826",
              spec="0.172547", gap="1.8375", l2r="0.101646")


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- DP.1/DP.2 exact algebra
def symbolic_checks():
    lam, t, s0, s1, q = sp.symbols("lambda t s0 s1 q", real=True)
    h, phi = sp.symbols("h phi")
    # beta: lam (h + t phi) - h = lam t ( ((lam-1)/(lam t)) h + phi )
    e1 = sp.expand(lam * (h + t * phi) - h - lam * t * (((lam - 1) / (lam * t)) * h + phi))
    # alpha: h - lam (h + t phi) = lam t ( ((1/lam - 1)/t) h - phi )
    e2 = sp.expand(h - lam * (h + t * phi) - lam * t * (((1 / lam - 1) / t) * h - phi))
    sol_b = sp.solve(sp.Eq((1 - 1 / lam) / t, s0), lam)
    sol_a = sp.solve(sp.Eq((1 / lam - 1) / t, s1), lam)
    ok_sol = (len(sol_b) == 1 and sp.simplify(sol_b[0] - 1 / (1 - t * s0)) == 0
              and len(sol_a) == 1 and sp.simplify(sol_a[0] - 1 / (1 + t * s1)) == 0)
    # monotonicity: (1 - 1/lam)/t increasing in lam (t > 0), (1/lam - 1)/t decreasing
    ok_mono = sp.simplify(sp.diff((1 - 1 / lam) / t, lam) - 1 / (lam ** 2 * t)) == 0
    theta = sp.log((1 + t * s1) / (1 - t * s0))
    ser = sp.series(theta, t, 0, 3).removeO()
    ok_ser = sp.simplify(ser - (t * (s0 + s1) + t ** 2 * (s0 ** 2 - s1 ** 2) / 2)) == 0
    # (s0 h + phi) + (s1 h - phi) = (s0 + s1) h
    ok_sum = sp.expand((s0 * h + phi) + (s1 * h - phi) - (s0 + s1) * h) == 0
    # DP.2: ratio of Theta(h + t q phi, h) to Theta(h + t phi, h) -> |q| (q > 0 directly; q < 0 via phi -> -phi)
    qp = sp.symbols("qp", positive=True)
    a0, a1 = sp.symbols("a0 a1", positive=True)
    th = lambda tt, x0, x1: sp.log((1 + tt * x1) / (1 - tt * x0))
    lim_pos = sp.limit(th(t * qp, a0, a1) / th(t, a0, a1), t, 0, "+")
    lim_neg = sp.limit(th(t * qp, a1, a0) / th(t, a0, a1), t, 0, "+")     # mu < 0: direction -phi, |q| = qp
    # T(h + t phi) = rho h + t mu phi = rho (h + t (mu/rho) phi): projective invariance
    rho, mu = sp.symbols("rho mu", positive=True)
    ok_proj = sp.expand(rho * h + t * mu * phi - rho * (h + t * (mu / rho) * phi)) == 0
    bstar = 7 * (5 + sp.sqrt(5)) / 30
    ok_bstar = sp.Rational(1688415, 10 ** 6) < bstar < sp.Rational(1688416, 10 ** 6)
    check(e1 == 0 and e2 == 0 and ok_sol and ok_mono,
          "DP.1 (sympy): l u - v = l t(((l-1)/(l t)) h + phi), v - l u = l t(((1/l-1)/t) h - phi) for u = h + t phi, v = h; "
          "so beta = 1/(1 - t s0), alpha = 1/(1 + t s1) exactly (monotone in l)")
    check(ok_ser and ok_sum,
          "DP.1 (sympy): Theta_K(h + t phi, h) = log((1+t s1)/(1-t s0)) = t (s0+s1) + t^2 (s0^2-s1^2)/2 + O(t^3); "
          "(s0 h + phi) + (s1 h - phi) = (s0+s1) h, so s0 + s1 >= osc(phi/h) > 0 for phi not in R h")
    check(sp.simplify(lim_pos - qp) == 0 and sp.simplify(lim_neg - qp) == 0 and ok_proj,
          "DP.2 (sympy): Theta(h + t q phi, h)/Theta(h + t phi, h) -> |q| as t -> 0+ for q > 0 and q < 0; "
          "T(h + t phi) = rho (h + t (mu/rho) phi), so tau(T; K) >= |mu|/rho")
    check(ok_bstar, "CF.2: b* = 7(5+sqrt5)/30 in (1.688415, 1.688416); every a > b* has Lt C_a in C_a and h in C_(b*)")


# ----------------------------------------------------------------- DP.3/DP.4 certified numbers
def consequences():
    ctx.prec = 128
    r_lo, r_hi = arb(ci.R_LO), arb(ci.R_HI)
    lam_lo, lam_up = arb(ci.LAM_LO), arb(ci.LAM_UP)
    mu_w = arb(ci.MU_MAX)
    q_w = mu_w / r_hi
    q1 = lam_lo / r_hi
    q1_up = lam_up / r_lo
    check(bool(q_w > arb(TARGET["q_w"])),
          f"DP.3 (explicit witness): CI.4's real C^1 eigenfunction v = sum mu^n K^n g at mu = {ci.MU_MAX} gives "
          f"tau(Lt^n; K)^(1/n) >= mu/r >= {q_w.lower().str(9, radius=False)} > {TARGET['q_w']}")
    check(bool(q1 > arb(TARGET["q1"])) and bool(q1_up < arb(TARGET["q1_up"])),
          f"DP.3 (CM.3: every real |mu| < rho(L_1) is a C^1 eigenvalue): cap rho(L_1)/r in [{q1.lower().str(9, radius=False)}, "
          f"{q1_up.upper().str(9, radius=False)}], so every CF.1-form per-level rate on a C^1-open cone is > {TARGET['q1']}")
    # diameter floors vs single-pair floors
    floors = {}
    ok = True
    parts = []
    for (cone, n), fl in PAIR_FLOORS.items():
        x = q1 ** n
        dfl = 4 * x.atanh()
        floors[n] = dfl
        ok = ok and bool(dfl < arb(fl))
        parts.append(f"n={n}: 4 artanh(q^{n}) < {dfl.upper().str(6, radius=False)} < {fl} ({cone})")
    check(ok, "DP.3 consistency: every diam_Theta Lt^n K >= 4 artanh(q^n); " + "; ".join(parts) +
              " (the single-pair floors of DL.3/DM.3/DJ.2 lie above, as they must)")
    rat = {}
    for name, rv in RATES.items():
        rat[name] = arb(rv) / q1
    ok = bool(rat["DM.5 (C_2, n=9)"] < arb(TARGET["ratio_dm"])) and bool(rat["DK.4 (C_100, n=9)"] < arb(TARGET["ratio_dk"]))
    ok = ok and all(bool(v > 1) for v in rat.values())
    check(ok, "DP.3: stated rates / cap: " + ", ".join(f"{k} {RATES[k]} -> < {v.upper().str(5, radius=False)}"
                                                        for k, v in rat.items()) +
              f"; so the C_2 route (DM.5) is within {TARGET['ratio_dm']} of what any C^1-open cone can give")
    mu_ratio = arb(MU_DM5) / lam_lo
    check(bool(mu_ratio > 1) and bool(mu_ratio < arb("1.1470")),
          f"DP.3: Birkhoff's eigenvalue bound |mu| <= tau^(1/n) r can at best reach rho(L_1) >= {ci.LAM_LO} (CJ: exactly the "
          f"C^1 radius); DM.5's {MU_DM5} is {mu_ratio.upper().str(5, radius=False)} times it")
    # DP.4: spectral rate, C^2 picture
    l2 = arb(cl.L2_ABS, float(cl.L2_ERR))
    spec_lo, spec_up = l2 / r_hi, l2 / r_lo
    gap = q1 / spec_up
    l2lo, l2up = arb(cl.CW[2][0]), arb(cl.CW[2][1])
    l2r = l2up / r_lo
    check(bool(spec_lo > arb(TARGET["spec"])) and bool(gap > arb(TARGET["gap"])),
          f"DP.3/DP.4: CN's exact level-sum rate |l_2|/r in [{spec_lo.lower().str(8, radius=False)}, "
          f"{spec_up.upper().str(8, radius=False)}]; cap / spectral > {gap.lower().str(6, radius=False)}: "
          f"no CF.1-form bound on a C^1-open cone reaches it")
    check(bool(l2r < arb(TARGET["l2r"])) and bool(l2up < l2 - arb("1e-6")),
          f"DP.4: rho(L_2) <= {cl.CW[2][1]} < |l_2|, rho(L_2)/r < {TARGET['l2r']}: on cones with V_h in C^2 the eigenfunction "
          f"floor is exactly |l_2|/r (CM.4: spec(C^2) = {{r, l_2}} u {{|z| <= rho(L_2)}})")
    # observation (certified arithmetic on DM.4's enclosures as stated)
    d7lo, d7hi = map(arb, DIAM_C2[7])
    d9lo, d9hi = map(arb, DIAM_C2[9])
    d5lo, d5hi = map(arb, DIAM_C2[5])
    rr97 = ((d9lo / d7hi).sqrt(), (d9hi / d7lo).sqrt())
    rr75 = ((d7lo / d5hi).sqrt(), (d7hi / d5lo).sqrt())
    c9 = ((d9lo / 4).tanh() / q1_up ** 9, (d9hi / 4).tanh() / q1 ** 9)
    print(f"     observation: two-level diameter ratios on C_2 (DM.4's enclosures): (D7/D5)^(1/2) in "
          f"[{rr75[0].lower().str(6, radius=False)}, {rr75[1].upper().str(6, radius=False)}], (D9/D7)^(1/2) in "
          f"[{rr97[0].lower().str(6, radius=False)}, {rr97[1].upper().str(6, radius=False)}] (cap {q1.lower().str(7, radius=False)}); "
          f"tanh(D9/4)/q^9 in [{c9[0].lower().str(6, radius=False)}, {c9[1].upper().str(6, radius=False)}]")
    cf = float(c9[1].upper())
    qf = float(q1.upper())
    n33 = math.ceil(math.log(cf) / math.log(0.33 / qf))
    n35 = math.ceil(math.log(cf) / math.log(0.35 / qf))
    print(f"     observation (floats, extrapolation, no claim): if tanh(D_n/4) ~ {cf:.3f} q^n persisted, rate 0.35 needs "
          f"n >= {n35} and rate 0.33 needs n >= {n33} (2^n words; DM's n = 9 is 512 words, ~33 min)")
    q2 = float(spec_up.upper())
    for k, bv in B_DO.items():
        truth = float(CS_LIM) * q2 ** k
        print(f"     observation: DO.3's B({k}) < {bv} vs CN-scale deviation ~{truth:.3e} (ratio ~{float(bv) / truth:.2e})")
    RESULT.update(q_w=q_w, q1=q1, q1_up=q1_up, floors=floors, rat=rat, mu_ratio=mu_ratio, spec=(spec_lo, spec_up),
                  gap=gap, l2r=l2r, rr97=rr97, c9=c9)


def main() -> int:
    t0 = time.time()
    print("Theorem DP: a cap for the cone-diameter route: every CF.1-form per-level rate on a C^1-open cone is "
          ">= rho(L_1)/r > 0.3170757. Not thrust.", flush=True)
    symbolic_checks()
    ci.collatz_two_sided()          # CI.1 re-run: rho(L_1) in [0.1343549, 0.1343551]
    cl.collatz_two_sided(2)         # CL.3 re-run: rho(L_2) in [0.0430702, 0.0430704]
    consequences()
    print(f"elapsed {time.time() - t0:.1f}s")
    nf = len(FAILS) + len(ci.FAILS) + len(cl.FAILS)
    if nf:
        print(f"{nf} FAILED")
        return 1
    print("all Theorem DP checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
