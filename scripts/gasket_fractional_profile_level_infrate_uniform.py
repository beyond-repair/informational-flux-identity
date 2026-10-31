#!/usr/bin/env python3
"""Theorem DA: a uniform-in-k rate for inf R_n^(k) -> rho(L_k) (CU/CV/CY/CZ's open item), with the Riesz-projection
constants Gamma_+-^(k) tracked in k: Gamma_+-^(k) = ((2k+3)/12) kappa_-+ + O(delta_k), delta_k/k -> 0 explicit.

Setting exactly as in Theorems CT-CZ (GASKET.md).  L_k f = u_1 f o sigma_1 + u_2 f o sigma_2 on X = [sqrt5/2, 5/2],
h_n = L_k^n 1, R_n^(k) = h_(n+1)/h_n.  CY's frame on H^oo(D), D = D(9/5, 4/5) (X inside D):
    L_k (h_k f) = 3^(-k-1) h_k Lt_k f,   Lt_k = C + E_k,   C f = f o sigma_1,   E_k f = a_k f o sigma_2,   h_k = P0 P1^k.
So with g = 1/h_k (holomorphic on D', zero-free P0, P1; CY.1):   h_n = 3^(-(k+1)n) h_k F_n,  F_n = Lt_k^n g,  and
    R_n = 3^(-k-1) F_(n+1)/F_n   (the factor h_k cancels).
CY.3: for k >= 24, Lt_k has two eigenvalues lambda_0 = 3^(k+1) rho(L_k), lambda_1 = 3^(k+1) nu_2^(k) in |z| >= 9/40, real and
simple, and none else there.  So F_n = lambda_0^n a + lambda_1^n b + e_n with a = Pi_0 g, b = Pi_1 g (Riesz projections) and
e_n = (2 pi i)^(-1) oint_(|z| = 9/40) z^n (z - Lt_k)^(-1) g dz.  B/A in CT.4/CU.4 is b/a (L_k's projections are
h_k Pi_j h_k^(-1)), so Gamma_+ = max_X b/a, Gamma_- = max_X (-b/a).

Left eigenfunctionals (exact algebra; CZ's Schur complement).  W = ker ell, Q f = f - ell(f) kappa, T = (lambda_1 - C - Q E)|_W^(-1);
W0 = ker delta, Q0 f = f - f(s*), T0 = (lambda_0 - C - Q0 E)|_W0^(-1).  Then
    ell~(f) = ell(f) + ell(E T Q f)   satisfies ell~ Lt = lambda_1 ell~, ell~(kappa) = 1,
    delta~(f) = f(s*) + delta(E T0 Q0 f)  satisfies delta~ Lt = lambda_0 delta~, delta~(1) = 1,
right eigenvectors phi_1 = kappa + w, phi_0 = 1 + w0 (CZ), and Pi_1 f = ell~(f) phi_1/ell~(phi_1), Pi_0 f = delta~(f) phi_0/delta~(phi_0).
Exactly g(s*) = 1 and ell(g) = g'(s*) = -(k Phi'(s*) + Psi'(s*)) = -(2k+3)/12 (Phi'(s*) = 1/6, Psi'(s*) = 1/4: geometric sums
along sigma_1'(s*) = -1/3).  Hence, on X,
    b/a = m (kappa + w)/(1 + w0),   m = [ell~(g)/ell~(phi_1)]/[delta~(g)/delta~(phi_0)] = -(2k+3)/12 + small,
and Gamma_+- = ((2k+3)/12) kappa_-+ + O(delta_k), kappa_+ = max_X kappa = kappa(5/2), kappa_- = -min_X kappa = -kappa(sqrt5/2),
and kappa_+ = 3 kappa_- exactly (sigma_1(5/2) = sqrt5/2).

Theorem DA (checked here; exit 1 on any failure).
  DA.1 [proved; arb] kappa' > 0 on X (Taylor arithmetic on a cover of X, Cauchy tail bounds), so kappa_+ = kappa(5/2),
       kappa_- = -kappa(sqrt5/2), enclosed; Phi'(s*) = 1/6, Psi'(s*) = 1/4 (exact geometric sums), checked by Taylor
       arithmetic: (P0 P1^k)'(s*) = (2k+3)/12 for k = 0, 1.
  DA.2 [proved from CY.1-CY.3, CZ.2, DA.1; every integer k >= K_DA = 29] |Gamma_+-^(k) - ((2k+3)/12) kappa_-+| <= delta_k, with
       delta_k explicit per k for K_DA <= k <= K1 = 200 and by monotone majorants for k >= K1; delta_k/(2k+3) -> 0
       (e.g. delta_40 <= 0.226, delta_60 <= 0.00969, delta_100 <= 2.1e-5).  For k = 24..28 the same bounds give
       delta_k >= Gamma_+ (not covered).
  DA.3 [proved; every k >= 29, every n >= n0(k)] CT.5's two-sided bound
           rho(1-q)|q|^n (Gamma_+- - theta_n) <= rho(L_k) - inf_X R_n^(k) <= rho(1-q)|q|^n (Gamma_+- + theta_n)
       (+ for even n, - for odd n), theta_n decreasing in n, theta_(n0(k)) < min Gamma_+- - delta_k, n0(k) <= 16 + k/3.
  DA.4 [corollary] uniform in k: for every k >= 29 and n >= 16 + k/3,
           (rho - inf R_n)/(rho (1-q)|q|^n) = ((2k+3)/12) kappa_-+ (1 + r),  |r| <= rel_k,  rel_k -> 0 (<= 1% for k >= 53);
       hence lim_k lim_(n even) (rho - inf R_n)/(k rho |q|^n) = (2/9) kappa_-, lim_k lim_(n odd) ... = (2/9) kappa_+ = (2/3) kappa_-.
Kept limits: 13 <= k <= 28 not covered (CY.3's localisation starts at 24, DA's bounds at 29; CV covers k <= 12); n < n0(k)
not covered and n0(k) grows linearly in k (sup_D |1/h_k| grows like sup|1/P1|^k = 1.0756^k); delta_k is not the exact
second-order rate; Gamma_+- are located only to +- delta_k (crude for k < 40).
Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import math
import os
import sys
import time
from fractions import Fraction as Fr

from flint import arb, arb_series, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_ratio_rate as cz  # noqa: E402

cy = cz.cy
FAILS: list[str] = []
RESULT: dict = {}
PREC = 160
K0 = 24
K1 = int(os.environ.get("DA_K1", "200"))     # per-k certification for K0..K1, monotone majorants beyond
RIN = Fr(9, 40)                               # CY.3's inner radius
NA, NB = int(os.environ.get("DA_NA", "16")), int(os.environ.get("DA_NB", "3"))    # n0(k) <= NA + k/NB (checked)
NXCOV = 64                                    # cover of X for kappa'
REPORT = (24, 25, 30, 40, 50, 60, 80, 100, 150, 200)

_a = cy._a
upper = cy.upper
lower = cy.lower
up_str = cy.up_str


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


def pm(x) -> arb:
    """the ball [-x, x] for x >= 0 (upper bound of x)."""
    return arb(0, upper(abs(arb(x))))


# ----------------------------------------------------------------------------------------- DA.1 kappa on X
def da1() -> bool:
    ctx.cap = 2
    steps = RESULT["steps"]
    c, RP = _a(cy.C0), _a(cy.RP)
    lo = arb(5).sqrt() / 2
    hi = arb(5) / 2
    ok = bool(lo - (c - _a(cy.R0)) > 0) and bool(c + _a(cy.R0) - hi > 0)       # X inside D
    # kappa' > 0 on X: cover [sqrt5/2, 5/2] by NXCOV balls
    dmin = None
    for i in range(NXCOV):
        a0 = lo + (hi - lo) * i / NXCOV
        a1 = lo + (hi - lo) * (i + 1) / NXCOV
        mid = (a0 + a1) / 2
        rad = upper((a1 - a0) / 2 + abs(a1 - a0) * arb("1e-30"))
        s0 = arb(mid.mid(), upper(rad + mid.rad()))
        d0 = lower(RP - abs(s0 - c))
        ok &= bool(d0 > 0)
        _, _, kp = cz.series_products(arb_series([s0, 1]), d0, arb(1), steps)
        dk = kp[1]
        dmin = dk if dmin is None else (dk if float(lower(dk)) < float(lower(dmin)) else dmin)
        ok &= bool(dk > 0)
    # endpoint values
    vals = {}
    for nm, x in (("lo", lo), ("hi", hi)):
        d0 = lower(RP - abs(x - c))
        _, _, kp = cz.series_products(arb_series([x, 1]), d0, arb(1), steps)
        vals[nm] = kp[0]
    kap_p, kap_m = vals["hi"], -vals["lo"]
    ok &= bool(kap_p > 0) and bool(kap_m > 0) and float(kap_p.rad()) < 1e-25 and float(kap_m.rad()) < 1e-25
    # exact: sigma_1(5/2) = sqrt(5/4) = sqrt5/2 and kappa o sigma_1 = -kappa/3, so kappa_+ = 3 kappa_- (checked to 1e-25)
    ok &= Fr(15, 4) - Fr(5, 2) == Fr(5, 4) and bool(abs(kap_p - 3 * kap_m) < arb("1e-25"))
    # Phi'(s*) = 1/6, Psi'(s*) = 1/4: (P0 P1^k)'(s*) = (2k+3)/12 for k = 0, 1 (and P0(s*) = P1(s*) = 1)
    ss = arb(3) / 2
    d0 = lower(RP - abs(ss - c))
    p0, p1, _ = cz.series_products(arb_series([ss, 1]), d0, arb(1), steps)
    ok &= bool(abs(p0[0] - 1) < arb("1e-25")) and bool(abs(p1[0] - 1) < arb("1e-25"))
    ok &= bool(abs(p0[1] - arb(1) / 4) < arb("1e-25")) and bool(abs(p1[1] - arb(1) / 6) < arb("1e-25"))
    # exact: Phi' = sum_j (2/9)(-1/3)^j... = 1/6, Psi' = (log wt_1)'(3/2) sum (-1/3)^j = (1/3)(3/4) = 1/4
    ok &= Fr(2, 9) * Fr(3, 4) == Fr(1, 6) and (Fr(1, 2) - Fr(2, 3) + Fr(1, 2)) * Fr(3, 4) == Fr(1, 4)
    ctx.cap = 10
    RESULT.update(kap_p=kap_p, kap_m=kap_m, KX=upper(kap_p), dkmin=dmin)
    check(ok, f"DA.1: X inside D; kappa' >= {lower(dmin).str(4, radius=False)} > 0 on X ({NXCOV} Taylor cells), so "
              f"kappa_+ = kappa(5/2) = {kap_p.str(15, radius=False)}, kappa_- = -kappa(sqrt5/2) = {kap_m.str(15, radius=False)} "
              f"(radii < 1e-25), kappa_+ = 3 kappa_- exactly (sigma_1(5/2) = sqrt5/2, kappa o sigma_1 = -kappa/3); P0(s*) = P1(s*) = 1, P0'(s*) = Psi'(s*) = 1/4, P1'(s*) = Phi'(s*) = 1/6, so g = 1/(P0 P1^k) has "
              f"g(s*) = 1 and ell(g) = g'(s*) = -(2k+3)/12 exactly")
    return ok


# ----------------------------------------------------------------------------------------- arcs: sup |g|
def arc_g() -> None:
    """per-arc upper bounds of |1/P0| and |1/P1| on dD (CY's 4096 arcs; max modulus: g = 1/(P0 P1^k) holomorphic on D')."""
    i0, i1 = [], []
    for i in range(cy.NARC):
        z = cy.arc(i, cy.NARC)
        p0, p1, _, _ = cy.orbit_products(z)
        i0.append(upper(abs(1 / p0)))
        i1.append(upper(abs(1 / p1)))
    RESULT.update(i0=i0, i1=i1, M0=max(i0, key=float), M1=max(i1, key=float))


def Mg(k: int) -> arb:
    m = arb(0)
    for x, y in zip(RESULT["i0"], RESULT["i1"]):
        m = max(m, upper(x * y ** k), key=float)
    return m


def S_of(r):
    return cz.S_of(r)


def core(k: int, frozen: dict | None = None) -> dict:
    """All DA quantities at integer k >= K0.  frozen=None: actual per-k bounds (CZ's bounds + per-arc sups).
    frozen=dict: monotone majorants valid for every k' >= frozen['K'] (values of the nonincreasing factors at K,
    explicit k-dependence elsewhere: G* B*^k, Lam_k, M0 M1^k, GM BMa^k)."""
    G, B = RESULT["Gst"], RESULT["Bst"]
    Kk, beta = _a(RESULT["KM"]), RESULT["beta"]
    gb = G * B ** k
    Lam = k * abs(RESULT["b1"]) + abs(RESULT["g1"]) + RESULT["SP"] / (2 * RESULT["t"])
    if frozen is None:
        b = cz.bounds(k)
        assert b is not None
        eta, etaK = cz.eta_ref(k)
        R1, eps, r2, lam0_hi = b["R1"], b["eps"], b["r2"], b["lam0_hi"]
        mg = Mg(k)
    else:
        eta = _a(RESULT["GM"]) * RESULT["aBmax"] ** k                   # >= eta'_k (per-arc max <= GM aBmax^k)
        etaK = eta * RESULT["aKmax"]
        R1, r2, lam0_hi = frozen["R1"], frozen["r2"], frozen["lam0_hi"]  # nonincreasing in k (CZ.4's monotone factors)
        eps = eta + Kk * gb * Lam
        mg = RESULT["M0"] * RESULT["M1"] ** k
    qek = etaK + Kk * gb * Lam * Kk
    ok = bool(R1 * eps < 1)
    T1 = R1 / (1 - R1 * eps)                                 # ||T|| on W
    w1 = T1 * qek                                            # ||w||, phi_1 = kappa + w
    eps0 = eta + gb
    R0 = beta / (arb(5) / 6) + (1 + beta) * RESULT["S_half"]
    ok &= bool(R0 * eps0 < 1)
    T0 = R0 / (1 - R0 * eps0)                                # ||T0|| on W0
    w0 = T0 * eps0                                           # ||w0||, phi_0 = 1 + w0
    m0 = -arb(2 * k + 3) / 12
    Qg = mg + abs(m0) * Kk                                   # ||Q g|| = ||g - ell(g) kappa||
    Q0g = mg + 1                                             # ||Q0 g|| = ||g - 1||
    E1 = gb * Lam * T1 * Qg                                  # |ell~(g) - ell(g)|
    E2 = gb * Lam * T1 * w1                                  # |ell~(phi_1) - 1|
    E3 = gb * T0 * Q0g                                       # |delta~(g) - 1|
    E4 = gb * T0 * w0                                        # |delta~(phi_0) - 1|
    Es = [upper(x) for x in (E1, E2, E3, E4, w0)]
    ok &= all(bool(x < arb(1) / 2) for x in Es[1:])
    Lam1 = (m0 + pm(Es[0])) / (1 + pm(Es[1]))                # ell~(g)/ell~(phi_1)
    Lam0 = (1 + pm(Es[2])) / (1 + pm(Es[3]))                 # delta~(g)/delta~(phi_0)
    m = Lam1 / Lam0
    dm = upper(abs(m - m0))
    mabs = upper(abs(m))
    KX = RESULT["KX"]
    delta = upper(dm * KX + mabs * (w1 + KX * w0) / (1 - w0))
    mina = lower(lower(Lam0) * (1 - w0))
    supb = upper(upper(abs(Lam1)) * (KX + w1))
    maxG = upper(mabs * (KX + w1) / (1 - w0))
    ok &= bool(mina > 0)
    # resolvent on |z| = 9/40: (Res) with |z - 1| >= 31/40, |z + 1/3| >= 13/120; Neumann with ||E_k|| <= eta
    r = _a(RIN)
    Rc = 1 / _a(1 - RIN) + beta / _a(Fr(1, 3) - RIN) + RESULT["Q"] * RESULT["S_rin"]
    ok &= bool(eta * Rc < 1)
    Kin = Rc / (1 - eta * Rc)
    tau1 = upper(r * Kin * mg)                               # ||e_n||/lambda_0^n <= tau1 r^n   (lambda_0 > 1)
    q_hi = upper(arb(1) / 3 + r2)                            # |q| = |lambda_1|/lambda_0 <= |lambda_1|
    l1_lo = lower(arb(1) / 3 - r2)                           # |lambda_1| >= 1/3 - r2
    q_lo = lower(l1_lo / lam0_hi)
    ok &= bool(q_hi < 1) and bool(r < l1_lo) and bool(q_lo > 0)
    gam_p = (2 * k + 3) * RESULT["kap_m"] / 12
    gam_m = (2 * k + 3) * RESULT["kap_p"] / 12
    gmin = lower(min(lower(gam_p - delta), lower(gam_m - delta), key=float))
    return dict(k=k, ok=ok, gb=upper(gb), eta=upper(eta), Lam=Lam, R1=R1, eps=eps, r2=r2, lam0_hi=lam0_hi, mg=mg, T1=upper(T1),
                T0=upper(T0), w1=upper(w1), w0=upper(w0), E=Es, m=m, dm=dm, delta=delta, mina=mina, supb=supb, maxG=maxG,
                Kin=upper(Kin), tau1=tau1, q_hi=q_hi, q_lo=q_lo, l1_lo=l1_lo, gam_p=gam_p, gam_m=gam_m, gmin=gmin)


def theta(b: dict, n: int) -> arb:
    """CT.5's theta_n with DA's constants: |(rho - R_n)/(rho(1-q)|q|^n) - (-1)^n b/a| <= theta_n on X."""
    r = _a(RIN)
    tn = b["tau1"] * r ** n
    tn1 = tn * r
    epsn = (b["q_hi"] ** n * b["supb"] + tn) / b["mina"]
    # Delta_n = (tau_n + tau_(n+1))/((1 - q)|q|^n min a), 1 - q >= 1 + q_lo; tau_n/|q|^n <= tau1 (r lambda_0/|lambda_1|)^n...
    # use |q|^n lambda_0^n = |lambda_1|^n and tau_n = ||e_n||/lambda_0^n: tau_n/|q|^n <= r Kin mg (r/|lambda_1|)^n
    ratio = b["tau1"] * (r / b["l1_lo"]) ** n
    Dn = ratio * (1 + r) / ((1 + b["q_lo"]) * b["mina"])
    if not bool(epsn < 1):
        return arb("inf")
    return upper((b["maxG"] * epsn + Dn) / (1 - epsn))


def n0_of(b: dict, nmax: int = 2000) -> int | None:
    for n in range(0, nmax):
        th = theta(b, n)
        if th.is_finite() and bool(th < b["gmin"]):
            return n
    return None


def da2_da3() -> bool:
    """Per-k certification for k = K0..K1.  K_DA = the first k with min Gamma_+- - delta_k > theta_n for some n (below it
    delta_k >= Gamma: kept limit); from K_DA on every k must pass."""
    ok = True
    tab = {}
    kda = None
    fails_below = []
    for k in range(K0, K1 + 1):
        b = core(k)
        n0 = n0_of(b) if b["ok"] else None
        b["n0"] = n0
        good = b["ok"] and n0 is not None and n0 <= NA + k / NB
        if good:
            b["theta_n0"] = theta(b, n0)      # theta_n decreasing in n: bases q_hi, r, r/|lambda_1| < 1 (checked in core)
            if kda is None:
                kda = k
            tab[k] = b
        elif kda is None:
            fails_below.append((k, b))
        else:
            ok = False
            print(f"  failure at k = {k}: ok={b['ok']}, n0={n0}", flush=True)
            break
    ok &= kda is not None
    RESULT.update(tab=tab, KDA=kda, below=fails_below)
    rows = []
    for k, b in fails_below:
        rows.append(f"k={k} (not covered): delta_k <= {up_str(b['delta'], 4)} vs Gamma_+ ~ {b['gam_p'].str(5, radius=False)}")
    for k in sorted(set(REPORT) | ({kda} if kda else set())):
        if k in tab:
            b = tab[k]
            rows.append(f"k={k}: Gamma_+ in {b['gam_p'].str(10, radius=False)} +- {up_str(b['delta'], 3)}, Gamma_- in "
                        f"{b['gam_m'].str(10, radius=False)} +- {up_str(b['delta'], 3)} (delta_k/(2k+3) <= "
                        f"{up_str(b['delta'] / (2 * k + 3), 3)}); min a >= {lower(b['mina']).str(6, radius=False)}, "
                        f"sup|b| <= {up_str(b['supb'], 5)}, sup_D|g| <= {up_str(b['mg'], 4)}, tau_1 <= {up_str(b['tau1'], 4)}; "
                        f"n0 = {b['n0']}, theta_n0 <= {up_str(b['theta_n0'], 4)}")
    check(ok, f"DA.2/DA.3 (per k, k = {kda}..{K1}): all Neumann conditions hold; |Gamma_+- - ((2k+3)/12) kappa_-+| <= delta_k; "
              f"theta_n < min Gamma_+- - delta_k from n0(k) <= {NA} + k/{NB}:\n      " + "\n      ".join(rows))
    return ok


def monotone_ok(x: arb, alpha: arb, beta: arb, K: int, p: int = 1) -> bool:
    """(alpha k + beta)^p x^k nonincreasing for k >= K (alpha, beta >= 0): x (1 + alpha/(alpha K + beta))^p <= 1."""
    return bool(x * (1 + alpha / (alpha * K + beta)) ** p <= 1)


def da_tail() -> bool:
    """k >= K1: monotone majorants.  Every factor frozen at K1 is nonincreasing in k (CZ.4's monotone factors: R1, r2,
    lam0_hi).  The k-dependent pieces are sums of terms c k^p x^k with x in {B*, aBmax, B* M1, B* aBmax/..}; each is
    nonincreasing for k >= K1 by monotone_ok.  Then E2, E3, E4, E1/|m0|, w1, w0, delta_k/(2k+3), tau1 r^(n0(k)) and
    q_hi^(n0(k)) sup|b| (with n0(k) = NA + ceil(k/NB) >= NA + k/NB) are nonincreasing, so the K1 values bound all k >= K1."""
    ok = cz.monotone_factors(K1)
    B, M1, aB = RESULT["Bst"], RESULT["M1"], RESULT["aBmax"]
    b1 = abs(RESULT["b1"])
    lam0 = abs(RESULT["g1"]) + RESULT["SP"] / (2 * RESULT["t"])
    r = _a(RIN)
    fr = dict(K=K1, R1=RESULT["tab"][K1]["R1"], r2=RESULT["tab"][K1]["r2"], lam0_hi=RESULT["tab"][K1]["lam0_hi"])
    bt = core(K1, frozen=fr)
    ok &= bt["ok"]
    # terms and their bases (p = power of a linear factor in k): gb Lam (x = B*, p=1); gb Lam^2 (p=2); gb Lam mg (B* M1, p=1);
    # gb Lam (2k+3) (p=2); eta Lam (aBmax, p=1); gb mg (B* M1, p=0); E1/|m0|, w1 |m0| (p=3 at most)
    for x in (B, B * M1, aB):
        ok &= monotone_ok(x, arb(1), arb(0), K1, p=3)
    # n-dependent pieces at n = n0(k) = NA + ceil(k/NB): (r^(1/NB) M1)^k, (q_hi^(1/NB))^k (2k+3), ((r/l1_lo)^(1/NB) M1)^k
    xr = r ** (arb(1) / NB) * M1
    xq = bt["q_hi"] ** (arb(1) / NB)
    xd = (r / bt["l1_lo"]) ** (arb(1) / NB) * M1
    ok &= bool(xr < 1) and bool(xd < 1) and monotone_ok(xq, arb(1), arb(0), K1, p=2)
    # majorant values at K1 (n = NA + K1/NB, real exponent: the bound uses n >= NA + k/NB)
    nK = NA + arb(K1) / NB
    epsn = (bt["q_hi"] ** nK * bt["supb"] + bt["tau1"] * r ** nK) / bt["mina"]
    Dn = bt["tau1"] * (r / bt["l1_lo"]) ** nK * (1 + r) / ((1 + bt["q_lo"]) * bt["mina"])
    th_rel = upper((bt["maxG"] / (2 * K1 + 3) * epsn + Dn / (2 * K1 + 3)) / (1 - epsn))
    del_rel = upper(bt["delta"] / (2 * K1 + 3))
    kmin = min(RESULT["kap_m"], RESULT["kap_p"], key=float)
    ok &= bool(epsn < arb(1) / 2) and bool(12 * (th_rel + del_rel) < lower(kmin))
    RESULT.update(tail=bt, th_rel=th_rel, del_rel=del_rel)
    check(ok, f"DA tail (every k >= {K1}): monotone majorants (bases B* = {B.str(6, radius=False)}, B* sup|1/P1| = "
              f"{(B * M1).str(6, radius=False)}, sup|B| = {aB.str(6, radius=False)}; at n = {NA} + k/{NB}: "
              f"(9/40)^(1/{NB}) sup|1/P1| = {xr.str(5, radius=False)}, ((9/40)/|lambda_1|)^(1/{NB}) sup|1/P1| = "
              f"{xd.str(5, radius=False)} < 1); sup_dD|1/P0| <= {up_str(RESULT['M0'], 6)}, sup_dD|1/P1| <= {up_str(M1, 6)}; "
              f"for every k >= {K1} and n >= {NA} + k/{NB}: delta_k/(2k+3) <= {up_str(del_rel, 3)}, "
              f"theta_n/(2k+3) <= {up_str(th_rel, 3)}, both nonincreasing in k, 12 (sum) < min kappa_+-")
    return ok


def da4() -> bool:
    """Uniform statement: for every k >= K_DA and n >= NA + k/NB, the normalised gap lies within rel_k of
    ((2k+3)/12) kappa_-+, rel_k = 12 (delta_k + theta_(n_u(k)))/((2k+3) kappa_-) (kappa_- <= kappa_+), n_u = ceil(NA + k/NB);
    rel_k is certified per k up to K1 and bounded by the tail majorant beyond; rel_k -> 0."""
    tab = RESULT["tab"]
    ok = True
    kmin = min(RESULT["kap_m"], RESULT["kap_p"], key=float)
    rel = {}
    for k, b in tab.items():
        nu = math.ceil(NA + k / NB)
        ok &= nu >= b["n0"]
        rel[k] = upper(12 * (b["delta"] + theta(b, nu)) / ((2 * k + 3) * kmin))
    rel_tail = upper(12 * (RESULT["th_rel"] + RESULT["del_rel"]) / kmin)
    ks = sorted(rel)
    # first K with rel_k <= 1/100 for every k >= K (tail included)
    k1pc = None
    if bool(rel_tail <= arb("0.01")):
        for k in reversed(ks):
            if bool(rel[k] <= arb("0.01")):
                k1pc = k
            else:
                break
    ok &= bool(rel[ks[0]] < 1) and k1pc is not None and all(bool(rel[k] < 1) for k in ks)
    lim_p = 2 * RESULT["kap_m"] / 9
    lim_m = 2 * RESULT["kap_p"] / 9
    RESULT.update(rel=rel, rel_tail=rel_tail, K1pc=k1pc, lim_p=lim_p, lim_m=lim_m)
    shown = [k for k in (ks[0], 30, 35, 40, 50, 60, 80, 100, K1) if k in rel]
    check(ok, f"DA.4: for every k >= {RESULT['KDA']} and every n >= {NA} + k/{NB}, (rho - inf R_n^(k))/(rho (1-q)|q|^n) = "
              f"((2k+3)/12) kappa_-+ (1 + r), |r| <= rel_k: " + ", ".join(f"k={k}: {up_str(rel[k], 3)}" for k in shown)
              + f", k >= {K1}: {up_str(rel_tail, 3)}; rel_k <= 1% for every k >= {k1pc}; so, uniformly, lim_k lim_(n even) "
              f"(rho - inf R_n)/(k rho |q|^n) = (2/9) kappa_- = {lim_p.str(12, radius=False)} and lim_k lim_(n odd) = (2/9) kappa_+ = "
              f"{lim_m.str(12, radius=False)}")
    return ok


def observation() -> None:
    """Floats only (no claim): CV.4's Gamma_+- for k = 5, 8, 12 against ((2k+3)/12) kappa_-+."""
    cv = {5: (0.416558, 1.129942), 8: (0.583993, 1.607694), 12: (0.824925, 2.33819)}
    kp, km = float(RESULT["kap_p"].mid()), float(RESULT["kap_m"].mid())
    rows = [f"k={k}: CV.4 ({gp}, {gm}) vs ({(2*k+3)/12*km:.6f}, {(2*k+3)/12*kp:.6f})" for k, (gp, gm) in cv.items()]
    print("obs  (floats, no claim, outside DA's range) " + "; ".join(rows))


def main() -> int:
    print("Theorem DA: uniform-in-k rate of inf R_n^(k) -> rho(L_k); Gamma_+-^(k) = ((2k+3)/12) kappa_-+ + o(k). Not thrust.",
          flush=True)
    t0 = time.time()
    FAILS.clear()
    RESULT.clear()
    old = ctx.prec
    ctx.prec = PREC
    try:
        cy.FAILS.clear()
        cy.RESULT.clear()
        cz.FAILS.clear()
        cz.RESULT.clear()
        ok = cy.cy1_disk() and cy.cy1_sups() and cy.cy2() and cy.cy3() and cy.cy4()
        check(ok and cy.RESULT["K0"] == K0, f"CY.1-CY.4 re-run (imported): K0 = {cy.RESULT.get('K0')}")
        cz.RESULT.update({k: cy.RESULT[k] for k in ("cn", "dNC", "beta", "KM", "BM", "GM")})
        cz.RESULT["steps"] = cy.RESULT["entry_steps"]
        ok = cz.cz1() and ok
        old_scan = cz.K_SCAN
        cz.K_SCAN = 60
        try:
            ok = cz.cz2() and ok
        finally:
            cz.K_SCAN = old_scan
        ok = cz.monotone_factors(K0) and ok
        RESULT.update({k: cz.RESULT[k] for k in ("Gst", "Bst", "b1", "g1", "SP", "t", "KM", "GM", "beta", "S_half",
                                                  "aBmax", "aKmax", "steps")})
        RESULT["Q"] = cy.RESULT["Q"]
        RESULT["S_rin"] = cy.RESULT["S"]
        ok = da1() and ok
        arc_g()
        if os.environ.get("DA_SETUP_ONLY") == "1":
            return 0
        ok = da2_da3() and ok
        if ok:
            ok = da_tail() and ok
            ok = da4() and ok
            observation()
    finally:
        ctx.prec = old
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem DA checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
