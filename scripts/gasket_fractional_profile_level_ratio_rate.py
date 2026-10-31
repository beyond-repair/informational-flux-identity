#!/usr/bin/env python3
"""Theorem CZ: the exact rate and the sign of q^(k) + 1/3 (CY's kept limit), and the exact first correction to
3^(k+1) rho(L_k) - 1 (CX's c_k -> 1 with an explicit error), for every integer k >= K0 = 24.

Setting exactly as in Theorem CY (GASKET.md; its CY.1-CY.4 are imported and re-run): on H^oo(D), D = D(9/5, 4/5),
    Lt_k = 3^(k+1) h_k^(-1) L_k h_k = C + E_k,   C f = f o sigma_1,   E_k f = a_k (f o sigma_2),   a_k = G B^k,
lambda_0 = 3^(k+1) rho(L_k), lambda_1 = 3^(k+1) nu_2^(k) (CY.3: real, simple, |lambda_0 - 1| < 1/2,
|lambda_1 + 1/3| <= r_2(k)), q^(k) = lambda_1/lambda_0.  s* = 3/2, t = sigma_2(s*) = sqrt(21)/2.
Left eigenfunctionals of C: delta(f) = f(s*) (eigenvalue 1, right eigenvector 1) and ell(f) = f'(s*) (eigenvalue -1/3,
right eigenvector kappa, kappa'(s*) = 1, kappa(s*) = 0).

Schur complement (exact algebra).  With W = ker ell, Q f = f - ell(f) kappa: if (lambda_1 - C - Q E_k) is invertible on
W, the eigenvector is kappa + w, w = (lambda_1 - C - Q E_k)|_W^(-1) Q E_k kappa, and
    lambda_1 = -1/3 + ell(E_k kappa) + ell(E_k w).
Likewise with W0 = ker delta, Q0 f = f - f(s*):  lambda_0 = 1 + delta(E_k 1) + delta(E_k w0) = 1 + a_k(s*) + delta(E_k w0).
First-order terms (exact): ell(E_k f) = a_k'(s*) f(t) + a_k(s*) f'(t)/(2t), a_k(s*) = G* B*^k, a_k'(s*) = G* B*^k (k b1 + g1),
b1 = B'(s*)/B*, g1 = G'(s*)/G*.  So
    lambda_1 + 1/3 = G* B*^k (k alpha' + gamma' - 1/3) + rho_1,   alpha' = b1 kappa(t),  gamma' = g1 kappa(t) + kappa'(t)/(2t) + 1/3,
    lambda_0 - 1   = G* B*^k (1 + rho_0'),
and q + 1/3 = N_k/lambda_0 with N_k = lambda_1 + 1/3 + (lambda_0 - 1)/3 = G* B*^k (k alpha' + gamma') + rho_1 + G* B*^k rho_0'/3.

Bounds (norms on H^oo(D)): ||E_k|| <= eta_k (CY), |f(t)| <= ||f||, |f'(t)| <= SP ||f||, SP = R/(R^2 - |t - c|^2) (Schwarz-Pick),
so |ell(E_k f)| <= G* B*^k Lam_k ||f||, Lam_k = k|b1| + |g1| + SP/(2t);  ||Q E_k||_W, ||Q E_k kappa||/||kappa|| <= eps_k = eta_k + K G* B*^k Lam_k;
(z - C)^(-1) on W: w = w(s*) 1 + g, g in V_2, ||g|| <= 2||w||: <= R_1(z) = 1/|z - 1| + 2 S(|z|);
(z - C)^(-1) on W0: f = f'(s*) kappa + g, ||g|| <= (1 + beta)||f||: <= R_0(z) = beta/|z + 1/3| + (1 + beta) S(|z|);
||Q0 E_k||_W0, ||Q0 E_k 1|| <= eps0_k = eta_k + G* B*^k.  Hence (Neumann series)
    |rho_1| <= G* B*^k Lam_k K R_1 eps_k/(1 - R_1 eps_k),     |rho_0'| <= w0_k = R_0 eps0_k/(1 - R_0 eps0_k),
with R_1 evaluated at |z| >= 1/3 - r_2(k), |z - 1| >= 4/3 - r_2(k), and R_0 at |z + 1/3| >= 5/6, |z| >= 1/2.

Theorem CZ (checked here; exit 1 on any failure).
  CZ.1 [proved; arb, 160 bits, Taylor arithmetic with explicit product tails and Cauchy bounds for their derivatives]
       B*, G*, b1, g1, kappa(t), kappa'(t), SP; alpha* = G* alpha' and gamma* = G* gamma' enclosed.
  CZ.2 [proved from CY.1-CY.3 and CZ.1, every integer k >= 24] the two expansions above with the explicit bounds:
       c_k = (3^(k+1) rho(L_k) - 1)/(G* B*^k) satisfies |c_k - 1| <= w0_k; |3^(k+1) nu_2 + 1/3 - G* B*^k (k alpha' + gamma' - 1/3)|
       <= G* B*^k Rem1_k.
  CZ.3 [proved from CZ.2] |q^(k) + 1/3 - (alpha* k + gamma*) B*^k| <= Err_k B*^k for every k >= 24, Err_k/k -> 0, so
       (q^(k) + 1/3)/(k B*^k) -> alpha* < 0: exact rate k B*^k.
  CZ.4 [proved from CZ.2 for k >= K1 (bound monotone in k), CU's routine for 24 <= k < K1] q^(k) < -1/3 for every k >= 24;
       with CY.6 (k = 8..23), CV.5 (k = 4..12, as stated): q^(k) < -1/3 for every k >= 4.
  CZ.5 [proved from CZ.2, every k >= KM] |q^(k)| strictly decreasing for k >= KM (two-sided bounds, monotone in k).
Stage 1 untouched; no Psi_info, W, 0.08, 0.23 or 30 microN/kW enters.  Functions on a graph, not a flux.  Not thrust.
"""
from __future__ import annotations

import os
import sys
import time
from fractions import Fraction as Fr

from flint import arb, arb_series, ctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gasket_fractional_profile_level_ratio_k as cy  # noqa: E402

FAILS: list[str] = []
RESULT: dict = {}
PREC = 160
K0 = 24
K_SCAN = 150                       # explicit per-k check range for CZ.4/CZ.5 below the monotone-tail threshold
# Neumann power 24 (CY.5's 40 underflows to nan in CR.2's Frobenius bound at k >= 24) and rho_in = 0.9 |nu_2|
M_POW_CZ, RIN_FRAC_CZ = int(os.environ.get('CZ_MPOW', '24')), float(os.environ.get('CZ_RINF', '0.9'))
REPORT = (24, 25, 26, 30, 40, 50, 60)
# CZ.4/CZ.5 gap: CU's routine (imported; as in CY.5, with the eigenvalue circles of radius GRAD_REL |nu_2| instead of
# 1e-20, Neumann power 24 and rho_in = 0.9 |nu_2|) for k = K0..KM.  Float seeds from CZ.2's first-order formula
# (floats only; Krawczyk certifies rho and nu_2; rho_in's validity is certified by the Neumann bound).
GRAD_REL = float(os.environ.get('CZ_GRAD', '1e-7'))


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


_a = cy._a
upper = cy.upper
up_str = cy.up_str


def ball(x) -> arb:
    return arb(0, upper(abs(arb(x))))


# ----------------------------------------------------------------------------------------- CZ.1 Taylor enclosures
def tail_eps(fac: Fr, steps: int) -> arb:
    """Uniform bound e with |prod_(j >= J) f_j - 1| <= e on D' (CY.1(a): |x_J - s*| <= (1/4)(7/20)^(J - steps) on D')."""
    dJ = _a(cy.DLT) * _a(cy.THETA) ** (cy.J - steps)
    t = _a(Fr(7, 6) * fac / (1 - cy.THETA)) * dJ
    return upper(t.exp() - 1)


def series_products(s: arb_series, d0: arb, ds: arb, steps: int):
    """P0, P1, kappa as order-2 Taylor series at s0 = s[0], s = s0 + ds*x + ... (s0 in D').  Head: J exact factors in
    Taylor arithmetic.  Tail T (holomorphic on D', |T - 1| <= e): T(s0) in 1 + [+-e], |T'(s0)| <= e/d0, d0 = dist(s0, dD');
    as a series in x the tail is 1 + [+-e] + [+-e/d0]|ds| x."""
    p0 = arb_series([1])
    p1 = arb_series([1])
    kp = s - arb(3) / 2
    x = s
    for _ in range(cy.J):
        y = (arb(15) / 4 - x).sqrt()
        p0 *= 3 * (x + arb(1) / 2) / (2 * x * (arb(7) / 2 - x))
        p1 *= 3 / (2 * y)
        kp *= 3 / (y + arb(3) / 2)
        x = y
    out = []
    for p, fac in ((p0, cy.FACT["P0"]), (p1, cy.FACT["P1"]), (kp, cy.FACT["K"])):
        e = tail_eps(fac, steps)
        out.append(p * arb_series([1 + arb(0, e), arb(0, upper(e / d0 * abs(ds)))]))
    return out


def cz1() -> bool:
    ctx.cap = 2
    steps = RESULT["steps"]
    c, R, RP = _a(cy.C0), _a(cy.R0), _a(cy.RP)
    ss = arb(3) / 2
    t = (arb(15) / 4 + ss).sqrt()                 # sigma_2(s*) = sqrt(21)/2
    dss = RP - abs(ss - c)
    dt = RP - abs(t - c)
    ok = bool(dss > 0) and bool(dt > 0)
    s = arb_series([ss, 1])
    y = (arb(15) / 4 + s).sqrt()                  # sigma_2(s*) + x/(2t) + ...
    p0, p1, kps = series_products(s, dss, arb(1), steps)
    q0, q1, _ = series_products(y, dt, 1 / (2 * t), steps)
    Bs = 3 / (2 * y) * q1 / p1
    Gs = 3 * (s - arb(1) / 2) / (2 * s * (s + arb(7) / 2)) * q0 / p0
    _, _, kt = series_products(arb_series([t, 1]), dt, arb(1), steps)
    Bst, b1 = Bs[0], Bs[1] / Bs[0]
    Gst, g1 = Gs[0], Gs[1] / Gs[0]
    kap_t, dkap_t = kt[0], kt[1]
    # consistency: kappa(s*) = 0, kappa'(s*) = 1; B*, G* as in CX.1
    ok &= bool(abs(kps[0]) < arb("1e-25")) and bool(abs(kps[1] - 1) < arb("1e-25"))
    ok &= bool(abs(Bst - arb("0.77737602763")) < arb("1e-11")) and bool(abs(Gst - arb("0.28706524176")) < arb("1e-11"))
    SP = R / (R * R - (t - c) ** 2)
    alpha_p = b1 * kap_t
    gamma_p = g1 * kap_t + dkap_t / (2 * t) + arb(1) / 3
    alpha_s, gamma_s = Gst * alpha_p, Gst * gamma_p
    ok &= bool(alpha_p < 0) and bool(gamma_p > 0) and bool(kap_t > 0)
    for v in (Bst, Gst, b1, g1, kap_t, dkap_t, alpha_s, gamma_s):
        ok &= float(v.rad()) < 1e-25
    RESULT.update(t=t, Bst=Bst, Gst=Gst, b1=b1, g1=g1, kap_t=kap_t, dkap_t=dkap_t, SP=SP,
                  alpha_p=alpha_p, gamma_p=gamma_p, alpha_s=alpha_s, gamma_s=gamma_s)
    check(ok, f"CZ.1: t = sigma_2(3/2) = {t.str(12, radius=False)}; B* = {Bst.str(15, radius=False)}, "
              f"G* = {Gst.str(15, radius=False)}, b1 = B'/B(s*) = {b1.str(15, radius=False)}, g1 = G'/G(s*) = "
              f"{g1.str(15, radius=False)}, kappa(t) = {kap_t.str(15, radius=False)}, kappa'(t) = {dkap_t.str(15, radius=False)}, "
              f"SP = {SP.str(8, radius=False)} (all radii < 1e-25); alpha* = G* b1 kappa(t) = {alpha_s.str(15, radius=False)}, "
              f"gamma* = G* gamma' = {gamma_s.str(15, radius=False)}")
    return ok


# ----------------------------------------------------------------------------------------- CZ.2 remainders
def S_of(r: arb) -> arb:
    """CY.2's S(r) = sum c_n r^(-n-1) (NC terms + CY.2's geometric tail), valid as a bound for |z| >= r, (7/20)^2 < r."""
    cn, dNC = RESULT["cn"], RESULT["dNC"]
    S = arb(0)
    for n in range(cy.NC + 1):
        S += cn[n] / r ** (n + 1)
    e = _a(abs(cy.SST - cy.C0))
    kq = (_a(cy.R0) * dNC / (_a(cy.R0) ** 2 - e * (e + dNC))) ** 2
    th2 = _a(cy.THETA) ** 2
    assert bool(th2 < r)
    return upper(S + kq * r ** (-cy.NC - 1) * (th2 / r) / (1 - th2 / r))


def arc_data() -> None:
    """Per-arc upper bounds on dD (CY.1(b)'s 4096 arcs, maximum modulus): |G|, |B| and |kappa o sigma_2|.  Then
    eta'_k = max_i |G|_i |B|_i^k >= sup_D |a_k| and etaK_k = max_i |G|_i |B|_i^k |kappa o sigma_2|_i >= ||E_k kappa||;
    both are nonincreasing in k (every |B|_i < 1)."""
    aG, aB, aK = [], [], []
    for i in range(cy.NARC):
        z = cy.arc(i, cy.NARC)
        p0, p1, _, _ = cy.orbit_products(z)
        y = cy.sig2(z)
        q0, q1, kq, _ = cy.orbit_products(y)
        aB.append(upper(abs(3 / (2 * y) * q1 / p1)))
        aG.append(upper(abs(3 * cy.wt2(z) * q0 / p0)))
        aK.append(upper(abs(kq)))
    RESULT.update(aG=aG, aB=aB, aK=aK, aBmax=max(aB, key=float), aKmax=max(aK, key=float))


def eta_ref(k: int):
    aG, aB, aK = RESULT["aG"], RESULT["aB"], RESULT["aK"]
    e = ek = arb(0)
    for g, b, kk in zip(aG, aB, aK):
        v = g * b ** k
        e = max(e, upper(v), key=float)
        ek = max(ek, upper(v * kk), key=float)
    return e, ek


def bounds(k: int) -> dict | None:
    """All CZ.2/CZ.3 quantities at integer k >= K0 (upper bounds unless stated).  None if a Neumann condition fails.
    Localisation of lambda_1: r2 = min(CY.3's r_2(k), G* B*^k (k|alpha'| + |gamma' - 1/3| + Rem1^(0))), Rem1^(0) the
    bound computed with CY.3's r_2(k) (a valid first pass); every factor is nonincreasing in k."""
    G, B = RESULT["Gst"], RESULT["Bst"]
    K, beta = _a(RESULT["KM"]), RESULT["beta"]
    gb = G * B ** k
    eta, etaK = eta_ref(k)
    _, _, r2cy = cy.r2(k)
    Lam = k * abs(RESULT["b1"]) + abs(RESULT["g1"]) + RESULT["SP"] / (2 * RESULT["t"])
    eps = eta + K * gb * Lam                                # >= ||Q E_k|_W||
    qek = etaK + K * gb * Lam * K                           # >= ||Q E_k kappa||
    al, ga = RESULT["alpha_p"], RESULT["gamma_p"]
    r2 = r2cy
    for _ in range(2):
        R1 = 1 / (arb(4) / 3 - r2) + 2 * S_of(lowerb(arb(1) / 3 - r2))
        if not bool(R1 * eps < 1):
            return None
        w1 = R1 * qek / (1 - R1 * eps)                      # ||w|| bound
        rem1 = upper(Lam * w1)                              # |rho_1| <= gb * rem1
        r2 = min(r2cy, upper(gb * (k * abs(al) + abs(ga - arb(1) / 3) + rem1)), key=float)
    eps0 = eta + gb
    R0 = beta / (arb(5) / 6) + (1 + beta) * RESULT["S_half"]
    if not bool(R0 * eps0 < 1):
        return None
    w0 = R0 * eps0 / (1 - R0 * eps0)                        # |c_k - 1| = |rho_0'| <= w0
    main = k * al + ga                                      # N_k = gb (main + rhoN), |rhoN| <= remN
    remN = rem1 + w0 / 3
    lam0_hi = 1 + gb * (1 + w0)                             # 1 < lambda_0 <= lam0_hi
    # |q + 1/3 - gb main| <= gb [remN + (|main| + remN)(lam0_hi - 1)]   (q + 1/3 = N_k/lambda_0)
    err = remN + (abs(main) + remN) * (lam0_hi - 1)
    return dict(k=k, gb=upper(gb), eta=eta, r2=r2, Lam=Lam, eps=eps, R1=R1, rem1=upper(rem1), w0=upper(w0),
                main=main, remN=upper(remN), lam0_hi=upper(lam0_hi), err=upper(err),
                Gerr=upper(G * err))


def lowerb(x: arb) -> arb:
    return cy.lower(x)


def cz2() -> bool:
    RESULT["S_half"] = S_of(arb(1) / 2)
    arc_data()
    ok = bool(RESULT["aBmax"] < 1)
    rows = []
    tab = {}
    for k in range(K0, K_SCAN + 1):
        b = bounds(k)
        if b is None:
            ok = False
            break
        tab[k] = b
    RESULT["tab"] = tab
    for k in REPORT:
        b = tab[k]
        rows.append(f"k={k}: R_1 eps_k <= {up_str(b['R1'] * b['eps'], 4)}, |c_k - 1| <= {up_str(b['w0'], 4)}, "
                    f"Rem1_k <= {up_str(b['rem1'], 4)}")
    check(ok, f"CZ.2: per-arc sup_dD |B| <= {up_str(RESULT['aBmax'], 6)} < 1, sup_dD |kappa o sigma_2| <= "
              f"{up_str(RESULT['aKmax'], 6)}, eta'_24 = {up_str(tab[24]['eta'], 5)}; S(1/2) <= {up_str(RESULT['S_half'], 6)}; Neumann conditions R_1 eps_k < 1 and R_0 eps0_k < 1 hold for "
              f"k = {K0}..{K_SCAN} (and beyond: every factor is nonincreasing in k, see CZ.4); lambda_0 = 1 + G* B*^k (1 + rho_0'), "
              f"|rho_0'| <= w0_k; lambda_1 + 1/3 = G* B*^k (k alpha' + gamma' - 1/3) + rho_1, |rho_1| <= G* B*^k Rem1_k:\n      "
              + "\n      ".join(rows))
    return ok


def monotone_factors(k: int) -> bool:
    """For k' >= k: eta'_k, etaK_k, r_2(k) nonincreasing (every |B|_i < 1; CY); B*^k Lam_k nonincreasing iff
    B* (1 + |b1|/Lam_k) < 1; B*^k (k|alpha'| + |gamma' - 1/3|) likewise; Lam_k/k and G* B*^k decrease.  Then eps_k, qek_k,
    the localisation r2, R_1, eps0_k, w0_k, G* B*^k Rem1_k, Rem1_k/k and remN_k/k are nonincreasing for every k' >= k."""
    B, b1 = RESULT["Bst"], abs(RESULT["b1"])
    Lam = k * b1 + abs(RESULT["g1"]) + RESULT["SP"] / (2 * RESULT["t"])
    al, c = abs(RESULT["alpha_p"]), abs(RESULT["gamma_p"] - arb(1) / 3)
    return (bool(B * (1 + b1 / Lam) < 1) and bool(B * (1 + al / (k * al + c)) < 1)
            and bool(abs(RESULT["g1"]) + RESULT["SP"] / (2 * RESULT["t"]) >= 0))


def cz3() -> bool:
    tab = RESULT["tab"]
    ok = monotone_factors(K0)
    rows = []
    for k in REPORT:
        b = tab[k]
        rows.append(f"k={k}: (alpha* k + gamma*) = {(RESULT['alpha_s'] * k + RESULT['gamma_s']).str(8, radius=False)}, "
                    f"|q + 1/3 - (alpha* k + gamma*) B*^k| <= {up_str(b['Gerr'], 4)} B*^k")
    # Err_k/k -> 0: Err_k/k <= (remN_k/k)(1 + ...) with remN_k/k nonincreasing and -> 0 (eps_k -> 0)
    b60, b400 = tab[60], tab[K_SCAN]
    ok &= bool(b400["Gerr"] / K_SCAN < b60["Gerr"] / 60)
    RESULT["err_over_k_scan"] = upper(b400["Gerr"] / K_SCAN)
    check(ok, f"CZ.3: for every k >= {K0}, |q^(k) + 1/3 - (alpha* k + gamma*) B*^k| <= Err_k B*^k with Err_k/k nonincreasing "
              f"and -> 0 (Err_k/k <= {up_str(RESULT['err_over_k_scan'], 3)} at k = {K_SCAN}); hence (q^(k) + 1/3)/(k B*^k) -> "
              f"alpha* = {RESULT['alpha_s'].str(12, radius=False)}:\n      " + "\n      ".join(rows))
    return ok


def cz4_bound() -> bool:
    """Sign from the bounds: N_k/(G* B*^k) <= main + remN < 0.  For k >= K1: main/k = alpha' + gamma'/k <= alpha' + gamma'/K1
    (gamma' > 0), remN/k nonincreasing, so main + remN < 0 at K1 in the form (alpha' + gamma'/K1) + remN_K1/K1 < 0
    implies it for every k >= K1."""
    tab = RESULT["tab"]
    al, ga = RESULT["alpha_p"], RESULT["gamma_p"]
    K1 = None
    for k in range(K0, K_SCAN + 1):
        if bool(al + ga / k + tab[k]["remN"] / k < 0) and monotone_factors(k):
            K1 = k
            break
    ok = K1 is not None
    first_direct = next((k for k in range(K0, K_SCAN + 1) if bool(tab[k]["main"] + tab[k]["remN"] < 0)), None)
    RESULT.update(K1=K1, K1_direct=first_direct)
    check(ok, f"CZ.4 (bounds): q^(k) + 1/3 < 0 for every k >= K1 = {K1} (uniform argument; per-k bound first negative at "
              f"k = {first_direct}): N_k/(G* B*^k) <= k alpha' + gamma' + remN_k < 0")
    return ok


def seeds(k: int):
    """Float seeds (rho, nu_2) from CZ.2: 3^(k+1) rho ~ 1 + G* B*^k, 3^(k+1) nu_2 ~ -1/3 + G* B*^k (k alpha' + gamma' - 1/3)."""
    gb = float(RESULT["Gst"].mid()) * float(RESULT["Bst"].mid()) ** k
    lam0 = 1 + gb
    lam1 = -1 / 3 + gb * (k * float(RESULT["alpha_p"].mid()) + float(RESULT["gamma_p"].mid()) - 1 / 3)
    return repr(lam0 / 3 ** (k + 1)), repr(lam1 / 3 ** (k + 1))


def cz4_cu(ks) -> bool:
    """CU's routine (imported, run unchanged; as in CY.5) for the gap k = K0..K1-1."""
    import gasket_fractional_profile_level_infrate_k as cu    # noqa: E402
    old = (cu.N, cu.cr.M_POW)
    cu.N, cu.cr.M_POW = cy.N_CU, M_POW_CZ
    out = {}
    ok = True
    try:
        for k in ks:
            r, nu = seeds(k)
            rf, nf = float(r), abs(float(nu))
            rin = Fr(1, round(1 / (RIN_FRAC_CZ * nf)))
            cu.CONF[k] = dict(seeds=[r, nu], rho_in=rin, proj=[Fr(1, round(3 / rf)), Fr(1, round(2 / (nf - float(rin))))],
                              cl3=None, qcap="0.36")
            cu.configure(k)
            grad = Fr(1, round(1 / (GRAD_REL * nf)))
            cu.G_RAD = [grad, grad]
            cu.RESULT.clear()
            nf0 = len(cu.FAILS) + len(cu.cr.FAILS)
            t0 = time.time()
            print(f"--- CZ.4, k = {k}: CU.1-CU.3 at N = {cu.N}, m = {cu.cr.M_POW}, rho_in = {rin}", flush=True)
            dc = cu.disk_constants()
            eps, beta = cu.tail_constants(dc, cu.N)
            M = cu.build_matrix(cu.N)
            sp = cu.spectral_split(M)
            good = sp is not None
            rho = nuv = q = None
            if good:
                pairs, T, alpha = sp
                _, rho, nuv, q, good = cu.certify(pairs, T, alpha, eps, beta)
            good = bool(good) and len(cu.FAILS) + len(cu.cr.FAILS) == nf0
            ctx.prec = PREC
            if good:
                good = bool(q + arb(1) / 3 < 0)
                out[k] = dict(rho=rho, nu=nuv, q=q)
            if good:     # consistency with CY.3 and CZ.3 (q + 1/3 within Err_k B*^k of (alpha* k + gamma*) B*^k)
                b = RESULT["tab"][k]
                pred = (RESULT["alpha_s"] * k + RESULT["gamma_s"]) * RESULT["Bst"] ** k
                good = bool(abs(q + arb(1) / 3 - pred) <= b["Gerr"] * RESULT["Bst"] ** k)
            check(good, f"CZ.4 (k = {k}, {time.time() - t0:.0f}s): CU's routine (circles of radius {float(grad):.2e}): two "
                        f"eigenvalues of L_k in |z| > {float(rin):.4e}, real, simple; rho = {rho.str(14) if q is not None else '?'}, "
                        f"q = {q.str(12) if q is not None else '?'}, q + 1/3 < 0, inside CZ.3's band")
            ok &= good
            del cu.CONF[k]
    finally:
        cu.N, cu.cr.M_POW = old
        ctx.prec = PREC
    RESULT["cu"] = out
    return ok


def cz5() -> bool:
    """|q^(k+1) + 1/3| < |q^(k) + 1/3| from two-sided bounds: |q + 1/3| in [gb (|main| - remN)/lam0_hi, gb (|main| + remN)].
    Sufficient: B* m_k (1 + d_(k+1)) lam0_hi_k < 1 - d_k, m_k = |main_(k+1)|/|main_k| = 1 + 1/(k - c) (c = gamma'/|alpha'|),
    d_k = remN_k/|main_k|.  d_k is nonincreasing (remN/k nonincreasing, |main|/k nondecreasing as gamma' > 0), m_k and
    lam0_hi_k decrease, so B* m_k (1 + d_k) lam0_hi_k < 1 - d_k at KM gives the sufficient condition for every k >= KM."""
    tab = RESULT["tab"]
    B = RESULT["Bst"]
    KM = None
    for k in range(K0, K_SCAN):
        b, b2 = tab[k], tab[k + 1]
        mk = abs(b2["main"]) / abs(b["main"])
        dk = b["remN"] / abs(b["main"])
        if bool(b["main"] < 0) and bool(B * mk * (1 + dk) * b["lam0_hi"] < 1 - dk) and monotone_factors(k):
            KM = k
            break
    ok = KM is not None
    RESULT["KM_mono"] = KM
    check(ok, f"CZ.5: |q^(k)| = 1/3 + |q^(k) + 1/3| strictly decreases for every k >= {KM}")
    return ok


# CY.5's certified |q^(23)| (used as stated) and CY.6: |q^(k)| strictly decreasing on 8..23, > 1/3
Q23_CY = ("0.33622604", "8.3e-9")


def cz6() -> bool:
    """Gap closure: CU's q^(k), k = K0..max(K1, KM): sign (CZ.4) and |q^(23)| > |q^(24)| > ... > |q^(KM)| (CZ.5)."""
    out = RESULT["cu"]
    ks = list(range(K0, RESULT["KM_mono"]))
    ok = sorted(out) == ks
    aq = {23: arb(Q23_CY[0]) + arb(0, float(Q23_CY[1]))}
    for k in ks:
        aq[k] = abs(out[k]["q"])
    seq = sorted(aq)
    dec = ok and all(bool(aq[seq[i + 1]] < aq[seq[i]]) for i in range(len(seq) - 1))
    # link to CZ.5: |q^(KM) + 1/3| <= G* B*^KM (|main| + remN) < |q^(ks[-1]) + 1/3| (CU), when ks[-1] = KM - 1
    KM = RESULT["KM_mono"]
    b = RESULT["tab"][KM]
    hiKM = upper(b["gb"] * (abs(b["main"]) + b["remN"]))
    link = ks[-1] == KM - 1 and bool(hiKM < abs(out[ks[-1]]["q"] + arb(1) / 3)) and bool(out[ks[-1]]["q"] + arb(1) / 3 < 0)
    dec = dec and link
    RESULT["hiKM"] = hiKM
    above = ok and all(bool(aq[k] > arb(1) / 3) for k in seq)
    RESULT.update(aq=aq, dec=dec, above=above)
    check(ok and dec and above, f"CZ.6: with CY.5's |q^(23)| (as stated), |q^(k)| strictly decreases for k = 23..{ks[-1]} "
                                f"and stays > 1/3, and |q^({KM}) + 1/3| <= {up_str(hiKM, 5)} (CZ.2 bound) < |q^({ks[-1]}) + 1/3|, "
                                f"so with CZ.5 |q^(k)| strictly decreases for every k >= 23 (with CY.6: every k >= 8): " + ", ".join(f"{k}:{aq[k].str(10, radius=False)}" for k in seq))
    return ok and dec and above


def main() -> int:
    print("Theorem CZ: exact rate k B*^k and sign of q^(k) + 1/3; |c_k - 1| explicit. Not thrust.", flush=True)
    t0 = time.time()
    FAILS.clear()
    RESULT.clear()
    old = ctx.prec
    ctx.prec = PREC
    try:
        cy.FAILS.clear()
        cy.RESULT.clear()
        ok = cy.cy1_disk() and cy.cy1_sups() and cy.cy2() and cy.cy3() and cy.cy4()
        check(ok and cy.RESULT["K0"] == K0, f"CY.1-CY.4 re-run (imported): K0 = {cy.RESULT.get('K0')}")
        RESULT.update({k: cy.RESULT[k] for k in ("cn", "dNC", "beta", "KM", "BM", "GM")})
        RESULT["steps"] = cy.RESULT["entry_steps"]
        ok = cz1() and ok
        ok = cz2() and ok
        ok = cz3() and ok
        ok = cz4_bound() and ok
        ok = cz5() and ok
        if os.environ.get("CZ_SKIP_CU") != "1" and RESULT.get("K1") and RESULT.get("KM_mono"):
            ks = list(range(K0, RESULT["KM_mono"]))
            ok = cz4_cu(ks) and ok
            ok = cz6() and ok
    finally:
        ctx.prec = old
        ctx.cap = 10
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem CZ checks passed; functions on a graph, not a flux; nothing here is thrust")
    return 0


if __name__ == "__main__":
    sys.exit(main())
