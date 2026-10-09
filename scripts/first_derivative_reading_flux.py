#!/usr/bin/env python3
"""Theorem U witness: first-derivative O(3)-covariant readings (non-polynomial allowed).

Family (Assumption U):
    T^{ij} = f(s) delta^{ij} + g(s) d^i Psi d^j Psi,   s = |grad Psi|^2,
with f, g in C^1([0, infinity)). This is the general O(3)-covariant local reading
that depends only on first derivatives of Psi. Affine (f, g) recover Theorem O's
quadratic/trace pieces; f, g may be non-polynomial (e.g. exp(-s)).

On-shell (static massless C6, Lap Psi = 0 in vacuum):
    d_j T^{ij} = (f' + g/2) d_i s + g' (d_j s) E^i E^j.
Conserved for every harmonic field iff g' = 0 and f' + g/2 = 0, i.e.
    T = kappa * delta + c * Q
(Q = Maxwell electrostatic stress / eps0). Then G = 0 on every vacuum surface
around an isolated device. Otherwise G(S) = -int_outside div T, moves with S,
and tends to 0 as the surface recedes. Removing the Maxwell piece leaves either
0 or that same surface-dependent remainder.

Part 0 (exact, Fractions): divergence identity for polynomial f, g of low degree;
conservation criterion on a harmonic polynomial.
Part 1 (numeric FD): same identity for non-polynomial f, g = exp(-s), 1/(1+s).
Part 2 (quadrature): Maxwell G ~ 0; three non-conserved readings (trace,
exp(-s) delta, E⊗E/(1+s)) give three values on three vacuum spheres; exterior
integral match; dilation decay.

Exits nonzero if any identity or recorded figure moves. Nothing here is thrust.
"""
from __future__ import annotations

import math
import random
import sys
from fractions import Fraction

import numpy as np

FAIL = []


def check(ok, msg):
    print(("ok   " if ok else "FAIL ") + msg)
    if not ok:
        FAIL.append(msg)


# ------------------------------------------------- exact polynomial helpers --
def padd(p, q, s=1):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + s * v
    return {k: v for k, v in r.items() if v != 0}


def pmul(p, q):
    r = {}
    for k1, v1 in p.items():
        for k2, v2 in q.items():
            k = tuple(a + b for a, b in zip(k1, k2))
            r[k] = r.get(k, 0) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}


def pscale(p, s):
    return {k: v * s for k, v in p.items() if v * s != 0}


def pdiff(p, i):
    r = {}
    for k, v in p.items():
        if k[i]:
            kk = list(k)
            kk[i] -= 1
            r[tuple(kk)] = r.get(tuple(kk), 0) + v * k[i]
    return {k: v for k, v in r.items() if v != 0}


def peval(p, xyz):
    t = Fraction(0)
    for k, v in p.items():
        m = v
        for a, e in enumerate(k):
            m *= Fraction(xyz[a]) ** e
        t += m
    return t


def lap(p):
    r = {}
    for a in range(3):
        r = padd(r, pdiff(pdiff(p, a), a))
    return r


def grad2(psi):
    r = {}
    for a in range(3):
        g = pdiff(psi, a)
        r = padd(r, pmul(g, g))
    return r


def ppow(p, n):
    if n == 0:
        return {(0, 0, 0): Fraction(1)}
    r = {(0, 0, 0): Fraction(1)}
    for _ in range(n):
        r = pmul(r, p)
    return r


def poly_f_of_s(coeffs, s):
    """f = sum coeffs[k] * s^k."""
    r = {}
    for k, c in enumerate(coeffs):
        if c:
            r = padd(r, pscale(ppow(s, k), c))
    return r


def T_poly(psi, fcoe, gcoe, i, j):
    s = grad2(psi)
    f = poly_f_of_s(fcoe, s)
    g = poly_f_of_s(gcoe, s)
    gi, gj = pdiff(psi, i), pdiff(psi, j)
    t = pscale(pmul(g, pmul(gi, gj)), 1)
    if i == j:
        t = padd(t, f)
    return t


def div_poly_direct(psi, fcoe, gcoe, i):
    d = {}
    for j in range(3):
        d = padd(d, pdiff(T_poly(psi, fcoe, gcoe, i, j), j))
    return d


def div_poly_formula(psi, fcoe, gcoe, i):
    """(f' + g/2) d_i s + g' (d_j s) E^i E^j + g E^i Lap."""
    s = grad2(psi)
    # f' coeffs: k*a_k for s^{k-1}
    fp = [Fraction(k) * fcoe[k] for k in range(1, len(fcoe))]
    gp = [Fraction(k) * gcoe[k] for k in range(1, len(gcoe))]
    f = poly_f_of_s(fcoe, s)
    g = poly_f_of_s(gcoe, s)
    fp_s = poly_f_of_s(fp, s) if fp else {}
    gp_s = poly_f_of_s(gp, s) if gp else {}
    ds_i = pdiff(s, i)
    E_i = pdiff(psi, i)
    L = lap(psi)
    # (f' + g/2) d_i s
    r = pmul(padd(fp_s, pscale(g, Fraction(1, 2))), ds_i)
    # g' (d_j s) E^i E^j
    for j in range(3):
        term = pmul(gp_s, pmul(pdiff(s, j), pmul(E_i, pdiff(psi, j))))
        r = padd(r, term)
    # g E^i Lap
    r = padd(r, pmul(g, pmul(E_i, L)))
    return r


print("=== Part 0: exact divergence identity (polynomial f, g) ===")
rng = random.Random(20261009)
monos = [(x, y, z) for x in range(4) for y in range(4) for z in range(4) if x + y + z <= 3]
TRIALS = 40
n_ok = 0
for _ in range(TRIALS):
    fcoe = [Fraction(rng.randint(-3, 3)) for _ in range(3)]  # a0 + a1 s + a2 s^2
    gcoe = [Fraction(rng.randint(-3, 3)) for _ in range(2)]  # b0 + b1 s
    psi = {m: Fraction(rng.randint(-4, 4)) for m in rng.sample(monos, 6)}
    psi = {k: v for k, v in psi.items() if v}
    ok = True
    for i in range(3):
        if div_poly_direct(psi, fcoe, gcoe, i) != div_poly_formula(psi, fcoe, gcoe, i):
            ok = False
            break
    n_ok += ok
check(n_ok == TRIALS,
      f"d_j T^ij = (f'+g/2) d_i s + g'(d_j s) E^i E^j + g E^i Lap as polynomials: {n_ok}/{TRIALS}")

# Conservation criterion: on harmonic Psi_h, div vanishes for all i iff
# T is kappa delta + c Q, i.e. f(s) = kappa - (c/2) s, g(s) = c (affine).
PSI_H = {(3, 0, 0): 1, (1, 2, 0): -3, (0, 1, 1): 2, (1, 1, 1): 1,
         (4, 0, 0): 1, (2, 2, 0): -6, (0, 4, 0): 1}
check(lap(PSI_H) == {}, "Psi_h of Theorem Q is harmonic")

# Maxwell Q: f = -s/2, g = 1
fQ = [Fraction(0), Fraction(-1, 2)]
gQ = [Fraction(1)]
check(all(div_poly_direct(PSI_H, fQ, gQ, i) == {} for i in range(3)),
      "Maxwell Q (f=-s/2, g=1): div vanishes on Psi_h")
# kappa delta: f = 1, g = 0
check(all(div_poly_direct(PSI_H, [Fraction(1)], [Fraction(0)], i) == {} for i in range(3)),
      "kappa delta (f=1, g=0): div vanishes on Psi_h (flux of delta is 0)")
# Non-conserved affine: f = s, g = 0 (trace)
div_tr = [div_poly_direct(PSI_H, [Fraction(0), Fraction(1)], [Fraction(0)], i) for i in range(3)]
check(any(d != {} for d in div_tr), "trace (f=s, g=0): div does not vanish on Psi_h")
# Non-conserved quadratic: f = s^2, g = 0
div_s2 = [div_poly_direct(PSI_H, [Fraction(0), Fraction(0), Fraction(1)], [Fraction(0)], i)
          for i in range(3)]
check(any(d != {} for d in div_s2), "f=s^2, g=0: div does not vanish on Psi_h")
# Affine conserved family: f = k - (c/2)s, g = c
cons_ok = 0
for k in range(-2, 3):
    for c in range(-2, 3):
        fcoe = [Fraction(k), Fraction(-c, 2)]
        gcoe = [Fraction(c)]
        if all(div_poly_direct(PSI_H, fcoe, gcoe, i) == {} for i in range(3)):
            cons_ok += 1
check(cons_ok == 25, f"all 25 affine Maxwell+kappa readings conserved on Psi_h: {cons_ok}/25")

# Spot-check: non-Maxwell affine f=0,g=1 (pure E⊗E) has nonzero div
PT = (1, Fraction(1, 2), -1)
div_ee = tuple(peval(div_poly_direct(PSI_H, [Fraction(0)], [Fraction(1)], i), PT) for i in range(3))
print(f"     div(E⊗E) on Psi_h at (1,1/2,-1) = {tuple(str(v) for v in div_ee)}")
check(div_ee != (0, 0, 0), "pure E⊗E is not on-shell conserved")


# -------------------------------- Part 1: FD check for non-polynomial f, g --
print("\n=== Part 1: finite-difference divergence for non-polynomial f, g ===")


def num_grad_hess_lap(fn, x, h=1e-5):
    """fn: R^3 -> R. Returns grad, Hessian, Laplacian at x."""
    g = np.zeros(3)
    H = np.zeros((3, 3))
    for i in range(3):
        ep = np.zeros(3)
        ep[i] = h
        g[i] = (fn(x + ep) - fn(x - ep)) / (2 * h)
    for i in range(3):
        for j in range(3):
            ei, ej = np.zeros(3), np.zeros(3)
            ei[i] = h
            ej[j] = h
            H[i, j] = (fn(x + ei + ej) - fn(x + ei - ej) - fn(x - ei + ej) + fn(x - ei - ej)) / (4 * h * h)
    return g, H, np.trace(H)


def T_val(E, f, g):
    s = float(np.dot(E, E))
    return f(s) * np.eye(3) + g(s) * np.outer(E, E)


def div_T_fd(fn, x, f, g, h=1e-5):
    """Finite-difference divergence of T[grad fn] at x."""
    out = np.zeros(3)
    for j in range(3):
        ep = np.zeros(3)
        ep[j] = h
        Ep, _, _ = num_grad_hess_lap(fn, x + ep, h=h * 0.5)
        Em, _, _ = num_grad_hess_lap(fn, x - ep, h=h * 0.5)
        Tp = T_val(Ep, f, g)
        Tm = T_val(Em, f, g)
        out += (Tp[:, j] - Tm[:, j]) / (2 * h)
    return out


def div_T_formula_num(E, H, Lap, f, g, fp, gp):
    s = float(np.dot(E, E))
    ds = 2 * H @ E  # d_i s = 2 E^k H_{ik}
    return (fp(s) + g(s) / 2) * ds + gp(s) * np.dot(ds, E) * E + g(s) * E * Lap


# Harmonic test field: Re(x+iy)^3 = x^3 - 3 x y^2
def harm_fn(x):
    return float(x[0] ** 3 - 3 * x[0] * x[1] ** 2 + 2 * x[1] * x[2])


# Non-polynomial pairs
NP_PAIRS = [
    ("exp-delta", lambda s: math.exp(-s), lambda s: 0.0,
     lambda s: -math.exp(-s), lambda s: 0.0),
    ("soft-E⊗E", lambda s: 0.0, lambda s: 1.0 / (1.0 + s),
     lambda s: 0.0, lambda s: -1.0 / (1.0 + s) ** 2),
    ("exp-Maxwell-ish", lambda s: -0.5 * (1.0 - math.exp(-s)), lambda s: math.exp(-s),
     lambda s: -0.5 * math.exp(-s), lambda s: -math.exp(-s)),
]

pts = [np.array([0.7, -0.3, 0.4]), np.array([1.2, 0.5, -0.8]), np.array([-0.4, 0.9, 0.2])]
for name, f, g, fp, gp in NP_PAIRS:
    max_rel = 0.0
    for x in pts:
        E, H, Lap = num_grad_hess_lap(harm_fn, x)
        # analytic grad/hess of harm_fn for Lap check
        # Lap(x^3-3xy^2+2yz) = 6x - 6x = 0
        fd = div_T_fd(harm_fn, x, f, g)
        form = div_T_formula_num(E, H, Lap, f, g, fp, gp)
        scale = max(1e-8, np.linalg.norm(form), np.linalg.norm(fd))
        max_rel = max(max_rel, np.linalg.norm(fd - form) / scale)
    check(max_rel < 5e-4, f"{name}: FD div matches formula on harmonic field (max rel {max_rel:.2e})")

# exp-Maxwell-ish is NOT conserved (f' + g/2 = -0.5 e^{-s} + 0.5 e^{-s} = 0, but g' != 0)
# Wait: f' + g/2 = -0.5 exp(-s) + 0.5 exp(-s) = 0, but g' = -exp(-s) != 0.
# So the g' term survives. Good non-conserved example that almost looks like Maxwell.
x0 = pts[0]
E, H, Lap = num_grad_hess_lap(harm_fn, x0)
form = div_T_formula_num(E, H, Lap, NP_PAIRS[2][1], NP_PAIRS[2][2], NP_PAIRS[2][3], NP_PAIRS[2][4])
nd = float(np.linalg.norm(form))
check(nd > 1e-3,
      "exp-weighted E⊗E (f=-(1-e^{-s})/2, g=e^{-s}) is not on-shell conserved (|div|=%.4f at test point)" % nd)


# ----------------------------------------------------- Part 2: quadrature ----
print("\n=== Part 2: quadrature on Coulomb cluster (static massless C6) ===")
NT, NP = 128, 256
u, wu = np.polynomial.legendre.leggauss(NT)
phi = 2 * np.pi * np.arange(NP) / NP
U, PH = np.meshgrid(u, phi, indexing="ij")
S = np.sqrt(1 - U ** 2)
NRM = np.stack([S * np.cos(PH), S * np.sin(PH), U], axis=-1)
W = wu[:, None] * np.full(NP, 2 * np.pi / NP)[None, :]

CLUSTER = [(3.0, (-0.6, 0.1, 0.0)), (-1.0, (0.5, 0.35, -0.2)), (2.0, (0.1, -0.45, 0.3))]


def fields(pts):
    g = np.zeros(pts.shape)
    H = np.zeros(pts.shape + (3,))
    for q, c in CLUSTER:
        dd = pts - np.asarray(c)
        r = np.linalg.norm(dd, axis=-1)
        g += (-q / r ** 3)[..., None] * dd
        H += q * (3 * dd[..., :, None] * dd[..., None, :] / r[..., None, None] ** 5
                  - np.eye(3) / r[..., None, None] ** 3)
    return g, H


def flux_and_abs(center, R, f, g):
    pts = np.asarray(center) + R * NRM
    E, H = fields(pts)
    s = np.sum(E * E, axis=-1)
    # T · n = f(s) n + g(s) (E·n) E
    En = np.sum(E * NRM, axis=-1)
    Tn = f(s)[..., None] * NRM + g(s)[..., None] * En[..., None] * E
    dA = R ** 2 * W
    G = np.einsum("ij,ijk->k", dA, Tn)
    abs_x = np.einsum("ij,ij->", dA, np.abs(Tn[..., 0]))
    return G, abs_x, E, H, s


def onshell_div(E, H, f, g, fp, gp):
    """Vector field (f'+g/2) ds + g' (ds·E) E  (Lap=0)."""
    s = np.sum(E * E, axis=-1)
    ds = 2 * np.einsum("...ik,...k->...i", H, E)
    return ((fp(s) + 0.5 * g(s))[..., None] * ds
            + (gp(s) * np.sum(ds * E, axis=-1))[..., None] * E)


NR = 96
tr, wr = np.polynomial.legendre.leggauss(NR)
tr, wr = 0.5 * (tr + 1), 0.5 * wr


def exterior_int(center, R, f, g, fp, gp):
    out = np.zeros(3)
    for t, w in zip(tr, wr):
        r = R / t
        E, H = fields(np.asarray(center) + r * NRM)
        d = onshell_div(E, H, f, g, fp, gp)
        jac = r ** 2 * (R / t ** 2) * w
        out -= jac * np.einsum("ij,ijk->k", W, d)
    return out


def f_id(s):
    return s


def g0(s):
    return 0.0 * s


def fp_id(s):
    return np.ones_like(s)


def gp0(s):
    return 0.0 * s


def f_exp(s):
    return np.exp(-s)


def fp_exp(s):
    return -np.exp(-s)


def g_soft(s):
    return 1.0 / (1.0 + s)


def gp_soft(s):
    return -1.0 / (1.0 + s) ** 2


def f0(s):
    return 0.0 * s


def fp0(s):
    return 0.0 * s


def f_max(s):
    return -0.5 * s


def g_max(s):
    return np.ones_like(s)


def fp_max(s):
    return -0.5 * np.ones_like(s)


READINGS = [
    ("Maxwell Q", f_max, g_max, fp_max, gp0, True),
    ("trace |E|^2 delta", f_id, g0, fp_id, gp0, False),
    ("exp(-s) delta", f_exp, g0, fp_exp, gp0, False),
    ("E⊗E/(1+s)", f0, g_soft, fp0, gp_soft, False),
]

SURF = [("c=0 R=2", (0.0, 0.0, 0.0), 2.0),
        ("c=(0.6,0.2,-0.3) R=2.4", (0.6, 0.2, -0.3), 2.4),
        ("c=(-0.9,0,0.5) R=3", (-0.9, 0.0, 0.5), 3.0)]

# First pass: compute and print; record Maxwell ~0 and non-conserved values
REC = {}
for rn, f, g, fp, gp, conserved in READINGS:
    vals = []
    for name, cen, R in SURF:
        G, abs_x, E, H, s = flux_and_abs(cen, R, f, g)
        ext = exterior_int(cen, R, f, g, fp, gp)
        rel = np.max(np.abs(G - ext)) / max(1e-30, np.max(np.abs(G)), np.max(np.abs(ext)))
        print(f"     {rn} sphere {name}: G = {np.array2string(G, precision=10)}  "
              f"(abs row-x {abs_x:.6f}; vs exterior rel {rel:.1e})")
        if conserved:
            check(np.max(np.abs(G)) < 1e-11 and abs_x > 1,
                  f"{rn} sphere {name}: |G| = {np.max(np.abs(G)):.1e} vs abs row-x {abs_x:.6f}")
        else:
            check(rel < 1e-9, f"{rn} sphere {name}: G = -int_outside div to {rel:.1e} relative")
            vals.append(tuple(float(x) for x in G))
    if not conserved:
        REC[rn] = vals
        gap = min(np.max(np.abs(np.array(vals[a]) - np.array(vals[b])))
                  for a in range(3) for b in range(a + 1, 3))
        check(gap > 0.01, f"{rn}: same device, three spheres, three values (smallest pairwise gap {gap:.4f})")

# Maxwell removal: T = Q + exp(-s) delta; residual flux = flux of exp(-s) delta
print("\n=== Maxwell removal leaves the non-polynomial remainder ===")
for s_idx, (name, cen, R) in enumerate(SURF):
    Gq, _, _, _, _ = flux_and_abs(cen, R, f_max, g_max)
    Ge, _, _, _, _ = flux_and_abs(cen, R, f_exp, g0)
    Gsum, _, _, _, _ = flux_and_abs(cen, R,
                                    lambda s: f_max(s) + f_exp(s),
                                    g_max)
    check(np.max(np.abs(Gsum - (Gq + Ge))) < 1e-12,
          f"linearity of flux on {name}")
    check(np.max(np.abs(Gsum - Ge)) < 1e-11,
          f"after Maxwell removal on {name}: residual = exp(-s) delta flux "
          f"(max |G_Q|={np.max(np.abs(Gq)):.1e})")

# Dilation: non-conserved first-derivative fluxes are O(R^{-2}) -> 0
# (near s=0: exp(-s) = 1 - s + O(s^2) so flux(exp delta) ~ -flux(trace);
#  E⊗E/(1+s) ~ E⊗E = Q + (1/2) trace, Maxwell vanishes, so ~ (1/2) trace).
print("\n=== Dilation: non-conserved first-derivative fluxes -> 0 ===")
C1 = np.array((0.3, -0.2, 0.1))
# Trace flux P on the same dilated spheres supplies the leading coefficient.
trace_scaled = []
exp_scaled = []
soft_scaled = []
mags_trace = []
for R in (4.0, 8.0, 16.0, 32.0, 64.0):
    Gt, _, _, _, _ = flux_and_abs(R * C1, R, f_id, g0)
    Ge, _, _, _, _ = flux_and_abs(R * C1, R, f_exp, g0)
    Gs, _, _, _, _ = flux_and_abs(R * C1, R, f0, g_soft)
    trace_scaled.append(R ** 2 * Gt)
    exp_scaled.append(R ** 2 * Ge)
    soft_scaled.append(R ** 2 * Gs)
    mags_trace.append(float(np.linalg.norm(Gt)))
    print(f"     R={R:>4}: |G_trace|={np.linalg.norm(Gt):.6e}, "
          f"|G_exp|={np.linalg.norm(Ge):.6e}, |G_soft|={np.linalg.norm(Gs):.6e}")
# Errors to the R=64 value should roughly halve when R doubles (O(1/R) remainder).
lim_tr = trace_scaled[-1]
lim_exp = exp_scaled[-1]
lim_soft = soft_scaled[-1]
err_tr = [np.linalg.norm(v - lim_tr) for v in trace_scaled[:-1]]
err_exp = [np.linalg.norm(v - lim_exp) for v in exp_scaled[:-1]]
# exp ~ -trace at leading order
check(np.linalg.norm(lim_tr) > 50 and np.linalg.norm(lim_exp + lim_tr) / np.linalg.norm(lim_tr) < 5e-2,
      f"dilation: R^2 G_exp -> -R^2 G_trace leading "
      f"(||lim_tr||={np.linalg.norm(lim_tr):.3f}, ||lim_exp+lim_tr||/||lim_tr||="
      f"{np.linalg.norm(lim_exp + lim_tr) / np.linalg.norm(lim_tr):.2e})")
check(np.linalg.norm(lim_soft - 0.5 * lim_tr) / np.linalg.norm(lim_tr) < 5e-2,
      f"dilation: R^2 G_soft -> (1/2) R^2 G_trace "
      f"(rel {(np.linalg.norm(lim_soft - 0.5 * lim_tr) / np.linalg.norm(lim_tr)):.2e})")
check(all(mags_trace[k + 1] < 0.55 * mags_trace[k] for k in range(4)),
      "dilation: |G_trace| on R S_1 decreases by about 1/4 each doubling (O(R^{-2}))")
check(mags_trace[-1] < 0.05,
      f"dilation: |G_trace|(R=64) = {mags_trace[-1]:.3e} -> 0")
# Record Richardson for docs
rich_tr = 2 * trace_scaled[-1] - trace_scaled[-2]
print(f"     Richardson 2 s_tr(64)-s_tr(32) = {np.array2string(rich_tr, precision=6)}")
print(f"     lim_exp + lim_tr rel to lim_tr: "
      f"{np.linalg.norm(lim_exp + lim_tr) / np.linalg.norm(lim_tr):.3e}")

# Locked recorded values (README / CLAIM_STATUS)
LOCKED = {
    "trace |E|^2 delta": [
        (-18.2952748241, -6.7816526690, 6.2312529880),
        (-35.0468798202, -11.1667620113, 14.4560275587),
        (4.2397653498, -2.8258079986, -3.1167387541),
    ],
    "exp(-s) delta": [
        (5.5717749687, 2.6920037157, -2.2351171899),
        (10.5241700683, 4.1814004818, -4.7864665372),
        (-3.2630881595, 2.1437905472, 2.4770827176),
    ],
    "E⊗E/(1+s)": [
        (-0.1238875858, -0.3611663762, 0.2007774801),
        (-0.6315002605, -0.6853815368, 0.4967991642),
        (0.9929314854, -0.6340327425, -0.8125889061),
    ],
}
print("\n=== Recorded non-conserved fluxes (lock) ===")
for rn, vals in REC.items():
    for s_idx, (got, exp) in enumerate(zip(vals, LOCKED[rn])):
        err = np.max(np.abs(np.array(got) - np.array(exp)))
        check(err < 1e-8, f"{rn} sphere[{s_idx}]: recorded value (err {err:.1e})")
        print(f"     {rn}[{s_idx}] = {np.array2string(np.array(got), precision=10)}")

if FAIL:
    print(f"\n{len(FAIL)} check(s) failed")
    sys.exit(1)
print("\nall Theorem U checks passed; nothing here is thrust")
