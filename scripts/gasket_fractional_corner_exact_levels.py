#!/usr/bin/env python3
"""Theorem EP: levels 2..6 by exact residue sums, and the full sign pattern at alpha = 0.85 (N = 1).  EO kept "at 0.7-0.85 level 1 only
(levels 2..5 and the number of sign changes there undecided); no closed form for n >= 2".  Every level n is a finite object (G_n has
(3^(n+1)+3)/2 vertices): r_n is a rational function with integer coefficients, so EO.1's partial-fraction route works at every level.

  EP.1 (THEOREM, exact; every alpha in (0, 1) for the integrals, every alpha in (gamma, 1) for Y; AQ (1), BL, DZ.3, EB as stated).
       AQ (1) (u_0 = 3 + s, u_(k+1)(s) = beta(s) u_k(p(s)) - e(s)) gives r_n = B_n/A_n with A_n, B_n in Z[s], computed exactly
       (python-flint fmpz_poly, gcd-reduced); for n = 1..6: deg A_n = 2^n, deg B_n = 2^n - 1, A_n has 2^n simple roots, all real and
       < 0 (certified isolation, arb_fmpz_poly_complex_roots; gcd(A_n, A_n') = 1).  The same holds for the denominator of
       r_n eta_1(p^(n-1)(s)) = B_n/(2 A_n (p^(n-1)(s) + 5)) after gcd reduction (BL: eta_1(z) = 1/(2(z + 5)), EO.1).  For a rational
       N/D with deg N < deg D and simple roots z_j of D off [0, oo), partial fractions N/D = sum_j c_j/(s - z_j), c_j = N(z_j)/D'(z_j),
       and int_0^oo s^-a/(s + w) ds = pi w^-a / sin(pi a) (w > 0; here w = -z_j > 0) give
         J_n = int_0^oo s^-a r_n = pi/sin(pi a) sum_j c_j (-z_j)^-a,   H_n = int_0^oo s^-a r_n eta_1(p^(n-1)(s)) likewise,
       om_n = 1 - delta_1^(n) = H_n/J_n, D'_n = J_n/rho^n (BL at k = n, N = 1; x = 5^n s), and (EB; nu < 0, DZ.3)
         Y_n - Y_oo = rho^(2n) (om_n - om) - om |nu| (rho^n + |nu|),
       decreasing in om = 1 - delta_1 and in |nu|.  n = 1 is EO.1.  Only EF.1's om and |nu| are certified quantities; no eps_n, eta_n,
       L, U.
  EP.2 (THEOREM, certified; arb, 128 bits; EF.1 point cells 4000/32000 per unit, EN's certify_point_en; residue sums at 512 bits):
       alpha = 0.85, N = 1: Y_n - Y_oo < 0 for n = 1, 2, 3 and > 0 for n = 4, 5, 6; EL.1 upper at k0 = 5 (1600 cells/unit) gives
       EJ.2 (a)'s onset 5 (Y_n > Y_oo for every n >= 5; EJ.4: 6 at k0 = 6).  So Y_n - Y_oo changes sign exactly once for all n >= 1,
       between levels 3 and 4 (EO.3: at least once in [1, 6]).
  EP.3 (THEOREM, certified; point cells 1000/8000 per unit): alpha = 0.7 and 0.8: Y_n - Y_oo > 0 for n = 1..5, and EJ.2 (a)'s onset
       (EL.1 upper, k0 in {5, 6}, 1600 cells/unit) is 5 (EJ.4: 6): Y_n > Y_oo for every n >= 1 at 0.7 and 0.8 -- no sign change at all (N = 1).
  EP.4 (OBSERVED / attempts to break): r_n from the level-n gasket Laplacian (numpy, n = 1..3) equals B_n/A_n and every root of A_n is
       minus an eigenvalue of L_n; D'_n exact lies inside EN.1's certified bracket at k0 = n (n = 1..5, 0.85); n = 1 equals EO.1's
       closed form; direct float trapezoid Y_n - Y_oo (EO.4's level_quad_float) lies inside the certified balls; the tail route
       (EJ.2 (b) with EN.1's D_n^hi, EJ.2 (a)) leaves n = 2, 3 and n = 4 undecided at 0.85 (not a failure: one-sided criteria).

Functions on a graph, not a flux.  Nothing here is thrust.
"""
from __future__ import annotations

import sys
import time

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import gasket_fractional_corner_two_sided_rate as ej  # noqa: E402
import gasket_fractional_corner_exact_0905_n3 as el  # noqa: E402
import gasket_fractional_corner_exact_091 as em  # noqa: E402
import gasket_fractional_corner_exact_088 as en  # noqa: E402
import gasket_fractional_corner_level_one as eo  # noqa: E402

FAILS: list[str] = []
NMAX = 6
A85 = "0.85"
GRID85 = (4000, 32000)
SIDE = {"0.7": (1000, 8000), "0.8": (1000, 8000)}
NCELL = 1600
SIGN85 = {1: -1, 2: -1, 3: -1, 4: 1, 5: 1, 6: 1}            # EP.2
EXPECT85 = (3, 4)                                          # EP.2: last negative level, first positive level
ONSET85 = 5                                                # EP.2: EJ.2 (a) onset with EL.1 at k0 = 5
EJ4_ONSET = 6                                              # EJ.4 at 0.7, 0.8, 0.85 (k0 = 6)
EO3 = {"0.7": 6.395764e-02, "0.8": 3.912902e-03, "0.85": -1.189998e-02}   # EO.3's level-1 bound (lower / upper end)
OM_EXACT85 = {1: 0.0907457636, 2: 0.0844691702, 3: 0.0783754360, 4: 0.0733690985, 5: 0.0695975186, 6: 0.0668440463}   # EP.2 (rounded), CI


def check(cond: bool, msg: str) -> None:
    print(("   ok   " if cond else "   FAIL ") + msg, flush=True)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- floats (numpy only; used by the CI test)
def gasket_laplacian(n: int):
    """The level-n gasket graph G_n (lattice coordinates, corners (0,0), (2^n,0), (0,2^n)); returns (L, v), v = e_(0,0) - e_(2^n,0)."""
    m = 2**n
    tris = [((0, 0), (m, 0), (0, m))]
    for _ in range(n):
        nxt = []
        for a, b, c in tris:
            ab = ((a[0] + b[0]) // 2, (a[1] + b[1]) // 2)
            bc = ((b[0] + c[0]) // 2, (b[1] + c[1]) // 2)
            ca = ((c[0] + a[0]) // 2, (c[1] + a[1]) // 2)
            nxt += [(a, ab, ca), (ab, b, bc), (ca, bc, c)]
        tris = nxt
    idx: dict = {}
    for t in tris:
        for p in t:
            idx.setdefault(p, len(idx))
    L = np.zeros((len(idx), len(idx)))
    for a, b, c in tris:
        for p, q in ((a, b), (b, c), (c, a)):
            i, j = idx[p], idx[q]
            L[i, i] += 1
            L[j, j] += 1
            L[i, j] -= 1
            L[j, i] -= 1
    v = np.zeros(len(idx))
    v[idx[(0, 0)]], v[idx[(m, 0)]] = 1.0, -1.0
    return L, v


def r_graph(n: int, s: float) -> float:
    L, v = gasket_laplacian(n)
    return 0.5 * float(v @ np.linalg.solve(L + s * np.eye(len(v)), v))


def r_recursion(n: int, s: float) -> float:
    """AQ (1) in floats: u_0 = 3 + s, u_(k+1)(s) = beta(s) u_k(p(s)) - e(s), r_n = 1/u_n."""
    orb = [s]
    for _ in range(n):
        orb.append(orb[-1] * (orb[-1] + 5))
    u = 3 + orb[n]
    for j in range(n - 1, -1, -1):
        t = orb[j]
        u = (t + 6) / ((t + 2) * (t + 5)) * u - 2 * t / (t + 2)
    return 1 / u


def gap_float(al: float, n: int, nu_abs: float, om: float, om_n: float) -> float:
    """EP.1: Y_n - Y_oo = rho^(2n) (om_n - om) - om |nu| (rho^n + |nu|)."""
    rho = 5**al / 3
    return rho ** (2 * n) * (om_n - om) - om * nu_abs * (rho**n + nu_abs)


# ----------------------------------------------------------------- exact (python-flint)
def rn_exact(n: int):
    """AQ (1) exactly in Z[s]: returns (B_n, A_n), r_n = B_n/A_n, gcd-reduced."""
    from flint import fmpz_poly
    s = fmpz_poly([0, 1])
    p = s * (s + 5)
    A, B = s + 3, fmpz_poly([1])
    for _ in range(n):
        Ap, Bp = A(p), B(p)
        A, B = (s + 6) * Ap - 2 * s * (s + 5) * Bp, (s + 2) * (s + 5) * Bp
        g = A.gcd(B)
        A, B = A // g, B // g
    if A.leading_coefficient() < 0:
        A, B = -A, -B
    return B, A


def h_exact(n: int):
    """The integrand of H_n: r_n eta_1(p^(n-1)(s)) = B_n / (2 A_n (p^(n-1)(s) + 5)), gcd-reduced."""
    from flint import fmpz_poly
    s = fmpz_poly([0, 1])
    B, A = rn_exact(n)
    f = s
    for _ in range(n - 1):
        f = f * (f + 5)
    num, den = B, 2 * A * (f + 5)
    g = num.gcd(den)
    num, den = num // g, den // g
    if den.leading_coefficient() < 0:
        num, den = -num, -den
    return num, den


def mellin_exact(num, den, a_s, prec=512):
    """int_0^oo s^-a num/den ds by EP.1's residue sum; asserts deg num < deg den, simple roots, all real and < 0.
    Returns (value as arb at the caller's precision, roots, all residues > 0)."""
    from flint import acb_poly, arb, ctx
    old = ctx.prec
    ctx.prec = prec
    try:
        info = {"deg": den.degree(), "degnum": num.degree(), "simple": den.gcd(den.derivative()).degree() == 0}
        roots = [r for r, _ in den.complex_roots()]
        mults = [m for _, m in den.complex_roots()]
        info["nroots"] = len(roots)
        info["real_neg"] = all(r.imag == 0 and bool(r.real < 0) for r in roots) and all(m == 1 for m in mults)
        if not (info["simple"] and info["real_neg"] and info["degnum"] < info["deg"] and info["nroots"] == info["deg"]):
            return None, info
        N = acb_poly([int(c) for c in num.coeffs()])
        dD = acb_poly([int(c) for c in den.derivative().coeffs()])
        al = arb(a_s)
        tot = arb(0)
        pos = True
        for r in roots:
            z = r.real
            c = (N(r) / dD(r)).real
            pos = pos and bool(c > 0)
            tot += c * (-z) ** (-al)
        info["res_pos"] = pos
        val = tot * arb.pi() / (arb.pi() * al).sin()
        info["roots"] = [float(r.real.mid()) for r in roots]
    finally:
        ctx.prec = old
    return arb(val), info


def level_exact(a_s, n):
    """EP.1 in arb: J_n, H_n, om_n, D'_n."""
    arb, _ = ej._arb()
    B, A = rn_exact(n)
    J, iJ = mellin_exact(B, A, a_s)
    hn, hd = h_exact(n)
    H, iH = mellin_exact(hn, hd, a_s)
    rho = arb(5) ** arb(a_s) / 3
    ok = J is not None and H is not None
    return {"J": J, "H": H, "om": H / J if ok else None, "D": J / rho**n if ok else None, "iJ": iJ, "iH": iH, "ok": ok,
            "degA": A.degree(), "degB": B.degree()}


def gap_arb(P, Lv, n):
    """EP.1: Y_n - Y_oo as a ball from EF.1's om, |nu| (P) and the exact om_n."""
    om, nab, rho = P["om"][1], -P["nu"], P["rho"]
    return rho ** (2 * n) * (Lv["om"] - om) - om * nab * (rho**n + nab)


def sgn(b) -> int:
    return 1 if bool(b > 0) else (-1 if bool(b < 0) else 0)


def main() -> int:
    t0 = time.time()
    print(f"Theorem EP: levels 2..{NMAX} by exact residue sums; full sign pattern at 0.85 (N = 1). Not thrust. ({em.WORKERS} workers)",
          flush=True)
    arb, ctx = ej._arb()
    ctx.prec = 128
    import sympy as sp
    R, nu, om, omn = sp.symbols("R nu om om_n", real=True)          # R = rho^n
    Yn = R * (R * (omn - om) - (-nu * om))                           # EB: X_n = rho^n (om_n - om), C_delta = -nu om, Y_n = rho^n (X_n - C_delta)
    check(sp.simplify(Yn - nu**2 * om - (R**2 * (omn - om) - om * (-nu) * (R + (-nu)))) == 0,
          "EP.1 (sympy): Y_n - Y_oo = rho^(2n) (om_n - om) - om |nu| (rho^n + |nu|) for nu < 0 (EB, DZ.3)")
    # ---- EP.1: exact rational r_n, structure of the denominators
    for n in range(1, NMAX + 1):
        B, A = rn_exact(n)
        hn, hd = h_exact(n)
        _, iJ = mellin_exact(B, A, "0.85")
        _, iH = mellin_exact(hn, hd, "0.85")
        check(A.degree() == 2**n and B.degree() == 2**n - 1 and iJ["simple"] and iJ["real_neg"] and iJ.get("res_pos", False)
              and iH["simple"] and iH["real_neg"] and hd.degree() == 2**n and hn.degree() < hd.degree(),
              f"EP.1 n = {n}: r_n = B_n/A_n, deg A_n = {A.degree()}, deg B_n = {B.degree()}; A_n and the H_n-denominator (deg "
              f"{hd.degree()}) have simple real negative roots; residues of r_n > 0")
    B1, A1 = rn_exact(1)
    check(list(B1.coeffs()) == [5, 1] and list(A1.coeffs()) == [9, 7, 1], "EP.1 n = 1: r_1 = (s + 5)/(s^2 + 7s + 9) (EO.1)")
    # ---- EP.4 (break): G_n's Laplacian, n = 1..3
    for n in (1, 2, 3):
        B, A = rn_exact(n)
        Lg, v = gasket_laplacian(n)
        ev = np.linalg.eigvalsh(Lg)
        rts = mellin_exact(B, A, "0.85")[1]["roots"]
        bf = [float(c) for c in B.coeffs()]
        af = [float(c) for c in A.coeffs()]
        dif = max(abs(r_graph(n, s) - np.polyval(bf[::-1], s) / np.polyval(af[::-1], s)) for s in (1e-3, 0.3, 1.0, 4.0, 50.0))
        check(len(v) == (3 ** (n + 1) + 3) // 2 and dif < 1e-12 and all(np.min(np.abs(ev + z)) < 1e-8 for z in rts),
              f"EP.4 (break) n = {n}: G_n ({len(v)} vertices) Laplacian gives r_n = B_n/A_n (max diff {dif:.1e}); every root of A_n is "
              f"-(an eigenvalue of L_n)")
    # ---- EP.2 at 0.85
    t1 = time.time()
    P = en.certify_point_en(A85, *GRID85)
    check(P["ok"], f"EF.1 routines valid at alpha = {A85} ({GRID85[0]}/{GRID85[1]})")
    print(f"   alpha = {A85} ({GRID85[0]}/{GRID85[1]}): |nu| in {ej.fmt(-P['nu'], 8)}, 1-delta_1 in {ej.fmt(P['om'][1], 8)}, "
          f"D'_oo in {ej.fmt(P['Dp'], 8)} ({time.time() - t1:.1f}s)", flush=True)
    Q1 = eo.level1_arb(A85)
    lev, gaps = {}, {}
    for n in range(1, NMAX + 1):
        lev[n] = level_exact(A85, n)
        gaps[n] = gap_arb(P, lev[n], n)
        print(f"     n = {n}: om_n = {float(lev[n]['om'].mid()):.10f}, D'_n = {float(lev[n]['D'].mid()):.8f}, Y_n - Y_oo in "
              f"[{float(gaps[n].lower()):+.6e}, {float(gaps[n].upper()):+.6e}]", flush=True)
        check(lev[n]["ok"] and sgn(gaps[n]) == SIGN85[n] and abs(float(lev[n]["om"].mid()) - OM_EXACT85[n]) < 1e-10, f"EP.2 alpha = {A85}, n = {n}: Y_n {'<' if SIGN85[n] < 0 else '>'} Y_oo")
    check(bool(abs(lev[1]["om"] - Q1["om1"]) < arb("1e-25")) and bool(abs(lev[1]["D"] - Q1["D1"]) < arb("1e-25")),
          "EP.4 (break): level 1 by the residue sum equals EO.1's closed form (om_1, D'_1)")
    check(float(gaps[1].upper()) <= EO3["0.85"] + 2e-4, f"EP.2 n = 1 consistent with EO.3 (upper end {EO3['0.85']:+.6e} on 1000/8000)")
    t1 = time.time()
    up = {}
    for k in range(1, 6):
        up.update(el.certify_upper_el(A85, k, NCELL, Ns=(1,)))
    low = em.certify_lower_em(A85, tuple(range(1, 6)), (1,), NCELL)
    on = ej.onset_upper(P, {5: up[5]}, 1)
    print(f"   upper (EL.1, k0 = 5, {NCELL}/unit): U(5) <= {float(up[5]['U'].upper()):.6f}, kappa_hi(1, 5) <= "
          f"{float(up[5]['kappa'][1]['hi'].upper()):.6f}; Phi({on[0]}) >= {float(on[2].lower()):.4g} ({time.time() - t1:.1f}s)", flush=True)
    check(up[5]["cond"] and on[0] == ONSET85, f"EJ.2 (a) at k0 = 5: Y_n > Y_oo for every n >= {ONSET85} at alpha = {A85} (EJ.4: {EJ4_ONSET})")
    neg = [n for n in range(1, NMAX + 1) if sgn(gaps[n]) < 0]
    if neg == [1, 2, 3] and all(sgn(gaps[n]) > 0 for n in range(4, NMAX + 1)) and on[0] <= 5:
        print(f"     -> alpha = {A85}, N = 1: Y_n < Y_oo for n = 1..3, Y_n > Y_oo for every n >= 4: exactly one sign change for all n >= 1, "
              f"between {EXPECT85[0]} and {EXPECT85[1]}", flush=True)
    check(neg == [1, 2, 3] and on[0] <= 5, "EP.2: the full sign pattern at 0.85 (N = 1): exactly one sign change, 3 -> 4")
    # ---- EP.4 (break) at 0.85: EN.1's brackets, tail route, floats
    for n in range(1, 6):
        dlo, dhi = en.dn_bounds(P, n, low["L"][n], up[n]["U"])
        check(bool(dlo < lev[n]["D"]) and bool(lev[n]["D"] < dhi),
              f"EP.4 (break): exact D'_{n} = {float(lev[n]['D'].mid()):.8f} inside EN.1's certified [{float(dlo.lower()):.6f}, "
              f"{float(dhi.upper()):.6f}] (k0 = {n})")
    tail = {}
    for n in (2, 3):
        vs = []
        for k in range(1, n + 1):
            dhi = en.dn_bounds(P, n, low["L"][k], up[k]["U"])[1]
            vs.append(en.lower_side(P, low["L"][k], low["kappa"][(1, k)], n, dhi))
        tail[n] = min(v for v in vs if v is not None)
    on4 = ej.onset_upper(P, {4: up[4]}, 1)
    print(f"     (OBSERVED) tail route at 0.85: EJ.2 (b) with EN.1's D_n^hi, best k0 <= n: n = 2: {tail[2]:+.4g}, n = 3: {tail[3]:+.4g} "
          f"(undecided); EJ.2 (a) at k0 = 4: onset {on4[0]} (level 4 undecided)", flush=True)
    check(tail[2] > 0 and tail[3] > 0 and on4[0] > 4, f"EP.4: the one-sided tail route does not reach n = 2, 3 (EJ.2 (b)) or n = 4 (EJ.2 (a) at k0 = 4: "
          f"first certified level {on4[0]}) at 0.85; only EP.1 decides them (not a failure: one-sided criteria)")
    ef = ej._mods()[0]
    c = ef.critical_floats(0.85, Ns=(1,))
    for n in range(1, NMAX + 1):
        fv = eo.gap_n_float(0.85, n, -c["nu"], c["om"][1])
        omq = eo.level_quad_float(0.85, n)[1]
        check(float(gaps[n].lower()) - 1e-6 < fv < float(gaps[n].upper()) + 1e-6 and abs(omq - float(lev[n]["om"].mid())) < 1e-8,
              f"EP.4 (OBSERVED): n = {n}: trapezoid om_n {omq:.10f} = exact; float Y_n - Y_oo {fv:+.4e} inside the certified ball")
    for n in (3, 4):
        print(f"     margin at n = {n}: distance of the ball from 0 {min(abs(float(gaps[n].lower())), abs(float(gaps[n].upper()))):.3e} "
              f"against its width {float(gaps[n].upper()) - float(gaps[n].lower()):.3e}", flush=True)
    # ---- EP.3 at 0.7, 0.8
    for a_s, grid in SIDE.items():
        t1 = time.time()
        Pa = en.certify_point_en(a_s, *grid)
        check(Pa["ok"], f"EF.1 routines valid at alpha = {a_s} ({grid[0]}/{grid[1]})")
        ga = {n: gap_arb(Pa, level_exact(a_s, n), n) for n in range(1, 6)}
        upa = {}
        for k in (5, 6):
            upa.update(el.certify_upper_el(a_s, k, NCELL, Ns=(1,)))
        ona = ej.onset_upper(Pa, upa, 1)
        print(f"   alpha = {a_s} ({grid[0]}/{grid[1]}): |nu| in {ej.fmt(-Pa['nu'], 8)}, 1-delta_1 in {ej.fmt(Pa['om'][1], 8)}; Y_n - Y_oo >= "
              + ", ".join(f"{float(ga[n].lower()):+.4e}" for n in range(1, 6)) + f" (n = 1..5); onset {ona[0]} (k0 = {ona[1]}) "
              f"({time.time() - t1:.1f}s)", flush=True)
        check(all(sgn(ga[n]) > 0 for n in range(1, 6)) and ona[0] <= 6 and float(ga[1].lower()) >= EO3[a_s] - 2e-4,
              f"EP.3 alpha = {a_s}: Y_n > Y_oo for n = 1..5 and (EJ.2 (a)) every n >= {ona[0]}: Y_n > Y_oo for every n >= 1, no sign change")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all EP checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
