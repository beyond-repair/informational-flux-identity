#!/usr/bin/env python3
"""Theorem BC: the centred (action-reaction) gradient stencil on the frozen 0.45 gasket mesh.

Setting: Assumption BA (README), with the corner gradient replaced by Assumption BC:
  g^GG_x = mean of the gradients of the affine interpolants of u on the level-n cells containing x
  (Green-Gauss vertex average).  On the triangular lattice this IS the centred six-neighbour stencil
  (1/3) sum_e e (u(x+e) - u(x)) of Kept Failure BA.1's torus witness (sympy identity below).

Theorem BC (checked here; exit 1 on any failure):
  (1) reflection law, every S3-equivariant linear gradient rule that kills constants: H_kc = R_ck H_ck
      (R_ck = reflection swapping corners c,k), H_cc = s (P_c - Pbar).  So the Maxwell piece of the pair
      force is exactly its chord-parallel (central) component, the residual is the self force plus the
      chord-normal (non-central) component, and
          kappa_n = -(4/R_n^2) (s + 2 sqrt3 h_perp),   h_perp = H_{c c'} . nu  (nu: unit normal to the chord,
      pointing to the third corner).  Action-reaction (phi = phi^M) <=> s = 0 and h_perp = 0.
  (2) no centred stencil anywhere: at every vertex of V_n (n = 1..6 checked, proof for all n) the
      neighbour offsets lie in a closed half-plane and are not closed under negation; corners lie in one
      cell, every other vertex in exactly two.  The isometry group of V_n is D3, which has no half-turn.
  (3) the Green-Gauss rule coincides with the one-sided corner-cell rule at the three corners, where all
      the charge sits (Theorem H), so kappa^GG_n = kappa_n exactly: certified kappa^GG_1..4(0.45) reproduce
      BA.3, and by Theorem BB 2^-n kappa^GG_n -> K(0.45) = -0.85464690512...
  (4) cell rule components: H_{cc'} = 2^n(-E_n (P_c-Pbar) + Y_n (P_c'-P_c'')), so h_par = 2^(n-1)(E_n+Y_n),
      h_perp = (sqrt3/6) 2^n (E_n - 3Y_n), s = 2^n (sigma_n - E_n); certified non-centrality
      tan theta_n = (E_n - 3Y_n)/(sqrt3 (E_n + Y_n)) < 0 at n = 1..4.
  (5) exact small-s ratio: lim_{s->0} e_n(s)/y_n(s) = 2 * 3^n/(3^n + 1) (telescoping product).
      Conjecture BC.2: E_n/Y_n -> 2 for every alpha, so tan theta_n -> -1/(3 sqrt3).

Not a field solver.  Nothing here is thrust.
"""
from __future__ import annotations

import math
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv, mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import mesh_maxwell_residual as ba  # noqa: E402
import mesh_maxwell_limit as bb  # noqa: E402

FAILS: list[str] = []

# BA.3 certified values (README), used as recorded figures
KAPPA_BA = {1: "-4.0427396046867839149", 2: "-6.6139165202761608886", 3: "-10.738239541819965987",
            4: "-18.526339511158477727"}


def check(cond: bool, msg: str) -> None:
    print(("ok   " if cond else "FAIL ") + msg)
    if not cond:
        FAILS.append(msg)


# ----------------------------------------------------------------- (A) symbolic
def symbolic_checks():
    r3 = sp.sqrt(3)
    P = [sp.Matrix([0, 0]), sp.Matrix([1, 0]), sp.Matrix([sp.Rational(1, 2), r3 / 2])]
    Pb = (P[0] + P[1] + P[2]) / 3

    def refl_fixing(c):
        d = P[c] - Pb
        d = d / sp.sqrt(d.dot(d))
        return sp.simplify(2 * d * d.T - sp.eye(2))

    rot = sp.Matrix([[-sp.Rational(1, 2), -r3 / 2], [r3 / 2, -sp.Rational(1, 2)]])
    ok = all(sp.simplify(rot * (P[c] - Pb) + Pb - P[(c + 1) % 3]) == sp.zeros(2, 1) for c in range(3))
    h1, h2, s, a, b = sp.symbols("h1 h2 s a b", real=True)
    # H from S3 equivariance H_{sigma c, sigma k} = R_sigma H_{ck}, starting from a generic H_{01}
    H = {(0, 1): sp.Matrix([h1, h2])}
    H[(0, 2)] = refl_fixing(0) * H[(0, 1)]
    for ck in [(0, 1), (0, 2)]:
        c, k = ck
        H[((c + 1) % 3, (k + 1) % 3)] = rot * H[ck]
        H[((c + 2) % 3, (k + 2) % 3)] = rot * rot * H[ck]
    for c in range(3):
        H[(c, c)] = s * (P[c] - Pb)          # stabiliser of c is the reflection fixing c
    law = True
    par_ok = True
    for c in range(3):
        for k in range(3):
            if c == k:
                continue
            R = refl_fixing(3 - c - k)       # swaps c and k
            law = law and sp.simplify(H[(k, c)] - R * H[(c, k)]) == sp.zeros(2, 1)
            ch = (P[k] - P[c]) / sp.sqrt((P[k] - P[c]).dot(P[k] - P[c]))
            anti = (H[(c, k)] - H[(k, c)]) / 2
            par_ok = par_ok and sp.simplify(anti - ch * ch.dot(H[(c, k)])) == sp.zeros(2, 1)
    check(ok and law, "S3 equivariance forces H_kc = R_ck H_ck (reflection swapping c,k), H_cc = s(P_c - Pbar)")
    check(par_ok, "antisymmetric (Maxwell) part of H_ck = chord-parallel projection: Maxwell piece = central force")
    x = [a, b, -a - b]
    Q = sum((x[c] ** 2 * (P[c] - Pb) for c in range(3)), sp.zeros(2, 1))
    cross = sum((x[c] * x[k] * H[(c, k)] for c in range(3) for k in range(3) if c != k), sp.zeros(2, 1))
    selfs = sum((x[c] ** 2 * H[(c, c)] for c in range(3)), sp.zeros(2, 1))
    check(sp.simplify(cross - 2 * r3 * h2 * Q) == sp.zeros(2, 1) and sp.simplify(selfs - s * Q) == sp.zeros(2, 1),
          "sum_{c!=k} x_c x_k H_ck = 2 sqrt3 h_perp Q(x), sum_c x_c^2 H_cc = s Q(x)  =>  kappa = -(4/R^2)(s + 2 sqrt3 h_perp)")
    # cell rule: H_{01} = 2^n(-E A_0 + Y B_0)
    E, Y = sp.symbols("E Y", positive=True)
    Hc = -E * (P[0] - Pb) + Y * (P[1] - P[2])
    check(sp.simplify(Hc[0] - (E + Y) / 2) == 0 and sp.simplify(Hc[1] - r3 / 6 * (E - 3 * Y)) == 0,
          "cell rule: h_par = 2^(n-1)(E_n+Y_n), h_perp = (sqrt3/6) 2^n (E_n - 3Y_n); with s = 2^n(sigma-E) this is BB.1")
    check(sp.simplify(-(s + 2 * r3 * Hc[1]) - (3 * Y - E - s)) == 0, "kappa^cross = (4/R^2) 2^n (3Y - E) recovered (BB.1)")
    # Green-Gauss average of the six cells at a triangular-lattice vertex = centred stencil (1/3) sum e du
    Eo = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    vec = lambda p, q: sp.Matrix([p + sp.Rational(1, 2) * q, q * r3 / 2])
    u = sp.symbols("u0:6")
    uc = sp.Symbol("uc")
    avg = sp.zeros(2, 1)
    for i in range(6):
        j = (i + 1) % 6
        M = sp.Matrix([[vec(*Eo[i])[0], vec(*Eo[i])[1]], [vec(*Eo[j])[0], vec(*Eo[j])[1]]])
        avg += M.inv() * sp.Matrix([u[i] - uc, u[j] - uc]) / 6
    cen = sum(((u[i] - uc) * vec(*Eo[i]) for i in range(6)), sp.zeros(2, 1)) / 3
    check(sp.simplify(avg - cen) == sp.zeros(2, 1),
          "triangular lattice: Green-Gauss average of the 6 cell gradients = BA.1's centred stencil (1/3) sum_e e du (odd)")
    # telescoping product at s -> 0: t_m(0+) = 2/(3^m+1), prod_m 3/(3 - t_m) = 2 3^n/(3^n+1)
    tel = True
    t = Fraction(1)
    prod = Fraction(1)
    for m in range(12):
        tel = tel and t == Fraction(2, 3**m + 1)
        prod *= Fraction(3) / (3 - t)
        tel = tel and prod == Fraction(2 * 3 ** (m + 1), 3 ** (m + 1) + 1)
        t = t / (3 - t)                       # t_{m+1}(0) = 2 t_m/(6 - 2 t_m)
    check(tel, "t_m(0+) = 2/(3^m+1) and prod_{m<n} 3/(3-t_m(0+)) = 2*3^n/(3^n+1), n = 1..12 (exact)")
    ss = sp.Symbol("s", positive=True)
    lim_ok = True
    for n in (1, 2, 3):
        r, tt, g, y = bb.rec_all(n, ss)
        lim = sp.limit(sp.simplify(ss * g / y), ss, 0)
        lim_ok = lim_ok and sp.simplify(lim - sp.Rational(2 * 3**n, 3**n + 1)) == 0
    check(lim_ok, "lim_{s->0} e_n(s)/y_n(s) = 2*3^n/(3^n+1) from the BB recursions, n = 1,2,3 (sympy limit)")


def ratio_product_check():
    """e_n/y_n = prod_k (6+z_k)/(6+z_k-2 t_{n-1-k}(z_{k+1})), z_k = p^k(s): exact rationals vs the recursions."""
    ok = True
    for n in (1, 2, 3, 4):
        for s in (Fraction(1, 7), Fraction(2), Fraction(1, 1000)):
            r, t, g, y = bb.rec_all(n, s)
            z = bb.orbit(s, n)
            pr = Fraction(1)
            for k in range(n):
                tk = bb.rec_all(n - 1 - k, z[k + 1])[1] if n - 1 - k > 0 else Fraction(1)
                pr *= (6 + z[k]) / (6 + z[k] - 2 * tk)
            ok = ok and pr == s * g / y
    check(ok, "exact product formula e_n/y_n = prod_k (6+z_k)/(6+z_k-2t_{n-1-k}(z_{k+1})) at n = 1..4, three s")


# ----------------------------------------------------------------- (B) combinatorics of V_n
def cells_of(n):
    P, adj, C, lat = ba.gasket_setup(n)
    h = Fraction(1, 2**n)
    where = {lat[i]: i for i in range(len(lat))}
    cells = []
    for i, (p, q) in enumerate(lat):
        j, k = where.get((p + h, q)), where.get((p, q + h))
        if j is not None and k is not None and j in adj[i] and k in adj[i] and k in adj[j]:
            cells.append((i, j, k))
    return P, adj, C, lat, cells


def combinatorics(levels=range(1, 7)):
    ok_cells = ok_half = ok_odd = True
    for n in levels:
        P, adj, C, lat, cells = cells_of(n)
        N = len(adj)
        cnt = [0] * N
        for cl in cells:
            for v in cl:
                cnt[v] += 1
        ok_cells = ok_cells and len(cells) == 3**n and all(cnt[v] == (1 if v in C else 2) for v in range(N))
        ok_cells = ok_cells and 3 * len(cells) == len(set(tuple(sorted(e)) for i in range(N) for e in [(i, j) for j in adj[i]]))
        for v in range(N):
            off = {(lat[w][0] - lat[v][0], lat[w][1] - lat[v][1]) for w in adj[v]}
            ok_odd = ok_odd and any((-a, -b) not in off for a, b in off)
            # closed half-plane: some direction among the six lattice normals with all offsets on one side
            vecs = [np.array([float(a) + float(b) / 2, float(b) * math.sqrt(3) / 2]) for a, b in off]
            dirs = [np.array([math.cos(t), math.sin(t)]) for t in np.arange(12) * math.pi / 6]
            ok_half = ok_half and any(all(d @ w >= -1e-12 for w in vecs) for d in dirs) and np.linalg.norm(sum(vecs)) > 1e-9
    check(ok_cells, "V_n, n = 1..6: 3^n cells (upward triangles), each edge in one cell; corners in 1 cell, all other vertices in exactly 2")
    check(ok_odd, "V_n, n = 1..6: no vertex has a neighbour set closed under negation (no odd stencil exists anywhere)")
    check(ok_half, "V_n, n = 1..6: at every vertex the neighbour offsets lie in a closed half-plane and do not sum to 0")
    # no half-turn symmetry: point reflection about the midpoint of two corners, or the centroid, is not a symmetry
    P, adj, C, lat, cells = cells_of(3)
    S = {(round(x, 9), round(y, 9)) for x, y in P}
    centres = [(P[C[0]] + P[C[1]]) / 2, P[C].mean(0)] + [P[c] for c in C]
    check(all(any((round(2 * m[0] - x, 9), round(2 * m[1] - y, 9)) not in S for x, y in P) for m in centres),
          "no half-turn (about a corner, a corner-pair midpoint, or the centroid) maps V_3 to itself (isometry group D3)")


# ----------------------------------------------------------------- (C) intervals
def certify_level(n, alpha=Fraction(45, 100), dps=120):
    iv.dps = dps
    P, adj, C, lat, cells = cells_of(n)
    g = ba.Graph(adj)
    mpoly = ba.krylov_minpoly(g)
    K = ba.krylov_columns(g, len(mpoly) - 1)
    ann = ba.verify_annihilates(K, mpoly)
    roots = ba.isolate_roots(mpoly)
    pos = [ba.euclid(t) for t in lat]
    incident = {v: [cl for cl in cells if v in cl] for v in range(len(adj))}

    def gg(u, x):
        """Green-Gauss vertex gradient: mean of cell gradients over the level-n cells containing x."""
        gx = gy = iv.mpf(0)
        cl = incident[x]
        for cell in cl:
            o = [w for w in cell if w != x]
            ev = [(pos[w][0] - pos[x][0], pos[w][1] - pos[x][1]) for w in o]
            a, b = ba.cell_gradient(ev, [u[o[0]] - u[x], u[o[1]] - u[x]])
            gx, gy = gx + a, gy + b
        return (gx / len(cl), gy / len(cl))

    al = ba.ivq(alpha)
    data = [iv.mpf(1), iv.mpf(-0.5), iv.mpf(-0.5)]
    out = ba.maxwell_reading(K, roots, al, C, data, gg)
    cen = [sum((pos[c][s] for c in C), iv.mpf(0)) / 3 for s in (0, 1)]
    Q = tuple(sum(((data[k] - 0) ** 2 * (pos[C[k]][s] - cen[s]) for k in range(3)), iv.mpf(0)) for s in (0, 1))  # dbar = 0
    kap = out["D"][0] / Q[0]
    # pair kernel H and Green entries
    G = ba.matfun(K, ba.interp_coeffs(roots, al, -1))
    N = len(adj)
    Hm = {(c, k): gg([G[j][k] for j in range(N)], c) for c in C for k in C}
    a0, a1, a2 = C
    pb = min(adj[a0], key=lambda j: np.linalg.norm(P[j] - P[a1]))
    Y = G[pb][a1] - G[pb][a2]
    E = G[pb][a1] + G[pb][a2] - G[a0][a1] - G[a0][a2]
    R = G[a0][a0] + G[a1][a1] - 2 * G[a0][a1]
    return dict(ann=ann, D=out["D"], M=out["M"], kap=kap, H=Hm, C=C, pos=pos, cen=cen, Y=Y, E=E, R=R, q=out["q"],
                kself=out["self"][0] / Q[0], incident=incident)


def interval_checks(levels=(1, 2, 3, 4)):
    r3 = iv.sqrt(iv.mpf(3))
    rows = {}
    for n in levels:
        t = time.time()
        r = certify_level(n)
        C, pos, H = r["C"], r["pos"], r["H"]
        a0, a1, a2 = C
        check(r["ann"], f"level {n}: Krylov polynomial annihilates L_n exactly")
        check(all(len(r["incident"][c]) == 1 for c in C),
              f"level {n}: each corner lies in exactly one cell, so g^GG_c is BA's one-sided cell gradient")
        kap = r["kap"]
        rec = iv.mpf(KAPPA_BA[n])
        check(abs(kap - rec).b < 1e-18 and kap.b < 0,
              f"level {n}: Green-Gauss (centred) rule kappa^GG_n in {ba.s(kap, 20)} = BA.3's kappa_n != 0 ({time.time()-t:.1f}s)")
        check(ba.contains0(r["M"][0]) and ba.contains0(r["M"][1]), f"level {n}: Maxwell piece still sums to 0")
        # reflection law H_{a1 a0} = R H_{a0 a1}, chord a0->a1
        ch = [pos[a1][s] - pos[a0][s] for s in (0, 1)]
        ln = iv.sqrt(ch[0] ** 2 + ch[1] ** 2)
        ch = [ch[0] / ln, ch[1] / ln]
        nu = [-ch[1], ch[0]]  # unit normal; a2 lies on this side
        side = (pos[a2][0] - pos[a0][0]) * nu[0] + (pos[a2][1] - pos[a0][1]) * nu[1]
        assert side.a > 0
        h = H[(a0, a1)]
        hpar = h[0] * ch[0] + h[1] * ch[1]
        hperp = h[0] * nu[0] + h[1] * nu[1]
        Rh = [h[s] - 2 * hpar * ch[s] for s in (0, 1)]  # reflection across the perpendicular bisector
        hk = H[(a1, a0)]
        check(ba.contains0(hk[0] - Rh[0]) and ba.contains0(hk[1] - Rh[1]),
              f"level {n}: H_(c'c) = R H_(cc') in intervals (reflection law; not antisymmetry)")
        two = iv.mpf(2) ** n
        check(ba.contains0(hpar - two * (r["E"] + r["Y"]) / 2) and ba.contains0(hperp - two * r3 / 6 * (r["E"] - 3 * r["Y"])),
              f"level {n}: h_par = 2^(n-1)(E+Y), h_perp = (sqrt3/6)2^n(E-3Y)")
        hs = H[(a0, a0)]
        A0 = [pos[a0][s] - r["cen"][s] for s in (0, 1)]
        sfac = hs[0] / A0[0]
        check(ba.contains0(hs[1] - sfac * A0[1]) and sfac.a > 0, f"level {n}: self kernel H_cc = s (P_c - Pbar), s > 0 certified")
        kform = -(4 / r["R"] ** 2) * (sfac + 2 * r3 * hperp)
        check(not ba.disjoint(kform, kap), f"level {n}: kappa_n = -(4/R^2)(s + 2 sqrt3 h_perp) inside the enclosures")
        tan = hperp / hpar
        ratio = r["E"] / r["Y"]
        check(tan.b < 0 and hpar.a > 0, f"level {n}: non-central pair force certified: tan theta_n in {ba.s(tan, 16)} < 0")
        print(f"   level {n}: E_n/Y_n in {ba.s(ratio, 16)}, s = {ba.s(sfac, 12)}, h_par = {ba.s(hpar, 12)}, h_perp = {ba.s(hperp, 12)}")
        rows[n] = (tan, ratio, sfac, hpar, hperp)
    return rows


# ----------------------------------------------------------------- (D) recursion table, Conjecture BC.2
def ratio_table(levels=(1, 4, 10, 20, 30), alphas=("0.45",)):
    out = {}
    for al in alphas:
        for n in levels:
            q = bb.quad_level(n, al)
            out[(al, n)] = q["E"] / q["Y"]
    return out


def main(cert_levels=(1, 2, 3, 4)) -> int:
    t0 = time.time()
    print("Theorem BC: the centred (action-reaction) stencil on the 0.45 gasket mesh. Not thrust.")
    symbolic_checks()
    ratio_product_check()
    combinatorics()
    rows = interval_checks(cert_levels)
    # Kept Failure BC.1 summary
    check(all(rows[n][0].b < 0 for n in rows) and all(rows[n][2].a > 0 for n in rows),
          "Kept Failure BC.1: the centred-stencil surrogate leaves a non-central cross force and a self force at every certified level")
    if 2 in rows and 1 in rows:
        check(rows[2][0].a > rows[1][0].b, "tan theta_2 > tan theta_1 (non-centrality shrinks from level 1 to 2)")
    mp.dps = 30
    tab = ratio_table()
    for (al, n), v in tab.items():
        dev = (2 - v) * mp.mpf(3) ** n
        print(f"   alpha={al} n={n:2d}: E_n/Y_n = {mp.nstr(v, 15)}, (2 - E/Y) 3^n = {mp.nstr(dev, 8)}, "
              f"tan theta = {mp.nstr((v - 3) / (mp.sqrt(3) * (v + 1)), 12)}")
    if 1 in rows:
        check(abs(float(tab[("0.45", 1)]) - float(rows[1][1].mid)) < 1e-9, "recursion quadrature reproduces certified E_1/Y_1")
    if 4 in rows:
        check(abs(float(tab[("0.45", 4)]) - float(rows[4][1].mid)) < 1e-9, "recursion quadrature reproduces certified E_4/Y_4")
    d20 = float((2 - tab[("0.45", 20)]) * 3**20)
    d30 = float((2 - tab[("0.45", 30)]) * 3**30)
    check(abs(float(tab[("0.45", 30)]) - 2) < 1e-12 and abs(d20 - d30) < 1e-4 * abs(d30) and d30 > 0,
          f"Conjecture BC.2 evidence: E_n/Y_n -> 2 from below with (2-E/Y)3^n -> {d30:.6f} (floats, not certified)")
    print(f"   conjectured limit tan theta = -1/(3 sqrt3) = {-1/(3*math.sqrt(3)):.12f}, angle {math.degrees(math.atan(1/(3*math.sqrt(3)))):.6f} deg")
    print(f"elapsed {time.time() - t0:.1f}s")
    if FAILS:
        print(f"{len(FAILS)} FAILED")
        return 1
    print("all Theorem BC checks passed; nothing here is thrust")
    return 0


if __name__ == "__main__":
    lv = tuple(range(1, int(sys.argv[1]) + 1)) if len(sys.argv) > 1 else (1, 2, 3, 4)
    sys.exit(main(cert_levels=lv))
