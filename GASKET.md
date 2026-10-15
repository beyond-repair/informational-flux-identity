# Gasket corner currents

Sequel to the rectangle identity. Same ledger: no selected \(W\), no thrust, no continuum force.

The graph is the combinatorial Sierpinski gasket used in [sierpinski-geometry-045](https://github.com/beyond-repair/sierpinski-geometry-045): \(L = D - A\), levels with \(N(0)=3\), \(N(1)=6\), \(N(n)=(3^{n+1}+3)/2\). Their `net_flux` weights three corner currents by \((P_c - \mathrm{centroid})\). This note computes those currents for a harmonic function and separates them from a net source.

```text
experimental_validation     = false
thrust_validated            = false
energy_extraction_validated = false
```

## Assumptions

1. **ESTABLISHED, this graph.** Undirected edges, symmetric weights \(w_{ij}=w_{ji}\). The combinatorial case is \(w_{ij}=1\) on edges. The outer corners at level \(n\ge 1\) each have two neighbors, both on the boundary edges of the big triangle.
2. **ESTABLISHED, linear algebra.** \((Lu)_i = \sum_{j\sim i} w_{ij}(u_i-u_j)\).
3. **PROPOSED only as their diagnostic.** The audit vector is \(F = \sum_c I_c (P_c - \bar P)\), with \(I_c = \sum_{j\sim c}(u_c-u_j)\) in the combinatorial case. Not a momentum flux of a field theory.
4. **ASSUMPTION, spectral calculus, used in Theorems H, J–L, X, AA–AD.** For the combinatorial Laplacian \(L\succeq 0\) with \(L\mathbf{1}=0\), set \(L^{\alpha}=V\mathrm{diag}(\lambda^{\alpha})V^{T}\) in an eigenbasis of \(L\), with the convention \(0^{\alpha}:=0\) for \(\alpha>0\). This is the same functional calculus as in the related repo's `fractional_laplacian`.
5. **ASSUMPTION, geometric weights, used only in Theorem I.** On the level-\(n\) combinatorial edge set of `build_gasket(n)`, set \(w_{ij}=1/\|P_i-P_j\|^2\). Every edge of that build has Euclidean length \(2^{-n}\), so \(w_{ij}=4^n\) uniformly and \(L_{\mathrm{geom}}=4^n L_{\mathrm{comb}}\).
6. **ASSUMPTION, uniform multiplicative weights, used only in Theorem S.** On the same level-\(n\) combinatorial edge set, set \(w_{ij}=\lambda^n\) for a fixed \(\lambda>0\). The combinatorial case is \(\lambda=1\); geometric \(1/d^2\) is \(\lambda=4\); Kigami resistance renormalization on this finest-edge build is \(\lambda=5/3\).
7. **ASSUMPTION, multi-scale edges, used in Theorems W, X.** Vertices are the level-\(n\) gasket points of `build_gasket(n)`. Edges are the three sides of every triangular cell at every generation \(g=0,\dots,n\) (so coarse edges are retained alongside finer ones). Generation-\(g\) edges get weight \(w=\lambda^g\). The finest-only build of Theorems G–S is the generation-\(n\) slice alone.
8. **CITED THEOREM, used only in Theorem AA.** Löwner–Heinz operator monotonicity: \(0\le M\le L\) and \(0\le\alpha\le1\) imply \(M^{\alpha}\le L^{\alpha}\).
9. **NOT ASSUMED.** A unit source on the interior, any value of \(W\), a nonzero continuum limit of the fractional audit vector, or a continuum kernel \(1/|x-y|^{d+\alpha}\).

## Theorem F. Cut identity

**THEOREM.** For any vertex subset \(S\) and any function \(u\),

\[
\sum_{i\in S}(Lu)_i = \sum_{i\in S,\, j\notin S} w_{ij}(u_i-u_j).
\]

**Proof.** Expand the left-hand side into oriented edge terms. An edge with both ends in \(S\) contributes \(w(u_i-u_j)+w(u_j-u_i)=0\). An edge with one end in \(S\) contributes once. Edges outside \(S\) do not appear.

**Corollary.** Take \(S\) to be every vertex. The cut is empty, so \(\sum_v (Lu)_v = 0\). Split the vertices into the three corners \(C\) and the interior. Then

\[
\sum_{c\in C} I_c = -\sum_{i\notin C}(Lu)_i,
\]

where \(I_c=(Lu)_c\). If \(u\) is harmonic on the interior, the right-hand side is \(0\), so the three corner currents sum to \(0\). Their vector \(F\) can still be nonzero, because it weights the currents by three different vectors. A nonzero \(\|F\|\) is not a net source.

## Theorem G. Harmonic corner currents scale by \(3/5\)

**THEOREM.** Let \(u\) be harmonic on the interior of the combinatorial gasket at level \(n\ge 1\), with corner values \(a,b,c\). Let \(I_a(n), I_b(n), I_c(n)\) be the combinatorial corner currents. Then

\[
I_a(n) = \Bigl(\frac{3}{5}\Bigr)^{n-1} I_a(1),
\]

and the same factor multiplies \(I_b\) and \(I_c\). In particular the triple remains neutral at every level, and \(F(n) = (3/5)^{n-1} F(1)\). Since the three corner positions do not move, \(\|F(n)\|\to 0\) as \(n\to\infty\).

### Level-1 cell

**DERIVED.** On the level-1 graph the midpoints \(m_{ab}, m_{bc}, m_{ca}\) are harmonic, so each has degree 4 and

\[
\begin{aligned}
4m_{ab}-m_{bc}-m_{ca} &= a+b,\\
4m_{bc}-m_{ca}-m_{ab} &= b+c,\\
4m_{ca}-m_{ab}-m_{bc} &= c+a.
\end{aligned}
\]

Adding these gives \(m_{ab}+m_{bc}+m_{ca}=a+b+c\). Substitute into the first equation:

\[
4m_{ab}-(a+b+c-m_{ab})=a+b,
\qquad
5m_{ab}=2a+2b+c.
\]

The other midpoints are the cyclic copies. Therefore

\[
m_{ab}=\frac{2a+2b+c}{5},\quad
m_{bc}=\frac{2b+2c+a}{5},\quad
m_{ca}=\frac{2c+2a+b}{5}.
\]

Each corner has those two midpoints as its only neighbors, so

\[
I_a(1)=(a-m_{ab})+(a-m_{ca}),
\]

and cyclic. For \((a,b,c)=(1,-1/2,0)\) this is exactly

\[
I(1)=\Bigl(\frac{3}{2},\,-\frac{6}{5},\,-\frac{3}{10}\Bigr).
\]

The sum is \(0\). The script checks both the closed form and a direct solve.

### One refinement

**DERIVED.** At level \(n\ge 1\) the two neighbors of an outer corner lie on boundary edges, so each is in exactly one cell. Refining that cell, the new neighbor of corner \(A\) on the edge toward the old neighbor \(p\), with the other old neighbor \(q\), is the midpoint of segment \(A p\) in the cell \(\{A,p,q\}\). The same three-midpoint solve gives

\[
u(q_{Ap})=\frac{2a+2u(p)+u(q)}{5}.
\]

Hence

\[
a-u(q_{Ap})=\frac{3a-2u(p)-u(q)}{5},
\]

and the companion edge contributes \((3a-2u(q)-u(p))/5\). Add them:

\[
I_a(n+1)=\frac{3}{5}\bigl(2a-u(p)-u(q)\bigr)=\frac{3}{5}I_a(n).
\]

The same holds at the other corners. Induction from level 1 is the theorem.

The corner positions are fixed. The barycenter of the vertex set equals the centroid of the three corners at every level, because the point set is invariant under rotation by \(2\pi/3\) about that centroid and no vertex sits on the centroid. The script's build has numerical mean error \(0\) at levels 2, 3, and 4; the identity \(F(n)=(3/5)^{n-1}F(1)\) uses that geometric fact. \(F\) is linear in the current triple, and the triple scales by \(3/5\), so \(F\) does too.

## Explicit dipole for data \((1,-1/2,0)\)

**DERIVED RESULT.** From the level-1 triple and Theorem G,

\[
I_a(n)=\frac{3^n}{2\cdot 5^{n-1}},\qquad
I(n)=I_a(n)\Bigl(1,\,-\frac{4}{5},\,-\frac{1}{5}\Bigr).
\]

With the centroid \(c=(1/2,\sqrt{3}/6)\),

\[
F_x=-\frac{9}{10}I_a(n),\qquad
F_y=-\frac{\sqrt{3}}{10}I_a(n),
\]

so

\[
\|F(n)\|=\Bigl(\frac{3}{5}\Bigr)^n\frac{\sqrt{21}}{2}.
\]

Checked values: \(\|F(2)\|=9\sqrt{21}/50\), and \(\|F(n+1)\|^2/\|F(n)\|^2=9/25\). The script checks the current triple through level 4. The norm formula is the algebra above, not a fit.

So at level 2 the audit vector has length \(9\sqrt{21}/50\approx 0.82486\) while the current sum is the integer \(0\). Refinement multiplies that length by \(3/5\) and never produces a net source.

## Theorem H. Fractional corner currents stay neutral; the \(3/5\) law does not inherit

**ASSUMPTION.** Spectral calculus as in item 4 above. Fix \(\alpha>0\). Write \(I_c^{\alpha}=(L^{\alpha}u)_c\) for the fractional corner currents, and write \(I_c^{L}=(Lu)_c\) for the combinatorial corner currents of the same function \(u\).

**THEOREM H (neutrality).** Let \(u\) satisfy \((L^{\alpha}u)_i=0\) on every interior vertex, with arbitrary corner values. Then

\[
\sum_{c\in C} I_c^{\alpha}=0.
\]

**Proof.** Spectral calculus gives \(L^{\alpha}\mathbf{1}=0\), so \(\mathbf{1}^{T}L^{\alpha}=0\) and \(\sum_v(L^{\alpha}u)_v=0\). Interior terms vanish by the Dirichlet condition, leaving the corner sum.

The cut form of Theorem F does **not** apply to \(L^{\alpha}\): that operator is dense, so \(\sum_{i\in S}(L^{\alpha}u)_i\) is not a boundary cut of edge fluxes. Only the global sum (\(S=\) all vertices) survives, and that is enough for neutrality.

**DERIVED (S3 / Schur).** The gasket automorphisms act as \(S_3\) on the three corners. For fixed level and \(\alpha\), the maps from corner values to the triples \((I_c^{\alpha})\) and \((I_c^{L})\) (both built from the \(L^{\alpha}\)-Dirichlet extension) are linear, kill constants, and intertwine that action. The plane \(x+y+z=0\) is an irreducible \(S_3\)-representation, so by Schur's lemma the two maps are scalar multiples of each other. In particular \(\sum_c I_c^{L}=0\) as well. This is why the related audit's combinatorial currents sum to roundoff even though \(Lu\) is not zero on the interior: Theorem F's interior vanishing is not the reason.

**KEPT FAILURE H.1 (no \(3/5\) inheritance).** The cell-by-cell identity \(I(n+1)=(3/5)I(n)\) of Theorem G uses locality of \(L\) inside one cell. It does not apply to \(L^{\alpha}\). Concretely, for \(\alpha=0.45\) and boundary data \((1,-1/2,0)\), the audit-style norms \(\|F\|\) built from combinatorial currents of the \(L^{0.45}\)-Dirichlet solution at levels \(1,2,3,4\) are

\[
\|F_{\mathrm{comb}}\|\approx(1.439599,\;1.165909,\;1.032272,\;0.956206),
\]

with successive ratios \(\approx(0.810,\;0.885,\;0.926)\), not the constant \(3/5\). The fractional-current norms give ratios \(\approx(0.850,\;0.894,\;0.928)\). The same script's \(\alpha=1\) control recovers ratio \(0.6\) and the exact formula \(\|F(n)\|=(3/5)^n\sqrt{21}/2\). Level 2 matches the related repo's lock \(\|F\|\approx 1.165909\).

So fractional Dirichlet data still produce a **neutral** dipole, and that dipole is **not** governed by the combinatorial refinement factor \(3/5\).

## Theorem I. Geometric \(1/d^2\) currents scale by \(12/5\); \(\|F\|\) diverges

**ASSUMPTION.** Item 5 above: the level-\(n\) edge set is the combinatorial gasket, every edge has length \(2^{-n}\), and \(w_{ij}=1/d^2=4^n\).

**THEOREM I.** Let \(u\) be harmonic on the interior for \(L_{\mathrm{geom}}=4^n L_{\mathrm{comb}}\) with fixed corner values \(a,b,c\). Write \(I_c^{\mathrm{geom}}(n)=(L_{\mathrm{geom}}u)_c\) and \(F_{\mathrm{geom}}(n)\) for the usual audit vector built from those currents. Then:

1. \(u\) coincides with the combinatorial harmonic extension of Theorem G (the scalar \(4^n\) cancels in the interior equations).
2. \(I_c^{\mathrm{geom}}(n)=4^n I_c^{\mathrm{comb}}(n)=4^n\bigl(3/5\bigr)^{n-1}I_c^{\mathrm{comb}}(1)\).
3. The triple stays neutral: \(\sum_c I_c^{\mathrm{geom}}(n)=0\).
4. The refinement ratio is constant:
   \[
   I^{\mathrm{geom}}(n+1)=\frac{12}{5}I^{\mathrm{geom}}(n),\qquad
   F_{\mathrm{geom}}(n+1)=\frac{12}{5}F_{\mathrm{geom}}(n).
   \]
5. For corner data \((1,-1/2,0)\),
   \[
   \|F_{\mathrm{geom}}(n)\|=\Bigl(\frac{12}{5}\Bigr)^n\frac{\sqrt{21}}{2}.
   \]
   In particular \(\|F_{\mathrm{geom}}(n)\|\to+\infty\) as \(n\to\infty\). There is no finite refinement limit.

**Proof.** Uniform edge length \(2^{-n}\) gives \(L_{\mathrm{geom}}=4^n L_{\mathrm{comb}}\). Interior harmonicity is equivalent for the two operators, so the potentials agree. Corner currents therefore pick up exactly the factor \(4^n\). Substitute Theorem G's factor \((3/5)^{n-1}\). The one-step ratio is \(4\cdot(3/5)=12/5\). The norm formula is \(4^n\) times the combinatorial formula \(\|F_{\mathrm{comb}}(n)\|=(3/5)^n\sqrt{21}/2\) already derived for this dipole. Neutrality is Theorem F (or the combinatorial neutrality times \(4^n\)).

**KEPT FAILURE I.1 (no finite nonzero geometric audit limit).** Under the assumption above, geometric weights \(1/d^2\) do **not** produce a finite nonzero refinement limit for \(\|F\|\). The audit norm grows by exactly \(12/5\) at every step. The script checks the closed form through level 5 and the constant ratio \(12/5\). This rejects the reading of the OPEN line that hoped \(1/d^2\) alone would stabilize a finite dipole. It does not speak to \(L^{0.45}\). Theorem J excludes divergence for that operator; zero versus a positive finite limit remains OPEN.

**NOT THIS.** Renormalizing currents by \(4^{-n}\) recovers the combinatorial dipole and its decay to \(0\) (Theorem G). A different edge set (e.g. multi-scale edges) is a different operator; Theorem I does not apply to it. Divergence of \(\|F_{\mathrm{geom}}\|\) is not thrust.

## Theorem S. Uniform multiplicative weights: trichotomy of the harmonic audit limit

**ASSUMPTION.** Item 6 above: on the level-\(n\) combinatorial edge set, \(w_{ij}=\lambda^n\) for fixed \(\lambda>0\). Corner data \((1,-1/2,0)\). The audit vector uses the weighted corner currents \(I_c=(L_\lambda u)_c\).

**THEOREM S.** Let \(u\) be harmonic on the interior for \(L_\lambda=\lambda^n L_{\mathrm{comb}}\). Write \(F_\lambda(n)\) for the audit vector built from those currents. Then:

1. \(u\) coincides with the combinatorial harmonic extension of Theorem G (the scalar \(\lambda^n\) cancels in the interior equations).
2. \(I_c^{(\lambda)}(n)=\lambda^n I_c^{\mathrm{comb}}(n)=\lambda^n(3/5)^{n-1}I_c^{\mathrm{comb}}(1)\).
3. The triple stays neutral, and
   \[
   \|F_\lambda(n)\|
     = \Bigl(\lambda\cdot\frac{3}{5}\Bigr)^n\frac{\sqrt{21}}{2}.
   \]
4. Trichotomy of the refinement limit:
   - if \(\lambda<5/3\), then \(\|F_\lambda(n)\|\to 0\);
   - if \(\lambda=5/3\) (Kigami resistance weight on this build), then \(\|F_\lambda(n)\|=\sqrt{21}/2\) for every \(n\ge 1\);
   - if \(\lambda>5/3\), then \(\|F_\lambda(n)\|\to+\infty\).

In particular Theorems G (\(\lambda=1\)) and I (\(\lambda=4\)) are the two extreme cases of the same formula, and resistance weights are the unique uniform multiplicative choice that produces a finite nonzero harmonic audit limit on this build.

**Proof.** Uniform weight \(\lambda^n\) gives \(L_\lambda=\lambda^n L_{\mathrm{comb}}\). Interior harmonicity is equivalent for the two operators, so the potentials agree. Corner currents pick up exactly \(\lambda^n\). Substitute Theorem G's factor \((3/5)^{n-1}\) and the dipole norm \(\|F_{\mathrm{comb}}(n)\|=(3/5)^n\sqrt{21}/2\). The one-step ratio is \(\lambda\cdot 3/5\), which is \(1\) precisely at \(\lambda=5/3\). Neutrality is Theorem F times the scalar.

**KEPT FAILURE S.1 (no fractional freeze from resistance weights).** The law "resistance weights \(\lambda=5/3\) freeze the \(\alpha=0.45\) fractional audit at a finite nonzero constant by the same cancellation that makes the harmonic audit constantly \(\sqrt{21}/2\)" is false as an inheritance. Spectral calculus gives \(L_{\mathrm{res}}=(5/3)^n L_{\mathrm{comb}}\) and therefore \(L_{\mathrm{res}}^{\alpha}=(5/3)^{n\alpha}L_{\mathrm{comb}}^{\alpha}\), so the \(L_{\mathrm{res}}^{\alpha}\)-Dirichlet solution coincides with the combinatorial fractional solution of Theorems H–L. Weighted corner currents are then exactly \((5/3)^n\) times the combinatorial ones, and
\[
\|F_{\mathrm{res}}(n,0.45)\|
  = \Bigl(\frac{5}{3}\Bigr)^n\|F(n,0.45)\|.
\]
Kept Failure H.1 already says the combinatorial ratios are not \(3/5\), so the product is not level-independent. Concretely the resistance-scaled norms through level \(7\) are strictly increasing
\[
\approx(2.399332,\; 3.238635,\; 4.779036,\; 7.378129,\; 11.698719,\; 18.863790,\; 30.750521),
\]
already past \(30\) at level \(7\), while the harmonic resistance value is the constant \(\sqrt{21}/2\approx 2.291288\). The script checks the closed form of Theorem S through level \(5\) for \(\lambda\in\{1,\,5/3,\,4,\,3/2,\,2\}\), and checks live through level \(4\) that the fractional resistance Dirichlet solution matches the combinatorial one with \(\|F_{\mathrm{res}}\|=(5/3)^n\|F\|\). This rejects that particular route to a finite nonzero \(\alpha=0.45\) limit. It does not decide whether the combinatorial \(\|F(n,0.45)\|\) tends to \(0\) or to a positive finite number; that question remains OPEN.

**NOT THIS.** Multi-scale edge sets (edges retained at every generation) are a different graph; Theorem S does not apply to them. That graph is settled for the harmonic audit in Theorem W. A nonzero resistance-weighted audit norm is not thrust.


## Theorem W. Multi-scale multiplicative weights: harmonic audit is a geometric series

**ASSUMPTION.** Item 7 above: vertices of the level-\(n\) gasket; edges of every generation \(g=0,\dots,n\) with weight \(w=\lambda^g\). Corner data \((1,-1/2,0)\). Write \(L^{(g)}\) for the Laplacian of the generation-\(g\) edges alone (unit weight), so the multi-scale operator is
\[
L_{\mathrm{ms}}=\sum_{g=0}^{n}\lambda^{g}L^{(g)}.
\]
The audit vector uses the multi-scale corner currents \(I_c=(L_{\mathrm{ms}}u)_c\).

**THEOREM W.** Let \(u\) be the combinatorial harmonic extension of Theorem G (finest-edge, unit weight). Then:

1. \(u\) is also harmonic for \(L_{\mathrm{ms}}\) at every \(\lambda>0\) (nested harmonicity: \(L^{(g)}u=0\) on every fine-interior vertex, for each generation \(g\)).
2. The multi-scale corner currents are
   \[
   I_{\mathrm{ms}}(n)=\Bigl(\sum_{k=0}^{n}\bigl(\lambda\cdot\tfrac{3}{5}\bigr)^{k}\Bigr)I(0),
   \]
   where \(I(0)=(5/2,-2,-1/2)\) is the level-0 triangle current of Theorem G.
3. The audit norm is therefore
   \[
   \|F_{\mathrm{ms},\lambda}(n)\|
     = \frac{\sqrt{21}}{2}\sum_{k=0}^{n}\Bigl(\lambda\cdot\frac{3}{5}\Bigr)^{k}.
   \]
4. **Trichotomy.** Write \(r=\lambda\cdot 3/5\).
   - if \(\lambda<5/3\) (\(r<1\)), then \(\|F_{\mathrm{ms},\lambda}(n)\|\to\dfrac{\sqrt{21}}{2(1-r)}=\dfrac{5\sqrt{21}}{2(5-3\lambda)}\) (finite, nonzero);
   - if \(\lambda=5/3\) (\(r=1\)), then \(\|F_{\mathrm{ms},\lambda}(n)\|=(n+1)\sqrt{21}/2\to+\infty\);
   - if \(\lambda>5/3\) (\(r>1\)), then \(\|F_{\mathrm{ms},\lambda}(n)\|\to+\infty\) (successive increments grow by exactly \(r\)).

In particular, multi-scale resistance weights diverge linearly, while multi-scale combinatorial weights (\(\lambda=1\)) freeze at \(5\sqrt{21}/4\). This is the opposite of the finest-only trichotomy of Theorem S, where resistance was the unique finite nonzero freeze.

**Proof.** The restriction of the combinatorial harmonic \(u\) to the generation-\(g\) vertex set is the level-\(g\) harmonic extension (Kigami nested harmonicity / decimation). Hence \(L^{(g)}u=0\) at every vertex that is interior for the fine graph (corners excluded). Summing with coefficients \(\lambda^g\) gives \(L_{\mathrm{ms}}u=0\) on the fine interior, so \(u\) solves the multi-scale Dirichlet problem. Corner currents add by linearity:
\[
I_{\mathrm{ms}}=\sum_{g=0}^{n}\lambda^{g}I^{(g)},
\]
and Theorem G supplies \(I^{(g)}=(3/5)^{g}I(0)\). The audit map is linear in the current triple and sends \(I(0)\) to a vector of norm \(\sqrt{21}/2\). The geometric series and its trichotomy are elementary.

**KEPT FAILURE W.1 (no multi-scale resistance freeze).** The law "resistance weights \(\lambda=5/3\) on the multi-scale edge set freeze the harmonic audit at the finite nonzero constant \(\sqrt{21}/2\) by the same cancellation as Theorem S" is false. Under Theorem W the multi-scale resistance audit is exactly \((n+1)\sqrt{21}/2\), already past \(16\) at level \(6\), while the finest-only resistance value stays \(\sqrt{21}/2\). Multi-scale resistance diverges; the finite freeze on this graph sits at \(\lambda<5/3\) (e.g. \(\lambda=1\) gives \(5\sqrt{21}/4\)). The script checks the closed form through level \(5\) for \(\lambda\in\{1,\,1.2,\,3/2,\,5/3,\,2,\,3,\,4\}\), nested annihilation of each \(L^{(g)}\) on the combinatorial harmonic through level \(4\), linear growth at \(\lambda=5/3\) through level \(6\), and the contrast with Theorem S. This rejects that particular multi-scale route to a finite nonzero harmonic dipole. It does not speak to \(L^{0.45}\) on either graph; zero versus positive at \(\alpha=0.45\) on the combinatorial finest-edge build remains OPEN.

**NOT THIS.** A nonzero multi-scale audit norm is not thrust. Spectral fractional operators on this graph are taken up in Theorem X.

## Theorem X. Multi-scale spectral \(L^{\alpha}\): neutrality; nested harmonicity is special to \(\alpha=1\)

**ASSUMPTION.** Item 7 (multi-scale edges, \(L_{\mathrm{ms}}=\sum_{g=0}^{n}\lambda^{g}L^{(g)}\)) together with the spectral calculus of item 4 applied to \(L_{\mathrm{ms}}\succeq 0\). Corner data \((1,-1/2,0)\). Write \(u_{\mathrm{ms},\alpha}\) for the \(L_{\mathrm{ms}}^{\alpha}\)-Dirichlet minimizer on the fine interior. The multi-scale audit uses the \(L_{\mathrm{ms}}\)-currents \(I_c=(L_{\mathrm{ms}}u)_c\); the fractional-current audit uses \(I_c^{\alpha}=(L_{\mathrm{ms}}^{\alpha}u)_c\).

**THEOREM X.** On every level \(n\ge 1\) and every \(\lambda>0\):

1. Both current triples are neutral (\(\sum I_c=0=\sum I_c^{\alpha}\)) by the cut identity (Theorem F) applied to \(L_{\mathrm{ms}}\) and to \(L_{\mathrm{ms}}^{\alpha}\) (which annihilates constants). With these corner data both triples lie on the Schur line \((1,-4/5,-1/5)\) by the gasket's order-3 symmetry (same argument as Theorems H and J).
2. At \(\alpha=1\), \(u_{\mathrm{ms},1}\) equals the combinatorial harmonic of Theorem G, and the multi-scale audit recovers Theorem W exactly.
3. For \(\alpha\in(0,1)\), nested harmonicity fails: the combinatorial harmonic \(u_H\) satisfies \(L_{\mathrm{ms}}u_H=0\) on the fine interior (Theorem W) but \(L_{\mathrm{ms}}^{\alpha}u_H\not=0\) there, so \(u_{\mathrm{ms},\alpha}\neq u_H\). Nested annihilation of each generation block is special to the linear operator \(\alpha=1\); spectral powers do not inherit it.

**Proof of (1)–(2).** Cut identity and \(L_{\mathrm{ms}}^{\alpha}\mathbf{1}=0\) give neutrality. The \(S_3\) action fixing the centroid cycles the corners and preserves both \(L_{\mathrm{ms}}\) and its functional calculus, so the current triple of any invariant Dirichlet problem with data \((1,-1/2,0)\) is a multiple of \((1,-4/5,-1/5)\). At \(\alpha=1\) the operator is \(L_{\mathrm{ms}}\) itself, so Theorem W applies.

**Proof of (3), witnessed.** If \(L_{\mathrm{ms}}^{\alpha}u_H\) vanished on the fine interior, uniqueness for the \(L_{\mathrm{ms}}^{\alpha}\)-Dirichlet problem would force \(u_{\mathrm{ms},\alpha}=u_H\). The witness records \(\|(L_{\mathrm{ms}}^{0.45}u_H)_{\mathrm{int}}\|_{\infty}\in(0.18,0.45)\) and \(\|u_{\mathrm{ms},0.45}-u_H\|_{\infty}\in(0.09,0.66)\) through levels \(1\)–\(6\) at \(\lambda=1\), so the identification fails. (At \(\alpha=1\) the same residual is at most \(10^{-12}\).)

**KEPT FAILURE X.1 (no multi-scale \(\lambda=1\) freeze for fractional).** The law "because multi-scale \(\lambda=1\) freezes the harmonic audit at the finite nonzero constant \(5\sqrt{21}/4\) (Theorem W), the spectral \(L^{0.45}\) audit on the same multi-scale graph also freezes at that constant" is false. Through levels \(1\)–\(6\) the \(L_{\mathrm{ms}}\)-audit norms of \(u_{\mathrm{ms},0.45}\) at \(\lambda=1\) are
\[
\|F_{\mathrm{ms},L}\|
  \approx(3.749328,\;5.059027,\;6.363924,\;7.716285,\;9.127123,\;10.587688),
\]
strictly increasing, already past \(5\sqrt{21}/4\approx 5.728220\) at level \(3\), and still climbing at level \(6\). The fractional-current norms
\[
\|F_{\mathrm{ms},\alpha}\|
  \approx(1.554895,\;1.754213,\;1.924315,\;2.087851,\;2.250028,\;2.409734)
\]
and the energies \(Q_{0.45}\) likewise increase. By contrast the harmonic multi-scale values at the same levels stay strictly below the freeze and approach it from below. The script also records that multi-scale resistance \(\lambda=5/3\) fractional audits diverge even faster (\(\|F_{\mathrm{ms},L}\|\approx 38.1\) already at level \(5\)). This rejects inheriting the Theorem W freeze (or any multi-scale resistance freeze) for \(L^{0.45}\). It does not decide the finest-edge \(\alpha=0.45\) zero-versus-positive question.

**KEPT FAILURE X.2 (no coincidence with finest-only fractional).** The law "multi-scale \(L_{\mathrm{ms}}^{0.45}\) coincides with finest-edge \(L^{0.45}\) up to a level-independent constant (same audit sequence up to scaling)" is false. On levels \(1\)–\(5\) the finest-edge audit decreases (Theorem J record), while the multi-scale \(L_{\mathrm{ms}}\)-audit increases; their ratio
\[
\|F_{\mathrm{ms},L}\|/\|F_{\mathrm{fin}}\|
  \approx(2.604,\;4.339,\;6.165,\;8.070,\;10.033)
\]
is strictly increasing and already past \(10\) at level \(5\). Pointwise \(\|u_{\mathrm{ms},0.45}-u_{\mathrm{fin},0.45}\|_{\infty}\) grows from \(0.02\) to \(0.19\) on those levels. Opposite monotonicity rules out any constant factor relating the two sequences.

**NOT THIS.** Theorem X does not evaluate a continuum fractional kernel, does not freeze or kill the finest-edge \(\alpha=0.45\) audit, and is not thrust. Zero versus positive at \(\alpha=0.45\) on the combinatorial finest-edge build remains OPEN.

## Theorem J. The fractional audit norm does not diverge

**ASSUMPTION.** Spectral calculus as in item 4. Fix \(\alpha\in(0,1)\). The audit norm \(\|F(n)\|\) is the one in the related repo's `net_flux`: combinatorial corner currents of the \(L^{\alpha}\)-Dirichlet solution, with the vertex-set centroid. Corner data are \((1,-1/2,0)\). This is not the norm built from fractional currents \(I_c^{\alpha}=(L^{\alpha}u)_c\); that companion norm is written \(\|F_{\mathrm{frac}}\|\).

**LEMMA (degree).** For \(n=0\) the graph is a triangle, so every vertex has degree \(2\). For \(n\ge 1\) the graph is three copies of the level-\((n-1)\) graph, glued pairwise at the three edge midpoints and nowhere else. By induction, each outer corner of a copy has degree \(2\) inside that copy, and every other vertex of a copy has degree at most \(4\) inside that copy. A vertex that is not a glued midpoint lies in exactly one copy, so its degree in the union is the degree inside that copy. A glued midpoint is an outer corner of exactly two copies, so its degree is \(2+2=4\). The three outer corners of the big triangle lie in one copy each and have degree \(2\).

**LEMMA (off-diagonal sign).** For \(0<\alpha<1\),
\[
\lambda^{\alpha}
  = \frac{\alpha}{\Gamma(1-\alpha)}
    \int_0^{\infty}(1-e^{-\lambda s})\,s^{-1-\alpha}\,ds.
\]
The identity follows by substituting \(u=\lambda s\) and evaluating \(\int_0^{\infty}(1-e^{-u})u^{-1-\alpha}\,du=\Gamma(1-\alpha)/\alpha\). Therefore
\[
L^{\alpha}
  = \frac{\alpha}{\Gamma(1-\alpha)}
    \int_0^{\infty}(I-e^{-sL})\,s^{-1-\alpha}\,ds.
\]
Let \(M=4I-L\). The degree lemma gives \(M_{ii}\ge 0\) and \(M_{ij}=-L_{ij}\ge 0\), with \(M_{cp}=1\) when \(p\sim c\). Hence \(e^{-sL}=e^{-4s}e^{sM}\) is entrywise nonnegative. The graph is connected and \(M_{ij}>0\) whenever \(i\sim j\), so for every pair \(i,j\) some power \(M^k\) has \((M^k)_{ij}>0\). The series for \(e^{sM}\) therefore has \((e^{sM})_{ij}>0\) for every \(s>0\), and \((e^{-sL})_{ij}>0\). Every off-diagonal entry of \(L^{\alpha}\) is therefore strictly negative, and \(L^{\alpha}\mathbf{1}=0\).

In particular, for each corner neighbor \(p\),
\[
-(L^{\alpha})_{cp}
  \ge \alpha\,4^{\alpha-1},
\]
because \((e^{-sL})_{cp}\ge s e^{-4s}\) and the integral collapses to that constant.

**LEMMA (maximum principle).** A solution of \((L^{\alpha}u)_i=0\) off the three corners attains its maximum and its minimum on the corners. Indeed \((L^{\alpha}u)_i=\sum_{j\neq i} w_{ij}(u_i-u_j)\) with \(w_{ij}=-(L^{\alpha})_{ij}>0\), so a maximum at an interior vertex forces \(u\) constant. For the data \((1,-1/2,0)\), \(-1/2\le u\le 1\) at every vertex.

**THEOREM J.** For every \(n\ge 1\) and every \(\alpha\in(0,1)\),
\[
\|F(n,\alpha)\| \le \frac{3}{5}\sqrt{21}.
\]
In particular \(\|F(n,0.45)\|\) does not tend to \(+\infty\).

**Proof.** Each corner has two neighbors, and each neighbor differs from the corner value by at most \(3/2\), so every combinatorial corner current satisfies \(|I_c|\le 3\). The corner with value \(1\) has nonnegative gaps, so its current \(I_a\) satisfies \(0\le I_a\le 3\). By the Schur argument in Theorem H, the current triple is the scalar multiple of \((1,-4/5,-1/5)\) determined by \(I_a\). The same centroid computation as in the harmonic dipole then gives
\[
\|F\| = I_a\cdot\frac{\sqrt{21}}{5} \le \frac{3}{5}\sqrt{21}.
\]
(The vertex barycenter coincides with the corner centroid by the rotational symmetry already used for Theorem G.)

The bound is uniform in the level. It does not decide whether the \(\alpha=0.45\) limit is \(0\) or a positive finite number.

## Theorem K. Above \(\log 3/\log 5\) the audit norm tends to \(0\)

**THEOREM K.** Let \(\alpha\in(\log 3/\log 5,\, 1)\) and let \(u_H\) be the combinatorial harmonic extension of the same corner data (Theorem G). Write \(Q_{\alpha}(v)=v^{T}L^{\alpha}v\). Then the fractional minimizer \(u\) satisfies
\[
Q_{\alpha}(u)\le Q_{\alpha}(u_H)
  \le \Bigl(\frac{7}{2}\Bigr)^{\alpha} 3^{1-\alpha}\,(3\cdot 5^{-\alpha})^{n},
\]
and the audit norm obeys
\[
\|F(n,\alpha)\|
  \le C(\alpha)\,(3\cdot 5^{-\alpha})^{n/2},
\]
where
\[
C(\alpha)
  = \frac{\sqrt{21}}{5}
    \sqrt{\frac{2}{\alpha\,4^{\alpha-1}}
      \Bigl(\frac{7}{2}\Bigr)^{\alpha} 3^{1-\alpha}}.
\]
Since \(3\cdot 5^{-\alpha}<1\), one has \(\|F(n,\alpha)\|\to 0\) as \(n\to\infty\). The same comparison gives \(\|F_{\mathrm{frac}}(n,\alpha)\|\to 0\).

**Proof.** The minimizer of \(Q_{\alpha}\) with fixed corners is the Dirichlet solution, so \(Q_{\alpha}(u)\le Q_{\alpha}(u_H)\). For the harmonic extension, Theorem G supplies the graph energy
\[
Q_1(u_H)=u_H^{T} L u_H = \frac{7}{2}\Bigl(\frac{3}{5}\Bigr)^{n}.
\]
(The pairing \(Q_1=\sum_c u_c I_c^{H}\) and \(I_a^{H}=(5/2)(3/5)^{n}\) are the dipole formulas from Theorem G.) Let \(\mu\) be the spectral measure of \(u_H\) on the positive eigenspace of \(L\), so \(Q_1=\int\lambda\,d\mu\), \(Q_{\alpha}=\int\lambda^{\alpha}\,d\mu\), and \(\int 1\,d\mu\le\|u_H\|^2\le N(n)\). The maximum principle for \(\alpha=1\) gives \(|u_H|\le 1\), and \(N(n)=(3^{n+1}+3)/2\le 3^{n+1}\). Hölder's inequality with exponents \(1/\alpha\) and \(1/(1-\alpha)\) yields
\[
Q_{\alpha}(u_H)\le Q_1(u_H)^{\alpha}\, N(n)^{1-\alpha}
  \le \Bigl(\frac{7}{2}\Bigr)^{\alpha} 3^{1-\alpha}\,(3\cdot 5^{-\alpha})^{n}.
\]
For the combinatorial current at the corner of value \(1\), the two incident edges are part of the Dirichlet form \(Q_{\alpha}(u)=\sum_{i<j} w_{ij}(u_i-u_j)^2\). With \(w_{cp}\ge \alpha\,4^{\alpha-1}\) and gaps \(g_1,g_2\ge 0\),
\[
Q_{\alpha}(u)
  \ge \alpha\,4^{\alpha-1}\,(g_1^2+g_2^2)
  \ge \frac{\alpha\,4^{\alpha-1}}{2}\, I_a^2,
\]
where the last step is \(g_1^2+g_2^2\ge (g_1+g_2)^2/2\). Therefore \(I_a\le \sqrt{2 Q_{\alpha}(u)/(\alpha 4^{\alpha-1})}\), and multiplying by \(\sqrt{21}/5\) is the audit norm. The hypothesis \(\alpha>\log 3/\log 5\) is exactly \(3\cdot 5^{-\alpha}<1\).

Fractional currents are the Schur multiple of the same vector with pairing \(Q_{\alpha}(u)=(7/5)I_a^{\alpha}\), so \(\|F_{\mathrm{frac}}\|=Q_{\alpha}(u)\sqrt{21}/7\) and the energy bound sends that norm to \(0\) as well.

**NOT THIS.** At \(\alpha=0.45<\log 3/\log 5\) the factor \(3\cdot 5^{-0.45}>1\), so Theorem K does not apply. On that side the harmonic comparison energy is not even a decaying majorant: the witness records \(Q_{0.45}(u_H)\) larger at level \(6\) than at level \(2\). Infinity is still ruled out by Theorem J. Zero versus a positive finite limit remains open at \(\alpha=0.45\).

The threshold \(\log 3/\log 5=d_h/d_w\) is the exponent where this comparison changes regime. The argument does not prove it is sharp. \(\alpha=1\) is already settled by Theorem G, which is stronger than the \(\alpha\to 1\) case of the estimate above; the integral representation was stated for \(\alpha<1\).

## Kept failures at \(\alpha=0.45\)

Audit norms computed for data \((1,-1/2,0)\), \(\alpha=0.45\), levels \(1\) through \(7\) (level \(8\) is recorded below when the witness was run with that argument). Level \(2\) reproduces the lock \(1.165909\).

\[
\begin{align*}
\|F\| &\approx (1.439599390,\; 1.165908607,\; 1.032271686,\; 0.956205578,\\
&\qquad 0.909692382,\; 0.880108995,\; 0.860817772).
\end{align*}
\]

Successive ratios
\[
(0.809884066,\; 0.885379591,\; 0.926311931,\; 0.951356489,\; 0.967479791,\; 0.978080870).
\]

Corner currents stay on the Schur line \((1,-4/5,-1/5)\) and sum to roundoff. The companion fractional-current norms on the same levels are approximately
\[
(1.00281550,\; 0.85240280,\; 0.76184056,\; 0.70681096,\; 0.67258683,\; 0.65073578,\; 0.63647516).
\]

**KEPT FAILURE J.1 (no uniform ratio \(0.95\)).** The law "\(\|F(n+1)\|/\ \|F(n)\|\le 0.95\) for every \(n\), hence \(\|F(n)\|\to 0\) geometrically" is false. The ratio from level \(4\) to level \(5\) is already \(0.95136\), and the ratio from level \(6\) to level \(7\) is \(0.97808\). A ratio bound of \(0.95\) is not available. This does not by itself decide the limit: ratios may tend to \(1\) with either a vanishing or a positive limit.

**KEPT FAILURE J.2 (no \(5^{-0.45}\) contraction of the gaps).** Let \(\delta_n=1-\|F(n)\|/\|F(n-1)\|\) for \(n\ge 2\). The law "\(\delta_{n+1}\le 5^{-0.45}\,\delta_n\) for every \(n\)" would have made \(\sum\delta_n<\infty\) and, with the monotone decrease seen on these levels, would have forced a strictly positive limit. It is false. Here \(5^{-0.45}\approx 0.4847\), while
\[
\frac{\delta_3}{\delta_2},\frac{\delta_4}{\delta_3},\frac{\delta_5}{\delta_4},\frac{\delta_6}{\delta_5},\frac{\delta_7}{\delta_6}
  \approx (0.6029,\; 0.6429,\; 0.6601,\; 0.6685,\; 0.6740).
\]
The contraction has already exceeded \(5^{-0.45}\) by \(\delta_6/\delta_5\). Rejecting the contraction rejects that particular sufficient condition for a positive limit. It does not prove the limit is \(0\).

## Theorem L. Dirichlet energy controls the audit norm for every \(\alpha\in(0,1)\)

**ASSUMPTION.** Spectral calculus as in item 4. The audit norm \(\|F\|\) and the fractional-current norm \(\|F_{\mathrm{frac}}\|\) are as in Theorem J. Write \(Q_{\alpha}(u)=u^{T}L^{\alpha}u\) for the \(L^{\alpha}\)-Dirichlet minimizer \(u\) with corner data \((1,-1/2,0)\), and write \(w_{\min}(\alpha)=\alpha\,4^{\alpha-1}\) for the neighbor-weight lower bound of Theorem J.

**THEOREM L.** For every \(\alpha\in(0,1)\) and every level \(n\ge 1\),

\[
\|F_{\mathrm{frac}}(n,\alpha)\|
  \ge w_{\min}(\alpha)\,\|F(n,\alpha)\|,
\qquad
\|F(n,\alpha)\|
  \le \kappa(\alpha)\,\sqrt{Q_{\alpha}(u_{n})},
\]
where
\[
\kappa(\alpha)
  = \frac{\sqrt{21}}{5}
    \sqrt{\frac{2}{w_{\min}(\alpha)}}.
\]
In particular \(Q_{\alpha}(u_{n})\to 0\) implies \(\|F_{\mathrm{frac}}(n,\alpha)\|\to 0\) and \(\|F(n,\alpha)\|\to 0\). Equivalently, a strictly positive \(\liminf_{n}\|F(n,\alpha)\|\) forces
\[
\liminf_{n} Q_{\alpha}(u_{n})>0
\quad\text{and}\quad
\liminf_{n}\|F_{\mathrm{frac}}(n,\alpha)\|>0.
\]

**Proof.** Off-diagonal signs and the bound \(w_{cp}\ge w_{\min}(\alpha)\) are the integral lemmas of Theorem J. The maximum principle gives \(u\le 1\), so the corner of value \(1\) has nonnegative gaps to every other vertex. Keeping only the two combinatorial neighbors of that corner,
\[
I_{a}^{\alpha}
  = \sum_{j\neq a} w_{aj}(1-u_{j})
  \ge w_{\min}(\alpha)\bigl((1-u_{p})+(1-u_{q})\bigr)
  = w_{\min}(\alpha)\,I_{a}.
\]
Schur's lemma (Theorem H) puts both current triples on the line \((1,-4/5,-1/5)\), and the centroid computation of Theorem J converts the inequality of scalars into \(\|F_{\mathrm{frac}}\|\ge w_{\min}\|F\|\).

For the energy upper bound on \(\|F\|\), the same two edges sit inside the Dirichlet form:
\[
Q_{\alpha}(u)
  \ge w_{\min}(\alpha)\,(g_{1}^{2}+g_{2}^{2})
  \ge \frac{w_{\min}(\alpha)}{2}\,I_{a}^{2},
\]
so \(I_{a}\le\sqrt{2Q_{\alpha}/w_{\min}}\) and \(\|F\|=I_{a}\sqrt{21}/5\). The pairing \(Q_{\alpha}=(7/5)I_{a}^{\alpha}\) (Theorem K) identifies vanishing of \(Q_{\alpha}\) with vanishing of \(\|F_{\mathrm{frac}}\|\), and the weight comparison passes that vanishing to \(\|F\|\).

**NOT THIS.** The comparison does not force \(Q_{0.45}(u_{n})\to 0\). Theorem K's harmonic majorant still grows at \(\alpha=0.45\), and the actual minimizer energy is a different sequence. Zero versus a positive finite audit-norm limit remains OPEN; Theorem L only says that a positive audit-norm limit would need a positive energy liminf, and that energy collapse would settle the limit at \(0\). Theorem AA later proves the energy does not collapse at \(\alpha=0.45\) (Kept Failure AA.1), so that route is closed.

**KEPT FAILURE J.3 (no uniform \(0.95\) decay of \(Q_{\alpha}\)).** The law "\(Q_{0.45}(u_{n+1})/Q_{0.45}(u_{n})\le 0.95\) for every \(n\), hence \(Q\to 0\) and therefore \(\|F\|\to 0\) by Theorem L" is false. The minimizer energies at levels \(1\) through \(6\) are
\[
Q_{0.45}(u_{n})
  \approx(1.531826,\;1.302067,\;1.163731,\;1.079672,\;1.027393,\;0.994015),
\]
with successive ratios
\[
(0.850010,\;0.893757,\;0.927768,\;0.951580,\;0.967512).
\]
The ratio from level \(4\) to level \(5\) is already \(0.95158>0.95\). Rejecting the ratio bound rejects that route to \(\|F\|\to 0\). It does not prove a positive limit.

**KEPT FAILURE J.4 (no \(\log 3/\log 5\) contraction of the gaps).** With \(\delta_n=1-\|F(n)\|/\|F(n-1)\|\) as in Kept Failure J.2, the weaker law "\(\delta_{n+1}/\delta_n<\log 3/\log 5\approx 0.682606\) for every \(n\)" would also have made \(\sum\delta_n<\infty\) and forced a strictly positive limit. It held through level \(9\) (\(\delta_9/\delta_8\approx 0.680954\)). It is false at level \(10\). The audit norms at levels \(8\), \(9\), \(10\) are
\[
\|F(n)\|\approx(0.848024495085,\;0.839442326601,\;0.833639705732),
\]
so \(\delta_{10}/\delta_9\approx 0.683038>\log 3/\log 5\). Equivalently, \(\|F(10)\|\) sits below the kill value \(0.833643371744\) recorded with level \(9\), by \(3.7\cdot 10^{-6}\). Level \(10\) has \(N=88575\) vertices and was computed without a dense eigenbasis (`scripts/gasket_fractional_matrix_free.py`): the \(1581\) eigenpairs with \(\lambda\le 0.05\) applied exactly (eigen-residual \(\le 4.0\cdot 10^{-13}\)), \(\lambda^{0.45}\) on the rest by a degree-\(150\) Chebyshev expansion (scalar error \(9.4\cdot 10^{-14}\)), interior residual \(5.0\cdot 10^{-14}\), corner-current sum \(-4.2\cdot 10^{-13}\). The same code reproduces the dense level-\(8\) norm and an independent level-\(9\) run to twelve digits, and a second level-\(10\) run with a different split (\(941\) exact modes below \(\lambda=0.03\), degree-\(250\) Chebyshev) returns the same \(\|F(10)\|\) to twelve digits. The margin is about seven orders of magnitude above every recorded error. Witness: `scripts/gasket_fractional_level10.py`, which replays the recorded values and runs a live level-\(2\) lock. Rejecting this law closes one more route to a positive limit. It does not prove the limit is \(0\). The minimizer energy also keeps falling, \(Q_{0.45}(u_{9})\approx 0.948091\) and \(Q_{0.45}(u_{10})\approx 0.941537\), which is a record, not energy collapse.

## Theorem AA. Below \(\alpha=1/2\) the fractional energy has a uniform floor

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\); \(u_n\) is the \(L^{\alpha}\)-Dirichlet minimizer at level \(n\), \(Q_{\alpha}(u_n)=u_n^{T}L^{\alpha}u_n\), and \(\|F_{\mathrm{frac}}\|\) is the fractional-current norm of Theorem J. **CITED THEOREM (Löwner–Heinz).** If \(0\le M\le L\) as symmetric matrices and \(0\le\alpha\le 1\), then \(M^{\alpha}\le L^{\alpha}\). This is the standard operator monotonicity of \(t\mapsto t^{\alpha}\), used here as an imported result.

**THEOREM AA.** For every \(\alpha\in(0,1/2)\) and every level \(n\ge 1\),
\[
\frac{9/4}{C(\alpha)}
  \le Q_{\alpha}(u_n)
  \le \frac{7}{6}\,8^{\alpha},
\qquad
C(\alpha)=4^{-\alpha}\Bigl(8+\frac{4}{\pi}\,B\bigl(\tfrac12-\alpha,\tfrac32\bigr)\Bigr),
\]
and therefore, by the pairing \(\|F_{\mathrm{frac}}\|=Q_{\alpha}\sqrt{21}/7\) of Theorem K,
\[
\frac{9\sqrt{21}}{28\,C(\alpha)}
  \le \|F_{\mathrm{frac}}(n,\alpha)\|
  \le \frac{\sqrt{21}}{6}\,8^{\alpha}.
\]
At \(\alpha=0.45\): \(C(0.45)\approx 17.532293\), so \(Q_{0.45}(u_n)\ge 0.128335\) and \(\|F_{\mathrm{frac}}(n,0.45)\|\ge 0.084015\) at every level, with ceilings \(2.973975\) and \(1.946953\). In particular \(\|F_{\mathrm{frac}}(n,0.45)\|\) neither tends to \(0\) nor to \(+\infty\).

**Proof.** *Subgraph comparison.* The side from corner \(a=(0,0)\) to corner \(b=(1,0)\) of the level-\(n\) build is a path \(P_m\) with \(m=2^n\) unit edges. Let \(M\) be the Laplacian of those \(m\) edges alone, extended by zero. Then \(L-M\) is the Laplacian of the remaining edges, so \(0\le M\le L\), and Löwner–Heinz gives \(M^{\alpha}\le L^{\alpha}\). Since \(M\) is block diagonal, \(M^{\alpha}=P_m^{\alpha}\oplus 0\) with the same convention \(0^{\alpha}=0\). Hence \(Q_{\alpha}(u)\ge u_{\mathrm{side}}^{T}P_m^{\alpha}u_{\mathrm{side}}\).

*Duality on the path.* Put \(v=e_0-e_m\perp\mathbf 1\). Cauchy–Schwarz in the \(P_m^{\alpha}\) inner product gives \((u^{T}v)^2\le (u^{T}P_m^{\alpha}u)\,(v^{T}(P_m^{\alpha})^{+}v)\). With \(u_a-u_b=3/2\),
\[
Q_{\alpha}(u)\ge \frac{9/4}{R_{\alpha}(m)},
\qquad
R_{\alpha}(m)=v^{T}(P_m^{\alpha})^{+}v .
\]

*Closed form.* The Neumann path has eigenpairs \(\lambda_k=4\sin^2\theta_k\), \(\theta_k=\pi k/(2(m+1))\), \(\phi_k(j)=\sqrt{2/(m+1)}\cos\bigl(\pi k(2j+1)/(2(m+1))\bigr)\), \(k=1,\dots,m\). Then \(\phi_k(m)=(-1)^k\phi_k(0)\), so only odd \(k\) contribute and
\[
R_{\alpha}(m)=\frac{8}{m+1}\sum_{k\ \mathrm{odd}\le m}\bigl(2\sin\theta_k\bigr)^{-2\alpha}\cos^2\theta_k .
\]

*Uniform bound.* Let \(f(\theta)=(2\sin\theta)^{-2\alpha}\cos^2\theta\), decreasing on \((0,\pi/2)\), and \(\Delta=\pi/(m+1)\). Odd \(k=2j+1\) gives \(\theta_k=(j+\tfrac12)\Delta<\pi/2\), and \(R_{\alpha}(m)=(8/\pi)\sum_j f(\theta_{2j+1})\Delta\). For \(j\ge 1\), \(f(\theta_{2j+1})\Delta\le\int_{(j-1/2)\Delta}^{(j+1/2)\Delta}f\), and these intervals are disjoint inside \((0,\pi/2)\). For \(j=0\), \(\sin x\ge 2x/\pi\) gives \(f(\theta_1)\Delta\le\pi\,4^{-\alpha}(m+1)^{2\alpha-1}\le\pi\,4^{-\alpha}\), using \(2\alpha<1\). The integral converges exactly when \(2\alpha<1\):
\[
\int_0^{\pi/2}f=4^{-\alpha}\,\tfrac12\,B\bigl(\tfrac12-\alpha,\tfrac32\bigr).
\]
So \(R_{\alpha}(m)\le C(\alpha)\) for every \(m\), which is the lower bound. (The Riemann sums also converge, so \(R_{\alpha}(m)\to 4^{-\alpha}(4/\pi)B(\tfrac12-\alpha,\tfrac32)\approx 13.245199\) at \(\alpha=0.45\).)

*Ceiling.* The competitor \(w\) with corner data \((1,-1/2,0)\) and value \(1/6\) at every other vertex differs from a constant by \(z=\tfrac56e_a-\tfrac23e_b-\tfrac16e_c\). Gershgorin with maximal degree \(4\) gives \(\lambda_{\max}(L)\le 8\), so \(Q_{\alpha}(u_n)\le Q_{\alpha}(w)=z^{T}L^{\alpha}z\le 8^{\alpha}\|z\|^2=\tfrac76 8^{\alpha}\), because \(u_n\) is the Dirichlet minimizer. The pairing \(Q_{\alpha}=(7/5)I_a^{\alpha}\) and \(\|F_{\mathrm{frac}}\|=I_a^{\alpha}\sqrt{21}/5\) convert both bounds.

**Witness, computed by `scripts/gasket_fractional_energy_floor.py`.** The closed-form \(R_{\alpha}(2^n)\) matches a dense pseudo-inverse for \(n=0,\dots,7\) and \(\alpha\in\{0.10,0.25,0.45,0.49,0.9\}\) to \(1.3\cdot10^{-11}\). For \(\alpha\in\{0.10,0.25,0.45,0.49\}\), \(R_{\alpha}(2^n)\) increases through \(n=22\) and stays below both the Riemann-sum limit and \(C(\alpha)\); at \(\alpha=0.45\), \(R(2^{22})\approx 10.525442\). Live on levels \(1\)–\(6\) at \(\alpha=0.45\): the side is a path of \(2^n\) gasket edges, \(\lambda_{\min}(L^{\alpha}-M^{\alpha})\ge -5.1\cdot10^{-16}\), and
\[
Q_{0.45}(u_n)\approx(1.531826,\;1.302067,\;1.163731,\;1.079672,\;1.027393,\;0.994015)
\]
dominates the path energy, which dominates \((9/4)/R(2^n)\approx(1.125000,\;0.868689,\;0.701170,\;0.587393,\;0.507267,\;0.448821)\), which dominates the floor \(0.128335\). \(\|F_{\mathrm{frac}}\|=Q\sqrt{21}/7\) to \(10^{-10}\), and \(\lambda_{\max}(L)\le 6\). The recorded level-\(9\) and level-\(10\) energies \(0.948091\) and \(0.941537\) lie inside \([0.128335,\,2.973975]\). Level \(4\) also clears the floor at \(\alpha\in\{0.10,0.25,0.49\}\). Control: at \(\alpha=1/2\), \(R(2^n)\) gains about \(3.53\) every four levels from \(n=10\) to \(n=22\) (logarithmic growth), so this route needs \(\alpha<1/2\). The script exits nonzero if any of these checks fails.

**KEPT FAILURE AA.1 (no energy collapse at \(\alpha=0.45\)).** Theorem L left one route to \(\|F(n,0.45)\|\to 0\): energy collapse \(Q_{0.45}(u_n)\to 0\). Theorem AA proves that route false. \(Q_{0.45}(u_n)\ge 0.128335\) at every level, so no subsequence collapses. The same floor makes "\(\|F_{\mathrm{frac}}(n,0.45)\|\to 0\)" false. Kept Failure J.3 rejected one ratio law for \(Q\). AA.1 rejects every collapse. If the combinatorial audit norm does tend to \(0\), the vanishing must come from the two-neighbor share \(I_a/I_a^{\alpha}\) of the fractional current, with the current itself bounded below.

**NOT THIS.** Theorem AA does not bound the combinatorial audit norm \(\|F(n,0.45)\|\) from below. Theorem L gives only \(\|F_{\mathrm{frac}}\|\ge w_{\min}\|F\|\), which is the wrong direction. Zero versus positive for \(\|F(n,0.45)\|\) (Conjecture J.1) remains OPEN. The window \(\alpha\in[1/2,\log 3/\log 5]\) is open for \(Q_{\alpha}\) and \(\|F_{\mathrm{frac}}\|\) as well, because the path has spectral dimension \(1\) where the gasket has \(2\log 3/\log 5\). A positive fractional energy or current is not a force and not thrust.

## Theorem AB. The combinatorial audit limit is the two-neighbor fractional share

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\); \(u_n\) is the \(L^{\alpha}\)-Dirichlet minimizer at level \(n\). Write \(I_a=(Lu_n)_a\) for the combinatorial corner current at the value-\(1\) corner, \(I_a^{\alpha}=(L^{\alpha}u_n)_a\) for the fractional corner current, \(\sigma_n=I_a/I_a^{\alpha}=\|F\|/\|F_{\mathrm{frac}}\|\), and \(w_{\min}(\alpha)=\alpha\,4^{\alpha-1}\) as in Theorem L. Let \(p,q\) be the two combinatorial neighbors of that corner.

**THEOREM AB.** For every \(\alpha\in(0,1)\) and every level \(n\ge 1\):

1. **Equal neighbor weights.** The altitude reflection through the value-\(1\) corner is a graph automorphism, so \(-(L^{\alpha})_{ap}=-(L^{\alpha})_{aq}=:w_n(\alpha)\).
2. **Factorization.** Write \(I_{\mathrm{loc}}=w_n I_a\) for the two-neighbor piece of the fractional current and \(\rho_n=I_{\mathrm{loc}}/I_a^{\alpha}\in(0,1]\) for its share. Then
   \[
   \sigma_n = \frac{\rho_n}{w_n}.
   \]
3. **Sandwich.** The off-diagonal lower bound of Theorem J and Gershgorin \(\lambda_{\max}(L)\le 8\) give \(w_{\min}(\alpha)\le w_n(\alpha)\le 8^{\alpha}\), and therefore
   \[
   \frac{\rho_n}{8^{\alpha}}
     \le \sigma_n
     \le \frac{\rho_n}{w_{\min}(\alpha)}
     \le \frac{1}{w_{\min}(\alpha)}.
   \]
4. **Equivalence after AA.** For \(\alpha\in(0,1/2)\), Theorem AA supplies a level-uniform floor \(\|F_{\mathrm{frac}}(n,\alpha)\|\ge c(\alpha)>0\). Combined with the sandwich and \(w_n\ge w_{\min}>0\),
   \[
   \|F(n,\alpha)\|\to 0
     \quad\Longleftrightarrow\quad
   \sigma_n\to 0
     \quad\Longleftrightarrow\quad
   \rho_n\to 0,
   \]
   and \(\liminf_n\|F(n,\alpha)\|>0\) if and only if \(\liminf_n\rho_n>0\).

**Proof.** Reflection through the altitude from the value-\(1\) corner swaps its two incident edges and fixes the corner, so it conjugates \(L\) and therefore \(L^{\alpha}\), giving equal neighbor weights. Then \(I_{\mathrm{loc}}=w_n\bigl((u_a-u_p)+(u_a-u_q)\bigr)=w_n I_a\), and dividing by \(I_a^{\alpha}\) is the factorization. The lower bound \(w_n\ge w_{\min}\) is the integral estimate of Theorem J. For the upper bound, \(L^{\alpha}\mathbf{1}=0\) with negative off-diagonals gives \((L^{\alpha})_{aa}=\sum_{j\neq a}w_{aj}\ge w_n\), while \(e_a^{T}L^{\alpha}e_a\le\lambda_{\max}(L)^{\alpha}\le 8^{\alpha}\). The equivalence claim multiplies the sandwich by the AA floor on \(\|F_{\mathrm{frac}}\|=I_a^{\alpha}\sqrt{21}/5\) (and uses \(w_n\ge w_{\min}\) to pass \(\sigma\to 0\) to \(\rho\to 0\) and back).

**Witness, computed by `scripts/gasket_fractional_share.py`.** Through level \(6\) at \(\alpha=0.45\): neighbor weights agree to \(10^{-15}\); \(\sigma=\rho/w_n\) to \(10^{-12}\); \(w_{\min}\le w_n\le 8^{\alpha}\) holds; the sandwich and Schur audit identities hold; currents stay neutral; level \(2\) recovers the lock \(\|F\|\approx 1.165909\). Recorded shares:
\[
\begin{align*}
w_n &\approx(0.389568,\;0.344325,\;0.334438,\;0.332620,\;0.332311,\;0.332260),\\
\rho_n &\approx(0.559247,\;0.470964,\;0.453153,\;0.449983,\;0.449460,\;0.449376),\\
\sigma_n &\approx(1.435558,\;1.367791,\;1.354971,\;1.352845,\;1.352528,\;1.352483).
\end{align*}
\]
Successive \(\rho\) ratios are
\[
(0.842140,\;0.962182,\;0.993003,\;0.998838,\;0.999814).
\]
The script exits nonzero if any check fails. Passing `--through 7` recomputes one level higher; the default witness does not.

**KEPT FAILURE AB.1 (no geometric share collapse).** The law "\(\rho_{n+1}/\rho_n\le 0.95\) for every \(n\), hence \(\rho_n\to 0\) geometrically and therefore \(\|F(n,0.45)\|\to 0\) by Theorem AB plus Theorem AA" is false. The ratio from level \(4\) to level \(5\) is already \(0.99884\), and from level \(5\) to level \(6\) it is \(0.99981\). Rejecting the ratio bound rejects that particular route to zero after energy collapse was closed. It does not prove a positive limit. On the same range \(\sigma_n>1.35\) for every \(n\ge 3\) (combinatorial current exceeds fractional); that observation is not a limit theorem.

**NOT THIS.** Theorem AB does not evaluate \(\lim\rho_n\). Zero versus positive for \(\|F(n,0.45)\|\) remains OPEN; the remaining route is exactly whether the two-neighbor fractional share tends to \(0\) or stays positive. Not thrust.

## Theorem AC. The energy floor holds for every \(\alpha<\log 3/\log 5\)

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\); \(u_n\) is the \(L^{\alpha}\)-Dirichlet minimizer at level \(n\) and \(Q_{\alpha}(u_n)=u_n^{T}L^{\alpha}u_n\). Write \(\gamma=\log 3/\log 5\approx 0.682606\), \(G_k\) for the level-\(k\) graph (unit edges, \(N(k)=(3^{k+1}+3)/2\) vertices), \(E(f)=f^{T}Lf\), and \(R(x,y)=\sup_f (f(x)-f(y))^2/E(f)\) for effective resistance. **ESTABLISHED, standard network facts used here.** \(R\) is a metric; Rayleigh monotonicity (deleting edges cannot lower \(R\)); and the corner-to-corner resistance of \(G_k\) is \((2/3)(5/3)^k\) (the witness checks it for \(k\le 5\)). No Löwner–Heinz is needed.

**THEOREM AC.** For every \(\alpha\in(0,\gamma)\) and every level \(n\ge 1\),
\[
Q_{\alpha}(u_n)\ \ge\ \frac{9/4}{C_{AC}(\alpha)},
\qquad
C_{AC}(\alpha)=\frac{2\gamma}{\gamma-\alpha}\,40^{\alpha/\gamma},
\]
and therefore \(\|F_{\mathrm{frac}}(n,\alpha)\|=Q_{\alpha}\sqrt{21}/7\ge 9\sqrt{21}/(28\,C_{AC}(\alpha))\). Together with Theorem K this makes \(\gamma\) the exact threshold for a level-uniform energy floor: a floor holds for every \(\alpha<\gamma\) and fails for every \(\alpha\in(\gamma,1)\). The endpoint \(\alpha=\gamma\) is decided on the collapsing side by Theorem AD. Item 4 of Theorem AB, previously stated for \(\alpha<1/2\), therefore holds for every \(\alpha<\gamma\): \(\|F(n,\alpha)\|\to 0\) if and only if \(\rho_n\to 0\).

At \(\alpha=0.45\), \(C_{AC}\approx 66.790994\) and the floor is \(0.033687\), weaker than AA's \(0.128335\). The gain is the window \(\alpha\in[1/2,\gamma)\), where AA's path route fails: floors \(0.020184\) at \(\alpha=0.50\), \(0.005319\) at \(0.60\), \(0.000109\) at \(0.68\).

**Proof.** *Cells.* \(G_n\) is a union of \(3^{n-k}\) copies of \(G_k\) (level-\(k\) cells). For a vertex \(x\) of a cell \(K\), choose a unit triangle \(T_0\ni x\) inside \(K\) and let \(T_j\) be the level-\(j\) cell containing it, so \(T_k=K\). Each \(T_j\) contains exactly one corner \(c_{j+1}\) of \(T_{j+1}\); with \(c_0=x\), both \(c_j\) and \(c_{j+1}\) are corners of \(T_j\). The triangle inequality and Rayleigh monotonicity give \(R(x,c_k)\le\sum_{j<k}(2/3)(5/3)^j=(5/3)^k-1\), so any two vertices of \(K\) satisfy \(R(x,y)\le R_k:=(8/3)(5/3)^k\), measured in \(G_n\).

*Low modes are flat.* Let \(g\) satisfy \(E(g)\le s\|g\|^2\). For \(x,y\) in one level-\(k\) cell \(K\), \(g(x)^2\le 2g(y)^2+2(g(x)-g(y))^2\le 2g(y)^2+2R_k s\|g\|^2\). Averaging over \(y\in K\) gives \(\max g^2\le K_s\|g\|^2\) with \(K_s=2/N(k)+2R_k s\).

*Spectral mass.* Put \(v=e_a-e_b\perp\mathbf 1\) and let \(\mu\) be its spectral measure for \(L\) (total mass \(\|v\|^2=2\), all on \(\lambda>0\)). For the spectral projection \(g=P_{(0,s]}v\), \(\|g\|^2=\langle g,v\rangle=g(a)-g(b)\le 2\sqrt{K_s}\,\|g\|\), so \(\mu((0,s])=\|g\|^2\le 4K_s\). For \(s\in[5^{-(k+1)},5^{-k})\) with \(0\le k\le n\), \(2/N(k)\le(4/3)3^{-k}\) and \(2R_k s\le(16/3)3^{-k}\), so \(\mu((0,s])\le(80/3)3^{-k}\le 80\,s^{\gamma}\), because \(3^{-k}=(5^{-k})^{\gamma}\le(5s)^{\gamma}=3s^{\gamma}\).

*Gap.* For \(f\perp\mathbf 1\), \(\sum_{x,y}(f(x)-f(y))^2=2N\|f\|^2\le N^2 R_n E(f)\), so \(\lambda_1(G_n)\ge 2/(N R_n)\ge 5^{-n}/4>5^{-(n+1)}\), using \(N\le 3^{n+1}\). Hence \(\mu((0,s])=0\) for \(s<5^{-(n+1)}\), and \(\mu((0,s])\le\min(2,\,80s^{\gamma})\) for every \(s>0\).

*Layer cake.* \(R_{\alpha}(n):=v^{T}(L^{\alpha})^{+}v=\int\lambda^{-\alpha}\,d\mu=\alpha\int_0^{\infty}s^{-\alpha-1}\mu((0,s])\,ds\) (Tonelli). Split at \(s_*=40^{-1/\gamma}\), where \(80s_*^{\gamma}=2\):
\[
R_{\alpha}(n)\le\frac{80\alpha\,s_*^{\gamma-\alpha}}{\gamma-\alpha}+2s_*^{-\alpha}
=\frac{2\gamma}{\gamma-\alpha}\,s_*^{-\alpha}=C_{AC}(\alpha),
\]
finite exactly because \(\alpha<\gamma\).

*Duality.* As in Theorem AA, Cauchy–Schwarz in the \(L^{\alpha}\) inner product with \(u_n^{T}v=u_a-u_b=3/2\) gives \(Q_{\alpha}(u_n)\ge(9/4)/R_{\alpha}(n)\ge(9/4)/C_{AC}(\alpha)\). The pairing of Theorem K converts it to \(\|F_{\mathrm{frac}}\|\). \(\square\)

**Witness, computed by `scripts/gasket_fractional_capacity_floor.py`.** Corner resistance of \(G_k\) equals \((2/3)(5/3)^k\) for \(k=0,\dots,5\), and the resistance diameter stays below \((8/3)(5/3)^k\) (on these levels it equals the corner resistance). Through level \(6\): \(\lambda_1\ge 5^{-n}/4\) (level \(6\): \(1.154\cdot10^{-3}\ge 1.6\cdot10^{-5}\)); \(\max_{s<1}\mu((0,s])/s^{\gamma}\le 0.98\), far inside the proved constant \(80\); total mass \(2\). For \(\alpha\in\{0.45,0.50,0.55,0.60,0.65,0.68\}\) on levels \(1\)–\(6\), the live \(R_{\alpha}(n)\) stays below \(C_{AC}(\alpha)\) and \(Q_{\alpha}(u_n)\ge(9/4)/R_{\alpha}(n)\ge(9/4)/C_{AC}(\alpha)\), with \(\|F_{\mathrm{frac}}\|=Q\sqrt{21}/7\) to \(10^{-9}\). At level \(6\), \(R_{\alpha}(6)\approx(2.347,\;2.521,\;2.753,\;3.061,\;3.470,\;3.777)\) and \(Q_{\alpha}(u_6)\approx(0.994015,\;0.925445,\;0.847484,\;0.762196,\;0.672392,\;0.617709)\). The recorded level-\(9\) and level-\(10\) energies at \(0.45\) clear the AC floor. The script exits nonzero if any check fails.

**KEPT FAILURE AC.1 (no energy floor at every \(\alpha\)).** The law "for every \(\alpha\in(0,1)\) there is \(c(\alpha)>0\) with \(Q_{\alpha}(u_n)\ge c(\alpha)\) at every level" is false. For \(\alpha>\gamma\), Theorem K gives \(Q_{\alpha}(u_n)\le(7/2)^{\alpha}3^{1-\alpha}(3\cdot5^{-\alpha})^n\to 0\). Live witness, levels \(1\)–\(6\): \(Q_{0.75}(u_n)\approx(1.822739,\;1.304966,\;0.975341,\;0.756383,\;0.603760,\;0.492589)\) and \(Q_{0.90}(u_n)\approx(1.984971,\;1.282422,\;0.846205,\;0.568710,\;0.387631,\;0.266891)\), strictly decreasing and under the Theorem K bound. So the floor proved here cannot be pushed past \(\gamma\), and AB's equivalence cannot be carried above \(\gamma\) by an energy floor.

**NOT THIS.** Theorem AC is a statement about the fractional energy and \(\|F_{\mathrm{frac}}\|\). It does not bound the combinatorial audit norm \(\|F(n,0.45)\|\) or the share \(\rho_n\) from below, so zero versus positive at \(\alpha=0.45\) remains OPEN. The constant \(80\) is crude (the witness sees \(0.98\)); no sharp constant is claimed. Not thrust.

## Theorem AD. At \(\alpha=\log 3/\log 5\) the energy collapses like \(1/n\)

**ASSUMPTION.** Spectral calculus as in item 4, notation of Theorem AC: \(\gamma=\log 3/\log 5\), \(v=e_a-e_b\), \(\mu\) its spectral measure for \(L\), \(R_{\alpha}(n)=v^{T}(L^{\alpha})^{+}v=\int\lambda^{-\alpha}\,d\mu\). Corner data \((1,-1/2,0)\). **CITED, standard.** The Stieltjes formula \(\lambda^{-\alpha}=\frac{\sin\pi\alpha}{\pi}\int_0^{\infty}s^{-\alpha}(\lambda+s)^{-1}\,ds\) for \(\lambda>0\), \(0<\alpha<1\), and the variational formula \(v^{T}(L+s)^{-1}v=\sup_f\langle f,v\rangle^2/(E(f)+s\|f\|^2)\) for \(s>0\).

**THEOREM AD.** (i) For every \(\alpha\in(0,1]\) and every level \(n\ge 1\), the corner Schur complement of \(L^{\alpha}\) is \(\kappa_n(3I-J)\), and
\[
Q_{\alpha}(u_n)=\frac{7/3}{R_{\alpha}(n)}\quad\text{exactly}.
\]
(ii) At the endpoint \(\alpha=\gamma\), with \(c_{AD}=\frac{\sin\pi\gamma}{\pi}\cdot\frac{3}{25}\ln 5\approx 0.051634\),
\[
c_{AD}\,n\ \le\ R_{\gamma}(n)\ \le\ 80\bigl((n+1)\ln 3+1\bigr),
\qquad
\frac{7/3}{80((n+1)\ln 3+1)}\ \le\ Q_{\gamma}(u_n)\ \le\ \frac{7/3}{c_{AD}\,n}.
\]
So \(Q_{\gamma}(u_n)\to 0\) at rate exactly order \(1/n\), \(\|F_{\mathrm{frac}}(n,\gamma)\|=Q_{\gamma}\sqrt{21}/7\to 0\), and by Theorem L the combinatorial audit norm obeys \(\|F(n,\gamma)\|\le\kappa(\gamma)\sqrt{(7/3)/(c_{AD}n)}\to 0\). Together with Theorems AC and K: a level-uniform energy floor holds exactly for \(\alpha<\gamma\), and \(\|F(n,\alpha)\|\to 0\) for every \(\alpha\ge\gamma\).

**Proof.** *(i)* The level-\(n\) graph has the dihedral symmetry of the triangle, which permutes the corners in every way and commutes with \(L\) and hence with \(L^{\alpha}\). The corner Schur complement \(S\) (the energy form of the \(L^{\alpha}\)-Dirichlet minimizer as a function of its corner values) therefore commutes with all \(3\times 3\) permutations, so \(S=xI+yJ\); \(S\mathbf 1=0\) forces \(S=\kappa(3I-J)\), with \(\kappa>0\) because the kernel of \(L^{\alpha}\) is the constants. Then \(Q_{\alpha}(u_n)=\kappa\sum_{i<j}(b_i-b_j)^2=\kappa(9/4+1/4+1)=(7/2)\kappa\). For \(v\) supported on the corners and \(\perp\mathbf 1\), any \(w\) with \(L^{\alpha}w=v\) is \(\alpha\)-harmonic off the corners, so \(R_{\alpha}=w^{T}L^{\alpha}w=v^{T}S^{+}v=|v|^2/(3\kappa)=2/(3\kappa)\). Eliminating \(\kappa\) gives \(Q_{\alpha}R_{\alpha}=7/3\). (This is sharper than the duality constant \(9/4\) used in AA and AC.)

*(ii), lower bound on \(R_{\gamma}\).* Put \(R^{s}=v^{T}(L+s)^{-1}v=\int(\lambda+s)^{-1}d\mu\) (no zero-mode issue since \(v\perp\mathbf 1\)). By Stieltjes and Tonelli, \(R_{\gamma}(n)=\frac{\sin\pi\gamma}{\pi}\int_0^{\infty}s^{-\gamma}R^{s}\,ds\). For \(0\le k\le n-1\), let \(f\) be the harmonic function on the level-\(k\) corner cell at \(a\) with \(f(a)=1\) and \(f=0\) at its other two corners, extended by \(0\). The cell meets the rest of \(G_n\) only at those two corners, where \(f=0\), so \(E(f)\) is the cell energy: the cell is electrically a \(\Delta\) of three resistors \((5/3)^k\) (corner resistance \((2/3)(5/3)^k\)), so \(E(f)=2(3/5)^k\). Also \(f(b)=0\), \(\langle f,v\rangle=1\), and \(0\le f\le 1\) gives \(\|f\|^2\le N(k)\le 3^{k+1}\). For \(s\in[5^{-(k+1)},5^{-k})\), \(E(f)+s\|f\|^2\le 2(3/5)^k+3(3/5)^k\), so \(R^{s}\ge(5/3)^k/5=5^{k(1-\gamma)}/5\ge(5s)^{\gamma-1}/5\). Hence
\[
R_{\gamma}(n)\ge\frac{\sin\pi\gamma}{\pi}\int_{5^{-n}}^{1}s^{-\gamma}\,\frac{(5s)^{\gamma-1}}{5}\,ds=\frac{\sin\pi\gamma}{\pi}\,5^{\gamma-2}\,n\ln 5=c_{AD}\,n,
\]
using \(5^{\gamma}=3\).

*Upper bound on \(R_{\gamma}\).* Theorem AC's layer cake at \(\alpha=\gamma\), with \(\mu((0,s])=0\) below \(5^{-(n+1)}\), \(\le 80s^{\gamma}\) up to \(s_*=40^{-1/\gamma}\) and \(\le 2\) above: if \(s_*\ge 5^{-(n+1)}\), \(R_{\gamma}\le 80\gamma\ln(s_*5^{n+1})+2s_*^{-\gamma}=80((n+1)\ln 3-\ln 40+1)\); otherwise \(3^{n+1}<40\) and \(R_{\gamma}\le 2\cdot 5^{(n+1)\gamma}=2\cdot 3^{n+1}<80\). Both are at most \(80((n+1)\ln 3+1)\).

*Conclusion.* Insert both bounds into (i); Theorem K's pairing gives \(\|F_{\mathrm{frac}}\|\), and Theorem L's \(\|F\|\le\kappa(\gamma)\sqrt{Q_{\gamma}}\) gives the audit norm. \(\square\)

**Witness, computed by `scripts/gasket_fractional_endpoint.py`.** Corner Schur complement equals \(\kappa(3I-J)\) and \(Q_{\alpha}R_{\alpha}=7/3\) to \(10^{-12}\) on levels \(1\)–\(5\) for \(\alpha\in\{0.45,0.6,\gamma,0.75,0.9,1\}\) (it also reproduces AC's recorded pairs, e.g. \((7/3)/2.347\approx 0.994\)). The corner bump energy is \(2(3/5)^k\) for \(k\le 5\). On levels \(1\)–\(7\), \(R^{s}\) stays above \((5s)^{\gamma-1}/5\) on \([5^{-n},1)\) (worst ratio about \(5.38\), so \(c_{AD}\) is crude). Stieltjes quadrature reproduces \(R_{\gamma}(3)=2.270069\). Live at \(\alpha=\gamma\), levels \(1\)–\(7\): \(R_{\gamma}(n)\approx(1.330561,\,1.780895,\,2.270069,\,2.776998,\,3.290658,\,3.806650,\,4.323426)\), increments rising toward about \(0.517\) per level; \(Q_{\gamma}(u_n)\approx(1.753647,\,1.310202,\,1.027869,\,0.840236,\,0.709078,\,0.612962,\,0.539695)\), strictly decreasing inside the proved sandwich; \(\|F(n,\gamma)\|\approx(1.407584,\,1.004514,\,0.778089,\,0.634142,\,0.534810,\,0.462255,\,0.406991)\) under \(\kappa(\gamma)\sqrt{Q_{\gamma}}\). The script exits nonzero if any check fails.

**KEPT FAILURE AD.1 (the AC floor does not reach the endpoint).** The law "Theorem AC's level-uniform energy floor extends to \(\alpha=\gamma\)" is false: by Theorem AD, \(Q_{\gamma}(u_n)\le 45.19/n\). Theorem K could not have decided this: at \(\alpha=\gamma\) its majorant is the constant \((7/2)^{\gamma}3^{1-\gamma}\approx 3.332879\) because \(3\cdot 5^{-\gamma}=1\), and the live harmonic comparison energy \(Q_{\gamma}(u_H)\approx(1.7601,\,1.3691,\,1.1702,\,1.0768,\,1.0358,\,1.0185,\,1.0115)\) levels off near \(1\) instead of decaying, while the minimizer energy keeps falling. The AC floor itself degenerates (\(0.000109\) at \(\alpha=0.68\), \(C_{AC}\to\infty\) as \(\alpha\to\gamma\)). So the endpoint belongs to the collapsing side.

**NOT THIS.** Theorem AD is about \(\alpha=\gamma\) and the exact pairing \(Q_{\alpha}R_{\alpha}=7/3\). It does not bound \(\|F(n,0.45)\|\) or \(\rho_n\) from below; at \(\alpha=0.45\) the energy has a floor (AA, AC) and zero versus positive for the combinatorial audit norm remains OPEN. The rate \(1/n\) is for the energy; the audit-norm rate \(n^{-1/2}\) is only an upper bound. No sharp constants are claimed. Not thrust.

## Theorem AE. Degree–gap factorization of the two-neighbor share

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\). Theorems AB and AD are used as stated; this note does not re-prove them.

**THEOREM AE.** For every \(\alpha\in(0,1]\) and every level \(n\ge 1\), with the \(L^{\alpha}\)-Dirichlet minimizer \(u\):

1. **Sharp weight ceiling.** \((L^{\alpha})_{aa}=e_a^{T}L^{\alpha}e_a\le\lambda_{\max}(L)^{\alpha}\le 8^{\alpha}\). By AB's equal neighbor weights, \((L^{\alpha})_{aa}\ge 2w_n\), so
   \[
   w_n(\alpha)\le 8^{\alpha}/2.
   \]
   This sharpens AB's published ceiling \(w_n\le 8^{\alpha}\).

2. **Degree–gap factorization.** Write
   \[
   \lambda_n:=\frac{2w_n}{(L^{\alpha})_{aa}}\in(0,1]
   \]
   for the two-neighbor share of the fractional corner degree, \(W_{\mathrm{far}}:=(L^{\alpha})_{aa}-2w_n\), \(\bar g_{\mathrm{nbr}}:=I_a/2\), and (when \(W_{\mathrm{far}}>0\)) \(\bar g_{\mathrm{far}}:=I_{\mathrm{far}}/W_{\mathrm{far}}\) with \(I_{\mathrm{far}}:=I_a^{\alpha}-I_{\mathrm{loc}}\) and \(\gamma_n:=\bar g_{\mathrm{far}}/\bar g_{\mathrm{nbr}}\). Then
   \[
   \rho_n=\frac{\lambda_n}{\lambda_n+(1-\lambda_n)\gamma_n},
   \]
   and the far share is \(\phi_n:=1-\rho_n=(1-\lambda_n)\gamma_n\big/\bigl(\lambda_n+(1-\lambda_n)\gamma_n\bigr)\).

3. **Resistance identity.** By Theorem AD's pairing \(Q_{\alpha}R_{\alpha}=7/3\) and the Schur identity \(Q_{\alpha}=(7/5)I_a^{\alpha}\),
   \[
   I_a^{\alpha}=\frac{5}{3R_{\alpha}(n)},
   \]
   and therefore
   \[
   \rho_n=\tfrac35\,w_n\,I_a\,R_{\alpha}(n)=\frac{\sqrt{21}}{7}\,w_n\,\|F(n,\alpha)\|\,R_{\alpha}(n).
   \]

4. **Split of the vanishing criterion.** For \(\alpha\in(0,\gamma)\) with \(\gamma=\log3/\log5\), Theorems AB and AC give \(\|F\|\to0\) if and only if \(\rho_n\to0\). By the factorization, \(\rho_n\to0\) if and only if \(\bigl((1-\lambda_n)/\lambda_n\bigr)\gamma_n\to+\infty\).

**Proof.** (1) The Rayleigh quotient \(e_a^{T}L^{\alpha}e_a\le\lambda_{\max}^{\alpha}\le8^{\alpha}\) is the same bound AB used for a single neighbor weight; equal neighbor weights upgrade it to \(2w_n\le(L^{\alpha})_{aa}\). (2) \(I_{\mathrm{loc}}=2w_n\bar g_{\mathrm{nbr}}=\lambda_n\cdot\mathrm{diag}\cdot\bar g_{\mathrm{nbr}}\) and \(I_a^{\alpha}=\mathrm{diag}\bigl(\lambda_n\bar g_{\mathrm{nbr}}+(1-\lambda_n)\bar g_{\mathrm{far}}\bigr)\), so the ratio is the displayed Möbius form in \((\lambda_n,\gamma_n)\). (3) Eliminate \(\kappa\) from AD's \(Q=(7/2)\kappa\) and \(R=2/(3\kappa)\), then use \(Q=(7/5)I_a^{\alpha}\) and \(\|F\|=I_a\sqrt{21}/5\). (4) Algebra from (2) together with AB+AC.

**Witness, computed by `scripts/gasket_fractional_degree_share.py`.** Through level \(6\) at \(\alpha=0.45\): neighbor weights agree to \(10^{-15}\); \(w_n\le8^{\alpha}/2\); the factorization \(\rho=\lambda/(\lambda+(1-\lambda)\gamma)\) and the AD identities \(QR=7/3\), \(I_a^{\alpha}=5/(3R)\), \(\rho=(3/5)w_n I_a R\) hold to \(10^{-10}\); currents stay neutral; level \(2\) recovers the lock \(\|F\|\approx1.165909\). Recorded degree and gap ratios:
\[
\begin{aligned}
\lambda_n&\approx(0.646893,\;0.553135,\;0.533143,\;0.529482,\;0.528861,\;0.528759),\\
\gamma_n&\approx(1.443835,\;1.390438,\;1.378096,\;1.375485,\;1.374963,\;1.374864).
\end{aligned}
\]
Successive \(\lambda\) ratios are
\[
(0.855064,\;0.963857,\;0.993134,\;0.998827,\;0.999807).
\]
The script exits nonzero if any check fails.

**KEPT FAILURE AE.1 (no geometric neighbor-degree wipeout).** The law "\(\lambda_{n+1}/\lambda_n\le0.95\) for every \(n\), hence \(\lambda_n\to0\) geometrically, and therefore (with any uniform positive lower bound on \(\gamma_n\)) \(\rho_n\to0\) and \(\|F(n,0.45)\|\to0\) by Theorems AB and AA" is false. The ratio from level \(3\) to level \(4\) is already \(0.99314\), and from level \(4\) to level \(5\) it is \(0.99883\). Rejecting geometric wipeout of the neighbor-degree share rejects that particular route to zero after AB.1 closed geometric wipeout of \(\rho\) itself. It does not prove a positive limit, and it does not bound \(\gamma_n\).

**NOT THIS.** Theorem AE factors \(\rho_n\) into a degree share and a gap ratio and ties it to AD's resistance. Finite-level \(\lambda_n\approx0.529\) and \(\gamma_n\approx1.375\) are not limit theorems. Zero versus positive for \(\|F(n,0.45)\|\) remains OPEN; the remaining route is whether \(\bigl((1-\lambda_n)/\lambda_n\bigr)\gamma_n\to\infty\). Not thrust.

## Theorem AF. The neighbor-degree share has a uniform floor

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\). Notation of Theorem AE: \(w_n=-(L^{\alpha})_{ap}=-(L^{\alpha})_{aq}\), \(\lambda_n=2w_n/(L^{\alpha})_{aa}\), \(\gamma_n=\bar g_{\mathrm{far}}/\bar g_{\mathrm{nbr}}\), and \(w_{\min}(\alpha)=\alpha\,4^{\alpha-1}\) from Theorem J (at \(\alpha=1\), \(w_n=1=w_{\min}(1)\) exactly). Theorems J, AB, AC, and AE are used as stated; this note does not re-prove them.

**THEOREM AF.** For every \(\alpha\in(0,1]\) and every level \(n\ge 1\):

1. **Uniform floor.**
   \[
   \lambda_n(\alpha)\ \ge\ \lambda_{\min}(\alpha):=\frac{2w_{\min}(\alpha)}{8^{\alpha}}=\frac{2\alpha\,4^{\alpha-1}}{8^{\alpha}}=\alpha\,2^{-1-\alpha}>0.
   \]
2. **Jensen sharpening.** \((L^{\alpha})_{aa}\le (L_{aa})^{\alpha}=2^{\alpha}\). Therefore
   \[
   \lambda_n(\alpha)\ \ge\ \lambda_{J}(\alpha):=\frac{2w_{\min}(\alpha)}{2^{\alpha}}=\alpha\,2^{\alpha-1}\ \ge\ \lambda_{\min}(\alpha),
   \qquad
   w_n(\alpha)\le 2^{\alpha-1},
   \]
   which sharpens AE's ceiling \(w_n\le 8^{\alpha}/2\) (equality \(w_n=1\) at \(\alpha=1\)).
3. **One-sided vanishing split.** For \(\alpha\in(0,\gamma)\), \(\gamma=\log3/\log5\), put \(C_J(\alpha)=1/\lambda_J(\alpha)-1\in[0,\infty)\). Then
   \[
   \rho_n\ \ge\ \frac{1}{1+C_J(\alpha)\,\gamma_n},
   \qquad\text{hence}\qquad
   \|F(n,\alpha)\|\to0\ \Longrightarrow\ \gamma_n\to+\infty .
   \]
   The same holds with the cruder \(C=1/\lambda_{\min}-1\).

**Proof.** (1) AE item 1: \((L^{\alpha})_{aa}=e_a^{T}L^{\alpha}e_a\le\lambda_{\max}(L)^{\alpha}\le 8^{\alpha}\), and \(w_n\ge w_{\min}\) (Theorem J), so \(\lambda_n=2w_n/(L^{\alpha})_{aa}\ge 2w_{\min}/8^{\alpha}\). The closed form is \(2\alpha\,2^{2\alpha-2}2^{-3\alpha}=\alpha\,2^{-1-\alpha}\). (2) Let \(\mu_a=\sum_k v_k(a)^2\,\delta_{\lambda_k}\) be the spectral measure of the unit vector \(e_a\) for \(L\); it is a probability measure on \([0,8]\) with mean \(\int t\,d\mu_a=L_{aa}=\deg a=2\) (degree lemma of Theorem J). For \(\alpha\in(0,1]\), \(t\mapsto t^{\alpha}\) is concave on \([0,\infty)\), so Jensen gives \((L^{\alpha})_{aa}=\int t^{\alpha}\,d\mu_a\le 2^{\alpha}\). Then \(\lambda_n\ge 2w_{\min}/2^{\alpha}=\alpha\,2^{\alpha-1}\), and \(2w_n\le(L^{\alpha})_{aa}\) (AE) gives \(w_n\le2^{\alpha-1}\). Since \(2^{\alpha}\le 8^{\alpha}\), \(\lambda_J\ge\lambda_{\min}\). (3) For \(\alpha<1\), all off-diagonal entries of \(L^{\alpha}\) are strictly negative (Theorem J) and \(N>3\), so \(W_{\mathrm{far}}>0\) and \(\gamma_n\) is defined; the maximum principle (\(u\le1\)) gives \(I_{\mathrm{far}}\ge0\), so \(\gamma_n\ge0\). From \(\lambda_J\le\lambda_n\le1\), \(0\le(1-\lambda_n)/\lambda_n\le C_J\). Dividing AE's factorization by \(\lambda_n\),
\[
\rho_n=\frac{1}{1+\frac{1-\lambda_n}{\lambda_n}\gamma_n}\ \ge\ \frac{1}{1+C_J\gamma_n}.
\]
By AB+AC (AE item 4), \(\|F\|\to0\) forces \(\rho_n\to0\), hence \(C_J\gamma_n\to+\infty\), hence \(\gamma_n\to+\infty\) (\(C_J\) is a fixed finite number). \(\square\)

**NOT CLAIMED (converse).** \(\gamma_n\to\infty\Rightarrow\|F\|\to0\) is not proved. It would follow if, in addition, \(\limsup_n\lambda_n\le\lambda_{\max}<1\) uniformly, because then \((1-\lambda_n)/\lambda_n\ge(1-\lambda_{\max})/\lambda_{\max}>0\). That \(\limsup\) is not proved here; if \(\lambda_n\to1\), \(\gamma_n\to\infty\) alone does not force \(\rho_n\to0\).

**Witness, computed by `scripts/gasket_fractional_lambda_floor.py`.** At \(\alpha=0.45\):
\[
\lambda_{\min}(0.45)=0.45\cdot2^{-1.45}\approx0.164710,\qquad
\lambda_J(0.45)=0.45\cdot2^{-0.55}\approx0.307359,\qquad
C_J\approx2.253524,\quad C\approx5.071290 .
\]
Through level \(6\): \((L^{\alpha})_{aa}\approx(1.204427,\,1.244994,\,1.254590,\,1.256396,\,1.256704,\,1.256755)\le2^{0.45}\approx1.366040\); \(w_{\min}\le w_n\le2^{-0.55}\); \(\lambda_n\ge\lambda_J\ge\lambda_{\min}\); \(\rho_n\ge1/(1+C_J\gamma_n)\); AE's factorization still holds to \(10^{-12}\); the recorded AE table of \(\lambda_n,\gamma_n\) is reproduced. The floors \(\lambda_n\ge\lambda_J\ge\lambda_{\min}\) and \((L^{\alpha})_{aa}\le2^{\alpha}\) also hold live through level \(5\) at \(\alpha\in\{0.10,0.30,0.60,0.68,0.80,0.95\}\); at \(\alpha=1\), \(\lambda_n=1\) exactly; at level \(3\), \((L^{0.45})_{ii}\le(L_{ii})^{0.45}\) at every vertex. Successive \(\gamma\) ratios are
\[
(0.963017,\;0.991123,\;0.998105,\;0.999620,\;0.999928).
\]
The script exits nonzero if any check fails.

**KEPT FAILURE AF.1 (no geometric gap-ratio blowup).** The law "\(\gamma_{n+1}/\gamma_n\ge1.05\) for every \(n\), hence \(\gamma_n\to\infty\) geometrically, and therefore \(\|F(n,0.45)\|\to0\) by Theorems AF and AE" is false. Every live ratio through level \(6\) is below \(1\), starting at \(0.96302\) from level \(1\) to level \(2\); \(\gamma_n\) decreases. Since AF makes \(\gamma_n\to\infty\) necessary for \(\|F\|\to0\), this rejects the geometric version of the only remaining route to zero. It does not prove a positive limit: \(\gamma_n\) could still turn and drift to \(\infty\) slowly beyond the computed range, and nothing here bounds \(\gamma_n\) above.

**OBSERVATION, not a limit theorem.** On levels \(1\)–\(6\) at \(\alpha=0.45\), \(\gamma_n\in(1,2)\) is decreasing toward \(\approx1.375\), with ratios rising toward \(1\) from below, and \(\lambda_n\) is decreasing toward \(\approx0.529\), well above \(\lambda_J\approx0.307\).

**NOT THIS.** Theorem AF bounds \(\lambda_n\) below and so reduces \(\|F\|\to0\) to a necessary condition on \(\gamma_n\) alone. The missing ceiling \(\limsup\lambda_n<1\) is supplied by Theorem AG below. Theorem AF itself does not bound \(\gamma_n\) above and does not bound the combinatorial \(\|F(n,0.45)\|\) below by a positive constant. Zero versus positive for \(\|F(n,0.45)\|\) remains OPEN; Conjecture J.1 is untouched. Not thrust.

## Theorem AG. The neighbor-degree share has a uniform ceiling below 1

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\). Notation of Theorems AE–AF: \(w_n=-(L^{\alpha})_{ap}=-(L^{\alpha})_{aq}\), \(\lambda_n=2w_n/(L^{\alpha})_{aa}\), \(W_{\mathrm{far}}=(L^{\alpha})_{aa}-2w_n\), \(\gamma_n=\bar g_{\mathrm{far}}/\bar g_{\mathrm{nbr}}\). The integral representation and \(M=4I-L\) heat-kernel comparison are those of Theorem J. Theorems J, AB, AC, AE, and AF are used as stated; this note does not re-prove them.

**THEOREM AG.** For every \(\alpha\in(0,1)\) and every level \(n\ge 1\):

1. **Far-weight floor.** The three vertices at combinatorial distance \(2\) from the value-\(1\) corner contribute
   \[
   W_{\mathrm{far}}(n,\alpha)\ \ge\ W_{\mathrm{far}}^{\min}(\alpha):=\frac{2\alpha(1-\alpha)}{4^{2-\alpha}}>0.
   \]
2. **Uniform ceiling.** With AF's Jensen bound \((L^{\alpha})_{aa}\le 2^{\alpha}\),
   \[
   \lambda_n(\alpha)\ \le\ \lambda_{\mathrm{ceil}}(\alpha):=1-\frac{\alpha(1-\alpha)}{2^{3-\alpha}}<1.
   \]
   At \(\alpha=0.45\), \(\lambda_{\mathrm{ceil}}\approx 0.957738\).
3. **Two-sided vanishing split.** For \(\alpha\in(0,\gamma)\) with \(\gamma=\log3/\log5\), Theorems AF and AB+AC give
   \[
   \|F(n,\alpha)\|\to 0\ \Longleftrightarrow\ \gamma_n\to+\infty.
   \]
   Equivalently, with \(c_*(\alpha)=(1-\lambda_{\mathrm{ceil}})/\lambda_{\mathrm{ceil}}>0\) and AF's \(C_J=1/\lambda_J-1\),
   \[
   \frac{1}{1+C_J\gamma_n}\ \le\ \rho_n\ \le\ \frac{1}{1+c_*\gamma_n}.
   \]

**Proof.** (1) Let \(a\) be the value-\(1\) corner. On every level \(n\ge 1\) the combinatorial neighbourhood of \(a\) is the same six-vertex patch: neighbors \(p,q\) (degree \(4\)) and three distance-\(2\) vertices \(f_1,f_2,f_3\) that are not neighbors of \(a\), with two length-\(2\) walks from \(a\) to the middle far vertex and one each to the other two. The matrix \(M=4I-L\) is entrywise nonnegative (degree lemma of Theorem J), and the path counts \((M^2)_{af}\) equal \(1,2,1\) independently of \(n\). Since \(e^{-sL}=e^{-4s}e^{sM}\) and \(e^{sM}\ge (s^2/2)M^2\) entrywise on the off-diagonal,
\[
(e^{-sL})_{af}\ \ge\ \frac{s^2}{2}\,(M^2)_{af}\,e^{-4s}.
\]
Theorem J's integral representation therefore yields
\[
-(L^{\alpha})_{af}\ \ge\ (M^2)_{af}\cdot\frac{\alpha(1-\alpha)}{2\cdot 4^{2-\alpha}}
\]
for each of the three far vertices. Summing gives \(W_{\mathrm{far}}\ge 2\alpha(1-\alpha)/4^{2-\alpha}\). (2) \(\lambda_n=1-W_{\mathrm{far}}/(L^{\alpha})_{aa}\le 1-W_{\mathrm{far}}^{\min}/2^{\alpha}\), and \(W_{\mathrm{far}}^{\min}/2^{\alpha}=\alpha(1-\alpha)/2^{3-\alpha}\). (3) AF supplies \(\lambda_J\le\lambda_n\le\lambda_{\mathrm{ceil}}<1\), so
\[
c_*\ \le\ \frac{1-\lambda_n}{\lambda_n}\ \le\ C_J.
\]
Dividing AE's factorization by \(\lambda_n\) sandwiches \(\rho_n\) between \(1/(1+C_J\gamma_n)\) and \(1/(1+c_*\gamma_n)\). Hence \(\rho_n\to 0\) if and only if \(\gamma_n\to+\infty\). Combined with AB+AC (\(\|F\|\to 0\iff\rho_n\to 0\)) and AF's one-sided half, that is the two-sided claim. \(\square\)

**Witness, computed by `scripts/gasket_fractional_lambda_ceiling.py`.** At \(\alpha=0.45\):
\[
W_{\mathrm{far}}^{\min}\approx 0.057731,\qquad
\lambda_{\mathrm{ceil}}\approx 0.957738,\qquad
c_*\approx 0.044127,\qquad
C_J\approx 2.253524.
\]
Through level \(6\): \(W_{\mathrm{far}}\ge W_{\mathrm{far}}^{\min}\); \(\lambda_J\le\lambda_n\le\lambda_{\mathrm{ceil}}\); \(c_*\le(1-\lambda_n)/\lambda_n\le C_J\); the sandwich \(1/(1+C_J\gamma_n)\le\rho_n\le 1/(1+c_*\gamma_n)\) holds; AE's factorization and the recorded AE/AF table of \(\lambda_n,\gamma_n\) are reproduced; the path counts \((M^2)_{af}\in\{1,2,1\}\) and the local-patch isolation of the corner are checked through level \(4\); the off-diagonal heat-kernel lower bounds hold live at level \(3\); the ceiling \(\lambda_n\le\lambda_{\mathrm{ceil}}\) and floor \(W_{\mathrm{far}}\ge W_{\mathrm{far}}^{\min}\) also hold through level \(5\) at \(\alpha\in\{0.10,0.30,0.60,0.68,0.80,0.95\}\). Successive \(\lambda\) ratios are
\[
(0.855064,\;0.963857,\;0.993134,\;0.998827,\;0.999807).
\]
The script exits nonzero if any check fails.

**KEPT FAILURE AG.1 (no geometric climb of \(\lambda\) toward \(1\)).** The law "\(\lambda_{n+1}/\lambda_n\ge 1.01\) for every \(n\), hence \(\lambda_n\to 1\) and AF's converse obstruction (\(\lambda\to 1\) blocking \(\gamma\to\infty\Rightarrow\rho\to 0\)) materializes" is false. Every live ratio through level \(6\) is below \(1\), starting at \(0.85506\) from level \(1\) to level \(2\); \(\lambda_n\) decreases. Theorem AG already closes that obstruction by a uniform ceiling \(\lambda_n\le\lambda_{\mathrm{ceil}}<1\); the failure rejects the geometric version of the same escape. It does not prove a positive limit: \(\gamma_n\) could still drift to \(\infty\) slowly, and nothing here bounds \(\gamma_n\) above.

**OBSERVATION, not a limit theorem.** On levels \(1\)–\(6\) at \(\alpha=0.45\), \(\lambda_n\) decreases toward \(\approx 0.529\) (well below \(\lambda_{\mathrm{ceil}}\approx 0.958\) and above \(\lambda_J\approx 0.307\)), and \(\gamma_n\) decreases toward \(\approx 1.375\).

**NOT THIS.** Theorem AG bounds \(\lambda_n\) above by a constant strictly less than \(1\) and so upgrades AF's one-sided vanishing split to a two-sided equivalence with \(\gamma_n\to\infty\). Theorem AG itself does not bound \(\gamma_n\) above and does not bound the combinatorial \(\|F(n,0.45)\|\) below by a positive constant; both are supplied by Theorem AH below. Conjecture J.1 is untouched. Not thrust.

## Theorem AH. The far-to-neighbor gap ratio has a uniform ceiling; the audit norm does not vanish below \(\log 3/\log 5\)

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\). Notation of Theorems AE–AG: \(a\) is the value-\(1\) corner, \(p,q\) its two neighbors (degree \(4\)), \(\mathrm{Far}=V\setminus\{a,p,q\}\) (every vertex at combinatorial distance \(\ge2\) from \(a\), the corners \(b,c\) included), \(w_{xy}=-(L^{\alpha})_{xy}>0\) for \(x\ne y\), \(H_t=e^{-tL}=e^{-4t}e^{tM}\) with \(M=4I-L\ge0\) entrywise, \(I_a=(Lu)_a\), \(I_a^{\alpha}=(L^{\alpha}u)_a\), \(I_{\mathrm{far}}=\sum_{j\in\mathrm{Far}}w_{aj}(u_a-u_j)\), \(W_{\mathrm{far}}=\sum_{j\in\mathrm{Far}}w_{aj}\), \(\gamma_n=\bar g_{\mathrm{far}}/\bar g_{\mathrm{nbr}}=2I_{\mathrm{far}}/(W_{\mathrm{far}}I_a)\), \(\sigma_n=I_a/I_a^{\alpha}=\|F\|/\|F_{\mathrm{frac}}\|\). The integral representation, the degree lemma, and the maximum principle are those of Theorem J. Theorems AA, AB, AC, AE, AF, and AG are used as stated; this note does not re-prove them. Fix \(s=1/10\) and put
\[
K(\alpha)=\frac{e^{8s}\,\alpha}{(2-\alpha)(1-2^{-\alpha})},\qquad
\kappa(\alpha)=\frac{s\,e^{-4s}\,2^{-1-\alpha}}{1+K(\alpha)}.
\]

**THEOREM AH.** For every \(\alpha\in(0,1)\) and every level \(n\ge 1\):

1. **Kernel comparison.** For every \(j\in\mathrm{Far}\) and \(x\in\{p,q\}\),
   \[
   w_{xj}\ \ge\ \kappa(\alpha)\,w_{aj}.
   \]
2. **One-step Harnack at the neighbors.**
   \[
   I_{\mathrm{far}}\ \le\ \frac{4^{\alpha}}{2\kappa(\alpha)}\,I_a .
   \]
3. **Uniform gap-ratio ceiling.** With AG's \(W_{\mathrm{far}}\ge W_{\mathrm{far}}^{\min}(\alpha)=2\alpha(1-\alpha)/4^{2-\alpha}\),
   \[
   \gamma_n\ \le\ \frac{4^{\alpha}}{\kappa(\alpha)\,W_{\mathrm{far}}}\ \le\ \Gamma_{AH}(\alpha):=\frac{4^{\alpha}}{\kappa(\alpha)\,W_{\mathrm{far}}^{\min}(\alpha)}<\infty .
   \]
4. **Share floors.** With AF's \(w_{\min}\le w_n\le 2^{\alpha-1}\),
   \[
   \sigma_n\ \ge\ \sigma_{\min}(\alpha):=\Bigl(2^{\alpha-1}+\frac{4^{\alpha}}{2\kappa(\alpha)}\Bigr)^{-1},\qquad
   \rho_n\ \ge\ \rho_{\min}(\alpha):=\frac{w_{\min}}{w_{\min}+4^{\alpha}/(2\kappa(\alpha))}.
   \]
5. **Positive audit floor below \(\log3/\log5\).** For \(\alpha\in(0,\gamma)\), \(\gamma=\log3/\log5\), Theorem AC gives
   \[
   \|F(n,\alpha)\|\ \ge\ \sigma_{\min}(\alpha)\,\|F_{\mathrm{frac}}(n,\alpha)\|\ \ge\ \sigma_{\min}(\alpha)\,\frac{9\sqrt{21}}{28\,C_{AC}(\alpha)}\ >\ 0
   \]
   at every level. At \(\alpha=0.45\), with Theorem AA's \(\|F_{\mathrm{frac}}\|\ge 0.084015\),
   \[
   \|F(n,0.45)\|\ \ge\ \sigma_{\min}(0.45)\cdot 0.084015\ \approx\ 6.442\cdot10^{-4}\qquad\text{for every }n\ge1 .
   \]
   Hence \(\|F(n,0.45)\|\not\to 0\): together with Theorem J, \(6.442\cdot10^{-4}\le\|F(n,0.45)\|\le\tfrac35\sqrt{21}\) at every level. By AG (3) this is the same statement as \(\gamma_n\not\to\infty\), now proved directly by (3).

**Proof.** (1) Let \(c_{\alpha}=\alpha/\Gamma(1-\alpha)\), so \(w_{xy}=c_{\alpha}\int_0^{\infty}H_t(x,y)\,t^{-1-\alpha}\,dt\) for \(x\ne y\) (Theorem J). Fix \(j\in\mathrm{Far}\). *Small times.* Since \(j\ne a\) and \(j\not\sim a\), \((M^0)_{aj}=(M^1)_{aj}=0\), so \(H_{\tau}(a,j)=e^{-4\tau}S(\tau)\) with \(S(\tau)=\sum_{k\ge2}\tau^k(M^k)_{aj}/k!\) nondecreasing. For \(\tau\le s\le\tau'\le 2s\), \(\tau^k\le(\tau/s)^2s^k\) for \(k\ge2\) gives \(H_{\tau}(a,j)\le(\tau/s)^2S(s)\le(\tau/s)^2e^{8s}H_{\tau'}(a,j)\). With \(m=\min_{[s,2s]}H_{\cdot}(a,j)\),
\[
\int_0^{s}H_{\tau}(a,j)\tau^{-1-\alpha}d\tau\le\frac{e^{8s}m\,s^{-\alpha}}{2-\alpha}
= K(\alpha)\,\frac{m\,s^{-\alpha}(1-2^{-\alpha})}{\alpha}\le K(\alpha)\int_s^{2s}H_{\tau}(a,j)\tau^{-1-\alpha}d\tau,
\]
so \(w_{aj}\le c_{\alpha}(1+K)\int_s^{\infty}H_{\tau}(a,j)\tau^{-1-\alpha}d\tau\). *Semigroup.* \(H_{\tau+s}(p,j)\ge H_s(p,a)H_{\tau}(a,j)\) entrywise, and \(H_s(p,a)\ge e^{-4s}sM_{pa}=s e^{-4s}\). Therefore
\[
w_{pj}\ge c_{\alpha}\int_0^{\infty}H_{\tau+s}(p,j)(\tau+s)^{-1-\alpha}d\tau
\ge c_{\alpha}\,s e^{-4s}\int_s^{\infty}H_{\tau}(a,j)(2\tau)^{-1-\alpha}d\tau
\ge\frac{s e^{-4s}2^{-1-\alpha}}{1+K}\,w_{aj},
\]
using \(\tau+s\le2\tau\) for \(\tau\ge s\). The same holds with \(q\) in place of \(p\).

(2) Put \(v=u_a-u\). By the maximum principle \(v\ge0\), \(v_a=0\), and \(L^{\alpha}\mathbf 1=0\) makes \(v\) \(L^{\alpha}\)-harmonic at the interior vertex \(p\): \((L^{\alpha})_{pp}\,v_p=\sum_{j\ne p}w_{pj}v_j\). Dropping \(j=a\) (\(v_a=0\)) and \(j=q\) (a nonnegative term) and applying (1),
\[
(L^{\alpha})_{pp}\,v_p\ \ge\ \sum_{j\in\mathrm{Far}}w_{pj}v_j\ \ge\ \kappa\sum_{j\in\mathrm{Far}}w_{aj}v_j=\kappa\,I_{\mathrm{far}} .
\]
Jensen as in AF, now at a degree-\(4\) vertex, gives \((L^{\alpha})_{pp}\le(L_{pp})^{\alpha}=4^{\alpha}\). The same holds at \(q\). Adding, \(I_a=v_p+v_q\ge 2\kappa I_{\mathrm{far}}/4^{\alpha}\).

(3) \(\gamma_n=2I_{\mathrm{far}}/(W_{\mathrm{far}}I_a)\le 4^{\alpha}/(\kappa W_{\mathrm{far}})\), then AG (1).

(4) \(I_a^{\alpha}=\sum_{j\ne a}w_{aj}v_j=w_nI_a+I_{\mathrm{far}}\le\bigl(w_n+4^{\alpha}/(2\kappa)\bigr)I_a\) and \(w_n\le2^{\alpha-1}\) (AF). For \(\rho_n=w_nI_a/I_a^{\alpha}\ge w_n/(w_n+4^{\alpha}/(2\kappa))\), the right side increases in \(w_n\ge w_{\min}\).

(5) \(\|F\|=(\sqrt{21}/5)I_a\) and \(\|F_{\mathrm{frac}}\|=(\sqrt{21}/5)I_a^{\alpha}\) (Schur line, Theorems J and AB), so \(\|F\|=\sigma_n\|F_{\mathrm{frac}}\|\ge\sigma_{\min}\|F_{\mathrm{frac}}\|\). Insert AC's floor \(\|F_{\mathrm{frac}}\|\ge 9\sqrt{21}/(28C_{AC})\) for \(\alpha<\gamma\), or AA's \(0.084015\) at \(\alpha=0.45\). \(\square\)

Items (1)–(4) hold for every \(\alpha\in(0,1)\), including \(\alpha>\gamma\); there \(\|F_{\mathrm{frac}}\|\to0\) (Theorem K) and (5) gives nothing, consistent with \(\|F\|\to0\). At \(\alpha=\gamma\), \(\rho_n\ge\rho_{\min}\) is consistent with AD because \(R_{\gamma}(n)\to\infty\).

**Witness, computed by `scripts/gasket_fractional_gamma_ceiling.py`.** At \(\alpha=0.45\):
\[
K\approx 2.411299,\quad
\kappa\approx 7.1923\cdot10^{-3},\quad
\frac{4^{\alpha}}{2\kappa}\approx 129.726,\quad
\Gamma_{AH}\approx 4494.13,\quad
\sigma_{\min}\approx 7.6682\cdot10^{-3},\quad
\rho_{\min}\approx 1.6157\cdot10^{-3},
\]
and the audit floor is \(\approx 6.4424\cdot10^{-4}\) (AA route) or \(\approx 1.6911\cdot10^{-4}\) (AC route). The heat-kernel steps are checked by quadrature at levels \(1\)–\(3\) (\(\alpha=0.45\)) and level \(2\) (\(\alpha=0.10,\,0.90\)), with \(H_t\) from the nonnegative \(M\)-series for \(t\le4\): the subordination integral reproduces \(w_{aj}\); \(\int_0^s/\int_s^{2s}\le0.658<K\) at \(0.45\); the pointwise small-time bound; \(H_s(p,a)\approx0.078557\ge se^{-4s}\approx0.067032\); and the chain \(w_{pj}\ge c_{\alpha}H_s(p,a)\int_s^{\infty}H_{\tau}(a,j)(\tau+s)^{-1-\alpha}\ge\kappa w_{aj}\) for every far \(j\). Through level \(6\) at \(\alpha=0.45\): \(v\ge0\), \(v\) harmonic at \(p,q\), \((L^{\alpha})_{pp}\le4^{\alpha}\), \(\min_jw_{pj}/w_{aj}\ge\kappa\), \(I_{\mathrm{far}}\le(4^{\alpha}/(2\kappa))I_a\), \(\gamma_n\le4^{\alpha}/(\kappa W_{\mathrm{far}})\le\Gamma_{AH}\), \(\sigma_n\ge\sigma_{\min}\), \(\rho_n\ge\rho_{\min}\), and \(\|F\|\ge\sigma_{\min}\|F_{\mathrm{frac}}\|\ge6.44\cdot10^{-4}\); the recorded AE table of \(\gamma_n\) and the level-\(2\) lock \(1.165909\) are reproduced; the same items hold through level \(5\) at \(\alpha\in\{0.10,0.30,0.50,0.60,0.68,0.80,0.95\}\), with the AC-route audit floor below \(\gamma\). Live values:
\[
I_{\mathrm{far}}/I_a\approx(0.307026,\;0.386781,\;0.403586,\;0.406564,\;0.407045,\;0.407121),\qquad
\min_{j}w_{pj}/w_{aj}\approx(1.134125,\;1.052430,\;1.014906,\;1.003499,\;1.000745,\;1.000152).
\]
The script exits nonzero if any check fails.

**KEPT FAILURE AH.1 (Theorem J's gap cap alone gives no floor).** The shortcut "every far gap satisfies \(u_a-u_j\le 3/2\) (maximum principle, Theorem J), so \(I_{\mathrm{far}}\le\tfrac32W_{\mathrm{far}}\) and \(I_a=(I_a^{\alpha}-I_{\mathrm{far}})/w_n\ge(I_a^{\alpha}-\tfrac32W_{\mathrm{far}})/w_n>0\) at every level" is false. Live at \(\alpha=0.45\), \(I_a^{\alpha}<\tfrac32W_{\mathrm{far}}\) at levels \(3\) through \(6\), and the shortcut bound is
\[
(1.171102,\;0.277444,\;-0.141539,\;-0.347367,\;-0.464237,\;-0.536750),
\]
negative from level \(3\). A sup-norm bound on the far gaps cannot see that \(L^{\alpha}\)-harmonicity at \(p,q\) ties the far gaps to the neighbor gaps; Theorem AH (1)–(2) supplies that tie.

**OBSERVATION, not a theorem.** Through level \(6\) the live kernel ratio \(\min_jw_{pj}/w_{aj}\) stays \(\ge1\) and decreases toward \(1\), so the proved \(\kappa\approx7.2\cdot10^{-3}\) is crude by about two orders of magnitude, and the proved audit floor \(6.44\cdot10^{-4}\) is far below the live \(\|F(6,0.45)\|\approx0.880109\). No sharp constant is claimed. Theorem AI below removes \(\kappa\) from the Harnack step by comparing the *sum* of the two neighbor kernels with \(2w_{aj}\).

**NOT THIS.** Theorem AH proves that the combinatorial audit norm \(\|F(n,\alpha)\|\) has a level-uniform positive floor for every \(\alpha<\log3/\log5\), in particular at \(\alpha=0.45\), so it does not tend to \(0\). It does not prove that \(\|F(n,0.45)\|\) converges, does not prove that it decreases, and does not prove Conjecture J.1's value \(L\ge0.70\). A positive floor on a graph audit norm is not thrust, not a force, and not a pass of any experimental criterion. Not thrust.

## Theorem AI. A two-neighbor heat-kernel sum sharpens the neighbor Harnack

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\), used only for the numerical floors. Notation of Theorem AH: \(a\) is the value-\(1\) corner (\(\deg a=2\)), \(p,q\) its two neighbors (\(\deg p=\deg q=4\), \(p\sim q\)), \(\mathrm{Far}=V\setminus\{a,p,q\}\), \(w_{xy}=-(L^{\alpha})_{xy}>0\) for \(x\ne y\), \(H_t=e^{-tL}\), \(c_{\alpha}=\alpha/\Gamma(1-\alpha)\), so \(w_{xy}=c_{\alpha}\int_0^{\infty}H_t(x,y)\,t^{-1-\alpha}\,dt\) for \(x\ne y\) (Theorem J). Theorems AA, AC, AF, AG, and AH are used as stated; this note does not re-prove them.

**First-step identity.** Row \(a\) of \(L\) is \(2e_a-e_p-e_q\), so the backward equation \(\partial_tH_t=-LH_t\) at \(a\) reads, for every \(j\ne a\) (where \(H_0(a,j)=0\)),
\[
\partial_tH_t(a,j)+2H_t(a,j)=H_t(p,j)+H_t(q,j),\qquad
H_t(a,j)=\int_0^te^{-2(t-s)}\bigl(H_s(p,j)+H_s(q,j)\bigr)\,ds .
\]

**THEOREM AI.** For every \(\alpha\in(0,1)\) and every level \(n\ge1\):

1. **Two-neighbor kernel sum.** For every \(j\in\mathrm{Far}\),
   \[
   w_{pj}+w_{qj}\ \ge\ 2\,w_{aj}.
   \]
2. **Sharpened one-step Harnack at the neighbors.**
   \[
   I_{\mathrm{far}}\ \le\ \frac{4^{\alpha}}{2}\,I_a .
   \]
   This is Theorem AH (2) with \(\kappa\) replaced by \(1\).
3. **Sharpened gap-ratio ceiling.** With AG's \(W_{\mathrm{far}}\ge W_{\mathrm{far}}^{\min}=2\alpha(1-\alpha)/4^{2-\alpha}\),
   \[
   \gamma_n\ \le\ \frac{4^{\alpha}}{W_{\mathrm{far}}}\ \le\ \Gamma_{AI}(\alpha):=\frac{4^{\alpha}}{W_{\mathrm{far}}^{\min}}=\frac{8}{\alpha(1-\alpha)} .
   \]
4. **Sharpened share floors.** With AF's \(w_{\min}=\alpha4^{\alpha-1}\le w_n\le2^{\alpha-1}\),
   \[
   \sigma_n\ \ge\ \sigma_{AI}(\alpha):=\Bigl(2^{\alpha-1}+\frac{4^{\alpha}}{2}\Bigr)^{-1}=\frac{2}{2^{\alpha}(1+2^{\alpha})},\qquad
   \rho_n\ \ge\ \rho_{AI}(\alpha):=\frac{w_{\min}}{w_{\min}+4^{\alpha}/2}=\frac{\alpha}{\alpha+2}.
   \]
5. **Sharpened audit floor below \(\log3/\log5\).** For \(\alpha\in(0,\gamma)\), \(\gamma=\log3/\log5\),
   \[
   \|F(n,\alpha)\|\ \ge\ \sigma_{AI}(\alpha)\,\|F_{\mathrm{frac}}(n,\alpha)\|\ \ge\ \sigma_{AI}(\alpha)\,\frac{9\sqrt{21}}{28\,C_{AC}(\alpha)}\ >\ 0
   \]
   at every level. At \(\alpha=0.45\), with Theorem AA's \(\|F_{\mathrm{frac}}\|\ge0.084015\),
   \[
   \|F(n,0.45)\|\ \ge\ \sigma_{AI}(0.45)\cdot0.084015\ \approx\ 0.051988\qquad\text{for every }n\ge1,
   \]
   about \(80.7\) times Theorem AH's \(6.442\cdot10^{-4}\).

**Proof.** (1) Fix \(j\in\mathrm{Far}\); then \(j\notin\{a,p,q\}\), so \(w_{aj}\), \(w_{pj}\), \(w_{qj}\) are all off-diagonal and have the subordination form. Insert the first-step identity and exchange the order of integration (Tonelli; every integrand is nonnegative):
\[
w_{aj}=c_{\alpha}\int_0^{\infty}\bigl(H_s(p,j)+H_s(q,j)\bigr)\,K_{\alpha}(s)\,ds,\qquad
K_{\alpha}(s)=\int_s^{\infty}e^{-2(t-s)}\,t^{-1-\alpha}\,dt .
\]
For \(t\ge s\), \(t^{-1-\alpha}\le s^{-1-\alpha}\), so \(K_{\alpha}(s)\le s^{-1-\alpha}\int_0^{\infty}e^{-2u}\,du=s^{-1-\alpha}/2\). Hence \(w_{aj}\le\tfrac12c_{\alpha}\int_0^{\infty}(H_s(p,j)+H_s(q,j))s^{-1-\alpha}ds=\tfrac12(w_{pj}+w_{qj})\).

(2) As in AH (2), \(v=u_a-u\ge0\) (maximum principle), \(v_a=0\), and \(v\) is \(L^{\alpha}\)-harmonic at the interior vertices \(p\) and \(q\). Dropping \(j=a\) and the nonnegative cross terms (\(j=q\) at \(p\), \(j=p\) at \(q\)),
\[
(L^{\alpha})_{pp}v_p\ge\sum_{j\in\mathrm{Far}}w_{pj}v_j,\qquad (L^{\alpha})_{qq}v_q\ge\sum_{j\in\mathrm{Far}}w_{qj}v_j .
\]
Jensen gives \((L^{\alpha})_{pp},(L^{\alpha})_{qq}\le4^{\alpha}\), and \(I_a=(Lu)_a=v_p+v_q\). Add the two inequalities *before* comparing kernels, then apply (1) with \(v_j\ge0\):
\[
4^{\alpha}I_a\ \ge\ \sum_{j\in\mathrm{Far}}(w_{pj}+w_{qj})\,v_j\ \ge\ 2\sum_{j\in\mathrm{Far}}w_{aj}v_j=2I_{\mathrm{far}} .
\]

(3)–(5) are AH (3)–(5) verbatim with \(4^{\alpha}/(2\kappa)\) replaced by \(4^{\alpha}/2\): \(\gamma_n=2I_{\mathrm{far}}/(W_{\mathrm{far}}I_a)\le4^{\alpha}/W_{\mathrm{far}}\); \(I_a^{\alpha}=w_nI_a+I_{\mathrm{far}}\le(2^{\alpha-1}+4^{\alpha}/2)I_a\); \(\rho_n\ge w_n/(w_n+4^{\alpha}/2)\) increases in \(w_n\ge w_{\min}\); and \(\|F\|=\sigma_n\|F_{\mathrm{frac}}\|\) with AC's or AA's floor. The closed forms use \(4^{\alpha}\cdot4^{2-\alpha}=16\) and \(w_{\min}/(4^{\alpha}/2)=\alpha/2\). \(\square\)

Item (1) uses only \(\deg a=2\) and the two edges \(a\sim p\), \(a\sim q\); it is a statement about the *sum* of the two neighbor kernels, not about each one. The pointwise comparison \(w_{pj}\ge w_{aj}\) is not proved here (see the observation). For \(\alpha>\gamma\), items (1)–(4) still hold and (5) gives nothing, consistent with Theorem K.

**Witness, computed by `scripts/gasket_fractional_harnack_sum.py`.** At \(\alpha=0.45\):
\[
2^{\alpha-1}\approx0.683020,\quad
\frac{4^{\alpha}}{2}\approx0.933033,\quad
\Gamma_{AI}\approx32.3232\ (\Gamma_{AH}\approx4494.13=\Gamma_{AI}/\kappa),\quad
\sigma_{AI}\approx0.618792,\quad
\rho_{AI}=\tfrac{9}{49}\approx0.183673,
\]
and the audit floor is \(\approx0.051988\) (AA route) or \(\approx0.013646\) (AC route). At levels \(1\)–\(3\) (\(\alpha=0.45\)) and level \(2\) (\(\alpha=0.10,\,0.90\)) the script checks that row \(a\) of \(L\) is \(2e_a-e_p-e_q\), the backward equation at \(a\), the Duhamel formula for \(H_t(a,j)\) by quadrature, \(0<K_{\alpha}(s)\le s^{-1-\alpha}/2\) on \(s\in[10^{-30},10^{12}]\), and that the weighted integral \(c_{\alpha}\int(H_s(p,j)+H_s(q,j))K_{\alpha}(s)\,ds\) reproduces \(w_{aj}\) on Far (relative error \(<10^{-6}\)). Through level \(6\) at \(\alpha=0.45\): \(w_{pj}+w_{qj}\ge2w_{aj}\) for every far \(j\), \(v\ge0\), \(v\) harmonic at \(p,q\), the two dropped-term inequalities, \((L^{\alpha})_{pp}\le4^{\alpha}\), \(I_{\mathrm{far}}\le(4^{\alpha}/2)I_a\), \(\gamma_n\le4^{\alpha}/W_{\mathrm{far}}\le\Gamma_{AI}\), \(\sigma_n\ge\sigma_{AI}\), \(\rho_n\ge\rho_{AI}\), and \(\|F\|\ge\sigma_{AI}\|F_{\mathrm{frac}}\|\ge0.051988\); the recorded AE table of \(\gamma_n\) and the level-\(2\) lock \(1.165909\) are reproduced; the same items hold through level \(5\) at \(\alpha\in\{0.10,0.30,0.50,0.60,0.68,0.80,0.95\}\), with the AC-route floor below \(\gamma\). Live values at \(0.45\):
\[
\min_{j\in\mathrm{Far}}\frac{w_{pj}+w_{qj}}{2w_{aj}}\approx(2.002493,\;1.153701,\;1.034036,\;1.007326,\;1.001513,\;1.000306),\qquad
\frac{4^{\alpha}}{W_{\mathrm{far}}}\approx(4.3877,\;3.3542,\;3.1860,\;3.1566,\;3.1517,\;3.1509).
\]
The script exits nonzero if any check fails.

**KEPT FAILURE AI.1 (\(\kappa\to1\) does not reach Conjecture J.1).** The claim "replacing \(\kappa\) by \(1\) in Theorem AH and combining with Theorem AA's floor yields \(\|F(n,0.45)\|\ge0.70\), Conjecture J.1's value" is false. The chain gives \(\sigma_{AI}(0.45)\cdot0.084015\approx0.051988\). The route cannot be repaired by the Harnack constant alone: with constant \(0\) it would give only \(2^{1-\alpha}\cdot0.084015\approx0.1230\). Even with the live \(\|F_{\mathrm{frac}}(n)\|\) in place of AA's floor, \(\sigma_{AI}\|F_{\mathrm{frac}}(n)\|\approx(0.620534,\;0.527460,\;0.471420,\;0.437369,\;0.416191,\;0.402670)\) on levels \(1\)–\(6\), all below \(0.70\). Reaching \(0.70\) through \(\|F\|=\sigma_n\|F_{\mathrm{frac}}\|\) needs \(\sigma_n>1\) uniformly (live \(\sigma_n\approx1.352\)) together with a fractional-current floor much nearer the live \(\|F_{\mathrm{frac}}(6)\|\approx0.650736\); neither is proved.

**OBSERVATION, not a theorem.** Through level \(6\) at \(0.45\) (level \(5\) at the other \(\alpha\)), the pointwise ratio \(\min_jw_{pj}/w_{aj}\) stays \(\ge1\) and decreases toward \(1\) (\(1.000152\) at level \(6\)), so pointwise \(\kappa=1\) holds on the witness range but is not proved; and since the pointwise version would give exactly the same constant \(4^{\alpha}/2\), Theorem AI loses nothing by proving only the sum. The minimum sum ratio decreases toward \(1\), so the constant \(2\) in (1) is nearly attained on the range. The Harnack constant is still not sharp: live \(I_{\mathrm{far}}/I_a\approx0.407121\) at level \(6\) against \(4^{\alpha}/2\approx0.933033\); live \(\gamma_6\approx1.374864\) against \(\Gamma_{AI}\approx32.32\); and the proved floor \(0.051988\) is far below the live \(\|F(6,0.45)\|\approx0.880109\). No sharp constant is claimed.

**NOT THIS.** Theorem AI improves the constants of Theorem AH (\(\Gamma\) by the factor \(1/\kappa\approx139\), the audit floor at \(0.45\) from \(6.44\cdot10^{-4}\) to \(0.051988\)). It does not prove that \(\|F(n,0.45)\|\) converges, does not prove that it decreases, and does not prove Conjecture J.1's value \(L\ge0.70\) (Kept Failure AI.1). A positive floor on a graph audit norm is not thrust, not a force, and not a pass of any experimental criterion. Not thrust.

## Theorem AJ. Keeping the neighbor-edge weight sharpens the Harnack constant below \(4^{\alpha}/2\)

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\), used only for the numerical floors. Notation of Theorems AH–AI: \(a\) is the value-\(1\) corner (\(\deg a=2\)), \(p,q\) its two neighbors (\(\deg p=\deg q=4\), \(p\sim q\)), \(\mathrm{Far}=V\setminus\{a,p,q\}\), \(w_{xy}=-(L^{\alpha})_{xy}>0\) for \(x\ne y\), \(v=u_a-u\ge0\), \(I_a=(Lu)_a=v_p+v_q\), \(I_{\mathrm{far}}=\sum_{j\in\mathrm{Far}}w_{aj}v_j\), \(w_{pq}=-(L^{\alpha})_{pq}\). Theorems J, AA, AC, AF, AG, and AI are used as stated; this note does not re-prove them.

**Exact neighbor identity.** Since \(v\) is \(L^{\alpha}\)-harmonic at the interior vertices \(p\) and \(q\) and \(v_a=0\),
\[
(L^{\alpha})_{pp}v_p=w_{pq}v_q+\sum_{j\in\mathrm{Far}}w_{pj}v_j,\qquad
(L^{\alpha})_{qq}v_q=w_{pq}v_p+\sum_{j\in\mathrm{Far}}w_{qj}v_j.
\]
Adding and rearranging,
\[
\sum_{j\in\mathrm{Far}}(w_{pj}+w_{qj})\,v_j=(L^{\alpha})_{pp}v_p+(L^{\alpha})_{qq}v_q-w_{pq}(v_p+v_q).
\]
Theorem AI drops the nonnegative cross terms \(w_{pq}v_q\) and \(w_{pq}v_p\) before comparing kernels; here they are kept.

**THEOREM AJ.** For every \(\alpha\in(0,1)\) and every level \(n\ge1\):

1. **Neighbor-edge weight floor.** On the edge \(p\sim q\),
   \[
   w_{pq}\ \ge\ w_{\min}(\alpha):=\alpha\,4^{\alpha-1}
   \]
   (Theorem J: \(H_t(p,q)\ge te^{-4t}\) via \(M=4I-L\ge0\) entrywise with \(M_{pq}=1\), then subordination).
2. **Sharpened one-step Harnack at the neighbors.**
   \[
   I_{\mathrm{far}}\ \le\ \frac{4^{\alpha}-w_{pq}}{2}\,I_a\ \le\ H_{AJ}(\alpha)\,I_a,\qquad
   H_{AJ}(\alpha):=\frac{4^{\alpha}-\alpha\,4^{\alpha-1}}{2}=\frac{4^{\alpha}}{2}\Bigl(1-\frac{\alpha}{4}\Bigr).
   \]
   Strictly sharper than Theorem AI's \(4^{\alpha}/2\).
3. **Sharpened gap-ratio ceiling.** With AG's \(W_{\mathrm{far}}\ge W_{\mathrm{far}}^{\min}=2\alpha(1-\alpha)/4^{2-\alpha}\),
   \[
   \gamma_n\ \le\ \frac{4^{\alpha}-\alpha\,4^{\alpha-1}}{W_{\mathrm{far}}}\ \le\ \Gamma_{AJ}(\alpha):=\frac{2(4-\alpha)}{\alpha(1-\alpha)}=\Gamma_{AI}(\alpha)\Bigl(1-\frac{\alpha}{4}\Bigr).
   \]
4. **Sharpened share floors.** With AF's \(w_{\min}\le w_n\le2^{\alpha-1}\),
   \[
   \sigma_n\ \ge\ \sigma_{AJ}(\alpha):=\bigl(2^{\alpha-1}+H_{AJ}(\alpha)\bigr)^{-1}=\frac{2}{2^{\alpha}+4^{\alpha}-\alpha\,4^{\alpha-1}},\qquad
   \rho_n\ \ge\ \rho_{AJ}(\alpha):=\frac{w_{\min}}{w_{\min}+H_{AJ}(\alpha)}=\frac{2\alpha}{\alpha+4}.
   \]
5. **Sharpened audit floor below \(\log3/\log5\).** For \(\alpha\in(0,\gamma)\), \(\gamma=\log3/\log5\),
   \[
   \|F(n,\alpha)\|\ \ge\ \sigma_{AJ}(\alpha)\,\|F_{\mathrm{frac}}(n,\alpha)\|\ \ge\ \sigma_{AJ}(\alpha)\,\frac{9\sqrt{21}}{28\,C_{AC}(\alpha)}\ >\ 0
   \]
   at every level. At \(\alpha=0.45\), with Theorem AA's \(\|F_{\mathrm{frac}}\|\ge0.084015\),
   \[
   \|F(n,0.45)\|\ \ge\ \sigma_{AJ}(0.45)\cdot0.084015\ \approx\ 0.055599\qquad\text{for every }n\ge1,
   \]
   about \(1.069\) times Theorem AI's \(0.051988\).

**Proof.** (1) is Theorem J on the edge \(p\sim q\): \(H_t=e^{-4t}e^{tM}\) with \(M=4I-L\ge0\) entrywise and \(M_{pq}=1\), so \(H_t(p,q)\ge te^{-4t}\); subordination then gives \(w_{pq}\ge c_{\alpha}\int_0^{\infty}te^{-4t}\,t^{-1-\alpha}\,dt=\alpha\,4^{\alpha-1}\).

(2) From the exact identity, AI's sum \(w_{pj}+w_{qj}\ge2w_{aj}\) with \(v_j\ge0\), and Jensen \((L^{\alpha})_{pp},(L^{\alpha})_{qq}\le4^{\alpha}\),
\[
2I_{\mathrm{far}}\le\sum_{j\in\mathrm{Far}}(w_{pj}+w_{qj})\,v_j=(L^{\alpha})_{pp}v_p+(L^{\alpha})_{qq}v_q-w_{pq}I_a\le(4^{\alpha}-w_{pq})I_a.
\]
Item (1) upgrades the live \(w_{pq}\) to the uniform \(H_{AJ}\).

(3)–(5) are AI (3)–(5) with \(4^{\alpha}/2\) replaced by \(H_{AJ}\): \(\gamma_n=2I_{\mathrm{far}}/(W_{\mathrm{far}}I_a)\le(4^{\alpha}-\alpha4^{\alpha-1})/W_{\mathrm{far}}\); \(I_a^{\alpha}=w_nI_a+I_{\mathrm{far}}\le(2^{\alpha-1}+H_{AJ})I_a\); \(\rho_n\ge w_n/(w_n+H_{AJ})\) increases in \(w_n\ge w_{\min}\); and \(\|F\|=\sigma_n\|F_{\mathrm{frac}}\|\) with AC's or AA's floor. The closed forms use \(\Gamma_{AI}=8/(\alpha(1-\alpha))\) and \(H_{AJ}=(4^{\alpha}/2)(1-\alpha/4)\), so \(\Gamma_{AJ}=\Gamma_{AI}(1-\alpha/4)\) and \(w_{\min}/H_{AJ}=2\alpha/(4-\alpha)\), whence \(\rho_{AJ}=2\alpha/(\alpha+4)\). \(\square\)

For \(\alpha>\gamma\), items (1)–(4) still hold and (5) gives nothing, consistent with Theorem K. The uniform floor on \(w_{pq}\) does not vanish as \(n\to\infty\) (same constant as Theorem J's \(w_{\min}\) on every edge).

**Witness, computed by `scripts/gasket_fractional_wpq_harnack.py`.** At \(\alpha=0.45\):
\[
w_{\min}\approx0.209932,\quad
H_{AJ}\approx0.828067\ (4^{\alpha}/2\approx0.933033),\quad
\Gamma_{AJ}\approx28.6869\ (\Gamma_{AI}\approx32.3232),\quad
\sigma_{AJ}\approx0.661775,\quad
\rho_{AJ}=\tfrac{9}{44.5}\approx0.202247,
\]
and the audit floor is \(\approx0.055599\) (AA route) or \(\approx0.014594\) (AC route). At levels \(1\)–\(3\) (\(\alpha=0.45\)) and level \(2\) (\(\alpha=0.10,\,0.90\)) the script checks \(p\sim q\), \(M=4I-L\ge0\) with \(M_{pq}=1\), \(H_t(p,q)\ge te^{-4t}\) on a grid of \(t\), and that \(\int te^{-4t}\,t^{-1-\alpha}\,dt\) recovers \(\alpha4^{\alpha-1}\). Through level \(6\) at \(\alpha=0.45\): the exact neighbor identity, \(w_{pj}+w_{qj}\ge2w_{aj}\), \(w_{pq}\ge w_{\min}\), \(v\ge0\), \(v\) harmonic at \(p,q\), \(I_{\mathrm{far}}\le((4^{\alpha}-w_{pq})/2)I_a\le H_{AJ}I_a\), \(\gamma_n\le\Gamma_{AJ}\), \(\sigma_n\ge\sigma_{AJ}\), \(\rho_n\ge\rho_{AJ}\), and \(\|F\|\ge\sigma_{AJ}\|F_{\mathrm{frac}}\|\ge0.055599\); the recorded AE table of \(\gamma_n\) and the level-\(2\) lock \(1.165909\) are reproduced; the same items hold through level \(5\) at \(\alpha\in\{0.10,0.30,0.50,0.60,0.68,0.80,0.95\}\), with the AC-route floor below \(\gamma\). Live values at \(0.45\):
\[
w_{pq}\approx(0.371367,\;0.320597,\;0.310397,\;0.308566,\;0.308257,\;0.308206),\qquad
\frac{4^{\alpha}-w_{pq}}{2}\approx(0.747349,\;0.772734,\;0.777834,\;0.778750,\;0.778904,\;0.778930).
\]
The script exits nonzero if any check fails.

**KEPT FAILURE AJ.1 (keeping \(w_{pq}\) does not reach Conjecture J.1).** The claim "keeping the neighbor-edge cross term in the Harnack identity and combining with Theorem AA's floor yields \(\|F(n,0.45)\|\ge0.70\), Conjecture J.1's value" is false. The chain gives \(\sigma_{AJ}(0.45)\cdot0.084015\approx0.055599\). The route cannot be repaired by the Harnack constant alone: with constant \(0\) it would give only \(2^{1-\alpha}\cdot0.084015\approx0.1230\). Even with the live \(\|F_{\mathrm{frac}}(n)\|\) in place of AA's floor, \(\sigma_{AJ}\|F_{\mathrm{frac}}(n)\|\approx(0.663639,\;0.564099,\;0.504167,\;0.467750,\;0.445101,\;0.430641)\) on levels \(1\)–\(6\), all below \(0.70\). Reaching \(0.70\) through \(\|F\|=\sigma_n\|F_{\mathrm{frac}}\|\) still needs \(\sigma_n\) nearer the live \(\approx1.352\) together with a fractional-current floor nearer the live \(\|F_{\mathrm{frac}}(6)\|\approx0.650736\); neither is proved.

**OBSERVATION, not a theorem.** Through level \(6\) at \(0.45\), live \(w_{pq}\) decreases toward \(\approx0.308\) (well above \(w_{\min}\approx0.210\)), so the uniform improvement from keeping the cross term does not vanish on the range. The live constant \((4^{\alpha}-w_{pq})/2\approx0.779\) still sits far above live \(I_{\mathrm{far}}/I_a\approx0.407\); the proved floor \(0.055599\) is far below live \(\|F(6,0.45)\|\approx0.880109\). No sharp constant is claimed.

**NOT THIS.** Theorem AJ improves the constants of Theorem AI (\(H\) and \(\Gamma\) by the factor \(1-\alpha/4=0.8875\) at \(0.45\), the audit floor from \(0.051988\) to \(0.055599\)). It does not prove that \(\|F(n,0.45)\|\) converges, does not prove that it decreases, and does not prove Conjecture J.1's value \(L\ge0.70\) (Kept Failure AJ.1). A positive floor on a graph audit norm is not thrust, not a force, and not a pass of any experimental criterion. Not thrust.

## Conjecture J.1

**CONJECTURE, not a theorem.** For \(\alpha=0.45\) and these corner data, \(\|F(n)\|\) decreases for every \(n\ge 1\) and
\[
\lim_{n\to\infty}\|F(n)\| = L \quad\text{with}\quad L\ge 0.70.
\]
The only evidence is the computed sequence above: the norm is still decreasing at level \(7\), and the gap ratios \(\delta_{n+1}/\delta_n\) are increasing but still near \(0.67<1\). A geometric tail with ratio \(0.67\) would leave \(L\) near \(0.82\); a later rise of that ratio toward \(1\) can still push the limit to \(0\). Theorem AH now excludes \(0\): \(\|F(n,0.45)\|\ge 6.44\cdot10^{-4}\) at every level. Theorem AI sharpens that floor to \(\approx0.051988\); Theorem AJ sharpens it further to \(\approx0.055599\) by keeping the neighbor-edge weight \(w_{pq}\), still far below \(0.70\) (Kept Failures AI.1, AJ.1). Neither proves that the limit exists, that the norm decreases, or that \(L\ge 0.70\); those remain conjecture.

**Falsifier.** The first level \(n\) with \(\|F(n)\|<0.70\), or the first \(n\ge 7\) with \(\|F(n+1)\|>\|F(n)\|+10^{-8}\). Level \(8\) is the first level not required by the default witness. Through level \(10\) neither falsifier has fired: \(\|F(10)\|\approx 0.833640\), still decreasing.

A second, separate numerical guess, also not a theorem: at \(\alpha=0.10\), levels \(1\) through \(6\) give
\[
\|F\|\approx(1.50450860,\; 1.44384308,\; 1.42364897,\; 1.41577938,\; 1.41266163,\; 1.41143248),
\]
and the successive drops after level \(3\) contract by a factor \(<1/2\). **CONJECTURE J.2.** \(\|F(n,0.10)\|\ge 1.410\) for every \(n\). Falsified by the first level with \(\|F\|<1.410\). The \(\alpha\to 0\) profile \(u\equiv (a+b+c)/3\) off the corners has audit norm \(\sqrt{7/3}\approx 1.527525\), independent of \(n\ge 1\); \(\alpha=0.10\) is a different operator, and this sentence does not identify its limit with \(\sqrt{7/3}\).


## What this does not say

- **NOT A LIMIT THEOREM AT \(\alpha=0.45\).** Levels \(1\) through \(7\) still have increasing ratios, now past \(0.97\). Theorem J excludes \(+\infty\). Theorem K gives \(\|F\|\to 0\) only for \(\alpha>\log 3/\log 5\). Theorem L says energy collapse would force \(\|F\|\to 0\) at every \(\alpha\in(0,1)\), including \(0.45\), but does not prove energy collapse. Neither decides zero versus positive at \(\alpha=0.45\); Theorem AH does (not zero). Conjecture J.1 is not a theorem. Kept Failures J.1, J.2, J.3, and J.4 reject four sufficient conditions that would have closed it. Level \(10\) is recorded, not a limit. Theorem S settles the harmonic audit under uniform multiplicative weights (finite nonzero only at resistance \(\lambda=5/3\)); Kept Failure S.1 rejects inheriting that freeze for \(L^{0.45}\). Theorem W settles the multi-scale harmonic audit (finite for \(\lambda<5/3\); resistance diverges); Kept Failure W.1 rejects a multi-scale resistance freeze. Theorem X takes spectral \(L^{\alpha}\) on the multi-scale graph: neutrality and Schur hold, nested harmonicity is special to \(\alpha=1\); Kept Failures X.1–X.2 reject inheriting the \(\lambda=1\) freeze and reject coincidence with finest-edge fractional. Theorem AA gives the level-uniform floors \(Q_{0.45}\ge 0.128335\) and \(\|F_{\mathrm{frac}}\|\ge 0.084015\) (a floor exists at every \(\alpha<1/2\)); Kept Failure AA.1 closes the energy-collapse route to zero. Theorem AB factors \(\sigma_n=\|F\|/\|F_{\mathrm{frac}}\|=\rho_n/w_n\) and, with AA, equates \(\|F\|\to 0\) with \(\rho_n\to 0\); Kept Failure AB.1 rejects geometric collapse of \(\rho\). Theorem AC extends the energy floor to every \(\alpha<\log 3/\log 5\) by cell resistance and spectral mass, so AB's equivalence holds on that whole range; Kept Failure AC.1 (via Theorem K) rejects a floor above it. Theorem AD proves the exact pairing \(Q_{\alpha}R_{\alpha}=7/3\) and settles the endpoint: \(Q_{\gamma}(u_n)=\Theta(1/n)\) and \(\|F(n,\gamma)\|\to 0\); Kept Failure AD.1 rejects extending the AC floor to \(\gamma\). Theorem AE factors \(\rho_n=\lambda_n/(\lambda_n+(1-\lambda_n)\gamma_n)\) with sharp ceiling \(w_n\le8^{\alpha}/2\) and the resistance identity \(\rho_n=(\sqrt{21}/7)w_n\|F\|R_{\alpha}\); Kept Failure AE.1 rejects geometric wipeout of the neighbor-degree share \(\lambda_n\). Theorem AF floors that share, \(\lambda_n\ge\alpha\,2^{-1-\alpha}\) (Jensen: \(\ge\alpha\,2^{\alpha-1}\)), so for \(\alpha<\log3/\log5\) vanishing of \(\|F\|\) requires \(\gamma_n\to\infty\) (one-sided); Kept Failure AF.1 rejects geometric blowup of \(\gamma_n\). Theorem AG ceilings the same share, \(\lambda_n\le 1-\alpha(1-\alpha)/2^{3-\alpha}<1\), so \(\|F\|\to0\) if and only if \(\gamma_n\to\infty\) (two-sided); Kept Failure AG.1 rejects geometric climb of \(\lambda\) toward \(1\). Theorem AH compares the jump kernels of \(p\) and \(a\) (\(w_{pj}\ge\kappa(\alpha)w_{aj}\) off the corner patch) and applies one-step \(L^{\alpha}\)-harmonicity at the two neighbors, so \(\gamma_n\le\Gamma_{AH}(\alpha)<\infty\) and \(\|F(n,\alpha)\|\ge\sigma_{\min}(\alpha)\|F_{\mathrm{frac}}\|>0\) at every level for every \(\alpha<\log3/\log5\) (\(\ge 6.44\cdot10^{-4}\) at \(0.45\)); Kept Failure AH.1 rejects the max-principle shortcut to the same floor. Theorem AI replaces \(\kappa\) by \(1\) through the first-step identity at the degree-\(2\) corner (\(w_{pj}+w_{qj}\ge2w_{aj}\) on Far), so \(I_{\mathrm{far}}\le(4^{\alpha}/2)I_a\), \(\gamma_n\le8/(\alpha(1-\alpha))\), and \(\|F(n,0.45)\|\ge0.051988\) at every level; Kept Failure AI.1 rejects reading that as Conjecture J.1's \(0.70\). Theorem AJ keeps the neighbor-edge weight \(w_{pq}\ge\alpha4^{\alpha-1}\) in the Harnack identity, so \(I_{\mathrm{far}}\le H_{AJ}I_a\) with \(H_{AJ}=(4^{\alpha}-\alpha4^{\alpha-1})/2\), \(\gamma_n\le2(4-\alpha)/(\alpha(1-\alpha))\), and \(\|F(n,0.45)\|\ge0.055599\); Kept Failure AJ.1 rejects reading that as \(0.70\). Neither floor is a limit theorem: convergence of \(\|F(n,0.45)\|\) and Conjecture J.1 stay open.
- **NOT A RESIDUAL.** Their tilt diagnostic sets \((Lu)_i=1\) on the interior. Theorem F then says the corner currents sum to \(-(N-3)\), the source that was inserted. That sum grows like \(3^n\). It is not a force left after the source is removed.
- **FAILED as a continuum thrust for the combinatorial harmonic case.** Under Theorem G, \(\|F(n)\|\to 0\) for every fixed corner triple. The combinatorial harmonic limit does not keep a finite audit vector.
- **OPEN.** Whether \(\|F(n)\|\) for \(L^{0.45}\) converges, and to what value. It does not tend to \(0\) (Theorem AH: \(\|F(n,0.45)\|\ge 6.44\cdot10^{-4}\) at every level; Theorem AI sharpens this to \(\ge0.051988\); Theorem AJ to \(\ge0.055599\)). Divergence is excluded (Theorem J). Neutrality is settled (Theorem H). For \(\alpha>\log 3/\log 5\) the same audit norm tends to \(0\) (Theorem K). Energy controls the audit norm at every \(\alpha\in(0,1)\) (Theorem L), so energy collapse would force \(\|F\|\to 0\), and a positive audit-norm liminf would force a positive energy liminf; the energy does not vanish at \(\alpha=0.45\) (Theorem AA: \(Q_{0.45}(u_n)\ge 0.128335\) and \(\|F_{\mathrm{frac}}\|\ge 0.084015\) at every level, for the same reason at every \(\alpha<1/2\)), so energy collapse cannot send \(\|F\|\) to \(0\) (Kept Failure AA.1); Theorem AB equates vanishing of \(\|F\|\) with vanishing of the two-neighbor fractional share \(\rho_n\) (Kept Failure AB.1 rejects geometric collapse of \(\rho\)); Theorem AC proves the energy floor for every \(\alpha<\log 3/\log 5\), the exact threshold by Theorem K (Kept Failure AC.1), and Theorem AD settles \(\alpha=\log 3/\log 5\) itself: \(Q_{\gamma}(u_n)=\Theta(1/n)\to 0\) and \(\|F(n,\gamma)\|\to 0\) (Kept Failure AD.1); Theorem AE splits \(\rho_n\to0\) into degree-share and gap-ratio blowup and rejects geometric \(\lambda\) wipeout (Kept Failure AE.1); Theorem AF floors \(\lambda_n\) uniformly, so \(\|F\|\to0\) would force the far-to-neighbor gap ratio \(\gamma_n\to\infty\), and Kept Failure AF.1 rejects geometric \(\gamma\) blowup (live \(\gamma_n\) decreases through level \(6\)); Theorem AG ceilings \(\lambda_n\) uniformly below \(1\), so \(\|F\|\to0\) if and only if \(\gamma_n\to\infty\), and Kept Failure AG.1 rejects geometric climb of \(\lambda\) toward \(1\); Theorem AH ceilings \(\gamma_n\le\Gamma_{AH}(\alpha)\) by a kernel comparison and one-step Harnack at the neighbors, so \(\|F\|\) does not vanish for any \(\alpha<\log3/\log5\), and Kept Failure AH.1 rejects the max-principle shortcut; Theorem AI sharpens the Harnack constant from \(4^{\alpha}/(2\kappa)\) to \(4^{\alpha}/2\) by a two-neighbor heat-kernel sum, giving \(\gamma_n\le8/(\alpha(1-\alpha))\) and \(\|F(n,0.45)\|\ge0.051988\), and Kept Failure AI.1 rejects reaching \(0.70\) that way; Theorem AJ keeps \(w_{pq}\) to push the floor to \(0.055599\), and Kept Failure AJ.1 rejects reaching \(0.70\) that way; whether \(\|F(n,0.45)\|\) converges is OPEN. Geometric weights \(1/d^2\) on this build are settled (Theorem I): that audit norm diverges. Uniform multiplicative weights are settled for the harmonic case (Theorem S): finite nonzero only at the resistance value \(\lambda=5/3\); resistance weights do not freeze the fractional \(\alpha=0.45\) audit (Kept Failure S.1). Multi-scale multiplicative weights are settled for the harmonic case (Theorem W): finite nonzero for every \(\lambda<5/3\), and resistance diverges like \((n+1)\sqrt{21}/2\) (Kept Failure W.1). Spectral \(L^{\alpha}\) on that multi-scale graph is settled structurally by Theorem X (neutrality, Schur, \(\alpha=1\) recovers W, no nested fractional harmonicity); the \(\lambda=1\) freeze does not transfer to \(\alpha=0.45\) (Kept Failure X.1), and multi-scale fractional does not coincide with finest-edge fractional (Kept Failure X.2). None of this is thrust.
- **OPEN.** The Stage 2 question from the rectangle note, whether an informational tensor on the \(0.45\) mesh has nonzero integrated divergence after the Maxwell piece is removed. Theorem F says the harmonic combinatorial gasket does not supply that divergence. Under the Hessian reading of the rank-2 symbol, Theorem M in the README reduces the question to \(\Delta\Psi_{\mathrm{info}}\) on the enclosing surface (zero for static massless \(A_0\) in vacuum, Kept Failure M.1). Under the quadratic reading, Theorem N makes it a volume integral of \(\partial_i\Psi\,\Delta\Psi\), which for static massless \(A_0\) is the Maxwell piece itself and vanishes for an isolated device (Kept Failure N.1). Theorem O extends this to every two-derivative reading: the only survivor is a surface-dependent trace flux whose all-space total is zero (Kept Failure O.1). Theorem P extends it to every linear reading of any order: zero for shift-invariant ones (Kept Failure P.1), and a surface-dependent \(-\int_V E\) for the undifferentiated term (Kept Failure P.2). Theorem Q extends it to every shift-invariant polynomial reading of any order and degree: zero on every vacuum surface if the divergence vanishes on harmonic fields, otherwise surface-dependent (Kept Failure Q.1). Theorem R drops shift invariance and keeps the same dichotomy; the undifferentiated \(\delta^{ij}\Psi^2\) term does not even decay under dilation (Kept Failure R.1). Theorem T takes the same polynomial family under screening and still finds only surface-dependent mass remainders (Kept Failure T.1). Theorem U covers every first-derivative \(O(3)\)-covariant reading \(f(s)\delta+g(s)E\otimes E\) with \(C^1\) coefficients (non-polynomial allowed): conserved ones are Maxwell plus \(\kappa\delta\) and give \(\mathcal{G}=0\); all others are surface-dependent and decay like \(R^{-2}\) (Kept Failure U.1). Theorem V covers the \(\Psi+\nabla\Psi\) isotropic extension \(f(\Psi,s)\delta+g(\Psi,s)E\otimes E\): conserved ones remain Maxwell plus \(\kappa\delta\) with \(\mathcal{G}=0\); all others are surface-dependent (Kept Failure V.1). Theorem Y allows explicit position dependence and keeps the same conserved class (Kept Failure Y.1). Theorem Z allows the pseudotensor \(h\,\varepsilon\cdot\nabla\Psi\): the only new conserved piece is a curl with \(\mathcal{G}=0\) on every closed surface (Kept Failure Z.1). Higher-derivative non-polynomial readings and the mesh solve itself remain OPEN.

## Reproduce

```bash
python3 scripts/gasket_corner_current.py
python3 scripts/gasket_fractional_currents.py
python3 scripts/gasket_geometric_currents.py
python3 scripts/gasket_uniform_weights.py
python3 scripts/gasket_multiscale_weights.py
python3 scripts/gasket_multiscale_fractional.py
python3 scripts/gasket_fractional_limit.py
python3 scripts/gasket_fractional_energy_floor.py
python3 scripts/gasket_fractional_share.py
python3 scripts/gasket_fractional_capacity_floor.py
python3 scripts/gasket_fractional_endpoint.py
python3 scripts/gasket_fractional_degree_share.py
python3 scripts/gasket_fractional_lambda_floor.py
python3 scripts/gasket_fractional_lambda_ceiling.py
python3 scripts/gasket_fractional_gamma_ceiling.py
python3 scripts/gasket_fractional_harnack_sum.py
python3 scripts/gasket_fractional_wpq_harnack.py
python3 scripts/gasket_fractional_level9.py
python3 scripts/gasket_fractional_level10.py
```

The first script exits nonzero unless the level-1 triple is \((3/2,-6/5,-3/10)\) and levels 2, 3, 4 match the factor \((3/5)^{n-1}\), including one arbitrary corner triple at levels 1 and 2. The second exits nonzero unless fractional currents at \(\alpha=0.45\) stay neutral through level 4, the level-2 audit norm matches \(1.165909\) within \(5\cdot 10^{-3}\), the refinement ratios differ from \(3/5\) by more than \(0.15\), and the \(\alpha=1\) control recovers Theorem G. The third exits nonzero unless geometric \(1/d^2\) currents match \(4^n\) times the combinatorial ones through level 5, the audit norm equals \((12/5)^n\sqrt{21}/2\), and every successive ratio equals \(12/5\). The fourth exits nonzero unless uniform multiplicative weights recover the closed form \((\lambda\cdot 3/5)^n\sqrt{21}/2\) through level 5 (including the resistance constant \(\sqrt{21}/2\)), and the resistance-scaled fractional audit at \(\alpha=0.45\) is not constant (Kept Failure S.1). The fifth exits nonzero unless multi-scale weights recover \((\sqrt{21}/2)\sum_{k=0}^{n}(\lambda\cdot 3/5)^k\) through level 5, nested harmonicity holds, and resistance multi-scale grows like \((n+1)\sqrt{21}/2\) (Kept Failure W.1). The sixth exits nonzero unless multi-scale spectral \(L^{0.45}\) stays neutral and on the Schur line through level 6, recovers Theorem W at \(\alpha=1\), shows nested harmonicity fails for \(\alpha=0.45\), rejects the \(\lambda=1\) freeze inheritance (Kept Failure X.1), and rejects coincidence with finest-edge fractional (Kept Failure X.2). The seventh exits nonzero unless, through level \(7\), the \(\alpha=0.45\) audit norm matches the level-2 lock \(1.165909\), stays neutral and inside the Theorem J cap, breaks both the \(0.95\) ratio law and the \(5^{-0.45}\) gap contraction, obeys the Theorem L energy sandwich (\(\|F_{\mathrm{frac}}\|\ge w_{\min}\|F\|\) and \(\|F\|\le\kappa\sqrt{Q}\)), breaks the \(0.95\) decay law for \(Q_{0.45}\), and the \(\alpha=0.9\) norm stays under the explicit Theorem K bound. `python3 scripts/gasket_fractional_limit.py 8` recomputes one level higher; the default witness does not. The eighth script exits nonzero unless the closed-form path resistance matches a dense pseudo-inverse, stays below \(C(\alpha)\) through \(m=2^{22}\) for \(\alpha<1/2\), and Löwner–Heinz, path duality, the energy floor and ceiling, and the pairing \(\|F_{\mathrm{frac}}\|=Q\sqrt{21}/7\) hold live on levels \(1\)–\(6\) (Theorem AA, Kept Failure AA.1). The ninth script exits nonzero unless altitude reflection equalizes the two neighbor weights, \(\sigma=\rho/w_n\) holds with \(w_{\min}\le w_n\le 8^{\alpha}\), the sandwich and Schur identities hold through level \(6\), and \(\rho_{5}/\rho_{4}>0.99\) (Theorem AB, Kept Failure AB.1). The tenth script exits nonzero unless the cell corner resistance is \((2/3)(5/3)^k\) with diameter under \((8/3)(5/3)^k\), \(\lambda_1\ge 5^{-n}/4\), the spectral mass of \(e_a-e_b\) stays under \(80s^{\gamma}\), and \(R_{\alpha}(n)\le C_{AC}(\alpha)\) with \(Q_{\alpha}(u_n)\ge(9/4)/R_{\alpha}(n)\) through level \(6\) for \(\alpha\in[0.45,0.68]\), and \(Q_{0.75}\), \(Q_{0.9}\) fall under the Theorem K bound (Theorem AC, Kept Failure AC.1). The eleventh script exits nonzero unless the corner Schur complement is \(\kappa(3I-J)\) with \(Q_{\alpha}R_{\alpha}=7/3\), the corner bump energy is \(2(3/5)^k\), \(R^{s}\ge(5s)^{\gamma-1}/5\) on \([5^{-n},1)\) through level \(7\), the Stieltjes quadrature matches \(R_{\gamma}(3)\), and \(R_{\gamma}(n)\), \(Q_{\gamma}(u_n)\), \(\|F(n,\gamma)\|\) obey the Theorem AD sandwich through level \(7\) while the Theorem K majorant stays constant (Theorem AD, Kept Failure AD.1). The twelfth script exits nonzero unless \(w_n\le8^{\alpha}/2\), the factorization \(\rho=\lambda/(\lambda+(1-\lambda)\gamma)\) and the AD resistance identities hold through level \(6\), and \(\lambda_4/\lambda_3>0.99\) (Theorem AE, Kept Failure AE.1). The thirteenth script exits nonzero unless \((L^{\alpha})_{aa}\le2^{\alpha}\le8^{\alpha}\), \(\lambda_n\ge\alpha\,2^{\alpha-1}\ge\lambda_{\min}(\alpha)\), \(w_n\le2^{\alpha-1}\), and \(\rho_n\ge1/(1+C_J\gamma_n)\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at six other \(\alpha\)), AE's factorization still holds, and every \(\gamma_{n+1}/\gamma_n<1\) (Theorem AF, Kept Failure AF.1). The fourteenth script exits nonzero unless \(W_{\mathrm{far}}\ge 2\alpha(1-\alpha)/4^{2-\alpha}\), \(\lambda_n\le 1-\alpha(1-\alpha)/2^{3-\alpha}\), the sandwich \(1/(1+C_J\gamma_n)\le\rho_n\le 1/(1+c_*\gamma_n)\), and the level-independent distance-\(2\) path counts hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at six other \(\alpha\)), and every \(\lambda_{n+1}/\lambda_n<1\) (Theorem AG, Kept Failure AG.1). The fifteenth script exits nonzero unless the subordination quadrature, the small-time bound \(\int_0^s\le K\int_s^{2s}\), \(H_s(p,a)\ge se^{-4s}\), and the kernel chain \(w_{pj}\ge\kappa w_{aj}\) hold at levels \(1\)–\(3\), and \(v=u_a-u\ge0\), \((L^{\alpha})_{pp}\le4^{\alpha}\), \(I_{\mathrm{far}}\le(4^{\alpha}/(2\kappa))I_a\), \(\gamma_n\le\Gamma_{AH}(\alpha)\), \(\sigma_n\ge\sigma_{\min}\), and \(\|F\|\ge\sigma_{\min}\|F_{\mathrm{frac}}\|\ge6.44\cdot10^{-4}\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the max-principle shortcut bound is negative at levels \(3\)–\(6\) (Theorem AH, Kept Failure AH.1). The sixteenth script exits nonzero unless the backward equation and Duhamel formula at the degree-\(2\) corner, \(0<K_{\alpha}(s)\le s^{-1-\alpha}/2\), and the weighted subordination identity for \(w_{aj}\) hold at levels \(1\)–\(3\), and \(w_{pj}+w_{qj}\ge2w_{aj}\) on Far, \(I_{\mathrm{far}}\le(4^{\alpha}/2)I_a\), \(\gamma_n\le\Gamma_{AI}(\alpha)\), \(\sigma_n\ge\sigma_{AI}\), \(\rho_n\ge\rho_{AI}\), and \(\|F\|\ge\sigma_{AI}\|F_{\mathrm{frac}}\|\ge0.051988\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the \(\kappa\to1\) chain stays below \(0.70\) (Theorem AI, Kept Failure AI.1). The seventeenth script exits nonzero unless \(H_t(p,q)\ge te^{-4t}\), \(w_{pq}\ge\alpha4^{\alpha-1}\), the exact neighbor identity with cross term, \(I_{\mathrm{far}}\le((4^{\alpha}-w_{pq})/2)I_a\le H_{AJ}I_a\), \(\gamma_n\le\Gamma_{AJ}\), \(\sigma_n\ge\sigma_{AJ}\), \(\rho_n\ge\rho_{AJ}\), and \(\|F\|\ge\sigma_{AJ}\|F_{\mathrm{frac}}\|\ge0.055599\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the \(w_{pq}\) chain stays below \(0.70\) (Theorem AJ, Kept Failure AJ.1). The level-\(9\) and level-\(10\) scripts replay recorded matrix-free values against a live level-\(2\) lock; the level-\(10\) script exits nonzero unless \(\delta_{10}/\delta_9>\log 3/\log 5\) (Kept Failure J.4), the norm still decreases above \(0.70\), and the Schur line, Theorem J cap, and Theorem L sandwich hold. To recompute level \(10\) itself (about half an hour; needs SciPy): `python3 scripts/gasket_fractional_matrix_free.py --level 10`.
