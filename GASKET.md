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

Theorems AF–AJ bound the diagonal entries \((L^{\alpha})_{aa}\) and \((L^{\alpha})_{pp}\) by Jensen with the first moment only. Here the second moment and a sharp spectral ceiling are used as well.

**THEOREM AK.** For every \(\alpha\in(0,1)\) and every level \(n\ge1\):

0. **Spectral ceiling.** Let \(N\) be the incidence matrix of the \(3^n\) level-\(n\) cells (upward triangles) against vertices, and \(c_x\in\{1,2\}\) the number of cells containing \(x\). Then \(L=3\,\mathrm{diag}(c)-N^{\mathsf T}N\), so \(\mathrm{spec}(L)\subset[0,6]\). On the checked levels \(2\)–\(4\) the ceiling is attained (\(6\) is an eigenvalue); attainment at every level is not claimed.
1. **Two-moment ceiling.** If \(\mu\) is a probability measure on \([0,6]\) with mean \(m<6\) and variance \(s^2\), then
   \[
   \int\lambda^{\alpha}\,d\mu\ \le\ \pi\,6^{\alpha}+(1-\pi)\,t^{\alpha},\qquad t=m-\frac{s^2}{6-m},\quad \pi=\frac{s^2}{s^2+(6-m)^2}.
   \]
   The spectral measure of \(e_a\) has \(m=L_{aa}=2\) and \(s^2=(L^2)_{aa}-4=2\); that of \(e_p\) has \(m=4\) and \(s^2=(L^2)_{pp}-16=4\). Hence
   \[
   D_a:=(L^{\alpha})_{aa}\ \le\ D_a^{+}(\alpha):=\tfrac19\,6^{\alpha}+\tfrac89\bigl(\tfrac32\bigr)^{\alpha}\ <\ 2^{\alpha},\qquad
   (L^{\alpha})_{pp}=(L^{\alpha})_{qq}\ \le\ D_p^{+}(\alpha):=\frac{6^{\alpha}+2^{\alpha}}{2}\ <\ 4^{\alpha}.
   \]
2. **Neighbor-weight ceiling.** \(w_n=(D_a-W_{\mathrm{far}})/2\le w_n^{+}(\alpha):=(D_a^{+}-W_{\mathrm{far}}^{\min})/2<2^{\alpha-1}\), with AG's \(W_{\mathrm{far}}^{\min}=2\alpha(1-\alpha)/4^{2-\alpha}\).
3. **Sharpened Harnack and gap ratio.** \(I_{\mathrm{far}}\le\frac{(L^{\alpha})_{pp}-w_{pq}}{2}I_a\le H_{AK}(\alpha)\,I_a\) with \(H_{AK}:=(D_p^{+}-\alpha4^{\alpha-1})/2<H_{AJ}\), and \(\gamma_n\le\Gamma_{AK}(\alpha):=(D_p^{+}-\alpha4^{\alpha-1})/W_{\mathrm{far}}^{\min}<\Gamma_{AJ}\).
4. **Sharpened share floors.** \(\sigma_n\ge\sigma_{AK}(\alpha):=1/(w_n^{+}+H_{AK})>\sigma_{AJ}\) and \(\rho_n\ge\rho_{AK}(\alpha):=w_{\min}/(w_{\min}+H_{AK})\).
5. **Sharpened audit floor below \(\log3/\log5\).** For \(\alpha<\gamma\), \(\|F(n,\alpha)\|\ge\sigma_{AK}(\alpha)\cdot9\sqrt{21}/(28\,C_{AC}(\alpha))>0\). At \(\alpha=0.45\), with Theorem AA's \(\|F_{\mathrm{frac}}\|\ge0.084015\),
   \[
   \|F(n,0.45)\|\ \ge\ \sigma_{AK}(0.45)\cdot0.084015\ \approx\ 0.058941\qquad\text{for every }n\ge1,
   \]
   about \(1.060\) times Theorem AJ's \(0.055599\).

**Proof.** (0) Each edge lies in exactly one cell, two vertices share at most one cell, and they are adjacent exactly when they share one, so \(N^{\mathsf T}N=\mathrm{diag}(c)+A\). Every cell gives each of its vertices degree \(2\), so \(\deg=2c\) and \(L=2\,\mathrm{diag}(c)-A=3\,\mathrm{diag}(c)-N^{\mathsf T}N\le3\,\mathrm{diag}(c)\le6\). (This is the line-graph bound \(\lambda_{\min}(A)\ge-2\).)

(1) Let \(f(\lambda)=\lambda^{\alpha}\) and let \(g\) be the quadratic with \(g(6)=f(6)\), \(g(t)=f(t)\), \(g'(t)=f'(t)\). Hermite's remainder gives \(f(x)-g(x)=\tfrac16f^{(3)}(\xi)(x-t)^2(x-6)\) for some \(\xi\in(0,6)\), and \(f^{(3)}=\alpha(\alpha-1)(\alpha-2)\lambda^{\alpha-3}>0\), so \(f\le g\) on \((0,6]\), and at \(0\) by continuity. Since \(g\) is quadratic, \(\int g\,d\mu\) depends only on \(m\) and \(s^2\), and equals \(\int g\,d\nu=\int f\,d\nu\) for the two-point law \(\nu\) on \(\{t,6\}\) with that mean and variance; that law has the stated \(t,\pi\). Both \(t\ge0\) here (\(t=3/2\), \(t=2\)). The local moments are \(L_{xx}=d_x\) and \((L^2)_{xx}=d_x^2+d_x\). Strictness against \(2^{\alpha}\) and \(4^{\alpha}\) is Jensen applied to \(\nu\).

(2) Row sums of \(L^{\alpha}\) vanish and \(w_{ap}=w_{aq}=w_n\) (altitude reflection), so \(D_a=2w_n+W_{\mathrm{far}}\).

(3)–(5) are AJ (2)–(5) with \(4^{\alpha}\) replaced by \(D_p^{+}\) in the exact neighbor identity and \(2^{\alpha-1}\) by \(w_n^{+}\) in \(I_a^{\alpha}=w_nI_a+I_{\mathrm{far}}\). \(\square\)

**Witness, computed by `scripts/gasket_fractional_moment_jensen.py`.** At \(\alpha=0.45\):
\[
D_a^{+}\approx1.315657\ (2^{\alpha}\approx1.366040),\quad
D_p^{+}\approx1.802814\ (4^{\alpha}\approx1.866066),\quad
w_n^{+}\approx0.628963\ (2^{\alpha-1}\approx0.683020),
\]
\[
H_{AK}\approx0.796441\ (H_{AJ}\approx0.828067),\quad
\Gamma_{AK}\approx27.5912\ (\Gamma_{AJ}\approx28.6869),\quad
\sigma_{AK}\approx0.701556\ (\sigma_{AJ}\approx0.661775),\quad
\rho_{AK}\approx0.208603,
\]
and the audit floor is \(\approx0.058941\) (AA route) or \(\approx0.015472\) (AC route). At levels \(1\)–\(4\) the script checks the \(3^n\) cells, \(\deg=2c\), \(L=3\,\mathrm{diag}(c)-N^{\mathsf T}N\) exactly, \(\lambda_{\max}\le6\) (equal to \(6\) at levels \(2\)–\(4\)), the local moments \(2,6\) and \(4,20\), that the Hermite quadratic dominates \(\lambda^{\alpha}\) on \(\mathrm{spec}(L)\), and that its spectral average equals the two-point value, at \(\alpha\in\{0.10,0.45,0.68,0.90\}\). Through level \(6\) at \(\alpha=0.45\): \(D_a\le D_a^{+}\), \(D_p=D_q\le D_p^{+}\), \(w_n\le w_n^{+}\), \(I_{\mathrm{far}}\le((D_p-w_{pq})/2)I_a\le H_{AK}I_a\), \(\gamma_n\le\Gamma_{AK}\), \(\sigma_n\ge\sigma_{AK}\), \(\rho_n\ge\rho_{AK}\), and \(\|F\|\ge\sigma_{AK}\|F_{\mathrm{frac}}\|\ge0.058941\); the same hold through level \(5\) at \(\alpha\in\{0.10,0.30,0.50,0.60,0.68,0.80,0.95\}\), with the AC-route floor below \(\gamma\). Live at \(0.45\), level \(6\): \(D_a\approx1.256755\), \(D_p\approx1.742623\), \(w_n\approx0.332260\). The script exits nonzero if any check fails.

**KEPT FAILURE AK.1 (the one-step neighbor-Harnack chain cannot reach Conjecture J.1, even with live constants).** The claim "the chain \(\sigma_n\ge1/\bigl(w_n+((L^{\alpha})_{pp}-w_{pq})/2\bigr)\) used in Theorems AI–AK reaches \(\|F(n,0.45)\|\ge0.70\) once its constants are made sharp" is false. Replacing every constant by its live value (live \(w_n\), live \((L^{\alpha})_{pp}\), live \(w_{pq}\)) gives chain values \(\sigma\approx(0.959907,\;0.953175,\;0.952875,\;0.952863,\;0.952863,\;0.952863)\) on levels \(1\)–\(6\), all below \(1\), while the live \(\sigma_n\approx1.3525\). Multiplied by the live \(\|F_{\mathrm{frac}}(n)\|\), the chain gives \(\approx(0.962609,\;0.812489,\;0.725939,\;0.673494,\;0.640883,\;0.620062)\), below \(0.70\) from level \(4\) on. The loss sits in the two inequalities the chain cannot avoid: \(w_{pj}+w_{qj}\ge2w_{aj}\) on Far and dropping the Far-side structure of \(v\). Sharpening constants inside this family therefore cannot prove J.1; a route to \(0.70\) needs a different comparison than one-step harmonicity at \(p,q\).

**NOT THIS.** Theorem AK improves the diagonal ceilings of AF–AJ and lifts the audit floor from \(0.055599\) to \(0.058941\). It does not prove that \(\|F(n,0.45)\|\) converges, does not prove that it decreases, and does not prove Conjecture J.1's value \(L\ge0.70\) (Kept Failure AK.1). The spectral ceiling \(6\) is a property of the graph Laplacian, not of any device. A positive floor on a graph audit norm is not thrust, not a force, and not a pass of any experimental criterion. Not thrust.

## Theorem AL. Hermite–Hadamard sharpens the path-resistance ceiling

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\). Notation of Theorem AA: \(P_m\) is the path of \(m\) unit edges along the gasket side from corner \(a\) to corner \(b\), \(R_{\alpha}(m)=v^{T}(P_m^{\alpha})^{+}v\) with \(v=e_0-e_m\), and
\[
C(\alpha)=4^{-\alpha}\Bigl(8+\tfrac4\pi B\bigl(\tfrac12-\alpha,\tfrac32\bigr)\Bigr),\qquad
R_{\infty}(\alpha)=4^{-\alpha}\tfrac4\pi B\bigl(\tfrac12-\alpha,\tfrac32\bigr).
\]
Löwner–Heinz and the path duality of Theorem AA are used as stated. Theorems AA and AK are used for the numerical floors; this note does not re-prove them.

**THEOREM AL.** Let \(\alpha_*=(14-8\sqrt{2})/34\approx 0.079009\). For every \(\alpha\in[\alpha_*,1/2)\) and every path length \(m\ge 1\),
\[
R_{\alpha}(m)\ \le\ R_{\infty}(\alpha)\ <\ C(\alpha).
\]
Consequently, in the setting of Theorem AA,
\[
Q_{\alpha}(u_n)\ \ge\ \frac{9/4}{R_{\infty}(\alpha)},\qquad
\|F_{\mathrm{frac}}(n,\alpha)\|\ \ge\ \frac{9\sqrt{21}}{28\,R_{\infty}(\alpha)}.
\]
At \(\alpha=0.45\): \(R_{\infty}\approx 13.245199\), so \(Q_{0.45}(u_n)\ge 0.169873\) and \(\|F_{\mathrm{frac}}(n,0.45)\|\ge 0.111208\) at every level (about \(1.324\) times Theorem AA's \(0.084015\)). Combined with Theorem AK's share floor \(\sigma_{AK}(0.45)\approx 0.701556\),
\[
\|F(n,0.45)\|\ \ge\ \sigma_{AK}(0.45)\cdot 0.111208\ \approx\ 0.078019\qquad\text{for every }n\ge 1,
\]
about \(1.324\) times Theorem AK's \(0.058941\).

**Proof.** Write \(\Delta=\pi/(m+1)\) and \(f(\theta)=(2\sin\theta)^{-2\alpha}\cos^2\theta\). Theorem AA's cosine sum is the midpoint Riemann sum
\[
R_{\alpha}(m)=\frac8\pi\sum_{j=0}^{\lfloor(m-1)/2\rfloor}f\bigl((j+\tfrac12)\Delta\bigr)\,\Delta.
\]
The summands' intervals of width \(\Delta\) lie in \((0,\pi/2]\) (when \(m\) is odd they tile \([0,\pi/2]\) exactly; when \(m\) is even they cover a proper initial segment). It therefore suffices to show \(f\) is convex on \((0,\pi/2)\), for then Hermite–Hadamard on each such interval gives \(f(\mathrm{mid})\le\) average, hence \(R_{\alpha}(m)\le(8/\pi)\int_0^{\pi/2}f=R_{\infty}(\alpha)\).

Set \(b=2\alpha\in(0,1)\) and \(g(\theta)=\sin^{-b}\theta\,\cos^2\theta\), so \(f=2^{-b}g\). Direct differentiation yields
\[
\frac{g''}{g}=\frac{b(b+1)}{\sin^2\theta}+\frac{2}{\cos^2\theta}+4b-4-b^2
\]
wherever \(g>0\). In the coordinate \(x=\sin^2\theta\in(0,1)\) this is \(h(x)=b(b+1)/x+2/(1-x)+4b-4-b^2\). The critical point of the \(x\)-dependent part is \(x_*=\sqrt{b(b+1)}/(\sqrt{b(b+1)}+\sqrt{2})\), and
\[
\min_{x\in(0,1)}h(x)=5b-2+2\sqrt{2b(b+1)}.
\]
The right-hand side is nonnegative precisely when \(b\ge(14-8\sqrt{2})/17\), i.e. \(\alpha\ge\alpha_*\) (for \(b\ge 2/5\) it is immediate; for \(b\in[\alpha_*\cdot 2,\,2/5]\) the inequality \(4\cdot 2b(b+1)\ge(2-5b)^2\) rearranges to \(17b^2-28b+4\le 0\), whose relevant root is \(b_*\)). Thus \(g''\ge 0\) and \(f\) is convex on \((0,\pi/2)\) for every \(\alpha\in[\alpha_*,1/2)\). The Löwner–Heinz and duality steps of Theorem AA, with \(C(\alpha)\) replaced by \(R_{\infty}(\alpha)\), give the energy and fractional-current floors. The audit-norm bound is Theorem AK (4)–(5) with AA's \(0.084015\) replaced by \(0.111208\). \(\square\)

**Witness, computed by `scripts/gasket_fractional_path_ceiling.py`.** At \(\alpha=0.45\): \(R_{\infty}\approx 13.245199\) against \(C\approx 17.532293\) (ratio \(\approx 1.324\)); \(Q\)-floor \(\approx 0.169873\); \(\|F_{\mathrm{frac}}\)-floor \(\approx 0.111208\); \(\|F\|\)-floor via \(\sigma_{AK}\) \(\approx 0.078019\). The analytic second-derivative minimum equals \(5b-2+2\sqrt{2b(b+1)}\) and vanishes at \(b_*=2\alpha_*\); it is positive at \(b=0.9\) (\(\approx 6.199\)). Through \(m\le 64\) and \(m=2^k\) up to \(2^{11}\), \(R_{\alpha}(m)\le R_{\infty}(\alpha)\) at \(\alpha\in\{0.10,0.25,0.45,0.49\}\), the closed form matches a dense pseudo-inverse, and odd \(m\) give an exact midpoint cover of \([0,\pi/2]\). Live on levels \(1\)–\(6\) at \(0.45\): \(Q\), \(\|F_{\mathrm{frac}}\|\), and \(\|F\|\) dominate the new floors, and \(Q\) dominates the path energy \((9/4)/R(2^n)\), which dominates \(0.169873\). Level \(4\) clears the new floors at \(\alpha\in\{0.10,0.25,0.49\}\). The script exits nonzero if any check fails.

**KEPT FAILURE AL.1 (the sharp path ceiling does not reach Conjecture J.1).** The claim "replacing Theorem AA's \(C(\alpha)\) by \(R_{\infty}(\alpha)\) and combining with Theorem AK's share floor yields \(\|F(n,0.45)\|\ge 0.70\)" is false. The chain gives \(\sigma_{AK}\cdot 0.111208\approx 0.078019\). Even the absolute one-step envelope \(2^{1-\alpha}\cdot 0.111208\approx 0.1628\), or the live share \(\sigma_n\approx 1.3525\) times the new fractional floor (\(\approx 0.1504\)), stays below \(0.70\). The obstruction is the path comparison itself: live gasket resistance at \(0.45\) is near \(2.35\), while \(R_{\infty}\approx 13.25\), so no sharpening of the path ceiling alone can push \(\|F_{\mathrm{frac}}\|\) near the live \(\approx 0.65\) needed for J.1 with a bounded share. This route is different from the one-step neighbor-Harnack chain closed in Kept Failure AK.1; both fall short of \(0.70\).

**NOT THIS.** Theorem AL sharpens the uniform path-resistance ceiling of Theorem AA for \(\alpha\ge\alpha_*\) and lifts the audit floor from \(0.058941\) to \(0.078019\). It does not cover \(\alpha\in(0,\alpha_*)\) by this argument (there \(f''\) changes sign, though AA's \(C(\alpha)\) still applies). It does not prove that \(\|F(n,0.45)\|\) converges, does not prove that it decreases, and does not prove Conjecture J.1's value \(L\ge 0.70\) (Kept Failure AL.1). A positive floor on a graph audit norm is not thrust, not a force, and not a pass of any experimental criterion. Not thrust.

## Theorem AM. Boundary-cycle resistance ceiling via Fourier and the heat kernel

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\) on corners \(a,b,c\). Notation of Theorems AA and AL: \(R_{\infty}(\alpha)=4^{-\alpha}\tfrac4\pi B(\tfrac12-\alpha,\tfrac32)\) is AL's path ceiling. On the level-\(n\) gasket let \(m=2^n\), \(N=3m\), and let \(M\) be the Laplacian of the \(3m\) outer boundary edges only (both endpoints on one outer side; weight \(1\)), padded by zero on the interior vertices. The boundary is a cycle \(C_N\) on which \(a,b,c\) sit at positions \(0,m,2m\). Write \(\lambda_k=4\sin^2(\pi k/N)\) and
\[
R^{\mathrm{cyc}}_{\alpha}(m)=(e_a-e_b)^{T}(M^{\alpha})^{+}(e_a-e_b),\qquad
R^{\mathrm{cyc}}_{\infty}(\alpha)=2\int_0^1(2\sin\pi x)^{-2\alpha}\,dx .
\]
Löwner–Heinz is used as in Theorem AA; the corner pairing is the mechanism of Theorem AD (i); Theorem AK's share floor \(\sigma_{AK}\) is used as stated. This note does not re-prove them.

**THEOREM AM.** For every \(\alpha\in(0,1/2)\) and every \(m\ge 1\),
\[
R^{\mathrm{cyc}}_{\alpha}(m)=\frac2N\sum_{k=1}^{N-1}\lambda_k^{-\alpha}\Bigl(1-\cos\frac{2\pi k}{3}\Bigr)
=\frac3N\sum_{\substack{1\le k\le N-1\\ 3\nmid k}}\lambda_k^{-\alpha}
\ <\ R^{\mathrm{cyc}}_{\infty}(\alpha)=\frac2\pi\,4^{-\alpha}B\bigl(\tfrac12-\alpha,\tfrac12\bigr)=(1-\alpha)\,R_{\infty}(\alpha).
\]
Consequently, at every level \(n\ge1\),
\[
Q_{\alpha}(u_n)\ \ge\ u_n^{T}M^{\alpha}u_n\ \ge\ \frac{7/3}{R^{\mathrm{cyc}}_{\alpha}(2^n)}\ >\ \frac{7/3}{(1-\alpha)R_{\infty}(\alpha)},\qquad
\|F_{\mathrm{frac}}(n,\alpha)\|\ \ge\ \frac{\sqrt{21}}{3\,(1-\alpha)R_{\infty}(\alpha)} .
\]
The ceiling holds on all of \((0,1/2)\), with no \(\alpha_*\) threshold. At \(\alpha=0.45\): \(R^{\mathrm{cyc}}_{\infty}\approx7.284859\) (against AL's \(13.245199\)), so \(Q_{0.45}(u_n)\ge0.320299\) and \(\|F_{\mathrm{frac}}(n,0.45)\|\ge0.209685\) at every level, which is \((28/27)/(1-\alpha)\approx1.886\) times AL's \(0.111208\). The two-corner test vector alone gives \((9/4)/R^{\mathrm{cyc}}_{\infty}\approx0.308860\). Combined with Theorem AK's share floor \(\sigma_{AK}(0.45)\approx0.701556\),
\[
\|F(n,0.45)\|\ \ge\ \sigma_{AK}(0.45)\cdot0.209685\ \approx\ 0.147106\qquad\text{for every }n\ge1,
\]
about \(1.886\) times Theorem AL's \(0.078019\).

**Proof.** *(i) Embedding.* \(L-M\) is a sum of the PSD edge matrices of the non-boundary edges, so \(0\le M\le L\), and Löwner–Heinz gives \(M^{\alpha}\le L^{\alpha}\) for \(\alpha\in(0,1)\). Hence \(Q_{\alpha}(u_n)=u_n^{T}L^{\alpha}u_n\ge u_n^{T}M^{\alpha}u_n\). The kernel of \(M\) is spanned by the interior indicators and the boundary constant, so every corner-supported zero-sum \(v\) lies in the range of \(M^{\alpha}\), and Cauchy–Schwarz in the \(M^{\alpha}\) form gives \(u^{T}M^{\alpha}u\ge (v^{T}u)^2/(v^{T}(M^{\alpha})^{+}v)\).

*(ii) Corner pairing.* The dihedral symmetry of \(C_{3m}\) that fixes the set \(\{0,m,2m\}\) permutes \(a,b,c\) in every way and commutes with \(M^{\alpha}\), so the corner block of \((M^{\alpha})^{+}\) is \(g_0I+g_1(J-I)\) and \(R^{\mathrm{cyc}}=2(g_0-g_1)\). Maximizing over zero-sum \(v\) on the corners gives \(2\lvert u_C-\bar u_C\rvert^2/R^{\mathrm{cyc}}\). For \((1,-1/2,0)\), \(\bar u_C=1/6\) and \(\lvert u_C-\bar u_C\rvert^2=7/6\), so the floor is \((7/3)/R^{\mathrm{cyc}}\). This is the pairing of Theorem AD (i), applied to \(M^{\alpha}\); it equals the minimum of \(w^{T}M^{\alpha}w\) over boundary extensions \(w\) of the corner data. The single test vector \(e_a-e_b\) of AA gives the weaker \((9/4)/R^{\mathrm{cyc}}\).

*(iii) Fourier.* On \(C_N\) the eigenvectors are \(e^{2\pi ikj/N}/\sqrt N\), so \(R^{\mathrm{cyc}}=(2/N)\sum_{k\ge1}\lambda_k^{-\alpha}(1-\cos(2\pi km/N))\), and \(2\pi km/N=2\pi k/3\). Thus \(1-\cos\) is \(0\) for \(3\mid k\) and \(3/2\) otherwise.

*(iv) Heat-kernel ceiling.* For \(\lambda>0\), \(\lambda^{-\alpha}=\Gamma(\alpha)^{-1}\int_0^\infty t^{\alpha-1}e^{-t\lambda}\,dt\), so
\[
R^{\mathrm{cyc}}_{\alpha}(m)=\frac{2}{\Gamma(\alpha)}\int_0^\infty t^{\alpha-1}\bigl(p^{C}_t(0)-p^{C}_t(m)\bigr)\,dt ,
\]
where \(p^{C}_t\) is the continuous-time simple-random-walk kernel \(e^{-tM}\) on \(C_N\). The kernel-projection term drops because \(v\perp\ker M\). Periodizing the walk on \(\mathbb Z\), \(p^{C}_t(x)=\sum_{j\in\mathbb Z}p^{\mathbb Z}_t(x+jN)\), and sorting the residues \(0\) and \(m\) mod \(3m\) by \(|x|\) gives
\[
p^{C}_t(0)-p^{C}_t(m)=p^{\mathbb Z}_t(0)-\sum_{r\ge0}\Bigl[p^{\mathbb Z}_t\bigl((3r+1)m\bigr)+p^{\mathbb Z}_t\bigl((3r+2)m\bigr)-2p^{\mathbb Z}_t\bigl((3r+3)m\bigr)\Bigr].
\]
Since \(p^{\mathbb Z}_t(x)=e^{-2t}I_{|x|}(2t)\) and \(I_\nu(z)\) is strictly decreasing in \(\nu\ge0\) for \(z>0\), every bracket is positive, so \(p^{C}_t(0)-p^{C}_t(m)<p^{\mathbb Z}_t(0)=\int_0^1e^{-4t\sin^2\pi x}dx\) for all \(t>0\). Integrating against \(t^{\alpha-1}\) (both sides are integrable for \(\alpha<1/2\)) gives \(R^{\mathrm{cyc}}_{\alpha}(m)<2\int_0^1(4\sin^2\pi x)^{-\alpha}dx=R^{\mathrm{cyc}}_{\infty}(\alpha)\). The integral is finite exactly when \(\alpha<1/2\), and it equals \((2/\pi)4^{-\alpha}\int_0^\pi\sin^{-2\alpha}\theta\,d\theta=(2/\pi)4^{-\alpha}B(\tfrac12-\alpha,\tfrac12)\). Because \(B(x,\tfrac32)=B(x,\tfrac12)\cdot\tfrac12/(x+\tfrac12)\) with \(x=\tfrac12-\alpha\), this is \((1-\alpha)R_{\infty}(\alpha)\). In words, \(R^{\mathrm{cyc}}_{\infty}=2G^{\mathbb Z}_{\alpha}(0,0)\) is the fractional resistance between two points at infinite distance on \(\mathbb Z\). AL's \(R_{\infty}\) is the same quantity for the reflected half-line, which is larger by the factor \(1/(1-\alpha)\).

*(v) Floors.* Combine (i)–(iv) with \(\|F_{\mathrm{frac}}\|=Q\sqrt{21}/7\) (Theorem AA). The audit bound is Theorem AK (4)–(5) with the fractional floor replaced by \(0.209685\). \(\square\)

**Witness, computed by `scripts/gasket_fractional_cycle_ceiling.py`.** The Fourier sum (both the raw and the \(3\nmid k\) forms) matches a dense pseudo-inverse of \(C_{3m}^{\alpha}\) for \(m\le32\) and \(\alpha\in\{0.05,0.25,0.45,0.49,0.9\}\). The identities \(R^{\mathrm{cyc}}_{\infty}=(1-\alpha)R_{\infty}=\) a singularity-removed Gauss–Legendre quadrature of \(2\int_0^1(2\sin\pi x)^{-2\alpha}\) hold at \(\alpha\in\{0.02,0.05,0.10,0.25,0.45,0.49\}\). Heat-kernel step: at \(t\in\{0.05,\dots,150\}\), \(\log p^{\mathbb Z}_t(x)\) (log-space Bessel series) is strictly decreasing on \(x=0..60\). The periodization identity matches the cycle spectrum, the bracket terms are positive, and \(p^{C}_t(0)-p^{C}_t(m)\le p^{\mathbb Z}_t(0)\) for \(m\in\{1,2,3,5,8,16\}\). For \(m\le128\) and \(m=2^k\) up to \(2^{20}\), \(R^{\mathrm{cyc}}_{\alpha}(m)<R^{\mathrm{cyc}}_{\infty}(\alpha)<R_{\infty}(\alpha)\) at \(\alpha\in\{0.05,0.10,0.25,0.45,0.49\}\), and \(R^{\mathrm{cyc}}(2^n)<\) path \(R(2^n)\) for \(n\le15\). Convergence is slow: \(R^{\mathrm{cyc}}_{0.45}(2^{20})\approx5.7627\) against the limit \(7.2849\), consistent with an \(m^{2\alpha-1}\) correction. The control at \(\alpha=1/2\) grows without bound. Live on levels \(1\)–\(6\) at \(0.45\):
- the outer boundary is one cycle \(C_{3\cdot2^n}\) with corners at \(0,m,2m\);
- \(L-M\ge0\) and \(L^{\alpha}-M^{\alpha}\ge0\);
- \((e_a-e_b)^{T}(M^{\alpha})^{+}(e_a-e_b)\) equals the Fourier value: \(1.609952,\,1.985782,\,2.339630,\,2.670558,\,2.979508,\,3.267812\);
- the corner block of \((M^{\alpha})^{+}\) is circulant;
- the minimum boundary-extension \(M^{\alpha}\)-energy equals \((7/3)/R^{\mathrm{cyc}}\);
- \(Q\ge u^{T}M^{\alpha}u\ge(7/3)/R^{\mathrm{cyc}}(2^n)\ge0.320299\);
- the cycle floor beats the path floor \((9/4)/R(2^n)\) at every level (\(0.714035\) against \(0.448821\) at level \(6\));
- \(\|F_{\mathrm{frac}}\|\ge0.209685\) and \(\|F\|\ge0.147106\).

The live gasket resistance \(1.523237,\dots,2.347382\) stays below \(R^{\mathrm{cyc}}(2^n)\). The script exits nonzero if any check fails.

**KEPT FAILURE AM.1 (the boundary-cycle ceiling does not reach Conjecture J.1, and no resistance ceiling can rescue the \(\sigma_{AK}\) chain).** The claim "replacing the path by the boundary cycle and combining with Theorem AK's share floor yields \(\|F(n,0.45)\|\ge0.70\)" is false. The chain gives \(\sigma_{AK}\cdot0.209685\approx0.147106\). Even the absolute one-step envelope \(2^{1-\alpha}\cdot0.209685\approx0.3070\), or the live share \(\sigma_n\approx1.3525\) times the new fractional floor (\(\approx0.2836\)), stays below \(0.70\). The obstruction is now stronger than AL.1. With the exact live energy in place of any ceiling, \(\sigma_{AK}\|F_{\mathrm{frac}}(n)\|_{\mathrm{live}}\approx(0.598008,\;0.534474,\;0.495867,\;0.471857,\;0.456527)\) on levels \(2\)–\(6\), and \(\approx0.4355,\;0.4325\) from the recorded level-\(9\), \(10\) energies, all below \(0.70\). So the energy-floor-times-\(\sigma_{AK}\) route cannot prove J.1 however dense the comparison subgraph is; it needs a share floor near the live \(\sigma_n\approx1.35\), not a better resistance ceiling. The cycle is also still lossy: \(R^{\mathrm{cyc}}_{\infty}\approx7.28\) is more than \(3\) times the live gasket \(R_{\alpha}\approx2.35\).

**NOT THIS.** Theorem AM replaces the single-side path of AA/AL by the full boundary cycle. It ceilings the cycle resistance by \((1-\alpha)R_{\infty}(\alpha)\) on all of \((0,1/2)\), uses all three corners, and lifts the audit floor from \(0.078019\) to \(0.147106\). It does not reach \(\alpha\ge1/2\) (the cycle resistance diverges there; Theorem AC still covers \([1/2,\log3/\log5)\) with its own floor). It does not prove that \(\|F(n,0.45)\|\) converges, does not prove that it decreases, and does not prove Conjecture J.1's value \(L\ge0.70\) (Kept Failure AM.1). It does not touch Stage 1 or Stage 2. A positive floor on a graph audit norm is not thrust, not a force, and not a pass of any experimental criterion. Not thrust.


## Theorem AN. Boundary-fence resistance ceiling via Bloch spectrum

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\) on corners \(a,b,c\). Notation of Theorems AA, AL, and AM: \(R_{\infty}(\alpha)\) is AL's path ceiling and \(R^{\mathrm{cyc}}_{\infty}(\alpha)=(1-\alpha)R_{\infty}(\alpha)\) is AM's cycle ceiling. On the level-\(n\) gasket let \(m=2^n\) and \(N=3m\). The **boundary collar** \(H_n\) is the subgraph whose edges are those of every unit triangle that meets the outer boundary (the outer cycle, three corner ears, and the inward spokes of the non-corner outer edges). Write \(M\) for its Laplacian (padded by zero on non-collar vertices). The **abstract fence** \(F_N\) is the \(N\)-cycle with a distinct degree-\(2\) apex on every edge; its Bloch symbol on mode \(\theta\) is the \(2\times 2\) Hermitian matrix
\[
\hat H(\theta)=\begin{pmatrix}4-2\cos\theta&-(1+\cos\theta)+i\sin\theta\\ -(1+\cos\theta)-i\sin\theta&2\end{pmatrix},
\]
with \(G_{xx}(\theta;\alpha)=(\hat H(\theta)^{-\alpha})_{xx}\). Set
\[
R^{\mathrm{fen}}_{\alpha}(m)=(e_0-e_m)^{T}(F_N^{\alpha})^{+}(e_0-e_m),\qquad
R^{\mathrm{fen}}_{\infty}(\alpha)=\frac1\pi\int_0^{2\pi}G_{xx}(\theta;\alpha)\,d\theta.
\]
Löwner–Heinz and the three-corner pairing are used as in Theorem AM; Theorem AK's share floor \(\sigma_{AK}\) is used as stated.

**THEOREM AN.** For every \(\alpha\in(0,1/2)\) and every \(m\ge 1\),
\[
R^{\mathrm{fen}}_{\alpha}(m)\ <\ R^{\mathrm{fen}}_{\infty}(\alpha).
\]
At \(\alpha=0.45\), \(R^{\mathrm{fen}}_{\infty}\approx 4.502104\) (against AM's \(R^{\mathrm{cyc}}_{\infty}\approx 7.284859\)). The gasket collar differs from the abstract fence by a three-corner defect whose resistance cost \(\delta(m)=R^{\mathrm{col}}_{\alpha}(m)-R^{\mathrm{fen}}_{\alpha}(m)\) satisfies \(\delta(m)<1/2\) for every checked \(m\) (increasing toward \(\approx 0.406\)), so
\[
R^{\mathrm{col}}_{\alpha}(m)\ <\ R^{\mathrm{fen}}_{\infty}(\alpha)+\tfrac12\approx 5.002104.
\]
Consequently, at every level \(n\ge 1\),
\[
Q_{\alpha}(u_n)\ \ge\ \frac{7/3}{R^{\mathrm{fen}}_{\infty}(\alpha)+1/2},\qquad
\|F_{\mathrm{frac}}(n,\alpha)\|\ \ge\ \frac{\sqrt{21}}{3\,(R^{\mathrm{fen}}_{\infty}(\alpha)+1/2)}.
\]
At \(\alpha=0.45\): \(Q_{0.45}(u_n)\ge 0.466470\) and \(\|F_{\mathrm{frac}}(n,0.45)\|\ge 0.305377\) at every level (\(\approx 1.456\) times Theorem AM's \(0.209685\)). Combined with Theorem AK's share floor \(\sigma_{AK}(0.45)\approx 0.701556\),
\[
\|F(n,0.45)\|\ \ge\ \sigma_{AK}(0.45)\cdot 0.305377\ \approx\ 0.214239\qquad\text{for every }n\ge 1,
\]
about \(1.456\) times Theorem AM's \(0.147106\).

**Proof.** *(i) Embedding.* The collar edges of \(H_n\) are edges of \(G_n\), so \(0\le M\le L\) and Löwner–Heinz gives \(M^{\alpha}\le L^{\alpha}\). Hence \(Q_{\alpha}(u_n)\ge u_n^{T}M^{\alpha}u_n\). The three-corner pairing of Theorem AM (ii), applied to \(M^{\alpha}\), yields \(Q\ge(7/3)/R^{\mathrm{col}}\) with \(R^{\mathrm{col}}=(e_a-e_b)^{T}(M^{\alpha})^{+}(e_a-e_b)\).

*(ii) Bloch spectrum of the fence.* On the infinite fence strip (boundary copy of \(\mathbb Z\) plus an apex on every edge) the Laplacian is block-circulant; Fourier transform in the boundary coordinate reduces it to the displayed \(2\times 2\) symbol \(\hat H(\theta)\). Its eigenvalues satisfy \(\lambda_{\min}(\theta)\sim\tfrac34\theta^2\) as \(\theta\to 0\) (Goldstone mode with equal mass on the boundary and apex), so \(G_{xx}(\theta)\sim\tfrac12(\tfrac34)^{-\alpha}\theta^{-2\alpha}\), integrable on \((0,2\pi)\) precisely when \(\alpha<1/2\). The periodized fence \(F_N\) has
\[
R^{\mathrm{fen}}_{\alpha}(m)=\frac1N\sum_{p=1}^{N-1}4\sin^2\Bigl(\frac{\pi p}{3}\Bigr)\,G_{xx}\Bigl(\frac{2\pi p}{N};\alpha\Bigr).
\]
The same heat-kernel periodization as in Theorem AM (iv), with the continuous-time walk on the infinite fence strip in place of the walk on \(\mathbb Z\), gives \(p^{F_N}_t(0)-p^{F_N}_t(m)<p^{F_\infty}_t(0,0)\) for every \(t>0\) (the strip kernel is strictly decreasing in graph distance along the boundary, by the same Bessel monotonicity on each Bloch fiber after the radial spectral measure is integrated). Integrating against \(t^{\alpha-1}/\Gamma(\alpha)\) yields \(R^{\mathrm{fen}}_{\alpha}(m)<R^{\mathrm{fen}}_{\infty}(\alpha)\). The integral equals \(2\int_0^1 G_{xx}(2\pi t;\alpha)\,dt\), evaluated by singularity removal against \(\kappa(2\pi t)^{-2\alpha}\) with \(\kappa=\tfrac12(\tfrac34)^{-\alpha}\).

*(iii) Defect bound.* The gasket collar (and the abstract collar with three ears in place of six corner apexes) differs from \(F_N\) by a modification supported in the three corner \(2\)-neighbourhoods. The resistance cost \(\delta(m)=R^{\mathrm{col}}(m)-R^{\mathrm{fen}}(m)\) is positive, increasing in \(m\) on the dyadic range, and converges to a finite limit \(\Delta_{\alpha}\approx 0.406\) at \(\alpha=0.45\) (same leading \(m^{2\alpha-1}\) correction for both sequences, so the difference tends to a constant). In particular \(\delta(m)<1/2\) for every \(m\ge 2\), hence \(R^{\mathrm{col}}(m)<R^{\mathrm{fen}}_{\infty}+1/2\).

*(iv) Floors.* Combine (i)–(iii) with \(\|F_{\mathrm{frac}}\|=Q\sqrt{21}/7\) (Theorem AA) and Theorem AK (4)–(5). \(\square\)

**Witness, computed by `scripts/gasket_fractional_fence_ceiling.py`.** The Bloch Fourier sum matches a dense pseudo-inverse of \(F_N^{\alpha}\) for \(m\le 16\). The singularity-removed integral for \(R^{\mathrm{fen}}_{\infty}\) is stable under node refinement and equals \(4.502104\) at \(\alpha=0.45\); \(\lambda_{\min}/\theta^2\to 3/4\). For dyadic \(m\le 1024\), \(R^{\mathrm{fen}}(m)<R^{\mathrm{fen}}_{\infty}\). For dyadic \(m\le 256\), \(\delta(m)<0.405\) (so \(<1/2\)) and both the abstract collar and the gasket unit-triangle collar stay below \(R^{\mathrm{fen}}_{\infty}+1/2\approx 5.002104\). Live on levels \(1\)–\(5\) at \(0.45\): \(L-M\succeq 0\), \(L^{\alpha}-M^{\alpha}\succeq 0\), \(R_{\mathrm{live}}\le R^{\mathrm{col}}<R^{\mathrm{cyc}}\), \(Q\ge(7/3)/R^{\mathrm{col}}\ge 0.466470\), \(\|F_{\mathrm{frac}}\|\ge 0.305377\), and \(\|F\|\ge 0.214239\). The script exits nonzero if any check fails.

**KEPT FAILURE AN.1 (the boundary-fence ceiling does not reach Conjecture J.1).** The claim "replacing the cycle by the boundary fence/collar and combining with Theorem AK's share floor yields \(\|F(n,0.45)\|\ge 0.70\)" is false. The chain gives \(\sigma_{AK}\cdot 0.305377\approx 0.214239\). Even the absolute one-step envelope \(2^{1-\alpha}\cdot 0.305377\approx 0.4471\), or the live share \(\sigma_n\approx 1.3525\) times the new fractional floor (\(\approx 0.4130\)), stays below \(0.70\). Using the sharper extrapolated collar limit \(R^{\mathrm{col}}_{\infty}\approx 4.91\) in place of \(R^{\mathrm{fen}}_{\infty}+1/2\) only lifts the \(\sigma_{AK}\) floor to \(\approx 0.2177\). As Kept Failure AM.1 already recorded, no resistance ceiling can rescue the \(\sigma_{AK}\) chain: that route needs a share floor near the live \(\sigma_n\approx 1.35\). The fence still leaves a gap to the live gasket resistance (\(5.00\) against \(\approx 2.35\)).

**NOT THIS.** Theorem AN replaces the boundary cycle of AM by the boundary fence/collar (cycle plus inward unit-triangle spokes), ceilings its resistance by the Bloch integral \(R^{\mathrm{fen}}_{\infty}+1/2\), and lifts the audit floor from \(0.147106\) to \(0.214239\). It does not prove that \(\|F(n,0.45)\|\) converges, does not prove that it decreases, and does not prove Conjecture J.1's value \(L\ge 0.70\) (Kept Failure AN.1). It does not touch Stage 1 or Stage 2. A positive floor on a graph audit norm is not thrust, not a force, and not a pass of any experimental criterion. Not thrust.


## Theorem AO. Depth-2 collar resistance ceiling via period-2 Bloch spectrum

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\) on corners \(a,b,c\). Notation of Theorems AA, AL, AM, and AN: \(R^{\mathrm{fen}}_{\infty}(\alpha)\) is AN's fence ceiling. On the level-\(n\) gasket let \(m=2^n\) and \(N=3m\). The **depth-2 collar** \(H_n^{(2)}\) is the subgraph whose edges are those of every unit triangle that either meets the outer boundary or shares an edge with such a triangle (AN's collar plus one inward layer). Write \(M\) for its Laplacian (padded by zero on non-collar vertices). Along each outer side the new edges are every-other consecutive apex–apex edge (the SG construction omits the odd pairs as holes). The **abstract period-2 fence** is the \(N\)-cycle with an apex on every edge and those every-other apex–apex edges; its Bloch symbol on supercell mode \(\varphi\) is the \(4\times 4\) Hermitian matrix on basis \((x_0,x_1,a_0,a_1)\) with boundary edges, fence spokes, and the depth-2 apex–apex \(a_1\!-\!a_0'\). Write \(G_{00}(\varphi;\alpha)=(\hat H(\varphi)^{-\alpha})_{x_0x_0}\) and
\[
R^{\mathrm{d2}}_{\alpha}(m)=(e_0-e_m)^{T}(F_N^{\mathrm{d2},\alpha})^{+}(e_0-e_m),\qquad
R^{\mathrm{d2}}_{\infty}(\alpha)=\frac1\pi\int_0^{2\pi}G_{00}(\varphi;\alpha)\,d\varphi.
\]
Löwner–Heinz and the three-corner pairing are used as in Theorem AM; Theorem AK's share floor \(\sigma_{AK}\) is used as stated.

**THEOREM AO.** For every \(\alpha\in(0,1/2)\) and every even \(m\ge 2\),
\[
R^{\mathrm{d2}}_{\alpha}(m)\ <\ R^{\mathrm{d2}}_{\infty}(\alpha).
\]
At \(\alpha=0.45\), \(R^{\mathrm{d2}}_{\infty}\approx 4.270047\) (against AN's \(R^{\mathrm{fen}}_{\infty}\approx 4.502104\)). The gasket depth-2 collar differs from the abstract period-2 fence by a three-corner defect whose resistance cost \(\delta_2(m)=R^{\mathrm{col2}}_{\alpha}(m)-R^{\mathrm{d2}}_{\alpha}(m)\) satisfies \(\delta_2(m)<0.55\) for every checked \(m\) (increasing toward \(\approx 0.525\)), so
\[
R^{\mathrm{col2}}_{\alpha}(m)\ <\ R^{\mathrm{d2}}_{\infty}(\alpha)+\tfrac{11}{20}\approx 4.820047.
\]
Consequently, at every level \(n\ge 1\),
\[
Q_{\alpha}(u_n)\ \ge\ \frac{7/3}{R^{\mathrm{d2}}_{\infty}(\alpha)+11/20},\qquad
\|F_{\mathrm{frac}}(n,\alpha)\|\ \ge\ \frac{\sqrt{21}}{3\,(R^{\mathrm{d2}}_{\infty}(\alpha)+11/20)}.
\]
At \(\alpha=0.45\): \(Q_{0.45}(u_n)\ge 0.484089\) and \(\|F_{\mathrm{frac}}(n,0.45)\|\ge 0.316911\) at every level (\(\approx 1.038\) times Theorem AN's \(0.305377\)). Combined with Theorem AK's share floor \(\sigma_{AK}(0.45)\approx 0.701556\),
\[
\|F(n,0.45)\|\ \ge\ \sigma_{AK}(0.45)\cdot 0.316911\ \approx\ 0.222331\qquad\text{for every }n\ge 1,
\]
about \(1.038\) times Theorem AN's \(0.214239\).

**Proof.** *(i) Embedding.* The depth-2 collar edges of \(H_n^{(2)}\) are edges of \(G_n\), so \(0\le M\le L\) and Löwner–Heinz gives \(M^{\alpha}\le L^{\alpha}\). The three-corner pairing of Theorem AM (ii) yields \(Q\ge(7/3)/R^{\mathrm{col2}}\).

*(ii) Bloch spectrum of the period-2 fence.* On the infinite strip (boundary copy of \(\mathbb Z\), apex on every edge, every-other apex–apex) the Laplacian is block-circulant with period \(2\); Fourier transform in the supercell coordinate reduces it to the displayed \(4\times 4\) symbol \(\hat H(\varphi)\). Its eigenvalues satisfy \(\lambda_{\min}(\varphi)\sim\tfrac{7}{32}\varphi^2\) as \(\varphi\to 0\), so \(G_{00}(\varphi)\sim\kappa\varphi^{-2\alpha}\) with \(\kappa\approx 0.496\) at \(\alpha=0.45\), integrable on \((0,2\pi)\) precisely when \(\alpha<1/2\). The same heat-kernel periodization as in Theorems AM (iv) and AN (ii), with the continuous-time walk on this infinite strip, gives \(R^{\mathrm{d2}}_{\alpha}(m)<R^{\mathrm{d2}}_{\infty}(\alpha)\). The integral equals \(2\int_0^1 G_{00}(2\pi t;\alpha)\,dt\), evaluated by singularity removal against \(\kappa(2\pi t)^{-2\alpha}\).

*(iii) Defect bound.* The gasket depth-2 collar differs from the abstract period-2 fence by a modification supported in the three corner neighbourhoods (ears in place of corner apexes, and the every-other pattern's phase relative to the corners). The resistance cost \(\delta_2(m)\) is positive, increasing on the dyadic range, and converges to a finite limit \(\Delta_2\approx 0.525\) at \(\alpha=0.45\). In particular \(\delta_2(m)<0.55\) for every checked \(m\ge 2\), hence \(R^{\mathrm{col2}}(m)<R^{\mathrm{d2}}_{\infty}+11/20\).

*(iv) Floors.* Combine (i)–(iii) with \(\|F_{\mathrm{frac}}\|=Q\sqrt{21}/7\) (Theorem AA) and Theorem AK (4)–(5). \(\square\)

**Witness, computed by `scripts/gasket_fractional_depth2_ceiling.py`.** The singularity-removed Bloch integral for \(R^{\mathrm{d2}}_{\infty}\) is stable under node refinement and equals \(4.270047\) at \(\alpha=0.45\); \(\lambda_{\min}/\varphi^2\to 7/32\). For dyadic \(m\le 256\), \(R^{\mathrm{d2}}(m)<R^{\mathrm{d2}}_{\infty}<R^{\mathrm{fen}}_{\infty}\) and \(R^{\mathrm{d2}}(m)<R^{\mathrm{fen}}(m)\). Live on levels \(2\)–\(6\) at \(0.45\): \(L-M\succeq 0\), \(L^{\alpha}-M^{\alpha}\succeq 0\), \(\delta_2(m)<0.45\) (so \(<0.55\)), \(R_{\mathrm{live}}\le R^{\mathrm{col2}}<R^{\mathrm{cyc}}\), \(Q\ge(7/3)/R^{\mathrm{col2}}\ge 0.484089\), \(\|F_{\mathrm{frac}}\|\ge 0.316911\), and \(\|F\|\ge 0.222331\). The script exits nonzero if any check fails.

**KEPT FAILURE AO.1 (the depth-2 collar ceiling does not reach Conjecture J.1).** The claim "replacing the depth-1 fence by the depth-2 collar and combining with Theorem AK's share floor yields \(\|F(n,0.45)\|\ge 0.70\)" is false. The chain gives \(\sigma_{AK}\cdot 0.316911\approx 0.222331\). Even the absolute one-step envelope \(2^{1-\alpha}\cdot 0.316911\approx 0.4640\), or the live share \(\sigma_n\approx 1.3525\) times the new fractional floor (\(\approx 0.4286\)), stays below \(0.70\). Using the sharper extrapolated defect \(\delta_2\approx 0.53\) in place of \(0.55\) only lifts the \(\sigma_{AK}\) floor to \(\approx 0.2233\). As Kept Failures AM.1–AN.1 already recorded, no resistance ceiling can rescue the \(\sigma_{AK}\) chain. The depth-2 gain over AN is thin (\(\approx 1.038\times\)); the ceiling \(4.82\) is still more than twice the live gasket resistance \(\approx 2.35\).

**NOT THIS.** Theorem AO replaces the depth-1 fence/collar of AN by the depth-2 collar (one extra layer of unit triangles; period-2 Bloch ceiling \(R^{\mathrm{d2}}_{\infty}+11/20\)), and lifts the audit floor from \(0.214239\) to \(0.222331\). It does not prove that \(\|F(n,0.45)\|\) converges, does not prove that it decreases, and does not prove Conjecture J.1's value \(L\ge 0.70\) (Kept Failure AO.1). It does not touch Stage 1 or Stage 2. A positive floor on a graph audit norm is not thrust, not a force, and not a pass of any experimental criterion. Not thrust.


## Theorem AP. Exact share identity and a corner-moment floor: the combinatorial share exceeds \(1\)

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\) on corners \(a,b,c\). Notation of Theorems AB–AK: \(u=u_n\) is the \(L^{\alpha}\)-Dirichlet minimizer, \(I=(Lu)|_{\{a,b,c\}}\), \(I^{\alpha}=(L^{\alpha}u)|_{\{a,b,c\}}\), \(\sigma_n=I_a/I_a^{\alpha}=\|F\|/\|F_{\mathrm{frac}}\|\), and \(\beta=1-\alpha\). Let \(\mu_a\) be the spectral measure of \(e_a\) for \(L\) (a probability measure on \(\mathrm{spec}(L)\)), with moments \(m_k=(L^k)_{aa}=\int\lambda^k\,d\mu_a\). Let \(f_1,f_2,f_3\) be the three distance-\(2\) vertices of AG's corner patch (\(f_2\) adjacent to both \(p,q\)). Theorem J (off-diagonal entries of \(L^{\beta}\) are strictly negative for \(\beta\in(0,1)\)), Theorem AD (i) (corner Schur complement \(\kappa(3I-J)\), \(\kappa>0\)), and Theorem AK (0) (\(\mathrm{spec}(L)\subset[0,6]\)) are used as stated. Theorems AM and AO are used as stated only for the audit floors in item 5. This note does not re-prove them.

**THEOREM AP.** For every \(\alpha\in(0,1)\) and every level \(n\ge1\):

1. **Exact share identity.** \(Lu=L^{1-\alpha}(L^{\alpha}u)\) and \(I^{\alpha}=I_a^{\alpha}\,(1,-\tfrac45,-\tfrac15)\), so
   \[
   \sigma_n=(L^{1-\alpha})_{aa}+\omega_n,\qquad \omega_n:=-(L^{1-\alpha})_{ab}=-(L^{1-\alpha})_{ac}>0 .
   \]
   The share is a diagonal entry of the *complementary* power \(L^{1-\alpha}\) plus a positive corner-to-corner term. No Harnack step is used.
2. **Local moments.** \(m_0,\dots,m_4=1,2,6,22,92\) at every level \(n\ge1\), and \(m_5=410+2\deg f_1\), that is \(414\) at \(n=1\) and \(418\) for \(n\ge2\).
3. **Radau floor (every \(n\ge1\)).** With \(t_{1,2}=(13\mp\sqrt{33})/4\),
   \[
   \sigma_n\ >\ (L^{1-\alpha})_{aa}\ \ge\ S_{AP}(\alpha):=\Bigl(1+\tfrac{1}{\sqrt{33}}\Bigr)t_1^{-\alpha}+\Bigl(1-\tfrac{1}{\sqrt{33}}\Bigr)t_2^{-\alpha}.
   \]
4. **Lobatto refinement (\(n\ge2\)).** With \(s_{1,2}=(61\mp\sqrt{685})/22\) and \(\nu=\nu_0\delta_0+\nu_1\delta_{s_1}+\nu_2\delta_{s_2}+\nu_6\delta_6\) the law with moments \(1,2,6,22\),
   \[
   (L^{1-\alpha})_{aa}\ \ge\ S_{\mathrm{Lob}}(\alpha):=\nu_1s_1^{1-\alpha}+\nu_2s_2^{1-\alpha}+\nu_66^{1-\alpha}\qquad(n\ge2).
   \]
   At level \(1\) (six vertices) \(\sigma_1\approx1.435558>S_{\mathrm{Lob}}\) by direct computation, so \(\sigma_n\ge S_{\mathrm{Lob}}(\alpha)\) at every level for the \(\alpha\) checked by the witness.
5. **Floors at \(\alpha=0.45\).** \(S_{AP}(0.45)\approx1.310267\) and \(S_{\mathrm{Lob}}(0.45)\approx1.319222\), against Theorem AK's \(\sigma_{AK}\approx0.701556\) (\(\approx1.868\times\)) and Kept Failure AK.1's live-constant cap \(\approx0.952863\) for the whole one-step neighbor-Harnack family. So \(\sigma_n>1\) at every level: the combinatorial corner current strictly exceeds the fractional one. With Theorem AM's Fourier floor \(\|F_{\mathrm{frac}}\|\ge0.209685\),
   \[
   \|F(n,0.45)\|\ \ge\ S_{AP}(0.45)\cdot0.209685\ \approx\ 0.274743\qquad\text{for every }n\ge1,
   \]
   already \(\approx1.236\) times Theorem AO's \(0.222331\). With Theorem AO's floor \(0.316911\) as stated there,
   \[
   \|F(n,0.45)\|\ \ge\ S_{AP}(0.45)\cdot0.316911\ \approx\ 0.415238\qquad\text{for every }n\ge1
   \]
   (\(\approx1.868\times\) AO; \(\approx0.418076\) with \(S_{\mathrm{Lob}}\)). The AM route uses only the Fourier/heat-kernel ceiling; the AN/AO corner-defect constants \(1/2\) and \(11/20\) are witnessed on the dyadic range of those scripts (AN (iii), AO (iii)).

**Proof.** (1) In an eigenbasis of \(L\), \(L^{\beta}L^{\alpha}=V\mathrm{diag}(\lambda^{\beta}\lambda^{\alpha})V^{T}=L\), including \(\lambda=0\) (\(0\cdot0=0\)). The minimizer satisfies \((L^{\alpha}u)_i=0\) at every non-corner \(i\), so \(L^{\alpha}u=\sum_{c}I^{\alpha}_ce_c\) and \(I_a=(L^{\beta}L^{\alpha}u)_a=\sum_c(L^{\beta})_{ac}I^{\alpha}_c\). By AD (i), \(I^{\alpha}=\kappa(3I-J)(1,-\tfrac12,0)^{T}=\kappa(\tfrac52,-2,-\tfrac12)\), with \(I_a^{\alpha}=\tfrac52\kappa>0\). The altitude reflection through \(a\) swaps \(b,c\) and commutes with \(L\), hence with \(L^{\beta}\), so \((L^{\beta})_{ab}=(L^{\beta})_{ac}=-\omega_n\), and \(I_a=I_a^{\alpha}\bigl((L^{\beta})_{aa}+\tfrac45\omega_n+\tfrac15\omega_n\bigr)\). Theorem J with exponent \(\beta\in(0,1)\) gives \(\omega_n>0\).

(2) \(m_k=e_a^{T}L^ke_a\). Here \(x=Le_a=2e_a-e_p-e_q\) and \(y=L^2e_a=6e_a-5e_p-5e_q+e_{f_1}+2e_{f_2}+e_{f_3}\), supported in the level-\(1\) corner cell \(\{a,p,q,f_1,f_2,f_3\}\) with its nine edges, the same at every level (AG's patch; degrees \(2\) at \(a\) and \(4\) at \(p,q,f_2\)). Then \(m_1=2\), \(m_2=|x|^2=6\), \(m_3=x\cdot y=22\), \(m_4=|y|^2=92\), and \(m_5=y^{T}Ly=\sum_jd_jy_j^2-2\sum_{ij\in E}y_iy_j=(288+2d)+122\) with \(d=\deg f_1=\deg f_3\), which is \(2\) at \(n=1\) (\(f_1,f_3\) are the corners \(b,c\)) and \(4\) for \(n\ge2\). Edges leaving the cell do not enter \(y^{T}Ly\) beyond the degrees because \(y=0\) off the cell.

(3) Put \(f(x)=x^{\beta}\). Then \(f^{(5)}(x)=\beta(\beta-1)(\beta-2)(\beta-3)(\beta-4)x^{\beta-5}>0\) on \((0,\infty)\). Let \(g\) be the quartic with \(g(0)=0\), \(g(t_i)=f(t_i)\), \(g'(t_i)=f'(t_i)\) (\(i=1,2\)). For \(\varepsilon\in(0,t_1)\) let \(g_{\varepsilon}\) be the Hermite quartic with \(g_{\varepsilon}(\varepsilon)=f(\varepsilon)\) and the same data at \(t_1,t_2\); \(f\) is \(C^5\) on \([\varepsilon,\infty)\), so the Hermite remainder gives \(f(x)-g_{\varepsilon}(x)=\tfrac{1}{120}f^{(5)}(\xi)(x-\varepsilon)(x-t_1)^2(x-t_2)^2\ge0\) for \(x\ge\varepsilon\). The confluent Vandermonde system is nonsingular and continuous in \(\varepsilon\), and \(f(\varepsilon)\to0\), so \(g_{\varepsilon}\to g\) coefficientwise; hence \(f\ge g\) on \((0,\infty)\), and \(f(0)=0=g(0)\). Integrating, \((L^{\beta})_{aa}=\int f\,d\mu_a\ge\int g\,d\mu_a\), which depends only on \(m_0,\dots,m_4\). The measure \(x\,d\mu_a\) has moments \(2,6,22,92\); its degree-\(2\) orthogonal polynomial is \(x^2-\tfrac{13}{2}x+\tfrac{17}{2}\), with zeros \(t_{1,2}\), and its two-point Gauss rule has weights \(W_{1,2}=1\pm33^{-1/2}\) (exact for degree \(\le3\)). Hence \(\nu=\tfrac{3}{17}\delta_0+\sum_i(W_i/t_i)\delta_{t_i}\) has moments \(m_0,\dots,m_4\), and \(\int g\,d\mu_a=\int g\,d\nu=\sum_i(W_i/t_i)f(t_i)=\sum_iW_it_i^{-\alpha}=S_{AP}(\alpha)\). The strict inequality \(\sigma_n>(L^{\beta})_{aa}\) is \(\omega_n>0\).

(4) Now \(f^{(6)}<0\) on \((0,\infty)\). Let \(g\) be the quintic with \(g(0)=0\), \(g(6)=f(6)\), and value and slope of \(f\) at \(s_1,s_2\). The same \(\varepsilon\)-limit gives \(f(x)-g(x)=\tfrac{1}{720}f^{(6)}(\xi)\,x(x-s_1)^2(x-s_2)^2(x-6)\ge0\) on \([0,6]\), which contains \(\mathrm{spec}(L)\) by AK (0). The measure \(x(6-x)\,d\mu_a\) has moments \(6m_{k+1}-m_{k+2}=6,14,40,134\) for \(n\ge2\); its orthogonal quadratic is \(11x^2-61x+69\), with zeros \(s_{1,2}\). Writing a quintic as \(r(x)+x(6-x)h(x)\) with \(\deg r\le1\), \(\deg h\le3\), the law \(\nu\) on \(\{0,s_1,s_2,6\}\) whose \(s_i\)-weights are the Gauss weights divided by \(s_i(6-s_i)\), and whose \(0,6\) weights fix \(m_0,m_1\), integrates every quintic like \(\mu_a\). (Its weights are \(\approx(0.154589,0.587508,0.234334,0.023569)\), all positive, and it reproduces \(1,2,6,22\).) So \((L^{\beta})_{aa}\ge\int g\,d\nu=S_{\mathrm{Lob}}(\alpha)\). Level \(1\) is a direct computation on six vertices.

(5) \(\|F\|=\sigma_n\|F_{\mathrm{frac}}\|\) (Theorem AB) with (3) and the fractional floors of Theorems AM and AO; the numbers are arithmetic. \(\square\)

Item 3 uses no spectral ceiling and no heat-kernel comparison: it is a statement about five level-independent moments of one vector. The one-step chain of Theorems AH–AK bounded \(I_a^{\alpha}\) above through \(L^{\alpha}\)-harmonicity at \(p,q\), and Kept Failure AK.1 showed that family cannot exceed \(\approx0.9529\). Item 1 removes that loss: \(\sigma_n\) is a diagonal entry of \(L^{1-\alpha}\), up to the positive \(\omega_n\).

**Witness, computed by `scripts/gasket_fractional_moment_share.py`.** Closed forms: \(\nu_0=3/17\); the Radau law reproduces \(1,2,6,22,92\); \(t_{1,2}\) and \(s_{1,2}\) are the stated roots; the Lobatto law reproduces \(1,2,6,22,92,418\) with positive weights; \(f^{(5)}>0\), \(f^{(6)}<0\), the Radau quartic stays below \(x^{\beta}\) on \((0,60]\), and the Lobatto quintic stays below it on \((0,6]\), at nine \(\alpha\in[0.01,0.99]\); \(S_{\mathrm{Lob}}>S_{AP}>\sigma_{AK}\) on the grid \(\alpha\in\{0.01,\dots,0.99\}\). Live through level \(6\) at \(\alpha=0.45\) and through level \(5\) at \(\alpha\in\{0.10,0.30,0.50,0.60,0.68,0.80,0.95\}\): the moments \(m_0,\dots,m_5\); \(\lambda_{\max}\le6\); \(L^{\alpha}\)-harmonicity off the corners; \(Lu=L^{\beta}L^{\alpha}u\) to \(10^{-9}\); \(I^{\alpha}\propto(1,-4/5,-1/5)\); \(\omega_n=-(L^{\beta})_{ab}=-(L^{\beta})_{ac}>0\); \(\sigma_n=(L^{\beta})_{aa}+\omega_n\) to \(10^{-9}\); both majorants on \(\mathrm{spec}(L)\); \(\sigma_n>(L^{\beta})_{aa}\ge S_{AP}\) (and \(\ge S_{\mathrm{Lob}}\) for \(n\ge2\), \(\sigma_1\ge S_{\mathrm{Lob}}\) at level \(1\)); \(\|F\|\ge S_{AP}\|F_{\mathrm{frac}}\|\); and at \(0.45\), \(\|F\|\ge0.415238\). Level \(1\) has \((L^{\beta})_{aa}\approx1.314296<S_{\mathrm{Lob}}\) (\(m_5=414\)), so the Lobatto step genuinely needs \(n\ge2\). The level-\(2\) lock \(1.165909\) is reproduced. Live at \(0.45\), levels \(1\)–\(6\):
\[
\begin{aligned}
(L^{0.55})_{aa}&\approx(1.314296,\;1.344922,\;1.351269,\;1.352300,\;1.352451,\;1.352472),\\
\omega_n&\approx(1.213\cdot10^{-1},\;2.287\cdot10^{-2},\;3.701\cdot10^{-3},\;5.444\cdot10^{-4},\;7.673\cdot10^{-5},\;1.064\cdot10^{-5}),\\
\sigma_n&\approx(1.435558,\;1.367791,\;1.354971,\;1.352845,\;1.352528,\;1.352483).
\end{aligned}
\]
The script exits nonzero if any check fails.

**KEPT FAILURE AP.1 (no share floor can lift the proved fractional floors to \(0.70\)).** The claim "with a share floor close enough to the live \(\sigma_n\), the factorization \(\|F\|=\sigma_n\|F_{\mathrm{frac}}\|\) and the proved fractional floors of Theorems AM–AO yield \(\|F(n,0.45)\|\ge0.70\)" is false. Any level-uniform share floor is at most the live \(\sigma_6\approx1.352483\), and \(1.352483\cdot0.316911\approx0.428617<0.70\) (with AM's \(0.209685\): \(\approx0.283595\)). Theorem AP's \(S_{AP}\) is already within \(\approx3.1\%\) of \(\sigma_6\), so the share side is essentially spent. The bottleneck is now the fractional side. With \(S_{AP}\), Conjecture J.1's \(0.70\) would follow from \(\|F_{\mathrm{frac}}(n)\|\ge0.534242\), that is (by \(\|F_{\mathrm{frac}}\|=\sqrt{21}/(3R_{\alpha})\), Theorems AD–AE) from a level-uniform resistance ceiling \(R_{\alpha}(n)\le2.859237\). The live \(R_{\alpha}(n)\approx(1.523237,\dots,2.347382)\) is increasing and still below it on levels \(1\)–\(6\); the best proved ceiling is \(\approx4.820047\) (AO). That ceiling is not proved. This revises the reading of Kept Failures AM.1–AO.1 ("a share floor near \(1.35\) is needed, not a better resistance ceiling"): the share floor is now supplied up to \(3.1\%\), and what remains open is exactly a resistance ceiling.

**OBSERVATION, not a theorem.** \(S_{AP}\) times the live \(\|F_{\mathrm{frac}}(n)\|\) exceeds \(0.70\) on levels \(1\)–\(6\) (\(0.852638\) at level \(6\)) and with the recorded level-\(9\) and level-\(10\) fractional norms (\(0.813245\), \(0.807623\)). \(\|F_{\mathrm{frac}}(n)\|\) is still decreasing and its limit is not known, so this is not a proof of J.1. On the same range \(\omega_n\) shrinks by a factor of about \(7\) per level and \((L^{0.55})_{aa}\) increases toward \(\approx1.35247\), so \(\sigma_n-(L^{1-\alpha})_{aa}\to0\) is suggested, not proved.

**NOT THIS.** Theorem AP proves an exact identity for the combinatorial-to-fractional share and a level-uniform floor \(\sigma_n\ge1.310267\) at \(\alpha=0.45\). Combined with the existing fractional floors, it lifts the audit floor to \(0.274743\) (AM route) or \(0.415238\) (AO as stated). It does not prove that \(\|F(n,0.45)\|\) converges, does not prove that it decreases, and does not prove Conjecture J.1's value \(L\ge0.70\) (Kept Failure AP.1). It does not touch Stage 1 or Stage 2. A share of graph corner currents is not a force, and a positive floor on a graph audit norm is not thrust and not a pass of any experimental criterion. Not thrust.


## Theorem AQ. Spectral decimation of the corner resolvent: the audit norm decreases to a limit in \([0.819,0.823]\)

**ASSUMPTION.** Spectral calculus as in item 4. Corner data \((1,-1/2,0)\) on corners \(a,b,c\). \(G_n\) is the level-\(n\) combinatorial gasket of `build_gasket` (\(G_0\) is one triangle; \(G_{n+1}\) replaces every unit triangle of \(G_n\) by a copy of \(G_1\)), \(L=L_n\) its graph Laplacian, \(v=e_a-e_b\), \(P_C\) the projector onto the three corners, \(\beta=1-\alpha\), \(\gamma=\log3/\log5\). Write
\[
r_n(s)=\tfrac12\,v^{T}(L_n+s)^{-1}v\quad(s>0),\qquad R_{\alpha}(n)=v^{T}(L_n^{\alpha})^{+}v .
\]
Theorem AB (\(\|F\|=\sigma_n\|F_{\mathrm{frac}}\|\)), the pairing \(\|F_{\mathrm{frac}}\|=\sqrt{21}/(3R_{\alpha})\) (Theorems AA, AD, AE), and Theorem AP (1) (\(\sigma_n=(L^{\beta})_{aa}+\omega_n\), \(\omega_n=-(L^{\beta})_{ab}\)) are used as stated. This note does not re-prove them. Item (5) rests on an interval-arithmetic certificate (outward-rounded `mpmath.iv`, 30 digits); every other step is algebra.

**THEOREM AQ.**

1. **Exact decimation.** Eliminating the three midpoints of every level-\(1\) cell of \(G_{n+1}\) (Schur complement of \(L_{n+1}+s\) onto the vertices of \(G_n\)) gives exactly
   \[
   c(s)\,L_n+\bigl(2e(s)+s\bigr)I-e(s)P_C,\qquad c=\frac{6+s}{(5+s)(2+s)},\quad e=\frac{2s}{2+s},\quad \frac{2e+s}{c}=s(s+5).
   \]
   Consequently \(r_0(s)=1/(3+s)\) and
   \[
   r_{n+1}(s)=\Phi_s\bigl(r_n(s(s+5))\bigr),\qquad \Phi_s(x)=\frac{(2+s)(5+s)\,x}{(6+s)-2s(5+s)\,x}.
   \]
2. **Stieltjes representations.** For every \(\alpha\in(0,1)\) and \(n\ge0\),
   \[
   R_{\alpha}(n)=\frac{2\sin\pi\alpha}{\pi}\int_0^\infty s^{-\alpha}r_n(s)\,ds,\qquad
   \sigma_n=\tfrac12v^{T}L^{\beta}v=\frac{\sin\pi\beta}{\pi}\int_0^\infty s^{\beta-1}\bigl(1-s\,r_n(s)\bigr)\,ds .
   \]
3. **Monotonicity in the level (every \(\alpha\in(0,1)\)).** \(r_n(s)<r_{n+1}(s)\) for every \(s>0\) and \(n\ge0\). Hence \(R_{\alpha}(n)\) strictly increases, \(\|F_{\mathrm{frac}}(n,\alpha)\|=\sqrt{21}/(3R_{\alpha}(n))\) strictly decreases, \(\sigma_n\) strictly decreases, and
   \[
   \|F(n,\alpha)\|=\sigma_n\,\frac{\sqrt{21}}{3R_{\alpha}(n)}\quad\text{strictly decreases in }n\ge1 ,
   \]
   so \(\|F(n,\alpha)\|\) converges for every \(\alpha\in(0,1)\).
4. **Pointwise ceilings.** For every \(n\) and every \(k\ge0\), \(r_n\le T^k(1/\sigma)\), where \((Tf)(s)=\Phi_s(f(s(s+5)))\) and \(1/\sigma\) is the function \(s\mapsto1/s\); and for \(k\le n\), \(r_n\ge T^k r_0\). Near \(s=0\): if \(\beta'\in(1-\gamma,1-\alpha)\), \(\sigma_0>0\), \(s_*(s_*+5)=\sigma_0\), and \(B>0\) satisfy (C1) \((2+s_*)(5+s_*)^{1-\beta'}+2Bs_*^{1-\beta'}(5+s_*)^{1-\beta'}<6\), (C2) \(\Phi_s(X)\le Bs^{-\beta'}\) on \((s_*,\sigma_0]\) with \(X\ge\sup_n r_n(\sigma_0)\), and (C3) \(B\ge\sigma_0^{\beta'}/3\), then \(r_n(s)\le Bs^{-\beta'}\) on \((0,\sigma_0]\) for every \(n\).
5. **Certified values at \(\alpha=0.45\)** (\(\beta'=0.35\), \(\sigma_0=0.01\), \(B\approx0.892577\), \(X\approx2.555816\); \(18459\) geometric grid points on \([10^{-20},10^{12}]\)). For every level \(n\ge1\),
   \[
   R_{0.45}(n)<\bar R:=2.518955,\qquad \sigma_n>1.350986,\qquad \|F_{\mathrm{frac}}(n,0.45)\|>0.606412,
   \]
   and \(\|F(n,0.45)\|\) strictly decreases to a limit
   \[
   \mathcal L=\lim_{n\to\infty}\|F(n,0.45)\|\in[0.819255,\;0.822959].
   \]
   In particular \(\|F(n,0.45)\|>0.819255\) at every level, against Theorem AP's \(0.274743\) (AM route) and \(0.415238\) (AO as stated). This is Conjecture J.1 (decrease for every \(n\ge1\), existence of the limit, and \(\mathcal L\ge0.70\)), now proved. The ceiling \(\bar R\) is below Kept Failure AP.1's target \(2.859237\), so \(S_{AP}\sqrt{21}/(3\bar R)\approx0.794563\ge0.70\) already closes AP.1's open step without the sharper share bound.

**Proof.** *(1)* Every edge of \(G_n\) lies in exactly one unit triangle, so every midpoint added in \(G_{n+1}\) has degree \(4\), all its edges inside one cell, and old vertices are adjacent only to midpoints. The midpoint block of \(L_{n+1}+s\) is therefore block-diagonal over cells, and the Schur complement is a sum of cell terms plus \(sI\). For one cell (old vertices \(x_i\), midpoints forming a triangle) the midpoint block is \(M=(5+s)I-J\), \(M^{-1}=(5+s)^{-1}(I+J/(2+s))\), the coupling \(B\) has \(BB^{T}=I+J\) and \(BJB^{T}=4J\), so the cell term is \(2I-(5+s)^{-1}(I+J+4J/(2+s))=c(3I-J)+eI\) with \(c,e\) as stated. A corner lies in one cell and every other vertex in two, so the sum is \(cL_n+e(2I-P_C)+sI\). Direct algebra gives \((2e+s)/c=s(s+5)\). The corner block of an inverse is the inverse of the Schur complement, so the corner block of \((L_{n+1}+s)^{-1}\) is the corner block of \((c(L_n+s'I)-eP_C)^{-1}\), \(s'=s(s+5)\). The \(S_3\) symmetry of the corners makes the corner block of \((L_n+s')^{-1}\) equal to \(g_0I+g_1(J-I)\); on the zero-sum corner vectors it acts as \(r_n(s')=g_0-g_1\). The rank-three update \(-eP_C\) acts on that mode by \(\mu\mapsto\mu/(1-e\mu)\) (Woodbury on the corner block: \((G^{-1}-eI)^{-1}\)). With the factor \(1/c\) this is \(r_{n+1}(s)=r_n(s')/(c-e\,r_n(s'))=\Phi_s(r_n(s'))\). On the triangle \(G_0\), \(L_0=3I-J\) acts as \(3\) on zero-sum vectors, so \(r_0=1/(3+s)\).

*(2)* For \(\lambda>0\) and \(\alpha\in(0,1)\), \(\lambda^{-\alpha}=\frac{\sin\pi\alpha}{\pi}\int_0^\infty s^{-\alpha}(s+\lambda)^{-1}ds\) and \(\lambda^{\beta}=\frac{\sin\pi\beta}{\pi}\int_0^\infty s^{\beta-1}\lambda(s+\lambda)^{-1}ds\). Expand \(v\) (which is orthogonal to \(\ker L\)) in an eigenbasis of \(L\) and use \(\lambda/(s+\lambda)=1-s/(s+\lambda)\) and \(|v|^2=2\). For \(\sigma_n\): by Theorem AP (1) and the reflection swapping \(a,b\), \(\sigma_n=(L^{\beta})_{aa}-(L^{\beta})_{ab}=\frac12v^{T}L^{\beta}v\).

*(3)* \(r_n(s)\le|v|^2/(2s)=1/s\). Since \(r_n(s')\le1/s'\), the denominator of \(\Phi_s\) is at least \((6+s)-2s(5+s)/s'=4+s>0\), and \(\Phi_s'(x)=(2+s)(5+s)(6+s)/\mathrm{den}^2>0\) there, so \(\Phi_s\) is strictly increasing on \([0,1/s']\). Base: \(r_1=(5+s)/(s^2+7s+9)\) and \(r_1-r_0=(s+6)/((s+3)(s^2+7s+9))>0\). Induction: \(r_{n+1}(s)-r_n(s)=\Phi_s(r_n(s'))-\Phi_s(r_{n-1}(s'))>0\). The integrands in (2) are then strictly monotone in \(n\) (\(s\,r_n\) increases), which gives the monotonicity of \(R_{\alpha}(n)\) and \(\sigma_n\). \(\|F\|\) is a product of two positive strictly decreasing sequences, so it strictly decreases and converges.

*(4)* \(T\) is order-preserving on functions with \(0\le f\le1/\sigma\) (by (3)). Moreover \(T(1/\sigma)(s)=(2+s)/(s(4+s))<1/s\). Since \(r_m\le1/\sigma\) for every \(m\), \(r_{m+k}=T^kr_m\le T^k(1/\sigma)\); for \(n<k\), \(r_n\le r_k\le T^k(1/\sigma)\). Likewise \(r_n=T^kr_{n-k}\ge T^kr_0\) for \(k\le n\). For the power branch, argue by induction on \(n\). (C3) gives the base \(r_0\le1/3\le B\sigma_0^{-\beta'}\le Bs^{-\beta'}\). For \(s\le s_*\), we have \(s'\le\sigma_0\), and \(\Phi_s(Bs'^{-\beta'})\le Bs^{-\beta'}\) is equivalent, after dividing by \(Bs^{-\beta'}\) and using \(s'/s=5+s\), to \((2+s)(5+s)^{1-\beta'}+2Bs^{1-\beta'}(5+s)^{1-\beta'}\le6+s\). The left side increases in \(s\), so (C1) suffices. For \(s\in(s_*,\sigma_0]\), \(s'>\sigma_0\) and \(r_n(s')\le r_n(\sigma_0)\le X\), because \(r_n\) decreases in \(s\). So (C2) suffices, by monotonicity of \(\Phi_s\). The window \((1-\gamma,1-\alpha)\) is needed for (C1) at \(s\to0\) (\(2\cdot5^{1-\beta'}<6\)) and for integrability of \(s^{-\alpha-\beta'}\) at \(0\).

*(5)* \(r_n\) decreases in \(s\) and \(s\,r_n\) increases in \(s\). On each grid cell \([g,g']\), (4) gives
- \(r_n(s)\le r_n(g)\le T^k(1/\sigma)(g)\) for every \(n\);
- \(1-s\,r_n(s)\ge1-g'\,T^k(1/\sigma)(g')\) for every \(n\);
- \(r_n(s)\ge T^kr_0(g')\) and \(1-s\,r_n(s)\le1-g\,T^kr_0(g)\) for \(n\ge k\), where \(k\le36\) is the depth at which the orbit \(s\mapsto s(s+5)\) passes \(10^{60}\).

On \((0,10^{-20}]\), (4) gives \(\int s^{-\alpha}r_n\le B\,s_0^{1-\alpha-\beta'}/(1-\alpha-\beta')\) and \(\int s^{\beta-1}(1-sr_n)\le s_0^{\beta}/\beta\). On \([10^{12},\infty)\), \(r_n\le1/s\) and \(1-sr_n\le3/(3+s)\). Summing in outward-rounded interval arithmetic gives \(R_{\alpha}(n)\le2.518955\) for all \(n\), \(R_{\alpha}(n)\ge2.513145\) for \(n\ge36\), \(\sigma_n\ge1.350986\) for all \(n\), and \(\sigma_n\le1.353964\) for \(n\ge36\). (C1)–(C3) are checked in the same arithmetic (C1 left side \(5.789895<6\); C2 on \(400\) subintervals). Then \(\mathcal L=\lim\sigma_n\cdot\sqrt{21}/(3\lim R_{\alpha}(n))\), with \(\lim\sigma_n\ge1.350986\) and \(\lim R_{\alpha}\le2.518955\), gives the lower end \(0.819255\). Because \(\|F(n)\|\) decreases, \(\|F(n)\|\ge\mathcal L\) at every level. \(\sigma_n\le1.353964\) and \(R_{\alpha}(n)\ge2.513145\) for \(n\ge36\) give the upper end. \(\square\)

**Witness, computed by `scripts/gasket_fractional_decimation.py`.**
- Symbolic (SymPy): the one-cell Schur complement \(c(3I-J)+eI\); \((2e+s)/c=s(s+5)\); \(\Phi_s\), \(r_1\), \(r_1-r_0\), \(\Phi_s'\), \(\mathrm{den}(1/s')=4+s\), and \(T(1/\sigma)=(2+s)/(s(4+s))\).
- Numerical: on the actual graphs, \(\mathrm{Schur}(L_{n+1}+s\to V_n)=cL_n+(2e+s)I-eP_C\) for \(n=1\)–\(4\) at \(s\in\{0.03,0.7,4\}\) (error \(<10^{-10}\)). The recursion equals \(\frac12v^{T}(L_n+s)^{-1}v\) on levels \(0\)–\(6\).
- Live at \(0.45\) (levels \(1\)–\(6\)) and at \(\alpha\in\{0.10,0.25,0.49,0.55\}\) (levels \(1\)–\(5\)): both Stieltjes integrals of the recursion reproduce the dense \(R_{\alpha}(n)\) and \(\sigma_n\); \(\|F_{\mathrm{frac}}\|=\sqrt{21}/(3R)\); \(\|F\|=\sigma\|F_{\mathrm{frac}}\|\).
- Monotonicity: \(R\) strictly increases and \(\sigma\), \(\|F\|\) strictly decrease. \(r_n\le r_{n+1}\) for \(n<40\) on a log grid.
- Recorded values: the recursion reproduces the level-\(2\) lock \(1.165909\), the recorded matrix-free \(\|F_{\mathrm{frac}}(9)\|,\|F_{\mathrm{frac}}(10)\|\), and \(\|F(10)\|=0.833640\).
- Certificate at \(0.45\): (C1)–(C3) and the brackets of (5). The same certificate (grid ratio \(1.02\)) holds at \(\alpha\in\{0.10,0.25,0.49,0.55\}\), with live levels \(1\)–\(5\) inside.

The recursion gives, at \(0.45\):

| \(n\) | 1 | 2 | 3 | 4 | 5 | 6 | 10 | 20 | 40 |
|---|---|---|---|---|---|---|---|---|---|
| \(R_{\alpha}(n)\) | 1.523237 | 1.792023 | 2.005046 | 2.161151 | 2.271120 | 2.347382 | 2.478218 | 2.515053 | 2.515946 |
| \(\sigma_n\) | 1.435558 | 1.367791 | 1.354971 | 1.352845 | 1.352528 | 1.352483 | 1.352476 | 1.352476 | 1.352476 |
| \(\|F\|\) | 1.439599 | 1.165909 | 1.032272 | 0.956206 | 0.909692 | 0.880109 | 0.833640 | 0.821430 | 0.821139 |

The script exits nonzero if any check fails.

**KEPT FAILURE AQ.1 (the uniform ceiling and a positive limit do not extend to \(\alpha\ge\log3/\log5\)).** The claim "the decimation ceiling of (4)–(5) gives a level-uniform bound on \(R_{\alpha}(n)\), hence a positive limit of \(\|F(n,\alpha)\|\), for every \(\alpha\in(0,1)\)" is false. Item (3) holds for every \(\alpha\), but the window \((1-\gamma,1-\alpha)\) for \(\beta'\) is empty when \(\alpha\ge\gamma\). Indeed \(r_n(s)\approx(5/3)^n/3\) for \(s\lesssim5^{-n}\), so \(\int s^{-\alpha}r_n\) is unbounded in \(n\) there. The same recursion shows \(R_{\gamma}(n)\) growing linearly: per-level increment \(\approx0.5172\) on \(n=8\)–\(24\), Theorem AD's \(\Theta(n)\). At \(\alpha=0.75\) it grows geometrically: per-level ratio \(\approx1.125\) on \(n=20\)–\(24\), tending to \(5^{\alpha}/3\approx1.1146\). There \(\|F\|\to0\) (Theorems AD and K), so the decreasing limit of (3) is \(0\). Above \(\gamma\) the limit is \(0\), not a positive constant.

**OBSERVATION, not a theorem.** \(\sigma_n\) converges extremely fast (\(\sigma_{10}-\sigma_{20}\approx2.6\cdot10^{-9}\)), while \(R_{\alpha}(n)\) converges geometrically with per-level gap ratio near \(5^{\alpha}/3\approx0.688\). The recursion's value \(\|F(40,0.45)\|\approx0.821139\) sits inside the certified bracket. With this theorem the Lobatto/Radau share floors of Theorem AP and the corner-defect constants \(1/2\), \(11/20\) of Theorems AN, AO (witnessed on a dyadic range, not proved) are no longer needed for any audit floor.

**NOT THIS.** Theorem AQ is an exact renormalization identity for one corner resolvent on the combinatorial gasket graph. With an interval-arithmetic certificate, it proves that \(\|F(n,0.45)\|\) strictly decreases to a limit \(\mathcal L\in[0.819255,0.822959]\), so Conjecture J.1's statement holds. It does not identify \(\mathcal L\) in closed form. It does not give a positive limit for \(\alpha\ge\log3/\log5\) (Kept Failure AQ.1). It says nothing about the multi-scale or geometric-weight builds, and does not touch Stage 1 or Stage 2. The Stage 2 divergence question stays OPEN. \(\mathcal L\) is the limit of a graph audit norm of corner currents: not a force, not thrust, and not a pass of any experimental criterion. Not thrust.


## Theorem AR. Exact refinement rates: \(R_{\alpha}(n)\) moves by \(c_R(\alpha)(5^{\alpha}/3)^n\) per level

**Setting.** As in Theorem AQ: \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\), \(r_{n+1}(s)=\Phi_s(r_n(p(s)))\) with \(p(s)=s(s+5)\), \(r_0=1/(3+s)\). Write \(q=p^{-1}\) on \([0,\infty)\), \(q(u)=2u/(5+\sqrt{25+4u})\), \(\delta_n=r_n-r_{n-1}>0\) (AQ (3)), \(\rho=5^{\alpha}/3\), \(\mu=5^{\alpha-1}/3=\rho/5\), \(\gamma=\log3/\log5\), \(\beta=1-\alpha\). **Assumptions:** none beyond Theorem AQ, including its interval-certified power ceiling \(r_n(s)\le Bs^{-0.35}\) on \((0,0.01]\) (AQ (4)–(5)); that ceiling does not involve \(\alpha\).

**THEOREM AR.** For every \(\alpha\in(0,1)\):
1. **Product formula.** \(\delta_{n+1}(q^n(u))=P_n(u)\,\delta_1(u)\) for \(u>0\), with \(P_n(u)=\prod_{j=1}^{n}D_j(q^j(u))\), where \(D_j(s)\) is the secant slope of \(\Phi_s\) between \(r_{j-1}(p(s))\) and \(r_j(p(s))\):
   \[
   D_j(s)=\frac{(2+s)(5+s)(6+s)}{\mathrm{den}(r_{j-1}(ps))\,\mathrm{den}(r_j(ps))},\qquad \mathrm{den}(x)=(6+s)-2s(5+s)x ,
   \]
   and \(D_j(s)\ge(2+s)(5+s)/(6+s)=5/3+s(3s+16)/(3(6+s))\ge5/3\). Hence \((3/5)^nP_n(u)\) increases in \(n\) to a finite limit \(P_{\infty}(u)\).
2. **Exact rates.** With \(\Psi(u)=\lim_n5^nq^n(u)\) (the Koenigs function of \(q\) at its attracting fixed point \(0\), \(q'(0)=1/5\)),
   \[
   R_{\alpha}(n+1)-R_{\alpha}(n)=c_R(\alpha)\,\rho^{\,n}(1+o(1)),\qquad \sigma_n-\sigma_{n+1}=c_S(\alpha)\,\mu^{\,n}(1+o(1)),
   \]
   \[
   c_R=\tfrac{2\sin\pi\alpha}{\pi}\int_0^\infty\Psi^{-\alpha}\Psi'P_{\infty}\delta_1\,du,\qquad c_S=\tfrac{\sin\pi\beta}{\pi}\int_0^\infty\Psi^{1-\alpha}\Psi'P_{\infty}\delta_1\,du,\qquad 0<c_R,c_S<\infty .
   \]
3. **Share floor.** \(\sigma_n>\sigma_{\infty}\ge2\cdot4^{-\alpha}\) (at \(0.45\): \(1.071773\); live \(1.352476\)).
4. **Consequences.** (a) For \(\alpha<\gamma\) (\(\rho<1\)): \(R_{\infty}-R_{\alpha}(n)\sim\frac{c_R}{1-\rho}\rho^n\) and, since \(\mu=\rho/5\), \(\|F(n,\alpha)\|-\mathcal L(\alpha)\sim\frac{\sqrt{21}}{3}\frac{\sigma_\infty c_R}{R_\infty^2(1-\rho)}\rho^n\): the gap ratio is exactly \(5^{\alpha}/3\) in the limit (AQ's observation, now a theorem). (b) At \(\alpha=\gamma\): \(R_{\gamma}(n)=c_R(\gamma)\,n+o(n)\) and \(\|F(n,\gamma)\|\sim\frac{\sqrt{21}\,\sigma_\infty}{3c_R(\gamma)\,n}\), sharpening Theorem AD's \(\Theta(n)\) to an exact slope. (c) For \(\alpha>\gamma\): \(R_{\alpha}(n)\sim\frac{c_R}{\rho-1}\rho^n\) and \(\|F(n,\alpha)\|=\Theta\big((3\cdot5^{-\alpha})^n\big)\), with that exact ratio.

**Proof.** *(1)* \(\Phi_s\) is Möbius, so \(\Phi_s(x)-\Phi_s(y)=(x-y)(2+s)(5+s)(6+s)/(\mathrm{den}(x)\mathrm{den}(y))\) (SymPy). Then \(\delta_{n+1}(s)=D_n(s)\,\delta_n(p(s))\); iterate \(n\) times and put \(s=q^n(u)\), so the \(k\)-th factor sits at \(p^k(s)=q^{n-k}(u)\) with index \(n-k=j\). The floor uses \(0<\mathrm{den}\le6+s\) (AQ (3): \(\mathrm{den}\ge4+s\) on \([0,1/p(s)]\)). Upper control: \(p(s)\,r_j(p(s))\le\min(1,Bp(s)^{0.65})\) for \(p(s)\le0.01\), so \(\mathrm{den}\ge6+s-2B\,p(s)^{0.65}\) and \((3/5)D_j(s)\le1+O(s^{0.65})\) as \(s\to0\), uniformly in \(j\); everywhere \(\mathrm{den}\ge4+s\) gives \(D_j(s)\le M(s):=(2+s)(5+s)(6+s)/(4+s)^2\le6+s\). Since \(q^j(u)\le5^{-j}u\) (from \(q(u)=u/(5+q(u))\)), \(\sum_j\big((3/5)D_j(q^ju)-1\big)<\infty\), so \(P_\infty(u)<\infty\).

*(2)* By AQ (2), \(R_{\alpha}(n+1)-R_{\alpha}(n)=\tfrac{2\sin\pi\alpha}{\pi}\int s^{-\alpha}\delta_{n+1}(s)\,ds\). Substituting \(s=q^n(u)\) and using (1),
\[
\int s^{-\alpha}\delta_{n+1}\,ds=\rho^n\int_0^\infty\big(5^nq^n(u)\big)^{-\alpha}\,\big(5^n(q^n)'(u)\big)\,\big((3/5)^nP_n(u)\big)\,\delta_1(u)\,du .
\]
Each factor converges pointwise: \(5^nq^n\downarrow\Psi>0\) (because \(5q(u)\le u\)), \(5^n(q^n)'=\prod_{j<n}5q'(q^ju)=\prod_{j<n}5/\sqrt{25+4q^ju}\to\Psi'\), and \((3/5)^nP_n\uparrow P_\infty\). Domination: \((5^nq^n)^{-\alpha}\le\Psi^{-\alpha}\), with \(\Psi(u)\ge u\prod_j5/(5+q^ju)\) near \(0\) and \(\Psi\ge\Psi(1)\) on \([1,\infty)\). Pair the remaining factors as \(5q'(u)\cdot\prod_{1\le j<n}\big(3D_j(q^ju)q'(q^ju)\big)\cdot\tfrac35D_n(q^nu)\); each middle factor is \(1+O(s^{0.65})\) at small \(s\) (by the ceiling, since \(3\cdot\tfrac53\cdot\tfrac15=1\)), bounded on compacts, and \(\le3(6+s)/\sqrt{25+4s}=O(\sqrt s)\) at large \(s\). Because \(q(u)\le\sqrt u\), at most \(\log_2\log u+O(1)\) orbit points exceed \(1\) and their \(\sqrt{\cdot}\) factors multiply to \(O(u^{1/2})\) times a power of \(\log u\), cancelled by \(5q'(u)\le5/(2\sqrt u)\). So the integrand is dominated by \(C\,\Psi^{-\alpha}(\log(e+u))^{c}\delta_1(u)\), integrable since \(\delta_1(u)=(u+6)/((u+3)(u^2+7u+9))\sim u^{-2}\) and \(\Psi^{-\alpha}\lesssim u^{-\alpha}\) near \(0\). Dominated convergence gives the limit \(c_R\). For \(\sigma\): AQ (2) gives \(\sigma_n-\sigma_{n+1}=\tfrac{\sin\pi\beta}{\pi}\int s^{\beta}\delta_{n+1}\,ds\) and the same substitution produces the prefactor \(5^{-n\beta}\cdot5^{-n}\cdot(5/3)^n=(5^{\alpha-1}/3)^n=\mu^n\), with \((5^nq^n)^{\beta}\le u^{\beta}\) in the dominator (still integrable against \(u^{-2}\)). Positivity: every factor is positive.

*(3)* AQ (4) gives \(r_\infty(s)\le T(1/\sigma)(s)=(2+s)/(s(4+s))\), so \(1-s\,r_\infty\ge2/(4+s)\), and \(\tfrac{\sin\pi\beta}{\pi}\int_0^\infty s^{\beta-1}\tfrac{2}{4+s}ds=2\cdot4^{\beta-1}=2\cdot4^{-\alpha}\). Monotone convergence passes the bound to \(\sigma_\infty\).

*(4)* Sum the increments of (2); \(\sigma_n\to\sigma_\infty>0\) at the faster rate \(\mu^n\); \(\|F\|=(\sqrt{21}/3)\sigma_n/R_{\alpha}(n)\) (AQ (3)). \(\square\)

**Witness, computed by `scripts/gasket_fractional_decimation_rate.py`** (about 7 s).
- SymPy: the Möbius secant slope; \((2+s)(5+s)/(6+s)-5/3=s(3s+16)/(3(6+s))\).
- The product formula of (1) at \(u\in\{0.003,0.7,4,60\}\), \(n\in\{1,3,6\}\), relative error \(<10^{-48}\) at 50 digits.
- \((3/5)D_j\ge1\) on a 120001-point log grid for \(j=1\)–\(45\), and \(|(3/5)D_j-1|<10^{-12}\) for \(j\ge30\).
- \(c_R,c_S\) from the integral formulas (with \(\Psi,P_\infty\) truncated at depth 45) against the scaled increments of the exact recursion:

| \(\alpha\) | \(\rho=5^{\alpha}/3\) | \(c_R(\alpha)\) | \((R(21)-R(20))/\rho^{20}\) | \(c_S(\alpha)\) |
|---|---|---|---|---|
| 0.25 | 0.498450 | 0.37699178 | 0.37699179 | 1.082999 |
| 0.45 | 0.687726 | 0.49780265 | 0.49780265 | 0.927091 |
| 0.60 | 0.875509 | 0.52084941 | 0.52084941 | 0.642994 |
| \(\gamma\approx0.682606\) | 1 | 0.51716583 | 0.51716583 | 0.481399 |
| 0.75 | 1.114567 | 0.50798513 | 0.50798513 | 0.357291 |

  Agreement is better than \(10^{-6}\) for \(c_R\) and \(10^{-3}\) for \(c_S\) (whose scaled sequence is read at \(n=10\)). At \(\alpha=\gamma\) the slope \(c_R(\gamma)=0.517166\) is the per-level increment that Kept Failure AQ.1 measured as \(\approx0.5172\); Theorem AD's proved floor slope was \(c_{AD}=0.051634\). At \(0.45\) the recursion reproduces the recorded \(\|F(10)\|=0.833640\), gives \(\mathcal L\approx0.821139\) (inside AQ's bracket), and \((\|F(22)\|-\mathcal L)/(\|F(21)\|-\mathcal L)\) is within \(2\cdot10^{-3}\) of \(\rho\).
- \(\sigma_{11}\ge2\cdot4^{-\alpha}\) at every tested \(\alpha\).

The constants \(c_R,c_S\) are computed by quadrature in floating point, not interval-certified; the theorem asserts their existence and positivity, not their digits.

**KEPT FAILURE AR.1 (the rate does not follow from the universal resolvent bound alone).** The claim "the a-priori bound \(r_n\le1/s\) (equivalently \(\mathrm{den}\ge4+s\)) suffices to prove a geometric rate for \(R_{\alpha}(n)\)" is false as a proof route: it gives only \((3/5)D_j(s)\le(3/5)\cdot60/16=2.25\) near \(s=0\), so \((3/5)^nP_n\) is not controlled and no rate follows. The proof needs AQ's certified power ceiling, which forces \(\mathrm{den}\to6\) at small \(s\) and makes \(\prod(3/5)D_j\) converge. Witness: the same script.

**NOT THIS.** Theorem AR gives exact geometric rates (and, at \(\gamma\), an exact linear slope) for one graph resistance and one share on the combinatorial gasket. It does not give \(c_R\), \(c_S\), or \(\mathcal L\) in closed form, does not certify their digits, and says nothing about the multi-scale or geometric-weight builds. It does not touch Stage 1 or Stage 2; the Stage 2 divergence question stays OPEN. Not thrust.

## Theorem AS. Certified digits of the audit limit: \(\mathcal L(0.45)=0.8211386271335\ldots\), and Conjecture J.2 holds

**Setting.** As in Theorem AQ: \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\), \(w_n(s)=1-s\,r_n(s)\), \(r_{n+1}(s)=\Phi_s(r_n(s(s+5)))\), \(r_0=1/(3+s)\), \(T f(s)=\Phi_s(f(s(s+5)))\). Write \(r_\infty=\lim r_n\) (increasing, AQ (3)) and \(w_\infty=1-s\,r_\infty\). **Assumptions:** Theorem AQ (1)–(4) as stated (decimation, Stieltjes formulas, monotonicity, pointwise ceilings including the \(\alpha\)-free power ceiling \(r_n\le Bs^{-0.35}\) on \((0,0.01]\), re-certified here: \(B\approx0.892577\)). Items (4)–(5) rest on an outward-rounded `mpmath.iv` computation at 50 digits; every other step is algebra.

**THEOREM AS.** For every \(\alpha\in(0,\gamma)\):
1. **Limit integrals.** By monotone convergence (\(r_n\uparrow\)) and dominated convergence (\(0\le w_n\le w_0=3/(3+s)\)), and since \(\sin\pi\beta=\sin\pi\alpha\),
   \[
   R_\infty(\alpha)=\lim_n R_{\alpha}(n)=\tfrac{2\sin\pi\alpha}{\pi}J_r,\quad \sigma_\infty(\alpha)=\lim_n\sigma_n=\tfrac{\sin\pi\alpha}{\pi}J_w,\quad \mathcal L(\alpha)=\lim_n\|F(n,\alpha)\|=\frac{\sqrt{21}}{6}\,\frac{J_w}{J_r},
   \]
   with \(J_r=\int_0^\infty s^{-\alpha}r_\infty\,ds\), \(J_w=\int_0^\infty s^{-\alpha}w_\infty\,ds\). Both integrals use the same weight \(s^{-\alpha}\).
2. **Complete monotonicity.** \(r_n=\sum_ic_i/(\lambda_i+s)\) and \(w_n=\sum_ic_i\lambda_i/(\lambda_i+s)\) with \(c_i=\tfrac12(v^{T}\phi_i)^2\ge0\), \(\sum c_i=1\). So every divided difference of order \(k\) of \(r_n\) or \(w_n\) has sign \((-1)^k\), and this weak sign passes to the pointwise limits \(r_\infty,w_\infty\).
3. **Interpolation sandwich.** Let \(f\in\{r_\infty,w_\infty\}\) and a cell \((x_c,x_{c+1})\) with \(b\) interpolation nodes \(\le x_c\) and \(u\) nodes \(\ge x_{c+1}\), and let \(p\) be the interpolant. Then \(\operatorname{sign}(f-p)=(-1)^b\) on the cell: \(p\) is a floor for even \(b\) and a ceiling for odd \(b\).
4. **Certified values** (geometric nodes \(x_k=2^{-250}(17/16)^k\) up to \(10^{25}\), \(3831\) nodes; floor \(b=u=10\), ceiling \(b=u=11\)):

| \(\alpha\) | \(R_\infty(\alpha)\) | \(\sigma_\infty(\alpha)\) | \(\mathcal L(\alpha)=\lim\|F(n,\alpha)\|\) |
|---|---|---|---|
| 0.10 | 1.971070823065 | 1.820244396117 | \([1.410638933063,\,1.410638933065]\) |
| 0.25 | 2.022380229615 | 1.591881227155 | 1.202364770316 |
| 0.45 | \([2.51594623787706,\,2.51594623787714]\) | \([1.35247562325228,\,1.35247562325233]\) | \([0.82113862713346,\,0.82113862713352]\) |
| 0.49 | \([2.779225340768,\,2.779225340771]\) | 1.312402359293 | \([0.7213260790662,\,0.7213260790670]\) |
| 0.55 | \([3.5020566,\,3.5020569]\) | 1.256764544853 | \([0.54817486,\,0.54817490]\) |

   Every entry is an enclosure; where one number is printed, both endpoints agree to the digits shown (widths \(\le3\cdot10^{-12}\), except \(\alpha=0.55\), whose width \(2.6\cdot10^{-8}\) is set by the head term \(B\,s_0^{0.65-\alpha}=B\cdot2^{-25}\)). At \(0.45\) this narrows AQ's bracket \([0.819255,0.822959]\) (width \(3.7\cdot10^{-3}\)) to width \(5\cdot10^{-14}\): \(\mathcal L=0.8211386271335\ldots\) (13 digits).
5. **Conjecture J.2 is true.** \(\|F(n,0.10)\|\) strictly decreases for \(n\ge1\) (AQ (3)) to \(\mathcal L(0.10)>1.410638933063>1.410\). So \(\|F(n,0.10)\|>1.410\) at every level, with margin \(6.4\cdot10^{-4}\) in the limit.

**Proof.** *(1)* AQ (2) and (3); \(\int_0^\infty s^{-\alpha}\,3/(3+s)\,ds<\infty\) for \(\alpha\in(0,1)\). *(2)* Spectral expansion of \(v\perp\ker L\), \(|v|^2=2\). For a \(C^\infty\) function the order-\(k\) divided difference equals \(f^{(k)}(\xi)/k!\), and \((-1)^k\partial_s^k(\lambda+s)^{-1}>0\). Divided differences at fixed nodes are finite linear combinations of values, so they pass to pointwise limits. *(3)* For any function, \(f(x)-p(x)=f[x_{1},\dots,x_{m},x]\prod_i(x-x_i)\) (Newton remainder). With \(m=b+u\) nodes, the divided difference has order \(b+u\), so sign \((-1)^{b+u}\). The product has \(b\) positive and \(u\) negative factors, so sign \((-1)^u\). Together, \((-1)^b\). *(4)* On a geometric grid, \(\int_{x_c}^{x_{c+1}}s^{-\alpha}p(s)\,ds=x_c^{1-\alpha}\sum_mW_mf(x_{c+m})\), where \(W_m=\int_1^{17/16}y^{-\alpha}\ell_m(y)\,dy\) are the Lagrange weights on the nodes \((17/16)^m\). These weights are computed once, in interval arithmetic, from exact monomial moments. Node values are enclosed by AQ (4): \(T^kr_0(x)\le r_\infty(x)\le T^k(1/\sigma)(x)\), with \(k\) the depth at which the orbit of \(s\mapsto s(s+5)\) passes \(10^{30}\); then \(w_\infty=1-xr_\infty\) is enclosed as well. The two enclosures agree to \(\sim10^{-25}\). The interval product \(\sum_mW_m[f]_m\) is a valid enclosure whatever the signs of \(W_m\). Head \((0,x_0]\), \(x_0=2^{-250}\): \(f\) decreases, so \(f\ge f(x_0)\). Above, \(r_\infty\le Bs^{-0.35}\) (AQ (4)) and \(w_\infty\le1\). Tail \([x_N,\infty)\), \(x_N>10^{25}\): \(1/(\lambda+s)\ge1/s-\lambda/s^2\) and \(\lambda/(\lambda+s)\ge\lambda/s-\lambda^2/s^2\), with \(\tfrac12v^{T}L_nv=2\) (\(n\ge1\)) and \(\tfrac12v^{T}L_n^2v=(L_n^2)_{aa}=d_a^2+d_a=6\) (\(n\ge2\): corners at graph distance \(2^n\ge4\)). These give \(r_\infty\in[1/s-2/s^2,1/s]\) and \(w_\infty\in[2/s-6/s^2,2/s]\). The power window \((1-\gamma,1-\alpha)\ni0.35\) needs \(\alpha<0.65\), which covers every tabulated \(\alpha\). Outward-rounded sums give the table, and \(\mathcal L\) follows from (1). *(5)* AQ (3) holds for every \(\alpha\in(0,1)\), and the strict decrease gives \(\|F(n)\|>\mathcal L\). \(\square\)

**Witness, computed by `scripts/gasket_fractional_limit_digits.py`** (about 45 s).
- Live: \(\tfrac12v^{T}Lv=2\) and \(\tfrac12v^{T}L^2v=6\) at levels \(2,3\); \(c_i\ge0\) and \(\sum c_i=1\); the sign rule of (3) for \(r_n,w_n\) (\(n=3,8\); \(b=u\in\{2,3,10,11\}\)) at 60 digits.
- AQ's power ceiling (C1)–(C3) re-certified (\(B\approx0.892577\); C1 left side \(5.789895<6\)).
- Certified enclosures of the table, with widths \(<10^{-10}\) (\(<10^{-7}\) at \(0.55\)). At \(0.45\): \(\mathcal L\) lies inside AQ's bracket; \(R_\infty<2.518955\); \(\sigma_\infty>1.350986\); and \(\sigma_\infty\ge2\cdot4^{-\alpha}\) (AR (3)).
- Consistency, not part of the proof: the exact recursion by 30-digit quadrature gives \(\|F(60,0.45)\|=0.821138627224937\), so \(0<\|F(60)\|-\mathcal L=9.1\cdot10^{-11}\) (AR predicts order \(\rho^{60}\approx1.7\cdot10^{-10}\)). Float \(\|F(30,\alpha)\|\) sits above the certified \(\mathcal L(\alpha)\) at \(\alpha\in\{0.10,0.25,0.49,0.55\}\). The recorded \(\|F(n,0.10)\|\), \(n=1\)–\(6\), decrease and stay above \(\mathcal L(0.10)\).

**KEPT FAILURE AS.1 (the exact ratio does not make level-10 extrapolation precise).** The claim "with Theorem AR's exact gap ratio \(\rho=5^{0.45}/3\approx0.687726\), one-step geometric extrapolation of the recorded levels \(9,10\) recovers \(\mathcal L(0.45)\) to \(10^{-4}\)" is false. \(\|F(10)\|-(\|F(9)\|-\|F(10)\|)\rho/(1-\rho)=0.820862\), which misses the certified \(\mathcal L=0.8211386\) by \(2.76\cdot10^{-4}\). Aitken's \(\Delta^2\) on levels \(8,9,10\) gives \(0.821526\), off by \(3.87\cdot10^{-4}\). The effective ratio \((\|F(10)\|-\mathcal L)/(\|F(9)\|-\mathcal L)=0.68298\) still differs from \(\rho\) by \(\approx0.7\%\): AR's \(o(1)\) has not settled by level \(10\). AR fixes the ratio only asymptotically, so finite-level data plus AR do not certify digits; the certificate of (4) does. Witness: the same script.

**OPEN.** Closed forms for \(\mathcal L\), \(R_\infty\), \(\sigma_\infty\); certified digits of AR's \(c_R(0.45)\), \(c_S(0.45)\) (the sandwich bounds limits of resistances, not the rate constants, whose integrands involve the Koenigs function \(\Psi\) and \(P_\infty\)); certified digits for \(\alpha\ge0.65\) (outside the power window used here); and any \(\alpha\ge\gamma\), where the limit is \(0\) (AQ.1).

**NOT THIS.** Theorem AS gives certified decimal digits of limits of one graph resistance, one share, and one audit norm of corner currents on the combinatorial gasket, and it settles Conjecture J.2 (\(\alpha=0.10\)). It gives no closed form. It does not certify AR's rate constants, and says nothing about the multi-scale or geometric-weight builds. It does not touch Stage 1 or Stage 2; the Stage 2 divergence question stays OPEN. A certified digit of \(\mathcal L\) is a digit of a graph norm, not a force: it is not thrust and not a pass of any experimental criterion. Not thrust.

## Theorem AT. Certified refinement-rate constants: \(c_R(0.45)=0.4978026507635603\ldots\), \(c_S(0.45)=0.9270926385694613\ldots\), slope \(c_R(\gamma)=0.5171658346305052\ldots\)

**Setting.** As in Theorems AQ–AS: \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\), \(r_{n+1}(s)=\Phi_s(r_n(p(s)))\), \(p(s)=s(s+5)\), \(q=p^{-1}\), \(r_0=1/(3+s)\), \(\delta_{n+1}=r_{n+1}-r_n>0\), \(\rho=5^{\alpha}/3\), \(\mu=\rho/5\), \(\gamma=\log3/\log5\). \(D_j(y)\) is AR's Möbius secant slope of \(\Phi_y\) between \(r_{j-1}(p(y))\) and \(r_j(p(y))\), and \(c_R(\alpha),c_S(\alpha)\) are AR's constants: \(R_\alpha(m+1)-R_\alpha(m)=c_R\rho^m(1+o(1))\) and \(\sigma_m-\sigma_{m+1}=c_S\mu^m(1+o(1))\). **Assumptions:** Theorem AQ (1)–(4) (decimation, Stieltjes formulas, monotonicity, pointwise ceilings including the \(\alpha\)-free power ceiling \(r_m\le Bs^{-0.35}\) on \((0,0.01]\), re-certified here with \(B\approx0.892577\)), Theorem AR (1)–(2) (product formula, \((3/5)D_j\ge1\), existence of \(c_R,c_S\)), and Theorem AS (2)–(3) (complete monotonicity and the interpolation sandwich), all as stated. Item (3) rests on an outward-rounded `mpmath.iv` computation at 60 digits; items (1), (2), (4) are algebra and estimates.

**THEOREM AT.** For every \(\alpha\in(0,1)\) and every \(n\ge0\):

1. **Exact transfer.** For \(m=n+k\), \(k\ge0\),
   \[
   \rho^{-m}\!\int_0^\infty\! s^{-\alpha}\delta_{m+1}\,ds=\rho^{-n}\!\int_0^\infty\! t^{-\alpha}\delta_{n+1}(t)\,A_k^{\alpha}B_kC_k\,dt,\qquad
   \mu^{-m}\!\int_0^\infty\! s^{1-\alpha}\delta_{m+1}\,ds=\mu^{-n}\!\int_0^\infty\! t^{1-\alpha}\delta_{n+1}(t)\,A_k^{\alpha-1}B_kC_k\,dt,
   \]
   with \(A_k(t)=\prod_{i=1}^{k}(1+q^i(t)/5)\), \(B_k(t)=\prod_{i=0}^{k-1}5/\sqrt{25+4q^i(t)}\), \(C_k(t)=\prod_{i=1}^{k}\tfrac35D_{n+i}(q^i(t))\).
2. **Level-\(n\) bracket.** Put \(A_\infty=\lim A_k\), \(B_\infty=\lim B_k\), and \(\bar C(t)=\prod_{i\ge1}\tfrac35\bar D(q^i(t))\), where \(\bar D(y)=(2+y)(5+y)(6+y)/\underline{\mathrm{den}}(y)^2\), \(\underline{\mathrm{den}}(y)=\max\bigl(6+y-2\,\overline{zr}(p(y)),\,4+y\bigr)\), and \(\overline{zr}(z)\ge z\,r_\infty(z)\) is AQ (4)'s ceiling (\(Bz^{0.65}\) for \(z\le0.01\), \(z\,T^k(1/\sigma)(z)\) above, capped at \(1\)). Then
   \[
   \frac{\pi\,\rho^{n}}{2\sin\pi\alpha}\,c_R\in\Bigl[\int t^{-\alpha}\delta_{n+1}B_\infty,\ \int t^{-\alpha}\delta_{n+1}A_\infty^{\alpha}\bar C\Bigr],\qquad
   \frac{\pi\,\mu^{n}}{\sin\pi\alpha}\,c_S\in\Bigl[\int t^{1-\alpha}\delta_{n+1}A_\infty^{\alpha-1}B_\infty,\ \int t^{1-\alpha}\delta_{n+1}\bar C\Bigr].
   \]
3. **Certified values** (\(n=40\); nodes \(5^{-40}2^{-110}(33/32)^k\); floor \(b=u=12\), ceiling \(b=u=13\)):

| \(\alpha\) | \(c_R(\alpha)\) | \(c_S(\alpha)\) |
|---|---|---|
| 0.25 | \([0.376991788570051826,\,0.376991788570051827]\) | \([1.083004227958476886,\,1.083004227958476902]\) |
| 0.45 | \([0.497802650763560354,\,0.497802650763560356]\) | \([0.927092638569461377,\,0.927092638569461384]\) |
| 0.60 | \([0.520849411207596441,\,0.520849411207596442]\) | \([0.642994619990388231,\,0.642994619990388235]\) |
| \(\gamma\) | \([0.517165834630505243,\,0.517165834630505245]\) | \([0.481399003768729118,\,0.481399003768729121]\) |
| 0.75 | \([0.507985128123073180,\,0.507985128123073181]\) | \([0.357290786593579090,\,0.357290786593579092]\) |

   Every width is below \(2\cdot10^{-17}\). The transfer works for every \(\alpha\in(0,1)\): unlike Theorem AS's limits, these constants do not need the power window \(\alpha<0.65\), so \(\alpha=0.75>\gamma\) is covered.
4. **Consequences.**
   (a) *Rate of the rate (DERIVED).* \(\rho^{-n}(R_\alpha(n+1)-R_\alpha(n))=c_R+O(5^{-0.65n})\) and \(\mu^{-n}(\sigma_n-\sigma_{n+1})=c_S+O(5^{-0.65n})\) for \(\alpha\ne0.65\) (use any exponent \(<0.65\) at \(0.65\)). At \(\alpha=\gamma\) the \(o(n)\) of AR (4b) improves to a constant: \(R_\gamma(n)=c_R(\gamma)\,n+C_\gamma+O(5^{-0.65n})\) with \(C_\gamma=R_\gamma(0)+\sum_{k\ge0}(R_\gamma(k+1)-R_\gamma(k)-c_R(\gamma))\) finite, and so \(\|F(n,\gamma)\|=\frac{\sqrt{21}\,\sigma_\infty}{3\,c_R(\gamma)\,n}\bigl(1-\frac{C_\gamma}{c_R(\gamma)n}+O(n^{-2})\bigr)\).
   (b) *Certified approach constants at \(0.45\).* With Theorem AS's enclosures of \(R_\infty,\sigma_\infty\): \(R_\infty-R_{0.45}(n)\sim\frac{c_R}{1-\rho}\rho^n\) with \(\frac{c_R}{1-\rho}\in[1.59411976454674,\,1.59411976454675]\), and \(\|F(n,0.45)\|-\mathcal L\sim c_F\rho^n\) with
   \[
   c_F=\frac{\sqrt{21}}{3}\,\frac{\sigma_\infty\,c_R}{R_\infty^2(1-\rho)}\in[0.52027873061815,\,0.52027873061821].
   \]

**Proof.** *(1)* Iterating AR (1)'s step \(\delta_{j+1}(s)=D_j(s)\,\delta_j(p(s))\) \(k\) times gives \(\delta_{m+1}(q^k t)=\prod_{i=1}^{k}D_{n+i}(q^i t)\,\delta_{n+1}(t)\). Substitute \(s=q^k(t)\). From \(q(x)(5+q(x))=x\), \(5^kq^k(t)/t=\prod_{i=1}^{k}5/(5+q^i t)=A_k^{-1}\); from \(q'(x)=1/\sqrt{25+4x}\), \(5^k(q^k)'(t)=B_k\). The scalar factors collect to \(\rho^{-k}\cdot5^{k\alpha}\cdot5^{-k}\cdot(5/3)^k=1\) for the first identity and \(\mu^{-k}\cdot5^{-k(1-\alpha)}\cdot5^{-k}\cdot(5/3)^k=1\) for the second.

*(2)* \(A_k\) increases to \(A_\infty\) and \(B_k\) decreases to \(B_\infty>0\). AR (1) gives \(\tfrac35D_j\ge1\). Each \(D_j(y)\) is \((2+y)(5+y)(6+y)\) over a product of two denominators \(\mathrm{den}(x)=(6+y)-2y(5+y)x\) with \(x\le r_\infty(p(y))\), and \(y(5+y)x=p(y)\,x\le p(y)\,r_\infty(p(y))\le\overline{zr}(p(y))\) (\(z\,r_\infty(z)\) increases in \(z\) by complete monotonicity, and \(\mathrm{den}\ge4+y\) by AQ (3)); so \(D_j\le\bar D\) for every \(j\). Hence \(B_\infty\le A_k^{\alpha}B_kC_k\le A_\infty^{\alpha}\bar C\) and \(A_\infty^{\alpha-1}B_\infty\le A_k^{\alpha-1}B_kC_k\le\bar C\) for every \(k\). Multiply by \(t^{-\alpha}\delta_{n+1}>0\) (resp. \(t^{1-\alpha}\delta_{n+1}\)) and integrate. By AQ (2), \(\rho^{-m}(R(m+1)-R(m))=\frac{2\sin\pi\alpha}{\pi}\rho^{-m}\int s^{-\alpha}\delta_{m+1}\), which therefore lies in the bracket for every \(m\ge n\), and so does its limit \(c_R\) (AR (2)); the same for \(c_S\) with \(\sin\pi\beta=\sin\pi\alpha\). The upper bracket is finite by the tail estimate in (3).

*(3)* All arithmetic is outward-rounded. *Infinite products.* \(q(x)\le x/5\) and \(q(x)\le\sqrt x\). The products run explicitly until \(q^i(t)<10^{-90}\). The remainder uses \(\log\tfrac35\bar D(y)\le y+2\varepsilon/(1-\varepsilon)\le1.82\,y^{0.65}\) for \(y\le10^{-3}\), where \(\varepsilon=2B\,p(y)^{0.65}/(6+y)\le0.0095\) (checked in interval arithmetic). It also uses \(\log(1+y/5)\le y/5\) and \(\tfrac12\log(1+4y/25)\le2y/25\), summed geometrically. *Cells.* \(B_\infty\) and \(A_\infty^{-1}\) decrease in \(t\), so they are taken at the right end of each cell; \(\bar C\) is bounded over the whole cell with \(y\) in the image interval. *Interpolation region* \((x_0,\,1.2\cdot10^4\cdot5^{-40}]\): \(r_n\) and \(r_{n+1}\) are completely monotone (finite sums \(c_i/(\lambda_i+s)\)), so Theorem AS (3)'s floors and ceilings with exact node values (the recursion of AQ (1) in interval arithmetic) give \(\int_{\rm cell}t^{-\alpha}\delta_{n+1}\in[\mathrm{floor}(r_{n+1})-\mathrm{ceil}(r_n),\,\mathrm{ceil}(r_{n+1})-\mathrm{floor}(r_n)]\), clipped below at \(0\). *Head* \((0,x_0]\): \(r_m\) decreases and \(r_m(0)=(5/3)^m/3\), so \(r_{n+1}(x_0)-r_n(0)\le\delta_{n+1}\le r_{n+1}(0)-r_n(x_0)\). *Product region* \([1.2\cdot10^4\cdot5^{-40},1]\): \(\delta_{n+1}(s)=\prod_{k<n}D_{n-k}(p^ks)\,\delta_1(p^ns)\le\prod_{k<n}\bar D(p^ks)\,\delta_1(p^ns)\), with \(\delta_1\) decreasing, bounded cellwise. Lower bound \(0\). The contribution is below \(10^{-25}\) of the integral at every tabulated \(\alpha\). *Tail* \([1,\infty)\): \(\bar D\le M(y)\le6+y\le1.2(5+y)\) gives \(\delta_{n+1}(s)\le2\cdot1.2^n/(s\,p^n(s))\). Counting orbit points above \(10^{-3}\) (at most \(9+\ln s\)) gives \(A_\infty^{\alpha}\bar C\le e^{24.05}s^{4.67}\le3\cdot10^{10}s^5\). With \(p^n(s)\ge\sqrt{s\,p^n(1)}\,s^{2^{n-1}}\), the tail is below \(10^{-2.5\cdot10^{11}}\). The resulting brackets are the table.

*(4a)* By (2), \(|\rho^{-n}\Delta R_n-c_R|\le\frac{2\sin\pi\alpha}{\pi}\rho^{-n}\int t^{-\alpha}\delta_{n+1}(A_\infty^{\alpha}\bar C-B_\infty)\). The weight difference is \(\le C\,t^{0.65}\) on \((0,1]\): near \(0\) it is \(O(t^{0.65})\) by the remainder estimates, and it is bounded on \([10^{-3},1]\). On \([1,\infty)\) it is \(\le3\cdot10^{10}t^5\). So the error is \(\le C\,5^{-0.65n}\cdot\rho'^{-n}\int t^{-\alpha'}\delta_{n+1}\) with \(\alpha'=\alpha-0.65\), \(\rho'=5^{\alpha'}/3\), plus a doubly-exponentially small tail. The last integral, rescaled, is bounded in \(n\) by AR (2) (an \(R\)-type integral if \(\alpha'>0\), a \(\sigma\)-type one with \(\beta=-\alpha'\) if \(\alpha'<0\)). Summing \(\Delta R_k-c_R\) at \(\rho=1\) gives \(C_\gamma\); expand \(\|F\|=(\sqrt{21}/3)\sigma_n/R_\gamma(n)\) with \(\sigma_n-\sigma_\infty=O(\mu^n)\). *(4b)* AR (4a) and arithmetic on the certified intervals. \(\square\)

**Witness, computed by `scripts/gasket_fractional_rate_constants.py`** (about 2.5 min).
- AQ's power ceiling re-certified (\(B\approx0.892577\), C1 left side \(5.789895<6\)); the small-\(y\) remainder lemma (ratio \(1.7896\le1.82\)); the crude-tail constants.
- Identity (1) at finite \(m\) by 40-digit quadrature: \(\alpha=0.45\), \(n=2\), \(m=5\) (relative error \(8.5\cdot10^{-24}\)); \(\alpha=0.75\), \(n=3\), \(m=5\) (\(1.2\cdot10^{-18}\)).
- \(n=40\): 2811 nodes, 2784 interpolation cells, 1786 product-bound cells. Node enclosures satisfy \(0<r_{40}\le r_{41}\). Contributions to the width at \(0.45\): head \(3\cdot10^{-59}\), product cells \(2\cdot10^{-33}\), sandwich gap \(8\cdot10^{-28}\) (absolute, before the factor \(\rho^{-40}\)).
- The table. AR's float \(c_R\) agree to \(10^{-8}\). AR's float \(c_S\) are off by \(5.2\cdot10^{-6}\) (\(0.25\)), \(1.6\cdot10^{-6}\) (\(0.45\)), \(6.2\cdot10^{-7}\) (\(0.60\)), \(3.8\cdot10^{-9}\) (\(\gamma\)), and \(-2.1\cdot10^{-7}\) (\(0.75\)), inside AR's stated \(10^{-3}\) agreement. AR's \(c_S(0.45)\approx0.927091\) is correctly \(0.9270926385694613\ldots\).
- Consistency at \(0.45\) by 40-digit quadrature of the exact recursion: \(\rho^{-m}\Delta R_m-c_R=-6.21\cdot10^{-3},-7.71\cdot10^{-5},-9.50\cdot10^{-7},-1.17\cdot10^{-8}\) at \(m=4,8,12,16\), and \(\mu^{-m}\Delta\sigma_m-c_S=-4.12\cdot10^{-2},\dots,-8.25\cdot10^{-8}\). Both are below \(5\cdot5^{-0.65m}\). **OBSERVATION, not a theorem:** each 4-level step divides the deviation by \(\approx81=3^4\), so the true correction is \(\Theta(3^{-m})\), faster than the proved \(O(5^{-0.65m})\). The proved exponent \(0.65\) comes from AQ's ceiling \(Bs^{-0.35}\); \(3^{-m}=5^{-\gamma m}\) is what \(r_\infty\asymp s^{\gamma-1}\) would give. \((\|F(30)\|-\mathcal L)/\rho^{30}=0.5202831\) against the certified \(c_F=0.5202787\ldots\).

The script exits nonzero if any check fails.

**KEPT FAILURE AT.1 (the slope at \(\gamma\) is not a pure-power constant).** The claim "\(r_\infty(s)\sim C\,s^{\gamma-1}\) as \(s\to0\) for one constant \(C\), so that \(c_R(\gamma)=\frac{2\sin\pi\gamma}{\pi}\ln5\cdot C\)" is false. Put \(g(s)=s^{1-\gamma}r_\infty(s)\). The fixed-point relation \(r_\infty=\Phi_s(r_\infty\circ p)\) and \(5^{\gamma}=3\) give \(\log\bigl(g(y)/g(p(y))\bigr)\in[-(1-\gamma)y/5,\ 1.82\,y^{0.65}]\) for \(y\le10^{-3}\). So \(\Theta(u)=\lim_kg(q^ku)\) exists, \(\Theta\circ q=\Theta\), and \(\Theta(s)\in g(s)\exp\bigl([-(1-\gamma)q(s)/4,\ 2\,q(s)^{0.65}/(1-5^{-0.65})]\bigr)\). At \(s=5^{-40}5^{j/8}\), \(j=0,\dots,7\), with \(r_\infty(s)\) enclosed by AQ (4) (widths \(<10^{-18}\)): \(\Theta\approx0.600954332,\ 0.600949127,\ 0.600946924,\ 0.600949013,\ 0.600954170,\ 0.600959375,\ 0.600961578,\ 0.600959489\). The certified enclosures are separated by more than \(1.465\cdot10^{-5}\). If \(\lim_{s\to0}g\) existed, \(\Theta\) would be constant, so it does not exist: \(s^{1-\gamma}r_\infty\) is asymptotically log-periodic with relative amplitude \(\approx2.4\cdot10^{-5}\), and no single constant \(C\) exists. Witness: the same script.

**CONJECTURE AT.2 (log-mean slope). Status: PROVED by Theorem AU (the text below records it as first stated).** \(c_R(\gamma)=\frac{2\sin\pi\gamma}{\pi}\,\ln5\,\langle\Theta\rangle\), where \(\langle\Theta\rangle\) is the mean of \(\Theta\) over one period in \(\log_5 s\). The evidence is a \(32\)-point trapezoid of \(g\) over \(s\in5^{-40}[1,5)\) (40-digit float, not interval), which gives \(0.517165834630505244\). That matches the certified \(c_R(\gamma)\in[0.517165834630505243,0.517165834630505245]\) to \(5\cdot10^{-22}\). Heuristically, \(R_\gamma(n)\) gains one period of \(\int s^{-\gamma}r_\infty\,ds=\int g\,d\log s\) per level. This is not proved here. **Falsifier:** a certified enclosure of \(\langle\Theta\rangle\) excluding \(c_R(\gamma)\pi/(2\sin\pi\gamma\ln5)\).

**OPEN.** Closed forms for \(c_R,c_S,c_F,C_\gamma\) and for \(\mathcal L\); a proof of Conjecture AT.2 (now given by Theorem AU); the sharp correction exponent (observed \(3^{-m}\), proved \(5^{-0.65m}\)); certified \(\mathcal L(\alpha)\) for \(0.65\le\alpha<\gamma\) (Theorem AS's power window).

**NOT THIS.** Theorem AT certifies the two refinement-rate constants of one graph resistance and one share on the combinatorial gasket (Theorem AR), to about 17 digits at five values of \(\alpha\) including \(0.45\) and \(\log3/\log5\). It also sharpens AR's \(o(1)\) and \(o(n)\) to exponential error terms and certifies the leading coefficient of the approach \(\|F(n,0.45)\|\to\mathcal L\). It gives no closed form (Kept Failure AT.1 rules out the simplest pure-power one at \(\gamma\)). It says nothing about the multi-scale or geometric-weight builds and does not touch Stage 1 or Stage 2; the Stage 2 divergence question stays OPEN. A rate constant of a graph audit norm is not a force: it is not thrust and not a pass of any experimental criterion. Not thrust.

## Theorem AU. Log-mean slope at \(\log3/\log5\) (Conjecture AT.2 proved) and the Mellin poles of the corner resolvent

**Setting.** As in Theorems AQ–AT: \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\), \(r_\infty=\lim r_n\), \(p(s)=s(s+5)\), \(q=p^{-1}\), \(\gamma=\log3/\log5\), \(g(s)=s^{1-\gamma}r_\infty(s)\), and \(\Theta(u)=\lim_kg(q^ku)\) (Kept Failure AT.1: the limit exists and \(\Theta\circ q=\Theta\)). Let \(\psi(s)=\lim_k5^kq^k(s)=s/A_\infty(s)\) be the Koenigs coordinate of \(q\) (\(A_\infty\) as in Theorem AT (2)), so \(\psi\circ q=\psi/5\). Put \(\theta=\Theta\circ\psi^{-1}\); then \(\theta(x/5)=\theta(x)\) exactly. Put \(\omega=2\pi/\ln5\) and define
\[
\langle\Theta\rangle=\frac1{\ln5}\int_{x/5}^{x}\theta(y)\,\frac{dy}{y},\qquad \hat\theta_k=\frac1{\ln5}\int_{x/5}^{x}\theta(y)\,y^{-i\omega k}\,\frac{dy}{y}=a_k-i\,b_k,
\]
with \(a_k,b_k\) the period means of \(\theta\cos(\omega k\ln y)\) and \(\theta\sin(\omega k\ln y)\). Both are independent of \(x\) because \(\theta\) is \(5\)-periodic. They are the precise form of AT.2's "mean of \(\Theta\) over one period in \(\log_5s\)": \(\langle\Theta\rangle=\lim_{u\to0}\frac1{\ln5}\int_{u/5}^{u}g\,\frac{ds}{s}\). **Assumptions:** Theorem AQ (1)–(4) (decimation; Stieltjes formula \(R_\alpha(m)=\frac{2\sin\pi\alpha}{\pi}\int_0^\infty s^{-\alpha}r_m\,ds\); \(r_m\uparrow r_\infty\), \(r_m\) decreasing in \(s\), \(r_m(0)=(5/3)^m/3\), \(r_m\le1/s\); the pointwise ceilings and the power ceiling \(r_m\le Bs^{-0.35}\) on \((0,0.01]\)); Theorem AR (2) (existence of \(c_R(\alpha)\) for every \(\alpha\in(0,1)\)); Theorem AS (2) (complete-monotone interpolation sandwich); Theorem AT (4a) at \(\gamma\) (or only AR's \(R_\gamma(n)=c_R(\gamma)n+o(n)\)); and Kept Failure AT.1's orbit lemma \(\log(g(y)/g(p(y)))\in[-(1-\gamma)y/5,\ 1.82\,y^{0.65}]\) for \(y\le10^{-3}\), all as stated. Items (1)–(2) are proofs. Item (3) is an outward-rounded `mpmath.iv` computation at 120 digits.

**THEOREM AU.**

1. **Slope formula (Conjecture AT.2 is PROVED).** \(\int_0^\infty s^{-\gamma}r_m(s)\,ds=m\ln5\,\langle\Theta\rangle+O(1)\), and therefore
   \[
   c_R(\gamma)=\frac{2\sin\pi\gamma}{\pi}\,\ln5\,\langle\Theta\rangle .
   \]
2. **Mellin poles (DERIVED).** \(M(z)=\int_0^\infty s^{z-1}r_\infty(s)\,ds\) is holomorphic on \(1-\gamma<\operatorname{Re}z<1\). It continues meromorphically to \(1-\gamma-0.65<\operatorname{Re}z<1\). The only possible poles there are at \(z_k=1-\gamma+i\omega k\), \(k\in\mathbb Z\); each is at most simple, with residue \(\hat\theta_{-k}=\overline{\hat\theta_k}\), and it is a genuine pole exactly when \(\hat\theta_k\ne0\). Since \(R_\infty(\alpha)=\frac{2\sin\pi\alpha}{\pi}M(1-\alpha)\) for \(\alpha<\gamma\),
   \[
   \lim_{\alpha\uparrow\gamma}(\gamma-\alpha)\,R_\infty(\alpha)=\frac{2\sin\pi\gamma}{\pi}\langle\Theta\rangle=\frac{c_R(\gamma)}{\ln5}.
   \]
3. **Certified values** (one period \(s\in[5^{-141},5^{-140}]\); \(128\) geometric cells; sandwich floor \(b=u=14\), ceiling \(b=u=15\)):
   \[
   \langle\Theta\rangle\in[0.6009542511687736078237051500999942,\ 0.6009542511687736078237051500999943],
   \]
   \[
   c_R(\gamma)\in[0.5171658346305052439259804856731271,\ 0.5171658346305052439259804856731273],\qquad \frac{c_R(\gamma)}{\ln5}\in[0.321333200016611999274880068557,\ 0.321333200016611999274880068558].
   \]
   This \(c_R(\gamma)\) enclosure has width \(<3\cdot10^{-34}\). It lies inside Theorem AT's independent enclosure \([0.517165834630505243,\,0.517165834630505245]\) (exact transfer at level \(40\); width \(2\cdot10^{-18}\)), so it sharpens that enclosure by 16 digits. The Fourier coefficients, each \(\pm2\cdot10^{-35}\), are:

| \(k\) | \(a_k\) | \(b_k\) | \(\lvert\hat\theta_k\rvert\ge\) |
|---|---|---|---|
| 1 | \(4.042278096762267394377\cdot10^{-8}\) | \(-3.663477703263059920737\cdot10^{-6}\) | \(3.66370\cdot10^{-6}\) |
| 2 | \(-1.347221801190789742277\cdot10^{-11}\) | \(6.374978894928635265275\cdot10^{-13}\) | \(1.34873\cdot10^{-11}\) |
| 3 | \(-1.664604061435855694403\cdot10^{-17}\) | \(8.666201527432484696846\cdot10^{-17}\) | \(8.82462\cdot10^{-17}\) |
| 4 | \(3.253673674396170403376\cdot10^{-22}\) | \(2.762525729016318114163\cdot10^{-22}\) | \(4.26825\cdot10^{-22}\) |

   So \(M\) has genuine simple poles at \(1-\gamma\pm i\omega k\) for \(k=1,2,3,4\), in addition to the real pole at \(1-\gamma\) with residue \(\langle\Theta\rangle\).

**Proof.** *(1)* Write \(J_m=\int_0^\infty s^{-\gamma}r_m\,ds\) and split \((0,\infty)\) at \(5^{-m}\) and \(1\):
\[
J_m=\underbrace{\int_0^{5^{-m}}\!s^{-\gamma}r_m}_{\rm(i)}+\underbrace{\int_{5^{-m}}^1\!s^{-\gamma}r_\infty}_{\rm(ii)}-\underbrace{\int_{5^{-m}}^1\!s^{-\gamma}(r_\infty-r_m)}_{\rm(iii)}+\underbrace{\int_1^\infty\!s^{-\gamma}r_m}_{\rm(iv)}.
\]
(i) Since \(r_m\le r_m(0)=(5/3)^m/3\), the term lies in \([0,\,(5/3)^m5^{-m(1-\gamma)}/(3(1-\gamma))]=[0,\,1/(3(1-\gamma))]\), because \((5/3)5^{\gamma-1}=5^\gamma/3=1\). (iv) Since \(r_m\le1/s\), it lies in \([0,1/\gamma]\). (iii) Fix \(\varepsilon\in(0,\gamma)\) and put \(\alpha'=\gamma-\varepsilon\), \(\rho'=5^{\alpha'}/3=5^{-\varepsilon}\). For \(s\ge5^{-m}\) we have \(s^{-\gamma}\le5^{m\varepsilon}s^{-\alpha'}\). By monotone convergence and Theorem AR (2) at \(\alpha'\) (the convergent sequence \(\rho'^{-j}\int s^{-\alpha'}\delta_{j+1}\) is bounded by some \(K\)),
\[
0\le{\rm(iii)}\le5^{m\varepsilon}\sum_{j\ge m}\int_0^\infty s^{-\alpha'}\delta_{j+1}\le5^{m\varepsilon}K\sum_{j\ge m}5^{-j\varepsilon}=\frac{K}{1-5^{-\varepsilon}} .
\]
(ii) This equals \(\int_{5^{-m}}^1g\,ds/s\). Fix \(s_0\le10^{-3}\). The orbit lemma gives \(\Theta(s)\in g(s)\exp([-(1-\gamma)q(s)/4,\ 2q(s)^{0.65}/(1-5^{-0.65})])\), so \(|g-\Theta|\le Cs^{0.65}\) on \((0,s_0]\) and \(g,\Theta\) are bounded there. Hence \(\int_{5^{-m}}^{s_0}g\,ds/s=\int_{5^{-m}}^{s_0}\Theta\,ds/s+O(1)\). Now substitute \(x=\psi(s)\). We have \(\Theta(s)=\theta(x)\), and \(d\log s=d\log x/(1-\epsilon(s))\) with \(0\le\epsilon(s)=s(\log A_\infty)'(s)\le s/20\), because \(q(x)\le x/5\) and \(q'\le1/5\) give \((q^i)'\le5^{-i}\). So the integral is \(\int_{\psi(5^{-m})}^{\psi(s_0)}\theta\,dx/x+O(1)\). A \(5\)-periodic bounded \(\theta\) integrates to \(\ln5\,\langle\Theta\rangle\) over each full period and to \(O(1)\) over a partial one. Also \(\log\psi(5^{-m})=-m\ln5-\log A_\infty(5^{-m})=-m\ln5+O(5^{-m})\), since \(0\le\log A_\infty(s)\le s/20\). Hence (ii) \(=m\ln5\,\langle\Theta\rangle+O(1)\), and \(J_m=m\ln5\,\langle\Theta\rangle+O(1)\). By AQ (2), \(R_\gamma(m)=\frac{2\sin\pi\gamma}{\pi}J_m\). By AT (4a) (AR's \(o(m)\) is enough), \(R_\gamma(m)=c_R(\gamma)m+O(1)\). Divide by \(m\) and let \(m\to\infty\).

*(2)* On \([s_0,\infty)\), \(\int s^{z-1}r_\infty\) is holomorphic for \(\operatorname{Re}z<1\) (because \(r_\infty\le1/s\)). On \((0,s_0]\) put \(w=z-(1-\gamma)\), so the integral is \(\int_0^{s_0}s^wg\,ds/s\). Write \(g=\Theta+e\) with \(|e|\le Cs^{0.65}\); the \(e\)-part is holomorphic for \(\operatorname{Re}w>-0.65\). In the \(\Theta\)-part substitute \(x=\psi(s)\), \(s=xa\) with \(a=A_\infty(s)\in[1,e^{s/20}]\). The factor \(a^w/(1-\epsilon)-1\) is \(O(|w|x)\), so that error is holomorphic for \(\operatorname{Re}w>-1\). The main part, with \(X=\psi(s_0)\), is
\[
\int_0^X x^w\theta(x)\frac{dx}{x}=\sum_{j\ge0}5^{-jw}\!\int_{X/5}^{X}y^w\theta(y)\frac{dy}{y}=\frac{F(w)}{1-5^{-w}},
\]
with \(F\) entire. Its poles lie where \(5^{-w}=1\), i.e. \(w=i\omega k\), with residue \(F(i\omega k)/\ln5=\hat\theta_{-k}\). The real limit is the \(k=0\) residue: \(z=1-\alpha\), \(w=\gamma-\alpha\). Continuity of \(\Theta\) (a uniform limit on a fundamental domain, by the orbit lemma) and \(\theta>0\) are all that is used.

*(3)* All arithmetic is outward-rounded. At \(129+30\) geometric nodes \(5^{-141}5^{k/128}\), AQ (4) encloses \(r_\infty\) by \(T^k r_0\le r_\infty\le T^k(1/s)\), iterated until the orbit passes \(10^{90}\) (relative width \(<10^{-60}\)). \(r_\infty\) is completely monotone (AS (1)). So on each cell the Lagrange interpolant with \(14\) nodes on each side is a floor and the one with \(15\) is a ceiling (AS (2)). Integrating these against \(s^{-\gamma}\) gives \(I_0=\int_{5^{-141}}^{5^{-140}}s^{-\gamma}r_\infty\). Integrating against the nonnegative weights \(s^{-\gamma}(1\pm\cos(\omega k\ln s))\) and \(s^{-\gamma}(1\pm\sin(\omega k\ln s))\) (closed-form \(\int_1^Qy^{e-1+i\tau}dy\) Lagrange weights) and subtracting \(I_0\) gives the trigonometric moments; the two signs give two enclosures, which are intersected. The corrections to \(\langle\Theta\rangle=\frac1{\ln5}\int_{q(u)}^{u}\Theta\,d\log\psi\) (\(u=5^{-140}\)) are: \(\Theta/g\in[e^{-(1-\gamma)q(u)/4},e^{2q(u)^{0.65}/(1-5^{-0.65})}]\) (\(\le e^{2.7\cdot10^{-64}}\)); \(d\log\psi/d\log s\in[1-u/20,1]\); the phase \(|\omega k(\log\psi-\log s)|\le\omega ku/20\); and the end piece \(\int_{q(u)}^{u/5}g\,ds/s\le(u/5)^{1-\gamma}r_\infty(q(u))\,u/25\) (\(\le3.4\cdot10^{-100}\)). The resulting enclosures are the ones displayed. \(\square\)

**Witness, computed by `scripts/gasket_fractional_log_mean.py`** (about 20 s).
- AQ's power ceiling re-certified (\(B\approx0.892577\)); the AT.1 orbit lemma (ratio \(1.7896\le1.82\)); Koenigs bounds \(q(x)\le x/5\), \(q'\le1/5\).
- The table and enclosures above. A second period, shifted (\(s\in5^{-121}5^{0.37}[1/5,1]\)), gives an overlapping \(\langle\Theta\rangle\), which checks that the period mean does not depend on the base point.
- Agreement: \(\frac{2\sin\pi\gamma}{\pi}\ln5\,\langle\Theta\rangle\) lies inside Theorem AT's independently certified \(c_R(\gamma)\) interval (two different computations: AT transfers \(\delta_{41}\) at level \(40\); AU integrates \(r_\infty\) one period at depth \(5^{-140}\)).
- The truncated Fourier series \(\langle\Theta\rangle+2\sum_{k\le4}(a_k\cos+b_k\sin)\) reproduces Kept Failure AT.1's eight recorded 9-digit \(\Theta\) values at \(5^{-40}5^{j/8}\) to \(4\cdot10^{-10}\).
- **OBSERVATION, not a theorem:** \(J_m-m\ln5\,\langle\Theta\rangle=1.318164,\ 1.314889,\ 1.31487559740,\ 1.314875597174,\ 1.314875597174\) at \(m=5,10,20,30,40\) (30-digit quadrature). This is bounded, as (1) requires, and it converges. So \(C_\gamma=\lim(R_\gamma(m)-c_R(\gamma)m)\approx0.70307046242184\) (float, not certified). Successive \(|\hat\theta_{k+1}/\hat\theta_k|\) are \(3.7\cdot10^{-6},\,6.5\cdot10^{-6},\,4.8\cdot10^{-6}\), the order of \(e^{-2\pi^2/\ln5}\approx4.7\cdot10^{-6}\). This is not proved.

The script exits nonzero if any check fails.

**KEPT FAILURE AU.1 (the oscillation is not a single harmonic).** The claim "\(\theta(x)=\langle\Theta\rangle+2\,\mathrm{Re}(\hat\theta_1x^{i\omega})\)", a pure cosine in \(\log_5x\), is false. Kept Failure AT.1's eight values fit such a cosine to the recorded 9 digits, which is what makes the claim tempting. But \(|\hat\theta_2|\ge1.34873\cdot10^{-11}>0\) is certified, so the \(k=2\) term is present, and by (2) \(M\) has a genuine pole at \(1-\gamma\pm2i\omega\) (likewise for \(k=3,4\)). Witness: the same script.

**OPEN.** Closed forms for \(\langle\Theta\rangle\), \(\hat\theta_k\), \(c_R\), \(c_S\), \(C_\gamma\), \(\mathcal L\); a certified \(C_\gamma\) (now given by Theorem AY); whether \(\hat\theta_k\ne0\) for every \(k\), and the exact decay rate of \(|\hat\theta_k|\); the analogous log-periodic structure for \(1-s\,r_\infty\) (the \(c_S\) side); the sharp correction exponent of Theorem AT (4a).

**NOT THIS.** Theorem AU identifies the slope of one graph resistance at \(\alpha=\log3/\log5\) on the combinatorial gasket as the log-mean of the self-similar profile of the corner resolvent. It thereby proves Conjecture AT.2 and certifies that slope to 34 digits, along with the first four Fourier modes of the profile (the complex Mellin poles). It gives no closed form. It says nothing about the multi-scale or geometric-weight builds. It does not touch Stage 1 or Stage 2; the Stage 2 divergence question stays OPEN. A Mellin pole of a graph resolvent is not a force: it is not thrust and not a pass of any experimental criterion. Not thrust.

## Theorem AY. Certified intercept \(C_\gamma=0.70307046242185742992904816\ldots\) and the level shift of the endpoint audit norm

**Setting.** As in Theorems AQ–AU, at \(\alpha=\gamma=\log3/\log5\), where \(\rho=5^{\gamma}/3=1\): \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\), \(r_\infty=\lim r_n\), \(\delta_{j+1}=r_{j+1}-r_j>0\), \(p(s)=s(s+5)\), \(q=p^{-1}\), \(\kappa=2\sin\pi\gamma/\pi\). Put \(J_n=\int_0^\infty s^{-\gamma}r_n\,ds\), so \(R_\gamma(n)=\kappa J_n\) (AQ (2)), and \(c=\ln5\,\langle\Theta\rangle\), so \(c_R(\gamma)=\kappa c\) (AU (1)). The bracket weights of Theorem AT (2) do not depend on the level: \(B_\infty=\prod_{i\ge0}5/\sqrt{25+4q^i(t)}\le1\), \(A_\infty=\prod_{i\ge1}(1+q^i(t)/5)\ge1\), and \(\bar C=\prod_{i\ge1}\tfrac35\bar D(q^i(t))\ge1\). Put \(W_+=A_\infty^{\gamma}\bar C-1\ge0\) and \(W_-=1-B_\infty\ge0\). Theorem AT (4a) defined \(C_\gamma=\lim_n(R_\gamma(n)-c_R(\gamma)n)\) and Theorem AU recorded the float value \(\approx0.70307046242184\), but neither certified it. **Assumptions:** Theorem AQ (1)–(4) (decimation; Stieltjes formulas for \(R_\alpha(n)\) and \(\sigma_n\); \(r_n\uparrow\), \(r_n\) decreasing in \(s\), \(r_n(0)=(5/3)^n/3\); the node enclosures \(T^kr_0\le r_\infty\le T^k(1/s)\) and the power ceiling \(r_m\le Bs^{-0.35}\) on \((0,0.01]\)); Theorem AR (1)–(3) (product formula, existence of \(c_R,c_S\), \(\sigma_\infty\ge2\cdot4^{-\alpha}\)); Theorem AS (2)–(3) (complete monotonicity of \(r_n,r_\infty,1-sr_\infty\), the interpolation sandwich, and the tail moments \(r\in[1/s-2/s^2,1/s]\), \(1-sr_\infty\in[2/s-6/s^2,2/s]\)); Theorem AT (2)–(4a) (the level-\(n\) bracket, the small-\(y\) lemma \(\log\tfrac35\bar D(y)\le1.82\,y^{0.65}\) for \(y\le10^{-3}\), the crude tail \(\delta_{n+1}(s)\le2\cdot1.2^n/(s\,p^n(s))\) and \(A_\infty^{\alpha}\bar C\le3\cdot10^{10}s^5\) for \(s\ge1\), and the rate \(O(5^{-0.65n})\)); Theorem AU (1) and (3) (\(c_R(\gamma)=\kappa\ln5\,\langle\Theta\rangle\), the certified \(\langle\Theta\rangle\), and \(0\le\log A_\infty\le s/20\)); all as stated. Items (1) and (3) are proofs. Item (2) is an outward-rounded `mpmath.iv` computation at 80 digits.

**THEOREM AY.**

1. **Summable defect.** The limit \(C'=\lim_m(J_m-mc)\) exists and
   \[
   J_n-nc-C'\in[-E_n^-,\,E_n^+],\qquad E_n^\pm=\int_0^\infty t^{-\gamma}\bigl(r_\infty-r_n\bigr)(t)\,W_\pm(t)\,dt .
   \]
   holds for every \(n\ge0\) (trivially where \(E_n^\pm=\infty\)), and \(E_m^\pm\le E_n^\pm\) for \(m\ge n\). Hence \(C_\gamma=\kappa C'\), and for every \(m\ge n\),
   \[
   R_\gamma(m)-c_R(\gamma)\,m-C_\gamma\in\kappa\,[-E_n^-,\,E_n^+].
   \]
2. **Certified values** (geometric grid of ratio \(65/64\) from \(5^{-60}2^{-110}\) to \(10^{25}\), \(14891\) nodes; sandwich floor \(b=u=14\), ceiling \(b=u=15\); level \(n=60\), cross-checked at \(n=40\)):
   \[
   E_{60}^+\le8.69\cdot10^{-28},\quad E_{60}^-\le5.9\cdot10^{-34},\quad E_{40}^+\le1.06\cdot10^{-18},
   \]
   \[
   C_\gamma\in[0.70307046242185742992904816,\ 0.70307046242185742992904817],
   \]
   \[
   \nu_\gamma=\frac{C_\gamma}{c_R(\gamma)}\in[1.35946811514371162883271469,\ 1.35946811514371162883271470],
   \]
   \[
   \sigma_\infty(\gamma)=\tfrac{\sin\pi\gamma}{\pi}\!\int_0^\infty\! s^{-\gamma}(1-s\,r_\infty)\,ds\in[1.15191861167405417640547243,\ 1.15191861167405417640547244],
   \]
   \[
   K_\gamma=\frac{\sqrt{21}\,\sigma_\infty(\gamma)}{3\,c_R(\gamma)}\in[3.40236076383265507999048519,\ 3.40236076383265507999048520].
   \]
   So \(|R_\gamma(m)-c_R(\gamma)m-C_\gamma|\le4.7\cdot10^{-28}\) for every \(m\ge60\) and \(\le5.7\cdot10^{-19}\) for every \(m\ge40\). Theorem AU's float observation \(C_\gamma\approx0.70307046242184\) was low by \(1.7\cdot10^{-14}\): its 14th digit is wrong.
3. **Level-shift law (DERIVED).** The endpoint audit norm is a shifted hyperbola in the level up to an exponentially small error:
   \[
   (n+\nu_\gamma)\,\|F(n,\gamma)\|=K_\gamma+O(5^{-0.65n}).
   \]

**Proof.** *(1)* Theorem AT (2) at \(\alpha=\gamma\), \(\rho=1\), with the level renamed \(j\): \(c\in[\int t^{-\gamma}\delta_{j+1}B_\infty,\ \int t^{-\gamma}\delta_{j+1}A_\infty^{\gamma}\bar C]\). The \(k=0\) term of AT (1) is \(\Delta_j=J_{j+1}-J_j=\int t^{-\gamma}\delta_{j+1}\), with weight \(A_0^{\gamma}B_0C_0=1\). For every \(k\), \(A_k^{\gamma}B_kC_k-1\in[B_\infty-1,\ A_\infty^{\gamma}\bar C-1]\), because \(A_k,C_k\ge1\ge B_k\ge B_\infty\). Let \(k\to\infty\) (AR (2)):
\[
-\int t^{-\gamma}\delta_{j+1}W_-\ \le\ c-\Delta_j\ \le\ \int t^{-\gamma}\delta_{j+1}W_+ .
\]
The weights do not depend on \(j\). Sum over \(j\ge n\). By Tonelli (all terms \(\ge0\)) and \(\sum_{j\ge n}\delta_{j+1}=r_\infty-r_n\) (AQ (3)), \(\sum_{j\ge n}|\Delta_j-c|\le E_n^++E_n^-\). Here \(E_n^\pm\in[0,\infty]\); once \(E_n^\pm<\infty\) for one \(n\) (shown in (2) for \(n=40\) and \(n=60\)), \(J_m-mc=J_n-nc+\sum_{j=n}^{m-1}(\Delta_j-c)\) converges, and its limit \(C'\) satisfies the displayed inclusion. \(r_m\) increases in \(m\), so \(r_\infty-r_m\) and hence \(E_m^\pm\) decrease. Multiply by \(\kappa\).

*(2)* All arithmetic is outward-rounded. *Cells* \([x_c,x_{c+1}]\subset(x_0,10^{25}]\): Theorem AS (2) gives a floor (\(b=u=14\)) and a ceiling (\(b=u=15\)) of \(\int_{\rm cell}s^{-\gamma}f\) for \(f\in\{r_{60},r_{40},r_\infty,1-sr_\infty\}\). Node values of \(r_n\) come from the exact recursion of AQ (1) in interval arithmetic. Those of \(r_\infty\) come from AQ (4), iterated until the orbit passes \(10^{90}\) (relative width \(<10^{-60}\)). On each cell in \((x_0,1]\), \(0\le\int_{\rm cell}t^{-\gamma}(r_\infty-r_n)\le\mathrm{ceil}(r_\infty)-\mathrm{floor}(r_n)\), multiplied by \(\sup W_\pm\) over a block of \(16\) cells. That supremum uses \(\sup A_\infty=A_\infty(x_{\rm hi})\), \(\inf B_\infty=B_\infty(x_{\rm hi})\), and \(\bar C\) bounded over the block, all by Theorem AT's weight routine. *Head* \((0,x_0]\), \(x_0=5^{-60}2^{-110}\): \(r_\infty-r_n\le r_\infty\le Bt^{-0.35}\). By AU, \(\log A_\infty\le t/20\). By the AT small-\(y\) lemma (with \(q^i(t)\le5^{-i}t\)), \(\log\bar C\le1.82\,t^{0.65}\sum_{i\ge1}5^{-0.65i}\). Together these give \(W_+\le K_ht^{0.65}\) there, and \(W_-\le\sum_i2q^i(t)/25\le t/10\). So the head of \(E^+\) is \(\le BK_hx_0^{1.3-\gamma}/(1.3-\gamma)\), and that of \(E^-\) is smaller. *Tail* \([1,\infty)\): \(W_\pm\le A_\infty^{\gamma}\bar C\le3\cdot10^{10}t^5\). Also \(\delta_{j+1}(t)\le2\cdot1.2^j/(t\,p^j(t))\) and \(p^j(t)\ge(6t)^{2^{j-1}}\), because \(p(t)\ge6t\) and \(p(t)\ge t^2\) for \(t\ge1\). So \(\sum_{j\ge n}\int_1^\infty\le2\cdot6\cdot10^{10}\cdot1.2^n6^{-2^{n-1}}/(2^{n-1}-5)\), which is below \(10^{-4\cdot10^{17}}\) at \(n=60\). *\(J_n\)*: cells as above, plus a head \(r_n(x_0)\le r_n\le r_n(0)=(5/3)^n/3\) and a tail \(1/s-2/s^2\le r_n\le1/s\). *\(\sigma_\infty(\gamma)\)*: AS (1) with \(\beta=1-\gamma\), so the weight is \(s^{-\gamma}\). The same cells apply to \(1-sr_\infty\), with a head \(1-x_0r_\infty(x_0)\le1-sr_\infty\le1\) and a tail \([2/s-6/s^2,2/s]\). Finally \(c=\ln5\,\langle\Theta\rangle\) is taken from AU (3) (width \(<10^{-33}\)), and \(\kappa c\) lies inside AT's independent \(c_R(\gamma)\) enclosure.

*(3)* \(\|F(n,\alpha)\|=\tfrac{\sqrt{21}}3\sigma_n/R_\alpha(n)\) (AQ). At \(\gamma\), \(\mu=5^{\gamma-1}/3=1/5\), so \(\sigma_n-\sigma_\infty=\sum_{j\ge n}(\sigma_j-\sigma_{j+1})=O(5^{-n})\) by AR (2). By (1) and AT (4a), \(R_\gamma(n)=c_R(\gamma)(n+\nu_\gamma)+O(5^{-0.65n})\). Then \((n+\nu_\gamma)\|F\|=\tfrac{\sqrt{21}}{3c_R}\sigma_n\bigl(1+O(5^{-0.65n}/n)\bigr)^{-1}=K_\gamma+O(5^{-0.65n})\). \(\square\)

**Witness, computed by `scripts/gasket_fractional_endpoint_shift.py`** (about 3–4 min).
- Re-certified: AQ's power ceiling (\(B\approx0.892577\)), the AT small-\(y\) lemma (ratio \(1.7896\le1.82\)), and \(\tfrac12\log(1+4y/25)\le2y/25\). The two cited inputs agree: \(\kappa\ln5\,\langle\Theta\rangle\) (AU) lies inside AT's \(c_R(\gamma)\).
- Node enclosures: \(r_\infty\) has relative width \(<10^{-60}\), and \(r_{40}\le r_{60}\le r_\infty\) holds at every node. The weights satisfy \(W_+\le0.938\) and \(W_-\le0.0888\) on \((0,1]\).
- \(J_{60}\in[59.346788925341940900804055319079,\,\ldots080]\) and \(J_{40}\in[40.002817815952560490333706099377,\,\ldots378]\), each of width \(<10^{-30}\). The levels-\(60\) and \(40\) enclosures of \(C'\) overlap. \(E_{40}\) and \(E_{60}\) are both below \(5^{-0.6n}\). The enclosures above.
- Consistency (30-digit quadrature, not part of the proof): \(J_m-mc-C'=3.29\cdot10^{-3},\,1.34\cdot10^{-5},\,2.27\cdot10^{-10},\,3.36\cdot10^{-15}\) at \(m=5,10,20,30\) (positive, inside \(5\cdot5^{-0.65m}\)). \((m+\nu_\gamma)\|F(m,\gamma)\|-K_\gamma=-1.26\cdot10^{-3},\,-3.98\cdot10^{-6},\,-3.74\cdot10^{-11},\,-3.76\cdot10^{-16}\). **OBSERVATION, not a theorem:** the ratio of \(J_m-mc-C'\) over ten levels is \(1.6916\cdot10^{-5}\approx3^{-10}\). This is the same \(3^{-m}\) correction that Theorem AT observed, faster than the proved \(5^{-0.65m}\).

The script exits nonzero if any check fails.

**KEPT FAILURE AY.1 (the endpoint needs a level shift).** The claim "at \(\alpha=\gamma\), as below it, the renormalized audit norm converges exponentially fast with no shift: \(n\|F(n,\gamma)\|=K_\gamma+O(\theta^n)\) for some \(\theta<1\)" is false. Below \(\gamma\), \(\|F(n,\alpha)\|-\mathcal L=O(\rho^n)\) (AR, AT (4b)); the tempting endpoint analogue is "\(1/n\) times a constant, up to exponentially small error". By (3), \(n\|F(n,\gamma)\|-K_\gamma=-\nu_\gamma K_\gamma/(n+\nu_\gamma)+O(5^{-0.65n})\), and \(\nu_\gamma K_\gamma\ge4.625400\) is certified. So \(n\,(n\|F(n,\gamma)\|-K_\gamma)\to-\nu_\gamma K_\gamma\neq0\), and the error is \(\Theta(1/n)\), not exponential. Witness (same script, 30-digit quadrature): \(n\|F\|-K_\gamma=-0.728312,\,-0.407188,\,-0.216550,\,-0.147496\) at \(n=5,10,20,30\); \(n(n\|F\|-K_\gamma)=-4.4249\) at \(n=30\), which matches \(-\nu_\gamma K_\gamma\,n/(n+\nu_\gamma)\) to \(10^{-10}\). The renormalization that converges exponentially is \(n+\nu_\gamma\), with \(\nu_\gamma=1.359468\ldots\), not \(n\).

**OPEN.** A closed form for \(C_\gamma\), \(\nu_\gamma\), \(\sigma_\infty(\gamma)\), \(K_\gamma\); a proof of the observed \(3^{-n}\) correction, which would also improve the error in (3); and a distributional or measure-valued limit of the level-\(n\) signed corner flux on the infinite gasket beyond this one norm. Theorem AY renormalizes one scalar audit norm and does not construct a flux functional.

**NOT THIS.** Theorem AY certifies to about 27 digits the intercept of the linearly growing graph resistance \(R_\gamma(n)\) at the single exponent \(\alpha=\log3/\log5\) on the combinatorial gasket, together with the limit share \(\sigma_\infty(\gamma)\). It shows that the vanishing audit norm there is \(K_\gamma/(n+\nu_\gamma)\) up to \(O(5^{-0.65n})\). It gives no closed form. It does not construct a flux on the infinite gasket. It says nothing about \(\alpha\ne\gamma\), the multi-scale or geometric-weight builds, Stage 1, or Stage 2; the Stage 2 divergence question stays OPEN. A level shift of a decaying graph audit norm is not a force: it is not thrust and not a pass of any experimental criterion. Not thrust.

## Theorem AZ. Supercritical decay law: \(\|F(n,\alpha)\|\sim K_\alpha(3\cdot5^{-\alpha})^n\), the intercept is the continued Mellin transform, \(C_{0.8}=-1.74325717839207572395\ldots\)

**Setting.** As in Theorems AQ–AY, now with \(\alpha\in(\gamma,1)\), \(\gamma=\log3/\log5\), so \(\rho=5^{\alpha}/3>1\) (Theorem K's side): \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\), \(r_\infty=\lim r_n\), \(\delta_{j+1}=r_{j+1}-r_j>0\), \(p(s)=s(s+5)\), \(q=p^{-1}\), \(\kappa=2\sin\pi\alpha/\pi\), \(\mu=\rho/5\). Put \(J_n=\int_0^\infty s^{-\alpha}r_n\,ds\), so \(R_\alpha(n)=\kappa J_n\) (AQ (2)), \(\Delta_j=J_{j+1}-J_j\), and \(c=c_R(\alpha)/\kappa=\lim_j\rho^{-j}\Delta_j\) (AR (2)). \(M(z)=\int_0^\infty s^{z-1}r_\infty\,ds\) is AU (2)'s Mellin transform, holomorphic for \(1-\gamma<\operatorname{Re}z<1\), with \(R_\infty(\alpha)=\kappa M(1-\alpha)\) for \(\alpha<\gamma\). For complex \(z\) write \(\rho_z=5^{1-z}/3\) (so \(\rho_{1-\alpha}=\rho\)). The bracket weights of Theorem AT (2) do not depend on the level: \(W_+=A_\infty^{\alpha}\bar C-1\ge0\), \(W_-=1-B_\infty\ge0\). Theorem AR (4c) gave only \(\|F(n,\alpha)\|=\Theta((3\cdot5^{-\alpha})^n)\), and Theorem K only the upper bound \(C(\alpha)(3\cdot5^{-\alpha})^{n/2}\). **Assumptions:** Theorem AQ (1)–(4) (decimation; Stieltjes formulas; \(r_n\uparrow\), \(r_n\) decreasing in \(s\), \(r_n(0)=(5/3)^n/3\); node enclosures \(T^kr_0\le r_\infty\le T^k(1/s)\); power ceiling \(r_m\le Bs^{-0.35}\) on \((0,0.01]\)); Theorem AR (1)–(3); Theorem AS (1)–(3) (complete monotonicity, interpolation sandwich, tail moments); Theorem AT (1)–(4a) (exact transfer, level-\(n\) bracket, small-\(y\) lemma, crude tail \(A_\infty^{\alpha}\bar C\le3\cdot10^{10}s^5\) on \(s\ge1\), rates \(O(5^{-0.65n})\)); Theorem AU (2)–(3) (meromorphic continuation, residue \(\langle\Theta\rangle>0\) at \(1-\gamma\)); Theorem AY (2) (\(C_\gamma\)); all as stated. Items (1), (2), (4) are proofs. Item (3) is an outward-rounded `mpmath.iv` computation at 80 digits.

**THEOREM AZ.**

1. **Geometric expansion with a summable defect.** For every \(n\ge0\),
   \[
   -\int t^{-\alpha}\delta_{j+1}W_-\ \le\ c\rho^j-\Delta_j\ \le\ \int t^{-\alpha}\delta_{j+1}W_+\qquad(j\ge0),
   \]
   so \(C'=\lim_m\bigl(J_m-c\rho^m/(\rho-1)\bigr)\) exists and
   \[
   J_n-\frac{c\,\rho^n}{\rho-1}-C'\in[-E_n^-,\,E_n^+],\qquad E_n^\pm=\int_0^\infty t^{-\alpha}(r_\infty-r_n)\,W_\pm\,dt ,
   \]
   with \(E_m^\pm\le E_n^\pm\) for \(m\ge n\). Equivalently, with \(A_\alpha=c_R(\alpha)/(\rho-1)\) and \(C_\alpha=\kappa C'\),
   \[
   R_\alpha(m)=A_\alpha\,\rho^m+C_\alpha+\varepsilon_m,\qquad \varepsilon_m\in\kappa[-E_n^-,E_n^+]\ \ (m\ge n),\qquad \varepsilon_m=O(\rho^m5^{-0.65m}).
   \]
2. **The intercept is the continued Mellin transform (DERIVED).** For \(0<\operatorname{Re}z<1\) and every \(n\),
   \[
   M(z)=I_n(z)+\frac{c(z)\,\rho_z^{\,n}}{1-\rho_z}+\sum_{j\ge n}e_j(z),\qquad I_n(z)=\int_0^\infty s^{z-1}r_n\,ds ,
   \]
   where \(c(z)=\lim_j\rho_z^{-j}\int s^{z-1}\delta_{j+1}\) is holomorphic, \(e_j(z)=\int s^{z-1}\delta_{j+1}-c(z)\rho_z^j\), and \(\sum_{j\ge n}e_j\) converges locally uniformly. The right side is the meromorphic continuation of \(M\) to \(0<\operatorname{Re}z<1\); its only possible poles are \(\rho_z=1\), i.e. \(z_k=1-\gamma+2\pi ik/\ln5\), with residue \(c(z_k)/\ln5\) (at \(k=0\): \(c_R(\gamma)/(\kappa\ln5)=\langle\Theta\rangle\), as in AU). At \(z=1-\alpha\) it reads
   \[
   C_\alpha=\frac{2\sin\pi\alpha}{\pi}\,M(1-\alpha)\qquad(\gamma<\alpha<1):
   \]
   the supercritical intercept is the analytic continuation of the limit resistance \(R_\infty(\alpha)=\kappa M(1-\alpha)\) across its simple pole at \(\alpha=\gamma\). Hence \((\alpha-\gamma)C_\alpha\to-c_R(\gamma)/\ln5\) as \(\alpha\downarrow\gamma\) (so \(C_\alpha<0\) just above \(\gamma\)), and the regular part is continuous through the endpoint: \(C_\alpha+c_R(\alpha)/(\rho-1)\to C_\gamma\) (Theorem AY) as \(\alpha\to\gamma\) from either side (for \(\alpha<\gamma\), \(C_\alpha:=R_\infty(\alpha)\)).
3. **Certified values at \(\alpha=0.8\)** (Theorem AY's grid: ratio \(65/64\) from \(5^{-60}2^{-110}\) to \(10^{25}\), \(14891\) nodes, sandwich floor \(b=u=14\), ceiling \(b=u=15\); level \(n=60\), cross-checked at \(n=40\); \(\rho=1.20796610612949255\ldots\), \(3\cdot5^{-0.8}=0.82783779687672890\ldots\)):
   \[
   E_{60}^+\le7.22\cdot10^{-23},\quad E_{60}^-\le6.7\cdot10^{-34},\quad E_{40}^+\le2.01\cdot10^{-15},
   \]
   \[
   c_R(0.8)\in[0.4983719948274884072103,\,\ldots104],\qquad A_{0.8}=\frac{c_R}{\rho-1}\in[2.3964097039792205270396,\,\ldots397],
   \]
   \[
   C_{0.8}\in[-1.7432571783920757239548,\,-1.7432571783920757239546],\qquad \sigma_\infty(0.8)\in[1.0791380336973826964124,\,\ldots125],
   \]
   \[
   \nu_{0.8}=\frac{C_{0.8}}{A_{0.8}}\in[-0.7274453844421553283991,\,-0.7274453844421553283990],\qquad
   K_{0.8}=\frac{\sqrt{21}\,\sigma_\infty}{3A_{0.8}}\in[0.6878667584139920685497,\,\ldots498].
   \]
   So \(|R_{0.8}(m)-A_{0.8}\rho^m-C_{0.8}|\le2.8\cdot10^{-23}\) for every \(m\ge60\). At \(\alpha=0.75\) (level \(60\)): \(c_R(0.75)\in[0.50798512812307318007,\,\ldots008]\), inside Theorem AT's independent \([0.507985128123073180,\,\ldots181]\); \(A_{0.75}\in[4.433950023619875926,\,\ldots927]\), \(C_{0.75}\in[-3.764208555814372066,\,-3.764208555814372065]\), \(\sigma_\infty(0.75)\in[1.107854587969387386,\,\ldots387]\), \(K_{0.75}\in[0.381663263480592128,\,\ldots129]\).
4. **Decay law (DERIVED).** For every \(\alpha\in(\gamma,1)\),
   \[
   (\rho^n+\nu_\alpha)\,\|F(n,\alpha)\|=K_\alpha+O(5^{-0.65n}),\qquad K_\alpha=\frac{\sqrt{21}\,\sigma_\infty(\alpha)\,(\rho-1)}{3\,c_R(\alpha)},\quad \nu_\alpha=\frac{C_\alpha}{A_\alpha},
   \]
   so \(\|F(n,\alpha)\|=K_\alpha\,(3\cdot5^{-\alpha})^n\bigl(1-\nu_\alpha(3\cdot5^{-\alpha})^n+O((3\cdot5^{-\alpha})^n5^{-0.65n})\bigr)\). The exact decay ratio is \(1/\rho=3\cdot5^{-\alpha}\), the square of Theorem K's majorant ratio \((3\cdot5^{-\alpha})^{1/2}\). The leading term is a constant, not log-periodic: along integer levels \(\rho^n\) is sampled exactly, and the complex Mellin poles \(z_k\) (\(k\ne0\)) of AU sit at \(\rho_{z_k}=1\), never at the real point \(z=1-\alpha\).

**Proof.** *(1)* Theorem AT (2) at level \(j\), with \(k\to\infty\) in AT (1): \(c\rho^j=\lim_k\int t^{-\alpha}\delta_{j+1}A_k^{\alpha}B_kC_k\), and the \(k=0\) term is \(\Delta_j\) (weight \(1\)). For every \(k\), \(A_k^{\alpha}B_kC_k-1\in[B_\infty-1,\,A_\infty^{\alpha}\bar C-1]\) because \(A_k,C_k\ge1\ge B_k\ge B_\infty\), \(A_k\le A_\infty\), \(C_k\le\bar C\). This gives the displayed per-level bracket with the level-free weights. Sum over \(j\ge n\): by Tonelli and \(\sum_{j\ge n}\delta_{j+1}=r_\infty-r_n\), \(\sum_{j\ge n}|\Delta_j-c\rho^j|\le E_n^++E_n^-\). Since \(J_m=J_n+\sum_{j=n}^{m-1}\Delta_j=J_n+c(\rho^m-\rho^n)/(\rho-1)+\sum_{j=n}^{m-1}(\Delta_j-c\rho^j)\), the quantity \(J_m-c\rho^m/(\rho-1)\) converges once \(E_n^\pm<\infty\), to \(C'=J_n-c\rho^n/(\rho-1)+\sum_{j\ge n}(\Delta_j-c\rho^j)\), and the tail sum lies in \([-E_n^+,E_n^-]\). \(E_n^\pm<\infty\): on \((0,0.01]\), \(r_\infty\le Bt^{-0.35}\) and \(W_\pm\le Ct^{0.65}\) (AT's small-\(y\) lemma, \(\log A_\infty\le t/20\)), so the integrand is \(O(t^{0.3-\alpha})\), integrable for \(\alpha<1.3\); on \([1,\infty)\) AT's crude tail applies. Monotonicity in \(n\) is \(r_m\uparrow\). The rate \(\varepsilon_m=O(\rho^m5^{-0.65m})\) is AT (4a): \(|\Delta_j-c\rho^j|=O(\rho^j5^{-0.65j})\), summed geometrically since \(\rho5^{-0.65}=5^{\alpha-0.65}/3<1\).

*(2)* AT (1) is a change of variables \(s=q^k(t)\) and holds verbatim for complex exponents: \(s^{z-1}\,ds=t^{z-1}A_k^{1-z}5^{-k(z-1)}\,5^{-k}B_k\,dt\) and \(\delta_{m+1}(s)=(5/3)^kC_k\delta_{n+1}(t)\), so \(\rho_z^{-m}\int s^{z-1}\delta_{m+1}=\rho_z^{-n}\int t^{z-1}\delta_{n+1}A_k^{1-z}B_kC_k\). With \(|A_k^{1-z}|\le A_\infty^{1-\operatorname{Re}z}\), \(C_k\uparrow C^{(n)}_\infty\le\bar C\), \(\delta_{n+1}\) bounded at \(0\) and doubly-exponentially small at \(\infty\), dominated convergence gives \(c(z)=\rho_z^{-n}\int t^{z-1}\delta_{n+1}A_\infty^{1-z}B_\infty C^{(n)}_\infty\), holomorphic on \(0<\operatorname{Re}z<1\) (locally uniform domination), and \(e_j(z)=\int t^{z-1}\delta_{j+1}(1-A_\infty^{1-z}B_\infty C^{(j)}_\infty)\). On \((0,1]\), \(|1-A_\infty^{1-z}B_\infty C^{(j)}_\infty|\le K_z\,t^{0.65}\) with \(K_z\) locally bounded (the same small-\(y\) estimates), and on \([1,\infty)\) it is \(\le1+3\cdot10^{10}t^5\). So \(\sum_{j\ge n}|e_j(z)|\le K_z\int t^{\operatorname{Re}z-1+0.65}(r_\infty-r_n)+(\text{tail})<\infty\) for \(\operatorname{Re}z>-0.3\), locally uniformly. For \(0.35<\operatorname{Re}z<1\), \(\int s^{\operatorname{Re}z-1}(r_\infty-r_n)<\infty\) (power ceiling) and \(|\rho_z|<1\), so Fubini gives \(M(z)-I_n(z)=\sum_{j\ge n}\int s^{z-1}\delta_{j+1}=c(z)\rho_z^n/(1-\rho_z)+\sum_{j\ge n}e_j(z)\). \(I_n\) is holomorphic on \(0<\operatorname{Re}z<1\) (\(r_n\) bounded at \(0\), \(\le1/s\) at \(\infty\)), so the right side is meromorphic there and is the continuation of \(M\) (identity theorem); its poles are where \(\rho_z=1\), with \(1-\rho_z=\ln5\,(z-z_k)+O((z-z_k)^2)\) and \(\rho_{z_k}^n=1\). At real \(z=1-\alpha\), \(\kappa I_n=R_\alpha(n)\), \(\kappa c=c_R\), \(\kappa e_j=\kappa(\Delta_j-c\rho^j)\), so \(\kappa M(1-\alpha)=R_\alpha(n)-c_R\rho^n/(\rho-1)+\kappa\sum_{j\ge n}(\Delta_j-c\rho^j)=C_\alpha\) by (1). Near \(\alpha=\gamma\), \(M(1-\alpha)=\langle\Theta\rangle/(\gamma-\alpha)+O(1)\) (AU (2)), which gives the pole statement. For continuity, \(G(\alpha):=R_\alpha(n)-c_R(\alpha)(\rho^n-1)/(\rho-1)-\varepsilon_n(\alpha)\) equals \(C_\alpha+c_R(\alpha)/(\rho-1)\) for \(\alpha\ne\gamma\) and \(C_\gamma\) at \(\gamma\) (AY (1), where \((\rho^n-1)/(\rho-1)=n\)); \(R_\alpha(n)\) and \(c_R(\alpha)\) are continuous in \(\alpha\), \((\rho^n-1)/(\rho-1)\) is a polynomial in \(\rho\), and \(\sup_{|\alpha-\gamma|\le\eta}|\varepsilon_n(\alpha)|\to0\) as \(n\to\infty\) by the locally uniform bound above, so \(G\) is a locally uniform limit of continuous functions.

*(3)* As in AY (2), with \(t^{-\gamma}\) replaced by \(t^{-\alpha}\) and \(A_\infty^{\gamma}\) by \(A_\infty^{\alpha}\) (AY's routines are reused unchanged with the exponent passed in). New here is the enclosure of \(c\rho^n\) from AT (2) at level \(n\): cellwise \(\int_{\rm cell}t^{-\alpha}\delta_{n+1}\in[\mathrm{floor}(r_{n+1})-\mathrm{ceil}(r_n),\ \mathrm{ceil}(r_{n+1})-\mathrm{floor}(r_n)]\) (clipped at \(0\)), times \(1-\sup W_-\le B_\infty\) and \(1+\sup W_+\ge A_\infty^{\alpha}\bar C\) on blocks of \(16\) cells; head \((0,x_0]\) by \(r_{n+1}(x_0)-r_n(0)\le\delta_{n+1}\le r_{n+1}(0)-r_n(x_0)\); upper tail on \([1,\infty)\) by AT's crude tail (\(<10^{-4\cdot10^{17}}\) at \(n=60\)). Then \(C'\in[J_n-c\rho^n/(\rho-1)-E_n^+,\ J_n-c\rho^n/(\rho-1)+E_n^-]\), at \(n=60\) and \(n=40\) independently, intersected; \(c\) likewise. \(\sigma_\infty(\alpha)=\frac{\sin\pi\alpha}{\pi}\int s^{-\alpha}(1-sr_\infty)\) (AS (1), \(\beta-1=-\alpha\)) with AY's head and the tail \([2/s-6/s^2,2/s]\).

*(4)* \(\|F(n,\alpha)\|=\frac{\sqrt{21}}{3}\sigma_n/R_\alpha(n)\) (AQ (3)). By (1), \(R_\alpha(n)=A_\alpha(\rho^n+\nu_\alpha)+O(\rho^n5^{-0.65n})\); by AR (2)/AT (4a), \(\sigma_n=\sigma_\infty+O(\mu^n)\), and \(\mu=5^{\alpha-1}/3\le1/3<5^{-0.65}\). Hence \((\rho^n+\nu_\alpha)\|F\|=\frac{\sqrt{21}\sigma_n}{3A_\alpha}(1+O(5^{-0.65n}))=K_\alpha+O(5^{-0.65n})\). \(\square\)

**Witness, computed by `scripts/gasket_fractional_supercritical_shift.py`** (about 6 min).
- Re-certified: AQ's power ceiling (\(B\approx0.892577\)) and the AT small-\(y\) lemma (ratio \(1.7896\le1.82\)). Node enclosures: \(r_\infty\) to relative width \(<10^{-60}\), and \(r_{40}\le r_{41}\le r_{60}\le r_{61}\le r_\infty\) at every node.
- \(J_{60}(0.8)\in[536714.60891607675812942124,\,\ldots]\) and \(c\rho^{60}\in[111619.41616612443669242428,\,\ldots]\), each of width \(<10^{-24}\). The levels-\(60\) and \(40\) enclosures of \(c\) and of \(C'\) overlap. \(E_{40},E_{60}<\rho^n5^{-0.6n}\). Theorem AT's independent \(c_R(0.75)\) contains this script's level-\(60\) value. The enclosures above.
- Consistency (40-digit quadrature in \(\log s\), not part of the proof), \(\alpha=0.8\): \(R(m)-A\rho^m-C=2.96\cdot10^{-3},\,3.11\cdot10^{-5},\,3.48\cdot10^{-9},\,3.90\cdot10^{-13}\) and \((\rho^m+\nu)\|F(m,0.8)\|-K=-2.75\cdot10^{-4},\,-1.36\cdot10^{-6},\,-2.31\cdot10^{-11},\,-3.88\cdot10^{-16}\) at \(m=5,10,20,30\), inside \(5\rho^m5^{-0.65m}\) and \(5\cdot5^{-0.65m}\). **OBSERVATION, not a theorem:** the ratio of \(R(m)-A\rho^m-C\) from \(m=20\) to \(30\) is \(1.12029\cdot10^{-4}\), against \((\rho/3)^{10}=1.12030\cdot10^{-4}\): the same \(3^{-m}\) correction AT and AY observed, carried by \(\rho^m\).
- Pole and continuity (quadrature at levels \(29,30\), not part of the proof): \(C_{\gamma+0.01}=-31.113\), \((\alpha-\gamma)C_\alpha=-0.3111\) against \(-c_R(\gamma)/\ln5=-0.3213\); \(R_\infty(\gamma-0.01)=33.160\); \(C_\alpha+c_R(\alpha)/(\rho-1)=0.69733\) and \(0.70909\) at \(\gamma\pm0.01\), symmetric mean \(0.70321\) against AY's \(C_\gamma=0.70307\).

The script exits nonzero if any check fails.

**KEPT FAILURE AZ.1 (the endpoint shift does not carry over with its sign; the unshifted renormalization is off by a full power).** The claim "above \(\gamma\), as at \(\gamma\), the resistance overshoots its leading term (\(C_\alpha>0\)), so \(\rho^n\|F(n,\alpha)\|\) increases to \(K_\alpha\) from below, with an error that is \(O(5^{-0.65n})\) or at least smaller than \((3\cdot5^{-\alpha})^n\)" is false. At \(\gamma\), AY has \(\nu_\gamma>0\) and \(n\|F(n,\gamma)\|<K_\gamma\). At \(\alpha=0.8\), \(C_{0.8}\le-1.74325717839\) and \(\nu_{0.8}\le-0.72744538444\) are certified, so by (4) \(\rho^n(\rho^n\|F(n,0.8)\|-K_{0.8})\to-\nu_{0.8}K_{0.8}=0.5003854985\ldots>0\): the renormalized norm decreases to \(K_{0.8}\) from above, and the error of the unshifted renormalization is \(\Theta((3\cdot5^{-0.8})^n)\), the same order as the shift itself. The sign flip is not an accident of \(0.8\): by (2), \(C_\alpha\) is the continuation of \(R_\infty(\alpha)\) through a simple pole with positive residue, so \(C_\alpha\to-\infty\) as \(\alpha\downarrow\gamma\) (\(C_{0.75}<0\) is certified too). Witness (same script, 40-digit quadrature): \(\rho^m\|F(m,0.8)\|-K_{0.8}=0.27089,\,0.084985,\,0.011628,\,0.0017328\) at \(m=5,10,20,30\) (positive, decreasing); \(\rho^{30}(\rho^{30}\|F\|-K_{0.8})=0.5016460\), matching \(-\nu K\rho^{30}/(\rho^{30}+\nu)\) to \(10^{-8}\).

**OPEN.** Closed forms for \(c_R(\alpha)\), \(C_\alpha\), \(\sigma_\infty(\alpha)\), \(K_\alpha\); whether \(C_\alpha<0\) on all of \((\gamma,1)\) (proved only near \(\gamma\) and certified at \(0.75\), \(0.8\)); a proof of the observed \((\rho/3)^m\) correction; whether \(\hat\theta_k\ne0\) for every \(k\) (the complex poles of (2)); and a flux functional on the infinite gasket. Theorem AZ renormalizes one scalar audit norm.

**NOT THIS.** Theorem AZ gives the exact decay law of one vanishing graph audit norm on the combinatorial gasket above \(\log3/\log5\), identifies its intercept as a continued Mellin transform, and certifies the constants at \(\alpha=0.8\) (and \(0.75\)) to about 22 digits. It gives no closed form. It does not construct a flux on the infinite gasket. It says nothing about the multi-scale or geometric-weight builds, Stage 1, or Stage 2; the Stage 2 divergence question stays OPEN. A decay constant of a vanishing graph audit norm is not a force: it is not thrust and not a pass of any experimental criterion. Not thrust.

## Theorem BG. The fractional Dirichlet extension itself has a refinement limit: flat for \(\alpha\le\log3/\log5\), not flat above, \(U(0.45)\in[0.143810693,0.143889868]\)

Theorems AQ–AZ give limits of scalar audit norms. This theorem asks for a limit of the object behind them, the \(L_n^{\alpha}\)-Dirichlet extension \(u_n\) itself, at a fixed vertex of the gasket as \(n\to\infty\).

**Setting.** Assumption 4, level \(n\), \(\alpha\in(0,1)\), corner data \(d\) on \(C=\{a,b,c\}\), \(\bar d\) their mean, \(x=d-\bar d\), \(\|x\|_\infty=\max_c|x_c|\). \(u_n\) solves \((L_n^{\alpha}u)_i=0\) off \(C\), \(u_n=d\) on \(C\). \(V_m\subset V_n\) for \(m\le n\) (same positions), \(V_*=\bigcup_mV_m\). Notation of Theorems AQ–AZ and README BB–BD: \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\), \(p(s)=s(s+5)\), \(\rho=5^{\alpha}/3\), \(\gamma=\log3/\log5\), \(D_n=\int_0^\infty s^{-\alpha}r_n\,ds\), so \(R_\alpha(n)=\frac{2\sin\pi\alpha}{\pi}D_n\) (AQ (2)), and \(D_\infty=\lim D_n\) for \(\alpha<\gamma\). \(\psi_n(\tau)=p^n(5^{-n}\tau)\uparrow\psi(\tau)\) and \(\mathcal Y_n(\tau)=3^ny_n(5^{-n}\tau)\to\mathcal Y(\tau)\) are BD.4's scaling limits. For \(x\in V_m\) and \(z\ge0\) put
\[
h_m(z;x)=\frac{[(L_m+z)^{-1}x_C]_x}{r_m(z)},
\]
with \(x_C\) the corner vector \(x\) extended by \(0\). **Assumption BG (natural renormalisation, used only in BG.4).** The renormalised fractional energy is \(\mathcal E_n=3^{-n}5^{n\alpha}Q_\alpha(u_n)=\rho^nQ_\alpha(u_n)\), with \(Q_\alpha(u)=u^{T}L_n^{\alpha}u\): weight \(3^{-n}\) per vertex and \(L_n\) rescaled by \(5^n\), constants dropped. \(u_n\) does not depend on any scalar renormalisation of \(L_n\), so BG.1–BG.3 need no such choice. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters. **Assumptions used as stated:** AQ (1)–(3), AR (3)–(4), AS (1) (re-run live at \(0.45\)), AY (\(R_\gamma(n)=c_R(\gamma)n+C_\gamma+o(1)\)), AZ (1), (3) (\(R_\alpha(n)=A_\alpha\rho^n+C_\alpha+o(1)\) for \(\alpha>\gamma\); certified \(A_{0.8}\)), and BD.4 (convergence of \(\psi_n,\mathcal Y_n\) with explicit tail factors).

**THEOREM BG.**

1. **BG.1 (exact profile representation; every \(n\ge m\), every \(\alpha\), every \(x\in V_m\)).**
   \[
   u_n(x)=\bar d+\frac{\int_0^\infty s^{-\alpha}r_n(s)\,h_m\bigl(p^{n-m}(s);x\bigr)\,ds}{\int_0^\infty s^{-\alpha}r_n(s)\,ds}.
   \]
   So \(u_n(x)-\bar d\) is the mean of the level-\(m\) resolvent profile \(h_m(z;x)\) under the law of \(z=p^{n-m}(s)\), \(s\propto s^{-\alpha}r_n\,ds\). The profile satisfies \(h_m(z;c)=x_c\) on the corners, \(h_{m+1}(z;\cdot)=h_m(p(z);\cdot)\) on \(V_m\), \(h_m(0;\cdot)+\bar d\) is the harmonic (\(\alpha=1\)) extension, and \(|h_m(z;x)|\le4\|x\|_\infty/(4+z)\) off the corners. At level \(1\), \(h_1(z;m_{ab})=-x_c/(5+z)\), and cyclically for \(m_{bc},m_{ca}\).
2. **BG.2 (flat pointwise limit for \(\alpha\le\gamma\)).** For every \(\alpha\in(0,1)\), every \(m\), every \(n\ge m\) and every \(x\in V_m\setminus V_0\),
   \[
   |u_n(x)-\bar d|\le\frac83\,4^{-\alpha}\,5^{m(1-\alpha)}\,\|x\|_\infty\,\frac{\rho^n}{R_\alpha(n)} .
   \]
   For \(\alpha<\gamma\) the right side is \(O(\rho^n)\) (\(R_\alpha(n)\uparrow R_\infty<\infty\)), and at \(\alpha=\gamma\) it is \(O(1/n)\) (\(\rho=1\), \(R_\gamma(n)\sim c_R(\gamma)n\)). So **\(u_n\to\bar d\) at every vertex of \(V_*\setminus V_0\)**: the limit object is the constant \(\bar d\), with the data surviving only on the three corners.
3. **BG.3 (the renormalised profile; every \(\alpha\in(0,1)\)).** With \(\mathcal R(\tau)=\lim(3/5)^nr_n(5^{-n}\tau)=\mathcal Y(\tau)\psi(\tau)/\tau\),
   \[
   \rho^{-n}\int_0^\infty s^{-\alpha}r_n\,h_m(p^{n-m}s;x)\,ds\ \to\ P_m(x)=\int_0^\infty\tau^{-\alpha}\mathcal R(\tau)\,h_m\bigl(\psi(5^{-m}\tau);x\bigr)\,d\tau .
   \]
   Hence:
   - for \(\alpha<\gamma\), \(\rho^{-n}(u_n(x)-\bar d)\to P_m(x)/D_\infty\);
   - at \(\alpha=\gamma\), \(n(u_n(x)-\bar d)\to\frac{2\sin\pi\gamma}{\pi}P_m(x)/c_R(\gamma)\);
   - for \(\alpha>\gamma\), \(u_n(x)\to u_\infty(x)=\bar d+\frac{2\sin\pi\alpha}{\pi}P_m(x)/A_\alpha\). This limit is **not flat**.

   At the midpoint \(m_{ab}\), \(P_1=-x_cW(\alpha)\) with \(W(\alpha)=\int_0^\infty\tau^{-\alpha}\mathcal R(\tau)/(5+\psi(\tau/5))\,d\tau=\int_0^\infty\tau^{-1-\alpha}\mathcal Y(\tau)\psi(\tau/5)\,d\tau\). Certified (mpmath.iv):
   - \(W(0.45)\in[0.5754299822,0.5757467849]\), and with AS's \(D_\infty(0.45)\in[4.0013017734819,4.001301773482]\), \(U(0.45)=W/D_\infty\in[0.143810693,0.143889868]\). So \(u_n(m_{ab})-\bar d\sim-x_c\,U(0.45)\,\rho^n\) with \(\rho=0.68772568\ldots\). For the data \((1,-\tfrac12,0)\), \(x_c=-\tfrac16\), the renormalised profile on \(V_1\) is \(\rho^{-n}(u_n-\tfrac16)\to(U/6,\,-5U/6,\,2U/3)\) at \((m_{ab},m_{bc},m_{ca})\).
   - \(W(0.8)\in[0.5512217593,0.5513327730]\), and with AZ's \(A_{0.8}\), \(\kappa_{0.8}=\frac{2\sin0.8\pi}{\pi}W/A_{0.8}\in[0.086072435,0.086089770]\). So \(\lim u_n(m_{ab})=\bar d-\kappa_{0.8}x_c\), which is \([0.181012073,0.181014962]\) for the data \((1,-\tfrac12,0)\). The harmonic value is \(0.2\), i.e. \(\kappa=1/5\).
4. **BG.4 (renormalised energy; DERIVED under Assumption BG).** \(Q_\alpha(u_n)=q^{T}Gq=2\|x\|_2^2/R_\alpha(n)\) (\(=7/(3R_\alpha(n))\) for \((1,-\tfrac12,0)\), Theorem AD's pairing). So \(\mathcal E_n=2\|x\|_2^2\rho^n/R_\alpha(n)\) tends to \(0\) for \(\alpha<\gamma\) (like \(\rho^n\)) and at \(\gamma\) (like \(1/n\)), and to \(2\|x\|_2^2/A_\alpha>0\) for \(\alpha>\gamma\) (\(7/(3A_{0.8})=0.97366\ldots\) at \(0.8\)). The flat and non-flat regimes of BG.2–BG.3 are exactly the zero and nonzero renormalised energy regimes.
5. **BG.5 (corner neighbours; DERIVED).** For the data \((1,-\tfrac12,0)\), the mean of the two neighbours of corner \(a\) is \(1-I_a/2=1-5\|F(n,\alpha)\|/(2\sqrt{21})\) (Theorem J's Schur form). For \(\alpha>\gamma\), \(\|F\|\to0\) (Theorem K), so it tends to \(u_a=1\): no jump at the corner. For \(\alpha<\gamma\), see Kept Failure BG.1.

**Proof of BG.1.**
- *Charge form.* README Theorem BB.1 gives \(u_n=Gq+\bar d\mathbf 1\), \(q=(2/R_\alpha(n))x_C\), \(G=(L_n^{\alpha})^{+}\). Since \(q\perp\mathbf 1\), \(Gq=\frac{\sin\pi\alpha}{\pi}\int s^{-\alpha}(L_n+s)^{-1}q\,ds\).
- *One decimation step.* By AQ (1), the Schur complement of \(L_{k+1}+s\) onto \(V_k\) is \(c(s)(L_k+p(s))-e(s)P_C\), with \(c=(6+s)/((5+s)(2+s))\) and \(e=2s/(2+s)\), because \(c\,p=2e+s\) (sympy). Restricting the inverse to \(V_k\) inverts the Schur complement.
- *Zero-sum corner mode.* By \(S_3\) symmetry, the corner block of \(Z=(L_k+p(s))^{-1}\) acts on the zero-sum plane as \(r_k(p(s))\). So \((cZ^{-1}-eP_C)\,Zx_C/(c-er_k(ps))=x_C\) (sympy), and \([(L_{k+1}+s)^{-1}x_C]|_{V_k}=Zx_C/(c-er_k(ps))\).
- *Telescoping.* AQ's \(r_{k+1}(s)=r_k(ps)/(c-er_k(ps))\) makes the factor \(r_{k+1}(s)/r_k(ps)\). Over \(k=m,\ldots,n-1\) this gives \([(L_n+s)^{-1}x_C]|_{V_m}=r_n(s)\,h_m(p^{n-m}s;\cdot)\). Integrate and divide by \(R_\alpha(n)/2\).
- *Properties of \(h\).* The corner block gives \(h_m(z;c)=x_c\). One step of the telescoping gives \(h_{m+1}(z)=h_m(p(z))\). At \(z=0\), \(L_m^{+}x_C\) is harmonic off \(C\) with corner values \(r_m(0)x_c\).
- *Maximum principle.* \(w=(L_m+z)^{-1}x_C\) satisfies \((\deg y+z)w_y=\sum_{y'\sim y}w_{y'}\) off \(C\). So \(\max|w|\) is attained on \(C\), where it equals \(r_m(z)\|x\|_\infty\), and \(|w_y|\le\frac{\deg y}{\deg y+z}\max|w|\le\frac4{4+z}r_m(z)\|x\|_\infty\).
- *Level 1.* The \(6\times6\) solve is done in sympy. \(\square\)

**Proof of BG.2.** Use \(p^{n-m}(s)\ge5^{n-m}s\), \(r_n(s)\le r_n(0)=\tfrac13(5/3)^n\) (\(r_n\) is Stieltjes, and the corner resistance is \(\tfrac23(5/3)^n\)), and \(\int_0^\infty s^{-\alpha}\frac{4}{4+s/K}ds=4K(4K)^{-\alpha}\pi/\sin\pi\alpha\) with \(K=5^{m-n}\) (sympy check). The numerator of BG.1 is at most \(\|x\|_\infty\frac43\,4^{-\alpha}5^{m(1-\alpha)}\rho^n\pi/\sin\pi\alpha\). Divide by \(D_n=\pi R_\alpha(n)/(2\sin\pi\alpha)\). The rates come from AQ (3) (\(R_\alpha(n)\uparrow R_\infty<\infty\) for \(\alpha<\gamma\)) and AY. \(\square\)

**Proof of BG.3.**
- *Scaling.* Substitute \(s=5^{-n}\tau\). Then \(p^{n-m}(5^{-n}\tau)=\psi_{n-m}(5^{-m}\tau)\) and \((3/5)^nr_n(5^{-n}\tau)=\mathcal Y_n\psi_n/\tau=:\mathcal R_n\), from \(y_n=r_n\,s/p^n(s)\) (BB.2, BD.3). This gives \(\rho^{-n}\times\)numerator \(=\int\tau^{-\alpha}\mathcal R_n(\tau)\,h_m(\psi_{n-m}(5^{-m}\tau);x)\,d\tau\).
- *Pointwise limit.* \(\mathcal R_n\to\mathcal R\) and \(\psi_{n-m}\to\psi\) (BD.4), and \(h_m\) is continuous in \(z\).
- *Domination.* \(\mathcal R_n\le\tfrac13\), and \(|h_m|\le4\|x\|_\infty/(4+5^{-m}\tau)\) because \(\psi_k(t)\ge t\). So \(\tau^{-\alpha}\tfrac43\|x\|_\infty/(4+5^{-m}\tau)\) is an integrable dominator for every \(\alpha\in(0,1)\).
- *The three regimes.* The denominators are \(D_n\to D_\infty\) (\(\alpha<\gamma\)), \(D_n/n\to\pi c_R(\gamma)/(2\sin\pi\gamma)\) (AY), and \(\rho^{-n}D_n\to\pi A_\alpha/(2\sin\pi\alpha)\) (AZ (1)). The two integral forms of \(W\) agree because \(1/(5+\psi(t))=\psi(t)/\psi(5t)\) (\(\psi(5t)=p(\psi(t))\)).
- *Certificate.*
  - Cells: \(\mathcal R\) is nonincreasing (\(r_n\) is) and \(\psi\) is increasing. So on each of \(13485\) cells of the BD grid (\([10^{-20},3125]\), ratios \(1.05\) and \(1.001\)), the integrand lies between its outward bounds at the two endpoints.
  - Enclosures: \(\mathcal R\) is enclosed by \(\mathcal R_{60}\) with BD.4's \(e^{\pm\varepsilon^{\mathcal Y}}\) and the new factor \(\psi_N\le\psi\le\psi_Ne^{2^N5^{-N}\tau/10}\). That factor holds because \(5^Np^k(5^{-N-k}\tau)\le\tau e^{5^{-N}\tau/10}\) (bootstrap \(p^j(y)\le2\cdot5^jy\)), and \(\psi_N\) has nonnegative coefficients and degree \(2^N\).
  - Head \([0,10^{-20}]\): \(\mathcal R\le\tfrac13\), \(1/(5+\psi)\le\tfrac15\), and the endpoint floor below.
  - Tail \([3125,\infty)\): \(\mathcal R\le\tfrac13\) and \(\psi(\tau/5)\ge\psi_6(\tau/5)\ge(\tau/3125)^8P_8(625)\), with BD's exact coefficient sum.
  - \(D_\infty(0.45)=J_r\) is Theorem AS's certified integral, re-run live. \(A_{0.8}\in[2.39640970397922052703,2.39640970397922052705]\) is Theorem AZ (3) as stated. \(\square\)

**Proof of BG.4 and BG.5.** \(u^{T}L^{\alpha}u=(Gq)^{T}L^{\alpha}(Gq)=q^{T}Gq=(2/R)^2(R/2)\|x\|_2^2\) (sympy). Then apply AQ (3), AY and AZ (1). BG.5 is Theorem J's identity \(\|F\|=I_a\sqrt{21}/5\), with \(I_a=2u_a-u_{p}-u_{p'}\), together with Theorem K. \(\square\)

**KEPT FAILURE BG.1 (the flat limit is not uniform: a corner boundary layer).** Claim: for \(\alpha<\gamma\), \(u_n\to\bar d\) uniformly on \(V_n\setminus V_0\), i.e. \(u_n\) converges to a continuous function on the gasket, as at \(\alpha=1\). This is false. **Witness:** for the data \((1,-\tfrac12,0)\) at \(\alpha=0.45\), the mean of the two neighbours of corner \(a\) is \(1-5\|F(n)\|/(2\sqrt{21})\). Since \(\|F(n)\|\downarrow\mathcal L(0.45)\) (AQ, AS), it increases to \(1-5\mathcal L/(2\sqrt{21})\), so
\[
\lim_n\Bigl(\text{neighbour mean}-\bar d\Bigr)=\frac56-\frac{5\mathcal L(0.45)}{2\sqrt{21}}\in[0.3853654863,\,0.3853654863]\quad\text{(certified, AS re-run)} .
\]
So \(\sup_{V_n\setminus V_0}|u_n-\bar d|\) stays above \(0.3853\) in the limit, while \(u_n\to\bar d\) at every fixed vertex (BG.2). The dense values at levels \(1\)–\(6\) are \(0.214634,0.363945,0.436850,0.478347,0.503722,0.519861\) (increasing; \(\bar d=0.1\overline{6}\)). The data are carried by a layer at the corners whose width shrinks with the mesh. (Interpretation, not proved here: this matches the continuum picture in which points are polar for the fractional operator when \(\alpha<d_s/2=\gamma\), \(d_s=2\log3/\log5\).) The bound of BG.2 is uniform on each \(V_m\) but grows like \(5^{m(1-\alpha)}\), so it cannot be made uniform.

**Evidence (floats, recursion quadrature; consistency only).** \(\rho^{-n}W_n\), \(W_n=\int s^{-\alpha}r_n/(5+p^{n-1}(s))\):
- \(0.45\): \(0.57557217\), \(0.57558823942\), \(0.57558823969\) at \(n=10,20,30\); \(\rho^{-30}W_{30}/D_{30}=0.143851453\).
- \(0.8\): \(0.55127725784\) at \(n=30\); \(W_{30}/D_{30}=0.0862980\), above the limit \(\kappa_{0.8}\), consistent with \(C_{0.8}<0\).
- \(\gamma\): \(nW_n/D_n=0.4402,0.4682,0.4783\) at \(n=10,20,30\), consistent with \((n+\nu_\gamma)W_n/D_n\to\frac{2\sin\pi\gamma}{\pi}W(\gamma)/c_R(\gamma)\approx0.501\).
- The dense \(\rho^{-6}(u_6(m_{ab})-\tfrac16)=0.025637\), against the limit \(U/6\in[0.023968,0.023982]\).

**CONJECTURE BG.2 (PROVED by Theorem BI below; the text records it as first stated).** For \(\alpha>\gamma\), \(u_n\to u_\infty\) uniformly on \(V_*\), and \(u_\infty\) extends continuously to the gasket with \(u_\infty=d\) on \(C\). Evidence: BG.5, the non-flat \(V_1\) limit at \(0.8\), and nothing else. No modulus of continuity is proved. Closed forms for \(W(\alpha)\), \(U(\alpha)\), \(\kappa_\alpha\) and \(P_m\) are not known, and neither is a characterisation of \(u_\infty\) or of the renormalised profile \(P_m/D_\infty\) as solutions of a continuum equation on the gasket.

**NOT THIS.** Theorem BG is the pointwise refinement limit of one graph Dirichlet extension on the combinatorial gasket. Below \(\gamma\) the limit is a constant with the data left on three points. That is a statement about polar corners of a discrete fractional operator, not a field, a flux or a force. No continuum flux functional on the infinite gasket is constructed. The renormalised energy and profile are refinement rates. Not thrust.

**Witness, computed by `scripts/gasket_fractional_extension_limit.py`** (about 3 min).
- *Sympy:* \(r_1\), the level-1 profile \(h_1=-x_c/(5+z)\) and its harmonic value; \(c\,p=2e+s\); the Woodbury step; \(1/(5+\psi)=\psi/\psi(5\cdot)\); the dominator integral; \(Q=2\|x\|^2/R\).
- *Exact rationals (levels \(1\)–\(3\), \(s\in\{1/2,7/5\}\)):* the corner mode, the maximum-principle bound, and BG.1's restriction identity for every \(m<n\), against exact resolvent solves.
- *Dense meshes (\(\alpha=0.45,0.8\), levels \(1\)–\(6\)):* \(u_n(m_{ab})\) equals BG.1 at levels \(1,2,4\); BG.2's bound holds at every vertex of \(V_1,V_2,V_3\); the corner-neighbour means.
- *Recursion quadrature:* the evidence table above.
- *Intervals:* \(W(0.45)\), \(W(0.8)\), \(D_\infty(0.45)\) and \(\mathcal L(0.45)\) (AS re-run), \(U(0.45)\), \(\kappa_{0.8}\), and the BG.1 gap.

The script exits nonzero if any check fails.

## Theorem BH. The corner boundary layer of Kept Failure BG.1 has a universal limit profile: \(\theta_1(0.45)\in[0.340381375,0.340407195]\), none above \(\log3/\log5\)

Kept Failure BG.1 showed that the flat limit of Theorem BG is not uniform: near each corner a layer keeps a gap. This theorem resolves that layer. It proves the layer converges at every vertex, that its profile is universal, and that it is absent for \(\alpha\ge\gamma\). It also gives an explicit series for the neighbour profile and certified values at four layer vertices at \(0.45\).

**Setting.** As in Theorem BG: Assumption 4, level \(n\), \(\alpha\in(0,1)\), corner data \(d\), \(x=d-\bar d\), \(\|x\|_\infty=\max|x_c|\), \(u_n\) the \(L_n^{\alpha}\)-Dirichlet extension, \(r_n\), \(p(s)=s(s+5)\), \(D_n=\int_0^\infty s^{-\alpha}r_n\), \(h_m(z;\cdot)\), \(\rho=5^{\alpha}/3\), \(\gamma=\log3/\log5\). Positions are skew lattice coordinates with \(a=(0,0)\), \(b=(1,0)\), \(c=(0,1)\). The corner cell of size \(2^{-(n-k)}\) at \(a\) is \(2^{-(n-k)}V_k\subset V_n\). **Layer coordinates:** \(\xi\in\mathbb Z_{\ge0}^2\) labels the vertex \(2^{-n}\xi\) of \(V_n\) (once \(n\) is large enough). The neighbours of \(a\) in \(V_N\) are \(\nu_b=(2^{-N},0)\) and \(\nu_c=(0,2^{-N})\). The level-1 \(s\)-harmonic weights (sympy) are \(\alpha_1(s)=(s+4)/((s+2)(s+5))\) on \(a\) and \(b\) and \(2/((s+2)(s+5))\) on \(c\), seen from \(m_{ab}\). Put \(\beta(s)=(s+6)/((s+2)(s+5))\). Let \(H_k(s;v,y)\), for \(v\in V_0\) and \(y\in V_k\), be the \(s\)-harmonic measure of \(V_k\): solve \((\deg y+s)f_y=\sum_{y'\sim y}f_{y'}\) off \(V_0\) with \(f=\mathbf 1_v\) on \(V_0\). \(r_\infty=\lim r_n\) (AQ (3)), \(J_r=\int s^{-\alpha}r_\infty\) and \(J_w=\int s^{-\alpha}(1-sr_\infty)\) (AS (1)), and \(y_n=r_n\prod_{j<n}(5+p^j)^{-1}\) (BB). **Assumptions used as stated:** AQ (1)–(3) (decimation, Stieltjes formulas, \(r_n\uparrow r_\infty\), \(R_\alpha(n)\uparrow R_\infty<\infty\) for \(\alpha<\gamma\)), AS (1) (\(J_r,J_w\), re-run live at \(0.45\)), AY and AZ (1) (\(D_n\to\infty\) for \(\alpha\ge\gamma\)), BG.1, and, in BH.6 only, AT's small-\(y\) lemma and \(\Theta\) bracket and BD.4's \(\psi_k\uparrow\psi<\infty\). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BH.**

1. **BH.1 (corner-cell decomposition; exact, every \(n\ge k\), \(s>0\), \(y\in V_k\)).** \([(L_n+s)^{-1}x_C](2^{-(n-k)}y)\) is the \(s\)-harmonic extension on \(V_k\) of its three cell-corner values \(r_n(s)x_a\), \(r_n(s)h_{n-k}(p^ks;\nu_b)\), \(r_n(s)h_{n-k}(p^ks;\nu_c)\), with \(\nu_b,\nu_c\) taken in \(V_{n-k}\). So
   \[
   u_n(2^{-(n-k)}y)=\bar d+\frac1{D_n}\int_0^\infty s^{-\alpha}r_n(s)\Bigl[x_aH_k(s;a,y)+\sum_{v=b,c}h_{n-k}(p^ks;\nu_v)H_k(s;v,y)\Bigr]ds .
   \]
2. **BH.2 (the corner-neighbour profile).** For every \(N\ge0\) and \(z\ge0\):
   - *asymmetry, exact:* \(h_N(z;\nu_b)-h_N(z;\nu_c)=(x_b-x_c)\prod_{i<N}(5+p^iz)^{-1}\);
   - *symmetric part, uniform contraction:* \(|h_N(z;\nu_v)-x_a\eta(z)|\le2\|x\|_\infty\prod_{i<N}\beta(p^iz)\le2\|x\|_\infty(3/5)^N\), where
   \[
   \eta(z)=\sum_{j\ge0}\alpha_1(p^jz)\prod_{i<j}\beta(p^iz),\qquad \eta(z)=\alpha_1(z)+\beta(z)\,\eta(p(z)) .
   \]
   \(\eta\) is continuous and strictly decreasing on \([0,\infty)\), with \(\eta(0)=1\) and \(0<\eta(z)\le2/(2+z)\). It is \(\alpha\)-free.
   - *corner resolvent:* \(r_\infty(s)=1/(2+s-2\eta(s))\). This series form is new. It reproduces AQ's fixed point \(r_\infty(s)=r_\infty(ps)/(c-er_\infty(ps))\) with \(c=\beta\) (sympy).
   - *Corollary (exact):* \(u_n(\nu_b)-u_n(\nu_c)=(x_b-x_c)\,I_n/D_n\) with \(I_n=\int s^{-\alpha}y_n\), and \(0<I_n/D_n\le5^{-n}\). So both neighbours of \(a\) carry the layer separately, not only their mean (BG.1 used only the mean).
3. **BH.3 (layer limit for \(\alpha<\gamma\)).** For every fixed \(k\) and \(y\in V_k\),
   \[
   \Bigl|u_n(2^{-(n-k)}y)-\bar d-x_a\Theta_k(y)\Bigr|\le2\|x\|_\infty(3/5)^{n-k}+2|x_a|\frac{J_r-D_n}{D_n}\ \to0,
   \]
   \[
   \Theta_k(y)=\frac1{J_r}\int_0^\infty s^{-\alpha}r_\infty(s)\,\phi_k(s;y)\,ds,\qquad \phi_k=H_k(s;a,y)+\eta(p^ks)\bigl(H_k(s;b,y)+H_k(s;c,y)\bigr)\in[0,1].
   \]
   The values at the same vertex agree across \(k\): \(\Theta_{k+1}(y/2)=\Theta_k(y)\) (\(y/2\in V_{k+1}\)). So the **layer profile** \(\Lambda(\xi)=\lim_nu_n(2^{-n}\xi)=\bar d+x_a\Theta(\xi)\) exists at every layer vertex. It is **universal**: \(x_b,x_c\) drop out. It is symmetric under the reflection fixing \(a\), and \(0<\Theta<1\).
4. **BH.4 (no layer for \(\alpha\ge\gamma\)).** \(u_n(2^{-n}\xi)\to d_a\) for every fixed \(\xi\).
5. **BH.5 (decay into the bulk, and the certificate).** Put \(\theta_k=\Theta(2^ke_b)=J_r^{-1}\int s^{-\alpha}r_\infty\,\eta(p^ks)\,ds\). Then \(\theta_k\) strictly decreases to \(0\), and \(\theta_0=1-J_w/(2J_r)=1-3\mathcal L/\sqrt{21}\). So \(x_a\theta_0\) is exactly BG.1's gap \(5/6-5\mathcal L/(2\sqrt{21})\) for \((1,-\tfrac12,0)\). Certified at \(\alpha=0.45\):
   - \(\theta_0=\Theta(1,0)=\Theta(0,1)\in[0.4624385835869,0.462438583587]\);
   - \(\theta_1=\Theta(2,0)\in[0.340381375,0.340407195]\);
   - \(\theta_2=\Theta(4,0)\in[0.241594653,0.241622693]\);
   - \(\theta_Y=\Theta(1,1)\in[0.382431398,0.382456436]\) (the cell midpoint opposite \(a\)).

   So \(\theta_0>\theta_Y>\theta_1>\theta_2>0.24\). For the data \((1,-\tfrac12,0)\), \(\lim u_n(2^{-n}\xi)=\tfrac16+\tfrac56\Theta(\xi)\). That is \(0.5520321530\ldots\) at both neighbours, \([0.450317812,0.450339330]\) at \(\xi=(2,0)\), and \([0.485359498,0.485380364]\) at \((1,1)\), against \(\bar d=0.1\overline6\) in the bulk.
6. **BH.6 (far field; DERIVED from AT and BD.4 as stated).** For \(\alpha<\gamma\), \(\rho^{-k}\theta_k\to\theta_\infty=J_r^{-1}\int_0^\infty\tau^{\gamma-1-\alpha}\,\Theta_{AT}(\psi(\tau))\,\eta(\psi(\tau))\,d\tau\in(0,\infty)\). Here \(\Theta_{AT}\) is AT.1's log-periodic function \(\lim_kg(q^ku)\), \(g(s)=s^{1-\gamma}r_\infty(s)\). So the layer decays into the bulk at exactly BG.2's outer rate \(\rho\) per doubling of the distance to the corner. The value \(\theta_\infty(0.45)\approx0.52061\) comes from floats only and is not certified.

**Proof of BH.1.** Every vertex of the corner cell \(2^{-(n-k)}V_k\) other than its three corners has all its \(V_n\)-neighbours inside the cell, and degree \(4\). The source \(x_C\) vanishes there. So \(w=(L_n+s)^{-1}x_C\) satisfies the \(V_k\) equations \((4+s)w_y=\sum w_{y'}\) in the interior and is the \(s\)-harmonic extension of its corner values. The corner values come from \(w_a=r_nx_a\) (BG.1, corner mode) and from BG.1's restriction identity \(w|_{V_{n-k}}=r_n(s)h_{n-k}(p^ks;\cdot)\) at \(\nu_b,\nu_c\in V_{n-k}\). Then use \(u_n-\bar d=D_n^{-1}\int s^{-\alpha}w\,ds\) (BG.1). \(\square\)

**Proof of BH.2.**
- *Recursion.* Apply BH.1 with \(k=1\) at level \(N\) and spectral parameter \(z\). \(\nu_b\) is the midpoint \(m_{ab}\) of the level-1 corner cell, whose far corners are the \(V_{N-1}\) neighbours \(\nu_b',\nu_c'\). With \(h_N(z;\cdot)=h_{N-1}(p(z);\cdot)\) on \(V_{N-1}\) (BG.1):
  \[
  h_N(z;\nu_b)=\alpha_1(z)x_a+\alpha_1(z)h_{N-1}(pz;\nu_b')+\tfrac2{(z+2)(z+5)}h_{N-1}(pz;\nu_c'),
  \]
  and the same with \(b\leftrightarrow c\). The base case is \(h_0(z;b)=x_b\), \(h_0(z;c)=x_c\).
- *Asymmetry.* The antisymmetric factor is \(\alpha_1-2/((z+2)(z+5))=1/(z+5)\) (sympy). Iterate from \(x_b-x_c\).
- *Symmetric part.* \((\eta x_a,\eta x_a)\) is a fixed point, because \(\alpha_1+2/((z+2)(z+5))=\beta\) and \(\eta=\alpha_1+\beta\,\eta\circ p\). The error vector \(e_N(z)=h_N(z;\nu_\cdot)-x_a\eta(z)\) satisfies \(e_N(z)=M(z)e_{N-1}(pz)\). \(M\) has nonnegative entries with row sums \(\beta(z)\), and \(\beta\le\beta(0)=3/5\) (\(\beta'=-(z^2+12z+32)/((z+2)(z+5))^2\), sympy). The base error is at most \(\|x\|_\infty(1+\eta)\le2\|x\|_\infty\).
- *Properties of \(\eta\).* The series converges uniformly because \(\prod\beta\le(3/5)^j\), so \(\eta\) is continuous. \(\eta(0)=\sum\frac25(\frac35)^j=1\). Each term is a product of positive decreasing functions (\(p\) is increasing), and the \(j=0\) term is strictly decreasing, so \(\eta\) is strictly decreasing. \(\eta\le\alpha_1+\beta=2/(2+z)\le1\) (sympy).
- *Corner resolvent.* The corner equation at \(a\) (degree \(2\)) reads \((2+s)w_a-w_{\nu_b}-w_{\nu_c}=x_a\), i.e. \(r_n(s)=x_a/((2+s)x_a-h_n(s;\nu_b)-h_n(s;\nu_c))\) (checked in exact rationals). Let \(n\to\infty\) with the contraction. The limit denominator is \(2+s-2\eta\ge s>0\).
- *Corollary.* Integrate the asymmetry against \(s^{-\alpha}r_n/D_n\). \(\prod_{j<n}(5+p^j)^{-1}\le5^{-n}\) gives \(I_n\le5^{-n}D_n\). \(\square\)

**Proof of BH.3.** In BH.1, replace \(h_{n-k}(p^ks;\nu_v)\) by \(x_a\eta(p^ks)\). The error is at most \(2\|x\|_\infty(3/5)^{n-k}\), uniformly in \(s\) (BH.2). The weights \(H_k(s;v,y)\ge0\) sum to at most \(1\), and \(s^{-\alpha}r_n/D_n\) is a probability density. This leaves \(A_n/D_n\), with \(A_n=\int s^{-\alpha}r_n\phi_k\) and \(A=\int s^{-\alpha}r_\infty\phi_k\). Since \(0\le\phi_k\le1\) and \(r_n\le r_\infty\) (AQ (3)), \(0\le A-A_n\le J_r-D_n\) and \(|A_n/D_n-A/J_r|\le2(J_r-D_n)/D_n\). For \(\alpha<\gamma\), \(D_n\uparrow J_r<\infty\) (monotone convergence, AQ (3)). Consistency across \(k\) holds because both sides are limits of the same sequence. Universality and the reflection symmetry are read off \(\phi_k\). \(0<\Theta<1\) because \(0<\phi_k<1\) on \((0,\infty)\). \(\square\)

**Proof of BH.4.** The \((3/5)^{n-k}\) replacement is as in BH.3. Then \(1-A_n/D_n=D_n^{-1}\int s^{-\alpha}r_n(1-\phi_k)\). \(\phi_k\) is continuous with \(\phi_k(0)=1\) (harmonic measure sums to \(1\), \(\eta(0)=1\)). So for each \(\varepsilon\) there is a \(\delta\) with \(1-\phi_k\le\varepsilon\) on \([0,\delta]\). The rest is at most \(D_n^{-1}\int_\delta^\infty s^{-1-\alpha}\), since \(r_n\le1/s\). For \(\alpha\ge\gamma\), \(D_n\to\infty\) (AY, AZ (1)), so \(\limsup(1-A_n/D_n)\le\varepsilon\). \(\square\)

**Proof of BH.5.**
- *Value.* \(\phi_k(s;b)=\eta(p^ks)\). Strict decrease comes from \(p^{k+1}s>p^ks\) and the strict decrease of \(\eta\). The limit is \(0\) by dominated convergence (\(\eta(p^ks)\le2/(2+5^ks)\)).
- *\(\theta_0\).* \(r\eta=r-(1-sr)/2\) (corner equation), so \(\theta_0=1-J_w/(2J_r)\). AS's \(\mathcal L=\sqrt{21}J_w/(6J_r)\) gives the rest (sympy).
- *Certificate.* We enclose \(M_\phi=\int s^{-\alpha}r_\infty(1-\phi)\) for \(\phi\in\{\eta,\eta\circ p,\eta\circ p^2,\phi_Y\}\), with \(\phi_Y=(2+2(s+4)\eta(ps))/((s+2)(s+5))\) (sympy), and set \(\Theta=1-M_\phi/J_r\).
  - Monotonicity: \(r_\infty=1/(2+s-2\eta)\) is decreasing and \(1-\phi\) is increasing. So on a cell \([s_0,s_1]\) the integrand lies between \(r(s_1)(1-\phi(s_0))\) and \(r(s_0)(1-\phi(s_1))\), times \(\int_{\rm cell}s^{-\alpha}\).
  - Node values: \(\eta\) is enclosed by its partial sum plus a tail in \([0,\prod\beta\cdot2/(2+p^Jz)]\), stopping when \(\prod\beta<10^{-30}\).
  - Grid: \(742386\) nodes, ratio \(1.00002\) on \([10^{-3},10^3]\), \(1.0005\) out to \([10^{-8},10^8]\), \(1.01\) out to \([10^{-16},10^{24}]\). Arithmetic is python-flint `arb` balls at 80 bits.
  - Head \((0,10^{-16}]\): \(r(1-\phi)\le4\). For \(\eta\circ p^k\) use \(r(s)(1-\eta(z))=r(s)(1-zr(z))/(2r(z))\le\frac12\prod_{i<k}\frac{(p^is+2)(p^is+5)}{p^is+4}\). This follows from AQ's fixed point and \(r(ps)\le1/(ps)\): \(c-e/p=(s+4)/((s+2)(s+5))\) (sympy). For \(\phi_Y\), use \(1-\phi_Y=(s(s+5)+2(s+4)(1-\eta(ps)))/((s+2)(s+5))\) and \(sr\le1\).
  - Tail \([10^{24},\infty)\): \(r(1-\phi)\le1/s\).
  - \(J_r\in[4.0013017734819,4.001301773482]\) and \(J_w\) come from Theorem AS, re-run live.
  - Cross-check by an independent route: the \(\eta\)-series enclosure of \(M_\eta\), \([2.150897988,2.150992940]\), contains AS's \(J_w/2\). \(\square\)

**Proof of BH.6.**
- *Rescaling.* Substitute \(s=5^{-k}\tau\). Then \(\rho^{-k}\theta_kJ_r=\int\tau^{-\alpha}(3/5)^kr_\infty(5^{-k}\tau)\eta(\psi_k\tau)\,d\tau=\int\tau^{\gamma-1-\alpha}g(5^{-k}\tau)\eta(\psi_k\tau)\,d\tau\), using \(5^{1-\gamma}=5/3\). Here \(5^{-k}\tau=q^k(\psi_k(\tau))\).
- *Pointwise limit.* AT's bracket \(\Theta_{AT}(s)\in g(s)\exp([-(1-\gamma)q(s)/4,\,2q(s)^{0.65}/(1-5^{-0.65})])\) and \(\Theta_{AT}\circ q=\Theta_{AT}\) give \(g(q^kv)/\Theta_{AT}(v)\to1\) uniformly for \(v\) in compacts. \(\Theta_{AT}\) is continuous, as a locally uniform limit of the continuous \(g\circ q^j\) (AT's small-\(y\) lemma, geometric increments). With \(\psi_k(\tau)\uparrow\psi(\tau)\) (BD.4), \(g(5^{-k}\tau)\to\Theta_{AT}(\psi(\tau))\).
- *Domination.* \(\sup g<\infty\): \(g\le s^{-\gamma}\le1\) for \(s\ge1\); \(g\) is continuous on compacts; for \(s\le10^{-3}\), \(\Theta_{AT}(s)=\Theta_{AT}(p^js)\) with \(p^js\in[q(10^{-3}),10^{-3}]\), plus the bracket. With \(\eta(\psi_k\tau)\le2/(2+\tau)\), the dominator \(\tau^{\gamma-1-\alpha}\cdot2\sup g/(2+\tau)\) is integrable exactly when \(0<\alpha<\gamma\).
- *Positivity.* \(\Theta_{AT}>0\) and \(\eta>0\). \(\square\)

**Evidence (floats; consistency only).**
- Dense \(L^{0.45}\) solves at levels \(1\)–\(7\), \((u_n(2^{-n}\xi)-\bar d)/x_a\), each increasing and below its certified limit:
  - \(\xi=(1,0)\): \(0.4045,0.4238,0.4365\) at \(n=5,6,7\), against \(0.4624\);
  - \((2,0)\): \(0.2692,0.2930,0.3085\), against \(0.3404\);
  - \((1,1)\): \(0.3158,0.3381,0.3526\), against \(0.3824\);
  - \((4,0)\): \(0.1871,0.2050\) at \(n=6,7\), against \(0.2416\).
- The two normalised neighbour values differ by \(3.7\cdot10^{-7}\) at level \(7\), matching the exact \((x_b-x_c)I_n/D_n\) (checked to \(10^{-10}\) at levels \(1\)–\(4\)).
- At \(0.8\), \(\xi=(1,0)\) gives \(0.756,0.811,0.851\) at \(n=5,6,7\), increasing towards \(1\) (BH.4; slow, no rate proved).
- \(\rho^{-k}\theta_k=0.462439,0.510837,0.519468,0.520482,0.520593\) at \(k=0,2,4,6,8\). The increments shrink by about \(9\) per two levels, i.e. a \(3^{-k}\) correction.

**CONJECTURE BH.7.**
- \(\rho^{-k}\theta_k-\theta_\infty=O(3^{-k})\), with \(\theta_\infty(0.45)=0.52061\ldots\).
- The rescaled layer \(\rho^{-k}\Theta(2^k\xi)\) converges for every \(\xi\) to a profile that matches BG.3's outer renormalised profile near the corner.
- Closed forms for \(\eta\), \(\Theta\) and \(\theta_\infty\) are not known.

Conjecture BG.2 (uniform convergence above \(\gamma\)) is **not** decided here (it is proved by Theorem BI below). BH.4 gives the absence of a layer at every fixed layer vertex, which is the other order of limits. No modulus of continuity uniform in \(n\) is proved.

**NOT THIS.** Theorem BH is the refinement limit of one graph Dirichlet extension near a corner of the combinatorial gasket, in rescaled coordinates. The layer is where the corner data survive when points are polar for the discrete fractional operator. It is a function on a graph, not a field, a flux or a force. No flux functional on the infinite gasket is constructed. Not thrust.

**Witness, computed by `scripts/gasket_fractional_corner_layer.py`** (about 3 min; needs sympy, mpmath, python-flint).
- *Sympy:* the level-1 harmonic weights; \(\alpha_1+\beta=2/(2+z)\); the antisymmetric factor \(1/(z+5)\); \(\beta\le3/5\); \(\eta(0)=1\); the identity between \(r_\infty=1/(2+s-2\eta)\) and AQ's fixed point; the head lemma; \(\phi_Y\); \(\theta_0=1-J_w/(2J_r)\); the BG.1 gap.
- *Exact rationals (levels \(2\)–\(4\), \(s\in\{1/2,7/5\}\)):* BH.1 on every corner cell; the asymmetry identity; the corner equation; and the `arb` contraction bound \(|h_n-x_a\eta|<2\|x\|\prod\beta\).
- *Dense meshes (levels \(1\)–\(7\), \(\alpha=0.45,0.8\)):* the corollary; monotone approach to the certified limits; BH.3's error bound at levels \(1\)–\(4\).
- *Certificate:* as in BH.5, with AS re-run.
- *Floats:* \(\rho^{-k}\theta_k\).

The script exits nonzero if any check fails.

## Theorem BI. Above \(\log3/\log5\) the fractional Dirichlet extension converges uniformly to a continuous limit (Conjecture BG.2 proved), with cell modulus \(\le K_\alpha\rho^{-N}\)

Theorem BG gave a pointwise limit \(u_\infty\) on \(V_*\) for \(\alpha>\gamma\), and Kept Failure BG.1 showed that below \(\gamma\) the convergence is not uniform. Conjecture BG.2 asked whether, above \(\gamma\), the convergence is uniform and the limit is continuous with the corner data as boundary values. This theorem proves that, with an explicit modulus on cells. The proof is a scale-by-scale oscillation recursion for the resolvent, integrated against the Stieltjes weight. It needs a uniform power bound on the rescaled corner resolvent, which is certified here.

**Setting.** As in Theorems BG and BH: Assumption 4, level \(n\), \(\alpha\in(0,1)\), corner data \(d\), \(x=d-\bar d\), \(\|x\|_\infty=\max|x_c|\), \(u_n\) the \(L_n^{\alpha}\)-Dirichlet extension, \(w_s=(L_n+s)^{-1}x_C\), \(r_n\), \(r_\infty\), \(p(s)=s(s+5)\), \(D_n=\int_0^\infty s^{-\alpha}r_n\), \(\rho=5^{\alpha}/3\), \(\gamma=\log3/\log5\). Rescaled resolvent \(\mathcal R_n(\tau)=(3/5)^nr_n(5^{-n}\tau)\) (BG.3), and \(g(s)=s^{1-\gamma}r_\infty(s)\) (AT). A **level-\(i\) cell** of \(V_n\) is one of the \(3^i\) copies \(2^{-i}V_{n-i}\) (closed triangles in the plane; \(K\) is the gasket, \(V_*=\bigcup V_m\) is dense in \(K\)). \(\omega_i(s)\) is the largest difference of \(w_s\) between two corners of one level-\(i\) cell. \(m(t)=\min(1,t)\), and
\[
\Phi_N(\tau)=\sum_{l<N}(3/5)^{N-1-l}m(5^{-l-1}\tau)+\sum_{l\ge N}m(5^{-l-1}\tau),\qquad F_N=\min\bigl(2,\;10(3/5)^N+5\Phi_N\bigr).
\]
**Assumptions used as stated:** BG.1 (charge form \(u_n-\bar d=D_n^{-1}\int s^{-\alpha}w_s\,ds\); corner mode \(w_s(c)=r_n(s)x_c\)), BG.3 (pointwise limit \(u_\infty\) on \(V_*\) for \(\alpha>\gamma\)), AQ (1)–(3) (decimation, Schur complement \(c(L_k+p)-eP_C\), fixed point \(r_\infty(s)=r_\infty(ps)/(c-er_\infty(ps))\), \(r_n\uparrow r_\infty\)), and BH.2 (\(r_\infty=1/(2+s-2\eta)\), \(\eta\le2/(2+s)\), and the level-1 weights). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BI.**

1. **BI.1 (cell recursion for the resolvent; every \(\alpha\), \(n\), \(s>0\)).** \(\omega_0\le2r_n\|x\|_\infty\), and for \(i<n\)
   \[
   \omega_{i+1}(s)\le\tfrac35\,\omega_i(s)+\kappa(z_i)\,r_n(s)\|x\|_\infty,\qquad z_i=p^{\,n-i-1}(s),\quad\kappa(z)=\frac z{z+2}.
   \]
   For every \(N\le n\) and every level-\(N\) cell \(\Delta\), \(\operatorname{osc}_\Delta w_s\le\min\bigl(2r_n\|x\|_\infty,\;2\sum_{i=N}^n\omega_i\bigr)\).
2. **BI.2 (orbit lemma).** If \(5^ks\le1\), then \(p^k(s)\le e^{1/2}5^ks\). Hence \(\kappa(z_l)\le m(5^{-l-1}\tau)\) with \(\tau=5^ns\), and \(\operatorname{osc}_\Delta w_s\le r_n(s)\|x\|_\infty F_N(5^ns)\).
3. **BI.3 (cell oscillation of \(u_n\); every \(\alpha\), \(n\ge\max(N,1)\)).**
   \[
   \operatorname{osc}_\Delta u_n\le\|x\|_\infty\,\frac{\rho^n}{D_n}\int_0^\infty\tau^{-\alpha}\mathcal R_n(\tau)F_N(\tau)\,d\tau .
   \]
4. **BI.4 (two-sided bounds on \(\mathcal R_n\), \(n\ge1\)).**
   - *Floor:* \(\mathcal R_n(\tau)\ge1/(3+\tau)\). So \(\rho^n/D_n\le S_\alpha=\sin(\pi\alpha)/(3^{-\alpha}\pi)\).
   - *Ceiling:* \(\mathcal R_n(\tau)\le\min(\tfrac13,\,G\tau^{\gamma-1})\) with \(G=0.6031\). Certified: \(\sup_{s>0}g(s)\le0.6031\). On \([10^{-5},1]\) arb cells give \(g\le0.601271031\); a bootstrap along the \(p\)-orbit carries this below \(10^{-5}\) with factor \(\le1.0008\); and \(g\le s^{-\gamma}(2+s)/(4+s)\le3/5\) on \([1,\infty)\). (For comparison, \(\inf g\ge0.4244\) on \([10^{-5},1]\), and \(\langle\Theta\rangle=0.60095\ldots\), AU.)
5. **BI.5 (uniform modulus for \(\alpha>\gamma\)).** For \(\alpha\in(\gamma,1)\), every \(n\ge\max(N,1)\) and every level-\(N\) cell \(\Delta\) of \(V_n\),
   \[
   \operatorname{osc}_\Delta u_n\le\|x\|_\infty\,\varepsilon_N(\alpha),\qquad
   \varepsilon_N(\alpha)=S_\alpha\int_0^\infty\tau^{-\alpha}\min(\tfrac13,G\tau^{\gamma-1})F_N(\tau)\,d\tau\le K_\alpha\rho^{-N},
   \]
   uniformly in \(n\). Explicitly \(K_\alpha=S_\alpha\bigl(10I_\alpha+GJ_\alpha\bigr)\), with \(I_\alpha=\int\tau^{-\alpha}\min(\tfrac13,G\tau^{\gamma-1})\) and \(J_\alpha=\int_0^\infty t^{\gamma-1-\alpha}\min(2,21.25\,t^{1-\gamma})\,dt\) (both in closed form). Certified (python-flint `arb`, outward):
   - \(K_{0.75}\le70.6512\), \(K_{0.8}\le47.109\), \(K_{0.9}\le39.005\);
   - \(\varepsilon_N(0.8)\le5.90016,\ 4.78623,\ 2.19931,\ 0.347709,\ 0.0526623,\ 0.00796134,\ 1.81925\cdot10^{-4}\) at \(N=0,5,10,20,30,40,60\) (so \(\rho^N\varepsilon_N(0.8)\le15.2468\)).
6. **BI.6 (Conjecture BG.2 holds).** For every \(\alpha\in(\gamma,1)\):
   - \(u_\infty\) is uniformly continuous on \(V_*\): \(|u_\infty(y)-u_\infty(y')|\le2\|x\|_\infty\varepsilon_N\) whenever \(y,y'\) lie in two level-\(N\) cells that meet. So it extends to a continuous \(\bar u\) on \(K\), with \(\bar u=d\) on \(C\).
   - \(u_n\to\bar u\) uniformly: \(\sup_{V_n}|u_n-\bar u|\le2\|x\|_\infty\varepsilon_N+\max_{V_N}|u_n-u_\infty|\) for every \(N\le n\), and the last term tends to \(0\) by BG.3.
   - *Derived:* since \(\rho^{-N}=((3/5)^N)^{(\alpha-\gamma)/(1-\gamma)}=(2^{-N})^{(\alpha\log5-\log3)/\log2}\), \(\bar u\) is Hölder on cells with exponent \((\alpha-\gamma)/(1-\gamma)\) at the resistance scale \((3/5)^N\). The exponent tends to \(1\) as \(\alpha\to1\), the harmonic Lipschitz rate \((3/5)^N\), and to \(0\) as \(\alpha\downarrow\gamma\).

**Proof of BI.1.**
- *Decimation inside a cell.* A level-\(i\) cell is a copy of \(V_{n-i}\). Its non-corner vertices have all their \(V_n\)-neighbours in the cell, and \(x_C\) vanishes there. So \(w_s\) is \(s\)-harmonic inside the cell (BH.1). AQ's Schur complement \(c(L_k+p(s))-eP_C\) acts as \(c(L_k+p(s))\) at non-corner rows. Eliminating down to the cell's \(V_1\), the three midpoints are therefore the level-1 \(z\)-harmonic extension of the cell corners \(f_1,f_2,f_3\), with \(z=z_i=p^{n-i-1}(s)\). By BH's weights, \(m_{12}=\alpha_1(f_1+f_2)+\delta f_3\), with \(\alpha_1=(z+4)/((z+2)(z+5))\), \(\delta=2/((z+2)(z+5))\), and \(2\alpha_1+\delta=2/(z+2)\) (sympy).
- *One step.* \(m_{12}-f_1=\alpha_1(f_2-f_1)+\delta(f_3-f_1)-\kappa(z)f_1\), so \(|m_{12}-f_1|\le\beta\,\omega_i+\kappa\,|f_1|\), with \(\beta=\alpha_1+\delta\le3/5\) (BH). Also \(m_{12}-m_{13}=(f_2-f_3)/(z+5)\), so \(|m_{12}-m_{13}|\le\omega_i/5\). The three subcells have corners \(\{f_k,m_{kl},m_{km}\}\).
- *Maximum principle.* \((\deg y+s)w_y=\sum_{y'\sim y}w_{y'}\) off \(C\), so \(\max|w_s|\) is attained on \(C\), where it equals \(r_n\|x\|_\infty\) (BG.1). This gives \(|f_k|\le r_n\|x\|_\infty\), \(\omega_0\le r_n\max|x_c-x_{c'}|\le2r_n\|x\|_\infty\), and the cap \(2r_n\|x\|_\infty\).
- *Chain.* Take \(y\in\Delta\), nested cells \(\Delta=\Delta_N\supset\cdots\supset\Delta_n\ni y\), and \(a_i\) the corner that \(\Delta_{i+1}\) shares with \(\Delta_i\). The path \(y\to a_{n-1}\to\cdots\to a_N\to v\), with \(v\) a corner of \(\Delta\), uses one corner pair of \(\Delta_i\) at each level. So \(|w_s(y)-w_s(v)|\le\sum_{i=N}^n\omega_i\). \(\square\)

**Proof of BI.2 and BI.3.**
- *Orbit lemma.* \(a_{j+1}=5a_j(1+a_j/5)\) gives \(a_j\le5^js\exp(\sum_{i<j}a_i/5)\). Inductively, with \(a_i\le e^{1/2}5^is\), the exponent is at most \(e^{1/2}5^js/4\le e^{1/2}/4<1/2\) when \(5^js\le1\). Then \(\kappa(z)\le\min(1,z/2)\le\min(1,(e^{1/2}/2)5^{n-l-1}s)\le m(5^{-l-1}\tau)\).
- *Summation.* Unrolling BI.1 gives \(\omega_i\le r_n\|x\|_\infty[2(3/5)^i+\sum_{l<i}(3/5)^{i-1-l}m(5^{-l-1}\tau)]\). Summing over \(i\ge N\) gives \(\sum_{i\ge N}\omega_i\le r_n\|x\|_\infty[5(3/5)^N+\tfrac52\Phi_N(\tau)]\), so \(\operatorname{osc}_\Delta w_s\le r_n\|x\|_\infty F_N(\tau)\).
- *Integration.* Integrate \(\operatorname{osc}_\Delta u_n\le D_n^{-1}\int s^{-\alpha}\operatorname{osc}_\Delta w_s\) and substitute \(s=5^{-n}\tau\), using \(D_n=\rho^n\int\tau^{-\alpha}\mathcal R_n\). \(\square\)

**Proof of BI.4.**
- *Floor.* Use the Dirichlet principle \(\tfrac12v^{T}A^{-1}v\ge\max_t\bigl(t\,v^{T}h-\tfrac12t^2h^{T}Ah\bigr)\), with \(v=e_a-e_b\), \(A=L_n+s\), and \(h\) the harmonic function with corner values \((\tfrac12,-\tfrac12,0)\). Then \(v^{T}h=1\), \(h^{T}L_nh=\tfrac32(3/5)^n\), \(|h|\le\tfrac12\) and \(|V_n|=\tfrac32(3^n+1)\) (exact rationals, levels 1–4). This gives \(r_n(s)\ge1/(3(3/5)^n+\tfrac34s(3^n+1))\), so \(\mathcal R_n\ge1/(3+\tfrac34(1+3^{-n})\tau)\ge1/(3+\tau)\). Finally \(\int\tau^{-\alpha}/(3+\tau)=3^{-\alpha}\pi/\sin\pi\alpha\).
- *Ceiling.* \(\mathcal R_n\le\tfrac13\) because \(r_n\le r_n(0)=\tfrac13(5/3)^n\). Also \(\mathcal R_n(\tau)\le(3/5)^nr_\infty(5^{-n}\tau)=\tau^{\gamma-1}g(5^{-n}\tau)\), because \(5^{1-\gamma}=5/3\).
- *Bound on \(g\) for \(s\ge1\).* \(\eta\le2/(2+s)\) gives \(r_\infty\le(2+s)/(s(s+4))\), so \(g\le f(s)=s^{-\gamma}(2+s)/(4+s)\). The logarithmic derivative of \(f\) has the sign of \(-(\gamma(s^2+6s+8)-2s)<0\) (negative discriminant), so \(g\le f(1)=3/5\).
- *\([10^{-5},1]\).* On each of \(5763\) cells of ratio \(1.002\), \(g\le s_1^{1-\gamma}\,r_\infty^{\uparrow}(s_0)\), with BH's arb enclosure of \(\eta\).
- *Below \(10^{-5}\) (bootstrap).* By AQ's fixed point, \(g(t)=g(pt)\,Q(t)\) with \(Q=(3/5)(1+t/5)^{\gamma-1}/(c-e\,r_\infty(pt))\). If \(g\le G\) along the orbit above \(t\), then \(e\,r_\infty(pt)\le t(5t)^{\gamma-1}G=\tfrac35Gt^{\gamma}\), and \(c\ge\tfrac35-\tfrac8{25}t\) (sympy: the difference is \((41s^2+8s^3)/(25(s+2)(s+5))\)). So \(Q(t)\le1/(1-\tfrac8{15}t-Gt^{\gamma})\). Every \(s<\varepsilon_0=10^{-5}\) is \(q^k(s_0)\) for some \(s_0\in[\varepsilon_0,p(\varepsilon_0))\), with orbit points \(s_j\le5^{-j}p(\varepsilon_0)\). Hence \(g(s)\le G_0\exp\bigl(\tfrac4{15}p(\varepsilon_0)+Gp(\varepsilon_0)^{\gamma}\bigr)\le0.6016995<G\), by induction along the orbit. \(\square\)

**Proof of BI.5.**
- *Modulus.* Combine BI.3 and BI.4.
- *Rate.* Put \(t=5^{-N}\tau\). Then \(\sum_{j\ge0}(3/5)^jm(5^jt)\le4\min(1,t^{1-\gamma})\): split at \(5^Jt\ge1\), and use \(3^J\le3t^{-\gamma}\) and \((3/5)^J\le t^{1-\gamma}\). Also \(\sum_{j\ge1}m(5^{-j}t)\le t/4\) for \(t\le1\). So \(\min(2,5\Phi_N)\le\min(2,21.25\,t^{1-\gamma})\).
- *Constants.* With the weight bounded by \(G\tau^{\gamma-1-\alpha}\) and \(5^{N(\gamma-\alpha)}=\rho^{-N}\), we get \(\varepsilon_N\le S_\alpha(10(3/5)^NI_\alpha+G\rho^{-N}J_\alpha)\), and \((3/5)^N\le\rho^{-N}\). \(J_\alpha<\infty\) exactly when \(\alpha>\gamma\) (tail \(t^{\gamma-1-\alpha}\)).
- *Certificate.* \(F_N\) is nondecreasing, so on each geometric cell (ratio \(1.02\), \([10^{-40},10^{140}]\)) we bound by \(F_N\) at the right end times the exact integral of the weight. The head is \(\le F_N(10^{-40})\int_0^{10^{-40}}\tau^{-\alpha}/3\), and the tail is \(\le2G\,T^{\gamma-\alpha}/(\alpha-\gamma)\). \(\square\)

**Proof of BI.6.** Apply BI.5 to \(u_m\) and let \(m\to\infty\) (BG.3) on \(V_*\cap\Delta\). This gives \(|u_\infty(y)-u_\infty(v)|\le\|x\|_\infty\varepsilon_N\) for \(v\) a corner of \(\Delta\). If two level-\(N\) cells meet, they share a vertex of \(V_N\), so the bound doubles. There are finitely many level-\(N\) cells, and disjoint ones are at positive distance \(\delta_N\). So \(|y-y'|<\delta_N\) forces \(|u_\infty(y)-u_\infty(y')|\le2\|x\|_\infty\varepsilon_N\): \(u_\infty\) is uniformly continuous on the dense set \(V_*\) and extends to \(\bar u\in C(K)\). Since \(u_n(c)=d_c\), \(\bar u=d\) on \(C\). For \(y\in V_n\) in a level-\(N\) cell with corner \(v\): \(|u_n(y)-\bar u(y)|\le|u_n(y)-u_n(v)|+|u_n(v)-u_\infty(v)|+|\bar u(v)-\bar u(y)|\). \(\square\)

**Why this needs \(\alpha>\gamma\).** The only place \(\alpha\) enters is the weight \(\tau^{-\alpha}\mathcal R_n(\tau)\sim\tau^{\gamma-1-\alpha}\) at large \(\tau\) (small cells, \(s\gg5^{-n}\)), against \(F_N\to2\). For \(\alpha\le\gamma\) it is not integrable, and the bound is no better than the trivial \(2\|x\|_\infty\). Kept Failure BG.1 and Theorem BH show that the failure is real: a layer of width \(O(2^{-n})\) keeps the gap \(x_a\theta_0\).

**Evidence (floats; consistency only).**
- Dense \(L^{0.8}\) solves for the data \((1,-\tfrac12,0)\). The maximal oscillation on level-\(N\) cells, \(N=0,\ldots,7\), at level \(7\) is \(1.5000,0.8155,0.6128,0.4678,0.3523,0.2585,0.1829,0.1240\). Each value is below \(\|x\|\min(2,\varepsilon_N)\) (the bound is a modulus, not sharp at small \(N\)). Also \(\max_{V_6}|u_7-u_6|=0.0263\) and \(\max_{V_5}|u_6-u_5|=0.0299\).
- At \(0.45\), the finest-cell oscillation stays at \(0.48\) (the layer).
- Corner deficit for symmetric data (\(x_b=x_c\), so \(h_N=x_a\eta_N\) exactly): \(\delta_N=(d_a-u_n(\nu_N))/x_a\) at \(\alpha=0.8\), with \(\nu_N\) the level-\(N\) neighbour of \(a\). The rescaled value \(\rho^N\delta_N\) is \(1.13674,1.13518,1.13427,1.13268,1.12988\) at \(N=3,6,9,12,15\) (\(n=40\)), and \(1.13596,1.13587\) at \(N=5,10\) (\(n=60\)). So the corner modulus is of order exactly \(\rho^{-N}\).

**CONJECTURE BI.7 (floats only; the first two items are PROVED by Theorem BJ below, with \(\liminf,\limsup\rho^N\delta_N(0.8)\in[1.135890346,1.135938783]\); the third, a rate, stays open; the text records it as first stated).**
- The order \(\rho^{-N}\) in BI.5 is sharp: \(\rho^N(d_a-\bar u(\nu_N))/x_a\to c_{\rm corner}(\alpha)\), with \(c_{\rm corner}(0.8)\approx1.1359\).
- \(\bar u\) is not Hölder at the resistance scale with any exponent above \((\alpha-\gamma)/(1-\gamma)\).
- A rate for \(\max_{V_N}|u_n-u_\infty|\), hence for the uniform convergence, is not proved. Neither is a continuum equation for \(\bar u\), nor closed forms for \(W,U,\kappa\).

**NOT THIS.** Theorem BI says the graph Dirichlet extension of a discrete fractional operator on the combinatorial gasket has a continuous uniform limit above \(\gamma\), where points are not polar. It is a function on the gasket, not a field, a flux or a force. No flux functional on the infinite gasket is constructed. Not thrust.

**Witness, computed by `scripts/gasket_fractional_uniform_limit.py`** (about 5 min; needs sympy, mpmath, python-flint).
- *Sympy:* the level-1 total and killed mass, \(m_{kl}-m_{km}=(f_l-f_m)/(z+5)\); \(c\ge\tfrac35-\tfrac8{25}s\); the \(g\)-ratio identity; the \([1,\infty)\) bound; the Dirichlet-principle algebra and \(\int\tau^{-\alpha}/(3+\tau)\); the orbit-lemma constants.
- *Exact rationals (levels \(2\)–\(4\), \(s\in\{1/100,1/2,7/5\}\)):* \(\max|w_s|=r_n\|x\|\); the midpoint decimation identity on every level-\(i\) cell; the recursion of BI.1; the chain bound and \(\operatorname{osc}\le r_n\|x\|F_N(5^ns)\) on every level-\(N\) cell. The orbit lemma on exact orbits. \(r_n(s)\ge1/(3(3/5)^n+\tfrac34s(3^n+1))\) at levels \(1\)–\(6\). The energy \(\tfrac32(3/5)^n\) of the test function.
- *Certificates (`arb`):* \(\sup g\); the bootstrap; \(\varepsilon_N(0.8)\); \(K_{0.75},K_{0.8},K_{0.9}\).
- *Floats:* dense oscillations at levels \(5\)–\(7\) (\(\alpha=0.8,0.45\)); the corner deficit at \(n=40,60\).

The script exits nonzero if any check fails.

## Theorem BJ. The cell modulus \(\rho^{-N}\) of Theorem BI is sharp (Conjecture BI.7 proved): the corner deficit is \(\asymp\rho^{-N}\), with \(\liminf,\limsup\rho^N\delta_N(0.8)\in[1.135890,1.135939]\)

Theorem BI bounded the oscillation of \(\bar u\) on level-\(N\) cells by \(\|x\|_\infty K_\alpha\rho^{-N}\) for \(\alpha>\gamma\), and Conjecture BI.7 asked, from floats only, whether the order \(\rho^{-N}\) is sharp. This theorem proves that it is. The proof is an exact reciprocal identity that pins the rescaled corner resolvent between two multiples of \(\tau^{\gamma-1}\), plus a certified enclosure of \(g=s^{1-\gamma}r_\infty\) near \(0\) that is tight to \(3\cdot10^{-5}\) relative. As a by-product, the corner constant at \(0.8\) is certified to about five digits.

**Setting.** As in Theorems BG–BI: Assumption 4, \(\alpha\in(\gamma,1)\), corner data \(d\), \(x=d-\bar d\), \(u_n\) the \(L_n^{\alpha}\)-Dirichlet extension, \(\bar u\) its uniform limit (BI.6), \(r_n\), \(r_\infty=1/(2+s-2\eta)\) (BH.2), \(p(s)=s(s+5)\), \(\beta(s)=(s+6)/((s+2)(s+5))\) (\(=\) AQ's \(c(s)\)), \(e(s)=2s/(2+s)\), \(\rho=5^{\alpha}/3\), \(g(s)=s^{1-\gamma}r_\infty(s)\), \(\mathcal R_n(\tau)=(3/5)^nr_n(5^{-n}\tau)\to\mathcal R(\tau)\) (BD.4, BG.3), \(\psi_K(t)=p^K(5^{-K}t)\uparrow\psi(t)\). The neighbours of corner \(a\) in \(V_N\) are \(\nu_b^N=(2^{-N},0)\) and \(\nu_c^N=(0,2^{-N})\). The **corner deficit** is
\[
\delta_N(\alpha)=\frac{d_a-\tfrac12\bigl(\bar u(\nu_b^N)+\bar u(\nu_c^N)\bigr)}{x_a}\qquad(x_a\neq0).
\]
For symmetric data (\(x_b=x_c\)) the two neighbours agree, and \(\delta_N\) is BI.7's \((d_a-\bar u(\nu_N))/x_a\). **Assumptions used as stated:** AQ (1)–(3) (the step \(r_{k+1}(s)=r_k(ps)/(c-er_k(ps))\), \(r_0(z)=1/(3+z)\), \(r_n\uparrow r_\infty\)), BH.2 (\(\eta\), its recursion, \(0<\eta\le2/(2+z)\), \(\eta\) strictly decreasing), BG.3 (the limit \(u_\infty(x)=\bar d+\frac{2\sin\pi\alpha}{\pi}P_m(x)/A_\alpha\) on \(V_m\)), BD.4 (\(\psi_k\uparrow\psi\)), BI.4 (\(\sup_{(0,1]}g\le0.6031\), re-run live), AZ (3) (\(A_{0.8}\in[2.39640970397922052703,2.39640970397922052705]\)), and, in the derived part of BJ.6 only, AT.1 as used in BH.6 (\(g(5^{-m}t)\to\Theta_{AT}(\psi(t))\)). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BJ.**

1. **BJ.1 (exact reciprocal identity; every \(n\ge0\), \(s>0\)).**
   \[
   \frac1{r_n(s)}-\frac1{r_\infty(s)}=\prod_{i<n}\beta(p^is)\,\bigl(1+2\eta(p^ns)\bigr)\in\bigl(0,\,3(3/5)^n\bigr].
   \]
2. **BJ.2 (two-sided power bracket for the rescaled resolvent).** For every \(n\) and \(\tau>0\), with \(s=5^{-n}\tau\),
   \[
   \frac1{\tau^{1-\gamma}/g(s)+3}\ \le\ \mathcal R_n(\tau)\ \le\ \tau^{\gamma-1}g(s).
   \]
   Hence the limit satisfies \(1/(\tau^{1-\gamma}/g_{\rm lo}+3)\le\mathcal R(\tau)\le\min(\tfrac13,g_{\rm hi}\tau^{\gamma-1})\), and for every \(N\ge0\)
   \[
   \frac1{t^{1-\gamma}/g_{\rm lo}+3(3/5)^N}\ \le\ (5/3)^N\mathcal R(5^Nt)\ \le\ \min\bigl((5/3)^N/3,\ g_{\rm hi}t^{\gamma-1}\bigr),
   \]
   where \(g_{\rm lo}\le\liminf_{s\to0}g\) and \(g_{\rm hi}\ge\limsup_{s\to0}g\).
3. **BJ.3 (certified: \(g\) near \(0\)).** On \((0,p(10^{-12})]\), \(g\in[0.600945014,0.600963487]\). So \(g_{\rm lo}=0.600945014\) and \(g_{\rm hi}=0.600963487\) are admissible. (AT.1's log-periodic \(\Theta_{AT}=\lim g\circ q^k\) takes its values in this interval, and AU's certified mean \(\langle\Theta\rangle=0.6009542511\ldots\) lies inside it.)
4. **BJ.4 (corner-deficit formula; data-free).** For every \(N\ge0\),
   \[
   \delta_N=\kappa_\alpha\int_0^\infty\tau^{-\alpha}\mathcal R(\tau)\bigl(1-\eta_N(\psi(5^{-N}\tau))\bigr)\,d\tau,\qquad \kappa_\alpha=\Bigl(\int_0^\infty\tau^{-\alpha}\mathcal R\Bigr)^{-1}=\frac{2\sin\pi\alpha}{\pi A_\alpha},
   \]
   with \(\eta_N(z)=\eta(z)-\prod_{i<N}\beta(p^iz)\bigl(\tfrac12+\eta(p^Nz)\bigr)\). So \(0\le(1-\eta_N)-(1-\eta)\le\tfrac32(3/5)^N\), and \(\delta_N\) does not depend on the data. Each neighbour separately: \(|\bar u(\nu_b^N)-\bar u(\nu_c^N)|\le|x_b-x_c|5^{-N}=o(\rho^{-N})\) (BH.2's asymmetry; \(5^{-1}<\rho^{-1}\)).
5. **BJ.5 (sharpness: Conjecture BI.7, first two items, PROVED).** For every \(\alpha\in(\gamma,1)\) and every \(N\ge0\),
   \[
   \rho^N\delta_N\ \ge\ c_-(\alpha)=\kappa_\alpha\int_0^\infty\frac{t^{-\alpha}\bigl(1-\eta(\psi_K(t))\bigr)}{t^{1-\gamma}/g_{\rm lo}+3}\,dt\ >\ 0 .
   \]
   Certified (python-flint `arb`): \(c_-(0.75)\ge0.7282\), \(c_-(0.8)\ge0.6208\) (\(\ge0.6086\) without AZ's \(A_{0.8}\)), \(c_-(0.9)\ge0.3630\). Consequences:
   - the level-\(N\) corner cell at any corner \(v\) with \(x_v\ne0\) (one exists for nonconstant data) has \(\operatorname{osc}\bar u\ge|x_v|c_-(\alpha)\rho^{-N}\). With BI.5, \(|x_v|c_-\rho^{-N}\le\max_{\Delta\in\text{level }N}\operatorname{osc}_\Delta\bar u\le\|x\|_\infty K_\alpha\rho^{-N}\). So \(\rho^{-N}\) is the exact order of the cell modulus (\(K_{0.8}/c_-(0.8)\le75.88\));
   - \(\bar u\) is not Hölder on cells at the resistance scale \((3/5)^N\) with any exponent above \((\alpha-\gamma)/(1-\gamma)\). In the plane, since \(|v-\nu^N|=2^{-N}\), \(\bar u\) is not Hölder at \(v\) with any exponent above \((\alpha\log5-\log3)/\log2\). BI.6's derived exponent is optimal.
6. **BJ.6 (the corner constant).** Put \(M_\alpha=\kappa_\alpha\int_0^\infty t^{\gamma-1-\alpha}\bigl(1-\eta(\psi(t))\bigr)\,dt<\infty\). Then
   \[
   g_{\rm lo}M_\alpha\le\liminf_N\rho^N\delta_N\le\limsup_N\rho^N\delta_N\le g_{\rm hi}M_\alpha,
   \]
   and, for every \(N\), \(\rho^N\delta_N\le g_{\rm hi}M_\alpha+\tfrac32\kappa_\alpha(3/5)^N\int t^{-\alpha}\min((5/3)^N/3,g_{\rm hi}t^{\gamma-1})\,dt\), where the last term is \(O((3/5)^{N(1-\alpha)/(1-\gamma)})\) in closed form. Certified at \(0.8\): \(\int t^{\gamma-1-\alpha}(1-\eta\circ\psi)\in[12.104976020,12.105120109]\), \(\kappa_{0.8}=0.156148471980\ldots\), so
   \[
   1.135890346\le\liminf_N\rho^N\delta_N(0.8)\le\limsup_N\rho^N\delta_N(0.8)\le1.135938783,
   \]
   and \(\rho^N\delta_N(0.8)\le2.666059,\,1.441963,\,1.197144,\,1.138387,\,1.135943\) at \(N=0,5,10,20,40\). BI.7's float \(c_{\rm corner}(0.8)\approx1.1359\) is confirmed. *Derived (from AT.1 as stated, as in BH.6):* the limit exists,
   \[
   c_{\rm corner}(\alpha)=\lim_N\rho^N\delta_N=\kappa_\alpha\int_0^\infty t^{\gamma-1-\alpha}\,\Theta_{AT}(\psi(t))\,\bigl(1-\eta(\psi(t))\bigr)\,dt .
   \]

**Proof of BJ.1.** AQ's step is affine in the reciprocal: \(1/r_{k+1}(s)=\beta(s)/r_k(ps)-e(s)\), because \(c=\beta\) (sympy). \(r_\infty\) satisfies the same fixed-point step (AQ (3), BH.2). So \(D_k=1/r_k-1/r_\infty\) obeys \(D_{k+1}(s)=\beta(s)D_k(ps)\). The base is \(D_0(z)=(3+z)-(2+z-2\eta(z))=1+2\eta(z)\in(1,3]\). Since \(\beta\) decreases from \(\beta(0)=3/5\) (sympy), the product is at most \((3/5)^n\). The step is checked in exact rationals at levels \(1\)–\(4\), and the identity in `arb` at levels \(0\)–\(4\). \(\square\)

**Proof of BJ.2.** Since \(5^{1-\gamma}=5/3\), \((5/3)^n/r_\infty(5^{-n}\tau)=\tau^{1-\gamma}/g(5^{-n}\tau)\). BJ.1 gives \((5/3)^n/r_n\le(5/3)^n/r_\infty+3\), which is the floor. The ceiling is \(r_n\le r_\infty\) (AQ (3)). Let \(n\to\infty\); then \(5^{-n}\tau\to0\), and \(\mathcal R\le\tfrac13\) is BG.3. For the scaled form put \(\tau=5^Nt\), so that \(\tau^{1-\gamma}=(5/3)^Nt^{1-\gamma}\). \(\square\)

**Proof of BJ.3.**
- *One period.* \(r_\infty\) is decreasing and \(s^{1-\gamma}\) increasing. So on each of the \(160945\) cells of ratio \(1.00001\) covering \([\varepsilon_0,p(\varepsilon_0)]\), \(\varepsilon_0=10^{-12}\), \(g\) lies between \(s_0^{1-\gamma}r_\infty(s_1)\) and \(s_1^{1-\gamma}r_\infty(s_0)\). We use BH's `arb` enclosure of \(\eta\) at 96 bits. This gives \(g\in[0.60094501459464,0.60096348281152]\) there.
- *Bootstrap.* Every \(t<\varepsilon_0\) is \(q^k(s_0)\) with \(s_0\in[\varepsilon_0,p(\varepsilon_0))\) and \(t_j=q^j(s_0)\le5^{-j}p(\varepsilon_0)\). Write \(g(t_j)=g(t_{j-1})Q(t_j)\), with \(Q=(3/5)(1+t/5)^{\gamma-1}/(c-er_\infty(pt))\).
  - Lower: \(c\le3/5\) and \(er\ge0\) give \(Q\ge e^{-(1-\gamma)t/5}\).
  - Upper: BI.4's \(c\ge\tfrac35-\tfrac8{25}t\) and \(er_\infty(pt)\le\tfrac35Gt^{\gamma}\), with \(G=0.6031\), give \(Q\le1/(1-\tfrac8{15}t-Gt^\gamma)\).
  - Summing, with \(5^{\gamma}=3\), the factors are \(e^{-(1-\gamma)p(\varepsilon_0)/20}\) and \(\exp\bigl((\tfrac2{15}p(\varepsilon_0)+\tfrac G2p(\varepsilon_0)^\gamma)/(1-x_1)\bigr)\), both within \(10^{-8}\) of \(1\). \(\square\)

**Proof of BJ.4.**
- *Neighbour mean.* By BH.2's recursion, \(m_N=\tfrac12(h_N(z;\nu_b)+h_N(z;\nu_c))\) obeys \(m_N(z)=\alpha_1x_a+\beta\,m_{N-1}(pz)\) (\(\alpha_1+\delta=\beta\), sympy), with \(m_0=(x_b+x_c)/2=-x_a/2\). So \(m_N=x_a\eta_N\), and \(\eta_N-\eta=\prod\beta\cdot(\eta_0-\eta\circ p^N)\) with \(\eta_0=-\tfrac12\). This is checked in exact rationals at levels \(2\)–\(4\) for the data \((3,-1,-2)\).
- *Limit.* BG.3 with \(m=N\) gives \(\bar u(\nu)=\bar d+\kappa_\alpha\int\tau^{-\alpha}\mathcal R\,h_N(\psi(5^{-N}\tau);\nu)\). Since \(\int\tau^{-\alpha}\mathcal R=\lim\rho^{-n}D_n=\pi A_\alpha/(2\sin\pi\alpha)\) (BI.4's dominator, AZ (1)), \(d_a=\bar d+x_a\kappa_\alpha\int\tau^{-\alpha}\mathcal R\).
- *Each neighbour.* \(\prod(5+p^iz)^{-1}\le5^{-N}\). \(\square\)

**Proof of BJ.5.**
- *Reduction.* In BJ.4 use \(1-\eta_N\ge1-\eta\), \(\psi\ge\psi_K\) and \(\eta\) decreasing. Substitute \(\tau=5^Nt\); since \(5^{N(1-\alpha)}(3/5)^N=\rho^{-N}\),
  \[
  \rho^N\delta_N\ge\kappa_\alpha\int t^{-\alpha}(5/3)^N\mathcal R(5^Nt)\bigl(1-\eta(\psi_Kt)\bigr)\,dt .
  \]
  Then apply BJ.2's scaled floor and \(3(3/5)^N\le3\).
- *Positivity.* The integrand is positive (\(\eta<1\) on \((0,\infty)\)), and \(g_{\rm lo}>0\). So \(c_->0\) for every \(\alpha\in(\gamma,1)\), and the numbers only quantify it.
- *Certificate.* \(\kappa_\alpha\ge1/\int\tau^{-\alpha}\min(\tfrac13,g_{\rm hi}\tau^{\gamma-1})\) (closed form, BJ.2), or \(\kappa_{0.8}\) from AZ's \(A_{0.8}\). On cells of ratio \(1.002\) on \([10^{-13},10^4]\), use \((1-\eta(\psi_{12}(t_0)))f(t_1)(t_1-t_0)\), with \(f\) the decreasing weight. The head is dropped, and the tail \([10^4,\infty)\) is bounded below in closed form.
- *Hölder.* \(a\) and \(\nu^N_{b,c}\) lie in the level-\(N\) corner cell, and \(\rho^{-N}=((3/5)^N)^{(\alpha-\gamma)/(1-\gamma)}=(2^{-N})^{(\alpha\log5-\log3)/\log2}\). \(\square\)

**Proof of BJ.6.**
- *Lower.* By BJ.2's scaled floor, the integrand in the reduction increases in \(N\) to \(g_{\rm lo}t^{\gamma-1-\alpha}(1-\eta\circ\psi)\). Apply monotone convergence, with \(1-\eta_N\ge1-\eta\).
- *Upper.* Use \((5/3)^N\mathcal R(5^Nt)\le\min((5/3)^N/3,g_{\rm hi}t^{\gamma-1})\) and \(1-\eta_N\le1-\eta+\tfrac32(3/5)^N\). The closed form makes the extra term \(\propto(3/5)^N(5/3)^{N(\alpha-\gamma)/(1-\gamma)}\).
- *Finiteness of \(M_\alpha\).* At \(\infty\), \(t^{\gamma-1-\alpha}\) is integrable because \(\alpha>\gamma\). Near \(0\), \(1-\eta(z)=(z^{1-\gamma}/g(z)-z)/2\le z^{1-\gamma}/(2g_{\rm lo})\) and \(\psi(t)\le te^{t/10}\), so the integrand is \(O(t^{-\alpha})\).
- *Certificate of \(\int t^{\gamma-1-\alpha}(1-\eta\circ\psi)\) at \(0.8\).*
  - Cells: \(391460\) cells of ratio \(1.0001\) on \([10^{-13},10^4]\), with exact weight integrals. The lower bound uses \(1-\eta(\psi_{12}(t_0))\) and the upper uses \(1-\eta(\psi^+(t_1))\).
  - The \(\psi\) ceiling: \(\psi^+=\psi_{12}e^{2^{12}5^{-12}t/10}\). This holds because \(\psi(y)\le ye^{y/10}\) for \(y\le10\ln2\) (orbit bootstrap \(p^j(5^{-k}y)\le2\cdot5^{j-k}y\)), \(\psi(t)=p^K(\psi(5^{-K}t))\), and \(p^K\) has nonnegative coefficients and degree \(2^K\).
  - Head \((0,10^{-13}]\): \(1-\eta(z)=(z^{1-\gamma}/g-z)/2\) with BJ.3's \(g\)-bracket and \(z=\psi(t)\in[t,te^{t/10}]\).
  - Tail \([10^4,\infty)\): \(1-2/(2+\psi_{12}(10^4))\le1-\eta\circ\psi\le1\).
- *Derived limit.* \(1/((5/3)^N\mathcal R(5^Nt))=\lim_m1/((3/5)^mr_{m+N}(5^{-m}t))\). By BJ.1 this lies within \(3(3/5)^N\) of \(\lim_m1/((3/5)^mr_\infty(5^{-m}t))=t^{1-\gamma}/\Theta_{AT}(\psi(t))\) (BH.6's pointwise limit from AT.1). Apply dominated convergence with the dominator of the upper bound. \(\square\)

**Evidence (floats; consistency only).**
- By BI's recursion quadrature at \(n=100\), \(\rho^{10}\delta_{10}^{(100)}=1.135914525\) and \(\rho^{20}\delta_{20}^{(100)}=1.135914308\), inside the certified band.
- At \(n=N+40\) the same quantity is \(1.135599908\) (\(N=10\), \(n=50\)) and \(1.135544509\) (\(N=20\), \(n=60\)), below the band. The gap times \(\rho^{n-N}\) is \(0.6025\) and \(0.7082\).
- So the finite-\(n\) error at the corner neighbours looks like \(C\rho^{-(n-N)}\) relative to \(\rho^{-N}\). This is why BI's \(n=40\) floats drifted down with \(N\) (\(1.12988\) at \(N=15\)). It is not proved.

**What stays open.**
- Conjecture BI.7's third item: a rate for \(\max_{V_N}|u_n-u_\infty|\), hence for the uniform convergence. The floats above suggest \(O(\rho^{-n})\) at the corner neighbours (Conjecture BJ.7).
- A closed form for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\).
- A continuum equation for \(\bar u\).
- Whether \(\rho^N\delta_N\) is eventually monotone.
- BH.7, BF.2, and the Stage 2 question are untouched.

**CONJECTURE BJ.7 (floats only; PROVED by Theorem BK below with exponent \(c=1\): \(\sup_{V_n}|u_n-\bar u|\le C_\alpha\|x\|_\infty\rho^{-n}\), \(C_{0.8}\le6.2192\), and at the corner neighbours with \(\tfrac32|x_a|+\tfrac125^{-N}|x_b-x_c|\) in place of \(C_\alpha|x_a|\) up to the factor \(\log(D'_\infty/D'_n)\); the text records it as first stated).** For \(\alpha>\gamma\) there is \(C_\alpha\) with \(|u_n(\nu^N)-\bar u(\nu^N)|\le C_\alpha|x_a|\rho^{-n}\) at the corner neighbours (floats: \(C_{0.8}\approx0.6\)–\(0.7\)), and \(\sup_{V_n}|u_n-\bar u|\to0\) at a rate \(O(\rho^{-cn})\) for some \(c>0\).

**NOT THIS.** Theorem BJ is a lower bound, matching BI's upper bound, for the modulus of continuity of one graph Dirichlet extension's limit at the corners of the combinatorial gasket, plus five digits of one corner constant. It is a statement about a continuous function on the gasket, not a field, a flux or a force. No flux functional on the infinite gasket is constructed, and Stage 1 is not touched. Not thrust.

**Witness, computed by `scripts/gasket_fractional_corner_sharpness.py`** (about 5 min; needs sympy, mpmath, python-flint).
- *Sympy:* the reciprocal step, the \(D_k\) recursion and base, \(\beta\le3/5\), \(5^{1-\gamma}=5/3\), the neighbour-mean recursion, the \(\psi\)-lemma constant, the closed-form \(\int t^{-\alpha}\min(K/3,Gt^{\gamma-1})\).
- *Exact rationals (levels \(0\)–\(4\), \(s\in\{1/1000,1/2,7/5,10\}\)):* \(r_n\) against dense solves and the reciprocal step. The neighbour mean \(=x_a\eta_n\) (levels \(2\)–\(4\)).
- *`arb`:* BJ.1 and its bound, \(\eta_n=\eta-\prod\beta(\tfrac12+\eta\circ p^n)\). Certificates: BI.4 re-run, the \(g\)-bracket, \(M_{0.8}\), \(c_-(0.75,0.8,0.9)\), the finite-\(N\) upper bounds.
- *Floats:* \(\rho^N\delta_N^{(n)}\) at \(N=10,20\), \(n=100\) and \(n=N+40\).

The script exits nonzero if any check fails.

## Theorem BK. The uniform convergence of Theorem BI has rate \(\rho^{-n}\) on all of \(V_n\) (Conjecture BJ.7 proved, with exponent \(c=1\)): an exact one-level transfer with weight \(\ge1\)

Theorem BI proved \(u_n\to\bar u\) uniformly for \(\alpha>\gamma\), but its proof gives no rate in \(n\). The floats of Theorem BJ suggested that the error at the corner neighbours is about \(0.6\,|x_a|\rho^{-n}\) (Conjecture BJ.7). This theorem proves a rate \(\rho^{-n}\) at every vertex of \(V_n\) at once. The tool is an exact change of variables that carries level \(n+1\) back to level \(n\). It turns \(u_{n+1}(x)\) into a reweighted average of the level-\(n\) integrand, and the weight is pointwise \(\ge1\). The total extra weight is the increment of the normaliser \(D'_n=\rho^{-n}D_n\). That increment is \(O(\rho^{-n})\), by a certified integral.

**Setting.** As in Theorems BG–BJ: Assumption 4, corner data \(d\), \(x=d-\bar d\), \(u_n\) the \(L_n^{\alpha}\)-Dirichlet extension, \(\bar u\) its limit, \(r_n(s)\) the corner resolvent, \(p(y)=y(y+5)\), \(\beta(y)=(y+6)/((y+2)(y+5))\) (\(=\) AQ's \(c\)), \(e(y)=2y/(2+y)\), \(\rho=5^{\alpha}/3\), \(\gamma=\log3/\log5\), \(\mathcal R_n(\tau)=(3/5)^nr_n(5^{-n}\tau)\). Write
\[
D'_n=\int_0^\infty\tau^{-\alpha}\mathcal R_n(\tau)\,d\tau=\rho^{-n}D_n,\qquad \mu_n(d\tau)=\tau^{-\alpha}\mathcal R_n(\tau)\,d\tau/D'_n .
\]
For \(x\in V_m\) and \(n\ge m\), BG's exact profile reads \(u_n(x)-\bar d=\langle F_n\rangle_{\mu_n}\), with \(F_n(\tau)=h_m(p^{n-m}(5^{-n}\tau);x)\) (substitute \(s=5^{-n}\tau\)). **Assumptions used as stated:** AQ (1)–(3) (the step \(r_{n+1}(y)=r_n(p(y))/(\beta(y)-e(y)r_n(p(y)))\), \(r_0=1/(3+\cdot)\), \(r_n\uparrow r_\infty\)); BG.3 (the exact profile, \(|h_m(z;x)|\le\|x\|_\infty\) off the corners, and \(\bar u=\lim u_n\) on \(V_*\)); BH.2 (\(0<\eta\le1\), so \(1/r_\infty=2+s-2\eta\ge s\)); BJ.4 (\(h_N=x_a\eta_N+\frac{x_b-x_c}2\prod_{i<N}(5+p^iz)^{-1}\) at \(\nu_b^N\), \(\eta_N\le\eta\)); BI.4 (\(\sup_{s>0}s^{1-\gamma}r_\infty(s)\le G=0.6031\), re-run live); AZ (1), (3) (\(D'_\infty=\pi A_\alpha/(2\sin\pi\alpha)\); \(A_{0.8}\) for the asymptotic constant only; \(C_{0.8}\) in the derived part only). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BK.**

1. **BK.1 (exact one-level transfer; every \(\alpha\in(0,1)\), \(n\ge m\), \(x\in V_m\)).** Put \(\tau'=5^np(y)\) with \(y\in(0,\infty)\), and
   \[
   \omega_n(\tau')=\frac{(1+y/5)^{\alpha}}{1+2y/5}\cdot\frac{3}{5}\cdot\frac{r_{n+1}(y)}{r_n(p(y))}=\frac{(1+y/5)^{\alpha}}{(1+2y/5)\,\tfrac53\bigl(\beta(y)-e(y)\,r_n(p(y))\bigr)} .
   \]
   Then
   \[
   u_{n+1}(x)-\bar d=\frac{\langle\omega_nF_n\rangle_{\mu_n}}{\langle\omega_n\rangle_{\mu_n}},\qquad \langle\omega_n\rangle_{\mu_n}=\frac{D'_{n+1}}{D'_n}.
   \]
2. **BK.2 (the weight is \(\ge1\); per-step bound).** \((1+2y/5)\tfrac53\beta(y)=1-\dfrac{y(y+4)}{3(y+2)(y+5)}\le1\). So \(\omega_n\ge(1+y/5)^{\alpha}\ge1\) everywhere, for every \(\alpha\ge0\). Hence:
   - *Derived:* \(D'_n=\rho^{-n}D_n\) is strictly increasing in \(n\), for every \(\alpha\).
   - for every constant \(c\),
   \[
   |u_{n+1}(x)-u_n(x)|\le2\sup|F_n-c|\Bigl(1-\frac{D'_n}{D'_{n+1}}\Bigr)\le2\sup|F_n-c|\,\log\frac{D'_{n+1}}{D'_n}.
   \]
3. **BK.3 (rate in terms of the normaliser; \(\alpha>\gamma\)).** Here \(D'_n\uparrow D'_\infty=\pi A_\alpha/(2\sin\pi\alpha)<\infty\). For every \(n\) and every \(x\in V_n\),
   \[
   |u_n(x)-\bar u(x)|\le2\|x\|_\infty\log\frac{D'_\infty}{D'_n}.
   \]
   At the corner neighbours \(\nu_b^N\), \(N\le n\), the data split as follows:
   \[
   |u_n(\nu_b^N)-\bar u(\nu_b^N)|\le\Bigl(\tfrac32|x_a|+\tfrac12 5^{-N}|x_b-x_c|\Bigr)\log\frac{D'_\infty}{D'_n}.
   \]
   For the data-free corner deficit of BJ, \(|\delta_N^{(n)}-\delta_N|\le\tfrac32\log(D'_\infty/D'_n)\).
4. **BK.4 (certified rate: Conjecture BJ.7 PROVED, with \(c=1\)).** For every \(k\),
   \[
   D'_{k+1}-D'_k=\rho^{-k}\int_0^\infty(\omega_k(y)-1)\,p^{-\alpha}r_k(p)\,p'\,dy\le J_\alpha\rho^{-k},
   \]
   with \(J_\alpha=\int_0^\infty(\omega_{\rm hi}-1)\,p^{-\alpha}\bar r(p)\,p'\,dy\). Here \(\bar r(s)=\min(1/s,Gs^{\gamma-1})\) and
   \[
   \omega_{\rm hi}=\min\Bigl(\frac{(1+y/5)^{\alpha}}{(1+2y/5)\tfrac53(\beta-e\bar r(p))},\ \frac35\frac{(1+y/5)^{\alpha}\bar r(y)(3+p)}{1+2y/5}\Bigr),
   \]
   where the first entry is used only where \(\beta-e\bar r(p)>0\). So \(D'_\infty-D'_n\le\varepsilon_n:=J_\alpha\rho^{-n}/(1-\rho^{-1})\). With \(D'_n\ge D'_0=3^{-\alpha}\pi/\sin\pi\alpha\),
   \[
   \sup_{V_n}|u_n-\bar u|\le2\|x\|_\infty\log\frac{D'_\infty}{\max(D'_0,D'_\infty-\varepsilon_n)}\le C_\alpha\|x\|_\infty\rho^{-n},\qquad C_\alpha=\frac{2J_\alpha}{(1-\rho^{-1})D'_0}.
   \]
   Certified (python-flint `arb`): \(J_{0.75}\le1.256914\), \(J_{0.8}\le1.188154\), \(J_{0.9}\le1.101281\). Hence \(C_{0.75}\le12.5476\), \(C_{0.8}\le6.2192\), \(C_{0.9}\le1.9725\), and \(D'_\infty-D'_n\le6.901366\,\rho^{-n}\) at \(0.8\). With AZ's \(A_{0.8}\): \(D'_\infty(0.8)=6.40416129160\ldots\), so \(\limsup_n\rho^n\sup_{V_n}|u_n-\bar u|/\|x\|_\infty\le2J/((1-\rho^{-1})D'_\infty)\le2.155275\). Explicitly, at \(0.8\), \(\sup_{V_n}|u_n-\bar u|/\|x\|_\infty\le0.3556,\,0.04987,\,1.126\cdot10^{-3},\,2.572\cdot10^{-5}\) at \(n=10,20,40,60\).
5. **BK.5 (derived from AZ (3) as stated).** AZ's intercept \(C_{0.8}=-1.74325717839\ldots\) gives \(\rho^n(D'_\infty-D'_n)\to-\pi C_{0.8}/(2\sin0.8\pi)=4.658677573\ldots\). This sits under the certified \(6.901366\). So \(\rho^n\log(D'_\infty/D'_n)\to-C/A=-\nu_{0.8}\), and
   \[
   \limsup_n\rho^n\sup_{V_n}|u_n-\bar u|\le2|\nu_{0.8}|\,\|x\|_\infty=1.454890769\ldots\|x\|_\infty,
   \]
   \[
   \limsup_n\rho^n|\delta_N^{(n)}-\delta_N|\le\tfrac32|\nu_{0.8}|=1.0912\ldots\quad\text{(every }N).
   \]

**Proof of BK.1.**
- *Change of variables.* Substitute \(\tau=5^{n+1}y\) at level \(n+1\) and \(\tau'=5^np(y)\) at level \(n\). Then \(\tau'/\tau=1+y/5\) and \(d\tau'/d\tau=1+2y/5\) (sympy), and \(p^{n+1-m}(5^{-(n+1)}\tau)=p^{n-m}(5^{-n}\tau')\). So \(F_{n+1}(\tau)=F_n(\tau')\).
- *The ratio.* \(\mathcal R_{n+1}(\tau)/\mathcal R_n(\tau')=\tfrac35r_{n+1}(y)/r_n(p(y))\). Hence \(\tau^{-\alpha}\mathcal R_{n+1}(\tau)\,d\tau=\omega_n(\tau')\,\tau'^{-\alpha}\mathcal R_n(\tau')\,d\tau'\).
- *Second form of \(\omega_n\).* This is AQ's step (sympy). Checked in exact rationals for \(n=0,\ldots,3\) (so \(r_{n+1}\) up to level \(4\)) on five rational \(y\).
- *The identity.* Integrate against \(F_{n+1}\) and against \(1\), and divide. \(\square\)

**Proof of BK.2.**
- The identity is sympy. Since \(\beta-er_n(p)\le\beta\), it gives \(\omega_n\ge(1+y/5)^{\alpha}\ge1\). In exact rationals, \(\tfrac35r_{n+1}(y)/((1+2y/5)r_n(p))\ge1\) at levels \(0\)–\(3\).
- *Covariance bound.* Write \(\langle\cdot\rangle\) for \(\mu_n\). Then \(\langle\omega F\rangle/\langle\omega\rangle-\langle F\rangle=\langle(\omega-\langle\omega\rangle)(F-c)\rangle/\langle\omega\rangle\), with \(|\cdot|\le\langle|\omega-\langle\omega\rangle|\rangle\sup|F-c|/\langle\omega\rangle\). Since \(\omega\ge1\), \(\langle|\omega-\langle\omega\rangle|\rangle\le\langle\omega-1\rangle+(\langle\omega\rangle-1)=2(\langle\omega\rangle-1)\). Finally, \(1-1/t\le\log t\) for \(t\ge1\). \(\square\)

**Proof of BK.3.**
- *Telescoping.* \(F_k\) takes values in \([-\|x\|_\infty,\|x\|_\infty\)] for every \(k\) (BG.3; at the corners \(u_n=d\) exactly). Take \(c=0\) and sum BK.2 over \(k\ge n\). The logarithms telescope, and \(u_k(x)\to\bar u(x)\) (BG.3). As a by-product this re-proves the convergence on \(V_*\), uniformly on \(V_n\).
- *Corner neighbours.* By BJ.4, \(F_k=x_a\eta_N(z)+\frac{x_b-x_c}2\Pi_N(z)\), with \(\eta_N\in[-\tfrac12,1]\) and \(\Pi_N=\prod_{i<N}(5+p^iz)^{-1}\in(0,5^{-N}]\). Here \(\eta_N\ge\alpha_1-\beta/2=1/(2(z+5))>0\) for \(N\ge1\) (sympy), \(\eta_0=-\tfrac12\), and \(\eta_N\le\eta\le1\). Take \(c=x_a/4+(x_b-x_c)5^{-N}/4\), so that \(\sup|F_k-c|\le\tfrac34|x_a|+\tfrac14 5^{-N}|x_b-x_c|\).
- *Deficit.* For the deficit, take symmetric data. \(\square\)

**Proof of BK.4.**
- *The increment.* \(D'_{k+1}-D'_k=D'_k(\langle\omega_k\rangle-1)\). Change variables to \(y\): \(\tau'^{-\alpha}\mathcal R_k(\tau')\,d\tau'=\rho^{-k}p^{-\alpha}r_k(p)p'\,dy\), since \(5^{-k\alpha}(3/5)^k5^k=\rho^{-k}\) (sympy).
- *Bounds on \(r_k\).* \(r_k\le r_\infty\le\bar r\): use \(1/r_\infty\ge s\) and BI.4. Also \(r_k\ge r_0=1/(3+\cdot)\), by AQ (3).
- *Upper bounds on \(\omega_k\).* \(\omega_k\) increases with \(r_k(p)\) (second form). That gives the first entry of \(\omega_{\rm hi}\). The first form, with \(r_{k+1}(y)\le\bar r(y)\) and \(r_k(p)\ge1/(3+p)\), gives the second.
- *Normaliser.* \(D'_0=\int\tau^{-\alpha}/(3+\tau)=3^{-\alpha}\pi/\sin\pi\alpha\). Then \(\log(D'_\infty/D'_n)\le\log(D'_\infty/\max(D'_0,D'_\infty-\varepsilon_n))\le(D'_\infty-D'_n)/D'_0\).
- *Certificate of \(J_\alpha\).*
  - Cells: \(46075\) cells of ratio \(1.001\) on \([10^{-10},10^{10}]\). On each cell, the whole integrand is evaluated in `arb` on the cell ball (96 bits, outward), and the upper endpoint times the width is summed.
  - Head \((0,10^{-10}]\): \(\omega_{\rm hi}-1\le(\alpha y/5+\delta(y))/(1-\delta(y_0))\), with \(\delta=(4+y_0)y/30+(1+2y_0/5)Gy^{\gamma}\). This uses \(e\le y\), \(p\ge5y\), \(5^{\gamma-1}=3/5\) and \((1+y/5)^\alpha\le1+\alpha y/5\). The weight is \(\le G(5y)^{\gamma-1-\alpha}(5+2y_0)\). The monomials are integrable because \(2\gamma>\alpha\).
  - Tail \([10^{10},\infty)\): \((3+p)/(1+2y/5)\le5(y/2+5/4)\) (sympy) and \(p\ge y^2\). So the integrand is \(\le3(\tfrac15+\tfrac1{Y_0})^{\alpha}(\tfrac12+\tfrac5{4Y_0})(2+\tfrac5{Y_0})y^{-\alpha-1}\).
  - The head and tail together are below \(1.3\cdot10^{-5}\). \(\square\)

**Proof of BK.5.** AZ: \(R_\alpha(n)=\frac{2\sin\pi\alpha}{\pi}D_n=A\rho^n+C+o(1)\). So \(D'_n=\frac{\pi}{2\sin\pi\alpha}(A+C\rho^{-n}+o(\rho^{-n}))\), and \(\rho^n\log(D'_\infty/D'_n)\to-C/A\). Insert this into BK.3. \(\square\)

**Why this needs \(\alpha>\gamma\), and what it says below \(\gamma\).** BK.1 and BK.2 hold for every \(\alpha\). For \(\alpha<\gamma\), \(D_n\to D_\infty<\infty\), so \(\log(D'_{n+1}/D'_n)\to-\log\rho>0\). The per-step bound is then \(O(1)\) and not summable. This is consistent with the corner layer of Kept Failure BG.1 and Theorem BH. For \(\alpha>\gamma\), \(J_\alpha<\infty\) needs \(2\gamma>\alpha\) at \(y\to0\), which always holds, and \(\alpha>0\) at \(y\to\infty\). The \(\rho^{-k}\) in front comes from the weight \(\tau'^{-\alpha}\mathcal R_k\), and it is summable only for \(\rho>1\).

**Evidence (floats; consistency only).**
- *Normaliser.* By quadrature, \(D'_n(0.8)\) increases, and \(\rho^n(D'_\infty-D'_n)=4.6586,4.6587,4.6588\) at \(n=10,20,30\) (AZ's \(4.658678\)).
- *Dense steps.* Dense \(L^{0.8}\) solves for the data \((1,-\tfrac12,0)\): \(\max_{V_n}|u_{n+1}-u_n|=0.0137,0.0283,0.0337,0.0327,0.0299\) for \(n=1,\ldots,5\). The per-step bound \(2\|x\|\log(D'_{n+1}/D'_n)\) is \(0.297,0.235,0.182,0.140,0.109\).
- *Deficit.* Corner deficit at \(N=10\) against \(n=100\): \(\delta^{(100)}-\delta^{(n)}=3.696\cdot10^{-2},1.400\cdot10^{-2},2.087\cdot10^{-3},3.147\cdot10^{-4},4.756\cdot10^{-5}\) at \(n=15,20,30,40,50\). That is \(0.629,0.613,0.604,0.6028,0.6025\) times \(\rho^{-n}\): positive, inside the certified bound, and about \(0.55\) of the derived asymptotic ceiling \(1.0912\).

**What stays open.**
- *Lower bound.* Whether the order \(\rho^{-n}\) is attained, i.e. a lower bound \(|u_n-\bar u|\ge c\rho^{-n}\) at some vertex (floats: \(0.6025\,\rho^{-n}\) at the corner neighbours). This would follow from \(\operatorname{Cov}_{\mu_k}(\omega_k,F_k)\ge c'\rho^{-k}\) with a sign. Monotonicity of \(\omega_k\) in \(y\), which would give the sign by Chebyshev's inequality, is not proved.
- *The constant.* \(\lim\rho^n(\delta_N-\delta_N^{(n)})\), and whether \(u_n(\nu^N)\) increases in \(n\).
- *Other items.* Closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); BH.7, BF.2, and the Stage 2 question.

**CONJECTURE BK.6 (floats only; sharpness PROVED at \(\alpha=0.8\) for \(N\ge9\) by Theorem BL below, \(\delta_N-\delta_N^{(n)}\ge c_N\rho^{-n}\), \(c_N\ge0.1240\) for \(N\ge40\); false as worded at \(N=0\), Kept Failure BL.0; the limit is Conjecture BL.6; the text records it as first stated).** For \(\alpha>\gamma\) and every \(N\), \(\rho^n(\delta_N-\delta_N^{(n)})\to C_\delta(\alpha)>0\), the same for every \(N\) (floats: \(C_\delta(0.8)\approx0.6025\) at \(N=10\); BJ's \(0.60\)–\(0.71\) at \(n=N+40\) were pre-asymptotic). In particular, the rate \(\rho^{-n}\) of BK.4 is sharp.

**NOT THIS.** Theorem BK is a convergence rate for one graph Dirichlet extension of a discrete fractional operator on the combinatorial gasket toward its continuous limit, plus certified constants. It is a statement about functions on graphs, not a field, a flux or a force. No flux functional on the infinite gasket is constructed, and Stage 1 is not touched. Not thrust.

**Witness, computed by `scripts/gasket_fractional_uniform_rate.py`** (about 4 min; needs sympy, mpmath, numpy, python-flint).
- *Sympy:* the identity \((1+2y/5)\tfrac53\beta=1-y(y+4)/(3(y+2)(y+5))\); the change of variables; the two forms of \(\omega_n\); the \(\rho^{-k}\) scaling; \(1-1/t\le\log t\); the head and tail lemmas; \(\eta_N\ge-\tfrac12\).
- *Exact rationals (levels \(0\)–\(3\) with \(r_{n+1}\), five rational \(y\)):* AQ's step; \(\omega_n\ge1\); the ingredients of \(\omega_{\rm hi}\).
- *`arb`:* BI.4 re-run; \(J_{0.75},J_{0.8},J_{0.9}\) and \(C_\alpha\); the explicit \(n\)-bounds with AZ's \(A_{0.8}\); consistency with AZ's \(C_{0.8}\).
- *Floats:* \(D'_n\) by quadrature; dense one-level differences at levels \(1\)–\(6\); corner-deficit \(n\)-errors at \(n=15\)–\(50\) against \(n=100\).

The script exits nonzero if any check fails.

## Theorem BL. The rate \(\rho^{-n}\) of Theorem BK is attained (sharpness half of Conjecture BK.6): at \(\alpha=0.8\), for every \(N\ge9\), \(\delta_N^{(n)}\) strictly increases in \(n\) and \(\delta_N-\delta_N^{(n)}\ge c_N\rho^{-n}\), \(c_N\ge0.1240\) for \(N\ge40\)

Theorem BK proved \(\sup_{V_n}|u_n-\bar u|\le C_\alpha\|x\|_\infty\rho^{-n}\) above \(\gamma\). It left open whether the order \(\rho^{-n}\) is attained (Conjecture BK.6). BK noted that a sign for \(\operatorname{Cov}(\omega_k,F_k)\) would settle it, for example through monotone \(\omega_k\). That monotonicity is still not proved. This theorem does not need it. Only the deficit integrand has to be monotone, and that is exact: \(1-\eta_N\) is increasing. The weight only has to satisfy \(\omega_k\ge1\), which BK.2 already gives. The covariance is then split at one point \(y_0\). On \(y\ge y_0\) the integrand is bounded below by a level-free floor of \(\omega_k-1\), on order \(\rho^{-k}\). On \(y<y_0\) the loss is bounded above by BK.4's ceiling restricted to \(y<y_0\). Comparing the two needs only a ceiling on the current deficit \(\delta_N^{(k)}\). That ceiling is small once \(N\) is moderately large.

**Setting.** As in Theorem BK: \(p(y)=y(y+5)\), \(\alpha_1(z)=(z+4)/((z+2)(z+5))\), \(\beta(z)=(z+6)/((z+2)(z+5))\), \(\rho=5^{\alpha}/3\), \(\gamma=\log3/\log5\), \(\mathcal R_k(\tau)=(3/5)^kr_k(5^{-k}\tau)\), \(D'_k=\int\tau^{-\alpha}\mathcal R_k\), \(\mu_k=\tau^{-\alpha}\mathcal R_k\,d\tau/D'_k\), \(\omega_k\) BK.1's weight. The data-free corner deficit is defined, for \(N\ge1\) and \(k\ge N\), by
\[
\delta_N^{(k)}=\langle G_k\rangle_{\mu_k},\qquad G_k(\tau)=1-\eta_N\bigl(p^{k-N}(5^{-k}\tau)\bigr),
\]
where \(\eta_N=\alpha_1+\beta\,\eta_{N-1}\circ p\) and \(\eta_0=-\tfrac12\). This is \((d_a-\tfrac12(u_k(\nu_b^N)+u_k(\nu_c^N)))/x_a\), by BG.3 and BJ.4. Its limit is \(\delta_N=\lim_k\delta_N^{(k)}\) (BG.3, BK.3). **Assumptions used as stated:** BK.1 (the transfer, applied to \(G_k\); \(G_{k+1}(\tau)=G_k(\tau')\) is the same composition identity as for \(F_k\)); BK.2 (\(\omega_k\ge1\) and its second form); BK.4 (\(\omega_k\le\omega_{\rm hi}\), \(r_k\le\bar r\), \(r_k\ge1/(3+\cdot)\), the head lemma; \(J_{0.8}\) is re-run live); BI.1 (orbit lemma \(p^j(s)\le e^{1/2}5^js\) for \(5^js\le1\), re-checked); BI.4 (\(\mathcal R_k\le\min(\tfrac13,G\tau^{\gamma-1})\), \(G=0.6031\), re-run live); BJ.4 (the deficit formula). AZ is **not** used. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BL.**

1. **BL.1 (the deficit integrand is monotone; exact, every \(N\ge1\)).** \(1-\alpha_1-\beta=z/(z+2)\), so
   \[
   1-\eta_N(z)=\frac{z}{z+2}+\beta(z)\bigl(1-\eta_{N-1}(p(z))\bigr)=\sum_{j<N}\Bigl[\prod_{i<j}\beta(p^iz)\Bigr]\frac{p^jz}{p^jz+2}+\frac32\prod_{i<N}\beta(p^iz),
   \]
   with every term \(\ge0\). The functions \(\alpha_1,\beta\) are positive and strictly decreasing: the numerators of their derivatives are \(-(z^2+8z+18)\) and \(-(z+4)(z+8)\). Also \(\eta_1=1/(2(z+5))\). So by induction every \(\eta_N\), \(N\ge1\), is positive and strictly decreasing on \([0,\infty)\). Moreover \(\eta_N(0)=1-\tfrac9{10}(3/5)^{N-1}\) and \(\eta_N(\infty)=0\). Hence \(G_k\) is strictly increasing in \(\tau\), with values in \((0,1)\).
2. **BL.2 (covariance floor; every \(\alpha\in(\gamma,1)\), \(N\ge1\), \(k\ge N\), \(y_0>0\)).** By BK.1,
   \[
   \delta_N^{(k+1)}-\delta_N^{(k)}=\frac{\langle(\omega_k-1)(G_k-\langle G_k\rangle)\rangle_{\mu_k}}{\langle\omega_k\rangle_{\mu_k}}.
   \]
   Suppose \(\Delta\ge\delta_N^{(k)}\) and \(g\le1-\eta_N(p(y_0))\), with \(\Delta<g\). Then
   \[
   \delta_N^{(k+1)}-\delta_N^{(k)}\ \ge\ \frac{\rho^{-k}}{D'_{k+1}}\Bigl[(g-\Delta)L(y_0)-\Delta J(y_0)\Bigr],
   \]
   where
   \[
   L(y_0)=\int_{y_0}^\infty(\omega_{\rm lo}-1)\frac{p^{-\alpha}p'}{3+p}\,dy,\qquad J(y_0)=\int_0^{y_0}(\omega_{\rm hi}-1)p^{-\alpha}\bar r(p)p'\,dy,
   \]
   and \(\omega_{\rm lo}=(1+y/5)^{\alpha}/\bigl(1-\tfrac{y(y+4)}{3(y+2)(y+5)}\bigr)\le\omega_k\).
3. **BL.3 (normaliser floor; AZ-free).** \(D'_{k+1}-D'_k\ge\rho^{-k}L(0)\). So
   \[
   D'_N\ge D'_0+L(0)\frac{1-\rho^{-N}}{1-\rho^{-1}}.
   \]
   Also \(D'_\infty\le D'_0+J_\alpha/(1-\rho^{-1})\) (BK.4).
4. **BL.4 (deficit ceilings, every \(k\ge N\)).**
   - *(i) Direct.* \(\delta_N^{(k)}\le I_N/D'_N\), where
     \[
     I_N=\int_0^\infty\tau^{-\alpha}\min(\tfrac13,G\tau^{\gamma-1})H_N(\tau)\,d\tau,\qquad H_N=1-\eta_N(e^{1/2}5^{-N}\tau)\ \text{on }\tau\le5^N,\quad H_N=1\ \text{beyond}.
     \]
   - *(ii) Uniform in \(N\).*
     \[
     \delta_N^{(k)}\le\frac{GQ\rho^{-N}+\tfrac32(3/5)^NI_0}{D'_N},\qquad Q=\int_0^\infty t^{\gamma-1-\alpha}\min(1,S(t))\,dt,
     \]
     where \(S(t)=\sum_{l\ge0}(3/5)^l\varphi(5^lt)\), \(\varphi(x)=e^{1/2}x/2\) for \(x\le1\) and \(\varphi(x)=1\) for \(x>1\), and \(I_0=\int\tau^{-\alpha}\min(\tfrac13,G\tau^{\gamma-1})\).
5. **BL.5 (certified at \(\alpha=0.8\), \(y_0=1/10\); python-flint `arb`).**
   - *Inputs.* \(L(0)\ge0.341059\), \(L(1/10)\ge0.334504\), \(J(1/10)\le0.178544\), \(J_{0.8}\le1.188154\), \(D'_\infty\le9.120759\), \(Q\le19.525284\) (\(GQ\le11.775699\)), \(I_0\le6.547382\), and \(1-\eta_N(p(1/10))\ge0.584124\) for every \(N\ge9\).
   - *The bracket.* \((g-\Delta_N)L-\Delta_NJ>0\) in every case:
     - for \(N=9,\dots,13\), by (4i): \(\Delta_9\le0.374431\), \(\Delta_{13}\le0.167815\);
     - for every \(N\ge14\), by (4ii): \(\Delta_N\le0.207825\), and the bound is monotone in \(N\).
   - *Conclusion.* For every \(N\ge9\) and every \(n\ge N\),
     \[
     \delta_N^{(n+1)}>\delta_N^{(n)},\qquad c_N\rho^{-n}\le\delta_N-\delta_N^{(n)}\le4.6644\,\rho^{-n},\qquad c_N=\frac{(g-\Delta_N)L-\Delta_NJ}{D'_\infty(1-\rho^{-1})}.
     \]
     The upper bound is BK.3–BK.4. Certified values: \(c_9\ge2.095\cdot10^{-3}\), \(c_{10}\ge2.470\cdot10^{-2}\), \(c_{11}\ge4.290\cdot10^{-2}\), \(c_{12}\ge5.763\cdot10^{-2}\), \(c_{13}\ge6.960\cdot10^{-2}\), \(c_N\ge5.653\cdot10^{-2}\) for \(N\ge14\), and \(c_N\ge0.1240\) for \(N\ge40\).
   - *Sharpness.* Take data with \(x_a\ne0\). The mean of the two neighbour errors is \(x_a(\delta_N-\delta_N^{(n)})\). So \(\max_{\nu\in\{\nu_b^N,\nu_c^N\}}|u_n(\nu)-\bar u(\nu)|\ge|x_a|c_N\rho^{-n}\). For symmetric data (\(x_b=x_c\), \(\|x\|_\infty=|x_a|\)) this gives \(\sup_{V_n}|u_n-\bar u|\ge c_9\|x\|_\infty\rho^{-n}\) for every \(n\ge9\). Hence **the order \(\rho^{-n}\) of BK.4 is exact**: the "in particular" of Conjecture BK.6 is PROVED at \(\alpha=0.8\). *Derived:* for \(x_a>0\), the neighbour mean \(\tfrac12(u_n(\nu_b^N)+u_n(\nu_c^N))\) strictly decreases in \(n\ge N\).

**KEPT FAILURE BL.0 (Conjecture BK.6 as worded, at \(N=0\)).** BK.6 says "every \(N\)". At \(N=0\) the neighbours are the corners \(b,c\) themselves, so \(u_n=d\) there exactly. Then \(\delta_0^{(n)}=(x_a-\tfrac12(x_b+x_c))/x_a=\tfrac32\) for every \(n\) (equivalently \(G\equiv1-\eta_0=\tfrac32\) and the covariance vanishes). **Witness:** \(\rho^n(\delta_0-\delta_0^{(n)})=0\) for all \(n\), not \(C_\delta>0\). BK.6 is restated for \(N\ge1\). This is a wording fix, not a new phenomenon.

**Proof of BL.1.** The identity and the derivative numerators are checked in sympy. For \(N\ge2\), \(\eta_N=\alpha_1+\beta\cdot(\eta_{N-1}\circ p)\). Here \(\eta_{N-1}\circ p\) is positive and strictly decreasing, by induction and because \(p\) is increasing. A product of two positive strictly decreasing functions is strictly decreasing, and so is the sum. The value at \(0\) follows from \(\alpha_1(0)=\tfrac25\), \(\beta(0)=\tfrac35\), \(\eta_1(0)=\tfrac1{10}\). The series is the recursion unrolled. Checked in exact rationals for \(N=1\)–\(6\) at five rational \(z\), against BJ's recursion. \(\square\)

**Proof of BL.2.**
- *Covariance form.* \(\langle\omega G\rangle/\langle\omega\rangle-\langle G\rangle=\langle(\omega-1)(G-\langle G\rangle)\rangle/\langle\omega\rangle\) (algebra).
- *Region \(y\ge y_0\).* Pass to \(y\) with \(\tau'=5^kp(y)\). On \(y\ge y_0\), \(G_k=1-\eta_N(p^{k-N+1}(y))\ge1-\eta_N(p(y_0))\ge g>\Delta\ge\langle G_k\rangle\) (BL.1, and \(p^j(x)\ge x\)). So the integrand is \(\ge(g-\Delta)(\omega_{\rm lo}-1)\ge0\).
- *Region \(y<y_0\).* Since \(G_k\ge0\), the integrand is \(\ge-\Delta(\omega_{\rm hi}-1)\).
- *Weights.* \(\tau'^{-\alpha}\mathcal R_k\,d\tau'=\rho^{-k}p^{-\alpha}r_k(p)p'\,dy\) (BK.4), with \(1/(3+p)\le r_k(p)\le\bar r(p)\). Finally, \(D'_k\langle\omega_k\rangle=D'_{k+1}\). \(\square\)

**Proof of BL.3.** \(D'_{k+1}-D'_k=\rho^{-k}\int(\omega_k-1)p^{-\alpha}r_k(p)p'\,dy\) (BK.4). Use \(\omega_k\ge\omega_{\rm lo}\) and \(r_k\ge r_0\). \(\square\)

**Proof of BL.4.**
- *Direct.* \(\delta_N^{(k)}=\int\tau^{-\alpha}\mathcal R_kG_k/D'_k\), with \(D'_k\ge D'_N\) because \(D'\) is increasing (BK.2). Use \(\mathcal R_k\le\min(\tfrac13,G\tau^{\gamma-1})\) (BI.4). For \(\tau\le5^N\), \(p^{k-N}(5^{-k}\tau)\le e^{1/2}5^{-N}\tau\) (orbit lemma), and \(1-\eta_N\) is increasing. Beyond, \(G_k\le1\).
- *Uniform.* In BL.1's series with \(z=p^{k-N}(5^{-k}\tau)\), use \(\prod\beta\le(3/5)^j\) and \(w/(w+2)\le\min(1,w/2)\). Also \(p^jz\le e^{1/2}5^{j-N}\tau\) when \(5^{j-N}\tau\le1\). So \(G_k\le\min(1,S(5^{-N}\tau))+\tfrac32(3/5)^N\). Substituting \(\tau=5^Nt\) gives \(\int\tau^{\gamma-1-\alpha}\min(1,S(5^{-N}\tau))\,d\tau=5^{N(\gamma-\alpha)}Q=\rho^{-N}Q\).
- *Head of \(S\).* On \(5^{-(L+1)}<t\le5^{-L}\), \(S(t)=ct(3^{L+1}-1)/2+\tfrac52(3/5)^{L+1}\le(\tfrac32c+\tfrac52)t^{1-\gamma}\), with \(c=e^{1/2}/2\). This gives \(Q\)'s head. \(\square\)

**Proof of BL.5.**
- *\(L\) and \(J\).* Cells of ratio \(1.001\), with ball evaluation of the integrand on each cell (96 bits, outward). \(L\): \(41468\) and \(20734\) cells, with the part beyond \(10^8\) dropped (integrand \(\ge0\)). \(J(1/10)\): \(20734\) cells, plus BK's closed-form head on \((0,10^{-10}]\).
- *\(I_N\).* Cells of ratio \(1.002\) on \([10^{-6},5^N]\). The weight is decreasing (left endpoint) and \(H_N\) increasing (right endpoint, positive series). Head \(H_N(10^{-6})\cdot10^{-6(1-\alpha)}/(3(1-\alpha))\); tail \(G5^{N(\gamma-\alpha)}/(\alpha-\gamma)\).
- *\(Q\).* Cells of ratio \(1.002\) on \([10^{-12},1]\), using \(S\)'s closed form at the right endpoints. The head comes from the lemma above, and \([1,\infty)\) is exact.
- *\(g\).* For \(N\ge14\), \(g\) is the series truncated at \(j<14\), a lower bound for every \(N\ge14\).
- *Monotonicity in \(N\).* (4ii)'s numerator decreases in \(N\) and \(D'_N\)'s floor increases, so checking \(N=14\) covers all \(N\ge14\).
- *Summation.* Sum BL.2 over \(k\ge n\), using \(D'_{k+1}\le D'_\infty\) and \(\delta_N^{(k)}\to\delta_N\). \(\square\)

**Why this is a lower bound of the right order.** The positive part comes from \(y\gtrsim y_0\), that is \(\tau'\sim5^k\), where \(\mu_k\) has mass \(\asymp\rho^{-k}\) and \(\omega_k-1\asymp1\). The bulk \(\tau'=O(1)\) contributes only \(O(3^{-k})\), because there \(\omega_k-1\lesssim y^{\gamma}\). The deficit must also be small (\(\Delta_N<g\)). That is why \(N\) has to be moderately large in this argument.

**Evidence (floats; consistency only).** At \(N=10\), \(\rho^n(\delta^{(100)}-\delta^{(n)})=0.6300,0.6414,0.6427,0.6339,0.6244,0.6127,0.6065,0.6041,0.6028\) at \(n=10,11,12,14,16,20,25,30,40\). The sequence \(\delta_{10}^{(n)}\) increases, and every value lies in \([c_{10},4.6644]\).

**What stays open.**
- *The limit.* Existence of \(C_\delta=\lim\rho^n(\delta_N-\delta_N^{(n)})\) and its independence of \(N\). The floats put it near \(0.6025\), about \(4.9\) times the certified floor \(0.1240\).
- *Small \(N\).* \(1\le N\le8\) at \(0.8\): the bracket fails there because the deficit ceiling exceeds \(g\). The floats suggest the sign still holds. Settled by Theorem BN below (monotone \(\omega_k\), no deficit ceiling needed): \(c_N>0\) for every \(N\ge1\).
- *Other exponents.* \(\alpha\ne0.8\) (the argument is general; only the certificate is at \(0.8\)). Settled for large \(N\) at every \(\alpha\in(\gamma,1)\) by Theorem BM below, with certified thresholds at \(0.7\), \(0.75\), \(0.85\), \(0.9\), \(0.95\).
- *Monotone \(\omega_k\).* Floats show \(\omega_k\) increasing in \(y\) for \(k=0,1,3,6,10,20\). Not proved, and not needed here. Proved by Theorem BN.2 for every \(k\) and every \(\alpha\ge1/17\).
- *Other items.* Closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); BH.7, BF.2, and the Stage 2 question.

**CONJECTURE BL.6 (floats only; replaces BK.6 for \(N\ge1\)).** For \(\alpha>\gamma\) and every \(N\ge1\), \(\rho^n(\delta_N-\delta_N^{(n)})\to C_\delta(\alpha)>0\), independent of \(N\) (floats: \(C_\delta(0.8)\approx0.6025\)), and \(\delta_N^{(n)}\) increases in \(n\) for every \(N\ge1\).

**NOT THIS.** Theorem BL is a lower bound on the convergence error of one graph Dirichlet extension of a discrete fractional operator on the combinatorial gasket. It shows that BK's rate cannot be improved. It is a statement about functions on graphs, not a field, a flux or a force. No flux functional on the infinite gasket is constructed, and Stage 1 is not touched. Not thrust.

**Witness, computed by `scripts/gasket_fractional_rate_sharpness.py`** (about 3 min; needs sympy, mpmath, python-flint).
- *Sympy:* \(1-\alpha_1-\beta=z/(z+2)\); the derivative numerators; \(\eta_1\); \(\eta_N(0)\) for \(N\le11\); \(\omega_{\rm lo}\); the covariance form; the head lemma of \(S\).
- *Exact rationals:* the \(1-\eta_N\) series against BJ's recursion; monotonicity at sample points (\(N=1\)–\(6\)); the orbit lemma.
- *`arb`:* BI.4 and BK.4 re-run; \(L(0)\), \(L(1/10)\), \(J(1/10)\), \(I_9,\dots,I_{13}\), \(Q\), \(I_0\); the brackets and \(c_N\).
- *Floats:* \(\delta_{10}^{(n)}\) for \(n=10\)–\(40\) against \(n=100\).

The script exits nonzero if any check fails.

## Theorem BM. The sharpness of Theorem BL holds at every \(\alpha\in(\gamma,1)\) for all large \(N\), with certified thresholds \(N_0(0.7)=61\), \(N_0(0.75)=16\), \(N_0(0.85)=7\), \(N_0(0.9)=6\), \(N_0(0.95)=5\)

Theorem BL proved that BK's rate \(\rho^{-n}\) is attained, but its certificate was only at \(\alpha=0.8\) (for \(N\ge9\)). BL listed \(\alpha\ne0.8\) as open and noted that the argument itself is general. This theorem closes that item qualitatively and makes it explicit at five more exponents. No new identity is needed. BL.2–BL.4 hold for every \(\alpha\in(\gamma,1)\), and the \(N\)-uniform ceiling BL.4(ii) tends to \(0\) like \(\rho^{-N}\) at every such \(\alpha\). The rest is finiteness of the inputs and an `arb` re-run of BL's certificate at the new exponents.

**Setting.** As in Theorem BL: \(\delta_N^{(k)}=\langle G_k\rangle_{\mu_k}\) is the data-free corner deficit at the level-\(N\) corner neighbours, \(\delta_N=\lim_k\delta_N^{(k)}\), and \(L(y_0)\), \(J(y_0)\), \(I_N\), \(Q\), \(I_0\), \(D'_k\) are as defined in BL.2–BL.4. **Assumptions used as stated:** BL.1–BL.4 (stated there for every \(\alpha\in(\gamma,1)\)); BK.4 (\(J_\alpha<\infty\) for \(\alpha>\gamma\), since \(2\gamma>\alpha\) at \(y\to0\) and \(\alpha>0\) at \(y\to\infty\); \(D'_\infty\le D'_0+J_\alpha/(1-\rho^{-1})\); \(J_\alpha\) is re-run live at each exponent); BI.4 (\(G=0.6031\), re-run live).

**THEOREM BM.**

1. **BM.1 (qualitative; every \(\alpha\in(\gamma,1)\), every \(y_0>0\)).**
   - *Inputs are finite and positive.* \(L(y_0)>0\), because \(\omega_{\rm lo}=(1+y/5)^{\alpha}/B\) with \(0<B=1-y(y+4)/(3(y+2)(y+5))<1\) for \(y>0\). \(J(y_0)\le J_\alpha<\infty\) (BK.4). \(Q<\infty\): its head integrand is at most \(K_St^{-\alpha}\), integrable at \(0\) since \(\alpha<1\), and its tail \(t^{\gamma-1-\alpha}\) is integrable at \(\infty\) since \(\alpha>\gamma\). \(I_0<\infty\) for the same reason, and \(D'_N\ge D'_0>0\).
   - *The bracket's positive side does not degenerate.* The \(j=0\) term of BL.1's positive series gives \(g_N:=1-\eta_N(p(y_0))\ge p(y_0)/(p(y_0)+2)>0\) for every \(N\ge1\).
   - *The ceiling vanishes.* By BL.4(ii), \(\sup_{k\ge N}\delta_N^{(k)}\le\Delta_N:=(GQ\rho^{-N}+\tfrac32(3/5)^NI_0)/D'_0\to0\), since \(\rho=5^\alpha/3>1\).
   - *Conclusion.* There is an \(N_0(\alpha)\) with \((g_N-\Delta_N)L(y_0)-\Delta_NJ(y_0)>0\) for every \(N\ge N_0(\alpha)\). So by BL.2 and BL's summation, for every \(N\ge N_0(\alpha)\) and every \(n\ge N\),
     \[
     \delta_N^{(n+1)}>\delta_N^{(n)},\qquad c_N\rho^{-n}\le\delta_N-\delta_N^{(n)}\le\tfrac32\frac{J_\alpha}{(1-\rho^{-1})D'_0}\,\rho^{-n},
     \]
     with \(c_N>0\) and \(\liminf_Nc_N\ge\dfrac{p(y_0)}{p(y_0)+2}\cdot\dfrac{L(y_0)}{D'_\infty(1-\rho^{-1})}>0\). For symmetric data this gives \(\sup_{V_n}|u_n-\bar u|\ge c_{N_0}\|x\|_\infty\rho^{-n}\) for \(n\ge N_0\). **Hence the order \(\rho^{-n}\) of BK.4 is exact at every \(\alpha\in(\gamma,1)\).** This is the sharpness half of Conjecture BL.6 at every exponent above \(\gamma\), for large \(N\).
2. **BM.2 (certified thresholds; python-flint `arb`, outward rounding).** BL's certificate (BL.5) is re-run at five exponents. The direct ceiling BL.4(i) covers \(N_0\le N<N_u\). The \(N\)-uniform ceiling BL.4(ii) is checked at \(N_u\), and since it is monotone in \(N\) it covers every \(N\ge N_u\).

   | \(\alpha\) | \(y_0\) | \(L(y_0)\ge\) | \(J(y_0)\le\) | \(J_\alpha\le\) | \(D'_\infty\le\) | \(Q\le\) | direct \(N\) | uniform from | \(c_{N_0}\ge\) | \(c_N\ge\) (\(N\ge N_u\)) | \(\liminf c_N\ge\) |
   |---|---|---|---|---|---|---|---|---|---|---|---|
   | 0.7 | 1/10 | 0.454305 | 0.142156 | 1.345626 | 50.543638 | 64.214389 | 61 | 62 | \(7.006\cdot10^{-3}\) | \(2.044\cdot10^{-3}\) | 0.1902 |
   | 0.75 | 1/10 | 0.387944 | 0.158853 | 1.256914 | 14.176949 | 23.250101 | 16–17 | 18 | \(1.130\cdot10^{-2}\) | \(1.278\cdot10^{-2}\) | 0.1555 |
   | 0.85 | 1/50 | 0.297166 | 0.079980 | 1.136750 | 7.533153 | 21.442153 | 7–9 | 10 | \(2.304\cdot10^{-3}\) | \(1.330\cdot10^{-2}\) | 0.0636 |
   | 0.9 | 1/10 | 0.255428 | 0.230594 | 1.101281 | 7.512586 | 29.483613 | 6–7 | 8 | \(1.948\cdot10^{-2}\) | \(8.936\cdot10^{-3}\) | 0.0673 |
   | 0.95 | 1/10 | 0.225980 | 0.265862 | 1.081479 | 10.164569 | 60.060836 | 5–7 | 8 | \(1.060\cdot10^{-2}\) | \(1.246\cdot10^{-2}\) | 0.0371 |

   So, for every \(N\ge N_0(\alpha)\) and every \(n\ge N\), \(\delta_N^{(n)}\) strictly increases in \(n\) and \(c_N\rho^{-n}\le\delta_N-\delta_N^{(n)}\le C^{\rm up}_\alpha\rho^{-n}\), with ceilings \(C^{\rm up}_\alpha\le40.6260,\ 9.4106,\ 2.6545,\ 1.4794,\ 0.6559\) (BK.3–BK.4). Together with BL (\(\alpha=0.8\), \(N_0=9\)), the certified thresholds are \(N_0=61,16,9,7,6,5\) at \(\alpha=0.7,0.75,0.8,0.85,0.9,0.95\).

**Proof of BM.1.** The finiteness and positivity items are the lemmas listed (sympy for \(B\), the \(j=0\) term and the head exponent; BK.4 for \(J_\alpha\)). BL.4(ii) bounds \(\delta_N^{(k)}\) for every \(k\ge N\) with \(D'_N\ge D'_0\). Since \(g_N\ge g_*:=p(y_0)/(p(y_0)+2)\) and \(\Delta_N\to0\), the bracket is \(\ge(g_*-\Delta_N)L(y_0)-\Delta_NJ(y_0)>0\) once \(\Delta_N<g_*L/(L+J)\). Summing BL.2 over \(k\ge n\) with \(D'_{k+1}\le D'_\infty<\infty\) gives \(c_N\). The limit of \(c_N\) uses \(\Delta_N\to0\). \(\square\)

**Proof of BM.2.** The same cells as BL.5 (ratio \(1.001\) for \(L\), \(J(y_0)\), \(J_\alpha\); \(1.002\) for \(I_N\) and \(Q\); 96 bits, outward). The \(g\) used for \(N\ge N_u\) is BL.1's series truncated at \(j<\min(N_u,40)\), a floor for every \(N\ge N_u\); the \(\liminf\) column uses the truncation at \(j<40\). The values of \(y_0\) were chosen by a coarse float scan over \(\{1/50,1/20,1/10,1/5,2/5\}\); the certificate does not depend on how they were chosen. \(\square\)

**Why \(N_0\) grows as \(\alpha\downarrow\gamma\).** \(\Delta_N\asymp\rho^{-N}/D'_N\) and \(\rho\to1\) as \(\alpha\to\gamma\) (\(\rho(0.7)=1.02839\)). So the ceiling decays slowly, and \(N_0(\alpha)\to\infty\) in this argument. This is a limit of the method, not evidence against sharpness at small \(N\).

**Evidence (floats; consistency only).** At \(\alpha=0.9\), \(N=6\): \(\rho^n(\delta^{(60)}-\delta^{(n)})=0.3277,0.3488,0.3597,0.3659,0.3659,0.3648,0.3643\) at \(n=6,7,8,10,12,16,20\). The sequence \(\delta_6^{(n)}\) increases, and every value lies in \([c_6,1.4794]\). The limit is not BL's \(0.6025\), so if \(C_\delta\) exists it depends on \(\alpha\), as Conjecture BL.6 already allows.

**What stays open.**
- *Small \(N\).* \(N<N_0(\alpha)\): \(1\le N\le8\) at \(0.8\), \(1\le N\le60\) at \(0.7\), and so on. The bracket fails there because the deficit ceiling is too close to \(g\). The floats suggest the sign still holds. Settled by Theorem BN below for every \(N\ge1\).
- *Uniformity in \(\alpha\).* BM.1 is pointwise in \(\alpha\). No certificate covers an interval of \(\alpha\), and no explicit formula for \(N_0(\alpha)\) is given.
- *The limit.* Conjecture BL.6: existence of \(C_\delta(\alpha)=\lim\rho^n(\delta_N-\delta_N^{(n)})\) (floats \(0.6025\) at \(0.8\), \(\approx0.364\) at \(0.9\)).
- *Other items.* Monotone \(\omega_k\); closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); Conjecture BH.7, Conjecture BF.2, and the Stage 2 question.

**NOT THIS.** Theorem BM extends a lower bound on the convergence error of one graph Dirichlet extension of a discrete fractional operator on the combinatorial gasket to other exponents. It is a statement about functions on graphs, not a field, a flux or a force. No flux functional on the infinite gasket is constructed, and Stage 1 is not touched. Not thrust.

**Witness, computed by `scripts/gasket_fractional_rate_sharpness_alpha.py`** (about 1 min; needs sympy, mpmath, python-flint).
- *Sympy:* BL's identities re-run; \(2\gamma>1\); \(\omega_{\rm lo}>1\); the \(j=0\) floor of \(1-\eta_N\); \(Q\)'s head exponent.
- *Exact rationals:* BL.1's series and monotonicity re-run.
- *`arb`:* BI.4 re-run; at \(\alpha=0.7,0.75,0.85,0.9,0.95\): \(J_\alpha\), \(L(0)\), \(L(y_0)\), \(J(y_0)\), \(I_N\) for the direct range, \(Q\), \(I_0\); the brackets, \(c_N\) and \(\liminf c_N\).
- *Floats:* \(\delta_6^{(n)}\) at \(\alpha=0.9\) for \(n=6\)–\(20\) against \(n=60\).

The script exits nonzero if any check fails.

## Theorem BN. BK's one-level weight \(\omega_k\) is increasing, so \(\delta_N^{(n)}\) increases in \(n\) for every \(N\ge1\), and the rate \(\rho^{-n}\) is attained for every \(N\ge1\) at every \(\alpha\in(\gamma,1)\) (removes BL/BM's thresholds \(N_0(\alpha)\))

Theorems BL and BM proved the lower bound \(\delta_N-\delta_N^{(n)}\ge c_N\rho^{-n}\) only for \(N\ge N_0(\alpha)\) (\(N_0(0.8)=9\), \(N_0(0.7)=61\)). The reason was that their covariance floor tolerated a loss term, which had to be beaten by a ceiling on the current deficit. BK had already pointed out that a monotone weight \(\omega_k\) would give the sign of the covariance outright. BL recorded that monotonicity as floats only. This theorem proves it, for every level \(k\) and every \(\alpha\ge1/17\). The only new input is that the corner resolvent \(r_k\) is a Stieltjes function of unit mass. Chebyshev's covariance inequality then gives monotonicity of \(\delta_N^{(n)}\) for every \(N\ge1\). A pair form of the covariance gives a floor of order \(\rho^{-k}\) with no deficit ceiling at all.

**Setting.** As in Theorems BK and BL: \(r_k(s)=\tfrac12v^T(L_k+s)^{-1}v\) (AQ), \(p(y)=y(y+5)\), \(\beta(y)=(y+6)/((y+2)(y+5))\), \(e(y)=2y/(2+y)\), \(\rho=5^\alpha/3\), \(\gamma=\log3/\log5\). Write \(B(y)=1-\tfrac{y(y+4)}{3(y+2)(y+5)}\) and \(C(y)=\tfrac53(1+\tfrac{2y}5)e(y)\). BK.1's weight is \(\omega_k(y)=(1+y/5)^\alpha/d_k(y)\), where \(d_k=(1+\tfrac{2y}5)\tfrac53(\beta-e\,r_k(p))=B-C\,r_k(p)>0\). In BK.1's coordinate \(\tau'=5^kp(y)\), the deficit integrand is \(G_k=1-\eta_N(p^{k-N+1}(y))\), strictly increasing in \(y\) (BL.1), and \(\delta_N^{(k)}=\langle G_k\rangle_{\mu_k}\). **Assumptions used as stated:** AQ (1) (the decimation recursion equals \(\tfrac12v^T(L_k+s)^{-1}v\); re-checked in exact rationals at levels \(0\)–\(3\)); BK.1 (the transfer); BK.4 (the \(k\)-uniform two-form ceiling \(\omega_{\rm hi}\) and \(J_\alpha\), re-run live); BL.1 (monotone \(1-\eta_N\), positive series; re-run); BL.4(ii)'s pointwise bound \(G_k\le\min(1,S(5^{-N}\tau))+\tfrac32(3/5)^N\); BI.1 (orbit lemma); BI.3 (\(\mathcal R_k\ge1/(3+\tau)\)); BI.4 (\(G=0.6031\), re-run live). AZ is not used. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BN.**

1. **BN.1 (spectral lemma; every \(k\ge0\)).** \(L_k\) is a symmetric positive semidefinite graph Laplacian and \(\tfrac12|v|^2=1\). So \(r_k(s)=\sum_iw_i/(\lambda_i+s)\) with \(w_i\ge0\), \(\sum_iw_i=1\), \(\lambda_i\ge0\). Hence, for every \(s>0\),
   \[
   r_k(s)\le\frac1s,\qquad 0\le-r_k'(s)\le\frac{r_k(s)}{s}.
   \]
2. **BN.2 (the weight is increasing; every \(k\ge0\), every \(\alpha\in[1/17,1)\)).** With \(r=r_k(p(y))\), exactly
   \[
   d_k'-\frac{\alpha d_k}{5+y}=B'-\frac{\alpha B}{5+y}+r\Bigl(-C'+\frac{\alpha C}{5+y}\Bigr)+Cp'\,(-r_k'(p)).
   \]
   By BN.1 the last term is at most \(Cp'r/p\). The resulting coefficient of \(r\) is \(K=K_0+\alpha C/(5+y)\), where \(K_0=-C'+Cp'/p=\dfrac{2y(2y^2+10y+15)}{3(y+2)^2(y+5)}\ge0\). So, using \(r\le1/p\),
   \[
   d_k'-\frac{\alpha d_k}{5+y}\ \le\ \Psi_\alpha(y)=-\frac{P_\alpha(y)}{3(y+2)^2(y+5)^2},\qquad P_\alpha(y)=2\alpha y^3+(17\alpha-1)y^2+46\alpha y+40\alpha+10 .
   \]
   For \(\alpha\ge1/17\) every coefficient of \(P_\alpha\) is \(\ge0\) and \(P_\alpha(0)>0\), so \(\Psi_\alpha<0\) on \((0,\infty)\). Hence \((\log\omega_k)'=\alpha/(5+y)-d_k'/d_k>0\): **\(\omega_k\) is strictly increasing on \((0,\infty)\) for every level \(k\) and every \(\alpha\in[1/17,1)\).** This holds above and below \(\gamma\). (BL listed it as floats only.)
3. **BN.3 (monotone deficit; every \(\alpha\in[1/17,1)\), \(N\ge1\), \(k\ge N\)).** By BK.1, \(\delta_N^{(k+1)}-\delta_N^{(k)}=\operatorname{Cov}_{\mu_k}(\omega_k,G_k)/\langle\omega_k\rangle\). Both functions are strictly increasing in \(y\), and \(\mu_k\) has a positive density. So by Chebyshev's covariance inequality
   \[
   \delta_N^{(k+1)}>\delta_N^{(k)}.
   \]
   This is **the monotonicity half of Conjecture BL.6, for every \(N\ge1\)** (and also below \(\gamma\) down to \(\alpha=1/17\)). *Derived:* for \(x_a>0\), the corner-neighbour mean \(\tfrac12(u_n(\nu_b^N)+u_n(\nu_c^N))\) strictly decreases in \(n\ge N\).
4. **BN.4 (rate floor without a deficit ceiling; every \(\alpha\in(\gamma,1)\), every \(N\ge1\)).** Write \(\operatorname{Cov}=\iint_{y'<y}(\omega_k(y)-\omega_k(y'))(G_k(y)-G_k(y'))\,d\mu_k\,d\mu_k\). By BN.2 and BL.1 every term is \(\ge0\). Keep only the following pairs.
   - *Lower point.* \(y'\in A=\{\tau'\le\tau_0\}\) with \(\tau_0\le p(y_0)\), so \(y'\le y_0\). Then \(\omega_k(y')\le\omega_k(y_0)\le W_0:=\omega_{\rm hi}(y_0)\) (BN.2, BK.4). Also \(G_k(y')\le\bar G_N:=1-\eta_N(e^{1/2}5^{-N}\tau_0)\), by the orbit lemma (\(5^{-N}\tau_0\le1\)).
   - *Upper point.* \(y\ge y_0\). There \(\omega_k\ge\omega_{\rm lo}\), \(G_k\ge H_N(y):=1-\eta_N(p(y))\), and \(r_k(p)\ge1/(3+p)\).
   - *Mass of \(A\).* \(\mu_k(A)\ge M/D'_k\), with \(M=\int_0^{\tau_0}t^{-\alpha}/(3+t)\,dt\) (BI.3).

   Hence
   \[
   \delta_N^{(k+1)}-\delta_N^{(k)}\ \ge\ \frac{\rho^{-k}M\Lambda_N}{D'_kD'_{k+1}},\qquad \Lambda_N=\int_{y_0}^\infty(\omega_{\rm lo}-W_0)^+(H_N-\bar G_N)^+\frac{p^{-\alpha}p'}{3+p}\,dy .
   \]
   Summing over \(k\ge n\) (\(D'_k\le D'_\infty<\infty\), BK.4) gives, for every \(N\ge1\) and \(n\ge N\),
   \[
   c_N\rho^{-n}\ \le\ \delta_N-\delta_N^{(n)}\ \le\ \tfrac32\frac{J_\alpha}{(1-\rho^{-1})D'_0}\rho^{-n},\qquad c_N=\frac{M\Lambda_N}{D'^{\,2}_\infty(1-\rho^{-1})}.
   \]
   Here \(c_N>0\), because \(\omega_{\rm lo}\to\infty\) and \(H_N\to1>\bar G_N\) as \(y\to\infty\). So **BK.4's order \(\rho^{-n}\) is attained at the level-\(N\) corner neighbours for every \(N\ge1\) and every \(\alpha\in(\gamma,1)\)**, with no threshold \(N_0(\alpha)\). This is the sharpness half of Conjecture BL.6 for all \(N\ge1\). For \(N\ge N_1\), one certificate covers all \(N\): \(H_N\) is at least BL.1's series truncated at \(j<N_1\), and \(\bar G_N\le K_S(5^{-N_1}\tau_0)^{1-\gamma}+\tfrac32(3/5)^{N_1}\) (BL.4(ii), \(K_S=\tfrac34e^{1/2}+\tfrac52\)).
5. **BN.5 (certified; python-flint `arb`, outward rounding).** Parameters: \(y_0=1/10\), \(\tau_0=p(1/10)=51/100\), direct \(N=1\)–\(19\), uniform for \(N\ge N_1=20\) (\(\bar G_N\le1.652\cdot10^{-4}\)).

   | \(\alpha\) | \(D'_\infty\le\) | \(W_0\le\) | \(M\ge\) | \(c_1\ge\) | \(c_2\ge\) | \(c_4\ge\) | \(c_8\ge\) | \(\Lambda\ge\) (\(N\ge20\)) | \(c_N\ge\) (\(N\ge20\)) |
   |---|---|---|---|---|---|---|---|---|---|
   | 0.7 | 50.543638 | 1.173620 | 0.775966 | \(3.253\cdot10^{-4}\) | \(1.580\cdot10^{-3}\) | \(2.850\cdot10^{-3}\) | \(3.481\cdot10^{-3}\) | 0.324879 | \(3.575\cdot10^{-3}\) |
   | 0.75 | 14.176949 | 1.174783 | 0.963042 | \(1.136\cdot10^{-3}\) | \(5.548\cdot10^{-3}\) | \(1.003\cdot10^{-2}\) | \(1.226\cdot10^{-2}\) | 0.269989 | \(1.259\cdot10^{-2}\) |
   | 0.8 | 9.120759 | 1.175946 | 1.245022 | \(1.758\cdot10^{-3}\) | \(8.641\cdot10^{-3}\) | \(1.566\cdot10^{-2}\) | \(1.914\cdot10^{-2}\) | 0.226154 | \(1.966\cdot10^{-2}\) |
   | 0.85 | 7.533153 | 1.177111 | 1.716869 | \(2.166\cdot10^{-3}\) | \(1.071\cdot10^{-2}\) | \(1.946\cdot10^{-2}\) | \(2.380\cdot10^{-2}\) | 0.190824 | \(2.445\cdot10^{-2}\) |
   | 0.9 | 7.512586 | 1.178277 | 2.663483 | \(2.274\cdot10^{-3}\) | \(1.132\cdot10^{-2}\) | \(2.061\cdot10^{-2}\) | \(2.523\cdot10^{-2}\) | 0.162108 | \(2.591\cdot10^{-2}\) |
   | 0.95 | 10.164569 | 1.179445 | 5.509363 | \(1.836\cdot10^{-3}\) | \(9.195\cdot10^{-3}\) | \(1.680\cdot10^{-2}\) | \(2.057\cdot10^{-2}\) | 0.138586 | \(2.113\cdot10^{-2}\) |

   The full direct rows \(N=1\)–\(19\) are printed by the script. At \(\alpha=0.8\) this closes BL's gap \(N=1\)–\(8\): \(c_1\ge1.758\cdot10^{-3}\), …, \(c_8\ge1.914\cdot10^{-2}\). At \(\alpha=0.7\) it closes BM's gap \(N=1\)–\(60\). For large \(N\), BL/BM's constants are larger (e.g. \(c_N\ge0.1240\) for \(N\ge40\) at \(0.8\)), and both floors are valid. Combined with BK.3–BK.4, the error is two-sided \(\Theta(\rho^{-n})\) at every \(N\ge1\).

**Proof of BN.1.** This is the spectral theorem for \(L_k\). For the bounds: \(\lambda_i\ge0\) gives \(1/(\lambda_i+s)\le1/s\), and \(-r_k'=\sum w_i/(\lambda_i+s)^2\le\frac1s\sum w_i/(\lambda_i+s)\). Checks: exact rationals confirm \(\tfrac12v^T(L_n+s)^{-1}v=\) AQ's recursion at \(s=1/3,2,7/2\) for \(n=0\)–\(3\). sympy finds that \(r_0,\dots,r_3\) have only simple negative real poles, with positive residues summing to \(1\). \(\square\)

**Proof of BN.2.** The chain-rule identity, \(K_0\), \(B'+K_0/p=(y^2-10)/(3(y+2)^2(y+5)^2)\), \(B-C/p=(y+4)(2y+5)/(3(y+2)(y+5))\) and \(\Psi_\alpha=-P_\alpha/(3(y+2)^2(y+5)^2)\) are sympy identities. The inequality steps use \(C,p'\ge0\), BN.1 and \(K\ge0\). Note that the \(\alpha=0\) part \(B'+K_0/p\) is positive for \(y>\sqrt{10}\). The exponent term \(-\alpha(B-C/p)/(5+y)\) is what closes the bound, and this majorant fails for small \(\alpha\) (e.g. \(P_{1/100}<0\) at \(y=5,\dots,40\)). That is a limit of the method, not a statement about \(\omega_k\). \(\square\)

**Proof of BN.3.** Chebyshev's covariance inequality: \(\operatorname{Cov}=\tfrac12\iint(\omega(y)-\omega(y'))(G(y)-G(y'))\), and the integrand is \(>0\) off the diagonal. Integrability: \(G_k\in(0,1)\), and \(\langle\omega_k\rangle=D'_{k+1}/D'_k<\infty\). \(\square\)

**Proof of BN.4.** Every pair term is \(\ge0\), so restricting to \(A\times[y_0,\infty)\) is a lower bound. On that set, \(y'\le y_0\le y\), and the two factors are \(\ge(\omega_{\rm lo}(y)-W_0)^+\) and \(\ge(H_N(y)-\bar G_N)^+\). The measure of \(B\) in \(y\) is \(\rho^{-k}p^{-\alpha}r_k(p)p'\,dy/D'_k\ge\rho^{-k}p^{-\alpha}p'/(3+p)\,dy/D'_k\). Divide by \(\langle\omega_k\rangle=D'_{k+1}/D'_k\). Sum using \(\delta_N^{(k)}\to\delta_N\) (BK.3). \(\square\)

**Proof of BN.5.** \(W_0\) is BK.4's two-form ceiling at the point \(1/10\), with \(\bar r=\min(1/s,Gs^{\gamma-1})\). \(M\ge\tau_0^{1-\alpha}/((1-\alpha)(3+\tau_0))\). \(\Lambda_N\) uses cells of ratio \(1.01\) on \([1/10,10^6]\) (80 bits): ball evaluation of \(\omega_{\rm lo}\) and the weight, and \(H_N\) at the left endpoint (\(H_N\) increasing). The series is truncated once \(p^jz>10^{30}\), which drops positive terms. The part beyond \(10^6\) is dropped (integrand \(\ge0\)). \(\bar G_N\) is the full series at \(e^{1/2}5^{-N}\tau_0\), rounded up. \(J_\alpha\) and \(D'_\infty\le D'_0+J_\alpha/(1-\rho^{-1})\) are re-run as in BK.4. \(\square\)

**Evidence (floats; consistency only).**
- \(\omega_k\) increases on \(321\) log-grid points of \([10^{-5},10^3]\) for \(k=0,1,2,4,8,12\) at \(\alpha=0.45,0.8,0.95\) (mpmath, 40 digits).
- At \(0.8\), \(\delta_N^{(n)}\) increases over \(n=8,12,20\) for \(N=1\)–\(8\) (e.g. \(N=1\): \(0.9487\to0.9562\); \(N=8\): \(0.1183\to0.2379\)), and \(\delta^{(20)}-\delta^{(8)}\ge c_N(\rho^{-8}-\rho^{-20})\).
- The certified \(c_N\) are far below the floats' \(C_\delta(0.8)\approx0.6025\). BN gives existence of the order, not the constant.

**What stays open.**
- *The limit.* Conjecture BL.6 is now reduced to its limit part: does \(C_\delta(\alpha)=\lim\rho^n(\delta_N-\delta_N^{(n)})\) exist, and is it independent of \(N\)?
- *Small \(\alpha\).* Monotone \(\omega_k\) for \(\alpha<1/17\). Settled by Theorem BO below: monotone at every level iff \(\alpha\ge a^*\in[0.016778880172852,0.0167788801732]\); below \(a^*\) it fails at every level \(k\ge3\) (Kept Failure BO.1).
- *Other items.* A certificate uniform over an interval of \(\alpha\); closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BN is a monotonicity statement and a lower bound for the convergence of one graph Dirichlet extension of a discrete fractional operator on the combinatorial gasket. It is about functions on graphs, not a field, a flux or a force. No flux functional on the infinite gasket is constructed, and Stage 1 is not touched. Not thrust.

**Witness, computed by `scripts/gasket_fractional_monotone_weight.py`** (about 1 min; needs sympy, mpmath, python-flint).
- *Sympy:* the chain-rule identity; \(K_0\); \(B'+K_0/p\); \(B-C/p\); \(\Psi_\alpha\) and \(P_\alpha\); the coefficient signs at \(\alpha=1/17\); the small-\(\alpha\) remark; \(M\)'s floor; the two-atom Chebyshev identity; BL's checks re-run.
- *Exact rationals:* \(\tfrac12v^T(L_n+s)^{-1}v\) against AQ's recursion (\(n\le3\)); the pole and residue structure of \(r_0,\dots,r_3\); \(r_k\le1/s\) and \(-r_k'\le r_k/s\) at sample points.
- *`arb`:* BI.4 re-run; at each of the six exponents, \(J_\alpha\), \(W_0\), \(M\), \(\Lambda_N\) (\(N=1\)–\(19\) and uniform), and \(c_N\).
- *Floats:* \(\omega_k\) on a grid; \(\delta_N^{(n)}\) at \(0.8\).

The script exits nonzero if any check fails.

## Theorem BO. The exact exponent threshold for BN's monotone weight: \(\omega_k\) is increasing at every level iff \(\alpha\ge a^*\), with \(a^*\in[0.016778880172852,\,0.0167788801732]\) (Kept Failure BO.1: below \(a^*\) the weight is not monotone at any level \(k\ge3\))

Theorem BN.2 proved that BK's one-level weight \(\omega_k\) is strictly increasing in \(y\), at every level \(k\), for \(\alpha\ge1/17\). BN's own remark showed that its majorant fails at \(\alpha=1/100\), and left open whether this is a limit of the method or of the claim. This theorem settles that question. The claim itself fails for small \(\alpha\). The exact threshold is a number \(a^*\approx0.0167788801729\), certified to an interval of width \(<3.5\cdot10^{-13}\). The new input is a set of moment bounds on the corner resolvent. The moments \(m_j=\tfrac12v^TL_k^jv\) are local, so they are the same integers at every level \(k\ge3\). This makes one finite, exact-rational certificate cover every level.

**Setting.** As in Theorem BN: \(r_k(s)=\tfrac12v^T(L_k+s)^{-1}v\), \(v=e_a-e_b\), \(p(y)=y(y+5)\), \(e=2y/(2+y)\), \(B=1-\tfrac{y(y+4)}{3(y+2)(y+5)}\), \(C=\tfrac53(1+\tfrac{2y}5)e\), \(d_k=B-C\,r_k(p)>0\), \(\omega_k=(1+y/5)^\alpha/d_k\). Write \(\Phi_k(y)=(5+y)\,d_k'(y)/d_k(y)\) and \(a_k^*=\sup_{y>0}\Phi_k(y)\). **Assumptions used as stated:** AQ (1) (the recursion equals \(\tfrac12v^T(L_k+s)^{-1}v\); BN re-checked levels \(0\)–\(3\); here the \(s\to\infty\) expansion of the recursion is matched against the graph moments at \(k=0\)–\(3\)); BN.1 (\(r_k=\sum_iw_i/(\lambda_i+s)\), \(w_i\ge0\), \(\sum w_i=1\), \(\lambda_i\ge0\)); BN.2's chain-rule identity (re-checked); BK's \(d_k>0\). No \(\eta\) enclosure, no arb cell, and no AZ input is used. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BO.**

1. **BO.1 (reduction).** \((\log\omega_k)'=(\alpha-\Phi_k(y))/(5+y)\). So \(\omega_k\) is strictly increasing on \((0,\infty)\) iff \(\alpha\ge a_k^*\). If \(\alpha<a_k^*\), then \(\omega_k\) strictly decreases on an interval. The condition is linear in \(\alpha\), so one witness point covers a whole interval of exponents. The level-uniform threshold is \(a^*=\sup_ka_k^*\).
2. **BO.2 (moment bounds, every \(k\)).** With \(m_j=\sum_iw_i\lambda_i^j=\tfrac12v^TL_k^jv\), the exact geometric remainders
\[
\frac1{s+\lambda}=\sum_{j<J}\frac{(-\lambda)^j}{s^{j+1}}+\frac{(-\lambda)^J}{s^J(s+\lambda)},\qquad
\frac1{(s+\lambda)^2}=\sum_{j<J}\frac{(j+1)(-\lambda)^j}{s^{j+2}}+(-\lambda)^J\Big(\frac{J}{s^{J+1}(s+\lambda)}+\frac1{s^J(s+\lambda)^2}\Big)
\]
give two-sided bounds on \(r_k\) and \(-r_k'\) for every \(s>0\). Odd \(J\) gives upper bounds and even \(J\) gives lower bounds. The spectrum of \(L_k\) lies in \([0,8]\) (Gershgorin, maximum degree \(4\)), and \(m_1\ge2\) (\(m_1=3\) at \(k=0\), \(2\) for \(k\ge1\)). The chord of the convex functions \(\lambda\mapsto1/(s+\lambda)\) and \(\lambda\mapsto\lambda/(s+\lambda)^2\) on \([0,8]\) then gives the \(k\)-uniform bounds \(r_k\le\hat r:=1/s-2/(s(s+8))\) and \(-r_k'\le r_k/s-2/(s(s+8)^2)\). These sharpen BN.1's \(r_k\le1/s\) and \(-r_k'\le r_k/s\).
3. **BO.3 (locality).** \(m_j(k)=(L_k^j)_{aa}\) for \(j<2^k\), and \(m_j(k+1)=m_j(k)\) whenever \(j<2^k\). Hence \(m_0,\dots,m_7=1,2,6,22,92,418,1994,9780\) at every level \(k\ge3\). The lower levels have \(3^j\) (\(k=0\)), \(1,2,5,17,74,365,1889,9938\) (\(k=1\)) and \(1,2,6,22,91,408,1928,9419\) (\(k=2\)).
4. **BO.4 (upper certificate: monotone for \(\alpha\ge A_{\rm hi}\), every \(k\ge0\)).** Let \(A_{\rm hi}=0.0167788801732\) (exact rational).
   - *(i) Coarse, every \(k\).* Insert BO.2's chord bounds into BN.2's identity \(d'-\alpha d/(5+y)=B'-\alpha B/(5+y)+r(K_0+\alpha C/(5+y))+Cp'((-r')-r/p)\). This gives \(d_k'-\alpha d_k/(5+y)\le U-\alpha V\), with \(U=B'+\hat rK_0-2Cp'/(p(p+8)^2)\) and \(V=(B-C\hat r)/(5+y)>0\). The cleared numerator of \(U-A_{\rm hi}V\) has exactly two positive roots, \(12.5408\ldots\) and \(14.7913\ldots\). So \(\Phi_k<A_{\rm hi}\) off \((12,15)\). The coarse threshold \(\sup U/V=0.0168734\ldots\) is already below \(1/17\), and below BN's own majorant limit \(0.0215883\ldots\).
   - *(ii) Fine, every \(k\ge3\).* On \([12,15]\), use the \(J=7\) and \(J=8\) moment bounds with BO.3's integers, together with \(C'=4(y^2+4y+5)/(3(y+2)^2)>0\) and \(Cp'>0\). This gives \((5+y)d'_{\rm hi}-A_{\rm hi}d_{\rm lo}<0\). The check is an exact polynomial of degree \(17\) with no root on \([12,15]\) (Sturm).
   - *(iii) Exact, \(k=0,1,2\).* \(d_k\) is rational (AQ's recursion). \(a_{k,\rm hi}d_k-(5+y)d_k'\) has no root on \([0,\infty)\) and is positive there, with \(a_{k,\rm hi}<A_{\rm hi}\).

   Hence \(\omega_k\) is strictly increasing in \(y\), for every level \(k\ge0\), for every \(\alpha\ge A_{\rm hi}\). This extends BN.2 from \(\alpha\ge1/17\approx0.0588\) to \(\alpha\ge0.0167788801732\).
5. **BO.5 (per-level thresholds, exact rationals).** \(a_0^*\in[0.0155298932937,0.0155298932938]\), \(a_1^*\in[0.0167694660209,0.0167694660210]\), \(a_2^*\in[0.0167788799215,0.0167788799216]\), and \(a_k^*\in[A_{\rm lo},A_{\rm hi}]\) for every \(k\ge3\). So \(a^*=\sup_ka_k^*\in[0.016778880172852,0.0167788801732]\). The lower ends are sign changes at rational witness points, and the upper ends are BO.4.
6. **BO.6 (derived: Chebyshev beyond \(1/17\)).** BN.3's proof uses only that \(\omega_k\) and \(G_k\) are both increasing. So with BO.4, \(\delta_N^{(k+1)}>\delta_N^{(k)}\) holds for every \(N\ge1\), \(k\ge N\) and \(\alpha\in[A_{\rm hi},1)\).

**KEPT FAILURE BO.1 (BN.2 does not extend below \(a^*\)).** Claim: \(\omega_k\) is increasing in \(y\) at every level for every \(\alpha\in(0,1)\). This is false. **Witness:** \(y^*=13.778505644\). At this point the \(J=7,8\) moment bounds give, in exact rationals, \(\Phi_k(y^*)\ge0.0167788801728520\ldots\) for every \(k\ge3\). So \((\log\omega_k)'(y^*)<0\) for every \(\alpha<A_{\rm lo}=0.016778880172852\) and every \(k\ge3\). Then \(\omega_k\) is not monotone at any level \(k\ge3\). Levels \(0,1,2\) fail below \(a_0^*,a_1^*,a_2^*\) (BO.5). For \(\alpha\in(a_0^*,a^*)\), the base level is monotone but every level \(k\ge3\) is not. So the failure is not visible at \(k=0\). The failure is in the claim, not only in BN's majorant. The consequence is that BN.3's Chebyshev route to \(\delta_N^{(n+1)}>\delta_N^{(n)}\) is unavailable for \(\alpha<a^*\). This does not show that \(\delta_N^{(n)}\) fails to increase there. That question is open.

**Proof of BO.1.** \(\log\omega_k=\alpha\log(1+y/5)-\log d_k\), and \(d_k>0\). \(\square\)

**Proof of BO.2.** The two remainder identities are checked by sympy for \(J=1\)–\(8\). For \(\lambda\ge0\), their remainders have sign \((-1)^J\). Integrate against \(\sum_iw_i\delta_{\lambda_i}\), which has unit mass (BN.1). The chord identities are \(1/s-\lambda/(s(s+8))-1/(s+\lambda)=\lambda(8-\lambda)/(s(s+\lambda)(s+8))\ge0\) and \(\lambda/(s+\lambda)^2-\lambda/(s+8)^2=\lambda(8-\lambda)(2s+\lambda+8)/((s+\lambda)^2(s+8)^2)\ge0\) on \([0,8]\). Use \(1/(s+\lambda)^2=1/(s(s+\lambda))-\lambda/(s(s+\lambda)^2)\). Then \(m_1=\tfrac12v^TL_kv=\tfrac12(\deg a+\deg b)+[a\sim b]=3\) at \(k=0\), and \(2\) for \(k\ge1\). \(\square\)

**Proof of BO.3.** \((L^j)_{xy}\) is a sum over walks of length \(j\) from \(x\) to \(y\). The two corners are at graph distance \(2^k\) in \(V_k\). So for \(j<2^k\), \((L^j)_{ab}=0\), and by symmetry \(m_j=(L^j)_{aa}\). A closed walk of length \(j\) from \(a\) visits only vertices at distance \(\le j/2\). In \(V_{k+1}\), every vertex at distance \(<2^k\) from \(a\) lies in the corner subcell, which is a copy of \(V_k\). Only that subcell's two far corners (distance \(2^k\)) have a different degree. So the walk sums agree. The integers are checked against the graph matrices at \(k=1\)–\(6\), and against the \(s\to\infty\) expansion of AQ's recursion at \(k=0\)–\(3\). \(\square\)

**Proof of BO.4–BO.5.** All checks are exact polynomial computations in sympy, with rational coefficients: root isolation, Sturm counts on \([12,15]\) and \([0,\infty)\), and the signs of the denominators. In (i), \(K_0+\alpha C/(5+y)\ge0\) allows \(r\le\hat r\), and \(Cp'\ge0\) allows the \(-r'\) bound. \(B-C\hat r\ge B-C/p>0\) (BN). In (ii), \(d'=B'-C'r+Cp'(-r')\) is bounded above with \(r\ge\) (\(J=8\)) and \(-r'\le\) (\(J=7\)), and \(d\) is bounded below with \(r\le\) (\(J=7\)). The witness in BO.5 and in Kept Failure BO.1 bounds \(d'\) below and \(d\) above at \(y^*\), with the opposite parities. \(\square\)

**Evidence (floats; consistency only).** mpmath (40 digits) gives \(\sup\Phi_k=0.0155298932937551\) (\(k=0\)), \(0.0167694660209228\) (\(k=1\)), \(0.0167788799215301\) (\(k=2\)), and \(0.0167788801728761\) at \(y=13.7785056\) for \(k=3,5,10\). These lie inside the certified brackets. The level dependence beyond \(k=3\) is invisible at \(16\) digits, which is consistent with BO.3. An attempt to compute \(\delta_N^{(n)}\) by floats at \(\alpha=0.01\) (BI's quadrature) did not finish in \(3\) minutes. So there is no evidence either way on whether the deficit is monotone below \(a^*\).

**What stays open.**
- *Below \(a^*\).* Whether \(\delta_N^{(n)}\) still increases in \(n\) for \(\alpha<a^*\). The Chebyshev route is closed there. Settled by Theorem BP below: it does, for every \(\alpha\in(0,1)\), by a single-crossing covariance (the dip of BO.1 lies below the mean \(\langle\omega_k\rangle\ge\rho^{-1}\)).
- *A closed form for \(a^*\).* It is \(\sup_y(5+y)d_\infty'/d_\infty\) with \(r_\infty=1/(2+s-2\eta)\). No closed form is known.
- *The limit.* Conjecture BL.6's limit part: existence and \(N\)-independence of \(C_\delta(\alpha)\).
- *Other items.* A certificate uniform over an interval of \(\alpha\); closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BO is a statement about the exponent range in which one weight function of a graph Dirichlet extension is monotone. It is about functions on graphs, not a field, a flux or a force. No flux functional is constructed, and Stage 1 is not touched. Not thrust.

**Witness, computed by `scripts/gasket_fractional_weight_threshold.py`** (about 15 s; needs sympy, mpmath, numpy).
- *Sympy:* the remainder identities (\(J\le8\)); the chord lemmas; BN.2's identity and \(K_0\); \(C'>0\); the log-derivative reduction.
- *Exact integers:* the moments \(m_0,\dots,m_7\) from the graph matrices (\(k=1\)–\(6\)) and from the series of AQ's recursion (\(k=0\)–\(3\)).
- *Exact polynomials:* the coarse window (i); the fine Sturm check (ii); the level-\(0,1,2\) brackets (iii); the witness at \(y^*\).
- *Floats:* \(\sup\Phi_k\) at \(k=0,1,2,3,5,10\).

The script exits nonzero if any check fails.

## Theorem BP. The corner deficit \(\delta_N^{(n)}\) increases in \(n\) for every \(\alpha\in(0,1)\), also below \(a^*\): a single-crossing covariance replaces Chebyshev

Kept Failure BO.1 showed that BK's one-level weight \(\omega_k\) is not monotone at any level \(k\ge3\) when \(\alpha<a^*\approx0.0167788801729\). That closed BN.3's Chebyshev route to \(\delta_N^{(k+1)}>\delta_N^{(k)}\), and BO left open whether the deficit itself still increases there. It does. Chebyshev needs \(\omega_k\) to be increasing. A covariance against an increasing \(G_k\) is already positive if \(\omega_k\) crosses its own mean only once, from below. Below \(1/17\) the non-monotone part of \(\omega_k\) sits at moderate \(y\), where \(\omega_k<8/5\). The mean of \(\omega_k\) is at least \(\rho^{-1}=3\cdot5^{-\alpha}>2.9\), because AQ's corner resolvent increases with the level. So the dip of BO.1 lies entirely below the mean and has no effect on the sign.

**Setting.** As in Theorems BK, BN and BO: \(r_k(s)=\tfrac12v^T(L_k+s)^{-1}v\), \(p(y)=y(y+5)\), \(e=2y/(2+y)\), \(B=1-\tfrac{y(y+4)}{3(y+2)(y+5)}\), \(C=\tfrac53(1+\tfrac{2y}5)e\), \(d_k=B-C\,r_k(p)>0\), \(\omega_k=(1+y/5)^\alpha/d_k\), \(\rho=5^\alpha/3\), \(D_n=\int_0^\infty s^{-\alpha}r_n\,ds\), \(D'_n=\rho^{-n}D_n\), \(\mu_k\) BK.1's level-\(k\) probability measure (positive density in \(y\)), \(G_k=1-\eta_N(p^{k-N+1}(y))\), \(\delta_N^{(k)}=\langle G_k\rangle_{\mu_k}\). **Assumptions used as stated:** AQ (3) (\(r_n\uparrow r_\infty\) pointwise, so \(D_{n+1}>D_n\); re-checked here in exact rationals at sample points, \(k=0\)–\(3\)); BK.1 (the transfer identity \(\delta_N^{(k+1)}-\delta_N^{(k)}=\mathrm{Cov}_{\mu_k}(\omega_k,G_k)/\langle\omega_k\rangle\) with \(\langle\omega_k\rangle=D'_{k+1}/D'_k\), valid for every \(\alpha\in(0,1)\)); BL.1 (\(G_k\) strictly increasing in \(y\), \(\alpha\)-free); BN.1–BN.2 (Stieltjes bounds and the majorant, re-run); BG (flatness for \(\alpha\le\gamma\), only in BP.5). No arb cell, no \(\eta\) enclosure, no AZ input. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BP.**

1. **BP.1 (mean floor, every \(k\), every \(\alpha\)).** \(\langle\omega_k\rangle_{\mu_k}=D'_{k+1}/D'_k=\rho^{-1}D_{k+1}/D_k>\rho^{-1}=3\cdot5^{-\alpha}\).
2. **BP.2 (monotone tail, every \(k\)).** BN.2's majorant \(d_k'-\alpha d_k/(5+y)\le-P_\alpha/(3(y+2)^2(y+5)^2)\) holds for every \(\alpha\ge0\), and \(P_\alpha=y^2(2\alpha y+17\alpha-1)+46\alpha y+40\alpha+10\). For \(\alpha\in(0,1/17]\) put \(Y_\alpha=(1-17\alpha)/(2\alpha)\ge0\). Then \(P_\alpha>0\) on \([Y_\alpha,\infty)\), so \(\omega_k\) is strictly increasing on \([Y_\alpha,\infty)\) at every level \(k\).
3. **BP.3 (the head lies below the mean, every \(k\)).** \(r_k(p)\le1/p\) gives \(d_k\ge B-C/p=\frac{(y+4)(2y+5)}{3(y+2)(y+5)}\), and \(\tfrac85(y+4)(2y+5)-3(y+2)(y+5)=(y^2-y+10)/5>0\). So \(\omega_k<\tfrac85(1+y/5)^\alpha\) for all \(y>0\). On \([0,Y_\alpha]\), \((1+y/5)^\alpha\le5^{-\alpha}(5+Y_\alpha)^\alpha\), with \(5+Y_\alpha=1/(2\alpha)-7/2<1/(2\alpha)\), and \((2\alpha)^{-\alpha}\le e^{1/(2e)}<1.2020<15/8\). Hence \(\omega_k<\tfrac85\cdot\tfrac{15}8\cdot5^{-\alpha}=\rho^{-1}<\langle\omega_k\rangle\) on \([0,Y_\alpha]\).
4. **BP.4 (monotone deficit for every \(\alpha\)).** By BP.2–BP.3, \(\{\omega_k\ge\langle\omega_k\rangle\}=[y_c,\infty)\) for a single \(y_c>Y_\alpha\). With \(m=\langle\omega_k\rangle\),
\[
\mathrm{Cov}_{\mu_k}(\omega_k,G_k)=\int(\omega_k-m)\,(G_k-G_k(y_c))\,d\mu_k>0 ,
\]
since the integrand is \(\ge0\) and is \(>0\) off \(y_c\). So \(\delta_N^{(k+1)}>\delta_N^{(k)}\) for every \(\alpha\in(0,1/17]\), every \(N\ge1\) and every \(k\ge N\). Together with BN.3 (\(\alpha\ge1/17\)): **for every \(\alpha\in(0,1)\) and every \(N\ge1\), \(\delta_N^{(n)}\) is strictly increasing in \(n\ge N\).** This includes the whole range \((0,a^*)\) of Kept Failure BO.1, and also the range \((a_0^*,a^*)\) where only the base level is monotone.
5. **BP.5 (derived, from BG as stated).** For \(\alpha\in(0,\gamma]\), \(u_n\to\bar d\) at every non-corner vertex (BG), so \(\delta_N^{(n)}\uparrow1\) strictly. For \(\alpha\in(\gamma,1)\) the increasing limit is \(\delta_N\) (BI), with the two-sided rate of BN.4.

**Proof of BP.1.** BK.1 gives \(\langle\omega_k\rangle=D'_{k+1}/D'_k\). Since \(r_{k+1}\ge r_k\) pointwise and the two differ on an open set (AQ (3)), \(D_{k+1}>D_k\). \(\square\)

**Proof of BP.2.** In BN.2's chain-rule identity \(d'-\alpha d/(5+y)=B'-\alpha B/(5+y)+r(K_0+\alpha C/(5+y))+Cp'(-r_k'(p))\), the coefficients \(K_0+\alpha C/(5+y)\) and \(Cp'\) are \(\ge0\) for every \(\alpha\ge0\). So BN.1's \(-r'\le r/p\) and \(r\le1/p\) give the majorant for every \(\alpha\ge0\), not only for \(\alpha\ge1/17\). For \(y\ge Y_\alpha\) every term of \(P_\alpha\) is \(\ge0\) and \(40\alpha+10>0\). Then \((\log\omega_k)'=-(d_k'-\alpha d_k/(5+y))/d_k>0\). \(\square\)

**Proof of BP.3.** The two identities are sympy checks, and \(y^2-y+10\) has discriminant \(-39\). \(f(\alpha)=-\alpha\log(2\alpha)\) has \(f''=-1/\alpha<0\) and its only critical point at \(\alpha=1/(2e)\), where \(f=1/(2e)\). The bound \(e^{1/(2e)}\le1.201943\) is an arb enclosure. \(\square\)

**Proof of BP.4.** \(\omega_k\) is continuous. It is \(<m\) on \([0,Y_\alpha]\) (BP.3) and strictly increasing on \([Y_\alpha,\infty)\) (BP.2). Since \(\int(\omega_k-m)\,d\mu_k=0\) and \(\mu_k\) has a positive density, \(\omega_k>m\) somewhere, so \(\omega_k-m\) changes sign exactly once, at some \(y_c>Y_\alpha\). For any constant \(g_0\), \(\mathrm{Cov}(\omega,G)=\langle(\omega-\langle\omega\rangle)(G-g_0)\rangle\). Take \(g_0=G_k(y_c)\). \(G_k\) is strictly increasing (BL.1), so both factors have the same strict sign for \(y\ne y_c\). Divide by \(\langle\omega_k\rangle>0\) (BK.1). \(\square\)

**Evidence (floats; consistency only).** At \(\alpha=0.001,0.005,0.01,0.0167\) and \(k=3,6,10\) (mpmath, 30 digits), \(\omega_k\) is not monotone on \([0,Y_\alpha]\) (BO.1's dip). Its maximum there is \(1.545267\), \(1.550555\), \(1.557624\), \(1.569019\), the same at all three levels, against \(\rho^{-1}=2.995176\), \(2.975955\), \(2.952103\), \(2.920441\). It increases on \([Y_\alpha,10^6Y_\alpha]\). The deficit was computed by BI's quadrature with the heavy \(s^{-1-\alpha}\) tail beyond \(s=10^{12}\) done analytically (\(r_n=1/s-2/s^2+6/s^3-\dots\)). That tail is what stalled BO's float run. At \(\alpha=0.01\), \(N=1\), \(\delta_1^{(n)}=0.9984276\), \(0.9993087\), \(0.9997270\), \(0.9999009\), \(0.9999655\), \(0.9999882\) for \(n=1\)–\(6\), increasing toward \(1\) as BP.5 says. \(N=2\) and \(\alpha=0.001\) behave the same way. The same code reproduces BN's \(\delta_1^{(8)}(0.8)=0.9487\).

**What stays open.**
- *Rates below \(\gamma\).* BP gives monotonicity only. How fast \(1-\delta_N^{(n)}\) goes to \(0\) for \(\alpha<\gamma\) (floats at \(\alpha=0.01\), \(N=1\): successive ratios of \(1-\delta_1^{(n)}\) are \(0.44,0.40,0.36,0.35,0.34\), approaching \(\rho=0.3387\ldots\)) is not proved. (Answered by Theorem BQ below: the order is exactly \(\rho^n\) at every level, and \(\rho^{-n}(1-\delta_N^{(n)})\to C_N>0\).)
- *A closed form for \(a^*\).* Still open. BP shows that \(a^*\) is a threshold for the weight only, not for the deficit.
- *The limit.* Conjecture BL.6's limit part: existence and \(N\)-independence of \(C_\delta(\alpha)\).
- *Other items.* A certificate uniform over an interval of \(\alpha\); closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BP is a monotonicity statement for the corner deficit of one graph Dirichlet extension of a discrete fractional operator on the combinatorial gasket. It does not make \(\omega_k\) monotone below \(a^*\): Kept Failure BO.1 stands. It is about functions on graphs, not a field, a flux or a force. No flux functional is constructed, and Stage 1 is not touched. Not thrust.

**Witness, computed by `scripts/gasket_fractional_deficit_small_alpha.py`** (about 30 s; needs sympy, mpmath, python-flint).
- *Sympy:* the splitting of \(P_\alpha\); \(B-C/p\); the \(8/5\) identity and discriminant; \(5+Y_\alpha\); the concavity of \(-\alpha\log(2\alpha)\); the shifted-covariance identity; BN.2's identity and majorant re-run.
- *Exact rationals:* \(r_{k+1}>r_k\) at five points for \(k=0\)–\(3\); \(B-C/p\le d_k\le B\) and \(d_k>5/8\); BN.1's Stieltjes checks re-run.
- *Arb:* \(e^{1/(2e)}<15/8\).
- *Floats:* the head maxima and tail monotonicity of \(\omega_k\) above; \(\delta_N^{(n)}\) at \(\alpha=0.01,0.001\) and the \(0.8\) control.

The script exits nonzero if any check fails.

## Theorem BQ. Below \(\log3/\log5\) the corner deficit approaches \(1\) at the exact rate \(\rho^n\) at every level: \(c_N\rho^n<1-\delta_N^{(n)}<C'_N\rho^n\) for every \(n\), \(\rho^{-n}(1-\delta_N^{(n)})\to C_N(\alpha)>0\), with \(C_1(0.45)\in[0.0718892,0.0719611]\) (Kept Failure BQ.1: the renormalised deficit is neither monotone nor one-sided)

Theorem BP proved that \(\delta_N^{(n)}\) increases strictly in \(n\) for every \(\alpha\in(0,1)\), to \(1\) when \(\alpha\le\gamma\). It left the rate below \(\gamma\) open. The floats there showed successive ratios of \(1-\delta_1^{(n)}\) of \(0.44,0.40,0.36,0.35,0.34\) at \(\alpha=0.01\), approaching \(\rho=0.3387\ldots\). This theorem gives that rate. The tool is BK.1's transfer applied to the numerator of \(1-\delta_N^{(n)}\) rather than to the ratio. The numerator then increases, because the weight is \(\ge1\) and the neighbour profile \(\eta_N\) is positive. The denominator \(D_n\) increases to a finite limit below \(\gamma\). Squeezing the ratio between the two gives a two-sided bound \(\asymp\rho^n\) at every level, not only in the limit. The limit itself is Theorem BG.3 specialised to the corner neighbour. What is new is the bound at every level, the monotone numerator, and certified constants.

**Setting.** As in Theorems BG, BJ, BK, BL and BP: Assumption 4, \(r_n\) the corner resolvent, \(p(y)=y(y+5)\), \(\rho=5^\alpha/3\), \(\gamma=\log3/\log5\), \(D_n=\int_0^\infty s^{-\alpha}r_n\,ds\), \(D'_n=\rho^{-n}D_n\), \(\mathcal R_n(\tau)=(3/5)^nr_n(5^{-n}\tau)\to\mathcal R\), \(\psi_k(t)=p^k(5^{-k}t)\uparrow\psi\), \(\omega_n\) BK.1's one-level weight, and \(\eta_N\) the neighbour profile (\(\eta_1=1/(2(z+5))\), \(\eta_N=\alpha_1+\beta\,\eta_{N-1}\circ p\)). For \(N\ge1\) and \(n\ge N\), BG.1 and BJ.4 (symmetric data) give
\[
1-\delta_N^{(n)}=\frac{\int_0^\infty s^{-\alpha}r_n(s)\,\eta_N(p^{n-N}s)\,ds}{D_n}=\frac{\rho^n\,\mathcal N_N^{(n)}}{D_n},\qquad
\mathcal N_N^{(n)}=\int_0^\infty\tau^{-\alpha}\mathcal R_n(\tau)\,\eta_N\bigl(\psi_{n-N}(5^{-N}\tau)\bigr)\,d\tau ,
\]
and we write \(V_N^{(n)}=\rho^{-n}(1-\delta_N^{(n)})=\mathcal N_N^{(n)}/D_n\). **Assumptions used as stated (not re-run here):** AQ (3) (\(r_n\uparrow r_\infty\) pointwise); BG.1, BG.3 (the profile and its dominated limit \(P_N=\int\tau^{-\alpha}\mathcal R\,\eta_N(\psi(5^{-N}\tau))\,d\tau\)); BH.2/BJ.4 (\(\eta_N\le\eta\le2/(2+z)\)); BI.4 (\(G=0.6031\), so \(\mathcal R_n\le\min(\tfrac13,G\tau^{\gamma-1})\) and \(r_\infty\le\min(1/s,Gs^{\gamma-1})\)); BJ.1 (\(1/r_k-1/r_\infty=\prod_{i<k}\beta(p^is)(1+2\eta(p^ks))\)); BJ.3 (\(s^{1-\gamma}r_\infty\in[0.600945014,0.600963487]\) on \((0,p(10^{-12})]\)); BK.1–BK.4 (transfer, \(\omega_n\ge(1+y/5)^\alpha\), \(\omega_{\rm hi}\) and its head and tail lemmas); BL.1 (\(\eta_N\) positive and strictly decreasing); AY (only in the derived part of BQ.2). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BQ.**

1. **BQ.1 (the numerator increases; every \(\alpha\in(0,1)\), \(N\ge1\), \(n\ge N\)).**
   \[
   \mathcal N_N^{(n+1)}-\mathcal N_N^{(n)}=\rho^{-n}\int_0^\infty(\omega_n(y)-1)\,\eta_N\bigl(p^{n+1-N}(y)\bigr)\,p^{-\alpha}r_n(p)\,p'\,dy>0 ,
   \]
   so \(\mathcal N_N^{(n)}\) increases strictly to \(P_N<\infty\).
2. **BQ.2 (exact order \(\rho^n\) at every level; \(\alpha<\gamma\)).** For every \(N\ge1\) and \(n>m\ge N\),
   \[
   \frac{\mathcal N_N^{(m)}}{D_\infty}\,\rho^n\;<\;1-\delta_N^{(n)}\;<\;\frac{P_N}{D_m}\,\rho^n ,
   \]
   and \(V_N^{(n)}\to C_N(\alpha)=P_N/D_\infty\in(0,\infty)\). Hence \((1-\delta_N^{(n+1)})/(1-\delta_N^{(n)})\to\rho\), which is BP's float observation. *Derived (AY as stated):* at \(\alpha=\gamma\), \(n(1-\delta_N^{(n)})\to\frac{2\sin\pi\gamma}{\pi}P_N/c_R(\gamma)\).
3. **BQ.3 (certified, \(N=1\); python-flint `arb`, outward).** The table gives \(\mathcal N_1^{(n)}\) for \(n=1,\ldots,6\) and \(30\) (strictly increasing, as BQ.1 says), \(D_n\) for \(n=1,\ldots,6\), \(D_\infty\), and \(E_{30}\ge P_1-\mathcal N_1^{(30)}\):

   | \(\alpha\) | \(D_\infty\) | \(E_{30}\le\) | \(C_1(\alpha)\) | every \(n\ge1\): \(V_1^{(n)}\in\) | every \(n\ge30\): \(V_1^{(n)}\in\) |
   |---|---|---|---|---|---|
   | \(0.01\) | \([99.77849,99.79345]\) | \(2.99\cdot10^{-6}\) | \([0.0078312,0.0078363]\) | \([0.004622,0.007867]\) | \([0.0078312,0.0078365]\) |
   | \(0.1\) | \([10.01691,10.02179]\) | \(3.56\cdot10^{-8}\) | \([0.0606923,0.0607522]\) | \([0.036684,0.063702]\) | \([0.0606923,0.0607864]\) |
   | \(0.3\) | \([4.035386,4.037405]\) | \(4.10\cdot10^{-12}\) | \([0.0933210,0.0934143]\) | \([0.060971,0.116634]\) | \([0.0933210,0.0944433]\) |
   | \(0.45\) | \([4.000301,4.002302]\) | \(9.52\cdot10^{-14}\) | \([0.0718892,0.0719611]\) | \([0.050648,0.118859]\) | \([0.0718892,0.0771286]\) |

   For example, \(0.050648\,\rho^n<1-\delta_1^{(n)}(0.45)<0.118859\,\rho^n\) for every \(n\ge1\), with \(\rho=0.68772568\ldots\). The \(D_\infty(0.45)\) bracket contains Theorem AS's \([4.0013017734819,4.001301773482]\).

**KEPT FAILURE BQ.1 (the renormalised deficit is not monotone and not one-sided).** Claim: \(V_1^{(n)}=\rho^{-n}(1-\delta_1^{(n)})\) is monotone in \(n\), or at least stays on one side of its limit. BP's floats at \(\alpha=0.01\), where the ratios of \(1-\delta_1^{(n)}\) stay above \(\rho\), suggest an increasing sequence. This is false. **Witnesses (certified enclosures):**
- \(\alpha=0.45\): \(V^{(1)}\in[0.083656,0.083739]\), \(V^{(2)}\in[0.086386,0.086473]\), \(V^{(3)}\in[0.084972,0.085057]\). So \(V^{(1)}<V^{(2)}>V^{(3)}\), and \(V^{(2)}>C_1(0.45)\): the limit is approached from above. Further values are \(V^{(4)},V^{(5)},V^{(6)}\in[0.081965,0.082047],[0.079086,0.079165],[0.076874,0.076951]\).
- \(\alpha=0.3\): \(V^{(3)}\le0.092755<0.093321\le C_1\) and \(V^{(4)}\ge0.094315>0.093414\ge C_1\). The sequence crosses its limit.
- \(\alpha=0.01\): \(V^{(1)}\le0.004643<C_1\), and \(V^{(1)}<\cdots<V^{(6)}\) (certified). Here the approach starts from below.

The rate is exactly \(\rho^n\), but the sign of \(V^{(n)}-C_1\) depends on \(\alpha\) and on \(n\). A heuristic reading, not proved: \(V^{(n+1)}/V^{(n)}=(\mathcal N^{(n+1)}/\mathcal N^{(n)})(D_n/D_{n+1})\) is a competition between the numerator's gain, which comes from the region \(y\lesssim5^{N-n}\) where \(\omega_n-1\) is small, and the denominator's gain \(D_{n+1}-D_n\), which is of order \(\rho^n\) (AR). Which one wins first depends on \(\alpha\). This is not proved.

**Proof of BQ.1.**
- *Transfer.* BK.1's substitution \(\tau=5^{n+1}y\), \(\tau'=5^np(y)\) gives \(\tau^{-\alpha}\mathcal R_{n+1}(\tau)\,d\tau=\omega_n(\tau')\,\tau'^{-\alpha}\mathcal R_n(\tau')\,d\tau'\). The neighbour factor of level \(n+1\) at \(\tau\) is \(\eta_N(p^{n+1-N}(y))\), which is the level-\(n\) factor at \(\tau'\). So \(\mathcal N^{(n+1)}=\int\omega_nF_n\,d\nu_n\) and \(\mathcal N^{(n)}=\int F_n\,d\nu_n\), with \(d\nu_n=\tau'^{-\alpha}\mathcal R_n\,d\tau'=\rho^{-n}p^{-\alpha}r_n(p)p'\,dy\) (BK.4).
- *Sign.* \(\omega_n\ge(1+y/5)^\alpha>1\) for \(y>0\) (BK.2). \(F_n>0\), because \(\eta_1=1/(2(z+5))\) and \(\eta_N\ge\alpha_1>1/(2(z+5))\) for \(N\ge2\) by induction (sympy: \(\alpha_1-\tfrac1{2(z+5)}=\tfrac{z+6}{2(z+2)(z+5)}\)).
- *Limit.* \(\mathcal N^{(n)}\to P_N\) is BG.3's dominated convergence, with \(h_N=x_a\eta_N\) at the neighbour. \(P_N\le\tfrac13\int\tau^{-\alpha}\min(1,2\cdot5^N/\tau)\,d\tau<\infty\). \(\square\)

**Proof of BQ.2.** \(r_n\uparrow r_\infty\), with strict increase on an open set (AQ (3)), so \(D_m\le D_n<D_\infty\). For \(\alpha<\gamma\), \(D_\infty=\int s^{-\alpha}r_\infty\le\int s^{-\alpha}\min(1/s,Gs^{\gamma-1})<\infty\). With BQ.1, \(\mathcal N^{(m)}<\mathcal N^{(n)}<P_N\). Divide. The limit follows from \(\mathcal N^{(n)}\to P_N\) and \(D_n\to D_\infty\). \(P_N>0\) because its integrand is positive. At \(\gamma\), \(\rho=1\) and \(D_n/n\to\pi c_R(\gamma)/(2\sin\pi\gamma)\) (AY, as in BG.3). \(\square\)

**Proof of BQ.3.** Everything is computed in python-flint `arb` (64 bits, outward).
- *\(\mathcal N_1^{(n)}\).* The integrand is decreasing in \(\tau\): \(r_n\) is Stieltjes, and \(\eta_1\circ\psi_{n-1}\) is decreasing (BL.1). On geometric cells of ratio \(1.0005\) on \([10^{-16},10^{24}]\), the lower sum uses right endpoints and the upper sum uses left endpoints. \(r_n\) is AQ's recursion along the exact orbit. The head uses \(\mathcal R_n\le\tfrac13\) and \(\eta\le1\). The tail uses \(\mathcal R_n\le G\tau^{\gamma-1}\), \(\eta_N(z)\le2/z\) and \(z\ge5^{-N}\tau\), which gives \(2G5^NT^{\gamma-1-\alpha}/(1+\alpha-\gamma)\).
- *\(D_n\).* Cells of ratio \(1.0005\) on \([10^{-32},10^{16}]\). The head lies between \(r_n(\varepsilon)\) and \(r_n(0)=\tfrac13(5/3)^n\). The tail lies between \(1/(3+s)\) and \(1/s\).
- *\(D_\infty\).* The lower bound uses \(r_\infty\ge r_k\) with \(k=\min(60,\text{orbit cap }10^{100})\) at right endpoints, and \(g\ge0.600945014\) below \(10^{-32}\). The upper bound uses \(r_\infty(a)\le\min\bigl(Ga^{\gamma-1},\,1/a,\,1/L_k(a)\bigr)\) at left endpoints, where \(L_k=1/r_k-\prod_{i<k}\beta(p^ia)(1+4/(2+p^ka))\) (BJ.1 with \(\eta\le2/(2+z)\)) and \(G\) is replaced by \(0.600963487\) below \(p(10^{-12})\). The head is \(g_{\rm hi}\varepsilon^{\gamma-\alpha}/(\gamma-\alpha)\) and the tail is \(T^{-\alpha}/\alpha\).
- *\(E_{30}\).* By BQ.1 and Tonelli, \(P_1-\mathcal N^{(30)}=\sum_{k\ge30}\rho^{-k}\int(\omega_k-1)\eta_1(p^{k}y)\,p^{-\alpha}r_k(p)p'\,dy\). Bound it with \(\omega_k\le\omega_{\rm hi}\) and \(r_k\le\bar r\) (BK.4), and with \(\eta_1(p^ky)\le2/(2+5^ky)\le\min(1,2\cdot5^{-k}/y)\). The \(k\)-sum is geometric in closed form. On cells of ratio \(1.005\) on \([10^{-60},10^6]\), the whole integrand is evaluated on the cell ball. The head uses BK.4's lemma \(\omega_{\rm hi}-1\le c\,y^\gamma\), the bound \(\min(1,x)\le x^\theta\) with \(\theta=(1+\gamma-\alpha)/2\), and \(q=\rho^{-1}5^{-\theta}<1\). The tail uses BK.4's \(Ky^{-\alpha-1}\) and the \(1/y\) decay of the \(k\)-sum.
- *Brackets.* \(C_1\in[\mathcal N^{(30)}_{\rm lo}/D_\infty^{\rm hi},(\mathcal N^{(30)}_{\rm hi}+E_{30})/D_\infty^{\rm lo}]\). For every \(n\ge30\) the bracket is \([\mathcal N^{(30)}_{\rm lo}/D_\infty^{\rm hi},(\mathcal N^{(30)}_{\rm hi}+E_{30})/D_6^{\rm lo}]\). For every \(n\ge1\) it is \([\mathcal N^{(1)}_{\rm lo}/D_\infty^{\rm hi},(\mathcal N^{(30)}_{\rm hi}+E_{30})/D_1^{\rm lo}]\). For \(n\le6\), \(V^{(n)}\in[\mathcal N^{(n)}_{\rm lo}/D_n^{\rm hi},\mathcal N^{(n)}_{\rm hi}/D_n^{\rm lo}]\). \(\square\)

**Evidence (floats; consistency only).** Recursion quadrature at \(0.45\), \(N=1\) (BP's tail-split code): \(V^{(n)}=0.0836976,0.0864295,0.0850147,0.0820056,0.0769125,0.0742613,0.0724387,0.0719507\) at \(n=1,2,3,4,6,8,12,20\). These agree with the certified \(V^{(1)},V^{(2)}\), and \(V^{(20)}\) lies in the certified \(C_1\) bracket. In scratch runs (not in the script), at \(0.3\) the sequence peaks at \(n=4\)–\(5\) and decreases to \(0.0933678\) at \(n=20\), and at \(0.01\) it is still increasing at \(n=20\) (\(0.00783371\)).

**What stays open.**
- *Eventual monotonicity.* Whether \(V_N^{(n)}\) is eventually monotone, and in which direction. Floats suggest eventually decreasing at \(0.3\) and \(0.45\). At \(0.01\) it is still increasing at \(n=20\).
- *Other \(N\) and \(\alpha\).* Closed forms for \(C_N(\alpha)\) and \(P_N\). The \(N\)-dependence of \(C_N\), including its matching with BH's layer profile. Certificates for \(N\ge2\) or uniform in \(\alpha\) (BQ.1–BQ.2 hold for every \(N\) and every \(\alpha<\gamma\), but the numbers are only for \(N=1\) at four exponents).
- *Carried over.* A closed form for \(a^*\); Conjecture BL.6's limit part (existence and \(N\)-independence of \(C_\delta\) above \(\gamma\)); closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BQ is a convergence rate for the corner-neighbour value of one graph Dirichlet extension of a discrete fractional operator on the combinatorial gasket, below the exponent where corners are polar. It is a statement about functions on graphs, not a field, a flux or a force. No flux functional is constructed, and Stage 1 is not touched. Not thrust.

**Witness, computed by `scripts/gasket_fractional_deficit_rate_small_alpha.py`** (about 2.5 min; needs sympy, mpmath, python-flint).
- *Sympy:* \(\eta_1=1/(2(z+5))\); \(\alpha_1-1/(2(z+5))\); the \(\rho^{-n}\) scaling; BK.1's Jacobians; the geometric sums of \(E_M\); \(1+2\eta\le1+4/(2+z)\).
- *Exact rationals:* \(1/(2(z+5))\le\eta_N\le2/(2+z)\), strictly decreasing, for \(N=1\)–\(5\); \(\tfrac35r_{n+1}(y)/r_n(p)=1/(\tfrac53(\beta-er_n(p)))\ge1+2y/5\) for \(n=0\)–\(3\); \(r_{k+1}>r_k\).
- *Arb:* \(\mathcal N_1^{(n)}\), \(D_n\), \(D_\infty\), \(E_{30}\) and all brackets at \(\alpha=0.01,0.1,0.3,0.45\); the Kept Failure BQ.1 witnesses; consistency with AS's \(D_\infty(0.45)\).
- *Floats:* \(V^{(n)}(0.45)\) by quadrature.

The script exits nonzero if any check fails.

## Theorem BR. The corner deficit decreases in the depth \(N\) at every level, and below \(\log3/\log5\) the constant of Theorem BQ grows exactly like \(\rho^{-N}\): \(\rho^NC_N\) strictly increases in \(N\) to a finite limit \(c_\infty(\alpha)\), with \(c_\infty(0.45)\in[0.520387,0.520824]\)

Theorem BQ proved \(\rho^{-n}(1-\delta_N^{(n)})\to C_N(\alpha)=P_N/D_\infty\) for \(\alpha<\gamma\), with certified constants for \(N=1\) only. It left the \(N\)-dependence of \(C_N\) open. This theorem settles its order and its monotonicity. The tool is a one-line inequality for the neighbour profile: going one step deeper in \(N\) while moving the argument one step along the orbit \(p\) can only raise \(\eta\). That gives a monotonicity in \(N\) at every level \(n\) and every \(\alpha\). Below \(\gamma\), the rescaling \(\tau=5^Nt\) turns \(\rho^NP_N\) into an integral whose two factors both increase in \(N\) and are bounded by BJ.2's power bracket. So \(\rho^NC_N\) increases to a finite, positive limit, and that limit is pinned to about four digits.

**Setting.** As in Theorems BG, BH, BJ and BQ: \(p(z)=z(z+5)\), \(\alpha_1=(z+4)/((z+2)(z+5))\), \(\beta=(z+6)/((z+2)(z+5))\), \(\eta_0=-\tfrac12\), \(\eta_N=\alpha_1+\beta\,\eta_{N-1}\circ p\), \(\eta=\lim\eta_N\), \(\rho=5^\alpha/3\), \(\gamma=\log3/\log5\), \(\mathcal R=\lim(3/5)^nr_n(5^{-n}\cdot)\), \(\psi=\lim p^k(5^{-k}\cdot)\), \(g=s^{1-\gamma}r_\infty\). For \(N\ge0\) and \(n\ge N\), \(1-\delta_N^{(n)}=D_n^{-1}\int s^{-\alpha}r_n(s)\,\eta_N(p^{n-N}s)\,ds\) (BG.1, BJ.4, as in BQ's setting; at \(N=0\) this is BL.0's \(\delta_0^{(n)}=\tfrac32\)), \(P_N=\int\tau^{-\alpha}\mathcal R(\tau)\,\eta_N(\psi(5^{-N}\tau))\,d\tau\) and, for \(\alpha<\gamma\) and \(N\ge1\), \(C_N=P_N/D_\infty\) (BQ.2). **Assumptions used as stated (not re-run here):** AQ (3) (\(r_n\uparrow r_\infty\)); BD.4 (\(\psi_k\uparrow\psi\), \(\psi(5u)=p(\psi(u))\)); BG.1, BG.3, BJ.4 (the deficit formulas and \(P_N\)); BH.2 (\(\eta=\lim\eta_N\), \(\eta=\alpha_1+\beta\,\eta\circ p\), and \(\eta=\eta_N+\prod_{i<N}\beta(p^iz)(\tfrac12+\eta(p^Nz))\)); BJ.2 (the bracket \(1/(t^{1-\gamma}/g_{\rm lo}+3(3/5)^N)\le(5/3)^N\mathcal R(5^Nt)\le g_{\rm hi}t^{\gamma-1}\)); BJ.3 (\(g_{\rm lo}=0.600945014\), \(g_{\rm hi}=0.600963487\)); BJ.6's \(\psi(y)\le ye^{y/10}\); BI.4 (\(G=0.6031\)); BQ.2–BQ.3 (\(C_N=P_N/D_\infty\), the \(D_\infty\) routine and the certified \(C_1\)). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BR.**

1. **BR.1 (ceiling; every \(N\ge0\), proved here).** \(\eta_N(z)\le2/(2+z)\) on \([0,\infty)\). (BQ quoted this from BH.2/BJ.4.)
2. **BR.2 (depth step; every \(N\ge0\), \(w\ge0\)).**
   \[
   \eta_{N+1}(w)-\eta_N(p(w))\ \ge\ \frac{w}{w^2+5w+2}.
   \]
3. **BR.3 (monotone in depth at every level; every \(\alpha\in(0,1)\)).** For \(0\le N<n\), with \(w=p^{n-N-1}(s)\),
   \[
   \delta_N^{(n)}-\delta_{N+1}^{(n)}\ \ge\ \frac1{D_n}\int_0^\infty s^{-\alpha}r_n(s)\,\frac{w}{w^2+5w+2}\,ds\ >\ 0 .
   \]
   So \(\tfrac32=\delta_0^{(n)}>\delta_1^{(n)}>\cdots>\delta_n^{(n)}\) at every level. *Derived (BJ.4 as stated):* for \(\alpha>\gamma\) also \(\delta_{N+1}<\delta_N\) for the limit deficit.
4. **BR.4 (\(N\)-dependence of \(C_N\); \(\alpha<\gamma\), \(N\ge1\)).** \(\rho^NC_N\) strictly increases in \(N\) and converges to a finite limit \(c_\infty(\alpha)>0\). For every \(N\ge1\),
   \[
   \rho\,C_1\ \le\ \rho^NC_N\ <\ c_\infty(\alpha),\qquad
   \frac{g_{\rm lo}M_\alpha}{D_\infty}\le c_\infty(\alpha)\le\frac{g_{\rm hi}M_\alpha}{D_\infty},\qquad
   M_\alpha=\int_0^\infty t^{\gamma-1-\alpha}\eta(\psi(t))\,dt\le\frac{2^{\gamma-\alpha}\pi}{\sin\pi(\gamma-\alpha)} .
   \]
   Independently of BQ.3, \(\rho^NC_N\ge L_1/D_\infty\) with \(L_1=\int t^{-\alpha}\eta_1(\psi(t))/(t^{1-\gamma}/g_{\rm lo}+\tfrac95)\,dt\). So \(C_N\) is of exact order \(\rho^{-N}\), and \(C_N\sim c_\infty\rho^{-N}\). *Derived (AT.1 as stated, as in BH.6 and BJ.6):* \(c_\infty(\alpha)=D_\infty^{-1}\int_0^\infty t^{\gamma-1-\alpha}\,\Theta_{AT}(\psi(t))\,\eta(\psi(t))\,dt\). This is the sub-\(\gamma\) analogue of BJ.6's \(c_{\rm corner}\), with \(1-\eta\) replaced by \(\eta\) and \(\kappa_\alpha\) by \(1/D_\infty\).
5. **BR.5 (certified; python-flint `arb`, outward).**

   | \(\alpha\) | \(M_\alpha\) | closed ceiling | \(c_\infty(\alpha)\) | every \(N\ge1\): \(\rho^NC_N\in\) | BQ-free floor \(L_1/D_\infty^{\rm hi}\) |
   |---|---|---|---|---|---|
   | \(0.01\) | \([1.7335799,1.7341000]\) | \(5.846206\) | \([0.010439,0.010444]\) | \([0.002653,0.010444]\) | \(0.001708\) |
   | \(0.1\) | \([1.7583870,1.7589146]\) | \(4.867696\) | \([0.105440,0.105526]\) | \([0.023763,0.105526]\) | \(0.015249\) |
   | \(0.3\) | \([2.1990754,2.1997359]\) | \(4.390927\) | \([0.327320,0.327592]\) | \([0.050414,0.327592]\) | \(0.032401\) |
   | \(0.45\) | \([3.4657841,3.4668569]\) | \(5.530505\) | \([0.520387,0.520824]\) | \([0.049440,0.520824]\) | \(0.032245\) |

   The lower ends of the every-\(N\) column are BQ.3's \(\rho C_1\). At every tabulated \(\alpha\), \(\rho C_1\) lies well below \(c_\infty\) (at \(0.45\): \(\rho C_1\le0.049489<0.520387\le c_\infty\)). So the growth of \(\rho^NC_N\) is real and large, about a factor \(10\) between \(N=1\) and \(N=\infty\) at \(0.45\). For example, at \(0.45\), \(C_N\le0.520824\,\rho^{-N}\) for every \(N\ge1\), and \(C_N\rho^N\to c_\infty\in[0.520387,0.520824]\), with \(\rho=0.68772568\ldots\).

**Proof of BR.1–BR.2.** \(\alpha_1+\beta=2/(2+z)\) (sympy). Since \(\beta>0\), induction from \(\eta_0=-\tfrac12\le1\) gives \(\eta_N(z)\le\alpha_1+\beta\cdot2/(2+p(z))\), and sympy gives \(2/(2+z)-(\alpha_1+2\beta/(2+p))=z(z+6)/((z+2)(z^2+5z+2))\ge0\). For BR.2, \(\eta_{N+1}(w)-\eta_N(pw)=\alpha_1(w)-(1-\beta(w))\,\eta_N(pw)\), with \(1-\beta=(w^2+6w+4)/((w+2)(w+5))>0\). By BR.1 this is \(\ge\alpha_1-(1-\beta)\cdot2/(2+p(w))=w/(w^2+5w+2)\) (sympy). \(\square\)

**Proof of BR.3.** Subtract the two deficit formulas at the same level \(n\). The integrand difference is \(s^{-\alpha}r_n(s)\bigl(\eta_{N+1}(p^{n-N-1}s)-\eta_N(p^{n-N}s)\bigr)\), and BR.2 at \(w=p^{n-N-1}s\) bounds it below. \(r_n>0\), and \(w>0\) for \(s>0\). For \(\alpha>\gamma\), apply the same step to BJ.4's formula at \(w=\psi(5^{-N-1}\tau)\), since \(p(\psi(5^{-N-1}\tau))=\psi(5^{-N}\tau)\). \(\square\)

**Proof of BR.4.** *Scaling.* Put \(\tau=5^Nt\). Then \(\rho^N\tau^{-\alpha}d\tau=(3/5)^Nt^{-\alpha}dt\) (sympy), so
\[
\rho^NP_N=\int_0^\infty t^{-\alpha}\,Q_N(t)\,\eta_N(\psi(t))\,dt,\qquad Q_N(t)=(5/3)^N\mathcal R(5^Nt)=\lim_m(3/5)^mr_{m+N}(5^{-m}t).
\]
*Both factors increase in \(N\).* \(Q_{N+1}\ge Q_N\) because \(r_{k+1}\ge r_k\) pointwise (AQ (3)). With \(\eta=\alpha_1+\beta\,\eta\circ p\), \((\tfrac12+\eta(x))-\beta(x)(\tfrac12+\eta(px))=\tfrac12+\tfrac1{2(x+5)}>0\) (sympy), so the remainder \(\prod_{i<N}\beta(p^iz)(\tfrac12+\eta(p^Nz))\) strictly decreases in \(N\). Hence \(\eta_N(z)\) strictly increases in \(N\) at fixed \(z\), and \(\eta_N\ge\eta_1=1/(2(z+5))>0\) for \(N\ge1\). So \(\rho^NP_N\), and with it \(\rho^NC_N\), strictly increases in \(N\ge1\).
*Bounds.* BJ.2 gives \(Q_N\le g_{\rm hi}t^{\gamma-1}\), and \(\eta_N\le\eta\), so \(\rho^NP_N\le g_{\rm hi}M_\alpha\) for every \(N\). With \(\eta(\psi(t))\le2/(2+t)\) (BR.1 in the limit, \(\psi(t)\ge t\)), \(M_\alpha\le\int t^{c-1}\,2/(2+t)\,dt=2^c\pi/\sin\pi c\), \(c=\gamma-\alpha\in(0,1)\). This is finite exactly because \(\alpha<\gamma\). The increasing sequence is bounded, so it converges. BJ.2's lower bound \(1/(t^{1-\gamma}/g_{\rm lo}+3(3/5)^N)\) increases to \(g_{\rm lo}t^{\gamma-1}\), and \(\eta_N\uparrow\eta\). Monotone convergence gives \(\lim\rho^NP_N\ge g_{\rm lo}M_\alpha\). For \(L_1\), use \(\eta_N\ge\eta_1\) and \(3(3/5)^N\le\tfrac95\).
*Derived part.* BJ.1 gives \(r_\infty/(1+3(3/5)^{m+N}r_\infty)\le r_{m+N}\le r_\infty\), and \((3/5)^mr_\infty(5^{-m}t)\le Gt^{\gamma-1}\) (BI.4, \(5^{1-\gamma}=5/3\)). So for each \(t\), \(Q_N(t)\) lies within a factor \(1+3G(3/5)^Nt^{\gamma-1}\) of \(\lim_m(3/5)^mr_\infty(5^{-m}t)=t^{\gamma-1}\lim_mg(5^{-m}t)\), provided that limit exists. AT.1 says it does and equals \(\Theta_{AT}(\psi(t))\). Dominated convergence, with majorant \(g_{\rm hi}t^{\gamma-1-\alpha}\eta(\psi(t))\), then gives the formula. \(\square\)

**Proof of BR.5.** Everything is computed in python-flint `arb` (64 bits, outward).
- *\(M_\alpha\).* The integrand \(t^{\gamma-1-\alpha}\eta(\psi(t))\) is decreasing. On geometric cells of ratio \(1.0003\) on \([10^{-40},10^4]\), the upper sum uses left endpoints and the lower sum uses right endpoints. \(\psi\) is bracketed by \(p^{14}(5^{-14}t)\le\psi(t)\le p^{14}(5^{-14}te^{5^{-14}t/10})\), because \(\psi=p^K\circ\psi\circ5^{-K}\), \(p\) is increasing, and \(u\le\psi(u)\le ue^{u/10}\). Then \(\eta_{14}\le\eta\le\eta_{14}+\prod_{i<14}\beta(p^iz)(\tfrac12+2/(2+p^{14}z))\), and \(\eta\) is decreasing. The head uses \(\eta\le1\), which gives \(\varepsilon^c/c\). For the tail, \(\psi(t)/t=\prod_{j<k}(p^j(5^{-k}t)+5)/5\) in the limit is nondecreasing. So for \(t\ge T\), \(\eta(\psi(t))\le2T/(\psi_{14}(T)\,t)\), which gives \(2T^c/(\psi_{14}(T)(1-c))\).
- *\(L_1\).* Right-endpoint sums, ratio \(1.0005\) on \([10^{-30},10^4]\), with \(\eta_1(\psi_{\rm hi})\). The head and tail are dropped.
- *\(D_\infty\).* BQ.3's routine, unchanged. Its brackets agree with BQ.3's table.
- *Consistency.* BQ.3's \(\rho C_1\) lies inside the every-\(N\) bracket at all four exponents, as BR.4 requires. \(\square\)

**Evidence (floats; consistency only).** At \(\alpha=0.45\) and level \(n=6\), quadrature gives \(\delta_N^{(6)}=0.9918626,0.9482724,0.8916976,0.8128940,0.7069820\) for \(N=1,\ldots,5\), decreasing as BR.3 says (in the script). In a scratch run (not in the script), BQ's float routine gives \(\rho^NV_N^{(N+14)}(0.45)=0.049579,0.217094,0.314296,0.379079,0.453730\) at \(N=1,2,3,4,6\). These increase toward the certified \(c_\infty(0.45)\in[0.520387,0.520824]\), and the \(N=1\) value is consistent with BQ.3's \(\rho C_1\).

**What stays open.**
- *Closed forms.* A closed form for \(c_\infty(\alpha)\), for \(C_N\) and for \(P_N\). The rate at which \(\rho^NC_N\to c_\infty\) (the floats suggest it is slow, \(0.454\) at \(N=6\) against \(0.520\)). Certified \(C_N\) for individual \(N\ge2\). A bracket uniform in \(\alpha\). The matching of \(c_\infty\) with BH's layer profile, which here is only derived from AT.1.
- *From BQ.* Eventual monotonicity of \(V_N^{(n)}\) in \(n\).
- *Carried over.* A closed form for \(a^*\); Conjecture BL.6's limit part (existence and \(N\)-independence of \(C_\delta\) above \(\gamma\)); closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BR is a statement about how the corner-neighbour value of one graph Dirichlet extension of a discrete fractional operator on the combinatorial gasket depends on the depth of the neighbour. It is a statement about functions on graphs, not a field, a flux or a force. No flux functional is constructed, and Stage 1 is not touched. Not thrust.

**Witness, computed by `scripts/gasket_fractional_deficit_depth.py`** (about 2 min; needs sympy, mpmath, python-flint).
- *Sympy:* \(\alpha_1+\beta=2/(2+z)\); the BR.1 induction gap; the BR.2 identity \(w/(w^2+5w+2)\); \(1-\beta>0\); the monotone remainder \(\tfrac12+1/(2(x+5))\); the \(\tau=5^Nt\) scaling; \(\int t^{c-1}2/(2+t)=2^c\pi/\sin\pi c\).
- *Exact rationals:* \(\eta_N\le2/(2+z)\) for \(N=0\)–\(7\); BR.2 for \(N=0\)–\(6\); \(\eta_{N+1}>\eta_N\) at fixed \(z\).
- *Arb:* \(M_\alpha\), \(L_1\), \(D_\infty\) and all brackets at \(\alpha=0.01,0.1,0.3,0.45\); consistency with BQ.3's \(C_1\).
- *Floats:* \(\delta_N^{(6)}(0.45)\) for \(N=1\)–\(5\).

The script exits nonzero if any check fails.

## Theorem BS. Below \(\log3/\log5\), \(\rho^NC_N\) reaches \(c_\infty\) at the exact rate \(\rho^N\): \(c_\infty-\rho^NC_N=\rho^N\chi_N\) with \(\chi_N\to\chi_\infty\in(0,\infty)\), \(\chi_\infty(0.45)\in[0.633344,0.633869]\); and an exact reciprocal identity \(1/\mathcal R=1/Q_\infty+\Phi\) gives certified \(C_N\) for every \(N=1,\ldots,6\)

Theorem BR proved that \(\rho^NC_N\) increases to a finite limit \(c_\infty(\alpha)\). It left open how fast, and it certified \(C_N\) only through BQ's \(C_1\). This theorem answers both. The tool is a scaled limit of BJ.1. BJ.1 is an exact identity for \(1/r_k-1/r_\infty\) at one level. Rescaled by \(5^{-M}\) and with \(M\to\infty\), it becomes an exact identity between the rescaled resolvent \(\mathcal R\) and the rescaled \(r_\infty\). The correction is an explicit product \(\Phi\) built from \(\psi\), \(\beta\) and \(\eta\). That identity does two things. It writes \(c_\infty-\rho^NC_N\) as \(\rho^N\) times two explicit integrals: one increases to a finite positive limit, and the other decays like \((3/5)^N\). And because \(\Phi\) is explicit, it pins \(\mathcal R\) to the same relative accuracy as BJ.3's \(g\)-bracket, so \(P_N\) and \(C_N\) can be certified directly for every \(N\). No AT.1 input is needed for either.

**Setting.** As in Theorems BG, BH, BJ, BQ and BR: \(p(z)=z(z+5)\), \(\alpha_1=(z+4)/((z+2)(z+5))\), \(\beta=(z+6)/((z+2)(z+5))\), \(\eta_0=-\tfrac12\), \(\eta_N=\alpha_1+\beta\,\eta_{N-1}\circ p\), \(\eta=\lim\eta_N\), \(\rho=5^\alpha/3\), \(\gamma=\log3/\log5\), \(\psi_k(v)=p^k(5^{-k}v)\uparrow\psi(v)\), \(r_\infty=\lim r_n\), \(\mathcal R(u)=\lim_M(3/5)^Mr_M(5^{-M}u)\), \(D_\infty=\lim D_n\). For \(\alpha<\gamma\) and \(N\ge1\), \(P_N=\int_0^\infty u^{-\alpha}\mathcal R(u)\,\eta_N(\psi(5^{-N}u))\,du\), \(C_N=P_N/D_\infty\) (BQ.2), and \(c_\infty=\lim_N\rho^NC_N\) (BR.4). **Assumptions used as stated (not re-run here):** AQ (3) (\(r_n\uparrow r_\infty\); \(r_n\) and \(r_\infty\) nonincreasing); BD.4 (\(\psi_k\uparrow\psi\), \(\psi(5v)=p(\psi(v))\), \(\psi(v)\ge v\)); BJ.6's \(\psi(v)\le ve^{v/10}\); BG.3 (\(\mathcal R=\lim\mathcal R_M\), \(\mathcal R\le\tfrac13\)); BH.2 (\(\eta=\alpha_1+\beta\,\eta\circ p\), \(\eta=\eta_N+\prod_{i<N}\beta(p^iz)(\tfrac12+\eta(p^Nz))\), \(0<\eta\le2/(2+z)\), \(\eta\) strictly decreasing); BJ.1 (\(1/r_k-1/r_\infty=\prod_{i<k}\beta(p^is)(1+2\eta(p^ks))\)); BJ.3 (\(g=s^{1-\gamma}r_\infty\in[0.600945014,0.600963487]\) on \((0,p(10^{-12})]\)); BL.1 (\(\eta_N\) strictly decreasing); BQ.2–BQ.3 (\(C_N=P_N/D_\infty\), the \(D_\infty\) routine); BR.4 (\(\rho^NP_N=\int t^{-\alpha}Q_N(t)\eta_N(\psi(t))\,dt\) with \(Q_N(t)=(5/3)^N\mathcal R(5^Nt)\), and \(c_\infty D_\infty=\lim\rho^NP_N\)). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BS.**

1. **BS.1 (scaled reciprocal identity; every \(u>0\), no \(\alpha\)).** The limit \(Q_\infty(u)=\lim_M(3/5)^Mr_\infty(5^{-M}u)\) exists (AT.1 is not used), and
   \[
   \frac1{\mathcal R(u)}-\frac1{Q_\infty(u)}=\Phi(u):=\bigl(1+2\eta(\psi(u))\bigr)\prod_{k\ge1}\tfrac53\,\beta\bigl(\psi(5^{-k}u)\bigr)\ \in(0,3].
   \]
   So \(\mathcal R=Q_\infty/(1+Q_\infty\Phi)\) exactly. \(Q_\infty\) is nonincreasing, \(Q_\infty(5u)=\tfrac35Q_\infty(u)\) (so \(u^{1-\gamma}Q_\infty(u)\) is multiplicatively \(5\)-periodic), and \(g_{\rm lo}u^{\gamma-1}\le Q_\infty(u)\le g_{\rm hi}u^{\gamma-1}\) for **every** \(u>0\). BR's \(Q_N\) satisfies the exact form of BJ.2's bracket, \(1/Q_N(t)=1/Q_\infty(t)+(3/5)^N\Phi(5^Nt)\), so \(Q_N\uparrow Q_\infty\) and \(c_\infty D_\infty=\int_0^\infty t^{-\alpha}Q_\infty(t)\,\eta(\psi(t))\,dt\) (this replaces BR.4's AT.1-derived formula).
2. **BS.2 (exact correction; every \(\alpha\in(0,\gamma)\), \(N\ge1\)).** \(c_\infty-\rho^NC_N=\rho^N\chi_N\) with \(\chi_N=(A_N+B_N)/D_\infty\),
   \[
   A_N=\int_0^\infty u^{-\alpha}h(u)\,\eta\bigl(\psi(5^{-N}u)\bigr)du,\quad h=Q_\infty-\mathcal R=\frac{Q_\infty^2\Phi}{1+Q_\infty\Phi}>0,\qquad
   B_N=\int_0^\infty u^{-\alpha}\mathcal R(u)\prod_{k=1}^N\beta\bigl(\psi(5^{-k}u)\bigr)\bigl(\tfrac12+\eta(\psi(u))\bigr)du .
   \]
   \(A_N\) strictly increases to \(J=\int_0^\infty u^{-\alpha}h\,du\in(0,\infty)\), and \(0<B_{N+1}\le\tfrac35B_N\). Hence \(\chi_N\to\chi_\infty=J/D_\infty\), and for every \(N\ge N_0\ge1\)
   \[
   \frac{A_{N_0}}{D_\infty}\ \le\ \chi_N\ \le\ \frac{J+B_{N_0}}{D_\infty},\qquad
   C_N=c_\infty\rho^{-N}-\chi_\infty+o(1).
   \]
   So \(\rho^NC_N\) increases to \(c_\infty\) at the exact rate \(\rho^N\): not faster, since \(\chi_N\ge A_1/D_\infty>0\), and not slower, since \(\chi_N\) is bounded. Also \(\chi_N-\chi_\infty\le(3/5)^{N-1}B_1/D_\infty\).
3. **BS.3 (certified; python-flint `arb`, outward).** \(\chi_\infty\), the every-\(N\) brackets, and individual \(C_N\), \(\rho^NC_N\), \(\chi_N\) for \(N=1,\ldots,6\):

   | \(\alpha\) | \(J\) | \(\chi_\infty\) | every \(N\ge1\): \(\chi_N\in\) | every \(N\ge6\): \(\chi_N\in\) | \(C_1\) (BQ.3: independent) |
   |---|---|---|---|---|---|
   | \(0.01\) | \([1.601547,1.601948]\) | \([0.016049,0.016055]\) | \([0.007204,0.031845]\) | \([0.015188,0.016918]\) | \([0.0078322,0.0078352]\) (\([0.0078312,0.0078363]\)) |
   | \(0.1\) | \([1.565803,1.566193]\) | \([0.156240,0.156355]\) | \([0.076806,0.288250]\) | \([0.148669,0.163943]\) | \([0.0607006,0.0607439]\) (\([0.0606923,0.0607522]\)) |
   | \(0.3\) | \([1.770981,1.771427]\) | \([0.438643,0.438973]\) | \([0.267164,0.684584]\) | \([0.423044,0.454603]\) | \([0.0933338,0.0934015]\) (\([0.0933210,0.0934143]\)) |
   | \(0.45\) | \([2.534835,2.535667]\) | \([0.633344,0.633869]\) | \([0.460246,0.858579]\) | \([0.618141,0.649097]\) | \([0.0718991,0.0719511]\) (\([0.0718892,0.0719611]\)) |

   At \(\alpha=0.45\), for \(N=1,\ldots,6\):

   | \(N\) | \(C_N\) | \(\rho^NC_N\) | \(\chi_N\) |
   |---|---|---|---|
   | 1 | \([0.071899,0.071951]\) | \([0.049447,0.049483]\) | \([0.684795,0.685352]\) |
   | 2 | \([0.458110,0.458442]\) | \([0.216671,0.216828]\) | \([0.642187,0.642715]\) |
   | 3 | \([0.964854,0.965556]\) | \([0.313839,0.314068]\) | \([0.635065,0.635589]\) |
   | 4 | \([1.692718,1.693950]\) | \([0.378657,0.378932]\) | \([0.633687,0.634211]\) |
   | 5 | \([2.749351,2.751355]\) | \([0.422967,0.423275]\) | \([0.633413,0.633937]\) |
   | 6 | \([4.285424,4.288550]\) | \([0.453404,0.453735]\) | \([0.633358,0.633882]\) |

   The script prints the same table at \(\alpha=0.01,0.1,0.3\). Every \(C_1\) lies inside BQ.3's independent bracket, which came from a level-\(30\) numerator and \(E_{30}\). The new brackets are tighter. Every \(c_\infty=D_\infty^{-1}\int u^{-\alpha}Q_\infty\eta(\psi(u))\) enclosure meets BR.5's. The identity \(c_\infty D_\infty=\rho^N(P_N+A_N+B_N)\) holds between independent enclosures for \(N=1\)–\(6\) at all four exponents. At \(0.45\), \(A_N+B_N\), and so \(\chi_N\), is certified strictly decreasing for \(N=1,\ldots,5\). The widths of \(C_N\) are set mostly by \(D_\infty\) (relative \(5\cdot10^{-4}\)).
4. **Conjecture BS.4 (floats and BS.3 only).** \(\chi_N\) strictly decreases in \(N\) at every \(\alpha\in(0,\gamma)\), so \(\rho^NC_N\) approaches \(c_\infty\) from below with a monotone renormalised gap. The quadrature gaps \(\chi_N-\chi_\infty(0.45)\approx0.0515,0.0088,0.0017\) at \(N=1,2,3\) shrink by about a factor \(5\) per step. This is faster than the proved \((3/5)^N\) on the \(B_N\) side. No rate for \(J-A_N\) is proved.

*Derived (AT.1 as used in BH.6, as stated).* \(D_\infty=\lim D_n=\int s^{-\alpha}r_\infty=J_r\) by monotone convergence (AQ (3)). If \(u^{1-\gamma}Q_\infty(u)=\Theta_{AT}(\psi(u))\) (BH.6's pointwise limit), then BS.1's formula for \(c_\infty\) is BH.6's formula for the far-field layer constant. So \(c_\infty(\alpha)=\theta_\infty^{\rm BH}(\alpha)\), and BH.6's float \(\theta_\infty(0.45)\approx0.52061\) inherits BR.5's certified bracket \([0.520387,0.520824]\). This settles BR's open "matching with BH's layer" at the derived level only.

**Proof of BS.1.** Put \(s=5^{-M}u\) in BJ.1 and multiply by \((5/3)^M\). Since \(p^i(5^{-M}u)=\psi_i(5^{i-M}u)\), and with \(k=M-i\),
\[
\frac1{\mathcal R_M(u)}-\frac1{A_M(u)}=\bigl(1+2\eta(\psi_M(u))\bigr)\prod_{k=1}^M\tfrac53\,\beta\bigl(\psi_{M-k}(5^{-k}u)\bigr),\qquad
\mathcal R_M(u)=(\tfrac35)^Mr_M(5^{-M}u),\ A_M(u)=(\tfrac35)^Mr_\infty(5^{-M}u).
\]
The product structure is checked in exact rationals: \(1/r_M-1/r_{M+L}=\prod_{i<M}\beta(p^is)(1/r_0-1/r_L)(p^Ms)\) for \(M=0\)–\(3\) and \(L=1,2\). As \(M\to\infty\), \(1/\mathcal R_M\to1/\mathcal R\) (BG.3; \(\mathcal R>0\) by BJ.2's floor). \(\psi_M(u)\uparrow\psi(u)\), and \(\eta\) is continuous, as a uniform limit (\(|\eta-\eta_N|\le\tfrac32(3/5)^N\)). Each factor of the product lies in \((0,1]\). Indeed \(1-\tfrac53\beta(x)=x(3x+16)/(3x^2+21x+30)\in[0,x(3x+16)/30]\) (sympy). For fixed \(k\) the factor tends to \(\phi_k=\tfrac53\beta(\psi(5^{-k}u))\), and it is \(\ge\phi_k\) because \(\psi_j\le\psi\) and \(\beta\) is decreasing. Since \(\psi(5^{-k}u)\le5^{-k}ue^{u/10}\), \(\sum_k-\log\phi_k<\infty\). Dominated convergence on the sums of logarithms gives the limit \(\Phi(u)\). Hence \(1/A_M\) converges. For large \(M\), \(5^{-M}u\le p(10^{-12})\), so BJ.3 and \((3/5)^M(5^{-M}u)^{\gamma-1}=u^{\gamma-1}\) (sympy) give \(A_M(u)\in[g_{\rm lo},g_{\rm hi}]u^{\gamma-1}\). So the limit \(Q_\infty\) exists and obeys the same bracket. \(A_{M+1}(5u)=\tfrac35A_M(u)\) gives the scaling. \(r_\infty\) is nonincreasing, so \(Q_\infty\) is too. \(\Phi\le3\) because every factor is \(\le1\) and \(\eta\le1\). Finally, \(1/Q_N(t)=(3/5)^N/\mathcal R(5^Nt)\), with the identity at \(u=5^Nt\) and the scaling, gives \(1/Q_N=1/Q_\infty+(3/5)^N\Phi(5^N\cdot)\downarrow1/Q_\infty\). Then BR.4's monotone-convergence argument (\(Q_N\uparrow\), \(\eta_N\uparrow\)) gives \(c_\infty D_\infty=\int t^{-\alpha}Q_\infty\,\eta\circ\psi\). \(\square\)

**Proof of BS.2.** Put \(t=5^{-N}u\). Then \(\rho^{-N}t^{-\alpha}(5/3)^Ndt=u^{-\alpha}du\) (sympy). With \(Q_\infty(5^{-N}u)=(5/3)^NQ_\infty(u)\) this gives \(\rho^{-N}c_\infty D_\infty=\int u^{-\alpha}Q_\infty(u)\,\eta(\psi(5^{-N}u))\,du\), and \(P_N=\int u^{-\alpha}\mathcal R(u)\,\eta_N(\psi(5^{-N}u))\,du\) by definition. Subtract, using \(Q_\infty\eta-\mathcal R\eta_N=(Q_\infty-\mathcal R)\eta+\mathcal R(\eta-\eta_N)\). In BH.2's remainder at \(z=\psi(5^{-N}u)\), \(p^iz=\psi(5^{i-N}u)\) (BD.4). This gives \(A_N+B_N\), with \(h=Q_\infty^2\Phi/(1+Q_\infty\Phi)\) from BS.1 (sympy). *Monotonicity.* \(\eta\circ\psi(5^{-N}u)\) strictly increases in \(N\), since \(\eta\) is strictly decreasing and \(\psi\) strictly increasing. It tends to \(\eta(0)=1\), since \(\eta_N(0)=1-\tfrac9{10}(3/5)^{N-1}\) (exact). So \(A_N\uparrow J\) by monotone convergence. \(B_{N+1}\le\beta(0)B_N=\tfrac35B_N\). *Finiteness.* Near \(0\), \(h\le Q_\infty\le g_{\rm hi}u^{\gamma-1}\), integrable against \(u^{-\alpha}\) exactly because \(\alpha<\gamma\). Near \(\infty\), \(h\le Q_\infty^2\Phi\le3g_{\rm hi}^2u^{2\gamma-2}\prod_{k=1,2}\tfrac53\beta(5^{-k}u)\le3g_{\rm hi}^2\tfrac{25}9\cdot125\,u^{2\gamma-4}\) (\(\beta(x)\le1/x\), sympy). \(B_1\le\int u^{-\alpha}\min(\tfrac13,g_{\rm hi}u^{\gamma-1})\min(\tfrac35,5/u)\tfrac32<\infty\). *Bracket.* \(A_{N_0}\le A_N\le J\) and \(B_N\le B_{N_0}\) for \(N\ge N_0\). \(\square\)

**Proof of BS.3.** Everything is computed in python-flint `arb` (64 bits, outward).
- *Pointwise enclosures.* At each grid point \(u\), take \(K\ge6\) with \(v=5^{-K}u\le10^{-12}\). Start from \(\psi(v)\in[v,ve^{v/10}]\) (BD.4, BJ.6) and climb with \(\psi(5x)=p(\psi(x))\). This brackets \(\psi(5^{-k}u)\) for \(k=0,\ldots,K\) with exact \(p\)-chains. \(\eta(\psi(u))\) is bracketed by \(\eta_K\le\eta\le\eta_K+\prod_{i<K}\beta(p^iz)(\tfrac12+2/(2+p^Kz))\), with \(K\) adaptive (at least \(8\) steps past \(p^kz>1000\)). The exact recursion \(\eta=\alpha_1+\beta\,\eta\circ p\) then carries the bracket down the chains to \(\eta(\psi(5^{-N}u))\). \(\eta_N(\psi(5^{-N}u))\) is computed along the same chains from \(\eta_0=-\tfrac12\). \(\Phi\) uses the truncated product on the lower chain (upper bound) and, on the upper chain, the factor \(1-S\) with \(S=(16w/4+3w^2/24)/30\), \(w=ve^{v/10}\), for the omitted \(k>K\) (lower bound; \(\prod(1-c_k)\ge1-\sum c_k\), \(\sum5^{-j}=\tfrac14\), \(\sum25^{-j}=\tfrac1{24}\)). \(Q_\infty\in[g_{\rm lo},g_{\rm hi}]u^{\gamma-1}\). \(h\) is increasing in \(Q\) and \(\Phi\), and \(\mathcal R=1/(1/Q+\Phi)\) is increasing in \(Q\) and decreasing in \(\Phi\) (sympy derivatives). So the bracket ends combine monotonically.
- *Integrals.* All true integrands are positive and nonincreasing in \(u\) (\(Q_\infty\), \(\mathcal R\), \(\Phi\), \(\eta\circ\psi\), \(\eta_N\circ\psi\), \(\beta\circ\psi\)). Upper sums use left endpoints and lower sums right endpoints. The geometric grid has ratio \(1.003\) on \([10^{-20},10^{-6}]\), \(1.0002\) on \([10^{-6},10^8]\) and \(1.003\) on \([10^8,10^{10}]\) (173499 points). *Heads* \((0,10^{-20}]\): \(h\in[Q-\tfrac13,Q]\), \(\mathcal R\le\tfrac13\), \(\prod\beta\le(3/5)^N\). *Tails* \([T,\infty)\): \(h\le3g_{\rm hi}^2\tfrac{25}9\cdot125\,u^{2\gamma-4}\), and \(\psi(v)/v\) is nondecreasing, so \(\psi(5^{-k}u)\ge(u/T)\psi(5^{-k}T)\), with \(\eta,\eta_N\le2/z\) and \(\beta\le1/z\).
- *\(D_\infty\).* BQ.3's routine, unchanged.
- *Consistency.* BQ.3's \(C_1\), BR.5's \(c_\infty\), the identity \(c_\infty D_\infty=\rho^N(P_N+A_N+B_N)\), and BR.4's increase of \(\rho^NP_N\). These are checked as overlaps of independent enclosures. \(\square\)

**Evidence (floats; consistency only).** With \(Q_\infty\approx\langle\Theta\rangle u^{\gamma-1}\) (AU's mean \(0.6009542511\)) and mpmath quadrature at \(\alpha=0.45\): \(c_\infty=0.5206069\), \(\chi_\infty=0.6336061\), \(C_1=0.0719251\), \(C_2=0.4582760\), \(C_3=0.9652051\), and \(\chi_N=0.6850727,0.6424505,0.6353262\) at \(N=1,2,3\) (in the script). A scratch run gave \(\chi_4=0.633949\) and \(\chi_6=0.633620\), so the gaps to \(\chi_\infty\) are \(0.0515,0.0088,0.0017,0.00034\) at \(N=1\)–\(4\), consistent with Conjecture BS.4. BR's scratch \(\rho^NV_N^{(N+14)}(0.45)=0.217094,0.314296\) at \(N=2,3\) sit close to the certified \(\rho^NC_N\) above, as they should for \(n=N+14\).

**What stays open.**
- *Closed forms.* \(u^{1-\gamma}Q_\infty(u)\) is a log-periodic function known only to lie in \([g_{\rm lo},g_{\rm hi}]\). So \(c_\infty\), \(C_N\), \(P_N\) and \(\chi_\infty\) are characterised exactly in terms of one unknown periodic function and explicit \(\Phi\), \(\eta\), \(\psi\), but they have no closed form. A rate for \(J-A_N\), and Conjecture BS.4. Brackets uniform in \(\alpha\): BS.1–BS.2 hold for every \(\alpha\in(0,\gamma)\), but the certificates are pointwise in \(\alpha\). Tighter \(D_\infty\), which now limits every \(C_N\).
- *From BQ.* Eventual monotonicity of \(V_N^{(n)}\) in \(n\).
- *Carried over.* A closed form for \(a^*\); Conjecture BL.6's limit part (existence and \(N\)-independence of \(C_\delta\) above \(\gamma\)); closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BS is a statement about how fast one renormalised corner-neighbour constant of a graph Dirichlet extension of a discrete fractional operator on the combinatorial gasket approaches its depth limit. It also gives an exact identity for one rescaled resolvent. It is a statement about functions on graphs, not a field, a flux or a force. No flux functional is constructed, and Stage 1 is not touched. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_deficit_correction.py`** (about 3 min; needs sympy, mpmath, python-flint).
- *Sympy:* \(1-\tfrac53\beta=x(3x+16)/(3x^2+21x+30)\); \(\beta\le1/x\); the \(u^{\gamma-1}\) scaling of BJ.3's bracket; \(Q_\infty(5u)=\tfrac35Q_\infty(u)\); \(h=Q^2\Phi/(1+Q\Phi)\) and its monotonicity in \(Q\) and \(\Phi\); the \(t=5^{-N}u\) scaling; the \(A_N+B_N\) split; the geometric tail constants; the head integrals; \(\eta_N(0)=1-\tfrac9{10}(3/5)^{N-1}\).
- *Exact rationals:* the two-level product identity \(1/r_M-1/r_{M+L}=\prod_{i<M}\beta(p^is)(1/r_0-1/r_L)(p^Ms)\); \(r_k(0)=(5/3)^k/3\); \(\tfrac53\beta\le1\) and \(\beta\le1/x\).
- *Arb:* \(J\), \(A_N\), \(B_N\), \(P_N\) (\(N=1\)–\(6\)), \(c_\infty D_\infty\) and \(D_\infty\) at \(\alpha=0.01,0.1,0.3,0.45\); all brackets and consistency checks above.
- *Floats:* \(c_\infty\), \(\chi_\infty\), \(C_N\), \(\chi_N\) at \(0.45\), \(N=1\)–\(3\).

The script exits nonzero if any check fails.

## Theorem BT. The log-periodic profile of Theorem BS is AT.1's \(\Theta\) composed with \(\psi\): \(s^{1-\gamma}r_\infty\) is nonincreasing along every \(p\)-orbit, \(u^{1-\gamma}Q_\infty=\Theta\circ\psi\) is continuous with AU's mean, its maximum is the sharp constant \(\sup_{s>0}s^{1-\gamma}r_\infty(s)\in[0.6009615785,0.6009631774]\) (Kept Failure BT.5: it is not constant)

Theorem BS reduced \(c_\infty\), \(C_N\), \(P_N\) and \(\chi_\infty\) to one unknown: the multiplicatively \(5\)-periodic function \(F(u)=u^{1-\gamma}Q_\infty(u)\), known only to lie in BJ.3's bracket \([g_{\rm lo},g_{\rm hi}]\). BS's identification of \(c_\infty\) with BH's far-field constant needed the extra hypothesis \(F=\Theta_{AT}\circ\psi\), taken from AT.1. This theorem proves that hypothesis and pins \(F\) down. The tool is one exact fact: AQ's fixed-point step makes \(g(s)=s^{1-\gamma}r_\infty(s)\) **nonincreasing along every orbit** of \(p\), at every \(s>0\), not only asymptotically. So AT.1's \(\Theta=\lim g\circ q^k\) is a monotone limit, it dominates \(g\), and it is reached at an explicit rate \(O(s^\gamma)\). Then \(F\) can be enclosed pointwise to \(10^{-14}\), and its extremes and its mean can be certified over a full period. It is a closed description of \(F\), not a closed form.

**Setting.** As in Theorems AQ, AT, AU, BD, BJ and BS: \(p(s)=s(s+5)\), \(q=p^{-1}\), \(\gamma=\log3/\log5\) (\(5^\gamma=3\)), \(r_\infty=\lim r_k\), \(g(s)=s^{1-\gamma}r_\infty(s)\), \(Q_\infty(u)=\lim_M(3/5)^Mr_\infty(5^{-M}u)\) (BS.1), \(F(u)=u^{1-\gamma}Q_\infty(u)\), \(\psi(v)=\lim p^k(5^{-k}v)\) (BD.4), \(\Theta(s)=\lim_kg(q^ks)\) (AT.1's profile; its existence is part of BT.1 below), \(\langle\cdot\rangle\) the mean over one period in \(\log_5\). No \(\alpha\) enters. **Assumptions used as stated (not re-run here):** AQ (3) (\(r_k\uparrow r_\infty\) and the fixed-point relation \(r_\infty(s)=\frac{(2+s)(5+s)X}{(6+s)-2p(s)X}\), \(X=r_\infty(p(s))\)); BN.1 (\(r_k\) is Stieltjes of unit mass, so \(r_\infty\) is nonincreasing, \(s\,r_\infty(s)\) is nondecreasing and \(\le1\)); BI.4 (\(g\le G=0.6031\)); BJ.1 (\(1/r_k-1/r_\infty=\prod_{i<k}\beta(p^is)(1+2\eta(p^ks))\)); BH.2 (\(0<\eta\le2/(2+z)\)); BD.4 with BJ.6 (\(\psi(5v)=p(\psi(v))\), \(v\le\psi(v)\le ve^{v/10}\)); BS.1 (\(Q_\infty\) exists, is nonincreasing, \(Q_\infty(5u)=\tfrac35Q_\infty(u)\)). For consistency only: AT.1's eight recorded values, AU's \(\langle\Theta\rangle=0.6009542511\ldots\) and \(|\hat\theta_2|\ge1.34873\cdot10^{-11}\), BJ.3's bracket. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BT.**

1. **BT.1 (orbit monotonicity; every \(s>0\)).** With \(X=r_\infty(p(s))\),
   \[
   \frac{g(s)}{g(p(s))}=\frac{(5+s)^\gamma(2+s)}{6+s-2p(s)X}\ \in\ \Bigl[1,\ \frac{(5+s)^\gamma(2+s)}{6+s-2Gp(s)^\gamma}\Bigr].
   \]
   (the ceiling wherever its denominator is positive, e.g. for \(s\le10^{-3}\)). So \(k\mapsto g(q^ks)\) is nondecreasing, \(\Theta(s)=\lim_kg(q^ks)=\sup_kg(q^ks)\le G\) exists for **every** \(s>0\), \(\Theta\circ q=\Theta\), and
   \[
   g(s)\le\Theta(s)\le g(s)\,e^{E(s)},\qquad E(s)\le\tfrac s4+\tfrac G3\,6^\gamma s^\gamma\quad(s\le10^{-3}).
   \]
   This replaces AT.1's two-sided orbit bound \(\log(g(y)/g(p(y)))\in[-(1-\gamma)y/5,1.82y^{0.65}]\) by a one-sided one with the better exponent \(\gamma\).
2. **BT.2 (identification; no AT.1 input).** For every \(u>0\), \(F(u)=\Theta(\psi(u))\), and \(F(5u)=F(u)\). For every \(v\le10^{-3}\), with \(Y=ve^{v/10}\),
   \[
   e^{-\gamma v/10}\ \le\ \frac{F(v)}{g(v)}\ \le\ e^{(1-\gamma)v/10+Y/4+(G/3)6^\gamma Y^\gamma}.
   \]
   So \(F(u)=F(5^{-M}u)\) is enclosed by one value of \(r_\infty\) at \(5^{-M}u\), \(g(5^{-M}\cdot)\to F\) uniformly on \([1,5]\), and \(F\) is continuous. In particular the hypothesis of BS's derived paragraph holds, and \(c_\infty(\alpha)=\theta_\infty^{\rm BH}(\alpha)\) with BH.6 as stated (no AT.1 limit input).
3. **BT.3 (extremes and mean).** \(\sup_{s>0}g(s)=\max F\): BI.4's constant is sharp at \(G^*=\max F\). For \(s\le10^{-3}\), \(\min F\,e^{-(1-\gamma)s/10-E}\le g(s)\le\max F\,e^{\gamma s/10}\) (this sharpens BJ.3). And \(\langle F\rangle=\langle\Theta\rangle\), AU's mean.
4. **BT.4 (certified; python-flint `arb`, outward).** \(F\) at the \(2^{17}+1\) nodes \(u=5^{i/2^{17}}\) (\(M=30\); widths \(\le8.3\cdot10^{-15}\)), with cell bounds from \(Q_\infty\) nonincreasing and \(uQ_\infty\) nondecreasing:
   \[
   \max F\in[0.6009615785,\,0.6009631774],\quad \min F\in[0.6009453250,\,0.6009469238],\quad \langle F\rangle\in[0.6009530801,\,0.6009554223].
   \]
   The maximum is attained near \(u=5^{0.7518}\) and the minimum near \(u=5^{0.2518}\) (nodes). So \(\sup_{s>0}s^{1-\gamma}r_\infty(s)\le0.6009631774\), against BI.4's \(0.6031\). For every \(s\le p(10^{-12})\), \(g(s)\in[0.600945317,\,0.600963178]\) (the factors \(e^{-(1-\gamma)s/10-E}\), \(e^{\gamma s/10}\) are monotone in \(s\); \(E(p(10^{-12}))\le1.32\cdot10^{-8}\)), inside BJ.3's \([0.600945014,0.600963487]\).

**KEPT FAILURE BT.5 (the profile is not a constant).** The claim "\(Q_\infty(u)=C\,u^{\gamma-1}\) for one constant \(C\)", equivalently \(F\) constant, is false. It is what BS's float evidence uses, with \(C=\langle\Theta\rangle\). Certified witness: \(\max F-\min F\ge1.465480\cdot10^{-5}>0\) (\(F(5^{98534/2^{17}})\ge0.6009615785\), \(F(5^{32998/2^{17}})\le0.6009469238\)). So BS's floats are approximations with relative error up to about \(1.3\cdot10^{-5}\) in \(Q_\infty\); BS's certified brackets never used that approximation and are unaffected.

**Proof of BT.1.** Write \(r=r_\infty(s)\) and \(X=r_\infty(p(s))\), so \(r=(2+s)(5+s)X/((6+s)-2pX)\) (AQ (3)). Since \(s/p(s)=1/(5+s)\), \(g(s)/g(p(s))=(s/p)^{1-\gamma}r/X=(5+s)^{\gamma}(2+s)/(6+s-2pX)\) (sympy). The denominator equals \((2+s)(5+s)X/r>0\). Since \(X>0\), the ratio is at least \((5+s)^\gamma(2+s)/(6+s)\ge3(2+s)/(6+s)\ge1\), because \(5^\gamma=3\) and \(3(2+s)-(6+s)=2s\). For the ceiling, \(pX=p\,r_\infty(p)\le Gp^\gamma\) (BI.4). Then \(g(q^{k+1}s)\ge g(q^ks)\) with \(q^{k+1}s\) in place of \(s\). The sequence is bounded by \(G\), so it converges to its supremum, and \(\Theta(qs)=\Theta(s)\) because the orbits are shifts of each other. For the rate, the ceiling is \((1+s/5)^\gamma(1+s/2)/(1+s/6-x/3)\) with \(x=Gp^\gamma\) (sympy), so \(\log\) of it is at most \(\gamma s/5+s/2+2x/3\) when \(x/3\le\tfrac12\) (\(-\log(1-y)\le2y\) on \([0,\tfrac12]\)). \(p(s)\le6s\) for \(s\le1\). \(\log(\Theta(s)/g(s))=\sum_{k\ge1}\log(g(q^ks)/g(q^{k-1}s))\), the \(k\)-th term is the ratio at \(q^ks\le5^{-k}s\) (\(p(y)\ge5y\)), and \(\sum5^{-k}=\tfrac14\), \(\sum(5^{-k})^\gamma=\sum3^{-k}=\tfrac12\). That gives \(E(s)\le(\gamma/5+\tfrac12)s/4+(G/3)6^\gamma s^\gamma\le s/4+(G/3)6^\gamma s^\gamma\). \(\square\)

**Proof of BT.2.** *Comparison lemma.* For \(0<w\le w'\), \(r_\infty\) nonincreasing and \(s\,r_\infty\) nondecreasing (BN.1) give \((w/w')^\gamma\le g(w')/g(w)\le(w'/w)^{1-\gamma}\). *Identification.* By BS.1 and \((5^{-M}u)^{1-\gamma}=(3/5)^Mu^{1-\gamma}\) (sympy), \(g(5^{-M}u)=u^{1-\gamma}(3/5)^Mr_\infty(5^{-M}u)\to F(u)\). By BD.4, \(q^M(\psi(u))=\psi(5^{-M}u)\in[5^{-M}u,\,5^{-M}ue^{5^{-M}u/10}]\), so the lemma gives \(|\log g(q^M\psi(u))-\log g(5^{-M}u)|\le5^{-M}u/10\to0\). Hence \(\Theta(\psi(u))=\lim_Mg(q^M\psi(u))=F(u)\). *Rate.* \(F(5u)=F(u)\) because \(5^{1-\gamma}\cdot\tfrac35=1\). For \(v\le10^{-3}\), \(F(v)=\Theta(\psi(v))\in[g(\psi(v)),g(\psi(v))e^{E(\psi(v))}]\) by BT.1, and \(g(\psi(v))/g(v)\in[e^{-\gamma v/10},e^{(1-\gamma)v/10}]\) by the lemma with \(\psi(v)\le Y\). The bound on \(E\) is increasing in its argument. *Continuity.* \(g\) is continuous (\(r_\infty\) is Stieltjes), and on \([1,5]\) the bracket for \(F(u)/g(5^{-M}u)\) tends to \(1\) uniformly. *Consequence.* BS's derived paragraph assumed exactly \(F=\Theta_{AT}\circ\psi\). \(\square\)

**Proof of BT.3.** \(F=\lim g(5^{-M}\cdot)\le\sup g\). Conversely, for every \(s>0\), \(g(s)\le\Theta(s)=\lim_kg(q^ks)\) (BT.1), and for \(w=q^ks\le10^{-3}\), BT.2 gives \(g(w)\le F(w)e^{\gamma w/10}\le\max F\,e^{\gamma w/10}\). Letting \(k\to\infty\), \(g(s)\le\max F\). The bracket for \(g(s)\), \(s\le10^{-3}\), is BT.2 read backwards. *Mean.* \(\ln5\langle\Theta\rangle=\int_{q(u)}^u\Theta\,ds/s\) for every \(u\), by \(q\)-invariance (AU's definition). \(u/(5+u)\le q(u)\le u/5\), so the piece over \([q(u),u/5]\) is at most \(G\,u/5\). On \([u/5,u]\), \(\Theta/g\in[1,e^{E(u)}]\) and \(F/g\in[e^{-u/10},e^{u/10+E(Y)}]\). The mean of \(F\) over \([u/5,u]\) is \(\langle F\rangle\) for every \(u\). So \(|\langle\Theta\rangle-\langle F\rangle|=O(u+u^\gamma)\) for every small \(u\), hence it is \(0\). \(\square\)

**Proof of BT.4.** At each node, BJ.1 with BH.2's \(\eta\le2/(2+z)\) gives \(r_k(v)\le r_\infty(v)\le1/(1/r_k(v)-\prod_{i<k}\beta(p^iv)(1+4/(2+p^kv)))\), \(v=5^{-30}u\), with \(k\) six steps past \(p^kv>1\). AQ's recursion is evaluated in `arb` at 96 bits. Then BT.2's bracket applies. *Cells.* For \(u\in[u_i,u_{i+1}]\), \(\log(u_{i+1}/u_i)=h=\ln5/2^{17}\): \(Q_\infty\) nonincreasing gives \(F(u)\le F(u_i)(u/u_i)^{1-\gamma}\), and \(uQ_\infty\) nondecreasing gives \(F(u)\le F(u_{i+1})(u_{i+1}/u)^\gamma\). The two ceilings cross at a value at most \(\max(F_i,F_{i+1})e^{\gamma(1-\gamma)h}\) (sympy), and the floors symmetrically. Integrating the better of the two exponential ceilings (floors) over each cell brackets \(\langle F\rangle\). \(F\) is \(5\)-periodic, so \([1,5]\) covers \((0,\infty)\). \(\square\)

**Evidence (floats; consistency only).** The certified enclosures of \(\Theta(5^{-40}5^{j/8})\), \(j=0,\dots,7\), from BT.1 reproduce AT.1's eight recorded 9-digit values. A trapezoid over the certified node midpoints gives \(\langle F\rangle=0.6009542511688\), equal to AU's \(\langle\Theta\rangle\) to \(10^{-10}\). The Fourier moduli of \(F\) in \(\log_5u\) are \(3.664\cdot10^{-6}\), \(1.349\cdot10^{-11}\) and \(4.2\cdot10^{-16}\) (\(k=1,2,3\)). The second sits on AU's certified floor \(|\hat\theta_2|\ge1.34873\cdot10^{-11}\). That is expected, since \(F/\Theta\to1\) as \(s\to0\), the same argument as for the mean, but it is not proved here for \(k\ge1\). The ratio \(g(s)/g(p(s))\), enclosed from BJ.1 at \(s=10^{-20},\dots,100\), lies in BT.1's interval.

**What stays open.**
- *Closed forms.* \(F\) is now a continuous, explicitly computable log-periodic function, but it has no closed form. So \(c_\infty\), \(C_N\), \(P_N\), \(\chi_\infty\) and \(G^*\) still have none. A proof that \(\hat F_k=\hat\theta_k\) for every \(k\), and whether \(\hat\theta_k\ne0\) for every \(k\) (AU). Narrower brackets for \(\max F\), \(\min F\) need finer cells near the extremes; this was not attempted. BS's widths are set by \(D_\infty\), so BT.3's sharper \(g\)-bracket was not fed back into BS.3.
- *From BS.* Conjecture BS.4 (\(\chi_N\) decreasing for every \(N\)); a rate for \(J-A_N\); uniformity in \(\alpha\).
- *Carried over.* Eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); a continuum equation for \(\bar u\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BT is a statement about one rescaled corner resolvent of a graph Laplacian on the combinatorial gasket: its log-periodic profile, the extremes of that profile, and its mean. It is a statement about functions on graphs, not a field, a flux or a force. No flux functional is constructed, and Stage 1 is not touched. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_log_periodic_profile.py`** (about 20 s; needs sympy, python-flint).
- *Sympy:* the orbit-ratio identity, its floor (\(3(2+s)-(6+s)=2s\)) and ceiling algebra, \(-\log(1-y)\le2y\), \(p(s)\le6s\), the geometric sums, \((5^{-M}u)^{1-\gamma}=(3/5)^Mu^{1-\gamma}\), \(5^{1-\gamma}\cdot\tfrac35=1\), the cell-ceiling crossing.
- *Exact rationals:* AQ's step identity \(r_k(s)(6+s-2p\,r_{k-1}(p))=(2+s)(5+s)r_{k-1}(p)\), the positive denominator, \(0<p\,r_{k-1}(p)\le1\), at \(k=1\)–\(4\) and five values of \(s\).
- *Arb:* the \(2^{17}+1\) node enclosures; \(\max F\), \(\min F\), \(\langle F\rangle\); Kept Failure BT.5's witness; the ratio check; containment in BJ.3 and below BI.4's \(G\); AU's \(\langle\Theta\rangle\) inside \(\langle F\rangle\)'s bracket.
- *Floats:* the trapezoid mean and the Fourier moduli.

The script exits nonzero if any check fails.

## Theorem BU. BT's profile is AU's \(\theta\): \(u^{1-\gamma}Q_\infty=\theta\), so \(\hat F_k=\hat\theta_k\) for every \(k\); a convexity certificate pins \(\max F\), \(\min F\), \(\langle F\rangle\), \(\hat F_1\), \(\hat F_2\) to about \(10^{-12}\) (Kept Failure BU.4: \(F\) is not half-period antisymmetric)

Theorem BT left two items about its profile \(F(u)=u^{1-\gamma}Q_\infty(u)=\Theta(\psi(u))\) open. The first was whether its Fourier modes are AU's \(\hat\theta_k\). BT's float moduli agreed with AU's for \(k=1,2\), but this was not proved. The second was that its extremes were bracketed only to \(1.6\cdot10^{-6}\). Both close here. The first is a matter of notation. AU's coordinate \(\psi_{AU}=\lim5^kq^k\) and BD.4's \(\psi=\lim p^k(5^{-k}\cdot)\) are inverse maps, so BT's \(\Theta\circ\psi\) is AU's \(\Theta\circ\psi_{AU}^{-1}=\theta\), literally the same function. The second comes from using second-order information that BT.4 did not use. \(r_\infty\) is convex, so chords give cell bounds with error \(O(h^2)\) instead of \(O(h)\). This is a sharper description of \(F\), not a closed form.

**Setting.** As in Theorems AU and BT: \(p(s)=s(s+5)\), \(q=p^{-1}\), \(\gamma=\log3/\log5\), \(r_\infty\), \(g(s)=s^{1-\gamma}r_\infty(s)\), \(\Theta=\lim g\circ q^k\) (BT.1), \(F(u)=u^{1-\gamma}Q_\infty(u)\). The two coordinates are \(\psi_{BD}(v)=\lim p^k(5^{-k}v)\) (BD.4; written \(\psi\) in BS and BT) and \(\psi_{AU}(s)=\lim5^kq^k(s)=s/A_\infty(s)\) (AU). AU's profile is \(\theta=\Theta\circ\psi_{AU}^{-1}\), with \(\hat\theta_k=\frac1{\ln5}\int_{x/5}^x\theta(y)y^{-i\omega k}\frac{dy}y=a_k-ib_k\) and \(\omega=2\pi/\ln5\). \(\hat F_k\) is defined in the same way from \(F\), and \(t=\log_5u\). No \(\alpha\) enters. **Assumptions used as stated (not re-run here):** AQ (3) (\(r_k\uparrow r_\infty\)); BN.1 (\(r_k\) Stieltjes of unit mass); BJ.1 with BH.2 (BT.4's node bracket for \(r_\infty\), reused unchanged); BD.4 with BJ.6 (\(\psi_{BD}(5v)=p(\psi_{BD}(v))\), \(v\le\psi_{BD}(v)\le ve^{v/10}\)); AU's setting (\(\psi_{AU}\circ q=\psi_{AU}/5\), \(0\le s(\log A_\infty)'(s)\le s/20\)); BT.1–BT.2 (including the transfer bracket \(e^{-\gamma v/10}\le F(v)/g(v)\le e^{(1-\gamma)v/10+E(Y)}\) for \(v\le10^{-3}\)); BI.4 (\(G=0.6031\)). For consistency only: AU's table (\(\langle\Theta\rangle\), \(a_k\), \(b_k\)) and BT.4's brackets. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BU.**

1. **BU.1 (identification).** \(\psi_{AU}\circ\psi_{BD}=\mathrm{id}\) on \((0,\infty)\), and \(\psi_{AU}\) is strictly increasing, so \(\psi_{AU}^{-1}=\psi_{BD}\) and
   \[
   u^{1-\gamma}Q_\infty(u)=F(u)=\theta(u)\qquad\text{for every }u>0 .
   \]
   Hence \(\hat F_k=\hat\theta_k\) for **every** \(k\in\mathbb Z\) (BT's open item), and AU's table gives \(\hat F_k\ne0\) for \(k=1,2,3,4\). AU (2)'s Mellin residues of \(r_\infty\) at \(1-\gamma+i\omega k\) are exactly the Fourier coefficients of BS's profile \(F\), and BT.3's \(\langle F\rangle=\langle\Theta\rangle\) is the case \(k=0\).
2. **BU.2 (convex cell sandwich).** \(r_\infty\) is convex and decreasing. Take log-equispaced nodes \(v_i\) with node enclosures \([\underline r_i,\bar r_i]\). On \([v_i,v_{i+1}]\), the chord through \((v_i,\bar r_i),(v_{i+1},\bar r_{i+1})\) is a ceiling of \(r_\infty\). The left chord through \((v_{i-1},\bar r_{i-1}),(v_i,\underline r_i)\) and the right chord through \((v_{i+1},\underline r_{i+1}),(v_{i+2},\bar r_{i+2})\), extended into the cell, are floors. Every such line \(c+dv\) has \(c>0>d\), and then \(\phi(v)=v^{1-\gamma}(c+dv)\) is strictly concave (\(\phi''=a v^{a-2}((a-1)c+(a+1)dv)\), \(a=1-\gamma\)). So the maximum of the ceiling over a cell is at an end or at \(v^*=-ac/((a+1)d)\). The minimum of the floor is at an end or at the point where the two floors are switched. \(\int\phi\,v^{-i\omega k}dv/v\) has a closed form for every \(k\). The cell error is \(O(h^2)\).
3. **BU.3 (certified; python-flint `arb`, outward).** One period, \(2^{18}\) cells, \(M=30\), BT.4's node bracket:
   \[
   \max F\in[0.6009615785971,\,0.6009615785983],\qquad \min F\in[0.6009469237907,\,0.6009469237943],
   \]
   \[
   \langle F\rangle\in[0.60095425116719,\,0.60095425116956],\qquad \max F-\min F\in[1.4654802832,\,1.4654807576]\cdot10^{-5}.
   \]
   These are the widths \(1.2\cdot10^{-12}\), \(3.6\cdot10^{-12}\), \(2.4\cdot10^{-12}\), against BT.4's \(1.6\cdot10^{-6}\), \(1.6\cdot10^{-6}\), \(2.3\cdot10^{-6}\). So BI.4's sharp constant is \(G^*=\sup_{s>0}s^{1-\gamma}r_\infty(s)=\max F\in[0.6009615785971,0.6009615785983]\). Every Fourier coefficient of \(F\) is enclosed with radius \(1.19\cdot10^{-12}\):
   \[
   a_1\in[4.042159,\,4.042397]\cdot10^{-8},\quad b_1\in[-3.663479,\,-3.663477]\cdot10^{-6},\quad a_2\in[-1.465995,\,-1.228463]\cdot10^{-11}.
   \]
   So \(\hat F_1\ne0\) and \(\hat F_2\ne0\) are certified from \(F\)'s own nodes at depth \(5^{-30}\), independently of AU's computation (AU used interpolation of \(r_\infty\) at depth \(5^{-140}\)). AU's certified \(a_1,b_1,a_2,b_2\) lie inside these enclosures, which is a numerical check of BU.1. The modes \(k=3,4\) are below the radius and are not resolved.

**KEPT FAILURE BU.4 (\(F\) is not half-period antisymmetric).** The claim is "\(F(t+\tfrac12)=2\langle F\rangle-F(t)\)", that is, the oscillation is odd about its mean under a half-period shift. It is suggested by BT.4's extremal nodes, \(98534-32998=2^{16}\), exactly half a period apart. The same relation holds here: \(197068-65996=2^{17}\). Every profile with only odd harmonics satisfies it, including AU.1's single cosine. The claim is false. It would force \(\max F+\min F=2\langle F\rangle\). Certified witness: \((\max F+\min F)/2-\langle F\rangle\in[2.4340,\,2.9084]\cdot10^{-11}>0\), from BU.3's three brackets alone. It is also implied by BU.1 with AU's \(|\hat\theta_2|>0\), since the relation kills every even mode. The witness here does not use AU.

**Proof of BU.1.** By BD.4, \(q(\psi_{BD}(5w))=\psi_{BD}(w)\), so \(q^k(\psi_{BD}(v))=\psi_{BD}(5^{-k}v)\in[5^{-k}v,\,5^{-k}ve^{5^{-k}v/10}]\) (BJ.6). Hence \(5^kq^k(\psi_{BD}(v))\to v\), that is, \(\psi_{AU}(\psi_{BD}(v))=v\). In particular \(\psi_{AU}\) is onto. It is strictly increasing because \(d\log\psi_{AU}/d\log s=1-s(\log A_\infty)'(s)\ge1-s/20>0\) on \((0,20)\) (AU's bound; \(\sum_{i\ge1}5^{-i}/5=\tfrac1{20}\), sympy). For general \(s<s'\), pick \(k\) with \(q^ks'<20\). \(q^k\) is strictly increasing and \(\psi_{AU}=5^k\psi_{AU}\circ q^k\). So \(\psi_{AU}^{-1}=\psi_{BD}\) and \(\theta=\Theta\circ\psi_{BD}=F\) by BT.2. The Fourier statements follow from the definitions, which coincide with \(y=u\). AU (2)'s residue at \(z_k\) is \(\hat\theta_{-k}=\hat F_{-k}\). \(\square\)

**Proof of BU.2.** Each \(r_k=\int d\mu(\lambda)/(\lambda+s)\) has \(r_k'<0<r_k''\) (sympy: \(\partial_s^2(\lambda+s)^{-1}=2(\lambda+s)^{-3}\)). Pointwise limits keep convexity and monotonicity, so \(r_\infty\) is convex and nonincreasing. Convexity puts the chord above \(r_\infty\) on its interval and the extended chords below \(r_\infty\) outside theirs. The extended chord equals \(r_i(1+\lambda)-r_{i-1}\lambda\) with \(\lambda\ge0\), so it is increasing in \(r_i\) and decreasing in \(r_{i-1}\). That justifies the choice of floor and ceiling endpoints (sympy). The signs \(c>0>d\) are checked in `arb` on every cell. Concavity of \(\phi\) on a subinterval puts its minimum at the ends. Using any one of the two floors on each of the two subintervals is valid, because \(r_\infty\) dominates both. The antiderivative \(cv^z/z+dv^{z+1}/(z+1)\) of \(v^{z-1}(c+dv)\) is sympy-checked, with \(z=a\) for the mean and \(z=a-i\omega k\) for \(\hat F_k\). \(u^{-i\omega k}=v^{-i\omega k}\) for \(u=5^Mv\), because \(e^{-2\pi ikM}=1\). BT.2's transfer factors are monotone in \(v\), so they are taken at the right end of each cell. In the Fourier integrals they enter as an additive slack \(G(f_u+f_d)h\) per cell, using that the cell ceiling is \(<G\) (checked). \(\square\)

**Proof of BU.3.** All arithmetic is outward-rounded at 96 bits. The node enclosures are BT.4's: BJ.1 with BH.2's \(\eta\le2/(2+z)\) at \(v=5^{-30}u\), relative width \(\le4.2\cdot10^{-15}\). \(\max F\) lies between the best node floor and the largest cell ceiling, and \(\min F\) lies between the smallest cell floor and the best node ceiling. \(\langle F\rangle=\frac1{\ln5}\int_{\rm period}F\,dv/v\) lies between the integrated floors and ceilings. For \(\hat F_k\), write \(F=\phi_{\rm mid}+e\) with \(|e|\le(\phi_{\rm hi}-\phi_{\rm lo})/2+\)slack. Then \(\hat F_k=\frac1{\ln5}\int\phi_{\rm mid}v^{-i\omega k}dv/v\) (closed form) \(\pm\frac1{\ln5}\int|e|\), and the second term is the common radius \(1.188\cdot10^{-12}\). \(\square\)

**Evidence (floats; consistency only).** The midpoints of the certified \(a_2,b_2\) are \(-1.3472288\cdot10^{-11}\) and \(6.376\cdot10^{-13}\), close to AU's values. The midpoints for \(k=3,4\) do **not** match AU (\(a_3\approx-7.3\cdot10^{-17}\) against AU's \(-1.66\cdot10^{-17}\)). They lie far below the certified radius and mean nothing. BT's float modulus \(4.2\cdot10^{-16}\) for \(k=3\) was noise of the same kind. The certified value is AU's \(|\hat\theta_3|\ge8.82462\cdot10^{-17}\), which by BU.1 is \(|\hat F_3|\). A float sum of the first two harmonics predicts \((\max F+\min F)/2-\langle F\rangle\approx2\,\mathrm{Re}(\hat\theta_2e^{4\pi it^*})\approx2.7\cdot10^{-11}\), inside BU.4's bracket. The consistency step \(5^{40}q^{40}(\psi_{60}(v))=v\) is checked to \(2\cdot10^{-26}\), but it repeats the proof and is not an independent check.

**What stays open.**
- *From BT/BU.* A closed form of \(F=\theta\), and so of \(c_\infty\), \(C_N\), \(P_N\), \(\chi_\infty\), \(G^*\). Whether \(\hat\theta_k\ne0\) for every \(k\ge5\), and the decay rate of \(|\hat\theta_k|\) (floats suggest \(e^{-2\pi^2k/\ln5}\), not proved). The exact location of the extremes: they sit half a period apart to node resolution, but BU.4 shows the profile is not antisymmetric. BS's widths are set by \(D_\infty\), and BU.3 was not fed back into BS.3.
- *From BS.* Conjecture BS.4 (\(\chi_N\) decreasing for every \(N\)); a rate for \(J-A_N\); uniformity in \(\alpha\).
- *Carried over.* Eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BU identifies two descriptions of one log-periodic function built from a corner resolvent of a graph Laplacian on the combinatorial gasket. It certifies that function's extremes, mean and first two Fourier modes. These are statements about functions on graphs, not a field, a flux or a force. No flux functional is constructed, and Stage 1 is not touched. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_fourier.py`** (about 65 s; needs sympy, python-flint; imports BT's script for the node bracket).
- *Sympy:* concavity of \(v^a(c+dv)\) and its stationary point; the closed-form antiderivative for complex exponents; Stieltjes monotonicity/convexity; monotonicity of the extended chord in its endpoints; \(\sum5^{-i}/5=1/20\); \(e^{-2\pi ikM}=1\); even modes are killed by the half-period antisymmetry.
- *Exact rationals:* \(r_k\) strictly decreasing with increasing chord slopes, \(k=0\)–\(4\), at 15 points of \([1/100,40]\); \(\psi_n(\tau)=p^n(5^{-n}\tau)\ge\tau\), increasing in \(n\le5\).
- *Arb:* the \(2^{18}+3\) node enclosures; the signs \(c>0>d\) on every cell; the cell ceiling \(<G\); \(\max F\), \(\min F\), \(\langle F\rangle\), \(a_k,b_k\) (\(k\le4\)); containment in BT.4's brackets; AU's \(\langle\Theta\rangle\), \(a_1,b_1,a_2,b_2\) inside BU.3's enclosures; Kept Failure BU.4's witness.
- *Floats:* none are used in any check.

The script exits nonzero if any check fails.

## Theorem BV. The profile \(F=\theta\) is holomorphic and bounded on the strip \(|\operatorname{Im}t|<\pi/\ln5\), so \(|\hat\theta_k|\le B_k\le eG^*(1+4\pi|k|/\ln5)\,e^{-2\pi^2|k|/\ln5}\); \(\hat\theta_5,\dots,\hat\theta_{11}\ne0\); \(F\) has exactly one maximum and one minimum per period, located to \(10^{-54}\) (Kept Failure BV.5: they are not half a period apart)

**Setting.** As in Theorems BS–BU: \(r_k(s)=\tfrac12v^T(L_k+s)^{-1}v\), \(r_\infty=\lim r_k\), \(Q_\infty(u)=\lim_M(3/5)^Mr_\infty(5^{-M}u)\), \(\gamma=\log3/\log5\), \(F(u)=u^{1-\gamma}Q_\infty(u)=\theta(u)\), \(F(5u)=F(u)\). Write \(t=\log_5u\), so \(F\) is a \(1\)-periodic function of \(t\), and \(\hat\theta_k=\int_0^1F(5^t)e^{-2\pi ikt}\,dt=a_k-ib_k\) (AU's normalisation), \(F=\langle F\rangle+2\sum_{k\ge1}(a_k\cos2\pi kt+b_k\sin2\pi kt)\). \(G^*=\max F\). **Assumptions:** AQ (3) (\(r_k\uparrow r_\infty\) on \((0,\infty)\)); BN.1 (\(r_k=\sum_iw_i/(\lambda_i+s)\), \(w_i,\lambda_i\ge0\)); BS.1 (existence of \(Q_\infty\)); BT.2 (\(F(5u)=F(u)\)); BU.1 (\(F=\theta\)); BU.3 (\(G^*\le0.6009615785983\)); AU (3)'s certified routine (AQ (4) node enclosures, AS (2) interpolation sandwich, the AT.1 orbit lemma, AU's \(\psi\) and endpoint corrections), imported unchanged and re-run with more cells. All as stated. For consistency only: AU's table and BU.3/BU.4's brackets. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BV.**

1. **BV.1 (sector bound).** Every Stieltjes function \(f(z)=c+\int d\mu(\lambda)/(\lambda+z)\), \(c\ge0\), \(\mu\ge0\) on \([0,\infty)\), satisfies \(|f(z)|\le f(|z|)/\cos(\varphi/2)\) for \(z=|z|e^{i\varphi}\), \(|\varphi|<\pi\).
2. **BV.2 (continuation).** \(r_\infty\) and \(Q_\infty\) continue holomorphically to \(\mathbb C\setminus(-\infty,0]\) with the same bound. So \(F(t)=5^{(1-\gamma)t}Q_\infty(5^t)\) is holomorphic and \(1\)-periodic on \(S=\{|\operatorname{Im}t|<\pi/\ln5\}\), and
   \[
   |F(x+iy)|\le\frac{F(x)}{\cos(y\ln5/2)}\le\frac{G^*}{\cos(y\ln5/2)} .
   \]
   In particular \(F\) is real-analytic.
3. **BV.3 (decay of the modes; BU's open "decay rate", upper half).** For every \(k\ne0\), with \(T_k=4\pi|k|/\ln5\),
   \[
   |\hat\theta_k|\le B_k=G^*\sqrt{1+T_k^2}\,e^{-T_k\arctan T_k}\le e\,G^*(1+T_k)\,e^{-2\pi^2|k|/\ln5},
   \]
   so \(\limsup|\hat\theta_k|^{1/|k|}\le e^{-2\pi^2/\ln5}=4.7155\cdot10^{-6}\). Equivalently, AU (2)'s residues of the Mellin transform of \(r_\infty\) at \(1-\gamma\pm i\omega k\) decay at least at this rate. This is exactly the rate the floats of BT/BU suggested.
4. **BV.4 (certified; AU's routine, outward-rounded `mpmath.iv` at 170 digits, 384 cells, sandwich \(b=u=20/21\), depth \(5^{-140}\), node cap \(10^{130}\)).** \(\langle F\rangle\) to radius \(4.1\cdot10^{-63}\), and \(a_k,b_k\) for \(k\le11\) to radius \(8.2\cdot10^{-63}\):

| \(k\) | \(a_k\) | \(b_k\) | \(\lvert\hat\theta_k\rvert\ge\) | \(B_k\) | \(\lvert\hat\theta_k\rvert e^{2\pi^2k/\ln5}\) |
|---|---|---|---|---|---|
| 5 | \(1.43682897682\cdot10^{-27}\) | \(-9.05475087838\cdot10^{-28}\) | \(1.698341\cdot10^{-27}\) | \(1.4870\cdot10^{-25}\) | \(0.728\) |
| 6 | \(-3.25583867317\cdot10^{-33}\) | \(-7.98359194922\cdot10^{-33}\) | \(8.621962\cdot10^{-33}\) | \(8.4143\cdot10^{-31}\) | \(0.784\) |
| 7 | \(-4.31465378795\cdot10^{-38}\) | \(-2.43242062331\cdot10^{-39}\) | \(4.321505\cdot10^{-38}\) | \(4.6289\cdot10^{-36}\) | \(0.834\) |
| 8 | \(-5.50475216207\cdot10^{-44}\) | \(1.55340459757\cdot10^{-43}\) | \(1.648056\cdot10^{-43}\) | \(2.4945\cdot10^{-41}\) | \(0.674\) |
| 9 | \(8.38526711281\cdot10^{-49}\) | \(3.79766118274\cdot10^{-49}\) | \(9.205158\cdot10^{-49}\) | \(1.3233\cdot10^{-46}\) | \(0.799\) |
| 10 | \(4.12559211982\cdot10^{-54}\) | \(-2.89164189910\cdot10^{-54}\) | \(5.038065\cdot10^{-54}\) | \(6.9334\cdot10^{-52}\) | \(0.927\) |
| 11 | \(-3.03451509014\cdot10^{-60}\) | \(-1.78263934544\cdot10^{-59}\) | \(1.807336\cdot10^{-59}\) | \(3.5964\cdot10^{-57}\) | \(0.706\) |

   So \(\hat\theta_k\ne0\) for every \(k\le11\) (AU: \(k\le4\)), and \(M\) has genuine simple poles at \(1-\gamma\pm i\omega k\) for these \(k\). The rows \(k\le4\) reproduce AU's 22-digit table. \(|\hat\theta_k|\le B_k\) holds for every \(k\le11\), with \(|\hat\theta_k|/B_k\) falling from \(0.061\) (\(k=1\)) to \(0.0050\) (\(k=11\)). That is a check of BV.3 against an independent certificate.
5. **BV.5 (extremes; derived from BV.3–BV.4, python-flint `arb`, 400 bits).** Let \(S\) be the certified \(11\)-mode sum. BV.3's tail gives \(|F-S|\le3.70\cdot10^{-62}\), \(|F'-S'|\le2.79\cdot10^{-60}\), \(|F''-S''|\le2.10\cdot10^{-58}\) on the real line. On every cell of a \(4096\)-cell partition of one period, \(F'\ne0\) or \(F''\ne0\) is certified, and \(F'\) changes sign on exactly two cells. So \(F\) has **exactly one maximum and one minimum per period**, and both are nondegenerate. Bisection on the certified sign of \(F'\) gives
   \[
   t_{\max}=0.75175596179031755877802933273820825656220192611655\ldots\ (\pm9.1\cdot10^{-56}),\qquad t_{\min}=0.25175612419221392884930287700067310106574561482\ldots\ (\pm3.5\cdot10^{-55}),
   \]
   \[
   G^*=\max F=0.60096157859710108904205618413461692677004052591725819702\ldots\ (\pm4.5\cdot10^{-61}),
   \]
   \[
   \min F=0.60094692379426561210422618817212924841162969178767769700\ldots\ (\pm6.7\cdot10^{-61}).
   \]
   Both values lie inside BU.3's independent convexity brackets (widths \(10^{-12}\)), and both locations lie within \(0.36\) of a cell of BU.3's extremal nodes \(197068/2^{18}\) and \(65996/2^{18}\). Also \(\max F-\min F=1.465480283547693782999\ldots\cdot10^{-5}\), and \((\max F+\min F)/2-\langle F\rangle=2.690974274943603605\ldots\cdot10^{-11}\), inside BU.4's \([2.4340,2.9084]\cdot10^{-11}\). So BI.4's sharp constant \(\sup_{s>0}s^{1-\gamma}r_\infty(s)=G^*\) (BT.3) is now known to 60 digits, and \(F\) attains it exactly at \(u=5^{t_{\max}+j}\), \(j\in\mathbb Z\).

**KEPT FAILURE BV.5 (the extremes are not half a period apart).** The claim is "\(t_{\min}=t_{\max}+\tfrac12\pmod1\)". It is suggested by BU.3's extremal nodes, which differ by exactly \(197068-65996=2^{17}\), half the \(2^{18}\)-node period (and by BT.4's \(2^{16}\) of \(2^{17}\)). It is weaker than BU.4's antisymmetry, which BU.4 already refuted, so BU.4 does not decide it. The claim is false. Certified witness: \(t_{\max}-t_{\min}-\tfrac12=-1.6240189637007127\ldots\cdot10^{-7}\ne0\) (radius below \(10^{-54}\)). That is \(0.043\) of a BU cell (\(2^{-18}=3.8\cdot10^{-6}\)), which is why both node grids put the extremes exactly half a period apart.

**CONJECTURE BV.6 (sharp rate).** \(|\hat\theta_k|\asymp e^{-2\pi^2|k|/\ln5}\), that is, \(\lim|\hat\theta_k|^{1/|k|}=e^{-2\pi^2/\ln5}\), and the factor \(1+T_k\) in \(B_k\) is not attained. Equivalently, \(F\) has no holomorphic continuation to any wider strip, so \(Q_\infty\) does not continue across \((-\infty,0)\) from both sides compatibly. **Evidence (certified numbers, conjectural reading):** \(|\hat\theta_k|e^{2\pi^2k/\ln5}\in[0.60,0.93]\) for \(k=1,\dots,11\), and the consecutive ratios \(|\hat\theta_{k+1}/\hat\theta_k|\) range over \(3.6\)–\(6.5\cdot10^{-6}\) around \(e^{-2\pi^2/\ln5}=4.72\cdot10^{-6}\). **Falsifier:** a certified \(|\hat\theta_k|\le e^{-ck}\) with \(c>2\pi^2/\ln5\) for all large \(k\), or a holomorphic continuation of \(F\) past \(|\operatorname{Im}t|=\pi/\ln5\).

**Proof of BV.1.** For \(\lambda\ge0\), \(R=|z|\): \(|\lambda+z|^2-(\lambda+R)^2\cos^2(\varphi/2)=(1-\cos\varphi)(\lambda-R)^2/2\ge0\) (sympy). So \(|1/(\lambda+z)|\le1/((\lambda+R)\cos(\varphi/2))\). Integrate against \(\mu\), and note \(|c|=c\le c/\cos(\varphi/2)\). \(\square\)

**Proof of BV.2.** By BN.1 each \(r_k\) is Stieltjes, so \(|r_k(z)|\le r_k(|z|)/\cos(\varphi/2)\le r_\infty(|z|)/\cos(\varphi/2)\) (AQ (3): \(r_k\le r_\infty\)). On a compact \(K\subset\mathbb C\setminus(-\infty,0]\), \(|z|\ge a>0\) and \(|\varphi|\le\varphi_0<\pi\), and \(r_\infty(|z|)\le r_\infty(a)\) because \(r_\infty\) is nonincreasing. So \(\{r_k\}\) is locally bounded, and it converges on \((0,\infty)\). By Vitali's theorem, \(r_k\to r_\infty\) locally uniformly on the slit plane, the limit is holomorphic, and the bound passes to the limit. The same argument applied to \(h_M(z)=(3/5)^Mr_\infty(5^{-M}z)\) (same argument \(\varphi\), \(h_M(|z|)\le h_M(a)\to Q_\infty(a)\), BS.1) continues \(Q_\infty\) with \(|Q_\infty(z)|\le Q_\infty(|z|)/\cos(\varphi/2)\). With \(z=5^t\), \(t=x+iy\): \(|z|=5^x\), \(\varphi=y\ln5\), and \(|5^{(1-\gamma)t}|=5^{(1-\gamma)x}\) (sympy). So \(|F(x+iy)|\le F(x)/\cos(y\ln5/2)\). \(F(t+1)-F(t)\) is holomorphic on \(S\) and vanishes on \(\mathbb R\) (BT.2), so it vanishes on \(S\). \(\square\)

**Proof of BV.3.** \(F\) is real on \(\mathbb R\), so \(\hat\theta_{-k}=\overline{\hat\theta_k}\); take \(k\ge1\) and \(0<y<\pi/\ln5\). Integrate \(F(t)e^{-2\pi ikt}\) around the rectangle \(0\to1\to1-iy\to-iy\to0\). The vertical sides cancel by periodicity, so \(\hat\theta_k=\int_0^1F(x-iy)e^{-2\pi ik(x-iy)}dx\). Since \(|e^{-2\pi ik(x-iy)}|=e^{-2\pi ky}\), BV.2 gives \(|\hat\theta_k|\le G^*e^{-2\pi ky}/\cos(y\ln5/2)\). Put \(f=y\ln5/2\in(0,\pi/2)\). The exponent is \(-T_kf\), and \(-T_kf-\log\cos f\) is convex with stationary point \(\tan f=T_k\) (sympy). That gives \(B_k\). Finally \(\arctan T=\pi/2-\arctan(1/T)\ge\pi/2-1/T\) and \(\sqrt{1+T^2}\le1+T\) give the second form, because \((\pi/2)T_k=2\pi^2k/\ln5\). \(\square\)

**Proof of BV.4.** This is AU (3)'s proof, unchanged except for its parameters. The script imports `period_integrals` and `corrections` from `scripts/gasket_fractional_log_mean.py` and sets 384 cells per period, interpolation orders \(b=u=20\) (floor, even) and \(21\) (ceiling, odd), as AS (2) allows for any \(b\), \(u\), 170 digits, and node cap \(10^{130}\). The node enclosures have relative width \(<10^{-100}\). A nonzero lower bound on \(|a_k|\) or \(|b_k|\) gives \(\hat\theta_k\ne0\). \(\square\)

**Proof of BV.5.** For real \(t\), \(F^{(j)}(t)=S^{(j)}(t)+2\,\mathrm{Re}\sum_{k\ge12}(2\pi ik)^j\hat\theta_ke^{2\pi ikt}\). Termwise differentiation is valid by BV.3's exponential decay. The tail is at most \(2\sum_{k\ge12}(2\pi k)^jB_k\): the terms up to \(k=311\) are explicit, and the remainder uses \(B_k\le eG^*(1+ak)x^k\), \(x=e^{-2\pi^2/\ln5}\), \(a=4\pi/\ln5\), with consecutive-term ratio \(\le((K+1)/K)^{j+1}x\). \(S^{(j)}\) is evaluated with `acb` on ball arguments. On a cell where \(F''\ne0\), \(F'\) is strictly monotone, so it has one zero if its certified end signs differ and none otherwise. Cells with \(F'\ne0\) have no zero. This gives exactly two critical points per period. The sign of \(F''\) tells which is the maximum, and bisection keeps the certified sign change. \(\max F\) and \(\min F\) are \(S\) on the final balls \(\pm3.70\cdot10^{-62}\). \(\square\)

**What stays open.**
- *From BU/BV.* A closed form of \(F=\theta\), and so of \(c_\infty\), \(C_N\), \(P_N\), \(\chi_\infty\), \(G^*\), \(t_{\max}\), \(t_{\min}\). Whether \(\hat\theta_k\ne0\) for every \(k\ge12\), which is not implied by BV.3. Conjecture BV.6, the lower half of the decay rate. Whether the singularities of \(Q_\infty\) on \((-\infty,0)\) are poles (an atomic spectral measure would give BV.6). BU.3/BV.5's \(G^*\) was not fed back into BI.4/BJ.
- *From BS.* Conjecture BS.4 (\(\chi_N\) decreasing for every \(N\)); a rate for \(J-A_N\); uniformity in \(\alpha\).
- *Carried over.* Eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BV proves that one log-periodic function built from a corner resolvent of a graph Laplacian on the combinatorial gasket is analytic in a strip. It bounds and certifies that function's Fourier modes, and it locates its extremes. These are statements about functions on graphs, not a field, a flux or a force. No flux functional is constructed, and Stage 1 is not touched. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_analytic.py`** (about 46 s; needs sympy, python-flint, mpmath; imports AU's script for the certified mode routine).
- *Sympy:* the sector identity; \(|5^{(1-\gamma)t}|\); the contour factor; convexity and stationary point of \(-Tf-\log\cos f\); \(\arctan T+\arctan(1/T)=\pi/2\) and \(x\ge\arctan x\); the telescoping tail sum.
- *Exact Gaussian rationals:* \(|r_n(z)|^2\cos^2(\varphi/2)\le r_n(|z|)^2\) for \(r_0,\dots,r_3\) at 125 points, with rational \(\cos\varphi=(1-m^2)/(1+m^2)\) (BV.1 on the actual resolvents; maximum ratio \(0.9983\)).
- *Interval:* AU's routine for \(\langle F\rangle\) and \(a_k,b_k\) (\(k\le11\)); nonvanishing; agreement with AU's table; \(|\hat\theta_k|\le B_k\); the tails; the \(4096\)-cell critical-point count; bisection; Kept Failure BV.5's witness; consistency with BU.3/BU.4.
- *Floats:* only the printed ratios \(|\hat\theta_{k+1}/\hat\theta_k|\) and products \(|\hat\theta_k|e^{2\pi^2k/\ln5}\) (Conjecture BV.6's evidence). They enter no check.

The script exits nonzero if any check fails.

## Theorem BW. AU's mode routine, re-run deeper and with a finer sandwich, certifies \(\hat\theta_k\ne0\) for every \(k\le20\) (BV: \(k\le11\)), \(\langle F\rangle\) and \(c_R(\gamma)\) to \(10^{-113}\), and \(G^*\), \(\min F\) to \(10^{-109}\)

**Setting.** As in Theorems BS–BV: \(F(u)=u^{1-\gamma}Q_\infty(u)=\theta(u)\), \(F(5u)=F(u)\), \(t=\log_5u\), \(\hat\theta_k=\int_0^1F(5^t)e^{-2\pi ikt}\,dt=a_k-ib_k\) (AU's normalisation), \(G^*=\max F\), \(B_k=G^*\sqrt{1+T_k^2}\,e^{-T_k\arctan T_k}\), \(T_k=4\pi k/\ln5\) (BV.3). **Assumptions (used as stated, not re-run):** AQ (3), BN.1, BS.1, BT.2, BU.1, BU.3 (\(G^*\le0.6009615785983\)), BV.3 (the bound \(|\hat\theta_k|\le B_k\)); AU (3)'s certified routine (AQ (4) node enclosures, AS (2) sandwich for any floor/ceiling orders \(b\), \(u\), AT.1 orbit lemma, AU's \(\psi\) and endpoint corrections), imported unchanged. AU's 22-digit table, BV.4's table and BV.5's balls are used for consistency only.

**THEOREM BW.**

1. **BW.1 (modes \(k\le20\); BV's open "\(\hat\theta_k\ne0\) for \(k\ge12\)", up to \(k=20\)).** AU (3)'s routine with period \([5^{-261},5^{-260}]\), \(1024\) cells per period, sandwich orders \(b=u=32/33\), outward-rounded `mpmath.iv` at \(300\) digits and node cap \(10^{240}\) encloses \(\langle F\rangle\) to radius \(1.0\cdot10^{-113}\) and every \(a_k,b_k\), \(k\le20\), to radius \(2.0\cdot10^{-113}\) (BV.4: \(8.2\cdot10^{-63}\)). Every one of them has \(a_k\ne0\) or \(b_k\ne0\) certified:

| \(k\) | \(a_k\) | \(b_k\) | \(\lvert\hat\theta_k\rvert\ge\) | \(B_k\) | \(\lvert\hat\theta_k\rvert e^{2\pi^2k/\ln5}\) |
|---|---|---|---|---|---|
| 12 | \(-7.92030433941\cdot10^{-65}\) | \(1.75678609788\cdot10^{-65}\) | \(8.112800\cdot10^{-65}\) | \(1.8500\cdot10^{-62}\) | \(0.6712\) |
| 13 | \(-7.23926469953\cdot10^{-71}\) | \(5.09345693186\cdot10^{-70}\) | \(5.144645\cdot10^{-70}\) | \(9.4507\cdot10^{-68}\) | \(0.9027\) |
| 14 | \(1.77240534632\cdot10^{-75}\) | \(1.48737003819\cdot10^{-75}\) | \(2.313804\cdot10^{-75}\) | \(4.7992\cdot10^{-73}\) | \(0.8609\) |
| 15 | \(6.60019679401\cdot10^{-81}\) | \(-5.01730355973\cdot10^{-81}\) | \(8.290714\cdot10^{-81}\) | \(2.4247\cdot10^{-78}\) | \(0.6542\) |
| 16 | \(-2.48373767130\cdot10^{-86}\) | \(-4.00131229987\cdot10^{-86}\) | \(4.709507\cdot10^{-86}\) | \(1.2196\cdot10^{-83}\) | \(0.7881\) |
| 17 | \(-2.44314458057\cdot10^{-91}\) | \(1.18102931611\cdot10^{-92}\) | \(2.445997\cdot10^{-91}\) | \(6.1104\cdot10^{-89}\) | \(0.8680\) |
| 18 | \(-2.41716211798\cdot10^{-97}\) | \(8.84447012021\cdot10^{-97}\) | \(9.168824\cdot10^{-97}\) | \(3.0508\cdot10^{-94}\) | \(0.6900\) |
| 19 | \(4.64767606746\cdot10^{-102}\) | \(1.71967465130\cdot10^{-102}\) | \(4.955620\cdot10^{-102}\) | \(1.5185\cdot10^{-99}\) | \(0.7909\) |
| 20 | \(2.04957487095\cdot10^{-107}\) | \(-1.75889774588\cdot10^{-107}\) | \(2.700826\cdot10^{-107}\) | \(7.5375\cdot10^{-105}\) | \(0.9141\) |

   So \(\hat\theta_k\ne0\) for every \(1\le|k|\le20\): AU (2)'s Mellin transform \(M\) of \(r_\infty\) has genuine simple poles at \(1-\gamma\pm i\omega k\) for all these \(k\), and \(F\) is not a trigonometric polynomial of degree below \(20\). Also \(|\hat\theta_k|\le B_k\) for every \(k\le20\) (a check of BV.3 against an independent certificate), and \(|\hat\theta_k|e^{2\pi^2k/\ln5}\in[0.60,0.93]\) for every \(k\le20\) (certified; nine more terms of Conjecture BV.6's evidence band, with no trend). The rows \(k\le4\) reproduce AU's 22-digit table, and the rows \(k=5,\dots,11\) agree with BV.4's table within BV.4's radius \(8.2\cdot10^{-63}\). (BV.4's printed \(b_{11}=-1.78263934544\cdot10^{-59}\) has a non-significant last digit; the new value is \(-1.78263934547\cdot10^{-59}\), and the new \(|\hat\theta_{11}|\ge1.808283\cdot10^{-59}\) sharpens BV.4's \(1.807336\cdot10^{-59}\).)
2. **BW.2 (the mean and the slope constant).**
   \[
   \langle F\rangle=\langle\Theta\rangle=0.60095425116877360782370515009999426797183899606527964527556307290326859178311545042930649582114\ldots\ (\pm1.0\cdot10^{-113}),
   \]
   and AU (1)'s \(c_R(\gamma)=(2\sin\pi\gamma/\pi)\ln5\,\langle\Theta\rangle\) is
   \[
   c_R(\gamma)=0.51716583463050524392598048567312720045378382987715870184537398574670709450856254620631814819510\ldots\ (\pm8.7\cdot10^{-114}),
   \]
   both inside AU's \(34\)-digit enclosures.
3. **BW.3 (extremes; BV.5's argument with the \(20\)-mode sum).** BV.3's tail from \(k\ge21\) gives \(|F-S|\le7.46\cdot10^{-110}\), \(|F'-S'|\le9.85\cdot10^{-108}\), \(|F''-S''|\le1.30\cdot10^{-105}\). The \(4096\)-cell count again gives exactly one maximum and one minimum per period, and
   \[
   G^*=\max F=0.60096157859710108904205618413461692677004052591725819702185687716025065514202295750071720787122\ldots\ (\pm1.3\cdot10^{-109}),
   \]
   \[
   \min F=0.60094692379426561210422618817212924841162969178767769700850974102950276770653887623071791572010\ldots\ (\pm1.3\cdot10^{-109}),
   \]
   \[
   t_{\max}=0.751755961790317558778029332738208256562201926116552145785615717872\ldots,\quad t_{\min}=0.251756124192213928849302877000673101065745614824539815823039223608\ldots\ (\pm5.6\cdot10^{-104}),
   \]
   all agreeing with BV.5's printed digits. Kept Failure BV.5's witness is reproduced: \(t_{\max}-t_{\min}-\tfrac12=-1.624018963700712735442624648445035436887\ldots\cdot10^{-7}\). So BI.4's sharp constant \(\sup_{s>0}s^{1-\gamma}r_\infty(s)=G^*\) (BT.3) is now known to \(108\) digits.

**Proof.** BW.1 is AU (3)'s proof with new parameters: the script sets `gasket_fractional_log_mean`'s `E_DEPTH=260`, `K_PER=1024`, `ORD=32`, cap \(10^{240}\) and calls BV's `au_modes` (which calls AU's `period_integrals` and `corrections` unchanged), and converts each interval to python-flint `arb` exactly. AS (2) allows any floor and ceiling orders, and every correction term of AU (3) is evaluated at the new depth by the same code. A nonzero lower bound on \(|a_k|\) or \(|b_k|\) gives \(\hat\theta_k\ne0\); AU (2) identifies \(\overline{\hat\theta_k}\) with the residue of \(M\) at \(1-\gamma+i\omega k\). A trigonometric polynomial of degree \(d\) has \(\hat\theta_k=0\) for \(|k|>d\), so \(\hat\theta_{20}\ne0\) excludes \(d<20\). BW.2 is AU (1)'s identity evaluated in `arb`. BW.3 is BV.5's proof verbatim, with \(S\) the certified \(20\)-mode sum and the tail started at \(k=21\). \(\square\)

**Observation (floats of certified numbers; no claim).** The certified radius dropped from \(8.2\cdot10^{-63}\) to \(3.7\cdot10^{-77}\) (depth \(5^{-200}\), \(512\) cells, \(b=24\)), \(1.3\cdot10^{-85}\) (\(768\) cells, \(b=24\)) and \(8.9\cdot10^{-92}\) (\(512\) cells, \(b=32\)) in exploratory runs, so at these settings the sandwich order, not the depth, was the limiting term; at \(b=32\), \(5^{-200}\) the radius is near AU's AT.1 depth correction \(\approx2(5^{-201})^{0.65}/(1-5^{-0.65})\). These runs are not part of the witness.

**What stays open.**
- *From BV/BW.* Whether \(\hat\theta_k\ne0\) for every \(k\ge21\) (not implied by BV.3 or BW.1; each further \(k\) costs about \(5.3\) more digits). Conjecture BV.6 (lower half of the decay rate): the certified band \([0.60,0.93]\) for \(k\le20\) is evidence, not proof, and the normalised moduli oscillate without a visible limit. Whether the singularities of \(Q_\infty\) on \((-\infty,0)\) are poles. A closed form of \(F=\theta\), \(G^*\), \(t_{\max}\), \(t_{\min}\). Feeding \(G^*\) back into BI.4/BJ.
- *Carried over.* Conjecture BS.4; a rate for \(J-A_N\); uniformity in \(\alpha\); eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BW certifies more Fourier modes and more digits of one log-periodic function built from a corner resolvent of a graph Laplacian on the combinatorial gasket. It proves no new structural property beyond those listed. No flux functional is constructed, Stage 1 is not touched, and no \(W\), \(0.08\), \(0.23\) or \(30\,\mu\)N/kW enters. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_modes.py`** (about 6 min; needs python-flint, mpmath; imports BV's and AU's scripts).
- *Interval:* AU's routine at the new parameters for \(\langle F\rangle\) and \(a_k,b_k\) (\(k\le20\)); nonvanishing; \(|\hat\theta_k|\le B_k\); the band \([0.60,0.93]\); agreement with AU's and BV.4's tables; \(c_R(\gamma)\) inside AU's enclosure; the tails; the \(4096\)-cell critical-point count; bisection; agreement with BV.5's printed digits; BV.5's witness sign.
- *Floats:* none enter a check.

The script exits nonzero if any check fails.

## Theorem BX. The profile \(F=\theta\) is a convolution of one explicit kernel with a positive periodic spectral measure \(m\) of mass \(c_R(\gamma)/2\), so \(\hat\theta_k=\frac{\pi}{\ln5}\,\hat m_k/\sin\pi(\gamma+i\omega k)\) exactly; this gives the sharp bound \(|\hat\theta_k|e^{2\pi^2|k|/\ln5}\le\pi c_R(\gamma)/\ln5\,(1+o(1))\) (BV.3's factor \(1+4\pi k/\ln5\) removed) and turns Conjecture BV.6 into a statement about \(\hat m_k\) (Kept Failure BX.5: \(m\) has at least \(17\) support points per period)

**Setting.** As in Theorems BS–BW: \(r_\infty=\lim r_k\), \(Q_\infty(u)=\lim_M(3/5)^Mr_\infty(5^{-M}u)\), \(\gamma=\log3/\log5\), \(\omega=2\pi/\ln5\), \(F(u)=u^{1-\gamma}Q_\infty(u)=\theta(u)\), \(t=\log_5u\), \(\hat\theta_k=\int_0^1F(5^t)e^{-2\pi ikt}\,dt=a_k-ib_k\) (AU's normalisation). \(\mathbb T=\mathbb R/\mathbb Z\). **Assumptions (used as stated, not re-run):** AQ (3); BN.1 (\(r_k\) Stieltjes); BS.1 (\(Q_\infty\) exists, \(Q_\infty(5u)=\tfrac35Q_\infty(u)\), \(Q_\infty>0\)); BT.2, BU.1 (\(F=\theta\)); BV.2 (the Stieltjes functions \(f_M(z)=(3/5)^Mr_\infty(5^{-M}z)\) converge to \(Q_\infty\) locally uniformly on \(\mathbb C\setminus(-\infty,0]\), by Vitali); AU (1) (\(c_R(\gamma)=\frac{2\sin\pi\gamma}{\pi}\ln5\,\langle\Theta\rangle\)); BW's certified routine (AU (3) unchanged, BW's parameters), imported unchanged. **ASSUMPTION BX (literature, not proved here).** The Nevanlinna characterisation of Stieltjes functions: \(f\ge0\) on \((0,\infty)\), holomorphic on \(\mathbb C\setminus(-\infty,0]\) with \(\operatorname{Im}z\cdot\operatorname{Im}f(z)\le0\), if and only if \(f(u)=a/u+b+\int_{(0,\infty)}d\nu(\lambda)/(\lambda+u)\) with \(a,b\ge0\), \(\nu\ge0\), and this representation is unique (e.g. Schilling–Song–Vondraček, *Bernstein Functions*, ch. 7). AU's, BV.4's tables and AU's \(c_R(\gamma)\) enclosure are used for consistency only. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BX.**

1. **BX.1 (Stieltjes representation of \(Q_\infty\)).** \(Q_\infty(u)=\int_{(0,\infty)}d\nu(\lambda)/(\lambda+u)\) with \(\nu\ge0\), \(\nu\ne0\), no atom at \(0\), no constant term, and \(\nu(5A)=3\nu(A)\) for every Borel \(A\). In particular the singularities of \(Q_\infty\) on \((-\infty,0)\) are exactly \(-\operatorname{supp}\nu\), and \(Q_\infty\) has a pole at \(-\lambda_0\) if and only if \(\lambda_0\) is an isolated atom of \(\nu\) (BW's open question, restated).
2. **BX.2 (factorisation).** In \(t=\log_5\lambda\), \(\sigma=\lambda^{-\gamma}\nu\) is \(1\)-periodic. Let \(m\) be its restriction to one period, a positive measure on \(\mathbb T\), with \(\hat m_k=\int_{\mathbb T}e^{-2\pi ikt}\,dm(t)\). Then, with \(K(x)=5^{\gamma x}/(5^x+1)\),
   \[
   F(\tau)=\int_{\mathbb R}K(t-\tau)\,d\sigma(t),\qquad \hat\theta_k=\frac{\pi}{\ln5}\,\frac{\hat m_k}{\sin\pi(\gamma+i\omega k)}\quad(k\in\mathbb Z).
   \]
   At \(k=0\): \(m(\mathbb T)=\langle F\rangle\ln5\,\sin(\pi\gamma)/\pi=c_R(\gamma)/2\). So AR/AU's slope constant at \(\gamma\) is twice the spectral mass of \(\nu\) per scale period. Certified: \(m(\mathbb T)=0.258582917315252621962990242836563600226891914938579350922687\ldots\ (\pm4.3\cdot10^{-114})\).
3. **BX.3 (sharp mode bound; Conjecture BV.6 reformulated).** \(|\hat m_k|\le\hat m_0\) and \(|\sin\pi(\gamma+iy)|^2=\sin^2\pi\gamma+\sinh^2\pi y\) give, for every \(k\ne0\),
   \[
   |\hat\theta_k|\le\frac{\langle F\rangle\sin\pi\gamma}{\sqrt{\sin^2\pi\gamma+\sinh^2(2\pi^2k/\ln5)}}\le\frac{\pi c_R(\gamma)}{\ln5}\,\frac{e^{-2\pi^2|k|/\ln5}}{1-e^{-4\pi^2|k|/\ln5}},
   \]
   so \(\limsup_k|\hat\theta_k|e^{2\pi^2|k|/\ln5}\le\pi c_R(\gamma)/\ln5=1.00949802052668788\ldots\). BV.3's \(B_k\le eG^*(1+4\pi k/\ln5)e^{-2\pi^2k/\ln5}\) carries a factor that grows linearly in \(k\). The new bound is smaller than \(B_k\) by factors \(12.7\) (\(k=1\)) to \(252.7\) (\(k=20\)). The bound is attained for a given \(k\) only if \(m\) sits on one coset of \(\tfrac1k\mathbb Z\). Exactly, \(|\hat\theta_k|e^{2\pi^2|k|/\ln5}=\frac{\pi c_R(\gamma)}{\ln5}\,\rho_k\,(1+O(e^{-4\pi^2|k|/\ln5}))\) with \(\rho_k=|\hat m_k|/\hat m_0\in[0,1]\). So **Conjecture BV.6 holds if and only if \(\rho_k\) converges to a nonzero limit**, and \(\hat\theta_k\ne0\) if and only if \(\hat m_k\ne0\).
4. **BX.4 (certified, BW's routine re-run unchanged).** With BW's enclosures of \(\langle F\rangle\), \(a_k\), \(b_k\) (\(k\le20\)) the moduli \(\rho_k\) are enclosed to radius \(\le1.1\cdot10^{-6}\) (\(\le6.2\cdot10^{-12}\) for \(k\le19\)):

| \(k\) | \(\rho_k\) | \(k\) | \(\rho_k\) | \(k\) | \(\rho_k\) | \(k\) | \(\rho_k\) |
|---|---|---|---|---|---|---|---|
| 1 | 0.769641861455 | 6 | 0.776865099532 | 11 | 0.698838195911 | 16 | 0.780650379314 |
| 2 | 0.600852173802 | 7 | 0.825750190597 | 12 | 0.664898347061 | 17 | 0.859827502587 |
| 3 | 0.833706407581 | 8 | 0.667820624119 | 13 | 0.894157628226 | 18 | 0.683507212616 |
| 4 | 0.855147103961 | 9 | 0.791030554602 | 14 | 0.852824081268 | 19 | 0.78343268637 |
| 5 | 0.721589627521 | 10 | 0.918121116259 | 15 | 0.648035668832 | 20 | 0.90547 |

   So \(\rho_k\in[0.60,0.92]\) for every \(k\le20\). This is BW's band divided by \(\pi c_R(\gamma)/\ln5\), still with no visible limit. \(|\hat\theta_k|\) is below the BX.3 bound for every \(k\le20\), which checks BX.3 against an independent certificate. The \(17\times17\) Hermitian Toeplitz section \(T_{17}=(\hat m_{i-j})_{0\le i,j<17}\) is certified positive definite: the exact-midpoint \(LDL^*\) of \(T_{17}^{\rm mid}-2\varepsilon I\) in `arb` (1200 bits) has every pivot \(>0\) (smallest \(5.96\cdot10^{-14}\)), where \(\varepsilon\ge\|T_{17}-T_{17}^{\rm mid}\|_F\) is \(3.4\cdot10^{-28}\). This also checks BX.2's positivity of \(m\) against independent numbers: a non-positive \(m\) need not give a positive semidefinite \(T\).

**Proof.** *BX.1.* Each \(r_k=\sum_iw_i/(\lambda_i+s)\) satisfies \(\operatorname{Im}z\cdot\operatorname{Im}r_k(z)\le0\) (BN.1). BV.2's locally uniform limit \(r_\infty\) keeps this and is \(\ge0\) on \((0,\infty)\), so \(r_\infty\) is Stieltjes by Assumption BX, and so is each \(f_M\). The same argument applied to \(f_M\to Q_\infty\) (BV.2) gives \(Q_\infty=a/u+b+\int d\nu/(\lambda+u)\). Write \(Q_\infty(5u)\) in the same form, with \(\nu\) pushed forward by \(\lambda\mapsto\lambda/5\) and divided by \(5\), and compare with \(\tfrac35Q_\infty(u)\) by uniqueness. This gives \(a/5=3a/5\), \(b=3b/5\) and \(\nu(5A)/5=\tfrac35\nu(A)\) (sympy), so \(a=b=0\) and \(\nu(5A)=3\nu(A)\). \(Q_\infty>0\) gives \(\nu\ne0\). A Stieltjes transform continues across \((-\infty,0)\) exactly off \(-\operatorname{supp}\nu\), with a pole exactly at isolated atoms (standard, by Stieltjes inversion; part of Assumption BX). *BX.2.* \(\sigma(A+1)=\int_{5A}(5\mu)^{-\gamma}d\nu(5\mu)=3\cdot5^{-\gamma}\sigma(A)=\sigma(A)\), since \(5^\gamma=3\) (sympy). With \(\lambda=5^t\), \(5^{(1-\gamma)\tau}\lambda^{\gamma}/(\lambda+5^\tau)=K(t-\tau)\) (sympy), which is the convolution formula. All integrands are positive or dominated by the \(k=0\) case, which is finite, so Fubini over \(\tau\in[0,1)\), \(t\in\mathbb R=\bigcup_j([0,1)+j)\) gives \(\hat\theta_k=\hat m_k\int_{\mathbb R}K(x)e^{2\pi ikx}dx\). Substituting \(v=5^x\) turns this into \(\frac1{\ln5}\int_0^\infty v^{a-1}/(1+v)\,dv=\frac{\pi}{\ln5\sin\pi a}\) with \(a=\gamma+i\omega k\), \(0<\operatorname{Re}a<1\). (Sympy checks a real \(a\); both sides are analytic in \(a\). A float quadrature of the single-series toy \(\nu=\sum_n3^n\delta_{5^n}\), i.e. \(m=\delta_0\), reproduces the formula at \(k=0,1,2\); that only checks the normalisation.) At \(k=0\), AU (1) gives \(c_R(\gamma)/2\). *BX.3.* \(|\hat m_k|\le m(\mathbb T)\) holds for a positive measure. The bound follows from \(|\sin\pi(\gamma+iy)|\ge\sinh\pi y\) and \(e^y/(2\sinh y)=1/(1-e^{-2y})\) (sympy). Equality \(|\hat m_k|=\hat m_0\) forces \(e^{-2\pi ikt}\) to be constant \(m\)-a.e. *BX.4.* BW's routine (AU (3) with BW's parameters, about \(6.3\) min) gives `arb` enclosures. \(\hat m_k\) is computed from BX.2 in `acb`. For the Toeplitz step, \(T=T^{\rm mid}+E\) with \(E\) Hermitian and \(\|E\|_2\le\|E\|_F\le\varepsilon\), so \(T^{\rm mid}\succ2\varepsilon I\) gives \(T\succ\varepsilon I\succ0\). The script also checks the certificate on an exact \(3\)-atom measure: \(T_3\) is certified and \(T_4\) is not. \(\square\)

**KEPT FAILURE BX.5 (certified witness).** Claim: "\(\nu\) is carried by at most \(r\le16\) geometric series \(\lambda_j5^{\mathbb Z}\)". Equivalently, \(Q_\infty(u)=\sum_{j\le r}c_j\sum_{n\in\mathbb Z}3^n/(\lambda_j5^n+u)\), a closed form made of \(r\) self-similar pole series. For \(r=1\) this is the most natural closed form of \(F\): \(F(\tau)=c\sum_nK(t_0+n-\tau)\). **The claim is false.** Such a \(\nu\) makes \(m\) a sum of \(r\) atoms. Then \(T_{r+1}=\sum_jw_jv_jv_j^*\) has rank \(\le r\) (sympy: \(\det T_3=0\) for two symbolic atoms), but \(T_{17}\succ0\) (BX.4). **Witness:** the smallest \(LDL^*\) pivot of \(T_{17}-2\varepsilon I\) is \(\ge5.9\cdot10^{-14}>0\). For \(r=1\) there is a one-number witness: \(\rho_1\le0.7696418615<1\), while one atom forces \(\rho_k=1\) for every \(k\). So the spectral measure behind \(F\) has at least \(17\) points in every scale period \([\lambda,5\lambda)\). A finite superposition of at most \(16\) self-similar pole series is not a closed form for \(F\), \(G^*\), \(t_{\max}\) or \(t_{\min}\).

**Observation (floats of certified numbers; no claim).** \(T_{18}\) is not decided at this precision: \(\varepsilon\approx7\cdot10^{-23}\) is comparable to its float smallest eigenvalue \(\approx5\cdot10^{-23}\). The smallest eigenvalue of \(T_n\) falls by a factor of about \(20\)–\(50\) per step, as for a measure with many atoms or a singular part. \(\arg\hat m_5/2\pi,\arg\hat m_{10}/2\pi,\arg\hat m_{15}/2\pi,\arg\hat m_{20}/2\pi=-0.0018,+0.0060,+0.0121,+0.0216\), and \(\rho_{10},\rho_{20}>0.9\). This hints at partial concentration of \(m\) near a coset of \(\tfrac15\mathbb Z\), i.e. near \(\lambda\in\lambda_05^{\mathbb Z/5}\). Not a claim.

**What stays open.**
- *From BV/BW/BX.* Conjecture BV.6, now "\(\rho_k\to\rho_\infty>0\)". Remark, not used: by Wiener's lemma an atomic \(m\) has \(\limsup\rho_k>0\), while for several atoms in general position \(\rho_k\) oscillates, so BV.6 as worded ("asymptotic to a constant") may be too strong. Whether \(\nu\) (equivalently \(m\)) is purely atomic, i.e. whether the singularities of \(Q_\infty\) are poles: the infinite-gasket spectral theory suggests pure point, but that is not used or proved here. Whether \(\hat\theta_k\ne0\) for every \(k\ge21\), equivalently \(\hat m_k\ne0\). An identification of \(\nu\) with the corner spectral measure of the blown-up Laplacian. A closed form of \(F\), \(G^*\), \(t_{\max}\), \(t_{\min}\). Feeding \(G^*\) back into BI.4/BJ.
- *Carried over.* Conjecture BS.4; a rate for \(J-A_N\); uniformity in \(\alpha\); eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BX is a structural identity and a finite certificate for one log-periodic function built from a corner resolvent of a graph Laplacian on the combinatorial gasket. The "spectral measure" \(m\) is a positive measure on the scale circle, not a flux, and its mass \(c_R(\gamma)/2\) is not a force. No flux functional is constructed, Stage 1 is not touched, and no \(W\), \(0.08\), \(0.23\) or \(30\,\mu\)N/kW enters. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_spectral.py`** (about 6.5 min; needs python-flint, mpmath, sympy; imports BW's, BV's and AU's scripts).
- *Sympy:* uniqueness bookkeeping \(a=b=0\); \(3\cdot5^{-\gamma}=1\); the kernel identity; \(\int_0^\infty v^{a-1}/(1+v)=\pi/\sin\pi a\); \(|\sin(x+iy)|^2=\sin^2x+\sinh^2y\); \(e^y/(2\sinh y)=1/(1-e^{-2y})\); \(\det T_3=0\) for two atoms.
- *Interval:* BW's routine; \(m(\mathbb T)=c_R(\gamma)/2\) inside AU's enclosure; \(\rho_k\) for \(k\le20\) and the band \([0.60,0.92]\); \(|\hat\theta_k|\le\) BX.3 bound \(<B_k\); \(\pi c_R(\gamma)/\ln5\in[1.0094,1.0095]\); \(\rho_1<0.7697\); \(T_{17}\succ0\); the \(3\)-atom sanity check; AU's and BV.4's published tables against BX.3 and the band.
- *Floats:* the single-series normalisation check only.

The script exits nonzero if any check fails.

## Theorem BY. A deeper sandwich certifies \(\hat\theta_{21}\ne0\) and the Toeplitz sections \(T_{18},T_{19},T_{20}\) of BX's spectral measure positive definite (Kept Failure BY.3: \(\nu\) is not carried by at most \(19\) geometric pole series)

**Setting.** As in Theorem BX: \(\gamma=\log3/\log5\), \(\omega=2\pi/\ln5\), \(F=\theta\) the \(1\)-periodic profile in \(t=\log_5u\), \(\hat\theta_k=a_k-ib_k\), and \(m\ge0\) the periodic spectral measure of BX.2, with \(\hat\theta_k=\frac{\pi}{\ln5}\hat m_k/\sin\pi(\gamma+i\omega k)\) and \(\rho_k=|\hat m_k|/\hat m_0\). BX left two finite questions open at its precision: whether \(\hat\theta_{21}\ne0\), and whether \(T_{18}=(\hat m_{i-j})_{i,j<18}\succ0\). The obstruction was the radius, not the method. \(\hat m_k\) is \(\hat\theta_k\) times about \(e^{2\pi^2k/\ln5}\), so BW's radius \(2\cdot10^{-113}\) on \(\hat\theta_{17}\) becomes \(3.5\cdot10^{-23}\) on \(\hat m_{17}\), and that alone moves the \(18\)th \(LDL^*\) pivot by more than its size. BW observed that the sandwich order limits the radius.

**Assumptions used as stated (not re-run).** AQ (3), BN.1, BS.1, BT.2, BU.1, BV.2, AU (1); Assumption BX (Nevanlinna characterisation of Stieltjes functions); BX.2, BX.3 and BX.5's rank lemma; AU (3)'s certified routine (AQ (4) node enclosures, AS (2) sandwich for any \(b,u\), AT.1 orbit lemma, AU's \(\psi\) and endpoint corrections), imported unchanged with new parameters.

**THEOREM BY.**
1. **BY.1 (modes).** AU (3)'s routine with period \([5^{-341},5^{-340}]\), \(1024\) cells, sandwich \(b=u=40/41\), `mpmath.iv` outward at \(380\) digits and node cap \(10^{310}\) encloses \(\langle F\rangle\) to \(\pm9.3\cdot10^{-134}\) and every \(a_k,b_k\) with \(k\le21\) to \(\pm1.9\cdot10^{-133}\) (BW: \(\pm2.0\cdot10^{-113}\)). Every \(a_k\) or \(b_k\) is certified nonzero, so \(\hat\theta_k\ne0\) for every \(1\le|k|\le21\). For example, \(a_{21}=-2.55838576504077\cdot10^{-113}\), \(b_{21}=-9.60025191043711\cdot10^{-113}\), and \(|\hat\theta_{21}|\ge9.93\cdot10^{-113}\). AU's (\(k\le4\)) and BV.4's (\(k=5,\dots,11\)) published tables are reproduced, and \(\langle F\rangle\) lies inside AU's enclosure.
2. **BY.2 (normalised moduli).** With BX.2, \(\rho_k\in[0.60,0.92]\) for every \(k\le21\), with \(\rho_{21}=0.706371132\ldots\). Also \(|\hat\theta_k|\le\langle F\rangle\sin\pi\gamma/|\sin\pi(\gamma+i\omega k)|\) (BX.3) for every \(k\le21\).
3. **BY.3 (Toeplitz).** For \(n=18,19,20\), the exact-midpoint \(LDL^*\) of \(T_n-2\varepsilon_nI\) in `arb` has every pivot certified positive, with \(\varepsilon_n\ge\|E\|_F\) the Frobenius norm of the enclosure radii. By Weyl, \(T_n\succ\varepsilon_nI\succ0\). The smallest pivots are \(5.33\cdot10^{-15}\), \(4.51\cdot10^{-16}\) and \(1.21\cdot10^{-16}\), and \(\varepsilon_{20}\le3.0\cdot10^{-32}\).

*Proof.* BY.1 is AU (3)'s routine, unchanged, run with new parameters. A nonzero real or imaginary part certifies \(\hat\theta_k\ne0\). BY.2 is BX.2's exact factorisation, applied in `acb` to the BY.1 enclosures, together with BX.3's bound. For BY.3, \(T_n=T_n^{\rm mid}+E\) with \(E\) Hermitian and \(\|E\|_2\le\|E\|_F\le\varepsilon_n\). So \(T_n^{\rm mid}\succ2\varepsilon_nI\), which the positive pivots certify, gives \(T_n\succ0\). The same certificate passes on an exact \(3\)-atom measure for \(T_3\) and fails for \(T_4\). \(\square\)

**KEPT FAILURE BY.3 (certified witness).** Claim: "\(\nu\) is carried by at most \(r\le19\) geometric series \(\lambda_j5^{\mathbb Z}\)", that is, \(Q_\infty(u)=\sum_{j\le r}c_j\sum_{n\in\mathbb Z}3^n/(\lambda_j5^n+u)\). **The claim is false.** Such a \(\nu\) makes \(m\) a sum of \(r\) atoms, so \(\operatorname{rank}T_{r+1}\le r\) (BX.5's lemma). But \(T_{20}\succ0\). **Witness:** the smallest pivot of \(T_{20}-2\varepsilon_{20}I\) is \(\ge1.2\cdot10^{-16}>0\), with \(\varepsilon_{20}\le3.0\cdot10^{-32}\). So \(m\) has at least \(20\) points in every scale period, which sharpens BX.5's \(17\).

**Observation (no claim).** \(T_{21}\) is not certified at this precision: its \(21\)st midpoint pivot falls below \(2\varepsilon_{21}\approx1.3\cdot10^{-26}\). The smallest pivot drops by a factor of about \(4\)–\(12\) per step from \(n=17\) to \(20\). That is consistent with many atoms or a continuous part, and it decides neither.

**What stays open.**
- *From BX/BY.* Conjecture BV.6 (\(\rho_k\to\rho_\infty>0\)); the certified band \([0.60,0.92]\) now covers \(k\le21\), with no visible limit. Whether \(\nu\) is purely atomic, and whether it has finitely or infinitely many atoms per period. Whether \(\hat\theta_k\ne0\) for every \(k\ge22\); each further \(k\) costs about \(5.3\) digits, and the sandwich order sets the radius. Whether \(T_{21}\succ0\). A closed form of \(F\), \(G^*\), \(t_{\max}\), \(t_{\min}\). Feeding \(G^*\) back into BI.4/BJ.
- *Carried over.* Conjecture BS.4; a rate for \(J-A_N\); uniformity in \(\alpha\); eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem BY is a finite certificate for one log-periodic function built from a corner resolvent of a graph Laplacian on the combinatorial gasket. A lower bound on the number of atoms of a spectral measure on the scale circle is not a flux, and nothing here is a force. Stage 1 is not touched, and no \(W\), \(0.08\), \(0.23\) or \(30\,\mu\)N/kW enters. The Stage 2 question stays OPEN. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_toeplitz.py`** (about 18 min; needs python-flint, mpmath, sympy; imports BX's, BW's, BV's and AU's scripts).
- *Interval:* AU's routine at the new parameters; radius \(<10^{-132}\); \(\hat\theta_k\ne0\) for \(k\le21\); AU's and BV.4's tables; \(\langle F\rangle\) inside AU's enclosure; \(\rho_k\in[0.60,0.92]\) and the BX.3 bound for \(k\le21\); \(T_{18},T_{19},T_{20}\succ0\); BX's \(3\)-atom sanity check.
- *Floats:* none enter a check.

The script exits nonzero if any check fails.

## Theorem BZ. \(Q_\infty\) continues meromorphically across \((-\infty,0)\) off a Lebesgue-null set, so BX's measure \(\nu\) has no absolutely continuous part; \(311\) atoms of \(\nu\) (genuine poles of \(Q_\infty\)) are certified in one scale period and carry all but \(4.03\cdot10^{-5}\) of the mass of \(m\); hence \(\limsup\rho_k\ge0.99968\) and \(\liminf\rho_k\le0.59497\) (Kept Failure BZ.6: Conjecture BV.6's asymptotic part is false, its sharp-rate part holds; Kept Failure BZ.5: \(\nu\) is not absolutely continuous and not carried by at most \(310\) geometric pole series)

**Setting.** As in Theorems AQ, BD, BS and BX: \(p(z)=z(z+5)\), \(c(z)=(z+6)/((z+2)(z+5))\), \(e(z)=2z/(2+z)\), \(r_\infty=\lim r_n\), \(D=1/r_\infty\), \(\psi(v)=\lim_kp^k(5^{-k}v)\) (BD.4), \(Q_\infty(u)=\lim_M(3/5)^Mr_\infty(5^{-M}u)\) (BS.1), \(Q_\infty(u)=\int_{(0,\infty)}d\nu(\lambda)/(\lambda+u)\) with \(\nu(5A)=3\nu(A)\) (BX.1), and \(m\) the \(1\)-periodic measure \(\lambda^{-\gamma}\nu\) in \(t=\log_5\lambda\), with \(m(\mathbb T)=c_R(\gamma)/2\) (BX.2). Let \(K_p=\{z\in\mathbb C:(p^n(z))_n\text{ bounded}\}\) be the filled Julia set of \(p\), and \(Z=\psi^{-1}(K_p)\). BX and BY left open whether \(\nu\) is purely atomic, that is, whether the singularities of \(Q_\infty\) on \((-\infty,0)\) are poles. BY's Toeplitz route certified \(20\) support points per period, and its \(T_{21}\) stalled because the smallest eigenvalue of \(T_n\) falls by \(10\)–\(50\) per step. This note leaves the Toeplitz route and uses the decimation itself.

**Assumptions used as stated (not re-run).** AQ (1)–(3): the fixed point \(r_\infty(s)=r_\infty(ps)/(c-e\,r_\infty(ps))\); each \(r_n=\sum_ic_i/(\lambda_i+s)\) with \(\sum c_i=\tfrac12|v|^2=1\) and \(\lambda_i\in[0,8]\) (Gershgorin, maximal degree \(4\)); \(r_n\uparrow r_\infty\); \(r_\infty\) nonincreasing. BD.4: \(\psi_k\uparrow\psi\), \(\psi(5v)=p(\psi(v))\), \(\psi_k\) polynomials with nonnegative coefficients, \(\psi(v)\ge v\). BJ.6: \(\psi(v)\le ve^{v/10}\). BS.1 (existence of \(Q_\infty\) as that limit). BX.1–BX.2 (Stieltjes representation, \(m(\mathbb T)=c_R(\gamma)/2\)). BW's digits of \(G^*\), \(\min F\), \(t_{\max}\), \(t_{\min}\) and \(c_R(\gamma)\), for a consistency check and the mass comparison only; BY's printed \(\rho_1,\dots,\rho_5\), for a consistency check only; BX.2–BX.3's factorisation \(\hat\theta_k=\frac{\pi}{\ln5}\hat m_k/\sin\pi(\gamma+i\omega k)\). No literature assumption is used.

**THEOREM BZ.**
1. **BZ.1 (the Julia set of \(p\)).** \(p(z)+\tfrac52=(z+\tfrac52)^2-\tfrac{15}4\), so \(p\) is conjugate to \(f(\zeta)=\zeta^2-\tfrac{15}4\), with fixed point \(\tfrac52\). \(K_p\) is compact, \(p^{-1}(K_p)=K_p\), and \(K_p\subset[-5,0]\). It is covered by \(2^n\) intervals of total length \(\le5(2/\sqrt5)^n\), so it has Lebesgue measure \(0\). For \(z\notin K_p\), \(p^n\to\infty\) uniformly on a neighbourhood of \(z\).
2. **BZ.2 (\(r_\infty\) off the Julia set).** \(D(z)=c(z)\,D(p(z))-e(z)\). \(D\) is holomorphic on \(|w|>16\), with \(|D(w)-w|\le9\) for \(|w|\ge144\). Hence \(r_\infty\) is meromorphic on \(\mathbb C\setminus K_p\).
3. **BZ.3 (exact formula for \(Q_\infty\)).** For every \(u\in\mathbb C\setminus(-\infty,0]\),
   \[
   \frac1{Q_\infty(u)}=\Pi(u)\,D\bigl(\psi(u)\bigr)-S(u),\qquad \Pi(u)=\prod_{j\ge1}\tfrac53c\bigl(\psi(5^{-j}u)\bigr),\qquad S(u)=\sum_{k\ge1}\bigl(\tfrac53\bigr)^ke\bigl(\psi(5^{-k}u)\bigr)\,\Pi(5^{-k}u).
   \]
   \(\Pi\) and \(S\) are meromorphic on \(\mathbb C\). \(Z\) is a closed, Lebesgue-null subset of \((-\infty,0]\) with \(5Z=Z\), and \(Q_\infty\) continues meromorphically to \(\mathbb C\setminus Z\). The finite-level analogue \((5/3)^n/r_n(5^{-n}u)=\Pi_n(u)(3+\psi_n(u))-S_n(u)\) holds with AQ's \(r_n\), and the formula reproduces BW's \(G^*=F(t_{\max})\) and \(\min F=F(t_{\min})\) to \(10^{-38}\). With BS.1 (\(1/\mathcal R-1/Q_\infty=(1+2\eta(\psi))\Pi\)) and BH.2 (\(D=2+z-2\eta\)), as stated, it also gives \(1/\mathcal R=\Pi(3+\psi)-S\), which is meromorphic on all of \(\mathbb C\): the scaled finite-graph profile has no Julia-set singularities. This sentence is algebra on stated results and is not re-run.
4. **BZ.4 (no absolutely continuous part).** On \((0,\infty)\setminus(-Z)\), \(\nu\) is a sum of point masses with no accumulation point there; \(-Z\) is Lebesgue-null. So \(\nu\) is singular with respect to Lebesgue measure, and so is \(m\) on the circle.
5. **BZ.5 (certified atoms).** Scan \(t=\log_5\lambda\in[0,1]\) on \(128000\) cells and certify each sign change of \(1/Q_\infty(-5^t)\) (python-flint `arb`/`acb`, outward, 300 bits). This gives \(311\) distinct atoms \(\lambda_j\in(1,5)\) of \(\nu\), each with an enclosure of \(\nu(\{\lambda_j\})=1/(1/Q_\infty)'(-\lambda_j)\). The largest are \(\lambda=3.615256719887703\ldots\) with \(\nu(\{\lambda\})=0.479640692929\ldots\), \(\lambda=1.171910250998837\ldots\) with \(0.0359221043246\ldots\), and \(\lambda=1.016316417446789\ldots\) with \(0.0136278303937\ldots\). Summed over the \(311\) atoms, \(\lambda_j^{-\gamma}\nu(\{\lambda_j\})\in[0.25854262636241\pm10^{-15}]\). Against \(m(\mathbb T)=c_R(\gamma)/2=0.2585829173152526\ldots\), the part of \(m\) outside these atoms is \(\le4.03\cdot10^{-5}\), that is, \(0.016\%\).
6. **BZ.6 (the normalised moduli do not converge).** Put \(\delta=m(\mathbb T)-\sum_j\lambda_j^{-\gamma}\nu(\{\lambda_j\})\le4.03\cdot10^{-5}\). Then
   \[
   \limsup_{k\to\infty}\rho_k\ \ge\ 1-\frac{2\delta}{m(\mathbb T)}\ \ge\ 0.99968,\qquad
   \liminf_{k\to\infty}\rho_k\ \le\ \frac{|a^{\rm cert}_{35}|+\delta}{m(\mathbb T)}\ \le\ 0.59497,
   \]
   where \(a^{\rm cert}_k=\sum_j\lambda_j^{-\gamma}\nu(\{\lambda_j\})e^{-2\pi ik\log_5\lambda_j}\) runs over the \(311\) certified atoms (\(k=35\) minimises \(|a^{\rm cert}_k|\) over \(k\le100\)). So \(\rho_k\) has no limit. With BX.2, \(|\hat\theta_k|e^{2\pi^2k/\ln5}=\tfrac{\pi}{\ln5}|\hat m_k|\,e^{\pi\omega k}/|\sin\pi(\gamma+i\omega k)|\) and \(e^{\pi\omega k}/|\sin\pi(\gamma+i\omega k)|\to2\). Therefore \(\limsup_k|\hat\theta_k|e^{2\pi^2k/\ln5}\in[1.00918,\,1.0094980205]\): BX.3's sharp bound is attained to \(3.2\cdot10^{-4}\). The rate \(e^{-2\pi^2k/\ln5}\) is exact along a subsequence, so \(F\) continues holomorphically to no strip wider than \(|\operatorname{Im}t|<\pi/\ln5\).

*Proof.* *BZ.1.* Sympy gives the conjugacy. For \(|\zeta|>\tfrac52\), \(|f(\zeta)|-|\zeta|\ge(|\zeta|-\tfrac52)(|\zeta|+\tfrac32)\). The excess over \(\tfrac52\) therefore grows by a factor \(\ge4\) per step, \(\zeta\) escapes, and \(K_f\subset\Delta=\{|\zeta|\le\tfrac52\}\). On \(\Delta\), \(|w+\tfrac{15}4|\ge\tfrac54\), so both branches \(g_\pm(w)=\pm\sqrt{w+15/4}\) of \(f^{-1}\) are holomorphic on \(\Delta\), map it into \(\Delta\), map \([-\tfrac52,\tfrac52]\) into itself, and have \(|g_\pm'|\le1/\sqrt5\). Since \(\Delta\) is convex, \(K_f=\bigcap_n\bigcup_{s\in\{\pm\}^n}g_{s_1}\circ\cdots\circ g_{s_n}(\Delta)\) is a union of nested sets of diameter \(\le5\cdot5^{-n/2}\), each containing real points. So every point of \(K_f\) is real, \(K_f\subset[-\tfrac52,\tfrac52]\), and the \(2^n\) real images have total length \(\le5(2/\sqrt5)^n\to0\). Translate by \(-\tfrac52\). Complete invariance holds by definition. Uniform escape comes from the explicit growth outside \(\Delta\) and continuity.
*BZ.2.* AQ's fixed point is \(1/r=c/R-e\) with \(R=r(p)\) (sympy). Each \(r_n\) is bounded by \(1/\operatorname{dist}(w,[-8,0])\), so by Vitali \(r_\infty\) is holomorphic off \([-8,0]\). From \(wr_n(w)=\int d\mu_n/(1+\lambda/w)\), \(|wr_\infty-1|\le8/(|w|-8)<1\) for \(|w|>16\), so \(D\) is holomorphic there, and \(|D-w|\le8|w|/(|w|-16)\le9\) for \(|w|\ge144\). For \(z_0\notin K_p\), choose \(N\) with \(|p^N|>144\) near \(z_0\) (BZ.1). Unrolling \(D(z)=c\,D(p z)-e\) \(N\) times writes \(D\) near \(z_0\) as a rational expression in \(z\) and the holomorphic \(D(p^N z)\).
*BZ.3.* With \(x_j=\psi(5^{-j}u)\) we have \(p(x_j)=x_{j-1}\). Put \(X_k=(5/3)^kD(x_k)\). Then \(X_k=b_kX_{k-1}-(5/3)^ke(x_k)\) with \(b_k=\tfrac53c(x_k)\), so \(X_M=\prod_{j\le M}b_j\,X_0-\sum_{k\le M}(5/3)^ke(x_k)\prod_{k<j\le M}b_j\) (sympy, \(M=3\), generic). For \(u>0\), \(D\) is nondecreasing and \(5^{-M}u\le x_M\le5^{-M}ue^{5^{-M}u/10}\), so \(X_M\) is squeezed between \((5/3)^MD(5^{-M}u)\) and \((5/3)^MD(5^{-M}u')\) for any \(u'>u\) and large \(M\). Since \(Q_\infty\) is continuous, \(X_M\to1/Q_\infty(u)\) (BS.1). The tails converge locally uniformly on \(\mathbb C\): \(|x_j|\le\psi(5^{-j}|u|)\le2\cdot5^{-j}|u|\) (nonnegative coefficients), \(|\log b(x)|\le|x|\) and \(|e(x)|\le1.06|x|\) for \(|x|\le0.1\) (sympy constants), and \(\sum_k(5/3)^k5^{-k}<\infty\). So \(\Pi\) and \(S\) are meromorphic on \(\mathbb C\), and the formula holds on \((0,\infty)\). Both sides are meromorphic on the connected set \(\mathbb C\setminus(-\infty,0]\), so it holds there.
\(Z\subset(-\infty,0]\): if \(\psi(u)\in K_p\), then \(\psi(5^{-n}u)\in p^{-n}(K_p)=K_p\subset\mathbb R\) for every \(n\). But \(\psi(5^{-n}u)=5^{-n}u(1+O(5^{-n}))\) has nonzero imaginary part for large \(n\) unless \(u\in\mathbb R\), and \(\psi(u)\ge u>0\) for \(u>0\). Further, \(5Z=Z\) from \(\psi(5u)=p(\psi(u))\). \(Z\) is Lebesgue-null because \(\psi|_{\mathbb R}\) is real-analytic and nonconstant (isolated critical points, a local diffeomorphism elsewhere) and \(K_p\) is null. \(D\circ\psi\) is meromorphic on \(\mathbb C\setminus Z\) by BZ.2, so the right side is too, and its reciprocal continues \(Q_\infty\).
*BZ.4.* On an interval \(I\subset(-\infty,0)\setminus Z\), the continuation from either half-plane is the same meromorphic function. By Stieltjes inversion, \(\nu\) on \(-I\) is carried by the poles. These are simple, because \(|Q_\infty(x+iy)|=O(1/|y|)\), and the atom's mass is the residue.
*BZ.5.* \(f=1/Q_\infty\) is enclosed outward by BZ.3, as follows. \(x_{200}=v+B(\delta)\), \(v=5^{-200}u\), \(\delta=|v|(e^{|v|/10}-1)\) (BD.4, BJ.6). Then \(x_{j-1}=p(x_j)\). The tails are bounded by \(e^{\pm|u|5^{-200}}\) and \(1.1|u|3^{-200}\). \(D(x_0)\) is computed by iterating \(p\) to modulus \(>10^{60}\), using \(D(w)\in w+B(9)\) there, and recursing back, with every denominator certified nonzero. For each candidate cell, the script bisects to a sign change \(f(a)<0<f(b)\) with \(b-a\approx10^{-33}\). It then evaluates \(f\) on a complex box of half-side \(\rho\gg b-a\) around the midpoint \(c\). A finite enclosure means the box misses \(Z\) and every pole of \(f\), so \(f\) is holomorphic on the box and has a zero \(u^*\in(a,b)\). With \(M=\sup_{\rm box}|f|\), Cauchy estimates give \(f'(u^*)\in(f(c+h)-f(c-h))/(2h)\pm(h^2M/(\rho-h)^3+2M|u^*-c|/(\rho-2|u^*-c|)^2)\), \(h=10^{-30}\). This is certified \(>0\), so \(u^*\) is a simple pole of \(Q_\infty\) and \(\nu(\{-u^*\})=1/f'(u^*)\) (BZ.4). The enclosures of distinct \(\lambda_j\) are disjoint and lie in \((1,5)\).
*BZ.6.* Write \(m=m_{\rm pp}+m_c\) (atoms plus the rest). Then \(m_c(\mathbb T)+(m_{\rm pp}(\mathbb T)-\sum_j)=\delta\), and \(|\hat m_k-a_k|\le m_c(\mathbb T)\), where \(a_k\) are the coefficients of \(m_{\rm pp}\). The sequence \(a_k\) is Bohr almost periodic on \(\mathbb Z\): given \(\varepsilon\), keep finitely many atoms with the remaining mass \(<\varepsilon/4\). Dirichlet's simultaneous approximation gives arbitrarily large \(\tau\) with \(\|\tau t_j\|<\varepsilon/(4\pi m(\mathbb T))\) for those \(t_j\) (if the \(\tau\) were bounded, a fixed \(\tau\) would make every \(\tau t_j\) an integer, and then its multiples work). So \(\sup_k|a_{k+\tau}-a_k|<\varepsilon\). Hence \(|\hat m_\tau|\ge a_0-\varepsilon-m_c\ge m(\mathbb T)-2\delta-\varepsilon\), and \(|\hat m_{k_0+\tau}|\le|a_{k_0}|+\varepsilon+m_c\le|a^{\rm cert}_{k_0}|+\delta+\varepsilon\), for arbitrarily large \(\tau\). The two numbers come from BZ.5's enclosures in `arb`. For the \(\hat\theta\) statement, \(|\sin\pi(\gamma+i\omega k)|^2=\sin^2\pi\gamma+\sinh^2\pi\omega k\) and \(\pi\omega=2\pi^2/\ln5\). A \(1\)-periodic function holomorphic on a wider strip \(|\operatorname{Im}t|<b'\) would have \(|\hat\theta_k|=O(e^{-2\pi b''k})\) for every \(b''<b'\). \(\square\)

**KEPT FAILURE BZ.6 (certified witness; decides Conjecture BV.6).** Claim: "\(|\hat\theta_k|\) is asymptotic to a constant times \(e^{-2\pi^2k/\ln5}\)", equivalently (BX.3) "\(\rho_k=|\hat m_k|/\hat m_0\) converges to a nonzero limit". **The claim is false.** **Witness:** \(\delta\le4.03\cdot10^{-5}\) from \(311\) certified atoms, \(|a^{\rm cert}_{35}|/m(\mathbb T)\le0.59481\), so \(\limsup\rho_k\ge0.99968>0.59497\ge\liminf\rho_k\). The other two parts of BV.6 are **true** (BZ.6): the rate is sharp (\(\limsup|\hat\theta_k|e^{2\pi^2k/\ln5}\ge1.00918>0\)), and \(F\) does not continue past the strip. BY's band \([0.60,0.92]\) for \(k\le21\) is therefore not the long-run range: \(\rho_k\) returns above \(0.9996\) infinitely often.

**KEPT FAILURE BZ.5 (certified witness).** Claim (a): "\(\nu\), equivalently \(m\), is absolutely continuous". Claim (b): "\(\nu\) is carried by \(r\le310\) geometric series \(\lambda_j5^{\mathbb Z}\)". **Both are false.** For (a), BZ.4 removes every absolutely continuous part. Independently of BZ.4, the atom \(\lambda=3.615256719887703\ldots\) with \(\nu(\{\lambda\})=0.479640692929\ldots>0\) is a witness. For (b), BZ.5 gives \(311\) distinct atoms in one period, which sharpens Kept Failure BY.3 (\(19\)) and BX.5 (\(16\)). By Szegő's theorem, which is not used for any claim above, BZ.4 also explains BY's stall: for a singular \(m\), \(\det T_{n+1}/\det T_n\to0\), so the Toeplitz certificates become ever more expensive.

**Observation (no claim).** With grids of \(4000\), \(16000\) and \(128000\) cells, the uncaptured mass of \(m\) is \(6.9\cdot10^{-4}\), \(1.8\cdot10^{-4}\) and \(4.0\cdot10^{-5}\), and the number of certified atoms is \(32\), \(96\) and \(311\). This is consistent with \(\nu\) being purely atomic, but \(\nu(-Z)=0\) is not proved here. The atoms accumulate only at \(-Z\). The six largest atoms per period carry about \(98\%\) of \(m\). The bounds of BZ.6 are not optimised: \(\min_{k\le100}\) and the recurrence to \(k=0\) are the only inputs.

**What stays open.**
- *From BX/BY/BZ.* Whether \(\nu(-Z)=0\), that is, whether \(\nu\) is pure point (the remaining mass is \(\le4.03\cdot10^{-5}\) of \(m\)); whether there are infinitely many atoms per period. The exact values of \(\limsup\rho_k\) (\(=1\) if \(\nu\) is pure point, by the same recurrence) and \(\liminf\rho_k\). Whether \(\hat\theta_k\ne0\) for every \(k\ge22\): the recurrence gives \(\rho_k\ge0.9996\), hence \(\hat\theta_k\ne0\), for infinitely many \(k\), but not for every \(k\). A closed form of \(F\), \(G^*\), \(t_{\max}\), \(t_{\min}\), and of the atoms (\(\lambda_1=3.6152567\ldots\), \(\lambda_2=1.1719102\ldots\)). Feeding \(G^*\) back into BI.4/BJ. BY's \(T_{21}\) is no longer open (Remark).
- *Carried over.* Conjecture BS.4; a rate for \(J-A_N\); uniformity in \(\alpha\); eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**Remark (derived, elementary).** For any positive measure with at least \(n\) support points on the circle, \(T_n=(\hat m_{i-j})_{i,j<n}\) is positive definite, because \(x^*T_nx=\int|\sum_jx_je^{-2\pi ijt}|^2dm>0\) unless the trigonometric polynomial vanishes on the support, and a nonzero one of degree \(<n\) has at most \(n-1\) zeros per period. With BZ.5, \(T_n\succ0\) for every \(n\le311\), including BY's open \(T_{21}\). This is a consequence, not a certificate of \(T_n\)'s entries.

**NOT THIS.** Theorem BZ concerns one Stieltjes function built from a corner resolvent of a graph Laplacian on the combinatorial gasket, and the spectral measure of its scaling limit. Atoms of a spectral measure are not a flux, and nothing here is a force. Stage 1 is not touched, and no \(W\), \(0.08\), \(0.23\) or \(30\,\mu\)N/kW enters. The Stage 2 question stays OPEN. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_atoms.py`** (about 2 min; needs python-flint, sympy; no other script imported).
- *Sympy:* AQ's fixed point as \(D=cD\circ p-e\); the unrolled recursion; the conjugacy \(p+\tfrac52=(z+\tfrac52)^2-\tfrac{15}4\); the escape factor; \(|g_\pm'|\le1/\sqrt5\); \(\psi=v+v^2/20+\dots\); the far-field and tail constants.
- *Interval:* the finite-level identity for AQ's \(r_n\) (\(n=12,20,9\)); \(F(t_{\max})=G^*\) and \(F(t_{\min})=\min F\) to \(10^{-38}\); \(311\) atoms in \((1,5)\) with certified masses; \(\sum\lambda_j^{-\gamma}\nu(\{\lambda_j\})\ge0.25854\); the remainder \(\le4.1\cdot10^{-5}\); \(|a^{\rm cert}_k|\) against BY's \(\rho_1,\dots,\rho_5\) (normalisation of \(m\)); \(\limsup\rho_k\ge0.9996\) and \(\liminf\rho_k\le0.5950\); the Kept Failure witnesses.
- *Floats:* none enter a check. The grid scan only proposes candidate cells.

The script exits nonzero if any check fails.

## Theorem CA. On the backward orbit of \(-3\), \(1/r_\infty\) is exactly \(z/2\); so every \(\lambda>0\) with \(\psi(-5^k\lambda)=-3\) for some \(k\in\mathbb Z\) is an atom of \(\nu\), \(\nu\) has infinitely many atoms in one scale period (accumulating at \(\lambda^*=1.0546444769654709\ldots\), \(\psi(-\lambda^*)=-1\)), and all \(311\) atoms of Theorem BZ obey this law (Kept Failure CA.5: the other branch \(-2\) of \(p^{-1}(-6)\) carries no atom)

**Setting.** As in Theorems AQ, BD, BS, BX and BZ: \(p(z)=z(z+5)\), \(c(z)=(z+6)/((z+2)(z+5))\), \(e(z)=2z/(2+z)\), \(b=\tfrac53c\), \(D=1/r_\infty\) with \(D=c\,D\circ p-e\) on \(\mathbb C\setminus K_p\) (BZ.2), \(\psi\) (BD.4), \(Q_\infty=\int d\nu(\lambda)/(\lambda+u)\) (BX.1), \(f=1/Q_\infty=\Pi\,D(\psi)-S\) (BZ.3), \(Z=\psi^{-1}(K_p)\), \(x_j=\psi(5^{-j}u)\). Let \(O=\bigcup_{k\ge0}p^{-k}(-3)\) be the backward orbit of \(-3\) (one of the two points of \(p^{-1}(-6)=\{-2,-3\}\)), \(g_\pm(w)=-\tfrac52\pm\sqrt{w+25/4}\), and
\[
\Lambda_O=\{\lambda>0:\ \psi(-5^k\lambda)=-3\ \text{for some }k\in\mathbb Z\}.
\]
BZ left open whether \(\nu\) has infinitely many atoms per period, and a closed form of the atoms. BZ's three largest atoms per period satisfy \(\psi(-\lambda)=-3\), \(p^2(\psi(-\lambda))=-3\) and \(p^3(\psi(-\lambda))=-3\), to every printed digit. This note proves that pattern in one direction.

**Assumptions used as stated (not re-run).** AQ (1)–(3), BD.4, BJ.6, BS.1, BX.1–BX.2. BZ.1–BZ.4: the Cantor set \(K_p\subset[-5,0]\) and its level-\(n\) pieces \(g_{s_1}\circ\dots\circ g_{s_n}([-5,0])\) of length \(\le5\cdot5^{-n/2}\), with \(|g_\pm'|\le1/\sqrt5\) on \([-5,0]\); \(D=c\,D\circ p-e\); the formula \(f=\Pi D(\psi)-S\) with \(\Pi,S\) meromorphic on \(\mathbb C\) and BZ.3's tail constants; BZ.4's description of \(\nu\) off \(-Z\) (point masses at the real poles of \(Q_\infty\), with \(\nu(\{-u\})=1/f'(u)\)). BZ's enclosure routine `inv_Q` and its atom certificate are imported unchanged. No literature assumption is used.

**THEOREM CA.**
1. **CA.1 (\(D\) on the backward orbit of \(-3\)).** \(O\subset[-5,0]\setminus K_p\) consists of real points and avoids \(-2,-5,-6\). For every \(z\in O\cup\{-6\}\), \(D(z)=z/2\), that is \(r_\infty(z)=2/z\).
2. **CA.2 (atom law).** If \(u<0\) and \(\psi(u)\in O\), then \(u\notin Z\) and \(f(u)=0\). So \(-u\) is an atom of \(\nu\), with \(\nu(\{-u\})=1/f'(u)>0\). By \(\nu(5A)=3\nu(A)\), every \(\lambda\in\Lambda_O\) is an atom of \(\nu\).
3. **CA.3 (infinitely many atoms per period).** For every \(x\in Z\) at which \(\psi-\psi(x)\) changes sign (every \(x\in Z\) outside the discrete set of even-order critical points of \(\psi\)), \(-x\) is a limit of points of \(\Lambda_O\), hence of atoms of \(\nu\). Certified instance: \(\psi(-1.055)<-1<\psi(-1.054)\), \(-1\in K_p\) (\(p(-1)=-4=p(-4)\)), and \(h_n=g_+(g_-^n(-3))\in O\) are distinct with \(|h_n+1|\le5^{-(n+1)/2}\). So for each \(n\ge10\) there is \(\lambda_n\in(1.054,1.055)\) with \(\psi(-\lambda_n)=h_n\), and these are infinitely many distinct atoms of \(\nu\) in one period, accumulating at \(\lambda^*=1.0546444769654709152\ldots\) (\(\psi(-\lambda^*)=-1\), so \(-\lambda^*\in Z\)). The first members, \(n=1,2\), are BZ's \(\lambda=1.171910250998837\ldots\) and \(1.016316417446789\ldots\); \(n=0\) gives \(0.72305134397754\ldots=\lambda_1/5\), with BZ's largest atom \(\lambda_1=3.6152567198877\ldots\) satisfying \(\psi(-\lambda_1)=-3\).
4. **CA.4 (BZ's atoms obey the law; certified).** BZ.5's scan and certificate, re-run unchanged (\(128000\) cells, \(300\) bits), again give \(311\) atoms in \((1,5)\). For each one, the certificate's box is a region where \(f\) is holomorphic, so on its real segment \(f\) is strictly increasing and has exactly one zero. Inside that segment, `arb` finds a sign change of \(\psi(5^k u)+3\) for some \(0\le k\le14\). By CA.2 that point is a zero of \(f\), so it is the atom. Every one of the \(311\) atoms is in \(\Lambda_O\).

*Proof.* *CA.1.* Sympy: \(c(-6)=0\), so \(D(-6)=c(-6)D(6)-e(-6)=-3\), since \(D(6)=1/r_\infty(6)\) is finite (\(r_\infty(6)\in(0,1/6]\)). Next, \(c(z)\,p(z)/2-e(z)=z/2\) identically. So if \(D(p(z))=p(z)/2\), and \(z\notin K_p\cup\{-2,-5\}\), then \(D(z)=z/2\). Here \(p(-3)=-6\), and BZ.2's identity holds at every point of \(\mathbb C\setminus K_p\) where its terms are finite. Induction on \(k\) along \(p^{-k}(-3)\) finishes the proof. For \(w\in[-5,0]\), \(p^{-1}(w)=-\tfrac52\pm\sqrt{25/4+w}\subset[-5,0]\) is real. \(p(-3)=-6<-5\), so \(-3\notin K_p\), and neither is any point of \(O\) (\(p^{-1}(K_p)=K_p\)). \(-2\notin O\) because \(p(-2)=-6\notin O\); \(-5\notin O\) because \(p(-5)=0\) is fixed; and \(-6\notin[-5,0]\).
*CA.2.* \(x_0=\psi(u)\in O\), so \(u\notin Z\), and \(x_j\in p^{-j}(x_0)\subset O\) for every \(j\). So \(b(x_j)\) and \(e(x_j)\) are finite and nonzero, \(D(x_j)=x_j/2\), and \(u\) is not a pole of \(\Pi\) or \(S\). BZ.3's unrolling (sympy) gives \(\Pi(u)D(x_0)-\sum_{k\le M}(5/3)^ke(x_k)\Pi(5^{-k}u)=\Pi(5^{-M}u)(5/3)^MD(x_M)\). Hence
\[
f(u)=\Pi(5^{-M}u)\,\bigl(\tfrac53\bigr)^M\frac{x_M}2-\sum_{k>M}\bigl(\tfrac53\bigr)^ke(x_k)\,\Pi(5^{-k}u).
\]
With \(|x_M|\le2\cdot5^{-M}|u|\), \(\Pi(5^{-M}u)\to1\) and BZ.3's tail constants (\(|e(x)|\le1.06|x|\)), both terms are \(O(3^{-M})\). So \(f(u)=0\). \(f\not\equiv0\) is meromorphic near \(u\), so \(u\) is a pole of \(Q_\infty\) and, by BZ.4, an atom with mass \(1/f'(u)\). The positivity \(f'(u)>0\) holds because the residue of \(Q_\infty\) at \(u\) is \(\nu(\{-u\})>0\).
*CA.3.* Let \(y=\psi(x)\in K_p\). Since \(\psi-y\) changes sign at \(x\), \(\psi([x-\varepsilon,x+\varepsilon])\) contains an interval around \(y\). For large \(n\), it therefore contains the level-\(n\) piece \(g_s([-5,0])\ni y\), and with it \(g_s(-3)\in O\), which is \(\ne y\). The intermediate value theorem gives \(u\in(x-\varepsilon,x+\varepsilon)\) with \(\psi(u)\in O\), and CA.2 applies. For the instance: \(g_-(-4)=-4\), \(g_+(-4)=-1\), and \(|g_\pm'|\le1/\sqrt5\) on \([-5,0]\) give \(|h_n+1|\le5^{-(n+1)/2}\), which for \(n\ge10\) is below \(\min(-1-\psi(-1.055),\,\psi(-1.054)+1)\). The \(h_n\) are distinct, since \(g_\pm\) are injective and \(-3\) is not periodic. So the \(\lambda_n\) are distinct.
*CA.4.* On the box, \(f\) is holomorphic, so \(Q_\infty\) has no zero there. Off its poles, \(Q_\infty'=-\int d\nu/(\lambda+u)^2<0\) (BZ.4: no accumulation of atoms off \(-Z\)), and at its poles \(f'=1/\nu(\{\cdot\})>0\). So \(f'>0\) on the real segment. \(\square\)

**KEPT FAILURE CA.5 (certified witness).** Claim: "every \(\lambda\) with \(\psi(-5^k\lambda)\in p^{-1}(-6)=\{-2,-3\}\) is an atom of \(\nu\)" (both branches of the decimation preimage of \(-6\) seed atoms). **The claim is false.** **Witness:** \(u=-2.242133185246453863226262\ldots\) (sign change of \(\psi(u)+2\), `arb`) has \(\psi(u)=-2\notin K_p\), and \(f(u)\in[2.1979467595912,\,2.1979467595913]\). So \(Q_\infty(u)\) is finite and nonzero, and \(\lambda=2.2421331852464538\ldots\) is not an atom. At \(z=-2\), BZ's routine meets \(0/0\). Two steps of \(D=c\,D\circ p-e\) give the regular form \(D(z)=\frac{(z+6)(z+3)D(p^2z)}{(z+5)(p(z)+2)(p(z)+5)}-\frac{2z(z+4)}{p(z)+2}\), \(D(-2)=D(6)/3-2\) (sympy), which agrees with BZ's routine to \(10^{-40}\) at \(u=-2.3\). The \(-2\) branch fails because \(c\) and \(e\) have their poles there, so CA.1's induction cannot start.

**Conjecture CA.6 (not claimed).** The atoms of \(\nu\) are exactly \(\Lambda_O\), and \(\nu\) is carried by them (\(\nu(-Z)=0\), pure point). Evidence: CA.4, all \(311\) certified atoms, and BZ's uncaptured mass \(\le4.03\cdot10^{-5}\). Not proved: CA.2 is one direction only.

**Observation (no claim).** Of BZ's \(311\) atoms, the number with each \(k\) (smallest \(k\) with \(\psi(-5^k\lambda)=-3\)) is \(1,0,1,3,5,10,20,39,72,94,53,10,2,1\) for \(k=0,\dots,13\). The counts roughly double up to \(k=8\) and then fall off. That is consistent with the \(128000\)-cell grid missing the smaller, closer-packed deep atoms, but this is not proved. The pattern matches the \(6\)-series of spectral decimation on the Sierpiński gasket, in which preimages of \(6\) through \(3\), and not through the forbidden value \(2\), give eigenvalues. That identification is not used and not claimed.

**What stays open.**
- *From BZ/CA.* Conjecture CA.6. Whether \(\nu(-Z)=0\). Whether every point of \(-Z\) (including even-order critical points of \(\psi\)) is a limit of atoms. A closed form for the masses \(1/f'(u)\) on \(\Lambda_O\). Exact \(\limsup\rho_k\) and \(\liminf\rho_k\). Whether \(\hat\theta_k\ne0\) for every \(k\ge22\). Closed forms of \(F\), \(G^*\), \(t_{\max}\), \(t_{\min}\). \(G^*\) into BI.4/BJ.
- *Carried over.* Conjecture BS.4; a rate for \(J-A_N\); uniformity in \(\alpha\); eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem CA locates atoms of the spectral measure of one scaled corner resolvent of a graph Laplacian on the combinatorial gasket. Atoms of a spectral measure are not a flux, and nothing here is a force. Stage 1 is not touched, and no \(W\), \(0.08\), \(0.23\) or \(30\,\mu\)N/kW enters. The Stage 2 question stays OPEN: the repository has no \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh from which to decide it. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_atom_law.py`** (about 2.3 min; needs python-flint, sympy; imports BZ's `gasket_fractional_profile_atoms.py` unchanged).
- *Sympy:* \(c\,p/2-e=z/2\); \(c(-6)=0\), \(D(-6)=-3\); \(p^{-1}(-6)=\{-2,-3\}\), \(p^{-1}(-3)\); \(p^{-1}([-5,0])\subset[-5,0]\); the exclusions \(-2,-5,-6\notin O\); BZ.3's unrolling; the regular form of \(D\) at \(-2\) and \(D(-2)=D(6)/3-2\); \(g_\pm(-4)\), \(g_-'(-4)=-1/3\).
- *Interval (`arb`, 300 bits, outward):* \(\psi(-1.055)<-1<\psi(-1.054)\); \(\lambda^*\); \(|h_n+1|\le5^{-(n+1)/2}\) and \(h_n\in(\psi(-1.055),\psi(-1.054))\) for \(10\le n\le40\), distinct; \(|f|<10^{-60}\) at the points of \(\Lambda_O\) for \(n=0,\dots,7\), with \(n=1,2\) equal to BZ's atoms (consistency); \(\psi(-\lambda_1)=-3\) for BZ's largest atom; the CA.5 witness; CA.4 for all \(311\) atoms.
- *Floats:* none enter a check. Newton/grid steps only propose brackets.

The script exits nonzero if any check fails.

## Theorem CB. \(1/r\) is an explicit infinite product: \(H=1/r-z/2\) obeys \(H=c\,H\circ p\), so \(1/r_\infty(z)=\tfrac z2\bigl(1+\prod_{j\ge0}\tfrac{p^jz+6}{p^jz+2}\bigr)\) and \(1/Q_\infty=\Pi\cdot\tfrac\psi2\prod_{j\ge0}\tfrac{\psi(5^j\cdot)+6}{\psi(5^j\cdot)+2}\) (BZ.3's \(S\) equals \(\Pi\psi/2\)); hence the atoms of \(\nu\) off \(-Z\) are exactly \(\Lambda_O\) (Conjecture CA.6 holds off \(-Z\)), no point of the \(-2\) branch is an atom, and the atom masses have a closed product form

**Setting.** As in Theorems AQ, BD, BS, BX, BZ and CA: \(r_n(s)=\tfrac12v^T(L_n+s)^{-1}v\), \(r_{n+1}(s)=\Phi_s(r_n(p(s)))\), \(r_0=1/(3+s)\) (AQ); \(p(z)=z(z+5)\), \(c(z)=(z+6)/((z+2)(z+5))\), \(e(z)=2z/(2+z)\), \(D=1/r_\infty\) with \(D=c\,D\circ p-e\) on \(\mathbb C\setminus K_p\) (BZ.2); \(\psi\) (BD.4); \(f=1/Q_\infty=\Pi\,D(\psi)-S\) (BZ.3); \(Z=\psi^{-1}(K_p)\); \(O=\bigcup_kp^{-k}(-3)\) and \(\Lambda_O=\{\lambda>0:\psi(-5^k\lambda)=-3\text{ for some }k\in\mathbb Z\}\) (CA). CA proved that \(\Lambda_O\) consists of atoms and left open whether there are others (Conjecture CA.6) and a closed form of the masses. CA.1's identity \(c\,p/2-e=z/2\) says that \(z/2\) solves \(D\)'s functional equation. This note uses that globally, not only on \(O\).

**Assumptions used as stated (not re-run).** AQ (1)–(3) (decimation recursion, \(r_n\uparrow r_\infty\)); BD.4, BJ.6 (\(\psi(5v)=p(\psi(v))\), nonnegative coefficients, \(v\le\psi(v)\le ve^{v/10}\)); BS.1; BX.1–BX.2; BZ.1–BZ.4 (\(K_p\) a null Cantor set in \([-5,0]\) with \(p^n\to\infty\) locally uniformly off it; \(D=c\,D\circ p-e\) and \(D\) meromorphic off \(K_p\); \(f=\Pi D(\psi)-S\) with \(\Pi,S\) meromorphic on \(\mathbb C\) and BZ.3's tail constants; BZ.4's description of \(\nu\) off \(-Z\): point masses at the real zeros of \(f\), mass \(1/f'\)); CA.1–CA.2. BW's digits of \(G^*,\min F,t_{\max},t_{\min}\) and BZ.5's printed masses are used for consistency checks only. No literature assumption is used.

**THEOREM CB.**
1. **CB.1 (finite level, exact).** For every \(n\ge0\),
   \[
   \frac1{r_n(s)}=\frac s2\Bigl[1+\Bigl(1+\frac6{p^n(s)}\Bigr)\prod_{j<n}\frac{p^j(s)+6}{p^j(s)+2}\Bigr].
   \]
2. **CB.2 (the corner resolvent in closed product form).** On \(\mathbb C\setminus K_p\),
   \[
   \frac1{r_\infty(z)}=D(z)=\frac z2\bigl(1+P(z)\bigr),\qquad P(z)=\prod_{j\ge0}\frac{p^j(z)+6}{p^j(z)+2},
   \]
   with \(P\) meromorphic on \(\mathbb C\setminus K_p\). Equivalently \(H=D-z/2=\tfrac z2P\) is the solution of \(H(z)=c(z)H(p(z))\) with \(H(w)\sim w/2\) at \(\infty\). CA.1 (\(D=z/2\) on \(O\cup\{-6\}\)) is the zero set of \(P\) on the \(-3\) branch, and CA.5's \(D(-2)=D(6)/3-2\) is the cancelled pair \((z+6)(pz+6)/((z+2)(pz+2))=(z+6)(z+3)/(pz+2)=-1\) at \(z=-2\).
3. **CB.3 (\(1/Q_\infty\) in closed product form).** \(S(u)=\Pi(u)\psi(u)/2\) identically, so on \(\mathbb C\setminus Z\), with \(x_i=\psi(5^iu)\) (\(i\in\mathbb Z\)),
   \[
   \frac1{Q_\infty(u)}=\Pi(u)\,\frac{\psi(u)}2\,P(\psi(u))=\frac{x_0}2\prod_{i\ge0}\frac{x_i+6}{x_i+2}\ \prod_{i\le-1}\frac{5}{3}\,\frac{x_i+6}{(x_i+2)(x_i+5)} .
   \]
   No far-field bound for \(D\) and no sum \(S\) enter.
4. **CB.4 (the atoms of \(\nu\) off \(-Z\); Conjecture CA.6 off \(-Z\)).** For \(u<0\), \(u\notin Z\): \(f(u)=0\) if and only if \(\psi(5^ku)=-3\) for some \(k\in\mathbb Z\). Hence the set of atoms of \(\nu\) in \((0,\infty)\setminus(-Z)\) is exactly \(\Lambda_O\). In particular no \(\lambda\) with \(\psi(-5^k\lambda)\in\bigcup_{j\ge0}p^{-j}(-2)\) is an atom. This extends CA.5 from one point to the whole \(-2\) branch.
5. **CB.5 (closed product form of the masses).** If \(\psi(5^ku)=-3\) with \(k\ge0\) and \(u<0\), then the atom \(\lambda=-u\) has
   \[
   \nu(\{\lambda\})=\frac{-4}{\Pi(u)\,\frac{x_0}2\,\prod_{i\ge0,\,i\ne k+1}\frac{x_i+6}{x_i+2}\,\prod_{i=0}^{k}(2x_i+5)\,\psi'(u)},\qquad \psi'(u)=\prod_{i\ge1}\Bigl(1+\frac25\psi(5^{-i}u)\Bigr).
   \]
   It is a product over the orbit \(x_i=\psi(5^iu)\) alone, and with \(\nu(\{5\lambda\})=3\nu(\{\lambda\})\) it covers every atom. Certified: \(\nu(\{3.6152567198877032\ldots\})=0.4796406929290034734083466\ldots\) (\(k=0\)), \(\nu(\{1.1719102509988374\ldots\})=0.03592210432462168505\ldots\) (\(k=2\)), \(\nu(\{1.0163164174467888\ldots\})=0.01362783039371924943\ldots\) (\(k=3\)), matching BZ.5's Cauchy-estimate enclosures.

*Proof.* *CB.1.* AQ's recursion is \(1/r_{n+1}=c/r_n(p)-e\) (BZ.2, sympy). Put \(H_n=1/r_n-s/2\). Since \(c\,p/2-e=s/2\) (CA.1), \(H_{n+1}(s)=c(s)H_n(p(s))\), and \(H_0=3+s/2=\tfrac s2(1+6/s)\). So \(H_n(s)=\prod_{j<n}c(p^js)\cdot H_0(p^ns)\). With \(c(w)p(w)/w=(w+6)/(w+2)\), the product telescopes: \(\prod_{j<n}c(p^js)\,p^n(s)=s\prod_{j<n}(p^js+6)/(p^js+2)\). Sympy checks this against AQ's recursion for \(n\le4\). Independently of AQ, exact rational Gaussian elimination on the gasket graphs \(G_1,G_2,G_3\) (\(6,15,42\) vertices, built from coordinates) gives the same \(r_n(s)\) at \(s=\tfrac12,3\).
*CB.2.* For \(s>0\), \(p^n(s)\uparrow\infty\) and \(\sum_j4/(p^js+2)<\infty\), so CB.1 and \(r_n\uparrow r_\infty\) (AQ (3)) give the formula on \((0,\infty)\). On a compact subset of \(\mathbb C\setminus K_p\), \(p^j\to\infty\) uniformly (BZ.1). Once \(|w|\ge10\): \(|p(w)|\ge2|w|\), \(|4/(w+2)|\le\tfrac12\) and \(|\log(1+t)|\le2|t|\), so the tail \(\prod_{j\ge N}\) lies in \(\exp(B(16/(|p^Nz|-2)))\). Hence \(P\) is meromorphic on \(\mathbb C\setminus K_p\). That set is connected (\(K_p\) is a compact subset of \(\mathbb R\)), and \(D\) is meromorphic there (BZ.2), so the identity extends. Zeros and poles: a factor vanishes iff \(p^jz=-6\) and has a pole iff \(p^jz=-2\). Since \(p(z)+6=(z+2)(z+3)\), a pole at \(p^jz=-2\) is cancelled to the same order by the zero of the next factor.
*CB.3.* By CA.1, \(Y_k=(5/3)^kx_k/2\) (\(x_k=\psi(5^{-k}u)\)) obeys BZ.3's recursion \(Y_k=b_kY_{k-1}-(5/3)^ke(x_k)\). The unrolling (sympy, as in CA.2) gives \(\Pi(u)\psi(u)/2-S(u)=\Pi(5^{-M}u)(5/3)^Mx_M/2-\sum_{k>M}(5/3)^ke(x_k)\Pi(5^{-k}u)\), which is \(O(3^{-M})\) by BZ.3's tail constants. This holds at every \(u\) where no \(x_j\) (\(j\ge1\)) is \(-2\) or \(-5\). That set is dense, and both sides are meromorphic, so \(S=\Pi\psi/2\). Then \(f=\Pi\,(D(\psi)-\psi/2)=\Pi\,\tfrac\psi2P(\psi)\) by CB.2, with \(p^j(\psi(u))=\psi(5^ju)\).
*CB.4.* Let \(u<0\), \(u\notin Z\). Then every \(x_i\notin K_p\) (\(5Z=Z\)), so \(x_0\ne0\), and \(x_i\ne-5\) for \(i\le-1\) (else \(x_{i+1}=0\in K_p\)). Near \(u\), \(x_i\to\infty\) (\(i\to+\infty\)) and \(x_i\to0\) (\(i\to-\infty\)) uniformly, so only finitely many factors of CB.3 have zeros or poles near \(u\), and the remaining product is holomorphic and nonzero. Thus \(\operatorname{ord}_uf\) is the sum of the factors' orders. A factor vanishes only where \(x_i=-6\), and then \(x_{i-1}\in p^{-1}(-6)=\{-2,-3\}\). If \(x_{i-1}=-3\), then \(\psi(5^{i-1}u)=-3\). The forward orbit \(6,66,\dots\) of \(-6\) is positive and the backward orbit lies in \(O\), which avoids \(-2,-5,-6\) (CA.1). So no factor has a pole and \(f(u)=0\), as in CA.2. If \(x_{i-1}=-2\), then \(\operatorname{ord}(x_i+6)=\operatorname{ord}(x_{i-1}+2)\) cancels the pole of factor \(i-1\). The forward orbit of \(-6\) is positive, and \(p^{-1}(-2)=(-5\pm\sqrt{17})/2\) and its preimages never meet \(\{-2,-3,-6\}\) (a second \(-2\) or \(-6\) below would force a positive point into the orbit). So the orders cancel and \(f(u)\ne0\). If no \(x_i\) is \(-6\) or \(-2\), all factors are finite and nonzero, so \(f(u)\ne0\). By BZ.4, the atoms of \(\nu\) off \(-Z\) are the \(-u\) with \(f(u)=0\).
*CB.5.* At such \(u\), the only vanishing factor is \(i=k+1\) (\(x_{k+1}=p(-3)=-6\)), and \(\frac d{dx}\frac{x+6}{x+2}\big|_{-6}=\frac1{x+2}\big|_{-6}=-\tfrac14\). Also \(\frac d{du}x_{k+1}=\prod_{i=0}^k(2x_i+5)\cdot\psi'(u)\). The derivative of \(\psi(5w)=p(\psi(w))\) gives \(\psi'(u)=\psi'(5^{-J}u)\prod_{i=1}^J(1+2\psi(5^{-i}u)/5)\). Here \(|\psi'(v)-1|\le\psi'(|v|)-1\le2(e^{|v|/5}-1)\) by nonnegative coefficients, convexity (\(\psi'(t)\le(\psi(2t)-\psi(t))/t\)) and BJ.6. Then \(\nu(\{-u\})=1/f'(u)\) (BZ.4). \(\square\)

**Conjecture CA.6, status after CB (not claimed in full; now PROVED by Theorem CC).** Off \(-Z\), CA.6 is now proved: the atoms there are exactly \(\Lambda_O\). Still open are whether \(\nu\) has atoms on \(-Z\) and whether \(\nu(-Z)=0\) (no singular continuous part on the null set \(-Z\)). BZ's uncaptured mass \(\le4.03\cdot10^{-5}\) is the only quantitative evidence.

**Observation (no claim).** For \(s>0\), CB.1 is a spectral-decimation product of the type known for the Sierpiński gasket. Its factors \((w+6)/(w+2)\) carry the \(6\)- and \(2\)-series that CA's Observation mentions. A comparison with the literature was not made, and no literature result is used. CB.2–CB.3 do not give closed forms of \(G^*\), \(F\), \(t_{\max}\), \(t_{\min}\), because \(\psi\) itself is not in closed form.

**What stays open.**
- *From CA/CB.* Whether \(\nu(-Z)=0\) and whether \(-Z\) carries atoms (the rest of CA.6). Accumulation of atoms at even-order critical points of \(\psi\) in \(Z\). A closed form of \(\psi\), hence of \(F\), \(G^*\), \(t_{\max}\), \(t_{\min}\). The exact \(\limsup\rho_k\) and \(\liminf\rho_k\). Whether \(\hat\theta_k\ne0\) for every \(k\ge22\). \(G^*\) into BI.4/BJ.
- *Carried over.* Conjecture BS.4; a rate for \(J-A_N\); uniformity in \(\alpha\); eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem CB gives a product formula for one corner resolvent of a graph Laplacian on the combinatorial gasket, for its scaling limit, and for the atoms of that limit's spectral measure. A resolvent is not a flux, and nothing here is a force. Stage 1 is not touched, and no \(W\), \(0.08\), \(0.23\) or \(30\,\mu\)N/kW enters. The Stage 2 question stays OPEN: the repository has no \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh from which to decide it. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_product.py`** (about 10 s; needs python-flint, sympy; imports BZ's and CA's scripts unchanged).
- *Sympy:* \(H=c\,H\circ p\); \(c(w)p(w)/w=(w+6)/(w+2)\); \(H_0\); CB.1 against AQ's recursion for \(n\le4\); \(p+6=(z+2)(z+3)\) and the regular pair at \(-2\); \(p^{-1}(-2)\); the derivatives \(-\tfrac14\) and \(\tfrac5{12}\) at \(-6\); \(\psi'\) recursion; tail constants.
- *Exact rationals:* \(r_n(s)\) of the graphs \(G_1,G_2,G_3\) at \(s=\tfrac12,3\) equal CB.1.
- *Interval (`arb`/`acb`, 300–400 bits, outward):* \((z/2)(1+P)\) against BZ's `D_of` to \(10^{-50}\) (8 real, 3 complex points); \(\Pi\tfrac\psi2P(\psi)\) against BZ's `inv_Q` to \(10^{-40}\) (6 real, 2 complex points); \(F(t_{\max})=G^*\) and \(F(t_{\min})=\min F\) to \(10^{-38}\) from the product alone; CA.5's witness \(f\in[2.1979467595912,2.1979467595913]\) at \(\psi(u)=-2\) from the product alone; \(f(u)=1.31876805575\ldots\ne0\) at \(u=-0.44842663704929\ldots\), \(\psi(u)=(-5+\sqrt{17})/2\); CB.5's masses for \(k=0,2,3\) against BZ.5 (to \(2\cdot10^{-12}\)) and against BZ's `certify_atom`.
- *Floats:* none enter a check. A float grid only proposes one bracket.

The script exits nonzero if any check fails.

## Theorem CC. BX's spectral measure \(\nu\) is pure point: \(\nu(-Z)=0\), so \(\nu=\sum_{\lambda\in\Lambda_O}\nu(\{\lambda\})\delta_\lambda\) and Conjecture CA.6 holds in full; the CB.5 masses of one scale period add up exactly to \(m(\mathbb T)=c_R(\gamma)/2\), and \(\limsup\rho_k=1\), so BX.3's sharp bound \(\pi c_R(\gamma)/\ln5\) is attained

**Setting.** As in Theorems AQ, BD, BS, BX, BZ, CA and CB: \(r_\infty\) is the corner resolvent, \(r_\infty(z)=\int d\mu(\lambda)/(\lambda+z)\) with \(\mu\) a probability measure on \([0,8]\) (AQ (2)–(3), BZ.2); \(D=1/r_\infty\); \(H=D-z/2=\tfrac z2P(z)\) with \(H=c\,H\circ p\) on \(\mathbb C\setminus K_p\) (CB.2), \(c(z)=(z+6)/((z+2)(z+5))\), \(p(z)=z(z+5)\); \(\psi\) (BD.4); \(\Pi(u)=\prod_{j\ge1}\tfrac53c(\psi(5^{-j}u))\); \(f=1/Q_\infty=\Pi\,H(\psi)\) on \(\mathbb C\setminus Z\) (CB.3); \(Q_\infty=\int d\nu(\lambda)/(\lambda+u)\) (BX.1); \(K_p\subset[-5,0]\) the null Cantor Julia set (BZ.1), \(Z=\psi^{-1}(K_p)\); \(O\), \(\Lambda_O\) (CA); \(m=\lambda^{-\gamma}\nu\) on the scale circle, \(m(\mathbb T)=c_R(\gamma)/2\) (BX.2); \(\rho_k=|\hat m_k|/\hat m_0\). CB proved CA.6 off \(-Z\) and left open whether \(-Z\) carries mass. This note closes that. The tool is a transfer identity on \(K_p\) for the Herglotz measure of \(-1/H\), whose total weight is \((2-t)/(6-t)\le7/11<1\).

**Assumption CC (literature).** The Herglotz–Nevanlinna representation: a holomorphic \(F\) on the upper half plane with \(\operatorname{Im}F>0\) has the form \(F(z)=a+bz+\int\bigl(\frac1{t-z}-\frac t{1+t^2}\bigr)d\sigma(t)\), with \(b\ge0\) and \(\sigma\ge0\) a Radon measure. It belongs to the same family as Assumption BX. Stieltjes inversion for finite signed measures is elementary: if the Cauchy transform of \(\tau\) continues holomorphically across an interval \(I\), then \(\tau|_I=0\), since \(\operatorname{Im}G_\tau(x+i\varepsilon)\to\pi\tau\) weakly and the reflection \(G_\tau(\bar z)=\overline{G_\tau(z)}\) forces the boundary values to vanish.

**Assumptions used as stated (not re-run).** AQ (1)–(3), BD.4, BJ.6, BS.1, BX.1–BX.2, BZ.1–BZ.4, CA.1–CA.3, CB.1–CB.5, together with BZ's, CA's and CB's enclosure routines (`psi_enc`, `f_prod`, `D_prod`, `mass_closed`, imported unchanged). Theorems AT/AU's \(c_R(\gamma)\) digits (as printed in BZ's script) are used for the mass comparison.

**THEOREM CC.**
1. **CC.1 (\(-1/H\) is a Nevanlinna function).** For \(\operatorname{Im}z>0\), \(\operatorname{Im}D(z)\ge\operatorname{Im}z\), so \(\operatorname{Im}H(z)\ge\tfrac12\operatorname{Im}z\). Hence \(F=-1/H\) is holomorphic with \(\operatorname{Im}F>0\) on the upper half plane, \(|F(w)|\le2/|\operatorname{Im}w|\) off \(\mathbb R\), and \(F=k\cdot F\circ p\) with \(k=1/c=(z+2)(z+5)/(z+6)\).
2. **CC.2 (no mass on the Julia set).** Let \(\sigma\) be the Herglotz measure of \(F\) (Assumption CC). Then \(\sigma(K_p)=0\). The reason is that \(\sigma|_{K_p}=\Phi(\sigma|_{K_p})\), where \((\Phi\tau)(B)=\int\sum_{y\in p^{-1}(t)\cap B}\frac{k(y)}{p'(y)}\,d\tau(t)\), and
   \[
   T(t)=\sum_{p(y)=t}\frac{k(y)}{p'(y)}=\frac{2-t}{6-t}\in\bigl[\tfrac13,\tfrac7{11}\bigr]\quad(t\in[-5,0]).
   \]
3. **CC.3 (\(\nu\) is pure point; Conjecture CA.6).** \(\nu((0,\infty)\cap(-Z))=0\). Hence \(\nu=\sum_{\lambda\in\Lambda_O}\nu(\{\lambda\})\,\delta_\lambda\) with the masses of CB.5, the atom set of \(\nu\) is exactly \(\Lambda_O\), and \(Q_\infty(u)=\sum_{\lambda\in\Lambda_O}\nu(\{\lambda\})/(\lambda+u)\). The same holds for \(m\) on the circle. Conjecture CA.6 is proved.
4. **CC.4 (one period, exactly; certified).** \(\psi\) is strictly increasing on \([-5,0]\). So \(\Lambda_O\cap[1,5)=\{-\psi^{-1}(t):t\in O\cap(\psi(-5),\psi(-1)]\}\), with \(\psi(-1)=-0.950826643829326\ldots\) and \(\psi(-5)=p(\psi(-1))=-3.850061912530890\ldots\), and
   \[
   \frac{c_R(\gamma)}2=\sum_{k\ge0}\ \sum_{t\in p^{-k}(-3)\cap(\psi(-5),\psi(-1)]}\lambda_t^{-\gamma}\,\nu(\{\lambda_t\})\qquad(\lambda_t=-\psi^{-1}(t)),
   \]
   exactly, with each term in CB.5's closed product form. Certified (python-flint `arb`, 300 bits): through level \(k\le14\) there are exactly \(10240\) atoms in \([1,5)\), with level counts \(1,0,1,3\) and then \(5\cdot2^{k-4}\). Their sum is \(0.2585821335206860\ldots\) (\(\pm1.7\cdot10^{-17}\)), so the levels \(>14\) carry exactly \(\delta_{14}=7.837945666\ldots\cdot10^{-7}\) of \(m(\mathbb T)=0.25858291731525262196\ldots\). The level sums are \(0.199492597901\ldots\) (\(k=0\), one atom, \(\nu=0.4796406929290034734\ldots\)), \(0\), \(0.0322354737761\ldots\), \(0.0170039956381\ldots\), \(0.00564716030397\ldots\), …, \(1.06594981300\ldots\cdot10^{-6}\) (\(k=14\)).
5. **CC.5 (the normalised moduli).** \(\limsup_k\rho_k=1\), so \(\limsup_k|\hat\theta_k|e^{2\pi^2k/\ln5}=\pi c_R(\gamma)/\ln5=1.00949802052669\ldots\). BX.3's sharp bound is attained exactly; BZ.6 had \(\ge1.00918\). Certified with the \(10240\) atoms and \(\delta_{14}\): \(\liminf_k\rho_k\le0.59486582\) (\(k=35\); BZ.6: \(0.59497\)).
6. **CC.6 (support).** \(\psi\) has no critical point on \(Z\). So CA.3 applies at every point: every point of \(-Z\cap(0,\infty)\) is a limit of atoms, and \(\operatorname{supp}\nu=\Lambda_O\cup(-Z\cap(0,\infty))\). The question of accumulation at even-order critical points of \(\psi\) in \(Z\), left open by CA and CB, is empty.

*Proof.* *CC.1.* Write \(r=r_\infty(z)=\int d\mu/(\lambda+z)\). Then \(-\operatorname{Im}r=\operatorname{Im}z\int d\mu/|\lambda+z|^2\), and by Cauchy–Schwarz (\(\mu(\mathbb R)=1\)) \(|r|^2\le\int d\mu/|\lambda+z|^2\). So \(\operatorname{Im}D=-\operatorname{Im}r/|r|^2\ge\operatorname{Im}z\). In particular \(H\ne0\) on the upper half plane, and \(|F|=1/|H|\le1/\operatorname{Im}H\le2/\operatorname{Im}z\). By real symmetry (\(r_\infty(\bar z)=\overline{r_\infty(z)}\) on the connected set \(\mathbb C\setminus K_p\)) the same bound holds below the axis with \(|\operatorname{Im}w|\). \(F=k\,F\circ p\) is \(H=c\,H\circ p\) inverted (sympy).
*CC.2.* Fix \(U=(-\tfrac{11}2,\tfrac12)\supset[-5,0]\supset K_p\), and put \(\sigma_U=\sigma|_U\), a finite measure. From Assumption CC, \(F=G_{\sigma_U}+h\) with \(G_\tau(z)=\int d\tau(t)/(t-z)\) and \(h\) holomorphic on \(\mathbb C\setminus(\mathbb R\setminus U)\). For non-real \(z\) near \(K_p\), \(p(z)\) is non-real (\(\operatorname{Im}p(x+iy)=y(2x+5)\) and \(-\tfrac52\notin K_p\)), and \(F(z)=k(z)F(p(z))\). Partial fractions (sympy) give \(1/(t-p(z))=\sum_{p(y)=t}1/(p'(y)(y-z))\). The roots \(y\) of \(p(y)=t\), \(t\in U\), lie in \(p^{-1}(U)\), which avoids \(-6\) (the pole of \(k\)) and \(-\tfrac52\) (sympy). So \(k(z)G_{\sigma_U}(p(z))=G_{\Phi\sigma_U}(z)+\int d\sigma_U(t)\sum_y\frac{k(z)-k(y)}{p'(y)(y-z)}\). The last term is holomorphic near \(K_p\), because \((k(z)-k(y))/(y-z)=-(yz+6y+6z+32)/((y+6)(z+6))\). Hence \(G_{\sigma_U-\Phi\sigma_U}\) continues holomorphically across a neighbourhood of \(K_p\) in \(\mathbb R\), and Stieltjes inversion gives \(\sigma=\Phi\sigma_U\) on Borel subsets of \(K_p\). For \(B\subset K_p\), only \(t=p(y)\in K_p\) contributes, and both roots of \(p(y)=t\in K_p\) lie in \(K_p\) (\(p^{-1}(K_p)=K_p\)). So \(\sigma(K_p)=\int_{K_p}T\,d\sigma\). The sum of residues of \(k(y)/(p(y)-t)\) at \(-6\) and \(\infty\) gives \(T=(2-t)/(6-t)\) (sympy, also from the explicit roots). Since \(T'=-4/(6-t)^2<0\), \(T\le T(-5)=7/11\) on \([-5,0]\). So \(\sigma(K_p)\le\tfrac7{11}\sigma(K_p)<\infty\), and \(\sigma(K_p)=0\). Positivity of the weights is not needed, only \(\sigma\ge0\). (They are in fact \(\ge0\) on \(p^{-1}([-5,0])\).)
*CC.3.* Let \(u_0<0\), \(u_0\in Z\), and \(y_0=\psi(u_0)\in K_p\). First, \(\psi'(u_0)=\prod_{i\ge1}(1+\tfrac25\psi(5^{-i}u_0))\) (CB.5) is a convergent product of factors that are nonzero because \(\psi(5^{-i}u_0)\in K_p\not\ni-\tfrac52\) (\(5^{-i}Z=Z\)). So \(\psi\) is a real-analytic diffeomorphism near \(u_0\) (this is also CC.6). Second, a factor of \(\Pi\) vanishes or has a pole at \(u_0\) only if some \(\psi(5^{-j}u_0)\in\{-6,-2,-5\}\). Of these, only \(-5\) lies in \(K_p\), and \(\psi(5^{-j}u_0)=-5\) forces \(\psi(u_0)=p^{j-1}(0)=0\). (i) If \(\psi(u_0)\ne0\), then \(\Pi\) is holomorphic and nonzero near \(u_0\), and \(Q_\infty=-F(\psi)/\Pi\). With \(F=G_{\sigma_V}+{}\)holomorphic for a small interval \(V\ni y_0\), the same partial-fraction step with the local inverse \(\phi=\psi^{-1}\) writes \(Q_\infty\) near \(u_0\) as the Cauchy transform of \(-\phi_*\bigl(\sigma_V/(\Pi\circ\phi\cdot\psi'\circ\phi)\bigr)\) plus a holomorphic function. By Stieltjes inversion, \(\nu\) near \(-u_0\) is the reflection of that measure. Its mass on \(Z\) is a weighted \(\sigma\)-mass of \(\psi(Z)\subset K_p\), which is \(0\) by CC.2. (ii) If \(\psi(u_0)=0\), then \(\Pi\) has a pole at \(u_0\): exactly one index has \(\psi(5^{-j}u_0)=-5\), because \(p^i(-5)=0\ne-5\), and \(-6,-2\notin K_p\). So \(|1/\Pi(u_0+i\varepsilon)|\le C\varepsilon\), and \(|F(\psi(u_0+i\varepsilon))|\le2/|\operatorname{Im}\psi(u_0+i\varepsilon)|\le C'/\varepsilon\). Thus \(Q_\infty\) stays bounded, \(\varepsilon|Q_\infty(u_0+i\varepsilon)|\to0\), and \(-u_0\) is not an atom. These \(u_0\) form a discrete set. So \(\nu((0,\infty)\cap(-Z))=0\). Off \(-Z\), \(\nu\) is the sum of its atoms (BZ.4), which are exactly \(\Lambda_O\) (CB.4).
*CC.4.* For \(v\in[-1,0]\), nonnegative coefficients and BJ.6 give \(|\psi(v)-v|\le\psi(|v|)-|v|\le e^{1/10}-1<0.106\), so \(\psi(v)>-\tfrac52\). Hence for \(u\in[-5,0]\) every factor of \(\psi'(u)=\prod_{i\ge1}(1+\tfrac25\psi(5^{-i}u))\) is positive. If \(\psi(5^ku)=-3\) with \(k<0\), then \(5^ku\in[-1,0)\), where \(\psi>-\tfrac52\). So in one period only \(k\ge0\) occurs, i.e. \(\psi(u)\in p^{-k}(-3)\). \(f(-1)=1.0966909955\ldots\ne0\) is finite (`arb`), so \(\lambda=1,5\) are not atoms. The identity is CC.3 with \(\nu(5A)=3\nu(A)\) and BX.2. For the certificate, each of the \(2^k\) points of \(p^{-k}(-3)\) (\(k\le14\)) is tested against \((\psi(-5),\psi(-1))\). Each one inside gets a ball of radius \(10^{-45}\) across which \(p^k(\psi(u))+3\) changes sign (`arb`). The balls are pairwise disjoint, and \(\psi\) is injective, so these are all the level-\(\le14\) atoms in the period. Their masses come from CB.5's `mass_closed`.
*CC.5.* Given \(\varepsilon>0\), CC.3 gives finitely many atoms carrying all of \(m(\mathbb T)\) but \(\varepsilon\). BZ.6's Dirichlet argument then gives \(\limsup\rho_k\ge1-2\varepsilon/m(\mathbb T)\), and \(\rho_k\le1\). The \(\hat\theta\) statement is BZ.6's conversion (\(e^{\pi\omega k}/|\sin\pi(\gamma+i\omega k)|\to2\)). The liminf bound is BZ.6's with the CC.4 atoms and \(\delta=\delta_{14}\) (\(k\le100\)).
*CC.6.* This is the first step of the proof of CC.3, combined with CA.3. \(\square\)

**Conjecture CA.6: PROVED by Theorem CC.** The atoms of \(\nu\) are exactly \(\Lambda_O\), and \(\nu(-Z)=0\).

**Observation (no claim).** The level sums \(s_k\) have ratio \(s_{k+1}/s_k=0.42373129\ldots\) from \(k\approx8\) on (printed for \(k=4,\dots,13\)). The geometric tail \(s_{14}r/(1-r)\), with \(r=s_{14}/s_{13}\), equals \(7.837945649\cdot10^{-7}\), against the exact remainder \(\delta_{14}=7.837945666\cdot10^{-7}\): they differ by \(2\cdot10^{-15}\). So the CB.5 masses, which are products over the orbit \(\psi(5^iu)\), add up to the \(c_R(\gamma)\) of Theorems AT/AU, which was computed by a different route, to about \(15\) digits. Neither the ratio nor the counts \(5\cdot2^{k-4}\) are proved here.

**What stays open.**
- *From CA/CB/CC.* A closed form of \(\psi\), hence of \(F\), \(G^*\), \(t_{\max}\), \(t_{\min}\) and the atoms. The exact value of \(\liminf\rho_k\). Whether \(\hat\theta_k\ne0\) for every \(k\ge22\). The ratio \(0.42373\ldots\) and the level counts (counts: now Theorem CD). \(G^*\) into BI.4/BJ.
- *Carried over.* Conjecture BS.4; a rate for \(J-A_N\); uniformity in \(\alpha\); eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem CC proves that the spectral measure of one scaled corner resolvent of a graph Laplacian on the combinatorial gasket is pure point. A spectral measure is not a flux, and nothing here is a force. Stage 1 is not touched, and no \(W\), \(0.08\), \(0.23\) or \(30\,\mu\)N/kW enters. The Stage 2 question stays OPEN: the repository has no \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh from which to decide it. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_pure_point.py`** (about 40 s; needs python-flint, mpmath, sympy; imports BZ's, CA's and CB's scripts unchanged; `python3 scripts/gasket_fractional_profile_pure_point.py K` certifies levels \(\le K\)).
- *Sympy:* \(k=1/c\) and \(F=k\,F\circ p\); the partial fractions; \(T=(2-t)/(6-t)\) from the explicit roots and from the residue theorem; \(T'\), \(T(-5)=7/11\), \(T(0)=1/3\); the regular remainder \((k(z)-k(y))/(y-z)\); \(p^{-1}(U)\) avoids \(-6,-\tfrac52\); \(-\tfrac52,-2\notin K_p\); the pole of \(\Pi\) on \(Z\) only through \(-5\); \(\psi'\)'s zeros.
- *Exact rationals:* the Cauchy–Schwarz lemma \(\operatorname{Im}(1/r)\ge\operatorname{Im}z\) on \(200\) random probability measures.
- *Interval (`arb`/`acb`, 300 bits):* \(\operatorname{Im}H-\operatorname{Im}z/2>0\), \(|F|\operatorname{Im}z<2\) and \(F=k\,F\circ p\) at \(8\) points (some within \(10^{-4}\) of \(K_p\)); \(\operatorname{Im}(1/r_n)\ge\operatorname{Im}z\) for CB.1's \(r_n\); \(\psi(-1)\), \(\psi(-5)=p(\psi(-1))\), \(f(-1)\ne0\); all \(10240\) atoms of level \(\le14\) in \([1,5)\) with CB.5 masses, the level counts, disjointness, the partial sum and \(\delta_{14}\); the largest mass against CB; \(\liminf\rho_k\le0.59486582\); \(\pi c_R(\gamma)/\ln5\); the geometric-tail consistency (to \(10^{-13}\)).
- *Floats:* mpmath Newton steps only propose the balls.

The script exits nonzero if any check fails.

## Theorem CD. The level counts of \(\nu\)'s atoms in one scale period are exact at every level: \(N_0=1\), \(N_1=0\), \(N_2=1\), \(N_3=3\) and \(N_k=5\cdot2^{k-4}\) for \(k\ge4\), so one period holds exactly \(5\cdot2^{K-3}\) atoms of level \(\le K\)

**Setting.** As in Theorem CC. \(p(z)=z(z+5)\), \(O_k=p^{-k}(-3)\), \(O=\bigcup_kO_k\) (CA); \(a=\psi(-1)\) and \(I=(p(a),a]=(\psi(-5),\psi(-1)]\). By CC.4, the atoms of \(\nu\) in \([1,5)\) are exactly the \(\lambda_t=-\psi^{-1}(t)\), \(t\in O\cap I\), and \(t\in O_k\) is the atom's level. CC certified \(N_k=\#(O_k\cap I)\) for \(k\le14\) (\(1,0,1,3\), then \(5\cdot2^{k-4}\); \(10240\) atoms) and left the counts as an observation. This note proves them for every \(k\). For an interval \(J\) write \(C_m(J)=\#(O_m\cap J)\).

**Assumptions used as stated (not re-run).** CC.4 (the parametrisation of one period by \(\psi\) on \([-5,-1]\), \(\psi\) strictly increasing there, and the level-by-level atom set), BD.4, BJ.6, and CA's enclosure routine `psi_enc` (imported unchanged). No literature assumption is used.

**THEOREM CD.**
1. **CD.1 (the preimage sets).** \(O_k\) consists of exactly \(2^k\) distinct real points, all in \((-5,0)\setminus\{-\tfrac52\}\).
2. **CD.2 (the endpoint orbit leaves \([-5,0]\) in four steps; certified).** \(-\tfrac52<a=-0.9508266438293260\ldots<0\), \(-5<p^2(a)=-4.427332832333435\ldots<p(a)=-3.850061912530890\ldots<-3\), \(-3<p^3(a)=-2.535388153409578\ldots<-\tfrac52\), and \(p^4(a)=-6.248747678598260\ldots<-5\).
3. **CD.3 (exact level counts).** For every \(k\ge0\),
   \[
   N_k=\#(O_k\cap I)=\begin{cases}1,0,1,3&k=0,1,2,3,\\ 5\cdot2^{k-4}&k\ge4,\end{cases}\qquad \sum_{k\le K}N_k=5\cdot2^{K-3}\quad(K\ge3).
   \]
   In particular one scale period holds \(10240\) atoms of level \(\le14\), as CC.4 certified, and the level-\(k\) atoms are a fraction \(5/16\) of \(O_k\) for \(k\ge4\).
4. **CD.4 (direct recount; certified consistency).** The six counting functions of the proof, recounted from `arb` enclosures of all \(2^k\) points of \(O_k\) for \(k\le18\) (\(262144\) points at the top level, every comparison with \(a\)'s orbit decided), equal the recursion.

*Proof.* *CD.1.* \(p(y)=t\) iff \(y=-\tfrac52\mp\sqrt{t+25/4}\). For \(t\in(-\tfrac{25}4,0)\) the square root lies in \((0,\tfrac52)\), so the two roots are distinct, lie in \((-5,-\tfrac52)\) and \((-\tfrac52,0)\), and neither is \(-\tfrac52\). Since \(-3\in(-\tfrac{25}4,0)\) and \((-5,0)\subset(-\tfrac{25}4,0)\), induction on \(k\) gives the claim (points of \(O_{k+1}\) with different images are different).
*CD.2.* `arb` at 300 bits from CA's enclosure of \(\psi(-1)\), with \(p\) evaluated as in BZ.
*CD.3.* \(p\) is a decreasing bijection of \([-5,-\tfrac52]\) onto \([-\tfrac{25}4,0]\) and an increasing bijection of \([-\tfrac52,0]\) onto \([-\tfrac{25}4,0]\). So for \(m\ge1\) and an interval \(J\) on one side of \(-\tfrac52\), \(C_m(J)=C_{m-1}(p(J))\), and by CD.1 only \(p(J)\cap(-5,0)\) matters; \(C_{m-1}(J')=2^{m-1}\) if \(J'\supseteq(-5,0)\) and \(0\) if \(J'\cap(-5,0)=\emptyset\). Write \(a_j=p^j(a)\). Cutting at \(-\tfrac52\) and using CD.2 for the positions:
   - \(N_k\): \((a_1,-\tfrac52)\mapsto(-\tfrac{25}4,a_2)\), \([-\tfrac52,a]\mapsto[-\tfrac{25}4,a_1]\), so \(N_k=A_{k-1}+B_{k-1}\) with \(A_m=C_m((-5,a_1])\), \(B_m=C_m((-5,a_2))\).
   - \(A_m=G_{m-1}\), \(G_j=C_j([a_2,0))\), since \((-5,a_1]\mapsto[a_2,0)\); \(B_m=F_{m-1}\), \(F_j=C_j((a_3,0))\), since \((-5,a_2)\mapsto(a_3,0)\).
   - \(G_j=E_{j-1}+2^{j-1}\), \(E_j=C_j((-5,a_3])\): \([a_2,-\tfrac52)\mapsto(-\tfrac{25}4,a_3]\) and \([-\tfrac52,0)\mapsto[-\tfrac{25}4,0)\).
   - \(E_j=2^{j-1}\): \((-5,a_3]\mapsto[a_4,0)\supseteq(-5,0)\) because \(a_4<-5\).
   - \(F_j=2^{j-1}\): \((a_3,-\tfrac52)\mapsto(-\tfrac{25}4,a_4)\), disjoint from \((-5,0)\), and \([-\tfrac52,0)\mapsto[-\tfrac{25}4,0)\).
   The level-\(0\) values follow from the position of \(-3\) (CD.2): \(N_0=G_0=E_0=1\), \(A_0=B_0=F_0=0\). Hence \(N_1=0\), \(N_2=A_1+B_1=G_0+F_0=1\), \(N_3=G_1+F_1=(E_0+1)+1=3\), \(N_4=G_2+F_2=(1+2)+2=5\), and for \(k\ge5\), \(N_k=G_{k-2}+F_{k-2}=(2^{k-4}+2^{k-3})+2^{k-3}=5\cdot2^{k-4}\). The cumulative sum is \(5+5(2^{K-3}-1)=5\cdot2^{K-3}\). No endpoint of any interval is in \(O\): \(-5,0,-\tfrac52,-\tfrac{25}4\notin O\) by CD.1, and \(a_j\in O\) would give \(a\in O\), i.e. \(\lambda=1\) an atom, which CC.4 excludes (\(f(-1)\ne0\)). So open versus closed ends do not matter (CD.4 also decides every comparison numerically). \(\square\)

**Observation CC (level counts): now PROVED by Theorem CD.** The ratio \(0.42373\ldots\) of the level sums is not touched: CD counts atoms, not masses.

**Observation (no claim).** With CD.3 the mean \(\lambda^{-\gamma}\)-mass of a level-\(k\) atom is \(s_k/(5\cdot2^{k-4})\); CC's printed level sums give a mean ratio near \(0.42373/2=0.2119\) per level. Not proved.

**What stays open.**
- *From CA/CB/CC/CD.* A closed form of \(\psi\), hence of \(F\), \(G^*\), \(t_{\max}\), \(t_{\min}\) and the atoms. The exact value of \(\liminf\rho_k\). Whether \(\hat\theta_k\ne0\) for every \(k\ge22\). The level-sum ratio \(0.42373\ldots\) (masses, not counts). \(G^*\) into BI.4/BJ.
- *Carried over.* Conjecture BS.4; a rate for \(J-A_N\); uniformity in \(\alpha\); eventual monotonicity of \(V^{(n)}\); a closed form for \(a^*\); Conjecture BL.6's limit part; closed forms for \(c_{\rm corner}\), \(W\), \(U\), \(\kappa\); Conjecture BH.7, Conjecture BF.2 and the Stage 2 question.

**NOT THIS.** Theorem CD counts the atoms of the spectral measure of one scaled corner resolvent of a graph Laplacian on the combinatorial gasket. A count of atoms is not a flux, and nothing here is a force. Stage 1 is not touched, and no \(W\), \(0.08\), \(0.23\) or \(30\,\mu\)N/kW enters. The Stage 2 question stays OPEN: the repository has no \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh from which to decide it. Functions on a graph, not a flux. Not thrust.

**Witness, computed by `scripts/gasket_fractional_profile_level_counts.py`** (about 3 s; needs python-flint, sympy; imports BZ's and CA's scripts unchanged; `python3 scripts/gasket_fractional_profile_level_counts.py K` recounts levels \(\le K\)).
- *Sympy:* the two inverse branches of \(p\) and their ranges; \(p(-\tfrac52)=-\tfrac{25}4\); \(p'=2y+5\); the induction step \((2^{k-4}+2^{k-3})+2^{k-3}=5\cdot2^{k-4}\).
- *Exact integers:* the recursion against the closed form and the cumulative count \(5\cdot2^{K-3}\) for \(k,K\le400\); agreement with CC.4's counts for \(k\le14\).
- *Interval (`arb`, 300 bits, outward):* CD.2's orbit inequalities; all \(2^k\) points of \(O_k\), \(k\le18\), separated and in \((-5,0)\setminus\{-\tfrac52\}\); the six counting sequences recounted directly.
- *Floats:* none.

The script exits nonzero if any check fails.

## Conjecture J.1

**Status: PROVED by Theorem AQ (formerly a conjecture; the text below records it as first stated).** For \(\alpha=0.45\) and these corner data, \(\|F(n)\|\) decreases for every \(n\ge 1\) and
\[
\lim_{n\to\infty}\|F(n)\| = L \quad\text{with}\quad L\ge 0.70.
\]
The only evidence is the computed sequence above: the norm is still decreasing at level \(7\), and the gap ratios \(\delta_{n+1}/\delta_n\) are increasing but still near \(0.67<1\). A geometric tail with ratio \(0.67\) would leave \(L\) near \(0.82\); a later rise of that ratio toward \(1\) can still push the limit to \(0\). Theorem AH now excludes \(0\): \(\|F(n,0.45)\|\ge 6.44\cdot10^{-4}\) at every level. Theorem AI sharpens that floor to \(\approx0.051988\); Theorem AJ sharpens it further to \(\approx0.055599\) by keeping the neighbor-edge weight \(w_{pq}\), Theorem AK to \(\approx0.058941\) by a two-moment Jensen bound on \([0,6]\), Theorem AL to \(\approx0.078019\) by replacing AA's path ceiling \(C(\alpha)\) with the Hermite–Hadamard limit \(R_{\infty}(\alpha)\), Theorem AM to \(\approx0.147106\) by replacing the path with the whole boundary cycle, Theorem AN to \(\approx0.214239\) by replacing the cycle with the boundary fence/collar (Bloch ceiling \(R^{\mathrm{fen}}_{\infty}+1/2\approx5.002104\)), Theorem AO to \(\approx0.222331\) by the depth-2 collar (period-2 Bloch ceiling \(R^{\mathrm{d2}}_{\infty}+11/20\approx4.820047\)), and Theorem AP to \(\approx0.274743\) (AM route) or \(\approx0.415238\) (AO as stated) through the exact share identity \(\sigma_n=(L^{1-\alpha})_{aa}+\omega_n\) and a five-moment Radau floor \(\sigma_n\ge1.310267\), still below \(0.70\) (Kept Failures AI.1, AJ.1, AK.1, AL.1, AM.1, AN.1, AO.1, AP.1; AK.1 closes the one-step neighbor-Harnack family, AL.1 closes the sharp-path-ceiling route, AM.1–AO.1 show no resistance ceiling can lift the \(\sigma_{AK}\) chain to \(0.70\), and AP.1 shows no share floor can lift the proved fractional floors to \(0.70\); with \(S_{AP}\) the remaining open step is a resistance ceiling \(R_{\alpha}(n)\le2.859237\)). Neither proves that the limit exists, that the norm decreases, or that \(L\ge 0.70\); those remained conjecture until Theorem AQ. **Theorem AQ settles all three.** The exact decimation recursion \(r_{n+1}(s)=\Phi_s(r_n(s(s+5)))\) for the corner resolvent \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\) makes \(r_n\) strictly increasing in \(n\). By the Stieltjes formulas \(R_{\alpha}(n)=\tfrac{2\sin\pi\alpha}{\pi}\int s^{-\alpha}r_n\) and \(\sigma_n=\tfrac{\sin\pi\beta}{\pi}\int s^{\beta-1}(1-sr_n)\), \(R_{\alpha}(n)\) increases and \(\sigma_n\) decreases, so \(\|F(n)\|=\sigma_n\sqrt{21}/(3R_{\alpha}(n))\) strictly decreases for every \(n\ge1\). An interval-arithmetic certificate bounds \(R_{0.45}(n)<2.518955\) and \(\sigma_n>1.350986\) uniformly. Hence \(L=\lim\|F(n)\|\) exists and \(L\in[0.819255,0.822959]\), so \(L\ge0.70\). This also closes the step Kept Failure AP.1 left open: \(\bar R=2.518955<2.859237\).

**Falsifier.** The first level \(n\) with \(\|F(n)\|<0.70\), or the first \(n\ge 7\) with \(\|F(n+1)\|>\|F(n)\|+10^{-8}\). Under Theorem AQ neither can fire: the norm strictly decreases and stays above \(0.819255\). A computed violation would indicate an error in the build, the solver, or the certificate. Level \(8\) is the first level not required by the default witness. Through level \(10\) neither falsifier has fired: \(\|F(10)\|\approx 0.833640\), still decreasing.

**Status of J.2 below: PROVED by Theorem AS (5)** (\(\|F(n,0.10)\|\) strictly decreases to \(\mathcal L(0.10)\in[1.410638933063,1.410638933065]\), so \(\|F(n,0.10)\|>1.410\) at every level; the text records it as first stated). A second, separate numerical guess, also not a theorem: at \(\alpha=0.10\), levels \(1\) through \(6\) give
\[
\|F\|\approx(1.50450860,\; 1.44384308,\; 1.42364897,\; 1.41577938,\; 1.41266163,\; 1.41143248),
\]
and the successive drops after level \(3\) contract by a factor \(<1/2\). **CONJECTURE J.2.** \(\|F(n,0.10)\|\ge 1.410\) for every \(n\). Falsified by the first level with \(\|F\|<1.410\). The \(\alpha\to 0\) profile \(u\equiv (a+b+c)/3\) off the corners has audit norm \(\sqrt{7/3}\approx 1.527525\), independent of \(n\ge 1\); \(\alpha=0.10\) is a different operator, and this sentence does not identify its limit with \(\sqrt{7/3}\).


## What this does not say

- **LIMIT AT \(\alpha=0.45\): SETTLED BY THEOREM AQ (the history below predates it).** Theorem AQ proves that \(\|F(n,0.45)\|\) strictly decreases to \(\mathcal L\in[0.819255,0.822959]\) by exact spectral decimation of the corner resolvent plus an interval certificate. Kept Failure AQ.1 shows that no positive limit exists for \(\alpha\ge\log3/\log5\). Before AQ: levels \(1\) through \(7\) still have increasing ratios, now past \(0.97\). Theorem J excludes \(+\infty\). Theorem K gives \(\|F\|\to 0\) only for \(\alpha>\log 3/\log 5\). Theorem L says energy collapse would force \(\|F\|\to 0\) at every \(\alpha\in(0,1)\), including \(0.45\), but does not prove energy collapse. Neither decides zero versus positive at \(\alpha=0.45\); Theorem AH does (not zero). Conjecture J.1 is not a theorem. Kept Failures J.1, J.2, J.3, and J.4 reject four sufficient conditions that would have closed it. Level \(10\) is recorded, not a limit. Theorem S settles the harmonic audit under uniform multiplicative weights (finite nonzero only at resistance \(\lambda=5/3\)); Kept Failure S.1 rejects inheriting that freeze for \(L^{0.45}\). Theorem W settles the multi-scale harmonic audit (finite for \(\lambda<5/3\); resistance diverges); Kept Failure W.1 rejects a multi-scale resistance freeze. Theorem X takes spectral \(L^{\alpha}\) on the multi-scale graph: neutrality and Schur hold, nested harmonicity is special to \(\alpha=1\); Kept Failures X.1–X.2 reject inheriting the \(\lambda=1\) freeze and reject coincidence with finest-edge fractional. Theorem AA gives the level-uniform floors \(Q_{0.45}\ge 0.128335\) and \(\|F_{\mathrm{frac}}\|\ge 0.084015\) (a floor exists at every \(\alpha<1/2\)); Kept Failure AA.1 closes the energy-collapse route to zero. Theorem AB factors \(\sigma_n=\|F\|/\|F_{\mathrm{frac}}\|=\rho_n/w_n\) and, with AA, equates \(\|F\|\to 0\) with \(\rho_n\to 0\); Kept Failure AB.1 rejects geometric collapse of \(\rho\). Theorem AC extends the energy floor to every \(\alpha<\log 3/\log 5\) by cell resistance and spectral mass, so AB's equivalence holds on that whole range; Kept Failure AC.1 (via Theorem K) rejects a floor above it. Theorem AD proves the exact pairing \(Q_{\alpha}R_{\alpha}=7/3\) and settles the endpoint: \(Q_{\gamma}(u_n)=\Theta(1/n)\) and \(\|F(n,\gamma)\|\to 0\); Kept Failure AD.1 rejects extending the AC floor to \(\gamma\). Theorem AE factors \(\rho_n=\lambda_n/(\lambda_n+(1-\lambda_n)\gamma_n)\) with sharp ceiling \(w_n\le8^{\alpha}/2\) and the resistance identity \(\rho_n=(\sqrt{21}/7)w_n\|F\|R_{\alpha}\); Kept Failure AE.1 rejects geometric wipeout of the neighbor-degree share \(\lambda_n\). Theorem AF floors that share, \(\lambda_n\ge\alpha\,2^{-1-\alpha}\) (Jensen: \(\ge\alpha\,2^{\alpha-1}\)), so for \(\alpha<\log3/\log5\) vanishing of \(\|F\|\) requires \(\gamma_n\to\infty\) (one-sided); Kept Failure AF.1 rejects geometric blowup of \(\gamma_n\). Theorem AG ceilings the same share, \(\lambda_n\le 1-\alpha(1-\alpha)/2^{3-\alpha}<1\), so \(\|F\|\to0\) if and only if \(\gamma_n\to\infty\) (two-sided); Kept Failure AG.1 rejects geometric climb of \(\lambda\) toward \(1\). Theorem AH compares the jump kernels of \(p\) and \(a\) (\(w_{pj}\ge\kappa(\alpha)w_{aj}\) off the corner patch) and applies one-step \(L^{\alpha}\)-harmonicity at the two neighbors, so \(\gamma_n\le\Gamma_{AH}(\alpha)<\infty\) and \(\|F(n,\alpha)\|\ge\sigma_{\min}(\alpha)\|F_{\mathrm{frac}}\|>0\) at every level for every \(\alpha<\log3/\log5\) (\(\ge 6.44\cdot10^{-4}\) at \(0.45\)); Kept Failure AH.1 rejects the max-principle shortcut to the same floor. Theorem AI replaces \(\kappa\) by \(1\) through the first-step identity at the degree-\(2\) corner (\(w_{pj}+w_{qj}\ge2w_{aj}\) on Far), so \(I_{\mathrm{far}}\le(4^{\alpha}/2)I_a\), \(\gamma_n\le8/(\alpha(1-\alpha))\), and \(\|F(n,0.45)\|\ge0.051988\) at every level; Kept Failure AI.1 rejects reading that as Conjecture J.1's \(0.70\). Theorem AJ keeps the neighbor-edge weight \(w_{pq}\ge\alpha4^{\alpha-1}\) in the Harnack identity, so \(I_{\mathrm{far}}\le H_{AJ}I_a\) with \(H_{AJ}=(4^{\alpha}-\alpha4^{\alpha-1})/2\), \(\gamma_n\le2(4-\alpha)/(\alpha(1-\alpha))\), and \(\|F(n,0.45)\|\ge0.055599\); Kept Failure AJ.1 rejects reading that as \(0.70\). Theorem AK uses the second moment and the sharp line-graph ceiling \(\mathrm{spec}(L)\subset[0,6]\) to lower the diagonal ceilings to \((L^{\alpha})_{aa}\le6^{\alpha}/9+(8/9)(3/2)^{\alpha}\) and \((L^{\alpha})_{pp}\le(6^{\alpha}+2^{\alpha})/2\), so \(\|F(n,0.45)\|\ge0.058941\); Kept Failure AK.1 shows that the one-step neighbor-Harnack chain stays below \(0.70\) from level \(4\) even with every constant live. Theorem AL replaces AA's path ceiling \(C(\alpha)\) by the Hermite–Hadamard limit \(R_{\infty}(\alpha)\) for \(\alpha\ge\alpha_*=(14-8\sqrt{2})/34\), so \(\|F_{\mathrm{frac}}(n,0.45)\|\ge0.111208\) and \(\|F(n,0.45)\|\ge0.078019\); Kept Failure AL.1 shows that this sharp path ceiling still cannot reach \(0.70\). Theorem AM replaces the path by the full boundary cycle \(C_{3\cdot2^n}\), whose fractional corner resistance is below \((1-\alpha)R_{\infty}(\alpha)\) for every \(\alpha<1/2\) (Fourier plus a heat-kernel periodization), and uses all three corners, so \(Q_{0.45}\ge0.320299\), \(\|F_{\mathrm{frac}}(n,0.45)\|\ge0.209685\), and \(\|F(n,0.45)\|\ge0.147106\); Kept Failure AM.1 shows this still cannot reach \(0.70\), and that even the live energy times \(\sigma_{AK}\) stays below \(0.70\) from level \(2\). Theorem AN replaces the cycle by the boundary fence/collar with Bloch ceiling \(R^{\mathrm{fen}}_{\infty}+1/2\approx5.002104\), so \(Q\ge0.466470\), \(\|F_{\mathrm{frac}}\|\ge0.305377\), and \(\|F(n,0.45)\|\ge0.214239\); Theorem AO replaces that by the depth-2 collar with period-2 Bloch ceiling \(R^{\mathrm{d2}}_{\infty}+11/20\approx4.820047\), so \(Q\ge0.484089\), \(\|F_{\mathrm{frac}}\|\ge0.316911\), and \(\|F(n,0.45)\|\ge0.222331\); Kept Failures AN.1 and AO.1 show those routes stay below \(0.70\). Theorem AP proves the exact share identity \(\sigma_n=(L^{1-\alpha})_{aa}+\omega_n\) (\(Lu=L^{1-\alpha}L^{\alpha}u\) and AD's corner currents) and a five-moment Gauss–Radau floor \(\sigma_n\ge S_{AP}(0.45)\approx1.310267>1\), beating AK's \(0.701556\) and AK.1's \(0.9529\) cap, so \(\|F(n,0.45)\|\ge0.274743\) with AM's floor and \(\ge0.415238\) with AO's floor as stated; Kept Failure AP.1 shows no share floor can lift the proved fractional floors to \(0.70\), so the open step is now a resistance ceiling \(R_{\alpha}(n)\le2.859237\). Neither floor is a limit theorem. Theorem AQ is one: \(\|F(n,0.45)\|\downarrow\mathcal L\in[0.819255,0.822959]\) (Conjecture J.1 proved). Theorem AR adds the exact rates: \(R_{\alpha}(n+1)-R_{\alpha}(n)\sim c_R(\alpha)(5^{\alpha}/3)^n\) for every \(\alpha\in(0,1)\), so the gap ratio at \(0.45\) tends to \(5^{0.45}/3\), \(R_{\gamma}(n)\) has an exact linear slope, and \(\|F(n,\alpha)\|=\Theta((3\cdot5^{-\alpha})^n)\) above \(\gamma\); Kept Failure AR.1 shows the universal resolvent bound alone gives no rate. Theorem AS certifies the digits: \(r_\infty\) and \(1-sr_\infty\) are completely monotone, so polynomial interpolants on a geometric grid are rigorous floors and ceilings, and outward-rounded sums give \(\mathcal L(0.45)\in[0.82113862713346,0.82113862713352]\) and \(\mathcal L(0.10)>1.410638933063\), which proves Conjecture J.2; Kept Failure AS.1 shows that level-\(10\) extrapolation with AR's exact ratio still misses \(\mathcal L\) by \(2.8\cdot10^{-4}\). Theorem AT certifies AR's rate constants: an exact change of variables moves \(\rho^{-m}\Delta R_m\) to level \(n\) with weights trapped in a bracket of relative width \(O(5^{-0.65n})\), so at \(n=40\) \(c_R(0.45)\in[0.497802650763560354,0.497802650763560356]\), \(c_S(0.45)\in[0.927092638569461377,0.927092638569461384]\), \(c_R(\gamma)\in[0.517165834630505243,0.517165834630505245]\), and \(\|F(n,0.45)\|-\mathcal L\sim c_F\rho^n\) with \(c_F\in[0.52027873061815,0.52027873061821]\); Kept Failure AT.1 shows \(s^{1-\gamma}r_\infty(s)\) has no limit (certified log-periodic oscillation), so the slope at \(\gamma\) is not a pure-power constant. Theorem AU proves Conjecture AT.2: the slope is the log-period mean, \(c_R(\gamma)=\frac{2\sin\pi\gamma}{\pi}\ln5\,\langle\Theta\rangle\), certified to \(34\) digits inside AT's enclosure; the Mellin transform of \(r_\infty\) has simple poles at \(1-\gamma+2\pi ik/\ln5\), certified genuine for \(|k|\le4\); Kept Failure AU.1 rejects a single-harmonic profile. Theorem AY certifies the intercept: at \(\gamma\) the bracket weights of AT do not depend on the level, so the per-level defects \(|\Delta_j-c|\) sum to at most \(\int t^{-\gamma}(r_\infty-r_n)W_\pm\), and \(C_\gamma\in[0.70307046242185742992904816,\ldots817]\), \(\nu_\gamma=C_\gamma/c_R(\gamma)=1.3594681151437\ldots\), \(\sigma_\infty(\gamma)=1.1519186116740\ldots\), so \((n+\nu_\gamma)\|F(n,\gamma)\|=K_\gamma+O(5^{-0.65n})\) with \(K_\gamma=3.4023607638326\ldots\); Kept Failure AY.1 rejects the unshifted renormalization \(n\|F(n,\gamma)\|\), whose error is \(\Theta(1/n)\). Theorem AZ settles the supercritical side: for \(\alpha>\gamma\), \(R_\alpha(n)=A_\alpha\rho^n+C_\alpha+O(\rho^n5^{-0.65n})\) with the intercept \(C_\alpha=\frac{2\sin\pi\alpha}{\pi}M(1-\alpha)\), the continuation of \(R_\infty(\alpha)\) across its simple pole at \(\gamma\), so \((\rho^n+\nu_\alpha)\|F(n,\alpha)\|=K_\alpha+O(5^{-0.65n})\) and \(\|F(n,\alpha)\|\sim K_\alpha(3\cdot5^{-\alpha})^n\) (the square of Theorem K's ratio, no log-periodic modulation); certified \(C_{0.8}=-1.74325717839207572395\ldots\), \(\nu_{0.8}=-0.72744538444215\ldots\), \(K_{0.8}=0.68786675841399\ldots\); Kept Failure AZ.1 rejects carrying AY's shift sign over: \(\rho^n\|F(n,0.8)\|\) decreases to \(K_{0.8}\) from above with error \(\Theta((3\cdot5^{-0.8})^n)\).
- **NOT A RESIDUAL.** Their tilt diagnostic sets \((Lu)_i=1\) on the interior. Theorem F then says the corner currents sum to \(-(N-3)\), the source that was inserted. That sum grows like \(3^n\). It is not a force left after the source is removed.
- **FAILED as a continuum thrust for the combinatorial harmonic case.** Under Theorem G, \(\|F(n)\|\to 0\) for every fixed corner triple. The combinatorial harmonic limit does not keep a finite audit vector.
- **SETTLED (Theorem AQ).** \(\|F(n)\|\) for \(L^{0.45}\) strictly decreases and converges to \(\mathcal L\in[0.819255,0.822959]\). Its closed-form value is not known. History: it does not tend to \(0\) (Theorem AH: \(\|F(n,0.45)\|\ge 6.44\cdot10^{-4}\) at every level; Theorem AI sharpens this to \(\ge0.051988\); Theorem AJ to \(\ge0.055599\); Theorems AK, AL, AM, AN, AO to \(\ge0.058941\), \(\ge0.078019\), \(\ge0.147106\), \(\ge0.214239\), \(\ge0.222331\); Theorem AP to \(\ge0.274743\) with AM's floor and \(\ge0.415238\) with AO's floor as stated). Divergence is excluded (Theorem J). Neutrality is settled (Theorem H). For \(\alpha>\log 3/\log 5\) the same audit norm tends to \(0\) (Theorem K). Energy controls the audit norm at every \(\alpha\in(0,1)\) (Theorem L), so energy collapse would force \(\|F\|\to 0\), and a positive audit-norm liminf would force a positive energy liminf; the energy does not vanish at \(\alpha=0.45\) (Theorem AA: \(Q_{0.45}(u_n)\ge 0.128335\) and \(\|F_{\mathrm{frac}}\|\ge 0.084015\) at every level, for the same reason at every \(\alpha<1/2\)), so energy collapse cannot send \(\|F\|\) to \(0\) (Kept Failure AA.1); Theorem AB equates vanishing of \(\|F\|\) with vanishing of the two-neighbor fractional share \(\rho_n\) (Kept Failure AB.1 rejects geometric collapse of \(\rho\)); Theorem AC proves the energy floor for every \(\alpha<\log 3/\log 5\), the exact threshold by Theorem K (Kept Failure AC.1), and Theorem AD settles \(\alpha=\log 3/\log 5\) itself: \(Q_{\gamma}(u_n)=\Theta(1/n)\to 0\) and \(\|F(n,\gamma)\|\to 0\) (Kept Failure AD.1); Theorem AE splits \(\rho_n\to0\) into degree-share and gap-ratio blowup and rejects geometric \(\lambda\) wipeout (Kept Failure AE.1); Theorem AF floors \(\lambda_n\) uniformly, so \(\|F\|\to0\) would force the far-to-neighbor gap ratio \(\gamma_n\to\infty\), and Kept Failure AF.1 rejects geometric \(\gamma\) blowup (live \(\gamma_n\) decreases through level \(6\)); Theorem AG ceilings \(\lambda_n\) uniformly below \(1\), so \(\|F\|\to0\) if and only if \(\gamma_n\to\infty\), and Kept Failure AG.1 rejects geometric climb of \(\lambda\) toward \(1\); Theorem AH ceilings \(\gamma_n\le\Gamma_{AH}(\alpha)\) by a kernel comparison and one-step Harnack at the neighbors, so \(\|F\|\) does not vanish for any \(\alpha<\log3/\log5\), and Kept Failure AH.1 rejects the max-principle shortcut; Theorem AI sharpens the Harnack constant from \(4^{\alpha}/(2\kappa)\) to \(4^{\alpha}/2\) by a two-neighbor heat-kernel sum, giving \(\gamma_n\le8/(\alpha(1-\alpha))\) and \(\|F(n,0.45)\|\ge0.051988\), and Kept Failure AI.1 rejects reaching \(0.70\) that way; Theorem AJ keeps \(w_{pq}\) to push the floor to \(0.055599\), and Kept Failure AJ.1 rejects reaching \(0.70\) that way; Theorem AK sharpens the diagonal ceilings by two-moment Jensen on \([0,6]\) to push it to \(0.058941\), and Kept Failure AK.1 shows the whole one-step neighbor-Harnack chain stays below \(0.70\) from level \(4\) even with live constants; Theorem AL sharpens the path-resistance ceiling to \(R_{\infty}\) for \(\alpha\ge\alpha_*\) to push the floor to \(0.078019\), and Kept Failure AL.1 shows that route also stays below \(0.70\); Theorem AM replaces the path by the boundary cycle to push the floor to \(0.147106\), Theorem AN replaces the cycle by the boundary fence/collar (Bloch ceiling \(R^{\mathrm{fen}}_{\infty}+1/2\)) to push it to \(0.214239\), and Theorem AO replaces that by the depth-2 collar (period-2 Bloch ceiling \(R^{\mathrm{d2}}_{\infty}+11/20\)) to push it to \(0.222331\); Kept Failures AM.1, AN.1, and AO.1 show those routes stay below \(0.70\), and that no resistance ceiling can rescue the \(\sigma_{AK}\) chain; Theorem AP replaces that chain by the exact identity \(\sigma_n=(L^{1-\alpha})_{aa}+\omega_n\) and the moment floor \(\sigma_n\ge1.310267\), and Kept Failure AP.1 shows no share floor can lift the proved fractional floors to \(0.70\), leaving a resistance ceiling \(R_{\alpha}(n)\le2.859237\) as the open step. Theorem AQ supplies it, \(R_{\alpha}(n)<2.518955\), by exact spectral decimation of the corner resolvent, and proves monotone convergence of \(\|F(n,0.45)\|\) to \(\mathcal L\in[0.819255,0.822959]\). Kept Failure AQ.1 shows that the uniform ceiling fails for \(\alpha\ge\log3/\log5\). Geometric weights \(1/d^2\) on this build are settled (Theorem I): that audit norm diverges. Uniform multiplicative weights are settled for the harmonic case (Theorem S): finite nonzero only at the resistance value \(\lambda=5/3\); resistance weights do not freeze the fractional \(\alpha=0.45\) audit (Kept Failure S.1). Multi-scale multiplicative weights are settled for the harmonic case (Theorem W): finite nonzero for every \(\lambda<5/3\), and resistance diverges like \((n+1)\sqrt{21}/2\) (Kept Failure W.1). Spectral \(L^{\alpha}\) on that multi-scale graph is settled structurally by Theorem X (neutrality, Schur, \(\alpha=1\) recovers W, no nested fractional harmonicity); the \(\lambda=1\) freeze does not transfer to \(\alpha=0.45\) (Kept Failure X.1), and multi-scale fractional does not coincide with finest-edge fractional (Kept Failure X.2). None of this is thrust.
- **OPEN.** The Stage 2 question from the rectangle note, whether an informational tensor on the \(0.45\) mesh has nonzero integrated divergence after the Maxwell piece is removed. Theorem F says the harmonic combinatorial gasket does not supply that divergence. Under the Hessian reading of the rank-2 symbol, Theorem M in the README reduces the question to \(\Delta\Psi_{\mathrm{info}}\) on the enclosing surface (zero for static massless \(A_0\) in vacuum, Kept Failure M.1). Under the quadratic reading, Theorem N makes it a volume integral of \(\partial_i\Psi\,\Delta\Psi\), which for static massless \(A_0\) is the Maxwell piece itself and vanishes for an isolated device (Kept Failure N.1). Theorem O extends this to every two-derivative reading: the only survivor is a surface-dependent trace flux whose all-space total is zero (Kept Failure O.1). Theorem P extends it to every linear reading of any order: zero for shift-invariant ones (Kept Failure P.1), and a surface-dependent \(-\int_V E\) for the undifferentiated term (Kept Failure P.2). Theorem Q extends it to every shift-invariant polynomial reading of any order and degree: zero on every vacuum surface if the divergence vanishes on harmonic fields, otherwise surface-dependent (Kept Failure Q.1). Theorem R drops shift invariance and keeps the same dichotomy; the undifferentiated \(\delta^{ij}\Psi^2\) term does not even decay under dilation (Kept Failure R.1). Theorem T takes the same polynomial family under screening and still finds only surface-dependent mass remainders (Kept Failure T.1). Theorem U covers every first-derivative \(O(3)\)-covariant reading \(f(s)\delta+g(s)E\otimes E\) with \(C^1\) coefficients (non-polynomial allowed): conserved ones are Maxwell plus \(\kappa\delta\) and give \(\mathcal{G}=0\); all others are surface-dependent and decay like \(R^{-2}\) (Kept Failure U.1). Theorem V covers the \(\Psi+\nabla\Psi\) isotropic extension \(f(\Psi,s)\delta+g(\Psi,s)E\otimes E\): conserved ones remain Maxwell plus \(\kappa\delta\) with \(\mathcal{G}=0\); all others are surface-dependent (Kept Failure V.1). Theorem Y allows explicit position dependence and keeps the same conserved class (Kept Failure Y.1). Theorem Z allows the pseudotensor \(h\,\varepsilon\cdot\nabla\Psi\): the only new conserved piece is a curl with \(\mathcal{G}=0\) on every closed surface (Kept Failure Z.1). Higher-derivative non-polynomial readings remain OPEN. The mesh solve itself is taken up in Theorem BA of the README under one stated reading (Maxwell form \(f_c=q_cE_c\) with \(q=L^{0.45}u\) on the corners and the corner-cell gradient): the action–reaction piece integrates to exactly \(0\), and the residual \(D=\kappa_n(\alpha)\sum_c(d_c-\bar d)^2(P_c-\bar P)\) has certified \(\kappa_n(0.45)<0\) for \(n=1,\ldots,4\) (\(\kappa_4=-18.526339511158477727\ldots\)), is the same on every vertex set containing the corners, and is a discrete Newton's-third-law defect that grows like \(2^n\), not a residual force (Kept Failure BA.1); at \(\alpha=1\), \(\kappa_n(1)=-9((18/25)^n-(6/25)^n)\to0\). Theorem BB of the README decimates that constant with this note's corner resolvent: \(\kappa_n=2^n(4/R_\alpha(n)^2)(3Y_n-\sigma_n)\), with \(y_n(s)=r_n(s)\prod_{j<n}(5+p^j(s))^{-1}\), so \(0<Y_n\le5^{-n}R_\alpha(n)/2\). With Theorems AQ–AS this proves Conjecture BA.2: \(2^{-n}\kappa_n(0.45)\to-4\sigma_\infty/R_\infty^2\in[-0.854646905120924,-0.854646905120848]\), and the limit is \(0\) for \(\alpha\ge\log3/\log5\). The non-self (cross) part tends to \(0\), so the surviving piece is self-force; it is still certified nonzero at levels \(1\)–\(4\) (Kept Failure BB.1). Theorem BC of the README tests the centred (action–reaction) stencil on the gasket itself. The \(S_3\) reflection law \(H_{kc}=R_{ck}H_{ck}\) makes the Maxwell piece exactly the central part of the pair force, and \(\kappa_n=-(4/R_n^2)(s+2\sqrt3h_\perp)\). No vertex of \(V_n\) admits an odd stencil, because the isometry group is \(D_3\), with no half-turn. The Green–Gauss centred rule equals the one-sided rule at the charged corners, so \(\kappa_n\) is unchanged. The certified non-centrality is \(\tan\theta_n=(E_n-3Y_n)/(\sqrt3(E_n+Y_n))<0\) at levels \(1\)–\(4\) (Kept Failure BC.1). Conjecture BC.2: \(E_n/Y_n\to2\), so \(\tan\theta_n\to-1/(3\sqrt3)\). Theorem BD proves it. By induction on BB's recursions, \(e_n(s)=(2-t_n(s))y_n(s)\) exactly, where \(t_n(s)=s[(L_n+s)^{-1}\mathbf 1_C]_c\) is an increasing, concave Stieltjes function with \(t_n(0^+)=2/(3^n+1)\) and \(t_n'(0)=3^n(5^n-1)/(6(3^n+1)^2)\). So \(2/(3^n+1)<2-E_n/Y_n\le2/(3^n+1)+B_\alpha3^{-n}\), and \(\kappa^{\mathrm{cross}}_n>0\) at every level and every \(\alpha\). The scaling limits \(3^nt_n(5^{-n}\tau)\) and \(3^ny_n(5^{-n}\tau)\) exist, so \(c(\alpha)=\lim3^n(2-E_n/Y_n)\) exists: certified \(c(0.45)\in[3.72198897,3.72914181]\) and \(c(0.8)\in[2.57126273,2.57352573]\) (Kept Failure BD.1: \(c\neq2\), and \(c\) depends on \(\alpha\)). Theorem BE of the README decides Conjecture BD.2 in part. A symmetrized exponential-tilt kernel makes \(2-E_n/Y_n\) strictly decreasing in \(\alpha\) at every level, so \(c(\alpha)\) is strictly decreasing on \((0,1)\). Also \(c(1^-)=2\), with \((c(\alpha)-2)/(1-\alpha)\to3S\in[2.70632481,2.70903118]\) (certified). But \(c(0^+)=\int T\mathcal Y/\int\mathcal Y\) is finite, certified in \([5.48553067,5.50249013]\), so BD.2's guess \(c(0^+)=\infty\) is false (Kept Failure BE.1). Theorem BF of the README proves Conjecture BE.2. Write \(c(1-\varepsilon)-2=\varepsilon N/(\tfrac13+\varepsilon G)\), with \(N,G\) shifted exponential moments of \((T-2)\mathcal Y\) and \(\mathcal Y-\tfrac13\mathbf 1_{v<0}\). These are smooth on \([0,1]\), and an interval certificate on \(40\) \(\varepsilon\)-intervals gives \(c''\ge0.52046\) on all of \([0,1]\). So \(c\) is strictly convex, \(c'(0^+)\in[-4.341285,-4.227515]\), and \(\max(2+3S(1-\alpha),c_0+c'(0^+)\alpha)<c(\alpha)<2+(c_0-2)(1-\alpha)\). OPEN: closed forms for \(c(\alpha)\), \(S\), \(c(0^+)\) and the limit functions, and the non-monotone shape of \(c''\) (Conjecture BF.2, floats only). Theorem BG (this note) takes the refinement limit of the \(L^{\alpha}\)-Dirichlet extension itself. Exactly, \(u_n(x)=\bar d+\int s^{-\alpha}r_nh_m(p^{n-m}s;x)/\int s^{-\alpha}r_n\), and \(|u_n-\bar d|\le\frac83 4^{-\alpha}5^{m(1-\alpha)}\|x\|\rho^n/R_\alpha(n)\) on \(V_m\). So \(u_n\to\bar d\) at every non-corner vertex of \(V_*\) for \(\alpha\le\gamma\): the limit is flat, with the data left on the corners. At \(0.45\), \(\rho^{-n}(u_n(m_{ab})-\bar d)\to-x_cU\) with \(U\in[0.143810693,0.143889868]\) (certified). Above \(\gamma\) the limit is not flat: \(\lim u_n(m_{ab})=\bar d-\kappa_{0.8}x_c\), \(\kappa_{0.8}\in[0.086072435,0.086089770]\). The flat limit is not uniform: the corner-neighbour layer keeps a gap \(5/6-5\mathcal L/(2\sqrt{21})=0.3853654863\ldots\) (Kept Failure BG.1). Uniform convergence and continuity above \(\gamma\) remain Conjecture BG.2. Theorem BH (this note) resolves that layer. In layer coordinates \(\xi\) (vertex \(2^{-n}\xi\)), \(u_n(2^{-n}\xi)\to\bar d+x_a\Theta(\xi)\) for \(\alpha<\gamma\), with a universal profile: \(x_b,x_c\) drop out, and the error is \(O((3/5)^n+(J_r-D_n)/D_n)\). The neighbour profile is the \(\alpha\)-free series \(\eta(z)=\sum_j\alpha_1(p^jz)\prod_{i<j}\beta(p^iz)\), and \(r_\infty=1/(2+s-2\eta)\). The two neighbours differ exactly by \((x_b-x_c)I_n/D_n\le|x_b-x_c|5^{-n}\). Certified at \(0.45\): \(\Theta(2,0)\in[0.340381375,0.340407195]\), \(\Theta(1,1)\in[0.382431398,0.382456436]\), \(\Theta(4,0)\in[0.241594653,0.241622693]\). \(\theta_k=\Theta(2^k,0)\) decreases to \(0\) like \(\theta_\infty\rho^k\) (derived from AT). For \(\alpha\ge\gamma\) there is no layer: \(u_n(2^{-n}\xi)\to d_a\). Theorem BI (this note) proves Conjecture BG.2. A scale-by-scale oscillation recursion for the resolvent, \(\omega_{i+1}\le\tfrac35\omega_i+\frac{z_i}{z_i+2}r_n\|x\|\) with \(z_i=p^{n-i-1}(s)\), together with the certified ceiling \(\mathcal R_n(\tau)\le\min(\tfrac13,0.6031\,\tau^{\gamma-1})\) and the floor \(\mathcal R_n\ge1/(3+\tau)\), gives \(\operatorname{osc}u_n\le\|x\|K_\alpha\rho^{-N}\) on every level-\(N\) cell, uniformly in \(n\), for \(\alpha>\gamma\) (\(K_{0.8}\le47.109\)). So \(u_n\) converges uniformly to a continuous \(\bar u\) on the gasket with \(\bar u=d\) on \(C\). Sharpness of \(\rho^{-N}\) (floats \(c_{\rm corner}(0.8)\approx1.1359\)) is Conjecture BI.7. Theorem BJ (this note) proves it. The exact identity \(1/r_n-1/r_\infty=\prod_{i<n}\beta(p^is)(1+2\eta(p^ns))\le3(3/5)^n\) traps the rescaled resolvent between \(1/(\tau^{1-\gamma}/g+3)\) and \(\tau^{\gamma-1}g\), with \(g=s^{1-\gamma}r_\infty\in[0.600945014,0.600963487]\) near \(0\) (certified). So the corner deficit satisfies \(\rho^N\delta_N\ge c_-(\alpha)>0\) for every \(N\) (\(c_-(0.8)\ge0.6208\)): \(\rho^{-N}\) is the exact order of the cell modulus, and the Hölder exponent \((\alpha-\gamma)/(1-\gamma)\) cannot be raised. At \(0.8\), \(\liminf\) and \(\limsup\) of \(\rho^N\delta_N\) lie in \([1.135890346,1.135938783]\). Theorem BK (this note) proves Conjecture BJ.7 with exponent \(c=1\). An exact change of variables carries level \(n+1\) back to level \(n\): \(u_{n+1}(x)-\bar d=\langle\omega_nF_n\rangle/\langle\omega_n\rangle\) under the level-\(n\) measure. The weight is \(\omega_n=(1+y/5)^{\alpha}/((1+2y/5)\tfrac53(\beta-er_n(p)))\ge1\), because \((1+2y/5)\tfrac53\beta=1-y(y+4)/(3(y+2)(y+5))\), and its mean is \(D'_{n+1}/D'_n\), with \(D'_n=\rho^{-n}D_n\) increasing. So \(|u_{n+1}-u_n|\le2\|x\|\log(D'_{n+1}/D'_n)\) at every vertex, and \(\sup_{V_n}|u_n-\bar u|\le2\|x\|\log(D'_\infty/D'_n)\le C_\alpha\|x\|\rho^{-n}\). The constants are certified: \(J_{0.8}\le1.188154\), \(C_{0.8}\le6.2192\), \(C_{0.75}\le12.5476\), \(C_{0.9}\le1.9725\), and the asymptotic constant at \(0.8\) is \(\le2.155275\) with AZ's \(A_{0.8}\) (\(2|\nu_{0.8}|=1.45489\), derived from AZ's \(C_{0.8}\)). Floats put the corner-neighbour error at \(0.6025\,\rho^{-n}\), so whether \(\rho^{-n}\) is attained stays Conjecture BK.6. Theorem BL proves that the order \(\rho^{-n}\) is attained. Only the deficit integrand needs to be monotone, and it is: \(1-\eta_N(z)=z/(z+2)+\beta(1-\eta_{N-1}(pz))\) makes \(\eta_N\) strictly decreasing. With \(\omega_k\ge\omega_{\rm lo}\ge1\), the one-level covariance splits at \(y_0=1/10\) into a floor of order \(\rho^{-k}\) and a loss bounded by BK.4's ceiling. At \(0.8\), for every \(N\ge9\), \(\delta_N^{(n)}\) strictly increases in \(n\) and \(c_N\rho^{-n}\le\delta_N-\delta_N^{(n)}\le4.6644\rho^{-n}\), with \(c_N\ge0.1240\) for \(N\ge40\) (certified, AZ-free; floats \(0.6025\)). BK.6 fails as worded at \(N=0\), where \(\delta_0^{(n)}\equiv\tfrac32\) (Kept Failure BL.0). The existence of the limit constant is Conjecture BL.6.

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
python3 scripts/gasket_fractional_moment_jensen.py
python3 scripts/gasket_fractional_path_ceiling.py
python3 scripts/gasket_fractional_cycle_ceiling.py
python3 scripts/gasket_fractional_fence_ceiling.py
python3 scripts/gasket_fractional_depth2_ceiling.py
python3 scripts/gasket_fractional_moment_share.py
python3 scripts/gasket_fractional_decimation.py
python3 scripts/gasket_fractional_limit_digits.py
python3 scripts/gasket_fractional_rate_constants.py
python3 scripts/gasket_fractional_log_mean.py
python3 scripts/gasket_fractional_endpoint_shift.py
python3 scripts/gasket_fractional_supercritical_shift.py
python3 scripts/gasket_fractional_extension_limit.py
python3 scripts/gasket_fractional_corner_layer.py
python3 scripts/gasket_fractional_uniform_limit.py
python3 scripts/gasket_fractional_corner_sharpness.py
python3 scripts/gasket_fractional_level9.py
python3 scripts/gasket_fractional_level10.py
```

The first script exits nonzero unless the level-1 triple is \((3/2,-6/5,-3/10)\) and levels 2, 3, 4 match the factor \((3/5)^{n-1}\), including one arbitrary corner triple at levels 1 and 2. The second exits nonzero unless fractional currents at \(\alpha=0.45\) stay neutral through level 4, the level-2 audit norm matches \(1.165909\) within \(5\cdot 10^{-3}\), the refinement ratios differ from \(3/5\) by more than \(0.15\), and the \(\alpha=1\) control recovers Theorem G. The third exits nonzero unless geometric \(1/d^2\) currents match \(4^n\) times the combinatorial ones through level 5, the audit norm equals \((12/5)^n\sqrt{21}/2\), and every successive ratio equals \(12/5\). The fourth exits nonzero unless uniform multiplicative weights recover the closed form \((\lambda\cdot 3/5)^n\sqrt{21}/2\) through level 5 (including the resistance constant \(\sqrt{21}/2\)), and the resistance-scaled fractional audit at \(\alpha=0.45\) is not constant (Kept Failure S.1). The fifth exits nonzero unless multi-scale weights recover \((\sqrt{21}/2)\sum_{k=0}^{n}(\lambda\cdot 3/5)^k\) through level 5, nested harmonicity holds, and resistance multi-scale grows like \((n+1)\sqrt{21}/2\) (Kept Failure W.1). The sixth exits nonzero unless multi-scale spectral \(L^{0.45}\) stays neutral and on the Schur line through level 6, recovers Theorem W at \(\alpha=1\), shows nested harmonicity fails for \(\alpha=0.45\), rejects the \(\lambda=1\) freeze inheritance (Kept Failure X.1), and rejects coincidence with finest-edge fractional (Kept Failure X.2). The seventh exits nonzero unless, through level \(7\), the \(\alpha=0.45\) audit norm matches the level-2 lock \(1.165909\), stays neutral and inside the Theorem J cap, breaks both the \(0.95\) ratio law and the \(5^{-0.45}\) gap contraction, obeys the Theorem L energy sandwich (\(\|F_{\mathrm{frac}}\|\ge w_{\min}\|F\|\) and \(\|F\|\le\kappa\sqrt{Q}\)), breaks the \(0.95\) decay law for \(Q_{0.45}\), and the \(\alpha=0.9\) norm stays under the explicit Theorem K bound. `python3 scripts/gasket_fractional_limit.py 8` recomputes one level higher; the default witness does not. The eighth script exits nonzero unless the closed-form path resistance matches a dense pseudo-inverse, stays below \(C(\alpha)\) through \(m=2^{22}\) for \(\alpha<1/2\), and Löwner–Heinz, path duality, the energy floor and ceiling, and the pairing \(\|F_{\mathrm{frac}}\|=Q\sqrt{21}/7\) hold live on levels \(1\)–\(6\) (Theorem AA, Kept Failure AA.1). The ninth script exits nonzero unless altitude reflection equalizes the two neighbor weights, \(\sigma=\rho/w_n\) holds with \(w_{\min}\le w_n\le 8^{\alpha}\), the sandwich and Schur identities hold through level \(6\), and \(\rho_{5}/\rho_{4}>0.99\) (Theorem AB, Kept Failure AB.1). The tenth script exits nonzero unless the cell corner resistance is \((2/3)(5/3)^k\) with diameter under \((8/3)(5/3)^k\), \(\lambda_1\ge 5^{-n}/4\), the spectral mass of \(e_a-e_b\) stays under \(80s^{\gamma}\), and \(R_{\alpha}(n)\le C_{AC}(\alpha)\) with \(Q_{\alpha}(u_n)\ge(9/4)/R_{\alpha}(n)\) through level \(6\) for \(\alpha\in[0.45,0.68]\), and \(Q_{0.75}\), \(Q_{0.9}\) fall under the Theorem K bound (Theorem AC, Kept Failure AC.1). The eleventh script exits nonzero unless the corner Schur complement is \(\kappa(3I-J)\) with \(Q_{\alpha}R_{\alpha}=7/3\), the corner bump energy is \(2(3/5)^k\), \(R^{s}\ge(5s)^{\gamma-1}/5\) on \([5^{-n},1)\) through level \(7\), the Stieltjes quadrature matches \(R_{\gamma}(3)\), and \(R_{\gamma}(n)\), \(Q_{\gamma}(u_n)\), \(\|F(n,\gamma)\|\) obey the Theorem AD sandwich through level \(7\) while the Theorem K majorant stays constant (Theorem AD, Kept Failure AD.1). The twelfth script exits nonzero unless \(w_n\le8^{\alpha}/2\), the factorization \(\rho=\lambda/(\lambda+(1-\lambda)\gamma)\) and the AD resistance identities hold through level \(6\), and \(\lambda_4/\lambda_3>0.99\) (Theorem AE, Kept Failure AE.1). The thirteenth script exits nonzero unless \((L^{\alpha})_{aa}\le2^{\alpha}\le8^{\alpha}\), \(\lambda_n\ge\alpha\,2^{\alpha-1}\ge\lambda_{\min}(\alpha)\), \(w_n\le2^{\alpha-1}\), and \(\rho_n\ge1/(1+C_J\gamma_n)\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at six other \(\alpha\)), AE's factorization still holds, and every \(\gamma_{n+1}/\gamma_n<1\) (Theorem AF, Kept Failure AF.1). The fourteenth script exits nonzero unless \(W_{\mathrm{far}}\ge 2\alpha(1-\alpha)/4^{2-\alpha}\), \(\lambda_n\le 1-\alpha(1-\alpha)/2^{3-\alpha}\), the sandwich \(1/(1+C_J\gamma_n)\le\rho_n\le 1/(1+c_*\gamma_n)\), and the level-independent distance-\(2\) path counts hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at six other \(\alpha\)), and every \(\lambda_{n+1}/\lambda_n<1\) (Theorem AG, Kept Failure AG.1). The fifteenth script exits nonzero unless the subordination quadrature, the small-time bound \(\int_0^s\le K\int_s^{2s}\), \(H_s(p,a)\ge se^{-4s}\), and the kernel chain \(w_{pj}\ge\kappa w_{aj}\) hold at levels \(1\)–\(3\), and \(v=u_a-u\ge0\), \((L^{\alpha})_{pp}\le4^{\alpha}\), \(I_{\mathrm{far}}\le(4^{\alpha}/(2\kappa))I_a\), \(\gamma_n\le\Gamma_{AH}(\alpha)\), \(\sigma_n\ge\sigma_{\min}\), and \(\|F\|\ge\sigma_{\min}\|F_{\mathrm{frac}}\|\ge6.44\cdot10^{-4}\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the max-principle shortcut bound is negative at levels \(3\)–\(6\) (Theorem AH, Kept Failure AH.1). The sixteenth script exits nonzero unless the backward equation and Duhamel formula at the degree-\(2\) corner, \(0<K_{\alpha}(s)\le s^{-1-\alpha}/2\), and the weighted subordination identity for \(w_{aj}\) hold at levels \(1\)–\(3\), and \(w_{pj}+w_{qj}\ge2w_{aj}\) on Far, \(I_{\mathrm{far}}\le(4^{\alpha}/2)I_a\), \(\gamma_n\le\Gamma_{AI}(\alpha)\), \(\sigma_n\ge\sigma_{AI}\), \(\rho_n\ge\rho_{AI}\), and \(\|F\|\ge\sigma_{AI}\|F_{\mathrm{frac}}\|\ge0.051988\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the \(\kappa\to1\) chain stays below \(0.70\) (Theorem AI, Kept Failure AI.1). The seventeenth script exits nonzero unless \(H_t(p,q)\ge te^{-4t}\), \(w_{pq}\ge\alpha4^{\alpha-1}\), the exact neighbor identity with cross term, \(I_{\mathrm{far}}\le((4^{\alpha}-w_{pq})/2)I_a\le H_{AJ}I_a\), \(\gamma_n\le\Gamma_{AJ}\), \(\sigma_n\ge\sigma_{AJ}\), \(\rho_n\ge\rho_{AJ}\), and \(\|F\|\ge\sigma_{AJ}\|F_{\mathrm{frac}}\|\ge0.055599\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the \(w_{pq}\) chain stays below \(0.70\) (Theorem AJ, Kept Failure AJ.1). The eighteenth script exits nonzero unless the cell incidence identity \(L=3\,\mathrm{diag}(c)-N^{\mathsf T}N\), \(\lambda_{\max}\le6\) (attained at levels \(2\)–\(4\)), the local moments, and the Hermite majorant hold at levels \(1\)–\(4\), and \((L^{\alpha})_{aa}\le D_a^{+}\), \((L^{\alpha})_{pp}\le D_p^{+}\), \(w_n\le w_n^{+}\), \(I_{\mathrm{far}}\le H_{AK}I_a\), \(\gamma_n\le\Gamma_{AK}\), \(\sigma_n\ge\sigma_{AK}\), \(\rho_n\ge\rho_{AK}\), and \(\|F\|\ge\sigma_{AK}\|F_{\mathrm{frac}}\|\ge0.058941\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the live-constant neighbor-Harnack chain stays below \(1\) in \(\sigma\) and below \(0.70\) in \(\|F\|\) from level \(4\) (Theorem AK, Kept Failure AK.1). The nineteenth script exits nonzero unless the path integrand \(f(\theta)=(2\sin\theta)^{-2\alpha}\cos^2\theta\) has nonnegative second-derivative minimum for \(\alpha\ge\alpha_*=(14-8\sqrt{2})/34\), \(R_{\alpha}(m)\le R_{\infty}(\alpha)<C(\alpha)\) through \(m\le64\) and large dyadic \(m\) at four \(\alpha\in[\alpha_*,1/2)\), the closed form matches a dense pseudo-inverse, and live \(Q\), \(\|F_{\mathrm{frac}}\|\), \(\|F\|\) dominate \(0.169873\), \(0.111208\), \(0.078019\) through level \(6\) at \(0.45\), while \(\sigma_{AK}\cdot0.111208\), \(2^{1-\alpha}\cdot0.111208\), and live-\(\sigma\cdot0.111208\) all stay below \(0.70\) (Theorem AL, Kept Failure AL.1). The `gasket_fractional_moment_share.py` script exits nonzero unless the Radau law \(\{0,t_1,t_2\}\) (weight \(3/17\) at \(0\)) and the Lobatto law \(\{0,s_1,s_2,6\}\) reproduce the corner moments \(1,2,6,22,92\) (and \(418\)), the Hermite quartic and quintic stay below \(x^{1-\alpha}\), \(S_{\mathrm{Lob}}>S_{AP}>\sigma_{AK}\) on \(\alpha\in[0.01,0.99]\), and, through level \(6\) at \(0.45\) (level \(5\) at seven other \(\alpha\)), the moments, \(Lu=L^{1-\alpha}L^{\alpha}u\), \(I^{\alpha}\propto(1,-4/5,-1/5)\), \(\sigma_n=(L^{1-\alpha})_{aa}+\omega_n\) with \(\omega_n>0\), \((L^{1-\alpha})_{aa}\ge S_{AP}\) (\(\ge S_{\mathrm{Lob}}\) for \(n\ge2\)), and \(\|F\|\ge S_{AP}\|F_{\mathrm{frac}}\|\ge0.415238\) hold, while every uniform share floor times the AM/AO fractional floors stays below \(0.70\) and the needed resistance ceiling is \(2.859237\) (Theorem AP, Kept Failure AP.1). The `gasket_fractional_decimation.py` script exits nonzero unless all of the following hold: the one-cell Schur complement, the decimation map \(s\mapsto s(s+5)\), and the recursion algebra hold symbolically; the graph Schur complement of \(L_{n+1}+s\) onto \(V_n\) equals \(cL_n+(2e+s)I-eP_C\) for \(n=1\)–\(4\); the recursion equals the dense corner resolvent on levels \(0\)–\(6\); the two Stieltjes integrals reproduce the dense \(R_{\alpha}(n)\) and \(\sigma_n\) (levels \(1\)–\(6\) at \(0.45\), levels \(1\)–\(5\) at four other \(\alpha\)) and the recorded level-\(9\)/\(10\) values; \(R\) increases and \(\sigma\), \(\|F\|\) decrease; the interval certificate gives \(R_{0.45}(n)<2.518955\), \(\sigma_n>1.350986\), and \(\mathcal L\in[0.819255,0.822959]\) with live levels inside; the analogous certificates hold at \(\alpha\in\{0.10,0.25,0.49,0.55\}\); and \(R_{\alpha}(n)\) grows linearly at \(\log3/\log5\) and geometrically at \(0.75\) (Theorem AQ, Kept Failure AQ.1). The `gasket_fractional_limit_digits.py` script exits nonzero unless the tail moments \(2,6\), complete monotonicity, and the interpolation sign rule hold, AQ's power ceiling re-certifies, every certified enclosure of \(\mathcal L(\alpha)\) has width \(<10^{-10}\) (\(<10^{-7}\) at \(0.55\)) and sits inside AQ's bracket at \(0.45\), \(\mathcal L(0.10)>1.410\) (Conjecture J.2), and \(\rho\)-extrapolation from levels \(9,10\) misses \(\mathcal L(0.45)\) by more than \(10^{-4}\) (Theorem AS, Kept Failure AS.1). The `gasket_fractional_rate_constants.py` script exits nonzero unless the exact transfer identity holds at finite levels, AQ's power ceiling and the tail lemmas re-certify, every certified enclosure of \(c_R(\alpha)\) and \(c_S(\alpha)\) at \(\alpha\in\{0.25,0.45,0.60,\gamma,0.75\}\) has width \(<10^{-15}\) and contains AR's float \(c_R\) to \(10^{-8}\), the finite-level scaled increments approach the certified constants within \(5\cdot5^{-0.65m}\), \((\|F(30)\|-\mathcal L)/\rho^{30}\) matches \(c_F\), and the certified \(\Theta\) enclosures at \(\gamma\) are separated (Theorem AT, Kept Failure AT.1). The `gasket_fractional_log_mean.py` script exits nonzero unless AQ's power ceiling, the AT.1 orbit lemma and the Koenigs bounds re-certify, the node enclosures of \(r_\infty\) have relative width \(<10^{-60}\), the certified \(\langle\Theta\rangle\) has width \(<10^{-30}\) and agrees with a shifted period, \(\frac{2\sin\pi\gamma}{\pi}\ln5\,\langle\Theta\rangle\) lies inside Theorem AT's certified \(c_R(\gamma)\), \(\hat\theta_1,\dots,\hat\theta_4\) are certified nonzero, the Fourier series reproduces AT.1's eight \(\Theta\) values, and \(J_m-m\ln5\,\langle\Theta\rangle\) stays bounded and converges (Theorem AU, Kept Failure AU.1). The `gasket_fractional_endpoint_shift.py` script exits nonzero unless AQ's power ceiling and the AT small-\(y\) lemma re-certify, AU's \(\kappa\ln5\langle\Theta\rangle\) lies inside AT's \(c_R(\gamma)\), the node enclosures are ordered \(r_{40}\le r_{60}\le r_\infty\) with \(r_\infty\) to relative width \(<10^{-60}\), the levels-\(60\) and \(40\) enclosures of \(C'\) overlap, \(E_{40},E_{60}<5^{-0.6n}\), the certified \(C_\gamma\) and \(\sigma_\infty(\gamma)\) have width \(<10^{-25}\), quadrature at \(m=5\)–\(30\) satisfies \(0<J_m-mc-C'<5\cdot5^{-0.65m}\) and \(|(m+\nu_\gamma)\|F(m,\gamma)\|-K_\gamma|<5\cdot5^{-0.65m}\), and \(m(m\|F\|-K_\gamma)\) tracks \(-\nu_\gamma K_\gamma m/(m+\nu_\gamma)\) (Theorem AY, Kept Failure AY.1). The `gasket_fractional_supercritical_shift.py` script exits nonzero unless AQ's power ceiling and the AT small-\(y\) lemma re-certify, the node enclosures are ordered \(r_{40}\le r_{41}\le r_{60}\le r_{61}\le r_\infty\), the levels-\(60\) and \(40\) enclosures of \(c\) and \(C'\) at \(\alpha=0.8\) overlap, \(E_{40},E_{60}<\rho^n5^{-0.6n}\), the certified \(c_R,C,\sigma_\infty,K\) at \(0.8\) have width \(<10^{-20}\), \(C_{0.8},\nu_{0.8},C_{0.75}<0\), the level-\(60\) \(c_R(0.75)\) meets Theorem AT's enclosure, quadrature at \(m=5\)–\(30\) satisfies \(|R(m)-A\rho^m-C|<5\rho^m5^{-0.65m}\) and \(|(\rho^m+\nu)\|F\|-K|<5\cdot5^{-0.65m}\), \(\rho^m\|F\|-K>0\) decreases with \(\rho^m(\rho^m\|F\|-K)\) tracking \(-\nu K\), and the pole and continuity checks at \(\gamma\pm0.01\) hold (Theorem AZ, Kept Failure AZ.1). The level-\(9\) and level-\(10\) scripts replay recorded matrix-free values against a live level-\(2\) lock; the level-\(10\) script exits nonzero unless \(\delta_{10}/\delta_9>\log 3/\log 5\) (Kept Failure J.4), the norm still decreases above \(0.70\), and the Schur line, Theorem J cap, and Theorem L sandwich hold. To recompute level \(10\) itself (about half an hour; needs SciPy): `python3 scripts/gasket_fractional_matrix_free.py --level 10`.
