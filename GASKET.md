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

- **LIMIT AT \(\alpha=0.45\): SETTLED BY THEOREM AQ (the history below predates it).** Theorem AQ proves that \(\|F(n,0.45)\|\) strictly decreases to \(\mathcal L\in[0.819255,0.822959]\) by exact spectral decimation of the corner resolvent plus an interval certificate. Kept Failure AQ.1 shows that no positive limit exists for \(\alpha\ge\log3/\log5\). Before AQ: levels \(1\) through \(7\) still have increasing ratios, now past \(0.97\). Theorem J excludes \(+\infty\). Theorem K gives \(\|F\|\to 0\) only for \(\alpha>\log 3/\log 5\). Theorem L says energy collapse would force \(\|F\|\to 0\) at every \(\alpha\in(0,1)\), including \(0.45\), but does not prove energy collapse. Neither decides zero versus positive at \(\alpha=0.45\); Theorem AH does (not zero). Conjecture J.1 is not a theorem. Kept Failures J.1, J.2, J.3, and J.4 reject four sufficient conditions that would have closed it. Level \(10\) is recorded, not a limit. Theorem S settles the harmonic audit under uniform multiplicative weights (finite nonzero only at resistance \(\lambda=5/3\)); Kept Failure S.1 rejects inheriting that freeze for \(L^{0.45}\). Theorem W settles the multi-scale harmonic audit (finite for \(\lambda<5/3\); resistance diverges); Kept Failure W.1 rejects a multi-scale resistance freeze. Theorem X takes spectral \(L^{\alpha}\) on the multi-scale graph: neutrality and Schur hold, nested harmonicity is special to \(\alpha=1\); Kept Failures X.1–X.2 reject inheriting the \(\lambda=1\) freeze and reject coincidence with finest-edge fractional. Theorem AA gives the level-uniform floors \(Q_{0.45}\ge 0.128335\) and \(\|F_{\mathrm{frac}}\|\ge 0.084015\) (a floor exists at every \(\alpha<1/2\)); Kept Failure AA.1 closes the energy-collapse route to zero. Theorem AB factors \(\sigma_n=\|F\|/\|F_{\mathrm{frac}}\|=\rho_n/w_n\) and, with AA, equates \(\|F\|\to 0\) with \(\rho_n\to 0\); Kept Failure AB.1 rejects geometric collapse of \(\rho\). Theorem AC extends the energy floor to every \(\alpha<\log 3/\log 5\) by cell resistance and spectral mass, so AB's equivalence holds on that whole range; Kept Failure AC.1 (via Theorem K) rejects a floor above it. Theorem AD proves the exact pairing \(Q_{\alpha}R_{\alpha}=7/3\) and settles the endpoint: \(Q_{\gamma}(u_n)=\Theta(1/n)\) and \(\|F(n,\gamma)\|\to 0\); Kept Failure AD.1 rejects extending the AC floor to \(\gamma\). Theorem AE factors \(\rho_n=\lambda_n/(\lambda_n+(1-\lambda_n)\gamma_n)\) with sharp ceiling \(w_n\le8^{\alpha}/2\) and the resistance identity \(\rho_n=(\sqrt{21}/7)w_n\|F\|R_{\alpha}\); Kept Failure AE.1 rejects geometric wipeout of the neighbor-degree share \(\lambda_n\). Theorem AF floors that share, \(\lambda_n\ge\alpha\,2^{-1-\alpha}\) (Jensen: \(\ge\alpha\,2^{\alpha-1}\)), so for \(\alpha<\log3/\log5\) vanishing of \(\|F\|\) requires \(\gamma_n\to\infty\) (one-sided); Kept Failure AF.1 rejects geometric blowup of \(\gamma_n\). Theorem AG ceilings the same share, \(\lambda_n\le 1-\alpha(1-\alpha)/2^{3-\alpha}<1\), so \(\|F\|\to0\) if and only if \(\gamma_n\to\infty\) (two-sided); Kept Failure AG.1 rejects geometric climb of \(\lambda\) toward \(1\). Theorem AH compares the jump kernels of \(p\) and \(a\) (\(w_{pj}\ge\kappa(\alpha)w_{aj}\) off the corner patch) and applies one-step \(L^{\alpha}\)-harmonicity at the two neighbors, so \(\gamma_n\le\Gamma_{AH}(\alpha)<\infty\) and \(\|F(n,\alpha)\|\ge\sigma_{\min}(\alpha)\|F_{\mathrm{frac}}\|>0\) at every level for every \(\alpha<\log3/\log5\) (\(\ge 6.44\cdot10^{-4}\) at \(0.45\)); Kept Failure AH.1 rejects the max-principle shortcut to the same floor. Theorem AI replaces \(\kappa\) by \(1\) through the first-step identity at the degree-\(2\) corner (\(w_{pj}+w_{qj}\ge2w_{aj}\) on Far), so \(I_{\mathrm{far}}\le(4^{\alpha}/2)I_a\), \(\gamma_n\le8/(\alpha(1-\alpha))\), and \(\|F(n,0.45)\|\ge0.051988\) at every level; Kept Failure AI.1 rejects reading that as Conjecture J.1's \(0.70\). Theorem AJ keeps the neighbor-edge weight \(w_{pq}\ge\alpha4^{\alpha-1}\) in the Harnack identity, so \(I_{\mathrm{far}}\le H_{AJ}I_a\) with \(H_{AJ}=(4^{\alpha}-\alpha4^{\alpha-1})/2\), \(\gamma_n\le2(4-\alpha)/(\alpha(1-\alpha))\), and \(\|F(n,0.45)\|\ge0.055599\); Kept Failure AJ.1 rejects reading that as \(0.70\). Theorem AK uses the second moment and the sharp line-graph ceiling \(\mathrm{spec}(L)\subset[0,6]\) to lower the diagonal ceilings to \((L^{\alpha})_{aa}\le6^{\alpha}/9+(8/9)(3/2)^{\alpha}\) and \((L^{\alpha})_{pp}\le(6^{\alpha}+2^{\alpha})/2\), so \(\|F(n,0.45)\|\ge0.058941\); Kept Failure AK.1 shows that the one-step neighbor-Harnack chain stays below \(0.70\) from level \(4\) even with every constant live. Theorem AL replaces AA's path ceiling \(C(\alpha)\) by the Hermite–Hadamard limit \(R_{\infty}(\alpha)\) for \(\alpha\ge\alpha_*=(14-8\sqrt{2})/34\), so \(\|F_{\mathrm{frac}}(n,0.45)\|\ge0.111208\) and \(\|F(n,0.45)\|\ge0.078019\); Kept Failure AL.1 shows that this sharp path ceiling still cannot reach \(0.70\). Theorem AM replaces the path by the full boundary cycle \(C_{3\cdot2^n}\), whose fractional corner resistance is below \((1-\alpha)R_{\infty}(\alpha)\) for every \(\alpha<1/2\) (Fourier plus a heat-kernel periodization), and uses all three corners, so \(Q_{0.45}\ge0.320299\), \(\|F_{\mathrm{frac}}(n,0.45)\|\ge0.209685\), and \(\|F(n,0.45)\|\ge0.147106\); Kept Failure AM.1 shows this still cannot reach \(0.70\), and that even the live energy times \(\sigma_{AK}\) stays below \(0.70\) from level \(2\). Theorem AN replaces the cycle by the boundary fence/collar with Bloch ceiling \(R^{\mathrm{fen}}_{\infty}+1/2\approx5.002104\), so \(Q\ge0.466470\), \(\|F_{\mathrm{frac}}\|\ge0.305377\), and \(\|F(n,0.45)\|\ge0.214239\); Theorem AO replaces that by the depth-2 collar with period-2 Bloch ceiling \(R^{\mathrm{d2}}_{\infty}+11/20\approx4.820047\), so \(Q\ge0.484089\), \(\|F_{\mathrm{frac}}\|\ge0.316911\), and \(\|F(n,0.45)\|\ge0.222331\); Kept Failures AN.1 and AO.1 show those routes stay below \(0.70\). Theorem AP proves the exact share identity \(\sigma_n=(L^{1-\alpha})_{aa}+\omega_n\) (\(Lu=L^{1-\alpha}L^{\alpha}u\) and AD's corner currents) and a five-moment Gauss–Radau floor \(\sigma_n\ge S_{AP}(0.45)\approx1.310267>1\), beating AK's \(0.701556\) and AK.1's \(0.9529\) cap, so \(\|F(n,0.45)\|\ge0.274743\) with AM's floor and \(\ge0.415238\) with AO's floor as stated; Kept Failure AP.1 shows no share floor can lift the proved fractional floors to \(0.70\), so the open step is now a resistance ceiling \(R_{\alpha}(n)\le2.859237\). Neither floor is a limit theorem. Theorem AQ is one: \(\|F(n,0.45)\|\downarrow\mathcal L\in[0.819255,0.822959]\) (Conjecture J.1 proved). Theorem AR adds the exact rates: \(R_{\alpha}(n+1)-R_{\alpha}(n)\sim c_R(\alpha)(5^{\alpha}/3)^n\) for every \(\alpha\in(0,1)\), so the gap ratio at \(0.45\) tends to \(5^{0.45}/3\), \(R_{\gamma}(n)\) has an exact linear slope, and \(\|F(n,\alpha)\|=\Theta((3\cdot5^{-\alpha})^n)\) above \(\gamma\); Kept Failure AR.1 shows the universal resolvent bound alone gives no rate. Theorem AS certifies the digits: \(r_\infty\) and \(1-sr_\infty\) are completely monotone, so polynomial interpolants on a geometric grid are rigorous floors and ceilings, and outward-rounded sums give \(\mathcal L(0.45)\in[0.82113862713346,0.82113862713352]\) and \(\mathcal L(0.10)>1.410638933063\), which proves Conjecture J.2; Kept Failure AS.1 shows that level-\(10\) extrapolation with AR's exact ratio still misses \(\mathcal L\) by \(2.8\cdot10^{-4}\). Theorem AT certifies AR's rate constants: an exact change of variables moves \(\rho^{-m}\Delta R_m\) to level \(n\) with weights trapped in a bracket of relative width \(O(5^{-0.65n})\), so at \(n=40\) \(c_R(0.45)\in[0.497802650763560354,0.497802650763560356]\), \(c_S(0.45)\in[0.927092638569461377,0.927092638569461384]\), \(c_R(\gamma)\in[0.517165834630505243,0.517165834630505245]\), and \(\|F(n,0.45)\|-\mathcal L\sim c_F\rho^n\) with \(c_F\in[0.52027873061815,0.52027873061821]\); Kept Failure AT.1 shows \(s^{1-\gamma}r_\infty(s)\) has no limit (certified log-periodic oscillation), so the slope at \(\gamma\) is not a pure-power constant. Theorem AU proves Conjecture AT.2: the slope is the log-period mean, \(c_R(\gamma)=\frac{2\sin\pi\gamma}{\pi}\ln5\,\langle\Theta\rangle\), certified to \(34\) digits inside AT's enclosure; the Mellin transform of \(r_\infty\) has simple poles at \(1-\gamma+2\pi ik/\ln5\), certified genuine for \(|k|\le4\); Kept Failure AU.1 rejects a single-harmonic profile. Theorem AY certifies the intercept: at \(\gamma\) the bracket weights of AT do not depend on the level, so the per-level defects \(|\Delta_j-c|\) sum to at most \(\int t^{-\gamma}(r_\infty-r_n)W_\pm\), and \(C_\gamma\in[0.70307046242185742992904816,\ldots817]\), \(\nu_\gamma=C_\gamma/c_R(\gamma)=1.3594681151437\ldots\), \(\sigma_\infty(\gamma)=1.1519186116740\ldots\), so \((n+\nu_\gamma)\|F(n,\gamma)\|=K_\gamma+O(5^{-0.65n})\) with \(K_\gamma=3.4023607638326\ldots\); Kept Failure AY.1 rejects the unshifted renormalization \(n\|F(n,\gamma)\|\), whose error is \(\Theta(1/n)\).
- **NOT A RESIDUAL.** Their tilt diagnostic sets \((Lu)_i=1\) on the interior. Theorem F then says the corner currents sum to \(-(N-3)\), the source that was inserted. That sum grows like \(3^n\). It is not a force left after the source is removed.
- **FAILED as a continuum thrust for the combinatorial harmonic case.** Under Theorem G, \(\|F(n)\|\to 0\) for every fixed corner triple. The combinatorial harmonic limit does not keep a finite audit vector.
- **SETTLED (Theorem AQ).** \(\|F(n)\|\) for \(L^{0.45}\) strictly decreases and converges to \(\mathcal L\in[0.819255,0.822959]\). Its closed-form value is not known. History: it does not tend to \(0\) (Theorem AH: \(\|F(n,0.45)\|\ge 6.44\cdot10^{-4}\) at every level; Theorem AI sharpens this to \(\ge0.051988\); Theorem AJ to \(\ge0.055599\); Theorems AK, AL, AM, AN, AO to \(\ge0.058941\), \(\ge0.078019\), \(\ge0.147106\), \(\ge0.214239\), \(\ge0.222331\); Theorem AP to \(\ge0.274743\) with AM's floor and \(\ge0.415238\) with AO's floor as stated). Divergence is excluded (Theorem J). Neutrality is settled (Theorem H). For \(\alpha>\log 3/\log 5\) the same audit norm tends to \(0\) (Theorem K). Energy controls the audit norm at every \(\alpha\in(0,1)\) (Theorem L), so energy collapse would force \(\|F\|\to 0\), and a positive audit-norm liminf would force a positive energy liminf; the energy does not vanish at \(\alpha=0.45\) (Theorem AA: \(Q_{0.45}(u_n)\ge 0.128335\) and \(\|F_{\mathrm{frac}}\|\ge 0.084015\) at every level, for the same reason at every \(\alpha<1/2\)), so energy collapse cannot send \(\|F\|\) to \(0\) (Kept Failure AA.1); Theorem AB equates vanishing of \(\|F\|\) with vanishing of the two-neighbor fractional share \(\rho_n\) (Kept Failure AB.1 rejects geometric collapse of \(\rho\)); Theorem AC proves the energy floor for every \(\alpha<\log 3/\log 5\), the exact threshold by Theorem K (Kept Failure AC.1), and Theorem AD settles \(\alpha=\log 3/\log 5\) itself: \(Q_{\gamma}(u_n)=\Theta(1/n)\to 0\) and \(\|F(n,\gamma)\|\to 0\) (Kept Failure AD.1); Theorem AE splits \(\rho_n\to0\) into degree-share and gap-ratio blowup and rejects geometric \(\lambda\) wipeout (Kept Failure AE.1); Theorem AF floors \(\lambda_n\) uniformly, so \(\|F\|\to0\) would force the far-to-neighbor gap ratio \(\gamma_n\to\infty\), and Kept Failure AF.1 rejects geometric \(\gamma\) blowup (live \(\gamma_n\) decreases through level \(6\)); Theorem AG ceilings \(\lambda_n\) uniformly below \(1\), so \(\|F\|\to0\) if and only if \(\gamma_n\to\infty\), and Kept Failure AG.1 rejects geometric climb of \(\lambda\) toward \(1\); Theorem AH ceilings \(\gamma_n\le\Gamma_{AH}(\alpha)\) by a kernel comparison and one-step Harnack at the neighbors, so \(\|F\|\) does not vanish for any \(\alpha<\log3/\log5\), and Kept Failure AH.1 rejects the max-principle shortcut; Theorem AI sharpens the Harnack constant from \(4^{\alpha}/(2\kappa)\) to \(4^{\alpha}/2\) by a two-neighbor heat-kernel sum, giving \(\gamma_n\le8/(\alpha(1-\alpha))\) and \(\|F(n,0.45)\|\ge0.051988\), and Kept Failure AI.1 rejects reaching \(0.70\) that way; Theorem AJ keeps \(w_{pq}\) to push the floor to \(0.055599\), and Kept Failure AJ.1 rejects reaching \(0.70\) that way; Theorem AK sharpens the diagonal ceilings by two-moment Jensen on \([0,6]\) to push it to \(0.058941\), and Kept Failure AK.1 shows the whole one-step neighbor-Harnack chain stays below \(0.70\) from level \(4\) even with live constants; Theorem AL sharpens the path-resistance ceiling to \(R_{\infty}\) for \(\alpha\ge\alpha_*\) to push the floor to \(0.078019\), and Kept Failure AL.1 shows that route also stays below \(0.70\); Theorem AM replaces the path by the boundary cycle to push the floor to \(0.147106\), Theorem AN replaces the cycle by the boundary fence/collar (Bloch ceiling \(R^{\mathrm{fen}}_{\infty}+1/2\)) to push it to \(0.214239\), and Theorem AO replaces that by the depth-2 collar (period-2 Bloch ceiling \(R^{\mathrm{d2}}_{\infty}+11/20\)) to push it to \(0.222331\); Kept Failures AM.1, AN.1, and AO.1 show those routes stay below \(0.70\), and that no resistance ceiling can rescue the \(\sigma_{AK}\) chain; Theorem AP replaces that chain by the exact identity \(\sigma_n=(L^{1-\alpha})_{aa}+\omega_n\) and the moment floor \(\sigma_n\ge1.310267\), and Kept Failure AP.1 shows no share floor can lift the proved fractional floors to \(0.70\), leaving a resistance ceiling \(R_{\alpha}(n)\le2.859237\) as the open step. Theorem AQ supplies it, \(R_{\alpha}(n)<2.518955\), by exact spectral decimation of the corner resolvent, and proves monotone convergence of \(\|F(n,0.45)\|\) to \(\mathcal L\in[0.819255,0.822959]\). Kept Failure AQ.1 shows that the uniform ceiling fails for \(\alpha\ge\log3/\log5\). Geometric weights \(1/d^2\) on this build are settled (Theorem I): that audit norm diverges. Uniform multiplicative weights are settled for the harmonic case (Theorem S): finite nonzero only at the resistance value \(\lambda=5/3\); resistance weights do not freeze the fractional \(\alpha=0.45\) audit (Kept Failure S.1). Multi-scale multiplicative weights are settled for the harmonic case (Theorem W): finite nonzero for every \(\lambda<5/3\), and resistance diverges like \((n+1)\sqrt{21}/2\) (Kept Failure W.1). Spectral \(L^{\alpha}\) on that multi-scale graph is settled structurally by Theorem X (neutrality, Schur, \(\alpha=1\) recovers W, no nested fractional harmonicity); the \(\lambda=1\) freeze does not transfer to \(\alpha=0.45\) (Kept Failure X.1), and multi-scale fractional does not coincide with finest-edge fractional (Kept Failure X.2). None of this is thrust.
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
python3 scripts/gasket_fractional_level9.py
python3 scripts/gasket_fractional_level10.py
```

The first script exits nonzero unless the level-1 triple is \((3/2,-6/5,-3/10)\) and levels 2, 3, 4 match the factor \((3/5)^{n-1}\), including one arbitrary corner triple at levels 1 and 2. The second exits nonzero unless fractional currents at \(\alpha=0.45\) stay neutral through level 4, the level-2 audit norm matches \(1.165909\) within \(5\cdot 10^{-3}\), the refinement ratios differ from \(3/5\) by more than \(0.15\), and the \(\alpha=1\) control recovers Theorem G. The third exits nonzero unless geometric \(1/d^2\) currents match \(4^n\) times the combinatorial ones through level 5, the audit norm equals \((12/5)^n\sqrt{21}/2\), and every successive ratio equals \(12/5\). The fourth exits nonzero unless uniform multiplicative weights recover the closed form \((\lambda\cdot 3/5)^n\sqrt{21}/2\) through level 5 (including the resistance constant \(\sqrt{21}/2\)), and the resistance-scaled fractional audit at \(\alpha=0.45\) is not constant (Kept Failure S.1). The fifth exits nonzero unless multi-scale weights recover \((\sqrt{21}/2)\sum_{k=0}^{n}(\lambda\cdot 3/5)^k\) through level 5, nested harmonicity holds, and resistance multi-scale grows like \((n+1)\sqrt{21}/2\) (Kept Failure W.1). The sixth exits nonzero unless multi-scale spectral \(L^{0.45}\) stays neutral and on the Schur line through level 6, recovers Theorem W at \(\alpha=1\), shows nested harmonicity fails for \(\alpha=0.45\), rejects the \(\lambda=1\) freeze inheritance (Kept Failure X.1), and rejects coincidence with finest-edge fractional (Kept Failure X.2). The seventh exits nonzero unless, through level \(7\), the \(\alpha=0.45\) audit norm matches the level-2 lock \(1.165909\), stays neutral and inside the Theorem J cap, breaks both the \(0.95\) ratio law and the \(5^{-0.45}\) gap contraction, obeys the Theorem L energy sandwich (\(\|F_{\mathrm{frac}}\|\ge w_{\min}\|F\|\) and \(\|F\|\le\kappa\sqrt{Q}\)), breaks the \(0.95\) decay law for \(Q_{0.45}\), and the \(\alpha=0.9\) norm stays under the explicit Theorem K bound. `python3 scripts/gasket_fractional_limit.py 8` recomputes one level higher; the default witness does not. The eighth script exits nonzero unless the closed-form path resistance matches a dense pseudo-inverse, stays below \(C(\alpha)\) through \(m=2^{22}\) for \(\alpha<1/2\), and Löwner–Heinz, path duality, the energy floor and ceiling, and the pairing \(\|F_{\mathrm{frac}}\|=Q\sqrt{21}/7\) hold live on levels \(1\)–\(6\) (Theorem AA, Kept Failure AA.1). The ninth script exits nonzero unless altitude reflection equalizes the two neighbor weights, \(\sigma=\rho/w_n\) holds with \(w_{\min}\le w_n\le 8^{\alpha}\), the sandwich and Schur identities hold through level \(6\), and \(\rho_{5}/\rho_{4}>0.99\) (Theorem AB, Kept Failure AB.1). The tenth script exits nonzero unless the cell corner resistance is \((2/3)(5/3)^k\) with diameter under \((8/3)(5/3)^k\), \(\lambda_1\ge 5^{-n}/4\), the spectral mass of \(e_a-e_b\) stays under \(80s^{\gamma}\), and \(R_{\alpha}(n)\le C_{AC}(\alpha)\) with \(Q_{\alpha}(u_n)\ge(9/4)/R_{\alpha}(n)\) through level \(6\) for \(\alpha\in[0.45,0.68]\), and \(Q_{0.75}\), \(Q_{0.9}\) fall under the Theorem K bound (Theorem AC, Kept Failure AC.1). The eleventh script exits nonzero unless the corner Schur complement is \(\kappa(3I-J)\) with \(Q_{\alpha}R_{\alpha}=7/3\), the corner bump energy is \(2(3/5)^k\), \(R^{s}\ge(5s)^{\gamma-1}/5\) on \([5^{-n},1)\) through level \(7\), the Stieltjes quadrature matches \(R_{\gamma}(3)\), and \(R_{\gamma}(n)\), \(Q_{\gamma}(u_n)\), \(\|F(n,\gamma)\|\) obey the Theorem AD sandwich through level \(7\) while the Theorem K majorant stays constant (Theorem AD, Kept Failure AD.1). The twelfth script exits nonzero unless \(w_n\le8^{\alpha}/2\), the factorization \(\rho=\lambda/(\lambda+(1-\lambda)\gamma)\) and the AD resistance identities hold through level \(6\), and \(\lambda_4/\lambda_3>0.99\) (Theorem AE, Kept Failure AE.1). The thirteenth script exits nonzero unless \((L^{\alpha})_{aa}\le2^{\alpha}\le8^{\alpha}\), \(\lambda_n\ge\alpha\,2^{\alpha-1}\ge\lambda_{\min}(\alpha)\), \(w_n\le2^{\alpha-1}\), and \(\rho_n\ge1/(1+C_J\gamma_n)\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at six other \(\alpha\)), AE's factorization still holds, and every \(\gamma_{n+1}/\gamma_n<1\) (Theorem AF, Kept Failure AF.1). The fourteenth script exits nonzero unless \(W_{\mathrm{far}}\ge 2\alpha(1-\alpha)/4^{2-\alpha}\), \(\lambda_n\le 1-\alpha(1-\alpha)/2^{3-\alpha}\), the sandwich \(1/(1+C_J\gamma_n)\le\rho_n\le 1/(1+c_*\gamma_n)\), and the level-independent distance-\(2\) path counts hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at six other \(\alpha\)), and every \(\lambda_{n+1}/\lambda_n<1\) (Theorem AG, Kept Failure AG.1). The fifteenth script exits nonzero unless the subordination quadrature, the small-time bound \(\int_0^s\le K\int_s^{2s}\), \(H_s(p,a)\ge se^{-4s}\), and the kernel chain \(w_{pj}\ge\kappa w_{aj}\) hold at levels \(1\)–\(3\), and \(v=u_a-u\ge0\), \((L^{\alpha})_{pp}\le4^{\alpha}\), \(I_{\mathrm{far}}\le(4^{\alpha}/(2\kappa))I_a\), \(\gamma_n\le\Gamma_{AH}(\alpha)\), \(\sigma_n\ge\sigma_{\min}\), and \(\|F\|\ge\sigma_{\min}\|F_{\mathrm{frac}}\|\ge6.44\cdot10^{-4}\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the max-principle shortcut bound is negative at levels \(3\)–\(6\) (Theorem AH, Kept Failure AH.1). The sixteenth script exits nonzero unless the backward equation and Duhamel formula at the degree-\(2\) corner, \(0<K_{\alpha}(s)\le s^{-1-\alpha}/2\), and the weighted subordination identity for \(w_{aj}\) hold at levels \(1\)–\(3\), and \(w_{pj}+w_{qj}\ge2w_{aj}\) on Far, \(I_{\mathrm{far}}\le(4^{\alpha}/2)I_a\), \(\gamma_n\le\Gamma_{AI}(\alpha)\), \(\sigma_n\ge\sigma_{AI}\), \(\rho_n\ge\rho_{AI}\), and \(\|F\|\ge\sigma_{AI}\|F_{\mathrm{frac}}\|\ge0.051988\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the \(\kappa\to1\) chain stays below \(0.70\) (Theorem AI, Kept Failure AI.1). The seventeenth script exits nonzero unless \(H_t(p,q)\ge te^{-4t}\), \(w_{pq}\ge\alpha4^{\alpha-1}\), the exact neighbor identity with cross term, \(I_{\mathrm{far}}\le((4^{\alpha}-w_{pq})/2)I_a\le H_{AJ}I_a\), \(\gamma_n\le\Gamma_{AJ}\), \(\sigma_n\ge\sigma_{AJ}\), \(\rho_n\ge\rho_{AJ}\), and \(\|F\|\ge\sigma_{AJ}\|F_{\mathrm{frac}}\|\ge0.055599\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the \(w_{pq}\) chain stays below \(0.70\) (Theorem AJ, Kept Failure AJ.1). The eighteenth script exits nonzero unless the cell incidence identity \(L=3\,\mathrm{diag}(c)-N^{\mathsf T}N\), \(\lambda_{\max}\le6\) (attained at levels \(2\)–\(4\)), the local moments, and the Hermite majorant hold at levels \(1\)–\(4\), and \((L^{\alpha})_{aa}\le D_a^{+}\), \((L^{\alpha})_{pp}\le D_p^{+}\), \(w_n\le w_n^{+}\), \(I_{\mathrm{far}}\le H_{AK}I_a\), \(\gamma_n\le\Gamma_{AK}\), \(\sigma_n\ge\sigma_{AK}\), \(\rho_n\ge\rho_{AK}\), and \(\|F\|\ge\sigma_{AK}\|F_{\mathrm{frac}}\|\ge0.058941\) hold through level \(6\) at \(\alpha=0.45\) (and through level \(5\) at seven other \(\alpha\)), and the live-constant neighbor-Harnack chain stays below \(1\) in \(\sigma\) and below \(0.70\) in \(\|F\|\) from level \(4\) (Theorem AK, Kept Failure AK.1). The nineteenth script exits nonzero unless the path integrand \(f(\theta)=(2\sin\theta)^{-2\alpha}\cos^2\theta\) has nonnegative second-derivative minimum for \(\alpha\ge\alpha_*=(14-8\sqrt{2})/34\), \(R_{\alpha}(m)\le R_{\infty}(\alpha)<C(\alpha)\) through \(m\le64\) and large dyadic \(m\) at four \(\alpha\in[\alpha_*,1/2)\), the closed form matches a dense pseudo-inverse, and live \(Q\), \(\|F_{\mathrm{frac}}\|\), \(\|F\|\) dominate \(0.169873\), \(0.111208\), \(0.078019\) through level \(6\) at \(0.45\), while \(\sigma_{AK}\cdot0.111208\), \(2^{1-\alpha}\cdot0.111208\), and live-\(\sigma\cdot0.111208\) all stay below \(0.70\) (Theorem AL, Kept Failure AL.1). The `gasket_fractional_moment_share.py` script exits nonzero unless the Radau law \(\{0,t_1,t_2\}\) (weight \(3/17\) at \(0\)) and the Lobatto law \(\{0,s_1,s_2,6\}\) reproduce the corner moments \(1,2,6,22,92\) (and \(418\)), the Hermite quartic and quintic stay below \(x^{1-\alpha}\), \(S_{\mathrm{Lob}}>S_{AP}>\sigma_{AK}\) on \(\alpha\in[0.01,0.99]\), and, through level \(6\) at \(0.45\) (level \(5\) at seven other \(\alpha\)), the moments, \(Lu=L^{1-\alpha}L^{\alpha}u\), \(I^{\alpha}\propto(1,-4/5,-1/5)\), \(\sigma_n=(L^{1-\alpha})_{aa}+\omega_n\) with \(\omega_n>0\), \((L^{1-\alpha})_{aa}\ge S_{AP}\) (\(\ge S_{\mathrm{Lob}}\) for \(n\ge2\)), and \(\|F\|\ge S_{AP}\|F_{\mathrm{frac}}\|\ge0.415238\) hold, while every uniform share floor times the AM/AO fractional floors stays below \(0.70\) and the needed resistance ceiling is \(2.859237\) (Theorem AP, Kept Failure AP.1). The `gasket_fractional_decimation.py` script exits nonzero unless all of the following hold: the one-cell Schur complement, the decimation map \(s\mapsto s(s+5)\), and the recursion algebra hold symbolically; the graph Schur complement of \(L_{n+1}+s\) onto \(V_n\) equals \(cL_n+(2e+s)I-eP_C\) for \(n=1\)–\(4\); the recursion equals the dense corner resolvent on levels \(0\)–\(6\); the two Stieltjes integrals reproduce the dense \(R_{\alpha}(n)\) and \(\sigma_n\) (levels \(1\)–\(6\) at \(0.45\), levels \(1\)–\(5\) at four other \(\alpha\)) and the recorded level-\(9\)/\(10\) values; \(R\) increases and \(\sigma\), \(\|F\|\) decrease; the interval certificate gives \(R_{0.45}(n)<2.518955\), \(\sigma_n>1.350986\), and \(\mathcal L\in[0.819255,0.822959]\) with live levels inside; the analogous certificates hold at \(\alpha\in\{0.10,0.25,0.49,0.55\}\); and \(R_{\alpha}(n)\) grows linearly at \(\log3/\log5\) and geometrically at \(0.75\) (Theorem AQ, Kept Failure AQ.1). The `gasket_fractional_limit_digits.py` script exits nonzero unless the tail moments \(2,6\), complete monotonicity, and the interpolation sign rule hold, AQ's power ceiling re-certifies, every certified enclosure of \(\mathcal L(\alpha)\) has width \(<10^{-10}\) (\(<10^{-7}\) at \(0.55\)) and sits inside AQ's bracket at \(0.45\), \(\mathcal L(0.10)>1.410\) (Conjecture J.2), and \(\rho\)-extrapolation from levels \(9,10\) misses \(\mathcal L(0.45)\) by more than \(10^{-4}\) (Theorem AS, Kept Failure AS.1). The `gasket_fractional_rate_constants.py` script exits nonzero unless the exact transfer identity holds at finite levels, AQ's power ceiling and the tail lemmas re-certify, every certified enclosure of \(c_R(\alpha)\) and \(c_S(\alpha)\) at \(\alpha\in\{0.25,0.45,0.60,\gamma,0.75\}\) has width \(<10^{-15}\) and contains AR's float \(c_R\) to \(10^{-8}\), the finite-level scaled increments approach the certified constants within \(5\cdot5^{-0.65m}\), \((\|F(30)\|-\mathcal L)/\rho^{30}\) matches \(c_F\), and the certified \(\Theta\) enclosures at \(\gamma\) are separated (Theorem AT, Kept Failure AT.1). The `gasket_fractional_log_mean.py` script exits nonzero unless AQ's power ceiling, the AT.1 orbit lemma and the Koenigs bounds re-certify, the node enclosures of \(r_\infty\) have relative width \(<10^{-60}\), the certified \(\langle\Theta\rangle\) has width \(<10^{-30}\) and agrees with a shifted period, \(\frac{2\sin\pi\gamma}{\pi}\ln5\,\langle\Theta\rangle\) lies inside Theorem AT's certified \(c_R(\gamma)\), \(\hat\theta_1,\dots,\hat\theta_4\) are certified nonzero, the Fourier series reproduces AT.1's eight \(\Theta\) values, and \(J_m-m\ln5\,\langle\Theta\rangle\) stays bounded and converges (Theorem AU, Kept Failure AU.1). The `gasket_fractional_endpoint_shift.py` script exits nonzero unless AQ's power ceiling and the AT small-\(y\) lemma re-certify, AU's \(\kappa\ln5\langle\Theta\rangle\) lies inside AT's \(c_R(\gamma)\), the node enclosures are ordered \(r_{40}\le r_{60}\le r_\infty\) with \(r_\infty\) to relative width \(<10^{-60}\), the levels-\(60\) and \(40\) enclosures of \(C'\) overlap, \(E_{40},E_{60}<5^{-0.6n}\), the certified \(C_\gamma\) and \(\sigma_\infty(\gamma)\) have width \(<10^{-25}\), quadrature at \(m=5\)–\(30\) satisfies \(0<J_m-mc-C'<5\cdot5^{-0.65m}\) and \(|(m+\nu_\gamma)\|F(m,\gamma)\|-K_\gamma|<5\cdot5^{-0.65m}\), and \(m(m\|F\|-K_\gamma)\) tracks \(-\nu_\gamma K_\gamma m/(m+\nu_\gamma)\) (Theorem AY, Kept Failure AY.1). The level-\(9\) and level-\(10\) scripts replay recorded matrix-free values against a live level-\(2\) lock; the level-\(10\) script exits nonzero unless \(\delta_{10}/\delta_9>\log 3/\log 5\) (Kept Failure J.4), the norm still decreases above \(0.70\), and the Schur line, Theorem J cap, and Theorem L sandwich hold. To recompute level \(10\) itself (about half an hour; needs SciPy): `python3 scripts/gasket_fractional_matrix_free.py --level 10`.
