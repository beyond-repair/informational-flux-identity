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
4. **ASSUMPTION, spectral calculus, used in Theorems H, J–L, X.** For the combinatorial Laplacian \(L\succeq 0\) with \(L\mathbf{1}=0\), set \(L^{\alpha}=V\mathrm{diag}(\lambda^{\alpha})V^{T}\) in an eigenbasis of \(L\), with the convention \(0^{\alpha}:=0\) for \(\alpha>0\). This is the same functional calculus as in the related repo's `fractional_laplacian`.
5. **ASSUMPTION, geometric weights, used only in Theorem I.** On the level-\(n\) combinatorial edge set of `build_gasket(n)`, set \(w_{ij}=1/\|P_i-P_j\|^2\). Every edge of that build has Euclidean length \(2^{-n}\), so \(w_{ij}=4^n\) uniformly and \(L_{\mathrm{geom}}=4^n L_{\mathrm{comb}}\).
6. **ASSUMPTION, uniform multiplicative weights, used only in Theorem S.** On the same level-\(n\) combinatorial edge set, set \(w_{ij}=\lambda^n\) for a fixed \(\lambda>0\). The combinatorial case is \(\lambda=1\); geometric \(1/d^2\) is \(\lambda=4\); Kigami resistance renormalization on this finest-edge build is \(\lambda=5/3\).
7. **ASSUMPTION, multi-scale edges, used in Theorems W, X.** Vertices are the level-\(n\) gasket points of `build_gasket(n)`. Edges are the three sides of every triangular cell at every generation \(g=0,\dots,n\) (so coarse edges are retained alongside finer ones). Generation-\(g\) edges get weight \(w=\lambda^g\). The finest-only build of Theorems G–S is the generation-\(n\) slice alone.
8. **NOT ASSUMED.** A unit source on the interior, any value of \(W\), a nonzero continuum limit of the fractional audit vector, or a continuum kernel \(1/|x-y|^{d+\alpha}\).

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

**NOT THIS.** The comparison does not force \(Q_{0.45}(u_{n})\to 0\). Theorem K's harmonic majorant still grows at \(\alpha=0.45\), and the actual minimizer energy is a different sequence. Zero versus a positive finite audit-norm limit remains OPEN; Theorem L only says that a positive audit-norm limit would need a positive energy liminf, and that energy collapse would settle the limit at \(0\).

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

## Conjecture J.1

**CONJECTURE, not a theorem.** For \(\alpha=0.45\) and these corner data, \(\|F(n)\|\) decreases for every \(n\ge 1\) and
\[
\lim_{n\to\infty}\|F(n)\| = L \quad\text{with}\quad L\ge 0.70.
\]
The only evidence is the computed sequence above: the norm is still decreasing at level \(7\), and the gap ratios \(\delta_{n+1}/\delta_n\) are increasing but still near \(0.67<1\). A geometric tail with ratio \(0.67\) would leave \(L\) near \(0.82\); a later rise of that ratio toward \(1\) can still push the limit to \(0\). Nothing proved here excludes \(0\).

**Falsifier.** The first level \(n\) with \(\|F(n)\|<0.70\), or the first \(n\ge 7\) with \(\|F(n+1)\|>\|F(n)\|+10^{-8}\). Level \(8\) is the first level not required by the default witness. Through level \(10\) neither falsifier has fired: \(\|F(10)\|\approx 0.833640\), still decreasing.

A second, separate numerical guess, also not a theorem: at \(\alpha=0.10\), levels \(1\) through \(6\) give
\[
\|F\|\approx(1.50450860,\; 1.44384308,\; 1.42364897,\; 1.41577938,\; 1.41266163,\; 1.41143248),
\]
and the successive drops after level \(3\) contract by a factor \(<1/2\). **CONJECTURE J.2.** \(\|F(n,0.10)\|\ge 1.410\) for every \(n\). Falsified by the first level with \(\|F\|<1.410\). The \(\alpha\to 0\) profile \(u\equiv (a+b+c)/3\) off the corners has audit norm \(\sqrt{7/3}\approx 1.527525\), independent of \(n\ge 1\); \(\alpha=0.10\) is a different operator, and this sentence does not identify its limit with \(\sqrt{7/3}\).


## What this does not say

- **NOT A LIMIT THEOREM AT \(\alpha=0.45\).** Levels \(1\) through \(7\) still have increasing ratios, now past \(0.97\). Theorem J excludes \(+\infty\). Theorem K gives \(\|F\|\to 0\) only for \(\alpha>\log 3/\log 5\). Theorem L says energy collapse would force \(\|F\|\to 0\) at every \(\alpha\in(0,1)\), including \(0.45\), but does not prove energy collapse. Neither decides zero versus positive at \(\alpha=0.45\). Conjecture J.1 is not a theorem. Kept Failures J.1, J.2, J.3, and J.4 reject four sufficient conditions that would have closed it. Level \(10\) is recorded, not a limit. Theorem S settles the harmonic audit under uniform multiplicative weights (finite nonzero only at resistance \(\lambda=5/3\)); Kept Failure S.1 rejects inheriting that freeze for \(L^{0.45}\). Theorem W settles the multi-scale harmonic audit (finite for \(\lambda<5/3\); resistance diverges); Kept Failure W.1 rejects a multi-scale resistance freeze. Theorem X takes spectral \(L^{\alpha}\) on the multi-scale graph: neutrality and Schur hold, nested harmonicity is special to \(\alpha=1\); Kept Failures X.1–X.2 reject inheriting the \(\lambda=1\) freeze and reject coincidence with finest-edge fractional.
- **NOT A RESIDUAL.** Their tilt diagnostic sets \((Lu)_i=1\) on the interior. Theorem F then says the corner currents sum to \(-(N-3)\), the source that was inserted. That sum grows like \(3^n\). It is not a force left after the source is removed.
- **FAILED as a continuum thrust for the combinatorial harmonic case.** Under Theorem G, \(\|F(n)\|\to 0\) for every fixed corner triple. The combinatorial harmonic limit does not keep a finite audit vector.
- **OPEN.** Whether \(\|F(n)\|\) for \(L^{0.45}\) tends to \(0\) or to a positive finite limit. Divergence is excluded (Theorem J). Neutrality is settled (Theorem H). For \(\alpha>\log 3/\log 5\) the same audit norm tends to \(0\) (Theorem K). Energy controls the audit norm at every \(\alpha\in(0,1)\) (Theorem L), so energy collapse would force \(\|F\|\to 0\), and a positive audit-norm liminf would force a positive energy liminf; whether either vanishes is OPEN. Geometric weights \(1/d^2\) on this build are settled (Theorem I): that audit norm diverges. Uniform multiplicative weights are settled for the harmonic case (Theorem S): finite nonzero only at the resistance value \(\lambda=5/3\); resistance weights do not freeze the fractional \(\alpha=0.45\) audit (Kept Failure S.1). Multi-scale multiplicative weights are settled for the harmonic case (Theorem W): finite nonzero for every \(\lambda<5/3\), and resistance diverges like \((n+1)\sqrt{21}/2\) (Kept Failure W.1). Spectral \(L^{\alpha}\) on that multi-scale graph is settled structurally by Theorem X (neutrality, Schur, \(\alpha=1\) recovers W, no nested fractional harmonicity); the \(\lambda=1\) freeze does not transfer to \(\alpha=0.45\) (Kept Failure X.1), and multi-scale fractional does not coincide with finest-edge fractional (Kept Failure X.2). None of this is thrust.
- **OPEN.** The Stage 2 question from the rectangle note, whether an informational tensor on the \(0.45\) mesh has nonzero integrated divergence after the Maxwell piece is removed. Theorem F says the harmonic combinatorial gasket does not supply that divergence. Under the Hessian reading of the rank-2 symbol, Theorem M in the README reduces the question to \(\Delta\Psi_{\mathrm{info}}\) on the enclosing surface (zero for static massless \(A_0\) in vacuum, Kept Failure M.1). Under the quadratic reading, Theorem N makes it a volume integral of \(\partial_i\Psi\,\Delta\Psi\), which for static massless \(A_0\) is the Maxwell piece itself and vanishes for an isolated device (Kept Failure N.1). Theorem O extends this to every two-derivative reading: the only survivor is a surface-dependent trace flux whose all-space total is zero (Kept Failure O.1). Theorem P extends it to every linear reading of any order: zero for shift-invariant ones (Kept Failure P.1), and a surface-dependent \(-\int_V E\) for the undifferentiated term (Kept Failure P.2). Theorem Q extends it to every shift-invariant polynomial reading of any order and degree: zero on every vacuum surface if the divergence vanishes on harmonic fields, otherwise surface-dependent (Kept Failure Q.1). Theorem R drops shift invariance and keeps the same dichotomy; the undifferentiated \(\delta^{ij}\Psi^2\) term does not even decay under dilation (Kept Failure R.1). Theorem T takes the same polynomial family under screening and still finds only surface-dependent mass remainders (Kept Failure T.1). Theorem U covers every first-derivative \(O(3)\)-covariant reading \(f(s)\delta+g(s)E\otimes E\) with \(C^1\) coefficients (non-polynomial allowed): conserved ones are Maxwell plus \(\kappa\delta\) and give \(\mathcal{G}=0\); all others are surface-dependent and decay like \(R^{-2}\) (Kept Failure U.1). Theorem V covers the \(\Psi+\nabla\Psi\) isotropic extension \(f(\Psi,s)\delta+g(\Psi,s)E\otimes E\): conserved ones remain Maxwell plus \(\kappa\delta\) with \(\mathcal{G}=0\); all others are surface-dependent (Kept Failure V.1).

## Reproduce

```bash
python3 scripts/gasket_corner_current.py
python3 scripts/gasket_fractional_currents.py
python3 scripts/gasket_geometric_currents.py
python3 scripts/gasket_uniform_weights.py
python3 scripts/gasket_multiscale_weights.py
python3 scripts/gasket_multiscale_fractional.py
python3 scripts/gasket_fractional_limit.py
python3 scripts/gasket_fractional_level9.py
python3 scripts/gasket_fractional_level10.py
```

The first script exits nonzero unless the level-1 triple is \((3/2,-6/5,-3/10)\) and levels 2, 3, 4 match the factor \((3/5)^{n-1}\), including one arbitrary corner triple at levels 1 and 2. The second exits nonzero unless fractional currents at \(\alpha=0.45\) stay neutral through level 4, the level-2 audit norm matches \(1.165909\) within \(5\cdot 10^{-3}\), the refinement ratios differ from \(3/5\) by more than \(0.15\), and the \(\alpha=1\) control recovers Theorem G. The third exits nonzero unless geometric \(1/d^2\) currents match \(4^n\) times the combinatorial ones through level 5, the audit norm equals \((12/5)^n\sqrt{21}/2\), and every successive ratio equals \(12/5\). The fourth exits nonzero unless uniform multiplicative weights recover the closed form \((\lambda\cdot 3/5)^n\sqrt{21}/2\) through level 5 (including the resistance constant \(\sqrt{21}/2\)), and the resistance-scaled fractional audit at \(\alpha=0.45\) is not constant (Kept Failure S.1). The fifth exits nonzero unless multi-scale weights recover \((\sqrt{21}/2)\sum_{k=0}^{n}(\lambda\cdot 3/5)^k\) through level 5, nested harmonicity holds, and resistance multi-scale grows like \((n+1)\sqrt{21}/2\) (Kept Failure W.1). The sixth exits nonzero unless multi-scale spectral \(L^{0.45}\) stays neutral and on the Schur line through level 6, recovers Theorem W at \(\alpha=1\), shows nested harmonicity fails for \(\alpha=0.45\), rejects the \(\lambda=1\) freeze inheritance (Kept Failure X.1), and rejects coincidence with finest-edge fractional (Kept Failure X.2). The seventh exits nonzero unless, through level \(7\), the \(\alpha=0.45\) audit norm matches the level-2 lock \(1.165909\), stays neutral and inside the Theorem J cap, breaks both the \(0.95\) ratio law and the \(5^{-0.45}\) gap contraction, obeys the Theorem L energy sandwich (\(\|F_{\mathrm{frac}}\|\ge w_{\min}\|F\|\) and \(\|F\|\le\kappa\sqrt{Q}\)), breaks the \(0.95\) decay law for \(Q_{0.45}\), and the \(\alpha=0.9\) norm stays under the explicit Theorem K bound. `python3 scripts/gasket_fractional_limit.py 8` recomputes one level higher; the default witness does not. The level-\(9\) and level-\(10\) scripts replay recorded matrix-free values against a live level-\(2\) lock; the level-\(10\) script exits nonzero unless \(\delta_{10}/\delta_9>\log 3/\log 5\) (Kept Failure J.4), the norm still decreases above \(0.70\), and the Schur line, Theorem J cap, and Theorem L sandwich hold. To recompute level \(10\) itself (about half an hour; needs SciPy): `python3 scripts/gasket_fractional_matrix_free.py --level 10`.
