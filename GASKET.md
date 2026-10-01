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
4. **NOT ASSUMED.** Fractional powers \(L^{\alpha}\), a unit source on the interior, geometric weights \(1/d^2\), or any value of \(W\).

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

## What this does not say

- **NOT THIS THEOREM.** The fractional audit \(L^{0.45}\) in `gasket_flux_audit.py`. That operator is nonlocal, so the cell-by-cell extension step does not apply. Their recorded \(\|F\|\approx 1.165909\) at level 2, \(\alpha=0.45\), boundary \((1,-1/2,0)\), is a different number. A float check of that same solve gives corner-current sum at roundoff and \(\|F\|\) matching their test, which is consistent with Theorem F's neutrality but is not the \(3/5\) law.
- **NOT A RESIDUAL.** Their tilt diagnostic sets \((Lu)_i=1\) on the interior. Theorem F then says the corner currents sum to \(-(N-3)\), the source that was inserted. That sum grows like \(3^n\). It is not a force left after the source is removed.
- **FAILED as a continuum thrust.** \(\|F(n)\|\to 0\) for every fixed corner triple. The limit does not keep a finite audit vector.
- **OPEN.** Whether a fractional or weighted operator has a refinement limit with a nonzero neutral dipole. Not decided here.
- **OPEN.** The Stage 2 question from the rectangle note, whether an informational tensor on the \(0.45\) mesh has nonzero integrated divergence after the Maxwell piece is removed. Theorem F says the harmonic combinatorial gasket does not supply that divergence.

## Reproduce

```bash
python3 scripts/gasket_corner_current.py
```

The script exits nonzero unless the level-1 triple is \((3/2,-6/5,-3/10)\) and levels 2, 3, 4 match the factor \((3/5)^{n-1}\), including one arbitrary corner triple at levels 1 and 2.
