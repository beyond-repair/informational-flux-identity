<div align="center">

```
╔══════════════════════════════════════════════════════════════╗
║   ATOMIC DREAM LABS  ·  BEYOND-REPAIR                        ║
╚══════════════════════════════════════════════════════════════╝
```

# Informational Flux Identity

### Closed-surface identity. Signed flux equals summed divergence. No selected W. No thrust.

[![Lifecycle](https://img.shields.io/badge/●_RESEARCH-a855f7?style=for-the-badge&labelColor=0f0f23)](https://github.com/beyond-repair/ADL-Governance)
[![Claim](https://img.shields.io/badge/Claim_≤1-22c55e?style=for-the-badge&labelColor=0f0f23)](https://github.com/beyond-repair/ADL-Governance/blob/main/docs/CLAIM_VALIDATION.md)
[![Governance](https://img.shields.io/badge/ADL--Governance-7c3aed?style=for-the-badge&labelColor=0f0f23)](https://github.com/beyond-repair/ADL-Governance)

```
LIFECYCLE   RESEARCH
CLAIM       ≤ 1   mathematical identity
NOT CLAIMED thrust · continuum limit · selected W
```

</div>

---
## ▌ STATUS

Classification follows [ADL-Governance](https://github.com/beyond-repair/ADL-Governance). A README facelift does not raise claim level. Physics and pharmacology stay at the evidenced cap. CI green is not experimental validation.

---

## ▌ PRESERVED BODY

# Informational flux identity

Closed-surface constraint on the frozen Coherence Drive residual. No value of \(W\) is selected. No thrust is predicted. No continuum limit is taken.

Parent chain, quoted from [coherence-drive/docs/MATH_THEORY_CLOSURE.md](https://github.com/beyond-repair/coherence-drive/blob/main/docs/MATH_THEORY_CLOSURE.md) and left frozen:

\[
T_{\mathrm{eff}}^{ij} = T_{\mathrm{EM}}^{ij} + W(n)\,\chi_{\mathrm{vac}}\,(\nabla\Psi_{\mathrm{info}})^{ij},
\qquad
\mathcal{G}_i = \oint (\nabla\Psi_{\mathrm{info}})^{ij}\, n_j\, dA,
\qquad
\Delta F_i(n) = W(n)\,\chi_{\mathrm{vac}}\,\mathcal{G}_i.
\]

This note does not replace that chain. It says what \(\mathcal{G}_i\) is allowed to be.

```text
experimental_validation     = false
thrust_validated            = false
energy_extraction_validated = false
W_0.08_derived              = false
```

## Object

On a closed rectangle, the signed boundary flux of a rank-2 array, and the condition under which that flux can be a nonzero \(\Delta F\).

## Assumptions

Labeled before use.

1. **ESTABLISHED.** On a finite rectangle partitioned into unit cells, face fluxes are ordinary integers. No continuum measure is assumed.
2. **ESTABLISHED (Riemannian geometry).** On a smooth pseudo-Riemannian manifold with the Levi-Civita connection, the Einstein tensor is divergenceless: \(\nabla^\mu G_{\mu\nu} = 0\). This is the contracted second Bianchi identity. It is not re-derived here.
3. **PROPOSED, their registered model, constant \(W\).** \(G_{\mu\nu} = 8\pi G\,(T_{\mu\nu} + W\, T^{\mathrm{info}}_{\mu\nu})\) with \(W\) a constant scalar, as registered in [CFTv3.3](https://github.com/beyond-repair/CFTv3.3-IQG-Unified-Framework). Not established as physics.
4. **PROPOSED, frozen residual.** \(\Delta F_i\) is defined by the display above. \(\chi_{\mathrm{vac}}\) and \((\nabla\Psi_{\mathrm{info}})^{ij}\) are not derived in this note.
5. **NOT ASSUMED.** \(W = 0.08\), \(W(n) = 0.08\,e^{0.23(n-3)}\), the factor \(10^3\) in the informational fork protocol, entanglement equilibrium, or a laboratory force.

## Theorem A. Rectangle telescoping

**THEOREM.** Let \(F^x_{x,y}\) be given for \(x = 0,\ldots,n_x\) and \(y = 0,\ldots,n_y-1\), and \(F^y_{x,y}\) for \(x = 0,\ldots,n_x-1\) and \(y = 0,\ldots,n_y\). Define the cell divergence

\[
(\mathrm{div}\,F)_{x,y}
  = \bigl(F^x_{x+1,y} - F^x_{x,y}\bigr)
  + \bigl(F^y_{x,y+1} - F^y_{x,y}\bigr)
\]

for \(x = 0,\ldots,n_x-1\) and \(y = 0,\ldots,n_y-1\). Then

\[
\sum_{x,y} (\mathrm{div}\,F)_{x,y}
  = \sum_y F^x_{n_x,y}
  - \sum_y F^x_{0,y}
  + \sum_x F^y_{x,n_y}
  - \sum_x F^y_{x,0}.
\]

The right-hand side is the outward signed flux. The identity is algebraic. It does not depend on symmetry of the array, on which face carries the largest absolute flux, or on any field equation.

**Proof.** Sum the \(F^x\) differences in \(x\). For each fixed \(y\) the sum collapses to \(F^x_{n_x,y} - F^x_{0,y}\). Sum the \(F^y\) differences in \(y\). For each fixed \(x\) the sum collapses to \(F^y_{x,n_y} - F^y_{x,0}\). Adding these is the claim.

**Corollary.** If \((\mathrm{div}\,F)_{x,y} = 0\) at every cell, the outward signed flux is \(0\). An interior rearrangement of a divergence-free array cannot change the outward flux, because every interior face is added once and subtracted once.

## Theorem B. A divergence-free witness with aft-heavy absolute flux

**DERIVED RESULT, from Theorem A plus the explicit array below.** There exists a divergence-free integer flux on the \(16\times 16\) rectangle whose absolute flux on the right face is \(349/366\) of the total absolute boundary flux, while the signed outward flux is \(0\).

Construct \(P_{x,y}\) on the \(17\times 17\) vertex grid by zeros except

\[
\begin{aligned}
P_{15,\cdot} &= (0,-4,-4,3,2,3,0,3,-2,0,3,-3,-2,-3,0,4,-3),\\
P_{16,\cdot} &= (-2,1,16,27,-28,-22,20,27,-15,-11,23,-5,-14,20,-15,-6,9).
\end{aligned}
\]

Set \(F^x_{x,y} = P_{x,y+1} - P_{x,y}\) and \(F^y_{x,y} = P_{x,y} - P_{x+1,y}\).

**Cancellation.** Substitute into the cell divergence:

\[
\begin{aligned}
F^x_{x+1,y} - F^x_{x,y}
  &= \bigl(P_{x+1,y+1}-P_{x+1,y}\bigr) - \bigl(P_{x,y+1}-P_{x,y}\bigr),\\
F^y_{x,y+1} - F^y_{x,y}
  &= \bigl(P_{x,y+1}-P_{x+1,y+1}\bigr) - \bigl(P_{x,y}-P_{x+1,y}\bigr).
\end{aligned}
\]

Each of \(P_{x+1,y+1}\), \(P_{x+1,y}\), \(P_{x,y+1}\), and \(P_{x,y}\) appears once with each sign, so the sum is \(0\). So \(\mathrm{div}\,F = 0\) at every cell, and Theorem A gives signed flux \(0\).

**Witness, computed by `scripts/flux_identity.py`.** Signed outward pieces: left \(0\), right \(+11\), bottom \(-2\), top \(-9\). Sum \(0\). Absolute pieces: left \(0\), right \(349\), bottom \(2\), top \(15\). Sum \(366\). Right-face share \(349/366\).

So a face can hold about \(95\%\) of the absolute flux while the net force, identified with the signed flux, is exactly zero. Localizing absolute stress on an aft face does not produce \(\Delta F\). This is the discrete form of the obstruction in front of the topological-pinch percentage. The percentage \(92\%\) is not used and is not confirmed here. The \(349/366\) figure is a witness that localization and net flux are different quantities, not a prediction for the \(0.45\) mesh.

**Source witness.** Add \(17\) to the single boundary entry \(F^x_{16,3}\). Theorem A still applies to the modified array, which is no longer divergence-free. The script records \(\sum\mathrm{div} = 17\) and signed flux \(17\). A nonzero net appears only when the divergence sum is nonzero, and the two integers agree.

## Theorem C. What a nonzero \(\mathcal{G}\) would mean

**DERIVED RESULT, inside the frozen definition.** Suppose the informational integrand is represented, on a rectangle, by a face-flux array \(F\), and \(\mathcal{G}\) is its outward signed flux. By Theorem A,

\[
\mathcal{G} = \sum_{\mathrm{cells}} \mathrm{div}\,F.
\]

Therefore \(\Delta F_i = W(n)\,\chi_{\mathrm{vac}}\,\mathcal{G}_i\) is zero whenever that informational array is divergence-free, for any constant \(W(n)\) and any constant \(\chi_{\mathrm{vac}}\), including the symbolic values left unfixed by the freeze. Asymmetry of the rectangle is never consulted.

Equivalently: a numerically nonzero \(\mathcal{G}\) on a closed surface around the device is a measurement of \(\sum\mathrm{div}\,F\), not of how the absolute flux is split between faces. Stage 2 can still find \(\sum\mathrm{div}\,F \neq 0\). This note does not forbid that. It says the quantity that has to come out nonzero is the integrated divergence, after the classical Maxwell piece has been subtracted, which is the unresolved item already named in the coherence-drive ledger ("genuine residual force after classical subtraction").

## Theorem D. Constant-\(W\) field equation cancels a lone informational flux

**DERIVED RESULT, from assumptions 2 and 3.** Assume \(G_{\mu\nu} = 8\pi G\,(T_{\mu\nu} + W T^{\mathrm{info}}_{\mu\nu})\) with \(W\) constant, and \(\nabla^\mu G_{\mu\nu} = 0\). Then

\[
\nabla^\mu\bigl(T_{\mu\nu} + W\, T^{\mathrm{info}}_{\mu\nu}\bigr) = 0.
\]

On any region where the integral form of this equation holds, the outward flux of \(T + W T^{\mathrm{info}}\) is zero. A nonzero flux of the informational piece is cancelled by an opposite flux of \(T\). The informational piece alone is not the momentum ledger of an isolated system.

This does not say the informational flux is numerically zero. It says that publishing that flux as a net force, while leaving out the flux of \(T\), is an incomplete ledger. That is the same gate momentum-closure already states as policy ("classical conservation closes first"). Here it is the integral consequence of the registered field equation plus Bianchi.

If assumption 3 is false, the cancellation is not claimed. A theory that adds \(W T^{\mathrm{info}}\) and also drops \(\nabla^\mu G_{\mu\nu} = 0\) is a different theory, and this theorem does not apply to it.

## Theorem E. A scalar \(W(x)\) still reduces to the boundary

**DERIVED RESULT, one-dimensional model of a variable weight.** This is not their frozen \(W(n)\), which is constant at fixed \(n\). It is the smallest variable-weight extension, checked because local \(W(x)\) is listed as open in [-ware-constant-derivation](https://github.com/beyond-repair/-ware-constant-derivation).

Let \(s_0,\ldots,s_{n-1}\) and \(w_0,\ldots,w_n\) be integers. Then

\[
\sum_{i=0}^{n-1} s_i\,(w_{i+1}-w_i)
  = s_{n-1} w_n - s_0 w_0 - \sum_{i=1}^{n-1} (s_i - s_{i-1})\, w_i.
\]

**Proof.** Expand the left-hand side: \(-s_0 w_0 + \sum_{i=1}^{n-1} w_i\,(s_{i-1}-s_i) + s_{n-1} w_n\). Reindex the middle sum. That is the right-hand side.

**Corollaries, still inside this identity.**

- If \(s_i = c\) for every \(i\) and \(w_0 = w_n\), the sum is \(0\). The script checks \(c = 4\) and \(w = (3,8,-1,6,3)\).
- If \(s_i = c\) for every \(i\), the sum equals \(c\,(w_n - w_0)\). The script checks \(c = 4\) and \(w = (3,8,-1,6,9)\), where the sum is \(24\). That integer is the boundary term. It is not an interior source.
- In one dimension, \(s_i - s_{i-1} = 0\) for every \(i\) means \(s\) is constant. The only array that is both constant and zero at the endpoints is \(s = 0\). A nontrivial divergence-free flux with absolute flux piled on one face needs at least two dimensions. That witness is Theorem B, not this line.

A first attempt at the second corollary claimed that a nonconstant array with zero endpoints pairs to zero against an arbitrary \(w\). That claim is false. The script rejected it (the volume pairing was not zero), and it is not part of the theorem. The identity that survived is the summation by parts above.

So a spatially varying weight does not by itself create an interior source. A nonzero sum is either the endpoint term \(c\,(w_n-w_0)\) or a nonzero edge-difference of \(s\), and both sit in the identity.

## Theorem M. Under the Hessian reading, \(\mathcal{G}\) sees only \(\Delta\Psi\) on the surface

**ASSUMPTION M (reading, not a change to the freeze).** The frozen chain writes \((\nabla\Psi_{\mathrm{info}})^{ij}\) without saying which rank-2 object it is. This section takes the Hessian reading \((\nabla\Psi_{\mathrm{info}})^{ij}=\partial^i\partial^j\Psi_{\mathrm{info}}\). The quadratic reading \(\partial^i\Psi\,\partial^j\Psi-\tfrac12\delta^{ij}|\nabla\Psi|^2\) is a different object and is not covered here. Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), or \(\kappa\) enters.

**THEOREM M (continuum).** Let \(V\subset\mathbb{R}^d\) be a bounded Lipschitz domain with boundary \(S\) and outward normal \(n\), and let \(\Psi\) be \(C^3\) on an open neighbourhood \(U\) of \(S\). Nothing is assumed about \(\Psi\) inside \(V\) away from \(U\): sources, singularities, and an asymmetric device are all allowed there. Then
\[
\mathcal{G}_i=\oint_S \partial_i\partial_j\Psi\,n_j\,dA=\oint_S \Delta\Psi\,n_i\,dA .
\]

**Proof.** Take a smooth \(\chi\) with \(\chi=1\) near \(S\) and \(\operatorname{supp}\chi\subset U\), and set \(\Phi=\chi\Psi\), extended by \(0\). Then \(\Phi\) is \(C^3\) on a neighbourhood of \(\bar V\) and agrees with \(\Psi\) to third order on \(S\). The divergence theorem twice, with \(\partial_j\partial_j\partial_i=\partial_i\Delta\), gives \(\oint_S\partial_i\partial_j\Phi\,n_j=\int_V\partial_i\Delta\Phi=\oint_S\Delta\Phi\,n_i\). Replace \(\Phi\) by \(\Psi\) on \(S\).

**Corollaries, inside Assumption M.**

1. If \(\Delta\Psi=0\) on \(S\), then \(\mathcal{G}=0\), for every interior source distribution and every interior asymmetry. The integrated divergence of the Hessian is the boundary integral of \(\Delta\Psi\), so a source inside \(V\) is invisible to \(\mathcal{G}\) unless it reaches the surface.
2. **KEPT FAILURE M.1.** Under C6 (\(\Psi_{\mathrm{info}}=A_0\)) with static massless fields (electrostatics), \(\Delta A_0=-\rho/\varepsilon_0=0\) on any closed surface drawn in vacuum around the device. So \(\mathcal{G}=0\) and \(\Delta F=W\chi_{\mathrm{vac}}\mathcal{G}=0\) for every \(W\) and every \(\kappa\) in C5. The route "Hessian reading plus C5 plus static massless C6" cannot give a nonzero residual on a vacuum surface. This is a failure of that route, not of Stage 2.
3. If \(\Psi\) obeys a screened (Proca-type) equation \((\Delta-m^2)\Psi=0\) on \(S\), then \(\mathcal{G}_i=m^2\oint_S\Psi\,n_i\,dA\). So on a source-free surface where the field equation is screened rather than Laplace, the surviving \(\mathcal{G}\) is the mass term itself. No \(m\) is selected and no number is claimed.

**THEOREM M (discrete, exact on integers).** On the \(N_x\times N_y\) cell rectangle of Theorem A, take any potential \(P\) on the vertices \([-1,N_x+1]\times[-1,N_y+1]\), set \(g^i=D_iP\) (forward difference), and let row \(i\) of the array be the face fluxes of \(g^i\): \(F^x=g^i_{x,y}-g^i_{x-1,y}\), \(F^y=g^i_{x,y}-g^i_{x,y-1}\). With \(\Delta\) the five-point Laplacian,
\[
\mathcal{G}_x=\sum_{y=0}^{N_y-1}\bigl(\Delta P_{N_x,y}-\Delta P_{0,y}\bigr),\qquad
\mathcal{G}_y=\sum_{x=0}^{N_x-1}\bigl(\Delta P_{x,N_y}-\Delta P_{x,0}\bigr).
\]
**Proof.** The cell divergence of row \(i\) is \(\Delta g^i=D_i\Delta P\), because constant-coefficient differences commute. Theorem A turns the signed flux into the cell sum, and the sum of \(D_i\Delta P\) telescopes in direction \(i\).

**Witness, computed by `scripts/hessian_flux_reduction.py`** on an \(8\times 8\) rectangle. The identity holds on \(200\) seeded random integer potentials. A lopsided interior lump (values \(7,11,-2\) in the left half and \(13,5\) in the right half, \(\Delta P=0\) on every boundary layer) has summed absolute cell divergence \(398\) in row \(x\) and \(422\) in row \(y\), and \(\mathcal{G}=(0,0)\) exactly. The potential \(P=x^3\) has \(\Delta P=6x\) and gives \(\mathcal{G}=(384,0)=(6N_xN_y,0)\). A column with \(\Delta P=P=1\) on the right boundary layer gives \(\mathcal{G}_x=8\), the discrete \(m^2\oint\Psi\,n_x\) with \(m^2=1\). The script exits nonzero if any of these integers moves.

**NOT THIS.** Theorem M does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh, does not decide which rank-2 reading the freeze intends, and does not say the quadratic reading gives zero. It narrows the Stage 2 question under the Hessian reading to one quantity: \(\Delta\Psi_{\mathrm{info}}\) on the enclosing surface. A nonzero \(\mathcal{G}\) is not thrust.

## Theorem N. Under the quadratic reading, \(\mathcal{G}\) is a volume integral of \(\partial_i\Psi\,\Delta\Psi\)

**ASSUMPTION N (reading, not a change to the freeze).** This section takes the quadratic reading \((\nabla\Psi_{\mathrm{info}})^{ij}=Q^{ij}:=\partial^i\Psi\,\partial^j\Psi-\tfrac12\delta^{ij}|\nabla\Psi|^2\). It is the second reading named in Theorem M and is a different object from the Hessian. \(\chi_{\mathrm{vac}}\) stays outside the integral, as the frozen chain writes it. Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters.

**THEOREM N (continuum).** For \(\Psi\in C^2\), \(\partial_jQ^{ij}=\partial_i\Psi\,\Delta\Psi\). So for a bounded Lipschitz domain \(V\subset\mathbb{R}^d\) with boundary \(S\) and \(\Psi\in C^2(\bar V)\),
\[
\mathcal{G}_i=\oint_S Q^{ij}n_j\,dA=\int_V \partial_i\Psi\,\Delta\Psi\,dV .
\]
If \(\Delta\Psi=0\) on the shell between two nested closed surfaces, \(\mathcal{G}\) is the same on both.

**Proof.** \(\partial_j(\partial_i\Psi\,\partial_j\Psi)=\partial_i\partial_j\Psi\,\partial_j\Psi+\partial_i\Psi\,\Delta\Psi\) and \(\partial_i(\tfrac12|\nabla\Psi|^2)=\partial_j\Psi\,\partial_i\partial_j\Psi\). Subtract, then apply the divergence theorem; on the shell the integrand is \(0\).

**Corollaries, inside Assumption N.**

1. **KEPT FAILURE N.1.** Under C6 (\(\Psi_{\mathrm{info}}=A_0\)) with static massless fields, \(\partial_iA_0=-E_i\), so \(\varepsilon_0Q^{ij}=\varepsilon_0(E^iE^j-\tfrac12\delta^{ij}|E|^2)\), which is exactly the Maxwell electrostatic stress. Under this route the informational tensor *is* the Maxwell piece divided by \(\varepsilon_0\); after the Maxwell piece is removed, nothing is left, identically. Without removing it, \(\Delta A_0=-\rho/\varepsilon_0\) gives \(\mathcal{G}_i=\varepsilon_0^{-1}\int_V\rho\,E_i\,dV\), the ordinary Coulomb force on the enclosed charge. For an isolated device (bounded, compactly supported, Hölder \(\rho\) inside \(S\), no charge outside), that self-force is \(\varepsilon_0^{-1}\iint\rho(x)\rho(y)\,(x-y)/(4\pi\varepsilon_0|x-y|^3)\,dx\,dy=0\): the integrand is absolutely integrable and odd under \(x\leftrightarrow y\). So \(\mathcal{G}=0\) and \(\Delta F=0\) for every \(W\) and every \(\kappa\), for every interior asymmetry. With charges outside \(S\), \(\mathcal{G}\) is the classical Coulomb force on the inside, whose reaction sits on the outside charges. The route "quadratic reading plus C5 plus static massless C6" cannot give a residual. This is a failure of that route, not of Stage 2.
2. **Screened (Proca-type) equation.** If \((\Delta-m^2)\Psi=-\rho/\varepsilon_0\), then \(\partial_i\Psi\,\Delta\Psi=\rho E_i/\varepsilon_0+\tfrac12m^2\partial_i(\Psi^2)\), so
\[
\mathcal{G}_i=\varepsilon_0^{-1}\int_V\rho\,E_i\,dV+\tfrac{m^2}{2}\oint_S\Psi^2\,n_i\,dA .
\]
The tensor \(Q^{ij}-\tfrac12m^2\Psi^2\delta^{ij}\) is divergence-free wherever \((\Delta-m^2)\Psi=0\). For an isolated source the Yukawa kernel is radial, so the self-force term is \(0\) by the same odd-integrand argument, and \(\mathcal{G}_i=\tfrac{m^2}{2}\oint_S\Psi^2n_i\,dA\). That term is nonzero for an asymmetric source, depends on which surface is drawn, and is \(O(e^{-2mR})\) on a sphere of radius \(R\) around the source, so it tends to \(0\) as the surface recedes. It is the flux of the mass term that the quadratic reading leaves out, not an interior source. No \(m\) is selected and no number is claimed.

**Witness, computed by `scripts/quadratic_flux_reduction.py`.** Exact (rational arithmetic) on the box \([0,2]\times[-1,1]\times[0,3]\): \(\oint Q\,n=\int\partial_i\Psi\,\Delta\Psi\) for \(60/60\) seeded random cubic integer polynomials, and \(\Psi=x^3\) gives \(\mathcal{G}_x=432\). The harmonic \(\Psi=x^3-3xy^2+2yz\) gives \(\mathcal{G}=(0,0,0)\) exactly while row \(x\) carries \(907/5\) and \(173/5\) on the two \(x\)-faces and \(-90\), \(-126\) on the \(y\)-faces (total absolute \(432\)). Quadrature (\(96\) Gauss–Legendre by \(192\) azimuthal nodes, \(4\pi\varepsilon_0=1\)) on the sphere \(R=2\) around the asymmetric point-charge cluster \(q=(3,-1,2)\): \(|\mathcal{G}|=3.0\times10^{-14}\) against \(\oint|Q^{xj}n_j|=15.131871\). Point charges are allowed here because the sphere lies in vacuum and replacing each charge by a small uniform ball leaves the exterior field unchanged. Adding a charge \(1.5\) at \((4,1,-0.5)\), outside the sphere, gives \(\mathcal{G}=4\pi F_{\mathrm{Coulomb}}\) to \(10^{-10}\) relative, with \(\mathcal{G}_x=-2.9694050402\). With Yukawa potentials at \(m=0.7\) and no outside charge, \(\mathcal{G}=\tfrac{m^2}{2}\oint\Psi^2n\) to \(10^{-10}\) relative at \(R=2,3,5\), with \(|\mathcal{G}|=1.524401,\ 0.3233750,\ 0.01713515\). The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem N does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which rank-2 reading the freeze intends. Together with Theorem M it covers the two readings named so far; under both, static massless C6 on a vacuum surface gives \(\mathcal{G}=0\) for an isolated device. A nonzero \(\mathcal{G}\) under screening is a surface-dependent mass-term flux. None of this is thrust.

## Theorem O. Every two-derivative reading reduces to a trace flux

**ASSUMPTION O (family of readings, not a change to the freeze).** Theorems M and N each took one reading of the rank-2 symbol \((\nabla\Psi_{\mathrm{info}})^{ij}\). This section takes all of them at once inside one class: \(T^{ij}\) is a local polynomial in the derivatives of \(\Psi\) at the same point, with constant coefficients, unchanged by \(\Psi\to\Psi+\text{const}\), covariant under rotations and reflections (\(O(d)\)), and carrying exactly two derivatives in each term. \(\chi_{\mathrm{vac}}\) stays outside the integral. Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters.

**THEOREM O (classification and divergence).** Every such \(T\) has the form
\[
T^{ij}=a\,\partial_i\partial_j\Psi+b\,\delta_{ij}\Delta\Psi+c\,\partial_i\Psi\,\partial_j\Psi+d\,\delta_{ij}|\nabla\Psi|^2,
\]
and for \(\Psi\in C^3\),
\[
\partial_jT^{ij}=(a+b)\,\partial_i\Delta\Psi+c\,\partial_i\Psi\,\Delta\Psi+\Bigl(\tfrac c2+d\Bigr)\partial_i|\nabla\Psi|^2 .
\]
The Hessian reading is \((1,0,0,0)\) and the quadratic reading of Theorem N is \((0,0,1,-\tfrac12)\).

**Proof.** A polynomial unchanged by every shift of \(\Psi\) does not depend on undifferentiated \(\Psi\), so each factor carries at least one derivative. With two derivatives in total, a term is either one factor \(\partial_k\partial_l\Psi\) or two factors \(\partial_k\Psi\,\partial_l\Psi\), contracted to the free indices \(i,j\) by an \(O(d)\)-invariant rank-4 tensor. Those are spanned by \(\delta_{ij}\delta_{kl}\), \(\delta_{ik}\delta_{jl}\), \(\delta_{il}\delta_{jk}\) (first fundamental theorem for \(O(d)\)), which gives the four terms. For the divergence, \(\partial_j\partial_i\partial_j\Psi=\partial_i\Delta\Psi\) and \(\partial_j(\partial_i\Psi\,\partial_j\Psi)=\tfrac12\partial_i|\nabla\Psi|^2+\partial_i\Psi\,\Delta\Psi\).

**Corollaries, inside Assumption O, \(d=3\).** Take C6 (\(\Psi_{\mathrm{info}}=A_0\)) with static massless fields and an isolated device: \(\rho\in C^\infty_c\), \(\varepsilon_0\Delta A_0=-\rho\), \(A_0\to0\) at infinity, and \(S\) a closed Lipschitz surface in vacuum enclosing \(\operatorname{supp}\rho\). Write \(E=-\nabla A_0\).

1. **THEOREM O.1 (trace flux).** \(\mathcal{G}_i(S)=\bigl(\tfrac c2+d\bigr)\oint_S|E|^2n_i\,dA\). *Proof.* The \(a\) and \(b\) terms give \((a+b)\oint_S\Delta A_0\,n_i=0\) (Theorem M; \(\Delta A_0=0\) on a vacuum surface). Split \(c\,\partial_i\Psi\partial_j\Psi+d\,\delta_{ij}|\nabla\Psi|^2=c\,Q^{ij}+(\tfrac c2+d)\delta_{ij}|\nabla\Psi|^2\); the \(Q\) flux is the Coulomb self-force over \(\varepsilon_0\), which is \(0\) (Kept Failure N.1). The Maxwell piece \(\varepsilon_0Q\) contributes \(0\) on every such \(S\), so removing it, with any coefficient, leaves \(\mathcal{G}\) unchanged.
2. **THEOREM O.2 (zero over all space).** \(\partial_jT^{ij}\) is absolutely integrable on \(\mathbb{R}^3\) and \(\int_{\mathbb{R}^3}\partial_jT^{ij}\,dV=0\). On every such \(S\), the integrated divergence inside equals minus the integrated divergence outside: \(\mathcal{G}_i(S)=-\int_{\mathbb{R}^3\setminus V}\partial_jT^{ij}\,dV\). *Proof.* The \(a+b\) and \(c\,\partial_i\Psi\Delta\Psi\) terms have compact support; outside the device the remaining term is \((\tfrac c2+d)\partial_i|E|^2=O(r^{-5})\). On the sphere \(S_R\), \(|E|=O(R^{-2})\), so \(|\mathcal{G}(S_R)|=O(R^{-2})\to0\).
3. **Asymptotics.** For the dilations \(S_R=R\,S_1\) of a fixed vacuum surface \(S_1\) around a point inside the device, \(E=Q_{\mathrm{tot}}x/(4\pi\varepsilon_0|x|^3)+O(|x|^{-3})\) gives \(\mathcal{G}_i(S_R)=(\tfrac c2+d)\,R^{-2}\,\tfrac{Q_{\mathrm{tot}}^2}{(4\pi\varepsilon_0)^2}\oint_{S_1}n_i\,|y|^{-4}\,dA_y+O(R^{-3})\). The leading coefficient vanishes on spheres centred at the expansion point and is generally nonzero otherwise.
4. **KEPT FAILURE O.1.** No reading in the family, combined with C5 and static massless C6, gives a net source for an isolated device. If \(c+2d=0\) (this includes the Hessian, \(\partial\partial\Psi-\delta\Delta\Psi\), and \(Q\)), \(\mathcal{G}=0\) on every vacuum surface. If \(c+2d\neq0\) (for example \(\partial_i\Psi\partial_j\Psi\) alone, its trace-free part, or \(\delta_{ij}|\nabla\Psi|^2\)), a finite surface does give \(\mathcal{G}\neq0\), but the value changes when the surface is moved around the same device, it is exactly the opposite of the exterior's integrated divergence (Theorem O.2), and it decays like \(R^{-2}\). So \(\Delta F=W\chi_{\mathrm{vac}}\mathcal{G}\) would depend on where the audit surface is drawn. Interpretive step, labeled as such: a force on the device is a property of the device and its field, not of an audit surface, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This is a failure of the route "two-derivative reading plus C5 plus static massless C6", not of Stage 2.

**Witness, computed by `scripts/general_reading_flux.py`.** The linear system for \(O(3)\)-invariant rank-4 tensors (two generic rotations, \(-I\), and one reflection) has a null space of dimension \(3\). Exact (rational arithmetic) on the box \([0,2]\times[-1,1]\times[0,3]\): the divergence formula holds as a polynomial identity, and \(\oint T\,n\) equals the volume integral of it, for \(50/50\) seeded random integer \((a,b,c,d)\) and cubic integer \(\Psi\). Quadrature (\(128\) Gauss–Legendre by \(256\) azimuthal nodes, \(4\pi\varepsilon_0=1\)) around the asymmetric cluster \(q=(3,-1,2)\) of Theorem N, on three vacuum spheres (centre \(0\), \(R=2\); centre \((0.6,0.2,-0.3)\), \(R=2.4\); centre \((-0.9,0,0.5)\), \(R=3\)): Hessian flux at most \(1.6\times10^{-13}\), Maxwell flux at most \(9.3\times10^{-14}\), and \(\mathcal{G}=(\tfrac c2+d)P\) for six readings, with trace flux \(P=\oint|\nabla\Psi|^2n\) equal to \((-18.295274824,\,-6.781652669,\,6.231252988)\), \((-35.04687982,\,-11.166762011,\,14.456027559)\), and \((4.23976535,\,-2.825807999,\,-3.116738754)\). Same device, three surfaces, three different values. Under dilation of the unit sphere centred at \((0.3,-0.2,0.1)\), \(R^2P\) at \(R=4,8,16,32,64\) approaches \(Q_{\mathrm{tot}}^2\oint_{S_1}n/|y|^4=(-105.632225,\,70.421484,\,-35.210742)\), with errors halving each doubling; the Richardson value \(2s(64)-s(32)\) is within \(1.2\times10^{-3}\) relative. The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem O does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh, does not decide which reading the freeze intends, and does not cover readings with more than two derivatives, position-dependent coefficients, undifferentiated \(\Psi\) (such as the screened mass term of Theorem N), or any \(\Psi_{\mathrm{info}}\neq A_0\). A nonzero \(\mathcal{G}\) on a finite surface is not thrust.

## Theorem P. Every linear reading, of any order, reduces to one surface scalar

**ASSUMPTION P (family of readings, not a change to the freeze).** Theorem O stops at two derivatives. This section takes every reading that is linear in \(\Psi\), with any finite number of derivatives:
\[
T^{ij}=p(\Delta)\,\partial^i\partial^j\Psi+\delta^{ij}\,q(\Delta)\,\Psi,\qquad p(s)=\sum_k p_ks^k,\quad q(s)=\sum_k q_ks^k,
\]
with real constant coefficients. The Hessian is \(p=1,q=0\); Theorem O's linear part is \(p=a,\ q=bs\). The reading is shift-invariant exactly when \(q(0)=0\). Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters.

**DERIVED (why this is the whole linear family).** A local, constant-coefficient, linear reading is \(T^{ij}=\sum_m A^{ij}{}_{k_1\cdots k_m}\partial_{k_1}\cdots\partial_{k_m}\Psi\), and \(A\) may be taken symmetric in the \(k\)'s. Covariance under \(O(d)\) makes each \(A\) an isotropic tensor, hence a combination of products of Kronecker deltas, and \(-I\in O(d)\) kills odd ranks. Contracting a product of deltas with symmetric derivatives leaves only \(\delta^{ij}\Delta^{m/2}\Psi\) and \(\partial^i\partial^j\Delta^{m/2-1}\Psi\). The witness checks the dimension count \(1,0,2,0,2\) for \(m=0,\dots,4\) in \(d=3\).

**THEOREM P (continuum).** Let \(V\subset\mathbb{R}^d\) be a bounded Lipschitz domain with boundary \(S\) and outward normal \(n\), let \(D=\max(\deg p,\deg q)\), and let \(\Psi\) be \(C^{2D+3}\) on an open neighbourhood \(U\) of \(S\). Nothing is assumed about \(\Psi\) inside \(V\) away from \(U\). With \(r(s)=s\,p(s)+q(s)\),
\[
\partial_jT^{ij}=\partial_i\,r(\Delta)\Psi\ \text{ on } U,\qquad
\mathcal{G}_i=\oint_S T^{ij}n_j\,dA=\oint_S r(\Delta)\Psi\;n_i\,dA .
\]

**Proof.** Constant-coefficient operators commute, so \(p(\Delta)\partial_i\partial_j\Psi=\partial_i\partial_j\Phi\) with \(\Phi=p(\Delta)\Psi\in C^3(U)\). Theorem M applied to \(\Phi\) gives \(\oint_S\partial_i\partial_j\Phi\,n_j=\oint_S\Delta\Phi\,n_i=\oint_S\Delta p(\Delta)\Psi\,n_i\). The \(\delta\) term contributes \(\oint_S q(\Delta)\Psi\,n_i\). Add. The divergence formula is the same commutation.

**Corollaries, inside Assumption P.**

1. If \(\Delta\Psi=0\) on \(U\), then \(\Delta^k\Psi=0\) there for every \(k\ge1\), so \(\mathcal{G}_i=q(0)\oint_S\Psi\,n_i\,dA\).
2. **KEPT FAILURE P.1.** Every shift-invariant linear reading (\(q(0)=0\)), of any order, combined with C5 and static massless C6 (\(\Psi=A_0\), \(\Delta A_0=0\) in vacuum), gives \(\mathcal{G}=0\) on every vacuum surface, for every interior source distribution, every \(W\), and every \(\kappa\). This extends Kept Failure M.1 from the Hessian to all linear readings. Adding derivatives does not open a residual.
3. **Screened equation.** If \((\Delta-m^2)\Psi=0\) on \(U\), then \(\Delta^k\Psi=m^{2k}\Psi\) and \(\mathcal{G}_i=r(m^2)\oint_S\Psi\,n_i\,dA\): the screened surface term of Theorem M, multiplied by \(r(m^2)/m^2\). No \(m\) is selected and no number is claimed.
4. **KEPT FAILURE P.2 (the undifferentiated term).** If \(q(0)=e\neq0\), static massless C6 gives \(\mathcal{G}_i=e\oint_S A_0\,n_i\,dA=-e\int_V E_i\,dV\). On a sphere centred at \(x_0\) that encloses the whole device, only the dipole term of the exterior multipole expansion survives the integral, so \(\mathcal{G}=e\,p(x_0)/(3\varepsilon_0)\), where \(p(x_0)=\int\rho(x)(x-x_0)\,dx=p-Q_{\mathrm{tot}}x_0\). For a charged device this moves with the centre, and under dilation of an off-centre sphere it grows linearly in \(R\). For a neutral device every enclosing sphere gives the same \(e\,p/(3\varepsilon_0)\), but an ellipsoid around the same charges gives a different value (the classical shape dependence of \(\int E\) over a region that contains the charges). The tensor \(T\) itself also changes by \(eC\delta^{ij}\) under the gauge shift \(\Psi\mapsto\Psi+C\), even though \(\mathcal{G}\) does not. Interpretive step, labeled as such and the same as in Kept Failure O.1: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This is a failure of the route "linear reading plus C5 plus static massless C6", not of Stage 2.

**Witness, computed by `scripts/linear_reading_flux.py`.** The linear system for \(O(3)\)-invariant tensors \(A^{ij}{}_{k_1\cdots k_m}\), symmetric in the \(k\)'s, has null spaces of dimension \(1,0,2,0,2\) for \(m=0,\dots,4\). Exact (rational arithmetic) on the box \([0,2]\times[-1,1]\times[0,3]\): the divergence formula holds as a polynomial identity, and \(\oint T\,n=\oint r(\Delta)\Psi\,n\), for \(40/40\) seeded random integer \(p,q\) of degree at most \(2\) (up to six derivatives) and integer \(\Psi\) of degree at most \(7\). On the harmonic \(\Psi=x^3-3xy^2+2yz+x^5-10x^3y^2+5xy^4\), \(20\) random shift-invariant readings give \(\mathcal{G}=(0,0,0)\) exactly, while row \(x\) of the reading \(p=1+2s-3s^2\), \(q=s+4s^2\) carries \(1584\) in absolute face flux; the undifferentiated reading \(\delta^{ij}\Psi\) gives \(\mathcal{G}=(80,36,0)\). Quadrature (\(128\) Gauss–Legendre by \(256\) azimuthal nodes, \(4\pi\varepsilon_0=1\), so \(\oint\Psi n=\tfrac{4\pi}{3}p(x_0)\) on enclosing spheres): for the charged cluster \(q=(3,-1,2)\) of Theorems N and O, the three spheres of Theorem O give \(\oint\Psi n=(-8.796459430,\,-3.979350695,\,3.351032164)\), \((-18.849555922,\,-7.330382858,\,8.377580410)\), and \((6.283185307,\,-3.979350695,\,-5.026548246)\); under dilation of the unit sphere centred at \((0.3,-0.2,0.1)\) the magnitude is \(30.584623,\,54.991733,\,104.759332,\,204.860214\) at \(R=4,8,16,32\). For the neutral cluster \(q=(3,-1,-2)\) at the same positions, every enclosing sphere gives \((-10.471975512,\,3.560471674,\,-1.675516082)\), the ellipsoid with semi-axes \((2,3,4)\) at the origin gives \((-15.196767609,\,3.257898425,\,-1.061936758)\) (confirmed independently, to \(10^{-9}\), as the sphere value minus a radial volume integral of \(\nabla\Psi\) over the shell between them), and the ellipsoid \((3,1.6,1.8)\) at \((0.1,0,-0.2)\) gives \((-6.204064424,\,4.570118967,\,-1.883253718)\). Same device, three surfaces, three values. The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem P does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends. It does not cover nonlinear readings with more than two derivatives (for example \(\partial^i\Psi\,\partial^j\Delta\Psi\) or quartic terms), position-dependent coefficients, time-dependent fields, or any \(\Psi_{\mathrm{info}}\neq A_0\). A nonzero \(\mathcal{G}\) on a finite surface is not thrust.

## Theorem Q. Every shift-invariant polynomial reading, of any order and degree, is zero or surface-dependent

**ASSUMPTION Q (family of readings, not a change to the freeze).** Theorem O stops at two derivatives and Theorem P at linear readings. This section takes every reading that is a finite sum of monomials, each a product of \(m\ge1\) factors \(\partial_{k_1}\cdots\partial_{k_s}\Psi\) with every \(s\ge1\) (so the reading is unchanged by \(\Psi\to\Psi+\text{const}\)), contracted to \(i,j\) by constant \(O(3)\)-invariant tensors. Write \(D\) for the total number of derivatives in a monomial. Any \(m\), any \(D\). Examples: \(\partial^i\Psi\,\partial^j\Delta\Psi\) (named in Theorem P), \(|\nabla\Psi|^2\partial^i\Psi\partial^j\Psi\), \(|\nabla\Psi|^2\partial^i\partial^j\Psi\). Theorems M, N (massless), O, and P with \(q(0)=0\) are special cases. Take C6 (\(\Psi_{\mathrm{info}}=A_0\)) with static massless fields and an isolated device: all charge inside an open ball \(B\), \(\varepsilon_0\Delta A_0=-\rho\), \(A_0\to0\) at infinity, \(d=3\). A vacuum surface is a closed Lipschitz surface outside \(\bar B\) enclosing it. Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters.

Call the reading **on-shell conserved** if \(\partial_jT^{ij}[h]=0\) for every harmonic polynomial \(h\). Because the coefficients are constant and the \(K\)-th Taylor polynomial of a harmonic function is harmonic, this is the same as \(\partial_jT^{ij}=0\) wherever \(\Delta\Psi=0\).

**THEOREM Q.** For every such reading and every vacuum surface \(S\) bounding \(V\):

1. **Q.1 (exterior formula).** \(\partial_jT^{ij}\) is absolutely integrable outside \(V\), \(\mathcal{G}(S_R)\to0\) on spheres \(S_R\) of radius \(R\to\infty\), and
\[
\mathcal{G}_i(S)=-\int_{\mathbb{R}^3\setminus V}\partial_jT^{ij}\,dV .
\]
2. **Q.2.** If the reading is on-shell conserved, \(\mathcal{G}(S)=0\) on every vacuum surface, for every interior asymmetry, every \(W\), and every \(\kappa\).
3. **Q.3.** If it is not, there are an isolated device made of finitely many point charges in \(B\) and two vacuum surfaces \(S_1,S_2\) around it with \(\mathcal{G}(S_1)\neq\mathcal{G}(S_2)\).

So no reading in the family gives a nonzero \(\mathcal{G}\) that is the same on every audit surface around the device. Removing the Maxwell piece changes nothing: \(\varepsilon_0Q\) is itself an on-shell conserved member (\(m=2\), \(D=2\)), and its flux on every vacuum surface is \(0\) (Kept Failure N.1).

**Proof of Q.1.** Outside a ball containing \(\bar B\), every \(s\)-th derivative of \(A_0\) is \(O(r^{-(s+1)})\) (exterior multipole expansion), so a monomial with \(m\) factors and \(D\) derivatives is \(O(r^{-(D+m)})\) and its divergence is \(O(r^{-(D+m+1)})\). If \(m=1\), the coefficient is an \(O(3)\)-invariant tensor of rank \(D+2\); \(-I\in O(3)\) kills odd ranks and \(D\ge1\), so \(D\ge2\). If \(m\ge2\), \(D\ge m\ge2\). Either way \(D+m\ge3\). So the flux through \(S_R\) is \(O(R^{2-D-m})\to0\) and the divergence is integrable at infinity. The divergence theorem on the region between \(S\) and \(S_R\) (where \(A_0\) is smooth) gives \(\mathcal{G}(S_R)-\mathcal{G}(S)=\int\partial_jT^{ij}\); let \(R\to\infty\).

**Proof of Q.2.** \(\Delta A_0=0\) outside \(V\), so the integrand in Q.1 is \(0\).

**Lemma Q.0 (point charges reach every harmonic jet).** Fix \(y\notin\bar B\) and \(K\ge0\). Every \(K\)-jet at \(y\) of a harmonic function is the \(K\)-jet at \(y\) of \(\sum_a c_a|z-x_a|^{-1}\) for finitely many \(x_a\in B\). *Proof.* Let \(H_K\) be the harmonic polynomials of degree at most \(K\) (the harmonic \(K\)-jets). The jet of \(|z-x|^{-1}\) at \(y\) lies in \(H_K\). If these jets, over \(x\in B\), spanned a proper subspace, some linear functional \(\ell\neq0\) on \(H_K\) would vanish on all of them. Extend \(\ell\) to all polynomials of degree at most \(K\); it has the form \(f\mapsto(P(\partial)f)(y)\) for a polynomial \(P\) of degree at most \(K\). Then \(w\mapsto(P(\partial)|\cdot|^{-1})(w)\) vanishes for \(w\in y-B\), is real-analytic on the connected set \(\mathbb{R}^3\setminus\{0\}\), and so vanishes there. The distribution \(P(\partial)|w|^{-1}\) is then supported at \(0\), so its Fourier transform \(P(i\xi)\,4\pi/|\xi|^2\) is a polynomial, so \(|\xi|^2\) divides \(P(i\xi)\) and \(P(\partial)=P_1(\partial)\Delta\). Then \(\ell(h)=(P_1(\partial)\Delta h)(y)=0\) for every \(h\in H_K\), a contradiction.

**Proof of Q.3.** Let \(K\) be one more than the highest derivative order in \(T\), so \(\partial_jT^{ij}\) at a point depends only on the \(K\)-jet of \(\Psi\) there. Not on-shell conserved gives a harmonic \(h\) and a component \(i\) with \((\partial_jT^{ij}[h])_i\neq0\) at some point; by constant coefficients, move that point to any \(y\notin\bar B\). Lemma Q.0 gives point charges in \(B\) whose potential has the same \(K\)-jet at \(y\), so \(c=\partial_jT^{ij}[A_0](y)\neq0\), and by continuity this component has the sign of \(c\) on a ball \(U\ni y\) with \(\bar U\cap\bar B=\emptyset\). Let \(S_1\) be the sphere about the centre of \(B\) through \(y\), and \(S_2\) the same surface pushed outward by a small smooth bump supported in a cap around \(y\), so the region between them is a nonempty open set inside \(U\). Both are vacuum surfaces around the device, and \(\mathcal{G}_i(S_2)-\mathcal{G}_i(S_1)\) is the integral of that component over the region between them, which is not \(0\).

**KEPT FAILURE Q.1.** No shift-invariant polynomial reading, of any order and any degree, combined with C5 and static massless C6, gives a net source for an isolated device. Either the reading is on-shell conserved and \(\mathcal{G}=0\) on every vacuum surface, or \(\mathcal{G}\) changes when the surface is moved around the same device, equals minus the exterior's integrated divergence (Q.1), and falls off at least like \(R^{2-D-m}\). Interpretive step, labeled as such and the same as in Kept Failures O.1 and P.2: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This closes the nonlinear higher-derivative gap left by Theorems O and P inside static massless C6. It is a failure of that route, not of Stage 2.

**Witness, computed by `scripts/nonlinear_reading_flux.py`.** Four readings: \(T_1=\partial_i\partial_k\Psi\,\partial_j\partial_k\Psi-\tfrac12\delta_{ij}|\partial\partial\Psi|^2\) (the Maxwell stress of \(\nabla\Psi\); \(m=2\), \(D=4\)), \(T_2=\partial_i\Psi\,\partial_j\Delta\Psi\), \(T_3=|\nabla\Psi|^2\partial_i\Psi\partial_j\Psi\) (\(m=4\), \(D=4\)), and \(T_4=|\nabla\Psi|^2\partial_i\partial_j\Psi\) (\(m=3\), \(D=4\)). Exact (rational arithmetic): the four divergence formulas hold as polynomial identities for \(25/25\) seeded random integer \(\Psi\) of degree at most \(5\); on the harmonic \(\Psi_h=x^3-3xy^2+2yz+xyz+x^4-6x^2y^2+y^4\), \(\partial_jT_1^{ij}\) and \(\partial_jT_2^{ij}\) vanish identically, while at \((1,\tfrac12,-1)\) \(\partial_jT_3^{ij}=(\tfrac{858843}{32},\tfrac{327369}{8},-\tfrac{970097}{128})\) and \(\partial_jT_4^{ij}=(\tfrac{19927}{8},-\tfrac{45683}{4},1024)\). Lemma Q.0, exact: the jets of \(|z-x|^{-1}\) at \(y=(3,\tfrac12,-1)\) for rational \(x\) in the unit ball are harmonic Taylor polynomials and have rank \(16\) for \(K=3\) and \(25\) for \(K=4\), the dimensions of \(H_3\) and \(H_4\). Quadrature (\(128\) Gauss–Legendre by \(256\) azimuthal nodes, \(4\pi\varepsilon_0=1\)) around the cluster \(q=(3,-1,2)\) of Theorems N and O, on the three spheres of Theorem O: \(T_1\) gives \(|\mathcal{G}|\le4.0\times10^{-14}\) against absolute row-\(x\) flux up to \(12.120152\); \(T_3\) gives \((-37.0759432675,\,-9.2196731911,\,10.0472161178)\), \((-95.5534175465,\,-21.3721686011,\,33.8903975982)\), \((1.7271019689,\,-1.2018086361,\,-1.0826483846)\), and \(T_4\) gives \((-21.7317583295,\,-5.4809207986,\,5.9790206965)\), \((-50.8676131315,\,-11.4164465132,\,18.1281743757)\), \((1.6845711842,\,-1.1109562809,\,-1.0556329878)\). Same device, three surfaces, three values. Each agrees with minus the exterior integral of the on-shell divergence (radial Gauss–Legendre in \(R/r\), \(96\) nodes) to \(1.1\times10^{-14}\) relative. Under dilation of the unit sphere centred at \((0.3,-0.2,0.1)\), \(R^6\mathcal{G}_{T_3}\) and \(R^5\mathcal{G}_{T_4}\) at \(R=4,\dots,64\) approach the monopole coefficients \((-5438.951969,\,3625.967979,\,-1812.98399)\) and \((-1443.514392,\,962.342928,\,-481.171464)\), errors roughly halving per doubling, Richardson values within \(5\times10^{-3}\) relative. The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem Q does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends. It does not cover readings with undifferentiated \(\Psi\) beyond Theorem P's linear term (for example \(\Psi^2\delta^{ij}\) or a screened mass term), non-polynomial functions of the derivatives, position-dependent coefficients, time-dependent or radiating fields, \(d\neq3\), or any \(\Psi_{\mathrm{info}}\neq A_0\). A nonzero \(\mathcal{G}\) on a finite surface is not thrust.

## Theorem R. Dropping shift invariance: every polynomial reading is still zero or surface-dependent

**ASSUMPTION R (family of readings, not a change to the freeze).** Theorem Q required every factor to carry a derivative, so the reading ignores \(\Psi\to\Psi+\text{const}\). Here that is dropped. A reading is a finite sum of monomials, each a product of \(m\ge0\) factors \(\partial_{k_1}\cdots\partial_{k_s}\Psi\) with every \(s\ge0\), contracted to \(i,j\) by constant \(O(3)\)-invariant tensors; \(D\) is the total number of derivatives in a monomial. Undifferentiated factors are allowed. Examples: \(\delta^{ij}\Psi^2\) (the shape of a mass term, read here under the massless equation), \(\Psi\,\partial^i\partial^j\Psi-\tfrac12\delta^{ij}|\nabla\Psi|^2\), \(\Psi^2\partial^i\partial^j\Psi\). Theorem Q and Theorem P's undifferentiated term \(e\,\delta^{ij}\Psi\) are special cases. C6 (\(\Psi_{\mathrm{info}}=A_0\)), static massless fields, an isolated device inside an open ball \(B\), \(A_0\to0\) at infinity, \(d=3\), and vacuum surfaces are exactly as in Assumption Q. On-shell conserved means \(\partial_jT^{ij}[h]=0\) for every harmonic \(h\), as in Theorem Q. Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters.

Write \(T=\sum_{m,D}T_{m,D}\) for the pieces with \(m\) factors and \(D\) derivatives. The coefficient tensor of a monomial has rank \(D+2\), and \(-I\in O(3)\) kills odd ranks, so \(D\) is even. The only pieces with \(D+m\le2\) are then the constant \(c\,\delta^{ij}\) (\(m=0\)), the linear term \(e\,\delta^{ij}\Psi\) (\(m=1\), \(D=0\)), and the quadratic term \(f\,\delta^{ij}\Psi^2\) (\(m=2\), \(D=0\)). Every other piece has \(D+m\ge3\).

**Lemma R.0 (conservation is graded).** \(T\) is on-shell conserved if and only if every \(T_{m,D}\) is. *Proof.* If \(h\) is harmonic, so is \(\lambda h(\mu x)\) for every \(\lambda\) and \(\mu>0\). With constant coefficients, \(\partial_jT^{ij}[\lambda h(\mu\,\cdot)](x)=\sum_{m,D}\lambda^m\mu^{D+1}\,(\partial_jT_{m,D}^{ij}[h])(\mu x)\). Evaluate at \(x/\mu\): a polynomial in \(\lambda\) and \(\mu\) with distinct exponent pairs vanishes for all \(\lambda\) and all \(\mu>0\) only if every coefficient does.

**THEOREM R.** For every such reading and every vacuum surface \(S\) bounding \(V\):

1. **R.1.** If \(e=f=0\), Theorem Q.1 holds verbatim: \(\partial_jT^{ij}\) is absolutely integrable outside \(V\), \(\mathcal{G}(S_R)\to0\) on large spheres, and \(\mathcal{G}_i(S)=-\int_{\mathbb{R}^3\setminus V}\partial_jT^{ij}\,dV\).
2. **R.2.** If the reading is on-shell conserved, then \(e=f=0\) and \(\mathcal{G}(S)=0\) on every vacuum surface, for every interior asymmetry, every \(W\), and every \(\kappa\).
3. **R.3.** If it is not, there are an isolated device made of finitely many point charges in \(B\) and two vacuum surfaces \(S_1,S_2\) around it with \(\mathcal{G}(S_1)\neq\mathcal{G}(S_2)\).
4. **R.4 (the quadratic term does not decay).** For \(T=\delta^{ij}\Psi^2\), let \(S_1\) be the unit sphere centred at \(C\) with \(0<a=|C|<1\), and \(RS_1\) its dilate about the origin of the device. With \(Q\) the total charge and \(4\pi\varepsilon_0=1\),
\[
\lim_{R\to\infty}\mathcal{G}(RS_1)=2\pi Q^2\Big(\frac1a-\frac{1+a^2}{2a^2}\ln\frac{1+a}{1-a}\Big)\frac{C}{a},
\]
which is nonzero whenever \(Q\neq0\) (the bracket is \(-\tfrac43a+O(a^3)\) and negative on \((0,1)\)). It is \(0\) for centred spheres. Shape and placement of the audit surface fix its value at every scale.

So no reading in the family, shift-invariant or not, gives a nonzero \(\mathcal{G}\) that is the same on every audit surface around the device.

**Proof of R.1.** The proof of Q.1 used shift invariance only to get \(D+m\ge3\). With \(e=f=0\) and \(c\,\delta^{ij}\) contributing neither divergence nor flux (\(\oint n=0\)), every remaining piece has \(D+m\ge3\), and the same \(O(r^{-(D+m)})\) bound and divergence theorem apply.

**Proof of R.2.** By Lemma R.0, \(e\,\delta^{ij}\Psi\) and \(f\,\delta^{ij}\Psi^2\) are each on-shell conserved. On the harmonic \(h=x\) their divergences are \((e,0,0)\) and \((2fx,0,0)\), so \(e=f=0\). Then R.1 applies and its integrand vanishes in vacuum.

**Proof of R.3.** The proof of Q.3 never used shift invariance. Lemma Q.0 matches the full \(K\)-jet at \(y\), value included, so the same point charges, sphere, and outward bump give \(\mathcal{G}_i(S_2)\neq\mathcal{G}_i(S_1)\).

**Proof of R.4.** \(\mathcal{G}(RS_1)=\oint_{S_1}\big(R\,A_0(Ry)\big)^2\,n\,dA(y)\), and \(R\,A_0(Ry)\to Q/|y|\) uniformly on \(S_1\) because \(|y|\ge1-a>0\) there. For large \(R\), \(RS_1\) encloses \(B\), since its centre sits at distance \(Ra\) from the device and its radius is \(R\). On \(S_1\), \(y=C+n\) and \(|y|^2=1+a^2+2a\cos\theta\) with \(\theta\) the angle from \(C\); by symmetry the limit points along \(C\), with component \(2\pi Q^2\int_{-1}^{1}\frac{u\,du}{1+a^2+2au}\), which is the bracket above.

**KEPT FAILURE R.1.** Dropping shift invariance does not open a net source. No \(O(3)\)-covariant, constant-coefficient polynomial reading of any order and degree, undifferentiated factors included, combined with C5 and static massless C6, gives a surface-independent nonzero \(\mathcal{G}\) for an isolated device. The on-shell conserved ones, shift-invariant or not, give \(0\) on every vacuum surface; undifferentiated factors can appear in a conserved reading (for example \(U_1\) below), but only in pieces with \(D+m\ge3\), whose flux still vanishes. Every other reading changes when the audit surface is moved around the same device. The two pieces that escape Q.1's decay make it worse, not better: \(e\,\delta^{ij}\Psi\) grows like \(R\) under dilation (Kept Failure P.2), and \(f\,\delta^{ij}\Psi^2\) tends to a nonzero limit set by the surface's shape and placement (R.4). Interpretive step, labeled as such and the same as in Kept Failures O.1, P.2, and Q.1: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This closes the non-shift-invariant polynomial gap left by Theorem Q inside static massless C6. It is a failure of that route, not of Stage 2.

**Witness, computed by `scripts/nonshift_reading_flux.py`.** Three readings: \(U_1=\Psi\,\partial_i\partial_j\Psi-\tfrac12\delta_{ij}|\nabla\Psi|^2\) (\(m=2\), \(D=2\)), \(U_2=\delta_{ij}\Psi^2\) (\(m=2\), \(D=0\)), and \(U_3=\Psi^2\partial_i\partial_j\Psi\) (\(m=3\), \(D=2\)). Exact (rational arithmetic): \(\partial_jU_1^{ij}=\Psi\,\partial_i\Delta\Psi\), \(\partial_jU_2^{ij}=2\Psi\,\partial_i\Psi\), and \(\partial_jU_3^{ij}=\Psi\,\partial_i|\nabla\Psi|^2+\Psi^2\partial_i\Delta\Psi\) hold as polynomial identities for \(25/25\) seeded random integer \(\Psi\) of degree at most \(5\), and \(U_1[\Psi+1]-U_1[\Psi]=\partial_i\partial_j\Psi\) for all \(25\), so \(U_1\) is not shift-invariant. On the harmonic \(\Psi_h\) of Theorem Q, \(\partial_jU_1^{ij}\) vanishes identically, while at \((1,\tfrac12,-1)\) \(\partial_jU_2^{ij}=(-\tfrac{297}{32},\tfrac{621}{16},-\tfrac{81}{16})\) and \(\partial_jU_3^{ij}=(-\tfrac{3051}{4},-\tfrac{3591}{8},\tfrac{7155}{64})\). On \(h=x\), \(e\,\delta\Psi+f\,\delta\Psi^2\) is conserved for only \((e,f)=(0,0)\) among the \(49\) integer pairs in \([-3,3]^2\). Quadrature (\(128\) Gauss–Legendre by \(256\) azimuthal nodes, \(4\pi\varepsilon_0=1\)) around the cluster \(q=(3,-1,2)\) of Theorems N, O, and Q, on the same three spheres: \(U_1\) gives \(|\mathcal{G}|\le2.6\times10^{-13}\) against absolute row-\(x\) flux up to \(46.594480\); \(U_2\) gives \((-35.4343595363,\,-15.3746788139,\,13.1778307095)\), \((-69.8470404027,\,-25.997156695,\,30.5781414053)\), \((17.0626744284,\,-11.0562983204,\,-13.460682377)\), and \(U_3\) gives \((-39.1942920352,\,-13.866363922,\,12.9663781189)\), \((-77.2995955976,\,-23.6858511431,\,31.324181094)\), \((5.9282613979,\,-4.0210589129,\,-4.2759207583)\). Same device, three surfaces, three values. Each agrees with minus the exterior integral of the on-shell divergence (radial Gauss–Legendre in \(R/r\), \(96\) nodes, angular first) to \(1.2\times10^{-13}\) relative; for \(U_2\), whose divergence is only conditionally integrable, this is a numerical agreement on spheres and is not claimed by R.1. Under dilation of the unit sphere centred at \((0.3,-0.2,0.1)\), \(\mathcal{G}_{U_2}\) itself at \(R=4,\dots,64\) tends to the R.4 limit \((-42.690617678,\,28.460411786,\,-14.230205893)\) (closed form and monopole quadrature agree to \(10^{-9}\)), errors halving per doubling, Richardson value within \(2.2\times10^{-4}\) relative; \(R^3\mathcal{G}_{U_3}\) tends to \((-491.610036,\,327.740024,\,-163.870012)\), Richardson within \(1.7\times10^{-3}\). The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem R does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends. It reads \(\delta^{ij}\Psi^2\) under the massless equation only; the screened field equation (\(\Delta\Psi=m^2\Psi\), the case left at \(\tfrac{m^2}{2}\oint\Psi^2n\) in Kept Failure N.1) is taken up in Theorem T. Nor does R cover non-polynomial functions of \(\Psi\) or its derivatives, position-dependent coefficients, pseudotensor (\(\varepsilon^{ijk}\)) contractions, time-dependent or radiating fields, \(d\neq3\), or any \(\Psi_{\mathrm{info}}\neq A_0\). A nonzero \(\mathcal{G}\) on a finite surface is not thrust.

## Theorem T. Screened polynomial readings leave only surface-dependent mass flux

**ASSUMPTION T (family of readings, screened field equation; not a change to the freeze).** Theorems Q and R covered every \(O(3)\)-covariant, constant-coefficient polynomial reading of any order and degree (shift-invariant, and with undifferentiated \(\Psi\) factors allowed) under the *massless* equation. Here the reading family is the same, but the field equation is screened:
\[
(\Delta-m^2)\Psi=-\rho/\varepsilon_0,\qquad m>0\text{ fixed}.
\]
C6 still takes \(\Psi_{\mathrm{info}}=A_0\); the device is isolated (\(\rho\in C^\infty_c\) supported inside an open ball \(B\)); vacuum surfaces lie outside \(\operatorname{supp}\rho\); \(d=3\); and \(A_0\to0\) at infinity through the Yukawa kernel \(e^{-mr}/r\). On-shell conserved means \(\partial_jT^{ij}=0\) whenever \((\Delta-m^2)\Psi=0\). Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) is selected as a design target (\(m\) is only a fixed parameter of the equation under test).

**THEOREM T.** For every such reading and every vacuum surface \(S\) bounding \(V\):

1. **T.1 (Yukawa decay).** On the sphere \(S_R\) of radius \(R\) centred at the origin of the device, \(\mathcal{G}(S_R)\to0\) as \(R\to\infty\). Every monomial built from \(\Psi\) and its derivatives is \(O(e^{-cmR}R^{K})\) for some \(c\ge m\) and some \(K\), so the surface integral vanishes exponentially.
2. **T.2 (conserved readings vanish).** If the reading is on-shell conserved, then \(\mathcal{G}(S)=0\) on every vacuum surface, for every interior asymmetry, every \(W\), and every \(\kappa\). *Proof.* In vacuum \(\partial_jT^{ij}=0\), so \(\mathcal{G}(S)=\mathcal{G}(S_R)\) for every large sphere enclosing \(S\); T.1 forces the common value to \(0\).
3. **T.3 (otherwise: mass remainder, surface-dependent, decaying).** If it is not on-shell conserved, then on every vacuum surface
\[
\mathcal{G}_i(S)=-\int_{\mathbb{R}^3\setminus V}\partial_jT^{ij}\Big|_{\Delta\Psi=m^2\Psi}\,dV.
\]
The integrand is a polynomial in \(\Psi\) and its derivatives with \(\Delta\Psi\) rewritten as \(m^2\Psi\) (a *mass/screening remainder*). The value depends on which surface is drawn, and \(\mathcal{G}(S)\to0\) as the surface recedes (T.1). After removing any classical Yukawa self-force piece coming from the \(\rho\)-supported expansion \(\Delta\Psi=m^2\Psi-\rho/\varepsilon_0\) (odd under \(x\leftrightarrow y\) for the symmetric Yukawa kernel, hence \(0\) for an isolated device, as in Kept Failure N.1's screened corollary), what remains is exactly this mass flux.
4. **T.4 (massless loopholes die).** The two pieces that escaped decay under the massless equation now decay under screening: \(\oint\Psi\,n\) (Kept Failure P.2 grew like \(R\)) and \(\oint\Psi^2\,n\) (Theorem R.4 tended to a nonzero shape-dependent limit) are both \(O(e^{-cmR})\) and tend to \(0\).

So no reading in the family, under C5 and screened C6, gives a nonzero \(\mathcal{G}\) that is the same on every audit surface around the device.

**Proof of T.1.** Each factor \(\partial^k\Psi\) of a Yukawa field is bounded by \(C_{k,\varepsilon}e^{-(m-\varepsilon)|x|}\) at infinity. A monomial of \(m_\ast\) factors is \(O(e^{-c|x|}|x|^{K})\) with \(c\ge m\); times area \(R^2\) this tends to \(0\).

**Proof of T.2.** Divergence theorem on the vacuum shell between \(S\) and \(S_R\), then T.1.

**Proof of T.3.** Same shell identity with the on-shell substitution \(\Delta\Psi=m^2\Psi\) in the exterior; the resulting integrand is not identically zero (else the reading would be on-shell conserved), so two surfaces with a vacuum region between them generally disagree. Concrete remainders for the readings of Theorems M–R include \(\partial_j(\partial_i\partial_j\Psi)=m^2\partial_i\Psi\), \(\partial_jQ^{ij}=m^2\Psi\partial_i\Psi\), \(\partial_j(\delta^{ij}\Psi^2)=2\Psi\partial_i\Psi\), and \(\partial_j(\delta^{ij}\Psi)=\partial_i\Psi\); each integrates to a surface flux built from powers of \(\Psi\) and \(|\nabla\Psi|^2\).

**Proof of T.4.** Immediate from the Yukawa pointwise bounds; the witness records the decay of both fluxes on centred and off-centre spheres.

**KEPT FAILURE T.1.** Screening does not open a net source inside the polynomial family. No \(O(3)\)-covariant, constant-coefficient polynomial reading of any order and degree, undifferentiated factors included, combined with C5 and screened C6, gives a surface-independent nonzero \(\mathcal{G}\) for an isolated device. On-shell conserved readings give \(0\) on every vacuum surface (T.2). Every other reading equals a mass/screening remainder that moves when the audit surface is moved and tends to \(0\) as the surface recedes (T.3); the massless escape routes (linear growth of \(\oint\Psi\,n\), nonzero dilation limit of \(\oint\Psi^2\,n\)) are closed by Yukawa decay (T.4). Interpretive step, labeled as such and the same as in Kept Failures O.1, P.2, Q.1, and R.1: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This closes the screened gap left by Theorem R inside polynomial readings. It is a failure of that route, not of Stage 2.

**Witness, computed by `scripts/screened_reading_flux.py`.** Seven readings: Hessian \(H=\partial_i\partial_j\Psi\), quadratic \(Q\), trace \(\operatorname{Tr}=\delta_{ij}|\nabla\Psi|^2\), \(U_1=\Psi\,\partial_i\partial_j\Psi-\tfrac12\delta_{ij}|\nabla\Psi|^2\), \(U_2=\delta_{ij}\Psi^2\), \(L_1=\delta_{ij}\Psi\), and the identically conserved \(C_0=\partial_i\partial_j\Psi-\delta_{ij}\Delta\Psi\). Exact (rational arithmetic): the off-shell divergence formulae and the on-shell remainders after \(\Delta\to m^2\) hold for \(25/25\) seeded random integer \(\Psi\) of degree at most \(5\); \(C_0\) is identically divergence-free; and \(\partial_j\bigl(Q^{ij}-\tfrac{m^2}{2}\delta^{ij}\Psi^2\bigr)=\partial_i\Psi\,(\Delta\Psi-m^2\Psi)\) as a polynomial identity. Quadrature (\(128\) Gauss–Legendre by \(256\) azimuthal nodes, \(4\pi\varepsilon_0=1\), Yukawa \(m=0.7\)) around the cluster \(q=(3,-1,2)\) of Theorems N–R, on the three vacuum spheres of Theorem O: \(C_0\) and \(Q-\tfrac{m^2}{2}U_2\) give \(|\mathcal{G}|\le8.9\times10^{-14}\); \(Q\) and \(U_1\) each equal \(\tfrac{m^2}{2}\oint\Psi^2n\) to \(10^{-13}\) relative, with \(\mathcal{G}_Q=(-1.3467625605,\,-0.5331301852,\,0.4751852306)\), \((-2.4534612442,\,-0.8178402483,\,1.0318819706)\), \((0.2122600081,\,-0.1608186359,\,-0.1554593837)\); \(H=m^2\oint\Psi\,n\) matches \(L_1\) to \(10^{-13}\). Same device, three surfaces, three values for every non-conserved reading. On centred spheres \(R=2,3,5\), \(|\mathcal{G}_Q|=1.524401,\ 0.3233750,\ 0.01713515\) (Theorem N lock), \(|\mathcal{G}_{U_2}|=6.222045,\ 1.319898,\ 0.06993938\), \(|\mathcal{G}_{L_1}|=6.158927,\ 3.950475,\ 1.414125\), all strictly decreasing. Under dilation of the off-centre unit sphere of Theorem R.4, \(|\mathcal{G}_{U_2}|\) at \(R=2,4,8,16\) is \(21.18818,\ 2.275332,\ 0.04147253,\ 2.111048\times10^{-5}\) (decays to \(0\); the massless R.4 limit is killed). The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem T does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends. First-derivative non-polynomial readings of the form \(f(s)\delta+g(s)E\otimes E\) are taken up in Theorem U. Theorem T does not cover non-polynomial dependence on undifferentiated \(\Psi\) or on higher derivatives, nonlocal kernels, position-dependent coefficients, pseudotensor (\(\varepsilon^{ijk}\)) contractions, time-dependent or radiating fields, \(d\neq3\), or any \(\Psi_{\mathrm{info}}\neq A_0\). No value of \(m\) is selected. A nonzero \(\mathcal{G}\) on a finite surface is not thrust.

## Theorem U. First-derivative readings (non-polynomial allowed) are still zero or surface-dependent

**ASSUMPTION U (family of readings, not a change to the freeze).** Theorems Q–T covered every \(O(3)\)-covariant, constant-coefficient *polynomial* reading. Here the family is every reading that depends only on first derivatives of \(\Psi\) and is \(O(3)\)-covariant:
\[
T^{ij}=f(s)\,\delta^{ij}+g(s)\,\partial^i\Psi\,\partial^j\Psi,\qquad s=|\nabla\Psi|^2,
\]
with \(f,g\in C^1([0,\infty))\). By the representation theory of isotropic tensor functions of one vector, this is the general local form of that class. Affine \((f,g)\) recover Theorem O's quadratic and trace pieces (including the Maxwell stress \(Q\), \(f=-s/2\), \(g=1\)). Non-polynomial examples such as \(f(s)=e^{-s}\), \(g=0\) and \(f=0\), \(g=1/(1+s)\) are allowed. C6 (\(\Psi_{\mathrm{info}}=A_0\)), static massless fields, an isolated device inside an open ball \(B\), \(A_0\to0\) at infinity, \(d=3\), and vacuum surfaces are as in Assumption Q. On-shell conserved means \(\partial_jT^{ij}=0\) whenever \(\Delta\Psi=0\). Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters.

**THEOREM U.** For every such reading and every vacuum surface \(S\) bounding \(V\):

1. **U.1 (divergence).** Writing \(E=\nabla\Psi\),
   \[
   \partial_jT^{ij}=(f'+g/2)\,\partial_is+g'(\partial_js)\,E^iE^j+g\,E^i\Delta\Psi.
   \]
2. **U.2 (conservation criterion).** The reading is on-shell conserved if and only if \(g'\equiv0\) and \(f'+g/2\equiv0\) on \([0,\infty)\), i.e.
   \[
   T=\kappa\,\delta^{ij}+c\,Q^{ij}
   \]
   for constants \(\kappa,c\) (a constant multiple of the identity plus a multiple of the Maxwell stress).
3. **U.3 (conserved readings vanish).** If the reading is on-shell conserved, then \(\mathcal{G}(S)=0\) on every vacuum surface around an isolated device: \(\oint\kappa\,\delta\,n=0\), and the Maxwell flux vanishes by Kept Failure N.1.
4. **U.4 (otherwise: surface-dependent, decaying).** If it is not on-shell conserved, then \(\partial_jT^{ij}\) is absolutely integrable outside \(V\), \(\mathcal{G}(S_R)\to0\) on large spheres (because \(f(s)=f(0)+O(s)\) and \(\oint f(0)\,n=0\), while the \(O(s)\) remainder is \(O(R^{-4})\) on \(S_R\)), and
   \[
   \mathcal{G}_i(S)=-\int_{\mathbb{R}^3\setminus V}\partial_jT^{ij}\,dV.
   \]
   The value depends on which surface is drawn. After the Maxwell piece is removed from any reading \(T=cQ+R\), the residual flux is exactly \(\mathcal{G}[R]\); if \(R\) is not of Maxwell+\(\kappa\delta\) form that residual still moves with the surface.

So no reading in the family gives a nonzero \(\mathcal{G}\) that is the same on every audit surface around the device. Removing the Maxwell piece does not create one.

**Proof of U.1.** Differentiate \(T^{ij}=f(s)\delta^{ij}+g(s)E^iE^j\) and use \(\partial_is=2E^k\partial_iE_k\).

**Proof of U.2.** On-shell, \(\Delta\Psi=0\), so the formula reduces to \((f'+g/2)\partial_is+g'(\partial_js)E^iE^j\). If \(g'=0\) and \(f'+g/2=0\), this vanishes identically. Conversely, the harmonic polynomial \(x^2-y^2\) realizes, at points with \(s>0\) and \(\nabla s\neq0\), independent combinations of those two coefficients, so both must vanish on \((0,\infty)\); continuity extends to \(0\). Integrating \(g'=0\) and \(f'=-g/2\) gives the stated form.

**Proof of U.3.** Immediate from U.2 and Kept Failure N.1.

**Proof of U.4.** Write \(f(s)=f(0)+s\tilde f(s)\) with \(\tilde f\) continuous. The \(f(0)\delta\) flux vanishes. The remainder and the \(gE\otimes E\) term are \(O(R^{-4})\) on \(S_R\) (since \(g\) is bounded near \(0\)), so the flux is \(O(R^{-2})\to0\). The divergence is \(O(R^{-5})\) and integrable. The shell identity between \(S\) and \(S_R\), then \(R\to\infty\), gives the exterior formula. Surface dependence is the bump argument of Theorem Q.3 applied to a harmonic jet where the on-shell divergence is nonzero (U.2).

**KEPT FAILURE U.1.** Allowing non-polynomial functions of first derivatives does not open a net source. No \(O(3)\)-covariant \(C^1\) reading of the form \(f(s)\delta+g(s)E\otimes E\), combined with C5 and static massless C6, gives a surface-independent nonzero \(\mathcal{G}\) for an isolated device. The on-shell conserved ones are exactly Maxwell plus a constant multiple of \(\delta^{ij}\) and give \(\mathcal{G}=0\) on every vacuum surface (U.3). Every other reading, polynomial or not, equals minus the exterior integrated divergence, moves when the audit surface is moved, and decays like \(R^{-2}\) (U.4). After Maxwell removal the residual is either \(0\) or still surface-dependent. Interpretive step, labeled as such and the same as in Kept Failures O.1, P.2, Q.1, R.1, and T.1: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This closes the first-derivative non-polynomial gap left open by Theorems Q–T. It is a failure of that route, not of Stage 2.

**Witness, computed by `scripts/first_derivative_reading_flux.py`.** Exact (rational arithmetic): the divergence formula holds as a polynomial identity for \(40/40\) seeded random integer \(\Psi\) of degree at most \(3\) and random polynomial \((f,g)\) of degrees \(\le2\) and \(\le1\); all \(25\) affine Maxwell+\(\kappa\delta\) readings are divergence-free on the harmonic \(\Psi_h\) of Theorem Q; trace (\(f=s\), \(g=0\)), \(f=s^2\), and pure \(E\otimes E\) are not. Finite differences on the harmonic field \(x^3-3xy^2+2yz\): non-polynomial pairs \((e^{-s},0)\), \((0,1/(1+s))\), and \((-(1-e^{-s})/2,\,e^{-s})\) match the formula to relative \(10^{-6}\); the third is not conserved (\(|\mathrm{div}|\approx0.091\) at a test point) even though \(f'+g/2=0\). Quadrature (\(128\) Gauss–Legendre by \(256\) azimuthal nodes, \(4\pi\varepsilon_0=1\)) around the cluster \(q=(3,-1,2)\) of Theorems N–T, on the three vacuum spheres of Theorem O: Maxwell \(Q\) gives \(|\mathcal{G}|\le9.3\times10^{-14}\) against absolute row-\(x\) flux up to \(15.712\); trace reproduces Theorem O's \(P=\oint|E|^2n\); \(e^{-s}\delta\) gives \((5.5717749687,\,2.6920037157,\,-2.2351171899)\), \((10.5241700683,\,4.1814004818,\,-4.7864665372)\), \((-3.2630881595,\,2.1437905472,\,2.4770827176)\); \(E\otimes E/(1+s)\) gives \((-0.1238875858,\,-0.3611663762,\,0.2007774801)\), \((-0.6315002605,\,-0.6853815368,\,0.4967991642)\), \((0.9929314854,\,-0.6340327425,\,-0.8125889061)\). Same device, three surfaces, three values. Each non-conserved flux agrees with minus the exterior integral of the on-shell divergence to \(10^{-13}\) relative. On \(T=Q+e^{-s}\delta\), Maxwell removal leaves exactly the \(e^{-s}\delta\) flux. Under dilation of the off-centre unit sphere of Theorem R.4, \(R^2\mathcal{G}\) for trace, \(e^{-s}\delta\), and \(E\otimes E/(1+s)\) at \(R=4,\dots,64\) tend to a common nonzero leading coefficient (up to the factors \(-1\) and \(1/2\) respectively), so \(\mathcal{G}\to0\) like \(R^{-2}\). The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem U does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends. Dependence on undifferentiated \(\Psi\) inside the first-derivative isotropic family is taken up in Theorem V. Theorem U does not cover second or higher derivatives (beyond what Theorems Q–T already closed for polynomials), nonlocal kernels, position-dependent coefficients, pseudotensor contractions, time-dependent or radiating fields, \(d\neq3\), or any \(\Psi_{\mathrm{info}}\neq A_0\). A nonzero \(\mathcal{G}\) on a finite surface is not thrust.

## Theorem V. \(\Psi\)-dependent first-derivative isotropic readings stay Maxwell or surface-dependent

**ASSUMPTION V (family of readings, not a change to the freeze).** Theorem U covered every first-derivative \(O(3)\)-covariant reading that depends only on \(s=|\nabla\Psi|^2\). Here undifferentiated \(\Psi\) is allowed as well:
\[
T^{ij}=f(\Psi,s)\,\delta^{ij}+g(\Psi,s)\,\partial^i\Psi\,\partial^j\Psi,\qquad s=|\nabla\Psi|^2,
\]
with \(f,g\in C^1(\mathbb{R}\times[0,\infty))\). By the representation theory of isotropic tensor functions of one scalar and one vector, this is the general local form of that class. Theorem U is the \(\Psi\)-independent slice. C6 (\(\Psi_{\mathrm{info}}=A_0\)), static massless fields, an isolated device inside an open ball \(B\), \(A_0\to0\) at infinity, \(d=3\), and vacuum surfaces are as in Assumption Q/U. On-shell conserved means \(\partial_jT^{ij}=0\) whenever \(\Delta\Psi=0\). Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters.

**THEOREM V.** For every such reading and every vacuum surface \(S\) bounding \(V\):

1. **V.1 (divergence).** Writing \(E=\nabla\Psi\),
\[
\partial_jT^{ij}=(f_\Psi+s\,g_\Psi)\,E^i+(f_s+g/2)\,\partial_is+g_s\,E^i(E\cdot\nabla s)+g\,E^i\Delta\Psi.
\]
(The identity \(\partial_j(E^i)\,E^j=\tfrac12\partial_is\) is used in the \(g\) term.)

2. **V.2 (conservation criterion).** The reading is on-shell conserved for all harmonic \(\Psi\) iff
\[
g_s\equiv0,\qquad f_s+g/2\equiv0,\qquad f_\Psi+s\,g_\Psi\equiv0
\]
on \(\mathbb{R}\times[0,\infty)\). Integrating: \(g=g(\Psi)\) only, then \(f=-(g(\Psi)/2)\,s+h(\Psi)\), then \(f_\Psi+s\,g_\Psi=h'+(g'/2)\,s\equiv0\) for all \(s\) forces \(g'\equiv0\) and \(h'\equiv0\), so \(g\equiv c\), \(h\equiv\kappa\), and
\[
T=\kappa\,\delta+c\,Q
\]
(exactly Maxwell plus a constant multiple of \(\delta\) — the same conserved class as Theorem U). \(\Psi\)-dependence does not enlarge the conserved class.

3. **V.3 (conserved readings vanish).** If the reading is on-shell conserved, then \(\mathcal{G}(S)=0\) on every vacuum surface around an isolated device (as U.3 / N.1).

4. **V.4 (otherwise: surface-dependent).** If it is not on-shell conserved, then
\[
\mathcal{G}_i(S)=-\int_{\mathbb{R}^3\setminus V}\partial_jT^{ij}\,dV.
\]
The value depends on which surface is drawn. After Maxwell removal from \(T=cQ+R\), the residual is \(\mathcal{G}[R]\); if \(R\) is not of Maxwell+\(\kappa\delta\) form that residual still moves with the surface. For typical \(C^1\) growth near infinity with \(A_0\sim1/R\): pieces already treated in U (functions of \(s\) alone near \(0\)) are \(O(R^{-2})\to0\); the Taylor piece \(e^{-\Psi^2}=1-\Psi^2+O(\Psi^4)\) recovers the nonzero but surface-dependent R.4 limit of \(\delta\Psi^2\); remainders such as the pseudo-Maxwell \(g=\Psi\), \(f=-(\Psi/2)s\) still decay like \(R^{-2}\).

So no reading in the family gives a nonzero \(\mathcal{G}\) that is the same on every audit surface around the device. Removing the Maxwell piece does not create one.

**Proof of V.1.** Differentiate \(T^{ij}=f(\Psi,s)\delta^{ij}+g(\Psi,s)E^iE^j\). The chain rule gives \(f_\Psi\partial_j\Psi\,\delta^{ij}+f_s\partial_js\,\delta^{ij}\) plus the \(g\) terms. Expand \(\partial_j(gE^iE^j)=g_\Psi E^i(E\cdot E)+g_s(\partial_js)E^iE^j+g(\partial_jE^i)E^j+gE^i\Delta\Psi\). Use \(\partial_jE^i\,E^j=\tfrac12\partial_is\) and \(E\cdot E=s\) to rearrange into the displayed formula.

**Proof of V.2.** On-shell, \(\Delta\Psi=0\), so the formula reduces to the three coefficient blocks in V.2. If all three vanish, the divergence vanishes identically. Conversely, harmonic polynomials realize independent combinations of those blocks (e.g. \(x^2-y^2\) and jets with independent \(\Psi\), \(s\), \(\nabla s\)), so each coefficient must vanish. The integration argument in the statement then forces \(g\equiv c\) and \(h\equiv\kappa\).

**Proof of V.3.** Immediate from V.2 and Kept Failure N.1.

**Proof of V.4.** Shell identity between \(S\) and \(S_R\), then \(R\to\infty\) where applicable, gives the exterior formula (the same argument as Q.1 / U.4, with the R.4 caveat for even undifferentiated powers). Surface dependence is the bump argument of Theorem Q.3 applied where the on-shell divergence is nonzero (V.2).

**KEPT FAILURE V.1.** Allowing undifferentiated \(\Psi\) inside the first-derivative isotropic family does not open a surface-independent nonzero \(\mathcal{G}\) for an isolated device under C5 and static massless C6. Conserved readings remain exactly Maxwell plus \(\kappa\delta\) and give \(\mathcal{G}=0\) on every vacuum surface (V.3). Every other reading equals minus the exterior integrated divergence and moves when the audit surface is moved (V.4); after Maxwell removal the residual is either \(0\) or still surface-dependent. \(\Psi\)-dependence does not enlarge the conserved class. Interpretive step, labeled as such and the same as in Kept Failures O.1, P.2, Q.1, R.1, T.1, and U.1: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This closes the \(\Psi+\nabla\Psi\) isotropic non-polynomial first-order gap left open by Theorem U. It is a failure of that route, not of Stage 2.

**Witness, computed by `scripts/psi_gradient_reading_flux.py`.** Exact (rational arithmetic): the divergence formula holds as a polynomial identity for \(40/40\) seeded random integer \(\Psi\) of degree at most \(3\) and random polynomial \((f,g)\) in \((\Psi,s)\); all \(25\) affine Maxwell+\(\kappa\delta\) readings are divergence-free on the harmonic \(\Psi_h\) of Theorem Q; pseudo-Maxwell (\(g=\Psi\), \(f=-(\Psi/2)s\)) has on-shell \(\mathrm{div}=(s/2)E\) exactly; pure \(f=\Psi\) (\(g=0\)) has \(\mathrm{div}=E\). Finite differences on harmonic and nonharmonic fields: non-polynomial pairs \((e^{-\Psi^2},0)\), \((0,1/(1+\Psi^2))\), and \((-(s/2)e^{-\Psi^2},\,e^{-\Psi^2})\) match the formula to relative \(10^{-6}\)–\(10^{-8}\); the third satisfies \(f_s+g/2=0\) and \(g_s=0\) but not \(f_\Psi+s\,g_\Psi=0\) (\(|\mathrm{div}|\approx1.268\) at a test point). Quadrature (\(128\) Gauss–Legendre by \(256\) azimuthal nodes, \(4\pi\varepsilon_0=1\)) around the cluster \(q=(3,-1,2)\) of Theorems N–U, on the three vacuum spheres of Theorem O: Maxwell \(Q\) gives \(|\mathcal{G}|\le9.3\times10^{-14}\) against absolute row-\(x\) flux up to \(15.712\); \(e^{-\Psi^2}\delta\) gives \((0.8677036002,\,0.4738897793,\,-0.3661840665)\), \((3.8753254205,\,1.7138521206,\,-1.7879281576)\), \((-2.7289157145,\,1.6232717555,\,2.2585882375)\); \(E\otimes E/(1+\Psi^2)\) gives \((-0.3438719407,\,-0.1605731743,\,0.1357384023)\), \((-0.7373169450,\,-0.2964568398,\,0.3352190949)\), \((0.2564455823,\,-0.1551684989,\,-0.1982562974)\); pseudo-Maxwell \(\Psi\) gives \((-4.8241970412,\,-1.6431402777,\,1.5606256234)\), \((-9.6387393650,\,-2.7964529504,\,3.8240961681)\), \((0.7309160876,\,-0.5067393424,\,-0.5251162281)\). Same device, three surfaces, three values. Each non-conserved flux agrees with minus the exterior integral of the on-shell divergence to \(10^{-12}\) relative. On \(T=Q+e^{-\Psi^2}\delta\), Maxwell removal leaves exactly the \(e^{-\Psi^2}\delta\) flux. Under dilation of the off-centre unit sphere of Theorem R.4, \(e^{-\Psi^2}\delta\to-(\text{R.4 limit of }\Psi^2\delta)\) (nonzero, surface-dependent); pseudo-Maxwell and \(E\otimes E/(1+\Psi^2)\) decay like \(R^{-2}\to0\). The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem V does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends. Position-dependent coefficients in the same isotropic first-derivative family are taken up in Theorem Y. Theorem V does not cover dependence on second or higher derivatives non-polynomially, nonlocal kernels, pseudotensor contractions, time-dependent or radiating fields, \(d\neq3\), or any \(\Psi_{\mathrm{info}}\neq A_0\). A nonzero \(\mathcal{G}\) on a finite surface is not thrust.


## Theorem Y. Position-dependent first-derivative isotropic readings stay Maxwell or surface-dependent

**ASSUMPTION Y (family of readings, not a change to the freeze).** Theorem V covered every first-derivative isotropic reading \(T^{ij}=f(\Psi,s)\,\delta^{ij}+g(\Psi,s)\,\partial^i\Psi\,\partial^j\Psi\) with no explicit \(x\). Here the free coefficient functions may depend on position as well:
\[
T^{ij}=f(x,\Psi,s)\,\delta^{ij}+g(x,\Psi,s)\,\partial^i\Psi\,\partial^j\Psi,\qquad s=|\nabla\Psi|^2,
\]
with \(f,g\in C^1\) on an open set of \(\mathbb{R}^3\times\mathbb{R}\times[0,\infty)\). Theorem V is the \(x\)-independent slice. C6 (\(\Psi_{\mathrm{info}}=A_0\)), static massless fields, an isolated device inside an open ball \(B\), \(A_0\to0\) at infinity, \(d=3\), and vacuum surfaces are as in Assumption Q/U/V. On-shell conserved means \(\partial_jT^{ij}=0\) whenever \(\Delta\Psi=0\). Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters. This is about the informational reading's own coefficients, not the C5 prefactor \(W(x)\) of Theorem E.

**THEOREM Y.** For every such reading and every vacuum surface \(S\) bounding \(V\):

1. **Y.1 (divergence).** Writing \(E=\nabla\Psi\) and \(\nabla_x\) for derivatives that hit only the explicit \(x\)-dependence of \(f\) and \(g\),
\[
\partial_jT^{ij}=\partial_{x_i}f+(\nabla_x g\cdot E)\,E^i+(f_\Psi+s\,g_\Psi)\,E^i+(f_s+g/2)\,\partial_is+g_s\,E^i(E\cdot\nabla s)+g\,E^i\Delta\Psi.
\]
(The V terms are recovered when \(\nabla_x f=\nabla_x g=0\).)

2. **Y.2 (conservation criterion).** The reading is on-shell conserved for all harmonic \(\Psi\) iff the three V coefficient identities hold pointwise in \(x\),
\[
g_s\equiv0,\qquad f_s+g/2\equiv0,\qquad f_\Psi+s\,g_\Psi\equiv0,
\]
and the remaining explicit-gradient identity
\[
\nabla_x h-\tfrac12(\nabla_x g)\,|E|^2+(\nabla_x g\cdot E)\,E\equiv0\qquad\text{for all }E\in\mathbb{R}^3
\]
holds, where the V integration has already forced \(g=g(x)\), \(f=-(g(x)/2)\,s+h(x)\). The \(E=0\) case of that identity forces \(\nabla_x h=0\); then \(E=\nabla_x g\) forces \(\nabla_x g=0\). So \(g\equiv c\), \(h\equiv\kappa\) are constants, and
\[
T=\kappa\,\delta+c\,Q
\]
(exactly Maxwell plus a constant multiple of \(\delta\) — the same conserved class as Theorems U and V). Position dependence does not enlarge the conserved class.

3. **Y.3 (conserved readings vanish).** If the reading is on-shell conserved, then \(\mathcal{G}(S)=0\) on every vacuum surface around an isolated device (as V.3 / N.1).

4. **Y.4 (otherwise: surface-dependent).** If it is not on-shell conserved, then whenever the surface term at infinity vanishes,
\[
\mathcal{G}_i(S)=-\int_{\mathbb{R}^3\setminus V}\partial_jT^{ij}\,dV,
\]
and in every case the value depends on which surface is drawn. After Maxwell removal from \(T=cQ+R\), the residual is \(\mathcal{G}[R]\); if \(R\) is not of Maxwell+\(\kappa\delta\) form that residual still moves with the surface. Special cases: \(a(x)Q\) and \(a(x)\delta\) with nonconstant \(a\) are never on-shell conserved; for \(a(x)=x\) one has \(\operatorname{div}(a\delta)=e_x\) (so the exterior integral diverges) and the exact identity \(\mathcal{G}(S)=\mathrm{Vol}(V)\,e_x\) on any ball; for \(a(x)=1/(1+|x|^2)\) the reading \(a\delta\) has a nonzero but surface-dependent dilation limit (an R.4 sibling), while \(aQ\) and \(a(x)=e^{-|x|^2}\) decay.

So no reading in the family gives a nonzero \(\mathcal{G}\) that is the same on every audit surface around the device. Removing the Maxwell piece does not create one.

**Proof of Y.1.** Differentiate \(T^{ij}=f(x,\Psi,s)\delta^{ij}+g(x,\Psi,s)E^iE^j\). The chain rule on \(f\) contributes \(\partial_{x_i}f+f_\Psi E^i+f_s\partial_is\). Expanding \(\partial_j(gE^iE^j)\) adds the explicit term \((\nabla_x g\cdot E)E^i\) to the V expansion already proved in V.1.

**Proof of Y.2.** On-shell, \(\Delta\Psi=0\). At a fixed point \(x\), harmonic jets realise any value \(\Psi\), any gradient \(E\), and any traceless symmetric Hessian \(H\) independently, and \(\nabla s=2HE\). The terms carrying \(H\) are \((f_s+g/2)\,2HE+g_s\,E\,(E\cdot2HE)\); varying \(H\) over traceless matrices at fixed \(E\) forces \(g_s\equiv0\) and \(f_s+g/2\equiv0\) pointwise in \(x\) (as in V.2). Integrating in \(s\) gives \(g=g(x,\Psi)\), \(f=-(g/2)s+h(x,\Psi)\). What remains has no \(H\): \(\nabla_x h-\tfrac12 s\,\nabla_x g+(\nabla_x g\cdot E)E+(h_\Psi+\tfrac12 s\,g_\Psi)E\equiv0\). The first three terms are even in \(E\) and the last is odd, so each part vanishes separately. The odd part gives \(h_\Psi+\tfrac12 s\,g_\Psi=0\) for all \(s\ge0\), so \(h_\Psi=g_\Psi=0\) (the third V identity) and \(g=g(x)\), \(h=h(x)\). The even part is the explicit-gradient identity of Y.2; \(E=0\) gives \(\nabla h=0\), and then \(E\perp\nabla g\) with \(E\neq0\) gives \(-\tfrac12|E|^2\nabla g=0\), so \(\nabla g=0\). Hence \(g\equiv c\), \(h\equiv\kappa\).

**Proof of Y.3.** Immediate from Y.2 and Kept Failure N.1.

**Proof of Y.4.** Shell identity between \(S\) and \(S_R\), then \(R\to\infty\) where the integrand decays, gives the exterior formula (as Q.1 / V.4). When it does not decay (e.g. \(a(x)=x\) times \(\delta\)), the divergence theorem on the enclosed volume still yields an exact surface-dependent expression. Surface dependence is the bump argument of Theorem Q.3 wherever the on-shell divergence is nonzero (Y.2).

**KEPT FAILURE Y.1.** Allowing explicit position dependence inside the first-derivative isotropic family does not open a surface-independent nonzero \(\mathcal{G}\) for an isolated device under C5 and static massless C6. Conserved readings remain exactly Maxwell plus \(\kappa\delta\) with constant coefficients and give \(\mathcal{G}=0\) on every vacuum surface (Y.3). Every other reading moves when the audit surface is moved (Y.4); after Maxwell removal the residual is either \(0\) or still surface-dependent. Position dependence does not enlarge the conserved class. Interpretive step, labeled as such and the same as in Kept Failures O.1, P.2, Q.1, R.1, T.1, U.1, and V.1: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This closes the position-dependent first-derivative isotropic gap left open by Theorem V. It is a failure of that route, not of Stage 2. Theorem E (scalar \(W(x)\) reduces to the boundary) constrains the C5 prefactor, not this reading family.

**Witness, computed by `scripts/position_dependent_reading_flux.py`.** Exact (rational arithmetic): the divergence formulae \(\operatorname{div}(aQ)=(E\cdot\nabla a)E-\tfrac12 s\,\nabla a+a\,E\Delta\Psi\) and \(\operatorname{div}(a\delta)=\nabla a\) hold as polynomial identities for \(40/40\) seeded random integer \(\Psi\) of degree at most \(3\) and random polynomial \(a,b\) of degree at most \(2\); all \(25\) affine constant Maxwell+\(\kappa\delta\) readings are divergence-free on the harmonic \(\Psi_h\) of Theorem Q; \(a(x)=x\) times \(Q\) has on-shell \(\mathrm{div}=(E\cdot e_x)E-\tfrac12 s\,e_x\) exactly; \(a(x)=x\) times \(\delta\) has \(\mathrm{div}=e_x\); \(a=x\,Q+b=y\,\delta\) matches the sum. Finite differences on harmonic and nonharmonic fields: non-polynomial pairs with \(a=1/(1+|x|^2)\) and \(a=e^{-|x|^2}\) times \(Q\) or \(\delta\), and the mixed \(a(x)\Psi\,\delta\), match the Y.1 formula to relative \(10^{-6}\)–\(10^{-11}\); soft and Gaussian \(aQ\) are not on-shell conserved (\(|\mathrm{div}|\approx1.717\) and \(2.481\) at a test point). Quadrature (\(128\) Gauss–Legendre by \(256\) azimuthal nodes, \(4\pi\varepsilon_0=1\)) around the cluster \(q=(3,-1,2)\) of Theorems N–V, on the three vacuum spheres of Theorem O: Maxwell \(Q\) gives \(|\mathcal{G}|\le9.3\times10^{-14}\) against absolute row-\(x\) flux up to \(15.712\); soft \(aQ\) gives \((0,0,0)\), \((-0.6364258354,\,-0.2102975036,\,0.3359147973)\), \((0.2414780083,\,0.0154592946,\,-0.1309581046)\); soft \(a\delta\) gives \((0,0,0)\), \((-1.5237956654,\,-0.5079318885,\,0.7618978327)\), \((2.0678503489,\,0,\,-1.1488057494)\); Gaussian \(aQ\) gives \((0,0,0)\), \((-0.2716189600,\,-0.0922300380,\,0.1508885882)\), \((0.0196246838,\,0.0046628728,\,-0.0104428429)\); Gaussian \(a\delta\) gives \((0,0,0)\), \((-0.3611976657,\,-0.1203992219,\,0.1805988329)\), \((0.1381385175,\,0,\,-0.0767436208)\); \(x\,Q\) gives \((17.5779417342,\,-0.1917692645,\,0.1652142724)\), \((16.3914344267,\,-0.6150475666,\,0.7051558669)\), \((11.4162693907,\,0.0429072747,\,0.0671943269)\); \(x\,\delta\) gives \((33.5103216383,\,0,\,0)\), \((57.9058357910,\,0,\,0)\), \((113.0973355292,\,0,\,0)\). Same device, three surfaces, three values. Each decaying non-conserved flux agrees with minus the exterior integral of the on-shell divergence (absolute \(10^{-13}\) or better; symmetry zeros on the centred sphere checked absolutely). For \(x\,\delta\) the exterior integral diverges, and the exact volume identity \(\mathcal{G}=(4\pi/3)R^3 e_x\) holds on every sphere. On \(T=Q+a(x)\delta\) with soft \(a\), Maxwell removal leaves exactly the soft \(a\delta\) flux. Under dilation of an off-centre sphere: soft \(aQ\) and both Gaussian readings decay to \(0\); soft \(a\delta\) tends to a nonzero limit aligned with \(-\hat c\) whose direction depends on the offset (surface-dependent, an R.4 sibling); \(x\,Q\) decays slowly. The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem Y does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends. It does not cover dependence on second or higher derivatives non-polynomially, nonlocal kernels, pseudotensor contractions (taken up in Theorem Z), time-dependent or radiating fields, \(d\neq3\), or any \(\Psi_{\mathrm{info}}\neq A_0\). It does not reopen Theorem E: a variable C5 weight \(W(x)\) is a different object from a position-dependent reading coefficient. A nonzero \(\mathcal{G}\) on a finite surface is not thrust.



## Theorem Z. Pseudotensor (hemitropic) first-derivative readings stay Maxwell, a boundary-only curl, or surface-dependent

**ASSUMPTION Z (family of readings, not a change to the freeze).** Theorems U, V, and Y used only \(\delta^{ij}\) and \(\partial^i\Psi\partial^j\Psi\), which is the general form under the full orthogonal group \(O(3)\). If the reading is only required to be covariant under rotations \(SO(3)\), the Levi-Civita symbol is allowed, and the general second-order tensor function of one scalar and one vector gains one more term:
\[
T^{ij}=f(x,\Psi,s)\,\delta^{ij}+g(x,\Psi,s)\,\partial^i\Psi\,\partial^j\Psi+h(x,\Psi,s)\,\varepsilon^{ijk}\partial_k\Psi,\qquad s=|\nabla\Psi|^2,
\]
with \(f,g,h\in C^1\) on an open set of \(\mathbb{R}^3\times\mathbb{R}\times[0,\infty)\). Theorem Y is the \(h=0\) slice. \(\mathcal{G}_i=\oint_S T^{ij}n_j\,dA\), so the new piece contributes \(\oint h\,(n\times\nabla\Psi)_i\,dA\). C6 (\(\Psi_{\mathrm{info}}=A_0\)), static massless fields, an isolated device inside an open ball \(B\), \(A_0\to0\) at infinity, \(d=3\), and vacuum surfaces are as in Assumption Q/U/V/Y. On-shell conserved means \(\partial_jT^{ij}=0\) whenever \(\Delta\Psi=0\). Stage 1 is not edited; no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters.

**THEOREM Z.** For every such reading:

1. **Z.1 (divergence).** With \(E=\nabla\Psi\), \(H=\nabla\nabla\Psi\), and \(\nabla_x\) hitting only the explicit \(x\)-dependence,
\[
\partial_j\big(h\,\varepsilon^{ijk}E_k\big)=\big((\nabla_x h+h_\Psi E+h_s\nabla s)\times E\big)^i=(\nabla_x h\times E)^i+2h_s\,(HE\times E)^i,
\]
with no \(\Delta\Psi\) term, because \(\varepsilon^{ijk}\partial_j\partial_k\Psi=0\) identically. The full divergence is the Y.1 expression plus this.

2. **Z.2 (conservation criterion).** The reading is on-shell conserved iff \(h_s\equiv0\), \(\nabla_x h\equiv0\), and the Y.2 conditions hold. So
\[
T=\kappa\,\delta+c\,Q+h(\Psi)\,\varepsilon\cdot\nabla\Psi=\kappa\,\delta+c\,Q+\varepsilon\cdot\nabla K(\Psi),\qquad K'=h,
\]
with constants \(\kappa,c\) and any \(h\in C^1(\mathbb{R})\). The pseudotensor enlarges the conserved class, but only by a piece that is divergence-free off shell too: \(\partial_j(\varepsilon^{ijk}\partial_kK)=0\) for every \(C^2\) field.

3. **Z.3 (the new conserved piece vanishes on every closed surface).** For every closed Lipschitz surface \(S\) near which \(\Psi\) is \(C^2\), whatever lies inside (charges, the device, a surface that splits the device),
\[
\oint_S h(\Psi)\,(n\times\nabla\Psi)\,dA=\oint_S n\times\nabla K(\Psi)\,dA=0.
\]
So conserved readings give \(\mathcal{G}=0\) on every vacuum surface around an isolated device, as in Y.3.

4. **Z.4 (otherwise: surface-dependent).** If the reading is not on-shell conserved, \(\mathcal{G}\) moves with the audit surface (bump argument of Theorem Q.3), and \(\mathcal{G}_i(S)=-\int_{\mathbb{R}^3\setminus V}\partial_jT^{ij}\,dV\) whenever the surface term at infinity vanishes. For \(h=h(x)\) alone, Stokes on \(K=h\Psi\) gives the exact reduction \(\oint h\,(n\times\nabla\Psi)=-\oint\Psi\,(n\times\nabla h)\); for \(h=x_m\) this is \(e_m\times\oint\Psi n\), the flux of the undifferentiated reading of Kept Failure P.2, and on a sphere of centre \(c\) enclosing the device it equals \(\tfrac{4\pi}{3}\,e_m\times p_c\) with \(p_c=\sum q\,(c_q-c)\) (units \(4\pi\varepsilon_0=1\)) at every radius. That depends on the centre whenever the net charge is nonzero.

So no reading in the family gives a nonzero \(\mathcal{G}\) that is the same on every audit surface around the device. Removing the Maxwell piece leaves \(0\), a boundary-only curl that is \(0\), or a surface-dependent remainder.

**Proof of Z.1.** \(\partial_j(h\varepsilon^{ijk}E_k)=\varepsilon^{ijk}(\partial_jh)E_k+h\,\varepsilon^{ijk}\partial_j\partial_k\Psi\). The second term is zero by antisymmetry. The chain rule gives \(\partial_jh=\partial_{x_j}h+h_\Psi E_j+h_s\partial_js\), \(\nabla s=2HE\), and \(E\times E=0\) removes \(h_\Psi\).

**Proof of Z.2.** On shell, at a fixed point \(x\), harmonic jets \(\Psi_0+E\cdot(y-x)+\tfrac12(y-x)\cdot H(y-x)\) realise any \(\Psi\), any \(E\), and any traceless symmetric \(H\) independently. For \(E\neq0\), \(H\mapsto HE\) maps traceless symmetric matrices onto \(\mathbb{R}^3\): \(H=Ew^{\top}+wE^{\top}\) with \(w\perp E\) is traceless with \(HE=s\,w\), and \(H=EE^{\top}-\tfrac s3 I\) gives \(HE=\tfrac{2s}{3}E\). Write \(v=HE\). The terms carrying \(H\) are \(2(f_s+g/2)\,v+2g_s(E\cdot v)E+2h_s\,v\times E\), and the rest does not depend on \(H\), so this linear map in \(v\) must vanish. Taking \(v\perp E\), \(v\neq0\), the vectors \(v\) and \(v\times E\) are orthogonal and nonzero, so \(f_s+g/2=0\) and \(h_s=0\). Taking \(v=E\) gives \(g_s=0\). These hold for \(s>0\) and extend to \(s=0\) by continuity. So \(g=g(x,\Psi)\), \(f=-(g/2)s+k(x,\Psi)\), \(h=h(x,\Psi)\), and what is left is
\[
\nabla_xk-\tfrac12s\,\nabla_xg+(\nabla_xg\cdot E)E\;+\;(k_\Psi+\tfrac12s\,g_\Psi)E+\nabla_xh\times E\equiv0.
\]
The first three terms are even in \(E\) and the last two odd, so each group vanishes. The even group is the explicit-gradient identity of Y.2, which forces \(\nabla_xk=\nabla_xg=0\). In the odd group, the component along \(E\) gives \(k_\Psi+\tfrac12s\,g_\Psi=0\) for all \(s\), so \(k_\Psi=g_\Psi=0\), and the component across \(E\) gives \(\nabla_xh\times E=0\) for every direction of \(E\), so \(\nabla_xh=0\). Hence \(k=\kappa\), \(g=c\), \(h=h(\Psi)\). Conversely each of the three pieces is conserved on shell (Y.2 and Z.1).

**Proof of Z.3.** Put \(F=K(\Psi)\,e_i\), which is \(C^2\) on a neighbourhood \(U\) of \(S\). Then \((n\times\nabla K)_i=\varepsilon^{ijk}n_j\partial_kK=n\cdot(\nabla\times F)\), since \((\nabla\times F)_j=\varepsilon^{jli}\partial_lK=\varepsilon^{ijl}\partial_lK\). Take a cutoff \(\chi\in C_c^\infty(U)\) with \(\chi=1\) near \(S\). Then \(\nabla\times(\chi F)\) is \(C^1\) on all of \(\mathbb{R}^3\), equals \(\nabla\times F\) on \(S\), and has zero divergence, so the divergence theorem on \(V\) gives \(\oint_S n\cdot\nabla\times F=\int_V\nabla\cdot\nabla\times(\chi F)=0\). Nothing about \(\Psi\) inside \(V\) away from \(U\) is used. The Y.3 part is Kept Failure N.1.

**Proof of Z.4.** Surface dependence is the bump argument of Theorem Q.3 wherever the on-shell divergence is nonzero (Z.2), and the exterior formula is the shell identity of Q.1 / V.4 when the integrand decays. For \(h=h(x)\): \(\nabla(h\Psi)=\Psi\nabla h+h\nabla\Psi\), and Z.3 applied to \(h\Psi\) gives the reduction. For a sphere of radius \(R\) about \(c\) containing a point charge at \(a\), only the \(\ell=1\) term survives in \(\oint n/|y-a|\,dA=\tfrac{4\pi}{3}(a-c)\); sum over charges.

**KEPT FAILURE Z.1.** Allowing pseudotensor (Levi-Civita) contractions in the first-derivative family does not open a surface-independent nonzero \(\mathcal{G}\) for an isolated device under C5 and static massless C6. The conserved class does grow, by \(\varepsilon\cdot\nabla K(\Psi)\), but that piece is a curl and its flux is \(0\) on every closed surface (Z.3). Every other reading moves when the audit surface is moved (Z.4). The position-dependent \(x_m\,\varepsilon\cdot\nabla\Psi\) is the P.2 flux \(\oint\Psi n\) rotated by \(e_m\times\), so it adds nothing beyond P.2. Interpretive step, labeled as such and the same as in Kept Failures O.1, P.2, Q.1, R.1, T.1, U.1, V.1, and Y.1: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This closes the pseudotensor first-derivative gap left open by Theorem Y. It is a failure of that route, not of Stage 2.

**Witness, computed by `scripts/hemitropic_reading_flux.py`.** Exact (sympy): the Z.1 formula holds for \(30/30\) seeded random integer \(\Psi\) of degree at most \(3\) and random polynomial \(h(x,\Psi,s)\); \(h(\Psi)\in\{\Psi^2,\sin\Psi,e^{-\Psi},1/(1+\Psi^2)\}\) times \(\varepsilon\cdot\nabla\Psi\) is divergence-free on a non-harmonic \(\Psi\); on a harmonic cubic at \((1,\tfrac12,-1)\), \(h=s\) gives on-shell divergence \(2HE\times E=(-265,\tfrac{517}{16},-\tfrac{1021}{4})\), \(h=x\) gives \(e_x\times E=(0,\tfrac12,-6)\), and \(h=\Psi s\) gives \(2\Psi\,HE\times E\neq0\); \(H\mapsto HE\) has rank \(3\) for \(20/20\) random integer \(E\neq0\), and the two explicit matrices of the proof check exactly. Quadrature (\(128\) Gauss–Legendre by \(256\) azimuthal nodes, \(4\pi\varepsilon_0=1\)) around the cluster \(q=(3,-1,2)\) of Theorems N–Y on its three vacuum spheres: Maxwell \(Q\), \(\Psi\,\varepsilon\cdot E\), \(e^{-\Psi}\varepsilon\cdot E\), and \(Q+\Psi\,\varepsilon\cdot E\) give \(|\mathcal{G}|\le1.5\times10^{-13}\) against absolute row-\(x\) flux up to \(22.18\). On a sphere of radius \(0.3\) around the \(q=3\) charge alone (not a vacuum surface for the device), \(\Psi\,\varepsilon\cdot E\) still gives \(6.7\times10^{-14}\) while Maxwell gives the force on that charge, \((-36.084231,\,56.288828,\,-32.298463)\). The non-conserved \(s\,\varepsilon\cdot E\) gives \((0.0396736069,\,0.1242382122,\,0.2473772397)\), \((0.2264880458,\,0.5105596407,\,0.9389671858)\), \((0.0058596383,\,0.0104726951,\,-0.0085588517)\) on the three spheres, each equal to minus the exterior integral of \(2HE\times E\) to \(4\times10^{-14}\), and decays under dilation. \(x\,\varepsilon\cdot E\) gives \((0,\,-3.3510321638,\,-3.9793506945)\), \((0,\,-8.3775804096,\,-7.3303828584)\), \((0,\,5.0265482457,\,-3.9793506945)\), matching \(e_x\times\oint\Psi n\) and \(\tfrac{4\pi}{3}e_x\times p_c\) to \(5\times10^{-12}\) on these spheres and at \(R=50\) and \(R=200\), so it does not decay and its value is set by the centre. On \(Q+s\,\varepsilon\cdot E\), Maxwell removal leaves exactly the \(s\,\varepsilon\cdot E\) flux. Same device, three surfaces, three values. The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem Z does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends, or whether a parity-odd reading is physically admissible. It does not cover second or higher derivatives non-polynomially, nonlocal kernels, time-dependent or radiating fields, \(d\neq3\), or any \(\Psi_{\mathrm{info}}\neq A_0\). A nonzero \(\mathcal{G}\) on a finite surface is not thrust.


## Theorem AV. Nonlocal (convolution and fractional) readings: zero only for the local kernel, otherwise surface-dependent, never a surface-independent nonzero value

**ASSUMPTION AV (family of readings, not a change to the freeze).** C6 (\(\Psi_{\mathrm{info}}=A_0\)), static massless fields, \(d=3\), \(\Delta\Psi=-\rho/\varepsilon_0\) with \(\rho\) a finite charge measure of compact support (point charges allowed), and \(\Psi\) the Newtonian potential, so \(\Psi\to0\) at infinity. A nonlocal reading replaces \(\Psi\) by \(\Phi_a=K_a\Psi\), where \(K_a\) is translation-invariant and is one of two kinds. (i) Convolution with a finite Borel measure \(K\), radial when rotation covariance is needed. Examples are the heat kernel \(e^{t\Delta}\), the Yukawa smearing \(\mu^2(\mu^2-\Delta)^{-1}\), and ball or shell averages. (ii) The fractional multiplier \((-\Delta)^{s/2}\) with \(0<s<2\). The stress is the symmetric Maxwell form
\[
T[\Phi_1,\Phi_2]^{ij}=\tfrac12\big(\partial^i\Phi_1\partial^j\Phi_2+\partial^i\Phi_2\partial^j\Phi_1\big)-\tfrac12\delta^{ij}\,\nabla\Phi_1\cdot\nabla\Phi_2 .
\]
\(T[\Phi,\Phi]\) is the Maxwell form of \(\Phi\), and \(T[\Psi,\Psi]=Q\) is Theorem N's reading. The effective sources are \(\sigma_a:=-\Delta\Phi_a\): \(K*\rho/\varepsilon_0\) in case (i) and \((-\Delta)^{s/2}\rho/\varepsilon_0\) in case (ii). Unless \(K\) has compact support, \(\sigma_a\) is not compactly supported. In both cases \(\Phi_a=\mathcal N*\sigma_a\) with \(\mathcal N=1/(4\pi|x|)\): use Tonelli in case (i) and the Fourier multiplier \(|k|^{s}\) in case (ii). Write \(U=\mathbb{R}^3\setminus\operatorname{supp}\rho\). A vacuum surface is a closed Lipschitz surface \(S\subset U\) that encloses \(\operatorname{supp}\rho\), and \(\mathcal{G}(S)=\oint_ST^{ij}n_j\,dA\). The regularity hypothesis is **(H)**: \(\Phi_a\in C^2(U)\), \(|\nabla\Phi_a(x)|=o(|x|^{-1})\), and \(\sigma_a\nabla\Phi_b\in L^1(\{|x|>R_0\})\). For a nonlocal reading, on-shell conserved means \(\partial_jT^{ij}=0\) at every point of \(U\) for every admissible source. A pointwise condition like \(\Delta\Psi(x)=0\) no longer fixes \(\sigma(x)\). Stage 1 is not edited, and no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), or \(m\) enters.

**THEOREM AV.**

1. **AV.1 (divergence; exterior and interior formulas).** Off shell,
\[
\partial_jT[\Phi_1,\Phi_2]^{ij}=\tfrac12\big(\partial^i\Phi_1\,\Delta\Phi_2+\partial^i\Phi_2\,\Delta\Phi_1\big)=-\tfrac12\big(\sigma_2\,\partial^i\Phi_1+\sigma_1\,\partial^i\Phi_2\big).
\]
Under (H), every vacuum surface satisfies
\[
\mathcal{G}(S)=-\int_{\mathrm{ext}\,S}\partial_jT^{ij}\,dV=\tfrac12\int_{\mathrm{ext}\,S}\big(\sigma_2\nabla\Phi_1+\sigma_1\nabla\Phi_2\big)\,dV .
\]
Suppose instead that \(\sigma_a\in L^1\cap L^\infty\cap C^1(\mathbb{R}^3)\). This holds for any finite measure \(K\) acting on \(\rho\in C_c^\infty\), and for the heat kernel acting on point charges. Then \(\int_{\mathbb{R}^3}\partial_jT^{ij}\,dV=0\) (Newton's third law for the effective sources, i.e. the \(k=0\) Fourier mode of \(\operatorname{div}T\) vanishes). So \(\mathcal{G}(S)=\int_V\partial_jT^{ij}\,dV=-\int_{\mathrm{ext}\,S}\partial_jT^{ij}\,dV\), whether or not the decay part of (H) holds. For \(T[\Phi,\Phi]\), \(\varepsilon_0\mathcal{G}(S)\) is the electrostatic force that the effective charge \(\varepsilon_0\sigma\) outside \(S\) exerts on the effective charge inside \(S\). It is an internal force of the smeared distribution, split by \(S\).

2. **AV.2 (no surface-independent nonzero value).** Let \(S_k\) be vacuum surfaces whose exteriors lie in \(\{|x|>R_k\}\) with \(R_k\to\infty\). Under (H), \(\mathcal{G}(S_k)\to0\). So if \(\mathcal{G}\) has the same value on every vacuum surface around the device, that value is \(0\). The proof uses only \(T\in C^1(U)\), \(\oint_{|x|=R}|T|\,dA\to0\), and \(\operatorname{div}T\in L^1\) near infinity. So AV.2 holds for every reading with that decay, local or not, linear or not, radial or not.

3. **AV.3 (conservation criterion for convolution readings).** Let \(K\) be a finite radial Borel measure, \(\Phi=K*\Psi\), and take the reading \(T[\Phi,\Phi]\) or the cross reading \(T[\Psi,\Phi]\). Over sources \(\rho\in C_c^\infty\) the following are equivalent: (a) \(K=c\,\delta\); (b) \(\operatorname{div}T=0\) on \(U\) for every \(\rho\); (c) \(\mathcal{G}(S)=0\) on every vacuum surface for every \(\rho\). Under (a) the reading is \(c^2Q\) or \(cQ\). For \((-\Delta)^{s/2}\) with \(0<s<2\), (b) and (c) fail: a single point charge already gives \(\operatorname{div}T\neq0\) at every vacuum point (AV.4).

4. **AV.4 (one radial source, exact).** Let \(\Phi\) be radial about \(a\) and \(C^2\) for \(r=|x-a|>r_s\), with \(\operatorname{supp}\rho\subset\bar B_{r_s}(a)\), (H) holding, and \(m(r):=-r^2\Phi'(r)\), so \(\sigma=m'/r^2\). Take a sphere \(S\) of radius \(R\) and centre \(c\), with \(b=|c-a|\) and \(r_s+b<R\). Then the reading \(T[\Phi,\Phi]\) gives
\[
\mathcal{G}(S)=\pi\,\hat e_b\int_{R-b}^{R+b}\frac{m(r)\,m'(r)\,(1-u_0(r)^2)}{r^2}\,dr,\qquad u_0=\frac{r^2+b^2-R^2}{2rb},\quad \hat e_b=\frac{c-a}{b},
\]
and \(\mathcal{G}=0\) for a centred sphere (\(b=0\)). For a point charge \(q\) (units \(4\pi\varepsilon_0=1\), \(\Psi=q/r\)) the cases are as follows. Maxwell has \(m\equiv q\), so \(\mathcal{G}=0\) on every sphere. The heat kernel gives \(m=q\big(\operatorname{erf}\tfrac{r}{2\sqrt t}-\tfrac{r}{\sqrt{\pi t}}e^{-r^2/4t}\big)\), and Yukawa gives \(m=q\big(1-(1+\mu r)e^{-\mu r}\big)\). In both, \(m'>0\), so \(0<\mathcal{G}\cdot\hat e_b\le\pi q^2\big(1-m(R-b)/q\big)/(R-b)^2\) on every off-centre sphere. The ceiling decays like a Gaussian or exponential tail. The fractional reading gives \(\Phi=qC_sr^{-1-s}\) with \(C_s=\tfrac{2}{\pi}\Gamma(s)\sin\tfrac{\pi s}{2}\), so \(\mathcal{G}\cdot\hat e_b=-\pi s(1+s)^2C_s^2q^2\int_{R-b}^{R+b}r^{-3-2s}(1-u_0^2)\,dr=-\pi s(1+s)^2C_s^2q^2\,\tfrac{4b}{3}R^{-3-2s}\big(1+O(b^2/R^2)\big)\). In closed form,
\[
s=1:\ \ \mathcal{G}=-\frac{64}{3\pi}\,\frac{q^2R^3b}{(R^2-b^2)^4}\,\hat e_b,\qquad s=\tfrac12:\ \ \mathcal{G}=-\frac35\,\frac{q^2\,b\,(5R^2-b^2)}{(R^2-b^2)^3}\,\hat e_b .
\]

So in this family, removing the Maxwell piece leaves \(0\) only for the local kernel. Otherwise it leaves a remainder that moves with the audit surface and tends to \(0\) at infinity. No nonlocal reading in the family gives a nonzero \(\mathcal{G}\) that is the same on every audit surface around an isolated static device.

**Proof of AV.1.** We have \(\partial_j(\partial^i\Phi_1\partial^j\Phi_2)=\partial^i\Phi_1\Delta\Phi_2+\partial_j\partial^i\Phi_1\,\partial^j\Phi_2\), and the same holds with \(1\leftrightarrow2\). The \(\delta\) term contributes \(-\tfrac12\partial^i(\nabla\Phi_1\cdot\nabla\Phi_2)\), which cancels the two Hessian terms. Under (H), \(T\in C^1\) between \(S\) and the sphere \(S_R\), and the divergence theorem gives \(\mathcal{G}(S)=\oint_{S_R}T n-\int_{\mathrm{between}}\operatorname{div}T\). Here \(\oint_{S_R}|T|\le4\pi R^2\sup_{S_R}|\nabla\Phi_1||\nabla\Phi_2|=o(1)\), and \(\operatorname{div}T\in L^1\) near infinity, so letting \(R\to\infty\) gives the exterior formula. For the interior form, write \(\int\sigma_2\nabla\Phi_1=\iint\sigma_2(x)\sigma_1(y)\nabla\mathcal N(x-y)\,dx\,dy\). This converges absolutely because \(\int|\sigma_1(y)||x-y|^{-2}dy\le4\pi\|\sigma_1\|_\infty+\|\sigma_1\|_1\). Swap \(x\leftrightarrow y\) and use \(\nabla\mathcal N(-z)=-\nabla\mathcal N(z)\) to get \(\int\sigma_1\nabla\Phi_2=-\int\sigma_2\nabla\Phi_1\), so \(\int_{\mathbb{R}^3}\operatorname{div}T=0\). In Fourier terms this is \(\widehat{\operatorname{div}T}(0)=0\): the integrand \(k\,|\hat\rho|^2\,\mathrm{Re}(\hat K_1\overline{\hat K_2})/k^2\) is odd. A \(C^1\) density with \(L^1\cap L^\infty\) bounds has a \(C^2\) Newtonian potential, so the divergence theorem on \(V\) gives \(\mathcal{G}(S)=\int_V\operatorname{div}T\). For the force reading, the force density on \(\varepsilon_0\sigma\) is \(-\varepsilon_0\sigma\nabla\Phi=\varepsilon_0\operatorname{div}T\). The self-force of the inside part is \(0\) by the same antisymmetry.

**Proof of AV.2.** By AV.1, \(\mathcal{G}(S_k)=-\int_{\mathrm{ext}\,S_k}\operatorname{div}T\). The integrand is dominated by \(|\operatorname{div}T|\in L^1(\{|x|>R_0\})\), and the domains shrink to the empty set. Only the shell identity and that decay were used.

**Proof of AV.3.** (a)\(\Rightarrow\)(b): \(\Phi=c\Psi\) and \(\Delta\Psi=0\) on \(U\). (a)\(\Rightarrow\)(c) is Kept Failure N.1. For (b)\(\Rightarrow\)(a) and (c)\(\Rightarrow\)(a), suppose \(K\neq c\delta\). Then some \(f\in C_c(\mathbb{R}^3\setminus\{0\})\), supported in \(\{|x|>2\varepsilon_1\}\), has \(\int f\,dK\neq0\). Let \(\rho_\varepsilon=\varepsilon_0\varphi_\varepsilon\) with \(\varphi_\varepsilon\ge0\) smooth, radial, supported in \(B_\varepsilon\), and of unit mass. The effective source \(\sigma_\varepsilon=K*\varphi_\varepsilon\) is smooth, radial, and in \(L^1\cap L^\infty\). Also \(\int f\sigma_\varepsilon=\int(f*\check\varphi_\varepsilon)\,dK\to\int f\,dK\), because \(f*\check\varphi_\varepsilon\to f\) uniformly and \(|K|\) is finite. So for some \(\varepsilon<\varepsilon_1\), \(\sigma_\varepsilon(x_1)\neq0\) at a point with \(|x_1|>2\varepsilon_1>\varepsilon\). Newton's shell theorem gives \(\nabla\Phi=-M(r)\,x/(4\pi r^3)\) with \(M(r)=\int_{B_r}\sigma_\varepsilon\), and \(\nabla\Psi=-x/(4\pi r^3)\) for \(r>\varepsilon\). For the cross reading, \(\operatorname{div}T=-\tfrac12\sigma\nabla\Psi\neq0\) at \(x_1\). For \(T[\Phi,\Phi]\), \(\operatorname{div}T=\sigma M\,x/(4\pi r^3)\). On an interval \(I\ni|x_1|\) where \(\sigma\neq0\), \(M'=4\pi r^2\sigma\neq0\), so \(M\) vanishes at most once on \(I\), and \(\operatorname{div}T(x_0)\neq0\) at some \(x_0\) with \(|x_0|\in I\). That contradicts (b). For (c), let \(S_1\) be the sphere \(|x|=r_1\) with \(\varepsilon<r_1<|x_0|\), placed so that \(|x_0|=r_1+\delta/2\). Then \(\mathcal{G}(S_1)=\tfrac12\oint\Phi'^2n=0\) by symmetry. Let \(S_2=\partial(B_{r_1}\cup B_\delta(x_0))\), a vacuum surface. Then \(\mathcal{G}(S_2)-\mathcal{G}(S_1)=\int_{B_\delta(x_0)\setminus B_{r_1}}\operatorname{div}T\). This is nonzero for small \(\delta\), because one component of \(\operatorname{div}T\) keeps a sign near \(x_0\). That contradicts (c). The fractional statement is AV.4: \(\sigma=-s(1+s)C_sq\,r^{-3-s}\) and \(m=(1+s)C_sq\,r^{-s}\) are both nonzero for every \(r>0\).

**Proof of AV.4.** Here \(\nabla\Phi=-(m/r^2)\hat y\) and \(\Delta\Phi=-m'/r^2\), with \(y=x-a\), so \(\operatorname{div}T=\nabla\Phi\,\Delta\Phi=(mm'/r^4)\,\hat y\). In the exterior formula, shells with \(r>R+b\) are full and contribute \(0\). Shells with \(r<R-b\) lie inside \(S\). For \(R-b<r<R+b\), the exterior part of the shell is the cap \(\hat y\cdot\hat e_b<u_0\), and \(\int_{\mathrm{cap}}\hat y\,dA=\hat e_b\,2\pi r^2\int_{-1}^{u_0}u\,du=-\hat e_b\,\pi r^2(1-u_0^2)\). The kernel formulas follow from \(e^{t\Delta}(q/r)=q\operatorname{erf}(r/2\sqrt t)/r\) (heat equation), \((\mu^2-\Delta)\Phi=\mu^2q/r\) for \(\Phi=q(1-e^{-\mu r})/r\), and the radial inverse transform \(\Phi(r)=\tfrac{2q}{\pi r}\int_0^\infty k^{s-1}\sin(kr)\,dk=\tfrac{2q}{\pi}\Gamma(s)\sin\tfrac{\pi s}{2}\,r^{-1-s}\) of \(\hat\Phi=4\pi q|k|^{s-2}\). The last integral is Dirichlet's for \(s<1\) and is Abel-summed for \(1\le s<2\). It matches the Riesz constant \(2^s\Gamma(\tfrac{1+s}2)/(\sqrt\pi\,\Gamma(1-\tfrac s2))\). The ceiling uses \(0<m<q\) and \(\int_{R-b}^{R+b}m'/r^2\le(q-m(R-b))/(R-b)^2\). The asymptotic form substitutes \(r=R+bv\), giving \(u_0=v+O(b/R)\) and \(\int_{-1}^1(1-v^2)\,dv=4/3\). The two closed forms are exact integrals of rational functions.

**KEPT FAILURE AV.1.** Replacing \(\Psi\) by a translation-invariant nonlocal reading of it does not open a surface-independent nonzero \(\mathcal{G}\) for an isolated static device under C5 and static massless C6. The readings tested are the heat kernel, Yukawa smearing, a general finite radial measure, and the fractional Laplacian. The effective source \(\sigma=K*\rho/\varepsilon_0\) leaks into the vacuum, so conservation on the vacuum region fails for every kernel except \(c\delta\) (AV.3). \(\mathcal{G}\) is then nonzero on some vacuum surfaces, even for a single point charge, but it moves with the surface. A centred sphere gives \(0\), and an off-centre sphere gives a value whose direction is set by the offset. Under dilation \(\mathcal{G}\) tends to \(0\): like the kernel's tail for heat or Yukawa, and like \(R^{-3-2s}\) for \((-\Delta)^{s/2}\) (AV.2, AV.4). Newton's third law for the effective charge makes the total \(\int_{\mathbb{R}^3}\operatorname{div}T\) exactly \(0\), so the remainder is an internal force of the smeared charge, booked across \(S\). Interpretive step, labeled as such and the same as in Kept Failures O.1, P.2, Q.1, R.1, T.1, U.1, V.1, Y.1, and Z.1: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent \(\mathcal{G}\) cannot be read as one. This closes the nonlocal (convolution and fractional Maxwell-form) gap left open by Theorems V, Y, and Z. It is a failure of that route, not of Stage 2.

**Witness, computed by `scripts/nonlocal_reading_flux.py`.**
- Exact (sympy): the AV.1 divergence identity holds for \(30/30\) seeded random integer cubic pairs \((\Phi_1,\Phi_2)\).
- Exact: the heat, Yukawa, and fractional kernel formulas, \(m\), and \(\sigma=m'/r^2\) are verified, as is \(m'>0\) for both positive kernels.
- \(C_s\) equals the Riesz constant at \(s=\tfrac14,\tfrac12,1,\tfrac32,\tfrac74\) to \(10^{-25}\), and the \(s=\tfrac12\) Dirichlet integral equals \(\sqrt{\pi/2}\) to \(10^{-20}\).
- Exact closed forms: \(\int_{R-b}^{R+b}r^{-5}(1-u_0^2)\,dr=4R^3b/(3(R^2-b^2)^4)\). For a unit charge at \((R,b)=(2,\tfrac12),(2,1),(3,1)\), \(s=1\) gives \(\mathcal{G}\cdot\hat e_b=-\tfrac{65536}{151875\pi},-\tfrac{512}{243\pi},-\tfrac{9}{64\pi}\), and \(s=\tfrac12\) gives \(-\tfrac{632}{5625},-\tfrac{19}{45},-\tfrac{33}{640}\). The \(s=1\) series is \(-(64b/3\pi)R^{-5}+O(R^{-7})\).
- One unit charge, 1D formula against a \(128\times256\) surface quadrature, agreeing to \(1.2\times10^{-14}\), on spheres \((R,b)=(2,0),(2,\tfrac12),(2,1),(3,1)\):
  - Maxwell: \(0\) on all four.
  - Heat kernel, \(t=\tfrac14\): \(0,\ 0.1163628425531,\ 0.4313643009608,\ 0.01228902003482\).
  - Yukawa, \(\mu=1\): \(0,\ 0.08832298677601,\ 0.2021386130048,\ 0.06782121850442\).
  - \((-\Delta)^{1/2}\) and \((-\Delta)^{1/4}\): negative off centre, equal to the exact values above.
- The positive-kernel ceiling holds at \(R=2,3,6\). The \((-\Delta)^{1/4}\) leading term \(-\tfrac94\cdot\tfrac{4b}{3}R^{-4}\) is matched to ratios \(1.0286,1.0070,1.0018\) at \(R=10,20,40\).
- Coulomb cluster \(q=(3,-1,2)\) of Theorems N–Z, on its three vacuum spheres:
  - Heat kernel: \((1.9216953999,0.6370531893,-0.6085148589)\), \((3.539183942,0.8309375122,-1.314652441)\), \((0.0831547402,0.1220272808,-0.0600820331)\).
  - Yukawa: \((1.4344502109,0.5915081429,-0.5173454202)\), \((2.5199205479,0.8488218105,-1.0655074819)\), \((-0.3265550523,0.2917522815,0.2468140837)\).
  - \((-\Delta)^{1/2}\): \((-2.2374241597,-0.7720884765,0.7284680158)\), \((-4.7064016253,-1.2940127131,1.8225234138)\), \((0.1814520415,-0.1797628144,-0.1251353099)\).
  - Each value equals \(\int_{\mathrm{ext}}\sigma\nabla\Phi\) to \(10^{-10}\). The Maxwell piece on the same spheres is at most \(1.0\times10^{-13}\), so removing it leaves exactly these values.
  - The heat-kernel values also equal the interior integral \(-\int_V\sigma\nabla\Phi\) to \(1.5\times10^{-13}\), which checks Newton's third law \(\int_{\mathbb{R}^3}\operatorname{div}T=0\).
  - The cross reading \(T[\Psi,e^{\Delta/4}\Psi]\) gives \((1.0653944,0.33417684,-0.32771438)\), \((2.14639883,0.49534143,-0.78857187)\), \((0.03635805,0.05977247,-0.02752305)\), equal to \(\tfrac12\int_{\mathrm{ext}}\sigma\nabla\Psi\).
- Dilation about \((0.6,0.2,-0.3)\):
  - Heat: \(|\mathcal{G}|=4.10\times10^{-1},2.19\times10^{-3},1.82\times10^{-6}\) at \(R=3,4,5\).
  - Yukawa: \(1.25\times10^{-1},4.19\times10^{-4},9.22\times10^{-9}\) at \(R=5,10,20\).
  - \((-\Delta)^{1/2}\): \(4.96\times10^{-2},1.30\times10^{-3},3.88\times10^{-5}\) at \(R=5,10,20\).

Same device, three surfaces, three values, all tending to \(0\). The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem AV does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends. AV.3's criterion is proved for the Maxwell-form (self and cross) readings of finite radial measures. Other nonlocal tensors, such as nonlocal Hessian or higher-derivative nonlinear functionals, are covered only by the general no-go AV.2, which needs the \(o(R^{-2})\) surface decay. Nonlocal undifferentiated pieces such as \(\delta^{ij}K*\Psi\) or \(\delta^{ij}(K*\Psi)^2\) are the P.2/R.1 situation and are not reopened here. Smoothing kernels that are not finite measures, such as \((-\Delta)^{-\tau}\) with \(\tau\ge\tfrac12\) acting on a charged device, violate the decay in (H) and are not covered. Kernels with non-\(C^2\) output on the vacuum region (shell averages acting on point charges) are covered by AV.3 only through smooth sources. Time-dependent or radiating fields, \(d\neq3\), and any \(\Psi_{\mathrm{info}}\neq A_0\) remain OPEN. A nonzero \(\mathcal{G}\) on a finite surface is not thrust.


## Theorem AW. Slow-tail smoothing readings \((-\Delta)^{-\tau}\Psi\): a surface-independent value is always \(0\), and the marginal readings give a radius-free but centre-dependent remainder

**ASSUMPTION AW (family of readings, not a change to the freeze).** These are the C6 and static massless hypotheses of Assumption AV, with \(d=3\), units \(4\pi\varepsilon_0=1\), and a compactly supported charge distribution \(\tilde\rho:=-\Delta\Psi\). Point charges and point dipoles are allowed, so a unit point charge has \(\Psi=1/r\). The reading is the slow-tail smoothing \(\Phi=(-\Delta)^{-\tau}\Psi\), \(\tau>0\), in the Maxwell form \(T[\Phi,\Phi]\) of Assumption AV, or the cross form \(T[\Psi,\Phi]\). We have \(\operatorname{div}T[\Phi,\Phi]=\nabla\Phi\,\Delta\Phi=-\sigma\nabla\Phi\) with effective source \(\sigma:=-\Delta\Phi\), so \(\hat\sigma=|k|^{-2\tau}\hat{\tilde\rho}\) on \(k\neq0\). A vacuum surface and \(\mathcal{G}(S)\) are as in Assumption AV. \(\Omega\) is the unbounded component of \(\mathbb{R}^3\setminus\operatorname{supp}\tilde\rho\). Stage 1 is not edited, and no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), \(m\), \(0.08\), \(0.23\), or the \(30\,\mu\)N/kW target enters.

**(a) Existence. FAILED as a literal definition; DERIVED replacement.** Let \(\ell\) be the order of the first nonzero multipole (\(\ell=0\) charged, \(\ell=1\) neutral dipole), so \(\hat{\tilde\rho}(k)\asymp|k|^{\ell}\) near \(0\). The literal multiplier \(\hat\Phi=|k|^{-2-2\tau}\hat{\tilde\rho}\) is locally integrable at \(k=0\) iff \(\ell-2\tau>-1\), i.e. \(\tau<(1+\ell)/2\). The real-space Riesz integral \(I_{2\tau}*\Psi\), with \(\Psi\asymp r^{-1-\ell}\), converges under the same condition. So for a charged device the literal \((-\Delta)^{-\tau}A_0\) **does not exist for \(\tau\ge\tfrac12\)**, and for a dipole it does not exist for \(\tau\ge1\). That is the regime AV left open, and the literal reading is rejected there. Note also that \(\sigma=I_{2\tau}\tilde\rho\) alone exists literally only for \(2\tau<3\) when the device is charged. The replacement used below is the canonical finite-part continuation \(\Phi=\mathcal K_\tau*\tilde\rho\), where \(\mathcal K_\tau\) is the radial fundamental solution of \((-\Delta)^{1+\tau}\). For a unit charge,
\[
\Phi=C_\tau\,r^{2\tau-1},\qquad C_\tau=\frac{\Gamma(\tfrac12-\tau)}{4^{\tau}\sqrt\pi\,\Gamma(1+\tau)},
\]
with the log kernel \(\Phi=-\tfrac{2}{\pi}\ln r\) at \(\tau=\tfrac12\). This is the analytic continuation in \(\tau\), and it equals the literal reading wherever the literal reading exists. At \(\tau=1\) it gives \(\Phi=-r/2\) and \(\sigma=\Psi\) exactly. The continuation is fixed only up to polynomials. AW.2 is proved for every such choice when \(\tau\notin\mathbb Z\), and for the canonical kernel when \(\tau\in\mathbb Z\). AW.1 and AW.3 are stated for \(0<\tau<\tfrac32\), which avoids the \(r^{2j}\ln r\) kernels at \(\tau=\tfrac32,\tfrac52,\dots\).

**THEOREM AW.**

1. **AW.1 (one charge, off-centre sphere, no decay hypothesis).** Let a point charge \(q\) sit at \(a\), and let \(S\) be the sphere of radius \(R\) about \(c\), with \(b=|c-a|<R\). Then AV.4 holds without (H):
\[
\mathcal{G}(S)=D_\tau\,q^2\,J_\tau(R,b)\,\hat e_b,\qquad D_\tau=2\pi\tau(2\tau-1)^2C_\tau^2>0,\qquad J_\tau=\int_{R-b}^{R+b}r^{4\tau-3}(1-u_0^2)\,dr ,
\]
with \(u_0\) and \(\hat e_b\) as in AV.4. Here \(D_\tau\) is continuous at \(\tau=\tfrac12\), with value \(4/\pi\). Exact values: \(D_{1/4}=1\), \(D_{3/4}=\tfrac43\), \(D_1=\tfrac\pi2\), \(D_{5/4}=\tfrac{16}5\). Then \(J_{1/4}=\tfrac{4b}{3(R^2-b^2)}\), \(J_{3/4}=\tfrac{4b}{3}\), and \(J_{5/4}=\tfrac{4b(5R^2-b^2)}{15}\). Under dilation \(J_\tau=\tfrac{4b}3R^{4\tau-3}(1+O(b^2/R^2))\). A centred sphere gives \(0\).
2. **AW.2 (no surface-independent nonzero value, for every \(\tau>0\)).** For the reading \(T[\Phi,\Phi]\) or \(T[\Psi,\Phi]\), if \(\mathcal{G}\) takes the same value on every vacuum surface, that value is \(0\). Moreover, if \(\tau\notin\mathbb Z\), the self reading \(T[\Phi,\Phi]\) is surface-dependent for every nonzero device. This holds for any choice of the polynomial ambiguity. For \(\tau\notin\mathbb Z\), the cross reading \(T[\Psi,\Phi]\) is surface-dependent for every device whose potential \(\Psi\) does not vanish identically outside the device. If \(\tau\in\mathbb Z\) and the canonical kernel is used, surface independence holds iff \(\sigma\) has compact support, i.e. \(\tilde\rho=(-\Delta)^{\tau}\varphi\) with \(\varphi\) compactly supported (or, for the cross reading, iff \(\Psi\equiv0\) on \(\Omega\)), and then \(\mathcal{G}\equiv0\). No decay of \(\nabla\Phi\) is assumed anywhere.
3. **AW.3 (dilation; the marginal readings).** Let \(Q\neq0\) be the total charge and \(x_Q\) the centre of charge, and let \(0<\tau<\tfrac32\). On spheres \(S_{c,R}\),
\[
\mathcal{G}(S_{c,R})=D_\tau Q^2J_\tau(R,|c-x_Q|)\,\hat e+O(R^{4\tau-4})=\tfrac43D_\tau Q^2\,(c-x_Q)\,R^{4\tau-3}\big(1+O(R^{-2})\big)+O(R^{4\tau-4}) .
\]
So \(\mathcal{G}\to0\) for \(\tau<\tfrac34\). This includes \(\tfrac12\le\tau<\tfrac34\), where (H) fails and \(\oint|T|\,dA\sim R^{4\tau-2}\) does not tend to \(0\). For \(\tau>\tfrac34\), \(\mathcal{G}\) diverges whenever \(c\neq x_Q\). At the **marginal** value \(\tau=\tfrac34\), one charge gives exactly
\[
\mathcal{G}(S_{c,R})=\tfrac{16}{9}\,q^2\,(c-a)\qquad\text{for every }R>|c-a| ,
\]
and a charged cluster gives \(\tfrac{16}9Q^2(c-x_Q)+O(R^{-1})\). **(d) Neutral dipole \(p\) at \(0\).** Centred spheres give \(\mathcal{G}=0\) for every \(\tau\) by parity. By homogeneity \(\mathcal{G}(S_{c,R})=R^{4\tau-5}M_\tau c+O(|c|^3R^{4\tau-7})\), so the marginal value moves to \(\tau=\tfrac54\), where exactly
\[
\mathcal{G}(S_{c,R})=\tfrac{16}{75}\big(|p|^2c-8\,(p\cdot c)\,p\big)\qquad\text{for every }R>|c| .
\]
For example, at \(\tau=\tfrac34\) with \(c\parallel p\), \(\mathcal{G}=-\tfrac{4}{15}|p|^2c/(R^2-|c|^2)\).

So removing the Maxwell piece leaves a slow-tail remainder that never has the same nonzero value on every surface. Below the marginal \(\tau\) it vanishes under dilation. At the marginal \(\tau\) it is independent of the radius but linear in the offset of the audit centre, and it is \(0\) on spheres centred at \(x_Q\) (or at the dipole). Above the marginal \(\tau\) it grows without bound.

**Proof of (a).** In polar coordinates, \(\int_{|k|<1}|k|^{-2-2\tau}|k|^{\ell}d^3k\propto\int_0^1k^{\ell-2\tau}dk\), which is finite iff \(\ell-2\tau>-1\). The Riesz integrand is \(|x-y|^{2\tau-3}|y|^{-1-\ell}\), which is integrable at infinity iff \(2\tau-4-\ell<-3\). For the kernel, \(\Delta r^{p}=p(p+1)r^{p-2}\) gives \(-\Delta(C_\tau r^{2\tau-1})=-2\tau(2\tau-1)C_\tau r^{2\tau-3}\). This equals \(4\pi c_{2\tau}r^{2\tau-3}\), where \(c_\alpha=\Gamma(\tfrac{3-\alpha}2)/(2^\alpha\pi^{3/2}\Gamma(\tfrac\alpha2))\) is the Riesz constant, so \(\sigma=(-\Delta)^{-\tau}(4\pi\delta)\). The poles of \(\Gamma(\tfrac12-\tau)\) at half-integer \(\tau\) are the log cases. At \(\tau=\tfrac12\), \((2\tau-1)C_\tau\to-2/\pi\) reproduces \(-\tfrac2\pi\ln r\).

**Proof of AW.1.** \(\Phi\) is radial about \(a\), so on the centred sphere \(S_0\) of radius \(R'>R+b\) we have \(Tn=\tfrac12\Phi'^2n\), and \(\oint_{S_0}n\,dA=0\). The region between \(S\) and \(S_0\) is a compact Lipschitz region in the vacuum, where \(T\in C^1\). The divergence theorem then gives \(\mathcal{G}(S)=\mathcal{G}(S_0)-\int_{\mathrm{between}}\operatorname{div}T=-\int_{\mathrm{between}}\operatorname{div}T\). This is the AV.4 cap integral: shells with \(r>R+b\) lie wholly inside \(S_0\) and outside \(S\), and they are full, so they contribute \(0\). No limit \(R'\to\infty\) is taken, so no decay is used. With \(m=-r^2\Phi'=-qC_\tau(2\tau-1)r^{2\tau}\), we get \(mm'/r^2=q^2C_\tau^2(2\tau-1)^2\,2\tau\,r^{4\tau-3}\). The closed forms of \(J\) are integrals of \(r^{\beta}\) times the quartic \(4r^2b^2-(r^2+b^2-R^2)^2\), over \(4r^2b^2\).

**Proof of AW.2.** (i) *Surface independence forces \(\operatorname{div}T=0\) on \(\Omega\).* Take \(x_0\in\Omega\) with \(|x_0|>r_s\), where \(\operatorname{supp}\tilde\rho\subset B_{r_s}\). Let \(D_1=B(0,|x_0|)\) and \(D_2=D_1\cup B(x_0,\varepsilon)\). Both \(\partial D_1\) and \(\partial D_2\) are vacuum surfaces, and \(\mathcal{G}(\partial D_2)-\mathcal{G}(\partial D_1)=\int_{B(x_0,\varepsilon)\setminus D_1}\operatorname{div}T\). If this is \(0\) for all small \(\varepsilon\), divide by the volume and let \(\varepsilon\to0\) to get \(\operatorname{div}T(x_0)=0\). (ii) *Analyticity.* \(\mathcal K_\tau\) is real-analytic on \(\mathbb{R}^3\setminus\{0\}\), so \(\Phi=\mathcal K_\tau*\tilde\rho+P\) (\(P\) a polynomial) and \(\sigma\) are real-analytic on the connected open set \(\Omega\). The real-analytic functions on \(\Omega\) form an integral domain. So \(\sigma\,\partial_i\Phi\equiv0\) on the open set \(\{|x|>r_s\}\) gives \(\sigma\equiv0\) or \(\nabla\Phi\equiv0\) there. The second gives \(\sigma=-\Delta\Phi=0\) as well. Either way \(\sigma\equiv0\) on \(\Omega\), so \(\operatorname{supp}\sigma\) is compact. (iii) *Fourier, \(\tau\notin\mathbb Z\).* \(\hat\sigma\) and \(\hat{\tilde\rho}\) are entire, and \(\hat{\tilde\rho}(k)=|k|^{2\tau}\hat\sigma(k)\) for \(k\neq0\). Polynomials only add terms supported at \(k=0\). On a line \(k=te\), suppose \(g(t)=\hat\sigma(te)\not\equiv0\). Write \(g=t^mh\) with \(h(0)\neq0\). Then \(\hat{\tilde\rho}(te)/(t^mh(t))\) is meromorphic at \(0\) and equals \(|t|^{2\tau}\) for \(t\neq0\). That is impossible unless \(2\tau\) is an even integer: \(t^{2\tau}\) is not meromorphic for \(2\tau\notin\mathbb Z\), and \(|t|^{2\tau}\) has different one-sided Laurent series for odd \(2\tau\). So \(\hat\sigma\equiv0\), hence \(\tilde\rho=0\). (iv) *\(\tau\in\mathbb Z\), canonical kernel.* \(\sigma=\mathcal L_\tau*\tilde\rho\), where \(\mathcal L_\tau=-\Delta\mathcal K_\tau\) is the fundamental solution of \((-\Delta)^\tau\), so \((-\Delta)^\tau\sigma=\tilde\rho\). Then \(\Phi=\mathcal K_\tau*(-\Delta)^\tau\sigma=\big((-\Delta)^\tau\mathcal K_\tau\big)*\sigma=\mathcal N*\sigma\), because \((-\Delta)^\tau\mathcal K_\tau\) is radial, homogeneous of degree \(-1\), and has \(-\Delta\) equal to \(\delta\), so it is \(\mathcal N\). With \(\sigma\) compactly supported, \(|\nabla\Phi|=O(r^{-2})\), so \(\oint_{S_R}|T|=O(R^{-2})\). A surface-independent value is therefore \(\lim_R\mathcal{G}(S_R)=0\). (v) *Cross reading.* On \(\Omega\), \(\Delta\Psi=0\), so \(\operatorname{div}T[\Psi,\Phi]=-\tfrac12\sigma\nabla\Psi\). Steps (i) and (ii) give \(\sigma\equiv0\) on \(\Omega\) or \(\nabla\Psi\equiv0\) on \(\Omega\). In the second case \(T[\Psi,\Phi]=0\) on \(\Omega\) and \(\mathcal{G}\equiv0\). In the first case, (iii) gives \(\tilde\rho=0\) for \(\tau\notin\mathbb Z\), and (iv) gives \(\mathcal{G}=\lim O(R^2\cdot R^{-2}\cdot R^{-2})=0\) for \(\tau\in\mathbb Z\). Conversely, if \(\operatorname{div}T=0\) on \(\Omega\), any two vacuum surfaces bound, together with a large sphere, regions in \(\Omega\), because the exterior of a vacuum surface is connected and unbounded. So \(\mathcal{G}\) is constant.

**Proof of AW.3.** About \(x_Q\) the dipole moment vanishes. Taylor's theorem for \(\nabla\mathcal K_\tau(x-y)\) in \(y\) gives \(\nabla\Phi=Q\nabla\mathcal K_\tau(x-x_Q)+E_2\) with \(|E_2|\le C r^{2\tau-4}\) on \(S_{c,R}\). Then \(T=T_{\mathrm{mono}}+O(R^{2\tau-2}R^{2\tau-4})\), and integrating over an area of \(4\pi R^2\) gives an error \(O(R^{4\tau-4})\). The monopole term is AW.1 with \(b=|c-x_Q|\). For the dipole, \(\Phi=-p\cdot\nabla\mathcal K_\tau\) and \(\operatorname{div}T\) is homogeneous of degree \(4\tau-7\) and odd. Rescaling gives \(\mathcal{G}(S_{\lambda c,\lambda R})=\lambda^{4\tau-4}\mathcal{G}(S_{c,R})\), and \(\mathcal{G}(S_{c,R})\) is odd in \(c\), because \(\operatorname{div}T\) is odd, and smooth in \(c\) for \(|c|<R\). Expanding gives \(R^{4\tau-5}M_\tau c+O(|c|^3R^{4\tau-7})\). The exact \(\tau=\tfrac54\) formula is the cap integral of the proof of AW.1, done symbolically, with \(\mathcal{G}=0\) on the centred sphere. The map \(p\mapsto\mathcal{G}\) is quadratic and equivariant under the \(O(2)\) fixing \(c\), so it has the form \(\alpha p_\parallel^2\hat c+\beta|p_\perp|^2\hat c+\gamma p_\parallel p_\perp\). The three exact cases \(p\parallel c\), \(p\perp c\), and \(p\) at \(45^\circ\) fix \(\alpha,\beta,\gamma\), and they give \(\tfrac{16}{75}(|p|^2c-8(p\cdot c)p)\). For the one-charge marginal formula, \(J_{3/4}=4b/3\) holds identically.

**KEPT FAILURE AW.1.** Slow-tail smoothing does not open a surface-independent nonzero \(\mathcal{G}\) for an isolated static device under C5 and static massless C6. The literal \((-\Delta)^{-\tau}A_0\) does not exist for \(\tau\ge(1+\ell)/2\). Its finite-part continuation exists for every \(\tau>0\), but any value that is the same on every audit surface is \(0\) (AW.2). For non-integer \(\tau\) the self reading is surface-dependent for every nonzero device, with no decay needed: the map \(|k|^{2\tau}\) is not analytic, so the smeared charge always leaks into the vacuum. The tempting cases are the marginal readings, \(\tau=\tfrac34\) for a charged device and \(\tau=\tfrac54\) for a neutral dipole. There \(\mathcal{G}\) is finite and nonzero. For one charge or one dipole it is exactly the same on every sphere of every radius with a given centre, and a charged cluster approaches the same value under dilation: \(\tfrac{16}9Q^2(c-x_Q)\), or \(\tfrac{16}{75}(|p|^2c-8(p\cdot c)p)\). But it is linear in where the audit sphere is centred, and it vanishes on spheres centred on the device. Above the marginal \(\tau\) it grows like \(R^{4\tau-3}\) (or \(R^{4\tau-5}\)). The integer-\(\tau\) loophole is not physical. A nonzero constant value needs a polynomial added to \(\Phi\), e.g. \(\Phi+E_0\cdot x\), which gives \(\mathcal{G}=-E_0\int\sigma\). That is an external uniform background field, not an isolated device, so it falls outside Assumption AW. Interpretive step, labeled as such and the same as in Kept Failures O.1–AV.1: a force on the device does not depend on where the audit surface is drawn, so a centre-dependent \(\mathcal{G}\) cannot be read as one. This closes the slow-tail gap left open by AV. It is a failure of that route, not of Stage 2.

**Witness, computed by `scripts/slow_tail_reading_flux.py`.**
- Exact (sympy and mpmath at 40 digits):
  - \(-\Delta(C_\tau r^{2\tau-1})=4\pi c_{2\tau}r^{2\tau-3}\) at \(\tau=\tfrac14,\tfrac34,1,\tfrac54,\tfrac75\) to \(10^{-30}\). The log kernel matches at \(\tau=\tfrac12\), and \(\sigma=\Psi\) at \(\tau=1\).
  - The literal reading exists or diverges as stated in (a), at \((\tau,\ell)=(\tfrac14,0),(\tfrac12,0),(\tfrac34,0),(\tfrac34,1),(1,1),(\tfrac54,1)\).
  - \(D_\tau\) takes the four exact values, and its limit at \(\tfrac12\) is \(4/\pi\). The three closed forms of \(J\) and the \(\tfrac{4b}3R^{4\tau-3}\) asymptotics at \(\tau=\tfrac14,1,\tfrac54\) are verified.
  - The dipole cap integral at \(\tau=\tfrac54\) for \(p=e_z,e_x,e_x+e_z\) gives \(-\tfrac{112}{75}b\), \(\tfrac{16}{75}b\), and \((-\tfrac{128}{75}b,0,-\tfrac{96}{75}b)\), equal to \(\tfrac{16}{75}(|p|^2b-8(p\cdot b)p)\) for every \(R\). At \(\tau=\tfrac34\) it gives \(-\tfrac4{15}b/(R^2-b^2)\).
- One unit charge, \(160\times320\) surface quadrature of \(T\) against \(D_\tau J_\tau\), agreeing to \(10^{-9}\) relative:
  - Values at \(\tau=0.25,0.5,0.6,0.75,1,1.25\) and \((R,b)=(2,0),(2,\tfrac12),(3,1),(10,1)\), e.g. \(0.435470590524\) at \((\tfrac12;2,\tfrac12)\), \(6.14125617949\) at \((1;3,1)\), and \(425.813333333\) at \((\tfrac54;10,1)\).
  - At \(\tau=\tfrac34\): \(\tfrac{16}9b\) on radii \(2,3,10\).
- (H) fails but \(\mathcal{G}\to0\), at \(\tau=0.6\): \(\oint|Tn|=5.4546,7.1888,9.4828,12.5116\), growing like \(R^{0.4}\), while \(\mathcal{G}=0.436616,0.287590,0.189662,0.125118\) at \(R=10,20,40,80\), decaying like \(R^{-0.6}\).
- Coulomb cluster \(q=(3,-1,2)\) on the three vacuum spheres of AV:
  - \(\tau=\tfrac34\): \((14.654685,6.421449,-5.478966)\), \((31.398911,11.734110,-13.766968)\), \((-10.613254,7.039006,8.410023)\).
  - \(\tau=1\): \((34.301914,15.313051,-12.967614)\), \((84.494665,32.406204,-37.385644)\), \((-37.522389,24.004630,29.909171)\).
  - The cross reading \(T[\Psi,(-\Delta)^{-3/4}\Psi]\) gives \((5.05502,2.143256,-1.854604)\), \((9.698608,3.50335,-4.19912)\), \((-1.975102,1.324669,1.54438)\).
  - The Maxwell piece is at most \(6.1\times10^{-14}\), so removing it leaves these values.
  - At \(\tau=\tfrac34\) about \(c=(0.6,0.2,-0.3)\), \(\mathcal{G}\to\tfrac{16}9\cdot16\,(c-x_Q)=(32,12.444444,-14.222222)\), with error \(3.6\times10^{-2},9.0\times10^{-3},2.2\times10^{-3},5.6\times10^{-4}\) at \(R=10,20,40,80\). Centred at \(x_Q\) it gives \(|\mathcal{G}|=5.7\times10^{-3}\to8.9\times10^{-5}\).
  - At \(\tau=1\), \(|\mathcal{G}|=753.29,1507.62,3015.76\) at \(R=20,40,80\), growing like \(R\).
- Dipole:
  - Centred spheres give \(0\) at \(\tau=0.75,1,1.25,1.4\).
  - At \(\tau=\tfrac54\), \(p=(0.3,-0.5,0.8)\): \(\mathcal{G}=(0.28330667,-0.18645333,0.49066667)\) at \(c=(0.4,0.7,-0.2)\) on both \(R=2.5\) and \(R=30\), and twice that at \(2c\), \(R=7\), matching the closed form to \(10^{-11}\).
  - Ratios under \(R\to2R\) are \(0.4996\) at \(\tau=1\) and \(1.5161\) at \(\tau=1.4\), against \(2^{4\tau-5}\).

The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem AW does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh. It does not decide which reading the freeze intends, and it does not say whether the finite-part continuation is physically meaningful where the literal reading diverges. The quantitative parts AW.1 and AW.3 are proved for \(0<\tau<\tfrac32\), point charges, and a point dipole. The \(r^{2j}\ln r\) kernels (\(\tau=\tfrac32,\tfrac52,\dots\)) are covered only by AW.2. **CONJECTURE AW.4:** a pure \(2^\ell\)-pole has \(\mathcal{G}\sim R^{4\tau-3-2\ell}M_{\tau,\ell}c\) with marginal \(\tau_\ell=(3+2\ell)/4\). *Proved by Theorem AX below* (exact at the marginal value for every \(\ell\), with a closed form for \(M\)). Kept Failure AX.1 shows that the radius-free value need not vanish on device-centred spheres once the next moment is present. Non-Maxwell slow-tail tensors, time-dependent fields, \(d\neq3\), and \(\Psi_{\mathrm{info}}\neq A_0\) remain OPEN. A nonzero \(\mathcal{G}\) on a finite surface, or a radius-independent one, is not thrust.


## Theorem AX. Marginal slow-tail readings of a \(2^\ell\)-pole: Conjecture AW.4 holds exactly, but the radius-free value need not vanish on device-centred spheres

**ASSUMPTION AX (same family as Assumption AW, not a change to the freeze).** These are the hypotheses of Assumption AW, with the finite-part kernel \(\mathcal K_\tau=C_\tau r^{2\tau-1}\) in Maxwell form \(T[\Phi,\Phi]\), \(\operatorname{div}T=f:=\nabla\Phi\,\Delta\Phi\). A *point source of symbol \(P\)* is the compactly supported distribution \(\tilde\rho=4\pi P(\partial)\delta\), where \(P\) is a real homogeneous polynomial of degree \(\ell\). It has \(\Psi=P(\partial)\tfrac1r\) and \(\Phi=P(\partial)\mathcal K_\tau\). The dipole of AW is \(P=-p\cdot\xi\). A cluster \(\sum_iq_i\delta_{a_i}\) has \(\Phi=\sum_nP_n(\partial)\mathcal K_\tau\) with \(P_n(\xi)=\sum_iq_i(-a_i\cdot\xi)^n/n!\), so "first nonzero multipole of order \(\ell\)" means \(P_0=\dots=P_{\ell-1}=0\neq P_\ell\) (Cartesian moments, as in AW's existence criterion). Write \(P=\sum_j|\xi|^{2j}H_{\ell-2j}\) with \(H_m\) harmonic, and \(S_{c,R}\) for the sphere of centre \(c\) and radius \(R>|c|\) about a source at \(0\). Stage 1 is not edited, and no value of \(W\), \(\chi_{\mathrm{vac}}\), \(\kappa\), \(m\), \(0.08\), \(0.23\), or the \(30\,\mu\)N/kW target enters.

**THEOREM AX.**

1. **AX.1 (exact sphere formula, every \(\tau\)).** For a point source of symbol \(P\) and every \(\tau>0\) with \(\tau\notin\tfrac12+\mathbb Z\),
\[
\mathcal{G}(S_{c,R})=\int_{S^2}f(\omega)\,\frac{\rho(\omega)^{k}-R^{k}}{k}\,d\omega,\qquad k=4\tau-2-2\ell,\qquad \rho(\omega)=\omega\cdot c+\sqrt{R^2-|c|^2+(\omega\cdot c)^2},
\]
(\(\ln(\rho/R)\) at \(k=0\)). Hence \(\mathcal{G}=R^{4\tau-3-2\ell}M_\tau c+O(|c|^3R^{4\tau-5-2\ell})\) with \(M_\tau=\oint_{S^2}f\otimes n\,dA\). This is the scaling of Conjecture AW.4. Wherever the literal reading exists (\(\tau<(1+\ell)/2\)), the exponent is \(<-1\), so \(\mathcal{G}\to0\) under dilation.
2. **AX.2 (the marginal value; proves Conjecture AW.4).** At \(\tau_\ell=(3+2\ell)/4\), i.e. \(k=1\),
\[
\mathcal{G}(S_{c,R})=M\,c\qquad\text{exactly, for every }R>|c| ,
\]
so the value is radius-free, exactly linear in the audit centre, and \(0\) on spheres centred at the source. The literal \((-\Delta)^{-\tau}A_0\) never exists there: \(\tau_\ell-(1+\ell)/2=\tfrac14\) for every \(\ell\). \(\tau_\ell\) is never an integer or a half-integer, so \(\mathcal K_{\tau_\ell}\) is the pure power (no \(r^{2j}\ln r\), unique homogeneous choice).
3. **AX.3 (closed form, all \(\ell\)).** \(M=\sum_j\mathcal M_{\ell-2j}[H_{\ell-2j}]\) with no cross terms between harmonic degrees, and
\[
\mathcal M_m[H]=\frac{8}{\pi(2m+3)^2}\Big(\frac N2\,I-\frac{4}{4m^2-1}\,B\Big),\qquad N=\oint_{S^2}H^2,\quad B=\oint_{S^2}\nabla H\otimes\nabla H,\qquad \operatorname{tr}\mathcal M_m[H]=-\frac{4N}{\pi(2m-1)(2m+3)} .
\]
It recovers AW: \(\tfrac{16}9q^2I\) (\(m=0\)) and \(\tfrac{16}{75}(|p|^2I-8pp^{\mathsf T})\) (\(m=1\)). **Quadrupole**, \(P=\xi\cdot A\xi\), \(A_0=A-\tfrac13(\operatorname{tr}A)I\), at \(\tau_2=\tfrac74\):
\[
M=\tfrac{16}{81}(\operatorname{tr}A)^2I+\tfrac{32}{735}\operatorname{tr}(A_0^2)\,I-\tfrac{512}{2205}A_0^2 ,
\]
e.g. \(\operatorname{diag}(64,64,-1472)/2205\) for \(2z^2-x^2-y^2\). The trace part behaves exactly like a charge: \(P=|\xi|^2\), whose Coulomb field vanishes outside the point (\(\Psi=-4\pi\delta\)), gives the unit-charge value \(\tfrac{16}9I\). So the marginal value is not a functional of the exterior Coulomb field.
4. **AX.4 (when \(M\neq0\)).** For \(\ell\ge1\), \(M\neq0\) whenever \(P\) has no \(|\xi|^\ell\) part, in particular if \(P\) is harmonic or \(\ell\) is odd: every trace in AX.3 with \(m\ge1\) is negative, and only \(m=0\) is positive. \(M\neq0\) also for every nonzero quadrupole. \(M\) can be singular: for the harmonic quadrupole \(x^2+ty^2-(1+t)z^2\), \(t=(\sqrt{69}-13)/10\), \(M_{zz}=0\). \(M\) can vanish: for \(P_*=H_4+\tfrac{4\sqrt{11}}{385}|\xi|^4\), \(H_4=\sum\xi_i^4-\tfrac35|\xi|^4\), \(M=0\), so \(\mathcal{G}=0\) on every sphere enclosing the point, for every centre and radius. On cubes, however, \(\mathcal{G}(\lambda Q)=\lambda\,\mathcal{G}(Q)\neq0\).
5. **AX.5 (next moment).** For \(P_\ell+P_{\ell+1}\) at \(\tau_\ell\), exactly,
\[
\mathcal{G}(S_{c,R})=M c+g_0+\tfrac12\ln\!\big(1-|c|^2/R^2\big)\,v+\frac{M'c}{R^2-|c|^2},\qquad g_0=\oint_{S^2}T[\Phi_\ell,\Phi_{\ell+1}]_{\mathrm{sym}}\,n ,
\]
with \(v=\oint_{S^2}f_{\ell,\ell+1}\) and \(M'=\oint_{S^2}f_{\ell+1,\ell+1}\otimes n\). Translating the device by \(a\) sends \(g_0\mapsto g_0-Ma\). **DERIVED** (multipole series plus the exact pair formulas, verified numerically below): a compactly supported device with first nonzero moment \(P_\ell\) has \(\mathcal{G}(S_{c,R})=Mc+g_0+O(R^{-2})\). When \(M\) is invertible this is \(M(c-x_*)\) with a reading-dependent centre \(x_*=-M^{-1}g_0\).

**Proof of AX.1–AX.2.** \(\Phi\) is homogeneous of degree \(2\tau-1-\ell\) with parity \((-1)^\ell\), so \(T\) is even and \(f\) is odd and homogeneous of degree \(4\tau-5-2\ell\), real-analytic off \(0\). On \(S_{0,R}\) the integrand \(Tn\) is odd, so \(\mathcal{G}(S_{0,R})=0\). Both balls \(B_{c,R}\) and \(B_{0,R}\) contain \(0\) in their interior and are star-shaped about it. The divergence theorem on \(B_{c,R}\setminus B_{0,R}\) and \(B_{0,R}\setminus B_{c,R}\), which avoid \(0\), gives \(\mathcal{G}(S_{c,R})=\int_{S^2}\int_R^{\rho(\omega)}f(r\omega)r^2\,dr\,d\omega\). Homogeneity gives the formula. Next, \(\rho(\omega)\rho(-\omega)=R^2-|c|^2\) and \(\rho(\omega)-\rho(-\omega)=2\omega\cdot c\). Since \(f\) is odd, only the odd part of \(\rho^k\) survives, and that is \(kR^{k-1}\omega\cdot c+O(|c|^3R^{k-3})\). It is exactly \(\omega\cdot c\) when \(k=1\). The existence statements are Theorem AW (a) with \(\tau=\tau_\ell\).

**Proof of AX.3.** For harmonic \(H\) of degree \(m\), Hobson's formula gives \(H(\partial)g(r)=H(x)(r^{-1}\tfrac d{dr})^mg\), so \(\Phi=A\,H(x)r^{s}\) with \(s=\tfrac12-m\) and \(A=C_{\tau_m}c_m\), \(c_m=\prod_{k<m}(m+\tfrac12-2k)\). Then \(\Delta\Phi=As(m+\tfrac32)Hr^{s-2}\). On \(S^2\), \(f_i=\kappa H(\partial_iH+sHx_i)\) with \(\kappa=A^2s(m+\tfrac32)\). The ball identity \(\oint gx_j=\frac1{d+2}\oint\partial_jg\) (\(g\) homogeneous of degree \(d\)), plus orthogonality of \(H\) to the degree-\((m-2)\) harmonics \(\partial_{ij}H\), gives \(\oint H\partial_iHx_j=B_{ij}/(2m+1)\) and \(\oint H^2x_ix_j=(2B_{ij}/(2m+1)+N\delta_{ij})/(2m+3)\). So \(M=\frac{\kappa}{2m+3}\big(\frac4{2m+1}B+(\tfrac12-m)NI\big)\). Writing \(\Gamma(1+\tau_m)=\Gamma(\tfrac12-\tau_m)\prod_{j=0}^{m+1}(j-\tfrac{2m+1}4)\) gives \(C_{\tau_m}c_m=(-1)^m4\sqrt2/(\sqrt\pi(2m-1)(2m+3))\), hence \(\kappa=-8/(\pi(2m-1)(2m+3))\), which is the stated form. For non-harmonic \(P\), \(-\Delta\mathcal K_\tau=\mathcal K_{\tau-1}\) gives \(|\xi|^{2j}H(\partial)\mathcal K_{\tau_\ell}=(-1)^jH(\partial)\mathcal K_{\tau_{\ell-2j}}\). So each piece is a harmonic \(2^{m}\)-pole at its own marginal value. For two pieces of degrees \(m\) and \(m'\) with \(|m-m'|\ge4\), every cross integral is a product of spherical harmonics of different degree and vanishes. For \(m'=m-2\), the two surviving integrals are \(W/(2m-1)\) and \(W/((2m+1)(2m-1))\), with \(W=\oint H_{m-2}\partial_{ij}H_m\). Their coefficient is \(\lambda_{m'}(2m+1+s_m)+\lambda_ms_{m'}=(\tfrac52-m)(m+\tfrac32)[(m-\tfrac12)+(\tfrac12-m)]=0\), where \(\lambda_m=s_m(m+\tfrac32)\).

**Proof of AX.4–AX.5.** Without an \(|\xi|^\ell\) part, every \(m\) in the decomposition is \(\ge1\), so \(\operatorname{tr}M\) is a sum of negative terms. For \(\ell=2\), rotate \(A_0\) to diagonal form, so \(M\) is diagonal. The Gröbner basis of the three diagonal entries contains \(\beta^4\) and \(2a_1^2-35\beta^2\), so \(M=0\) forces \(P=0\). The singular and vanishing witnesses are exact sympy evaluations. In AX.5, the pair \((\ell,\ell+1)\) has even \(f\) of degree \(-3\), so \(\int_R^\rho r^{-1}dr=\ln(\rho/R)\), and the even part of \(\ln\rho\) is \(\tfrac12\ln(R^2-|c|^2)\). The pair \((\ell+1,\ell+1)\) has odd \(f\) of degree \(-4\), and the odd part of \(-1/\rho\) is \(\omega\cdot c/(R^2-|c|^2)\). The centred-sphere value of the \((\ell,\ell+1)\) pair is \(g_0\), which is radius-free by homogeneity (degree \(-2\) times area \(R^2\)).

**KEPT FAILURE AX.1.** Conjecture AW.4 is true. But the reading AW.1 drew from it, that the radius-free marginal value is linear in the audit centre and \(0\) on spheres centred on the device, is not a property of the marginal readings in general. **Witness:** the point source \(\tilde\rho=4\pi\big(P_*(\partial)+\partial_x^5\big)\delta\) at \(\tau=\tfrac{11}4\). Its exterior Coulomb field is that of the cubic hexadecapole \(H_4\) plus \(\partial_x^5\tfrac1r\), because the \(|\xi|^4\) part gives \(\Delta^2\tfrac1r\), which is supported at the point. On every sphere,
\[
\mathcal{G}(S_{c,R})=g_0+\frac{M'c}{R^2-|c|^2},\qquad g_0=\Big(\frac{64\sqrt{11}}{8085}-\frac{5888}{266805}\Big)e_x\approx0.004185497113\,e_x .
\]
So every sphere centred on the device gives the same nonzero \(g_0\) at every radius, and spheres of every centre tend to \(g_0\). Centred cubes give a different radius-free value, \((-0.0238957422,0,0)\). Off-centre dilated cubes give \(\lambda\,\mathcal{G}[P_*]+B+C/\lambda\), which grows without bound. So \(\mathcal{G}\) still moves with the surface, and AW.2 stands. What fails is the sphere-only argument: "not a force" must now be read off non-spherical surfaces, since on spheres the marginal reading can mimic a centre-free, radius-free vector. Interpretive step, labeled as such and the same as in Kept Failures O.1–AW.1: a force on the device does not depend on the audit surface, so a value that differs between spheres and cubes cannot be read as one. Not thrust.

**Witness, computed by `scripts/multipole_marginal_flux.py`.**
- Exact (sympy):
  - \(\tau_\ell-(1+\ell)/2=\tfrac14\). \(\tau_\ell\notin\tfrac12\mathbb Z\) for \(\ell<40\). \(-\Delta\mathcal K_\tau=\mathcal K_{\tau-1}\) at \(\tau=\tfrac74,\tfrac94,\tfrac{11}4,\tfrac{13}4\). \(C_{\tau_\ell}c_\ell\) and \(\kappa_\ell=-8/(\pi(2\ell-1)(2\ell+3))\) for \(\ell=0,\dots,8\).
  - The sphere integral \(M\) equals the AX.3 closed form for random harmonic \(H\) of degree \(1\)–\(5\), and recovers \(\tfrac{16}9I\) and \(\tfrac{16}{75}(|p|^2I-8pp^{\mathsf T})\).
  - Quadrupole: \(\operatorname{diag}(64,64,-1472)/2205\) for \(2z^2-x^2-y^2\) and \(\operatorname{diag}(-\tfrac{16}{441},-\tfrac{16}{441},\tfrac{16}{735})\) for \(xy\). The general six-coefficient formula holds, \(P=|\xi|^2\) gives \(\tfrac{16}9I\), and the Gröbner basis contains \(\beta^4\).
  - Harmonic degrees decouple at \(\ell=3,4,5\).
  - Singular quadrupole \(M_{zz}=0\). \(M[H_4]=-\tfrac{256}{121275}I\) and \(M[P_*]=0\).
  - \(g_0\) and \(v=0\) of the Kept Failure witness. An axial quadrupole plus octupole \(x^3-3xy^2+z^3\) gives \(g_0=(0,0,-0.6965986395)\).
- Quadrature (\(160\times320\) nodes):
  - \(\mathcal{G}=Mc\) to \(2\times10^{-13}\) relative on \(R=1.5,3,6,40\) for \(2z^2-x^2-y^2\), \(x^2\), a random harmonic quadrupole, a random harmonic octupole, and a non-harmonic octupole. For example, \(\mathcal{G}=(0.011609977,-0.008707483,-0.333786848)\) at \(c=(0.4,-0.3,0.5)\) for \(2z^2-x^2-y^2\).
  - Off the marginal value (\(\tau=\tfrac85,\tfrac{19}{10}\)), the AX.1 integral matches to \(10^{-13}\), and \(R\to2R\) ratios are \(0.6595,1.5159\) against \(2^{4\tau-7}\).
  - The singular quadrupole gives \(0\) on spheres centred on the \(z\)-axis.
  - \(P_*\) gives \(|\mathcal{G}|\le2\times10^{-16}\) on three off-centre spheres, and \((-0.00031154,0.00020248,-0.00059395)\times\lambda\) on cubes.
  - The Kept Failure witness gives \(g_0\) on centred spheres \(R=1,3,30\) and the exact AX.5 formula to \(10^{-14}\) on off-centre spheres. At centre \((-3,2,1)\), \(R=600\), it gives \((0.00418760,1.7\times10^{-6},8.5\times10^{-7})\). Centred cubes give \(-0.0238957422\,e_x\) at every size.
- Point-charge clusters at \(\tau=\tfrac74\):
  - The linear quadrupole \((1,-2,1)\) gives \(\mathcal{G}\to Mc\) with error \(2.0\times10^{-5},5.0\times10^{-6},1.2\times10^{-6}\) at \(R=5,10,20\).
  - An \(8\)-charge neutral, dipole-free device with an octupole gives \(\mathcal{G}\to Mc+g_0\), with \(g_0=(0.00483962,-0.00249034,0.00205322)\) and error \(\sim R^{-2}\). Spheres centred at \(x_*=(-0.195235,0.512464,-0.007802)\) give \(|\mathcal{G}|\to0\), and spheres centred at the origin do not.

The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem AX does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh, and it does not choose a reading. The device asymptotics in AX.5 are DERIVED (proof sketch plus two numerical devices), not a full theorem with explicit constants. Non-marginal \(\tau\) is covered only by AX.1 for pure point sources. **OPEN:** a classification of the symbols with \(M=0\) for even \(\ell\ge4\) (they need an \(|\xi|^\ell\) part), a full theorem with constants for AX.5 on extended devices, the \(r^{2j}\ln r\) kernels; non-Maxwell slow-tail tensors; time-dependent fields; \(d\neq3\); and \(\Psi_{\mathrm{info}}\neq A_0\). A radius-free \(\mathcal{G}\) on spheres, even a centre-free one, is not thrust.

## Theorem BA. The Maxwell-form reading on the frozen \(0.45\) mesh: the action–reaction piece integrates to \(0\), the residual does not, and the residual is a mesh-symmetry defect

This is the deferred Stage 2 question, asked on the actual discrete mesh rather than in the continuum. It is one precise reading, labeled as an assumption. The freeze does not choose it.

**ASSUMPTION BA (mesh reading, not a change to the freeze).**

1. **Mesh.** The level-\(n\) combinatorial gasket \(V_n\) of `build_gasket(n)` (as in [GASKET.md](GASKET.md)), with corners \(C=\{c_0,c_1,c_2\}\), positions \(P_i\in\mathbb R^2\), integer Laplacian \(L_n\), and \(A=L_n^{\alpha}\) by spectral calculus with \(0^{\alpha}=0\) (GASKET Assumption 4), at \(\alpha=0.45\).
2. **Field and charge.** \(\Psi=u\) is the \(A\)-Dirichlet extension of corner data \(d\in\mathbb R^3\), so \((Au)_i=0\) for \(i\notin C\). The discrete charge is \(q=Au\). It is supported on \(C\) and neutral (Theorem H).
3. **Field at a corner.** \(g_c\) is the gradient of the affine interpolant of \(u\) on the unique level-\(n\) cell at \(c\). The corner has exactly the two neighbours of that cell. Set \(E_c=-g_c\).
4. **Discrete divergence.** For the quadratic reading \(T=\nabla\Psi\otimes\nabla\Psi-\tfrac12|\nabla\Psi|^2\delta\) the continuum divergence is \(\partial_jT^{ij}=\partial_i\Psi\,\Delta\Psi=\rho E_i\) (Theorem N). On the mesh it is \(f_i:=q_iE_i=-q_ig_i\), which is \(0\) off \(C\). The integrated divergence over a vertex set \(S\) is \(\mathcal G(S):=\sum_{i\in S}f_i\), and \(D:=\mathcal G(V_n)\).
5. **Maxwell piece.** Let \(G=A^{+}\), the symmetric Green kernel with \(G\mathbf 1=0\). Then \(u=Gq+\bar u\mathbf 1\). The pair force of \(k\) on \(c\) is \(\varphi_{ck}:=-q_cq_kH_{ck}\), with \(H_{ck}:=g_c[Ge_k]\). The Maxwell (action–reaction) piece is \(\varphi^{M}_{ck}=\tfrac12(\varphi_{ck}-\varphi_{kc})\). The residual is \(\varphi^{R}_{ck}=\tfrac12(\varphi_{ck}+\varphi_{kc})\), self terms \(\varphi_{cc}\) included.

Stage 1 is not edited. No value of \(W\), \(\chi_{\mathrm{vac}}\), \(0.08\), \(0.23\), or the \(30\,\mu\)N/kW target enters.

**THEOREM BA.**

1. **BA.1 (exact identities; every \(n\), every \(\alpha>0\), every \(d\), every linear gradient rule that kills constants).** \(\sum_{c,k}\varphi^{M}_{ck}=0\). For every \(S\), \(\sum_{c\in S}\sum_k\varphi^{M}_{ck}=\sum_{c\in S\cap C,\,k\in C\setminus S}\varphi^{M}_{ck}\) (cut form, as in Theorem F). So the Maxwell piece has integrated divergence \(0\) over the whole mesh, and the integrated residual is \(\sum\varphi^{R}=D\). Also \(\mathcal G(S)=\sum_{c\in S\cap C}f_c\), which depends on \(S\) only through \(S\cap C\). In particular \(\mathcal G(S)=D\) for every \(S\supseteq C\).
2. **BA.2 (one scalar per level, by Schur).** \(D_n(d)=\kappa_n(\alpha)\,Q(d)\), where \(Q(d):=\sum_c(d_c-\bar d)^2(P_c-\bar P)\) and \(|Q(d)|=|d-\bar d\mathbf 1|^2/(2\sqrt3)\). So \(D\neq0\) for every nonconstant \(d\) if and only if \(\kappa_n(\alpha)\neq0\).
3. **BA.3 (certified at \(\alpha=0.45\)).** \(\kappa_1=-4.0427396046867839149\ldots\), \(\kappa_2=-6.6139165202761608886\ldots\), \(\kappa_3=-10.738239541819965987\ldots\), \(\kappa_4=-18.526339511158477727\ldots\). The enclosure widths are \(\le3\cdot10^{-37}\) (mpmath.iv). So at levels \(1\)–\(4\) of the \(0.45\) mesh, the integrated divergence after the Maxwell piece is removed is **nonzero for every nonconstant corner data**. At level \(4\) with \(d=(1,-\tfrac12,-\tfrac12)\), \(D=(6.9473773166844291477\ldots,\;4.0110701639496550916\ldots)\). It is the same for every vertex set containing all three corners. It is surface-dependent otherwise: the certified sets with \(S\cap C=\{c_0\}\) or \(\{c_0,c_1\}\) give other values. The self (diagonal) part alone is \(\kappa^{\mathrm{self}}_4=-18.530320559786\ldots\), so the symmetric cross part is only \(+0.00398\).
4. **BA.4 (\(\alpha=1\) control, exact).** \(\kappa_n(1)=-9\big((18/25)^n-(6/25)^n\big)\) for every \(n\ge1\), e.g. \(-108/25\), \(-2592/625\), \(-50544/15625\). The harmonic mesh has a nonzero residual at each level, and it tends to \(0\).
5. **BA.5 (translation-invariant control).** Take a periodic lattice whose Green kernel is translation-invariant and even, \(G_{xy}=\mathsf g(x-y)\) with \(\mathsf g(-z)=\mathsf g(z)\), and a gradient stencil that is odd, \(S(z):=\mathrm{grad}[\mathsf g(\cdot+z)](0)\) with \(S(-z)=-S(z)\). Then \(H_{ck}=-H_{kc}\) and \(H_{cc}=0\). So \(\varphi=\varphi^{M}\) and \(D=0\) for every device and every data set.

**Proof of BA.1.** \(\varphi^{M}\) is antisymmetric, so its sum over all ordered pairs is \(0\), and pairs with both ends in \(S\) cancel. \(g_c\) is linear and kills constants, so \(f_c=-q_cg_c[Gq]=\sum_k\varphi_{ck}\). Hence \(\sum\varphi^{R}=\sum\varphi-\sum\varphi^{M}=D\). Off \(C\), \(q_i=0\), so \(f_i=0\).

**Proof of BA.2.** Replacing \(d\) by \(d+t\mathbf 1\) replaces \(u\) by \(u+t\mathbf 1\): \(A\mathbf 1=0\), and the Dirichlet extension is unique because \(A_{II}\succ0\). This leaves \(q\) and \(g\) unchanged, so \(D\) is a quadratic form in \(x=d-\bar d\mathbf 1\in V=\{\sum x_c=0\}\). Each symmetry \(\sigma\in S_3\) of the gasket is an isometry \(T_\sigma\) with linear part \(R_\sigma\). It fixes \(\bar P\) and permutes vertices and level-\(n\) cells. Since \(u\circ T_\sigma^{-1}\) is the extension of \(d\circ\sigma^{-1}\), \(q\) moves as a scalar and \(g_{\sigma c}[u\circ T_\sigma^{-1}]=R_\sigma g_c[u]\). So \(D(\sigma x)=R_\sigma D(x)\). The plane carries the standard representation \(V\) of \(S_3\), and \(\mathrm{Sym}^2V\cong\mathbf 1\oplus V\) (character \((3,1,0)\)). So \(\mathrm{Hom}_{S_3}(\mathrm{Sym}^2V,V)\) is one-dimensional. \(Q\) is a nonzero element, and \(|Q|^2=|x|^4/12\) is a polynomial identity (sympy).

**Proof of BA.3.** Every step is exact or outward-rounded:
- An integer Krylov sequence gives a monic \(m\in\mathbb Z[x]\) with \(m(L_n)e_j=0\) for every basis vector \(e_j\), checked in exact integers. So \(\mathrm{spec}(L_n)\) lies in the roots of \(m\). \(m\) is squarefree, of degree \(4,8,18,36\) at \(n=1,\ldots,4\).
- Its roots are all real. sympy isolates them in rational intervals of width \(10^{-70}\).
- \(A\) and \(A^{+}\) are \(p_\pm(L_n)\), where \(p_\pm\) interpolates \(x^{\pm\alpha}\) (and \(0\) at \(0\)) on those roots. This is exact on the spectrum. The divided differences are taken in mpmath.iv at \(120\) digits and multiplied into the exact integer columns \(L_n^te_j\).
- The Dirichlet system is solved with a Krawczyk enclosure (\(\|I-RA_{II}\|_\infty<1\)).
- The script then checks four things on the enclosures: \((Au)_I\ni0\), \(\sum q\ni0\), \(\sum\varphi^{M}\ni0\), and \(\sum\varphi^{R}\cap D\neq\emptyset\). The last one is a nontrivial check of \(u=Gq+\mathrm{const}\).

**Proof of BA.4.** At \(\alpha=1\), \(q_c=(Lu)_c\) uses only the two neighbours in the corner cell. The values on that cell are \(M^n(d_c,d_{c'},d_{c''})\) with the harmonic-extension matrix \(M=\begin{psmallmatrix}1&0&0\\2/5&2/5&1/5\\2/5&1/5&2/5\end{psmallmatrix}\) and \(h=2^{-n}\). Diagonalizing \(M\) gives \(D_x/Q_x\) as an exponential polynomial in \(n\), and sympy reduces it to the stated closed form symbolically. Exact rational solves at \(n=1,2,3\) agree.

**Proof of BA.5.** \(S(-z)=\sum_e c_e\mathsf g(-z+e)=\sum_e c_e\mathsf g(z-e)=\sum_e c_{-e}\mathsf g(z+e)\). With an odd stencil, \(c_{-e}=-c_e\), so this is \(-S(z)\). Then \(H_{ck}=S(P_c-P_k)\).

**KEPT FAILURE BA.1 (a nonzero residual on the mesh, but a symmetry defect of the discrete pair force, not a Stage 2 signal).** BA.3 answers the deferred question under Assumption BA: yes, \(\sum\mathrm{div}\neq0\) after the Maxwell piece is removed. Reading that as an informational residual fails, for three reasons.

**Witness (i).** Take the same operator \(L^{0.45}\) and the same one-sided cell-gradient rule on the translation-invariant \(6\times6\) triangular torus, which has no boundary and no corners. The device is \((0,0),(2,0),(0,2)\). It gives the certified \(D=(2.8976275758051\ldots,\;1.6729460609024\ldots)\) for \(d=(1,-\tfrac12,-\tfrac12)\) and \((2.2840764921345\ldots,\;1.2310526376196\ldots)\) for \((1,-\tfrac12,0)\). The centred six-neighbour gradient on the same torus gives \(D\ni0\) with width below \(10^{-75}\), also for a scalene device, as BA.5 requires. So the nonzero value is produced by a stencil that is not odd: the discrete pair force breaks Newton's third law. At a gasket corner no odd stencil exists, because the corner has two neighbours in a \(60^\circ\) wedge.

**Witness (ii).** On the gasket the residual is almost all self-force: \(\kappa^{\mathrm{self}}_4=-18.53032\) against \(\kappa_4=-18.52634\).

**Witness (iii).** There is no refinement-stable value. \(|\kappa_n|\) strictly increases at every certified level and in floats through level \(6\), roughly like \(2^n\): \(2^{-n}\kappa_n=-2.0214,-1.6535,-1.3423,-1.1579,-1.0488,-0.9818\). At \(\alpha=1\) it decays instead (BA.4).

Interpretive step, labeled as such and the same as in Kept Failures O.1–AX.1: a quantity that appears on a boundaryless translation-invariant mesh as soon as the stencil changes, and that grows like the inverse mesh spacing, cannot be read as a force on the device. Not thrust.

**CONJECTURE BA.2.** \(2^{-n}\kappa_n(0.45)\) converges to a nonzero limit. Aitken on float levels \(5\)–\(7\) gives about \(-0.86\); this is not certified. The observed mechanism: the corner jump \(u_{c_0}-\tfrac12(u_p+u_r)\) for \(d=(1,-\tfrac12,-\tfrac12)\) runs \(0.942,0.763,0.676,0.626,0.596,0.576,0.564\) for \(n=1,\ldots,7\), apparently toward a positive limit, and \(q_{c_0}\) also stays positive, so \(g_c\sim2^n\). Falsified by a certified \(2^{-n}\kappa_n\) sequence tending to \(0\) or to \(-\infty\). **Status: proved for every \(\alpha<\log3/\log5\) by Theorem BB.** The limit is \(-4\sigma_\infty/R_\infty^2\), and \(K(0.45)\in[-0.854646905120924,-0.854646905120848]\), not \(\approx-0.86\). It is \(0\) for \(\alpha\ge\log3/\log5\).

**NOT THIS.** Theorem BA does not decide which reading the freeze intends. Other gradient rules keep BA.1 and BA.2 (the identities and the Schur form) but change \(\kappa_n\). BA.3 is certified at levels \(1\)–\(4\) only; levels \(5\)–\(7\) are floats. Using the combinatorial current \(Lu\) as the charge is not covered. The device is the Dirichlet corners, so the mesh has no "vacuum surface" around the device short of containing \(C\). A nonzero \(D\) is not thrust.

**Witness, computed by `scripts/mesh_maxwell_residual.py`.** Exact: Krylov minimal polynomials annihilate \(L_n\) (\(n\le4\)) and the torus Laplacian. \(|Q|^2=|x|^4/12\). \(\kappa_n(1)=-9((18/25)^n-(6/25)^n)\) symbolically, plus exact rational solves at \(n=1,2,3\). Intervals (mpmath.iv, \(120\) digits, Krawczyk): \(\kappa_1,\ldots,\kappa_4\) above, all certified \(<0\) and strictly increasing in modulus. For each level the script checks interior harmonicity, neutrality, \(\sum\varphi^{M}\ni0\), residual total \(=D\), \(D\parallel Q\), and that sets with one or two corners give values other than \(D\). It also certifies the torus one-sided \(D\neq0\) and torus centred \(D\ni0\) (two data sets plus a scalene device). Floats at levels \(1\)–\(6\) check the Schur form with a second data set to \(10^{-12}\) relative and match the certificates to \(10^{-9}\). The script exits nonzero if any of these figures moves.

## Theorem BB. The mesh residual constant decimates exactly: Conjecture BA.2 holds, the limit is \(-4\sigma_\infty/R_\infty^2\), and all of it is self-force

Theorem BA left one open number: how \(\kappa_n(0.45)\) behaves under refinement. This theorem expresses \(\kappa_n\) through the corner resolvent of [GASKET.md](GASKET.md) Theorems AQ–AS, for which exact decimation is already proved.

**SETTING.** Assumption BA, level \(n\ge1\), \(\alpha\in(0,1)\), \(G=(L_n^{\alpha})^{+}\). The notation is that of GASKET Theorems AQ–AS: \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\) with \(v=e_a-e_b\), \(R_n:=R_\alpha(n)=v^{T}Gv\), \(\sigma_n=\tfrac12v^{T}L_n^{1-\alpha}v\), \(p(s)=s(s+5)\), and \(\gamma=\log3/\log5\). For a corner \(c\) with the other corners \(c',c''\), let \(p',p''\) be its two neighbours, \(p'\) on the edge towards \(c'\). Set
\[
Y_n=G_{p'c'}-G_{p'c''},\qquad E_n=G_{p'c'}+G_{p'c''}-G_{cc'}-G_{cc''}.
\]
Split \(D\) into the self part \(\sum_c\varphi_{cc}=\kappa^{\mathrm{self}}_nQ(d)\) and the cross part \(\sum_{c\neq k}\varphi_{ck}=\kappa^{\mathrm{cross}}_nQ(d)\). Both have the Schur form, by BA.2's argument. **Assumptions:** GASKET Theorems AQ (1)–(3) (decimation, Stieltjes formulas, monotonicity in \(n\), \(p(s)r_n(p(s))\le1\)), AR (3)–(4) (\(\sigma_\infty\ge2\cdot4^{-\alpha}\); \(R_\infty<\infty\) for \(\alpha<\gamma\), \(R_n\to\infty\) for \(\alpha\ge\gamma\)), and AS (1), (4) (limit integrals and the interpolation sandwich, re-run live at \(0.45\)), all as stated, plus Theorem J's entrywise positivity of \((L_n+s)^{-1}\). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BB.**

1. **BB.1 (exact decomposition; every \(n\ge1\), every \(\alpha\in(0,1)\)).** \(q_c=2x_c/R_n\) with \(x=d-\bar d\mathbf 1\), and
   \[
   \kappa_n=2^n\frac{4}{R_n^2}\,(3Y_n-\sigma_n),\qquad
   \kappa^{\mathrm{self}}_n=-2^n\frac{4}{R_n^2}\,(\sigma_n-E_n),\qquad
   \kappa^{\mathrm{cross}}_n=2^n\frac{4}{R_n^2}\,(3Y_n-E_n).
   \]
2. **BB.2 (tangential decimation).** \(y_n(s):=[(L_n+s)^{-1}]_{p'c'}-[(L_n+s)^{-1}]_{p'c''}=r_n(s)\prod_{j=0}^{n-1}\bigl(5+p^{j}(s)\bigr)^{-1}\). Hence \(Y_n=\frac{\sin\pi\alpha}{\pi}\int_0^\infty s^{-\alpha}y_n\,ds\) satisfies \(0<Y_n\le5^{-n}R_n/2\).
3. **BB.3 (far-field decimation).** The resolvent version of \(E_n\) is \(e_n(s)=s\,g_n(s)\), where \(g_n(s)=[(L_n+s)^{-1}]_{cc'}\). It obeys
   \[
   g_{n+1}(s)=\frac{(2+s)(5+s)(6+s)}{\bigl(6+s-2t_n(ps)\bigr)\bigl(6+s-2p(s)r_n(ps)\bigr)}\,g_n(p(s)),\qquad g_0=\frac1{s(3+s)},
   \]
   with \(t_n(s)=s[(L_n+s)^{-1}\mathbf 1_C]_c\in(0,1]\), \(t_{n+1}(s)=(2+s)t_n(ps)/(6+s-2t_n(ps))\), \(t_0=1\). Consequently \(g_n(s)\le(3/20)^n/s^2\), \(e_n\le\tfrac13\), and
   \[
   0<E_n\le\frac{\sin\pi\alpha}{\pi}\Bigl(\frac{3}{20}\Bigr)^{n(1-\alpha)}\Bigl(\frac1{3(1-\alpha)}+\frac1\alpha\Bigr).
   \]
4. **BB.4 (Conjecture BA.2 holds; every \(\alpha<\gamma\)).**
   \[
   \lim_{n\to\infty}2^{-n}\kappa_n(\alpha)=K(\alpha):=-\frac{4\sigma_\infty(\alpha)}{R_\infty(\alpha)^2}\le-\frac{8\cdot4^{-\alpha}}{R_\infty(\alpha)^2}<0 .
   \]
   At \(\alpha=0.45\) the certificate gives \(K\in[-0.854646905120924,\,-0.854646905120848]\), i.e. \(K(0.45)=-0.85464690512\ldots\), from \(R_\infty\in[2.51594623787707,2.51594623787714]\) and \(\sigma_\infty\in[1.35247562325228,1.35247562325232]\). For \(\alpha\ge\gamma\) the limit is \(0\), because \(R_n\to\infty\). Then BB.1 gives \(\kappa_n(\gamma)\sim-4\sigma_\infty2^n/(c_R(\gamma)n)^2\), and, for \(\alpha>\gamma\), \(\kappa_n\sim-4\sigma_\infty\bigl((\rho-1)/c_R\bigr)^2(2/\rho^2)^n\) with \(\rho=5^{\alpha}/3\). At \(\alpha=1\) that is BA.4's \(-9(18/25)^n\).
5. **BB.5 (the growth is all self-force; every \(\alpha<\alpha_*:=1-\log2/\log(20/3)\approx0.6346\)).** \(\kappa^{\mathrm{cross}}_n\to0\), with
   \[
   |\kappa^{\mathrm{cross}}_n|\le\frac{4}{R_1^2}\Bigl[\tfrac32(2/5)^nR_\infty+2^nE_n\Bigr]=O\bigl((2\cdot(3/20)^{1-\alpha})^n\bigr).
   \]
   At \(0.45\) the base is \(2\cdot0.15^{0.55}\approx0.704\). So \(2^{-n}\kappa^{\mathrm{self}}_n\to K\) as well, and the self term carries both the \(2^n\) growth and the whole nonzero limit.

**DERIVED (rate, not certified).** Since \(\sigma_n\) and \(Y_n\) converge faster than \(R_n\) (AR (2): \(\mu^n\), and \(5^{-n}\)), AR (4a) gives \(2^{-n}\kappa_n-K\sim-\frac{8\sigma_\infty c_R}{R_\infty^3(1-\rho)}\rho^n\). With AT's \(c_R(0.45)\) the constant is \(-1.083020\), and the recursion at \(n=40\) reproduces it to six digits.

**Proof of BB.1.**
- *Charge.* By \(S_3\) symmetry, the corner block of \(G\) acts on the zero-sum plane as \(R_n/2\), since \(v^{T}Gv=R_n\) and \(|v|^2=2\). From \(u=Gq+\bar u\mathbf 1\) with \(q\) zero-sum, \(x=(R_n/2)q\).
- *Combinatorial current.* \(LG=L^{1-\alpha}\) on \(\mathbf 1^{\perp}\), so \(I:=(Lu)_C=(L^{1-\alpha}q)_C\). On the zero-sum plane this block is \(\sigma_n\) (AQ (2)), so \(I_c=\sigma_nq_c\).
- *Gradient.* The corner cell is equilateral with side \(h=2^{-n}\). Split the gradient along the altitude and its normal: \(g_c=2^n\bigl[I_c(P_c-\bar P)+(u_{p'}-u_{p''})(P_{c'}-P_{c''})\bigr]\) (BA's unit triangle).
- *Tangential increment.* The reflection fixing \(c\) gives \(G_{p''c'}=G_{p'c''}\) and \(G_{p'c}=G_{p''c}\). So \(u_{p'}-u_{p''}=Y_n(q_{c'}-q_{c''})\).
- *Assembly.* \(D=-\sum_cq_cg_c=-2^n(4/R_n^2)\bigl[\sigma_nQ+Y_n\sum_cx_c(x_{c'}-x_{c''})(P_{c'}-P_{c''})\bigr]\), and the last sum is \(-3Q\) (sympy).
- *Self term.* \(g_c[Ge_c]=2^n\cdot2(G_{cc}-G_{p'c})(P_c-\bar P)\). Applying the current functional to \(Ge_k\) gives \(2(G_{cc}-G_{p'c})\) for \(k=c\) and \(-E_n\) for \(k\neq c\). Hence \(\sigma_n=2(G_{cc}-G_{p'c})+E_n\), which gives \(\kappa^{\mathrm{self}}_n\). The cross part is the difference. \(\square\)

**Proof of BB.2.** \(\lambda^{-\alpha}=\frac{\sin\pi\alpha}{\pi}\int_0^\infty s^{-\alpha}(s+\lambda)^{-1}ds\) for \(\lambda>0\). The coefficients of \(Y_n\) and \(E_n\) sum to \(0\), so the kernel of \(L_n\) drops out. By AQ (1), eliminating the midpoints of \(G_{n+1}\) leaves \(c\,Z^{-1}-eP_C\) on \(V_n\), with \(Z=(L_n+p(s))^{-1}\). Woodbury on the corner modes multiplies a zero-sum corner column by \(1/(c-e\,r_n(ps))\) and a \(\mathbf 1_C\) column by \(1/(c-e\,u_n(ps))\), where \(u_n=g_0+2g_1\) is the \(\mathbf 1_C\) eigenvalue of the corner block. Midpoint values are \(M^{-1}Bw_{\mathrm{old}}\) with \(M^{-1}=(5+s)^{-1}(I+J/(2+s))\). The new neighbour \(p'\) is the midpoint of \(c\) and its level-\(n\) neighbour \(b'\) towards \(c'\). For the column \(e_{c'}-e_{c''}\), reflection gives \(w_c=0\) and \(w_{e'}=-w_{b'}\), so \(w_{p'}=w_{b'}/(5+s)\) (sympy). Therefore \(y_{n+1}(s)=y_n(ps)/((5+s)(c-er_n(ps)))\). Since \(r_{n+1}(s)=r_n(ps)/(c-er_n(ps))\), \(y/r\) gains a factor \(1/(5+s)\) per level. At level \(0\), \(b'=c'\) and \(y_0=r_0\). \(\square\)

**Proof of BB.3.**
- *Corner row.* The corner row of \((L_n+s)w=e_{c'}+e_{c''}\) reads \((2+s)w_c-w_{p'}-w_{p''}=0\) with \(w_{p'}=w_{p''}\). So \(e_n=w_{p'}-w_c=s\,g_n\).
- *Recursion.* \(3g=u-r\), and both modes take their Woodbury factor. This gives \(g_{n+1}=c\,g_n(ps)/((c-eu')(c-er'))\), and the stated \(F_n\) follows by algebra (sympy).
- *Bounds on \(t\) and \(pr\).* Entrywise positivity gives \(0<(L+s)^{-1}\mathbf 1_C\le(L+s)^{-1}\mathbf 1=\mathbf 1/s\), so \(t_n\le1\). AQ (3) gives \(p\,r_n(p)\le1\).
- *Ceiling.* Both denominators are therefore \(\ge4+s\), and \(F_n\le(2+s)(5+s)(6+s)/(4+s)^2\le\tfrac3{20}(5+s)^2\), since \(3(5+s)(4+s)^2-20(6+s)(2+s)=3s^3+19s^2+8s\).
- *Iteration.* \(\prod_{j<n}(5+p^js)=p^n(s)/s\) gives \(g_n(s)\le(3/20)^np^n(s)/(s^2(3+p^n(s)))\le(3/20)^n/s^2\). Also \(e_n=s\,g_n\le t_n/3\le1/3\).
- *Integral.* Split \(\int s^{-\alpha}e_n\) at \(s_1=(3/20)^n\): use \(e_n\le\tfrac13\) below and \(e_n\le(3/20)^n/s\) above. \(\square\)

**Proof of BB.4.** By BB.1, \(2^{-n}\kappa_n=-4\sigma_n/R_n^2+12Y_n/R_n^2\). By BB.2 and AQ (3) (\(R_n\ge R_1>0\)), \(12Y_n/R_n^2\le6\cdot5^{-n}/R_1\to0\). Next, \(\sigma_n\downarrow\sigma_\infty\ge2\cdot4^{-\alpha}\) and \(R_n\uparrow R_\infty\), which is finite for \(\alpha<\gamma\) and infinite for \(\alpha\ge\gamma\) (AQ (3), AR (3)–(4), AS (1)). At \(0.45\), AS's outward-rounded sandwich (re-run by the witness) encloses \(R_\infty\) and \(\sigma_\infty\), and interval division encloses \(K\). The \(\gamma\) and \(\alpha>\gamma\) asymptotics substitute AR (4b), (4c) and AZ's \(R_\alpha(n)=\frac{c_R}{\rho-1}\rho^n+C_\alpha+o(1)\) into BB.1. \(\square\)

**Proof of BB.5.** \(|\kappa^{\mathrm{cross}}_n|\le2^n(4/R_n^2)(3Y_n+E_n)\). Then apply BB.2, BB.3, and \(R_1\le R_n\le R_\infty\). The base \(2(3/20)^{1-\alpha}\) is \(<1\) exactly when \(\alpha<\alpha_*\). \(\square\)

**DERIVED (no centred corner stencil).** Every vertex of \(V_n\) lies in the closed \(60^\circ\) wedge at a corner \(c\). So a gradient stencil at \(c\) cannot be odd in the sense of BA.5. Any rule supported on the corner cell that is exact on affine functions is the one-sided cell gradient. The only gasket analogue of the torus symmetrisation is to drop the self term \(\varphi_{cc}\). BB.5 shows that dropping it removes the growth and the limit. Kept Failure BB.1 shows that it does not remove the residual at any finite level.

**KEPT FAILURE BB.1 (self-force removal is not an exact conservation law on the gasket).** Claim: with the self pair force removed, the integrated divergence after the Maxwell piece is exactly \(0\) on the gasket, as the centred stencil gives on the torus (BA.5). This is false. The certified cross constants at \(0.45\) (mpmath.iv, the BA pipeline) are
\[
\kappa^{\mathrm{cross}}_{1..4}=0.48882868502562873746\ldots,\;0.086933485121556483865\ldots,\;0.017575791192880243766\ldots,\;0.0039810486280556649983\ldots,
\]
all \(>0\). The recursion gives \(1.2757\cdot10^{-6}\) at \(n=10\) and \(3.0737\cdot10^{-12}\) at \(n=20\). It is nonzero at every computed level and tends to \(0\). A second failure: the decreasing modulus of BA's float table was not evidence of a zero limit. \(2^{-n}\kappa_n\) increases towards \(K\) from below: \(-0.880867405864\) at \(n=10\), \(-0.855253882504\) at \(n=20\), \(-0.854647244939\) at \(n=40\). BA's Aitken estimate "\(\approx-0.86\)" was \(0.005\) too low.

**CONJECTURE BB.2.** \(\kappa^{\mathrm{cross}}_n(\alpha)\sim c\,(2\cdot5^{\alpha-1}/3)^n\). On \(n=10\)–\(20\) the per-level ratio is \(0.2743\), against \(2\cdot5^{-0.55}/3=0.2750\). Both \(Y_n\) and \(E_n\) decay like \((5^{\alpha-1}/3)^n\) in floats. Only the cruder BB.3 rate is proved.

**NOT THIS.** Theorem BB is a statement about one discrete constant of one stated reading (Assumption BA). It does not choose the reading. \(K(0.45)\) is the limit of a renormalised (\(\times2^{-n}\)) discrete self-force constant. The unrenormalised \(\kappa_n\) still diverges like \(2^n\), and no continuum \(D\) is constructed. Per BB.5 the surviving piece is self-force, which is the interpretive conclusion of Kept Failure BA.1, now with a proof. A nonzero \(K\) is not a residual force on the device. Not thrust.

**Witness, computed by `scripts/mesh_maxwell_limit.py`** (about 2–3 min).
- *Sympy:* the Woodbury mode factors and \(t\)-recursion; \(F_n\) and its ceiling polynomial; the midpoint value \(w_{b'}/(5+s)\); the level-\(0\) data; the \(-3Q\) gradient identity; and \(\alpha=1\) consistency (BB.1 with \(R=(2/3)(5/3)^n\), \(Y=5^{-n}R/2\), \(\sigma=1\) reproduces BA.4).
- *Exact rationals:* \(r_n\), \(y_n\), \(g_n\), \(t_n\), and \(e_n=s\,g_n\) from the recursions equal exact resolvent solves at levels \(1\)–\(3\), \(s\in\{1/2,3,7/5\}\).
- *Dense floats:* BB.1 matches BA's \(\kappa_n\) and a direct self-force sum at levels \(1\)–\(6\) (\(\alpha=0.45\)) and \(1\)–\(4\) (\(\alpha=0.25,0.6\)). It also checks \(q_{c}=2/R\) and the bounds of BB.2–BB.3.
- *Recursion quadrature:* agrees with the dense mesh at \(n=1,4,6\), and continues to \(n=10,20,40\).
- *Intervals:* \(\kappa^{\mathrm{cross}}_{1..4}>0\), and AS's sandwich re-run gives \(K(0.45)\).
- *\(\alpha=0.75\):* \(2^{-n}\kappa_n=-0.0507,-0.0036,-0.00036\) at \(n=10,20,30\), tending to \(0\).

The script exits nonzero if any check fails.

## Theorem BC. The centred stencil on the gasket: action–reaction is centrality, no vertex admits an odd stencil, and the centred surrogate leaves \(\kappa_n\) unchanged

Kept Failure BA.1 traced the mesh residual to a stencil that is not odd: on the translation-invariant torus the centred stencil gives \(D=0\). This theorem asks the same question on the gasket itself. The centred rule is stated as an assumption. The result explains what replaces the torus mechanism.

**ASSUMPTION BC (centred reading, not a change to the freeze).** Assumption BA (1), (2), (4), (5), with item (3) replaced by the Green–Gauss vertex gradient
\[
g^{\mathrm{GG}}_x[u]=\frac1{|\mathcal C(x)|}\sum_{K\in\mathcal C(x)}\nabla\bigl(u|_K\bigr)_{\mathrm{affine}},
\]
the mean of the affine-interpolant gradients over the set \(\mathcal C(x)\) of level-\(n\) cells containing \(x\). On the triangular lattice this rule is exactly the centred six-neighbour stencil \(\tfrac13\sum_e e\,(u(x+e)-u(x))\) of Kept Failure BA.1's torus witness (sympy identity). So it is the centred stencil, written in a form that exists on any triangulated mesh. Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BC.**

1. **BC.1 (reflection law; every \(n\), every \(\alpha\), every \(S_3\)-equivariant linear gradient rule that kills constants).** For corners \(c\neq k\), let \(R_{ck}\) be the reflection that swaps them. Then \(H_{kc}=R_{ck}H_{ck}\) and \(H_{cc}=s\,(P_c-\bar P)\) for one scalar \(s\). Hence the Maxwell piece \(\varphi^{M}_{ck}\) is exactly the chord-parallel (central) component of the pair force \(\varphi_{ck}\). The residual is the self force plus the chord-normal (non-central) component. With \(h_\perp:=H_{cc'}\cdot\nu\), where \(\nu\) is the unit normal to the chord \(cc'\) pointing towards \(c''\), and with \(q_c=2x_c/R_n\) (BB.1),
   \[
   \kappa_n=-\frac{4}{R_n^2}\bigl(s+2\sqrt3\,h_\perp\bigr),\qquad\kappa^{\mathrm{self}}_n=-\frac{4s}{R_n^2},\qquad\kappa^{\mathrm{cross}}_n=-\frac{8\sqrt3\,h_\perp}{R_n^2}.
   \]
   The pair force respects action–reaction (\(\varphi=\varphi^{M}\)) if and only if \(s=0\) and \(h_\perp=0\), that is, no self force and central cross forces.
2. **BC.2 (no centred stencil anywhere on the gasket).** Every level-\(n\) cell is an upward triangle. A corner lies in one cell and every other vertex in exactly two. At every vertex of \(V_n\), the neighbour offsets lie in a closed half-plane, do not sum to \(0\), and are not closed under \(e\mapsto-e\). So no vertex, corners or junctions, carries an odd stencil in the sense of BA.5. The isometry group of \(V_n\) is \(D_3\), which contains no half-turn.
3. **BC.3 (the centred surrogate changes nothing where the charge is).** At each corner \(|\mathcal C(c)|=1\), so \(g^{\mathrm{GG}}_c\) is BA's one-sided cell gradient. All of the charge sits on the corners (Theorem H). Therefore \(D^{\mathrm{GG}}=D\) and \(\kappa^{\mathrm{GG}}_n=\kappa_n\) for every \(n\) and \(\alpha\). In particular \(\kappa^{\mathrm{GG}}_{1..4}(0.45)=-4.0427396046867839149\ldots\), \(-6.6139165202761608886\ldots\), \(-10.738239541819965987\ldots\), \(-18.526339511158477727\ldots\) (re-certified through the Green–Gauss implementation). By Theorem BB, \(\kappa^{\mathrm{GG}}_n\) grows like \(2^n\), with \(2^{-n}\kappa^{\mathrm{GG}}_n\to K(0.45)\in[-0.854646905120924,-0.854646905120848]\).
4. **BC.4 (the corner-cell pair kernel, explicitly).** \(H_{cc'}=2^n\bigl(-E_n(P_c-\bar P)+Y_n(P_{c'}-P_{c''})\bigr)\), so \(h_\parallel=2^{n-1}(E_n+Y_n)\), \(h_\perp=\tfrac{\sqrt3}{6}2^n(E_n-3Y_n)\), and \(s=2^n(\sigma_n-E_n)\). This reproduces BB.1. The non-centrality angle obeys \(\tan\theta_n=h_\perp/h_\parallel=(E_n-3Y_n)/(\sqrt3(E_n+Y_n))\). Certified at \(0.45\) (mpmath.iv, \(120\) digits): \(\tan\theta_{1..4}=-0.3917585210483279\ldots\), \(-0.277452251033554\ldots\), \(-0.2247945906847502\ldots\), \(-0.2038773766517652\ldots\), with \(E_n/Y_n=1.383015302328912\ldots\), \(1.701677898838603\ldots\), \(1.879032444489115\ldots\), \(1.956117962609032\ldots\), and \(s=2.62859242109\ldots\), \(5.37968631401\ldots\), \(10.8101556491\ldots\), \(21.6368071437\ldots>0\).
5. **BC.5 (exact small-\(s\) ratio).** In BB's resolvent notation, \(e_n(s)/y_n(s)=\prod_{k<n}(6+z_k)/(6+z_k-2t_{n-1-k}(z_{k+1}))\) with \(z_k=p^k(s)\). As \(s\to0\), \(t_m(0^+)=2/(3^m+1)\), and the product telescopes to \(\lim_{s\to0}e_n/y_n=2\cdot3^n/(3^n+1)\).

**Proof of BC.1.** By BA.2's argument, \(g_{\sigma c}[w\circ T_\sigma^{-1}]=R_\sigma g_c[w]\) and \(G\) is \(S_3\)-invariant. So \(H_{\sigma c,\sigma k}=R_\sigma H_{ck}\). Take \(\sigma=(c\,k)\): this gives \(H_{kc}=R_{ck}H_{ck}\). The stabiliser of \(c\) is the reflection fixing \(c\), so \(H_{cc}\) lies on its axis \(P_c-\bar P\). Then \(\tfrac12(\varphi_{ck}-\varphi_{kc})=-q_cq_k\tfrac12(I-R_{ck})H_{ck}\), and \(\tfrac12(I-R_{ck})\) is the projection onto the chord. The sums \(\sum_{c\neq k}x_cx_kH_{ck}=2\sqrt3\,h_\perp Q(x)\) and \(\sum_cx_c^2H_{cc}=sQ(x)\) are sympy identities for a generic \(H_{cc'}\) propagated by \(S_3\). \(\square\)

**Proof of BC.2.** The cells are the images of the unit triangle under the maps \(x\mapsto x/2+b\), so they are all upward translates. A vertex is a corner of at most two cells, in two different roles (bottom-left, bottom-right, top). The offset pairs of the roles are \(\{(h,0),(h/2,h\sqrt3/2)\}\), \(\{(-h,0),(-h/2,h\sqrt3/2)\}\), and \(\{(\pm h/2,-h\sqrt3/2)\}\). No union of two of them is closed under negation, and each union lies in a closed half-plane. The script also checks this at every vertex for \(n=1,\ldots,6\). An isometry of \(V_n\) preserves the edges (pairs at distance \(h\)), hence degrees. The corners are the only degree-\(2\) vertices, so it permutes them and lies in \(D_3\). \(\square\)

**Proof of BC.3–BC.5.** BC.3 is BC.2 at the corners plus \(\operatorname{supp}q=C\). BC.4: substitute BB.1's gradient split into \(g_c[Ge_{c'}]\). The current functional gives \(-E_n\) (BB, self term) and the tangential increment gives \(Y_n\). BC.5: BB's recursions give \(e_{n+1}(s)/y_{n+1}(s)=\frac{6+s}{6+s-2t_n(ps)}\,e_n(ps)/y_n(ps)\), with \(e_0/y_0=1\). At \(s=0\), \(t_{m+1}=t_m/(3-t_m)\), so \(\tau_m:=1/t_m\) satisfies \(\tau_{m+1}=3\tau_m-1\), \(\tau_m=(3^m+1)/2\). Each factor is \(3/(3-t_m)=3\tau_m/\tau_{m+1}\), so the product telescopes to \(3^n\tau_0/\tau_n=2\cdot3^n/(3^n+1)\). Exact rationals at \(n\le4\) and sympy limits at \(n\le3\) confirm it. \(\square\)

**DERIVED (what replaces the torus mechanism).** On the torus, the half-turn about the midpoint of \(P_c,P_k\) is a mesh symmetry that swaps \(c\) and \(k\) with linear part \(-I\). A rule equivariant under it, which is exactly an odd stencil, gets \(H_{kc}=-H_{ck}\). The half-turn about \(c\) gives \(H_{cc}=0\). That is BA.5. On the gasket the only isometry swapping two corners is a reflection, so symmetry makes only the chord-parallel part antisymmetric. Newton's third law becomes a condition (\(s=0\), \(h_\perp=0\)), not a consequence. The non-translation-invariant structure that produces the residual is therefore the \(D_3\) symmetry group, which has no inversion, together with the degree-\(2\) wedge corners where the charge sits. The residual decimates like \(2^n\) with the finite renormalised limit \(K\) (BB.4). The non-central cross part tends to \(0\) (BB.5).

**KEPT FAILURE BC.1 (a centred stencil does not restore Newton's third law on the gasket).** Claim: the centred stencil that gives \(D=0\) on the torus gives \(D=0\) on the frozen \(0.45\) gasket after the Maxwell piece is removed. This is false. The only form of that stencil that exists on the gasket (Green–Gauss) coincides with the one-sided rule at every charged vertex. **Witness:** certified \(\kappa^{\mathrm{GG}}_{1..4}(0.45)=\kappa_{1..4}\neq0\) (BC.3). The cross pair force is certified non-central, \(\tan\theta_{1..4}<0\), and the self-force coefficient is certified \(s>0\) at levels \(1\)–\(4\) (BC.4). The torus cancellation used a half-turn symmetry that the gasket does not have (BC.2).

**CONJECTURE BC.2.** For every \(\alpha\in(0,1)\), \(E_n/Y_n\to2\) from below, with \(2-E_n/Y_n\sim c(\alpha)3^{-n}\). Hence \(\tan\theta_n\to-1/(3\sqrt3)=-0.19245008973\ldots\), a non-centrality angle of \(10.893395^\circ\) that does not depend on \(\alpha\). Evidence (recursion quadrature, floats, not certified): at \(0.45\), \(E_n/Y_n=1.99993691129\) at \(n=10\), \(1.99999999893\) at \(n=20\), and \(1.99999999999998\) at \(n=30\), with \((2-E_n/Y_n)3^n\to3.7255634\). At \(\alpha=0.25,0.6,0.8\) and \(n=16\), \((2-E/Y)3^n\approx4.4720,3.2076,2.5724\). The pointwise ratio \(e_n(5^{-n}t)/y_n(5^{-n}t)\) also tends to \(2\) for \(t\le10^3\). Numerically \(e_n/y_n\le2\cdot3^n/(3^n+1)\) for every \(s\) (\(n\le11\)). The factorwise bound that would prove this fails, so it remains a conjecture. A proof of \(E_n<3Y_n\) for every \(n\) would extend Kept Failure BB.1 (\(\kappa^{\mathrm{cross}}_n>0\)) from levels \(1\)–\(4\) to all levels.

**NOT THIS.** Theorem BC does not choose the reading. A stencil designed to satisfy \(s=-2\sqrt3h_\perp\), which is \(D=0\) without action–reaction, or \(s=h_\perp=0\), would be a fit to the desired answer. It is not attempted. Doubling the mesh through a corner (the bow-tie, where the corner becomes an interior vertex with an odd neighbour set) would change the frozen mesh, so it is not used. A non-central residual is not a force on the device. Not thrust.

**Witness, computed by `scripts/mesh_centred_stencil.py`** (about 1.5 min).
- *Sympy:* the reflection law and the chord-parallel form of the Maxwell piece for a generic equivariant kernel; the \(\kappa\) formula; the cell-rule components; Green–Gauss equals centred on the triangular lattice; the telescoping product; \(\lim_{s\to0}e_n/y_n\) at \(n\le3\).
- *Exact:* the product formula at \(n\le4\); the cell and vertex combinatorics and the half-plane property at \(n\le6\); no half-turn symmetry.
- *Intervals (BA pipeline with the Green–Gauss rule, levels \(1\)–\(4\)):* \(\kappa^{\mathrm{GG}}_n=\) BA.3; the Maxwell sum \(\ni0\); \(H_{c'c}-R H_{cc'}\ni0\); the BC.4 components; the \(\kappa\) formula; \(\tan\theta_n<0\); \(s>0\).
- *Recursion quadrature:* the \(E_n/Y_n\) table to \(n=30\).

The script exits nonzero if any check fails.

## Theorem BD. The non-centrality ratio tends to 2: Conjecture BC.2 holds for every \(\alpha\), with an \(\alpha\)-dependent \(3^{-n}\) constant

Theorem BC left open whether the corner-cell non-centrality \(\tan\theta_n=(E_n-3Y_n)/(\sqrt3(E_n+Y_n))\) has an \(\alpha\)-independent limit. This theorem proves it does, from one exact identity between the two resolvent entries.

**SETTING.** Assumption BA/BC, nothing new. BB's notation: for \(s>0\), \(y_n(s)=[(L_n+s)^{-1}]_{p'c'}-[(L_n+s)^{-1}]_{p'c''}\), \(e_n(s)=s\,g_n(s)\), \(t_n(s)=s[(L_n+s)^{-1}\mathbf 1_C]_c\), \(p(s)=s(s+5)\), and \(Y_n,E_n=\frac{\sin\pi\alpha}{\pi}\int_0^\infty s^{-\alpha}(y_n,e_n)\,ds\). Write \(d_n:=2-E_n/Y_n\). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BD.**

1. **BD.1 (exact identity; every \(n\), every \(s>0\)).** \(e_n(s)=(2-t_n(s))\,y_n(s)\). Hence \(d_n=\langle t_n\rangle_n\), the mean of \(t_n\) under the probability measure proportional to \(s^{-\alpha}y_n(s)\,ds\).
2. **BD.2 (the corner capture \(t_n\)).** \(t_n(s)=\tfrac13\sum_\lambda\|P_\lambda\mathbf 1_C\|^2\,s/(\lambda+s)\). So \(t_n\) is increasing and concave, with \(0<t_n<1\) for \(n\ge1\), \(t_n(0^+)=a_n=2/(3^n+1)\), and \(t_n'(0)=b_n=3^n(5^n-1)/(6(3^n+1)^2)\le\tfrac16(5/3)^n\). Consequences for every \(n\ge1\) and \(\alpha\in(0,1)\):
   - \(1<E_n/Y_n<2\cdot3^n/(3^n+1)\). Pointwise, \(e_n(s)/y_n(s)\le2\cdot3^n/(3^n+1)\) for every \(s\); BC observed this numerically for \(n\le11\).
   - \(3Y_n-E_n=Y_n(1+d_n)>0\), so \(\kappa^{\mathrm{cross}}_n=2^n(4/R_n^2)Y_n(1+d_n)>0\) **at every level and every \(\alpha\)**. Kept Failure BB.1, certified at levels \(1\)–\(4\), holds for all \(n\).
   - \(\tan\theta_n=-\frac1{3\sqrt3}-\frac{4d_n}{3\sqrt3(3-d_n)}<-\frac1{3\sqrt3}\): the limit is approached from below.
3. **BD.3 (Conjecture BC.2; two-sided rate; every \(\alpha\in(0,1)\), \(n\ge2\)).**
   \[
   \frac2{3^n+1}<2-\frac{E_n}{Y_n}\le\frac2{3^n+1}+B_\alpha3^{-n},\qquad B_\alpha=\frac{(1-\alpha)J_\alpha e^{0.1}(3+e^{0.1})}{18},
   \]
   with \(J_\alpha=\int_0^\infty\frac{\tau^{1-\alpha}d\tau}{(1+\tau/25)(1+\tau/125)}=\frac{125}{4}\frac{\pi}{\sin\pi(2-\alpha)}\bigl(25^{1-\alpha}-125^{1-\alpha}\bigr)\). At \(0.45\), \(B=115.196453544\ldots\). So \(E_n/Y_n\to2\) and \(\tan\theta_n\to-1/(3\sqrt3)=-0.19245008973\ldots\) for every \(\alpha\). Conjecture BC.2 is proved.
4. **BD.4 (the \(3^{-n}\) constant exists; every \(\alpha\in(0,1)\)).** Put \(T_n(\tau)=3^nt_n(5^{-n}\tau)\) and \(\mathcal Y_n(\tau)=3^ny_n(5^{-n}\tau)\). Then \(T_n\uparrow T\le2+\tau/6\), and \(\mathcal Y_n\to\mathcal Y>0\) pointwise, so
   \[
   c(\alpha):=\lim_{n\to\infty}3^n\Bigl(2-\frac{E_n}{Y_n}\Bigr)=\frac{\int_0^\infty\tau^{-\alpha}T(\tau)\mathcal Y(\tau)\,d\tau}{\int_0^\infty\tau^{-\alpha}\mathcal Y(\tau)\,d\tau},\qquad2<c(\alpha)\le2+B_\alpha ,
   \]
   and \(3^n(\tan\theta_n+\tfrac1{3\sqrt3})\to-4c(\alpha)/(9\sqrt3)\). Certified (mpmath.iv): \(c(0.45)\in[3.72198897,3.72914181]\), which contains BC's float \(3.7255634\), and \(c(0.8)\in[2.57126273,2.57352573]\).

**Proof of BD.1.** BC.5 gives \(Q_{n+1}(s)=\frac{6+s}{6+s-2t}Q_n(ps)\) for \(Q_n=e_n/y_n\), and BB gives \(t_{n+1}(s)=\frac{(2+s)t}{6+s-2t}\), with \(t=t_n(ps)\). If \(Q_n=2-t_n\), then \(Q_{n+1}(s)=\frac{(6+s)(2-t)}{6+s-2t}=2-\frac{(2+s)t}{6+s-2t}=2-t_{n+1}(s)\) (sympy). The base case is \(Q_0=1=2-t_0\). Exact rationals at \(n\le8\), and exact resolvent solves on the meshes at levels \(1\)–\(3\), confirm the identity independently of the recursions. Integrate against \(s^{-\alpha}\). \(\square\)

**Proof of BD.2.** \(\mathbf 1_C\) and \(L_n\) are \(S_3\)-invariant, and \(S_3\) is transitive on \(C\). So \(P_\lambda\mathbf 1_C\) takes one value on \(C\), namely \(\tfrac13\langle\mathbf 1_C,P_\lambda\mathbf 1_C\rangle\ge0\). This gives the Stieltjes form, monotonicity and concavity. The \(\lambda=0\) term is \(3/N_n=2/(3^n+1)\), and \(s/(\lambda+s)<1\) for \(\lambda>0\) gives \(t_n<1\). Differentiating BB's \(t\)-recursion at \(0\) gives \(b_{n+1}=(60b_n+4a_n-2a_n^2)/(6-2a_n)^2\) with \(b_0=0\), and the closed form satisfies it (sympy, and exact to \(n=29\)). It equals \(\tfrac13\langle w,L^+w\rangle\), \(w=\mathbf 1_C-3/N\), on the dense meshes. Then \(t_n\ge a_n\) gives \(e_n/y_n\le2-a_n\), and \(t_n<1\) gives \(E_n>Y_n\). The cross part is BB.1 with \(3Y-E=Y(1+d)\), and the \(\tan\theta\) form is algebra (sympy). \(\square\)

**Proof of BD.3.** *Lower bound:* \(t_n>a_n\) on \(s>0\). *Upper bound:* concavity gives \(\langle t_n\rangle\le a_n+b_n\langle s\rangle\), with \(\langle s\rangle=I_1/I_0\) and \(I_k=\int s^{k-\alpha}y_n\).
- *Upper bound on \(I_1\).* \(y_n=r_n(s)\,s/p^n(s)\) (BB.2), \(r_n(s)\le r_n(0)=\tfrac13(5/3)^n\), and \(s/p^n(s)=5^{-n}\prod_{j<n}(1+p^j(s)/5)^{-1}\le5^{-n}/((1+5^{n-2}s)(1+5^{n-3}s))\) for \(n\ge2\), because \(p^j(s)\ge5^js\). Substituting \(s=5^{-n}\tau\) gives \(I_1\le\tfrac13 3^{-n}5^{-n(2-\alpha)}J_\alpha\).
- *Lower bound on \(I_0\).* \(\Phi_s(r)\ge\frac{(2+s)(5+s)}{6+s}r\ge\tfrac53r\) gives \(r_n(s)\ge(5/3)^n/(3+p^n(s))\). For \(s\le5^{-n}\), induction gives \(p^n(s)\le e^{0.1}5^ns\). So \(y_n\ge3^{-n}/(e^{0.1}(3+e^{0.1}))\) there, and \(I_0\ge3^{-n}5^{-n(1-\alpha)}/((1-\alpha)e^{0.1}(3+e^{0.1}))\).
- Combining these with \(b_n\le\tfrac16(5/3)^n\) gives \(3^nb_n\langle s\rangle\le B_\alpha\). \(J_\alpha\) is the partial-fraction continuation of \(\int t^{\mu-1}/(c+t)=c^{\mu-1}\pi/\sin\pi\mu\) (sympy for the fractions; quadrature agrees). \(\square\)

**Proof of BD.4.**
- *\(T\).* From the \(t\)-recursion, \(T_{n+1}(\tau)=\frac{3(2+s)}{6+s-2t}T_n(\tau(1+\delta))\), with \(s=5^{-n-1}\tau\), \(\delta=s/5\), \(5^np(s)=\tau(1+\delta)\). Since \(3(2+s)\ge6+s-2t\) and \(T_n\) is increasing, \(T_{n+1}\ge T_n\). Concavity gives \(T_n\le3^na_n+3^nb_n5^{-n}\tau\le2+\tau/6\). With \(T_n(\tau(1+\delta))\le(1+\delta)T_n(\tau)\), the increment satisfies \(\log(T_{n+1}/T_n)\le\delta+s/2+3^{-n}(1+\tau(1+\delta)/12)\), which is summable.
- *\(\mathcal Y\).* \(y_{n+1}(s)=\frac{(2+s)y_n(ps)}{6+s-2p(s)r_n(ps)}\) (sympy). So \(\mathcal Y_{n+1}(\tau)=F_n\,\mathcal Y_n(\tau(1+\delta))\) with \(0\le\log F_n\le s/2+\tau3^{-n}(1+\delta)/6\), using \(p\,r_n(p)\le\min(1,p\,r_n(0))\).
- *Log-Lipschitz step.* \(s\,r_n(s)\) is increasing. \(p^j\) has nonnegative coefficients and degree \(2^j\), and \(p^j(x)\le5^{j-n}\psi_n(\tau)\) with \(\psi_n(\tau)=p^n(5^{-n}\tau)\uparrow\), \(\psi_n\le\bar\Psi(\tau)=p^k(ue^{u/10})\), \(u=5^{-k}\tau\le1\). These give \(0\le\log\mathcal Y_n(\tau)-\log\mathcal Y_n(\tau(1+\delta))\le\delta+\bar\Psi(\tau)\eta_ne^{\eta_n}/45\), \(\eta_n=2^n\delta\).
- *Limits.* Both increments are summable, so \(\mathcal Y_n\to\mathcal Y\in(0,\infty)\), with explicit tail factors \(e^{\pm\varepsilon^{\mathcal Y}_N(\tau)}\), \(e^{\varepsilon^T_N(\tau)}\).
- *Dominated convergence.* \(\mathcal Y_n\le\tfrac13/((1+\tau/25)(1+\tau/125))\) dominates, and \(\mathcal Y_n\ge e^{-0.1}/(3+e^{0.1})\) on \(\tau\le1\) keeps the denominator away from \(0\).
- *\(c>2\).* \(T\ge2\), and \(T\ge T_1=3(10+\tau)/(20+\tau)>2\) for \(\tau>10\).
- *Certificate.* \(T_{60},\mathcal Y_{60}\) are computed in interval arithmetic on a geometric grid of \(13485\) cells over \([10^{-20},3125]\). On each cell, \(T\) increasing, \(\mathcal Y\) decreasing and \(\varepsilon\) increasing give outward bounds. The head \([0,10^{-20}]\) uses \(0\le T\le2+\tau/6\), \(\mathcal Y\le\tfrac13\). The tail \([3125,\infty)\) uses \(\mathcal Y\le\tau/(3\psi_6(\tau))\) and \(\psi_6(\tau)\ge(\tau/3125)^8P_8\) (nonnegative coefficients, exact sympy polynomial). \(\square\)

**KEPT FAILURE BD.1 (the \(3^{-n}\) constant is not the \(s\to0\) value and is not universal).** Two guesses fail. (i) The correction is the small-\(s\) one, \(2-E_n/Y_n\sim t_n(0^+)=2\cdot3^{-n}\), i.e. \(c=2\). (ii) Like the limit \(-1/(3\sqrt3)\), the constant does not depend on \(\alpha\). Both are false. **Witness:** certified \(c(0.45)\in[3.72198897,3.72914181]\) and \(c(0.8)\in[2.57126273,2.57352573]\), disjoint and both \(>2\). The proved bound \(c(\alpha)>2\) holds for every \(\alpha\). The measure \(s^{-\alpha}y_n\,ds\) lives on the scale \(s\sim5^{-n}\), where \(t_n\) has already grown from \(2\cdot3^{-n}\) to \(3^{-n}T(\tau)\), with \(T(5)\approx2.81\) and \(T(50)\approx8.79\). Floats: \(c(0.25)\approx4.4719725\), \(c(0.45)\approx3.7255634\), \(c(0.8)\approx2.5724294\).

**CONJECTURE BD.2.** \(c(\alpha)\) decreases on \((0,1)\), with \(c(1^-)=2\) (the weight \(\tau^{-\alpha}\) concentrates at \(\tau=0\), where \(T=2\)) and \(c(0^+)=\infty\) (the bound \(B_\alpha\) diverges like \(1/\alpha\)). Only the three float values support it. No closed form for \(T\), \(\mathcal Y\) or \(c\) is known. Theorem BE decides it in part: the decrease and \(c(1^-)=2\) are proved, and \(c(0^+)=\infty\) is false (Kept Failure BE.1).

**NOT THIS.** Theorem BD is the asymptotics of one ratio of two discrete Green-kernel entries on the combinatorial gasket. It fixes the non-centrality angle of the corner-cell pair force under one stated reading (Assumption BA/BC) at \(\arctan(1/(3\sqrt3))=10.893395^\circ\) in the limit. It does not choose the reading. A positive cross part at every level means the self-force-removed residual is never \(0\); it does not mean a force on a device. The constant \(c(\alpha)\) is a refinement rate, not a physical coefficient. Not thrust.

**Witness, computed by `scripts/mesh_noncentrality_limit.py`** (about 4–5 min).
- *Sympy:* the induction step of BD.1; the \(y\)-recursion; the \(t'(0)\) recursion and the closed form of \(b_n\); the \(\tan\theta\) form; the \(J_\alpha\) partial fractions.
- *Exact:* BD.1 for \(n\le8\) at four \(s\); BD.1 against exact resolvent solves at levels \(1\)–\(3\); \(a_n,b_n\) for \(n\le29\).
- *Dense meshes:* nonnegative Stieltjes weights of \(t_n\) and \(b_n=\tfrac13\langle w,L^+w\rangle\) at levels \(1\)–\(4\); \(2-E_n/Y_n=\langle t_n\rangle\) at levels \(1,2,4\).
- *Intervals:* BC's certified \(E_n/Y_n\) at levels \(1\)–\(3\) lie in the BD.2/BD.3 bracket, with \(\tan\theta_n<-1/(3\sqrt3)\); \(B_{0.45}\); \(T_n\) monotone at three \(\tau\); \(c(0.45)\) and \(c(0.8)\).
- *Recursion quadrature:* \(3^nd_n\) at \(n\le30\) for \(\alpha=0.25,0.45,0.8\).

The script exits nonzero if any check fails.

## Theorem BE. The \(3^{-n}\) constant \(c(\alpha)\) decreases strictly from a finite value to \(2\): Conjecture BD.2 decided in part

Theorem BD left open Conjecture BD.2: does \(c(\alpha)=\lim3^n(2-E_n/Y_n)\) decrease on \((0,1)\), with \(c(1^-)=2\) and \(c(0^+)=\infty\)? Monotonicity and the endpoint \(1\) are true. The endpoint \(0\) is false.

**SETTING.** Assumption BA/BC and the notation of Theorems BB–BD, nothing new: \(d_n(\alpha)=2-E_n/Y_n=\langle t_n\rangle_n\) under \(s^{-\alpha}y_n(s)\,ds\) (BD.1), and \(T,\mathcal Y\) are the scaling limits of BD.4, so \(c(\alpha)=\int\tau^{-\alpha}T\mathcal Y/\int\tau^{-\alpha}\mathcal Y\). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BE.**

1. **BE.1 (monotone in \(\alpha\) at every level).** For \(0<\alpha<\beta<1\) and every \(n\ge1\), \(d_n(\alpha)>d_n(\beta)\). So \(E_n/Y_n\) is strictly increasing in \(\alpha\), and so is \(\tan\theta_n\): the non-centrality \(|\theta_n|\) shrinks as \(\alpha\) grows.
2. **BE.2 (monotone limit).** \(c(\alpha)\) is strictly decreasing on \((0,1)\).
3. **BE.3 (\(c(1^-)=2\), with slope).** \(c(\alpha)\le2+B_\alpha\) and \(B_\alpha=O(1-\alpha)\), since \(J_1=\tfrac{125}4\ln5\) is finite. More precisely,
   \[
   \frac{c(\alpha)-2}{1-\alpha}\to3S,\qquad S=\int_0^\infty\frac{T(\tau)-2}{\tau}\,\mathcal Y(\tau)\,d\tau,\qquad 3S\in[2.70632481,\,2.70903118]\ \text{(certified)}.
   \]
4. **BE.4 (\(c(0^+)\) is finite).** \(c(0^+)=c_0:=\int_0^\infty T\mathcal Y\,d\tau\big/\int_0^\infty\mathcal Y\,d\tau\), certified \(c_0\in[5.48553067,\,5.50249013]\). By BE.2, \(c(\alpha)<c_0\) for every \(\alpha\in(0,1)\).

**Proof of BE.1 and BE.2.** Let \(m\) be a positive measure on \((0,\infty)\), \(f\) nondecreasing, and \(\langle f\rangle_\alpha=\int x^{-\alpha}f\,dm/Z_\alpha\). Then
\[
Z_\alpha Z_\beta\bigl(\langle f\rangle_\alpha-\langle f\rangle_\beta\bigr)=\tfrac12\iint\bigl(f(x)-f(y)\bigr)\bigl(x^{-\alpha}y^{-\beta}-x^{-\beta}y^{-\alpha}\bigr)\,m(dx)\,m(dy).
\]
The kernel equals \((xy)^{-\alpha}(y^{\alpha-\beta}-x^{\alpha-\beta})\), which has the sign of \(x-y\) when \(\alpha<\beta\) (sympy). So the integrand is \(\ge0\). It is \(>0\) on a set of positive \(m\otimes m\) measure unless \(f\) is \(m\)-a.e. constant. This is the monotone-likelihood-ratio ordering of the tilts \(x^{-\alpha}m\). For BE.1 take \(m=y_n\,ds\) (\(y_n>0\), BB.2) and \(f=t_n\). By BD.2 \(t_n\) is a Stieltjes function, strictly increasing from \(a_n<1\) to \(1\). For BE.2 take \(m=\mathcal Y\,d\tau\) (\(\mathcal Y>0\)) and \(f=T\). \(T\) is nondecreasing as a limit of the nondecreasing \(T_n\). It is not constant, because \(2\le T\le2+\tau/6\) gives \(T(0^+)=2\), while \(T\ge T_1=2+(\tau-10)/(20+\tau)>2\) for \(\tau>10\). The integrals converge by BD.4's dominations. The \(\tan\theta_n\) statement follows from BD.2's form \(\tan\theta_n=-\frac1{3\sqrt3}-\frac{4d_n}{3\sqrt3(3-d_n)}\), which is decreasing in \(d_n\in(0,1)\). \(\square\)

**Proof of BE.3.**
- *First bound.* \(B_\alpha=(1-\alpha)J_\alpha e^{0.1}(3+e^{0.1})/18\), and \(J_\alpha\to\tfrac{125}4\ln5\) as \(\alpha\to1\) (sympy limit). So \(c(\alpha)-2\le B_\alpha\to0\).
- *Slope.* Write \(c-2=N(\alpha)/D(\alpha)\), with \(N=\int\tau^{-\alpha}(T-2)\mathcal Y\) and \(D=\int\tau^{-\alpha}\mathcal Y\).
- *Numerator.* \(\tau^{-\alpha}\le\tau^{-1}\) on \((0,1]\) and \(\le1\) on \([1,\infty)\). Also \((T-2)/\tau\le1/6\) and \((T-2)\mathcal Y\le\tau\mathcal Y/6\), which is integrable (BD.4). So dominated convergence gives \(N\to S\).
- *Denominator.* \(\mathcal Y\le\tau/(3\psi(\tau))\le\tfrac13\), because \(\psi_n(\tau)\ge\tau\). And for \(\tau\le1\), \(\mathcal Y\ge e^{-\tau/10}/(3+\tau e^{\tau/10})\ge\tfrac13-C\tau\), by BD.3's lower bound on \(r_n\) and \(\psi_n(\tau)\le\tau e^{\tau/10}\). Hence \((1-\alpha)\int_0^1\tau^{-\alpha}\mathcal Y\to\tfrac13\). The part over \([1,\infty)\) stays bounded, so \((1-\alpha)D\to\tfrac13\).
- *Certificate.* \(T\) is concave, as a limit of the concave \(T_n\), with \(T(0^+)=2\). So \((T-2)/\tau\) is nonincreasing, and \(\mathcal Y\) is nonincreasing. On each cell the integrand therefore lies between its outward interval bounds at the two endpoints, computed from \(T_{60},\mathcal Y_{60}\) with BD.4's tail factors \(T_{60}\le T\le T_{60}e^{\varepsilon^T}\) and \(e^{\pm\varepsilon^{\mathcal Y}}\). This uses \(12903\) cells on \([10^{-8},3125]\). The head contributes \(\le10^{-8}/18\) (\((T-2)/\tau\le\tfrac16\), \(\mathcal Y\le\tfrac13\)). The tail uses \(\mathcal Y\le\tau/(3\psi_6(\tau))\) and BD's exact polynomial bound \(\psi_6(\tau)\ge(\tau/3125)^8P_8\). Result: \(S\in[0.9021082713,0.9030103931]\). \(\square\)

**Proof of BE.4.** For \(\alpha\le\tfrac12\), \(\tau^{-\alpha}\le1+\tau^{-1/2}\). \(T\mathcal Y\le(2+\tau/6)\mathcal Y\), and \(\mathcal Y\le\tau/(3\psi_6(\tau))=O(\tau^{-7})\) at infinity. So dominated convergence gives the numerator and the denominator of \(c(\alpha)\) as \(\alpha\to0^+\), and the limit denominator is \(>0\). The certificate is BD.4's enclosure run at \(\alpha=0\) (\(13485\) cells, same head and tail bounds). \(\square\)

**KEPT FAILURE BE.1 (\(c(0^+)=\infty\) is false).** Conjecture BD.2 guessed \(c(0^+)=\infty\) because the proved bound \(B_\alpha\) diverges like \(1/\alpha\) (sympy: \(\alpha J_\alpha\to3125\)). That divergence belongs to the bound, not to \(c\). BD.3 dominates \(\mathcal Y\) by the polynomial \((1+\tau/25)^{-1}(1+\tau/125)^{-1}\), whose first moment diverges at \(\alpha=0\). But the true \(\mathcal Y\le\tau/(3\psi(\tau))\) decays faster than any power. **Witness:** certified \(c(0^+)\in[5.48553067,5.50249013]\), finite and above the certified \(c(0.45)\) and a coarse certified \(c(0.25)\in[4.416733,4.537510]\). The monotone direction of BD.2 and \(c(1^-)=2\) survive (BE.2, BE.3).

**CONJECTURE BE.2.** \(c\) is convex on \((0,1)\). Only floats support it: the secant slopes of \(3^{20}d_{20}(\alpha)\) over \(\alpha=0.1,0.25,0.45,0.65,0.8,0.95\) are \(-4.0095,-3.7321,-3.4136,-3.1361,-2.8983\) (increasing), and the slope at \(1\) is \(-3S\approx-2.708\). No closed form for \(T\), \(\mathcal Y\), \(c\), \(S\) or \(c_0\) is known. Theorem BF proves the convexity (strict, with \(c''>0\) on \([0,1]\)); the closed forms remain open.

**NOT THIS.** Theorem BE is about how one refinement-rate constant of a ratio of two discrete Green-kernel entries on the combinatorial gasket depends on the fractional exponent. It does not choose the reading (Assumption BA/BC). It does not select \(\alpha\), and it does not change the limit angle \(\arctan(1/(3\sqrt3))\). Its monotonicity says only that the finite-level non-centrality is smaller at larger \(\alpha\). \(c\), \(S\) and \(c_0\) are rates, not physical coefficients. Not thrust.

**Witness, computed by `scripts/mesh_noncentrality_monotone.py`** (about 4–5 min).
- *Sympy:* the kernel factorization and the symmetrization; \(J_1=\tfrac{125}4\ln5\), \((1-\alpha)J_\alpha\to0\), \(\alpha J_\alpha\to3125\); \(T_1-2=(\tau-10)/(20+\tau)\).
- *Recursion quadrature (floats):* \(d_n(\alpha)\) strictly decreasing over six \(\alpha\) at levels \(1,2,4\); \(3^{20}d_{20}(\alpha)\) decreasing, with the secant slopes above.
- *Intervals:* \(S\) (hence \(3S\)); \(c(0^+)\); a coarse \(c(0.25)\), ordered between \(c(0.45)\) and \(c(0^+)\).

The script exits nonzero if any check fails.

## Theorem BF. The \(3^{-n}\) constant \(c(\alpha)\) is strictly convex: Conjecture BE.2 proved

Theorem BE left open Conjecture BE.2: is \(c(\alpha)=\int\tau^{-\alpha}T\mathcal Y/\int\tau^{-\alpha}\mathcal Y\) convex on \((0,1)\)? It is, strictly, and \(c''>0\) up to and including both endpoints.

**SETTING.** Assumption BA/BC and the notation of Theorems BB–BE, nothing new. \(T,\mathcal Y\) are BD.4's scaling limits. Put \(u=\log\tau\), \(\varepsilon=1-\alpha\), \(P=(T-2)\mathcal Y\ge0\), and pick any shift \(s_0\), with \(v=u-s_0\). Stage 1 is not edited. No \(W\), \(0.08\), \(0.23\), or \(30\,\mu\)N/kW enters.

**THEOREM BF.**

1. **BF.1 (regular representation).** For every \(s_0\) and \(\varepsilon\in(0,1)\),
   \[
   c(1-\varepsilon)-2=F(\varepsilon):=\frac{\varepsilon N(\varepsilon)}{\tfrac13+\varepsilon G(\varepsilon)},\qquad N(\varepsilon)=\int_{\mathbb R}e^{\varepsilon v}P\,dv,\quad G(\varepsilon)=\int_{\mathbb R}e^{\varepsilon v}\bigl(\mathcal Y-\tfrac13\mathbf 1_{v<0}\bigr)\,dv .
   \]
   \(N\) and \(G\) are smooth on a neighbourhood of \([0,1]\), so \(c\) extends smoothly to \([0,1]\), with \(c(1)=2\), \(c(0)=c_0\) (BE.4) and \(c''(\alpha)=F''(1-\alpha)\).
2. **BF.2 (strict convexity).** \(c''(\alpha)>0\) for every \(\alpha\in[0,1]\). The certified lower bound is \(c''\ge0.52046\) on all of \([0,1]\). So \(c\) is strictly convex on \((0,1)\), which proves Conjecture BE.2.
3. **BF.3 (two-sided linear bounds).** \(c'(1^-)=-3S\) (BE.3), and \(c'(0^+)=-F'(1)\in[-4.341285,\,-4.227515]\) (certified). Convexity then gives, for every \(\alpha\in(0,1)\),
   \[
   \max\bigl(2+3S(1-\alpha),\;c_0+c'(0^+)\,\alpha\bigr)<c(\alpha)<2+(c_0-2)(1-\alpha).
   \]
   With BE's enclosures this reads \(\max(2+2.70632481(1-\alpha),\,5.48553067-4.341285\,\alpha)<c(\alpha)<2+3.50249013(1-\alpha)\). For example, at \(\alpha=0.45\) it gives \(3.5319<c<3.9264\), which contains BD's certified \(c(0.45)\). It sharpens BD.3's \(c\le2+B_\alpha\), since \(B_{0.45}=115.2\).

**Proof of BF.1.** \(\int\tau^{-\alpha}f\,d\tau=\int e^{\varepsilon u}f\,du\). Shifting \(u=s_0+v\) multiplies the numerator and the denominator by the same \(e^{\varepsilon s_0}\). In the denominator, \(\int_{v<0}e^{\varepsilon v}\tfrac13\,dv=1/(3\varepsilon)\) (sympy). Multiplying the numerator and the denominator by \(\varepsilon\) gives the formula. *Smoothness:*
- as \(v\to-\infty\): \(0\le P\le\tau/18\) (\(T-2\le\tau/6\), \(\mathcal Y\le\tfrac13\)) and \(|\mathcal Y-\tfrac13|\le\tfrac{e^{0.1}+0.3}9\tau\) for \(\tau\le1\). The second bound follows from BE.3's \(\mathcal Y\ge e^{-\tau/10}/(3+\tau e^{\tau/10})\), since \(\tfrac13-\) that \(=(\tau e^x+3(1-e^{-x}))/(3(3+\tau e^x))\le(\tau e^{0.1}+3x)/9\), with \(x=\tau/10\) (sympy for the identity);
- as \(v\to+\infty\): \(\mathcal Y\le\tau/(3\psi_6(\tau))\) decays faster than any power (BD.4).

So \(e^{\varepsilon v}|v|^k\) times either integrand is integrable uniformly for \(\varepsilon\) near \([0,1]\). Differentiation under the integral is then legitimate, and \(h=\tfrac13+\varepsilon G=\varepsilon e^{-\varepsilon s_0}\int e^{\varepsilon u}\mathcal Y>0\) on \((0,1]\), with \(h(0)=\tfrac13\). \(\square\)

**Proof of BF.2 (certificate).**
- *Formula.* With \(\Phi=\varepsilon N\) and \(h=\tfrac13+\varepsilon G\), \(F''=(\Phi''-2F'h'-Fh'')/h\) and \(F'=(\Phi'-Fh')/h\) (sympy), where \(\Phi''=2N'+\varepsilon N''\), \(h''=2G'+\varepsilon G''\), and \(N^{(k)},G^{(k)}\) carry an extra \(v^k\).
- *Cover.* \([0,1]\) is covered by \(40\) closed intervals \(E_j=[j/40,(j+1)/40]\). On each one, \(F''\) is evaluated in mpmath.iv with \(\varepsilon\) itself an interval. The shift is \(s_0=2.25+\tfrac14\lfloor5\varepsilon_{\mathrm{mid}}+\tfrac12\rfloor\in\{2.25,\ldots,3.5\}\). Centring \(v\) near the mass of \(e^{\varepsilon u}P\) is what makes the interval evaluation tight. In a float prototype of the same cell scheme, the unshifted choice \(s_0=0\) gives an enclosure of \(F''(1)\) that contains \(0\) (method note, not a claim).
- *Cells.* There are \(14025\) cells of width \(1/500\) on \(u\in[-20,8.05]\). On each cell \(T\in[T_{60}(\tau_a),T_{60}(\tau_b)e^{\varepsilon^T}]\) and \(\mathcal Y\in[\mathcal Y_{60}(\tau_b)e^{-\varepsilon^{\mathcal Y}},\mathcal Y_{60}(\tau_a)e^{\varepsilon^{\mathcal Y}}]\), using \(T\) increasing, \(\mathcal Y\) decreasing, and BD.4's tail factors. The factors \(e^{\varepsilon v}\) and \(v^k\) are enclosed over the box (cell \(\times E_j\)).
- *Head* \(u<-20\): \(e^{\varepsilon v}\le1\) and \(\int_{-\infty}^{-m}|v|^ke^v\,dv=e^{-m}\sum_{j\le k}\tfrac{k!}{j!}m^j\) (sympy), with the two bounds above.
- *Tail* \(\tau\ge3125\): \(e^{\varepsilon v}\le\tau e^{-s_0}\), \(0<v^k\le(\log\tau)^k\le k!\,\tau\), \(T-2\le\tau/6\), and BD's \(\mathcal Y\le\tau^{-7}3125^8/(3P_8)\).
- *Result.* Every \(E_j\) gives a certified \(F''>0\), with smallest lower bound \(0.52046\) (on \(E_{39}\), i.e. \(\alpha\in[0,0.025]\)). \(\square\)

*Identity behind it, not a proof of sign:* under \(\nu_\alpha\propto\tau^{-\alpha}\mathcal Y\,d\tau\), \(c'(\alpha)=-\mathrm{Cov}(T,u)\) and \(c''(\alpha)=\kappa(T,u,u)=\mathbb E[(T-c)(u-\bar u)^2]\), the third joint cumulant (sympy). The first is \(<0\) by BE's monotone-likelihood-ratio argument. No such principle fixes the sign of the second, and that is why BF.2 is an interval certificate.

**Proof of BF.3.** \(c\) is \(C^2\) on \([0,1]\) with \(c''>0\). So it lies strictly above its tangents at \(0\) and at \(1\), and strictly below the chord. \(c'(1^-)=-F'(0)=-3S\) is consistent with BE.3: the point enclosure \(F'(0)\in[2.702670,2.713194]\) meets BE's \(3S\). \(c'(0^+)=-F'(1)\) is a point enclosure at \(\varepsilon=1\), with \(s_0=3.5\). \(\square\)

**Consistency (point enclosures, same cells).** \(2+F(1)\in[5.472917,5.515498]\) meets BE.4's \(c_0\). \(2+F(0.55)\in[3.719089,3.732225]\) meets BD's \(c(0.45)\).

**CONJECTURE BF.2 (floats only).** \(c''\) is not monotone. It rises from \(c''(1^-)\approx1.531\) to about \(1.593\) near \(\alpha\approx0.46\), then falls to \(c''(0^+)\approx1.557\). So \(c'''\) changes sign once. The certified enclosures of \(c''\) are too wide (\(\approx[0.52,2.6]\)) to decide this. No closed form for \(T\), \(\mathcal Y\), \(c\), \(S\), \(c_0\), \(c'(0^+)\) or \(c''\) is known.

**NOT THIS.** Theorem BF is the shape of one refinement-rate constant of a ratio of discrete Green-kernel entries on the combinatorial gasket, as a function of the fractional exponent. It does not choose the reading (Assumption BA/BC), does not select \(\alpha\), and does not change the limit angle \(\arctan(1/(3\sqrt3))\). \(c\) and its derivatives are rates, not physical coefficients. Not thrust.

**Witness, computed by `scripts/mesh_noncentrality_convex.py`** (about 2 min).
- *Sympy:* the \(F',F''\) quotient formulas; \(\varepsilon\int_{v<0}e^{\varepsilon v}/3=\tfrac13\); \(c'=-\mathrm{Cov}\), \(c''=\kappa(T,u,u)\); the head integrals for \(k\le3\); the \(|\mathcal Y-\tfrac13|\) identity.
- *Intervals:* \(F''>0\) on all \(40\) \(\varepsilon\)-intervals; point enclosures of \(F'(0)\), \(F(1)\), \(F'(1)\) and \(F(0.55)\), checked against BD and BE.
- *Floats (display only):* \(c''\) at five \(\alpha\) (Conjecture BF.2).

The script exits nonzero if any check fails.

## What was tried and does not follow

- **FAILED as a derivation of thrust.** Nothing above produces a nonzero \(\Delta F\) for a divergence-free informational tensor. The mathematics does not fail at a hidden algebraic step. The net is zero because the sum telescopes.
- **FAILED as a derivation of \(0.08\) or \(0.23\).** Those numerals never enter.
- **FAILED as a reading of the informational fork protocol.** The inequality \(T_{\mathrm{Red}} > T_{\mathrm{CIS}}\times 10^3\) is not a stress tensor and supplies no \(\mathrm{div}\,F\). It is not used.
- **NOT DERIVED.** Einstein gravity from entanglement. The essay in [-Entanglement-and-Emergence](https://github.com/beyond-repair/-Entanglement-and-Emergence) remains an essay. The only entanglement-facing consequence here is conditional: if an entanglement stress is divergenceless, Theorem A says its closed-surface signed flux is zero. Whether any concrete entanglement stress is divergenceless is **OPEN**.
- **OPEN.** Whether \((\nabla\Psi_{\mathrm{info}})^{ij}\) on the frozen \(0.45\) mesh has \(\sum\mathrm{div} \neq 0\) after the Maxwell piece is removed. That is a Stage 2 computation, not a symbolic gap. Under the Hessian reading, Theorem M reduces it to whether \(\Delta\Psi_{\mathrm{info}}\neq 0\) on the enclosing surface; with static massless \(\Psi=A_0\) on a vacuum surface it is \(0\) (Kept Failure M.1). Under the quadratic reading, Theorem N makes \(\mathcal{G}\) the volume integral of \(\partial_i\Psi\,\Delta\Psi\); with static massless \(\Psi=A_0\) that tensor is the Maxwell stress over \(\varepsilon_0\), so nothing remains after the Maxwell piece is removed, and an isolated device gives \(\mathcal{G}=0\) (Kept Failure N.1). Theorem O covers every shift-invariant, rotation-covariant, constant-coefficient reading with two derivatives: with static massless \(\Psi=A_0\), \(\mathcal{G}\) is either \(0\) or a trace flux \((\tfrac c2+d)\oint|E|^2n\) that moves with the surface and sums to \(0\) over all space (Kept Failure O.1). Theorem P covers every linear, rotation-covariant, constant-coefficient reading of any order: shift-invariant ones give \(\mathcal{G}=0\) on a vacuum surface (Kept Failure P.1), and the one undifferentiated term \(\delta^{ij}\Psi\) gives \(-\int_V E\), which depends on the surface (Kept Failure P.2). Theorem Q covers every shift-invariant polynomial reading of any order and degree: with static massless \(\Psi=A_0\) and an isolated device, \(\mathcal{G}\) is either \(0\) on every vacuum surface or moves with the surface (Kept Failure Q.1). Theorem R drops shift invariance: undifferentiated factors such as \(\Psi^2\delta^{ij}\) leave the same dichotomy, zero on every vacuum surface or surface-dependent (Kept Failure R.1). Theorem T takes the same polynomial family under a screened (Proca / Yukawa) equation: on-shell conserved readings still give \(\mathcal{G}=0\) on every vacuum surface, and every other reading is a surface-dependent mass remainder that decays exponentially (Kept Failure T.1); the massless escape routes of P.2 and R.4 are closed by Yukawa decay. Theorem U covers every first-derivative \(O(3)\)-covariant reading \(f(s)\delta+g(s)E\otimes E\) with \(f,g\in C^1([0,\infty))\), polynomial or not: on-shell conserved ones are exactly Maxwell plus \(\kappa\delta\) and give \(\mathcal{G}=0\); all others are surface-dependent and decay like \(R^{-2}\) (Kept Failure U.1); Maxwell removal leaves either \(0\) or that same surface-dependent remainder. Theorem V closes the \(\Psi+\nabla\Psi\) isotropic first-order gap: allowing undifferentiated \(\Psi\) in \(f(\Psi,s)\delta+g(\Psi,s)E\otimes E\) does not enlarge the conserved class beyond Maxwell+\(\kappa\delta\) (Kept Failure V.1). Theorem Y closes the position-dependent extension of that family: \(f(x,\Psi,s)\delta+g(x,\Psi,s)E\otimes E\) is on-shell conserved iff Maxwell plus constant \(\kappa\delta\) (Kept Failure Y.1). Theorem Z allows the pseudotensor term \(h\,\varepsilon^{ijk}\partial_k\Psi\): the conserved class grows only by \(\varepsilon\cdot\nabla K(\Psi)\), whose flux is \(0\) on every closed surface, and every other reading is surface-dependent (Kept Failure Z.1). Higher-derivative non-polynomial readings, nonlocal kernels, time-dependent fields, and any \(\Psi_{\mathrm{info}}\) not equal to \(A_0\) remain OPEN, as does the \(0.45\) mesh. Theorem AV takes up nonlocal kernels. For translation-invariant convolution readings (finite radial measure, e.g. heat or Yukawa smearing) and fractional readings \((-\Delta)^{s/2}\) in Maxwell form, the effective source leaks into the vacuum. Conservation on the vacuum region holds iff \(K=c\delta\). Otherwise \(\mathcal{G}\) is nonzero on some vacuum surfaces even for one point charge, but it moves with the surface and tends to \(0\) under dilation, and Newton's third law makes its total \(0\) (Kept Failure AV.1). Any reading with \(o(R^{-2})\) surface decay has no surface-independent nonzero value (AV.2). Theorem AW takes up the slow-tail smoothing readings \((-\Delta)^{-\tau}\Psi\) outside that decay class. The literal reading does not exist for \(\tau\ge(1+\ell)/2\). For its finite-part continuation, a surface-independent \(\mathcal{G}\) is \(0\) for every \(\tau>0\), and for non-integer \(\tau\) the reading is surface-dependent for every nonzero device. At the marginal \(\tau=\tfrac34\) (charged) and \(\tfrac54\) (dipole), \(\mathcal{G}\) does not depend on the radius but is linear in the audit centre, e.g. \(\tfrac{16}9Q^2(c-x_Q)\) (Kept Failure AW.1). Theorem AX proves Conjecture AW.4 for every \(\ell\): at \(\tau_\ell=(3+2\ell)/4\) a point \(2^\ell\)-pole gives exactly \(\mathcal{G}=Mc\) on every enclosing sphere, with \(M=\sum_j\mathcal M_{\ell-2j}[H_{\ell-2j}]\) in closed form. The literal reading never exists there (\(\tau_\ell-(1+\ell)/2=\tfrac14\)). But with the next moment present, the marginal value on spheres can be a nonzero vector that does not depend on the centre or the radius, while cubes give other values (Kept Failure AX.1). Non-Maxwell slow-tail tensors, time-dependent fields, and \(\Psi_{\mathrm{info}}\neq A_0\) remain OPEN. Theorem BA takes the question to the actual \(0.45\) mesh under one stated reading (Maxwell form \(f_c=q_cE_c\), \(q=L^{0.45}u\), corner-cell gradient, Green kernel \(A^+\)): the action–reaction piece has integrated divergence exactly \(0\), and the residual is \(D=\kappa_n(\alpha)\sum_c(d_c-\bar d)^2(P_c-\bar P)\) with certified \(\kappa_1,\ldots,\kappa_4<0\) (e.g. \(\kappa_4(0.45)=-18.526339511158477727\ldots\)), so \(\sum\mathrm{div}\neq0\) after Maxwell removal, the same on every vertex set containing the three corners. But it is a broken-Newton's-third-law defect of the discrete pair force, not an informational residual: the same rule on a translation-invariant torus gives a nonzero value with a one-sided stencil and exactly \(0\) with a centred one, it is almost all self-force, and it grows like \(2^n\) (Kept Failure BA.1). Theorem BB gives the exact decomposition \(\kappa_n=2^n(4/R_n^2)(3Y_n-\sigma_n)\) with \(0<Y_n\le5^{-n}R_n/2\). This proves Conjecture BA.2: \(2^{-n}\kappa_n(0.45)\to-4\sigma_\infty/R_\infty^2\in[-0.854646905120924,-0.854646905120848]\), and the limit is \(0\) for \(\alpha\ge\log3/\log5\). The cross (non-self) pair part tends to \(0\), so the growth and the limit are pure self-force. Removing the self force still leaves a certified nonzero residual at levels \(1\)–\(4\) (Kept Failure BB.1). Theorem BC tests the centred stencil on the gasket itself. By the \(S_3\) reflection law \(H_{kc}=R_{ck}H_{ck}\), the Maxwell piece is exactly the central (chord-parallel) part of the pair force, and \(\kappa_n=-(4/R_n^2)(s+2\sqrt3h_\perp)\) is self force plus non-central force. No vertex of the gasket admits an odd stencil, because its isometry group \(D_3\) has no half-turn. The Green–Gauss form of the centred stencil equals the one-sided rule at the charged corners, so \(\kappa^{\mathrm{GG}}_n=\kappa_n\): nonzero (certified at levels \(1\)–\(4\)), growing like \(2^n\) with limit \(K\). The certified non-centrality is \(\tan\theta_4=-0.2038773766517652\ldots\) (Kept Failure BC.1). Theorem BD proves Conjecture BC.2. The exact identity \(e_n=(2-t_n)y_n\) gives \(2/(3^n+1)<2-E_n/Y_n\le2/(3^n+1)+B_\alpha3^{-n}\), so \(\tan\theta_n\to-1/(3\sqrt3)\) from below for every \(\alpha\), and \(\kappa^{\mathrm{cross}}_n>0\) at every level. The \(3^{-n}\) constant exists and depends on \(\alpha\): certified \(c(0.45)\in[3.72198897,3.72914181]\) and \(c(0.8)\in[2.57126273,2.57352573]\) (Kept Failure BD.1). Theorem BE shows \(c(\alpha)\) is strictly decreasing, with \(c(1^-)=2\) and slope \(3S\in[2.70632481,2.70903118]\), while \(c(0^+)\) is finite, certified in \([5.48553067,5.50249013]\) (Kept Failure BE.1). Theorem BF proves that \(c\) is strictly convex, with certified \(c''\ge0.52046\) on \([0,1]\) (Conjecture BE.2), and \(c'(0^+)\in[-4.341285,-4.227515]\).
- **OPEN.** Distributional flux on the infinite gasket. Theorem A is the finite rectangle. It does not pass to a limit that has not been constructed. On the combinatorial gasket, before Theorem BG, only scalar audit norms had limits: \(\|F(n,\alpha)\|\to\mathcal L(\alpha)\) below \(\log3/\log5\) (Theorems AQ–AT), and at \(\gamma=\log3/\log5\) the renormalization \((n+\nu_\gamma)\|F(n,\gamma)\|\to K_\gamma\) with \(\nu_\gamma=1.3594681\ldots\) (Theorem AY; Kept Failure AY.1: the unshifted \(n\|F\|\) converges only at rate \(1/n\)). Above \(\gamma\), \((\rho^n+\nu_\alpha)\|F(n,\alpha)\|\to K_\alpha\) with \(\rho=5^\alpha/3\), e.g. \(\nu_{0.8}=-0.7274453\ldots\), \(K_{0.8}=0.6878667\ldots\) (Theorem AZ; Kept Failure AZ.1: the shift has the opposite sign to the endpoint's). No flux functional on the infinite gasket is constructed. The first non-scalar limit is [GASKET.md](GASKET.md) Theorem BG, the \(L^{\alpha}\)-Dirichlet extension itself. For \(\alpha\le\log3/\log5\) it tends to the constant \(\bar d\) at every non-corner vertex of \(V_*\), at rate \(\rho^n\) below \(\gamma\), with certified renormalised midpoint profile \(U(0.45)\in[0.143810693,0.143889868]\). The convergence is not uniform: a corner layer keeps the gap \(0.3853654863\ldots\) (Kept Failure BG.1). Above \(\gamma\) the limit is not flat (\(\kappa_{0.8}\in[0.086072435,0.086089770]\)). That limit is a function on \(V_*\), not a flux. [GASKET.md](GASKET.md) Theorem BH resolves the corner layer. For \(\alpha<\gamma\), \(u_n(2^{-n}\xi)\to\bar d+x_a\Theta(\xi)\) at every fixed layer vertex \(\xi\), with a universal profile built from the \(\alpha\)-free neighbour series \(\eta\) (\(r_\infty=1/(2+s-2\eta)\)). Certified \(\Theta(2,0)\in[0.340381375,0.340407195]\) at \(0.45\). The layer decays into the bulk like \(\rho^k\) per doubling, and it is absent for \(\alpha\ge\gamma\). It is a function on a graph, not a flux. [GASKET.md](GASKET.md) Theorem BI proves Conjecture BG.2. For \(\alpha>\gamma\), the oscillation of \(u_n\) on every level-\(N\) cell is at most \(\|x\|K_\alpha\rho^{-N}\), uniformly in \(n\) (certified \(K_{0.8}\le47.109\), from a cell recursion for the resolvent and the certified ceiling \(\sup s^{1-\gamma}r_\infty\le0.6031\)). So \(u_n\) converges uniformly to a continuous function on the gasket with the corner data as boundary values. That is a continuous function, not a flux. [GASKET.md](GASKET.md) Theorem BJ proves that this modulus is sharp (Conjecture BI.7). The exact identity \(1/r_n-1/r_\infty=\prod_{i<n}\beta(p^is)(1+2\eta(p^ns))\le3(3/5)^n\) and a certified bracket \(s^{1-\gamma}r_\infty(s)\in[0.600945014,0.600963487]\) near \(0\) give a corner deficit \(\rho^N\delta_N\ge c_-(\alpha)>0\) at every \(N\) (\(c_-(0.8)\ge0.6208\)). So the Hölder exponent \((\alpha-\gamma)/(1-\gamma)\) is optimal, and \(\liminf,\limsup\rho^N\delta_N(0.8)\in[1.135890346,1.135938783]\). It is a modulus of a continuous function, not a flux. [GASKET.md](GASKET.md) Theorem BK proves Conjecture BJ.7 with exponent \(c=1\). An exact change of variables carries level \(n+1\) back to level \(n\): \(u_{n+1}(x)-\bar d=\langle\omega_nF_n\rangle/\langle\omega_n\rangle\) under the level-\(n\) measure. The weight is \(\omega_n=(1+y/5)^{\alpha}/((1+2y/5)\tfrac53(\beta-er_n(p)))\ge1\), because \((1+2y/5)\tfrac53\beta=1-y(y+4)/(3(y+2)(y+5))\), and its mean is \(D'_{n+1}/D'_n\), with \(D'_n=\rho^{-n}D_n\) increasing. So \(|u_{n+1}-u_n|\le2\|x\|\log(D'_{n+1}/D'_n)\) at every vertex, and \(\sup_{V_n}|u_n-\bar u|\le2\|x\|\log(D'_\infty/D'_n)\le C_\alpha\|x\|\rho^{-n}\). The constants are certified: \(J_{0.8}\le1.188154\), \(C_{0.8}\le6.2192\), \(C_{0.75}\le12.5476\), \(C_{0.9}\le1.9725\), and the asymptotic constant at \(0.8\) is \(\le2.155275\) with AZ's \(A_{0.8}\) (\(2|\nu_{0.8}|=1.45489\), derived from AZ's \(C_{0.8}\)). Floats put the corner-neighbour error at \(0.6025\,\rho^{-n}\), so whether \(\rho^{-n}\) is attained stays Conjecture BK.6. It is a convergence rate for functions on graphs, not a flux. [GASKET.md](GASKET.md) Theorem BL proves that this order is attained. The deficit integrand is exactly monotone (\(1-\eta_N(z)=z/(z+2)+\beta(1-\eta_{N-1}(pz))\)), and \(\omega_k\ge1\). So each one-level covariance has a floor of order \(\rho^{-k}\) once the deficit is small. At \(0.8\), for every \(N\ge9\), \(\delta_N^{(n)}\) strictly increases in \(n\) and \(c_N\rho^{-n}\le\delta_N-\delta_N^{(n)}\le4.6644\rho^{-n}\), with \(c_N\ge0.1240\) for \(N\ge40\) (certified, AZ-free; floats \(0.6025\)). BK.6 fails as worded at \(N=0\), where \(\delta_0^{(n)}\equiv\tfrac32\) (Kept Failure BL.0). The limit constant remains Conjecture BL.6. It is a lower bound for functions on graphs, not a flux. [GASKET.md](GASKET.md) Theorem BM extends BL to every \(\alpha\in(\gamma,1)\). BL's ingredients are finite and positive there, and the \(N\)-uniform deficit ceiling tends to \(0\) like \(\rho^{-N}\). So for all large \(N\), \(\delta_N^{(n)}\) increases in \(n\) and \(\delta_N-\delta_N^{(n)}\ge c_N\rho^{-n}\). BK's order \(\rho^{-n}\) is therefore exact at every exponent above \(\gamma\). Certified thresholds: \(N\ge61,16,7,6,5\) at \(\alpha=0.7,0.75,0.85,0.9,0.95\) (python-flint arb; e.g. \(\liminf_Nc_N(0.75)\ge0.1555\)). Small \(N\) and the limit constant stay open. It is a lower bound for functions on graphs, not a flux. [GASKET.md](GASKET.md) Theorem BN removes the thresholds. The corner resolvent \(r_k\) is a Stieltjes function of unit mass, so \(r_k\le1/s\) and \(-r_k'\le r_k/s\). That makes BK's weight \(\omega_k\) strictly increasing at every level for every \(\alpha\ge1/17\): \(d_k'-\alpha d_k/(5+y)\le-P_\alpha/(3(y+2)^2(y+5)^2)\), \(P_\alpha=2\alpha y^3+(17\alpha-1)y^2+46\alpha y+40\alpha+10\). By Chebyshev's covariance inequality, \(\delta_N^{(n)}\) then increases in \(n\) for every \(N\ge1\). A pair form of the covariance gives \(\delta_N-\delta_N^{(n)}\ge c_N\rho^{-n}\) for every \(N\ge1\) at every \(\alpha\in(\gamma,1)\), with no deficit ceiling. Certified at \(0.7\)–\(0.95\), e.g. \(c_1(0.8)\ge1.758\cdot10^{-3}\), \(c_8(0.8)\ge1.914\cdot10^{-2}\). Only the limit constant of Conjecture BL.6 stays open. It is a statement about functions on graphs, not a flux. [GASKET.md](GASKET.md) Theorem BO finds the exact exponent range of BN's monotone weight. The corner-resolvent moments \(m_j=\tfrac12v^TL_k^jv\) are local, so they are the same integers \(1,2,6,22,92,418,1994,9780\) at every level \(k\ge3\). Exact geometric remainders then bound \(r_k\) and \(-r_k'\) on both sides. So \(\omega_k\) is increasing at every level iff \(\alpha\ge a^*\), with \(a^*\in[0.016778880172852,0.0167788801732]\) (exact-rational polynomial certificates, Sturm). This extends BN from \(1/17\) and with it BN's monotone deficit. Kept Failure BO.1: for \(\alpha<a^*\), \(\omega_k\) is not monotone at any level \(k\ge3\) (witness \(y^*=13.778505644\)), so BN's Chebyshev route is closed there. Whether the deficit still increases there is open. It is a statement about functions on graphs, not a flux. [GASKET.md](GASKET.md) Theorem BP answers that question: yes. Chebyshev is replaced by a single-crossing covariance. Below \(1/17\), \(\omega_k\) is strictly increasing on \([Y_\alpha,\infty)\), \(Y_\alpha=(1-17\alpha)/(2\alpha)\) (BN's majorant holds for every \(\alpha\ge0\)), and \(\omega_k<\tfrac85(1+y/5)^\alpha<3\cdot5^{-\alpha}=\rho^{-1}\) on \([0,Y_\alpha]\), while its mean \(D'_{k+1}/D'_k\) is \(>\rho^{-1}\) because the corner resolvent increases with the level. So \(\omega_k\) crosses its mean once and \(\mathrm{Cov}(\omega_k,G_k)>0\). Hence \(\delta_N^{(n)}\) strictly increases in \(n\) for every \(\alpha\in(0,1)\) and every \(N\ge1\), to \(1\) for \(\alpha\le\gamma\). The dip of BO.1 lies below the mean and does not affect the sign. It is a statement about functions on graphs, not a flux. [GASKET.md](GASKET.md) Theorem BQ gives the rate below \(\gamma\). BK.1's transfer applied to the numerator of \(1-\delta_N^{(n)}\) shows that the numerator increases strictly in \(n\), for every \(\alpha\), because \(\omega_n\ge1\) and \(\eta_N>0\). Below \(\gamma\), \(D_n\) increases to a finite \(D_\infty\). So \(1-\delta_N^{(n)}\) is of exact order \(\rho^n\) at every level, not only in the limit, and \(\rho^{-n}(1-\delta_N^{(n)})\to C_N(\alpha)=P_N/D_\infty>0\). Certified (python-flint arb): \(0.050648\,\rho^n<1-\delta_1^{(n)}(0.45)<0.118859\,\rho^n\) for every \(n\ge1\), \(C_1(0.45)\in[0.0718892,0.0719611]\), \(C_1(0.01)\in[0.0078312,0.0078363]\). Kept Failure BQ.1: \(\rho^{-n}(1-\delta_1^{(n)})\) is neither monotone nor one-sided. At \(0.45\) it rises from \(n=1\) to \(n=2\), then falls toward its limit from above; at \(0.3\) it crosses its limit between \(n=3\) and \(n=4\). It is a statement about functions on graphs, not a flux. [GASKET.md](GASKET.md) Theorem BR gives the depth dependence. Going one step deeper in \(N\) while moving the argument one step along \(p\) can only raise the neighbour profile: \(\eta_{N+1}(w)-\eta_N(p(w))\ge w/(w^2+5w+2)\), from the ceiling \(\eta_N\le2/(2+z)\), now proved by induction. So \(\delta_N^{(n)}\) strictly decreases in \(N\) at every level \(n\) and every \(\alpha\in(0,1)\). Below \(\gamma\), the rescaling \(\tau=5^Nt\) writes \(\rho^NP_N\) as an integral whose two factors both increase in \(N\) and are bounded by BJ.2's power bracket. So \(\rho^NC_N\) strictly increases in \(N\) to a finite limit \(c_\infty(\alpha)\in[g_{\rm lo},g_{\rm hi}]\cdot M_\alpha/D_\infty\), \(M_\alpha=\int t^{\gamma-1-\alpha}\eta(\psi(t))\,dt\le2^{\gamma-\alpha}\pi/\sin\pi(\gamma-\alpha)\), and \(C_N\) has exact order \(\rho^{-N}\). Certified (python-flint arb): \(c_\infty(0.45)\in[0.520387,0.520824]\) and \(0.049440\le\rho^NC_N(0.45)\le0.520824\) for every \(N\ge1\); \(c_\infty(0.01)\in[0.010439,0.010444]\). It is a statement about functions on graphs, not a flux. [GASKET.md](GASKET.md) Theorem BS gives the rate. Rescaling BJ.1 by \(5^{-M}\) and letting \(M\to\infty\) gives an exact identity, \(1/\mathcal R=1/Q_\infty+\Phi\). Here \(Q_\infty(u)=\lim(3/5)^Mr_\infty(5^{-M}u)\) exists without AT.1 and lies in \([g_{\rm lo},g_{\rm hi}]u^{\gamma-1}\) for every \(u\), and \(\Phi=(1+2\eta\circ\psi)\prod_{k\ge1}\tfrac53\beta(\psi(5^{-k}u))\) is explicit. So \(c_\infty-\rho^NC_N=\rho^N\chi_N\), \(\chi_N=(A_N+B_N)/D_\infty\), where \(A_N\) increases to \(J=\int u^{-\alpha}(Q_\infty-\mathcal R)\) and \(B_N\le(3/5)^{N-1}B_1\). Hence \(\chi_N\to\chi_\infty=J/D_\infty>0\), \(C_N=c_\infty\rho^{-N}-\chi_\infty+o(1)\), and the rate \(\rho^N\) is exact. Certified (python-flint arb): \(\chi_\infty(0.45)\in[0.633344,0.633869]\), and \(\chi_N(0.45)\in[0.618141,0.649097]\) for every \(N\ge6\). Individual \(C_N(0.45)\) for \(N=1,\ldots,6\), e.g. \(C_2\in[0.458110,0.458442]\) and \(C_6\in[4.285424,4.288550]\). \(C_1\) agrees with BQ.3's independent bracket at all four exponents, and the new brackets are tighter. Whether \(\chi_N\) decreases in \(N\) is Conjecture BS.4 (certified for \(N=1\)–\(5\) at \(0.45\)). It is a statement about functions on graphs, not a flux. Theorem BT pins the one unknown in that identity: \(u^{1-\gamma}Q_\infty(u)=\Theta_{AT}(\psi(u))\), because \(s^{1-\gamma}r_\infty\) is nonincreasing along every \(p\)-orbit; the profile is continuous, has AU's mean, is certified non-constant (Kept Failure BT.5: amplitude \(\ge1.465\cdot10^{-5}\)), and its maximum is the sharp constant \(\sup_s s^{1-\gamma}r_\infty(s)\in[0.6009615785,0.6009631774]\). Theorem BU shows that profile is AU's \(\theta\) itself. AU's Koenigs coordinate \(\lim5^kq^k\) is the inverse of BD.4's \(\psi\), so \(\hat F_k=\hat\theta_k\) for every mode. Convexity of \(r_\infty\) then gives second-order cell bounds: \(\max F\in[0.6009615785971,0.6009615785983]\), \(\min F\in[0.6009469237907,0.6009469237943]\), \(\langle F\rangle\) to \(2.4\cdot10^{-12}\), and \(\hat F_1,\hat F_2\ne0\) certified from \(F\)'s own nodes, matching AU. Kept Failure BU.4: although the extremes sit half a period apart at node resolution, \(F\) is not half-period antisymmetric (\((\max F+\min F)/2-\langle F\rangle\ge2.43\cdot10^{-11}\)). [GASKET.md](GASKET.md) Theorem BV shows that \(F=\theta\) is holomorphic on the strip \(|\operatorname{Im}\log_5u|<\pi/\ln5\). Every Stieltjes function obeys \(|f(z)|\le f(|z|)/\cos(\arg z/2)\), and that bound passes by Vitali to \(r_\infty\) and \(Q_\infty\). So \(|\hat\theta_k|\le G^*\sqrt{1+T_k^2}e^{-T_k\arctan T_k}\le eG^*(1+T_k)e^{-2\pi^2|k|/\ln5}\), \(T_k=4\pi|k|/\ln5\), which is the decay rate the floats suggested (upper half; the matching lower rate is Conjecture BV.6). AU's routine with more cells certifies \(\hat\theta_k\ne0\) for every \(k\le11\), to \(\pm8.2\cdot10^{-63}\). With the tail bound, \(F\) has exactly one maximum and one minimum per period, and \(G^*=\max F=0.600961578597101089042056184\ldots\) and \(\min F=0.600946923794265612104226188\ldots\) are certified to \(10^{-60}\). Kept Failure BV.5: the extremes are not half a period apart (\(t_{\max}-t_{\min}-\tfrac12=-1.624019\cdot10^{-7}\)), although both node grids put them exactly \(2^{17}\) of \(2^{18}\) cells apart. Functions on a graph, not a flux. Not thrust.
- **OPEN.** The renormalized local propagator for \(W(x)\), already open in the Ware derivation ledger. Theorem E constrains the force ledger of a variable weight. It does not construct \(Z_{\mathrm{ren}}\).

## Not a consequence

No laboratory force, no energy-extraction law, no selected \(W\), no confirmation of a pinch percentage, and no change to the Stage 1 freeze. The symbolic chain may stay written as \(\Delta F = W\chi\mathcal{G}\). This note fixes what \(\mathcal{G}\) has to equal.

## Reproduce

```bash
python3 scripts/flux_identity.py
python3 scripts/hessian_flux_reduction.py
python3 scripts/quadratic_flux_reduction.py
python3 scripts/general_reading_flux.py
python3 scripts/linear_reading_flux.py
python3 scripts/nonlinear_reading_flux.py
python3 scripts/nonshift_reading_flux.py
python3 scripts/screened_reading_flux.py
python3 scripts/first_derivative_reading_flux.py
python3 scripts/psi_gradient_reading_flux.py
python3 scripts/position_dependent_reading_flux.py
python3 scripts/hemitropic_reading_flux.py
python3 scripts/nonlocal_reading_flux.py
python3 scripts/slow_tail_reading_flux.py
python3 scripts/multipole_marginal_flux.py
python3 scripts/mesh_maxwell_residual.py
python3 scripts/mesh_maxwell_limit.py
python3 scripts/mesh_centred_stencil.py
python3 scripts/mesh_noncentrality_limit.py
python3 scripts/mesh_noncentrality_monotone.py
python3 scripts/mesh_noncentrality_convex.py
python3 scripts/gasket_fractional_extension_limit.py
python3 scripts/gasket_fractional_corner_layer.py
python3 scripts/gasket_fractional_uniform_limit.py
python3 scripts/gasket_fractional_corner_sharpness.py
python3 scripts/gasket_fractional_uniform_rate.py
python3 scripts/gasket_fractional_rate_sharpness.py
python3 scripts/gasket_fractional_rate_sharpness_alpha.py
python3 scripts/gasket_fractional_monotone_weight.py
python3 scripts/gasket_fractional_weight_threshold.py
python3 scripts/gasket_fractional_deficit_small_alpha.py
python3 scripts/gasket_fractional_deficit_rate_small_alpha.py
python3 scripts/gasket_fractional_deficit_depth.py
python3 scripts/gasket_fractional_deficit_correction.py
python3 scripts/gasket_fractional_log_periodic_profile.py
python3 scripts/gasket_fractional_profile_fourier.py
python3 scripts/gasket_fractional_profile_analytic.py
```

The script writes `witness.json` and exits nonzero if the integers quoted above move.

## Run it

Python 3.10+ and NumPy. SciPy is only needed for the optional matrix-free audit.

```bash
git clone https://github.com/beyond-repair/informational-flux-identity.git
cd informational-flux-identity
python3 -m venv .venv && . .venv/bin/activate
pip install -e ".[test]"        # numpy, plus pytest and scipy for the tests

flux-identity witness           # Theorem B witness: net 0, right-face share 349/366; +17 source gives net 17
flux-identity reproduce         # runs every reproduction script below in a scratch copy; prints 51/51 passed
flux-identity reproduce --matrix-free-level 5   # also runs the SciPy matrix-free audit at level 5
flux-identity list              # what each script checks
python -m pytest                # identity, CLI, reproduction, and matrix-free vs dense checks
```

`flux-identity check FILE.json` applies Theorem A to your own integer face-flux array. Give either `{"fx": [[...]], "fy": [[...]]}` with `fx` of shape `(nx+1, ny)` and `fy` of shape `(nx, ny+1)`, or `{"potential": [[...]]}` on the `(nx+1) x (ny+1)` vertex grid (always divergence-free). It prints the summed divergence, the signed and absolute flux on each face, and whether the two totals agree. Use `-` to read from stdin and `--json` for machine output.

Exit codes: `0` ok, `1` a recorded figure moved or a reproduction script failed, `2` bad input. Without installing, `python3 -m scripts witness` works from the checkout root, and the original `python3 scripts/<name>.py` commands are unchanged. `flux-identity reproduce` runs on a temporary copy, so it never rewrites the committed `witness.json`. The level-9 and level-10 scripts check recorded digits from long matrix-free runs; they do not redo those runs. None of this output is thrust.

## Sequel

On the combinatorial gasket, harmonic corner currents are a neutral dipole and shrink by exactly \(3/5\) at each refinement. For the spectral fractional operator, corner currents stay neutral (Theorem H) but do not inherit the factor \(3/5\) (Kept Failure H.1). The audit norm at \(\alpha=0.45\) does not tend to infinity (Theorem J) and does not tend to \(0\) (Theorem AH, below); it strictly decreases to a limit \(\mathcal L\in[0.819255,0.822959]\) (Theorem AQ, below, proves Conjecture J.1; Theorem AR gives the exact approach rate \((5^{\alpha}/3)^n\), the exact linear slope of \(R_{\gamma}(n)\) at \(\alpha=\log3/\log5\), and the exact decay ratio \(3\cdot5^{-\alpha}\) above it; Kept Failures J.1–J.4 had rejected four earlier ways of closing it; Theorem AS certifies \(\mathcal L=0.8211386271335\ldots\) to 13 digits by a complete-monotonicity interpolation sandwich and proves Conjecture J.2, \(\|F(n,0.10)\|>1.410\) at every level; Kept Failure AS.1 shows that level-10 extrapolation with the exact ratio misses \(\mathcal L\) by \(2.8\cdot10^{-4}\); Theorem AT certifies AR's rate constants to about 17 digits by an exact change of variables to level \(40\), e.g. \(c_R(0.45)=0.4978026507635603\ldots\), \(c_S(0.45)=0.9270926385694613\ldots\), slope \(c_R(\log3/\log5)=0.5171658346305052\ldots\), and \(\|F(n,0.45)\|-\mathcal L\sim0.5202787306\,(5^{0.45}/3)^n\); Kept Failure AT.1 certifies that \(s^{1-\gamma}r_\infty(s)\) oscillates log-periodically, so that slope is not a pure-power constant; Theorem AU proves Conjecture AT.2, \(c_R(\gamma)=\frac{2\sin\pi\gamma}{\pi}\ln5\,\langle\Theta\rangle\) with \(\langle\Theta\rangle\) the log-period mean of that oscillation, certifies \(\langle\Theta\rangle=0.60095425116877360782370515009999\ldots\) and so \(c_R(\gamma)=0.51716583463050524392598048567312\ldots\) to 34 digits (inside AT's independent enclosure), and shows that the Mellin transform of \(r_\infty\) has simple poles at \(1-\gamma+2\pi ik/\ln5\) with residues the Fourier coefficients of the oscillation, certified nonzero for \(|k|\le4\); Kept Failure AU.1: the oscillation is not a single harmonic; Theorem AY certifies the intercept of the linear growth at \(\gamma\), \(R_\gamma(n)=c_R(\gamma)n+C_\gamma+\varepsilon_n\) with \(C_\gamma=0.70307046242185742992904816\ldots\) and \(|\varepsilon_n|\le4.7\cdot10^{-28}\) for \(n\ge60\), because AT's bracket weights do not depend on the level and the per-level defects therefore sum to one integral of \(r_\infty-r_n\); with the certified limit share \(\sigma_\infty(\gamma)=1.1519186116740\ldots\) this gives \((n+\nu_\gamma)\|F(n,\gamma)\|=K_\gamma+O(5^{-0.65n})\), \(\nu_\gamma=1.3594681151437\ldots\), \(K_\gamma=3.4023607638326\ldots\); Kept Failure AY.1: without the level shift \(\nu_\gamma\) the error is \(\Theta(1/n)\), not exponential; Theorem AZ gives the supercritical law: for \(\alpha>\gamma\), \(R_\alpha(n)=A_\alpha\rho^n+C_\alpha+O(\rho^n5^{-0.65n})\), where the intercept \(C_\alpha=\frac{2\sin\pi\alpha}{\pi}M(1-\alpha)\) is the Mellin transform of \(r_\infty\) continued across its simple pole at \(\gamma\), so \((\rho^n+\nu_\alpha)\|F(n,\alpha)\|=K_\alpha+O(5^{-0.65n})\) and \(\|F(n,\alpha)\|\sim K_\alpha(3\cdot5^{-\alpha})^n\), with certified \(C_{0.8}=-1.74325717839207572395\ldots\), \(\nu_{0.8}=-0.72744538444215\ldots\), \(K_{0.8}=0.68786675841399\ldots\); Kept Failure AZ.1: \(C_\alpha<0\), so \(\rho^n\|F(n,0.8)\|\) decreases to \(K_{0.8}\) from above, opposite to the endpoint, with error \(\Theta((3\cdot5^{-0.8})^n)\) if unshifted). Dirichlet energy controls that audit norm at every \(\alpha\in(0,1)\) (Theorem L), so energy collapse would force the limit to \(0\). Theorem AA proves that collapse does not happen for any \(\alpha<1/2\). A boundary path is a subgraph, so Löwner–Heinz plus a closed-form path resistance gives \(Q_{0.45}(u_n)\ge 0.128335\) and fractional-current norm \(\|F_{\mathrm{frac}}\|\ge 0.084015\) at every level. That closes the energy-collapse route to zero (Kept Failure AA.1). Theorem AB factors the combinatorial-to-fractional share as \(\sigma_n=\rho_n/w_n\) and, with AA, equates \(\|F\|\to 0\) with vanishing of the two-neighbor fractional share \(\rho_n\); Kept Failure AB.1 rejects geometric collapse of \(\rho\). Theorem AC extends the energy floor from \(\alpha<1/2\) to every \(\alpha<\log 3/\log 5\), using cell resistances and the spectral mass of \(e_a-e_b\) instead of a path, so AB's equivalence holds on that whole range; by Theorem K no floor exists above it (Kept Failure AC.1), and at \(\alpha=\log 3/\log 5\) itself Theorem AD proves the exact pairing \(Q_{\alpha}R_{\alpha}=7/3\) and \(c/n\le Q\le C/n\), so the energy and the audit norm tend to \(0\) there (Kept Failure AD.1: the floor does not reach the endpoint). Theorem AE factors \(\rho_n=\lambda_n/(\lambda_n+(1-\lambda_n)\gamma_n)\) with sharp ceiling \(w_n\le8^{\alpha}/2\) and the resistance identity tying \(\rho\) to AD's \(R_{\alpha}\); Kept Failure AE.1 rejects geometric wipeout of the neighbor-degree share \(\lambda_n\). Theorem AF floors that share uniformly (\(\lambda_n\ge\alpha\,2^{-1-\alpha}\), sharpened by Jensen to \(\alpha\,2^{\alpha-1}\)), so below \(\log3/\log5\) a vanishing audit norm would need the gap ratio \(\gamma_n\to\infty\); Kept Failure AF.1 rejects geometric blowup of \(\gamma_n\), which decreases on the computed levels. Theorem AG ceilings the same share (\(\lambda_n\le 1-\alpha(1-\alpha)/2^{3-\alpha}<1\)), so vanishing is equivalent to \(\gamma_n\to\infty\); Kept Failure AG.1 rejects geometric climb of \(\lambda\) toward \(1\). Theorem AH closes the gap ratio: a heat-kernel comparison gives \(w_{pj}\ge\kappa(\alpha)w_{aj}\) for every vertex \(j\) off the corner patch, and \(L^{\alpha}\)-harmonicity of \(u_a-u\ge0\) at the two corner neighbors then gives \(I_{\mathrm{far}}\le(4^{\alpha}/(2\kappa))I_a\), so \(\gamma_n\le\Gamma_{AH}(\alpha)<\infty\) and \(\|F(n,\alpha)\|\ge\sigma_{\min}(\alpha)\|F_{\mathrm{frac}}(n,\alpha)\|>0\) at every level for every \(\alpha<\log3/\log5\); at \(\alpha=0.45\), \(\|F(n,0.45)\|\ge 6.44\cdot10^{-4}\), so the combinatorial audit norm does not vanish. Kept Failure AH.1 rejects the shortcut through Theorem J's \(3/2\) gap cap, which gives no floor from level \(3\). Theorem AI sharpens AH by a first-step heat-kernel identity at the degree-\(2\) corner: \(w_{pj}+w_{qj}\ge2w_{aj}\) on Far, so the neighbor Harnack holds with \(\kappa\) replaced by \(1\), \(I_{\mathrm{far}}\le(4^{\alpha}/2)I_a\), \(\gamma_n\le8/(\alpha(1-\alpha))\), and \(\|F(n,0.45)\|\ge0.051988\) at every level; Kept Failure AI.1 rejects reading this as Conjecture J.1's \(0.70\). Theorem AJ keeps the neighbor-edge weight \(w_{pq}\ge\alpha4^{\alpha-1}\) in the exact harmonicity identity at \(p,q\), so \(I_{\mathrm{far}}\le((4^{\alpha}-w_{pq})/2)I_a\le H_{AJ}I_a\) with \(H_{AJ}=(4^{\alpha}-\alpha4^{\alpha-1})/2\), \(\gamma_n\le2(4-\alpha)/(\alpha(1-\alpha))\), and \(\|F(n,0.45)\|\ge0.055599\); Kept Failure AJ.1 rejects reading this as \(0.70\). Theorem AK writes \(L=3\,\mathrm{diag}(c)-N^{\mathsf T}N\) through the cell incidence matrix, so \(\mathrm{spec}(L)\subset[0,6]\), and a two-moment Jensen bound on \([0,6]\) lowers the diagonal ceilings to \((L^{\alpha})_{aa}\le6^{\alpha}/9+(8/9)(3/2)^{\alpha}\) and \((L^{\alpha})_{pp}\le(6^{\alpha}+2^{\alpha})/2\), giving \(\|F(n,0.45)\|\ge0.058941\); Kept Failure AK.1 shows the one-step neighbor-Harnack chain stays below \(0.70\) from level \(4\) even with every constant at its live value, so that family cannot prove J.1. Theorem AL replaces AA's path ceiling \(C(\alpha)\) by the Hermite–Hadamard limit \(R_{\infty}(\alpha)\) for \(\alpha\ge\alpha_*=(14-8\sqrt{2})/34\), giving \(\|F_{\mathrm{frac}}\|\ge0.111208\) and \(\|F(n,0.45)\|\ge0.078019\); Kept Failure AL.1 shows the sharp path ceiling still cannot reach \(0.70\). Theorem AM replaces the single side by the whole boundary cycle \(C_{3\cdot2^n}\): Fourier plus a heat-kernel periodization bound its fractional corner resistance by \((1-\alpha)R_{\infty}(\alpha)\) for every \(\alpha<1/2\), and the cycle's corner symmetry gives the three-corner pairing \(Q\ge(7/3)/R^{\mathrm{cyc}}\), so \(\|F_{\mathrm{frac}}\|\ge0.209685\) and \(\|F(n,0.45)\|\ge0.147106\); Kept Failure AM.1 shows this still cannot reach \(0.70\), and that even the live energy times \(\sigma_{AK}\) stays below \(0.70\) from level \(2\), so a better resistance ceiling alone cannot prove J.1. Theorem AN replaces the cycle by the boundary fence/collar (Bloch ceiling \(R^{\mathrm{fen}}_{\infty}+1/2\approx5.002104\)), giving \(\|F_{\mathrm{frac}}\|\ge0.305377\) and \(\|F(n,0.45)\|\ge0.214239\); Kept Failure AN.1 shows this still cannot reach \(0.70\). Theorem AO thickens that strip by one layer (depth-2 collar; period-2 Bloch ceiling \(R^{\mathrm{d2}}_{\infty}+11/20\approx4.820047\)), giving \(\|F_{\mathrm{frac}}\|\ge0.316911\) and \(\|F(n,0.45)\|\ge0.222331\); Kept Failure AO.1 shows the gain is thin (\(\approx1.038\times\) AN) and still far below \(0.70\). Theorem AP leaves the Harnack family entirely: \(Lu=L^{1-\alpha}(L^{\alpha}u)\) and the corner currents of Theorem AD give the exact share identity \(\sigma_n=(L^{1-\alpha})_{aa}+\omega_n\) with \(\omega_n>0\), and a Gauss–Radau bound from the level-independent corner moments \(1,2,6,22,92\) gives \(\sigma_n\ge1.310267>1\) at \(\alpha=0.45\) (against AK's \(0.701556\)), so \(\|F(n,0.45)\|\ge0.274743\) with AM's floor and \(\ge0.415238\) with AO's floor as stated; Kept Failure AP.1 shows no share floor can lift the proved fractional floors to \(0.70\), so what is left is a resistance ceiling \(R_{\alpha}(n)\le2.859237\). Those floors are crude (live \(\|F(6)\|\approx0.880\)). Theorem AQ replaces them with a limit theorem. Eliminating the midpoints of every level-\(1\) cell gives an exact one-dimensional recursion for the corner resolvent, \(r_{n+1}(s)=\Phi_s(r_n(s(s+5)))\) with \(r_n(s)=\tfrac12v^{T}(L_n+s)^{-1}v\). Through Stieltjes integrals for \(R_{\alpha}(n)\) and \(\sigma_n\), this makes \(R_{\alpha}(n)\) strictly increasing and \(\sigma_n\) strictly decreasing at every \(\alpha\in(0,1)\), so \(\|F(n,\alpha)\|\) strictly decreases and converges. An interval-arithmetic certificate at \(0.45\) gives \(R_{\alpha}(n)<2.518955\), \(\sigma_n>1.350986\), and \(\lim\|F(n,0.45)\|\in[0.819255,0.822959]\), so Conjecture J.1 holds. Kept Failure AQ.1: there is no uniform ceiling and no positive limit for \(\alpha\ge\log3/\log5\). For \(\alpha>\log 3/\log 5\) that same audit norm tends to \(0\) (Theorem K). For geometric weights \(1/d^2\) on the same finest-edge build, Theorem I gives exact growth by \(12/5\) so \(\|F\|\to\infty\) (Kept Failure I.1). Uniform multiplicative weights \(w=\lambda^n\) give the harmonic trichotomy of Theorem S (finite nonzero only at the resistance value \(\lambda=5/3\)); resistance weights do not freeze the fractional \(\alpha=0.45\) audit (Kept Failure S.1). Multi-scale multiplicative weights (edges at every generation) give the geometric-series trichotomy of Theorem W (finite nonzero for every \(\lambda<5/3\); resistance diverges like \((n+1)\sqrt{21}/2\), Kept Failure W.1). Spectral \(L^{\alpha}\) on that multi-scale graph (Theorem X) keeps neutrality and the Schur line, recovers W at \(\alpha=1\), and shows nested harmonicity is special to \(\alpha=1\); the \(\lambda=1\) freeze does not transfer to \(\alpha=0.45\) (Kept Failure X.1), and multi-scale fractional does not coincide with finest-edge fractional (Kept Failure X.2). None of these norms is thrust. Proofs and witnesses are in [GASKET.md](GASKET.md). Run `python3 scripts/gasket_corner_current.py`, `python3 scripts/gasket_fractional_currents.py`, `python3 scripts/gasket_geometric_currents.py`, `python3 scripts/gasket_uniform_weights.py`, `python3 scripts/gasket_multiscale_weights.py`, `python3 scripts/gasket_multiscale_fractional.py`, `python3 scripts/gasket_fractional_limit.py`, `python3 scripts/gasket_fractional_energy_floor.py`, `python3 scripts/gasket_fractional_share.py`, `python3 scripts/gasket_fractional_capacity_floor.py`, `python3 scripts/gasket_fractional_endpoint.py`, `python3 scripts/gasket_fractional_degree_share.py`, `python3 scripts/gasket_fractional_lambda_floor.py`, `python3 scripts/gasket_fractional_lambda_ceiling.py`, `python3 scripts/gasket_fractional_gamma_ceiling.py`, `python3 scripts/gasket_fractional_harnack_sum.py`, `python3 scripts/gasket_fractional_wpq_harnack.py`, `python3 scripts/gasket_fractional_moment_jensen.py`, `python3 scripts/gasket_fractional_path_ceiling.py`, and `python3 scripts/gasket_fractional_cycle_ceiling.py`, and `python3 scripts/gasket_fractional_fence_ceiling.py`, `python3 scripts/gasket_fractional_depth2_ceiling.py`, `python3 scripts/gasket_fractional_moment_share.py`, `python3 scripts/gasket_fractional_decimation.py`, `python3 scripts/gasket_fractional_limit_digits.py`, and `python3 scripts/gasket_fractional_rate_constants.py`.

---

<div align="center">

**REWRITE · BUILD · TRANSCEND**

**William (Brian) Ware** · [Atomic Dream Labs](https://github.com/beyond-repair)  
Governing source: [ADL-Governance](https://github.com/beyond-repair/ADL-Governance) · [Claim levels 0–5](https://github.com/beyond-repair/ADL-Governance/blob/main/docs/CLAIM_VALIDATION.md)

</div>
