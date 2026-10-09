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

**NOT THIS.** Theorem V does not compute \(\Psi_{\mathrm{info}}\) on the frozen \(0.45\) mesh and does not decide which reading the freeze intends. It does not cover dependence on second or higher derivatives non-polynomially, nonlocal kernels, position-dependent coefficients, pseudotensor contractions, time-dependent or radiating fields, \(d\neq3\), or any \(\Psi_{\mathrm{info}}\neq A_0\). A nonzero \(\mathcal{G}\) on a finite surface is not thrust.


## What was tried and does not follow

- **FAILED as a derivation of thrust.** Nothing above produces a nonzero \(\Delta F\) for a divergence-free informational tensor. The mathematics does not fail at a hidden algebraic step. The net is zero because the sum telescopes.
- **FAILED as a derivation of \(0.08\) or \(0.23\).** Those numerals never enter.
- **FAILED as a reading of the informational fork protocol.** The inequality \(T_{\mathrm{Red}} > T_{\mathrm{CIS}}\times 10^3\) is not a stress tensor and supplies no \(\mathrm{div}\,F\). It is not used.
- **NOT DERIVED.** Einstein gravity from entanglement. The essay in [-Entanglement-and-Emergence](https://github.com/beyond-repair/-Entanglement-and-Emergence) remains an essay. The only entanglement-facing consequence here is conditional: if an entanglement stress is divergenceless, Theorem A says its closed-surface signed flux is zero. Whether any concrete entanglement stress is divergenceless is **OPEN**.
- **OPEN.** Whether \((\nabla\Psi_{\mathrm{info}})^{ij}\) on the frozen \(0.45\) mesh has \(\sum\mathrm{div} \neq 0\) after the Maxwell piece is removed. That is a Stage 2 computation, not a symbolic gap. Under the Hessian reading, Theorem M reduces it to whether \(\Delta\Psi_{\mathrm{info}}\neq 0\) on the enclosing surface; with static massless \(\Psi=A_0\) on a vacuum surface it is \(0\) (Kept Failure M.1). Under the quadratic reading, Theorem N makes \(\mathcal{G}\) the volume integral of \(\partial_i\Psi\,\Delta\Psi\); with static massless \(\Psi=A_0\) that tensor is the Maxwell stress over \(\varepsilon_0\), so nothing remains after the Maxwell piece is removed, and an isolated device gives \(\mathcal{G}=0\) (Kept Failure N.1). Theorem O covers every shift-invariant, rotation-covariant, constant-coefficient reading with two derivatives: with static massless \(\Psi=A_0\), \(\mathcal{G}\) is either \(0\) or a trace flux \((\tfrac c2+d)\oint|E|^2n\) that moves with the surface and sums to \(0\) over all space (Kept Failure O.1). Theorem P covers every linear, rotation-covariant, constant-coefficient reading of any order: shift-invariant ones give \(\mathcal{G}=0\) on a vacuum surface (Kept Failure P.1), and the one undifferentiated term \(\delta^{ij}\Psi\) gives \(-\int_V E\), which depends on the surface (Kept Failure P.2). Theorem Q covers every shift-invariant polynomial reading of any order and degree: with static massless \(\Psi=A_0\) and an isolated device, \(\mathcal{G}\) is either \(0\) on every vacuum surface or moves with the surface (Kept Failure Q.1). Theorem R drops shift invariance: undifferentiated factors such as \(\Psi^2\delta^{ij}\) leave the same dichotomy, zero on every vacuum surface or surface-dependent (Kept Failure R.1). Theorem T takes the same polynomial family under a screened (Proca / Yukawa) equation: on-shell conserved readings still give \(\mathcal{G}=0\) on every vacuum surface, and every other reading is a surface-dependent mass remainder that decays exponentially (Kept Failure T.1); the massless escape routes of P.2 and R.4 are closed by Yukawa decay. Theorem U covers every first-derivative \(O(3)\)-covariant reading \(f(s)\delta+g(s)E\otimes E\) with \(f,g\in C^1([0,\infty))\), polynomial or not: on-shell conserved ones are exactly Maxwell plus \(\kappa\delta\) and give \(\mathcal{G}=0\); all others are surface-dependent and decay like \(R^{-2}\) (Kept Failure U.1); Maxwell removal leaves either \(0\) or that same surface-dependent remainder. Theorem V closes the \(\Psi+\nabla\Psi\) isotropic first-order gap: allowing undifferentiated \(\Psi\) in \(f(\Psi,s)\delta+g(\Psi,s)E\otimes E\) does not enlarge the conserved class beyond Maxwell+\(\kappa\delta\) (Kept Failure V.1). Higher-derivative non-polynomial readings, nonlocal kernels, position-dependent coefficients, time-dependent fields, and any \(\Psi_{\mathrm{info}}\) not equal to \(A_0\) remain OPEN, as does the \(0.45\) mesh.
- **OPEN.** Distributional flux on the infinite gasket. Theorem A is the finite rectangle. It does not pass to a limit that has not been constructed.
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
flux-identity reproduce         # runs every reproduction script below in a scratch copy; prints 18/18 passed
flux-identity reproduce --matrix-free-level 5   # also runs the SciPy matrix-free audit at level 5
flux-identity list              # what each script checks
python -m pytest                # identity, CLI, reproduction, and matrix-free vs dense checks
```

`flux-identity check FILE.json` applies Theorem A to your own integer face-flux array. Give either `{"fx": [[...]], "fy": [[...]]}` with `fx` of shape `(nx+1, ny)` and `fy` of shape `(nx, ny+1)`, or `{"potential": [[...]]}` on the `(nx+1) x (ny+1)` vertex grid (always divergence-free). It prints the summed divergence, the signed and absolute flux on each face, and whether the two totals agree. Use `-` to read from stdin and `--json` for machine output.

Exit codes: `0` ok, `1` a recorded figure moved or a reproduction script failed, `2` bad input. Without installing, `python3 -m scripts witness` works from the checkout root, and the original `python3 scripts/<name>.py` commands are unchanged. `flux-identity reproduce` runs on a temporary copy, so it never rewrites the committed `witness.json`. The level-9 and level-10 scripts check recorded digits from long matrix-free runs; they do not redo those runs. None of this output is thrust.

## Sequel

On the combinatorial gasket, harmonic corner currents are a neutral dipole and shrink by exactly \(3/5\) at each refinement. For the spectral fractional operator, corner currents stay neutral (Theorem H) but do not inherit the factor \(3/5\) (Kept Failure H.1). The audit norm at \(\alpha=0.45\) does not tend to infinity (Theorem J); whether it tends to \(0\) or to a positive finite limit is OPEN (Conjecture J.1 is not a theorem; Kept Failures J.1–J.4 reject four ways of closing it). Dirichlet energy controls that audit norm at every \(\alpha\in(0,1)\) (Theorem L), so energy collapse would force the limit to \(0\), but collapse is not proved at \(\alpha=0.45\). For \(\alpha>\log 3/\log 5\) that same audit norm tends to \(0\) (Theorem K). For geometric weights \(1/d^2\) on the same finest-edge build, Theorem I gives exact growth by \(12/5\) so \(\|F\|\to\infty\) (Kept Failure I.1). Uniform multiplicative weights \(w=\lambda^n\) give the harmonic trichotomy of Theorem S (finite nonzero only at the resistance value \(\lambda=5/3\)); resistance weights do not freeze the fractional \(\alpha=0.45\) audit (Kept Failure S.1). Multi-scale multiplicative weights (edges at every generation) give the geometric-series trichotomy of Theorem W (finite nonzero for every \(\lambda<5/3\); resistance diverges like \((n+1)\sqrt{21}/2\), Kept Failure W.1). None of these norms is thrust. Proofs and witnesses are in [GASKET.md](GASKET.md). Run `python3 scripts/gasket_corner_current.py`, `python3 scripts/gasket_fractional_currents.py`, `python3 scripts/gasket_geometric_currents.py`, `python3 scripts/gasket_uniform_weights.py`, `python3 scripts/gasket_multiscale_weights.py`, and `python3 scripts/gasket_fractional_limit.py`.

---

<div align="center">

**REWRITE · BUILD · TRANSCEND**

**William (Brian) Ware** · [Atomic Dream Labs](https://github.com/beyond-repair)  
Governing source: [ADL-Governance](https://github.com/beyond-repair/ADL-Governance) · [Claim levels 0–5](https://github.com/beyond-repair/ADL-Governance/blob/main/docs/CLAIM_VALIDATION.md)

</div>
