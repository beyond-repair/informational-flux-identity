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

## What was tried and does not follow

- **FAILED as a derivation of thrust.** Nothing above produces a nonzero \(\Delta F\) for a divergence-free informational tensor. The mathematics does not fail at a hidden algebraic step. The net is zero because the sum telescopes.
- **FAILED as a derivation of \(0.08\) or \(0.23\).** Those numerals never enter.
- **FAILED as a reading of the informational fork protocol.** The inequality \(T_{\mathrm{Red}} > T_{\mathrm{CIS}}\times 10^3\) is not a stress tensor and supplies no \(\mathrm{div}\,F\). It is not used.
- **NOT DERIVED.** Einstein gravity from entanglement. The essay in [-Entanglement-and-Emergence](https://github.com/beyond-repair/-Entanglement-and-Emergence) remains an essay. The only entanglement-facing consequence here is conditional: if an entanglement stress is divergenceless, Theorem A says its closed-surface signed flux is zero. Whether any concrete entanglement stress is divergenceless is **OPEN**.
- **OPEN.** Whether \((\nabla\Psi_{\mathrm{info}})^{ij}\) on the frozen \(0.45\) mesh has \(\sum\mathrm{div} \neq 0\) after the Maxwell piece is removed. That is a Stage 2 computation, not a symbolic gap. Under the Hessian reading, Theorem M reduces it to whether \(\Delta\Psi_{\mathrm{info}}\neq 0\) on the enclosing surface; with static massless \(\Psi=A_0\) on a vacuum surface it is \(0\) (Kept Failure M.1). Under the quadratic reading, Theorem N makes \(\mathcal{G}\) the volume integral of \(\partial_i\Psi\,\Delta\Psi\); with static massless \(\Psi=A_0\) that tensor is the Maxwell stress over \(\varepsilon_0\), so nothing remains after the Maxwell piece is removed, and an isolated device gives \(\mathcal{G}=0\) (Kept Failure N.1). Theorem O covers every shift-invariant, rotation-covariant, constant-coefficient reading with two derivatives: with static massless \(\Psi=A_0\), \(\mathcal{G}\) is either \(0\) or a trace flux \((\tfrac c2+d)\oint|E|^2n\) that moves with the surface and sums to \(0\) over all space (Kept Failure O.1). Theorem P covers every linear, rotation-covariant, constant-coefficient reading of any order: shift-invariant ones give \(\mathcal{G}=0\) on a vacuum surface (Kept Failure P.1), and the one undifferentiated term \(\delta^{ij}\Psi\) gives \(-\int_V E\), which depends on the surface (Kept Failure P.2). Nonlinear readings with more than two derivatives, position-dependent coefficients, and any \(\Psi_{\mathrm{info}}\) not equal to \(A_0\) are not covered.
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
flux-identity reproduce         # runs every reproduction script below in a scratch copy; prints 9/9 passed
flux-identity reproduce --matrix-free-level 5   # also runs the SciPy matrix-free audit at level 5
flux-identity list              # what each script checks
python -m pytest                # identity, CLI, reproduction, and matrix-free vs dense checks
```

`flux-identity check FILE.json` applies Theorem A to your own integer face-flux array. Give either `{"fx": [[...]], "fy": [[...]]}` with `fx` of shape `(nx+1, ny)` and `fy` of shape `(nx, ny+1)`, or `{"potential": [[...]]}` on the `(nx+1) x (ny+1)` vertex grid (always divergence-free). It prints the summed divergence, the signed and absolute flux on each face, and whether the two totals agree. Use `-` to read from stdin and `--json` for machine output.

Exit codes: `0` ok, `1` a recorded figure moved or a reproduction script failed, `2` bad input. Without installing, `python3 -m scripts witness` works from the checkout root, and the original `python3 scripts/<name>.py` commands are unchanged. `flux-identity reproduce` runs on a temporary copy, so it never rewrites the committed `witness.json`. The level-9 and level-10 scripts check recorded digits from long matrix-free runs; they do not redo those runs. None of this output is thrust.

## Sequel

On the combinatorial gasket, harmonic corner currents are a neutral dipole and shrink by exactly \(3/5\) at each refinement. For the spectral fractional operator, corner currents stay neutral (Theorem H) but do not inherit the factor \(3/5\) (Kept Failure H.1). The audit norm at \(\alpha=0.45\) does not tend to infinity (Theorem J); whether it tends to \(0\) or to a positive finite limit is OPEN (Conjecture J.1 is not a theorem; Kept Failures J.1–J.3 reject three ways of closing it). Dirichlet energy controls that audit norm at every \(\alpha\in(0,1)\) (Theorem L), so energy collapse would force the limit to \(0\), but collapse is not proved at \(\alpha=0.45\). For \(\alpha>\log 3/\log 5\) that same audit norm tends to \(0\) (Theorem K). For geometric weights \(1/d^2\) on the same finest-edge build, Theorem I gives exact growth by \(12/5\) so \(\|F\|\to\infty\) (Kept Failure I.1). None of these norms is thrust. Proofs and witnesses are in [GASKET.md](GASKET.md). Run `python3 scripts/gasket_corner_current.py`, `python3 scripts/gasket_fractional_currents.py`, `python3 scripts/gasket_geometric_currents.py`, and `python3 scripts/gasket_fractional_limit.py`.

---

<div align="center">

**REWRITE · BUILD · TRANSCEND**

**William (Brian) Ware** · [Atomic Dream Labs](https://github.com/beyond-repair)  
Governing source: [ADL-Governance](https://github.com/beyond-repair/ADL-Governance) · [Claim levels 0–5](https://github.com/beyond-repair/ADL-Governance/blob/main/docs/CLAIM_VALIDATION.md)

</div>
