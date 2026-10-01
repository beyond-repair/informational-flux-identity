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

## What was tried and does not follow

- **FAILED as a derivation of thrust.** Nothing above produces a nonzero \(\Delta F\) for a divergence-free informational tensor. The mathematics does not fail at a hidden algebraic step. The net is zero because the sum telescopes.
- **FAILED as a derivation of \(0.08\) or \(0.23\).** Those numerals never enter.
- **FAILED as a reading of the informational fork protocol.** The inequality \(T_{\mathrm{Red}} > T_{\mathrm{CIS}}\times 10^3\) is not a stress tensor and supplies no \(\mathrm{div}\,F\). It is not used.
- **NOT DERIVED.** Einstein gravity from entanglement. The essay in [-Entanglement-and-Emergence](https://github.com/beyond-repair/-Entanglement-and-Emergence) remains an essay. The only entanglement-facing consequence here is conditional: if an entanglement stress is divergenceless, Theorem A says its closed-surface signed flux is zero. Whether any concrete entanglement stress is divergenceless is **OPEN**.
- **OPEN.** Whether \((\nabla\Psi_{\mathrm{info}})^{ij}\) on the frozen \(0.45\) mesh has \(\sum\mathrm{div} \neq 0\) after the Maxwell piece is removed. That is a Stage 2 computation, not a symbolic gap.
- **OPEN.** Distributional flux on the infinite gasket. Theorem A is the finite rectangle. It does not pass to a limit that has not been constructed.
- **OPEN.** The renormalized local propagator for \(W(x)\), already open in the Ware derivation ledger. Theorem E constrains the force ledger of a variable weight. It does not construct \(Z_{\mathrm{ren}}\).

## Not a consequence

No laboratory force, no energy-extraction law, no selected \(W\), no confirmation of a pinch percentage, and no change to the Stage 1 freeze. The symbolic chain may stay written as \(\Delta F = W\chi\mathcal{G}\). This note fixes what \(\mathcal{G}\) has to equal.

## Reproduce

```bash
python3 scripts/flux_identity.py
```

The script writes `witness.json` and exits nonzero if the integers quoted above move.

## Sequel

On the combinatorial gasket, harmonic corner currents are a neutral dipole and shrink by exactly \(3/5\) at each refinement. For the spectral fractional operator \(L^{0.45}\), corner currents stay neutral (Theorem H) but the refinement ratios are not \(3/5\) (Kept Failure H.1); whether the fractional audit norm has a nonzero limit remains OPEN. Proofs and witnesses are in [GASKET.md](GASKET.md). Run `python3 scripts/gasket_corner_current.py` and `python3 scripts/gasket_fractional_currents.py`.

---

<div align="center">

**REWRITE · BUILD · TRANSCEND**

**William (Brian) Ware** · [Atomic Dream Labs](https://github.com/beyond-repair)  
Governing source: [ADL-Governance](https://github.com/beyond-repair/ADL-Governance) · [Claim levels 0–5](https://github.com/beyond-repair/ADL-Governance/blob/main/docs/CLAIM_VALIDATION.md)

</div>
