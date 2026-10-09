## Theorem V. Psi-dependent first-derivative isotropic readings stay Maxwell or surface-dependent

**ASSUMPTION V (family of readings, not a change to the freeze).** Theorem U covered every first-derivative O(3)-covariant reading that depends only on s=|grad Psi|^2. Here undifferentiated Psi is allowed as well:

T^{ij} = f(Psi,s) delta^{ij} + g(Psi,s) partial^i Psi partial^j Psi,  s = |grad Psi|^2,

with f,g in C^1(R x [0, infinity)). By the representation theory of isotropic tensor functions of one scalar and one vector, this is the general local form of that class. Theorem U is the Psi-independent slice. C6 (Psi_info = A_0), static massless fields, an isolated device inside an open ball B, A_0 -> 0 at infinity, d=3, and vacuum surfaces are as in Assumption Q/U. On-shell conserved means partial_j T^{ij} = 0 whenever Delta Psi = 0. Stage 1 is not edited; no value of W, chi_vac, kappa, or m enters.

**THEOREM V.** For every such reading and every vacuum surface S bounding V:

1. **V.1 (divergence).** Writing E = grad Psi,

   partial_j T^{ij} = (f_Psi + s g_Psi) E^i + (f_s + g/2) partial_i s + g_s E^i (E · grad s) + g E^i Delta Psi.

   (The identity partial_j(E^i) E^j = (1/2) partial_i s is used in the g term.)

2. **V.2 (conservation criterion).** The reading is on-shell conserved for all harmonic Psi iff

   g_s ≡ 0,  f_s + g/2 ≡ 0,  f_Psi + s g_Psi ≡ 0

   on R x [0, infinity). Integrating: g = g(Psi) only, then f = -(g(Psi)/2) s + h(Psi), then f_Psi + s g_Psi = h' + (g'/2) s ≡ 0 for all s forces g' ≡ 0 and h' ≡ 0, so g ≡ c, h ≡ kappa, and

   T = kappa delta + c Q

   (exactly Maxwell plus a constant multiple of delta — the same conserved class as Theorem U). Psi-dependence does not enlarge the conserved class.

3. **V.3 (conserved readings vanish).** If the reading is on-shell conserved, then G(S) = 0 on every vacuum surface around an isolated device (as U.3 / N.1).

4. **V.4 (otherwise: surface-dependent).** If it is not on-shell conserved, then

   G_i(S) = -int_{outside V} partial_j T^{ij} dV.

   The value depends on which surface is drawn. After Maxwell removal from T = cQ + R, the residual is G[R]; if R is not of Maxwell+kappa delta form that residual still moves with the surface. For typical C^1 growth near infinity with A_0 ~ 1/R: pieces already treated in U (functions of s alone near 0) are O(R^{-2}) -> 0; the Taylor piece e^{-Psi^2} = 1 - Psi^2 + O(Psi^4) recovers the nonzero but surface-dependent R.4 limit of delta Psi^2; remainders such as the pseudo-Maxwell g=Psi, f=-(Psi/2)s still decay like R^{-2}.

So no reading in the family gives a nonzero G that is the same on every audit surface around the device. Removing the Maxwell piece does not create one.

**Proof of V.1.** Differentiate T^{ij} = f(Psi,s) delta^{ij} + g(Psi,s) E^i E^j. The chain rule gives f_Psi partial_j Psi delta^{ij} + f_s partial_j s delta^{ij} plus the g terms. Expand partial_j(g E^i E^j) = g_Psi E^i (E·E) + g_s (partial_j s) E^i E^j + g (partial_j E^i) E^j + g E^i Delta Psi. Use partial_j E^i E^j = (1/2) partial_i s and E·E = s to rearrange into the displayed formula.

**Proof of V.2.** On-shell, Delta Psi = 0, so the formula reduces to the three coefficient blocks in V.2. If all three vanish, the divergence vanishes identically. Conversely, harmonic polynomials realize independent combinations of those blocks (e.g. x^2-y^2 and jets with independent Psi, s, grad s), so each coefficient must vanish. The integration argument in the statement then forces g ≡ c and h ≡ kappa.

**Proof of V.3.** Immediate from V.2 and Kept Failure N.1.

**Proof of V.4.** Shell identity between S and S_R, then R -> infinity where applicable, gives the exterior formula (the same argument as Q.1 / U.4, with the R.4 caveat for even undifferentiated powers). Surface dependence is the bump argument of Theorem Q.3 applied where the on-shell divergence is nonzero (V.2).

**KEPT FAILURE V.1.** Allowing undifferentiated Psi inside the first-derivative isotropic family does not open a surface-independent nonzero G for an isolated device under C5 and static massless C6. Conserved readings remain exactly Maxwell plus kappa delta and give G = 0 on every vacuum surface (V.3). Every other reading equals minus the exterior integrated divergence and moves when the audit surface is moved (V.4); after Maxwell removal the residual is either 0 or still surface-dependent. Psi-dependence does not enlarge the conserved class. Interpretive step, labeled as such and the same as in Kept Failures O.1, P.2, Q.1, R.1, T.1, and U.1: a force on the device does not depend on where the audit surface is drawn, so a surface-dependent G cannot be read as one. This closes the Psi+grad Psi isotropic non-polynomial first-order gap left open by Theorem U. It is a failure of that route, not of Stage 2.

**Witness, computed by `scripts/psi_gradient_reading_flux.py`.** Exact (rational arithmetic): the divergence formula holds as a polynomial identity for 40/40 seeded random integer Psi of degree at most 3 and random polynomial (f,g) in (Psi,s); all 25 affine Maxwell+kappa delta readings are divergence-free on the harmonic Psi_h of Theorem Q; pseudo-Maxwell (g=Psi, f=-(Psi/2)s) has on-shell div = (s/2) E exactly; pure f=Psi (g=0) has div = E. Finite differences on harmonic and nonharmonic fields: non-polynomial pairs (e^{-Psi^2},0), (0,1/(1+Psi^2)), and (-(s/2)e^{-Psi^2}, e^{-Psi^2}) match the formula to relative 10^{-6}–10^{-8}; the third satisfies f_s+g/2=0 and g_s=0 but not f_Psi+s g_Psi=0 (|div|≈1.268 at a test point). Quadrature (128 Gauss–Legendre by 256 azimuthal nodes, 4 pi eps0 = 1) around the cluster q=(3,-1,2) of Theorems N–U, on the three vacuum spheres of Theorem O: Maxwell Q gives |G|≤9.3×10^{-14} against absolute row-x flux up to 15.712; e^{-Psi^2} delta gives (0.8677036002, 0.4738897793, -0.3661840665), (3.8753254205, 1.7138521206, -1.7879281576), (-2.7289157145, 1.6232717555, 2.2585882375); E⊗E/(1+Psi^2) gives (-0.3438719407, -0.1605731743, 0.1357384023), (-0.7373169450, -0.2964568398, 0.3352190949), (0.2564455823, -0.1551684989, -0.1982562974); pseudo-Maxwell Psi gives (-4.8241970412, -1.6431402777, 1.5606256234), (-9.6387393650, -2.7964529504, 3.8240961681), (0.7309160876, -0.5067393424, -0.5251162281). Same device, three surfaces, three values. Each non-conserved flux agrees with minus the exterior integral of the on-shell divergence to 10^{-12} relative. On T=Q+e^{-Psi^2} delta, Maxwell removal leaves exactly the e^{-Psi^2} delta flux. Under dilation of the off-centre unit sphere of Theorem R.4, e^{-Psi^2} delta → −(R.4 limit of Psi^2 delta) (nonzero, surface-dependent); pseudo-Maxwell and E⊗E/(1+Psi^2) decay like R^{-2}→0. The script exits nonzero if any of these figures moves.

**NOT THIS.** Theorem V does not compute Psi_info on the frozen 0.45 mesh and does not decide which reading the freeze intends. It does not cover dependence on second or higher derivatives non-polynomially, nonlocal kernels, position-dependent coefficients, pseudotensor contractions, time-dependent or radiating fields, d≠3, or any Psi_info ≠ A_0. A nonzero G on a finite surface is not thrust.
