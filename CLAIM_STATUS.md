# Claim status

**Sweep:** 169 (2026-10-02); gasket limit note same day, claim ceiling unchanged
**Classification:** RESEARCH
**Claim ceiling:** 1 (finite-rectangle algebraic identity plus an explicit integer witness)

| Statement | State |
|-----------|--------|
| Signed outward flux equals summed cell divergence on a finite rectangle | VERIFIED (Theorem A proof in README; `tests/test_flux_identity.py`) |
| Explicit 16x16 divergence-free array with signed net 0 and absolute right-face flux 349/366 | VERIFIED (local ledger match against `witness.json`) |
| Adding 17 to one boundary edge makes div sum and signed net both 17 | VERIFIED (script asserts; not re-run in CI this pass beyond unit tests of the identity) |
| 1D summation-by-parts identity | VERIFIED in script asserts |
| Laboratory thrust, selected W, continuum limit, gasket force | NOT CLAIMED |
| Fractional audit norm, \(\alpha=0.45\), zero vs positive | OPEN; divergence excluded by Theorem J; energy controls the norm (Theorem L; `scripts/gasket_fractional_limit.py`); log 3/log 5 gap contraction rejected at level 10 (Kept Failure J.4; `scripts/gasket_fractional_level10.py`); resistance weights do not freeze it (Kept Failure S.1; `scripts/gasket_uniform_weights.py`). Not thrust. |
| Fractional audit norm for \(\alpha>\log 3/\log 5\) | VERIFIED tends to 0 (Theorem K; same script). Not thrust. |
| Geometric \(1/d^2\) harmonic audit norm on finest-edge build | VERIFIED diverge as \((12/5)^n\sqrt{21}/2\) (Theorem I; `scripts/gasket_geometric_currents.py`) |
| Uniform multiplicative weights \(w=\lambda^n\) on finest-edge build (harmonic) | VERIFIED trichotomy: \(\to 0\) if \(\lambda<5/3\), constant \(\sqrt{21}/2\) if \(\lambda=5/3\), \(\to\infty\) if \(\lambda>5/3\) (Theorem S; `scripts/gasket_uniform_weights.py`). Resistance weights do not freeze the fractional \(\alpha=0.45\) audit (Kept Failure S.1). Not thrust. |
| Hessian reading \((\nabla\Psi)^{ij}=\partial^i\partial^j\Psi\): \(\mathcal{G}_i=\oint\Delta\Psi\,n_i\); zero for static massless \(A_0\) on a vacuum surface | VERIFIED (Theorem M, Kept Failure M.1; continuum proof in README, exact integer witness `scripts/hessian_flux_reduction.py`). Hessian reading only; 0.45 mesh solve OPEN. Not thrust. |
| Quadratic reading \((\nabla\Psi)^{ij}=\partial^i\Psi\partial^j\Psi-\tfrac12\delta^{ij}|\nabla\Psi|^2\): \(\mathcal{G}_i=\int_V\partial_i\Psi\,\Delta\Psi\); for static massless \(A_0\) it is the Maxwell stress over \(\varepsilon_0\) and an isolated device gives 0; screened case leaves only \(\tfrac{m^2}{2}\oint\Psi^2n\) | VERIFIED (Theorem N, Kept Failure N.1; continuum proof in README, exact rational and quadrature witness `scripts/quadratic_flux_reduction.py`). 0.45 mesh solve and other readings OPEN. Not thrust. |
| Every shift-invariant, \(O(d)\)-covariant, constant-coefficient two-derivative reading: \(\mathcal{G}=(\tfrac c2+d)\oint|E|^2n\) for static massless \(A_0\), zero if \(c+2d=0\), otherwise surface-dependent | VERIFIED (Theorem O, Kept Failure O.1; `scripts/general_reading_flux.py`). Not thrust. |
| Every linear, \(O(d)\)-covariant, constant-coefficient reading of any order: \(\mathcal{G}=\oint r(\Delta)\Psi\,n\); zero for shift-invariant readings with static massless \(A_0\) on a vacuum surface; the undifferentiated term gives a surface-dependent \(-\int_V E\) | VERIFIED (Theorem P, Kept Failures P.1, P.2; `scripts/linear_reading_flux.py`). Nonlinear higher-derivative readings and the 0.45 mesh solve OPEN. Not thrust. |
| Every shift-invariant, \(O(3)\)-covariant, constant-coefficient polynomial reading of any order and degree, static massless \(A_0\), isolated device: \(\mathcal{G}(S)=-\int_{\text{outside }S}\operatorname{div}T\); zero on every vacuum surface if on-shell conserved, otherwise surface-dependent | VERIFIED (Theorem Q, Lemma Q.0, Kept Failure Q.1; `scripts/nonlinear_reading_flux.py`). Non-shift-invariant nonlinear readings and the 0.45 mesh solve OPEN. Not thrust. |
| Same family without shift invariance (undifferentiated factors such as \(\delta^{ij}\Psi^2\) allowed): on-shell conserved readings have no \(\delta\Psi\) or \(\delta\Psi^2\) part and give \(\mathcal{G}=0\); all others are surface-dependent; \(\delta^{ij}\Psi^2\) has a nonzero, shape-dependent limit under dilation | VERIFIED (Theorem R, Lemma R.0, Kept Failure R.1; `scripts/nonshift_reading_flux.py`). Non-polynomial readings and the 0.45 mesh solve OPEN. Not thrust. |
| Same polynomial family under screened (Proca / Yukawa) equation \((\Delta-m^2)\Psi=-\rho/\varepsilon_0\): on-shell conserved readings give \(\mathcal{G}=0\) on every vacuum surface; all others equal a surface-dependent mass/screening remainder that decays exponentially; massless P.2/R.4 loopholes die | VERIFIED (Theorem T, Kept Failure T.1; `scripts/screened_reading_flux.py`). Non-polynomial readings and the 0.45 mesh solve OPEN. Not thrust. |
| Installable `flux-identity` CLI (`witness`, `check`, `reproduce`, `list`) and test suite | RUNS (pytest; `reproduce` executes every reproduction script; CI command `python -m unittest tests/test_flux_identity.py` now collects 3 tests). Tooling only; claim ceiling unchanged |

Parent index remains `coherence-drive` (RESEARCH). This repository does not replace it.
