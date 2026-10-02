# Claim status

**Sweep:** 169 (2026-10-02)
**Classification:** RESEARCH
**Claim ceiling:** 1 (finite-rectangle algebraic identity plus an explicit integer witness)

| Statement | State |
|-----------|--------|
| Signed outward flux equals summed cell divergence on a finite rectangle | VERIFIED (Theorem A proof in README; `tests/test_flux_identity.py`) |
| Explicit 16x16 divergence-free array with signed net 0 and absolute right-face flux 349/366 | VERIFIED (local ledger match against `witness.json`) |
| Adding 17 to one boundary edge makes div sum and signed net both 17 | VERIFIED (script asserts; not re-run in CI this pass beyond unit tests of the identity) |
| 1D summation-by-parts identity | VERIFIED in script asserts |
| Laboratory thrust, selected W, continuum limit, gasket force | NOT CLAIMED |
| Fractional gasket audit norm limit | OPEN (GASKET.md; not elevated) |
| Geometric \(1/d^2\) harmonic audit norm on finest-edge build | VERIFIED diverge as \((12/5)^n\sqrt{21}/2\) (Theorem I; `scripts/gasket_geometric_currents.py`) |

Parent index remains `coherence-drive` (RESEARCH). This repository does not replace it.
