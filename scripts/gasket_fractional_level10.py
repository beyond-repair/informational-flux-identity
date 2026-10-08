#!/usr/bin/env python3
"""Level-10 audit-norm record for spectral L^0.45: Kept Failure J.4.

Not a theorem about the limit. Not thrust. Does not unfreeze Stage 1.

Operator, data, and norm are those of gasket_fractional_limit.py: spectral
L^α of the combinatorial gasket Laplacian (0^α := 0), Dirichlet corner data
(1, -1/2, 0), α = 0.45, audit norm = net_flux of the combinatorial corner
currents.

Level 10 (N = 88575, 177147 edges) was computed with
scripts/gasket_fractional_matrix_free.py --level 10 --cut 0.05 --deg 150:
every eigenpair with λ ≤ 0.05 applied exactly (1581 modes, max eigen-residual
4.0e-13), λ^0.45 on the orthogonal complement as a degree-150 Chebyshev
expansion on [0.05, 8] with scalar error 9.4e-14, CG converged in 209
iterations. Interior residual of L^0.45 u 5.0e-14, corner-current sum
-4.2e-13, u in [-1/2, 1], Schur-line error 8.3e-13, centroid error 3.2e-13.
The same code and parameters reproduce the dense level-8 norm
(0.848024495086 vs 0.848024495085) and the independent level-9 matrix-free
value (cut 0.19, degree 100) to all 12 recorded digits (0.839442326601).
An independent level-10 rerun with different spectral split
(--cut 0.03 --deg 250: 941 exact modes, eigen-residual 3.7e-13, Chebyshev
scalar error 2.2e-14, interior residual 5.0e-14) gave the same
||F(10)|| = 0.833639705732 to all 12 digits.

KEPT FAILURE J.4. The contraction guess recorded with level 9,
"δ_{n+1}/δ_n < log 3 / log 5 for every n", where δ_n = 1 - ||F(n)||/||F(n-1)||,
is false: δ_10/δ_9 = 0.683038 > log 3 / log 5 = 0.682606. Equivalently,
||F(10)|| = 0.833639705732 is below the stated kill value 0.833643371744.
The margin, 3.7e-6, is about seven orders of magnitude above every recorded
numerical error. This rejects one sufficient condition for a positive limit.
It does not prove the limit is 0, and it does not prove it is positive.
Conjecture J.1 (decrease, limit at least 0.70) survives at level 10.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

F7 = 0.860817772
F8 = 0.848024495085
F9 = 0.839442326601
F10 = 0.833639705732
I10 = (0.909575490, -0.727660390, -0.181915100)  # numpy-printed digits from the run log
SUM_I10 = -4.234e-13
RES10 = 4.991e-14
EIG_RES10 = 3.982e-13
CHEB_ERR10 = 9.392e-14
Q9 = 0.948090519636
Q10 = 0.941536882637
FF9 = 0.620670938843
FF10 = 0.616380576325
F10_RERUN = 0.833639705732  # cut 0.03, degree 250
LEVEL10_KILL = 0.833643371744  # stated in gasket_fractional_level9.py
ALPHA = 0.45
ALPHA_C = math.log(3.0) / math.log(5.0)


def main() -> None:
    path = Path(__file__).with_name("gasket_fractional_limit.py")
    spec = importlib.util.spec_from_file_location("gasket_fractional_limit", path)
    gfl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gfl)
    row2 = gfl.audit_level(2, [ALPHA])
    n2 = row2["by_alpha"][ALPHA]["nF"]
    assert abs(n2 - 1.165908607) < 5e-6, n2

    delta8 = 1.0 - F8 / F7
    delta9 = 1.0 - F9 / F8
    delta10 = 1.0 - F10 / F9
    kill = F9 * (1.0 - ALPHA_C * delta9)
    assert abs(kill - LEVEL10_KILL) < 1e-12
    assert abs(F10 - F10_RERUN) < 1e-12

    # Level 9 survived the guess; level 10 breaks it.
    assert delta9 / delta8 < ALPHA_C
    contraction10 = delta10 / delta9
    assert contraction10 > ALPHA_C
    assert F10 <= LEVEL10_KILL
    margin = LEVEL10_KILL - F10
    num_err = max(RES10, EIG_RES10, CHEB_ERR10, abs(SUM_I10))
    assert margin > 1e6 * num_err, (margin, num_err)

    # Conjecture J.1 still stands at this level (not a theorem).
    assert F10 < F9 < F8 < F7
    assert F10 > 0.70

    # Schur line (Theorem H) and Theorem J cap at level 10.
    ia, ib, ic = I10
    assert abs(ia + ib + ic) < 1e-8
    assert abs(ib + 0.8 * ia) < 1e-8 and abs(ic + 0.2 * ia) < 1e-8
    assert abs(F10 - ia * math.sqrt(21.0) / 5.0) < 1e-8
    assert F10 <= 0.6 * math.sqrt(21.0)

    # Theorem L sandwich at level 10.
    w_min = ALPHA * 4.0 ** (ALPHA - 1.0)
    kappa = math.sqrt(21.0) / 5.0 * math.sqrt(2.0 / w_min)
    assert FF10 >= w_min * F10
    assert F10 <= kappa * math.sqrt(Q10)

    print(f"level2_lock {n2:.9f}")
    print(f"F9 {F9:.12f}")
    print(f"F10 {F10:.12f}")
    print(f"ratio_10_over_9 {F10 / F9:.12f}")
    print(f"delta9_over_delta8 {delta9 / delta8:.12f}")
    print(f"delta10_over_delta9 {contraction10:.12f}")
    print(f"alpha_c {ALPHA_C:.12f}")
    print(f"level10_kill {LEVEL10_KILL:.12f} margin {margin:.3e}")
    print(f"Q9 {Q9:.12f} Q10 {Q10:.12f} ratio {Q10 / Q9:.12f}")
    print("kept_failure_J4_contraction_guess_REJECTED_at_level_10", True)
    print("conjecture_J1_still_open_above_0.70", True)
    print("not_a_theorem_about_the_limit", True)
    print("not_thrust", True)


if __name__ == "__main__":
    main()
