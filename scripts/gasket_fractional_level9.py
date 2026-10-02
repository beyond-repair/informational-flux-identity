#!/usr/bin/env python3
"""Level 8 and level 9 audit-norm record for spectral L^0.45.

Not a theorem about the limit. Not thrust. Does not unfreeze Stage 1.

Operator and norm are those of gasket_fractional_limit.py: the spectral
power L^α of the combinatorial gasket Laplacian (0^α := 0), Dirichlet data
(1, -1/2, 0), α = 0.45, audit norm equal to the net_flux of the combinatorial
corner currents. This file does not re-prove Theorems J or K and does not
re-prove Kept Failures J.1 or J.2.

Level 8 was recomputed with that script's dense eigenpath (N = 9843):
  ||F(8)|| = 0.848024495085
  corner currents (0.925270580926, -0.740216464748, -0.185054116189),
  sum -1.14e-11
  fractional-current norm 0.627016409174
These match the previously computed 9-digit values
(0.848024495 and (0.925270581, -0.740216465, -0.185054116)).

Level 9 (N = 29526) is the same spectral operator, applied without forming
the dense eigenbasis: every eigenpair with λ ≤ 0.190556309 was applied
exactly (1147 modes, max eigen-residual 5.6e-13) and λ^0.45 on the
orthogonal complement was a degree-100 Chebyshev expansion on [λ_cut, 8]
with scalar error 4.5e-14. The same matrix-free code reproduced the dense
level-8 norm to 1e-12 (0.848024495086 versus 0.848024495085) before level 9
was accepted. The level-9 interior residual of L^0.45 u was 2.0e-14, the
corner-current sum was 4.4e-13, and u stayed in [-1/2, 1].

A contraction guess, not a theorem: δ_{n+1}/δ_n stays below log 3 / log 5,
where δ_n = 1 - ||F(n)||/||F(n-1)||. It is a sufficient condition for a
positive limit only if it holds for every n. At level 9 it still holds.
That does not prove the limit is positive, and it does not prove the limit
is 0. Conjecture J.1 (decrease, limit at least 0.70) also still stands at
this single level.
"""

from __future__ import annotations

import math
import importlib.util
from pathlib import Path

# Recorded from the dense witness at level 8 and the matrix-free spectral
# solve at level 9. Twelve digits; the level-9 norm agreed with the Schur
# identity Ia * sqrt(21) / 5 to about 5e-14.
F7 = 0.860817772
F8 = 0.848024495085
F9 = 0.839442326601
I9 = (0.915906667428, -0.732725333943, -0.183181333485)
SUM_I9 = 4.425349e-13
RES9 = 2.039e-14
# Killing value supplied with the level-8 data: ||F(9)|| at or below this
# would have pushed δ9/δ8 to log 3 / log 5 or above.
LEVEL9_KILL = 0.839421509
ALPHA_C = math.log(3.0) / math.log(5.0)


def main() -> None:
    # Live lock: same operator, level 2, before the recorded level-9 checks.
    path = Path(__file__).with_name("gasket_fractional_limit.py")
    spec = importlib.util.spec_from_file_location("gasket_fractional_limit", path)
    gfl = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gfl)
    row2 = gfl.audit_level(2, [0.45])
    n2 = row2["by_alpha"][0.45]["nF"]
    assert abs(n2 - 1.165908607) < 5e-6, n2

    assert F9 > LEVEL9_KILL
    assert F9 > 0.70
    assert F9 < F8 < F7

    ratio = F9 / F8
    delta9 = 1.0 - ratio
    delta8 = 1.0 - F8 / F7
    contraction = delta9 / delta8
    assert contraction < ALPHA_C

    # Schur line and the centroid identity used by the audit norm.
    ia, ib, ic = I9
    assert abs((ia + ib + ic)) < 1e-12
    assert abs(ib + 0.8 * ia) < 1e-12
    assert abs(ic + 0.2 * ia) < 1e-12
    assert abs(F9 - ia * math.sqrt(21.0) / 5.0) < 1e-12
    assert abs(SUM_I9) < 1e-9
    assert RES9 < 1e-8

    # Level-10 value that would falsify the contraction guess. Same rule as
    # the level-9 threshold: ||F(10)|| ≤ ||F(9)|| * (1 - α_c δ_9).
    level10_kill = F9 * (1.0 - ALPHA_C * delta9)
    assert abs(level10_kill - 0.833643371744) < 1e-12

    print(f"level2_lock {n2:.9f}")
    print(f"F8 {F8:.12f}")
    print(f"F9 {F9:.12f}")
    print(f"ratio_9_over_8 {ratio:.12f}")
    print(f"I9 {ia:.12f} {ib:.12f} {ic:.12f}")
    print(f"sum_I9 {SUM_I9:.3e}")
    print(f"delta9_over_delta8 {contraction:.12f}")
    print(f"alpha_c {ALPHA_C:.12f}")
    print("contraction_guess_at_level_9 SURVIVES")
    print(f"level10_kill_if_F10_at_or_below {level10_kill:.12f}")
    print("conjecture_J1_still_open_above_0.70", True)
    print("not_a_theorem_about_the_limit", True)
    print("not_thrust", True)


if __name__ == "__main__":
    main()
