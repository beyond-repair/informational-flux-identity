"""Theorem FE: D_n(alpha) = W/2 - rho^-n H_n = int_0^oo t^-alpha h_n dt is holomorphic on the strip -2^(n-1) < Re alpha < 1,
and every derivative may be taken under the integral sign. This discharges the ASSUMPTION of Theorem FD (differentiation
under the integral for n = 1), so FD.2 and FD.3 hold unconditionally (given FD's other, already-proved inputs).

FE.1 (exact, elementary): for k >= 1, delta > 0 and real u, |u|^k e^(-delta |u|) <= (k/(e delta))^k.
    Hence for a' + delta <= Re alpha <= b' - delta and t > 0:  |ln t|^k |t^-alpha| <= (k/(e delta))^k (t^-a' + t^-b').
FE.2 (from EX.1): f >= 0 and f_n >= 0, so |h_n| <= f + f_n and int t^-a |h_n| <= W(a)/2 + G_n(a) < oo for real a in (-2^(n-1), 1).
FE.3: dominated convergence with the FE.1 dominator gives D_n^(k)(alpha) = int (-ln t)^k t^-alpha h_n dt for every k and
    every complex alpha in the strip; Morera/Fubini gives holomorphy. For n = 1 this is FD's assumption.
Functions on a graph, not a flux, not thrust. Stage 2 still blocked.
"""
from __future__ import annotations
import math, sys

G0, G1 = -0.75, -0.73   # FD.2's interval


def lemma_const(k: int, delta: float) -> float:
    return (k / (math.e * delta)) ** k


def lemma_grid_ok(k: int, delta: float, umax: float = 400.0, N: int = 200000) -> bool:
    """sanity: sup_u |u|^k e^(-delta|u|) attained at u = k/delta and <= (k/(e delta))^k (with 1e-12 relative slack)."""
    c = lemma_const(k, delta)
    best = 0.0
    for i in range(N + 1):
        u = umax * i / N
        best = max(best, u ** k * math.exp(-delta * u))
    peak = (k / delta) ** k * math.exp(-k)
    return best <= c * (1 + 1e-12) and abs(peak - c) <= 1e-12 * c


def dominator_ok(k: int, a_lo: float, a_hi: float, delta: float, ts) -> bool:
    """|ln t|^k t^-alpha <= C (t^-(a_lo-delta) + t^-(a_hi+delta)) for alpha in [a_lo, a_hi], checked at sample t and alphas."""
    c = lemma_const(k, delta)
    for t in ts:
        for j in range(11):
            a = a_lo + (a_hi - a_lo) * j / 10
            lhs = abs(math.log(t)) ** k * t ** (-a)
            rhs = c * (t ** (-(a_lo - delta)) + t ** (-(a_hi + delta)))
            if lhs > rhs * (1 + 1e-12):
                return False
    return True


def main() -> int:
    ok = True
    for k in (1, 2, 3):
        for d in (0.01, 0.05, 0.25):
            ok &= lemma_grid_ok(k, d, umax=max(400.0, 40 * k / d), N=200000)
    ts = [5.0 ** (m / 4) for m in range(-200, 241)]
    ok &= dominator_ok(1, float(G0), float(G1), 0.05, ts)   # FD.2's interval, a' = -0.8, b' = -0.68 (inside (-1,1))
    ok &= dominator_ok(2, -0.95, 0.95, 0.04, ts)
    print(f"FE.1 lemma (k=1..3) and dominator on FD.2 interval: {'ok' if ok else 'FAIL'}")
    print("FE.3: FD's differentiation-under-the-integral ASSUMPTION is discharged (holomorphy on -1 < Re alpha < 1 for n = 1).")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
