#!/usr/bin/env python3
"""One-shot patcher: fold docs/THEOREM_V.md into README.md and update GASKET.md / counts.

Idempotent. Run from repo root: python3 scripts/apply_theorem_v_docs.py
"""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
GASKET = ROOT / "GASKET.md"
THEOREM_V = ROOT / "docs" / "THEOREM_V.md"


def main() -> None:
    v_body = THEOREM_V.read_text(encoding="utf-8").strip() + "\n"
    # Prefer LaTeX-styled section title matching README house style.
    if v_body.startswith("## Theorem V. Psi-dependent"):
        v_body = (
            "## Theorem V. $\\Psi$-dependent first-derivative isotropic readings "
            "stay Maxwell or surface-dependent\n\n"
            + v_body.split("\n", 1)[1]
        )
        # The docs file is plain-ASCII; keep it as the README body content
        # (readable; full LaTeX polish lives in the prose already written).
        v_body = THEOREM_V.read_text(encoding="utf-8").strip() + "\n"

    text = README.read_text(encoding="utf-8")
    if "## Theorem V." not in text:
        needle = (
            "**NOT THIS.** Theorem U does not compute $\\Psi_{\\mathrm{info}}$ "
            "on the frozen $0.45$ mesh"
        )
        # README uses \( \) not $ $
        needle = (
            "**NOT THIS.** Theorem U does not compute "
            "\\(\\Psi_{\\mathrm{info}}\\) on the frozen "
            "\\(0.45\\) mesh and does not decide which reading the freeze intends. "
            "It does not cover dependence on undifferentiated "
            "\\(\\Psi\\), on second or higher derivatives (beyond what Theorems Q–T "
            "already closed for polynomials), nonlocal kernels, position-dependent "
            "coefficients, pseudotensor contractions, time-dependent or radiating "
            "fields, "
            "\\(d\\neq 3\\), or any "
            "\\(\\Psi_{\\mathrm{info}}\\neq A_0\\). "
            "A nonzero "
            "\\(\\mathcal{G}\\) on a finite surface is not thrust.\n"
        )
        if needle not in text:
            raise SystemExit("README: Theorem U NOT THIS block not found")
        replacement = (
            "**NOT THIS.** Theorem U does not compute "
            "\\(\\Psi_{\\mathrm{info}}\\) on the frozen "
            "\\(0.45\\) mesh and does not decide which reading the freeze intends. "
            "Dependence on undifferentiated "
            "\\(\\Psi\\) inside the first-derivative isotropic family is taken up in "
            "Theorem V. Theorem U does not cover second or higher derivatives "
            "(beyond what Theorems Q–T already closed for polynomials), nonlocal "
            "kernels, position-dependent coefficients, pseudotensor contractions, "
            "time-dependent or radiating fields, "
            "\\(d\\neq 3\\), or any "
            "\\(\\Psi_{\\mathrm{info}}\\neq A_0\\). "
            "A nonzero "
            "\\(\\mathcal{G}\\) on a finite surface is not thrust.\n\n"
            + v_body
            + "\n"
        )
        text = text.replace(needle, replacement, 1)

    old_open = (
        "Theorem U covers every first-derivative "
        "\\(O(3)\\)-covariant reading "
        "\\(f(s)\\delta+g(s)E\\otimes E\\) with "
        "\\(f,g\\in C^1([0,\\infty))\\), polynomial or not: on-shell conserved ones "
        "are exactly Maxwell plus "
        "\\(\\kappa\\delta\\) and give "
        "\\(\\mathcal{G}=0\\); all others are surface-dependent and decay like "
        "\\(R^{-2}\\) (Kept Failure U.1); Maxwell removal leaves either "
        "\\(0\\) or that same surface-dependent remainder. Readings that depend on "
        "undifferentiated "
        "\\(\\Psi\\) or on higher derivatives *non-polynomially*, nonlocal kernels, "
        "position-dependent coefficients, time-dependent fields, and any "
        "\\(\\Psi_{\\mathrm{info}}\\) not equal to "
        "\\(A_0\\) are not covered."
    )
    new_open = (
        "Theorem U covers every first-derivative "
        "\\(O(3)\\)-covariant reading "
        "\\(f(s)\\delta+g(s)E\\otimes E\\) with "
        "\\(f,g\\in C^1([0,\\infty))\\), polynomial or not: on-shell conserved ones "
        "are exactly Maxwell plus "
        "\\(\\kappa\\delta\\) and give "
        "\\(\\mathcal{G}=0\\); all others are surface-dependent and decay like "
        "\\(R^{-2}\\) (Kept Failure U.1); Maxwell removal leaves either "
        "\\(0\\) or that same surface-dependent remainder. Theorem V closes the "
        "\\(\\Psi+\\nabla\\Psi\\) isotropic first-order gap: allowing undifferentiated "
        "\\(\\Psi\\) in "
        "\\(f(\\Psi,s)\\delta+g(\\Psi,s)E\\otimes E\\) does not enlarge the conserved "
        "class beyond Maxwell+"
        "\\(\\kappa\\delta\\) (Kept Failure V.1). Higher-derivative non-polynomial "
        "readings, nonlocal kernels, position-dependent coefficients, "
        "time-dependent fields, and any "
        "\\(\\Psi_{\\mathrm{info}}\\) not equal to "
        "\\(A_0\\) remain OPEN, as does the "
        "\\(0.45\\) mesh."
    )
    if old_open in text:
        text = text.replace(old_open, new_open, 1)

    if "psi_gradient_reading_flux.py" not in text:
        text = text.replace(
            "python3 scripts/first_derivative_reading_flux.py\n```",
            "python3 scripts/first_derivative_reading_flux.py\n"
            "python3 scripts/psi_gradient_reading_flux.py\n```",
            1,
        )

    text = text.replace("prints 15/15 passed", "prints 16/16 passed")

    README.write_text(text, encoding="utf-8")

    g = GASKET.read_text(encoding="utf-8")
    g_old = (
        "Theorem U covers every first-derivative "
        "\\(O(3)\\)-covariant reading "
        "\\(f(s)\\delta+g(s)E\\otimes E\\) with "
        "\\(C^1\\) coefficients (non-polynomial allowed): conserved ones are "
        "Maxwell plus "
        "\\(\\kappa\\delta\\) and give "
        "\\(\\mathcal{G}=0\\); all others are surface-dependent and decay like "
        "\\(R^{-2}\\) (Kept Failure U.1)."
    )
    g_new = (
        g_old
        + " Theorem V covers the "
        "\\(\\Psi+\\nabla\\Psi\\) isotropic extension "
        "\\(f(\\Psi,s)\\delta+g(\\Psi,s)E\\otimes E\\): conserved ones remain "
        "Maxwell plus "
        "\\(\\kappa\\delta\\) with "
        "\\(\\mathcal{G}=0\\); all others are surface-dependent "
        "(Kept Failure V.1)."
    )
    if "Kept Failure V.1" not in g:
        if g_old not in g:
            raise SystemExit("GASKET: Stage 2 OPEN Theorem U sentence not found")
        g = g.replace(g_old, g_new, 1)
        GASKET.write_text(g, encoding="utf-8")

    print("patched README.md and GASKET.md for Theorem V")


if __name__ == "__main__":
    main()
