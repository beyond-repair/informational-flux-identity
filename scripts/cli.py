"""Command-line entry point: `flux-identity` / `python -m scripts`.

Subcommands
  witness    recompute the Theorem B witness and the +17 source witness
  check      apply Theorem A to your own integer face-flux array (JSON)
  reproduce  run the repository's reproduction scripts in a scratch copy
  list       list the reproduction scripts and what each one checks

Exit codes: 0 ok, 1 a recorded figure moved or a script failed, 2 bad input.
Nothing here selects W, predicts thrust, or takes a continuum limit.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import flux_identity as fi  # noqa: E402

__version__ = "0.1.0"

EXPECTED_SIGNED = {"left": 0, "right": 11, "bottom": -2, "top": -9}
EXPECTED_ABSOLUTE = {"left": 0, "right": 349, "bottom": 2, "top": 15}

# (script, extra args, what it checks). Order matches README / GASKET.md.
SCRIPTS = [
    ("flux_identity.py", [], "Theorems A, B, E: rectangle witness 349/366, +17 source, 1D summation by parts"),
    ("hessian_flux_reduction.py", [], "Theorem M: Hessian reading reduces G to the surface Laplacian; Kept Failure M.1"),
    ("quadratic_flux_reduction.py", [], "Theorem N: quadratic reading gives G = int d_i Psi Lap Psi; Kept Failure N.1"),
    ("general_reading_flux.py", [], "Theorem O: every two-derivative reading reduces to a trace flux; Kept Failure O.1"),
    ("linear_reading_flux.py", [], "Theorem P: every linear reading of any order reduces to oint r(Lap) Psi n; Kept Failures P.1, P.2"),
    ("nonlinear_reading_flux.py", [], "Theorem Q: every shift-invariant polynomial reading is zero or surface-dependent; Kept Failure Q.1"),
    ("nonshift_reading_flux.py", [], "Theorem R: non-shift-invariant polynomial readings are still zero or surface-dependent; Kept Failure R.1"),
    ("screened_reading_flux.py", [], "Theorem T: screened polynomial readings leave only surface-dependent mass flux; Kept Failure T.1"),
    ("first_derivative_reading_flux.py", [], "Theorem U: first-derivative readings (non-polynomial allowed) are zero or surface-dependent; Kept Failure U.1"),
    ("gasket_corner_current.py", [], "harmonic gasket corner currents: neutral dipole scaling by 3/5"),
    ("gasket_fractional_currents.py", [], "Theorem H: fractional corner currents neutral; Kept Failure H.1 (no 3/5)"),
    ("gasket_geometric_currents.py", [], "Theorem I: geometric 1/d^2 currents grow by 12/5"),
    ("gasket_uniform_weights.py", [], "Theorem S: uniform multiplicative weights trichotomy; Kept Failure S.1"),
    ("gasket_fractional_limit.py", [], "Theorems J, K, L: fractional audit norm bounded, decays above log3/log5"),
    ("gasket_fractional_level9.py", [], "recorded level-9 matrix-free audit norm (checks recorded digits)"),
    ("gasket_fractional_level10.py", [], "Kept Failure J.4: log3/log5 gap contraction rejected at level 10 (recorded digits)"),
]
MATRIX_FREE = "gasket_fractional_matrix_free.py"


class InputError(ValueError):
    pass


def _int_array(value, name: str) -> np.ndarray:
    try:
        arr = np.array(value)
    except Exception as exc:  # ragged lists etc.
        raise InputError(f"{name}: not a rectangular array ({exc})") from exc
    if arr.ndim != 2 or arr.size == 0:
        raise InputError(f"{name}: expected a non-empty 2D array, got shape {arr.shape}")
    if arr.dtype == bool or not np.issubdtype(arr.dtype, np.number):
        raise InputError(f"{name}: entries must be integers")
    if not np.issubdtype(arr.dtype, np.integer):
        if not np.all(np.isfinite(arr)) or not np.all(arr == np.round(arr)):
            raise InputError(f"{name}: entries must be integers (Theorem A is stated for integer face fluxes)")
    return arr.astype(np.int64)


def load_fluxes(data: dict) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(data, dict):
        raise InputError("input must be a JSON object with 'fx' and 'fy', or 'potential'")
    if "potential" in data:
        p = _int_array(data["potential"], "potential")
        if p.shape[0] < 2 or p.shape[1] < 2:
            raise InputError("potential must be at least 2x2 vertices")
        return fi.fluxes(p)
    if "fx" not in data or "fy" not in data:
        raise InputError("input must contain 'fx' and 'fy', or 'potential'")
    fx = _int_array(data["fx"], "fx")
    fy = _int_array(data["fy"], "fy")
    nx, ny = fx.shape[0] - 1, fx.shape[1]
    if nx < 1 or fy.shape != (nx, ny + 1):
        raise InputError(
            f"shape mismatch: fx is {fx.shape}, so fy must be ({nx}, {ny + 1}); got {fy.shape}"
        )
    return fx, fy


def cmd_witness(args) -> int:
    fx, fy = fi.fluxes(fi.potential())
    div_free = fi.ledger(fx, fy)
    fx_src = fx.copy()
    fx_src[-1, 3] += 17
    sourced = fi.ledger(fx_src, fy)
    ok = (
        div_free["div_max_abs"] == 0
        and div_free["net"] == 0
        and div_free["signed"] == EXPECTED_SIGNED
        and div_free["absolute"] == EXPECTED_ABSOLUTE
        and div_free["absolute_total"] == 366
        and sourced["div_sum"] == 17
        and sourced["net"] == 17
    )
    if args.json:
        print(json.dumps({"divergence_free": div_free, "boundary_source_17": sourced, "matches_readme": ok}, indent=2))
    else:
        print("Theorem B witness, 16x16 rectangle")
        print(f"  max |div| = {div_free['div_max_abs']}   sum div = {div_free['div_sum']}   signed net = {div_free['net']}")
        print(f"  signed   left {div_free['signed']['left']}  right {div_free['signed']['right']}  bottom {div_free['signed']['bottom']}  top {div_free['signed']['top']}")
        print(f"  absolute left {div_free['absolute']['left']}  right {div_free['absolute']['right']}  bottom {div_free['absolute']['bottom']}  top {div_free['absolute']['top']}  total {div_free['absolute_total']}")
        print(f"  right-face share of absolute flux = {div_free['absolute']['right']}/{div_free['absolute_total']}")
        print("Source witness, +17 on one right-face edge")
        print(f"  sum div = {sourced['div_sum']}   signed net = {sourced['net']}")
        print("matches README:", "yes" if ok else "NO")
        print("Not thrust. Signed flux is the summed divergence; absolute flux on one face is not a net force.")
    return 0 if ok else 1


def cmd_check(args) -> int:
    try:
        if args.file == "-":
            data = json.load(sys.stdin)
        else:
            data = json.loads(Path(args.file).read_text(encoding="utf-8"))
        fx, fy = load_fluxes(data)
    except (OSError, json.JSONDecodeError, InputError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    led = fi.ledger(fx, fy)
    identity_holds = led["div_sum"] == led["net"]
    result = {
        "nx": int(fy.shape[0]),
        "ny": int(fx.shape[1]),
        **led,
        "divergence_free": led["div_max_abs"] == 0,
        "theorem_a_holds": identity_holds,
    }
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        print(f"{result['nx']}x{result['ny']} rectangle")
        print(f"  sum div = {led['div_sum']}   signed outward flux = {led['net']}   (Theorem A: equal -> {'yes' if identity_holds else 'NO'})")
        print(f"  divergence-free: {'yes' if result['divergence_free'] else 'no'} (max |div| = {led['div_max_abs']})")
        print("  signed   " + "  ".join(f"{k} {v}" for k, v in led["signed"].items()))
        print("  absolute " + "  ".join(f"{k} {v}" for k, v in led["absolute"].items()) + f"  total {led['absolute_total']}")
    return 0 if identity_holds else 1


def cmd_list(args) -> int:
    for name, extra, what in SCRIPTS:
        print(f"{name:34s} {what}")
    print(f"{MATRIX_FREE:34s} matrix-free L^0.45 audit at --level N (needs scipy; only with reproduce --matrix-free-level)")
    return 0


def cmd_reproduce(args) -> int:
    selected = SCRIPTS
    if args.only:
        names = set(args.only)
        unknown = names - {s[0] for s in SCRIPTS} - {MATRIX_FREE}
        if unknown:
            print(f"error: unknown script(s): {', '.join(sorted(unknown))}", file=sys.stderr)
            return 2
        selected = [s for s in SCRIPTS if s[0] in names]
    runs = [(name, extra) for name, extra, _ in selected]
    if args.matrix_free_level is not None:
        if not 1 <= args.matrix_free_level <= 10:
            print("error: --matrix-free-level must be between 1 and 10", file=sys.stderr)
            return 2
        runs.append((MATRIX_FREE, ["--level", str(args.matrix_free_level)]))
    failures = 0
    with tempfile.TemporaryDirectory(prefix="flux-identity-") as tmp:
        # Scratch copy: flux_identity.py writes ../witness.json next to its dir,
        # so running a copy never touches the checkout or site-packages.
        work = Path(tmp) / "scripts"
        shutil.copytree(HERE, work, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        for name, extra in runs:
            start = time.perf_counter()
            proc = subprocess.run(
                [sys.executable, str(work / name), *extra],
                cwd=tmp,
                capture_output=True,
                text=True,
            )
            secs = time.perf_counter() - start
            status = "PASS" if proc.returncode == 0 else f"FAIL (exit {proc.returncode})"
            print(f"{status:16s} {name} {' '.join(extra)}  [{secs:.1f}s]".rstrip(), flush=True)
            if args.verbose or proc.returncode != 0:
                out = (proc.stdout + proc.stderr).rstrip()
                if out:
                    print("\n".join("    " + line for line in out.splitlines()[-40:]))
            if proc.returncode != 0:
                failures += 1
    print(f"{len(runs) - failures}/{len(runs)} reproduction scripts passed")
    return 0 if failures == 0 else 1


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="flux-identity",
        description="Closed-rectangle flux identity: signed outward flux equals summed divergence. Not thrust.",
    )
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    w = sub.add_parser("witness", help="recompute the Theorem B and +17 source witnesses")
    w.add_argument("--json", action="store_true")
    w.set_defaults(func=cmd_witness)

    c = sub.add_parser("check", help="apply Theorem A to an integer face-flux array in JSON")
    c.add_argument("file", help="JSON file with {'fx': [[..]], 'fy': [[..]]} or {'potential': [[..]]}; '-' for stdin")
    c.add_argument("--json", action="store_true")
    c.set_defaults(func=cmd_check)

    r = sub.add_parser("reproduce", help="run the reproduction scripts in a scratch copy")
    r.add_argument("--only", nargs="+", metavar="SCRIPT", help="run only these scripts (names from `list`)")
    r.add_argument("--matrix-free-level", type=int, metavar="N", help="also run the matrix-free audit at level N (needs scipy)")
    r.add_argument("-v", "--verbose", action="store_true", help="show each script's output")
    r.set_defaults(func=cmd_reproduce)

    ls = sub.add_parser("list", help="list the reproduction scripts")
    ls.set_defaults(func=cmd_list)
    return ap


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
