"""CLI and reproduction checks. Runs under pytest or `python -m unittest discover tests`."""

from __future__ import annotations

import io
import json
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import cli  # noqa: E402

try:
    import scipy  # noqa: F401

    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


def run(argv, stdin=None):
    out, err = io.StringIO(), io.StringIO()
    old_stdin = sys.stdin
    if stdin is not None:
        sys.stdin = io.StringIO(stdin)
    try:
        with redirect_stdout(out), redirect_stderr(err):
            code = cli.main(argv)
    finally:
        sys.stdin = old_stdin
    return code, out.getvalue(), err.getvalue()


class WitnessTests(unittest.TestCase):
    def test_witness_text(self):
        code, out, _ = run(["witness"])
        self.assertEqual(code, 0)
        self.assertIn("349/366", out)
        self.assertIn("matches README: yes", out)

    def test_witness_json_matches_committed_witness(self):
        code, out, _ = run(["witness", "--json"])
        self.assertEqual(code, 0)
        data = json.loads(out)
        committed = json.loads((ROOT / "witness.json").read_text(encoding="utf-8"))
        self.assertEqual(data["divergence_free"], committed["divergence_free"])
        self.assertEqual(data["boundary_source_17"], committed["boundary_source_17"])
        self.assertTrue(data["matches_readme"])


class CheckTests(unittest.TestCase):
    def _write(self, obj) -> str:
        tmp = tempfile.NamedTemporaryFile("w", suffix=".json", delete=False)
        json.dump(obj, tmp)
        tmp.close()
        self.addCleanup(Path(tmp.name).unlink)
        return tmp.name

    def test_potential_input_is_divergence_free(self):
        path = self._write({"potential": [[0, 1, 2], [3, 5, 8], [1, -4, 2]]})
        code, out, _ = run(["check", path, "--json"])
        self.assertEqual(code, 0)
        data = json.loads(out)
        self.assertTrue(data["divergence_free"])
        self.assertEqual(data["net"], 0)
        self.assertTrue(data["theorem_a_holds"])

    def test_face_flux_input_reports_source(self):
        fx = [[0, 0], [0, 0], [5, 0]]  # (nx+1, ny) with nx=2, ny=2
        fy = [[0, 0, 0], [0, 0, 0]]  # (nx, ny+1)
        code, out, _ = run(["check", "-", "--json"], stdin=json.dumps({"fx": fx, "fy": fy}))
        self.assertEqual(code, 0)
        data = json.loads(out)
        self.assertEqual(data["div_sum"], 5)
        self.assertEqual(data["net"], 5)
        self.assertEqual(data["signed"]["right"], 5)
        self.assertFalse(data["divergence_free"])

    def test_bad_inputs_exit_2(self):
        cases = [
            {"fx": [[1, 2], [3, 4.5]], "fy": [[1, 2, 3]]},  # non-integer
            {"fx": [[1, 2], [3, 4]], "fy": [[1, 2]]},  # shape mismatch
            {"nothing": 1},
            [1, 2, 3],
        ]
        for case in cases:
            code, _, err = run(["check", "-"], stdin=json.dumps(case))
            self.assertEqual(code, 2, case)
            self.assertIn("error:", err)
        code, _, _ = run(["check", "-"], stdin="not json")
        self.assertEqual(code, 2)


class ReproduceTests(unittest.TestCase):
    def test_reproduce_all_scripts_pass_without_touching_checkout(self):
        before = (ROOT / "witness.json").read_bytes()
        code, out, _ = run(["reproduce"])
        self.assertEqual(code, 0, out)
        self.assertIn(f"{len(cli.SCRIPTS)}/{len(cli.SCRIPTS)} reproduction scripts passed", out)
        self.assertEqual((ROOT / "witness.json").read_bytes(), before)

    def test_reproduce_unknown_script_exit_2(self):
        code, _, err = run(["reproduce", "--only", "nope.py"])
        self.assertEqual(code, 2)
        self.assertIn("unknown", err)

    @unittest.skipUnless(HAVE_SCIPY, "scipy not installed")
    def test_matrix_free_matches_dense_at_level_5(self):
        import gasket_fractional_limit as gfl

        dense = gfl.audit_level(5, [0.45])["by_alpha"][0.45]["nF"]
        proc = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / cli.MATRIX_FREE), "--level", "5"],
            capture_output=True,
            text=True,
            check=True,
        )
        record = [line for line in proc.stdout.splitlines() if line.startswith("RECORD")][0]
        f_value = float(record.split()[1].split("=")[1])
        self.assertAlmostEqual(f_value, dense, places=9)


class EntryPointTests(unittest.TestCase):
    def test_python_m_scripts(self):
        proc = subprocess.run([sys.executable, "-m", "scripts", "list"], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("flux_identity.py", proc.stdout)


if __name__ == "__main__":
    unittest.main()
