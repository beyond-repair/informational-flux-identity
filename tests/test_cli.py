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


class SlowTailReadingTests(unittest.TestCase):
    def test_theorem_aw_marginal_values(self):
        """Theorem AW: tau=3/4 one-charge G=(16/9)q^2 b on every radius; dipole tau=5/4 closed form."""
        proc = subprocess.run([sys.executable, str(ROOT / "scripts" / "slow_tail_reading_flux.py")], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout[-3000:])
        self.assertIn("J_(3/4) = 4b/3 for every R>b", proc.stdout)
        self.assertIn("all Theorem AW checks passed; nothing here is thrust", proc.stdout)
        self.assertNotIn("FAIL", proc.stdout)

class MultipoleMarginalTests(unittest.TestCase):
    def test_theorem_ax_marginal_multipoles(self):
        """Theorem AX: G=Mc exactly at tau_l=(3+2l)/4; closed form for M; Kept Failure AX.1 witness."""
        proc = subprocess.run([sys.executable, str(ROOT / "scripts" / "multipole_marginal_flux.py")], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout[-3000:])
        self.assertIn("axial quadrupole 2z^2-x^2-y^2: M = diag(64, 64, -1472)/2205", proc.stdout)
        self.assertIn("all Theorem AX checks passed; nothing here is thrust", proc.stdout)
        self.assertNotIn("FAIL", proc.stdout)


class MeshMaxwellResidualTests(unittest.TestCase):
    def test_theorem_ba_mesh_residual(self):
        """Theorem BA: certified kappa_1..3(0.45) < 0, Maxwell piece sums to 0, torus control; Kept Failure BA.1."""
        proc = subprocess.run([sys.executable, str(ROOT / "scripts" / "mesh_maxwell_residual.py"), "3"], capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stdout[-3000:])
        self.assertIn("kappa_n(0.45) in [-10.7382395418199659874", proc.stdout)
        self.assertIn("torus one-sided cell stencil: D != 0 certified (Kept Failure BA.1)", proc.stdout)
        self.assertIn("all Theorem BA checks passed; nothing here is thrust", proc.stdout)
        self.assertNotIn("FAIL", proc.stdout)
        sys.path.insert(0, str(ROOT / "scripts"))
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("mesh_maxwell_residual.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("18.526339511158477727", (ROOT / doc).read_text())


class MeshMaxwellLimitTests(unittest.TestCase):
    def test_theorem_bb_kappa_decimation(self):
        """Theorem BB: exact decomposition of kappa_n, decimation products, recorded limit K(0.45) (fast float check).

        The interval certificate of K (Theorem AS sandwich) and of the cross parts runs in `reproduce` (about 2-3 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import mesh_maxwell_limit as bb
        import mesh_maxwell_residual as ba_mod

        bb.FAILS.clear()
        bb.symbolic_checks()
        bb.exact_checks()
        self.assertEqual(bb.FAILS, [])
        for n in (1, 2, 3, 4):
            d = bb.dense_level(n)
            k, ks, kc = bb.kap(d)
            self.assertLess(abs(2**n * k - ba_mod.float_kappa(n)[0]), 1e-9 * 2**n)
            self.assertLess(abs(2**n * ks - d["kself"]), 1e-9 * 2**n)
            self.assertTrue(0 < d["Y"] <= 5.0**-n * d["R"] / 2)
            self.assertTrue(0 < d["E"] <= bb.ebound(n, 0.45))
            self.assertGreater(kc, 0)                                   # Kept Failure BB.1
        K = -0.854646905120886
        q = bb.quad_level(20, "0.45")
        k20 = float(bb.kap(q)[0])
        rho = 5**0.45 / 3
        self.assertLess(abs(k20 - K + 1.08302 * rho**20), 1e-5)        # 2^-n kappa_n -> K at rate rho^n
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("mesh_maxwell_limit.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("0.8546469051208", (ROOT / doc).read_text())


class MeshCentredStencilTests(unittest.TestCase):
    def test_theorem_bc_centred_stencil(self):
        """Theorem BC: reflection law, no odd stencil on V_n, Green-Gauss rule = BA at levels 1-2; Kept Failure BC.1.

        Levels 3-4 and the E_n/Y_n recursion table run in `reproduce` (about 1.5 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import mesh_centred_stencil as bc

        bc.FAILS.clear()
        bc.symbolic_checks()
        bc.ratio_product_check()
        bc.combinatorics(levels=range(1, 5))
        rows = bc.interval_checks((1, 2))
        self.assertEqual(bc.FAILS, [])
        self.assertLess(abs(float(rows[1][0].mid) + 0.3917585210483279), 1e-13)
        self.assertLess(abs(float(rows[2][0].mid) + 0.277452251033554), 1e-13)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("mesh_centred_stencil.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Kept Failure BC.1", (ROOT / doc).read_text())
        self.assertIn("0.2038773766517652", (ROOT / "README.md").read_text())


class EndpointShiftTests(unittest.TestCase):
    def test_theorem_ay_level_shift(self):
        """Theorem AY: recorded C_gamma, nu_gamma, K_gamma reproduce the endpoint audit norm (fast float check).

        The full interval certificate runs in `reproduce` (about 3-4 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_endpoint_shift as ay
        from mpmath import mp, mpf

        mp.dps = 30
        g = mp.log(3) / mp.log(5)
        kap = 2 * mp.sin(mp.pi * g) / mp.pi
        c = mp.log(5) * mpf(ay.AU_THETA[0])
        Cg, nu, K = mpf("0.70307046242185742992904816"), mpf("1.35946811514371162883271469"), mpf("3.40236076383265507999048519")
        self.assertLess(abs(Cg / (kap * c) - nu), mpf(10) ** -25)
        for m in (10, 20):
            Jm, sm = ay.float_level(m, g)
            F = mp.sqrt(21) / 3 * sm / (kap * Jm)
            self.assertLess(abs(kap * (Jm - m * c) - Cg), 5 * 5 ** (-0.65 * m))       # R_gamma(m) - c_R m - C_gamma
            self.assertLess(abs((m + nu) * F - K), 5 * 5 ** (-0.65 * m))           # level-shift law
            self.assertLess(m * F - K, -0.2)                                       # Kept Failure AY.1: no shift is off by ~nu K/m
        self.assertLess(ay.tail_log10(60), -1e17)
        names = [s for s, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_endpoint_shift.py", names)
        for doc in ("GASKET.md", "README.md", "CLAIM_STATUS.md"):
            self.assertIn("0.70307046242185742992904816", (ROOT / doc).read_text())


class SupercriticalShiftTests(unittest.TestCase):
    def test_theorem_az_supercritical_law(self):
        """Theorem AZ: recorded A, C, nu, K at alpha=0.8 reproduce R_0.8(m) and the audit norm (fast float check).

        The full interval certificate runs in `reproduce` (about 6 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_supercritical_shift as az
        from mpmath import mp, mpf

        mp.dps = 40
        a = mpf("0.8")
        rho = mpf(5) ** a / 3
        cR, A = mpf("0.4983719948274884072103"), mpf("2.3964097039792205270396")
        C, nu, K = mpf("-1.7432571783920757239547"), mpf("-0.7274453844421553283990"), mpf("0.6878667584139920685497")
        self.assertLess(abs(cR / (rho - 1) - A), mpf(10) ** -21)
        self.assertLess(abs(C / A - nu), mpf(10) ** -21)
        for m in (10, 20):
            R, s = az.float_RS(m, a)
            F = mp.sqrt(21) / 3 * s / R
            self.assertLess(abs(R - A * rho ** m - C), 5 * rho ** m * 5 ** (-0.65 * m))   # Theorem AZ (1)
            self.assertLess(abs((rho ** m + nu) * F - K), 5 * 5 ** (-0.65 * m))           # Theorem AZ (4)
            self.assertGreater(rho ** m * (rho ** m * F - K), mpf("0.5"))               # Kept Failure AZ.1: approach from above
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_supercritical_shift.py", names)
        for doc in ("GASKET.md", "README.md", "CLAIM_STATUS.md"):
            self.assertIn("1.74325717839207572395", (ROOT / doc).read_text())


class EntryPointTests(unittest.TestCase):
    def test_python_m_scripts(self):
        proc = subprocess.run([sys.executable, "-m", "scripts", "list"], cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertIn("flux_identity.py", proc.stdout)


if __name__ == "__main__":
    unittest.main()
