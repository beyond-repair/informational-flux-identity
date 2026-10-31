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


class MeshNoncentralityLimitTests(unittest.TestCase):
    def test_theorem_bd_noncentrality_limit(self):
        """Theorem BD: e_n=(2-t_n)y_n, b_n closed form, bracket at levels 1-2, coarse c(0.45) enclosure (fast).

        The fine c(0.45), c(0.8) enclosures and the recursion table run in `reproduce` (about 4-5 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import mesh_noncentrality_limit as bd
        from mpmath import iv

        bd.FAILS.clear()
        bd.symbolic_checks()
        bd.exact_checks()
        bd.dense_checks(levels=(1, 2, 3))
        self.assertEqual(bd.FAILS, [])
        iv.dps = 30
        B = bd.B_const(iv.mpf("0.45"), iv)
        self.assertLess(abs(float(B.mid) - 115.196453544), 1e-6)
        for n in (2, 4):
            mt = float(bd.quad_mean_t(n))
            d = bd.bb.dense_level(n)
            self.assertLess(abs(d["E"] / d["Y"] - (2 - mt)), 1e-8)
            self.assertTrue(2 / (3**n + 1) < mt <= 2 / (3**n + 1) + float(B.b) / 3**n)
        c, _ = bd.c_bracket("0.45", N=40, tau0=1e-12, tau1=625, r2=1.01, M=6, k=8)
        self.assertTrue(2 < c.a < 3.7255634 < c.b < 2 + float(B.a))
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("mesh_noncentrality_limit.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Kept Failure BD.1", (ROOT / doc).read_text())
            self.assertIn("3.72198897", (ROOT / doc).read_text())


class MeshNoncentralityMonotoneTests(unittest.TestCase):
    def test_theorem_be_monotone_constant(self):
        """Theorem BE: symmetrised kernel, J_a limits, d_n decreasing in alpha at level 2, coarse S and c(0+) (fast).

        The fine S and c(0+) enclosures run in `reproduce` (about 4-5 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import mesh_noncentrality_monotone as be

        be.FAILS.clear()
        be.symbolic_checks()
        be.level_monotone(levels=(2,))
        self.assertEqual(be.FAILS, [])
        S, _ = be.S_bracket(N=40, r2=1.01)
        self.assertTrue(S.a < 2.70632481 / 3 and 2.70903118 / 3 < S.b and S.b - S.a < 0.02)
        c0, _ = be.bd.c_bracket("0", N=40, tau0=1e-12, tau1=625, r2=1.01, M=6, k=8)
        self.assertTrue(3.73 < c0.a < 5.48553067 and 5.50249013 < c0.b < 6)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("mesh_noncentrality_monotone.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Kept Failure BE.1", (ROOT / doc).read_text())
            self.assertIn("5.48553067", (ROOT / doc).read_text())


class MeshNoncentralityConvexTests(unittest.TestCase):
    def test_theorem_bf_convex_constant(self):
        """Theorem BF: sympy identities, coarse cells (1/100, level 40): point enclosures of F'(0), F''(0), F'(1)
        and F'' > 0 on the first eps-interval (fast).  The full 40-interval certificate runs in `reproduce` (about 2 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import mesh_noncentrality_convex as bf
        from mpmath import iv

        bf.FAILS.clear()
        bf.symbolic_checks()
        self.assertEqual(bf.FAILS, [])
        cells, _ = bf.build_cells(den=100, N=40)
        tail = bf.tail_consts()

        def fd(E, s0):
            N, G = bf.moments(cells, E, s0[0], s0[1], 3, den=100, tail=tail)
            return bf.F_derivs(E, N, G)

        r0 = fd(iv.mpf(0), (9, 4))
        self.assertTrue(r0[1].a < 2.70632481 and 2.70903118 < r0[1].b and r0[1].b - r0[1].a < 0.1)
        self.assertTrue(0 < r0[2].a < 1.5308 < r0[2].b)
        r1 = fd(iv.mpf(1), (14, 4))
        self.assertTrue(r1[1].a < 4.227515 and 4.341285 < r1[1].b and r1[2].a > 0)
        self.assertGreater(fd(iv.mpf([0, 0.025]), (9, 4))[2].a, 0)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("mesh_noncentrality_convex.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BF", (ROOT / doc).read_text())
            self.assertIn("0.52046", (ROOT / doc).read_text())


class FractionalExtensionLimitTests(unittest.TestCase):
    def test_theorem_bg_extension_limit(self):
        """Theorem BG: sympy identities, exact decimation of the corner-data resolvent at levels 1-2, dense u_n at
        levels 1-3 against the profile formula and the flat-limit bound, coarse W(0.45) enclosure (fast).
        The fine certificate (13485 cells, level 60) and AS's D_inf run in `reproduce` (about 3 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_extension_limit as bg

        bg.FAILS.clear()
        bg.symbolic_checks()
        bg.exact_checks(levels=(1, 2))
        bg.dense_checks(levels=(1, 2, 3), alphas=(0.45,), quad_levels=(1, 2))
        self.assertEqual(bg.FAILS, [])
        cells = bg.build_cells(N=40, tau0=1e-12, r2=1.01)
        W, _ = bg.W_bracket(cells, "0.45")
        self.assertTrue(0.5754299822 < W.b and W.a < 0.5757467849 and W.a < 0.5755882397 < W.b and W.b - W.a < 0.01)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_extension_limit.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Kept Failure BG.1", (ROOT / doc).read_text())
            self.assertIn("0.143810693", (ROOT / doc).read_text())


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


class FractionalCornerLayerTests(unittest.TestCase):
    def test_theorem_bh_corner_layer(self):
        """Theorem BH: sympy identities, exact corner-cell decomposition and neighbour recursion at levels 2-3,
        coarse arb enclosures of the layer integrals (fast). The fine certificate (742386 nodes) and AS's J_r, J_w
        run in `reproduce` (about 3 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_corner_layer as bh

        bh.FAILS.clear()
        bh.symbolic_checks()
        bh.exact_checks(levels=(2, 3))
        self.assertEqual(bh.FAILS, [])
        M, _ = bh.M_bracket("0.45", pts=bh.grid(rf=1.002, rm=1.02, rc=1.2))
        Jr = 4.0013017734819
        th = [1 - m / Jr for m in M]
        fine = [(0.4624385835869, 0.462438583587), (0.340381375, 0.340407195), (0.241594653, 0.241622693), (0.382431398, 0.382456436)]
        for x, (lo, hi) in zip(th, fine):
            self.assertTrue(x.a < lo and hi < x.b and x.b - x.a < 0.02)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_corner_layer.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BH", (ROOT / doc).read_text())
            self.assertIn("0.340381375", (ROOT / doc).read_text())


class FractionalUniformLimitTests(unittest.TestCase):
    def test_theorem_bi_uniform_limit(self):
        """Theorem BI: sympy identities, exact cell recursion and cell bound at levels 2-3, coarse arb certificate of
        sup g and the bootstrap, coarse eps_N(0.8) (upper bounds above the fine ones) and K_0.8 (fast). The fine
        certificate and the dense/float evidence run in `reproduce` (about 5 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        from fractions import Fraction

        import gasket_fractional_uniform_limit as bi
        from flint import arb

        bi.FAILS.clear()
        bi.symbolic_checks()
        bi.exact_checks(levels=(2, 3), svals=(Fraction(1, 100), Fraction(7, 5)))
        self.assertEqual(bi.FAILS, [])
        G0, _, _ = bi.sup_g(ratio="1.005")
        G = arb("0.6031")
        ok, _ = bi.bootstrap(G0, G)
        self.assertTrue(ok and G0.upper() < G.lower())
        e40 = float(bi.eps_bound("0.8", 40, G, lo="1e-30", hi="1e120", ratio="1.2").upper())
        self.assertTrue(0.00796134 <= e40 < 0.05)
        K = float(bi.K_alpha("0.8", G).upper())
        self.assertTrue(47.0 < K <= 47.109)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_uniform_limit.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BI", (ROOT / doc).read_text())
            self.assertIn("47.109", (ROOT / doc).read_text())


class FractionalCornerSharpnessTests(unittest.TestCase):
    def test_theorem_bj_corner_sharpness(self):
        """Theorem BJ: sympy identities, exact reciprocal identity at levels 0-3, coarse g-bracket (contains the fine
        one), coarse corner-constant bracket at 0.8 and a coarse N-uniform floor c_-(0.8). The fine certificate and
        the floats run in `reproduce` (about 4 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_corner_sharpness as bj
        from flint import arb

        bj.FAILS.clear()
        bj.symbolic_checks()
        bj.exact_checks(levels=(0, 1, 2, 3))
        self.assertEqual(bj.FAILS, [])
        lo0, hi0, _ = bj.g_range(ratio="1.0002")
        glo, ghi, ok = bj.g_bootstrap(lo0, hi0, arb("0.6031"))
        self.assertTrue(ok and 0.6009 < glo.lower() <= 0.600945014 and 0.600963487 <= ghi.upper() < 0.6011)
        a = arb("0.8")
        A = arb(bj.A08[0]).union(arb(bj.A08[1]))
        kap = 2 * (arb.pi() * a).sin() / (arb.pi() * A)
        Mlo, Mhi, _ = bj.M_bracket("0.8", glo, ghi, ratio="1.003")
        clo, chi = float((kap * glo * Mlo).lower()), float((kap * ghi * Mhi).upper())
        self.assertTrue(1.13 < clo <= 1.135890346 and 1.135938783 <= chi < 1.142)
        B = bj.B_lower("0.8", glo, ratio="1.02")
        self.assertTrue(float((kap * B).lower()) > 0.55)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_corner_sharpness.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BJ", (ROOT / doc).read_text())
            self.assertIn("1.135890", (ROOT / doc).read_text())


class FractionalUniformRateTests(unittest.TestCase):
    def test_theorem_bk_uniform_rate(self):
        """Theorem BK: sympy identities, exact rationals at n = 0-2, a coarse J_0.8 bound (contains the fine one) and the
        resulting C_0.8. The fine certificate and the floats run in `reproduce` (about 4 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_uniform_rate as bk
        from flint import arb

        bk.FAILS.clear()
        bk.symbolic_checks()
        bk.exact_checks(levels=(0, 1, 2))
        self.assertEqual(bk.FAILS, [])
        J, _, _, _ = bk.J_bound("0.8", arb("0.6031"), ratio="1.01")
        self.assertTrue(1.188154 <= float(J.upper()) < 1.25)
        _, _, _, C = bk.constants("0.8", J)
        self.assertTrue(6.2 < float(C.upper()) < 6.6)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_uniform_rate.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BK", (ROOT / doc).read_text())
            self.assertIn("6.2192", (ROOT / doc).read_text())


class FractionalRateSharpnessTests(unittest.TestCase):
    def test_theorem_bl_rate_sharpness(self):
        """Theorem BL: sympy identities, exact rationals for the eta_N series and monotonicity, and a coarse version of the
        N-uniform certificate (coarse cells are still rigorous, only looser). The fine certificate, the direct N = 9..13
        cases and the floats run in `reproduce` (about 3 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_rate_sharpness as bl
        from flint import arb

        bl.FAILS.clear()
        bl.symbolic_checks()
        bl.exact_checks()
        self.assertEqual(bl.FAILS, [])
        G = arb("0.6031")
        a = arb("0.8")
        rho = arb(5) ** a / 3
        Ly, _ = bl.L_lower("0.8", "0.1", ratio="1.01")
        L0, _ = bl.L_lower("0.8", "1e-10", ratio="1.01")
        Jy, _ = bl.J_upper("0.8", G, "0.1", ratio="1.01")
        Q = bl.Q_upper("0.8", ratio="1.02")
        self.assertTrue(0.32 < float(Ly.lower()) <= 0.334504 and 0.178544 <= float(Jy.upper()) < 0.19 and 19.525 <= float(Q.upper()) < 19.6)
        D0 = arb(3) ** (-a) * arb.pi() / (arb.pi() * a).sin()
        N = 14
        Dl = D0 + L0 * (1 - rho ** (-N)) / (1 - 1 / rho)
        I0 = __import__("gasket_fractional_corner_sharpness").I_closed("0.8", arb(1), G)
        Du = (G * Q * rho ** (-N) + arb(3) / 2 * (arb(3) / 5) ** N * I0) / Dl
        p0 = arb("0.1") * (arb("0.1") + 5)
        _, g = bl.one_minus_eta_series(p0, N)
        br = (g - Du) * Ly - Du * Jy
        self.assertTrue(float(br.lower()) > 0.05)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_rate_sharpness.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BL", (ROOT / doc).read_text())
            self.assertIn("0.1240", (ROOT / doc).read_text())


class FractionalRateSharpnessAlphaTests(unittest.TestCase):
    def test_theorem_bm_rate_sharpness_alpha(self):
        """Theorem BM: sympy lemmas and a coarse (still rigorous, looser) N-uniform certificate at alpha = 0.9, y0 = 1/10,
        N = 9. The fine certificates at five exponents and the floats run in `reproduce` (about 1 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_rate_sharpness_alpha as bm
        import gasket_fractional_rate_sharpness as bl
        import gasket_fractional_uniform_rate as bk
        from flint import arb

        bm.FAILS.clear()
        bm.symbolic_checks()
        self.assertEqual(bm.FAILS, [])
        G = arb("0.6031")
        al = "0.9"
        a = arb(al)
        rho = arb(5) ** a / 3
        Ly, _ = bl.L_lower(al, "0.1", ratio="1.01")
        L0, _ = bl.L_lower(al, "1e-10", ratio="1.01")
        Jy, _ = bl.J_upper(al, G, "0.1", ratio="1.01")
        Q = bl.Q_upper(al, ratio="1.02")
        _, D0, _, _ = bk.constants(al, arb(1))
        N = 9
        Dl = D0 + L0 * (1 - rho ** (-N)) / (1 - 1 / rho)
        I0 = __import__("gasket_fractional_corner_sharpness").I_closed(al, arb(1), G)
        Du = (G * Q * rho ** (-N) + arb(3) / 2 * (arb(3) / 5) ** N * I0) / Dl
        p0 = arb("0.1") * (arb("0.1") + 5)
        _, g = bl.one_minus_eta_series(p0, N)
        br = (g - Du) * Ly - Du * Jy
        self.assertTrue(float(br.lower()) > 0)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_rate_sharpness_alpha.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BM", (ROOT / doc).read_text())


class FractionalMonotoneWeightTests(unittest.TestCase):
    def test_theorem_bn_monotone_weight(self):
        """Theorem BN: sympy lemmas (monotone omega_k for alpha >= 1/17) and a coarse (still rigorous, looser) floor at
        alpha = 0.8, N = 1. The fine certificates at six exponents and the floats run in `reproduce` (about 1 min)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_monotone_weight as bn
        import gasket_fractional_rate_sharpness as bl
        from flint import arb

        bn.FAILS.clear()
        bn.symbolic_checks()
        self.assertEqual(bn.FAILS, [])
        G = arb("0.6031")
        y0 = arb("0.1")
        tau0 = y0 * (y0 + 5)
        W0 = bn.omega_hi_point("0.8", G, y0)
        self.assertTrue(float(W0.upper()) < 1.18)
        N = 1
        Gbar = arb(bl.one_minus_eta_series(arb(0.5).exp() * arb(5) ** (-N) * tau0, N)[0].upper())
        Lam = bn.Lambda_lower("0.8", W0, Gbar, lambda z: bn.series_lo(z, N), yhi="1e4", ratio="1.05")
        self.assertTrue(float(Lam.lower()) > 0)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_monotone_weight.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BN", (ROOT / doc).read_text())


class FractionalWeightThresholdTests(unittest.TestCase):
    def test_theorem_bo_weight_threshold(self):
        """Theorem BO: sympy lemmas, the local moments, and the exact-rational witness of Kept Failure BO.1 at y*.
        The full upper certificate (Sturm checks, about 15 s) runs in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_weight_threshold as bo

        bo.FAILS.clear()
        bo.symbolic_checks()
        bo.moment_checks()
        bo.lower_checks()
        self.assertEqual(bo.FAILS, [])
        self.assertTrue(bo.A_LO < bo.A_HI < bo.sp.Rational(1, 17))
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_weight_threshold.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BO", (ROOT / doc).read_text())


class FractionalDeficitSmallAlphaTests(unittest.TestCase):
    def test_theorem_bp_deficit_small_alpha(self):
        """Theorem BP: sympy lemmas and exact-rational checks of the single-crossing argument.
        The float evidence (about 30 s) runs in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_deficit_small_alpha as bp

        bp.FAILS.clear()
        bp.symbolic_checks()
        bp.exact_checks()
        self.assertEqual(bp.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_deficit_small_alpha.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BP", (ROOT / doc).read_text())


class FractionalDeficitRateSmallAlphaTests(unittest.TestCase):
    def test_theorem_bq_deficit_rate_small_alpha(self):
        """Theorem BQ: sympy lemmas and exact-rational checks of the monotone numerator.
        The arb certificates (about 2.5 min) run in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_deficit_rate_small_alpha as bq

        bq.FAILS.clear()
        bq.symbolic_checks()
        bq.exact_checks()
        self.assertEqual(bq.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_deficit_rate_small_alpha.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BQ", (ROOT / doc).read_text())


class FractionalDeficitDepthTests(unittest.TestCase):
    def test_theorem_br_deficit_depth(self):
        """Theorem BR: sympy lemmas (ceiling induction, depth step, monotone remainder, scaling) and exact rationals.
        The arb certificates (about 2 min) run in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_deficit_depth as br

        br.FAILS.clear()
        br.symbolic_checks()
        br.exact_checks()
        self.assertEqual(br.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_deficit_depth.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BR", (ROOT / doc).read_text())


class FractionalDeficitCorrectionTests(unittest.TestCase):
    def test_theorem_bs_deficit_correction(self):
        """Theorem BS: sympy lemmas (Phi factors, h monotonicity, scalings, A_N + B_N split) and exact rationals (two-level
        product identity behind the scaled BJ.1). The arb certificates (about 3 min) run in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_deficit_correction as bs

        bs.FAILS.clear()
        bs.symbolic_checks()
        bs.exact_checks()
        self.assertEqual(bs.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_deficit_correction.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BS", (ROOT / doc).read_text())


class FractionalLogPeriodicProfileTests(unittest.TestCase):
    def test_theorem_bt_log_periodic_profile(self):
        """Theorem BT: sympy lemmas (orbit ratio identity and bounds, scalings, cell-bound crossing) and exact rationals.
        The arb certificate (about 20 s) runs in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_log_periodic_profile as bt

        bt.FAILS.clear()
        bt.symbolic_checks()
        bt.exact_checks()
        self.assertEqual(bt.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_log_periodic_profile.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BT", (ROOT / doc).read_text())


class FractionalProfileFourierTests(unittest.TestCase):
    def test_theorem_bu_profile_fourier(self):
        """Theorem BU: sympy lemmas (concavity of v^a(c+dv), closed-form Mellin antiderivative, Stieltjes convexity, extended-chord
        monotonicity, psi_AU derivative bound, even-mode antisymmetry) and exact rationals. The arb certificate (about 65 s) runs
        in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_fourier as bu

        bu.FAILS.clear()
        bu.symbolic_checks()
        bu.exact_checks()
        self.assertEqual(bu.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_fourier.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BU", (ROOT / doc).read_text())


class FractionalProfileAnalyticTests(unittest.TestCase):
    def test_theorem_bv_profile_analytic(self):
        """Theorem BV: sympy lemmas (Stieltjes sector identity, contour factor, optimal strip height, arctan bound, tail sum),
        exact Gaussian-rational sector checks for r_0..r_3, and B_k >= AU's published |theta_k| for k <= 4 (arb). The AU-routine
        re-run and the extremes certificate (about 46 s) run in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_analytic as bv
        from flint import arb, ctx

        bv.FAILS.clear()
        bv.symbolic_checks()
        bv.exact_checks()
        self.assertEqual(bv.FAILS, [])
        ctx.prec = 200
        G = arb(bv.G_STAR_HI)
        x = (-2 * arb.pi() ** 2 / arb(5).log()).exp()
        for k in range(1, 5):
            mod = (arb(bv.AU_TABLE_A[k]) ** 2 + arb(bv.AU_TABLE_B[k]) ** 2).sqrt()
            self.assertTrue(bool(mod < bv.B_bound(k, G)))
            T = 4 * arb.pi() * k / arb(5).log()
            self.assertTrue(bool(bv.B_bound(k, G) <= arb(1).exp() * G * (1 + T) * x ** k))
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_analytic.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BV", (ROOT / doc).read_text())


class FractionalProfileModesTests(unittest.TestCase):
    def test_theorem_bw_profile_modes(self):
        """Theorem BW: the published tables it is checked against (AU k <= 4, BV.4 k = 5..11) obey BV.3's B_k and the
        [0.60, 0.93] normalised band (arb, fast); the deep AU-routine re-run (k <= 20, about 6 min) runs in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_modes as bw
        from flint import arb, ctx

        ctx.prec = 200
        G = arb(bw.bv.G_STAR_HI)
        x = (-2 * arb.pi() ** 2 / arb(5).log()).exp()
        table = {k: (bw.bv.AU_TABLE_A[k], bw.bv.AU_TABLE_B[k]) for k in range(1, 5)}
        table.update(bw.BV_TABLE)
        for k, (a, b) in table.items():
            mod = (arb(a) ** 2 + arb(b) ** 2).sqrt()
            self.assertTrue(bool(mod < bw.bv.B_bound(k, G)))
            c = mod / x ** k
            self.assertTrue(bool(c > arb("0.60")) and bool(c < arb("0.93")))
        self.assertEqual(bw.KMAX, 20)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_modes.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BW", (ROOT / doc).read_text())


class FractionalProfileSpectralTests(unittest.TestCase):
    def test_theorem_bx_profile_spectral(self):
        """Theorem BX: sympy lemmas (uniqueness bookkeeping, 3*5^-gamma=1, kernel identity, Mellin integral, |sin(x+iy)|^2,
        Toeplitz rank for two atoms), the BX.3 bound and the rho_k band on AU's and BV.4's published tables, and the Toeplitz
        certificate on an exact 3-atom measure (arb, fast). BW's routine re-run and T_17 > 0 (about 6.5 min) run in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_spectral as bx

        bx.FAILS.clear()
        bx.symbolic_checks()
        self.assertTrue(bx.table_checks())
        bx.toeplitz_toy_checks()
        self.assertEqual(bx.FAILS, [])
        self.assertEqual(bx.N_TOEP, 17)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_spectral.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BX", (ROOT / doc).read_text())


class FractionalProfileToeplitzTests(unittest.TestCase):
    def test_theorem_by_profile_toeplitz(self):
        """Theorem BY: parameters, the Toeplitz certificate's 3-atom sanity check (arb, fast) and doc wiring; the deep
        AU-routine re-run (k <= 21, sandwich 40/41, about 18 min) and T_18..T_20 > 0 run in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_toeplitz as by

        self.assertEqual(by.KMAX, 21)
        self.assertEqual(by.N_TOEP, (18, 19, 20))
        self.assertGreater(by.ORD, by.bw.ORD)
        by.bx.FAILS.clear()
        by.bx.toeplitz_toy_checks()
        self.assertEqual(by.bx.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_toeplitz.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BY", (ROOT / doc).read_text())


class FractionalProfileAtomsTests(unittest.TestCase):
    def test_theorem_bz_profile_atoms(self):
        """Theorem BZ: symbolic checks, the finite-level identity, the G*/min F consistency of the exact formula, one
        certified atom (the largest, lam = 3.6152567...) and doc wiring; the 128000-cell scan (311 atoms, about 2 min)
        and BZ.6 run in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_atoms as bz
        from flint import arb

        bz.FAILS.clear()
        bz.symbolic_checks()
        bz.finite_level_check()
        bz.profile_checks()
        self.assertEqual(bz.FAILS, [])
        five = arb(5)
        lo, hi = (-five ** (arb(3195) / 4000)).mid(), (-five ** (arb(3194) / 4000)).mid()
        lam, w = bz.certify_atom(lo, hi)
        self.assertTrue(lam > arb(bz.LAM_BIG[0]) and lam < arb(bz.LAM_BIG[1]))
        self.assertTrue(w > arb("0.4796") and w < arb("0.4797"))
        self.assertEqual(bz.N_SCAN, 128000)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_atoms.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem BZ", (ROOT / doc).read_text())


class FractionalProfileAtomLawTests(unittest.TestCase):
    def test_theorem_ca_profile_atom_law(self):
        """Theorem CA: symbolic checks, the accumulation certificate at lam*, the Kept Failure CA.5 witness and doc
        wiring; the identification of all 311 BZ atoms (BZ.5 scan re-run, about 2.3 min) runs in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_atom_law as ca

        ca.FAILS.clear()
        ca.symbolic_checks()
        ca.accumulation()
        ca.minus_two_witness()
        self.assertEqual(ca.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_atom_law.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CA", (ROOT / doc).read_text())


class FractionalProfileProductTests(unittest.TestCase):
    def test_theorem_cb_profile_product(self):
        """Theorem CB: symbolic checks, exact graph resolvents G_1-G_3 against the finite-level product, the product
        formula against BZ's D_of/inv_Q, the -2 branch witnesses and the closed-form atom masses (about 10 s)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_product as cb

        cb.FAILS.clear()
        cb.symbolic_checks()
        cb.graph_checks()
        cb.consistency()
        cb.branch_checks()
        cb.mass_checks()
        self.assertEqual(cb.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_product.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CB", (ROOT / doc).read_text())


class FractionalProfilePurePointTests(unittest.TestCase):
    def test_theorem_cc_profile_pure_point(self):
        """Theorem CC: symbolic transfer identity T = (2-t)/(6-t), the Cauchy-Schwarz lemma, the Nevanlinna checks and
        a level <= 6 certificate of the period sum (fast); the level <= 14 certificate (about 40 s) runs in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_pure_point as cc

        cc.FAILS.clear()
        cc.symbolic_checks()
        cc.cauchy_schwarz_exact()
        cc.nevanlinna_checks()
        found, tot, delta, mT, sums = cc.mass_identity(cc.K_FAST)
        self.assertEqual(len(found), 40)
        cc.rho_bounds(found, delta, mT, kmax=20)
        self.assertEqual(cc.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_pure_point.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CC", (ROOT / doc).read_text())


class FractionalProfileLevelCountsTests(unittest.TestCase):
    def test_theorem_cd_profile_level_counts(self):
        """Theorem CD: the endpoint orbit leaves [-5,0] in four steps, the exact recursion gives N_k = 1,0,1,3 then
        5*2^(k-4), and a direct arb recount of all 2^k preimages (k <= K_FAST) agrees; k <= 18 runs in `reproduce`."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_counts as cd

        cd.FAILS.clear()
        cd.symbolic_checks()
        orb = cd.endpoint_orbit()
        seq = cd.recursion_checks()
        self.assertEqual(seq["N"][:8], [1, 0, 1, 3, 5, 10, 20, 40])
        self.assertEqual(sum(seq["N"][:15]), 10240)
        lev = cd.preimage_sets(cd.K_FAST)
        cd.cd1_certified(lev)
        cd.direct_recount(lev, orb, cd.K_FAST)
        self.assertEqual(cd.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_counts.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CD", (ROOT / doc).read_text())


class FractionalProfileLevelRatioTests(unittest.TestCase):
    def test_theorem_ce_profile_level_ratio(self):
        """Theorem CE: sympy identities, the sigma law on O_k (k <= 6), the 5+11 level-4 branches, the exact rational
        Collatz-Wielandt certificate r in [lo, hi] (width 2e-38) and CC's level sums through level K_FAST against r."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_ratio as ce

        ce.FAILS.clear()
        ce.symbolic_checks()
        ce.cone_constants()
        ce.sigma_law(6)
        ce.branches4()
        lo, hi = ce.certificate()
        self.assertTrue(hi - lo <= 2 * ce.EPS)
        self.assertTrue(abs(float(lo) - 0.42373128696260794) < 1e-15)
        found = ce.level_sum_consistency(ce.K_FAST, lo, hi)
        ce.mass_factorisation([f for f in found if f[2] <= 5][:6])
        ce.naive_failure(10, float((lo + hi) / 2))
        self.assertEqual(ce.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_ratio.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CE", (ROOT / doc).read_text())


class FractionalProfileLevelRateTests(unittest.TestCase):
    def test_theorem_cf_profile_level_rate(self):
        """Theorem CF: sympy proof steps of Birkhoff's contraction, the s-variable cone constants, the certified diameter
        and contraction factor (arb), CF.3's bound against arb iterates for m < M_FAST, and the gap bound."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_ratio as ce
        import gasket_fractional_profile_level_rate as cf

        cf.FAILS.clear()
        cf.birkhoff_lemma()
        cf.symbolic_cone()
        consts = cf.cone_constants()
        self.assertTrue(float(consts["theta"]) < 0.5442 and float(consts["rate"]) < 0.8855)
        ce.FAILS.clear()
        lo, hi = ce.certificate()
        self.assertEqual(ce.FAILS, [])
        self.assertTrue(cf.pointwise_consistency(cf.M_FAST, consts, lo, hi))
        cf.eigen_gap(consts, hi)
        self.assertEqual(cf.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_rate.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CF", (ROOT / doc).read_text())


class FractionalProfileLevelSpectrumTests(unittest.TestCase):
    def test_theorem_cg_profile_level_spectrum(self):
        """Theorem CG: disk constants and eps_N (arb), the 81x81 Taylor matrix, Pellet counts and sign changes on its arb
        characteristic polynomial, certified resolvent bounds on the circles around r, l_2, l_3 and on |z| = 3/100
        (the |z| = 1/1000 circle runs in `reproduce`), and CE.4's [lo, hi] inside the disk around r."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_ratio as ce
        import gasket_fractional_profile_level_spectrum as cg

        cg.FAILS.clear()
        cg.symbolic_checks()
        dc = cg.disk_constants()
        eps, beta = cg.tail_constants(dc, cg.N_TRUNC)
        self.assertTrue(float(eps) < 2e-17)
        M = cg.build_matrix(cg.N_TRUNC)
        P, mids = cg.charpoly_counts(M, cg.N_TRUNC)
        deltas = cg.deltas_for(cg.N_TRUNC)
        self.assertTrue(cg.sign_changes(P, mids, deltas))
        self.assertTrue(abs(float(mids[1]) + 0.0731138175124) < 1e-12)
        self.assertTrue(abs(float(mids[1]) / float(mids[0]) + 0.1725476021) < 1e-9)
        cg.transfer(M, cg.N_TRUNC, eps, beta, mids, deltas, radii=cg.RADII_FAST)
        ce.FAILS.clear()
        lo, hi = ce.certificate()
        self.assertEqual(ce.FAILS, [])
        cg.summary(mids, deltas, lo, hi, sum((M[i, i] for i in range(cg.N_TRUNC + 1)), cg.arb(0)))
        self.assertEqual(cg.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_spectrum.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CG", (ROOT / doc).read_text())


class FractionalProfileLevelRegularityTests(unittest.TestCase):
    def test_theorem_ch_profile_level_regularity(self):
        """Theorem CH: exact CH.1 identities, Collatz-Wielandt bounds for rho(L_k), k = 1..4 (arb), the comparisons with
        CG's r, l_2, l_3, and Kept Failure CH.4's certified majorant rho(A) < 9.18 (every |mu| <= 0.1089 a C^1 eigenvalue)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_regularity as ch

        ch.FAILS.clear()
        ch.symbolic_checks()
        lams = {k: ch.collatz_wielandt(k) for k in (1, 2, 3, 4)}
        self.assertTrue(abs(lams[1] - 0.134355) < 1e-5 and abs(lams[2] - 0.0430703) < 1e-6)
        rho = ch.k_majorant()
        self.assertTrue(9.17 < float(rho.mid()) < 9.18)
        self.assertTrue(float(ch.MU_MAX) > 0.0731138175124 and float(ch.LAM[2]) < 0.0731138175124)
        ch.consequences(lams)
        self.assertTrue(5.0 < ch.k_growth_observation() < 9.18)
        self.assertEqual(ch.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_regularity.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CH", (ROOT / doc).read_text())


class FractionalProfileLevelEssentialTests(unittest.TestCase):
    def test_theorem_ci_profile_level_essential(self):
        """Theorem CI: exact CI.2 identities, two-sided Collatz-Wielandt rho(L_1) in [0.1343549, 0.1343551] (arb), the
        certified majorant rho(A) < 7.44297805367 of the sharp right inverse K (every |mu| <= 0.1343547 a C^1 eigenvalue),
        and the comparisons with CH's 0.1089, 0.1345 and CG's r."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_essential as ci

        ci.FAILS.clear()
        ci.symbolic_checks()
        lam_f, coef = ci.collatz_two_sided()
        self.assertTrue(abs(lam_f - 0.134354983) < 1e-8)
        self.assertEqual(len(coef), ci.DEG_H + 1)
        rho = ci.k_majorant(coef)
        self.assertTrue(7.4429 < float(rho.mid()) < 7.44297805367)
        self.assertTrue(float(ci.MU_MAX) * float(rho.mid()) < 1)
        self.assertTrue(float(ci.MU_CH) < float(ci.MU_MAX) < float(ci.LAM_LO) < float(ci.LAM_UP) < float(ci.LAM1_CH))
        ci.consequences()
        self.assertEqual(ci.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_essential.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CI", (ROOT / doc).read_text())


class FractionalProfileLevelDiskTests(unittest.TestCase):
    def test_theorem_cj_profile_level_disk(self):
        """Theorem CJ: Conjecture CI.6: exact (log u_i)', the certified uniform cone for L_1 (B <= 2.4121, diameter
        < 13.31), the gap-weight constant Q <= 0.8050 < 1, and the consequences spec(C^1) = {r} u {|z| <= rho(L_1)}."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_disk as cj

        cj.FAILS.clear()
        cj.symbolic_checks()
        k = cj.constants()
        self.assertTrue(float(k["B"].mid()) < 2.4121)
        self.assertTrue(float(k["Delta"].mid()) < 13.31)
        self.assertTrue(float(k["Q"].mid()) < 0.8050)
        self.assertTrue(float(k["Q"].upper()) < 1)
        cj.consequences()
        self.assertEqual(cj.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_disk.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CJ", (ROOT / doc).read_text())


class FractionalProfileLevelFredholmTests(unittest.TestCase):
    def test_theorem_ck_profile_level_fredholm(self):
        """Theorem CK: Assumption CH discharged: exact Faa di Bruno split (Lt f)^(k) = M_k f^(k) + lower order (k = 1..4),
        disjoint branch images and positive weights (arb), finite one-step Lasota-Yorke constants, re-run CI.1 / CH.2
        radii, and the now unconditional consequences (spec(C^1) = {r} u {|z| <= rho(L_1)})."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_fredholm as ck

        ck.FAILS.clear()
        coeffs = ck.symbolic_checks()
        self.assertEqual(len(coeffs), 2 * ck.KMAX)
        consts = ck.certified_constants(coeffs)
        self.assertTrue(0.28166 < float(consts[1][0].mid()) < 0.28167)   # sup L_1 1
        self.assertTrue(all(float(consts[k][0].upper()) < 0.3 for k in consts))
        ck.consequences()
        self.assertTrue(0.14 < ck.norm_identity_observation(nmax=5) < 0.16)
        self.assertEqual(ck.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_fredholm.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CK", (ROOT / doc).read_text())


class FractionalProfileLevelCalkinTests(unittest.TestCase):
    def test_theorem_cl_profile_level_calkin(self):
        """Theorem CL: r_ess(C^k) = rho(L_k) exactly (Calkin norm), two-sided arb Collatz-Wielandt for rho(L_2..4),
        C^1 Fredholm essential spectrum = closed disk, jet-zero eigenfunctions never reach the circle."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_calkin as cl

        cl.FAILS.clear()
        cl.symbolic_checks()
        floats = {k: cl.collatz_two_sided(k) for k in (2, 3, 4)}
        for k, v in floats.items():
            self.assertTrue(float(cl.CW[k][0]) < v < float(cl.CW[k][1]))
        Dw, Ds = cl.distortion_constants(n_cover=1000)
        self.assertTrue(1.68 < Dw < 1.70 and 0.72 < Ds < 0.73)
        cl.consequences()
        out, maxdl = cl.float_observations()
        self.assertTrue(out["witness_disjoint"] and out["witness_err"] < 1e-12)
        self.assertTrue(maxdl < Dw)
        self.assertEqual(cl.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_calkin.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CL", (ROOT / doc).read_text())


class FractionalProfileLevelFillTests(unittest.TestCase):
    def test_theorem_cm_profile_level_fill(self):
        """Theorem CM: every |mu| < rho(L_k) is a C^k eigenvalue of infinite multiplicity; sigma_e(C^k) = closed disk;
        Harnack constants certified, exact telescoping identities, float growth rates of the explicit series."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_fill as cm

        cm.FAILS.clear()
        cm.cl.FAILS.clear()
        cm.symbolic_checks(kmax=3)
        Ck = cm.harnack_constants(n_cover=500)
        self.assertTrue(28.0 < Ck[1] < 28.2 and 76.0 < Ck[2] < 76.5 and Ck[4] < 1e3)
        cm.consequences()
        out = cm.float_observations(ks=(2,))
        lam, res = out[2]
        for mu, ratio in res.items():
            self.assertAlmostEqual(ratio[2], abs(mu) / lam, delta=0.02 * abs(mu) / lam)
        self.assertTrue(out[("N", 2)] <= 1 + 1e-9)
        self.assertTrue(out["eig_err"] < 1e-12)
        self.assertTrue(all(out["harnack"][k] <= Ck[k] for k in (1, 2, 3, 4)))
        self.assertEqual(cm.FAILS, [])
        self.assertEqual(cm.cl.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_fill.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CM", (ROOT / doc).read_text())


class FractionalProfileLevelCoefficientTests(unittest.TestCase):
    def test_theorem_cn_profile_level_coefficient(self):
        """Theorem CN: B != 0 in CG.4. Holomorphic continuation of Lt^4 G (box checks, Krawczyk on 256 arcs), certified
        resolvent circles, CC's level sums, the enclosures of a_1, a_2, B and the explicit rate constant."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_coefficient as cn

        cn.FAILS.clear()
        cn.tail_constants()
        self.assertTrue(cn.box_checks(nx=44, ny=6))
        S = cn.circle_sup(K=256)
        self.assertTrue(0.0167 < float(S.mid()) < 0.021)
        eig = cn.spectral_layout()
        res, cg = cn.resolvent_constants()
        ec = cn.expansion_constants(S, res)
        sums = cn.level_data(13)
        a1, a2, B = cn.coefficients(sums, ec, eig)
        self.assertTrue(-8.25 < float(B.mid()) < -8.17 and float(a2.mid()) < 0 < float(a1.mid()))
        C, rem = cn.exact_rate(sums, ec, eig, a1, a2, B)
        self.assertTrue(float(C.mid()) < 15)
        self.assertEqual(cn.FAILS, [])
        self.assertEqual(cg.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_coefficient.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CN", (ROOT / doc).read_text())

    def test_theorem_co_profile_level_third(self):
        """Theorem CO: a_3 != 0 in CN.2. Continuation of Lt^4 G to a larger disk (1024 arcs), certified Taylor coefficients
        (2048 points, aliasing bound), Krawczyk eigenvectors of M_80, Riesz-projection error bounds, a_1..a_3, B, B_3."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_third as co

        co.FAILS.clear()
        co.cn.FAILS.clear()
        S = co.cn.circle_sup(K=256)
        self.assertTrue(co.big_disk(K=1024))
        S_rho = co.circle_sup_rho(K=1024)
        self.assertTrue(0.0249 < float(S_rho.mid()) < 0.03)
        f = co.taylor_coeffs(S_rho, M=2048)
        self.assertIsNotNone(f)
        a, c0 = co.projection_coeffs(f, S)
        self.assertEqual(sorted(a), ["l2", "l3", "r"])
        B, B3 = co.compare_and_conclude(a)
        a3 = a["l3"]
        self.assertFalse(a3.contains(0))
        self.assertTrue(-2.37e-5 < float(a3.mid()) < -2.35e-5 and float(a3.rad()) < 1e-7)
        self.assertTrue(-8.2108 < float(B.mid()) < -8.2106 and -1580 < float(B3.mid()) < -1560)
        self.assertEqual(co.FAILS, [])
        self.assertEqual(co.cn.FAILS, [])
        self.assertEqual(co.cg.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_third.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CO", (ROOT / doc).read_text())

    def test_theorem_cp_profile_level_threshold(self):
        """Theorem CP: explicit n_0 in CO.5. CC's level sums give sign(T_n) = (-1)^n for 4 <= n <= 12; the explicit tail
        bound E_1(12) < B_3(q_3 - 1) covers n >= 12; E_2(15)/|w|^15 < (B^2/2)(1 - q_2^2) gives the exact rate q_2^2/|q_3|."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_threshold as cp

        cp.FAILS.clear()
        cp.cn.FAILS.clear()
        k = cp.constants()
        sums = cp.cn.level_data(13)
        lo, hi, vals = cp.data_part(k, sums)
        self.assertTrue(1500 < float(lo.mid()) and float(hi.mid()) < 1750)
        P = cp.tail_coeffs(k)
        out = cp.conclude(k, P, lo, hi)
        self.assertTrue(float(out["E1"].mid()) < 270 and 1360 < float(out["c"].mid()) and float(out["C"].mid()) < 1907)
        self.assertTrue(float(out["rel"].mid()) < 5.1)
        self.assertTrue(0.7391 < float(abs(k["w"]).mid()) < 0.7392)
        self.assertEqual(cp.FAILS, [])
        self.assertEqual(cp.cn.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_threshold.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CP", (ROOT / doc).read_text())

    def test_theorem_cq_profile_level_fourth(self):
        """Theorem CQ: the fourth eigenvalue. Pellet counts at 1/1000 and 1/20000, certified resolvent circles G_4 and
        |z| = 1/20000, Krawczyk enclosure of l_4 (N = 140), a_4 != 0 by the second resolvent identity, exact rate of e_k."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_fourth as cq

        for mod in (cq, cq.cn, cq.cg, cq.co):
            mod.FAILS.clear()
        cq.cn.tail_constants()
        S = cq.cn.circle_sup(K=1024)
        self.assertTrue(cq.pellet(cq.cg.build_matrix(cq.N1)))
        M80 = cq.cg.build_matrix(cq.N1, prec=320)
        circ, eps80, beta80, dc = cq.counting(M80)
        self.assertTrue(circ["G4"][0] and circ["in"][0])
        loc = cq.localise(cq.cg.build_matrix(cq.N2, prec=320), M80, dc)
        self.assertIsNotNone(loc)
        l4 = loc["l4"]
        self.assertTrue(-2.1645e-4 < float(l4.mid()) < -2.1644e-4 and float(l4.rad()) < 1e-17)
        self.assertTrue(cq.co.big_disk(K=1024))
        S_rho = cq.co.circle_sup_rho(K=1024)
        f = cq.co.taylor_coeffs(S_rho, M=2048, Nc=cq.N2)
        self.assertIsNotNone(f)
        a4, c0, ev = cq.a4_coeff(f, S, loc)
        self.assertFalse(a4.contains(0))
        self.assertTrue(-2.21e-7 < float(a4.mid()) < -2.20e-7 and float(a4.rad()) < 1e-10)
        eps0p, k0 = cq.rate(a4, c0, ev, S, circ["in"][2], l4)
        self.assertIsNotNone(k0)
        self.assertLessEqual(k0, 30)
        for mod in (cq, cq.cn, cq.cg, cq.co):
            self.assertEqual(mod.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_fourth.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CQ", (ROOT / doc).read_text())

    def test_theorem_cs_profile_level_sharp(self):
        """Theorem CS: sharp constants for the level-sum rate. CP.2's two-sided bound makes D_n strictly decreasing to
        L = B(q_2 - 1) from D_4; the sharp C_G* in CF.3 is |log(s_5/(r s_4))| theta^(-4/5)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_sharp as cs

        for mod in (cs, cs.cp, cs.cp.cn, cs.cf):
            mod.FAILS.clear()
        self.assertEqual(cs.main(), 0)
        res = cs.RESULT
        self.assertTrue(9.6273 < float(res["L"].mid()) < 9.6275 and float(res["L"].rad()) < 1e-4)
        self.assertTrue(14.5950 < float(res["D4"].mid()) < 14.5952 and float(res["D4"].rad()) < 1e-8)
        self.assertTrue(float(res["contr"].mid()) < 0.33)
        self.assertTrue(0.0210512 < float(res["CG"].mid()) < 0.0210514 and float(res["CG"].rad()) < 1e-10)
        self.assertTrue(float(res["base"].mid()) < 0.2)
        D = res["D"]
        for n in range(4, 12):
            self.assertTrue(bool(D[n] > D[n + 1]))
        for mod in (cs, cs.cp, cs.cp.cn, cs.cf):
            self.assertEqual(mod.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_sharp.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CS", (ROOT / doc).read_text())

    def test_theorem_cz_profile_level_ratio_rate(self):
        """Theorem CZ: (q^(k) + 1/3)/(k B*^k) -> alpha* = -0.04916130431; q^(k) < -1/3 for k >= 24; |q^(k)| decreasing for k >= 23."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_ratio_rate as cz

        self.assertEqual(cz.main(), 0)
        self.assertEqual(cz.FAILS, [])
        r = cz.RESULT
        self.assertAlmostEqual(float(r["alpha_s"].mid()), -0.0491613043141717, places=14)
        self.assertAlmostEqual(float(r["gamma_s"].mid()), 0.179336821291891, places=14)
        self.assertAlmostEqual(float(r["Bst"].mid()), 0.777376027633332, places=14)
        self.assertLess(float(r["alpha_s"].rad()), 1e-25)
        self.assertEqual(r["K1"], 24)
        self.assertEqual(r["KM_mono"], 30)
        tab = r["tab"]
        self.assertLess(float(tab[24]["Gerr"].mid()), 0.416)
        self.assertLess(float(tab[40]["Gerr"].mid()), 0.0248)
        self.assertLess(float(tab[40]["w0"].mid()), 0.00683)
        self.assertLess(float(tab[60]["w0"].mid()), 0.000274)
        # first-order formula agrees with CU's certified values (CZ.6) and with CY's float collocation at k = 30, 40
        self.assertEqual(sorted(r["cu"]), list(range(24, 30)))
        self.assertTrue(r["dec"] and r["above"])
        self.assertAlmostEqual(float(r["cu"][24]["q"].mid()), -0.3357008, places=6)
        for k, fl in ((30, -6.783e-4), (40, -7.542e-5)):
            pred = float(((r["alpha_s"] * k + r["gamma_s"]) * r["Bst"] ** k).mid())
            self.assertLess(abs(pred - fl) / abs(fl), 2e-3)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_ratio_rate.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CZ", (ROOT / doc).read_text())

    def test_theorem_da_profile_level_infrate_uniform(self):
        """Theorem DA: Gamma_+-^(k) = ((2k+3)/12) kappa_-+ + O(delta_k) and CT.5's two-sided inf R_n bound for every k >= 29,
        n >= n0(k) <= 16 + k/3 (local only, about 95 s; needs python-flint)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_infrate_uniform as da

        self.assertEqual(da.main(), 0)
        self.assertEqual(da.FAILS, [])
        r = da.RESULT
        self.assertAlmostEqual(float(r["kap_m"].mid()), 0.370860305428934, places=14)
        self.assertAlmostEqual(float(r["kap_p"].mid()), 3 * float(r["kap_m"].mid()), places=13)
        self.assertLess(float(r["kap_m"].rad()), 1e-25)
        self.assertEqual(r["KDA"], 29)
        tab = r["tab"]
        self.assertEqual(sorted(tab), list(range(29, da.K1 + 1)))
        self.assertLess(float(tab[40]["delta"].mid()), 0.227)
        self.assertLess(float(tab[60]["delta"].mid()), 0.0097)
        for k, b in tab.items():
            self.assertLessEqual(b["n0"], 16 + k / 3)
        self.assertLess(float(r["rel"][40].mid()), 0.09)
        self.assertLessEqual(r["K1pc"], 53)
        self.assertLess(float(r["rel_tail"].mid()), 1e-7)
        self.assertAlmostEqual(float(r["lim_p"].mid()), 0.0824134012064, places=12)
        # float collocation (GASKET.md observation) lies inside DA.2's bands
        for k, gp, gm in ((30, 1.95186, 5.84369), (40, 2.56698, 7.69921), (60, 3.80143, 11.40425), (80, 5.03752, 15.11257)):
            b = tab[k]
            self.assertLess(abs(gp - float(b["gam_p"].mid())), float(b["delta"].mid()) + 1e-5)
            self.assertLess(abs(gm - float(b["gam_m"].mid())), float(b["delta"].mid()) + 1e-5)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_infrate_uniform.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem DA", (ROOT / doc).read_text())

    def test_theorem_db_profile_level_infrate_mid(self):
        """Theorem DB: CU.1-CU.5 certified on 3^(k+1) L_k for every k = 13..28, CT.5's two-sided inf R_n bound from n0(k);
        DB.6/DB.7 (local only, about 27 min on 8 processes; needs python-flint)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_infrate_mid as db

        self.assertEqual(db.KS, list(range(13, 29)))
        self.assertEqual(db.main(), 0)
        self.assertEqual(db.FAILS, [])
        r = db.RESULT
        tab = r["tab"]
        self.assertEqual(sorted(tab), list(range(13, 29)))
        self.assertEqual([db.n_of(k) for k in range(14, 20)], [200, 210, 220, 230, 240, 250])
        for k, b in tab.items():
            self.assertTrue(all(float(e.upper()) < 1e-12 for e in b["errs"]))  # CU.4's own threshold, unchanged
            self.assertGreater(float(b["Gp"].lower()), 0)
            self.assertGreater(float(b["Gm"].lower()), 0)
        self.assertAlmostEqual(float(tab[13]["q"].mid()), -0.349942106203, places=11)
        self.assertAlmostEqual(float(tab[13]["Gp"].mid()), 0.887577, places=5)
        self.assertAlmostEqual(float(tab[28]["Gp"].mid()), 1.82873, places=4)
        self.assertEqual(r["n0s"], dict(zip(range(13, 29), [48, 51, 55, 58, 62, 65, 68, 72, 75, 78, 82, 85, 88, 92, 95, 100])))
        self.assertEqual(r["K_2pc"], 19)
        self.assertTrue(r["inc"])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_infrate_mid.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem DB", (ROOT / doc).read_text())

    def test_theorem_dc_profile_level_birkhoff(self):
        """Theorem DC: CF.3's Birkhoff-route constant explicit: general-f K(beta) and, for the level sums,
        Theta_a(f_0, Lt f_0) <= 0.28154093, C_G^B <= 0.74544057 (local only, about 15 min; needs python-flint)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_birkhoff as dc

        for mod in (dc, dc.co, dc.cn, dc.cf, dc.co.cg):
            mod.FAILS.clear()
        self.assertEqual(dc.main(), 0)
        for mod in (dc, dc.co, dc.cn, dc.cf, dc.co.cg):
            self.assertEqual(mod.FAILS, [])
        r = dc.RESULT
        tab = r["tab"]
        self.assertLess(float(tab["0"]["K"].upper()), 1.44879)
        self.assertLess(float(tab["b*"]["K"].upper()), 4.95456)
        self.assertEqual(tab["100"]["j0"], 2)
        self.assertEqual(tab["1e4"]["j0"], 8)
        self.assertLess(float(r["b0"].upper()), 0.773338)
        self.assertLess(float(r["ThetaG"].upper()), 0.28154093)
        self.assertLess(float(r["CB"].upper()), 0.74544057)
        self.assertGreater(float(r["logR_lo"].lower()), 0.0633864)
        self.assertAlmostEqual(float(r["CGs"].mid()), 0.0210512781, places=10)
        self.assertLess(float(r["worst"].upper()), 0.05)
        self.assertTrue(35.3 < float(r["ratioC"].mid()) < 35.5)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_birkhoff.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem DC", (ROOT / doc).read_text())

    def test_theorem_dd_profile_level_cone(self):
        """Theorem DD: a sharper diameter for CF.2's cone C_(47/2) under Lt^5: V_D <= 0.73438, theta_D <= 0.42672,
        Birkhoff rate theta_D^(1/5) < 0.84339, C_D <= 1.09984, floor V_low >= 0.51383 (about 25 s; needs python-flint)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_cone as dd

        dd.FAILS.clear()
        dd.cf.FAILS.clear()
        self.assertEqual(dd.main(obs=False), 0)
        self.assertEqual(dd.FAILS, [])
        r = dd.RESULT
        self.assertLess(float(r["V"].upper()), 0.73438)
        self.assertLess(float(r["Delta"].upper()), 1.82349)
        self.assertLess(float(r["theta"].upper()), 0.42672)
        self.assertLess(float(r["rate"].upper()), 0.84339)
        self.assertLess(float(r["CD"].upper()), 1.09984)
        self.assertGreater(float(r["Vlow"].lower()), 0.51383)
        self.assertLess(float(r["worst"].upper()), 0.05)
        self.assertLess(float(r["theta"].upper()), float(r["theta_CF"].lower()))
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_cone.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem DD", (ROOT / doc).read_text())

    def test_theorem_de_profile_level_joint(self):
        """Theorem DE: the joint (ratio-of-sums) oscillation for CF.2's cone C_(47/2) under Lt^5: W_E <= 0.29040, floor
        0.29009, theta_E <= 0.15990, Birkhoff rate theta_E^(1/5) < 0.69306, |mu| < 0.29367 (about 30 s; needs python-flint)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_joint as de

        de.FAILS.clear()
        de.cf.FAILS.clear()
        de.dd.FAILS.clear()
        self.assertEqual(de.main(obs=False), 0)
        self.assertEqual(de.FAILS, [])
        r = de.RESULT
        self.assertLess(float(r["W"].upper()), 0.29040)
        self.assertLess(float(r["Delta"].upper()), 0.64513)
        self.assertLess(float(r["theta"].upper()), 0.15990)
        self.assertLess(float(r["rate"].upper()), 0.69306)
        self.assertLess(float(r["CE"].upper()), 5.28923)
        self.assertGreater(float(r["Wlow"].lower()), 0.29009)
        self.assertLess(float(r["Wlow"].upper()), float(r["W"].lower()))
        self.assertLess(float(r["worst"].upper()), 0.05)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_joint.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem DE", (ROOT / doc).read_text())

    def test_theorem_df_profile_level_block(self):
        """Theorem DF: DE's joint programme on the cone C_100 under Lt^8: W_F <= 0.04451, floor 0.04297, theta_F <= 0.02958,
        Birkhoff rate theta_F^(1/8) < 0.64398, |mu| < 0.27288, Theta_G' < 0.11477 (about 80 s; needs python-flint)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_block as df

        df.FAILS.clear()
        df.de.FAILS.clear()
        df.cf.FAILS.clear()
        df.dd.FAILS.clear()
        self.assertEqual(df.main(), 0)
        self.assertEqual(df.FAILS, [])
        r = df.RESULT
        self.assertLess(float(r["W"].upper()), 0.04451)
        self.assertLess(float(r["Delta"].upper()), 0.11835)
        self.assertLess(float(r["theta"].upper()), 0.02958)
        self.assertLess(float(r["rate"].upper()), 0.64398)
        self.assertLess(float(r["mu"].upper()), 0.27288)
        self.assertLess(float(r["TG"].upper()), 0.11477)
        self.assertLess(float(r["CF"].upper()), 14.531)
        self.assertGreater(float(r["Wlow"].lower()), 0.04297)
        self.assertLess(float(r["Wlow"].upper()), float(r["W"].lower()))
        self.assertEqual(r["worse"], [9, 10, 11, 19])
        self.assertLess(float(r["worst"].upper()), 0.2)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_block.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem DF", (ROOT / doc).read_text())

    def test_theorem_cy_profile_level_ratio_k(self):
        """Theorem CY: q^(k) = nu_2^(k)/rho(L_k) -> -1/3 (|q + 1/3| <= 4.83*0.8529^k + 0.0958*0.8352^k, k >= 24); CU's routine k = 13..23."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_ratio_k as cy
        from fractions import Fraction as Fr

        self.assertEqual(cy.main(), 0)
        self.assertEqual(cy.FAILS, [])
        r = cy.RESULT
        self.assertEqual((cy.C0, cy.R0, cy.RIN, cy.R1), (Fr(9, 5), Fr(4, 5), Fr(9, 40), Fr(1, 2)))
        self.assertEqual((r["BM"], r["GM"], r["KM"]), (Fr(8529, 10000), Fr(3572, 10000), Fr(12426, 10000)))
        self.assertEqual(r["K0"], 24)
        self.assertLess(float(r["A"].mid()), 91.4)
        self.assertLess(float(r["C2"].mid()), 4.828)
        self.assertLess(float(r["r2K0"].mid()), 1 / 3 - 9 / 40)
        self.assertLess(float(r["qbound"][40][1].mid()), 2.44e-3)
        self.assertLess(float(r["qbound"][60][1].mid()), 9.5e-5)
        self.assertLess(float(r["word_worst"].mid()), 0.8)
        # CY.5: every k = 13..23 certified by CU's routine; CY.6: |q| decreasing on 8..23 and > 1/3
        self.assertEqual(sorted(r["cu"]), list(range(13, 24)))
        self.assertTrue(r["dec"] and r["above"])
        self.assertAlmostEqual(float(r["cu"][13]["q"].mid()), -0.349942106203, places=11)
        self.assertAlmostEqual(float(r["cu"][23]["q"].mid()), -0.33622604, places=7)
        # Koenigs function: kappa o sigma_1 = -kappa/3 at a sample point, kappa'(s*) = 1 (difference quotient)
        cy.ctx.prec = 160
        s = cy.acb(2)
        k1 = cy.orbit_products(cy.sig1(s))[2]
        k0 = cy.orbit_products(s)[2]
        self.assertLess(float(abs(k1 + k0 / 3).mid()), 1e-30)
        h = cy.acb(cy.arb(10) ** -20)
        d = cy.orbit_products(cy.acb(3) / 2 + h)[2] / h
        self.assertLess(float(abs(d - 1).mid()), 1e-15)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_ratio_k.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CY", (ROOT / doc).read_text())

    def test_theorem_cx_profile_level_rate_k(self):
        """Theorem CX: (3^(k+1) rho(L_k) - 1)/(G* B*^k) -> 1, B* = B(3/2) = 0.7773760276; 3^k rho(L_k) strictly decreasing."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_rate_k as cx
        from fractions import Fraction as Fr

        self.assertEqual(cx.main(), 0)
        self.assertEqual(cx.FAILS, [])
        r = cx.RESULT
        self.assertAlmostEqual(float(r["Bstar"].mid()), 0.7773760276333319, places=14)
        self.assertAlmostEqual(float(r["Gstar"].mid()), 0.2870652417603890, places=14)
        self.assertLess(float(r["Bstar"].rad()), 1e-20)
        self.assertEqual((cx.K, cx.DELTA, cx.TAU), (40, Fr(2, 100), Fr(8, 1000)))
        self.assertLessEqual(r["LAM_B"], Fr(213, 1000))
        self.assertLessEqual(r["LAM_G"], Fr(204, 1000))
        br = r["bracket"]
        self.assertEqual(sorted(br), list(range(cx.K + 1)))
        for k, (lo, hi) in br.items():
            self.assertLess(lo, hi)
        self.assertLess(br[13][1] - br[13][0], Fr(6, 1000))
        self.assertLess(br[39][1] - br[39][0], Fr(2, 10000))
        self.assertTrue(Fr(98, 100) <= br[cx.K][0] and br[cx.K][1] <= Fr(102, 100))
        self.assertLess(float(r["worst_ratio"].mid()), 0.8136)
        self.assertLess(float(r["mono"].mid()), 0.8092)
        self.assertLess(float(r["cx3_up"].mid()), 0.9901)
        self.assertGreater(float(r["cx3_lo"].mid()), 1.0104)
        self.assertTrue(r["consistency"])
        # theta0 = (3 - sqrt5)/2 is the contraction factor towards s* = 3/2, attained at s = 5/2 (sigma_1(5/2) = sqrt5/2)
        for s, strict in ((cx.arb(2), True), (cx.arb(5) / 2, False)):
            d = float((abs(cx.sig1(s) - cx.arb(3) / 2) - r["theta0"] * abs(s - cx.arb(3) / 2)).mid())
            self.assertLess(d, -1e-3) if strict else self.assertLess(abs(d), 1e-30)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_rate_k.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CX", (ROOT / doc).read_text())

    def test_theorem_cw_profile_level_uniform_k(self):
        """Theorem CW: 0.1999*0.5999^k <= 3^(k+1) rho(L_k) - 1 <= 0.2873*0.8352^k for every k >= 0 (3^k rho(L_k) -> 1/3)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_uniform_k as cw
        from fractions import Fraction as Fr

        self.assertEqual(cw.main(), 0)
        self.assertEqual(cw.FAILS, [])
        self.assertLessEqual(4 * cw.LO * cw.LO, 5)            # X' contains X = [sqrt5/2, 5/2]
        self.assertEqual(cw.HI, Fr(5, 2))
        self.assertEqual(cw.BOUNDS, dict(B_lo=Fr(5999, 10000), B_hi=Fr(8352, 10000),
                                         G_lo=Fr(1999, 10000), G_hi=Fr(2873, 10000)))
        r = cw.RESULT
        self.assertGreater(r["cover"], 100)
        self.assertTrue(r["consistency"])
        self.assertGreater(float(r["Bmin"].mid()), 0.5999)
        self.assertLess(float(r["Bmax"].mid()), 0.8352)
        self.assertGreater(float(r["Gmin"].mid()), 0.1999)
        self.assertLess(float(r["Gmax"].mid()), 0.2873)
        # exact single-branch identity at the fixed point s* = 3/2 of sigma_1: 2 sigma_1(s*) = 3, wt_1(s*) = 1/3
        s = cw.arb(3) / 2
        self.assertLess(float(abs(2 * cw.sig1(s) - 3).mid()), 1e-30)
        self.assertLess(float(abs(cw.wt1(s) - cw.arb(1) / 3).mid()), 1e-30)
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_uniform_k.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CW", (ROOT / doc).read_text())

    def test_theorem_cv_profile_level_infrate_k5(self):
        """Theorem CV: the exact rate of inf R_n^(k) -> rho(L_k) for k = 5..12; Kept Failure CV.6: |q^(k)| not monotone."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_infrate_k5 as cv

        for mod in (cv, cv.cu, cv.cu.cr):
            mod.FAILS.clear()
        n_old, m_old = cv.cu.N, cv.cu.cr.M_POW
        self.assertEqual(cv.main(), 0)
        self.assertEqual((cv.cu.N, cv.cu.cr.M_POW), (n_old, m_old))   # CU's globals restored
        self.assertFalse(any(k >= 5 for k in cv.cu.CONF))
        expect = {  # k: (rho, nu_2, |q|, n_0)
            5: (1.4792559239634748e-3, -5.161358132104635e-4, 0.3489158332, 16),
            8: (5.27169494750854e-5, -1.89167039184835e-5, 0.3588353292, 31),
            9: (1.74329038013139e-5, -6.2429809833751e-6, 0.3581148072, 34),
            12: (6.359658135257e-7, -2.240022214049e-7, 0.3522236835, 47),
        }
        for k, (rho, nu, aq, n0) in expect.items():
            res = cv.ALL[k]
            self.assertEqual(res["k"], k)
            self.assertAlmostEqual(float(res["rho"].mid()) / rho, 1.0, delta=1e-12)
            self.assertLess(float(res["rho"].rad()), 1e-19)
            self.assertAlmostEqual(float(res["nu"].mid()) / nu, 1.0, delta=1e-12)
            self.assertLess(float(res["nu"].rad()), 1e-19)
            self.assertAlmostEqual(abs(float(res["q"].mid())), aq, delta=1e-9)
            self.assertLessEqual(res["n0"], n0)
            self.assertGreater(float(res["ev"][0].lower()), 0)
            self.assertGreater(float(res["od"][0].lower()), 0)
        aq = cv.ALL["aq"]
        self.assertTrue(all(bool(aq[k + 1] > aq[k]) for k in range(1, 8)))
        self.assertTrue(all(bool(aq[k + 1] < aq[k]) for k in range(8, 12)))
        self.assertGreater(float(cv.ALL["gap"].lower()), 7e-4)
        t3 = cv.ALL["t3"]
        self.assertTrue(all(float(t3[k].lower()) > 1 / 3 for k in t3))
        for mod in (cv, cv.cu, cv.cu.cr):
            self.assertEqual(mod.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_infrate_k5.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CV", (ROOT / doc).read_text())

    def test_theorem_cu_profile_level_infrate_k(self):
        """Theorem CU: the exact rate of inf R_n^(k) -> rho(L_k) for k = 2, 3, 4. Two certified eigenvalues of L_k on H^2(D)
        outside rho_in; rho(L_k) - inf R_n between explicit multiples of |nu_2/rho|^n for every n >= n_0(k)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_infrate_k as cu

        for mod in (cu, cu.cr):
            mod.FAILS.clear()
        self.assertEqual(cu.main(), 0)
        expect = {  # k: (rho, nu_2, q, Gamma_+, Gamma_-, c_+, c_-, n_0)
            2: (0.04307031122, -0.01247174356, -0.2895670639, 0.2546354, 0.6796090, 0.01414298, 0.03774688, 7),
            3: (0.01391925477, -0.00444017237, -0.3189949783, 0.3087743, 0.8310494, 0.00566892, 0.01525759, 9),
            4: (0.00452646909, -0.00152826859, -0.3376293000, 0.3625752, 0.9801676, 0.002195297, 0.00593466, 12),
        }
        lo_cl3 = {2: 0.0430702, 3: 0.0139192, 4: 0.0045264}
        for k, (rho, nu, q, gp, gm, cp, cm, n0) in expect.items():
            res = cu.ALL[k]
            self.assertEqual(res["k"], k)
            self.assertAlmostEqual(float(res["rho"].mid()), rho, delta=1e-10)
            self.assertLess(float(res["rho"].rad()), 1e-19)
            self.assertTrue(lo_cl3[k] < float(res["rho"].mid()) < lo_cl3[k] + 2e-7)
            self.assertAlmostEqual(float(res["nu"].mid()), nu, delta=1e-10)
            self.assertLess(float(res["nu"].rad()), 1e-19)
            self.assertAlmostEqual(float(res["q"].mid()), q, delta=1e-9)
            self.assertAlmostEqual(float(res["Gp"].mid()), gp, delta=1e-6)
            self.assertAlmostEqual(float(res["Gm"].mid()), gm, delta=1e-6)
            self.assertLess(float(res["Gp"].rad()), 1e-6)
            self.assertLess(float(res["Gm"].rad()), 1e-6)
            self.assertAlmostEqual(float(res["cp"].mid()), cp, delta=1e-7)
            self.assertAlmostEqual(float(res["cm"].mid()), cm, delta=1e-7)
            self.assertLessEqual(res["n0"], n0)
            self.assertGreater(float(res["ev"][0].lower()), 0)
            self.assertGreater(float(res["od"][0].lower()), 0)
        for mod in (cu, cu.cr):
            self.assertEqual(mod.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_infrate_k.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CU", (ROOT / doc).read_text())

    def test_theorem_ct_profile_level_infrate(self):
        """Theorem CT: the exact rate of inf R_n -> rho(L_1). Two certified eigenvalues of L_1 on H^2(D) outside 1/100;
        rho(L_1) - inf R_n between explicit multiples of |nu_2/rho|^n for every n >= 10."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_infrate as ct

        for mod in (ct, ct.cr):
            mod.FAILS.clear()
        self.assertEqual(ct.main(), 0)
        res = ct.RESULT
        rho, nu, q = res["rho"], res["nu"], res["q"]
        self.assertTrue(0.1343549 < float(rho.mid()) < 0.1343551 and float(rho.rad()) < 1e-19)
        self.assertTrue(-0.03279338 < float(nu.mid()) < -0.03279337 and float(nu.rad()) < 1e-19)
        self.assertTrue(-0.24409 < float(q.mid()) < -0.24407)
        self.assertTrue(0.20006 < float(res["Gp"].mid()) < 0.20007 and float(res["Gp"].rad()) < 1e-6)
        self.assertTrue(0.52490 < float(res["Gm"].mid()) < 0.52491 and float(res["Gm"].rad()) < 1e-6)
        self.assertTrue(0.03344 < float(res["cp"].mid()) < 0.03345)
        self.assertTrue(0.08773 < float(res["cm"].mid()) < 0.08774)
        self.assertLessEqual(res["n0"], 10)
        for mod in (ct, ct.cr):
            self.assertEqual(mod.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_infrate.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CT", (ROOT / doc).read_text())

    def test_theorem_cr_profile_level_fifth(self):
        """Theorem CR: the fifth eigenvalue. Krawczyk eigenpairs of M_300 for r, l_2..l_5, a certified spectral remainder,
        resolvent circles (CG.3's homotopy): exactly five eigenvalues in |z| > 1e-8; a_5 != 0 and the exact rate of e'_k.
        Slow (about 9 minutes)."""
        sys.path.insert(0, str(ROOT / "scripts"))
        import gasket_fractional_profile_level_fifth as cr

        for mod in (cr, cr.cn, cr.cg, cr.co):
            mod.FAILS.clear()
        self.assertEqual(cr.N, 300)
        self.assertEqual(cr.main(), 0)
        res = cr.RESULT
        l5, a5 = res["l5"], res["a5"]
        self.assertTrue(-4.5502e-8 < float(l5.mid()) < -4.5501e-8 and float(l5.rad()) < 1e-20)
        self.assertFalse(a5.contains(0))
        self.assertTrue(-2.133e-11 < float(a5.mid()) < -2.131e-11 and float(a5.rad()) < 1e-15)
        self.assertIsNotNone(res["k0"])
        self.assertLessEqual(res["k0"], 50)
        for mod in (cr, cr.cn, cr.cg, cr.co):
            self.assertEqual(mod.FAILS, [])
        names = [n for n, _, _ in __import__("cli").SCRIPTS]
        self.assertIn("gasket_fractional_profile_level_fifth.py", names)
        for doc in ("README.md", "GASKET.md", "CLAIM_STATUS.md"):
            self.assertIn("Theorem CR", (ROOT / doc).read_text())
