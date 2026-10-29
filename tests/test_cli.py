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
