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
    ("psi_gradient_reading_flux.py", [], "Theorem V: Psi-dependent first-derivative isotropic readings stay Maxwell or surface-dependent; Kept Failure V.1"),
    ("position_dependent_reading_flux.py", [], "Theorem Y: position-dependent first-derivative isotropic readings stay Maxwell or surface-dependent; Kept Failure Y.1"),
    ("hemitropic_reading_flux.py", [], "Theorem Z: pseudotensor first-derivative readings stay Maxwell, a boundary-only curl, or surface-dependent; Kept Failure Z.1"),
    ("nonlocal_reading_flux.py", [], "Theorem AV: nonlocal (convolution / fractional) Maxwell-form readings: conserved on vacuum iff K=c delta; otherwise surface-dependent, ->0 under dilation, total 0 (Newton); exact G=-(64/3pi)q^2R^3b/(R^2-b^2)^4 for (-Lap)^(1/2); Kept Failure AV.1"),
    ("slow_tail_reading_flux.py", [], "Theorem AW: slow-tail (-Lap)^(-tau) readings: literal reading diverges for tau>=(1+l)/2; surface-independent G is 0 for every tau, surface-dependent for non-integer tau; marginal tau=3/4 gives (16/9)Q^2(c-x_Q), dipole tau=5/4 (16/75)(|p|^2c-8(p.c)p); Kept Failure AW.1"),
    ("multipole_marginal_flux.py", [], "Theorem AX: proves Conjecture AW.4: marginal tau_l=(3+2l)/4 gives G=Mc exactly on every enclosing sphere, M closed form for all l (quadrupole M=(16/81)(trA)^2I+(32/735)tr(A0^2)I-(512/2205)A0^2); literal reading never exists there; Kept Failure AX.1: P*+d_x^5 gives a centre- and radius-free nonzero sphere value, cubes differ"),
    ("mesh_maxwell_residual.py", [], "Theorem BA: Maxwell-form reading on the 0.45 mesh: action-reaction piece integrates to 0 exactly; residual D=kappa_n Q(d), certified kappa_1..4(0.45)<0 (kappa_4=-18.5263395111584777...), kappa_n(1)=-9((18/25)^n-(6/25)^n); Kept Failure BA.1: torus control (one-sided stencil nonzero, centred 0), self-force dominated, grows like 2^n"),
    ("mesh_maxwell_limit.py", [], "Theorem BB: proves Conjecture BA.2: kappa_n=2^n(4/R_n^2)(3Y_n-sigma_n), 0<Y_n<=5^-n R_n/2 (decimation product), so 2^-n kappa_n(0.45) -> -4 sigma_inf/R_inf^2 in [-0.854646905120924,-0.854646905120848]; limit 0 for alpha>=log3/log5; cross (non-self) part -> 0; Kept Failure BB.1: cross part >0 certified at levels 1-4"),
    ("mesh_centred_stencil.py", [], "Theorem BC: centred stencil on the 0.45 gasket: S3 reflection law H_kc=R_ck H_ck makes the Maxwell piece the central pair force, kappa_n=-(4/R^2)(s+2 sqrt3 h_perp); no vertex admits an odd stencil (D3, no half-turn); Green-Gauss centred rule = one-sided at corners so kappa unchanged (certified levels 1-4); Kept Failure BC.1: certified non-centrality tan theta_4=-0.20387737665...; Conjecture BC.2 E_n/Y_n->2"),
    ("mesh_noncentrality_limit.py", [], "Theorem BD: proves Conjecture BC.2: exact e_n=(2-t_n)y_n with t_n increasing concave Stieltjes, t_n(0+)=2/(3^n+1), t_n'(0)=3^n(5^n-1)/(6(3^n+1)^2); 2/(3^n+1)<2-E_n/Y_n<=2/(3^n+1)+B_a 3^-n so tan theta_n -> -1/(3 sqrt3) for every alpha; kappa_cross>0 at every level; c(alpha)=lim 3^n(2-E_n/Y_n) exists, certified c(0.45) in [3.72198897,3.72914181], c(0.8) in [2.57126273,2.57352573] (Kept Failure BD.1)"),
    ("mesh_noncentrality_monotone.py", [], "Theorem BE: Conjecture BD.2 decided in part: d_n(alpha)=2-E_n/Y_n strictly decreasing in alpha at every level (symmetrised exponential-tilt kernel), so c(alpha) strictly decreasing on (0,1); c(1-)=2 with certified slope 3S in [2.70632481,2.70903118]; Kept Failure BE.1: c(0+)=int T Y/int Y is finite, certified in [5.48553067,5.50249013], not infinite"),
    ("mesh_noncentrality_convex.py", [], "Theorem BF: proves Conjecture BE.2: c(1-eps)-2=eps N/(1/3+eps G) with shifted exponential moments smooth on [0,1]; interval certificate on 40 eps-intervals gives c''>=0.52046 on [0,1], so c(alpha) strictly convex; c'(0+) in [-4.341285,-4.227515]; tangent/chord bounds max(2+3S(1-a), c0+c'(0+)a) < c(a) < 2+(c0-2)(1-a)"),
    ("gasket_fractional_extension_limit.py", [], "Theorem BG: the L^alpha Dirichlet extension itself has a refinement limit: exact profile u_n(x)=dbar+int s^-a r_n h_m(p^(n-m)s;x)/int s^-a r_n; |u_n-dbar|<=(8/3)4^-a 5^(m(1-a))||x|| rho^n/R_a(n) on V_m, so flat (u->dbar) for alpha<=log3/log5; rho^-n(u_n(m_ab)-dbar)->-x_c U, U(0.45) in [0.143810693,0.143889868]; non-flat for alpha>log3/log5 (kappa_0.8 in [0.086072435,0.086089770]); Kept Failure BG.1: not uniform (corner boundary layer 0.3853654863)"),
    ("gasket_fractional_corner_layer.py", [], "Theorem BH: corner boundary layer of BG.1 has a limit profile: corner-cell decomposition; neighbour recursion h_N(z)=alpha1+beta h_(N-1)(pz) contracts (beta<=3/5) to x_a eta(z), eta=sum alpha1(p^j z)prod beta(p^i z), r_inf=1/(2+s-2eta); exact asymmetry u_n(nu_b)-u_n(nu_c)=(x_b-x_c)I_n/D_n; for alpha<gamma u_n(2^-n xi)->dbar+x_a Theta(xi) (universal), for alpha>=gamma ->d_a; certified theta_1(0.45) in [0.340381375,0.340407195], theta_2 in [0.241594653,0.241622693], theta_Y in [0.382431398,0.382456436]; rho^-k theta_k converges (from AT)"),
    ("gasket_fractional_uniform_limit.py", [], "Theorem BI: proves Conjecture BG.2: cell recursion omega_(i+1)<=(3/5)omega_i+z_i/(z_i+2) r_n||x|| for the resolvent, orbit lemma, R_n>=1/(3+tau), certified sup s^(1-gamma) r_inf<=0.6031; for alpha>gamma osc of u_n on level-N cells <=||x|| K_alpha rho^-N uniformly in n (K_0.8<=47.109), so u_n->continuous limit uniformly, =d on corners"),
    ("gasket_fractional_corner_sharpness.py", [], "Theorem BJ: proves Conjecture BI.7 (sharpness): exact 1/r_n-1/r_inf=prod beta(p^i s)(1+2eta(p^n s))<=3(3/5)^n, so 1/(tau^(1-gamma)/g+3)<=R_n<=tau^(gamma-1)g; certified g in [0.600945014,0.600963487] near 0; corner deficit rho^N delta_N>=c_-(alpha)>0 for every N (c_-(0.8)>=0.6208), so the BI modulus rho^-N is exact and the Hoelder exponent (alpha-gamma)/(1-gamma) optimal; liminf/limsup rho^N delta_N(0.8) in [1.135890346,1.135938783]"),
    ("gasket_fractional_uniform_rate.py", [], "Theorem BK: proves Conjecture BJ.7 (rate, c=1): exact one-level transfer u_(n+1)-dbar=<omega F_n>/<omega> under the level-n measure, omega=(1+y/5)^a/((1+2y/5)(5/3)(beta-e r_n(p)))>=1, <omega>=D'_(n+1)/D'_n; so sup_(V_n)|u_n-u-bar|<=2||x|| log(D'_inf/D'_n)<=C_a||x|| rho^-n; certified J_0.8<=1.188154, C_0.8<=6.2192, C_0.75<=12.5476, C_0.9<=1.9725"),
    ("gasket_fractional_rate_sharpness.py", [], "Theorem BL: sharpness half of Conjecture BK.6: 1-eta_N(z)=z/(z+2)+beta(1-eta_(N-1)(pz)) so eta_N strictly decreasing and the deficit integrand G_k increasing; increment delta^(k+1)-delta^(k)=<(omega-1)(G-<G>)>/<omega> >= rho^-k[(g-Delta)L(y0)-Delta J(y0)]/D'_(k+1) with omega>=omega_lo; certified at 0.8 (y0=1/10): positive for N=9..13 (direct) and all N>=14 (uniform), so delta_N^(n) strictly increases in n>=N and delta_N-delta_N^(n)>=c_N rho^-n (c_N>=0.1240 for N>=40): BK's rate rho^-n is exact. Not thrust."),
    ("gasket_fractional_rate_sharpness_alpha.py", [], "Theorem BM: BL's sharpness at every alpha in (gamma,1) for large N: L(y0)>0, J(y0)<=J_a<inf, Q<inf, 1-eta_N(p(y0))>=p(y0)/(p(y0)+2), and the uniform deficit ceiling ->0 like rho^-N, so delta_N-delta_N^(n)>=c_N rho^-n for N>=N_0(alpha); certified N_0=61,16,7,6,5 at alpha=0.7,0.75,0.85,0.9,0.95. Not thrust."),
    ("gasket_fractional_monotone_weight.py", [], "Theorem BN: BK's weight omega_k is increasing for every level and alpha>=1/17 (r_k Stieltjes of unit mass: r<=1/s, -r'<=r/s; d_k'-a d_k/(5+y)<=-P_a/(3(y+2)^2(y+5)^2)), so delta_N^(n) increases in n for every N>=1 (Chebyshev) and, via a pair form of the covariance, delta_N-delta_N^(n)>=c_N rho^-n for EVERY N>=1 at every alpha in (gamma,1); certified at 0.7-0.95, e.g. c_1(0.8)>=1.758e-3; not thrust"),
    ("gasket_fractional_weight_threshold.py", [], "Theorem BO: exact exponent threshold of BN's monotone weight: omega_k increasing at every level iff alpha >= a* in [0.016778880172852, 0.0167788801732] (local moments m_j=1,2,6,22,92,418,1994,9780 for k>=3, exact geometric remainders, exact-rational Sturm certificates); Kept Failure BO.1: below a* omega_k is not monotone at any level k>=3 (witness y*=13.778505644); not thrust"),
    ("gasket_fractional_deficit_small_alpha.py", [], "Theorem BP: delta_N^(n) strictly increases in n for EVERY alpha in (0,1), also below a* where omega_k is not monotone (BO.1): BN's majorant holds for every alpha>=0 so omega_k increases on [Y_a,inf), Y_a=(1-17a)/(2a); omega_k<(8/5)(1+y/5)^a<3*5^-a=rho^-1 on [0,Y_a]; <omega_k>=rho^-1 D_(k+1)/D_k>rho^-1 (AQ); single crossing gives Cov(omega_k,G_k)>0. Not thrust."),
    ("gasket_fractional_deficit_rate_small_alpha.py", [], "Theorem BQ: below gamma, 1-delta_N^(n) is of exact order rho^n at every level: the numerator calN^(n) increases (BK.1 transfer, omega>=1, eta_N>0) and D_n increases to D_inf, so calN^(m) rho^n/D_inf < 1-delta_N^(n) < P_N rho^n/D_m and rho^-n(1-delta_N^(n)) -> C_N>0; certified C_1(0.45) in [0.0718892,0.0719611]; Kept Failure BQ.1: rho^-n(1-delta) not monotone, not one-sided. Not thrust."),
    ("gasket_fractional_deficit_depth.py", [], "Theorem BR: depth dependence of the corner deficit: eta_N<=2/(2+z) (induction) and eta_(N+1)(w)-eta_N(pw)>=w/(w^2+5w+2), so delta_N^(n) strictly decreases in N at every level and every alpha; below gamma rho^N C_N strictly increases in N (Q_N=(5/3)^N R(5^N t) and eta_N both increase) to c_inf in [g_lo,g_hi] M_a/D_inf; certified c_inf(0.45) in [0.520387,0.520824]. Not thrust."),
    ("gasket_fractional_deficit_correction.py", [], "Theorem BS: below gamma, rho^N C_N -> c_inf at the exact rate rho^N: scaled BJ.1 gives 1/R = 1/Q_inf + Phi exactly (Q_inf = lim (3/5)^M r_inf(5^-M u), AT.1-free; Phi = (1+2 eta o psi) prod (5/3) beta(psi(5^-k u))), so c_inf - rho^N C_N = rho^N (A_N + B_N)/D_inf with A_N increasing to J, B_N <= (3/5)^(N-1) B_1; certified chi_inf(0.45) in [0.633344,0.633869] and C_N(0.45), N=1..6. Not thrust."),
    ("gasket_fractional_log_periodic_profile.py", [], "Theorem BT: the log-periodic profile F=u^(1-gamma)Q_inf of BS: g=s^(1-gamma)r_inf is nonincreasing along p-orbits (1<=g(s)/g(ps), every s), so Theta=lim g o q^k exists everywhere, F=Theta o psi (AT.1-free), F continuous, <F>=<Theta>, sup g=max F in [0.6009615785,0.6009631774]; Kept Failure BT.5: F not constant (max-min>=1.465e-5). Not thrust"),
    ("gasket_fractional_profile_fourier.py", [], "Theorem BU: BT's profile F=u^(1-gamma)Q_inf is AU's theta (psi_AU=lim 5^k q^k inverts BD.4's psi), so F_k=theta_k for every k; convexity of r_inf gives O(h^2) chord cell bounds: max F in [0.6009615785971,0.6009615785983], min F in [0.6009469237907,0.6009469237943], <F> to 2.4e-12, F_1, F_2 != 0 from F's own nodes (match AU); Kept Failure BU.4: F(t+1/2)=2<F>-F(t) false, (max+min)/2-<F> >= 2.43e-11"),
    ("gasket_fractional_profile_analytic.py", [], "Theorem BV: F=theta is holomorphic on |Im log_5 u|<pi/ln5 (Stieltjes sector bound + Vitali), so |theta_k|<=B_k<=e G*(1+4 pi k/ln5)e^(-2 pi^2 k/ln5); AU's routine with 384 cells certifies theta_1..theta_11 != 0 (+-8.2e-63); exactly one max and one min per period, G*=0.600961578597101089042056184..., min F=0.600946923794265612104226188... (+-1e-60); Kept Failure BV.5: t_max-t_min-1/2=-1.624019e-7 != 0"),
    ("gasket_fractional_profile_modes.py", [], "Theorem BW: AU's certified mode routine re-run at depth 5^-260 with 1024 cells and sandwich b=u=32/33 (300 digits, about 6 min) certifies theta_1..theta_20 != 0 (+-2.0e-113; BV: k<=11), |theta_k|<=B_k and |theta_k|e^(2 pi^2 k/ln5) in [0.60,0.93] for k<=20, <F> and c_R(gamma) to 1e-113, and with BV.3's tail G*=0.600961578597101089042056184... and min F to 1.3e-109, t_max, t_min to 5.6e-104. Not thrust"),
    ("gasket_corner_current.py", [], "harmonic gasket corner currents: neutral dipole scaling by 3/5"),
    ("gasket_fractional_currents.py", [], "Theorem H: fractional corner currents neutral; Kept Failure H.1 (no 3/5)"),
    ("gasket_geometric_currents.py", [], "Theorem I: geometric 1/d^2 currents grow by 12/5"),
    ("gasket_uniform_weights.py", [], "Theorem S: uniform multiplicative weights trichotomy; Kept Failure S.1"),
    ("gasket_multiscale_weights.py", [], "Theorem W: multi-scale multiplicative weights geometric-series trichotomy; Kept Failure W.1"),
    ("gasket_multiscale_fractional.py", [], "Theorem X: multi-scale spectral L^alpha neutrality; Kept Failures X.1-X.2"),
    ("gasket_fractional_limit.py", [], "Theorems J, K, L: fractional audit norm bounded, decays above log3/log5"),
    ("gasket_fractional_energy_floor.py", [], "Theorem AA: energy floor for L^alpha below alpha=1/2 (Loewner-Heinz path bound); Kept Failure AA.1"),
    ("gasket_fractional_share.py", [], "Theorem AB: combinatorial/fractional share sigma=rho/w_n; Kept Failure AB.1"),
    ("gasket_fractional_capacity_floor.py", [], "Theorem AC: energy floor for every alpha<log3/log5 (cell resistance + spectral mass); Kept Failure AC.1"),
    ("gasket_fractional_endpoint.py", [], "Theorem AD: exact pairing Q*R_alpha=7/3; at alpha=log3/log5 energy and audit norm -> 0 like 1/n; Kept Failure AD.1"),
    ("gasket_fractional_degree_share.py", [], "Theorem AE: degree-gap factorization of rho; sharp w_n<=8^a/2; Kept Failure AE.1"),
    ("gasket_fractional_lambda_floor.py", [], "Theorem AF: uniform floor lambda_n>=a 2^(-1-a) (Jensen a 2^(a-1)); ||F||->0 forces gamma_n->inf; Kept Failure AF.1"),
    ("gasket_fractional_lambda_ceiling.py", [], "Theorem AG: uniform ceiling lambda_n<=1-a(1-a)/2^(3-a); ||F||->0 iff gamma_n->inf; Kept Failure AG.1"),
    ("gasket_fractional_gamma_ceiling.py", [], "Theorem AH: kernel comparison + neighbor Harnack ceil gamma_n; ||F(n,a)||>0 floor for a<log3/log5 (6.44e-4 at 0.45); Kept Failure AH.1"),
    ("gasket_fractional_harnack_sum.py", [], "Theorem AI: two-neighbor heat-kernel sum w_pj+w_qj>=2w_aj; Harnack 4^a/2; ||F(n,0.45)||>=0.051988; Kept Failure AI.1"),
    ("gasket_fractional_wpq_harnack.py", [], "Theorem AJ: keep w_pq cross term; H_AJ=(4^a-a 4^(a-1))/2; ||F(n,0.45)||>=0.055599; Kept Failure AJ.1"),
    ("gasket_fractional_moment_jensen.py", [], "Theorem AK: spec(L) in [0,6] + two-moment Jensen; ||F(n,0.45)||>=0.058941; Kept Failure AK.1: live neighbor-Harnack chain < 0.70"),
    ("gasket_fractional_path_ceiling.py", [], "Theorem AL: Hermite-Hadamard path ceiling R_inf; ||F_frac||>=0.111208, ||F(n,0.45)||>=0.078019; Kept Failure AL.1: sharp path ceiling < 0.70"),
    ("gasket_fractional_cycle_ceiling.py", [], "Theorem AM: boundary-cycle resistance ceiling via Fourier; ||F_frac||>=0.209685, ||F(n,0.45)||>=0.147106; Kept Failure AM.1"),
    ("gasket_fractional_fence_ceiling.py", [], "Theorem AN: boundary-fence Bloch ceiling; ||F_frac||>=0.305377, ||F(n,0.45)||>=0.214239; Kept Failure AN.1"),
    ("gasket_fractional_moment_share.py", [], "Theorem AP: exact share sigma=(L^(1-a))_aa+omega; Radau corner-moment floor sigma>=1.310267; ||F(n,0.45)||>=0.274743 (AM) / 0.415238 (AO); Kept Failure AP.1"),
    ("gasket_fractional_decimation.py", [], "Theorem AQ: exact spectral decimation r_{n+1}(s)=Phi_s(r_n(s(s+5))); R_a(n) up, sigma_n down; ||F(n,0.45)|| decreases to L in [0.819255,0.822959] (Conjecture J.1); Kept Failure AQ.1"),
    ("gasket_fractional_decimation_rate.py", [], "Theorem AR: exact rates R(n+1)-R(n)~c_R(5^a/3)^n, sigma gap ~c_S(5^(a-1)/3)^n; slope c_R(gamma)=0.517166 at log3/log5; sigma_inf>=2*4^-a; Kept Failure AR.1"),
    ("gasket_fractional_limit_digits.py", [], "Theorem AS: complete-monotone interpolation sandwich certifies L(0.45) in [0.82113862713346,0.82113862713352]; L(0.10)>1.410 proves Conjecture J.2; Kept Failure AS.1"),
    ("gasket_fractional_rate_constants.py", [], "Theorem AT: certified rate constants c_R(0.45)=0.49780265076356035..., c_S(0.45)=0.92709263856946138..., c_R(gamma)=0.51716583463050524...; c_F for ||F||-L; Kept Failure AT.1 (log-periodic s^(1-gamma) r_inf)"),
    ("gasket_fractional_log_mean.py", [], "Theorem AU: proves Conjecture AT.2, c_R(gamma)=(2 sin(pi gamma)/pi) ln5 <Theta>; certified <Theta>=0.600954251168773607823705150099994..., c_R(gamma) to 34 digits; Mellin poles 1-gamma+2pi i k/ln5, theta_hat_1..4 nonzero; Kept Failure AU.1"),
    ("gasket_fractional_endpoint_shift.py", [], "Theorem AY: certified intercept C_gamma=0.703070462421857429929048..., R_gamma(m)=c_R m+C_gamma+O(5e-28) for m>=60; level shift nu_gamma=1.359468115..., sigma_inf(gamma)=1.151918611674..., (n+nu)||F(n,gamma)||->K_gamma=3.402360763832... exponentially; Kept Failure AY.1 (unshifted n||F|| converges only like 1/n)"),
    ("gasket_fractional_supercritical_shift.py", [], "Theorem AZ: alpha>log3/log5: R_a(n)=A rho^n+C_a+O(rho^n 5^-0.65n), intercept C_a=(2 sin pi a/pi) M(1-a) (Mellin transform of r_inf continued across the pole at gamma); (rho^n+nu)||F(n,a)||->K_a; certified C_0.8=-1.743257178392075723954..., nu_0.8=-0.727445384442..., K_0.8=0.687866758413992...; Kept Failure AZ.1 (shift sign opposite to the endpoint)"),
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
