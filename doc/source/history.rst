.. _history:

History and changes
===============================================================================

For more details, view the commit log
in the git repository https://github.com/fredrik-johansson/arb

Old releases of the code can be accessed from
https://github.com/fredrik-johansson/arb/releases

Old versions of the documentation
-------------------------------------------------------------------------------

* http://arblib.org/arb-2.23.0.pdf
* http://arblib.org/arb-2.22.0.pdf
* http://arblib.org/arb-2.21.0.pdf
* http://arblib.org/arb-2.20.0.pdf
* http://arblib.org/arb-2.19.0.pdf
* http://arblib.org/arb-2.18.1.pdf
* http://arblib.org/arb-2.18.0.pdf
* http://arblib.org/arb-2.17.0.pdf
* http://arblib.org/arb-2.16.0.pdf
* http://arblib.org/arb-2.15.0.pdf
* http://arblib.org/arb-2.14.0.pdf
* http://arblib.org/arb-2.13.0.pdf
* http://arblib.org/arb-2.12.0.pdf
* http://arblib.org/arb-2.11.1.pdf
* http://arblib.org/arb-2.11.0.pdf
* http://arblib.org/arb-2.10.0.pdf
* http://arblib.org/arb-2.9.0.pdf
* http://arblib.org/arb-2.8.1.pdf
* http://arblib.org/arb-2.8.0.pdf
* http://arblib.org/arb-2.7.0.pdf
* http://arblib.org/arb-2.6.0.pdf
* http://arblib.org/arb-2.5.0.pdf
* http://arblib.org/arb-2.4.0.pdf
* http://arblib.org/arb-2.3.0.pdf

2022-06-29 -- version 2.23.0
-------------------------------------------------------------------------------

* Performance

  * Multithreaded numerical integration.
  * Multithreaded binary splitting computation of mathematical constants.
  * Multithreaded computation of Bernoulli numbers.
  * Multithreaded computation of Euler numbers.
  * Multithreaded refinement of Riemann zeta zeros.
  * Multithreaded complex_plot example program.
  * Multithreaded elementary functions.
  * Multithreaded computation of Hilbert class polynomials.
  * Improved multithreaded partition function.
  * Use FLINT's FFT multiplication instead of GMP in appropriate ranges.
  * New, faster algorithm for elementary functions between roughly 10^3 and 10^6 digits.
  * Faster computation of log using Newton-like iteration instead of using MPFR.
  * Faster computation of atan using Newton-like iteration instead of the bit-burst algorithm.
  * Fix performance bug in atan() leading to quadratic running time with large arguments in high precision.
  * Optimized high-precision complex squaring.
  * Added internal function arb_flint_get_num_available_threads() to improve tuning for multithreaded algorithms
  * Fixed performance bug making erf() slower at high precision with multiple threads.

* Features

  * Implemented the Lerch transcendent (acb_dirichlet_lerch_phi()).
  * fpwrap wrapper for Lerch transcendent (contributed by Valentin Boettcher).
  * Added a rudimentary module for Gaussian integers (fmpzi.h).
  * Added zeta_zeros example program (contributed by D.H.J. Polymath).
  * Added functions for simultaneous high-precision computation of logarithms
    of primes and arctangents for primitive angles.
  * Added bernoulli, class_poly, functions_benchmark example programs
    for benchmarking use.
  * Multiplying a signed number by an infinity yields an infinity instead of
    [0 +/- inf] (contributed by Erik Postma).

* Miscellaneous

  * Deprecated doubles version of partition function.
  * Fix crash in erf on some systems including mips64el (reported by Julien Puydt).
  * Fixed MINGW64 build (contributed by Massoud Mazar).
  * Avoid deprecated FLINT function n_gcd_full.
  * Documentation fixes.


2022-01-25 -- version 2.22.1
-------------------------------------------------------------------------------

* Fixed bugs causing some hypergeometric functions hang or crash for
  some input on various non-x86 architectures.
* Fixed a minor bug in acb_hypgeom_m (NaN result sometimes only set the
  real part to NaN).

2022-01-15 -- version 2.22.0
-------------------------------------------------------------------------------

* Special functions

  * Use numerical integration in some cases to compute the hypergeometric
    functions 0F1, 1F1, U, 2F1, incomplete gamma and beta, modified Bessel, etc.
    with real parameters and argument, improving performance and
    accuracy when the parameters are large.
  * Much faster computation of Bernoulli numbers using hybrid numerical-modular
    algorithm (modular code adapted from bernmm by David Harvey).
  * Faster computation of Euler numbers using hybrid algorithm; added
    arb_fmpz_euler_number_ui.
  * Added inverse error function (arb_hypgeom_erfinv, arb_hypgeom_erfcinv).
  * New (faster, more accurate) implementations of real error functions
    (arb_hypgeom_erf, arb_hypgeom_erfc) and trigonometric integrals
    (arb_hypgeom_si, arb_hypgeom_ci).
  * Added acb_dirichlet_l_fmpq and acb_dirichlet_l_fmpq_afe:
    reduced-complexity evaluation of L-functions at rational points.
  * Added functions for computing primorials (arb_primorial_ui, arb_primorial_nth_ui).
  * New, highly optimized internal code for real hypergeometric series
    (arb_hypgeom_sum_fmpq_arb, etc.; currently only used in some functions).
  * Fix arb_fpwrap_double_hypgeom_2f1 which computed the wrong thing.

* Core arithmetic and functions

  * Faster implementation of arb_ui_pow_ui.
  * Added arb_fma_si, arb_fma_fmpz.
  * Added arf_equal_ui, arf_equal_d.
  * Added arf_get_str.
  * Use arb-based printing code instead of MPFR in arf_printd and mag_printd
    so that large exponents work.
  * Fixed bug in arb_get_str that caused loss of precision
    when printing more than about 10^6 digits.
  * Allow negative exponents in mag_pow_fmpz.
  * Added the double_interval module for fast machine-precision interval
    arithmetic (experimental, intended for internal use).

2021-10-20 -- version 2.21.1
-------------------------------------------------------------------------------

* Fixed 32-bit test failures for arb_hypgeom_gamma_fmpq.
* Added pow function to the fpwrap module.
* Added missing header file includes.
* Do not encode the library version in the SONAME on Android (contributed by Andreas Enge).

2021-09-25 -- version 2.21.0
-------------------------------------------------------------------------------

* Experimental new arb_fpwrap module: accurate floating-point wrappers
  of Arb mathematical functions (supersedes the external arbcmath.h).
* Fixed memory leak in arf_load_file (reported by Dave Platt).
* New and faster gamma function code.
* Most gamma function internals are now located in the arb_hypgeom and
  acb_hypgeom modules. The user-facing functions (arb_gamma, etc.) are still
  available under the old names for compatibility. The internal
  algorithms for rising factorials (binary splitting, etc.) have been moved
  without aliases.
* Added arb_fma, arb_fma_arf, arb_fma_ui (like addmul, but take a separate input and output).
* Slightly faster internal Bernoulli number generation for small n.
* Better enclosures for acb_barnes_g at negative reals.
* Added Graeffe transforms (arb_poly_graeffe_transform, acb_poly_graeffe_transform)
  (contributed by Matthias Gessinger).
* Fixed conflict with musl libc (reported by Gonzalo Tornaría).
* Added acb_add_error_arb (contributed by Albin Ahlbäck).

2021-07-25 -- version 2.20.0
-------------------------------------------------------------------------------

* Flint 2.8 support.
* Change arb_get_str with ARB_STR_NO_RADIUS: [+/- 1.20e-15] now prints as 0e-14.
* Uniformly distributed random number functions arf_urandom, arb_urandom
  (contributed by Albin Ahlbäck).
* Use quasilinear algorithm in arb_gamma_fmpq for all small fractions.
* Added derivative of Weierstrass elliptic function (acb_elliptic_p_prime)
  (contributed by Daniel Schultz).
* Added dot products with integer coefficients: arb_dot_fmpz, arb_dot_siui,
  arb_dot_uiui, arb_dot_si, arb_dot_ui, acb_dot_fmpz, acb_dot_siui,
  acb_dot_uiui, acb_dot_si, acb_dot_ui.
* Faster arb_fmpz_poly_evaluate_arb and arb_fmpz_poly_evaluate_acb.
* Explicitly guarantee that roots are isolated in arb_fmpz_poly_complex_roots
  (could previously theoretically fail when using the deflation hack).
* Use GNUInstallDirs in CMakeLists.txt to support standard GNU installation
  directories (contributed by Michael Orlitzky).
* Fixed bug for aliased multiplication of window matrices (contributed by
  David Berghaus).
* Documentation fixes (contributed by Joel Dahne, Hanno Rein).

2020-12-06 -- version 2.19.0
-------------------------------------------------------------------------------

* Significant improvements to the implementation of Platt's algorithm for
  computing Riemann zeta function zeros at large height (contributed by
  p15-git-acc).
* Better criterion for selecting asymptotic expansion of incomplete gamma
  function (contributed by p15-git-acc).
* Multithreaded acb_dft for power-of-two lengths (contributed by p15-git-acc).
* Added acb_csc_pi, arb_csc_pi (contributed by p15-git-acc).
* Fixed segfault in acb_mat_eig_simple_rump when called with L non-NULL
  and R NULL (contributed by p15-git-acc).
* Fixed bug in acb_real_abs (contributed by Joel Dahne).
* Changed several functions to more consistently return infinities
  instead of NaNs where reasonable (contributed by p15-git-acc).
* Added Fransen-Robinson as an integral example (contributed by p15-git-acc).
* Cleaned up makefile (contributed by p15-git-acc).
* Fixed several typos and some omitted functions in the documentation
  (contributed by Joel-Dahne, p15-git-acc).


2020-06-25 -- version 2.18.1
-------------------------------------------------------------------------------

* Support MinGW64.
* Added version numbers (__ARB_VERSION, __ARB_RELEASE, ARB_VERSION) to arb.h.

2020-06-09 -- version 2.18.0
-------------------------------------------------------------------------------

* General

  * Flint 2.6 support.
  * Several build system improvements (contributed by Isuru Fernando).
  * Changed arf_get_mpfr to return an MPFR underflow/overflow result
    (rounding to 0 or infinity with the right sign and MPFR overflow flags)
    instead of throwing flint_abort() if the exponent is out of bounds for MPFR.
  * Documentation and type corrections (contributed by Joel Dahne).

* Arithmetic

  * The number of iterations per precision level in arb_fmpz_poly_complex_roots
    has been tweaked to avoid extreme slowdown for some polynomials with
    closely clustered roots.
  * Added arb_contains_interior, acb_contains_interior.

* Special functions

  * Fixed unsafe shifts causing Dirichlet characters for certain moduli
    exceeding 32 bits to crash.
  * Added acb_agm for computing the arithmetic-geometric mean of two complex
    numbers.
  * acb_elliptic_rj now uses a slow fallback algorithm in cases where Carlson's
    algorithm is not known to be valid. This fixes instances where
    acb_elliptic_pi, acb_elliptic_pi_inc and acb_elliptic_rj previously ended
    up on the wrong branch. Users should be cautioned that the new version can
    give worse enclosures and sometimes fails to converge in some cases where
    the old algorithm did (the pi flag for acb_elliptic_pi_inc is useful as a
    workaround).
  * Optimized some special cases in acb_hurwitz_zeta.

2019-10-16 -- version 2.17.0
-------------------------------------------------------------------------------

* General

  * Added exact serialization methods (arb_dump_str, arb_load_str, arb_dump_file,
    arb_load_file, arf_dump_str, arf_load_str, arf_dump_file, arf_load_file,
    mag_dump_str, mag_load_str, mag_dump_file, mag_load_file)
    (contributed by Julian Rüth).
  * Removed many obsolete fmpr methods and de-inlined several helper functions
    to slightly improve compile time and library size.
  * Fixed a namespace clash for an internal function (contributed by Julian Rüth).
  * Added the helper function arb_sgn_nonzero.
  * Added the helper function acb_rel_one_accuracy_bits.

* Riemann zeta function

  * Added a function for efficiently computing individual zeros of the Riemann
    zeta function using Turing's method (acb_dirichlet_zeta_zero)
    (contributed by D.H.J. Polymath).
  * Added a function for counting zeros of the Riemann zeta function up to
    given height using Turing's method (acb_dirichlet_zeta_nzeros)
    (contributed by D.H.J. Polymath).
  * Added the Backlund S function (acb_dirichlet_backlund_s).
  * Added a function for computing Gram points (acb_dirichlet_gram_point).
  * Added acb_dirichlet_zeta_deriv_bound for quickly bounding the derivative
    of the Riemann zeta function.
  * Fast multi-evaluation of the Riemann zeta function using Platt's algorithm
    (acb_dirichlet_platt_multieval) (contributed by D.H.J. Polymath).

* Other special functions

  * Improved the algorithm in acb_hypgeom_u to estimate precision loss
    more accurately.
  * Implemented Coulomb wave functions (acb_hypgeom_coulomb,
    acb_hypgeom_coulomb_series and other functions).
  * Faster algorithm for Catalan's constant.
  * Added acb_modular_theta_series.
  * Added arb_poly_sinc_pi_series (contributed by D.H.J. Polymath).
  * Improved tuning in acb_hypgeom_pfq_series_sum for higher derivatives
    at high precision (reported by Mark Watkins).


2018-12-07 -- version 2.16.0
-------------------------------------------------------------------------------

* Linear algebra and arithmetic

  * Added acb_mat_approx_eig_qr for approximate computation of eigenvalues
    and eigenvectors of complex matrices.
  * Added acb_mat_eig_enclosure_rump implementing Rump's algorithm for
    certification of eigenvalue-eigenvector pairs as well as clusters.
  * Added acb_mat_eig_simple_rump for certified diagonalization of matrices
    with simple eigenvalues.
  * Added acb_mat_eig_simple_vdhoeven_mourrain, acb_mat_eig_simple for fast
    certified diagonalization of matrices with simple eigenvalues.
  * Added acb_mat_eig_multiple_rump, acb_mat_eig_multiple for certified
    computation of eigenvalues with possible overlap.
  * Added acb_mat_eig_global_enclosure for fast global inclusion of eigenvalues
    without isolation.
  * Added arb_mat_companion, acb_mat_companion for constructing companion
    matrices.
  * Added several arb_mat and acb_mat helper functions: indeterminate,
    is_exact, is_zero, is_finite, is_triu, is_tril, is_diag, diag_prod.
  * Added arb_mat_approx_inv, acb_mat_approx_inv.
  * Optimized arb_mat_mul_block by using arb_dot when the blocks are small.
  * Added acb_get_mid.
  * Updated hilbert_matrix example program.


2018-10-25 -- version 2.15.1
-------------------------------------------------------------------------------

* Fixed precision issue leading to spurious NaN results in incomplete elliptic integrals

2018-09-18 -- version 2.15.0
-------------------------------------------------------------------------------

* Arithmetic

  * Added arb_dot and acb_dot for efficient evaluation of dot products.
  * Added arb_approx_dot and acb_approx_dot for efficient evaluation of dot products without error bounds.
  * Converted loops to arb_dot and acb_dot in the arb_poly and acb_poly methods mullow_classical, inv_series, div_series, exp_series_basecase, sin_cos_series_basecase, sinh_cosh_series_basecase, evaluate_rectangular, evaluate2_rectangular, revert_series_lagrange_fast. Also changed the algorithm cutoffs for mullow, exp_series, sin_cos_series, sinh_cosh_series.
  * Converted loops to arb_dot and acb_dot in the arb_mat and acb_mat methods mul_classical, mul_threaded, solve_tril, solve_triu, charpoly. Also changed the algorithm cutoffs for mul, solve_tril, solve_triu.
  * Converted loops to arb_approx_dot and acb_approx_dot in the arb_mat and acb_mat methods approx_solve_tril, approx_solve_triu. Also changed the algorithm cutoffs.
  * Added arb_mat_approx_mul and acb_mat_approx_mul for matrix multiplication without error bounds.

* Miscellaneous

  * Added arb_hypgeom_airy_zero for computing zeros of Airy functions.
  * Added arb_hypgeom_dilog wrapper.
  * Optimized arb_const_pi and arb_const_log2 by using a static table at low precision, giving a small speedup and avoiding common recomputation when starting threads.
  * Optimized mag_set_ui_2exp_si.
  * Remove obsolete and unused function _arb_vec_dot.
  * Converted some inline functions to ordinary functions to reduce library size.
  * Fixed acb_dirichlet_stieltjes to use the integration algorithm also when a != 1.
  * Fixed test failure for acb_dirichlet_stieltjes on ARM64 (reported by Gianfranco Costamagna). Special thanks to Julien Puydt for assistance with debugging.
  * Fixed crash in acb_dft_bluestein with zero length (reported by Gianfranco Costamagna).

2018-07-22 -- version 2.14.0
-------------------------------------------------------------------------------

* Linear algebra

  * Faster and more accurate real matrix multiplication using block decomposition, scaling, and multiplying via FLINT integer matrices in combination with safe use of doubles for radius matrix multiplications.
  * Faster and more accurate complex matrix multiplication by reordering and taking advantage of real matrix multiplication.
  * The new multiplication algorithm methods (arb_mat_mul_block, acb_mat_mul_reorder) are used automatically by the main multiplication methods.
  * Faster and more accurate LU factorization by using a block recursive algorithm that takes advantage of matrix multiplication. Added separate algorithm methods (arb/acb)_mat_lu_(recursive/classical) with an automatic algorithm choice in the default lu methods.
  * Added methods (arb/acb)_mat_solve_(tril/triu) (and variants) for solving upper or lower triangular systems using a block recursive algorithm taking advantage of matrix multiplication.
  * Improved linear solving and inverse for large well-conditioned matrices by using a preconditioning algorithm. Added separate solving algorithm methods (arb/acb)_mat_solve_(lu/precond) with an automatic algorithm choice in the default solve methods (contributed by anonymous user arbguest).
  * Improved determinants using a preconditioning algorithm. Added separate determinant algorithm methods (arb/acb)_mat_det_(lu/precond) with an automatic algorithm choice in the default det methods.
  * Added automatic detection of triangular matrices in arb_mat_det and acb_mat_det.
  * Added arb_mat_solve_preapprox which allows certifying a precomputed approximate solution (contributed by anonymous user arbguest).
  * Added methods for constructing various useful test matrices: arb_mat_ones, arb_mat_hilbert, arb_mat_pascal, arb_mat_stirling, arb_mat_dct, acb_mat_ones, acb_mat_dft.
  * Added support for window matrices (arb/acb_mat_window_init/clear).
  * Changed random test matrix generation (arb/acb_mat_randtest) to produce sparse matrices with higher probability.
  * Added acb_mat_conjugate and acb_mat_conjugate_transpose.

* Arithmetic and elementary functions

  * Improved arb_sin_cos, arb_sin and arb_cos to produce more accurate enclosures for wide input intervals. The working precision is also reduced automatically based on the accuracy of the input to improve efficiency.
  * Improved arb_sinh_cosh, arb_sinh and arb_cosh to produce more accurate enclosures for wide input intervals. The working precision is also reduced automatically based on the accuracy of the input to improve efficiency.
  * Improved arb_exp_invexp and arb_expm1 to produce more accurate enclosures for wide input intervals. The working precision is also reduced automatically based on the accuracy of the input to improve efficiency.
  * Improved acb_rsqrt to produce more accurate enclosures for wide intervals.
  * Made mag_add_ui_lower public.
  * Added mag_sinh, mag_cosh, mag_sinh_lower, mag_cosh_lower.
  * Fixed minor precision loss near -1 in arb_log_hypot and acb_log.
  * Return imaginary numbers with exact zero real part when possible in acb_acos and acb_acosh (contributed by Ralf Stephan).
  * Improved special cases in arb_set_interval_arf (reported by Marc Mezzarobba).

* Special functions

  * Added a function for computing isolated generalized Stieltjes constants (acb_dirichlet_stieltjes).
  * Added scaled versions of Bessel functions (acb_hypgeom_bessel_i_scaled, acb_hypgeom_bessel_k_scaled).
  * The interface for the internal methods computing Bessel functions (i_asymp, k_asymp, etc.) has been changed to accommodate computing scaled versions.
  * Added Riemann xi function (acb_dirichlet_xi) (contributed by D.H.J Polymath).
  * Fixed infinite error bounds in the Riemann zeta function when evaluating at a ball containing zero centered in the left plane (contributed by D.H.J Polymath).
  * Fixed precision loss in Airy functions with huge input and high precision.
  * Legendre functions of the first kind (legendre_p): handle inexact integer a+b-c in 2F1 better (contributed by Joel Dahne).

* Example programs and documentation

  * Added more color functions to complex_plot.c.
  * Added more example integrals suggested by Nicolas Brisebarre and Bruno Salvy to integrals.c
  * Changed Sphinx style and redesigned the documentation front page.
  * Miscellaneous documentation cleanups.
  * Added documentation page about contributing.

* Other

  * Fixed a crash on some systems when calling acb_dft methods with a length of zero.
  * Fixed issue with setting rpath in configure (contributed by Vincent Delecroix).


2018-02-23 -- version 2.13.0
-------------------------------------------------------------------------------

* Major bugs

  * Fixed rounding direction in arb_get_abs_lbound_arf() which in some cases
    could result in an invalid lower bound being returned, and added forgotten
    test code for this and related functions (reported by deinst). Although
    this bug could lead to incorrect results, it probably had limited impact in
    practice (explaining why it was not caught indirectly by other test code)
    since a single rounding in the wrong direction in this operation generally
    will be dwarfed by multiple roundings in the correct direction in
    surrounding operations.

* Important notes about bounds

  * Many functions have been modified to compute tighter enclosures
    when the input balls are wide. In most cases the bounds should be
    improved, but there may be some regressions. Bug reports about any
    significant regressions are welcome.
  * Division by zero in arb_div() has been changed to return [NaN +/- inf]
    instead of [+/- inf]. This change might be reverted in the future if it
    proves to be too inconvenient. In either case, users should only assume
    that division by zero produces something non-finite, and user code that
    depends on division by zero to produce [0 +/- inf] should be modified to
    handle zero-containing denominators as a separate case.

* Improvements to arithmetic and elementary functions

  * Faster implementation of acb_get_mag_lower().
  * Optimized arb_get_mag_lower(), arb_get_mag_lower_nonnegative().
  * Added arb_set_interval_mag() and arb_set_interval_neg_pos_mag() for
    constructing an arb_t from a pair of mag_t endpoints.
  * Added mag_const_pi_lower(), mag_atan(), mag_atan_lower().
  * Added mag_div_lower(), mag_inv(), mag_inv_lower().
  * Added mag_sqrt_lower() and mag_rsqrt_lower().
  * Added mag_log(), mag_log_lower(), mag_neg_log(), mag_neg_log_lower().
  * Added mag_exp_lower(), mag_expinv_lower() and tweaked mag_exp().
  * Added mag_pow_fmpz_lower(), mag_get_fmpz(), mag_get_fmpz_lower().
  * Improved arb_exp() for wide input.
  * Improved arb_log() for wide input.
  * Improved arb_sqrt() for wide input.
  * Improved arb_rsqrt() for wide input.
  * Improved arb_div() for wide input.
  * Improved arb_atan() for wide input and slightly optimized arb_atan2()
    for input spanning multiple signs.
  * Improved acb_rsqrt() for wide input and improved stability of this
    function generally in the left half plane.
  * Added arb_log_hypot() and improved acb_log() for wide input.
  * Slightly optimized trigonometric functions (acb_sin(), acb_sin_pi(),
    acb_cos(), acb_cos_pi(), acb_sin_cos(), acb_sin_cos_pi()) for pure real or
    imaginary input.

* Special functions

  * Slightly improved bounds for gamma function (arb_gamma(), acb_gamma(),
    arb_rgamma(), acb_rgamma()) for wide input.
  * Improved bounds for Airy functions for wide input.
  * Simplifications to code for computing Gauss period minimal polynomials
    (contributed by Jean-Pierre Flori).
  * Optimized arb_hypgeom_legendre_p_ui() further by avoiding divisions in the
    basecase recurrence and computing the prefactor more quickly in the
    asymptotic series (contributed by Marc Mezzarobba).
  * Small further optimization of arb_hypgeom_legendre_p_ui_root()
    (contributed by Marc Mezzarobba).
  * Improved derivative bounds for Legendre polynomials (contributed by
    Marc Mezzarobba).

* Numerical integration

  * Increased default quadrature deg_limit at low precision to improve
    performance for integration of functions without singularities near the
    path.
  * Added several more integrals to examples/integrals.c
  * Added utility functions acb_real_abs(), acb_real_sgn(),
    acb_real_heaviside(), acb_real_floor(), acb_real_ceil(), acb_real_min(),
    acb_real_max(), acb_real_sqrtpos(), useful for numerical integration.
  * Added utility functions acb_sqrt_analytic(), acb_rsqrt_analytic(),
    acb_log_analytic(), acb_pow_analytic() with branch cut detection, useful
    for numerical integration.

* Build system and compatibility issues

  * Removed -Wl flag from Makefile.subdirs to fix "-r and -pie may not be used
    together" compilation error on some newer Linux distributions (reported
    by many users).
  * Fixed broken test code for l_vec_hurwitz which resulted in spurious
    failures on 32-bit systems (originally reported by Thierry Monteil on
    Sage trac).
  * Avoid using deprecated MPFR function mpfr_root() with MPFR
    versions >= 4.0.0.
  * Remark: the recently released MPFR 4.0.0 has a bug in mpfr_div() leading
    to test failures in Arb (though not affecting correctness of Arb itself).
    Users should make sure to install the patched version MPFR 4.0.1.
  * Added missing C++ include guards in arb_fmpz_poly.h and dlog.h (reported
    by Marc Mezzarobba).
  * Fixed Travis builds on Mac OS again (contributed by Isuru Fernando).
  * Added missing declaration for arb_bell_ui() (reported by numsys).

2017-11-29 - version 2.12.0
-------------------------------------------------------------------------------

* Numerical integration

  * Added a new function (acb_calc_integrate) for rigorous numerical
    integration using adaptive subdivision and Gauss-Legendre quadrature. This
    largely obsoletes the old integration code using Taylor series.
  * Added new integrals.c example program (old example program moved to
    integrals_taylor.c).

* Discrete Fourier transforms

  * Added acb_dft module with various FFT algorithm implementations, including
    top level O(n log n) acb_dft and acb_dft_inverse functions
    (contributed by Pascal Molin).

* Legendre polynomials

  * Added arb_hypgeom_legendre_p_ui for fast and accurate evaluation of
    Legendre polynomials. This is also used automatically by the Legendre
    functions, where it is substantially faster and gives better error
    bounds than the generic algorithm.
  * Added arb_hypgeom_legendre_p_ui_root for fast computation of Legendre
    polynomial roots and Gauss-Legendre quadrature nodes (used internally
    by the new integration code).
  * Added arb_hypgeom_central_bin_ui for fast computation of central
    binomial coefficients (used internally for Legendre polynomials).

* Dirichlet L-functions and zeta functions

  * Fixed a bug in the Riemann zeta function involving a too small error
    bound in the implementation of the Riemann-Siegel formula for inexact
    input. This bug could result in a too small enclosure when evaluating the
    Riemann zeta function at an argument of large imaginary height without
    also computing derivatives, if the input interval was very wide.
  * Add acb_dirichlet_zeta_jet; also made computation of the first derivative
    of Riemann zeta function use the Riemann-Siegel formula where appropriate.
  * Added acb_dirichlet_l_vec_hurwitz for fast simultaneous evaluation of
    Dirichlet L-functions for multiple characters using Hurwitz zeta function
    and FFT (contributed by Pascal Molin).
  * Simplified interface for using hurwitz_precomp functions.
  * Added lcentral.c example program (contributed by Pascal Molin).
  * Improved error bounds when evaluating Dirichlet L-functions using
    Euler product.

* Elementary functions

  * Faster custom implementation of sin, cos at 4600 bits and above
    instead of using MPFR (30-40% asymptotic improvement, up to a factor
    two speedup).
  * Faster code for exp between 4600 and 19000 bits.
  * Improved error bounds for acb_atan by using derivative.
  * Improved error bounds for arb_sinh_cosh, arb_sinh and arb_cosh when
    the input has a small midpoint and large radius.
  * Added reciprocal trigonometric and hyperbolic functions (arb_sec, arb_csc,
    arb_sech, arb_csch, acb_sec, acb_csc, acb_sech, acb_csch).
  * Changed the interface of _acb_vec_unit_roots to take an extra length
    parameter (compatibility-breaking change).
  * Improved arb_pow and acb_pow with an inexact base and a negative integer
    or negative half-integer exponent; the inverse is now computed before
    performing binary exponentiation in this case to avoid spurious blow-up.

* Elliptic functions

  * Improved Jacobi theta functions to reduce the argument modulo the lattice
    parameter, greatly improving speed and numerical stability for large input.
  * Optimized arb_agm by using a final series expansion and using special code
    for wide intervals.
  * Optimized acb_agm1 by using a final series expansion and using special code
    for positive real input.
  * Optimized derivative of AGM for high precision by using a central
    difference instead of a forward difference.
  * Optimized acb_elliptic_rf and acb_elliptic_rj for high precision by using
    a variable length series expansion.

* Other

  * Fixed incorrect handling of subnormals in arf_set_d.
  * Added mag_bin_uiui for bounding binomial coefficients.
  * Added mag_set_d_lower, mag_sqrt_lower, mag_set_d_2exp_fmpz_lower.
  * Implemented multithreaded complex matrix multiplication.
  * Optimized arb_rel_accuracy_bits by adding fast path.
  * Fixed a spurious floating-point exception (division by zero) in the
    t-gauss_period_minpoly test program triggered by new code optimizations
    in recent versions of GCC that are unsafe together with FLINT inline
    assembly functions (a workaround was added to the test code, and a proper
    fix for the assembly code has been added to FLINT).

2017-07-10 - version 2.11.1
-------------------------------------------------------------------------------

* Avoid use of a function that was unavailable in the latest public FLINT release

2017-07-09 - version 2.11.0
-------------------------------------------------------------------------------

* Special functions

  * Added the Lambert W function (arb_lambertw, acb_lambertw, arb_poly_lambertw_series, acb_poly_lambertw_series). All complex branches and evaluation of derivatives are supported.
  * Added the acb_expm1 method, complementing arb_expm1.
  * Added arb_sinc_pi, acb_sinc_pi.
  * Optimized handling of more special cases in the Hurwitz zeta function.

* Polynomials

  * Added the arb_fmpz_poly module to provide Arb methods for FLINT integer polynomials.
  * Added methods for evaluating an fmpz_poly at arb_t and acb_t arguments.
  * Added arb_fmpz_poly_complex_roots for computing the real and complex roots of an integer polynomial, turning the functionality previously available in the poly_roots.c example program into a proper library function.
  * Added a method (arb_fmpz_poly_gauss_period_minpoly) for constructing minimal polynomials of Gaussian periods.
  * Added arb_poly_product_roots_complex for constructing a real polynomial from complex conjugate roots.

* Miscellaneous

  * Fixed test code in the dirichlet module for 32-bit systems (contributed by Pascal Molin).
  * Use flint_abort() instead of abort() (contributed by Tommy Hofmann).
  * Fixed the static library install path (contributed by François Bissey).
  * Made arb_nonnegative_part() a publicly documented method.
  * Arb now requires FLINT version 2.5 or later.

2017-02-27 - version 2.10.0
-------------------------------------------------------------------------------

* General

  * Changed a large number of methods from inline functions to normal
    functions, substantially reducing the size of the built library.
  * Fixed a few minor memory leaks (missing clear() calls).

* Basic arithmetic

  * Added arb_is_int_2exp_si and acb_is_int_2exp_si.
  * Added arf_sosq for computing x^2+y^2 of floating-point numbers.
  * Improved error bounds for complex square roots in the left half plane.
  * Improved error bounds for complex reciprocal (acb_inv) and division.
  * Added the internal helper mag_get_d_log2_approx as a public method.

* Elliptic functions and integrals

  * New module acb_elliptic.h for elliptic functions and integrals.
  * Added complete elliptic integral of the third kind.
  * Added Legendre incomplete elliptic integrals (first, second, third kinds).
  * Added Carlson symmetric incomplete elliptic integrals (RF, RC, RG, RJ, RD).
  * Added Weierstrass elliptic zeta and sigma functions.
  * Added inverse Weierstrass elliptic p-function.
  * Added utility functions for computing the Weierstrass invariants and lattice roots.
  * Improved computation of derivatives of Jacobi theta functions by
    using modular transformations, and added a main evaluation function
    (acb_modular_theta_jet).
  * Improved detection of pure real or pure imaginary parts in various cases
    of evaluating theta and modular functions.

* Other special functions

  * New, far more efficient implementation of the dilogarithm function (acb_polylog with s = 2).
  * Fixed an issue in the Hurwitz zeta function leading to unreasonable
    slowdown for certain complex input.
  * Added add acb_poly_exp_pi_i_series.
  * Added arb_poly_log1p_series, acb_poly_log1p_series.

2016-12-02 - version 2.9.0
-------------------------------------------------------------------------------

* License

  * Changed license from GPL to LGPL.

* Build system and compatibility

  * Fixed FLINT includes to use flint/foo.h instead of foo.h, simplifying compilation on many systems.
  * Added another alias for the dynamic library to fix make check on certain systems (contributed by Andreas Enge).
  * Travis CI support (contributed by Isuru Fernando).
  * Added support for ARB_TEST_MULTIPLIER environment variable to control the number of test iterations.
  * Support building with CMake (contributed by Isuru Fernando).
  * Support building with MSVC on Windows (contributed by Isuru Fernando).
  * Fixed unsafe use of FLINT_ABS for slong -> ulong conversion in arf.h,
    which caused failures on MIPS and ARM systems.

* Basic arithmetic and methods

  * Fixed mag_addmul(x,x,x) with x having a mantissa of all ones. This could
    produce a non-normalized mag_t value, potentially leading to
    incorrect results in arb and acb level arithmetic. This bug was caught by
    new test code, and fortunately would have been hard to trigger accidentally.
  * Added fasth paths for error bound calculations in arb_sqrt and arb_div, speeding up these operations significantly at low precision
  * Added support for round-to-nearest in all arf methods.
  * Added fprint methods (contributed by Alex Griffing).
  * Added acb_printn and acb_fprintn methods to match arb_printn.
  * Added arb_equal_si and acb_equal_si.
  * Added arb_can_round_mpfr.
  * Added arb_get_ubound_arf, arb_get_lbound_arf (contributed by Tommy Hofmann).
  * Added sign function (arb_sgn).
  * Added complex sign functions (acb_sgn, acb_csgn).
  * Rewrote arb_contains_fmpq to make the test exact.
  * Optimized mag_get_fmpq.
  * Optimized arf_get_fmpz and added more robust test code.
  * Rewrote arb_get_unique_fmpz and arb_get_interval_fmpz_2exp, reducing overhead, making them more robust with huge exponents, and documenting their behavior more carefully.
  * Optimized arb_union.
  * Optimized arf_is_int, arf_is_int_2exp_si and changed these from inline to normal functions.
  * Added mag_const_pi, mag_sub, mag_expinv.
  * Optimized binary-to-decimal conversion for huge exponents by using exponential function instead of binary powering.
  * Added arb_intersection (contributed by Alex Griffing).
  * Added arb_min, arb_max (contributed by Alex Griffing).
  * Fixed a bug in arb_log and in test code on 64-bit Windows due to unsafe use of MPFR which only uses 32-bit exponents on Win64.
  * Improved some test functions to reduce the chance of reporting spurious failures.
  * Added squaring functions (arb_sqr, acb_sqr) (contributed by Ricky Farr).
  * Added arf_frexp.
  * Added arf_cmp_si, arf_cmp_ui, arf_cmp_d.
  * Added methods to count allocated bytes (arb_allocated_bytes, _arb_vec_allocated_bytes, etc.).
  * Added methods to predict memory usage for large vectors (_arb/_acb_vec_estimate_allocated_bytes).
  * Changed clear() methods from inline to normal functions, giving 8% faster compilation and 25% smaller libarb.so.
  * Added acb_unit_root and _acb_vec_unit_roots (contributed by Pascal Molin).

* Polynomials

  * Added sinh and cosh functions of power series (arb/acb_poly_sinh/cosh_series and sinh_cosh_series).
  * Use basecase series inversion algorithm to improve speed and error bounds in arb/acb_poly_inv_series.
  * Added functions for fast polynomial Taylor shift (arb_poly_taylor_shift, acb_poly_taylor_shift and variants).
  * Fast handling of special cases in polynomial composition.
  * Added acb_poly scalar mul and div convenience methods (contributed by Alex Griffing).
  * Added set_trunc, set_trunc_round convenience methods.
  * Added add_series, sub_series methods for truncating addition.
  * Added polynomial is_zero, is_one, is_x, valuation convenience methods.
  * Added hack to arb_poly_mullow and acb_poly_mullow to avoid overhead when doing an in-place multiplication with length at most 2.
  * Added binomial and Borel transform methods for acb_poly.

* Matrices

  * Added Cholesky decomposition plus solving and inverse
    for positive definite matrices (arb_mat_cho, arb_mat_spd_solve, arb_mat_spd_inv
    and related methods) (contributed by Alex Griffing).
  * Added LDL decomposition and inverse and solving based on LDL decomposition
    for real matrices (arb_mat_ldl, arb_mat_solve_ldl_precomp, arb_mat_inv_ldl_precomp)
    (contributed by Alex Griffing).
  * Improved the entrywise error bounds in matrix exponential computation
    to preserve sparsity and give exact entries where possible in many cases
    (contributed by Alex Griffing).
  * Added public functions for computing the truncated matrix exponential
    Taylor series (arb_mat_exp_taylor_sum, acb_mat_exp_taylor_sum).
  * Added functions related to sparsity structure (arb_mat_entrywise_is_zero,
    arb_mat_count_is_zero, etc.) (contributed by Alex Griffing).
  * Entrywise multiplication (arb_mat_mul_entrywise, acb_mat_mul_entrywise)
    (contributed by Alex Griffing).
  * Added is_empty and is_square convenience methods (contributed by Alex Griffing).
  * Added the bool_mat helper module for matrices over the boolean semiring (contributed by Alex Griffing).
  * Added Frobenius norm computation (contributed by Alex Griffing).

* Miscellaneous special functions

  * Added evaluation of Bernoulli polynomials (arb_bernoulli_poly_ui, acb_bernoulli_poly_ui).
  * Added convenience function for evaluation of huge Bernoulli numbers (arb_bernoulli_fmpz).
  * Added Euler numbers (arb_euler_number_ui, arb_euler_number_fmpz).
  * Added fast approximate partition function (arb_partitions_fmpz/ui).
  * Optimized partition function for n < 1000 by using recurrence for the low 64 bits.
  * Improved the worst-case error bound in arb_atan.
  * Added arb_log_base_ui.
  * Added complex sinc function (acb_sinc).
  * Special handling of z = 1 when computing polylogarithms.
  * Fixed agm(-1,-1) to output 0 instead of indeterminate.
  * Made working precision in arb_gamma and acb_gamma more sensitive to the input accuracy.

* Hypergeometric functions

  * Compute erf and erfc without cancellation problems for large or complex z.
  * Avoid re-computing the square root of pi in several places.
  * Added generalized hypergeometric function (acb_hypgeom_pfq).
  * Implement binary splitting and rectangular splitting for evaluation of hypergeometric series with a power series parameter, greatly speeding up Y_n, K_n and other functions at high precision, as well as speeding up high-order parameter derivatives.
  * Use binary splitting more aggressively in acb_hypgeom_pfq_sum to reduce error bound inflation.
  * Asymptotic expansions of hypergeometric functions: more accurate parameter selection, and better handling of terminating cases.
  * Tweaked algorithm selection and working precision in acb_hypgeom_m.
  * Avoid dividing by the denominator of the next term in acb_hypgeom_sum, which would lead to a division by zero when evaluating hypergeometric polynomials.
  * Fixed a bug in hypergeometric series evaluation resulting in near-integers not being skipped in some cases, leading to unnecessary loss of precision.
  * Added series expansions of Airy functions (acb_hypgeom_airy_series, acb_hypgeom_airy_jet).
  * Fixed a case where Airy functions accidentally chose the worst algorithm instead of the best one.
  * Added functions for computing erf, erfc, erfi of power series in the acb_hypgeom module.
  * Added series expansion of the logarithmic integral (acb_hypgeom_li_series).
  * Added Fresnel integrals (acb_hypgeom_fresnel, acb_hypgeom_fresnel_series).
  * Added the lower incomplete gamma function (acb_hypgeom_gamma_lower) (contributed by Alex Griffing).
  * Added series expansion of the lower incomplete gamma function (acb_hypgeom_gamma_lower_series) (contributed by Alex Griffing).
  * Added support for computing the regularized incomplete gamma functions.
  * Use slightly sharper error bound for analytic continuation of 2F1.
  * Added support for computing finite limits of 2F1 with inexact parameters differing by integers.
  * Added the incomplete beta function (acb_hypgeom_beta_lower, acb_hypgeom_beta_lower_series)
  * Improved acb_hypgeom_u to use a division-avoiding algorithm for small polynomial cases.
  * Added arb_hypgeom module, wrapping the complex hypergeometric functions for more convenient use with the arb_t type.

* Dirichlet L-functions and Riemann zeta function

  * New module dirichlet for working algebraically with Dirichlet groups and characters (contributed by Pascal Molin).
  * New module acb_dirichlet for numerical evaluation of Dirichlet characters and L-functions (contributed by Pascal Molin).
  * Efficient representation and manipulation of Dirichlet characters using the Conrey representation (contributed by Pascal Molin).
  * New module dlog for word-size discrete logarithm evaluation, used to support algorithms on Dirichlet characters (contributed by Pascal Molin).
  * Methods for properties, evaluation, iteration, pairing, lift, lowering etc. of Dirichlet characters (contributed by Pascal Molin).
  * Added acb_dirichlet_roots methods for fast evaluation of many roots of unity (contributed by Pascal Molin).
  * Added acb_dirichlet_hurwitz_precomp methods for fast multi-evaluation of the Hurwitz zeta function for many parameter values.
  * Added methods for computing Gauss, Jacobi and theta sums over Dirichlet characters (contributed by Pascal Molin).
  * Added methods (acb_dirichlet_l, acb_dirichlet_l_jet, acb_dirichlet_l_series) for evaluation of Dirichlet L-functions and their derivatives.
  * Implemented multiple algorithms for evaluation of Dirichlet L-functions depending on the argument (Hurwitz zeta function decomposition, Euler product, functional equation).
  * Added methods (acb_dirichlet_hardy_z, acb_dirichlet_hardy_z_series, etc.) for computing the Hardy Z-function corresponding to a Dirichlet L-function.
  * Added fast bound for Hurwitz zeta function (mag_hurwitz_zeta_uiui).
  * Improved parameter selection in Hurwitz zeta function to target relative
    instead of absolute error for large positive s.
  * Improved parameter selection in Hurwitz zeta function to avoid computing
    unnecessary Bernoulli numbers for large imaginary s.
  * Added Dirichlet eta function (acb_dirichlet_eta).
  * Implemented the Riemann-Siegel formula for faster evaluation of the Riemann zeta function at large height.
  * Added smooth-index algorithm for the main sum when evaluating the Riemann zeta function, avoiding the high memory usage of the full sieving algorithm when the number of terms gets huge.
  * Improved tuning for using the Euler product when computing the Riemann zeta function.

* Example programs

  * Added logistic map example program.
  * Added lvalue example program.
  * Improved poly_roots in several ways: identify roots that are exactly real,
    automatically perform squarefree factorization, use power hack, and
    allow specifying a product of polynomials as input on the command line.

* Housekeeping

  * New section in the documentation giving an introduction to ball arithmetic and using the library.
  * Tidied, documented and added test code for the fmpz_extras module.
  * Added proper documentation and test code for many helper methods.
  * Removed the obsolete fmprb module entirely.
  * Documented more algorithms and formulas.
  * Clarified integer overflow issues and use of ARF_PREC_EXACT in the documentation.
  * Added .gitignore file.
  * Miscellaneous improvements to the documentation.

2015-12-31 - version 2.8.1
-------------------------------------------------------------------------------

* Fixed 32-bit test failure for the Laguerre function.
* Made the Laguerre function indeterminate at negative integer orders, to be consistent with the test code.

2015-12-29 - version 2.8.0
-------------------------------------------------------------------------------

* Compatibility and build system

  * Windows64 support (contributed by Bill Hart).
  * Fixed a bug that broke basic arithmetic on targets where FLINT uses fallback code instead of assembly code, such as PPC64 (contributed by Jeroen Demeyer).
  * Fixed configure to use EXTRA_SHARED_FLAGS/LDFLAGS, and other build system fixes (contributed by Tommy Hofmann, Bill Hart).
  * Added soname versioning (contributed by Julien Puydt).
  * Fixed test code on MinGW (contributed by Hrvoje Abraham).
  * Miscellaneous fixes to simplify interfacing Arb from Julia.

* Arithmetic and elementary functions

  * Fixed arf_get_d to handle underflow/overflow correctly and to support round-to-nearest.
  * Added more complex inverse hyperbolic functions (acb_asin, acb_acos, acb_asinh, acb_acosh, acb_atanh).
  * Added arb_contains_int and acb_contains_int for testing whether an interval contains any integer.
  * Added acb_quadratic_roots_fmpz.
  * Improved arb_sinh to use a more accurate formula for x < 0.
  * Added sinc function (arb_sinc) (contributed by Alex Griffing).
  * Fixed bug in arb_exp affecting convergence for huge input.
  * Faster implementation of arb_div_2expm1_ui.
  * Added mag_root, mag_geom_series.
  * Improved and added test code for arb_add_error functions.
  * Changed arb_pow and acb_pow to make pow(0,positive) = 0 instead of nan.
  * Improved acb_sqrt to return finite output for finite input straddling the branch cut.
  * Improved arb_set_interval_arf so that [inf,inf] = inf instead of an infinite interval.
  * Added computation of Bell numbers (arb_bell_fmpz).
  * Added arb_power_sum_vec for computing power sums using Bernoulli numbers.
  * Added computation of the Fujiwara root bound for acb_poly.
  * Added code to identify all the real roots of a real polynomial (acb_poly_validate_real_roots).
  * Added several convenient assignment functions, including arb_set_d, acb_set_d, acb_set_d_d, acb_set_fmpz_fmpz (contributed by Ricky Farr).
  * Added many accessor functions (_arb/acb_vec_entry_ptr, arb_get_mid/rad_arb, acb_real/imag_ptr, arb_mid/rad_ptr, acb_get_real/imag).
  * Added missing functions acb_add_si, acb_sub_si.
  * Renamed arb_root to arb_root_ui (keeping alias) and added acb_root_ui.

* Special functions

  * Implemented the Gauss hypergeometric function 2F1 and its regularized version.
  * Fixed two bugs in acb_hypgeom_pfq_series_direct discovered while implementing 2F1. In rare cases, these could lead to incorrect values for functions depending on parameter derivatives of hypergeometric series.

    * The first bug involved incorrect handling of negative integer parameters. The bug only affected 2F1 and higher functions; it did not affect correctness of any previously implemented functions that relied on acb_hypgeom_pfq_series_direct (such as Bessel Y and K functions of integer order).
    * The second bug involved a too small bound being computed for the sum of a geometric series. The geometric series bound is nearly tight for 2F1, and the incorrect version caused immediate test failures for that function. Theoretically, this bug affected correctness of some previously-implemented functions that relied on acb_hypgeom_pfq_series_direct (such as Bessel Y and K functions of integer order), but since the geometric bound is not as tight in those cases, those functions were still reliable in practice (no failing test case has been found).

  * Implemented Airy functions and their derivatives (acb_hypgeom_airy).
  * Implemented the confluent hypergeometric function 0F1 (acb_hypgeom_0f1).
  * Implemented associated Legendre functions P and Q.
  * Implemented Chebyshev, Jacobi, Gegenbauer, Laguerre, Hermite functions.
  * Implemented spherical harmonics.
  * Added function for computing Bessel J and Y functions simultaneously.
  * Added rising factorials for non-integer n (arb_rising, acb_rising).
  * Made rising factorials use gamma function for large integer n.
  * Faster algorithm for theta constants and Dedekind eta function at very high precision.
  * Fixed erf to give finite values instead of +/-inf for big imaginary input.
  * Improved acb_zeta (and arb_zeta) to automatically use fast code for integer zeta values.
  * Added double factorial (arb_doublefac_ui).
  * Added code for generating Hilbert class polynomials (acb_modular_hilbert_class_poly).

* Matrices

  * Added faster matrix squaring (arb/acb_mat_sqr) (contributed by Alex Griffing).
  * Added matrix trace (arb/acb_mat_trace) (contributed by Alex Griffing).
  * Added arb/acb_mat_set_round_fmpz_mat, acb_mat_set(_round)_arb_mat (contributed by Tommy Hofmann).
  * Added arb/acb_mat_transpose (contributed by Tommy Hofmann).
  * Added comparison methods arb/acb_mat_eq/ne (contributed by Tommy Hofmann).

* Other

  * Added complex_plot example program.
  * Added Airy functions to real_roots example program.
  * Other minor patches were contributed by Alexander Kobel, Marc Mezzarobba, Julien Puydt.
  * Removed obsolete file config.h.

2015-07-14 - version 2.7.0
-------------------------------------------------------------------------------

* Hypergeometric functions

  * Implemented Bessel I and Y functions (acb_hypgeom_bessel_i, acb_hypgeom_bessel_y).
  * Fixed bug in Bessel K function giving the wrong branch for negative real arguments.
  * Added code for evaluating complex hypergeometric series binary splitting.
  * Added code for evaluating complex hypergeometric series using fast multipoint evaluation.

* Gamma related functions

  * Implemented the Barnes G-function and its continuous logarithm (acb_barnes_g, acb_log_barnes_g).
  * Implemented the generalized polygamma function (acb_polygamma).
  * Implemented the reflection formula for the logarithmic gamma function (acb_lgamma, acb_poly_lgamma_series).
  * Implemented the digamma function of power series (arb_poly_digamma_series, acb_poly_digamma_series).
  * Improved acb_poly_zeta_series to produce exact zero imaginary parts in most cases when the result should be real-valued.
  * Made the real logarithmic gamma function (arb_lgamma, arb_poly_lgamma_series) abort more quickly for negative input.

* Elementary functions

  * Added arb_exp_expinv and acb_exp_expinv functions for simultaneously computing exp(x), exp(-x).
  * Improved acb_tan, acb_tan_pi, acb_cot and acb_cot_pi for input with large imaginary parts.
  * Added complex hyperbolic functions (acb_sinh, acb_cosh, acb_sinh_cosh, acb_tanh, acb_coth).
  * Added acb_log_sin_pi for computing the logarithmic sine function without branch cuts away from the real line.
  * Added arb_poly_cot_pi_series, acb_poly_cot_pi_series.
  * Added arf_root and improved speed of arb_root.
  * Tuned algorithm selection in arb_pow_fmpq.

* Other

  * Added documentation for arb and acb vector functions.

2015-04-19 - version 2.6.0
-------------------------------------------------------------------------------

* Special functions

  * Added the Bessel K function.
  * Added the confluent hypergeometric functions M and U.
  * Added exponential, trigonometric and logarithmic integrals ei, si, shi, ci, chi, li.
  * Added the complete elliptic integral of the second kind E.
  * Added support for computing hypergeometric functions with power series as parameters.
  * Fixed special cases in Bessel J function returning useless output.
  * Fixed precision of zeta function accidentally being capped at 7000 digits (bug in 2.5).
  * Special-cased real input in the gamma functions for complex types.
  * Fixed exp of huge numbers outputting unnecessarily useless intervals.
  * Fixed broken code in erf that sometimes gave useless output.
  * Made selection of number of terms in hypergeometric series more robust.

* Polynomials and power series.

  * Added sin_pi, cos_pi and sin_cos_pi for real and complex power series.
  * Speeded up series reciprocal and division for length = 2.
  * Added add_si methods for polynomials.
  * Made inv_series and div_series with zero input produce indeterminates instead of aborting.
  * Added arb_poly_majorant, acb_poly_majorant.

* Basic functions

  * Added comparison methods arb_eq, arb_ne, arb_lt, arb_le, arb_gt, arb_ge, acb_eq, acb_ne.
  * Added acb_rel_accuracy_bits and improved the real version.
  * Fixed precision of constants like pi behaving more nondeterministically than necessary.
  * Fixed arf_get_mag_lower(nan) to output 0 instead of inf.

* Other

  * Removed call to fmpq_dedekind_sum which only exists in the git version of flint.
  * Fixed a test code bug that could cause crashes on some systems.
  * Added fix for static build on OS X (thanks Marcello Seri).
  * Miscellaneous corrections to the documentation.

2015-01-28 - version 2.5.0
-------------------------------------------------------------------------------

* String conversion

  * Added arb_set_str.
  * Added arb_get_str and arb_printn for pretty-printed rigorous decimal output.
  * Added helper functions for binary to decimal conversion.

* Core arithmetic

  * Improved speed of division when using GMP instead of MPIR.
  * Improved complex division with a small denominator.
  * Removed a little bit of overhead for complex squaring.

* Special functions

  * Faster code for atan at very high precision, used instead of mpfr_atan.
  * Optimized elementary functions slightly for small input.
  * Added modified error functions erfc and erfi.
  * Added the generalized exponential integral.
  * Added the upper incomplete gamma function.
  * Implemented the complete elliptic integral of the first kind.
  * Implemented the arithmetic-geometric mean of complex numbers.
  * Optimized arb_digamma for small integers.
  * Made mag_log_ui, mag_binpow_uiui and mag_polylog_tail proper functions.
  * Added pow, agm, erf, elliptic_k, elliptic_p as functions of complex power series.
  * Added incomplete gamma function of complex power series.
  * Improved code for bounding complex rising factorials (the old code could
    potentially have given wrong results in degenerate cases).
  * Added arb_sqrt1pm1, arb_atanh, arb_asinh, arb_atanh.
  * Added arb_log1p, acb_log1p, acb_atan.
  * Added arb_hurwitz_zeta.
  * Improved parameter selection in the Hurwitz zeta function to try to
    avoid stalling when given enormous input.
  * Optimized sqrt and rsqrt of power series when given a binomial as input.
  * Made arb_bernoulli_ui(2^64-2) not crash.
  * Fixed rgamma of negative integers returning indeterminate.

* Polynomials and matrices

  * Added characteristic polynomial computation for real and complex matrices.
  * Added polynomial set_round methods.
  * Added is_real methods for more types.
  * Added more get_unique_fmpz methods.
  * Added code for generating Swinnerton-Dyer polynomials.
  * Improved error bounding in det() and exp() of complex matrices to
    recognize when the result is real-valued.
  * Changed polynomial divrem to return success/fail instead of aborting on divide by zero.

* Miscellaneous

  * Added logo to documentation.
  * Made inlined functions build as part of the library.
  * Silenced a clang warning.
  * Made _acb_vec_sort_pretty a library function.

2014-11-15 - version 2.4.0
-------------------------------------------------------------------------------

* Arithmetic and core functions

  * Made evaluation of sin, cos and exp at medium precision faster using the sqrt trick.
  * Optimized arb_sinh and arb_sinh_cosh.
  * Optimized complex division with a small denominator.
  * Optimized cubing of complex numbers.
  * Added floor and ceil functions for the arf and arb types.
  * Added acb_poly powering functions.
  * Added acb_exp_pi_i.
  * Added functions for evaluation of Chebyshev polynomials.
  * Fixed arb_div to output nan for input containing nan.

* Added a module acb_hypgeom for hypergeometric functions

  * Evaluation of the generalized hypergeometric function in convergent cases.
  * Evaluation of confluent hypergeometric functions using asymptotic expansions.
  * The Bessel function of the first kind for complex input.
  * The error function for complex input.

* Added a module acb_modular for modular forms and elliptic functions

  * Support for working with modular transformations.
  * Mapping a point to the fundamental domain.
  * Evaluation of Jacobi theta functions and their series expansions.
  * The Dedekind eta function.
  * The j-invariant and the modular lambda and delta function.
  * Eisenstein series.
  * The Weierstrass elliptic function and its series expansion.

* Miscellaneous

  * Fixed mag_print printing a too large exponent.
  * Fixed printd methods to use a fallback instead of aborting when printing numbers too large for MPFR.
  * Added version number string (arb_version).
  * Various additions to the documentation.

2014-09-25 - version 2.3.0
-------------------------------------------------------------------------------

* Removed most of the legacy (Arb 1.x) modules.
* Updated build scripts, hopefully fixing various issues.
* New implementations of arb_sin, arb_cos, arb_sin_cos, arb_atan, arb_log, arb_exp, arb_expm1, much faster up to a few thousand bits.
* Ported the bit-burst code for high-precision exponentials to the arb type.
* Speeded up arb_log_ui_from_prev.
* Added mag_exp, mag_expm1, mag_exp_tail, mag_pow_fmpz.
* Improved various mag functions.
* Added arb_get/set_interval_mpfr, arb_get_interval_arf, and improved arb_set_interval_arf.
* Improved arf_get_fmpz.
* Prettier printing of complex numbers with negative imaginary part.
* Changed some frequently-used functions from inline to non-inline to reduce code size.

2014-08-01 - version 2.2.0
-------------------------------------------------------------------------------

* Added functions for computing polylogarithms and order expansions
  of polylogarithms, with support for real and complex s, z.
* Added a missing cast affecting C++ compatibility.
* Generalized powsum functions to allow a geometric factor.
* Improved powsum functions slightly when the exponent is an integer.
* Faster arb_log_ui_from_prev.
* Added mag_sqrt and mag_rsqrt functions.
* Fixed various minor bugs and added missing tests and documentation entries.

2014-06-20 - version 2.1.0
-------------------------------------------------------------------------------

* Ported most of the remaining functions to the new arb/acb types,
  including:

  * Elementary functions (log, atan, etc.).
  * Hypergeometric series summation.
  * The gamma function.
  * The Riemann zeta function and related functions.
  * Bernoulli numbers.
  * The partition function.
  * The calculus modules (rigorous real root isolation, rigorous numerical integration of complex-valued functions).
  * Example programs.

* Added several missing utility functions to the arf and mag modules.

2014-05-27 - version 2.0.0
-------------------------------------------------------------------------------

* New modules mag, arf, arb, arb_poly, arb_mat, acb, acb_poly,
  acb_mat for higher-performance ball arithmetic.

* Poly_roots2 and hilbert_matrix2 example programs.

* Vector dot product and norm functions (contributed by Abhinav Baid).

2014-05-03 - version 1.1.0
-------------------------------------------------------------------------------

* Faster and more accurate error bounds for polynomial multiplication
  (error bounds are now always as good as with classical multiplication,
  and multiplying high-degree polynomials with approximately equal
  coefficients now has proper quasilinear complexity).

* Faster and much less memory-hungry exponentials at very high precision.

* Improved the partition function to support n bigger than a single word,
  and enabled the possibility to use two threads for the computation.

* Fixed a bug in floating-point arithmetic that caused a too small bound
  for the rounding error to be reported when the result of an inexact
  operation was rounded up to a power of two (this bug did
  not affect the correctness of ball arithmetic, because operations on
  ball midpoints always round down).

* Minor optimizations to floating-point arithmetic.

* Improved argument reduction of the digamma function and short series
  expansions of the rising factorial.

* Removed the holonomic module for now, as it did not really do anything
  very useful.

2013-12-21 - version 1.0.0
-------------------------------------------------------------------------------

* New example programs directory

  * poly_roots example program.
  * real_roots example program.
  * pi_digits example program.
  * hilbert_matrix example program.
  * keiper_li example program.

* New fmprb_calc module for calculus with real functions

  * Bisection-based root isolation.
  * Asymptotically fast Newton root refinement.

* New fmpcb_calc module for calculus with complex functions

  * Numerical integration using Taylor series.

* Scalar functions

  * Simplified fmprb_const_euler using published error bound.
  * Added fmprb_inv.
  * Added fmprb_trim, fmpcb_trim.
  * Added fmpcb_rsqrt (complex reciprocal square root).
  * Fixed bug in fmprb_sqrtpos with nonfinite input.
  * Slightly improved fmprb powering code.
  * Added various functions for bounding fmprs by powers of two.
  * Added fmpr_is_int.

* Polynomials and power series

  * Implemented scaling to speed up blockwise multiplication.
  * Slightly faster basecase power series exponentials.
  * Improved sin/cos/tan/exp for short power series.
  * Added complex sqrt_series, rsqrt_series.
  * Implemented the Riemann-Siegel Z and theta functions for real power series.
  * Added fmprb_poly_pow_series, fmprb_poly_pow_ui and related methods.
  * Added fmprb/fmpcb_poly_contains_fmpz_poly.
  * Faster composition by monomials.
  * Implemented Borel transform and binomial transform for real power series.

* Matrices

  * Implemented matrix exponentials.
  * Multithreaded fmprb_mat_mul.
  * Added matrix infinity norm functions.
  * Added some more matrix-scalar functions.
  * Added matrix contains and overlaps methods.

* Zeta function evaluation

  * Multithreaded power sum evaluation.
  * Faster parameter selection when computing many derivatives.
  * Implemented binary splitting to speed up computing many derivatives.

* Miscellaneous

  * Corrections for C++ compatibility (contributed by Jonathan Bober).
  * Several minor bugfixes and test code enhancements.

2013-08-07 - version 0.7
-------------------------------------------------------------------------------

* Floating-point and ball functions

  * Documented, added test code, and fixed bugs for various operations involving a ball containing an infinity or NaN.
  * Added reciprocal square root functions (fmpr_rsqrt, fmprb_rsqrt) based on mpfr_rec_sqrt.
  * Faster high-precision division by not computing an explicit remainder.
  * Slightly faster computation of pi by using new reciprocal square root and division code.
  * Added an fmpr function for approximate division to speed up certain radius operations.
  * Added fmpr_set_d for conversion from double.
  * Allow use of doubles to optionally compute the partition function faster but without an error bound.
  * Bypass mpfr overflow when computing the exponential function to extremely high precision (approximately 1 billion digits).
  * Made fmprb_exp faster for large numbers at extremely high precision by skipping the log(2) removal.
  * Made fmpcb_lgamma faster at high precision by speeding up the argument reduction branch computation.
  * Added fmprb_asin, fmprb_acos.
  * Added various other utility functions to the fmprb module.
  * Added a function for computing the Glaisher constant.
  * Optimized evaluation of the Riemann zeta function at high precision.

* Polynomials and power series

  * Made squaring of polynomials faster than generic multiplication.
  * Implemented power series reversion (various algorithms) for the fmprb_poly type.
  * Added many fmprb_poly utility functions (shifting, truncating, setting/getting coefficients, etc.).
  * Improved power series division when either operand is short
  * Improved power series logarithm when the input is short.
  * Improved power series exponential to use the basecase algorithm for short input regardless of the output size.
  * Added power series square root and reciprocal square root.
  * Added atan, tan, sin, cos, sin_cos, asin, acos fmprb_poly power series functions.
  * Added Newton iteration macros to simplify various functions.
  * Added gamma functions of real and complex power series ([fmprb/fmpcb]_poly_[gamma/rgamma/lgamma]_series).
  * Added wrappers for computing the Hurwitz zeta function of a power series ([fmprb/fmpcb]_poly_zeta_series).
  * Implemented sieving and other optimizations to improve performance for evaluating the zeta function of a short power series.
  * Improved power series composition when the inner series is linear.
  * Added many fmpcb_poly versions of nearly all fmprb_poly functions.
  * Improved speed and stability of series composition/reversion by balancing the power table exponents.

* Other

  * Added support for freeing all cached data by calling flint_cleanup().
  * Introduced fmprb_ptr, fmprb_srcptr, fmpcb_ptr, fmpcb_srcptr typedefs for cleaner function signatures.
  * Various bug fixes and general cleanup.

2013-05-31 - version 0.6
-------------------------------------------------------------------------------

* Made fast polynomial multiplication over the reals numerically stable by using a blockwise algorithm.
* Disabled default use of the Gauss formula for multiplication of complex polynomials, to improve numerical stability.
* Added division and remainder for complex polynomials.
* Added fast multipoint evaluation and interpolation for complex polynomials.
* Added missing fmprb_poly_sub and fmpcb_poly_sub functions.
* Faster exponentials (fmprb_exp and dependent functions) at low precision, using precomputation.
* Rewrote fmpr_add and fmpr_sub using mpn level code, improving efficiency at low precision.
* Ported the partition function implementation from flint (using ball arithmetic
  in all steps of the calculation to guarantee correctness).
* Ported algorithm for computing the cosine minimal polynomial from flint (using
  ball arithmetic to guarantee correctness).
* Support using GMP instead of MPIR.
* Only use thread-local storage when enabled in flint.
* Slightly faster error bounding for the zeta function.
* Added some other helper functions.

2013-03-28 - version 0.5
-------------------------------------------------------------------------------

* Arithmetic and elementary functions

  * Added fmpr_get_fmpz, fmpr_get_si.
  * Fixed accuracy problem with fmprb_div_2expm1.
  * Special-cased squaring of complex numbers.
  * Added various fmpcb convenience functions (addmul_ui, etc).
  * Optimized fmpr_cmp_2exp_si and fmpr_cmpabs_2exp_si, and added test code for comparison functions.
  * Added fmprb_atan2, also fixing a bug in fmpcb_arg.
  * Added fmprb_sin_pi, cos_pi, sin_cos_pi, etc.
  * Added fmprb_sin_pi_fmpq (etc.) using algebraic methods for fast evaluation of roots of unity.
  * Faster fmprb_poly_evaluate and evaluate_fmpcb using rectangular splitting.
  * Added fmprb_poly_evaluate2, evaluate2_fmpcb for simultaneously evaluating the derivative.
  * Added fmprb_poly root polishing code using near-optimal Newton steps (experimental).
  * Added fmpr_root, fmprb_root (currently based on MPFR).
  * Added fmpr_min, fmpr_max.
  * Added fmprb_set_interval_fmpr, fmprb_union.
  * Added fmpr_bits, fmprb_bits, fmpcb_bits for obtaining the mantissa width.
  * Added fmprb_hypot.
  * Added complex square roots.
  * Improved fmprb_log to slightly improve speed, and properly support huge arguments.
  * Fixed exp, cosh, sinh to work with huge arguments.
  * Added fmprb_expm1.
  * Fixed sin, cos, atan to work with huge arguments.
  * Improved fmprb_pow and fmpcb_pow, including automatic detection of small integer and half-integer exponents.
  * Added many more elementary functions: fmprb_tan/cot/tanh/coth, fmpcb_tan/cot, and pi versions.
  * Added fmprb const_e, const_log2, const_log10, const_catalan.
  * Fixed ball containment/overlap checking to work operate efficiently and correctly with huge exponents.
  * Strengthened test code for many core operations.

* Special functions

  * Reorganized zeta function related code.
  * Faster evaluation of the Riemann zeta function via sieving.
  * Documented and improved efficiency of the zeta constant binary splitting code.
  * Calculate error bound in Borwein's algorithm with fmprs instead of using doubles.
  * Optimized divisions in zeta evaluation via the Euler product.
  * Use functional equation for Riemann zeta function of a negative argument.
  * Compute single Bernoulli numbers using ball arithmetic instead of relying on the floating-point code in flint.
  * Initial code for evaluating the gamma function using its Taylor series.
  * Much faster rising factorials at high precision, using difference polynomials.
  * Much faster gamma function at high precision.
  * Added complex gamma function, log gamma function, and other versions.
  * Added fmprb_agm (real arithmetic-geometric mean).
  * Added fmprb_gamma_fmpq, supporting rapid computation of gamma(p/q) for q = 1,2,3,4,6.
  * Added real and complex digamma function.
  * Fixed unnecessary recomputation of Bernoulli numbers.
  * Optimized computation of Euler's constant, and added proper error bounds.
  * Avoid reliance on doubles in the hypergeometric series tail bound.
  * Cleaned up factorials and binomials, computing factorials via gamma.

* Other

  * Added an fmpz_extras module to collect various internal fmpz helper functions.
  * Fixed detection of flint header files.
  * Fixed various other small bugs.

2013-01-26 - version 0.4
-------------------------------------------------------------------------------

* Much faster fmpr_mul, fmprb_mul and set_round, resulting in general speed improvements.
* Code for computing the complex Hurwitz zeta function with derivatives.
* Fixed and documented error bounds for hypergeometric series.
* Better algorithm for series evaluation of the gamma function at a rational point.
* Much faster generation of Bernoulli numbers.
* Complex log, exp, pow, trigonometric functions (currently based on MPFR).
* Complex nth roots via Newton iteration.
* Added code for arithmetic on fmpcb_polys.
* Code for computing Khinchin's constant.
* Code for rising factorials of polynomials or power series
* Faster sin_cos.
* Better div_2expm1.
* Many other new helper functions.
* Improved thread safety.
* More test code for core operations.

2012-11-07 - version 0.3
-------------------------------------------------------------------------------

* Converted documentation to Sphinx.
* New module fmpcb for ball interval arithmetic over the complex numbers

  * Conversions, utility functions and arithmetic operations.

* New module fmpcb_mat for matrices over the complex numbers

  * Conversions, utility functions and arithmetic operations.
  * Multiplication, LU decomposition, solving, inverse and determinant.

* New module fmpcb_poly for polynomials over the complex numbers

  * Root isolation for complex polynomials.

* New module fmpz_holonomic for functions/sequences
  defined by linear differential/difference equations
  with polynomial coefficients

  * Functions for creating various special sequences and functions.
  * Some closure properties for sequences.
  * Taylor series expansion for differential equations.
  * Computing the nth entry of a sequence using binary splitting.
  * Computing the nth entry mod p using fast multipoint evaluation.

* Generic binary splitting code with automatic error bounding is now
  used for evaluating hypergeometric series.
* Matrix powering.
* Various other helper functions.

2012-09-29 - version 0.2
-------------------------------------------------------------------------------

* Code for computing the gamma function (Karatsuba, Stirling's series).
* Rising factorials.
* Fast exp_series using Newton iteration.
* Improved multiplication of small polynomials by using classical multiplication.
* Implemented error propagation for square roots.
* Polynomial division (Newton-based).
* Polynomial evaluation (Horner) and composition (divide-and-conquer).
* Product trees, fast multipoint evaluation and interpolation (various algorithms).
* Power series composition (Horner, Brent-Kung).
* Added the fmprb_mat module for matrices of balls of real numbers.
* Matrix multiplication.
* Interval-aware LU decomposition, solving, inverse and determinant.
* Many helper functions and small bugfixes.

2012-09-14 - version 0.1
-------------------------------------------------------------------------------

* 2012-08-05 - Began simplified rewrite.
* 2012-04-05 - Experimental ball and polynomial code (first commit).

