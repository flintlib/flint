---
title: Arb 0.1 to 2.8
summary: Releases of Arb from 2012 to 2015: from the first ball arithmetic code to the mag, arf, arb and acb modules of Arb 2.0, power series, hypergeometric and other special functions.
order: 2
---

Arb, the library for arbitrary-precision ball arithmetic built on
FLINT 2 and merged into FLINT 3.0 in 2023, started with version 0.1 in
September 2012. The release notes below cover the versions up to 2.8.1
(December 2015), among them Arb 2.0 (May 2014), which introduced the
modules `mag`, `arf`, `arb` and `acb`; the later versions are in
[](history/merged/arb-2-9-to-2-23).

## Arb 2.8.1 (2015-12-31) {#arb-2-8-1}

* Fixed 32-bit test failure for the Laguerre function.
* Made the Laguerre function indeterminate at negative integer orders, to be consistent with the test code.

## Arb 2.8.0 (2015-12-29) {#arb-2-8-0}

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

## Arb 2.7.0 (2015-07-14) {#arb-2-7-0}

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

## Arb 2.6.0 (2015-04-19) {#arb-2-6-0}

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

## Arb 2.5.0 (2015-01-28) {#arb-2-5-0}

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

## Arb 2.4.0 (2014-11-15) {#arb-2-4-0}

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

## Arb 2.3.0 (2014-09-25) {#arb-2-3-0}

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

## Arb 2.2.0 (2014-08-01) {#arb-2-2-0}

* Added functions for computing polylogarithms and order expansions
  of polylogarithms, with support for real and complex s, z.
* Added a missing cast affecting C++ compatibility.
* Generalized powsum functions to allow a geometric factor.
* Improved powsum functions slightly when the exponent is an integer.
* Faster arb_log_ui_from_prev.
* Added mag_sqrt and mag_rsqrt functions.
* Fixed various minor bugs and added missing tests and documentation entries.

## Arb 2.1.0 (2014-06-20) {#arb-2-1-0}

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

## Arb 2.0.0 (2014-05-27) {#arb-2-0-0}

* New modules mag, arf, arb, arb_poly, arb_mat, acb, acb_poly,
  acb_mat for higher-performance ball arithmetic.

* Poly_roots2 and hilbert_matrix2 example programs.

* Vector dot product and norm functions (contributed by Abhinav Baid).

## Arb 1.1.0 (2014-05-03) {#arb-1-1-0}

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

## Arb 1.0.0 (2013-12-21) {#arb-1-0-0}

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

## Arb 0.7 (2013-08-07) {#arb-0-7}

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

## Arb 0.6 (2013-05-31) {#arb-0-6}

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

## Arb 0.5 (2013-03-28) {#arb-0-5}

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

## Arb 0.4 (2013-01-26) {#arb-0-4}

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

## Arb 0.3 (2012-11-07) {#arb-0-3}

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

## Arb 0.2 (2012-09-29) {#arb-0-2}

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

## Arb 0.1 (2012-09-14) {#arb-0-1}

* 2012-08-05 - Began simplified rewrite.
* 2012-04-05 - Experimental ball and polynomial code (first commit).
