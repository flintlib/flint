.. _history:

History and changes
===============================================================================

FLINT version history
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

2023-11-10 -- FLINT 3.0.1
-------------------------------------------------------------------------------

* Build issues

  * Fix LIBS2 order for static linking (Tomás Oliveira e Silva).
  * Fix substitution of version number for older autotools (Albin Ahlbäck).
  * Fix use of AC_SEARCH_LIBS to find cblas_dgemm (Gonzalo Tornaría).
  * Add FlexiBLAS as a cblas option (Mahrud Sayrafi).
  * Don't use deprecated PythonInterp in CMake build (Mahrud Sayrafi).
  * Fix setting version numbers and strings in CMake build (Mahrud Sayrafi).
  * Only link with NTL for the tests on CMake (Mahrud Sayrafi).

* Bugs

  * Fix bug in nmod32 on 32-bit systems.
  * Fix missing modulus assignment in nmod_poly_mat_window_init (Vincent Neiger).
  * Fix tmp allocation size in _fmpz_set_str_basecase.
  * Fix rare arithmetic bug and memory leak in n_factor_ecm_select_curve.

* Other

  * Some corrections to the documentation.

2023-10-20 -- FLINT 3.0.0
-------------------------------------------------------------------------------

Merged libraries and reorganisation
...................................

* The following libraries have been merged into FLINT:

  * Arb 2.23 (arbitrary-precision ball arithmetic)
  * Calcium 0.4 (exact real and complex arithmetic)
  * Antic 0.2.5 (number fields, binary quadratic forms)

* Arb, Calcium and Antic will no longer be maintained as separate
  libraries. Users upgrading to FLINT 3.0 should ensure that they
  no longer link to the
  old Arb, Calcium or Antic library files or include any header files
  from those libraries which may be incompatible.

* The FLINT 3.0 API is largely backwards-compatible with FLINT 2.9,
  Arb 2.23, Calcium 0.4 and Antic 0.2.5, except for changes to
  rarely-used and internal functions documented below.
  However, the following changes to the handling of header files
  are likely to require (trivial) patches in many downstream codebases:

  * Header files belonging to Arb, Calcium and Antic
    now appear in the ``flint/`` subdirectory.
    For example, instead of ``#include "arb.h"``, it is necessary to
    ``#include "flint/arb.h"`` unless ``<INCLUDE_DIR>/flint`` has
    been added to the include path.

  * Most header files no longer include their implicit dependencies.
    For example, ``fmpz_poly.h`` no longer includes ``fmpz.h``.
    Code that used functions from the ``fmpz`` module but only
    included ``fmpz_poly.h`` may thus now need to include ``fmpz.h``
    explicitly. Likewise, many inclusions of system libraries like
    ``stdlib.h`` have been removed.

* The following people helped with the merge: Fredrik Johansson,
  Isuru Fernando, Albin Ahlbäck.

* FLINT 3.0 has a new build system based on Autotools,
  contributed by Albin Ahlbäck.
  Among other improvements, parallel builds are much faster
  and it is possible to build individual targets.
  Additional build system and CI improvements have been made by Marc Mezzarobba,
  Max Horn, Edgar Costa, Alex Best, Andreas Enge, and others.

* It is now necessary to run ``bootstrap.sh`` to generate the
  ``configure`` script in order to build FLINT from the git repository.

* Some ``configure`` options have changed: for
  example, ``--reentrant`` is now ``--enable-reentrant``.

* The root directory has been cleaned up by moving all source code
  into the ``src`` directory. This should not affect any users.

* The NTL interface has been moved to a single header file. The ``--with-ntl``
  build flag is now only needed to build the test code for this interface.

* The C++ interface (flintxx) has been removed. This interface is now
  maintained in the separate repository https://github.com/flintlib/flintxx
  (Edgar Costa).


Generic rings
..................

* The new ``gr`` module supports generic programming.
  It provides wrappers for most builtin FLINT types and allows
  constructing generic structures (polynomials, matrices, etc.) over
  arbitrary base rings.
  The following modules are available:

  * ``gr_generic`` (various generic algorithms)
  * ``gr_mat`` (matrices with generic elements)
  * ``gr_mpoly`` (multivariate polynomials with generic elements)
  * ``gr_poly`` (univariate polynomials with generic elements)
  * ``gr_special`` (special functions for generic elements)
  * ``gr_vec`` (vectors with generic elements)

  This feature is experimental: it is highly likely that some interfaces
  will change in a future FLINT release.

* There is also a Python wrapper (``flint_ctypes``) included with FLINT
  available in the ``src/python`` directory.
  Unlike other third-party FLINT wrappers available currently, this wrapper
  uses the ``gr`` interface to wrap (nearly) all FLINT types
  at once. This wrapper is not officially supported and will likely be
  deprecated in the future, but it can be useful for experimenting with FLINT.

* The generics system supports certain representations that do not have
  dedicated FLINT modules, for example 8-bit and 32-bit nmods.

Small-prime FFT
..................

* The new ``fft_small`` module implements FFTs modulo word-size
  primes and multiplication based on such FFTs.
  This module requires AVX2 or NEON vector instructions and will not be
  built on targets that do not support them.
  The small-prime  FFT speeds up the following functions for huge input,
  sometimes by a factor 2x to 10x:

  * ``flint_mpn_mul`` and variants, and indirectly any function based on
    FLINT's integer multiplication for large inputs. For example,
    ``fmpz_mul`` and ``arb_mul`` are faster, but ``fmpz_gcd`` is
    currently unaffected since it calls GMP.
  * ``nmod_poly_mul`` and variants, and indirectly any function based on
    ``nmod_poly`` multiplication.
  * ``fmpz_poly_mul`` and variants, and indirectly any function based on
    ``fmpz_poly`` multiplication.
  * Division functions for ``fmpz`` and ``arb``, which now use Newton
    iteration instead of calling GMP for huge input.
  * ``fmpz_mod`` arithmetic.
  * Radix conversion functions like ``fmpz_get_str``, ``fmpz_set_str``
    and ``arb_get_str``.

* The FFT was contributed by Daniel Schultz, with final integration
  work and adaptations for other FLINT functions (Newton iteration
  implementations, etc.) done by Fredrik Johansson.

Other changes
..................

* Changed the order of the ``alloc`` and ``length`` fields in ``arb_poly_t``,
  ``acb_poly_t`` and ``ca_poly_t`` to match the FLINT types.
* Added ``fmpzi`` division, norm and GCD functions (gcd_shortest by Daniel Schultz).
* Added an ``acf`` type for complex floating-point numbers.
* Added error handling to ``dirichlet_group_init``.
* Increased the prime factor limit in ``dirichlet_group_init`` from 1e12 to 1e16.
* Added ``arb_nonnegative_abs`` (Erik Postma).
* Fixed ``arb_pow`` for x just barely containing 0, y > 0 (Erik Postma).
* Improved precision handling in ``arb_gamma`` for huge input.
* Faster ``arb_contains_arf``, ``arb_overlaps``, ``arb_gt``, ``arb_lt``.
* Changed the argument order of ``_fmpz_mod_poly_mullow`` and
  ``_fmpz_mod_poly_div_series``.
* Changed the call signature of many ``_fmpz_mod_poly`` methods to
  take a context object as input instead of the raw modulus.
* Support test coverage reports (``--enable-coverage``).
* Added ``fmpz_poly_randtest_irreducible``.
* Improved tuning for various ``nmod_poly`` functions.
* Most Newton polynomial division and square root functions now use
  the Karp-Markstein algorithm.
* Replaced ``count_leading_zeros`` and ``count_trailing_zeros`` macros with ``flint_clz`` and ``flint_ctz``.
* Fixed ``nmod_poly_compose`` which was not using an asymptotically fast algorithm.
* Various functions in the ``nmod``, ``fmpz_mod``, ``fq`` modules and
  elsewhere have been rewritten to use algorithms in the generics
  module. In many cases the corresponding type-specific algorithm
  implementation has been removed entirely
  (for example, ``nmod_poly_divrem_newton`` no longer exists).
* Fixed ``fmpz_mod_poly_factor_squarefree``, ``nmod_poly_factor_squarefree``
  and ``fq_*_poly_factor_squarefree`` sometimes returning non-monic factors.
  Among other consequences, this could lead to functions like ``fq_poly_roots``
  returning incorrect roots
* Fixed several bugs in the ``fq_default`` modules (Tommy Hofmann).
* Fixed stack overflow in ``mpoly_divrem_ideal`` functions.
* Handle the the zero polynomial correctly in ``nmod_poly_shift_left`` (Vincent Neiger).
* Fixed handling of permutations in ``invert_cols`` matrix methods (Vincent Neiger).
* Added ``nmod_mat_permute_rows`` (Vincent Neiger).
* Fixed bug in ``mpoly_monomial_halves`` (Daniel Schultz).
* Fixed overflow bug in ``fmpz_mod_mpoly_divrem_ideal`` (Daniel Schultz).
* Optimized ``fmpz_addmul``, ``fmpz_addmul_ui``, ``fmpz_submul``, ``fmpz_submul_ui`` for small arguments.
* Fixed demotion bug in ``fmpz_addmul_si`` and ``fmpz_submul_si``.
* Optimized ``fmpq_cmp``, ``fmpq_cmp_ui``, ``fmpq_cmp_si``, ``fmpq_cmp_fmpz`` for small arguments.
* Optimized ``fmpz_poly_resultant_modular`` by using a tighter bound.
* Allow lll to work with rank deficient Z basis (Daniel Schultz).
* Added ``fmpq_mat_can_solve_dixon`` (William Hart).
* Inlined ``n_gcd`` (Albin Ahlbäck).
* Fixed fallback code for ``sub_ddmmss`` when given signed arguments.
* Many documentation fixes (Håvard Damm-Johnsen, Joel Dahne, Albin Ahlbäck, David Einstein, Alex Best, and others).
* Code simplifications (Vincent Neiger).
* Fixed several type signatures (Ricardo Buring).
* Fixed several memory leaks (Ricardo Buring).
* Fixed ``fmpz_poly_factor_squarefree`` crashing when given the zero polynomial.
* Added ``arb_minmax`` (Joel Dahne).
* Added ``_push_term_ffmpz`` functions to mpoly types (David Einstein).
* Added functions for printing nmod vectors (Vincent Neiger).
* Added ``nmod_poly_is_monic`` (Vincent Neiger).
* Fixed threaded Arb functions to use the thread pool (Albin Ahlbäck).
* Removed nmod_poly mpz functions (Ricardo Buring).
* Fixed file handling in qsieve (Michiel de Wilde, Oscar Benjamin).
* Free memory in case of failure in ``fq_zech_ctx_init`` (Claus Fieker).
* Fixed corrupted output in ``fmpz_or``.
* Added several ``nmod_poly_mat`` utility functions (Vincent Neiger).

List of additions
.................

* FLINT 3.0 includes all functions in FLINT 2.9, Arb 2.23, Calcium 0.4
  and Antic 0.2.5 except those listed under "list of removals".
  On top of this, the following functions have been added.
  This list is incomplete; many internal functions and functions
  starting with an underscore have been omitted.
* ``mpn_mul_default_mpn_ctx``, ``_nmod_poly_mul_mid_default_mpn_ctx``, ``_fmpz_poly_mul_mid_default_mpn_ctx`` and many internal functions in the new ``fft_small`` module
* ``acb_poly_nth_derivative, arb_div_arf_newton, arb_div_newton, arb_fmpz_divapprox, arb_nint, arb_poly_nth_derivative, arb_rsqrt_arf, arb_rsqrt_arf_newton, arb_sqrt_arf_newton, arb_sqrt_newton, arb_trunc, arb_minmax``
* ``ca_set_fmpzi``
* ``flint_aligned_alloc, flint_aligned_free``
* ``flint_get_num_available_threads``
* ``flint_mpn_add_inplace_c, flint_mpn_cmp_ui_2exp, flint_mpn_mul_large, flint_mpn_nbits``
* ``fmpz_get_str_bsplit_threaded``
* ``fmpz_mat_equal_col, fmpz_mat_equal_row, fmpz_neg_ui_array``
* ``fmpz_poly_randtest_irreducible``
* ``fmpz_poly_q_evaluate_fmpq, fmpz_poly_q_scalar_div_fmpq, fmpz_poly_q_scalar_div_fmpz, fmpz_poly_q_scalar_mul_fmpq, fmpz_poly_q_scalar_mul_fmpz``
* ``fmpz_ui_pow_ui``
* ``fmpzi_set_qqbar``
* ``get_default_mpn_ctx``
* ``gr_abs, gr_acos, gr_acos_pi, gr_acosh, gr_acot, gr_acot_pi, gr_acoth, gr_acsc, gr_acsc_pi, gr_acsch, gr_add, gr_add_fmpq, gr_add_fmpz, gr_add_other, gr_add_si, gr_add_ui, gr_addmul, gr_addmul_fmpq, gr_addmul_fmpz, gr_addmul_other, gr_addmul_si, gr_addmul_ui, gr_agm, gr_agm1, gr_airy, gr_airy_ai, gr_airy_ai_prime, gr_airy_ai_prime_zero, gr_airy_ai_zero, gr_airy_bi, gr_airy_bi_prime, gr_airy_bi_prime_zero, gr_airy_bi_zero, gr_asec, gr_asec_pi, gr_asech, gr_asin, gr_asin_pi, gr_asinh, gr_atan, gr_atan2, gr_atan_pi, gr_atanh, gr_barnes_g, gr_bellnum_fmpz, gr_bellnum_ui, gr_bellnum_vec, gr_bernoulli_fmpz, gr_bernoulli_ui, gr_bernoulli_vec, gr_bernpoly_ui, gr_bessel_i, gr_bessel_i_scaled, gr_bessel_j, gr_bessel_j_y, gr_bessel_k, gr_bessel_k_scaled, gr_bessel_y, gr_beta, gr_beta_lower, gr_bin, gr_bin_ui, gr_bin_ui_vec, gr_bin_uiui, gr_bin_vec, gr_carlson_rc, gr_carlson_rd, gr_carlson_rf, gr_carlson_rg, gr_carlson_rj, gr_catalan, gr_ceil, gr_chebyshev_t, gr_chebyshev_t_fmpz, gr_chebyshev_u, gr_chebyshev_u_fmpz, gr_clear, gr_cmp, gr_cmp_other, gr_cmpabs, gr_cmpabs_other, gr_conj, gr_cos, gr_cos_integral, gr_cos_pi, gr_cosh, gr_cosh_integral, gr_cot, gr_cot_pi, gr_coth, gr_coulomb, gr_coulomb_f, gr_coulomb_g, gr_coulomb_hneg, gr_coulomb_hpos, gr_csc, gr_csc_pi, gr_csch, gr_csgn, gr_ctx_ca_get_option, gr_ctx_ca_set_option, gr_ctx_clear, gr_ctx_cmp_coercion, gr_ctx_data_as_ptr, gr_ctx_data_ptr, gr_ctx_fmpz_mod_set_primality, gr_ctx_fq_degree, gr_ctx_fq_order, gr_ctx_fq_prime, gr_ctx_get_real_prec, gr_ctx_get_str, gr_ctx_has_real_prec, gr_ctx_init_complex_acb, gr_ctx_init_complex_algebraic_ca, gr_ctx_init_complex_ca, gr_ctx_init_complex_float_acf, gr_ctx_init_complex_qqbar, gr_ctx_init_dirichlet_group, gr_ctx_init_fmpq, gr_ctx_init_fmpz, gr_ctx_init_fmpz_mod, gr_ctx_init_fmpz_poly, gr_ctx_init_fmpzi, gr_ctx_init_fq, gr_ctx_init_fq_nmod, gr_ctx_init_fq_zech, gr_ctx_init_gr_series, gr_ctx_init_gr_series_mod, gr_ctx_init_matrix_domain, gr_ctx_init_matrix_ring, gr_ctx_init_matrix_space, gr_ctx_init_gr_mpoly, gr_ctx_init_nf, gr_ctx_init_nf_fmpq_poly, gr_ctx_init_nmod, gr_ctx_init_nmod8, gr_ctx_init_nmod32, gr_ctx_init_perm, gr_ctx_init_gr_poly, gr_ctx_init_psl2z, gr_ctx_init_random, gr_ctx_init_real_algebraic_ca, gr_ctx_init_real_arb, gr_ctx_init_real_ca, gr_ctx_init_real_float_arf, gr_ctx_init_real_qqbar, gr_ctx_init_vector_gr_vec, gr_ctx_init_vector_space_gr_vec, gr_ctx_is_algebraically_closed, gr_ctx_is_canonical, gr_ctx_is_commutative_ring, gr_ctx_is_exact, gr_ctx_is_field, gr_ctx_is_finite, gr_ctx_is_finite_characteristic, gr_ctx_is_integral_domain, gr_ctx_is_multiplicative_group, gr_ctx_is_ordered_ring, gr_ctx_is_ring, gr_ctx_is_threadsafe, gr_ctx_is_unique_factorization_domain, gr_ctx_matrix_is_fixed_size, gr_ctx_print, gr_ctx_println, gr_ctx_set_real_prec, gr_ctx_sizeof_ctx, gr_ctx_sizeof_elem, gr_ctx_vector_gr_vec_is_fixed_size, gr_ctx_write, gr_dedekind_eta, gr_dedekind_eta_q, gr_digamma, gr_dilog, gr_dirichlet_beta, gr_dirichlet_chi_fmpz, gr_dirichlet_chi_vec, gr_dirichlet_eta, gr_dirichlet_hardy_theta, gr_dirichlet_hardy_z, gr_dirichlet_l, gr_div, gr_div_fmpq, gr_div_fmpz, gr_div_other, gr_div_si, gr_div_ui, gr_divexact, gr_divexact_fmpq, gr_divexact_fmpz, gr_divexact_other, gr_divexact_si, gr_divexact_ui, gr_divides, gr_dot_other, gr_doublefac, gr_doublefac_ui, gr_eisenstein_e, gr_eisenstein_g, gr_eisenstein_g_vec, gr_elliptic_e, gr_elliptic_e_inc, gr_elliptic_f, gr_elliptic_invariants, gr_elliptic_k, gr_elliptic_pi, gr_elliptic_pi_inc, gr_elliptic_roots, gr_equal, gr_erf, gr_erfc, gr_erfcinv, gr_erfcx, gr_erfi, gr_erfinv, gr_euclidean_div, gr_euclidean_divrem, gr_euclidean_rem, gr_euler, gr_eulernum_fmpz, gr_eulernum_ui, gr_eulernum_vec, gr_eulerpoly_ui, gr_evaluate_fmpz_mpoly_iter, gr_exp, gr_exp10, gr_exp2, gr_exp_integral, gr_exp_integral_ei, gr_exp_pi_i, gr_expm1, gr_fac, gr_fac_fmpz, gr_fac_ui, gr_fac_vec, gr_factor, gr_falling, gr_falling_ui, gr_fib_fmpz, gr_fib_ui, gr_fib_vec, gr_floor, gr_fmms, gr_fmpz_mpoly_evaluate, gr_fmpz_mpoly_evaluate_horner, gr_fmpz_poly_evaluate, gr_fmpz_poly_evaluate_horner, gr_fmpz_poly_evaluate_rectangular, gr_fq_frobenius, gr_fq_is_primitive, gr_fq_multiplicative_order, gr_fq_norm, gr_fq_pth_root, gr_fq_trace, gr_fresnel, gr_fresnel_c, gr_fresnel_s, gr_gamma, gr_gamma_fmpq, gr_gamma_fmpz, gr_gamma_lower, gr_gamma_upper, gr_gcd, gr_gegenbauer_c, gr_gen, gr_generic_acot, gr_generic_acoth, gr_generic_acsc, gr_generic_acsch, gr_generic_add_fmpq, gr_generic_add_fmpz, gr_generic_add_other, gr_generic_add_si, gr_generic_add_ui, gr_generic_addmul, gr_generic_addmul_fmpq, gr_generic_addmul_fmpz, gr_generic_addmul_other, gr_generic_addmul_si, gr_generic_addmul_ui, gr_generic_asec, gr_generic_asech, gr_generic_asin, gr_generic_asinh, gr_generic_atan, gr_generic_atanh, gr_generic_bellnum_fmpz, gr_generic_bellnum_ui, gr_generic_bellnum_vec, gr_generic_bernoulli_fmpz, gr_generic_bernoulli_ui, gr_generic_bernoulli_vec, gr_generic_beta, gr_generic_bin, gr_generic_bin_ui, gr_generic_bin_ui_vec, gr_generic_bin_uiui, gr_generic_bin_vec, gr_generic_chebyshev_t2_fmpz, gr_generic_chebyshev_t_fmpz, gr_generic_chebyshev_u2_fmpz, gr_generic_chebyshev_u_fmpz, gr_generic_cmp, gr_generic_cmp_other, gr_generic_cmpabs, gr_generic_cmpabs_other, gr_generic_cos, gr_generic_ctx_clear, gr_generic_ctx_predicate, gr_generic_ctx_predicate_false, gr_generic_ctx_predicate_true, gr_generic_div_fmpq, gr_generic_div_fmpz, gr_generic_div_other, gr_generic_div_si, gr_generic_div_ui, gr_generic_divexact, gr_generic_doublefac, gr_generic_doublefac_ui, gr_generic_erfcx, gr_generic_eulernum_fmpz, gr_generic_eulernum_ui, gr_generic_eulernum_vec, gr_generic_exp, gr_generic_exp10, gr_generic_exp2, gr_generic_expm1, gr_generic_fac, gr_generic_fac_fmpz, gr_generic_fac_ui, gr_generic_fac_vec, gr_generic_falling, gr_generic_falling_ui, gr_generic_fib2_fmpz, gr_generic_fib_fmpz, gr_generic_fib_ui, gr_generic_fib_vec, gr_generic_get_fmpz_2exp_fmpz, gr_generic_harmonic, gr_generic_harmonic_ui, gr_generic_hilbert_class_poly, gr_generic_inv, gr_generic_is_invertible, gr_generic_is_neg_one, gr_generic_is_one, gr_generic_is_square, gr_generic_is_zero, gr_generic_log, gr_generic_log10, gr_generic_log1p, gr_generic_log2, gr_generic_mul_2exp_fmpz, gr_generic_mul_2exp_si, gr_generic_mul_fmpq, gr_generic_mul_fmpz, gr_generic_mul_other, gr_generic_mul_si, gr_generic_mul_two, gr_generic_mul_ui, gr_generic_mul_ui_via_ZZ, gr_generic_neg_one, gr_generic_other_add, gr_generic_other_add_vec, gr_generic_other_div, gr_generic_other_div_vec, gr_generic_other_divexact_vec, gr_generic_other_mul, gr_generic_other_mul_vec, gr_generic_other_pow, gr_generic_other_pow_vec, gr_generic_other_sub, gr_generic_other_sub_vec, gr_generic_partitions_fmpz, gr_generic_partitions_ui, gr_generic_partitions_vec, gr_generic_pow_fmpq, gr_generic_pow_fmpz, gr_generic_pow_fmpz_binexp, gr_generic_pow_other, gr_generic_pow_si, gr_generic_pow_ui, gr_generic_pow_ui_binexp, gr_generic_randtest_not_zero, gr_generic_rfac, gr_generic_rfac_fmpz, gr_generic_rfac_ui, gr_generic_rfac_vec, gr_generic_rising, gr_generic_rising_ui, gr_generic_rsqrt, gr_generic_scalar_add_vec, gr_generic_scalar_div_vec, gr_generic_scalar_divexact_vec, gr_generic_scalar_mul_vec, gr_generic_scalar_other_add_vec, gr_generic_scalar_other_div_vec, gr_generic_scalar_other_divexact_vec, gr_generic_scalar_other_mul_vec, gr_generic_scalar_other_pow_vec, gr_generic_scalar_other_sub_vec, gr_generic_scalar_pow_vec, gr_generic_scalar_sub_vec, gr_generic_set_fmpq, gr_generic_set_fmpz_2exp_fmpz, gr_generic_set_other, gr_generic_set_shallow, gr_generic_sin, gr_generic_sin_cos, gr_generic_sqr, gr_generic_sqrt, gr_generic_stirling_s1_ui_vec, gr_generic_stirling_s1_uiui, gr_generic_stirling_s1u_ui_vec, gr_generic_stirling_s1u_uiui, gr_generic_stirling_s2_ui_vec, gr_generic_stirling_s2_uiui, gr_generic_sub_fmpq, gr_generic_sub_fmpz, gr_generic_sub_other, gr_generic_sub_si, gr_generic_sub_ui, gr_generic_submul, gr_generic_submul_fmpq, gr_generic_submul_fmpz, gr_generic_submul_other, gr_generic_submul_si, gr_generic_submul_ui, gr_generic_tan, gr_generic_vec_add, gr_generic_vec_add_other, gr_generic_vec_add_scalar, gr_generic_vec_add_scalar_fmpq, gr_generic_vec_add_scalar_fmpz, gr_generic_vec_add_scalar_other, gr_generic_vec_add_scalar_si, gr_generic_vec_add_scalar_ui, gr_generic_vec_clear, gr_generic_vec_div, gr_generic_vec_div_other, gr_generic_vec_div_scalar, gr_generic_vec_div_scalar_fmpq, gr_generic_vec_div_scalar_fmpz, gr_generic_vec_div_scalar_other, gr_generic_vec_div_scalar_si, gr_generic_vec_div_scalar_ui, gr_generic_vec_divexact, gr_generic_vec_divexact_other, gr_generic_vec_divexact_scalar, gr_generic_vec_divexact_scalar_fmpq, gr_generic_vec_divexact_scalar_fmpz, gr_generic_vec_divexact_scalar_other, gr_generic_vec_divexact_scalar_si, gr_generic_vec_divexact_scalar_ui, gr_generic_vec_dot, gr_generic_vec_dot_fmpz, gr_generic_vec_dot_rev, gr_generic_vec_dot_si, gr_generic_vec_dot_ui, gr_generic_vec_equal, gr_generic_vec_init, gr_generic_vec_is_zero, gr_generic_vec_mul, gr_generic_vec_mul_other, gr_generic_vec_mul_scalar, gr_generic_vec_mul_scalar_2exp_si, gr_generic_vec_mul_scalar_fmpq, gr_generic_vec_mul_scalar_fmpz, gr_generic_vec_mul_scalar_other, gr_generic_vec_mul_scalar_si, gr_generic_vec_mul_scalar_ui, gr_generic_vec_neg, gr_generic_vec_normalise, gr_generic_vec_normalise_weak, gr_generic_vec_pow, gr_generic_vec_pow_other, gr_generic_vec_pow_scalar, gr_generic_vec_pow_scalar_fmpq, gr_generic_vec_pow_scalar_fmpz, gr_generic_vec_pow_scalar_other, gr_generic_vec_pow_scalar_si, gr_generic_vec_pow_scalar_ui, gr_generic_vec_reciprocals, gr_generic_vec_scalar_addmul, gr_generic_vec_scalar_addmul_si, gr_generic_vec_scalar_submul, gr_generic_vec_scalar_submul_si, gr_generic_vec_set, gr_generic_vec_set_powers, gr_generic_vec_sub, gr_generic_vec_sub_other, gr_generic_vec_sub_scalar, gr_generic_vec_sub_scalar_fmpq, gr_generic_vec_sub_scalar_fmpz, gr_generic_vec_sub_scalar_other, gr_generic_vec_sub_scalar_si, gr_generic_vec_sub_scalar_ui, gr_generic_vec_swap, gr_generic_vec_zero, gr_generic_write_n, gr_get_d, gr_get_fmpq, gr_get_fmpz, gr_get_fmpz_2exp_fmpz, gr_get_si, gr_get_str, gr_get_str_n, gr_get_ui, gr_glaisher, gr_harmonic, gr_harmonic_ui, gr_heap_clear, gr_heap_clear_vec, gr_heap_init, gr_heap_init_vec, gr_hermite_h, gr_hilbert_class_poly, gr_hurwitz_zeta, gr_hypgeom_0f1, gr_hypgeom_1f1, gr_hypgeom_2f1, gr_hypgeom_pfq, gr_hypgeom_u, gr_i, gr_im, gr_init, gr_inv, gr_is_invertible, gr_is_neg_one, gr_is_one, gr_is_square, gr_is_zero, gr_jacobi_p, gr_jacobi_theta, gr_jacobi_theta_1, gr_jacobi_theta_2, gr_jacobi_theta_3, gr_jacobi_theta_4, gr_khinchin, gr_laguerre_l, gr_lambertw, gr_lambertw_fmpz, gr_lcm, gr_legendre_p, gr_legendre_p_root_ui, gr_legendre_q, gr_lerch_phi, gr_lgamma, gr_log, gr_log10, gr_log1p, gr_log2, gr_log_barnes_g, gr_log_integral, gr_log_pi_i, gr_mat_add, gr_mat_add_scalar, gr_mat_addmul_scalar, gr_mat_adjugate, gr_mat_adjugate_charpoly, gr_mat_adjugate_cofactor, gr_mat_apply_row_similarity, gr_mat_charpoly, gr_mat_charpoly_berkowitz, gr_mat_charpoly_danilevsky, gr_mat_charpoly_faddeev, gr_mat_charpoly_faddeev_bsgs, gr_mat_charpoly_from_hessenberg, gr_mat_charpoly_gauss, gr_mat_charpoly_householder, gr_mat_clear, gr_mat_concat_horizontal, gr_mat_concat_vertical, gr_mat_det, gr_mat_det_berkowitz, gr_mat_det_cofactor, gr_mat_det_fflu, gr_mat_det_generic, gr_mat_det_generic_field, gr_mat_det_generic_integral_domain, gr_mat_det_lu, gr_mat_diag_mul, gr_mat_diagonalization, gr_mat_diagonalization_generic, gr_mat_diagonalization_precomp, gr_mat_div_scalar, gr_mat_eigenvalues, gr_mat_eigenvalues_other, gr_mat_entry_ptr, gr_mat_entry_srcptr, gr_mat_equal, gr_mat_exp, gr_mat_exp_jordan, gr_mat_fflu, gr_mat_find_nonzero_pivot, gr_mat_find_nonzero_pivot_generic, gr_mat_find_nonzero_pivot_large_abs, gr_mat_gr_poly_evaluate, gr_mat_hadamard, gr_mat_hessenberg, gr_mat_hessenberg_gauss, gr_mat_hessenberg_householder, gr_mat_hilbert, gr_mat_init, gr_mat_init_set, gr_mat_inv, gr_mat_invert_cols, gr_mat_invert_rows, gr_mat_is_diagonal, gr_mat_is_empty, gr_mat_is_hessenberg, gr_mat_is_lower_triangular, gr_mat_is_neg_one, gr_mat_is_one, gr_mat_is_scalar, gr_mat_is_square, gr_mat_is_upper_triangular, gr_mat_is_zero, gr_mat_jordan_blocks, gr_mat_jordan_form, gr_mat_jordan_transformation, gr_mat_log, gr_mat_log_jordan, gr_mat_lu, gr_mat_lu_classical, gr_mat_lu_recursive, gr_mat_minpoly_field, gr_mat_mul, gr_mat_mul_classical, gr_mat_mul_diag, gr_mat_mul_generic, gr_mat_mul_scalar, gr_mat_mul_strassen, gr_mat_neg, gr_mat_nonsingular_solve, gr_mat_nonsingular_solve_den, gr_mat_nonsingular_solve_den_fflu, gr_mat_nonsingular_solve_fflu, gr_mat_nonsingular_solve_fflu_precomp, gr_mat_nonsingular_solve_lu, gr_mat_nonsingular_solve_lu_precomp, gr_mat_nonsingular_solve_tril, gr_mat_nonsingular_solve_tril_classical, gr_mat_nonsingular_solve_tril_recursive, gr_mat_nonsingular_solve_triu, gr_mat_nonsingular_solve_triu_classical, gr_mat_nonsingular_solve_triu_recursive, gr_mat_nullspace, gr_mat_one, gr_mat_ones, gr_mat_pascal, gr_mat_print, gr_mat_randops, gr_mat_randpermdiag, gr_mat_randrank, gr_mat_randtest, gr_mat_rank, gr_mat_rank_fflu, gr_mat_rank_lu, gr_mat_reduce_row, gr_mat_rref, gr_mat_rref_den, gr_mat_rref_den_fflu, gr_mat_rref_fflu, gr_mat_rref_lu, gr_mat_set, gr_mat_set_fmpq, gr_mat_set_fmpq_mat, gr_mat_set_fmpz, gr_mat_set_fmpz_mat, gr_mat_set_jordan_blocks, gr_mat_set_scalar, gr_mat_set_si, gr_mat_set_ui, gr_mat_solve_field, gr_mat_sqr, gr_mat_stirling, gr_mat_sub, gr_mat_sub_scalar, gr_mat_submul_scalar, gr_mat_swap, gr_mat_swap_cols, gr_mat_swap_entrywise, gr_mat_swap_rows, gr_mat_trace, gr_mat_trace_prod2, gr_mat_transpose, gr_mat_transpose_resize, gr_mat_window_clear, gr_mat_window_init, gr_mat_write, gr_mat_zero, gr_method_tab_init, gr_modular_delta, gr_modular_j, gr_modular_lambda, gr_mpoly_add, gr_mpoly_assert_canonical, gr_mpoly_clear, gr_mpoly_combine_like_terms, gr_mpoly_equal, gr_mpoly_fit_bits, gr_mpoly_fit_length, gr_mpoly_fit_length_fit_bits, gr_mpoly_fit_length_reset_bits, gr_mpoly_gen, gr_mpoly_get_coeff_scalar_fmpz, gr_mpoly_get_coeff_scalar_ui, gr_mpoly_init, gr_mpoly_init2, gr_mpoly_init3, gr_mpoly_is_canonical, gr_mpoly_is_gen, gr_mpoly_is_one, gr_mpoly_is_zero, gr_mpoly_mul, gr_mpoly_mul_fmpq, gr_mpoly_mul_fmpz, gr_mpoly_mul_johnson, gr_mpoly_mul_monomial, gr_mpoly_mul_scalar, gr_mpoly_mul_si, gr_mpoly_mul_ui, gr_mpoly_neg, gr_mpoly_one, gr_mpoly_print_pretty, gr_mpoly_push_term_scalar_fmpz, gr_mpoly_push_term_scalar_ui, gr_mpoly_randtest_bits, gr_mpoly_randtest_bound, gr_mpoly_set, gr_mpoly_set_coeff_fmpq_fmpz, gr_mpoly_set_coeff_fmpq_ui, gr_mpoly_set_coeff_fmpz_fmpz, gr_mpoly_set_coeff_fmpz_ui, gr_mpoly_set_coeff_scalar_fmpz, gr_mpoly_set_coeff_scalar_ui, gr_mpoly_set_coeff_si_fmpz, gr_mpoly_set_coeff_si_ui, gr_mpoly_set_coeff_ui_fmpz, gr_mpoly_set_coeff_ui_ui, gr_mpoly_set_fmpq, gr_mpoly_set_fmpz, gr_mpoly_set_scalar, gr_mpoly_set_si, gr_mpoly_set_ui, gr_mpoly_sort_terms, gr_mpoly_sub, gr_mpoly_swap, gr_mpoly_write_pretty, gr_mpoly_zero, gr_mul, gr_mul_2exp_fmpz, gr_mul_2exp_si, gr_mul_fmpq, gr_mul_fmpz, gr_mul_other, gr_mul_si, gr_mul_two, gr_mul_ui, gr_neg, gr_neg_one, gr_nint, gr_not_equal, gr_not_implemented, gr_not_in_domain, gr_one, gr_other_add, gr_other_div, gr_other_divexact, gr_other_mul, gr_other_pow, gr_other_sub, gr_partitions_fmpz, gr_partitions_ui, gr_partitions_vec, gr_pi, gr_poly_acos_series, gr_poly_acosh_series, gr_poly_add, gr_poly_add_series, gr_poly_asin_series, gr_poly_asinh_series, gr_poly_atan_series, gr_poly_atanh_series, gr_poly_clear, gr_poly_compose, gr_poly_compose_divconquer, gr_poly_compose_horner, gr_poly_compose_series, gr_poly_compose_series_brent_kung, gr_poly_compose_series_divconquer, gr_poly_compose_series_horner, gr_poly_derivative, gr_poly_div, gr_poly_div_basecase, gr_poly_div_divconquer, gr_poly_div_newton, gr_poly_div_series, gr_poly_div_series_basecase, gr_poly_div_series_invmul, gr_poly_div_series_newton, gr_poly_divrem, gr_poly_divrem_basecase, gr_poly_divrem_divconquer, gr_poly_divrem_newton, gr_poly_entry_ptr, gr_poly_equal, gr_poly_evaluate, gr_poly_evaluate_horner, gr_poly_evaluate_other, gr_poly_evaluate_other_horner, gr_poly_evaluate_other_rectangular, gr_poly_evaluate_rectangular, gr_poly_evaluate_vec_fast, gr_poly_evaluate_vec_iter, gr_poly_exp_series, gr_poly_exp_series_basecase, gr_poly_exp_series_basecase_mul, gr_poly_exp_series_newton, gr_poly_factor_squarefree, gr_poly_fit_length, gr_poly_gcd, gr_poly_gcd_euclidean, gr_poly_gcd_hgcd, gr_poly_gen, gr_poly_get_coeff_scalar, gr_poly_init, gr_poly_init2, gr_poly_integral, gr_poly_inv, gr_poly_inv_series, gr_poly_inv_series_basecase, gr_poly_inv_series_newton, gr_poly_is_gen, gr_poly_is_monic, gr_poly_is_one, gr_poly_is_zero, gr_poly_length, gr_poly_log1p_series, gr_poly_log_series, gr_poly_make_monic, gr_poly_mul, gr_poly_mul_scalar, gr_poly_mullow, gr_poly_neg, gr_poly_neg_one, gr_poly_nth_derivative, gr_poly_one, gr_poly_pow_fmpz, gr_poly_pow_series_fmpq_recurrence, gr_poly_pow_series_ui, gr_poly_pow_series_ui_binexp, gr_poly_pow_ui, gr_poly_pow_ui_binexp, gr_poly_print, gr_poly_randtest, gr_poly_rem, gr_poly_resultant, gr_poly_resultant_euclidean, gr_poly_resultant_hgcd, gr_poly_resultant_small, gr_poly_resultant_sylvester, gr_poly_reverse, gr_poly_roots, gr_poly_roots_other, gr_poly_rsqrt_series, gr_poly_rsqrt_series_basecase, gr_poly_rsqrt_series_miller, gr_poly_rsqrt_series_newton, gr_poly_set, gr_poly_set_coeff_fmpq, gr_poly_set_coeff_fmpz, gr_poly_set_coeff_scalar, gr_poly_set_coeff_si, gr_poly_set_coeff_ui, gr_poly_set_fmpq, gr_poly_set_fmpq_poly, gr_poly_set_fmpz, gr_poly_set_fmpz_poly, gr_poly_set_gr_poly_other, gr_poly_set_scalar, gr_poly_set_si, gr_poly_set_ui, gr_poly_shift_left, gr_poly_shift_right, gr_poly_sin_cos_series_basecase, gr_poly_sin_cos_series_tangent, gr_poly_sqrt_series, gr_poly_sqrt_series_basecase, gr_poly_sqrt_series_miller, gr_poly_sqrt_series_newton, gr_poly_squarefree_part, gr_poly_sub, gr_poly_sub_series, gr_poly_swap, gr_poly_tan_series, gr_poly_tan_series_basecase, gr_poly_tan_series_newton, gr_poly_taylor_shift, gr_poly_taylor_shift_convolution, gr_poly_taylor_shift_divconquer, gr_poly_taylor_shift_horner, gr_poly_truncate, gr_poly_write, gr_poly_xgcd_euclidean, gr_poly_xgcd_hgcd, gr_poly_zero, gr_polygamma, gr_polylog, gr_pow, gr_pow_fmpq, gr_pow_fmpz, gr_pow_other, gr_pow_si, gr_pow_ui, gr_print, gr_println, gr_randtest, gr_randtest_not_zero, gr_randtest_small, gr_re, gr_rfac, gr_rfac_fmpz, gr_rfac_ui, gr_rfac_vec, gr_rgamma, gr_riemann_xi, gr_rising, gr_rising_ui, gr_rising_ui_forward, gr_rsqrt, gr_sec, gr_sec_pi, gr_sech, gr_series_acos, gr_series_acosh, gr_series_add, gr_series_agm1, gr_series_airy, gr_series_airy_ai, gr_series_airy_ai_prime, gr_series_airy_bi, gr_series_airy_bi_prime, gr_series_asin, gr_series_asinh, gr_series_atan, gr_series_atanh, gr_series_beta_lower, gr_series_clear, gr_series_cos_integral, gr_series_cosh_integral, gr_series_digamma, gr_series_dirichlet_hardy_theta, gr_series_dirichlet_hardy_z, gr_series_dirichlet_l, gr_series_div, gr_series_elliptic_k, gr_series_equal, gr_series_erf, gr_series_erfc, gr_series_erfi, gr_series_exp, gr_series_exp_integral_ei, gr_series_fresnel, gr_series_fresnel_c, gr_series_fresnel_s, gr_series_gamma, gr_series_gamma_lower, gr_series_gamma_upper, gr_series_gen, gr_series_hurwitz_zeta, gr_series_hypgeom_pfq, gr_series_init, gr_series_inv, gr_series_is_one, gr_series_is_zero, gr_series_jacobi_theta, gr_series_jacobi_theta_1, gr_series_jacobi_theta_2, gr_series_jacobi_theta_3, gr_series_jacobi_theta_4, gr_series_lgamma, gr_series_log, gr_series_log_integral, gr_series_make_exact, gr_series_mul, gr_series_neg, gr_series_one, gr_series_polylog, gr_series_randtest, gr_series_rgamma, gr_series_rsqrt, gr_series_set, gr_series_set_fmpq, gr_series_set_fmpz, gr_series_set_gr_poly, gr_series_set_scalar, gr_series_set_si, gr_series_set_ui, gr_series_sin_integral, gr_series_sinh_integral, gr_series_sqrt, gr_series_sub, gr_series_swap, gr_series_tan, gr_series_weierstrass_p, gr_series_write, gr_series_zero, gr_set, gr_set_d, gr_set_fmpq, gr_set_fmpz, gr_set_fmpz_2exp_fmpz, gr_set_other, gr_set_shallow, gr_set_si, gr_set_str, gr_set_ui, gr_sgn, gr_sin, gr_sin_cos, gr_sin_cos_pi, gr_sin_integral, gr_sin_pi, gr_sinc, gr_sinc_pi, gr_sinh, gr_sinh_cosh, gr_sinh_integral, gr_spherical_y_si, gr_sqr, gr_sqrt, gr_stieltjes, gr_stirling_s1_ui_vec, gr_stirling_s1_uiui, gr_stirling_s1u_ui_vec, gr_stirling_s1u_uiui, gr_stirling_s2_ui_vec, gr_stirling_s2_uiui, gr_stream_init_file, gr_stream_init_str, gr_stream_write, gr_stream_write_fmpz, gr_stream_write_free, gr_stream_write_si, gr_stream_write_ui, gr_sub, gr_sub_fmpq, gr_sub_fmpz, gr_sub_other, gr_sub_si, gr_sub_ui, gr_submul, gr_submul_fmpq, gr_submul_fmpz, gr_submul_other, gr_submul_si, gr_submul_ui, gr_swap, gr_swap2, gr_tan, gr_tan_pi, gr_tanh, gr_test_add_aliasing, gr_test_add_associative, gr_test_add_commutative, gr_test_add_type_variants, gr_test_addmul_submul, gr_test_addmul_type_variants, gr_test_binary_op_aliasing, gr_test_binary_op_associative, gr_test_binary_op_commutative, gr_test_binary_op_left_distributive, gr_test_binary_op_right_distributive, gr_test_binary_op_type_variants, gr_test_complex_parts, gr_test_div_right_distributive, gr_test_div_then_mul, gr_test_div_type_variants, gr_test_divexact, gr_test_divexact_type_variants, gr_test_equal, gr_test_field, gr_test_get_fmpq, gr_test_get_fmpz, gr_test_get_fmpz_2exp_fmpz, gr_test_get_si, gr_test_get_ui, gr_test_init_clear, gr_test_integral_domain, gr_test_inv_involution, gr_test_inv_multiplication, gr_test_iter, gr_test_mat_mul_classical_associative, gr_test_mul_2exp_fmpz, gr_test_mul_2exp_si, gr_test_mul_aliasing, gr_test_mul_associative, gr_test_mul_commutative, gr_test_mul_left_distributive, gr_test_mul_right_distributive, gr_test_mul_then_div, gr_test_mul_type_variants, gr_test_multiplicative_group, gr_test_neg, gr_test_one, gr_test_ordered_ring_cmp, gr_test_ordered_ring_cmpabs, gr_test_pow_fmpz_exponent_addition, gr_test_pow_ui_aliasing, gr_test_pow_ui_base_multiplication, gr_test_pow_ui_base_scalar_multiplication, gr_test_pow_ui_exponent_addition, gr_test_randtest_not_zero, gr_test_ring, gr_test_rsqrt, gr_test_set_fmpq, gr_test_set_fmpz, gr_test_set_si, gr_test_set_ui, gr_test_sqrt, gr_test_sub_aliasing, gr_test_sub_equal_neg_add, gr_test_sub_type_variants, gr_test_submul_type_variants, gr_test_swap, gr_test_vec_add, gr_test_vec_binary_op, gr_test_vec_div, gr_test_vec_divexact, gr_test_vec_dot, gr_test_vec_mul, gr_test_vec_pow, gr_test_vec_sub, gr_test_zero_one, gr_trunc, gr_vec_append, gr_vec_clear, gr_vec_entry_ptr, gr_vec_entry_srcptr, gr_vec_fit_length, gr_vec_init, gr_vec_length, gr_vec_print, gr_vec_set, gr_vec_set_length, gr_vec_write, gr_weierstrass_p, gr_weierstrass_p_inv, gr_weierstrass_p_prime, gr_weierstrass_sigma, gr_weierstrass_zeta, gr_write, gr_write_n, gr_zero, gr_zeta, gr_zeta_nzeros, gr_zeta_ui, gr_zeta_zero, gr_zeta_zero_vec``
* ``gr_pos_inf, gr_neg_inf, gr_uinf, gr_undefined, gr_unknown, gr_arg, gr_ctx_init_complex_extended_ca, gr_poly_divexact_basecase_bidirectional, gr_poly_divexact_bidirectional, gr_poly_divexact_basecase, gr_poly_is_scalar, gr_poly_div_series_divconquer, gr_poly_divexact_series_basecase``
* ``nmod_mat_fprint_pretty, nmod_mat_print, nmod_mat_fprint, nmod_poly_is_monic``
* ``nmod_poly_mat_set_trunc, nmod_poly_mat_truncate, nmod_poly_mat_shift_left, nmod_poly_mat_shift_right, nmod_poly_mat_get_coeff_mat, nmod_poly_mat_set_coeff_mat, nmod_poly_mat_set_nmod_mat, nmod_poly_mat_equal_nmod_mat, nmod_poly_mat_degree``
* ``qqbar_set_fmpzi``
* ``fmpq_mpoly_push_term_fmpq_ffmpz, fmpq_mpoly_push_term_fmpz_ffmpz, fmpq_mpoly_push_term_ui_ffmpz, fmpq_mpoly_push_term_si_ffmpz, fmpz_mod_mpoly_push_term_fmpz_ffmpz, fmpz_mod_mpoly_push_term_ui_ffmpz, fmpz_mod_mpoly_push_term_si_ffmpz, fmpz_mpoly_push_term_fmpz_ffmpz, fmpz_mpoly_push_term_ui_ffmpz, fmpz_mpoly_push_term_si_ffmpz, fq_nmod_mpoly_push_term_fq_nmod_ffmpz, nmod_mpoly_push_term_ui_ffmpz, fmpq_mpoly_push_term_fmpz_ffmpz, fmpq_mpoly_push_term_fmpq_ffmpz, fmpq_mpoly_push_term_ui_ffmpz, fmpq_mpoly_push_term_si_ffmpz, fmpq_mpoly_push_term_fmpq_ffmpz``

List of removals
................

* The following functions that were present in FLINT 2.9, Arb 2.23 or
  Calcium 0.4 have been removed, deprecated, or replaced.
  Most are algorithms obsoleted by new gr implementations,
  functions dealing with removed types (fmpr) or GMP types (mpz, etc.),
  and internal functions that are no longer needed.
* ``__fmpz_clear, __fmpz_eq, __fmpz_gt, __fmpz_gte, __fmpz_init, __fmpz_init_set, __fmpz_init_set_ui, __fmpz_lt, __fmpz_lte, __fmpz_neg, __fmpz_neq, __fmpz_set_si, __fmpz_set_ui``
* ``__fmpz_mod_poly_div_divconquer, __fmpz_mod_poly_divrem_divconquer, __fq_nmod_poly_divrem_divconquer, __fq_poly_divrem_divconquer, __fq_zech_poly_divrem_divconquer``
* ``__nmod_poly_div_divconquer, __nmod_poly_divrem_divconquer, __nmod_poly_invsqrt_series_prealloc``
* ``_acb_poly_compose_axnc, _acb_poly_compose_divconquer, _acb_poly_compose_horner, _acb_poly_compose_series_brent_kung, _acb_poly_compose_series_horner, _acb_poly_sin_cos_series_basecase, _acb_poly_sin_cos_series_tangent, _acb_poly_taylor_shift_convolution, _acb_poly_taylor_shift_divconquer, _acb_poly_taylor_shift_horner``
* ``acb_rising_ui_bs, acb_rising_ui_rs, acb_rising_ui_rec``
* ``_arb_poly_compose_axnc, _arb_poly_compose_divconquer, _arb_poly_compose_horner, _arb_poly_compose_series_brent_kung, _arb_poly_compose_series_horner, _arb_poly_sin_cos_series_basecase, _arb_poly_sin_cos_series_tangent, _arb_poly_taylor_shift_convolution, _arb_poly_taylor_shift_divconquer, _arb_poly_taylor_shift_horner``
* ``arb_rising_ui_bs, arb_rising_ui_rs, arb_rising_ui_rec, arb_rising2_ui_bs, arb_rising2_ui_rs, arb_rising2_ui``
* ``_arith_bernoulli_number_vec_zeta, _arith_bernoulli_number_zeta, _arith_cos_minpoly, _arith_euler_number_zeta, _arith_number_of_partitions_mpfr``
* ``_ca_poly_atan_series, _ca_poly_compose_divconquer, _ca_poly_compose_horner``
* ``_fmpq_poly_set_array_mpq``
* ``_fmpr_add_1x1, _fmpr_add_eps, _fmpr_add_mpn, _fmpr_mul_1x1, _fmpr_mul_mpn, _fmpr_normalise_naive, _fmpr_set_round, _fmpr_set_round_mpn``
* ``_fmpz_deprecated_multi_crt_local_size, _fmpz_deprecated_multi_crt_run, _fmpz_deprecated_multi_crt_run_p, _fmpz_mod_poly_compose_divconquer, _fmpz_mod_poly_compose_divconquer_recursive, _fmpz_mod_poly_compose_horner, _fmpz_mod_poly_div_basecase, _fmpz_mod_poly_div_divconquer, _fmpz_mod_poly_div_divconquer_recursive, _fmpz_mod_poly_div_newton, _fmpz_mod_poly_divrem_divconquer, _fmpz_mod_poly_divrem_divconquer_recursive, _fmpz_mod_poly_gcd_cofactors, _fmpz_mod_poly_gcd_euclidean, _fmpz_mod_poly_gcd_hgcd, _fmpz_mod_poly_hgcd_recursive, _fmpz_mod_poly_hgcd_recursive_iter, _fmpz_mod_poly_hgcd_res, _fmpz_mod_poly_xgcd_euclidean, _fmpz_mod_poly_xgcd_hgcd, _fmpz_poly_evaluate_mpfr``
* ``_fmpz_ui_pow_ui, _fmpz_vec_get_mpf_vec``
* ``_fq_nmod_poly_compose_divconquer, _fq_nmod_poly_compose_horner, _fq_nmod_poly_div_basecase, _fq_nmod_poly_divrem_basecase, _fq_nmod_poly_divrem_divconquer, _fq_nmod_poly_divrem_divconquer_recursive, _fq_nmod_poly_gcd_euclidean, _fq_nmod_poly_gcd_hgcd, _fq_nmod_poly_hgcd, _fq_nmod_poly_hgcd_recursive, _fq_nmod_poly_hgcd_recursive_iter, _fq_nmod_poly_xgcd_euclidean``
* ``_fq_poly_compose_divconquer, _fq_poly_compose_horner, _fq_poly_div_basecase, _fq_poly_divrem_basecase, _fq_poly_divrem_divconquer, _fq_poly_divrem_divconquer_recursive, _fq_poly_gcd_euclidean, _fq_poly_gcd_hgcd, _fq_poly_hgcd, _fq_poly_hgcd_recursive, _fq_poly_hgcd_recursive_iter, _fq_poly_xgcd_euclidean``
* ``_fq_zech_poly_compose_divconquer, _fq_zech_poly_compose_horner, _fq_zech_poly_div_basecase, _fq_zech_poly_divrem_basecase, _fq_zech_poly_divrem_divconquer, _fq_zech_poly_divrem_divconquer_recursive, _fq_zech_poly_gcd_euclidean, _fq_zech_poly_gcd_hgcd, _fq_zech_poly_hgcd, _fq_zech_poly_hgcd_recursive, _fq_zech_poly_hgcd_recursive_iter, _fq_zech_poly_xgcd_euclidean``
* ``_nmod_mat_set_mod``
* ``_nmod_poly_compose_divconquer, _nmod_poly_compose_series_brent_kung, _nmod_poly_compose_series_divconquer, _nmod_poly_compose_series_horner, _nmod_poly_div_basecase, _nmod_poly_div_basecase_1, _nmod_poly_div_basecase_2, _nmod_poly_div_basecase_3, _nmod_poly_div_divconquer, _nmod_poly_div_divconquer_recursive, _nmod_poly_div_newton, _nmod_poly_divrem_basecase_1, _nmod_poly_divrem_basecase_2, _nmod_poly_divrem_basecase_3, _nmod_poly_divrem_divconquer, _nmod_poly_divrem_divconquer_recursive, _nmod_poly_divrem_newton, _nmod_poly_divrem_q0, _nmod_poly_divrem_q1, _nmod_poly_exp_series_basecase, _nmod_poly_exp_series_monomial_ui, _nmod_poly_exp_series_newton, _nmod_poly_hgcd_recursive, _nmod_poly_hgcd_recursive_iter, _nmod_poly_hgcd_res, _nmod_poly_integral_offset, _nmod_poly_log_series_monomial_ui, _nmod_poly_rem_basecase, _nmod_poly_rem_basecase_1, _nmod_poly_rem_basecase_2, _nmod_poly_rem_basecase_3``
* ``acb_poly_compose_divconquer, acb_poly_compose_horner, acb_poly_compose_series_brent_kung, acb_poly_compose_series_horner, acb_poly_sin_cos_series_basecase, acb_poly_sin_cos_series_tangent, acb_poly_taylor_shift_convolution, acb_poly_taylor_shift_divconquer, acb_poly_taylor_shift_horner``
* ``arb_flint_get_num_available_threads``
* ``arb_poly_compose_divconquer, arb_poly_compose_horner, arb_poly_compose_series_brent_kung, arb_poly_compose_series_horner, arb_poly_sin_cos_series_basecase, arb_poly_sin_cos_series_tangent, arb_poly_taylor_shift_convolution, arb_poly_taylor_shift_divconquer, arb_poly_taylor_shift_horner``
* ``arb_test_multiplier``
* ``arb_thread_pool_num_available``
* ``arf_get_fmpr, arf_set_fmpr``
* ``arith_cos_minpoly, arith_number_of_partitions_mpfr``
* ``ca_mat_transpose_resize, ca_poly_atan_series, ca_poly_compose_divconquer, ca_poly_compose_horner, calcium_test_multiplier``
* ``cos_minpoly, cos_pi_pq``
* ``fmpq_poly_evaluate_mpq, fmpq_poly_evaluate_mpz, fmpq_poly_get_coeff_mpq, fmpq_poly_scalar_div_mpq, fmpq_poly_scalar_div_mpz, fmpq_poly_scalar_mul_mpq, fmpq_poly_scalar_mul_mpz, fmpq_poly_set_array_mpq, fmpq_poly_set_coeff_mpq, fmpq_poly_set_coeff_mpz, fmpq_poly_set_mpq, fmpq_poly_set_mpz``
* ``fmpr_add, fmpr_add_fmpz, fmpr_add_naive, fmpr_add_si, fmpr_add_ui, fmpr_addmul, fmpr_addmul_fmpz, fmpr_addmul_si, fmpr_addmul_ui, fmpr_check_ulp, fmpr_cmp, fmpr_cmp_2exp_si, fmpr_cmpabs, fmpr_cmpabs_2exp_si, fmpr_cmpabs_ui, fmpr_div, fmpr_div_fmpz, fmpr_div_si, fmpr_div_ui, fmpr_exp, fmpr_expm1, fmpr_fmpz_div, fmpr_fmpz_div_fmpz, fmpr_get_d, fmpr_get_fmpq, fmpr_get_fmpz, fmpr_get_fmpz_2exp, fmpr_get_fmpz_fixed_fmpz, fmpr_get_fmpz_fixed_si, fmpr_get_mpfr, fmpr_get_si, fmpr_log, fmpr_log1p, fmpr_mul, fmpr_mul_fmpz, fmpr_mul_naive, fmpr_mul_si, fmpr_mul_ui, fmpr_pow_sloppy_fmpz, fmpr_pow_sloppy_si, fmpr_pow_sloppy_ui, fmpr_print, fmpr_printd, fmpr_randtest, fmpr_randtest_not_zero, fmpr_randtest_special, fmpr_root, fmpr_rsqrt, fmpr_set_d, fmpr_set_fmpq, fmpr_set_fmpz_2exp, fmpr_set_mpfr, fmpr_set_round_ui_2exp_fmpz, fmpr_set_round_uiui_2exp_fmpz, fmpr_si_div, fmpr_sqrt, fmpr_sub, fmpr_sub_fmpz, fmpr_sub_si, fmpr_sub_ui, fmpr_submul, fmpr_submul_fmpz, fmpr_submul_si, fmpr_submul_ui, fmpr_ui_div, fmpr_ulp``
* ``fmpz_deprecated_multi_crt, fmpz_deprecated_multi_crt_clear, fmpz_deprecated_multi_crt_init, fmpz_deprecated_multi_crt_precomp, fmpz_deprecated_multi_crt_precomp_p, fmpz_deprecated_multi_crt_precompute, fmpz_deprecated_multi_crt_precompute_p``
* ``fmpz_mat_col_equal, fmpz_mat_get_mpf_mat, fmpz_mat_row_equal``
* ``fmpz_mod_ctx_get_modulus_mpz_read_only``
* ``fmpz_mod_poly_compose_divconquer, fmpz_mod_poly_compose_horner, fmpz_mod_poly_div_basecase, fmpz_mod_poly_div_divconquer, fmpz_mod_poly_div_newton, fmpz_mod_poly_divrem_divconquer, fmpz_mod_poly_gcd_euclidean, fmpz_mod_poly_gcd_hgcd, fmpz_mod_poly_get_coeff_mpz, fmpz_mod_poly_set_coeff_mpz, fmpz_mod_poly_xgcd_euclidean, fmpz_mod_poly_xgcd_hgcd``
* ``fmpz_poly_evaluate_mpfr, fmpz_poly_evaluate_mpq, fmpz_poly_get_coeff_mpz, fmpz_poly_q_evaluate, fmpz_poly_q_scalar_div_mpq, fmpz_poly_q_scalar_div_mpz, fmpz_poly_q_scalar_mul_mpq, fmpz_poly_q_scalar_mul_mpz, fmpz_poly_scalar_divexact_mpz, fmpz_poly_scalar_fdiv_mpz, fmpz_poly_scalar_mul_mpz, fmpz_poly_set_coeff_mpz, fmpz_poly_set_mpz``
* ``fq_nmod_poly_compose_divconquer, fq_nmod_poly_compose_horner, fq_nmod_poly_divrem_basecase, fq_nmod_poly_divrem_divconquer, fq_nmod_poly_gcd_euclidean, fq_nmod_poly_gcd_hgcd, fq_nmod_poly_xgcd_euclidean, fq_poly_compose_divconquer, fq_poly_compose_horner, fq_poly_divrem_basecase, fq_poly_divrem_divconquer, fq_poly_gcd_euclidean, fq_poly_gcd_hgcd, fq_poly_xgcd_euclidean, fq_zech_poly_compose_divconquer, fq_zech_poly_compose_horner, fq_zech_poly_divrem_basecase, fq_zech_poly_divrem_divconquer, fq_zech_poly_gcd_euclidean, fq_zech_poly_gcd_hgcd, fq_zech_poly_xgcd_euclidean``
* ``mag_get_fmpr, mag_set_fmpr``
* ``mpfr_cos_pi_pq, mpfr_zeta_inv_euler_product``
* ``nmod_poly_compose_divconquer, nmod_poly_compose_series_brent_kung, nmod_poly_compose_series_divconquer, nmod_poly_compose_series_horner, nmod_poly_div_basecase, nmod_poly_div_divconquer, nmod_poly_div_newton, nmod_poly_divrem_divconquer, nmod_poly_divrem_newton, nmod_poly_exp_series_basecase, nmod_poly_exp_series_monomial_ui, nmod_poly_factor_get_nmod_poly, nmod_poly_log_series_monomial_ui, nmod_poly_rem_basecase, nmod_poly_set_fmpz_poly, sinh_cosh_divk_precomp``
* ``_nmod_poly_powmod_mpz_binexp, nmod_poly_powmod_mpz_binexp, _nmod_poly_powmod_mpz_binexp_preinv, nmod_poly_powmod_mpz_binexp_preinv, _nmod_poly_powmod_mpz_binexp, nmod_poly_powmod_mpz_binexp, _nmod_poly_powmod_mpz_binexp_preinv, nmod_poly_powmod_mpz_binexp_preinv``

2022-06-24 -- FLINT 2.9.0
-------------------------------------------------------------------------------

* Add fmpz_mod_poly_div function
* Add _flint_get_memory function
* Add Eulerian polynomials
* Support "multivariate" polynomials with zero variables
* Improve Stirling numbers of both kinds
* Speed up numerous fmpz functions for small inputs
* Improve Bell numbers
* Speedups to nmod arithmetic
* Improve nmod_mat LU decomposition
* Fully separate nmod module from nmod_vec
* Speed up Hermite polynomials
* Add n-th derivative for Z[x] and Q[x]
* Improve fq_default module (nmod is now used where optimal)
* Add sqrt functions for numerous polynomial/series modules and finite fields
* Add FFT matrix multiplication
* Improve CI
* Improve LLL for general use
* Add matrix-vector products over Q
* Add can_solve function for fmpq_mat, handling non-square/singular matrices
* Document fmpz_mod_vec module
* Fix and document qadic_sqrt function
* Add parallel programming helpers

2022-04-25 -- FLINT 2.8.5
-------------------------------------------------------------------------------

* Fix a serious bug in LLL

2021-11-17 -- FLINT 2.8.4
-------------------------------------------------------------------------------

* Fix a serious bug in fmpz_mod_poly_xgcd for polynomials of large length
* Fix an assertion failure in fmpz_mat_solve_fflu (only relevant if asserts enabled)
* Fix some bugs on 32 bit machines
* Work around a compiler bug on Apple M1
* Fix bug in nmod_mpoly_factor (returned 0 for some factorisations)
* Fix some documentation build errors and some doc formatting issues

2021-11-03 -- FLINT 2.8.3
-------------------------------------------------------------------------------

* Fix a serious bug in nmod_poly_xgcd_hgcd, nmod_poly_xgcd, fmpz_poly_xgcd_modular, fmpz_poly_xgcd,
  fmpq_poly_xgcd for polynomials of length >= 340.
* Fix some copyright assignments
* Fix some documentation errors

2021-10-15 -- FLINT 2.8.2
-------------------------------------------------------------------------------

* Fix an issue with --disable-dependency-tracking on some distributions

2021-10-01 -- FLINT 2.8.1
-------------------------------------------------------------------------------

* Numerous bug fixes
* Adjust soname on android
* Allow disabling of dependency tracking

2021-07-23 -- FLINT 2.8.0
-------------------------------------------------------------------------------

* New fq_default module which combines existing finite fields
* Speedups for linear algebra when using BLAS and/or threading
* New series expansions with coefficients in QQ
* Faster CRT
* New fmpz_mod_mpoly module
* Polynomial factoring improvements over ZZ
* Fixed bugs in gmpcompat on Windows
* Add fmpz_mat_can_solve_fflu and fmpz_mat_can_solve
* Cleanup of the nmod_poly and nmod_poly_factor code
* Implement nmod_mat_det_howell
* Add fmpz_mod_poly_divides, fmpz_divides, n_divides, nmod_poly_divides
* Interface for multiplying matrices by vectors and arrays
* Nearest Euclidean division
* Subresultant GCD
* XGCD over ZZ with canonical Bezout coefficients
* Add fmpz_mpoly resultant and discriminant
* Add deprecations list
* Add FLINT_SGN macro
* Speedups for series computations
* Switch to GitHub Actions for CI
* Improve Taylor shift
* Numerous bug fixes and speedups

2021-01-18 -- FLINT 2.7.1
-------------------------------------------------------------------------------

* Fix build bug due to missing test files
* Fix bug in fmpz_mod_poly_factor when there are more than five factors
* Fix issue when using MPIR 3.0.0 on Win64 with command line build
* Fix bug in fmpz_mod_poly_div_series
* Fix some broken asserts
* Support standard GNU installation directories in CMake build
* Fix stack overflow with ICC

2020-12-18 -- FLINT 2.7.0
-------------------------------------------------------------------------------

* Multivariate factorisation
* Square root and square testing for finite fields
* Square root and square testing for multivariates
* Zassenhaus factoring speedups (incl. degree pruning)
* Fast factorisation of cubic univariate polynomials
* Add context objects to fmpz_mod_poly functions
* Use BLAS for matrix multiplication over Z/nZ (small n)
* Linear solving for non-square/singular matrices (can_solve)
* Speed up factorisation over Z/nZ (for multiprecision n)

2020-08-12 -- FLINT 2.6.3
-------------------------------------------------------------------------------

* Fix a bug in generator of finite field in characteristic 2
* Allow Flint to work with GMP 6.1.2 and 6.2.0 interchangeably
* Fix some old license headers

2020-07-31 -- FLINT 2.6.2
-------------------------------------------------------------------------------

* Fix for choice of generator in an fq finite field of degree one
* Fix an incorrectly written test

2020-07-23 -- FLINT 2.6.1
-------------------------------------------------------------------------------

* Fix issues on Debian major architectures
* Fix numerous reported bugs (mpoly, fq_poly, mpn_mul_1, mod 2 code, etc.)

2020-06-05 -- FLINT 2.6.0
-------------------------------------------------------------------------------

* multivariate polynomials over most standard rings (sparse distributed)
* APR-CL primality proving
* elliptic curve integer factoring
* minpoly and charpoly
* improved quadratic sieve for integer factoring
* embeddings of finite fields
* pollard rho integer factoring
* p+1 integer factoring
* best of breed smooth integer factoring routine
* best of breed general integer factoring routine
* howell and strong echelon form
* large speedups for solve and hence inverse over Z and Q
* randprime and nextprime functions
* pernet-stein HNF improvements
* moller-granlund precomputed inverses
* resultant_modular_div
* fibonacci polynomials
* exception mechanism/flint_abort
* sqrt of series and polynomials
* division of series over Z
* power sums
* improved base cases of various power series functions
* ability to switch memory allocators
* fast recurrence for Hermite polys
* shifted Legendre polynomials
* Laguerre polynomials
* Gegenbauer polys
* sphinx documentation
* van hoeij with gradual feeding implementation of polynomial factoring over Z
* perfect power detection
* divisibility testing for polynomials
* fast block based memory manager for bundling fmpz allocations
* uniform random generation
* CMake build system
* linear algebra speedups when everything can be kept in longs
* nmod module for integers mod (small) n
* fmpz_mod_mat module for matrices over integers mod multiprecision n
* kronecker product (tensor product)
* random primitive polys (for finite fields)
* thread pool implementation
* threading of FFT for integer and polynomial multiplication over Z
* threading of quadratic sieve for integer factoring
* improved threading of factoring of polynomials mod p
* threading for multivariate polynomial multiplication, division and GCD
* threaded multiplication of matrices mod p
* Berlekamp-Massey (nmod)
* fmpz_mod module for integers mod multiprecision n
* Pohlig-Hellman (discrete log)
* farey_neighbours
* remove openMP option
* additional integer division variants
* speed up mpn_mulmod_preinv
* fft precaching
* cyclotomic polynomial detection
* polynomial root finding over finite fields
* GMP 6.2 support
* MPIR 3.0.0 support
* many small speedups and additional convenience functions added

2015-08-13 -- FLINT 2.5.2
-------------------------------------------------------------------------------

* Fix further issues with soname versioning and ldconfig
* Fix a bug when using GMP instead of MPIR.

2015-08-12 -- FLINT 2.5.1
-------------------------------------------------------------------------------

* Fix some build bugs related to soname versioning and ldconfig
* Fix issue with Windows MSVC build

2015-08-07 -- FLINT 2.5.0
-------------------------------------------------------------------------------

* LLL (rational, Nguyen-Stehle, from Gram matrix, with_removal, Storjohann/ULLL)
* Hermite normal form (naive, xgcd, Domich-Kannan-Trotter, Kannan-Bachem, Pernet-Stein)
* Smith normal form (diagonal, Kannan-Bachem, Iliopoulos)
* Paterson-Stockmeyer algorithm
* modular resultant
* hgcd resultant
* polynomial discriminant
* multithreaded multimodular Taylor shift
* multithreaded Brent-Kung composition
* multithreaded Kaltofen-Shoup distinct degree factorisation
* multiplication based reduced row echelon form
* place inline functions in library for foreign function interfaces
* Primality tests for large integers (Pocklington, Morrison)
* Probable prime tests for large integers (Lucas, Baillie-PSW, strong-psp, Brillhart-Lehmer-Selfridge)
* CRT for large integers
* Dixon algorithm for nullspace
* Brent-Kung composition in irreducibility and distinct degree factorisation
* floating point QR decomposition
* Schwarz-Rutishauser Gram-Schmidt algorithm
* Ogita-Rump-Oishi dot product
* matrix window functions
* MSVC support (Brian Gladman)
* fast cube/nth-root (Newton, Kahan, magic number, Chebyshev)
* Bodrato matrix squaring
* matrix concatenation functions
* matrix content
* faster n_gcd
* faster n_sqrtmod and fmpz_sqrtmod
* additional functions for returning factor of modulus in polys over Z/nZ
* Hadamard matrix construction
* series addition/subtraction
* faster prime_pi bounds
* speedup creation of sparse polynomials
* speedup n_isprime n_nextprime
* speedup n_isprime_pocklington
* speedups to fmpq_poly and fmpz_poly arithmetic
* speedup polynomial irreducibility testing over Z/pZ
* speedup of rank computation over ZZ
* made CPimport compile time dependency only
* teach flint_printf/sprintf about explicit width format specifiers
* support relative paths in configure
* library soname versioning
* ARM64 patches
* Support MSYS2
* Progress towards supporting MIPS64
* Fix a serious bug in fmpz_poly_signature

????-??-?? -- FLINT 2.4.5
-------------------------------------------------------------------------------

* fixed a severe bug in flint's fmpz_poly_gcd_heuristic, reported by
  Anton Mellit.

????-??-?? -- FLINT 2.4.4
-------------------------------------------------------------------------------

* fixed a severe bug in flint's primality code (n_is_prime() affecting n_factor())

2014-04-01 -- FLINT 2.4.3
-------------------------------------------------------------------------------

* Fix a linker issue on Mac OSX.

2014-03-11 -- FLINT 2.4.2
-------------------------------------------------------------------------------

* Fix bug in ARM assembly

2012-11-20 -- FLINT 2.4
-------------------------------------------------------------------------------

* C++ expressions template wrapper
* Fast factorisation of polynomials over Z/nZ
* improved p-adics
* polynomials/matrices over p-adics
* qadics
* Finite fields (small and large F_q), polynomials/matrices over F_q
* Finite fields with Zech logarithm representation
* Fast factorisation of polynomials over F_q
* Faster Brent-Kung modular composition
* New prime sieving code
* Lambert-W function
* Precomputed inverses for polynomials and large integers
* Williams' P+1 integer factoring algorithm
* Harvey's KS2/KS4 polynomial multiplication  
* Faster primality testing up to 64 bits
* Support for Cygwin64 and MinGW64
* Support for Clang
* Support for GMP
* Support for Boehm-Demers-Weiser GC
* Support for flint extension modules

2012-07-01 -- FLINT 2.3
-------------------------------------------------------------------------------

* general

  * many changes to the build system
  * added NTL interface
  * switched to custom memory allocation functions flint_malloc etc
  * in addition to the entries below, fixed a large number of memory leaks,
    problems with the test code, and bugs in corner cases of various functions
  * added _fmpz_cleanup_mpz_content as an alternative to _fmpz_cleanup
  * support MinGW32
  * support Cygwin
  * bugfix on ia64
  * support sparc32/sparc64
  * support OSX
  * support Solaris, NetBSD, OpenBSD, etc (if bash, GNU Make present)

* ulong_extras

  * implemented the improved Lehman algorithm
  * added n_jacobi_unsigned to allow n > WORD_MAX
  * fixed n_sqrtmod for n > WORD_MAX
  * fixed bug causing n_sqrtmod to hang
  * added sublinear algorithm for computing factorials mod p
  * added n_sqrtmod_primepow, n_sqrtmodn and associated functions for
    computing square roots modulo composite integers
  * fixed bugs in n_is_prime_pocklington
  * fixed UWORD_MAX case in powmod and powmod2
  * fixed problems with the random number generators
  * fixed rare bug in n_mod_precomp
  * fixed rare bug in n_is_prime_pseudosquare

* long_extras

  * added z_sizeinbase

* qsieve

  * new module implementing a quadratic sieve for numbers up to two limbs

* fft

  * new module providing an efficient Schoenhage-Strassen FFT

* longlong

  * added assembly code for ia64 and ARM
  * fixed bug in fallback version of add_sssaaaaaa

* fmpz

  * added fmpz_fib_ui
  * added double precision natural logarithm
  * added fmpz_val2 for 2-valuation
  * added mul_2exp, div_2exp, cdiv_q_2exp, tdiv_q_2exp, fdiv_r, fdiv_r_2exp,
    tdiv_ui, mul_tdiv_q_2exp
  * added get_d/set_d
  * added fmpz_divisible, divisible_si
  * optimised fmpz_powm and fmpz_powm_ui
  * added clog, clog_ui, flog, flog_ui for computing logarithms
  * added abs_lbound_ui_2exp, ubound_ui_2exp
  * added fmpz_rfac_ui and fmpz_rfac_uiui for rising factorials
  * added functions to obtain read-only fmpz_t's from mpz_t's
  * added fmpz_init_set, init_set_ui
  * added fmpz_gcdinv
  * added fmpz_is_square
  * added fmpz_tstbit, setbit, clrbit, complement, combit, and, or, xor, popcnt
  * added a sign flag for CRT instead of using separate methods
  * fixed bugs in fmpz_sqrtmod
  * fixed a bug in fmpz_bit_unpack that could cause corruption of the global
    fmpz array when compiled in single mode
  * fixed a bug in fmpz_sub_ui that could cause memory corruption

* fmpz_vec

  * added functions for obtaining the largest absolute value coefficient
  * added functions for computing the sum and product of an integer vector
  * made max_bits much faster
  * added _fmpz_vec_mod_fmpz
  * made randtest produce sparse output

* fmpz_poly

  * added fmpz_poly_sqr, fmpz_poly_sqrlow for squaring a polynomial
  * added fmpz_poly_lcm
  * made multipoint interpolation faster by using the Newton basis
  * added a function for fast division by a linear polynomial
  * added power series composition (classical and Brent-Kung)
  * added power series reversion (classical, Newton, fast Lagrange)
  * added a function for obtaining the largest absolute value coefficient
  * fixed quadratic memory usage and stack overflow when performing
    unbalanced division or pseudo division using the divconquer algorithm
  * fixed a bug in poly_zero_coeffs
  * fixed a bug in xgcd_modular
  * allowing +/-1 in the constant term of power series inversion
  * fixed aliasing bug in divrem
  * added restartable Hensel lifting and associated utility functions
  * fixed rem, which used to only call the basecase algorithm
  * fixed pseudo_divrem, which used to only call the basecase algorithm
  * implemented Schoenhage-Strassen multiplication (mul_SS, mullow_SS)
    and enabled this by default
  * fixed a bug in the heuristic GCD algorithm
  * added functions for Newton basis conversion
  * added functions for fast Taylor shift
  * added fmpz_poly_sqrt implementing a basecase algorithm
  * added scalar mul_2exp, fdiv_2exp, tdiv_2exp
  * made randtest produce sparse output
  * added fmpz_poly_equal_fmpz
  * improved performance by always using basecase multiplication
    when one polynomial is short
  * improved algorithm selection for fmpz_poly_gcd
  * fixed several bugs in gcd_modular
  * improved performance of gcd_modular

* fmpz_poly_factor

  * new module for factorisation of fmpz_polys
  * added a naive implementation of the Zassenhaus algorithm

* fmpz_mod_poly

  * new module for polynomials modulo over Z/nZ for arbitrary-precision n
  * multiplication, powering
  * classical and divconquer division
  * series inversion
  * Euclidean GCD and XGCD
  * invmod
  * radix conversion
  * divconquer composition
  * GCD and division functions that test invertibility of the
    leading coefficient

* fmpz_mat

  * added det_divisor for computing a random divisor of the determinant
  * faster determinant computation using divisor trick
  * faster determinant computation by using multimodular updates
  * fixed n x 0 x m product not zeroing the result
  * various interface improvements
  * faster implementation of Cramer's rule for multiple right hand sides
  * added fmpz_mat_fread and read
  * added multi CRT/mod functions
  * added trace

* fmpz_poly_mat

  * fixed n x 0 x m product not zeroing the result
  * added inverse
  * added rank computation
  * added reduced row echelon form and nullspace computation
  * added more utility functions
  * added squaring and exponentiation
  * added balanced product of a sequence of matrices
  * added truncate, mullow, sqrlow, pow_trunc
  * added trace

* fmpz_factor

  * new module providing interface for integer factorisation
  * fast expansion of a factored integer

* fmpq

  * cleaned up and improved performance of rational reconstruction code
  * allow separate numerator and denominator bounds for rational
    reconstruction
  * added continued fraction expansion
  * added functions for summation using binary splitting
  * added fmpq_swap
  * added fmpq_print, fmpq_get_str
  * added fmpq_pow_si
  * added functions to obtain read-only fmpq_t's from mpq_t's
  * added fmpq_cmp

* fmpq_mat

  * fixed n x 0 x m product not zeroing the result
  * added fmpq_mat_transpose
  * added trace

* fmpq_poly

  * improved speed of multipoint interpolation using _fmpz_poly_div_root
  * fmpq_poly: added power series composition (classical and Brent-Kung)
  * fmpq_poly: added power series reversion (classical, Newton, fast Lagrange)
  * fixed bug wherein set_array_mpq modified the input
  * added gcd, xgcd, lcm, resultant
  * added fmpq_poly_set_fmpq
  * added fmpq_poly_get_slice, fmpq_poly_reverse
  * fixed aliasing bug in divrem
  * changed some functions to use FLINT scalar types instead of MPIR data types
  * added fmpq_poly_get_numerator

* nmod_poly

  * implemented the half gcd algorithm for subquadratic gcd and xgcd
  * added multipoint evaluation and interpolation
  * added asymptotically fast multipoint evaluation and interpolation
  * added a function for forming the product of linear factors
  * added a function for fast division by a linear polynomial
  * added power series composition (classical and Brent-Kung)
  * added power series reversion (classical, Newton, fast Lagrange)
  * added nmod_poly_mulmod, powmod and related functions
    (ported from flint1)
  * added squarefree, irreducibility tests (ported from flint1)
  * added Berlekamp and Cantor-Zassenhaus factoring (ported from flint1)
  * fixed quadratic memory usage and stack overflow when performing
    unbalanced division using the divconquer algorithm
  * added compose_series_divconquer
  * added resultant
  * fixed aliasing bug in divrem
  * added rem functions
  * added divrem_q0, q1 for special cases of division
  * added functions for fast Taylor shift
  * added nmod_poly_sqrt
  * made fread read the modulus from the file
  * made randtest produce sparse output
  * fixed bug in xgcd_euclidean with scalar inputs

* nmod_vec

  * added functions and macros for computing dot products
  * made randtest produce sparse output

* nmod_mat

  * added addmul/submul functions
  * asymptotically fast solving of triangular systems
  * asymptotically fast LUP decomposition
  * asymptotically fast determinant and rank computation
  * asymptotically fast reduced row echelon form and nullspace
  * asymptotically fast nonsingular solving
  * asymptotically fast inverse
  * tidied some interfaces
  * fixed n x 0 x m product not zeroing the result
  * added trace
  * made multiplication faster for tiny moduli by means of bit packing

* nmod_poly_mat

  * new module for matrices over Z/nZ[x], with similar
    functionality as the fmpz_poly_mat module
  * determinant, rank, solving, reduced echelon form, nullspace
  * fraction-free Gaussian elimination
  * multiplication using bit packing
  * multiplication using evaluation-interpolation
  * determinant using evaluation-interpolation

* padic

  * restructured and improved much of the code
  * added padic_log
  * improved log and exp using rectangular splitting
  * added asymptotically fast log and exp based on binary splitting

* perm

  * added the perm module for permutation matrices
  * computing the parity of a permutation
  * inverting a permutation

* arith

  * added generation of cyclotomic polynomials
  * added functions for evaluating Dedekind sums
  * fast computation of the partition function
  * added a function for factoring a Hardy-Ramanujan-Rademacher
    type exponential sum
  * added Chebyshev polynomials T and U
  * added computation of the minimal polynomial of cos(2pi/n)
  * added asymptotically fast high-precision approximation of zeta(n)
  * added asymptotically fast computation of Euler's constant
  * added new algorithms and functions for computing Bell numbers
  * fast computation of pi (adapting code written by Hanhong Xue)
  * added functions for computing the number of sum of squares
    representations of an integer
  * renamed functions to have an arith prefix

2011-06-04 -- FLINT 2.2
-------------------------------------------------------------------------------

* fmpq (multiprecision rational numbers)

  * Basic arithmetic functions
  * Utility functions
  * Rational reconstruction
  * Functions for enumerating the rationals

* fmpq_mat (matrices over Q)

  * Basic arithmetic functions
  * Utility functions
  * Fast multiplication
  * Classical and fraction-free reduced row echelon form
  * Determinants
  * Fast non-singular solving

* fmpz_poly_mat (matrices over Z[x]

  * Basic arithmetic functions
  * Utility functions
  * Fraction-free row reduction and determinants
  * Fast determinants (experimental)

* fmpz_mat

  * Added more utility functions (scalar multiplication, etc)
  * Added Dixon's p-adic algorithm (used by fast nonsingular rational system
    solving)
  * Added reduced row echelon form
  * Added conversions between fmpz_mat and nmod_mat
  * Added CRT functions for fmpz_mats
  * Faster matrix multiplication for small to medium dimensions

* longlong.h

  * Added x86 assembly macros for accumulating sums of two limb operands

* nmod_mat

  * Sped up arithmetic for moduli close to FLINT_BITS

* arith

  * Changed interface of various functions to use new fmpq type

* fmpz

  * Added fmpz_set_ui_mod
  * Inlined fmpz_neg, fmpz_set_si, fmpz_set_ui for better performance
  * Added fmpz_lcm
  * Small performance improvement to fmpz_CRT_ui

* fmpz_vec

  * Added _fmpz_vec_lcm

* fmpz_poly_q (rational functions over Q, modeled as quotients of fmpz_polys)

  * Basic arithmetic functions
  * Conversion and IO functions
  * Evaluation

* padic (p-adic numbers -- experimental)

  * Basic arithmetic
  * Data conversion and IO
  * Inverse and square root using Newton iteration
  * Teichmuller lifts (not optimised)
  * p-adic exponential function (not optimised)

* fmpz_poly

  * Added fmpz_poly_gcd_modular (and fmpz_poly_gcd wrapper)
  * Added fmpz_poly_xgcd_modular (and fmpz_poly_xgcd wrapper)
  * Added conversions between fmpz_poly and nmod_poly
  * Added CRT functions
  * Added multipoint evaluation and interpolation

* nmod_poly

  * Added nmod_poly_xgcd_euclidean (and nmod_poly_xgcd wrapper)
  * nmod_poly_gcd wrapper

* mpn_extras

  * Added MPN_NORM and MPN_SWAP macros.
  * Added mpn_gcd_full to remove some of the restrictions from the usual mpn_gcd

* build fixes

  * fixed make install to create nonexistent dirs (reported by Serge Torres)
  * -L use /usr instead of /usr/local by default (reported by Serge Torres)
  * guards for system headers because of flint's use of ulong

2011-03-09 -- FLINT 2.1
-------------------------------------------------------------------------------

* fmpz

  * Simplified interface for fast multimodular reduction and CRT reconstruction
  * Fixed segmentation fault in fmpz_multi_mod_ui when the input exceeds the product of the moduli
  * Added simple incremental CRT functions (fmpz_CRT_ui, fmpz_CRT_ui_unsigned) to complement the existing fast ones
  * Added example programs for CRT
  * Added random number generators designed for testing modular code (fmpz_randtest_mod, fmpz_randtest_mod_signed)
  * Added fmpz_fdiv_ui for remainder on division by an ulong
  * Added fmpz_bin_uiui for computing binomial coefficients
  * Added fmpz_mul2_uiui and fmpz_divexact2_uiui for multiplying or dividing an fmpz by a pair of ulongs (efficiently when their product fits in a single limb)

* fmpz_mat

  * Added utility functions for basic arithmetic and creating unit matrices
  * Added multimodular determinant computation (certified or heuristic)
  * Added support for computing right nullspaces (fmpz_mat_kernel). Fast only for small matrices.
  * Some internal code cleanup and various small fixes

* nmod_mat

  * Faster Gaussian elimination for small moduli
  * Faster determinants
  * Faster matrix inversion and nonsingular solving

* nmod_poly

  * Added nmod_poly_integral for computing integrals
  * Added fast square root and inverse square root of power series
  * Added fast transcendental functions of power series (log, exp, sin, cos, tan, sinh, cosh, tanh, asin, atan, asinh, atanh)
  * Made nmod_poly_inv_series_newton more memory efficient

* fmpq_poly

  * Added fmpq_poly_integral for computing integrals
  * Added fast transcendental functions of power series (log, exp, sin, cos, tan, sinh, cosh, tanh, asin, atan, asinh, atanh)

* arith

  * Made computation of vectors of Bernoulli numbers faster
  * Added fast computation of single Bernoulli numbers
  * Added a separate function for computing denominators of Bernoulli numbers
  * Added fast computation of Bell numbers (vector and single)
  * Added fast computation of Euler numbers (vector and single)
  * Added fast computation of Euler polynomials
  * Added fast computation of Swinnerton-Dyer polynomials
  * Added fast computation of Legendre polynomials
  * Added fast vector computation of the partition function
  * Added fast vector computation of Landau's function

* ulong_extras

  * Added a function for computing factorials mod n

* build system

  * Added support for building static and shared libraries
  * All object files and test/profile/example binaries now build in separate build directory

* documentation

  * Large number of corrections

2011-01-16 -- FLINT 2.0
-------------------------------------------------------------------------------

N.B: FLINT 2 is  a complete rewrite of flint from scratch
It includes the following modules:

* ulong_extras: (word sized integers and modular arithmetic)

  * random numbers (randint, randbits, randprime, randint)
  * powering
  * reverse binary
  * mod, divrem, mulmod all with precomputed inverses
  * gcd, invgcd, xgcd
  * jacobi symbols
  * addmod, submod, invmod, powmod
  * prime sieve, nextprime, prime-pi, nth-prime
  * primality testing (small, binary search, Pocklington-Lehmer, Pseudosquare)
  * probably prime tests (strong base-a, Fermat, Fibonacci, BPSW, Lucas)
  * sqrt, sqrtrem, is-square, perfect-power (2,3,5)
  * remove, is-squarefree
  * factorisation (trial-range, trial, power (2,3,5), one-line, SQUFOF)
  * Moebius mu, Euler phi
   
* fmpz: (memory managed multiple precision integers)

  * memory management (init, clear)
  * random numbers (randbits, randm)
  * conversion to and from long, ulong, doubles, mpz's, strings
  * read/write to file, stdin, stdout
  * sizeinbase, bits, size, sgn, swap, set, zero
  * cmp, cmp-ui, cmpabs, equal, is-zero, is-one
  * neg, abs, add, add-ui, sub, sub-ui, mul, mul-si, mul-ui, mul-2exp
  * addmul, addmul-ui, submul, submul-ui
  * cdiv-q, cdiv-q-si, cdiv-q-ui
  * fdiv-q, fdiv-q-si, fdiv-q-ui, fdiv-qr, fdiv-q-2exp
  * tdiv-q, tdiv-q-si
  * divexact, divexact-si, divexact-ui
  * mod, mod-ui
  * powering
  * sqrt, sqrt-rem
  * factorial
  * gcd, invmod
  * bit-pack, bit-unpack
  * multimodular reduction, CRT 

* fmpz_vec: (vectors over fmpz's)

  * memory management (init, clear)
  * random vectors
  * max-bits, max-limbs
  * read/write to file/stdin/stdout
  * set, swap, zero, neg
  * equal, is-zero
  * sort
  * add, sub
  * scalar multiplication by fmpz, ulong, long, 2exp
  * exact division by fmpz, long, ulong
  * fdiv-q by fmpz, long, ulong, 2exp
  * tdiv-q by fmpq, long, ulong
  * addmul by fmpz, long, long by 2exp
  * submul by fmpz, long, long by 2exp
  * Gaussian content

* fmpz_poly: (polys over fmpz's)

  * memory management (init, realloc, fit-length, clear)
  * random polys
  * normalise, set-length, truncate
  * length, degree, max-limbs, max-bits
  * set, set-si, set-ui, set-fmpz, set-str
  * get-str, get-str-pretty
  * zero, one, zero-coeffs
  * swap, reverse
  * get/set coeffs from fmpz, long, ulong
  * get-coeff-ptr, lead
  * equal, is-zero
  * add, sub
  * scalar multiplication by fmpz, long, ulong
  * scalar addmul/submul by fmpz
  * scalar fdiv by fmpz, long, ulong
  * scalar tdiv by fmpz, long, ulong
  * scalar divexact by fmpz, long, ulong
  * bit pack, bit unpack
  * multiplication (classical, karatsuba, KS)
  * mullow (classical, karatsuba, KS)
  * mulhigh (classical, karatsuba)
  * middle product (classical)
  * powering (small, binary exponent, binomial, multinomial, addition chains)
  * truncated powering (binary exponent)
  * shift left/right
  * euclidean norm
  * gcd (subresultant)
  * resultant
  * content, primitive part
  * divrem (basecase, divide-and-conquer)
  * div (basecase, divide-and-conquer)
  * rem (basecase)
  * invert series (basecase, Newton)
  * div series
  * pseudo divrem (basecase, divide-and-conquer, Cohen)
  * rem (Cohen)
  * div
  * evaluate (Horner) at fmpz, mpq, a mod n
  * composition (Horner, divide-and-conquer)
  * signature
  * read/write to file/stdin/stdout

* fmpq_poly: (polynomials over Q stored as poly over fmpz with fmpz denominator)

  * memory management (init, realloc, fit-length, clear)
  * random polys
  * set-length, canonicalise, normalise, truncate
  * is-canonical, length, degree
  * reference to numerator, denominator
  * set, set-si, set-ui, set-fmpz, set-mpz, set-mpq
  * set-array-mpq, set-str
  * get-str, get-str-pretty
  * zero, neg, swap
  * invert
  * set coefficient to mpq, long, ulong, fmpz, mpz
  * get coefficient as mpq
  * equal, cmp, is-one, is-zero
  * add, sub
  * scalar multiplication by long, ulong, fmpz, mpq
  * scalar division by fmpz, long, ulong, mpq
  * multiplication, mullow
  * powering 
  * shift left/right
  * divrem, div, rem
  * invert series (Newton iteration)
  * divide series
  * derivative
  * evaluate at fmpz, mpq
  * composition, scale by constant
  * content, primitive part
  * make-monic, is-monic
  * is-squarefree
  * read/write to file/stdin/stdout

* nmod_vec: (vectors over Z/nZ for n fitting in a machine word)

  * memory management (init/clear)
  * macros for efficient reduction of 1, 2 and 3 limb integers mod n
  * macro for addmul mod n
  * add/sub/neg individual coefficients mod n
  * random vectors
  * set, zero, swap
  * reduce, max-bits
  * equal
  * add, sub, neg
  * scalar multiplication by a value reduced mod n
  * scalar addmul by a value reduced mod n

* nmod_poly: (polynomials over Z/nZ for n fitting in a machine word)

  * memory management (init, realloc, fit-length, clear)
  * random polys
  * normalise, truncate
  * length, degree, modulus, max-bits
  * set, swap, zero, reverse
  * get/set coefficients as ulongs, strings
  * read/write to file, stdin, stdout
  * equal, is-zero
  * shift left/right
  * add, sub, neg
  * scalar multiplication by a value reduced mod n
  * make-monic
  * bit pack, bit unpack
  * multiplication (classical, KS)
  * mullow (classical, KS)
  * mulhigh (classical)
  * powering (binary exponent)
  * pow-trunc (binary exponent)
  * divrem (basecase, divide-and-conquer, Newton iteration)
  * div (basecase, divide-and-conquer, Newton iteration)
  * invert series (basecase, Newton iteration)
  * divide series (Newton iteration)
  * derivative
  * evaluation at a value taken mod n
  * composition (Horner, divide-and-conquer)
  * gcd (euclidean)

* fmpz_mat: (matrices over fmpz's)

  * memory management (init, clear)
  * random matrices (bits, integer relations, simultaneous diophantine equations
    NTRU-like, ajtai, permutation of rows and cols of a diagonal matrix, random
    of given rank, random of given determinant, random elementary operations)
  * set, init-set, swap, entry pointer
  * write to file or stdout
  * equal
  * transpose
  * multiplication (classical, multimodular)
  * inverse
  * determinant
  * row reduce (Gaussian and Gauss-Jordan fraction-free elimination)
  * rank
  * solve Ax = b, solve AX = B
  * fraction free LU decomposition

* nmod_mat: (matrices over Z/nZ for n fitting in a machine word)

  * memory management (init, clear)
  * random matrices (uniform, full, permutations of diagonal matrix, random of 
    given rank, random elementary operations)
  * set, equal
  * print to stdout
  * add
  * transpose
  * multiplication (classical, Strassen, A*B^T)
  * row reduction (Gaussian and Gauss-Jordan)
  * determinant
  * rank
  * solve (Ax = b, AX = B, solve with precomputed LU)
  * invert

* arith: (arithmetic functions)

  * Bernoulli numbers
  * Bernoulli polynomials
  * primorials (product of primes up to n)
  * harmonic numbers
  * Stirling numbers
  * Euler phi function
  * Moebius mu function
  * Sigma (sum of powers of divisors)
  * Ramanujan tau function

* examples: (example programs)

  * compute coefficients of q-series of Delta function

* mpfr_vec: (vectors over mpfr reals)

  * memory management (init, clear)
  * add
  * set, zero
  * scalar multiplication by mpfr, 2exp
  * scalar product

* mpfr_mat: (matrices over mpfr reals)

  * memory management (init, clear)

2010-12-24 -- FLINT 1.6.0
-------------------------------------------------------------------------------

* Bugs:

  * Fixed a memory leak in mpz_poly_to_string_pretty
  * Fixed a bug inherited from an old version of fpLLL 
  * Makefile to respect CC and CXX
  * Fixed bug in F_mpz_set_si
  * Fixed bug in F_mpz_equal
  * Most for loops to C90 standard (for easier MSVC porting)
  * Better Cygwin support
  * Fixed a bug in zmod_poly_resultant
  * Fixed bug in F_mpz_mul_KS/2
  * Fixed bug in tinyQS
  * Worked around some known bugs in older GMP/MPIR's

* Major new functionality

  * F_mpz_poly_factor_zassenhaus 
  * F_mpz_poly_factor (incl. fmpz_poly_factor wrapper) using new vH-N approach
    (see the paper of van Hoeij and Novocin and the paper of van Hoeij, Novocin 
    and Hart)
  * Implementation of new CLD bounds function for polynomial factors
    (see the paper of van Hoeij, Novocin and Hart
  * Restartable Hensel lifting
  * Heuristic LLL implementations using doubles and mpfr
  * LLL implementations optimised for knapsack lattices
  * New (probably subquadratic) LLL implementation (ULLL)
  * zmod_poly_factor_cantor_zassenhaus
  * New F_mpz_mod_poly module for polynomials over Z/pZ for multiprec. p

* Some of the other new functions added

  * F_mpz

   * F_mpz_gcd
   * F_mpz_smod
   * F_mpz_mod_preinv
   * F_mpz_fdiv_qr
   * F_mpz_get/set_mpfr/2exp
   * F_mpz_sscanf
   * F_mpz_set_d

  * F_mpz_poly:

   * read F_mpz_poly to_string/from_string/fprint/print/fread/pretty
   * F_mpz_poly_to/from_zmod_poly
   * F_mpz_poly_scalar_div_exact
   * F_mpz_poly_smod
   * F_mpz_poly_derivative, F_mpz_poly_content, F_mpz_poly_eval_horner_d/2exp
   * F_mpz_poly_scalar_abs
   * F_mpz_poly_set_d_2exp
   * F_mpz_poly_div/divrem
   * F_mpz_poly_gcd
   * F_mpz_poly_is_squarefree
   * F_mpz_poly_factor_squarefree
   * F_mpz_poly_mul_trunc_left
   * F_mpz_poly_pseudo_div
   * F_mpz_poly_set_coeff
   * F_mpz_poly_pow_ui
   * Inflation/deflation trick for factorisation

 * zmod_poly:

   * Inflation/deflation trick for factorisation

  * mpz_mat:

   * mpz_mat_from_string/to_string/fprint/fread/pretty

  * mpq_mat:

   * mpq_mat_init/clear
   * Gramm-schmidt Orthogonalisation

  * F_mpz_mat:

   * F_mpz_mat_print/fprint/fread/pretty
   * F_mpz_mat_mul_classical
   * F_mpz_mat_max_bits/2
   * F_mpz_mat_scalar_mul/div_2exp
   * F_mpz_mat_col_equal
   * F_mpz_mat_smod
   * F_mpz_vec_scalar_product/norm
   * F_mpz_vec_add/submul_ui/si/F_mpz/2exp

  * zmod_mat:

   * classical multiplication
   * strassen multiplication
   * scalar multiplication
   * zmod_mat_equal
   * zmod_mat_add/sub
   * zmod_mat_addmul_classical

  * d_mat:

   * d_vec_norm, d_vec_scalar_product

  * mpfr_mat:

   * mpfr_vec_scalar_product/norm

2009-09-22 -- FLINT 1.5.0
-------------------------------------------------------------------------------

* Added multimodular reduction and CRT to F_mpz module
* Fixed some bugs in F_mpz module and numerous bugs in test code
* Added zmod_poly_compose
* Added zmod_poly_evaluate
* Added functions for reduced evaluation and composition to fmpz_poly module (contributed by Burcin Erocal)
* Fixed bugs in the primality tests in long_extras
* Removed all polynomial multimodular multiplication functions
* Added new thetaproduct code used in the 1 trillion triangles computation
* Fixed a severe bug in the fmpz_poly_pseudo_div function reported by Sebastian Pancratz
* Added fmpz_comb_temp_init/clear functions
* Fixed a normalisation buglet in fmpz_poly_pack_bytes
* Added F_mpz_pow_ui function (contributed by Andy Novocin)
* Fixed a severe bug in the FFT reported by William Stein and Mariah Lennox (fix contributed by David Harvey)
* Removed some memory leaks from F_mpz test code
* Fixed bug in zmod_poly_evaluate test code

2009-07-06 -- FLINT 1.4.0
-------------------------------------------------------------------------------

* Sped up zmod_poly division in case where operands are the same length
* Sped up zmod_poly division in case where operands have lengths differing by 1
* Fixed a bug in zmod_poly_gcd for polynomials of zero length
* Sped up zmod_poly_gcd considerably (both euclidean and half gcd)
* Sped up zmod_poly_gcd_invert and zmod_poly_xgcd considerably
* Made zmod_poly_gcd_invert and zmod_poly_xgcd asymptotically fast
* Made zmod_poly_resultant asymptotically fast
* Added optimised zmod_poly_rem function
* Fixed a divide by zero bug in zmod_poly_factor_berlekamp 
* Added test code for z_factor_tinyQS and z_factor_HOLF
* Fixed many bugs in the z_factor code, tinyQS and mpQS
* Corrected numerous typos in the documentation and added missing descriptions
* Added F_mpz_cmp function
* Added documentation to the manual for the new F_mpz module

2009-06-09 -- FLINT 1.3.0
-------------------------------------------------------------------------------

* Added new code for checking 2nd, 3rd and 5th roots
* Fixed a bug in z_factor
* Connected quadratic sieve for factoring large ulongs
* Added one line factor algorithm
* constructed best of breed factor algorithm
* Fixed termination conditions for z_intcuberoot and z_intfifthroot which were broken
* Added some code for special cases which cause infinite loops in cuberoot and fifthroot
* Went back to ceil(pow(n, 0.33333333)) and ceil(pow(n, 0.2)) for initial guesses in cube and fifth root functions as these were about 50% faster than sqrt(n) and sqrt(sqrt(n)) respectively.
* Added test code for z_intfifthroot
* Added test code for z_factor_235power
* Fixed some uninitialised data found by valgrind in intcuberoot and intfifthroot
* Fixed multiply defined PRIME_COUNT in long_extras-test
* Got rid of gotos in some functions in long_extras
* Knocked optimisation level back to -O2 because it miscompiles on sage.math
* Changed tables to use uint64_t's instead of ulongs which are not 64 bits on a 32 bit machine
* Only checked MAX_HOLF on 64 bit machine
* Changed MAX_SQUFOF to WORD(-1)
* Check constant 0x3FFFFFFFFUL only on a 64 bit machine
* Fixed a bug in z_oddprime_lt_4096 on 32 bit machines
* Fixed some TLS issues with Cygwin
* Fixed some typos in makefile
* Fixed a wrong path in fmpz.c

2009-04-18 -- FLINT 1.2.5
-------------------------------------------------------------------------------

* Upgraded to zn_poly-0.9 to avoid a serious error in squaring of large polynomials over Z/nZ

2009-04-04 -- FLINT 1.2.4
-------------------------------------------------------------------------------

* Defined THREAD to be blank on Apple CC and __thread for thread local storage on other gcc's (where it's defined)
* #undef ulong in profiler.h where time.h and other system time headers are included (both reported by M. Abshoff)

2009-03-31 -- FLINT 1.2.3
-------------------------------------------------------------------------------

* Fixed bugs in all fmpz_poly evaluation functions, identified by Burcin Erocal. 

2009-03-20 -- FLINT 1.2.2
-------------------------------------------------------------------------------

* Fixed a memory leak in zmod_poly_factor
* Fixed zmod_poly-profile build

2009-03-14 -- FLINT 1.2.1
-------------------------------------------------------------------------------

* Removed some FLINT 2.0 code which was interfering with the build of the NTL-interface
* Removed an omp.h from fmpz_poly.c.

2009-03-10 -- FLINT 1.2.0
-------------------------------------------------------------------------------

* made memory manager, fmpz and fmpz_poly threadsafe
* Code for running tests in parallel (not activated)
* Sped up fmpz_poly_scalar_div_ui/si when scalar is 1/-1
* Parallelise _fmpz_poly_mul_modular 
* fmpz_mul_modular_packed to pack coefficients to the byte before running _fmpz_poly_mul_modular 
* fmpz_poly_pseudo_rem_cohen (not documented)
* special case for leading coeff 1/-1 in fmpz_poly_pseudo_divrem_basecase
* removed a memory allocation bug which caused a massive slowdown in fmpz_poly_pseudo_divrem_basecase
* fmpz_poly_pseudo_rem_basecase (not documented)
* fmpz_poly_pseudo_rem (not asymptotically fast)
* fmpz_poly_signature (not asymptotically fast)
* basic fmpz_poly_is_squarefree function
* included zn_poly in trunk and made FLINT build zn_poly as part of its build process
* switched to using zn_poly for polynomial multiplication, newton inversion, scalar multiplication in zmod_poly
* Integer cube root of word sized integers
* Fibonacci pseudoprime test
* BSPW probable prime test
* n - 1 primality test
* Complete implementation of z_issquarefree
* Significantly improved the thetaproduct example program. 
* Fixed bug in fmpz_poly_byte_pack which is triggered when trying to pack into fields a multiple of 8 bytes (could cause a segfault)
* Fixed a bug in fmpz_poly_pseudo_divrem (relied on an uninitialised poly to have length 0)
* Fixed bug in fmpz_multi_CRT_ui (could segfault)
* Fixed bug in fmpz_multi_mod_ui (could segfault)
* Fixed memory leak in zmod_poly_factor_squarefree 
* Fixed memory leak in zmod_poly_from_string

2009-03-01 -- FLINT 1.1.3
-------------------------------------------------------------------------------

* Inserted some missing return values in zmod_poly test code.

2009-03-01 -- FLINT 1.1.2
-------------------------------------------------------------------------------

* Fixed some memory allocation slowdowns and bugs in fmpz_poly division and pseudo division functions (reported by William Stein).
  
2009-02-11 -- FLINT 1.1.1
-------------------------------------------------------------------------------

* Fixed bugs in _fmpz_poly_scalar_mul_fmpz, fmpz_poly_gcd_heuristic and fmpz_poly_gcd_subresultant and fixed bugs in test__fmpz_poly_scalar_div_fmpz, test_fmpz_poly_scalar_div_fmpz and test_fmpz_poly_scalar_div_mpz.

2008-12-21 -- FLINT 1.1.0
-------------------------------------------------------------------------------

Some of the following features were previewed in FLINT 1.0.11.

* integer gcd (this just wraps the GMP gcd code)
* polynomial content
* primitive part
* convert to and from FLINT and NTL integers and polynomials
* get a coefficient of a polynomial efficiently as a read only mpz_t
* print polynomials in a prettified format with a specified variable 
* Sped up integer multiplication
* Convert to and from zmod_polys from fmpz_polys
* Chinese remainder for fmpz_polys
* Leading coeff macro
* Euclidean norm of polynomials
* Exact division testing of polynomials
* Polynomial GCD (subresultant, heuristic, modular)
* Modular inversion of polynomials
* Resultant
* XGCD (Pohst-Zassenhaus)
* Multimodular polynomial multiplication
* Rewritten karatsuba_trunc function
* Rewritten division functions
* Polynomial derivative
* Polynomial evaluation
* Polynomial composition
* Addition and subtraction of zmod_polys
* Sped up multiplication of zmod_polys
* Extended multiplication of zmod_polys to allow up to 63 bit moduli
* zmod_poly subpolynomials
* zmod_poly reverse
* Truncated multiplication for zmod_polys (left, right, classical and KS)
* Scalar multiplication
* Division for zmod_polys (divrem and div, classical, divide and conquer and newton)
* Series inversion for zmod_polys
* Series division for zmod_polys
* Resultant for zmod_polys
* GCD for zmod_polys including half-gcd
* Inversion modulo a polynomial for zmod_polys
* XGCD for zmod_polys
* Squarefree factorisation for zmod_polys
* Berlekamp factorisation for zmod_polys
* Irreducibility testing for zmod_polys
* Derivative for zmod_polys
* Middle product for zmod_polys (sped up newton division)
* addmod, submod and divmod for ulongs
* Sped up limb sized integer square root
* Partial factoring of ulongs
* z_randbits
* Pocklington-Lehmer primality testing
* BSPW pseudo-primality testing
* Fermat pseudo-primality testing
* Fast Legendre symbol computation
* Chinese remainder for fmpzs
* Square root with remainder for fmpzs
* Left and right shift for fmpzs
* Reduction modulo a ulong for fmpzs
* Montgomery redc, mulmod, divmod and mod for fmpzs
* Multimodular reduction and CRT for fmpzs
* fmpz_mulmod and fmpz_divmod
* fmpz_invert for inversion modulo an fmpz
* Dramatically sped up gcd for small fmpzs
* Computation of 1D, 2D and some 3D theta functions
* Example program for multiplying theta functions
* Test code now times test functions
* Quick and dirty timing function for profiler
* Tiny quadratic sieve for small one and two limb integers
* Completely rewritten self initialising multiple polynomial quadratic sieve
* Build fix for 64 bit OSX dylibs (reported by Michael Abshoff)

2008-12-25 -- FLINT 1.0.21
-------------------------------------------------------------------------------

* Fixed the Christmas bug reported by Michael Abshoff which causes a test failure in fmpz_poly_gcd_modular and a hang in fmpz_poly_invmod_modular on 32 bit machines 

2008-12-13 -- FLINT 1.0.20
-------------------------------------------------------------------------------

* Fixed some bugs in conversion of zmod_poly's to and from strings

2008-12-12 -- FLINT 1.0.19
-------------------------------------------------------------------------------

* Fixed a bug in z_remove_precomp

2008-12-05 -- FLINT 1.0.18
-------------------------------------------------------------------------------

* Fixed another bug in the fmpz_poly_set_coeff_* functions which resulted in dirty coefficients

2008-11-30 -- FLINT 1.0.17
-------------------------------------------------------------------------------

* Fixed a segfault caused by left shifting of polynomials with zero limbs allocated in division and pseudo division functions.
* Fixed a bound used in fmpz_gcd_modular to use a proven bound
* Fixed a bug in fmpz_poly-profile where the top bit of random coefficients of n bits was always set

2008-10-22 -- FLINT 1.0.16
-------------------------------------------------------------------------------

* Fixed a segfault when trying to truncate a polynomial to an longer length than it currently is, with the function fmpz_poly_truncate (reported by Craig Citro)

2008-10-15 -- FLINT 1.0.15
-------------------------------------------------------------------------------

* Fixed a bug which causes a segfault when setting a coefficient of the zero polynomial to zero
* Fixed build bug in longlong.h on s390 platform

2008-09-23 -- FLINT 1.0.14
-------------------------------------------------------------------------------

* Update long_extras and test code for the sake of new quadratic sieve (new functions in long_extras remain undocumented)
* Removed many bugs from tinyQS and mpQS and joined them into a single program for factoring integers

2008-07-13 -- FLINT 1.0.13
-------------------------------------------------------------------------------

* Fixed memory leaks and dirty memory issues in test code for numerous modules.
* Removed further issues with cache prefetching in mpn_extras.c

2008-07-11 -- FLINT 1.0.12
-------------------------------------------------------------------------------

* Removed some Opteron tuning flags which cause illegal instruction errors on Pentium4
* Fixed numerous memory leaks in fmpz_poly test code
* Fixed memory leak in fmpz_poly_power_trunc_n
* Fixed major memory leaks in fmpz_poly_xgcd_modular
* Rewrote __fmpz_poly_mul_karatrunc_recursive and _fmpz_poly_mul_karatsuba_trunc to "prove code" and got rid of some dirty memory issues
* Fixed some potential illegal memory accesses to do with cache prefetching in fmpz_poly.c

2008-07-09 -- FLINT 1.0.11
-------------------------------------------------------------------------------

* Fixed a bug in z_ll_mod_precomp on ia64 (reported by Michael Abshoff and William Stein)

2008-06-16 -- FLINT 1.0.10
-------------------------------------------------------------------------------

* integer gcd (this just wraps the GMP gcd code)
* polynomial content
* convert to and from FLINT and NTL integers and polynomials
* get a coefficient of a polynomial efficiently as a read only mpz_t
* print polynomials in a prettified format with a specified variable         

2008-03-11 -- FLINT 1.0.9
-------------------------------------------------------------------------------

* Fixed a memory allocation bug in fmpz_poly_power

2008-02-15 -- FLINT 1.0.8
-------------------------------------------------------------------------------

* Fixed a bug in fmpz_poly_right_shift (reported by Kiran Kedlaya)

2008-01-22 -- FLINT 1.0.7
-------------------------------------------------------------------------------

* Made F_mpn_mul binary compatible with the way mpn_mul *operates* in practice.

2008-01-17 -- FLINT 1.0.6
-------------------------------------------------------------------------------

* Fixed an issue with FLINT_BIT_COUNT on certain machines (probably due to arithmetic shift issues)

2008-01-05 -- FLINT 1.0.5
-------------------------------------------------------------------------------

* Fixed some inline issues which cause problems because of the C99 inline rules (reported by David Harvey). 
* Fixed a makefile issue reported (and solved) by David Harvey when *not* linking against NTL.

2008-01-04 -- FLINT 1.0.4
-------------------------------------------------------------------------------

* Fixed a bug in the bernoulli_zmod example program and associated polynomial zmod code which caused memory corruption.
* Fixed a bug in the fmpz-test code which manifested on 32 bit machines, reported by David Harvey.   
* Fixed some bugs in the pari profiling code.     

2007-12-16 -- FLINT 1.0.3
-------------------------------------------------------------------------------

* Fixed a bug in the polynomial memory management code which caused a segfault
* Fixed a bug in the pseudo division code which caused a block overrun

2007-12-10 -- FLINT 1.0.2
-------------------------------------------------------------------------------

* Rewrote tuning code for integer multiplication functions, making it more robust and fixing a bug
  which showed up on 32 bit machines (reported by Michael Abshoff and Jaap Spies). Factored the tuning
  code out into a number of macros.

2007-12-07 -- FLINT 1.0.1
-------------------------------------------------------------------------------

* Fixed a bug in _fmpz_poly_maxbits1 on 32 bit machines, reported by Michael Abshoff and Carl Witty
* Removed some instances of u_int64_t and replaced them with uint64_t, reported by Michael Abshoff
* Replaced sys/types.h with stdint.h
* Added FLINT macros to documentation
* Corrected numerous typos in documentation   

2007-12-02 -- FLINT 1.0
-------------------------------------------------------------------------------

* First version of FLINT, includes fmpz_poly, fmpz and mpQS


Antic version history
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

2021-06-24 -- Antic 0.2.5
-------------------------------------------------------------------------------

* TODO: list changes here

2021-04-15 -- Antic 0.2.4
-------------------------------------------------------------------------------

* TODO: list changes here

2020-12-11 -- Antic 0.2.3
-------------------------------------------------------------------------------

* TODO: list changes here

2020-06-30 -- Antic 0.2.2
-------------------------------------------------------------------------------

* TODO: list changes here

2020-06-16 -- Antic 0.2.1
-------------------------------------------------------------------------------

* TODO: list changes here

2019-02-12 -- Antic 0.2
-------------------------------------------------------------------------------

* Many bug fixes, standalone build system, continuous integration.

2013-05-12 -- Antic 0.1
-------------------------------------------------------------------------------

* First version of antic, including a qfb module for (positive definite) binary quadratic forms.


Calcium version history
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

2021-05-28 -- Calcium 0.4
-------------------------------------------------------------------------------

* Algebraic numbers

  * Fixed bug in special-casing of roots of unity in qqbar_root_ui.
  * Fixed qqbar_randtest with bits == 1.
  * Faster qqbar_cmp_re for nearby reals.
  * Faster qqbar polynomial evaluation and powering using linear algebra.
  * Improved qqbar_abs, qqbar_abs2 to produce cleaner enclosures.
  * Use a slightly better method to detect real numbers in qqbar_sgn_im.
  * Added qqbar_hash.
  * Added qqbar_get_fmpq, qqbar_get_fmpz.
  * Added qqbar_pow_fmpq, qqbar_pow_fmpz, qqbar_pow_si.
  * Added qqbar_numerator, qqbar_denominator.

* Basic arithmetic and elementary functions

  * Improved ca_condense_field: automatically demote to a simple number field
    when the only used extension number is algebraic.
  * Improved multivariate field arithmetic to automatically remove algebraic
    or redundant monomial factors from denominators.
  * Added ca_pow_si_arithmetic (guaranteed in-field exponentiation).
  * Added polynomial evaluation functions (ca_fmpz_poly_evaluate,
    ca_fmpq_poly_evaluate, ca_fmpz_poly_evaluate, ca_fmpz_mpoly_q_evaluate).
  * Added several helper functions (ca_is_special, ca_is_qq_elem,
    ca_is_qq_elem_zero, ca_is_qq_elem_one, ca_is_qq_elem_integer,
    ca_is_nf_elem, ca_is_cyclotomic_nf_elem, ca_is_generic_elem).
  * Added ca_rewrite_complex_normal_form.
  * Perform direct complex conjugation in cyclotomic fields.
  * Use ca_get_acb_raw instead of ca_get_acb when printing to avoid expensive
    recomputations.
  * Added alternative algorithms for various basic functions.
  * Deep complex conjugation.
  * Use complex conjugation in is_real, is_imaginary, is_negative_real.
  * Added functions for unsafe inversion for internal use.
  * Significantly stronger zero testing in mixed algebraic-transcendental fields.
  * Added ca_arg.
  * Added special case for testing equality between number field elements
    and rationals.
  * Added ca_sin_cos, ca_sin, ca_cos, ca_tan and variants.
  * Added ca_atan, ca_asin, ca_acos and variants.
  * Added ca_csgn.
  * Improved ca_get_acb and ca_get_acb_accurate_parts to fall back on exact
    zero tests when direct numerical evaluation does not give precise enclosures.
  * Added ca_get_decimal_str.
  * More automatic simplifications of logarithms (simplify logarithms of
    exponentials, square roots and powers raised to integer powers).
  * More automatic simplifications of square roots (simplify square roots of
    exponentials, square roots and powers raised to integer powers).
  * Improved order comparisons (ca_check_ge etc.) to handle special values
    and to fall back on strong equality tests.
  * Fixed a test failure in the ca_mat module.

* Polynomials

  * Added ca_poly_inv_series, ca_poly_div_series (power series division).
  * Added ca_poly_exp_series (power series exponential).
  * Added ca_poly_log_series (power series logarithm).
  * Added ca_poly_atan_series (power series arctangent).

* Other

  * Added fmpz_mpoly_q_used_vars.
  * Remove useless rpath line from configure (reported by Julien Puydt).
  * Added missing declaration of fexpr_hash.
  * Fixed crashes on OS X in Python interface (contributed by deinst).
  * Fixed memory leaks in Python string conversions (contributed by deinst).
  * Reserve I, E for symbolic expressions in Python interface.


2021-04-23 -- Calcium 0.3
-------------------------------------------------------------------------------

* Symbolic expressions

  * Added the fexpr module for flat-packed unevaluated symbolic expressions.
  * LaTeX output.
  * Basic manipulation (construction, replacement, accessing subexpressions).
  * Numerical evaluation with Arb.
  * Expanded normal form.
  * Conversion methods for other types.
  * Enable LaTeX rendering of objects in Jupyter notebooks.

* Algebraic numbers

  * Fix a major performance issue (slow root refinement) that made Calcium as a whole far slower than necessary.
  * Added qqbar_cmp_root_order; sort polynomial roots consistently by default.
  * Added qqbar_get_quadratic.
  * Added qqbar_equal_fmpq_poly_val and use it to speed up checking guessed values.
  * Conversion of qqbar_t to and from symbolic expression (qqbar_set_fexpr, qqbar_get_fexpr_repr, qqbar_get_fexpr_root_nearest, qqbar_get_fexpr_root_indexed, qqbar_get_fexpr_formula).
  * Fixed bugs in qqbar_cmpabs_re, cmpabs_im.
  * Optimized qqbar_cmp_im and qqbar_cmpabs_im for conjugates with mirror symmetry.
  * Added qqbar_pow (taking a qqbar exponent).
  * Special-case roots of unity in qqbar_pow_ui, qqbar_root_ui, qqbar_abs and qqbar_abs2.
  * Wrapped qqbar in Python.

* Polynomials

  * Added several utility functions.
  * Optimized polynomial multiplication with rational entries.
  * Fast polynomial multiplication over number fields.

* Matrices

  * Fast matrix multiplication over number fields.
  * Right kernel (ca_mat_right_kernel).
  * Matrix diagonalization (ca_mat_diagonalization).
  * Jordan normal form  (ca_mat_jordan_form, ca_mat_jordan_transformation, ca_mat_jordan_blocks).
  * Matrix exponential (ca_mat_exp).
  * Matrix logarithm (ca_mat_log).
  * Polynomial evaluation (ca_mat_ca_poly_evaluate).
  * Cofactor expansion algorithm for determinant and adjugate (ca_mat_adjugate_cofactor).
  * Added several utility functions.
  * Improved algorithm selection in ca_mat_inv.
  * Solving using the adjugate matrix.
  * Danilevsky characteristic polynomial algorithm (ca_mat_charpoly_danilevsky).

* Field elements

  * Use factoring in ca_sqrt to enable more simplifications.
  * Simplify square roots and logarithms of negative real numbers.
  * Optimized ca_sub.
  * Conversion of ca_t to and from symbolic expressions (ca_set_fexpr, ca_get_fexpr).
  * Added function for assigning elements between context objects (ca_transfer).
  * Fixed a possible memory corruption bug when Vieta's formulas are used.
  * Optimized constructing square roots of rational numbers.

* Other

  * Added demonstration notebook to documentation.
  * Fixed OSX compatibility in Python wrapper (contributed by deinst).
  * Fixed bug in calcium_write_acb.
  * Fixed bug in fmpz_mpoly_vec_set_primitive_unique (contributed by gbunting).


2020-10-16 -- Calcium 0.2
-------------------------------------------------------------------------------

* Basic arithmetic and expression simplification

  * Use GrÃ¶bner basis for reduction ideals,  making simplification much more robust.
  * Compute all linear relations with LLL simultaneously instead of piecemeal.
  * Make monomial ordering configurable (default is lex as before).
  * Use Vieta's formulas to simplify expressions involving conjugate algebraic numbers.
  * Denest exponentials of symbolic logarithms.
  * Denest logarithms of symbolic powers and square roots.
  * Denest powers of symbolic powers.
  * Simplify exponentials that evaluate to roots of unity.
  * Simplify logarithms of roots of unity.
  * Improve ideal reduction to avoid some unnecessary GCD computations.

* Python wrapper

  * Calcium now includes a minimal ctypes-based Python wrapper for testing.

* New ca_mat module for matrices

  * Mostly using naive basecase algorithms.
  * Matrix arithmetic, basic manipulation.
  * Construction of special matrices (Hilbert, Pascal, Stirling, DFT).
  * LU factorization.
  * Fraction-free LU decomposition.
  * Nonsingular solving and inverse.
  * Reduced row echelon form.
  * Rank.
  * Trace and determinant.
  * Characteristic polynomial.
  * Computation of eigenvalues with multiplicities.

* New ca_poly module for polynomials

  * Mostly using naive basecase algorithms.
  * Polynomial arithmetic, basic manipulation.
  * Polynomial division.
  * Evaluation and composition.
  * Derivative and integral.
  * GCD (Euclidean algorithm).
  * Squarefree factorization.
  * Computation of roots with multiplicities.
  * Construction from given roots.

* New ca_vec module for vectors.

  * Memory management and basic scalar operations.

* Bug fixes

  * Fix bug in powering number field elements.
  * Fix bug in qqbar_log_pi_i.
  * Fix aliasing bug in ca_pow.

* New basic functions

  * Conversion from double: ca_set_d, ca_set_d_d.
  * Special functions: ca_erf, ca_erfi, ca_erfc, with algebraic relations.
  * Special functions: ca_gamma (incomplete simplification algorithms).

* New utils_flint module for Flint utilities

  * Vectors of multivariate polynomials.
  * Construction of elementary symmetric polynomials.
  * GrÃ¶bner basis computation (naive Buchberger algorithm).

* Documentation and presentation

  * Various improvements to the documentation.
  * DFT example program.


2020-09-08 -- Calcium 0.1
-------------------------------------------------------------------------------

* Initial test release.
* ca module (exact real and complex numbers).
* fmpz_mpoly_q module (multivariate rational functions over Q).
* qqbar module (algebraic numbers represented by minimal polynomials).
* Example programs.


Arb version history
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

2022-06-29 -- Arb 2.23.0
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


2022-01-25 -- Arb 2.22.1
-------------------------------------------------------------------------------

* Fixed bugs causing some hypergeometric functions hang or crash for
  some input on various non-x86 architectures.
* Fixed a minor bug in acb_hypgeom_m (NaN result sometimes only set the
  real part to NaN).

2022-01-15 -- Arb 2.22.0
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

2021-10-20 -- Arb 2.21.1
-------------------------------------------------------------------------------

* Fixed 32-bit test failures for arb_hypgeom_gamma_fmpq.
* Added pow function to the fpwrap module.
* Added missing header file includes.
* Do not encode the library version in the SONAME on Android (contributed by Andreas Enge).

2021-09-25 -- Arb 2.21.0
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

2021-07-25 -- Arb 2.20.0
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

2020-12-06 -- Arb 2.19.0
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


2020-06-25 -- Arb 2.18.1
-------------------------------------------------------------------------------

* Support MinGW64.
* Added version numbers (__ARB_VERSION, __ARB_RELEASE, ARB_VERSION) to arb.h.

2020-06-09 -- Arb 2.18.0
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

2019-10-16 -- Arb 2.17.0
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


2018-12-07 -- Arb 2.16.0
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


2018-10-25 -- Arb 2.15.1
-------------------------------------------------------------------------------

* Fixed precision issue leading to spurious NaN results in incomplete elliptic integrals

2018-09-18 -- Arb 2.15.0
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

2018-07-22 -- Arb 2.14.0
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


2018-02-23 -- Arb 2.13.0
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

2017-11-29 -- Arb 2.12.0
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

2017-07-10 -- Arb 2.11.1
-------------------------------------------------------------------------------

* Avoid use of a function that was unavailable in the latest public FLINT release

2017-07-09 -- Arb 2.11.0
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

2017-02-27 -- Arb 2.10.0
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

2016-12-02 -- Arb 2.9.0
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

2015-12-31 -- Arb 2.8.1
-------------------------------------------------------------------------------

* Fixed 32-bit test failure for the Laguerre function.
* Made the Laguerre function indeterminate at negative integer orders, to be consistent with the test code.

2015-12-29 -- Arb 2.8.0
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

2015-07-14 -- Arb 2.7.0
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

2015-04-19 -- Arb 2.6.0
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

2015-01-28 -- Arb 2.5.0
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

2014-11-15 -- Arb 2.4.0
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

2014-09-25 -- Arb 2.3.0
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

2014-08-01 -- Arb 2.2.0
-------------------------------------------------------------------------------

* Added functions for computing polylogarithms and order expansions
  of polylogarithms, with support for real and complex s, z.
* Added a missing cast affecting C++ compatibility.
* Generalized powsum functions to allow a geometric factor.
* Improved powsum functions slightly when the exponent is an integer.
* Faster arb_log_ui_from_prev.
* Added mag_sqrt and mag_rsqrt functions.
* Fixed various minor bugs and added missing tests and documentation entries.

2014-06-20 -- Arb 2.1.0
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

2014-05-27 -- Arb 2.0.0
-------------------------------------------------------------------------------

* New modules mag, arf, arb, arb_poly, arb_mat, acb, acb_poly,
  acb_mat for higher-performance ball arithmetic.

* Poly_roots2 and hilbert_matrix2 example programs.

* Vector dot product and norm functions (contributed by Abhinav Baid).

2014-05-03 -- Arb 1.1.0
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

2013-12-21 -- Arb 1.0.0
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

2013-08-07 -- Arb 0.7
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

2013-05-31 -- Arb 0.6
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

2013-03-28 -- Arb 0.5
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

2013-01-26 -- Arb 0.4
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

2012-11-07 -- Arb 0.3
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

2012-09-29 -- Arb 0.2
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

2012-09-14 -- Arb 0.1
-------------------------------------------------------------------------------

* 2012-08-05 - Began simplified rewrite.
* 2012-04-05 - Experimental ball and polynomial code (first commit).

