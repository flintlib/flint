---
title: FLINT 1
summary: Releases from 2007 to 2010: FLINT 1.0 to 1.6, the original library for polynomial and integer arithmetic.
order: 2
---

FLINT 1, the original version of the library, was released from
version 1.0 (December 2007), which contained `fmpz_poly`, `fmpz` and the
quadratic sieve `mpQS`, to 1.6.0 (December 2010). It was replaced by
FLINT 2, a complete rewrite; see [](history/older/flint-2).

## FLINT 1.6.0 (2010-12-24) {#flint-1-6-0}

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

## FLINT 1.5.0 (2009-09-22) {#flint-1-5-0}

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

## FLINT 1.4.0 (2009-07-06) {#flint-1-4-0}

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

## FLINT 1.3.0 (2009-06-09) {#flint-1-3-0}

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

## FLINT 1.2.5 (2009-04-18) {#flint-1-2-5}

* Upgraded to zn_poly-0.9 to avoid a serious error in squaring of large polynomials over Z/nZ

## FLINT 1.2.4 (2009-04-04) {#flint-1-2-4}

* Defined THREAD to be blank on Apple CC and __thread for thread local storage on other gcc's (where it's defined)
* #undef ulong in profiler.h where time.h and other system time headers are included (both reported by M. Abshoff)

## FLINT 1.2.3 (2009-03-31) {#flint-1-2-3}

* Fixed bugs in all fmpz_poly evaluation functions, identified by Burcin Erocal.

## FLINT 1.2.2 (2009-03-20) {#flint-1-2-2}

* Fixed a memory leak in zmod_poly_factor
* Fixed zmod_poly-profile build

## FLINT 1.2.1 (2009-03-14) {#flint-1-2-1}

* Removed some FLINT 2.0 code which was interfering with the build of the NTL-interface
* Removed an omp.h from fmpz_poly.c.

## FLINT 1.2.0 (2009-03-10) {#flint-1-2-0}

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

## FLINT 1.1.3 (2009-03-01) {#flint-1-1-3}

* Inserted some missing return values in zmod_poly test code.

## FLINT 1.1.2 (2009-03-01) {#flint-1-1-2}

* Fixed some memory allocation slowdowns and bugs in fmpz_poly division and pseudo division functions (reported by William Stein).

## FLINT 1.1.1 (2009-02-11) {#flint-1-1-1}

* Fixed bugs in _fmpz_poly_scalar_mul_fmpz, fmpz_poly_gcd_heuristic and fmpz_poly_gcd_subresultant and fixed bugs in test__fmpz_poly_scalar_div_fmpz, test_fmpz_poly_scalar_div_fmpz and test_fmpz_poly_scalar_div_mpz.

## FLINT 1.1.0 (2008-12-21) {#flint-1-1-0}

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

## FLINT 1.0.21 (2008-12-25) {#flint-1-0-21}

* Fixed the Christmas bug reported by Michael Abshoff which causes a test failure in fmpz_poly_gcd_modular and a hang in fmpz_poly_invmod_modular on 32 bit machines

## FLINT 1.0.20 (2008-12-13) {#flint-1-0-20}

* Fixed some bugs in conversion of zmod_poly's to and from strings

## FLINT 1.0.19 (2008-12-12) {#flint-1-0-19}

* Fixed a bug in z_remove_precomp

## FLINT 1.0.18 (2008-12-05) {#flint-1-0-18}

* Fixed another bug in the fmpz_poly_set_coeff_* functions which resulted in dirty coefficients

## FLINT 1.0.17 (2008-11-30) {#flint-1-0-17}

* Fixed a segfault caused by left shifting of polynomials with zero limbs allocated in division and pseudo division functions.
* Fixed a bound used in fmpz_gcd_modular to use a proven bound
* Fixed a bug in fmpz_poly-profile where the top bit of random coefficients of n bits was always set

## FLINT 1.0.16 (2008-10-22) {#flint-1-0-16}

* Fixed a segfault when trying to truncate a polynomial to an longer length than it currently is, with the function fmpz_poly_truncate (reported by Craig Citro)

## FLINT 1.0.15 (2008-10-15) {#flint-1-0-15}

* Fixed a bug which causes a segfault when setting a coefficient of the zero polynomial to zero
* Fixed build bug in longlong.h on s390 platform

## FLINT 1.0.14 (2008-09-23) {#flint-1-0-14}

* Update long_extras and test code for the sake of new quadratic sieve (new functions in long_extras remain undocumented)
* Removed many bugs from tinyQS and mpQS and joined them into a single program for factoring integers

## FLINT 1.0.13 (2008-07-13) {#flint-1-0-13}

* Fixed memory leaks and dirty memory issues in test code for numerous modules.
* Removed further issues with cache prefetching in mpn_extras.c

## FLINT 1.0.12 (2008-07-11) {#flint-1-0-12}

* Removed some Opteron tuning flags which cause illegal instruction errors on Pentium4
* Fixed numerous memory leaks in fmpz_poly test code
* Fixed memory leak in fmpz_poly_power_trunc_n
* Fixed major memory leaks in fmpz_poly_xgcd_modular
* Rewrote __fmpz_poly_mul_karatrunc_recursive and _fmpz_poly_mul_karatsuba_trunc to "prove code" and got rid of some dirty memory issues
* Fixed some potential illegal memory accesses to do with cache prefetching in fmpz_poly.c

## FLINT 1.0.11 (2008-07-09) {#flint-1-0-11}

* Fixed a bug in z_ll_mod_precomp on ia64 (reported by Michael Abshoff and William Stein)

## FLINT 1.0.10 (2008-06-16) {#flint-1-0-10}

* integer gcd (this just wraps the GMP gcd code)
* polynomial content
* convert to and from FLINT and NTL integers and polynomials
* get a coefficient of a polynomial efficiently as a read only mpz_t
* print polynomials in a prettified format with a specified variable

## FLINT 1.0.9 (2008-03-11) {#flint-1-0-9}

* Fixed a memory allocation bug in fmpz_poly_power

## FLINT 1.0.8 (2008-02-15) {#flint-1-0-8}

* Fixed a bug in fmpz_poly_right_shift (reported by Kiran Kedlaya)

## FLINT 1.0.7 (2008-01-22) {#flint-1-0-7}

* Made F_mpn_mul binary compatible with the way mpn_mul *operates* in practice.

## FLINT 1.0.6 (2008-01-17) {#flint-1-0-6}

* Fixed an issue with FLINT_BIT_COUNT on certain machines (probably due to arithmetic shift issues)

## FLINT 1.0.5 (2008-01-05) {#flint-1-0-5}

* Fixed some inline issues which cause problems because of the C99 inline rules (reported by David Harvey).
* Fixed a makefile issue reported (and solved) by David Harvey when *not* linking against NTL.

## FLINT 1.0.4 (2008-01-04) {#flint-1-0-4}

* Fixed a bug in the bernoulli_zmod example program and associated polynomial zmod code which caused memory corruption.
* Fixed a bug in the fmpz-test code which manifested on 32 bit machines, reported by David Harvey.
* Fixed some bugs in the pari profiling code.

## FLINT 1.0.3 (2007-12-16) {#flint-1-0-3}

* Fixed a bug in the polynomial memory management code which caused a segfault
* Fixed a bug in the pseudo division code which caused a block overrun

## FLINT 1.0.2 (2007-12-10) {#flint-1-0-2}

* Rewrote tuning code for integer multiplication functions, making it more robust and fixing a bug
  which showed up on 32 bit machines (reported by Michael Abshoff and Jaap Spies). Factored the tuning
  code out into a number of macros.

## FLINT 1.0.1 (2007-12-07) {#flint-1-0-1}

* Fixed a bug in _fmpz_poly_maxbits1 on 32 bit machines, reported by Michael Abshoff and Carl Witty
* Removed some instances of u_int64_t and replaced them with uint64_t, reported by Michael Abshoff
* Replaced sys/types.h with stdint.h
* Added FLINT macros to documentation
* Corrected numerous typos in documentation

## FLINT 1.0 (2007-12-02) {#flint-1-0}

* First version of FLINT, includes fmpz_poly, fmpz and mpQS
