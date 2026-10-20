---
title: FLINT 2
summary: Releases from 2011 to 2022: the rewrite of FLINT 1 as FLINT 2 with fmpz, fmpz_poly, nmod_poly, matrices, p-adic numbers and finite fields; later multivariate polynomials and their factorisation, APR-CL primality proving, ECM and the quadratic sieve, and the fq_default finite fields.
order: 1
---

FLINT 2.0, released in January 2011, was a complete rewrite of FLINT 1.
The FLINT 2 series ended with version 2.9.0 in June 2022, the last
release before FLINT 3.0. FLINT 2.6.0 (June 2020) added multivariate
polynomials over most standard rings, APR-CL primality proving and
elliptic curve factoring.

## FLINT 2.9.0 (2022-06-24) {#flint-2-9-0}

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

## FLINT 2.8.5 (2022-04-25) {#flint-2-8-5}

* Fix a serious bug in LLL

## FLINT 2.8.4 (2021-11-17) {#flint-2-8-4}

* Fix a serious bug in fmpz_mod_poly_xgcd for polynomials of large length
* Fix an assertion failure in fmpz_mat_solve_fflu (only relevant if asserts enabled)
* Fix some bugs on 32 bit machines
* Work around a compiler bug on Apple M1
* Fix bug in nmod_mpoly_factor (returned 0 for some factorisations)
* Fix some documentation build errors and some doc formatting issues

## FLINT 2.8.3 (2021-11-03) {#flint-2-8-3}

* Fix a serious bug in nmod_poly_xgcd_hgcd, nmod_poly_xgcd, fmpz_poly_xgcd_modular, fmpz_poly_xgcd,
  fmpq_poly_xgcd for polynomials of length >= 340.
* Fix some copyright assignments
* Fix some documentation errors

## FLINT 2.8.2 (2021-10-15) {#flint-2-8-2}

* Fix an issue with --disable-dependency-tracking on some distributions

## FLINT 2.8.1 (2021-10-01) {#flint-2-8-1}

* Numerous bug fixes
* Adjust soname on android
* Allow disabling of dependency tracking

## FLINT 2.8.0 (2021-07-23) {#flint-2-8-0}

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

## FLINT 2.7.1 (2021-01-18) {#flint-2-7-1}

* Fix build bug due to missing test files
* Fix bug in fmpz_mod_poly_factor when there are more than five factors
* Fix issue when using MPIR 3.0.0 on Win64 with command line build
* Fix bug in fmpz_mod_poly_div_series
* Fix some broken asserts
* Support standard GNU installation directories in CMake build
* Fix stack overflow with ICC

## FLINT 2.7.0 (2020-12-18) {#flint-2-7-0}

* Multivariate factorisation
* Square root and square testing for finite fields
* Square root and square testing for multivariates
* Zassenhaus factoring speedups (incl. degree pruning)
* Fast factorisation of cubic univariate polynomials
* Add context objects to fmpz_mod_poly functions
* Use BLAS for matrix multiplication over Z/nZ (small n)
* Linear solving for non-square/singular matrices (can_solve)
* Speed up factorisation over Z/nZ (for multiprecision n)

## FLINT 2.6.3 (2020-08-12) {#flint-2-6-3}

* Fix a bug in generator of finite field in characteristic 2
* Allow FLINT to work with GMP 6.1.2 and 6.2.0 interchangeably
* Fix some old license headers

## FLINT 2.6.2 (2020-07-31) {#flint-2-6-2}

* Fix for choice of generator in an fq finite field of degree one
* Fix an incorrectly written test

## FLINT 2.6.1 (2020-07-23) {#flint-2-6-1}

* Fix issues on Debian major architectures
* Fix numerous reported bugs (mpoly, fq_poly, mpn_mul_1, mod 2 code, etc.)

## FLINT 2.6.0 (2020-06-05) {#flint-2-6-0}

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

## FLINT 2.5.2 (2015-08-13) {#flint-2-5-2}

* Fix further issues with soname versioning and ldconfig
* Fix a bug when using GMP instead of MPIR.

## FLINT 2.5.1 (2015-08-12) {#flint-2-5-1}

* Fix some build bugs related to soname versioning and ldconfig
* Fix issue with Windows MSVC build

## FLINT 2.5.0 (2015-08-07) {#flint-2-5-0}

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

## FLINT 2.4.5 {#flint-2-4-5}

* fixed a severe bug in flint's fmpz_poly_gcd_heuristic, reported by
  Anton Mellit.

## FLINT 2.4.4 {#flint-2-4-4}

* fixed a severe bug in flint's primality code (n_is_prime() affecting n_factor())

## FLINT 2.4.3 (2014-04-01) {#flint-2-4-3}

* Fix a linker issue on Mac OSX.

## FLINT 2.4.2 (2014-03-11) {#flint-2-4-2}

* Fix bug in ARM assembly

## FLINT 2.4 (2012-11-20) {#flint-2-4}

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

## FLINT 2.3 (2012-07-01) {#flint-2-3}

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

## FLINT 2.2 (2011-06-04) {#flint-2-2}

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

## FLINT 2.1 (2011-03-09) {#flint-2-1}

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

## FLINT 2.0 (2011-01-16) {#flint-2-0}

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
