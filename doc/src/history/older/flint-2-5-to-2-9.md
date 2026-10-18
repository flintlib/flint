---
title: FLINT 2.5 to 2.9
summary: Releases from 2015 to 2022: multivariate polynomials and their factorisation, APR-CL primality proving, ECM and the quadratic sieve, the fq_default finite fields, the last FLINT 2 release 2.9.0.
order: 1
---

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
