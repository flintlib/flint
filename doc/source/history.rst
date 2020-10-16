.. _history:

History and changes
===============================================================================

For more details, view the commit log
in the git repository https://github.com/fredrik-johansson/calcium

Old releases of the code can be accessed from
https://github.com/fredrik-johansson/calcium/releases

2020-10-16 - version 0.2
-------------------------------------------------------------------------------

* Basic arithmetic and expression simplification

  * Use Gröbner basis for reduction ideals,  making simplification much more robust.
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
  * Gröbner basis computation (naive Buchberger algorithm).

* Documentation and presentation

  * Various improvements to the documentation.
  * DFT example program.


2020-09-08 - version 0.1
-------------------------------------------------------------------------------

* Initial test release.
* ca module (exact real and complex numbers).
* fmpz_mpoly_q module (multivariate rational functions over Q).
* qqbar module (algebraic numbers represented by minimal polynomials).
* Example programs.


.. raw:: latex

    \newpage

