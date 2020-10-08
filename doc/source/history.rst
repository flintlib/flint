.. _history:

History and changes
===============================================================================

For more details, view the commit log
in the git repository https://github.com/fredrik-johansson/calcium

Old releases of the code can be accessed from
https://github.com/fredrik-johansson/calcium/releases

Future - version 0.2-git
-------------------------------------------------------------------------------

* Simplification and basic arithmetic

  * Compute Gröbner bases for reduction ideals,  making simplification much more robust.
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

  * TODO: matrix functions.

* New ca_poly module for polynomials

  * TODO: polynomial functions.

* New ca_vec module for vectors.

  * TODO: vector functions

* Bug fixes

  * Fix bug in powering number field elements.
  * Fix bug in qqbar_log_pi_i.
  * Fix aliasing bug in ca_pow.

* New basic functions

  * Conversion from double: ca_set_d, ca_set_d_d.
  * Special functions: ca_erf, ca_erfi, ca_erfc, with algebraic relations.
  * Special functions: ca_gamma (incomplete simplification algorithms).

* New utils_flint module for Flint utilities

  * TODO: utility methods.

* Documentation and presentation

  * Various improvements to the documentation.
  * DFT example program.


2012-09-08 - version 0.1
-------------------------------------------------------------------------------

* Initial test release.
* ca module (exact real and complex numbers).
* fmpz_mpoly_q module (multivariate rational functions over Q).
* qqbar module (algebraic numbers represented by minimal polynomials).
* Example programs.


.. raw:: latex

    \newpage

