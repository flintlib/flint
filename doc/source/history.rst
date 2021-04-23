.. _history:

History and changes
===============================================================================

For more details, view the commit log
in the git repository https://github.com/fredrik-johansson/calcium

Old releases of the code can be accessed from
https://github.com/fredrik-johansson/calcium/releases

2021-04-23 - version 0.3
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

