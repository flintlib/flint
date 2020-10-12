.. _ca-mat:

**ca_mat.h** -- matrices over the real and complex numbers
===============================================================================

A :type:`ca_mat_t` represents a dense matrix over the real or
complex numbers,
implemented as an array of entries of type :type:`ca_struct`.
The dimension (number of rows and columns) of a matrix is fixed at
initialization, and the user must ensure that inputs and outputs to
an operation have compatible dimensions. The number of rows or columns
in a matrix can be zero.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: ca_mat_struct

.. type:: ca_mat_t

    Contains a pointer to a flat array of the entries (*entries*), an array of
    pointers to the start of each row (*rows*), and the number of rows (*r*)
    and columns (*c*).

    A *ca_mat_t* is defined as an array of length one of type
    *ca_mat_struct*, permitting a *ca_mat_t* to
    be passed by reference.

.. macro:: ca_mat_entry(mat, i, j)

    Macro giving a pointer to the entry at row *i* and column *j*.

.. macro:: ca_mat_nrows(mat)

    Returns the number of rows of the matrix.

.. macro:: ca_mat_ncols(mat)

    Returns the number of columns of the matrix.

.. function:: ca_ptr ca_mat_entry_ptr(ca_mat_t mat, slong i, slong j)

    Returns a pointer to the entry at row *i* and column *j*.
    Equivalent to :macro:`ca_mat_entry` but implemented as a function.

Memory management
-------------------------------------------------------------------------------

.. function:: void ca_mat_init(ca_mat_t mat, slong r, slong c, ca_ctx_t ctx)

    Initializes the matrix, setting it to the zero matrix with *r* rows
    and *c* columns.

.. function:: void ca_mat_clear(ca_mat_t mat, ca_ctx_t ctx)

    Clears the matrix, deallocating all entries.

.. function:: void ca_mat_swap(ca_mat_t mat1, ca_mat_t mat2, ca_ctx_t ctx)

    Efficiently swaps *mat1* and *mat2*.

.. function:: void ca_mat_window_init(ca_mat_t window, const ca_mat_t mat, slong r1, slong c1, slong r2, slong c2, ca_ctx_t ctx)

    Initializes *window* to a window matrix into the submatrix of *mat*
    starting at the corner at row *r1* and column *c1* (inclusive) and ending
    at row *r2* and column *c2* (exclusive).

.. function:: void ca_mat_window_clear(ca_mat_t window, ca_ctx_t ctx)

    Frees the window matrix.

Conversions
-------------------------------------------------------------------------------

.. function:: void ca_mat_set(ca_mat_t dest, const ca_mat_t src, ca_ctx_t ctx)
              void ca_mat_set_fmpz_mat(ca_mat_t dest, const fmpz_mat_t src, ca_ctx_t ctx)
              void ca_mat_set_fmpq_mat(ca_mat_t dest, const fmpq_mat_t src, ca_ctx_t ctx)

    Sets *dest* to *src*. The operands must have identical dimensions.

Random generation
-------------------------------------------------------------------------------

.. function:: void ca_mat_randtest(ca_mat_t mat, flint_rand_t state, slong depth, slong bits, ca_ctx_t ctx)

    Sets *mat* to a random matrix with entries having complexity up to
    *depth* and *bits* (see :func:`ca_randtest`).

.. function:: void ca_mat_randtest_rational(ca_mat_t mat, flint_rand_t state, slong bits, ca_ctx_t ctx)

    Sets *mat* to a random rational matrix with entries up to *bits* bits in size.


Input and output
-------------------------------------------------------------------------------

.. function:: void ca_mat_print(const ca_mat_t mat, ca_ctx_t ctx)

    Prints *mat* to standard output. The entries are printed on separate lines.

.. function:: void ca_mat_printn(const ca_mat_t mat, slong digits, ca_ctx_t ctx)

    Prints a decimal representation of *mat* with precision specified by *digits*.
    The entries are comma-separated with square brackets and comma separation
    for the rows.

Special matrices
-------------------------------------------------------------------------------

.. function:: void ca_mat_zero(ca_mat_t mat, ca_ctx_t ctx)

    Sets all entries in *mat* to zero.

.. function:: void ca_mat_one(ca_mat_t mat, ca_ctx_t ctx)

    Sets the entries on the main diagonal of *mat* to one, and
    all other entries to zero.

.. function:: void ca_mat_ones(ca_mat_t mat, ca_ctx_t ctx)

    Sets all entries in *mat* to one.

.. function:: void ca_mat_pascal(ca_mat_t mat, int triangular, ca_ctx_t ctx)

    Sets *mat* to a Pascal matrix, whose entries are binomial coefficients.
    If *triangular* is 0, constructs a full symmetric matrix
    with the rows of Pascal's triangle as successive antidiagonals.
    If *triangular* is 1, constructs the upper triangular matrix with
    the rows of Pascal's triangle as columns, and if *triangular* is -1,
    constructs the lower triangular matrix with the rows of Pascal's
    triangle as rows.

.. function:: void ca_mat_stirling(ca_mat_t mat, int kind, ca_ctx_t ctx)

    Sets *mat* to a Stirling matrix, whose entries are Stirling numbers.
    If *kind* is 0, the entries are set to the unsigned Stirling numbers
    of the first kind. If *kind* is 1, the entries are set to the signed
    Stirling numbers of the first kind. If *kind* is 2, the entries are
    set to the Stirling numbers of the second kind.

.. function:: void ca_mat_hilbert(ca_mat_t mat, ca_ctx_t ctx)

    Sets *mat* to the Hilbert matrix, which has entries `A_{i,j} = 1/(i+j+1)`.

.. function:: void ca_mat_dft(ca_mat_t mat, int type, ca_ctx_t ctx)

    Sets *mat* to the DFT (discrete Fourier transform) matrix of order *n*
    where *n* is the smallest dimension of *mat* (if *mat* is not square,
    the matrix is extended periodically along the larger dimension).
    The *type* parameter selects between four different versions
    of the DFT matrix (in which `\omega = e^{2\pi i/n}`):

    * Type 0 -- entries `A_{j,k} = \omega^{-jk}`
    * Type 1 -- entries `A_{j,k} = \omega^{jk} / n`
    * Type 2 -- entries `A_{j,k} = \omega^{-jk} / \sqrt{n}`
    * Type 3 -- entries `A_{j,k} = \omega^{jk} / \sqrt{n}`

    The type 0 and 1 matrices are inverse pairs, and similarly for the
    type 2 and 3 matrices.

Comparisons and properties
-------------------------------------------------------------------------------

.. function:: truth_t ca_mat_check_equal(const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)

    Compares *A* and *B* for equality.

.. function:: truth_t ca_mat_check_is_zero(const ca_mat_t A, ca_ctx_t ctx)

    Tests if *A* is the zero matrix.

.. function:: truth_t ca_mat_check_is_one(const ca_mat_t A, ca_ctx_t ctx)

    Tests if *A* has ones on the main diagonal and zeros elsewhere.

Conjugate and transpose
-------------------------------------------------------------------------------

.. function:: void ca_mat_transpose(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)

    Sets *res* to the transpose of *A*.

.. function:: void ca_mat_conjugate(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)

    Sets *res* to the entrywise complex conjugate of *A*.

.. function:: void ca_mat_conjugate_transpose(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)

    Sets *res* to the conjugate transpose (Hermitian transpose) of *A*.

Arithmetic
-------------------------------------------------------------------------------

.. function:: void ca_mat_neg(ca_mat_t res, const ca_mat_t A, ca_ctx_t ctx)

    Sets *res* to the negation of *A*.

.. function:: void ca_mat_add(ca_mat_t res, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)

    Sets *res* to the sum of *A* and *B*.

.. function:: void ca_mat_sub(ca_mat_t res, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)

    Sets *res* to the difference of *A* and *B*.

.. function:: void ca_mat_mul(ca_mat_t C, const ca_mat_t A, const ca_mat_t B, ca_ctx_t ctx)

    Sets *res* to the matrix product of *A* and *B*.

.. function:: void ca_mat_mul_si(ca_mat_t B, const ca_mat_t A, slong c, ca_ctx_t ctx)
              void ca_mat_mul_fmpz(ca_mat_t B, const ca_mat_t A, const fmpz_t c, ca_ctx_t ctx)
              void ca_mat_mul_fmpq(ca_mat_t B, const ca_mat_t A, const fmpq_t c, ca_ctx_t ctx)
              void ca_mat_mul_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)

    Sets *B* to *A* multiplied by the scalar *c*.

.. function:: void ca_mat_div_si(ca_mat_t B, const ca_mat_t A, slong c, ca_ctx_t ctx)
              void ca_mat_div_fmpz(ca_mat_t B, const ca_mat_t A, const fmpz_t c, ca_ctx_t ctx)
              void ca_mat_div_fmpq(ca_mat_t B, const ca_mat_t A, const fmpq_t c, ca_ctx_t ctx)
              void ca_mat_div_ca(ca_mat_t B, const ca_mat_t A, const ca_t c, ca_ctx_t ctx)

    Sets *B* to *A* divided by the scalar *c*.


Trace
-------------------------------------------------------------------------------

.. function:: void ca_mat_trace(ca_t trace, const ca_mat_t mat, ca_ctx_t ctx)

    Sets *trace* to the sum of the entries on the main diagonal of *mat*.


Gaussian elimination and LU factorization
-------------------------------------------------------------------------------

.. function:: truth_t ca_mat_find_pivot(slong * pivot_row, const ca_mat_t mat, slong start_row, slong end_row, slong column, ca_ctx_t ctx)

    Attempts to find a nonzero entry in *mat* with column index *column*
    and row index between *start_row* (inclusive) and *end_row* (exclusive).

    If the return value is ``T_TRUE``, such an element exists,
    and *pivot_row* is set to the row index.
    If the return value is ``T_FALSE``, no such element exists
    (all entries in this part of the column are zero).
    If the return value is ``T_UNKNOWN``, it is unknown whether such
    an element exists (zero certification failed).

.. function:: truth_t ca_mat_nonsingular_lu(slong * P, ca_mat_t LU, const ca_mat_t A, ca_ctx_t ctx)

    Given a square matrix *A*, attempts to prove invertibility of *A* and
    compute a LU decomposition `PLU = A`
    using Gaussian elimination with partial pivoting.
    The input and output matrices can be the same, performing the
    decomposition in-place.
    Entry `i` in the permutation vector *P* is set to the row index in
    the input matrix corresponding to row `i` in the output matrix.

    * If *A* can be proved to be invertible/nonsingular, returns ``T_TRUE`` and sets *P* and *LU* to a LU decomposition `A = PLU`.

    * If *A* can be proved to be singular, returns ``T_FALSE``.

    * If *A* cannot be proved to be either singular or nonsingular, returns ``T_UNKNOWN``.

    When the return value is ``T_FALSE`` or ``T_UNKNOWN``, the values of
    *P* and *LU* are arbitrary.

.. function:: truth_t ca_mat_nonsingular_fflu(slong * P, ca_mat_t LU, ca_t det, const ca_mat_t A, ca_ctx_t ctx)

    Similar to :func:`ca_mat_nonsingular_lu`, but computes a fraction-free
    LU decomposition using the Bareiss algorithm.
    The denominator is written to *det*.
    Note that despite being "fraction-free", this algorithm may
    introduce fractions due to incomplete symbolic simplifications.


Determinant
-------------------------------------------------------------------------------

.. function:: void ca_mat_det_berkowitz(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
              int ca_mat_det_lu(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
              int ca_mat_det_bareiss(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
              void ca_mat_det_cofactor(ca_t det, const ca_mat_t A, ca_ctx_t ctx)
              void ca_mat_det(ca_t det, const ca_mat_t A, ca_ctx_t ctx)

    Sets *det* to the determinant of the square matrix *A*.
    Various algorithms are available:

    * The *berkowitz* version uses the division-free Berkowitz algorithm
      performing `O(n^4)` operations. Since no zero tests are required, it
      is guaranteed to succeed.

    * The *cofactor* version performs cofactor expansion. This is currently
      only supported for matrices up to size 4.

    * The *lu* and *bareiss* versions use rational LU decomposition
      and fraction-free LU decomposition (Bareiss algorithm) respectively,
      requiring `O(n^3)` operations. These algorithms can fail if zero
      certification fails (see :func:`ca_mat_nonsingular_lu`); they
      return 1 for success and 0 for failure.
      Note that the Bareiss algorithm, despite being "fraction-free",
      may introduce fractions due to incomplete symbolic simplifications.

    The default function chooses an algorithm automatically.
    It will, in addition, recognize trivially rational and integer
    matrices and evaluate those determinants using
    :type:`fmpq_mat_t` or :type:`fmpz_mat_t`.

    The various algorithms can produce different symbolic
    forms of the same determinant. Which algorithm performs better
    depends strongly and sometimes
    unpredictably on the structure of the matrix and the fields of


Characteristic polynomial
-------------------------------------------------------------------------------

.. function:: void _ca_mat_charpoly(ca_ptr cp, const ca_mat_t mat, ca_ctx_t ctx)

.. function:: void ca_mat_charpoly(ca_poly_t cp, const ca_mat_t mat, ca_ctx_t ctx)

    Sets *poly* to the characteristic polynomial of *mat* which must be
    a square matrix. If the matrix has *n* rows, the underscore method
    requires space for `n + 1` output coefficients.
    Employs the division-free Berkowitz algorithm using
    `O(n^4)` operations.

.. function:: int ca_mat_companion(ca_mat_t mat, const ca_poly_t poly, ca_ctx_t ctx)

    Sets *mat* to the companion matrix of *poly*.
    This function verifies that the leading coefficient of *poly*
    is provably nonzero and that the output matrix has the right size,
    returning 1 on success.
    It returns 0 if the leading coefficient of *poly* cannot be
    proved nonzero or if the size of the output matrix does not match.


Eigenvalues and eigenvectors
-------------------------------------------------------------------------------

.. function:: int ca_mat_eigenvalues(ca_vec_t lambda, ulong * exp, ca_mat_t mat, ca_ctx_t ctx)

    Attempts to compute all complex eigenvalues of the given matrix *mat*.
    On success, returns 1 and sets *lambda* to the distinct eigenvalues
    with corresponding multiplicities in *exp*.
    The eigenvalues are returned in arbitrary order.
    On failure, returns 0 and leaves the values in *lambda* and *exp*
    arbitrary.

    This function effectively computes the characteristic polynomial
    and then calls :type:`ca_poly_roots`.
