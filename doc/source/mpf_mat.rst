.. _mpf-mat:

**mpf_mat.h** -- matrices of MPF floating-point numbers
===============================================================================



Memory management
--------------------------------------------------------------------------------


.. function:: void mpf_mat_init(mpf_mat_t mat, slong rows, slong cols, flint_bitcnt_t prec)

    Initialises a matrix with the given number of rows and columns and the
    given precision for use. The precision is at least the precision of the
    entries.

.. function:: void mpf_mat_clear(mpf_mat_t mat)
 
    Clears the given matrix.


Basic assignment and manipulation
--------------------------------------------------------------------------------


.. function:: void mpf_mat_set(mpf_mat_t mat1, const mpf_mat_t mat2)

    Sets ``mat1`` to a copy of ``mat2``. The dimensions of 
    ``mat1`` and ``mat2`` must be the same.

.. function:: void mpf_mat_swap(mpf_mat_t mat1, mpf_mat_t mat2)

    Swaps two matrices. The dimensions of ``mat1`` and ``mat2`` 
    are allowed to be different.

.. function:: mpf * mpf_mat_entry(const mpf_mat_t * mat, slong i, slong j)

    Returns a reference to the entry of ``mat`` at row `i` and column `j`.
    Both `i` and `j` must not exceed the dimensions of the matrix.
    The return value can be used to either retrieve or set the given entry.

.. function:: void mpf_mat_zero(mpf_mat_t mat)

    Sets all entries of ``mat`` to 0.

.. function:: void mpf_mat_one(mpf_mat_t mat)

    Sets ``mat`` to the unit matrix, having ones on the main diagonal
    and zeroes elsewhere. If ``mat`` is nonsquare, it is set to the
    truncation of a unit matrix.


Random matrix generation
--------------------------------------------------------------------------------


.. function:: void mpf_mat_randtest(mpf_mat_t mat, flint_ranmpf_t state, flint_bitcnt_t bits)

    Sets the entries of ``mat`` to random numbers in the 
    interval `[0, 1)` with ``bits`` significant bits in the mantissa or less if
    their precision is smaller.


Input and output
--------------------------------------------------------------------------------


.. function:: void mpf_mat_print(const mpf_mat_t mat)

    Prints the given matrix to the stream ``stdout``.


Comparison
--------------------------------------------------------------------------------


.. function:: int mpf_mat_equal(const mpf_mat_t mat1, const mpf_mat_t mat2)

    Returns a non-zero value if ``mat1`` and ``mat2`` have 
    the same dimensions and entries, and zero otherwise.
    
.. function:: int mpf_mat_approx_equal(const mpf_mat_t mat1, const mpf_mat_t mat2, flint_bitcnt_t bits)

    Returns a non-zero value if ``mat1`` and ``mat2`` have 
    the same dimensions and the first ``bits`` bits of their entries
    are equal, and zero otherwise.

.. function:: int mpf_mat_is_zero(const mpf_mat_t mat)

    Returns a non-zero value if all entries ``mat`` are zero, and
    otherwise returns zero.

.. function:: int mpf_mat_is_empty(const mpf_mat_t mat)

    Returns a non-zero value if the number of rows or the number of
    columns in ``mat`` is zero, and otherwise returns
    zero.

.. function:: int mpf_mat_is_square(const mpf_mat_t mat)

    Returns a non-zero value if the number of rows is equal to the
    number of columns in ``mat``, and otherwise returns zero.


Matrix multiplication
--------------------------------------------------------------------------------


.. function:: void mpf_mat_mul(mpf_mat_t C, const mpf_mat_t A, const mpf_mat_t B)

    Sets ``C`` to the matrix product `C = A B`. The matrices must have
    compatible dimensions for matrix multiplication (an exception is raised
    otherwise). Aliasing is allowed.


Gram-Schmidt Orthogonalisation and QR Decomposition
--------------------------------------------------------------------------------


.. function:: void mpf_mat_gso(mpf_mat_t B, const mpf_mat_t A)

    Takes a subset of `R^m` `S = {a_1, a_2, \ldots ,a_n}` (as the columns of
    a `m x n` matrix ``A``) and generates an orthonormal set
    `S' = {b_1, b_2, \ldots ,b_n}` (as the columns of the `m x n` matrix 
    ``B``) that spans the same subspace of `R^m` as `S`.

    This uses an algorithm of Schwarz-Rutishauser. See pp. 9 of
    https://people.inf.ethz.ch/gander/papers/qrneu.pdf
    
.. function:: void mpf_mat_qr(mpf_mat_t Q, mpf_mat_t R, const mpf_mat_t A)

    Computes the `QR` decomposition of a matrix ``A`` using the Gram-Schmidt
    process. (Sets ``Q`` and ``R`` such that `A = QR` where ``R`` is
    an upper triangular matrix and ``Q`` is an orthogonal matrix.)

    This uses an algorithm of Schwarz-Rutishauser. See pp. 9 of
    https://people.inf.ethz.ch/gander/papers/qrneu.pdf
