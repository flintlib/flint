.. _ca-mat:

**ca_mat.h** -- matrices over the real and complex numbers
===============================================================================

An :type:`ca_mat_t` represents a dense matrix over the real numbers,
implemented as an array of entries of type :type:`ca_struct`.
The dimension (number of rows and columns) of a matrix is fixed at
initialization, and the user must ensure that inputs and outputs to
an operation have compatible dimensions. The number of rows or columns
in a matrix can be zero.

Types, macros and constants
-------------------------------------------------------------------------------

.. type:: ca_mat_struct

.. type:: ca_mat_t

    Contains a pointer to a flat array of the entries (entries), an array of
    pointers to the start of each row (rows), and the number of rows (r)
    and columns (c).

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
