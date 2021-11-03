.. _padic-mat:

**padic_mat.h** -- matrices over p-adic numbers
===============================================================================

Module documentation
--------------------------------------------------------------------------------

We represent a matrix over `\mathbf{Q}_p` as a product `p^v M`, 
where `p` is a prime number, `v \in \mathbf{Z}` and `M` a matrix 
over `\mathbf{Z}`.
We say this matrix is in *canonical form* if either `M` is zero, 
in which case we choose `v = 0`, too, or if `M` contains at least one 
`p`-adic unit.
We say this matrix is *reduced* modulo `p^N` if it is canonical 
form and if all coefficients of `M` lie in the range `[0, p^{N-v})`.


Macros
--------------------------------------------------------------------------------


.. function:: fmpz_mat_struct * padic_mat(const padic_mat_t A)

    Returns a pointer to the unit part of the matrix, which 
    is a matrix over `\mathbf{Z}`.

    The return value can be used as an argument to 
    the functions in the ``fmpz_mat`` module.

.. function:: fmpz * padic_mat_entry(const padic_mat_t A, slong i, slong j)

    Returns a pointer to unit part of the entry in position `(i, j)`. 
    Note that this is not necessarily a unit.

    The return value can be used as an argument to 
    the functions in the ``fmpz`` module.

.. function:: slong padic_mat_val(const padic_mat_t A)

    Allow access (as L-value or R-value) to ``val`` field of `A`.
    This function is implemented as a macro.

.. function:: slong padic_mat_prec(const padic_mat_t A)

    Allow access (as L-value or R-value) to ``prec`` field of `A`.
    This function is implemented as a macro.

.. function:: slong padic_mat_get_val(const padic_mat_t A)

    Returns the valuation of the matrix.

.. function:: slong padic_mat_get_prec(const padic_mat_t A)

    Returns the `p`-adic precision of the matrix.

.. function:: slong padic_mat_nrows(const padic_mat_t A)

    Returns the number of rows of the matrix `A`.

.. function:: slong padic_mat_ncols(const padic_mat_t A)

    Returns the number of columns of the matrix `A`.


Memory management
--------------------------------------------------------------------------------


.. function:: void padic_mat_init(padic_mat_t A, slong r, slong c)

    Initialises the matrix `A` as a zero matrix with the specified numbers 
    of rows and columns and precision ``PADIC_DEFAULT_PREC``.

.. function:: void padic_mat_init2(padic_mat_t A, slong r, slong c, slong prec)

    Initialises the matrix `A` as a zero matrix with the specified numbers 
    of rows and columns and the given precision.

.. function:: void padic_mat_clear(padic_mat_t A)

    Clears the matrix `A`.

.. function:: void _padic_mat_canonicalise(padic_mat_t A, const padic_ctx_t ctx)

    Ensures that the matrix `A` is in canonical form.

.. function:: void _padic_mat_reduce(padic_mat_t A, const padic_ctx_t ctx)

    Ensures that the matrix `A` is reduced modulo `p^N`, 
    assuming that it is in canonical form already.

.. function:: void padic_mat_reduce(padic_mat_t A, const padic_ctx_t ctx)

    Ensures that the matrix `A` is reduced modulo `p^N`, 
    without assuming that it is necessarily in canonical form.

.. function:: int padic_mat_is_empty(const padic_mat_t A)

    Returns whether the matrix `A` is empty, that is, 
    whether it has zero rows or zero columns.

.. function:: int padic_mat_is_square(const padic_mat_t A)

    Returns whether the matrix `A` is square.

.. function:: int padic_mat_is_canonical(const padic_mat_t A, const fmpz_t p)

    Returns whether the matrix `A` is in canonical form.


Basic assignment
--------------------------------------------------------------------------------


.. function:: void padic_mat_set(padic_mat_t B, const padic_mat_t A)

    Sets `B` to a copy of `A`, respecting the precision of `B`.

.. function:: void padic_mat_swap(padic_mat_t A, padic_mat_t B)

    Swaps the two matrices `A` and `B`.  This is done efficiently by 
    swapping pointers.

.. function:: void padic_mat_swap_entrywise(padic_mat_t mat1, padic_mat_t mat2)

    Swaps two matrices by swapping the individual entries rather than swapping
    the contents of the structs.

.. function:: void padic_mat_zero(padic_mat_t A)

    Sets the matrix `A` to zero.

.. function:: void padic_mat_one(padic_mat_t A)

    Sets the matrix `A` to the identity matrix.  If the precision 
    is negative then the matrix will be the zero matrix.


Conversions
--------------------------------------------------------------------------------


.. function:: void padic_mat_set_fmpq_mat(padic_mat_t B, const fmpq_mat_t A, const padic_ctx_t ctx)

    Sets the `p`-adic matrix `B` to the rational matrix `A`, reduced 
    according to the given context.

.. function:: void padic_mat_get_fmpq_mat(fmpq_mat_t B, const padic_mat_t A, const padic_ctx_t ctx)

    Sets the rational matrix `B` to the `p`-adic matrices `A`;  
    no reduction takes place.


Entries
--------------------------------------------------------------------------------

Because of the choice of the data structure, representing the matrix 
as `p^v M`, setting an entry of the matrix might lead to changes in 
all entries in the matrix `M`.  Also, a specific entry is not readily 
available as a `p`-adic number.  Thus, there are separate functions 
available for getting and setting entries.

.. function:: void padic_mat_get_entry_padic(padic_t rop, const padic_mat_t op, slong i, slong j, const padic_ctx_t ctx)

    Sets ``rop`` to the entry in position `(i, j)` in the matrix ``op``.

.. function:: void padic_mat_set_entry_padic(padic_mat_t rop, slong i, slong j, const padic_t op, const padic_ctx_t ctx)

    Sets the entry in position `(i, j)` in the matrix to ``rop``.


Comparison
--------------------------------------------------------------------------------


.. function:: int padic_mat_equal(const padic_mat_t A, const padic_mat_t B)

    Returns whether the two matrices `A` and `B` are equal.

.. function:: int padic_mat_is_zero(const padic_mat_t A)

    Returns whether the matrix `A` is zero.


Input and output
--------------------------------------------------------------------------------


.. function:: int padic_mat_fprint(FILE * file, const padic_mat_t A, const padic_ctx_t ctx)

    Prints a simple representation of the matrix `A` to the 
    output stream ``file``.  The format is the number of rows, 
    a space, the number of columns, two spaces, followed by a list 
    of all the entries, one row after the other.

    In the current implementation, always returns `1`.

.. function:: int padic_mat_fprint_pretty(FILE * file, const padic_mat_t A, const padic_ctx_t ctx)

    Prints a *pretty* representation of the matrix `A` 
    to the output stream ``file``.  

    In the current implementation, always returns `1`.

.. function:: int padic_mat_print(const padic_mat_t A, const padic_ctx_t ctx)
              int padic_mat_print_pretty(const padic_mat_t A, const padic_ctx_t ctx)


Random matrix generation
--------------------------------------------------------------------------------


.. function:: void padic_mat_randtest(padic_mat_t A, flint_rand_t state, const padic_ctx_t ctx)

    Sets `A` to a random matrix.

    The valuation will be in the range `[- \lceil N/10\rceil, N)`, 
    `[N - \lceil -N/10\rceil, N)`, or `[-10, 0)` as `N` is positive, 
    negative or zero.


Transpose
--------------------------------------------------------------------------------


.. function:: void padic_mat_transpose(padic_mat_t B, const padic_mat_t A)

    Sets `B` to `A^t`.


Addition and subtraction
--------------------------------------------------------------------------------


.. function:: void _padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

    Sets `C` to the exact sum `A + B`, ensuring that the result is in 
    canonical form.

.. function:: void padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

    Sets `C` to the sum `A + B` modulo `p^N`.

.. function:: void _padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

    Sets `C` to the exact difference `A - B`, ensuring that the result is in 
    canonical form.

.. function:: void padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

    Sets `C` to `A - B`, ensuring that the result is reduced.

.. function:: void _padic_mat_neg(padic_mat_t B, const padic_mat_t A)

    Sets `B` to `-A` in canonical form.

.. function:: void padic_mat_neg(padic_mat_t B, const padic_mat_t A, const padic_ctx_t ctx)

    Sets `B` to `-A`, ensuring the result is reduced.


Scalar operations
--------------------------------------------------------------------------------


.. function:: void _padic_mat_scalar_mul_padic(padic_mat_t B, const padic_mat_t A, const padic_t c, const padic_ctx_t ctx)

    Sets `B` to `c A`, ensuring that the result is in canonical form.

.. function:: void padic_mat_scalar_mul_padic(padic_mat_t B, const padic_mat_t A, const padic_t c, const padic_ctx_t ctx)

    Sets `B` to `c A`, ensuring that the result is reduced.

.. function:: void _padic_mat_scalar_mul_fmpz(padic_mat_t B, const padic_mat_t A, const fmpz_t c, const padic_ctx_t ctx)

    Sets `B` to `c A`, ensuring that the result is in canonical form.

.. function:: void padic_mat_scalar_mul_fmpz(padic_mat_t B, const padic_mat_t A, const fmpz_t c, const padic_ctx_t ctx)

    Sets `B` to `c A`, ensuring that the result is reduced.

.. function:: void padic_mat_scalar_div_fmpz(padic_mat_t B, const padic_mat_t A, const fmpz_t c, const padic_ctx_t ctx)

    Sets `B` to `c^{-1} A`, assuming that `c \neq 0`.  
    Ensures that the result `B` is reduced.


Multiplication
--------------------------------------------------------------------------------


.. function:: void _padic_mat_mul(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

    Sets `C` to the product `A B` of the two matrices `A` and `B`, 
    ensuring that `C` is in canonical form.

.. function:: void padic_mat_mul(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

    Sets `C` to the product `A B` of the two matrices `A` and `B`, 
    ensuring that `C` is reduced.
