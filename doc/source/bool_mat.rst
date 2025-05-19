.. _bool-mat:

**bool_mat.h** -- matrices over booleans
===============================================================================

A :type:`bool_mat_t` represents a dense matrix over the boolean
semiring `\langle \left\{0, 1\right\}, \vee, \wedge \rangle`,
implemented as an array of entries of type ``int``.

The dimension (number of rows and columns) of a matrix is fixed at
initialization, and the user must ensure that inputs and outputs to
an operation have compatible dimensions. The number of rows or columns
in a matrix can be zero.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: bool_mat_struct

.. type:: bool_mat_t

    Contains a pointer to a flat array of the entries (entries), an array of
    pointers to the start of each row (rows), and the number of rows (r)
    and columns (c).

    An *bool_mat_t* is defined as an array of length one of type
    *bool_mat_struct*, permitting an *bool_mat_t* to
    be passed by reference.

.. function:: int bool_mat_get_entry(const bool_mat_t mat, slong i, slong j)

    Returns the entry of matrix *mat* at row *i* and column *j*.

.. function:: void bool_mat_set_entry(bool_mat_t mat, slong i, slong j, int x)

    Sets the entry of matrix *mat* at row *i* and column *j* to *x*.

.. macro:: bool_mat_nrows(mat)

    Returns the number of rows of the matrix.

.. macro:: bool_mat_ncols(mat)

    Returns the number of columns of the matrix.


Memory management
-------------------------------------------------------------------------------

.. function:: void bool_mat_init(bool_mat_t mat, slong r, slong c)

    Initializes the matrix, setting it to the zero matrix with *r* rows
    and *c* columns.

.. function:: void bool_mat_clear(bool_mat_t mat)

    Clears the matrix, deallocating all entries.

.. function:: int bool_mat_is_empty(const bool_mat_t mat)

    Returns nonzero iff the number of rows or the number of columns in *mat*
    is zero. Note that this does not depend on the entry values of *mat*.

.. function:: int bool_mat_is_square(const bool_mat_t mat)

    Returns nonzero iff the number of rows is equal to the number of columns in *mat*.


Conversions
-------------------------------------------------------------------------------

.. function:: void bool_mat_set(bool_mat_t dest, const bool_mat_t src)

    Sets *dest* to *src*. The operands must have identical dimensions.


Input and output
-------------------------------------------------------------------------------

.. function:: void bool_mat_print(const bool_mat_t mat)

    Prints each entry in the matrix.

.. function:: void bool_mat_fprint(FILE * file, const bool_mat_t mat)

    Prints each entry in the matrix to the stream *file*.


Value comparisons
-------------------------------------------------------------------------------

.. function:: int bool_mat_equal(const bool_mat_t mat1, const bool_mat_t mat2)

    Returns nonzero iff the matrices have the same dimensions
    and identical entries.

.. function:: int bool_mat_any(const bool_mat_t mat)

    Returns nonzero iff *mat* has a nonzero entry.

.. function:: int bool_mat_all(const bool_mat_t mat)

    Returns nonzero iff all entries of *mat* are nonzero.

.. function:: int bool_mat_is_diagonal(const bool_mat_t A)

    Returns nonzero iff `i \ne j \implies \bar{A_{ij}}`.

.. function:: int bool_mat_is_lower_triangular(const bool_mat_t A)

    Returns nonzero iff `i < j \implies \bar{A_{ij}}`.

.. function:: int bool_mat_is_transitive(const bool_mat_t mat)

    Returns nonzero iff `A_{ij} \wedge A_{jk} \implies A_{ik}`.

.. function:: int bool_mat_is_nilpotent(const bool_mat_t A)

    Returns nonzero iff some positive matrix power of `A` is zero.


Random generation
-------------------------------------------------------------------------------

.. function:: void bool_mat_randtest(bool_mat_t mat, flint_rand_t state)

    Sets *mat* to a random matrix.

.. function:: void bool_mat_randtest_diagonal(bool_mat_t mat, flint_rand_t state)

    Sets *mat* to a random diagonal matrix.

.. function:: void bool_mat_randtest_nilpotent(bool_mat_t mat, flint_rand_t state)

    Sets *mat* to a random nilpotent matrix.


Special matrices
-------------------------------------------------------------------------------

.. function:: void bool_mat_zero(bool_mat_t mat)

    Sets all entries in mat to zero.

.. function:: void bool_mat_one(bool_mat_t mat)

    Sets the entries on the main diagonal to ones,
    and all other entries to zero.

.. function:: void bool_mat_directed_path(bool_mat_t A)

    Sets `A_{ij}` to `j = i + 1`.
    Requires that `A` is a square matrix.

.. function:: void bool_mat_directed_cycle(bool_mat_t A)

    Sets `A_{ij}` to `j = (i + 1) \mod n`
    where `n` is the order of the square matrix `A`.


Transpose
-------------------------------------------------------------------------------

.. function:: void bool_mat_transpose(bool_mat_t dest, const bool_mat_t src)

    Sets *dest* to the transpose of *src*. The operands must have
    compatible dimensions. Aliasing is allowed for square matrices.


Arithmetic
-------------------------------------------------------------------------------

.. function:: void bool_mat_complement(bool_mat_t B, const bool_mat_t A)

    Sets *B* to the logical complement of *A*.
    That is `B_{ij}` is set to `\bar{A_{ij}}`.
    The operands must have the same dimensions.

.. function:: void bool_mat_add(bool_mat_t res, const bool_mat_t mat1, const bool_mat_t mat2)

    Sets *res* to the sum of *mat1* and *mat2*.
    The operands must have the same dimensions.

.. function:: void bool_mat_mul(bool_mat_t res, const bool_mat_t mat1, const bool_mat_t mat2)

    Sets *res* to the matrix product of *mat1* and *mat2*.
    The operands must have compatible dimensions for matrix multiplication.

.. function:: void bool_mat_mul_entrywise(bool_mat_t res, const bool_mat_t mat1, const bool_mat_t mat2)

    Sets *res* to the entrywise product of *mat1* and *mat2*.
    The operands must have the same dimensions.

.. function:: void bool_mat_sqr(bool_mat_t B, const bool_mat_t A)

   Sets *B* to the matrix square of *A*.
   The operands must both be square with the same dimensions.

.. function:: void bool_mat_pow_ui(bool_mat_t B, const bool_mat_t A, ulong exp)

    Sets *B* to *A* raised to the power *exp*.
    Requires that *A* is a square matrix.


Special functions
-------------------------------------------------------------------------------

.. function:: int bool_mat_trace(const bool_mat_t mat)

    Returns the trace of the matrix, i.e. the sum of entries on the
    main diagonal of *mat*. The matrix is required to be square.
    The sum is in the boolean semiring, so this function returns nonzero iff
    any entry on the diagonal of *mat* is nonzero.

.. function:: slong bool_mat_nilpotency_degree(const bool_mat_t A)

    Returns the nilpotency degree of the `n \times n` matrix *A*.
    It returns the smallest positive `k` such that `A^k = 0`.
    If no such `k` exists then the function returns `-1` if `n` is positive,
    and otherwise it returns `0`.

.. function:: void bool_mat_transitive_closure(bool_mat_t B, const bool_mat_t A)

    Sets *B* to the transitive closure `\sum_{k=1}^\infty A^k`.
    The matrix *A* is required to be square.

.. function:: slong bool_mat_get_strongly_connected_components(slong * p, const bool_mat_t A)

    Partitions the `n` row and column indices of the `n \times n` matrix *A*
    according to the strongly connected components (SCC) of the graph
    for which *A* is the adjacency matrix.
    If the graph has `k` SCCs then the function returns `k`,
    and for each vertex `i \in [0, n-1]`,
    `p_i` is set to the index of the SCC to which the vertex belongs.
    The SCCs themselves can be considered as nodes in a directed acyclic
    graph (DAG), and the SCCs are indexed in postorder with respect to that DAG.

.. function:: slong bool_mat_all_pairs_longest_walk(fmpz_mat_t B, const bool_mat_t A)

    Sets `B_{ij}` to the length of the longest walk with endpoint vertices
    `i` and `j` in the graph whose adjacency matrix is *A*.
    The matrix *A* must be square.  Empty walks with zero length
    which begin and end at the same vertex are allowed.  If `j` is not
    reachable from `i` then no walk from `i` to `j` exists and `B_{ij}`
    is set to the special value `-1`.
    If arbitrarily long walks from `i` to `j` exist then `B_{ij}`
    is set to the special value `-2`.

    The function returns `-2` if any entry of `B_{ij}` is `-2`,
    and otherwise it returns the maximum entry in `B`, except if `A` is empty
    in which case `-1` is returned.
    Note that the returned value is one less than
    that of :func:`nilpotency_degree`.

    This function can help quantify entrywise errors in a truncated evaluation
    of a matrix power series.  If *A* is an indicator matrix with the same
    sparsity pattern as a matrix `M` over the real or complex numbers,
    and if `B_{ij}` does not take the special value `-2`, then the tail
    `\left[ \sum_{k=N}^\infty a_k M^k \right]_{ij}`
    vanishes when `N > B_{ij}`.
