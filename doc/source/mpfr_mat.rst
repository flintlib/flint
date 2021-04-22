.. _mpfr-mat:

**mpfr_mat.h** -- matrices of MPFR floating-point numbers
===============================================================================



Memory management
--------------------------------------------------------------------------------


.. function:: void mpfr_mat_init(mpfr_mat_t mat, slong rows, slong cols, mpfr_prec_t prec)

    Initialises a matrix with the given number of rows and columns and the
    given precision for use. The precision is the exact precision of the
    entries.

.. function:: void mpfr_mat_clear(mpfr_mat_t mat)
 
    Clears the given matrix.


Basic manipulation
--------------------------------------------------------------------------------


.. function:: __mpfr_struct * mpfr_mat_entry(const mpfr_mat_t mat, slong i, slong j)

    Return a reference to the entry at row `i` and column `j` of the given
    matrix. The values `i` and `j` must be within the bounds for the matrix.
    The reference can be used to either return or set the given entry.

.. function:: void mpfr_mat_swap(mpfr_mat_t mat1, mpfr_mat_t mat2)

    Efficiently swap matrices ``mat1`` and ``mat2``.

.. function:: void mpfr_mat_swap_entrywise(mpfr_mat_t mat1, mpfr_mat_t mat2)

    Swaps two matrices by swapping the individual entries rather than swapping
    the contents of the structs.

.. function:: void mpfr_mat_set(mpfr_mat_t mat1, const mpfr_mat_t mat2)

    Set ``mat1`` to the value of ``mat2``.

.. function:: void mpfr_mat_zero(mpfr_mat_t mat)

    Set ``mat`` to the zero matrix.


Comparison
--------------------------------------------------------------------------------


.. function:: int mpfr_mat_equal(const mpfr_mat_t mat1, const mpfr_mat_t mat2)

    Return `1` if the two given matrices are equal, otherwise return `0`.


Randomisation
--------------------------------------------------------------------------------


.. function:: void mpfr_mat_randtest(mpfr_mat_t mat, flint_rand_t state)

    Generate a random matrix with random number of rows and columns and random
    entries for use in test code.


Basic arithmetic
--------------------------------------------------------------------------------


.. function:: void mpfr_mat_mul_classical(mpfr_mat_t C, const mpfr_mat_t A, const mpfr_mat_t B, mpfr_rnd_t rnd)

    Set `C` to the product of `A` and `B` with the given rounding mode, using
    the classical algorithm.
