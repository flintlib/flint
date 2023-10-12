.. _fq-zech-vec:

**fq_zech_vec.h** -- vectors over finite fields (Zech logarithm representation)
===============================================================================

Memory management
--------------------------------------------------------------------------------


.. function:: fq_zech_struct * _fq_zech_vec_init(slong len, const fq_zech_ctx_t ctx)

    Returns an initialised vector of ``fq_zech``'s of given length.

.. function:: void _fq_zech_vec_clear(fq_zech_struct * vec, slong len, const fq_zech_ctx_t ctx)

    Clears the entries of ``(vec, len)`` and frees the space allocated
    for ``vec``.


Randomisation
--------------------------------------------------------------------------------


.. function:: void _fq_zech_vec_randtest(fq_zech_struct * f, flint_rand_t state, slong len, const fq_zech_ctx_t ctx)

    Sets the entries of a vector of the given length to elements of
    the finite field.


Input and output
--------------------------------------------------------------------------------


.. function:: int _fq_zech_vec_fprint(FILE * file, const fq_zech_struct * vec, slong len, const fq_zech_ctx_t ctx)

    Prints the vector of given length to the stream ``file``. The
    format is the length followed by two spaces, then a space separated
    list of coefficients. If the length is zero, only `0` is printed.

    In case of success, returns a positive value.  In case of failure,
    returns a non-positive value.

.. function:: int _fq_zech_vec_print(const fq_zech_struct * vec, slong len, const fq_zech_ctx_t ctx)

    Prints the vector of given length to ``stdout``.

    For further details, see ``_fq_zech_vec_fprint()``.


Assignment and basic manipulation
--------------------------------------------------------------------------------


.. function:: void _fq_zech_vec_set(fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len2, const fq_zech_ctx_t ctx)

    Makes a copy of ``(vec2, len2)`` into ``vec1``.

.. function:: void _fq_zech_vec_swap(fq_zech_struct * vec1, fq_zech_struct * vec2, slong len2, const fq_zech_ctx_t ctx)

    Swaps the elements in ``(vec1, len2)`` and ``(vec2, len2)``.

.. function:: void _fq_zech_vec_zero(fq_zech_struct * vec, slong len, const fq_zech_ctx_t ctx)

    Zeros the entries of ``(vec, len)``.

.. function:: void _fq_zech_vec_neg(fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len2, const fq_zech_ctx_t ctx)

    Negates ``(vec2, len2)`` and places it into ``vec1``.


Comparison
--------------------------------------------------------------------------------


.. function:: int _fq_zech_vec_equal(const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len, const fq_zech_ctx_t ctx)

    Compares two vectors of the given length and returns `1` if they are
    equal, otherwise returns `0`.

.. function:: int _fq_zech_vec_is_zero(const fq_zech_struct * vec, slong len, const fq_zech_ctx_t ctx)

    Returns `1` if ``(vec, len)`` is zero, and `0` otherwise.


Addition and subtraction
--------------------------------------------------------------------------------


.. function:: void _fq_zech_vec_add(fq_zech_struct * res, const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len2, const fq_zech_ctx_t ctx)

    Sets ``(res, len2)`` to the sum of ``(vec1, len2)``
    and ``(vec2, len2)``.

.. function:: void _fq_zech_vec_sub(fq_zech_struct * res, const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len2, const fq_zech_ctx_t ctx)

    Sets ``(res, len2)`` to ``(vec1, len2)`` minus ``(vec2, len2)``.


Scalar multiplication and division
--------------------------------------------------------------------------------

.. function:: void _fq_zech_vec_scalar_addmul_fq_zech(fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len2, const fq_zech_t c, const fq_zech_ctx_t ctx)

    Adds ``(vec2, len2)`` times `c` to ``(vec1, len2)``, where
    `c` is a ``fq_zech_t``.

.. function:: void _fq_zech_vec_scalar_submul_fq_zech(fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len2, const fq_zech_t c, const fq_zech_ctx_t ctx)

    Subtracts ``(vec2, len2)`` times `c` from ``(vec1, len2)``,
    where `c` is a ``fq_zech_t``.


Dot products
--------------------------------------------------------------------------------


.. function:: void _fq_zech_vec_dot(fq_zech_t res, const fq_zech_struct * vec1, const fq_zech_struct * vec2, slong len2, const fq_zech_ctx_t ctx)

    Sets ``res`` to the dot product of (``vec1``, ``len``)
    and (``vec2``, ``len``).
