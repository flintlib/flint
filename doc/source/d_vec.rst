.. _d-vec:

**d_vec.h** -- double precision vectors
===============================================================================


Memory management
--------------------------------------------------------------------------------


.. function:: double * _d_vec_init(slong len)

    Returns an initialised vector of ``double``\s of given length. The
    entries are not zeroed.

.. function:: void _d_vec_clear(double * vec)

    Frees the space allocated for ``vec``.


Randomisation
--------------------------------------------------------------------------------


.. function:: void _d_vec_randtest(double * f, flint_rand_t state, slong len, slong minexp, slong maxexp)

    Sets the entries of a vector of the given length to random signed numbers
    with exponents between ``minexp`` and ``maxexp`` or zero.


Assignment and basic manipulation
--------------------------------------------------------------------------------


.. function:: void _d_vec_set(double * vec1, const double * vec2, slong len2)

    Makes a copy of ``(vec2, len2)`` into ``vec1``.

.. function:: void _d_vec_zero(double * vec, slong len)

    Zeros the entries of ``(vec, len)``.


Comparison
--------------------------------------------------------------------------------


.. function:: int _d_vec_equal(const double * vec1, const double * vec2, slong len)

    Compares two vectors of the given length and returns `1` if they are 
    equal, otherwise returns `0`.

.. function:: int _d_vec_is_zero(const double * vec, slong len)

    Returns `1` if ``(vec, len)`` is zero, and `0` otherwise.
    
.. function:: int _d_vec_is_approx_zero(const double * vec, slong len, double eps)

    Returns `1` if the entries of ``(vec, len)`` are zero to within
    ``eps``, and `0` otherwise.

.. function:: int _d_vec_approx_equal(const double * vec1, const double * vec2, slong len, double eps)

    Compares two vectors of the given length and returns `1` if their entries 
    are within ``eps`` of each other, otherwise returns `0`.


Addition and subtraction
--------------------------------------------------------------------------------


.. function:: void _d_vec_add(double * res, const double * vec1, const double * vec2, slong len2)

    Sets ``(res, len2)`` to the sum of ``(vec1, len2)`` 
    and ``(vec2, len2)``.

.. function:: void _d_vec_sub(double * res, const double * vec1, const double * vec2, slong len2)

    Sets ``(res, len2)`` to ``(vec1, len2)`` minus ``(vec2, len2)``.


Dot product and norm
--------------------------------------------------------------------------------


.. function:: double _d_vec_dot(const double * vec1, const double * vec2, slong len2)

    Returns the dot product of ``(vec1, len2)`` 
    and ``(vec2, len2)``.
    
.. function:: double _d_vec_norm(const double * vec, slong len)

    Returns the square of the Euclidean norm of ``(vec, len)``.
    
.. function:: double _d_vec_dot_heuristic(const double * vec1, const double * vec2, slong len2, double * err)

    Returns the dot product of ``(vec1, len2)``
    and ``(vec2, len2)`` by adding up the positive and negative products,
    and doing a single subtraction of the two sums at the end. ``err`` is a
    pointer to a double in which an error bound for the operation will be
    stored.

.. function:: double _d_vec_dot_thrice(const double * vec1, const double * vec2, slong len2, double * err)

    Returns the dot product of ``(vec1, len2)``
    and ``(vec2, len2)`` using error-free floating point sums and products
    to compute the dot product with three times (thrice) the working precision.
    ``err`` is a pointer to a double in which an error bound for the
    operation will be stored.

    This implements the algorithm of Ogita-Rump-Oishi. See
    http://www.ti3.tuhh.de/paper/rump/OgRuOi05.pdf.
