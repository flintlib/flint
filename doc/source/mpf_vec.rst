.. _mpf-vec:

**mpf_vec.h** -- vectors of MPF floating-point numbers
===============================================================================


Memory management
--------------------------------------------------------------------------------


.. function:: mpf * _mpf_vec_init(slong len)

    Returns a vector of the given length of initialised ``mpf``'s
    with at least the given precision.
 
.. function:: void _mpf_vec_clear(mpf * vec, slong len)

    Clears the given vector.


Randomisation
--------------------------------------------------------------------------------


.. function:: void _mpf_vec_randtest(mpf * f, flint_rand_t state, slong len, flint_bitcnt_t bits)

    Sets the entries of a vector of the given length to random numbers in the 
    interval `[0, 1)` with ``bits`` significant bits in the mantissa or less if
    their precision is smaller.


Assignment and basic manipulation
--------------------------------------------------------------------------------


.. function:: void _mpf_vec_zero(mpf * vec, slong len)

    Zeros the vector ``(vec, len)``.

.. function:: void _mpf_vec_set(mpf * vec1, const mpf * vec2, slong len2)

    Copies the vector ``vec2`` of the given length into ``vec1``. 
    A check is made to ensure ``vec1`` and ``vec2`` are different.
    

Comparison
--------------------------------------------------------------------------------


.. function:: int _mpf_vec_equal(const mpf * vec1, const mpf * vec2, slong len)

    Compares two vectors of the given length and returns `1` if they are 
    equal, otherwise returns `0`.

.. function:: int _mpf_vec_is_zero(const mpf * vec, slong len)

    Returns `1` if ``(vec, len)`` is zero, and `0` otherwise.
    
.. function:: int _mpf_vec_approx_equal(const mpf * vec1, const mpf * vec2, slong len, flint_bitcnt_t bits)

    Compares two vectors of the given length and returns `1` if the first
    ``bits`` bits of their entries are equal, otherwise returns `0`.
    

Addition and subtraction
--------------------------------------------------------------------------------


.. function:: void _mpf_vec_add(mpf * res, const mpf * vec1, const mpf * vec2, slong len2)

    Adds the given vectors of the given length together and stores the 
    result in ``res``.
    
.. function:: void _mpf_vec_sub(mpf * res, const mpf * vec1, const mpf * vec2, slong len2)

    Sets ``(res, len2)`` to ``(vec1, len2)`` minus ``(vec2, len2)``.


Scalar multiplication
--------------------------------------------------------------------------------


.. function:: void _mpf_vec_scalar_mul_mpf(mpf * res, const mpf * vec, slong len, mpf_t c)

    Multiplies the vector with given length by the scalar `c` and 
    sets ``res`` to the result.

.. function:: void _mpf_vec_scalar_mul_2exp(mpf * res, const mpf * vec, slong len, flint_bitcnt_t exp)

    Multiplies the given vector of the given length by ``2^exp``.


Dot product and norm
--------------------------------------------------------------------------------


.. function:: void _mpf_vec_dot(mpf_t res, const mpf * vec1, const mpf * vec2, slong len2)

    Sets ``res`` to the dot product of ``(vec1, len2)`` with 
    ``(vec2, len2)``.
    
.. function:: void _mpf_vec_norm(mpf_t res, mpf * vec, slong len)

    Sets ``res`` to the square of the Euclidean norm of 
    ``(vec, len)``.

.. function:: int _mpf_vec_dot2(mpf_t res, const mpf * vec1, const mpf * vec2, slong len2, flint_bitcnt_t prec)

    Sets ``res`` to the dot product of ``(vec1, len2)`` with
    ``(vec2, len2)``. The temporary variable used has its precision
    set to be at least ``prec`` bits. Returns 0 if a probable
    cancellation is detected, and otherwise returns a non-zero value.

.. function:: void _mpf_vec_norm2(mpf_t res, mpf * vec, slong len, flint_bitcnt_t prec)

    Sets ``res`` to the square of the Euclidean norm of 
    ``(vec, len)``. The temporary variable used has its precision
    set to be at least ``prec`` bits.
