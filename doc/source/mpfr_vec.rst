.. _mpfr-vec:

**mpfr_vec.h** -- vectors of MPFR floating-point numbers
===============================================================================


Memory management
--------------------------------------------------------------------------------


.. function:: flint_mpfr * _mpfr_vec_init(slong len, flint_bitcnt_t prec)

    Returns a vector of the given length of initialised ``mpfr``'s 
    with the given exact precision.
 
.. function:: void _mpfr_vec_clear(flint_mpfr * vec, slong len)

    Clears the given vector.


Arithmetic
--------------------------------------------------------------------------------


.. function:: void _mpfr_vec_zero(flint_mpfr * vec, slong len)

    Zeros the vector ``(vec, len)``.

.. function:: void _mpfr_vec_set(flint_mpfr * vec1, const flint_mpfr * vec2, slong len)

    Copies the vector ``vec2`` of the given length into ``vec1``. 
    No check is made to ensure ``vec1`` and ``vec2`` are different.

.. function:: void _mpfr_vec_add(flint_mpfr * res, const flint_mpfr * vec1, const flint_mpfr * vec2, slong len)

    Adds the given vectors of the given length together and stores the 
    result in ``res``.

.. function:: void _mpfr_vec_scalar_mul_mpfr(flint_mpfr * res, const flint_mpfr * vec, slong len, mpfr_t c)

    Multiplies the vector with given length by the scalar `c` and 
    sets ``res`` to the result.

.. function:: void _mpfr_vec_scalar_mul_2exp(flint_mpfr * res, const flint_mpfr * vec, slong len, flint_bitcnt_t exp)

    Multiplies the given vector of the given length by ``2^exp``.

.. function:: void _mpfr_vec_scalar_product(mpfr_t res, const flint_mpfr * vec1, const flint_mpfr * vec2, slong len)

   Sets ``res`` to the scalar product of ``(vec1, len)`` with 
    ``(vec2, len)``. Assumes ``len > 0``.
