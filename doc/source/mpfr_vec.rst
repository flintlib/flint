.. _mpfr-vec:

**mpfr_vec.h** -- vectors of MPFR floating-point numbers
===============================================================================


Memory management
--------------------------------------------------------------------------------


.. function:: mpfr * _mpfr_vec_init(slong len, mp_bitcnt_t prec)

    Returns a vector of the given length of initialised ``mpfr``'s 
    with the given exact precision.
 
.. function:: void _mpfr_vec_clear(mpfr * vec, slong len)

    Clears the given vector.


Arithmetic
--------------------------------------------------------------------------------


.. function:: void _mpfr_vec_zero(mpfr * vec, slong len)

    Zeros the vector ``(vec, len)``.

.. function:: void _mpfr_vec_set(mpfr * vec1, const mpfr * vec2, slong len)

    Copies the vector ``vec2`` of the given length into ``vec1``. 
    No check is made to ensure ``vec1`` and ``vec2`` are different.

.. function:: void _mpfr_vec_add(mpfr * res, const mpfr * vec1, const mpfr * vec2, slong len)

    Adds the given vectors of the given length together and stores the 
    result in ``res``.

.. function:: void _mpfr_vec_scalar_mul_mpfr(mpfr * res, const mpfr * vec, slong len, mpfr_t c)

    Multiplies the vector with given length by the scalar `c` and 
    sets ``res`` to the result.

.. function:: void _mpfr_vec_scalar_mul_2exp(mpfr * res, const mpfr * vec, slong len, mp_bitcnt_t exp)

    Multiplies the given vector of the given length by ``2^exp``.

.. function:: void _mpfr_vec_scalar_product(mpfr_t res, const mpfr * vec1, const mpfr * vec2, slong len)

   Sets ``res`` to the scalar product of ``(vec1, len)`` with 
    ``(vec2, len)``. Assumes ``len > 0``.
