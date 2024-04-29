.. _mpfr-vec:

**mpfr_vec.h** -- vectors of MPFR floating-point numbers
===============================================================================

.. note::

    This module is deprecated, and will probably removed completely in the
    future.


Memory management
--------------------------------------------------------------------------------


.. function:: mpfr_ptr _mpfr_vec_init(slong len, flint_bitcnt_t prec)

    Returns a vector of the given length of initialised ``mpfr``'s 
    with the given exact precision.
 
.. function:: void _mpfr_vec_clear(mpfr_ptr vec, slong len)

    Clears the given vector.


Arithmetic
--------------------------------------------------------------------------------


.. function:: void _mpfr_vec_zero(mpfr_ptr vec, slong len)

    Zeros the vector ``(vec, len)``.

.. function:: void _mpfr_vec_set(mpfr_ptr vec1, mpfr_srcptr vec2, slong len)

    Copies the vector ``vec2`` of the given length into ``vec1``. 
    No check is made to ensure ``vec1`` and ``vec2`` are different.
