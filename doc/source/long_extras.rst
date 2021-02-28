.. _long-extras:

**long_extras.h** -- support functions for signed word arithmetic
===============================================================================


Properties
--------------------------------------------------------------------------------


.. function:: size_t z_sizeinbase(slong n, int b)

    Returns the number of digits in the base `b` representation 
    of the absolute value of the integer `n`.

    Assumes that `b \geq 2`.


Checked Arithmetic
--------------------------------------------------------------------------------

.. function:: int z_mul_checked(slong * a, slong b, slong c)

    Set `*a` to `b` times `c` and return `1` if the product overflowed. Otherwise, return `0`.


Random functions 
--------------------------------------------------------------------------------


.. function:: mp_limb_signed_t z_randtest(flint_rand_t state)

    Returns a pseudo random number with a random number of bits, from 
    `0` to ``FLINT_BITS``.  The probability of the special values `0`, 
    `\pm 1`, ``COEFF_MAX``, ``COEFF_MIN``, ``WORD_MAX`` and 
    ``WORD_MIN`` is increased.

    This random function is mainly used for testing purposes.

.. function:: mp_limb_signed_t z_randtest_not_zero(flint_rand_t state)

    As for ``z_randtest(state)``, but does not return `0`.

.. function:: mp_limb_t z_randint(flint_rand_t state, mp_limb_t limit)

    Returns a pseudo random number of absolute value less than 
    ``limit``.  If ``limit`` is zero or exceeds ``WORD_MAX``, 
    it is interpreted as ``WORD_MAX``.


Modular arithmetic
--------------------------------------------------------------------------------

.. function:: int z_kronecker(slong a, slong n)

    Return the Kronecker symbol `\left(\frac{a}{n}\right)` for any `a` and any `n`.
