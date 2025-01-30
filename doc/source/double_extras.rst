.. _double-extras:

**double_extras.h** -- support functions for double arithmetic
===============================================================================

Random functions 
--------------------------------------------------------------------------------


.. function:: double d_randtest(flint_rand_t state)

    Returns a random number in the interval `[0.5, 1)`.

.. function:: double d_randtest_signed(flint_rand_t state, slong minexp, slong maxexp)

    Returns a random signed number with exponent between ``minexp`` and
    ``maxexp`` or zero.

.. function:: double d_randtest_special(flint_rand_t state, slong minexp, slong maxexp)

    Returns a random signed number with exponent between ``minexp`` and
    ``maxexp``, zero, ``D_NAN`` or `\pm`\ ``D_INF``.



Arithmetic
--------------------------------------------------------------------------------


.. function:: double d_polyval(const double * poly, int len, double x)

    Uses Horner's rule to evaluate the polynomial defined by the given
    ``len`` coefficients. Requires that ``len`` is nonzero.


.. function:: double d_mul_2exp_inrange(double x, int i)
              double d_mul_2exp_inrange2(double x, int i)
              double d_mul_2exp(double x, int i)

    Returns `x \cdot 2^i`.

    The *inrange* version requires that `2^i` is in the normal exponent
    range. The *inrange2* version additionally requires that both
    `x` and `x \cdot 2^i` are in the normal exponent range,
    and in particular also assumes that `x \ne 0`.


Special functions
--------------------------------------------------------------------------------


.. function:: double d_lambertw(double x)

    Computes the principal branch of the Lambert W function, solving
    the equation `x = W(x) \exp(W(x))`. If `x < -1/e`, the solution is
    complex, and NaN is returned.

    Depending on the magnitude of `x`, we start from a piecewise rational
    approximation or a zeroth-order truncation of the asymptotic expansion
    at infinity, and perform 0, 1 or 2 iterations with Halley's
    method to obtain full accuracy.

    A test of `10^7` random inputs showed a maximum relative error smaller
    than 0.95 times ``DBL_EPSILON`` (`2^{-52}`) for positive `x`.
    Accuracy for negative `x` is slightly worse, and can grow to
    about 10 times ``DBL_EPSILON`` close to `-1/e`.
    However, accuracy may be worse depending on compiler flags and
    the accuracy of the system libm functions.

.. function:: int d_is_nan(double x)

    Returns a nonzero integral value if ``x`` is ``D_NAN``, and otherwise
    returns 0.

.. function:: double d_log2(double x)

    Returns the base 2 logarithm of ``x`` provided ``x`` is positive. If
    a domain or pole error occurs, the appropriate error value is returned.

