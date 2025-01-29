.. _double_interval:

**double_interval.h** -- double-precision interval arithmetic and helpers
===============================================================================

This module provides helper functions for computing fast enclosures
using ``double`` arithmetic.


Types, macros and constants
-------------------------------------------------------------------------------

.. type:: di_t

    Holds two ``double`` endpoints ``a`` and ``b`` representing
    the extended real interval `[a, b]`. We generally assume that
    `a \le b` and that neither endpoint is NaN.

Basic manipulation
-------------------------------------------------------------------------------

.. function:: di_t di_interval(double a, double b)

    Returns the interval `[a, b]`. We require that the endpoints
    are ordered and not NaN.

.. function:: di_t arb_get_di(const arb_t x)

    Returns the ball *x* converted to a double-precision interval.

.. function:: void arb_set_di(arb_t res, di_t x, slong prec)

    Sets the ball *res* to the double-precision interval *x*,
    rounded to *prec* bits.

.. function:: void di_print(di_t x)

    Prints *x* to standard output. This simply prints decimal
    representations of the floating-point endpoints; the
    decimals are not guaranteed to be rounded outward.

.. function:: double d_randtest2(flint_rand_t state)

    Returns a random non-NaN ``double`` with any exponent.
    The value can be infinite or subnormal.

.. function:: di_t di_randtest(flint_rand_t state)

    Returns an interval with random endpoints.

Arithmetic
-------------------------------------------------------------------------------

.. function:: di_t di_neg(di_t x)

    Returns the exact negation of *x*.

Fast arithmetic
-------------------------------------------------------------------------------

The following methods perform fast but sloppy interval arithmetic:
we manipulate the endpoints with default rounding and then add
or subtract generic perturbations regardless of whether the
operations were exact.
It is currently assumed that the CPU rounding mode is to nearest.

.. function:: di_t di_fast_add(di_t x, di_t y)
              di_t di_fast_sub(di_t x, di_t y)
              di_t di_fast_mul(di_t x, di_t y)
              di_t di_fast_div(di_t x, di_t y)

    Returns the sum, difference, product or quotient of *x* and *y*.
    Division by zero is currently defined to return `[-\infty, +\infty]`.

.. function:: di_t di_fast_sqr(di_t x)

    Returns the square of *x*. The output is clamped to
    be nonnegative.

.. function:: di_t di_fast_add_d(di_t x, double y)
              di_t di_fast_sub_d(di_t x, double y)
              di_t di_fast_mul_d(di_t x, double y)
              di_t di_fast_div_d(di_t x, double y)

    Arithmetic with an exact ``double`` operand.

.. function:: di_t di_fast_log_nonnegative(di_t x)

    Returns an enclosure of `\log(x)`. The lower endpoint of *x*
    is rounded up to 0 if it is negative.

.. function:: di_t di_fast_mid(di_t x)

    Returns an enclosure of the midpoint of *x*.

.. function:: double di_fast_ubound_radius(di_t x)

    Returns an upper bound for the radius of *x*.
