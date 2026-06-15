.. _radix-padic:

**radix_padic.h** -- p-adic numbers using radix representation
===============================================================================

This module implements multiprecision `p`-adic numbers for word-size `p`
using base `p^e` limbs building on the :ref:`radix <radix>` module.
By default the limb radix is chosen as large as possible,
e.g. `p = 7` will use radix `7^{22}` internally on a 64-bit machine.
The limb radix is mainly an implementation detail: most user-facing
functions measure the precision in digits.

This module is designed to use the :ref:`generics <gr>` interface.
As such, the ring is represented by a :type:`gr_ctx_t` context object,
methods return status flags (``GR_SUCCESS``, ``GR_UNABLE``, ``GR_DOMAIN``),
and one can use generic structures such as :type:`gr_poly_t` for
polynomials and :type:`gr_mat_t` for matrices.

Representation of numbers
--------------------------------------------------------------------------------

A radix `p`-adic number is represented as a ball `u p^v + O(p^N)`.
The unit `u` is stored as a :type:`radix_integer_t` in the limb radix `B = p^e`.
The valuation `v` and the accuracy `N` are measured in powers of `p` (digits),
exactly as in the ``padic`` module.
We support exact elements (without error term).

A nonzero unit `u` is always canonicalised to be non-divisible by `p`.
As a special case, we always set `v = 0` when `u = 0`.

The functions in this module support but are not currently
optimized for `p = 2`.

Precision
--------------------------------------------------------------------------------

We support both *absolute* and *relative* precision.
If the relative precision is *r* and the absolute precision is *a*,
`u p^v` will be truncated to order `O(p^{v + r})` and `O(p^{a})`.
Basic operations track exactly whether a truncation has occurred and will
not introduce an `O(x^N)` term as long as the result is exact.
One can set `a = +\infty` to work with relative precision only
and `r = +\infty` to work with absolute precision only.
If both are set to `+\infty`, we effectively restrict all arithmetic to exact
operations in the subring `\mathbb{Z}[1/p] \subset \mathbb{Q}_p`.

.. macro:: RADIX_PADIC_EXACT

    Special value of `N` used to mark an exact element.

.. macro:: RADIX_PADIC_ERR_MAX

    Upper bound for admissible `N`: `O(x^N)` will automatically be
    clamped to this bound. Also, negative `N` underflowing
    the negation of this value will trigger a ``GR_UNABLE`` status flag.

.. macro:: RADIX_PADIC_PREC_INF

     Special precision value representing an infinite precision `+\infty`.

Context objects
--------------------------------------------------------------------------------

.. function:: int gr_ctx_init_radix_padic(gr_ctx_t ctx, ulong p, slong prec_rel, slong prec_abs, int flags)

    Initialize *ctx* to represent the field of `p`-adic numbers
    with default relative precision *prec_rel* and absolute precision
    *prec_abs* (either or both can be set to ``RADIX_PADIC_PREC_INF``).
    It is required (but not checked) that `p` is a prime number.
    The following flags are supported.

.. macro:: RADIX_PADIC_SIGNED

    Allow signed units for exact elements.

.. macro:: RADIX_PADIC_NO_ERROR

    Floating-point mode: never track error (not implemented).

.. macro:: RADIX_PADIC_DECIMAL

    Print the unit as a decimal integer.

.. macro:: RADIX_PADIC_TEST_LIMITS

    When set, randtest generates `v` and `N` values near the representability limits.

.. function:: void radix_padic_ctx_clear(gr_ctx_t ctx)
              int radix_padic_ctx_write(gr_stream_t out, gr_ctx_t ctx)

    Implementations of standard ``gr_ctx`` methods.

.. macro:: RADIX_PADIC_CTX_RADIX(ctx)

    Macro accessing the :type:`radix_t` context for the internal limb
    arithmetic of a ``radix_padic`` context object.

.. macro:: RADIX_PADIC_CTX_PREC_ABS(ctx)
           RADIX_PADIC_CTX_PREC_REL(ctx)
           RADIX_PADIC_CTX_FLAGS(ctx)
           RADIX_PADIC_CTX_SIGNED(ctx)
           RADIX_PADIC_CTX_DECIMAL(ctx)

    Macros accessing settings and flags of a ``radix_padic`` context object.

Elements
--------------------------------------------------------------------------------

.. type:: radix_padic_struct
          radix_padic_t

    Represents a `p`-adic number.

.. macro:: RADIX_PADIC_UNIT(x)
           RADIX_PADIC_VAL(x)
           RADIX_PADIC_N(x)

    Macros accessing the components of a :type:`radix_padic_t`.

The following methods mainly implement the :ref:`generics <gr>` interface.

.. function:: void radix_padic_init(radix_padic_t res, gr_ctx_t ctx)
              void radix_padic_clear(radix_padic_t res, gr_ctx_t ctx)
              void radix_padic_swap(radix_padic_t x, radix_padic_t y, gr_ctx_t ctx)
              void radix_padic_set_shallow(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)

.. function:: slong radix_padic_get_error(const radix_padic_t x, gr_ctx_t ctx)
              truth_t radix_padic_is_exact(const radix_padic_t x, gr_ctx_t ctx)

.. function:: int radix_padic_randtest(radix_padic_t res, flint_rand_t state, gr_ctx_t ctx)
              int radix_padic_write(gr_stream_t out, const radix_padic_t x, gr_ctx_t ctx)

.. function:: truth_t radix_padic_is_zero(const radix_padic_t x, gr_ctx_t ctx)
              truth_t radix_padic_is_one(const radix_padic_t x, gr_ctx_t ctx)
              truth_t radix_padic_is_neg_one(const radix_padic_t x, gr_ctx_t ctx)
              truth_t radix_padic_equal(const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)

.. function:: int radix_padic_zero(radix_padic_t res, gr_ctx_t ctx)
              int radix_padic_one(radix_padic_t res, gr_ctx_t ctx)
              int radix_padic_set(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
              int radix_padic_set_ui(radix_padic_t res, ulong x, gr_ctx_t ctx)
              int radix_padic_set_si(radix_padic_t res, slong x, gr_ctx_t ctx)
              int radix_padic_set_fmpz(radix_padic_t res, const fmpz_t x, gr_ctx_t ctx)
              int radix_padic_exact_set_ui(radix_padic_t res, ulong x, gr_ctx_t ctx)
              int radix_padic_exact_set_si(radix_padic_t res, slong x, gr_ctx_t ctx)
              int radix_padic_exact_set_fmpz(radix_padic_t res, const fmpz_t x, gr_ctx_t ctx)
              int radix_padic_get_fmpz(fmpz_t res, const radix_padic_t x, gr_ctx_t ctx)
              int radix_padic_neg(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
              int radix_padic_add(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
              int radix_padic_sub(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
              int radix_padic_mul(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
              int radix_padic_inv(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
              int radix_padic_div(radix_padic_t res, const radix_padic_t x, const radix_padic_t y, gr_ctx_t ctx)
              truth_t radix_padic_is_invertible(const radix_padic_t x, gr_ctx_t ctx)
              int radix_padic_sqrt(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
              int radix_padic_rsqrt(radix_padic_t res, const radix_padic_t x, gr_ctx_t ctx)
              truth_t radix_padic_is_square(const radix_padic_t x, gr_ctx_t ctx)

