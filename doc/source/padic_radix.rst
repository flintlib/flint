.. _padic-radix:

**padic_radix.h** -- p-adic numbers using radix representation
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

.. macro:: PADIC_RADIX_EXACT

    Special value of `N` used to mark an exact element.

.. macro:: PADIC_RADIX_ERR_MAX

    Upper bound for admissible `N`: `O(x^N)` will automatically be
    clamped to this bound. Also, negative `N` underflowing
    the negation of this value will trigger a ``GR_UNABLE`` status flag.

.. macro:: PADIC_RADIX_PREC_INF

     Special precision value representing an infinite precision `+\infty`.

Context objects
--------------------------------------------------------------------------------

.. function:: int gr_ctx_init_padic_radix(gr_ctx_t ctx, ulong p, slong prec_rel, slong prec_abs, int flags)

    Initialize *ctx* to represent the field of `p`-adic numbers
    with default relative precision *prec_rel* and absolute precision
    *prec_abs* (either or both can be set to ``PADIC_RADIX_PREC_INF``).
    It is required (but not checked) that `p` is a prime number.
    The following flags are supported.

.. macro:: PADIC_RADIX_SIGNED

    Allow signed units for exact elements.

.. macro:: PADIC_RADIX_NO_ERROR

    Floating-point mode: never track error (not implemented).

.. macro:: PADIC_RADIX_DECIMAL

    Print the unit as a decimal integer.

.. macro:: PADIC_RADIX_TEST_LIMITS

    When set, randtest generates `v` and `N` values near the representability limits.

.. function:: void padic_radix_ctx_clear(gr_ctx_t ctx)
              int padic_radix_ctx_write(gr_stream_t out, gr_ctx_t ctx)

    Implementations of standard ``gr_ctx`` methods.

.. macro:: PADIC_RADIX_CTX_RADIX(ctx)

    Macro accessing the :type:`radix_t` context for the internal limb
    arithmetic of a ``padic_radix`` context object.

.. macro:: PADIC_RADIX_CTX_PREC_ABS(ctx)
           PADIC_RADIX_CTX_PREC_REL(ctx)
           PADIC_RADIX_CTX_FLAGS(ctx)
           PADIC_RADIX_CTX_SIGNED(ctx)
           PADIC_RADIX_CTX_DECIMAL(ctx)

    Macros accessing settings and flags of a ``padic_radix`` context object.

Elements
--------------------------------------------------------------------------------

.. type:: padic_radix_struct
          padic_radix_t

    Represents a `p`-adic number.

.. macro:: PADIC_RADIX_UNIT(x)
           PADIC_RADIX_VAL(x)
           PADIC_RADIX_N(x)

    Macros accessing the components of a :type:`padic_radix_t`.

The following methods mainly implement the :ref:`generics <gr>` interface.

.. function:: void padic_radix_init(padic_radix_t res, gr_ctx_t ctx)
              void padic_radix_clear(padic_radix_t res, gr_ctx_t ctx)
              void padic_radix_swap(padic_radix_t x, padic_radix_t y, gr_ctx_t ctx)
              void padic_radix_set_shallow(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)

.. function:: slong padic_radix_get_error(const padic_radix_t x, gr_ctx_t ctx)
              truth_t padic_radix_is_exact(const padic_radix_t x, gr_ctx_t ctx)

.. function:: int padic_radix_randtest(padic_radix_t res, flint_rand_t state, gr_ctx_t ctx)
              int padic_radix_write(gr_stream_t out, const padic_radix_t x, gr_ctx_t ctx)

.. function:: truth_t padic_radix_is_zero(const padic_radix_t x, gr_ctx_t ctx)
              truth_t padic_radix_is_one(const padic_radix_t x, gr_ctx_t ctx)
              truth_t padic_radix_is_neg_one(const padic_radix_t x, gr_ctx_t ctx)
              truth_t padic_radix_equal(const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)

.. function:: int padic_radix_zero(padic_radix_t res, gr_ctx_t ctx)
              int padic_radix_one(padic_radix_t res, gr_ctx_t ctx)
              int padic_radix_set(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
              int padic_radix_set_ui(padic_radix_t res, ulong x, gr_ctx_t ctx)
              int padic_radix_set_si(padic_radix_t res, slong x, gr_ctx_t ctx)
              int padic_radix_set_fmpz(padic_radix_t res, const fmpz_t x, gr_ctx_t ctx)
              int padic_radix_exact_set_ui(padic_radix_t res, ulong x, gr_ctx_t ctx)
              int padic_radix_exact_set_si(padic_radix_t res, slong x, gr_ctx_t ctx)
              int padic_radix_exact_set_fmpz(padic_radix_t res, const fmpz_t x, gr_ctx_t ctx)
              int padic_radix_get_fmpz(fmpz_t res, const padic_radix_t x, gr_ctx_t ctx)
              int padic_radix_neg(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
              int padic_radix_add(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
              int padic_radix_sub(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
              int padic_radix_mul(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
              int padic_radix_inv(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
              int padic_radix_div(padic_radix_t res, const padic_radix_t x, const padic_radix_t y, gr_ctx_t ctx)
              truth_t padic_radix_is_invertible(const padic_radix_t x, gr_ctx_t ctx)
              int padic_radix_sqrt(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
              int padic_radix_rsqrt(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
              truth_t padic_radix_is_square(const padic_radix_t x, gr_ctx_t ctx)

Dot products
--------------------------------------------------------------------------------

.. function:: int padic_radix_dot_strided_delayed(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, slong stride1, const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
              int padic_radix_dot_strided_naive(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, slong stride1, const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)
              int padic_radix_dot_strided(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, slong stride1, const padic_radix_struct * vec2, slong stride2, slong len, gr_ctx_t ctx)

    General strided dot product. The *naive* algorithm just calls :func:`padic_radix_mul`
    and :func:`padic_radix_add` in a loop.
    The *delayed* algorithm computes a global precision and sums the limb cross-products
    from small operands into multiple accumulators, allowing shifts and divisions by
    the digit or limb radix to be delayed until the final recombination step.
    The default algorithm chooses a method automatically.

.. function:: int padic_radix_dot(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, const padic_radix_struct * vec2, slong len, gr_ctx_t ctx)
              int padic_radix_dot_rev(padic_radix_t res, const padic_radix_t initial, int subtract, const padic_radix_struct * vec1, const padic_radix_struct * vec2, slong len, gr_ctx_t ctx)

    Wrappers for :func:`padic_radix_dot_strided`.

Transcendental functions
--------------------------------------------------------------------------------

.. function:: slong _padic_radix_exp_bound(slong v, slong N, ulong p)

    Returns the number of terms `n` such that the tail of the exponential
    series for `x = p^v u` has valuation at least `N`, i.e. the smallest `n`
    with `nv - v_p(n!) \ge N`:

    .. math ::

        n = \max(1, \left \lceil \frac{N(p-1) - 1}{v(p-1) - 1} \right \rceil)

    Assumes (not checked) that `v` is large enough for the series to converge.

.. function:: int _padic_radix_exp_rectangular(radix_integer_t rop, const radix_integer_t u, slong v, slong N, const radix_t radix)
              int _padic_radix_exp_balanced(radix_integer_t rop, const radix_integer_t u, slong v, slong N, const radix_t radix)
              int _padic_radix_exp(radix_integer_t rop, const radix_integer_t u, slong v, slong N, const radix_t radix)

    Given a unit `u` and a valuation `v` such that `v \ge 1` if `p \ne 2`
    and `v \ge 2` if `p = 2` (not checked), sets *rop* to `\exp(u p^v)`
    truncated to precision `p^N`.

    The *rectangular* algorithm uses rectangular splitting.
    The *balanced* algorithm uses the `p`-adic bit-burst algorithm with
    binary splitting. The default function chooses an algorithm automatically.

    The *rectangular* algorithm can return ``GR_UNABLE`` if `N` is too large
    for the implementation; the *balanced* algorithm and the default
    algorithm always succeed and return ``GR_SUCCESS`` (given valid input).

.. function:: int padic_radix_exp_rectangular(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
              int padic_radix_exp_balanced(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)
              int padic_radix_exp(padic_radix_t res, const padic_radix_t x, gr_ctx_t ctx)

    Sets *res* to the `p`-adic exponential function `\exp(x)`.
    Returns ``GR_DOMAIN`` if the series does not converge
    and ``GR_UNABLE`` if this convergence cannot be determined or if the
    precision is too large for the implementation.

