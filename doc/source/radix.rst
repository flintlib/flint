.. _radix:

**radix.h** -- multiprecision arithmetic in general radix
===============================================================================

This module implements multiprecision arithmetic in
a general radix `2 \le b < \beta` where `\beta = 2^{\mathtt{FLINT\_BITS}}` is the
machine word radix.

The user selects both a *digit radix* `b` and a *limb radix* `B = b^e`
where `B < \beta`, allowing several digits to be packed into each
machine word. A *limb* is a chunk of `e` consecutive digits stored
in one word.
Efficiency is typically maximized by choosing the largest limb radix
that fits in a machine word. For example, for decimal arithmetic one should
typically select `B = 10^{19}` on a 64-bit machine and `B = 10^9` on a 32-bit
machine. A smaller limb radix may however be more efficient in some
circumstances.

Arithmetic in radix `B` is generally slower than arithmetic in the
machine word radix `\beta`, but has the advantage that
base-related operations such as reduction modulo `b^n` and string input/output
in base `b` can be done in time `O(n)` rather than `O(M(n))` or `O(M(n) \log n)`.
The two principal applications are decimal arithmetic (with `b = 10`)
and `p`-adic arithmetic (with `b = p` a small odd prime number).
This module can also be useful for experimenting with multi-limb algorithms
where debugging and discovering corner cases is easier in a small base like
10 than with full-word limbs.

The functions in this module are not optimized for radix `2^e`, where
the machine word radix is applicable (see GMP's ``mpn`` and FLINT's
``flint_mpn`` routines).

Radix objects
--------------------------------------------------------------------------------

.. type :: radix_struct
           radix_t

    A context object defines the radix `B = b^e` and stores additional
    data such as various precomputed powers `b^k` and their inverses used
    for fast division and modular reduction.

.. function:: void radix_init(radix_t radix, ulong b, unsigned int e)

    Initializes context object representing radix `B = b^e` given
    `b \ge 2` and an exponent `e \ge 1` such that `B < \beta`.
    Optionally, the user can input `e = 0`, in which case `e` will
    automatically be set to the largest exponent such that `b^e < \beta`.

.. function:: void radix_clear(radix_t radix)

    Frees any memory allocated by the context object.

.. function:: void radix_init_randtest(radix_t radix, flint_rand_t state)

    Initializes to a randomly chosen radix.

.. function:: ulong radix_digit_radix(const radix_t radix)
              ulong radix_limb_radix(const radix_t radix)
              ulong radix_limb_exponent(const radix_t radix)

    Respectively returns `b`, `B` and `e`.

Low-level natural number arithmetic
--------------------------------------------------------------------------------

A limb has type :type:`ulong`. The type :type:`nn_ptr` denotes a pointer
to an array of writable limbs and :type:`nn_srcptr` denotes a pointer
to an array of readonly limbs.
We use the notation *(x, n)* for the natural number in the range `[0, B^n-1]`
represented by an array of *n* limbs starting at address *x*.
Except where otherwise noted, the following rules apply:

* Input limbs must be normalised to the range `[0, B - 1]`, and output limbs
  will also be normalised to this range.
* Aliasing is permitted as long as aliased operands start at the same address.

.. function:: void radix_rand_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)

    Sets *(res, n)* to a uniformly random integer in `[0, B^n-1]`.

.. function:: void radix_rand_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)

    Sets *(res, rn)* to a uniformly random integer in `[0, b^n-1]`
    where `rn = \lceil n / e \rceil`.

.. function:: void radix_randtest_limbs(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)

    Sets *(res, n)* to a non-uniformly random integer in `[0, B^n-1]`.
    This will produce long strings of limbs with the value 0 or `B-1`
    with high probability.

.. function:: void radix_randtest_digits(nn_ptr res, flint_rand_t state, slong n, const radix_t radix)

    Sets *(res, rn)* to a non-uniformly random integer in `[0, b^n-1]`
    where `rn = \lceil n / e \rceil`.

.. function:: ulong radix_neg(nn_ptr res, nn_srcptr x, slong n, const radix_t radix)

    Sets *(res, n)* to the `B`:s complement negation of *(x, n)*,
    i.e., `(-x) \bmod B^n`, returning borrow.

.. function:: ulong radix_add(nn_ptr res, nn_srcptr x, slong xn, nn_srcptr y, slong yn, const radix_t radix)
              ulong radix_sub(nn_ptr res, nn_srcptr x, slong xn, nn_srcptr y, slong yn, const radix_t radix)

    Sets *(res, xn)* to the sum or difference of *(x, xn)* and *(y, yn)*, returning
    carry or borrow. Requires `xn \ge yn \ge 1`.

.. function:: void radix_mul(nn_ptr res, nn_srcptr x, slong xn, nn_srcptr y, slong yn, const radix_t radix)

    Sets *(res, xn + yn)* to the product of *(x, xn)* and *(y, yn)*.
    Requires `xn \ge yn \ge 1`. Does not allow *res* to be aliased with
    either *x* or *y*.

    This function is currently a wrapper of :func:`radix_mulmid`.

.. function:: void radix_sqr(nn_ptr res, nn_srcptr x, slong xn, const radix_t radix)

    Sets *(res, 2 xn)* to the square of *(x, xn)*. Requires `xn \ge 1`. Does
    not allow *res* to be aliased with *x*.

.. function:: void radix_mulmid(nn_ptr res, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong lo, slong hi, const radix_t radix)
              void radix_mulmid_classical(nn_ptr res, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong lo, slong hi, const radix_t radix)
              void radix_mulmid_fft_small(nn_ptr res, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong lo, slong hi, const radix_t radix)
              void radix_mulmid_KS(nn_ptr res, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong lo, slong hi, const radix_t radix)
              void radix_mulmid_naive(nn_ptr res, nn_srcptr x, slong xn, nn_srcptr y, slong yn, slong lo, slong hi, const radix_t radix)

    Short, truncated or middle product.
    Requires `xn \ge yn \ge 1` and `0 \le lo < hi \le xn + yn`.
    Does not allow *res* to be aliased with either *x* or *y*.

    Viewing as `x` and `y` as polynomials `X, Y \in \mathbb{Z}[T]` evaluated
    at `T = B`, we implicitly compute the polynomial product `Z = XY` and 
    extract the slice `Z_{lo:hi} = \sum_{k=lo}^{hi-1} ([T^k] Z) T^{k-lo}`,
    and write `Z_{lo:hi}(B) \bmod B^{hi-lo}` to *(res, hi - lo)*.
    We omit carries from lower order terms and we discard carries that
    would go into higher order terms.

    In the special case `lo = 0`, this sets *(res, hi)* to
    the low product `xy \bmod B^{hi}`.

    In the special case `hi = xn + yn`, this sets *(res, hi - lo)*
    to an approximate high product `r` satisfying
    `xy \ge r B^{lo} \ge xy - \min(xn,yn,lo) B^{lo+1}`.

    In the special case `lo = 0`, `hi = xn + yn`, this computes the full
    product.

    Various multiplication algorithms are implemented:

    * *classical* is an efficient implementation of the schoolbook algorithm,
      good when either *yn* or `hi - lo` is small.
    * *fft_small* multiplies using FFT. This is only available if FLINT is
      built with support for the *fft_small* module.
    * *KS* denotes Kronecker substitution, converting to a larger integer
      multiplication in the machine word radix.
    * *naive* is an unoptimized reference implementation which does a polynomial
      multiplication using ``fmpz_poly``.

.. function:: ulong radix_divrem_1(nn_ptr res, nn_srcptr x, slong xn, ulong d, const radix_t radix)

    Sets *(res, xn)* to the quotient of *(x, xn)* divided by *d*, returning the
    remainder. Requires `1 \le d \le B - 1`.

.. function:: void radix_divexact_1(nn_ptr res, nn_srcptr x, slong xn, ulong d, const radix_t radix)

    Sets *(res, xn)* to the quotient of *(x, xn)* divided by *d* assuming that
    *(x, xn)* is an exact multiple of *d*, with undefined behavior otherwise.
    Requires `1 \le d \le B - 1`.

.. function:: int radix_invmod_bn(nn_ptr res, nn_srcptr x, slong xn, slong n, const radix_t radix)

    If *(x, xn)* is invertible modulo `B^n`, sets *(res, n)* to the inverse
    and returns 1. If *(x, xn)* is not invertible, returns 0.

.. function:: int radix_cmp_bn_half(nn_srcptr x, slong n, const radix_t radix)

    Returns -1, 0 or 1 according to whether *(x, n)* is less, equal
    or greater than `\lfloor B^n / 2 \rfloor`. We require *n* to be positive.

Radix conversion
--------------------------------------------------------------------------------

.. function:: slong radix_get_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
              slong radix_get_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
              slong radix_get_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)

    Convert *(a, an)* to the machine word radix, writing the output limbs
    to *res* and returning the exact size of the result in the machine word radix.
    Leading zero limbs are omitted from the output size and may or may not be written.
    Requires that *res* has enough space for the largest representable integer
    with *an* limbs in the given radix, plus one extra limb.

.. function:: slong radix_set_mpn_basecase(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
              slong radix_set_mpn_divconquer(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)
              slong radix_set_mpn(nn_ptr res, nn_srcptr a, slong an, const radix_t radix)

    Convert *(a, an)* from the machine word radix, writing the output limbs
    to *res* and returning the exact size of the result in the target radix.
    Leading zero limbs are omitted from the output size and may or may not be written.
    Requires that *res* has space for at least ``radix_set_mpn_need_alloc(an, radix)``
    limbs.

.. function:: slong radix_set_mpn_need_alloc(slong n, const radix_t radix)

    Return a number of output limbs for which :func:`radix_set_mpn` is safe to call
    with input of length `n`.

String conversion
--------------------------------------------------------------------------------

.. function:: char * radix_get_str_decimal(char * res, nn_srcptr x, slong n, int negative, const radix_t radix)

    Returns a null-terminated string representing the integer *(x, n)*
    in base 10 (e.g. ``-12345``). If `n = 0`, the string ``0`` is returned.
    Otherwise, the most significant limb is assumed to be zero; a leading
    minus sign will be inserted if the flag *negative* is 1.

    If the input radix is of the form `10^e`, this function runs in `O(n)`.
    Otherwise, a radix conversion to decimal will be performed internally,
    which runs in `O(M(n) \log n)`.

    If the given pointer *res* is *NULL*, a new string will be allocated
    with :func:`flint_malloc` and returned. Otherwise, *res* will be used,
    and is assumed to have sufficient space to store the result including
    ter.

.. function:: char * radix_get_str_sum(char * res, nn_srcptr x, slong n, int negative, int ascending, const radix_t radix)

    Returns *(x, n)* represented in expression form as a sum over its limbs.
    For example, 12345678 in radix `10^3` will be printed as
    ``12 * 10^6 + 345 * 10^3 + 678`` if *negative* and *ascending*
    flags are both set to 0,
    and as ``-(678 + 345 * 10^3 + 12 * 10^6)`` if the *negative* and
    *ascending* flags are both set to 1.

    If the given pointer *res* is *NULL*, a new string will be allocated
    with :func:`flint_malloc` and returned. Otherwise, *res* will be used,
    and is assumed to have sufficient space to store the result including
    null terminator.

Memory-managed integers
--------------------------------------------------------------------------------

.. type :: radix_integer_struct
           radix_integer_t

    Represents a signed multi-limb integer (with respect to a separate radix
    object), in a format similar to ``mpz``.
    There are three slots: ``d`` (pointer to limbs), ``alloc``, and
    ``size`` whose absolute value stores the normalized limb count (negated to
    indicate a negative value).

.. function:: void gr_ctx_init_radix_integer(gr_ctx_t ctx, ulong b, unsigned int e)

    Initialize a generic ring for working with integers in radix `b^e`
    (where `e = 0` can be provided to select the largest admissible limb radix
    automatically). Elements have type :type:`radix_integer_t`.

.. function:: void radix_integer_init(radix_integer_t res, const radix_t radix)
              void radix_integer_clear(radix_integer_t res, const radix_t radix)
              nn_ptr radix_integer_fit_limbs(radix_integer_t res, slong nlimbs, const radix_t radix)
              void radix_integer_zero(radix_integer_t res, const radix_t radix)
              void radix_integer_rand_limbs(radix_integer_t res, flint_rand_t state, slong n, const radix_t radix)
              void radix_integer_randtest_limbs(radix_integer_t res, flint_rand_t state, slong max_limbs, const radix_t radix)
              void radix_integer_one(radix_integer_t res, const radix_t radix)
              void radix_integer_neg_one(radix_integer_t res, const radix_t radix)
              int radix_integer_is_zero(const radix_integer_t x, const radix_t radix)
              int radix_integer_is_one(const radix_integer_t x, const radix_t radix)
              int radix_integer_is_neg_one(const radix_integer_t x, const radix_t radix)
              int radix_integer_equal(const radix_integer_t x, const radix_integer_t y, const radix_t radix)
              int radix_integer_cmp(const radix_integer_t x, const radix_integer_t y, const radix_t radix)
              int radix_integer_cmpabs(const radix_integer_t x, const radix_integer_t y, const radix_t radix)
              void radix_integer_set(radix_integer_t res, const radix_integer_t x, const radix_t radix)
              void radix_integer_set_ui(radix_integer_t res, ulong x, const radix_t radix)
              void radix_integer_set_si(radix_integer_t res, slong x, const radix_t radix)
              void radix_integer_set_fmpz(radix_integer_t res, const fmpz_t x, const radix_t radix)
              void radix_integer_get_fmpz(fmpz_t res, const radix_integer_t x, const radix_t radix)
              void radix_integer_neg(radix_integer_t res, const radix_integer_t x, const radix_t radix)
              void radix_integer_abs(radix_integer_t res, const radix_integer_t x, const radix_t radix)
              int radix_integer_sgn(const radix_integer_t x, const radix_t radix)
              void radix_integer_add(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix)
              void radix_integer_sub(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix)
              void radix_integer_mul(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, const radix_t radix)

    Direct interface for working with :type:`radix_integer_t`. These can be used
    as an alternative to the generic interface.

.. function:: int radix_integer_is_normalised(const radix_integer_t x, const radix_t radix)

    Returns whether all limbs of *x* are reduced modulo the limb radix
    and the most significant limb is nonzero. If the result is 0,
    operations on *x* have undefined behavior.

.. function:: slong radix_integer_size(const radix_integer_t x, const radix_t radix)

    Number of limbs in *x*, ignoring the sign.

.. function:: ulong radix_integer_get_limb(const radix_integer_t x, slong index, const radix_t radix)

    Returns the limb at position *index* in *x*.
    The index must be nonnegative. The sign of *x* is ignored.

.. function:: void radix_integer_set_limb(radix_integer_t res, const radix_integer_t x, slong index, ulong c, const radix_t radix)

    Sets *res* to *x* with the limb at position *index* replaced by the given
    limb value *c*. The index must be nonnegative and *c* must be smaller
    than the limb radix.
    The sign of *x* is preserved. If *x* is zero and *c* is nonzero,
    the sign of the result will be positive.

.. function:: void radix_integer_lshift_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)

    Sets *res* to *x* multiplied by `B^n`.

.. function:: void radix_integer_rshift_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)

    Sets *res* to *x* divided by `B^n`, truncating the result (rounding to zero).

.. function:: slong radix_integer_valuation_limbs(const radix_integer_t x, const radix_t radix)

    Returns the largest *n* such that `B^n` divides *x* exactly.
    In the special case `x = 0`, returns zero.

.. function:: void radix_integer_trunc_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)

    Sets *res* to *x* with the absolute value reduced modulo `B^n`, preserving
    the sign of *x*.

.. function:: void radix_integer_mod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)

    Sets *res* to *x* reduced modulo `B^n`, returning the unique
    unsigned representative in `[0, B^n-1]`.

.. function:: void radix_integer_smod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)

    Sets *res* to *x* reduced modulo `B^n`, returning the unique
    signed representative in `[-\lfloor B^n/2 \rfloor, \lfloor B^n/2 \rfloor]`
    if *B* is odd and in `[-B^n/2, B^n/2-1]` if *B* is even.

.. function:: void radix_integer_mullow_limbs(radix_integer_t res, const radix_integer_t x, const radix_integer_t y, slong n, const radix_t radix)

    Sets *res* to the product of *x* and *y* with the absolute value
    reduced modulo `B^n`.

.. function:: int radix_integer_invmod_limbs(radix_integer_t res, const radix_integer_t x, slong n, const radix_t radix)

    If *x* is invertible modulo `B^n`, sets *res* to the inverse
    and returns 1. The sign of the input is preserved.
    If *x* is not invertible, returns 0.

