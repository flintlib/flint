.. _decimal:

**decimal.h** -- decimal floating-point and ball arithmetic
===============================================================================

This module implements arbitrary-precision decimal floating-point
arithmetic (:type:`decfloat_t`), a bounded-precision nonnegative radius
type (:type:`decmag_t`), decimal ball arithmetic (:type:`decball_t`),
and the corresponding complex types (:type:`deccfloat_t` and
:type:`deccball_t`, pairs of real and imaginary parts),
building on the :ref:`radix <radix>` module with base `B = 10^e` limbs
(by default `e = 19` on a 64-bit machine and `e = 9` on a 32-bit machine).

This module is designed to use the :ref:`generics <gr>` interface.
Each of the four types is represented by a
:type:`gr_ctx_t` context object which stores the radix, the default
precision (measured in decimal digits), the default rounding mode
(separately for real and imaginary parts, for complex floating-point
numbers), the radius precision, exponent limits and various flags.
Methods return status flags (``GR_SUCCESS``, ``GR_UNABLE``, ``GR_DOMAIN``),
and one can use generic structures such as :type:`gr_poly_t` for
polynomials and :type:`gr_mat_t` for matrices.

Representation of numbers
--------------------------------------------------------------------------------

A finite nonzero decimal floating-point number is represented as

.. math::

    x = (-1)^s \cdot M \cdot B^v, \quad B = 10^e,

where the mantissa `M` is a nonnegative integer stored as
:type:`radix_integer_t` limbs in the limb radix `B`, the sign `s`
is stored in the sign of the :type:`radix_integer_t`, and
`v` is an :type:`fmpz_t` exponent which can be arbitrarily large.
Nonzero mantissas are kept canonical by stripping zero limbs
at both ends (`M_0 \ne 0` and `M_{n-1} \ne 0`), so that equal
values have identical representations.

Precision is measured in decimal digits, independently of the limb
size: a number has `\operatorname{sd}(x) = \operatorname{digits}(M) - \nu_{10}(M_0)`
*significant digits*, where the trailing zero digits of the lowest limb
do not count (they are an artefact of the limb alignment). All rounding
operations produce results with `\operatorname{sd}(x) \le \mathrm{prec}`
in the requested rounding mode, i.e. exactly the results one would
obtain with the classical
definition of correctly rounded decimal arithmetic to `\mathrm{prec}`
significant digits.

**Limb-aligned versus digit-aligned exponents.**
Counting the exponent in limbs (rather than in digits, with the mantissa
normalised to a most or least significant digit as in ``arf`` and
``padic_radix``) means that:

* addition never requires a digit shift: operands are aligned by a limb
  offset, which is free;
* rounding to `\mathrm{prec}` digits after any operation costs `O(1)`:
  whole limbs are dropped and at most `e - 1` digits are cleared in a
  single boundary limb using a precomputed division;
* the digit exponent of a number is still available in `O(1)` since the
  digit count of the top limb is `O(1)`, so comparisons and the
  user-facing semantics (scientific exponent, ulp, rounding to a number
  of significant digits) do not depend on the internal alignment.

The price is that a `\mathrm{prec}`-digit number may occupy one limb
more than `\lceil \mathrm{prec}/e \rceil` since its significant digits
need not be aligned with a limb boundary; for example, with `e = 19`,
`1/3` to 30 digits occupies 2 limbs while `(1/3) \cdot 10^5` occupies 3.
This costs a constant factor (at most 2 for one-limb precision) in the
multiplication of very short numbers, and is negligible for long numbers.
Setting `e = 1` recovers a digit-aligned representation (this is mainly
useful for testing).

Special values `+\infty`, `-\infty` and NaN are encoded by a zero
mantissa with a code in the exponent field; whether they are admitted
is controlled by context flags. Zero is unsigned.

Contexts
--------------------------------------------------------------------------------

.. type:: decimal_ctx_struct

    Context data shared by the real and complex floating-point and ball
    rings. Its fields
    are accessed through the macros ``DECIMAL_CTX_RADIX``,
    ``DECIMAL_CTX_PREC``, ``DECIMAL_CTX_RND``, ``DECIMAL_CTX_RND_IM``,
    ``DECIMAL_CTX_FLAGS``,
    ``DECIMAL_CTX_RAD_PREC``, ``DECIMAL_CTX_EMIN``, ``DECIMAL_CTX_EMAX``,
    ``DECIMAL_CTX_E`` (the limb exponent `e`) and ``DECIMAL_CTX_B``
    (the limb radix `B`). The macros ``DECIMAL_CTX_WHICH``,
    ``DECIMAL_CTX_IS_BALL`` and ``DECIMAL_CTX_IS_COMPLEX`` identify
    the type of the elements.

.. function:: void gr_ctx_init_decfloat(gr_ctx_t ctx, slong prec, int flags)
              void gr_ctx_init_decball(gr_ctx_t ctx, slong prec, int flags)
              void gr_ctx_init_deccfloat(gr_ctx_t ctx, slong prec, int flags)
              void gr_ctx_init_deccball(gr_ctx_t ctx, slong prec, int flags)

    Initializes *ctx* to the ring of decimal floating-point numbers,
    decimal balls, complex decimal floating-point numbers, respectively
    complex decimal balls, with default precision *prec* digits
    (at least 1) and the given *flags*. The default rounding mode is
    ``DECIMAL_RND_NEAR`` (for both real and imaginary parts)
    and the radius precision is
    ``DECMAG_DEFAULT_PREC`` digits. Exponents are unbounded.
    The largest possible limb radix is used.

    Passing ``DECIMAL_PREC_EXACT`` as the precision gives the exact ring
    `\mathbb{Z}[1/10]` of decimal fractions: no rounding is ever performed,
    and inexact operations (division by 3, square roots of non-squares)
    return ``GR_UNABLE``. Additions of operands whose exponents are
    astronomically far apart also return ``GR_UNABLE`` in this mode
    rather than allocating huge mantissas.

.. function:: void _gr_ctx_init_decimal(gr_ctx_t ctx, decimal_ctx_which which, unsigned int e, slong prec, int rnd, int flags)

    Generalized initialization where *which* is ``DECIMAL_CTX_FLOAT``,
    ``DECIMAL_CTX_BALL``, ``DECIMAL_CTX_CFLOAT`` or ``DECIMAL_CTX_CBALL``,
    *e* is the limb exponent (`B = 10^e`,
    or `e = 0` to select the largest possible), and *rnd* is the
    default rounding mode. Elements of contexts with different
    limb exponents are not interchangeable (but can be converted
    with ``gr_set_other``).

.. function:: void gr_ctx_init_decfloat_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec)
              void gr_ctx_init_decball_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec)
              void gr_ctx_init_deccfloat_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec)
              void gr_ctx_init_deccball_randtest(gr_ctx_t ctx, flint_rand_t state, slong max_prec)

    Initializes a random context (random limb exponent, precision,
    rounding modes, radius precision, flags and exponent limits).

.. function:: void decimal_ctx_set_prec(gr_ctx_t ctx, slong prec)
              void decimal_ctx_set_rnd(gr_ctx_t ctx, int rnd)
              void decimal_ctx_set_rnd_im(gr_ctx_t ctx, int rnd)
              void decimal_ctx_set_rad_prec(gr_ctx_t ctx, slong rad_prec)
              void decimal_ctx_set_exp_limits(gr_ctx_t ctx, slong emin, slong emax)
              void decimal_ctx_set_flags(gr_ctx_t ctx, int flags)

    Changes the default precision, the default rounding mode, the
    radius precision (clamped to `[1, R]` where ``R = DECMAG_MAX_PREC``
    is 9 on a 64-bit machine and 4 on a 32-bit machine), the exponent
    limits and the
    flags of a context. Exponent limits apply to the scientific
    exponent `E` of a number `x`, `10^E \le |x| < 10^{E+1}`;
    ``WORD_MIN`` and ``WORD_MAX`` denote the absence of a limit.
    All settings can be changed freely at any time: see
    *Operands and results* below.

    In a complex floating-point context, real parts are rounded in the
    mode set by ``decimal_ctx_set_rnd``, which also resets the mode for
    imaginary parts, and ``decimal_ctx_set_rnd_im`` subsequently
    changes the mode for imaginary parts alone. The imaginary rounding
    mode is ignored in the other contexts.

    Through the generics interface, ``gr_ctx_has_real_prec`` is true and
    ``gr_ctx_set_real_prec`` and ``gr_ctx_get_real_prec`` are supported
    with the precision measured in *bits*: setting `p` bits selects
    `\lceil p \log_{10} 2 \rceil` digits, and a context with `d` digits
    reports `\lfloor d \log_2 10 \rfloor - 1` bits (directed rounding to
    `d` digits has a relative error up to `10^{1-d}`, slightly more than
    that of `d \log_2 10` bits). This is a compromise which lets the
    generic algorithms taking a precision in bits (series expansions,
    linear algebra with adaptive precision, etc.) work with decimal
    contexts; those algorithms assume radix-2 semantics such as exact
    representability of `2^{-p}` and error bounds in powers of two, so
    their results need not be precisely what they would be for binary
    types (they remain correct, since all decimal operations return
    correctly rounded results or rigorous enclosures).

Operands and results
................................................................................

Operands of an operation need not have been produced by the context (or
by a context with the same settings) performing the operation. The
requirements are only that the operands have the same limb radix as the
context (elements of contexts with different limb exponents `e` are
converted with ``gr_set_other``) and that they are elements of the ring
the context represents:

* midpoints and floating-point numbers may have any number of digits
  (the context only rounds *results* to its precision);
* radii may have any radius precision: the arithmetic on radii rounds
  operands to the radius precision of the context (in the direction
  which preserves the bound) and normalizes its results to it, so
  changing the radius precision of a context is always safe;
* the exponent limits apply to results only, not to operands (an
  operation can take an operand with an exponent out of bounds and
  return a result within the bounds);
* infinities and NaN are elements of the floating-point rings created
  with the ``DECIMAL_ALLOW_INF`` and ``DECIMAL_ALLOW_NAN`` flags, and of
  nothing else. Passing an infinity to a floating-point
  context without ``DECIMAL_ALLOW_INF`` returns ``GR_DOMAIN`` and passing
  a NaN to one without ``DECIMAL_ALLOW_NAN`` returns ``GR_UNABLE``, even
  where the result would be finite (`1/\infty`, `\infty^0`, comparisons);
  predicates returning :type:`truth_t` (which have no error channel)
  treat NaN as unknown and infinities as ordinary values. Balls
  represent real (or complex) numbers only: their midpoints are always
  finite, and an operation whose result would be infinite or undefined
  (`\log 0`, `\Gamma(0)`, `1/0`, a midpoint overflowing the exponent
  limits) returns ``GR_DOMAIN`` or ``GR_UNABLE``. A ball with an infinite
  radius is the whole real line `(-\infty, \infty)`, a perfectly good
  enclosure of a real number, and is returned whenever nothing better
  is known (for example when an intermediate ``arb`` computation
  overflows).

Functions on radii (:type:`decmag_t`) take a context argument but are
not ``gr`` methods; all their names carry a leading underscore. The same
convention is used for predicates and comparisons which return an
``int`` or ``slong`` result directly instead of following the ``gr``
conventions (a :type:`truth_t`, or a status code with the result in an
output argument): ``_decfloat_cmp`` versus ``decfloat_cmp``, or
``_decball_contains``, ``_decball_is_positive``, ``_decfloat_is_finite``.

.. function:: slong decimal_ctx_get_prec(gr_ctx_t ctx)
              int decimal_ctx_get_rnd(gr_ctx_t ctx)
              int decimal_ctx_get_rnd_im(gr_ctx_t ctx)
              slong decimal_ctx_get_rad_prec(gr_ctx_t ctx)
              int decimal_ctx_get_flags(gr_ctx_t ctx)
              slong decimal_ctx_get_limb_digits(gr_ctx_t ctx)
              void decimal_ctx_get_exp_limits(slong * emin, slong * emax, gr_ctx_t ctx)

    Accessors for the corresponding settings of a context; the limb
    digits are the number `e` such that the internal limb radix is
    `10^e`.

Rounding modes
................................................................................

.. macro:: DECIMAL_RND_DOWN
           DECIMAL_RND_UP
           DECIMAL_RND_FLOOR
           DECIMAL_RND_CEIL
           DECIMAL_RND_NEAR
           DECIMAL_RND_NEAR_AWAY
           DECIMAL_RND_NEAR_ZERO

    Rounding toward zero, away from zero, toward `-\infty`, toward
    `+\infty`, to nearest with ties to even, to nearest with ties away
    from zero, and to nearest with ties toward zero. All operations
    are correctly rounded in the selected mode.

.. macro:: DECIMAL_PREC_EXACT

    Special precision value indicating that no rounding should be
    performed.

Flags
................................................................................

.. macro:: DECIMAL_ALLOW_INF
           DECIMAL_ALLOW_NAN

    Admit infinities, respectively NaN, as floating-point numbers
    (the flags have no effect in ball contexts).
    Without these flags, operations that would produce such values
    (overflow when exponent limits are set, `1/0`, `\infty - \infty`)
    return ``GR_DOMAIN`` or ``GR_UNABLE``. TODO: a flag making a ball
    context represent intervals of extended real numbers, with infinite
    endpoints and an indeterminate value, as ``arb`` does.

.. macro:: DECIMAL_ALLOW_UNDERFLOW

    Flush results whose exponent falls below the lower exponent limit
    to zero. Without this flag such results return ``GR_UNABLE``.
    In the ball ring, a flushed midpoint is accompanied by an increase
    of the radius by `10^{e_{\min}}` so that the ball remains valid.

.. macro:: DECIMAL_SLOPPY_RADIUS

    In the ball ring, bound the rounding error of each midpoint operation
    by one ulp (half an ulp in the nearest rounding modes) instead of
    computing a tight upper bound (the exact error rounded up to the
    radius precision), which is the default. The tight bound costs
    essentially nothing extra, so this flag mainly exists for testing.

.. macro:: DECIMAL_WRITE_SCIENTIFIC

    Print numbers in scientific notation (except for numbers with
    scientific exponent zero, such as ``1`` and ``-2.5``, so that
    generators of polynomial rings print as expected).

Floating-point numbers
--------------------------------------------------------------------------------

.. type:: decfloat_struct
          decfloat_t

    A decimal floating-point number. The mantissa is the
    ``m`` field (a :type:`radix_integer_struct`) and the limb exponent is
    the ``exp`` field (an :type:`fmpz`).

.. macro:: DECFLOAT_IS_ZERO(x)
           DECFLOAT_IS_POS_INF(x)
           DECFLOAT_IS_NEG_INF(x)
           DECFLOAT_IS_INF(x)
           DECFLOAT_IS_NAN(x)
           DECFLOAT_IS_SPECIAL(x)
           DECFLOAT_IS_FINITE(x)

    Tests for special values.

Memory management and assignment
................................................................................

.. function:: void decfloat_init(decfloat_t res, gr_ctx_t ctx)
              void decfloat_clear(decfloat_t res, gr_ctx_t ctx)
              void decfloat_swap(decfloat_t x, decfloat_t y, gr_ctx_t ctx)
              void decfloat_set_shallow(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)

.. function:: int decfloat_zero(decfloat_t res, gr_ctx_t ctx)
              int decfloat_one(decfloat_t res, gr_ctx_t ctx)
              int decfloat_neg_one(decfloat_t res, gr_ctx_t ctx)
              int decfloat_pos_inf(decfloat_t res, gr_ctx_t ctx)
              int decfloat_neg_inf(decfloat_t res, gr_ctx_t ctx)
              int decfloat_nan(decfloat_t res, gr_ctx_t ctx)

    Constants. The special values return ``GR_DOMAIN`` (infinities)
    or ``GR_UNABLE`` (NaN) if not admitted by the context.

.. function:: int decfloat_set(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_set_si(decfloat_t res, slong x, gr_ctx_t ctx)
              int decfloat_set_ui(decfloat_t res, ulong x, gr_ctx_t ctx)
              int decfloat_set_fmpz(decfloat_t res, const fmpz_t x, gr_ctx_t ctx)
              int decfloat_set_fmpq(decfloat_t res, const fmpq_t x, gr_ctx_t ctx)
              int decfloat_set_d(decfloat_t res, double x, gr_ctx_t ctx)
              int decfloat_set_str(decfloat_t res, const char * s, gr_ctx_t ctx)
              int decfloat_set_arf(decfloat_t res, const arf_t x, gr_ctx_t ctx)
              int decfloat_set_fmpz_10exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
              int decfloat_set_fmpz_10exp_si(decfloat_t res, const fmpz_t m, slong e, gr_ctx_t ctx)
              int decfloat_set_si_10exp_si(decfloat_t res, slong m, slong e, gr_ctx_t ctx)
              int decfloat_set_fmpz_2exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
              int decfloat_set_other(decfloat_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)

    Sets *res* to the given value rounded to the context precision
    in the context rounding mode. Strings may be plain decimal literals
    (``"-12.5e-3"``, ``"inf"``, ``"nan"``) or arithmetic expressions
    (``"1/3"``), the latter being evaluated with the generic string
    parser. Binary values `m 2^k` are converted with exact integer
    arithmetic on a scaled truncation of the value (the digits needed
    plus guard digits, and a sticky bit) when the sizes are moderate,
    and otherwise (huge `|k|`, or a huge mantissa) by scaling with a
    power of ten in ``arb`` arithmetic combined with Ziv's strategy,
    after checking whether the value is an exact decimal number of at
    most *prec* digits (which Ziv's strategy could not decide); this is
    correctly rounded in all cases, with a cost essentially linear in
    the size of `k` for huge exponents.
    ``set_other`` accepts integers, rationals, ``arf``, ``arb`` and
    ``acb`` values (through their midpoints, with ``GR_UNABLE`` for an
    inexact ball when the target precision is exact), elements of all
    the decimal rings (balls likewise through their midpoints; complex
    values must have a zero imaginary part or ``GR_DOMAIN`` is returned),
    algebraic numbers (``qqbar``, see :func:`decfloat_set_qqbar`),
    and elements of other rings through a conversion to ``acb`` at a
    somewhat higher precision, which is not correctly rounded.

.. function:: int decfloat_set_round_qqbar(decfloat_t res, const qqbar_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_qqbar(decfloat_t res, const qqbar_t x, gr_ctx_t ctx)
              int decball_set_qqbar(decball_t res, const qqbar_t x, gr_ctx_t ctx)
              int deccfloat_set_qqbar(deccfloat_t res, const qqbar_t x, gr_ctx_t ctx)
              int deccball_set_qqbar(deccball_t res, const qqbar_t x, gr_ctx_t ctx)

    Converts the algebraic number *x*; these are also used by ``set_other``
    for elements of the ``qqbar`` rings. The prototypes are available when
    ``qqbar.h`` is included (before or after ``decimal.h``). The real
    and imaginary parts are converted separately: floating-point results
    are correctly rounded (the imaginary part of a ``deccfloat`` in the
    imaginary rounding mode), and a ball part is exact (has zero radius)
    whenever the part is a decimal number of at most the context precision
    that is within the exponent range; for example, `\sqrt{2} + 0.0023 i`
    gives an exact imaginary part at a precision of two or more digits.
    The real conversions return ``GR_DOMAIN`` for a nonreal *x*.

    The conversion first rounds a numerical enclosure obtained by refining
    the enclosure of *x*. Only when this cannot decide the result
    (a floating-point part lying on a rounding boundary, or a ball part
    whose enclosure contains a decimal number of at most the precision)
    is it checked exactly whether the part equals that decimal number:
    for the real part by subtracting the rational number, for the
    imaginary part by subtracting its product with `i` (a resultant
    computation), and then testing whether the real or imaginary part
    of the difference vanishes. The real part of a real algebraic
    number of degree greater than one is irrational and needs no check.
    With an exact target precision, a part must be an exact decimal
    number (computed exactly if *x* is not rational), otherwise
    ``GR_UNABLE`` is returned.

.. function:: int decfloat_set_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_round_si(decfloat_t res, slong x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_round_ui(decfloat_t res, ulong x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_round_fmpz(decfloat_t res, const fmpz_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_round_fmpq(decfloat_t res, const fmpq_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_round_d(decfloat_t res, double x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_round_str(decfloat_t res, const char * s, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_round_fmpz_10exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_round_fmpz_2exp_fmpz(decfloat_t res, const fmpz_t m, const fmpz_t e, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_set_round_arf(decfloat_t res, const arf_t x, slong prec, int rnd, gr_ctx_t ctx)

    Versions with an explicit precision and rounding mode.
    Passing ``DECIMAL_PREC_EXACT`` sets the exact value.

.. function:: int decfloat_set_round_info(decfloat_t res, const decfloat_t x, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)

    Rounds *x*, additionally reporting whether the result is inexact
    in *info* and, if *err* is not *NULL*, an upper bound for the absolute
    rounding error in *err*. Both pointers are optional.

.. function:: int decfloat_set_round_fmpq_reference(decfloat_t res, const fmpq_t x, slong prec, int rnd, gr_ctx_t ctx)

    Reference implementation of correctly rounded conversion using
    plain :type:`fmpz` arithmetic, used as an oracle by the test code.

Conversions
................................................................................

.. function:: int decfloat_get_fmpz(fmpz_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_get_fmpq(fmpq_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_get_si(slong * res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_get_ui(ulong * res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_get_d(double * res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_get_fmpz_10exp_fmpz(fmpz_t m, fmpz_t e, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_get_arf(arf_t res, const decfloat_t x, slong prec_bits, int rnd, gr_ctx_t ctx)
              int decfloat_get_arb(arb_t res, const decfloat_t x, slong prec_bits, gr_ctx_t ctx)
              int decfloat_get_fmpz_fixed_si(fmpz_t res, const decfloat_t x, slong e, int rnd, gr_ctx_t ctx)

    Conversions out. The integer conversions return ``GR_DOMAIN`` for
    non-integers and ``GR_UNABLE`` for integers with more than
    ``DECIMAL_CONV_DIGITS_LIMIT`` digits. The decomposition
    `x = m \cdot 10^e` is returned with `m` not divisible by 10.
    Doubles are correctly rounded to nearest; values outside
    the normal range of doubles give ``GR_UNABLE``.
    Conversion to ``arf`` is correctly rounded to *prec_bits* bits
    in the ``arf`` rounding mode *rnd* (``ARF_PREC_EXACT`` is allowed
    and succeeds when the value is a dyadic rational), and conversion
    to ``arb`` gives an enclosure with *prec_bits* bits. The latter is
    exact for a dyadic value whose exact conversion is no more expensive
    than the target precision; otherwise only the leading
    `O(\mathrm{prec\_bits})` digits of the mantissa are used (with an
    error bound for the rest), so that the cost is independent of the
    length of the mantissa and of the size of the exponent.
    The last function computes `x / 10^e` rounded to an integer.

.. function:: ulong decfloat_get_digit_si(const decfloat_t x, slong k, gr_ctx_t ctx)
              ulong decfloat_get_limb_si(const decfloat_t x, slong k, gr_ctx_t ctx)
              int decfloat_set_digit_si(decfloat_t res, const decfloat_t x, slong k, ulong c, gr_ctx_t ctx)
              int decfloat_set_limb_si(decfloat_t res, const decfloat_t x, slong k, ulong c, gr_ctx_t ctx)

    The digit `c \in \{0, \ldots, 9\}` of `|x|` at the absolute position
    `10^k`, respectively the limb `0 \le c < B` at position `B^k` (zero
    outside the mantissa and for special values), and *x* with that digit
    or limb replaced by *c*. Setting is exact (no rounding to the context
    precision, although the exponent limits are applied); the sign of *x*
    is kept, and setting a digit of zero gives a positive number.
    ``GR_DOMAIN`` is returned for an out-of-range *c* or a special value
    of *x*, and ``GR_UNABLE`` if the result would have more than
    ``DECIMAL_CONV_DIGITS_LIMIT`` digits.

.. function:: int decfloat_get_radix_integer(radix_integer_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_get_radix_integer_Bexp_fmpz(radix_integer_t m, fmpz_t e, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_set_radix_integer(decfloat_t res, const radix_integer_t m, gr_ctx_t ctx)
              int decfloat_set_radix_integer_Bexp_fmpz(decfloat_t res, const radix_integer_t m, const fmpz_t e, gr_ctx_t ctx)

    Conversions between *x* and integers in the limb radix of the
    context (``DECIMAL_CTX_RADIX(ctx)``), in time linear in the number of
    limbs. The first function requires an integer *x* (``GR_DOMAIN``
    otherwise); the second gives the exact representation `x = m B^e`
    with *m* the canonical mantissa (no zero limbs at either end). The
    set functions round `m` or `m B^e` to the context precision.

.. function:: char * decfloat_get_str(const decfloat_t x, gr_ctx_t ctx)
              char * decfloat_get_str_sci(const decfloat_t x, gr_ctx_t ctx)
              int decfloat_write(gr_stream_t out, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_write_sci(gr_stream_t out, const decfloat_t x, gr_ctx_t ctx)

    Prints *x* with all its significant digits, either in positional
    notation (when `-6 \le E \le 20` and the context does not select
    scientific notation, or when `E = 0`) or in scientific notation
    (``1.2345e-30``).

Properties and comparisons
................................................................................

.. function:: truth_t decfloat_is_zero(const decfloat_t x, gr_ctx_t ctx)
              truth_t decfloat_is_one(const decfloat_t x, gr_ctx_t ctx)
              truth_t decfloat_is_neg_one(const decfloat_t x, gr_ctx_t ctx)
              truth_t decfloat_equal(const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              truth_t decfloat_is_integer(const decfloat_t x, gr_ctx_t ctx)
              int decfloat_cmp(int * res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int decfloat_cmpabs(int * res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int decfloat_sgn(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)

    Predicates and comparisons following the ``gr`` conventions:
    predicates return ``T_UNKNOWN`` for NaN, comparisons return
    ``GR_UNABLE`` for NaN (and ``GR_DOMAIN`` for infinities when the
    context does not allow them). A number is an integer iff its limb
    exponent is nonnegative; infinities are not integers.

.. function:: int _decfloat_is_special(const decfloat_t x)
              int _decfloat_is_finite(const decfloat_t x)
              int _decfloat_is_nan(const decfloat_t x)
              int _decfloat_is_inf(const decfloat_t x)
              int _decfloat_is_pos_inf(const decfloat_t x)
              int _decfloat_is_neg_inf(const decfloat_t x)
              int _decfloat_is_int(const decfloat_t x, gr_ctx_t ctx)
              int _decfloat_sgn(const decfloat_t x, gr_ctx_t ctx)
              int _decfloat_cmp(const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int _decfloat_cmpabs(const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int _decfloat_cmp_si(const decfloat_t x, slong c, gr_ctx_t ctx)
              int _decfloat_cmp_ui(const decfloat_t x, ulong c, gr_ctx_t ctx)
              int _decfloat_cmpabs_ui(const decfloat_t x, ulong c, gr_ctx_t ctx)

    Raw versions returning the result directly. Comparisons involving
    NaN return 0.

.. function:: slong decfloat_digits(const decfloat_t x, gr_ctx_t ctx)
              slong decfloat_limbs(const decfloat_t x, gr_ctx_t ctx)
              slong _decfloat_mant_digits(const decfloat_t x, gr_ctx_t ctx)
              slong _decfloat_val10_clamped(const decfloat_t x, gr_ctx_t ctx)
              void decfloat_get_sci_exp(fmpz_t E, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_get_sci_exp_si(slong * E, const decfloat_t x, gr_ctx_t ctx)

    Writing `|x| = M \cdot 10^v` with `M` an integer not divisible by 10:
    the number of digits of `M` (the smallest precision at which *x* is
    representable; decimal numbers do not track significance, so this is
    a property of the value), the number of limbs of the mantissa,
    the number of digits of the limb mantissa including the trailing zeros
    of its lowest limb, the exponent `v` (clamped to
    `\pm` ``DECFLOAT_EXP_CLAMP``), and the scientific exponent `E` with
    `10^E \le |x| < 10^{E+1}` (the last function returns 0 if `E` does
    not fit in a word). The digit counts are zero for special values; the
    value must be nonzero and finite for the exponent functions.

Arithmetic
................................................................................

.. function:: int decfloat_neg(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_abs(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_add(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int decfloat_sub(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int decfloat_mul(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int decfloat_sqr(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_div(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int decfloat_inv(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_sqrt(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_rsqrt(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_add_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx)
              int decfloat_add_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx)
              int decfloat_add_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx)
              int decfloat_sub_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx)
              int decfloat_sub_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx)
              int decfloat_sub_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx)
              int decfloat_mul_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx)
              int decfloat_mul_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx)
              int decfloat_mul_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx)
              int decfloat_div_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx)
              int decfloat_div_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx)
              int decfloat_div_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx)
              int decfloat_mul_two(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_mul_10exp_si(decfloat_t res, const decfloat_t x, slong e, gr_ctx_t ctx)
              int decfloat_mul_10exp_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t e, gr_ctx_t ctx)
              int decfloat_mul_2exp_si(decfloat_t res, const decfloat_t x, slong e, gr_ctx_t ctx)
              int decfloat_mul_2exp_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t e, gr_ctx_t ctx)
              int decfloat_floor(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_ceil(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_trunc(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_nint(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_pow_ui(decfloat_t res, const decfloat_t x, ulong y, gr_ctx_t ctx)
              int decfloat_pow_si(decfloat_t res, const decfloat_t x, slong y, gr_ctx_t ctx)
              int decfloat_pow_fmpz(decfloat_t res, const decfloat_t x, const fmpz_t y, gr_ctx_t ctx)

    Arithmetic operations, correctly rounded to the context precision
    in the context rounding mode. Integer powers are computed exactly
    and then rounded once when the exact result has at most 100000
    digits (always for powers of ten), and correctly rounded via
    ``arb`` otherwise. Multiplication by `10^e` is exact
    for any *e* (a digit shift is only needed when *e* is not a multiple
    of the limb exponent). The integer rounding functions first
    round to an integer exactly and then round the integer to the context
    precision. Division by zero returns ``GR_DOMAIN`` unless infinities
    are admitted; square roots of negative numbers return ``GR_DOMAIN``.

    Multiplication by `2^e` is exact followed by one rounding for
    `|e| < 2^{16}`, and otherwise correctly rounded through the binary
    conversion (so that huge `|e|` costs `O(\log |e|)` arithmetic).

    In the exact ring, division and square root succeed exactly when
    the result is a decimal fraction (denominator of the form `2^i 5^j`).
    Exact results with more than ``DECIMAL_CONV_DIGITS_LIMIT`` digits
    are not computed: products, integer powers (estimated before
    computing) and multiplications by `2^e` then return ``GR_UNABLE``.

.. function:: int decfloat_neg_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_abs_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_add_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_sub_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_mul_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_sqr_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_div_round(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_inv_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_sqrt_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_rsqrt_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_mul_10exp_si_round(decfloat_t res, const decfloat_t x, slong e, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_mul_10exp_fmpz_round(decfloat_t res, const decfloat_t x, const fmpz_t e, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_mul_2exp_fmpz_round(decfloat_t res, const decfloat_t x, const fmpz_t e, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_floor_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_ceil_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_trunc_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)
              int decfloat_nint_round(decfloat_t res, const decfloat_t x, slong prec, int rnd, gr_ctx_t ctx)

    Versions with an explicit precision and rounding mode, allowing
    variable-precision computations without modifying the context.

.. function:: int _decfloat_add(decfloat_t res, const decfloat_t x, const decfloat_t y, int negate_y, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
              int _decfloat_mul(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
              int _decfloat_div(decfloat_t res, const decfloat_t x, const decfloat_t y, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
              int _decfloat_sqrt(decfloat_t res, const decfloat_t x, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
              int _decfloat_mul_10exp(decfloat_t res, const decfloat_t x, const fmpz_t e, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
              int _decfloat_round_to_int(decfloat_t res, const decfloat_t x, int int_rnd, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)

    Underlying implementations which additionally report rounding
    information (see :func:`decfloat_set_round_info`). These are used
    by the ball arithmetic.

    Addition aligns the operands by a limb offset and adds them
    exactly, unless the smaller operand lies entirely below the
    rounding horizon of the larger one, in which case the sum is
    computed as the larger operand plus a sticky tail (this makes
    additions with astronomically different exponents cost `O(\mathrm{prec})`).
    Multiplication first computes a truncated high product with a
    known error bound and rounds both ends of the resulting interval;
    if the two roundings agree the result is correctly rounded, and
    otherwise (with probability roughly `n / B`) the full product is
    computed. Division computes a truncated quotient with at least
    `\mathrm{prec} + 1` digits by Euclidean division and uses the remainder
    as a sticky flag; square roots are computed likewise with
    :func:`radix_sqrtrem`.

.. function:: int _decfloat_set_round_limbs(decfloat_t res, nn_srcptr d, slong n, int negative, const fmpz_t exp, int eps, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)
              slong _decimal_round_mantissa(nn_ptr d, slong n, int negative, int eps, slong prec, int rnd, slong * newn, decimal_rounding_info * info, decmag_ptr err, const fmpz_t exp, gr_ctx_t ctx)

    The rounding core: sets *res* to `(-1)^{\mathrm{negative}} (d, n) B^{\mathrm{exp}}`
    (plus an unrepresented positive tail smaller than one unit
    of the lowest limb if *eps* is set) rounded to *prec* digits.
    Limbs of *d* may be zero at either end. The buffer must have room
    for `n + 1` limbs (a carry may propagate out of the top) and is
    modified. The second function performs the in-place rounding,
    returning the number of low limbs to discard.

.. type:: decimal_rounding_info

    Structure with fields ``inexact``, ``increased`` (whether the
    magnitude was rounded up), ``underflow`` and ``overflow``.

Correct rounding
................................................................................

The following operations on ``decfloat`` are guaranteed to be correctly
rounded (to the target precision in the target rounding mode, with the
result also being exact whenever the mathematical result is
representable):

* conversions from integers, rationals, ``arf``, ``arb`` midpoints,
  doubles and strings (plain literals), and to ``arf``;
* ``add``, ``sub``, ``mul``, ``sqr``, ``div``, ``inv``, ``sqrt``,
  ``rsqrt``, ``mul_2exp``, ``mul_10exp`` (exact), ``floor``, ``ceil``,
  ``trunc``, ``nint`` (exact), ``neg``, ``abs``;
* ``pow`` (all cases), ``pow_ui``, ``pow_si``, ``pow_fmpz``;
* every function listed under *Elementary and special functions* below,
  subject to the caveat that Ziv's strategy is used for values that are
  neither exact nor covered by the special-argument handling, and
  returns ``GR_UNABLE`` in the (practically nonexistent) event that the
  rounding has not been decided at 64 times the working precision.

Operations reached through the generic fallbacks of the ``gr`` interface
are *not* correctly rounded, since they are composed of several rounded
operations. Examples include: ``gr_vec_dot`` and hence all matrix and
polynomial arithmetic (each product and sum is rounded separately),
``gr_pow_fmpq``, ``gr_hypot``, ``gr_fac``, ``gr_bin``, ``gr_sqrt_ui``
and other functions built
generically from the basic operations, power series functions, and
``gr_set_str`` applied to an arithmetic expression such as ``"1/3 + 2/7"``
(each operation in the expression is rounded).

Elementary and special functions
................................................................................

.. function:: int decfloat_pi(decfloat_t res, gr_ctx_t ctx)

    Constants, correctly rounded.

.. function:: int decfloat_exp(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_expm1(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_log(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_log1p(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_sin(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_cos(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_sin_cos(decfloat_t res1, decfloat_t res2, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_tan(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_asin(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_acos(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_atan(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_atan2(decfloat_t res, const decfloat_t y, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_sinh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_cosh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_sinh_cosh(decfloat_t res1, decfloat_t res2, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_tanh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_asinh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_acosh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_atanh(decfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int decfloat_pow(decfloat_t res, const decfloat_t x, const decfloat_t y, gr_ctx_t ctx)

    The most common elementary functions, correctly rounded. Versions for
    the other types (``decball_exp``, ``deccfloat_exp``, ``deccball_exp``,
    ...) are documented with those types.

Other elementary and special functions are only available through the
generic-ring interface: with a decimal context, :func:`gr_gamma`,
:func:`gr_bessel_j` etc. dispatch to implementations handling all four
decimal types, which are correctly rounded (componentwise for complex
floats) in the floating-point contexts and give enclosures in the ball
contexts. The arguments and flags have the same meaning as for the
corresponding ``gr`` methods. The following methods are implemented
(for all four types unless marked with *complex*):

* Constants: ``pi``, ``euler``, ``catalan``, ``khinchin``, ``glaisher``.
* Exponentials and logarithms: ``exp``, ``expm1``, ``exp2``, ``exp10``,
  ``exp_pi_i`` (complex), ``log``, ``log1p``, ``log2``, ``log10``,
  ``log_pi_i`` (complex), ``pow``.
* Trigonometric functions: ``sin``, ``cos``, ``sin_cos``, ``tan``,
  ``cot``, ``sec``, ``csc``, ``sinc``, ``sin_pi``, ``cos_pi``,
  ``sin_cos_pi``, ``tan_pi``, ``cot_pi``, ``sec_pi``, ``csc_pi``,
  ``sinc_pi``, ``asin``, ``acos``, ``atan``, ``atan2`` (real), ``acot``,
  ``asec``, ``acsc``, ``asin_pi``, ``acos_pi``, ``atan_pi``,
  ``acot_pi``, ``asec_pi``, ``acsc_pi``.
* Hyperbolic functions: ``sinh``, ``cosh``, ``sinh_cosh``, ``tanh``,
  ``coth``, ``sech``, ``csch``, ``asinh``, ``acosh``, ``atanh``,
  ``acoth``, ``asech``, ``acsch``.
* Lambert W function: ``lambertw``, ``lambertw_fmpz``.
* Gamma and related functions: ``gamma``, ``gamma_fmpz``,
  ``gamma_fmpq``, ``rgamma``, ``lgamma``, ``digamma``, ``polygamma``
  (complex), ``fac_ui``, ``fac_fmpz``, ``rising``, ``rising_ui``,
  ``barnes_g``, ``log_barnes_g``.
* Zeta and related functions: ``zeta``, ``hurwitz_zeta``,
  ``dirichlet_eta`` (complex), ``riemann_xi`` (complex), ``polylog``,
  ``dilog``, ``lerch_phi`` (complex).
* Error functions and integrals: ``erf``, ``erfc``, ``erfi``,
  ``erfinv``, ``erfcinv``, ``fresnel``, ``fresnel_s``, ``fresnel_c``,
  ``gamma_upper``, ``gamma_lower``, ``beta_lower``, ``exp_integral``,
  ``exp_integral_ei``, ``sin_integral``, ``cos_integral``,
  ``sinh_integral``, ``cosh_integral``, ``log_integral``.
* Bessel and Airy functions: ``bessel_j``, ``bessel_y``, ``bessel_i``,
  ``bessel_k``, ``bessel_i_scaled``, ``bessel_k_scaled``, ``airy``,
  ``airy_ai``, ``airy_bi``, ``airy_ai_prime``, ``airy_bi_prime``,
  ``coulomb_f``, ``coulomb_g``.
* Orthogonal polynomials: ``chebyshev_t``, ``chebyshev_u``,
  ``jacobi_p``, ``gegenbauer_c``, ``laguerre_l``, ``hermite_h``,
  ``legendre_p``, ``legendre_q``.
* Hypergeometric functions: ``hypgeom_0f1``, ``hypgeom_1f1``,
  ``hypgeom_u``, ``hypgeom_2f1``.
* Arithmetic-geometric mean and elliptic functions: ``agm``, ``agm1``,
  ``elliptic_k``, ``elliptic_e``, ``elliptic_pi``, ``elliptic_f``,
  ``elliptic_e_inc`` (complex).
* Modular forms: ``modular_j``, ``modular_lambda``, ``modular_delta``,
  ``dedekind_eta`` (complex).

.. function:: int decfloat_via_arb(decfloat_t res, const decfloat_t x, int (*func)(arb_t, const arb_t, slong), gr_ctx_t ctx)
              int _decfloat_arb_gr(decfloat_ptr * res, slong nres, decfloat_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx)

    The machinery behind the real functions above: evaluates a user-supplied
    ``arb`` function (which must return a ``gr`` status code), or the
    ``gr`` method *method* of the ``arb`` ring applied to *nargs*
    arguments (with an integer flag if *has_flag* is set) producing
    *nres* results, at increasing precision (Ziv's strategy) until the
    rounding of every result is decided. ``GR_UNABLE`` is returned if
    this does not happen within 64 times the working precision; for the
    elementary functions (exponentials, logarithms, powers, trigonometric
    and hyperbolic functions and their inverses), whose cost at high
    precision is quasi-linear, the limit is increased by twice the
    number of digits of the arguments, so that arguments with more digits
    than the target precision (which matter near cancellation points,
    e.g. `\log_2(1 + 10^{-10^6})`) are resolved in time proportional to
    their length. Points outside the real domain of the inverse
    functions and logarithms (`\mathrm{acosh}(1/2)`, `\log(-1)`) are
    detected before any evaluation and give ``GR_DOMAIN``.

    Ziv's strategy cannot decide the rounding when the exact result lies
    extremely close to a representable number. This happens
    systematically in three situations, which the functions above
    therefore treat separately:

    * Exact values: `\Gamma(n) = (n-1)!`, `1/\Gamma(n)`, `\zeta(-n)`,
      `\zeta(0) = -1/2`, `\log_2(2^n)`, `\log_{10}(10^n)`,
      `\sin(\pi k/4)`, `\tan(\pi k/4)`, `\cot(\pi k/4)`, `\csc(\pi k/4)`,
      `\sec(\pi k/4)`, `\mathrm{asin}(1/2)/\pi = 1/6` and the other
      rational values of the inverse trigonometric functions divided by
      `\pi`, `\mathrm{acos}(1) = \mathrm{acosh}(1) = 0`,
      `\log \Gamma(1) = \log \Gamma(2) = 0`, `G(n)`, `\log G(n)` for
      `n = 1, 2, 3`, `f(0)` for functions whose ``arb`` implementation may
      not return an exact zero or one, powers `x^n`, `x^{p/q}` and
      `2^x`, `10^x` with exact results, and `\mathrm{atan2}` on the axes.
      (``arb`` does not detect most of these.) These are detected from
      the sign, the scientific exponent, the exponent `v` above and
      individual digits of the argument, in time independent of its
      length and exponent (for instance `\sin(\pi x)` depends on three
      digits, and all trivial zeros `\zeta(-2k)` and poles of `\Gamma`
      are found for any `k`); the argument is only converted to an
      integer or rational number when it is small or a genuine candidate
      (`\log_2(x)` first checks the digit count, the last digit and a
      residue modulo a word-size prime).
    * Tiny arguments: when the leading Taylor term `S` (`x`, `1`, `\pm
      1/x`, `x - 1` for `\log`, `1/(x-1)` for `\zeta` near 1, `-1/2`
      for `\zeta` near 0, `1/2` or `\operatorname{sgn}(x)/2` for
      `\mathrm{acos}(x)/\pi`, `\mathrm{acot}(x)/\pi` near 0 and
      `\mathrm{atan}(x)/\pi`, `\mathrm{asec}(x)/\pi` near infinity,
      `x` for the Barnes G-function) is an exact
      decimal and the remaining tail, of known sign, lies more than
      ``prec`` `+ 4e` digits below it, the result is `S` rounded with the
      tail acting as a sticky bit. For example `\sin(10^{-10^9})`,
      `\Gamma(10^{-10^9}) = 10^{10^9} - \gamma + \ldots`,
      `\zeta(1 + 10^{-10^9})` and `\cot(10^{-10^9})` are all correctly
      rounded at negligible cost.
    * Large arguments: `\tanh`, `\coth`, `\mathrm{erf}`, `\mathrm{erfc}`,
      `\mathrm{expm1}` and `\zeta` converge exponentially fast to
      `\pm 1`, `2`, `-1` or `1`, and `\mathrm{acot}`, `\mathrm{acoth}`,
      `\mathrm{acsc}`, `\mathrm{acsch}` behave like `1/x`; these are
      handled in the same way.

    The remaining special functions (Bessel functions, hypergeometric
    functions, ...) rely on Ziv's strategy alone and can in principle
    return ``GR_UNABLE`` for arguments at which the exact value is a
    representable number not detected by ``arb``.

    Some arguments are rejected with ``GR_UNABLE`` before any evaluation
    (for balls as well), since the ``arb``/``acb`` algorithms would need
    time or memory far beyond what the result is worth: for the
    non-elementary functions, an argument with a part of magnitude at
    least `10^{10^6}` (and for the Barnes G-function also at most
    `10^{-10^6}`, and at least `10^{50}`); for the polylogarithm, the
    Hurwitz zeta function and the polygamma function, an order `s` with
    `|\operatorname{Re}(s)| \ge 10^{10}` or
    `|\operatorname{Im}(s)| \ge 10^{7}` (there is no Riemann-Siegel type
    formula); and for the zeta-type functions (`\zeta`,
    Hurwitz `\zeta`, `\eta`, `\beta`, `\operatorname{Li}_s`, `\xi`)
    in the complex contexts, `|\operatorname{Re}(s)| \ge 10^{30}`,
    `|\operatorname{Im}(s)| \ge 10^{15}` (the number of terms grows with
    `|\operatorname{Im}(s)|`) and, for `\xi`, an imaginary part whose
    cancellation (`\xi(s)` decays like `e^{-\pi |\operatorname{Im}(s)|/4}`)
    exceeds the maximal working precision.

Random generation
................................................................................

.. function:: int decfloat_randtest(decfloat_t res, flint_rand_t state, gr_ctx_t ctx)
              int decfloat_randtest_special(decfloat_t res, flint_rand_t state, gr_ctx_t ctx)

    Generates a random element rounded to the context precision,
    respectively a random element which may have more significant
    digits than the context precision.

Radius bounds
--------------------------------------------------------------------------------

.. type:: decmag_struct
          decmag_t

    A nonnegative number `m \cdot 10^{k}` with the mantissa `m` stored
    in a single word, `1 \le m < 10^{R}` where ``R = DECMAG_MAX_PREC``
    (9 on a 64-bit machine, 4 on a 32-bit machine), and an
    :type:`fmpz` exponent `k` measured in digits, or the special values
    zero (`m = 0`) and `+\infty` (`m = ` ``UWORD_MAX``). The radius
    precision *rad_prec* is a context setting between 1 and `R` digits.
    The results of all operations are *normalized*, with exactly
    *rad_prec* digits (`10^{\mathrm{rad\_prec}-1} \le m < 10^{\mathrm{rad\_prec}}`,
    tested by the macro ``DECMAG_IS_NORMALIZED(x, ctx)``), so that
    products of mantissas fit in a word and all operations reduce to
    a few word operations plus a small exponent adjustment.
    Operands may have any number of digits (up to `R`): they are first
    rounded to *rad_prec* digits in the direction preserving the bound
    being computed, which costs two comparisons per operand in the
    common case where nothing needs to be done. Comparisons are
    comparisons of values, independently of the number of digits.

    The functions below compute upper bounds unless their name ends in
    ``_lower``. Upper bounds are tight: the result is the exact value
    rounded up to *rad_prec* digits, except where noted.

.. function:: void _decmag_init(decmag_t res, gr_ctx_t ctx)
              void _decmag_clear(decmag_t res, gr_ctx_t ctx)
              void _decmag_swap(decmag_t x, decmag_t y, gr_ctx_t ctx)
              void _decmag_zero(decmag_t res, gr_ctx_t ctx)
              void _decmag_one(decmag_t res, gr_ctx_t ctx)
              void _decmag_inf(decmag_t res, gr_ctx_t ctx)
              void _decmag_set(decmag_t res, const decmag_t x, gr_ctx_t ctx)
              void _decmag_set_lower(decmag_t res, const decmag_t x, gr_ctx_t ctx)
              void _decmag_set_round(decmag_t res, const decmag_t x, slong prec, gr_ctx_t ctx)
              int _decmag_is_zero(const decmag_t x, gr_ctx_t ctx)
              int _decmag_is_inf(const decmag_t x, gr_ctx_t ctx)
              int _decmag_is_finite(const decmag_t x, gr_ctx_t ctx)
              int _decmag_equal(const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              int _decmag_cmp(const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              int _decmag_cmp_10exp_si(const decmag_t x, slong e, gr_ctx_t ctx)
              int _decmag_is_10exp(const decmag_t x, gr_ctx_t ctx)
              slong _decmag_digits(const decmag_t x)
              void _decmag_get_sci_exp(fmpz_t E, const decmag_t x)

    Basic operations. ``set`` rounds *x* up to the radius precision of
    the context (a plain copy for a normalized *x*), ``set_lower`` rounds
    down, and ``set_round`` rounds up to *prec* digits (clamped to
    `[1, R]`) regardless of the context setting. ``digits`` returns the
    number of digits of the mantissa (0 for zero and infinity) and
    ``get_sci_exp`` the exponent `E` with `10^E \le x < 10^{E+1}`
    (*x* finite and nonzero).

.. function:: void _decmag_set_ui(decmag_t res, ulong x, gr_ctx_t ctx)
              void _decmag_set_ui_lower(decmag_t res, ulong x, gr_ctx_t ctx)
              void _decmag_set_fmpz(decmag_t res, const fmpz_t x, gr_ctx_t ctx)
              void _decmag_set_fmpz_lower(decmag_t res, const fmpz_t x, gr_ctx_t ctx)
              int _decmag_set_d(decmag_t res, double x, gr_ctx_t ctx)
              void _decmag_set_10exp_si(decmag_t res, slong e, gr_ctx_t ctx)
              void _decmag_set_10exp_fmpz(decmag_t res, const fmpz_t e, gr_ctx_t ctx)
              void _decmag_set_ui_10exp_si(decmag_t res, ulong m, slong e, gr_ctx_t ctx)
              void _decmag_set_ui_10exp_fmpz(decmag_t res, ulong m, const fmpz_t e, gr_ctx_t ctx)
              void _decmag_set_ui_10exp_fmpz_lower(decmag_t res, ulong m, const fmpz_t e, gr_ctx_t ctx)
              void _decmag_set_uiui_10exp_fmpz(decmag_t res, ulong hi, ulong lo, ulong hi_radix, int sticky, const fmpz_t e, gr_ctx_t ctx)
              void _decmag_set_decfloat(decmag_t res, const decfloat_t x, gr_ctx_t ctx)
              void _decmag_set_decfloat_lower(decmag_t res, const decfloat_t x, gr_ctx_t ctx)
              void _decmag_set_mag(decmag_t res, const mag_t x, gr_ctx_t ctx)
              int _decmag_get_decfloat(decfloat_t res, const decmag_t x, gr_ctx_t ctx)
              void _decmag_get_fmpq(fmpq_t res, const decmag_t x, gr_ctx_t ctx)
              int _decmag_get_d(double * res, const decmag_t x, gr_ctx_t ctx)
              void _decmag_get_mag(mag_t res, const decmag_t x, gr_ctx_t ctx)

    Conversions. The absolute value is taken of signed inputs.
    ``set_uiui_10exp_fmpz`` sets *res* to an upper bound for
    `(\mathrm{hi} \cdot \mathrm{hi\_radix} + \mathrm{lo} + \varepsilon) 10^e`
    where *hi_radix* is a power of ten and `0 \le \varepsilon < 1`
    with `\varepsilon > 0` iff *sticky* is set.
    The conversion to :type:`decfloat_t` is exact. ``set_mag`` gives
    the tightest upper bound at the radius precision when the binary
    exponent is at most a few thousand (using exact integer arithmetic
    on a scaled truncation), and a bound within a relative error of
    about `10^{-12}` (rounded up to the radius precision) otherwise,
    computed from a fixed-point approximation of `\log_{10} 2` and a
    double-precision power of ten with a safety margin; this covers
    arbitrarily large exponents, in time essentially independent of
    their size.

.. function:: void _decmag_set_ulp(decmag_t res, const decfloat_t x, slong prec, gr_ctx_t ctx)

    Sets *res* to `10^{E - \mathrm{prec} + 1}` where `E` is the scientific
    exponent of *x*, i.e. the unit in the last place of *x* rounded
    to *prec* digits. Requires *x* to be nonzero and finite.

.. function:: void _decmag_add(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_add_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_sub_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_mul(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_mul_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_addmul(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_div(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_div_lower(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_inv(decmag_t res, const decmag_t x, gr_ctx_t ctx)
              void _decmag_sqrt(decmag_t res, const decmag_t x, gr_ctx_t ctx)
              void _decmag_sqrt_lower(decmag_t res, const decmag_t x, gr_ctx_t ctx)
              void _decmag_rsqrt(decmag_t res, const decmag_t x, gr_ctx_t ctx)
              void _decmag_mul_ui(decmag_t res, const decmag_t x, ulong y, gr_ctx_t ctx)
              void _decmag_div_ui(decmag_t res, const decmag_t x, ulong y, gr_ctx_t ctx)
              void _decmag_mul_10exp_si(decmag_t res, const decmag_t x, slong e, gr_ctx_t ctx)
              void _decmag_mul_10exp_fmpz(decmag_t res, const decmag_t x, const fmpz_t e, gr_ctx_t ctx)
              void _decmag_max(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_min(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)
              void _decmag_pow_ui(decmag_t res, const decmag_t x, ulong e, gr_ctx_t ctx)
              void _decmag_hypot(decmag_t res, const decmag_t x, const decmag_t y, gr_ctx_t ctx)

    Arithmetic with upper (or lower) bound rounding. ``sub_lower``
    computes a lower bound for `\max(x - y, 0)`. ``max`` rounds up and
    ``min`` rounds down to the radius precision.

.. function:: char * _decmag_get_str(const decmag_t x, gr_ctx_t ctx)
              int _decmag_write(gr_stream_t out, const decmag_t x, gr_ctx_t ctx)
              void _decmag_randtest(decmag_t res, flint_rand_t state, gr_ctx_t ctx)
              void _decmag_randtest_special(decmag_t res, flint_rand_t state, gr_ctx_t ctx)

    A radius is printed with all its digits, in positional notation when
    its scientific exponent `E` satisfies `-3 \le E \le 5` and in
    scientific notation otherwise (or always, with
    ``DECIMAL_WRITE_SCIENTIFIC``).

Balls
--------------------------------------------------------------------------------

.. type:: decball_struct
          decball_t

    A ball `[m \pm r]` consisting of a midpoint of type
    :type:`decfloat_struct` (field ``mid``) and a radius of type
    :type:`decmag_struct` (field ``rad``), representing the set
    of real numbers `\{ x : |x - m| \le r \}`. The midpoint is always
    finite; an infinite radius denotes the whole real line.
    Arithmetic operations on balls are guaranteed to produce balls
    containing the exact result of the operation applied to any points
    of the input balls, with midpoints computed to the context precision
    in the context rounding mode.

    The ball ring is presented to the generics interface as the
    (inexact) field of real numbers: predicates return ``T_UNKNOWN`` when
    the balls do not decide the answer, ``cmp`` returns ``GR_UNABLE`` for
    overlapping balls, division by a ball containing zero returns
    ``GR_UNABLE`` (``GR_DOMAIN`` for an exact zero), and any operation
    whose result would be infinite or undefined returns ``GR_DOMAIN`` or
    ``GR_UNABLE``.

.. macro:: DECBALL_MIDREF(x)
           DECBALL_RADREF(x)

    Pointers to the midpoint and radius.

.. function:: decfloat_ptr decball_midref(decball_t x)
              decmag_ptr decball_radref(decball_t x)

    Function versions of the above, for use through foreign function
    interfaces.

.. function:: void decball_init(decball_t res, gr_ctx_t ctx)
              void decball_clear(decball_t res, gr_ctx_t ctx)
              void decball_swap(decball_t x, decball_t y, gr_ctx_t ctx)
              void decball_set_shallow(decball_t res, const decball_t x, gr_ctx_t ctx)
              int _decball_is_exact(const decball_t x, gr_ctx_t ctx)
              int _decball_is_finite(const decball_t x, gr_ctx_t ctx)

.. function:: int decball_zero(decball_t res, gr_ctx_t ctx)
              int decball_one(decball_t res, gr_ctx_t ctx)
              int decball_neg_one(decball_t res, gr_ctx_t ctx)
              int decball_zero_pm_inf(decball_t res, gr_ctx_t ctx)

    Constants; ``zero_pm_inf`` is the whole real line `[0 \pm \infty]`.
    The ``gr`` methods ``pos_inf``, ``neg_inf``, ``uinf``, ``undefined``
    and ``unknown`` return ``GR_DOMAIN``.

.. function:: int decball_set(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_set_round(decball_t res, const decball_t x, slong prec, gr_ctx_t ctx)
              int decball_set_round2(decball_t res, const decball_t x, slong prec, slong rad_prec, gr_ctx_t ctx)
              int decball_set_decfloat(decball_t res, const decfloat_t x, gr_ctx_t ctx)
              int decball_set_si(decball_t res, slong x, gr_ctx_t ctx)
              int decball_set_ui(decball_t res, ulong x, gr_ctx_t ctx)
              int decball_set_fmpz(decball_t res, const fmpz_t x, gr_ctx_t ctx)
              int decball_set_fmpq(decball_t res, const fmpq_t x, gr_ctx_t ctx)
              int decball_set_d(decball_t res, double x, gr_ctx_t ctx)
              int decball_set_str(decball_t res, const char * s, gr_ctx_t ctx)
              int decball_set_fmpz_10exp_fmpz(decball_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
              int decball_set_other(decball_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int decball_set_arb(decball_t res, const arb_t x, gr_ctx_t ctx)
              int decball_get_arb(arb_t res, const decball_t x, slong prec_bits, gr_ctx_t ctx)
              int decball_get_fmpz(fmpz_t res, const decball_t x, gr_ctx_t ctx)
              int decball_get_fmpq(fmpq_t res, const decball_t x, gr_ctx_t ctx)
              int decball_get_d(double * res, const decball_t x, gr_ctx_t ctx)
              int decball_get_mid(decfloat_t res, const decball_t x, gr_ctx_t ctx)
              void decball_get_rad(decmag_t res, const decball_t x, gr_ctx_t ctx)
              void decball_get_abs_ubound(decmag_t res, const decball_t x, gr_ctx_t ctx)
              void decball_get_abs_lbound(decmag_t res, const decball_t x, gr_ctx_t ctx)

    Assignment and conversions. Values are rounded to the context
    precision with the rounding error added to the radius, and radii
    are rounded up to the radius precision of the context; ``set_round``
    rounds the midpoint to *prec* digits instead and ``set_round2``
    additionally rounds the radius up to *rad_prec* digits (both
    preserve the enclosure). An ``arb`` ball with an infinite or NaN
    midpoint is not a real number and converts with ``GR_DOMAIN``
    respectively ``GR_UNABLE``; one with an infinite radius converts
    to `[0 \pm \infty]`. Plain decimal
    literals are parsed exactly (then rounded); any other string goes
    through the generic expression parser, so ``"[mid +/- rad]"`` (the
    output format), ``"mid +/- rad"``, ``"+/- rad"`` and arithmetic
    expressions are accepted, and the same syntax works inside strings
    for polynomials, matrices and other objects over decimal balls. Note
    that a radius smaller than `10^{emin}` can only be parsed when
    underflow is allowed in the context. Since ``inf`` is not an element
    of the ring, the printed forms of balls with an infinite radius
    (``"[m +/- inf]"``, ``"[+/- inf]"``) are parsed through
    ``set_interval_mid_inf`` rather than by evaluating ``inf``. Algebraic
    numbers are converted with :func:`decball_set_qqbar`, which gives an
    exact ball for an exactly representable value. The conversions to and from
    ``arb`` are rigorous (the radius accounts for all conversion errors)
    and handle arbitrarily large exponents by scaling with a power of
    ten. The conversion from ``arb`` is also tight: the midpoint is
    correctly rounded and the exact rounding error, bounded to the
    radius precision, is added to the radius, and an exact midpoint
    that is exactly representable at the context precision gives an
    exact ball even when its exponent is huge; when scaling in ``arb``
    is used, the working precision is increased until the scaling error
    is negligible compared to the radius (unless this would require an
    unreasonable precision).
    Integer conversion succeeds for exact integer balls, returns
    ``GR_DOMAIN`` for balls containing no integer, and ``GR_UNABLE``
    otherwise.

.. function:: int decball_set_interval_mid_rad(decball_t res, const decball_t m, const decball_t r, gr_ctx_t ctx)
              int decball_set_interval_mid_inf(decball_t res, const decball_t m, gr_ctx_t ctx)

    Sets *res* to the ball *m* with an upper bound for `|r|` added to its
    radius, respectively with an infinite radius. These implement the
    ``+/-`` operator of the generic expression parser.

.. function:: int decball_mid(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_rad(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_shell(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_lower(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_upper(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_abs_lower(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_abs_upper(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_add_rad(decball_t res, const decball_t x, const decball_t r, gr_ctx_t ctx)
              int decball_set_interval(decball_t res, const decball_t lo, const decball_t hi, gr_ctx_t ctx)

    Components and bounds of a ball `x = [m \pm r]`, as balls.
    ``mid`` gives the exact ball `m` and ``rad`` the exact ball `r`
    (midpoints and radii are always representable), and
    ``shell`` gives `[0 \pm r]`. ``lower`` and ``upper`` give the endpoints
    `m - r` and `m + r` rounded to the context precision toward
    `-\infty` and `+\infty` respectively, and ``abs_lower`` and
    ``abs_upper`` give lower and upper bounds for `|x|`
    (`\max(|m| - r, 0)` rounded toward zero and `|m| + r` rounded up),
    all as exact balls, so that the results enclose the true bounds.
    An infinite radius gives ``GR_DOMAIN`` (the endpoints are not real
    numbers) except for ``abs_lower``, which gives zero, and ``rad``,
    which cannot represent an infinite radius as a ball.
    ``add_rad`` sets *res* to `[x \pm |r|]` (the same as
    ``set_interval_mid_rad``), and ``set_interval`` sets *res* to
    the smallest ball containing both *lo* and *hi* (their convex hull,
    in either order), with the midpoint rounded to the context precision.

.. function:: truth_t decball_is_zero(const decball_t x, gr_ctx_t ctx)
              truth_t decball_is_one(const decball_t x, gr_ctx_t ctx)
              truth_t decball_is_neg_one(const decball_t x, gr_ctx_t ctx)
              truth_t decball_equal(const decball_t x, const decball_t y, gr_ctx_t ctx)
              truth_t decball_is_integer(const decball_t x, gr_ctx_t ctx)
              int _decball_contains_zero(const decball_t x, gr_ctx_t ctx)
              int _decball_contains_decfloat(const decball_t x, const decfloat_t y, gr_ctx_t ctx)
              int _decball_contains_fmpq(const decball_t x, const fmpq_t y, gr_ctx_t ctx)
              int _decball_contains_fmpz(const decball_t x, const fmpz_t y, gr_ctx_t ctx)
              int _decball_contains_si(const decball_t x, slong y, gr_ctx_t ctx)
              int _decball_contains(const decball_t x, const decball_t y, gr_ctx_t ctx)
              int _decball_overlaps(const decball_t x, const decball_t y, gr_ctx_t ctx)
              int _decball_is_positive(const decball_t x, gr_ctx_t ctx)
              int _decball_is_negative(const decball_t x, gr_ctx_t ctx)
              int _decball_is_nonnegative(const decball_t x, gr_ctx_t ctx)
              int _decball_is_nonpositive(const decball_t x, gr_ctx_t ctx)
              int _decball_contains_negative(const decball_t x, gr_ctx_t ctx)
              int _decball_contains_positive(const decball_t x, gr_ctx_t ctx)
              int _decball_contains_nonnegative(const decball_t x, gr_ctx_t ctx)
              int _decball_contains_nonpositive(const decball_t x, gr_ctx_t ctx)
              int decball_cmp(int * res, const decball_t x, const decball_t y, gr_ctx_t ctx)
              int decball_cmpabs(int * res, const decball_t x, const decball_t y, gr_ctx_t ctx)
              int decball_sgn(decball_t res, const decball_t x, gr_ctx_t ctx)

    Predicates. Containment and overlap tests are decided exactly
    (using rational arithmetic in the rare cases where cheap bounds do
    not suffice), except that ``overlaps`` conservatively returns 1 and
    ``contains`` conservatively returns 0 when the exponents are too
    large for exact rational arithmetic.

.. function:: int decball_neg(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_abs(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_add(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx)
              int decball_sub(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx)
              int decball_mul(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx)
              int decball_sqr(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_div(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx)
              int decball_inv(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_sqrt(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_rsqrt(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_add_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx)
              int decball_add_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx)
              int decball_add_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx)
              int decball_sub_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx)
              int decball_sub_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx)
              int decball_sub_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx)
              int decball_mul_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx)
              int decball_mul_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx)
              int decball_mul_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx)
              int decball_div_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx)
              int decball_div_si(decball_t res, const decball_t x, slong y, gr_ctx_t ctx)
              int decball_div_fmpz(decball_t res, const decball_t x, const fmpz_t y, gr_ctx_t ctx)
              int decball_mul_two(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_mul_10exp_si(decball_t res, const decball_t x, slong e, gr_ctx_t ctx)
              int decball_mul_10exp_fmpz(decball_t res, const decball_t x, const fmpz_t e, gr_ctx_t ctx)
              int decball_mul_2exp_si(decball_t res, const decball_t x, slong e, gr_ctx_t ctx)
              int decball_mul_2exp_fmpz(decball_t res, const decball_t x, const fmpz_t e, gr_ctx_t ctx)
              int decball_floor(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_ceil(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_trunc(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_nint(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_pow_ui(decball_t res, const decball_t x, ulong y, gr_ctx_t ctx)

    Ball arithmetic. The radius of a product is
    `|x_m| y_r + |y_m| x_r + x_r y_r` plus the rounding error of the
    midpoint, the radius of a quotient is
    `(|x_m| y_r + |y_m| x_r) / (|y_m| (|y_m| - y_r))`, and so on,
    all evaluated with upper-bound rounding. Square roots of balls
    containing negative numbers return ``GR_UNABLE`` (``GR_DOMAIN`` if
    the whole ball is negative). Rounding to an integer returns an
    exact integer when the ball contains no integer boundary and
    otherwise `[\operatorname{round}(m) \pm (r + 1)]`.

.. function:: int decball_add_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx)
              int decball_sub_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx)
              int decball_mul_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx)
              int decball_div_round(decball_t res, const decball_t x, const decball_t y, slong prec, gr_ctx_t ctx)
              int decball_inv_round(decball_t res, const decball_t x, slong prec, gr_ctx_t ctx)
              int decball_sqrt_round(decball_t res, const decball_t x, slong prec, gr_ctx_t ctx)

    Versions with an explicit midpoint precision.

.. function:: int decball_add_error_decmag(decball_t res, const decmag_t err, gr_ctx_t ctx)
              int decball_add_error_10exp_si(decball_t res, slong e, gr_ctx_t ctx)
              int decball_add_error_decfloat(decball_t res, const decfloat_t err, gr_ctx_t ctx)
              int decball_trim(decball_t res, const decball_t x, gr_ctx_t ctx)

    Increases the radius, respectively rounds the midpoint to the number
    of digits justified by the radius.

.. function:: slong decball_rel_accuracy_digits(const decball_t x, gr_ctx_t ctx)

    Returns the number of correct significant digits of the midpoint,
    computed as the difference between the scientific exponents of the
    midpoint and the radius (the exponent of a zero midpoint is taken as
    zero). Returns ``DECIMAL_PREC_EXACT`` for an exact ball and
    ``-DECIMAL_PREC_EXACT`` for a ball with an infinite radius.

.. function:: int decball_pi(decball_t res, gr_ctx_t ctx)
              int decball_exp(decball_t res, const decball_t x, gr_ctx_t ctx)
              int decball_pow(decball_t res, const decball_t x, const decball_t y, gr_ctx_t ctx)
              int decball_via_arb(decball_t res, const decball_t x, int (*func)(arb_t, const arb_t, slong), gr_ctx_t ctx)
              int _decball_arb_gr(decball_ptr * res, slong nres, decball_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx)

    Every function listed for ``decfloat`` under *Elementary and special
    functions* is also available for balls, with the same name and
    signature (``decball_`` instead of ``decfloat_``), and so are the
    other special functions through the generic-ring methods. They are computed
    by a rigorous roundtrip through ``arb``: the input balls are
    converted to ``arb`` balls, the ``gr`` method of the ``arb`` ring is
    evaluated at a precision corresponding to the context precision,
    and the results are converted back, so that the output balls contain
    the true values for every point of the input balls. The two last
    functions expose the mechanism for a user-supplied ``arb`` function
    or an arbitrary ``gr`` method of the ``arb`` ring.

.. function:: char * decball_get_str(const decball_t x, gr_ctx_t ctx)
              int decball_write(gr_stream_t out, const decball_t x, gr_ctx_t ctx)
              int decball_randtest(decball_t res, flint_rand_t state, gr_ctx_t ctx)

    Output in the form ``[mid +/- rad]`` (or just ``mid`` for exact
    balls), and random generation of finite balls.

Complex floating-point numbers
--------------------------------------------------------------------------------

.. type:: deccfloat_struct
          deccfloat_t

    A complex number `a + bi` stored as a pair of real (field ``re``) and
    imaginary (field ``im``) parts of type :type:`decfloat_struct`, in the
    same way as ``acf``. The real and imaginary parts are rounded
    separately to the context precision, the real part in the mode
    ``DECIMAL_CTX_RND`` and the imaginary part in the mode
    ``DECIMAL_CTX_RND_IM`` (which is the same mode unless changed
    with ``decimal_ctx_set_rnd_im``).

    Special values are per component: with ``DECIMAL_ALLOW_INF`` and
    ``DECIMAL_ALLOW_NAN`` a number can be, say, `\infty + 2i` or
    `1 + \mathrm{NaN} \cdot i`, and a number is finite if both parts are.
    A real number is one with a zero imaginary part; there are no signed
    zeros.

    The ring is presented to the generics interface in the same way as the
    real floating-point ring (the exact ring `\mathbb{Z}[1/10][i]` when
    the precision is exact and no special values are admitted, and
    otherwise a non-ring), with
    ``GR_METHOD_I``, ``GR_METHOD_CONJ``, ``GR_METHOD_RE``, ``GR_METHOD_IM``,
    ``GR_METHOD_ABS``, ``GR_METHOD_SGN``, ``GR_METHOD_CSGN`` and
    ``GR_METHOD_ARG`` implemented. Comparison (``cmp``) is only defined
    for real numbers (``GR_DOMAIN`` otherwise), while ``cmpabs`` compares
    the exactly computed squared moduli. Rounding to integers
    (``floor``, ``ceil``, ``trunc``, ``nint``) likewise requires real
    input.

.. macro:: DECCFLOAT_REALREF(z)
           DECCFLOAT_IMAGREF(z)

    Pointers to the real and imaginary parts.

.. function:: decfloat_ptr deccfloat_realref(deccfloat_t z)
              decfloat_ptr deccfloat_imagref(deccfloat_t z)

    Function versions of the above, for use through foreign function
    interfaces.

.. function:: void deccfloat_init(deccfloat_t res, gr_ctx_t ctx)
              void deccfloat_clear(deccfloat_t res, gr_ctx_t ctx)
              void deccfloat_swap(deccfloat_t x, deccfloat_t y, gr_ctx_t ctx)
              void deccfloat_set_shallow(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int _deccfloat_is_real(const deccfloat_t x)
              int _deccfloat_is_imaginary(const deccfloat_t x)
              int _deccfloat_is_finite(const deccfloat_t x)
              int _deccfloat_is_nan(const deccfloat_t x)

.. function:: int deccfloat_zero(deccfloat_t res, gr_ctx_t ctx)
              int deccfloat_one(deccfloat_t res, gr_ctx_t ctx)
              int deccfloat_neg_one(deccfloat_t res, gr_ctx_t ctx)
              int deccfloat_i(deccfloat_t res, gr_ctx_t ctx)
              int deccfloat_nan(deccfloat_t res, gr_ctx_t ctx)
              int deccfloat_pos_inf(deccfloat_t res, gr_ctx_t ctx)
              int deccfloat_neg_inf(deccfloat_t res, gr_ctx_t ctx)

    Constants. The infinities are the real infinities `\pm \infty + 0i`,
    and NaN has both parts indeterminate.

.. function:: int deccfloat_set(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_set_round(deccfloat_t res, const deccfloat_t x, slong prec, int rnd, int rnd_im, gr_ctx_t ctx)
              int deccfloat_set_decfloat(deccfloat_t res, const decfloat_t x, gr_ctx_t ctx)
              int deccfloat_set_decfloat_decfloat(deccfloat_t res, const decfloat_t re, const decfloat_t im, gr_ctx_t ctx)
              int deccfloat_set_si(deccfloat_t res, slong x, gr_ctx_t ctx)
              int deccfloat_set_ui(deccfloat_t res, ulong x, gr_ctx_t ctx)
              int deccfloat_set_fmpz(deccfloat_t res, const fmpz_t x, gr_ctx_t ctx)
              int deccfloat_set_fmpq(deccfloat_t res, const fmpq_t x, gr_ctx_t ctx)
              int deccfloat_set_d(deccfloat_t res, double x, gr_ctx_t ctx)
              int deccfloat_set_str(deccfloat_t res, const char * s, gr_ctx_t ctx)
              int deccfloat_set_other(deccfloat_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int deccfloat_set_fmpz_10exp_fmpz(deccfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
              int deccfloat_set_fmpz_2exp_fmpz(deccfloat_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
              int deccfloat_set_acf(deccfloat_t res, const acf_t x, gr_ctx_t ctx)
              int deccfloat_set_acb(deccfloat_t res, const acb_t x, gr_ctx_t ctx)
              int deccfloat_get_acf(acf_t res, const deccfloat_t x, slong prec_bits, int rnd, gr_ctx_t ctx)
              int deccfloat_get_acb(acb_t res, const deccfloat_t x, slong prec_bits, gr_ctx_t ctx)
              int deccfloat_get_fmpz(fmpz_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_get_fmpq(fmpq_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_get_si(slong * res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_get_ui(ulong * res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_get_d(double * res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_get_re(decfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_get_im(decfloat_t res, const deccfloat_t x, gr_ctx_t ctx)

    Assignment and conversions, with the parts rounded to the context
    precision in the respective rounding modes. ``set_other`` accepts
    elements of all the decimal rings (a ball is converted to its
    midpoint, or returns ``GR_UNABLE`` when inexact and the target
    precision is exact), of ``acf``,
    ``acb``, ``arf``, ``arb``, ``qqbar`` (correctly rounded, see
    :func:`deccfloat_set_qqbar`) and of the rings accepted by
    ``decfloat_set_other``; ``deccfloat`` and ``deccball`` are conversely
    accepted by ``gr_set_other`` in the ``fmpz``, ``fmpq``, ``arf``, ``arb``,
    ``acf`` and ``acb`` rings and in all the decimal rings. The conversions
    to real numbers or integers return ``GR_DOMAIN`` for a nonzero
    imaginary part. Strings of the form ``"a"``, ``"b*I"``,
    ``"(a + b*I)"`` and ``"(a - b*I)"`` with plain decimal literals *a*
    and *b* (the output format) are parsed and rounded directly; any
    other string, such as ``"1+2*I"``, ``"I^2"`` or ``"(3+4*I)^(1/2)"``,
    goes through the generic expression parser with ``I`` bound to the
    imaginary unit and every operation rounded (with an exact result in
    the examples above).

.. function:: char * deccfloat_get_str(const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_write(gr_stream_t out, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_randtest(deccfloat_t res, flint_rand_t state, gr_ctx_t ctx)
              int deccfloat_randtest_special(deccfloat_t res, flint_rand_t state, gr_ctx_t ctx)

    Output in the form ``a`` for real numbers, ``b*I`` for imaginary
    numbers and ``(a + b*I)`` or ``(a - b*I)`` otherwise, where the parts
    are printed like real floating-point numbers, and random generation.

.. function:: truth_t deccfloat_is_zero(const deccfloat_t x, gr_ctx_t ctx)
              truth_t deccfloat_is_one(const deccfloat_t x, gr_ctx_t ctx)
              truth_t deccfloat_is_neg_one(const deccfloat_t x, gr_ctx_t ctx)
              truth_t deccfloat_equal(const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
              truth_t deccfloat_is_integer(const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_cmp(int * res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
              int deccfloat_cmpabs(int * res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)

    Predicates and comparisons. Equality is decided exactly (``T_UNKNOWN``
    only for NaN). ``cmp`` returns ``GR_DOMAIN`` unless both numbers are
    real, and ``cmpabs`` compares `|x|^2` and `|y|^2` exactly.

.. function:: int deccfloat_neg(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_conj(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_re(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_im(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_abs(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_arg(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_sgn(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_csgn(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_add(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
              int deccfloat_sub(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
              int deccfloat_mul(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
              int deccfloat_sqr(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_div(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
              int deccfloat_inv(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_sqrt(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_rsqrt(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_add_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx)
              int deccfloat_add_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx)
              int deccfloat_add_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccfloat_sub_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx)
              int deccfloat_sub_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx)
              int deccfloat_sub_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccfloat_mul_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx)
              int deccfloat_mul_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx)
              int deccfloat_mul_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccfloat_div_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx)
              int deccfloat_div_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx)
              int deccfloat_div_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccfloat_mul_decfloat(deccfloat_t res, const deccfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int deccfloat_div_decfloat(deccfloat_t res, const deccfloat_t x, const decfloat_t y, gr_ctx_t ctx)
              int deccfloat_mul_i(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_div_i(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_mul_two(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_mul_10exp_si(deccfloat_t res, const deccfloat_t x, slong e, gr_ctx_t ctx)
              int deccfloat_mul_10exp_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t e, gr_ctx_t ctx)
              int deccfloat_mul_2exp_si(deccfloat_t res, const deccfloat_t x, slong e, gr_ctx_t ctx)
              int deccfloat_mul_2exp_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t e, gr_ctx_t ctx)
              int deccfloat_floor(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_ceil(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_trunc(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_nint(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccfloat_pow_ui(deccfloat_t res, const deccfloat_t x, ulong y, gr_ctx_t ctx)
              int deccfloat_pow_si(deccfloat_t res, const deccfloat_t x, slong y, gr_ctx_t ctx)
              int deccfloat_pow_fmpz(deccfloat_t res, const deccfloat_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccfloat_pow(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
              int deccfloat_vec_dot(deccfloat_t res, const deccfloat_t initial, int subtract, deccfloat_srcptr vec1, deccfloat_srcptr vec2, slong len, gr_ctx_t ctx)
              int deccfloat_vec_dot_rev(deccfloat_t res, const deccfloat_t initial, int subtract, deccfloat_srcptr vec1, deccfloat_srcptr vec2, slong len, gr_ctx_t ctx)

    Arithmetic. Every operation is *componentwise correctly rounded*:
    the real and imaginary parts of the exact result of the operation
    on the exact inputs are each rounded once, to the context precision
    in the respective rounding modes. In particular, a product is computed
    as two exact dot products `x_r y_r - x_i y_i` and `x_r y_i + x_i y_r`
    each rounded once (the ``vec_dot`` methods, on the other hand, round
    after each term); a quotient `x / y` is computed by forming the exact
    numerators `x \bar{y}` and denominator `|y|^2` followed by a single
    correctly rounded division for each part (with direct divisions
    when `y` is real or imaginary); `\sqrt{z}`, `1/\sqrt{z}`, `|z|`,
    `\operatorname{sgn}(z) = z/|z|` and `\arg z` are exact whenever the
    result is representable (as in `\sqrt{3 + 4i} = 2 + i`) and otherwise
    computed by Ziv's strategy on top of ``acb`` as described below. The
    absolute value is computed as the correctly rounded square root of
    the exact `x_r^2 + x_i^2`, and `\arg z` as ``decfloat_atan2``.
    Integer powers are exact when the exact result has at most about
    `10^5` digits (as in `(1+i)^{1000000}`, which is a power of `2` times
    a power of `i`), and otherwise correctly rounded through ``acb``;
    general powers `x^y` reduce to integer powers, square roots, real
    powers and the exponential function as appropriate.

    Zeros in the input are treated as exact zeros: a zero factor
    annihilates its term even when the other factor is infinite, so that
    `i \cdot (\infty + 0i) = \infty i` rather than
    `\mathrm{NaN} + \infty i`, and dividing a number with a zero part by a
    real number keeps that part zero (`i/0 = \infty i` rather than
    `\mathrm{NaN} + \infty i` when infinities are admitted).

.. function:: int _decfloat_dot2(decfloat_t res, const decfloat_t x1, const decfloat_t y1, const decfloat_t x2, const decfloat_t y2, int subtract, slong prec, int rnd, decimal_rounding_info * info, decmag_ptr err, gr_ctx_t ctx)

    Sets *res* to `x_1 y_1 \pm x_2 y_2` correctly rounded (with the
    exact-zero convention above), optionally recording the rounding
    error. This is the primitive behind complex multiplication and
    division for both floating-point numbers and balls.

.. function:: int _deccfloat_mul_exact(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, gr_ctx_t ctx)
              int _deccfloat_add_exact(deccfloat_t res, const deccfloat_t x, const deccfloat_t y, int subtract, gr_ctx_t ctx)
              int _deccfloat_inv_exact(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int _deccfloat_pow_ui_exact(deccfloat_t res, const deccfloat_t x, ulong n, gr_ctx_t ctx)
              int _deccfloat_sqrt_exact(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)

    Exact operations, returning ``GR_UNABLE`` when the result is not
    representable (or, for additions and powers, would be absurdly large).

Complex balls
--------------------------------------------------------------------------------

.. type:: deccball_struct
          deccball_t

    A complex ball stored as a pair of real balls (fields ``re`` and
    ``im`` of type :type:`decball_struct`), representing the rectangle
    `\{ a + bi : |a - m_r| \le r_r, |b - m_i| \le r_i \}` as in ``acb``.
    Operations produce balls containing the exact result for every point
    of the input balls, with midpoints computed to the context precision
    in the context rounding mode (the same mode for both parts).

    The ring is presented to the generics interface as an inexact,
    algebraically closed field: predicates return ``T_UNKNOWN`` when the
    balls do not decide the answer, ``cmp`` is only defined for real
    balls, and division by a ball containing zero returns ``GR_UNABLE``.

.. macro:: DECCBALL_REALREF(z)
           DECCBALL_IMAGREF(z)

    Pointers to the real and imaginary parts.

.. function:: decball_ptr deccball_realref(deccball_t z)
              decball_ptr deccball_imagref(deccball_t z)

    Function versions of the above, for use through foreign function
    interfaces.

.. function:: void deccball_init(deccball_t res, gr_ctx_t ctx)
              void deccball_clear(deccball_t res, gr_ctx_t ctx)
              void deccball_swap(deccball_t x, deccball_t y, gr_ctx_t ctx)
              void deccball_set_shallow(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int _deccball_is_exact(const deccball_t x, gr_ctx_t ctx)
              int _deccball_is_finite(const deccball_t x, gr_ctx_t ctx)
              int _deccball_is_real(const deccball_t x, gr_ctx_t ctx)

.. function:: int deccball_zero(deccball_t res, gr_ctx_t ctx)
              int deccball_one(deccball_t res, gr_ctx_t ctx)
              int deccball_neg_one(deccball_t res, gr_ctx_t ctx)
              int deccball_i(deccball_t res, gr_ctx_t ctx)

.. function:: int deccball_set(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_set_round(deccball_t res, const deccball_t x, slong prec, gr_ctx_t ctx)
              int deccball_set_round2(deccball_t res, const deccball_t x, slong prec, slong rad_prec, gr_ctx_t ctx)
              int deccball_set_decball(deccball_t res, const decball_t x, gr_ctx_t ctx)
              int deccball_set_decball_decball(deccball_t res, const decball_t re, const decball_t im, gr_ctx_t ctx)
              int deccball_set_decfloat(deccball_t res, const decfloat_t x, gr_ctx_t ctx)
              int deccball_set_deccfloat(deccball_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccball_set_si(deccball_t res, slong x, gr_ctx_t ctx)
              int deccball_set_ui(deccball_t res, ulong x, gr_ctx_t ctx)
              int deccball_set_fmpz(deccball_t res, const fmpz_t x, gr_ctx_t ctx)
              int deccball_set_fmpq(deccball_t res, const fmpq_t x, gr_ctx_t ctx)
              int deccball_set_d(deccball_t res, double x, gr_ctx_t ctx)
              int deccball_set_str(deccball_t res, const char * s, gr_ctx_t ctx)
              int deccball_set_other(deccball_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
              int deccball_set_fmpz_10exp_fmpz(deccball_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
              int deccball_set_fmpz_2exp_fmpz(deccball_t res, const fmpz_t m, const fmpz_t e, gr_ctx_t ctx)
              int deccball_set_acb(deccball_t res, const acb_t x, gr_ctx_t ctx)
              int deccball_get_acb(acb_t res, const deccball_t x, slong prec_bits, gr_ctx_t ctx)
              int deccball_get_fmpz(fmpz_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_get_fmpq(fmpq_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_get_si(slong * res, const deccball_t x, gr_ctx_t ctx)
              int deccball_get_ui(ulong * res, const deccball_t x, gr_ctx_t ctx)
              int deccball_get_d(double * res, const deccball_t x, gr_ctx_t ctx)
              int deccball_get_mid(deccfloat_t res, const deccball_t x, gr_ctx_t ctx)

    Assignment and conversions, following ``decball``. Strings in the
    output format (see below) and expressions such as
    ``"[1 + 2*I +/- 0.01]"``, ``"[1 + 2*I +/- (0.01 + 0.02*I)]"``
    or ``"([1 +/- 0.1] + [2 +/- 0.2]*I) / (3 - I)"`` are accepted, also
    inside strings for polynomials and matrices over complex balls.

.. function:: int deccball_set_interval_mid_rad(deccball_t res, const deccball_t m, const deccball_t r, gr_ctx_t ctx)
              int deccball_set_interval_mid_inf(deccball_t res, const deccball_t m, gr_ctx_t ctx)

    Sets *res* to the ball *m* with the radius of the real part increased
    by an upper bound for `|\operatorname{Re}(r)|` and the radius of the
    imaginary part by an upper bound for `|\operatorname{Im}(r)|` (as in
    ``acb``, so that a real *r* only widens the real part), respectively
    with an infinite radius of the real part. These implement
    the ``+/-`` operator of the generic expression parser.

.. function:: int deccball_mid(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_shell(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_add_rad(deccball_t res, const deccball_t x, const deccball_t r, gr_ctx_t ctx)

    The exact midpoint, the ball with the same radii centered at zero, and
    the same as ``set_interval_mid_rad``, respectively.

.. function:: char * deccball_get_str(const deccball_t x, gr_ctx_t ctx)
              int deccball_write(gr_stream_t out, const deccball_t x, gr_ctx_t ctx)
              int deccball_randtest(deccball_t res, flint_rand_t state, gr_ctx_t ctx)

    Output in the form ``a``, ``b*I``, ``(a + b*I)`` or ``(a - b*I)`` where
    the parts are printed as real balls (``[mid +/- rad]``, or just
    ``mid`` when exact), and random generation.

.. function:: truth_t deccball_is_zero(const deccball_t x, gr_ctx_t ctx)
              truth_t deccball_is_one(const deccball_t x, gr_ctx_t ctx)
              truth_t deccball_is_neg_one(const deccball_t x, gr_ctx_t ctx)
              truth_t deccball_equal(const deccball_t x, const deccball_t y, gr_ctx_t ctx)
              truth_t deccball_is_integer(const deccball_t x, gr_ctx_t ctx)
              int _deccball_contains_zero(const deccball_t x, gr_ctx_t ctx)
              int _deccball_contains(const deccball_t x, const deccball_t y, gr_ctx_t ctx)
              int _deccball_overlaps(const deccball_t x, const deccball_t y, gr_ctx_t ctx)
              int _deccball_contains_deccfloat(const deccball_t x, const deccfloat_t y, gr_ctx_t ctx)
              int deccball_cmp(int * res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
              int deccball_cmpabs(int * res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)

    Predicates, applied to both parts.

.. function:: int deccball_neg(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_conj(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_re(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_im(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_abs(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_arg(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_sgn(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_csgn(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_add(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
              int deccball_sub(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
              int deccball_mul(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
              int deccball_sqr(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_div(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)
              int deccball_inv(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_sqrt(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_rsqrt(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_add_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx)
              int deccball_add_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx)
              int deccball_add_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccball_sub_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx)
              int deccball_sub_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx)
              int deccball_sub_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccball_mul_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx)
              int deccball_mul_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx)
              int deccball_mul_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccball_div_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx)
              int deccball_div_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx)
              int deccball_div_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccball_mul_decball(deccball_t res, const deccball_t x, const decball_t y, gr_ctx_t ctx)
              int deccball_div_decball(deccball_t res, const deccball_t x, const decball_t y, gr_ctx_t ctx)
              int deccball_mul_i(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_div_i(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_mul_two(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_mul_10exp_si(deccball_t res, const deccball_t x, slong e, gr_ctx_t ctx)
              int deccball_mul_10exp_fmpz(deccball_t res, const deccball_t x, const fmpz_t e, gr_ctx_t ctx)
              int deccball_mul_2exp_si(deccball_t res, const deccball_t x, slong e, gr_ctx_t ctx)
              int deccball_mul_2exp_fmpz(deccball_t res, const deccball_t x, const fmpz_t e, gr_ctx_t ctx)
              int deccball_floor(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_ceil(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_trunc(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_nint(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              int deccball_pow_ui(deccball_t res, const deccball_t x, ulong y, gr_ctx_t ctx)
              int deccball_pow_si(deccball_t res, const deccball_t x, slong y, gr_ctx_t ctx)
              int deccball_pow_fmpz(deccball_t res, const deccball_t x, const fmpz_t y, gr_ctx_t ctx)
              int deccball_pow(deccball_t res, const deccball_t x, const deccball_t y, gr_ctx_t ctx)

    Ball arithmetic. Products are computed with :func:`_decball_dot2`,
    which forms the exact midpoint products so that the radius only
    accounts for the input radii and a single rounding per part; the same
    holds for quotients when the modulus of the divisor is bounded away
    from zero by the componentwise bounds, and otherwise (when the
    rectangle `|y|^2` computed from the parts contains zero but `y` does
    not) the quotient is computed through ``acb``. Square roots, absolute
    values, arguments and signs of non-real balls go through ``acb`` as
    the transcendental functions do.

.. function:: int _decball_dot2(decball_t res, const decball_t x1, const decball_t y1, const decball_t x2, const decball_t y2, int subtract, slong prec, gr_ctx_t ctx)

    Sets *res* to a ball containing `x_1 y_1 \pm x_2 y_2`, with the
    midpoint computed as the correctly rounded exact product sum of the
    midpoints.

.. function:: int deccball_add_error_decmag(deccball_t res, const decmag_t err, gr_ctx_t ctx)
              int deccball_trim(deccball_t res, const deccball_t x, gr_ctx_t ctx)
              slong deccball_rel_accuracy_digits(const deccball_t x, gr_ctx_t ctx)

    Increases both radii, trims both parts, respectively returns the
    number of correct significant digits computed from the larger
    midpoint and the larger radius.

Complex functions
--------------------------------------------------------------------------------

.. function:: int deccfloat_exp(deccfloat_t res, const deccfloat_t x, gr_ctx_t ctx)
              int deccball_exp(deccball_t res, const deccball_t x, gr_ctx_t ctx)

    Every function listed for ``decfloat`` under *Elementary and special
    functions* is available for complex floating-point numbers and complex
    balls with the same name and signature (``deccfloat_`` or
    ``deccball_`` in place of ``decfloat_``), except for ``atan2``
    (use ``arg``). All the other special functions, including the
    complex-only ones marked in the list there, are available through the
    generic-ring methods, and the ``gr`` methods ``GR_METHOD_I``,
    ``GR_METHOD_CONJ`` and so on are implemented.

    For floating-point numbers, the result of a function is
    *componentwise correctly rounded* in the sense defined above:
    the real and the imaginary part of the exact value are each rounded
    once, to the context precision in the respective rounding modes.
    Each function is evaluated by the first applicable of the following
    methods:

    * *Real and imaginary arguments.* For a real argument `t`, the real
      ``decfloat`` function is called (inheriting its exact cases and
      special-argument handling); if it fails because the value is not
      real, the standard reductions are applied: `\log(-t) = \log t + \pi
      i`, `\operatorname{acos} t = \pm i\, \operatorname{acosh} t` and
      `\pi - i \operatorname{acosh}(-t)` for `t > 1` and `t < -1`,
      `\operatorname{acosh} t = i \operatorname{acos} t` for `-1 \le t <
      1` and `\operatorname{acosh}(-t) + \pi i` for `t < -1`,
      `\exp(\pi i t) = \cos(\pi t) + i \sin(\pi t)`, and `\sqrt{-t} = i
      \sqrt{t}`. For an imaginary argument `it`, the corresponding real
      function is used where one exists: `\sin(it) = i \sinh t`,
      `\cos(it) = \cosh t`, `\tan(it) = i \tanh t`, `\cot(it) = -i \coth
      t`, `\sec`, `\csc`, `\sinh`, `\cosh`, `\tanh`, `\coth`,
      `\operatorname{sech}`, `\operatorname{csch}`,
      `\operatorname{asin}(it) = i \operatorname{asinh} t`,
      `\operatorname{atan}(it) = i \operatorname{atanh} t` (for `|t| <
      1`), `\operatorname{asinh}`, `\operatorname{atanh}`,
      `\operatorname{erf}(it) = i \operatorname{erfi} t`,
      `\operatorname{erfi}(it) = i \operatorname{erf} t`,
      `\operatorname{Si}(it) = i \operatorname{Shi} t` and
      `\operatorname{Shi}(it) = i \operatorname{Si} t`, together with
      `\exp(it) = \cos t + i \sin t`. For functions of several arguments
      the real function is used when all arguments are real.
      Consequently, all the exact cases and tiny and large argument
      handling of the real functions carry over to real and imaginary
      arguments: `\sqrt{-4} = 2i`, `\log(-1) = \pi i`, `\Gamma(5) = 24`,
      `\sin(10^{-10^9} i) = i \sinh(10^{-10^9})` and `\tanh(10^{30} + 0i)
      = 1` are all correctly rounded at negligible cost.
    * *Tiny arguments.* For `\exp`, `\operatorname{expm1}`, `\log(1 +
      z)`, `\log` near `1`, `\sin`, `\cos`, `\tan`, `\sinh`, `\cosh`,
      `\tanh`, `\operatorname{asin}`, `\operatorname{atan}`,
      `\operatorname{asinh}`, `\operatorname{atanh}`,
      `\operatorname{sinc}`, `\sec`, `\operatorname{sech}`, `\cot`,
      `\csc`, `\coth`, `\operatorname{csch}`, `\Gamma`, `1/\Gamma`,
      `\psi`, `\zeta` near `1`, `W`, `\operatorname{Li}_2`, `\cos(\pi
      z)`, `\operatorname{sinc}(\pi z)`, `\operatorname{Si}`,
      `\operatorname{Shi}` and, through `1/z`, `\operatorname{acot}`,
      `\operatorname{acoth}`, `\operatorname{acsc}` and
      `\operatorname{acsch}` at large arguments, a complex version of the
      tiny-argument mechanism of the real functions is used. The leading
      part `S` (`z`, `1`, `\pm 1/z`, `z - 1` or `1/(z-1)`) is computed
      exactly, and the sign of each component of the remainder
      `f(z) - S` is determined from the first terms of the Taylor
      expansion, whose powers `z^k` are computed in ball arithmetic
      (exactly when representable); a component of `S` which is zero is
      replaced by the first nonzero Taylor term. Each component of the
      result is then `S` rounded with the remainder acting as a sticky
      bit, provided that the remainder lies more than ``prec`` `+ 4e`
      digits below `S`, which is the case whenever the rounding could
      otherwise not be decided. This works for arbitrarily unbalanced
      arguments; for example, `\sin(10^{-10^9} + 10^{-2 \cdot 10^9} i)`
      rounded toward zero has both parts one unit in the last place
      below `10^{-10^9}` and `10^{-2 \cdot 10^9}` respectively, while
      `\exp` of the same argument is `1 + 10^{-2 \cdot 10^9} i` when
      rounded down and has both parts one unit in the last place above
      these values when rounded up. Components whose
      rounding is not decided by this mechanism fall through to the
      next method, with the decided components kept.
    * *Ziv's strategy through* ``acb``. The argument is converted to an
      ``acb`` ball, the corresponding ``gr`` method of the ``acb`` ring is
      evaluated at increasing precision until the rounding of each
      (remaining) component is decided, and the parts are rounded
      separately. ``GR_UNABLE`` is returned if this does not happen at
      64 times the working precision.

    The last method is the only one used for the special functions
    beyond those listed, and it is also all that the ball versions need:
    for complex balls, every function is computed by a rigorous
    roundtrip through ``acb`` (except that the componentwise operations
    of the previous section and integer powers are computed directly),
    so that the output balls contain the true values for every point of
    the input balls. The ball versions inherit the limitations of ``acb``;
    for instance, `\exp(10^{1000} + i/2)` overflows to an indeterminate
    ball, while the floating-point version returns the correctly rounded
    result.

    The functions above are correctly rounded whenever they succeed. The
    caveats are the same as for the real functions: Ziv's strategy
    returns ``GR_UNABLE`` when the exact value of a component is a
    representable number or extremely close to one and the situation is
    not detected by the reductions above (``acb`` itself computes many
    exact values of the special functions at non-real arguments exactly,
    such as `\operatorname{Li}_0(1+i) = -1 + i`, `T_2(i) = -3` and
    `\zeta(0, 1+i) = -1/2 - i`, in which case Ziv's strategy succeeds).
    Known examples are components converging exponentially fast to a
    representable number, as `\operatorname{erf}(300 + i)` whose real
    part is `1 - O(10^{-39000})` (`\tanh`, `\coth`, `\operatorname{expm1}`
    and `\zeta` at large real parts are handled by ``acb`` itself, but
    `\operatorname{erf}` and `\operatorname{erfc}` at large arguments are
    only handled for real arguments), and functions for which the
    ``acb`` ring provides no method for non-real arguments
    (``sec_pi``, ``asin_pi``, ``acos_pi``, ``atan_pi``, ``acot_pi``,
    ``asec_pi``, ``acsc_pi``, ``erfinv``, ``erfcinv``), which return
    ``GR_UNABLE`` in those cases. Operations
    reached through the generic fallbacks of the ``gr`` interface (vector
    dot products and hence polynomial and matrix arithmetic, ``gr_hypot``,
    ``gr_pow_fmpq``, expressions parsed by ``gr_set_str``, and so on) are
    not correctly rounded, as for the real types.

.. function:: int _deccfloat_acb_gr(deccfloat_ptr * res, slong nres, deccfloat_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx)
              int _deccball_acb_gr(deccball_ptr * res, slong nres, deccball_srcptr * args, slong nargs, int has_flag, int flag, int method, gr_ctx_t ctx)

    The machinery behind the complex functions: evaluates the ``gr``
    method *method* of the ``acb`` ring applied to *nargs* arguments
    (with an integer flag if *has_flag* is set) producing *nres* results,
    with Ziv's strategy for floating-point numbers and a single rigorous
    evaluation for balls, corresponding to :func:`_decfloat_arb_gr` and
    :func:`_decball_arb_gr`.
