.. _mp_real:

**mp_real.h** -- low-level limb-radix multiple-precision real arithmetic
===============================================================================

This module provides multiple-precision real arithmetic in the radix
`B = 2^{\mathrm{FLINT\_BITS}}`, with low overhead and limb granularity
throughout, intended as a backend for arbitrary-precision numerical
algorithms: binary splitting, AGM iterations, Newton iterations and the
evaluation of elementary functions and constants.  It works at two
levels:

* :type:`mp_real_t`, a ball (a midpoint with a limb mantissa and a
  radix-`B` exponent, and a rigorous radius), with arithmetic, roots,
  the arithmetic-geometric mean, constants and hypergeometric series.
  These functions start with ``mp_real_``.

* Fixed-point numbers stored in limb arrays, with the elementary
  functions, Newton inverses and square roots, and verified floors of
  constants.  These functions start with ``_mp_real_``, as does any
  function with preconditions it does not check.

The module is mainly optimized for 64-bit systems.  With 32-bit limbs,
some generated straight-line and register implementations are disabled
and evaluation goes through generic code paths.

Conventions
-------------------------------------------------------------------------------

**Fixed-point numbers.**  A fixed-point number `(x, n)` is an unsigned
`n`-limb fraction ``x[0], ..., x[n-1]`` representing
`\sum_i x_i B^{i - n}`, so `0 \le x < 1`, with unit in the last place
(ulp) `B^{-n}`.  Outputs of size `n + 1` additionally carry an integer
(units) limb at index `n`, and outputs of size `n + 2` two integral
limbs.

**Error bounds.**  The fixed-point functions with an ulp-accurate
contract take a pointer ``ulong * err`` right after their outputs.
They write a bound on `|\mathrm{res} - f(x)|` in ulps of the output
(for functions with two outputs, one bound valid for each), and
*err* may be ``NULL`` when the caller does not need the bound.  Each
function documents the largest value it can write, which callers can
use to size their guard limbs.  Functions built on ball arithmetic
write the bound they compute at run time, usually well below that
maximum; the table-based and series functions write their a priori
bound.  Unless stated otherwise the bounds are two-sided.

The generated per-size kernels, the Newton inverses and roots (whose
bounds are relative to the output) and the exact floors of constants
take no *err*.

**Precision.**  All functions take their precision `n` in limbs.

Ball arithmetic
-------------------------------------------------------------------------------

Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. type:: mp_real_struct

.. type:: mp_real_t

    A ``mp_real_struct`` has fields ``d`` (the mantissa limbs),
    ``alloc``, ``size``, ``negative``, ``exp`` and ``err``, and
    represents the ball

    .. math::

        (-1)^{\mathrm{negative}} (d, \mathrm{size}) \, B^{\mathrm{exp} - \mathrm{size}}
        \;\pm\; \mathrm{err} \, B^{\mathrm{exp} - \mathrm{size}}.

    The radius is a count of ulps of the mantissa's bottom limb (for a
    zero mantissa, ``size == 0``, of `B^{\mathrm{exp}}`).  The top limb
    of the mantissa is nonzero, so
    `B^{\mathrm{exp} - 1} \le |\mathrm{mid}| < B^{\mathrm{exp}}`.
    A ``mp_real_t`` is defined as an array of length one of type
    ``mp_real_struct``, permitting a ``mp_real_t`` to be passed by
    reference.

    The count ``err`` is a single limb, `0 \le \mathrm{err} < B`, and
    `0` means exact.  The radius has no scale of its own:

    * a bound below one ulp is installed by padding the mantissa with
      zero limbs down to the bound's scale (balls like `1 - \epsilon`
      carry their zero limbs explicitly);
    * a bound of `B` ulps or more truncates the mantissa by as many
      limbs as needed, each dropped run of limbs being worth one more
      unit at the new bottom limb.

    The mantissa length thus tracks the number of accurate limbs, with
    at most one limb of noise retained.  Exact values carry no low zero
    limbs, so exact small integers stay small.

.. type:: mp_real_ptr

    Alias for ``mp_real_struct *``, used for vectors of balls.

.. type:: mp_real_srcptr

    Alias for ``const mp_real_struct *``, used for vectors of balls
    when passed as constant input to functions.

Exponent safety and working precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The arithmetic is rigorous except that it does not check for exponent
overflow or underflow; such checks would be overhead on every
operation.  Instead, a caller can check that its operands lie in a safe
range.

.. macro:: MP_REAL_EXP_MAX

    The largest safe limb exponent, `\mathrm{WORD\_MAX} / (4 \cdot
    \mathrm{FLINT\_BITS})`.  Since several functions form the bit
    exponent ``FLINT_BITS * exp``, this bounds the bit exponents by
    `\mathrm{WORD\_MAX} / 4`, which is `2^{61}` on 64-bit systems.

.. function:: int mp_real_exp_is_safe(slong e)
              int mp_real_is_safe(const mp_real_t x)

    Returns whether `|e| \le` :macro:`MP_REAL_EXP_MAX`, respectively
    whether the exponent of *x* is safe.  If all ball operands of an
    operation are safe and every bit-exponent argument `e` (of
    :func:`mp_real_mul_2exp_si`, :func:`mp_real_add_error_2exp_si`,
    :func:`mp_real_add_rel_error_2exp_si`, :func:`_mp_real_set_mpn_2exp`)
    satisfies `|e| \le \mathrm{FLINT\_BITS} \cdot` :macro:`MP_REAL_EXP_MAX`,
    the operation computes without overflow.  Its result can lie outside
    the safe range, so a chain of operations that may grow exponents
    should check again.

.. macro:: MP_REAL_PREC_GUARD

    The number of guard limbs added by the functions below, currently 1.

.. function:: slong mp_real_prec_limbs(slong prec_limbs)
              slong mp_real_prec_bits(slong prec_bits)

    Returns a working precision `n` (in limbs) such that the basic
    operations (:func:`mp_real_add`, :func:`mp_real_sub`,
    :func:`mp_real_mul`, :func:`mp_real_div`, :func:`mp_real_sqrt`,
    :func:`mp_real_rsqrt`, :func:`mp_real_mul_ui`,
    :func:`mp_real_div_ui`, :func:`mp_real_addmul_ui`) at precision `n`
    return a relative radius of at most `4 \varepsilon` on exact
    operands and at most `16 \varepsilon` on operands of relative
    radius at most `\varepsilon`, where `\varepsilon =
    B^{-\mathrm{prec\_limbs}}` resp. `2^{-\mathrm{prec\_bits}}`
    (sums of operands of like sign, and differences without
    cancellation).  The result is ``prec_limbs`` resp.
    `\lceil \mathrm{prec\_bits} / \mathrm{FLINT\_BITS} \rceil` plus
    :macro:`MP_REAL_PREC_GUARD`; without a guard limb a product of
    exact operands can lose a whole limb.

Memory management and assignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: void mp_real_init(mp_real_t x)
              void mp_real_clear(mp_real_t x)

    Initializes *x* to exact zero, respectively frees its memory.

.. function:: void mp_real_fit_length(mp_real_t x, slong len)
              void _mp_real_grow(mp_real_t x, slong len)

    Ensures that *x* has room for *len* mantissa limbs (the inline
    function calls the second one when it has to reallocate).

.. function:: void mp_real_swap(mp_real_t x, mp_real_t y)

    Swaps *x* and *y* efficiently.

.. function:: void mp_real_zero(mp_real_t x)
              void mp_real_set(mp_real_t res, const mp_real_t x)
              void mp_real_set_ui(mp_real_t x, ulong c)
              void mp_real_set_si(mp_real_t x, slong c)

    Sets the ball to zero, a copy of *x*, or the exact integer *c*.

.. function:: void mp_real_set_fmpz(mp_real_t x, const fmpz_t f)

    Sets *x* to the integer *f* exactly.  The limbs of *f* above its
    valuation (the low zero limbs) are copied directly into the mantissa,
    without temporaries.

.. function:: void mp_real_set_trunc(mp_real_t res, const mp_real_t x, slong n)

    Sets *res* to *x* rounded to its top `n \ge 1` limbs, for any input
    (a copy if *x* has at most *n* limbs).  The dropped tail and the old
    radius are each below one ulp of the new bottom limb, so the new
    radius is `(\text{tail} \ne 0) + (\text{radius} \ne 0)` ulps; the
    result is exact if nothing nonzero is dropped from an exact input.
    *res* may alias *x*.

.. function:: void _mp_real_set_mpn_2exp(mp_real_t x, nn_srcptr p, slong len, slong e)

    Sets *x* exactly to `(p, \mathrm{len}) \cdot 2^e` for the unsigned
    integer `(p, \mathrm{len})` (which may have zero top limbs).  The bit
    part of *e* costs one mpn shift, after which the exponent is
    limb-aligned.

.. function:: void mp_real_neg(mp_real_t res, const mp_real_t x)

    Sets *res* to `-x`.  In place this is a single bit flip.

.. function:: void mp_real_mul_2exp_si(mp_real_t res, const mp_real_t x, slong e)

    Sets *res* to `x \cdot 2^e` exactly (a limb shift of the exponent
    and at most one mpn shift of the mantissa).

Predicates, magnitudes and the radius
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: int mp_real_is_zero(const mp_real_t x)

    Returns whether *x* is exactly zero (zero midpoint and radius).

.. function:: slong mp_real_abs_bound_lt_2exp_si(const mp_real_t x)

    Returns an `e` with `|t| < 2^e` for every `t` in the ball.

.. function:: slong mp_real_rel_radius_lt_2exp_si(const mp_real_t x)

    Returns an `e` such that the radius of *x* is below `2^e` times its
    midpoint: ``-WORD_MAX / 2`` for an exact *x* and ``WORD_MAX / 2``
    for a zero midpoint with a nonzero radius.

.. function:: void mp_real_add_error_2exp_si(mp_real_t x, slong e)
              void mp_real_add_rel_error_2exp_si(mp_real_t x, slong e)

    Adds `2^e` to the radius of *x*, respectively enlarges *x* to contain
    `x (1 + t)` for every `|t| \le 2^e`.

.. function:: void _mp_real_add_error_ulps(mp_real_t x, double v)
              void _mp_real_add_error_ulps_at(mp_real_t x, double v, slong anchor)

    Adds `v` ulps of the mantissa's bottom limb, respectively `v B^{\mathrm{anchor}}`,
    to the radius of *x*.  The first form depends on the current
    normalization of the mantissa (an exact import strips low zero limbs,
    so ulps at the current anchor can mean a much coarser scale than
    intended); the second states the scale explicitly.

Arithmetic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These functions truncate their outputs to about `n` limbs, and fewer
when an operand is accurate to fewer; outputs may alias inputs.

.. function:: void mp_real_add(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
              void mp_real_sub(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)

    Sets *res* to `a + b` resp. `a - b`.

.. function:: void mp_real_addmul_ui(mp_real_t res, const mp_real_t x, const mp_real_t y, ulong c, slong n)
              void mp_real_submul_ui(mp_real_t res, const mp_real_t x, const mp_real_t y, ulong c, slong n)

    Sets *res* to `x + y c` resp. `x - y c` for a word `c`.  In place
    (*res* aliasing *x*) and with `y c` inside the window of *x*, this is
    one ``mpn_addmul_1`` resp. ``mpn_submul_1`` over the limbs of *y*,
    independent of the size of *x*.  These are the primitives of a
    rectangular splitting or of an interpolation formula summed into one
    accumulator.

.. function:: void mp_real_mul(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
              void mp_real_mul_ui(mp_real_t res, const mp_real_t a, ulong c, slong n)

    Sets *res* to `a b` resp. `a c`.

.. function:: void mp_real_mul_complex(mp_real_t rr, mp_real_t ri, const mp_real_t ar, const mp_real_t ai, const mp_real_t br, const mp_real_t bi, slong n)

    Sets `rr + i \, ri` to `(ar + i \, ai)(br + i \, bi)`, to about `n`
    limbs relative to the larger part of the result (both parts kept
    down to the same limb, the precision of a complex number being that
    of its modulus), by one transform-sharing complex product instead of
    four.

.. function:: void _mp_real_submul_bounded(mp_real_t res, const mp_real_t a, const mp_real_t b, const mp_real_t c, slong E, slong n)

    Sets *res* to `a - b c` GIVEN the bound `|a - b c| < B^E` from the
    caller (the residual of an iteration, whose bound comes from the
    radius of the previous approximation): only the window of the
    product that reaches the result is computed.  The result is wrong if
    the bound does not hold.

.. function:: void mp_real_div(mp_real_t res, const mp_real_t a, const mp_real_t b, slong n)
              void mp_real_div_ui(mp_real_t res, const mp_real_t a, ulong c, slong n)

    Sets *res* to `a / b` resp. `a / c`.  The divisor ball must have a
    relative radius below `2^{-30}` (checked; wider divisors throw), and
    `c \ne 0`.

Roots and the AGM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: void mp_real_sqrt(mp_real_t res, const mp_real_t x, slong n)
              void mp_real_rsqrt(mp_real_t res, const mp_real_t x, slong n)

    Sets *res* to `\sqrt{x}` resp. `1/\sqrt{x}` for `x > 0` of relative
    radius below `2^{-30}` (checked).  A ball with a negative midpoint or
    a zero midpoint throws, except that :func:`mp_real_sqrt` of an exact
    zero gives zero.  The exponent is aligned, the
    root computed by :func:`_mp_real_sqrt_newton` resp.
    :func:`_mp_real_rsqrt_newton` (or ``flint_mpn_sqrtrem`` below a
    cutoff) and the radius bounded through the derivative.  *res* may
    alias *x* (through a temporary).

.. function:: void mp_real_rsqrt_ui(mp_real_t res, ulong c, slong n)

    Sets *res* to `1/\sqrt{c}` for a word `c \ge 1` (throws for `c = 0`),
    by :func:`_mp_real_rsqrt_ui_newton`.

.. macro:: MP_REAL_ROOT_K_MAX

    The bound on `k` for the functions below: `2^{40}` on 64-bit
    machines and `2^{27}` on 32-bit machines.

.. function:: void mp_real_root_ui(mp_real_t res, const mp_real_t x, ulong k, slong n)
              void mp_real_rroot_ui(mp_real_t res, const mp_real_t x, ulong k, slong n)
              void _mp_real_root_ui_order(mp_real_t res, const mp_real_t x, ulong k, slong n, int r, int recip)

    Sets *res* to `x^{1/k}` respectively `x^{-1/k}` for `x > 0` and
    `1 \le k <` :macro:`MP_REAL_ROOT_K_MAX`.  For `k = 2` these call
    :func:`mp_real_sqrt` and :func:`mp_real_rsqrt`.

    For other `k` the roots come from a high-order iteration written in
    ``mp_real_t`` arithmetic, which handles the exponents and the
    propagation of the arithmetic errors.  With `z` any approximation
    of `x^{-1/k}` and `u = 1 - x z^k`, exactly

    .. math::

        x^{-1/k} = z (1 - u)^{-1/k}, \qquad
        x^{1/k} = S (1 - u)^{-(k-1)/k}, \quad S = x z^{k-1}.

    The binomial series `(1-u)^{-b} = \sum_j c_j u^j` truncated after
    `u^{r-1}` has a tail below `2 c_r |u|^r` for `|u| \le 1/2` (its
    coefficients decrease, `b < 1`).  So `z` accurate to a fraction
    `1/r` of the precision, obtained by a recursive call whose radius is
    then discarded (only the size of the final `u` matters, and that is
    computed rigorously), gives the root to the full precision from one
    evaluation of `S`, `u` and the series.  The recursion starts from a
    double (through the logarithm, as the reduced operand can reach
    `2^k`).

    The series is summed over a common denominator `D`:

    .. math::

        \mathrm{base} \cdot \Bigl(1 + \frac{W}{D}\Bigr), \qquad
        W = \sum_{1 \le j < r} (D c_j) \, u^j,

    with `D` the lcm of the denominators of the `c_j` in lowest terms
    when the coefficients are words (a power of two for the square
    roots, so the division is a shift), and `k^{r-1} (r-1)!` otherwise.
    The coefficient tables are built once per root and cached per
    thread.  The powers of `u` come by squaring, each at the precision
    it contributes at; the division acts on `W` (a fraction `1 - 1/r` of
    the precision), and one full product applies the base, with the tail
    added to the radius.

    The order `r` of the steps is 3 for `k < 8` and 4 above; the
    third function takes it as an argument.

.. function:: void mp_real_agm(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n)
              void _mp_real_agm_order(mp_real_t res, const mp_real_t x, const mp_real_t y, slong n, int m)

    Sets *res* to the arithmetic-geometric mean of `x, y \ge 0` to about
    `n` limbs.  The iteration `a' = (a+b)/2`, `b' = \sqrt{ab}` runs in
    ball arithmetic until `z = (a-b)/(a+b)` satisfies
    `z^2 < 2^{-p/m}`, and is finished by the series

    .. math::

        \operatorname{agm}(a, b) = \frac{(a+b)/2}{{}_2F_1(1/2, 1/2; 1; z^2)},
        \qquad
        \frac{1}{{}_2F_1(1/2, 1/2; 1; x)} = 1 - \sum_{j \ge 1} c_j x^j.

    The coefficients `c_j` are dyadic, positive (Kaluza's theorem: the
    coefficients of `{}_2F_1` are log-convex) and sum to 1 (the
    reciprocal vanishes at `x = 1`), so the truncation after `x^{m-1}`
    has a tail below `x^m` at every order `m`.  The terms are integers
    over one power of two, and the powers of `x` come by squaring at
    the precision each contributes at.  The order `m` is at most 16
    (10 on 32-bit machines); the default depends on the precision, and
    the second function takes it as an argument.

    When `a - b` no longer determines its own size, or `a + b` or `ab`
    is too wide for the division and the square root, the result is the
    enclosure `(a+b)/2 \pm |a-b|/2` (valid since
    `b \le \operatorname{agm}(a, b) \le a`).

    The working set is four numbers of `n` limbs (the finish reuses the
    storage of two of them, and frees each power of `x` after its last
    use); an input aliased by *res* is taken over rather than copied.

Conversions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: void mp_real_get_arb(arb_t res, const mp_real_t x)

    Sets *res* to the ball *x*, losslessly.

.. function:: int mp_real_unique_integer(mp_real_t res, const mp_real_t x)

    If the ball *x* contains exactly one integer, sets *res* to that
    integer (exactly) and returns 1; otherwise returns 0 and leaves
    *res* unchanged.  This is :func:`arb_get_unique_fmpz` with the result
    kept as an :type:`mp_real_t`, free of the size limits of
    :type:`fmpz_t` (the result of a huge computation such as a value of
    the partition function can exceed the `2^{31}` limbs of an ``mpz``).
    With a fractional part of `f` limbs and a radius below one ulp, the
    candidates follow from the carry of the fraction plus the radius and
    the borrow of the fraction less the radius, which are decided by
    scanning the fraction from its top limb until a limb is neither all
    zeros nor all ones; the integer part is then moved down (and
    incremented if needed), in place when *res* aliases *x*, without a
    reallocation.

.. function:: void mp_real_set_dfloat(mp_real_t res, const double * x, slong n, double rad)
              void mp_real_get_dfloat(double * res, double * rad, slong n, const mp_real_t x)

    Conversions to and from the double expansions of the dfloat module,
    given as arrays of `n \le 8` doubles (the midpoint components) and a
    double radius, without going through :type:`arb_t` (about half the
    cost).  :func:`mp_real_set_dfloat` sets *res* to the ball
    `\sum_i x_i \pm \mathrm{rad}`: the sum is exact (accumulated in a
    two's complement frame spanning the components' exponents) and the
    radius is installed as an ulp count, rounded up when it does not fit
    exactly.  Nonfinite components or radius throw (mp_real has no whole
    line).  :func:`mp_real_get_dfloat` rounds greedily, each component the
    double nearest to the remainder (so `|x_{i+1}| \le
    \operatorname{ulp}(x_i)/2`), and sets *rad* (if not ``NULL``) to an
    upper bound for the remainder, the truncated tail of the mantissa and
    the radius of *x*.  Components that would be subnormal are not
    produced (their value goes to the radius), and a value beyond the
    double range gives dfloat's whole line: zero components and an
    infinite radius.

.. function:: void mp_real_print(const mp_real_t x)

    Prints the fields of *x* and its value.

.. function:: void _mp_real_get_fixed(nn_ptr y, ulong * err, const mp_real_t x, slong n)

    Writes a ball known to lie in `[0, 1)` as the fixed-point number
    `(y, n)`, truncating, and sets *err* (if not ``NULL``) to a rigorous
    bound on the error in ulps `B^{-n}`, rounded up; ``UWORD_MAX`` means
    that no bound below `B - 1` ulps is known.  A negative midpoint
    within the radius is clamped to zero and absorbed into the bound.

.. function:: int _mp_real_get_fixed_floor(nn_ptr y, slong n, const mp_real_t x)

    Writes `y = \lfloor x B^n \rfloor` for a ball known to lie in
    `[0, 1)` and returns 1 if the radius determines that floor uniquely;
    returns 0 otherwise (the caller retries at a higher precision).

Implementation notes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Radius.**  Rounding a bound up to a whole count at the bottom limb
costs up to one ulp per operation.  The mantissa deliberately retains at
most one limb of noise, since interval radii over correlated errors grow
faster than true errors.  The bottom limb of an inexact ball is its
noise floor: no ball holds limbs below it, so a sum takes nothing below
an inexact operand's floor (the discarded tail of the other operand
costs one unit there).

**Bounds.**  The bounds are composed in limb arithmetic.  Internally a
bound is a two-limb count at a limb anchor (``umul_ppmm`` of a radius by
a magnitude, ``add_ssaaaa`` of two counts); a sum aligns the lower anchor
to the higher by dividing its count by `B`, rounding up.  The count is
reduced to a single limb, rounding up and moving the anchor, when it is
installed on a ball.  The divisions and roots bound through doubles,
converted to a count on installation.  Magnitudes are bounded through
the leading `\mathrm{FLINT\_BITS}` bits of a mantissa (across a limb
boundary): bounding through the top limb alone would overstate by a
factor 2 for a top limb of 1, the working range of any iteration near 1,
and with the radius held as a whole count of ulps such pessimism turns
directly into truncated limbs.

**Products.**  Products truncate their operands to `n + 2` limbs
beforehand; the dropped tail, below one unit of the lowest kept limb,
enters the radius through the other operand's magnitude.  A long exact
operand thus costs no more than the precision asks.  Products and
quotients of inexact operands keep size + 2 limbs.  Short products are
formed in full, longer ones by windowed middle products
(:func:`flint_mpn_mulmid`); a square always goes through
``flint_mpn_sqr`` in full.

**Sums.**  Sums and differences are formed in the destination buffer
without a temporary.  A difference is taken from the larger operand,
decided up front from the exponents and top limbs, so that no negation
pass follows.  An in-place sum ``x := x + y`` with `y` inside `x`'s
window costs the ``mpn_add`` of `y`'s limbs plus a constant, whatever
the size of `x`: the carry or borrow stops when it runs out, and the
normalization and the installation of the radius are constant time.
A `y` reaching below an inexact `x`'s bottom limb is cut there, so `x`
never moves.  Only a `y` below an exact `x` (the contiguous mantissa
must move) or an `x` longer than the precision asks (a truncation)
cost a pass over `x`.

``mp_real_addmul_ui`` and ``mp_real_submul_ui`` work in place when
`y c` lands inside `x`'s window (and, for a difference, cannot flip its
sign); `y` may reach below an inexact `x`, its tail entering through the
high word of its top dropped limb times `c`.  Other shapes form `y c`
and add it.

**Division by a word.**  ``mp_real_div_ui`` divides by
``mpn_divrem_1`` over an appended zero limb (for an exact dividend,
enough zero limbs for a quotient of `n + 2` limbs), keeping the
quotient's own precision; the inherited radius and the truncation are
bounded one limb below the dividend's anchor.  A chain of divisions by
a common denominator therefore loses no accuracy relative to the
quotients.

**Bounded residuals.**  ``_mp_real_submul_bounded`` reads the sign of
`a - b c` off the top limb of the product window; the bound passed by
the caller must be rigorous.

**Complex products.**  ``mp_real_mul_complex`` forms
`(a_r + i a_i)(b_r + i b_i)` by one call to
:func:`flint_mpn_mul_complex` over the four parts (its high half only,
by :func:`flint_mpn_mulhigh_n_complex`, for long operands).  Both parts
are kept down to the same limb: an accumulation of `\prod \exp(i x_k)`
needs its small imaginary part to the absolute precision of the frame.
Each operand's parts are likewise windowed to a common bottom limb.

**Divisions and roots** require operand balls of relative radius below
`2^{-30}` (checked: a wider divisor is a usage error whose quotient
would be pure noise).  ``mp_real_rsqrt`` and ``mp_real_sqrt`` align the
exponent, call :func:`_mp_real_rsqrt_newton` resp.
:func:`_mp_real_sqrt_newton` (or ``flint_mpn_sqrtrem`` below its cutoff)
and bound the result through the derivative.

Constants
-------------------------------------------------------------------------------

.. function:: void mp_real_const_pi4(mp_real_t res, slong n, int cache)
              void mp_real_const_log2(mp_real_t res, slong n, int cache)
              void mp_real_const_euler(mp_real_t res, slong n, int cache)
              void mp_real_const_e(mp_real_t res, slong n, int cache)
              void mp_real_const_log10(mp_real_t res, slong n, int cache)
              void mp_real_const_catalan(mp_real_t res, slong n, int cache)
              void mp_real_const_zeta3(mp_real_t res, slong n, int cache)
              void mp_real_const_zeta5(mp_real_t res, slong n, int cache)
              void mp_real_const_gamma_1_3(mp_real_t res, slong n, int cache)
              void mp_real_const_gamma_1_4(mp_real_t res, slong n, int cache)
              void mp_real_const_2_div_pi(mp_real_t res, slong n, int cache)
              void mp_real_const_inv_log2_frac(mp_real_t res, slong n, int cache)
              void _mp_real_const_pi4(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_log2(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_euler(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_e(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_log10(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_catalan(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_zeta3(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_zeta5(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_gamma_1_3(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_gamma_1_4(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_2_div_pi(nn_ptr res, ulong * err, slong n, int cache)
              void _mp_real_const_inv_log2_frac(nn_ptr res, ulong * err, slong n, int cache)

    The constants `c = \pi/4`, `\log 2`, Euler's constant `\gamma`,
    `e`, `\log 10`, Catalan's constant `G`, `\zeta(3)`, `\zeta(5)`,
    `\Gamma(1/3)`, `\Gamma(1/4)`, `2/\pi` (for argument reduction
    by `\pi/2`) and `1/\log 2 - 1` (for the quotient of the reduction
    by `\log 2`, whose top limbs :func:`mp_real_exp_bits` and the
    ``nfloat`` exponential read from the static table).

    The ball versions set *res* to a ball for `c` accurate to about `n`
    limbs.  The limb versions set *res* to EXACTLY
    `\lfloor c B^n \rfloor` -- `n` limbs for the constants below 1
    (`\pi/4`, `\log 2`, `\gamma`, `G`, `2/\pi`, `1/\log 2 - 1`), `n + 1` limbs with a
    units limb for the others -- and set *err* (if not ``NULL``) to 1, a
    bound on `c - \lfloor c B^n \rfloor B^{-n}` in ulps.  A floor is
    computed from a ball at a few guard limbs, retried with more guard
    limbs until the rigorous radius determines it uniquely.

    With *cache* = 0 the value is computed from scratch and not kept.
    With *cache* = 1 it comes from a shared cache (one per constant,
    read by every thread; see `Shared tables`_) holding the floor of
    the constant at the largest precision computed so far.  Floors
    nest (the top limbs of `\lfloor c B^N \rfloor` are
    `\lfloor c B^n \rfloor` for `n \le N`), so a shorter request is a
    copy.  A longer one recomputes the entry at `\max(n + 5, 1.5 N)`
    limbs, so that a sequence of growing requests costs a constant
    factor over the last one.  The ball version with *cache* = 1
    encloses the cached floor with a radius of one ulp.  `\pi/4`,
    `\log 2`, `2/\pi` and `1/\log 2 - 1` come from static 64-limb
    tables (generated by ``dev/gen_mp_real_elem_tables.py``) up to 64
    limbs, with *cache* = 1 as well as internally (where they are read
    in place without copying), and from the cache beyond.

.. function:: void _mp_real_const_clear_cache(void)

    Frees the shared caches of the constants; no other thread may be
    using them.  The caches are also freed by
    :func:`flint_cleanup_master`.

.. type:: mp_real_const_func

.. function:: void _mp_real_const_arb(arb_t res, mp_real_const_func f, slong prec)

    Sets *res* to the constant computed by the ball function *f* (with
    *cache* = 0), rounded to *prec* bits: the entry point of the arb
    wrappers (:func:`arb_const_pi`, :func:`arb_const_log2`, ...), which
    keep their own caches.

.. function:: void _mp_real_const_euler_tune(mp_real_t res, slong n, int set)

    Euler's constant with a forced choice of the logarithms of `\log m`
    (see below): *set* = 0, 1, 2, 3 selects the sets `\{2\}`,
    `\{2, 3\}`, `\{2, 3, 5\}`, `\{2, 3, 5, 7\}`, and `-1` the automatic
    choice.

Algorithms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Binary splitting.**  All the series below are summed by binary
splitting in ``mp_real_t`` ball arithmetic.  Blocks of terms at the
leaves accumulate iteratively over exact mpn integers with a backward
recurrence; the tree keeps exact integers until they outgrow the target
precision, and truncated balls afterwards.

**Pi and log 2** are hypergeometric series summed by
:func:`mp_real_hypgeom_series`: `\pi` by the Chudnovsky series (see the
example under :func:`mp_real_hypgeom_series_int64`), and `\log 2` by the
series of [Zun2025]_ with 11.9 bits per term,

.. math::

    \log 2 = \frac{1}{2160} \sum_{k \ge 0} (1497 + 1794 k)
        \prod_{j=1}^{k} \frac{j (2j - 1)}{216 (6j + 1)(6j + 5)}.

At small precision they are read from static tables or the cache.

**Euler's constant** comes from a static 3456-bit table at low
precision and otherwise from the Brent-McMillan formula ([BM1980]_)
with the error bound of [BJ2013]_:

.. math::

    \gamma = \frac{S_0}{I_0} - \frac{K_0}{I_0^2} - \log m
    + O(e^{-8m}), \quad
    I_0 = \sum_{k \ge 0} \left(\frac{m^k}{k!}\right)^2, \quad
    S_0 = \sum_{k \ge 0} \left(\frac{m^k}{k!}\right)^2 H_k,

with `K_0 = I_0(2m) K_0(2m)` from its asymptotic series and
`8m \ge b \log 2` for `b` bits.  `S_0` and `I_0` come from one binary
splitting over `4.97 m` terms in dual numbers: with `\varepsilon^2 = 0`,

.. math::

    F(\varepsilon) = \sum_{k \ge 0} \frac{m^{2k}}
        {\prod_{j=1}^k (j + \varepsilon)^2},
    \quad I_0 = F(0), \quad S_0 = -\tfrac{1}{2} F'(0).

The splitting of `F` carries a scalar `P`, `D = D_0 + D_1 \varepsilon
= \prod (k + 1 + \varepsilon)` (so that `Q = D^2`) and
`T = T_0 + T_1 \varepsilon`, merged as

.. math::

    T = P_1 T_2 + D_2^2 T_1, \qquad D = D_1 D_2, \qquad P = P_1 P_2,

five full-size and five half-size products per merge.  At the end, with
`X = T_0 + D_0^2`,

.. math::

    \frac{S_0}{I_0} = \frac{2 T_0 D_1 - T_1 D_0}{2 D_0 X}, \qquad
    \frac{K_0}{I_0^2} = \frac{D_0^4 T_K}{X^2 Q_K},

where `K_0` takes a second splitting over `2m` terms at half the
precision.  The main sum runs first and leaves only the quotient and two
half-precision values, so that the peak memory is that of the main sum.

The parameter `m` is rounded up so that `\log m` is a combination of a
few fast logarithms, over the primes `\{2\}`, `\{2, 3\}`, `\{2, 3, 5\}`
or `\{2, 3, 5, 7\}`:

* `\{2, 3\}`: `\log(9/8)` and `\log(256/243)`;
* `\{2, 3, 5\}`: `\log(81/80)`, `\log(32805/32768)`, `\log(25/24)`;
* `\{2, 3, 5, 7\}`: `\log(2401/2400)`, `\log(4375/4374)`,
  `\log(225/224)`, `\log(64/63)`,

each by Zuniga's series (see `Machin terms as Zuniga series`_).  The
splittings cost about `m (61 + \log_2 o)` for `o` the odd part of `m`
(the power of two in `m` being free), so an `m` a little larger with a
small odd part can be faster.  Every `m` in `[m_0, 1.5 m_0]` that is
smooth over `\{2, 3, 5, 7\}` is scored by this cost plus that of its
logarithms, and the best is taken.

**Other constants.**  Except for `\log 10`, each is a hypergeometric
series through :func:`mp_real_hypgeom_series`, taken from the formula
files of y-cruncher (https://github.com/Mysticial/y-cruncher-Formulas);
of the files available for each constant, the one measured fastest in
this binary splitting was chosen.

* `e = \sum_{k \ge 0} 1/k!`.
* `\log 10 = \log 2 + \log 5`, from the three-prime Machin-type set of
  :func:`_mp_real_log_primes_vec` (three atanh terms by Zuniga's
  series, on up to three threads).
* `G`: Pilehrood's short series (2010).
* `\zeta(3)`: Zuniga's 2023-vi series (2.05 bits per term).
* `\zeta(5)`: Zhi-Wei Sun's identity (2025),
  `\zeta(5) = (3S + 56 \pi^2 \zeta(3))/540` with `S` a
  `{}_3F_2`-type series.
* `\Gamma(1/3) = (810^{1/4} \pi / X)^{1/3}` with `X` Guillera's 2023
  series.
* `\Gamma(1/4) = (\pi^6 / (322 S^4))^{1/8}` with `S` Ebisu's 2016
  lemniscate series; the eighth root is three square roots.  From
  `2^{16}` limbs on, the AGM formula

  .. math::

      \Gamma(1/4) = \sqrt{\frac{(2\pi)^{3/2}}{\operatorname{agm}(1, \sqrt 2)}}

  is used instead, with the AGM computed in parallel with `\pi`: it
  costs `O(M(n) \log n)` against the `O(M(n) \log^2 n)` of the series.

The independent parts of `\zeta(5)` and of the two `\Gamma` values run
on separate threads.

Elementary functions
-------------------------------------------------------------------------------

.. function:: void mp_real_sin_cos_bits(mp_real_t res1, mp_real_t res2, const mp_real_t x, slong prec)

    Sets *res1* and *res2* to balls containing `\sin x` and `\cos x`
    for every `x` in the ball *x*, to a relative accuracy of about
    `2^{-\mathrm{prec}}`.  Either output may be ``NULL`` and either may
    alias *x*; the numbers of working and output limbs are chosen
    internally.

    If *x* is inexact, *prec* is first lowered to a few bits beyond the
    accuracy its radius allows (relative for `|x| < 1`, where the sine
    keeps the relative accuracy of `x`, absolute otherwise); the
    midpoint is evaluated exactly as given and the radius is added to
    both outputs at the end (`|\sin'|, |\cos'| \le 1`).  A ball
    within `(-1, 1)` whose radius exceeds a quarter of its midpoint
    gives `[0 \pm 2^e]` and `1 \pm 2^{2e-1}` for `|x| < 2^e`; a
    wider ball elsewhere gives `[0 \pm 1]`, as do arguments with
    `|x| \ge 2^{\max(65536, 4 \mathrm{prec})}`.

    Evaluation goes by the number of leading zero bits `z` of the
    midpoint (`|x| < 2^{-z}`):

    * `2z \ge \mathrm{prec} + 3`: `\sin x = x`, `\cos x = 1`, with the
      Taylor remainders in the radius.
    * `4z \ge \mathrm{prec} + 3`: `x - x^3/6` and `1 - x^2/2`.
    * Otherwise, for `|x| < 1`: the fixed-point kernel on `x` at
      `\mathrm{prec} + z` bits plus guard bits, so that the sine keeps
      its relative accuracy.  From a size-dependent number of leading
      zero bits on (from about 20 at a few limbs to thousands at 65536
      limbs), the argument goes to a series directly instead, whose
      cost falls with the argument's size: the sine and `1 - \cos` of
      :func:`_mp_real_series_tapered` resp. :func:`_mp_real_series_rs`
      below 32 zero bits, :func:`_mp_real_sin_cos_reduced` beyond.
    * `|x| \ge 1`: reduction mod `\pi/2` with a quotient `q` that
      approximates `x \cdot 2/\pi` from below (so
      `t = x - q \pi/2 \ge 0` and no correction is needed), a fold
      `v = \pi/2 - t` when `t > \pi/4`, and the kernel on `v`, with the
      outputs swapped and negated according to `q \bmod 4`.  The
      reduction adds at most 3 ulps to the kernel's bound.  Near a zero
      of the sine or cosine the loss of significance shows in the
      output's radius; the reduction is not repeated at a higher
      precision.

    The kernel is the dense-table kernel at 4 to 32 limbs on 64-bit
    machines (`Dense-table kernels at four to ten limbs`_, `Dense-table
    kernels with run-time tables, eleven to thirty-two limbs`_),
    :func:`_mp_real_sin_cos_bitwise_rs` otherwise up to 600 limbs,
    :func:`_mp_real_sin_cos_diophantine` (32 primes) up to 65536 limbs
    and :func:`_mp_real_sin_cos_notab` above.

.. function:: int mp_real_tan_bits(mp_real_t res, const mp_real_t x, slong prec)

    Sets *res* to a ball containing `\tan y` for every `y` in the ball
    *x*, to a relative accuracy of about `2^{-\mathrm{prec}}`, and returns
    1; or returns 0 with *res* set to zero when the ball may contain a
    pole (or `|x| \ge 2^{\max(65536, 4 \mathrm{prec})}`, where the
    reduction is not attempted).  *res* may alias *x*.

    The midpoint `m` is evaluated exactly as given, as by
    :func:`mp_real_sin_cos_bits`: for `|m| < 1` as `m` itself when
    `2z \ge \mathrm{prec} + 3` for its leading zero bits `z` (beyond,
    the kernel's series for small arguments is faster than a Taylor
    polynomial in ball arithmetic), else by the tan kernel on `|m|` at
    `z` bits more precision; otherwise after the reduction
    `|m| = a \pi/2 + s v`, `v \in [0, \pi/4]`, as `\tan m = s \tan v`
    (`a` even) or `-s / \tan v` (`a` odd, by one Newton reciprocal of the
    kernel's fixed-point output, :func:`_mp_real_inv_newton`).  The
    reduction is absolute, and near a zero or pole of the tangent `v` is
    small: an absolute error is then a relative one of the result, so
    the kernel runs at the leading zero bits of `v` more precision, the
    reduction being repeated at those limbs when the first pass had
    fewer.  An exact `m` is never a
    multiple of `\pi/2`, and by the irrationality measure of `\pi` a
    `b`-bit `m` lies about `2^{-7.2 b}` away from one at least; the
    repetition is capped at `8b + 1024` more bits (status 0 beyond).
    For an inexact `x` no repetition is made once `v` falls below the
    radius.

    The kernel on `v \in [0, 1)` is the tangent series
    (:func:`_mp_real_series_rs_tan`) where its tables reach, else from a
    size-dependent number of leading zero bits the reduced sine and
    `1 - \cos` of :func:`_mp_real_sin_cos_reduced` and a division, else
    on 64-bit machines up to 22 limbs the sine and cosine kernels and a
    division (:func:`_mp_real_tan_sc`), else
    :func:`_mp_real_tan_bitwise_rs` up to 600 limbs,
    :func:`_mp_real_tan_diophantine` up to 65536 limbs and the sine and
    cosine of :func:`_mp_real_sin_cos_notab` and a division above.  It
    costs about as much as :func:`mp_real_sin_cos_bits` at every size.

    The radius `r` of `x` enters through a pole test: with `T`
    enclosing `\tan m` and `U = r (1 + |T|) \ge r / |\cos m|`, `U < 1`
    means that `|\cos|` (which is 1-Lipschitz) stays above
    `|\cos m| - r > 0` over the ball, which then contains no pole, and

    .. math::

        |\tan y - \tan m| \le \frac{r}{(|\cos m| - r)^2}
                          \le \frac{r (1 + T^2)}{(1 - U)^2},

    evaluated in two-limb balls.  Otherwise the status is 0; the test is
    conservative by a factor below `\sqrt 2`, where the output would
    have almost no accuracy anyway.  *prec* is first lowered to a few
    bits beyond the accuracy the radius allows, as for
    :func:`mp_real_sin_cos_bits`.

.. function:: void mp_real_sin_cos_pi_bits(mp_real_t res1, mp_real_t res2, const mp_real_t x, slong prec)
              int mp_real_tan_pi_bits(mp_real_t res, const mp_real_t x, slong prec)

    The same as :func:`mp_real_sin_cos_bits` and :func:`mp_real_tan_bits`
    for `\sin(\pi y)`, `\cos(\pi y)` and `\tan(\pi y)`, with an exact
    argument reduction: the midpoint is dyadic, and its units bit and
    fraction bits give

    .. math::

        |m| = \frac{a}{2} + s t, \qquad a \in \{0, 1, 2, 3\}, \;
        s = \pm 1, \; t \in [0, 1/4]

    exactly, after which `v = \pi t` enters the same formulas as after the
    reduction mod `\pi/2`.  `t = 0` gives the exact values `0`, `\pm 1`
    (and a pole of the tangent for odd `a`, status 0), `t = 1/4` the
    exact `\tan(\pi/4) = 1`.  Otherwise `t` has `z` leading zero bits,
    and `v = (\pi/4) \cdot 4t` goes to the kernel at `z` more bits of
    absolute precision, or for `2z \ge \mathrm{prec} + 3` to the
    first-order values `\sin v = \tan v = v`, `\cos v = 1` with the
    Taylor remainders in the radius: the outputs keep their relative
    accuracy near every zero and pole, at any magnitude of the argument
    (whose reduction costs nothing).  `t` is kept in a buffer the size of
    the argument's mantissa and the results are written into the outputs
    directly, with no temporary balls.  The radius `r` of `x` is a radius `\pi r` of the argument:
    added to the sine and cosine, and through the pole test for the
    tangent.  A radius `\pi r \ge 1/4` gives `[0 \pm 1]` for the sine and
    cosine.

.. function:: void mp_real_sin_cos_pi_ui_div_ui(mp_real_t res1, mp_real_t res2, ulong p, ulong q, slong n)

    Sets *res1* and *res2* to balls containing `\sin(\pi p/q)` and
    `\cos(\pi p/q)` for words `p \ge 0` and `q \ge 1` (throws for
    `q = 0`), to a relative accuracy of about `B^{-n}`.  Either output
    may be ``NULL``; the two must be distinct.  The rational values
    `0, \pm 1, \pm 1/2` are exact.

    The fraction is reduced to lowest terms and folded by the symmetries
    to `\pi p/q = \pi a \pm \pi t/q` with `0 \le 2t \le q`; the
    nonzero values of `\sin(\pi t/q) \ge 2t/q` and
    `\cos(\pi t/q) \ge (q - 2t)/q` are at least `1/q`, so
    :func:`_mp_real_sin_cos_pi_ui_div_ui` at
    `n + \lceil (\operatorname{bits}(q) + 16) / \mathrm{FLINT\_BITS} \rceil`
    limbs gives the relative accuracy, its bound becoming the radius.
    The values `0, \pm 1, \pm 1/2` are set directly (only the other
    output, if any, is computed).

.. function:: int mp_real_tan_pi_ui_div_ui(mp_real_t res, ulong p, ulong q, slong n)

    Sets *res* to a ball containing `\tan(\pi p/q)` for words `p \ge 0`
    and `q \ge 1` (throws for `q = 0`), to a relative accuracy of about
    `B^{-n}`, and returns 1; or returns 0 with *res* set to zero at a pole
    (`p/q` a half-integer).  The values `0` and `\pm 1` are exact.

    With the fraction in lowest terms, `\tan(\pi p/q) = \pm \tan(\pi t/q)`,
    `0 < 2t < q`, is computed at the working precision of
    :func:`mp_real_sin_cos_pi_ui_div_ui`: by the closed forms
    `\sqrt 3`, `1/\sqrt 3`, `\sqrt 2 \mp 1`, `2 \mp \sqrt 3` for
    `q = 3, 6, 8, 12` and `\sqrt{5 \mp 2\sqrt 5}`, `\sqrt{1 \mp 2/\sqrt 5}`
    for `q = 5, 10` (from 12 limbs); by the Chebyshev sine and cosine and
    one division for the denominators where they beat the kernel for sine
    and cosine; else by the tan kernel on `\pi t/q`, or on its complement
    `\pi/2 - \pi t/q` for `4t > q` followed by one reciprocal.

.. function:: void mp_real_exp_bits(mp_real_t res, const mp_real_t x, slong prec)
              int mp_real_log_bits(mp_real_t res, const mp_real_t x, slong prec)
              void mp_real_atan_bits(mp_real_t res, const mp_real_t x, slong prec)

    Set *res* to a ball containing `\exp(y)`, `\log(y)` resp.
    `\operatorname{atan}(y)` for every `y` in the ball *x*, to a
    relative accuracy of about `2^{-\mathrm{prec}}`, in the manner of
    :func:`mp_real_sin_cos_bits`: the midpoint `m` is evaluated
    exactly as given, *prec* is first lowered to a few bits beyond what
    the radius `\rho` allows, the numbers of working and output limbs
    are chosen internally, and *res* may alias *x*.  The results keep
    their relative accuracy near `\exp(y) = 1`, `\log(y) = 0` and
    `\operatorname{atan}(y) = 0`.

    The radius enters as the factor `[1 \pm (e^{\rho} - 1)]` (exp), as
    `\rho / (m - \rho)` (log; the mean value theorem) and as `\rho`, or
    `\rho / (|m| - \rho)^2` when the ball stays beyond 1 in absolute
    value (atan).  A ball too wide for two bits gives `(0, e^{m + \rho}]`
    enclosed as `[0 \pm e^{m + \rho}]` for exp.

    :func:`mp_real_exp_bits` throws for `|y| \ge 2^{\mathrm{FLINT\_BITS}
    - 5}`, where the result's exponent would leave the safe range.
    :func:`mp_real_log_bits` returns 1, or 0 with *res* set to zero when
    the ball is not strictly positive (decided exactly, by the midpoint's
    mantissa against the radius count).

    In all three functions, arguments with many leading zero bits (`m`
    for exp and atan, `1/|m|` for large atan arguments, `|m - 1|` for
    log) go to a reduced series from a size-dependent threshold on,
    since its cost falls with the argument's size where that of the
    table-based kernels does not.  The kernels are the bitwise ones up
    to 600 limbs; beyond, exp uses :func:`_mp_real_exp_diophantine` up
    to 65536 limbs and :func:`_mp_real_exp_notab` above, and log and
    atan use the Newton-Taylor functions :func:`_mp_real_neglog_newton`
    and :func:`_mp_real_atan_newton`.

    **exp.**  For small `|m|`: `1 + m`, `+ m^2/2`, `+ m^3/6` with the
    next term in the radius; the series of :func:`_mp_real_exp_reduced`
    resp. the alternating series of `\exp(-|m|)` by
    :func:`_mp_real_series_rs`; `m \in (0, 1)` the kernel directly.
    Otherwise `|m| = q \log 2 + t`, `t \in [0, \log 2)`, with `q` from
    two limbs of `1/\log 2` (at most one too low, then corrected) and
    one multiply-subtract by `\log 2`; for `m < 0` the quotient is
    rounded up and the fraction complemented, so that
    `\exp(m) = 2^{-q} \exp(t)` directly.

    **log.**  With `m = 2^E u`, `u \in [1/2, 1)`: near 1, `D`,
    `D - D^2/2`, `+ D^3/3` for `D = m - 1`, the rest in the radius;
    then `2 \operatorname{atanh}(s / (2 \pm s))`, `s = |m - 1|`, by one
    fixed-point division and a reduced series; otherwise

    .. math::

        \log m = (E - 1) \log 2 + \operatorname{log1p}(2u - 1)
        \quad\text{resp.}\quad
        \log m = E \log 2 - (-\log u),

    with the log1p kernel (dense-table at 4 to 32 limbs on 64-bit
    machines, bitwise otherwise) resp. :func:`_mp_real_neglog_newton`,
    combined in fixed point against the floor of `\log 2`.  Near 1 each
    form cancels up to `z + 1` bits, which the working precision
    absorbs.

    **atan.**  `|m| < 1`: `m`, `m - m^3/3` with the next term in the
    radius, then a reduced series, then the kernel at
    `\mathrm{prec} + z` bits.  For `|m| \ge 1`:

    .. math::

        \operatorname{atan}(1) = \frac{\pi}{4}, \qquad
        \operatorname{atan}(m) = \frac{\pi}{4} +
            \operatorname{atan}\Bigl(\frac{s}{2 + s}\Bigr),
            \; s = m - 1 \;\; (1 < m < 2), \qquad
        \operatorname{atan}(m) = \frac{\pi}{2} -
            \operatorname{atan}\Bigl(\frac{1}{m}\Bigr) \;\; (m \ge 2),

    the ratio by one ``flint_mpn_divapprox_fraction`` and the reciprocal
    by one ``flint_mpn_invapprox``.  These results exceed `\pi/4`, so
    absolute accuracy suffices.

    At few limbs the four functions call the per-size kernels
    (the table-driven kernels of `Table-driven kernels at one to three
    limbs`_ at `n \le 3`, ``*_opt_<n>.c`` otherwise) directly, take their guard bits from the per-size
    reduction parameters (``MP_REAL_*_OPT_R`` in ``impl.h``) and read
    the argument's limbs in place.

    The radius bounds use double arithmetic that never leaves the
    normal range: scalings by powers of `B` go through
    ``d_mul_2exp_inrange`` with clamped exponents, so no call raises
    ``FE_UNDERFLOW``.

.. function:: int mp_real_log1p_bits(mp_real_t res, const mp_real_t x, slong prec)
              void mp_real_expm1_bits(mp_real_t res, const mp_real_t x, slong prec)
              void mp_real_sinh_cosh_bits(mp_real_t res1, mp_real_t res2, const mp_real_t x, slong prec)
              void mp_real_tanh_bits(mp_real_t res, const mp_real_t x, slong prec)
              int mp_real_asin_bits(mp_real_t res, const mp_real_t x, slong prec)
              int mp_real_acos_bits(mp_real_t res, const mp_real_t x, slong prec)
              void mp_real_asinh_bits(mp_real_t res, const mp_real_t x, slong prec)
              int mp_real_acosh_bits(mp_real_t res, const mp_real_t x, slong prec)
              int mp_real_atanh_bits(mp_real_t res, const mp_real_t x, slong prec)

    Set *res* (*res1* to `\sinh`, *res2* to `\cosh`, either may be
    ``NULL``) to a ball containing `\operatorname{log1p}(y) = \log(1+y)`,
    `\operatorname{expm1}(y) = e^y - 1`, `\sinh(y)`, `\cosh(y)`,
    `\tanh(y)`, `\operatorname{asin}(y)`, `\operatorname{acos}(y)`,
    `\operatorname{asinh}(y)`, `\operatorname{acosh}(y)` resp.
    `\operatorname{atanh}(y)` for every `y` in the ball *x*, to a
    relative accuracy of about `2^{-\mathrm{prec}}`, or less as the
    radius allows; the outputs may alias *x*.

    The functions are compositions of :func:`mp_real_exp_bits`,
    :func:`mp_real_log_bits` and :func:`mp_real_atan_bits` by formulas
    without cancellation, so that the results keep their relative
    accuracy at the zeros of the functions:
    `\operatorname{log1p}(y) = \log(1 + y)` with `1 + y` formed exactly
    (or to `\mathrm{prec} + z` bits for `|y| < 2^{-z}`);
    `\operatorname{expm1}` from `\exp` at `\mathrm{prec} + z` bits;
    `\sinh` and `\cosh` from `\operatorname{expm1}(\pm m)`;
    `\tanh(m) = \operatorname{expm1}(2m) / (\operatorname{expm1}(2m) + 2)`;
    `\operatorname{asin}(m) = \operatorname{atan}(m / \sqrt{(1-m)(1+m)})`
    and `\operatorname{acos}(m) = 2 \operatorname{atan}(\sqrt{(1-m)/(1+m)})`
    with the differences `1 \pm m` exact;
    `\operatorname{asinh}(m) = \operatorname{log1p}(|m| + m^2 / (1 + \sqrt{1 + m^2}))`;
    `\operatorname{acosh}(m) = \operatorname{log1p}(d + \sqrt{d (m + 1)})`
    with `d = m - 1` exact; and
    `\operatorname{atanh}(m) = \operatorname{log1p}(2|m| / (1 - |m|))/2`.
    Tiny and huge arguments are handled by the leading terms with a
    bounded remainder (`m`, `\log(2|m|)`, `\pm 1`), without forming
    `m^2` (whose exponent could leave the safe range); the values
    `\operatorname{asin}(\pm 1)`, `\operatorname{acos}(\pm 1)` and
    `\operatorname{acosh}(1)` are exact.

    The midpoint `m` is evaluated as given, and the radius `\rho` enters
    through a bound for the variation over the ball: for example
    `\rho \cosh(|m| + \rho)` for `\sinh`, `4 \rho e^{-2(|m| - \rho)}`
    for `\tanh`, `\rho / (|m| - \rho)` for large asinh arguments,
    `\rho / \sqrt{1 - |m| - \rho}` for asin and acos (or
    `2 \sqrt{1 - |m| + \rho}` for a ball reaching close to `\pm 1`), and
    `\rho / (1 - |m| - \rho)` for atanh.

    The functions returning ``int`` return 1, or 0 with *res* set to
    zero when the ball is not certified to lie in the domain: `1 + y > 0`
    for log1p, `|y| < 1` for atanh, and for asin, acos and acosh
    `|y| \le 1` resp. `y \ge 1` for an exact argument and strictly
    inside for an inexact one.  :func:`mp_real_expm1_bits` and
    :func:`mp_real_sinh_cosh_bits` throw like :func:`mp_real_exp_bits`
    for `y \ge 2^{\mathrm{FLINT\_BITS} - 5}` resp.
    `|y| \ge 2^{\mathrm{FLINT\_BITS} - 5}` (expm1 of large negative
    arguments gives `-1` with a tiny radius).

.. function:: void mp_real_exp_notab_log2(mp_real_t res, const mp_real_t x, slong n)
              void mp_real_exp_notab_squaring(mp_real_t res, const mp_real_t x, slong n)
              void mp_real_exp_agm(mp_real_t res, const mp_real_t x, slong n)

    Set *res* to a ball containing `\exp(y)` for every `y` in the ball
    *x*, to a relative accuracy of about `B^{-n}`, without building the
    tables of the bitwise and diophantine kernels: intended for a single
    exponential at huge precision, where the tables of
    :func:`mp_real_exp_bits` do not pay for themselves.  The ball
    handling (the radius factor, the lowered precision, the range
    `|y| < 2^{\mathrm{FLINT\_BITS} - 5}`) is that of
    :func:`mp_real_exp_bits`; any midpoint is accepted.

    :func:`mp_real_exp_notab_log2` is :func:`mp_real_exp_bits` with
    :func:`_mp_real_exp_notab` as the kernel: `m = k \log 2 + t`, so it
    computes (and caches) `\log 2` at the working precision; the fastest
    choice when `\log 2` is already cached.

    :func:`mp_real_exp_notab_squaring` needs no `\log 2`:
    `\exp(|m|) = \exp(|m| 2^{-s})^{2^s}` with `|m| 2^{-s} < 1/2`,
    :func:`_mp_real_exp_notab` on the reduced argument and `s` ball
    squarings (with `s + 16` guard bits for the doubling relative
    error), one reciprocal for `m < 0`.  It is 1.2-1.7 times faster than
    the reduction when `\log 2` must be computed.

    :func:`mp_real_exp_agm` takes one Newton-Taylor step on the AGM
    logarithm: `y_0 = \exp(m)` to about `P/(N+1)` limbs (`N = 8` Taylor
    terms, 12 from 20000 limbs) by the squarings, then with
    `s = y_0 2^j \ge 2^{p/2 + 16}`,
    `\log s = \pi / (2 \operatorname{agm}(1, 4/s)) + \varepsilon`,
    `|\varepsilon| \le 64 (\log s + 8)/s^2` (Borwein and Borwein, *Pi and
    the AGM*, Thm. 7.2), `\log y_0 = \log s - j \log 2`, and
    `\exp(m) = y_0 \exp(m - \log y_0)` with the correction by the Taylor
    series of the Newton steps of `-\log` and `\operatorname{atan}`
    (``newton.c``).  When `y_0` is already large (`j = 0`, as for the
    `\exp(C)` of the partition function) no `\log 2` is needed at all.
    Below about 30 limbs it calls the squarings.  It overtakes them at
    about 15000 limbs (1.0-1.2 times faster at 18000-180000 limbs).

The partition function
-------------------------------------------------------------------------------

.. function:: void mp_real_partitions_hrr(mp_real_t res, ulong nhi, ulong nlo)
              void _mp_real_partitions_hrr(mp_real_t res, const fmpz_t n, int flags)

    Sets *res* to the exact value of the number of partitions `p(n)`,
    `n = n_{lo} + n_{hi} 2^{\mathrm{FLINT\_BITS}}`.  The result is an
    exact :type:`mp_real_t`, so its size is not limited by that of an
    :type:`fmpz_t` (an ``mpz`` holds at most `2^{31}` limbs, exceeded from
    `n \approx 1.4 \cdot 10^{21}`).  For `n < 128` it is read from a
    table; for `n < 417` (on 64-bit machines, while `p(n) < 2^{64}`) it is
    computed by the pentagonal recurrence, which costs less than the
    formula below up to about `n = 500` so that no method gluing the
    two pays off; otherwise by the Hardy-Ramanujan-Rademacher formula

    .. math::

        p(n) = \sum_{k=1}^N A_k(n) \sqrt{\frac{3}{k}} \frac{4}{24n-1}
            \left(\cosh z_k - \frac{\sinh z_k}{z_k}\right) + R_N,
        \qquad z_k = \frac{\pi \sqrt{24n-1}}{6k},

    with `A_k(n)` from :func:`arith_hrr_expsum_factored`, `N` and the
    bound for `R_N` (:func:`partitions_rademacher_bound`) as in the
    previous arb implementation, and each term to the precision it
    contributes at.  The terms of at most about 200 bits are evaluated
    in batches of dfloat balls (double to quad-double, with the vector
    exponential, cosine, sine and square roots) when dfloat is
    supported; the others in mp_real balls, with `\exp(z_k)` as the
    `k`-th root of `\exp(z_1)` (computed without tables, by the
    squarings or the AGM Newton-Taylor step of :func:`mp_real_exp_agm`
    at huge precision).  The ball containing `p(n)` is resolved by
    :func:`mp_real_unique_integer`.

    With several threads (:func:`flint_set_num_threads`, from
    `n = 10^8`), the constants are computed once with the whole thread
    budget (the parallelism inside the binary splitting, the exponential
    and the multiplications), and the terms are then summed in ranges of
    `k` of about equal estimated cost taken from a queue, each into its
    own accumulator.  The peak memory is about that of computing `\pi` to
    the same precision (about 24 numbers of the size of `p(n)` on one
    thread at `n = 10^{16}` and `10^{17}`).

    The second function always uses the formula for `n \ge 2` (and gives
    `p(n)` for any *n*); *flags* is a combination of
    ``MP_REAL_PARTITIONS_NO_DFLOAT`` (all terms in mp_real) and
    ``MP_REAL_PARTITIONS_TEST`` (the rarely taken paths at any `n`:
    several tasks and threads, the division by an inverse of `24n - 1`,
    the dfloat square roots of products beyond `2^{53}`), for testing.

Series
-------------------------------------------------------------------------------

.. function:: void _mp_real_atan_frac_bsplit(mp_real_t res, nn_srcptr p, slong pn, nn_srcptr q, slong qn, slong n)
              void _mp_real_atanh_frac_bsplit(mp_real_t res, nn_srcptr p, slong pn, nn_srcptr q, slong qn, slong n)

    Sets *res* to `\operatorname{atan}(p/q)` resp.
    `\operatorname{atanh}(p/q)` for the mpn integers `0 \le p < q`
    (with `p/q` bounded away from 1, at most about 0.99), to about `n`
    limbs, by binary splitting of

    .. math::

        \frac{p}{q} \sum_{k \ge 0} \frac{s^k}{2k+1} \Bigl(\frac{p^2}{q^2}\Bigr)^k,
        \qquad s = -1 \text{ (atan)}, \; s = 1 \text{ (atanh)}.

    Over a range of terms the sum is `N / (D Q^{\mathrm{len}-1})` with
    `D = \prod (2k+1)`, `Q = q^2`.  Leaves of up to 32 terms are exact
    mpn Horner loops.  When `Q` is large compared to the `2k + 1`, the
    tree is balanced (`2^J` leaves of `L` terms), so that the pure
    powers `Q^{L 2^i}` and `P^{L 2^i}` needed by the merges

    .. math::

        N = N_1 D_2 Q^{b-m} + s^{m-a} P^{m-a} N_2 D_1

    can be tabulated by squarings; this saves a full-size product per
    node.  Otherwise the denominator `D Q^{b-a}` is carried through the
    tree.

Hypergeometric series (y-cruncher format)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A generic binary splitting backend for hypergeometric series in the
format of y-cruncher's ``SeriesHypergeometric`` (the formula files in
https://github.com/Mysticial/y-cruncher-Formulas).  The series is

.. math::

    S = \left(\frac{c_Q + c_P \sum_{k \ge 1} \frac{P(k)}{Q(k)}
        \prod_{j=1}^{k-1} \frac{R(j)}{Q(j)}}{c_D}\right)^{\pm 1}

for integer polynomials `P, Q, R` (coefficient `i` of degree `i`,
`Q(k) \ne 0` for `k \ge 1`) and nonzero integers `c_P, c_Q, c_D`
(``CoefficientP``, ``CoefficientQ``, ``CoefficientD`` and ``Power`` in
y-cruncher's files; the prefactors y-cruncher multiplies in, such as an
inverse square root, are left to the caller).  Unlike y-cruncher, the
coefficients may have any size.  The term and the ratios share the
denominator `Q` and the product runs to `k - 1`, so that most formulas
need only small polynomials.  The series must converge geometrically:
`\deg R \le \deg Q`, with `|\mathrm{lc}(R)| < |\mathrm{lc}(Q)|` when
the degrees agree.

.. type:: mp_real_hypgeom_int_struct

    A signed integer `(-1)^{\mathrm{neg}} (d, n)` given as a limb array
    ``d`` of length ``n`` (zero for ``n = 0``) and a sign ``neg``.

.. type:: mp_real_hypgeom_series_struct

    The series: ``power`` (1 or -1), ``coefP``, ``coefQ``, ``coefD``,
    and the polynomials as arrays ``P``, ``Q``, ``R`` of
    ``mp_real_hypgeom_int_struct`` with lengths ``Plen``, ``Qlen``,
    ``Rlen``.

.. function:: void mp_real_hypgeom_series(mp_real_t res, const mp_real_hypgeom_series_struct * s, slong n)

    Sets *res* to a ball containing `S`, with about
    ``FLINT_BITS * (n - 1)`` accurate bits (as elsewhere, `n` includes
    the caller's guard limbs).  Throws if the series is not
    geometrically convergent or its tail cannot be bounded.

.. function:: void mp_real_hypgeom_series_int64(mp_real_t res, int power, int64_t coefP, int64_t coefQ, int64_t coefD, const int64_t * P, slong Plen, const int64_t * Q, slong Qlen, const int64_t * R, slong Rlen, slong n)

    Convenience wrapper taking the coefficients as ``int64_t``, as they
    appear in the y-cruncher files.  For example, the Chudnovsky
    formula gives `\pi` as `S / \sqrt{10005}` with

    * ``power`` `= -1`, `c_P = 1`, `c_Q = 13591409`, `c_D = 4270934400`,
    * `P = [-67957045, -2100495856, 23608573992, -57896553024, 39250089648]`,
    * `Q = [0, 0, 0, -10939058860032000]`,
    * `R = [-5, 46, -108, 72]`,

    and `\operatorname{atan}(p/q)` is `S` with `c_P = c_Q = p`,
    `c_D = q`, `P = R = [p^2, -2p^2]`, `Q = [q^2, 2q^2]` (atanh:
    `P = R = [-p^2, 2p^2]`).

The implementation sums

.. math::

    T = T_1 Q_2 + R_1 T_2, \qquad Q = Q_1 Q_2, \qquad R = R_1 R_2

over a tree whose leaves of 24 to 32 terms are exact backward
recurrences in mpn arithmetic.  The polynomial values are computed by
Horner's rule in one word (a whole leaf at a time), two words, or
generic multi-word two's complement.  When `R` divides `P` (as for
`\pi`, `\log 2` and most formulas) the leaf uses
`T \leftarrow R(k)(T + A(k) Q_s)` with `A = P/R`.  A large constant
content of `Q` (or `R`), as `q^2` for `\operatorname{atan}(p/q)`, is
split off and supplied from a table of powers over a balanced tree.

The number of terms comes from a rigorous tail bound: an exact scan of
`\prod |R(j)/Q(j)|`, with a closed form beyond the point where every
polynomial is dominated by its leading term (or a Taylor-shifted bound
before that), and the bound is added to the result.  Below about 30
limbs the fixed setup cost (about 1--2 microseconds) dominates.

Machin terms as Zuniga series
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each term of a Machin-type formula for logarithms is the logarithm of a
ratio of smooth numbers,

.. math::

    2 \operatorname{atanh}(1/x) = \log\frac{u}{v}, \qquad
    \frac{u}{v} = \frac{x+1}{x-1} \text{ in lowest terms}

(so `u - v = 1` for odd `x` and `u - v = 2` for even `x`).  Zuniga's
family of Ramanujan-type series [Zun2025]_ gives, for rationals
`u/v > 1`,

.. math::

    \log\frac{u}{v} = \sum_{n=1}^{\infty} \rho^n
        \frac{\alpha n + \beta}{\gamma\, n (2n-1)}
        \frac{(1)_n (1/2)_n}{(1/6)_n (5/6)_n},
    \qquad
    \rho = \frac{(u-v)^6}{108\, u^2 v^2 (u+v)^2},

.. math::

    \alpha = -2 (u+v)(u^2 - 14uv + v^2)(u^2 + 4uv + v^2), \quad
    \beta = (u+v)^3 (u^2 - 8uv + v^2), \quad
    \gamma = 2 (u-v)^5.

For instance, `u/v = 3` gives Zuniga's series for `\log 3`, with
`\rho = 1/243` and `(\alpha n + \beta)/\gamma = 88n - 14`, and
`u/v = 2` gives the series with `\rho = 1/3888` used for `\log 2`.
Substituting `u = x + 1`, `v = x - 1` (the series is homogeneous in
`u, v`) gives the direct translation of a Machin term:

.. math::

    2 \operatorname{atanh}(1/x) = \frac{x}{4} \sum_{n=1}^{\infty}
        \rho^n \frac{2 (3x^2-1)(3x^2-4)\, n - x^2 (3x^2-5)}{n (2n-1)}
        \frac{(1)_n (1/2)_n}{(1/6)_n (5/6)_n},
    \qquad \rho = \frac{4}{27\, x^2 (x^2-1)^2}.

The arguments and the coefficient matrix of the Machin formula are
unchanged; only the series evaluating each term is replaced.  In the
format of :func:`mp_real_hypgeom_series`, with
`\mathrm{num}/\mathrm{den} = (u-v)^6 / (6 u^2 v^2 (u+v)^2)` in lowest
terms:

.. math::

    P &= [-(u+v)^2 (u^2 - 8uv + v^2),\; 2 (u^2 - 14uv + v^2)(u^2 + 4uv + v^2)], \\
    Q &= \mathrm{den} \cdot [5, -36, 36], \qquad
    R = \mathrm{num} \cdot [0, -1, 2], \\
    \frac{c_P}{c_D} &= -\frac{(u+v)\,\mathrm{num}}{2 (u-v)^5}, \qquad c_Q = 0.

**Cost.**  One term of the series gains `6 \log_2 x + \log_2(27/4)`
bits, as much as three terms of the Taylor series of
`\operatorname{atanh}(1/x)` and 2.75 bits more.  In binary splitting
the cost per bit of precision is governed by the number of bits by
which the numbers outgrow the precision gained per term:

* about `\log_2 (2n)` against `2 \log_2 x` for the Taylor series;
* about `2 \log_2 n + 1` against `6 \log_2 x + 2.75` for Zuniga's
  series.

This is two thirds as much relatively, and much less for small `x`.
In exchange the terms are three times longer, which costs more at low
precision and for large `x`.  Accordingly, Zuniga's series is used for
arguments of `b` bits from a precision of about `b^2/2` limbs.

Measured on the logarithms of the first *num* primes
(:func:`_mp_real_log_primes_vec`, 64-bit, one thread), the speedup of
Zuniga's series over the Taylor series for the same Machin formulas:

    ===========  ===========  ===========  ===========
    *num*        `10^6` bits  `10^7` bits  `10^8` bits
    ===========  ===========  ===========  ===========
    2, 4         1.45         1.36         1.49
    8            1.23         1.31         1.28
    13           1.13         1.15         1.19
    20           1.06         1.05
    32--48       1.03--1.05   1.06--1.07
    ===========  ===========  ===========  ===========

For 2--4 primes (arguments of 8--13 bits) Zuniga's series is faster
from the smallest precisions at which the tables are computed, and for
32--48 primes (arguments of 62--87 bits) from about 2000--3000 limbs.
For 2 and 3 primes the set of 4 primes is used; the best dedicated
sets, `x = 7, 17` and `x = 31, 49, 161`, are no faster with Zuniga's
series.

**Arctangents.**  The same series with `u, v = t \pm i d` gives
`\log(u/v) = 2i \operatorname{atan}(d/t)` [Zun2025b]_: in terms of
`s = u + v = 2t` and `p = uv = t^2 + d^2`, which are real,

.. math::

    \rho = -\frac{4 d^6}{27\, t^2 p^2}, \qquad
    \frac{\mathrm{num}}{\mathrm{den}} = -\frac{8 d^6}{3\, t^2 p^2}, \qquad
    \frac{c_P}{c_D} = \frac{t\, \mathrm{num}}{64\, d^5},

.. math::

    P = [-s^2 (s^2 - 10p),\; 2 (s^2 - 16p)(s^2 + 2p)],

an alternating series which for `\operatorname{atan}(1/x)` gains the
same `6 \log_2 x + \log_2(27/4)` bits per term as for
`\operatorname{atanh}(1/x)` (:func:`_mp_real_atan_ratio_zuniga`).  The
terms of the Machin formulas for the arguments of Gaussian primes (and
their followup terms `\operatorname{atan}(p/q)`) use it under the same
rule, from about `b^2/2` limbs for arguments of `b` bits.  Measured on
the angles of the first *num* Gaussian primes
(:func:`_mp_real_atan_gauss_vec`, 64-bit, one thread), the speedup over
the Taylor series:

    ===========  ===========  ===========  ===========
    *num*        `10^5` bits  `10^6` bits  `10^7` bits
    ===========  ===========  ===========  ===========
    4            1.9          1.9          1.7
    8            1.2          1.2
    13           1.3          1.2          1.2
    20           1.0          1.4
    32--48       0.9--1.2     1.0          1.0
    ===========  ===========  ===========  ===========

Fixed-point Newton inverses and roots
-------------------------------------------------------------------------------

The *newton* suffix flags that the results of these functions are not
ulp-accurate: their error bounds are relative to the output.

.. function:: void _mp_real_inv_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n)
              void _mp_real_inv_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

    Given `(a, an)` with `a_{an-1} \ne 0` representing a fixed-point number
    `a \in [1/B, 1)` with `an` fraction limbs, sets `(q, n+2)` to an
    approximation of `1/a \in (1, B]` with `n` fraction limbs and two
    integral limbs (the highest limb may be zero).  The absolute error is
    bounded by `4 B^{-n} / a`.  The basecase computes the truncated
    reciprocal of the top `\min(an, n + 1)` limbs of `a` with
    :func:`_flint_mpn_inv_basecase`.

    The main function runs a Newton iteration on middle products, each
    step taking one of two forms from an approximation `T` of `1/a`:

    .. math::

        T' = T + T (1 - aT) \quad \text{(second order, from } m \approx n/2 \text{ limbs)},

    .. math::

        T' = T (1 - e + e^2), \; e = aT - 1
        \quad \text{(third order, from } m \approx n/3 \text{ limbs)}.

    In the third-order form the band of `aT` is merely wider, `e^2`
    comes from its top `n - 2m` limbs, and one middle product applies
    the correction.  The second-order step costs about one full
    multiplication at `n` limbs (the residual being a single middle
    product), and a third-order step 1.1 to 1.5 times as much, for a
    recursion over a third of the precision instead of half; the order
    is therefore chosen by precision: third-order for `24 \le n \le 80`
    and `n \ge 1000`, second-order otherwise.

.. function:: void _mp_real_div_newton_invmul(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n)
              void _mp_real_div_newton(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n)

    Given a numerator `(b, bn)` with `bn \ge 1` fraction limbs representing
    `b \in [0, 1)` and a denominator `(a, an)` with `a_{an-1} \ne 0`
    representing `a \in [1/B, 1)`, sets `(q, n+2)` to an approximation of
    `b/a` with `n` fraction limbs and two integral limbs.  The absolute
    error is bounded by `4 B^{-n} / a`.  The *invmul* variant multiplies
    the numerator by :func:`_mp_real_inv_newton`; the main function
    performs a Karp-Markstein iteration.

.. function:: void _mp_real_rsqrt_ui_newton_basecase(nn_ptr res, ulong a, slong n)
              void _mp_real_rsqrt_ui_newton(nn_ptr res, ulong a, slong n)

    Sets `(res, n)` to the fraction limbs of an approximation of
    `1/\sqrt{a}`, requiring `2 \le a < B`.  The error is bounded by
    `2 B^{-n}`.  The main function uses the third-order steps of
    :func:`_mp_real_rsqrt_newton`, with the residual `a y^2 - 1` exact and
    short (the low limbs of `y^2` times `a`).

.. function:: void _mp_real_rsqrt_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n)
              void _mp_real_rsqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

    Given `(a, an)` representing `a \in [B^{-2}, 1)` with `an` fraction
    limbs (at least one of the two highest limbs must be nonzero), sets
    `(q, n+2)` to an approximation of `1/\sqrt{a} \in (1, B]` with `n`
    fraction limbs and two integral limbs.  The absolute error is bounded
    by `4 B^{-n} / \sqrt{a}`.  The basecase combines
    ``flint_mpn_sqrtrem`` and ``mpn_tdiv_qr``.

    The main function performs third-order Newton steps from
    `T \approx a^{-1/2}` at `m \approx n/3` limbs (recursively):

    .. math::

        u = 1 - a T^2, \qquad
        w = \frac{u}{2} + \frac{3 u^2}{8}, \qquad
        q = T + T w.

    `u` is formed as the band of the product at weights `B^{-(n+2)}` to
    `B^{-(m-2)}`, with the unit wrapped into a control limb (a short `a`
    can start the band below the product, whose missing limbs are zero),
    `u^2` from its top `n - 2m` limbs, `w` in units two limbs below the
    output before those guard limbs are dropped, and `T w` by a middle
    product.  The neglected tail `(5/16) u^3 T` is below the rounding for
    `3m \ge n + 4`.  :func:`_mp_real_sqrt_newton` takes its reciprocal
    root from this function, and so does ``flint_mpn_sqrtrem`` in its
    Newton regime.

.. function:: void _mp_real_sqrt_newton_rsqrtmul(nn_ptr q, nn_srcptr a, slong an, slong n)
              void _mp_real_sqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

    Input as for :func:`_mp_real_rsqrt_newton`; sets `(q, n+2)` to an
    approximation of `\sqrt{a} \in [1/B, 1)` (the computed value can
    round up to 1) with absolute error bounded by `4 B^{-n} / \sqrt{a}`.
    Note that the error is proportional to `1/\sqrt{a}` rather than to
    the output.  The main function performs a Karp-Markstein iteration.

Fixed-point elementary functions
-------------------------------------------------------------------------------

The following routines implement evaluation of elementary functions
following [HJ2024]_ and [Joh2014c]_.  The exponential and the
trigonometric functions come with three kinds of argument reduction,
which trade precomputation against the cost per evaluation:

* **bitwise** (BKM-style) reduction with tables of `\log(1 + 2^{-i})`
  resp. `\operatorname{atan}(2^{-i})`, `i \le r`: the fastest per call,
  with a table whose cost is linear in `r n`; used up to 600 limbs;
* **diophantine** (multi-prime) reduction by precomputed logarithms of
  primes resp. arguments of Gaussian primes: one Machin-type
  evaluation per prime, cheaper to precompute; used from 600 to 65536
  limbs;
* **no table** (*notab*): halvings of the argument and a bit-burst
  evaluation; used above 65536 limbs.

The logarithm and the arctangent use the bitwise reduction up to 600
limbs, and above that one Newton-Taylor step over the diophantine
exponential resp. sine and cosine (or the AGM, for the logarithm; see
`Newton-Taylor inverses and the AGM logarithm`_).

Table-driven kernels at one to three limbs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On 64-bit machines, the kernels at one to three fraction limbs (the
dispatchers ``_mp_real_exp_opt``, ``_mp_real_log1p_opt``,
``_mp_real_sin_cos_opt``, ``_mp_real_atan_opt`` in ``impl.h``, and through
them the elementary functions up to 192 bits, the tangent with one
division, and the native ``nfloat`` functions at one and two limbs) are
short table-driven evaluations in registers, defined as inline functions
in ``src/mp_real/small.h``, with static tables (``small_tables.c``,
generated by ``dev/gen_mp_real_small_tables.py``, every entry at three
limbs so that all sizes share them; the one-limb kernels read the top
limbs):

* exp: `v = j_1/2^7 + j_2/2^{14} (+ j_3/2^{21}) + r` with tables of
  `\exp(j/2^7)` and `\exp(j/2^{14}) - 1`, `\exp(j/2^{21}) - 1`
  (128 entries each), and a Taylor polynomial in `r`;
* log1p: a multiplicative reduction over three (one and two limbs)
  resp. four levels of 7 bits with short one-limb factors `1 - D_k[j]`
  (exact products) and tables of `-\log(1 - D_k[j])`, then the Taylor
  polynomial of `\log(1 + u)`, `u < 2^{-21}` resp. `2^{-28}`;
* sin and cos: the pairs `(\sin c, 1 - \cos c)` for `c = j/2^7` and
  `j/2^{14}` combined with the Taylor polynomials of the remainder by
  angle addition;
* atan: `\operatorname{atan}(t) = \operatorname{atan}(c) +
  \operatorname{atan}(w)`, `c = j/2^9` (513 table entries), `w = (t -
  c)/(1 + tc) < 2^{-9}` by a division in registers, and the odd series
  in `w`. For a ratio `t = Y/X`, the same division gives
  `w = (512 Y - jX)/(512 X + jY)` with an exact numerator and
  denominator, so that `\operatorname{atan}(Y/X)` costs one division
  (and one short division for `j`) rather than a division followed by
  the kernel.

The polynomials taper their precision (one limb for the innermost terms
on a scaled argument, then two, then three) and the products keep the
top limbs only; at one limb every product is a single high word (and
the arctangent's division one hardware division).  The scalar cores
(``_mp_real_small_exp_2`` and so on, on and to limbs passed in registers)
are shared with the ``nfloat`` fast paths.  Against the register bitwise
kernels they replace at one limb, they are 3.9 (exp), 7 (log1p), 2.6
(atan) and 3.2 (sin and cos) times faster.

.. function:: void _mp_real_small_exp_1_mpn(nn_ptr y, ulong * err, nn_srcptr v)
              void _mp_real_small_exp_2_mpn(nn_ptr y, ulong * err, nn_srcptr v)
              void _mp_real_small_exp_3_mpn(nn_ptr y, ulong * err, nn_srcptr v)
              void _mp_real_small_log1p_1_mpn(nn_ptr y, ulong * err, nn_srcptr v)
              void _mp_real_small_log1p_2_mpn(nn_ptr y, ulong * err, nn_srcptr v)
              void _mp_real_small_log1p_3_mpn(nn_ptr y, ulong * err, nn_srcptr v)
              void _mp_real_small_sin_cos_1_mpn(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v)
              void _mp_real_small_sin_cos_2_mpn(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v)
              void _mp_real_small_sin_cos_3_mpn(nn_ptr ys, nn_ptr yc, ulong * err, nn_srcptr v)
              void _mp_real_small_atan_1_mpn(nn_ptr y, ulong * err, nn_srcptr v)
              void _mp_real_small_atan_2_mpn(nn_ptr y, ulong * err, nn_srcptr v)
              void _mp_real_small_atan_3_mpn(nn_ptr y, ulong * err, nn_srcptr v)

    For `v \in [0, 1)` at `n = 1`, `2` resp. `3` fraction limbs, set *y* to
    `\exp(v)` (`n + 1` limbs, with an integral limb), `\log(1 + v)`
    (`n` limbs), `\sin(v)` and `\cos(v)` (*ys*, *yc*, `n + 1` limbs each)
    resp. `\operatorname{atan}(v)` (`n` limbs), and *err* to a bound for
    the error in ulps of `B^{-n}` (in either direction; between 8 and
    64). The contract is that of the other fixed-point kernels. Only
    available with 64-bit limbs.

.. function:: void _mp_real_small_exp_pm_2(ulong * y2, ulong * y1, ulong * y0, ulong * z2, ulong * z1, ulong * z0, ulong v1, ulong v0)
              void _mp_real_small_exp_pm_3(ulong * y3, ulong * y2, ulong * y1, ulong * y0, ulong * z3, ulong * z2, ulong * z1, ulong * z0, ulong v2, ulong v1, ulong v0)

    For `v \in [0, 1)` at two resp. three fraction limbs (passed most
    significant limb first), set *y* to `\exp(v)` and *z* to `\exp(-v)`,
    each with an integral limb and two resp. three fraction limbs, within
    ``MP_REAL_SMALL_EXP_ERR_2`` resp. ``MP_REAL_SMALL_EXP_ERR_3`` ulps.
    The two evaluations share the splitting of `v` and the even and odd
    parts of the Taylor polynomial (with tables of `\exp(-j/2^7)` and
    `1 - \exp(-j/2^{14})`, `1 - \exp(-j/2^{21})` for the negative
    side), at about 1.5 times the cost of one exponential. Inline, 64-bit
    only.

.. function:: void _mp_real_small_divq_3(ulong * w2, ulong * w1, ulong * w0, ulong n2, ulong n1, ulong n0, ulong d2, ulong d1, ulong d0)

    For fractions `N < D` at three limbs with `D` normalized (top bit
    set), set `w = \lfloor N B^3 / D \rfloor`, the quotient at three
    fraction limbs (three 3-by-2 divisions in registers). Inline, 64-bit
    only.

.. function:: void _mp_real_small_atan_ratio_2(ulong * r1, ulong * r0, ulong y1, ulong y0, ulong x1, ulong x0)
              void _mp_real_small_atan_ratio_3(ulong * r2, ulong * r1, ulong * r0, ulong y2, ulong y1, ulong y0, ulong x2, ulong x1, ulong x0)

    For fractions `X \ge 1/2` and `Y \le X` at two resp. three limbs
    (passed and returned most significant limb first), set the result to
    `\operatorname{atan}(Y/X)` at the same number of fraction limbs,
    within ``MP_REAL_SMALL_ATAN_ERR_2`` resp. ``MP_REAL_SMALL_ATAN_ERR_3``
    ulps (the bounds of the kernels for `\operatorname{atan}(v)`). Inline
    functions, 64-bit only.

Dense-table kernels at four to ten limbs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On 64-bit machines, the kernels at four to ten fraction limbs (the
``r = 0`` per-size implementations ``exp_opt_<n>.c``, ``log1p_opt_<n>.c``,
``atan_opt_<n>.c`` and ``sin_cos_opt_<n>.c``, reached through the same
dispatchers, so that they serve the elementary functions from 193 to 640
bits and the native ``nfloat`` functions at three to nine limbs) reduce
the argument seven bits per level against static tables instead of one
bit per step.  Every level's multiplier is exact and short -- a single
limb -- so that applying it costs one row of single-limb products, and
the remainder after `L` levels lies below `2^{-(7L - 1)}`, which is the
`r` of these kernels (the reduced radius; the error bounds keep the forms
`9r + 100`, `3r + 64`, `4r + 64` and `6r + 128` of the bitwise kernels,
which the dense kernels stay far below, typically within ten to twenty
ulps and within fifty for exp at `r = 55`):

* exp, `v \in [0, 1)`: `t \gets v - \log(2^s (1 + C_1[j_1]))` (`j_1`
  the top seven bits, `C_1[j_1] = \lfloor (e^{j_1/2^7}/2^s - 1) 2^{64}
  \rfloor / 2^{64}` with `s = 1` for `e^{j_1/2^7} \ge 2`), then
  `t \gets t - \log(1 + j_k 2^{-7k})` with `j_k = \lfloor t 2^{7k} \rfloor`
  (tables of the floors, so that `t` stays nonnegative; the overshoot of
  a level makes `j_{k+1}` up to 190 at the third level and 130 beyond),
  the series of `\exp(t)`, the factors `1 + j_k 2^{-7k}` and
  `1 + C_1[j_1]` by single-limb products, and the shift by `s`;
* log1p: the multiplicative reduction of the two- and three-limb
  kernels over all levels, `u \gets (1 + u)(1 - D) - 1` with
  `D_1 = \lfloor 2^{64} j/(2^7 + j)\rfloor / 2^{64}` tabulated and, from the
  second level on, `D = c - \lceil c^2 2^{64} \rceil / 2^{64}`,
  `c = j 2^{-7k}`, formed in place (`D \le c/(1 + c)` keeps `u \ge 0`),
  with tables of `-\log(1 - D)`; then the series of `\log(1 + u)` itself
  (no division);
* atan: the vector `(1, v)` rotated by the factors `1 - i c_k`,
  `c_1 = j_1/2^7`, `c_k = j_k 2^{-7k}`, with the index `j_k` from an
  estimate of `Y/X` by a one-limb reciprocal of `X` (lowered by `2^{-12}`
  so that `c_k \le Y/X`), tables of `\operatorname{atan}(c)`, and one
  reciprocal of `X` for the final ratio;
* sin and cos: on `w = v/2`, the angles `\operatorname{atan}(c_k)` are
  subtracted (`c_1 = \lfloor \tan(j_1/2^7) 2^{64} \rfloor / 2^{64}` and
  `c_k = j_k 2^{-7k}`, the latter sharing the atan tables), the pair
  `(\cos t, \sin t)` of the series is rotated by the factors
  `1 + i c_k`, which scale it by `\sqrt{N}`, `N = \prod (1 + c_k^2)`,
  and the double-angle formulas `\cos v = (X^2 - Y^2)/N`,
  `\sin v = 2XY/N` absorb the scale with one reciprocal of
  `N = X^2 + Y^2` and no square root.

The single-limb rows are emitted as straight-line products followed by
add and subtract chains (the ``longlong.h`` macros, up to eleven limbs);
at these sizes a level costs about a quarter of what
``mpn_addmul_1`` and its copies cost.  The reciprocals are the inline
``_flint_mpn_recip_small`` of ``mpn_extras/recip_small.h``.

The series that finishes a kernel is either a hand series built for
exactly `t < 2^{-r}` (``dev/gen_mp_real_exp_hand_mixed.py``,
``dev/gen_mp_real_trig_hand_mixed.py`` and the log1p series of
``dev/gen_mp_real_dense_bodies.py``) or the hardcoded rectangular
splitting of :func:`_mp_real_exp_rs`, :func:`_mp_real_atan_rs` and
:func:`_mp_real_sin_cos_rs`, whose tiers sit at `2^{-32}` and `2^{-64}`;
two values of `r` select the latter: `r = 32` (five levels, the
remainder below `2^{-34}`, for the tier of `2^{-32}`) and `r = 63` (nine
levels and a *top level* that takes the whole top limb of the remainder,
`t < (1 + 2^{-8}) 2^{-63}`, as its index `j \le 2` with the factor
`1 + j 2^{-64}` and a table of four values, leaving `t < 2^{-64}` but for
rare inputs, where the limb below is within a few units of `2^{64}` and
the series takes the tier of `2^{-32}`).  Measured by
``dev/tune_mp_real.py``, the hardcoded tiers win for exp at nine and ten
limbs (`r = 63`, 14 to 18 percent faster than the hand series at
`r = 62`) and for sin and cos at eight to ten limbs (`r = 32`, 4 to 9
percent faster than at `r = 41`); elsewhere the hand series built for the
exact radius stay ahead: a tier serves any argument below its threshold,
and from `2^{-32}` it evaluates more terms than a remainder below
`2^{-48}` needs, while reaching `2^{-64}` takes one more full level per
seven bits (and for atan and sin/cos rotations by two rows), which only
pays when the top level closes the gap of two bits after nine levels.
The log1p kernels keep their own series because the hardcoded atanh
needs a division by `2 + u`.

The tables (``dense_tables.c``, generated by
``dev/gen_mp_real_dense_tables.py``, which also checks the index and
final remainder bounds) store every value from the third level on at
ten fraction limbs, the smaller sizes reading the top limbs of an entry:
exp 7 levels and the top level, log1p 7 levels, atan 4 levels (shared
with sin and cos), 191 KB in all; the first two levels (and the first
level of sin and cos) come from the 21-limb tables that the run-time
kernels use (below), whose top ten limbs are the same values.  The bodies are
generated by ``dev/gen_mp_real_dense_bodies.py`` through
``dev/tune_mp_real.py`` (the functions ``exp``, ``log1p``, ``atan`` and
``sin_cos``, with `r = 7L - 1`, 32 or 63), straight-line code on named
limbs except for the series,
the reciprocal and the high products.  Against the bitwise per-size
kernels at four to seven limbs they are about 2 times faster for exp,
2.8 to 3.4 times for log1p, 1.8 to 2.2 times for atan and 1.6 to 1.8
times for sin and cos; at eight to ten limbs, against the generic
bitwise kernels, 2 to 2.5, 3.5 to 4.4, 2.7 to 3.8 and 1.5 to 1.75
times.  The tangent takes the sine and cosine kernels and one division
(:func:`_mp_real_tan_bitwise_rs`).

.. function:: void _mp_real_exp_opt_4(nn_ptr res, nn_srcptr x)
              void _mp_real_log1p_opt_4(nn_ptr res, nn_srcptr x)
              void _mp_real_atan_opt_4(nn_ptr res, nn_srcptr x)
              void _mp_real_sin_cos_opt_4(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)

    And likewise for `n = 5, \ldots, 10`: for `x \in [0, 1)` at `n`
    fraction limbs, set `(res, n + 1)` to `\exp(x)`, `(res, n)` to
    `\log(1 + x)` resp. `\operatorname{atan}(x)`, and `(ysin, n + 1)`,
    `(ycos, n + 1)` to `\sin(x)` and `\cos(x)` (each with an integral
    limb), within the bounds that the dispatchers report
    (``MP_REAL_*_OPT_R``, ``MP_REAL_SIN_COS_OPT_R`` in ``impl.h``).
    64-bit only.

Dense-table kernels with run-time tables, eleven to thirty-two limbs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Past the static tables (``MP_REAL_DENSE_W + 1 = 11`` up to
``MP_REAL_DENSE_RT_MAX = 32`` fraction limbs, so that they serve the
elementary functions from 641 to 2048 bits and the native ``nfloat``
functions at ten to thirty-one limbs), the dispatchers run the same four
dense reductions as generic loops (``dense_rt.c``) on tables partly
computed at run time: nine levels of seven bits and the top level, for
the series of the tier of `2^{-64}` (which the hardcoded routines of
``exp_rs.c`` and ``trig_rs.c`` cover up to 21, 22 and 35 limbs for exp,
sin and cos, and atan and atanh, and :func:`_mp_real_series_rs` beyond;
measured against the generic tapered rectangular splitting at the same
radius, the hardcoded routines stay 1.4 to 2.2 times faster up to about
20 limbs, and reach parity around 22 limbs for exp and sin/cos, where
they stop, and beyond 24 for atan; the deeper reduction that a tier of
`2^{-128}` would need costs more levels than the shorter series saves at
these sizes).  Every multiplier is one limb (the top level of log1p has
`D = c - c^2` as two limbs); log1p and atan end with one division
(:func:`flint_mpn_divapprox_fraction`: `u/(2 + u)` for
`2 \operatorname{atanh}`, resp. `Y/X`), sin and cos with one reciprocal
of `N` (:func:`_flint_mpn_recip_small` up to eleven limbs,
:func:`flint_mpn_invapprox` beyond).  The bounds reported are those of
the bitwise form with ``MP_REAL_DENSE_RT_R = 63``; the measured errors
stay within twenty ulps.

The kernel bodies are written once, inline in the size `n`, and
instantiated straight-line at 11 to 16 limbs (where add-with-carry
intrinsics exist: x86-64, or ``__builtin_addcll``): with a constant `n`
the carry chains unroll into add-with-carry sequences, and on builds
with the x86-64 ADX assembly every row `d \pm \lfloor s c / B \rfloor`
is one pass of ``mulx``, ``adcx`` and ``adox``
(``dev/gen_mp_real_dense_rt_asm.py``), about half the cost of
``mpn_mul_1`` and ``mpn_add_n`` or of the same row in C.  At 11 to 14
limbs log1p ends with the hand series of the static kernels at `r = 63`
(``dev/gen_mp_real_dense_rt_series.py``) instead of the division and the
atanh series.  Against the generic loops (the other sizes, and other
machines), measured: exp 1.12 to 1.21 times faster, log1p 1.11 to 1.47,
atan 1.19 to 1.29, sin and cos 1.11 to 1.24 (the unrolled chains in C
alone, without the rows in assembly, gave 3 to 10 percent); about 100 KB
of code for the 24 instances.  All variants produce the same bits.  At
17 to 21 limbs the same instances would gain 1.06 to 1.19 times for
about 85 KB more.

The tables hold floors (or lower bounds within one ulp), a kernel at
`n` limbs reading the top `n` limbs of an entry.  The first two levels
of every family and the first level of sin and cos are static, at
``MP_REAL_DENSE_RT_W = 32`` limbs (``dense_rt_tables.c``, emitted by
``dev/gen_mp_real_dense_rt_tables.py`` with mpmath, 215 KB): their
arguments are the largest, so they are the most expensive entries to
compute.  The static per-size kernels at four to ten limbs read their
first two levels from the same tables (their ten-limb tables,
``dense_tables.c``, start at the third level: 195 KB).  The deeper
levels are filled as they are needed, in two tiers: at 21 limbs for the
kernels up to 21 limbs and at 32 limbs beyond, so that a process working
at the smaller sizes fills only the cheaper entries.  The storage of a
family in a tier is a shared cache (`Shared tables`_) with a single
version, allocated on its first use, so that a reader finds it with one
acquire load; every entry has a state byte (empty, being written,
ready).  A kernel that meets an entry that is not ready computes it into
a buffer of its own at one more limb, uses it, and stores it in the
table if it wins the compare-and-swap of the state from empty to being
written, publishing it with a release store of the state ready; no
reader waits, and no arb is involved.  An entry is `\log(1 + c)`,
`-\log(1 - D)` or `\operatorname{atan}(c)` for a one-limb `c` (two
limbs for the top level of log1p), evaluated by the Horner scheme in
`c` (every step one or two passes of ``mpn_mul_1``, the denominators
gathered into running products of a limb, tapered as the powers of `c`
grow) and lowered by its error bound; it costs 0.5 to 3 microseconds
at 21 limbs (all 900 to 970 entries of a family 1 to 1.5 ms) and 1 to
5 microseconds at 32 limbs (2 to 3 ms).  The storage, 163, 152 and
158 KB for exp, log1p and atan (with sin and cos) in the first tier and
249, 231 and 241 KB in the second, lives until
:func:`flint_cleanup_master`, and only the pages holding written
entries are touched.

At 11 to 21 limbs the kernels are 1.9 to 2.5 times faster than the
bitwise kernels for exp, 1.7 to 2 times for log1p, 1.9 to 2.7 times for
atan and 1.4 to 1.7 times for sin and cos; at 22 to 32 limbs 1.5 to 1.7
times for exp, 1.5 to 1.8 times for log1p, 1.9 to 2.4 times for atan
and 1.1 to 1.2 times for sin and cos; the static kernels remain 1.15 to
1.75 times faster than these loops at ten limbs.  The first call at a
run-time size costs about 0.3 ms (most of it the constants at that
precision) instead of 0.5 to 6 ms for building the bitwise tables; the
first hundred calls on random arguments 1 to 2.3 ms at 11 to 21 limbs
(bitwise: 1.8 to 6.9 ms), but 1.8 to 2.6 ms at 22 to 32 limbs (bitwise:
0.9 to 2 ms), where the run-time entries cost about as much as a
bitwise evaluation and pay for themselves after about a thousand calls.

At equal tables, the dense kernels would stay ahead of the bitwise
kernels well beyond 32 limbs (measured: up to about 85 limbs for exp,
90 for log1p, 48 for sin and cos, and beyond 128 for atan, by 1.2 to
1.5 times at 64 limbs except for sin and cos), but every limb of width
adds 6.7 KB to the static levels, and the run-time entries at 64 limbs
cost 3 to 20 microseconds each, more than a bitwise evaluation, so that
a process would need several thousand calls to recover them.

.. function:: void _mp_real_exp_dense_rt(nn_ptr y, nn_srcptr v, slong n)
              void _mp_real_log1p_dense_rt(nn_ptr y, nn_srcptr v, slong n)
              void _mp_real_atan_dense_rt(nn_ptr y, nn_srcptr v, slong n)
              void _mp_real_sin_cos_dense_rt(nn_ptr ysin, nn_ptr ycos, nn_srcptr v, slong n)

    For `v \in [0, 1)` at `n` fraction limbs, `2 \le n \le`
    ``MP_REAL_DENSE_RT_MAX``, set `(y, n + 1)` to `\exp(v)`, `(y, n)` to
    `\log(1 + v)` resp. `\operatorname{atan}(v)`, and `(ysin, n + 1)`,
    `(ycos, n + 1)` to `\sin(v)` and `\cos(v)`, by the run-time dense
    tables (built on the first call).  64-bit only; the dispatchers call
    them from eleven limbs on.

.. function:: slong _mp_real_dense_rt_prepare(int which, slong n)

    Fills every entry of the run-time tables of a family
    (``MP_REAL_DENSE_RT_EXP``, ``MP_REAL_DENSE_RT_LOG1P`` or
    ``MP_REAL_DENSE_RT_ATAN``) that is not there yet (the kernels
    otherwise fill them as they need them), in the tier serving the
    kernels at `n` limbs, and returns its width: 21 for `n \le 21`, else
    ``MP_REAL_DENSE_RT_MAX``.


.. function:: void _flint_mpn_recip_small(nn_ptr I, nn_srcptr F, slong n, slong bits)

    For an `n`-limb fraction `F`, `4 \le n \le 11`, sets `(I, n)` to a
    lower bound of `1/(1 + F)` within six ulps (within `2^{-bits}` for a
    caller needing fewer bits, which may lower the order of the last
    step): a one-limb reciprocal by one two-by-one division, a cubic step
    at three limbs whose products are all by that limb, and one step at
    `n` limbs, of order two, three or four as `n` needs, on products by
    the three limbs (the residual has two zero top limbs).  An inline
    function in ``mpn_extras/recip_small.h`` for a compile-time `n`,
    64-bit only.  Inlined into a kernel it saves 20 to 25 percent of the
    arctangent at four to seven limbs against
    :func:`flint_mpn_divapprox_fraction`; as a general division with a
    run-time size, normalization and the floor-or-one-more contract it
    is not faster than the existing small divisions, which is why
    those keep their algorithms.

Series evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The series evaluation functions require an argument reduced below
`2^{-32}` (checked with a ``FLINT_ASSERT``) and dispatch internally on
the top limbs of `x`: hardcoded rectangular splitting (64-bit) for
`2^{-64} \le x < 2^{-32}` and `2^{-128} \le x < 2^{-64}` at small sizes,
and otherwise the tapered rectangular splitting of
:func:`_mp_real_series_rs` at the argument's actual number of leading
zero bits.  The hardcoded tiers are straight-line routines (generated
by ``dev/gen_mp_real_exp_rs_hard.py`` and
``dev/gen_mp_real_trig_rs_hard.py``) only at the sizes on the hot paths
of the elementary functions: below `2^{-64}`, to which the dense kernels
and the larger bitwise kernels reduce, exp from 9 limbs on, atan and atanh from 11, sin and cos from
11; below `2^{-32}` only sin and cos at 8 to 10 limbs (the per-size
kernels pinned to `r = 32`).  A generic engine (``series_hard.c``, about
12 KB of code against about 80 KB of straight-line routines it
replaces) serves the other sizes and the hyperbolic functions: it makes
the same calls in the same order from the same tuned parameters, so
its results are identical, at a cost of a few hundred instructions per
call (twice the time at four limbs, 10 to 15 percent at twenty).  The bounds written to *err* are 10 ulps for the exponential
and 15 for the others; they are one-sided (the computed result never
exceeds the true value) for :func:`_mp_real_exp_rs`, the hyperbolic
functions and :func:`_mp_real_atanh_rs`, two-sided for the alternating
functions.

.. function:: void _mp_real_exp_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)

    Sets `(res, n + 1)` to an approximation of `\exp((x, n))`,
    requiring `x < 2^{-32}`; *err* at most 10, one-sided.

.. function:: void _mp_real_sin_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
              void _mp_real_cos_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
              void _mp_real_sin_cos_rs(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n)
              void _mp_real_sinh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
              void _mp_real_cosh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
              void _mp_real_sinh_cosh_rs(nn_ptr ysinh, nn_ptr ycosh, ulong * err, nn_srcptr x, slong n)

    Set `(res, n + 1)` to an approximation of the respective function
    of `(x, n)`, requiring `x < 2^{-32}`.  The combined versions allow
    either output pointer to be *NULL*.  *err* at most 15.

.. function:: void _mp_real_atan_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)
              void _mp_real_atanh_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n)

    Sets `(res, n)` to an approximation of `\operatorname{atan}((x, n))`
    resp. `\operatorname{atanh}((x, n))`, requiring `x < 2^{-32}`;
    *err* at most 15.

.. function:: void _mp_real_series_rs(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r, int func)
              void _mp_real_series_rs_sin_cos(nn_ptr ysin, nn_ptr yg, nn_srcptr x, slong n, flint_bitcnt_t r, int hyperbolic)

    Tapered rectangular splitting for a fraction `0 \le x < 2^{-r}`,
    `r \ge 8`, of `n` limbs (declared in ``mp_real/impl.h``).  *func*
    selects

    * `\exp(x)` or `\exp(-x)` (``MP_REAL_SERIES_EXP``, ``_EXP_NEG``;
      `(res, n + 1)` with a units limb, within 5 ulps);
    * `\sin`, `\sinh`, `\operatorname{atan}`, `\operatorname{atanh}`,
      `1 - \cos` or `\cosh - 1` (``_SIN``, ``_SINH``, ``_ATAN``,
      ``_ATANH``, ``_COS``, ``_COSH``; `(res, n)`, within 2 ulps).

    The second function computes the sine and `1 - \cos` (resp.
    `\sinh` and `\cosh - 1`) sharing the powers, either output possibly
    ``NULL``.  Outputs may alias *x*.

    With `v = x` (exp) or `v = x^2`, and `s = r` resp. `2r` the bits
    each power of `v` gains, the series

    .. math::

        P = \sum_{k < N} a_k v^k, \qquad
        a_k = \frac{1}{k!}, \; \frac{1}{(2k+1)!}, \; \frac{1}{(2k+2)!}
        \text{ or } \frac{1}{2k+1}

    (alternating for the circular functions and `\exp(-x)`; the result
    is `P`, `x P` or `v P`) is split into blocks of
    `m \approx \sqrt{N/1.5}` (even) terms,

    .. math::

        S_b = \sum_{i < m} a_{mb+i} v^i + v^m S_{b+1},

    with the powers `v, \ldots, v^m` computed once and every term a
    single ``mpn_addmul_1`` (resp. ``submul_1``) by a one-limb
    coefficient.  The factorial coefficients chain integrally downwards,
    and one ``mpn_divrem_1`` closes each run that fills a limb;
    `1/(2k+1)` uses the integer numerators `D/(2k+1)` of groups of odd
    numbers with product `D`, and one rescaling by `D'/D` per group.

    **Tapering.**  Block `b` enters the sum multiplied by
    `v^{mb} < 2^{-smb}`, so it runs in a window of

    .. math::

        \min\Bigl(n, \Bigl\lceil \frac{\mathrm{FLINT\_BITS} \cdot n + g - smb}{\mathrm{FLINT\_BITS}} \Bigr\rceil\Bigr)

    limbs, `g = \operatorname{bits}(N) + 3`, at bit granularity (not
    only when `x` has whole zero limbs).  Within a block `v^i` skips its
    whole zero limbs, and the product of `S_{b+1}` by `v^m` reads only
    the significant limbs of both.  The cost is about `m` full products
    for the powers, `N/(3m)` full products for the block boundaries and
    the scalar passes over the shrinking windows: about `2 \sqrt{N/3}`
    full products in all, against `N/3` for the tapered Horner scheme of
    :func:`_mp_real_series_tapered`, which is faster only below about a
    dozen limbs.

    In the alternating families the sum is kept in two's complement
    across its units limb; with `m` even every block starts with a
    positive term, and a division following a negative term offsets by
    the divisor to stay unsigned.  Errors committed in block `b \ge 1`
    are damped by `2^{-\min(g, sm)}`, those of block 0 are dominated by
    the powers' truncations weighted by `a_i`, and the final product by
    `x` (or `v`) damps everything by `2^{-r}`: hence 2 ulps, except for
    exp, whose sum is the result.

.. function:: void _mp_real_series_tapered(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r, int func)

    Sets `(res, n)` to `f((x, n))` for `0 \le x < 2^{-r}` within 6
    ulps, *func* selecting `f = \tan` (``MP_REAL_SERIES_TAN``),
    `\operatorname{atan}`, `\operatorname{atanh}`, `\sin` or
    `1 - \cos` (``_ATAN``, ``_ATANH``, ``_SIN``, ``_COS``).  With
    `z = x^2`, Horner's rule evaluates

    .. math::

        f(x) = x \pm x^3 V_0 \;\; (\text{resp. } z V_0), \qquad
        V_k = c_k \pm z V_{k+1},

    over TAPERED widths: level `k` enters the result scaled by
    `x^{e_k} < 2^{-r e_k}`, `e_k = 2k + 3` (resp. `2k + 2`), so it runs
    at `\lceil (\mathrm{FLINT\_BITS} \cdot n + g - r e_k) /
    \mathrm{FLINT\_BITS} \rceil` limbs (`g` a few guard bits), one
    ``flint_mpn_mulhigh_n`` per level.  Unlike the table-based kernels,
    the cost falls with `r`, i.e. with the argument's leading zero bits.

    The coefficients (for the tangent `c_k = T_{k+2}/(2k+3)!`, `T_m`
    the tangent numbers) are static tables generated by
    ``dev/gen_mp_real_elem_tables.py`` (``elem_tables.c``, for 64- and
    32-bit limbs): `\lfloor c_k B^{W_k} \rfloor` at the widest width
    `W_k` the supported range needs, whose top limbs serve every smaller
    width (floors nest), with an upper bound for `\log_2 c_k` to count
    the levels.  The supported ranges are `n \le 80`, `r \ge 32` for the
    tangent, `n \le 13`, `r \ge 10` for atan and atanh (which share the
    table), and `n \le 16`, `r \ge 18` for the sine and cosine.

.. function:: void _mp_real_series_rs_tan(nn_ptr res, nn_srcptr x, slong n, flint_bitcnt_t r)
              int _mp_real_series_rs_tan_ok(slong n, flint_bitcnt_t r)

    Sets `(res, n)` to `\tan((x, n))` for `0 \le x < 2^{-r}`,
    `r \ge 32`, within 3 ulps (declared in ``mp_real/impl.h``;
    *res* may alias *x*), for `n \le 80` and beyond as long as
    :func:`_mp_real_series_rs_tan_ok` returns nonzero.

    With `\tan t = t + t^3 S`, `S = \sum_k c_k z^k`, the coefficients
    are not ratios of small integers, but the leading ones share small
    denominators.  With `Q_j` the lcm of the reduced denominators of
    `c_0, \ldots, c_k`, chunk `j` holds the terms whose `Q_j` fits
    `j + 1` limbs (on 64-bit machines, 12 terms in one limb, then 5--7
    more per limb, 77 terms in 12 chunks), with the integer numerators
    `N_k = c_k Q_j` in static tables.  One tapered rectangular splitting
    (as in :func:`_mp_real_series_rs`) runs over the chunked terms,
    descending, carrying `Q_j` times the partial sum with `j + 1` units
    limbs:

    * a term costs `j + 1` ``mpn_addmul_1`` passes;
    * crossing into chunk `j - 1` costs one division by `Q_j / Q_{j-1}`;
    * the end is one ``mpn_divrem_1`` by `Q_0`.

    A chunk is used while its width is at most a third of the window at
    its first term; the rest is the tapered Horner tail of
    :func:`_mp_real_series_tapered`, entering as the top block's
    continuation.  At small sizes the Horner scheme alone is faster and
    is used.

Bitwise argument reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions compute the elementary functions on the unit
interval using bitwise (BKM-style) argument reduction to bring the
argument below a tuned threshold `2^{-r}`, followed by series
evaluation and reconstruction.  On 32-bit systems they require
`n \ge 2`.

.. function:: void _mp_real_exp_reduced(nn_ptr y, ulong * err, nn_srcptr t, slong n, flint_bitcnt_t r, int alg)

    Sets `(y, n + 1)` (`n` fraction limbs and a units limb) to an
    approximation of `\exp(t)` for a reduced argument `(t, n)` with
    `t < 2^{-r}`, `r \ge 16`, independent of any particular argument
    reduction (algorithms 1 and 2 additionally require `r \ge 32`).
    *err* is the rigorous bound computed at run time, at most 96.
    *alg* selects the internal method:

    * 0: the tuned automatic choice;
    * 1: the direct rectangular-splitting series (:func:`_mp_real_exp_rs`,
      whose output is already in this format, *err* = 10);
    * 2: the sinh series plus a square root;
    * 3: one bit-burst step: the leading slice of *t*, on limb
      boundaries, evaluated by binary splitting, and the remainder by
      the sinh series at the doubled rate;
    * 4: the full bit-burst algorithm with slice lengths doubling, which
      is asymptotically quasi-optimal for very large `n`.

    Except for the direct series, the arithmetic around the series
    kernels and the exact binary splitting is ``mp_real_t`` ball
    arithmetic: the kernels are imported with their documented errors,
    the sinh reconstruction is a ball squaring, sum and square root,
    each burst slice's exact splitting output `T / (Q B^{Q_e})` becomes
    the factor `Q B^{Q_e} + T` by one ball addition, the products are
    ball products at the frame precision, and the finish is one ball
    division.  The bound comes out rigorous and is checked against the
    budget on export.

    The burst machinery works at limb granularity throughout: slice
    boundaries and splitting-tree truncation frames are limb counts,
    and dead low limbs of a slice fold into its frame, so sparse
    arguments (a few significant limbs over a deep frame, the shape of
    a double-precision input at high precision) shrink the whole
    computation.  The automatic thresholds can be recalibrated with
    ``tune/tune-exp-reduced``.

.. function:: void _mp_real_exp_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)

    Sets `(res, n + 1)` to an approximation of `\exp((x, n))` for any
    `0 \le x < 1`, with *err* at most `9 r + 100` for the reduction
    parameter `r` used.  The bound grows linearly with `r` because all
    work happens at the output precision; callers wanting sub-ulp
    accuracy should pad the precision by one limb themselves.

    The argument is reduced below `2^{-r}` (`r` internally clamped to
    `\mathrm{FLINT\_BITS} \, n - 16`) by subtracting in turn each
    logarithm `L_i = \log(1 + 2^{-i})`, `i = 0, 1, \ldots, r`, for
    which `L_i \le x`; the Taylor series is evaluated on the reduced
    argument; and the used factors `1 + 2^{-i}` are multiplied back in,
    each by a single shift-and-add.  When the reduced series is long,
    `\sinh` is evaluated instead (half the terms) and the exponential
    reconstructed as `\exp(t) = \sinh(t) + \sqrt{1 + \sinh(t)^2}`.
    The logarithm table is generated at run time and kept in a shared
    cache at the largest precision and index range requested so far by
    any thread (`Shared tables`_).

    Passing `r = 0` selects a tuned default (see `Tuning the reduction
    parameter`_); explicit values require `r \ge 32`, the contract of
    the residual series.  Small `r` minimizes table and reduction work,
    large `r` shortens the series; the tuned values grow from 64 at a
    dozen limbs to 768 at several hundred.

    With `r = 0` the sizes `n \le 10` go to the per-size kernels on
    64-bit machines (which also set *err*): the table-driven kernels of
    `Table-driven kernels at one to three limbs`_ for `n \le 3`, and
    for `4 \le n \le 10` the dense-table kernels (`Dense-table kernels
    at four to ten limbs`_, one source file per size, ``exp_opt_<n>.c``,
    emitted and tuned by ``dev/tune_mp_real.py``), which reduce seven
    bits per level rather than bitwise.

    Otherwise, from ``MP_REAL_EXP_BITWISE_VEC_MIN = 4`` limbs on 64-bit
    machines, the work is done by :func:`_mp_real_exp_bitwise_vec`,
    with signed digits and vector arithmetic in radix `2^{56}`; the
    description above is that of the remaining sizes (and of 32-bit
    machines).

.. function:: void _mp_real_exp_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n + 1)` to an approximation of `\exp((x, n))` for
    `0 \le x < 1` and `32 \le r \le 64 n - 16`, within the bound of
    :func:`_mp_real_exp_bitwise_rs` (measured: below 10 ulps).  64-bit
    only; :func:`_mp_real_exp_bitwise_rs` calls it from
    ``MP_REAL_EXP_BITWISE_VEC_MIN`` limbs.

    The reduction uses signed digits: with `a_i = \log(1 + 2^{-i})` and
    `b_i = -\log(1 - 2^{-i})`, it writes
    `x = d_0 a_0 + \sum_{i=1}^{r+1} (a_i \text{ if } d_i = 1, -b_i \text{ if } d_i = -1) + t`
    with `d_i \in \{-1, 0, 1\}`, choosing `d_i` by comparing `x` with
    `\pm 0.75 \cdot 2^{-i}`.  This keeps `|x| < 2 a_i` before step `i`
    for any approximation of `x` within about `0.2 \cdot 2^{-i}`, so the
    decisions need no exact comparison: they run on a signed 64-bit
    model of the remainder over windows of 52 indices, refreshed
    between windows from a head of `\lceil (r + 127) / 64 \rceil`
    ordinary limbs to which each window's digits are applied, and never
    branch on the full-precision operand.  About `r/3` digits are
    nonzero (the unsigned greedy reduction uses `r/2`).  All of them
    are then subtracted from the full operand in one pass of vector
    additions, and a negative remainder (in `(-2^{-r-1}, 0)`) receives
    one more `b_{r+1}`.  The reconstruction multiplies `\exp(t)` by the
    factors `1 + d_i 2^{-i}`, each one pass `y \leftarrow y \pm
    (y \gg i)`.

    The full-precision operands are held in radix `2^{56}`: `F = \lceil
    64 n / 56 \rceil + 1` digits of 56 bits in the 64-bit lanes of
    :type:`vec4n` (see :ref:`machine-vectors-vec4n-extra`), in signed
    carry-save form, so that every addition and every shift-and-add is
    a lanewise vector operation without carries.  The headroom of eight
    bits per lane admits 62 table entries per pass of the reduction and
    about 64 factors of the reconstruction between carry normalizations
    (one vector pass each).  The radix trades the digit count against
    that headroom: `2^{52}` and `2^{58}` measured within a few percent,
    `2^{60}` clearly slower (carry passes after every few factors).

    The tables are shared caches (`Shared tables`_) at the largest
    precision and `r` used so far by any thread, smaller requests
    reading the top digits and limbs of the rows: `a_i` and `b_i` in
    radix `2^{56}`, their top limbs for the head, and the window models;
    `b_i` is computed from the cached `a_i` (by `b_i = a_i + b_{2i}`,
    and for `i > r/2` by the direct fixed-point summation of
    `\sum_k 2^{-ik}/k`), so no other precomputation is needed.  At 600
    limbs and `r = 768` the radix tables take 8.5 MB.

    Measured (x86-64, AVX2), against the unsigned reduction with limb
    arithmetic at the best `r` of each: 1.2 times faster at 24 to 32
    limbs, 1.45 to 1.5 times at 48 to 300 limbs.  Since the reduction
    and the reconstruction are much cheaper, the tuned `r` is larger
    (192 at 64 limbs, 768 from 448 limbs), and so is the cost of the
    first call (the tables at that `r`): at 300 limbs about 20 ms
    instead of 10 ms.  Built for a target without native integer
    vectors (the generic backends) the function still breaks even.

.. function:: void _mp_real_log1p_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)

    Sets `(res, n)` to an approximation of `\log(1 + (x, n))` for any
    `0 \le x < 1`, with *err* at most `3 r + 64`, by the dual of the
    bitwise exponential reduction (the L-mode BKM recurrence).  A
    product `P` of factors `1 + 2^{-i}`, `i = 1, \ldots, r`, is built
    up greedily below `X = 1 + x`, each accepted factor costing one
    shift-and-add, and the residual is evaluated as

    .. math::

        \log\frac{X}{P} = 2 \operatorname{atanh}\frac{X - P}{X + P},

    where the numerator is the exact deficit maintained through the
    reduction and the single division fuses the normalization by `P`
    with the atanh transformation, so that the quotient satisfies the
    :func:`_mp_real_atanh_rs` contract directly (and the odd series needs
    half the terms of `\log(1 + u)`).  Finally the tabulated logarithms
    of the used factors are added.  The same cached table serves
    :func:`_mp_real_exp_bitwise_rs` and this function.  The significant
    length of `P` grows by `i` bits per accepted factor, so the
    per-factor work starts out at a single limb.

    Passing `r = 0` selects the per-size kernels for `n \le 10` on
    64-bit machines, as for :func:`_mp_real_exp_bitwise_rs`.  Explicit
    values require `r \ge 32`.

    From ``MP_REAL_LOG1P_BITWISE_VEC_MIN = 48`` limbs on 64-bit machines
    the work is done by :func:`_mp_real_log1p_bitwise_vec`.  (Below,
    the limb reduction above, whose work grows with the significant
    length of `P`, keeps up with it.)

.. function:: void _mp_real_log1p_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n)` to an approximation of `\log(1 + (x, n))` for
    `0 \le x < 1` and `32 \le r \le 64 n - 16`, within the bound of
    :func:`_mp_real_log1p_bitwise_rs` (measured: below 7 ulps), with the
    signed digits and the radix-`2^{56}` vector arithmetic of
    :func:`_mp_real_exp_bitwise_vec`.  64-bit only.

    The digits drive `Y = 1 + x` multiplicatively to 1:
    `Y \prod_{i=1}^{r} (1 + d_i 2^{-i}) = 1 + u` with `|u| < 2^{-r}`.
    With `w = -\log Y` these are exactly the dynamics of the exponential
    reduction, so the same thresholds apply, transferred to
    `u = Y - 1`: `d_i = -1` if `u \ge c + c^2/2`, `d_i = 1` if
    `u \le -(c - c^2/2)`, `c = 0.75 \cdot 2^{-i}` (equal to the
    thresholds on `w` up to `c^3/6`).  No step at `i = 0` is needed,
    since `\log 2 < 2 \log(3/2)`.  The decisions run on a 64-bit model
    of `u` per window of 52 indices, updated by
    `U \leftarrow U + d_i (2^{P-i} + (U \gg i))`, and between windows
    the factors are applied to a head of `Y` in ordinary limbs.  The
    full-precision `Y` then takes all factors in one series of vector
    passes (as in the exponential's reconstruction), and

    .. math::

        \log(1 + x) = \log(1 + u) - \sum_{d_i = 1} a_i + \sum_{d_i = -1} b_i,
        \qquad \log(1 + u) = \pm 2 \operatorname{atanh}\frac{|u|}{2 + u},

    one division and the atanh series for the residual, and one batched
    vector pass for the table sum (the cached tables of
    :func:`_mp_real_exp_bitwise_vec`).

    Measured (x86-64, AVX2) at the tuned `r` of each, against the limb
    reduction: 1.05 to 1.2 times faster from 64 to 600 limbs.  The
    reduction itself becomes 1.5 times faster; the atanh series, which
    the limb reduction already shortens to the same `r`, and the
    division dominate the remaining time.

.. function:: void _mp_real_sin_cos_reduced(nn_ptr ysin, nn_ptr yg, ulong * err, nn_srcptr t, slong n, flint_bitcnt_t r, int alg)

    Sets `(ysin, n)` and `(yg, n)` to approximations of `\sin(t)` and
    `g = 1 - \cos(t)` for a reduced argument `(t, n)` with
    `t < 2^{-r}`, `r \ge 16`, independent of any particular argument
    reduction (algorithms 1 and 2 additionally require `r \ge 32`).
    Both results are pure fractions (`\sin t < 2^{-r}`,
    `g < 2^{-2r-1}`), but both buffers must have room for `n + 1` limbs,
    the top limb being scratch.  *err* is the rigorous bound computed at
    run time for both outputs, at most 96.  Returning `g` rather than
    `\cos` preserves the information near the top: the tangent
    half-angle reconstruction of :func:`_mp_real_sin_cos_bitwise_rs` and
    :func:`_mp_real_tan_bitwise_rs` consumes exactly this pair, with
    `\cos t` never formed explicitly.

    *alg* selects the internal method:

    * 0: the tuned automatic choice;
    * 1: the direct sine and cosine rectangular-splitting series (the
      outputs of :func:`_mp_real_sin_cos_rs` converted in place,
      *err* = 15);
    * 2: the sine series plus a squaring and a square root,
      `g = 1 - \sqrt{1 - \sin^2 t}` (half the series terms);
    * 3: one bit-burst step: the leading slice of *t*, on limb
      boundaries, evaluated as a complex exponential by Gaussian binary
      splitting, and the remainder by the tuned series at the raised
      rate;
    * 4: the full bit-burst algorithm with slice lengths tripling, which
      is asymptotically quasi-optimal for very large `n`.

    The burst evaluates each slice factor
    `\exp(i x_k) = \cos x_k + i \sin x_k` by a JOINT truncated binary
    splitting of the sine and `1 - \cos` series over one shared exact
    denominator per slice (``_mp_real_sin_cos_sum_bs_powtab``).  Past a
    per-slice term threshold the cosine track is dropped from the tree
    and recovered from the slice's own denominator by one square root,

    .. math::

        \cos_k Q = \sqrt{(Q B^{Q_e})^2 - (\sin_k Q B^{Q_e})^2}.

    In ``mp_real_t`` arithmetic, the factors are assembled by ball
    additions, the complex numerator accumulates by one
    ``mp_real_mul_complex`` per slice and the real denominator by one
    product, and the finish is two ball divisions against the single
    accumulated denominator, with no full-precision square root per
    slice.  As in :func:`_mp_real_exp_reduced`, everything is
    limb-granular and dead low limbs of a slice fold into its frame
    (including the internal power tables, whose entries strip their
    dead low limbs into per-entry exponents).

.. function:: void _mp_real_exp_notab(nn_ptr y, ulong * err, nn_srcptr x, slong n)
              void _mp_real_exp_notab_r(nn_ptr y, ulong * err, nn_srcptr x, slong n, int r)

    Sets `(y, n + 1)` (`n` fraction limbs and a unit limb) to an
    approximation of `\exp(x)` for any `(x, n)` in `[0, 1)`, without
    table-based argument reduction.  With `z` the number of leading zero
    bits of `x`, the argument is halved `h = \max(0, r(n) - z)` times (an
    exact shift inside the internal guard limbs),
    :func:`_mp_real_exp_reduced` runs on the halved argument at depth
    `z + h`, and the result is squared `h` times, each squaring one
    ``sqrhigh`` at the working precision.  All intermediate values lie
    in `[1, e)`.  Relative errors double per squaring, so
    `h + \log_2 n + 8` guard bits (rounded up to limbs) keep the
    amplified component below one output ulp; *err* is at most 128.

    The reduction depth `r(n)` comes from a tuned table: `r = 32` and
    `r = 16` alternate through the rectangular-splitting range, and
    `r = 32` in the bit-burst regime.  The second function takes the
    depth `r` explicitly, for tuning.

.. function:: void _mp_real_sin_cos_notab(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n)
              void _mp_real_sin_cos_notab_r(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n, int r)

    Sets `(ysin, n + 1)` and `(ycos, n + 1)` (`n` fraction limbs and a
    unit limb each) to approximations of `\sin(x)` and `\cos(x)` for any
    `(x, n)` in `[0, 1)`, without table-based argument reduction.  As in
    :func:`_mp_real_exp_notab` the argument is halved
    `h = \max(0, r(n) - z)` times and :func:`_mp_real_sin_cos_reduced`
    runs at depth `z + h`; the angle is then doubled back on
    `g = 1 - \cos` alone,

    .. math::

        g(2\theta) = 2 g (2 - g),

    one ``sqrhigh`` per doubling with every intermediate a pure fraction
    (the chain ends at `g(x) \le 1 - \cos 1 < 0.46`, so its inputs stay
    below `g(1/2) < 0.13`).  The final sine comes from one square root,

    .. math::

        \sin x = \sqrt{2g - g^2},

    through ``mpn_sqrtrem`` at small sizes and
    :func:`_mp_real_sqrt_newton` above, the input taken at a limb
    position of matching parity so that the root placement is a limb
    copy.  The absolute error of `g` multiplies by at most 4 per
    doubling and the square root divides it by `2 \sin x`; since any
    path that doubles at all has `x \ge 2^{-z-1}`, the total
    amplification is bounded by `2^{2 r - z + 2}`, and
    `2h + z + \log_2 n + 8` guard bits keep the amplified component
    below one output ulp; *err* is at most 128.

    The tuned depth is `r = 32` across the rectangular-splitting range
    and `r = 24` in the bit-burst regime; the second function takes it
    explicitly.

.. function:: void _mp_real_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n, int r)

    Sets `(ysin, n+1)` and `(ycos, n+1)` to approximations of
    `\sin((x, n))` and `\cos((x, n))` for any `0 \le x < 1`; either
    output may be ``NULL``.  *err* is at most `6 r + 128`.  The greedy
    reduction with the cached angles `A_i = \operatorname{atan}(2^{-i})`
    is applied to `x/2` for `i = 1, \ldots, r`:

    .. math::

        \frac{x}{2} = \sum_{i \in \text{used}} A_i + t', \qquad t' < 2^{-r}.

    The halving is required by the windowed decision model of the
    reduction, which needs table entries below `2^{-i}`: `A_i`
    satisfies this, as `\log(1 + 2^{-i})` does, but the doubled angles
    `2 \operatorname{atan}(2^{-i}) \approx 2^{1-i}` of the underlying
    rotation identity do not.  Concavity of `\operatorname{atan}` gives
    `A_{i-1} < 2 A_i`, so each index is used at most once.

    The used indices drive the rotation `W = \prod (1 + i 2^{-i})`, two
    shifts and two add/subtracts per factor, and everything is then read
    off through the tangent half-angle reconstruction described under
    :func:`_mp_real_tan_bitwise_rs`: with `u = \tan t'`,
    `T = w_y + w_x u` and `D = w_x - w_y u` (whose ratio is
    `\tan(x/2)`, though the quotient is never formed),

    .. math::

        A = T D, \quad B = D^2, \quad C = T^2, \qquad
        \sin x = \frac{2A}{B + C}, \quad \cos x = 1 - \frac{2C}{B + C},

    each output one reciprocal division and one mulhigh.  Every divisor
    is arranged, by conditionally halving `T` and `D` together and then
    shifting by at most two more bits (the exponent folds into the final
    doubling), to be an `n`-limb value with its top bit set.

    Passing `r = 0` selects the per-size kernels for `n \le 10` on
    64-bit machines (as for :func:`_mp_real_exp_bitwise_rs`).  Beyond,
    from ``MP_REAL_TRIG_BITWISE_VEC_MIN`` limbs on 64-bit machines, the
    reduction and the rotation are those of
    :func:`_mp_real_tan_halfangle_vec`, whose residual enters through
    :func:`_mp_real_series_rs_tan` where its tables reach; otherwise
    (and on 32-bit machines) the sine and `1 - \cos` of the residual come
    from :func:`_mp_real_sin_cos_reduced`, as described under
    :func:`_mp_real_tan_bitwise_rs`.  Explicit values require `r \ge 32`.
    The angle table is shared with :func:`_mp_real_atan_bitwise_rs`.

.. function:: void _mp_real_atan_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)

    Sets `(res, n)` to an approximation of `\operatorname{atan}((x, n))`
    for any `0 \le x < 1`, with *err* at most `4 r + 64`, by the
    vectoring dual of the rotation above (as
    :func:`_mp_real_log1p_bitwise_rs` is the dual of
    :func:`_mp_real_exp_bitwise_rs`).  The vector `(X, Y) = (1, x)` is
    rotated towards the real axis by the factors `1 - i 2^{-i}`,
    `i = 1, \ldots, r`, applying a factor -- two shifts and two
    add/subtracts -- whenever it keeps `Y \ge 0`, that is whenever
    `Y \ge X 2^{-i}`.  This is a greedy reduction of the angle with
    steps `A_i`, so afterwards `\operatorname{atan}(Y/X) < 2^{-r}` and a
    single division yields a residual meeting the
    :func:`_mp_real_atan_rs` contract.  Because the angle is scale
    invariant, the growth of the vector needs no compensation.  The
    tabulated angles of the used factors are summed at the end; the
    total is below `\sum_{i \ge 1} A_i \approx 0.898 < 1`, so no
    rescaling is needed.  The decisions of the vectoring loop never
    consult the table (they compare `Y` against `X 2^{-i}`).

    Passing `r = 0` selects the per-size kernels for `n \le 10` on
    64-bit machines, as for :func:`_mp_real_exp_bitwise_rs`.
    Explicit values require `r \ge 32`.

    From ``MP_REAL_ATAN_BITWISE_VEC_MIN = 4`` limbs on 64-bit machines
    the work is done by :func:`_mp_real_atan_bitwise_vec`.

.. function:: void _mp_real_atan_bitwise_vec(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n)` to an approximation of `\operatorname{atan}((x, n))`
    for `0 \le x < 1` and `32 \le r \le 64 n - 16`, within the bound of
    :func:`_mp_real_atan_bitwise_rs` (measured: below 4 ulps), with the
    signed digits and the radix-`2^{56}` vector arithmetic of
    :func:`_mp_real_exp_bitwise_vec`.  64-bit only.

    The vector `(X, Y) = (1, x)` is rotated by factors
    `1 - d_i i 2^{-i}`, `d_i \in \{-1, 0, 1\}`, `i = 1, \ldots, r`,
    `d_i = \pm 1` when `\pm Y \ge 0.75 \cdot 2^{-i} X`.  On the angle
    `\theta` of `(X, Y)` this is the signed-digit reduction with steps
    `A_i` on both sides: it keeps `|\theta| < 2 A_i` before step `i`
    (by `A_i < 2 A_{i+1}`), for any decision threshold in
    `[A_i - 2 A_{i+1}, 2 A_{i+1}]` (the lower end is negative), which
    leaves a tolerance of about `0.23 \cdot 2^{-i}`; so the decisions
    run on 64-bit models of `X` (scaled by `2^{62}`) and `Y` (scaled by
    `2^{P}`) per window of 52 indices, the heads of `X` and `Y` in
    ordinary two's complement limbs being rotated between windows.  The
    full-precision vector is rotated in fused vector passes (both
    components from the old values, one pass per factor), and

    .. math::

        \operatorname{atan}(x) = \sum_{d_i = 1} A_i - \sum_{d_i = -1} A_i
            \pm \operatorname{atan}\frac{|Y|}{X},

    one division, the atan series and one batched vector pass over the
    radix-`2^{56}` copy of the angle table (a shared cache like the
    logarithms).

    Measured (x86-64, AVX2) at the tuned `r` of each, against the limb
    reduction: 1.3 times faster from 16 to 96 limbs, 1.4 to 1.5 times
    from 128 to 600 limbs.

.. function:: void _mp_real_tan_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)

    Sets `(res, n+1)` to an approximation of `\tan((x, n))` for any
    `0 \le x < 1`, with *err* at most `8 r + 256` (for `r = 0` at
    `n \le 10` on 64-bit machines, the bound of :func:`_mp_real_tan_sc`,
    below 1800).  Since `\tan(1) < 1.56` the result carries a unit limb.

    The angle `x/2` is reduced by the greedy rotation of
    :func:`_mp_real_sin_cos_bitwise_rs`, accumulating
    `W = \prod (1 + i 2^{-i})` over the accepted factors.  Since
    `\arg W = \sum A_i`, one has `\tan(\sum A_i) = w_y / w_x`, and the
    growth of `|W|` cancels in the ratio.  With `u = \tan(t')` the
    addition formula gives the half-angle tangent in a single division,

    .. math::

        t = \tan(x/2) = \frac{w_y + w_x u}{w_x - w_y u},

    whose denominator lies near `w_x \in (0.72, 1.17)`, and then

    .. math::

        \sin x = \frac{2t}{1 + t^2}, \quad
        \cos x = \frac{1 - t^2}{1 + t^2}, \quad
        \tan x = \frac{2t}{1 - t^2}.

    `t` itself is never divided out: with `T = w_y + w_x u`,
    `D = w_x - w_y u`, `A = TD`, `B = D^2` and `C = T^2` these are

    .. math::

        \sin x = \frac{2A}{B + C}, \quad
        \cos x = 1 - \frac{2C}{B + C}, \quad
        \tan x = \frac{2A}{B - C},

    one reciprocal division and one mulhigh per output.  In the limb
    path the residual enters through `\sin t'` and `\cos t' = 1 - g`
    without their quotient being formed,

    .. math::

        T = w_y + w_x \sin t' - w_y g, \qquad
        D = w_x - w_x g - w_y \sin t',

    four small multiplications; `\cos t'` cancels in every ratio just as
    `|W|` does.  The pair `(\sin t', g)` comes from
    :func:`_mp_real_sin_cos_reduced`.

    Passing `r = 0` at `n \le 10` on 64-bit machines gives
    :func:`_mp_real_tan_sc` instead.  From
    ``MP_REAL_TRIG_BITWISE_VEC_MIN = 4`` limbs on 64-bit machines
    (beyond those sizes for `r = 0`), the reduction
    and the rotation of this function and of
    :func:`_mp_real_sin_cos_bitwise_rs` are done by
    :func:`_mp_real_tan_halfangle_vec`.

.. function:: void _mp_real_tan_sc(nn_ptr y, ulong * err, nn_srcptr v, slong n)

    Sets `(y, n + 1)` to an approximation of `\tan((v, n))` for
    `0 \le v < 1` and `n \le 32`, from the sine and cosine kernels
    (``_mp_real_sin_cos_opt``: the table-driven kernels up to three
    limbs, the dense-table kernels beyond) and one division
    (``flint_mpn_tdiv_q``; a hardware division at one limb).  With both
    kernel outputs within `E` ulps and `\cos v > \cos 1 > 0.5403`, *err*
    is `\lfloor 19 E / 4 \rfloor + 2`.  The tangent kernel uses it up to
    22 limbs, where it is 1.1 to 1.85 times faster than the half-angle
    reconstruction (beyond, the wide tier of the run-time dense tables
    makes it 0.85 times as fast).  64-bit only; declared in
    ``mp_real/impl.h``.

.. function:: void _mp_real_tan_halfangle_vec(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan, nn_srcptr x, slong n, int r)

    Sets `(ysin, n + 1)`, `(ycos, n + 1)` and `(ytan, n + 1)` (each may
    be ``NULL``) to approximations of `\sin x`, `\cos x`, `\tan x` for
    `(x, n)` in `[0, 1)` and `32 \le r \le 64 n - 16`, within the bounds
    of :func:`_mp_real_sin_cos_bitwise_rs` and
    :func:`_mp_real_tan_bitwise_rs` (measured: below 25 resp. 45 ulps,
    where the limb path reaches 600 resp. 2100), with the signed digits and the
    radix-`2^{56}` vector arithmetic of :func:`_mp_real_exp_bitwise_vec`.
    64-bit only.  Declared in ``mp_real/impl.h``.

    The angle `x/2` (exact in radix `2^{56}`) is reduced as the
    exponential's argument, with the steps `\pm A_i` and the thresholds
    `\pm 0.75 \cdot 2^{-i}` (the tolerance of about `0.23 \cdot 2^{-i}`
    discussed under :func:`_mp_real_atan_bitwise_vec`), to
    `x/2 = \sum d_i A_i + t'`; a negative remainder, above
    `-2 A_{r+2}`, receives one more `A_r`, so that `0 \le t' < 2^{-r}`.
    The residual vector `(1, \tan t')` (from :func:`_mp_real_series_rs_tan`
    where the per-size tangent series applies) or `(\cos t', \sin t')`
    (from :func:`_mp_real_sin_cos_reduced`) is converted to radix
    `2^{56}` and rotated by the factors `1 + d_i i 2^{-i}`, in fused
    vector passes, to the pair `(D, T)` of the half-angle reconstruction
    described under :func:`_mp_real_tan_bitwise_rs`, which finishes as
    there.  (This replaces both the product `W` and the four
    multiplications combining it with the residual.)

    Measured (x86-64, AVX2) at the tuned `r` of each, against the limb
    reduction: sin and cos 1.1 times faster from 16 to 96 limbs, 1.1 to
    1.25 times from 128 to 600 limbs; tan 1.1 to 1.3 times.  The
    reconstruction's reciprocal and products and the residual series
    dominate the remaining time.

.. function:: void _mp_real_log1p_2mexp_ui_bs(nn_ptr res, ulong i, slong n)
              void _mp_real_atan_2mexp_ui_bs(nn_ptr res, ulong i, slong n)

    Sets `(res, n)` to a one-sided fixed-point approximation of
    `\log(1 + 2^{-i})` resp. `\operatorname{atan}(2^{-i})`, `i \ge 1`:
    the floor of the value scaled by `2^{\mathrm{FLINT\_BITS} \cdot n}`,
    or one ulp below it, never above.  These build the entries of the
    cached reduction tables (declared in ``mp_real/impl.h``) by binary
    splitting in ``mp_real_t`` ball arithmetic, which tracks all
    truncation errors and the series tail; the result is read off the
    ball as a lower bound.  Blocks of terms at the leaves are summed
    exactly by nonallocating mpn code, and when the whole series fits in
    one block the exact partial sum is floored by one division.  For
    large `i` the logarithm sums the alternating series
    `2^{-i} \sum (-1)^k 2^{-ik}/(k+1)`, otherwise
    `2\operatorname{atanh}(1/(2^{i+1}+1))`.

The shared logarithm and arctangent tables are built on demand in two
tiers.  The smallest indices (`i \le 6` for the logarithms, `i \le 3`
for the angles) are combined from `\log 2, \log 3, \ldots` resp.
Machin-type arctangent formulas, computed by binary splitting.  Larger
indices come from a fixed-point multi-summation in which one reciprocal
per odd series index serves every table entry at once.  Each cached
entry carries a guard limb below its value limbs, and entries are
one-sided (the exact floor or one ulp below).

Tuning the reduction parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reduction parameter selected by `r = 0` is tuned in two tiers.

**Small sizes** are tuned per (function, size) at arbitrary `r` by
``dev/tune_mp_real.py``.  It generates an out-of-tree source with one
fully specialized candidate per reduction parameter (the same bodies
production uses), builds it against the in-tree library, validates
every candidate against MPFR and times it, and selects the fastest
whose measured error stays within a margin of the documented budget.
With ``--emit`` it writes the production file
``src/mp_real/FUNC_opt_<n>.c``; ``--pin R`` at the shipped `r`
reproduces the shipped file byte for byte.  These are the dense-table
kernels at four to ten limbs (exp, log1p, atan and sin_cos), whose
candidates are `r = 7L - 1` for `L` levels (up to the levels the tables
provide).  The table-driven kernels below four limbs have no free
parameter.

**Large sizes** take `r` from a ladder of 32 and the multiples of 64,
via per-function crossover tables, from interleaved timings of every
ladder value at fixed sizes.  ``src/mp_real/tune/tune-bitwise-r``
regenerates the tables by binary searches for the crossovers (which
need a quiet machine).  The tables end at `r = 512` (exp, sin/cos),
448 (log1p) and 320 (atan).  :func:`_mp_real_sin_cos_bitwise_rs` and
:func:`_mp_real_tan_bitwise_rs` share one table.

On 64-bit machines the vector paths have separate tables, from their
own sweeps (x86-64, AVX2; the minimum of interleaved runs of every
ladder value at 25 to 27 sizes up to 600 limbs):

* :func:`_mp_real_exp_bitwise_vec`: 64 below 24 limbs, 96, 128 from 40,
  192 from 56, 256 from 80, 320 from 128, 448 from 160, 512 from 224,
  640 from 256, 768 from 448;
* :func:`_mp_real_log1p_bitwise_vec` (from 48 limbs): 96, 128 from 90,
  256 from 120, 320 from 240, 448 from 280, 512 from 380;
* :func:`_mp_real_atan_bitwise_vec`: 32 below 18 limbs, 48, 64 from 30,
  96 from 44, 192 from 88, 320 from 200, 384 from 240, 512 from 380;
* :func:`_mp_real_tan_halfangle_vec` (sin, cos and tan, beyond the
  specialized sizes): 48 below 28 limbs, 64, 96 from 32, 128 from 44,
  192 from 90, 256 from 120, 320 from 150, 384 from 230, 448 from 280,
  512 from 380, 640 from 420.

The cheaper reductions move the optimum to larger `r` than for the limb
paths, shortening the series, at the price of larger tables: the first
call at 600 limbs (computing the tables) takes 54 ms for log1p, 40 ms for
atan and 49 ms for sin and cos, against 32, 23 and 37 ms with the limb
paths.
Built without native integer vectors (the generic backends of
``machine_vectors.h``, e.g. for baseline x86-64), the vector paths
measure roughly even with the limb paths: atan 1.05 to 1.2 times
faster, log1p even, sin, cos and tan 0.85 to 1.1 times.

.. function:: int _mp_real_exp_bitwise_rs_default_r(slong n)
              int _mp_real_log1p_bitwise_rs_default_r(slong n)
              int _mp_real_atan_bitwise_rs_default_r(slong n)
              int _mp_real_trig_bitwise_rs_default_r(slong n)

    The reduction parameter that `r = 0` selects at size `n`: that of
    the per-size kernel where one exists (``MP_REAL_*_OPT_R`` in
    ``impl.h``: the depth of the table-driven kernels' remainder at
    `n \le 3`, the compile-time constant of the dense-table kernels at
    `4 \le n \le 10`), the tuned ladder value beyond.

``src/mp_real/profile/p-mp_real FUNC [nmax]`` prints, for
`n = 1, \ldots, 12` and then geometric steps of about `4/3`, the
precision, the selected `r` and per-call timings.  Each size is called
once before timing so that table precomputation stays out of the
measurement, and the timing loop cycles over an array of random inputs
so that the branchy reductions pay their real misprediction costs.
``profile/p-threads [FUNCS] [MINBITS] [MAXBITS] [THREADS] [SECONDS]``
times :func:`mp_real_exp_bits`, :func:`mp_real_sin_cos_bits`,
:func:`mp_real_log_bits` and :func:`mp_real_atan_bits` against arb at
large precisions for each thread count in the list (default
``1,2,4,8``), printing the ratio and each library's scaling over its
one-thread time; ``p-threads reduced ...`` instead times the algorithms
of :func:`_mp_real_exp_reduced` and :func:`_mp_real_sin_cos_reduced` at
several reduction depths per thread count, the data for retuning the
burst thresholds, which were tuned on one thread.

Progress reports
-------------------------------------------------------------------------------

.. function:: void mp_real_set_verbose(int level)

.. function:: int mp_real_get_verbose(void)

    Sets or reads the level of progress reports on ``stderr`` for huge
    computations: 0 is silent (the default, unless the environment
    variable ``FLINT_MP_REAL_VERBOSE`` sets another level, which is read
    once, when nothing has been set), 1 prints at most one line every
    two seconds from the long loops (the merges of the binary
    splittings, the slices of the bit-burst evaluations) plus the end
    of every phase that took at least that long, and 2 prints the start
    and end of every phase of an evaluation above a few thousand limbs
    with its size, thread budget and time, and (on Linux) the resident
    memory of the process and its high-water mark.  Each line carries the time
    since the first message and the id of the calling thread; the
    messages are serialized, so lines from different threads never
    interleave.

Diophantine argument reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions compute the elementary functions on the unit
interval by *diophantine* (multi-prime) argument reduction [HJ2024]_.
For the exponential, with `p_0 = 2, p_1 = 3, \ldots` the first primes,
integers `c_j` are chosen so that

.. math::

    t = x - \sum_j c_j \log p_j

is tiny, and then

.. math::

    \exp(x) = 2^{c_0} \, \frac{p}{q} \, \exp(t), \qquad
    p = \prod_{c_j > 0} p_j^{c_j}, \quad
    q = \prod_{c_j < 0} p_j^{-c_j}.

The work beyond `\exp(t)` is one long-by-short multiplication by `p`,
one long-by-short division by `q`, and the products `p` and `q`
themselves.  For the trigonometric functions the logarithms of primes
are replaced by `\pi/2` and the doubled arguments `2 \arg \pi_j` of the
first nonreal Gaussian primes `\pi_j = 1 + i, 1 + 2i, 2 + 3i, \ldots`,
and the rational correction by a Gaussian rational.

**Precomputed tables.**  The method rests on two kinds of precomputed
data:

* The values `\log p_j` (resp. `\arg \pi_j`) to the working precision,
  kept in a shared cache as exact floors (see
  :func:`_mp_real_log_primes_vec` under `Reduction tables`_ and `Shared
  tables`_).  They are computed together by one
  Machin-type formula per set of primes, a small number of fast
  `\operatorname{atanh}(1/x_j)` (resp. `\operatorname{atan}(1/x_j)`)
  series with large `x_j`, evaluated by Zuniga's series (see
  `Machin terms as Zuniga series`_).  The formulas were found with
  https://github.com/fredrik-johansson/machin.
* A table of integer relations (:func:`_mp_real_rel_table`), independent
  of the precision: row `i` holds coefficients `d_{ij}` with

  .. math::

      \sum_j d_{ij} \log p_j = \varepsilon_i, \qquad
      |\varepsilon_1| > |\varepsilon_2| > \cdots,

  with small `d_{ij}`.  Tables for common numbers of primes are
  precomputed and stored in the library.

Compared to the bitwise reductions, the precomputation is much cheaper
(one Machin-type evaluation per prime instead of a table of `r` series
at the full precision), at the price of a shallower reduction and hence
a slower evaluation.  This is why the bitwise kernels are used up to
600 limbs and the diophantine ones above.  ``profile/p-diophantine``
prints the trade-off across precisions.

**The descent.**  After the free step by `\log 2` (its power is a
shift), each row of the relation table subtracts the nearest integer
multiple of `\varepsilon_i` from the residual, adding the same multiple
of `d_i` to the coefficient vector.  The descent stops when the weight

.. math::

    \sum_{j > 0} |c_j| \log_2 p_j,

a proxy for the size of `p q`, would exceed a budget.  Larger budgets
give deeper reductions and larger products `p, q`, which have about
half the weight in bits each.  The residual is tracked in double
precision and resynchronized every eight rows from a short fixed-point
evaluation of `x - \sum c_j \log p_j`.  Rounding to nearest leaves `t`
of either sign, while the reduced series takes `t \ge 0`; a negative
`t` is lifted by adding one relation with a positive `\varepsilon`,
which leaves `0 < t \le \varepsilon`.

**Prime powers.**  `p` and `q` are built by vector exponentiation on the
bits of the exponent vector,

.. math::

    P(c) = P(\lfloor c / 2 \rfloor)^2 \prod_{j : c_j \text{ odd}} p_j,

one squaring per exponent bit with the odd-subset product a single limb,
so that the whole product costs about one squaring at the final size.
Products beyond the working precision are carried truncated.

.. function:: void _mp_real_exp_diophantine(nn_ptr y, ulong * err, nn_srcptr x, slong n)
              void _mp_real_exp_diophantine_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, slong num_primes, double max_weight)

    Sets `(y, n + 1)` to `\exp(x)` for `(x, n)` in `[0, 1)`, with *err*
    at most 3, the reduced `\exp(t)` coming from
    :func:`_mp_real_exp_reduced`.  The tunable variant takes the number
    of primes (2 to 64) and the weight budget of the descent; the
    default uses 13 primes and a budget equal to the precision in bits.

    The cached logarithms are one-sided floors, so the dot product
    `\sum c_j \log p_j` is short of the exact value by less than
    `\sum |c_j| + 1` working ulps, a relative error of `\exp(t)` of the
    same size.  The multiplication by `p` is exact and the division by
    `q` a floor, and the factor `2^{c_0} p / q < e` amplifies all
    errors by less than 3; one guard limb keeps the total below one
    output ulp.

.. function:: void _mp_real_sin_cos_diophantine(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n)
              void _mp_real_sin_cos_diophantine_tune(nn_ptr ysin, nn_ptr ycos, ulong * err, nn_srcptr x, slong n, slong num_primes, double max_weight)
              void _mp_real_tan_diophantine(nn_ptr res, ulong * err, nn_srcptr x, slong n)
              void _mp_real_tan_diophantine_tune(nn_ptr res, ulong * err, nn_srcptr x, slong n, slong num_primes, double max_weight)

    Set `(ysin, n + 1)`, `(ycos, n + 1)` (either may be ``NULL``) to
    the sine and cosine, respectively `(res, n + 1)` to the tangent, of
    `(x, n)` in `[0, 1)`, with *err* at most 4.  The reduction writes

    .. math::

        x = c_0 \frac{\pi}{2} + \sum_j c_j \, 2 \arg(\pi_j) + t,
        \qquad
        e^{ix} = i^{c_0} e^{it} \frac{A^2}{N},

    for the Gaussian integer `A = \prod \pi_j^{c_j}` (conjugates for
    negative `c_j`), whose norm `N = |A|^2` is a rational integer.  The
    normalization is one reciprocal of an integer and two
    multiplications, and the tangent is the ratio of the two parts of
    `e^{it} A^2` with no normalization at all.  Gaussian products beyond
    the working precision are carried truncated through the high complex
    products of ``mpn_extras``.  The reduced sine and cosine come from
    :func:`_mp_real_sin_cos_reduced`, whose cost dominates and which
    needs deeper reductions than the exponential to be efficient; the
    default therefore uses 32 Gaussian primes and a weight budget of
    four times the precision.

Newton-Taylor inverses and the AGM logarithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Beyond the range of the bitwise reductions, the logarithm and the
arctangent are computed from the forward functions by one Newton-Taylor
step.  From a starting value `t` correct to about `p_0` bits, the
forward function at the full precision `p` gives a residual `w` with
`|w| < 2^{-p_0}` and

.. math::

    f(x) = t \pm g(w)

exactly, `g` being the Taylor series of `\log(1 + w)` or
`\operatorname{atan}(w)`, so that `p / p_0` terms of `g` suffice.  Since
the `k`-th term has `k p_0` leading zero bits, the polynomial costs a
fraction of a multiplication.  The starting value comes from the
bitwise function up to 512 limbs and from the step itself above, so the
whole costs one forward evaluation plus about a tenth of one.

Everything around the forward function and the starting value is
``mp_real_t`` arithmetic: the residual (its cancellation and the
forward function's documented error carried by the balls), the series
with its powers by squaring at the precision each term contributes at,
the Taylor coefficients scaled to integers by the lcm of their
denominators and one ``mp_real_div_ui``, the tail added to the radius
from the ball's bound on `|w|`, and the sum with `t`.  The export
truncates to `n` limbs and checks the radius, which is tiny since the
forward errors enter through `w` only.  ``tune/tune-newton.c`` measures
the overhead over the forward functions and the effect of the number of
terms.

.. function:: void _mp_real_neglog_newton(nn_ptr y, ulong * err, nn_srcptr x, slong n)
              void _mp_real_neglog_newton_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, int forward, slong N)

    Sets `(y, n)` to `-\log(x)` for `(x, n)` in `[1/2, 1)`; *err* is
    the rigorous bound computed at run time, at most 2.  With
    `E = \exp(t)` from :func:`_mp_real_exp_diophantine` at `n + 1` limbs
    and `w = x E - 1`,

    .. math::

        -\log x = t - \log(1 + w).

    The tunable variant takes the forward algorithm (0 = diophantine,
    1 = bitwise, 2 = notab) and the number `N` of terms of the series
    (up to 16; the starting value is then computed at about `1/(N + 1)`
    of the precision), 0 selecting the default.

.. function:: void _mp_real_atan_newton(nn_ptr y, ulong * err, nn_srcptr x, slong n)
              void _mp_real_atan_newton_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, int forward, slong N)

    Sets `(y, n)` to `\operatorname{atan}(x)` for `(x, n)` in `[0, 1)`;
    *err* is the rigorous bound computed at run time, at most 2.  With
    `(s, c) = (\sin t, \cos t)` from :func:`_mp_real_sin_cos_diophantine`
    at `n + 1` limbs,

    .. math::

        \operatorname{atan} x = t + \operatorname{atan}(w), \qquad
        w = \frac{x c - s}{x s + c}.

    The numerator cancels to about the precision the starting value
    leaves, so the denominator and the division are only needed to that
    precision.  The series in `w^2` has up to 13 terms (the starting
    value at about `1/(2N)` of the precision); the parameters are as for
    the logarithm.

.. function:: void _mp_real_neglog_agm(nn_ptr y, ulong * err, nn_srcptr x, slong n)
              void _mp_real_neglog_agm_tune(nn_ptr y, ulong * err, nn_srcptr x, slong n, slong N)

    The same `-\log(x)` (*err* computed at run time, at most 2) by the
    Sasaki-Kanada formula

    .. math::

        \log(1/q) = \frac{\pi}{\operatorname{agm}(\theta_2(q)^2, \theta_3(q)^2)},

    exact for `0 < q < 1`.  With `r = x 2^{-e}` and `q = r^4` (so that
    `\theta_2(q) = 2 r S_2` needs no fourth root),

    .. math::

        -\log x = \frac{\pi/4}{\operatorname{agm}(4 r^2 S_2^2, \theta_3^2)} - e \log 2,
        \qquad
        \theta_3 = 1 + 2 \sum_{n \ge 1} q^{n^2}, \quad
        S_2 = 1 + \sum_{n \ge 1} q^{n(n+1)},

    with `\pi/4` and `\log 2` from the shared caches of
    :func:`_mp_real_const_pi4` and :func:`_mp_real_const_log2`.

    The exponents of `\theta_3` and `S_2` merged in order are
    `\lfloor m^2/4 \rfloor`, consecutive ones differing by
    `\lfloor m/2 \rfloor`, so both series cost one product per term and
    one per two terms for the powers `q^j` (a squaring for even `j`),
    each at the precision the term contributes at; the omitted terms are
    below `2 q^{k}` for the first omitted exponent `k`.  The first AGM
    step needs no square root, since
    `\theta_2^2 \theta_3^2 = (2 r S_2 \theta_3)^2`, and the pair becomes
    the theta functions at `q^{1/2}` by the doubling formulas.

    The shift is `e = \lceil p / (4 (N+1)^2) \rceil` for `N` terms of
    `\theta_3` (default 4).  The AGM spends about `\log_2 e` iterations
    in its slow phase, and each doubling of `N` saves two of them for
    about `3N` more products at decreasing precision.

    The AGM logarithm costs `O(M(n) \log n)`, against the
    `O(M(n) \log^2 n)` of the Newton-Taylor logarithm over the
    diophantine exponential; with warm caches the two cross over near
    `7 \cdot 10^5` limbs.  Without any precomputed tables, i.e. against
    the Newton-Taylor logarithm over :func:`_mp_real_exp_notab`, the
    AGM logarithm is faster from about 20000 limbs.

Sine and cosine of rational multiples of `\pi`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The values `\sin(\pi t/q)` and `\cos(\pi t/q)` are algebraic numbers
of degree at most `q`, and past a few hundred limbs it is cheaper to
refine them as roots of a polynomial than to evaluate the sine and
cosine of `\pi t/q`: the root iteration costs `O(M(n) \log q)`, the
elementary function `O(M(n) \log^2 n)` with a much larger constant.

**The Chebyshev iteration.** For `\varphi^* = \pi R/Q \in (0, \pi/2)`
in lowest terms, `x^* = \cos \varphi^*` is a simple root of the
Chebyshev polynomial `U_{Q-1}`, with `T_Q(x^*) = \varepsilon =
(-1)^R`.  Given an exact approximation `x_0 = \cos \varphi_0`,
`\varphi_0 = \varphi^* + \delta`, the order-doubling (Montgomery
ladder) evaluation of the pair `(T_k(x_0), T_{k+1}(x_0))` -- one
squaring and one multiplication per bit of `Q`, in the form
`V_k = 2 T_k` with `V_{2k} = V_k^2 - 2`,
`V_{2k+1} = V_k V_{k+1} - V_1` -- gives

.. math::

    C = T_Q(x_0) = \varepsilon \cos(Q\delta), \qquad
    A = x_0 T_Q(x_0) - T_{Q+1}(x_0) = \varepsilon s_0 \sin(Q \delta),
    \quad s_0 = \sin \varphi_0

(for odd `Q = 2k + 1`, `A = (V_k - V_{k+1})(V_k + V_{k+1})/4` replaces
the last doubling and the product by `x_0`).  With `a = \varepsilon A`,
`\sigma = \sin(Q\delta) = a / s_0`, `y = \sigma^2 = 1 - C^2` and
`b = 1/Q`, exactly

.. math::

    x^* = x_0 \cos \delta + s_0 \sin \delta
        = x_0 + x_0 G(y) + \frac{a}{Q} F(y),

    G(y) = {}_2F_1\left(-\tfrac{b}{2}, \tfrac{b}{2}; \tfrac12; y\right) - 1,
    \qquad
    F(y) = {}_2F_1\left(\tfrac{1-b}{2}, \tfrac{1+b}{2}; \tfrac32; y\right),

the hypergeometric forms of `\cos(b \operatorname{asin} \sigma)` and
`Q \sin(b \operatorname{asin} \sigma) / \sigma`, valid on asin's
principal branch, `|Q \delta| < \pi/2`.  No division or square root
enters, `s_0` appearing only through `a` and `y`.  The coefficients of
`G` (from `-1/(2Q^2)`) and of `F` (from 1) decrease in magnitude, so
for `|y| \le 1/2` the tails after `y^J` are below `|y|^J / Q^2`
resp. `2 |y|^J`.  In `z = y / (4Q^2)` both series have integer
coefficients over the common denominators `m! (2m-1)!!` resp.
`m! (2m+1)!!` for `m` terms, and share the powers of `z`.

With `x_0` correct to `e` bits, `\sigma \approx Q 2^{-e}` and summing
until the terms fall below the target makes a step of any order `r`:
`x_0` is needed to about `(p + (r - 1) \log_2 Q) / r` bits of the
target `p`.  A step costs the ladder at full precision, about
`(\log_2 Q - \log_2 r) (S + M)(n)` (the first `\log_2 r` doublings of
the short `x_0` are short products) and a few products for the series
(about `3 M(n)` at order 12).  The steps are written in ``mp_real_t``
arithmetic, which carries the rounding errors of the ladder (its
radius grows by about two bits per doubling, paid for in guard bits);
each level's radius is discarded after its midpoint is truncated to
`x_0`, but serves as the rigorous bound `|x_0 - x^*| \le \rho` from
which the step checks the branch condition (sufficiently,
`\rho \le \min(s^2/8, s/Q)` for `s = \sin \varphi^* \ge 2R/Q`, by
Jordan's inequality).  The recursion ends in the kernel below as soon
as the kernel is cheaper at the level's precision, so for large `Q`
the result is the kernel at `n/r` limbs followed by one step.  The
order is 12: 8, 12 and 16 are within a few percent of each other, 6
about 9% and 4 about 17% slower.

Of the targets `\cos(\pi t/q)` (denominator `q`) and
`\sin(\pi t/q) = \cos(\pi (q - 2t)/(2q))` (denominator `2q` in lowest
terms, `q` or `q/2` for even `q`), the one with the smaller
denominator is iterated, its conditioning (`\sin \varphi^*` small near
`\varphi^* = 0`) paid for in guard bits.  The other value is
`\sqrt{1 - x^{*2}}`, about `3 M(n)`.

**The complex alternative.** Newton-Taylor steps for the reciprocal
`Q`-th root `w = e^{-i \varphi^*}` of `\varepsilon`, as
`w = w_0 (1 - u)^{-1/Q}` with `u = 1 - \varepsilon w_0^Q` (the
complex analogue of :func:`mp_real_rroot_ui`), give both values at
once but need a complex squaring (two products) per bit of `Q` and a
complex product by the short `w_0` per set bit.  Written in the same
ball arithmetic (``tune/tune-sin-cos-pi.c``), it measured 1.4 to 2
times slower than the Chebyshev iteration including its square root,
from `Q = 5` up.

**Choice.**  The Chebyshev iteration is used for denominators up to a
bit length measured against the kernel (``tune/tune-sin-cos-pi.c``,
x86-64): none below 48 limbs, 3 bits at 64 limbs, 7 at 128, 9 at 256,
12 at 512, 17 at 1024, 20 at 2048, 25 at 4096, 30 at 8192, 40 at
16384, 43 at 32768, 48 at 65536, and any admissible denominator
(below `2^{\mathrm{FLINT\_BITS} - 8}`) beyond, where the kernel no longer
uses tables.  For comparison, the kernel costs about 16 `M(n)` at 64
limbs, 20 at 256, 40 at 1024, 50 at 4096 and 150 at 65536.

.. function:: void _mp_real_sin_cos_pi_ui_div_ui(nn_ptr ys, nn_ptr yc, ulong * err, ulong p, ulong q, slong n)
              void _mp_real_sin_cos_pi_ui_div_ui_tune(nn_ptr ys, nn_ptr yc, ulong * err, ulong p, ulong q, slong n, int alg, int r)

    Sets `(ys, n + 1)` and `(yc, n + 1)` (either may be ``NULL``) to
    `\sin(\pi p/q)` and `\cos(\pi p/q)`, requiring `q \ge 1` and
    `0 \le p \le q/2` (so that both lie in `[0, 1]`).  The fraction is
    reduced first; `q \le 2` gives the exact values, `q = 3, 4, 6` the
    quadratic ones from :func:`_mp_real_rsqrt_ui_newton` (*err* at most
    4).  For `q = 5, 10` (`(\sqrt 5 \pm 1)/4` and
    `\sqrt{(5 \mp \sqrt 5)/8}`), `q = 12` (`(\sqrt 6 \pm \sqrt 2)/4`)
    and `q = 8` (`\sqrt{2 \pm \sqrt 2}/2`, from one reciprocal square
    root of `2 + \sqrt 2`) the closed forms in ball arithmetic beat the
    kernel from about 12 limbs (24 for `q = 8`) and the Chebyshev
    iteration by 1.4 to 5 times.  Otherwise *err* is at most 4, computed
    at run time, for the Chebyshev iteration and the closed forms, and
    the kernel's bound plus 2 (at most `6r + 130` with `r \le 768`) for
    the kernel.

    The kernel path forms `v = (\pi/4) \cdot 4t/q` (or its complement
    `(\pi/4) \cdot 2(q - 2t)/q` for `4t > q`, the outputs swapping) from
    the cached floor of `\pi/4` by one ``mpn_mul_1`` and one
    ``mpn_divrem_1``, two ulps below the exact value, and calls the
    kernel of :func:`mp_real_sin_cos_bits`.  The tunable variant takes
    the algorithm (*alg* = 0: the tuned choice, 1: the kernel, 2: the
    Chebyshev iteration for `q < 2^{\mathrm{FLINT\_BITS} - 8}`, whose
    recursion still ends in the kernel, 3: the closed forms for
    `q = 5, 8, 10, 12`, the kernel otherwise) and the order *r* (0: the
    default; at most 16, 10 on 32-bit machines).

Reduction tables
-------------------------------------------------------------------------------

These tables are shared with arb (``arb_exp_arf_log_reduction``,
``arb_sin_cos_arf_atan_reduction``, ``arb_log_primes_vec_bsplit``,
``arb_atan_gauss_primes_vec_bsplit``).

.. type:: mp_real_rel_struct

.. function:: const mp_real_rel_struct * _mp_real_rel_table(int gaussian, slong num)
              int _mp_real_rel_table_is_cached(int gaussian, slong num)

    Returns the relation table for the first *num* primes
    (*gaussian* = 0: `\alpha_j = \log p_j`) or Gaussian primes
    (*gaussian* = 1: `\alpha_0 = \pi/2`, `\alpha_j = 2 \arg \pi_j`).
    The rows of ``d`` are integer relations

    .. math::

        \sum_j d_{ij} \alpha_j = \epsilon_i

    with `|\epsilon_i|` decreasing (an extra row starting with
    ``MP_REAL_REL_TERMINATOR`` ends the table), and the structure also
    carries the reciprocals `1/\epsilon_i`, the weights, and the primes.
    Tables for *num* = 2, 4, 6, 8, 10, 12, 13, 16, 20, 24, 32, 40, 48
    are precomputed; any other size up to ``MP_REAL_REL_MAX`` is
    generated on first use (``_arb_log_precompute_reductions``; seconds
    around 48 primes) and cached per thread.  The second function
    returns nonzero when a call would find the table ready.  Tables are
    freed by :func:`flint_cleanup`.

.. type:: mp_real_machin_struct

.. function:: const mp_real_machin_struct * _mp_real_machin_table(int gaussian, slong num)
              slong _mp_real_machin_table_max(int gaussian)
              void _mp_real_machin_get_x(fmpz_t q, const mp_real_machin_struct * tab, slong i)
              void _mp_real_machin_get_c(fmpz_t c, const mp_real_machin_struct * tab, slong i, slong j)
              void _mp_real_machin_get_c_row(fmpz * row, const mp_real_machin_struct * tab, slong i)

    Machin-type sets for the logarithms of the first primes and for the
    arguments of the first nonreal Gaussian primes,

    .. math::

        \log p_i = \frac{1}{\mathrm{den}} \sum_j c_{ij} \operatorname{atanh}\frac{1}{x_j},
        \qquad
        \arg \pi_i = \frac{1}{\mathrm{den}} \sum_j c_{ij} \operatorname{atan}\frac{1}{x_j},

    with square coefficient matrices.  The first function returns the
    best set for *num* values: the largest one with at most *num* terms,
    the remaining values being left to one followup series each, unless
    that would need more followups than a measured limit, in which case
    the next larger set is used.  Sets exist for every size from 4 (3
    for the Gaussian primes) to 32, and for 40 and 48; the second
    function gives the largest.  New sets can be generated with
    https://github.com/fredrik-johansson/machin.

    The arguments are stored as 128-bit values
    (``MP_REAL_MACHIN_X_LIMBS`` limbs each).  The coefficient matrices
    are not stored but reconstructed on first use from the arguments,
    the denominator and ``cbits`` (factoring each argument over the
    primes of the set and inverting the exponent matrix modulo one or
    two primes), and cached per thread.  The accessors return any
    argument or coefficient, or a whole row, as ``fmpz``
    (``_mp_real_machin_c_row_raw`` gives a row as the cached limbs and
    sign bytes).

    In this module the sets are evaluated in ``mp_real_t`` arithmetic
    and combined by exact mpn dot products into the cached logarithms of
    primes and angles of Gaussian primes of the diophantine reductions,
    whose exact floors are read off the balls (the first 13 values of
    each, up to 4608 bits, come from static tables).  The
    `\operatorname{atanh}` and `\operatorname{atan}` terms are
    evaluated by Zuniga's series from a precision of about `b^2/2` limbs
    for arguments of `b` bits, and by :func:`_mp_real_atanh_frac_bsplit`
    resp. :func:`_mp_real_atan_frac_bsplit` below.

.. macro:: MP_REAL_ATAN_GAUSS_MAX

.. var:: const signed char _mp_real_gaussian_primes[]

    The real and imaginary parts, consecutively, of the first
    :macro:`MP_REAL_ATAN_GAUSS_MAX` = 64 nonreal Gaussian primes in order
    of norm: `1 + i`, `1 + 2i`, `2 + 3i`, ...

.. function:: void _mp_real_atan_gauss_vec_arb(arb_ptr res, slong num, slong prec)

    Sets *res* to the angles `\pi/2` and `2 \arg \pi_j`, `j < \mathit{num}`,
    as arb balls of *prec* bits.

Multithreading
-------------------------------------------------------------------------------

The constants, the sets of logarithms and Gaussian-prime arguments,
:func:`mp_real_hypgeom_series` and the bit-burst evaluations behind the
elementary functions (:func:`_mp_real_exp_reduced` and
:func:`_mp_real_sin_cos_reduced`, hence the diophantine and table-free
kernels and everything above them at large precision, including the
Newton iterations for `\log` and `\operatorname{atan}`) use the threads of
:func:`flint_set_num_threads`.  The constants and series are identical
for any number of threads; the elementary functions at large precision
choose between a burst step and a series by the thread count (below),
so their rounding may differ.

* A set of series (the Machin terms of :func:`_mp_real_log_primes_vec`
  and :func:`_mp_real_atan_gauss_vec` together with their followup
  terms, or the logarithms of `\log m` for `\gamma`) runs on up to one
  thread per series.  The threads take the series costliest first
  (about `n / \log_2 x` terms for `1/x`) as they become free, so at
  most one series per thread is in flight, and threads left over when
  there are fewer series than threads go to the splitting of the
  series.  The Machin combination frees each series value as it reads
  it.

* A binary splitting (both sums for `\gamma`,
  :func:`_mp_real_atan_frac_bsplit` and the generic backend) splits the
  thread budget between the halves of a node, as
  :func:`flint_parallel_binary_splitting` does, with one difference.
  The halves run on two threads only if the node's exact values (as
  estimated from the bits per term) stay within ``MP_REAL_PAR_CAP = 4``
  times the working precision.  Above that, each half would hold
  full-size numbers, so both halves run one after the other with the
  whole budget, and the merge runs on two threads: its products fall
  into two independent groups, such as `\{T_1 Q_2, Q_1 Q_2\}` and
  `\{P_1 T_2, P_1 P_2\}`, and each group's multiplications can use its
  half of the threads for the FFT.  A thread started on a subtree gets
  its own leaf buffers and level temporaries, freed when it is done.

* The bit-burst loops use the two strategies of ``arb_exp_arf_bb``.
  The slices of the argument, and the series remainder that ends a
  cascade of single steps, are independent, so they run as tasks on up
  to one thread per task, and the numerator and denominator of the
  result are product trees whose halves fork.  The tasks are listed
  costliest first, and the threads' budgets for the parallelism inside
  a task (the FFT, the forks of the splitting) follow the estimated
  costs, so that a series remainder much costlier than the slices does
  not run on a single thread while the others idle.  A burst step's
  slice parallelizes while the series it replaces does not, so with
  `t` threads the cascade takes a step from about `1/t` of the serial
  crossover (``MP_REAL_EXP_BURST_TERMS``, ``MP_REAL_TRIG_BURST_TERMS``
  terms): the algorithm, and hence the exact rounding of the result,
  depends on the thread count, as in arb.

  Each slice's splitting integers and power table come to several
  times the working precision, and every factor is a full-precision
  number.  The tasks are therefore dealt in advance to one lane per
  thread (greedily by estimated cost), and each lane folds its factors
  into its own running numerator and denominator as they complete:
  the live full-precision numbers are then two per lane plus those
  being computed, independently of the number of slices, and the
  lanes' products are combined by product trees at the end.  The
  assignment depends only on the thread setting, not on timing or on
  how many threads the pool grants, so the result is the same on every
  run.  From ``MP_REAL_EXP_BURST_SERIAL_LIMBS`` (`10^9` bits) for the
  exponential and ``MP_REAL_TRIG_BURST_SERIAL_LIMBS`` (`4 \cdot 10^8`
  bits) for the sine and cosine, the switches of arb's bit-burst
  functions, there is a single lane, and the threads go to the binary splitting
  inside each slice, whose subtrees above a size threshold fork and
  whose merges run their independent products (`\{T_1 Q_2, Q_1 Q_2\}`
  and `\{x^m T_2\}` for the exponential; the two tracks of the
  trigonometric tree) on two threads.

  The slices keep the ratio of the growth of their splitting integers
  (the frame, 64 `D` bits per term, or 128 `D` for the `x^2` of the
  sine and cosine) to the decay of their terms at 2 (3 for the
  tripling ladder of the sine and cosine).  Below a limb of reduction
  (`r < 64`, as in the table-free kernels at `r = 24, 32`) the first
  slice is the top limb alone for that reason: a first slice of two
  (three) limbs has ratio `128 / r` (`384 / 2r`) and was both the
  costliest slice and the peak of the memory.

  The complex products of the trigonometric accumulation share
  transforms (:func:`mp_real_mul_complex`) and so need an FFT buffer
  about twice a real product's; fft_small keeps its buffers per
  thread, so at `10^9` bits the first accumulation raised the peak by
  about 1.2 GB.  Beyond the transform-storage budget of
  ``flint_fft_small_max_transformed_ring_size`` the complex product
  falls back to three real products by itself.

  Memory above the persistent FFT buffers, on the development VM at
  `10^8` bits, in multiples of the size of one number (arb in
  parentheses): sin_cos 20 (15) on one thread, 29 (46) on two.  What
  remains of the difference on one thread is mostly the `m` stored
  powers of the rectangular-splitting series that ends the cascade.
  A smaller block size, `m \approx \sqrt{N/3}` instead of
  `\sqrt{N/1.5}`, would store `\sqrt 2` times fewer (0.7 GB less at
  `10^9` bits) but made that series about 10% slower at `10^8` bits,
  so it is not used.

With the cap, the peak memory stays close to the one-thread peak plus
the FFT's own per-thread tables, which every thread that multiplies at
full size keeps until :func:`flint_cleanup`.

Internals
-------------------------------------------------------------------------------

The declarations used only inside the module (the per-size kernels,
the cached reduction tables and their accessors, the bound machinery of
the ball arithmetic, the from-scratch evaluations of the constants and
the thread helpers) are in ``src/mp_real/impl.h``, which the test,
tuning and profiling programs include as ``"mp_real/impl.h"``.

Shared tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The run-time tables of the module -- the bitwise tables of
`\log(1 + 2^{-i})` and `\operatorname{atan}(2^{-i})` and their
radix-`2^{56}` copies, the run-time levels of the dense tables, the
Machin-type constants `\log p_j` and `2 \arg \pi_j` of the diophantine
reductions, and the constants of :func:`mp_real_const_pi4` and its
siblings -- are process-wide shared caches (``flint_shared_cache_struct``
of ``thread_support.h``; see :ref:`threading`), not thread-local: when
several threads evaluate elementary functions, all of them read the
same tables, which are built once and compete for the cache hierarchy
once.  A table holds one read-only version at a time, built at the
largest precision and index range requested so far by any thread;
since the entries are exact floors, which nest, a request at a smaller
precision reads the top limbs (or radix digits) of the entries, so the
results of a function do not depend on which version served it.  A
thread whose request exceeds the current version builds a larger one
under the cache's lock -- the limb count at least one and a half times
the old version's, the index range as far as requested with the old
entries copied -- and publishes it with one release store; readers
never wait and never see a partial table.  The versions superseded
stay allocated, since other threads may still be reading them, until
:func:`flint_cleanup_master` frees every version of every cache; the
geometric growth bounds their total at a small multiple of the current
version.  The ``_clear`` functions of the tables and
:func:`_mp_real_const_clear_cache` free a cache at once, for the tuning
programs that measure cold builds; they must not run while another
thread may be using the table.  The caches that stay thread-local are
the small ones built per use of a Machin-type set or a relation table
(``machin_tab.c``, ``rel_tab.c``) and the root coefficients
(``root.c``).

None of this is reached at small precision: with the tuned defaults,
the elementary functions up to 8 limbs (on 32-bit machines too, where
the guard limbs of the trigonometric functions at 9 limbs reach the
tuned reduction parameter 64, beyond the static prefixes) and the four
static constants up to 64 limbs touch no run-time table at all,
and on 64-bit machines up to 10 limbs the only run-time storage touched
is the lazily filled dense tables, which the guard limbs of the
functions reach from 9 (two guard limbs: expm1, sinh, cosh, tanh,
asinh, and the nfloat exp10, log10 and log1p) or 10 limbs;
``t-static_regime.c`` checks this for every function of the module and
of ``nfloat``.

Code outside the library (the tests) reads the bitwise tables through
the exported ``_entry`` accessors, which return into the current
version.
