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
    radius below `2^{-30}` (checked).  The exponent is aligned, the
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

Conversions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: void mp_real_get_arb(arb_t res, const mp_real_t x)

    Sets *res* to the ball *x*, losslessly.

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

    The constants `c = \pi/4`, `\log 2`, Euler's constant `\gamma`,
    `e`, `\log 10`, Catalan's constant `G`, `\zeta(3)`, `\zeta(5)`,
    `\Gamma(1/3)`, `\Gamma(1/4)` and `2/\pi` (for argument reduction
    by `\pi/2`).

    The ball versions set *res* to a ball for `c` accurate to about `n`
    limbs.  The limb versions set *res* to EXACTLY
    `\lfloor c B^n \rfloor` -- `n` limbs for the constants below 1
    (`\pi/4`, `\log 2`, `\gamma`, `G`, `2/\pi`), `n + 1` limbs with a
    units limb for the others -- and set *err* (if not ``NULL``) to 1, a
    bound on `c - \lfloor c B^n \rfloor B^{-n}` in ulps.  A floor is
    computed from a ball at a few guard limbs, retried with more guard
    limbs until the rigorous radius determines it uniquely.

    With *cache* = 0 the value is computed from scratch and not kept.
    With *cache* = 1 it comes from a per-thread cache holding the floor
    of each constant at the largest precision computed so far.  Floors
    nest (the top limbs of `\lfloor c B^N \rfloor` are
    `\lfloor c B^n \rfloor` for `n \le N`), so a shorter request is a
    copy.  A longer one recomputes the entry at `\max(n + 5, 1.5 N)`
    limbs, so that a sequence of growing requests costs a constant
    factor over the last one.  The ball version with *cache* = 1
    encloses the cached floor with a radius of one ulp.  Internally,
    `\pi/4`, `\log 2` and `2/\pi` are read in place without copying:
    from static 64-limb tables (generated by
    ``dev/gen_mp_real_elem_tables.py``) up to 64 limbs, and from the
    cache beyond.

.. function:: void _mp_real_const_clear_cache(void)

    Frees this thread's cache of constants.  The cache is also freed by
    :func:`flint_cleanup`.

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

**π and log 2** are hypergeometric series summed by
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

    The kernel is :func:`_mp_real_sin_cos_bitwise_rs` up to 600 limbs,
    :func:`_mp_real_sin_cos_diophantine` (32 primes) up to 65536 limbs
    and :func:`_mp_real_sin_cos_notab` above.

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

    with the bitwise kernel resp. :func:`_mp_real_neglog_newton`,
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
    (``*_opt_<n>.c``) directly, take their guard bits from the per-size
    reduction parameters (``MP_REAL_*_OPT_R`` in ``impl.h``) and read
    the argument's limbs in place.

    The radius bounds use double arithmetic that never leaves the
    normal range: scalings by powers of `B` go through
    ``d_mul_2exp_inrange`` with clamped exponents, so no call raises
    ``FE_UNDERFLOW``.

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

Series evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The series evaluation functions require an argument reduced below
`2^{-32}` (checked with a ``FLINT_ASSERT``) and dispatch internally on
the top limbs of `x`: hardcoded straight-line routines (64-bit) for
`2^{-64} \le x < 2^{-32}` and `2^{-128} \le x < 2^{-64}` at small sizes,
and otherwise the tapered rectangular splitting of
:func:`_mp_real_series_rs` at the argument's actual number of leading
zero bits.  The bounds written to *err* are 10 ulps for the exponential
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
    The logarithm table is generated at run time and cached per thread
    at the largest precision and index range requested so far.

    Passing `r = 0` selects a tuned default (see `Tuning the reduction
    parameter`_); explicit values require `r \ge 32`, the contract of
    the residual series.  Small `r` minimizes table and reduction work,
    large `r` shortens the series; as a rule of thumb on 64-bit
    machines, `r = 32` is best up to about 8 limbs, `r = 64` up to
    about 32 limbs, and `r = 128` or `192` beyond.

    With `r = 0` the sizes `n \le 7` are fully specialized on 64-bit
    machines, one source file per size (``exp_opt_<n>.c``, emitted and
    tuned by ``dev/tune_mp_real.py``): the reduction parameter is a
    compile-time constant and the series is built for exactly that `r`,
    with the smallest number of terms `N` satisfying
    `N r + \log_2(N!) \ge 64 n`.

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

    Passing `r = 0` selects the fully specialized per-size
    implementations for `n \le 7` on 64-bit machines
    (``log1p_opt_<n>.c``: decisions and updates in straight-line masked
    carry chains, with a compile-time `r` and the atanh series built for
    it).  Explicit values require `r \ge 32`.

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
    doubling), to be an `n`-limb value with its top bit set.  For
    `n \le 4` the whole computation past the reduction runs in
    registers, except that one division.

    The residual's tangent `u` comes from the per-size kernels
    (``trig_opt_<n>.c``, `n \le 12`, selected by `r = 0` on 64-bit
    machines, each with a compile-time `r` and the tangent series built
    for it), then from :func:`_mp_real_series_rs_tan` up to 600 limbs.
    Beyond that (or for an explicit `r` whose series outgrows its
    tables), and on 32-bit machines, the sine and `1 - \cos` of the
    residual come from :func:`_mp_real_sin_cos_reduced`, as described
    under :func:`_mp_real_tan_bitwise_rs`.  Explicit values require
    `r \ge 32`.  The angle table is shared with
    :func:`_mp_real_atan_bitwise_rs`.

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

    Passing `r = 0` selects the fully specialized per-size
    implementations for `n \le 7` on 64-bit machines
    (``atan_opt_<n>.c``: the vectoring in straight-line masked borrow
    chains, a compile-time `r`, and the series built for it).
    Explicit values require `r \ge 32`.

.. function:: void _mp_real_tan_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)

    Sets `(res, n+1)` to an approximation of `\tan((x, n))` for any
    `0 \le x < 1`, with *err* at most `8 r + 256`.  Since
    `\tan(1) < 1.56` the result carries a unit limb.

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

    one reciprocal division and one mulhigh per output.  For sizes
    without a tabulated tangent series the residual enters through
    `\sin t'` and `\cos t' = 1 - g` without their quotient being formed,

    .. math::

        T = w_y + w_x \sin t' - w_y g, \qquad
        D = w_x - w_x g - w_y \sin t',

    four small multiplications; `\cos t'` cancels in every ratio just as
    `|W|` does.  The pair `(\sin t', g)` comes from
    :func:`_mp_real_sin_cos_reduced`.

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
reproduces the shipped file byte for byte.

**Large sizes** take `r` from a ladder of 32 and the multiples of 64,
via per-function crossover tables, from interleaved timings of every
ladder value at fixed sizes.  ``src/mp_real/tune/tune-bitwise-r``
regenerates the tables by binary searches for the crossovers (which
need a quiet machine).  The tables end at `r = 512` (exp, sin/cos),
448 (log1p) and 320 (atan).  :func:`_mp_real_sin_cos_bitwise_rs` and
:func:`_mp_real_tan_bitwise_rs` share one table.

.. function:: int _mp_real_exp_bitwise_rs_default_r(slong n)
              int _mp_real_log1p_bitwise_rs_default_r(slong n)
              int _mp_real_atan_bitwise_rs_default_r(slong n)
              int _mp_real_trig_bitwise_rs_default_r(slong n)

    The reduction parameter that `r = 0` selects at size `n`: the
    compile-time constant of the specialized per-size implementation
    where one exists, the tuned ladder value beyond (for the
    trigonometric functions, a separate table over the range of
    :func:`_mp_real_series_rs_tan`).

``src/mp_real/profile/p-mp_real FUNC [nmax]`` prints, for
`n = 1, \ldots, 12` and then geometric steps of about `4/3`, the
precision, the selected `r` and per-call timings.  Each size is called
once before timing so that table precomputation stays out of the
measurement, and the timing loop cycles over an array of random inputs
so that the branchy reductions pay their real misprediction costs.

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
  cached per thread as exact floors (see :func:`_mp_real_log_primes_vec`
  under `Reduction tables`_).  They are computed together by one
  Machin-type formula per set of primes, a small number of fast
  `\operatorname{atanh}(1/x_j)` (resp. `\operatorname{atan}(1/x_j)`)
  series with large `x_j`, the atanh terms by Zuniga's series (see
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

    with `\pi/4` and `\log 2` from the per-thread caches of
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
    `\operatorname{atan}` terms are evaluated by
    :func:`_mp_real_atan_frac_bsplit`, and the `\operatorname{atanh}`
    terms of the logarithms by Zuniga's series from a precision of about
    `b^2/2` limbs for arguments of `b` bits, by
    :func:`_mp_real_atanh_frac_bsplit` below.

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

The constants, the sets of logarithms and Gaussian-prime arguments and
:func:`mp_real_hypgeom_series` use the threads of
:func:`flint_set_num_threads`.  The results are identical for any number
of threads.

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

With the cap, the peak memory stays close to the one-thread peak plus
the FFT's own per-thread tables, which every thread that multiplies at
full size keeps until :func:`flint_cleanup`.

Internals
-------------------------------------------------------------------------------

The declarations used only inside the module (the per-size kernels,
the cached reduction tables and their accessors, the bound machinery of
the ball arithmetic, the from-scratch evaluations of the constants and
the thread helpers) are in ``src/mp_real/impl.h``, which the test,
tuning and profiling programs include as ``"mp_real/impl.h"``.  The
per-thread tables cannot be exported from a Windows DLL, so code outside
the library reads them through their exported ``_entry`` accessors.
