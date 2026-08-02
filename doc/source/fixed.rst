.. _fixed:

**fixed.h** -- fixed-point real arithmetic
===============================================================================

This module provides low-level, low-overhead functions for
fixed-point real arithmetic, intended as kernels for
arbitrary-precision numerical algorithms.

A fixed-point number `(x, n)` is an unsigned `n`-limb fraction
``x[0], ..., x[n-1]`` representing
`\sum_i x_i \, 2^{\mathrm{FLINT\_BITS} (i - n)}`, i.e. `0 \le x < 1`
with unit in the last place (ulp) `2^{-\mathrm{FLINT\_BITS}\, n}`.
Outputs of size `n + 1` additionally carry an integer (units) limb at
index `n`.

This module is mainly optimized for 64-bit systems. With 32-bit limbs, some
generated straight-line and register implementations are disabled and
evaluation goes through generic and fallback code paths.

Arithmetic
-------------------------------------------------------------------------------

Newton-based division and square root
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: void fixed_inv_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n)
              void fixed_inv_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

    Given `(a, an)` with `a_{an-1} \ne 0` representing a fixed-point number
    `a \in [1/B, 1)` with `an` fraction limbs, sets `(q, n+2)` to an
    approximation of `1/a \in (1, B]` with `n` fraction limbs and two
    integral limbs (the highest limb may be zero). The absolute error is
    bounded by `4 B^{-n} / a`. The *newton* suffix flags that the result
    is not ulp-accurate. The basecase divides by ``mpn_tdiv_qr``; the
    main function runs a Newton iteration on middle products, ported
    from :func:`radix_inv_approx`.

.. function:: void fixed_div_newton_invmul(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n)
              void fixed_div_newton(nn_ptr q, nn_srcptr b, slong bn, nn_srcptr a, slong an, slong n)

    Given a numerator `(b, bn)` with `bn \ge 1` fraction limbs representing
    `b \in [0, 1)` and a denominator `(a, an)` with `a_{an-1} \ne 0`
    representing `a \in [1/B, 1)`, sets `(q, n+2)` to an approximation of
    `b/a` with `n` fraction limbs and two integral limbs. The absolute
    error is bounded by `4 B^{-n} / a`. The *invmul* variant multiplies
    the numerator by :func:`fixed_inv_newton`; the main function performs
    a Karp-Markstein iteration, ported from :func:`radix_div_approx`.

.. function:: void fixed_rsqrt_ui_newton_basecase(nn_ptr res, ulong a, slong n)
              void fixed_rsqrt_ui_newton(nn_ptr res, ulong a, slong n)

    Sets `(res, n)` to the fraction limbs of an approximation of
    `1/\sqrt{a}`, requiring `2 \le a < B`. The error is bounded by
    `2 B^{-n}`. Ported from :func:`radix_rsqrt_1_approx`.

.. function:: void fixed_rsqrt_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n)
              void fixed_rsqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

    Given `(a, an)` representing `a \in [B^{-2}, 1)` with `an` fraction
    limbs (at least one of the two highest limbs must be nonzero), sets
    `(q, n+2)` to an approximation of `1/\sqrt{a} \in (1, B]` with `n`
    fraction limbs and two integral limbs. The absolute error is bounded
    by `4 B^{-n} / \sqrt{a}`. Ported from :func:`radix_rsqrt_approx`;
    the basecase combines ``mpn_sqrtrem`` and ``mpn_tdiv_qr``.

.. function:: void fixed_sqrt_newton_rsqrtmul(nn_ptr q, nn_srcptr a, slong an, slong n)
              void fixed_sqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

    Input as for :func:`fixed_rsqrt_newton`; sets `(q, n+2)` to an
    approximation of `\sqrt{a} \in [1/B, 1)` (the computed value can
    round up to 1) with absolute error bounded by `4 B^{-n} / \sqrt{a}`.
    Note that the error is proportional to `1/\sqrt{a}` rather than to
    the output. The main function performs a Karp-Markstein iteration,
    ported from :func:`radix_sqrt_approx`.

Elementary functions
-------------------------------------------------------------------------------

The following routines implement evaluation of elementary functions
following [HJ2024]_ and [Joh2014c]_.

Series evaluation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The series evaluation functions require an argument reduced below
`2^{-32}` (checked with a ``FLINT_ASSERT``) and dispatch internally on
the top limb of `x`: nonzero selects hardcoded straight-line routines
for the 32-bit reduction range, zero selects hardcoded or windowed
general routines exploiting a whole leading zero limb (or more).  All
evaluation uses rectangular splitting.  Error bounds are given in ulp
by the following macros; for :func:`fixed_exp_rs` and the hyperbolic
functions they are one-sided (the computed result never exceeds the
true value), for the alternating functions two-sided.

.. macro:: FIXED_EXP_RS_MAX_ERR(n)
           FIXED_SIN_RS_MAX_ERR(n)
           FIXED_COS_RS_MAX_ERR(n)
           FIXED_SIN_COS_RS_MAX_ERR(n)
           FIXED_SINH_RS_MAX_ERR(n)
           FIXED_COSH_RS_MAX_ERR(n)
           FIXED_SINH_COSH_RS_MAX_ERR(n)
           FIXED_ATAN_RS_MAX_ERR(n)
           FIXED_ATANH_RS_MAX_ERR(n)

    Bounds, in ulp, for the error of the corresponding function
    below.  These are small constants.

.. function:: void fixed_exp_rs(nn_ptr res, nn_srcptr x, slong n)

    Sets `(res, n + 1)` to an approximation of `\exp((x, n))`,
    requiring `x < 2^{-32}`.

.. function:: void fixed_sin_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_cos_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_sin_cos_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n)
              void fixed_sinh_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_cosh_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_sinh_cosh_rs(nn_ptr ysinh, nn_ptr ycosh, nn_srcptr x, slong n)

    Set `(res, n + 1)` to an approximation of the respective function
    of `(x, n)`, requiring `x < 2^{-32}`.  The combined versions allow
    either output pointer to be *NULL*.

.. function:: void fixed_atan_rs(nn_ptr res, nn_srcptr x, slong n)
              void fixed_atanh_rs(nn_ptr res, nn_srcptr x, slong n)

    Sets `(res, n)` to an approximation of `\operatorname{atan}((x, n))`
    resp. `\operatorname{atanh}((x, n))`, requiring `x < 2^{-32}`.

.. function:: void _fixed_exp_rs_fallback(nn_ptr res, nn_srcptr x, slong n)
              void _fixed_sin_cos_rs_fallback(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n, int alternating)
              void _fixed_atan_rs_fallback(nn_ptr res, nn_srcptr x, slong n, int alternating)

    As the corresponding functions above, requiring `x < 2^{-32}` (the number
    of terms is chosen from the actual leading zero bits of `x`). These are
    portable fallbacks which run at constant full precision with coefficients
    generated on the fly; they serve out-of-table sizes, 32-bit machines,
    and the test code.

Bitwise argument reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following functions compute the elementary functions on the unit interval
using bitwise (BKM-style) argument reduction to bring the argument below a
tuned threshold `2^{-r}` followed by series evaluation and reconstruction.
On 32-bit systems, it is assumes that `n \ge 2`.

.. macro:: FIXED_EXP_BITWISE_RS_MAX_ERR(n, r)

    Bound, in ulp, for the error of :func:`fixed_exp_bitwise_rs`,
    currently `9 r + 100`.  The bound grows linearly with `r` because
    all work happens at the output precision; callers wanting sub-ulp
    accuracy should pad the precision by one limb themselves.

.. function:: void fixed_exp_reduced(nn_ptr y, nn_srcptr t, slong wn, flint_bitcnt_t r, int alg)

    Sets `(y, wn + 1)` (``wn`` fraction limbs and a units limb) to an
    approximation of `\exp(t)` for a reduced argument `(t, wn)` with
    `t < 2^{-r}`, `r \ge 16`, independent of any particular argument
    reduction (algorithms 1 and 2 additionally require
    `r \ge 32`).  The error is at most ``FIXED_EXP_REDUCED_MAX_ERR``
    ulps.  *alg* selects the internal method: 0 the tuned automatic
    choice, 1 the direct rectangular-splitting series, 2 the sinh
    series plus a square root, 3 one bit-burst step (the leading
    slice of *t*, on limb boundaries, evaluated by binary splitting
    and the remainder by the sinh series at the doubled rate), 4 the
    full bit-burst algorithm with slice lengths doubling, which is
    asymptotically quasi-optimal for very large ``wn``.  The burst
    machinery works at limb granularity throughout -- slice
    boundaries, splitting-tree truncation frames and quotient
    placements are all limb counts, with no bit shifts on the path
    -- and dead low limbs of a slice fold into its frame, so sparse
    arguments (a few significant limbs over a deep frame, the shape
    of a double-precision input at high precision) shrink the whole
    computation.  The automatic thresholds can be recalibrated with
    ``tune/tune-exp-reduced``.

.. function:: void fixed_exp_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n + 1)` to an approximation of `\exp((x, n))` for any
    `0 \le x < 1`.  Passing `r = 0` selects a
    tuned default from a built-in table of crossovers (the largest
    tabulated value serving all larger `n`); the table can be
    regenerated for a given machine with
    ``src/fixed/tune/tune-bitwise-r.c``.  The argument is reduced below `2^{-r}` (where
    `r \ge 32` is a tuning parameter, internally clamped to
    `\mathrm{FLINT\_BITS} \, n - 16`) by subtracting in turn each
    logarithm `L_i = \log(1 + 2^{-i})`, `i = 0, 1, \ldots, r`, for
    which `L_i \le x`; the Taylor series is evaluated on the reduced
    argument; and the used factors `1 + 2^{-i}` are multiplied back
    in, each by a single shift-and-add.  When the reduced series is
    long, `\sinh` is evaluated instead (half the terms) and the
    exponential reconstructed as
    `\exp(t) = \sinh(t) + \sqrt{1 + \sinh(t)^2}`.

    The logarithm table is generated at runtime and cached per thread
    at the largest precision and index range requested so far.

    Explicit values require `r \ge 32`, the contract of the residual
    series.

    With `r = 0` the sizes `n \le 7` are fully specialized on 64-bit
    machines, one source file per size (``exp_opt_<n>.c``, emitted and
    tuned by ``dev/tune_fixed.py``): the reduction parameter is a
    compile-time constant and the series is built for exactly that
    `r`, with the
    smallest number of terms `N` satisfying
    `N r + \log_2(N!) \ge 64 n`. For `n = 1` and `n = 2` the
    series is hand-written rather than generated.

    Small `r` minimizes table and reduction work, large `r` shortens
    the series; as a rule of thumb on 64-bit machines, `r = 32` is
    best up to about 8 limbs, `r = 64` up to about 32 limbs, and
    `r = 128` or `192` beyond.  See ``profile/p-exp_bitwise_rs.c``.

.. macro:: FIXED_LOG1P_BITWISE_RS_MAX_ERR(n, r)

    Bound, in ulp, for the error of :func:`fixed_log1p_bitwise_rs`,
    currently `3 r + 64`; as for the exponential, all work happens at
    the output precision, so callers wanting sub-ulp accuracy should
    pad the precision by one limb themselves.

.. function:: void fixed_log1p_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n)` to an approximation of `\log(1 + (x, n))` for any
    `0 \le x < 1`, by the dual of the bitwise exponential reduction
    (the L-mode BKM recurrence): a product `P` of factors
    `1 + 2^{-i}`, `i = 1, \ldots, r`, is built up greedily below
    `X = 1 + x`, each accepted factor costing one shift-and-add; the
    residual is evaluated as
    `\log(X/P) = 2 \operatorname{atanh}((X - P)/(X + P))`, where the
    numerator is the exact deficit maintained through the reduction
    and the single division fuses the normalization by `P` with the
    atanh transformation, so that the quotient satisfies the
    :func:`fixed_atanh_rs` contract directly (and the odd series
    needs half the terms of `\log(1 + u)`); finally the tabulated
    logarithms of the used factors are added.  The same cached table
    serves :func:`fixed_exp_bitwise_rs` and this function.

    Passing `r = 0` selects the fully specialized per-size
    implementations for `n \le 7` on 64-bit machines
    (``log1p_opt_<n>.c``, emitted and tuned by ``dev/tune_fixed.py``;
    decisions and updates in straight-line masked carry chains, with a
    compile-time `r` and the atanh series built for it -- the
    hand-written `n \le 2` paths run at `r = 16` with custom series
    evaluation routines.  Explicit values
    require `r \ge 32`, the :func:`fixed_atanh_rs` contract.

    The significant length of `P` grows gradually (by `i` bits per
    accepted factor), so the per-factor work starts out at a single
    limb and `P` is only ever truncated once its exact length exceeds
    `n` limbs.

.. macro:: FIXED_SIN_COS_REDUCED_MAX_ERR

    Bound, in ulp, for the error of each output of
    :func:`fixed_sin_cos_reduced`, currently 96 -- a constant, like
    ``FIXED_EXP_REDUCED_MAX_ERR``, because the internal working
    precision carries guard limbs.

.. function:: void fixed_sin_cos_reduced(nn_ptr ysin, nn_ptr yg, nn_srcptr t, slong wn, flint_bitcnt_t r, int alg)

    Sets `(ysin, wn)` and `(yg, wn)` to approximations of `\sin(t)`
    and `g = 1 - \cos(t)` for a reduced argument `(t, wn)` with
    `t < 2^{-r}`, `r \ge 16`, independent of any particular argument
    reduction (algorithms 1 and 2 additionally require
    `r \ge 32`).  Both results are pure fractions
    (`\sin t < 2^{-r}`, `g < 2^{-2r-1}`), but both buffers must have
    room for `wn + 1` limbs, the top limb being scratch.  The error
    of each output is at most ``FIXED_SIN_COS_REDUCED_MAX_ERR``
    ulps.  Returning `g` rather than `\cos` preserves the
    information near the top: the tangent half-angle reconstruction
    of :func:`fixed_sin_cos_bitwise_rs` and
    :func:`fixed_tan_bitwise_rs` consumes exactly this pair, with
    `\cos t` never formed explicitly.

    *alg* selects the internal method: 0 the tuned automatic choice,
    1 the direct sine and cosine rectangular-splitting series, 2 the
    sine series plus a squaring and a square root
    (`g = 1 - \sqrt{1 - \sin^2 t}`, half the series terms), 3 one
    bit-burst step (the leading slice of *t*, on limb boundaries,
    evaluated as a complex exponential by Gaussian binary splitting
    and the remainder by the tuned series at the raised rate), 4 the
    full bit-burst algorithm with slice lengths tripling, which is
    asymptotically quasi-optimal for very large ``wn``.

    The burst evaluates each slice factor
    `\exp(i x_k) = (\cos x_k + i \sin x_k)` by a JOINT truncated
    binary splitting of the sine and `1 - \cos` series over one
    shared exact denominator per slice
    (``_fixed_sin_cos_sum_bs_powtab``), accumulates the complex
    numerator by windowed middle products (a Karatsuba 3-product
    update when the imaginary content dominates the window slack, a
    borrow-safe direct sum of nonnegative windows otherwise) and the
    real denominator by middle products in ping-pong buffers, and
    finishes with two balanced Newton divisions against the single
    accumulated denominator -- where the classical per-slice scheme
    spends a full-precision square root per slice.  Past a per-slice
    term threshold the cosine track is dropped from the tree and
    recovered from the slice's own denominator by one square root,
    `\cos_k Q = \sqrt{(Q B^{Q_e})^2 - (\sin_k Q B^{Q_e})^2}`.  As
    in :func:`fixed_exp_reduced`, everything is limb-granular and
    dead low limbs of a slice fold into its frame, so sparse
    arguments shrink the whole computation (including the internal
    power tables, whose entries strip their dead low limbs into
    per-entry exponents).  The automatic thresholds can be
    recalibrated with the tuner alongside the exp ones.

.. macro:: FIXED_EXP_NOTAB_MAX_ERR

    Bound, in ulp, for the error of :func:`fixed_exp_notab`,
    currently 128.

.. function:: void fixed_exp_notab(nn_ptr y, nn_srcptr x, slong n)

    Sets `(y, n + 1)` (`n` fraction limbs and a unit limb) to an
    approximation of `\exp(x)` for any `(x, n)` in `[0, 1)`,
    without table-based argument reduction: the number of leading
    zero bits `z` of `x` is inspected, the argument is halved
    `h = \max(0, r(n) - z)` times (an exact shift inside the
    internal guard limbs), :func:`fixed_exp_reduced` runs on the
    halved argument at depth `z + h`, and the result is squared `h`
    times, each squaring one ``sqrhigh`` at the working precision.
    All intermediate values lie in `[1, e)`.  Relative errors double
    per squaring, so `h + \log_2 n + 8` guard bits (rounded up to
    limbs) keep the amplified component below one output ulp; the
    error is at most ``FIXED_EXP_NOTAB_MAX_ERR`` ulps.

    The reduction depth `r(n)` comes from a tuned table: `r = 32`
    and `r = 16` alternate in windows through the
    rectangular-splitting range (the 16-windows are where the
    automatic dispatch inside :func:`fixed_exp_reduced` switches to
    the bit-burst path early and beats the series), and `r = 32`
    carries the whole bit-burst regime -- unlike arb's bit-burst
    exponential (16 squarings below `\sim 10^8` bits), the deeper
    reduction measured consistently ahead on the development
    machine, the shrunken first-level splitting trees outweighing
    the extra squarings.  The internal worker
    ``_fixed_exp_notab_r(y, x, n, r)`` takes the depth explicitly,
    for tuning; ``profile/p-fixed exp_notab [nmax]`` compares
    against :func:`arb_exp` without stopping while arb is ahead
    (the interesting regime is the asymptotic one: arb's table-based
    reduction wins below roughly `2^{19}` bits, after which the
    table-free path here measured `1.1\text{--}1.4\times` faster
    than arb's bit-burst exponential up to `8 \times 10^6` bits).

.. macro:: FIXED_SIN_COS_NOTAB_MAX_ERR

    Bound, in ulp, for the error of each output of
    :func:`fixed_sin_cos_notab`, currently 128.

.. function:: void fixed_sin_cos_notab(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n)

    Sets `(ysin, n + 1)` and `(ycos, n + 1)` (`n` fraction limbs
    and a unit limb each) to approximations of `\sin(x)` and
    `\cos(x)` for any `(x, n)` in `[0, 1)`, without table-based
    argument reduction.  As in :func:`fixed_exp_notab` the argument
    is halved `h = \max(0, r(n) - z)` times and
    :func:`fixed_sin_cos_reduced` runs at depth `z + h`; the angle
    is then doubled back on the cosine alone, working with
    `g = 1 - \cos` throughout:

    .. math::

        g(2\theta) = 2 g (2 - g),

    one ``sqrhigh`` per doubling with every intermediate a pure
    fraction (the chain ends at `g(x) \le 1 - \cos 1 < 0.46`, so
    its inputs stay below `g(1/2) < 0.13`).  The final sine comes
    from one square root, `\sin x = \sqrt{2g - g^2}`, through
    ``mpn_sqrtrem`` at small sizes and :func:`fixed_sqrt_newton`
    above a cutoff, the input taken at a limb position of matching
    parity so that the root placement is a limb copy.  The absolute
    error of `g` multiplies by at most 4 per doubling and the
    square root divides it by `2 \sin x`; since any path that
    doubles at all has `x \ge 2^{-z-1}`, the total amplification is
    bounded by `2^{2 r - z + 2}`, and `2h + z + \log_2 n + 8` guard
    bits keep the amplified component below one output ulp.  The
    error of each output is at most
    ``FIXED_SIN_COS_NOTAB_MAX_ERR`` ulps.

    The tuned depth is `r = 32` across the rectangular-splitting
    range and `r = 24` in the bit-burst regime (matching arb's
    bit-burst sine/cosine); the internal worker
    ``_fixed_sin_cos_notab_r`` takes it explicitly, and
    ``profile/p-fixed sin_cos_notab [nmax]`` compares against
    :func:`arb_sin_cos` (measured `1.2\text{--}1.8\times` faster
    from `6 \times 10^4` through `8 \times 10^6` bits on the
    development machine).

.. macro:: FIXED_SIN_COS_BITWISE_RS_MAX_ERR(n, r)

    Bound, in ulp, for the error of :func:`fixed_sin_cos_bitwise_rs`,
    currently `6 r + 128`; as elsewhere in this module all work
    happens at the output precision, so callers wanting sub-ulp
    accuracy should pad the precision by one limb themselves.

.. function:: void fixed_sin_cos_bitwise_rs(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n, int r)

    Sets `(ysin, n+1)` and `(ycos, n+1)` to approximations of
    `\sin((x, n))` and `\cos((x, n))` for any `0 \le x < 1`; either
    output may be ``NULL``.  The greedy reduction with the cached
    angles `A_i = \operatorname{atan}(2^{-i})` is applied to `x/2`
    for `i = 1, \ldots, r`,

    .. math ::

        x/2 = \sum_{i \in \text{used}} A_i + t', \qquad t' < 2^{-r}.

    The halved argument is not an accident of table sharing with
    :func:`fixed_atan_bitwise_rs`: the windowed decision model of
    the shared reduction requires table entries below `2^{-i}`,
    which `A_i` satisfies -- as `\log(1 + 2^{-i})` does -- but the
    doubled angles `2 \operatorname{atan}(2^{-i}) \approx 2^{1-i}`
    of the underlying rotation identity do not.  Concavity of
    `\operatorname{atan}` gives `A_{i-1} < 2 A_i`, so each index is
    used at most once.

    The used indices drive the rotation
    `W = \prod (1 + i 2^{-i})`, two shifts and two add/subtracts per
    factor (in pure registers for `n \le 7`), and everything is then
    read off through the tangent half-angle reconstruction described
    under :func:`fixed_tan_bitwise_rs`: with `T = w_y + w_x u` and
    `D = w_x - w_y u` (whose ratio is `\tan(x/2)`, though the quotient
    is never formed) and

    .. math ::

        A = T D, \qquad B = D^2, \qquad C = T^2,

    each output is one reciprocal division and one mulhigh,

    .. math ::

        \sin x = \frac{2A}{B + C}, \qquad \cos x = 1 - \frac{2C}{B + C}.

    Every divisor is arranged, by conditionally halving `T` and `D`
    together and then shifting by at most two more bits (the exponent
    folds into the final doubling), to be an `n`-limb value with its
    top bit set, as :func:`mpn_tdiv_qr` wants.  For `n \le 4` the whole
    computation past the reduction runs in registers, except that one
    division.

    Passing `r = 0` selects the fully specialized per-size
    implementations for `n \le 12` on 64-bit machines
    (``trig_opt_<n>.c``, emitted and tuned by ``dev/tune_fixed.py``,
    each with a compile-time `r` and the tangent series built for
    it).  Explicit values require `r \ge 32`, the contract of
    :func:`fixed_sin_cos_rs`, which supplies the residual there.  The
    angle table is shared with :func:`fixed_atan_bitwise_rs`.

.. macro:: FIXED_ATAN_BITWISE_RS_MAX_ERR(n, r)

    Bound, in ulp, for the error of :func:`fixed_atan_bitwise_rs`,
    currently `4 r + 64`.

.. function:: void fixed_atan_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n)` to an approximation of
    `\operatorname{atan}((x, n))` for any `0 \le x < 1`, by the
    vectoring dual of the rotation above (as
    :func:`fixed_log1p_bitwise_rs` is the dual of
    :func:`fixed_exp_bitwise_rs`).  The vector `(X, Y) = (1, x)` is
    rotated towards the real axis by the factors `1 - i 2^{-i}`,
    `i = 1, \ldots, r`, applying a factor -- two shifts and two
    add/subtracts on the unrotated components -- whenever it keeps
    `Y \ge 0`, that is whenever `Y \ge X 2^{-i}`.  This is a greedy
    on the angle with steps `A_i = \operatorname{atan}(2^{-i})`, so
    afterwards `\operatorname{atan}(Y/X) < 2^{-r}` and a single
    division yields a residual meeting the :func:`fixed_atan_rs`
    contract.  Because the angle is scale invariant, the growth of
    `|Z|` needs no compensation at all.  The tabulated angles of the
    used factors -- the same `A_i` cached for
    :func:`fixed_sin_cos_bitwise_rs`, which uses them against `x/2`
    -- are simply summed at the end; the total is below
    `\sum_{i \ge 1} A_i \approx 0.898 < 1`, so no rescaling is
    needed.

    Note that the decisions of the vectoring loop never consult the
    table (they compare `Y` against `X 2^{-i}`), so only the
    rotation direction constrains how the shared angles are scaled.

    Passing `r = 0` selects the fully specialized per-size
    implementations for `n \le 7` on 64-bit machines
    (``atan_opt_<n>.c``, emitted and tuned by ``dev/tune_fixed.py``:
    the vectoring in straight-line masked borrow chains, a
    compile-time `r`, and the series built for it).  Explicit values
    require `r \ge 32`.

.. macro:: FIXED_TAN_BITWISE_RS_MAX_ERR(n, r)

    Bound, in ulp, for the error of :func:`fixed_tan_bitwise_rs`,
    currently `8 r + 256`.

.. function:: void fixed_tan_bitwise_rs(nn_ptr res, nn_srcptr x, slong n, int r)

    Sets `(res, n+1)` to an approximation of `\tan((x, n))` for any
    `0 \le x < 1`.  Since `\tan(1) < 1.56` the result carries a unit
    limb.

    Shares the tangent half-angle path of
    :func:`fixed_sin_cos_bitwise_rs`: the angle `x/2` is reduced by
    the same greedy rotation, accumulating
    `W = \prod (1 + i 2^{-i})` over the accepted factors.  Since
    `\arg W = \sum A_i`, one has `\tan(\sum A_i) = w_y / w_x`, and
    because the tangent is a ratio the growth of `|W|` cancels.
    With `u = \tan(t')` the addition formula
    gives the half-angle tangent in a single division,

    .. math ::

        t = \tan(x/2) = \frac{w_y + w_x u}{w_x - w_y u},

    whose denominator lies near `w_x \in (0.72, 1.17)`, and then

    .. math ::

        \sin x = \frac{2t}{1 + t^2}, \quad
        \cos x = \frac{1 - t^2}{1 + t^2}, \quad
        \tan x = \frac{2t}{1 - t^2},

    where `t` itself is never divided out: with `T = w_y + w_x u`,
    `D = w_x - w_y u`, `A = TD`, `B = D^2` and `C = T^2` these are

    .. math ::

        \sin x = \frac{2A}{B + C}, \quad
        \cos x = 1 - \frac{2C}{B + C}, \quad
        \tan x = \frac{2A}{B - C},

    one reciprocal division and one mulhigh per output, and each of
    the three functions is obtained separately.  For sizes without a
    tabulated tangent series the residual contributes through `\sin`
    and `\cos` of `t'` without their quotient being formed: writing
    `\cos t' = 1 - g`,

    .. math ::

        T = w_y + w_x \sin t' - w_y g, \qquad
        D = w_x - w_x g - w_y \sin t',

    four small multiplications, and `\cos t'` cancels in every ratio
    just as `|W|` does.  The pair `(\sin t', g)` itself comes from
    :func:`fixed_sin_cos_reduced`, which owns the tuned choice
    between the direct series, the sine-plus-square-root variant and
    the bit-burst modes for large `n / r`.

.. function:: void fixed_log1p_2mexp_ui_bs(nn_ptr res, ulong i, slong n)
              void fixed_atan_2mexp_ui_bs(nn_ptr res, ulong i, slong n)

    Sets `(res, n)` to a one-sided fixed-point approximation of
    `\log(1 + 2^{-i})` resp. `\operatorname{atan}(2^{-i})`, `i \ge 1`:
    at most the true value, short by no more than a couple of ulps.
    These build the entries of the cached reduction tables (which call
    them with one guard limb) by binary splitting in mpn arithmetic:
    the split products are carried as truncated mantissas with
    limb-radix exponents, an iterative basecase accumulates blocks of
    terms with shifts and single-limb multiplications, and every
    truncation rounds numerators down and denominators up so that the
    final truncating division stays one-sided.  The logarithm sums
    `2\operatorname{atanh}(1/(2^{i+1}+1))`, carrying one `q^2` factor
    per TERM in the materialized denominator so that ranges compose
    without any `q`-power at the merges; the `q^2` multiplications
    all happen at the leaves, each a single ``mpn_mul_1`` whenever
    `q^2` fits in a limb.  At high precision the direct `\log(1+x)`
    series takes over, its larger term count offset by purely dyadic
    merges.

    The smallest indices of the tables (`i \le 6` for the logarithms,
    `i \le 3` for the angles) are instead combined from
    `\log 2, \log 3, \ldots` resp. Gauss-machin style atans, computed
    by generic binary splitting.


The reduction parameter selected by `r = 0` is tuned in two tiers.

**Small sizes** are tuned per (function, size) at arbitrary `r` by
``dev/tune_fixed.py``: it generates an out-of-tree source with one
fully specialized candidate per reduction parameter -- the same
bodies production uses -- builds it against the in-tree library,
validates every candidate against MPFR and times it, and selects the
fastest whose measured error stays within a margin of the documented
budget (near-ties resolve to the cleanest error, since sweep-sized
samples underestimate the maximum).  With ``--emit`` it writes the
production file ``src/fixed/FUNC_opt_<n>.c``; ``--pin R`` at the
shipped `r` reproduces the shipped file byte for byte.

**Large sizes** take `r` from a ladder of 32 and the multiples of 64,
where the general series evaluation is available in the library, via
per-function crossover tables.  ``src/fixed/tune/tune-bitwise-r``
regenerates the tables: for each consecutive pair of ladder values it
binary-searches the smallest `n` at which the larger parameter stops
losing (the optimum is nondecreasing in `n` to within noise), warming
the shared angle/logarithm tables before every measurement, and
prints the ``r_tab``/``n_tab`` pairs to paste into the dispatch
files.  The shipped tables run to `r = 768`.
:func:`fixed_sin_cos_bitwise_rs` and :func:`fixed_tan_bitwise_rs`
share the half-angle path and hence a single table, tuned on the
sine/cosine call with the tangent run printed as a cross-check.

The selection is queryable:

.. function:: int fixed_exp_bitwise_rs_default_r(slong n)
              int fixed_log1p_bitwise_rs_default_r(slong n)
              int fixed_atan_bitwise_rs_default_r(slong n)
              int fixed_trig_bitwise_rs_default_r(slong n)

    The reduction parameter that `r = 0` selects at size `n`: the
    compile-time constant of the specialized per-size implementation
    where one exists, the tuned ladder value beyond.

The shared logarithm and arctangent tables are built on demand in
two tiers: binary splitting for small `i` (through arb for now) and
a fixed-point multi-summation for large `i`, in which one reciprocal
per odd series index serves every table entry at once; each cached
entry carries a guard limb below its value limbs, and entries are
one-sided (the exact floor or one ulp below).

``src/fixed/profile/p-fixed FUNC [nmax]`` prints, for
`n = 1, \ldots, 12` and then geometric steps of about `4/3`, the
precision in bits and digits, the selected `r`, the per-call times of
arb and of the fixed function, and the speedup ratio, stopping at
``nmax`` or once arb has been at least as fast twice in a row.  Each
size is called once before timing so that table precomputation stays
out of the measurement, and the timing loop cycles over an array of
random inputs so that the branchy reductions pay their real
misprediction costs.

Verified constants
-------------------------------------------------------------------------------

.. function:: void fixed_const_pi_div_4(nn_ptr y, slong n)
              void fixed_const_log2(nn_ptr y, slong n)

    Set `(y, n)` to EXACTLY `\lfloor c \, B^n \rfloor` for
    `c = \pi/4` respectively `c = \log 2` -- a verified correct
    floor truncation, not merely an approximation within some ulp
    budget, in the same sense as the cached
    `\log(1 + 2^{-i})` and `\operatorname{atan}(2^{-i})` table
    entries.  Computed limbs are cached per thread and extended on
    demand: the cache stores the floor at the largest size requested
    so far (rounded up geometrically), and since floors nest, any
    shorter request is served as the top limbs of the cached entry.
    Each extension evaluates the constant through binary splitting
    in ``fball`` ball arithmetic (the Chudnovsky series for `\pi`,
    the hypergeometric `\log 2` series of [Zun2025]_ with `11.9`
    bits per term, the fastest known)
    with increasing guard precision until the ball's rigorous radius
    determines the floor uniquely, converted directly from the fball
    representation (limb-aligned mantissa split at the output grid,
    the radius checked to fit strictly inside the tail on both
    sides) without passing through arb.  The reconstruction scalars
    -- `D = 640320^2/12` with the extra `1/4` for `\pi/4`, and the
    denominator `2160` for `\log 2` -- are baked into `q(0)` at the
    leftmost leaf of the splitting tree, which scales the root `Q`
    while leaving `T` invariant, so no final scalar multiplication
    or bit shift remains.  The caches can be freed
    explicitly with ``_fixed_const_pi_div_4_clear`` /
    ``_fixed_const_log2_clear`` and are released by
    :func:`flint_cleanup`.

fball: semi-private ball arithmetic
-------------------------------------------------------------------------------

``fball`` is an internal mpf-like floating-point type with arb-like
ball semantics in the radix `B = 2^{\mathrm{FLINT\_BITS}}`,
intended as a backend for algorithms on scaled or growing numbers
(binary splitting summation, AGM iterations) where all cheap
operations should stay in mpn arithmetic and exponents are
limb-aligned so that shifts are limb copies.  It is DOCUMENTED BUT
SEMI-PRIVATE: the interface lives in ``fixed.h`` for use inside
FLINT and by expert callers, but makes no stability promises across
releases -- the representation, the error-normalization policy and
the set of operations are all subject to future revision.

A ball is `(-1)^{\mathrm{negative}} (d, \mathrm{size})
B^{\mathrm{exp} - \mathrm{size}}` with a top-normalized mantissa
(``size == 0`` encodes zero), plus a rigorous radius of ``err`` ulps
at the anchor `B^{\mathrm{exp} - \mathrm{size} + \mathrm{erra}}`:
the anchor offset ``erra`` keeps `1 \le \mathrm{err} < 2^{128}`
whenever nonzero, so radii far below one mantissa ulp (balls like
`1 - \epsilon`) neither underflow the double nor get flattened to a
whole ulp.  All bound composition rounds up, anchored pairs
`v \, B^a` absorb the exponent ranges a plain double cannot, every
compound bound is multiplied by a fudge factor covering the bound
arithmetic's own rounding, and normalization truncates low mantissa
limbs once the radius exceeds a threshold sitting a couple of limbs
above one ulp, so mantissa length tracks the number of accurate
limbs (with content deliberately retained slightly below the noise
floor, since interval radii over correlated errors grow faster than
true errors).

Operations take a precision `n` in limbs, the caller supplying 2--4
guard limbs; outputs are truncated to about `n` limbs and less when
an operand is accurate to fewer.  Division and the square roots
require operand balls of relative radius below `2^{-30}` (checked,
not assumed: wider divisors are a usage error whose mantissa would
be pure noise).  ``fball_get_fixed`` writes a ball known to lie in
`[0, 1)` as a truncated fixed-point fraction and returns a rigorous
bound in output ulps; ``fball_get_fixed_floor`` is its
verified-floor sibling, writing `\lfloor x B^n \rfloor` when the
radius determines that floor uniquely and reporting failure
otherwise; ``fball_get_arb`` converts losslessly to an arb ball.  ``fball_const_pi_chudnovsky`` and ``fball_const_log2``
evaluate the constants by binary splitting over exact-then-truncated
fballs: blocks of terms accumulate iteratively over exact mpn
integers with a backward recurrence (whole-limb factors on every
word size), the tree keeps exact integers until they outgrow the
target precision and balls afterwards, and the right spine skips
its unused `P` products.  A factored-`Q` variant of the `\pi`
splitting (powers of the constant `C` kept implicit) was measured
11--17% slower at every size and has been dropped.  See ``fixed.h``
for the full interface and per-function contracts.
