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
    is not ulp-accurate. The basecase computes the truncated reciprocal of
    the top `\min(an, n + 1)` limbs of `a` with
    :func:`_flint_mpn_inv_basecase`; the
    main function runs a Newton iteration on middle products, ported
    from :func:`radix_inv_approx`, each step taking either the
    second-order form `T + T(1 - aT)` from `m \approx n/2` limbs or the
    third-order form `T(1 - e + e^2)`, `e = aT - 1`, from
    `m \approx n/3` limbs (the band of `aT` merely wider, `e^2` from
    its top `n - 2m` limbs, and one middle product for the correction).
    The second-order step is already lean (about one full multiplication
    at `n` limbs, the residual being a single middle product), and a
    third-order step costs 1.1 to 1.5 times as much for recursion over a
    third of the precision instead of half, so the order is chosen by
    precision: third-order for `24 \le n \le 80` and `n \ge 1000`.
    Against the pure second-order iteration this measures 11--18% faster
    from 1000 to 16000 limbs, 2--3% beyond, and within noise elsewhere;
    both forms measure worst errors at half the bound.

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
    `2 B^{-n}`.  The main function uses the third-order steps of
    :func:`fixed_rsqrt_newton`, with the residual `a y^2 - 1` exact and
    short (the low limbs of `y^2` times `a`).

.. function:: void fixed_rsqrt_newton_basecase(nn_ptr q, nn_srcptr a, slong an, slong n)
              void fixed_rsqrt_newton(nn_ptr q, nn_srcptr a, slong an, slong n)

    Given `(a, an)` representing `a \in [B^{-2}, 1)` with `an` fraction
    limbs (at least one of the two highest limbs must be nonzero), sets
    `(q, n+2)` to an approximation of `1/\sqrt{a} \in (1, B]` with `n`
    fraction limbs and two integral limbs. The absolute error is bounded
    by `4 B^{-n} / \sqrt{a}`.  The basecase combines
    ``flint_mpn_sqrtrem`` and ``mpn_tdiv_qr``.

    The main function performs third-order Newton steps: from
    `T \approx a^{-1/2}` at `m \approx n/3` limbs (recursively),
    `u = 1 - a T^2` is formed as the band of the product at weights
    `B^{-(n+2)}` to `B^{-(m-2)}` with the unit wrapped into a control
    limb (the construction of :func:`fixed_inv_newton`; a short `a` can
    start the band below the product, whose missing limbs are zero),
    `u^2` from its top `n - 2m` limbs, the correction
    `w = (u \pm 3u^2/4)/2` in units two limbs below the output before
    those guard limbs are dropped, and `q = T \pm T w` by a middle
    product.  The tail `(5/16) u^3 T` is below the rounding for
    `3m \ge n + 4`; worst errors measure at half the bound.  Against
    the second-order iteration ported from :func:`radix_rsqrt_approx`
    (which this replaces) it is 8--34% faster at every size from 32 to
    262144 limbs: the recursion covers a third of the precision instead
    of half, for a step of about the same cost.  :func:`fixed_sqrt_newton`
    takes its reciprocal root from it, and so does ``flint_mpn_sqrtrem``
    in its Newton regime (8--19% faster from 4096 limbs on).

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
    boundaries and splitting-tree truncation frames are limb counts,
    with no bit shifts on the path -- and dead low limbs of a slice
    fold into its frame, so sparse arguments (a few significant limbs
    over a deep frame, the shape of a double-precision input at high
    precision) shrink the whole computation.  The automatic
    thresholds can be recalibrated with ``tune/tune-exp-reduced``.

    All the arithmetic around the two series kernels and the exact
    binary splitting is ``fball`` ball arithmetic
    (``_fball_exp_reduced``): the kernels are imported with their
    documented errors, the sinh reconstruction is a ball squaring,
    sum and square root, each burst slice's exact splitting output
    `T / (Q B^{Q_e})` becomes the factor `Q B^{Q_e} + T` by one ball
    addition, the products `\mathrm{NUM} \cdot f` and
    `\mathrm{DEN} \cdot Q` are ball products at the frame precision,
    the remainder of the one-step variant is the function itself at
    the doubled rate, and the finish is one ball division; the bound
    comes out rigorous (27 ulps for the cascade, 1 for the full
    burst) and is checked against the budget on export.  That is 200
    lines in place of the earlier driver's hand-written windows,
    exponents and informal error accounting (windowed middle
    products over limb-aligned frames), which measured the same
    within the noise at every size and rate (256 to 65536 limbs, `r`
    = 16 to 1024) and was dropped: the binary splitting kernel is
    85% of the time, the ball products 15%, the division 2%, and the
    driver's own bookkeeping is invisible.

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
    (``_fixed_sin_cos_sum_bs_powtab``), and, in ``fball`` arithmetic
    as :func:`fixed_exp_reduced`, assembles the factors by ball
    additions and one short product (past a per-slice term threshold
    the cosine track is dropped from the tree and recovered from the
    slice's own denominator by one square root,
    `\cos_k Q = \sqrt{(Q B^{Q_e})^2 - (\sin_k Q B^{Q_e})^2}`, a
    ball squaring, difference and square root), accumulates the
    complex numerator by one ``fball_mul_complex`` per slice and the
    real denominator by one product, and finishes with two ball
    divisions against the single accumulated denominator -- where
    the classical per-slice scheme spends a full-precision square
    root per slice.  Everything is limb-granular and dead low limbs
    of a slice fold into its frame, so sparse arguments shrink the
    whole computation (including the internal power tables, whose
    entries strip their dead low limbs into per-entry exponents).
    The earlier driver in hand-written mpn arithmetic measured the
    same phase for phase (once the complex product took the
    high-half path; with the exact product that phase was 20%
    behind, 5% overall) and was dropped.  At 65536 limbs and `r =
    32` the splitting kernels take 70% of the time, the complex
    accumulation 13%, the square roots of the sine-only slices 12%,
    the two divisions 2%.  The automatic thresholds can be
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
    the floor of the value scaled by `2^{\mathrm{FLINT\_BITS} \cdot n}`, or one ulp below it,
    never above.  These build the entries of the cached reduction
    tables (which call them with one guard limb) by binary splitting
    in ``fball`` ball arithmetic, which tracks all truncation errors
    and the series tail; the result is read off the ball as a lower
    bound.  Blocks of terms at the leaves are summed exactly by
    nonallocating mpn code (shifts and fused single-limb
    multiply-accumulates), and when the whole series fits in one such
    block the exact partial sum is simply floored by one division.
    The logarithm sums the alternating series
    `2^{-i} \sum (-1)^k 2^{-ik}/(k+1)` for `i \ge 30` (and for
    `i \ge 20` at moderate precision), otherwise
    `2\operatorname{atanh}(1/(2^{i+1}+1))`, carrying one `q^2` factor
    per term in the materialized denominator so that ranges compose
    without any `q`-power at the merges.

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

Diophantine argument reduction
-------------------------------------------------------------------------------

The following functions compute the elementary functions on the unit
interval by *diophantine* (multi-prime) argument reduction, the
fixed-point counterpart of ``arb_exp_arf_log_reduction`` and
``arb_sin_cos_arf_atan_reduction``: the argument is reduced by an
integer combination of the logarithms of the first primes, or of
`\pi/2` and the doubled arguments `2 \arg(\pi_j)` of the first nonreal
Gaussian primes `\pi_j`, found by a descent through a table of integer
relations between these values, and the correction on the value side
is a rational number (a product of prime powers) or a Gaussian
rational.  Compared to the bitwise reductions the precomputation is
much cheaper -- one Machin-type binary splitting per prime instead of
a table of `r` series -- at the price of a slower evaluation, so they
pay off for moderate numbers of evaluations at a given precision.
``profile/p-diophantine`` prints the trade-off across precisions.

.. function:: void fixed_exp_diophantine(nn_ptr y, nn_srcptr x, slong n)
              void _fixed_exp_diophantine_tune(nn_ptr y, nn_srcptr x, slong n, slong num_primes, double max_weight)

    Sets `(y, n + 1)` to `\exp(x)` for `(x, n)` in `[0, 1)`, within
    ``FIXED_EXP_DIOPHANTINE_MAX_ERR`` ulps.  Writes
    `x = c_0 \log 2 + \sum_j c_j \log p_j + t` with `t` tiny and
    evaluates `\exp(x) = 2^{c_0} (p / q) \exp(t)` with `p, q` products
    of prime powers (truncated to the working precision once they
    exceed it), the reduced `\exp(t)` coming from
    :func:`fixed_exp_reduced`.  The tunable variant takes the number
    of primes (2 to 64) and the weight budget of the descent
    (the sum of `|c_j| \log_2 p_j`; larger budgets give deeper
    reductions and larger products); the default uses 13 primes and
    a budget equal to the precision, as arb does.  The descent is
    bit-identical to arb's.

.. function:: void fixed_sin_cos_diophantine(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n)
              void _fixed_sin_cos_diophantine_tune(nn_ptr ysin, nn_ptr ycos, nn_srcptr x, slong n, slong num_primes, double max_weight)
              void fixed_tan_diophantine(nn_ptr res, nn_srcptr x, slong n)
              void _fixed_tan_diophantine_tune(nn_ptr res, nn_srcptr x, slong n, slong num_primes, double max_weight)

    Set `(ysin, n + 1)`, `(ycos, n + 1)` (either may be ``NULL``) to
    the sine and cosine, respectively `(res, n + 1)` to the tangent,
    of `(x, n)` in `[0, 1)`, within ``FIXED_SIN_COS_DIOPHANTINE_MAX_ERR``
    ulps.  Writes `x = c_0 \pi/2 + \sum_j c_j \, 2 \arg(\pi_j) + t` and
    evaluates `e^{ix} = i^{c_0} e^{it} A^2 / N` for the Gaussian
    integer `A = \prod \pi_j^{c_j}` (conjugates for negative `c_j`),
    whose norm `N = |A|^2` is a rational integer: the normalization is
    one reciprocal of an integer and two multiplications, and the
    tangent is the ratio of the two parts of `e^{it} A^2` with no
    normalization at all.  Gaussian products beyond the working
    precision are carried truncated through the high complex
    products of ``mpn_extras``.  The reduced sine and cosine come from
    :func:`fixed_sin_cos_reduced`, whose cost dominates; since it only
    matches the exponential's efficiency at reductions deeper than
    about `2^{-300}`, the default uses 32 Gaussian primes and a weight
    budget of four times the precision.

.. type:: fixed_rel_struct

.. function:: const fixed_rel_struct * fixed_rel_table(int gaussian, slong num)
              int fixed_rel_table_is_cached(int gaussian, slong num)

    Returns the relation table for the first *num* primes
    (*gaussian* = 0: `\alpha_j = \log p_j`) or Gaussian primes
    (*gaussian* = 1: `\alpha_0 = \pi/2`, `\alpha_j = 2 \arg \pi_j`):
    the rows of ``d`` are integer relations
    `\sum_j d_{ij} \alpha_j = \epsilon_i` with `|\epsilon_i|`
    decreasing (an extra row starting with ``FIXED_REL_TERMINATOR``
    ends the table), and the structure also carries the reciprocals
    `1/\epsilon_i`, the weights, and the primes.  Tables for
    *num* = 2, 4, 6, 8, 10, 12, 13, 16, 20, 24, 32, 40, 48 are
    precomputed; any other size up to ``FIXED_REL_MAX`` is generated
    on first use (``_arb_log_precompute_reductions``; seconds around
    48 primes) and cached per thread, which the second function
    predicts: it returns nonzero when a call would find the table
    ready.  Tables are freed by :func:`flint_cleanup`.

.. type:: fixed_machin_struct

.. function:: const fixed_machin_struct * fixed_machin_table(int gaussian, slong num)
              slong fixed_machin_table_max(int gaussian)
              void fixed_machin_get_x(fmpz_t q, const fixed_machin_struct * tab, slong i)
              void fixed_machin_get_c(fmpz_t c, const fixed_machin_struct * tab, slong i, slong j)
              void fixed_machin_get_c_row(fmpz * row, const fixed_machin_struct * tab, slong i)

    Machin-type sets for the logarithms of the first primes,
    `\log p_i = (1/\mathrm{den}) \sum_j c_{ij} \operatorname{atanh}(1/x_j)`,
    and for the arguments of the first nonreal Gaussian primes,
    `\arg \pi_i = (1/\mathrm{den}) \sum_j c_{ij} \operatorname{atan}(1/x_j)`,
    with square coefficient matrices.  The first function returns the
    best set for *num* values: the largest one with at most *num*
    terms, the remaining values being left to one followup series
    each, unless that would need more followups than the source file's
    measured limit, in which case the next larger set is used.  Sets
    exist for
    every size from 4 (3 for the Gaussian primes) to 32, and for 40
    and 48, on 32- and 64-bit systems alike; the second function
    gives the largest.  The
    arguments are stored as 128-bit values (``FIXED_MACHIN_X_LIMBS``
    limbs each); the coefficient matrices are not stored but
    reconstructed on first use from the arguments, the denominator
    and ``cbits`` (factoring each argument over the primes of the set
    and inverting the exponent matrix modulo one or two primes), and
    cached per thread.  The accessors return any argument or
    coefficient, or a whole row, as ``fmpz``
    (``fixed_machin_c_row_raw`` gives a row as the cached limbs and
    sign bytes).  These are the tables behind
    ``arb_log_primes_vec_bsplit`` and
    ``arb_atan_gauss_primes_vec_bsplit``, and -- evaluated in ``fball``
    arithmetic and combined by exact mpn dot products -- behind the
    cached logarithms of primes and angles of Gaussian primes of the
    diophantine reductions, whose exact floors are read off the balls
    (the first 13 values of each, up to 4608 bits, come directly from
    arb's static tables).  There the `\operatorname{atan}` terms are
    evaluated by :func:`fball_atan_frac_bsplit`, and the
    `\operatorname{atanh}` terms of the logarithms by Zuniga's series
    (below) from a precision of about `b^2/2` limbs for arguments of
    `b` bits, by :func:`fball_atan_frac_bsplit` below.  New sets can
    be generated with https://github.com/fredrik-johansson/machin.

Machin terms as Zuniga series
...............................................................................

Each term of a Machin-type formula for logarithms is a logarithm of a
ratio of smooth numbers: `2 \operatorname{atanh}(1/x) = \log(u/v)`
with `u/v = (x+1)/(x-1)` in lowest terms (so `u - v = 1` for odd `x`
and `u - v = 2` for even `x`).  Zuniga's family of Ramanujan-type
series [Zun2025]_ gives, for rationals `u/v > 1`,

.. math::

    \log\frac{u}{v} = \sum_{n=1}^{\infty} \rho^n
        \frac{\alpha n + \beta}{\gamma\, n (2n-1)}
        \frac{(1)_n (1/2)_n}{(1/6)_n (5/6)_n},
    \qquad
    \rho = \frac{(u-v)^6}{108\, u^2 v^2 (u+v)^2},

.. math::

    \alpha = -2 (u+v)(u^2 - 14uv + v^2)(u^2 + 4uv + v^2), \quad
    \beta = (u+v)^3 (u^2 - 8uv + v^2), \quad
    \gamma = 2 (u-v)^5

(for instance `u/v = 3` gives Zuniga's series for `\log 3`, with
`\rho = 1/243` and `(\alpha n + \beta)/\gamma = 88n - 14`, and
`u/v = 2` the series with `\rho = 1/3888` used for `\log 2`).
Substituting `u = x + 1`, `v = x - 1` (the series is homogeneous in
`u, v`) gives the direct translation of a Machin term:

.. math::

    2 \operatorname{atanh}(1/x) = \frac{x}{4} \sum_{n=1}^{\infty}
        \rho^n \frac{2 (3x^2-1)(3x^2-4)\, n - x^2 (3x^2-5)}{n (2n-1)}
        \frac{(1)_n (1/2)_n}{(1/6)_n (5/6)_n},
    \qquad \rho = \frac{4}{27\, x^2 (x^2-1)^2}.

The arguments and the coefficient matrix of the Machin formula are
unchanged; only the series evaluating each term is replaced.  In the
format of :func:`fball_hypgeom_series`, with
`\mathrm{num}/\mathrm{den} = (u-v)^6 / (6 u^2 v^2 (u+v)^2)` in lowest
terms: `P = [-(u+v)^2 (u^2 - 8uv + v^2),\; 2 (u^2 - 14uv + v^2)(u^2 +
4uv + v^2)]`, `Q = \mathrm{den} \cdot [5, -36, 36]`,
`R = \mathrm{num} \cdot [0, -1, 2]`,
`c_P / c_D = -(u+v)\,\mathrm{num} / (2 (u-v)^5)` and `c_Q = 0`.

One term of the series gains `6 \log_2 x + \log_2(27/4)` bits, as
much as three terms of the Taylor series of `\operatorname{atanh}(1/x)`
and 2.75 bits more.  In binary splitting the cost per bit of precision
is governed by the number of bits by which the numbers outgrow the
precision gained per term: about `\log_2 (2n)` against `2 \log_2 x`
for the Taylor series, and about `2 \log_2 n + 1` against
`6 \log_2 x + 2.75` for Zuniga's series -- two thirds as much
relatively, and much less for small `x`.  In exchange the terms are
three times longer, which costs more at low precision and for large
`x`.  Measured on the logarithms of the first *num* primes
(:func:`_fixed_log_primes_vec_fball`, 64-bit, one thread), Zuniga's
series are faster by a factor 1.45 (*num* = 2, 4), 1.23 (8), 1.13 (13),
1.06 (20) and 1.03--1.05 (32--48) at `10^6` bits, 1.36 (2, 4), 1.31
(8), 1.15 (13), 1.05 (20) and 1.06--1.07 (32--48) at `10^7` bits, and
1.49 (4), 1.28 (8) and 1.19 (13) at `10^8` bits.  For 2--4 primes
(arguments of 8--13 bits) they are faster from the smallest
precisions at which the tables are used, for 32--48 primes
(arguments of 62--87 bits) from about 2000--3000 limbs; the switch
at `b^2/2` limbs follows these crossovers.  (For 2 and 3 primes the
set of 4 primes is used; the best dedicated sets, `x = 7, 17` and
`x = 31, 49, 161`, are no faster with Zuniga's series.)

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

.. function:: void fixed_const_euler(nn_ptr y, slong n)
              void fball_const_euler(fball_t res, slong n)

    The same for Euler's constant `\gamma`: the first function sets
    `(y, n)` to exactly `\lfloor \gamma B^n \rfloor` (cached per
    thread as above, freed with ``_fixed_const_euler_clear``), the
    second gives a ball for `\gamma` accurate to about
    ``FLINT_BITS * (n - 1)`` bits (from arb's 3456-bit table below 53
    limbs).  The algorithm is that of :func:`arb_const_euler`, the
    Brent-McMillan formula ([BM1980]_) with the error bound of
    [BJ2013]_,

    .. math::

        \gamma = \frac{S_0}{I_0} - \frac{K_0}{I_0^2} - \log m
        + O(e^{-8m}), \quad
        I_0 = \sum_{k \ge 0} \left(\frac{m^k}{k!}\right)^2, \quad
        S_0 = \sum_{k \ge 0} \left(\frac{m^k}{k!}\right)^2 H_k,

    `K_0 = I_0(2m) K_0(2m)` from its asymptotic series, with
    `8m \ge b \log 2` for `b` bits, in ``fball`` arithmetic.  `S_0`
    and `I_0` come from one binary splitting over `4.97 m` terms in
    dual numbers: with `\varepsilon^2 = 0` and

    .. math::

        F(\varepsilon) = \sum_{k \ge 0} \frac{m^{2k}}
            {\prod_{j=1}^k (j + \varepsilon)^2},
        \quad I_0 = F(0), \quad S_0 = -\tfrac{1}{2} F'(0),

    the splitting of `F` carries a scalar `P`, `D = D_0 + D_1
    \varepsilon = \prod (k + 1 + \varepsilon)` (so that `Q = D^2`) and
    `T = T_0 + T_1 \varepsilon`, merged as `T = P_1 T_2 + D_2^2 T_1`,
    `D = D_1 D_2`, `P = P_1 P_2`: five full-size and five half-size
    products per merge against twelve mostly full-size ones for the
    quantities `P, Q, T, C, D, V` of the arb code (23% faster overall at
    `10^6` and `10^7` bits).  At the end `X = T_0 + D_0^2`,
    `S_0 / I_0 = (2 T_0 D_1 - T_1 D_0) / (2 D_0 X)` and
    `K_0 / I_0^2 = D_0^4 T_K / (X^2 Q_K)`.  `K_0` takes a second
    splitting over `2m` terms at half the precision; the leaves are
    exact ``fmpz`` products, and the power of two in `P = m^{2 \ell}`
    is not stored.  The final two divisions and ten products cost
    about 1.5% of the total, so fusing the divisions (which trades a
    half-precision division for extra full-precision products) is not
    done.  The main sum runs first, with nothing else live, and leaves
    only the quotient and the two half-precision values `Q_0^2`,
    `X^2`; `K_0` and `\log m` follow, each freeing its scratch space,
    and large right children are freed after their merge.  The peak
    heap use is below that of :func:`arb_const_euler` (54 MB against
    60 MB at `10^7` bits, where keeping all scratch space to the end
    took 141 MB).

    The parameter `m` is rounded up so that `\log m` is a combination
    of a few fast logarithms: over `\{2\}` (`m` a power of two;
    :func:`fball_const_log2`), `\{2, 3\}` (Zuniga's series for
    `\log(9/8)` and `\log(256/243)`), `\{2, 3, 5\}` (`81/80`,
    `32805/32768`, `25/24`) or `\{2, 3, 5, 7\}` (`2401/2400`,
    `4375/4374`, `225/224`, `64/63`), where arb uses `\{2, 3, 5\}` and
    three Machin terms.  Profiling showed that the choice of `m`
    matters more than the cost of the logarithms (6--7% of the total
    for each set): the splittings cost about
    `m (61 + \log_2 o)` for `o` the odd part of `m` (the power of two
    in `m` being free), so an `m` a little larger with a small odd
    part (`98304 = 2^{15} \cdot 3` rather than `87480 = 2^3 \cdot 3^7
    \cdot 5` for `m \ge 86644` at `10^6` bits) is up to 7% faster.
    Every `m` in `[m_0, 1.5 m_0]` smooth over `\{2, 3, 5, 7\}` is
    scored by this cost plus that of its logarithms, and the best is
    taken.  Measured over 32 precisions from `10^4` to `3 \cdot 10^6`
    bits against the best of the four sets at each precision, each
    with its own best `m`, the geometric mean overheads are 1.26 for
    `\{2\}` (up to 1.6 just above a power of two), 1.03 for
    `\{2, 3\}`, 1.02 for `\{2, 3, 5\}` and `\{2, 3, 5, 7\}`, and
    1.002 for the automatic choice (within 1%, the noise level, after
    the dual-number merge).  Against :func:`arb_const_euler` (one
    thread) this is 2.2 times faster at `10^4` bits, 1.65 times at
    `10^5` bits, 1.55 times at `10^6` bits and 1.7 times at `10^7`
    bits.

Further constants
-------------------------------------------------------------------------------

.. function:: void fball_const_e(fball_t res, slong n)
              void fball_const_log10(fball_t res, slong n)
              void fball_const_catalan(fball_t res, slong n)
              void fball_const_zeta3(fball_t res, slong n)
              void fball_const_zeta5(fball_t res, slong n)
              void fball_const_gamma_1_3(fball_t res, slong n)
              void fball_const_gamma_1_4(fball_t res, slong n)

    Balls for `e`, `\log 10`, Catalan's constant `G`, `\zeta(3)`,
    `\zeta(5)`, `\Gamma(1/3)` and `\Gamma(1/4)`, accurate to about
    ``FLINT_BITS * (n - 1)`` bits.  These have no caches of their own;
    the arb constants (:func:`arb_const_e`, :func:`arb_const_log10`,
    :func:`arb_const_catalan`, :func:`arb_const_apery`,
    :func:`arb_const_zeta5` and the `\Gamma(p/q)` constants of
    :func:`arb_hypgeom_gamma_fmpq`) are thin wrappers keeping theirs.

    Except for `\log 10`, each is a hypergeometric series through
    :func:`fball_hypgeom_series`, taken from the formula files of
    y-cruncher (https://github.com/Mysticial/y-cruncher-Formulas); of
    the files available for each constant, the one measured fastest in
    this binary splitting at `10^5` and `10^6` bits was chosen.

    * `e = \sum_{k \ge 0} 1/k!`.
    * `\log 10 = \log 2 + \log 5`, from the three-prime Machin-type
      set of :func:`_fixed_log_primes_vec_fball` (three atanh terms
      from Zuniga's series, on up to three threads), where arb used the
      same three terms as plain arctangent series.
    * `G`: Pilehrood's short series (2010), 10% ahead of Zuniga's 2023
      series and 45% ahead of Guillera's 2019 series.
    * `\zeta(3)`: Zuniga's 2023-vi series, the best of the twelve files
      (2.05 bits per term); the Amdeberhan-Zeilberger series of the old
      arb code is 20% slower.
    * `\zeta(5)`: Zhi-Wei Sun's identity (2025),
      `\zeta(5) = (3S + 56 \pi^2 \zeta(3))/540` with `S` a
      `{}_3F_2`-type series of cost 10.9, against 32.3 for Y. Zhao's
      single series (3.6 times slower in this splitting).
    * `\Gamma(1/3) = (810^{1/4} \pi / X)^{1/3}` with `X` Guillera's
      2023 series; Zuniga's 2024 and Brown's 2011 series are within
      10%.  The series alone is 4 times faster than the one the old arb
      code used.
    * `\Gamma(1/4) = (\pi^6 / (322 S^4))^{1/8}` with `S` Ebisu's 2016
      lemniscate series (Zuniga's 2023-x series is within 1%); the
      eighth root is three square roots.

    The independent parts of `\zeta(5)` and of the two `\Gamma` values
    run on separate threads (at `10^6` bits on two threads, 0.84 s
    becomes 0.68 s for `\zeta(5)` and 0.12 s becomes 0.087 s for
    `\Gamma(1/3)`).

    Against the old arb implementations at `10^6` bits (one thread, the
    whole constant from scratch, `\pi` and `\zeta(3)` included where
    they are needed):

    ===============  ========  ========  =======
    constant         arb       fixed     speedup
    ===============  ========  ========  =======
    `e`              0.022 s   0.016 s   1.4
    `\log 10`        0.174 s   0.085 s   2.0
    `G`              0.259 s   0.196 s   1.3
    `\zeta(3)`       0.144 s   0.121 s   1.2
    `\zeta(5)`       2.342 s   0.844 s   2.8
    `\Gamma(1/3)`    0.153 s   0.123 s   1.2
    `\Gamma(1/4)`    0.094 s   0.109 s   0.9
    ===============  ========  ========  =======

    `\Gamma(1/4)` is the exception: arb computed it from
    `\operatorname{agm}(1, \sqrt 2)`, which costs `O(M(n) \log n)`
    against the `O(M(n) \log^2 n)` of any of these series, so it stays
    ahead at high precision (1.3 s against 1.9 s at `10^7` bits when
    `\pi` is already cached, 1.9 s against 1.9 s when it is not).
    With :func:`fball_agm`, :func:`fball_const_gamma_1_4` takes
    `\Gamma(1/4) = \sqrt{(2\pi)^{3/2} / \operatorname{agm}(1, \sqrt 2)}`
    from `2^{16}` limbs, the AGM in parallel with `\pi`: equal to the
    series up to `4 \cdot 10^6` bits and a few percent ahead at `10^7`
    (1.8 s of AGM against 2.2 s of series alongside `\pi` on two
    threads), the gap growing with the precision.

.. function:: void fball_agm(fball_t res, const fball_t x, const fball_t y, slong n)
              void _fball_agm_order(fball_t res, const fball_t x, const fball_t y, slong n, int m)

    Sets *res* to the arithmetic-geometric mean of `x, y \ge 0` to about
    `n` limbs: the iteration `a' = (a+b)/2`, `b' = \sqrt{ab}` in ball
    arithmetic (the square roots are 72% of the time), finished by the
    series of `\operatorname{agm}(a, b) = ((a+b)/2) / {}_2F_1(1/2, 1/2;
    1; z^2)`, `z = (a-b)/(a+b)`, once `z^2 < 2^{-p/m}`.  The
    coefficients of `1/{}_2F_1(1/2, 1/2; 1; x) = 1 - \sum_{j \ge 1}
    c_j x^j` are dyadic, positive (Kaluza's theorem: the coefficients
    of `{}_2F_1` are log-convex) and sum to 1 (the reciprocal vanishes
    at `x = 1`), so the truncation after `x^{m-1}` has a tail below
    `x^m` at every order `m`; the terms are integers over one power of
    two and the powers of `x` come by squaring at the precision each
    contributes at.  The order `m` (up to 16; 10 on 32-bit limbs) is 8
    below 1024 limbs and 12 above by default, 5--9% faster than `m = 2`
    and 3--5% faster than arb's `m = 5` at high precision; the whole
    is 10--20% faster than :func:`arb_agm` from 64 to 65536 limbs for
    both a full-precision pair and a far-apart one.  When `a - b` no
    longer determines its own size, or `a + b` or `ab` is too wide for
    the division and the square root, the result is the enclosure
    `(a+b)/2 \pm |a-b|/2` (valid since `b \le \operatorname{agm}(a,
    b) \le a`).

.. function:: void fball_root_ui(fball_t res, const fball_t x, ulong k, slong n)
              void fball_rroot_ui(fball_t res, const fball_t x, ulong k, slong n)
              void _fball_root_ui_order(fball_t res, const fball_t x, ulong k, slong n, int r, int recip)

    Sets *res* to `x^{1/k}` respectively `x^{-1/k}` for `x > 0` and
    `1 \le k < 2^{40}` (`k = 2` goes to :func:`fball_sqrt` and
    :func:`fball_rsqrt`).  For other `k` (`k = 3` included, which
    :func:`fball_const_gamma_1_3` uses) the roots come
    from a high-order iteration written in ``fball`` arithmetic, which
    handles the exponents and the propagation of the arithmetic errors:
    with `z` any approximation of `x^{-1/k}` and `u = 1 - x z^k`,

    .. math::

        x^{-1/k} = z (1 - u)^{-1/k}, \qquad
        x^{1/k} = S (1 - u)^{-(k-1)/k}, \quad S = x z^{k-1},

    exactly, and the binomial series `(1-u)^{-b} = \sum_j c_j u^j`
    truncated after `u^{r-1}` has a tail below `2 c_r |u|^r` for
    `|u| \le 1/2` (its coefficients decrease, `b < 1`).  So `z`
    accurate to a fraction `1/r` of the precision, obtained by a
    recursive call whose radius is then discarded (only the size of the
    final `u` matters, and that is computed rigorously), gives the root
    to the full precision from one evaluation of `S`, `u` and the
    series `\sum_{j < r} c_j u^j`, taken over a common denominator
    `D` as `\mathrm{base} + \mathrm{base} \cdot W / D` with
    `W = \sum (D c_j) u^j`, `D` the lcm of the denominators of the
    `c_j` in lowest terms (a power of two for the square roots: 8 at
    order 3, 16 at order 4, and the division a shift) when the
    coefficients are words and `k^{r-1} (r-1)!` otherwise (divided
    out in one ``fball_div_ui`` when it fits a word, else one division
    per factor); the tables of coefficients are built once per root
    and cached per thread.  The powers of `u` come by squaring, one
    product each at the precision they contribute at, the integer
    coefficients by ``fball_mul_ui``, the division acts on `W` (a
    fraction `1 - 1/r` of the precision) and one full product applies
    the base, with the tail added to the radius; measured 3--8%
    faster than a chain of terms each carrying its own product,
    ``mul_ui`` and ``div_ui`` from 256 limbs on.  The recursion starts from a double (through
    the logarithm, as the reduced operand reaches `2^k`).
    The order `r` of the steps is 3 for `k < 8` and 4 above (the
    tuning entry point takes it as an argument; measured: within 5% of
    each other for `k = 3`, the fourth order 5--15% faster for `k = 7`
    at all sizes, 25% at `k = 100` and the largest sizes).

    Two dedicated fixed-point engines were written and measured
    against this iteration before being dropped: a general engine for
    the same iteration with its own floating mantissas and
    hand-assembled product windows (much longer; within a factor 1.3
    of it from a few hundred limbs on and slower above 4096 limbs,
    once :func:`fball_mul` truncated its operands to the precision of
    the product), and a cube root in the mold of
    :func:`fixed_sqrt_newton` (a reciprocal-root Newton step
    `T (1 + u/3)` and a Karp-Markstein finish), which was 1.5--2.5
    times faster below 128 limbs but 25--30% slower from 256 limbs
    on (0.50 ms against 0.37 at 1024 limbs, 11.1 ms against 7.9 at
    16384), the third-order step of the iteration beating a
    second-order Newton step once the products dominate.  Below 64
    limbs the per-operation constant of the ball arithmetic leaves
    the iteration 1.5--2 times behind :func:`arb_root_ui` (MPFR's
    ``mpfr_rootn_ui`` there).  Above, against :func:`arb_root_ui`
    (`k = 7`, one thread): 13 µs against 16 at 64 limbs, 0.34 ms
    against 0.94 at 1024, 7.6 ms against 40 at 16384, 34 ms against
    176 at 65536
    limbs; for `k = 100`, 0.67 ms against 7.9 at 1024 limbs and 15 ms
    against 355 at 16384.  At `10^6` bits the fifth root takes 6.9 ms
    against 22 ms and the seventh 8.3 ms against 31 ms.

Newton-Taylor inverses
-------------------------------------------------------------------------------

Beyond the range of the bitwise reductions, the logarithm and the
arctangent are computed from the forward functions by one
Newton-Taylor step (the fixed-point counterpart of ``arb_log_newton``
and ``arb_atan_newton``): from a starting value `t` correct to about
`p_0` bits, the forward function at the full precision `p` gives a
residual `w` with `|w| < 2^{-p_0}` and `f(x) = t \pm g(w)` exactly,
`g` being the Taylor series of `\log(1 + w)` or `\operatorname{atan}(w)`,
so that `p / p_0` terms of `g` suffice; since the `k`-th term has
`k p_0` leading zero bits, the polynomial costs a fraction of a
multiplication.  The starting value comes from the bitwise function
up to ``FIXED_NEWTON_CUTOFF`` limbs and from the step itself above,
so the whole costs one forward evaluation plus about a tenth of one.

.. function:: void fixed_neglog_newton(nn_ptr y, nn_srcptr x, slong n)
              void _fixed_neglog_newton_tune(nn_ptr y, nn_srcptr x, slong n, int forward, slong N)
              void _fball_neglog_newton(fball_t res, nn_srcptr x, slong n, int forward, slong N)

    Sets `(y, n)` to `-\log(x)` for `(x, n)` in `[1/2, 1)`, within
    ``FIXED_NEGLOG_NEWTON_MAX_ERR`` ulps.  With `E = \exp(t)` from
    :func:`fixed_exp_diophantine` at `n + 1` limbs and `w = x E - 1`,
    `-\log x = t - \log(1 + w)`.  The tunable variant takes the forward
    algorithm (0 = diophantine, 1 = bitwise) and the number `N` of
    terms of the series (up to 16; the starting value is then computed
    at about `1/(N + 1)` of the precision), 0 selecting the measured
    default; the ball-valued worker returns the result with its
    rigorous radius.

.. function:: void fixed_atan_newton(nn_ptr y, nn_srcptr x, slong n)
              void _fixed_atan_newton_tune(nn_ptr y, nn_srcptr x, slong n, int forward, slong N)
              void _fball_atan_newton(fball_t res, nn_srcptr x, slong n, int forward, slong N)

    Sets `(y, n)` to `\operatorname{atan}(x)` for `(x, n)` in `[0, 1)`,
    within ``FIXED_ATAN_NEWTON_MAX_ERR`` ulps.  With
    `(s, c) = (\sin t, \cos t)` from :func:`fixed_sin_cos_diophantine`
    at `n + 1` limbs and `w = (x c - s) / (x s + c)`,
    `\operatorname{atan} x = t + \operatorname{atan}(w)`.  The series
    in `w^2` (up to 13 terms, the starting value at about `1/(2N)`
    of the precision) and the parameters are as for the logarithm.

.. function:: void fixed_neglog_agm(nn_ptr y, nn_srcptr x, slong n)
              void _fixed_neglog_agm_tune(nn_ptr y, nn_srcptr x, slong n, slong N)
              void _fball_neglog_agm(fball_t res, nn_srcptr x, slong n, slong N)

    The same `-\log(x)` by the Sasaki-Kanada formula
    `\log(1/q) = \pi / \operatorname{agm}(\theta_2(q)^2, \theta_3(q)^2)`,
    exact for `0 < q < 1`: with `r = x 2^{-e}` and `q = r^4` (so that
    `\theta_2(q) = 2 r S_2` needs no fourth root),
    `-\log x = (\pi/4) / \operatorname{agm}(4 r^2 S_2^2, \theta_3^2) - e \log 2`,
    `\pi/4` and `\log 2` coming from the per-thread caches of
    :func:`fixed_const_pi_div_4` and :func:`fixed_const_log2`.  The
    exponents of `\theta_3 = 1 + 2 \sum q^{n^2}` and
    `S_2 = 1 + \sum q^{n(n+1)}` merged in order are `\lfloor m^2/4
    \rfloor`, consecutive ones differing by `\lfloor m/2 \rfloor`, so
    both series cost one product per term and one per two terms for
    the powers `q^j` (a squaring `(q^{j/2})^2` for even `j`, each to the
    precision of the term `q^{j^2}` that first uses it; `q^2` and `q^4`
    are squarings too), each at the precision the term contributes at
    (the squarings measured 3--5% at the default `N` and 10% at
    `N = 12`--16); the omitted terms are below `2 q^{k}` for the first omitted
    exponent `k`.  The first AGM step needs no square root:
    `\theta_2^2 \theta_3^2 = (2 r S_2 \theta_3)^2`, and the pair
    becomes the theta functions at `q^{1/2}` by the doubling formulas
    (3% saved).  `e = \lceil p / (4 (N+1)^2) \rceil` for `N` terms of
    `\theta_3`; the AGM spends about `\log_2 e` iterations in its slow
    phase, and each doubling of `N` saves two of them for about `3N`
    more products at decreasing precision, which measured flat from
    `N = 1` to 6 (default 4).  Measured against
    :func:`fixed_neglog_newton` with both caches warm (ms): 6.8
    against 3.5 at 1024 limbs, 46 against 26 at 4096, 250 against 190
    at 16384, 1370 against 1060 at 65536, 7.3 s against 6.5 s at
    262144 and 34.6 s against 36.0 s at `10^6` limbs (19 million
    digits): the `O(M(n) \log n)` of the AGM against the
    `O(M(n) \log^2 n)` of the diophantine exponential crosses over
    near `7 \cdot 10^5` limbs.  Over the table-free exponential
    (``forward = 2``, :func:`fixed_exp_notab`, itself 1.3--1.4 times
    slower than the diophantine one from 65536 limbs on, for the
    sparse starting values of the step as for dense arguments), the
    Newton logarithm takes 238 ms at 16384 limbs, 1.59 s at 65536,
    8.5 s at 262144 and 45.2 s at `10^6` limbs, so the AGM logarithm
    overtakes it near 20000 limbs (a million bits) and is 20--25%
    faster above: with no tables at all, the AGM is the method of
    choice from there on.  The ball-valued worker returns the
    result with its rigorous radius.

Everything around the forward function and the starting value is
``fball`` arithmetic: the residual (its cancellation and the forward
function's documented error carried by the balls; the arctangent's
denominator and division at the precision the cancelled numerator has
left), the series with its powers by squaring at the precision each
term contributes at, the Taylor coefficients scaled to integers by the
lcm of their denominators and one ``fball_div_ui``, the tail added to
the radius from the ball's bound on `|w|`, and the sum with `t`; the
export truncates to `n` limbs and checks the radius, which comes out
around `2^{-40}` ulps since the forward errors enter through `w`
only.  An earlier implementation in hand-written mpn arithmetic
(windowed middle products of the top limbs, a rectangular splitting of
the polynomial with ``mpn_addmul_1`` rows and its own error
accounting) measured the same at every size from 100 to 64000 limbs
-- the forward function is 90% of the time -- and was dropped.
``tune/tune-newton.c`` measures the overhead over the forward
functions (8--13% for the logarithm, 10--20% for the arctangent, whose
division is the difference) and the effect of `N`, which is flat to
within a few percent around the defaults.

Multithreading
-------------------------------------------------------------------------------

The constants, the sets of logarithms and Gaussian-prime arguments and
:func:`fball_hypgeom_series` use the threads of :func:`flint_set_num_threads`.
The results are identical for any number of threads.

* A set of series (the Machin terms of :func:`_fixed_log_primes_vec_fball`
  and :func:`_fixed_atan_gauss_vec_fball` together with their
  followup terms, or the two to four logarithms of `\log m` for
  `\gamma`) runs on up to one thread per series.  The threads take the
  series costliest first (about `n / \log_2 x` terms for `1/x`) as they
  become free, so at most one series per thread is in flight, and
  threads left over when there are fewer series than threads go to the
  splitting of the series.  The Machin combination frees each series
  value as it reads it.

* A binary splitting (`\pi`, `\log 2`, both sums for `\gamma`,
  :func:`fball_atan_frac_bsplit` and the generic backend) splits the
  thread budget between the halves of a node, as
  :func:`flint_parallel_binary_splitting` does, with one difference.
  The halves run on two threads only if the node's exact values (as
  estimated from the bits per term) stay within ``FIXED_PAR_CAP = 4``
  times the working precision.  Above that, each half would hold
  full-size numbers, so both halves run one after the other with the
  whole budget, and the merge runs on two threads: its products fall
  into two independent groups, such as `\{T_1 Q_2, Q_1 Q_2\}` and
  `\{P_1 T_2, P_1 P_2\}`.  Each group's multiplications can use its
  half of the threads for the FFT.  A thread started on a subtree
  gets its own leaf buffers and level temporaries, freed when it is
  done.

Forking at every level from the top, as arb does, is a few percent
faster (3% for `\gamma` at `10^7` bits on two threads) but holds one
full-size working set per thread.  With the cap, the peak memory stays
close to the one-thread peak plus FFT's own per-thread tables.  Peak
heap in MB at `10^7` bits (64-bit, two cores; times on 2 threads),
with ``arb`` for the corresponding arb functions:

    ===========  =====  ======  ======  ======  ======  =============
    constant     code   1 thr   2 thr   4 thr   8 thr   time 1 / 2 thr
    ===========  =====  ======  ======  ======  ======  =============
    `\pi`        arb    27.7    39.5    49.4    57.0    0.86 / 0.49 s
    `\pi`        fball  25.1    36.8    49.4    60.2    0.69 / 0.40 s
    `\log 2`     arb    32.1    52.6    70.8    94.5    1.76 / 0.97 s
    `\log 2`     fball  32.1    48.6    62.5    75.7    1.42 / 0.80 s
    `\gamma`     arb    59.8   100.6   142.3   240.3   30.3 / 16.1 s
    `\gamma`     fball  54.4    79.0    93.7   111.0   18.4 / 10.2 s
    13 logs      arb    48.0    61.5    96.6   174.6   5.03 / 2.87 s
    13 logs      fball  48.0    61.5    92.6   163.3   4.13 / 2.20 s
    13 args      arb    47.8    61.0    90.6   174.9   5.61 / 3.11 s
    13 args      fball  57.4    70.7    96.8   149.0   5.25 / 2.76 s
    ===========  =====  ======  ======  ======  ======  =============

Most of the growth with the number of threads is flint's FFT tables,
which every thread that multiplies at full size keeps until
:func:`flint_cleanup` (19 MB for each thread that multiplies at
full size here): 87 of the 111 MB for
`\gamma` on 8 threads.

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
of the mantissa's bottom limb, i.e. `\mathrm{err} \, B^{\mathrm{exp}
- \mathrm{size}}` (or `\mathrm{err} \, B^{\mathrm{exp}}` for a zero
mantissa), where ``err`` is a single limb, `0 \le \mathrm{err} < B`,
and `0` means exact.  The radius thus has no scale of its own: a
bound below one ulp is installed by padding the mantissa with zero
limbs down to the bound's scale (balls like `1 - \epsilon` carry
their zero limbs explicitly), and a bound of `B` ulps or more by
truncating the mantissa by as many limbs as needed (each dropped
run of limbs being worth one more unit at the new bottom limb), so
that mantissa length tracks the number of accurate limbs with at
most one limb of noise retained (deliberately, since interval radii
over correlated errors grow faster than true errors).  Rounding a
bound up to a whole count at the bottom limb costs up to one ulp
per operation; measured over the `k`-th root iteration and the
constants this loses at most a few bits against the earlier design
of a double-precision radius with its own anchor.

Operations take a precision `n` in limbs, the caller supplying 2--4
guard limbs; outputs are truncated to about `n` limbs and less when
an operand is accurate to fewer, and the products truncate their
operands to `n + 2` limbs beforehand (the dropped tail, below one
unit of the lowest kept limb, enters the radius through the other
operand's magnitude), so that a long exact operand costs no more
than the precision asks; ``fball_div_ui`` divides by a machine word
with ``mpn_divrem_1`` over one appended zero limb (two for an exact
dividend), keeping the quotient's own precision: a quotient
re-anchored at the dividend's bottom limb with its truncation error
there loses `\log_2 c` bits relative to its magnitude at every
division, which a chain of divisions by a common denominator `k^m`
compounded to nothing (found by the root iteration's high-order base
case at `k \approx 2^{40}`); the inherited radius and the truncation
are bounded one limb below the dividend's anchor instead.
Products of fewer than 64 kept limbs are formed in full (the windowed
middle product only pays from there on, measured), and short exact
operands whose product is kept whole take a path without windows,
scratch or bound arithmetic.  ``fball_submul_bounded`` computes
`a - b c` GIVEN a bound `|a - b c| < B^E` from the caller (the
residual of an iteration, whose bound comes from the radius of the
previous approximation): only the window of the product that reaches
the result is computed and the sign is read off its top limb; the
result is wrong if the bound is violated, so the bound must be
rigorous.  Whether this saves anything depends on
:func:`flint_mpn_mulmid` charging for the window rather than for the
full product, which at present it hardly does, so the operation is
within a few percent of a plain product and subtraction.

``fball_mul_complex`` forms `(a_r + i a_i)(b_r + i b_i)` by one
complex product of :func:`flint_mpn_mul_complex` over the four parts
(its high half only above the middle-product cutoff, by
:func:`flint_mpn_mulhigh_n_complex`: four forward and two inverse
transforms in place of the twelve of four products).  Its precision
is that of a complex number, relative to the larger part of the
result with both parts kept down to the same limb: an accumulation of
`\prod \exp(i x_k)` needs its small imaginary part to the absolute
precision of the frame, and asking for the imaginary part's own
relative precision instead was measured to lengthen the windows by
the ratio of the parts (2000 limbs of 16000 in the trigonometric bit
burst).  Each operand's parts are likewise windowed to a common
bottom keep + 2 limbs below its larger part.

The bounds are composed in limb arithmetic: internally a bound is a
128-bit count at a limb anchor (``umul_ppmm`` of a radius by a
magnitude, ``add_ssaaaa`` of two counts), a sum aligning the lower anchor
to the higher by dividing its count by `B` rounding up (one unit when
it sits two or more limbs lower), and the count is reduced to a
single limb, rounding up and moving the anchor, when it is installed
on a ball.  The divisions and roots keep bounding through doubles,
converted to a count on installation; ``fball_div_ui`` bounds the
radius' quotient through a double when that is exact enough.
Magnitudes are bounded through the leading 64 bits of a mantissa
(across a limb boundary), within `2^{-63}`: bounding through the top
limb alone, `|x| < (\mathrm{top} + 1) B^{\mathrm{exp} - 1}`,
overstates by a factor 2 for a top limb of 1 -- every value in
`[1, 2)`, the working set of any iteration near 1 -- and the
derivative bounds of the square roots taken from the blanket
`|x| \ge B^{\mathrm{exp} - 1}` by up to a limb and a half; with the
radius held as a whole count of ulps such pessimism turns directly
into truncated limbs (the bit-burst exponential lost a bit per slice
to it, 195 ulps against 27 now).  These
helpers are force-inlined (a three-word bound passed by value through
memory was measured to cost as much as the mantissa arithmetic of a
one-limb operation), the normalization and installation having inline
fast paths for the common cases.  Sums and differences are formed in
the destination buffer without a temporary (the operands are swapped
so that the result never aliases the second one, and in place when it
aliases the first); a difference is taken from the larger operand,
decided up front from the exponents and top limbs, so that no
negation pass follows, and is subtracted in reverse in place when the
result aliases the smaller; exact operands of at most two limbs whose
exponents are within two limbs are summed in a fixed window of
registers with no bound arithmetic at all.  An in-place sum
``x := x + y`` with `y` inside `x`'s window costs the ``mpn_add`` of
`y`'s limbs plus 55--75 cycles whatever the size of `x` (measured
with `x` of `10^5` limbs and `y` of 1 to 1000, at its bottom, middle
or top, exact or not): the carry or borrow stops when it runs out,
and the normalization and the installation of the radius are
constant time.  A `y` reaching below an inexact `x`'s bottom limb is
cut there (the radius of `x` is at least one ulp of that limb, and
the dropped tail costs one unit), so `x` never moves: about 220
cycles for a 10-limb `y` whether it reaches below or not (it was a
pass over `x` down and back up, 112000 cycles at `10^5` limbs, when
operations kept a pad of six limbs under the noise floor, a relic of
the double-radius design that bought no accuracy with limb-count
radii: products and quotients of inexact operands now keep size + 2
limbs, and the root iteration, the AGM, the logarithms and all the
constants measured the same accuracy without it).  Only a `y` below
an exact `x` (the contiguous mantissa must move) or an `x` longer
than the precision asks (a truncation) cost a pass over `x`.
``fball_addmul_ui`` and ``fball_submul_ui`` compute `x \pm y c` for a
word `c` in place by ``mpn_addmul_1`` resp. ``mpn_submul_1`` over `y`'s
limbs when `y c` lands inside `x`'s window (and, for a difference,
cannot flip its sign; `y` may reach below an inexact `x`, its tail
entering by the high word of its top dropped limb times `c`, within
2 units): 44--58 cycles plus the ``mpn_addmul_1`` for
exact operands, 100--200 for inexact ones, independent of the size of
`x` and about half of ``fball_mul_ui`` followed by ``fball_add``; other
shapes form `y c` and add it.  They are the primitives of a
rectangular splitting or of an interpolation formula summed into one
accumulator.  Measured, a one-limb sum costs
about 65 cycles in all (40 for short exact operands), a one-limb
product 60 (26 when exact and kept whole), a product or quotient by a
word 30--60, against 110--150 for the sums and 70 for the product with
the double radius on the same machine; a difference of 1024-limb
operands costs the same as a sum (it was 1.5 times more), and one
taken in place of its second operand 40% less than out of place (it
was 30% more).  Division and the square roots
require operand balls of relative radius below `2^{-30}` (checked,
not assumed: wider divisors are a usage error whose mantissa would
be pure noise).  ``fball_rsqrt`` and ``fball_sqrt`` align the exponent,
call :func:`fixed_rsqrt_newton` resp. :func:`fixed_sqrt_newton` (or
``flint_mpn_sqrtrem`` below its cutoff) and bound the result through
the derivative; the third-order reciprocal root inside both makes
them as fast as or faster than an earlier arrangement that routed
high precision through the root iteration of :func:`fball_rroot_ui`
and a Karp-Markstein finish in ball arithmetic (measured: square root
3--7% faster, reciprocal root equal to 4% faster at full-precision
inputs, both equal to 10% faster for `x = 2`), which was dropped.
A square always goes through ``flint_mpn_sqr`` in
full, also in the middle-product path and in
``fball_submul_bounded``, at two thirds of a product against the
0.9 of any window of a middle product.  ``fball_get_fixed`` writes a ball known to lie in
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
11--17% slower at every size and has been dropped (re-measured with
a balanced tree and a table of the powers `C^{L 2^i}`: 6--24% slower;
the powers of `C` are about half of `Q`, and splitting them off costs
a full-size product in the `T` merge for the smaller `Q` product it
saves).

``fball_atan_frac_bsplit`` computes `\operatorname{atan}(p/q)` or
`\operatorname{atanh}(p/q)` for mpn integers `0 \le p < q` (with `p/q`
bounded away from 1) by binary splitting of
`x \sum_k s^k (p^2/q^2)^k / (2k+1)`: over a range the sum is
`N / (D Q^{\mathrm{len}-1})` with `D = \prod (2k+1)`, `Q = q^2`.
Leaves of up to 32 terms are exact mpn Horner loops (one fused
``mpn_addmul_1`` / ``mpn_submul_1`` per term when `p^2 (2k+1)` fits a
limb, one ``mpn_mul_1`` for `(2k+1) q^2` when that fits); the tree is
balanced (`2^J` leaves of `L` terms), so the pure powers `Q^{L 2^i}`
(and `P^{L 2^i}`) needed by the merges
`N = N_1 D_2 Q^{b-m} + s^{m-a} P^{m-a} N_2 D_1` can be tabulated by
squarings.  This saves a full-size product per node over carrying the
denominator `D Q^{b-a}` through the tree, but the products by `D` are
only short when `q` is large, and the table costs a final `Q^{N}`:
the table is used when `Q` has at least four times the bits of the
`2k+1` and `n \ge 64` (8--12% faster at `10^6` bits for
`q \ge 2^{40}`), the carried denominator otherwise (3--7% faster for
small `q`).  Against ``arb_atan_frac_bsplit`` it is 5--30% faster.
See ``fixed.h`` for the full interface and per-function contracts.

Hypergeometric series (y-cruncher format)
-------------------------------------------------------------------------------

A generic binary splitting backend for hypergeometric series in the
format of y-cruncher's ``SeriesHypergeometric`` (the formula files in
https://github.com/Mysticial/y-cruncher-Formulas).  It replaces the
``hypgeom`` module, which has been removed: every constant that used it
is now computed here.  The series is

.. math::

    S = \left(\frac{c_Q + c_P \sum_{k \ge 1} \frac{P(k)}{Q(k)}
        \prod_{j=1}^{k-1} \frac{R(j)}{Q(j)}}{c_D}\right)^{\pm 1}

for integer polynomials `P, Q, R` (coefficient `i` of degree `i`,
`Q(k) \ne 0` for `k \ge 1`) and nonzero integers `c_P, c_Q, c_D`
(``CoefficientP``, ``CoefficientQ``, ``CoefficientD`` and ``Power`` in
y-cruncher's files; the prefactors y-cruncher multiplies in, such as an
inverse square root, are left to the caller).  Unlike y-cruncher, the
coefficients may have any size.  Unlike ``hypgeom`` (terms
`A(k)/B(k) \prod_{j=1}^{k} P(j)/Q(j)` with a separate term
denominator `B`), the term and the ratios share the denominator `Q` and
the product runs to `k - 1`, so that most formulas need only small
polynomials and y-cruncher's formula files translate directly.  The
series must converge geometrically: `\deg R \le \deg Q`, with
`|\mathrm{lc}(R)| < |\mathrm{lc}(Q)|` when the degrees agree.

.. type:: fixed_hypgeom_int_struct

    A signed integer `(-1)^{\mathrm{neg}} (d, n)` given as a limb array
    ``d`` of length ``n`` (zero for ``n = 0``) and a sign ``neg``.

.. type:: fixed_hypgeom_series_struct

    The series: ``power`` (1 or -1), ``coefP``, ``coefQ``, ``coefD``,
    and the polynomials as arrays ``P``, ``Q``, ``R`` of
    ``fixed_hypgeom_int_struct`` with lengths ``Plen``, ``Qlen``,
    ``Rlen``.

.. function:: void fball_hypgeom_series(fball_t res, const fixed_hypgeom_series_struct * s, slong n)

    Sets *res* to a ball containing `S`, with about
    ``FLINT_BITS * (n - 1)`` accurate bits (as elsewhere, `n` includes
    the caller's guard limbs).  Throws if the series is not
    geometrically convergent or its tail cannot be bounded.

.. function:: void fball_hypgeom_series_int64(fball_t res, int power, int64_t coefP, int64_t coefQ, int64_t coefD, const int64_t * P, slong Plen, const int64_t * Q, slong Qlen, const int64_t * R, slong Rlen, slong n)

    Convenience wrapper taking the coefficients as ``int64_t``, as they
    appear in the y-cruncher files.  For example, the Chudnovsky
    formula gives `\pi` as `S / \sqrt{10005}` with ``power = -1``,
    `c_P = 1`, `c_Q = 13591409`, `c_D = 4270934400`,
    `P = [-67957045, -2100495856, 23608573992, -57896553024,
    39250089648]`, `Q = [0, 0, 0, -10939058860032000]`,
    `R = [-5, 46, -108, 72]`, and `\operatorname{atan}(p/q)` is
    `S` with `c_P = c_Q = p`, `c_D = q`, `P = R = [p^2, -2p^2]`,
    `Q = [q^2, 2q^2]` (atanh: `P = R = [-p^2, 2p^2]`).

The implementation sums `T = T_1 Q_2 + R_1 T_2`, `Q = Q_1 Q_2`,
`R = R_1 R_2` over a tree whose leaves of 24 to 32 terms (about 24
limbs when that allows) are exact backward recurrences in mpn arithmetic, with the
polynomial values computed by Horner's rule in one word (a whole leaf
at a time, and a leaf specialized to one-word values), two words, or
generic multi-word two's complement; when `R` divides `P` (as for
`\pi`, `\log 2` and most formulas) the leaf uses `T \leftarrow R(k)(T
+ A(k) Q_s)` with `A = P/R`.  A large constant content of `Q` (or
`R`), as `q^2` for `\operatorname{atan}(p/q)`, is split off and
supplied from a table of powers over a balanced tree when it has at
least eight times the bits of the rest (the power-table strategy of
``fball_atan_frac_bsplit``; in the generic code it gains only a few
percent, as the product of the long `T_1` by the short `Q_2'` still
costs about half a full product under FFT multiplication).  The number of terms comes from a rigorous
tail bound: an exact scan of `\prod |R(j)/Q(j)|` with a closed form
beyond the point where every polynomial is dominated by its leading
term (or a Taylor-shifted bound before that), and the bound is added
to the result.  Measured against the special-purpose code (64-bit,
generic including the same final square root for `\pi`): `\pi` and
`\log 2` within about 5% from `n = 100` limbs and within 2--4% (often
faster) from `n = 1000`; `\operatorname{atan}(p/q)` and
`\operatorname{atanh}(p/q)` within 0--8% of ``fball_atan_frac_bsplit``
from `n = 1000` (up to 14% faster for `p > 1`, where the content of
`R` is not split off when small).  Below about `n = 30` the fixed
setup cost (about 1--2 microseconds) dominates.
