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
    just as `|W|` does.


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
