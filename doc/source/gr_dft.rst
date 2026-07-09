.. _gr-dft:

**gr_dft.h** -- fast Fourier transforms over generic rings
===============================================================================

This module implements discrete Fourier transforms of arbitrary length
over generic rings. Given a ring `R`, an integer `n \ge 1`, and a
principal `n`-th root of unity `w \in R` (meaning that
`\sum_{j=0}^{n-1} w^{jk} = 0` for every `k \not\equiv 0 \pmod n`;
over an integral domain, any primitive `n`-th root of unity qualifies),
the forward transform of `(x_0, \ldots, x_{n-1})` is

.. math::

    X_k = \sum_{j=0}^{n-1} x_j w^{jk}

and the inverse transform is

.. math::

    x_j = \frac{1}{n} \sum_{k=0}^{n-1} X_k w^{-jk}.

Over `\mathbb{C}`, taking `w = e^{-2\pi i/n}` (the default chosen by
:func:`gr_dft_default_root`) recovers the classical DFT with the same
convention as the :ref:`acb_dft <acb-dft>` module.

The fast algorithms rely on the principality of `w` (concretely,
identities such as `1 + w^{n/p} + w^{2n/p} + \cdots + w^{(p-1)n/p} = 0`
for the primes `p \mid n`); this condition is only partially checked
by the plan constructors, and using a non-principal root gives
undefined results.

Transforms are performed with respect to a precomputed *plan*
(:type:`gr_dft_pre_t`) which stores the algorithm parameters together
with a table of the roots of unity `w^0, w^1, \ldots, w^{n-1}` and,
optionally, tables for complex Karatsuba multiplication (see below).

Algorithms and flags
-------------------------------------------------------------------------------

The *alg* parameter of the plan constructors selects the algorithm:

* ``GR_DFT_ALG_AUTO`` -- choose automatically based on the length.
* ``GR_DFT_ALG_NAIVE`` -- direct `O(n^2)` evaluation.
* ``GR_DFT_ALG_CT`` -- iterative radix-2 Cooley-Tukey, performing
  `\tfrac{1}{2} n \log_2 n` multiplications by roots of unity.
* ``GR_DFT_ALG_BAILEY`` -- the four-step (Bailey) algorithm, which
  decomposes a transform of length `n = n_1 n_2` with
  `n_1, n_2 \approx \sqrt{n}` into column transforms of length `n_2`,
  pointwise twiddle multiplications, row transforms of length `n_1`,
  and a final transposition. This is cache-friendlier than plain
  Cooley-Tukey for transforms too large to fit in cache. The row and
  column transforms currently use radix-2 Cooley-Tukey sub-plans.
  This is the algorithm supporting multithreading (see below), and
  ``GR_DFT_ALG_AUTO`` selects it for large power-of-two lengths.
* ``GR_DFT_ALG_SPLIT`` -- recursive split-radix, performing about
  `\tfrac{1}{3} n \log_2 n` multiplications by nontrivial roots of
  unity when multiplication by `w^{n/4}` (playing the role of `-i`)
  is free, as in complex mode. In complex mode with Karatsuba products
  enabled, the number of real multiplications is `n \log_2 n - 3n + 4`,
  which meets the minimal
  counts of the Winograd small-FFT (WFTA) modules for
  `n \in \{4, 8, 16\}` (`0`, `4` and `20` real multiplications) and
  is within a few multiplications of the Heideman-Burrus lower bound
  `4n - 2\log_2^2 n - 2\log_2 n - 4` beyond that. This algorithm
  always produces natural ordering and works out of place internally.
* ``GR_DFT_ALG_PFA`` -- the Good-Thomas prime factor algorithm for a
  composite length `n = n_1 n_2` with `\gcd(n_1, n_2) = 1`. With the
  input index map `j = (n_2 j_1 + n_1 j_2) \bmod n` and the CRT output
  index map `k \equiv k_1 \pmod{n_1}`, `k \equiv k_2 \pmod{n_2}`,
  the transform becomes a pure two-dimensional DFT of size
  `n_1 \times n_2` with roots `w^{n_2}` and `w^{n_1}`, requiring *no*
  twiddle factors between the two stages. The sub-transforms use
  automatically chosen sub-plans; with ``GR_DFT_ALG_AUTO``, a general
  length is decomposed recursively into its coprime prime power
  components this way. Requires at least two distinct prime factors.
* ``GR_DFT_ALG_MIXED`` -- recursive mixed-radix Cooley-Tukey
  decimation in time, splitting off one prime factor `p` per
  recursion level. Transforms of prime length `p` are computed by a
  direct kernel using the identity `1 + v + \cdots + v^{p-1} = 0`
  (with `v = w^{n/p}`) to eliminate one root power, which requires
  `(p-1)(p-2)` multiplications by roots of unity instead of the naive
  `(p-1)^2`; for `p = 3`, for instance,
  `X_1 = (x_0 - x_2) + v (x_1 - x_2)` and
  `X_2 = (x_0 - x_1) + v (x_2 - x_1)`. For a prime power `n = p^e`
  with `p` large, a Bluestein sub-plan is used for the length-`p`
  transforms when possible (see below).
* ``GR_DFT_ALG_BLUESTEIN`` -- Bluestein's chirp-z algorithm for odd
  lengths. Writing `t = (n+1)/2` (so that `2t \equiv 1 \pmod n`),
  the identity `jk \equiv t (j^2 + k^2 - (j-k)^2) \pmod n` turns the
  transform into a cyclic convolution of the chirped input
  `w^{t j^2} x_j` with the even, `n`-periodic kernel `w^{-t j^2}`,
  which is embedded in a cyclic convolution of power-of-two length
  `\ge 2n - 1` and evaluated with a power-of-two sub-plan (two fast
  transforms plus a pointwise multiplication by the precomputed
  transformed kernel). All chirp factors are entries of the root
  table. The power-of-two sub-plan requires a root of unity for the
  convolution length, obtained with :func:`gr_dft_default_root`;
  over rings where this fails (e.g. finite fields without roots of
  the required order), Bluestein is unavailable and plans fall back
  to direct prime kernels.

The *flags* parameter is a bitwise combination of:

* ``GR_DFT_SCRAMBLED`` -- the forward transform leaves the output in
  scrambled order and the inverse transform expects its input in the
  same scrambled order. This saves the reordering passes, which is
  useful, for example, when computing convolutions
  (forward transform, pointwise multiplication, inverse transform)
  where the intermediate order is irrelevant. The scrambled order is
  bit-reversed for ``GR_DFT_ALG_NAIVE`` and ``GR_DFT_ALG_CT`` and
  a matrix transposition (a partial bit reversal) for
  ``GR_DFT_ALG_BAILEY``; the precise permutation can be queried with
  :func:`gr_dft_precomp_output_perm`. The flag is only supported for
  power-of-two lengths with these three algorithms, and is cleared in
  the plan (so that natural ordering is used) in all other cases.

Complex mode and Karatsuba multiplication
-------------------------------------------------------------------------------

When the ring is a complex extension `C = R[i]` of a real ring `R`, with
elements of `C` stored as a real and an imaginary part laid out
contiguously as two elements of `R` (as with ``acb`` over ``arb``),
plans can be constructed in *complex mode* with
:func:`gr_dft_precomp_init_karatsuba`, passing the real ring `R` as
*real_ctx* in addition to the ring `C` as *ctx*. It is the user's
responsibility to ensure that the element layout assumption holds; the
constructor verifies at least that the element size of *ctx* is twice
that of *real_ctx*, returning ``GR_UNABLE`` otherwise.

Complex mode always uses the standard root `w = e^{-2\pi i/n}` (there
is no variant taking a user-supplied root), so the special powers of
`w` are known structurally from their exponents, with no value
inspection: the quarter points `\pm 1`, `\pm i` are free rotations,
and the odd multiples of `n/8`, which are of the form `c(1 \pm i)`,
require only 2 real multiplications. Note that no comparison of
computed values could certify these forms over inexact rings such as
``acb``, where the enclosure of, say, `d + c` is a ball around zero;
the cheap multiplication formulas use at most the enclosure of `c`, so
enclosures of the exact products are still obtained over ball rings.
These special rotations are used at any precision. The root table is
computed as `w^j = \operatorname{exp}(-2\pi i j/n)` directly when
possible, giving tighter enclosures over inexact rings than a chain of
multiplications.

In addition, multiplications by the remaining (generic) roots of unity
can be performed with 3 instead of 4 multiplications in `R`, using the
complex Karatsuba formula (also attributed to Gauss) for multiplication
by a precomputed complex constant: for a root `w = c + di` and an
element `x = x_r + x_i i`,

.. math::

    k_1 = c (x_r + x_i), \quad
    \operatorname{Re}(wx) = k_1 - x_i (d + c), \quad
    \operatorname{Im}(wx) = k_1 + x_r (d - c),

where the values `c`, `d - c` and `d + c` are stored in the plan. This
trades one multiplication for three additions per root multiplication,
which is profitable when multiplications in `R` are expensive compared
to additions. The Karatsuba products are therefore only enabled over
exact rings (where the number of multiplications is what matters) and
over inexact rings at sufficiently high working precision, as reported
by :func:`gr_ctx_get_real_prec`; at low precision, a native complex
multiplication (a single multiplication in `C`) is used for the generic
roots instead, which is faster e.g. for ``acb`` at machine-word
precisions.

Drop-in transforms for acb vectors
-------------------------------------------------------------------------------

The module provides drop-in replacements for :func:`acb_dft` and
:func:`acb_dft_inverse` which compute correct enclosures of the DFT of
an ``acb`` vector, using either ball arithmetic or the fixed-point
contexts internally.

.. function:: void gr_dft_acb(acb_ptr w, acb_srcptr v, slong n, slong prec)
              void gr_dft_acb_inverse(acb_ptr w, acb_srcptr v, slong n, slong prec)
              int _gr_dft_acb(acb_ptr w, acb_srcptr v, slong n, int inverse, int which, slong prec)

    Sets `w_k = \sum_j v_j e^{-2\pi i jk/n}` (respectively the
    inverse with the `1/n` scaling), where the output balls contain
    the exact transform of every point of the input balls. *w* may
    alias *v*. The underscore variant selects the internal arithmetic
    (*which*: 0 automatic, 1 ``acb`` ball arithmetic, 2 fixed-point)
    and returns ``GR_UNABLE`` if the requested path cannot handle the
    input (e.g. non-finite midpoints or precision beyond the
    fixed-point range); in automatic mode such inputs are routed to
    the ball path.

    In the fixed-point path, only the input midpoints pass through the
    transform: they are scaled by an exact power of two, determined
    from the magnitude bound of the plan so that every intermediate
    value satisfies the `|t| < 1` requirement, and truncated to fixed
    point with error below 1 ulp per component. The output radii are
    the sum of two independently computed bounds: the fixed-point
    roundoff and root-table error from
    :func:`gr_dft_precomp_nfixed_bound`, and the effect of the input
    radii, which is not propagated through the computation but bounded
    directly by the DFT perturbation inequality
    `|\Delta X_k| \le \sum_j |\Delta x_j|` (the transform matrix
    and its unscaled inverse have unit-modulus entries), evaluated
    once with upward rounding. This is both sharper and cheaper than
    pushing intervals through the transform, and makes the radii of
    the input independent of the working precision of the kernel.

    The number of limbs is chosen exactly, with no over-provisioning:
    the plan is first constructed as a layout only
    (:func:`_gr_dft_precomp_init_layout`), whose magnitude and error
    bounds are independent of the working precision, and the smallest
    limb count for which the absolute roundoff error stays below
    `2^{-\mathrm{prec}}` times the largest input magnitude is
    computed from those bounds before any root of unity is evaluated;
    only then is the plan realized over the fixed-point contexts of
    the chosen precision. At *prec* = 64 this typically selects
    2 limbs (the guard amounting to some tens of bits, growing
    logarithmically with `n`). Output components at the scale of the
    input then carry roughly *prec* relative bits of accuracy. (As
    with any fixed-point scheme, components that are tiny due to
    cancellation have correspondingly fewer relative bits, though
    their enclosures remain correct.)

    At low precision the transform itself is fast enough that the
    surrounding linear work is significant, so it is streamlined and,
    for large lengths, threaded. The conversions between ``arf``
    midpoints and fixed-point components are performed directly on the
    mantissa limbs (:func:`_arf_get_integer_mpn` on input;
    :func:`_arf_set_mpn_fixed` on output, which also performs the
    single rounding to the target precision, so no separate rounding
    pass over the output is needed), the scratch vectors are treated
    as plain limb blocks rather than generic ring vectors, and the
    input is scanned in one pass that simultaneously checks
    finiteness, accumulates the radius bound and locates the largest
    midpoint by pointer. For the inverse, the `1/n` normalization is
    folded into the output exponent for free when `n` is a power of
    two, and otherwise performed as the single rounded division per
    component during output conversion. When `2n` components meet the
    internal threshold (8192) and multiple threads are available, the
    scan and both conversion loops run on the thread pool; the scan
    reduces over a fixed 16-chunk grid combined in canonical order, so
    results (including the upward-rounded radius sums) are bitwise
    independent of the number of threads.

Fixed-point contexts and error bounds
-------------------------------------------------------------------------------

For numerical DFTs, the module provides real and complex fixed-point
contexts as a faster alternative to ``arb`` and ``acb`` (following the
arithmetic used by the fixed-point matrix multiplication in the
``nfloat`` module). An element of the real context with `k`-limb
precision consists of `k + 1` contiguous limbs: a sign limb (0 or 1)
followed by `k` fraction limbs storing an absolute value in `[0, 1)`,
so that representable values `t` satisfy `|t| < 1`. An element of the
complex context consists of a real and an imaginary part laid out
contiguously, matching the layout assumption of complex mode.

Additions and subtractions are exact as long as the result remains in
range. Arithmetic operations whose exact result would reach magnitude
1 return ``GR_UNABLE``, leaving a wrapped value in the output; the
detection is branch-free (the carry is shifted into the status word)
and the flag is absorbed by the status accumulation that transforms
perform anyway, so it costs nothing to intercept at the end of a
large computation. A ``GR_SUCCESS`` status from a whole transform
therefore certifies a posteriori that no intermediate value
overflowed, which makes the guarantees of
``gr_dft_precomp_nfixed_bound`` unconditional: callers need not trust
their own compliance with the magnitude bounds, they can verify it.
(The constants ``one`` and ``neg_one`` clamp to the largest
representable magnitude `1 - \operatorname{ulp}`, a representational
choice covered by the root-table error budgets.)
Real multiplications compute the high product with
``flint_mpn_mulhigh_n``, with error at most 2 ulp, where
`\operatorname{ulp} = 2^{-64 k}`. The contexts install specialized
method tables with fully inlined arithmetic kernels for 1, 2, 3 and 4
limbs (using ``umul_ppmm`` and ``FLINT_MPN_MUL_2X2``, whose products
are exactly truncated with error below 1 ulp) and generic ``mpn``
based methods otherwise.

Complex multiplications use the schoolbook formula with 4 real
multiplications up to a limb-count cutoff (per-component error at most
4 ulp, or 2 ulp for the exact 1- and 2-limb kernels): unlike the
precomputed Karatsuba tables of complex mode (whose entries
`d \pm c` can reach `\sqrt{2}`), all its intermediate values are
bounded by the complex moduli of the operands, so plans over the
fixed-point contexts always disable those tables while still using the
free and cheap rotation classes. From the cutoff onward, the complex
multiplication instead uses the 3-multiplication complex Karatsuba
formula internally, with the sums `x_r + x_i` (of magnitude up to 2)
held as carry-extended values with an extra high limb: the middle
product is obtained from a single ``mulhigh`` on the fraction parts
plus conditional exact additions selected by the carries, and the
final combination runs in exact `(k+1)`-limb arithmetic before storing
back to `k` limbs. This costs 3 instead of 4 ``mulhigh`` calls per
complex product at a per-component error of at most 6 ulp. The default
cutoff is 22 limbs, the measured crossover on x86-64 (where the extra
additions and carry logic take that long to amortize against the saved
multiplication); ``profile/p-gr_dft_nfixed.c`` benchmarks this and the
other kernel-level implementation choices for retuning on other
machines.

The values `1`, `-1`, `i`, `-i` are not representable and are rounded
to magnitude `1 - \operatorname{ulp}`; the transforms never multiply
by these values explicitly. Consequently the contexts are not rings in
the strict sense: ``is_one`` and ``is_neg_one`` return ``T_FALSE``,
and integers of absolute value greater than 1 cannot be assigned. In
particular :func:`gr_dft_inverse_precomp`, which scales by `1/n`,
returns ``GR_DOMAIN`` over these contexts; use the unscaled inverse
:func:`_gr_dft_precomp_raw` and incorporate the scaling into the
caller's own normalization. For the same reason the Bluestein
algorithm is not available over fixed-point contexts (its kernel
normalization and convolution intermediates do not fit the
representation); prime and prime-power lengths automatically fall
back to the direct kernels of the mixed-radix algorithm.

The root table of a fixed-point plan is built from a primitive root
computed via ``arb`` at elevated internal precision, truncated to
fixed point with error at most 2 ulp; the powers up to an eighth of
the circle (or a quarter, or a half, depending on the divisibility of
`n`) are computed as
`w^j = w^{\lfloor j/2 \rfloor} w^{\lceil j/2 \rceil}` (with a
2-multiplication complex squaring for even `j`), and the remaining
entries are exact sign-flipped and part-swapped copies given by the
eight-fold symmetry of the unit circle. The per-entry errors satisfy
a linear recurrence with a closed-form bound of `O(j)` ulp, stored
with the plan; the root tables of sub-plans are exact strided copies
of the parent table and cost nothing to build, and canonical PFA
plans, whose algorithm performs no multiplications by roots at all,
build no table of their own (their ``roots`` field is ``NULL``; the
sub-plans construct their own small tables directly). (Canonical tables over
``acb`` are constructed analogously, from a running product at
elevated precision with the same symmetry fills, rather than by
per-entry exponentials.)

It is the caller's responsibility to scale the input so that every
input, intermediate and output value of the transform satisfies
`|t| < 1`. To this end, :func:`gr_dft_precomp_nfixed_bound` computes,
at plan-construction cost `O(\log n)` (using memoized affine
recursions over the plan structure, evaluated in machine doubles with
conservative fudge factors), a bound *peak* on the modulus of every
value appearing in the computation, together with a bound on the
output errors of the real and imaginary parts measured in ulp. The
computation is safe from overflow if and only if *peak* is less
than 1; since the bounds are linear in the input magnitude, a suitable
input scale is obtained as `c / \operatorname{peak}(1)` for some
margin `c < 1`. The peak includes a factor `\sqrt{2}` at every root
multiplication, accounting for the intermediate sums inside the
2-multiplication rotation classes at odd multiples of `n/8`.

.. function:: int gr_dft_ctx_init_nfixed(gr_ctx_t ctx, slong nlimbs)
              int gr_dft_ctx_init_nfixed_complex(gr_ctx_t ctx, slong nlimbs)

    Initializes *ctx* as the real respectively complex fixed-point
    context with *nlimbs*-limb precision, returning ``GR_UNABLE`` if
    *nlimbs* is not between 1 and 64. The complex context is intended
    to be passed as *ctx*, and the real context as *real_ctx*, to
    :func:`gr_dft_precomp_init_karatsuba`.

.. function:: int gr_dft_nfixed_set_arf(gr_ptr res, const arf_t x, gr_ctx_t ctx)
              void gr_dft_nfixed_get_arf(arf_t res, gr_srcptr x, gr_ctx_t ctx)

    Conversions between elements of a real fixed-point context and
    ``arf`` values. Setting truncates toward zero (error less than
    1 ulp) and saturates values with `|x| \ge 1` to magnitude
    `1 - \operatorname{ulp}`, returning ``GR_DOMAIN`` only for
    non-finite input; getting is exact. The real and imaginary parts
    of a complex element can be converted by treating them as two
    consecutive real elements.

.. function:: void gr_dft_precomp_nfixed_bound(double * peak, double * err_ulps, double in_mag, double in_err_ulps, const gr_dft_pre_t P)

    Given that the inputs of a transform have complex modulus at most
    *in_mag* and componentwise errors of at most *in_err_ulps* ulp,
    computes a bound *peak* on the modulus of all input, intermediate
    and output values, and a bound *err_ulps* on the errors of the
    real and imaginary parts of the output in ulp, for a forward or
    unscaled inverse transform with the plan *P* executed in
    fixed-point arithmetic. The transform is free of overflow iff
    *peak* is less than 1. The model assumes exact additions and
    ``mulhigh`` multiplications with at most 2 ulp error, uses the
    per-component complex multiplication error of the plan's actual
    context (2, 4 or 6 ulp; see above), and uses the root table error
    bound stored in the plan; for plans over other rings the bounds
    describe a hypothetical fixed-point execution with the worst-case
    constants.

.. function:: int _gr_dft_nfixed_roots(gr_ptr roots, ulong n, double * err_ulps, gr_ctx_t ctx)

    Fills *roots* with fixed-point approximations of
    `w^j, 0 \le j < n`, `w = e^{-2\pi i/n}`, over the complex
    fixed-point context *ctx*, writing a bound (in ulp) for the
    errors of the table entries to *err_ulps*. Used automatically by
    the plan constructors over fixed-point contexts.

.. function:: int gr_dft_acb_precomp_init(gr_dft_acb_pre_t Q, slong n, slong prec)
              void gr_dft_acb_precomp_clear(gr_dft_acb_pre_t Q)
              void gr_dft_acb_precomp(acb_ptr w, acb_srcptr v, const gr_dft_acb_pre_t Q, slong prec)
              void gr_dft_acb_inverse_precomp(acb_ptr w, acb_srcptr v, const gr_dft_acb_pre_t Q, slong prec)
              int _gr_dft_acb_precomp(acb_ptr w, acb_srcptr v, int inverse, const gr_dft_acb_pre_t Q, slong prec)

    Precomputed variants, amortizing the internal precomputation over
    repeated transforms of the same length: the object stores the
    realized plan together with the input-independent constants of the
    fixed-point scaling and error analysis (the power-of-two input
    scale itself is chosen per transform from the input magnitudes,
    so a single object serves inputs of arbitrary scale). The
    accuracy is governed by the precision given at initialization;
    the *prec* argument of the transforms only controls the rounding
    of the results. Inputs the fixed-point path cannot handle fall
    back to a one-shot ball transform. Since the plan construction is
    amortized, the automatic backend selection routes large prime
    factors to the ball path at a lower threshold than the one-shot
    interface.

Two-phase plan construction
-------------------------------------------------------------------------------

Internally, plans are built in two phases: a *layout* phase resolving
the algorithm and the complete decomposition (including the layouts of
all sub-plans) without touching any ring elements, and a *realize*
phase binding the ring and computing the root tables, the complex mode
tables and the convolution kernels. Because
:func:`gr_dft_precomp_nfixed_bound` only depends on the layout (the
root-table error of an unrealized layout is a worst-case estimate for
the doubling construction, replaced by the actual bound at
realization), cost and error bounds can be evaluated before choosing
an internal working precision, as done by the ``acb`` drop-in
transforms.

.. function:: int _gr_dft_precomp_init_layout(gr_dft_pre_t P, ulong n, int alg, int flags, int complex_mode)

    Initializes the layout of a plan of length *n* without computing
    any tables. *complex_mode* indicates whether the plan will be
    realized in complex mode (it affects the automatic algorithm
    selection). The layout can be queried with
    :func:`gr_dft_precomp_nfixed_bound` and must be freed with
    :func:`gr_dft_precomp_clear` (or realized).

.. function:: int _gr_dft_precomp_realize(gr_dft_pre_t P, gr_ctx_t real_ctx, gr_ctx_t ctx)

    Realizes a layout over the ring *ctx* with canonical roots,
    passing *real_ctx* (which may be ``NULL``) to enable complex mode,
    after which the plan is ready for transforms. Returns
    ``GR_UNABLE`` (clearing the plan) if the ring cannot provide the
    canonical root table.

Multithreading
-------------------------------------------------------------------------------

Plain Cooley-Tukey and complex-mode split-radix plans (without the
complex Karatsuba tables, which only arise at high precision where
transforms are compute-bound) store their twiddle factors packed per
stage instead of as a serial table of all powers, in the same total
memory: every table access during a transform is then a sequential
walk, where strided reads into a length-`n` table miss the cache once
per entry as soon as the stride exceeds a cache line (measured at
about 28% of a fixed-point transform of length `2^{22}` before the
change). Inverse transforms reuse the forward tables through the
identities `w^{-j r} = -\,\mathrm{stage}_s[h_m - j]` (a reversed
walk, with the negation folded into the butterfly) and, for
split-radix, reversed pair walks with free quarter rotations.

For power-of-two lengths the automatic algorithm selection uses
split-radix in complex mode (fewest multiplications) and plain Cooley-Tukey
otherwise (and whenever scrambled ordering is requested, which
split-radix does not support). The Bailey four-step algorithm remains
available explicitly; in measurements up to `n = 2^{20}` over ball,
fixed-point and word-size modular rings its cache blocking did not
outperform the plain transforms, serially or threaded.

Transforms use worker threads from the global FLINT thread pool when
the ring reports itself threadsafe: threads are requested at transform
time (respecting ``flint_set_num_threads``), or borrowed handles can
be attached with :func:`gr_dft_precomp_set_threads`. Each algorithm
distributes its own structure across as many threads as it can use:

* Cooley-Tukey distributes the top `\log_2(\mathrm{chunks})`
  butterfly passes -- the only levels whose butterflies genuinely span
  the array -- across the workers (one synchronization round per
  pass), after which the array decomposes into independent contiguous
  blocks and each worker runs all the remaining passes over its own
  blocks without further synchronization, in a depth-first blocked
  order (``GR_DFT_CT_BLOCK_BYTES``) so that cache-sized sub-blocks
  complete all their levels while resident; the serial transform uses
  the same order. (Sweeping every level across all workers instead is limited
  by the shared memory bandwidth, which does not grow with the thread
  count; this matters most over cheap rings such as word-size
  modular arithmetic.)

* Split-radix forks its independent sub-transforms (the half-length
  one against the two quarter-length ones, recursively partitioning
  the worker handles) and distributes the `n/4` independent
  iterations of each combine pass.

* Bailey's four-step algorithm distributes the work items of its
  three phases (`n_1` column transforms, `n_2` twiddled rows); the
  transpose is serial.

The only tuning parameter is the granularity: work is divided into at
most `n / \mathrm{serial\_block}` chunks per phase or pass, and
threaded recursion stops below that size
(``GR_DFT_SERIAL_BLOCK_DEFAULT`` = 8 elements by default, deliberately
small; users of cheap rings should raise it with
:func:`gr_dft_precomp_set_serial_block`). No algorithm switching is
performed based on thread availability, so serial and threaded
transforms of the same plan run the same arithmetic. The profile
program ``p-gr_dft_threads`` measures the parallel speedup of each
algorithm by thread count, and ``p-gr_dft_nfixed`` microbenchmarks
the fixed-point arithmetic primitives (carry-chain styles, signed
subtraction strategies, saturation costs, and the classical/Karatsuba
complex multiplication crossover) for machine tuning of the
implementation choices in ``nfixed.c``.

.. function:: void gr_dft_precomp_fprint(FILE * out, const gr_dft_pre_t P)
              void gr_dft_precomp_print(const gr_dft_pre_t P)

    Prints a diagnostic description of the plan to *out* (respectively
    to standard output): the algorithm and decomposition of each node
    (with one indented block per sub-plan), the state of the root and
    kernel tables, the fixed-point table error bounds, and the
    threading configuration. Works on layouts as well as realized
    plans.

Types and macros
-------------------------------------------------------------------------------

.. type:: gr_dft_pre_struct

.. type:: gr_dft_pre_t

    A precomputed plan for transforms of a fixed length over a fixed
    ring, storing the algorithm parameters, the table of roots of unity,
    optional complex Karatsuba multiplication tables, and sub-plans and auxiliary
    tables for the composite algorithms (four-step, prime factor,
    mixed-radix and Bluestein).

Plan construction
-------------------------------------------------------------------------------

.. function:: int gr_dft_default_root(gr_ptr w, ulong n, gr_ctx_t ctx)

    Sets *w* to a default principal *n*-th root of unity: `1` for
    `n = 1`, `-1` for `n = 2`, and `e^{-2 \pi i / n}` otherwise
    (computed via ``exp_pi_i`` or, failing that, via ``pi``, ``i``
    and ``exp``). Returns ``GR_UNABLE`` if the ring does not implement
    the required operations, in which case a root must be supplied
    manually to :func:`gr_dft_precomp_init_root`; this is typically
    the case over finite fields.

.. function:: int gr_dft_precomp_init_root(gr_dft_pre_t P, gr_srcptr w, ulong n, int alg, int flags, gr_ctx_t ctx)
              int gr_dft_precomp_init(gr_dft_pre_t P, ulong n, int alg, int flags, gr_ctx_t ctx)
              int gr_dft_precomp_init_karatsuba(gr_dft_pre_t P, ulong n, int alg, int flags, gr_ctx_t real_ctx, gr_ctx_t ctx)

    Initializes *P* with a plan for transforms of length `n \ge 1`
    over the ring *ctx*, using the given algorithm and flags.
    ``GR_DOMAIN`` is returned if the algorithm does not support the
    length: ``GR_DFT_ALG_CT``, ``GR_DFT_ALG_BAILEY`` and
    ``GR_DFT_ALG_SPLIT`` require a power of two, ``GR_DFT_ALG_PFA``
    requires at least two distinct prime factors, and
    ``GR_DFT_ALG_BLUESTEIN`` requires an odd `n \ge 3` (and a ring
    supporting :func:`gr_dft_default_root`, since a root of unity of
    power-of-two order is needed for the underlying convolution;
    ``GR_UNABLE`` is returned otherwise). With ``GR_DFT_ALG_AUTO``,
    any length is supported: a power of two selects among the
    power-of-two algorithms, a prime power selects
    ``GR_DFT_ALG_MIXED``, and other composite lengths select
    ``GR_DFT_ALG_PFA`` over the coprime prime power components. Over
    rings where Bluestein is unavailable, transforms of lengths with
    large prime factors `p` degrade gracefully to `O(p^2)` prime
    kernels.

    The *root* variant takes a principal *n*-th root of unity *w*,
    treated as generic; the other variants use
    :func:`gr_dft_default_root`. The *karatsuba* variant additionally
    takes the real ring *real_ctx* and enables complex mode (with
    Karatsuba multiplication where profitable) as described above. The constructors partially check that *w* is
    a principal root (verifying `w^{n/2} = -1` for even `n` and
    `w^{n/p} \ne 1` for odd primes `p \mid n`, when these conditions
    are decidable in the ring), returning ``GR_DOMAIN`` on a definite
    failure.

    A sanity check `w^{n/2} = -1` is performed; if this can be decided
    and is false, ``GR_DOMAIN`` is returned. If the return status is
    not ``GR_SUCCESS``, the plan is left in a cleared state and must
    not be used (calling :func:`gr_dft_precomp_clear` is harmless).

    The plan retains pointers to *ctx* and *real_ctx*, which must
    stay valid for the lifetime of the plan.

.. function:: void gr_dft_precomp_clear(gr_dft_pre_t P)

    Frees the data stored in the plan *P*.

.. function:: void gr_dft_precomp_set_threads(gr_dft_pre_t P, thread_pool_handle * threads, slong num_threads)

    Attaches the given worker threads to the plan and, recursively, to
    its sub-plans. The handles are borrowed for the lifetime of the
    plan: they are used by every subsequent transform, are *not* given
    back to the thread pool by :func:`gr_dft_precomp_clear`, and remain
    the responsibility of the caller (typically obtained with
    :func:`flint_request_threads` and released with
    :func:`flint_give_back_threads` after the plan has been cleared).
    Passing *NULL*, 0 detaches the workers, restoring the default
    behavior of requesting threads from the global thread pool during
    each transform. Attached workers are only used when the ring is
    certified thread-safe by :func:`gr_ctx_is_threadsafe` (see above).

.. function:: void gr_dft_precomp_set_serial_block(gr_dft_pre_t P, slong serial_block)

    Sets the threading granularity of the plan and, recursively, of its
    sub-plans (see the discussion above). Passing zero restores the
    default ``GR_DFT_SERIAL_BLOCK_DEFAULT``.

.. function:: void gr_dft_precomp_output_perm(ulong * perm, const gr_dft_pre_t P)

    Writes to *perm* the output permutation of the plan, defined so that
    entry *i* of the forward output contains the DFT coefficient
    `X_{\operatorname{perm}(i)}`. Similarly, the inverse transform
    expects the coefficient `X_{\operatorname{perm}(i)}` in entry *i*
    of its input. This is the identity permutation unless
    ``GR_DFT_SCRAMBLED`` is set and the algorithm supports scrambled
    ordering (other algorithms do not).

Transforms
-------------------------------------------------------------------------------

.. function:: int gr_dft_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_pre_t P, gr_ctx_t ctx)
              int gr_dft_inverse_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_pre_t P, gr_ctx_t ctx)

    Sets *res* to the forward respectively inverse DFT of *vec*, both of
    length `n`, using the plan *P*. Aliasing of *res* and *vec* is
    allowed. The inverse transform includes the scaling by `1/n` and
    returns ``GR_DOMAIN`` if `n` is not invertible in the ring.

.. function:: int gr_dft(gr_ptr res, gr_srcptr vec, ulong n, gr_ctx_t ctx)
              int gr_dft_inverse(gr_ptr res, gr_srcptr vec, ulong n, gr_ctx_t ctx)

    Convenience functions performing a single forward or inverse
    transform of length *n* with a temporary plan constructed with
    default parameters.

Internal functions
-------------------------------------------------------------------------------

.. function:: int _gr_dft_mul_root(gr_ptr res, gr_srcptr x, ulong e, int inverse, gr_ptr rtmp, const gr_dft_pre_t P)

    Sets *res* to `w^e x` (or `w^{-e} x` if *inverse* is set) where
    `0 \le e < n`, using the plain ring multiplication or the complex
    Karatsuba tables depending on the plan. *res* must not alias *x*. In
    complex Karatsuba
    mode, *rtmp* must point to an initialized temporary element of the
    real ring.

.. function:: void _gr_dft_bit_reverse(gr_ptr x, slong stride, int depth, gr_ctx_t ctx)

    Permutes the strided vector *x* of length `2^{depth}` in place by
    the bit reversal permutation.

.. function:: int _gr_dft_ct(gr_ptr x, slong stride, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx)

    In-place radix-2 Cooley-Tukey transform of the strided vector *x*.
    The inverse transform omits the scaling by `1/n`.

.. function:: int _gr_dft_bailey(gr_ptr x, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx)

    In-place four-step transform of the contiguous vector *x*.
    The inverse transform omits the scaling by `1/n`.

.. function:: int _gr_dft_split(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)

    Split-radix transform of *vec*, written to *res* in natural order.
    *res* must not alias *vec*. The inverse transform omits the
    scaling by `1/n`.

.. function:: int _gr_dft_mul_const(gr_ptr res, gr_srcptr x, gr_srcptr w, gr_srcptr wtab3, gr_ptr rtmp, const gr_dft_pre_t P)

    Sets *res* to `w x` for a fixed plan constant *w*, using the
    precomputed complex Karatsuba multiplication table entry *wtab3*
    (three elements `c, d - c, d + c` of the real ring) if it is not
    *NULL*, and a plain ring multiplication otherwise. *res* must not
    alias *x*.

.. function:: int _gr_dft_precomp_raw(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)

    Applies the plan without the `1/n` scaling of the inverse
    transform. Aliasing of *res* and *vec* is allowed.

.. function:: int _gr_dft_naive(gr_ptr res, gr_srcptr vec, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx)

    Naive `O(n^2)` transform. *res* must not alias *vec*.
    The inverse transform omits the scaling by `1/n`.

.. function:: int _gr_dft_mixed(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)

    Mixed-radix transform (``GR_DFT_ALG_MIXED``). *res* must not alias
    *vec*.

.. function:: int _gr_dft_pfa(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)

    Good-Thomas prime factor transform (``GR_DFT_ALG_PFA``). *res*
    must not alias *vec*.

.. function:: int _gr_dft_bluestein(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)

    Bluestein chirp-z transform (``GR_DFT_ALG_BLUESTEIN``). *res* must
    not alias *vec*. The inverse transform is computed as the forward
    transform of the cyclically reversed input.
