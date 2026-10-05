.. _dfloat:

**dfloat.h** -- floating-point expansions and balls on doubles
===============================================================================

This module provides real numbers represented as unevaluated sums of
one to four IEEE doubles (double, double-double, triple-double and
quad-double), rigorous ball arithmetic on top of them, and complex
numbers and complex balls with such parts.
The arithmetic is built entirely from double operations (with a
correctly rounded ``fma``), and is intended to be much faster than
:type:`arb_t` for precisions up to about 200 bits while remaining
rigorous.

Eight real formats are provided:

* ``d1``, ``d2``, ``d3``, ``d4``: a value `x = d_0 + d_1 + \ldots + d_{N-1}`
  of `N` doubles, roughly nonoverlapping (`|d_{i+1}|` is at most about
  `\operatorname{ulp}(d_i)`), carrying about `53 N` bits of precision.
  These are approximate types in the style of ``double``: results are
  rounded, overflow gives infinity, and infinities and NaN propagate.
  ``d1`` is a plain double, included for consistency.
* ``d1b``, ``d2b``, ``d3b``, ``d4b``: balls with a ``dN`` midpoint and a
  double radius, representing every real number `x` with
  `|x - \operatorname{mid}| \le \operatorname{rad}`.

and eight complex ones, ``d1c`` to ``d4c`` (complex numbers with
``dN`` real and imaginary parts) and ``d1cb`` to ``d4cb`` (rectangular
complex balls with ``dNb`` parts), described in the section *Complex
numbers* below.

The ball types implement the real numbers rigorously in the sense of
interval arithmetic: for every operation, the output ball contains the
exact image of the input balls, including the cases where the double
arithmetic underflows or overflows. The conventions for the whole
real line and for the domains of the functions are set out under
*Nonfinite values and domains* below.

Design
-------------------------------------------------------------------------------

**Error bounds.** The radii are loose (constant factors are not
optimized) but they are not a priori: every midpoint operation is
carried out with error-free transformations (TwoSum, TwoProd) that
produce a list of exact terms, the list is renormalized to `N`
components, and the magnitudes of the dropped terms are added to the
radius. The only a-priori bounds are for partial products of order
`\ge N` in multiplication (bounded by their own rounded magnitude),
for the Taylor remainder of the exponential function, and for the
rounding of the radius arithmetic itself, which is done in nearest
rounding with a single slack factor `1 + 2^{-45}` applied at the end.
As a consequence, exact inputs give exact outputs whenever the
floating-point operations on the midpoints were exact: `[1 \pm 0] +
[1 \pm 0] = [2 \pm 0]`, `[3 \pm 0] \cdot [3 \pm 0] = [9 \pm 0]`,
`(2^{60} + 1) \cdot 3` is exact in ``d2b``, and so on.
The elementary functions (``exp``, ``expm1``, ``sin``, ``cos``,
``log``, ``log1p``, ``atan``) are the exception to the tracking: they
are pure floating-point kernels with worst-case error bounds proved
once and for all (below), and the input radius is propagated
separately; the remaining functions (``tan``, ``sinh``, ``cosh``,
``tanh``, ``asin``, ``acos``, ``atan2``, ``pow``) are compositions of
these with the arithmetic.

**Error modes of the products.** Tracking the rounding errors of a
product exactly costs about as much as the product itself for
`N = 4` (the residuals of the last order, which the plain product
discards, must be computed). The ball products (``mul``, ``sqr`` and
the residual products inside ``div`` and ``sqrt``) therefore track
exactly only when both inputs are exact (radius zero), which is what
preserves exact results, and otherwise *bound* the rounding errors of
the last order by `u` times the accumulated values (``kernels.inc``,
``track = 2``), which is cheaper and only loses a constant factor on
a radius that is nonzero anyway. The flag ``DFLOAT_FAST`` selects the
bound always; it makes no difference in speed for inexact inputs and
loses exactness (`(2^{60} + 1) \cdot 3` then gets a radius of about
`2^{-46}`), so it exists mainly for experiments.

**Renormalization.** The kernels renormalize with a single bottom-up
VecSum, which is all the error analysis needs (*weak*
renormalization): the components are nonoverlapping in the sense that
each is the rounding error of a sum whose rounded value is the
previous one, but after a cancellation a component can be smaller than
the one following it, so that a result needing all `53 N` bits may
come out with a tiny nonzero radius, and the representation is not
canonical (``equal`` compares values, not components). The functions
``dN_canonicalise`` bring an expansion into *canonical* form exactly:
`|x_{i+1}| \le \operatorname{ulp}(x_i)/2`, zeros only at the end,
which is the form produced by rounding the value to nearest term by
term. A context created with the flag ``DFLOAT_STRONG`` canonicalises
the result of every arithmetic operation (*strong* renormalization).
For `N = 2` the weak form already is canonical; for `N = 3, 4` about
4-7% of the results of random operations are not, and the
canonicalisation costs a top-down TwoSum chain plus a check, with a
branchy repair pass in about `10^{-5}` to `10^{-4}` of the cases.
The penalty is measured below.

**Underflow and overflow.** The error-free product `a b = p + e`
requires that `a b` does not underflow. The multiplication kernels
check that all nonzero components are at least `2^{-480}` in absolute
value (so that no partial product falls below `2^{-960}`) and
otherwise take a slow path that scales the operands to unit exponent
(``scaled.c``); components that would become subnormal are flushed
into the radius. Division and square root use the scaled path when
the heads are outside `[2^{-150}, 2^{150}]`. Radius products that
would underflow are bumped up; a midpoint that overflows becomes the
whole real line. Subnormal numbers are handled correctly but may be
slow, as with any double arithmetic.

**Nonfinite values and domains.** The ball rings implement the field
of real numbers, with these conventions, which every function
follows (the ``dNb`` functions as well as the generic rings):

* A ball with a nonfinite radius or midpoint component (`\pm \infty`
  or NaN) represents the whole real line `W`: ``[0 +/- inf]``,
  ``[5 +/- inf]``, ``[inf +/- 0]``, ``[nan +/- 1]`` and ``[1, nan +/- 0]``
  are all `W`. Every function accepts any of these forms and treats
  it as `W`. The functions return `W` only in the normal form
  `[0 \pm \infty]` (all midpoint components `+0`, radius `+\infty`);
  any other result has finite midpoint components and a finite
  nonnegative radius. (``dNb_set`` copies its input as is; so do the
  data movement methods of the generic rings: ``swap``, ``vec_set``,
  ``vec_gather``, ``vec_scatter``.) `W` results from overflow, from
  nonfinite input and from functions that are unbounded on the input.

* Operations defined on all of `\mathbb{R}` (``add``, ``sub``, ``mul``,
  ``sqr``, ``neg``, ``abs``, ``mul_2exp_si``, ``exp``, ``expm1``,
  ``sin``, ``cos``, ``sin_cos``, ``atan``, ``sinh``, ``cosh``, ``tanh``,
  ``atan2``, ``sgn``, the dot products,
  integer powers with a nonnegative exponent, and division by a ball
  not containing zero) always succeed. On `W` they give `W`, or a ball
  containing the range where the function is bounded: `[0 \pm 1]` for
  ``sin``, ``cos``, ``tanh`` and ``sgn``, `[0 \pm \pi/2]` for ``atan``
  and `[0 \pm \pi]` for ``atan2`` (which, like :func:`arb_atan2`, has
  `\operatorname{atan2}(0, 0) = 0`). An exact zero times anything,
  including `W`, is exactly zero, and `x^0 = 1` exactly for every `x`.

* Operations defined on a subset `D` of `\mathbb{R}` (or
  `\mathbb{R}^2`) return ``GR_SUCCESS`` only if the input balls are
  contained in `D`, ``GR_DOMAIN`` only if they are disjoint from `D`,
  and ``GR_UNABLE`` otherwise (the output is then unspecified). Since
  `W` contains points on both sides of every boundary, an operation
  with a nonfinite argument returns ``GR_UNABLE`` unless the other
  argument alone decides. The domains: ``div`` and ``inv``, a nonzero
  divisor; ``sqrt``, `x \ge 0`; ``rsqrt`` and ``log``, `x > 0`;
  ``log1p``, `x > -1`; ``asin`` and ``acos``, `|x| \le 1`; ``tan``,
  `x \notin \pi/2 + \pi \mathbb{Z}` (never ``GR_DOMAIN``, since the
  poles are irrational); ``pow``, `x > 0`, or `x = 0` and `y \ge 0`,
  or `x < 0` and `y \in \mathbb{Z}`. For example:

      ==================================== ========================================
      operation                            status
      ==================================== ========================================
      `[0 \pm \infty] / 3`                 ``GR_SUCCESS``, `[0 \pm \infty]`
      `1 / 0`                              ``GR_DOMAIN``
      `[0 \pm \infty] / 0`                 ``GR_DOMAIN``
      `[0 \pm \infty] / [3 \pm 2]`         ``GR_SUCCESS``, `[0 \pm \infty]`
      `1 / [1 \pm 2]`                      ``GR_UNABLE``
      `1 / [0 \pm \infty]`                 ``GR_UNABLE``
      `\exp([0 \pm \infty])`               ``GR_SUCCESS``, `[0 \pm \infty]`
      `\sin([\infty \pm 0])`               ``GR_SUCCESS``, `[0 \pm 1]`
      `\sqrt{-1}`                          ``GR_DOMAIN``
      `\sqrt{[0 \pm 1]}`                   ``GR_UNABLE``
      `\log 0`, `\log [-2 \pm 1]`          ``GR_DOMAIN``
      `\log [0 \pm \infty]`                ``GR_UNABLE``
      `\operatorname{asin}(1)`             ``GR_SUCCESS``, `\pi/2`
      `\operatorname{asin}([2 \pm 0.5])`   ``GR_DOMAIN``
      `0^{-1/2}`                           ``GR_DOMAIN``
      `(-2)^{1/2}`                         ``GR_DOMAIN``
      `(-2)^{[0 \pm \infty]}`              ``GR_UNABLE``
      `[0 \pm \infty]^3`                   ``GR_SUCCESS``, `[0 \pm \infty]`
      ==================================== ========================================

* The predicates of the generic rings (``is_zero``, ``equal``, ``cmp``,
  ...) answer ``T_TRUE``/``T_FALSE`` (or ``GR_SUCCESS``) only when
  every point of the balls decides the question, and ``T_UNKNOWN``
  (``GR_UNABLE``) otherwise, in particular for `W`; they may answer
  ``T_UNKNOWN`` in borderline cases.

* Conversions: ``gr_set_d`` of `\pm \infty` or NaN returns
  ``GR_DOMAIN`` (they are not real numbers; ``dNb_set_d`` gives `W`);
  ``gr_set_other`` from an :type:`arb_t` that is not finite gives `W`;
  ``gr_get_fmpz`` (and ``get_si``, ``get_ui``) succeeds only for an
  exact integer, returns ``GR_DOMAIN`` for a ball that contains no
  integer (or an integer out of range) and ``GR_UNABLE`` otherwise;
  ``gr_get_d`` and ``dNb_get_d`` return 0 for `W`;
  ``gr_get_interval_mid_rad`` returns ``GR_UNABLE`` for `W` (whose
  radius is not a real number); every form of `W` prints as
  ``[+/- inf]``.

These rules have been checked by fuzzing against :type:`arb_t` (every
method of the generic rings and every ``dNb`` function, for all `N` and
flags, on random and special balls). The complex balls follow the
same rules part by part, on `\mathbb{C}` (see *Complex numbers*).

**Requirements.** Round-to-nearest IEEE double arithmetic with a
correctly rounded ``fma()``, and every operation on doubles rounded
once to double precision. Hardware FMA (x86-64-v3 or later, or any
64-bit ARM) is required for performance; without it ``fma()`` is a
slow library call but the results are still correct. The error-free
transformations break if the compiler contracts a separate multiply
and add into an ``fma``; the module is compiled with
``-ffp-contract=off``, as detected by ``configure`` and CMake, and its
sources also turn contraction off with a pragma for MSVC and clang
(``src/dfloat/fp_contract.h``). They also break with the excess
precision of x87 arithmetic (``FLT_EVAL_METHOD`` 2, the GCC default
for 32-bit x86): there ``configure`` compiles the module with
``-mfpmath=sse`` if the target has SSE2 (``CFLAGS`` with ``-msse2``,
or a ``-march`` that implies it). When either requirement cannot be
met (another compiler that does not take ``-ffp-contract=off``, or
32-bit x86 without SSE2), the module declines to work:
:func:`dfloat_is_supported` returns 0, :func:`gr_ctx_init_dfloat`
returns ``GR_UNABLE`` and the power sums return 0, so that their
callers fall back to other code.

The vector code is written against ``machine_vectors.h`` (``vec4d``
with masks, bitwise and integer-lane operations, fused multiply-adds
and table gathers) and is enabled with AVX2 and FMA and with NEON on
AArch64; elsewhere the vector functions use the scalar code, with
the same results. Defining ``DFLOAT_FORCE_SIMD`` when compiling the
module uses the vector code with whichever backend
``machine_vectors.h`` selects, including its generic (GNU vector
extension or ISO C) backends, and ``DFLOAT_NO_SIMD`` disables it.

.. function:: int dfloat_is_supported(void)

    Returns whether the module works on this target, i.e. whether it
    was compiled with doubles evaluated in double precision
    (``FLT_EVAL_METHOD`` 0). Where it returns 0, the functions
    documented here must not be called directly;
    :func:`gr_ctx_init_dfloat` returns ``GR_UNABLE``.

**Performance.** Nanoseconds per operation on one core of an Intel
Xeon (Skylake-AVX512 class, 2.8 GHz, a virtual machine; the absolute
numbers vary by up to 1.5x between runs on this machine, so only
ratios within a table are meaningful), with :type:`nfloat` (64, 128,
192, 256 bits, the closest fixed-precision software floats in FLINT)
and :type:`arb_t` (at `53 N` bits) for comparison
(``src/dfloat/profile/p-arith.c``). For the arithmetic, the scalar
columns are calls on independent operands through the generic ring;
``d2s``, ``d2bs`` etc. are the strong (canonicalising) contexts and
``d2bf``, ``d2bfs`` the fast ones. The vector columns are per element
of a vector of 1024 (``vec_add``, ``vec_mul``, and ``vec_dot``, a dot
product), which use SIMD four elements at a time for the ``d`` types:

    ========== ===== ===== ===== ===== ===== ======= ======= =======
    ring         add   mul   sqr   div  sqrt vec_add vec_mul vec_dot
    ========== ===== ===== ===== ===== ===== ======= ======= =======
    d1           1.4   1.3   1.4   1.8   1.9    0.27    0.25    0.34
    d1s          1.5   1.4   1.8   1.8   1.9    0.25    0.26    0.31
    d1b          2.8   4.8   3.2   8.8   6.7    0.75     1.6     2.1
    d1bs         2.3   4.1   3.5   8.8   7.2    0.77     1.6     2.0
    d1bf         2.4   4.3   3.5   8.6   6.9    0.77     1.4     1.9
    d1bfs        2.4   4.4   3.3   8.1   7.2    0.76     1.4     1.9
    nfloat64     6.4   4.5   3.5    23    21     3.2     1.6     2.0
    arb 53        24    20  17.4    55    76      21    18.8     8.6
    d2           3.3   3.0   3.5   8.3   7.4     1.1    0.94     1.5
    d2s          4.5   4.7   4.4   9.1   8.8     2.0     1.9     1.3
    d2b          5.5   8.8   6.7    30    26     3.1     4.7     6.2
    d2bs         8.1  11.4   7.9    36    33     4.3     5.2     5.7
    d2bf         6.4   9.3   7.2    32    29     3.1     4.1     5.5
    d2bfs        6.9  10.0   7.6    36    37     4.3     5.0     5.5
    nfloat128    6.7   5.5   3.9    44    88     3.8     2.7     4.8
    arb 106       36    22    20    80   114      34      21    10.7
    d3           7.4   8.8   8.1    25    39     3.8     3.5     4.0
    d3s         10.7  13.1  11.4    29    43     6.7     6.4     4.7
    d3b         10.3  15.7  11.4    64    80     4.8     7.4     8.8
    d3bs        13.8    20  15.5    74    84     8.1    11.2    10.2
    d3bf        11.7    20  14.3    76    84     4.8     7.6     8.3
    d3bfs       14.1  19.1  15.3    64    84     7.6     9.8     8.3
    nfloat192    7.2   7.2   5.0    84   118     4.4     6.7     8.1
    arb 159       41    33    32    80   130      40      32      25
    d4          14.5    20  14.3    59    69     5.2     8.3     8.6
    d4s         19.1    27    23    76    76     9.5    12.9     8.8
    d4b         17.4    28    20   114   118     8.6    11.0    13.3
    d4bs          24    38    29   118   126    13.3    16.0    13.8
    d4bf        17.2    30    20   114   118     8.8    11.2    12.9
    d4bfs         24    35    26   122   122    13.3    16.2    12.6
    nfloat256    7.2   8.3   6.4    99   118     5.2     8.3    10.5
    arb 212       41    37    35    99   153      41      35      26
    ========== ===== ===== ===== ===== ===== ======= ======= =======

Against ``nfloat``, the plain expansions win clearly for `N \le 2`
(the ``d1`` operations are single instructions) and the balls are
comparable for `N \le 2` while carrying a rigorous error bound; the
vector operations and dot products are 1.5-3x faster than ``nfloat``
for `N \le 3`. For `N = 3, 4`, ``nfloat`` addition and multiplication
(limb arithmetic) are faster than the expansion operations, by 2x
for the products at `N = 4`. Division and square root of expansions
are a long division and a Newton iteration with tapered precision
(each step at the number of components it needs; the digit `q_k =
r_0 / y_0` is a single division, whose head `q_k y_0` is cancelled
exactly against the remainder by a TwoProd and Sterbenz's lemma, so
that no accuracy is lost at the cancellation) and are 1.5-4x faster
than ``nfloat``. A Newton-Karp-Markstein division (a reciprocal at
half precision, one product and one correction) is slower at every
`N` (75 vs 49 ns at `N = 4` in the same run): the
long division's steps are cheap products of a single double with
the divisor, and it has no full-precision square. Precomputing
`1/y_0` and multiplying does not help either: the divisions in the
chain are not on the critical path of the products, and the
division unit is idle otherwise. The overhead of the balls over the
plain numbers is 1.2-1.7x for the products with the bounding mode
above, and 2-4x for division and square root (the exact residual
`x - q y` costs a product and a difference).

The penalty of strong renormalization is nil for `N = 1` (nothing to
do), a check for `N = 2` (1-2 ns), and a dependent TwoSum chain for
`N = 3, 4`: about `+5` ns and `+8` ns per scalar operation (+60% on
a plain addition, +30% on a ball addition, +25% on a ball product),
and up to +100% on the cheapest vector operations, whose results are
canonicalised one at a time after the SIMD block.

The elementary functions ``exp``, ``sin``, ``cos``, ``sin_cos``,
``log`` and ``atan`` are pure floating-point kernels with static
error bounds; scalar calls are bound by the latency of the dependent
chain of expansion operations, which is kept short by evaluating the
polynomial as four chains in `w^2` (or `w^4`) in the four SIMD lanes
at once (and, for the trigonometric functions, the table and
combination products in the lanes as well). For `N = 1` the kernels
are ordinary double code (one table level, ``fma`` polynomials) and
about as fast as the ``libm`` functions, whose near-correct rounding
is not attempted (the observed errors are `2^{-52}` to `2^{-51}`;
``libm`` on this machine: ``exp`` 4.6, ``sin`` 6.7, ``log`` 3.3,
``atan`` 4.8 ns against 5.6, 4.8, 3.4 and 5.7). The vector versions
``_dN_vec_exp`` etc. evaluate the same kernels on four arguments at a
time with AVX2, 3-5x faster per element. Scalar calls (through the
generic ring) against the vector functions (``_dN_vec_exp``,
``_dNb_vec_exp`` etc., per element of a vector of 1024), with
:type:`arb_t` (whose vector versions are loops over the scalar
function) for reference:

    ========== ===== ======= ===== ======= ===== ======= ===== ========
    ring         exp vec_exp   sin vec_sin   log vec_log  atan vec_atan
    ========== ===== ======= ===== ======= ===== ======= ===== ========
    d1           8.3     1.8   6.2     2.7   4.5     3.3   6.9      6.0
    d1b         15.3     3.5   7.2     3.4   6.9     6.0   8.8      8.3
    arb 53       176     179   198     202   153     149   179      183
    d2            80    13.6    76      25    59      23    74       22
    d2b           88    14.8    80      26    58      22    76       24
    arb 106      214     210   244     248   210     210   248      236
    d3           153      32   187      56   160      48   183       53
    d3b          156      34   168      58   172      49   176       50
    arb 159      278     305   305     305   252     248   320      305
    d4           305      62   320     111   305      92   366       99
    d4b          305      58   320     111   305      92   382      103
    arb 212      336     351   412     382   305     305   382      382
    ========== ===== ======= ===== ======= ===== ======= ===== ========

(nanoseconds per element). The ball versions cost about the same as
the plain ones for `N \ge 2` (the radius is a static relative bound
plus the propagated input radius) and up to 2x more for `N = 1`, where
the plain kernels are a few instructions. Against ``arb``, the scalar
functions are 20-40x faster for `N = 1`, 3x for `N = 2` and about as
fast for `N = 4`, the vector functions 50-100x, 10x and 4-6x.
``sin_cos`` costs about 10% more than ``sin`` alone; ``expm1`` and
``log1p`` cost the same as ``exp`` and ``log`` (6.2, 64, 118, 259 and
3.5, 43, 130, 259 ns for `N = 1, \ldots, 4`), and the compositions
the sum of their parts (they have no vector versions).

**Argument reduction by division?** For ``log``, the reduction
`z = m r - 1` with a tabulated reciprocal `r` costs a TwoProd and a
product by a double (a few ns) per level and needs no division;
two levels bring `|z|` down to `2^{-14}`. A division-based
reduction `z = (m - c) / (m + c)` (with `\log m = \log c + 2
\operatorname{atanh} z`) would halve the number of coefficients by
using only odd powers, but the even/odd split of the polynomial runs
the two parities in parallel lanes anyway, so nothing would be gained
for the price of a division (46 ns at `N = 4`). For ``atan``, one
division is unavoidable (`z = (x - c) / (1 + x c)`; the numerator is
exact); a second level would shorten the chains by two tapered steps,
which cost about 4 ns at `N = 4` (the schedule is `[4, 3, 2, 1]`: the
inner steps run at one or two components), against a second division
of 46 ns. Tapered chains make the degree of the polynomial almost
irrelevant to the latency; what matters is the number of operations
at full precision, which a table lookup keeps at one or two.

Context objects
-------------------------------------------------------------------------------

.. macro:: DFLOAT_BALL
           DFLOAT_STRONG
           DFLOAT_FAST
           DFLOAT_COMPLEX

    Flags for :func:`gr_ctx_init_dfloat`.

.. function:: int gr_ctx_init_dfloat(gr_ctx_t ctx, int n, int flags)

    Initializes *ctx* to the ring of `n`-term double expansions
    (an approximate floating-point ring) or, with ``DFLOAT_BALL`` set
    in *flags*, the ring of balls with such midpoints (a rigorous
    representation of the real numbers). With ``DFLOAT_STRONG`` set,
    the arithmetic operations (``add``, ``sub``, ``mul``, ``sqr``,
    ``div``, ``inv``, ``sqrt``, ``rsqrt``, ``exp``, ``log``, ``sin``,
    ``cos``, ``sin_cos``, ``atan`` and the vector operations)
    canonicalise their results (see Design). With ``DFLOAT_FAST`` set
    (balls only), the products always bound their rounding errors
    instead of tracking them exactly for exact inputs (see Design).
    With ``DFLOAT_COMPLEX`` set, the ring is that of the complex
    numbers ``dNc`` (plain) or complex balls ``dNcb`` (with
    ``DFLOAT_BALL``), with the other flags as for the real rings (see
    *Complex numbers*). Passing 0 or 1 as *flags* selects the plain or
    the ball real ring with the defaults.
    Returns ``GR_UNABLE`` if `n` is not between
    1 and :macro:`DFLOAT_MAX_N`, or if :func:`dfloat_is_supported`
    returns 0. The element size is `8 (n + \text{ball})` bytes, twice
    that for the complex rings. The real ball ring reports itself as a
    field, ordered, real, inexact and non-canonical, like
    :func:`gr_ctx_init_real_arb`, and the complex ball ring as an
    algebraically closed field like :func:`gr_ctx_init_complex_acb`;
    the precision reported by :func:`gr_ctx_get_real_prec` is `53 n`
    and cannot be changed.

    The rings implement only what has native code: the arithmetic,
    the elementary functions listed below and their vector versions.
    Special functions (``gamma``, ``zeta``, ``erf``, ...) have no
    methods, so that the generic code returns ``GR_UNABLE`` for them:
    evaluating them through :type:`arb_t` would only hide the absence
    of fast implementations. The constants ``pi``, ``euler`` and
    ``catalan`` are table lookups (5-term expansions; the balls get the
    dropped terms and the table residual in the radius).

Basic functions
-------------------------------------------------------------------------------

The functions below exist for each format: ``dN`` (plain) and ``dNb``
(ball) for `N = 1, 2, 3, 4`; in the text, ``dN_exp`` stands for
:func:`d1_exp`, ..., :func:`d4_exp` and ``dNb_exp`` for their ball
versions. Ball functions that can fail return a status flag; the others
return nothing. Aliasing of inputs and outputs is allowed everywhere.

.. function:: void d1_set_d(d1_t res, double x)
              void d2_set_d(d2_t res, double x)
              void d3_set_d(d3_t res, double x)
              void d4_set_d(d4_t res, double x)
              void d1b_set_d(d1b_t res, double x)
              void d2b_set_d(d2b_t res, double x)
              void d3b_set_d(d3b_t res, double x)
              void d4b_set_d(d4b_t res, double x)

    Sets *res* to the double *x* (exactly; a nonfinite *x* gives the
    whole real line in the ball version).

.. function:: void d1_set_dd(d1_t res, double x0, double x1)
              void d2_set_dd(d2_t res, double x0, double x1)
              void d3_set_dd(d3_t res, double x0, double x1)
              void d4_set_dd(d4_t res, double x0, double x1)
              void d1b_set_dd(d1b_t res, double x0, double x1)
              void d2b_set_dd(d2b_t res, double x0, double x1)
              void d3b_set_dd(d3b_t res, double x0, double x1)
              void d4b_set_dd(d4b_t res, double x0, double x1)

    Sets *res* to `x_0 + x_1`, renormalized (the ball version records
    anything that does not fit in the radius).

.. function:: void d1b_set_d_rad(d1b_t res, double x, double rad)
              void d2b_set_d_rad(d2b_t res, double x, double rad)
              void d3b_set_d_rad(d3b_t res, double x, double rad)
              void d4b_set_d_rad(d4b_t res, double x, double rad)

    Sets *res* to the ball `[x \pm \text{rad}]` (the whole line if *x* or
    *rad* is not finite or *rad* is negative).

.. function:: void d1_zero(d1_t res)
              void d2_zero(d2_t res)
              void d3_zero(d3_t res)
              void d4_zero(d4_t res)
              void d1_one(d1_t res)
              void d2_one(d2_t res)
              void d3_one(d3_t res)
              void d4_one(d4_t res)
              void d1b_zero(d1b_t res)
              void d2b_zero(d2b_t res)
              void d3b_zero(d3b_t res)
              void d4b_zero(d4b_t res)
              void d1b_one(d1b_t res)
              void d2b_one(d2b_t res)
              void d3b_one(d3b_t res)
              void d4b_one(d4b_t res)
              void d1b_indeterminate(d1b_t res)
              void d2b_indeterminate(d2b_t res)
              void d3b_indeterminate(d3b_t res)
              void d4b_indeterminate(d4b_t res)

    Constants; ``indeterminate`` is the whole real line in its normal
    form `[0 \pm \infty]`.

.. function:: void d1_set(d1_t res, const d1_t x)
              void d2_set(d2_t res, const d2_t x)
              void d3_set(d3_t res, const d3_t x)
              void d4_set(d4_t res, const d4_t x)
              void d1_neg(d1_t res, const d1_t x)
              void d2_neg(d2_t res, const d2_t x)
              void d3_neg(d3_t res, const d3_t x)
              void d4_neg(d4_t res, const d4_t x)
              void d1_abs(d1_t res, const d1_t x)
              void d2_abs(d2_t res, const d2_t x)
              void d3_abs(d3_t res, const d3_t x)
              void d4_abs(d4_t res, const d4_t x)
              void d1b_set(d1b_t res, const d1b_t x)
              void d2b_set(d2b_t res, const d2b_t x)
              void d3b_set(d3b_t res, const d3b_t x)
              void d4b_set(d4b_t res, const d4b_t x)
              void d1b_neg(d1b_t res, const d1b_t x)
              void d2b_neg(d2b_t res, const d2b_t x)
              void d3b_neg(d3b_t res, const d3b_t x)
              void d4b_neg(d4b_t res, const d4b_t x)
              void d1b_abs(d1b_t res, const d1b_t x)
              void d2b_abs(d2b_t res, const d2b_t x)
              void d3b_abs(d3b_t res, const d3b_t x)
              void d4b_abs(d4b_t res, const d4b_t x)

    Copy, negation and absolute value (``dNb_set`` copies any form of
    the whole line as is; the others return its normal form).

    Assignment, negation and absolute value (exact).

.. function:: void d1_add(d1_t res, const d1_t x, const d1_t y)
              void d2_add(d2_t res, const d2_t x, const d2_t y)
              void d3_add(d3_t res, const d3_t x, const d3_t y)
              void d4_add(d4_t res, const d4_t x, const d4_t y)
              void d1_sub(d1_t res, const d1_t x, const d1_t y)
              void d2_sub(d2_t res, const d2_t x, const d2_t y)
              void d3_sub(d3_t res, const d3_t x, const d3_t y)
              void d4_sub(d4_t res, const d4_t x, const d4_t y)
              void d1_mul(d1_t res, const d1_t x, const d1_t y)
              void d2_mul(d2_t res, const d2_t x, const d2_t y)
              void d3_mul(d3_t res, const d3_t x, const d3_t y)
              void d4_mul(d4_t res, const d4_t x, const d4_t y)
              void d1_sqr(d1_t res, const d1_t x)
              void d2_sqr(d2_t res, const d2_t x)
              void d3_sqr(d3_t res, const d3_t x)
              void d4_sqr(d4_t res, const d4_t x)
              void d1_mul_d(d1_t res, const d1_t x, double y)
              void d2_mul_d(d2_t res, const d2_t x, double y)
              void d3_mul_d(d3_t res, const d3_t x, double y)
              void d4_mul_d(d4_t res, const d4_t x, double y)
              void d1_div(d1_t res, const d1_t x, const d1_t y)
              void d2_div(d2_t res, const d2_t x, const d2_t y)
              void d3_div(d3_t res, const d3_t x, const d3_t y)
              void d4_div(d4_t res, const d4_t x, const d4_t y)
              void d1_inv(d1_t res, const d1_t x)
              void d2_inv(d2_t res, const d2_t x)
              void d3_inv(d3_t res, const d3_t x)
              void d4_inv(d4_t res, const d4_t x)
              void d1_sqrt(d1_t res, const d1_t x)
              void d2_sqrt(d2_t res, const d2_t x)
              void d3_sqrt(d3_t res, const d3_t x)
              void d4_sqrt(d4_t res, const d4_t x)
              void d1_rsqrt(d1_t res, const d1_t x)
              void d2_rsqrt(d2_t res, const d2_t x)
              void d3_rsqrt(d3_t res, const d3_t x)
              void d4_rsqrt(d4_t res, const d4_t x)
              void d1_mul_2exp_si(d1_t res, const d1_t x, slong e)
              void d2_mul_2exp_si(d2_t res, const d2_t x, slong e)
              void d3_mul_2exp_si(d3_t res, const d3_t x, slong e)
              void d4_mul_2exp_si(d4_t res, const d4_t x, slong e)

    Arithmetic on the plain types, with results accurate to a few ulps
    at `53 N` bits (the elementary functions, with static bounds, are
    documented below). The arithmetic operations produce *compact* weak
    expansions: a component that comes out exactly zero (after an
    exact cancellation) is squeezed out, so that the components always
    decrease geometrically; the tapered kernels (``div``, ``sqrt``,
    the elementary functions) rely on the positions of the components
    and lose accuracy on an expansion with an internal zero, such as
    `(2^{-86}, 0, 2^{-140}, 0)`, which the module therefore never
    produces (a user-built one is canonicalised by the elementary
    functions, but not by ``div``, ``sqrt`` and ``rsqrt``).

    ``rsqrt`` is the reciprocal square root `x^{-1/2}`, by the Newton
    iteration `r \leftarrow r + r (1 - x r^2)/2` from the double
    reciprocal root, doubling the number of correct components each
    step (the first step in doubles from the exact products `r_0^2`
    and `x_0 r_0^2`, the later ones with the level kernels); it costs
    about as much as ``sqrt`` (``sqrt`` is the corresponding iteration
    with a tapered division) and is accurate to about `8 u^N`.

.. function:: void d1b_add(d1b_t res, const d1b_t x, const d1b_t y)
              void d2b_add(d2b_t res, const d2b_t x, const d2b_t y)
              void d3b_add(d3b_t res, const d3b_t x, const d3b_t y)
              void d4b_add(d4b_t res, const d4b_t x, const d4b_t y)
              void d1b_sub(d1b_t res, const d1b_t x, const d1b_t y)
              void d2b_sub(d2b_t res, const d2b_t x, const d2b_t y)
              void d3b_sub(d3b_t res, const d3b_t x, const d3b_t y)
              void d4b_sub(d4b_t res, const d4b_t x, const d4b_t y)
              void d1b_mul(d1b_t res, const d1b_t x, const d1b_t y)
              void d2b_mul(d2b_t res, const d2b_t x, const d2b_t y)
              void d3b_mul(d3b_t res, const d3b_t x, const d3b_t y)
              void d4b_mul(d4b_t res, const d4b_t x, const d4b_t y)
              void d1b_sqr(d1b_t res, const d1b_t x)
              void d2b_sqr(d2b_t res, const d2b_t x)
              void d3b_sqr(d3b_t res, const d3b_t x)
              void d4b_sqr(d4b_t res, const d4b_t x)
              void d1b_mul_d(d1b_t res, const d1b_t x, double y)
              void d2b_mul_d(d2b_t res, const d2b_t x, double y)
              void d3b_mul_d(d3b_t res, const d3b_t x, double y)
              void d4b_mul_d(d4b_t res, const d4b_t x, double y)
              void d1b_mul_2exp_si(d1b_t res, const d1b_t x, slong e)
              void d2b_mul_2exp_si(d2b_t res, const d2b_t x, slong e)
              void d3b_mul_2exp_si(d3b_t res, const d3b_t x, slong e)
              void d4b_mul_2exp_si(d4b_t res, const d4b_t x, slong e)

    Rigorous ball arithmetic; see the design notes above for the
    error bounds and the exactness guarantee.

.. function:: int d1b_div(d1b_t res, const d1b_t x, const d1b_t y)
              int d2b_div(d2b_t res, const d2b_t x, const d2b_t y)
              int d3b_div(d3b_t res, const d3b_t x, const d3b_t y)
              int d4b_div(d4b_t res, const d4b_t x, const d4b_t y)
              int d1b_inv(d1b_t res, const d1b_t x)
              int d2b_inv(d2b_t res, const d2b_t x)
              int d3b_inv(d3b_t res, const d3b_t x)
              int d4b_inv(d4b_t res, const d4b_t x)

    Division. Returns ``GR_DOMAIN`` if the divisor is exactly zero
    (whatever the dividend), ``GR_UNABLE`` if it contains zero (in
    particular, if it is the whole line), and ``GR_SUCCESS`` otherwise. The bound used is
    `|x'/y' - q| \le (|x - q y| + r_x + |q| r_y) / (|y| - r_y)` for the
    computed quotient `q`, with `x - q y` enclosed by the tracked
    kernels.

.. function:: int d1b_sqrt(d1b_t res, const d1b_t x)
              int d2b_sqrt(d2b_t res, const d2b_t x)
              int d3b_sqrt(d3b_t res, const d3b_t x)
              int d4b_sqrt(d4b_t res, const d4b_t x)

    Square root. Returns ``GR_DOMAIN`` if *x* is entirely negative,
    ``GR_UNABLE`` if it contains negative numbers, and ``GR_SUCCESS``
    otherwise (with `\sqrt{[0 \pm 0]} = [0 \pm 0]`).

.. function:: int d1b_rsqrt(d1b_t res, const d1b_t x)
              int d2b_rsqrt(d2b_t res, const d2b_t x)
              int d3b_rsqrt(d3b_t res, const d3b_t x)
              int d4b_rsqrt(d4b_t res, const d4b_t x)

    Reciprocal square root. Returns ``GR_DOMAIN`` if *x* is entirely
    nonpositive, ``GR_UNABLE`` if it contains zero, and ``GR_SUCCESS``
    otherwise. The bound used is, for the computed `r` and `e' = 1 -
    x' r^2`, `x^{-1/2} = r (1 - e')^{-1/2}`, so that `|x'^{-1/2} - r|
    \le |r| (|e'|/2)(1 + |e'|)` for `|e'| \le 1/4`, with `|e'| \le
    |1 - x r^2| + r_x r^2` and `1 - x r^2` enclosed by the tracked
    kernels (two products and a subtraction; the ball version is
    therefore relatively more expensive than the ball ``sqrt``, whose
    residual is one product).

.. function:: void d1_canonicalise(d1_t res, const d1_t x)
              void d2_canonicalise(d2_t res, const d2_t x)
              void d3_canonicalise(d3_t res, const d3_t x)
              void d4_canonicalise(d4_t res, const d4_t x)
              void d1b_canonicalise(d1b_t res, const d1b_t x)
              void d2b_canonicalise(d2b_t res, const d2b_t x)
              void d3b_canonicalise(d3b_t res, const d3b_t x)
              void d4b_canonicalise(d4b_t res, const d4b_t x)

    Sets *res* to *x* (exactly; the radius is unchanged) with the
    components in canonical form: `|x_{i+1}| \le \operatorname{ulp}(x_i)/2`
    with zeros only at the end, for any *x* whose components overlap
    by at most a few bits, in particular any result of this module.
    A nonfinite head is copied as is by ``dN_canonicalise``;
    ``dNb_canonicalise`` gives the whole line in its normal form.

.. macro:: DFLOAT_EXP_EPS_1
           DFLOAT_EXP_EPS_2
           DFLOAT_EXP_EPS_3
           DFLOAT_EXP_EPS_4
           DFLOAT_EXP_EPS(n)

    The static relative error bounds of the exponential kernels,
    about `2^{-49.7}`, `2^{-100.8}`, `2^{-151.8}` and `2^{-202.7}`: for
    `-746 \le x \le 709.9`, ``dN_exp`` returns `\exp(x) (1 + \delta)` with
    `|\delta| \le` ``DFLOAT_EXP_EPS_N``, up to an additional absolute
    error of `2^{-1070}` when the result is subnormal. The bounds are
    derived in ``dev/dfloat_exp_bound.py``; the observed worst case
    on random arguments is 4-80 times smaller.

.. function:: void d1_exp(d1_t res, const d1_t x)
              void d2_exp(d2_t res, const d2_t x)
              void d3_exp(d3_t res, const d3_t x)
              void d4_exp(d4_t res, const d4_t x)

    The exponential function, with the guarantee above. The input is
    first canonicalised (exactly). With `k = \operatorname{round}(x_0 / \log 2)`,
    the head is reduced to `h = x_0 - k L_0` exactly (`L_0` is
    `\log 2` to 41 bits, so `k L_0` is exact), then to
    `s_0 = h - i/256 - j/65536` with `|s_0| \le 2^{-17}`, and the
    reduced argument `s = s_0 + (x_1 + \cdots + x_{N-1}) - k (L_1 +
    \cdots + L_N)` is assembled with the expansion kernels in an
    order that keeps every term significant. Then with
    `T = \exp(i/256) \exp(j/65536)` and `w = s^2`,
    `\exp(x) = 2^k (T + (T s) E(w) + (T w) O(w))` where
    `E(w) = \sum_j w^j/(2j+1)!` and `O(w) = \sum_j w^j/(2j+2)!` are
    the even and odd parts of the Taylor polynomial of ``expm1``
    (degree 3, 6, 8, 11), evaluated by Horner's rule with tapered
    precision (fewer components at the inner levels, the schedule
    chosen by ``dev/dfloat_exp_bound.py``). The scalar kernel splits each of
    `E` and `O` once more into its even and odd parts in `w`, so
    that the four chains `E_0(v), E_1(v), O_0(v), O_1(v)` in
    `v = w^2` have half the length and run in the four SIMD lanes at
    once; the products `(T s) E_0 + (T s w) E_1 + (T w) O_0 + (T
    w^2) O_1` are one lane-wise product and the final five-term sum
    a single kernel operation, all of which keeps the dependent chain
    short (the vector kernel, whose lanes are busy with four
    arguments, uses the two chains in `w`). The tables (179 + 257
    entries of 5 doubles, nearest-rounded expansions) fit in L1
    cache. The exactness of the reduction (the three two-sums) is
    proved in the source (Sterbenz's lemma and ulp counting) rather
    than checked; the canonical form of the input is checked at run
    time, and the function falls back to :type:`arb_t` if it fails
    after the canonicalisation (which cannot happen). For `N = 1` the
    kernel is a plain double computation (one table level `i/256`,
    the tail `k L_1` by a TwoProd, the polynomial of degree 6 by
    ``fma``, `2^k` by an exponent addition), with the same style of
    analysis. Nonfinite arguments and arguments beyond the range give
    the double result (``inf``, 0, ``nan``).

.. function:: void d1b_exp(d1b_t res, const d1b_t x)
              void d2b_exp(d2b_t res, const d2b_t x)
              void d3b_exp(d3b_t res, const d3b_t x)
              void d4b_exp(d4b_t res, const d4b_t x)

    The exponential function on balls: the kernel above is applied
    to the midpoint `m`, and the radius is
    `|y| (1 + 2 \varepsilon) ((e^r - 1) + 2 \varepsilon)` with
    `\varepsilon` the static bound and `e^r - 1` bounded by
    `r + r^2` for `r \le 1` and by a power of two otherwise (since
    `|\exp(x') - \exp(m)| \le \exp(m) (e^r - 1)` for `|x' - m| \le r`,
    and `\exp(m) \le |y| / (1 - \varepsilon)`). Where that bound
    overflows (a wide ball) but the upper end `h \ge m + r` of the
    ball does not, the result is `[0 \pm \exp(h)]` (``expm1``:
    `[-1 \pm \exp(h)]`), which also covers the midpoints below the
    double range: `\exp([-10^{20} \pm 10^5]) \subset [0 \pm 2^{-1074}]`.
    Arguments beyond the double range give the whole real line
    (overflow) or a ball of radius at most `2^{-1074}` around zero
    (underflow); a subnormal result gets the radius bumped by
    `2^{-1069}`. `\exp([0 \pm r]) =
    [1 \pm (e^r - 1)]`, exactly 1 for `r = 0`.

.. macro:: DFLOAT_TRIG_EPS_1
           DFLOAT_TRIG_EPS_2
           DFLOAT_TRIG_EPS_3
           DFLOAT_TRIG_EPS_4
           DFLOAT_TRIG_EPS(n)

    The static error bounds of the sine and cosine kernels, about
    `2^{-50.4}`, `2^{-99.1}`, `2^{-148.5}` and `2^{-200.1}`: for
    `|x| \le 2^{20} \pi/2 - 1`, ``dN_sin`` returns `\sin x` with an
    absolute error of at most `\varepsilon (|x| + \min(|x|, 1))` and
    ``dN_cos`` returns `\cos x` with an absolute error of at most
    `\varepsilon (|x| + 1)`, `\varepsilon =` ``DFLOAT_TRIG_EPS_N``.
    The bounds are absolute, scaled by `|x|`: the relative accuracy
    near the zeros of the functions (multiples of `\pi/2`) is not
    preserved, since the argument reduction is done with `\pi/2` to
    about `53 N + 48` bits only and the reduced argument keeps the
    absolute accuracy `u^N |x|` of the input (the ``sin`` bound is
    `2 \varepsilon |x|` for `|x| \le 1`, a relative one where
    `\sin x \approx x`). The bounds are derived in
    ``dev/dfloat_exp_bound.py`` with the model of the exponential;
    the observed worst case on random arguments is about 5-100
    times smaller.

.. function:: void d1_sin(d1_t res, const d1_t x)
              void d2_sin(d2_t res, const d2_t x)
              void d3_sin(d3_t res, const d3_t x)
              void d4_sin(d4_t res, const d4_t x)
              void d1_cos(d1_t res, const d1_t x)
              void d2_cos(d2_t res, const d2_t x)
              void d3_cos(d3_t res, const d3_t x)
              void d4_cos(d4_t res, const d4_t x)
              void d1_sin_cos(d1_t sn, d1_t cs, const d1_t x)
              void d2_sin_cos(d2_t sn, d2_t cs, const d2_t x)
              void d3_sin_cos(d3_t sn, d3_t cs, const d3_t x)
              void d4_sin_cos(d4_t sn, d4_t cs, const d4_t x)

    The sine and cosine, with the guarantee above (``sin_cos``
    computes both from one reduction, at about 10% more than one of
    them; the results are identical to those of ``sin`` and ``cos``).
    The input is canonicalised (exactly). With `k = \operatorname{round}(2
    x_0 / \pi)`, the head is reduced to `h = x_0 - k P_0` exactly
    (`P_0` is `\pi/2` to 33 bits, `|k| \le 2^{20}`), the reduced
    argument `t = h + (x_1 + \cdots) - k (P_1 + \cdots + P_N)` is
    assembled with the kernels and canonicalised (so that the tail of
    `t` is bounded by `u |t|` and not by `u |x|`, which matters for
    large `x`), and split as `t = i/64 + j/2^{13} + s` with
    `|s| \le 2^{-14}`. Then with `A = i/64 + j/2^{13}`, `S_A = \sin A`
    and `C_A = \cos A` from the products of two table levels (51 + 65
    entries, sign symmetry), and `w = s^2`,

    .. math::

        \sin t = S_A + (S_A w) Q(w) + (C_A s) P(w), \quad
        \cos t = C_A + (C_A w) Q(w) - (S_A s) P(w)

    where `\sin s = s P(s^2)` and `\cos s = 1 + s^2 Q(s^2)` are the
    Taylor polynomials (degrees 7/8, 7/8, 11/10, 13/12 for `N = 1,
    \ldots, 4`), and `\sin x`, `\cos x` are `\pm \sin t`, `\pm \cos t`
    by the quadrant `k \bmod 4`. As for the exponential, the scalar
    kernel evaluates `P_0, P_1, Q_0, Q_1` in `v = w^2` in the four
    SIMD lanes; the four table products, the four multipliers `S_A
    w`, `C_A s`, `C_A w`, `S_A s` and their products with `w` are
    also computed as one vector product each, so that each of the two
    outputs costs one more vector product and a five-term sum. For
    `|x| < 2^{-300}` the results are `x` and 1 (which keeps all
    intermediate quantities of the analysis normal); for `N = 1` the
    kernel is plain double code with one table level. Arguments
    beyond `2^{20} \pi/2 - 1` in absolute value (where `k P_0` would
    no longer be exact) are evaluated through :type:`arb_t`, with the
    same guarantee; nonfinite arguments give ``nan``.

.. function:: void d1b_sin(d1b_t res, const d1b_t x)
              void d2b_sin(d2b_t res, const d2b_t x)
              void d3b_sin(d3b_t res, const d3b_t x)
              void d4b_sin(d4b_t res, const d4b_t x)
              void d1b_cos(d1b_t res, const d1b_t x)
              void d2b_cos(d2b_t res, const d2b_t x)
              void d3b_cos(d3b_t res, const d3b_t x)
              void d4b_cos(d4b_t res, const d4b_t x)
              void d1b_sin_cos(d1b_t sn, d1b_t cs, const d1b_t x)
              void d2b_sin_cos(d2b_t sn, d2b_t cs, const d2b_t x)
              void d3b_sin_cos(d3b_t sn, d3b_t cs, const d3b_t x)
              void d4b_sin_cos(d4b_t sn, d4b_t cs, const d4b_t x)

    The sine and cosine on balls: the kernel is applied to the
    midpoint `m` and the radius is `\varepsilon (|m| + a) + r` with
    `a = \min(|m|, 1)` for the sine and 1 for the cosine (the
    functions are 1-Lipschitz), rounded up with the usual slack and
    bumped by `2^{-1069}` when tiny; a radius of 1 or more is replaced
    by the trivial enclosure `[0 \pm 1]`. Infinite radii, nonfinite
    midpoints and midpoints beyond the kernel's range go through
    :type:`arb_t`. `\sin([0 \pm 0]) = [0 \pm 0]` and `\cos([0 \pm 0]) =
    [1 \pm 0]` exactly.

.. macro:: DFLOAT_LOG_EPS_1
           DFLOAT_LOG_EPS_2
           DFLOAT_LOG_EPS_3
           DFLOAT_LOG_EPS_4
           DFLOAT_LOG_EPS(n)
           DFLOAT_ATAN_EPS_1
           DFLOAT_ATAN_EPS_2
           DFLOAT_ATAN_EPS_3
           DFLOAT_ATAN_EPS_4
           DFLOAT_ATAN_EPS(n)

    The static relative error bounds of the logarithm and arctangent
    kernels: `|` ``dN_log(x)`` `- \log x| \le \varepsilon |\log x|` for
    `2^{-1022} \le x < 2^{1022}` with `\varepsilon` about `2^{-49.3}`,
    `2^{-99.3}`, `2^{-150.2}` and `2^{-200.8}`, and `|` ``dN_atan(x)``
    `- \operatorname{atan} x| \le \varepsilon |\operatorname{atan} x|`
    for all finite `x` with `\varepsilon` about `2^{-49.8}`,
    `2^{-99.6}`, `2^{-150.5}` and `2^{-201.4}`. Both bounds are relative
    everywhere, in particular near `x = 1` for the logarithm and near
    `0` for the arctangent, where the reductions are exact. They are
    derived in ``dev/dfloat_exp_bound.py``, which also models the
    long division; the observed worst case on random arguments is
    10-30 times smaller.

.. function:: void d1_log(d1_t res, const d1_t x)
              void d2_log(d2_t res, const d2_t x)
              void d3_log(d3_t res, const d3_t x)
              void d4_log(d4_t res, const d4_t x)

    The natural logarithm, with the guarantee above. With `x = 2^e m`,
    `m \in [M_0, 2 M_0)`, `M_0 = 0.70703125` (the exponent and the
    table index are read off the bits of the head), the reduction is
    `z = m (r r') - 1` with two tabulated reciprocals: `r \approx 1/m`
    from 128 entries indexed by the leading bits of `m` (`r = 1` on
    the two intervals adjacent to 1, which keeps `\log x` relative
    near `x = 1`), giving `|m r - 1| \le 2^{-7}`, and `r' \approx 1/(1
    + j/2^{13})`, `j = \operatorname{round}(2^{13} (m r - 1))` from a
    double approximation, giving `|z| \le 2^{-14} (1 + 2^{-7})`. The
    product `r r'` is exact as a double-double, and `z = (h_0 - 1) +
    l_0 + (m_1 + \cdots) R_h + m R_l` with `m_0 R_h = h_0 + l_0` is
    computed with the head cancelled exactly and the other terms at
    full precision, so that its error is `u^N |z|` rather than `u^N`;
    no division is involved. The tables hold `L = -\log r` for the
    rounded entries, so that `\log x = e \log 2 + L + L' +
    \operatorname{log1p}(z)` exactly. With `w = z^2`,
    `\operatorname{log1p}(z) = z + w E(w) + z w O(w)` with `E(w) = -1/2
    - w/4 - \cdots`, `O(w) = 1/3 + w/5 + \cdots` (degree 8, 8, 12, 15
    in `z`), evaluated as four chains in `v = w^2` in the SIMD lanes,
    and the result is the five-term sum of `D = e \log 2 + L + L' + z`
    (assembled while the chains run) and the four products. Below
    `|z| < 2^{-200}`, `\operatorname{log1p}(z) = z`. For `N = 1` the
    kernel is plain double code with the first level only and a
    polynomial of degree 8 (as fast as ``libm``). Arguments outside
    `[2^{-1022}, 2^{1022})` (subnormal or huge) go through
    :type:`arb_t`; `x \le 0` gives ``-inf`` (zero) or ``nan``,
    ``inf`` gives ``inf``, and `\log 1 = 0` exactly.

.. function:: void d1_atan(d1_t res, const d1_t x)
              void d2_atan(d2_t res, const d2_t x)
              void d3_atan(d3_t res, const d3_t x)
              void d4_atan(d4_t res, const d4_t x)

    The arctangent, with the guarantee above. For `x \ge 0` (odd
    symmetry), with `c = i/64`: for `x \le 1`, `i = \operatorname{round}(64 x)`
    and `\operatorname{atan} x = \operatorname{atan} c +
    \operatorname{atan} z`, `z = (x - c)/(1 + x c)`; for `x > 1`,
    `i = \operatorname{round}(64/x)` and `\operatorname{atan} x =
    (\pi/2 - \operatorname{atan} c) - \operatorname{atan} z`, `z = (1 -
    c x)/(x + c)`; in both cases `|z| \le 2^{-7}`, the numerator's
    head is exact, and the division is the long division kernel. For
    `x \le 1/128` no division is done (`z = x`). With `w = z^2`,
    `\operatorname{atan} z = z + z w E(w)`, `E(w) = -1/3 + w/5 -
    \cdots` with 3, 7, 11, 14 coefficients, evaluated as four chains
    in `w^4` in the SIMD lanes (the multipliers `z w`, `z w^2`, `z
    w^3`, `z w^4` are one vector product), and the result is the
    five-term sum of `\operatorname{atan} c + z` and the four
    products. For `|x| < 2^{-200}` the result is `x`, for `|x| \ge
    2^{900}` it is `\pi/2 - 1/x_0`. For `N = 1` the kernel is plain
    double code with a double division. ``inf`` gives `\pm \pi/2`,
    ``nan`` gives ``nan``, and `\operatorname{atan} 0 = 0` exactly.

.. function:: int d1b_log(d1b_t res, const d1b_t x)
              int d2b_log(d2b_t res, const d2b_t x)
              int d3b_log(d3b_t res, const d3b_t x)
              int d4b_log(d4b_t res, const d4b_t x)
              void d1b_atan(d1b_t res, const d1b_t x)
              void d2b_atan(d2b_t res, const d2b_t x)
              void d3b_atan(d3b_t res, const d3b_t x)
              void d4b_atan(d4b_t res, const d4b_t x)

    The logarithm and arctangent on balls `[m \pm r]`: the kernels on
    `m` and the radii `\varepsilon |y| + r / (m - r)` (the derivative
    `1/t` is largest at the left end of the ball; a lower bound on
    `m - r` is computed from the head, its tail and the radius) and
    `\varepsilon |y| + r` (the derivative is at most 1), with the usual
    slack and underflow bump. ``dNb_log`` returns ``GR_DOMAIN`` if the
    ball is entirely nonpositive and ``GR_UNABLE`` if it is not
    contained in `(0, \infty)` otherwise. Infinite
    radii, nonfinite midpoints and midpoints beyond the kernels'
    ranges go through :type:`arb_t`.

.. macro:: DFLOAT_EXPM1_EPS_1
           DFLOAT_EXPM1_EPS_2
           DFLOAT_EXPM1_EPS_3
           DFLOAT_EXPM1_EPS_4
           DFLOAT_EXPM1_EPS(n)
           DFLOAT_LOG1P_EPS_1
           DFLOAT_LOG1P_EPS_2
           DFLOAT_LOG1P_EPS_3
           DFLOAT_LOG1P_EPS_4
           DFLOAT_LOG1P_EPS(n)

    The static relative error bounds of ``dN_expm1`` (for `-746 \le x
    \le 709.9`; about `2^{-48.9}`, `2^{-98.3}`, `2^{-149.2}`,
    `2^{-199.8}`) and ``dN_log1p`` (for `1 + x \in [2^{-1022},
    2^{1022})`; about `2^{-49.3}`, `2^{-99.1}`, `2^{-149.9}`,
    `2^{-200.5}`), relative everywhere and in particular near 0, where
    both are `x (1 + O(x))`.

.. function:: void d1_expm1(d1_t res, const d1_t x)
              void d2_expm1(d2_t res, const d2_t x)
              void d3_expm1(d3_t res, const d3_t x)
              void d4_expm1(d4_t res, const d4_t x)
              void d1_log1p(d1_t res, const d1_t x)
              void d2_log1p(d2_t res, const d2_t x)
              void d3_log1p(d3_t res, const d3_t x)
              void d4_log1p(d4_t res, const d4_t x)

    `e^x - 1` and `\log(1 + x)` with the guarantees above.
    ``expm1`` runs the exponential kernel for `|x| \le 0.34` (where
    `k = 0`) with the constant term `T` replaced by `T - 1 = E_1 + T_1
    E_2` from tables of `\exp(i/256) - 1` and `\exp(j/65536) - 1`, so
    that the result `C + (T s) E + (T w) O` is relative when `C` is
    small or zero (`C = 0` and `s = x` exactly for `|x| < 2^{-17}`);
    beyond 0.34 it is `\exp(x) - 1`, which loses nothing since `|e^x -
    1| \ge 0.29`. ``log1p`` carries `1 + x` exactly into the logarithm
    kernel for `|x| \le 1/2`: the head `s = \operatorname{fl}(1 + x_0)`
    with its rounding error `t` and the tail of `x` form an
    `(N + 1)`-term argument whose reduction `z = m (r r') - 1` is
    computed to `u^N |z|` as for ``log``; for `|x| > 1/2`, `1 + x` is
    rounded to `N` terms and `|\log(1 + x)| \ge 0.4` absorbs the
    rounding. Both are derived in ``dev/dfloat_exp_bound.py`` (for the non-SIMD
    build the ``expm1`` chains run untapered: the exponential's
    schedule tolerates errors that `s` scales down, which is not the
    case when the result is `s E + w O` itself).

.. function:: void d1_tan(d1_t res, const d1_t x)
              void d2_tan(d2_t res, const d2_t x)
              void d3_tan(d3_t res, const d3_t x)
              void d4_tan(d4_t res, const d4_t x)
              void d1_sinh(d1_t res, const d1_t x)
              void d2_sinh(d2_t res, const d2_t x)
              void d3_sinh(d3_t res, const d3_t x)
              void d4_sinh(d4_t res, const d4_t x)
              void d1_cosh(d1_t res, const d1_t x)
              void d2_cosh(d2_t res, const d2_t x)
              void d3_cosh(d3_t res, const d3_t x)
              void d4_cosh(d4_t res, const d4_t x)
              void d1_tanh(d1_t res, const d1_t x)
              void d2_tanh(d2_t res, const d2_t x)
              void d3_tanh(d3_t res, const d3_t x)
              void d4_tanh(d4_t res, const d4_t x)
              void d1_asin(d1_t res, const d1_t x)
              void d2_asin(d2_t res, const d2_t x)
              void d3_asin(d3_t res, const d3_t x)
              void d4_asin(d4_t res, const d4_t x)
              void d1_acos(d1_t res, const d1_t x)
              void d2_acos(d2_t res, const d2_t x)
              void d3_acos(d3_t res, const d3_t x)
              void d4_acos(d4_t res, const d4_t x)
              void d1_atan2(d1_t res, const d1_t y, const d1_t x)
              void d2_atan2(d2_t res, const d2_t y, const d2_t x)
              void d3_atan2(d3_t res, const d3_t y, const d3_t x)
              void d4_atan2(d4_t res, const d4_t y, const d4_t x)
              void d1_pow(d1_t res, const d1_t x, const d1_t y)
              void d2_pow(d2_t res, const d2_t x, const d2_t y)
              void d3_pow(d3_t res, const d3_t x, const d3_t y)
              void d4_pow(d4_t res, const d4_t x, const d4_t y)

    Compositions of the kernels with the arithmetic (the arguments
    are canonicalised first): `\tan x = \sin x / \cos x` (the absolute
    bounds of the parts make it relative away from the poles);
    `\sinh |x| = e (e + 2) / (2 (e + 1))` with `e = \operatorname{expm1}
    |x|` for `|x| < 1` and `(E - 1/E)/2`, `E = e^{|x|}`, beyond;
    `\cosh x = (E + 1/E)/2`; `\tanh |x| = e/(e + 2)` with `e =
    \operatorname{expm1}(2 |x|)` for `|x| < 1` and `1 - 2/(E + 1)`, `E =
    e^{2|x|}`, beyond (no cancellation in either branch); `\operatorname{asin} x =
    \operatorname{atan}(x / \sqrt{(1 - x)(1 + x)})` and `\operatorname{acos} x = 2
    \operatorname{atan} \sqrt{(1 - x)/(1 + x)}` (`1 - x`, `1 + x` are
    exact where they cancel, so both are relative up to `\pm 1`);
    `\operatorname{atan2}(y, x) = \operatorname{atan}(y/x)` adjusted by
    `\pm \pi` for `x < 0` and `\pm \pi/2` for `x = 0` (nonfinite
    arguments as for the double ``atan2``); `x^y = \exp(y \log x)` for
    `x > 0` and an integer `y` with `|y| < 2^{40}` by binary powering
    of `x` or `1/x` (any base), otherwise ``nan`` for `x < 0`. The
    relative errors are the sums of those of the parts (about
    `2^{-49}`, `2^{-98}`, `2^{-149}`, `2^{-199}`), except that the
    error of `\log x` is amplified by `|y \log x|` in ``pow`` (up to
    `2^{9.5}` at the end of the exponential's range) and that
    `\tan`, `\sinh`, `\tanh` near their zeros inherit the absolute
    errors of `\sin`, `\cos` and of `e^{|x|} - 1/e^{|x|}`. The ball
    versions are the same compositions on balls, and ``dNb_tan``,
    ``dNb_asin``, ``dNb_acos``, ``dNb_pow`` and ``dNb_log1p`` return a
    status following the domains above. Where the composition fails
    for a ball inside the domain, ``dNb_asin`` and ``dNb_acos`` (at and
    near `\pm 1`, where `1 - x^2` contains zero) evaluate the monotone
    functions at the exact endpoints through :type:`arb_t`, and
    ``dNb_tan`` (near a pole, where the cosine's absolute error covers
    zero) calls :func:`arb_tan`; ``dNb_pow`` with an exact integer
    exponent of `2^{40}` or more uses :func:`arb_pow_fmpz`. ``dNb_atan2`` always
    returns a finite enclosure: for `x` away from zero it is
    `\operatorname{atan}(y/x)` adjusted by `\pm \pi`, `\pi` when `x
    < 0` and `y = 0` exactly, or `[0 \pm \pi]` when `x < 0` and `y`
    otherwise contains zero (the branch cut); for `x`
    containing zero and `y` away from it, `\pm \pi/2 -
    \operatorname{atan}(x/y)`, an identity valid for all `x`; when
    both contain zero, `[0 \pm \pi]`.

.. function:: void d1b_add_error_d(d1b_t res, double err)
              void d2b_add_error_d(d2b_t res, double err)
              void d3b_add_error_d(d3b_t res, double err)
              void d4b_add_error_d(d4b_t res, double err)

    Adds `|\text{err}|` to the radius (a nonfinite *err* gives the
    whole line).

.. function:: int d1_is_zero(const d1_t x)
              int d2_is_zero(const d2_t x)
              int d3_is_zero(const d3_t x)
              int d4_is_zero(const d4_t x)
              int d1_is_one(const d1_t x)
              int d2_is_one(const d2_t x)
              int d3_is_one(const d3_t x)
              int d4_is_one(const d4_t x)
              int d1_is_finite(const d1_t x)
              int d2_is_finite(const d2_t x)
              int d3_is_finite(const d3_t x)
              int d4_is_finite(const d4_t x)
              int d1_is_nan(const d1_t x)
              int d2_is_nan(const d2_t x)
              int d3_is_nan(const d3_t x)
              int d4_is_nan(const d4_t x)
              int d1_equal(const d1_t x, const d1_t y)
              int d2_equal(const d2_t x, const d2_t y)
              int d3_equal(const d3_t x, const d3_t y)
              int d4_equal(const d4_t x, const d4_t y)
              int d1_cmp(const d1_t x, const d1_t y)
              int d2_cmp(const d2_t x, const d2_t y)
              int d3_cmp(const d3_t x, const d3_t y)
              int d4_cmp(const d4_t x, const d4_t y)
              int d1_sgn(const d1_t x)
              int d2_sgn(const d2_t x)
              int d3_sgn(const d3_t x)
              int d4_sgn(const d4_t x)

    Predicates and exact comparison on the plain types (``cmp`` compares
    the exact values).

.. function:: int d1b_is_exact(const d1b_t x)
              int d2b_is_exact(const d2b_t x)
              int d3b_is_exact(const d3b_t x)
              int d4b_is_exact(const d4b_t x)
              int d1b_is_finite(const d1b_t x)
              int d2b_is_finite(const d2b_t x)
              int d3b_is_finite(const d3b_t x)
              int d4b_is_finite(const d4b_t x)
              int d1b_is_zero(const d1b_t x)
              int d2b_is_zero(const d2b_t x)
              int d3b_is_zero(const d3b_t x)
              int d4b_is_zero(const d4b_t x)
              int d1b_contains_zero(const d1b_t x)
              int d2b_contains_zero(const d2b_t x)
              int d3b_contains_zero(const d3b_t x)
              int d4b_contains_zero(const d4b_t x)
              int d1b_contains_d(const d1b_t x, double c)
              int d2b_contains_d(const d2b_t x, double c)
              int d3b_contains_d(const d3b_t x, double c)
              int d4b_contains_d(const d4b_t x, double c)
              int d1b_contains(const d1b_t x, const d1b_t y)
              int d2b_contains(const d2b_t x, const d2b_t y)
              int d3b_contains(const d3b_t x, const d3b_t y)
              int d4b_contains(const d4b_t x, const d4b_t y)
              int d1b_overlaps(const d1b_t x, const d1b_t y)
              int d2b_overlaps(const d2b_t x, const d2b_t y)
              int d3b_overlaps(const d3b_t x, const d3b_t y)
              int d4b_overlaps(const d4b_t x, const d4b_t y)
              int d1b_equal(const d1b_t x, const d1b_t y)
              int d2b_equal(const d2b_t x, const d2b_t y)
              int d3b_equal(const d3b_t x, const d3b_t y)
              int d4b_equal(const d4b_t x, const d4b_t y)
              int d1b_cmp(int * res, const d1b_t x, const d1b_t y)
              int d2b_cmp(int * res, const d2b_t x, const d2b_t y)
              int d3b_cmp(int * res, const d3b_t x, const d3b_t y)
              int d4b_cmp(int * res, const d4b_t x, const d4b_t y)

    Ball predicates. ``is_exact`` is true for a finite ball of radius
    zero and ``is_finite`` for a ball that is not the whole line.
    ``is_zero``, ``equal``, ``contains`` and ``contains_d`` answer 1
    only when certain (``is_zero`` and ``equal`` only for exact balls,
    ``contains_d`` only for a finite *c*); ``contains_zero`` and
    ``overlaps`` answer 0 only when certain (a doubtful case counts as
    containing or overlapping). The whole line contains everything and
    overlaps everything, and is contained only in the whole line.
    ``cmp`` sets `-1, 0, 1` and returns ``GR_SUCCESS`` when the balls
    decide the comparison (or both are the same exact number) and
    ``GR_UNABLE`` otherwise.

.. function:: double d1_get_d(const d1_t x)
              double d2_get_d(const d2_t x)
              double d3_get_d(const d3_t x)
              double d4_get_d(const d4_t x)
              double d1b_get_d(const d1b_t x)
              double d2b_get_d(const d2b_t x)
              double d3b_get_d(const d3b_t x)
              double d4b_get_d(const d4b_t x)

    The leading component (the double nearest to the value, up to
    ties); 0 for the whole line.

.. function:: void d1b_get_mid_d(double * mid, const d1b_t x)
              void d2b_get_mid_d(double * mid, const d2b_t x)
              void d3b_get_mid_d(double * mid, const d3b_t x)
              void d4b_get_mid_d(double * mid, const d4b_t x)

    Sets *mid* to the `N` midpoint components (zeros for the whole
    line).

.. function:: void d1_get_arf(arf_t res, const d1_t x)
              void d2_get_arf(arf_t res, const d2_t x)
              void d3_get_arf(arf_t res, const d3_t x)
              void d4_get_arf(arf_t res, const d4_t x)
              void d1_set_arf(d1_t res, const arf_t x)
              void d2_set_arf(d2_t res, const arf_t x)
              void d3_set_arf(d3_t res, const arf_t x)
              void d4_set_arf(d4_t res, const arf_t x)
              void d1_get_arb(arb_t res, const d1_t x)
              void d2_get_arb(arb_t res, const d2_t x)
              void d3_get_arb(arb_t res, const d3_t x)
              void d4_get_arb(arb_t res, const d4_t x)
              void d1_set_arb(d1_t res, const arb_t x)
              void d2_set_arb(d2_t res, const arb_t x)
              void d3_set_arb(d3_t res, const arb_t x)
              void d4_set_arb(d4_t res, const arb_t x)
              void d1b_get_arb(arb_t res, const d1b_t x)
              void d2b_get_arb(arb_t res, const d2b_t x)
              void d3b_get_arb(arb_t res, const d3b_t x)
              void d4b_get_arb(arb_t res, const d4b_t x)
              void d1b_set_arb(d1b_t res, const arb_t x)
              void d2b_set_arb(d2b_t res, const arb_t x)
              void d3b_set_arb(d3b_t res, const arb_t x)
              void d4b_set_arb(d4b_t res, const arb_t x)
              void d1b_set_arf(d1b_t res, const arf_t x)
              void d2b_set_arf(d2b_t res, const arf_t x)
              void d3b_set_arf(d3b_t res, const arf_t x)
              void d4b_set_arf(d4b_t res, const arf_t x)

    Conversions. The ``get`` functions are exact (the radius is rounded
    up to a :type:`mag_t`); the ``set`` functions round the midpoint
    to `N` doubles by greedy rounding and, for balls, add the remainder
    to the radius. A nonfinite :type:`arb_t` or :type:`arf_t` gives the
    whole real line, and ``dNb_get_arb`` gives the whole line as
    `[0 \pm \infty]`.

.. function:: void d1_set_fmpz(d1_t res, const fmpz_t x)
              void d2_set_fmpz(d2_t res, const fmpz_t x)
              void d3_set_fmpz(d3_t res, const fmpz_t x)
              void d4_set_fmpz(d4_t res, const fmpz_t x)
              void d1_set_fmpq(d1_t res, const fmpq_t x)
              void d2_set_fmpq(d2_t res, const fmpq_t x)
              void d3_set_fmpq(d3_t res, const fmpq_t x)
              void d4_set_fmpq(d4_t res, const fmpq_t x)
              void d1b_set_fmpz(d1b_t res, const fmpz_t x)
              void d2b_set_fmpz(d2b_t res, const fmpz_t x)
              void d3b_set_fmpz(d3b_t res, const fmpz_t x)
              void d4b_set_fmpz(d4b_t res, const fmpz_t x)
              void d1b_set_fmpq(d1b_t res, const fmpq_t x)
              void d2b_set_fmpq(d2b_t res, const fmpq_t x)
              void d3b_set_fmpq(d3b_t res, const fmpq_t x)
              void d4b_set_fmpq(d4b_t res, const fmpq_t x)
              int d1_set_str(d1_t res, const char * s)
              int d2_set_str(d2_t res, const char * s)
              int d3_set_str(d3_t res, const char * s)
              int d4_set_str(d4_t res, const char * s)
              int d1b_set_str(d1b_t res, const char * s)
              int d2b_set_str(d2b_t res, const char * s)
              int d3b_set_str(d3b_t res, const char * s)
              int d4b_set_str(d4b_t res, const char * s)
              char * d1_get_str(const d1_t x, slong digits)
              char * d2_get_str(const d2_t x, slong digits)
              char * d3_get_str(const d3_t x, slong digits)
              char * d4_get_str(const d4_t x, slong digits)
              char * d1b_get_str(const d1b_t x, slong digits)
              char * d2b_get_str(const d2b_t x, slong digits)
              char * d3b_get_str(const d3b_t x, slong digits)
              char * d4b_get_str(const d4b_t x, slong digits)
              void d1_print(const d1_t x)
              void d2_print(const d2_t x)
              void d3_print(const d3_t x)
              void d4_print(const d4_t x)
              void d1b_print(const d1b_t x)
              void d2b_print(const d2b_t x)
              void d3b_print(const d3b_t x)
              void d4b_print(const d4b_t x)

    Conversions from integers, rationals and decimal strings (through
    :type:`arb_t`, rigorously for the balls), and printing through
    :func:`arb_get_str` (``digits = 0`` selects the full precision).

.. function:: void d1_randtest(d1_t res, flint_rand_t state)
              void d2_randtest(d2_t res, flint_rand_t state)
              void d3_randtest(d3_t res, flint_rand_t state)
              void d4_randtest(d4_t res, flint_rand_t state)
              void d1_randtest_special(d1_t res, flint_rand_t state)
              void d2_randtest_special(d2_t res, flint_rand_t state)
              void d3_randtest_special(d3_t res, flint_rand_t state)
              void d4_randtest_special(d4_t res, flint_rand_t state)
              void d1b_randtest(d1b_t res, flint_rand_t state)
              void d2b_randtest(d2b_t res, flint_rand_t state)
              void d3b_randtest(d3b_t res, flint_rand_t state)
              void d4b_randtest(d4b_t res, flint_rand_t state)
              void d1b_randtest_special(d1b_t res, flint_rand_t state)
              void d2b_randtest_special(d2b_t res, flint_rand_t state)
              void d3b_randtest_special(d3b_t res, flint_rand_t state)
              void d4b_randtest_special(d4b_t res, flint_rand_t state)

    Random values; the ``special`` versions cover the whole exponent
    range including subnormals and the ends of the double range, and
    radii from zero to infinity.

.. function:: void _d1_vec_add(d1_ptr res, d1_srcptr x, d1_srcptr y, slong len)
              void _d2_vec_add(d2_ptr res, d2_srcptr x, d2_srcptr y, slong len)
              void _d3_vec_add(d3_ptr res, d3_srcptr x, d3_srcptr y, slong len)
              void _d4_vec_add(d4_ptr res, d4_srcptr x, d4_srcptr y, slong len)
              void _d1_vec_sub(d1_ptr res, d1_srcptr x, d1_srcptr y, slong len)
              void _d2_vec_sub(d2_ptr res, d2_srcptr x, d2_srcptr y, slong len)
              void _d3_vec_sub(d3_ptr res, d3_srcptr x, d3_srcptr y, slong len)
              void _d4_vec_sub(d4_ptr res, d4_srcptr x, d4_srcptr y, slong len)
              void _d1_vec_mul(d1_ptr res, d1_srcptr x, d1_srcptr y, slong len)
              void _d2_vec_mul(d2_ptr res, d2_srcptr x, d2_srcptr y, slong len)
              void _d3_vec_mul(d3_ptr res, d3_srcptr x, d3_srcptr y, slong len)
              void _d4_vec_mul(d4_ptr res, d4_srcptr x, d4_srcptr y, slong len)
              void _d1_vec_mul_scalar(d1_ptr res, d1_srcptr x, slong len, d1_srcptr c)
              void _d2_vec_mul_scalar(d2_ptr res, d2_srcptr x, slong len, d2_srcptr c)
              void _d3_vec_mul_scalar(d3_ptr res, d3_srcptr x, slong len, d3_srcptr c)
              void _d4_vec_mul_scalar(d4_ptr res, d4_srcptr x, slong len, d4_srcptr c)
              void _d1_vec_addmul_scalar(d1_ptr res, d1_srcptr x, slong len, d1_srcptr c)
              void _d2_vec_addmul_scalar(d2_ptr res, d2_srcptr x, slong len, d2_srcptr c)
              void _d3_vec_addmul_scalar(d3_ptr res, d3_srcptr x, slong len, d3_srcptr c)
              void _d4_vec_addmul_scalar(d4_ptr res, d4_srcptr x, slong len, d4_srcptr c)
              void _d1_vec_submul_scalar(d1_ptr res, d1_srcptr x, slong len, d1_srcptr c)
              void _d2_vec_submul_scalar(d2_ptr res, d2_srcptr x, slong len, d2_srcptr c)
              void _d3_vec_submul_scalar(d3_ptr res, d3_srcptr x, slong len, d3_srcptr c)
              void _d4_vec_submul_scalar(d4_ptr res, d4_srcptr x, slong len, d4_srcptr c)
              void _d1_vec_dot(d1_t res, d1_srcptr initial, int subtract, d1_srcptr x, d1_srcptr y, slong len)
              void _d2_vec_dot(d2_t res, d2_srcptr initial, int subtract, d2_srcptr x, d2_srcptr y, slong len)
              void _d3_vec_dot(d3_t res, d3_srcptr initial, int subtract, d3_srcptr x, d3_srcptr y, slong len)
              void _d4_vec_dot(d4_t res, d4_srcptr initial, int subtract, d4_srcptr x, d4_srcptr y, slong len)
              void _d1_vec_dot_rev(d1_t res, d1_srcptr initial, int subtract, d1_srcptr x, d1_srcptr y, slong len)
              void _d2_vec_dot_rev(d2_t res, d2_srcptr initial, int subtract, d2_srcptr x, d2_srcptr y, slong len)
              void _d3_vec_dot_rev(d3_t res, d3_srcptr initial, int subtract, d3_srcptr x, d3_srcptr y, slong len)
              void _d4_vec_dot_rev(d4_t res, d4_srcptr initial, int subtract, d4_srcptr x, d4_srcptr y, slong len)
              void _d1_vec_div(d1_ptr res, d1_srcptr x, d1_srcptr y, slong len)
              void _d2_vec_div(d2_ptr res, d2_srcptr x, d2_srcptr y, slong len)
              void _d3_vec_div(d3_ptr res, d3_srcptr x, d3_srcptr y, slong len)
              void _d4_vec_div(d4_ptr res, d4_srcptr x, d4_srcptr y, slong len)
              void _d1_vec_sqrt(d1_ptr res, d1_srcptr x, slong len)
              void _d2_vec_sqrt(d2_ptr res, d2_srcptr x, slong len)
              void _d3_vec_sqrt(d3_ptr res, d3_srcptr x, slong len)
              void _d4_vec_sqrt(d4_ptr res, d4_srcptr x, slong len)
              void _d1_vec_rsqrt(d1_ptr res, d1_srcptr x, slong len)
              void _d2_vec_rsqrt(d2_ptr res, d2_srcptr x, slong len)
              void _d3_vec_rsqrt(d3_ptr res, d3_srcptr x, slong len)
              void _d4_vec_rsqrt(d4_ptr res, d4_srcptr x, slong len)
              void _d1_vec_exp(d1_ptr res, d1_srcptr x, slong len)
              void _d2_vec_exp(d2_ptr res, d2_srcptr x, slong len)
              void _d3_vec_exp(d3_ptr res, d3_srcptr x, slong len)
              void _d4_vec_exp(d4_ptr res, d4_srcptr x, slong len)
              void _d1_vec_sin(d1_ptr res, d1_srcptr x, slong len)
              void _d2_vec_sin(d2_ptr res, d2_srcptr x, slong len)
              void _d3_vec_sin(d3_ptr res, d3_srcptr x, slong len)
              void _d4_vec_sin(d4_ptr res, d4_srcptr x, slong len)
              void _d1_vec_cos(d1_ptr res, d1_srcptr x, slong len)
              void _d2_vec_cos(d2_ptr res, d2_srcptr x, slong len)
              void _d3_vec_cos(d3_ptr res, d3_srcptr x, slong len)
              void _d4_vec_cos(d4_ptr res, d4_srcptr x, slong len)
              void _d1_vec_sin_cos(d1_ptr sres, d1_ptr cres, d1_srcptr x, slong len)
              void _d2_vec_sin_cos(d2_ptr sres, d2_ptr cres, d2_srcptr x, slong len)
              void _d3_vec_sin_cos(d3_ptr sres, d3_ptr cres, d3_srcptr x, slong len)
              void _d4_vec_sin_cos(d4_ptr sres, d4_ptr cres, d4_srcptr x, slong len)
              void _d1_vec_log(d1_ptr res, d1_srcptr x, slong len)
              void _d2_vec_log(d2_ptr res, d2_srcptr x, slong len)
              void _d3_vec_log(d3_ptr res, d3_srcptr x, slong len)
              void _d4_vec_log(d4_ptr res, d4_srcptr x, slong len)
              void _d1_vec_atan(d1_ptr res, d1_srcptr x, slong len)
              void _d2_vec_atan(d2_ptr res, d2_srcptr x, slong len)
              void _d3_vec_atan(d3_ptr res, d3_srcptr x, slong len)
              void _d4_vec_atan(d4_ptr res, d4_srcptr x, slong len)

    Vector operations (also for ``dNb``; the ball ``_dNb_vec_log``,
    ``_dNb_vec_sqrt`` and ``_dNb_vec_rsqrt`` return the or of the
    statuses of the elements, whose results are unspecified where
    they fail), used by the generic ring. With SIMD (see the
    requirements above), the additions,
    products, dot products, divisions, square roots and elementary
    functions process four elements at a time; the elementary
    functions run the same kernels on the four lanes (with the
    two-way split of the polynomials), so that the static bounds
    apply to every element, and redo in scalar the elements that hit
    a special case (``vec_log`` and ``vec_atan`` for `N = 1` use the
    scalar kernels, which are as fast). The vector division is the
    long division on four lanes (1.3, 4.5 and 11 ns per element for
    `N = 2, 3, 4`, against 7.6, 23 and 45 ns scalar); the elements
    whose divisor is zero, nonfinite, or tiny enough for the residual
    products to underflow are redone in scalar. The vector square
    roots run the Newton iterations of ``sqrt`` and ``rsqrt`` on four
    lanes (a lane whose remainder is left with a zero head by an
    exact cancellation, or whose result is not canonical, is redone
    in scalar), and the ball versions compute the residual bounds of
    the scalar functions in the lanes. The ball products by a scalar
    (``vec_mul_scalar``, ``vec_addmul_scalar``, ``vec_submul_scalar``)
    run the SIMD product with the scalar broadcast (and the SIMD sum),
    with the same results as the scalar operations.

    The dot products keep four lane accumulators in SIMD registers and
    combine them at the end; ``vec_dot_rev`` (the kernel of the
    classical polynomial multiplication and of the power series
    basecases: inverse, division, square root, exponential, ...) runs
    the same loop with the blocks of *y* loaded in reverse lane order,
    so that it gives the same result as ``vec_dot`` on a reversed copy
    of *y*. Short dot products (below 8 terms for the balls, 10 for
    ``d2`` and ``d3``, 32 for ``d1``, where the SIMD setup and the final
    combination of the lanes do not pay off) take a lean scalar loop.
    Classical multiplication of two polynomials of length 256 costs, per
    coefficient product, 0.54, 3.4, 9.0, 18 ns for ``d1`` to ``d4`` and
    3.5, 9.2, 16, 25 ns for ``d1b`` to ``d4b`` (against 1.0, 12, 25, 61
    and 8.4, 9.7, 16, 27 ns with a scalar ``vec_dot_rev`` loop);
    the gains start at length 16, and the short
    lengths are within a few percent of the scalar loop. Nanoseconds per element, one
    run on a machine going about twice slower than in the tables
    above (the ``vec_div`` and ``vec_exp`` columns for reference):

    ====== ======== ============ ========= ============= ========== ============== =========== =============== =========== =============
    format ``sqrt`` ``vec_sqrt`` ``rsqrt`` ``vec_rsqrt`` ``b sqrt`` ``b vec_sqrt`` ``b rsqrt`` ``b vec_rsqrt`` ``vec_div`` ``b vec_exp``
    ====== ======== ============ ========= ============= ========== ============== =========== =============== =========== =============
    d1          2.1          0.9       3.4           1.5         10            3.1          18             5.2         3.1           7.9
    d2           12          3.7        10           3.1         48             16          57              18         6.4            35
    d3           79           19        88            30        135             38         196              56          13            86
    d4          122           34       150            45        227             53         312              89          28           150
    ====== ======== ============ ========= ============= ========== ============== =========== =============== =========== =============

Complex numbers
-------------------------------------------------------------------------------

The complex formats ``d1c`` ... ``d4c`` (types ``dNc_t``) are pairs of
``dN`` numbers, the real and the imaginary part, and the complex ball
formats ``d1cb`` ... ``d4cb`` (types ``dNcb_t``) are pairs of ``dNb``
balls: rectangular balls, as in :type:`acb_t`, representing every
`a + b i` with `a` in the real ball and `b` in the imaginary ball. The
parts are the structure members ``re`` and ``im`` (with the accessor
macros ``dNc_realref``, ``dNc_imagref`` etc.) and are stored
consecutively, so that an array of `n` complex numbers is also an
array of `2n` real numbers of the same format.

**Conventions.** Each part follows the conventions of the real balls:
a part with a nonfinite component or radius is the whole real line
`W`, so that `W + 2i`, `3 + W i` and `W + W i` are the sets `\mathbb{R}
+ 2i`, `3 + \mathbb{R} i` and `\mathbb{C}`, and results use the normal
form `[0 \pm \infty]` for such parts. The complex ball rings implement
the field `\mathbb{C}` with the status rules of the real rings:
operations defined on all of `\mathbb{C}` (``add``, ``sub``, ``mul``,
``sqr``, ``neg``, ``conj``, ``mul_2exp_si``, ``exp``, ``sin``,
``cos``, ``sin_cos``, ``sinh``, ``cosh``, ``sqrt``, ``abs``, ``arg``,
``sgn``, the dot products and integer powers with a nonnegative
exponent) always succeed; the others return ``GR_SUCCESS`` only for
inputs inside their domain, ``GR_DOMAIN`` only for inputs outside
it, and ``GR_UNABLE`` otherwise. The domains: ``div`` and ``inv``, a
nonzero divisor; ``rsqrt`` and ``log``, `z \ne 0`; ``tan``, `z \notin
\pi/2 + \pi \mathbb{Z}` and ``tanh``, `z \notin i (\pi/2 + \pi
\mathbb{Z})` (never ``GR_DOMAIN``, the poles being irrational);
``pow``, `x \ne 0`, or `x = 0` and `\operatorname{Re}(y) > 0`, or `y
\in \{0, 1, 2, \ldots\}`. The functions with branch cuts (``log``,
``sqrt``, ``rsqrt``, ``pow``, ``arg``) are the principal branches,
with the cut `(-\infty, 0]` continuous from above (`\log(-1) = \pi i`,
`\sqrt{-4} = 2i`), and a ball crossing the cut gets an enclosure of
the values on both sides. For example:

    ===================================== ===============================================
    operation                             status and result
    ===================================== ===============================================
    `W / 3`                               ``GR_SUCCESS``, `W`
    `1 / 0`, `W / 0`                      ``GR_DOMAIN``
    `1 / ([0 \pm 1] + [0 \pm 1] i)`       ``GR_UNABLE``
    `1 / (W + 2i)`                        ``GR_SUCCESS``, `[0 \pm 1/2] + [0 \pm 1/2] i`
    `\exp(W)`                             ``GR_SUCCESS``, `W + 0 i`
    `\exp(W i)`                           ``GR_SUCCESS``, `[0 \pm 1] + [0 \pm 1] i`
    `\sqrt{-4}`                           ``GR_SUCCESS``, `2 i`
    `\sqrt{[0 \pm 1] + [0 \pm 1] i}`      ``GR_SUCCESS`` (an enclosure)
    `\log 0`                              ``GR_DOMAIN``
    `\log(-1 + [0 \pm 1/4] i)`            ``GR_SUCCESS``, imaginary part `[0 \pm \pi]`
    `\log(W + 2i)`                        ``GR_SUCCESS``, real part `W`
    `\tan(\pi/2 + 10^{-3} i)`             ``GR_SUCCESS``
    `\tan([\pi/2 \pm 0.01])`              ``GR_UNABLE``
    `0^{2 + i}`, `0^0`                    ``GR_SUCCESS``, 0 and 1
    `0^{-1}`, `0^i`                       ``GR_DOMAIN``
    `[0 \pm 1/2]^{1/2}`                   ``GR_SUCCESS`` (an enclosure)
    ===================================== ===============================================

The plain complex types follow the plain real types: results are
approximate (normwise accurate to a few ulps where the functions are
well conditioned), nonfinite values propagate, and the generic ring
returns ``GR_DOMAIN`` only for an exact singularity (a zero divisor,
`1/\sqrt 0`, `0^y` with `\operatorname{Re}(y) \le 0`, `y \ne 0`); ``log
0`` is `-\infty` as for the real type.

The conversions of the generic rings: from the real types (and
:type:`arb_t`, ``fmpz``, ``fmpq``, the real ``dfloat`` rings) through
the real ring of the parts, with a zero imaginary part; from
:type:`acb_t` part by part (a nonfinite part gives `W`); between
complex formats part by part as between real formats; from any other
ring through :type:`acb_t`. To a real type (``get_fmpz``, ``get_d``,
``gr_set_other`` into a real ring or :type:`arb_t`), the imaginary part
must be exactly zero: ``GR_DOMAIN`` if it is not zero (or certainly
not, for a ball) and ``GR_UNABLE`` if a ball only contains zero. The
:type:`arb_t` and :type:`acb_t` rings convert from the complex
``dfloat`` rings in ``gr_set_other``. Elements print as in the
:type:`acb_t` ring: ``re``, ``im*I`` or ``(re + im*I)``.

**Implementation.** The complex code has no kernels of its own and is
compiled per format only where speed matters:

* The arithmetic (``add``, ``sub``, ``mul``, ``sqr``, ``mul_re``,
  ``conj``, ...) calls the real operations, which the compiler inlines
  where they are small (for `N = 1`). The product is `(a c - b d) + (a
  d + b c) i`, or two real products when either factor is exactly real
  or imaginary (which is also tighter for balls); the square is `(a +
  b)(a - b) + 2 a b i` for the plain types (a relative real part) and
  `(a^2 - b^2) + 2 a b i` for the balls (a smaller radius).
* The vector operations reuse the real vector operations: ``vec_add``
  and ``vec_sub`` are the real ones on `2 \cdot` ``len`` elements, and so
  are the products by a scalar that is exactly real; ``vec_exp`` works
  on blocks of 32 elements, with the real ``vec_exp`` of the real
  parts, ``vec_sin_cos`` of the imaginary parts and two ``vec_mul``.
* The dot products are two real dot products per block of 32 terms:
  with `X` the `2 \cdot` ``len`` real numbers of `x`, the real part is
  the real dot product of `X` with the copy `(\operatorname{Re} y_j,
  -\operatorname{Im} y_j)_j` and the imaginary part that of `X` with
  `(\operatorname{Im} y_j, \operatorname{Re} y_j)_j` (with `y` reversed
  for ``vec_dot_rev``), the result of a block being the initial value of
  the next; they run the SIMD real dot products and cost 4-5 times a
  real dot product per term (the four real products), 3-4 times less
  than a loop of complex products and sums. (For ``d1c`` a direct
  loop over the doubles is used instead, since the copies would cost
  more than the products.)
* Everything else (division, square roots, the elementary functions,
  the conversions, strings, random elements and the corresponding gr
  methods) is implemented once for all formats in
  ``src/dfloat/complex.c``, through a table of the real operations of
  the format (``_dfloat_ops_struct``: pointers to the real functions,
  including the ball products with a run-time error mode and the
  format's complex product), so that it exists once instead of eight
  times, compiled for size. Plain and ball formats share the code; the
  checks that only concern balls are made where the table says the
  format is a ball.

The functions of ``complex.c``: the quotient is `x \bar y / |y|^2`
after scaling `y` by a power of two (so that `|y|^2` neither
overflows nor underflows), or two real divisions for a real or
imaginary `y`; for balls, `|y|^2` gets the lower bound `(\min |a|)^2 +
(\min |b|)^2` where the squares of wide balls give less, and where the
formula still fails for a divisor away from zero (a wide rectangle,
or one with a whole-line part), the result is the crude enclosure
`|x/y| \le \max|x| / \min|y|` of both parts. `|z|` is the square root
of the scaled `a^2 + b^2`. The square root is `t + u i` with `t =
\sqrt{(|z| + a)/2}` and `u = b/(2t)` for `a \ge 0`, and `|u| =
\sqrt{(|z| - a)/2}` with the sign of `b` and `t = |b|/(2|u|)` for `a <
0` (neither form cancels); for balls these hold at every point as
long as the divisor does not vanish, i.e. away from the cut, a ball
with `a < 0` crossing the cut gets `u = [0 \pm |u|]`, and a ball
touching the cut or zero gets the enclosure `[0, s] + [-s, s] i`, `s =
\sqrt{\max |z|}`. `\exp(a + bi) = e^a (\cos b + i \sin b)`; `\log z =
\log|z| + i \operatorname{atan2}(b, a)` with `\log|z| =
\operatorname{log1p}((p - 1)(p + 1) + q^2)/2` (`p` the larger part)
where `|z|^2` is near 1, so that it stays relative, and `\log(|z
2^{-e}|^2)/2 + e \log 2` otherwise (for a wide ball where that fails,
the interval `[\log \min|z|, \log \max|z|]`); `\sin`, `\cos`, `\sinh`,
`\cosh` by the addition formulas with the real sine, cosine and the
real `\sinh`, `\cosh` from one exponential; `\tan z = (t + u i)/(1 -
t u i)` with `t = \tan a`, `u = \tanh b` (the denominator has real
part 1), or, where `\tan a` fails (a real part near a pole or `W`) and
`b` is away from zero, `i (1 - w)/(1 + w)` with `w = e^{2iz}` taken with
`|w| < 1` (by conjugation), or the enclosure `|\tan z| \le \coth |b|
< 1/|b| + 1`; `\tanh z = -i \tan(iz)`; `x^y` by binary powering for an
integer `y` with `|y| < 2^{40}` and `\exp(y \log x)` otherwise, with the
enclosure `|x^y| \le \max(1, |x|)^{\operatorname{Re} y} e^{\pi
|\operatorname{Im} y|}` of both parts for `x` containing zero and
`\operatorname{Re}(y) > 0`. The ball results are rigorous by
composition. With the crude enclosures, a complex ball function
returns ``GR_UNABLE`` for an input inside its domain only for 0.2-0.5%
of the fuzzer's inputs (balls whose lower bound on `|z|` is below
`2^{-1000}`, or tangents at the edge of what the double range can
decide).

**Performance.** Nanoseconds per call through the generic ring (one
core, the machine of the tables above; the absolute numbers vary by
about 1.5x between runs), against :type:`acb_t` at `53 N` bits:

    ======= ===== ===== ====== ===== ====== ===== ====== ======
    ring      mul   div   sqrt   exp    log   sin    tan    pow
    ======= ===== ===== ====== ===== ====== ===== ====== ======
    d1cb       35   118    113    41     99    94    162    182
    acb 53     86   360    512   652    766  1042   1362   1055
    d2cb       43   148    170   202    270   313    457    544
    acb 106   103   386    447   492    612   771   1329   1287
    d3cb       90   303    448   413    620   650    920   1174
    acb 159   120   576    670   687    724  1016   1433   1619
    d4cb      119   412    478   695    971   960   1232   1617
    acb 212   104   418    529   706    694   952   1546   1581
    ======= ===== ===== ====== ===== ====== ===== ====== ======

The balls are several times faster than :type:`acb_t` for `N \le 2`,
faster for `N = 3`, and on par for `N = 4`, where the four real
products of a complex product (with their error tracking) cost as
much as the corresponding :type:`acb_t` product; the elementary
functions are compositions of the real ones and cost about the sum of
their parts. The dot products, per term, for length 64 (``dN dot``
and ``dNb dot`` the real ones for comparison, ``loop`` a loop of
complex ``gr_mul`` and ``gr_add``):

    ====== ======= ======== ======== ======= ======== =========
    format dN dot  dNc dot  dNc loop dNb dot dNcb dot dNcb loop
    ====== ======= ======== ======== ======= ======== =========
    d1        0.36     1.06     3.99    1.69      7.8      23.5
    d2        1.83     8.30     32.1    4.86     21.6      65.0
    d3        5.06     22.1     69.7    9.88     36.7      95.7
    d4        10.0     42.8      144    16.5     61.7       197
    ====== ======= ======== ======== ======= ======== =========

**Code size.** The per-format complex code (the arithmetic, the
vectors, the dot products, their gr methods and the thin wrappers of
the functions of ``complex.c``, and the table of real operations)
adds about 27 KB of text per `N`, and ``complex.c`` 22 KB, together
about 13% of the real code.

Types and functions
...............................................................................

The functions exist for each format: ``dNc`` (plain, with the real
type ``dN``) and ``dNcb`` (balls, with the real type ``dNb``) for
`N = 1, 2, 3, 4`. Aliasing is allowed everywhere.

.. function:: void d1c_zero(d1c_t res)
              void d2c_zero(d2c_t res)
              void d3c_zero(d3c_t res)
              void d4c_zero(d4c_t res)
              void d1c_one(d1c_t res)
              void d2c_one(d2c_t res)
              void d3c_one(d3c_t res)
              void d4c_one(d4c_t res)
              void d1c_onei(d1c_t res)
              void d2c_onei(d2c_t res)
              void d3c_onei(d3c_t res)
              void d4c_onei(d4c_t res)
              void d1c_set(d1c_t res, const d1c_t x)
              void d2c_set(d2c_t res, const d2c_t x)
              void d3c_set(d3c_t res, const d3c_t x)
              void d4c_set(d4c_t res, const d4c_t x)
              void d1c_set_d(d1c_t res, double x)
              void d2c_set_d(d2c_t res, double x)
              void d3c_set_d(d3c_t res, double x)
              void d4c_set_d(d4c_t res, double x)
              void d1c_set_d_d(d1c_t res, double x, double y)
              void d2c_set_d_d(d2c_t res, double x, double y)
              void d3c_set_d_d(d3c_t res, double x, double y)
              void d4c_set_d_d(d4c_t res, double x, double y)
              void d1c_set_re_im(d1c_t res, const d1_t x, const d1_t y)
              void d2c_set_re_im(d2c_t res, const d2_t x, const d2_t y)
              void d3c_set_re_im(d3c_t res, const d3_t x, const d3_t y)
              void d4c_set_re_im(d4c_t res, const d4_t x, const d4_t y)
              void d1cb_zero(d1cb_t res)
              void d2cb_zero(d2cb_t res)
              void d3cb_zero(d3cb_t res)
              void d4cb_zero(d4cb_t res)
              void d1cb_one(d1cb_t res)
              void d2cb_one(d2cb_t res)
              void d3cb_one(d3cb_t res)
              void d4cb_one(d4cb_t res)
              void d1cb_onei(d1cb_t res)
              void d2cb_onei(d2cb_t res)
              void d3cb_onei(d3cb_t res)
              void d4cb_onei(d4cb_t res)
              void d1cb_indeterminate(d1cb_t res)
              void d2cb_indeterminate(d2cb_t res)
              void d3cb_indeterminate(d3cb_t res)
              void d4cb_indeterminate(d4cb_t res)
              void d1cb_set(d1cb_t res, const d1cb_t x)
              void d2cb_set(d2cb_t res, const d2cb_t x)
              void d3cb_set(d3cb_t res, const d3cb_t x)
              void d4cb_set(d4cb_t res, const d4cb_t x)
              void d1cb_set_d(d1cb_t res, double x)
              void d2cb_set_d(d2cb_t res, double x)
              void d3cb_set_d(d3cb_t res, double x)
              void d4cb_set_d(d4cb_t res, double x)
              void d1cb_set_d_d(d1cb_t res, double x, double y)
              void d2cb_set_d_d(d2cb_t res, double x, double y)
              void d3cb_set_d_d(d3cb_t res, double x, double y)
              void d4cb_set_d_d(d4cb_t res, double x, double y)
              void d1cb_set_re_im(d1cb_t res, const d1b_t x, const d1b_t y)
              void d2cb_set_re_im(d2cb_t res, const d2b_t x, const d2b_t y)
              void d3cb_set_re_im(d3cb_t res, const d3b_t x, const d3b_t y)
              void d4cb_set_re_im(d4cb_t res, const d4b_t x, const d4b_t y)

    Assignments (``indeterminate`` is `W + W i`; ``set`` and
    ``set_re_im`` copy the parts as they are).

.. function:: void d1c_neg(d1c_t res, const d1c_t x)
              void d2c_neg(d2c_t res, const d2c_t x)
              void d3c_neg(d3c_t res, const d3c_t x)
              void d4c_neg(d4c_t res, const d4c_t x)
              void d1c_conj(d1c_t res, const d1c_t x)
              void d2c_conj(d2c_t res, const d2c_t x)
              void d3c_conj(d3c_t res, const d3c_t x)
              void d4c_conj(d4c_t res, const d4c_t x)
              void d1c_add(d1c_t res, const d1c_t x, const d1c_t y)
              void d2c_add(d2c_t res, const d2c_t x, const d2c_t y)
              void d3c_add(d3c_t res, const d3c_t x, const d3c_t y)
              void d4c_add(d4c_t res, const d4c_t x, const d4c_t y)
              void d1c_sub(d1c_t res, const d1c_t x, const d1c_t y)
              void d2c_sub(d2c_t res, const d2c_t x, const d2c_t y)
              void d3c_sub(d3c_t res, const d3c_t x, const d3c_t y)
              void d4c_sub(d4c_t res, const d4c_t x, const d4c_t y)
              void d1c_mul(d1c_t res, const d1c_t x, const d1c_t y)
              void d2c_mul(d2c_t res, const d2c_t x, const d2c_t y)
              void d3c_mul(d3c_t res, const d3c_t x, const d3c_t y)
              void d4c_mul(d4c_t res, const d4c_t x, const d4c_t y)
              void d1c_sqr(d1c_t res, const d1c_t x)
              void d2c_sqr(d2c_t res, const d2c_t x)
              void d3c_sqr(d3c_t res, const d3c_t x)
              void d4c_sqr(d4c_t res, const d4c_t x)
              void d1c_mul_re(d1c_t res, const d1c_t x, const d1_t y)
              void d2c_mul_re(d2c_t res, const d2c_t x, const d2_t y)
              void d3c_mul_re(d3c_t res, const d3c_t x, const d3_t y)
              void d4c_mul_re(d4c_t res, const d4c_t x, const d4_t y)
              void d1c_mul_onei(d1c_t res, const d1c_t x)
              void d2c_mul_onei(d2c_t res, const d2c_t x)
              void d3c_mul_onei(d3c_t res, const d3c_t x)
              void d4c_mul_onei(d4c_t res, const d4c_t x)
              void d1c_mul_2exp_si(d1c_t res, const d1c_t x, slong e)
              void d2c_mul_2exp_si(d2c_t res, const d2c_t x, slong e)
              void d3c_mul_2exp_si(d3c_t res, const d3c_t x, slong e)
              void d4c_mul_2exp_si(d4c_t res, const d4c_t x, slong e)
              void d1c_div(d1c_t res, const d1c_t x, const d1c_t y)
              void d2c_div(d2c_t res, const d2c_t x, const d2c_t y)
              void d3c_div(d3c_t res, const d3c_t x, const d3c_t y)
              void d4c_div(d4c_t res, const d4c_t x, const d4c_t y)
              void d1c_inv(d1c_t res, const d1c_t x)
              void d2c_inv(d2c_t res, const d2c_t x)
              void d3c_inv(d3c_t res, const d3c_t x)
              void d4c_inv(d4c_t res, const d4c_t x)
              void d1c_sqrt(d1c_t res, const d1c_t x)
              void d2c_sqrt(d2c_t res, const d2c_t x)
              void d3c_sqrt(d3c_t res, const d3c_t x)
              void d4c_sqrt(d4c_t res, const d4c_t x)
              void d1c_rsqrt(d1c_t res, const d1c_t x)
              void d2c_rsqrt(d2c_t res, const d2c_t x)
              void d3c_rsqrt(d3c_t res, const d3c_t x)
              void d4c_rsqrt(d4c_t res, const d4c_t x)
              void d1c_abs(d1_t res, const d1c_t x)
              void d2c_abs(d2_t res, const d2c_t x)
              void d3c_abs(d3_t res, const d3c_t x)
              void d4c_abs(d4_t res, const d4c_t x)
              void d1c_arg(d1_t res, const d1c_t x)
              void d2c_arg(d2_t res, const d2c_t x)
              void d3c_arg(d3_t res, const d3c_t x)
              void d4c_arg(d4_t res, const d4c_t x)
              void d1cb_neg(d1cb_t res, const d1cb_t x)
              void d2cb_neg(d2cb_t res, const d2cb_t x)
              void d3cb_neg(d3cb_t res, const d3cb_t x)
              void d4cb_neg(d4cb_t res, const d4cb_t x)
              void d1cb_conj(d1cb_t res, const d1cb_t x)
              void d2cb_conj(d2cb_t res, const d2cb_t x)
              void d3cb_conj(d3cb_t res, const d3cb_t x)
              void d4cb_conj(d4cb_t res, const d4cb_t x)
              void d1cb_add(d1cb_t res, const d1cb_t x, const d1cb_t y)
              void d2cb_add(d2cb_t res, const d2cb_t x, const d2cb_t y)
              void d3cb_add(d3cb_t res, const d3cb_t x, const d3cb_t y)
              void d4cb_add(d4cb_t res, const d4cb_t x, const d4cb_t y)
              void d1cb_sub(d1cb_t res, const d1cb_t x, const d1cb_t y)
              void d2cb_sub(d2cb_t res, const d2cb_t x, const d2cb_t y)
              void d3cb_sub(d3cb_t res, const d3cb_t x, const d3cb_t y)
              void d4cb_sub(d4cb_t res, const d4cb_t x, const d4cb_t y)
              void d1cb_mul(d1cb_t res, const d1cb_t x, const d1cb_t y)
              void d2cb_mul(d2cb_t res, const d2cb_t x, const d2cb_t y)
              void d3cb_mul(d3cb_t res, const d3cb_t x, const d3cb_t y)
              void d4cb_mul(d4cb_t res, const d4cb_t x, const d4cb_t y)
              void d1cb_sqr(d1cb_t res, const d1cb_t x)
              void d2cb_sqr(d2cb_t res, const d2cb_t x)
              void d3cb_sqr(d3cb_t res, const d3cb_t x)
              void d4cb_sqr(d4cb_t res, const d4cb_t x)
              void d1cb_mul_re(d1cb_t res, const d1cb_t x, const d1b_t y)
              void d2cb_mul_re(d2cb_t res, const d2cb_t x, const d2b_t y)
              void d3cb_mul_re(d3cb_t res, const d3cb_t x, const d3b_t y)
              void d4cb_mul_re(d4cb_t res, const d4cb_t x, const d4b_t y)
              void d1cb_mul_onei(d1cb_t res, const d1cb_t x)
              void d2cb_mul_onei(d2cb_t res, const d2cb_t x)
              void d3cb_mul_onei(d3cb_t res, const d3cb_t x)
              void d4cb_mul_onei(d4cb_t res, const d4cb_t x)
              void d1cb_mul_2exp_si(d1cb_t res, const d1cb_t x, slong e)
              void d2cb_mul_2exp_si(d2cb_t res, const d2cb_t x, slong e)
              void d3cb_mul_2exp_si(d3cb_t res, const d3cb_t x, slong e)
              void d4cb_mul_2exp_si(d4cb_t res, const d4cb_t x, slong e)
              int d1cb_div(d1cb_t res, const d1cb_t x, const d1cb_t y)
              int d2cb_div(d2cb_t res, const d2cb_t x, const d2cb_t y)
              int d3cb_div(d3cb_t res, const d3cb_t x, const d3cb_t y)
              int d4cb_div(d4cb_t res, const d4cb_t x, const d4cb_t y)
              int d1cb_inv(d1cb_t res, const d1cb_t x)
              int d2cb_inv(d2cb_t res, const d2cb_t x)
              int d3cb_inv(d3cb_t res, const d3cb_t x)
              int d4cb_inv(d4cb_t res, const d4cb_t x)
              void d1cb_sqrt(d1cb_t res, const d1cb_t x)
              void d2cb_sqrt(d2cb_t res, const d2cb_t x)
              void d3cb_sqrt(d3cb_t res, const d3cb_t x)
              void d4cb_sqrt(d4cb_t res, const d4cb_t x)
              int d1cb_rsqrt(d1cb_t res, const d1cb_t x)
              int d2cb_rsqrt(d2cb_t res, const d2cb_t x)
              int d3cb_rsqrt(d3cb_t res, const d3cb_t x)
              int d4cb_rsqrt(d4cb_t res, const d4cb_t x)
              void d1cb_abs(d1b_t res, const d1cb_t x)
              void d2cb_abs(d2b_t res, const d2cb_t x)
              void d3cb_abs(d3b_t res, const d3cb_t x)
              void d4cb_abs(d4b_t res, const d4cb_t x)
              void d1cb_arg(d1b_t res, const d1cb_t x)
              void d2cb_arg(d2b_t res, const d2cb_t x)
              void d3cb_arg(d3b_t res, const d3cb_t x)
              void d4cb_arg(d4b_t res, const d4cb_t x)

    Arithmetic, with the conventions and algorithms above
    (``mul_re`` multiplies by a real number, ``mul_onei`` by `i`; ``arg``
    is `\operatorname{atan2}(b, a) \in (-\pi, \pi]`, 0 for 0). The ball
    functions that can fail return a status.

.. function:: int d1c_is_zero(const d1c_t x)
              int d2c_is_zero(const d2c_t x)
              int d3c_is_zero(const d3c_t x)
              int d4c_is_zero(const d4c_t x)
              int d1c_is_one(const d1c_t x)
              int d2c_is_one(const d2c_t x)
              int d3c_is_one(const d3c_t x)
              int d4c_is_one(const d4c_t x)
              int d1c_is_real(const d1c_t x)
              int d2c_is_real(const d2c_t x)
              int d3c_is_real(const d3c_t x)
              int d4c_is_real(const d4c_t x)
              int d1c_is_finite(const d1c_t x)
              int d2c_is_finite(const d2c_t x)
              int d3c_is_finite(const d3c_t x)
              int d4c_is_finite(const d4c_t x)
              int d1c_equal(const d1c_t x, const d1c_t y)
              int d2c_equal(const d2c_t x, const d2c_t y)
              int d3c_equal(const d3c_t x, const d3c_t y)
              int d4c_equal(const d4c_t x, const d4c_t y)
              int d1cb_is_zero(const d1cb_t x)
              int d2cb_is_zero(const d2cb_t x)
              int d3cb_is_zero(const d3cb_t x)
              int d4cb_is_zero(const d4cb_t x)
              int d1cb_is_one(const d1cb_t x)
              int d2cb_is_one(const d2cb_t x)
              int d3cb_is_one(const d3cb_t x)
              int d4cb_is_one(const d4cb_t x)
              int d1cb_is_real(const d1cb_t x)
              int d2cb_is_real(const d2cb_t x)
              int d3cb_is_real(const d3cb_t x)
              int d4cb_is_real(const d4cb_t x)
              int d1cb_is_exact(const d1cb_t x)
              int d2cb_is_exact(const d2cb_t x)
              int d3cb_is_exact(const d3cb_t x)
              int d4cb_is_exact(const d4cb_t x)
              int d1cb_is_finite(const d1cb_t x)
              int d2cb_is_finite(const d2cb_t x)
              int d3cb_is_finite(const d3cb_t x)
              int d4cb_is_finite(const d4cb_t x)
              int d1cb_contains_zero(const d1cb_t x)
              int d2cb_contains_zero(const d2cb_t x)
              int d3cb_contains_zero(const d3cb_t x)
              int d4cb_contains_zero(const d4cb_t x)
              int d1cb_contains(const d1cb_t x, const d1cb_t y)
              int d2cb_contains(const d2cb_t x, const d2cb_t y)
              int d3cb_contains(const d3cb_t x, const d3cb_t y)
              int d4cb_contains(const d4cb_t x, const d4cb_t y)
              int d1cb_overlaps(const d1cb_t x, const d1cb_t y)
              int d2cb_overlaps(const d2cb_t x, const d2cb_t y)
              int d3cb_overlaps(const d3cb_t x, const d3cb_t y)
              int d4cb_overlaps(const d4cb_t x, const d4cb_t y)
              int d1cb_equal(const d1cb_t x, const d1cb_t y)
              int d2cb_equal(const d2cb_t x, const d2cb_t y)
              int d3cb_equal(const d3cb_t x, const d3cb_t y)
              int d4cb_equal(const d4cb_t x, const d4cb_t y)

    Predicates, those of the real parts combined (``is_real``: the
    imaginary part is exactly zero); for balls, answers of 1 are
    certain except for ``contains_zero`` and ``overlaps``, whose answers
    of 0 are.

.. function:: void d1c_exp(d1c_t res, const d1c_t x)
              void d2c_exp(d2c_t res, const d2c_t x)
              void d3c_exp(d3c_t res, const d3c_t x)
              void d4c_exp(d4c_t res, const d4c_t x)
              void d1c_log(d1c_t res, const d1c_t x)
              void d2c_log(d2c_t res, const d2c_t x)
              void d3c_log(d3c_t res, const d3c_t x)
              void d4c_log(d4c_t res, const d4c_t x)
              void d1c_sin(d1c_t res, const d1c_t x)
              void d2c_sin(d2c_t res, const d2c_t x)
              void d3c_sin(d3c_t res, const d3c_t x)
              void d4c_sin(d4c_t res, const d4c_t x)
              void d1c_cos(d1c_t res, const d1c_t x)
              void d2c_cos(d2c_t res, const d2c_t x)
              void d3c_cos(d3c_t res, const d3c_t x)
              void d4c_cos(d4c_t res, const d4c_t x)
              void d1c_sin_cos(d1c_t sn, d1c_t cs, const d1c_t x)
              void d2c_sin_cos(d2c_t sn, d2c_t cs, const d2c_t x)
              void d3c_sin_cos(d3c_t sn, d3c_t cs, const d3c_t x)
              void d4c_sin_cos(d4c_t sn, d4c_t cs, const d4c_t x)
              void d1c_tan(d1c_t res, const d1c_t x)
              void d2c_tan(d2c_t res, const d2c_t x)
              void d3c_tan(d3c_t res, const d3c_t x)
              void d4c_tan(d4c_t res, const d4c_t x)
              void d1c_sinh(d1c_t res, const d1c_t x)
              void d2c_sinh(d2c_t res, const d2c_t x)
              void d3c_sinh(d3c_t res, const d3c_t x)
              void d4c_sinh(d4c_t res, const d4c_t x)
              void d1c_cosh(d1c_t res, const d1c_t x)
              void d2c_cosh(d2c_t res, const d2c_t x)
              void d3c_cosh(d3c_t res, const d3c_t x)
              void d4c_cosh(d4c_t res, const d4c_t x)
              void d1c_tanh(d1c_t res, const d1c_t x)
              void d2c_tanh(d2c_t res, const d2c_t x)
              void d3c_tanh(d3c_t res, const d3c_t x)
              void d4c_tanh(d4c_t res, const d4c_t x)
              void d1c_pow(d1c_t res, const d1c_t x, const d1c_t y)
              void d2c_pow(d2c_t res, const d2c_t x, const d2c_t y)
              void d3c_pow(d3c_t res, const d3c_t x, const d3c_t y)
              void d4c_pow(d4c_t res, const d4c_t x, const d4c_t y)
              void d1cb_exp(d1cb_t res, const d1cb_t x)
              void d2cb_exp(d2cb_t res, const d2cb_t x)
              void d3cb_exp(d3cb_t res, const d3cb_t x)
              void d4cb_exp(d4cb_t res, const d4cb_t x)
              int d1cb_log(d1cb_t res, const d1cb_t x)
              int d2cb_log(d2cb_t res, const d2cb_t x)
              int d3cb_log(d3cb_t res, const d3cb_t x)
              int d4cb_log(d4cb_t res, const d4cb_t x)
              void d1cb_sin(d1cb_t res, const d1cb_t x)
              void d2cb_sin(d2cb_t res, const d2cb_t x)
              void d3cb_sin(d3cb_t res, const d3cb_t x)
              void d4cb_sin(d4cb_t res, const d4cb_t x)
              void d1cb_cos(d1cb_t res, const d1cb_t x)
              void d2cb_cos(d2cb_t res, const d2cb_t x)
              void d3cb_cos(d3cb_t res, const d3cb_t x)
              void d4cb_cos(d4cb_t res, const d4cb_t x)
              void d1cb_sin_cos(d1cb_t sn, d1cb_t cs, const d1cb_t x)
              void d2cb_sin_cos(d2cb_t sn, d2cb_t cs, const d2cb_t x)
              void d3cb_sin_cos(d3cb_t sn, d3cb_t cs, const d3cb_t x)
              void d4cb_sin_cos(d4cb_t sn, d4cb_t cs, const d4cb_t x)
              int d1cb_tan(d1cb_t res, const d1cb_t x)
              int d2cb_tan(d2cb_t res, const d2cb_t x)
              int d3cb_tan(d3cb_t res, const d3cb_t x)
              int d4cb_tan(d4cb_t res, const d4cb_t x)
              void d1cb_sinh(d1cb_t res, const d1cb_t x)
              void d2cb_sinh(d2cb_t res, const d2cb_t x)
              void d3cb_sinh(d3cb_t res, const d3cb_t x)
              void d4cb_sinh(d4cb_t res, const d4cb_t x)
              void d1cb_cosh(d1cb_t res, const d1cb_t x)
              void d2cb_cosh(d2cb_t res, const d2cb_t x)
              void d3cb_cosh(d3cb_t res, const d3cb_t x)
              void d4cb_cosh(d4cb_t res, const d4cb_t x)
              int d1cb_tanh(d1cb_t res, const d1cb_t x)
              int d2cb_tanh(d2cb_t res, const d2cb_t x)
              int d3cb_tanh(d3cb_t res, const d3cb_t x)
              int d4cb_tanh(d4cb_t res, const d4cb_t x)
              int d1cb_pow(d1cb_t res, const d1cb_t x, const d1cb_t y)
              int d2cb_pow(d2cb_t res, const d2cb_t x, const d2cb_t y)
              int d3cb_pow(d3cb_t res, const d3cb_t x, const d3cb_t y)
              int d4cb_pow(d4cb_t res, const d4cb_t x, const d4cb_t y)

    The elementary functions (principal branches), as described
    above.

.. function:: void d1c_get_acb(acb_t res, const d1c_t x)
              void d2c_get_acb(acb_t res, const d2c_t x)
              void d3c_get_acb(acb_t res, const d3c_t x)
              void d4c_get_acb(acb_t res, const d4c_t x)
              void d1c_set_acb(d1c_t res, const acb_t x)
              void d2c_set_acb(d2c_t res, const acb_t x)
              void d3c_set_acb(d3c_t res, const acb_t x)
              void d4c_set_acb(d4c_t res, const acb_t x)
              char * d1c_get_str(const d1c_t x, slong digits)
              char * d2c_get_str(const d2c_t x, slong digits)
              char * d3c_get_str(const d3c_t x, slong digits)
              char * d4c_get_str(const d4c_t x, slong digits)
              void d1c_print(const d1c_t x)
              void d2c_print(const d2c_t x)
              void d3c_print(const d3c_t x)
              void d4c_print(const d4c_t x)
              void d1c_randtest(d1c_t res, flint_rand_t state)
              void d2c_randtest(d2c_t res, flint_rand_t state)
              void d3c_randtest(d3c_t res, flint_rand_t state)
              void d4c_randtest(d4c_t res, flint_rand_t state)
              void d1c_randtest_special(d1c_t res, flint_rand_t state)
              void d2c_randtest_special(d2c_t res, flint_rand_t state)
              void d3c_randtest_special(d3c_t res, flint_rand_t state)
              void d4c_randtest_special(d4c_t res, flint_rand_t state)
              void d1cb_get_acb(acb_t res, const d1cb_t x)
              void d2cb_get_acb(acb_t res, const d2cb_t x)
              void d3cb_get_acb(acb_t res, const d3cb_t x)
              void d4cb_get_acb(acb_t res, const d4cb_t x)
              void d1cb_set_acb(d1cb_t res, const acb_t x)
              void d2cb_set_acb(d2cb_t res, const acb_t x)
              void d3cb_set_acb(d3cb_t res, const acb_t x)
              void d4cb_set_acb(d4cb_t res, const acb_t x)
              char * d1cb_get_str(const d1cb_t x, slong digits)
              char * d2cb_get_str(const d2cb_t x, slong digits)
              char * d3cb_get_str(const d3cb_t x, slong digits)
              char * d4cb_get_str(const d4cb_t x, slong digits)
              void d1cb_print(const d1cb_t x)
              void d2cb_print(const d2cb_t x)
              void d3cb_print(const d3cb_t x)
              void d4cb_print(const d4cb_t x)
              void d1cb_randtest(d1cb_t res, flint_rand_t state)
              void d2cb_randtest(d2cb_t res, flint_rand_t state)
              void d3cb_randtest(d3cb_t res, flint_rand_t state)
              void d4cb_randtest(d4cb_t res, flint_rand_t state)
              void d1cb_randtest_special(d1cb_t res, flint_rand_t state)
              void d2cb_randtest_special(d2cb_t res, flint_rand_t state)
              void d3cb_randtest_special(d3cb_t res, flint_rand_t state)
              void d4cb_randtest_special(d4cb_t res, flint_rand_t state)

    Conversions part by part (as for the real types), printing and
    random elements (with real, imaginary and zero parts among them).

.. function:: void _d1c_vec_add(d1c_ptr res, d1c_srcptr x, d1c_srcptr y, slong len)
              void _d2c_vec_add(d2c_ptr res, d2c_srcptr x, d2c_srcptr y, slong len)
              void _d3c_vec_add(d3c_ptr res, d3c_srcptr x, d3c_srcptr y, slong len)
              void _d4c_vec_add(d4c_ptr res, d4c_srcptr x, d4c_srcptr y, slong len)
              void _d1c_vec_sub(d1c_ptr res, d1c_srcptr x, d1c_srcptr y, slong len)
              void _d2c_vec_sub(d2c_ptr res, d2c_srcptr x, d2c_srcptr y, slong len)
              void _d3c_vec_sub(d3c_ptr res, d3c_srcptr x, d3c_srcptr y, slong len)
              void _d4c_vec_sub(d4c_ptr res, d4c_srcptr x, d4c_srcptr y, slong len)
              void _d1c_vec_mul(d1c_ptr res, d1c_srcptr x, d1c_srcptr y, slong len)
              void _d2c_vec_mul(d2c_ptr res, d2c_srcptr x, d2c_srcptr y, slong len)
              void _d3c_vec_mul(d3c_ptr res, d3c_srcptr x, d3c_srcptr y, slong len)
              void _d4c_vec_mul(d4c_ptr res, d4c_srcptr x, d4c_srcptr y, slong len)
              void _d1c_vec_mul_scalar(d1c_ptr res, d1c_srcptr x, slong len, d1c_srcptr c)
              void _d2c_vec_mul_scalar(d2c_ptr res, d2c_srcptr x, slong len, d2c_srcptr c)
              void _d3c_vec_mul_scalar(d3c_ptr res, d3c_srcptr x, slong len, d3c_srcptr c)
              void _d4c_vec_mul_scalar(d4c_ptr res, d4c_srcptr x, slong len, d4c_srcptr c)
              void _d1c_vec_addmul_scalar(d1c_ptr res, d1c_srcptr x, slong len, d1c_srcptr c)
              void _d2c_vec_addmul_scalar(d2c_ptr res, d2c_srcptr x, slong len, d2c_srcptr c)
              void _d3c_vec_addmul_scalar(d3c_ptr res, d3c_srcptr x, slong len, d3c_srcptr c)
              void _d4c_vec_addmul_scalar(d4c_ptr res, d4c_srcptr x, slong len, d4c_srcptr c)
              void _d1c_vec_submul_scalar(d1c_ptr res, d1c_srcptr x, slong len, d1c_srcptr c)
              void _d2c_vec_submul_scalar(d2c_ptr res, d2c_srcptr x, slong len, d2c_srcptr c)
              void _d3c_vec_submul_scalar(d3c_ptr res, d3c_srcptr x, slong len, d3c_srcptr c)
              void _d4c_vec_submul_scalar(d4c_ptr res, d4c_srcptr x, slong len, d4c_srcptr c)
              void _d1c_vec_dot(d1c_t res, d1c_srcptr initial, int subtract, d1c_srcptr x, d1c_srcptr y, slong len)
              void _d2c_vec_dot(d2c_t res, d2c_srcptr initial, int subtract, d2c_srcptr x, d2c_srcptr y, slong len)
              void _d3c_vec_dot(d3c_t res, d3c_srcptr initial, int subtract, d3c_srcptr x, d3c_srcptr y, slong len)
              void _d4c_vec_dot(d4c_t res, d4c_srcptr initial, int subtract, d4c_srcptr x, d4c_srcptr y, slong len)
              void _d1c_vec_dot_rev(d1c_t res, d1c_srcptr initial, int subtract, d1c_srcptr x, d1c_srcptr y, slong len)
              void _d2c_vec_dot_rev(d2c_t res, d2c_srcptr initial, int subtract, d2c_srcptr x, d2c_srcptr y, slong len)
              void _d3c_vec_dot_rev(d3c_t res, d3c_srcptr initial, int subtract, d3c_srcptr x, d3c_srcptr y, slong len)
              void _d4c_vec_dot_rev(d4c_t res, d4c_srcptr initial, int subtract, d4c_srcptr x, d4c_srcptr y, slong len)
              void _d1c_vec_exp(d1c_ptr res, d1c_srcptr x, slong len)
              void _d2c_vec_exp(d2c_ptr res, d2c_srcptr x, slong len)
              void _d3c_vec_exp(d3c_ptr res, d3c_srcptr x, slong len)
              void _d4c_vec_exp(d4c_ptr res, d4c_srcptr x, slong len)

    Vector operations (also for ``dNcb``), through the real vector
    operations as described above.

Application: the sieved power sum
-------------------------------------------------------------------------------

.. function:: int gr_powsum_sieved(acb_t res, const acb_t s, ulong N, gr_ctx_t T, gr_ctx_t P, gr_ctx_t F, double dm, double da, ulong M)

    Sets *res* to an enclosure of `\sum_{n=1}^N n^{-s}`, as
    :func:`acb_dirichlet_powsum_sieved` with ``len = 1``, computed
    with generic rings through runtime dispatch: the terms live in
    the ball ring *T* (a real ring; complex values are pairs of real
    vectors), the phases `t \log p` in the ball ring *P*, which should
    carry about `\log_2(t \log N)` more bits than *T*, and, if *F* is
    not ``NULL``, the composite terms and the block sums in the plain
    (non-ball) ring *F*, whose product and sum must have relative
    errors of at most *dm* and *da* (relative to `|x| |y|` and to
    `|x| + |y|`), with the error of that part bounded a priori (see
    below); with ``F = NULL`` every term is a ball. Returns a *gr*
    status (``GR_UNABLE`` if the a priori bound cannot be made
    finite, if `N \ge 2^{50}`, or if some operation fails). The
    function has been run with the ``dfloat`` rings (the fast case)
    and with ``arb`` for *T* and *P* and ``arf`` for *F*. The table
    bound *M* (see below) selects the table of the full algorithm for
    `M \ge N/3`; otherwise the two-phase version with a table bound
    of `\max(M, \lceil N^{2/3} \rceil)`; ``M = 0`` (the default) is
    the full table for `N < 2048` and `\lceil N^{2/3} \rceil`
    otherwise.

    The algorithm is that of :func:`acb_dirichlet_powsum_sieved`
    (terms at composite `n` as products `f(d) f(n/d)` with `d` the
    smallest prime factor, a table of `f(k)` for odd `k \le N/3`, the
    even `n` by Horner's rule in `2^{-s}` over the dyadic ranges of
    the odd part), rearranged so that all the bulk work goes through
    vector kernels: the odd numbers are processed in blocks of 4096
    with a segmented sieve that carries, for each small prime, its
    next multiple and the cofactor of that multiple (no divisions);
    the primes of a block are evaluated together (``_gr_vec_log`` in
    *P*, the phase reduced modulo `2 \pi` as `(t/2\pi) \log p` minus a
    nearby integer, subtracted exactly with ``_gr_vec_sub``, rounded to
    *T* with ``_gr_vec_set_other``, then ``_gr_vec_sin_cos``,
    ``_gr_vec_exp`` and ``_gr_vec_mul``; on the critical line
    `\sigma = 1/2` the magnitudes `p^{-1/2}` come from
    ``_gr_vec_rsqrt`` of the exact `p` instead of the exponential,
    which is 2-3 times cheaper and a few bits more accurate); the
    composites of a block
    are ``_gr_vec_gather`` from the table (real and imaginary parts
    interleaved, so that a lookup touches one cache line), four
    ``_gr_vec_mul`` and ``_gr_vec_add``/``_gr_vec_sub``, and
    ``_gr_vec_scatter`` back (for blocks beyond `N/3` of the first
    ones, where all cofactors precede the block; the first block is
    done with scalar operations); the terms of each dyadic range are
    summed pairwise with ``_gr_vec_add``. Only the scalar bookkeeping
    (the constants, the integers subtracted from the phases via
    ``gr_get_d``/``gr_set_d``, the Horner recurrence over the
    `\log_2 N` buckets) goes through scalar dispatch.

    **Bounded memory: the two-phase version.** The table of `f(k)`,
    `k \le N/3`, takes `8 N / 3` bytes per double of *F*, which is
    what limited the full algorithm (like the ``acb`` version, which
    switches to :func:`acb_dirichlet_powsum_smooth` beyond 4 GB). With
    a table bound `M < N/3` (and `M \ge N^{2/3}`), the odd `k \le M`
    are processed as above (with the table of `f(k)`, `k \le M`),
    and for the odd `c > M` only the *roots* are: the primes and the
    composites with cofactor `c/d \le M`, `d` the smallest prime
    factor of `c`. Every `n \le N` whose odd part exceeds `M` is
    uniquely `n = m c` with `c` a root and `m \le N/c < N/M` an
    integer whose odd part has no prime factor above `d` (the chain
    `n, n/d_1, n/(d_1 d_2), \ldots` of divisions by the smallest
    prime factor reaches `c` as the last element above `M`), so that

    .. math ::

        \sum_{n \le N} f(n) = \sum_j f(2)^j \sum_{k \text{ odd} \le \min(M, N/2^j)} f(k)
            + \sum_{c \text{ root}} f(c) \, G(\lfloor N/c \rfloor, d_c)

    where `G(z, p)` is the sum of `f(m)` over the `m \le z` whose odd
    part has no prime factor above `p` (and `d_c = \infty` for a
    prime). With `Z = \lfloor N/(M + 1) \rfloor \le N^{1/3}`, *G*
    is a table of balls in *T* (computed directly, `f(m)` for
    `m \le Z`, then cumulative sums: a full row and one row for each
    odd prime `p < Z`, about `Z^2 / (2 \log Z)` entries). The second
    phase gathers `G(\lfloor N/c \rfloor, d_c)` per root (the
    quotient by a double division with an exact correction), and one
    complex ``_gr_vec_mul`` and a pairwise sum per block; the other
    composites (all but about a fraction `1.12 / \log(N/M)` of them)
    are never computed. With a plain *F*, *G* is used through its
    midpoints and its largest radius `r_G`; as `g(k) = k^{-\sigma_{lo}}`
    is completely multiplicative, the sum over the roots of `g(c)`
    times the sum of `g(m)` over the `m` of `G` is exactly the sum of
    `g(n)` over the `n \le N` with odd part above `M`, so the a priori
    bound is as tight as in the first phase (plus `r_G` times the sum
    of `g(c)` over the odd `c > M`).

    **Threads.** The second phase, which is all but `O(N^{2/3})` of the
    work, is split into chunks of a fixed number of odd `c` (at least
    `2^{18}`, at least 16 sieve segments, at most 65536 chunks), each
    with its own scratch and sieve state (the first odd multiples of
    the small primes from the start of the chunk), run with
    :func:`flint_parallel_do` on :func:`flint_get_num_threads` threads;
    the partial sums are added in the order of the chunks, so the
    result does not depend on the number of threads. The sieve runs
    over segments of about `\sqrt{N}/2` odd numbers (up to 128 blocks
    of 4096), so that the small primes are visited about once per
    `\sqrt N` numbers even when `\pi(\sqrt N)` exceeds the block
    length. On the 2-core test machine, 2 threads take 0.095 s instead
    of 0.194 s at `N = 10^7` and 7.5 s instead of 15.4 s at `N = 10^9`
    (``(1, 3)``); `N = 10^{10}` (`t \approx 6.3 \cdot 10^{20}`) takes
    76 s on 2 threads, 15.2 ns per term and thread, in 64 MB.

    Besides bounding the memory by `O(N^{2/3})`, this is faster: the
    table fits in the cache, and fewer products are done. At
    `N = 10^8` (one core), with ``(T, P) = (1, 3)``, ``(2, 4)`` and
    ``(3, 4)`` (balls): 24.6, 38.8 and 72 ns per term with the full
    table (262 MB and 517 MB of memory for the first two), 15.8, 24.3
    and 54.5 ns with `M = N^{2/3}` (11 MB and 14 MB), with the same
    radii to within a few percent; at `N = 10^9` (`t \approx 6.3
    \cdot 10^{18}`), 15.7 and 25.3 ns per term in 20 MB and 31 MB
    (the full table would take 2.6 GB and 5.3 GB); at
    `N = 3 \cdot 10^9` (`t \approx 5.7 \cdot 10^{19}`), 48 s for
    ``(1, 3)`` in 33 MB. For comparison,
    :func:`acb_dirichlet_powsum_smooth` takes about 520 ns per term
    and :func:`acb_dirichlet_powsum_sieved` about 420 ns (at `N =
    10^7`, at the working precision of
    :func:`acb_dirichlet_zeta_rs_r`).

    With a plain ring *F*, only the prime terms are balls, and the
    error of the rest is bounded a priori and added per dyadic range:
    with `\varepsilon` the largest relative radius among the prime
    terms (obtained from the balls with
    ``_gr_vec_get_interval_mid_rad``) and
    `\mu = \sqrt 2 (\delta_m + \delta_a + \delta_m \delta_a)` the
    relative error of a plain complex product, `n` with `\Omega(n)`
    prime factors has a relative error of at most `(1 +
    \varepsilon)^{\Omega} (1 + \mu)^{\Omega - 1} - 1`, `\Omega \le
    \log_3 N`; the pairwise sums of at most 4096 terms add `((1 +
    \delta_a)^{12} - 1)` times the sum of the magnitudes, which is
    bounded by an integral.

.. function:: int dfloat_powsum_sieved(acb_t res, const acb_t s, ulong N, slong prec)
              int _dfloat_powsum_sieved(acb_t res, const acb_t s, ulong N, int T, int P, ulong M)
              int _dfloat_powsum_sieved_ball(acb_t res, const acb_t s, ulong N, int T, int P, ulong M)

    The ``dfloat`` instantiations of :func:`gr_powsum_sieved`: terms
    in ``dTb`` balls, phases in ``dPb`` balls (a strong context, so
    that the phases are canonical after the cancellation of the
    integer part and their rounding to *T* is a truncation), and, for
    `T \le 2`, the composites and block sums in plain ``dT`` numbers,
    whose outputs are always canonical, with `\delta_m = \delta_a = u`
    for ``d1`` and `8 u^2`, `3 u^2` for ``d2`` (from the kernel model
    in ``dev/dfloat_exp_bound.py``); the ``_ball`` version keeps every term a
    ball (and is the only one for `T = 3`); its radii are about 10
    times smaller, and it is 20-30% slower. These return 1 on success
    and 0 (leaving *res* undefined) for unsupported arguments. The
    first function chooses `T` so that the terms have `b = \text{prec}
    + \log_2(N)/2 + 8` bits (the error grows as `\sqrt{N} u^T`, so
    the result has about *prec* bits after the binary point) and `P`
    so that the phase `t \log N` (with `t = \operatorname{Im}(s)`)
    leaves `b + 6` bits after the binary point (not necessarily `T`
    full doubles: at the height of the `10^{15}`-th zero, 132 bits
    take ``(3, 4)``); with `P = 4`, a phase short of that by up to 8
    bits is accepted (the result then has up to about 6 bits less than
    *prec*), and returns 0 if no pair `(T, P) \in \{(1, 2), (1, 3), (2, 3),
    (2, 4), (3, 4)\}` fits, or if `|\operatorname{Re}(s)| > 2^{10}`,
    `|t| \ge 10^{60}` or `N \ge 2^{50}`, or if the tables (with the
    default table bound) would take more than
    ``DFLOAT_POWSUM_MAX_BYTES`` (`4 \cdot 10^9`; `5 \cdot 10^8` on
    32-bit systems) bytes, which happens (64-bit)
    at `N \approx 7.5 \cdot 10^{12}` (`t \approx 3.5 \cdot 10^{26}`)
    for `T = 1` and `N \approx 2.8 \cdot 10^{12}` (`t \approx 5
    \cdot 10^{25}`) for `T = 2` (`T = 3` is limited to `t \lesssim
    10^{14}` by the phase precision `P \le 4`); the per-thread scratch
    (a few MB) is not counted. *M* is the table bound of :func:`gr_powsum_sieved`
    (0 for the default). It is used by
    :func:`acb_dirichlet_zeta_rs_r` for the main sum of the
    Riemann-Siegel formula.

    Performance at the height of the `10^{15}`-th zero (`t \approx
    2.1 \cdot 10^{14}`, `N = 5\,760\,732`, one core, same run):
    :func:`acb_dirichlet_powsum_sieved` at 179 bits takes 2.4 s;
    ``(1, 2)`` 0.13 s (radius `2^{-30}`), ``(2, 3)`` 0.27 s (radius
    `2^{-79}`), ``(3, 4)`` 0.60 s (radius `2^{-133}`). These are
    for the full table; the default two-phase version (see
    :func:`gr_powsum_sieved`) is another 25-30% faster, with the same
    radii (on a second machine: 0.076 s, 0.135 s and 0.32 s, against
    0.102 s, 0.184 s and 0.43 s with the full table). The generic
    version costs about 5% over a version templated on the element
    types (which it replaced): the dispatch is per vector, and the
    remaining per-element dispatch (the integer parts of the phases,
    the maximum of the prime radii) is a few nanoseconds per prime.
    About half of the time goes to the prime terms (`\pi(N) \approx
    4 \cdot 10^5` logarithms, sines and exponentials), a quarter to
    the composites, whose cost is dominated by the cache misses of
    the table lookups `f(n/d)` for large `d`, and the rest to the
    sieve and the sums. Inside :func:`acb_dirichlet_hardy_z` this
    makes Riemann-Siegel evaluations at the precisions used by zero
    isolation 17-18 times faster at the heights of the `10^{15}`-th
    to `10^{17}`-th zeros (0.14 s against 2.5 s at `10^{15}`, 1.3 s
    against 25 s at `10^{17}`), and 8-9 times faster at twice that
    precision; above what `(3, 4)` can deliver (about 140 bits after
    the binary point at `10^{15}`) the arb version is used. Isolating
    the `10^{15}`-th zero takes 3.6 s against 56 s with the arb
    version, the `10^{16}`-th 15 s against 259 s, and computing the
    `10^{15}`-th zero to 64 bits 5.6 s against 78 s.

    The refinement of a zero to the default precision of the
    ``zeta_zeros`` example (`64 + \log_2 n` bits) evaluates *Z* at the
    target precision plus 12 bits (see
    :func:`_acb_dirichlet_refine_hardy_z_zero`), which asks the main
    sum for about 84 bits after the binary point at any height, so
    ``(2, 4)`` up to about the `10^{17}`-th zero and ``(3, 4)`` beyond.
    ``zeta_zeros -count 1 -noplatt`` on two threads of the 2-core test
    machine: 2.3 s at `n = 10^{15}`, 8.5 s at `10^{16}`, 18.5 s at
    `10^{17}`, 81 s at `10^{18}`. The rest of the
    time is the isolation (Gram points and Turing's method, about 20-35
    evaluations at 20-25 bits after the binary point, ``(1, 2)``).

    At low heights the cost is a fixed setup of about 3-8 microseconds
    (the contexts, the constants, the scratch vectors, which are sized
    to the number of odd terms up to one block, and the error bounds,
    computed with ``mag_t``) plus the terms: against
    :func:`acb_dirichlet_powsum_sieved` at the working precision of
    :func:`acb_dirichlet_zeta_rs_r`, the arb version is up to 2-3
    times faster for `N < 10`, they break even around `N = 10`-`16`
    (`t \approx 600`-`1600`), and the dfloat version is 1.5-2 times
    faster at `N \approx 40`, 2-5 times at `N \approx 100`-`400`
    (`t \approx 10^5`-`10^6`) and 2-10 times at `N \approx 4 \cdot
    10^4` (`t = 10^{10}`), for 32 to 160 bits. So
    :func:`acb_dirichlet_zeta_rs_r` uses it only for `N \ge 16`.

    **Additions to the generic interface.** Writing the power sum
    against *gr* needed the following methods, all with generic
    fallbacks (elementwise loops) in ``gr_generic`` and overrides in
    the ``dfloat`` rings (and, where applicable, ``arb``/``acb``):

    * ``GR_METHOD_VEC_SQRT``, ``VEC_RSQRT``, ``VEC_EXP``, ``VEC_LOG``,
      ``VEC_SIN``, ``VEC_COS``, ``VEC_SIN_COS`` (:func:`_gr_vec_exp`
      etc.): elementwise square roots and elementary functions, so
      that a ring with SIMD kernels can evaluate them four at a
      time.
    * ``GR_METHOD_VEC_SET_OTHER`` (:func:`_gr_vec_set_other`):
      elementwise conversion from another ring, for the phase and
      logarithm vectors from *P* to *T* and the prime terms from *T*
      to *F* (the ``dfloat`` override is a copy, a truncation with
      the tail in the radius, or a renormalisation, by the lengths
      and the canonical flag of the source).
    * ``GR_METHOD_VEC_GATHER`` and ``GR_METHOD_VEC_SCATTER``
      (:func:`_gr_vec_gather`, :func:`_gr_vec_scatter`): indexed
      data movement, ``res[i] = vec[idx[i]]`` and ``vec[idx[i]] =
      src[i]``, for the table lookups and the placement of the prime
      and composite terms in a block (the ``dfloat`` override is a
      struct copy with prefetching).
    * ``GR_METHOD_GET_INTERVAL_MID_RAD`` and its vector form
      (:func:`gr_get_interval_mid_rad`,
      :func:`_gr_vec_get_interval_mid_rad`), the inverse of
      :func:`gr_set_interval_mid_rad`: the midpoint and the radius of
      a ball as exact elements of the ring (``GR_UNABLE`` for the
      whole line), from which the a priori
      analysis gets the relative radii of the prime terms with vector
      arithmetic and one ``gr_get_d`` per prime.
    * Conversions from the ``dfloat`` rings in ``gr_set_other`` of the
      ``arb`` and ``acb`` rings, so that a generic ball result can be
      returned as an ``acb_t``.

    Two further primitives would remove the last per-element scalar
    dispatch: a vector ``nint``/``frac`` (the integer parts of the
    phases are taken with ``gr_get_d`` and ``gr_set_d``, two calls per
    element, before one exact ``_gr_vec_sub``), and a vector
    reduction for the maximum of a vector of nonnegative numbers (the
    maximum of the relative radii is one ``gr_get_d`` per element).
    Both are a few percent of the running time.


.. function:: int _dfloat_platt_smk(acb_ptr S, const fmpz * smk_points, const arb_t t0, slong A, slong B, ulong J, slong K, int P, ulong M, slong prec)
              int _dfloat_platt_smk_dd(double * S5, const fmpz * smk_points, const arb_t t0, slong A, slong B, ulong J, slong K, int P, ulong M)

    The sums over `j` of Platt's multi-evaluation of the scaled Lambda
    function (:func:`acb_dirichlet_platt_multieval`): with the buckets
    `m < N = AB` of the `j \le J` with ``smk_points[m]`` `\le j <`
    ``smk_points[m+1]`` (as computed by the arb version),
    `S_{k N + m} = \sum_j z_j b_j^k` for `0 \le k < K`, where
    `z_j = j^{-1/2} e^{-i t_0 \log(j \sqrt{\pi})}` and
    `b_j = \log(j \sqrt{\pi})/(2\pi) - m/B`, `|b_j| \le 1/(2B)`.
    The second version omits the factor
    `c = e^{-i t_0 \log \sqrt{\pi}}` and writes compact balls, five
    doubles per entry (the real and imaginary midpoints as double-words
    and a common radius), 40 bytes instead of about 100 for an ``acb_t``.
    *M* is the bound of the table of terms (see the sieving below), or
    0 for the default; a small *M* exercises the splitting of the large
    `j`. They return 1 on success, 0 if *P* is not 3 or 4 or dfloat is
    not supported (the arb version is then used).

    **Accuracy.** The rest of the multi-evaluation amplifies an error in
    the moment `k` by about `2^{24 + 8.6 k}` with the heuristic
    parameters of the zero finder (`B = 4096`), while
    `|S_k| \lesssim W 2^{-13 k}`, `W` the sum of `j^{-1/2}` over the
    bucket: for a grid good to 75 bits after the binary point (what the
    default precision of ``zeta_zeros`` needs at any height), `S_0`
    needs triple-double accuracy, `S_1` and `S_2` double-double with
    sums of logarithmic depth, and the moments beyond `k \approx 13`
    double. So the terms `f(j) = j^{-s}` and `b_j` are ``d3b`` balls
    (phases in ``d{P}b`` balls), `S_0` is summed in ``d3b`` balls, and
    the moments `1 \le k < 14` are computed in plain double-double
    arithmetic and `k \ge 14` in plain double arithmetic, over all `j` of
    a block for one `k` at a time (vectorised), with the a priori bound

    .. math ::

        \beta^k \sum R + \sum Z \bigl(k \rho \beta^{k-1} + \beta^k (e_k + s_k (1 + e_k))\bigr)

    on each part of each bucket sum, where `Z` and `R` bound the
    midpoint and the radius of `z_j` (the rounding of the triple-double
    midpoint to a double-word included), `\rho` and `\beta` the radius
    and the magnitude of `b_j`, `e_k = 8 k u^2` (``DWTimesDW1``, at most
    `7u^2` per product) and, beyond the double-double moments, a
    rounding to double and `u` per double product, and `s_k` the depth
    of the summation tree times `4u^2` (``AccurateDWPlusDW``, at most
    `3u^2/(1 - 4u)` of the result) or `u`; the sums over a bucket use
    runs of 64 terms in 8 interleaved accumulators combined pairwise,
    then a pairwise tree over the runs, so the depth is about
    `11 + \log_2(L/64)` plus the number of blocks the bucket spans;
    plus an allowance of `2^{-1060}` per operation for underflow.

    **Sieving.** `f` is completely multiplicative and `\log` additive:
    `f(j) = f(d) f(j/d)`, `\log j = \log d + \log(j/d)` with `d` the
    smallest prime factor, from a table of `f(k)`, `\log k`
    (``d3b`` balls) for `k \le M`, by default `M = \min(J/2, 5 \cdot 10^6)`
    (500 MB; 100 MB on 32-bit systems). Beyond `M`, a composite `j` whose cofactor exceeds `M`
    is split as `j = xy` with `2 \le x, y \le M` when possible, from
    its factorization over the primes up to `\sqrt{J}` (sieved along)
    and the remaining factor `R` (1 or a prime): the prime powers go
    to `x` while it stays at most `M`, the rest to `y`, and `R` to
    whichever of them it fits in; the exact identity `xy = j` is
    checked. 1, the primes and the remaining composites (with a prime
    factor above `M`) are computed directly (as the prime terms of the
    power sum); at the `10^{18}`-th zero (`M/J \approx 1/42`), about
    80% of the `j` beyond `M` come from the table (the multi-evaluation
    takes 163 s on one thread, against 197 s with all of them computed
    directly). The
    `j` are processed in dyadic ranges `[2^i, 2^{i+1})` while
    `2^i \le M` (the cofactors of a range lie in the previous ones), then
    `(2^i, J]`, each in chunks of a fixed length (at least `2^{18}`)
    computed in parallel, in waves of two chunks per thread whose
    bucket sums are added in order (so that the result does not depend
    on the number of threads).

    **Performance.** At the `10^{15}`-th zero (`J = 7.75 \cdot 10^6`,
    `K = 37`), one thread: 0.42 microseconds per `j` (0.23 for the terms,
    0.19 for the moments), against 7.7 for the arb version at 228 bits.
    With :func:`acb_dirichlet_platt_multieval` (the transforms at
    `\text{prec} - \log_2 t_0 + 16` bits, the table rows generated per
    group of `k`, the convolutions in waves), the Platt method in
    ``zeta_zeros`` on two threads takes: 6.9 s for 300 zeros from
    `10^{12}`, 5.3 s for a zero at `10^{15}`, 9.4 s at `10^{16}`, 24 s
    at `10^{17}` and 109 s at `10^{18}` in 0.87 GB, with the same digits
    as with the arb sums and as the Riemann-Siegel method.

.. function:: void dfloat_set_dfloat(double * res, int nres, const double * x, int nx)
              void dfloat_ball_set_ball(double * res, int nres, const double * x, int nx)

    Conversion between formats of different length (the ball version
    takes and returns the radius at index ``n``); rounding to a
    smaller length is exact for the balls (the dropped terms go to
    the radius). The vector conversions :func:`_gr_vec_set_other`
    between ``dfloat`` contexts copy when the source has at most as
    many components, truncate (the tail going to the radius) when the
    source context is strong, and renormalise otherwise.

Internals
-------------------------------------------------------------------------------

The arithmetic is generated from ``src/dfloat/template.inc`` for each
`N`. The kernels in ``src/dfloat/kernels.inc`` are written against an
abstract element type and number of components and instantiated for
doubles and for ``vec4d`` vectors (one element per lane), at `N`
components and at every smaller number of components (for the
tapered algorithms, inline); they compute the error-free terms
order by order and renormalize with a bottom-up VecSum. The general
renormalization (``_renorm_terms``) additionally runs a top-down pass
that packs the residuals and falls back to a branchy VecSumErrBranch
when a gap would be left. The functions
``_dN_add_tracked``, ``_dN_mul_tracked`` etc. expose the tracked
kernels to the generic (runtime `n`) code in ``scaled.c``, which
implements the slow paths by scaling, and ``conv.c`` holds the
conversions. The tables for the elementary functions are generated
by ``dev/gen_dfloat_exp_tables.py``, and their error bounds and
tapering schedules by ``dev/dfloat_exp_bound.py``, which models
every kernel operation of the functions (the terms of each order,
the residuals of the TwoSum accumulations and the dropped terms,
the digits and remainders of the long division) with rigorous
bounds on component magnitudes; the kernels are compiled with the
option of a strong (canonical) renormalization pass
(``DFLOAT_KERNEL_STRONG``) for experiments, which the contexts do not
use since the separate pass is as fast.

The complex types are generated from ``src/dfloat/ctemplate.inc``
(the arithmetic, the vectors and the table of real operations, after
``template.inc`` in the same translation unit) and ``complex.c`` (the
rest, once for all formats).

**Code size.** Everything that depends on `N` is compiled four times,
and the kernels are fully unrolled, so the object code is large
(about 0.9 MB of text for the four types and their complex versions,
most of it the scalar and vector cores of the elementary functions and
the division and square root kernels). The small operations are
kernels (``_inline``, ``_impl``) that the public functions wrap and
that the rest of the module calls, so that no call within a shared
library goes through the PLT (a call to an exported function can be
interposed, so it is neither inlined nor made directly): the trivial
ones are always inlined, the arithmetic (and the division and square
roots) only for `N = 1`, and otherwise each is one static copy called
directly, its cost dominating that of a call. The following keep the
code from growing further, at no measurable cost in speed against
fully inlined code (operation by operation): the general
renormalization with a
run-time number of terms, used only by slow paths and conversions, is
implemented once for all `N` (``scaled.c``); the ball operations take
their error mode (automatic or fast context) as a run-time argument
of one out-of-line body (the products for `N \ge 2` as well), except
where specialising the mode in a hot loop is measurably faster (the
dot products, and ``vec_sqrt`` for `N = 1`); the division kernel for
`N \ge 3` is one copy shared by the plain division and the
elementary functions, and ``_dN_sqrt_approx`` and
``_dN_rsqrt_approx`` go through the full functions rather than
another copy of their kernels; the scalar fallbacks of the vector loops, the reversed dot
product, ``inv``, and the cores of ``expm1`` (`N \ge 2`) and
``log1p`` (`N \ge 3`) are shared out of line; the public elementary
functions are never inlined into their gr wrappers; and the
comparisons and predicates share one out-of-line bound on the
difference. Relative to inlining everything, this saves 36% of the
text and 30% of the object file with debug information.

**Benchmarking.** On Skylake-family processors, whose JCC erratum
microcode keeps jumps that cross or end on a 32-byte boundary out of
the decoded-instruction cache, the speed of the short hot loops (the
dot products, the short vector operations) varies by up to 25% with
the placement of the code, between builds of identical instructions.
Comparisons of such loops are best made with the dfloat module
assembled with ``-Wa,-mbranches-within-32B-boundaries``, which removes
the effect (and makes the text about 1% larger).
