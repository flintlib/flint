---
title: gr_dft.h
header: gr_dft.h
description: Discrete Fourier transforms of any length over generic rings, with rigorous fixed-point transforms of complex ball vectors.
category: Fast Fourier transforms
order: 3
---

This module computes discrete Fourier transforms (DFTs) of any length
$n \ge 1$ over any ring of the generic rings module `gr` that contains a
suitable $n$-th root of unity: complex balls (`acb`), exact algebraic
numbers, $\Z/m\Z$ and finite fields with a user-supplied root, and two
fixed-point contexts provided by the module itself. A transform is
applied through a precomputed *plan* (`gr_dft_pre_t`), which chooses
among radix-2 and split-radix Cooley–Tukey, the four-step algorithm,
the Good–Thomas prime factor algorithm, mixed-radix Cooley–Tukey and
Bluestein's algorithm, and stores the tables of roots of unity. On top of
the cyclic transforms, the module provides DFTs on products of cyclic
groups and on groups of Dirichlet characters, and functions for complex
ball vectors that compute rigorous enclosures in fixed-point arithmetic.
These back `acb_dft` and `acb_dirichlet_dft`.

Start with the guide [](gr_dft/getting-started), which also explains the
conventions of `gr` (contexts, status codes, element pointers) used
throughout this module. The guide [](gr_dft/performance) explains how a
plan chooses its algorithm and how to tune plans and threads. The
algorithms and their conditions on the ring are described in
[](algorithms/fft/dft-algorithms), and the error analysis of the
fixed-point transforms in [](algorithms/fft/fixed-point-dft).

## Conventions {#conventions}

Let $R$ be a ring, $n \ge 1$ and $w \in R$ with $w^n = 1$. The *forward
transform* of $x = (x_0, \ldots, x_{n-1}) \in R^n$ with respect to $w$ is
$X = (X_0, \ldots, X_{n-1})$, and the *inverse transform* of $X$ is $x$:

$$
X_k = \sum_{j=0}^{n-1} x_j w^{jk},
\qquad
x_j = \frac{1}{n} \sum_{k=0}^{n-1} X_k w^{-jk}.
$$

- The *default root* is $w = e^{-2\pi i/n}$, as in `acb_dft` and in most
  numerical libraries. Over rings without it, such as $\Z/m\Z$, the root
  must be supplied.
- The root must be a *principal* $n$-th root of unity: $w^n = 1$ and
  $w^{n/p} - 1$ is not a zero divisor for every prime $p \mid n$. Over a
  field this means that $w$ has order exactly $n$. The plan constructors
  check only part of this condition, and with a root that is not
  principal the fast algorithms do not compute the transform above.
- The inverse transform divides by $n$ at the end, so it fails with
  `GR_DOMAIN` when the ring cannot divide by $n$.
- Vectors are contiguous arrays of $n$ ring elements (`gr_ptr`). Unless
  stated otherwise, output and input may be the same array.
- Functions returning `int` return a `gr` status code: `GR_SUCCESS`, or
  `GR_DOMAIN` and `GR_UNABLE` (possibly combined with `|`) as described
  for each function. Functions declared with `WARN_UNUSED_RESULT` must
  have their status checked.

## Plans {#plans}

A plan holds everything that depends only on the length, the ring and the
root: the algorithm and its decomposition into sub-plans, the table of
roots of unity, and the auxiliary tables of the composite algorithms.
[](gr_dft/plans) explains how to construct plans and what can
go wrong.

::: type
gr_dft_pre_struct
gr_dft_pre_t

A plan for transforms of one length $n$ over one ring. `gr_dft_pre_t` is
an array of one `gr_dft_pre_struct`, passed by reference. The fields are
internal; `n` (the length) and `alg` (the resolved algorithm) may be
read.
:::

::: macro
GR_DFT_ALG_AUTO
GR_DFT_ALG_NAIVE
GR_DFT_ALG_CT
GR_DFT_ALG_BAILEY
GR_DFT_ALG_SPLIT
GR_DFT_ALG_PFA
GR_DFT_ALG_MIXED
GR_DFT_ALG_BLUESTEIN

Values 0 to 7 of the `alg` argument of the plan constructors: automatic
choice, the $O(n^2)$ definition, radix-2 Cooley–Tukey, the four-step
algorithm and split-radix (these three for powers of two), the
Good–Thomas prime factor algorithm (for lengths with at least two
distinct prime factors), mixed-radix Cooley–Tukey (any $n \ge 2$) and
Bluestein's algorithm (odd $n \ge 3$).

::: details
`GR_DFT_ALG_AUTO` chooses split-radix for powers of two in complex mode
(see `gr_dft_precomp_init_karatsuba`) and radix-2 Cooley–Tukey otherwise,
the prime factor algorithm when $n$ has at least two distinct prime
factors, and mixed-radix for powers of an odd prime $p$. For $p > 19$,
mixed-radix computes the transforms of length $p$ with a Bluestein
sub-plan when the ring allows it; see [](gr_dft/performance).
:::
:::

::: macro
GR_DFT_SCRAMBLED

Flag for the plan constructors: the forward transform leaves its output in
a permuted order and the inverse transform expects its input in the same
order, which saves the reordering passes. Supported only for powers of
two with `GR_DFT_ALG_NAIVE`, `GR_DFT_ALG_CT` and `GR_DFT_ALG_BAILEY`, and
with `GR_DFT_ALG_AUTO`, which then chooses radix-2; it is silently
cleared otherwise.

::: details
The permutation is bit reversal for the naive and radix-2 algorithms and
a matrix transposition for the four-step algorithm; it can be queried
with `gr_dft_precomp_output_perm`. The flag is useful for convolutions,
where the order of the transformed values does not matter. With
`GR_DFT_ALG_AUTO` it makes a complex-mode plan use radix-2 instead of
split-radix.
:::
:::

::: function
int gr_dft_default_root(gr_ptr w, ulong n, gr_ctx_t ctx)

Sets $w$ to the default $n$-th root of unity: 1 for $n \le 1$, $-1$ for
$n = 2$, and $e^{-2\pi i/n}$ otherwise. Returns `GR_DOMAIN` or `GR_UNABLE`
if the ring cannot provide $e^{-2\pi i/n}$, for example over $\Z/m\Z$,
finite fields, real numbers and `acf`.
:::

::: function
int gr_dft_precomp_init_root(gr_dft_pre_t P, gr_srcptr w, ulong n, int alg, int flags, gr_ctx_t ctx)
int gr_dft_precomp_init(gr_dft_pre_t P, ulong n, int alg, int flags, gr_ctx_t ctx)
int gr_dft_precomp_init_karatsuba(gr_dft_pre_t P, ulong n, int alg, int flags, gr_ctx_t real_ctx, gr_ctx_t ctx)

Initialise `P` with a plan for transforms of length $n \ge 1$ over `ctx`
with the algorithm `alg` and the flags `flags`. The first function uses
the principal $n$-th root of unity $w$, the other two the default root
$e^{-2\pi i/n}$. The third also enables *complex mode* when `real_ctx` is
not `NULL`: `ctx` must then be a complex ring whose elements are stored
as a real and an imaginary part from `real_ctx`, such as `acb` over
`arb`. On failure the plan is left cleared.

::: details
The functions return `GR_DOMAIN` if $n = 0$ or the algorithm does not
support $n$, and `gr_dft_precomp_init_root` also if a check of the root
fails: $w^{n/2} \ne -1$ for even $n$, or $w^{n/p} = 1$ for an odd prime
$p \mid n$. The other two return the status of `gr_dft_default_root` if
the ring has no default root, and `GR_UNABLE` if complex mode is
requested with mismatched element sizes. Explicitly requesting
`GR_DFT_ALG_BLUESTEIN` returns `GR_UNABLE` over rings without default
roots. The plan keeps pointers to `ctx` and `real_ctx`, which must stay
valid while it is used.
:::
:::

::: function
void gr_dft_precomp_clear(gr_dft_pre_t P)

Frees the memory used by `P`, including its sub-plans. It is harmless to
call it after a failed initialisation.
:::

::: function
void gr_dft_precomp_output_perm(ulong * perm, const gr_dft_pre_t P)

Writes to `perm` (an array of $n$ entries) the permutation by which entry
$i$ of the forward output holds $X_{\mathit{perm}[i]}$; the inverse
transform expects $X_{\mathit{perm}[i]}$ in entry $i$ of its input. This
is the identity unless the plan has `GR_DFT_SCRAMBLED` set.
:::

::: function
void gr_dft_precomp_fprint(FILE * out, const gr_dft_pre_t P)
void gr_dft_precomp_print(const gr_dft_pre_t P)

Print a description of the plan to `out` or to standard output: the
algorithm and decomposition of each node, with one indented block per
sub-plan, the tables present, and the threading settings if they are not
the defaults.

::: details
The output is meant for diagnostics and may change. It works for layouts
(see `_gr_dft_precomp_init_layout`) as well as for realised plans. The
"root err" figure is the error bound of a fixed-point root table, in
units in the last place; over other rings it is a worst-case estimate
with no meaning for the actual computation, and it is marked "(est.)"
whenever the plan does not store the full table, even if it is exact.
`gr_dft_precomp_fprint` is only declared when `stdio.h` is included
before the FLINT headers.
:::
:::

## Transforms {#transforms}

::: function
int gr_dft_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_pre_t P, gr_ctx_t ctx)
int gr_dft_inverse_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_pre_t P, gr_ctx_t ctx)

Set `res` to the forward or inverse transform of `vec` with the plan `P`.
Both vectors have length $n$ and may be the same. `ctx` must be the ring
of the plan. The inverse transform includes the division by $n$, and
returns `GR_DOMAIN` if the ring cannot perform it.
:::

::: function
int gr_dft(gr_ptr res, gr_srcptr vec, ulong n, gr_ctx_t ctx)
int gr_dft_inverse(gr_ptr res, gr_srcptr vec, ulong n, gr_ctx_t ctx)

Set `res` to the forward or inverse transform of `vec` with respect to
the default root, using a temporary plan from `gr_dft_precomp_init` with
`GR_DFT_ALG_AUTO`. Return `GR_DOMAIN` or `GR_UNABLE` if the ring has no
default root of order $n$.
:::

## Threads {#threads}

The radix-2, split-radix and four-step algorithms, and the product
transforms, distribute their work over threads of the global FLINT thread
pool when `flint_get_num_threads` is greater than 1 and
`gr_ctx_is_threadsafe` returns `T_TRUE` for the ring (and for the real
ring in complex mode). The guide [](gr_dft/performance) explains how the
work is divided.

::: function
void gr_dft_precomp_set_threads(gr_dft_pre_t P, thread_pool_handle * threads, slong num_threads)

Attaches `num_threads` worker threads to `P` and its sub-plans; every
later transform uses them instead of requesting threads from the pool.
The handles are borrowed: the caller obtains them, for example with
`flint_request_threads`, and gives them back with
`flint_give_back_threads` after the plan is cleared. Passing `NULL` and 0
detaches them.
:::

::: function
void gr_dft_precomp_set_serial_block(gr_dft_pre_t P, slong serial_block)

Sets the threading granularity of `P` and its sub-plans: each pass or
phase is divided into at most $n / \mathit{serial\_block}$ chunks, so a
transform with $n < 2 \cdot \mathit{serial\_block}$ runs serially. Zero
restores the default `GR_DFT_SERIAL_BLOCK_DEFAULT`.
:::

::: macro
GR_DFT_SERIAL_BLOCK_DEFAULT

The default threading granularity, 256 elements.
:::

::: macro
GR_DFT_CT_BLOCK_BYTES

The block size, 128 KiB, below which the radix-2 transform completes all
remaining levels of a block before moving to the next one. It is a
compile-time constant, to be at most about half the per-core cache.
:::

## Products of cyclic groups {#prod}

These functions compute the DFT on the group
$\Z/c_0\Z \times \cdots \times \Z/c_{r-1}\Z$, that is, a multidimensional
DFT of an array of $n = c_0 \cdots c_{r-1}$ entries stored in row-major
order (the first index varies slowest):

$$
X_{k_0, \ldots, k_{r-1}} = \sum_{j_0, \ldots, j_{r-1}} x_{j_0, \ldots, j_{r-1}}
\prod_{a=0}^{r-1} w_a^{j_a k_a},
$$

where $w_a$ is a principal $c_a$-th root of unity. It is computed as one
cyclic transform along each axis. [](gr_dft/prod) has an
example.

::: type
gr_dft_prod_pre_struct
gr_dft_prod_pre_t

A plan for product DFTs with fixed component lengths over one ring. It
holds one cyclic plan for each distinct component length.
:::

::: function
int gr_dft_prod_precomp_init(gr_dft_prod_pre_t P, const ulong * cyc, slong num, int flags, gr_ctx_t ctx)
int gr_dft_prod_precomp_init_root(gr_dft_prod_pre_t P, gr_srcptr w, ulong order, const ulong * cyc, slong num, int flags, gr_ctx_t ctx)

Initialise `P` for the product DFT with the `num` component lengths
`cyc`. The first function uses the default roots $w_a = e^{-2\pi i/c_a}$.
The second takes a root $w$ of the given order, which every component
length must divide, and uses $w_a = w^{\mathit{order}/c_a}$. The flags
are passed to the component plans. After any return value, `P` must be
cleared with `gr_dft_prod_precomp_clear`.

::: details
Each component plan is made with `GR_DFT_ALG_AUTO` by
`gr_dft_precomp_init` or `gr_dft_precomp_init_root`, and their status is
returned, so the conditions of those functions apply to each component.
`gr_dft_prod_precomp_init_root` returns `GR_DOMAIN` if some $c_a$ is zero
or does not divide `order`. With `num` = 0 the transform is the identity
on one entry.
:::
:::

::: function
void gr_dft_prod_precomp_clear(gr_dft_prod_pre_t P)

Frees the memory used by `P`.
:::

::: function
int gr_dft_prod_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_prod_pre_t P, gr_ctx_t ctx)
int gr_dft_prod_inverse_precomp(gr_ptr res, gr_srcptr vec, const gr_dft_prod_pre_t P, gr_ctx_t ctx)
int gr_dft_prod(gr_ptr res, gr_srcptr vec, const ulong * cyc, slong num, gr_ctx_t ctx)
int gr_dft_prod_inverse(gr_ptr res, gr_srcptr vec, const ulong * cyc, slong num, gr_ctx_t ctx)

Set `res` to the forward or inverse product DFT of `vec`, both of length
$n = c_0 \cdots c_{r-1}$, which may be the same. The last two use a
temporary plan with default roots. The inverse divides by $n$ once at the
end.
:::

::: function
void gr_dft_prod_precomp_set_serial_block(gr_dft_prod_pre_t P, slong serial_block)

Sets the threading granularity of the product transform: the lines of
one axis are distributed over at most $n / \mathit{serial\_block}$ chunks.
It does not change the component plans. Zero restores the default
`GR_DFT_SERIAL_BLOCK_DEFAULT`.
:::

## Dirichlet characters {#dirichlet}

The characters modulo $q$ form a product of cyclic groups, so the DFT on
this group is a product DFT over the orders $\varphi_k$ of the Conrey
generators (`G->P[k].phi.n` for a `dirichlet_group_t` `G`). For the
character $\chi_a$ with Conrey number $a$, the transform is
$w_a = \sum_x \overline{\chi_a(x)} v_x$, summed over $x$ coprime to $q$.

::: function
int gr_dft_dirichlet_index(gr_ptr w, gr_srcptr v, const dirichlet_group_t G, gr_ctx_t ctx)
int gr_dft_dirichlet(gr_ptr w, gr_srcptr v, const dirichlet_group_t G, gr_ctx_t ctx)

Set `w` to the DFT of `v` over the characters of `G`. In
`gr_dft_dirichlet_index` the vectors have `G->phi_q` entries in the
lexicographic order of Conrey logarithms (that of `dirichlet_char_next`).
In `gr_dft_dirichlet` they have `G->q` entries indexed by Conrey numbers;
entries at indices not coprime to $q$ are ignored in `v` and not written
in `w`. The vectors may be the same. Both use default roots, so they
need a ring with $e^{-2\pi i/\varphi_k}$ such as `acb`, and return
`GR_UNABLE` or `GR_DOMAIN` otherwise.
:::

## Complex ball vectors {#acb}

These functions compute transforms of `acb` vectors with the default
root $e^{-2\pi i/n}$. They use the fixed-point contexts internally with a
rigorous error bound, and ball arithmetic when the fixed-point path does
not apply. The outputs contain the exact transform of every vector of
complex numbers contained in the input balls. `acb_dft` and
`acb_dirichlet_dft` are thin wrappers around them. [](gr_dft/acb) explains how the method works and when it falls back.

::: type
gr_dft_acb_pre_struct
gr_dft_acb_pre_t
gr_dft_acb_prod_pre_struct
gr_dft_acb_prod_pre_t

Precomputed data for repeated transforms of `acb` vectors of one length,
or one product shape, at one working precision. `acb_dft_pre_t` is the
same type as `gr_dft_acb_pre_t`.
:::

::: function
void gr_dft_acb(acb_ptr w, acb_srcptr v, slong n, slong prec)
void gr_dft_acb_inverse(acb_ptr w, acb_srcptr v, slong n, slong prec)
int _gr_dft_acb(acb_ptr w, acb_srcptr v, slong n, int inverse, int which, slong prec)

Set `w` to the forward or inverse transform (including the division by
$n$) of `v`, both of length $n$, with working precision `prec`. The
vectors may be the same, and nothing is done for $n \le 0$.
`_gr_dft_acb` forces the method with `which` (0 automatic, 1 ball
arithmetic, 2 fixed point) and returns `GR_UNABLE` if the forced
fixed-point path cannot handle the input, in which case `w` is set to
indeterminate values.

::: details
Outputs at the scale of the largest input carry about `prec` accurate
bits. For $n < 0$, `_gr_dft_acb` returns `GR_DOMAIN`.
:::
:::

::: function
int gr_dft_acb_precomp_init(gr_dft_acb_pre_t Q, slong n, slong prec)
void gr_dft_acb_precomp_clear(gr_dft_acb_pre_t Q)
void gr_dft_acb_precomp(acb_ptr w, acb_srcptr v, const gr_dft_acb_pre_t Q, slong prec)
void gr_dft_acb_inverse_precomp(acb_ptr w, acb_srcptr v, const gr_dft_acb_pre_t Q, slong prec)
int _gr_dft_acb_precomp(acb_ptr w, acb_srcptr v, int inverse, const gr_dft_acb_pre_t Q, slong prec)

Precomputed versions of `gr_dft_acb` and `gr_dft_acb_inverse` for
repeated transforms of length $n \ge 1$. The precision given to
`gr_dft_acb_precomp_init` determines the accuracy; the `prec` of the
transforms only sets the rounding of the outputs. The initialisation
returns `GR_DOMAIN` for $n \le 0$; `Q` may be cleared after any return
status.

::: details
The initialisation chooses the fixed-point precision from the error bound
of the plan and builds the plan over the fixed-point contexts, or over
`acb` if that fails. A transform whose input the fixed-point path cannot
handle (non-finite or extremely large midpoints) falls back to a one-shot
ball transform. `_gr_dft_acb_precomp` returns the status of the
transform, which is `GR_SUCCESS` except in degenerate cases, where `w` is
set to indeterminate values.
:::
:::

::: function
void gr_dft_acb_prod(acb_ptr w, acb_srcptr v, const ulong * cyc, slong num, slong prec)
void gr_dft_acb_prod_inverse(acb_ptr w, acb_srcptr v, const ulong * cyc, slong num, slong prec)
int _gr_dft_acb_prod(acb_ptr w, acb_srcptr v, const ulong * cyc, slong num, int inverse, int which, slong prec)

The product DFT of [](#prod) for `acb` vectors, computed like
`gr_dft_acb`. `_gr_dft_acb_prod` returns `GR_DOMAIN` if a component length
is zero or the total length is too large.
:::

::: function
int gr_dft_acb_prod_precomp_init(gr_dft_acb_prod_pre_t Q, const ulong * cyc, slong num, slong prec)
void gr_dft_acb_prod_precomp_clear(gr_dft_acb_prod_pre_t Q)
void gr_dft_acb_prod_precomp(acb_ptr w, acb_srcptr v, const gr_dft_acb_prod_pre_t Q, slong prec)
void gr_dft_acb_prod_inverse_precomp(acb_ptr w, acb_srcptr v, const gr_dft_acb_prod_pre_t Q, slong prec)
int _gr_dft_acb_prod_precomp(acb_ptr w, acb_srcptr v, int inverse, const gr_dft_acb_prod_pre_t Q, slong prec)

Precomputed versions of `gr_dft_acb_prod` and `gr_dft_acb_prod_inverse`.
The initialisation returns `GR_DOMAIN` if a component length is zero or
the total length exceeds about `WORD_MAX`/8; `Q` may be cleared after any
return status.
:::

::: function
void gr_dft_acb_dirichlet_index(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)
void gr_dft_acb_dirichlet(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)

The transforms of `gr_dft_dirichlet_index` and `gr_dft_dirichlet` for
`acb` vectors, computed with `gr_dft_acb_prod`. They implement
`acb_dirichlet_dft_index` and `acb_dirichlet_dft`
([](acb_dirichlet#dft)).
:::

## Fixed-point contexts {#nfixed}

Two `gr` contexts for numerical transforms. An element of the real
context with $k$ limbs is a sign limb (0 or 1) followed by $k$ limbs
holding an absolute value in $[0, 1)$ with unit $2^{-k \bits}$ (one
*ulp*); an element of the complex context is a real and an imaginary
part. Additions are exact, and multiplications truncate toward zero.
Every value must stay below 1 in absolute value: an operation that
overflows returns `GR_UNABLE`. [](gr_dft/nfixed) shows how to
scale the input with the help of `gr_dft_precomp_nfixed_bound`.

::: function
int gr_dft_ctx_init_nfixed(gr_ctx_t ctx, slong nlimbs)
int gr_dft_ctx_init_nfixed_complex(gr_ctx_t ctx, slong nlimbs)

Initialise `ctx` as the real or complex fixed-point context with `nlimbs`
fraction limbs. Return `GR_UNABLE` unless
$1 \le \mathit{nlimbs} \le$ `GR_DFT_NFIXED_MAX_NLIMBS`. The contexts are
freed with `gr_ctx_clear`.

::: details
The contexts are not rings in the strict sense: 1 and $-1$ are rounded to
$\pm(1 - \mathrm{ulp})$, `gr_is_one` and `gr_is_neg_one` always return
`T_FALSE`, and `gr_set_si` and `gr_set_ui` return `GR_DOMAIN` for integers
other than 0 and $\pm 1$. Besides the basic operations, they support
`gr_div_ui` (truncating) and are thread-safe.
:::
:::

::: macro
GR_DFT_NFIXED_MAX_NLIMBS

The largest number of limbs of a fixed-point context, 2048.
:::

::: function
int gr_dft_nfixed_set_arf(gr_ptr res, const arf_t x, gr_ctx_t ctx)
void gr_dft_nfixed_get_arf(arf_t res, gr_srcptr x, gr_ctx_t ctx)

Convert between `arf` values and elements of a real fixed-point context.
Setting truncates toward zero, with error below 1 ulp, and saturates
values with $|x| \ge 1$ to $\pm(1 - \mathrm{ulp})$; it returns
`GR_DOMAIN` only for non-finite $x$. Getting is exact. The parts of a
complex element are converted by treating them as two consecutive real
elements.
:::

::: function
void gr_dft_precomp_nfixed_bound(double * peak, double * err_ulps, double in_mag, double in_err_ulps, const gr_dft_pre_t P)
void gr_dft_prod_precomp_nfixed_bound(double * peak, double * err_ulps, double in_mag, double in_err, const gr_dft_prod_pre_t P)

For a forward transform or an inverse transform without the division by
$n$, executed in fixed-point arithmetic with the plan `P`, given inputs
of complex modulus at most `in_mag` with errors of at most `in_err_ulps`
ulp in each part: set `peak` to a bound for the modulus of every input,
intermediate and output value, and `err_ulps` to a bound for the error of
each part of the output, in ulp. There is no overflow if `peak` is less
than 1.

::: details
The bounds assume a plan in complex mode, made by
`gr_dft_precomp_init_karatsuba` with the fixed-point real and complex
contexts. For a plan over another ring, or a layout, they describe a
hypothetical fixed-point execution with the largest multiplication error
of the model. Both bounds are linear in `in_mag` and `in_err_ulps`, and
the cost is $O(\log n)$. [](algorithms/fft/fixed-point-dft) derives
them.
:::
:::

## Two-phase plan construction {#layout}

The `acb` functions construct their plans in two steps: a *layout*, which
fixes the algorithm and the complete decomposition without any ring, so
that `gr_dft_precomp_nfixed_bound` can be evaluated before the precision
is chosen, and a *realisation* over the chosen contexts, which computes
the tables with the default roots.

::: function
int _gr_dft_precomp_init_layout(gr_dft_pre_t P, ulong n, int alg, int flags, int complex_mode)
int _gr_dft_precomp_realize(gr_dft_pre_t P, gr_ctx_struct * real_ctx, gr_ctx_t ctx)

Initialise `P` as a layout for length $n$ with the given algorithm and
flags (`complex_mode` tells the automatic choice whether the plan will be
realised in complex mode), and realise a layout over `ctx` with default
roots, in complex mode if `real_ctx` is not `NULL`. A layout can be
passed to `gr_dft_precomp_nfixed_bound`, `gr_dft_precomp_print` and
`gr_dft_precomp_clear`. On failure, both functions leave `P` cleared.
:::

::: function
int _gr_dft_prod_precomp_init_layout(gr_dft_prod_pre_t P, const ulong * cyc, slong num, int flags, int complex_mode)
int _gr_dft_prod_precomp_realize(gr_dft_prod_pre_t P, gr_ctx_struct * real_ctx, gr_ctx_t ctx)

The same for product plans. `P` must be cleared with
`gr_dft_prod_precomp_clear` after any return value of either function.
:::

## Internal functions {#internal}

These functions are used by the transforms and the tests. They take the
plan's root table directly and have restrictions that the public
functions handle.

::: function
int _gr_dft_precomp_raw(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)
int _gr_dft_prod_precomp_raw(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_prod_pre_t P, gr_ctx_t ctx)

Apply the plan, forward or inverse, without the division by $n$ of the
inverse transform. `res` and `vec` may be the same.
:::

::: function
int _gr_dft_naive(gr_ptr res, gr_srcptr vec, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx)
int _gr_dft_ct(gr_ptr x, slong stride, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx)
int _gr_dft_ct_threaded(gr_ptr x, slong stride, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx)
int _gr_dft_bailey(gr_ptr x, int inverse, int scrambled, const gr_dft_pre_t P, gr_ctx_t ctx)
int _gr_dft_split(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)
int _gr_dft_split_threaded(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)
int _gr_dft_mixed(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)
int _gr_dft_pfa(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)
int _gr_dft_bluestein(gr_ptr res, gr_srcptr vec, int inverse, const gr_dft_pre_t P, gr_ctx_t ctx)

The kernels of the algorithms, each for a plan of the matching algorithm,
without the division by $n$. The radix-2 and four-step kernels work in
place on `x` (the radix-2 kernels on the entries `x[i * stride]`); the
others write to `res`, which must not alias `vec`. The `_threaded`
versions use threads and may only be called on a top-level plan; the
others never use threads.

::: details
`scrambled` selects the bit-reversed order (the transposed order for
`_gr_dft_bailey`) and must be 0 unless $n$ is a power of two.
`_gr_dft_bluestein` computes the inverse
transform as the forward transform of the input $x_{(n - j) \bmod n}$.
:::
:::

::: function
int _gr_dft_mul_root(gr_ptr res, gr_srcptr x, ulong e, int inverse, gr_ptr rtmp, const gr_dft_pre_t P)
int _gr_dft_mul_const(gr_ptr res, gr_srcptr x, gr_srcptr w, gr_srcptr wtab3, gr_ptr rtmp, const gr_dft_pre_t P)

The first sets `res` to $w^e x$ (or $w^{-e} x$ if `inverse` is set),
for $0 \le e < n$ and a plan that stores the full root table. Radix-2
plans, and split-radix plans in complex mode, of length at least 4 and
without Karatsuba tables store per-stage tables instead. The second sets `res` to $w x$ for a plan constant $w$, using
the complex Karatsuba entry `wtab3` (three real elements $c$, $d - c$,
$d + c$ for $w = c + di$) if it is not `NULL`.
`res` must not alias `x`. In complex mode, `rtmp` is an initialised
element of the real ring used as scratch space.

::: details
In complex mode, `_gr_dft_mul_root` multiplies by the powers $\pm 1$,
$\pm i$ at no cost and by the odd multiples of $w^{n/8}$, which have the
form $c(1 \pm i)$, with two real multiplications. Other roots take three
real multiplications with the Karatsuba tables or one complex
multiplication without them.
:::
:::

::: macro
GR_DFT_ROOT_GENERIC
GR_DFT_ROOT_DMC_ZERO
GR_DFT_ROOT_DPC_ZERO
GR_DFT_ROOT_NEG_ONE
GR_DFT_ROOT_I
GR_DFT_ROOT_NEG_I

The cost classes of the roots $c + di$ in complex mode: generic, $d = c$,
$d = -c$, $-1$, $i$ and $-i$.
:::

::: function
void _gr_dft_bit_reverse(gr_ptr x, slong stride, int depth, gr_ctx_t ctx)

Permutes the entries `x[i * stride]`, $0 \le i < 2^{\mathit{depth}}$, by
bit reversal of the index, in place.
:::

::: function
int _gr_dft_ctx_is_nfixed(gr_ctx_t ctx)
int _gr_dft_ctx_is_nfixed_complex(gr_ctx_t ctx)

Return 1 if `ctx` is a real or complex fixed-point context, and 0
otherwise.
:::

::: function
int _gr_dft_nfixed_roots(gr_ptr roots, ulong n, double * err_ulps, gr_ctx_t ctx)
double _gr_dft_nfixed_root_err_bound(ulong n, double kmul)
double _gr_dft_nfixed_cmul_err_ulps(gr_ctx_t ctx)

The root table of a fixed-point plan and its error bounds. The first
function sets `roots` to approximations of $e^{-2\pi i j/n}$,
$0 \le j < n$, over the complex fixed-point context `ctx` and writes a
bound for their errors, in ulp, to `err_ulps`. The second bounds the
error of the entries of a table of length $n$ built with products whose
rounding error is at most `kmul` ulp in modulus. The third returns the
error bound per part of one complex multiplication in `ctx`: 2, 4 or 6
ulp.
:::

::: function
int _gr_dft_nfixed_cmul_schoolbook(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
int _gr_dft_nfixed_cmul_karatsuba(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)

Multiply two elements of a complex fixed-point context with four or with
three real multiplications, regardless of the limb count. They exist for
tuning the crossover between the two, 22 limbs, with the profile program
`p-gr_dft_nfixed`.
:::
