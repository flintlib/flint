---
title: fft_small.h
header: fft_small.h
description: Number-theoretic transforms modulo 50-bit primes computed in SIMD floating point, and the multimodular multiplication of integers and polynomials built on them.
category: Fast Fourier transforms
order: 2
---

This module multiplies large integers, and polynomials over $\Z/n\Z$ and
$\Z$, with number-theoretic transforms modulo a few 50-bit primes $p$ with
$p - 1$ divisible by a large power of two. The residues are stored as
integers in `double` variables and all modular arithmetic is done in
floating point with fused multiply-add, on vectors of four or eight
doubles.
Products whose coefficients exceed one prime are computed modulo up to
eight primes and reconstructed with the Chinese remainder theorem.

The module is compiled only on 64-bit targets with AVX2 (x86-64) or NEON
(ARM64), where `FLINT_HAVE_FFT_SMALL` is 1; on other targets the header
declares only `FLINT_FFT_SMALL_ALIGNMENT` and
`flint_fft_small_max_transformed_ring_size`, and FLINT multiplies large
integers with GMP's `mpn_mul` and, from 32000 limbs, the
[Schönhage–Strassen code](fft) instead.

Most programs use the module indirectly. `flint_mpn_mul` calls it for
operands of a few hundred limbs and more, and `nmod_poly_mul`,
`fmpz_poly_mul` and their truncated and middle-product variants call it
for long polynomials; the guide [](fft_small/usage) lists every caller and
its thresholds. The functions below are for code that wants to bypass
those dispatchers, reuse a transformed operand (for polynomials, or for
integers that are multiplied by a common factor), accumulate several
products in transform space ([](fft_small/transform-space)), or compute
middle products as transposed products.

::: summary algorithms/fft/multimodular :::

::: summary algorithms/fft/small-prime-fft :::

## Conventions {#conventions}

- $B = 2^{\bits}$, and a *limb* is one `ulong`. Integers are given as
  limb arrays, least significant limb first, as in GMP's `mpn` layer.
- Polynomials are coefficient arrays, constant term first. A *window*
  $[z_l, z_h)$ of a product $c = ab$ of lengths $a_n$ and $b_n$ means the
  coefficients $c_{z_l}, \ldots, c_{z_h - 1}$, where coefficients of
  index $a_n + b_n - 1$ and beyond are zero. Windows give full, low
  (`mullow`) and middle products with one function.
- Over $\Z/n\Z$, coefficients must be reduced modulo $n$.
- Functions with an `mpn_ctx_t` argument take their primes, tables and
  scratch memory from it. The functions with `default` in their name use
  the context returned by `get_default_mpn_ctx`, which is local to the
  calling thread. A context must not be used by two threads at the same
  time; the functions themselves spread their work over the threads of
  FLINT's thread pool (`flint_set_num_threads`).
- Unless stated otherwise, outputs must not overlap inputs.

## Contexts {#contexts}

::: type
mpn_ctx_struct
mpn_ctx_t

A multiplication context. It holds the transform contexts of
`MPN_CTX_NCRTS` = 8 primes, the Chinese remaindering data for the first
$k$ of them for every $k \le 8$, tables of powers of two modulo each
prime for converting integers, the table of packing profiles, the
constants of the two-prime pipeline, and a scratch buffer that grows on
demand and is kept between calls.

::: details
The primes are generated from the initial prime given to `mpn_ctx_init`.
FLINT always uses $p_0 = 63 \cdot 2^{44} + 1$, which yields

| $i$ | $p_i$ | form |
|---|---|---|
| 0 | 1108307720798209 | $63 \cdot 2^{44} + 1$ |
| 1 | 659706976665601 | $75 \cdot 2^{43} + 1$ |
| 2 | 1086317488242689 | $247 \cdot 2^{42} + 1$ |
| 3 | 910395627798529 | $207 \cdot 2^{42} + 1$ |
| 4 | 699289395265537 | $159 \cdot 2^{42} + 1$ |
| 5 | 1022545813831681 | $465 \cdot 2^{41} + 1$ |
| 6 | 1013749720809473 | $461 \cdot 2^{41} + 1$ |
| 7 | 868614185943041 | $395 \cdot 2^{41} + 1$ |

All have 50 bits. A transform of length $2^d$ modulo the first $k$
primes needs $2^d \mid p_i - 1$ for each of them, so $d \le 41$ in
general and $d \le 43$ for the two-prime pipeline: far beyond any
realistic size. The fields are described with the low-level functions
in [](#crt).
:::
:::

::: function
void mpn_ctx_init(mpn_ctx_t R, ulong p)
void mpn_ctx_clear(mpn_ctx_t R)

Initialises `R` from the initial prime `p` and frees it. FLINT itself
always uses $p = 63 \cdot 2^{44} + 1$ (`0x0003f00000000001`), the value
`get_default_mpn_ctx` passes.

::: details
The eight primes are $p_0 = p$ (or, if $p$ is not prime, the first prime
in the sequence below) and then the next primes of the sequence obtained
by repeatedly subtracting $2^{v + 1}$ from the previous value, where
$2^v \| p - 1$; when that would lose a bit, the bit size is kept and
the two-adic valuation is decreased by one. Each prime must pass
`fft_small_mulmod_satisfies_bounds` and have $2^4 \mid p - 1$, or
`sd_fft_ctx_init_prime` throws. Other starting values are meant only for
testing: the packing profiles, the two-prime pipeline and the shortcut
for moduli that are context primes all assume 50-bit primes.

Initialisation computes the twiddle factors for transforms of length up
to $2^{12}$, the cofactors $P/p_i$ for the Chinese remainder theorem, and
the tables $2^j \bmod p_i$ for $j < 256$. Longer twiddle tables are added
on first use.
:::
:::

::: function
mpn_ctx_struct * get_default_mpn_ctx(void)

Returns the default context of the calling thread, creating it on the
first call. `flint_cleanup` and `flint_cleanup_master` free it.
:::

::: function
void * mpn_ctx_fit_buffer(mpn_ctx_t R, ulong n)

Returns at least `n` bytes of scratch memory from the context, aligned to
4096 bytes. The memory belongs to the context and must not be freed. A
later call may move the buffer, so the pointer is valid only until the
next call on `R`, and only one request can be outstanding.

::: details
The buffer grows geometrically (by at least $1/16$) and is kept until
`mpn_ctx_clear`, so repeated multiplications on the same thread allocate
nothing. While a reservation made with `mpn_ctx_fit_buffer_reserve` is
live, requests of at most the reserved head size are served from the
buffer, and larger ones from a second buffer in the context with the same
rules.
:::
:::

::: function
void * mpn_ctx_fit_buffer_reserve(mpn_ctx_t R, ulong head, ulong tail)
void mpn_ctx_fit_buffer_release(mpn_ctx_t R)

Reserves the tail of the scratch buffer: the buffer is grown to hold
`head` plus `tail` bytes (each rounded up to 4096), and the returned
pointer to the last `tail` bytes stays valid until the release. Returns
`NULL` if a reservation is already live. Reserving, like any growth,
invalidates a pointer returned by `mpn_ctx_fit_buffer`.

::: details
A reservation is meant for a caller that keeps transformed operands
alive across many operations on one thread, such as an operand pool
(`fft_small_op_pool_fit` with `FFT_SMALL_OP_POOL_RESERVE`) and through it
the transformed integer rings of `gr`. While it is live, multiplications
on the same thread still work: requests that do not fit in the head are
served from the secondary buffer, at the cost of an allocation.
:::
:::

## Integer multiplication {#mpn}

The integers are cut into chunks of a fixed number of bits, the chunk
sequences are multiplied as polynomials modulo two to eight primes, and
the product is reassembled with carries. Up to 40000 limbs, two primes and
chunks of 40 to 45 bits are used; above, four to eight primes and chunks
of 84 to 192 bits, chosen from a table of profiles.

::: function
void mpn_ctx_mpn_mul(mpn_ctx_t R, ulong * z, const ulong * a, ulong an, const ulong * b, ulong bn)
void mpn_mul_default_mpn_ctx(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn)

Sets $(z, a_n + b_n)$ to the product of $(a, a_n)$ and $(b, b_n)$, using
`R` or the default context. Requires $a_n, b_n \ge 1$, in either order.
If `a == b` and $a_n = b_n$, a square is computed.

::: details
Operands with $\max(a_n, b_n) \le 40000$ (`FLINT_MPN_MUL_NP2_MAX_AN`) use
the single-threaded two-prime pipeline. Larger ones use a profile chosen
for the sizes and the number of threads available, and up to eight
threads. The number of products accumulating in one coefficient is
bounded using $b_n$, so for very unbalanced operands it is slightly better
to pass the longer one first. See [](fft_small/mpn_mul).
:::
:::

::: function
void _mpn_ctx_mpn_mul_range(mpn_ctx_t R, ulong * z, ulong lo, ulong hi, const ulong * a, ulong an, const ulong * b, ulong bn)

Writes to $(z, \mathit{hi} - \mathit{lo})$ the limbs $[\mathit{lo},
\mathit{hi})$ of the product of $(a, a_n)$ and $(b, b_n)$. The result is
exact for $\mathit{lo} = 0$; otherwise carries from below limb
$\mathit{lo}$ may be lost, making it smaller than the exact window by a
few units. Limbs at or above $a_n + b_n$ are zero.

::: details
Requires $a_n, b_n \ge 1$; nothing is written if
$\mathit{lo} \ge \mathit{hi}$. Precisely, the result is congruent modulo
$B^{\mathit{hi}-\mathit{lo}}$ to $\lfloor ab / B^{\mathit{lo}} \rfloor -
\delta$ with $0 \le \delta < \max(2, \min(a_n, b_n, \mathit{lo}))$.
Only the coefficients of the chunk product that reach the window are
reconstructed. A window shorter than the longer operand is computed as a
transposed product of the shorter operand and the window when that
transform is strictly shorter (not for squarings), with the same result;
otherwise, above 40000 limbs, the transform is cyclic of half the length
when the window allows it. See [](fft_small/mpn_mul#algorithm).
:::
:::

## Integer products with a common operand {#mpn-shared}

When several products share an operand, its transform can be computed
once. `_fft_small_mpn_mul_shared` forms two products $ax$ and $bx$, and the
`_fft_small_mpn_mul_precomp` functions keep the transform of a fixed
factor $x$ for any number of products. `flint_mpn_mul_shared` chooses
between the first and two ordinary products. See [](fft_small/mpn_shared).

::: function
void _fft_small_mpn_mul_shared(nn_ptr z1, nn_srcptr a, ulong an, nn_ptr z2, nn_srcptr b, ulong bn, nn_srcptr x, ulong xn)

Sets $(z_1, a_n + x_n)$ to $ax$ and $(z_2, b_n + x_n)$ to $bx$, computing
the transform of $x$ once: three forward and two inverse transforms
instead of four and two. Requires $a_n, b_n, x_n \ge 1$, in any order of
size. `z1` may overlap `a` and `x` but not `b`; `z2` may overlap any input
but not `z1`.

::: details
The transform of $x$ is made for the longer of the two products; each
partner is converted, transformed, multiplied pointwise and transformed
back at the truncation of its own product, so a shorter partner costs
less. Up to `FLINT_MPN_MUL_NP2_MAX_AN` limbs (for the longest operand)
the two-prime geometry is used and each product runs in one pass over its
data, single-threaded; above, a plan with at least three primes, threaded
over the primes. There is no size threshold: callers should check
`_flint_mpn_mul_shared_use_fft` first, as `flint_mpn_mul_shared` does,
since for short operands two calls of `flint_mpn_mul` are faster.
`_fft_small_mpn_mul_complex` uses this function for a full product of a
real or imaginary operand by a complex one.
:::
:::

::: type
fft_small_mpn_mul_precomp_struct
fft_small_mpn_mul_precomp_t

The transform of a fixed integer $x$ of `xn` limbs for products with
partners of up to `yn_max` limbs: the transform data `d`, `xstride`
doubles per prime, and the geometry, either two-prime (chunk size `bits`
and `depth`) or, if `wide` is set, given by the plan `P`. `owns` records
whether `d` was allocated by the initialisation.
:::

::: function
ulong _fft_small_mpn_mul_precomp_size(ulong xn, ulong yn_max)
void _fft_small_mpn_mul_precomp_init(fft_small_mpn_mul_precomp_t P, nn_srcptr x, ulong xn, ulong yn_max, double * storage)
void _fft_small_mpn_mul_precomp_clear(fft_small_mpn_mul_precomp_t P)

`_fft_small_mpn_mul_precomp_size` returns the number of doubles that the
transform of an integer of $x_n$ limbs needs for partners of up to
`yn_max` limbs. `_fft_small_mpn_mul_precomp_init` stores the transform of
$(x, x_n)$ in `storage`, which must hold that many doubles and be aligned
to 32 bytes, or in memory it allocates if `storage` is `NULL`.
`_fft_small_mpn_mul_precomp_clear` frees that memory, but not storage
provided by the caller. Require $x_n, \mathtt{yn\_max} \ge 1$.

::: details
The geometry is that of the product of $x_n$ by `yn_max` limbs: two primes
up to `FLINT_MPN_MUL_NP2_MAX_AN` limbs, at least three primes above. Only
the points read by the longest product are kept (its truncation, rather
than the whole transform of length $2^d$, which can be up to twice as
long). The operand $x$ is not referenced after initialisation.
:::
:::

::: function
void _fft_small_mpn_mul_precomp(nn_ptr z, const fft_small_mpn_mul_precomp_t P, nn_srcptr y, ulong yn)

Sets $(z, x_n + y_n)$ to $xy$, where $x$ is the integer whose transform
`P` holds. Requires $1 \le y_n \le$ `yn_max`; `z` may overlap `y`.

::: details
A product converts and transforms $y$ at the truncation of its own
length, multiplies pointwise by the stored transform, transforms back and
reconstructs, which saves the conversion and forward transform of $x$.
`P` is only read, so several threads may multiply with it at the same
time, each using the scratch buffer of its own default context.
:::
:::

## Polynomial multiplication {#poly}

Products of polynomials over $\Z/n\Z$ and $\Z$, as windows of the full
product, and transposed products over $\Z/n\Z$, which give middle
products at the cost of a product of the output length. The
coefficients of the integer product are computed modulo as many primes as
their size requires and reconstructed by the Chinese remainder theorem.

::: function
void _nmod_poly_mul_mid_mpn_ctx(ulong * z, ulong zl, ulong zh, const ulong * a, ulong an, const ulong * b, ulong bn, nmod_t mod, mpn_ctx_t R)
void _nmod_poly_mul_mid_default_mpn_ctx(nn_ptr z, slong zl, slong zh, nn_srcptr a, slong an, nn_srcptr b, slong bn, nmod_t mod)

Sets $(z, z_h - z_l)$ to the coefficients $[z_l, z_h)$ of the product of
$(a, a_n)$ and $(b, b_n)$ modulo $n$. Requires $a_n, b_n \ge 1$, in
either order, and reduced coefficients; any modulus $n \ge 1$ is
allowed.

::: details
Nothing is written if $z_l \ge z_h$, and coefficients from index
$a_n + b_n - 1$ on are zero. With $k$ the bit length of $n$, the number
of primes is the least whose product is at least $m 2^{2k}$, with
$m = \min(a_n, b_n)$: one for moduli up to about
$2^{(50 - \log_2 m)/2}$, two up to about $2^{(99 - \log_2 m)/2}$ and
three for all other 64-bit moduli. A single prime suffices also when $n$
is one of the context primes, or when $n$ is a prime of at least 20 bits
with $2^d \mid n - 1$ for the transform length $2^d$ and
$\max(a_n, b_n) \ge 100$; the transform is then computed modulo $n$
itself. Nothing depends on the order of the operands. See
[](fft_small/poly_mul).
:::
:::

::: function
void _nmod_poly_mul_tr_mpn_ctx(ulong * z, ulong zn, const ulong * a, ulong an, const ulong * v, ulong vn, nmod_t mod, mpn_ctx_t R)
void _nmod_poly_mul_tr_default_mpn_ctx(nn_ptr z, slong zn, nn_srcptr a, slong an, nn_srcptr v, slong vn, nmod_t mod)

Transposed product: sets $z_j = \sum_i a_i v_{i+j} \bmod n$ for
$0 \le j < z_n$, with $(v, v_n)$ read as padded with zeros to length
$a_n + z_n - 1$. This is the transpose of the multiplication by
$(a, a_n)$ of polynomials of length $z_n$; with $a$ the reversal of $f$
it is the window $[a_n - 1, a_n - 1 + z_n)$ of the product $fv$. Requires
$a_n, v_n, z_n \ge 1$, $v_n \le a_n + z_n - 1$ and reduced coefficients;
the output must not overlap the inputs.

::: details
The cost is that of a product of lengths $a_n$ and $z_n$, whereas
`_nmod_poly_mul_mid_mpn_ctx` computes the same window from a product of
lengths $a_n$ and $v_n$, unless a power of two allows a cyclic
convolution. The plan, and with it the primes and the transform modulo
$n$ itself, is that of the product of lengths $a_n$ and $z_n$; the
convolution engine runs in its transposed mode, sending $v$ through
`sd_ifft_trunc_tr`, multiplying pointwise by the transform of $a$ and
applying `sd_fft_trunc_tr`. No wrap-around is used. `nmod_poly_mulmid`
reaches this function through `_nmod_poly_mulmid_fft_small` when its
transform is shorter. See [](fft_small/poly_mul#transposed) and
[](algorithms/polynomials/middle-product#transposed-fft).
:::
:::

::: function
int _fmpz_poly_mul_mid_mpn_ctx(fmpz * z, ulong zl, ulong zh, const fmpz * a, ulong an, const fmpz * b, ulong bn, mpn_ctx_t R)
int _fmpz_poly_mul_mid_default_mpn_ctx(fmpz * z, slong zl, slong zh, const fmpz * a, slong an, const fmpz * b, slong bn)

Sets $(z, z_h - z_l)$ to the coefficients $[z_l, z_h)$ of the product of
$(a, a_n)$ and $(b, b_n)$ over $\Z$ and returns 1, if eight primes
suffice for the coefficient size; otherwise returns 0 and leaves `z`
unchanged. Requires $a_n, b_n \ge 1$; `z` must hold initialised `fmpz`
values.

::: details
With $\beta_a$ and $\beta_b$ the maximal bit sizes of the coefficients
of $a$ and $b$, the product of the primes must be at least
$b_n 2^{\beta_a + \beta_b + 1}$, so the function succeeds roughly when
$\beta_a + \beta_b + \log_2 b_n < 397$. Signed coefficients are lifted to
the symmetric range. See [](fft_small/poly_mul).
:::
:::

::: function
void _nmod_poly_divrem_mpn_ctx(ulong * q, ulong * r, const ulong * a, ulong an, const ulong * b, ulong bn, nmod_t mod, mpn_ctx_t R)

Sets $(q, a_n - b_n + 1)$ and $(r, b_n - 1)$ to the quotient and
remainder of $(a, a_n)$ divided by $(b, b_n)$ modulo $n$. Requires
$a_n \ge b_n \ge 2$ and an invertible leading coefficient of $b$.

::: details
The quotient is a middle product of the reversed dividend with a power
series inverse of the reversed divisor, computed by
`_nmod_poly_inv_series`; the remainder is $a - bq$ computed modulo
$x^N - 1$ with $N$ the least power of two that is at least $b_n - 1$ and
256, which is exact because $\deg r < b_n - 1 \le N$. FLINT's
`nmod_poly_divrem` does not call this function; see
[](fft_small/precomp) for the version with a precomputed divisor.
:::
:::

::: function
void _nmod_poly_mul_mid(ulong * z, slong zl, slong zh, const ulong * a, slong an, const ulong * b, slong bn, nmod_t mod)
void _nmod_poly_mul_mid_classical(ulong * z, slong zl, slong zh, const ulong * a, slong an, const ulong * b, slong bn, nmod_t mod)

Reference implementations of the window product without transforms, used
by the tests and profiles. `_nmod_poly_mul_mid_classical` evaluates each
coefficient as a sum reduced after every term;
`_nmod_poly_mul_mid` trims the operands to the window and splits a long
operand into pieces of the length of the short one before calling it.
Neither is used by FLINT's polynomial arithmetic, which has its own
`_nmod_poly_mulmid`.
:::

## Precomputed operands {#precomp}

When one operand is used in many products, its transform can be computed
once. The second operand is then converted and transformed, multiplied
pointwise and transformed back, which saves about a third of the work.

::: type
mul_precomp_struct

A transformed polynomial over $\Z/n\Z$: a plan for cyclic convolutions of
length $2^d$ and the transform of the operand modulo $x^{2^d} - 1$, with
its length `bn` and the truncation `btrunc` it was transformed for. If the
field `transposed` is set (by `_mul_precomp_init_tr`), it holds instead
the images of a polynomial $v$ through the transposed inverse transform,
for transposed products with $v$.
:::

::: function
int _mul_precomp_init(mul_precomp_struct * M, const ulong * b, ulong bn, ulong btrunc, ulong depth, ulong len_bound, ulong direct_len, nmod_t mod, mpn_ctx_t R)
void _mul_precomp_clear(mul_precomp_struct * M)

Initialises `M` with the transform of $(b, b_n)$ for convolutions of
length $2^d$, $d$ = `depth` $\ge 4$, in which at most `len_bound`
products accumulate in one coefficient, and returns 1; returns 0 if more
than eight primes would be needed. `btrunc` is the truncation of the
transforms `M` will be used with, an upper bound for the output lengths.

::: details
`len_bound` = $2^d$ is always enough; a smaller bound such as
$\min(a_n, b_n)$ for products of linear length at most $2^d$ may save a
prime. `direct_len` enables the transform modulo $n$ itself as in
`fft_small_plan_init_nmod` (0 disables it). An operand longer than $2^d$
is reduced modulo $x^{2^d} - 1$. See [](fft_small/precomp).
:::
:::

::: function
int _nmod_poly_mul_mid_precomp(ulong * z, ulong zl, ulong zh, const ulong * a, ulong an, const mul_precomp_struct * M, nmod_t mod, mpn_ctx_t R)

Sets $(z, z_h - z_l)$ to the coefficients $[z_l, z_h)$ of the product of
$(a, a_n)$ with the polynomial stored in `M` and returns 1, or returns 0
if the window needs a transform longer than $2^d$. Requires $a_n \ge 1$.

::: details
The window fits if $a_n + b_n - 1 \le 2^d$, or if the coefficients beyond
$2^d$ fold onto coefficients below $z_l$, that is, if
$z_h \le 2^d \le a_n + b_n - 1 \le z_l + 2^d$ and $a_n \le 2^d$. `M` is
only read, so several threads may use it at
the same time, each with its own context `R` (normally its default
context); `R` must have the primes of the context `M` was made with.
:::
:::

::: function
int _mul_precomp_init_tr(mul_precomp_struct * M, const ulong * v, ulong vn, ulong len_bound, ulong direct_len, nmod_t mod, mpn_ctx_t R)

Initialises `M` for the transposed products with $(v, v_n)$ computed by
`_nmod_poly_mul_tr_precomp`, in which at most `len_bound` products
accumulate in one coefficient, and returns 1; returns 0 if no plan is
available. `M` is freed with `_mul_precomp_clear`.

::: details
These transposed products are transposes of products of length at most
$v_n$, so one plan serves them all: that of the window $[0, v_n)$ of a
convolution of length $v_n$, whose transform length is $v_n$ rounded up
to a multiple of 256 (for $v_n > 128$) rather than a power of two. `M`
stores the images of $v$ through the transposed inverse transform
(`fft_small_ifft_tr_nmod`). `direct_len` enables the transform modulo $n$
itself as in `fft_small_plan_init_nmod` (0 disables it).
:::
:::

::: function
int _nmod_poly_mul_tr_precomp(ulong * z, ulong zn, const ulong * r, ulong rn, const mul_precomp_struct * M, nmod_t mod, mpn_ctx_t R)

Sets $z_j = \sum_i r_i v_{i+j} \bmod n$ for $0 \le j < z_n$, with
$(v, v_n)$ the polynomial of `M` made by `_mul_precomp_init_tr`, and
returns 1; returns 0, writing nothing, if $r_n + z_n - 1 > v_n$.
Requires $r_n, z_n \ge 1$, at most `len_bound` nonzero coefficients in
$r$ and reduced coefficients.

::: details
$r$ is transformed, multiplied pointwise by the stored images of $v$, and
the transposed forward transform gives the $z_n$ coefficients: two
transforms of the length of `M` instead of the three of
`_nmod_poly_mul_tr_mpn_ctx` (the convolution engine in its second
transposed mode). With $r$ the reversal of $a$, zero padded in front,
this is a window of the product $av$; this is how the evaluation and
extrapolation at geometric progressions use it. `M` is only read, as for
`_nmod_poly_mul_mid_precomp`.
:::
:::

::: type
nmod_poly_divrem_precomp_struct

A divisor $b$ prepared for division with remainder: the transform of the
reversed power series inverse of the reversed $b$, and the transform of
$b$ itself for the cyclic product that gives the remainder.
:::

::: function
void _nmod_poly_divrem_precomp_init(nmod_poly_divrem_precomp_struct * M, const ulong * b, ulong bn, ulong Bn, nmod_t mod, mpn_ctx_t R)
void _nmod_poly_divrem_precomp_clear(nmod_poly_divrem_precomp_struct * M)

Prepares the divisor $(b, b_n)$ for quotients of length at most $B_n$,
that is, for dividends of length at most $b_n + B_n - 1$. Requires
$b_n \ge 2$ and an invertible leading coefficient.
:::

::: function
int _nmod_poly_divrem_precomp(ulong * q, ulong * r, const ulong * a, ulong an, nmod_poly_divrem_precomp_struct * M, nmod_t mod, mpn_ctx_t R)

Sets $(q, a_n - b_n + 1)$ and $(r, b_n - 1)$ to the quotient and
remainder of $(a, a_n)$ by the divisor prepared in `M`, and returns 1.
Requires $a_n \ge b_n$. Returns 0 without computing anything if the
quotient length $a_n - b_n + 1$ exceeds the precision $B_n$ of the stored
inverse or needs a longer transform than the one prepared.
:::

## Plans {#plans}

A plan fixes everything that two transformed operands must agree on to be
combined pointwise: the primes, the transform length $2^d$, the bound on
the output coefficients, and the output window. Operands are transformed
for a plan, combined pointwise any number of times, transformed back once
and converted out. The guide [](fft_small/transform-space) explains how
to use plans and operands.

::: type
fft_small_plan_struct
fft_small_plan_t

A plan. The public fields are `np` (number of primes), `depth` ($d$),
`ztrunc` (the transform truncation, $2^d$ for cyclic plans), the window
`zl`, `zh` of a convolution of length `zn`, the coefficient bound
`bound_c` $\cdot 2^{\mathtt{bound\_e}}$, and for integer plans the chunk
size `bits`.

::: details
The other fields are the context `R`, the offset of the first prime in
it, the transform contexts in use (those of `R`, or a private context
`direct_fft` for a transform modulo $n$ itself), the per-prime stride of
the data in doubles (a multiple of 256), a sign flag, and the
normalisation factors `m[i]` $= (c_i 2^d)^{-1} \bmod p_i$, where $c_i$ is
the reduced Chinese remainder cofactor (1 for one prime). The two-prime
integer plans set `m[i]` to 1 and normalise during reconstruction.
:::
:::

::: function
int fft_small_plan_init_nmod(fft_small_plan_t P, mpn_ctx_t R, ulong zl, ulong zh, ulong zn, ulong xtrunc_max, ulong len_bound, ulong prod_bits, nmod_t mod, ulong direct_len)

Initialises `P` for computing the coefficients $[z_l, z_h)$ of
convolutions of length $z_n$ over $\Z/n\Z$ in which at most `len_bound`
products, each less than $2^{\mathtt{prod\_bits}}$, accumulate in one
coefficient. Returns 0 if more than eight primes would be needed.
Requires $z_l < z_h \le z_n$.

::: details
`xtrunc_max` bounds the operand lengths, rounded up to a multiple of 256.
For a product of lengths $a_n$ and $b_n$, take
`len_bound` $= \min(a_n, b_n)$ and `prod_bits` $= 2 \lceil \log_2 n
\rceil$; for a sum of $K$ products, `len_bound` $= K \min(a_n, b_n)$.
Unreduced inputs, such as sums of operands formed in transform space,
are accounted for by increasing `prod_bits`. If $z_n \le z_l + 2^d$ for a
power of two $2^d$ that is at least $z_h$ and the operand lengths, the
transform has length $2^d$ and the excess coefficients wrap around below
$z_l$. If `mod.n` is a context prime, one prime is used without Chinese
remaindering. If `direct_len` is at least 100, $n$ has at least 20 bits,
$2^d \mid n - 1$ and $n$ is prime, the transform is done modulo $n$
itself; pass 0 to disable this.
:::
:::

::: function
int fft_small_plan_init_nmod_cyclic(fft_small_plan_t P, mpn_ctx_t R, ulong depth, ulong zh, ulong len_bound, ulong prod_bits, nmod_t mod)

Initialises `P` for cyclic convolutions modulo $x^N - 1$, $N = 2^d$
with $d$ = `depth` $\ge 4$, writing the coefficients $[0, z_h)$, where
$z_h \le N$. `len_bound` and `prod_bits` are as for
`fft_small_plan_init_nmod`; a full cyclic product accumulates up to $N$
products. Returns 0 if more than eight primes would be needed.
:::

::: function
int fft_small_plan_init_mpn(fft_small_plan_t P, mpn_ctx_t R, ulong an_max, ulong bn_max, ulong len_bound, int is_signed)

Initialises `P` for products of nonnegative integers of at most
`an_max` and `bn_max` limbs, of which at most `len_bound` are added in
transform space. If `is_signed` is nonzero, the plan also allows
subtractions and the signed exports. Returns 1 on success and 0 if no
prime count and chunk size satisfy the bound.

::: details
For `len_bound` $\le 2$ and operands of at most 40000 limbs the
two-prime geometry is chosen if its exactness bound allows; otherwise a
profile with four to eight primes, scored by the estimated transform
cost. `is_signed` doubles the internal capacity, because the signed
exports reconstruct values in the symmetric range modulo the product of
the primes; callers must not double `len_bound` themselves. The pipeline
choice is made before the doubling, so a signed sum of two products still
gets the two-prime geometry.
:::
:::

::: function
int fft_small_plan_init_mpn_np(fft_small_plan_t P, mpn_ctx_t R, ulong an_max, ulong bn_max, ulong len_bound, int is_signed, ulong np)

Like `fft_small_plan_init_mpn`, with the regime forced for testing and
profiling: `np` = 0 chooses automatically, `np` = 2 forces the two-prime
geometry (and fails if its bound cannot be met), and any other value
forces the profile selection with four to eight primes.
:::

::: function
int fft_small_plan_init_mpn_cyclic(fft_small_plan_t P, mpn_ctx_t R, ulong * nn, ulong len_bound)

Initialises `P` for products of nonnegative integers of at most
$n_n$ = `*nn` limbs modulo $B^{n_n} - 1$. The length is rounded up to the
nearest one for which the chunks split $n_n$ limbs exactly into $2^d$
chunks, and the new value is written to `*nn`. Returns 0 if no geometry
exists.

::: details
Operands longer than $n_n$ limbs must be reduced by the caller. The
export `fft_small_export_mpn` returns a nonnegative value less than
`len_bound` $\cdot 2^{\bits n_n + 2 \mathtt{bits} + d}$, which it writes
exactly into $z_n$ limbs if $\bits z_n \ge \bits n_n + 2\,\mathtt{bits}
+ d + \lceil \log_2(\mathtt{len\_bound} + 1) \rceil$; the caller folds
the limbs above $n_n$ back to reduce modulo $B^{n_n} - 1$.
:::
:::

::: function
int fft_small_plan_init_mpn_window(fft_small_plan_t P, mpn_ctx_t R, ulong an_max, ulong bn_max, ulong len_bound, int is_signed, ulong lo, ulong hi)
ulong fft_small_plan_mpn_nslots(const fft_small_plan_t P, ulong an_max, ulong bn_max)
int fft_small_plan_is_cyclic(const fft_small_plan_t P)

`fft_small_plan_init_mpn_window` is `fft_small_plan_init_mpn` for products
of which only the limbs $[\mathit{lo}, \mathit{hi})$ are needed.
`fft_small_plan_mpn_nslots` returns the number of chunk coefficients
(slots) that the inverse transform computes for operands of at most
`an_max` and `bn_max` limbs, and `fft_small_plan_is_cyclic` whether slots
alias, that is, whether the transform is shorter than the product.

::: details
Only the chunk coefficients that can reach the window are needed. If a
power of two $2^d$ is at least the operand chunk counts and the highest
of them, and the coefficients past $2^d$ fold onto coefficients below the
lowest one, the transforms are cyclic of length $2^d$. The slots below
the window then hold no exact values, and exports must read windows
inside $[\mathit{lo}, \mathit{hi})$ with
`fft_small_export_mpn_signed_window`, which are within one unit of the
exact window. With $\mathit{lo} = 0$ and $\mathit{hi}$ at least the
product length this is the ordinary plan.
:::
:::

::: function
void fft_small_plan_clear(fft_small_plan_t P)

Frees the private transform context of a plan that transforms modulo
$n$ itself, and does nothing otherwise.
:::

## Transformed operands {#ops}

An operand holds the forward transforms of one polynomial or integer
modulo each prime of a plan. It is in the *primal* domain after a forward
transform and in the *product* domain after a pointwise product, which
carries the normalisation factor of the plan. Products can be added and
subtracted in transform space, but the two domains must not be mixed, and
a product cannot be multiplied again. The pointwise operations process the
first `P->ztrunc` points, rounded up to the granularity of the plan's
transforms (whole blocks of 256, or the full length $2^d$ below one
block), which are all the points the transforms produce and read.

::: macro
FFT_SMALL_OP_PRIMAL
FFT_SMALL_OP_PRODUCT

The two values of the `domain` field of an operand.
:::

::: type
fft_small_op_struct
fft_small_op_t

A transformed operand: `np` blocks of `stride` doubles at `data`, the
`depth` and prime `offset` of the plan it was made for, the input
truncation `itrunc`, the `domain`, and whether it owns its data.
:::

::: function
void fft_small_op_init(fft_small_op_t X, const fft_small_plan_t P)
void fft_small_op_init_borrowed(fft_small_op_t X, const fft_small_plan_t P, double * data)
ulong fft_small_op_sizeof_data(const fft_small_plan_t P)
void fft_small_op_clear(fft_small_op_t X)

Initialise an operand for the plan `P`, with storage allocated by the
operand or provided by the caller, and free it. The borrowed storage must
hold `fft_small_op_sizeof_data(P)` bytes, be aligned to 4096 bytes and
outlive the operand; `fft_small_op_clear` does not free it.

::: details
The data size is $\mathtt{np} \cdot \mathtt{stride} \cdot 8$ bytes rounded
up to 4096. For transforms shorter than one block the padding up to the
stride is zeroed by both initialisers. An operand can be used with any
plan that has the same primes, depth and stride, which is checked only by
assertions.
:::
:::

::: function
void fft_small_fft_nmod(fft_small_op_t X, const ulong * a, ulong an, ulong itrunc, nmod_t mod, const fft_small_plan_t P)

Sets `X` to the forward transform of $(a, a_n)$, a polynomial over
$\Z/n\Z$ with reduced coefficients, read as zero-padded to length
`itrunc`. Requires `itrunc` to be a multiple of 8 and at most
`P->ztrunc`; normally `itrunc` is $a_n$ rounded up to a multiple of 256,
or $2^d$ if the depth $d$ is less than 8.

::: details
If $a_n$ exceeds `itrunc`, the polynomial is reduced modulo
$x^{\mathtt{itrunc}} - 1$, which is meaningful only for cyclic plans with
`itrunc` $= 2^d$. The primes are transformed in parallel from
$a_n \ge 1500$.
:::
:::

::: function
void fft_small_fft_mpn(fft_small_op_t X, const ulong * a, ulong an, const fft_small_plan_t P)

Sets `X` to the forward transform of the integer $(a, a_n)$, cut into
chunks of `P->bits` bits. Requires an integer plan and $a_n \ge 1$, at
most the operand length the plan was made for.

::: details
Two-prime plans pack the chunks with plain bit slicing. Other plans use
the vectorised packing of their profile when one exists, and a scalar
conversion with tables of powers of two otherwise. The work is spread
over the primes from $a_n \ge 1024$.
:::
:::

::: function
void _fft_small_fft_mpn_trunc(fft_small_op_t X, const ulong * a, ulong an, ulong otrunc, const fft_small_plan_t P)

Like `fft_small_fft_mpn`, with the transforms truncated to their first
`otrunc` points, which suffices for products of at most `otrunc` slots.
Requires `otrunc` $\le$ `P->ztrunc`, a multiple of 256 or equal to
`P->ztrunc`; `fft_small_fft_mpn` is the case `otrunc` = `P->ztrunc`.
:::

::: function
void fft_small_op_mul(fft_small_op_t Z, const fft_small_op_t A, const fft_small_op_t B, const fft_small_plan_t P)
void fft_small_op_sqr(fft_small_op_t Z, const fft_small_op_t A, const fft_small_plan_t P)

Set `Z` to the pointwise product $AB$, respectively $A^2$, of operands in
the primal domain. `Z` is then in the product domain and may alias `A` or
`B`.
:::

::: function
void fft_small_op_addmul(fft_small_op_t Z, const fft_small_op_t A, const fft_small_op_t B, const fft_small_plan_t P)
void fft_small_op_submul(fft_small_op_t Z, const fft_small_op_t A, const fft_small_op_t B, const fft_small_plan_t P)
void fft_small_op_addsqr(fft_small_op_t Z, const fft_small_op_t A, const fft_small_plan_t P)
void fft_small_op_subsqr(fft_small_op_t Z, const fft_small_op_t A, const fft_small_plan_t P)

Add or subtract the pointwise product $AB$, respectively $A^2$, of
primal operands to `Z`, which must be in the product domain. The `sqr`
versions read `A` once and are preferable to `addmul` with $B = A$.

::: details
A subtraction is exact only if the exported coefficients stay
nonnegative, unless the plan is signed and a signed export is used: the
unsigned exports reconstruct values in $[0, P)$, with $P$ the product of
the primes.
:::
:::

::: function
void fft_small_op_add(fft_small_op_t Z, const fft_small_op_t A, const fft_small_op_t B, const fft_small_plan_t P)
void fft_small_op_sub(fft_small_op_t Z, const fft_small_op_t A, const fft_small_op_t B, const fft_small_plan_t P)
void fft_small_op_neg(fft_small_op_t Z, const fft_small_op_t A, const fft_small_plan_t P)

Set `Z` to $A + B$, $A - B$ and $-A$, for operands in the same domain;
`Z` gets that domain. Sums of primal operands represent sums of the
inputs, whose larger coefficients the plan's bound must allow for.
:::

::: function
void fft_small_ifft(fft_small_op_t X, const fft_small_plan_t P)

Applies the inverse transform in place to an operand in the product
domain, at the truncation of the plan. The normalisation was applied by
the pointwise products. The primes are processed in parallel from
`P->ztrunc` $\ge 3072$.
:::

::: function
void fft_small_ifft_trunc(fft_small_op_t X, ulong itrunc, ulong clear_to, const ulong * scale, const fft_small_plan_t P)

Inverse transform of the first `itrunc` slots in place, with the
truncation rounded up as for the plan. If `scale` is not `NULL`,
`scale[i]` is a normalisation factor modulo the $i$-th prime that was
not applied in the pointwise stage, and is applied after the transform.
The slots $[\mathtt{itrunc}, \mathtt{clear\_to})$ are set to zero, for
exports that read whole chunks. `X` is marked as being in the product
domain.
:::

::: function
void fft_small_ifft_tr_nmod(fft_small_op_t X, const ulong * a, ulong an, nmod_t mod, const fft_small_plan_t P)
void fft_small_fft_tr(fft_small_op_t X, ulong itrunc, const fft_small_plan_t P)

The transposes of `fft_small_ifft` and `fft_small_fft_nmod`, for
transposed products. `fft_small_ifft_tr_nmod` converts $(a, a_n)$, read
as zero-padded to `P->ztrunc` coefficients (so $a_n \le$ `P->ztrunc`), and
applies the transposed inverse transform; `X` is then in the primal
domain. `fft_small_fft_tr` applies the transposed forward transform with
truncation `itrunc` $\le$ `P->ztrunc` in place to an operand in the
product domain, whose first `itrunc` slots can then be exported.

::: details
For a plan of the window $[0, z_n)$ of a convolution of length
$a_n + z_n - 1$, with `B` the transform of $(a, a_n)$ made by
`fft_small_fft_nmod`, the calls

    fft_small_ifft_tr_nmod(X, v, vn, mod, P);
    fft_small_op_mul(X, X, B, P);
    fft_small_fft_tr(X, ztrunc, P);
    fft_small_export_nmod_range(z, X, 0, zn, mod, P);

compute the transposed product of `_nmod_poly_mul_tr_mpn_ctx`, where
`ztrunc` is the truncation for $z_n$ (`_fft_small_op_trunc`). The
normalisation is applied by the pointwise product, as for ordinary
products. The primes are processed in parallel as in
`fft_small_fft_nmod` and `fft_small_ifft`. See
[](fft_small/transform-space#transposed).
:::
:::

## Conversion out of transform space {#exports}

After `fft_small_ifft`, the slots hold the coefficients of the result
modulo each prime. The exports reconstruct them by the Chinese remainder
theorem. The `mpn` exports may overwrite the operand (the two-prime ones
always do), which is why they take it as non-`const`; copy it first if
it is needed afterwards.

::: function
void fft_small_export_nmod(ulong * z, const fft_small_op_t X, nmod_t mod, const fft_small_plan_t P)
void fft_small_export_nmod_range(ulong * z, const fft_small_op_t X, ulong zl, ulong zh, nmod_t mod, const fft_small_plan_t P)

Write the coefficients of the plan's window $[z_l, z_h)$, respectively
of the subwindow $[z_l, z_h) \subseteq$ [`P->zl`, `P->zh`), reduced modulo
$n$, to `z`, which receives $z_h - z_l$ coefficients. The exact integer
coefficients must lie in $[0, P)$ for the product $P$ of the plan's
primes, as guaranteed by its bound for sums of products.
:::

::: function
void fft_small_export_mpn(ulong * z, ulong zn, fft_small_op_t X, const fft_small_plan_t P)

Writes to $(z, z_n)$ the integer $\sum_k s_k 2^{k \cdot \mathtt{bits}}$
reduced modulo $B^{z_n}$, where $s_k \in [0, P)$ are the reconstructed
slots. For a product of integers this is the product, exact if $z_n$ is
at least its length.
:::

::: function
void fft_small_export_mpn_signed(ulong * z, ulong zn, int * sign, fft_small_op_t X, ulong nslots, const fft_small_plan_t P)
void fft_small_export_mpn_signed_trunc(ulong * z, ulong zn, int * sign, fft_small_op_t X, ulong nslots, ulong lo_limbs, const fft_small_plan_t P)

Reconstruct a signed integer from the first `nslots` slots, each lifted
to the symmetric range $(-P/2, P/2]$. The absolute value is written to
$(z, z_n)$, zero-padded, and its sign to `*sign`. The `trunc` version
writes the limbs from limb `lo_limbs` on, within $(-1.5, 0.5)$ units of
its lowest limb; with `lo_limbs` = 0 both are exact.

::: details
Requires a plan initialised with `is_signed` set and
$\bits (\mathtt{lo\_limbs} + z_n) \ge \mathtt{nslots} \cdot \mathtt{bits}
+ \bits c + 2$, where $c$ is the limb count of the Chinese remainder data
(`fft_small_export_mpn_signed_limbs` computes the least $z_n$ for
`lo_limbs` = 0). The truncated version processes only the slots that can
reach limb `lo_limbs`: the error is the truncation itself, below one
unit, plus the dropped slots, less than half a unit either way.
:::
:::

::: function
void fft_small_export_mpn_signed_low(ulong * z, ulong zn, fft_small_op_t X, ulong nslots, const fft_small_plan_t P)
void fft_small_export_mpn_signed_window(ulong * z, ulong zn, fft_small_op_t X, ulong nslots, ulong lo_limbs, const fft_small_plan_t P)

Write to $(z, z_n)$ the limbs $[\mathtt{lo\_limbs}, \mathtt{lo\_limbs} +
z_n)$ of the two's complement value of the signed reconstruction, without
resolving the sign; `_low` is the case `lo_limbs` = 0, which is the
residue modulo $B^{z_n}$ and exact. The value need not fit in the window.

::: details
Only the slots that start below limb $\mathtt{lo\_limbs} + z_n$ are read,
since carries and borrows only propagate upwards. For `lo_limbs` $> 0$
the slots below the first one that can reach limb `lo_limbs` are dropped
as in the truncated export, and the window is within one unit of the
exact one. This holds even if those slots hold arbitrary values, as the
aliased slots of a cyclic window plan do, because every lifted slot is
less than $P/2$ in absolute value.
:::
:::

::: function
ulong fft_small_export_mpn_signed_limbs(const fft_small_plan_t P, ulong nslots)
ulong fft_small_export_mpn_limbs(const fft_small_plan_t P, ulong len_bound, ulong nslots)

Return the number of limbs that the signed reconstruction of `nslots`
slots, respectively the unsigned one with up to `len_bound` accumulated
products per slot, writes. An export into a shorter window must go
through scratch memory, as `fft_small_stage_begin` arranges.
:::

## Operand pools and staged output {#pools}

Drivers that keep many transformed operands alive allocate their storage
in pools, and stage conversions out through scratch memory taken from the
same pool.

::: macro
FLINT_FFT_SMALL_ALIGNMENT

The alignment in bytes of transform data, 4096.
:::

::: variable
ulong flint_fft_small_max_transformed_ring_size

The largest amount of transform storage, in bytes, that a single
transformed ring or context may plan for: its expected number of
simultaneously live elements times the size of one element. Constructors
decline larger requests, and their callers fall back to other algorithms
(matrix multiplication then works in blocks). The default is
$4 \cdot 2^{30}$ on 64-bit machines and $2^{30}$ otherwise; it may be
changed. It is defined even without `FLINT_HAVE_FFT_SMALL`.
:::

::: macro
FFT_SMALL_OP_POOL_RESERVE
FFT_SMALL_OP_POOL_ALLOC

Flags for `fft_small_op_pool_fit`: take the block as a reservation in the
context's scratch buffer, respectively allocate it.
:::

::: type
fft_small_op_pool_struct
fft_small_op_pool_t

A block of equal slabs for operand data plus an optional staging region,
with a stack of free slabs protected by a mutex.
:::

::: function
void fft_small_op_pool_init(fft_small_op_pool_t pool)
int fft_small_op_pool_fit(fft_small_op_pool_t pool, mpn_ctx_t R, ulong slab_size, slong count, ulong stage_bytes, int flags)
void fft_small_op_pool_clear(fft_small_op_pool_t pool)

Initialise an empty pool, give it `count` slabs of `slab_size` bytes and
a staging region of `stage_bytes` bytes, and free it. `fit` returns 0,
leaving the pool empty, if the block can be neither reserved nor
allocated as `flags` allow.

::: details
`slab_size` must be a multiple of `FLINT_FFT_SMALL_ALIGNMENT`, normally
`fft_small_op_sizeof_data` of the plan. With `FFT_SMALL_OP_POOL_RESERVE`
the block is reserved in the scratch buffer of `R`
(`mpn_ctx_fit_buffer_reserve`), memory that is usually already mapped;
with `FFT_SMALL_OP_POOL_ALLOC` it is allocated if no reservation is
possible. A reserved pool must be cleared on the thread that owns `R`,
after every slab has been given back.
:::
:::

::: function
double * fft_small_op_pool_take(fft_small_op_pool_t pool)
void fft_small_op_pool_give(fft_small_op_pool_t pool, double * d)
int fft_small_op_pool_owns(const fft_small_op_pool_t pool, const void * d)

Take a free slab, for `fft_small_op_init_borrowed`, or `NULL` if none is
free; give a slab back; and test whether an address lies in the pool's
block. Taking and giving are thread safe.
:::

::: function
ulong * fft_small_op_pool_scratch_take(fft_small_op_pool_t pool, ulong limbs, int use_slabs)
void fft_small_op_pool_scratch_give(fft_small_op_pool_t pool, ulong * t)

Return scratch space for `limbs` limbs and give it back. The space is the
staging region if it is free and large enough, otherwise a free slab if
`use_slabs` is nonzero and one is large enough, otherwise a fresh
allocation; concurrent requests never share it.
:::

::: function
ulong * fft_small_stage_begin(ulong * z, ulong zn, ulong need, fft_small_op_pool_struct * pool, int use_slabs)
int fft_small_stage_end(ulong * z, ulong zn, ulong * buf, ulong len, fft_small_op_pool_struct * pool)

Staged conversion of a value whose reconstruction writes `need` limbs
into a window of $z_n$ limbs. `begin` returns `z` if $z_n \ge$ `need`
and scratch space from the pool otherwise (`pool` may be `NULL`, meaning
a plain allocation). After exporting there, and possibly adjusting the
value, `end` copies the normalised value of length `len` to $(z, z_n)$,
zero pads it, releases the scratch space, and returns 0 if it does not
fit (the window is then undefined).
:::

## Complex integer products {#complex}

The products behind the `flint_mpn_*_complex` functions of `mpn_extras`,
which multiply Gaussian integers $(a_r + i a_i)(b_r + i b_i)$ with parts
given as magnitudes and sign bits. They are used above 230 limbs (270 for
squares).

::: function
int _fft_small_mpn_mulmid_complex(nn_ptr zr, nn_ptr zi, slong lo, slong hi, nn_srcptr ar, slong arn, int ar_sgn, nn_srcptr ai, slong ain, int ai_sgn, nn_srcptr br, slong brn, int br_sgn, nn_srcptr bi, slong bin, int bi_sgn)
int _fft_small_mpn_sqrmid_complex(nn_ptr zr, nn_ptr zi, slong lo, slong hi, nn_srcptr ar, slong arn, int ar_sgn, nn_srcptr ai, slong ain, int ai_sgn)

Write the limbs $[\mathit{lo}, \mathit{hi})$ of the real and imaginary
parts of the product, respectively the square, as two's complement
windows of $\mathit{hi} - \mathit{lo}$ limbs. Each is within one unit of
$\lfloor x / B^{\mathit{lo}} \rfloor$ for the exact part $x$, and exact
(the residue modulo $B^{\mathit{hi}}$) for $\mathit{lo} = 0$. Return 0 if
no plan exists.

::: details
A part is given by a magnitude of the given length (0 for a zero part) and
a sign bit (1 for negative). The window plan is used, with cyclic
transforms when the window allows.
:::
:::

::: function
int _fft_small_mpn_mul_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, slong lo, slong hi, nn_srcptr ar, slong arn, int ar_sgn, nn_srcptr ai, slong ain, int ai_sgn, nn_srcptr br, slong brn, int br_sgn, nn_srcptr bi, slong bin, int bi_sgn)
int _fft_small_mpn_sqr_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, slong lo, slong hi, nn_srcptr ar, slong arn, int ar_sgn, nn_srcptr ai, slong ain, int ai_sgn)

Write the magnitudes of the limbs $[\mathit{lo}, \mathit{hi})$ of the real
and imaginary parts of the product, respectively the square, to `zr` and
`zi` ($\mathit{hi} - \mathit{lo}$ limbs each, enough to hold the parts
shifted right by $\mathit{lo}$ limbs), and their lengths, negated for
negative parts, to `*zr_len` and `*zi_len`. Exact for $\mathit{lo} = 0$
and otherwise within $(-1.5, 0.5)$ units of the lowest limb. Return 0 if
no plan exists.

::: details
Only the nonzero parts are transformed, so a real or imaginary operand
costs one transform for both output parts. Each output part is a sum of
at most two products, formed with its signs in transform space and
exported signed; an output part that is a single product with
$\mathit{lo} = 0$ is exported unsigned with its known sign, and if all
are, the plan is unsigned and packs more bits per chunk. Storage for five
transforms is needed for a product and three for a square. A full product
($\mathit{lo} = 0$, not a square) of an operand with a single nonzero part
by a complex one consists of two products with a common factor and is
passed to `_fft_small_mpn_mul_shared` when the outputs do not overlap the
operands in a way that function forbids.
:::
:::

## Multiplication modulo $2^N + 1$ {#negacyclic}

Multiplication of residues modulo $2^N + 1$ by a negacyclic convolution
of their digits, computed with the small-prime transforms. This is the
pointwise stage of the Schönhage–Strassen convolutions in the `fft`
module (`fft_convolution_sd_fft`), which `fmpz_poly` and `fmpzi`
multiplication use for large coefficients.

Residues occupy $n + 1$ limbs, $n = N/\bits$. An input may have top limb
0 or 1 with arbitrary other limbs and is read modulo $2^N + 1$. The output
is fully reduced: the unique representative in $[0, 2^N]$, with top limb 1
only for $2^N$.

::: type
sd_fft_mpn_mulmod_2expp1_ctx_struct
sd_fft_mpn_mulmod_2expp1_ctx_t
sd_fft_mpn_mulmod_2expp1_scratch_struct
sd_fft_mpn_mulmod_2expp1_scratch_t

A context for one ring $\Z/(2^N + 1)\Z$, holding $N$, the digit size $b$,
the digit count $m = N/b$, the number of primes and the twisting tables,
and a scratch object with the mutable state of one multiplication.
:::

::: function
void sd_fft_mpn_mulmod_2expp1_ctx_init(sd_fft_mpn_mulmod_2expp1_ctx_struct * C, mpn_ctx_struct * R, slong N, slong m)
void sd_fft_mpn_mulmod_2expp1_ctx_clear(sd_fft_mpn_mulmod_2expp1_ctx_struct * C)

Initialises `C` for multiplication modulo $2^N + 1$ with $m$ digits of
$b = N/m$ bits, taking the primes from `R`, normally
`get_default_mpn_ctx()`, and frees it. $N$ must be a positive multiple
of `FLINT_BITS`, and $m$ a power of two at least 16 dividing $N$, or 0 to
let the function choose. Throws if $N$ is not a positive multiple of
`FLINT_BITS` or no configuration exists.

::: details
The number of primes is the least $k \in [2, 8]$ whose product exceeds
$2^{2b + \log_2 m + 1}$ and for which digits of $b$ bits can be
converted: any $b \le 50$, or $b > 50$ with a vectorised conversion for
$k$ primes, which exists for $b \in \{64, 84, 88, 92, 112, 116, 120, 126,
128, 136, 140, 144, 160, 164, 168, 184, 188, 192\}$ and suitable $k$.
With $m = 0$ the digit count is chosen among the powers of two $m \ge 16$
with $m \le N/20$, $m \mid N$ and $20 \le b \le 200$ by an estimated cost
$w_k \, k\, m (\log_2 m + 4 + b/16)$, with weights $w_k$ measured per
prime count. For example $N = 2^{10}, \ldots, 2^{24}$ give $b = 64$ with
three primes, $N = 47104$ gives $b = 92$ with four primes, and
$N = 344064$ gives $b = 42$ with two primes. That $N$ is a multiple of
`FLINT_BITS` is not checked: other values give wrong results.
:::
:::

::: function
void sd_fft_mpn_mulmod_2expp1_scratch_init(sd_fft_mpn_mulmod_2expp1_scratch_struct * S, const sd_fft_mpn_mulmod_2expp1_ctx_struct * C)
void sd_fft_mpn_mulmod_2expp1_scratch_clear(sd_fft_mpn_mulmod_2expp1_scratch_struct * S)

Initialise and free the scratch object for multiplications with `C`. The
context is read-only during multiplication, so several threads can share
it, each with its own scratch object.
:::

::: function
void sd_fft_mpn_mulmod_2expp1(sd_fft_mpn_mulmod_2expp1_ctx_struct * C, nn_ptr z, nn_srcptr x, nn_srcptr y, sd_fft_mpn_mulmod_2expp1_scratch_struct * S)

Sets $z = xy \bmod 2^N + 1$. `z` may alias `x` or `y`; if `x == y` a
square is computed with one forward transform.
:::

::: function
slong sd_fft_mpn_mulmod_2expp1_transformed_size(const sd_fft_mpn_mulmod_2expp1_ctx_struct * C)
void sd_fft_mpn_mulmod_2expp1_transform(sd_fft_mpn_mulmod_2expp1_ctx_struct * C, double * F, nn_srcptr y, sd_fft_mpn_mulmod_2expp1_scratch_struct * S)
void sd_fft_mpn_mulmod_2expp1_mul_cached(sd_fft_mpn_mulmod_2expp1_ctx_struct * C, nn_ptr z, nn_srcptr x, const double * F, sd_fft_mpn_mulmod_2expp1_scratch_struct * S)

Multiplication by a fixed operand. `transform` stores the transform of $y$
in `F`, an array of `transformed_size(C)` doubles; `mul_cached` then sets $z = xy \bmod 2^N + 1$, like
`sd_fft_mpn_mulmod_2expp1` with one forward transform less.
:::

::: function
slong sd_fft_mpn_mulmod_2expp1_choose_np(mpn_ctx_struct * R, slong b, slong depth)

Returns the number of primes the context uses for digits of $b$ bits and
$m = 2^{\mathtt{depth}}$ digits, as described for the context, or $-1$ if
there is none.
:::

## Transforms modulo one prime {#sd-fft}

The number-theoretic transforms themselves. Data is an array of doubles
holding integers, transformed in place; the algorithm and the bounds are
described in [](algorithms/fft/small-prime-fft). These functions are internal, but
stable enough to be used directly; see [](fft_small/sd_fft).

::: macro
LG_BLK_SZ
BLK_SZ

The block size $256 = 2^8$ in doubles. Transforms longer than one block
work on whole blocks in their outer layers, and truncations are rounded
up to whole blocks.
:::

::: type
sd_fft_ctx_struct
sd_fft_ctx_t

The transform context of one prime $p$: `p` and `pinv` $= 1/p$ as
doubles, `mod` for integer arithmetic modulo $p$, a primitive root of
unity `primitive_2power_root` of order $2^v$ with $2^v \| p - 1$, and the
twiddle tables `w2tab` filled to `w2tab_depth`, guarded by a mutex.

::: details
The twiddle factors are the roots of unity $w_j$ (see `sd_fft_ctx_w`)
with even index, stored in $[-p/2, p/2]$ in tables of doubling length:
`w2tab[0]` holds $w_0$, `w2tab[1]` holds $w_2$, and `w2tab[k]` for
$k \ge 2$ holds the $2^{k-1}$ values $w_{2j}$ with $2^{k-1} \le j < 2^k$.
The tables up to `SD_FFT_CTX_W2TAB_INIT` are consecutive in memory.
Tables that have been filled never move, so readers need no lock.
:::
:::

::: macro
SD_FFT_CTX_W2TAB_INIT
SD_FFT_CTX_W2TAB_SIZE

The number of twiddle tables filled at initialisation (12, for transforms
up to length $2^{12}$) and the maximal number of tables (40).
:::

::: function
void sd_fft_ctx_init_prime(sd_fft_ctx_t Q, ulong pp)
void sd_fft_ctx_clear(sd_fft_ctx_t Q)

Initialise the context of the prime $p$ = `pp` and free it. Throws if $p$
fails `fft_small_mulmod_satisfies_bounds`, if $p - 1$ is odd or if
$2^4 \nmid p - 1$. Primality is assumed, not checked.
:::

::: function
void sd_fft_ctx_fit_depth(sd_fft_ctx_t Q, ulong depth)
void sd_fft_ctx_fit_depth_with_lock(sd_fft_ctx_t Q, ulong k)

Make sure that the twiddle tables cover transforms of length
$2^{\mathtt{depth}}$, filling the missing tables under the context's
mutex. Throw if $2^{\mathtt{depth}} \nmid p - 1$. `sd_fft_trunc` and
`sd_ifft_trunc` call this themselves.
:::

::: function
void sd_fft_trunc(sd_fft_ctx_t Q, double * d, ulong L, ulong itrunc, ulong otrunc)

Truncated forward transform of length $2^L$ in place: with $f$ the
polynomial with coefficients `d[0]`, ..., `d[itrunc - 1]`, sets
`d[sd_fft_ctx_trunc_index(L, i)]` to a value congruent to $f(w_i)$ for
$0 \le i < \mathtt{otrunc}$. The inputs must lie in $(-3p, 3p)$, and so
do the outputs.

::: details
`d` must be aligned to 32 bytes and have room for $2^L$ doubles. For
$L \le 8$ the entries from
`itrunc` on are set to zero and the whole transform is computed; for
larger $L$, `itrunc` and `otrunc` are rounded up to multiples of 256, the
entries from `itrunc` to the next multiple of 256 are set to zero, and
only the blocks below the rounded `otrunc` are written. Any
$1 \le \mathtt{itrunc}, \mathtt{otrunc} \le 2^L$ is allowed. The
evaluation points do not depend on $L$, so if both truncations are at
most $2^{L-1}$ the result is that of depth $L - 1$ at about the same
cost; recovering the input with `sd_ifft_trunc` requires
$\mathtt{otrunc} \ge \mathtt{itrunc}$. See [](fft_small/sd_fft).
:::
:::

::: function
void sd_ifft_trunc(sd_fft_ctx_t Q, double * d, ulong L, ulong trunc)

Truncated inverse transform of length $2^L$ in place: if `d[i]` holds
the value at position $i$ of the forward transform of a polynomial $g$ of
length at most `trunc`, for $i < \mathtt{trunc}$, sets `d[i]` to a value
congruent to $2^L g_i$ for $i < \mathtt{trunc}$. The inputs must lie in
$(-2p, 2p)$; the outputs lie in $(-3p, 3p)$.

::: details
For $L \le 8$ the whole inverse transform is computed and `trunc` is
ignored, so all $2^L$ values must be present. For larger $L$, `trunc` is
rounded up to a multiple of 256 and $g$ may have that length. The
outputs are not normalised; the multiplication drivers fold $2^{-L}$ into
the pointwise products.
:::
:::

::: function
void sd_fft_trunc_tr(sd_fft_ctx_t Q, double * d, ulong L, ulong itrunc, ulong otrunc)
void sd_ifft_trunc_tr(sd_fft_ctx_t Q, double * d, ulong L, ulong trunc)

The transposes of `sd_fft_trunc` and `sd_ifft_trunc` with the same
parameters, in place: for all vectors $a$ and $v$,
$\langle \mathtt{sd\_fft\_trunc}(a), v \rangle =
\langle a, \mathtt{sd\_fft\_trunc\_tr}(v) \rangle$, and likewise for the
inverse transform. `sd_fft_trunc_tr` reads the first `otrunc` entries,
indexed like the output of `sd_fft_trunc`, and sets
$d_t$ to a value congruent to
$\sum_{i < \mathtt{otrunc}} d[\mathtt{sd\_fft\_ctx\_trunc\_index}(L, i)]\, w_i^t$
for $0 \le t < \mathtt{itrunc}$. Each costs the same as the transform it
transposes.

::: details
`sd_fft_trunc_tr` accepts inputs in $(-2p, 2p)$ and produces outputs in
$(-2p, 2p)$, so it can follow `sd_fft_ctx_point_mul`; `sd_ifft_trunc_tr`
accepts and produces values in $(-3p, 3p)$, like `sd_fft_trunc`. The
truncations are rounded as for the untransposed transforms: for $L > 8$
to multiples of 256 (`sd_fft_trunc_tr` first zeroes the entries from
`otrunc` to the next multiple of 256), and for $L \le 8$
the full transform is computed (`sd_fft_trunc_tr` zeroes the entries from
`otrunc` on, and `sd_ifft_trunc_tr` ignores `trunc`). The transposed
inverse transform followed by a pointwise product and the transposed
forward transform computes a transposed convolution; see
[](fft_small/sd_fft#transposed) and
[](algorithms/fft/small-prime-fft#transposed).
:::
:::

::: function
void sd_fft_ctx_point_mul(const sd_fft_ctx_t Q, double * a, const double * b, ulong m_, ulong depth)
void sd_fft_ctx_point_sqr(const sd_fft_ctx_t Q, double * a, ulong m_, ulong depth)

Set $a_i \gets m a_i b_i$, respectively $m a_i^2$, modulo $p$ for
$i < 2^{\mathtt{depth}}$, with $0 \le m \le p$ and depth at least 4.
Inputs in $(-3p, 3p)$ give outputs in $(-3p/2, 3p/2)$, valid inputs for
`sd_ifft_trunc`.
:::

::: function
double sd_fft_ctx_w(const sd_fft_ctx_t Q, ulong j)
double sd_fft_ctx_w2(const sd_fft_ctx_t Q, ulong j)
double sd_fft_ctx_w2inv(const sd_fft_ctx_t Q, ulong j)

Return the root of unity $w_j$, $w_{2j}$, and $-w_{2j}^{-1}$
modulo $p$, as integers in $[-p/2, p/2]$. Here $w_0 = 1$ and, for
$2^{k-1} \le j < 2^k$, $w_j = \omega^{r}$ where $\omega$ is a primitive
$2^k$-th root of unity and $r$ is $j$ with its $k$ bits reversed; the
roots for different $k$ are chosen compatibly, $\omega_{2^k} =
\omega_{2^{k+1}}^2$. The tables must cover the index.

::: details
So $w_0, w_1, w_2, w_3, \ldots = 1, -1, \omega_4, -\omega_4, \omega_8,
\omega_8^5, \omega_8^3, \omega_8^7, \ldots$, and the first $2^L$ of
them are the $2^L$-th roots of unity in bit-reversed order. Since
$w_{2j+1} = -w_{2j}$, only even indices are stored.
:::
:::

::: function
ulong sd_fft_ctx_trunc_index(ulong L, ulong i)

The position of $f(w_i)$ in the output of a forward transform of length
$2^L$: $i$ itself for $L < 4$, and otherwise $i$ with the two 2-bit
fields of its low four bits exchanged, which transposes each group of 16
values as a $4 \times 4$ matrix.
:::

::: function
ulong sd_fft_ctx_blk_offset(ulong I)
double * sd_fft_ctx_blk_index(double * d, ulong I)
ulong sd_fft_ctx_data_size(ulong L)

The offset $256 I$ of block $I$, the address `d` $+ 256 I$, and the
number $2^L$ of doubles of a transform of length $2^L$.
:::

::: macro
SET_J_BITS_AND_J_R(j_bits, j_r, j)
SET_J_BITS_AND_J_MR(j_bits, j_mr, j)

Split a twiddle index $j$ into the table number `j_bits` (the bit length
of $j$) and the index in that table: `j_r` $= j - 2^{\mathtt{j\_bits} -
1}$ for $w_{2j}$, and `j_mr` $= 2^{\mathtt{j\_bits}} - 1 - j$, which
locates $-w_{2j}^{-1}$. Both are 0 for $j = 0$.
:::

::: function
int fft_small_mulmod_satisfies_bounds(ulong n)

Returns whether the floating-point modular multiplication modulo $n$
maps products in $(-2n^2, 2n^2)$ into $(-n, n)$ and products in
$(-4n^2, 4n^2)$ into $(-3n/2, 3n/2)$, by a bound that uses the actual
error of the double $1/n$. Every prime below $2^{49}$ passes, and about
one 50-bit prime in 900 fails.

::: summary algorithms/simd/fp-modular-arithmetic :::
:::

::: function
int fft_small_prime_supports_depth(ulong n, ulong depth)

Returns whether a transform of length $2^{\mathtt{depth}}$ can be computed
modulo $n$ itself: $2 < n \le 2^{50}$, $n$ passes
`fft_small_mulmod_satisfies_bounds`, depth is at most 40 and
$2^{\max(\mathtt{depth}, 4)} \mid n - 1$ (the context initialisation
`sd_fft_ctx_init_prime` needs $2^4 \mid n - 1$ for every depth).
Primality is assumed.
:::

## Chinese remaindering and integer conversion {#crt}

Data and helpers behind the conversions of the multiplication drivers.
They are internal.

::: type
crt_data_struct
crt_data_t

Chinese remaindering data for the first `nprimes` primes $p_0, \ldots,
p_{k-1}$ of a context, with product $P$: the cofactors $P / p_i$ and $P$,
each of `coeff_len` limbs (enough to hold $kP$), and the reduced
cofactors $(P / p_i) \bmod p_i$. `prime` is $p_{k-1}$.
:::

::: function
void crt_data_init(crt_data_t C, ulong prime, ulong coeff_len, ulong nprimes)
void crt_data_clear(crt_data_t C)

Allocate storage for the data, which the caller fills in, and free it.
:::

::: function
ulong * crt_data_co_prime(const crt_data_t C, ulong i)
ulong * _crt_data_co_prime(const crt_data_t C, ulong i, ulong n)
ulong * crt_data_prod_primes(const crt_data_t C)
ulong * crt_data_co_prime_red(const crt_data_t C, ulong i)

Pointers to the cofactor $P / p_i$ (`coeff_len` limbs; the second form
takes `coeff_len` as `n` for the compiler), to $P$ (`coeff_len` limbs),
and to the single limb $(P / p_i) \bmod p_i$.
:::

::: macro
MPN_CTX_NCRTS
MAX_NPROFILES
VEC_SZ
MPN_CTX_TWO_POWER_TAB_SIZE

The number of primes of a context (8), the capacity of its profile table
(20), the lane count of the vectors used in the conversions (4), and the
length of the tables of powers of two (256), which must grow with the
number of primes.
:::

::: type
profile_entry_struct
profile_entry_t
to_ffts_func

A packing profile: the number of primes `np`, the chunk size `bits`, the
largest second operand length `bn_bound` for which $\lceil 64 b_n /
\mathtt{bits} \rceil 2^{2 \mathtt{bits}} \le P$, and the vectorised
conversion `to_ffts` of type `to_ffts_func`, which reduces the chunks of
an integer modulo all primes at once.

::: details
The default context has the profiles (np, bits) = (4, 84), (4, 88),
(4, 92), (5, 112), (5, 116), (5, 120), (6, 136), (6, 140), (6, 144),
(7, 160), (7, 164), (7, 168), (8, 184), (8, 188), (8, 192); their bounds
are listed on [](algorithms/fft/multimodular#profiles).
:::
:::

::: type
from_ffts_func

The type of the windowed Chinese remainder reconstructions of the integer
drivers, which turn residues of chunk coefficients into limbs.
:::

::: function
to_ffts_func _mpn_ctx_to_ffts_func(ulong np, ulong bits)
from_ffts_func _mpn_ctx_from_ffts_func(ulong np)

Return the vectorised conversion of integers into chunks of `bits` bits
modulo `np` primes, or `NULL` if none is instantiated, and the
reconstruction for $3 \le$ `np` $\le 8$ primes. Used by the negacyclic
engine.
:::

::: function
void _convert_block(ulong * Xs, sd_fft_ctx_struct * Rffts, double * d, ulong dstride, ulong np, ulong Iv)
void _convert_block_1_mod(ulong * Xs, sd_fft_ctx_struct * Rffts, double * d, ulong dstride, ulong I, double p2, double p2inv)

For each prime $l <$ `np`, reduce block $I$ of the data `d` $+ l \cdot$
`dstride` to $[0, p_l)$ and store it as integers at `Xs` $+ 256 l$. The
second form converts only the first prime and reduces the results further
modulo the double `p2` (with `p2inv` $= 1/p_2$).
:::

## Plan internals {#plan-internals}

The steps of plan construction, used by drivers in other modules that
build their own plans.

::: function
int _fft_small_plan_set_bound(fft_small_plan_t P, ulong c, ulong e, ulong np_max)

Chooses the least number of primes $k \le$ `np_max` whose product is at
least $c 2^e$, starting at the first prime of `P->R`, and records the
bound. Returns 0 if there is none.
:::

::: function
void _fft_small_plan_set_window(fft_small_plan_t P, ulong zl, ulong zh, ulong zn, ulong xtrunc_max, ulong min_depth)
void _fft_small_plan_set_window_nmod(fft_small_plan_t P, ulong zl, ulong zh, ulong zn, ulong xtrunc_max)

Choose the depth and truncation for the window $[z_l, z_h)$ of a
convolution of length $z_n$ with operands of truncated length at most
`xtrunc_max`, and record the window. Require $z_l < z_h \le z_n$.

::: details
With $2^d \le z_n < 2^{d+1}$: if `xtrunc_max` $\le 2^d$, $z_h \le 2^d$
and $z_n \le z_l + 2^d$, the transform has length $2^{\max(d, d_0)}$
with $d_0$ = `min_depth` and wraps around. Otherwise the depth is
$\max(d_0, \lceil \log_2 (z_n + 4) \rceil)$ and the truncation $z_n$
rounded up to a multiple of 256 (or the full length below 256). The four
extra coefficients cover the bias digits of signed integer exports;
the `nmod` version omits them and uses $d_0 = 4$. Both set `bits` to 0.
Use `min_depth` = 8 for plans whose conversions work on whole blocks.
:::
:::

::: function
int _fft_small_plan_fit_window(fft_small_plan_t P, ulong zl, ulong zh, ulong zn, ulong xtrunc_max, ulong min_depth)

Like `_fft_small_plan_set_window` with the depth of `P` fixed. Returns 0
if the depth is below `min_depth` or the window needs a longer transform.

::: details
For depths below 8 the truncation is the full length $2^d$ and the
function accepts every window, including windows of convolutions longer
than $2^d$ that do not wrap below $z_l$; callers must check $z_n \le 2^d$
or the wrap-around condition themselves.
:::
:::

::: macro
_fft_small_len_trunc(x)
_fft_small_op_trunc(x, P)

The truncation of an operand or output of length $x$: $x$ rounded up to a
multiple of 256 or, below one block, the full transform length.
`_fft_small_len_trunc` anticipates it before a plan exists, using the
length $2^{\max(4, \lceil \log_2 x \rceil)}$ below 256;
`_fft_small_op_trunc` gives the one the plan `P` imposes, $2^d$ for a
depth $d < 8$.
:::

::: function
void _fft_small_plan_set_normalizers(fft_small_plan_t P)

Computes the stride and the normalisation factors `m[i]` from the primes
and the depth, which must be set.
:::

## Convolution engine {#engine}

The engine shared by the polynomial drivers of this and other modules.
Stage 1 converts and transforms the operands for each prime, multiplies
pointwise with the normalisation folded in, and transforms back; it runs
in parallel over the primes. Stage 2 reconstructs ranges of output
coefficients by the Chinese remainder theorem, in parallel over the
ranges. The coefficient type enters through callbacks.

::: type
fft_small_mod_func
fft_small_crt_func
fft_small_s2_finish_func

The callbacks: conversion of an operand into the double buffer of one
prime; reconstruction of the output coefficients
$[\mathit{start}, \mathit{stop})$ (relative to $z_l$) from the residues;
and an optional serial pass after stage 2 that receives the range
boundaries and the per-range values that the reconstruction left in its
`local` slots.
:::

::: macro
FFT_SMALL_CRT_LOCAL

The number of words (4) of per-range output, such as carries, a
reconstruction callback can leave for the finishing pass.
:::

::: type
fft_small_crt_range_struct

One stage 2 range: its end `stop_zi` and its `local` words.
:::

::: type
fft_small_conv_tuning

Three threshold predicates for a coefficient type: the number of threads
for stage 1, whether each stage 1 worker should transform the second
operand on an extra thread, and whether to request up to eight threads
for stage 2 (may be `NULL`).
:::

::: type
fft_small_conv_arg_struct

The arguments of one convolution: both operands with their lengths,
auxiliary values and truncations (or a precomputed transform `bfft` of
the second operand with its stride), the squaring flag, the `transposed`
flag, the size used by the tuning predicates, the bound `min_an_bn` on the
products per output coefficient, the callbacks with their `params`, and
the tuning table.

::: details
With `transposed` set, stage 1 computes the transpose, with respect to the
first operand, of the multiplication by the second: the first operand is
converted at the full truncation `P->ztrunc` and sent through
`sd_ifft_trunc_tr`, the pointwise product is unchanged, and
`sd_fft_trunc_tr` with output truncation `atrunc` replaces
`sd_ifft_trunc`. Squaring is not allowed in this mode. With `transposed`
equal to 2, the images of the first operand through `sd_ifft_trunc_tr`
are given in `bfft` (the first operand is ignored) and the second operand
is converted and transformed in the working buffer, so that the
convolution costs two transforms. `_nmod_poly_mul_tr_mpn_ctx` and the
transposed products of `radix_mulmid_fft_small` set the flag to 1,
`_nmod_poly_mul_tr_precomp` sets it to 2, and all other drivers to 0.
:::
:::

::: function
void _fft_small_conv(void * z, const fft_small_plan_t P, const fft_small_conv_arg_struct * A)

Computes the window of the plan of the convolution described by `A` into
`z`, using the scratch buffer of the plan's context.
:::

## Two-prime pipeline {#two-prime}

The integer multiplication used up to 40000 limbs. The chunks are below
both primes, so conversion is plain bit slicing, and reconstruction is
one Garner step [@Gar1959] with the normalisation folded into its
constants. See [](algorithms/fft/multimodular#two-prime).

::: macro
FLINT_MPN_MUL_NP2_MAX_AN

The largest operand length, 40000 limbs, for which `mpn_ctx_mpn_mul` and
`fft_small_plan_init_mpn` prefer the two-prime pipeline.
:::

::: function
void _mpn_ctx_mpn_mul_np2(mpn_ctx_t R, nn_ptr z, nn_srcptr a, ulong an, nn_srcptr b, ulong bn)
void _mpn_ctx_mpn_mullow_np2(mpn_ctx_t R, nn_ptr z, ulong zh, nn_srcptr a, ulong an, nn_srcptr b, ulong bn)
void _mpn_ctx_mpn_mul_window_np2(mpn_ctx_t R, nn_ptr z, ulong zl, ulong zh, nn_srcptr a, ulong an, nn_srcptr b, ulong bn)

The full product, its low $z_h$ limbs, and its limbs $[z_l, z_h)$, by
the two-prime pipeline. Require $a_n \ge b_n \ge 1$ and
$z_l < z_h \le a_n + b_n$. The window is exact for $z_l = 0$ and
otherwise smaller than the exact window by 0 or 1 units of its lowest
limb.

::: details
Requires a chunk size satisfying the bound of `_fft_small_np2_choose`
with `len_bound` = 1, which holds up to well beyond 40000 limbs (with
40-bit chunks, for $b_n$ up to about 300000 limbs).
`_mpn_ctx_mpn_mul_window_np2` computes a window with fewer chunk
coefficients than $a$ as a transposed product of the reversed chunks of
$b$ and the chunks of $a$ from a limb boundary on, when its transform is
strictly shorter than that of the full product (which it does not wrap
around), except for squarings; the result is the same, bit for bit.
:::
:::

::: function
int _fft_small_np2_choose(const mpn_ctx_struct * R, ulong an, ulong bn, ulong len_bound, ulong * bits, ulong * depth)

Chooses the least depth $d \in [8, 30]$ and then the largest chunk size
$b \in [40, 45]$ with $\lceil 64 a_n / b \rceil + \lceil 64 b_n / b
\rceil - 1 \le 2^d$ and `len_bound` $\cdot \lceil 64 b_n / b \rceil \cdot
2^{2b} < p_0 p_1$. Returns 0 if there is none.
:::

::: function
void _fft_small_np2_pack(double * d, double * d2, nn_srcptr x, ulong xn, ulong len, ulong itrunc, ulong bits)

Writes the first `len` chunks of `bits` $\le 45$ bits of $(x, x_n)$ as
doubles to both `d` and `d2`, followed by zeros up to `itrunc`.
:::

::: function
void _fft_small_np2_crt_recompose(const mpn_ctx_struct * R, nn_ptr z, ulong zl, ulong zh, double * d0, double * d1, ulong nchunks, ulong navail8, ulong bits, ulong depth)

Writes to `z` the limbs $[z_l, z_h)$ of $\sum_j v_j 2^{j \cdot
\mathtt{bits}}$ over `nchunks` chunk values $v_j$ given by their residues
modulo $p_0$ in `d0` and modulo $p_1$ in `d1`, which it overwrites.
Unless `depth` is `UWORD_MAX`, the residues are raw inverse transform
outputs and are multiplied by $2^{-\mathtt{depth}}$. `navail8` is the
number of entries available to the vectorised Garner step.
:::

## Word and limb helpers {#helpers}

Small inline helpers defined in this header.

::: function
ulong n_pow2(int k)
ulong n_min(ulong a, ulong b)
ulong n_max(ulong a, ulong b)
slong z_min(slong a, slong b)
slong z_max(slong a, slong b)

Return $2^k$, and the minimum and maximum of two unsigned and of two
signed words.
:::

::: function
ulong n_cdiv(ulong a, ulong b)
ulong n_round_up(ulong a, ulong b)
ulong n_round_down(ulong a, ulong b)

Return $\lceil a/b \rceil$, $b \lceil a/b \rceil$ and $b \lfloor a/b
\rfloor$. The first two compute $a + b - 1$ and are wrong if it
overflows.
:::

::: function
ulong n_nbits(ulong x)
ulong n_nbits_nz(ulong x)
ulong n_leading_zeros(ulong x)
ulong n_trailing_zeros(ulong x)
ulong n_next_pow2m1(ulong a)

The bit length of $x$ (0 for $x = 0$; `_nz` requires $x \ne 0$), the
numbers of leading and trailing zero bits (`FLINT_BITS` for $x = 0$), and
$2^{k} - 1$ for $a$ of bit length $k$.
:::

::: function
ulong n_clog2(ulong x)
ulong n_flog2(ulong x)

Return $\lceil \log_2 x \rceil$ and $\lfloor \log_2 x \rfloor$ for
$x \ge 1$; both are 0 for $x = 0$.
:::

::: function
ulong flint_mpn_nbits(const ulong * a, ulong an)
int flint_mpn_cmp_ui_2exp(const ulong * a, ulong an, ulong b, ulong e)

The bit length of $(a, a_n)$, which need not be normalised, and the sign
of $a - b 2^e$.
:::

::: function
unsigned char flint_mpn_add_inplace_c(ulong * z, ulong zn, ulong * a, ulong an, unsigned char cf)

Adds $(a, a_n)$ and the carry `cf` to $(z, z_n)$ in place, propagating
the carry through at most $z_n$ limbs, and returns the carry out.
Requires $z_n \ge a_n \ge 1$.
:::
