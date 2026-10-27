---
title: mpn_extras.h
header: mpn_extras.h
description: Arithmetic on natural numbers stored as arrays of limbs, extending and replacing GMP's mpn layer.
category: Integers
order: 5
---

This module provides arithmetic on natural numbers stored as arrays of
limbs, in the representation of GMP's `mpn` functions. Part of it fills
gaps in GMP's public interface (modular arithmetic with precomputed
inverses, multimodular reduction and Chinese remaindering, high, low and
middle products, Hensel division); most of the rest replaces GMP's
functions by faster ones: multiplication with hand-written assembly
basecases and the small-prime FFT of `fft_small`, division, exact
division, square roots and powering. It is the layer on which `fmpz`,
`mpn_mod`, `fft_small`, `radix` and the real and ball arithmetic of
`mp_real` and `arb` are built.

The functions are grouped by task. Most groups have one *dispatching*
function, which chooses the algorithm from the operand sizes, followed by
the individual algorithms, which are exposed for testing and tuning and
whose names usually begin with an underscore. The central operations
have pages of their own: [](mpn_extras/mul), [](mpn_extras/mulhigh),
[](mpn_extras/tdiv_qr) and [](mpn_extras/multi_crt).

## Conventions {#conventions}

Throughout, $B = 2^{\bits}$ is the limb radix. The notation $(x, n)$ (or
`{x, n}`) stands for the natural number
$x_0 + x_1 B + \cdots + x_{n-1} B^{n-1}$ held in the $n$ limbs
$x_0, \ldots, x_{n-1}$ starting at the pointer `x`; within formulas we
often write $x$ for this number. The types `mp_ptr`, `mp_srcptr`,
`nn_ptr` and `nn_srcptr` are all pointers to limbs (`ulong`), the last
two being FLINT's names for them.

- An array is *normalised* if its top limb is nonzero (or its length is
  zero). Unless stated otherwise inputs need not be normalised, but many
  functions require a nonzero top limb of a divisor or a modulus.
- A divisor or modulus $d$ of $n$ limbs is *normalised* in the sense of
  division if the most significant bit of $d_{n-1}$ is set. Functions
  that need this take $d 2^{\mathrm{norm}}$ together with the shift
  `norm`.
- Output arrays must not overlap input arrays unless stated otherwise.
  When aliasing is allowed it means *exact* aliasing (the same pointer).
- None of the requirements is checked, except by `FLINT_ASSERT` in builds
  configured with `--enable-assert`.
- Functions with `n` in their name (as in `flint_mpn_mul_n`) take operands
  of the same length $n$.

## Basic operations {#basic}

::: function
void flint_mpn_zero(mp_ptr xp, mp_size_t n)
void flint_mpn_store(mp_ptr xp, mp_size_t n, mp_limb_t y)

Set the $n$ limbs of `xp` to zero, respectively to $y$. These are inline
loops, cheaper than calls to `memset` or `mpn_zero` for the short arrays
typical of FLINT.
:::

::: function
void flint_mpn_copyi(mp_ptr xp, mp_srcptr yp, mp_size_t n)
void flint_mpn_copyd(mp_ptr xp, mp_srcptr yp, mp_size_t n)

Copy the $n$ limbs of `yp` to `xp`, in increasing respectively decreasing
order of addresses. The arrays may overlap if `xp` $\le$ `yp` (for
`flint_mpn_copyi`) or `xp` $\ge$ `yp` (for `flint_mpn_copyd`).
:::

::: function
int flint_mpn_zero_p(mp_srcptr x, mp_size_t xsize)

Returns 1 if all limbs of $(x, \mathit{xsize})$ are zero, and 0 otherwise.
:::

::: function
int flint_mpn_equal_p(mp_srcptr x, mp_srcptr y, mp_size_t xsize)

Returns 1 if the arrays $(x, \mathit{xsize})$ and $(y, \mathit{xsize})$
are equal limb by limb, and 0 otherwise.
:::

::: macro
MPN_NORM(a, an)

Decreases the length variable `an` until `an` is zero or `a[an - 1]` is
nonzero.
:::

::: macro
MPN_SWAP(a, an, b, bn)

Swaps the pointers `a`, `b` (of type `mp_ptr`) and the lengths `an`, `bn`.
:::

::: macro
BITS_TO_LIMBS(b)

The number of limbs $\lceil b / \bits \rceil$ needed for $b$ bits.
:::

## Conversion and debugging {#conversion}

::: function
char * flint_mpn_get_str(char * res, int base, mp_srcptr x, mp_size_t xn, int negative)

Returns the representation of $(x, \mathit{xn})$, or of its negation if
`negative` is 1, in the given base, which must be one supported by GMP
($2 \le$ `base` $\le 62$). If `res` is `NULL`, a new string is
allocated with `flint_malloc`; otherwise `res` must have room for all the
digits, one extra digit, a minus sign and the terminating null.

::: details
In base 10 and with `fft_small`, the conversion is done with the `radix`
module (`radix_set_mpn` into the radix $10^{19}$, then
`radix_get_str_decimal`) except for $40 < \mathit{xn} < 1400$, where
GMP's `mpz_get_str` is faster; other bases and builds without
`fft_small` use `mpz_get_str`. See [](algorithms/integers/radix-conversion).
:::
:::

::: function
double flint_mpn_get_d(mp_srcptr ptr, mp_size_t size, mp_size_t sign, long exp)

Returns the double nearest towards zero to $\pm (\mathit{ptr}, \mathit{size})
\cdot 2^{\mathit{exp}}$, where the sign is that of `sign`. Requires
`size` $\ge 0$ and a nonzero top limb if `size` $> 0$. Overflow gives
an infinity of the right sign and underflow (below the denormal range)
gives zero. This is a copy of GMP's internal `mpn_get_d`.
:::

::: function
void flint_mpn_debug(mp_srcptr x, mp_size_t xsize)

Prints $(x, \mathit{xsize})$ in decimal to `stdout`, followed by the
binary digits of each limb, least significant bit first.
:::

## Addition and subtraction {#addition}

The basic carry-propagating additions are GMP's `mpn_add_n`, `mpn_sub_n`
and so on. FLINT adds `mpn_add_nc` and `mpn_sub_nc` (with a carry in)
where GMP does not export them, and `flint_mpn_sumdiff_n` (documented with
the [`fft` module](fft)) for simultaneous sum and difference.

::: function
mp_limb_t flint_mpn_2add_n_inplace(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)

Sets $(\mathit{rp}, n)$ to $r + a + b$ modulo $B^n$ and returns the carry
(0, 1 or 2). Requires $n \ge 4$. Only available on x86-64 with the ADX
instructions, where it is an assembly routine running two carry chains in
parallel; `FLINT_HAVE_NATIVE_mpn_2add_n_inplace` is then defined.
:::

::: function
int flint_mpn_signed_sub_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t n)

Sets $(\mathit{res}, n)$ to $|x - y|$ and returns 0 if $x \ge y$ and 1 if
$x < y$. The output may alias either input.
:::

::: function
int flint_mpn_signed_sub_2(mp_ptr res, mp_srcptr x, mp_srcptr y)
int flint_mpn_signed_sub_3(mp_ptr res, mp_srcptr x, mp_srcptr y)
int flint_mpn_signed_sub_4(mp_ptr res, mp_srcptr x, mp_srcptr y)
int flint_mpn_signed_sub_5(mp_ptr res, mp_srcptr x, mp_srcptr y)
int flint_mpn_signed_sub_6(mp_ptr res, mp_srcptr x, mp_srcptr y)
int flint_mpn_signed_sub_7(mp_ptr res, mp_srcptr x, mp_srcptr y)
int flint_mpn_signed_sub_8(mp_ptr res, mp_srcptr x, mp_srcptr y)

`flint_mpn_signed_sub_n` for a fixed length from 2 to 8, with the
subtraction done by straight-line carry chains (the `NN_SUB` macros).
:::

::: function
void flint_mpn_signed_div2(mp_ptr res, mp_srcptr x, mp_size_t n)

Sets $(\mathit{res}, n)$ to $x / 2$ rounded towards $-\infty$, where $x$ is
read as an $n$-limb two's complement integer (an arithmetic shift right by
one bit). The output may alias the input.
:::

::: function
mp_size_t flint_mpn_fmms1(mp_ptr y, mp_limb_t a1, mp_srcptr x1, mp_limb_t a2, mp_srcptr x2, mp_size_t n)

Tries to compute $y = a_1 x_1 - a_2 x_2$, where $x_1$ and $x_2$ have $n \ge 1$
limbs (not necessarily normalised). If the result is nonnegative and fits
in $n$ limbs, it is written to $(y, n)$ and its normalised length is
returned; otherwise $-1$ is returned and $(y, n)$ holds garbage. $y$ may
alias $x_1$ but not $x_2$.
:::

## Modular addition {#addmod}

Additions modulo an $n$-limb modulus $m$; the operands must already be
reduced modulo $m$. These are inline functions used by `mpn_mod` and
`fmpz_mod`.

::: function
void flint_mpn_negmod_n(mp_ptr res, mp_srcptr x, mp_srcptr m, mp_size_t n)
void flint_mpn_addmod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)
void flint_mpn_submod_n(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m, mp_size_t n)

Set $(\mathit{res}, n)$ to $-x$, $x + y$ and $x - y$ modulo
$(m, n)$, for $x$ and $y$ reduced modulo $m$. The output may alias the
inputs.
:::

::: function
void flint_mpn_addmod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)
void flint_mpn_submod_n_m(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_size_t yn, mp_srcptr m, mp_size_t n)

Versions of `flint_mpn_addmod_n` and `flint_mpn_submod_n` where the second
operand has only `yn` $\le n$ limbs (and is reduced modulo $m$).
:::

::: function
void flint_mpn_negmod_2(mp_ptr res, mp_srcptr x, mp_srcptr m)
void flint_mpn_addmod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
void _flint_mpn_addmod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)
void flint_mpn_submod_2(mp_ptr res, mp_srcptr x, mp_srcptr y, mp_srcptr m)

Versions for two-limb moduli, written with two-limb arithmetic macros.
`_flint_mpn_addmod_2` requires $m < B^2/2$ (the top bit of `m[1]` clear),
so that the sum cannot overflow two limbs.
:::

## Multiplication {#multiplication}

[](mpn_extras/mul) explains the algorithm selection.

::: summary algorithms/integers/basecase-multiplication :::

::: function
mp_limb_t flint_mpn_mul(mp_ptr z, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)

Sets $(z, \mathit{xn} + \mathit{yn})$ to the product of
$(x, \mathit{xn})$ and $(y, \mathit{yn})$ and returns its top limb
$z_{xn + yn - 1}$. Requires $\mathit{xn} \ge \mathit{yn} \ge 1$; $z$
must not overlap the inputs.

::: details
For short operands this calls a hard-coded routine from a table
(assembly on x86-64 with ADX for $\mathit{xn} \le 16$ and on ARM64 for
$\mathit{yn} \le 15$; portable C for $\mathit{xn} \le 7$ elsewhere);
otherwise `_flint_mpn_mul` selects between `mpn_mul_1`, GMP's
`mpn_mul_basecase` (only without FLINT assembly),
`flint_mpn_mul_toom22`, GMP's `mpn_mul` and the FFT.
:::
:::

::: function
void flint_mpn_mul_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_size_t n)

Sets $(z, 2n)$ to the product of $(x, n)$ and $(y, n)$. Requires $n \ge 1$;
$z$ must not overlap the inputs.
:::

::: function
mp_limb_t flint_mpn_sqr(mp_ptr z, mp_srcptr x, mp_size_t n)

Sets $(z, 2n)$ to the square of $(x, n)$ and returns its top limb.
Requires $n \ge 1$; $z$ must not overlap $x$.
:::

::: function
void flint_mpn_mul_shared(mp_ptr z1, mp_srcptr a, mp_size_t an, mp_ptr z2, mp_srcptr b, mp_size_t bn, mp_srcptr x, mp_size_t xn)

Sets $(z_1, \mathit{an} + \mathit{xn})$ to $ax$ and
$(z_2, \mathit{bn} + \mathit{xn})$ to $bx$. All lengths must be at
least 1, in any order; $z_1$ and $z_2$ must overlap neither each other
nor any input. When both products are long enough, the transform of $x$
is computed only once ([](fft_small/mpn_shared)).
:::

::: function
int _flint_mpn_mul_shared_use_fft(mp_size_t an, mp_size_t bn, mp_size_t xn)

Returns nonzero if `flint_mpn_mul_shared` uses a shared transform for
these lengths, that is, if `fft_small` is available and
`_fft_small_mpn_mul_shared` is expected to beat two separate products.
Callers that need the looser aliasing rules of
`_fft_small_mpn_mul_shared` can call it directly when this returns
nonzero.
:::

::: function
mp_limb_t _flint_mpn_mul(mp_ptr r, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn)
void _flint_mpn_mul_n(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n)
mp_limb_t _flint_mpn_sqr(mp_ptr r, mp_srcptr x, mp_size_t n)

The non-inline part of `flint_mpn_mul`, `flint_mpn_mul_n` and
`flint_mpn_sqr`, used for sizes beyond the hard-coded routines. Same
requirements.
:::

::: function
void flint_mpn_mul_toom22(mp_ptr pp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn, mp_ptr scratch)

Karatsuba (Toom-22) multiplication, a port of GMP's `mpn_toom22_mul` that
recurses into `flint_mpn_mul`, so that FLINT's assembly basecases are
used at the leaves. Sets $(\mathit{pp}, \mathit{an} + \mathit{bn})$ to
the product, for $\mathit{an} \ge \mathit{bn}$ with $4 \mathit{an} <
5 \mathit{bn}$ (as required by GMP). `scratch` must have room for
$2 \mathit{an} + 2 \bits$ limbs, or be `NULL` for internal
allocation.
:::

::: function
void flint_mpn_mul_toom32(mp_ptr pp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn, mp_ptr scratch)

Toom-32 multiplication of unbalanced operands, a port of GMP's
`mpn_toom32_mul` (evaluation at $0, \pm 1, \infty$) built on
`flint_mpn_mul`. Requires $\mathit{bn} + 2 \le \mathit{an}$ and
$\mathit{an} + 6 \le 3 \mathit{bn}$ (GMP's conditions). `scratch` may be
`NULL`. It is not currently used by `flint_mpn_mul`.
:::

::: function
mp_limb_t flint_mpn_mul_2(mp_ptr rp, mp_srcptr ap, mp_size_t n, mp_srcptr bp)

Sets $(\mathit{rp}, n + 2)$ to the product of $(\mathit{ap}, n)$ and the
two-limb $(\mathit{bp}, 2)$ and returns the top limb. Only available on
ARM64 (`FLINT_HAVE_NATIVE_mpn_mul_2`), where it is an assembly routine
used by `flint_mpn_mul`.
:::

::: macro
FLINT_MPN_MUL_WITH_SPECIAL_CASES(z, x, xn, y, yn)

A statement that sets $(z, \mathit{xn} + \mathit{yn})$ to the product
for operands in either order, with the one- and two-limb balanced cases
inlined and squaring detected when `x == y`.
:::

::: macro
FLINT_MPN_MUL_2X1(r2, r1, r0, a1, a0, b0)
FLINT_MPN_MUL_2X2(r3, r2, r1, r0, a1, a0, b1, b0)
FLINT_MPN_MULLOW_2X2(r1, r0, a1, a0, b1, b0)
FLINT_MPN_MUL_3P2X2(r2, r1, r0, a1, a0, b1, b0)
FLINT_MPN_SQR_2X2(r3, r2, r1, r0, a1, a0)

Products of one- and two-limb operands in registers: the full
$2 \times 1$ and $2 \times 2$ products, the low two and the low three
limbs of a $2 \times 2$ product, and a $2$-limb square. Limbs are listed
most significant first.
:::

## High and low products {#short}

Given two $n$-limb integers, the *low product* is the product modulo
$B^n$, and a *high product* is an approximation of the top $n$ limbs
$\lfloor uv / B^n \rfloor$ of the full product. The low product is
uniquely defined; the high product is not, since the carry from the low
part is not computed. [](mpn_extras/mulhigh) states the
precise contracts.

::: summary algorithms/integers/short-products :::

::: function
mp_limb_t flint_mpn_mulhigh_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)

Sets $(\mathit{res}, n)$ to an approximation of
$\lfloor uv / B^n \rfloor$ and returns an approximation of limb $n - 1$
of the product. Together they form a lower bound for
$\lfloor uv / B^{n-1} \rfloor$ that is wrong by a few units of the
returned limb at most. Requires $n \ge 1$. The output may alias the
inputs only for $n$ at most `FLINT_MPN_MULHIGH_FUNC_TAB_WIDTH`.
:::

::: function
mp_limb_t flint_mpn_sqrhigh(mp_ptr res, mp_srcptr u, mp_size_t n)

The square $u^2$ with the contract of `flint_mpn_mulhigh_n`.
:::

::: type
mp_limb_pair_t

A pair of limbs, with fields `m1` and `m2`.
:::

::: function
mp_limb_pair_t flint_mpn_mulhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
mp_limb_pair_t flint_mpn_sqrhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_size_t n)

`flint_mpn_mulhigh_n` and `flint_mpn_sqrhigh` followed by normalisation:
if the top bit of the high product is clear, the $(n+1)$-limb result
(high limbs and returned limb) is shifted left by one bit. Returns the
low limb in `m1` and the shift (0 or 1) in `m2`. For normalised inputs
the product is at least $B^{2n}/4$, so the result is normalised. This is
the form needed for floating-point mantissas (`mp_real`, `arf`). The
output must not alias the inputs.
:::

::: function
mp_limb_t flint_mpn_mullow_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)

Sets $(\mathit{res}, n)$ to $uv \bmod B^n$ and returns limb $n$ of the
product, that is, $\lfloor uv / B^n \rfloor \bmod B$ (this extra limb
may be removed from the interface in the future). Requires $n \ge 1$;
`res` must not alias `u`.
:::

::: function
void flint_mpn_mul_or_mullow_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)

Writes $uv \bmod B^{n+1}$ to the low $n + 1$ limbs of `res`. The output
must have room for $2n$ limbs: the high limbs may be used as scratch
space, or the whole product may be written when this is fastest.
:::

::: function
void flint_mpn_mul_or_mulhigh_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)

Writes a high product (with the error of `flint_mpn_mulhigh_n`) as
$n + 1$ limbs to `res + n - 1`; the low $n - 1$ limbs of the $2n$-limb
output may be used as scratch space, or the full product may be
written.
:::

::: function
mp_limb_t _flint_mpn_mulhigh_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
mp_limb_t _flint_mpn_mulhigh_basecase(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
mp_limb_t _flint_mpn_mulhigh_n_recursive(mp_ptr r, mp_srcptr x, mp_srcptr y, mp_size_t n)
mp_limb_t _flint_mpn_mulhigh_n_mulders(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
mp_limb_t _flint_mpn_mulhigh_n_mul(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
mp_limb_t _flint_mpn_mulhigh_n_naive(mp_ptr rp, mp_srcptr up, mp_srcptr vp, mp_size_t n)

The algorithms behind `flint_mpn_mulhigh_n`, with its interface.

::: details
- `_flint_mpn_mulhigh_n` is the dispatcher used beyond the hard-coded
  sizes: the basecase up to `FLINT_MPN_MULHIGH_MULDERS_CUTOFF` (40)
  limbs, Mulders' algorithm up to `FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF`
  (with `fft_small`) or `FLINT_MPN_MULHIGH_MUL_CUTOFF` (2000), and
  above, a middle product with `fft_small` respectively a full product.
- `_flint_mpn_mulhigh_basecase` is the quadratic algorithm (an assembly
  routine on x86-64 with ADX, which requires $n \ge 6$, and on ARM64,
  which requires $n > 8$). On x86-64 with ADX the basecase used by the
  dispatcher is `_flint_mpn_mulhigh_n_recursive` up to 22 limbs: it
  splits off a block handled by the fastest hard-coded routine
  (`FLINT_MPN_MULHIGH_BEST_TAB_N` limbs) and recurses.
- `_flint_mpn_mulhigh_n_mulders` computes a rough high product with one
  extra limb using `_flint_mpn_mulhigh_n_mulders_recursive`.
- `_flint_mpn_mulhigh_n_mul` computes the full product; its result is the
  exact truncation.
- `_flint_mpn_mulhigh_n_naive` is a portable reference implementation
  of the basecase, against which the assembly routines are tested
  (they return identical results).
:::
:::

::: function
void _flint_mpn_mulhigh_n_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)
void _flint_mpn_sqrhigh_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_size_t n)

Rough high product by Mulders' recursive algorithm, as described by Harvey
and Zimmermann [@HZ2011]. Writes to $\mathit{rp}_n, \ldots,
\mathit{rp}_{2n-1}$ an approximation of the top $n$ limbs of the product
(respectively the square), which is a lower bound with an error of less
than $n$ units of $\mathit{rp}_n$. `rp` must have room for $2n$ limbs;
the low $n$ limbs are used as scratch space. The split points come from
a tuned table `flint_mpn_mulhigh_k_tab`.
:::

::: function
mp_limb_t _flint_mpn_sqrhigh(mp_ptr res, mp_srcptr u, mp_size_t n)
mp_limb_t _flint_mpn_sqrhigh_basecase(mp_ptr res, mp_srcptr u, mp_size_t n)
mp_limb_t _flint_mpn_sqrhigh_mulders(mp_ptr res, mp_srcptr u, mp_size_t n)
mp_limb_t _flint_mpn_sqrhigh_sqr(mp_ptr res, mp_srcptr u, mp_size_t n)

The algorithms behind `flint_mpn_sqrhigh`: dispatcher, quadratic basecase
(on x86-64 with ADX an assembly routine requiring $n \ge 8$; elsewhere
the multiplication basecase), Mulders' algorithm (from
`FLINT_MPN_SQRHIGH_MULDERS_CUTOFF`, 90 limbs) and a full square (from
`FLINT_MPN_SQRHIGH_SQR_CUTOFF`, 2000 limbs, or a middle product with
`fft_small` from `FLINT_MPN_SQRHIGH_FFT_SMALL_CUTOFF`).
:::

::: function
mp_limb_pair_t _flint_mpn_mulhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
mp_limb_pair_t _flint_mpn_sqrhigh_normalised(mp_ptr rp, mp_srcptr xp, mp_size_t n)

The non-inline part of `flint_mpn_mulhigh_normalised` and
`flint_mpn_sqrhigh_normalised`, used beyond the hard-coded sizes.
:::

::: function
mp_limb_t _flint_mpn_mullow_n(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
mp_limb_t flint_mpn_mullow_basecase(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
mp_limb_t _flint_mpn_mullow_n_mulders(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
mp_limb_t _flint_mpn_mullow_n_mul(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
void _flint_mpn_mullow_n_mulders_recursive(mp_ptr rp, mp_srcptr np, mp_srcptr mp, mp_size_t n)

The algorithms behind `flint_mpn_mullow_n`: dispatcher, quadratic
basecase (up to `FLINT_MPN_MULLOW_MULDERS_CUTOFF`, 50 limbs; an assembly
routine on x86-64 with ADX), Mulders' recursive short product (up to
`FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF` or `FLINT_MPN_MULHIGH_MUL_CUTOFF`),
and a full product (or an `fft_small` middle product). The low product is
exact in every case, and so is the returned limb $n$.
`_flint_mpn_mullow_n_mulders_recursive` writes the $n$ low limbs only and
uses `rp` beyond them as scratch space (it needs room for $2n$ limbs).
:::

::: function
mp_limb_t _flint_mpn_mulhigh_n_fft_small(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)
mp_limb_t _flint_mpn_mullow_n_fft_small(mp_ptr res, mp_srcptr u, mp_srcptr v, mp_size_t n)

High and low products computed by `flint_mpn_mulmid_fft_small`. Only
available with `fft_small`.
:::

## Middle products {#middle}

The *windowed middle product* extracts a chosen range of limbs of a
product. For $\mathit{an}, \mathit{bn} \ge 1$ and
$0 \le \mathit{zlo} < \mathit{zhi} \le \mathit{an} + \mathit{bn}$, it
writes $\mathit{zhi} - \mathit{zlo}$ limbs to $z$ approximating the limbs
$[\mathit{zlo}, \mathit{zhi})$ of $ab$, that is,
$\lfloor ab / B^{\mathit{zlo}} \rfloor \bmod B^{\mathit{zhi} - \mathit{zlo}}$.
It is a *lower approximation*: the partial products $a_p b_q$ with
$p + q < \mathit{zlo}$ may be dropped, so the result never exceeds the
exact window, and the deficit is the carry from below limb $\mathit{zlo}$,
less than $\min(\mathit{an}, \mathit{bn}, \mathit{zlo}) \cdot B$ units of
the lowest output limb (so only the two lowest output limbs can be
affected, unless a carry propagates). With $\mathit{zlo} = 0$ the window
is exact. See
[](mpn_extras/mulhigh) and [](algorithms/integers/short-products).

::: function
void flint_mpn_mulmid(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)

Sets $(z, \mathit{zhi} - \mathit{zlo})$ to the window
$[\mathit{zlo}, \mathit{zhi})$ of $ab$ with the contract above,
choosing the fastest of the methods below for the shape. The operands
may be given in either order; squaring is detected. $z$ must not overlap
the inputs.
:::

::: function
void flint_mpn_mulmid_classical(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
void flint_mpn_mulmid_via_mul(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
void flint_mpn_mulmid_via_mullow_n(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
void flint_mpn_mulmid_via_mulhigh_n(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
void flint_mpn_mulmid_via_n_padded(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)
void flint_mpn_mulmid_fft_small(mp_ptr z, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t zlo, mp_size_t zhi)

The methods behind `flint_mpn_mulmid`, with its interface and contract:
row-wise schoolbook multiplication in time
$O((\mathit{zhi} - \mathit{zlo}) \min(\mathit{an}, \mathit{bn}))$; a
full product (`flint_mpn_mul`); a balanced low product
(`flint_mpn_mullow_n`); a balanced high product
(`flint_mpn_mulhigh_n`); GMP's balanced middle product
(`flint_mpn_mulmid_n`); and the small-prime FFT, which only
reconstructs the coefficients that reach the window. Each accepts any
shape (padding internally) but is economical only in its own regime.
`flint_mpn_mulmid_via_n_padded` exists only if
`FLINT_HAVE_NATIVE_mpn_mulmid_n` is set, and `flint_mpn_mulmid_fft_small`
only with `fft_small`.

::: details
`flint_mpn_mulmid` computes full products with `flint_mpn_mul`. Otherwise
it truncates both operands to $\mathit{zhi}$ limbs and uses: for
$\min(\mathit{an}, \mathit{bn}, \mathit{zhi} - \mathit{zlo}) < 48$, the
low product if $\mathit{zlo} = 0$, $\mathit{zhi} \ge 8$ and both operands
have at least $7 \mathit{zhi}/8$ limbs, and the schoolbook method
otherwise; else `fft_small` if both operands reach the `fft_small`
multiplication threshold; else the low product when the window starts in
the bottom eighth of the product, the high product when it ends in the
top eighth (both if the operands are about as long as the window); else
GMP's middle product for the shape $2m \times m$ with window
$[m, 2m)$ (within 12.5%); and a full product otherwise.

`flint_mpn_mulmid_fft_small` (`_mpn_ctx_mpn_mul_range`) computes a window
shorter than the longer operand as a transposed product, whose cost is
that of a product of the shorter operand and the window, when its
transform is strictly shorter than that of the window of the full
product (not for squarings); the result is the same, bit for bit. For
the window $[n, 2n)$ of a $2n \times n$ product this is a transform of
about $2n$ instead of $3n$ limbs' worth of coefficients. See
[](algorithms/integers/short-products#transposed-fft).
:::
:::

::: function
void flint_mpn_mulmid_n(mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)

Balanced middle product of $(\mathit{ap}, 2n - 1)$ and
$(\mathit{bp}, n)$: writes the $n + 2$ limbs of
$\sum_{n-1 \le i + j \le 2n - 2} a_i b_j B^{i+j-n+1}$, that is, the
diagonals $n - 1$ to $2n - 2$ of the product. Compared with the
corresponding limbs of the full product only the carry from the lower
diagonals is missing, so the low two limbs serve as guard limbs. This is
GMP's
internal `mpn_mulmid_n`, and is only declared if
`FLINT_HAVE_NATIVE_mpn_mulmid_n` is set (when the build may call GMP
internals and GMP exports the symbol).
:::

::: function
void _flint_mpn_mulhigh_known_low(mp_ptr out, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn, mp_srcptr kl, mp_size_t kl_len, mp_size_t klo, mp_size_t khi, mp_ptr scratch)

Sets $(\mathit{out}, \mathit{khi} - \mathit{klo})$ to the exact limbs
$[\mathit{klo}, \mathit{khi})$ of the product $xy$, given that its low
$\mathit{klo}$ limbs are known and equal to $(\mathit{kl}, \mathit{kl\_len})$
(zero above `kl_len`). The window is computed from a windowed middle
product starting three guard limbs below $\mathit{klo}$, whose deficit is
less than $B^2$ units of its lowest limb, and the known low limbs tell
which carry was lost. `scratch` must have room for
$\mathit{khi}$ limbs. This is used by Hensel division and the 2-adic
square root, where the low part of a product is known in advance.
:::

## Complex multiplication {#complex}

Multiplication of Gaussian integers $a_r + a_i i$ whose real and imaginary
parts are given as magnitudes (limb arrays) with separate sign flags
(0 for nonnegative, 1 for negative). These are the kernels of `fmpzi`
multiplication and of the complex arithmetic of `mp_real`. The outputs
must not overlap each other or the inputs.

::: function
void flint_mpn_mul_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
void flint_mpn_sqr_complex(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)

Sets $z_r + z_i i = (a_r + a_i i)(b_r + b_i i)$, respectively
$(a_r + a_i i)^2$. Each part has its own length, at least 1 limb, not
necessarily normalised. A *signed length* is written for each output:
the magnitude occupies `|*zr_len|` limbs and a negative value means that
the part is negative. Limbs beyond the signed length are not
significant. Each output must have room for
$\max(\mathit{arn}, \mathit{ain}) + \max(\mathit{brn}, \mathit{bin}) + 1$
limbs ($2 \max(\mathit{arn}, \mathit{ain}) + 1$ for the square).

::: details
The algorithm is chosen from the shape: the schoolbook method (four real
products, or two squares and a product) when a part is much shorter than
its partner, Karatsuba's three-product formula when the parts of each
operand are of similar length, and a transformed method on `fft_small`
transforms for large balanced operands (from
`flint_mpn_mul_complex_fft_cutoff` and
`flint_mpn_sqr_complex_fft_cutoff` limbs). The two operands need not be
of similar size.
:::
:::

::: function
void flint_mpn_mul_complex_classical(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
void flint_mpn_mul_complex_karatsuba(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
int flint_mpn_mul_complex_fft_small(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
void flint_mpn_sqr_complex_classical(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)
void flint_mpn_sqr_complex_karatsuba(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)
int flint_mpn_sqr_complex_fft_small(nn_ptr zr, slong * zr_len, nn_ptr zi, slong * zi_len, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)

The individual algorithms behind the general functions, with the same
interface; all accept any shape. The `fft_small` versions return 0,
leaving the outputs untouched, when the method is unavailable or
refuses the operands, and 1 on success.
:::

::: function
void flint_mpn_mulhigh_n_complex(nn_ptr zr, int * zr_sgn, nn_ptr zi, int * zi_sgn, nn_srcptr ar, int ar_sgn, nn_srcptr ai, int ai_sgn, nn_srcptr br, int br_sgn, nn_srcptr bi, int bi_sgn, mp_size_t n)
void flint_mpn_sqrhigh_n_complex(nn_ptr zr, int * zr_sgn, nn_ptr zi, int * zi_sgn, nn_srcptr ar, int ar_sgn, nn_srcptr ai, int ai_sgn, mp_size_t n)

High products: all parts have $n$ limbs, and each output part receives
exactly $n + 1$ limbs, zero padded, approximating the limbs $[n, 2n]$
of the exact part, together with its sign.

::: details
Relative to the exact value, the error is below
$2 + 3(n + 4)/B$ units of the lowest returned limb for the product and
below $2 + 2(n + 4)/B$ for the square, hence below 3 units for any
practical $n$: each underlying `flint_mpn_mulhigh_n`, read as $n$ limbs,
errs by an amount in $(-1 - \varepsilon, \varepsilon)$ units with
$\varepsilon = (n + 4)/B$, and each output part combines at most three of
them. The transformed path stays within $(-1.5, 0.5)$ units.
:::
:::

::: function
void flint_mpn_mullow_complex(nn_ptr zr, nn_ptr zi, mp_size_t n, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
void flint_mpn_sqrlow_complex(nn_ptr zr, nn_ptr zi, mp_size_t n, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)

Low products: `zr` and `zi` receive exactly $n \ge 1$ limbs, the real and
imaginary parts of the product (respectively the square) modulo $B^n$ as
two's complement residues. No sign is returned: a negative part $x$ is
represented by $B^n + x$ when $|x| < B^n$. The inputs are as for
`flint_mpn_mul_complex`, and only their low $n$ limbs are read. $n$ may
exceed the length of the product, the residues then being the
sign-extended exact parts.

::: details
These are what a Kronecker unpacking of signed fields below limb $n$
needs. The algorithms are those of the full products with low products
in place of full ones: four or three real low products from
`flint_mpn_mulmid` with window start 0, combined modulo $B^n$, or above
the cutoffs `_fft_small_mpn_mulmid_complex` with the window $[0, n)$ on
the inputs truncated to $n$ limbs, which reconstructs only the CRT
coefficients below limb $n$. The transforms keep the length of the
product of the truncated inputs, since a cyclic wrap would land in the
low limbs; the saving is in the reconstruction.
:::
:::

::: function
void flint_mpn_mullow_complex_classical(nn_ptr zr, nn_ptr zi, mp_size_t n, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
void flint_mpn_mullow_complex_karatsuba(nn_ptr zr, nn_ptr zi, mp_size_t n, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
int flint_mpn_mullow_complex_fft_small(nn_ptr zr, nn_ptr zi, mp_size_t n, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
void flint_mpn_sqrlow_complex_classical(nn_ptr zr, nn_ptr zi, mp_size_t n, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)
void flint_mpn_sqrlow_complex_karatsuba(nn_ptr zr, nn_ptr zi, mp_size_t n, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)
int flint_mpn_sqrlow_complex_fft_small(nn_ptr zr, nn_ptr zi, mp_size_t n, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)

The individual algorithms behind the low products, with the conventions
of the full products' individual algorithms.
:::

::: function
void flint_mpn_mulmid_complex(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
void flint_mpn_sqrmid_complex(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)

Windowed products: `zr` and `zi` receive $\mathit{hi} - \mathit{lo}$
limbs each, $0 \le \mathit{lo} < \mathit{hi}$, the limbs
$[\mathit{lo}, \mathit{hi})$ of the real and imaginary parts in two's
complement: for each exact part $x$ the window is
$\lfloor x / B^{\mathit{lo}} \rfloor + e \bmod B^{\mathit{hi} - \mathit{lo}}$
with a small error $e$. For $\mathit{lo} = 0$ these are the low products
($e = 0$); on the transformed path $|e| \le 1$; otherwise
$|e| < 4 (m + 1) B$ with $m = \min(\max(\mathit{arn}, \mathit{ain}),
\max(\mathit{brn}, \mathit{bin}))$. Inputs are as for the low products
(only their low $\mathit{hi}$ limbs are read).

::: details
The two's complement window of a signed part is linear, so the real
windows of the schoolbook and Karatsuba schemes combine modulo
$B^{\mathit{hi} - \mathit{lo}}$, and it is what a Kronecker unpacking of
signed fields reads. To extract the polynomial coefficients
$[n_{lo}, n_{hi})$ of $b$-bit fields, take
$\mathit{hi} = \lceil n_{hi} b / \bits \rceil$ and $\mathit{lo}$ a few
guard limbs below bit $n_{lo} b - 1$; the unpacked fields are exact
provided that $|e| B^{\mathit{lo}}$ is below the distance from the low part
$\sum_{k < n_{lo}} r_k 2^{k b}$ to $\pm 2^{n_{lo} b - 1}$, which a
bound $|r_k| \le 2^{b - 1} - 2^{b - s}$ makes at least
$2^{n_{lo} b - s}$ (see `_fmpzi_poly_mulmid_KS`, which uses three
guard limbs).

The schoolbook and Karatsuba variants combine real windows from
`flint_mpn_mulmid`. The transformed variant splits the integers into
chunks and multiplies them as chunk sequences on `fft_small` transforms
(four forward and two inverse transforms; two and two for the square),
reconstructing only the coefficients that can reach limbs
$[\mathit{lo}, \mathit{hi})$. When a power of two $2^d$ admits it
(operands and the coefficients below the top of the window fitting in
$2^d$, and the coefficients beyond $2^d$ aliasing below the first one
reaching the window) the transforms are cyclic of length $2^d$. For the
shape of the Newton middle product ($2n$ by $n$ limbs, window $[n, 2n)$)
this gives a transform of about $2n$ chunks instead of $3n$ when $2n$
chunks lie a little below a power of two.
:::
:::

::: function
void flint_mpn_mulmid_complex_classical(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
void flint_mpn_mulmid_complex_karatsuba(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
int flint_mpn_mulmid_complex_fft_small(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn, nn_srcptr br, mp_size_t brn, int br_sgn, nn_srcptr bi, mp_size_t bin, int bi_sgn)
void flint_mpn_sqrmid_complex_classical(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)
void flint_mpn_sqrmid_complex_karatsuba(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)
int flint_mpn_sqrmid_complex_fft_small(nn_ptr zr, nn_ptr zi, mp_size_t lo, mp_size_t hi, nn_srcptr ar, mp_size_t arn, int ar_sgn, nn_srcptr ai, mp_size_t ain, int ai_sgn)

The individual algorithms behind the windowed products, with the error
bounds stated there.
:::

::: variable
slong flint_mpn_mul_complex_fft_cutoff
slong flint_mpn_sqr_complex_fft_cutoff
slong flint_mpn_mullow_complex_fft_cutoff
slong flint_mpn_sqrlow_complex_fft_cutoff
slong flint_mpn_mulmid_complex_fft_cutoff
slong flint_mpn_sqrmid_complex_fft_cutoff

Sizes in limbs from which the general complex functions use the
transformed path. They are machine dependent and writable, for tuning.
:::

::: variable
slong flint_mpn_mul_complex_sparse_fft_cutoff
slong flint_mpn_mulmid_complex_sparse_fft_cutoff

Sizes in limbs (of the shorter operand) from which the full and low
products, respectively the windows with $\mathit{lo} > 0$, of operands
with a zero part (a real or purely imaginary operand) use the transformed
path. That path transforms only the nonzero parts, so that a real
operand is transformed once for both output parts (three forward and two
inverse transforms for a real times a complex operand, against four and
two for two separate real products).
:::

## Packed vectors and polynomials {#packed}

The following routines work on arrays of integers that all occupy the same
number of limbs, as when the coefficients of a polynomial are packed into
a uniform number of limbs. Results are computed modulo $B^s$ for a fixed
number $s$ of output limbs, which makes them exact when the true results
(as unsigned integers, or as two's complement integers in the signed
versions) are known to fit in $s$ limbs. They are the kernels of
`fmpz_poly` and `fmpz_mat` multiplication for small coefficients.

::: type
flint_mpn_dot_func_t
flint_mpn_dot_strided_func_t

Function pointer types `void (*)(nn_ptr res, nn_srcptr a, nn_srcptr b, slong len)`
and `void (*)(nn_ptr res, nn_srcptr a, slong astride, nn_srcptr b, slong bstride, slong len)`
for the dot product kernels.
:::

::: variable
const flint_mpn_dot_strided_func_t flint_mpn_dot_strided_tab
const flint_mpn_dot_func_t flint_mpn_dot_tab
const flint_mpn_dot_func_t flint_mpn_dot_rev_tab

Tables of dot product kernels. The kernel at index
`[sgn][n1][n2][s - (n1 + n2 - 1)]` sets $(\mathit{res}, s)$ to
$\sum_{k=0}^{\mathit{len}-1} a_k b_k \bmod B^s$, where the $a_k$ have
$n_1$ limbs and the $b_k$ have $n_2$ limbs, unsigned for `sgn` = 0 and
two's complement for `sgn` = 1, for $1 \le n_2 \le n_1$ and
$s \in \{n_1 + n_2 - 1, n_1 + n_2, n_1 + n_2 + 1\}$. Entries are `NULL`
where no kernel exists (in particular for $n_1 < n_2$).

::: details
The table dimensions are `[2][FLINT_MPN_DOT_TAB_N + 1][FLINT_MPN_DOT_TAB_N + 1][3]`
for the strided kernels and
`[2][FLINT_MPN_DOT_DEDICATED_TAB_N + 1][FLINT_MPN_DOT_DEDICATED_TAB_N + 1][3]`
for the others, with `FLINT_MPN_DOT_TAB_N` = 4 and
`FLINT_MPN_DOT_DEDICATED_TAB_N` = 2. In `flint_mpn_dot_tab` the entries
are contiguous, $a_k = a[k n_1]$, $b_k = b[k n_2]$; in
`flint_mpn_dot_rev_tab` the second array is read backwards,
$b_k$ being the entry $\mathit{len} - 1 - k$; in the strided kernels the
entries are `astride` and `bstride` limbs apart (the strides may be
negative). The kernels (static functions such as `_flint_mpn_dot_2x2_5`)
are generated by `dev/gen_mpn_dot_rev.py`.
:::
:::

::: function
void _flint_mpn_dot_rev_generic(nn_ptr res, nn_srcptr a, slong n1, nn_srcptr b, slong n2, slong len, slong s, nn_ptr scratch)
void _flint_mpn_dot_rev_generic_sdiff(nn_ptr res, nn_srcptr a, slong n1, nn_srcptr b, slong n2, slong len, int sdiff, nn_ptr t)
void _flint_mpn_dot_rev_generic_signed(nn_ptr res, nn_srcptr a, slong n1, nn_srcptr b, slong n2, slong len, slong s, nn_ptr scratch)
void _flint_mpn_dot_rev_generic_signmag(nn_ptr res, nn_srcptr a, slong n1, nn_srcptr b, slong n2, slong len, slong s, nn_ptr t)

Reversed dot products for any $n_1 \ge n_2 \ge 1$, $\mathit{len} \ge 1$
and $s \in \{n_1 + n_2 - 1, n_1 + n_2, n_1 + n_2 + 1\}$, using
`flint_mpn_mul` and `mpn_add_n`. The scratch space must have
$n_1 + n_2$ limbs. The `sdiff` version takes $s - (n_1 + n_2 - 1)$, which
should be a compile-time constant. The signed version handles the signs
with branches and is competitive only when they are predictable. In the
sign-magnitude version each coefficient is a sign limb (0 or 1) followed
by the magnitude, so that the strides are $n_1 + 1$ and $n_2 + 1$; the
output is in two's complement, and the scratch space must have
$n_1 + n_2 + 2s$ limbs.
:::

::: type
flint_mpn_poly_mulmid_func_t
flint_mpn_poly_sqrmid_func_t

Function pointer types
`void (*)(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlo, slong nhi)`
and `void (*)(nn_ptr res, nn_srcptr f, slong flen, slong nlo, slong nhi)`
for the fixed-size polynomial multiplication kernels.
:::

::: variable
const flint_mpn_poly_mulmid_func_t flint_mpn_poly_mulmid_tab
const flint_mpn_poly_sqrmid_func_t flint_mpn_poly_sqrmid_tab
const flint_mpn_poly_mulmid_func_t flint_mpn_poly_mulmid_short_tab
const flint_mpn_poly_sqrmid_func_t flint_mpn_poly_sqrmid_short_tab
const flint_mpn_poly_mulmid_func_t flint_mpn_poly_mulmid_signmag_tab
const flint_mpn_poly_sqrmid_func_t flint_mpn_poly_sqrmid_signmag_tab

Instances of `_flint_mpn_poly_mulmid_classical` (and of squaring) for
fixed coefficient sizes, indexed like the dot product tables (the
sign-magnitude tables have no sign index). The coefficient loop and the
dot product code are fused, so that no function call is made per term or
per coefficient; in the squaring kernels the doubling of the cross terms
and the addition of the square term are done in registers.

::: details
The plain tables (dimensions `[2][4][4][3]`) cover $n_1 \le 3$. The
`short` tables (dimensions `[2][5][5][3]`) cover $n_1 = 4$ with
$s \le 8$. The sign-magnitude tables (`[5][5][3]`) take coefficients
with a sign limb. The kernels are static and generated by
`dev/gen_mpn_dot_rev.py`.
:::
:::

::: function
void _flint_mpn_poly_mulmid_classical(nn_ptr res, nn_srcptr f, slong flen, slong n1, nn_srcptr g, slong glen, slong n2, slong nlo, slong nhi, slong s, int sgn)

Sets $(\mathit{res}, (\mathit{nhi} - \mathit{nlo}) s)$ to the coefficients
$[\mathit{nlo}, \mathit{nhi})$ of the product of $(f, \mathit{flen}\, n_1)$
and $(g, \mathit{glen}\, n_2)$, viewed as polynomials with $n_1$-limb
respectively $n_2$-limb coefficients, computed modulo $B^s$ by the
classical algorithm. Requires $n_1 \ge n_2$ and
$s \in \{n_1 + n_2 - 1, n_1 + n_2, n_1 + n_2 + 1\}$. The coefficients are
unsigned for `sgn` = 0, two's complement for `sgn` = 1, and
sign-magnitude for `sgn` = 2 (a sign limb followed by $n_1$ respectively
$n_2$ magnitude limbs; the output is two's complement). The results are
exact provided that the requested coefficients fit in $s$ limbs.
Squaring is detected when `f == g` and `flen == glen`.

::: details
Each output coefficient is a dot product. For $n_1 \le 3$ the function
dispatches to the fused kernels of `flint_mpn_poly_mulmid_tab`; for
$n_1 = 3, 4$ and short lengths to those of
`flint_mpn_poly_mulmid_short_tab`; otherwise to loops over the dot
product kernels or `flint_mpn_mul` with `mpn_add_n`.
:::
:::

::: function
void _flint_mpn_poly_mul_classical(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlimbs, slong slimbs, int sgn)
void _flint_mpn_poly_mul_karatsuba(nn_ptr res, nn_srcptr f, slong flen, nn_srcptr g, slong glen, slong nlimbs, slong slimbs, slong cutoff, int norm, int sgn)

Set $(\mathit{res}, (\mathit{flen} + \mathit{glen} - 1)\, \mathit{slimbs})$
to the full product of $(f, \mathit{flen}\, \mathit{nlimbs})$ and
$(g, \mathit{glen}\, \mathit{nlimbs})$, viewed as polynomials with
`nlimbs`-limb coefficients (unsigned for `sgn` = 0, two's complement for
`sgn` = 1), by the classical algorithm respectively Karatsuba's. The
output coefficients are computed modulo $B^{\mathit{slimbs}}$, where
`slimbs` is one of $2\,\mathit{nlimbs} - 1$, $2\,\mathit{nlimbs}$,
$2\,\mathit{nlimbs} + 1$. Squaring is detected when `f == g` and
`flen == glen`.

::: details
The results are exact provided that all coefficients of the product
fit in `slimbs` limbs. For Karatsuba this must also hold for all
intermediate coefficients (the products of sums of the parts, which grow
by up to two bits per level of recursion), and `norm` must bound the
number of leading unused bits (excluding the sign bit when `sgn` is 1)
in the top limb of the input coefficients. When `norm` reaches zero in the
recursion, the coefficients of the sums are extended to
$\mathit{nlimbs} + 1$ limbs, so `slimbs` must then be at least
$2\,\mathit{nlimbs} + 1$. The recursion switches to the classical
algorithm when the shorter length is below `cutoff`; the classical
instances are resolved once per level, so that the base cases are called
directly. Unbalanced lengths are handled by splitting the longer operand
into pieces of the length of the shorter one at every level.
:::
:::

::: function
void _flint_mpn_poly_mulmid_karatsuba(nn_ptr res, nn_srcptr f, slong flen, slong nlimbs1, slong norm1, nn_srcptr g, slong glen, slong nlimbs2, slong norm2, slong nlo, slong nhi, slong slimbs, slong cutoff, int sgn)

Sets $(\mathit{res}, (\mathit{nhi} - \mathit{nlo})\, \mathit{slimbs})$ to
the coefficients $[\mathit{nlo}, \mathit{nhi})$ of the product of
$(f, \mathit{flen}\, \mathit{nlimbs1})$ and
$(g, \mathit{glen}\, \mathit{nlimbs2})$, with the conventions of
`_flint_mpn_poly_mul_karatsuba` except that the operands may have
different numbers of limbs, with headrooms `norm1` and `norm2`. `slimbs`
must be at least $\mathit{nlimbs1} + \mathit{nlimbs2} - 1$, and at least
$\mathit{nlimbs1} + \mathit{nlimbs2} + 1$ if either operand may get
extended in the recursion (when its headroom is smaller than the
recursion depth plus one). Squaring is detected when `f == g`,
`flen == glen` and `nlimbs1 == nlimbs2`.

::: details
The coefficients of the inputs that cannot contribute to the window are
discarded first, and the algorithm is chosen according to the shape:

- a full product if the window is the whole product;
- Mulders' short product for a low product ($\mathit{nlo} = 0$): a full
  product of the low parts, split at about $0.7 n$, plus two recursive
  short products (one when squaring); a high product is computed as the
  reversal of a low product of the reversed inputs;
- the Karatsuba middle product of Hanrot, Quercia and Zimmermann
  [@HQZ2004] for windows in which all coefficients of the shorter operand
  contribute (for example the $n$ coefficients $[n-1, 2n-1)$ of a
  $(2n-1) \times n$ product), computed in blocks of the length of the
  shorter operand, with the edges computed classically;
- otherwise a full product of the trimmed inputs clipped to the window,
  or classical multiplication for small windows.

The short products are used from about 4 times `cutoff`. The Karatsuba
middle product uses two's complement arithmetic internally (it involves
differences of parts of the longer operand), so `slimbs` must include a
sign bit even for unsigned inputs; it consumes one bit of headroom on
every branch of the recursion, so when the headroom is smaller than the
recursion depth the coefficients are extended early and a larger cutoff
is used.
:::
:::

::: macro
FLINT_MPN_POLY_MUL_UNSIGNED
FLINT_MPN_POLY_MUL_SIGNED
FLINT_MPN_POLY_MUL_SIGNMAG
FLINT_MPN_POLY_MUL_BIAS

The representations of packed coefficients accepted by
`_flint_mpn_poly_mulmid` (values 0 to 3): nonnegative integers
(zero-extended); two's complement integers; sign-magnitude, a sign limb
followed by the magnitude (classical algorithm only); and biased, $x +
2^b$ where $b$ is the bit bound of the operand, which the multiplication
treats as unsigned, leaving it to the caller to correct the outputs using
sliding window sums of the inputs. The bias correction costs
$O(\mathit{slimbs})$ per output coefficient but nothing per term, so it
wins for long dot products and with Karatsuba.
:::

::: type
mpn_poly_mul_params_struct
mpn_poly_mul_params_t

Packing parameters for `_flint_mpn_poly_mulmid`: the representation
`method` (one of the macros above), the numbers of limbs `nlimbs1` and
`nlimbs2` per coefficient of the two operands (magnitude limbs in the
sign-magnitude case), their headroom `norm1` and `norm2` (leading unused
bits, not counting the sign bit), and the number of output limbs
`slimbs`.
:::

::: function
void _flint_mpn_poly_mulmid_params(mpn_poly_mul_params_t P, slong len1, slong bits1, slong len2, slong bits2, slong nlo, slong nhi, int squaring)

Chooses packing parameters that are correct and (as far as the tuning
goes) optimal for computing the coefficients $[\mathit{nlo}, \mathit{nhi})$
of the product of polynomials of lengths `len1` and `len2` with
`_flint_mpn_poly_mulmid`, given bounds in the format of
`_fmpz_vec_max_bits`: the coefficients have at most $|\mathit{bits1}|$ and
$|\mathit{bits2}|$ bits, and a negative bound indicates that negative
coefficients may be present. `squaring` should be set if the operands
are the same.

::: details
For nonnegative coefficients the representation is always unsigned, with
the minimal number of limbs. Otherwise the choice between two's
complement, biased and sign-magnitude coefficients depends on the limb
counts, the lengths and whether Karatsuba will be used: the two's
complement kernels are somewhat slower than the unsigned ones and only
exist in the classical algorithm for small limb counts, the bias method
has an overhead per output coefficient, and when the sign bit would cost
a limb (coefficients of exactly $k \bits$ bits) the classical
sign-magnitude kernels remain the fastest option up to many times the
Karatsuba cutoff.
:::
:::

::: function
slong _flint_mpn_poly_mulmid_slimbs(slong minlen, slong nlimbs1, slong norm1, slong nlimbs2, slong norm2, int method, slong cutoff)

Returns the number of output limbs required by `_flint_mpn_poly_mulmid`
for the given packed operands when the dot products have at most `minlen`
terms and the Karatsuba cutoff is `cutoff`: the bound for the classical
algorithm, or, if Karatsuba is used, the bound including a sign bit and
the growth in the recursion, and $\mathit{nlimbs1} + \mathit{nlimbs2} + 1$
if an operand may get extended. Passing `WORD_MAX` as the cutoff gives
the classical bound.
:::

::: function
slong _flint_mpn_poly_mulmid_cutoff(slong nlimbs1, slong nlimbs2, int squaring)
slong _flint_mpn_poly_karatsuba_cutoff(slong n1, slong n2, int squaring)

The Karatsuba cutoff (shorter length) used by `_flint_mpn_poly_mulmid` for
coefficients of $n_1$ and $n_2$ limbs. `_flint_mpn_poly_karatsuba_cutoff`
is the tuned value, from 32 for one-limb products down to 4 for large
coefficients (64 down to 8 for squaring); `_flint_mpn_poly_mulmid_cutoff`
returns it unless overridden with the global
`_flint_mpn_poly_mulmid_force_cutoff`.
:::

::: function
int _flint_mpn_poly_karatsuba_may_extend(slong norm, slong minlen, slong cutoff)
int _flint_mpn_poly_mulmid_use_karatsuba(slong minlen, slong norm, int method, slong cutoff)

Whether the Karatsuba recursion may need to extend an operand with the
given headroom by a limb, and whether `_flint_mpn_poly_mulmid` uses
Karatsuba for dot products of `minlen` terms: never for sign-magnitude
coefficients, and with tight headroom only from 1.5 times the cutoff.
:::

::: variable
int _flint_mpn_poly_mulmid_force_method
slong _flint_mpn_poly_mulmid_force_cutoff

Overrides of the representation chosen by `_flint_mpn_poly_mulmid_params`
and of the Karatsuba cutoff, for tests and profiling ($-1$ means
automatic).
:::

::: function
void _flint_mpn_poly_mulmid(nn_ptr res, nn_srcptr f, slong flen, slong nlimbs1, slong norm1, nn_srcptr g, slong glen, slong nlimbs2, slong norm2, slong nlo, slong nhi, slong slimbs, int method)
void _flint_mpn_poly_mul(nn_ptr res, nn_srcptr f, slong flen, slong nlimbs1, slong norm1, nn_srcptr g, slong glen, slong nlimbs2, slong norm2, slong slimbs, int method)

Set $(\mathit{res}, (\mathit{nhi} - \mathit{nlo})\, \mathit{slimbs})$ to
the coefficients $[\mathit{nlo}, \mathit{nhi})$ of the product
(respectively the full product) of $(f, \mathit{flen})$ and
$(g, \mathit{glen})$, packed as described by `method` (a biased operand
is treated as unsigned), choosing between
`_flint_mpn_poly_mulmid_classical` (always for sign-magnitude
coefficients) and `_flint_mpn_poly_mulmid_karatsuba`, after discarding
the input coefficients that cannot contribute to the window. `slimbs`
must be at least `_flint_mpn_poly_mulmid_slimbs` (as computed by
`_flint_mpn_poly_mulmid_params`). Squaring is detected when `f == g`,
`flen == glen` and `nlimbs1 == nlimbs2`.
:::

## Division {#division}

Euclidean division and its relatives. Unless stated otherwise the dividend
$(a, \mathit{an})$ and the divisor $(b, \mathit{bn})$ satisfy
$\mathit{an} \ge \mathit{bn} \ge 1$ and $b_{bn-1} \ne 0$ (the top limb of
$a$ may be zero), and the quotient has $\mathit{an} - \mathit{bn} + 1$
limbs. [](mpn_extras/tdiv_qr) describes the algorithm
selection.

::: summary algorithms/integers/division :::

::: function
void flint_mpn_tdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
void flint_mpn_tdiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
void flint_mpn_tdiv_r(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

Euclidean division: sets $(q, \mathit{an} - \mathit{bn} + 1)$ to
$\lfloor a/b \rfloor$ and $(r, \mathit{bn})$ to $a - qb$. The quotient
must not overlap the inputs. As with GMP's `mpn_tdiv_qr`, $r$ may be equal
to $a$ (the remainder then overwrites the low $\mathit{bn}$ limbs of $a$)
but must not overlap $b$. These replace GMP's `mpn_tdiv_qr` and are
faster for all sizes.
:::

::: function
void flint_mpn_cdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
void flint_mpn_cdiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
void flint_mpn_cdiv_r(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

Ceiling division: sets $(q, \mathit{an} - \mathit{bn} + 2)$ to
$\lceil a/b \rceil$ and $(r, \mathit{bn})$ to $qb - a \in [0, b)$. The
quotient can overflow $\mathit{an} - \mathit{bn} + 1$ limbs (for example
for $a = B^2 - 1$, $b = B$), so one extra limb, which is 0 or 1, is
written. Aliasing as for `flint_mpn_tdiv_qr`.
:::

::: function
int flint_mpn_ndiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
void flint_mpn_ndiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
int flint_mpn_ndiv_r(mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

Division rounding to nearest, ties to even: sets
$(q, \mathit{an} - \mathit{bn} + 2)$ to the integer nearest to $a/b$ and
$(r, \mathit{bn})$ to $|a - qb| \le b/2$, and returns the sign ($-1$, 0
or 1) of $a - qb$. As for ceiling division, one extra quotient limb
(0 or 1) is written.
:::

::: function
int flint_mpn_div(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

If $b$ divides $a$, sets $(q, \mathit{an} - \mathit{bn} + 1)$ to the
quotient and returns 1; otherwise returns 0, the contents of $q$ being
undefined. This is a version of `flint_mpn_divides` based on
`flint_mpn_tdiv_qr`.
:::

::: function
int flint_mpn_divides(mp_ptr q, mp_srcptr array1, mp_size_t limbs1, mp_srcptr arrayg, mp_size_t limbsg, mp_ptr temp)

If $(\mathit{arrayg}, \mathit{limbsg})$ divides
$(\mathit{array1}, \mathit{limbs1})$, sets
$(q, \mathit{limbs1} - \mathit{limbsg} + 1)$ to the quotient and returns
1; otherwise returns 0. Requires $\mathit{limbs1} \ge \mathit{limbsg} > 0$
and a nonzero top limb of the divisor. `temp` must have room for
`limbsg` limbs.
:::

::: function
void flint_mpn_divapprox(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

Approximate quotient: sets $(q, \mathit{an} - \mathit{bn} + 1)$ to
$\lfloor a/b \rfloor$ or $\lfloor a/b \rfloor + 1$, like GMP's internal
`mpn_divappr_q`. No aliasing is permitted. This follows the algorithm
selection of `flint_mpn_tdiv_q`, so that it is never slower, but skips the
final corrections.

::: details
It uses `_mp_real_div_newton` with one guard limb, rounded up by more
than its error bound, instead of the Newton division, and does so already
instead of the divide and conquer or short division when both
$\mathit{bn}$ and $\mathit{an} - \mathit{bn} + 1$ are at least
`FLINT_MPN_DIVAPPROX_NEWTON_CUTOFF`.
:::
:::

::: function
void flint_mpn_divapprox_fraction(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t f)

Sets $(q, \mathit{an} - \mathit{bn} + 1 + f)$ to $\lfloor a B^f / b \rfloor$
or $\lfloor a B^f / b \rfloor + 1$, the approximate quotient with $f$
fraction limbs. A negative $f$ (with $\mathit{an} - \mathit{bn} + 1 + f
\ge 1$) gives the top limbs of the quotient.

::: details
The algorithm is that of `flint_mpn_divapprox` on the zero padded
numerator, so that short divisors with many fraction limbs cost
$O(qn \cdot \mathit{bn})$; the padded numerator is only formed for the
register-based, schoolbook, blockwise and divide and conquer algorithms,
not for the short and Newton divisions.
:::
:::

::: function
void flint_mpn_inv(mp_ptr q, mp_srcptr x, mp_size_t xn, mp_size_t n)

Correctly truncated reciprocal: sets $(q, n - \mathit{xn} + 2)$ to
$\lfloor B^n / x \rfloor$. Requires $x_{xn-1} \ne 0$ and
$n \ge \mathit{xn}$. The top output limb is nonzero only when $x$ is a
power of $B$.

::: details
One- and two-limb $x$ use a chain of 2/1 or 3/2 divisions of zero limbs.
Otherwise, when $\min(\mathit{xn}, n - \mathit{xn} + 2)$ is below
`FLINT_MPN_INV_NEWTON_CUTOFF` (quotients of at most $\mathit{xn} + 2$
limbs), `FLINT_MPN_INV_NEWTON_LONG_CUTOFF` (quotients of at most
$6 \mathit{xn}$ limbs) or `FLINT_MPN_INV_NEWTON_VERYLONG_CUTOFF` (longer
quotients), the quotient is computed by `flint_mpn_tdiv_q` (reading only
the top $O(n - \mathit{xn})$ limbs of $x$ for short quotients), and
above by `_mp_real_inv_newton` with three extra fraction limbs and the
certification and correction scheme of the Newton division, which avoids
forming the numerator $B^n$ in the common case.
:::
:::

::: function
void flint_mpn_invapprox(mp_ptr q, mp_srcptr x, mp_size_t xn, mp_size_t n)

Approximate reciprocal with the interface of `flint_mpn_inv`: sets
$(q, n - \mathit{xn} + 2)$ to $\lfloor B^n / x \rfloor$ or
$\lfloor B^n / x \rfloor + 1$. It uses the algorithms of `flint_mpn_inv`
(with the same cutoffs) without the correction steps, and above the
Newton cutoffs `_mp_real_inv_newton` with one guard limb instead of
three, so that it is never slower.
:::

::: function
void flint_mpn_divexact(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

Exact division: sets $(q, \mathit{an} - \mathit{bn} + 1)$ to $a/b$,
assuming that $b$ divides $a$ (otherwise the result is meaningless).

::: details
Single-limb divisors are handled by `mpn_divexact_1`. Otherwise the
division is done 2-adically: writing $b = 2^v B^k b'$ with $b'$ odd, the
quotient is $(a / (2^v B^k))\, b'^{-1} \bmod B^m$, where $m$ is the
number of quotient limbs ($\mathit{an} - \mathit{bn} + 1$, or one less
when the top limb of $a$ is below that of $b$). With
$d = \min(\mathit{bn}, m)$ significant divisor limbs, $d \le 7$ uses
register-resident code specialised for each $d$, with the shift done on
the fly and the limb inverse $b'^{-1} \bmod B$ computed inline. Larger
$d$ use ports of GMP's schoolbook and divide and conquer Hensel division
(`mpn_sbpi1_bdiv_q`, and `mpn_dcpi1_bdiv_q` from
`FLINT_MPN_DC_BDIV_Q_CUTOFF` limbs) with FLINT's multiplication, and the
Newton iteration `flint_mpn_bdiv_q` from $d \ge$
`FLINT_MPN_DIVEXACT_NEWTON_CUTOFF`, or already from $\mathit{bn} \ge$
`FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF` when $m \ge 4 \mathit{bn}$.
:::
:::

::: type
flint_mpn_divexact_preinv_struct
flint_mpn_divexact_preinv_t

Precomputed data for repeated exact division by a fixed divisor $b$: its
odd part $b'$ (with $b = 2^v B^k b'$), $v$ and $k$, and the 2-adic inverse
$b'^{-1} \bmod B^m$, where $m$ is one more than the number of limbs of
$b'$ unless a smaller precision was requested with
`flint_mpn_divexact_preinv_init_prec`.
:::

::: function
void flint_mpn_divexact_preinv_init(flint_mpn_divexact_preinv_t pre, mp_srcptr b, mp_size_t bn)
void flint_mpn_divexact_preinv_init_prec(flint_mpn_divexact_preinv_t pre, mp_srcptr b, mp_size_t bn, mp_size_t qn)
void flint_mpn_divexact_preinv_clear(flint_mpn_divexact_preinv_t pre)

Initialise `pre` for exact division by $(b, \mathit{bn})$ (with
$b_{bn-1} \ne 0$), and free it. The `init_prec` version computes the
inverse only to $\min(\mathit{bn}' + 1, \mathit{qn})$ limbs, which is
much cheaper when all quotients are known to have at most $\mathit{qn}$
limbs; longer quotients remain correct but are not accelerated.
:::

::: function
void flint_mpn_divexact_preinv(mp_ptr q, mp_srcptr a, mp_size_t an, const flint_mpn_divexact_preinv_t pre)

Exact division by the fixed divisor of `pre`: sets
$(q, \mathit{an} - \mathit{bn} + 1)$ to $a/b$, assuming that $b$ divides
$a$. Quotients of at most $\mathit{bn}' + 1$ limbs cost a single
`flint_mpn_mullow_n`; longer ones use block Hensel division with the
stored inverse (one low and one high product per block), except for
divisors $b'$ of at most 7 limbs, where the register-based code of
`flint_mpn_divexact` is faster. Only the low limbs of $a$ that can
influence the quotient are read.
:::

::: function
int flint_mpn_divisible(mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

Returns 1 if $b$ divides $a$ and 0 otherwise. Requires $\mathit{bn} \ge 1$
and $b_{bn-1} \ne 0$; $a$ may have zero top limbs and $\mathit{an}$ may be
zero or less than $\mathit{bn}$. Used by `fmpz_divisible`.

::: details
After the trivial cases and the 2-adic part ($b = 2^v B^k b'$ with $b'$
odd), one-limb $b'$ use the Hensel remainder `mpn_modexact_1_odd`. When
both the divisor and the quotient have at least 16 limbs, the operands
are first screened modulo $2^{48} - 1$ (which has nine small prime
factors; `flint_mpn_mod_2exp48m1`): a prime dividing $b$ but not $a$ is
revealed for about 60% of random pairs. Then, with
$m = \mathit{an} - \mathit{bn} + 1$: $b'$ of at most
`FLINT_MPN_DIVEXACT_SMALL_BN` limbs use `_flint_mpn_divisible_small`;
huge operands ($\min(m, \mathit{bn}) \ge$
`FLINT_MPN_DIVEXACT_NEWTON_CUTOFF`, or $\mathit{bn} \ge$
`FLINT_MPN_DIVEXACT_UNBALANCED_CUTOFF` and $m \ge 8 \mathit{bn}$) take
the candidate quotient from the Newton Hensel division and form only the
high part of its product with $b$ (`_flint_mpn_mulhigh_known_low`), the
low part being known; otherwise, as in GMP's `mpn_divisible_p`, a Hensel
division with remainder (ports of `mpn_sbpi1_bdiv_qr` and
`mpn_dcpi1_bdiv_qr`, the latter from `FLINT_MPN_DC_BDIV_QR_CUTOFF` limbs)
is done after padding $a$ so that its top limb is below that of $b'$, and
$b'$ divides $a$ if and only if the remainder equals $b'$.
:::
:::

::: function
int flint_mpn_divisible_1_odd(mp_srcptr x, mp_size_t xsize, mp_limb_t d)

Returns 1 if $(x, \mathit{xsize})$ is divisible by the odd limb $d \ge 3$
and 0 otherwise, using GMP's `mpn_modexact_1_odd` where available.
:::

::: function
mp_limb_t flint_mpn_mod_2exp48m1(mp_srcptr a, mp_size_t n)

Returns a value below $2^{50}$ congruent to $(a, n)$ modulo $2^{48} - 1$,
the quantity computed by GMP's internal `mpn_mod_34lsub1` (which is called
for inputs shorter than 48 limbs when available). Longer inputs use an
AVX2 kernel when compiled with AVX2 support, or a portable carry-chain
loop. Only available for 64-bit
limbs.
:::

## Division with precomputed inverses {#preinv}

Division by a normalised divisor whose inverse has been computed once, in
the style of [](algorithms/arithmetic/moller-granlund) for one-limb
divisors and its 3/2 analogue for longer ones. These functions are used
by `fmpz_mod` and `mpn_mod` for reductions modulo a fixed modulus.

::: function
mp_limb_t flint_mpn_preinv1(mp_limb_t d1, mp_limb_t d0)

Returns the 3/2 inverse $\lfloor (B^3 - 1) / (d_1 B + d_0) \rfloor - B$
of the leading two limbs $d_1, d_0$ of a divisor, for use by the `preinv1`
functions. Requires $d_1 \ge B/2$. It is computed inline from the limb
inverse of $d_1$, as GMP's `invert_pi1`.
:::

::: function
mp_limb_t flint_mpn_divrem_preinv1(mp_ptr q, mp_ptr a, mp_size_t m, mp_srcptr b, mp_size_t n, mp_limb_t dinv)

Divides $(a, m)$ by the normalised $(b, n)$, returning the top limb of
the quotient (0 or 1), writing the remaining $m - n$ quotient limbs to $q$
and the remainder in place to $(a, n)$. Requires $m \ge n \ge 2$ and
`dinv` = `flint_mpn_preinv1(b[n - 1], b[n - 2])`. This is
`_flint_mpn_divrem_preinv1` (schoolbook or divide and conquer division)
with scratch space allocated as needed.
:::

::: function
mp_limb_t flint_mpn_divrem_1_preinv(mp_ptr qp, mp_srcptr up, mp_size_t n, mp_limb_t d, mp_limb_t dinv, unsigned int norm)

Divides $(\mathit{up}, n)$ by the limb $d$, writing the quotient to
$(\mathit{qp}, n)$ and returning the remainder. Requires $n \ge 1$,
$d \ne 0$, `norm` the number of leading zero bits of $d$ and `dinv` =
`n_preinvert_limb(d)`. `qp` may alias `up`. This is equivalent to
`mpn_divrem_1(qp, 0, up, n, d)` and faster for small $n$; for large $n$
GMP's assembly is usually faster.
:::

::: function
mp_limb_t flint_mpn_divrem_2_1_preinv_norm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv)
mp_limb_t flint_mpn_divrem_2_1_preinv_unnorm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv, unsigned int norm)
mp_limb_t flint_mpn_divrem_3_1_preinv_norm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv)
mp_limb_t flint_mpn_divrem_3_1_preinv_unnorm(mp_ptr qp, mp_srcptr up, mp_limb_t d, mp_limb_t dinv, unsigned int norm)

`flint_mpn_divrem_1_preinv` for dividends of exactly 2 and 3 limbs. The
`norm` versions require a normalised divisor, the `unnorm` versions an
unnormalised one ($\mathit{norm} > 0$).
:::

::: function
void flint_mpn_preinvn(mp_ptr dinv, mp_srcptr d, mp_size_t n)

Sets $(\mathit{dinv}, n)$ to the precomputed inverse
$\lfloor B^{2n} / (d + 1) \rfloor - B^n$ of the normalised $n$-limb
integer $d$ (zero if $d = B^n - 1$), for use by the `preinvn`
functions ([multiword modular multiplication](algorithms/arithmetic/multiword-mulmod#inverse)).
:::

::: function
void flint_mpn_mod_preinvn(mp_ptr r, mp_srcptr a, mp_size_t m, mp_srcptr d, mp_size_t n, mp_srcptr dinv)

Sets $(r, n)$ to $a \bmod d$, where $(a, m)$ has $m \ge n$ limbs and the
normalised $(d, n)$ has the inverse `dinv` from `flint_mpn_preinvn`. No
aliasing of $r$ with $a$ is permitted. This is not always as fast as
ordinary division.
:::

::: function
mp_limb_t flint_mpn_divrem_preinvn(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t m, mp_srcptr d, mp_size_t n, mp_srcptr dinv)

Sets $(q, m - n)$ to the quotient and $(r, n)$ to the remainder of
$(a, m)$ by the normalised $(d, n)$ with inverse `dinv`, returning the
top quotient limb (0 or 1). Requires $m \ge n$; no aliasing is permitted.
This is used by `flint_mpn_tdiv_qr` for short divisors with long
dividends.
:::

::: macro
FLINT_MPN_UDIV_QR_2BY1(q, r, nh, nl, d, di)
FLINT_MPN_UDIV_QR_3BY2(q, r1, r0, n2, n1, n0, d1, d0, dinv)
FLINT_MPN_UDIV_QR_3BY2_HW(q, r1, r0, n2, n1, n0, d1, d0)

Möller–Granlund 2/1 division of $(n_h, n_l)$ by the normalised limb $d$
with inverse `di` = $\lfloor (B^2 - 1)/d \rfloor - B$ (requires
$n_h < d$), and 3/2 division of $(n_2, n_1, n_0)$ by the normalised
$(d_1, d_0)$ with `dinv` = `flint_mpn_preinv1(d1, d0)` (requires
$(n_2, n_1) < (d_1, d_0)$) [@GraMol2010]. The `HW` version does the 3/2
division without an inverse, by one hardware 2/1 division and at most two
corrections (a step of Knuth's algorithm D); it is only sensible where
hardware division is fast (`FLINT_PREINVERT_LIMB_USE_NATIVE`).
:::

::: macro
FLINT_MPN_INVERT_LIMB(v, d)

Sets $v$ to $\lfloor (B^2 - 1)/d \rfloor - B$ for a normalised limb $d$:
by a hardware division where that is fast, by GMP's assembly
`mpn_invert_limb` where available, and otherwise by a table-based Newton
iteration.
:::

## Hensel division and 2-adic square roots {#hensel}

Arithmetic modulo $B^n$, developing results from the least significant
limb upwards. These are ports of the `radix_invmod_bn`,
`radix_divmod_bn`, `radix_rsqrtmod_bn` and `radix_sqrtmod_bn` algorithms
of the `radix` module to the radix $B$, where the reductions modulo the
digit radix disappear.

::: function
void flint_mpn_binv(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)

Sets $(\mathit{res}, n)$ to $x^{-1} \bmod B^n$, for odd $x$. Newton
iteration on middle products; up to `FLINT_MPN_BINV_MULLOW_CUTOFF` limbs
full low products are used instead, and the cases $n \le 3$ are inlined.
:::

::: function
void flint_mpn_bdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t n)
void flint_mpn_bdiv_q(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t n)

Hensel division: for odd $(b, \mathit{bn})$, sets $(q, n)$ to the 2-adic
quotient with $qb \equiv a \pmod{B^n}$. If $r$ is not `NULL` (it is
`NULL` for `bdiv_q`), the remainder $(a - qb)/B^n \bmod B^{\mathit{bn}}$
is written to $(r, \mathit{bn})$, so that
$a \equiv qb + B^n r \pmod{B^{n + \mathit{bn}}}$; the division is exact
with an $n$-limb quotient if and only if $r = 0$ and $a$ has no nonzero
limbs at positions $\ge n + \mathit{bn}$. $q$ may alias $a$. Unlike GMP's
`mpn_bdiv_qr`, the number $n$ of quotient limbs is a parameter
independent of $\mathit{an}$ and $\mathit{bn}$.
:::

::: function
void flint_mpn_bdiv_qr_1(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_limb_t b, mp_size_t n)
void flint_mpn_bdiv_qr_classical(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t n)
void _flint_mpn_bdiv_qr_classical_preinv(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_srcptr binv, mp_size_t n)
void flint_mpn_bdiv_qr_karp_markstein(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_size_t n)

The algorithms behind `flint_mpn_bdiv_qr`, with its interface: a
one-limb divisor; the block algorithm running from the bottom up in blocks
of $\mathit{bn}$ limbs, each costing one low and one high
$\mathit{bn} \times \mathit{bn}$ product (used by the dispatcher when
$\mathit{bn} < n/2$; the `preinv` version takes
$b^{-1} \bmod B^{\mathit{bn}}$); and the Karp–Markstein form, which uses a
half-precision inverse from `flint_mpn_binv` and one refinement step
[@KM1997].
:::

::: function
int flint_mpn_brsqrt(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)
int flint_mpn_bsqrt(mp_ptr res, mp_srcptr x, mp_size_t xn, mp_size_t n)

If $x \equiv 1 \pmod 8$, set $(\mathit{res}, n)$ to the 2-adic reciprocal
square root $y$ (with $x y^2 \equiv 1 \pmod{B^n}$), respectively the
square root $s$ (with $s^2 \equiv x \pmod{B^n}$), that is congruent to 1
modulo 4, and return 1. Otherwise return 0 without writing anything (an
odd $x$ is a square modulo $B^n$ if and only if $x \equiv 1 \pmod 8$).

::: details
The reciprocal square root is computed by a Newton iteration in bit
precision (each step maps $c$ correct bits to $2c - 2$), using high
products with known low part above the basecase range and full low
products below `FLINT_MPN_BRSQRT_MULLOW_CUTOFF` limbs. The square root
adds one Karp–Markstein refinement to a reciprocal square root of
$\lceil n/2 \rceil + 1$ limbs.
:::
:::

## Square roots {#sqrt}

::: function
mp_size_t flint_mpn_sqrtrem(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)

Sets $(s, \lceil \mathit{an}/2 \rceil)$ to $\lfloor \sqrt a \rfloor$.
If $r$ is not `NULL`, the remainder $a - s^2 \le 2s$ is written to
$(r, \lceil \mathit{an}/2 \rceil + 1)$, zero padded, and its number of
limbs is returned; $r$ must have room for $\max(\mathit{an}, 2)$ limbs. If
$r$ is `NULL`, the return value is 0 if $a$ is a perfect square and 1
otherwise (GMP's convention). Requires $\mathit{an} \ge 1$ and
$a_{an-1} \ne 0$.

::: details
On 64-bit machines, inputs of up to four limbs have dedicated code.
Two-limb inputs use the hardware double-precision square root for the
initial approximation, followed for $a \ge 2^{108}$ by one Newton step
and a final adjustment. When the compiler reports a fast fused
multiply-add (`FP_FAST_FMA`), the Newton step is evaluated in floating
point: the residual $a - D^2$ of the integral approximation
$D = 2^{32} \sqrt{a_1}$ is computed exactly by an FMA from an exact
splitting of $a$, so that no 128/64-bit division is needed. Three- and
four-limb inputs use one fully unrolled step of the divide and conquer
square root of Zimmermann [@Zim1999], as GMP's `mpn_dc_sqrtrem`: the
two-limb code gives the root of the top two limbs, one 2/1 division gives
the remaining limb of the root (for three limbs the remaining 32 bits,
using an unbalanced split instead of GMP's zero padding), and the
remainder is corrected once without branching.

Larger inputs below `FLINT_MPN_SQRTREM_NEWTON_CUTOFF` limbs use
`_flint_mpn_sqrtrem_divconquer` on 64-bit machines (GMP's `mpn_sqrtrem`
on 32-bit machines): GMP's divide and conquer algorithm with the code
above as basecases, FLINT's squaring and the divide and conquer division
of `_flint_mpn_divrem_preinv1`. When no remainder is wanted, the last
step uses an approximate quotient with a guard limb
(`flint_mpn_divapprox`), and the remainder of the division is evaluated
(from the low part of a product, `flint_mpn_mulmid`) only if that is
inconclusive, in particular for perfect squares.

Above the cutoff, $a$ is viewed as a fixed-point number in $[B^{-2}, 1)$
and `_mp_real_sqrt_newton` is used with three extra fraction limbs, so
that the truncated root is certified whenever the first fraction limb
lies in $[2, B-2]$; otherwise it is corrected by $O(1)$ steps. With a
`NULL` remainder, the Newton path tests $a = s^2$ on the low two limbs
before comparing in full.
:::
:::

::: function
int flint_mpn_sqrt(mp_ptr s, mp_srcptr a, mp_size_t an)

Checked exact square root, the analogue of `flint_mpn_div`: if $a$ is a
perfect square, sets $(s, \lceil \mathit{an}/2 \rceil)$ to its square
root and returns 1; otherwise returns 0, the contents of $s$ being
undefined. Requires $a_{an-1} \ne 0$.

::: details
Nonsquares are usually rejected before any root is computed by quadratic
residue tests modulo 256 and modulo the factors 9, 5, 7, 13, 17, 97, 241
and 257 of $2^{48} - 1$, all obtained from one reduction modulo
$2^{48} - 1$ (`flint_mpn_mod_2exp48m1`); about one random nonsquare in
900 survives. On 32-bit builds the residues come from two `mpn_mod_1`
calls instead.
:::
:::

::: function
int flint_mpn_is_square(mp_srcptr a, mp_size_t an)

Returns 1 if $(a, \mathit{an})$ is a perfect square and 0 otherwise,
using the screening of `flint_mpn_sqrt` and allocating space for the root
only when the screening passes. Requires $a_{an-1} \ne 0$. Used by
`fmpz_is_square`.
:::

::: function
mp_size_t _flint_mpn_sqrtrem(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
mp_size_t _flint_mpn_sqrtrem_divconquer(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
void _flint_mpn_sqrtrem_newton(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
void _flint_mpn_sqrtrem_gmp(mp_ptr s, mp_ptr r, mp_srcptr a, mp_size_t an)
int _flint_mpn_sqrt(mp_ptr s, mp_srcptr a, mp_size_t an)
int _flint_mpn_is_square(mp_srcptr a, mp_size_t an)

The algorithms behind `flint_mpn_sqrtrem`, `flint_mpn_sqrt` and
`flint_mpn_is_square`, with their interfaces: the dispatcher, the divide
and conquer algorithm (requires $\mathit{an} \ge 5$ and needs no remainder
space beyond $\lceil \mathit{an}/2 \rceil + 1$ limbs), the Newton
algorithm, GMP's `mpn_sqrtrem`, and the checked root and square test
after the test modulo 256.
:::

## Powering {#powering}

::: function
mp_size_t flint_mpn_pow_bound_limbs(mp_srcptr x, mp_size_t xn, ulong e)

Returns an upper bound for the number of limbs of $x^e$, including one
limb of slack for the intermediate products of `flint_mpn_pow`. Requires
$x_{xn-1} \ne 0$. For small results the bound $e \cdot
\mathrm{bits}(x)$ is used; otherwise $e \log_2 x$ is evaluated in double
precision from the top two limbs and rounded up, which is essentially
tight also for small bases with large exponents.
:::

::: function
mp_size_t flint_mpn_pow(mp_ptr res, mp_srcptr x, mp_size_t xn, ulong e)

Sets $\mathit{res}$ to $x^e$ by left-to-right binary exponentiation and
returns the number of limbs of the result. `res` must have room for
`flint_mpn_pow_bound_limbs(x, xn, e)` limbs and must not overlap $x$.
Requires $\mathit{xn} \ge 1$ and $x_{xn-1} \ne 0$; $e = 0$ gives 1.

::: details
Low zero limbs of $x$ are stripped and restored as an offset in the
output; trailing zero bits are shifted out only when the base is short
(at most 8 limbs) and restored with one shift of the result. While the
power fits in one or two limbs, fixed-size arithmetic is used. Afterwards
the running power alternates between `res` and a scratch buffer without
copying: the parity of the remaining out-of-place operations is chosen so
that the final product lands in `res`, and multiplications by a
single-limb base are done in place with `mpn_mul_1`.
:::
:::

## The mpz interface {#mpz}

::: function
void flint_mpz_tdiv_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_tdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_tdiv_r(mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_fdiv_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_fdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_fdiv_r(mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_cdiv_qr(mpz_ptr q, mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_cdiv_q(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_cdiv_r(mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_mod(mpz_ptr r, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_divexact(mpz_ptr q, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_sqrtrem(mpz_ptr s, mpz_ptr r, mpz_srcptr a)
void flint_mpz_sqrt(mpz_ptr s, mpz_srcptr a)
void flint_mpz_gcd(mpz_ptr g, mpz_srcptr a, mpz_srcptr b)
void flint_mpz_gcdext(mpz_ptr g, mpz_ptr s, mpz_ptr t, mpz_srcptr a, mpz_srcptr b)
int flint_mpz_invert(mpz_ptr r, mpz_srcptr x, mpz_srcptr n)

Drop-in replacements for the GMP functions of the same names, with
identical semantics for signs, rounding, cofactors and aliasing, and the
limb-level work done by `flint_mpn_tdiv_qr`, `flint_mpn_divexact`,
`flint_mpn_sqrtrem`, `flint_mpn_gcd` and `_flint_mpn_gcdext_inplace`. They are
used by the `fmpz` division, square root, GCD and inversion functions
for large operands. In `flint_mpz_gcdext`, any of $g$, $s$, $t$ may be
`NULL` (not computed). In `flint_mpz_invert`, $n$ must be nonzero.
:::

## Modular multiplication {#mulmod}

Multiplication modulo a fixed $n$-limb modulus $d$. The functions take the
shifted modulus $d 2^{\mathrm{norm}}$ (normalised) and its precomputed
inverse. Unless stated otherwise, operands must be reduced modulo $d$ and
padded to $n$ limbs. The algorithms, including the preconditioned
methods, are described in [](algorithms/arithmetic/multiword-mulmod).

::: function
void flint_mpn_mulmod_preinv1(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_srcptr d, mp_limb_t dinv, ulong norm)

Sets $(r, n)$ to $ab \bmod d$, where `d` is the normalised modulus
$(d, n)$ with `dinv` = `flint_mpn_preinv1(d[n - 1], d[n - 2])` (for
$n = 1$, `flint_mpn_preinv1(d[0], 0)`, which equals
`n_preinvert_limb(d[0])`). If $a$, $b$ and $d$ have all been shifted left
by `norm` bits so that $d$ is normalised, the result carries the same
shift. The product is reduced by `_flint_mpn_divrem_preinv1`.
:::

::: function
void flint_mpn_mulmod_preinvn(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm)

Sets $(r, n)$ to $ab \bmod d$ for the normalised $(d, n)$ with inverse
`dinv` from `flint_mpn_preinvn`, with the shift convention of
`flint_mpn_mulmod_preinv1`.
:::

::: function
void flint_mpn_mulmod_preinvn_2(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_srcptr d, mp_srcptr dinv, ulong norm)

`flint_mpn_mulmod_preinvn` for $n = 2$, except that $a$, $b$ and the
result are *not* shifted: only `d` holds the shifted modulus.
:::

::: function
void flint_mpn_powmod_preinvn(mp_ptr res, mp_srcptr a, mp_srcptr e, mp_size_t en, mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm)

Sets $(\mathit{res}, n)$ to $a^e \bmod d$, where the exponent is the
nonnegative $(e, \mathit{en})$ (not necessarily normalised), with the
conventions of `flint_mpn_mulmod_preinvn` ($a$ and the result are
shifted, the exponent is not). Requires the unshifted modulus to be at
least 2. `res` must not alias $a$ or $e$.
:::

::: function
void flint_mpn_fmmamod_preinvn(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_srcptr e, mp_srcptr f, mp_size_t n, mp_srcptr d, mp_srcptr dinv, ulong norm)
void flint_mpn_fmmamod_preinvn_2(mp_ptr r, mp_srcptr a, mp_srcptr b, mp_srcptr e, mp_srcptr f, mp_srcptr d, mp_srcptr dinv, ulong norm)

Set $(r, n)$ to $ab + ef \bmod d$ with a single reduction, given `d`
holding $d 2^{\mathrm{norm}}$ and its inverse `dinv`; the operands and the
result are unshifted. The second version is for $n = 2$.
:::

::: function
int flint_mpn_mulmod_want_precond(mp_size_t n, slong num, ulong norm)

Given a modulus of $n$ limbs and shift `norm`, returns the best method for
`num` multiplications by a fixed operand, accounting for the cost of
the precomputation: `MPN_MULMOD_PRECOND_NONE` (use
`flint_mpn_mulmod_preinvn`), `MPN_MULMOD_PRECOND_SHOUP` (use
`flint_mpn_mulmod_precond_shoup`) or `MPN_MULMOD_PRECOND_MATRIX` (use
`flint_mpn_mulmod_precond_matrix`).
:::

::: function
void flint_mpn_mulmod_precond_shoup_precompute(mp_ptr apre, mp_srcptr a, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)
void flint_mpn_mulmod_precond_shoup(mp_ptr res, mp_srcptr a, mp_srcptr apre, mp_srcptr b, mp_size_t n, mp_srcptr d, ulong norm)

Multiplication by a fixed operand $0 \le a < d$ with Shoup's method
([](algorithms/arithmetic/shoup)): the precomputation writes the $n$-limb
quotient $\lfloor a B^n / d \rfloor$ to `apre` (the modulus given as
`dnormed` = $d 2^{\mathrm{norm}}$ with its inverse); the multiplication
sets $(\mathit{res}, n)$ to $ab \bmod d$ with a high product and a low
product, for $b$ reduced modulo $d$. Here `d` is the unshifted modulus.
:::

::: function
mp_size_t flint_mpn_mulmod_precond_matrix_alloc(mp_size_t n)
void flint_mpn_mulmod_precond_matrix_precompute(mp_ptr apre, mp_srcptr a, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)
void flint_mpn_mulmod_precond_matrix(mp_ptr rp, mp_srcptr apre, mp_srcptr b, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)
void flint_mpn_fmmamod_precond_matrix(mp_ptr rp, mp_srcptr a1pre, mp_srcptr b1, mp_srcptr a2pre, mp_srcptr b2, mp_size_t n, mp_srcptr dnormed, mp_srcptr dinv, ulong norm)

The matrix method of preconditioned multiplication
([multiword modular multiplication](algorithms/arithmetic/multiword-mulmod#matrix)): the precomputation
writes $a 2^{\mathrm{norm}} B^i \bmod d 2^{\mathrm{norm}}$ for
$0 \le i < n$ (plus one junk limb) to `apre`, which must have room for
`flint_mpn_mulmod_precond_matrix_alloc(n)` limbs; $ab \bmod d$ is then a
matrix-vector product with the limbs of $b$ followed by one short
reduction, in $O(n^2)$ operations. The `fmma` version computes
$a_1 b_1 + a_2 b_2 \bmod d$. Requires $n \ge 2$ and $b$ reduced modulo
$d$.
:::

::: function
int flint_mpn_mulmod_2expp1_basecase(mp_ptr xp, mp_srcptr yp, mp_srcptr zp, int c, flint_bitcnt_t b, mp_ptr tp)

Multiplication modulo $2^b + 1$ for the representation of the
Schönhage–Strassen FFT (see [](fft)): the inputs are $y + c_1 2^b$ and
$z + c_0 2^b$, where $(y, n)$ and $(z, n)$ with $n = \lceil b/\bits
\rceil$ are reduced modulo $2^b$ and `c` holds the top bits ($c_1 = 2$
bit, $c_0 = 1$ bit), each input being at most $2^b$. Sets $(x, n)$ and
returns the top bit of the product modulo $2^b + 1$, fully reduced. `tp`
must have room for $2n$ limbs.
:::

## Square roots modulo a prime {#sqrtmod}

The following take an $n$-limb odd modulus $d$ with $d_{n-1} \ne 0$ and an
operand already reduced to $[0, d)$. Primality of $d$ is assumed and never
checked. They are the implementation behind `fmpz_sqrtmod` and the square
roots of `mpn_mod` and `fmpz_mod`.

::: function
int flint_mpn_is_square_mod(nn_srcptr a, nn_srcptr d, mp_size_t n)

Returns 1 if $a$ is a square modulo the odd prime $d$, and 0 if it is
not. This evaluates the Jacobi symbol $(a/d)$ with GMP, in quasi-linear
time; a return value of 0 means that the symbol is $-1$, so it also
proves that $a$ is not a square for an odd $d$ that is not prime, while a
return value of 1 proves nothing in that case.
:::

::: function
int flint_mpn_sqrtmod(nn_ptr res, nn_srcptr a, nn_srcptr d, mp_size_t n)
int flint_mpn_sqrtmod_preinv(nn_ptr res, nn_srcptr a, nn_srcptr d, mp_size_t n, nn_srcptr dinv, flint_bitcnt_t norm)

Sets $(\mathit{res}, n)$ to a square root of $a$ modulo the odd prime $d$
and returns 1. Returns 0, setting `res` to zero, if the Jacobi symbol rules
a root out (which is conclusive whether or not $d$ is prime), and $-1$,
again setting `res` to zero, if the algorithm itself failed, which happens
only for $d$ not prime (or, with negligible probability, for a prime
whose least quadratic nonresidue exceeds $2^{20}$). The `preinv` version
takes the inverse `dinv` of $d 2^{\mathrm{norm}}$ from `flint_mpn_preinvn`
and $\mathrm{norm} = \mathrm{clz}(d_{n-1})$, as an `mpn_mod` or
`fmpz_mod` context holds them; `dinv` may be `NULL`.

::: details
For $d - 1 = q 2^s$ with $q$ odd, the root is a single exponentiation
when $s \le 2$ ($a^{(d+1)/4}$, or Atkin's formula); Tonelli–Shanks is used
for $3 \le s < 32$, and Müller's variant of the Cipolla–Lehmer algorithm,
whose cost does not depend on $s$, for $s \ge 32$. See
[](algorithms/number-theory/sqrt-mod-p). The arithmetic is
`flint_mpn_mulmod_preinvn` and `flint_mpn_powmod_preinvn` throughout.
:::
:::

## Greatest common divisor {#gcd}

`flint_mpn_gcd` and `flint_mpn_gcdext` compute the GCD and the extended
GCD of two nonzero values, given in either order, without modifying them.
`_flint_mpn_gcd_inplace` and `_flint_mpn_gcdext_inplace` are the same
algorithms with the interfaces of GMP's `mpn_gcd` and `mpn_gcdext`, which
restrict and destroy the operands but save copying them.
`flint_mpn_gcd_vec_1` and `flint_mpn_gcd_vec` compute the GCD of many
values. See [](algorithms/integers/gcd) for the algorithms.

::: function
mp_size_t flint_mpn_gcd(mp_ptr gp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)

Sets $(\mathit{gp}, \ell)$ to the greatest common divisor of
$(\mathit{ap}, \mathit{an})$ and $(\mathit{bp}, \mathit{bn})$ and returns
its length $\ell$. The operands must be nonzero with nonzero top limbs
($\mathit{an}, \mathit{bn} \ge 1$); they are not modified. `gp` needs
room for $\min(\mathit{an}, \mathit{bn})$ limbs; it may be equal to `ap`
or `bp` but must not otherwise overlap them.
:::

::: function
mp_size_t flint_mpn_gcdext(mp_ptr gp, mp_ptr sp, mp_size_t * sn, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)

Sets $(\mathit{gp}, \ell)$ to $g = \gcd(a, b)$, returning $\ell$, and
$(\mathit{sp}, |{*\mathit{sn}}|)$ with the sign of $*\mathit{sn}$ to the
cofactor $s$ of $a$ in $g = s a + t b$ that `mpz_gcdext` gives: normally
$|s| < b/(2g)$, and the conventions for the exceptional cases determine
$s$ uniquely. The other cofactor is $t = (g - s a)/b$. The operands must
be nonzero with nonzero top limbs; they are not modified. `gp` needs room
for $\min(\mathit{an}, \mathit{bn})$ limbs and `sp` for $\mathit{bn}$
limbs. The outputs must not overlap each other; each may be equal to
`ap` or `bp` (if that has the room) but must not otherwise overlap them.
When $\mathit{an} < \mathit{bn}$, $s$ is obtained from the cofactor of
$b$ by a multiplication and an exact division.
:::

::: macro
flint_mpn_gcd_full
flint_mpn_gcd_full2

Deprecated. `flint_mpn_gcd_full` is an alias of `flint_mpn_gcd`, and
`flint_mpn_gcd_full2(gp, ap, an, bp, bn, temp)` calls
`flint_mpn_gcd(gp, ap, an, bp, bn)`, ignoring the scratch space `temp`.
:::

::: function
mp_size_t _flint_mpn_gcd_inplace(mp_ptr gp, mp_ptr up, mp_size_t un, mp_ptr vp, mp_size_t vn)

The GCD with the interface of GMP's `mpn_gcd`: sets
$(\mathit{gp}, \ell)$ to the greatest common divisor of
$(\mathit{up}, \mathit{un})$ and $(\mathit{vp}, \mathit{vn})$ and returns
$\ell$. It is required that $\mathit{un} \ge \mathit{vn} > 0$, that the
top limb of $v$ is nonzero and that at least one of the operands is odd.
Both inputs are destroyed. `gp` needs room for $\mathit{vn}$ limbs; it
may be equal to `up` or `vp` but must not otherwise overlap them.
:::

::: function
mp_size_t _flint_mpn_gcdext_inplace(mp_ptr gp, mp_ptr sp, mp_size_t * sn, mp_ptr up, mp_size_t un, mp_ptr vp, mp_size_t vn)

The extended GCD with the interface of GMP's `mpn_gcdext`: sets
$(\mathit{gp}, \ell)$ to $g = \gcd(u, v)$, returning $\ell$, and
$(\mathit{sp}, |{*\mathit{sn}}|)$ with the sign of $*\mathit{sn}$ to the
cofactor $s$ of $u$ as for `flint_mpn_gcdext`. It is required that
$\mathit{un} \ge \mathit{vn} > 0$ and that the top limb of $v$ is
nonzero; the top limb of $u$ may be zero, so that $u$ can be padded to
the length of $v$ to obtain its cofactor directly (as for a modular
inverse). Both inputs are destroyed, and both input arrays must have room
for one more limb. `gp` needs room for $\mathit{vn}$ limbs and `sp` for
$\mathit{vn} + 1$ limbs; they must not overlap each other or the inputs.
:::

::: function
ulong flint_mpn_gcd_vec_1(mp_srcptr x, slong len, ulong g)

Returns the greatest common divisor of $g$ and the `len` limbs
$x_0, \ldots, x_{\mathit{len}-1}$, any of which (and $g$) may be zero;
the result is 0 only if all are zero.
:::

::: function
mp_size_t flint_mpn_gcd_vec(mp_ptr g, mp_size_t gn, const mp_srcptr * xp, const mp_size_t * xn, slong len)

Sets $(g, \ell)$ to the greatest common divisor of $(g, \mathit{gn})$
and the `len` values $(\mathit{xp}_i, \mathit{xn}_i)$, and returns its
length $\ell$ (0 if all are zero). The values are given by pointer, so
that they need not be contiguous or copied, and must be normalised (top
limb nonzero), with $\mathit{xn}_i = 0$ for zero; $\mathit{gn} = 0$ means
that $g$ is not used as input. The inputs are not modified. The output
array $g$ must have room for $\min(\mathit{gn}, \mathit{xn}_i)$ limbs,
the minimum taken over the nonzero inputs, and must not overlap the
values.
:::

## Trial division and removal of factors {#trial}

::: function
mp_size_t flint_mpn_remove_2exp(mp_ptr x, mp_size_t xsize, flint_bitcnt_t * bits)

Divides $(x, \mathit{xsize})$ in place by $2^v$, where $v$ is the number
of trailing zero bits, stores $v$ in `bits` and returns the new length.
$x$ must be nonzero.
:::

::: function
mp_size_t flint_mpn_remove_power_ascending(mp_ptr x, mp_size_t xsize, mp_ptr p, mp_size_t psize, ulong * exp)

Divides $(x, \mathit{xsize})$ in place by the largest power $p^v$ of
$(p, \mathit{psize})$ that divides it, stores $v$ in `exp` and returns the
new length. $x$ must be nonzero and $p > 2$. The divisibility is tested by
ascending squares $p, p^2, p^4, \ldots$, which is efficient for large
powers but has a high overhead, so it should not be used as the first
stage of trial division.
:::

::: function
int flint_mpn_factor_trial(mp_srcptr x, mp_size_t xsize, slong start, slong stop)

Searches for a factor of $(x, \mathit{xsize})$ among the primes
`flint_primes[i]` for `start` $\le i <$ `stop`. Returns the first such $i$
for which the prime divides $x$, and 0 if there is none. Requires
`start` $\ge 1$, so that 0 is never a valid index.
:::

::: function
int flint_mpn_factor_trial_tree(slong * factors, mp_srcptr x, mp_size_t xsize, slong num_primes)

Searches for factors of $(x, \mathit{xsize})$ among roughly the first
`num_primes` primes, by GCDs with products of primes organised in a
product tree; writes the indices in `flint_primes` of the primes found to
`factors` and returns their number. Requires $0 \le$ `num_primes`
$\le 3512$. If $x$ fits in a small `fmpz`, it is factored completely
instead. The product tree is computed once and cached globally (in a
thread-safe way). This is the trial division used by `fmpz_factor`.
:::

## Multimodular reduction and Chinese remaindering {#crt}

Simultaneous reduction of a multiprecision integer modulo a fixed list of
single-limb moduli (*multi-mod*) and the inverse reconstruction
(*multi-CRT*), with precomputed data stored in a `flint_mpn_crt_t`. The
moduli must be nonzero and, for Chinese remaindering, pairwise coprime;
they need not be prime. The `fmpz_comb_t` of the `fmpz` module is a thin
wrapper around this structure, and [](mpn_extras/multi_crt)
shows how to use it.

::: summary algorithms/integers/multimodular-reduction :::

::: type
flint_mpn_crt_struct
flint_mpn_crt_t

Precomputed data for multi-mod and multi-CRT with respect to the moduli
$m_0, \ldots, m_{k-1}$ with product $P$. The following fields are
public: `num_primes`, `primes` (a copy of the moduli), `prod` and
`prod_len` (the product $P$ as a limb array), and `tmp_limbs` (the size
of the workspace required by the conversion functions).
:::

::: type
flint_mpn_crt_leaf_struct

Reduction data for one leaf modulus (a product of consecutive moduli that
fits in a limb): the `nmod_t`, $B \bmod n$, $B^2 \bmod n$, the quotient
$\lfloor B^2/n \rfloor$ and a Barrett precomputation.
:::

::: macro
FLINT_MPN_CRT_MOD
FLINT_MPN_CRT_CRT

Flags (1 and 2) selecting the operations to precompute for.
:::

::: function
void flint_mpn_crt_init(flint_mpn_crt_t C, nn_srcptr primes, slong num_primes)
void flint_mpn_crt_init2(flint_mpn_crt_t C, nn_srcptr primes, slong num_primes, int flags)

Initialise `C` for the given moduli. The version with `flags` (a bitwise
or of `FLINT_MPN_CRT_MOD` and `FLINT_MPN_CRT_CRT`) only performs the
precomputations needed for the selected operations, which saves time and
memory when only one direction is needed; `flint_mpn_crt_init` selects
both. Throws an exception if a modulus is zero, or if Chinese remaindering
is selected and the moduli are not pairwise coprime. The precomputation
costs $O(M(N) \log N)$ for $N$ total bits and uses $O(N \log N)$ memory.
:::

::: function
void flint_mpn_crt_init_tuned(flint_mpn_crt_t C, nn_srcptr primes, slong num_primes, int flags, slong crt_chunk_bits, slong mod_base_bits, slong preinv_cutoff)

`flint_mpn_crt_init2` with explicit tuning parameters: the approximate
number of bits of a CRT basecase chunk, the number of bits at which the
modular reduction switches to its basecase, and the number of limbs from
which precomputed inverses are used for division. Intended for profiling.
:::

::: function
void flint_mpn_crt_clear(flint_mpn_crt_t C)

Frees the memory used by `C`.
:::

::: function
void flint_mpn_multi_mod(nn_ptr out, nn_srcptr x, slong xn, const flint_mpn_crt_t C, nn_ptr tmp)

Writes the residues of the nonnegative $(x, \mathit{xn})$ modulo the
moduli to `out[0], ..., out[num_primes - 1]`. The input may have any
length, including zero, and need not be normalised. `tmp` must have
room for `C->tmp_limbs` limbs, or be `NULL` for internal allocation.
Requires `C` to be initialised for `FLINT_MPN_CRT_MOD`.
:::

::: function
int flint_mpn_multi_crt(nn_ptr out, nn_srcptr res, const flint_mpn_crt_t C, int sign, nn_ptr tmp)

Reconstructs the integer $x$ with the residues
`res[0], ..., res[num_primes - 1]`, each reduced modulo its modulus, and
writes it as `C->prod_len` limbs, zero padded, to `out`. If `sign` is zero,
$x$ is the unique solution with $0 \le x < P$ and the return value is 0.
Otherwise $x$ is the unique solution with $-P < 2x \le P$, its absolute
value is written, and the return value is 1 if $x$ is negative. `tmp` is
as for `flint_mpn_multi_mod`.
:::

::: function
void flint_mpn_multi_mod_vec(nn_ptr out, slong out_stride, nn_srcptr x, slong xn, slong len, const flint_mpn_crt_t C, nn_ptr tmp)
void flint_mpn_multi_crt_vec(nn_ptr out, slong out_stride, int * negative, nn_srcptr res, slong res_stride, slong len, const flint_mpn_crt_t C, int sign, nn_ptr tmp)

Vector versions. In the mod version, `x` is a packed array of `len`
nonnegative integers of `xn` limbs each, and the residue of entry $i$
modulo modulus $l$ is written to `out[l * out_stride + i]` (so
`out_stride` must be at least `len`). In the CRT version, the residue of
entry $i$ modulo modulus $l$ is read from `res[l * res_stride + i]`, entry
$i$ of the output is written to `out + i * out_stride` (with
`out_stride` $\ge$ `C->prod_len`), and if `sign` is nonzero the sign flags
are written to `negative[i]` (`negative` may be `NULL` when `sign` is
zero). These are faster than loops over the entries when the entries are
small compared to $P$, since the precomputed tables are then traversed
once per block of entries.
:::

::: function
void flint_mpn_multi_mod_once(nn_ptr out, nn_srcptr x, slong xn, nn_srcptr primes, slong num_primes)
int flint_mpn_multi_crt_once(nn_ptr out, slong * outn, nn_ptr prod, nn_srcptr res, nn_srcptr primes, slong num_primes, int sign)

One-shot versions of `flint_mpn_multi_mod` and `flint_mpn_multi_crt`
that need no precomputed structure. In the CRT version, `out` must have
room for `num_primes` limbs; the length of the product $P$ is written to
`outn`, and if `prod` is not `NULL`, $P$ itself is written there (with
room for `num_primes` limbs). The return value is as for
`flint_mpn_multi_crt`.

::: details
The mod version creates a mod-only structure internally. The CRT version
instead traverses the subproduct tree depth first, computing subproducts
and cofactors on the fly and freeing them as soon as possible, so that
its memory use is a small multiple of the size of $P$ rather than
$O(N \log N)$, which matters for very large reconstructions. It is also
somewhat faster than initialising a structure for a single use. Large instances are
parallelised over the nodes of the tree when several threads are
available.
:::
:::

::: function
ulong flint_mpn_crt_mod_leaf(nn_srcptr a, slong an, const flint_mpn_crt_t C, slong j)

Returns the nonnegative $(a, \mathit{an})$ reduced modulo the $j$-th leaf
modulus of `C` (the product of one or more consecutive batched moduli;
leaf $j$ is modulus $j$ when all moduli exceed 32 bits). Requires
`an` $\le$ `C->mod_pow_limbs`. This exposes the basecase of
`flint_mpn_multi_mod`; with a single modulus it is a faster replacement
for `mpn_mod_1` when many reductions with the same modulus are needed.
:::

## Random numbers {#random}

::: function
void flint_mpn_urandomb(mp_ptr rp, flint_rand_t state, flint_bitcnt_t n)

Sets $(\mathit{rp}, \lceil n / \bits \rceil)$ to a uniformly random
integer in $[0, 2^n)$.
:::

::: function
void flint_mpn_urandomm(mp_ptr rp, flint_rand_t state, mp_srcptr xp, mp_size_t xn)

Sets $(\mathit{rp}, \mathit{xn})$ to a uniformly random integer in
$[0, x)$. Requires $x_{xn-1} \ne 0$. All $\mathit{xn}$ limbs are written
even if the value is shorter.
:::

::: function
void flint_mpn_rrandom(mp_ptr rp, flint_rand_t state, mp_size_t n)
void flint_mpn_rrandomb(mp_ptr rp, flint_rand_t state, flint_bitcnt_t n)

Set $(\mathit{rp}, n)$, respectively $(\mathit{rp}, \lceil n/\bits \rceil)$
with an $n$-bit value, to a random integer with long runs of zeros and
ones in its binary representation. Such numbers trigger corner cases
(carry propagation, near-overflow) more often than uniform ones and are
used in tests.
:::

## Internal functions {#internal}

### Division algorithms

::: function
void _flint_mpn_tdiv_qr(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
void _flint_mpn_tdiv_qr_newton(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
void _flint_mpn_tdiv_qr_unbalanced(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
void _flint_mpn_tdiv_qr_preinv(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, mp_srcptr binv, mp_size_t binvn)
void _flint_mpn_tdiv_qr_preinvn(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

The dispatcher behind `flint_mpn_tdiv_qr` and the individual algorithms
for large operands, with its interface; in all of them $r$ may be `NULL`.
The Newton version views $a$ and $b$ as fixed-point numbers and computes
$m + 2$ fraction limbs of $a/b$, where $m = \mathit{an} - \mathit{bn} + 1$;
since the error is below $4 B^{-2}$ at the integer scale, the integer
part is certified whenever the first fraction limb lies in $[2, B-2]$,
and the remainder then follows from a low product; otherwise the
candidate quotient is corrected by $O(1)$ steps using a full product. The
unbalanced version (requires $\mathit{an} > 2\mathit{bn}$ and
$\mathit{bn} \ge 3$) divides blockwise, sharing one approximate inverse
of $b$ between the $2\mathit{bn} \times \mathit{bn}$ blocks. The preinv
version requires $(\mathit{binv}, \mathit{binvn} + 2)$ to be the output of
`_mp_real_inv_newton(binv, b, bn, binvn)` with $\mathit{binvn} \ge m + 2$
and $\mathit{an} \ge m + 2$. The preinvn version uses
`flint_mpn_preinvn` and `flint_mpn_divrem_preinvn`.
:::

::: function
void _flint_mpn_tdiv_qr_small(mp_ptr q, mp_ptr r, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

Euclidean division for $1 \le \mathit{bn} \le$ `FLINT_MPN_DIV_SMALL_BN`
(7), with the conventions of `flint_mpn_tdiv_qr` ($r$ may be `NULL`; no
normalisation of $b$ is needed). There is one code path per divisor
length, in which the normalised divisor and the partial remainder are
kept in registers and the dividend is normalised on the fly; each
quotient limb takes a 3/2 division (2/1 for $\mathit{bn} = 1$) and, for
$\mathit{bn} \ge 3$, the subtraction of a $1 \times (\mathit{bn} - 2)$
product, as in GMP's `mpn_sbpi1_div_qr` but without function calls or
memory traffic. The cost is linear in $\mathit{an}$; beyond about
$3\mathit{bn}$ limbs the schoolbook division is faster for
$\mathit{bn} \ge 4$.
:::

::: macro
FLINT_MPN_DIV_SMALL_BN
FLINT_MPN_DIV_SMALL_SHAPE(an, bn, want_r)
FLINT_MPN_DIVREM_1_USE_HW(an, d)
FLINT_MPN_DIVEXACT_SMALL_BN

The largest divisor length (7) of `_flint_mpn_tdiv_qr_small`; whether it
is used for a given shape (with or without remainder); whether a one-limb
division uses a chain of hardware divisions (below
`FLINT_MPN_DIVREM_1_HW_CUTOFF` or `FLINT_MPN_DIVREM_1_NORM_HW_CUTOFF`
dividend limbs, only with `FLINT_PREINVERT_LIMB_USE_NATIVE`); and the
largest divisor length (7) of the register-based exact division.
:::

::: function
mp_limb_t _flint_mpn_divrem_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
mp_limb_t _flint_mpn_div_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
mp_limb_t _flint_mpn_divapprox_basecase_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv)
mp_limb_t _flint_mpn_divrem_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv, mp_ptr tp)
mp_limb_t _flint_mpn_divapprox_preinv1(mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn, mp_limb_t dinv, mp_ptr tp)
mp_limb_t _flint_mpn_divrem_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n, mp_limb_t dinv, mp_ptr tp)
mp_limb_t _flint_mpn_divapprox_n_divconquer_preinv1(mp_ptr qp, mp_ptr np, mp_srcptr dp, mp_size_t n, mp_limb_t dinv, mp_ptr tp)
void _flint_mpn_tdiv_qr_divconquer(mp_ptr Q, mp_ptr R, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn)
void _flint_mpn_tdiv_q_divconquer(mp_ptr Q, mp_srcptr A, mp_size_t An, mp_srcptr B, mp_size_t Bn)

Division with a normalised divisor and the 3/2 inverse `dinv` =
`flint_mpn_preinv1` of its top two limbs, ported from GMP's internal
`mpn_sbpi1_*`, `mpn_dcpi1_*`, `mpn_tdiv_qr` and `mpn_div_q` functions,
with FLINT's multiplication in the divide and conquer steps
[@BurZie1998].

::: details
The schoolbook functions require $\mathit{dn} \ge 3$ and
$\mathit{nn} \ge \mathit{dn}$. `_flint_mpn_divrem_basecase_preinv1`
divides $(\mathit{np}, \mathit{nn})$ in place by $(\mathit{dp},
\mathit{dn})$: the remainder is left in $(\mathit{np}, \mathit{dn})$, the
low $\mathit{nn} - \mathit{dn}$ limbs of the quotient are written to
`qp` and its top limb is returned. `_flint_mpn_div_basecase_preinv1`
computes the same quotient but develops the last $\mathit{dn} - 2$
quotient limbs with a truncated divisor (correcting the result when
needed), and `_flint_mpn_divapprox_basecase_preinv1` computes an
approximate quotient (correct or one too large) the same way; both
destroy $(\mathit{np}, \mathit{nn})$.

`_flint_mpn_divrem_preinv1` and `_flint_mpn_divapprox_preinv1` are the
analogues for $\mathit{dn} \ge 2$ with divide and conquer steps above
`FLINT_MPN_DIV_DC_CUTOFF` and `FLINT_MPN_DIVAPPR_DC_CUTOFF` limbs (ports
of `mpn_dcpi1_div_qr` and `mpn_dcpi1_divappr_q`); `tp` must have room for
$\mathit{dn}$ limbs. The balanced $2n/n$ steps ($n \ge 6$, `tp` of $n$
limbs) are `_flint_mpn_divrem_n_divconquer_preinv1` and
`_flint_mpn_divapprox_n_divconquer_preinv1`. The latter is not correct up
to one unit: each level of its recursion can make the quotient a few
units too large (never too small; at most 8 units in tests with
$n \le 200$), so it is used with a guard limb.

`_flint_mpn_tdiv_qr_divconquer` and `_flint_mpn_tdiv_q_divconquer` are
Euclidean divisions for any divisor with $\mathit{bn} \ge 2$ ($r$ may be
`NULL`). When the quotient is about as long as the divisor or longer, a
normalised copy is divided as above. For short quotients ($m$ limbs),
the top $2m$ limbs of $a$ are divided by the top $m$ limbs of $b$ and
the quotient (at most two too large) is corrected using its product with
the ignored divisor limbs. Quotients of $m \le \mathit{bn} + 2$ limbs
from `FLINT_MPN_DIVAPPROX_SHORT_CUTOFF` limbs use a short division in the
style of Mulders [@Mul2000; @HZ2011]: the quotient with one guard limb
is computed from the top $m + 1$ limbs of $a$ and $b$ only, recursively,
the correction for the rest of the divisor being needed only above the
current precision, i.e. as a high product.

The cutoffs are machine dependent and measured by
`src/mpn_extras/tune/tune-div.c` for `flint-mparam.h`; the divide and
conquer cutoffs must be at least 6.
:::
:::

::: function
void _flint_mpn_inv_basecase(mp_ptr q, mp_srcptr x, mp_size_t xn, mp_size_t n)

`flint_mpn_inv` without Newton inversion. It is used by
`_mp_real_inv_newton_basecase`, so that the mutual calls between the
Newton inversion and the division only ever go to lower precisions,
whatever the cutoffs.
:::

::: function
void _flint_mpn_divexact(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
void _flint_mpn_divexact_hensel(mp_ptr q, mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)

The part of `flint_mpn_divexact` for divisors of at least two limbs, and
the version that always uses the Newton Hensel division
`flint_mpn_bdiv_q`.
:::

::: function
int _flint_mpn_divisible(mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn)
int _flint_mpn_divisible_bdiv(mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, unsigned int v)
int _flint_mpn_divisible_small(mp_srcptr a, mp_size_t an, mp_srcptr b, mp_size_t bn, unsigned int v)

The parts of `flint_mpn_divisible`, for $b = 2^v b'$ with $b'$ odd, given
that $2^v$ divides $a$. The small version ($b'$ of at most
`FLINT_MPN_DIVEXACT_SMALL_BN` limbs, $\mathit{bn} \ge 2$, $v < 64$,
$a_{an-1} \ne 0$) streams the shifted $a$ through a Hensel remainder
computation with $b'$ and the running remainder held in registers; $b'$
divides $a$ if and only if the final remainder is zero.
:::

### Half-GCD

::: type
flint_mpn_hgcd_matrix1_struct
flint_mpn_hgcd_matrix_struct

A $2 \times 2$ matrix of words `u[i][j]`, and a $2 \times 2$ matrix of
limb arrays `p[i][j]` with entries of `n` limbs (normalised in the sense
that one of the top limbs is nonzero) and room for `alloc` limbs, as
GMP's `struct hgcd_matrix1` and `struct hgcd_matrix`.
:::

::: function
int _flint_mpn_hgcd2(mp_limb_t ah, mp_limb_t al, mp_limb_t bh, mp_limb_t bl, flint_mpn_hgcd_matrix1_struct * M)
mp_size_t _flint_mpn_hgcd_mul_matrix1_vector(const flint_mpn_hgcd_matrix1_struct * M, mp_ptr rp, mp_srcptr ap, mp_ptr bp, mp_size_t n)
mp_size_t _flint_mpn_hgcd_mul_matrix1_inverse_vector(const flint_mpn_hgcd_matrix1_struct * M, mp_ptr rp, mp_srcptr ap, mp_ptr bp, mp_size_t n)
void _flint_mpn_gcd_22(mp_ptr g, mp_limb_t u1, mp_limb_t u0, mp_limb_t v1, mp_limb_t v0)

The word-level half-GCD of the two-limb numbers $(a_h, a_l)$ and $(b_h,
b_l)$, returning 0 if no reduction is possible; the products $(r, b)
\gets (a, b) M$ and $(r, b) \gets M^{-1} (a, b)$ by such a matrix; and
the GCD of two odd two-limb numbers (GMP's `mpn_gcd_22` if available,
otherwise a port of GMP's C version). These are
ports of the GMP functions of the same names, except that the two
products form both outputs in a single pass, and allow `rp` to be equal
to `ap`.
:::

::: function
void _flint_mpn_hgcd_matrix_init(flint_mpn_hgcd_matrix_struct * M, mp_size_t n, mp_ptr p)
void _flint_mpn_hgcd_matrix_update_q(flint_mpn_hgcd_matrix_struct * M, mp_srcptr qp, mp_size_t qn, unsigned int col, mp_ptr tp)
void _flint_mpn_hgcd_matrix_mul_1(flint_mpn_hgcd_matrix_struct * M, const flint_mpn_hgcd_matrix1_struct * M1, mp_ptr tp)
void _flint_mpn_hgcd_matrix_mul(flint_mpn_hgcd_matrix_struct * M, const flint_mpn_hgcd_matrix_struct * M1, mp_ptr tp)
mp_size_t _flint_mpn_hgcd_matrix_adjust(const flint_mpn_hgcd_matrix_struct * M, mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t p, mp_ptr tp)
mp_size_t _flint_mpn_hgcd_matrix22_mul_strassen_itch(mp_size_t rn, mp_size_t mn)
void _flint_mpn_hgcd_matrix22_mul_strassen(mp_ptr r0, mp_ptr r1, mp_ptr r2, mp_ptr r3, mp_size_t rn, mp_srcptr m0, mp_srcptr m1, mp_srcptr m2, mp_srcptr m3, mp_size_t mn, mp_ptr tp)

The matrix operations of the half-GCD, as GMP's `mpn_hgcd_matrix_*`.
`_flint_mpn_hgcd_matrix_init` initialises $M$ to the identity for inputs
of $n$ limbs, using `FLINT_MPN_HGCD_MATRIX_INIT_ITCH(n)` limbs at `p`. The
product and the adjustment ($(a, b) \gets M^{-1} (a, b)$ for the low $p$
limbs of $a$ and $b$) are computed in transform space for large operands
([](algorithms/integers/gcd#use)). `_flint_mpn_hgcd_matrix22_mul_strassen` is GMP's
`mpn_matrix22_mul`, the classical or a Strassen-like product of $2 \times
2$ matrices, which the half-GCD uses below the transform-space range; it
is not meant as a general $2 \times 2$ matrix product.
:::

::: function
mp_size_t _flint_mpn_gcd_subdiv_step(mp_ptr ap, mp_ptr bp, mp_size_t n, mp_size_t s, _flint_mpn_gcd_subdiv_step_hook * hook, void * ctx, mp_ptr tp)
mp_size_t _flint_mpn_hgcd_step(mp_size_t n, mp_ptr ap, mp_ptr bp, mp_size_t s, flint_mpn_hgcd_matrix_struct * M, mp_ptr tp)
mp_size_t _flint_mpn_hgcd_itch(mp_size_t n)
mp_size_t _flint_mpn_hgcd(mp_ptr ap, mp_ptr bp, mp_size_t n, flint_mpn_hgcd_matrix_struct * M, mp_ptr tp)
mp_size_t _flint_mpn_hgcd_appr_itch(mp_size_t n)
int _flint_mpn_hgcd_appr(mp_ptr ap, mp_ptr bp, mp_size_t n, flint_mpn_hgcd_matrix_struct * M, mp_ptr tp)
mp_size_t _flint_mpn_hgcd_reduce_itch(mp_size_t n, mp_size_t p)
mp_size_t _flint_mpn_hgcd_reduce(flint_mpn_hgcd_matrix_struct * M, mp_ptr ap, mp_ptr bp, mp_size_t n, mp_size_t p, mp_ptr tp)
void _flint_mpn_gcdext_hook(void * p, mp_srcptr gp, mp_size_t gn, mp_srcptr qp, mp_size_t qn, int d)
mp_size_t _flint_mpn_gcdext_lehmer_n(mp_ptr gp, mp_ptr up, mp_size_t * usize, mp_ptr ap, mp_ptr bp, mp_size_t n, mp_ptr tp)

Ports of the GMP functions of the same names (without the leading
`_flint`). `_flint_mpn_hgcd` reduces $a, b$ of $n$ limbs (one of the top
limbs nonzero) until $|a - b|$ fits in $\lfloor n/2 \rfloor + 1$ limbs,
accumulating the reduction in $M$ (initialised for $n$ limbs) so that
$(a, b)^T$ is $M$ times the reduced pair, and returns the new size, or 0
if no reduction is possible; `tp` must have room for
`_flint_mpn_hgcd_itch(n)` limbs. `_flint_mpn_hgcd_appr` computes an
approximation of the same matrix from destroyed inputs.
:::

::: macro
FLINT_MPN_HGCD_CUTOFF
FLINT_MPN_HGCD_APPR_CUTOFF
FLINT_MPN_HGCD_REDUCE_CUTOFF
FLINT_MPN_GCD_DC_CUTOFF
FLINT_MPN_GCDEXT_DC_CUTOFF
FLINT_MPN_HGCD_MATRIX22_STRASSEN_CUTOFF

Cutoffs in limbs, with the meanings of GMP's `HGCD_THRESHOLD`,
`HGCD_APPR_THRESHOLD`, `HGCD_REDUCE_THRESHOLD`, `GCD_DC_THRESHOLD`,
`GCDEXT_DC_THRESHOLD` and `MATRIX22_STRASSEN_THRESHOLD`.
:::

### Fixed-size arithmetic macros

::: macro
NN_ADD_2(r, u, v)
NN_ADD_3(r, u, v)
NN_ADD_4(r, u, v)
NN_ADD_5(r, u, v)
NN_ADD_6(r, u, v)
NN_ADD_7(r, u, v)
NN_ADD_8(r, u, v)
NN_SUB_2(r, u, v)
NN_SUB_3(r, u, v)
NN_SUB_4(r, u, v)
NN_SUB_5(r, u, v)
NN_SUB_6(r, u, v)
NN_SUB_7(r, u, v)
NN_SUB_8(r, u, v)

$r \gets u + v$ and $r \gets u - v$ modulo $B^k$ for the $k$-limb arrays
$r$, $u$, $v$ ($2 \le k \le 8$), as straight-line carry chains (the
`add_sssaaaaaa` family of [`longlong.h`](longlong)), without carry in or
out.
:::

::: macro
NN_DOTREV_S3_1X1(s2, s1, s0, u, v, n)
NN_DOTREV_S3_1X1_HIGH(s2, s1, u, v, n)
NN_DOTREV_S3_A3_1X1(r2, r1, r0, s2, s1, s0, u, v, n)
NN_ADDMUL_S2_A2_1X1(r1, r0, s1, s0, x, y)

Reversed dot products of limbs into a three-limb accumulator:
$(s_2, s_1, s_0) = \sum_{i<n} u_i v_{n-1-i}$ (for $n \ge 2$); the same
summing only the high limbs of the products; the sum added to a given
three-limb value (for $n \ge 1$, the result fitting in three limbs); and
a two-limb $(r_1, r_0) = (s_1, s_0) + xy$ without carry out. They are the
building blocks of the generic basecases.
:::

::: macro
MPN_LEFT_SHIFT_HI(hi, lo, n)
MPN_RIGHT_SHIFT_LOW(hi, lo, n)

The high limb of the two-limb $(\mathit{hi}, \mathit{lo})$ shifted left by
$0 \le n < \bits$ bits, and the low limb of it shifted right.
:::

::: macro
eudiv_qrnd(q, r, n, d, t)
eudiv_qqrnnd(q1, q0, r0, n1, n0, d0, t)
eudiv_qrrnndd(q0, r1, r0, n1, n0, d1, d0, t)

Euclidean divisions of one or two limbs by one or two limbs where the
quotient is expected to be small (as in a continued fraction expansion):
the quotient is found by repeated subtraction for the first few values
before falling back to a division. They require $n > d > 0$; the
two-by-two version requires $d_1 \ne 0$. `t` is a prefix for the local
variables.
:::

### Assembly kernels and tables

::: type
flint_mpn_mul_func_t
flint_mpn_mul_func_n_t
flint_mpn_sqr_func_t
flint_mpn_mulhigh_normalised_func_t
flint_mpn_sqrhigh_normalised_func_t

Function pointer types of the hard-coded multiplication routines:
`mp_limb_t (*)(mp_ptr, mp_srcptr, mp_srcptr)` (fixed sizes),
`mp_limb_t (*)(mp_ptr, mp_srcptr, mp_srcptr, mp_size_t)` (one fixed size,
ARM64), `mp_limb_t (*)(mp_ptr, mp_srcptr)` (squaring), and the
normalised high products returning an `mp_limb_pair_t`.
:::

::: macro
FLINT_HAVE_MUL_FUNC(n, m)
FLINT_HAVE_MUL_N_FUNC(n)
FLINT_HAVE_SQR_FUNC(n)
FLINT_HAVE_MULLOW_FUNC(n)
FLINT_HAVE_MULHIGH_FUNC(n)
FLINT_HAVE_SQRHIGH_FUNC(n)
FLINT_HAVE_MULHIGH_NORMALISED_FUNC(n)
FLINT_HAVE_SQRHIGH_NORMALISED_FUNC(n)

Whether a hard-coded routine exists for the given sizes on the current
architecture. The tables of routines are `flint_mpn_mul_func_tab`
(indexed `[n][m]`), `flint_mpn_mul_n_func_tab`,
`flint_mpn_mul_func_n_tab` (ARM64), `flint_mpn_sqr_func_tab`,
`flint_mpn_mullow_func_tab`, `flint_mpn_mulhigh_func_tab`,
`flint_mpn_sqrhigh_func_tab`, `flint_mpn_mulhigh_normalised_func_tab` and
`flint_mpn_sqrhigh_normalised_func_tab`. See
[](algorithms/integers/basecase-multiplication) for the sizes covered.
:::

::: macro
FLINT_MPN_MUL_HARD(rp, xp, xn, yp, yn)
FLINT_MPN_MUL_N_HARD(rp, xp, yp, n)
FLINT_MPN_SQR_HARD(rp, xp, n)

Call the hard-coded routine for the given sizes; only valid when the
corresponding `FLINT_HAVE_*_FUNC` macro is true.
:::

::: function
mp_limb_t _flint_mpn_mulhigh_n_basecase2(mp_ptr rp, mp_srcptr xp, mp_srcptr yp, mp_size_t n)
mp_limb_t _flint_mpn_sqrhigh_basecase_even(mp_ptr rp, mp_srcptr xp, mp_size_t n)
mp_limb_t _flint_mpn_sqrhigh_basecase_odd(mp_ptr rp, mp_srcptr xp, mp_size_t n)

The basecase high product used by `_flint_mpn_mulhigh_n`
(`_flint_mpn_mulhigh_n_recursive` up to 22 limbs and
`_flint_mpn_mulhigh_basecase` above on x86-64 with ADX), and the two
assembly halves of `_flint_mpn_sqrhigh_basecase` on x86-64 with ADX, for
inputs of $2n$ and $2n + 1$ limbs.
:::
