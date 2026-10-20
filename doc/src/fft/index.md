---
title: fft.h
header: fft.h
description: Schönhage–Strassen FFTs over $\Z/(2^N + 1)\Z$ for multiplying large integers and polynomials with large integer coefficients.
category: Fast Fourier transforms
order: 1
---

This module implements the Schönhage–Strassen method [@SS1971]: fast
Fourier transforms over the rings $\Z/(2^N + 1)\Z$, in which $2$ is a
root of unity of order $2N$, so that every multiplication by a root of
unity is a shift. On top of the transforms it provides multiplication of
large integers (`flint_mpn_mul_fft_main`), convolutions of vectors of
large integers (`fft_convolution`, with a variant for a fixed operand)
and multiplication modulo $2^N + 1$ (`fft_mulmod_2expp1`). The
implementation, by William Hart, uses the $\sqrt 2$ trick to double the
transform length, van der Hoeven's truncated FFT [@vdH2004] to avoid
jumps in running time at powers of two, Bailey's matrix Fourier algorithm
[@Bai1990] for cache efficiency, and threads for large transforms.

Most users never call this module directly: `flint_mpn_mul`, `fmpz_mul`
and `fmpz_poly_mul` choose an algorithm themselves. FLINT uses the module
as follows.

- `flint_mpn_mul` (and so `fmpz_mul` and everything built on it) calls
  `flint_mpn_mul_fft_main` for operands of at least
  `FLINT_FFT_MUL_THRESHOLD` limbs (32000 on all tuned targets) only in
  builds *without* the `fft_small` module, that is, with 32-bit limbs
  or on processors without AVX2 or NEON. Where `fft_small` is
  available, which includes all common 64-bit desktop and server
  processors, it has replaced this module for integer multiplication:
  `mpn_mul_default_mpn_ctx` is used from a few hundred to 1000 limbs on,
  depending on the target.
- `fmpz_poly` multiplication (`_fmpz_poly_mulmid` and everything that
  dispatches through it: products, low, high and middle products,
  squaring) uses Schönhage–Strassen convolutions (`fft_convolution`) for
  polynomials with large coefficients. Without `fft_small` they are the
  main algorithm for large inputs. With `fft_small`, Kronecker
  substitution on top of `fft_small` wins almost everywhere, and the
  convolutions are used only for medium-size squarings, for huge
  products, with four or more threads, and for middle products in two
  bands where a cyclic transform saves a quarter of the work
  (`fft_convolution_choose_mid`). Their pointwise products are then
  computed with the negacyclic small-prime engine of `fft_small`
  (`fft_convolution_sd_fft`).
- `fmpz_poly_mul_SS_precache_init` and `fmpz_poly_mullow_SS_precache`
  use `fft_precache` and `fft_convolution_precache_sd_fft`, and the
  Gaussian integer polynomials of `fmpzi` use the same convolutions for
  huge products.
- `fmpz_mat_mul_fft` uses the transforms directly; `fmpz_mat_mul` itself
  calls `fmpz_mat_mul_fft_small` instead.

The guide [](fft/using) shows how to multiply integers and polynomials
with the module, including the memory layout that the transforms expect.

::: summary algorithms/fft/schonhage-strassen :::

## Conventions {#conventions}

Throughout, $B = 2^{\bits}$.

**Residues.** An element of $\Z/p\Z$, $p = 2^N + 1$ with
$N = \mathit{limbs} \cdot \bits$, is stored in $\mathit{limbs} + 1$
limbs as a signed integer in two's complement. Additions and
subtractions are carried out without reduction: the top limb absorbs the
carries and borrows, and only the value modulo $p$ matters. A residue is
*reduced* if it lies in $[0, p)$, that is, its top limb is $0$, or its
top limb is $1$ and all other limbs are $0$ (the value $2^N \equiv -1$).
`mpn_normmod_2expp1` reduces a residue.

**Transform parameters.** A transform of length $2n$,
$n = 2^{\mathit{depth}}$, works modulo $2^N + 1$ with $N = nw$ for a
parameter $w \ge 1$, and $N$ must be a multiple of $\bits$
($\mathit{limbs} = nw/\bits$). Then $\theta = 2^w$ satisfies
$\theta^n = 2^N \equiv -1$, so it is a root of unity of order exactly
$2n$, and $2$ is invertible. The *sqrt2* transforms have length $4n$ over
the same ring and use the root $\sqrt 2^{\,w}$, where
$\sqrt 2 = 2^{3N/4} - 2^{N/4}$; see
[The $\sqrt 2$ trick](algorithms/fft/schonhage-strassen#sqrt2).

**Coefficient arrays.** A vector of residues is an array `ii` of
pointers, each to a block of $\mathit{limbs} + 1$ limbs. The transforms
reorder data by swapping pointers: they exchange entries of `ii` with
each other and with the scratch pointers `*t1`, `*t2` and `*temp`, so
after a call `ii[i]` may point to a block that was passed as scratch.
The blocks should therefore be carved out of one allocation that is freed
as a whole. Distinct pointers must point to non-overlapping blocks.

**Order.** The forward transforms produce their output in bit-reversed
order, and the inverse transforms take their input in that order. Inverse
transforms do not divide by the transform length unless stated.

**Scratch space and threads.** The functions that may run threads
(`fft_convolution`, `fft_convolution_sd_fft` and the `_outer` and
`_inner` matrix Fourier functions) take arrays `t1`, `t2`, `s1` (or
`temp`) and `tt` with one entry per thread, that is, at least
`flint_get_num_threads()` entries at the time of the call. Each entry of
`t1`, `t2` and `s1` points to $\mathit{limbs} + 1$ limbs, each entry of
`tt` to $2(\mathit{limbs} + 1)$ limbs. The other functions use only the
first entry.

None of the preconditions below is checked, except by `FLINT_ASSERT` in
builds configured with `--enable-assert`.

## Integer multiplication {#mul}

::: function
void flint_mpn_mul_fft_main(mp_ptr r1, mp_srcptr i1, mp_size_t n1, mp_srcptr i2, mp_size_t n2)

Sets `(r1, n1 + n2)` to the product of `(i1, n1)` and `(i2, n2)`, choosing
the transform parameters from tuning tables and calling
`mul_truncate_sqrt2` or, for large products, `mul_mfa_truncate_sqrt2`.
Requires $n_1, n_2 \ge 1$, in either order. The output must not overlap
the inputs; `i1 == i2` gives a squaring. Intended for operands of
thousands of limbs and more; assert builds reject products shorter than
about 57 limbs. Uses threads for large products.

::: details
The parameter choice, the tuning tables and the sizes at which each
variant is used are given in
[Choice of parameters](algorithms/fft/schonhage-strassen#parameters).
:::
:::

::: function
void mul_truncate_sqrt2(mp_ptr r1, mp_srcptr i1, mp_size_t n1, mp_srcptr i2, mp_size_t n2, flint_bitcnt_t depth, flint_bitcnt_t w)

Sets `(r1, n1 + n2)` to the product of `(i1, n1)` and `(i2, n2)` with a
truncated sqrt2 convolution of length at most $4n$,
$n = 2^{\mathit{depth}}$, modulo $2^{nw} + 1$. The operands are split
into chunks of $b = \lfloor (nw - \mathit{depth} - 1)/2 \rfloor$ bits; if
this gives $j_1$ and $j_2$ chunks, requires $j_1 + j_2 - 1 \le 4n$.
Requires $nw$ to be a multiple of $\bits$. The pointwise products use
`flint_mpn_mulmod_2expp1_basecase`. Single-threaded except for the
splitting.
:::

::: function
void mul_mfa_truncate_sqrt2(mp_ptr r1, mp_srcptr i1, mp_size_t n1, mp_srcptr i2, mp_size_t n2, flint_bitcnt_t depth, flint_bitcnt_t w)

As `mul_truncate_sqrt2`, with the same requirements, but using the
cache-friendly matrix Fourier algorithm, threads, and
`fft_mulmod_2expp1` for the pointwise products. Any $w$ is allowed, but
if the limb count $nw/\bits$ exceeds `FFT_MULMOD_2EXPP1_CUTOFF` it must
be a value returned by `fft_adjust_limbs`. Requires
$\mathit{depth} \ge 2$.
:::

## Convolutions {#convolution}

::: function
void fft_convolution(mp_limb_t ** ii, mp_limb_t ** jj, slong depth, slong limbs, slong trunc, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)

Multiplies the polynomials with coefficients `ii` and `jj` modulo
$p = 2^N + 1$, $N = \mathit{limbs} \cdot \bits$, with a truncated sqrt2
transform of length $4n$, $n = 2^{\mathit{depth}}$, and writes the first
`trunc` coefficients of the product, reduced, to `ii`. Requires
$2n < \mathit{trunc} \le 4n$ (for $\mathit{depth} \le 6$ and an even
$w = N/n$, $0 < \mathit{trunc} \le 4n$ suffices; other values can crash
or hang) and that the product has at most `trunc`
coefficients; for $\mathit{trunc} = 4n$ the product is otherwise
computed modulo $x^{4n} - 1$. Requires $n$ to divide $N$ and, if `limbs`
exceeds `FFT_MULMOD_2EXPP1_CUTOFF`, `limbs` to be a value returned by
`fft_adjust_limbs`. Both arrays have $4n$ entries, and the entries
beyond each input must be zero. `jj` is overwritten unless `jj == ii`
(squaring). The scratch arrays need one entry per thread.
:::

::: function
void fft_convolution_sd_fft(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx, mp_limb_t ** ii, mp_limb_t ** jj, slong depth, slong limbs, slong trunc, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)

As `fft_convolution`, but if `sdctx` is not `NULL`, the pointwise
products modulo $2^N + 1$ are computed by `sd_fft_mpn_mulmod_2expp1` of
`fft_small` with this context, which must have been initialised for
$N = \mathit{limbs} \cdot \bits$ and the digit count returned by
`fft_adjust_limbs_sd_fft`. With `NULL` it is `fft_convolution`.
:::

::: function
slong fft_convolution_choose_mid(slong * trunc, slong len_out, slong nlo, slong nhi)

Plans a convolution of which only the coefficients with indices in
$[\mathit{nlo}, \mathit{nhi})$ of a product of length `len_out` are
needed, where $0 \le \mathit{nlo} < \mathit{nhi} \le \mathit{len\_out}$.
Returns the `depth` and sets `trunc` for `fft_convolution`: either the
truncated transform with $\mathit{trunc} = \mathit{len\_out}$, or a
cyclic transform of length $\mathit{trunc} = 4 \cdot 2^{\mathit{depth}}$
whose wrapped terms all fall below `nlo`, whichever processes fewer
points. The truncated choice violates $\mathit{trunc} > 2n$ when
$\mathit{len\_out} \le 2$, so callers handle that case separately.
:::

## Precomputed transforms {#precache}

::: function
void fft_precache(mp_limb_t ** jj, slong depth, slong limbs, slong trunc, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1)

Replaces the coefficients `jj` by their reduced truncated sqrt2
transform, for use as the fixed operand of `fft_convolution_precache`
with the same `depth`, `limbs` and `trunc`. The requirements are those of
`fft_convolution`. Single-threaded.
:::

::: function
void fft_convolution_precache(mp_limb_t ** ii, mp_limb_t ** jj, slong depth, slong limbs, slong trunc, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)

As `fft_convolution`, where `jj` has been transformed by `fft_precache`
with the same parameters. `jj` is not modified, so it can serve for any
number of products. Single-threaded.
:::

::: function
void fft_convolution_precache_sd_fft(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx, mp_limb_t ** ii, mp_limb_t ** jj, slong depth, slong limbs, slong trunc, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)

As `fft_convolution_precache`, with the optional pointwise engine of
`fft_convolution_sd_fft`, whose pointwise products use threads.
:::

## Multiplication modulo $2^N + 1$ {#mulmod}

::: function
void fft_mulmod_2expp1(mp_limb_t * r, mp_limb_t * i1, mp_limb_t * i2, mp_size_t n, mp_size_t w, mp_limb_t * tt)

Sets `r` to the product of the reduced residues `i1` and `i2` modulo
$2^{nw} + 1$, reduced, where $nw$ is a multiple of $\bits$. Up to
`FFT_MULMOD_2EXPP1_CUTOFF` limbs (128 on all targets) it uses
`flint_mpn_mulmod_2expp1_basecase`. Above, it uses a negacyclic FFT
convolution, which requires the limb count $nw/\bits$ to be a value
returned by `fft_adjust_limbs`; for other limb counts the result is
wrong. `r` may alias the inputs. `tt` must have room for $2nw/\bits$
limbs.
:::

::: function
slong fft_adjust_limbs(mp_size_t limbs)

Returns the least limb count $\ge$ `limbs` that `fft_mulmod_2expp1` can
handle: `limbs` itself if it is at most `FFT_MULMOD_2EXPP1_CUTOFF`, and
otherwise a value rounded up to a suitable multiple, at most the next
power of two.
:::

::: function
slong fft_adjust_limbs_sd_fft(slong tight, slong n_outer, slong * m_out)

Chooses the coefficient ring of a Schönhage–Strassen convolution with
transform parameter $n$ = `n_outer` and coefficients of at least `tight`
limbs, taking the negacyclic engine of `fft_small` into account. If a
ring $N' = mb$ suitable for that engine is found, sets `*m_out` to its
digit count $m$ and returns $N'/\bits$; otherwise sets `*m_out` to zero
and returns `fft_adjust_limbs(tight)`. Without `fft_small` it always does
the latter.
:::

::: macro
FLINT_FFT_SD_FFT_MIN_RING

The smallest ring size $N$ in bits (16384) for which
`fft_adjust_limbs_sd_fft` selects the negacyclic engine.
:::

::: function
void _fft_mulmod_2expp1(mp_limb_t * r1, mp_limb_t * i1, mp_limb_t * i2, mp_size_t r_limbs, flint_bitcnt_t depth, flint_bitcnt_t w)

Sets `r1` to the reduced product of `i1` and `i2` modulo $2^N + 1$,
$N = \mathit{r\_limbs} \cdot \bits$, with a negacyclic transform of
length $2n$, $n = 2^{\mathit{depth}}$, modulo $2^{nw} + 1$ combined with
a one-limb convolution. The inputs must be less than $2^N$ (top limbs
zero), and `depth` and `w` must be as chosen by `fft_mulmod_2expp1`:
$nw = N/n$, with $N/(2n)$ a multiple of $\bits$.
:::

::: function
void fft_naive_convolution_1(mp_limb_t * r, mp_limb_t * ii, mp_limb_t * jj, mp_size_t m)

Sets $r_k = \sum_{i + j = k} \mathit{ii}_i \mathit{jj}_j - \sum_{i + j = k + m} \mathit{ii}_i \mathit{jj}_j \bmod B$
for $0 \le k < m$: the negacyclic convolution of two vectors of $m$
limbs, with wrapping word arithmetic. `r` must not alias the inputs.
:::

## Splitting and combining {#split}

::: function
mp_size_t fft_split_limbs(mp_limb_t ** poly, mp_srcptr limbs, mp_size_t total_limbs, mp_size_t coeff_limbs, mp_size_t output_limbs)
mp_size_t fft_split_bits(mp_limb_t ** poly, mp_srcptr limbs, mp_size_t total_limbs, flint_bitcnt_t bits, mp_size_t output_limbs)

Split the integer `(limbs, total_limbs)` into consecutive chunks of
`coeff_limbs` limbs, respectively `bits` bits, starting with the least
significant, and write them to `poly[0]`, `poly[1]`, …, each zero-padded
to $\mathit{output\_limbs} + 1$ limbs. Return the number of chunks
written, $\lceil \mathit{total\_limbs}/\mathit{coeff\_limbs} \rceil$,
respectively $\lceil \bits \cdot \mathit{total\_limbs}/\mathit{bits} \rceil$;
later entries of `poly` are not touched. Requires `output_limbs` to be at
least the chunk size in limbs. May use threads.
:::

::: function
void fft_combine_limbs(mp_limb_t * res, mp_limb_t ** poly, slong length, mp_size_t coeff_limbs, mp_size_t output_limbs, mp_size_t total_limbs)
void fft_combine_bits(mp_limb_t * res, mp_limb_t ** poly, slong length, flint_bitcnt_t bits, mp_size_t output_limbs, mp_size_t total_limbs)

Evaluate the polynomial with the `length` coefficients `poly` at
$B^{\mathit{coeff\_limbs}}$, respectively $2^{\mathit{bits}}$, and set
`(res, total_limbs)` to the value, discarding the bits beyond
`total_limbs` limbs. Each coefficient must be less than
$B^{\mathit{output\_limbs}}$ and have $\mathit{output\_limbs} + 1$
readable limbs. `res` must be zero on entry: the coefficients are added
into it, but carries out of the window of each coefficient are lost; see
[](fft/split#combine).
:::

## Arithmetic modulo $2^N + 1$ {#arith}

These functions work on residues modulo $p = 2^N + 1$,
$N = \mathit{limbs} \cdot \bits$, in the representation described under
[](#conventions). All of them allow the output to alias the input.

::: function
void mpn_addmod_2expp1_1(mp_limb_t * r, mp_size_t limbs, mp_limb_signed_t c)

Adds the signed limb $c$ to the $(\mathit{limbs} + 1)$-limb two's
complement integer `r`, in place, without reduction. Called with
`r + k` and `limbs - k`, it adds $c B^k$. Inline; when no carry leaves
the first limb it costs a few instructions.
:::

::: function
void mpn_normmod_2expp1(mp_limb_t * t, mp_size_t limbs)

Reduces the residue `t`, whose top limb is a small signed value, to
$[0, p)$.
:::

::: function
void mpn_negmod_2expp1(mp_limb_t * z, const mp_limb_t * a, mp_size_t limbs)

Sets `z` to $-a \bmod p$ for a reduced residue `a`. The result is
reduced.
:::

::: function
void mpn_mul_2expmod_2expp1(mp_limb_t * t, mp_limb_t * i1, mp_size_t limbs, flint_bitcnt_t d)
void mpn_div_2expmod_2expp1(mp_limb_t * t, mp_limb_t * i1, mp_size_t limbs, flint_bitcnt_t d)

Set `t` to $i_1 \cdot 2^d$, respectively $i_1 / 2^d$, modulo $p$, for
$0 \le d < \bits$ and a residue $i_1$ whose top limb is small. The
result is not necessarily reduced.
:::

::: function
void fermat_to_mpz(mpz_t m, mp_limb_t * i, mp_size_t limbs)

Sets `m` to the signed integer represented by the $\mathit{limbs} + 1$
limbs of `i` in two's complement. The value is *not* reduced modulo $p$:
to obtain the residue in $[0, p)$, reduce `i` with `mpn_normmod_2expp1`
first. Used in tests.
:::

## Butterflies {#butterflies}

Here $p = 2^N + 1$ with $N = \mathit{limbs} \cdot \bits$, and the
parameters $i$ and $w$ specify the twiddle factor $2^{iw}$ or
$\sqrt 2^{\,iw}$. The inputs need not be reduced, and the outputs are
not reduced. Outputs must not alias inputs.

::: function
void fft_adjust(mp_limb_t * r, mp_limb_t * i1, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w)

Sets `r` to $2^{iw} i_1 \bmod p$, for $0 \le iw \le N$.
:::

::: function
void fft_adjust_sqrt2(mp_limb_t * r, mp_limb_t * i1, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp)

Sets `r` to $\sqrt 2^{\,iw} i_1 \bmod p$, for odd $i < 2n$ and odd $w$,
where $nw = N$. `temp` is scratch space of $\mathit{limbs} + 1$ limbs.
:::

::: function
void butterfly_lshB(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1, mp_limb_t * i2, mp_size_t limbs, mp_size_t x, mp_size_t y)

Sets $t = (i_1 + i_2) B^x$ and $u = (i_1 - i_2) B^y$ modulo $p$, for
$0 \le x, y < \mathit{limbs}$. The inputs are not modified.
:::

::: function
void butterfly_rshB(mp_limb_t * t, mp_limb_t * u, mp_limb_t * i1, mp_limb_t * i2, mp_size_t limbs, mp_size_t x, mp_size_t y)

Sets $t = i_1 / B^x + i_2 / B^y$ and $u = i_1 / B^x - i_2 / B^y$ modulo
$p$, for $0 \le x, y < \mathit{limbs}$. The inputs are overwritten when
$x \ne y$.
:::

::: function
void fft_butterfly(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1, mp_limb_t * i2, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w)

The forward butterfly: sets $s = i_1 + i_2$ and
$t = 2^{iw} (i_1 - i_2)$ modulo $p$, for $0 \le iw < N$. The inputs are
not modified.
:::

::: function
void ifft_butterfly(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1, mp_limb_t * i2, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w)

The inverse butterfly: sets $s = i_1 + 2^{-iw} i_2$ and
$t = i_1 - 2^{-iw} i_2$ modulo $p$, for $0 \le iw < N$. `i2` is
overwritten.
:::

::: function
void fft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1, mp_limb_t * i2, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp)

Sets $s = i_1 + i_2$ and $t = \sqrt 2^{\,iw} (i_1 - i_2)$ modulo $p$, for
odd $i < 2n$ and odd $w$, where $nw = N$. `temp` is scratch space of
$\mathit{limbs} + 1$ limbs.
:::

::: function
void ifft_butterfly_sqrt2(mp_limb_t * s, mp_limb_t * t, mp_limb_t * i1, mp_limb_t * i2, mp_size_t i, mp_size_t limbs, flint_bitcnt_t w, mp_limb_t * temp)

Sets $s = i_1 + \sqrt 2^{\,-iw} i_2$ and $t = i_1 - \sqrt 2^{\,-iw} i_2$
modulo $p$, for odd $i < 2n$ and odd $w$, where $nw = N$. `i2` is
overwritten; `temp` is scratch space of $\mathit{limbs} + 1$ limbs.
:::

## Radix-2 and truncated transforms {#transforms}

In this section `ii` holds residues $a_0, \ldots, a_{L-1}$ modulo
$p = 2^{nw} + 1$, $nw$ a multiple of $\bits$, with $L = 2n$, or $L = 4n$
for the sqrt2 functions. With $\omega = 2^w$, respectively
$\sqrt 2^{\,w}$, and $\mathrm{rev}$ the reversal of $\log_2 L$ bits, the
transform is $\hat a_k = \sum_j a_j \omega^{j \, \mathrm{rev}(k)}$, the
values of $a(x) = \sum_j a_j x^j$ at the powers of $\omega$ in
bit-reversed order. `t1`, `t2` and `temp` point to pointers to scratch
blocks of $\mathit{limbs} + 1$ limbs, which are swapped into `ii`.
Results are not reduced. See [](fft/transforms) and
[](algorithms/fft/truncated-fft).

::: function
void fft_radix2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2)
void ifft_radix2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2)

`fft_radix2` replaces $a$ by $\hat a$, the transform of length $2n$ with
root $2^w$. `ifft_radix2` is its inverse up to scaling: it replaces
$\hat a$ by $2n \cdot a$.
:::

::: function
void fft_truncate(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc)
void fft_truncate1(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc)

Truncated FFTs: compute only $\hat a_k$ for $k <$ `trunc`.
`fft_truncate` assumes $a_j = 0$ for $j \ge$ `trunc` and uses those
entries as scratch; `fft_truncate1` reads all $2n$ inputs. Requires
`trunc` even with $2 \le \mathit{trunc} \le 2n$.
:::

::: function
void ifft_truncate(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc)
void ifft_truncate1(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t trunc)

Inverse truncated FFTs. Given $\hat a_k$ for $k <$ `trunc`, they set
entry $j$ to $2n \cdot a_j$ for $j <$ `trunc`. `ifft_truncate` assumes
$a_j = 0$ for $j \ge$ `trunc` and uses those entries as scratch;
`ifft_truncate1` instead requires entry $j$ to hold $2n \cdot a_j$ for
$j \ge$ `trunc` on input, and overwrites it. Requires `trunc` even with
$2 \le \mathit{trunc} \le 2n$.
:::

::: function
void fft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t trunc)
void ifft_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t trunc)

As `fft_truncate` and `ifft_truncate` for the transform of length $4n$
with root $\sqrt 2^{\,w}$ over the same ring; the inverse returns
$4n \cdot a_j$. Requires `trunc` even with $2n < \mathit{trunc} \le 4n$.
:::

::: function
void fft_negacyclic(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp)
void ifft_negacyclic(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp)

The weighted transform of length $2n$ for negacyclic convolution.
`fft_negacyclic` multiplies $a_j$ by $\varphi^j$, $\varphi = \sqrt 2^{\,w}$,
and then applies `fft_radix2`, so that
$\hat a_k = a(\varphi^{2 \mathrm{rev}(k) + 1})$ are the values at the
roots of $x^{2n} + 1$. `ifft_negacyclic` applies `ifft_radix2` and then
multiplies entry $j$ by $\varphi^{-j}$, which gives $2n \cdot a_j$.
:::

## Matrix Fourier transforms {#mfa}

The matrix Fourier algorithm computes a transform of length $2n = n_1 n_2$
as $n_1$ column transforms of length $n_2$, twiddle factors, and $n_2$
row transforms of length $n_1$ [@Bai1990]. Here $n_1$ (`n1`, also called
`sqrt`) is a power of two with $2 \le n_1 \le n$; FLINT always uses
$n_1 = 2^{\lfloor \mathit{depth}/2 \rfloor}$. See
[](algorithms/fft/matrix-fourier).

::: function
void fft_butterfly_twiddle(mp_limb_t * u, mp_limb_t * v, mp_limb_t * s, mp_limb_t * t, mp_size_t limbs, flint_bitcnt_t b1, flint_bitcnt_t b2)
void ifft_butterfly_twiddle(mp_limb_t * u, mp_limb_t * v, mp_limb_t * s, mp_limb_t * t, mp_size_t limbs, flint_bitcnt_t b1, flint_bitcnt_t b2)

`fft_butterfly_twiddle` sets $u = 2^{b_1}(s + t)$ and $v = 2^{b_2}(s - t)$.
`ifft_butterfly_twiddle` sets $u = 2^{-b_1} s + 2^{-b_2} t$ and
$v = 2^{-b_1} s - 2^{-b_2} t$, overwriting `s` and `t`. Both work modulo
$2^N + 1$, $N = \mathit{limbs} \cdot \bits$, for $0 \le b_1, b_2 < 2N$.
They form the last layer of the forward and the first layer of the
inverse column transforms, with the twiddle factors of the matrix
Fourier algorithm merged into the shifts.
:::

::: function
void fft_radix2_twiddle(mp_limb_t ** ii, mp_size_t is, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs)
void ifft_radix2_twiddle(mp_limb_t ** ii, mp_size_t is, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs)

As `fft_radix2` and `ifft_radix2` on the entries `ii[0]`, `ii[is]`, …,
`ii[(2n - 1) is]`, with an extra twiddle: output $k$ of the forward
transform is multiplied by
$2^{\mathit{ws} \cdot c \cdot (r + \mathit{rs} \cdot \mathrm{rev}(k))}$,
and input $k$ of the inverse is divided by the same factor first.
:::

::: function
void fft_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc)
void ifft_truncate1_twiddle(mp_limb_t ** ii, mp_size_t is, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_size_t ws, mp_size_t r, mp_size_t c, mp_size_t rs, mp_size_t trunc)

The truncated versions of `fft_radix2_twiddle` and `ifft_radix2_twiddle`,
in the sense of `fft_truncate1` and `ifft_truncate1`.
:::

::: function
void fft_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
void ifft_mfa_truncate_sqrt2(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)

As `fft_truncate_sqrt2` and `ifft_truncate_sqrt2` (the inverse returns
$4n \cdot a_j$), computed with the matrix Fourier algorithm on each half
of length $2n$. Requires $2n < \mathit{trunc} \le 4n$ with `trunc` a
multiple of $2 n_1$. The transformed values are stored in a different
order than by `fft_truncate_sqrt2`, the order that
`ifft_mfa_truncate_sqrt2` expects. Single-threaded.
:::

::: function
void fft_mfa_truncate_sqrt2_outer(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)
void fft_mfa_truncate_sqrt2_inner(mp_limb_t ** ii, mp_limb_t ** jj, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc, mp_limb_t ** tt)
void ifft_mfa_truncate_sqrt2_outer(mp_limb_t ** ii, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** temp, mp_size_t n1, mp_size_t trunc)

The three stages of a threaded convolution with the matrix Fourier
algorithm, as used by `fft_convolution` and `mul_mfa_truncate_sqrt2`.
`fft_mfa_truncate_sqrt2_outer` performs the sqrt2 layer and the column
transforms of `ii`. `fft_mfa_truncate_sqrt2_inner` performs, one row at
a time, the row transforms of `ii` and `jj`, the pointwise products with
`fft_mulmod_2expp1` and the inverse row transform. `ifft_mfa_truncate_sqrt2_outer`
performs the inverse column transforms and the sqrt2 layer and divides
by $4n$, leaving reduced residues. The requirements are those of
`fft_mfa_truncate_sqrt2` and, for the pointwise products, those of
`fft_mulmod_2expp1`. The scratch arrays need one entry per thread;
`temp` is not used by the `_inner` function.
:::

## Internal functions {#internal}

::: function
mp_limb_t flint_mpn_sumdiff_n(mp_ptr s, mp_ptr d, mp_srcptr x, mp_srcptr y, mp_size_t n)

Sets `(s, n)` to $x + y$ and `(d, n)` to $x - y$, each modulo $B^n$, and
returns $2c + b$, where $c$ is the carry of the sum and $b$ the borrow of
the difference. The outputs may alias the inputs. Declared here for the
butterflies; defined in `mpn_extras`.
:::

::: macro
fft_sumdiff(t, u, r, s, n)

`flint_mpn_sumdiff_n(t, u, r, s, n)`, but 0 without any work if
$n = 0$.
:::

::: macro
SWAP_PTRS(xx, yy)

Exchanges the two `mp_ptr` variables `xx` and `yy`.
:::

::: macro
random_fermat(nn, state, limbs)

Sets the residue `nn` of $\mathit{limbs} + 1$ limbs to a random value for
testing, including the extreme value $2^N$ and unreduced values with top
limbs up to $\pm 1023$. Needs `mpn_extras.h`.
:::

::: function
void fft_convolution_basic(mp_limb_t ** ii, mp_limb_t ** jj, slong depth, slong limbs, slong trunc, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1, mp_limb_t ** tt)

As `fft_convolution`, but computing the transforms and the pointwise
products in separate single-threaded passes instead of the fused
threaded stages. Used to test `fft_convolution`.
:::

::: function
void fft_mfa_truncate_sqrt2_inner_sd_fft(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx, mp_limb_t ** ii, mp_limb_t ** jj, mp_size_t n, flint_bitcnt_t w, mp_limb_t ** t1, mp_limb_t ** t2, mp_limb_t ** s1, mp_size_t sqrt, mp_size_t trunc, mp_limb_t ** tt)

As `fft_mfa_truncate_sqrt2_inner`, with the pointwise products computed
by `sd_fft_mpn_mulmod_2expp1` with the context `sdctx`, initialised for
$N = nw$; each thread gets its own scratch for the context. Only defined
in builds with `fft_small`.
:::

::: function
void _fft_pointwise_sd_fft(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx, mp_limb_t ** ii, mp_limb_t ** jj, slong j0, slong j1, slong limbs)
void _fft_pointwise_sd_fft_ro(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx, mp_limb_t ** ii, mp_limb_t ** jj, slong j0, slong j1, slong limbs)

For $j_0 \le j < j_1$, reduce `ii[j]` (and, in the first function,
`jj[j]`) and set `ii[j]` to the product of `ii[j]` and `jj[j]` modulo
$2^N + 1$ with `sd_fft_mpn_mulmod_2expp1`, using threads. The `_ro`
version does not modify `jj`, which must already be reduced. Without
`fft_small` they do nothing; they are then never called.
:::

::: function
void _fft_pointwise_rows_sd_fft(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx, mp_limb_t ** ii, mp_limb_t ** jj, slong base, slong srt, slong nrows, slong rvbits, slong limbs)
void _fft_pointwise_rows_sd_fft_ro(struct sd_fft_mpn_mulmod_2expp1_ctx_struct * sdctx, mp_limb_t ** ii, mp_limb_t ** jj, slong base, slong srt, slong nrows, slong rvbits, slong limbs)

As `_fft_pointwise_sd_fft` and `_fft_pointwise_sd_fft_ro`, for the
indices $j = \mathit{base} + \mathrm{rev}(s) \cdot \mathit{srt} + u$ with
$0 \le s < \mathit{nrows}$ and $0 \le u < \mathit{srt}$, where
$\mathrm{rev}$ reverses `rvbits` bits: the rows of the truncated half of
a matrix Fourier transform.
:::
