---
title: The matrix Fourier algorithm
summary: Bailey's four-step algorithm computes an FFT of length $n_1 n_2$ as $n_1$ transforms of length $n_2$, twiddle factors and $n_2$ transforms of length $n_1$, so that each partial transform fits in cache; FLINT fuses it with the pointwise products of a convolution, truncates at whole rows and runs the partial transforms in threads.
complexity: The same number of butterflies as the radix-2 FFT; about two passes over the data per transform instead of $\log_2(L/C)$ for a cache of $C$ coefficients.
order: 30
---

A radix-2 FFT of length $L$ performs $\log_2 L$ layers of butterflies,
and each layer touches every coefficient. Once the data no longer fits in
cache, the top layers stream the whole vector through memory one after
another, and memory bandwidth rather than arithmetic determines the
running time. The *matrix Fourier algorithm*, also called the four-step
or six-step FFT, regroups the work into short transforms that each fit in
cache. It goes back to Gentleman and Sande [@GS1966] and was popularised
for hierarchical memory by Bailey [@Bai1990]; Arndt describes it under
this name [@Arn2010]. FLINT's `fft` module uses it for transforms of
length $4n$ with $n > 64$ (the functions `fft_mfa_truncate_sqrt2`,
`ifft_mfa_truncate_sqrt2` and their `_outer` and `_inner` parts), that
is, for `fft_convolution` with $\mathit{depth} > 6$ and for
`flint_mpn_mul_fft_main` from a total size of about 65000 limbs on
(32600 limbs per operand for balanced products).

## The decomposition {#decomposition}

Let $L = n_1 n_2$ and let $\omega$ be a root of unity of order $L$. View
the coefficients as a matrix with $n_2$ rows and $n_1$ columns stored by
rows: $a_j$ with $j = j_1 + n_1 j_2$ sits in row $j_2$ and column $j_1$.

::: theorem Four-step formula {#thm-four-step}
For $0 \le k_1 < n_1$ and $0 \le k_2 < n_2$,
$$
a\bigl(\omega^{k_2 + n_2 k_1}\bigr)
= \sum_{j_1 = 0}^{n_1 - 1} (\omega^{n_2})^{j_1 k_1} \,
  \omega^{j_1 k_2}
  \sum_{j_2 = 0}^{n_2 - 1} a_{j_1 + n_1 j_2} \, (\omega^{n_1})^{j_2 k_2}.
$$
:::

::: proof
Expand $\omega^{(j_1 + n_1 j_2)(k_2 + n_2 k_1)} = \omega^{j_1 k_2} \, \omega^{n_2 j_1 k_1} \, \omega^{n_1 j_2 k_2} \, \omega^{n_1 n_2 j_2 k_1}$
and use $\omega^{n_1 n_2} = 1$.
:::

The inner sum is a transform of length $n_2$ of column $j_1$ (root
$\omega^{n_1}$), evaluated at $k_2$. It is multiplied by the *twiddle
factor* $\omega^{j_1 k_2}$, and the outer sum is a transform of length
$n_1$ along row $k_2$ (root $\omega^{n_2}$). This gives the algorithm:

::: algorithm Matrix Fourier transform {#alg-mfa}
**Input:** $a_0, \ldots, a_{L-1}$ as an $n_2 \times n_1$ matrix by rows.

**Output:** the values $a(\omega^{k_2 + n_2 k_1})$, stored in row $k_2$,
column $k_1$.

1. For each column $j_1$: transform it with length $n_2$ and root
   $\omega^{n_1}$, and multiply its entry in row $k_2$ by
   $\omega^{j_1 k_2}$.
2. For each row $k_2$: transform it with length $n_1$ and root
   $\omega^{n_2}$.
:::

Bailey's six-step version adds transpositions so that the values come out
in natural order. A convolution does not need them: the pointwise
products do not care about the order, and the inverse transform performs
the steps in reverse (inverse row transforms, division by the twiddle
factors, inverse column transforms), which undoes any order. The
`gr_dft` module implements the same decomposition over general rings
([](algorithms/fft/dft-algorithms#four-step)).

In `fft.h`, $L = 2n$ is one half of a sqrt2 transform, $\omega = 2^w$,
$n_1 = 2^{\lfloor \mathit{depth}/2 \rfloor}$ and $n_2 = 2n/n_1$, so the
column transforms are slightly longer than the rows. The column transform
`fft_radix2_twiddle` produces its values in bit-reversed order; the
twiddle factor for position $k$ is then $\omega^{j_1 \mathrm{rev}(k)}$,
which is exactly what its parameters $(\mathit{ws}, r, c, \mathit{rs}) = (w, 0, j_1, 1)$
specify, and it is merged into the shifts of the last layer of
butterflies (`fft_butterfly_twiddle`), so it costs nothing extra. After
each column and row transform the entries are permuted from bit-reversed
to natural order by swapping pointers, so that the final layout is that
of [](#alg-mfa). The forward and inverse transforms of the two halves of
the sqrt2 transform are computed this way, with the first sqrt2 layer of
butterflies done column by column just before the column transforms of
the first half.

## Cache behaviour {#cache}

The coefficients are blocks of $N/\bits + 1$ limbs reached through an
array of pointers, so a column with stride $n_1$ is a set of separate
blocks, not a strided walk through one array. Suppose the cache holds
$C$ coefficients. A recursive radix-2 transform of length $L > C$ works
in cache once its subtransforms have length at most $C$, but each of the
$\log_2(L/C)$ layers above has to stream all $L$ coefficients between
memory and cache. [](#alg-mfa) needs only that a column of $n_2$ and a
row of $n_1$ coefficients fit, about $\sqrt{2L}$ and $\sqrt{L/2}$
coefficients: then step 1 and step 2 each move every coefficient once, so
the transform makes two passes over memory whatever its length.

For example, `flint_mpn_mul_fft_main` multiplies two integers of
$10^6$ limbs with $\mathit{depth} = 13$ and $N = 8192$ bits: each half
of each transform has $2n = 16384$ coefficients of about 1 KB, that is,
16 MB. The columns have $n_2 = 256$ coefficients (256 KB, which fits in
a typical L2 cache) and the rows $n_1 = 64$ (64 KB). A radix-2 transform
would stream the 16 MB four to six times, depending on the cache size.

## Fusion with the convolution {#fusion}

A convolution needs the forward transforms of both operands, the
pointwise products and the inverse transform. `fft_convolution` and
`mul_mfa_truncate_sqrt2` organise them in three stages so that the data
of each operand passes through the cache three times in all:

1. `fft_mfa_truncate_sqrt2_outer`: the sqrt2 layer and the column
   transforms of both halves, for each operand.
2. `fft_mfa_truncate_sqrt2_inner`: for each row, the row transforms of
   both operands, the pointwise products with `fft_mulmod_2expp1` (or
   `sd_fft_mpn_mulmod_2expp1`), and the inverse row transform, while the
   row is in cache.
3. `ifft_mfa_truncate_sqrt2_outer`: the inverse column transforms, the
   inverse sqrt2 layer and the division by $4n$, with the final reduction
   of each coefficient done while it is in cache.

`fft_convolution_sd_fft` departs from this when the rows are few compared
with the threads; see [](fft/convolution#algorithm). The precached
convolutions and `fft_convolution_basic` instead call the complete
transforms `fft_mfa_truncate_sqrt2` and `ifft_mfa_truncate_sqrt2` and
multiply in a separate pass.

## Truncation {#truncation}

The matrix layout combines with [truncation](algorithms/fft/truncated-fft)
at the granularity of rows. FLINT requires the number of points $t$ of
the sqrt2 transform to satisfy $2n < t \le 4n$, so the first half is
always complete, and rounds $t$ up to a multiple of $2 n_1$. In the
second half, the
$t_2 = (t - 2n)/n_1$ needed rows are the first $t_2$ outputs of the
column transforms in bit-reversed order, so the column transforms are
truncated transforms with full input (`fft_truncate1_twiddle`, since
the sqrt2 layer has filled the whole second half), and only the rows
$\mathrm{rev}(s)$, $0 \le s < t_2$, are transformed, multiplied and
transformed back. The inverse uses `ifft_truncate1_twiddle` on the
columns. Requiring $t_2$ to be even keeps the truncated column transforms
within FLINT's convention. The rounding costs at most $2 n_1 - 1 \approx 2\sqrt{2n}$
extra points.

## Threads {#threads}

All column transforms are independent, and so are all rows in the inner
stage. The `_outer` and `_inner` functions start up to
`flint_get_num_threads()` threads of the global thread pool, which take
blocks of 16 columns or rows from a shared counter protected by a mutex,
so faster threads take more blocks. Each thread needs its own scratch
coefficients, which is why these functions take arrays `t1`, `t2`,
`temp` and `tt` with one entry per thread. With the `fft_small` pointwise
engine, whose products dominate the work of a row, the inner stage
requests one thread per column or row rather than one per block of 16.

`fft_small` applies the same idea with a different granularity: its
transforms of length $2^k$ split recursively into $2^{k_1} \times 2^{k_2}$
blocks with $k_1 = \lfloor k/2 \rfloor$, on vectors of 256 doubles, and
truncate in units of such blocks; see [](algorithms/fft/small-prime-fft).
