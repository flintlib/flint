---
title: Truncated FFTs
summary: Van der Hoeven's truncated Fourier transform evaluates a polynomial of length $t$ at the first $t$ points of an FFT of length $L = 2^\ell \ge t$ in about $(t/2) \log_2 L$ butterflies, and its inverse recovers the polynomial from these values, which removes the jumps of radix-2 convolution at powers of two.
complexity: At most $(t/2)\log_2 L + L$ butterflies (or additions) for $t$ of $L$ points, forward and inverse.
order: 20
---

A convolution with radix-2 FFTs uses a transform length $L = 2^\ell$ at
least as large as the length $t$ of the product. Rounding $t$ up to a
power of two can double the work: a product of length $2^{10} + 1$ costs
as much as one of length $2^{11}$. (A middle product can use a shorter
length, since the wrap-around may fall outside the coefficients it needs;
see [](algorithms/polynomials/middle-product#fft).) The *truncated Fourier transform*
(TFT) of van der Hoeven [@vdH2004; @vdH2005] computes only $t$ of the $L$
values, and its inverse recovers a polynomial of length $t$ from them, in
time proportional to $t \log L$ rather than $L \log L$. Harvey made it
cache-friendly [@Har2009], Harvey and Roche [@HR2010] and Arnold
[@Arn2013] gave in-place variants.

FLINT uses truncated transforms in both of its integer FFTs: `fft.h`
implements van der Hoeven's recursion in `fft_truncate` and
`ifft_truncate` (and their sqrt2 and matrix Fourier variants), and
`fft_small` implements Harvey's version with separate input and output
truncation in `sd_fft_trunc` and `sd_ifft_trunc`, and their transposes
in `sd_fft_trunc_tr` and `sd_ifft_trunc_tr`
([](algorithms/fft/small-prime-fft#transposed)). The algorithms are
stated below in the form used by `fft.h`, with proofs.

## Notation {#notation}

Let $R$ be a commutative ring, $L = 2^\ell$, $m = L/2$, and let
$\omega \in R$ satisfy $\omega^m = -1$, with $2$ invertible in $R$. (In
`fft.h`, $R = \Z/(2^N + 1)\Z$ and $\omega = 2^w$ or $\sqrt 2^{\,w}$; see
[](algorithms/fft/schonhage-strassen#ring).) Let $\mathrm{rev}_L(k)$ be
the integer whose $\ell$ binary digits are those of $k$ reversed. For
$a = (a_0, \ldots, a_{L-1})$ and $a(x) = \sum_j a_j x^j$, define

$$
\hat a_k = a\bigl(\omega^{\mathrm{rev}_L(k)}\bigr), \qquad 0 \le k < L,
$$

the transform in bit-reversed order that `fft_radix2` computes by
radix-2 decimation in frequency
([](algorithms/fft/dft-algorithms#alg-dif)). The
values at positions below $t$ are the *first $t$ values*.

::: lemma Splitting {#lem-split}
For $0 \le i < m$ let $b_i = a_i + a_{m+i}$ and
$c_i = \omega^i (a_i - a_{m+i})$. With the transforms of length $m$ taken
with respect to $\omega^2$,
$$
\hat a_k = \hat b_k \quad \text{and} \quad \hat a_{m+k} = \hat c_k
\qquad (0 \le k < m).
$$
:::

::: proof
For $k < m$, $\mathrm{rev}_L(k) = 2\,\mathrm{rev}_m(k)$ and
$\mathrm{rev}_L(m + k) = 2\,\mathrm{rev}_m(k) + 1$. Let $e = \mathrm{rev}_m(k)$.
Since $\omega^{2em} = 1$, $a(\omega^{2e}) = \sum_{i<m} (a_i + a_{m+i}) \omega^{2ei} = b(\omega^{2e})$.
Since $\omega^{(2e+1)m} = -1$,
$a(\omega^{2e+1}) = \sum_{i<m} (a_i - a_{m+i}) \omega^{i(2e+1)} = \sum_{i<m} c_i (\omega^2)^{ei} = c(\omega^{2e})$.
:::

The splitting is one layer of decimation-in-frequency butterflies
[@GS1966]; applied recursively it gives `fft_radix2`. Truncation consists
in skipping all work whose results are not among the first $t$ values.

## The forward transform {#forward}

::: algorithm TFT with full input (`fft_truncate1`) {#alg-tft1}
**Input:** $a_0, \ldots, a_{L-1}$; an even $t$ with $2 \le t \le L$.

**Output:** $\hat a_0, \ldots, \hat a_{t-1}$.

1. If $t = L$, compute the full transform.
2. If $t \le m$: set $b_i \gets a_i + a_{m+i}$ for $i < m$ and return the
   first $t$ values of the transform of $b$ (recursively, with
   $\omega^2$).
3. Otherwise: compute $b$ and $c$ by $m$ butterflies, compute the full
   transform of $b$, and return it followed by the first $t - m$ values
   of the transform of $c$ (recursively).
:::

::: algorithm TFT with zero padding (`fft_truncate`) {#alg-tft}
**Input:** $a_0, \ldots, a_{t-1}$ (and $a_j = 0$ for $j \ge t$); an even
$t$ with $2 \le t \le L$.

**Output:** $\hat a_0, \ldots, \hat a_{t-1}$.

1. If $t = L$, compute the full transform.
2. If $t \le m$: then $b_i = a_i$, so return the first $t$ values of the
   transform of $(a_0, \ldots, a_{m-1})$, recursively with zero padding.
3. Otherwise: for $i < t - m$ compute $b_i, c_i$ by a butterfly; for
   $t - m \le i < m$, where $a_{m+i} = 0$, set $b_i \gets a_i$ and
   $c_i \gets \omega^i a_i$. Compute the full transform of $b$ and the
   first $t - m$ values of that of $c$ by [](#alg-tft1).
:::

Both are correct by [](#lem-split) and induction on $L$. In step 3 of
[](#alg-tft) the vector $c$ is in general not zero-padded, which is why
the recursion switches to [](#alg-tft1).

::: theorem Cost {#thm-cost}
[](#alg-tft1) and [](#alg-tft) use at most $(t/2) \log_2 L + L$
butterflies, additions or multiplications by powers of $\omega$.
:::

::: proof
Let $C(L, t)$ be the worst count, and $C(L, L) = (L/2)\log_2 L$ for the
full transform. If $t \le m$, the count is at most $m + C(m, t)$. If
$t > m$, it is at most $m + (m/2)\log_2 m + C(m, t - m)$. We show
$C(L, t) \le (t/2)\log_2 L + L$ by induction. The case $t = L$ is clear.
If $t \le m$, by induction
$C(L, t) \le m + (t/2)(\log_2 L - 1) + m \le (t/2)\log_2 L + L$. If
$t > m$, $C(L, t) \le m + (m/2)(\log_2 L - 1) + ((t - m)/2)(\log_2 L - 1) + m = (t/2)(\log_2 L - 1) + L$.
:::

The full transform costs $(L/2)\log_2 L$, so for $t$ just above $L/2$
the truncated transform costs about half. In practice the additive term
$L$ is smaller than the bound suggests, since [](#alg-tft) does not add
anything in step 2.

## The inverse transform {#inverse}

The inverse must recover $a_0, \ldots, a_{t-1}$ from
$\hat a_0, \ldots, \hat a_{t-1}$ when $a_j = 0$ for $j \ge t$. Since
the evaluation points $\omega^{\mathrm{rev}(k)}$, $k < t$, do not form a
subgroup, this is an interpolation problem, but the splitting lemma
reduces it to the same recursion if, at each level, the missing values
of one half are reconstructed from the other. The general step works
with a *mixed* vector that holds transform values at the positions below
$t$ and (scaled) coefficients at the positions from $t$ on.

::: algorithm Inverse TFT with mixed input (`ifft_truncate1`) {#alg-itft1}
**Input:** $x_k = \hat a_k$ for $k < t$ and $x_j = L a_j$ for
$t \le j < L$; an even $t$ with $2 \le t \le L$.

**Output:** $x_j = L a_j$ for $j < t$.

1. If $t = L$, apply the full inverse transform, which gives $L a$.
2. If $t \le m$:
   1. for $t \le i < m$: $x_i \gets (x_i + x_{m+i})/2$;
   2. apply the algorithm recursively to $x_0, \ldots, x_{m-1}$ with $t$
      (length $m$, root $\omega^2$);
   3. for $i < t$: $x_i \gets 2 x_i - x_{m+i}$.
3. If $t > m$:
   1. apply the full inverse transform of length $m$ to
      $x_0, \ldots, x_{m-1}$;
   2. for $t - m \le i < m$: $d \gets x_i - x_{m+i}$,
      $x_{m+i} \gets \omega^i d$, $x_i \gets x_i + d$;
   3. apply the algorithm recursively to $x_m, \ldots, x_{L-1}$ with
      $t - m$;
   4. for $i < t - m$:
      $(x_i, x_{m+i}) \gets (x_i + \omega^{-i} x_{m+i}, x_i - \omega^{-i} x_{m+i})$.
:::

::: theorem {#thm-itft1}
[](#alg-itft1) is correct: on return, $x_j = L a_j$ for $j < t$.
:::

::: proof
Induction on $L$; let $b$, $c$ be as in [](#lem-split), and note that the
inverse transform of length $m$ returns $m$ times its argument.

Case $t \le m$. After step 2.1, for $t \le i < m$,
$x_i = (L a_i + L a_{m+i})/2 = m b_i$, while for $k < t$,
$x_k = \hat a_k = \hat b_k$. This is valid input for length $m$, so after
step 2.2, $x_i = m b_i$ for $i < t$. Since $m + i \ge t$, $x_{m+i} = L a_{m+i}$
was given, and step 2.3 yields $2m b_i - L a_{m+i} = L a_i$.

Case $t > m$. The first $m$ values are $\hat b_0, \ldots, \hat b_{m-1}$,
so step 3.1 gives $x_i = m b_i$ for all $i < m$. For $t - m \le i < m$,
the position $m + i \ge t$ holds $L a_{m+i} = 2m a_{m+i}$, so
$d = m(a_i + a_{m+i}) - 2m a_{m+i} = m(a_i - a_{m+i})$; step 3.2 sets
$x_{m+i} = m c_i$ and $x_i = 2m a_i = L a_i$, which is final. Now the
second half holds $\hat c_k = \hat a_{m+k}$ at its positions $k < t - m$
and $m c_i$ at its positions $i \ge t - m$, valid input for length $m$,
so after step 3.3 it holds $m c_i$ for $i < t - m$. Finally, since
$\omega^{-i} c_i = a_i - a_{m+i}$, step 3.4 turns $(m b_i, m c_i)$ into
$(2m a_i, 2m a_{m+i})$.
:::

::: algorithm Inverse TFT with zero padding (`ifft_truncate`) {#alg-itft}
**Input:** $\hat a_0, \ldots, \hat a_{t-1}$, where $a_j = 0$ for $j \ge t$;
an even $t$ with $2 \le t \le L$.

**Output:** $L a_0, \ldots, L a_{t-1}$.

1. If $t = L$, apply the full inverse transform.
2. If $t \le m$: then $b = (a_0, \ldots, a_{m-1})$ with zero padding;
   recursively obtain $m a_i$ for $i < t$ and double them.
3. If $t > m$: apply the full inverse transform of length $m$, giving
   $m b_i$; for $t - m \le i < m$, where $a_{m+i} = 0$ and so
   $c_i = \omega^i b_i$, set $x_{m+i} \gets \omega^i x_i$; apply
   [](#alg-itft1) to the second half with $t - m$; apply the inverse
   butterflies for $i < t - m$ and double $x_i$ for $t - m \le i < m$.
:::

::: theorem Inversion {#thm-inverse}
If $a_j = 0$ for $j \ge t$, then $\hat a_0, \ldots, \hat a_{t-1}$
determine $a$, and [](#alg-itft) computes $L a_0, \ldots, L a_{t-1}$ from
them. Both inverse algorithms use at most $(t/2)\log_2 L + L$ butterflies,
additions, halvings and multiplications by powers of $\omega$.
:::

::: proof
[](#alg-itft) is [](#alg-itft1) with the known values $L a_j = 0$ for
$j \ge t$ substituted: in step 2, $x_i = (0 + 0)/2$ makes the first loop
trivial and $2 x_i - 0$ the last one a doubling; in step 3, $d = x_i$, so
$x_{m+i} = \omega^i x_i$ and the final value $x_i + d$ is a doubling,
which can be done at the end. Correctness follows from
[](#thm-itft1), and the algorithm shows that the $t$ values determine $a$,
since $L$ is invertible. For the cost, the loops of step 2 take $m$
operations in total, and those of step 3 take $(2m - t) + (t - m) = m$
operations besides the transform of length $m$, so the count satisfies
the recursion of [](#thm-cost).
:::

To multiply polynomials $f$ and $g$ whose product has length at most
$t$, one therefore computes the first $t$ values of $\hat f$ and
$\hat g$, multiplies them pointwise (these are the first $t$ values of
the transform of $fg$, because $\widehat{fg}_k = (fg)(\omega^{\mathrm{rev}(k)}) = \hat f_k \hat g_k$),
and applies the inverse TFT. This is `fft_convolution`: the product
costs about $3 \cdot (t/2) \log_2 L$ butterflies instead of
$3 \cdot (L/2) \log_2 L$, and its running time grows smoothly with $t$.

## Use in FLINT {#flint}

- `fft_truncate`, `fft_truncate1`, `ifft_truncate` and `ifft_truncate1`
  are [](#alg-tft), [](#alg-tft1), [](#alg-itft) and [](#alg-itft1) with
  $\omega = 2^w$ modulo $2^{nw} + 1$; the halving in step 2.1 of
  [](#alg-itft1) is `mpn_div_2expmod_2expp1` with $d = 1$. FLINT
  requires $t$ to be even, since its recursion bottoms out at length 2
  with $t = 2$; the callers round $t$ up.
- `fft_truncate_sqrt2` and `ifft_truncate_sqrt2` use length $4n$ with
  $\omega = \sqrt 2^{\,w}$ and require $t > 2n$: they always perform the
  first layer in the form of step 3 (with the
  [$\sqrt 2$ twiddles](algorithms/fft/schonhage-strassen#sqrt2)), a full
  transform of length $2n$ on the first half and a truncated one on the
  second.
- The matrix Fourier functions truncate at whole rows: the number of
  points is rounded up to a multiple of $2 n_1$, where $n_1$ is the row
  length, so the row transforms are complete and only the column
  transforms of the second half are truncated, by
  `fft_truncate1_twiddle` and `ifft_truncate1_twiddle`; see
  [](algorithms/fft/matrix-fourier#truncation).
- `fft_small` uses the more general form of Harvey [@Har2009]:
  `sd_fft_trunc` takes separate truncations for the input (coefficients
  known to be zero from `itrunc` on) and the output (values wanted below
  `otrunc`), and the inverse `sd_ifft_trunc` is a recursion in terms of
  the number $z$ of possibly nonzero coefficients, the number $n$ of known
  values and a flag $f$ asking for one more value, which [](#alg-itft1)
  specialises to $z = n = t$ and $f = 0$. There the truncation works in
  blocks of 256 points and is combined with a decomposition into rows and
  columns as in the [matrix Fourier algorithm](algorithms/fft/matrix-fourier).
  See [](algorithms/fft/small-prime-fft).
