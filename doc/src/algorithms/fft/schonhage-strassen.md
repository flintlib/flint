---
title: Schönhage–Strassen multiplication
summary: Multiplies large integers with FFTs over $\Z/(2^N + 1)\Z$, where the roots of unity are powers of 2 and every twiddle factor is a shift; FLINT adds the $\sqrt 2$ trick, truncation, negacyclic pointwise products and tuned parameters.
complexity: $O(n \log n \log \log n)$ bit operations for $n$-bit integers with full recursion; FLINT recurses one level, and its transforms cost $O(n \log n)$ bit operations.
order: 10
---

Schönhage and Strassen [@SS1971] multiply integers with a discrete
Fourier transform over the ring $\Z/(2^N + 1)\Z$. In this ring $2$ has
order $2N$, so the roots of unity needed by the FFT can be chosen as
powers of two, and multiplying by them is a shift with wraparound,
costing $O(N)$ bit operations instead of a multiplication. The pointwise
products are products in the same kind of ring and can be computed with
the same method, which gives $O(n \log n \log \log n)$ bit operations for
$n$-bit integers. Fürer [@Fur2009] and Harvey and van der Hoeven
[@HvdH2021] improved the asymptotic bound to $O(n \log n)$, but the
Schönhage–Strassen algorithm remains the standard exact integer FFT, as
in GMP [@GKZ2007] and in FLINT's `fft` module, written by William Hart.

The transforms of `fft.h` are truncated
([](algorithms/fft/truncated-fft)) and, for large lengths, computed with
the matrix Fourier algorithm ([](algorithms/fft/matrix-fourier)).

## The ring {#ring}

Let $n = 2^{\mathit{depth}}$, let $w \ge 1$, and let $N = nw$ and
$p = 2^N + 1$. Put $\theta = 2^w \in \Z/p\Z$.

::: lemma Roots of unity {#lem-root}
$\theta$ has order exactly $2n$, and for every integer $k$ that is not a
multiple of $2n$,
$$
\sum_{j=0}^{2n-1} \theta^{jk} = 0 \quad \text{in } \Z/p\Z.
$$
Moreover, $2$ is a unit in $\Z/p\Z$.
:::

::: proof
$\theta^n = 2^N = -1$, so $\theta^{2n} = 1$ and the order divides $2n$
but not $n$; as $2n$ is a power of two, the order is $2n$. For the sum,
we may reduce $k$ modulo $2n = 2^t$ and assume $0 < k < 2n$. Write
$k = 2^s k'$ with $k'$ odd, so $s < t$, and let $m = 2^{t - 1 - s}$.
Then $\theta^{mk} = \theta^{2^{t-1} k'} = \theta^{n k'} = (-1)^{k'} = -1$.
The indices $0, \ldots, 2n - 1$ fall into blocks of length $2m$, which
divides $2n$, and within a block the terms $j$ and $j + m$ cancel, since
$\theta^{(j+m)k} = -\theta^{jk}$. Finally, $2 \cdot 2^{2N - 1} = 2^{2N} = 1$.
:::

The lemma says that $\theta$ is a *principal* root of unity, which is
what the inversion formula needs in a ring that is not a field
[@Pol1971]: for $\hat a_k = \sum_j a_j \theta^{jk}$,

$$
\sum_{k=0}^{2n-1} \theta^{-ik} \hat a_k
= \sum_j a_j \sum_k \theta^{(j - i)k} = 2n \, a_i,
$$

and $2n$ is invertible. Hence pointwise products of transforms are the
transforms of cyclic convolutions of length $2n$, and the inverse
transform recovers these after a division by $2n$, which is a shift.

Multiplication by $\theta^i = 2^{iw}$ is a shift by $iw$ bits. Since
$2^{N + s} \equiv -2^s$, the bits shifted out at the top re-enter at the
bottom with a negative sign. FLINT splits the shift into whole limbs,
realised by reading or writing at an offset with the wrapped part
subtracted, and a remaining shift by fewer than $\bits$ bits. Residues
are kept in $N/\bits + 1$ limbs in two's complement, and sums and
differences are not reduced: the top limb absorbs the carries, and
`mpn_normmod_2expp1` reduces only where a reduced value is needed. A
butterfly $(a, b) \mapsto (a + b, 2^{iw}(a - b))$ is thus a combined
addition and subtraction (`flint_mpn_sumdiff_n`) followed by a shift: a
few passes over $N$ bits. With these butterflies the transforms are
radix-2 Cooley–Tukey FFTs ([](algorithms/fft/dft-algorithms#radix2)).

## The $\sqrt 2$ trick {#sqrt2}

A transform of length $2n$ with root $2^w$ needs a ring of $N = nw \ge n$
bits. The following identity, also used by GMP [@GKZ2007], gives a root
of unity of order $4n$ in the same ring.

::: lemma {#lem-sqrt2}
If $4 \mid N$, then $\sigma = 2^{3N/4} - 2^{N/4}$ satisfies
$\sigma^2 \equiv 2 \pmod{2^N + 1}$. Hence $\varphi = \sigma^w$ has
order $4n$ and $\varphi^2 = \theta$.
:::

::: proof
$\sigma^2 = 2^{3N/2} - 2 \cdot 2^N + 2^{N/2}$. Modulo $2^N + 1$ we have
$2^{3N/2} = 2^N 2^{N/2} \equiv -2^{N/2}$ and $2 \cdot 2^N \equiv -2$, so
$\sigma^2 \equiv -2^{N/2} + 2 + 2^{N/2} = 2$. Then
$\varphi^2 = \sigma^{2w} = 2^w = \theta$ and $\varphi^{2n} = \theta^n = -1$,
so $\varphi$ has order $4n$.
:::

Even powers $\varphi^{2j} = 2^{jw}$ are shifts. For odd $i$ and odd $w$,
$iw$ is odd and $\varphi^i = \sigma^{iw} = \sigma \cdot 2^{(iw - 1)/2}$.
Since $\sigma = 2^{N/4} (2^{N/2} - 1)$, multiplying by $\varphi^i$ costs a
shift by $(iw - 1)/2 + N/4$ bits, a shift by $N/2$ bits and a
subtraction; this is `fft_adjust_sqrt2` and `fft_butterfly_sqrt2`. For
even $w$, $\varphi = 2^{w/2}$ and nothing special is needed. In a
transform of length $4n$ by decimation in frequency, only the first layer
uses odd powers of $\varphi$; the two halves are transforms of length
$2n$ with root $\theta$. The trick therefore doubles the transform length
for a given ring, or halves the ring for a given number of points, at
the cost of a slightly more expensive first layer. All `_sqrt2` functions
of `fft.h` use it, with $N$ a multiple of $\bits$, hence of 4.

## Integer multiplication {#integers}

To multiply $x$ and $y$, FLINT cuts them into chunks of $b$ bits,
$x = \sum_{i < j_1} x_i 2^{bi}$ and $y = \sum_{i < j_2} y_i 2^{bi}$ with
$0 \le x_i, y_i < 2^b$ (`fft_split_bits`). The product $xy$ is the value
at $2^b$ of the polynomial $\sum_k c_k t^k$ with
$c_k = \sum_{i + j = k} x_i y_j$, which is computed by a convolution of
length $4n$ modulo $2^N + 1$ using the $\sqrt 2$ trick.

::: algorithm Integer multiplication (`mul_truncate_sqrt2`) {#alg-mul}
**Input:** $x$, $y$; parameters $\mathit{depth}$ and $w$ with
$N = nw$ a multiple of $\bits$, $n = 2^{\mathit{depth}}$.

**Output:** $xy$.

1. $b \gets \lfloor (N - \mathit{depth} - 1)/2 \rfloor$. Split $x$ and $y$
   into $j_1$ and $j_2$ chunks of $b$ bits; require
   $j_1 + j_2 - 1 \le 4n$.
2. $t \gets j_1 + j_2 - 1$, raised to $2n + 1$ if smaller, and rounded up
   to an even number.
3. Compute the truncated sqrt2 transforms of length $4n$ of both chunk
   vectors at the first $t$ points (`fft_truncate_sqrt2`).
4. Reduce the $t$ transformed values and multiply them pairwise modulo
   $2^N + 1$.
5. Apply the truncated inverse transform (`ifft_truncate_sqrt2`), divide
   by $4n$ and reduce; this gives $c_0, \ldots, c_{t - 1}$.
6. Return $\sum_k c_k 2^{bk}$ (`fft_combine_bits`).
:::

::: theorem {#thm-mul}
[](#alg-mul) returns $xy$.
:::

::: proof
Each $c_k$ is a sum of at most $\min(j_1, j_2)$ products $x_i y_j$. From
$j_1 + j_2 \le 4n + 1$ we get $\min(j_1, j_2) \le 2n$, so
$$
0 \le c_k \le 2n (2^b - 1)^2 < 2^{\mathit{depth} + 1 + 2b} \le 2^N < p,
$$
using $2b \le N - \mathit{depth} - 1$. Since the product polynomial has
$j_1 + j_2 - 1 \le t \le 4n$ coefficients, it equals the cyclic product
of length $4n$, and by [](algorithms/fft/truncated-fft#thm-inverse) the
truncated transforms at $t$ points determine it. The residues in
$[0, p)$ found in step 5 are therefore the integers $c_k$, and step 6
evaluates the product polynomial at $2^b$.
:::

The bound $2^{\mathit{depth} + 1 + 2b} \le 2^N$ means that almost half of
each transformed coefficient carries data: $N \approx 2b$. The product
has about $4n \cdot b \approx 2nN$ bits, so each transform handles about
twice the size of the product. `mul_mfa_truncate_sqrt2` performs the same
steps with the matrix Fourier algorithm, and `fft_convolution` performs
steps 3 to 5 on given coefficient vectors, also with signed
coefficients, which are then read from their residues in the symmetric
range $(-p/2, p/2)$. When only a window of the coefficients is wanted, as
in a middle product, `fft_convolution_choose_mid` shortens the cyclic
convolution so that the wrap-around misses the window
([](fft/convolution#mid)), as the small-prime FFT does
([](algorithms/fft/multimodular#windows)).

## Multiplication modulo $2^N + 1$ {#mulmod}

The pointwise products in step 4 are products modulo $2^N + 1$. Up to
`FFT_MULMOD_2EXPP1_CUTOFF` limbs (128 limbs, 8192 bits),
`flint_mpn_mulmod_2expp1_basecase` computes the full product with
`flint_mpn_mul_n` and subtracts its high half from its low half. Beyond
that, `fft_mulmod_2expp1` applies the method of Schönhage and Strassen
once more, now with a *negacyclic* convolution, which reduces modulo
$2^N + 1$ for free.

Write $N = 2n' b'$ with $n' = 2^{d'}$ and cut $x$ and $y$, both less than
$2^N$, into $2n'$ pieces $x_i$, $y_i$ of $b'$ bits. Since
$(2^{b'})^{2n'} = 2^N \equiv -1$, the product $xy$ modulo $2^N + 1$ is the
value at $2^{b'}$ of $X(t) Y(t)$ modulo $t^{2n'} + 1$, whose coefficients
are

$$
c_k = \sum_{i + j = k} x_i y_j - \sum_{i + j = k + 2n'} x_i y_j,
\qquad 0 \le k < 2n',
$$

so that $|c_k| < 2n' \, 2^{2b'}$. The negacyclic convolution is a cyclic
convolution of the weighted vectors $(\varphi'^i x_i)$ and
$(\varphi'^i y_i)$, where $\varphi'$ is a root of unity of order $4n'$,
followed by multiplication of coefficient $i$ by $\varphi'^{-i}$
[@Nus1980; @CF1994]. FLINT computes it modulo $p' = 2^{N'} + 1$ with
$N' = 2b' = n' w'$ and $\varphi' = \sqrt 2^{\,w'}$ (`fft_negacyclic`,
`ifft_negacyclic`). This ring is too small for the coefficients by
$d' + 2$ bits; the missing information comes from a second, cheap
convolution.

::: lemma Recovering the coefficients {#lem-crt}
Let $p' = 2^{N'} + 1$ with $\bits \le N'$, and let $B = 2^{\bits}$.
Suppose $|c| < (B/2 - 1) p'$, $v \in [0, p')$ with $v \equiv c \pmod{p'}$
and $r \equiv c \pmod{B}$. Then $c = v + q p'$, where $q$ is the integer
in $[-B/2, B/2)$ congruent to $r - v$ modulo $B$.
:::

::: proof
Since $B$ divides $2^{N'}$, $p' \equiv 1 \pmod{B}$. Write $c = v + q p'$
with an integer $q$. Reducing modulo $B$ gives
$q \equiv c - v \equiv r - v$. From $0 \le v < p'$ and the bound on $|c|$,
$-(B/2 - 1)p' - p' < q p' < (B/2 - 1)p'$, so $-B/2 < q < B/2 - 1$, and
$q$ is determined by its residue modulo $B$.
:::

The coefficients satisfy $|c_k| < 2n' \, 2^{N'} < (B/2 - 1) p'$ as long as
$4n' < B - 2$, which holds for any feasible transform. The residues
modulo $B$ are the negacyclic convolution of the lowest limbs of the
pieces in wrapping word arithmetic, `fft_naive_convolution_1`, which
costs $4n'^2$ word multiplications, negligible next to the rest. Without
this trick the inner ring would need $2b' + d' + 2$ bits, which is not of
the form $n' w'$ with a multiple of $\bits$. Finally,
$xy \bmod (2^N + 1)$ is assembled from the $c_k$ as in step 6 of
[](#alg-mul), with a correction for the negative ones.

::: algorithm Multiplication modulo $2^N + 1$ (`fft_mulmod_2expp1`) {#alg-mulmod}
**Input:** reduced $x$, $y$ modulo $2^N + 1$, with $N/\bits$ a value
returned by `fft_adjust_limbs`.

**Output:** $xy \bmod (2^N + 1)$, reduced.

1. If $x = 2^N$ or $y = 2^N$, return the negation of the other.
2. If $N \le 128 \bits$, return `flint_mpn_mulmod_2expp1_basecase`.
3. $d \gets \lceil \log_2 N \rceil$ and $d' \gets \lfloor d/2 \rfloor - k$
   with $k$ from the table `MULMOD_TAB`; $n' \gets 2^{d'}$,
   $b' \gets N/(2n')$, $w' \gets N/n'^2$.
4. Compute the negacyclic convolution of the $b'$-bit pieces modulo
   $2^{2b'} + 1$ with transforms of length $2n'$, using
   `flint_mpn_mulmod_2expp1_basecase` for its pointwise products, and the
   negacyclic convolution of the lowest limbs modulo $B$.
5. Recover the $c_k$ by [](#lem-crt) and return
   $\sum_k c_k 2^{b' k} \bmod (2^N + 1)$.
:::

`fft_adjust_limbs` makes step 3 well defined: it rounds the limb count up
so that $N$ is a multiple of $n'^2$, making $w'$ an integer, and $b'$ a
multiple of $\bits$. The pointwise products of step 4 are not split
further, so `fft_mulmod_2expp1` has one level of recursion; its inner
products of $2b'$ bits are full products by `flint_mpn_mul_n`. In builds
with `fft_small`, the convolutions of `fmpz_poly` and `fmpzi` use instead
`sd_fft_mpn_mulmod_2expp1`, a negacyclic convolution modulo word-size
primes in floating point, whenever `fft_adjust_limbs_sd_fft` finds a
suitable ring; it is about twice as fast for $N \ge 16384$. See
[](fft/mulmod#sd-fft) and, for the algorithm,
[](fft_small/mulmod_2expp1#algorithm).

## Choice of parameters {#parameters}

`flint_mpn_mul_fft_main` chooses $(\mathit{depth}, w)$ in two steps.

1. It tries $(\mathit{depth}, w)$ in the order
   $(6, 1), (6, 2), (7, 1), (7, 2), \ldots$ and stops at the first pair for
   which chunks of $b = \lfloor (nw - \mathit{depth} - 1)/2 \rfloor$ bits
   fit, $j_1 + j_2 - 1 \le 4n$.
2. If $\mathit{depth} < 11$, it shortens the transform by the factor $2^k$
   given by the table entry `FFT_TAB[depth - 6][w - 1]` and multiplies $w$
   by $4^k$, which keeps the capacity $4n \cdot b$ roughly the same with
   $2^k$ times fewer, larger coefficients. It then lowers $w$ in steps of
   $\max(1, 2^{6 - \mathit{depth}})$, keeping $nw$ a multiple of 64, as
   long as the product still fits, and calls `mul_truncate_sqrt2`.
   Otherwise, if $j_1 + j_2 - 1 \le 3n$, it halves the transform length and
   triples $w$, and then calls `mul_mfa_truncate_sqrt2`.

For small products, few large coefficients win because the pointwise
products, by Toom–Cook and basecase multiplication, are faster per bit
than additional transform layers; for large products, the transform must
stay long so that the coefficients and the pointwise products remain
small. The tables are tuned per processor with `src/fft/tune/tune-fft.c`.
For $\mathit{depth} = 6, \ldots, 10$, each with $w = 1, 2$:

| target | `FFT_TAB` |
|---|---|
| generic, ARM64, Apple M1, Skylake, Zen 3 | $(4, 4), (4, 3), (3, 2), (2, 1), (2, 1)$ |
| Ice Lake, Sapphire Rapids, Meteor Lake, Arrow Lake | $(4, 4), (4, 3), (3, 2), (2, 2), (2, 1)$ |
| Raptor Lake | $(4, 4), (4, 3), (3, 2), (2, 2), (1, 1)$ |
| Tiger Lake | $(4, 4), (4, 3), (3, 2), (2, 1), (1, 1)$ |
| Broadwell | $(4, 4), (4, 3), (3, 3), (3, 2), (2, 1)$ |
| Zen 4 | $(4, 4), (4, 3), (3, 2), (2, 2), (2, 0)$ |

With the generic table, balanced products of $\ell$ limbs each get the
following parameters ("first fit" is the pair found in step 1,
"truncate" and "MFA" stand for `mul_truncate_sqrt2` and
`mul_mfa_truncate_sqrt2`):

| $\ell$ | first fit | used | function | $4n$ | $N$ | $b$ | points | pointwise products |
|---:|---|---|---|---:|---:|---:|---:|---|
| 100 | (6, 2) | (2, 416) | truncate | 16 | 1664 | 830 | 15 | basecase |
| 1000 | (8, 2) | (6, 16) | truncate | 256 | 1024 | 508 | 251 | basecase |
| 10000 | (10, 1) | (8, 10) | truncate | 1024 | 2560 | 1275 | 1003 | basecase |
| 32000 | (10, 2) | (9, 8) | truncate | 2048 | 4096 | 2043 | 2005 | basecase |
| $10^5$ | (11, 2) | (11, 2) | MFA | 8192 | 4096 | 2042 | 6269 | basecase |
| $3 \cdot 10^5$ | (12, 2) | (11, 6) | MFA | 8192 | 12288 | 6138 | 6257 | 16 pieces of 768 bits |
| $10^6$ | (13, 1) | (13, 1) | MFA | 32768 | 8192 | 4089 | 31303 | basecase |
| $10^7$ | (15, 1) | (14, 3) | MFA | 65536 | 49152 | 24568 | 52101 | 32 pieces of 1536 bits |
| $10^8$ | (16, 2) | (15, 6) | MFA | 131072 | 196608 | 98296 | 130219 | 128 pieces of 1536 bits |

The pointwise products of `fft_mulmod_2expp1` use
$d' = \lfloor \lceil \log_2 N \rceil / 2 \rfloor - k$ with $k$ from
`MULMOD_TAB`, indexed by $\lceil \log_2 N \rceil - 12$ and clamped at the
last entry. The values lie between 1 and 4; the generic table has 19
entries, $4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1$. The
pieces then have between about $\sqrt N$ and $8 \sqrt N$ bits, as in the
last column above. In builds without `fft_small`, `flint_mpn_mul`
switches from GMP to `flint_mpn_mul_fft_main` at
`FLINT_FFT_MUL_THRESHOLD` $= 32000$ limbs.

## Complexity {#complexity}

In this section $n$ is the bit size of the operands. A transform of
length $L$ with $N$-bit coefficients costs $O(L \log L)$ butterflies of
$O(N)$ bit operations. With $LN \approx 4n$, the transforms cost $O(n \log n)$ bit operations,
and the $L$ pointwise products cost $L \cdot M(N)$, where $M(N)$ is the
cost of a multiplication modulo $2^N + 1$. Choosing $N \approx \sqrt n$
and applying the method recursively to the pointwise products gives the
bound $O(n \log n \log \log n)$ of [@SS1971], since each of the
$O(\log \log n)$ levels costs $O(n \log n)$. FLINT stops after one
negacyclic level and uses the best available multiplication for the
remaining products of a few thousand bits. In practice the cost per bit
grows slowly with $n$ and is dominated by the memory traffic of the
transforms, which the [matrix Fourier algorithm](algorithms/fft/matrix-fourier)
keeps in cache. FFTs modulo word-size primes in floating point are faster
on current processors ([](algorithms/fft/small-prime-fft));
[](algorithms/fft/choosing) compares the methods.
