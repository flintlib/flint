---
title: Multiplying large integers
functions: mpn_ctx_mpn_mul mpn_mul_default_mpn_ctx _mpn_ctx_mpn_mul_range
---

These functions multiply nonnegative integers given as limb arrays, like
GMP's `mpn_mul`. `flint_mpn_mul` calls `mpn_mul_default_mpn_ctx` from a
few hundred limbs on, so a direct call is only needed to bypass the
dispatch or to use a context of one's own; `_mpn_ctx_mpn_mul_range`
computes only a window of limbs, as `flint_mpn_mulhigh_n`,
`flint_mpn_mullow_n` and `flint_mpn_mulmid` do above their cutoffs.

## Contract {#contract}

- $a_n, b_n \ge 1$; the operands may come in either order, and the output
  must not overlap them. If `a == b` and $a_n = b_n$, the square is
  computed, saving one forward transform per prime.
- `mpn_ctx_mpn_mul` writes all $a_n + b_n$ limbs, the top one possibly
  zero.
- `_mpn_ctx_mpn_mul_range` writes $\mathit{hi} - \mathit{lo}$ limbs, or
  nothing if $\mathit{lo} \ge \mathit{hi}$. Limbs from $a_n + b_n$ on are
  zero. For $\mathit{lo} > 0$ it returns the exact window minus a small
  $\delta \ge 0$, modulo $B^{\mathit{hi} - \mathit{lo}}$ (so a deficit
  can borrow through limbs that are zero in the exact result), with
  $\delta < \max(2, \min(a_n, b_n, \mathit{lo}))$. In practice $\delta$
  is 0, 1 or 2. This is the usual contract of FLINT's high and middle
  products, whose callers allow for a few units of error in the lowest
  limb.
- The context is used for scratch memory and must not be in use on
  another thread. The default context is local to the calling thread.

## Algorithm {#algorithm}

The integers are cut into chunks of $\beta$ bits, $a = \sum_i a_i
2^{\beta i}$ and $b = \sum_j b_j 2^{\beta j}$, and the chunk sequences
are multiplied as polynomials. A coefficient of the product is a sum of
at most $\lceil 64 b_n / \beta \rceil$ products of chunks, so it is less
than $\lceil 64 b_n / \beta \rceil 2^{2 \beta}$; if the product $P$ of the
primes exceeds this, the coefficients are recovered exactly from their
residues. The product is $\sum_k c_k 2^{\beta k}$, assembled with carries.

| Operands | Primes | Chunk size $\beta$ | Threads |
|---|---|---|---|
| $\max(a_n, b_n) \le 40000$ | 2 | 40 to 45 bits | 1 |
| larger | 4 to 8 | 84 to 192 bits | up to 8 |

Up to 40000 limbs (`FLINT_MPN_MUL_NP2_MAX_AN`) the *two-prime pipeline*
is used: the chunks are smaller than both primes, so the conversion into
the transform is bit slicing with no reductions, and the reconstruction
is a single vectorised Garner step. It chooses the shortest transform
and then the largest chunk size in $[40, 45]$ that satisfies the bound.
For larger operands a *profile* $(k, \beta)$ is taken from a fixed table
of fifteen, with $k$ primes and chunks of $\beta \in [84, 192]$ bits,
reduced modulo all $k$ primes at once by vectorised code. Among the
profiles whose bound allows $b_n$, the one with the largest chunks for
each $k$ competes, and the one with the smallest estimated cost $(1 -
r/4)\, k\, d\, t$ wins, where $2^d$ is the transform length, $t$ the
truncated length and $r = t / 2^d$; with several threads only profiles
whose $k$ is a multiple of the thread count qualify. The resulting
choices for balanced products are, for example:

| $a_n = b_n$ (limbs) | 1000 | 40000 | 100000 | $10^6$ | $10^7$ |
|---|---|---|---|---|---|
| primes, chunk bits (one thread) | 2, 44 | 2, 41 | 5, 116 | 6, 136 | 7, 160 |
| transform length | $2^{12}$ | $2^{17}$ | $2^{17}$ | $2^{20}$ | $2^{23}$ |

The transforms are truncated to the product length, rounded up to a
multiple of 256. For a window, only the chunk coefficients that reach the
limbs $[\mathit{lo}, \mathit{hi})$ are reconstructed. Above 40000 limbs,
if a power of two $2^d$ is at least both chunk counts and the highest
needed coefficient, while the coefficients beyond $2^d$ fold onto
coefficients that lie entirely below limb $\mathit{lo}$, the transform is
cyclic of length $2^d$, which can halve it for a high product.

A window of $m$ chunk coefficients that is shorter than the longer
operand (in chunks) is computed, on both paths, as a *transposed product*
when its transform is strictly shorter: the coefficients
$c_{c_0 + j} = \sum_i r_i v_{i+j}$, $0 \le j < m$, with $r$ the chunks
of the shorter operand reversed and $v$ the chunks of the longer one
from index $c_0 - (r_n - 1)$ on. This is the transpose of a product of
lengths $r_n$ and $m$, computed with the transposed transforms
(`sd_ifft_trunc_tr` on $v$, the forward transform of $r$, a pointwise
product and `sd_fft_trunc_tr`), so its transform has length about
$r_n + m$ instead of the length of the full product: for the middle
product of $2n$ by $n$ limbs that a Newton iteration needs, about $2n$
instead of $3n$ limbs' worth of chunks. $c_0$ and $r_n$ are rounded so
that $v$ starts on a limb boundary and the reconstruction sees the same
coefficients at the same bit positions, so the result is the same, bit
for bit, as that of the window of the full product. Squarings do not use
it. The details, including the proofs of the bounds, are on the
algorithm pages; for the transposition principle see
[](algorithms/polynomials/middle-product#transposed-fft).

::: summary algorithms/fft/multimodular :::

::: summary algorithms/fft/small-prime-fft :::

## Memory and threads {#resources}

The transforms of both operands are kept for all primes: $2 k$ arrays of
$2^d$ doubles, taken from the context's scratch buffer, which is kept for
the next call. For two operands of 40000 limbs this is 4 MiB, about six
times the size of the product. Above 40000 limbs the work is spread over up to
eight threads of FLINT's thread pool, by prime for the transforms and by
ranges of coefficients for the reconstruction.

## Example {#example}

Multiplying $3^{300000}$ by $7^{100000}$, and computing the top four limbs
of the product alone:

```c example name=fft_small_mpn_mul requires=fft_small
#include <gmp.h>
#include <flint/flint.h>
#include <flint/mpn_extras.h>
#include <flint/fft_small.h>

int main(void)
{
    mpz_t a, b, c;
    nn_ptr z, w;
    slong an, bn, cn, lo, hi;

    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_ui_pow_ui(a, 3, 300000);
    mpz_ui_pow_ui(b, 7, 100000);
    an = mpz_size(a);
    bn = mpz_size(b);

    /* full product with the default context of this thread */
    z = flint_malloc((an + bn) * sizeof(ulong));
    mpn_mul_default_mpn_ctx(z, mpz_limbs_read(a), an, mpz_limbs_read(b), bn);

    mpz_mul(c, a, b);
    cn = mpz_size(c);
    flint_printf("%wd x %wd limbs: %s\n", an, bn,
        (mpn_cmp(z, mpz_limbs_read(c), cn) == 0
            && flint_mpn_zero_p(z + cn, an + bn - cn)) ? "equal to mpz_mul" : "WRONG");

    /* the top four limbs alone */
    lo = an + bn - 4;
    hi = an + bn;
    w = flint_malloc((hi - lo) * sizeof(ulong));
    _mpn_ctx_mpn_mul_range(get_default_mpn_ctx(), w, lo, hi,
                           mpz_limbs_read(a), an, mpz_limbs_read(b), bn);
    flint_printf("limbs [%wd, %wd): %s\n", lo, hi,
        mpn_cmp(w, z + lo, hi - lo) == 0 ? "exact" : "low by a few units");

    flint_free(z);
    flint_free(w);
    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    flint_cleanup_master();
    return 0;
}
```
```output
7430 x 4387 limbs: equal to mpz_mul
limbs [11813, 11817): exact
```
