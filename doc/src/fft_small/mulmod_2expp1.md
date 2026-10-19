---
title: Multiplication modulo 2^N + 1
functions: sd_fft_mpn_mulmod_2expp1_ctx_struct sd_fft_mpn_mulmod_2expp1_ctx_init sd_fft_mpn_mulmod_2expp1_scratch_init sd_fft_mpn_mulmod_2expp1 sd_fft_mpn_mulmod_2expp1_transformed_size sd_fft_mpn_mulmod_2expp1_choose_np
---

Schönhage–Strassen multiplication reduces a large product to many
products modulo $2^N + 1$ for a moderately large $N$
([](algorithms/fft/schonhage-strassen)). The `fft` module computes these
pointwise products with the functions here when the ring has at least
16384 bits (`FLINT_FFT_SD_FFT_MIN_RING`), choosing $N$ with
`fft_adjust_limbs_sd_fft` so that a good digit splitting exists. They
compute a product modulo $2^N + 1$ as a negacyclic convolution of $b$-bit
digits with the small-prime transforms.

## Usage {#usage}

1. Create a context for $N$ with `sd_fft_mpn_mulmod_2expp1_ctx_init`,
   passing `get_default_mpn_ctx()` and $m = 0$ to let it choose the digit
   count. $N$ must be a positive multiple of `FLINT_BITS`; otherwise an
   exception is raised.
2. Create one scratch object per thread with
   `sd_fft_mpn_mulmod_2expp1_scratch_init`. The context is read-only
   during multiplications and can be shared.
3. Multiply with `sd_fft_mpn_mulmod_2expp1`. Operands and result have
   $n + 1$ limbs, $n = N / \bits$; operands may have top limb 0 or 1 with
   arbitrary low limbs, and the result is fully reduced into $[0, 2^N]$.
   The result may alias an operand.
4. For many products with the same $y$, store its transform once with
   `sd_fft_mpn_mulmod_2expp1_transform` in an array of
   `sd_fft_mpn_mulmod_2expp1_transformed_size(C)` doubles, aligned to 32
   bytes, and multiply with `sd_fft_mpn_mulmod_2expp1_mul_cached`.

## Algorithm {#algorithm}

Write $x = \sum_{i<m} x_i 2^{bi}$ with digits $0 \le x_i < 2^b$, where
$N = mb$ and a top limb 1 contributes $2^N \equiv -1$, which is
subtracted from $x_0$. Since $2^{bm} \equiv -1$, the product is
$\sum_k z_k 2^{bk}$ with the negacyclic convolution

$$
z_k = \sum_{i + j = k} x_i y_j - \sum_{i + j = k + m} x_i y_j,
\qquad |z_k| < m 2^{2b} = 2^{2b + d},
$$

where $m = 2^d$. Modulo each prime $p$, with $\theta$ a primitive
$2m$-th root of unity ($\theta^m = -1$), the weighted sequences
$\theta^i x_i$ and $\theta^j y_j$ have the cyclic convolution
$\theta^k z_k$ [@CF1994], which a transform of length $m$ computes. The
residues are untwisted by $\theta^{-k}$, the factor $m^{-1}$ of the
inverse transform and the Chinese remainder normalisation being folded
into the same table. To use the unsigned reconstruction, $2^{2b + d}$ is
added to every coefficient, which makes them lie in $[0, 2^{2b + d + 1})$,
so the product of the primes must be at least $2^{2b + d + 1}$. After
reconstructing $\sum_k (z_k + 2^{2b+d}) 2^{bk}$ the bias polynomial is
subtracted, and the part of the signed result above $2^N$ is folded back
by subtraction.

The number of primes is the least $k \ge 2$ with enough capacity for
which digits of $b$ bits can be converted: digits of at most 50 bits are
converted directly, larger ones need one of the vectorised conversions of
the integer multiplication profiles. With $m = 0$ the context compares
the admissible digit counts $m \ge 16$ (with $20 \le b \le 200$) by an
estimate $w_k \, k \, m (d + 4 + b/16)$ with measured weights $w_k$ per
prime count, which favours two primes slightly and penalises eight.

| $N$ | $m$ | $b$ | primes |
|---|---|---|---|
| $2^{10}$ to $2^{24}$ | $N / 64$ | 64 | 3 |
| 10752 | 256 | 42 | 2 |
| 30720 | 256 | 120 | 6 |
| 47104 | 512 | 92 | 4 |
| 344064 | 8192 | 42 | 2 |

::: summary algorithms/fft/multimodular :::

## Example {#example}

```c example name=fft_small_mulmod_2expp1 requires=fft_small
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fft_small.h>

int main(void)
{
    slong N = 65536, n = N / FLINT_BITS, i;
    sd_fft_mpn_mulmod_2expp1_ctx_t C;
    sd_fft_mpn_mulmod_2expp1_scratch_t S;
    nn_ptr x = flint_calloc(n + 1, sizeof(ulong));
    nn_ptr y = flint_calloc(n + 1, sizeof(ulong));
    nn_ptr z = flint_calloc(n + 1, sizeof(ulong));
    fmpz_t X, Y, Z, M;

    sd_fft_mpn_mulmod_2expp1_ctx_init(C, get_default_mpn_ctx(), N, 0);
    sd_fft_mpn_mulmod_2expp1_scratch_init(S, C);
    flint_printf("N = %wd: %wd digits of %wd bits, %wd primes\n",
                 C->N, C->m, C->b, C->np);

    for (i = 0; i < n; i++)
    {
        x[i] = UWORD(0x0123456789abcdef) * (i + 1);
        y[i] = UWORD(0xfedcba9876543210) ^ i;
    }
    y[n] = 0;
    x[n] = 1;       /* x = 2^N + (low limbs), not reduced */

    sd_fft_mpn_mulmod_2expp1(C, z, x, y, S);

    fmpz_init(X); fmpz_init(Y); fmpz_init(Z); fmpz_init(M);
    fmpz_set_ui_array(X, x, n + 1);
    fmpz_set_ui_array(Y, y, n + 1);
    fmpz_one(M);
    fmpz_mul_2exp(M, M, N);
    fmpz_add_ui(M, M, 1);
    fmpz_mul(Z, X, Y);
    fmpz_mod(Z, Z, M);
    fmpz_set_ui_array(X, z, n + 1);
    flint_printf("x*y mod 2^N + 1: %s\n", fmpz_equal(X, Z) ? "agrees" : "WRONG");

    fmpz_clear(X); fmpz_clear(Y); fmpz_clear(Z); fmpz_clear(M);
    sd_fft_mpn_mulmod_2expp1_scratch_clear(S);
    sd_fft_mpn_mulmod_2expp1_ctx_clear(C);
    flint_free(x);
    flint_free(y);
    flint_free(z);
    flint_cleanup_master();
    return 0;
}
```
```output
N = 65536: 1024 digits of 64 bits, 3 primes
x*y mod 2^N + 1: agrees
```
