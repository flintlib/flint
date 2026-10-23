---
title: Multiplying limb arrays
functions: flint_mpn_mul flint_mpn_mul_n flint_mpn_sqr flint_mpn_mul_shared
---

`flint_mpn_mul`, `flint_mpn_mul_n` and `flint_mpn_sqr` are FLINT's
replacements for GMP's `mpn_mul`, `mpn_mul_n` and `mpn_sqr`, with the same
interface except for the order requirement of `flint_mpn_mul`
($\mathit{xn} \ge \mathit{yn}$) and the returned top limb. They are faster
than GMP's functions at both ends of the size range: for short operands
because of hand-written assembly routines for each size, and for long
operands because of the small-prime FFT of [`fft_small`](fft_small).
`flint_mpn_mul_shared` computes two products with a common factor.

::: summary algorithms/integers/multiplication :::

## Algorithm selection {#selection}

The functions are inline. On x86-64 with the ADX and BMI2 instructions
(Broadwell and later, Zen and later) and on ARM64, products of short
operands are computed by calling one routine from a table of
[hard-coded assembly kernels](algorithms/integers/basecase-multiplication),
with no further dispatch:

| Architecture | `flint_mpn_mul` | `flint_mpn_mul_n` | `flint_mpn_sqr` |
|---|---|---|---|
| x86-64 with ADX | $\mathit{xn} \le 16$ (any $\mathit{yn} \le \mathit{xn}$) | $n \le 16$ | $n \le 14$ |
| ARM64 | $\mathit{yn} \le 15$ (any $\mathit{xn}$) | $n \le 15$ | $n \le 9$ |
| other | $\mathit{xn} \le 7$, or $\mathit{xn} \le 14$ and $\mathit{yn} = 1$ | $n \le 7$ | none |

On machines without FLINT assembly the table holds portable C code. The
ARM64 kernels are indexed by the shorter length and loop over the longer
operand, so they serve any $\mathit{xn}$. Beyond the table, `_flint_mpn_mul` uses,
in order:

1. `mpn_mul_1` if $\mathit{yn} = 1$ (on ARM64 also `flint_mpn_mul_2` if
   $\mathit{yn} = 2$); on machines without FLINT assembly, GMP's
   `mpn_mul_basecase` for $\mathit{xn} \le 16$;
2. on x86-64 with ADX and on ARM64, `flint_mpn_mul_toom22` (Karatsuba
   recursing into the assembly kernels) if the operands are nearly
   balanced, $5\,\mathit{yn} \ge 4\,\mathit{xn}$, and $\mathit{yn}$ is at
   most 230 (x86-64) or 483 (ARM64);
3. GMP's `mpn_mul` (Toom–Cook) below the FFT range;
4. the FFT: `mpn_mul_default_mpn_ctx` from `fft_small` if available, and
   otherwise the Schönhage–Strassen FFT of the [`fft`](fft) module
   (`flint_mpn_mul_fft_main`, from `FLINT_FFT_MUL_THRESHOLD` limbs).

With `fft_small` and the threshold $T$ = `FLINT_FFT_SMALL_MUL_THRESHOLD`
(1000 limbs in the generic `flint-mparam.h`, between 240 and 1000 limbs
in the machine-specific ones), the FFT is
used when $\mathit{yn} \ge T$, and for unbalanced operands already when
$\mathit{xn} + \mathit{yn} > 2T$ and $\mathit{yn} > T/3$: the cost of
GMP's multiplication grows like $(\mathit{xn}/\mathit{yn}) M(\mathit{yn})$
while that of the FFT depends on $\mathit{xn} + \mathit{yn}$.
`_flint_mpn_mul_n` follows the same steps, and `_flint_mpn_sqr` uses GMP's
`mpn_sqr` below `FLINT_FFT_SMALL_SQR_THRESHOLD` and the FFT above.

## Shared operands {#shared}

`flint_mpn_mul_shared` computes $ax$ and $bx$. With `fft_small`, when
`_flint_mpn_mul_shared_use_fft` says so, it calls
`_fft_small_mpn_mul_shared`, which transforms $x$ once and uses the
transform for both products; otherwise it makes two calls to
`flint_mpn_mul`. With $T$ = `FLINT_MPN_MUL_SHARED_THRESHOLD` (by default
the `fft_small` multiplication threshold), the shared transform is used
when the longer product (by $\max(\mathit{an}, \mathit{bn})$) satisfies
the FFT rule of `flint_mpn_mul` above with threshold $T$ and the shorter
product satisfies it with threshold $T/2$: once the longer product is
transformed, the shorter one costs about its share of the transforms and
need not reach the FFT range on its own. This is used for
instance by the Chinese remaindering trees of
[](algorithms/integers/multimodular-reduction), where each node is
multiplied by two cofactors.

## Example {#example}

The product of two 2000-limb numbers (in the `fft_small` range) agrees
with GMP's, and a shared product agrees with two separate ones:

```c example name=mpn_mul_basic requires=fft_small
#include <gmp.h>
#include <flint/flint.h>
#include <flint/mpn_extras.h>

#define N 2000

int main(void)
{
    flint_rand_t state;
    mp_limb_t x[N], y[N], a[N / 2], z[2 * N], w[2 * N], z1[2 * N], z2[2 * N];

    flint_rand_init(state);
    flint_mpn_urandomb(x, state, N * FLINT_BITS);
    flint_mpn_urandomb(y, state, N * FLINT_BITS);
    flint_mpn_urandomb(a, state, N / 2 * FLINT_BITS);

    flint_mpn_mul_n(z, x, y, N);
    mpn_mul_n(w, x, y, N);
    flint_printf("mul_n agrees with GMP: %d\n", mpn_cmp(z, w, 2 * N) == 0);

    flint_mpn_sqr(z, x, N);
    mpn_sqr(w, x, N);
    flint_printf("sqr agrees with GMP: %d\n", mpn_cmp(z, w, 2 * N) == 0);

    /* z1 = y x and z2 = a x with one transform of x */
    flint_mpn_mul_shared(z1, y, N, z2, a, N / 2, x, N);
    flint_printf("shared transform: %d\n", _flint_mpn_mul_shared_use_fft(N, N / 2, N));
    flint_mpn_mul_n(z, y, x, N);
    flint_mpn_mul(w, x, N, a, N / 2);
    flint_printf("shared products agree: %d\n",
        mpn_cmp(z1, z, 2 * N) == 0 && mpn_cmp(z2, w, 3 * N / 2) == 0);

    flint_rand_clear(state);
    return 0;
}
```
```output
mul_n agrees with GMP: 1
sqr agrees with GMP: 1
shared transform: 1
shared products agree: 1
```
