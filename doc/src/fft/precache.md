---
title: Convolutions with a fixed operand
functions: fft_precache fft_convolution_precache fft_convolution_precache_sd_fft
---

When many polynomials are multiplied by the same polynomial $g$, its
transform needs to be computed only once. `fft_precache` replaces the
coefficients of $g$ by their transform, and `fft_convolution_precache`
then multiplies any number of polynomials by $g$ at the cost of one
forward transform, the pointwise products and one inverse transform each,
instead of two forward transforms, which saves about a third of the
work. The public functions `fmpz_poly_mul_SS_precache_init` and
`fmpz_poly_mullow_SS_precache` are built on these functions.

## Usage {#usage}

1. Choose `depth` and `limbs` as for `fft_convolution` (see
   [Parameters and layout](fft/convolution#layout)), for the longest
   product that will be computed, and set `trunc` to its length.
2. Store $g$ in an array `jj` of $4n$ coefficients, zero-padded, and
   call `fft_precache(jj, depth, limbs, trunc, t1, t2, s1)`. Afterwards
   `jj` holds the transform with every used entry reduced.
3. For each multiplicand, store it in a zero-padded array `ii` and call
   `fft_convolution_precache(ii, jj, depth, limbs, trunc2, t1, t2, s1, tt)`
   with the same `depth` and `limbs` and any $\mathit{trunc2} \le \mathit{trunc}$
   with $\mathit{trunc2} > 2n$ that is at least the length of this product.
   This works because the truncated transform at fewer points is a prefix
   of the one at more points. The product is in `ii`, reduced; `jj` is
   unchanged.

Both functions run in a single thread and use only the first entry of
each scratch array, so pointers to single scratch blocks can be passed
as in the example below. `fft_convolution_precache_sd_fft` with a
negacyclic context (see [](fft/convolution#algorithm)) computes the
pointwise products in threads, without modifying `jj`.

## Example {#example}

Three products with the same polynomial $g$ of length 200 with 20-bit
coefficients. The products have coefficients below $200 \cdot 2^{40}$,
and $N = 128$ is a multiple of $n = 128$.

```c example name=fft_precache
#include <flint/flint.h>
#include <flint/fft.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/mpn_extras.h>

int main(void)
{
    slong len = 200, len_out = 2 * len - 1;
    slong depth = FLINT_CLOG2(len_out) - 2, n = WORD(1) << depth;   /* n = 128 */
    slong limbs = 2, size = limbs + 1;   /* N = 128 bits > 200 * 2^40 */
    slong i, k;
    mp_limb_t ** ii, ** jj, * mem, * t1, * t2, * s1, * tt;
    fmpz_poly_t f, g, h;
    int ok = 1;

    ii = flint_malloc(8 * n * sizeof(mp_limb_t *));
    jj = ii + 4 * n;
    mem = flint_calloc((8 * n + 5) * size, sizeof(mp_limb_t));
    for (i = 0; i < 8 * n; i++)
        ii[i] = mem + i * size;
    t1 = mem + 8 * n * size; t2 = t1 + size; s1 = t2 + size; tt = s1 + size;

    fmpz_poly_init(f); fmpz_poly_init(g); fmpz_poly_init(h);

    /* the fixed operand g, transformed once */
    for (i = 0; i < len; i++)
    {
        fmpz_poly_set_coeff_ui(g, i, (i * 40503 + 7) % 1000003);
        jj[i][0] = (i * 40503 + 7) % 1000003;
    }
    fft_precache(jj, depth, limbs, len_out, &t1, &t2, &s1);

    for (k = 1; k <= 3; k++)
    {
        for (i = 0; i < 4 * n; i++)
            flint_mpn_zero(ii[i], size);
        for (i = 0; i < len; i++)
        {
            fmpz_poly_set_coeff_ui(f, i, (k * i * i + 1) % 999983);
            ii[i][0] = (k * i * i + 1) % 999983;
        }

        fft_convolution_precache(ii, jj, depth, limbs, len_out, &t1, &t2, &s1, &tt);

        fmpz_poly_mul(h, f, g);
        for (i = 0; i < len_out; i++)
            ok = ok && fmpz_equal_ui(h->coeffs + i, ii[i][0]) && ii[i][1] == 0 && ii[i][2] == 0;
        flint_printf("product %wd: coefficient 100 = %wu\n", k, ii[100][0]);
    }
    flint_printf("all equal: %d\n", ok);

    fmpz_poly_clear(f); fmpz_poly_clear(g); fmpz_poly_clear(h);
    flint_free(mem); flint_free(ii);
    return 0;
}
```
```output
product 1: coefficient 100 = 145641580639
product 2: coefficient 100 = 291233620886
product 3: coefficient 100 = 436825661133
all equal: 1
```
