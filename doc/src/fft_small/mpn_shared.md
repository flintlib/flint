---
title: Integer products with a common operand
functions: _fft_small_mpn_mul_shared fft_small_mpn_mul_precomp_struct _fft_small_mpn_mul_precomp_size _fft_small_mpn_mul_precomp
---

A product of two integers by `fft_small` costs, per prime, two
conversions and forward transforms, a pointwise product and an inverse
transform, followed by the reconstruction. When one factor $x$ occurs in
several products, its conversion and transform can be computed once.
`_fft_small_mpn_mul_shared` does this for two products $ax$ and $bx$, and
the `_fft_small_mpn_mul_precomp` functions keep the transform of $x$ for
any number of products, possibly on several threads.
`flint_mpn_mul_shared` in `mpn_extras` is the user-level entry point for
the first.

## Contract {#contract}

- All lengths are at least 1, in any order of size. The outputs receive
  $a_n + x_n$, $b_n + x_n$ and $x_n + y_n$ limbs, the top one possibly
  zero.
- `_fft_small_mpn_mul_shared(z1, a, an, z2, b, bn, x, xn)` writes `z1`
  before it reads `b`: `z1` may overlap `a` and `x` but not `b`, and `z2`
  may overlap any input but not `z1`. `flint_mpn_mul_shared` is stricter:
  its outputs must not overlap each other or the inputs.
- `_fft_small_mpn_mul_precomp_init(P, x, xn, yn_max, storage)` prepares
  products by partners of at most `yn_max` limbs; a longer partner is not
  allowed. `storage` is `NULL` or an array of
  `_fft_small_mpn_mul_precomp_size(xn, yn_max)` doubles aligned to 32
  bytes, which must outlive `P`. The integer $x$ itself is not referenced
  after initialisation.
- `_fft_small_mpn_mul_precomp(z, P, y, yn)` allows `z` to overlap `y`.
  `P` is only read, so several threads may use it at once; each product
  uses the scratch buffer of the calling thread's default context
  (`get_default_mpn_ctx`), which invalidates earlier pointers into that
  buffer.
- There are no size cutoffs: for short operands these functions are
  slower than separate products. `flint_mpn_mul_shared` checks
  `_flint_mpn_mul_shared_use_fft` first.

## Algorithm {#algorithm}

Both functions use the geometry that `flint_mpn_mul` would use for the
longest product, of $x$ by the longer partner (`yn_max` for the
precomputed version).

**Two primes** (all operands of at most `FLINT_MPN_MUL_NP2_MAX_AN` =
40000 limbs). The chunk size and depth are chosen as for a product of
these lengths ([](algorithms/fft/multimodular#two-prime)), and $x$ is
packed and transformed modulo both primes up to the truncation of the
longest product. Each product is then computed in one pass while its
data is in cache: the partner is packed, transformed modulo the first
prime, multiplied pointwise and transformed back, then the same modulo
the second prime, and the result is recomposed with one Garner step. The
partner's transforms are truncated to its own product, so a short partner
costs less than a long one.

**Several primes** (longer operands). The plan of a product of $x$ by the
longest partner with at least three primes is used. The transform of $x$
is computed with `fft_small_fft_mpn`; for each partner, the transform is
computed with `_fft_small_fft_mpn_trunc` at the truncation of its own
product, then the pointwise product (with the plan's normalisation) and
the inverse transform run prime by prime, threaded over the primes from
3072 transform points, and the result is reconstructed by
`fft_small_export_mpn`.

The precomputed version keeps only the points of the transform of $x$
that the longest product reads, per prime. Its products save the
conversion and transform of $x$, about 30% of a two-prime product.

Sharing pays as soon as one of the two products is in the FFT range: the
rule of `_flint_mpn_mul_shared_use_fft` is that the longer product passes
the FFT test of `flint_mpn_mul` with the threshold
$T$ = `FLINT_MPN_MUL_SHARED_THRESHOLD` and the shorter one with $T/2$.
$T$ is `FLINT_FFT_SMALL_MUL_THRESHOLD` (500 for the Intel targets from
Skylake on, 240 on Zen 3 and 4) unless `flint-mparam.h` sets it; the
`sapphirerapids` target sets 240.

::: summary algorithms/fft/multimodular :::

## Examples {#examples}

Two products with a common factor of 3000 limbs, compared with
`flint_mpn_mul`, and products by a fixed factor with a precomputed
transform:

```c example name=fft_small_mpn_shared requires=fft_small
#include <flint/flint.h>
#include <flint/mpn_extras.h>
#include <flint/fft_small.h>

int main(void)
{
    slong an = 2500, bn = 800, xn = 3000, i;
    nn_ptr a, b, x, z1, z2, r;
    fft_small_mpn_mul_precomp_t P;
    flint_rand_t state;
    int ok = 1;

    flint_rand_init(state);
    a = flint_malloc(an * sizeof(ulong));
    b = flint_malloc(bn * sizeof(ulong));
    x = flint_malloc(xn * sizeof(ulong));
    z1 = flint_malloc((an + xn) * sizeof(ulong));
    z2 = flint_malloc((bn + xn) * sizeof(ulong));
    r = flint_malloc((an + xn) * sizeof(ulong));
    for (i = 0; i < an; i++) a[i] = n_randlimb(state);
    for (i = 0; i < bn; i++) b[i] = n_randlimb(state);
    for (i = 0; i < xn; i++) x[i] = n_randlimb(state);

    /* z1 = a x and z2 = b x with one transform of x */
    flint_printf("shared transform used by flint_mpn_mul_shared: %d\n",
                 _flint_mpn_mul_shared_use_fft(an, bn, xn));
    _fft_small_mpn_mul_shared(z1, a, an, z2, b, bn, x, xn);

    flint_mpn_mul(r, x, xn, a, an);
    ok = ok && flint_mpn_equal_p(z1, r, an + xn);
    flint_mpn_mul(r, x, xn, b, bn);
    ok = ok && flint_mpn_equal_p(z2, r, bn + xn);
    flint_printf("two products: %s\n", ok ? "agree" : "WRONG");

    /* the transform of x kept for partners of up to an limbs */
    _fft_small_mpn_mul_precomp_init(P, x, xn, an, NULL);
    for (i = 1; i <= 3; i++)
    {
        slong yn = (an * i) / 3;

        _fft_small_mpn_mul_precomp(z1, P, a, yn);
        flint_mpn_mul(r, x, xn, a, yn);
        flint_printf("x times %wd limbs: %s\n", yn,
            flint_mpn_equal_p(z1, r, xn + yn) ? "agrees" : "WRONG");
    }
    _fft_small_mpn_mul_precomp_clear(P);

    flint_free(a);
    flint_free(b);
    flint_free(x);
    flint_free(z1);
    flint_free(z2);
    flint_free(r);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}
```
```output
shared transform used by flint_mpn_mul_shared: 1
two products: agree
x times 833 limbs: agrees
x times 1666 limbs: agrees
x times 2500 limbs: agrees
```
