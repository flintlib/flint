---
title: Multiplying polynomial matrices modulo n
functions: nmod_poly_mat_mul nmod_poly_mat_mul_classical nmod_poly_mat_mul_interpolate nmod_poly_mat_mulmid_fft_small
---

`nmod_poly_mat_mul` multiplies matrices whose entries are polynomials
over $\Z/n\Z$. Depending on the dimensions, the lengths of the entries
and the modulus, it computes the product entry by entry with polynomial
multiplications, as one large integer matrix product by Kronecker
substitution, as many small matrix products by evaluation and
interpolation, or with the small-prime FFT, transforming each entry only
once. All methods give the same result.

::: summary algorithms/linear-algebra/polynomial-matrices :::

## Choice of method {#dispatch}

Let $d$ be the smallest of the three dimensions, $\ell_A$ and $\ell_B$
the largest lengths of the entries of $A$ and $B$, and
$\ell = \min(\ell_A, \ell_B)$.

| Condition | Method |
|---|---|
| $d \ge 2$, $\ell \ge 2d$, $\ell \ge c$ | `nmod_poly_mat_mulmid_fft_small`, if it succeeds |
| $d < 10$ | `nmod_poly_mat_mul_classical` |
| $n$ prime, $\log_2 n > \beta/4$, $n \ge \ell_A + \ell_B - 1$, $d > 60 + \sqrt{\ell}$ | `nmod_poly_mat_mul_interpolate` |
| $\max(\ell_A, \ell_B) > 128$ | `nmod_poly_mat_mul_classical` |
| otherwise | `nmod_poly_mat_mul_KS` |

Here $c = 64$, or 32 if $n$ has more than 40 bits, and $c$ is doubled
for $d = 2$. The classical method pays off for small dimensions, where
there is little to share between the products of entries; Kronecker
substitution for short entries, where it produces integer matrices with
entries of a few hundred bits that `fmpz_mat_mul` multiplies efficiently
([](algorithms/linear-algebra/integer-matrix-multiplication)); evaluation
and interpolation for large matrices over a large prime field, where the
fast `nmod_mat_mul` kernels do the work; and the transformed method for
long entries.

## Window of coefficients {#mulmid}

`nmod_poly_mat_mulmid_fft_small(C, A, B, zl, zh)` computes only the
coefficients of degrees $z_l, \ldots, z_h - 1$ of each entry of $AB$, as
the entries of $C$ (shifted down by $z_l$). This is the middle product
used, for example, in Newton iterations for matrix power series, where
the low part of a product is known in advance. Its cost is that of
transforming the $mk + kn$ input entries and the $mn$ outputs once, plus
$mkn$ pointwise products, with transforms of the length needed for the
window; see [](algorithms/fft/multimodular#windows). It returns 0 when
`fft_small` is not available or no plan fits the accumulation bound; the
caller must then use another method.

## Example {#example}

Products of $4 \times 4$ matrices with entries of length 150 modulo a
62-bit prime, by all methods, and a window of coefficients:

```c example name=nmod_poly_mat_mul_methods
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

int main(void)
{
    nmod_poly_mat_t A, B, C, D;
    flint_rand_t state;
    ulong p = n_nextprime(UWORD(1) << 61, 1);
    slong i, j, d = 4, len = 150;
    int ok, all = 1;

    flint_rand_init(state);
    nmod_poly_mat_init(A, d, d, p);
    nmod_poly_mat_init(B, d, d, p);
    nmod_poly_mat_init(C, d, d, p);
    nmod_poly_mat_init(D, d, d, p);
    nmod_poly_mat_rand(A, state, len);
    nmod_poly_mat_rand(B, state, len);

    nmod_poly_mat_mul(C, A, B);
    flint_printf("degree of product: %wd\n", nmod_poly_mat_degree(C));

    nmod_poly_mat_mul_classical(D, A, B);
    flint_printf("classical:   %d\n", nmod_poly_mat_equal(C, D));
    nmod_poly_mat_mul_KS(D, A, B);
    flint_printf("KS:          %d\n", nmod_poly_mat_equal(C, D));
    nmod_poly_mat_mul_interpolate(D, A, B);
    flint_printf("interpolate: %d\n", nmod_poly_mat_equal(C, D));

    /* coefficients 100 to 199 of each entry */
    ok = nmod_poly_mat_mulmid_fft_small(D, A, B, 100, 200);
    if (ok)
    {
        for (i = 0; i < d; i++)
            for (j = 0; j < d; j++)
            {
                nmod_poly_t t;
                nmod_poly_init(t, p);
                nmod_poly_shift_right(t, nmod_poly_mat_entry(C, i, j), 100);
                nmod_poly_truncate(t, 100);
                all = all && nmod_poly_equal(t, nmod_poly_mat_entry(D, i, j));
                nmod_poly_clear(t);
            }
    }
    flint_printf("mulmid:      %d\n", all);

    nmod_poly_mat_clear(A);
    nmod_poly_mat_clear(B);
    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(D);
    flint_rand_clear(state);
    return 0;
}
```
```output
degree of product: 298
classical:   1
KS:          1
interpolate: 1
mulmid:      1
```
