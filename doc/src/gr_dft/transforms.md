---
title: Applying transforms
functions: gr_dft_precomp gr_dft_inverse_precomp gr_dft gr_dft_inverse
---

`gr_dft_precomp` and `gr_dft_inverse_precomp` apply a plan to a vector
of length $n$; `gr_dft` and `gr_dft_inverse` construct a temporary plan
with the default root for a single transform. All four return a `gr`
status code, which must be checked: a transform that does not return
`GR_SUCCESS` leaves meaningless values in the output.

## The inverse transform and division by n {#division}

The inverse transform is computed as the forward transform with $w^{-1}$,
followed by a division of every entry by $n$ with `_gr_vec_div_scalar_ui`.
The division is the ring's own operation, so its outcome depends on the
ring:

- in a field, or modulo $m$ with $\gcd(n, m) = 1$, it always succeeds;
- modulo $m$ with $\gcd(n, m) > 1$ it fails with `GR_DOMAIN`, and the
  transform cannot be inverted;
- in $\Z$ it succeeds exactly when every entry is divisible by $n$, which
  is the case when the input is the transform of an integer vector;
- in the fixed-point contexts it truncates with an error below 1 ulp.

Because the division is the last step, the status of the inverse
transform also tells whether the input was in the image of the forward
transform, in rings such as $\Z$. For $n = 1$ nothing is divided.

## Status codes {#status}

The transforms combine the status of every ring operation with `|` and
continue to the end even after a failure. Over the rings with exact
arithmetic that are typically used, transforms never fail except for the
division above. Over the fixed-point contexts, an overflow of the
representable range during the transform gives `GR_UNABLE`, so
`GR_SUCCESS` certifies that the input was scaled correctly (see
[](gr_dft/nfixed)). Over `acb`, transforms always succeed.

## Aliasing and memory {#memory}

Input and output may be the same vector for all plans. The radix-2 and
four-step algorithms work in place: they copy the input to the output
first if the two differ. The other algorithms write their output
separately from their input and use a temporary vector of length $n$ when
the two are the same. Besides this, the transforms allocate temporary
storage: a few elements for radix-2 and split-radix, $2p$ elements for
mixed-radix with largest prime factor $p$, $n$ elements for the
transposition of the four-step algorithm, $n + 2 \max(n_1, n_2)$
elements for the prime factor algorithm with $n = n_1 n_2$, and twice the
convolution length (between $2n$ and $4n$) for Bluestein's algorithm.

## One-shot transforms {#one-shot}

`gr_dft` and `gr_dft_inverse` construct a plan with
`gr_dft_precomp_init` and `GR_DFT_ALG_AUTO`, apply it and free it. They
need a ring with default roots, and return `GR_UNABLE` or `GR_DOMAIN`
otherwise; over $\Z/m\Z$ a plan with an explicit root is always needed.
Since the plan costs about as much as a few transforms and does not use
complex mode, repeated transforms over `acb` should use a plan from
`gr_dft_precomp_init_karatsuba`, or the functions of [](gr_dft/acb).

## Convolutions {#convolution}

If $w$ is a principal $n$-th root of unity, the transform turns cyclic
convolution into pointwise multiplication:

$$
\mathrm{DFT}(a \ast b) = \mathrm{DFT}(a) \cdot \mathrm{DFT}(b),
\qquad
(a \ast b)_k = \sum_{i + j \equiv k \pmod n} a_i b_j.
$$

The product of two polynomials with $\deg f + \deg g < n$ is therefore the
inverse transform of the pointwise product of the transforms of their
coefficient vectors, padded with zeros to length $n$. Since the order of
the transformed values does not matter here, a plan with
`GR_DFT_SCRAMBLED` can be used for powers of two. It saves the two
bit-reversal passes; the inverse transform then expects its input in the
same scrambled order, which the pointwise product preserves.

For the multiplication of polynomials over $\Z/p\Z$ with word-size $p$,
the specialised modules `fft_small` and `nmod_poly` are much faster than
generic transforms; see [](algorithms/polynomials/multiplication).

## Example {#example}

The product $(1 + 2x + 3x^2 + 4x^3)(5 + 6x + 7x^2)$ modulo the prime
$p = 998244353 = 119 \cdot 2^{23} + 1$, by a cyclic convolution of
length 8. Since 3 is a primitive root modulo $p$, the root
$3^{(p-1)/8}$ has order 8.

```c example name=gr_dft_convolution
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_vec.h>
#include <flint/gr_dft.h>

/* (1 + 2x + 3x^2 + 4x^3)(5 + 6x + 7x^2) modulo p = 119 * 2^23 + 1,
   by a cyclic convolution of length 8 */
int main(void)
{
    gr_ctx_t ctx;
    gr_dft_pre_t P;
    gr_ptr a, b, w;
    slong j, n = 8;
    ulong p = 998244353;
    slong sz;

    gr_ctx_init_nmod(ctx, p);
    sz = ctx->sizeof_elem;
    a = gr_heap_init_vec(n, ctx);
    b = gr_heap_init_vec(n, ctx);
    w = gr_heap_init(ctx);

    /* 3 is a primitive root, so 3^((p - 1)/8) has order 8 */
    GR_MUST_SUCCEED(gr_set_ui(w, 3, ctx));
    GR_MUST_SUCCEED(gr_pow_ui(w, w, (p - 1) / n, ctx));
    GR_MUST_SUCCEED(gr_dft_precomp_init_root(P, w, n, GR_DFT_ALG_AUTO,
                                             GR_DFT_SCRAMBLED, ctx));

    /* the remaining entries are zero after gr_heap_init_vec */
    for (j = 0; j < 4; j++)
        GR_MUST_SUCCEED(gr_set_si(GR_ENTRY(a, j, sz), j + 1, ctx));
    for (j = 0; j < 3; j++)
        GR_MUST_SUCCEED(gr_set_si(GR_ENTRY(b, j, sz), j + 5, ctx));

    GR_MUST_SUCCEED(gr_dft_precomp(a, a, P, ctx));
    GR_MUST_SUCCEED(gr_dft_precomp(b, b, P, ctx));
    GR_MUST_SUCCEED(_gr_vec_mul(a, a, b, n, ctx));
    GR_MUST_SUCCEED(gr_dft_inverse_precomp(a, a, P, ctx));

    _gr_vec_print(a, n, ctx);
    flint_printf("\n");

    gr_dft_precomp_clear(P);
    gr_heap_clear(w, ctx);
    gr_heap_clear_vec(a, n, ctx);
    gr_heap_clear_vec(b, n, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```
```output
[5, 16, 34, 52, 45, 28, 0, 0]
```
