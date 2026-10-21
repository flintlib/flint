---
title: acb_dft.h
header: acb_dft.h
description: Discrete Fourier transforms of complex ball vectors, as thin wrappers around gr_dft.
category: Special functions
order: 7
---

This module computes the discrete Fourier transform (DFT) of a vector of
complex balls,

$$
w_k = \sum_{j=0}^{n-1} v_j \, e^{-2\pi i jk/n}, \qquad 0 \le k < n,
$$

for any length $n \ge 1$, and its inverse
$v_j = \frac{1}{n} \sum_{k} w_k e^{2\pi i jk/n}$. The output balls
contain the exact transform of every complex vector in the input balls.

All functions are thin wrappers around the complex ball transforms of
[`gr_dft`](gr_dft): `acb_dft` calls `gr_dft_acb`, and
`acb_dft_pre_t` is the plan type `gr_dft_acb_pre_t`. These run the
transform on the midpoints in fixed-point arithmetic with an a priori
error bound, and fall back to ball arithmetic when fixed point does not
apply (non-finite or extremely large midpoints); see [](gr_dft/acb) for
the method and [](algorithms/fft/fixed-point-dft) for the error
analysis. The plan chooses among radix-2, split-radix, four-step,
Good–Thomas, mixed-radix and Bluestein algorithms according to the
factorisation of $n$ ([](algorithms/fft/dft-algorithms)), so that a
transform costs $O(n \log n)$ operations for every $n$.

Earlier versions of this module also provided product DFTs, cyclic
convolutions and direct access to the naive, CRT, cyclic, radix-2 and
Bluestein algorithms. These have been removed: product DFTs of complex
ball vectors are `gr_dft_acb_prod`, transforms on groups of Dirichlet
characters are `acb_dirichlet_dft`, and the individual algorithms,
transforms over other rings and convolutions are in [`gr_dft`](gr_dft).

::: summary algorithms/fft/dft-algorithms :::

## Conventions {#conventions}

- Input and output vectors have length $n$ and may be the same vector.
- With working precision `prec`, the outputs carry about `prec` accurate
  bits *relative to the largest input*: an output much smaller than the
  inputs has fewer accurate bits.
- A length $n \le 0$ is allowed and does nothing.

## Transforms {#transforms}

::: function
void acb_dft(acb_ptr w, acb_srcptr v, slong len, slong prec)
void acb_dft_inverse(acb_ptr w, acb_srcptr v, slong len, slong prec)

Set `w` to the DFT, respectively the inverse DFT (including the division
by $n$), of the vector `v` of length $n$ = `len`. Each call builds a
plan, so repeated transforms of one length are faster with the
precomputed versions below. These call `gr_dft_acb` and
`gr_dft_acb_inverse`.
:::

## Precomputed plans {#precomp}

::: type
acb_dft_pre_struct
acb_dft_pre_t

A precomputed transform of a fixed length, the same type as
`gr_dft_acb_pre_t`. The length is the field `n`. An `acb_dft_pre_t` is
an array of length 1, so that it is passed by reference.
:::

::: function
void acb_dft_precomp_init(acb_dft_pre_t pre, slong len, slong prec)
void acb_dft_precomp_clear(acb_dft_pre_t pre)

Initialise a plan for transforms of length `len`, and free it. The tables
of roots of unity and the fixed-point precision are chosen for the
working precision `prec`, so the plan should be built at the precision
of the transforms that will use it. For `len` $\le 0$ an empty plan is
made, for which the transforms do nothing. Throws an exception if the
plan cannot be built.
:::

::: function
void acb_dft_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)
void acb_dft_inverse_precomp(acb_ptr w, acb_srcptr v, const acb_dft_pre_t pre, slong prec)

Set `w` to the DFT, respectively the inverse DFT, of `v`, both of length
`pre->n`, using the plan. The accuracy is determined by the precision
given to `acb_dft_precomp_init`; the argument `prec` only sets the
rounding of the outputs. A plan is not modified by a transform and may
be shared by several threads.
:::

## Example {#example}

The DFT of $(1, 2, 3, 4, 5)$. The first output is the sum of the
inputs; its imaginary part is a tiny ball around zero, since the
fixed-point transform does not know that it is exactly zero. The inverse
transform recovers the input within the error bounds.

```c example name=acb_dft_basic
#include <flint/acb.h>
#include <flint/acb_dft.h>

int main(void)
{
    slong k, n = 5, prec = 64;
    acb_ptr v = _acb_vec_init(n), w = _acb_vec_init(n);

    for (k = 0; k < n; k++)
        acb_set_ui(v + k, k + 1);

    acb_dft(w, v, n, prec);
    for (k = 0; k < n; k++)
    {
        acb_printn(w + k, 10, ARB_STR_NO_RADIUS);
        flint_printf("\n");
    }

    acb_dft_inverse(w, w, n, prec);
    flint_printf("round trip contains input: %d\n",
                 _acb_vec_contains(w, v, n));

    _acb_vec_clear(v, n);
    _acb_vec_clear(w, n);
    return 0;
}
```
```output
15.00000000 + 0e-33*I
-2.500000000 + 3.440954801*I
-2.500000000 + 0.8122992406*I
-2.500000000 - 0.8122992406*I
-2.500000000 - 3.440954801*I
round trip contains input: 1
```
