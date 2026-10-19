---
title: High and low products
functions: flint_mpn_mulhigh_n flint_mpn_sqrhigh flint_mpn_mullow_n flint_mpn_mul_or_mullow_n flint_mpn_mul_or_mulhigh_n
---

Many algorithms need only part of a product. Newton iteration for the
reciprocal or the square root needs the top half of a product (or the
middle part, when the low part is known to cancel); floating-point
multiplication needs the top limbs of the product of two mantissas; Hensel
division and 2-adic arithmetic need the bottom half. A *short product*
computes such a part in less time than the full product: about half the
time for quadratic multiplication, and a fraction between $1/2$ and 1 for
Karatsuba and Toom–Cook multiplication. The descriptions below state
precisely what the functions compute; the algorithms and the proofs of
the error bounds are on the algorithm page.

::: summary algorithms/integers/short-products :::

## Low products {#low}

The low product $uv \bmod B^n$ of two $n$-limb integers is well defined,
and `flint_mpn_mullow_n` computes it exactly. The function also returns
limb $n$ of the product, $\lfloor uv/B^n \rfloor \bmod B$, which comes
almost for free from the carries; `flint_mpn_mul_or_mullow_n` stores the
same $n + 1$ limbs in an output of $2n$ limbs, which lets it compute the
full product when that is the fastest method (from
`FLINT_MPN_MULHIGH_MUL_CUTOFF` limbs without `fft_small`) or use the
high half as scratch space.

## High products {#high}

The top $n$ limbs $\lfloor uv / B^n \rfloor$ of a $2n$-limb product
depend on all the partial products through the carries, and cannot be
computed exactly faster than the full product in the worst case. A
*high product* therefore computes an approximation. We use three
notions [@HZ2011]:

- a *rough* high product sums enough partial products to produce $n$
  limbs of which the top $n - 1$ are essentially correct, the lowest one
  having an error of $O(n)$ units;
- a *precise* high product sums one more diagonal of partial products and
  produces $n + 1$ limbs, the top $n$ being essentially correct and the
  lowest one having an error of $O(n)$ units;
- the *exact* high product is $\lfloor uv/B^n \rfloor$.

`flint_mpn_mulhigh_n` and `flint_mpn_sqrhigh` compute precise high
products: the $n$ limbs written to `res` and the returned limb $r$ form
an $(n+1)$-limb number $h = \mathit{res} \cdot B + r$ with

$$
\left\lfloor \frac{uv}{B^{n-1}} \right\rfloor - 2n < h \le
\left\lfloor \frac{uv}{B^{n-1}} \right\rfloor .
$$

The result is never too large, so $(\mathit{res}, n)$ is either
$\lfloor uv/B^n \rfloor$ or one less, and it is the exact value unless
$r > B - 2n$, since only then can the missing amount carry into the high
limbs. For random inputs this almost never happens, so in practice one
knows the exact top limbs without further work; when it does happen, the
exact value can be found from a low product, as
`_flint_mpn_mulhigh_known_low` does.
Above `FLINT_MPN_MULHIGH_MUL_CUTOFF` limbs (without `fft_small`) the full
product is computed and the result is exact.

For floating-point arithmetic, `flint_mpn_mulhigh_normalised` also shifts
the result so that its top bit is set.

## Middle products {#middle}

The windowed middle product `flint_mpn_mulmid` generalises both: it
computes any range of limbs $[\mathit{zlo}, \mathit{zhi})$ of a product of
operands of any lengths, as a lower approximation with a deficit of
less than $\min(\mathit{an}, \mathit{bn}, \mathit{zlo}) \cdot B$ units of
its lowest limb, exact for $\mathit{zlo} = 0$. For $n$-limb operands, the low product is the window
$[0, n)$ and the high product is (up to the error bound) the window
$[n - 1, 2n)$; `flint_mpn_mulmid` calls `flint_mpn_mullow_n` and
`flint_mpn_mulhigh_n` for windows of these shapes. The balanced middle product of a $(2n-1)$-limb and
an $n$-limb operand (the window $[n, 2n)$), which is what a Newton step
for the reciprocal needs, can be computed in the time of an $n \times n$
product by the transposition principle [@HQZ2004]; FLINT has it from
GMP (`flint_mpn_mulmid_n`, where available) and, for large operands,
from `fft_small`, which computes such a window as a transposed product
of the length of the shorter operand and the window
([](algorithms/integers/short-products#transposed-fft)), or with a
cyclic convolution.

## Choice of algorithm {#selection}

On x86-64 with ADX, high products of up to 13 limbs, square high products
of up to 8 limbs and low products of up to 8 limbs are hard-coded
assembly routines (on ARM64, high products and square high products up
to 8 limbs). Beyond, high products use a quadratic basecase up to 40
limbs, Mulders' recursive algorithm up to `FLINT_MPN_MULHIGH_FFT_SMALL_CUTOFF`
(500 limbs by default), and above, a middle product computed with
`fft_small`; low products use the basecase up to 50 limbs, then Mulders'
algorithm and `fft_small` in the same way. Without `fft_small`, the
last step is a full product from 2000 limbs.

## Example {#example}

The high product of two 30-limb numbers, checked against the full
product:

```c example name=mpn_mulhigh_bound
#include <gmp.h>
#include <flint/flint.h>
#include <flint/mpn_extras.h>

#define N 30

int main(void)
{
    flint_rand_t state;
    mp_limb_t u[N], v[N], full[2 * N], hi[N + 1], lo[N];
    mp_limb_t diff[N + 1];
    int i, ok_bound = 1, ok_low = 1;

    flint_rand_init(state);

    for (i = 0; i < 1000; i++)
    {
        flint_mpn_rrandom(u, state, N);
        flint_mpn_rrandom(v, state, N);
        flint_mpn_mul_n(full, u, v, N);

        /* hi[0] is the returned limb, hi[1..N] the high limbs */
        hi[0] = flint_mpn_mulhigh_n(hi + 1, u, v, N);

        /* 0 <= floor(uv / B^(N-1)) - h < 2N */
        if (mpn_cmp(full + N - 1, hi, N + 1) < 0)
            ok_bound = 0;
        mpn_sub_n(diff, full + N - 1, hi, N + 1);
        if (!flint_mpn_zero_p(diff + 1, N) || diff[0] >= 2 * N)
            ok_bound = 0;

        /* the low product and its extra limb are exact */
        if (flint_mpn_mullow_n(lo, u, v, N) != full[N]
                || mpn_cmp(lo, full, N) != 0)
            ok_low = 0;
    }

    flint_printf("high product within bound: %d\n", ok_bound);
    flint_printf("low product exact: %d\n", ok_low);

    flint_rand_clear(state);
    return 0;
}
```
```output
high product within bound: 1
low product exact: 1
```
