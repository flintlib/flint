---
title: Division with a precomputed inverse
functions: n_preinvert_limb n_mod2_preinv n_divrem2_preinv n_div2_preinv n_ll_mod_preinv n_lll_mod_preinv
---

Division by a word $n$ costs tens of cycles in hardware. When $n$ is used
more than a few times, it is faster to compute a reciprocal of $n$ once
with `n_preinvert_limb` and divide with multiplications. The functions
below accept any nonzero $n$ and any dividend, with two exceptions:
`n_lll_mod_preinv` requires its top word to be reduced, and
`n_divrem_preinv_unnorm` requires $n$ not to be normalised.

## Usage {#usage}

The value returned by `n_preinvert_limb` belongs to $n$ and must be passed
together with it. It is the reciprocal of the normalised divisor
$n 2^s$, where $s$ is the number of leading zero bits of $n$. Functions such
as `n_mulmod_preinv` that require a normalised modulus take the same
inverse, so one call serves both kinds of functions.

`n_preinvert_limb_prenorm` computes the same value without the shift when
$n$ is already normalised. `n_divrem_preinv` and `n_divrem_preinv_unnorm`
take the shift $s$ as an argument, saving its recomputation when many
divisions share one $n$.

## Algorithm {#algorithm}

The dividend is shifted left by $s$ bits into two words, and the quotient is
estimated as the high word of $v u_1 + u$, where $v$ is the reciprocal and
$u_1$ the high word of the shifted dividend. The estimate is off by at most
one in each direction, which is fixed with one comparison and, rarely, a
second one. The remainder is shifted back at the end.

`n_ll_mod_preinv` first reduces the high word modulo $n$ (if it is not
already smaller than $n$) and then performs one division step;
`n_lll_mod_preinv` performs two steps.

::: summary algorithms/arithmetic/moller-granlund :::

## Example {#example}

```c example name=preinv_basic
#include <flint/flint.h>
#include <flint/ulong_extras.h>

int main(void)
{
    ulong n = 1000003;
    ulong ninv = n_preinvert_limb(n);   /* one division, done once */
    ulong a, q, r, hi, lo, wrong = 0;

    r = n_divrem2_preinv(&q, 987654321, n, ninv);
    flint_printf("987654321 = %wu * %wu + %wu\n", q, n, r);

    for (a = 0; a < 1000000; a++)
        wrong += (n_mod2_preinv(7919 * a, n, ninv) != (7919 * a) % n);
    flint_printf("disagreements with %%: %wu\n", wrong);

    /* reduce a full double-word product */
    umul_ppmm(hi, lo, UWORD_MAX / 3, UWORD_MAX / 5);
    r = n_ll_mod_preinv(hi, lo, n, ninv);
    flint_printf("same as n_mulmod2_preinv: %s\n",
        r == n_mulmod2_preinv(UWORD_MAX / 3, UWORD_MAX / 5, n, ninv) ? "yes" : "no");
    return 0;
}
```
```output
987654321 = 987 * 1000003 + 651360
disagreements with %: 0
same as n_mulmod2_preinv: yes
```
