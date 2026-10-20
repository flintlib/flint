---
title: Complete factorisation
functions: n_factor_t n_factor_init n_factor
---

`n_factor` is the function to use for factoring a word-size integer. It
stores the result in an `n_factor_t`, which has room for the at most 15
distinct prime factors of a 64-bit integer. The primes are not necessarily
in increasing order.

The structure has fixed size and needs no memory management:
`n_factor_init` only sets `num` to zero, and there is no clear function.
`n_factor` overwrites its output, so the initialisation is only needed
before functions that add to a factorisation, such as `n_factor_trial`
and `n_factor_insert`.

## Algorithm {#algorithm}

After removing the power of 2, `n_factor` trial-divides by the primes below
$2^{16}$, testing the cofactor for primality early on so that numbers with
a large prime factor are not trial-divided needlessly. Any remaining
composite has two or three prime factors above $2^{16}$. These are split by
Pollard–Brent rho, Hart's one-line method or SQUFOF, depending on the size,
after checking for perfect squares, cubes and fifth powers.

::: summary algorithms/number-theory/factoring-word :::

Integers of more than one word are factored by `fmpz_factor`; the methods
for them are compared in [](algorithms/number-theory/integer-factoring).

## Related functions {#related}

- `n_factor_partial` stops once the product of the factors found exceeds a
  limit; `n_factor_trial` only does trial division.
- `n_factor_evaluate` multiplies a factorisation back out, and
  `n_factor_insert` adds a prime power to one.

## Example {#example}

```c example name=factor_basic
#include <flint/flint.h>
#include <flint/ulong_extras.h>

int main(void)
{
    ulong values[] = { 360, 3215031751, 4294967295, 4294967291 };
    n_factor_t fac;
    int i, j;

    for (i = 0; i < 4; i++)
    {
        n_factor_init(&fac);
        n_factor(&fac, values[i], 1);

        flint_printf("%wu =", values[i]);
        for (j = 0; j < fac.num; j++)
        {
            flint_printf(" %wu", fac.p[j]);
            if (fac.exp[j] > 1)
                flint_printf("^%d", fac.exp[j]);
            if (j + 1 < fac.num)
                flint_printf(" *");
        }
        flint_printf("\n");
    }
    return 0;
}
```
```output
360 = 2^3 * 3^2 * 5
3215031751 = 151 * 751 * 28351
4294967295 = 3 * 5 * 17 * 257 * 65537
4294967291 = 4294967291
```
