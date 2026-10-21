---
title: Exponentiation with a precomputed inverse
functions: n_powmod2_preinv n_powmod2_ui_preinv n_powmod2_fmpz_preinv
---

These functions compute $a^e \bmod n$ for a reduced base $a < n$ and a
word-size modulus $n$, using the `ulong` inverse from `n_preinvert_limb`
(`n_powmod2` computes it itself). They differ in the type of the exponent:
`slong` for `n_powmod2_preinv` and `n_powmod2`, where a negative $e$
computes a power of $a^{-1} \bmod n$; `ulong` for the full unsigned range;
and `fmpz_t` for nonnegative exponents of any size. The functions with a
`double` inverse (`n_powmod_precomp` and its relatives) are documented separately.

By convention $a^0 = 1$ (also for $a = 0$), reduced modulo $n$, so the
result is 0 when $n = 1$.

## Normalised moduli {#normalised}

`n_powmod_ui_preinv` works like `n_mulmod_preinv`: the modulus is shifted
left until its top bit is set, and the base is shifted by the same number
of bits $s$. For any $n$,
`n_powmod_ui_preinv(a << s, e, n << s, ninv, s) >> s` equals
`n_powmod2_ui_preinv(a, e, n, ninv)`, where $s$ is the number of leading zero
bits of $n$; `n_powmod2_preinv` is implemented this way. Calling it directly
saves the shifts when the operands are already kept in shifted form.

## Algorithm {#algorithm}

Right-to-left binary exponentiation: the base is squared once for each bit
of the exponent below the top bit, and multiplied into the result once for
each set bit after the lowest. Each step is one modular multiplication with
a [Möller–Granlund](algorithms/arithmetic/moller-granlund) reduction, so
the cost is at most $2 \log_2 e$ modular multiplications.

## Example {#example}

```c example name=powmod_basic
#include <flint/flint.h>
#include <flint/ulong_extras.h>

int main(void)
{
    ulong n = 1000003, ninv = n_preinvert_limb(n);
    ulong x, norm;

    /* Fermat: 2^(n-1) = 1 mod n for prime n */
    flint_printf("2^(n-1) mod n = %wu\n", n_powmod2_preinv(2, n - 1, n, ninv));

    /* a negative exponent gives a modular inverse */
    x = n_powmod2_preinv(3, -1, n, ninv);
    flint_printf("3^-1 mod n = %wu, check: %wu\n", x, n_mulmod2_preinv(3, x, n, ninv));

    /* n_powmod_ui_preinv works with the shifted (normalised) modulus */
    norm = flint_clz(n);
    x = n_powmod_ui_preinv(UWORD(5) << norm, 1000, n << norm, ninv, norm) >> norm;
    flint_printf("5^1000 mod n = %wu = %wu\n", x, n_powmod2_ui_preinv(5, 1000, n, ninv));
    return 0;
}
```
```output
2^(n-1) mod n = 1
3^-1 mod n = 666669, check: 1
5^1000 mod n = 463979 = 463979
```

::: warning Shift with a ulong
Write `UWORD(5) << norm`, not `5 << norm`: the literal 5 is an `int`, and
shifting it by 31 bits or more is undefined behaviour, which in practice
gives wrong results.
:::
