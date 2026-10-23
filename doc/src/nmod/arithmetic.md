---
title: Basic arithmetic modulo n
summary: How to set up an `nmod_t`, which operations need reduced operands, how products are reduced, and how to divide when the divisor may not be invertible.
order: 1
---

An element of $\Z/n\Z$ is a `ulong` in $[0, n)$, and the modulus travels
with it as an `nmod_t`. The modulus is set up once with `nmod_init`, which
computes the precomputed inverse needed for fast reduction. The other
functions take the `nmod_t` by value, and the basic operations are
inline, so a
call such as `nmod_mul(a, b, mod)` compiles to a few multiplications with
no function call and no division.

## Usage {#usage}

1. Call `nmod_init(&mod, n)` once for each modulus. It computes a
   reciprocal of $n$, which costs about as much as one hardware division.
2. Reduce the inputs with `nmod_set_ui` or `nmod_set_si`. They accept any
   value; most other functions require reduced operands.
3. Compute with `nmod_add`, `nmod_sub`, `nmod_neg`, `nmod_mul`,
   `nmod_addmul`, `nmod_pow_ui`, `nmod_inv` and `nmod_div`. Every result
   is reduced again.

The reduced-operand requirement is not checked. Passing an unreduced value
to `nmod_add` or `nmod_mul` silently gives a wrong result. Use
`nmod_ui_add_ui`, `nmod_ui_mul_ui` and `nmod_ui_pow_ui` when the operands
may be unreduced.

Any modulus $1 \le n < 2^{\bits}$ is allowed. For $n = 1$ every result is
0. A few functions with a leading underscore are faster but restrict $n$:
`_nmod_add` and `_nmod_sub` need $n \le 2^{\bits - 1}$, and
`_nmod_mul_fullword` needs $n \ge 2^{\bits - 1}$.

## How it works {#algorithm}

`nmod_init` stores $n$, the number `norm` of leading zero bits of $n$, and
the reciprocal `ninv` of the normalised modulus $n 2^{\mathit{norm}}$.
Products are reduced with one step of
[Möller–Granlund division](algorithms/arithmetic/moller-granlund):
`nmod_mul` shifts one operand left by `norm` bits, multiplies, reduces the
two-word product modulo $n 2^{\mathit{norm}}$, and shifts the remainder
back. Addition and subtraction use a comparison and need no
multiplication.

`nmod_pow_ui` uses binary exponentiation. For odd $n$ and exponents of
at least $2^{11}$ it converts to
[Montgomery form](algorithms/arithmetic/montgomery), where each
multiplication has a shorter chain of dependent operations. The guide
[](nmod/redc) explains how to use Montgomery form directly for other long
chains of multiplications.

::: summary algorithms/arithmetic/moller-granlund :::

## Inverses and division {#division}

`nmod_inv` and `nmod_div` require the divisor to be invertible, that is,
coprime to $n$. Otherwise they throw an exception (`FLINT_IMPINV`). Each
call runs the extended Euclidean algorithm, which is much slower than a
multiplication, so divisions by the same value should be replaced by one
inversion and several multiplications.

When $c$ may share a factor with $n$, use `nmod_divides`. It decides
whether $ac \equiv b \pmod n$ has a solution $a$ and, if so, returns one.
A solution exists exactly when $\gcd(c, n)$ divides $b$.

## Examples {#examples}

Basic arithmetic modulo a prime:

```c example name=nmod_basic
#include <flint/flint.h>
#include <flint/nmod.h>

int main(void)
{
    nmod_t mod;
    ulong a, b, x;

    nmod_init(&mod, 1000003);           /* a prime */

    a = nmod_set_ui(123456789, mod);    /* reduce the inputs once */
    b = nmod_set_si(-5, mod);
    flint_printf("a = %wu, b = %wu\n", a, b);

    flint_printf("a + b = %wu\n", nmod_add(a, b, mod));
    flint_printf("a - b = %wu\n", nmod_sub(a, b, mod));
    flint_printf("a * b = %wu\n", nmod_mul(a, b, mod));
    flint_printf("a^(n-1) = %wu\n", nmod_pow_ui(a, mod.n - 1, mod));

    x = nmod_inv(a, mod);
    flint_printf("1/a = %wu, check: %wu\n", x, nmod_mul(a, x, mod));
    flint_printf("b/a = %wu\n", nmod_div(b, a, mod));
    return 0;
}
```
```output
a = 456420, b = 999998
a + b = 456415
a - b = 456425
a * b = 717909
a^(n-1) = 1
1/a = 301287, check: 1
b/a = 493571
```

Division by a non-invertible element with `nmod_divides`. Modulo 12, the
equation $8a = b$ is solvable exactly when $\gcd(8, 12) = 4$ divides $b$:

```c example name=nmod_divides
#include <flint/flint.h>
#include <flint/nmod.h>

int main(void)
{
    nmod_t mod;
    ulong b, x;

    nmod_init(&mod, 12);

    for (b = 0; b < 12; b++)
    {
        if (nmod_divides(&x, b, 8, mod))
            flint_printf("8 * %wu = %wu mod 12\n", x, b);
    }
    return 0;
}
```
```output
8 * 0 = 0 mod 12
8 * 11 = 4 mod 12
8 * 10 = 8 mod 12
```

Each solvable equation has four solutions here; for $b = 4$ they are 2, 5,
8 and 11, and `nmod_divides` returns one of them.
