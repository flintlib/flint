---
title: Multiplication with a precomputed inverse
functions: n_mulmod2_preinv n_mulmod2 n_mulmod_preinv
---

These functions compute $ab \bmod n$ for a word-size modulus $n$. They form
the inner loop of most modular arithmetic in FLINT; the `nmod` module wraps
them together with the modulus and its inverse in an `nmod_t`.

## Choosing a function {#choosing}

- `n_mulmod2_preinv` is the general-purpose function. Compute
  `ninv = n_preinvert_limb(n)` once per modulus.
- `n_mulmod2` computes `ninv` itself and is only worth using for a single
  multiplication.
- `n_mulmod_preinv` saves a few operations by requiring the modulus to be
  normalised (top bit set) and the operands reduced. To use it with any
  $n$, work with $n 2^s$ and operands shifted left by $s$ bits, where $s$ is
  the number of leading zero bits of $n$, and pass $s$ as `norm`. Results
  then also come out shifted by $s$ bits.

`n_mulmod_precomp` uses a `double` inverse instead and requires
$n < 2^{53}$. When one operand is fixed, `n_mulmod_shoup` is faster still;
see [](algorithms/arithmetic/choosing) for all options.

## Algorithm {#algorithm}

The product $ab$ is computed as a two-word integer and reduced with one
step of [Möller–Granlund division](algorithms/arithmetic/moller-granlund):
one full and one low multiplication by the precomputed inverse, followed by
at most two corrections. `n_mulmod_preinv` performs exactly this step,
since its operands are already shifted. `n_mulmod2_preinv` passes the
product to `n_ll_mod_preinv`, which shifts it to normalise the modulus;
if the operands are not reduced, the high word of the product may be at
least $n$, and a second division step reduces it first.

::: summary algorithms/arithmetic/moller-granlund :::

## Examples {#examples}

The most basic use:

```c example name=mulmod_basic
#include <flint/flint.h>
#include <flint/ulong_extras.h>

int main(void)
{
    ulong a = 12345678, b = 87654321, n = 111111111;
    ulong ninv = n_preinvert_limb(n);

    flint_printf("%wu * %wu mod %wu = %wu\n", a, b, n,
                 n_mulmod2_preinv(a, b, n, ninv));
    return 0;
}
```
```output
12345678 * 87654321 mod 111111111 = 23456790
```

The inverse pays off when the modulus is reused. This computes $1000!$
modulo a prime with 999 multiplications and checks the result against
`n_factorial_mod2_preinv`:

```c example name=mulmod_factorial
#include <flint/flint.h>
#include <flint/ulong_extras.h>

int main(void)
{
    ulong p = 4294967291;               /* the largest prime below 2^32 */
    ulong pinv = n_preinvert_limb(p);   /* computed once, used 999 times */
    ulong f = 1, k;

    for (k = 2; k <= 1000; k++)
        f = n_mulmod2_preinv(f, k, p, pinv);

    flint_printf("1000! mod p = %wu\n", f);
    flint_printf("check       = %wu\n", n_factorial_mod2_preinv(1000, p, pinv));
    return 0;
}
```
```output
1000! mod p = 444887038
check       = 444887038
```
