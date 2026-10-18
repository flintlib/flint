---
title: Using Montgomery arithmetic
summary: When Montgomery form pays off, how to convert values in and out, the rules for lazy residues and half-word arithmetic, and how the reduction works.
order: 3
---

Montgomery arithmetic represents $x \in \Z/n\Z$ by $x R \bmod n$, where $R$
is a power of two. In this form a modular product needs no division, only
a reduction modulo $R$ (a mask) and a division by $R$ (a shift). It
requires an odd modulus, and values must be converted into Montgomery
form and back. It is worth using for long sequences of multiplications,
such as exponentiation, primality tests and polynomial arithmetic over
$\Z/n\Z$; `nmod_pow_ui` and `nmod_sqrt` use it internally for this
reason. This guide describes the `nmod_redc` functions listed under
[](nmod/index#redc).

## Usage {#usage}

1. Initialise a context with `nmod_redc_ctx_init_ui` or, if an `nmod_t`
   is already available, `nmod_redc_ctx_init_nmod`. The modulus must be
   odd. The context needs no clearing.
2. Convert the operands with `nmod_redc_set_nmod` (for reduced $x$) or
   `nmod_redc_set_ui` (for any $x$).
3. Compute with `nmod_redc_add`, `nmod_redc_sub`, `nmod_redc_neg`,
   `nmod_redc_mul` and, for powers with exponent at least 1,
   `_nmod_redc_pow_ui`. Constants must also be converted: the Montgomery
   form of 1 is `nmod_redc_set_ui(1, ctx)`, which is not 1.
4. Convert results back with `nmod_redc_get_nmod`.

Addition, subtraction and negation are the ordinary operations modulo $n$,
because the map $x \mapsto xR$ is additive. Comparison with 0 and equality
tests also work directly on Montgomery forms in $[0, n)$.

## Lazy residues {#fast}

When $n < 2^{\bits - 2}$, which `nmod_redc_can_use_fast` tests, the `fast`
functions can be used instead. They accept and return residues in
$[0, 2n)$, which saves one comparison and subtraction per multiplication.
The rules are:

- `fast` and standard functions use the same context and the same
  Montgomery form, so values in $[0, n)$ can be passed to `fast` functions.
- Before passing a value from a `fast` function to a standard function,
  or comparing two values, reduce it to $[0, n)$ with
  `nmod_redc_fast_normalise`.
- `nmod_redc_get_nmod` accepts any word, so it can convert lazy residues
  directly.

## Half-word arithmetic {#half}

The `half` functions use $R = 2^{\bits/2}$, so that a product of two
residues and the reduction fit in one word. This avoids two-word
arithmetic, but they are not necessarily faster than the full-word
functions. They need their own context, initialised with
`nmod_redc_half_ctx_init_ui` or `nmod_redc_half_ctx_init_nmod`, and their
Montgomery form $x 2^{\bits/2} \bmod n$ differs from that of the full-word
functions, so the two kinds cannot be mixed. The restrictions are:

| Functions | Maximum $n$ (64-bit) | Maximum $n$ (32-bit) |
|---|---|---|
| `nmod_redc_` | $2^{64} - 1$ | $2^{32} - 1$ |
| `nmod_redc_fast_` | $2^{62} - 1$ | $2^{30} - 1$ |
| `nmod_redc_half_` | $2^{31} - 1$ | $2^{15} - 1$ |
| `nmod_redc_half_fast_` | $2^{30} - 1$ | $2^{14} - 1$ |

`nmod_redc_half_can_use_fast` tests the condition for the `half_fast`
functions.

## Algorithm {#algorithm}

A product $T = ab$ is reduced by adding the multiple $mn$ of $n$, with
$m = T n' \bmod R$ and $n' = -n^{-1} \bmod R$, which makes $T + mn$
divisible by $R$. The quotient $(T + mn)/R$ is congruent to $T R^{-1}$ and
lies in $[0, 2n)$. The standard functions use a variant with a subtraction
that also works for $n$ close to $2^{\bits}$; the `fast` functions skip the
final correction. The inverse $n'$ is computed once in the context by
`n_binvert`, and conversion into Montgomery form uses the `nmod_t` stored
in the context. The context functions are thin wrappers of low-level
functions such as `n_mulmod_redc`, `n_redc` and `n_to_redc_preinv`,
which take $n$ and $n'$ directly.

::: summary algorithms/arithmetic/montgomery :::

## Examples {#examples}

Repeated multiplication in Montgomery form, with standard and lazy
residues:

```c example name=nmod_redc_basic
#include <flint/flint.h>
#include <flint/nmod.h>

int main(void)
{
    nmod_redc_ctx_t ctx;
    ulong n = 1000003, a = 12345, k;
    ulong x, ar, xr;

    nmod_redc_ctx_init_ui(ctx, n);          /* n must be odd */

    /* a^1000 mod n, keeping the running product in Montgomery form */
    ar = nmod_redc_set_nmod(a, ctx);         /* a R mod n */
    xr = nmod_redc_set_ui(1, ctx);           /* 1 R mod n */
    for (k = 0; k < 1000; k++)
        xr = nmod_redc_mul(xr, ar, ctx);
    x = nmod_redc_get_nmod(xr, ctx);

    flint_printf("a^1000 = %wu\n", x);
    flint_printf("check  = %wu\n", nmod_pow_ui(a, 1000, ctx->mod));

    /* the same with lazy residues in [0, 2n) */
    if (nmod_redc_can_use_fast(ctx))
    {
        xr = nmod_redc_set_ui(1, ctx);
        for (k = 0; k < 1000; k++)
            xr = nmod_redc_fast_mul(xr, ar, ctx);
        xr = nmod_redc_fast_normalise(xr, ctx);
        flint_printf("fast   = %wu\n", nmod_redc_get_nmod(xr, ctx));
    }
    return 0;
}
```
```output
a^1000 = 294199
check  = 294199
fast   = 294199
```

A Fermat test $3^{n-1} \equiv 1 \pmod n$ with half-word arithmetic. The
moduli are below $2^{15}$, so the program also works on 32-bit machines.
32757 is composite, the other two are prime:

```c example name=nmod_redc_half
#include <flint/flint.h>
#include <flint/nmod.h>

int main(void)
{
    nmod_redc_ctx_t ctx;
    ulong ns[3] = { 32749, 32757, 30011 };
    ulong n, e, x, b;
    int i;

    for (i = 0; i < 3; i++)
    {
        n = ns[i];
        nmod_redc_half_ctx_init_ui(ctx, n);

        b = nmod_redc_half_set_ui(3, ctx);
        x = nmod_redc_half_set_ui(1, ctx);
        for (e = n - 1; e != 0; e >>= 1)       /* right-to-left binary powering */
        {
            if (e & 1)
                x = nmod_redc_half_mul(x, b, ctx);
            b = nmod_redc_half_mul(b, b, ctx);
        }

        flint_printf("n = %wu: 3^(n-1) = %wu\n", n, nmod_redc_half_get_nmod(x, ctx));
    }
    return 0;
}
```
```output
n = 32749: 3^(n-1) = 1
n = 32757: 3^(n-1) = 3840
n = 30011: 3^(n-1) = 1
```
