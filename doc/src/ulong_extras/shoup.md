---
title: Shoup multiplication
functions: n_mulmod_precomp_shoup n_mulmod_precomp_shoup_preinv n_mulmod_shoup n_mulmod_precomp_shoup_quo_rem n_mulmod_and_precomp_shoup
---

When many products $a b_i \bmod n$ share the same $a$ and $n$, Shoup's
method [@Sho2021] computes the quotient
$\hat a = \lfloor a 2^{\bits} / n \rfloor$ once, after which each product
needs one high and two low multiplications and a single
comparison, with no reduction step. This is the fastest way to multiply
modulo $n$ in FLINT when it applies: scaling vectors and matrices, and
multiplying by fixed roots of unity in number-theoretic transforms.

## Requirements {#requirements}

- The modulus must satisfy $n < 2^{\bits - 1}$. This is what allows the
  final correction to be a single comparison.
- The fixed operand must be reduced, $a < n$, so that $\hat a$ fits in a
  word.
- The other operand $b$ of `n_mulmod_shoup` may be any word.
  `n_mulmod_and_precomp_shoup` requires $b < n$.

## Computing the precomputed quotient {#precomp}

`n_mulmod_precomp_shoup` computes $\hat a$ with a two-word by one-word
hardware division. Its cost varies widely between processors: it is a
few multiplications' worth on recent x86-64, but tens of cycles on
older ones such as Skylake and Cascade Lake, and 64-bit ARM has no
two-word division instruction (FLINT then divides in software). When the
precomputed inverse of $n$ is at hand, as `mod.ninv` in an `nmod_t`,
`n_mulmod_precomp_shoup_preinv` and
`n_mulmod_precomp_shoup_quo_rem_preinv` give the same results with
multiplications only. On recent x86-64 processors the division is fast
and runs alongside multiplications, so the two versions are close and the
division can even win in code that is busy multiplying.

## Tables of powers {#powers}

To multiply by $a, a^2, a^3, \ldots$ in turn, each power needs its own
precomputed quotient. Computing it with `n_mulmod_precomp_shoup` would cost a
division per power. `n_mulmod_and_precomp_shoup` instead obtains the
quotient of $ab \bmod n$ from the quotient and remainder of $a$ and the
quotient of $b$, using one extra high product and two low products. The
quotient and remainder of $a$ come from `n_mulmod_precomp_shoup_quo_rem`.

::: summary algorithms/arithmetic/shoup :::

## Example {#example}

```c example name=shoup_basic
#include <flint/flint.h>
#include <flint/ulong_extras.h>

int main(void)
{
    ulong n = 1000003, a = 12345;      /* n < 2^(FLINT_BITS - 1), a < n */
    ulong a_pre = n_mulmod_precomp_shoup(a, n);
    ulong v[5] = { 1, 2, 3, 1000002, 4000000000 };
    ulong quo, rem, pw, pw_pre;
    int i;

    /* scale a vector by a; the entries need not be reduced */
    for (i = 0; i < 5; i++)
        v[i] = n_mulmod_shoup(a, v[i], a_pre, n);
    flint_printf("a*v = %wu %wu %wu %wu %wu\n", v[0], v[1], v[2], v[3], v[4]);

    /* powers of a, each with its precomputed quotient, without divisions */
    n_mulmod_precomp_shoup_quo_rem(&quo, &rem, a, n);
    pw = a;
    pw_pre = quo;
    flint_printf("powers:");
    for (i = 2; i <= 6; i++)
    {
        n_mulmod_and_precomp_shoup(&pw, &pw_pre, a, pw, quo, rem, pw_pre, n);
        flint_printf(" %wu", pw);
        if (pw_pre != n_mulmod_precomp_shoup(pw, n))
            flint_printf("(bad quotient)");
    }
    flint_printf("\n");
    return 0;
}
```
```output
a*v = 12345 24690 37035 987658 860447
powers: 398569 319545 771193 349025 700701
```
