---
title: Square roots modulo a prime
functions: nmod_sqrt _nmod_sqrt_small
---

`nmod_sqrt` is the function to use for square roots modulo a word-size
prime; `_nmod_sqrt_small` is the exhaustive search it falls back to for
small odd moduli. Of the two square roots $r$ and $n - r$ of a nonzero
square modulo a prime $n$, both functions return the one in $[0, n/2]$, so the
result does not depend on the algorithm that found it.

## Composite and even moduli {#composite}

`nmod_sqrt` does not check that $n$ is prime. For an arbitrary modulus
it behaves as follows.

- For $a = 0$ and $a = 1$ it returns 1 and sets `*res` to $a$, for every
  $n$.
- For even $n$ and $a \ge 2$ it returns 0, even when $a$ is a square.
- For odd $n < 128$ it calls `_nmod_sqrt_small`, which is exact for every
  odd $n$: it returns the least root, or 0 if there is none.
- For odd composite $n \ge 128$, a return value of 1 comes with a genuine
  root (its square is $a$ modulo $n$), but it need not be the least root,
  and a return value of 0 does not prove that $a$ is a non-square.

To find all square roots modulo a composite number, factor it and use
`n_sqrtmodn` from `ulong_extras`.

## Algorithm {#algorithm}

For odd $n < 128$, the function searches through
$t = 1, \ldots, (n-1)/2$, computing $t^2$ incrementally. Otherwise all
arithmetic is done in
[Montgomery form](algorithms/arithmetic/montgomery), with half-word,
full-word lazy or full-word arithmetic depending on the size of $n$, and
the method depends on $n$:

| Case | Method | Cost |
|---|---|---|
| $n \equiv 3 \pmod 4$ | $r = a^{(n+1)/4}$ | one exponentiation |
| $n \equiv 5 \pmod 8$ | Atkin's formula | one exponentiation |
| $n \equiv 1 \pmod 8$, $s < 12$ | Tonelli–Shanks | one exponentiation and $O(s^2)$ dependent multiplications |
| $n \equiv 1 \pmod 8$, $s \ge 12$ | Tonelli–Shanks by bits | one exponentiation and $O(s^2)$ independent multiplications |
| $n \equiv 1 \pmod 8$, large $s$ | Müller's Lucas sequence method | about one exponentiation and a few Jacobi symbols |

Here $s$ is the exponent of the largest power of 2 dividing $n - 1$.
Müller's method is used when $6s$ is at least the bit length of $n$ plus
120, which requires $n \ge 2^{25}$; it is typical of primes chosen for
number-theoretic transforms, such as $469762049 = 7 \cdot 2^{26} + 1$.

None of the methods tests whether $a$ is a square first. Each produces a
candidate root, which is accepted only if its square is $a$. For prime
$n$ this check is the residuosity test, and it is also what makes the
result trustworthy for composite $n$.

::: summary algorithms/number-theory/sqrt-mod-p :::

## Example {#example}

The four primes below exercise the four main code paths: $n \equiv 3
\pmod 4$, $n \equiv 5 \pmod 8$, Tonelli–Shanks with $s = 5$ (where 10 is
not a square), and Müller's method with $s = 26$.

```c example name=nmod_sqrt_basic
#include <flint/flint.h>
#include <flint/nmod.h>

int main(void)
{
    ulong primes[4] = { 1000003, 1000037, 1000033, 469762049 };
    nmod_t mod;
    ulong r, a;
    int i;

    for (i = 0; i < 4; i++)
    {
        nmod_init(&mod, primes[i]);
        a = nmod_set_ui(10, mod);

        if (nmod_sqrt(&r, a, mod))
            flint_printf("n = %wu: sqrt(10) = %wu, r^2 = %wu\n",
                         mod.n, r, nmod_mul(r, r, mod));
        else
            flint_printf("n = %wu: 10 is not a square\n", mod.n);
    }

    /* composite n: the root found is genuine */
    nmod_init(&mod, 135);
    if (nmod_sqrt(&r, 81, mod))
        flint_printf("n = 135: sqrt(81) = %wu\n", r);
    return 0;
}
```
```output
n = 1000003: sqrt(10) = 394215, r^2 = 10
n = 1000037: sqrt(10) = 122489, r^2 = 10
n = 1000033: 10 is not a square
n = 469762049: sqrt(10) = 115267370, r^2 = 10
n = 135: sqrt(81) = 54
```

Modulo 135, the square roots of 81 in $[0, 67]$ are 9, 36 and 54. The
function returns 54, a genuine root but not the least one.
