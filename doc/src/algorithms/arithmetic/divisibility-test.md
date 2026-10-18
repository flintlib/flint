---
title: Divisibility testing by multiplication
summary: Tests whether an odd $d$ divides $n$ with one multiplication and one comparison, using $d^{-1} \bmod 2^{\bits}$. Also gives the exact quotient for free.
complexity: One low product and one comparison per test; two words of precomputed data per divisor.
order: 30
---

Testing $d \mid n$ by computing $n \bmod d$ is wasteful when only a yes or
no answer is needed. Granlund and Montgomery [@GraMon1994, §9] observed
that for odd $d$ the test reduces to a multiplication modulo $B = 2^{\bits}$.
FLINT uses it for trial division in `n_is_prime` and `n_factor`, and exposes
it as `n_divisible_odd_gm`.

::: theorem {#thm-div}
Let $d$ be odd, $d' = d^{-1} \bmod B$ and $L = \lfloor (B - 1)/d \rfloor$.
For $0 \le n < B$,
$$
d \mid n \iff (n d' \bmod B) \le L.
$$
Moreover, if $d \mid n$ then $n/d = n d' \bmod B$.
:::

::: proof
Since $d'$ is odd, $\phi(n) = n d' \bmod B$ is a bijection of
$\{0, 1, \ldots, B - 1\}$. The multiples of $d$ in this range are $kd$ for
$k = 0, 1, \ldots, L$, and $\phi(kd) = k d d' \bmod B = k$. So $\phi$ maps the
$L + 1$ multiples of $d$ onto $\{0, \ldots, L\}$. Being a bijection, it maps
no other $n$ into this set. The second claim is $\phi(kd) = k$.
:::

For even $d = 2^e d_0$, test that the lowest $e$ bits of $n$ are zero and
that $d_0 \mid n$.

## Use in FLINT {#use}

`n_is_prime` stores the pairs $(d', L)$ for the 34 odd primes from 3 to 149
in a table. For $n < 2^{32}$ the 32 tests after 3 and 5 are written as a
branch-free loop, which compilers vectorise. Trial division in `n_factor`
uses the same test, and when it succeeds divides out the prime with the
exact quotient $n d' \bmod B$ from the theorem.

::: functions n_divisible_odd_gm n_is_prime n_factor :::
