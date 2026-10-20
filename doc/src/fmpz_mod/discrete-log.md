---
title: Discrete logarithms modulo a multiprecision prime
functions: fmpz_mod_discrete_log_pohlig_hellman_init fmpz_mod_discrete_log_pohlig_hellman_precompute_prime fmpz_mod_discrete_log_pohlig_hellman_primitive_root fmpz_mod_discrete_log_pohlig_hellman_run fmpz_next_smooth_prime
---

These functions solve $\alpha^x \equiv y \pmod p$ for a prime $p$ of any
size and a fixed primitive root $\alpha$, by the method of Pohlig and
Hellman [@PohHel1978]: the logarithm is computed modulo each prime power
$q^e$ dividing $p - 1$, one base-$q$ digit at a time with Shanks'
baby-step giant-step search in the subgroup of order $q$, and the
results are combined by the Chinese remainder theorem. This is fast when
every prime factor $q$ of $p - 1$ is small; the cost grows like
$\sqrt{q}$ for the largest $q$. The method and its analysis are
described in [](algorithms/number-theory/discrete-log); the word-size
version is `nmod_discrete_log_pohlig_hellman_run`.

## Usage {#usage}

1. Initialise an `fmpz_mod_discrete_log_pohlig_hellman_t` with
   `fmpz_mod_discrete_log_pohlig_hellman_init`.
2. Call `fmpz_mod_discrete_log_pohlig_hellman_precompute_prime` with the
   prime $p$. It factors $p - 1$, finds the least primitive root
   $\alpha$, which `fmpz_mod_discrete_log_pohlig_hellman_primitive_root`
   returns, and precomputes the tables. The return value estimates the
   number of multiplications modulo $p$ per logarithm, which a caller
   can use to decide whether $p$ is suitable.
3. Call `fmpz_mod_discrete_log_pohlig_hellman_run` for each nonzero
   residue $y$; the result is the logarithm in $[0, p - 1)$. The
   precomputed object is not modified, so several threads may share it.
4. Free the memory with `fmpz_mod_discrete_log_pohlig_hellman_clear`.

The function does not check that $p$ is prime or that $y$ is nonzero and
reduced. The precomputation throws an exception if $p - 1$ has a prime
factor that does not fit in a word, and the search throws one if it
finds no logarithm, which can only happen for invalid input. The
factorisation of $p - 1$ dominates the precomputation unless $p - 1$ is
known to be smooth; `fmpz_next_smooth_prime` provides primes for which
it is, with $p - 1$ a product of powers of 2, 5, 7, 11, 13, 17, 19 and 23.

The tables for a prime $q$ have at most 100 entries: a giant step is
about $\sqrt{q}$ for $q$ below about $10^4$ and is doubled until the
table is small enough for larger $q$, so the number of baby steps grows
from $\sqrt{q}$ to about $q/50$ for large $q$.

## Example {#example}

The program takes the next smooth prime after $10^{30}$ (a 128-bit prime;
the table of smooth primes has no entries between $2^{64}$ and $2^{127}$),
computes a power of the primitive root and recovers the exponent.

```c example name=fmpz_mod_discrete_log
#include <flint/fmpz.h>
#include <flint/fmpz_mod.h>

int main(void)
{
    fmpz_t p, x, y, e;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_discrete_log_pohlig_hellman_t L;
    double cost;

    fmpz_init(p);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(e);

    /* a prime p > 10^30 with p - 1 smooth */
    fmpz_set_str(x, "1000000000000000000000000000000", 10);
    if (!fmpz_next_smooth_prime(p, x))
        return 1;
    printf("p = "); fmpz_print(p); printf("\n");

    fmpz_mod_discrete_log_pohlig_hellman_init(L);
    cost = fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(L, p);
    printf("base = ");
    fmpz_print(fmpz_mod_discrete_log_pohlig_hellman_primitive_root(L));
    printf(", cost estimate %s 2000\n", cost < 2000 ? "<" : ">=");

    /* y = base^e mod p, then recover e */
    fmpz_mod_ctx_init(ctx, p);
    fmpz_set_str(e, "123456789012345678901234567890", 10);
    fmpz_mod_pow_fmpz(y, fmpz_mod_discrete_log_pohlig_hellman_primitive_root(L), e, ctx);
    fmpz_mod_discrete_log_pohlig_hellman_run(x, L, y);
    printf("log = "); fmpz_print(x); printf("\n");

    fmpz_mod_ctx_clear(ctx);
    fmpz_mod_discrete_log_pohlig_hellman_clear(L);
    fmpz_clear(p);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(e);
    return 0;
}
```

```output
p = 201315753943344863644036821430284967937
base = 3, cost estimate < 2000
log = 123456789012345678901234567890
```

::: summary algorithms/number-theory/discrete-log :::
