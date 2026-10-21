---
title: Factoring integers
functions: fmpz_factor fmpz_factor_no_trial fmpz_factor_smooth
---

`fmpz_factor` is the function to use for the complete factorisation of an
integer. It always succeeds, given enough time, and returns prime factors whose primality has been *proved*. Its
running time is governed by the second largest prime factor: integers
whose factors are all small are factored almost instantly whatever their
size, while a product of two primes of 45 digits each takes minutes with
the quadratic sieve. `fmpz_factor_smooth` instead looks only for factors
up to a given size, which is what one needs to compute the smooth part of
a number or to test whether it is smooth.

::: summary algorithms/number-theory/integer-factoring :::

## Strategy {#strategy}

`fmpz_factor` proceeds as follows.

1. Single-word inputs are factored by `n_factor`. Otherwise the power of 2
   is removed.
2. **Trial division** by the primes in `flint_primes`, in blocks of 1000
   primes starting from 3 (`flint_mpn_factor_trial`). Whenever a prime is
   found, its full power is removed and the next block of 1000 primes
   after it is searched, so that highly composite numbers such as
   factorials are factored by trial division alone. As soon as a block
   yields nothing, the cofactor goes to step 3. If it shrinks to one word,
   `n_factor` finishes.
3. `fmpz_factor_no_trial` on the cofactor $m$: if $m$ is prime (proved by
   `fmpz_is_prime`), it is a factor. If $m$ is a perfect power $r^k$
   (`fmpz_is_perfect_power`, see [](fmpz/roots)), $r$ is factored
   recursively.
4. Otherwise `_fmpz_factor_smooth` in adaptive mode with the bound
   $\max(b/3 - 17, 2)$ bits for a $b$-bit $m$: trial division by the
   first 3512 primes (product tree), then the **elliptic curve method**
   with parameters for factors of increasing size. Every time a factor is
   found, the bound is lowered to the value for the new cofactor, so that
   the search for factors stops once the [quadratic sieve](qsieve) is
   expected to finish the cofactor faster.
5. A remaining composite cofactor $c$ is split by:
   - for $c < 2^{128}$ (on 64-bit machines), Pollard–Brent rho in two-word
     arithmetic (`n_ll_factor_rho`, 2048 iterations up to 72 bits and
     8192 above), and from 118 bits a few ECM curves (`n_ll_factor_ecm`
     with $B_1 = 1000$, $B_2 = 10^5$);
   - for $c < 2^{72}$, SQUFOF (`n_ll_factor_SQUFOF`, 50000 or 150000
     iterations) if rho failed;
   - otherwise, or if these fail, the self-initialising quadratic sieve
     (`qsieve_factor`).
6. The pieces are factored recursively with `fmpz_factor_no_trial`.

The perfect power test before the sieve is necessary: the quadratic sieve
cannot factor perfect powers. When rho or SQUFOF returns a composite
factor sharing a prime with its cofactor (as $21$ and $35$ for $735$),
the split is discarded and the sieve is used instead, since the pieces
are combined without merging equal primes.

All prime factors are proved prime by `fmpz_is_prime`, which uses
$n \pm 1$ tests and falls back on
[APRCL or ECPP](algorithms/number-theory/primality-proving). For large
prime factors this proof can take longer than finding them.

The factors are *not* sorted: those found in steps 1–2 come first in
increasing order, and the others follow in the order in which they were
split off.

## Expected running times {#timings}

For a product of two random primes of equal size, `fmpz_factor` took the
following times on one core of the x86-64 machine used to write this
documentation: 8 ms at 30 digits, 27 ms at 40 digits, 0.3 to 0.5 s at 50
digits and 2.7 s at 60 digits, almost all of it in the quadratic sieve,
whose running time grows by a factor of 6 to 10 for every 10 digits in
this range. When one factor is much smaller than the others, the elliptic
curve method finds it in a time that depends mainly on the size of that
factor: factors of up to 20 digits are typically found in a fraction of a
second even in numbers of hundreds of digits. With several threads, the
elliptic curve stage runs the curves in parallel.

## Smooth parts {#smooth}

`fmpz_factor_smooth(factor, n, bits, proved)` stops after the elliptic
curve stage of step 4 with the given bound, and returns 0 if a composite
cofactor remains. With `proved` = 0 only probable primality is checked,
which is much faster for large factors. The bound is a size, not a
guarantee: a factor of about `bits` bits is found with probability at
least about $2/3$, and smaller factors with higher probability.

## Example {#example}

```c example name=fmpz_factor_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>

int main(void)
{
    fmpz_t n;
    fmpz_factor_t f;
    int done;

    fmpz_init(n);
    fmpz_factor_init(f);

    /* the Fermat number F7 = 2^128 + 1 */
    fmpz_ui_pow_ui(n, 2, 128);
    fmpz_add_ui(n, n, 1);
    fmpz_factor(f, n);
    fmpz_factor_print(f);
    flint_printf("\n");

    /* -(2^64 - 1)^2 * 3^5 */
    fmpz_set_ui(n, UWORD_MAX);
    fmpz_mul(n, n, n);
    fmpz_mul_ui(n, n, 243);
    fmpz_neg(n, n);
    fmpz_factor(f, n);
    fmpz_factor_print(f);
    flint_printf("\n");

    /* factors of 100! + 1 up to about 32 bits */
    fmpz_fac_ui(n, 100);
    fmpz_add_ui(n, n, 1);
    done = fmpz_factor_smooth(f, n, 32, 0);
    flint_printf("%wd factors, complete: %d\n", f->num, done);
    fmpz_print(f->p + 0); flint_printf(" ");
    fmpz_print(f->p + 1); flint_printf(" ");
    fmpz_print(f->p + 2); flint_printf("\n");

    fmpz_factor_clear(f);
    fmpz_clear(n);
    return 0;
}
```
```output
59649589127497217 * 5704689200685129054721
-1 * 3^7 * 5^2 * 17^2 * 257^2 * 641^2 * 65537^2 * 6700417^2
4 factors, complete: 0
101 14303 149239
```

The last entry of the smooth factorisation of $100! + 1$ is a composite
cofactor of 147 digits.
