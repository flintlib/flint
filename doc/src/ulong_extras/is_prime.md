---
title: Deterministic primality test
functions: n_is_prime n_is_prime_odd_no_trial
---

`n_is_prime` is the function to use for testing whether a word-size integer
is prime. Its answer is always correct: it is not a probabilistic test,
even though it is built from probable prime tests, because every base-2
strong pseudoprime below $2^{64}$ is known and handled.

`n_is_prime_odd_no_trial` performs the steps below except the trial
division (step 2). This pays off only when trial division is known to be
useless, for example for numbers produced by a sieve; for random input it
is much slower than `n_is_prime`.

## Algorithm {#algorithm}

1. Look up small $n$ in tables: a bit table for $n < 2^{15}$, and a table
   modulo 30 for $n < 2^{20}$ after removing factors 2, 3 and 5.
2. Trial-divide by the primes from 7 to 149 using
   [multiplication by inverses](algorithms/arithmetic/divisibility-test).
3. Below $2^{32}$: a strong probable prime test to base 2, then a lookup
   in the table of the 2314 base-2 strong pseudoprimes in this range.
4. Above $2^{32}$: a base-2 test, a check for a family of semiprimes that
   have few witnesses, and a second test with a base read from a hashed
   table of witnesses.

## Cost {#cost}

For random input most numbers are rejected by trial division, which is
cheap and uses no division instructions. A number that survives
trial division costs one strong probable prime test below $2^{32}$, and at
most two above. Composites usually fail the first test, so in practice the
second test is mostly paid by primes.

::: summary algorithms/number-theory/primality-word :::

## Example {#example}

```c example name=is_prime_basic
#include <flint/flint.h>
#include <flint/ulong_extras.h>

int main(void)
{
    ulong n, count = 0;
    ulong tests[] = { 2147483647, 3215031751, 4294967291, 4294967295 };
    int i;

    for (n = 0; n < 1000000; n++)
        count += n_is_prime(n);
    flint_printf("primes below 10^6: %wu (n_prime_pi gives %wu)\n",
                 count, n_prime_pi(999999));

    for (i = 0; i < 4; i++)
        flint_printf("%wu: %s\n", tests[i],
                     n_is_prime(tests[i]) ? "prime" : "composite");
    return 0;
}
```
```output
primes below 10^6: 78498 (n_prime_pi gives 78498)
2147483647: prime
3215031751: composite
4294967291: prime
4294967295: composite
```

The composite $3215031751 = 151 \cdot 751 \cdot 28351$ in this example is
a strong pseudoprime to the bases 2, 3, 5 and 7 simultaneously, so a test
with only these bases would wrongly report it prime.
