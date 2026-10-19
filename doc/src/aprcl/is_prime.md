---
title: The APRCL primality test
functions: aprcl_is_prime
---

`aprcl_is_prime` proves or disproves the primality of $n$ with the Jacobi
sum test. Unlike a probable prime test, its answer is always correct, but
no certificate is produced; use [`ecpp_prove`](ecpp/prove) if one is
needed, and above about 20000 bits, where ECPP is usually much faster. For general
use, `fmpz_is_prime` is preferable: it disposes of small inputs and of
inputs with easily factored $n \pm 1$ first, and it chooses between APRCL
and ECPP depending on the size and the number of threads.

::: summary algorithms/number-theory/primality-proving :::

## Cost {#cost}

The cost is dominated by the Jacobi sum tests, one for each pair $(p, q)$
with $q \mid s$ and $p \mid q - 1$, which consist of powerings in rings
$\Z[\zeta_{p^k}]/(n)$ of degree up to a few hundred. The number of pairs
and their degrees grow slowly with $n$ through the parameter $R$ (from
`aprcl_R_value`): $R = 5040$ up to 344 bits, $R = 720720$ up to 1566
bits. The tests are independent and are distributed over the threads set
with `flint_set_num_threads`, so APRCL scales well with threads. On one
core of the machine used for these timings, $10^{100} + 267$ took 29 ms
and $2^{521} - 1$ took 0.15 s.

## Example {#example}

```c example name=aprcl_is_prime_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/aprcl.h>

int main(void)
{
    fmpz_t n;
    fmpz_init(n);

    /* the first prime after 10^100 */
    fmpz_ui_pow_ui(n, 10, 100);
    fmpz_add_ui(n, n, 267);
    flint_printf("10^100 + 267: %d\n", aprcl_is_prime(n));

    /* the Mersenne prime 2^127 - 1 */
    fmpz_ui_pow_ui(n, 2, 127);
    fmpz_sub_ui(n, n, 1);
    flint_printf("2^127 - 1: %d\n", aprcl_is_prime(n));

    /* the Fermat number F7 is composite */
    fmpz_ui_pow_ui(n, 2, 128);
    fmpz_add_ui(n, n, 1);
    flint_printf("2^128 + 1: %d\n", aprcl_is_prime(n));

    fmpz_clear(n);
    return 0;
}
```
```output
10^100 + 267: 1
2^127 - 1: 1
2^128 + 1: 0
```
