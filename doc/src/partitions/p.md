---
title: Computing p(n)
functions: partitions_fmpz_fmpz
---

`partitions_fmpz_ui` and `partitions_fmpz_fmpz` compute one value of the
partition function exactly, in time roughly proportional to the size of
the result, $O(n^{1/2 + o(1)})$. They are the functions to use for a
single $p(n)$ with $n$ beyond a few thousand. For all values up to $n$,
`arith_number_of_partitions_vec` (power series inversion, $O(n^{3/2 +
o(1)})$ for the whole vector) is cheaper than $n$ separate calls; for
$p(n)$ modulo a word-size integer and all $n$ up to a bound,
`arith_number_of_partitions_nmod_vec`.

::: summary algorithms/number-theory/partitions :::

## Methods {#methods}

| Range of $n$ | Method |
|---|---|
| $n < 0$ | $p(n) = 0$ |
| $0 \le n < 128$ | table `partitions_lookup` |
| $128 \le n < 417$ | Euler's pentagonal recurrence in single words (64-bit machines), $p(416) < 2^{64}$ |
| $n \ge 417$ | Hardy–Ramanujan–Rademacher formula, `mp_real_partitions_hrr` |

The formula is evaluated as a ball whose radius covers the rounding
errors of every term and Rademacher's bound for the truncation error;
the function returns the unique integer in the ball
([](algorithms/floating-point/mp-real-balls#unique-integer)). It never returns an
incorrect value: if the ball contained no integer or several, it would
throw, which is not expected to happen. The terms are computed in parallel
from about $n = 10^8$ when threads are enabled with
`flint_set_num_threads`.

The limits are those of the data types: `partitions_fmpz_fmpz` throws for
$n \ge 2^{128}$, and the result must fit in an `fmpz`, whose size is
limited to about $2^{37}$ bits on 64-bit machines, which holds up to $n
\approx 1.4 \cdot 10^{21}$. The memory needed is a small multiple of the
size of the result, which is about $0.46 \sqrt n$ bytes.

## Example {#example}

The program prints $p(100)$ and $p(1000)$, computes $p(10^9)$ and
compares it with the leading term of the series. The relative difference
is of the order of $e^{-t/2}$, where $t = \pi\sqrt{24 n - 1}/6$, which for
$n = 10^9$ is about $10^{-17600}$; at 128 bits the ratio is 1 up to the
rounding error.

```c example name=partitions_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/arb.h>
#include <flint/partitions.h>

int main(void)
{
    fmpz_t p, n;
    arb_t lead, x;

    fmpz_init(p);
    fmpz_init(n);
    arb_init(lead);
    arb_init(x);

    partitions_fmpz_ui(p, 100);
    flint_printf("p(100) = ");
    fmpz_print(p);
    partitions_fmpz_ui(p, 1000);
    flint_printf("\np(1000) = ");
    fmpz_print(p);
    flint_printf("\n");

    fmpz_set_ui(n, 1000000000);
    partitions_fmpz_fmpz(p, n, 0);
    flint_printf("p(10^9) has %wu bits\n", fmpz_bits(p));

    partitions_leading_fmpz(lead, n, 128);
    arb_set_fmpz(x, p);
    arb_div(x, x, lead, 128);
    flint_printf("p(10^9) / leading term = ");
    arb_printn(x, 25, 0);
    flint_printf("\n");

    fmpz_clear(p);
    fmpz_clear(n);
    arb_clear(lead);
    arb_clear(x);
    return 0;
}
```
```output
p(100) = 190569292
p(1000) = 24061467864032622473692149727991
p(10^9) has 116993 bits
p(10^9) / leading term = [1.000000000000000000000000 +/- 9.48e-38]
```
