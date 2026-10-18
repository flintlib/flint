---
title: Continued fractions
functions: fmpq_get_cfrac fmpq_set_cfrac fmpq_cfrac_bound
---

Every rational number $x$ has a finite simple continued fraction

$$
x = c_0 + \cfrac{1}{c_1 + \cfrac{1}{c_2 + \cfrac{1}{\ddots +
\cfrac{1}{c_{k-1}}}}} = [c_0; c_1, \ldots, c_{k-1}]
$$

with $c_0 = \lfloor x \rfloor \in \Z$ and $c_1, \ldots, c_{k-1} \ge 1$.
The expansion is unique if the last term is required to be at least 2
(for $k \ge 2$); the other expansion replaces $c_{k-1}$ by $c_{k-1} - 1,
1$. `fmpq_get_cfrac` computes the terms of the shorter expansion,
`fmpq_set_cfrac` evaluates a continued fraction, and `fmpq_cfrac_bound`
bounds the length, so that the caller can allocate the output vector.

::: summary algorithms/number-theory/continued-fractions :::

## Computing the terms {#get}

`fmpq_get_cfrac(c, rem, x, n)` writes at most $n$ terms to $c$ (which must
have room for $n$ initialised `fmpz`), sets `rem` to the remainder $r$,
and returns the number $k$ of terms written, so that

$$
x = [c_0; c_1, \ldots, c_{k-2}, c_{k-1} + r], \qquad 0 \le r < 1.
$$

If $r = 0$, the expansion is complete. If $r \ne 0$, then $k = n$, and the
expansion continues with the terms of $1/r$: passing $1/r$ (computed with
`fmpq_inv`) as the next input generates $c_k, c_{k+1}, \ldots$, so that
the terms can be extracted in chunks of any size. For canonical $x$ and
$n > 0$ the remainder is canonical, and $k = n$ or $r = 0$; with $n \ge$
`fmpq_cfrac_bound(x)` the expansion is always complete.

The function uses a half-GCD algorithm
([](algorithms/number-theory/continued-fractions)): it finds the terms
from the leading bits of the numerator and denominator, recursively, and
costs $O(M(b) \log b)$ for $b$-bit inputs, against $O(b^2)$ for the
Euclidean algorithm in `fmpq_get_cfrac_naive`. The terms are written in
batches, so the time to find the first few terms of a huge fraction is
small.

With $n \le 0$ no term is written and `rem` is set to $1/x$, in canonical
form for $x \ne 0$ (for example $-2/3$ for $x = -3/2$) and the
non-canonical $1/0$ for $x = 0$, so that the convention of the remainder
holds trivially.

## Evaluating a continued fraction {#set}

`fmpq_set_cfrac(x, c, n)` sets $x = [c_0; c_1, \ldots, c_{n-1}]$, for $n
> 0$ and $c_1, \ldots, c_{n-1}$ positive (then the result is canonical).
The numerator and denominator are the first column of the matrix
product

$$
\begin{pmatrix} c_0 & 1 \\ 1 & 0 \end{pmatrix}
\begin{pmatrix} c_1 & 1 \\ 1 & 0 \end{pmatrix} \cdots
\begin{pmatrix} c_{n-1} & 1 \\ 1 & 0 \end{pmatrix},
$$

which the function evaluates by splitting it in halves recursively, with
the plain product below 32 factors. The balanced product has cost
$O(M(b) \log n)$ for a $b$-bit result, where multiplying the factors one
by one would cost $O(nb)$.

## Bounding the length {#bound}

`fmpq_cfrac_bound(x)` returns 1 for an integer and $\lfloor 1.4404\, b
\rfloor + 2$ otherwise, where $b$ is the bit length of the denominator.
It follows from the fact that a fraction whose shorter expansion has $k$
terms has denominator at least the Fibonacci number $F_{k+1}$.

## Example {#example}

The terms of $-103993/33102$, a convergent of $-\pi$, are computed in
two chunks, and the full expansion and a truncation are evaluated.

```c example name=fmpq_cfrac_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpq.h>

int main(void)
{
    fmpq_t x, r, y;
    fmpz * c;
    slong i, k, n, len;

    fmpq_init(x);
    fmpq_init(r);
    fmpq_init(y);

    fmpq_set_si(x, -103993, 33102);
    len = fmpq_cfrac_bound(x);
    c = _fmpz_vec_init(len);

    /* the first two terms and the remainder */
    k = fmpq_get_cfrac(c, r, x, 2);
    flint_printf("%wd terms, remainder ", k);
    fmpq_print(r);
    flint_printf("\n");

    /* continue with 1/r */
    fmpq_inv(r, r);
    n = k + fmpq_get_cfrac(c + k, r, r, len - k);
    flint_printf("terms:");
    for (i = 0; i < n; i++)
    {
        flint_printf(" ");
        fmpz_print(c + i);
    }
    flint_printf(" (bound %wd)\n", len);

    fmpq_set_cfrac(y, c, n);
    fmpq_print(y);
    flint_printf("\n");
    fmpq_set_cfrac(y, c, 4);
    fmpq_print(y);
    flint_printf("\n");

    _fmpz_vec_clear(c, len);
    fmpq_clear(x);
    fmpq_clear(r);
    fmpq_clear(y);
    return 0;
}
```
```output
2 terms, remainder 4687/28415
terms: -4 1 6 15 1 292 (bound 25)
-103993/33102
-333/106
```
