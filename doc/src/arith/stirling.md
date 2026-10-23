---
title: Stirling numbers
functions: arith_stirling_number_1u arith_stirling_number_1u_vec arith_stirling_matrix_1u
---

The Stirling numbers of the first kind connect the falling and rising
factorials with the powers of $x$, and those of the second kind the
powers with the falling factorials:

$$
x(x-1)\cdots(x-n+1) = \sum_{k=0}^n s(n,k)\, x^k, \qquad
x^n = \sum_{k=0}^n S(n,k)\, x(x-1)\cdots(x-k+1).
$$

The unsigned numbers $|s(n,k)| = (-1)^{n-k} s(n,k)$ are the coefficients
of the rising factorial $x(x+1)\cdots(x+n-1)$ and count the permutations
of $n$ elements with $k$ cycles; $S(n,k)$ counts the partitions of an
$n$-set into $k$ blocks [@Comtet1974]. All of them are zero for $k > n$,
and $s(n,n) = S(n,n) = 1$, $s(n,0) = S(n,0) = 0$ for $n \ge 1$.

The module offers three interfaces, which use different algorithms:

- `arith_stirling_number_1u`, `arith_stirling_number_1` and
  `arith_stirling_number_2` compute one value;
- the `_vec` functions compute a row $k = 0, \ldots, \mathit{klen} - 1$
  for fixed $n$;
- the `matrix` functions fill a whole table of rows $0 \le n < r$.

To compute many values, use the row or matrix functions: they cost about
as much as a few isolated values.

::: summary algorithms/number-theory/stirling-bell :::

## Isolated values {#single}

The values $s(n,k)$ with $k$ near 0 or near $n$ are cheap: $|s(n,1)| =
(n-1)!$, $S(n,2) = 2^{n-1} - 1$, $S(n, n-1) = |s(n, n-1)| = \binom{n}{2}$.
In general, the functions choose by $n$ and the ratio $k/n$:

- $|s(n,k)|$ is a coefficient of the rising factorial, which
  `arith_stirling_number_1u` expands by binary splitting truncated to the
  needed length: as $(1+x)(2+x)\cdots(n-1+x)$ to length $k$ for $k < n/2$,
  and as $(1+x)(1+2x)\cdots(1+(n-1)x)$ to length $n-k+1$ otherwise, with
  only the wanted coefficient extracted from the last product. For $n >
  140$ and $k > 0.87 n$ it uses instead the exponential generating
  function $|s(n,k)| = (n!/k!)\, [x^{n-k}] \bigl(-\log(1-x)/x\bigr)^k$.
- $S(n,k)$ is computed by `arith_stirling_number_2` from the triangular
  recurrence $S(n,k) = k S(n-1,k) + S(n-1,k-1)$ in one or two words for
  $n \le 43$ (on 64-bit machines), and otherwise from one of
  $$
  S(n,k) = \frac{1}{k!} \sum_{j=0}^{k} (-1)^{k-j} \binom{k}{j} j^n
  \quad\text{or}\quad
  S(n,k) = \frac{n!}{k!} [x^{n-k}] \left(\frac{e^x - 1}{x}\right)^k :
  $$
  the power sum exactly for small $k/n$ (below about $0.95
  e^{-0.00022 n}$, and at most 0.8), the power series for $k/n$ above
  about $0.92 + 0.005 \log n$, and the power sum modulo word-size primes
  with Chinese remaindering in between.

The signed version `arith_stirling_number_1` is the unsigned one with the
sign $(-1)^{n-k}$.

## Rows {#rows}

A full row of $|s(n,k)|$ is the coefficient list of the rising factorial,
which `arith_stirling_number_1u_vec` computes by binary splitting with
fast polynomial multiplication, truncated to `klen` coefficients. The
row of $S(n,k)$ comes from the triangular recurrence for $n \le 80$; for
larger $n$, it is the convolution
$$
S(n,k) = \sum_{j=0}^{k} \frac{j^n}{j!} \cdot \frac{(-1)^{k-j}}{(k-j)!},
$$
a polynomial product computed exactly over $\Z$ (with a common
denominator) for $\mathit{klen} < n/2$ and modulo word-size primes,
followed by Chinese remaindering, otherwise. `klen` may exceed $n + 1$;
the extra entries are zero.

## Matrices {#matrices}

`arith_stirling_matrix_1u`, `arith_stirling_matrix_1` and
`arith_stirling_matrix_2` are wrappers of `gr_mat_stirling`, which fills
row $n$ from row $n - 1$ with the recurrences
$$
|s(n,k)| = (n-1) |s(n-1,k)| + |s(n-1,k-1)|, \qquad S(n,k) = k S(n-1,k) +
S(n-1,k-1).
$$
The matrix may have any shape; its entries above the diagonal are zero.
By the defining identities, the square matrices of $s$ and $S$ are
inverses of each other.

## Example {#example}

```c example name=arith_stirling_basic
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_mat.h>
#include <flint/arith.h>

int main(void)
{
    fmpz_t s;
    fmpz * row;
    fmpz_mat_t A, B, C;
    slong k;

    fmpz_init(s);
    arith_stirling_number_1(s, 10, 3);
    flint_printf("s(10,3) = ");
    fmpz_print(s);
    arith_stirling_number_2(s, 10, 3);
    flint_printf(", S(10,3) = ");
    fmpz_print(s);
    flint_printf("\n");

    row = _fmpz_vec_init(7);
    arith_stirling_number_2_vec(row, 6, 7);
    flint_printf("S(6,k):");
    for (k = 0; k < 7; k++)
    {
        flint_printf(" ");
        fmpz_print(row + k);
    }
    flint_printf("\n");

    arith_stirling_number_2(s, 2000, 1000);
    flint_printf("S(2000,1000) has %wu bits\n", fmpz_bits(s));

    fmpz_mat_init(A, 8, 8);
    fmpz_mat_init(B, 8, 8);
    fmpz_mat_init(C, 8, 8);
    arith_stirling_matrix_1(A);
    arith_stirling_matrix_2(B);
    fmpz_mat_mul(C, A, B);
    flint_printf("s S = I: %d\n", fmpz_mat_is_one(C));

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
    _fmpz_vec_clear(row, 7);
    fmpz_clear(s);
    return 0;
}
```
```output
s(10,3) = -1172700, S(10,3) = 9330
S(6,k): 0 1 31 90 65 15 1
S(2000,1000) has 11144 bits
s S = I: 1
```
