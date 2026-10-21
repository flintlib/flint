---
title: Determinants of polynomial matrices
functions: fmpz_poly_mat_det fmpz_poly_mat_rank
---

`fmpz_poly_mat_det` computes the determinant of a square matrix with
entries in $\Z[x]$, a polynomial of degree at most $n(\ell - 1)$ for an
$n \times n$ matrix whose entries have length at most $\ell$.
`fmpz_poly_mat_rank` computes the rank
over $\Q(x)$ by fraction-free elimination.

::: summary algorithms/linear-algebra/polynomial-matrices :::

## Methods {#methods}

| Dimension | Method |
|---|---|
| $n \le 2$ | explicit formula |
| $3 \le n \le 14$ | `fmpz_poly_mat_det_fflu` |
| $n \ge 15$ | `fmpz_poly_mat_det_interpolate` |

`fmpz_poly_mat_det_fflu` runs fraction-free elimination
(`fmpz_poly_mat_fflu`) with rank check on a copy of the matrix; the last
pivot is the determinant up to the sign of the row permutation. All
intermediate entries are minors of the matrix, so their degrees are at
most $n(\ell - 1)$, and pivots of smallest length and coefficient size
are preferred ([](algorithms/linear-algebra/fraction-free)).

`fmpz_poly_mat_det_interpolate` evaluates the matrix at the
$L = n(\ell - 1) + 1$ integers $-\lfloor L/2 \rfloor, \ldots,
L - 1 - \lfloor L/2 \rfloor$, which are small so that the integer
matrices stay small; computes the $L$ determinants with `fmpz_mat_det`,
which uses multimodular methods for large matrices
([](algorithms/linear-algebra/multimodular#det)); and interpolates with
`fmpz_poly_interpolate_exact`. Its cost grows like $L$ times that of an
integer determinant, while fraction-free elimination performs $O(n^3)$
operations on polynomials of growing degree and coefficient size.

The rank is the number of pivots found by `fmpz_poly_mat_fflu`. A
determinant can also be computed modulo a prime with
`nmod_poly_mat_det`, and an integer determinant from the evaluation at
a point with `fmpz_poly_mat_evaluate_fmpz`.

## Example {#example}

The determinant of $xI - A$ for an integer matrix $A$ is its
characteristic polynomial, which can be compared with
`fmpz_mat_charpoly`; both determinant methods agree, and replacing a row
by a combination of the others makes the matrix singular:

```c example name=fmpz_poly_mat_det_charpoly
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly_mat.h>

int main(void)
{
    fmpz_mat_t A;
    fmpz_poly_mat_t M;
    fmpz_poly_t d1, d2, c;
    slong i, j, n = 6;

    fmpz_mat_init(A, n, n);
    fmpz_poly_mat_init(M, n, n);
    fmpz_poly_init(d1);
    fmpz_poly_init(d2);
    fmpz_poly_init(c);

    /* M = x I - A with A[i][j] = (i * j + 1) mod 5 - 2 */
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpz_set_si(fmpz_mat_entry(A, i, j), (i * j + 1) % 5 - 2);
            fmpz_poly_set_si(fmpz_poly_mat_entry(M, i, j),
                             -((i * j + 1) % 5 - 2));
        }
        fmpz_poly_set_coeff_si(fmpz_poly_mat_entry(M, i, i), 1, 1);
    }

    fmpz_poly_mat_det_fflu(d1, M);
    fmpz_poly_mat_det_interpolate(d2, M);
    fmpz_mat_charpoly(c, A);

    flint_printf("det = %{fmpz_poly}\n", d1);
    flint_printf("interpolation agrees: %d\n", fmpz_poly_equal(d1, d2));
    flint_printf("equals charpoly: %d\n", fmpz_poly_equal(d1, c));
    flint_printf("rank: %wd\n", fmpz_poly_mat_rank(M));

    /* last row := x * (row 0) + row 1, a dependent row */
    for (j = 0; j < n; j++)
    {
        fmpz_poly_shift_left(fmpz_poly_mat_entry(M, n - 1, j),
                             fmpz_poly_mat_entry(M, 0, j), 1);
        fmpz_poly_add(fmpz_poly_mat_entry(M, n - 1, j),
                      fmpz_poly_mat_entry(M, n - 1, j),
                      fmpz_poly_mat_entry(M, 1, j));
    }
    fmpz_poly_mat_det(d1, M);
    flint_printf("singular: det = %{fmpz_poly}, rank %wd\n",
                 d1, fmpz_poly_mat_rank(M));

    fmpz_mat_clear(A);
    fmpz_poly_mat_clear(M);
    fmpz_poly_clear(d1);
    fmpz_poly_clear(d2);
    fmpz_poly_clear(c);
    return 0;
}
```
```output
det = x^6 + 6 * x^5 - 10 * x^4 - 80 * x^3 + 25 * x^2 + 250 * x
interpolation agrees: 1
equals charpoly: 1
rank: 6
singular: det = 0, rank 5
```
