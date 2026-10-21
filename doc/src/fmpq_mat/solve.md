---
title: Solving rational linear systems
functions: fmpq_mat_solve fmpq_mat_can_solve fmpq_mat_inv
---

The solvers of this module find exact solutions of $AX = B$ over $\Q$
for rational or integer matrices $A$ and $B$, returning the solution as an
`fmpq_mat_t` in canonical form. All of them reduce the problem to integer
matrices and call the solvers of [](fmpz_mat), which work modulo primes
and reconstruct the rational solution; see [](fmpz_mat/solve).

::: summary algorithms/linear-algebra/multimodular :::

## Reduction to integer systems {#clearing}

Multiplying row $i$ of the augmented matrix $(A \mid B)$ by the least
common denominator $d_i$ of its entries gives an integer system
$A'X = B'$ with the same solutions, since
$A' = \operatorname{diag}(d) A$ and $B' = \operatorname{diag}(d) B$ with
$\operatorname{diag}(d)$ invertible. `fmpq_mat_get_fmpz_mat_rowwise_2`
computes $A'$ and $B'$. The integer solver returns $X$ as an integer
matrix and a denominator, which `fmpq_mat_set_fmpz_mat_div_fmpz` turns
into a rational matrix in canonical form. The variants whose names contain
`fmpz_mat` take $A$ and $B$ as integer matrices directly.

Clearing denominators row by row keeps the integer entries as small as
possible: for a matrix with a common structure of denominators, such as
the Hilbert matrix, the integers $d_i$ are the least common multiples of
the denominators in each row, not of the whole matrix.

## Choice of algorithm {#dispatch}

For a square $n \times n$ matrix $A$ and $k$ right-hand sides,
`fmpq_mat_solve` and `fmpq_mat_solve_fmpz_mat` use

- for $n \le 15$: `fmpz_mat_solve`, that is Cramer's rule for $n \le 3$
  and fraction-free elimination otherwise (the `fraction_free` variants);
- for $n \ge 16$ and $2k \ge n$: multimodular solving
  (`fmpz_mat_solve_multi_mod_den`, the `multi_mod` variants);
- otherwise: $p$-adic lifting (`fmpz_mat_solve_dixon_den`, the `dixon`
  variants).

All of them return 1 if $A$ is nonsingular and 0 otherwise, and in the
latter case leave $X$ unchanged. The modular solvers are output-sensitive:
their running time depends on the size of the solution rather than on a
worst-case bound.

`fmpq_mat_can_solve` accepts any $m \times n$ matrix $A$. It returns 1
with a solution if the system is consistent, and 0 if it is not; the
answer is proven, by exact verification of the solution or by an exact
certificate of inconsistency. For $m \le 15$ it uses fraction-free
elimination (`fmpq_mat_can_solve_fraction_free`), otherwise the modular
method of `fmpz_mat_can_solve`; `fmpq_mat_can_solve_dixon` and
`fmpq_mat_can_solve_multi_mod` force the solver used for the nonsingular
pivot block.

`fmpq_mat_inv` solves $A' X = \operatorname{diag}(d)$ with
`fmpq_mat_solve_fmpz_mat`, where $A' = \operatorname{diag}(d) A$, so that
$X = A^{-1}$; it uses explicit formulas for dimensions up to 2.

## Example {#example}

The $10 \times 10$ Hilbert matrix $H$, with entries $1/(i + j + 1)$, is
notoriously ill-conditioned but has an integral inverse. Solving
$Hx = (1, \ldots, 1)^T$ exactly, inverting $H$, and testing two
singular systems:

```c example name=fmpq_mat_solve_hilbert
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

int main(void)
{
    fmpq_mat_t H, B, X, Hinv, A, Y;
    slong i, n = 10;

    fmpq_mat_init(H, n, n);
    fmpq_mat_init(B, n, 1);
    fmpq_mat_init(X, n, 1);
    fmpq_mat_init(Hinv, n, n);

    fmpq_mat_hilbert_matrix(H);
    for (i = 0; i < n; i++)
        fmpq_one(fmpq_mat_entry(B, i, 0));

    flint_printf("solvable: %d\n", fmpq_mat_solve(X, H, B));
    flint_printf("x_0 = %{fmpq}, x_9 = %{fmpq}\n",
                 fmpq_mat_entry(X, 0, 0), fmpq_mat_entry(X, 9, 0));

    flint_printf("invertible: %d\n", fmpq_mat_inv(Hinv, H));
    flint_printf("inverse is integral: %d\n", fmpq_mat_is_integral(Hinv));
    flint_printf("entry (4, 4) of the inverse: %{fmpq}\n", fmpq_mat_entry(Hinv, 4, 4));

    /* a singular system: rows 0, 1 and their sum */
    fmpq_mat_init(A, 3, 2);
    fmpq_mat_init(Y, 2, 1);
    fmpq_mat_clear(B);
    fmpq_mat_init(B, 3, 1);
    fmpq_set_si(fmpq_mat_entry(A, 0, 0), 1, 2);
    fmpq_set_si(fmpq_mat_entry(A, 0, 1), 1, 3);
    fmpq_set_si(fmpq_mat_entry(A, 1, 0), 1, 4);
    fmpq_set_si(fmpq_mat_entry(A, 1, 1), -1, 5);
    fmpq_add(fmpq_mat_entry(A, 2, 0), fmpq_mat_entry(A, 0, 0), fmpq_mat_entry(A, 1, 0));
    fmpq_add(fmpq_mat_entry(A, 2, 1), fmpq_mat_entry(A, 0, 1), fmpq_mat_entry(A, 1, 1));
    fmpq_set_si(fmpq_mat_entry(B, 0, 0), 1, 1);
    fmpq_set_si(fmpq_mat_entry(B, 1, 0), 2, 1);
    fmpq_set_si(fmpq_mat_entry(B, 2, 0), 3, 1);

    flint_printf("consistent: %d, ", fmpq_mat_can_solve(Y, A, B));
    flint_printf("y = (%{fmpq}, %{fmpq})\n",
                 fmpq_mat_entry(Y, 0, 0), fmpq_mat_entry(Y, 1, 0));

    fmpq_set_si(fmpq_mat_entry(B, 2, 0), 4, 1);
    flint_printf("consistent: %d\n", fmpq_mat_can_solve(Y, A, B));

    fmpq_mat_clear(H);
    fmpq_mat_clear(B);
    fmpq_mat_clear(X);
    fmpq_mat_clear(Hinv);
    fmpq_mat_clear(A);
    fmpq_mat_clear(Y);
    return 0;
}
```
```output
solvable: 1
x_0 = -10, x_9 = 923780
invertible: 1
inverse is integral: 1
entry (4, 4) of the inverse: 176752976400
consistent: 1, y = (52 / 11, -45 / 11)
consistent: 0
```
