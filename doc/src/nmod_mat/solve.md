---
title: Solving systems of linear equations
functions: nmod_mat_solve nmod_mat_solve_vec nmod_mat_can_solve nmod_mat_can_solve_inner
---

These functions solve $AX = B$ over $\Z/p\Z$ for a prime $p$.
`nmod_mat_solve` and `nmod_mat_solve_vec` are for a square matrix $A$ and
report failure if $A$ is singular. `nmod_mat_can_solve` accepts any shape
and decides whether a solution exists. All of them compute an LU
decomposition of a copy of $A$ and finish with triangular solves, so
their cost is that of `nmod_mat_lu` plus two triangular solves with the
columns of $B$.

::: summary algorithms/linear-algebra/gaussian-elimination :::

## Square systems {#square}

`nmod_mat_solve(X, A, B)` decomposes a copy of $A$ with `nmod_mat_lu`
and `rank_check` set, so a singular $A$ is detected as soon as a column
without a pivot appears. It then permutes the rows of $B$ as in the
decomposition and solves $Ly = PB$ and $UX = y$ with
`nmod_mat_solve_tril` and `nmod_mat_solve_triu`. If $A$ is singular it
returns 0 without writing to $X$, even if the system has solutions; use
`nmod_mat_can_solve` for that case. Since $A$ and $B$ are copied before
$X$ is written, $X$ may be the same matrix as $A$ or $B$ when the
dimensions agree. If $A$ has no rows or $B$ has no columns, the
functions return 1.

`nmod_mat_solve_vec(x, A, b)` is the same for one right-hand side given
as an array, and allows `x == b`.

## Systems of any shape {#general}

`nmod_mat_can_solve(X, A, B)`, for an $m \times n$ matrix $A$, works as
follows.

1. Compute an LU decomposition $A = \Pi L U$ of rank $r$ without rank
   check, and the pivot columns $c_0 < \cdots < c_{r-1}$ of $U$.
2. Permute $B$ and solve $L_1 Y = (\Pi^{-1} B)_{1}$ for the first $r$
   rows, where $L_1$ is the top $r \times r$ block of $L$.
3. If $r < m$, check that the remaining $m - r$ rows satisfy
   $L_2 Y = (\Pi^{-1} B)_2$, with one matrix multiplication. If not, the
   system is inconsistent: set $X$ to zero and return 0.
4. Solve $U' Z = Y$, where $U'$ is the $r \times r$ upper triangular
   matrix formed by the pivot columns of $U$, and set row $c_i$ of $X$ to
   row $i$ of $Z$ and all other rows to zero.

The solution returned is thus the one whose entries for the free
(nonpivot) variables are zero. All solutions are obtained by adding
elements of the nullspace, see `nmod_mat_nullspace`.
`nmod_mat_can_solve_inner` also returns the rank, the permutation and the
pivot columns found in steps 1 and 2. The functions return 0 if the
dimensions of $X$, $A$ and $B$ do not match, and 1 with $X = 0$ if $A$ has
no rows or $B$ has no columns. $X$ may be aliased with $B$ but not with
$A$.

## Many right-hand sides {#reuse}

Putting all right-hand sides in the columns of one matrix $B$ is the most
efficient way to solve several systems with the same matrix: the
triangular solves then reduce to matrix multiplication. When the
right-hand sides become available one at a time, compute the LU
decomposition once and do the two triangular solves for each of them, as
in the second example below.

## Composite moduli {#composite}

For a composite modulus the functions may raise a `FLINT_IMPINV`
exception when a pivot is not invertible, and a return value of 0 does
not prove that no solution exists. Over $\Z/n\Z$, linear systems can be
solved with the Howell form [@StoMul1998], which `nmod_mat_howell_form`
computes; the module has no solver built on it.

## Examples {#examples}

A nonsingular system, then a singular $3 \times 4$ system with a
consistent and an inconsistent right-hand side:

```c example name=nmod_mat_solve_basic
#include <flint/flint.h>
#include <flint/nmod_mat.h>

int main(void)
{
    nmod_mat_t A, X, B, S;
    ulong x[3], b[3] = { 1, 2, 3 };
    slong rank, perm[3], pivots[3];
    int ok;

    /* a nonsingular system modulo 7 */
    nmod_mat_init(A, 3, 3, 7);
    nmod_mat_entry(A, 0, 0) = 2; nmod_mat_entry(A, 0, 1) = 1;
    nmod_mat_entry(A, 1, 1) = 3; nmod_mat_entry(A, 1, 2) = 1;
    nmod_mat_entry(A, 2, 0) = 1; nmod_mat_entry(A, 2, 2) = 4;

    ok = nmod_mat_solve_vec(x, A, b);
    flint_printf("solve_vec: %d, x = (%wu, %wu, %wu)\n", ok, x[0], x[1], x[2]);

    /* a singular 3 x 4 matrix: row 2 = row 0 + row 1 */
    nmod_mat_clear(A);
    nmod_mat_init(A, 3, 4, 7);
    nmod_mat_entry(A, 0, 0) = 1; nmod_mat_entry(A, 0, 1) = 2;
    nmod_mat_entry(A, 1, 2) = 1; nmod_mat_entry(A, 1, 3) = 5;
    nmod_mat_entry(A, 2, 0) = 1; nmod_mat_entry(A, 2, 1) = 2;
    nmod_mat_entry(A, 2, 2) = 1; nmod_mat_entry(A, 2, 3) = 5;

    nmod_mat_init(B, 3, 1, 7);
    nmod_mat_init(X, 4, 1, 7);
    nmod_mat_init(S, 3, 1, 7);

    /* consistent right-hand side: b2 = b0 + b1 */
    nmod_mat_entry(B, 0, 0) = 3;
    nmod_mat_entry(B, 1, 0) = 4;
    nmod_mat_entry(B, 2, 0) = 0;
    ok = nmod_mat_can_solve_inner(&rank, perm, pivots, X, A, B);
    nmod_mat_mul(S, A, X);
    flint_printf("can_solve: %d, rank %wd, pivots %wd %wd, X^T = (%wu %wu %wu %wu), AX == B: %d\n",
        ok, rank, pivots[0], pivots[1],
        nmod_mat_entry(X, 0, 0), nmod_mat_entry(X, 1, 0),
        nmod_mat_entry(X, 2, 0), nmod_mat_entry(X, 3, 0), nmod_mat_equal(S, B));

    /* inconsistent right-hand side */
    nmod_mat_entry(B, 2, 0) = 1;
    ok = nmod_mat_can_solve(X, A, B);
    flint_printf("can_solve: %d, X is zero: %d\n", ok, nmod_mat_is_zero(X));

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    nmod_mat_clear(X);
    nmod_mat_clear(S);
    return 0;
}
```
```output
solve_vec: 1, x = (0, 1, 6)
can_solve: 1, rank 2, pivots 0 2, X^T = (3 0 4 0), AX == B: 1
can_solve: 0, X is zero: 1
```

The free variables $x_1$ and $x_3$ (columns 1 and 3 have no pivot) are
set to zero. One LU decomposition reused for three right-hand sides:

```c example name=nmod_mat_solve_reuse_lu
#include <flint/flint.h>
#include <flint/nmod.h>
#include <flint/nmod_mat.h>

int main(void)
{
    slong n = 50, i, j, t, rank;
    slong * P = flint_malloc(n * sizeof(slong));
    nmod_mat_t A, LU, b, y, x, check;
    ulong p = 1000003;
    int all_ok = 1;

    nmod_mat_init(A, n, n, p);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            nmod_mat_entry(A, i, j) = nmod_pow_ui(i + 2, j + 1, A->mod);

    /* a Vandermonde-like matrix, nonsingular; factor it once */
    nmod_mat_init_set(LU, A);
    rank = nmod_mat_lu(P, LU, 1);
    flint_printf("rank = %wd\n", rank);

    nmod_mat_init(b, n, 1, p);
    nmod_mat_init(y, n, 1, p);
    nmod_mat_init(x, n, 1, p);
    nmod_mat_init(check, n, 1, p);

    /* then solve A x = b for several right-hand sides b */
    for (t = 0; t < 3 && rank == n; t++)
    {
        for (i = 0; i < n; i++)
            nmod_mat_entry(b, i, 0) = (t + 1) * (i + 1) % p;

        /* permute b, solve L y = b with unit diagonal, then U x = y */
        for (i = 0; i < n; i++)
            nmod_mat_entry(y, i, 0) = nmod_mat_entry(b, P[i], 0);
        nmod_mat_solve_tril(y, LU, y, 1);
        nmod_mat_solve_triu(x, LU, y, 0);

        nmod_mat_mul(check, A, x);
        all_ok &= nmod_mat_equal(check, b);
        flint_printf("b%wd: x[0] = %wu\n", t, nmod_mat_entry(x, 0, 0));
    }
    flint_printf("A x == b in all cases: %d\n", all_ok);

    nmod_mat_clear(A);
    nmod_mat_clear(LU);
    nmod_mat_clear(b);
    nmod_mat_clear(y);
    nmod_mat_clear(x);
    nmod_mat_clear(check);
    flint_free(P);
    return 0;
}
```
```output
rank = 50
b0: x[0] = 862864
b1: x[0] = 725725
b2: x[0] = 588586
A x == b in all cases: 1
```

The upper triangular solve reads only the upper triangular part of `LU`, and the
lower triangular solve with `unit` set reads only the entries below the
diagonal, so both use the packed output of `nmod_mat_lu` directly.
