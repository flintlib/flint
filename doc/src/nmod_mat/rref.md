---
title: Computing echelon forms and nullspaces
functions: nmod_mat_rref nmod_mat_nullspace nmod_mat_left_nullspace
---

The reduced row echelon form of a matrix over $\Z/p\Z$ is the unique
matrix row equivalent to it in which each nonzero row starts with a 1
(its pivot), the pivots move strictly to the right from row to row, the
zero rows come last, and each pivot is the only nonzero entry in its
column. `nmod_mat_rref` computes it in place, and `nmod_mat_nullspace`
and `nmod_mat_left_nullspace` read bases of the nullspaces off it. The
modulus must be prime.

::: summary algorithms/linear-algebra/gaussian-elimination :::

## Algorithm {#algorithm}

`nmod_mat_rref` does not use Gauss–Jordan elimination. For an
$m \times n$ matrix $A$ it

1. computes an LU decomposition with `nmod_mat_lu_with_pivots`, which
   gives the rank $r$, the row echelon form $U$ (the first $r$ rows) and
   the pivot columns $c_0 < \cdots < c_{r-1}$;
2. clears the entries of $L$;
3. if $r < n$, writes $U = (U_1 \mid U_2)$ up to a column permutation,
   with $U_1$ the $r \times r$ upper triangular matrix of pivot columns
   and $U_2$ the other $n - r$ columns, and replaces $U_2$ by
   $U_1^{-1} U_2$ with `nmod_mat_solve_triu`;
4. sets the pivot columns to the columns of the identity matrix.

The zero rows of the result come last. The cost is that of the LU
decomposition plus one triangular solve, $O(m n r^{\omega - 2})$ field
operations with the exponent $\omega$ of matrix multiplication
[@JPS2013]. A matrix with one row is just divided by its first nonzero
entry.

## Nullspaces {#nullspace}

Let $R$ be the reduced row echelon form of $A$, with pivot columns
$c_0, \ldots, c_{r-1}$ and free columns $f_0 < \cdots < f_{d-1}$, where
$d = n - r$. For each free column $f_i$, the vector $x$ with
$x_{f_i} = 1$, $x_{f_j} = 0$ for $j \ne i$ and $x_{c_k} = -R_{k, f_i}$
satisfies $Rx = 0$, hence $Ax = 0$, and these $d$ vectors form a basis of
the right nullspace. `nmod_mat_nullspace(X, A)` writes them to the first
$d$ columns of $X$ in this order, sets the remaining entries of $X$ to
zero and returns $d$. $X$ must have $n$ rows and at least $d$ columns.
Since $d$ is not known in advance, an $n \times n$ matrix is the usual
choice. This basis is the unique one whose vectors have the entries
$x_{f_j} = \delta_{ij}$ in the free positions, so it depends only on $A$.

`nmod_mat_left_nullspace(X, A)` applies `nmod_mat_nullspace` to the
transpose of $A$ and writes the basis vectors of $\{v : vA = 0\}$ to the
first $d = m - r$ rows of $X$, which must have at least $d$ rows and
exactly $m$ columns.

## Example {#example}

A $3 \times 5$ matrix of rank 2 modulo 13:

```c example name=nmod_mat_rref_nullspace
#include <flint/flint.h>
#include <flint/nmod_mat.h>

int main(void)
{
    /* a 3 x 5 matrix of rank 2 modulo 13: row 2 = 2 row 0 + row 1 */
    ulong a[3][5] = { { 2, 4, 1, 3, 0 },
                      { 1, 2, 5, 0, 7 },
                      { 5, 10, 7, 6, 7 } };
    nmod_mat_t A, R, X, Y, T;
    slong i, j, rank, nullity;

    nmod_mat_init(A, 3, 5, 13);
    for (i = 0; i < 3; i++)
        for (j = 0; j < 5; j++)
            nmod_mat_entry(A, i, j) = a[i][j];

    nmod_mat_init_set(R, A);
    rank = nmod_mat_rref(R);
    flint_printf("rank %wd, reduced row echelon form:\n", rank);
    nmod_mat_print_pretty(R);

    /* right nullspace: the first 'nullity' columns of X */
    nmod_mat_init(X, 5, 5, 13);
    nullity = nmod_mat_nullspace(X, A);
    flint_printf("nullity %wd, basis in the columns of:\n", nullity);
    nmod_mat_print_pretty(X);

    nmod_mat_init(T, 3, 5, 13);
    nmod_mat_mul(T, A, X);
    flint_printf("A X is zero: %d\n", nmod_mat_is_zero(T));

    /* left nullspace: the first rows of Y */
    nmod_mat_init(Y, 3, 3, 13);
    nullity = nmod_mat_left_nullspace(Y, A);
    flint_printf("left nullity %wd, first row (%wu %wu %wu)\n", nullity,
        nmod_mat_entry(Y, 0, 0), nmod_mat_entry(Y, 0, 1), nmod_mat_entry(Y, 0, 2));
    nmod_mat_mul(T, Y, A);
    flint_printf("Y A is zero: %d\n", nmod_mat_is_zero(T));

    nmod_mat_clear(A);
    nmod_mat_clear(R);
    nmod_mat_clear(X);
    nmod_mat_clear(Y);
    nmod_mat_clear(T);
    return 0;
}
```
```output
rank 2, reduced row echelon form:
<3 x 5 integer matrix mod 13>
[ 1  2  0  6  5]
[ 0  0  1  4  3]
[ 0  0  0  0  0]
nullity 3, basis in the columns of:
<5 x 5 integer matrix mod 13>
[11  7  8  0  0]
[ 1  0  0  0  0]
[ 0  9 10  0  0]
[ 0  1  0  0  0]
[ 0  0  1  0  0]
A X is zero: 1
left nullity 1, first row (11 12 1)
Y A is zero: 1
```

The pivot columns are 0 and 2, the free columns 1, 3 and 4. The first
basis vector of the nullspace is $(-2, 1, 0, 0, 0)$, read off column 1 of
$R$. The left nullspace is spanned by $(-2, -1, 1)$, since row 2 of $A$ is
twice row 0 plus row 1.
