---
title: Computing an LU decomposition
functions: nmod_mat_lu nmod_mat_lu_with_pivots
---

An LU decomposition with row pivoting writes an $m \times n$ matrix $A$
over a field $\Z/p\Z$ as a row permutation of a product $LU$, with $L$
lower triangular and $U$ in row echelon form. It is the basis of
`nmod_mat_solve`, `nmod_mat_inv`, `nmod_mat_det`, `nmod_mat_rank`,
`nmod_mat_rref` and `nmod_mat_nullspace`, and it is useful directly when
several systems with the same matrix must be solved. The functions
overwrite $A$ with both factors and return the rank.

::: summary algorithms/linear-algebra/gaussian-elimination :::

## Output {#output}

Let $r$ be the returned rank. On return:

- `P[i]` is the index of the row of the original $A$ that ended up in row
  $i$. `P` must have room for $m$ entries and need not be initialised.
- $U$ is the $r \times n$ matrix formed by the first $r$ rows of the
  output, keeping the entries in columns $j \ge i$ of row $i$ and
  replacing those in columns $j < i$ by zero. It is in row echelon form:
  the pivot (first nonzero entry) of row $i$ is in column $c_i$, with
  $c_0 < c_1 < \cdots < c_{r-1}$.
- $L$ is the $m \times r$ unit lower triangular matrix whose entry
  $(i, j)$ for $j < \min(i, r)$ is entry $(i, j)$ of the output; its
  diagonal entries are 1 and are not stored. The rows $i \ge r$ of the
  output contain $L$ in their first $r$ columns and zeros elsewhere.

Then row $i$ of $LU$ equals row `P[i]` of the original $A$. In matrix
terms, $A = \Pi L U$ where $\Pi$ is the permutation matrix with a 1 in
row `P[i]` and column $i$. `nmod_mat_permute_rows(A, P, NULL)` applies
the same permutation to another copy of $A$, as in the example below.

The pivots are chosen column by column: the pivot of column $c$ is the
first row, among those not yet used, with a nonzero entry in column $c$
after elimination. The output therefore does not depend on the algorithm,
and $c_0, \ldots, c_{r-1}$ is the column rank profile of $A$, the
lexicographically smallest set of $r$ linearly independent columns.
`nmod_mat_lu_with_pivots` returns these indices, followed by the indices
of the other columns.

## Checking the rank {#rank-check}

With a nonzero `rank_check`, the decomposition stops as soon as a column
without a pivot is found, and the function returns 0, leaving $A$ in an
undefined state. This is what `nmod_mat_solve` and `nmod_mat_inv` use to
reject singular matrices early. If $m \ge n$, the return value is then 0
exactly when the rank is less than $n$. If $m < n$ the result depends on
the algorithm: for example, `nmod_mat_lu_recursive` returns 0 for the
$2 \times 6$ matrix $(I_2 \mid 0)$, which has full rank 2. Use
`rank_check` only with $m \ge n$.

## Moduli {#moduli}

The algorithms divide by pivots, so the modulus should be prime. For a
composite modulus the functions still return a decomposition if every
pivot they meet happens to be invertible, but they raise a
`FLINT_IMPINV` exception when they meet a nonzero pivot that is not, and
the returned "rank" has no meaning. Use `nmod_mat_howell_form` or
`nmod_mat_det` for composite moduli.

## Algorithms {#algorithms}

`nmod_mat_lu_classical` is Gaussian elimination, with one reduction for
each update of an entry. `nmod_mat_lu_classical_delayed` performs the
same operations but accumulates the updates of each entry of the active
submatrix unreduced in one, two or three words, as in a dot product, and
reduces an entry only when its column becomes the pivot column or its row
the pivot row. `nmod_mat_lu_recursive` splits the columns into halves of
$\lfloor n/2 \rfloor$ and $\lceil n/2 \rceil$, decomposes the left half
recursively, updates the right half with a triangular solve
(`nmod_mat_solve_tril`) and a matrix product (`nmod_mat_submul`), and
decomposes the remaining block recursively. Its cost is dominated by
matrix multiplication, $O(m n r^{\omega - 2})$ operations with the
exponent $\omega$ of the multiplication algorithm used.

`nmod_mat_lu` uses the classical algorithm when $\min(m, n) \le 3$, or
$\le 7$ for $p \ge 2^{28}$. Otherwise it uses the delayed algorithm
when $\min(m, n)$ is below a crossover that depends on the bit size $b$
of $p$, and the recursive algorithm above it. The recursive algorithm
calls `nmod_mat_lu` for the halves, so the crossover also sets its base
case. The values were tuned on an AMD Zen 3 processor:

| $b$ | crossover |
|---|---|
| 1–2 | 64 |
| 3–7 | 244–332 |
| 8–11 | 792–800 |
| 12–24 | 400–424 |
| 25, 26, 27, 28 | 484, 1352, 1032, 260 |
| 29 | 68 |
| 30–32 | 56 |
| 33–60 | 112–168 |
| 61, 62 | 68, 64 |
| 63–64 | 104 |

The details, with proofs and complexity bounds, are in
[](algorithms/linear-algebra/gaussian-elimination).

## Example {#example}

The decomposition of a $4 \times 5$ matrix of rank 3 modulo 17, and a
check that $LU$ is $A$ with its rows permuted by `P`:

```c example name=nmod_mat_lu_unpack
#include <flint/flint.h>
#include <flint/nmod_mat.h>

int main(void)
{
    /* a 4 x 5 matrix of rank 3 modulo 17: row 3 = row 0 + row 1 */
    ulong a[4][5] = { { 0, 0, 3, 1, 4 },
                      { 2, 4, 1, 0, 5 },
                      { 1, 2, 7, 7, 1 },
                      { 2, 4, 4, 1, 9 } };
    nmod_mat_t A, LU, L, U, T;
    slong P[4], piv[5], i, j, r;

    nmod_mat_init(A, 4, 5, 17);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 5; j++)
            nmod_mat_entry(A, i, j) = a[i][j];

    nmod_mat_init_set(LU, A);
    r = nmod_mat_lu_with_pivots(P, piv, LU);

    flint_printf("rank %wd, P = [%wd %wd %wd %wd]\n", r, P[0], P[1], P[2], P[3]);
    flint_printf("pivot columns %wd %wd %wd, others %wd %wd\n",
                 piv[0], piv[1], piv[2], piv[3], piv[4]);
    nmod_mat_print_pretty(LU);

    /* unpack L (4 x r, unit diagonal) and U (r x 5) */
    nmod_mat_init(L, 4, r, 17);
    nmod_mat_init(U, r, 5, 17);
    for (i = 0; i < 4; i++)
        for (j = 0; j < r; j++)
            nmod_mat_entry(L, i, j) = (i == j) ? 1 :
                                      (j < i) ? nmod_mat_entry(LU, i, j) : 0;
    for (i = 0; i < r; i++)
        for (j = i; j < 5; j++)
            nmod_mat_entry(U, i, j) = nmod_mat_entry(LU, i, j);

    /* row i of L*U is row P[i] of A */
    nmod_mat_init(T, 4, 5, 17);
    nmod_mat_mul(T, L, U);
    nmod_mat_permute_rows(A, P, NULL);
    flint_printf("L*U == rows of A permuted by P: %d\n", nmod_mat_equal(T, A));

    nmod_mat_clear(A);
    nmod_mat_clear(LU);
    nmod_mat_clear(L);
    nmod_mat_clear(U);
    nmod_mat_clear(T);
    return 0;
}
```
```output
rank 3, P = [1 0 2 3]
pivot columns 0 2 3, others 1 4
<4 x 5 integer matrix mod 17>
[ 2  4  1  0  5]
[ 0  0  3  1  4]
[ 9  5  0  2  4]
[ 1  1  0  0  0]
L*U == rows of A permuted by P: 1
```

Row 0 of $A$ has a zero in column 0, so rows 0 and 1 are swapped. Column
1 has no pivot because it is twice column 0. In row 2 of the output, the
entries 9 and 5 belong to $L$ and the entries $0, 2, 4$ in columns 2 to 4
to $U$, whose pivot in this row is the 2 in column 3. Row 3 lies in the
span of the rows above, so it holds only its row of $L$, which says that
it is the sum of the first two rows of $U$.
