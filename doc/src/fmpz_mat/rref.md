---
title: Fraction-free elimination and echelon forms
functions: fmpz_mat_fflu fmpz_mat_rref
---

Gaussian elimination over $\Q$ applied to an integer matrix creates
fractions whose numerators and denominators grow quickly. Fraction-free
elimination avoids them: it multiplies rows by pivots instead of dividing
by them, and divides out the previous pivot exactly, so that every
intermediate entry is a minor of the input. `fmpz_mat_fflu` computes such
a decomposition, from which determinants, ranks and solutions follow;
`fmpz_mat_rref` computes the reduced row echelon form, either by
fraction-free elimination or by a multimodular method with an exact
certificate.

::: summary algorithms/linear-algebra/fraction-free :::

## Fraction-free LU decomposition {#fflu}

`fmpz_mat_fflu(B, den, perm, A, rank_check)` performs Bareiss' algorithm
[@Bar1968] on the $m \times n$ matrix $A$ and returns its rank $r$. Pivots
are chosen with `fmpz_mat_find_pivot_any`: in each column, the first row
at or below the current pivot row with a nonzero entry; rows are swapped
to bring it into place, and the swaps are applied to `perm` if it is not
`NULL`. When the $k$-th pivot $p_k$ has been chosen, every entry of the
remaining submatrix is updated as

$$
b_{ij} \gets \frac{p_k b_{ij} - b_{i c_k} b_{k j}}{p_{k-1}},
$$

with $p_{-1} = 1$, where $c_k$ is the pivot column; the division is
exact. After the elimination:

- the first $r$ rows of $B$, from their pivot columns on, form the
  fraction-free echelon form: the pivot in row $k$ is
  $p_k$, the determinant of the leading $(k+1) \times (k+1)$ submatrix of
  the row-permuted $A$ restricted to the first $k + 1$ pivot columns;
- below each pivot, $B$ keeps the entries of the pivot column as they were
  when the pivot was chosen (the multipliers of the fraction-free LU
  decomposition [@NakTurWil1997]); the other entries of the rows below the
  rank are zero;
- `den` is the last pivot $p_{r-1}$, or 1 if $r = 0$.

For a nonsingular square $A$, `den` is $\det(PA) = \pm \det(A)$, where
$P$ is the row permutation. With `rank_check` set, the function stops at
the first column without pivot and returns 0 with `den` = 0, which is
what the determinant and solving functions need. $B$ may be $A$.

The entries of $B$ are minors of $A$, so they have at most about
$n(b + \frac12 \log_2 n)$ bits for $b$-bit entries, and the cost is
$O(mnr)$ multiplications and exact divisions of such integers. For small
entries (at most 62 bits on a 64-bit machine) the updates use two-word
arithmetic and division by a precomputed inverse of the previous pivot
until an entry outgrows a word.

## Reduced row echelon form {#rref}

`fmpz_mat_rref(B, den, A)` sets ($B$, `den`) to the reduced row echelon
form of $A$: $B/\mathit{den}$ is in reduced row echelon form, so every
pivot entry of $B$ equals `den`, the other entries of the pivot columns
are zero, and the rows below the rank are zero. It returns the rank. The
denominator is not normalised: it is a nonzero $r \times r$ minor of $A$
for `fmpz_mat_rref_fflu`, and 1 or a divisor of such a minor for
`fmpz_mat_rref_mul`, so the two can return different but equivalent
results. Divide by the gcd of `den` and the entries for a canonical form,
or use `fmpq_mat_rref`.

For an $m \times n$ matrix, `fmpz_mat_rref` calls `fmpz_mat_rref_fflu`
if $m \le 3$, if $n \le 2$, if $m \le 20$ and $n > m$, or if
$20 < m \le 100$ and $n > m + \frac{m - 20}{80} m$; otherwise it calls
`fmpz_mat_rref_mul`. Wide
matrices thus go to fraction-free elimination, whose cost is dominated
by the elimination of the first $m$ columns, and tall or square ones to
the multimodular method.

### Fraction-free back substitution {#rref-fflu}

`fmpz_mat_rref_fflu` computes the fraction-free echelon form with
`fmpz_mat_fflu`. Writing its first $r$ rows, up to a permutation of the
columns, as $(U \mid V)$ with $U$ upper triangular, the reduced form is
$(I \mid U^{-1} V)$, and $\mathit{den} \cdot U^{-1} V$ is integral since
$\mathit{den} = \det U$ up to sign. Its columns are computed by back
substitution with exact divisions by the pivots, scaled to the common
denominator `den`; only the non-pivot columns need to be updated. This is
the fraction-free Gauss–Jordan elimination of [@NakTurWil1997]
restricted to the part $V$.

### Multimodular echelon form {#rref-mul}

`fmpz_mat_rref_mul` follows [@Stein2007]:

1. Reduce $A$ modulo a prime $p$ of 17 bits and compute its rank $r$,
   pivot rows and pivot columns modulo $p$ with
   `nmod_mat_lu_with_pivots`. If $r = n$, the reduced form is the identity
   matrix on top of zero rows, and $\mathit{den} = 1$.
2. Let $B_0$ be the $r \times r$ submatrix of $A$ formed by the pivot rows
   and columns, which is nonsingular modulo $p$ and hence over $\Q$, and
   $C$ the pivot rows restricted to the other columns. Solve
   $B_0 E_2 = \mathit{den} \cdot C$ with `fmpz_mat_solve`; the candidate is
   $\mathit{den}$ on the pivot positions and $E_2$ in the other columns.
3. Certify: check that the candidate is in reduced row echelon form, and
   that every non-pivot row of $A$ is annihilated by the nullspace basis
   defined by the candidate (one matrix product). Both checks are exact,
   so if they pass, the rows of $A$ lie in the span of the pivot rows,
   the rank is $r$ and the candidate is the reduced row echelon form
   (`fmpz_mat_rref_upper_certify_lu_mod_p`).
4. If a check fails, the rank modulo $p$ was below the rank of $A$; try
   the next prime.

The work is one solve with $n - r$ right-hand sides, done by the
output-sensitive modular solvers, and one product, instead of an
elimination with growing entries.

## Example {#example}

A fraction-free LU decomposition of a $4 \times 4$ matrix, whose pivots
are the leading principal minors $2, 5, 18, 82$, and the reduced row
echelon form of a $3 \times 5$ matrix of rank 2 computed in two ways:

```c example name=fmpz_mat_fflu_rref
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

int main(void)
{
    fmpz_mat_t A, B;
    fmpz_t den;
    slong i, j, rank, perm[4] = { 0, 1, 2, 3 };

    fmpz_init(den);
    fmpz_mat_init(A, 4, 4);
    fmpz_mat_init(B, 4, 4);

    /* tridiagonal with i + 2 on the diagonal, plus A[0][3] = 3 */
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            fmpz_set_si(fmpz_mat_entry(A, i, j),
                (i == j) ? i + 2 : (j == i + 1 || i == j + 1));
    fmpz_set_si(fmpz_mat_entry(A, 0, 3), 3);

    rank = fmpz_mat_fflu(B, den, perm, A, 0);
    flint_printf("rank %wd, den %{fmpz}\n", rank, den);
    fmpz_mat_print_pretty(B);
    flint_printf("\n");

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);

    /* row 2 = row 0 + 2 row 1 */
    fmpz_mat_init(A, 3, 5);
    fmpz_mat_init(B, 3, 5);
    for (j = 0; j < 5; j++)
    {
        fmpz_set_si(fmpz_mat_entry(A, 0, j), 2 * j + 2);
        fmpz_set_si(fmpz_mat_entry(A, 1, j), j * j - 3);
        fmpz_set_si(fmpz_mat_entry(A, 2, j), 2 * j + 2 + 2 * (j * j - 3));
    }

    rank = fmpz_mat_rref_fflu(B, den, A);
    flint_printf("rref_fflu: rank %wd, den %{fmpz}\n", rank, den);
    fmpz_mat_print_pretty(B);
    flint_printf("\n");

    rank = fmpz_mat_rref_mul(B, den, A);
    flint_printf("rref_mul: rank %wd, den %{fmpz}\n", rank, den);
    fmpz_mat_print_pretty(B);
    flint_printf("\n");

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_clear(den);
    return 0;
}
```
```output
rank 4, den 82
[[2 1 0 3]
[1 5 2 -3]
[0 2 18 8]
[0 0 5 82]
]
rref_fflu: rank 2, den 8
[[8 0 -16 -40 -72]
[0 8 20 36 56]
[0 0 0 0 0]
]
rref_mul: rank 2, den 8
[[8 0 -16 -40 -72]
[0 8 20 36 56]
[0 0 0 0 0]
]
```

Both echelon forms have the denominator 8, while the reduced fraction form
has denominator 2 since all entries are divisible by 4. In the
decomposition, the entries 1, 2 and 5 below the diagonal are the
multipliers.
