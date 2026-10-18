---
title: LU decompositions and linear systems
functions: gr_mat_lu gr_mat_fflu gr_mat_nonsingular_solve gr_mat_nonsingular_solve_den gr_mat_solve_field
---

Linear systems over a generic ring are solved by Gaussian elimination
in one of two forms. Over a field, `gr_mat_lu` computes an LU
decomposition with row pivoting, using inverses of the pivots. Over an
integral domain such as $\Z$ or $\Z[x]$, `gr_mat_fflu` computes a
fraction-free LU decomposition by Bareiss' algorithm, which divides
only exactly and keeps the entries as minors of the input. The solving
functions are built on these two decompositions and on triangular
solving.

::: summary algorithms/linear-algebra/gaussian-elimination :::

## LU decomposition {#lu}

`gr_mat_lu(rank, P, LU, A, rank_check, ctx)` computes $PA = LU$ for an
$m \times n$ matrix $A$ of rank $r$, where $P$ is a permutation of the
rows (the array `P` of length $m$ maps each row of $PA$ to its row in
$A$), $L$ is $m \times r$ lower triangular with ones on its diagonal
and $U$ is $r \times n$ in row echelon form. Both are stored in `LU`,
which has the shape of $A$: $U$ in the first $r$ rows on and above the
diagonal of pivots, $L$ below it without its unit diagonal, and zeros
elsewhere. For a nonsingular square matrix this is the usual
decomposition with $U$ upper triangular. `LU` may be the same matrix as
`A`.

The pivot in each column is the row found by
`gr_mat_find_nonzero_pivot`, a ring method: the first provably nonzero
entry by default, and the entry of largest absolute value for
floating-point rings. A column in which all candidates are provably
zero is skipped (it is not a pivot column); if some candidate cannot be
decided, the decomposition stops with `GR_UNABLE`, since the rank is
then unknown. The pivots are inverted with `gr_inv`; a pivot that is
nonzero but not invertible (in a ring that is not a field) stops the
decomposition with `GR_DOMAIN`, or `GR_UNABLE` if the ring pretends to
be a field.

If `rank_check` is nonzero, the function stops as soon as it has shown
that a square matrix is singular, returning `GR_SUCCESS` with `rank`
set to 0 and `LU` undefined; this saves work when only invertible
matrices are of interest. After any failure, `rank`, `P` and `LU` are
meaningless.

`gr_mat_lu_classical` eliminates one column at a time. `gr_mat_lu_recursive`
splits the columns in halves, decomposes the left half recursively,
updates the right half with a unit lower triangular solve and one
matrix product $A_{11} \gets A_{11} - A_{10} A_{01}$, and decomposes the
updated block, so that its cost is $O(n^\omega)$ with the exponent
$\omega$ of `gr_mat_mul` [@JPS2013]. The default method
`gr_mat_lu_generic` uses the classical algorithm when $A$ has fewer than
5 rows or columns and the recursive one otherwise; the recursive halves
are computed by `gr_mat_lu`, so a ring-specific method is used at every
level. `fmpz_mod` uses the same algorithms with cutoff 8, `mpn_mod` its
own decomposition with delayed reductions, and `nfloat` its own
floating-point algorithm.

## Fraction-free LU decomposition {#fflu}

`gr_mat_fflu(rank, P, LU, den, A, rank_check, ctx)` performs Bareiss'
elimination [@Bar1968]: at step $k$ with pivot $p_k$, every entry below
and to the right is replaced by $(p_k a_{ij} - a_{ik} a_{kj}) / p_{k-1}$,
where the division is exact by Sylvester's identity, since all
intermediate entries are minors of $A$
([](algorithms/linear-algebra/fraction-free)). The output has the same
packed format as for `gr_mat_lu`, but with the fraction-free factors,
and `den` is set to the last pivot, a minor of $A$ of order $r$ (for a
nonsingular square matrix, $\pm \det A$, the sign being that of the
permutation). The ring must be an integral domain with exact division
`gr_divexact`; over $\Z$ the entries stay bounded by Hadamard's bound
instead of growing exponentially as in the rational algorithm.

## Solving square systems {#nonsingular}

`gr_mat_nonsingular_solve(X, A, B, ctx)` solves $AX = B$ for a square $A$
and any number of right-hand sides (the columns of $B$). It uses
`gr_mat_nonsingular_solve_lu` (LU decomposition with `rank_check`, then
two triangular solves) over fields, rings pretending to be fields and
inexact rings, and `gr_mat_nonsingular_solve_fflu` (fraction-free LU,
then the solution times `den`, then an exact division by `den`) over
other rings. It returns `GR_DOMAIN` if $A$ is singular, even when the
system has solutions, and also, over a ring that is not a field, when
the solution has entries outside the ring. Over real and complex balls
this is Gaussian elimination in ball arithmetic; the preconditioned
solvers `arb_mat_solve` and `acb_mat_solve` usually give much tighter
enclosures for large systems
([](algorithms/ball-arithmetic/linear-systems)).

`gr_mat_nonsingular_solve_den(X, den, A, B, ctx)` solves the system over
the fraction field of an integral domain without leaving the ring: it
returns $X$ and $\mathit{den}$ with $AX = \mathit{den} \cdot B$, from the
fraction-free decomposition, so that $X / \mathit{den}$ is the solution.
The denominator is $\pm \det A$ and need not be minimal.

`gr_mat_inv` solves $AX = I$ with `gr_mat_nonsingular_solve`. To solve
many systems with the same matrix, compute the decomposition once and
call `gr_mat_nonsingular_solve_lu_precomp` or
`gr_mat_nonsingular_solve_fflu_precomp` for each right-hand side.

## General systems {#general}

`gr_mat_solve_field(X, A, B, ctx)` solves $AX = B$ for an $m \times n$
matrix $A$ of any rank over a field. It computes an LU decomposition,
solves with the pivot columns and sets the variables of the non-pivot
columns to zero, then checks that the remaining equations hold: if they
do not, the system is inconsistent and `GR_DOMAIN` is returned (with $X$
set to zero); if the check cannot be decided, `GR_UNABLE`. When the
system has several solutions, this particular one is returned; the
others differ from it by elements of the nullspace (`gr_mat_nullspace`).

## Example {#example}

The program solves $\begin{pmatrix} 2 & 1 \\ 1 & 3 \end{pmatrix} X =
\begin{pmatrix} 1 \\ 2 \end{pmatrix}$, whose solution $(1/5, 3/5)$ is not
integral, over $\Z$ and over $\mathbb{Q}$, and then a singular consistent and a
singular inconsistent system over $\mathbb{Q}$.

```c example name=gr_mat_solve_systems
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_mat.h>

/* print a column vector as (x0, x1, ...) */
static void
print_col(const gr_mat_t X, gr_ctx_t ctx)
{
    slong i;
    flint_printf("(");
    for (i = 0; i < gr_mat_nrows(X, ctx); i++)
    {
        gr_print(gr_mat_entry_srcptr(X, i, 0, ctx), ctx);
        if (i + 1 < gr_mat_nrows(X, ctx))
            flint_printf(", ");
    }
    flint_printf(")");
}

int main(void)
{
    gr_ctx_t ZZ, QQ;
    gr_mat_t A, B, X, S;
    gr_ptr den;
    int status;

    gr_ctx_init_fmpz(ZZ);
    gr_ctx_init_fmpq(QQ);

    /* A X = B over Z: the solution has denominators */
    gr_mat_init(A, 2, 2, ZZ);
    gr_mat_init(B, 2, 1, ZZ);
    gr_mat_init(X, 2, 1, ZZ);
    GR_MUST_SUCCEED(gr_mat_set_str(A, "[[2, 1], [1, 3]]", 0, ZZ));
    GR_MUST_SUCCEED(gr_mat_set_str(B, "[[1], [2]]", 0, ZZ));

    status = gr_mat_nonsingular_solve(X, A, B, ZZ);
    flint_printf("solve over Z:      status %d\n", status);

    den = gr_heap_init(ZZ);
    GR_MUST_SUCCEED(gr_mat_nonsingular_solve_den(X, den, A, B, ZZ));
    flint_printf("solve_den over Z:  X = "); print_col(X, ZZ);
    flint_printf(", den = "); gr_print(den, ZZ); flint_printf("\n");
    gr_heap_clear(den, ZZ);
    gr_mat_clear(X, ZZ);

    /* the same system over Q */
    gr_mat_init(S, 2, 1, QQ);
    {
        gr_mat_t AQ, BQ;
        gr_mat_init(AQ, 2, 2, QQ);
        gr_mat_init(BQ, 2, 1, QQ);
        GR_MUST_SUCCEED(gr_mat_set_gr_mat_other(AQ, A, ZZ, QQ));
        GR_MUST_SUCCEED(gr_mat_set_gr_mat_other(BQ, B, ZZ, QQ));
        GR_MUST_SUCCEED(gr_mat_nonsingular_solve(S, AQ, BQ, QQ));
        flint_printf("solve over Q:      X = "); print_col(S, QQ); flint_printf("\n");

        /* a singular but consistent system */
        GR_MUST_SUCCEED(gr_mat_set_str(AQ, "[[1, 2], [2, 4]]", 0, QQ));
        GR_MUST_SUCCEED(gr_mat_set_str(BQ, "[[3], [6]]", 0, QQ));
        status = gr_mat_nonsingular_solve(S, AQ, BQ, QQ);
        flint_printf("singular system:   nonsingular_solve status %d\n", status);
        status = gr_mat_solve_field(S, AQ, BQ, QQ);
        flint_printf("                   solve_field status %d, X = ", status);
        print_col(S, QQ); flint_printf("\n");

        /* an inconsistent system */
        GR_MUST_SUCCEED(gr_mat_set_str(BQ, "[[3], [7]]", 0, QQ));
        status = gr_mat_solve_field(S, AQ, BQ, QQ);
        flint_printf("inconsistent:      solve_field status %d\n", status);

        gr_mat_clear(AQ, QQ);
        gr_mat_clear(BQ, QQ);
    }

    gr_mat_clear(S, QQ);
    gr_mat_clear(A, ZZ);
    gr_mat_clear(B, ZZ);
    gr_ctx_clear(ZZ);
    gr_ctx_clear(QQ);
    return 0;
}
```
```output
solve over Z:      status 1
solve_den over Z:  X = (1, 3), den = 5
solve over Q:      X = (1/5, 3/5)
singular system:   nonsingular_solve status 1
                   solve_field status 0, X = (3, 0)
inconsistent:      solve_field status 1
```

Over $\Z$ the solution does not exist in the ring, so
`gr_mat_nonsingular_solve` returns `GR_DOMAIN` (1), while
`gr_mat_nonsingular_solve_den` returns the integral vector $(1, 3)$ with
the denominator $5 = \det A$.
