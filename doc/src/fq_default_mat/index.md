---
title: fq_default_mat.h
header: fq_default_mat.h
description: Dense matrices over any finite field $\mathbb{F}_q$, in the representation chosen by `fq_default`.
category: Finite fields
order: 2
---

This module implements dense matrices over a finite field given by an
`fq_default_ctx_t` (see [](fq_default)). A matrix is a union of the
matrix types of the five representations (`fq_zech_mat`, `fq_nmod_mat`,
`fq_mat`, `nmod_mat` and `fmpz_mod_mat`), and every function calls the
function of the representation of the context. The semantics are those
of [](fq_mat), where the algorithms are described: matrices are
initialised with fixed dimensions, outputs may alias inputs of the same
dimensions, and the linear algebra functions use Gaussian elimination
based on fast matrix multiplication. Over prime fields the `nmod_mat`
and `fmpz_mod_mat` functions are used, which have their own algorithms
(see [](nmod_mat)). The characteristic and minimal polynomials,
`fq_default_mat_charpoly` and `fq_default_mat_minpoly`, are declared in
`fq_default_poly.h` and documented with [](fq_default_poly).

::: summary algorithms/linear-algebra/gaussian-elimination :::

## Types {#types}

::: type
fq_default_mat_struct
fq_default_mat_t

A matrix over an `fq_default` field: a union of an `fq_mat_t`, an
`fq_nmod_mat_t`, an `fq_zech_mat_t`, an `nmod_mat_t` and an
`fmpz_mod_mat_t`, of which the member matching the context is used.
`fq_default_mat_t` is an array of one union.
:::

## Memory management {#memory}

::: function
void fq_default_mat_init(fq_default_mat_t mat, slong rows, slong cols, const fq_default_ctx_t ctx)

Initialises `mat` as a `rows` by `cols` zero matrix over the field.
:::

::: function
void fq_default_mat_init_set(fq_default_mat_t mat, const fq_default_mat_t src, const fq_default_ctx_t ctx)

Initialises `mat` as a copy of `src`.
:::

::: function
void fq_default_mat_clear(fq_default_mat_t mat, const fq_default_ctx_t ctx)

Frees the memory used by `mat`.
:::

::: function
void fq_default_mat_window_init(fq_default_mat_t window, const fq_default_mat_t mat, slong r1, slong c1, slong r2, slong c2, const fq_default_ctx_t ctx)
void fq_default_mat_window_clear(fq_default_mat_t window, const fq_default_ctx_t ctx)

Initialise `window` as the submatrix of rows $r_1 \le i < r_2$ and
columns $c_1 \le j < c_2$ of `mat`, sharing its entries, and clear such a
window without freeing the entries.
:::

## Entries and dimensions {#entries}

::: function
void fq_default_mat_entry(fq_default_t val, const fq_default_mat_t mat, slong i, slong j, const fq_default_ctx_t ctx)

Sets `val` to a copy of the entry in row $i$ and column $j$. (Unlike
`fq_mat_entry`, this does not return a pointer.)
:::

::: function
void fq_default_mat_entry_set(fq_default_mat_t mat, slong i, slong j, const fq_default_t x, const fq_default_ctx_t ctx)
void fq_default_mat_entry_set_fmpz(fq_default_mat_t mat, slong i, slong j, const fmpz_t x, const fq_default_ctx_t ctx)

Set the entry in row $i$ and column $j$ to $x$, respectively to the image
of the integer $x$ in $\mathbb{F}_p$.
:::

::: function
slong fq_default_mat_nrows(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
slong fq_default_mat_ncols(const fq_default_mat_t mat, const fq_default_ctx_t ctx)

Return the number of rows, respectively columns.
:::

## Assignment and permutations {#assignment}

::: function
void fq_default_mat_set(fq_default_mat_t mat1, const fq_default_mat_t mat2, const fq_default_ctx_t ctx)

Copies `mat2` into `mat1`, which must have the same dimensions.
:::

::: function
void fq_default_mat_swap(fq_default_mat_t mat1, fq_default_mat_t mat2, const fq_default_ctx_t ctx)

Swaps the two matrices efficiently; the dimensions may differ.
:::

::: function
void fq_default_mat_zero(fq_default_mat_t A, const fq_default_ctx_t ctx)
void fq_default_mat_one(fq_default_mat_t A, const fq_default_ctx_t ctx)

Set `A` to the zero matrix, respectively to the matrix with ones on the
main diagonal.
:::

::: function
void fq_default_mat_swap_rows(fq_default_mat_t mat, slong * perm, slong r, slong s, const fq_default_ctx_t ctx)
void fq_default_mat_swap_cols(fq_default_mat_t mat, slong * perm, slong r, slong s, const fq_default_ctx_t ctx)
void fq_default_mat_invert_rows(fq_default_mat_t mat, slong * perm, const fq_default_ctx_t ctx)
void fq_default_mat_invert_cols(fq_default_mat_t mat, slong * perm, const fq_default_ctx_t ctx)

Swap rows or columns $r$ and $s$, or reverse the order of the rows or
columns; if `perm` is not `NULL`, the same permutation is applied to it.
:::

## Conversion and concatenation {#conversion}

::: function
void fq_default_mat_set_nmod_mat(fq_default_mat_t mat1, const nmod_mat_t mat2, const fq_default_ctx_t ctx)
void fq_default_mat_set_fmpz_mod_mat(fq_default_mat_t mat1, const fmpz_mod_mat_t mat2, const fq_default_ctx_t ctx)
void fq_default_mat_set_fmpz_mat(fq_default_mat_t mat1, const fmpz_mat_t mat2, const fq_default_ctx_t ctx)

Set `mat1` to the image of the matrix `mat2` over $\mathbb{F}_p$, or over
$\Z$ with entries reduced modulo $p$. The modulus of `mat2` must be $p$
and the dimensions must agree.
:::

::: function
void fq_default_mat_concat_vertical(fq_default_mat_t res, const fq_default_mat_t mat1, const fq_default_mat_t mat2, const fq_default_ctx_t ctx)
void fq_default_mat_concat_horizontal(fq_default_mat_t res, const fq_default_mat_t mat1, const fq_default_mat_t mat2, const fq_default_ctx_t ctx)

Set `res` to the block matrix with `mat1` above, respectively to the left
of, `mat2`.
:::

## Comparison {#comparison}

::: function
int fq_default_mat_equal(const fq_default_mat_t mat1, const fq_default_mat_t mat2, const fq_default_ctx_t ctx)

Returns whether the matrices have the same dimensions and entries.
:::

::: function
int fq_default_mat_is_zero(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
int fq_default_mat_is_one(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
int fq_default_mat_is_empty(const fq_default_mat_t mat, const fq_default_ctx_t ctx)
int fq_default_mat_is_square(const fq_default_mat_t mat, const fq_default_ctx_t ctx)

Return whether the matrix is zero, the identity, empty (no rows or no
columns), respectively square.
:::

## Random matrices {#random}

::: function
void fq_default_mat_randtest(fq_default_mat_t mat, flint_rand_t state, const fq_default_ctx_t ctx)
void fq_default_mat_randrank(fq_default_mat_t mat, flint_rand_t state, slong rank, const fq_default_ctx_t ctx)
void fq_default_mat_randops(fq_default_mat_t mat, flint_rand_t state, slong count, const fq_default_ctx_t ctx)
void fq_default_mat_randtril(fq_default_mat_t mat, flint_rand_t state, int unit, const fq_default_ctx_t ctx)
void fq_default_mat_randtriu(fq_default_mat_t mat, flint_rand_t state, int unit, const fq_default_ctx_t ctx)

Random matrices for testing: random entries, a sparse matrix of the
given rank, `count` random elementary row and column operations (which
preserve the rank and the determinant), and random lower or upper
triangular matrices with unit (if `unit` is nonzero) or random nonzero
diagonal, as for `fq_mat`.
:::

## Input and output {#io}

::: function
int fq_default_mat_fprint(FILE * file, const fq_default_mat_t mat, const fq_default_ctx_t ctx)
int fq_default_mat_print(const fq_default_mat_t mat, const fq_default_ctx_t ctx)

Print the matrix to `file`, respectively `stdout`, in the plain format of
the representation (for `fq_mat`: the dimensions followed by the entries
in row-major order).
:::

::: function
int fq_default_mat_fprint_pretty(FILE * file, const fq_default_mat_t mat, const fq_default_ctx_t ctx)
int fq_default_mat_print_pretty(const fq_default_mat_t mat, const fq_default_ctx_t ctx)

Print the matrix as a bracketed list of bracketed rows to `file`,
respectively `stdout`.
:::

## Arithmetic {#arithmetic}

::: function
void fq_default_mat_transpose(fq_default_mat_t B, const fq_default_mat_t A, const fq_default_ctx_t ctx)

Sets $B$ to $A^T$.
:::

::: function
void fq_default_mat_add(fq_default_mat_t C, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)
void fq_default_mat_sub(fq_default_mat_t C, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)
void fq_default_mat_neg(fq_default_mat_t B, const fq_default_mat_t A, const fq_default_ctx_t ctx)

Set the output to $A + B$, $A - B$, respectively $-A$.
:::

::: function
void fq_default_mat_scalar_mul(fq_default_mat_t B, const fq_default_mat_t A, const fq_default_t c, const fq_default_ctx_t ctx)
void fq_default_mat_scalar_mul_fmpz(fq_default_mat_t B, const fq_default_mat_t A, const fmpz_t c, const fq_default_ctx_t ctx)
void fq_default_mat_scalar_mul_si(fq_default_mat_t B, const fq_default_mat_t A, const slong c, const fq_default_ctx_t ctx)
void fq_default_mat_scalar_mul_ui(fq_default_mat_t B, const fq_default_mat_t A, const ulong c, const fq_default_ctx_t ctx)

Set $B$ to $cA$.
:::

::: function
void fq_default_mat_mul(fq_default_mat_t C, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)

Sets $C = AB$, with Kronecker substitution for extension fields
(`fq_mat_mul` and its variants) and the fast multiplication of
`nmod_mat` or `fmpz_mod_mat` for prime fields.
:::

::: function
void fq_default_mat_submul(fq_default_mat_t D, const fq_default_mat_t C, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)

Sets $D = C - AB$.
:::

## Gaussian elimination {#elimination}

::: function
slong fq_default_mat_lu(slong * P, fq_default_mat_t A, int rank_check, const fq_default_ctx_t ctx)

Computes a generalised LU decomposition $PA = LU$ in place and returns
the rank, as `fq_mat_lu`; with nonzero `rank_check` it returns 0 as soon
as $A$ is found to be singular.
:::

::: function
slong fq_default_mat_rref(fq_default_mat_t B, const fq_default_mat_t A, const fq_default_ctx_t ctx)

Sets $B$ to the reduced row echelon form of $A$ and returns the rank.
:::

::: function
slong fq_default_mat_rank(const fq_default_mat_t A, const fq_default_ctx_t ctx)

Returns the rank of $A$.
:::

::: function
int fq_default_mat_inv(fq_default_mat_t B, fq_default_mat_t A, const fq_default_ctx_t ctx)

If the square matrix $A$ is invertible, sets $B = A^{-1}$ and returns 1;
otherwise returns 0.
:::

::: function
void fq_default_mat_solve_tril(fq_default_mat_t X, const fq_default_mat_t L, const fq_default_mat_t B, int unit, const fq_default_ctx_t ctx)
void fq_default_mat_solve_triu(fq_default_mat_t X, const fq_default_mat_t U, const fq_default_mat_t B, int unit, const fq_default_ctx_t ctx)

Set $X = L^{-1}B$, respectively $X = U^{-1}B$, for a lower, respectively
upper, triangular matrix with nonzero diagonal (assumed to be one if
`unit` is nonzero).
:::

::: function
int fq_default_mat_solve(fq_default_mat_t X, const fq_default_mat_t A, const fq_default_mat_t C, const fq_default_ctx_t ctx)

Solves $AX = C$ for a square matrix $A$; returns 0 if $A$ is singular.
:::

::: function
int fq_default_mat_can_solve(fq_default_mat_t X, const fq_default_mat_t A, const fq_default_mat_t B, const fq_default_ctx_t ctx)

Solves $AX = B$ for a matrix $A$ of any shape; returns 1 and sets $X$ to
a solution if one exists, and 0 otherwise.
:::

::: function
slong fq_default_mat_nullspace(fq_default_mat_t X, const fq_default_mat_t A, const fq_default_ctx_t ctx)

Sets the first columns of $X$ to a basis of the right nullspace of $A$
and returns its dimension; $X$ must be $n \times n$ for an $m \times n$
matrix $A$, and its other columns are set to zero.
:::

::: function
void fq_default_mat_similarity(fq_default_mat_t A, slong r, fq_default_t d, const fq_default_ctx_t ctx)

Replaces the square matrix $A$ by the similar matrix $P^{-1}AP$ described
for `fq_mat_similarity`.
:::

## Example {#example}

Solving a linear system over $\mathbb{F}_{2^{10}}$:

```c example name=fq_default_mat_solve
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fq_default_mat.h>

int main(void)
{
    fmpz_t p;
    fq_default_ctx_t K;
    fq_default_mat_t A, X, B, C;
    fq_default_t a;
    slong i, j;

    fmpz_init_set_ui(p, 2);
    fq_default_ctx_init(K, p, 10, "a");
    fq_default_mat_init(A, 3, 3, K);
    fq_default_mat_init(X, 3, 1, K);
    fq_default_mat_init(B, 3, 1, K);
    fq_default_mat_init(C, 3, 1, K);
    fq_default_init(a, K);

    /* A[i][j] = a^(i*j + i + j), B[i] = a^i */
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            fq_default_gen(a, K);
            fq_default_pow_ui(a, a, i * j + i + j, K);
            fq_default_mat_entry_set(A, i, j, a, K);
        }
        fq_default_gen(a, K);
        fq_default_pow_ui(a, a, i, K);
        fq_default_mat_entry_set(B, i, 0, a, K);
    }

    flint_printf("solvable: %d\n", fq_default_mat_solve(X, A, B, K));
    fq_default_mat_mul(C, A, X, K);
    flint_printf("A*X == B: %d\n", fq_default_mat_equal(C, B, K));
    flint_printf("rank: %wd\n", fq_default_mat_rank(A, K));

    fq_default_clear(a, K);
    fq_default_mat_clear(A, K);
    fq_default_mat_clear(X, K);
    fq_default_mat_clear(B, K);
    fq_default_mat_clear(C, K);
    fq_default_ctx_clear(K);
    fmpz_clear(p);
    return 0;
}
```
```output
solvable: 1
A*X == B: 1
rank: 3
```
