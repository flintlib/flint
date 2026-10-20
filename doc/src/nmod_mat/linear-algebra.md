---
title: Linear algebra modulo n
summary: How matrices are stored, how to work on blocks with windows, which functions need a prime modulus, and how to organise a computation so that it runs at the speed of matrix multiplication.
order: 1
---

## Storage and windows {#storage}

An `nmod_mat_t` holds a pointer to its first entry, its dimensions, a
row stride and the modulus. Entry $(i, j)$ is at offset
$i \cdot \mathit{stride} + j$, so each row is a contiguous array of
`ulong` that can be passed to the functions of [](nmod_vec), for example
`_nmod_vec_dot(nmod_mat_row_ptr(A, i), x, n, mod, params)`.

A *window*, created with `nmod_mat_window_init`, is an `nmod_mat_t` that
refers to a rectangular block of another matrix. It costs nothing to
create, needs no memory and can be used as input or output of any
function. Block algorithms, such as the recursive LU decomposition and
triangular solving of this module, are written with windows. The only
rules are that windows must be released with `nmod_mat_window_clear`
rather than `nmod_mat_clear`, and that the inputs and the output of a
function must either be the same matrix (where the function allows
aliasing) or not overlap at all. The first example below computes a Schur
complement in place with windows.

Entries must be kept in $[0, n)$. Functions that read entries assume this
and do not check it; for instance the multiplication kernels compute
bounds from $n - 1$, and an unreduced entry can make them return wrong
results without any warning. When filling a matrix from arbitrary
integers, reduce them with `nmod_set_ui`, `nmod_set_si` or
`fmpz_get_nmod`.

## Which functions need a prime modulus {#moduli}

The modulus $n$ may be any word $n \ge 1$, but many algorithms divide by
pivots and are only correct over a field. The table summarises the
situation; "prime" means that the result is only guaranteed for prime
$n$. None of these conditions is checked. With a composite modulus, the
functions marked "prime" raise a `FLINT_IMPINV` exception when they need
the inverse of a zero divisor, and may also return results that are
meaningless without any exception.

| Operation | Functions | Modulus |
|---|---|---|
| arithmetic, multiplication, powers, transposition | `nmod_mat_add`, `nmod_mat_mul`, `nmod_mat_pow`, ... | any |
| determinant | `nmod_mat_det`, `nmod_mat_det_howell` | any |
| characteristic polynomial | `nmod_mat_charpoly`, `nmod_mat_charpoly_berkowitz` | any |
| characteristic polynomial | `nmod_mat_charpoly_danilevsky` | any; returns 0 if it fails |
| Howell and strong echelon form | `nmod_mat_howell_form`, `nmod_mat_strong_echelon_form` | any |
| triangular solving | `nmod_mat_solve_tril`, `nmod_mat_solve_triu` | any, if the diagonal entries are invertible |
| LU decomposition, rank | `nmod_mat_lu`, `nmod_mat_rank` | prime |
| solving, inverse | `nmod_mat_solve`, `nmod_mat_can_solve`, `nmod_mat_inv` | prime |
| echelon forms, nullspaces | `nmod_mat_rref`, `nmod_mat_nullspace`, `nmod_mat_left_nullspace`, `nmod_mat_reduce_row` | prime |
| minimal polynomial | `nmod_mat_minpoly` | prime |

For a prime power or a general composite modulus, the Howell form plays
the role of the reduced row echelon form [@How1986; @StoMul1998]: two
matrices span the same submodule of $(\Z/n\Z)^c$ if and only if their
Howell forms agree. The second example below shows this. Alternatively,
one can work modulo the prime factors of $n$ and lift or combine the
results.

## Making computations fast {#performance}

Matrix multiplication modulo $n$ runs at a large fraction of the
arithmetic peak of the processor, through vector instructions, delayed
modular reduction, Strassen–Winograd and several threads; see
[](nmod_mat/mul). The other operations are fast in proportion to how
much of their work they hand to it.

- *Use matrices, not loops over vectors.* Solving $AX = B$ for 100
  right-hand sides at once costs far less than 100 calls with one
  right-hand side, because the triangular solves become matrix products.
  The same holds for products: `nmod_mat_mul` of a $k \times k$ matrix
  by a $k \times 100$ matrix is much faster than 100 matrix-vector
  products.
- *Prefer the high-level functions.* `nmod_mat_lu`, `nmod_mat_solve`,
  `nmod_mat_inv` and `nmod_mat_rref` switch to recursive algorithms built
  on matrix multiplication for large dimensions, and to elimination with
  delayed reduction for small ones.
- *Choose the modulus.* When the modulus can be chosen, as in
  multimodular algorithms, moduli up to $2^{26}$, up to $2^{32}$ and up
  to $2^{50}$ or $2^{52}$ are handled by successively more expensive
  vectorised kernels, and moduli above $2^{52}$ by a slower
  multimodular floating-point method or by the classical algorithm. The
  macro `NMOD_MAT_OPTIMAL_MODULUS_BITS` (`FLINT_BITS - 5`) predates the
  vectorised kernels and is no longer a good choice for them.
  Very small moduli $n \le 255$, in particular $n = 2$ and $n = 3$, have
  their own byte- and bit-packed kernels.
- *Use threads.* After `flint_set_num_threads(t)`, `nmod_mat_mul` and the
  functions built on it use up to $t$ threads for large matrices.

## Examples {#examples}

The Schur complement $S = D - C A^{-1} B$ of a $4 \times 4$ matrix
$M = \begin{pmatrix} A & B \\ C & D \end{pmatrix}$ with $2 \times 2$
blocks, computed in place with windows, and the identity
$\det M = \det A \det S$:

```c example name=nmod_mat_windows_schur
#include <flint/flint.h>
#include <flint/nmod_mat.h>

int main(void)
{
    nmod_mat_t M, A, B, C, D, S;
    slong i, j;

    /* M = [A B; C D] with 2 x 2 blocks, modulo 101 */
    nmod_mat_init(M, 4, 4, 101);
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            nmod_mat_entry(M, i, j) = (3 * i + j + (i == j)) % 101;

    nmod_mat_window_init(A, M, 0, 0, 2, 2);
    nmod_mat_window_init(B, M, 0, 2, 2, 4);
    nmod_mat_window_init(C, M, 2, 0, 4, 2);
    nmod_mat_window_init(D, M, 2, 2, 4, 4);

    /* Schur complement S = D - C A^-1 B, computed in place of D:
       first B <- A^-1 B, then D <- D - C B */
    nmod_mat_init(S, 2, 2, 101);
    if (nmod_mat_solve(B, A, B))
    {
        nmod_mat_submul(D, D, C, B);
        nmod_mat_set(S, D);
    }

    flint_printf("M after the update:\n");
    nmod_mat_print_pretty(M);
    flint_printf("det(M) = det(A) det(S): %wu\n",
                 nmod_mat_det(A) * nmod_mat_det(S) % 101);

    nmod_mat_window_clear(A);
    nmod_mat_window_clear(B);
    nmod_mat_window_clear(C);
    nmod_mat_window_clear(D);
    nmod_mat_clear(S);

    /* the determinant of the original matrix, for comparison */
    for (i = 0; i < 4; i++)
        for (j = 0; j < 4; j++)
            nmod_mat_entry(M, i, j) = (3 * i + j + (i == j)) % 101;
    flint_printf("det(M) directly:         %wu\n", nmod_mat_det(M));

    nmod_mat_clear(M);
    return 0;
}
```
```output
M after the update:
<4 x 4 integer matrix mod 101>
[  1   1  53  55]
[  3   5  50  49]
[  6   7  48  43]
[  9  10  44  38]
det(M) = det(A) det(S): 66
det(M) directly:         66
```

`nmod_mat_solve(B, A, B)` overwrites the window $B$ with $A^{-1} B$,
which the function allows, and `nmod_mat_submul(D, D, C, B)` updates the
window $D$ in place.

Two generating sets of the same submodule of $(\Z/12\Z)^2$ have the same
Howell form. Since $(6, 9) - (4, 6) = (2, 3)$, both modules are generated
by $(2, 3)$, and $6 \cdot (2, 3) = (0, 6)$ is the element that a row
echelon form would miss:

```c example name=nmod_mat_howell_form_mod12
#include <flint/flint.h>
#include <flint/nmod_mat.h>

int main(void)
{
    nmod_mat_t A, B;
    slong k;

    /* the rows (4, 6) and (6, 9) span a submodule of (Z/12Z)^2 */
    nmod_mat_init(A, 2, 2, 12);
    nmod_mat_entry(A, 0, 0) = 4; nmod_mat_entry(A, 0, 1) = 6;
    nmod_mat_entry(A, 1, 0) = 6; nmod_mat_entry(A, 1, 1) = 9;

    /* the same module with other generators: (2, 3) and (8, 0) */
    nmod_mat_init(B, 3, 2, 12);
    nmod_mat_entry(B, 0, 0) = 2; nmod_mat_entry(B, 0, 1) = 3;
    nmod_mat_entry(B, 1, 0) = 8; nmod_mat_entry(B, 1, 1) = 0;

    k = nmod_mat_howell_form(A);
    flint_printf("Howell form of A (%wd nonzero rows):\n", k);
    nmod_mat_print_pretty(A);

    k = nmod_mat_howell_form(B);
    flint_printf("Howell form of B (%wd nonzero rows):\n", k);
    nmod_mat_print_pretty(B);

    nmod_mat_clear(A);
    nmod_mat_clear(B);
    return 0;
}
```
```output
Howell form of A (2 nonzero rows):
<2 x 2 integer matrix mod 12>
[ 2  3]
[ 0  6]
Howell form of B (2 nonzero rows):
<3 x 2 integer matrix mod 12>
[ 2  3]
[ 0  6]
[ 0  0]
```

The Howell form functions need at least as many rows as columns, so a
module with fewer generators than the dimension must be padded with zero
rows.
