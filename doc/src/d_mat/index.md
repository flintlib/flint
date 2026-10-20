---
title: d_mat.h
header: d_mat.h
description: Dense matrices of double-precision floating-point numbers
category: Floating-point arithmetic
order: 8
---

This module provides dense matrices of `double` entries with basic
operations: memory management, comparison, transposition and classical
multiplication. Arithmetic uses the default rounding of the processor
(round to nearest). It is a small support module; for numerical linear
algebra at double or higher precision with more algorithms, use the
generic matrices `gr_mat_t` over an `nfloat` or `dfloat` context, or
`arb_mat_t` for rigorous results.

## Types and memory management {#memory}

::: type
d_mat_struct
d_mat_t

A matrix with `r` rows and `c` columns. The entries are stored in one
array `entries` in row-major order, and `rows` is an array of `r`
pointers to the starts of the rows. A `d_mat_t` is a length-one array
of `d_mat_struct`.
:::

::: function
void d_mat_init(d_mat_t mat, slong rows, slong cols)

Initialises `mat` as a `rows` by `cols` matrix with all entries zero.
Throws an exception if the number of entries overflows an `slong`.
:::

::: function
void d_mat_clear(d_mat_t mat)

Frees the memory used by `mat`.
:::

::: macro
d_mat_entry(mat, i, j)

The entry in row $i$ and column $j$ of `mat`, as an lvalue of type
`double`. The indices are not checked.
:::

::: function
slong d_mat_nrows(const d_mat_t mat)
slong d_mat_ncols(const d_mat_t mat)

Return the number of rows and columns of `mat`.
:::

## Assignment {#assignment}

::: function
void d_mat_set(d_mat_t mat1, const d_mat_t mat2)

Copies the entries of `mat2` to `mat1`, which must have the same
dimensions.
:::

::: function
void d_mat_swap_entrywise(d_mat_t mat1, d_mat_t mat2)

Swaps the entries of `mat1` and `mat2` one by one, rather than the
structs; the row pointers of both matrices stay valid. The matrices
must have the same dimensions.
:::

::: function
void d_mat_zero(d_mat_t mat)

Sets all entries of `mat` to zero.
:::

## Random generation and output {#random}

::: function
void d_mat_randtest(d_mat_t mat, flint_rand_t state, slong minexp, slong maxexp)

Sets the entries of `mat` to random values from `d_randtest_signed`:
signed numbers with exponents between `minexp` and `maxexp`, or zero.
:::

::: function
void d_mat_print(const d_mat_t mat)

Prints `mat` to standard output, one bracketed row per line, with the
entries in the `printf` format `%E`.
:::

## Comparison {#comparison}

::: function
int d_mat_equal(const d_mat_t mat1, const d_mat_t mat2)

Returns 1 if `mat1` and `mat2` have the same dimensions and equal
entries (as doubles), and 0 otherwise.
:::

::: function
int d_mat_approx_equal(const d_mat_t mat1, const d_mat_t mat2, double eps)

Returns 1 if `mat1` and `mat2` have the same dimensions and their
entries differ by at most `eps` in absolute value, and 0 otherwise.
:::

::: function
int d_mat_is_square(const d_mat_t mat)

Returns 1 if `mat` has as many rows as columns, and 0 otherwise.
:::

## Transpose and multiplication {#mul}

::: function
void d_mat_transpose(d_mat_t B, const d_mat_t A)

Sets $B$ to the transpose $A^T$. Throws an exception unless the
dimensions of $B$ are those of $A^T$. Aliasing is allowed (for square
matrices, by the dimension condition).
:::

::: function
void d_mat_mul_classical(d_mat_t C, const d_mat_t A, const d_mat_t B)

Sets $C = AB$ using the classical algorithm in double-precision
arithmetic. Throws an exception unless $C$ has as many rows as $A$ and
as many columns as $B$; the number of columns of $A$ must equal the
number of rows of $B$ (not checked). Aliasing of $C$ with $A$ or $B$
is allowed.

::: details
The function transposes $B$ and then forms the dot products of rows of
$A$ and rows of $B^T$ in blocks of 8 columns, so that the inner loop
runs over contiguous memory. Each entry is accumulated from partial
sums of at most 8 terms, so the result can differ in the last bits from
a straightforward loop. If $C$ is aliased with $A$ or $B$, the product
is computed into a temporary matrix and swapped into $C$, and the
dimensions of $C$ are not checked.
:::
:::
