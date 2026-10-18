---
title: padic_mat.h
header: padic_mat.h
description: Matrices over the $p$-adic numbers, stored as $p^v M$ with an integer matrix $M$ reduced modulo $p^{N - v}$.
category: p-adic numbers
order: 3
---

This module implements dense matrices over $\mathbf{Q}_p$ in the
fixed-precision model of `padic.h`. A matrix is stored as $p^v M$ with
$v \in \Z$ and $M$ an `fmpz_mat_t`, together with one precision $N$
for all entries; it represents the class of $p^v M$ modulo
$p^N \mathbf{Z}_p$ entrywise. Every function computes its result
exactly from the stored inputs and reduces it to the precision of the
output. The module provides conversions, entry access, addition, scalar
multiplication and matrix multiplication; multiplication is done over
$\Z$ by `fmpz_mat_mul`. Linear algebra over $\mathbf{Q}_p$ with error
tracking (determinants, solving) is available through the generic rings,
as `gr_mat` over `padic_radix`.

::: summary algorithms/p-adic/overview :::

## Representation {#representation}

::: type
padic_mat_struct
padic_mat_t

A matrix $p^v M$ over $\mathbf{Q}_p$. The fields are `mat` (an
`fmpz_mat_struct`), `val` ($v$) and `N` (the precision). A
`padic_mat_t` is an array of one `padic_mat_struct`.
:::

The matrix is *canonical* if some entry of $M$ is not divisible by $p$,
or if $M = 0$ and $v = 0$; then $v$ is the minimum valuation of the
entries. It is *reduced* if it is canonical and all entries of $M$ lie
in $[0, p^{N - v})$, so that a reduced matrix with $v \ge N$ is zero.
Functions require canonical inputs and produce reduced outputs unless
stated otherwise; dimensions must be compatible, which is not checked.
One precision for all entries means that an entry of valuation $w > v$
is stored to relative precision $N - w$ only.

::: function
fmpz_mat_struct * padic_mat(const padic_mat_t A)

Returns a pointer to the integer matrix $M$, for use with the functions
of `fmpz_mat.h`.
:::

::: function
fmpz * padic_mat_entry(const padic_mat_t A, slong i, slong j)

Returns a pointer to the entry $(i, j)$ of $M$, which is not necessarily
a unit.
:::

::: macro
padic_mat_val(A)
padic_mat_prec(A)

The valuation $v$ and the precision $N$ of $A$, as lvalues.
:::

::: function
slong padic_mat_get_val(const padic_mat_t A)
slong padic_mat_get_prec(const padic_mat_t A)

Return the valuation $v$ and the precision $N$ of $A$.
:::

::: function
slong padic_mat_nrows(const padic_mat_t A)
slong padic_mat_ncols(const padic_mat_t A)

Return the number of rows and columns of $A$.
:::

## Memory management {#memory}

::: function
void padic_mat_init(padic_mat_t A, slong r, slong c)
void padic_mat_init2(padic_mat_t A, slong r, slong c, slong prec)

Initialise $A$ as the zero matrix with $r$ rows and $c$ columns and
precision `PADIC_DEFAULT_PREC` (that is, 20), respectively `prec`.
:::

::: function
void padic_mat_clear(padic_mat_t A)

Clears $A$.
:::

::: function
void _padic_mat_canonicalise(padic_mat_t A, const padic_ctx_t ctx)

Brings $A$ into canonical form by dividing $M$ by the largest power of
$p$ dividing all entries; sets $v = 0$ if $M = 0$. No reduction.
:::

::: function
void _padic_mat_reduce(padic_mat_t A, const padic_ctx_t ctx)

Reduces the canonical matrix $A$ modulo $p^N$, setting it to zero if
$v \ge N$.
:::

::: function
void padic_mat_reduce(padic_mat_t A, const padic_ctx_t ctx)

Canonicalises and reduces $A$.
:::

::: function
int padic_mat_is_empty(const padic_mat_t A)
int padic_mat_is_square(const padic_mat_t A)

Return whether $A$ has zero rows or zero columns, respectively as many
rows as columns.
:::

::: function
int padic_mat_is_canonical(const padic_mat_t A, const padic_ctx_t ctx)

Returns whether $A$ is in canonical form.
:::

::: function
int padic_mat_is_reduced(const padic_mat_t A, const padic_ctx_t ctx)

Returns whether $A$ is reduced: empty, or zero with $v = 0$, or
canonical with $v < N$ and all entries less than $p^{N - v}$. The
entries are not checked to be nonnegative.
:::

## Assignment and conversion {#conversion}

::: function
void padic_mat_set(padic_mat_t B, const padic_mat_t A, const padic_ctx_t ctx)

Sets $B$ to $A$, reduced to the precision of $B$.
:::

::: function
void padic_mat_swap(padic_mat_t A, padic_mat_t B)

Swaps $A$ and $B$, including dimensions, valuations and precisions, in
constant time.
:::

::: function
void padic_mat_swap_entrywise(padic_mat_t mat1, padic_mat_t mat2)

Swaps the entries of the integer matrices of `mat1` and `mat2`, which
must have the same dimensions, one by one. The valuations and
precisions are not swapped, so the result is meaningful only if they
agree.
:::

::: function
void padic_mat_zero(padic_mat_t A)

Sets $A$ to the zero matrix.
:::

::: function
void padic_mat_one(padic_mat_t A)

Sets the square matrix $A$ to the identity matrix, or to zero if its
precision is not positive.
:::

::: function
void padic_mat_set_fmpq_mat(padic_mat_t B, const fmpq_mat_t A, const padic_ctx_t ctx)

Sets $B$ to the rational matrix $A$, reduced to the precision of $B$.
The valuation of $B$ is the minimum valuation of the entries of $A$.
:::

::: function
void padic_mat_get_fmpq_mat(fmpq_mat_t B, const padic_mat_t A, const padic_ctx_t ctx)

Sets the rational matrix $B$ to $p^v M$.
:::

## Entries {#entries}

Because all entries share one valuation, setting an entry may rescale the
whole matrix, and an entry is not directly available as a `padic_t`.

::: function
void padic_mat_get_entry_padic(padic_t rop, const padic_mat_t op, slong i, slong j, const padic_ctx_t ctx)

Sets `rop` to the entry $(i, j)$, that is, to $p^v M_{ij}$ reduced to
the precision of `rop`.
:::

::: function
void padic_mat_set_entry_padic(padic_mat_t rop, slong i, slong j, const padic_t op, const padic_ctx_t ctx)

Sets the entry $(i, j)$ of `rop` to `op` and restores canonical form.
If `op` has smaller valuation than `rop`, all entries are multiplied by
a power of $p$. The value is not reduced to the precision of `rop`; call
`padic_mat_reduce` afterwards if `op` may have a higher precision than
`rop`.
:::

## Comparison {#comparison}

::: function
int padic_mat_equal(const padic_mat_t A, const padic_mat_t B)

Returns whether the stored values are equal, that is, whether the
valuations and the integer matrices agree. The matrices must be
canonical; the precisions are not compared.
:::

::: function
int padic_mat_is_zero(const padic_mat_t A)

Returns whether $A$ is the zero matrix (with $v = 0$).
:::

## Input and output {#io}

::: function
int padic_mat_fprint(FILE * file, const padic_mat_t A, const padic_ctx_t ctx)
int padic_mat_print(const padic_mat_t A, const padic_ctx_t ctx)

Print the number of rows, a space, the number of columns, two spaces and
the entries row by row, separated by spaces. In the mode `PADIC_TERSE`
the entries are printed as rationals, in `PADIC_VAL_UNIT` in the form
`p^w*u`; the mode `PADIC_SERIES` is not supported and throws. Return 1.
:::

::: function
int padic_mat_fprint_pretty(FILE * file, const padic_mat_t A, const padic_ctx_t ctx)
int padic_mat_print_pretty(const padic_mat_t A, const padic_ctx_t ctx)

Print $A$ as a list of rows in square brackets, one row per line, each
entry in the printing mode of the context. Return 1.
:::

## Random generation {#random}

::: function
void padic_mat_randtest(padic_mat_t A, flint_rand_t state, const padic_ctx_t ctx)

Sets $A$ to a random reduced matrix, with valuation drawn as for
`padic_randtest`.
:::

## Transpose {#transpose}

::: function
void padic_mat_transpose(padic_mat_t B, const padic_mat_t A)

Sets $B$ to the transpose of $A$. The precision of $B$ is not applied.
Aliasing is allowed for square matrices.
:::

## Addition and subtraction {#addition}

::: function
void _padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

Sets $C$ to $A + B$, reduced, for nonempty matrices with
$\operatorname{ord}_p(A) \ge \operatorname{ord}_p(B)$, which is not
checked.
:::

::: function
void padic_mat_add(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

Sets $C$ to $A + B$, reduced.
:::

::: function
void _padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

Sets $C$ to $A - B$ in canonical form, for nonempty matrices. The result
is not reduced in general.
:::

::: function
void padic_mat_sub(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

Sets $C$ to $A - B$, reduced.
:::

::: function
void _padic_mat_neg(padic_mat_t B, const padic_mat_t A)

Sets $B$ to $-A$ in canonical form; the entries are negated and not
reduced.
:::

::: function
void padic_mat_neg(padic_mat_t B, const padic_mat_t A, const padic_ctx_t ctx)

Sets $B$ to $-A$, reduced.
:::

## Scalar operations {#scalar}

::: function
void _padic_mat_scalar_mul_padic(padic_mat_t B, const padic_mat_t A, const padic_t c, const padic_ctx_t ctx)
void _padic_mat_scalar_mul_fmpz(padic_mat_t B, const padic_mat_t A, const fmpz_t c, const padic_ctx_t ctx)

Set $B$ to $cA$ in canonical form, without reduction. For an integer $c$
its power of $p$ is moved into the valuation.
:::

::: function
void padic_mat_scalar_mul_padic(padic_mat_t B, const padic_mat_t A, const padic_t c, const padic_ctx_t ctx)
void padic_mat_scalar_mul_fmpz(padic_mat_t B, const padic_mat_t A, const fmpz_t c, const padic_ctx_t ctx)

Set $B$ to $cA$, reduced.
:::

::: function
void padic_mat_scalar_div_fmpz(padic_mat_t B, const padic_mat_t A, const fmpz_t c, const padic_ctx_t ctx)

Sets $B$ to $A/c$, reduced. The power of $p$ in $c$ is subtracted from
the valuation and the rest of $c$ is inverted modulo $p^{N - v}$.
Throws if $c = 0$.
:::

## Multiplication {#multiplication}

::: function
void padic_mat_mul(padic_mat_t C, const padic_mat_t A, const padic_mat_t B, const padic_ctx_t ctx)

Sets $C$ to $AB$, reduced. The integer matrices are multiplied with
`fmpz_mat_mul` ([](algorithms/linear-algebra/integer-matrix-multiplication)),
the valuations are added, and the product is canonicalised and reduced:
unlike for polynomials, the product of two canonical matrices may be
divisible by $p$. Allows aliasing.
:::
