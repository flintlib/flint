---
title: fmpz_mod_vec.h
header: fmpz_mod_vec.h
description: Vectors of residues modulo $n$ stored as `fmpz` arrays, with elementwise arithmetic, scalar multiplication and dot products with a single reduction.
category: Integers mod n
order: 4
---

This module provides operations on vectors over $\Z/n\Z$ whose entries
are canonical residues in $[0, n)$ stored as plain `fmpz` arrays, with
the modulus given by an `fmpz_mod_ctx_t`. Vectors are allocated,
cleared, copied and compared with the `fmpz_vec` functions
(`_fmpz_vec_init`, `_fmpz_vec_clear`, `_fmpz_vec_set`,
`_fmpz_vec_equal`, ...). The functions are the building blocks of
`fmpz_mod_poly` and `fmpz_mod_mat`.

All functions take a vector as a pointer and a length. Unless stated
otherwise, inputs must be canonical and the output is canonical, and the
output may alias an input of the same length (the loops are elementwise).

## Conversion {#conversion}

::: function
void _fmpz_mod_vec_set_fmpz_vec(fmpz * A, const fmpz * B, slong len, const fmpz_mod_ctx_t ctx)

Sets the entries of $(A, \mathit{len})$ to the entries of
$(B, \mathit{len})$ reduced modulo $n$; $B$ may contain any integers and
may alias $A$.

::: details
Long vectors, or vectors with large entries or a large modulus, are reduced in parallel
when FLINT uses several threads.

:::
:::

::: function
void _fmpz_mod_vec_get_fmpz_vec_smod(fmpz * res, const fmpz * vec, slong len, const fmpz_mod_ctx_t ctx)

Given canonical residues $(\mathit{vec}, \mathit{len})$, sets
$(\mathit{res}, \mathit{len})$ to the symmetric residues in
$(-n/2, n/2]$. This is the same as `_fmpz_vec_scalar_smod_fmpz` but
faster, since the input is known to be reduced: each entry costs one
comparison and at most one subtraction.
:::

## Arithmetic {#arithmetic}

::: function
void _fmpz_mod_vec_neg(fmpz * A, const fmpz * B, slong len, const fmpz_mod_ctx_t ctx)

Sets $(A, \mathit{len})$ to $-(B, \mathit{len})$.
:::

::: function
void _fmpz_mod_vec_add(fmpz * a, const fmpz * b, const fmpz * c, slong n, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_vec_sub(fmpz * a, const fmpz * b, const fmpz * c, slong n, const fmpz_mod_ctx_t ctx)

Set $(a, n)$ to $(b, n) + (c, n)$ and $(b, n) - (c, n)$, elementwise.
Here $n$ is the length of the vectors, not the modulus.
:::

::: function
void _fmpz_mod_vec_mul(fmpz * A, const fmpz * B, const fmpz * C, slong len, const fmpz_mod_ctx_t ctx)

Sets $(A, \mathit{len})$ to the elementwise product of
$(B, \mathit{len})$ and $(C, \mathit{len})$.
:::

## Scalar multiplication {#scalar}

::: function
void _fmpz_mod_vec_scalar_mul_fmpz_mod(fmpz * A, const fmpz * B, slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx)

Sets $(A, \mathit{len})$ to $c \cdot (B, \mathit{len})$, for canonical
$c$. For $c = 1$ this is a copy.
:::

::: function
void _fmpz_mod_vec_scalar_addmul_fmpz_mod(fmpz * A, const fmpz * B, slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx)

Sets $(A, \mathit{len})$ to $(A, \mathit{len}) + c \cdot (B, \mathit{len})$,
for canonical $c$. Each entry is updated with one exact `fmpz_addmul`
and one reduction.
:::

::: function
void _fmpz_mod_vec_scalar_div_fmpz_mod(fmpz * A, const fmpz * B, slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx)

Sets $(A, \mathit{len})$ to $c^{-1} \cdot (B, \mathit{len})$. Throws an
exception (from `fmpz_mod_inv`) if $c$ is not invertible modulo $n$.
:::

## Dot products {#dot}

::: function
void _fmpz_mod_vec_dot(fmpz_t d, const fmpz * A, const fmpz * B, slong len, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_vec_dot_rev(fmpz_t r, const fmpz * a, const fmpz * b, slong len, const fmpz_mod_ctx_t ctx)

Set the output to $\sum_{i} A_i B_i$, respectively
$\sum_i a_i b_{\mathit{len} - 1 - i}$, modulo $n$.

::: details
The sum is computed exactly with `_fmpz_vec_dot` (respectively
`_fmpz_vec_dot_general`) and reduced once, so the entries may in fact be
any integers.

:::
:::
