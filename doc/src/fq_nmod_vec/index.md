---
title: fq_nmod_vec.h
header: fq_nmod_vec.h
template-headers: fq_vec_templates
description: Vectors over finite fields $\mathbb{F}_q$ in the word-size characteristic.
category: Finite fields
order: 12
---

This module provides low-level functions on vectors over a finite field $\mathbb{F}_q$ in the representation of [](fq_nmod) (a prime $p < 2^{\bits}$, elements stored as `nmod_poly`): plain C arrays of `fq_nmod_struct` with their length. It is generated from the same template as [](fq_vec) and has the same functions, with `fq_nmod` in place of `fq`; each function links to the `fq_vec` version, which has the details.

`_fq_nmod_vec_dot` (a wrapper of `_gr_vec_dot`) accumulates the products of the representatives modulo $p$ and reduces modulo the defining polynomial once.

## Memory management {#memory}

::: function
fq_nmod_struct * _fq_nmod_vec_init(slong len, const fq_nmod_ctx_t ctx)

Returns a newly allocated vector of `len` elements, all initialised to
zero.

::: details
Generated from the same template as `_fq_vec_init`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_vec_clear(fq_nmod_struct * vec, slong len, const fq_nmod_ctx_t ctx)

Clears the `len` elements of `vec` and frees the array.

::: details
Generated from the same template as `_fq_vec_clear`, whose documentation has the details.
:::
:::

## Assignment and basic manipulation {#assignment}

::: function
void _fq_nmod_vec_set(fq_nmod_struct * v, const fq_nmod_struct * f, slong len, const fq_nmod_ctx_t ctx)

Copies `(f, len)` into `(v, len)`.

::: details
Generated from the same template as `_fq_vec_set`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_vec_swap(fq_nmod_struct * vec1, fq_nmod_struct * vec2, slong len2, const fq_nmod_ctx_t ctx)

Swaps the elements of `(vec1, len2)` and `(vec2, len2)` efficiently.

::: details
Generated from the same template as `_fq_vec_swap`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_vec_zero(fq_nmod_struct * v, slong len, const fq_nmod_ctx_t ctx)

Sets the `len` elements of `v` to zero.

::: details
Generated from the same template as `_fq_vec_zero`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_vec_neg(fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const fq_nmod_ctx_t ctx)

Sets `(vec1, len2)` to $-\mathit{vec2}$.

::: details
Generated from the same template as `_fq_vec_neg`, whose documentation has the details.
:::
:::

## Comparison {#comparison}

::: function
int _fq_nmod_vec_equal(const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len, const fq_nmod_ctx_t ctx)

Returns whether the two vectors of length `len` are equal.

::: details
Generated from the same template as `_fq_vec_equal`, whose documentation has the details.
:::
:::

::: function
int _fq_nmod_vec_is_zero(const fq_nmod_struct * vec, slong len, const fq_nmod_ctx_t ctx)

Returns whether all `len` elements of `vec` are zero.

::: details
Generated from the same template as `_fq_vec_is_zero`, whose documentation has the details.
:::
:::

## Arithmetic {#arithmetic}

::: function
void _fq_nmod_vec_add(fq_nmod_struct * res, const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const fq_nmod_ctx_t ctx)
void _fq_nmod_vec_sub(fq_nmod_struct * res, const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const fq_nmod_ctx_t ctx)

Set `(res, len2)` to $\mathit{vec1} + \mathit{vec2}$, respectively
$\mathit{vec1} - \mathit{vec2}$, elementwise.

::: details
Generated from the same template as `_fq_vec_add`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_vec_scalar_mul_fq_nmod(fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const fq_nmod_t x, const fq_nmod_ctx_t ctx)

Sets `(vec1, len2)` to $x \cdot \mathit{vec2}$.

::: details
Generated from the same template as `_fq_vec_scalar_mul_fq`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_vec_scalar_mul_fmpz(fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const fmpz_t x, const fq_nmod_ctx_t ctx)
void _fq_nmod_vec_scalar_mul_ui(fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const ulong x, const fq_nmod_ctx_t ctx)
void _fq_nmod_vec_scalar_mul_si(fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const slong x, const fq_nmod_ctx_t ctx)

Set `(vec1, len2)` to $\mathit{vec2}$ times the image of the integer $x$
in $\mathbb{F}_p$.

::: details
Generated from the same template as `_fq_vec_scalar_mul_fmpz`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_vec_scalar_addmul_fq_nmod(fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
void _fq_nmod_vec_scalar_submul_fq_nmod(fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const fq_nmod_t x, const fq_nmod_ctx_t ctx)

Add $x \cdot \mathit{vec2}$ to, respectively subtract it from,
`(vec1, len2)`.

::: details
Generated from the same template as `_fq_vec_scalar_addmul_fq`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_vec_dot(fq_nmod_t res, const fq_nmod_struct * vec1, const fq_nmod_struct * vec2, slong len2, const fq_nmod_ctx_t ctx)

Sets `res` to the dot product $\sum_{i < \mathit{len2}}
\mathit{vec1}_i \mathit{vec2}_i$, which is zero for $\mathit{len2} \le 0$.

::: details
Generated from the same template as `_fq_vec_dot`, whose documentation has the details.
:::
:::

## Random vectors {#random}

::: function
void _fq_nmod_vec_randtest(fq_nmod_struct * f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)

Sets the `len` elements of `f` to random elements for testing (with
probability 1/2 most of them are zero).

::: details
Generated from the same template as `_fq_vec_randtest`, whose documentation has the details.
:::
:::

## Input and output {#io}

::: function
int _fq_nmod_vec_fprint(FILE * file, const fq_nmod_struct * vec, slong len, const fq_nmod_ctx_t ctx)
int _fq_nmod_vec_print(const fq_nmod_struct * vec, slong len, const fq_nmod_ctx_t ctx)

Print the vector to `file`, respectively `stdout`: the length, then for a
nonempty vector two spaces and the elements in the format of `fq_nmod_fprint`
separated by single spaces. Return a positive value on success and a
non-positive value on failure.

::: details
Generated from the same template as `_fq_vec_fprint`, whose documentation has the details.
:::
:::

## Macros {#macros}

::: macro
FQ_NMOD_VEC_NORM(vec, i, ctx)

Decreases the `slong` variable `i` while `i` is positive and the element
`vec[i - 1]` is zero, so that `(vec, i)` becomes a normalised polynomial.
:::

::: macro
FQ_NMOD_VEC_SWAP(vec1, len1, vec2, len2)

Swaps the pointers `vec1` and `vec2` and the lengths `len1` and `len2`.
:::
