---
title: fq_vec.h
header: fq_vec.h
template-headers: fq_vec_templates
description: Vectors over finite fields $\mathbb{F}_q$ of arbitrary characteristic, stored as arrays of `fq_struct`.
category: Finite fields
order: 6
---

This module provides low-level functions on vectors over a finite field
$\mathbb{F}_q$ in the representation of [](fq): a vector is a plain
C array of `fq_struct` together with its length, and memory management
is up to the caller. The functions are used by `fq_poly` and `fq_mat`.
They are generated from the template header `fq_vec_templates.h`, which
also defines the same functions for `fq_nmod_vec` and `fq_zech_vec`.

All functions allow the output vector to be the same as an input vector,
but not to overlap partially.

## Memory management {#memory}

::: function
fq_struct * _fq_vec_init(slong len, const fq_ctx_t ctx)

Returns a newly allocated vector of `len` elements, all initialised to
zero.
:::

::: function
void _fq_vec_clear(fq_struct * vec, slong len, const fq_ctx_t ctx)

Clears the `len` elements of `vec` and frees the array.
:::

## Assignment and basic manipulation {#assignment}

::: function
void _fq_vec_set(fq_struct * v, const fq_struct * f, slong len, const fq_ctx_t ctx)

Copies `(f, len)` into `(v, len)`.
:::

::: function
void _fq_vec_swap(fq_struct * vec1, fq_struct * vec2, slong len2, const fq_ctx_t ctx)

Swaps the elements of `(vec1, len2)` and `(vec2, len2)` efficiently.
:::

::: function
void _fq_vec_zero(fq_struct * v, slong len, const fq_ctx_t ctx)

Sets the `len` elements of `v` to zero.
:::

::: function
void _fq_vec_neg(fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_ctx_t ctx)

Sets `(vec1, len2)` to $-\mathit{vec2}$.
:::

## Comparison {#comparison}

::: function
int _fq_vec_equal(const fq_struct * vec1, const fq_struct * vec2, slong len, const fq_ctx_t ctx)

Returns whether the two vectors of length `len` are equal.
:::

::: function
int _fq_vec_is_zero(const fq_struct * vec, slong len, const fq_ctx_t ctx)

Returns whether all `len` elements of `vec` are zero.
:::

## Arithmetic {#arithmetic}

::: function
void _fq_vec_add(fq_struct * res, const fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_ctx_t ctx)
void _fq_vec_sub(fq_struct * res, const fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_ctx_t ctx)

Set `(res, len2)` to $\mathit{vec1} + \mathit{vec2}$, respectively
$\mathit{vec1} - \mathit{vec2}$, elementwise.
:::

::: function
void _fq_vec_scalar_mul_fq(fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_t x, const fq_ctx_t ctx)

Sets `(vec1, len2)` to $x \cdot \mathit{vec2}$.
:::

::: function
void _fq_vec_scalar_mul_fmpz(fq_struct * vec1, const fq_struct * vec2, slong len2, const fmpz_t x, const fq_ctx_t ctx)
void _fq_vec_scalar_mul_ui(fq_struct * vec1, const fq_struct * vec2, slong len2, const ulong x, const fq_ctx_t ctx)
void _fq_vec_scalar_mul_si(fq_struct * vec1, const fq_struct * vec2, slong len2, const slong x, const fq_ctx_t ctx)

Set `(vec1, len2)` to $\mathit{vec2}$ times the image of the integer $x$
in $\mathbb{F}_p$.
:::

::: function
void _fq_vec_scalar_addmul_fq(fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_t x, const fq_ctx_t ctx)
void _fq_vec_scalar_submul_fq(fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_t x, const fq_ctx_t ctx)

Add $x \cdot \mathit{vec2}$ to, respectively subtract it from,
`(vec1, len2)`.
:::

::: function
void _fq_vec_dot(fq_t res, const fq_struct * vec1, const fq_struct * vec2, slong len2, const fq_ctx_t ctx)

Sets `res` to the dot product $\sum_{i < \mathit{len2}}
\mathit{vec1}_i \mathit{vec2}_i$, which is zero for $\mathit{len2} \le 0$.

::: details
This is a wrapper of the generic `_gr_vec_dot`. The `fq` implementation
accumulates the unreduced products of the representatives as one integer
polynomial and reduces it once at the end, modulo $p$ and modulo the
defining polynomial, so a dot product costs about as much as $\ell$
polynomial products and one reduction instead of $\ell$ reductions. The
`fq_nmod` version accumulates the products modulo $p$ and also reduces
modulo the defining polynomial once; the `fq_zech` version works with
discrete logarithms and the Zech table (see [](fq_zech_vec)).
:::
:::

## Random vectors {#random}

::: function
void _fq_vec_randtest(fq_struct * f, flint_rand_t state, slong len, const fq_ctx_t ctx)

Sets the `len` elements of `f` to random elements for testing (with
probability 1/2 most of them are zero).
:::

## Input and output {#io}

::: function
int _fq_vec_fprint(FILE * file, const fq_struct * vec, slong len, const fq_ctx_t ctx)
int _fq_vec_print(const fq_struct * vec, slong len, const fq_ctx_t ctx)

Print the vector to `file`, respectively `stdout`: the length, then for a
nonempty vector two spaces and the elements in the format of `fq_fprint`
separated by single spaces. Return a positive value on success and a
non-positive value on failure.
:::

## Macros {#macros}

::: macro
FQ_VEC_NORM(vec, i, ctx)

Decreases the `slong` variable `i` while `i` is positive and the element
`vec[i - 1]` is zero, so that `(vec, i)` becomes a normalised polynomial.
:::

::: macro
FQ_VEC_SWAP(vec1, len1, vec2, len2)

Swaps the pointers `vec1` and `vec2` and the lengths `len1` and `len2`.
:::
