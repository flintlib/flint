---
title: fmpq_vec.h
header: fmpq_vec.h
description: Vectors of rational numbers stored as arrays of fmpq.
category: Rational numbers
order: 2
---

This module provides a few functions on vectors of rational numbers,
stored as C arrays of `fmpq` of a given length. Each entry is a canonical
`fmpq` (see [](fmpq/canonical-form)). Polynomials and matrices over $\Q$
use other representations: `fmpq_poly` stores an integer polynomial with
a common denominator, and `fmpq_mat` a matrix of `fmpq`.

## Memory management {#memory}

::: function
fmpq * _fmpq_vec_init(slong len)

Allocates a vector of `len` entries and initialises them to zero.
:::

::: function
void _fmpq_vec_clear(fmpq * vec, slong len)

Clears the `len` entries and frees the vector.
:::

## Random generation {#random}

::: function
void _fmpq_vec_randtest(fmpq * f, flint_rand_t state, slong len, flint_bitcnt_t bits)

Sets the entries to random canonical fractions as given by
`fmpq_randtest`, with numerators and denominators of at most `bits` bits.
With probability $1/2$ the vector is made sparse, with a random density
of nonzero entries.
:::

::: function
void _fmpq_vec_randtest_uniq_sorted(fmpq * vec, flint_rand_t state, slong len, flint_bitcnt_t bits)

Sets the entries to distinct random fractions as given by
`fmpq_randtest`, sorted in increasing order. Throws if `bits` is less than
four times the bit length of `len`, to make collisions unlikely.
:::

## Heights {#heights}

::: function
void _fmpq_vec_max_height(fmpz_t height, const fmpq * vec, slong len)
flint_bitcnt_t _fmpq_vec_max_height_bits(const fmpq * vec, slong len)

Set `height` to, or return the bit length of, the maximum of the heights
$\max(|p|, q)$ of the entries $p/q$ (see `fmpq_height`). Both are zero for
an empty vector.
:::

## Comparison and sorting {#comparison}

::: function
int _fmpq_vec_equal(const fmpq * vec1, const fmpq * vec2, slong len)

Returns 1 if the vectors have equal entries, otherwise 0.
:::

::: function
void _fmpq_vec_sort(fmpq * vec, slong len)

Sorts the entries in increasing order, with `qsort` and `fmpq_cmp`.
:::

## Conversion {#conversion}

::: function
void _fmpq_vec_set_fmpz_vec(fmpq * res, const fmpz * vec, slong len)

Sets the entries of `res` to the integers in `vec`.
:::

::: function
void _fmpq_vec_get_fmpz_vec_fmpz(fmpz * num, fmpz_t den, const fmpq * a, slong len)

Writes the vector `a` with a common denominator: sets `den` to the least
common multiple of the denominators of the entries (1 for an empty
vector) and `num[i]` to $a_i \cdot \mathit{den}$.
:::

## Dot product {#dot}

::: function
void _fmpq_vec_dot(fmpq_t res, const fmpq * vec1, const fmpq * vec2, slong len)

Sets `res` to $\sum_{i < \mathit{len}} \mathit{vec1}_i \,
\mathit{vec2}_i$, accumulated with `fmpq_addmul`. `res` must not be an
entry of the input vectors.

::: details
Each step costs GCDs of the size of the partial sum. For long vectors it
is faster to clear denominators with `_fmpq_vec_get_fmpz_vec_fmpz` and
compute an integer dot product with `_fmpz_vec_dot`.
:::
:::

## Input and output {#io}

::: function
int _fmpq_vec_fprint(FILE * file, const fmpq * vec, slong len)
int _fmpq_vec_print(const fmpq * vec, slong len)

Print the vector to `file` or `stdout`: the length, two spaces, and the
entries separated by single spaces, in the format of `fmpq_print`. An
empty vector is printed as `0`. Return a positive value on success and a
non-positive value on failure.
:::
