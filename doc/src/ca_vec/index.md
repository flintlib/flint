---
title: ca_vec.h
header: ca_vec.h
description: Vectors of exact real and complex numbers: memory management, entrywise arithmetic and zero tests.
category: Exact real and complex numbers
order: 3
---

This module implements vectors of exact real and complex numbers of type
`ca_t` (see [](ca)). A `ca_vec_t` is a growable array with a length that
is set by the user: unlike a [](ca_poly) polynomial, it is never
normalised by removing zero entries. The module is small; it supplies
the memory management and the entrywise operations used by [](ca_poly)
and [](ca_mat).

As in other FLINT modules, functions with a leading underscore work on
plain arrays of type `ca_ptr` (or `ca_srcptr` for read-only input) with
an explicit length, and leave memory management to the caller. Every
function takes the context object `ctx` that the entries belong to.
Entries are arbitrary `ca_t` values, including special values such as
*Unknown*, and arithmetic follows the rules of the scalar functions in
[](ca).

## Types {#types}

::: type
ca_vec_struct
ca_vec_t

A vector of `ca_t` numbers. The structure has the fields `entries` (an
array of `alloc` initialised entries), `length` and `alloc`. A
`ca_vec_t` is an array of one structure, so it is passed by reference.
All `alloc` entries are initialised; those beyond `length` are zero.
:::

::: macro
ca_vec_entry(vec, i)

Expands to a pointer (`ca_ptr`) to entry `i` of `vec`. The index must
be nonnegative and less than the length of `vec`; this is not checked.
:::

::: function
ca_ptr ca_vec_entry_ptr(ca_vec_t vec, slong i)

Returns a pointer to entry `i` of `vec`, like the macro `ca_vec_entry`.
:::

## Memory management {#memory}

::: function
ca_ptr _ca_vec_init(slong len, ca_ctx_t ctx)

Allocates an array of `len` entries, each initialised to zero, and
returns a pointer to it. Free it with `_ca_vec_clear`.
:::

::: function
void ca_vec_init(ca_vec_t vec, slong len, ca_ctx_t ctx)

Initialises `vec` as a vector of length `len` with all entries zero. For
`len` = 0 nothing is allocated.
:::

::: function
void _ca_vec_clear(ca_ptr vec, slong len, ca_ctx_t ctx)

Clears the `len` entries of `vec` and frees the array itself.
:::

::: function
void ca_vec_clear(ca_vec_t vec, ca_ctx_t ctx)

Clears `vec`, freeing all its entries.
:::

::: function
void _ca_vec_swap(ca_ptr vec1, ca_ptr vec2, slong len, ca_ctx_t ctx)
void ca_vec_swap(ca_vec_t vec1, ca_vec_t vec2, ca_ctx_t ctx)

Swaps the contents of the two vectors, entry by entry for the
underscore version (which requires arrays of the same length `len`),
and by exchanging the structures in constant time for `ca_vec_swap`.
:::

## Length {#length}

::: function
slong ca_vec_length(const ca_vec_t vec, ca_ctx_t ctx)

Returns the length of `vec`.
:::

::: function
void _ca_vec_fit_length(ca_vec_t vec, slong len, ca_ctx_t ctx)

Ensures that `vec` has space for at least `len` entries, without
changing its length. If the array has to grow, it grows to at least
twice its previous size, and the new entries are initialised to zero.
:::

::: function
void ca_vec_set_length(ca_vec_t vec, slong len, ca_ctx_t ctx)

Sets the length of `vec` to `len`. A shorter vector is extended with
zeros; a longer vector is truncated, and the removed entries are set to
zero (their memory is kept for reuse).
:::

## Assignment and special vectors {#assignment}

::: function
void _ca_vec_set(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx)
void ca_vec_set(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx)

Sets `res` to a copy of `src`. The vector version also sets the length
of `res` to that of `src`.
:::

::: function
void _ca_vec_zero(ca_ptr res, slong len, ca_ctx_t ctx)
void ca_vec_zero(ca_vec_t res, slong len, ca_ctx_t ctx)

Sets the `len` entries of `res` to zero. The vector version also sets
the length of `res` to `len`.
:::

::: function
void _ca_vec_unknown(ca_ptr vec, slong len, ca_ctx_t ctx)
void _ca_vec_undefined(ca_ptr vec, slong len, ca_ctx_t ctx)

Sets the `len` entries of `vec` to the special value *Unknown* or
*Undefined*. Functions of [](ca_poly) and [](ca_mat) use these to signal
that a result could not be computed.
:::

::: function
void ca_vec_append(ca_vec_t vec, const ca_t f, ca_ctx_t ctx)

Appends a copy of `f` at the end of `vec`, increasing its length by one.
The amortised cost is constant thanks to the doubling strategy of
`_ca_vec_fit_length`.
:::

## Input and output {#io}

::: function
void ca_vec_print(const ca_vec_t vec, ca_ctx_t ctx)

Prints the line `ca_vec of length n:` followed by the entries of `vec`,
one per line, each printed with `ca_print` and indented by four spaces,
and a blank line.
:::

::: function
void ca_vec_printn(const ca_vec_t vec, slong digits, ca_ctx_t ctx)

Prints decimal approximations of the entries with `digits` significant
digits (see `ca_printn`), comma-separated in square brackets, followed
by a newline.
:::

## Arithmetic {#arithmetic}

The arithmetic functions apply the corresponding scalar function of
[](ca) to each entry. The underscore versions require all arrays to have
length `len` and allow the output to be the same array as an input.

::: function
void _ca_vec_neg(ca_ptr res, ca_srcptr src, slong len, ca_ctx_t ctx)
void ca_vec_neg(ca_vec_t res, const ca_vec_t src, ca_ctx_t ctx)

Sets `res` to the negation of `src`. The vector version sets the length
of `res` to that of `src`. Aliasing of `res` and `src` is allowed.
:::

::: function
void _ca_vec_add(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx)
void _ca_vec_sub(ca_ptr res, ca_srcptr vec1, ca_srcptr vec2, slong len, ca_ctx_t ctx)

Sets `res` to the entrywise sum or difference of `vec1` and `vec2`.
:::

::: function
void _ca_vec_scalar_mul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)
void _ca_vec_scalar_div_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)

Sets `res` to `src` multiplied or divided by the scalar `c`. Division
computes $c^{-1}$ once with `ca_inv` and multiplies by it; if $c = 0$
the entries become unsigned infinity or *Undefined* (for zero entries),
and *Unknown* if it cannot be decided whether $c = 0$.
:::

::: function
void _ca_vec_scalar_addmul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)
void _ca_vec_scalar_submul_ca(ca_ptr res, ca_srcptr src, slong len, const ca_t c, ca_ctx_t ctx)

Adds `src` multiplied by `c` to `res`, or subtracts it from `res`.
:::

## Comparisons and predicates {#predicates}

::: function
truth_t _ca_vec_check_is_zero(ca_srcptr vec, slong len, ca_ctx_t ctx)

Checks whether all `len` entries are zero, using `ca_check_is_zero`.
Returns `T_FALSE` as soon as one entry is proved nonzero, `T_TRUE` if
all entries are proved zero, and `T_UNKNOWN` otherwise. The empty vector
is zero.
:::

## Conversions for rational vectors {#rational}

The following helpers are used to pass vectors whose entries happen to
be rational numbers to faster algorithms over $\mathbb{Q}$ or $\Z$. They
inspect the representation only: a number that is rational but is
represented as an element of a larger field (see `ca_is_qq_elem`) is not
recognised.

::: function
int _ca_vec_is_fmpq_vec(ca_srcptr vec, slong len, ca_ctx_t ctx)

Returns whether every entry of `vec` is represented as an element of the
field $\mathbb{Q}$.
:::

::: function
int _ca_vec_fmpq_vec_is_fmpz_vec(ca_srcptr vec, slong len, ca_ctx_t ctx)

Assuming that every entry is represented as a rational number (see
`_ca_vec_is_fmpq_vec`), returns whether every entry is an integer.
:::

::: function
void _ca_vec_fmpq_vec_get_fmpz_vec_den(fmpz * c, fmpz_t den, ca_srcptr vec, slong len, ca_ctx_t ctx)

Assuming that every entry is represented as a rational number, sets
`den` to the least common multiple of the denominators and `c` to the
integer vector with $\mathit{vec}_i = c_i / \mathit{den}$. The array `c`
must have room for `len` entries.
:::

::: function
void _ca_vec_set_fmpz_vec_div_fmpz(ca_ptr res, const fmpz * v, const fmpz_t den, slong len, ca_ctx_t ctx)

Sets the entries of `res` to the rational numbers $v_i / \mathit{den}$.
The denominator must be nonzero.
:::
