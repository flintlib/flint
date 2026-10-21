---
title: gr_vec.h
header: gr_vec.h
description: Vectors over generic rings: memory management, elementwise arithmetic, dot products, sums and products, sorting and parsing.
category: Generic rings
order: 5
---

This module provides vectors over the rings of the generic interface
[`gr`](gr). Most functions act on plain arrays of `len` contiguous
elements, passed as a `gr_ptr` or `gr_srcptr` with an explicit length
and a leading underscore in the name; the type `gr_vec_t` is a resizable
vector that manages its own memory and length, used for example for the
factors returned by `gr_factor` and for lists of roots. The underscore
functions are the second level of the three levels of generic
operations ([](algorithms/generic/method-tables#levels)): each of them
dispatches to a method of the ring, whose default loops over the
elementwise methods, and rings with cheap elements override the
important ones (addition, scalar multiplication, dot products) with
specialised loops. Generic algorithms for polynomials and matrices are
written in terms of these functions, so their speed depends mostly on
the vector methods of the ring.

Unless stated otherwise, the output array may be the same as an input
array (complete aliasing) but must not overlap it partially, and a
length of 0 is allowed and does nothing. The guide [](gr/getting-started#vectors) shows the
basic use.

## Types and memory management {#types}

::: type
gr_vec_struct
gr_vec_t

A vector with fields `entries` (the array of elements), `alloc` (the
number of allocated, initialised entries) and `length`. A `gr_vec_t` is
an array of length one. The context of the entries is not stored and
must be passed to every function.
:::


::: function
void gr_vec_init(gr_vec_t vec, slong len, gr_ctx_t ctx)

Initialises `vec` to a vector of length `len` with all entries zero; `len` must be nonnegative.
:::

::: function
void gr_vec_clear(gr_vec_t vec, gr_ctx_t ctx)

Clears the vector, clearing all allocated entries and freeing the array.
:::

::: macro
GR_VEC_ENTRY(vec, i, sz)

The address of entry $i$ of the `gr_vec_t` `vec`, whose elements have
`sz` bytes. The index is not checked.
:::


::: function
gr_ptr gr_vec_entry_ptr(gr_vec_t vec, slong i, gr_ctx_t ctx)
gr_srcptr gr_vec_entry_srcptr(const gr_vec_t vec, slong i, gr_ctx_t ctx)

Return a pointer to entry $i$ of `vec`, indexed from zero. The index must be in bounds (not checked).
:::

::: function
slong gr_vec_length(const gr_vec_t vec, gr_ctx_t ctx)

Returns the length of `vec`.
:::

::: function
void gr_vec_fit_length(gr_vec_t vec, slong len, gr_ctx_t ctx)

Ensures that at least `len` entries are allocated, without changing the length; the allocation at least doubles when it grows. Allocated entries beyond the length are initialised (to zero).
:::

::: function
void gr_vec_set_length(gr_vec_t vec, slong len, gr_ctx_t ctx)

Sets the length to `len`, which must be nonnegative. New entries are zero; entries beyond the new length are set to zero, freeing the memory they own, but stay allocated.
:::

::: function
int gr_vec_set(gr_vec_t res, const gr_vec_t src, gr_ctx_t ctx)

Sets `res` to a copy of `src`, resizing `res`.
:::

::: function
void gr_vec_set_shallow(gr_vec_t res, const gr_vec_t src, gr_ctx_t ctx)

Resizes `res` to the length of `src` and copies the bytes of the entries of `src` into it, without copying owned data, so that both vectors share the memory owned by the entries. The previous contents of `res` are overwritten without being freed. For temporary views; `res` must not be cleared while `src` is in use, and not at all if the entries own memory.
:::

::: function
int gr_vec_append(gr_vec_t vec, gr_srcptr f, gr_ctx_t ctx)
void gr_vec_append_swap(gr_vec_t vec, gr_ptr f, gr_ctx_t ctx)

Append a copy of `f` to `vec`, respectively append `f` by swapping it into the new last entry (leaving `f` with the previous value of that entry, which is zero if it was never used). `f` must not point into `vec`, since appending may reallocate the array.
:::

::: function
void _gr_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx)
void _gr_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx)

Initialise the `len` elements of the array `vec` to zero, respectively clear them, freeing the memory they own but not the array itself. The array memory is allocated by the caller (see also `gr_heap_init_vec` and `GR_TMP_INIT_VEC`).
:::

::: function
void _gr_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx)

Swaps the entries of `vec1` and `vec2` elementwise.
:::

## Assignment and data movement {#assignment}


::: function
int _gr_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx)

Sets the `len` entries of `vec` to zero.
:::

::: function
int _gr_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)
int _gr_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)

Set `res` to a copy of `src`, respectively its negation, elementwise.
:::

::: function
void _gr_vec_set_shallow(gr_ptr dest, gr_srcptr src, slong len, gr_ctx_t ctx)

Copies the bytes of `len` elements from `src` to `dest` with `memcpy`, without copying owned data (see `gr_set_shallow`).
:::

::: function
int _gr_vec_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, slong len, gr_ctx_t ctx)

Sets `res[i]` to `x[i]` converted from the ring `x_ctx`, as `gr_set_other` elementwise (method `GR_METHOD_VEC_SET_OTHER`). If `x_ctx` is `ctx`, this is `_gr_vec_set`.
:::

::: function
int _gr_vec_gather(gr_ptr res, gr_srcptr vec, const slong * idx, slong len, gr_ctx_t ctx)
int _gr_vec_scatter(gr_ptr vec, const slong * idx, gr_srcptr src, slong len, gr_ctx_t ctx)

Indexed copies: `res[i] = vec[idx[i]]`, respectively `vec[idx[i]] = src[i]`, for $0 \le i < \mathit{len}$ (methods `GR_METHOD_VEC_GATHER` and `GR_METHOD_VEC_SCATTER`). The indices need not be distinct; when they repeat in a scatter, the last write wins. The defaults loop over `gr_set`; rings whose elements own no memory override them with plain copies.
:::

::: function
int _gr_vec_randtest(gr_ptr res, flint_rand_t state, slong len, gr_ctx_t ctx)

Sets each entry to zero with probability 1/2 and to a random element (`gr_randtest`) otherwise.
:::

::: function
int _gr_vec_step(gr_ptr vec, gr_srcptr start, gr_srcptr step, slong len, gr_ctx_t ctx)

Sets `vec` to the arithmetic progression $s, s + d, s + 2d, \ldots$ of length `len` with $s$ = `start` and $d$ = `step`, by repeated addition.
:::

::: function
int _gr_vec_set_powers(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)

Sets `res` to the powers $1, x, x^2, \ldots, x^{len-1}$.

::: details
The default multiplies by $x$ repeatedly, except in finite rings and rings with a real precision, where it computes even powers by squaring ($x^{2k} = (x^k)^2$, $x^{2k+1} = x^{k+1} x^k$), which reduces the error growth of floating-point and ball arithmetic.
:::
:::

::: function
int _gr_vec_reciprocals(gr_ptr res, slong len, gr_ctx_t ctx)

Sets `res` to the reciprocals $1, 1/2, 1/3, \ldots, 1/\mathit{len}$ of the first `len` positive integers. The default stops at the first failure, so in $\Z/n\Z$ it returns `GR_DOMAIN` if some $k \le \mathit{len}$ is not invertible.
:::

## Comparisons and normalisation {#comparisons}


::: function
truth_t _gr_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

Returns whether the two arrays are equal elementwise: `T_FALSE` if some pair is known to differ, `T_TRUE` if all are known equal, and `T_UNKNOWN` otherwise.
:::

::: function
truth_t _gr_vec_is_zero(gr_srcptr vec, slong len, gr_ctx_t ctx)

Returns whether all entries are zero, in the same three-valued sense.
:::

::: function
truth_t _gr_vec_contains(gr_srcptr vec, slong len, gr_srcptr x, gr_ctx_t ctx)
truth_t gr_vec_contains(const gr_vec_t vec, gr_srcptr x, gr_ctx_t ctx)

Return whether `x` equals some entry: `T_TRUE` if an entry is known equal, `T_FALSE` if all are known different, `T_UNKNOWN` otherwise.
:::

::: function
int _gr_vec_normalise(slong * res, gr_srcptr vec, slong len, gr_ctx_t ctx)

Sets `*res` to the length of `vec` with trailing zero entries removed: the largest $n \le \mathit{len}$ such that entry $n - 1$ is nonzero. Returns `GR_UNABLE` if `gr_is_zero` returns `T_UNKNOWN` for an entry that would decide the length; `*res` is then the length up to that entry.
:::

::: function
slong _gr_vec_normalise_weak(gr_srcptr vec, slong len, gr_ctx_t ctx)

Returns the length with trailing entries removed while they are known to be zero (`gr_is_zero` returns `T_TRUE`); the result is only an upper bound for the true length in inexact rings. Polynomials over inexact rings are normalised this way ([](algorithms/generic/inexact#weak-normalisation)).
:::

## Elementwise arithmetic {#arithmetic}

These functions apply a binary operation elementwise. In the scalar
versions one operand is the same for all entries. In the `_other`
versions one operand (a vector or a scalar) belongs to another ring,
`ctx2`, `ctx3` or `cctx`, and is converted as by `gr_set_other`. The
status is the OR of the elementwise statuses; the defaults continue
after a failure.


::: function
int _gr_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int _gr_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int _gr_vec_mul(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int _gr_vec_div(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int _gr_vec_divexact(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
int _gr_vec_pow(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)

Set `res[i]` to `src1[i]` op `src2[i]`.
:::

::: function
int _gr_vec_add_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int _gr_vec_sub_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int _gr_vec_mul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int _gr_vec_div_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int _gr_vec_divexact_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int _gr_vec_pow_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)

Set `vec1[i]` to `vec2[i]` op $c$.
:::

::: function
int _gr_vec_add_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int _gr_vec_sub_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int _gr_vec_mul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int _gr_vec_div_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int _gr_vec_divexact_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int _gr_vec_pow_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int _gr_vec_add_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int _gr_vec_sub_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int _gr_vec_mul_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int _gr_vec_div_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int _gr_vec_divexact_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int _gr_vec_pow_scalar_ui(gr_ptr vec1, gr_srcptr vec2, slong len, ulong c, gr_ctx_t ctx)
int _gr_vec_add_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int _gr_vec_sub_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int _gr_vec_mul_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int _gr_vec_div_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int _gr_vec_divexact_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int _gr_vec_pow_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)
int _gr_vec_add_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int _gr_vec_sub_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int _gr_vec_mul_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int _gr_vec_div_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int _gr_vec_divexact_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)
int _gr_vec_pow_scalar_fmpq(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpq_t c, gr_ctx_t ctx)

Set `vec1[i]` to `vec2[i]` op $c$ for a scalar $c$ of type `slong`, `ulong`, `fmpz_t` or `fmpq_t`.
:::

::: function
int _gr_scalar_add_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_sub_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_mul_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_div_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_pow_vec(gr_ptr vec1, gr_srcptr c, gr_srcptr vec2, slong len, gr_ctx_t ctx)

Set `vec1[i]` to $c$ op `vec2[i]`.
:::

::: function
int _gr_vec_add_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int _gr_vec_sub_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int _gr_vec_mul_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int _gr_vec_div_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int _gr_vec_divexact_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)
int _gr_vec_pow_other(gr_ptr vec1, gr_srcptr vec2, gr_srcptr vec3, gr_ctx_t ctx3, slong len, gr_ctx_t ctx)

Set `vec1[i]` to `vec2[i]` op `vec3[i]`, where `vec3` is an array over `ctx3`.
:::

::: function
int _gr_other_add_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int _gr_other_sub_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int _gr_other_mul_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int _gr_other_div_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int _gr_other_divexact_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)
int _gr_other_pow_vec(gr_ptr vec1, gr_srcptr vec2, gr_ctx_t ctx2, gr_srcptr vec3, slong len, gr_ctx_t ctx)

Set `vec1[i]` to `vec2[i]` op `vec3[i]`, where `vec2` is an array over `ctx2`.
:::

::: function
int _gr_vec_add_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int _gr_vec_sub_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int _gr_vec_mul_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int _gr_vec_div_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int _gr_vec_divexact_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)
int _gr_vec_pow_scalar_other(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t cctx, gr_ctx_t ctx)

Set `vec1[i]` to `vec2[i]` op $c$, where $c$ is an element of `cctx`.
:::

::: function
int _gr_scalar_other_add_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_other_sub_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_other_mul_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_other_div_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_other_divexact_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_scalar_other_pow_vec(gr_ptr vec1, gr_srcptr c, gr_ctx_t cctx, gr_srcptr vec2, slong len, gr_ctx_t ctx)

Set `vec1[i]` to $c$ op `vec2[i]`, where $c$ is an element of `cctx`.
:::

::: function
int _gr_vec_addmul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int _gr_vec_submul_scalar(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
int _gr_vec_addmul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int _gr_vec_submul_scalar_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
int _gr_vec_addmul_scalar_fmpz(gr_ptr vec1, gr_srcptr vec2, slong len, const fmpz_t c, gr_ctx_t ctx)

Set `vec1[i]` to `vec1[i]` $\pm$ `vec2[i]` $\cdot c$ (an *axpy* operation). These are the inner loops of classical polynomial division and Gaussian elimination, and rings with cheap elements override them.
:::

::: function
int _gr_vec_mul_scalar_2exp_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)

Sets `vec1[i]` to `vec2[i]` $\cdot 2^c$, as `gr_mul_2exp_si` elementwise.
:::

## Elementary functions {#elementary}


::: function
int _gr_vec_sqrt(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int _gr_vec_rsqrt(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int _gr_vec_exp(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int _gr_vec_log(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int _gr_vec_sin(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int _gr_vec_cos(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx)
int _gr_vec_sin_cos(gr_ptr res1, gr_ptr res2, gr_srcptr x, slong len, gr_ctx_t ctx)

Elementwise square roots, reciprocal square roots, exponentials, logarithms, sines and cosines (methods `GR_METHOD_VEC_SQRT` and so on). The defaults loop over the scalar functions of [`gr_special`](gr_special); rings with vectorised kernels (`dfloat`) override them.
:::

::: function
int _gr_vec_get_interval_mid_rad(gr_ptr m, gr_ptr r, gr_srcptr x, slong len, gr_ctx_t ctx)

Applies `gr_get_interval_mid_rad` elementwise (method `GR_METHOD_VEC_GET_INTERVAL_MID_RAD`).
:::

## Sums and products {#sum}


::: function
int _gr_vec_sum(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
int _gr_vec_product(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)

Set `res` to the sum, respectively the product, of the `len` entries; the empty sum is 0 and the empty product 1 (methods `GR_METHOD_VEC_SUM`, `GR_METHOD_VEC_PRODUCT`).

::: details
The defaults are `_gr_vec_sum_generic` and `_gr_vec_product_generic`. In finite rings they add or multiply from left to right. In other rings, where entries can grow, the sum of more than 100 entries uses binary splitting with basecase length 100, and the product of more than 20 entries binary splitting with basecase length 20, or, for more than 500 entries in a thread-safe ring, parallel binary splitting with basecase length 500. Balanced binary splitting computes a product of $n$ integers of $b$ bits each in $O(M(nb) \log n)$ time instead of the $O(n^2 b)$ of the left-to-right loop, and similarly sums of fractions; see [](algorithms/integers/binary-splitting).
:::
:::

::: function
int _gr_vec_sum_generic(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
int _gr_vec_product_generic(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)

The default implementations described above.
:::

::: function
int _gr_vec_sum_serial(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
int _gr_vec_product_serial(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)

Add or multiply the entries from left to right.
:::

::: function
int _gr_vec_sum_bsplit(gr_ptr res, gr_srcptr vec, slong len, slong basecase_cutoff, gr_ctx_t ctx)
int _gr_vec_product_bsplit(gr_ptr res, gr_srcptr vec, slong len, slong basecase_cutoff, gr_ctx_t ctx)

Binary splitting: split the array in two halves recursively, down to subarrays of fewer than `basecase_cutoff` entries, which are summed or multiplied with `_gr_vec_sum` or `_gr_vec_product`.
:::

::: function
int _gr_vec_sum_bsplit_parallel(gr_ptr res, gr_srcptr vec, slong len, slong basecase_cutoff, gr_ctx_t ctx)
int _gr_vec_product_bsplit_parallel(gr_ptr res, gr_srcptr vec, slong len, slong basecase_cutoff, gr_ctx_t ctx)

Binary splitting with the recursion distributed over threads with `flint_parallel_binary_splitting`; the ring must be thread-safe.
:::

::: function
int _gr_vec_sum_parallel(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)
int _gr_vec_product_parallel(gr_ptr res, gr_srcptr vec, slong len, gr_ctx_t ctx)

Split the array into one chunk per available thread, compute the chunks in parallel with `_gr_vec_sum` or `_gr_vec_product`, and combine the partial results. The number of threads is at most `flint_get_num_threads()` and at most $\mathit{len}/2$.
:::

## Dot products {#dot}

The function page [](gr_vec/dot) explains the arguments, aliasing and
the specialised implementations.


::: function
int _gr_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
int _gr_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)

Set `res` to $c + \sum_{i=0}^{n-1} a_i b_i$, respectively $c + \sum_{i=0}^{n-1} a_i b_{n-1-i}$, or $c - \sum$ if `subtract` is nonzero, where $a$ = `vec1`, $b$ = `vec2`, $n$ = `len` and $c$ = `initial`, which may be `NULL` for $c = 0$.
:::

::: function
int _gr_vec_dot_strided(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, slong stride1, gr_srcptr vec2, slong stride2, slong len, gr_ctx_t ctx)

Sets `res` to $c \pm \sum_{i=0}^{n-1} a_{i s_1} b_{i s_2}$ with strides $s_1$ = `stride1` and $s_2$ = `stride2` (in elements), which may be negative or zero.
:::

::: function
int _gr_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const slong * vec2, slong len, gr_ctx_t ctx)
int _gr_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const ulong * vec2, slong len, gr_ctx_t ctx)
int _gr_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const fmpz * vec2, slong len, gr_ctx_t ctx)

As `_gr_vec_dot`, with the second vector an array of `slong`, `ulong` or `fmpz` integers.
:::

## Sorting and permutations {#sorting}


::: function
int _gr_vec_sort(gr_ptr vec, slong len, gr_ctx_t ctx)
int gr_vec_sort(gr_vec_t dest, const gr_vec_t src, gr_ctx_t ctx)

Sort the entries in increasing order with respect to `gr_cmp`. `_gr_vec_sort` works in place. If a comparison fails, the two entries are treated as equal and the failure is returned; the array is then a permutation of the input but may not be sorted.
:::

::: function
void _gr_vec_permute(gr_ptr vec, slong * perm, slong len, gr_ctx_t ctx)
int gr_vec_permute(gr_vec_t dest, gr_vec_t src, slong * perm, gr_ctx_t ctx)

Apply the permutation `perm` of $\{0, \ldots, \mathit{len}-1\}$ to the entries: `dest[perm[i]] = src[i]`. The underscore version works in place by swapping along cycles and overwrites `perm`; `gr_vec_permute` works on a copy of `perm`.
:::

::: function
void _gr_vec_permute_inv(gr_ptr vec, slong * perm, slong len, gr_ctx_t ctx)
int gr_vec_permute_inv(gr_vec_t dest, gr_vec_t src, slong * perm, gr_ctx_t ctx)

Apply the inverse permutation: `dest[i] = src[perm[i]]`. The underscore version works in place and overwrites `perm`.
:::

::: function
void _gr_vec_shuffle(gr_ptr vec, flint_rand_t state, slong len, gr_ctx_t ctx)

Applies a uniformly random permutation to the entries.
:::

## Input and output {#io}


::: function
int _gr_vec_write(gr_stream_t out, gr_srcptr vec, slong len, gr_ctx_t ctx)
int gr_vec_write(gr_stream_t out, const gr_vec_t vec, gr_ctx_t ctx)
int _gr_vec_print(gr_srcptr vec, slong len, gr_ctx_t ctx)
int gr_vec_print(const gr_vec_t vec, gr_ctx_t ctx)

Write the vector to a stream, or print it to standard output, in the form `[a, b, c]` with each entry written by `gr_write`. The format specifier `%{gr*}` of `flint_printf` does the same for an array ([](flint/printf)).
:::

::: function
int gr_vec_set_str(gr_vec_t vec, const char * s, int resize, gr_ctx_t ctx)

Sets `vec` to the vector described by the string `s`, of the form `[expr1, expr2, ..., exprn]` (or `[]`), parsing each entry with `gr_set_str` over `ctx`.

::: details
Whitespace around the brackets and separators is ignored, and commas inside parentheses, brackets or braces of an entry belong to that entry. If `resize` is 1, `vec` is resized to the number of entries; if `resize` is 0, the function returns `GR_DOMAIN` when that number differs from the current length of `vec`. Returns `GR_UNABLE` if `s` is not a well-formed list, and propagates any failure from parsing an entry.
:::
:::

::: function
int gr_vec_str_count_entries(slong * count, const char * s, gr_ctx_t ctx)

Sets `*count` to the number of top-level entries of the bracketed list at the start of `s`, without evaluating them. Separators are counted only at the outermost nesting level of parentheses, brackets and braces, and characters after the closing bracket are ignored. Returns `GR_UNABLE` if the list is not well-formed.
:::

::: function
int _gr_vec_set_str(gr_ptr res, const char * s, slong len, gr_ctx_t ctx)

Evaluates the bracketed list at the start of `s` into the array `res` of `len` initialised elements, returning `GR_DOMAIN` if the number of entries differs from `len`. Characters after the closing bracket are ignored, so the function can read a row inside a matrix string.
:::

## Internal functions {#internal}


::: function
int _gr_vec_check_resize(gr_vec_t res, slong n, gr_ctx_t ctx)

For a vector context (`gr_ctx_init_vector_gr_vec` and `gr_ctx_init_vector_space_gr_vec`): resizes `res` to length $n$ if the context allows all lengths, and otherwise returns `GR_DOMAIN` unless $n$ is the fixed length. Here `ctx` is the vector context.
:::

::: function
int _gr_str_is_space(char c)
const char * _gr_str_find_sep(const char * s, const char * end, int * err)
int _gr_str_strip_brackets(const char ** content_begin, const char ** content_end, const char * s, char open, char close)
int _gr_str_blank(const char * s, const char * end)
void _gr_str_trim(const char ** s, const char ** end)

String helpers shared by `gr_vec_set_str` and `gr_mat_set_str`: test for a whitespace character; find the next top-level comma between `s` and `end`, tracking the nesting of `()`, `[]` and `{}` and setting `*err` on unbalanced input; find the content between a leading `open` and its matching `close` bracket; test whether a range is blank; and trim whitespace from both ends of a range.
:::
