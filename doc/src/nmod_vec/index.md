---
title: nmod_vec.h
header: nmod_vec.h
description: Vectors over $\Z/n\Z$ for word-size $n$, stored as plain arrays of `ulong`.
category: Integers mod n
order: 2
---

This module provides arithmetic on vectors over $\Z/n\Z$ for a word-size
modulus $n$. A vector is a plain array of `ulong` together with its length,
and the modulus is passed as an `nmod_t` set up with `nmod_init`. The
functions are the building blocks of `nmod_poly`, `nmod_mat` and the other
modules over $\Z/n\Z$. Dot products, the most heavily used operation, choose
among several strategies depending on the length, the size of $n$ and the
instruction set; see [](nmod_vec/dot). Multiplication by a scalar and
inversion of all entries are explained in [](nmod_vec/scalar) and
[](nmod_vec/invert).

## Conventions {#conventions}

- Unless stated otherwise, input entries must be reduced, that is, lie in
  $[0, n)$, and output entries are reduced. This is not checked.
- Lengths are passed separately. A vector of length 0 is valid.
- Element-wise functions allow the output to be the same array as an input,
  but not a partially overlapping one.

## Memory management {#memory}

::: function
nn_ptr _nmod_vec_init(slong len)

Returns an array of `len` limbs allocated with `flint_malloc`. The entries
are not initialised.
:::

::: function
void _nmod_vec_clear(nn_ptr vec)

Frees an array allocated with `_nmod_vec_init`.
:::

## Random vectors {#random}

::: function
void _nmod_vec_rand(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)

Sets the entries of `vec` to integers chosen uniformly at random in
$[0, n)$.
:::

::: function
void _nmod_vec_randtest(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)

Sets the entries of `vec` to values of `n_randtest` reduced modulo $n$.
With probability 1/2 the vector is made sparse. Intended for test code.
:::

## Basic manipulation and comparison {#basic}

::: function
void _nmod_vec_set(nn_ptr res, nn_srcptr vec, slong len)

Copies the entries of `vec` to `res`.
:::

::: function
void _nmod_vec_zero(nn_ptr vec, slong len)

Sets the entries of `vec` to zero.
:::

::: function
void _nmod_vec_swap(nn_ptr a, nn_ptr b, slong length)

Swaps the entries of `a` and `b` one by one.
:::

::: function
int _nmod_vec_equal(nn_srcptr vec, nn_srcptr vec2, slong len)
int _nmod_vec_is_zero(nn_srcptr vec, slong len)

`_nmod_vec_equal` returns 1 if the two vectors have the same entries and 0
otherwise. `_nmod_vec_is_zero` returns 1 if every entry is zero and 0
otherwise.
:::

::: macro
NMOD_VEC_NORM(vec, i)

Decreases the integer variable `i` while `vec[i - 1]` is zero, so that `i`
becomes the length of `vec` with trailing zeros removed. Used to normalise
polynomials.
:::

::: function
flint_bitcnt_t _nmod_vec_max_bits(nn_srcptr vec, slong len)

Returns the largest bit length of an entry of `vec`, or 0 if all entries
are zero. The entries need not be reduced.
:::

::: function
void _nmod_vec_reduce(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)

Sets each entry of `res` to the corresponding entry of `vec` reduced
modulo $n$. The entries of `vec` may be any words.
:::

## Printing {#printing}

::: function
void _nmod_vec_print_pretty(nn_srcptr vec, slong len, nmod_t mod)
int _nmod_vec_fprint_pretty(FILE * file, nn_srcptr vec, slong len, nmod_t mod)

Prints a header line, then the entries in brackets, each padded to the
number of decimal digits of $n$. `_nmod_vec_fprint_pretty` returns a
positive value on success and a nonpositive value on a write error. Like
all functions with a `FILE` argument, it is declared only if `<stdio.h>` is
included before the FLINT headers.

::: details
For example, `_nmod_vec_print_pretty` prints a vector of length 5 modulo 197
as

```
<length-5 integer vector mod 197>
[ 33 181 107   6   0]
```

For a vector of length 0 only the header line is printed.
:::
:::

::: function
int _nmod_vec_print(nn_srcptr vec, slong len, nmod_t mod)
int _nmod_vec_fprint(FILE * f, nn_srcptr vec, slong len, nmod_t mod)

The same as `_nmod_vec_fprint_pretty`, writing to `stdout` and `f`
respectively.
:::

## Addition and negation {#addition}

::: function
void _nmod_vec_add(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
void _nmod_vec_sub(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)

Sets `res` to the entry-wise sum or difference of `vec1` and `vec2`.
:::

::: function
void _nmod_vec_neg(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)

Sets `res` to the entry-wise negation of `vec`.
:::

## Scalar multiplication {#scalar}

::: function
void _nmod_vec_scalar_mul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
void _nmod_vec_scalar_mul_nmod_generic(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
void _nmod_vec_scalar_mul_nmod_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
void _nmod_vec_scalar_mul_nmod_redc(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)

Sets `res` to $c$ times `vec`. All require $c < n$. The first function
chooses among the others. `_shoup` requires $n < 2^{\bits - 1}$ and
`_redc` requires $n$ odd.

::: details
`_nmod_vec_scalar_mul_nmod` requires reduced entries, since it copies or
negates `vec` when $c = 1$ or $c = n - 1$. The `_generic`, `_shoup` and
`_redc` variants also accept entries of `vec` that are not reduced.
:::
:::

::: function
void _nmod_vec_scalar_addmul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
void _nmod_vec_scalar_addmul_nmod_generic(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
void _nmod_vec_scalar_addmul_nmod_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)

Adds $c$ times `vec` to `res`. All require $c < n$. The first function
chooses among the others. `_shoup` requires $n < 2^{\bits - 1}$.

::: details
The entries of `res` must always be reduced. `_nmod_vec_scalar_addmul_nmod`
and `_nmod_vec_scalar_addmul_nmod_generic` also need reduced entries of
`vec`; the former adds or subtracts `vec` when $c = 1$ or $c = n - 1$.
`_nmod_vec_scalar_addmul_nmod_shoup` accepts entries of `vec` that are not
reduced.
:::
:::

## Inversion {#inversion}

::: function
void _nmod_vec_invert(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
void _nmod_vec_invert_generic(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
void _nmod_vec_invert_shoup(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)
void _nmod_vec_invert_naive(nn_ptr res, nn_srcptr vec, ulong len, nmod_t mod)

Sets each entry of `res` to the inverse of the corresponding entry of
`vec`. Raises a `FLINT_IMPINV` exception if an entry is not invertible.
The first three use a single modular inversion. `_shoup` requires
$n < 2^{\bits - 1}$. `_naive` inverts each entry separately. All four
allow `res` to be the same array as `vec`.

::: summary algorithms/arithmetic/batch-inversion :::
:::

## Dot products {#dot}

A dot product $\sum_i a_i b_i \bmod n$ is computed by choosing a method
once with `_nmod_vec_dot_params` and then calling `_nmod_vec_dot` or one of
its variants any number of times. The sum is reduced modulo $n$ only once,
at the end.

::: summary algorithms/linear-algebra/dot-products :::

::: type
dot_params_t

The method chosen by `_nmod_vec_dot_params`, defined in `nmod_types.h`.
Its fields `method` (of the enumeration type `dot_method_t`) and
`pow2_precomp` are internal.
:::

::: function
dot_params_t _nmod_vec_dot_params(ulong len, nmod_t mod)
ulong _nmod_vec_dot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
ulong _nmod_vec_dot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
ulong _nmod_vec_dot_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod, dot_params_t params)
ulong _nmod_vec_dot_strided(nn_srcptr vec1, slong stride1, nn_srcptr vec2, slong stride2, slong len, nmod_t mod, dot_params_t params)

`_nmod_vec_dot_params` returns parameters for dot products of length at
most `len` modulo $n$. The others return
$\sum_{i=0}^{\mathit{len}-1} a_i b_i \bmod n$, where $a_i$ is `vec1[i]`
(`vec1[i * stride1]` for the strided version) and $b_i$ is `vec2[i]`,
`vec2[len - 1 - i]`, `vec2[i][offset]` or `vec2[i * stride2]`
respectively. `params` must come from `_nmod_vec_dot_params` with the same
modulus and a length at least `len`.
:::

::: macro
NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, params)

Sets `res` to $\sum_{i=0}^{\mathit{len}-1} e_1(i) \, e_2(i) \bmod n$, where
`expr1` and `expr2` are expressions in the loop variable `i` that evaluate
to reduced values. `params` is as for `_nmod_vec_dot`.
:::

::: function
int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod)
int _nmod_vec_dot_bound_limbs_from_params(slong len, nmod_t mod, dot_params_t params)

Returns the number of limbs (0 to 3) of $\mathit{len} \, (n - 1)^2$, the
largest possible unreduced dot product of length `len`. The second version
reads the answer off `params`, which must have been computed for this
`len` and `mod`.
:::

## Dot product kernels {#dot-kernels}

These are the methods behind `_nmod_vec_dot`, `_nmod_vec_dot_rev` and
`_nmod_vec_dot_ptr`. Each kernel comes in the same three forms and returns
the same value. It takes the same arguments, except that `params` is
omitted or, for some kernels, replaced by `pow2_precomp`. The kernels are
normally reached through those functions, which handle lengths up to 11
themselves. A direct caller must respect the conditions given for each
kernel, where the unreduced sum $\sum_i a_i b_i$ is at most
$\mathit{len} \, (n - 1)^2$.

::: function
ulong _nmod_vec_dot1(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot1_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot1_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod)

Dot product kernel that sums in one limb. Requires the unreduced sum
$\sum_i a_i b_i$ to be less than $2^{\bits}$, and on 64-bit machines with
AVX2 also $n \le 2^{32}$.
:::

::: function
ulong _nmod_vec_dot_pow2(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot_pow2_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot_pow2_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod)

Dot product kernel for $n$ a power of two: sums in one limb modulo
$2^{\bits}$, which $n$ divides. No condition on the length.
:::

::: function
ulong _nmod_vec_dot2_split(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
ulong _nmod_vec_dot2_split_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
ulong _nmod_vec_dot2_split_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod, ulong pow2_precomp)

Dot product kernel for 64-bit machines that keeps the unreduced sum as
$h \, 2^{56} + \ell$ in two words and folds $h$ in with `pow2_precomp`
$= 2^{56} \bmod n$. Requires $n \le 1515531528$ and a bound on the length;
see
[splitting at bit 56](algorithms/linear-algebra/dot-products#split).
:::

::: function
ulong _nmod_vec_dot2_half(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, ulong pow2_precomp)
ulong _nmod_vec_dot2_half_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot2_half_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod)

Dot product kernel for $n \le 2^{\bits/2}$ that adds one-limb products
into a two-limb sum. On
64-bit machines with AVX2, `_nmod_vec_dot2_half` instead splits each
product at bit 56, needs `pow2_precomp` $= 2^{56} \bmod n$ and has a length
limit; see [splitting each product](algorithms/linear-algebra/dot-products#half).
:::

::: function
ulong _nmod_vec_dot2(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot2_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot2_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod)

Dot product kernel that adds two-limb products into a two-limb sum.
Requires the unreduced sum $\sum_i a_i b_i$ to be less than $2^{2 \bits}$.
:::

::: function
ulong _nmod_vec_dot3_acc(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot3_acc_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot3_acc_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod)

Dot product kernel that adds blocks of 8 products in two limbs and the block sums in three limbs.
Requires $8 (n - 1)^2 < 2^{2 \bits}$, that is
$n \le 6521908912666391107 \approx 2^{62.5}$ on 64-bit machines.
:::

::: function
ulong _nmod_vec_dot3(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot3_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot3_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod)

Dot product kernel that adds each product into a three-limb sum. No
conditions.
:::

::: function
ulong _nmod_vec_dot_u52(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot_u52_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot_u52_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod)

AVX512-IFMA dot product kernel for $n \le 2^{52}$. Requires the unreduced
sum $\sum_i a_i b_i$ to be less than $2^{2 \bits}$. Calls
the corresponding `_nmod_vec_dot2` function when IFMA is not available or
the length is below `FLINT_NMOD_VEC_DOT_U52_MIN_LEN`.
:::

::: function
ulong _nmod_vec_dot_u64(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot_u64_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot_u64_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod)

AVX512-IFMA dot product kernel for any $n$ and any length. Calls the corresponding
`_nmod_vec_dot3` function when IFMA is not available.
:::

::: function
ulong _nmod_vec_dot_split_limbs(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot_split_limbs_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
ulong _nmod_vec_dot_split_limbs_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod)

AVX2 or AVX-512 dot product kernel for any $n$ and any length, on 64-bit machines.
Calls the corresponding `_nmod_vec_dot3` function when AVX2 is not
available.
:::

## Arithmetic without reduction {#nored}

These functions add $c$ times `vec` to `res` without any modular
reduction, with `res` holding entries of one, two or three limbs. They
serve delayed reduction in Gaussian elimination and polynomial
division. The caller must ensure that nothing overflows.

::: function
void _nmod_vec_nored_scalar_addmul_halflimb(nn_ptr res, nn_srcptr vec, slong len, ulong c)
void _nmod_vec_nored_scalar_addmul_halflimb_avx2(nn_ptr res, nn_srcptr vec, slong len, ulong c)

Sets $r_i \gets r_i + c \, v_i$ for $0 \le i < \mathit{len}$ without
reduction, where $v_i$ is `vec[i]` and $r_i$ is the single limb `res[i]`.
Requires $c, v_i < 2^{\bits/2}$ and $r_i + c v_i < 2^{\bits}$. The AVX2 version is only declared on 64-bit
machines with AVX2, where the first function uses it from length 16 on.
:::

::: function
void _nmod_vec_nored_ll_scalar_addmul_halflimb(nn_ptr res, nn_srcptr vec, slong len, ulong c)
void _nmod_vec_nored_ll_scalar_addmul(nn_ptr res, nn_srcptr vec, slong len, ulong c)

Sets $r_i \gets r_i + c \, v_i$ for $0 \le i < \mathit{len}$ without
reduction, where $v_i$ is `vec[i]` and $r_i$ is the two-limb integer
stored least significant limb first at `res + 2 i`. Requires $r_i + c v_i < 2^{2 \bits}$. The `_halflimb` version
also requires $c v_i < 2^{\bits}$.
:::

::: function
void _nmod_vec_nored_lll_scalar_addmul(nn_ptr res, nn_srcptr vec, slong len, ulong c)

Sets $r_i \gets r_i + c \, v_i$ for $0 \le i < \mathit{len}$ without
reduction, where $v_i$ is `vec[i]` and $r_i$ is the three-limb integer
stored least significant limb first at `res + 3 i`. Requires $r_i + c v_i < 2^{3 \bits}$.
:::
