---
title: d_vec.h
header: d_vec.h
description: Vectors of double-precision floating-point numbers
category: Floating-point arithmetic
order: 7
---

This module provides basic operations on vectors of `double`, given as
a pointer and a length. All arithmetic uses the default rounding of the
processor (round to nearest). Besides plain loops it offers two dot
products that return an error estimate: a cheap heuristic one and the
compensated "dot product in thrice the working precision" of Ogita,
Rump and Oishi [@ORO2005], whose result is as accurate as if computed
with about 159-bit arithmetic and then rounded. For vectors of
double-word or quad-double numbers, see the `_dN_vec_` functions of
`dfloat`.

Unless stated otherwise, output vectors may alias input vectors.

## Memory management {#memory}

::: function
double * _d_vec_init(slong len)

Returns a newly allocated vector of `len` doubles. The entries are not
initialised.
:::

::: function
void _d_vec_clear(double * vec)

Frees a vector allocated by `_d_vec_init`.
:::

## Random generation {#random}

::: function
void _d_vec_randtest(double * f, flint_rand_t state, slong len, slong minexp, slong maxexp)

Sets the `len` entries of `f` to random values from
`d_randtest_signed`: signed numbers with exponents between `minexp` and
`maxexp`, or zero.
:::

## Assignment {#assignment}

::: function
void _d_vec_set(double * vec1, const double * vec2, slong len2)

Copies the `len2` entries of `vec2` to `vec1`.
:::

::: function
void _d_vec_zero(double * vec, slong len)

Sets the `len` entries of `vec` to zero.
:::

## Comparison {#comparison}

::: function
int _d_vec_equal(const double * vec1, const double * vec2, slong len)

Returns 1 if the entries of the two vectors are equal as doubles
(so that $+0 = -0$, and a NaN entry makes the result 0 unless the
pointers coincide), and 0 otherwise.
:::

::: function
int _d_vec_approx_equal(const double * vec1, const double * vec2, slong len, double eps)

Returns 1 if $|\mathit{vec1}[i] - \mathit{vec2}[i]| \le \varepsilon$
for every $i$ (an absolute tolerance, with the difference rounded to
nearest), and 0 otherwise.
:::

::: function
int _d_vec_is_zero(const double * vec, slong len)

Returns 1 if every entry is zero, and 0 otherwise.
:::

::: function
int _d_vec_is_approx_zero(const double * vec, slong len, double eps)

Returns 1 if $|\mathit{vec}[i]| \le \varepsilon$ for every $i$, and 0
otherwise.
:::

## Arithmetic {#arithmetic}

::: function
void _d_vec_add(double * res, const double * vec1, const double * vec2, slong len2)
void _d_vec_sub(double * res, const double * vec1, const double * vec2, slong len2)

Set `res` to the entrywise sum or difference of the two vectors of
length `len2`.
:::

::: function
void _d_vec_mul_2exp(double * res, const double * x, slong len, int e)

Sets `res` to the vector `x` multiplied by $2^e$, which is exact unless
an entry overflows or becomes subnormal (see `d_mul_2exp`).
:::

## Dot products and norms {#dot}

::: function
double _d_vec_dot(const double * vec1, const double * vec2, slong len2)

Returns the dot product of the two vectors, accumulated from left to
right in double precision. The error is bounded by
$\gamma_{n} \sum_i |x_i y_i|$ with $\gamma_n = n u/(1 - n u)$,
$u = 2^{-53}$ and $n$ the length, barring overflow and underflow
[@MBD2018].
:::

::: function
double _d_vec_norm(const double * vec, slong len)

Returns the square $\sum_i x_i^2$ of the Euclidean norm (not the norm
itself), accumulated in double precision.
:::

::: function
double _d_vec_dot_heuristic(const double * vec1, const double * vec2, slong len2, double * err)

Returns the dot product computed by summing the nonnegative and the
negative products separately and subtracting the two sums at the end.
If `err` is not `NULL`, sets `*err` to $2 n \max(2^{e_p - 53}, 2^{e_q - 53})$,
where $2^{e_p}$ and $2^{e_q}$ bound the two partial sums and $n$ is the
length. This is an estimate of the error, not a rigorous bound: it
does not account for overflow or underflow, and the factor $2n$ is
heuristic for long vectors.
:::

::: function
double _d_vec_dot_thrice(const double * vec1, const double * vec2, slong len2, double * err)

Returns the dot product computed as if in three times the working
precision and then rounded to double, using error-free transformations
(the algorithm `Dot3` of Ogita, Rump and Oishi [@ORO2005]). The result
is accurate to nearly full double precision unless the dot product is
extremely ill-conditioned. If `err` is not `NULL`, sets `*err` to the
error bound of [@ORO2005] in terms of the returned
value $r$, evaluated in double precision:
$(\varepsilon + 2\gamma^2) |r| + \gamma^3 \|x\|_2 \|y\|_2$ with
$\varepsilon = 2^{-52}$ and $\gamma = (4n-2)\varepsilon / (1 - (4n-2)\varepsilon)$.

::: details
Each product $x_i y_i$ is split into $h_i + e_i$ exactly by Dekker's
product with Veltkamp splitting (the code does not use FMA), the sum of
the $h_i$ is accumulated with Knuth's TwoSum, keeping the $n - 1$
rounding errors, and the resulting $2n$ terms are summed again with a
cascade of TwoSums followed by an ordinary sum of the errors; see
[](algorithms/floating-point/double-word#eft) for the error-free
transformations. The products must not overflow, and products below
about $2^{-969}$ lose the exactness of the splitting. A temporary
vector of $2n$ doubles is allocated.
:::
:::
