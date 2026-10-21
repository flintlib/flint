---
title: long_extras.h
header: long_extras.h
description: Functions for signed words: digit counts, overflow-checked arithmetic, sign of a 2x2 determinant, random values and the Kronecker symbol.
category: Word-size integers
order: 2
---

This module collects a few functions on signed words (`slong`) that
complement the unsigned word functions of [](ulong_extras). The prefix
`z_` stands for a signed integer, as in the GMP conventions.

## Properties {#properties}

::: function
size_t z_sizeinbase(slong n, int b)

Returns the number of digits of $|n|$ in base $b \ge 2$, the sign not
counted. For $n = 0$ it returns 1. Unlike GMP's `mpz_sizeinbase`, the
result is exact for every base.
:::

## Checked arithmetic {#checked}

::: function
int z_add_checked(slong * a, slong b, slong c)
int z_mul_checked(slong * a, slong b, slong c)

Set `*a` to $b + c$ or $b c$ and return 0 if the result fits in an
`slong`; otherwise return 1, and `*a` holds the exact result reduced
modulo $2^{\bits}$ (into the range of `slong`). With GCC and Clang they
compile to the overflow-checking built-ins.
:::

::: function
int z_mat22_det_is_negative(slong m11, slong m12, slong m21, slong m22)

Returns 1 if the determinant $m_{11} m_{22} - m_{12} m_{21}$ is negative
and 0 otherwise. The determinant is computed exactly in two words, so
there is no overflow for any inputs.
:::

## Random functions {#random}

::: function
slong z_randtest(flint_rand_t state)

Returns a random value for testing. With probability $1/8$ it is one of
the special values $0$, $1$, $-1$, `COEFF_MAX`, `COEFF_MIN`, `WORD_MAX`
and `WORD_MIN`, chosen uniformly; otherwise it is a nonnegative random
value of between 0 and $\bits - 1$ bits, produced by `n_randbits`.
:::

::: function
slong z_randtest_not_zero(flint_rand_t state)

As `z_randtest`, but never returns 0.
:::

::: function
slong z_randint(flint_rand_t state, ulong limit)

Returns a random value $z$ with $|z| <$ `limit` and a random sign (so 0
is returned twice as often as other values). If `limit` is zero or larger
than `WORD_MAX`, it is replaced by `WORD_MAX`. The absolute value is
`n_randlimb(state) % limit`, which is close to uniform.
:::

## Number theory {#number-theory}

::: function
int z_kronecker(slong a, slong n)

Returns the Kronecker symbol $\left(\frac{a}{n}\right)$ for any $a$ and
$n$, including negative and even $n$ and $n = 0$ (where the symbol is 1
if $a = \pm 1$ and 0 otherwise). For odd positive $n$ this is the Jacobi
symbol, which `n_jacobi` also computes.
:::
