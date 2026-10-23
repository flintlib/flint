---
title: fmpz_mod_mpoly_q.h
header: fmpz_mod_mpoly_q.h
description: Multivariate rational functions over a prime field $\mathbb{F}_p$, as reduced fractions of `fmpz_mod_mpoly` polynomials with monic denominators.
category: Multivariate polynomials
order: 11
---

An `fmpz_mod_mpoly_q_t` represents an element of the rational function
field $\mathbb{F}_p(x_0, \ldots, x_{m-1})$ as a fraction $N/D$ of two
`fmpz_mod_mpoly_t`. It is the analogue over $\mathbb{F}_p$ of
`fmpz_mpoly_q` (see [](fmpz_mpoly_q)). The polynomial context
`fmpz_mod_mpoly_ctx_t` (see [](fmpz_mod_mpoly)) fixes the prime $p$, the
number of variables and the monomial ordering; the user must make sure
that $p$ is prime, otherwise the GCD computations may fail or throw.

A fraction is in *canonical form* if $\gcd(N, D) = 1$ and $D$ is monic
(in particular nonzero); zero is $0/1$. All functions assume canonical
inputs and produce canonical outputs, which is not checked; values built
directly through `fmpz_mod_mpoly_q_numref` and `fmpz_mod_mpoly_q_denref`
must be passed through `fmpz_mod_mpoly_q_canonicalise`. Arithmetic keeps
fractions reduced with as few GCD computations as possible, in the manner
of `fmpq` arithmetic: for example a sum $a/b + c/d$ first computes
$g = \gcd(b, d)$ and only needs to cancel common factors with $g$.

## Types and macros {#types}

::: type
fmpz_mod_mpoly_q_struct
fmpz_mod_mpoly_q_t

A pair of `fmpz_mod_mpoly_struct`, the numerator `num` and the
denominator `den`. An `fmpz_mod_mpoly_q_t` is an array of length one of
`fmpz_mod_mpoly_q_struct`.
:::

::: macro
fmpz_mod_mpoly_q_numref(x)
fmpz_mod_mpoly_q_denref(x)

Pointers to the numerator and the denominator of `x`, usable as
`fmpz_mod_mpoly_t`.
:::

## Memory management {#memory}

::: function
void fmpz_mod_mpoly_q_init(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_clear(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)

Initialise `res` to zero, and free it.
:::

## Assignment {#assignment}

::: function
void fmpz_mod_mpoly_q_swap(fmpz_mod_mpoly_q_t x, fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)

Swaps `x` and `y` efficiently.
:::

::: function
void fmpz_mod_mpoly_q_set(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_set_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_t x, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_set_si(fmpz_mod_mpoly_q_t res, slong x, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_q_set_fmpq(fmpz_mod_mpoly_q_t res, const fmpq_t x, const fmpz_mod_mpoly_ctx_t ctx)

Set `res` to `x` (integers are reduced modulo $p$). The `fmpq` version
returns 1 if the denominator of `x` is invertible modulo $p$ and 0
otherwise, in which case `res` is undefined.
:::

## Canonicalisation {#canonical}

::: function
void fmpz_mod_mpoly_q_canonicalise(fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

Puts `x` in canonical form by dividing numerator and denominator by
their GCD and making the denominator monic. The denominator must be
nonzero.
:::

::: function
int fmpz_mod_mpoly_q_is_canonical(const fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)

Returns whether `x` is in canonical form: the numerator and denominator
are canonical polynomials, the denominator is nonzero and monic, and the
GCD of numerator and denominator is 1. Used for testing.
:::

## Properties and special values {#properties}

::: function
int fmpz_mod_mpoly_q_is_zero(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_q_is_one(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_q_is_fmpz_mod(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

Return whether `x` is 0, 1, or a constant.
:::

::: function
void fmpz_mod_mpoly_q_zero(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_one(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_gen(fmpz_mod_mpoly_q_t res, slong i, const fmpz_mod_mpoly_ctx_t ctx)

Set `res` to 0, 1, or the variable of index $i$, $0 \le i < m$.
:::

## Input and output {#io}

The array `x` (or `vars`) gives the variable names, starting with variable
0; if it is `NULL`, the variables are named `x1`, `x2`, and so on.

::: function
void fmpz_mod_mpoly_q_print_pretty(const fmpz_mod_mpoly_q_t f, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
char * fmpz_mod_mpoly_q_get_str_pretty(const fmpz_mod_mpoly_q_t f, const char ** vars, const fmpz_mod_mpoly_ctx_t ctx)

Print `f` to standard output, or return it as a string that the caller
frees with `flint_free`, in the form `(N)/(D)` (or just `N` if $D = 1$).
:::

::: function
int fmpz_mod_mpoly_q_set_str_pretty(fmpz_mod_mpoly_q_t res, const char * s, const char ** vars, fmpz_mod_mpoly_ctx_t ctx)

Sets `res` to the rational function in the null-terminated string `s`,
written with integers, the variable names, parentheses, `+`, `-`, `*`, `/`
and `^` followed by an integer exponent. Returns 0 on success; on a parse
error or a division by zero, sets `res` to zero and returns $-1$.
:::

## Random generation {#random}

::: function
void fmpz_mod_mpoly_q_randtest(fmpz_mod_mpoly_q_t res, flint_rand_t state, slong length, slong exp_bound, const fmpz_mod_mpoly_ctx_t ctx)

Sets `res` to a random canonical fraction whose numerator and denominator
have at most `length` terms and exponents below `exp_bound`.
:::

## Comparison {#comparison}

::: function
int fmpz_mod_mpoly_q_equal(const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)

Returns whether `x` and `y` are equal.
:::

## Arithmetic {#arithmetic}

The functions that take an `fmpq` or divide by a scalar return 1 if the
denominator (resp. the divisor) is invertible modulo $p$ and 0 otherwise;
in the latter case `res` is undefined. The `fmpz_mod` versions require
the scalar to be reduced modulo $p$.

::: function
void fmpz_mod_mpoly_q_neg(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

Sets `res` to $-x$.
:::

::: function
void fmpz_mod_mpoly_q_add(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_q_add_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_add_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_add_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_add_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong c, const fmpz_mod_mpoly_ctx_t ctx)

Set `res` to $x + y$.
:::

::: function
void fmpz_mod_mpoly_q_sub(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_q_sub_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_sub_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_sub_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_sub_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong c, const fmpz_mod_mpoly_ctx_t ctx)

Set `res` to $x - y$.
:::

::: function
void fmpz_mod_mpoly_q_mul(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_q_mul_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_mul_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_mul_fmpz_mod(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_q_mul_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong c, const fmpz_mod_mpoly_ctx_t ctx)

Set `res` to $x y$. For $x = a/b$ and $y = c/d$, `fmpz_mod_mpoly_q_mul`
cancels $\gcd(a, d)$ and $\gcd(b, c)$ before multiplying.
:::

::: function
void fmpz_mod_mpoly_q_div(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_q_t y, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_q_div_fmpq(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpq_t y, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_q_div_fmpz(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_q_div_si(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, slong c, const fmpz_mod_mpoly_ctx_t ctx)

Set `res` to $x / y$. `fmpz_mod_mpoly_q_div` throws `FLINT_ERROR` if $y
= 0$.
:::

::: function
void fmpz_mod_mpoly_q_inv(fmpz_mod_mpoly_q_t res, const fmpz_mod_mpoly_q_t x, const fmpz_mod_mpoly_ctx_t ctx)

Sets `res` to $1/x$. Throws `FLINT_ERROR` if $x = 0$.
:::

## Internal functions {#internal}

::: function
void _fmpz_mod_mpoly_q_add(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den, const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den, const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_q_sub(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den, const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den, const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_q_mul(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den, const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den, const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_q_div(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den, const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den, const fmpz_mod_mpoly_t y_num, const fmpz_mod_mpoly_t y_den, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_q_add_fmpz_mod(fmpz_mod_mpoly_t res_num, fmpz_mod_mpoly_t res_den, const fmpz_mod_mpoly_t x_num, const fmpz_mod_mpoly_t x_den, const fmpz_t y, const fmpz_mod_mpoly_ctx_t ctx)

Arithmetic on canonical fractions given as separate numerators and
denominators; the outputs must not alias the inputs.
:::

::: function
void fmpz_mod_mpoly_gcd_assert_successful(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t x, const fmpz_mod_mpoly_t y, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_q_mpoly_divexact(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t x, const fmpz_mod_mpoly_t y, const fmpz_mod_mpoly_ctx_t ctx)

`fmpz_mod_mpoly_gcd_assert_successful` sets `res` to the GCD of `x` and
`y` and throws if `fmpz_mod_mpoly_gcd` fails. `_fmpz_mod_mpoly_q_mpoly_divexact`
sets `res` to $x/y$ for an exact division, dividing by a constant `y`
directly.
:::
