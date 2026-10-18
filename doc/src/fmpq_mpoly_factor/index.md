---
title: fmpq_mpoly_factor.h
header: fmpq_mpoly_factor.h
description: Factorisation of multivariate polynomials over the rational numbers, by reduction to the integers.
category: Multivariate polynomials
order: 6
---

This module factors polynomials in $\Q[x_0, \ldots, x_{n-1}]$ (see
[](fmpq_mpoly)). Since an `fmpq_mpoly_t` is stored as a rational content
times an integer polynomial, factoring over $\Q$ amounts to factoring the
integer part with [](fmpz_mpoly_factor) and multiplying the rational
content into the constant. The factorisation over $\Q$ is unique up to
multiplying the bases by nonzero rationals; the two normalising functions
below fix this freedom.

## Types {#types}

::: type
fmpq_mpoly_factor_struct
fmpq_mpoly_factor_t

A factored polynomial $c \prod_{i < \mathit{num}} p_i^{e_i}$ with a
rational constant: the fields are `constant` (an `fmpq_t`), `poly` (the
bases, of type `fmpq_mpoly_struct`), `exp` (the exponents, of type
`fmpz`), `num` and `alloc`.
:::

## Memory management {#memory}

::: function
void fmpq_mpoly_factor_init(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

Initialises `f` as the empty product with constant 1.
:::

::: function
void fmpq_mpoly_factor_realloc(fmpq_mpoly_factor_t f, slong alloc, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_factor_fit_length(fmpq_mpoly_factor_t f, slong len, const fmpq_mpoly_ctx_t ctx)

Change the allocation to exactly `alloc` factors (freeing the factors
beyond it, and everything if $\mathit{alloc} \le 0$), respectively ensure
room for at least `len` factors.
:::

::: function
void fmpq_mpoly_factor_clear(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

Frees the memory used by `f`.
:::

## Basic manipulation {#basic}

::: function
slong fmpq_mpoly_factor_length(const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

Returns the number of bases of `f`.
:::

::: function
void fmpq_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

Sets $c$ to the constant of `f`.
:::

::: function
void fmpq_mpoly_factor_get_base(fmpq_mpoly_t p, const fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_factor_swap_base(fmpq_mpoly_t p, fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)

Set $p$ to, respectively swap $p$ with, the base of index $i$.
:::

::: function
slong fmpq_mpoly_factor_get_exp_si(fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)

Returns the exponent of the base of index $i$, which must fit an `slong`.
:::

::: function
void fmpq_mpoly_factor_sort(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

Sorts the factors by increasing exponent, and factors with equal exponents
by `fmpq_mpoly_cmp`.
:::

::: function
int fmpq_mpoly_factor_expand(fmpq_mpoly_t A, const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to the product represented by `f`. Returns 0 if an exponent is
negative or a power cannot be computed, and 1 otherwise.
:::

::: function
int fmpq_mpoly_factor_make_monic(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_factor_make_integral(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)

Rescale the bases and adjust the constant so that the product is
unchanged: `make_monic` makes every base monic (leading coefficient 1 in
the monomial ordering), and `make_integral` moves the rational content of
every base into the constant, leaving the integer part (primitive, with
positive leading coefficient). Return 0 if a base is zero or a power of a
content cannot be computed, and 1 otherwise.
:::

## Factorisation {#factorisation}

The factoring functions return 1 for success and 0 for failure, in which
case `f` is unspecified; failure happens only in the circumstances
described for [](fmpz_mpoly_factor). They overwrite `f`. The
normalisation of the bases is not part of the specification: currently
they are integral, primitive and with positive leading coefficients, as
after `fmpq_mpoly_factor_make_integral`, but code that depends on a
normalisation should call one of the normalising functions.

::: function
int fmpq_mpoly_factor_squarefree(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Sets `f` to a factorisation of $A$ whose bases are squarefree and
pairwise coprime, computed with `fmpz_mpoly_factor_squarefree` on the
integer part of $A$. One multiplicity can occur in several bases; to get
one base for each multiplicity, sort with `fmpq_mpoly_factor_sort` and
multiply the bases with equal exponents.
:::

::: function
int fmpq_mpoly_factor(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Sets `f` to the factorisation of $A$ into a rational constant and
distinct irreducible polynomials, computed with `fmpz_mpoly_factor` on the
integer part of $A$. See [](fmpz_mpoly_factor/factor) for the methods.
:::

## Internal functions {#internal}

::: function
void _fmpq_mpoly_factor_swap_fmpz_mpoly_factor(fmpq_mpoly_factor_t f, fmpz_mpoly_factor_t g, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)

Sets `f` to $c$ times the factorisation $g$ over the integers, moving the
bases out of $g$ (which is left with unspecified bases). $g$ must be over
the integer context `ctx->zctx`.
:::
