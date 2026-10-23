---
title: fq_default_poly_factor.h
header: fq_default_poly_factor.h
description: Factorisation, irreducibility testing and root finding for polynomials over any finite field $\mathbb{F}_q$, in the representation chosen by `fq_default`.
category: Finite fields
order: 4
---

This module factors polynomials of [](fq_default_poly) into irreducible
factors, tests irreducibility and squarefreeness, and finds roots, for
any finite field given by an `fq_default_ctx_t`. Every function calls the
function of the representation of the context: those of
[](fq_poly_factor) (and its `fq_nmod` and `fq_zech` versions) for
extension fields, and those of `nmod_poly_factor` and
`fmpz_mod_poly_factor` for prime fields. The guide
[](fq_default/getting-started) has an example.

::: summary algorithms/polynomials/factoring-finite-fields :::

## Factorisation structure {#type}

::: type
fq_default_poly_factor_struct
fq_default_poly_factor_t

A factorisation: a union of the factorisation types of the five
representations, of which the member matching the context is used. The
factors and exponents are read with `fq_default_poly_factor_length`,
`fq_default_poly_factor_get_poly` and `fq_default_poly_factor_exp`.
:::

::: function
void fq_default_poly_factor_init(fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
void fq_default_poly_factor_clear(fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)

Initialise `fac` as an empty factorisation, respectively free its memory.
:::

::: function
void fq_default_poly_factor_realloc(fq_default_poly_factor_t fac, slong alloc, const fq_default_ctx_t ctx)
void fq_default_poly_factor_fit_length(fq_default_poly_factor_t fac, slong len, const fq_default_ctx_t ctx)

Change the allocation to `alloc` factors, respectively ensure room for at
least `len` factors.
:::

::: function
slong fq_default_poly_factor_length(fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
slong fq_default_poly_factor_exp(fq_default_poly_factor_t fac, slong i, const fq_default_ctx_t ctx)

Return the number of factors, respectively the exponent of the factor of
index $i$.
:::

::: function
void fq_default_poly_factor_get_poly(fq_default_poly_t poly, const fq_default_poly_factor_t fac, slong i, const fq_default_ctx_t ctx)

Sets `poly` to the factor of index $i$.
:::

::: function
void fq_default_poly_factor_set(fq_default_poly_factor_t res, const fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)

Sets `res` to a copy of `fac`.
:::

::: function
void fq_default_poly_factor_insert(fq_default_poly_factor_t fac, const fq_default_poly_t poly, slong exp, const fq_default_ctx_t ctx)
void fq_default_poly_factor_concat(fq_default_poly_factor_t res, const fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
void fq_default_poly_factor_pow(fq_default_poly_factor_t fac, slong exp, const fq_default_ctx_t ctx)

Multiply the factorisation by `poly^exp` (adding to the exponent if
`poly` is already a factor), multiply it by all factors of `fac`, and
multiply all exponents by `exp`, as for `fq_poly_factor`.
:::

::: function
void fq_default_poly_factor_print(const fq_default_poly_factor_t fac, const fq_default_ctx_t ctx)
void fq_default_poly_factor_print_pretty(const fq_default_poly_factor_t fac, const char * var, const fq_default_ctx_t ctx)

Print the factors and exponents, one factor per line, in the plain,
respectively readable, format (with the variable name `var`).
:::

## Squarefree and irreducible polynomials {#tests}

::: function
int fq_default_poly_is_squarefree(const fq_default_poly_t f, const fq_default_ctx_t ctx)

Returns whether $f$ is squarefree.
:::

::: function
void fq_default_poly_factor_squarefree(fq_default_poly_factor_t res, const fq_default_poly_t f, const fq_default_ctx_t ctx)

Sets `res` to a squarefree factorisation of $f$: pairwise coprime monic
squarefree factors with their multiplicities, whose product with the
leading coefficient of $f$ is $f$.
:::

::: function
int fq_default_poly_is_irreducible(const fq_default_poly_t f, const fq_default_ctx_t ctx)

Returns whether $f$ is irreducible. Whether constants are considered
irreducible depends on the representation; `fq_poly_is_irreducible`
returns 1 for all polynomials of length at most 2.
:::

## Factorisation {#factorisation}

::: function
void fq_default_poly_factor(fq_default_poly_factor_t result, fq_default_t leading_coeff, const fq_default_poly_t input, const fq_default_ctx_t ctx)

Factors `input` as `leading_coeff` times the product of the monic
irreducible factors in `result` raised to their exponents. The algorithm
is chosen by the underlying module (see `fq_poly_factor`).
:::

::: function
void fq_default_poly_factor_distinct_deg(fq_default_poly_factor_t res, const fq_default_poly_t poly, slong * const * degs, const fq_default_ctx_t ctx)

Distinct-degree factorisation of the squarefree polynomial `poly`, as
`fq_poly_factor_distinct_deg`: appends the products $g_i$ of the
irreducible factors of degree $d_i$ to `res` and writes the $d_i$ to
`*degs`, which must have room for $\deg(\mathit{poly})$ entries.
:::

::: function
void fq_default_poly_factor_equal_deg(fq_default_poly_factor_t factors, const fq_default_poly_t pol, slong d, const fq_default_ctx_t ctx)

Appends to `factors` the irreducible factors of `pol`, a monic squarefree
product of irreducible polynomials of degree $d$.
:::

::: function
void fq_default_poly_factor_split_single(fq_default_poly_t linfactor, const fq_default_poly_t input, const fq_default_ctx_t ctx)

Sets `linfactor` to a monic linear factor of `input`, which must split
into linear factors. Only for extension fields: for the types
`FQ_DEFAULT_NMOD` and `FQ_DEFAULT_FMPZ_MOD` it throws `FLINT_ERROR`.
:::

## Roots {#roots}

::: function
void fq_default_poly_roots(fq_default_poly_factor_t r, const fq_default_poly_t f, int with_multiplicity, const fq_default_ctx_t ctx)

Sets `r` to the roots of the nonzero polynomial $f$ in the field, stored
as monic linear factors $x - a$, with their multiplicities as exponents
if `with_multiplicity` is nonzero, and exponent 1 otherwise.
:::
