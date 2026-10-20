---
title: fq_zech_mpoly.h
header: fq_zech_mpoly.h
description: Sparse multivariate polynomials over small finite fields with coefficients in Zech logarithm representation; a partial port of `fq_nmod_mpoly` used mainly by the factoring code of `nmod_mpoly`.
category: Multivariate polynomials
order: 14
---

This module implements multivariate polynomials over a finite field
$\mathbb{F}_q$ whose elements are represented by Zech logarithms (see
[](fq_zech)), which makes multiplication in the field cheap but limits $q$
to fields small enough for the logarithm tables, which have $q$ entries. The representation and the conventions are those of
`fq_nmod_mpoly` (see [](fq_nmod_mpoly) and [](nmod_mpoly#conventions)):
sorted terms with packed exponent vectors, and a context
`fq_zech_mpoly_ctx_t` with the number of variables, the monomial ordering
and the field. Unlike `fq_nmod_mpoly`, the coefficients are stored as an
array of `fq_zech_struct`.

The module exists mainly to support the `medprime` factoring algorithms
of `nmod_mpoly_factor`, which map a polynomial over $\mathbb{F}_p$ to
$\mathbb{F}_{p^k}$ with $p^k \le 10^6$ and factor it there
([](algorithms/multivariate/modp-factoring)). Its interface is therefore
incomplete: arithmetic, division, evaluation in one variable,
derivatives, degrees and GCDs are implemented (with the algorithms of
[](algorithms/multivariate/multiplication),
[](algorithms/multivariate/division) and
[](algorithms/multivariate/gcd)), while parsing, most term access
functions, random generation, composition and deflation, for example,
are not available. GCDs and factorisations are computed by
converting to `fq_nmod_mpoly`.

Unless stated otherwise, each function has the semantics of the
`fq_nmod_mpoly` function with the same name, with `fq_nmod` replaced by
`fq_zech`.

## Types {#types}

::: type
fq_zech_mpoly_struct
fq_zech_mpoly_t

A multivariate polynomial over $\mathbb{F}_q$, with the fields `coeffs`
(an array of `fq_zech_struct`), `exps`, `alloc` (the number of terms
allocated, for both arrays), `length` and `bits`.
:::

::: type
fq_zech_mpoly_ctx_struct
fq_zech_mpoly_ctx_t

The context of a polynomial ring $\mathbb{F}_q[x_0, \ldots, x_{m-1}]$,
with the fields `minfo` (an `mpoly_ctx_t`) and `fqctx` (an
`fq_zech_ctx_t`).
:::

::: type
fq_zech_mpoly_univar_struct
fq_zech_mpoly_univar_t

A polynomial in a main variable with `fq_zech_mpoly_t` coefficients and
`fmpz` exponents, as `fq_nmod_mpoly_univar_t`.
:::

## Context objects {#context}

::: function
void fq_zech_mpoly_ctx_init_deg(fq_zech_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, ulong p, slong deg)

Initialises `ctx` for polynomials in `nvars` variables with the ordering
`ord` over $\mathbb{F}_{p^{\mathit{deg}}}$, defined by a Conway polynomial
if the database has one and a random primitive polynomial otherwise
(see `fq_zech_ctx_init_ui`); the generator is printed as `#`. This is the
only way to create a context.
:::

::: function
void fq_zech_mpoly_ctx_change_modulus(fq_zech_mpoly_ctx_t ctx, slong deg)

Replaces the field of `ctx` by $\mathbb{F}_{p^{\mathit{deg}}}$ with the
same characteristic. Polynomials of the context become invalid.
:::

::: function
slong fq_zech_mpoly_ctx_nvars(const fq_zech_mpoly_ctx_t ctx)
ordering_t fq_zech_mpoly_ctx_ord(const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_ctx_clear(fq_zech_mpoly_ctx_t ctx)

Return the number of variables and the ordering, and free the context.
:::

::: function
nmod_t fq_zech_ctx_mod(const fq_zech_ctx_t ctx)

Returns the `nmod_t` of the characteristic $p$ of the field.
:::

## Memory management {#memory}

::: function
void fq_zech_mpoly_init(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_init2(fq_zech_mpoly_t A, slong alloc, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_init3(fq_zech_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_clear(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)

Initialise `A` as zero (allocating nothing, space for `alloc` terms, or
space for `alloc` terms with `bits`-bit exponent fields), and free it.
:::

::: function
void fq_zech_mpoly_fit_length(fq_zech_mpoly_t A, slong length, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_fit_bits(fq_zech_mpoly_t A, slong bits, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_fit_length_reset_bits(fq_zech_mpoly_t A, slong len, flint_bitcnt_t bits, const fq_zech_mpoly_ctx_t ctx)

Ensure space for `length` terms; ensure exponent fields of at least `bits`
bits, repacking the exponents; ensure space for `len` terms with exactly
`bits` bits, discarding the exponent data.
:::

::: function
void fq_zech_mpoly_truncate(fq_zech_mpoly_t A, slong newlen, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_repack_bits_inplace(fq_zech_mpoly_t A, flint_bitcnt_t Abits, const fq_zech_mpoly_ctx_t ctx)

Keep only the `newlen` leading terms, and repack the exponents into
`Abits` bits (returning 0 if they do not fit).
:::

## Input and output {#io}

::: function
int fq_zech_mpoly_fprint_pretty(FILE * file, const fq_zech_mpoly_t A, const char ** x, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_print_pretty(const fq_zech_mpoly_t A, const char ** x, const fq_zech_mpoly_ctx_t ctx)

Print `A` to `file` or `stdout` with the variable names `x` (or `x1`,
`x2`, ... if `x` is `NULL`).
:::

## Basic manipulation and constants {#basic}

::: function
void fq_zech_mpoly_gen(fq_zech_mpoly_t A, slong var, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_set(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_swap(fq_zech_mpoly_t A, fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_equal(const fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)

The variable of index `var`, assignment, swap and equality test.
:::

::: function
int fq_zech_mpoly_is_fq_zech(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_set_fq_zech(fq_zech_mpoly_t A, const fq_zech_t c, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_set_ui(fq_zech_mpoly_t A, ulong c, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_zero(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_one(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_is_zero(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_is_one(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)

Test whether `A` is constant, set it to the constant $c$ (an integer is
reduced modulo $p$), to 0 or to 1, and test for 0 and 1.
:::

::: function
fq_zech_struct * fq_zech_mpoly_leadcoeff(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)

Returns a pointer to the leading coefficient of the nonzero `A`.
:::

## Degrees and coefficients {#degrees}

::: function
int fq_zech_mpoly_degrees_fit_si(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_degrees_si(slong * degs, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_degrees_fmpz(fmpz ** degs, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
slong fq_zech_mpoly_degree_si(const fq_zech_mpoly_t A, slong var, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_degree_fmpz(fmpz_t deg, const fq_zech_mpoly_t A, slong var, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_total_degree_fits_si(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
slong fq_zech_mpoly_total_degree_si(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_total_degree_fmpz(fmpz_t td, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)

Degrees in each variable, in one variable, and the total degree; $-1$ for
the zero polynomial.
:::

::: function
void fq_zech_mpoly_get_coeff_vars_ui(fq_zech_mpoly_t C, const fq_zech_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fq_zech_mpoly_ctx_t ctx)

Sets `C` to the coefficient of $\prod_i x_{v_i}^{e_i}$ in `A` viewed as a
polynomial in the variables `vars[i]` with coefficients in the others.
:::

## Container operations {#container}

::: function
slong fq_zech_mpoly_length(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_is_canonical(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_assert_canonical(const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)

Return the number of terms, whether `A` is canonical, or throw if it is
not.
:::

::: function
int fq_zech_mpoly_term_exp_fits_ui(const fq_zech_mpoly_t A, slong i, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_term_exp_fits_si(const fq_zech_mpoly_t A, slong i, const fq_zech_mpoly_ctx_t ctx)

Return whether the exponents of the term of index $i$ fit a `ulong` or
`slong`.
:::

::: function
void fq_zech_mpoly_sort_terms(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_combine_like_terms(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)

Sort the terms into the order of the context, and combine equal adjacent
monomials and remove zero coefficients.
:::

## Arithmetic {#arithmetic}

::: function
void fq_zech_mpoly_add(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_t C, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_sub(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_t C, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_sub_fq_zech(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_t C, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_neg(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)

Set `A` to $B + C$, $B - C$, $B - c$ or $-B$.
:::

::: function
void fq_zech_mpoly_scalar_mul_fq_zech(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_t c, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_scalar_addmul_fq_zech(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_t C, const fq_zech_t d, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_make_monic(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)

Set `A` to $c B$, to $B + d C$, or to `B` divided by its leading
coefficient (throwing if `B` is zero).
:::

::: function
void fq_zech_mpoly_derivative(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, slong var, const fq_zech_mpoly_ctx_t ctx)

Sets `A` to the derivative of `B` with respect to the variable of index
`var`.
:::

::: function
void fq_zech_mpoly_evaluate_one_fq_zech(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, slong var, const fq_zech_t val, const fq_zech_mpoly_ctx_t ctx)

Sets `A` to `B` with the variable of index `var` replaced by `val`.
:::

::: function
void fq_zech_mpoly_mul(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_t C, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_mul_johnson(fq_zech_mpoly_t poly1, const fq_zech_mpoly_t poly2, const fq_zech_mpoly_t poly3, const fq_zech_mpoly_ctx_t ctx)

Set `A` to $B C$ by Johnson's heap algorithm [@Joh1974].
:::

::: function
int fq_zech_mpoly_pow_ui(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, ulong k, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_pow_rmul(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, ulong k, const fq_zech_mpoly_ctx_t ctx)

Set `A` to $B^k$; `fq_zech_mpoly_pow_ui` returns 1 for success.
:::

::: function
int fq_zech_mpoly_divides(fq_zech_mpoly_t Q, const fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_divides_monagan_pearce(fq_zech_mpoly_t poly1, const fq_zech_mpoly_t poly2, const fq_zech_mpoly_t poly3, const fq_zech_mpoly_ctx_t ctx)

If `B` divides `A`, set `Q` (`poly1`) to the quotient and return 1;
otherwise return 0. Uses the heap algorithm of Monagan and Pearce
[@MonPea2011].
:::

::: function
void fq_zech_mpoly_divrem(fq_zech_mpoly_t Q, fq_zech_mpoly_t R, const fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_divrem_monagan_pearce(fq_zech_mpoly_t q, fq_zech_mpoly_t r, const fq_zech_mpoly_t poly2, const fq_zech_mpoly_t poly3, const fq_zech_mpoly_ctx_t ctx)

Set `Q` and `R` to the quotient and remainder of `A` divided by `B`.
:::

## Greatest common divisors {#gcd}

::: function
int fq_zech_mpoly_gcd(fq_zech_mpoly_t G, const fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_gcd_cofactors(fq_zech_mpoly_t G, fq_zech_mpoly_t Abar, fq_zech_mpoly_t Bbar, const fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctx)

Set `G` to the monic GCD of `A` and `B` (and the cofactors $A/G$, $B/G$),
computed by converting to `fq_nmod_mpoly` and calling `fq_nmod_mpoly_gcd`
or `fq_nmod_mpoly_gcd_cofactors`. Return 1 for success.
:::

## Univariate representation {#univar}

::: function
void fq_zech_mpoly_univar_init(fq_zech_mpoly_univar_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_univar_clear(fq_zech_mpoly_univar_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_univar_fit_length(fq_zech_mpoly_univar_t A, slong length, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_univar_swap(fq_zech_mpoly_univar_t A, fq_zech_mpoly_univar_t B, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_univar_print_pretty(const fq_zech_mpoly_univar_t A, const char ** x, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_univar_assert_canonical(fq_zech_mpoly_univar_t A, const fq_zech_mpoly_ctx_t ctx)

Memory management, printing and consistency check of the univariate
representation.
:::

::: function
void fq_zech_mpoly_to_univar(fq_zech_mpoly_univar_t A, const fq_zech_mpoly_t B, slong var, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_from_univar(fq_zech_mpoly_t A, const fq_zech_mpoly_univar_t B, slong var, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_from_univar_bits(fq_zech_mpoly_t A, flint_bitcnt_t Abits, const fq_zech_mpoly_univar_t B, slong var, const fq_zech_mpoly_ctx_t ctx)

Convert `B` to a polynomial in the variable of index `var`, and back
(with `Abits`-bit exponents in the last version).
:::

::: function
int fq_zech_mpoly_univar_degree_fits_si(const fq_zech_mpoly_univar_t A, const fq_zech_mpoly_ctx_t ctx)
slong fq_zech_mpoly_univar_length(const fq_zech_mpoly_univar_t A, const fq_zech_mpoly_ctx_t ctx)
slong fq_zech_mpoly_univar_get_term_exp_si(fq_zech_mpoly_univar_t A, slong i, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_univar_get_term_coeff(fq_zech_mpoly_t c, const fq_zech_mpoly_univar_t A, slong i, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_univar_swap_term_coeff(fq_zech_mpoly_t c, fq_zech_mpoly_univar_t A, slong i, const fq_zech_mpoly_ctx_t ctx)

Access the terms of the univariate representation.
:::

## Internal functions {#internal}

::: function
void _fq_zech_mpoly_fit_length(fq_zech_struct ** coeff, ulong ** exps, slong * alloc, slong len, slong N, const fq_zech_ctx_t fqctx)
void _fq_zech_mpoly_set_length(fq_zech_mpoly_t A, slong newlen, const fq_zech_mpoly_ctx_t ctx)

Ensure room for `len` coefficients and exponent vectors of `N` words (the
common allocation is `alloc`), and set the length.
:::

::: function
void _fq_zech_mpoly_radix_sort1(fq_zech_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask)
void _fq_zech_mpoly_radix_sort(fq_zech_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask)
slong _fq_zech_mpoly_add(fq_zech_struct * coeff1, ulong * exp1, fq_zech_struct * coeff2, const ulong * exp2, slong len2, fq_zech_struct * coeff3, const ulong * exp3, slong len3, slong N, const ulong * cmpmask, const fq_zech_ctx_t fqctx)
slong _fq_zech_mpoly_mul_johnson(fq_zech_struct ** coeff1, ulong ** exp1, slong * alloc, const fq_zech_struct * coeff2, const ulong * exp2, slong len2, const fq_zech_struct * coeff3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_zech_ctx_t fqctx)
void _fq_zech_mpoly_eval_all_fq_zech(fq_zech_t eval, const fq_zech_struct * Acoeffs, const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, fq_zech_struct * const * alphas, const mpoly_ctx_t mctx, const fq_zech_ctx_t fqctx)

Radix sort of terms, merge-addition and heap multiplication of term
arrays, and evaluation at a point.
:::

::: function
int _fq_zech_mpoly_get_nmod_mpoly(nmod_mpoly_t s, const nmod_mpoly_ctx_t sctx, const fq_zech_mpoly_t t, const fq_zech_mpoly_ctx_t tctx)
void _fq_zech_mpoly_set_nmod_mpoly(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t Actx, const nmod_mpoly_t B, const nmod_mpoly_ctx_t Bctx)

`_fq_zech_mpoly_get_nmod_mpoly` sets `s` to `t` if all coefficients of
`t` lie in the prime field $\mathbb{F}_p$ and returns 1, and returns 0
otherwise. `_fq_zech_mpoly_set_nmod_mpoly` maps a polynomial over
$\mathbb{F}_p$ into $\mathbb{F}_q$. The two contexts must have the same
variables and ordering.
:::

::: function
void _fq_zech_mpoly_get_fq_nmod_mpoly(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctxA, const fq_zech_mpoly_t B, const fq_zech_mpoly_ctx_t ctxB)
void _fq_zech_mpoly_set_fq_nmod_mpoly(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctxA, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctxB)

Convert between Zech logarithm and polynomial representation of the same
field (the `fq_nmod_ctx` of the Zech context) with the same variables and
ordering.
:::

::: function
void fq_zech_mpolyl_lead_coeff(fq_zech_mpoly_t c, const fq_zech_mpoly_t A, slong num_vars, const fq_zech_mpoly_ctx_t ctx)

Sets `c` to the leading coefficient of `A` viewed as a polynomial in its
first `num_vars` variables.
:::

::: type
fq_zech_mpolyu_struct
fq_zech_mpolyu_t
fq_zech_mpolyn_struct
fq_zech_mpolyn_t
fq_zech_mpolyun_struct
fq_zech_mpolyun_t
fq_zech_mpoly_geobucket_struct
fq_zech_mpoly_geobucket_t

Recursive representations and geobuckets, as for `fq_nmod_mpoly`. Only
`fq_zech_mpolyu_t` has functions.
:::

::: function
void fq_zech_mpolyu_init(fq_zech_mpolyu_t A, flint_bitcnt_t bits, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpolyu_clear(fq_zech_mpolyu_t A, const fq_zech_mpoly_ctx_t uctx)
void fq_zech_mpolyu_swap(fq_zech_mpolyu_t A, fq_zech_mpolyu_t B)
void fq_zech_mpolyu_zero(fq_zech_mpolyu_t A, const fq_zech_mpoly_ctx_t uctx)
void fq_zech_mpolyu_one(fq_zech_mpolyu_t A, const fq_zech_mpoly_ctx_t uctx)
int fq_zech_mpolyu_is_one(fq_zech_mpolyu_t A, const fq_zech_mpoly_ctx_t uctx)
int fq_zech_mpolyu_is_canonical(const fq_zech_mpolyu_t poly, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpolyu_fit_length(fq_zech_mpolyu_t A, slong length, const fq_zech_mpoly_ctx_t uctx)
void fq_zech_mpolyu_print_pretty(const fq_zech_mpolyu_t poly, const char ** x, const fq_zech_mpoly_ctx_t ctx)

Memory management and basic operations of `fq_zech_mpolyu_t`.
:::
