---
title: fq_nmod_mpoly.h
header: fq_nmod_mpoly.h
description: Sparse multivariate polynomials over finite fields $\mathbb{F}_q$ of word-size characteristic, with heap arithmetic, GCDs that fall back to extension fields for small $q$, square roots and resultants.
category: Multivariate polynomials
order: 12
---

This module implements multivariate polynomials over a finite field
$\mathbb{F}_q$, $q = p^d$, with a word-size prime $p$, the field being
given by an `fq_nmod_ctx_t` (see [](fq_nmod)). It has the representation
and the interface of `nmod_mpoly` (see [](nmod_mpoly), [](mpoly) and
[](algorithms/multivariate/representation)):
sorted terms with packed exponent vectors, and a context
`fq_nmod_mpoly_ctx_t` holding the number of variables, the monomial
ordering and the field. Each coefficient is stored inline as $d$ words,
the coordinates of the element on the basis $1, \theta, \ldots,
\theta^{d-1}$ of $\mathbb{F}_q = \mathbb{F}_p[\theta]/(f(\theta))$ (the
`n_fq` format), so no `fq_nmod_struct` is stored and no reference to a
coefficient as an `fq_nmod_t` is available; the coefficient functions copy
values in and out.

Multiplication uses Johnson's heap algorithm [@Joh1974]
([](algorithms/multivariate/multiplication)) and division the heap
algorithm of Monagan and Pearce [@MonPea2011]
([](algorithms/multivariate/division)). GCDs use Brown's, Zippel's and
Hensel-lifting algorithms ([](algorithms/multivariate/gcd),
[](algorithms/multivariate/hensel-lifting)) with evaluation points in
$\mathbb{F}_q$;
when $\mathbb{F}_q$ is too small, they work over extensions
$\mathbb{F}_{q^n}$ obtained by reducing one variable modulo an irreducible
polynomial over $\mathbb{F}_q$
([](algorithms/multivariate/modp-gcd)). Factorisation is in
`fq_nmod_mpoly_factor.h`. The module `fq_zech_mpoly` provides the same
polynomials with coefficients in Zech logarithm representation. The guide
[](fq_nmod_mpoly/getting-started) shows how to set up a context and work
with coefficients involving the generator of $\mathbb{F}_q$.

## Types and conventions {#conventions}

::: type
fq_nmod_mpoly_struct
fq_nmod_mpoly_t

A multivariate polynomial over $\mathbb{F}_q$, with the fields `coeffs`
($d$ words per coefficient), `exps`, `length`, `bits`, `coeffs_alloc` and
`exps_alloc` (allocations in words).
:::

::: type
fq_nmod_mpoly_ctx_struct
fq_nmod_mpoly_ctx_t

The context of a polynomial ring $\mathbb{F}_q[x_0, \ldots, x_{m-1}]$,
with the fields `minfo` (an `mpoly_ctx_t`: number of variables and
ordering) and `fqctx` (an `fq_nmod_ctx_t`, owned by the context).
:::

The conventions of [](nmod_mpoly#conventions) apply: all polynomials
passed to one function belong to the same context, variable 0 is the most
significant, inputs and outputs are canonical (sorted terms, nonzero
coefficients) except for the container operations, and outputs may alias
inputs unless stated otherwise. Since $\mathbb{F}_q$ is a field, division,
GCDs and factorisation have no extra preconditions.

## Context objects {#context}

::: function
void fq_nmod_mpoly_ctx_init(fq_nmod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, const fq_nmod_ctx_t fqctx)

Initialises `ctx` for polynomials in `nvars` variables with the ordering
`ord` (`ORD_LEX`, `ORD_DEGLEX` or `ORD_DEGREVLEX`) and coefficients in
the field of `fqctx`. The context makes its own copy of the field, so
`fqctx` may be cleared afterwards.
:::

::: function
void fq_nmod_mpoly_ctx_init_deg(fq_nmod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, ulong p, slong deg)

Initialises `ctx` with coefficients in $\mathbb{F}_{p^{\mathit{deg}}}$,
defined by the modulus chosen by `fq_nmod_ctx_init_ui` (a Conway
polynomial if the database has one, otherwise an irreducible polynomial of
minimal weight); the generator is printed as `#`.
:::

::: function
void fq_nmod_mpoly_ctx_init_rand(fq_nmod_mpoly_ctx_t ctx, flint_rand_t state, slong max_nvars, flint_bitcnt_t p_bits, slong deg_bound)

Initialises `ctx` with a random number of variables in $[0,$
`max_nvars`$]$, a random ordering, and a random field of characteristic a
prime of at most `p_bits` bits and degree at most `deg_bound`. Used for
testing.
:::

::: function
slong fq_nmod_mpoly_ctx_nvars(const fq_nmod_mpoly_ctx_t ctx)
ordering_t fq_nmod_mpoly_ctx_ord(const fq_nmod_mpoly_ctx_t ctx)

Return the number of variables and the ordering.
:::

::: function
void fq_nmod_mpoly_ctx_change_modulus(fq_nmod_mpoly_ctx_t ctx, slong deg)

Replaces the field of `ctx` by $\mathbb{F}_{p^{\mathit{deg}}}$ with the
same characteristic, as in `fq_nmod_mpoly_ctx_init_deg`. Polynomials of
the context become invalid. Used by the algorithms that try extension
fields of increasing degree.
:::

::: function
void fq_nmod_mpoly_ctx_clear(fq_nmod_mpoly_ctx_t ctx)

Frees the memory used by `ctx`.
:::

## Memory management {#memory}

::: function
void fq_nmod_mpoly_init(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_init2(fq_nmod_mpoly_t A, slong alloc, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_init3(fq_nmod_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)

Initialise `A` as the zero polynomial, allocating nothing, space for
`alloc` terms, or space for `alloc` terms with `bits` bits per exponent
field.
:::

::: function
void fq_nmod_mpoly_clear(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Frees the memory used by `A`.
:::

::: function
void fq_nmod_mpoly_fit_length(fq_nmod_mpoly_t A, slong length, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_realloc(fq_nmod_mpoly_t A, slong alloc, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_fit_length_fit_bits(fq_nmod_mpoly_t A, slong len, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_fit_length_reset_bits(fq_nmod_mpoly_t A, slong len, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)

Ensure space for `length` terms, change the allocation to `alloc` terms
(not less than the length), and ensure space with exponent fields of at
least (`fit_bits`, repacking) or exactly (`reset_bits`, discarding the
exponent data) `bits` bits.
:::

::: function
void fq_nmod_mpoly_truncate(fq_nmod_mpoly_t A, slong newlen, const fq_nmod_mpoly_ctx_t ctx)

Keeps only the `newlen` leading terms of `A` if it is longer.
:::

::: function
int fq_nmod_mpoly_repack_bits(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_repack_bits_inplace(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx)

Repack the exponents of `B` into `A` (resp. of `A` in place) with `Abits`
bits per field; return 0 if some exponent does not fit, setting `A` to
zero.
:::

## Input and output {#io}

The array `x` gives the variable names, starting with variable 0; if it is
`NULL`, the variables are named `x1`, `x2`, and so on. Coefficients are
printed as polynomials in the generator of the field, with the name stored
in the field context.

::: function
char * fq_nmod_mpoly_get_str_pretty(const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_fprint_pretty(FILE * file, const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_print_pretty(const fq_nmod_mpoly_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)

Return a string representing `A` (freed by the caller with
`flint_free`), or print it to `file` or `stdout`.
:::

::: function
int fq_nmod_mpoly_set_str_pretty(fq_nmod_mpoly_t A, const char * str, const char ** x, const fq_nmod_mpoly_ctx_t ctx)

Parses `str` with integers, the variables in `x`, `+`, `-`, `*`, exact
`/` and `^` followed by an integer exponent. The generator of the field is
not available as a symbol unless it is included in `x`. Returns 0 on
success, and $-1$ (setting `A` to zero) on failure.
:::

## Basic manipulation {#basic}

::: function
void fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_is_gen(const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to the variable of index `var`; return whether `A` is that
variable, or, for `var` $< 0$, any variable.
:::

::: function
void fq_nmod_mpoly_set(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_swap(fq_nmod_mpoly_t A, fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_equal(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Assignment, swap and equality test.
:::

## Constants {#constants}

::: function
int fq_nmod_mpoly_is_fq_nmod(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_get_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Return whether `A` is constant, and set `c` to the constant `A`
(throwing if it is not constant).
:::

::: function
void fq_nmod_mpoly_set_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_n_fq(fq_nmod_mpoly_t A, const ulong * c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_ui(fq_nmod_mpoly_t A, ulong c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_fmpz(fq_nmod_mpoly_t A, const fmpz_t c, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to the constant $c$, given as an `fq_nmod_t`, in `n_fq` format
($d$ words), or as an integer reduced modulo $p$.
:::

::: function
void fq_nmod_mpoly_set_fq_nmod_gen(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Sets `A` to the constant $\theta$, the generator of $\mathbb{F}_q$ over
$\mathbb{F}_p$ (as returned by `fq_nmod_gen`).
:::

::: function
void fq_nmod_mpoly_zero(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_one(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to 0 or 1.
:::

::: function
int fq_nmod_mpoly_equal_fq_nmod(const fq_nmod_mpoly_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_is_zero(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_is_one(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Return whether `A` equals the constant $c$, zero or one.
:::

::: function
ulong * fq_nmod_mpoly_get_nonzero_n_fq(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Returns a pointer to the coefficient ($d$ words) of `A`, which must be a
nonzero constant.
:::

## Degrees {#degrees}

The zero polynomial has degree $-1$ in every variable and total degree
$-1$.

::: function
int fq_nmod_mpoly_degrees_fit_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_degrees_si(slong * degs, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_degrees_fmpz(fmpz ** degs, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Set `degs` to the degrees of `A` in all variables, and test whether they
fit an `slong`.
:::

::: function
slong fq_nmod_mpoly_degree_si(const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_degree_fmpz(fmpz_t deg, const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)

The degree of `A` in the variable of index `var`.
:::

::: function
int fq_nmod_mpoly_total_degree_fits_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
slong fq_nmod_mpoly_total_degree_si(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_total_degree_fmpz(fmpz_t td, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

The total degree of `A`, and whether it fits an `slong`.
:::

::: function
void fq_nmod_mpoly_used_vars(int * used, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Sets `used[i]` to a nonzero value if the variable of index $i$ appears in
`A`, and to zero otherwise.
:::

## Coefficients {#coefficients}

::: function
void fq_nmod_mpoly_get_coeff_fq_nmod_monomial(fq_nmod_t c, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t M, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_coeff_fq_nmod_monomial(fq_nmod_mpoly_t A, const fq_nmod_t c, const fq_nmod_mpoly_t M, const fq_nmod_mpoly_ctx_t ctx)

Get or set the coefficient of the monomial `M` (throws if `M` is not a
monomial).
:::

::: function
void fq_nmod_mpoly_get_coeff_fq_nmod_fmpz(fq_nmod_t c, const fq_nmod_mpoly_t A, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_get_coeff_fq_nmod_ui(fq_nmod_t c, const fq_nmod_mpoly_t A, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A, const fq_nmod_t c, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_coeff_fq_nmod_ui(fq_nmod_mpoly_t A, const fq_nmod_t c, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx)

Get or set the coefficient of the monomial with exponent vector `exp`.
Setting a coefficient to zero removes the term. Each call costs time
linear in the length of `A`.
:::

::: function
void fq_nmod_mpoly_get_coeff_vars_ui(fq_nmod_mpoly_t C, const fq_nmod_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fq_nmod_mpoly_ctx_t ctx)

Sets `C` to the coefficient of $\prod_i x_{v_i}^{e_i}$ in `A` viewed as a
polynomial in the distinct variables $x_{v_i}$, $v_i$ = `vars[i]`, $e_i$ =
`exps[i]`, with coefficients in the other variables.
:::

::: function
int fq_nmod_mpoly_is_monic(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Returns whether `A` is nonzero with leading coefficient 1.
:::

## Conversion to univariate polynomials {#univariate-poly}

::: function
int fq_nmod_mpoly_is_fq_nmod_poly(const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_get_fq_nmod_poly(fq_nmod_poly_t A, const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_fq_nmod_poly(fq_nmod_mpoly_t A, const fq_nmod_poly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

Test whether `A` involves only the variable of index `var`; convert such a
polynomial to an `fq_nmod_poly` (returning 0 if it involves another
variable); and set `A` to the univariate polynomial `B` in that variable.
:::

## Comparison {#comparison}

::: function
int fq_nmod_mpoly_cmp(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Returns 1, $-1$ or 0 according as `A` is greater than, less than or equal
to `B` in a total order on polynomials that agrees with the monomial
ordering on monomials.
:::

## Container operations {#container}

These functions access terms by their index $i$, $0 \le i <$ length, and
can produce non-canonical polynomials, which `fq_nmod_mpoly_sort_terms`
and `fq_nmod_mpoly_combine_like_terms` repair. Out-of-range indices throw.

::: function
int fq_nmod_mpoly_is_canonical(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_assert_canonical(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Return whether `A` is canonical, or throw if it is not.
:::

::: function
slong fq_nmod_mpoly_length(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_resize(fq_nmod_mpoly_t A, slong new_length, const fq_nmod_mpoly_ctx_t ctx)

Return the number of terms, or change it to `new_length`, deleting terms
at the end or appending zero terms.
:::

::: function
void fq_nmod_mpoly_get_term_coeff_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_term_coeff_fq_nmod(fq_nmod_mpoly_t A, slong i, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)

Get or set the coefficient of the term of index $i$.
:::

::: function
int fq_nmod_mpoly_term_exp_fits_ui(const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_term_exp_fits_si(const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_get_term_exp_fmpz(fmpz ** exp, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_get_term_exp_ui(ulong * exp, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_get_term_exp_si(slong * exp, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
ulong fq_nmod_mpoly_get_term_var_exp_ui(const fq_nmod_mpoly_t A, slong i, slong var, const fq_nmod_mpoly_ctx_t ctx)
slong fq_nmod_mpoly_get_term_var_exp_si(const fq_nmod_mpoly_t A, slong i, slong var, const fq_nmod_mpoly_ctx_t ctx)

Test whether the exponents of the term of index $i$ fit a `ulong` or
`slong`, and get its exponent vector or the exponent of one variable; the
`ulong` and `slong` versions throw if a value does not fit.
:::

::: function
void fq_nmod_mpoly_set_term_exp_fmpz(fq_nmod_mpoly_t A, slong i, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_term_exp_ui(fq_nmod_mpoly_t A, slong i, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx)

Set the exponent vector of the term of index $i$.
:::

::: function
void fq_nmod_mpoly_get_term(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_get_term_monomial(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)

Set `M` to the term of index $i$, or to its monomial with coefficient 1.
:::

::: function
void fq_nmod_mpoly_push_term_fq_nmod_fmpz(fq_nmod_mpoly_t A, const fq_nmod_t c, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_push_term_fq_nmod_ffmpz(fq_nmod_mpoly_t A, const fq_nmod_t c, const fmpz * exp, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_push_term_fq_nmod_ui(fq_nmod_mpoly_t A, const fq_nmod_t c, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx)

Append a term with coefficient $c$ and exponent vector `exp` (pointers to
`fmpz`, an array of `fmpz`, or `ulong`) in constant amortised time.
:::

::: function
void fq_nmod_mpoly_sort_terms(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_combine_like_terms(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Sort the terms into the order of the context, and combine adjacent equal
monomials and remove zero coefficients.
:::

::: function
void fq_nmod_mpoly_reverse(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Sets `A` to `B` with its terms in reverse order.
:::

## Random generation {#random}

::: function
void fq_nmod_mpoly_randtest_bound(fq_nmod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_randtest_bounds(fq_nmod_mpoly_t A, flint_rand_t state, slong length, ulong * exp_bounds, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_randtest_bits(fq_nmod_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t exp_bits, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to a random polynomial with at most `length` terms and exponents
below `exp_bound`, below `exp_bounds[i]` for variable $i$, or of at most
`exp_bits` bits.
:::

## Addition and subtraction {#addition}

::: function
void fq_nmod_mpoly_add_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_add_n_fq(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const ulong * c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_sub_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_t C, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to $B + c$ or $B - c$.
:::

::: function
void fq_nmod_mpoly_add(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_sub(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to $B + C$ or $B - C$.
:::

## Scalar operations {#scalar}

::: function
void fq_nmod_mpoly_neg(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Sets `A` to $-B$.
:::

::: function
void fq_nmod_mpoly_scalar_mul_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_scalar_mul_n_fq(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const ulong * c, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to $c B$.
:::

::: function
void fq_nmod_mpoly_scalar_addmul_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_t e, const fq_nmod_mpoly_ctx_t ctx)

Sets `A` to $B + e C$.
:::

::: function
void fq_nmod_mpoly_make_monic(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Sets `A` to `B` divided by its leading coefficient; throws if `B` is zero.
:::

## Differentiation {#derivative}

::: function
void fq_nmod_mpoly_derivative(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

Sets `A` to the derivative of `B` with respect to the variable of index
`var`.
:::

## Evaluation and composition {#evaluation}

::: function
void fq_nmod_mpoly_evaluate_all_fq_nmod(fq_nmod_t ev, const fq_nmod_mpoly_t A, fq_nmod_struct * const * vals, const fq_nmod_mpoly_ctx_t ctx)

Sets `ev` to the value of `A` at the point `vals[0]`, `vals[1]`, ... (one
per variable).
:::

::: function
void fq_nmod_mpoly_evaluate_one_fq_nmod(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, slong var, const fq_nmod_t val, const fq_nmod_mpoly_ctx_t ctx)

Sets `A` to `B` with the variable of index `var` replaced by `val`.
:::

::: function
int fq_nmod_mpoly_compose_fq_nmod_poly(fq_nmod_poly_t A, const fq_nmod_mpoly_t B, fq_nmod_poly_struct * const * C, const fq_nmod_mpoly_ctx_t ctx)

Sets `A` to `B` evaluated at the univariate polynomials `C[i]`. Returns 1
for success and 0 if the result would be too large.
:::

::: function
int fq_nmod_mpoly_compose_fq_nmod_mpoly(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, fq_nmod_mpoly_struct * const * C, const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)
int fq_nmod_mpoly_compose_fq_nmod_mpoly_horner(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, fq_nmod_mpoly_struct * const * C, const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)
int fq_nmod_mpoly_compose_fq_nmod_mpoly_geobucket(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, fq_nmod_mpoly_struct * const * C, const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)

Set `A` to `B` evaluated at the polynomials `C[i]` of the context
`ctxAC`, one per variable of `ctxB`; `A` must not alias `B` or any
`C[i]`. Return 1 for success and 0 if the computation would be too large.
The main function uses a matrix transformation of the exponents when every
`C[i]` is zero or a monomial with coefficient 1, and otherwise a Horner
scheme or a geobucket sum.
:::

::: function
void fq_nmod_mpoly_compose_fq_nmod_mpoly_gen(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const slong * c, const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)

Sets `A` to `B` with the variable of index $i$ of `ctxB` replaced by the
variable of index `c[i]` of `ctxAC`, or by zero if `c[i]` is negative.
:::

## Multiplication {#multiplication}

::: function
void fq_nmod_mpoly_mul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t C, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_mul_johnson(fq_nmod_mpoly_t poly1, const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to $B C$ by Johnson's heap algorithm [@Joh1974].
:::

## Powering {#powering}

::: function
int fq_nmod_mpoly_pow_ui(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, ulong k, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_pow_fmpz(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fmpz_t k, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to $B^k$. `fq_nmod_mpoly_pow_fmpz` throws if $k < 0$ and
returns 0 if $k$ does not fit an `slong` and `B` has more than one term;
otherwise both return 1.
:::

::: function
void fq_nmod_mpoly_pow_rmul(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, ulong k, const fq_nmod_mpoly_ctx_t ctx)

Sets `A` to $B^k$ by repeated multiplication.
:::

## Division {#division}

::: function
int fq_nmod_mpoly_divides(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

If `B` divides `A`, sets `Q` to $A/B$ and returns 1; otherwise sets `Q`
to zero and returns 0. Throws `FLINT_DIVZERO` if `B` is zero.
:::

::: function
void fq_nmod_mpoly_divexact(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Sets `Q` to $A/B$, assuming that the division is exact.
:::

::: function
void fq_nmod_mpoly_div(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_divrem(fq_nmod_mpoly_t Q, fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Set `Q` and `R` to the quotient and remainder of `A` divided by `B`, with
no term of $R$ divisible by the leading monomial of $B$.
:::

::: function
void fq_nmod_mpoly_divrem_ideal(fq_nmod_mpoly_struct ** Q, fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A, fq_nmod_mpoly_struct * const * B, slong len, const fq_nmod_mpoly_ctx_t ctx)

Divides `A` by the `len` polynomials `B[i]`: $A = R + \sum_i Q_i B_i$
with no term of $R$ divisible by any leading monomial of a $B_i$.
:::

::: function
int fq_nmod_mpoly_divides_monagan_pearce(fq_nmod_mpoly_t poly1, const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_div_monagan_pearce(fq_nmod_mpoly_t q, const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_divrem_monagan_pearce(fq_nmod_mpoly_t q, fq_nmod_mpoly_t r, const fq_nmod_mpoly_t poly2, const fq_nmod_mpoly_t poly3, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_divrem_ideal_monagan_pearce(fq_nmod_mpoly_struct ** q, fq_nmod_mpoly_t r, const fq_nmod_mpoly_t poly2, fq_nmod_mpoly_struct * const * poly3, slong len, const fq_nmod_mpoly_ctx_t ctx)

The heap implementations [@MonPea2011] behind the division functions.
:::

## Square roots {#sqrt}

::: function
int fq_nmod_mpoly_sqrt(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_sqrt_heap(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_is_square(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

If `A` is a square, set `Q` to a square root and return 1; otherwise
return 0. In characteristic 2 the square root exists exactly when all
exponents are even, and it is obtained by halving the exponents and
taking the square roots $c^{2^{d-1}}$ of the coefficients.
:::

::: function
int fq_nmod_mpoly_quadratic_root(fq_nmod_mpoly_t Q, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

If $Q^2 + AQ = B$ has a solution, sets `Q` to one and returns 1;
otherwise returns 0. In characteristic 2 the equation is solved term by
term ([](algorithms/multivariate/modp-factoring#quadratic)).
:::

## Greatest common divisors {#gcd}

The GCD functions return 1 for success and 0 for failure, which only
happens when the exponents need more than one word per field even after
deflation.

::: function
void fq_nmod_mpoly_term_content(fq_nmod_mpoly_t M, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Sets `M` to the GCD of the terms of `A`: zero, or the monomial with
coefficient 1 whose exponents are the minima over the terms of `A`.
:::

::: function
int fq_nmod_mpoly_content_vars(fq_nmod_mpoly_t g, const fq_nmod_mpoly_t A, slong * vars, slong vars_length, const fq_nmod_mpoly_ctx_t ctx)

Sets `g` to the GCD of the coefficients of `A` viewed as a polynomial in
the variables `vars[i]`, and returns 1 for success.
:::

::: function
int fq_nmod_mpoly_gcd(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_gcd_cofactors(fq_nmod_mpoly_t G, fq_nmod_mpoly_t Abar, fq_nmod_mpoly_t Bbar, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Set `G` to the monic GCD of `A` and `B` (zero if both are zero), and
`fq_nmod_mpoly_gcd_cofactors` also $\bar{A} = A/G$ and $\bar{B} = B/G$.
The algorithm is chosen as for `nmod_mpoly_gcd`; when $\mathbb{F}_q$ is
too small, Brown's and Zippel's algorithms continue with images over
$\mathbb{F}_q[x]/(h) \cong \mathbb{F}_{q^n}$ for irreducible $h$ in one of
the variables.
:::

::: function
int fq_nmod_mpoly_gcd_brown(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_gcd_hensel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_gcd_zippel(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_gcd_zippel2(fq_nmod_mpoly_t G, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Set `G` to the monic GCD by one algorithm: Brown's dense interpolation,
Hensel lifting, Zippel's sparse interpolation, or sparse interpolation
with two main variables. Return 0 if the algorithm cannot be applied.
:::

::: function
void fq_nmod_mpoly_deflation(fmpz * shift, fmpz * stride, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_deflate(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fmpz * shift, const fmpz * stride, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_inflate(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fmpz * shift, const fmpz * stride, const fq_nmod_mpoly_ctx_t ctx)

As `nmod_mpoly_deflation`, `nmod_mpoly_deflate` and `nmod_mpoly_inflate`:
the shifts and strides of the exponents of each variable, and the
substitutions $e \mapsto (e - \mathit{shift})/\mathit{stride}$ and $e
\mapsto \mathit{shift} + \mathit{stride} \cdot e$.
:::

## Resultants and discriminants {#resultant}

::: function
int fq_nmod_mpoly_resultant(fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_discriminant(fq_nmod_mpoly_t R, const fq_nmod_mpoly_t A, slong var, const fq_nmod_mpoly_ctx_t ctx)

Set `R` to the resultant of `A` and `B`, or the discriminant of `A`, with
respect to the variable of index `var`, by the subresultant algorithm of
Ducos [@Duc2000]. Return 1 for success.
:::

## Univariate representation {#univar}

An `fq_nmod_mpoly_univar_t` holds a polynomial in a main variable with
`fq_nmod_mpoly_t` coefficients in the other variables, as in
[](nmod_mpoly#univar).

::: type
fq_nmod_mpoly_univar_struct
fq_nmod_mpoly_univar_t

The fields are `coeffs`, `exps` (`fmpz` exponents), `alloc` and
`length`.
:::

::: function
void fq_nmod_mpoly_univar_init(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_univar_clear(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_univar_fit_length(fq_nmod_mpoly_univar_t A, slong length, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_univar_swap(fq_nmod_mpoly_univar_t A, fq_nmod_mpoly_univar_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_univar_zero(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_univar_print_pretty(const fq_nmod_mpoly_univar_t A, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_univar_assert_canonical(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)

Memory management, printing and consistency check.
:::

::: function
void fq_nmod_mpoly_to_univar(fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_from_univar(fq_nmod_mpoly_t A, const fq_nmod_mpoly_univar_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

Convert `B` to a polynomial in the variable of index `var`, and back.
:::

::: function
int fq_nmod_mpoly_univar_degree_fits_si(const fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)
slong fq_nmod_mpoly_univar_length(const fq_nmod_mpoly_univar_t A, const fq_nmod_mpoly_ctx_t ctx)
slong fq_nmod_mpoly_univar_get_term_exp_si(fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_univar_get_term_coeff(fq_nmod_mpoly_t c, const fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_univar_swap_term_coeff(fq_nmod_mpoly_t c, fq_nmod_mpoly_univar_t A, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_univar_set_coeff_ui(fq_nmod_mpoly_univar_t A, ulong e, const fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)

Access the terms: whether the degree fits an `slong`, the number of terms,
the exponent and the coefficient of the term of index $i$, and setting the
coefficient of $x^e$.
:::

::: function
int fq_nmod_mpoly_univar_pseudo_gcd(fq_nmod_mpoly_univar_t Gx, const fq_nmod_mpoly_univar_t Ax, const fq_nmod_mpoly_univar_t Bx, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_univar_resultant(fq_nmod_mpoly_t R, const fq_nmod_mpoly_univar_t Ax, const fq_nmod_mpoly_univar_t Bx, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_univar_discriminant(fq_nmod_mpoly_t D, const fq_nmod_mpoly_univar_t Fx, const fq_nmod_mpoly_ctx_t ctx)

A GCD up to factors free of the main variable, the resultant and the
discriminant with respect to the main variable, by subresultants. Return
1 for success.
:::

## Internal functions {#internal}

These functions are used by the arithmetic, GCD and factorisation code of
this module and of `nmod_mpoly`. Their interfaces may change without
notice.

### Low-level helpers

::: function
void _fq_nmod_mpoly_fit_length(ulong ** coeffs, slong * coeffs_alloc, slong d, ulong ** exps, slong * exps_alloc, slong N, slong length)
void _fq_nmod_mpoly_set_length(fq_nmod_mpoly_t A, slong newlen, const fq_nmod_mpoly_ctx_t ctx)

Ensure room for `length` coefficients of `d` words and exponent vectors of
`N` words in the given arrays, and set the length field.
:::

::: function
ulong * _fq_nmod_mpoly_leadcoeff(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Returns a pointer to the leading coefficient ($d$ words) of the nonzero
`A`.
:::

::: function
void _fq_nmod_mpoly_radix_sort1(fq_nmod_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask, slong d)
void _fq_nmod_mpoly_radix_sort(fq_nmod_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask, slong d)
void _fq_nmod_mpoly_push_exp_ffmpz(fq_nmod_mpoly_t A, const fmpz * exp, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_push_exp_pfmpz(fq_nmod_mpoly_t A, fmpz * const * exp, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_push_exp_ui(fq_nmod_mpoly_t A, const ulong * exp, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_set_coeff_fq_nmod_fmpz(fq_nmod_mpoly_t A, const fq_nmod_t c, const fmpz * exp, const fq_nmod_mpoly_ctx_t ctx)

Radix sort of the terms `left` to `right - 1` (with `d` words per
coefficient), appending a term with undefined coefficient, and setting a
coefficient given an array of `fmpz` exponents.
:::

::: function
slong _fq_nmod_mpoly_add(ulong * coeff1, ulong * exp1, ulong * coeff2, const ulong * exp2, slong len2, ulong * coeff3, const ulong * exp3, slong len3, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx)
void _fq_nmod_mpoly_mul_johnson(fq_nmod_mpoly_t A, const ulong * Bcoeffs, const ulong * Bexps, slong Blen, const ulong * Ccoeffs, const ulong * Cexps, slong Clen, flint_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t ctx)
int _fq_nmod_mpoly_divides_monagan_pearce(fq_nmod_mpoly_t A, const ulong * coeff2, const ulong * exp2, slong len2, const ulong * coeff3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N, const ulong * cmpmask, const fq_nmod_ctx_t fqctx)
void _fq_nmod_mpoly_eval_all_fq_nmod(fq_nmod_t ev, const ulong * Acoeffs, const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, fq_nmod_struct * const * alphas, const mpoly_ctx_t mctx, const fq_nmod_ctx_t fqctx)
void _fq_nmod_mpoly_compose_mat(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fmpz_mat_t M, const fq_nmod_mpoly_ctx_t ctxB, const fq_nmod_mpoly_ctx_t ctxAC)

Kernels on term arrays: merge-addition, heap multiplication and exact
division (returning 0 if not exact), evaluation at a point, and
substitution of monomials by a matrix acting on exponent vectors.
:::

::: function
void _fq_nmod_mpoly_set_fq_nmod_poly(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_struct * Bcoeffs, slong Blen, slong var, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_to_fq_nmod_poly_deflate(fq_nmod_poly_t A, const fq_nmod_mpoly_t B, slong var, const ulong * Bshift, const ulong * Bstride, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_from_fq_nmod_poly_inflate(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_poly_t B, slong var, const ulong * Ashift, const ulong * Astride, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_from_univar(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpoly_univar_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

Conversions to and from univariate representations with given exponent
bits, shifts and strides.
:::

::: function
void mpoly_void_ring_init_fq_nmod_mpoly_ctx(mpoly_void_ring_t R, const fq_nmod_mpoly_ctx_t ctx)

Sets up the generic ring interface `R` used by the parser and the
subresultant code.
:::

::: function
void fq_nmod_mpoly_remainder_strongtest(const fq_nmod_mpoly_t r, const fq_nmod_mpoly_t g, const fq_nmod_mpoly_ctx_t ctx)

Throws if some term of `r` is divisible by the leading monomial of `g`.
:::

::: function
int fq_nmod_next(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx)
void fq_nmod_next_not_zero(fq_nmod_t alpha, const fq_nmod_ctx_t fqctx)

Step `alpha` through the elements of $\mathbb{F}_q$ in the lexicographic
order of their coordinates: `fq_nmod_next` returns 0 when it wraps around
to zero, `fq_nmod_next_not_zero` skips zero.
:::

::: type
fq_nmod_mpoly_geobucket_struct
fq_nmod_mpoly_geobucket_t

A geobucket [@Yan1998], as for `nmod_mpoly_geobucket_t`.
:::

::: function
void fq_nmod_mpoly_geobucket_init(fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_geobucket_clear(fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_geobucket_fit_length(fq_nmod_mpoly_geobucket_t B, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_geobucket_set(fq_nmod_mpoly_geobucket_t B, fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_geobucket_add(fq_nmod_mpoly_geobucket_t B, fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_geobucket_sub(fq_nmod_mpoly_geobucket_t B, fq_nmod_mpoly_t p, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_geobucket_empty(fq_nmod_mpoly_t p, fq_nmod_mpoly_geobucket_t B, const fq_nmod_mpoly_ctx_t ctx)

As the `nmod_mpoly_geobucket` functions: `set` swaps `p` into the
geobucket, `add` and `sub` add or subtract `p`, `empty` returns the sum in
`p`.
:::

### Embeddings into extension fields

Let $\mathbb{F}_q = \mathbb{F}_p[\theta]/(f(\theta))$ with $\deg f = m$,
and let $\mathbb{F}_{q^n} = \mathbb{F}_p[\varphi]/(g(\varphi))$ with $\deg
g = mn$. An embedding presents the large field also as
$\mathbb{F}_q[x]/(h(x))$ for an irreducible $h$ of degree $n$ over
$\mathbb{F}_q$, and stores the images of $\theta$ and $x$ in
$\mathbb{F}_p[\varphi]/(g)$ and of $\varphi$ in $\mathbb{F}_q[x]$, with
the matrices of the two $\mathbb{F}_p$-linear isomorphisms. There are $m$
such embeddings, one for each root of $f$ in the large field. They let
the GCD and factoring algorithms reduce a polynomial over $\mathbb{F}_q$
modulo $h$ in one variable and work over $\mathbb{F}_{q^n}$. The names
start with `bad_` because the construction is a simple, not particularly
fast one (see `fq_nmod_embed` for the general machinery).

::: type
bad_fq_nmod_embed_struct
bad_fq_nmod_embed_t
bad_fq_nmod_mpoly_embed_chooser_struct
bad_fq_nmod_mpoly_embed_chooser_t

An embedding, and a chooser that iterates over the embeddings for
extension degrees $n$ of increasing size.
:::

::: function
void bad_fq_nmod_embed_array_init(bad_fq_nmod_embed_struct * emb, const fq_nmod_ctx_t bigctx, const fq_nmod_ctx_t smallctx)
void bad_fq_nmod_embed_clear(bad_fq_nmod_embed_t emb)

Initialise the $m$ embeddings of the field of `smallctx` into the field of
`bigctx` (whose degree must be a multiple of $m$), and clear one.
:::

::: function
void bad_fq_nmod_embed_sm_to_lg(fq_nmod_t out, const fq_nmod_poly_t in, const bad_fq_nmod_embed_t emb)
void bad_fq_nmod_embed_lg_to_sm(fq_nmod_poly_t out, const fq_nmod_t in, const bad_fq_nmod_embed_t emb)
void bad_n_fq_embed_sm_to_lg(ulong * out_, const n_poly_t in_, const bad_fq_nmod_embed_t emb)
void bad_n_fq_embed_lg_to_sm(n_poly_t out_, const ulong * in_, const bad_fq_nmod_embed_t emb)
void bad_fq_nmod_embed_n_fq_sm_to_fq_nmod_lg(fq_nmod_t out, const n_poly_t in_, const bad_fq_nmod_embed_t emb)
void bad_fq_nmod_embed_fq_nmod_lg_to_n_fq_sm(n_poly_t out_, const fq_nmod_t in, const bad_fq_nmod_embed_t emb)

Map a polynomial in $x$ over $\mathbb{F}_q$ to its residue in the large
field (reduction modulo $h$), and an element of the large field to its
representative of degree less than $n$ in $x$, for the `fq_nmod` and
`n_fq` representations.
:::

::: function
void bad_fq_nmod_embed_sm_elem_to_lg(fq_nmod_t out, const fq_nmod_t in, const bad_fq_nmod_embed_t emb)
void bad_n_fq_embed_sm_elem_to_lg(ulong * out, const ulong * in, const bad_fq_nmod_embed_t emb)

Map an element of $\mathbb{F}_q$ into the large field.
:::

::: function
bad_fq_nmod_embed_struct * bad_fq_nmod_mpoly_embed_chooser_init(bad_fq_nmod_mpoly_embed_chooser_t embc, fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t randstate)
bad_fq_nmod_embed_struct * bad_fq_nmod_mpoly_embed_chooser_next(bad_fq_nmod_mpoly_embed_chooser_t embc, fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t randstate)
void bad_fq_nmod_mpoly_embed_chooser_clear(bad_fq_nmod_mpoly_embed_chooser_t embc, fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t randstate)

Start iterating over embeddings: `init` sets `ectx` to a context over an
extension of degree $n \ge 2$ (with $n m \log_2 p$ about 20) and returns
the first embedding; `next` returns the next one, moving to the next
degree (and changing `ectx`) after all $m$ embeddings of a degree have
been used; `clear` frees the chooser and `ectx`.
:::

### Recursive representations

::: type
fq_nmod_mpolyu_struct
fq_nmod_mpolyu_t
fq_nmod_mpolyn_struct
fq_nmod_mpolyn_t
fq_nmod_mpolyun_struct
fq_nmod_mpolyun_t

Sparse univariate polynomials with `fq_nmod_mpoly_t` coefficients (in a
lexicographic context with uniform exponent bits), multivariate
polynomials whose coefficients are univariate polynomials in the last
variable (stored as `n_fq_poly_struct`), and sparse univariate
polynomials with `fq_nmod_mpolyn_t` coefficients, as the corresponding
`nmod_mpoly` types.
:::

::: function
void fq_nmod_mpolyu_init(fq_nmod_mpolyu_t A, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_clear(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
void fq_nmod_mpolyu_swap(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t uctx)
void fq_nmod_mpolyu_zero(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
void fq_nmod_mpolyu_one(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
int fq_nmod_mpolyu_is_one(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx)
int fq_nmod_mpolyu_is_canonical(const fq_nmod_mpolyu_t poly, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_fit_length(fq_nmod_mpolyu_t A, slong length, const fq_nmod_mpoly_ctx_t uctx)
void fq_nmod_mpolyu_set(fq_nmod_mpolyu_t A, const fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t uctx)
void fq_nmod_mpolyu_print_pretty(const fq_nmod_mpolyu_t poly, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_degrees_si(slong * degs, const fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_repack_bits_inplace(fq_nmod_mpolyu_t A, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
fq_nmod_mpoly_struct * _fq_nmod_mpolyu_get_coeff(fq_nmod_mpolyu_t A, ulong pow, const fq_nmod_mpoly_ctx_t uctx)
ulong * fq_nmod_mpolyu_leadcoeff(const fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_shift_right(fq_nmod_mpolyu_t A, ulong s)
void fq_nmod_mpolyu_shift_left(fq_nmod_mpolyu_t A, ulong s)
void fq_nmod_mpolyu_scalar_mul_fq_nmod(fq_nmod_mpolyu_t A, fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyu_content_mpoly(fq_nmod_mpoly_t g, const fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_divexact_mpoly_inplace(fq_nmod_mpolyu_t A, fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_mul_mpoly(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_mul_mpoly_inplace(fq_nmod_mpolyu_t A, fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyuu_divides(fq_nmod_mpolyu_t Q, const fq_nmod_mpolyu_t A, const fq_nmod_mpolyu_t B, slong nmainvars, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_setform(fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx)

Operations on `fq_nmod_mpolyu_t`, as the `nmod_mpolyu` functions:
memory management, access to the coefficient of $x^{\mathit{pow}}$
(inserted if absent), the leading coefficient of the leading coefficient,
shifts of the main exponent, scaling, the GCD of the coefficients,
multiplication and exact division by a polynomial in the other variables,
exact division with `nmainvars` main variables, and copying the monomial
structure.
:::

::: function
void fq_nmod_mpoly_to_mpolyu_perm_deflate(fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t uctx, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride)
void fq_nmod_mpoly_from_mpolyu_perm_inflate(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx, const fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t uctx, const slong * perm, const ulong * shift, const ulong * stride)
void fq_nmod_mpoly_to_mpolyl_perm_deflate(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t lctx, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride)
void fq_nmod_mpoly_from_mpolyl_perm_inflate(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t lctx, const slong * perm, const ulong * shift, const ulong * stride)
void fq_nmod_mpoly_to_mpolyn_perm_deflate(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t nctx, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride)
void fq_nmod_mpoly_from_mpolyn_perm_inflate(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t ctx, const fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t nctx, const slong * perm, const ulong * shift, const ulong * stride)
void fq_nmod_mpoly_cvtto_mpolyn(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_cvtfrom_mpolyn(fq_nmod_mpoly_t A, fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_cvtto_mpolyun(fq_nmod_mpolyun_t A, const fq_nmod_mpolyu_t B, slong k, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyu_cvtfrom_mpolyun(fq_nmod_mpolyu_t A, fq_nmod_mpolyun_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

Conversions between the representations, with a permutation `perm` of the
variables, deflation by `shift` and `stride`, and the inverse operations
with `Abits`-bit exponents.
:::

::: function
void fq_nmod_mpolyl_lead_coeff(fq_nmod_mpoly_t c, const fq_nmod_mpoly_t A, slong num_vars, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyl_content(fq_nmod_mpoly_t g, const fq_nmod_mpoly_t A, slong num_vars, const fq_nmod_mpoly_ctx_t ctx)

For `A` viewed as a polynomial in its first `num_vars` variables: the
leading coefficient, and the GCD of the coefficients (returning 1 for
success).
:::

::: function
void fq_nmod_mpolyn_init(fq_nmod_mpolyn_t A, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_clear(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_swap(fq_nmod_mpolyn_t A, fq_nmod_mpolyn_t B)
void fq_nmod_mpolyn_zero(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_one(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyn_is_zero(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyn_is_canonical(const fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_fit_length(fq_nmod_mpolyn_t A, slong length, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_fit_bits(fq_nmod_mpolyn_t A, slong bits, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_set(fq_nmod_mpolyn_t A, const fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_print_pretty(const fq_nmod_mpolyn_t A, const char ** x_in, const fq_nmod_mpoly_ctx_t ctx)
ulong * fq_nmod_mpolyn_leadcoeff(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
n_poly_struct * fq_nmod_mpolyn_leadcoeff_poly(const fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
slong fq_nmod_mpolyn_lastdeg(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyn_is_nonzero_fq_nmod(const fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_scalar_mul_fq_nmod(fq_nmod_mpolyn_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_mul_poly(fq_nmod_mpolyn_t A, const fq_nmod_mpolyn_t B, const fq_nmod_poly_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_divexact_poly(fq_nmod_mpolyn_t A, const fq_nmod_mpolyn_t B, const fq_nmod_poly_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_content_poly(fq_nmod_poly_t g, fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ctx)

Operations on `fq_nmod_mpolyn_t`: memory management, the leading
coefficient (an element, or the univariate polynomial), the degree in the
last variable, the test for a nonzero constant, scaling, multiplication
and exact division by a univariate polynomial in the last variable, and
the GCD `g` of the univariate coefficients.
:::

::: function
void fq_nmod_mpolyun_init(fq_nmod_mpolyun_t A, flint_bitcnt_t bits, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_clear(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_swap(fq_nmod_mpolyun_t A, fq_nmod_mpolyun_t B)
void fq_nmod_mpolyun_zero(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_one(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyun_is_canonical(const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_fit_length(fq_nmod_mpolyun_t A, slong length, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_set(fq_nmod_mpolyun_t A, const fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_print_pretty(const fq_nmod_mpolyun_t poly, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
n_poly_struct * fq_nmod_mpolyun_leadcoeff_poly(const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
slong fq_nmod_mpolyun_lastdeg(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyun_is_nonzero_fq_nmod(const fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_scalar_mul_fq_nmod(fq_nmod_mpolyun_t A, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_shift_right(fq_nmod_mpolyun_t A, ulong s)
void fq_nmod_mpolyun_shift_left(fq_nmod_mpolyun_t A, ulong s)
void fq_nmod_mpolyun_mul_poly(fq_nmod_mpolyun_t A, const fq_nmod_mpolyun_t B, const fq_nmod_poly_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_divexact_poly(fq_nmod_mpolyun_t A, const fq_nmod_mpolyun_t B, const fq_nmod_poly_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_content_poly(fq_nmod_poly_t g, fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ctx)

The same operations for `fq_nmod_mpolyun_t`.
:::

### GCD algorithms

::: function
int _fq_nmod_mpoly_gcd_algo(fq_nmod_mpoly_t G, fq_nmod_mpoly_t Abar, fq_nmod_mpoly_t Bbar, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx, unsigned int algo)

Computes the GCD `G` and, if `Abar` and `Bbar` are not `NULL`, the
cofactors of the nonzero `A` and `B` with the algorithms in the bit mask
`algo` (the `MPOLY_GCD_USE_` constants). Returns 1 for success.
:::

::: function
int fq_nmod_mpolyn_gcd_brown_smprime(fq_nmod_mpolyn_t G, fq_nmod_mpolyn_t Abar, fq_nmod_mpolyn_t Bbar, fq_nmod_mpolyn_t A, fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyn_gcd_brown_lgprime(fq_nmod_mpolyn_t G, fq_nmod_mpolyn_t Abar, fq_nmod_mpolyn_t Bbar, fq_nmod_mpolyn_t A, fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

Brown's algorithm with evaluation points in $\mathbb{F}_q$, or with
reduction of the variable `var` modulo irreducible polynomials of
increasing degree over $\mathbb{F}_q$ (images over $\mathbb{F}_{q^n}$).
Return 1 for success, also setting the cofactors.
:::

::: function
int fq_nmod_mpolyu_gcdm_zippel(fq_nmod_mpolyu_t G, fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar, fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, fq_nmod_mpoly_ctx_t ctx, flint_rand_t randstate)
int fq_nmod_mpolyu_gcdp_zippel(fq_nmod_mpolyu_t G, fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar, fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, slong var, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t randstate)
int fq_nmod_mpolyu_gcdp_zippel_univar(fq_nmod_mpolyu_t G, fq_nmod_mpolyu_t Abar, fq_nmod_mpolyu_t Bbar, fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyu_gcdp_zippel_univar_no_cofactors(fq_nmod_mpolyu_t G, fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx)
nmod_gcds_ret_t fq_nmod_mpolyu_gcds_zippel(fq_nmod_mpolyu_t G, fq_nmod_mpolyu_t A, fq_nmod_mpolyu_t B, fq_nmod_mpolyu_t f, slong var, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t randstate, slong * degbound)

Zippel's sparse interpolation on `fq_nmod_mpolyu_t`: `gcdp` with points
in $\mathbb{F}_q$ (the `univar` versions for the base case of one main
variable, with or without cofactors), `gcds` solving for the coefficients
of the assumed form `f`, and `gcdm` trying `gcdp` first and then images
over extensions $\mathbb{F}_{q^n}$ chosen by an embedding chooser.
:::

::: function
int fq_nmod_mpolyl_gcd_zippel_smprime(fq_nmod_mpoly_t rG, const slong * rGdegs, fq_nmod_mpoly_t rAbar, fq_nmod_mpoly_t rBbar, const fq_nmod_mpoly_t A, const slong * Adegs, const fq_nmod_mpoly_t B, const slong * Bdegs, const fq_nmod_mpoly_t gamma, const slong * gammadegs, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyl_gcd_zippel_lgprime(fq_nmod_mpoly_t rG, const slong * rGdegs, fq_nmod_mpoly_t rAbar, fq_nmod_mpoly_t rBbar, const fq_nmod_mpoly_t A, const slong * Adegs, const fq_nmod_mpoly_t B, const slong * Bdegs, const fq_nmod_mpoly_t gamma, const slong * gammadegs, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyl_gcd_hensel_smprime(fq_nmod_mpoly_t G, slong Gdeg, fq_nmod_mpoly_t Abar, fq_nmod_mpoly_t Bbar, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)

Sparse interpolation with two main variables given degree bounds and the
leading coefficient of the GCD, over $\mathbb{F}_q$ or over extensions;
and the GCD by Hensel lifting over $\mathbb{F}_q$.
:::

### Interpolation

::: function
void _fq_nmod_mpoly_monomial_evals_cache(n_fq_poly_t E, const ulong * Aexps, flint_bitcnt_t Abits, slong Alen, const fq_nmod_struct * betas, slong start, slong stop, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_monomial_evals2_cache(n_fq_polyun_t E, const ulong * Aexps, flint_bitcnt_t Abits, slong Alen, const fq_nmod_struct * betas, slong m, const fq_nmod_mpoly_ctx_t ctx)
void n_fq_bpoly_eval_step_sep(n_fq_bpoly_t E, n_fq_polyun_t cur, const n_fq_polyun_t inc, const fq_nmod_mpoly_t A, const fq_nmod_ctx_t ctx)
void n_fq_polyun_zip_start(n_fq_polyun_t Z, n_fq_polyun_t H, slong req_images, const fq_nmod_ctx_t ctx)
int n_fq_polyu2n_add_zip_must_match(n_fq_polyun_t Z, const n_fq_bpoly_t A, slong cur_length, const fq_nmod_ctx_t ctx)
int n_fq_polyun_zip_solve(fq_nmod_mpoly_t A, n_fq_polyun_t Z, n_fq_polyun_t H, n_fq_polyun_t M, const fq_nmod_mpoly_ctx_t ctx)

Sparse interpolation helpers, as for `nmod_mpoly`: monomial evaluations
at powers of the evaluation points, the next bivariate image from the
current and incremental monomial values, and collection and solution of
the Vandermonde systems.
:::

::: function
void fq_nmod_mpolyn_interp_reduce_sm_poly(fq_nmod_poly_t E, const fq_nmod_mpolyn_t A, const fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_interp_lift_sm_poly(fq_nmod_mpolyn_t A, const fq_nmod_poly_t B, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyn_interp_crt_sm_poly(slong * lastdeg_, fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, const fq_nmod_poly_t A, const fq_nmod_poly_t modulus, const fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_interp_lift_sm_bpoly(fq_nmod_mpolyn_t F, n_fq_bpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyn_interp_crt_sm_bpoly(slong * lastdeg, fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, const n_fq_bpoly_t A, const n_fq_poly_t modulus, n_fq_poly_t alphapow, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_interp_reduce_sm_mpolyn(fq_nmod_mpolyn_t E, fq_nmod_mpolyn_t A, slong var, fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_interp_lift_sm_mpolyn(fq_nmod_mpolyn_t A, fq_nmod_mpolyn_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyn_interp_crt_sm_mpolyn(slong * lastdeg_, fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_mpolyn_t A, slong var, fq_nmod_poly_t modulus, const fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyn_interp_mcrt_sm_mpoly(slong * lastdeg_, fq_nmod_mpolyn_t F, fq_nmod_mpoly_t A, const n_fq_poly_t modulus, n_fq_poly_t alphapow, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_interp_reduce_sm_mpolyu(fq_nmod_mpolyu_t B, fq_nmod_mpolyun_t A, fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyn_interp_lift_sm_mpoly(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyun_interp_lift_sm_mpolyu(fq_nmod_mpolyun_t A, const fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpolyun_interp_crt_sm_mpolyu(slong * lastdeg, fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_mpolyu_t A, fq_nmod_poly_t modulus, fq_nmod_t alpha, const fq_nmod_mpoly_ctx_t ctx)

Interpolation steps of Brown's and Zippel's algorithms with points
$\alpha \in \mathbb{F}_q$: evaluation of the last variable (`reduce`),
first image (`lift`) and Chinese remaindering with a new image (`crt`,
returning whether the interpolant changed and updating `lastdeg`).
:::

::: function
void fq_nmod_mpolyn_interp_reduce_lg_poly(fq_nmod_poly_t E, const fq_nmod_mpoly_ctx_t ectx, fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb)
void fq_nmod_mpolyn_interp_lift_lg_poly(slong * lastdeg_, fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t B, const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb)
int fq_nmod_mpolyn_interp_crt_lg_poly(slong * lastdeg_, fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_poly_t modulus, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t A, const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb)
void fq_nmod_mpolyn_interp_lift_lg_bpoly(slong * lastdeg_, fq_nmod_mpolyn_t F, const fq_nmod_mpoly_ctx_t smctx, n_fq_bpoly_t A, const fq_nmod_mpoly_ctx_t lgctx, const bad_fq_nmod_embed_t emb)
int fq_nmod_mpolyn_interp_crt_lg_bpoly(slong * lastdeg, fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, n_fq_poly_t modulus, const fq_nmod_mpoly_ctx_t smctx, n_fq_bpoly_t A, const fq_nmod_mpoly_ctx_t lgctx, const bad_fq_nmod_embed_t emb)
void fq_nmod_mpolyn_interp_reduce_lg_mpolyn(fq_nmod_mpolyn_t E, const fq_nmod_mpoly_ctx_t ectx, fq_nmod_mpolyn_t A, slong var, const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb)
void fq_nmod_mpolyn_interp_lift_lg_mpolyn(slong * lastdeg_, fq_nmod_mpolyn_t A, slong var, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb)
int fq_nmod_mpolyn_interp_crt_lg_mpolyn(slong * lastdeg_, fq_nmod_mpolyn_t F, fq_nmod_mpolyn_t T, fq_nmod_poly_t modulus, slong var, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb)
void fq_nmod_mpolyn_interp_reduce_lg_mpoly(fq_nmod_mpoly_t A, fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb)
void fq_nmod_mpolyun_interp_reduce_lg_mpolyu(fq_nmod_mpolyu_t A, fq_nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ectx, const fq_nmod_mpoly_ctx_t ctx, const bad_fq_nmod_embed_t emb)
void fq_nmod_mpolyn_interp_lift_lg_mpoly(fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ctx,fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb)
void fq_nmod_mpolyun_interp_lift_lg_mpolyu(fq_nmod_mpolyun_t A, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t B, const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb)
int fq_nmod_mpolyun_interp_crt_lg_mpolyu(slong * lastdeg, fq_nmod_mpolyun_t F, fq_nmod_mpolyun_t T, fq_nmod_poly_t m, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ectx, const bad_fq_nmod_embed_t emb)
int fq_nmod_mpolyun_interp_mcrt_lg_mpolyu(slong * lastdeg, fq_nmod_mpolyun_t H, const fq_nmod_mpoly_ctx_t ctx, fq_nmod_poly_t m, fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ectx, bad_fq_nmod_embed_t emb)

The same steps for images over an extension $\mathbb{F}_{q^n}$ (context
`ectx`) obtained with the embedding `emb`: reduction modulo the
irreducible polynomial $h$ in the last variable, lifting of an image, and
Chinese remaindering in $\mathbb{F}_q[x]$ with the product `modulus` of
the polynomials used so far.
:::

::: function
void nmod_mpolyn_interp_reduce_lg_poly(fq_nmod_poly_t E, const fq_nmod_ctx_t fqctx, const nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_interp_lift_lg_poly(slong * lastdeg_, nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx, const fq_nmod_poly_t B, const fq_nmod_ctx_t fqctx)
int nmod_mpolyn_interp_crt_lg_poly(slong * lastdeg_, nmod_mpolyn_t F, nmod_mpolyn_t T, n_poly_t modulus, const nmod_mpoly_ctx_t ctx, fq_nmod_poly_t A, const fq_nmod_ctx_t fqctx)
void nmod_mpolyn_interp_lift_lg_bpoly(slong * lastdeg_, nmod_mpolyn_t F, const nmod_mpoly_ctx_t smctx, n_fq_bpoly_t A, const fq_nmod_mpoly_ctx_t lgctx)
int nmod_mpolyn_interp_crt_lg_bpoly(slong * lastdeg, nmod_mpolyn_t F, nmod_mpolyn_t T, n_fq_poly_t modulus, const nmod_mpoly_ctx_t smctx, n_fq_bpoly_t A, const fq_nmod_mpoly_ctx_t lgctx)
void nmod_mpolyn_interp_reduce_lg_mpolyn(fq_nmod_mpolyn_t E, fq_nmod_mpoly_ctx_t ectx, nmod_mpolyn_t A, slong var, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_interp_lift_lg_mpolyn(slong * lastdeg, nmod_mpolyn_t A, slong var, const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ectx)
int nmod_mpolyn_interp_crt_lg_mpolyn(slong * lastdeg_, nmod_mpolyn_t F, nmod_mpolyn_t T, n_poly_t modulus, slong var, const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyn_t A, const fq_nmod_mpoly_ctx_t ectx)
void nmod_mpolyn_interp_reduce_lg_mpoly(fq_nmod_mpoly_t A, nmod_mpolyn_t B, const fq_nmod_mpoly_ctx_t ffctx, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_interp_reduce_lg_mpolyu(fq_nmod_mpolyu_t A, nmod_mpolyun_t B, const fq_nmod_mpoly_ctx_t ffctx, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_interp_lift_lg_mpoly(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx, fq_nmod_mpoly_t Ap, const fq_nmod_mpoly_ctx_t ctxp)
void nmod_mpolyun_interp_lift_lg_mpolyu(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t Ap, const fq_nmod_mpoly_ctx_t ctxp)
int nmod_mpolyun_interp_crt_lg_mpolyu(slong * lastdeg, nmod_mpolyun_t F, nmod_mpolyun_t T, n_poly_t m, const nmod_mpoly_ctx_t ctx, fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ffctx)
int nmod_mpolyn_interp_mcrt_lg_mpoly(slong * lastdeg_, nmod_mpolyn_t H, const nmod_mpoly_ctx_t smctx, const n_poly_t m, const ulong * inv_m_eval,fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t lgctx)
int nmod_mpolyun_interp_mcrt_lg_mpolyu(slong * lastdeg, nmod_mpolyun_t H, const nmod_mpoly_ctx_t ctx, n_poly_t m, fq_nmod_mpolyu_t A, const fq_nmod_mpoly_ctx_t ctxp)

The `lgprime` interpolation steps of the `nmod_mpoly` GCD algorithms
(declared here because they use `fq_nmod_mpoly`): a polynomial over
$\mathbb{F}_p$ is reduced modulo an irreducible polynomial $m$ in its last
variable, giving a polynomial over $\mathbb{F}_p[x]/(m) \cong
\mathbb{F}_{p^{\deg m}}$ in the context `ectx` (or `ffctx`), and images
are combined by Chinese remaindering in $\mathbb{F}_p[x]$.
:::
