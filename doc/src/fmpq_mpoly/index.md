---
title: fmpq_mpoly.h
header: fmpq_mpoly.h
description: Sparse multivariate polynomials over the rational numbers, stored as a rational content times a primitive integer polynomial.
category: Multivariate polynomials
order: 5
---

This module implements sparse polynomials in $\Q[x_0, \ldots, x_{n-1}]$.
A polynomial is stored as a rational number, its *content*, times a
primitive integer polynomial with positive leading coefficient, an
`fmpz_mpoly_t` (see [](fmpz_mpoly)). Almost all operations reduce to the
corresponding operations on the integer part, followed by an adjustment of
the content, so the algorithms, the monomial orderings, the exponent
representation and the performance characteristics are those of
`fmpz_mpoly`; the guide [](fmpz_mpoly/getting-started) applies with
`fmpz` replaced by `fmpq`. Factorisation is in [](fmpq_mpoly_factor).

For polynomials over an arbitrary generic ring, including $\Q$, see
`gr_mpoly`; `fmpq_mpoly` is the specialised and faster implementation for
$\Q$.

## Types and conventions {#types}

::: type
fmpq_mpoly_struct
fmpq_mpoly_t

A polynomial $A = c \cdot Z$ with fields `content` (the rational $c$, an
`fmpq_t`) and `zpoly` (the integer polynomial $Z$, an `fmpz_mpoly_t`). The
representation is canonical when either both are zero, or $c \ne 0$ and
$Z$ is a canonical `fmpz_mpoly` whose coefficients have GCD 1 and whose
leading coefficient is positive. Then the sign of $A$ is carried by $c$.
The functions of this module return canonical polynomials when given
canonical input, except where noted.
:::

::: type
fmpq_mpoly_ctx_struct
fmpq_mpoly_ctx_t

The context (parent ring) of an `fmpq_mpoly`: the number of variables and
the monomial ordering. It contains a single field `zctx`, the
`fmpz_mpoly_ctx_t` that is used for the integer parts.
:::

Variables are numbered from 0; variable 0 is the most significant one in
the ordering. Exponent vectors passed as `ulong *`, `slong *` or arrays of
`fmpz` are indexed by variable number. Output polynomials may alias inputs
unless stated otherwise.

## Context object {#context}

::: function
void fmpq_mpoly_ctx_init(fmpq_mpoly_ctx_t ctx, slong nvars, const ordering_t ord)

Initialises a context for polynomials in `nvars` variables with the
monomial ordering `ord`, one of `ORD_LEX`, `ORD_DEGLEX` and
`ORD_DEGREVLEX`.
:::

::: function
void fmpq_mpoly_ctx_init_rand(fmpq_mpoly_ctx_t ctx, flint_rand_t state, slong max_nvars)

Initialises a context with a random number of variables between 0 and
`max_nvars` and a random ordering, for testing.
:::

::: function
slong fmpq_mpoly_ctx_nvars(const fmpq_mpoly_ctx_t ctx)
ordering_t fmpq_mpoly_ctx_ord(const fmpq_mpoly_ctx_t ctx)

Return the number of variables and the ordering of the context.
:::

::: function
void fmpq_mpoly_ctx_clear(fmpq_mpoly_ctx_t ctx)

Frees the memory used by the context.
:::

## Memory management {#memory}

::: function
void fmpq_mpoly_init(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_init2(fmpq_mpoly_t A, slong alloc, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_init3(fmpq_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fmpq_mpoly_ctx_t ctx)

Initialise $A$ to zero, in the last two versions with room for `alloc`
terms, and in the last with exponent fields of `bits` bits (otherwise at
least `MPOLY_MIN_BITS`).
:::

::: function
void fmpq_mpoly_fit_length(fmpq_mpoly_t A, slong len, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_fit_bits(fmpq_mpoly_t A, flint_bitcnt_t bits, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_realloc(fmpq_mpoly_t A, slong alloc, const fmpq_mpoly_ctx_t ctx)

Ensure room for `len` terms, ensure exponent fields of at least `bits`
bits (repacking the terms), and reallocate the integer part to exactly
`alloc` terms (its length must not exceed `alloc`).
:::

::: function
void fmpq_mpoly_clear(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Frees the memory used by $A$.
:::

## Input and output {#io}

The variable names `x` are indexed by variable number. If `x` is `NULL`,
the names `x1`, `x2`, ... are used.

::: function
char * fmpq_mpoly_get_str_pretty(const fmpq_mpoly_t A, const char ** x, const fmpq_mpoly_ctx_t ctx)

Returns a string representation of $A$, to be freed with `flint_free`.
:::

::: function
int fmpq_mpoly_fprint_pretty(FILE * file, const fmpq_mpoly_t A, const char ** x, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_print_pretty(const fmpq_mpoly_t A, const char ** x, const fmpq_mpoly_ctx_t ctx)

Print $A$ to `file` or to `stdout`, with coefficients such as `1/2*x^2`.
:::

::: function
int fmpq_mpoly_set_str_pretty(fmpq_mpoly_t A, const char * str, const char ** x, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to the polynomial given by the null-terminated string `str` and
returns 0, or returns $-1$ if parsing fails, in which case the value of
$A$ is unspecified. The syntax is that of `fmpz_mpoly_set_str_pretty`:
integers, the names in `x`, parentheses, `+`, `-`, `*`, `/` and `^`
followed by a nonnegative integer; division by a polynomial must be exact,
and division by integers gives rational coefficients.
:::

## Basic manipulation {#basic}

::: function
void fmpq_mpoly_gen(fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to the variable of index `var`.
:::

::: function
int fmpq_mpoly_is_gen(const fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)

If $\mathit{var} \ge 0$, returns whether $A$ is the variable of index
`var`; if $\mathit{var} < 0$, returns whether $A$ is any variable.
:::

::: function
void fmpq_mpoly_set(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_swap(fmpq_mpoly_t A, fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

Set $A$ to $B$, and swap $A$ and $B$ in constant time.
:::

::: function
int fmpq_mpoly_equal(const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

Returns whether $A = B$.
:::

## Constants {#constants}

::: function
int fmpq_mpoly_is_fmpq(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Returns whether $A$ is a constant.
:::

::: function
void fmpq_mpoly_get_fmpq(fmpq_t c, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Sets $c$ to the constant $A$; throws if $A$ is not constant.
:::

::: function
void fmpq_mpoly_set_fmpq(fmpq_mpoly_t A, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_set_fmpz(fmpq_mpoly_t A, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_set_ui(fmpq_mpoly_t A, ulong c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_set_si(fmpq_mpoly_t A, slong c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_zero(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_one(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Set $A$ to the constant $c$, to 0 or to 1.
:::

::: function
int fmpq_mpoly_equal_fmpq(const fmpq_mpoly_t A, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_equal_fmpz(const fmpq_mpoly_t A, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_equal_ui(const fmpq_mpoly_t A, ulong c, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_equal_si(const fmpq_mpoly_t A, slong c, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_is_zero(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_is_one(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Return whether $A$ is the constant $c$, 0 or 1.
:::

## Degrees {#degrees}

The degree of the zero polynomial is $-1$.

::: function
int fmpq_mpoly_degrees_fit_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Returns whether the degree of $A$ in each variable fits an `slong`.
:::

::: function
void fmpq_mpoly_degrees_fmpz(fmpz ** degs, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_degrees_si(slong * degs, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Set `degs[i]` to the degree of $A$ in the variable of index $i$.
:::

::: function
void fmpq_mpoly_degree_fmpz(fmpz_t deg, const fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)
slong fmpq_mpoly_degree_si(const fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)

Set or return the degree of $A$ in the variable of index `var`.
:::

::: function
int fmpq_mpoly_total_degree_fits_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_total_degree_fmpz(fmpz_t tdeg, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
slong fmpq_mpoly_total_degree_si(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Return whether the total degree of $A$ fits an `slong`, and set or return
the total degree.
:::

::: function
void fmpq_mpoly_used_vars(int * used, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Sets `used[i]` to a nonzero value if the variable of index $i$ occurs in
$A$, and to zero otherwise.
:::

## Coefficients {#coefficients}

::: function
void fmpq_mpoly_get_denominator(fmpz_t d, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Sets $d$ to the smallest positive integer such that $dA$ has integer
coefficients (the denominator of the content).
:::

::: function
int fmpq_mpoly_is_monic(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Returns whether the leading coefficient of $A$ is 1 (false for zero).
:::

::: function
void fmpq_mpoly_get_coeff_fmpq_monomial(fmpq_t c, const fmpq_mpoly_t A, const fmpq_mpoly_t M, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_set_coeff_fmpq_monomial(fmpq_mpoly_t A, const fmpq_t c, const fmpq_mpoly_t M, const fmpq_mpoly_ctx_t ctx)

Get or set the coefficient in $A$ of the monomial of the single-term
polynomial $M$ (whose coefficient is ignored). Throw if $M$ does not have
exactly one term.
:::

::: function
void fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t A, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_get_coeff_fmpq_ui(fmpq_t c, const fmpq_mpoly_t A, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
void _fmpq_mpoly_get_coeff_fmpq_fmpz(fmpq_t c, const fmpq_mpoly_t A, const fmpz * exp, const fmpq_mpoly_ctx_t ctx)

Set $c$ to the coefficient of the monomial with exponent vector `exp`
(zero if it does not occur). The underscore version takes the exponents
as a plain array of `fmpz`.
:::

::: function
void fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t A, const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_set_coeff_fmpq_ui(fmpq_mpoly_t A, const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
void _fmpq_mpoly_set_coeff_fmpq_fmpz(fmpq_mpoly_t A, const fmpq_t c, const fmpz * exp, const fmpq_mpoly_ctx_t ctx)

Set the coefficient of the monomial with exponent vector `exp` to $c$,
inserting or removing a term as needed.
:::

::: function
void fmpq_mpoly_get_coeff_vars_ui(fmpq_mpoly_t C, const fmpq_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fmpq_mpoly_ctx_t ctx)

Sets $C$ to the coefficient of $\prod_k x_{\mathit{vars}[k]}^{\mathit{exps}[k]}$
in $A$ viewed as a polynomial in the variables `vars` (a polynomial in the
other variables). The `length` variables in `vars` must be distinct, with
$0 < \mathit{length} \le n$.
:::

## Conversion {#conversion}

::: function
int fmpq_mpoly_is_fmpq_poly(const fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)

Returns whether $A$ involves no variable other than the one of index
`var`.
:::

::: function
int fmpq_mpoly_get_fmpq_poly(fmpq_poly_t A, const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to $B$ as a univariate polynomial in the variable of index
`var` and returns 1. Returns 0 if $B$ involves another variable or the
degree does not fit an `slong`.
:::

::: function
void fmpq_mpoly_set_fmpq_poly(fmpq_mpoly_t A, const fmpq_poly_t B, slong var, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to the univariate polynomial $B$ in the variable of index `var`.
:::

## Comparison {#comparison}

::: function
int fmpq_mpoly_cmp(const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

Returns $1$, $-1$ or $0$ according to a fixed total order on polynomials
that agrees with the monomial ordering for monomials with equal
coefficients.
:::

## Container operations {#container}

These functions give access to the terms of $A$ by index (from 0, in
decreasing monomial order) and to the internal representation. An index
outside $[0, \mathit{length})$ makes them throw. The mutating functions
may leave $A$ non-canonical; after setting up terms with distinct
exponents in any order, call `fmpq_mpoly_sort_terms` and
`fmpq_mpoly_combine_like_terms`, or `fmpq_mpoly_reduce` if the terms are
already sorted, nonzero and distinct.

::: function
fmpq * fmpq_mpoly_content_ref(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
fmpz_mpoly_struct * fmpq_mpoly_zpoly_ref(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
fmpz * fmpq_mpoly_zpoly_term_coeff_ref(fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

Return pointers to the content, to the integer part, and to the integer
coefficient of the term of index $i$ of the integer part.
:::

::: function
int fmpq_mpoly_is_canonical(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_assert_canonical(const fmpq_mpoly_t poly, const fmpq_mpoly_ctx_t ctx)

Return whether $A$ is in the canonical form described in the section on
types, respectively throw an exception describing the defect if `poly`
is not.
:::

::: function
slong fmpq_mpoly_length(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Returns the number of terms of $A$.
:::

::: function
void fmpq_mpoly_resize(fmpq_mpoly_t A, slong new_length, const fmpq_mpoly_ctx_t ctx)

Sets the length of $A$ to `new_length`, deleting terms at the end or
appending terms with zero coefficient and zero exponents.
:::

::: function
void fmpq_mpoly_get_term_coeff_fmpq(fmpq_t c, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_set_term_coeff_fmpq(fmpq_mpoly_t A, slong i, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)

Get or set the coefficient of the term of index $i$.
:::

::: function
int fmpq_mpoly_term_exp_fits_si(const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_term_exp_fits_ui(const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

Return whether every exponent of the term of index $i$ fits an `slong`,
respectively a `ulong`.
:::

::: function
void fmpq_mpoly_get_term_exp_fmpz(fmpz ** exps, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_get_term_exp_ui(ulong * exps, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_get_term_exp_si(slong * exps, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

Set `exps` to the exponent vector of the term of index $i$. The `ui` and
`si` versions throw if an exponent does not fit.
:::

::: function
ulong fmpq_mpoly_get_term_var_exp_ui(const fmpq_mpoly_t A, slong i, slong var, const fmpq_mpoly_ctx_t ctx)
slong fmpq_mpoly_get_term_var_exp_si(const fmpq_mpoly_t A, slong i, slong var, const fmpq_mpoly_ctx_t ctx)

Return the exponent of the variable of index `var` in the term of index
$i$; throw if it does not fit.
:::

::: function
void fmpq_mpoly_set_term_exp_fmpz(fmpq_mpoly_t A, slong i, fmpz * const * exps, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_set_term_exp_ui(fmpq_mpoly_t A, slong i, const ulong * exps, const fmpq_mpoly_ctx_t ctx)

Set the exponent vector of the term of index $i$, which may break the
ordering of the terms.
:::

::: function
void fmpq_mpoly_get_term(fmpq_mpoly_t M, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_get_term_monomial(fmpq_mpoly_t M, const fmpq_mpoly_t A, slong i, const fmpq_mpoly_ctx_t ctx)

Set $M$ to the term of index $i$, respectively its monomial with
coefficient 1.
:::

::: function
void fmpq_mpoly_push_term_fmpq_fmpz(fmpq_mpoly_t A, const fmpq_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_fmpq_ffmpz(fmpq_mpoly_t A, const fmpq_t c, const fmpz * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_fmpz_fmpz(fmpq_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_fmpz_ffmpz(fmpq_mpoly_t A, const fmpz_t c, const fmpz * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_ui_fmpz(fmpq_mpoly_t A, ulong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_ui_ffmpz(fmpq_mpoly_t A, ulong c, const fmpz * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_si_fmpz(fmpq_mpoly_t A, slong c, fmpz * const * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_si_ffmpz(fmpq_mpoly_t A, slong c, const fmpz * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_fmpq_ui(fmpq_mpoly_t A, const fmpq_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_fmpz_ui(fmpq_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_ui_ui(fmpq_mpoly_t A, ulong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_push_term_si_ui(fmpq_mpoly_t A, slong c, const ulong * exp, const fmpq_mpoly_ctx_t ctx)

Append a term with coefficient $c$ and exponent vector `exp` (an array of
pointers to `fmpz`, an array of `fmpz`, or an array of `ulong`). The
integer part and the denominator of the content are rescaled as needed, so
the cost per term is constant on average when the denominators of the
pushed coefficients are bounded. The result is not reduced; finish with
`fmpq_mpoly_reduce` (or sort and combine first).
:::

::: function
void fmpq_mpoly_reduce(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Moves the content and sign of the integer part into the content, so that
the integer part is primitive with positive leading coefficient. If the
terms of $A$ are nonzero, sorted and have distinct exponents, the result
is canonical.
:::

::: function
void fmpq_mpoly_reduce_easy(fmpq_mpoly_t A, slong easy_length, const fmpq_mpoly_ctx_t ctx)

A shortcut for `fmpq_mpoly_reduce` after an addition: if the integer part
has `easy_length` terms (no cancellation happened), it is assumed to be
primitive and only its sign is normalised; otherwise `fmpq_mpoly_reduce`
is called. The integer part must have at most `easy_length` terms.
:::

::: function
void fmpq_mpoly_sort_terms(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Sorts the terms of the integer part into decreasing monomial order,
without combining like terms, removing zero terms or reducing.
:::

::: function
void fmpq_mpoly_combine_like_terms(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Combines adjacent terms with equal exponents, removes zero terms and
calls `fmpq_mpoly_reduce`. If the terms were sorted, the result is
canonical.
:::

## Random generation {#random}

The parameter `coeff_bits` is only a suggestion for the size of the
coefficients, and the content is a random nonzero rational.

::: function
void fmpq_mpoly_randtest_bound(fmpq_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t coeff_bits, ulong exp_bound, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_randtest_bounds(fmpq_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t coeff_bits, ulong * exp_bounds, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_randtest_bits(fmpq_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t coeff_bits, flint_bitcnt_t exp_bits, const fmpq_mpoly_ctx_t ctx)

Set $A$ to a random polynomial with at most `length` terms and exponents
below `exp_bound`, below `exp_bounds[i]` for variable $i$, or of at most
`exp_bits` bits.
:::

## Addition and subtraction {#addition}

::: function
void fmpq_mpoly_add_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_add_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_add_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_add_si(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_sub_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_sub_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_sub_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_sub_si(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong c, const fmpq_mpoly_ctx_t ctx)

Set $A$ to $B + c$, respectively $B - c$.
:::

::: function
void fmpq_mpoly_add(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_sub(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx)

Set $A$ to $B + C$, respectively $B - C$. The integer parts are combined
with the cofactors of the GCD of the contents.
:::

## Scalar operations {#scalar}

Unless $B$ or the scalar is zero, these functions only change the
content, so they take constant time when $A$ and $B$ are aliased.

::: function
void fmpq_mpoly_neg(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to $-B$.
:::

::: function
void fmpq_mpoly_scalar_mul_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_scalar_mul_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_scalar_mul_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_scalar_mul_si(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong c, const fmpq_mpoly_ctx_t ctx)

Set $A$ to $cB$.
:::

::: function
void fmpq_mpoly_scalar_div_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_scalar_div_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_scalar_div_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong c, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_scalar_div_si(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong c, const fmpq_mpoly_ctx_t ctx)

Set $A$ to $B/c$ for nonzero $c$.
:::

::: function
void fmpq_mpoly_make_monic(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to $B$ divided by its leading coefficient; throws if $B$ is zero.
:::

## Differentiation and integration {#calculus}

::: function
void fmpq_mpoly_derivative(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to the derivative of $B$ with respect to the variable of index
`var`.
:::

::: function
void fmpq_mpoly_integral(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to the antiderivative of $B$ with respect to the variable of
index `var` that has no term free of that variable.
:::

## Evaluation and composition {#evaluation}

These functions return 1 for success and 0 if the result would be
unreasonably large (for example a huge power of a rational with a large
height).

::: function
int fmpq_mpoly_evaluate_all_fmpq(fmpq_t ev, const fmpq_mpoly_t A, fmpq * const * vals, const fmpq_mpoly_ctx_t ctx)

Sets `ev` to $A$ evaluated at $x_i = \mathit{vals}[i]$ for all $i$.
:::

::: function
int fmpq_mpoly_evaluate_one_fmpq(fmpq_mpoly_t A, const fmpq_mpoly_t B, slong var, const fmpq_t val, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to $B$ with the variable of index `var` replaced by `val`.
:::

::: function
int fmpq_mpoly_compose_fmpq_poly(fmpq_poly_t A, const fmpq_mpoly_t B, fmpq_poly_struct * const * C, const fmpq_mpoly_ctx_t ctxB)

Sets $A$ to $B(C_0, \ldots, C_{n-1})$ for univariate polynomials $C_i$.
:::

::: function
int fmpq_mpoly_compose_fmpq_mpoly(fmpq_mpoly_t A, const fmpq_mpoly_t B, fmpq_mpoly_struct * const * C, const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC)

Sets $A$ to $B(C_0, \ldots, C_{n-1})$, where $B$ is in `ctxB` and $A$
and the $C_i$ are in `ctxAC`. $A$ must not alias $B$ or any $C_i$.
:::

::: function
void fmpq_mpoly_compose_fmpq_mpoly_gen(fmpq_mpoly_t A, const fmpq_mpoly_t B, const slong * c, const fmpq_mpoly_ctx_t ctxB, const fmpq_mpoly_ctx_t ctxAC)

Sets $A$ to $B$ with the variable of index $i$ of `ctxB` replaced by the
variable of index `c[i]` of `ctxAC`, or by zero if `c[i]` is negative.
Each nonnegative `c[i]` must be less than the number of variables of
`ctxAC`.
:::

## Multiplication {#multiplication}

::: function
void fmpq_mpoly_mul(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_t C, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to $BC$, multiplying the integer parts with `fmpz_mpoly_mul`
(see [](fmpz_mpoly/mul)).
:::

## Powering {#powering}

::: function
int fmpq_mpoly_pow_fmpz(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz_t k, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_pow_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B, ulong k, const fmpq_mpoly_ctx_t ctx)

Set $A$ to $B^k$. Return 0 if $k$ is negative (in the `fmpz` version) or
the result would be unreasonably large, and 1 otherwise.
:::

## Division {#division}

Division by zero throws. Since $\Q$ is a field, division with remainder
is by the leading term of the divisor without any scaling: the
remainder has no term divisible by the leading monomial of the divisor.

::: function
int fmpq_mpoly_divides(fmpq_mpoly_t Q, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

If $B$ divides $A$, sets $Q = A/B$ and returns 1; otherwise sets $Q$ to
zero and returns 0. See [](fmpz_mpoly/division).
:::

::: function
void fmpq_mpoly_div(fmpq_mpoly_t Q, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_divrem(fmpq_mpoly_t Q, fmpq_mpoly_t R, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

Set $Q$ (and $R$) to the quotient (and remainder) of $A$ by $B$, with
$A = QB + R$ and no term of $R$ divisible by the leading monomial of $B$.
They are computed with `fmpz_mpoly_quasidiv` and `fmpz_mpoly_quasidivrem`
on the integer parts.
:::

::: function
void fmpq_mpoly_divrem_ideal(fmpq_mpoly_struct ** Q, fmpq_mpoly_t R, const fmpq_mpoly_t A, fmpq_mpoly_struct * const * B, slong len, const fmpq_mpoly_ctx_t ctx)

Generalised division by the `len` divisors $B_i$: sets the $Q_i$ and $R$
with $A = \sum_i Q_i B_i + R$ and no term of $R$ divisible by the leading
monomial of any $B_i$.
:::

## Square roots {#sqrt}

::: function
int fmpq_mpoly_sqrt(fmpq_mpoly_t Q, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

If $A$ is the square of a polynomial over $\Q$, sets $Q$ to the square
root with positive leading coefficient and returns 1; otherwise sets $Q$
to zero and returns 0.
:::

::: function
int fmpq_mpoly_is_square(const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Returns whether $A$ is the square of a polynomial over $\Q$.
:::

## Content and GCD {#gcd}

The GCD functions return the monic GCD (zero for two zero inputs),
computed with the function of the same name on the integer parts (see
[](fmpz_mpoly/gcd)). They return 1 for success and 0 for failure, which
happens only for unreasonably large inputs or when the requested algorithm
does not apply; $G$ is then unspecified.

::: function
void fmpq_mpoly_content(fmpq_t g, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Sets $g$ to the nonnegative content of $A$ (the absolute value of the
stored content).
:::

::: function
void fmpq_mpoly_term_content(fmpq_mpoly_t M, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Sets $M$ to the GCD of the terms of $A$: zero if $A$ is zero, and
otherwise the largest monomial dividing every term, with coefficient 1.
:::

::: function
int fmpq_mpoly_content_vars(fmpq_mpoly_t g, const fmpq_mpoly_t A, slong * vars, slong vars_length, const fmpq_mpoly_ctx_t ctx)

Sets $g$ to the GCD of the coefficients of $A$ viewed as a polynomial in
the variables `vars`, so that $g$ does not involve these variables.
Returns 1 for success and 0 for failure.
:::

::: function
int fmpq_mpoly_gcd(fmpq_mpoly_t G, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

Sets $G$ to the monic GCD of $A$ and $B$.
:::

::: function
int fmpq_mpoly_gcd_cofactors(fmpq_mpoly_t G, fmpq_mpoly_t Abar, fmpq_mpoly_t Bbar, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

Sets $G$ to the monic GCD of $A$ and $B$ and also $\bar A = A/G$ and
$\bar B = B/G$.
:::

::: function
int fmpq_mpoly_gcd_brown(fmpq_mpoly_t G, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_gcd_hensel(fmpq_mpoly_t G, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_gcd_subresultant(fmpq_mpoly_t G, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_gcd_zippel(fmpq_mpoly_t G, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
int fmpq_mpoly_gcd_zippel2(fmpq_mpoly_t G, const fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)

Set $G$ to the monic GCD of $A$ and $B$ using one particular algorithm.
The warning on `fmpz_mpoly_gcd_zippel` applies to `fmpq_mpoly_gcd_zippel`.
:::

::: function
int fmpq_mpoly_resultant(fmpq_mpoly_t R, const fmpq_mpoly_t A, const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)

Sets $R$ to the resultant of $A$ and $B$ with respect to the variable of
index `var`, computed with `fmpz_mpoly_resultant` on the integer parts and
multiplied by the appropriate powers of the contents. Returns 1 for
success and 0 for failure.
:::

::: function
int fmpq_mpoly_discriminant(fmpq_mpoly_t R, const fmpq_mpoly_t A, slong var, const fmpq_mpoly_ctx_t ctx)

Sets $R$ to the discriminant of $A$ with respect to the variable of index
`var`. Returns 1 for success and 0 for failure.
:::

::: function
void fmpq_mpoly_inflate(fmpq_mpoly_t A, const fmpq_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to $B$ with each exponent $e$ of variable $i$ replaced by
$\mathit{shift}[i] + \mathit{stride}[i] \cdot e$.
:::

## Univariate polynomials with multivariate coefficients {#univar}

An `fmpq_mpoly_univar_t` holds a polynomial in one main variable with
`fmpq_mpoly_t` coefficients in the other variables, which is useful for
rewriting an element of $\Q[x_0, \ldots, x_{n-1}]$ as an element of
$\Q[x_0, \ldots, \hat{x}_v, \ldots, x_{n-1}][x_v]$ and back.

::: type
fmpq_mpoly_univar_struct
fmpq_mpoly_univar_t

A sparse univariate polynomial with fields `coeffs` (nonzero polynomials
in the context of the multivariate polynomial, not involving the main
variable), `exps` (strictly decreasing nonnegative `fmpz` exponents),
`alloc` and `length`.
:::

::: function
void fmpq_mpoly_univar_init(fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_univar_clear(fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_univar_swap(fmpq_mpoly_univar_t A, fmpq_mpoly_univar_t B, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_univar_fit_length(fmpq_mpoly_univar_t A, slong length, const fmpq_mpoly_ctx_t ctx)

Initialise to zero, free, swap in constant time, and ensure room for
`length` terms.
:::

::: function
void fmpq_mpoly_univar_print_pretty(const fmpq_mpoly_univar_t A, const char ** x, const fmpq_mpoly_ctx_t ctx)

Prints $A$ as a sum of terms `(c)*X^e`.
:::

::: function
void fmpq_mpoly_univar_assert_canonical(fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)

Throws if the exponents of $A$ are not nonnegative and strictly
decreasing or a coefficient is not canonical.
:::

::: function
void fmpq_mpoly_to_univar(fmpq_mpoly_univar_t A, const fmpq_mpoly_t B, slong var, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to $B$ as a polynomial in the variable of index `var`.
:::

::: function
void fmpq_mpoly_from_univar(fmpq_mpoly_t A, const fmpq_mpoly_univar_t B, slong var, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_from_univar_bits(fmpq_mpoly_t A, flint_bitcnt_t Abits, const fmpq_mpoly_univar_t B, slong var, const fmpq_mpoly_ctx_t ctx)

Set $A$ to the polynomial represented by $B$, with the variable of index
`var` as main variable. The coefficients of $B$ must not involve that
variable. The second version packs the result into `Abits` bits, which
must suffice.
:::

::: function
int fmpq_mpoly_univar_degree_fits_si(const fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)
slong fmpq_mpoly_univar_length(const fmpq_mpoly_univar_t A, const fmpq_mpoly_ctx_t ctx)
slong fmpq_mpoly_univar_get_term_exp_si(fmpq_mpoly_univar_t A, slong i, const fmpq_mpoly_ctx_t ctx)

Return whether the degree of $A$ fits an `slong`, the number of terms,
and the exponent of the term of index $i$.
:::

::: function
void fmpq_mpoly_univar_get_term_coeff(fmpq_mpoly_t c, const fmpq_mpoly_univar_t A, slong i, const fmpq_mpoly_ctx_t ctx)
void fmpq_mpoly_univar_swap_term_coeff(fmpq_mpoly_t c, fmpq_mpoly_univar_t A, slong i, const fmpq_mpoly_ctx_t ctx)

Set $c$ to, respectively swap $c$ with, the coefficient of the term of
index $i$.
:::

## Internal functions {#internal}

::: function
void _fmpq_mpoly_push_rescale(fmpq_mpoly_t A, fmpq_t C, const fmpq_mpoly_ctx_t ctx)

Prepares appending a term with coefficient $C$: multiplies the integer
part by the numerator of the content (which becomes 1) and by the
denominator of $C \cdot \mathrm{den}(c)$, so that the new term can be
appended with the integer coefficient left in the numerator of $C$.
:::

::: function
void _fmpq_mpoly_make_monic_inplace(fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)

Makes the canonical nonzero $A$ monic by setting its content to the
inverse of the leading coefficient of the integer part.
:::

::: function
int _fmpq_mpoly_rescale(fmpq_t Acontent, fmpz * Acoeff, const fmpq_mpoly_t B, const fmpq * scales, const fmpq_mpoly_ctx_t ctx)

Writes $B(s_0 x_0, \ldots, s_{n-1} x_{n-1})$, for the rationals $s_i$ in
`scales`, as `Acontent` times the integer coefficients `Acoeff` on the
monomials of $B$, by clearing the denominators $b_i^{d_i}$ with $d_i$ the
maximal exponent of $x_i$ and taking out $a_i^{f_i}$ with $f_i$ the
minimal exponent. Returns 0 if the powers would be too large. Used by the
composition functions.
:::

::: function
int fmpq_mpoly_repack_bits(fmpq_mpoly_t A, const fmpq_mpoly_t B, flint_bitcnt_t Abits, const fmpq_mpoly_ctx_t ctx)

Sets $A$ to $B$ with exponents packed into `Abits` bits. Returns 0 if
the exponents do not fit, in which case $A$ is unspecified.
:::

::: function
void mpoly_void_ring_init_fmpq_mpoly_ctx(mpoly_void_ring_t R, const fmpq_mpoly_ctx_t ctx)

Fills the generic ring interface `R` (used by the expression parser) with
the operations of this module.
:::

::: function
void fmpq_mpoly_remainder_test(const fmpq_mpoly_t r, const fmpq_mpoly_t g, const fmpq_mpoly_ctx_t ctx)

Test helper: throws if some term of $r$ is divisible by the leading
monomial of $g$, using `fmpz_mpoly_remainder_strongtest` on the integer
parts.
:::
