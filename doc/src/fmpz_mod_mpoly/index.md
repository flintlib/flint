---
title: fmpz_mod_mpoly.h
header: fmpz_mod_mpoly.h
description: Sparse multivariate polynomials over $\Z/n\Z$ for a modulus $n$ of any size, with heap arithmetic, GCDs, square roots, resultants and naive Gröbner bases.
category: Multivariate polynomials
order: 9
---

This module implements multivariate polynomials over $\Z/n\Z$ for an
arbitrary modulus $n \ge 1$ given as an `fmpz`. It has the same
representation and interface as `nmod_mpoly` (see [](nmod_mpoly),
[](mpoly) and [](algorithms/multivariate/representation)): an
`fmpz_mod_mpoly_t` is a sorted array of terms, each an `fmpz`
coefficient in $[0, n)$ and a packed exponent vector, and the number
of variables, the monomial ordering and the modulus are held in an
`fmpz_mod_mpoly_ctx_t`, whose modulus is an `fmpz_mod_ctx_t` (see
[](fmpz_mod)). The functions taking scalars accept them as `fmpz`, `ulong`
or `slong`, which are reduced modulo $n$.

Multiplication uses Johnson's heap algorithm [@Joh1974] or Kronecker
substitution into `fmpz_mod_poly`
([](algorithms/multivariate/multiplication)), and division the heap
algorithm of Monagan and Pearce [@MonPea2011]
([](algorithms/multivariate/division)). The GCD algorithms (Brown,
Zippel, Hensel lifting; see [](algorithms/multivariate/gcd) and
[](algorithms/multivariate/hensel-lifting)) work with evaluation points
in $\Z/p\Z$; when the modulus fits a word, GCDs, square roots and factorisations are computed by converting to
`nmod_mpoly`, which also handles small primes
([](algorithms/multivariate/modp-gcd)). Factorisation is in
`fmpz_mod_mpoly_factor.h` and rational functions in `fmpz_mod_mpoly_q.h`.
The module also has a naive implementation of Buchberger's algorithm for
Gröbner bases [@BW1993], described with an example in
[](fmpz_mod_mpoly/groebner).

## Types and conventions {#conventions}

::: type
fmpz_mod_mpoly_struct
fmpz_mod_mpoly_t

A multivariate polynomial over $\Z/n\Z$, with the fields `coeffs` (an
array of `fmpz`), `exps` (the packed exponent vectors), `length`, `bits`
(bits per exponent field), `coeffs_alloc` and `exps_alloc`. An
`fmpz_mod_mpoly_t` is an array of length one of `fmpz_mod_mpoly_struct`.
:::

::: type
fmpz_mod_mpoly_ctx_struct
fmpz_mod_mpoly_ctx_t

The context of a polynomial ring $(\Z/n\Z)[x_0, \ldots, x_{m-1}]$, with
the fields `minfo` (an `mpoly_ctx_t` holding the number of variables and
the ordering) and `ffinfo` (an `fmpz_mod_ctx_t` for the modulus).
:::

The conventions of [](nmod_mpoly#conventions) apply: all polynomials passed
to one function belong to the same context; variable 0 is the most
significant; inputs and outputs are canonical (terms sorted in strictly
decreasing order, nonzero coefficients reduced into $[0, n)$), except for
the container operations; outputs may alias inputs unless stated
otherwise. Ring operations work for any modulus, division needs an
invertible leading coefficient of the divisor, and GCDs, resultants,
square roots and factorisation assume that $n$ is prime, which is not
checked.

## Context objects {#context}

::: function
void fmpz_mod_mpoly_ctx_init(fmpz_mod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, const fmpz_t modulus)

Initialises `ctx` for polynomials in `nvars` variables with the monomial
ordering `ord` (`ORD_LEX`, `ORD_DEGLEX` or `ORD_DEGREVLEX`) and
coefficients modulo `modulus`, which must be positive.
:::

::: function
void fmpz_mod_mpoly_ctx_init_rand(fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state, slong max_nvars, const fmpz_t modulus)
void fmpz_mod_mpoly_ctx_init_rand_bits(fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state, slong max_nvars, flint_bitcnt_t max_bits)
void fmpz_mod_mpoly_ctx_init_rand_bits_prime(fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state, slong max_nvars, flint_bitcnt_t max_bits)

Initialise `ctx` with a random number of variables in $[0,$
`max_nvars`$]$ and a random ordering, and the given modulus, a random
modulus of at most `max_bits` bits, or a random prime of at most
`max_bits` bits. Used for testing.
:::

::: function
slong fmpz_mod_mpoly_ctx_nvars(const fmpz_mod_mpoly_ctx_t ctx)
ordering_t fmpz_mod_mpoly_ctx_ord(const fmpz_mod_mpoly_ctx_t ctx)

Return the number of variables and the ordering of `ctx`.
:::

::: function
const fmpz * fmpz_mod_mpoly_ctx_modulus(const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_ctx_get_modulus(fmpz_t m, const fmpz_mod_mpoly_ctx_t ctx)

Return a pointer to the modulus of `ctx` (owned by the context), or set
`m` to it.
:::

::: function
void fmpz_mod_mpoly_ctx_clear(fmpz_mod_mpoly_ctx_t ctx)

Frees the memory used by `ctx`.
:::

## Memory management {#memory}

::: function
void fmpz_mod_mpoly_init(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_init2(fmpz_mod_mpoly_t A, slong alloc, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_init3(fmpz_mod_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx)

Initialise `A` as the zero polynomial, allocating nothing, space for
`alloc` terms, or space for `alloc` terms with `bits` bits per exponent
field.
:::

::: function
void fmpz_mod_mpoly_clear(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Frees the memory used by `A`.
:::

::: function
void fmpz_mod_mpoly_fit_length(fmpz_mod_mpoly_t A, slong length, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_fit_length_fit_bits(fmpz_mod_mpoly_t A, slong len, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_fit_length_reset_bits(fmpz_mod_mpoly_t A, slong len, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx)

Ensure space for `length` terms; the second also repacks the exponents to
at least `bits` bits, the third sets the field width to `bits` and leaves
the exponent data undefined.
:::

::: function
void fmpz_mod_mpoly_truncate(fmpz_mod_mpoly_t A, slong newlen, const fmpz_mod_mpoly_ctx_t ctx)

Keeps only the `newlen` leading terms of `A` if it is longer.
:::

::: function
int fmpz_mod_mpoly_repack_bits(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_repack_bits_inplace(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to `B` (resp. modify `A`) with exponents repacked into `Abits`
bits per field; return 0 if some exponent does not fit, setting `A` to
zero.
:::

## Input and output {#io}

The array `x` gives the variable names, starting with variable 0; if it is
`NULL`, the variables are named `x1`, `x2`, and so on.

::: function
char * fmpz_mod_mpoly_get_str_pretty(const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_fprint_pretty(FILE * file, const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_print_pretty(const fmpz_mod_mpoly_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)

Return a string representing `A`, which the caller frees with
`flint_free`, or print it to `file` or `stdout`. The printing functions
return a positive value on success.
:::

::: function
int fmpz_mod_mpoly_set_str_pretty(fmpz_mod_mpoly_t A, const char * str, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)

Parses `str` as in `nmod_mpoly_set_str_pretty`: integers, the variables
in `x`, `+`, `-`, `*`, exact `/` and `^`. Returns 0 on success, and $-1$
(setting `A` to zero) on failure.
:::

## Basic manipulation {#basic}

::: function
void fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_is_gen(const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to the variable of index `var`; return whether `A` is that
variable, or, for `var` $< 0$, any variable.
:::

::: function
void fmpz_mod_mpoly_set(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_swap(fmpz_mod_mpoly_t A, fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_equal(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

Assignment, swap and equality test.
:::

## Constants {#constants}

::: function
int fmpz_mod_mpoly_is_fmpz(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_get_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Return whether `A` is constant, and set `c` to the value of the constant
`A`, throwing if `A` is not constant.
:::

::: function
void fmpz_mod_mpoly_set_fmpz_mod(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_ui(fmpz_mod_mpoly_t A, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_si(fmpz_mod_mpoly_t A, slong c, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to the constant $c$; `fmpz_mod_mpoly_set_fmpz_mod` requires $c$
to be reduced, the other functions reduce it.
:::

::: function
void fmpz_mod_mpoly_zero(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_one(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to 0 or 1.
:::

::: function
int fmpz_mod_mpoly_equal_fmpz(const fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_equal_ui(const fmpz_mod_mpoly_t A, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_equal_si(const fmpz_mod_mpoly_t A, slong c, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_is_zero(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_is_one(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Return whether `A` equals the constant $c \bmod n$, zero or one.
:::

## Degrees {#degrees}

The zero polynomial has degree $-1$ in every variable and total degree
$-1$.

::: function
int fmpz_mod_mpoly_degrees_fit_si(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_degrees_si(slong * degs, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Set `degs` to the degrees of `A` in all variables, and test whether they
fit an `slong`.
:::

::: function
slong fmpz_mod_mpoly_degree_si(const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)

The degree of `A` in the variable of index `var`.
:::

::: function
int fmpz_mod_mpoly_total_degree_fits_si(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
slong fmpz_mod_mpoly_total_degree_si(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_total_degree_fmpz(fmpz_t td, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

The total degree of `A`, and whether it fits an `slong`.
:::

::: function
void fmpz_mod_mpoly_used_vars(int * used, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Sets `used[i]` to a nonzero value if the variable of index $i$ appears in
`A`, and to zero otherwise.
:::

## Coefficients {#coefficients}

::: function
void fmpz_mod_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_coeff_fmpz_monomial(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_ctx_t ctx)

Get or set the coefficient of the monomial `M`, which must be a monomial
(throws otherwise).
:::

::: function
void fmpz_mod_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mod_mpoly_t A, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_coeff_fmpz_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_coeff_ui_fmpz(fmpz_mod_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_coeff_si_fmpz(fmpz_mod_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_coeff_fmpz_ui(fmpz_mod_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_coeff_ui_ui(fmpz_mod_mpoly_t A, ulong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_coeff_si_ui(fmpz_mod_mpoly_t A, slong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)

Get or set the coefficient of the monomial with exponent vector `exp`,
given as `fmpz` pointers or as `ulong`. Setting a coefficient to zero
removes the term. Each call costs time linear in the length of `A`.
:::

::: function
void fmpz_mod_mpoly_get_coeff_vars_ui(fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fmpz_mod_mpoly_ctx_t ctx)

Sets `C` to the coefficient of $\prod_i x_{v_i}^{e_i}$ in `A` viewed as a
polynomial in the distinct variables $x_{v_i}$, $v_i$ = `vars[i]`, $e_i$ =
`exps[i]`, $0 \le i <$ `length`, with coefficients in the other variables.
:::

::: function
fmpz * fmpz_mod_mpoly_leadcoeff(fmpz_mod_mpoly_t A)

Returns a pointer to the leading coefficient of `A`, which must be
nonzero.
:::

## Conversion to univariate polynomials {#univariate-poly}

::: function
int fmpz_mod_mpoly_is_fmpz_mod_poly(const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_get_fmpz_mod_poly(fmpz_mod_poly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_fmpz_mod_poly(fmpz_mod_mpoly_t A, const fmpz_mod_poly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)

Test whether `A` involves only the variable of index `var`; convert such a
polynomial to an `fmpz_mod_poly` (returning 0 if it involves another
variable); and set `A` to the univariate polynomial `B` in that variable.
:::

## Comparison {#comparison}

::: function
int fmpz_mod_mpoly_cmp(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

Returns 1, $-1$ or 0 according as `A` is greater than, less than or equal
to `B` in a total order on polynomials that agrees with the monomial
ordering on monomials.
:::

## Container operations {#container}

These functions access terms by their index $i$, $0 \le i <$ length, and
can produce non-canonical polynomials, which `fmpz_mod_mpoly_sort_terms`
and `fmpz_mod_mpoly_combine_like_terms` repair. An index out of range
throws.

::: function
int fmpz_mod_mpoly_is_canonical(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_assert_canonical(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Return whether `A` is canonical (sorted, nonzero reduced coefficients,
valid exponents), or throw if it is not.
:::

::: function
slong fmpz_mod_mpoly_length(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_resize(fmpz_mod_mpoly_t A, slong new_length, const fmpz_mod_mpoly_ctx_t ctx)

Return the number of terms, or change it to `new_length`, deleting terms
at the end or appending zero terms.
:::

::: function
void fmpz_mod_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_term_coeff_fmpz(fmpz_mod_mpoly_t A, slong i, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_term_coeff_ui(fmpz_mod_mpoly_t A, slong i, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_term_coeff_si(fmpz_mod_mpoly_t A, slong i, slong c, const fmpz_mod_mpoly_ctx_t ctx)

Get or set the coefficient of the term of index $i$ (reduced modulo $n$).
:::

::: function
int fmpz_mod_mpoly_term_exp_fits_ui(const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_term_exp_fits_si(const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_get_term_exp_si(slong * exp, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
ulong fmpz_mod_mpoly_get_term_var_exp_ui(const fmpz_mod_mpoly_t A, slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx)
slong fmpz_mod_mpoly_get_term_var_exp_si(const fmpz_mod_mpoly_t A, slong i, slong var, const fmpz_mod_mpoly_ctx_t ctx)

Test whether the exponents of the term of index $i$ fit a `ulong` or
`slong`, and get its exponent vector or the exponent of one variable; the
`ulong` and `slong` versions throw if a value does not fit.
:::

::: function
void fmpz_mod_mpoly_set_term_exp_fmpz(fmpz_mod_mpoly_t A, slong i, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_term_exp_ui(fmpz_mod_mpoly_t A, slong i, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)

Set the exponent vector of the term of index $i$.
:::

::: function
void fmpz_mod_mpoly_get_term(fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_get_term_monomial(fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)

Set `M` to the term of index $i$, or to its monomial with coefficient 1.
:::

::: function
void fmpz_mod_mpoly_push_term_fmpz_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_push_term_fmpz_ffmpz(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_push_term_ui_fmpz(fmpz_mod_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_push_term_ui_ffmpz(fmpz_mod_mpoly_t A, ulong c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_push_term_si_fmpz(fmpz_mod_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_push_term_si_ffmpz(fmpz_mod_mpoly_t A, slong c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_push_term_fmpz_ui(fmpz_mod_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_push_term_ui_ui(fmpz_mod_mpoly_t A, ulong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_push_term_si_ui(fmpz_mod_mpoly_t A, slong c, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)

Append a term with coefficient $c \bmod n$ and exponent vector `exp` in
constant amortised time. The exponents are an array of pointers to `fmpz`
(`_fmpz`), an array of `fmpz` (`_ffmpz`) or of `ulong` (`_ui`).
:::

::: function
void fmpz_mod_mpoly_sort_terms(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_combine_like_terms(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Sort the terms into the order of the context (without combining), and
combine adjacent equal monomials and remove zero coefficients; together
they make `A` canonical.
:::

## Random generation {#random}

::: function
void fmpz_mod_mpoly_randtest_bound(fmpz_mod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_randtest_bounds(fmpz_mod_mpoly_t A, flint_rand_t state, slong length, ulong * exp_bounds, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_randtest_bits(fmpz_mod_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t exp_bits, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to a random polynomial with at most `length` terms and exponents
below `exp_bound`, below `exp_bounds[i]` for variable $i$, or of at most
`exp_bits` bits.
:::

## Addition and subtraction {#addition}

::: function
void fmpz_mod_mpoly_add_fmpz_mod(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_add_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_add_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_add_si(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_sub_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_sub_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_sub_si(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong c, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to $B + c$ or $B - c$. `fmpz_mod_mpoly_add_fmpz_mod` requires $c$
to be reduced.
:::

::: function
void fmpz_mod_mpoly_add(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_sub(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to $B + C$ or $B - C$.
:::

## Scalar operations {#scalar}

::: function
void fmpz_mod_mpoly_neg(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

Sets `A` to $-B$.
:::

::: function
void fmpz_mod_mpoly_scalar_mul_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_scalar_mul_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_scalar_mul_si(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to $c B$. The last version requires $c$ to be reduced and
invertible.
:::

::: function
void fmpz_mod_mpoly_scalar_addmul_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_t d, const fmpz_mod_mpoly_ctx_t ctx)

Sets `A` to $B + d C$.
:::

::: function
void fmpz_mod_mpoly_make_monic(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

Sets `A` to `B` divided by its leading coefficient; throws if `B` is zero
or the leading coefficient is not invertible.
:::

## Differentiation {#derivative}

::: function
void fmpz_mod_mpoly_derivative(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)

Sets `A` to the derivative of `B` with respect to the variable of index
`var`.
:::

## Evaluation and composition {#evaluation}

::: function
void fmpz_mod_mpoly_evaluate_all_fmpz(fmpz_t eval, const fmpz_mod_mpoly_t A, fmpz * const * alphas, const fmpz_mod_mpoly_ctx_t ctx)

Sets `eval` to the value of `A` at the point with coordinates
`alphas[0]`, `alphas[1]`, ... (one per variable).
:::

::: function
void fmpz_mod_mpoly_evaluate_one_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_t val, const fmpz_mod_mpoly_ctx_t ctx)

Sets `A` to `B` with the variable of index `var` replaced by `val`.
:::

::: function
int fmpz_mod_mpoly_compose_fmpz_mod_mpoly(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz_mod_mpoly_struct * const * C, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)
int fmpz_mod_mpoly_compose_fmpz_mod_mpoly_horner(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz_mod_mpoly_struct * const * C, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)
int fmpz_mod_mpoly_compose_fmpz_mod_mpoly_geobucket(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz_mod_mpoly_struct * const * C, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)

Set `A` to `B` evaluated at the polynomials `C[i]` of the context
`ctxAC`, one per variable of `ctxB`. `A` must not alias `B` or any
`C[i]`. Return 1 for success and 0 if the computation would be too large.
The main function transforms exponents by a matrix when every `C[i]` is
zero or a monomial with coefficient 1, and otherwise uses a Horner scheme
or a geobucket sum.
:::

::: function
void fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const slong * c, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)

Sets `A` to `B` with the variable of index $i$ of `ctxB` replaced by the
variable of index `c[i]` of `ctxAC`, or by zero if `c[i]` is negative.
:::

## Multiplication {#multiplication}

::: function
void fmpz_mod_mpoly_mul(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)

Sets `A` to $B C$, by Kronecker substitution if the product is dense
enough and otherwise by Johnson's heap algorithm.
:::

::: function
void fmpz_mod_mpoly_mul_johnson(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)

Sets `A` to $B C$ by Johnson's heap algorithm [@Joh1974].
:::

::: function
int fmpz_mod_mpoly_mul_dense(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t C, const fmpz_mod_mpoly_ctx_t ctx)

Tries to set `A` to $B C$ by Kronecker substitution into `fmpz_mod_poly`;
returns 0 if the dense representation would be too large.
:::

## Powering {#powering}

::: function
int fmpz_mod_mpoly_pow_ui(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong k, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_pow_fmpz(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_t k, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to $B^k$. `fmpz_mod_mpoly_pow_fmpz` throws if $k < 0$ and
returns 0 if $k$ does not fit an `slong` and `B` has more than one term;
otherwise both return 1.
:::

::: function
void fmpz_mod_mpoly_pow_rmul(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, ulong k, const fmpz_mod_mpoly_ctx_t ctx)

Sets `A` to $B^k$ by repeated multiplication.
:::

## Division {#division}

The division functions need the leading coefficient of the divisor to be
invertible; they throw `FLINT_DIVZERO` for a zero divisor.

::: function
int fmpz_mod_mpoly_divides(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

If `B` divides `A`, sets `Q` to $A/B$ and returns 1; otherwise sets `Q`
to zero and returns 0.
:::

::: function
void fmpz_mod_mpoly_divexact(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

Sets `Q` to $A/B$, assuming that the division is exact.
:::

::: function
int fmpz_mod_mpoly_divides_monagan_pearce(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_divides_dense(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

Variants of `fmpz_mod_mpoly_divides`: the heap algorithm of Monagan and
Pearce [@MonPea2011], and Kronecker substitution, which returns $-1$ if the
dense representation would be too large.
:::

::: function
void fmpz_mod_mpoly_div(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_divrem(fmpz_mod_mpoly_t Q, fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

Set `Q` and `R` to the quotient and remainder of `A` divided by `B`:
$A = QB + R$ with no term of $R$ divisible by the leading monomial of
$B$.
:::

::: function
void fmpz_mod_mpoly_divrem_ideal(fmpz_mod_mpoly_struct ** Q, fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, fmpz_mod_mpoly_struct * const * B, slong len, const fmpz_mod_mpoly_ctx_t ctx)

Divides `A` by the `len` polynomials `B[i]`: $A = R + \sum_i Q_i B_i$
with no term of $R$ divisible by any leading monomial of a $B_i$.
:::

::: function
void fmpz_mod_mpoly_div_monagan_pearce(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_divrem_monagan_pearce(fmpz_mod_mpoly_t Q, fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_divrem_ideal_monagan_pearce(fmpz_mod_mpoly_struct ** Q, fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, fmpz_mod_mpoly_struct * const * B, slong len, const fmpz_mod_mpoly_ctx_t ctx)

Heap implementations of the division functions.
:::

## Square roots {#sqrt}

::: function
int fmpz_mod_mpoly_sqrt(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_sqrt_heap(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_is_square(const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

If `A` is a square, set `Q` to a square root and return 1, otherwise
return 0. The modulus must be prime. A word-size modulus is handled by
`nmod_mpoly_sqrt`.
:::

::: function
int fmpz_mod_mpoly_quadratic_root(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

If $Q^2 + AQ = B$ has a solution, sets `Q` to one and returns 1;
otherwise returns 0. The modulus must be prime.
:::

## Greatest common divisors {#gcd}

These functions assume that the modulus is prime. They return 1 for
success and 0 for failure, which only happens when the exponents need
more than one word per field even after deflation.

::: function
void fmpz_mod_mpoly_term_content(fmpz_mod_mpoly_t M, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Sets `M` to the GCD of the terms of `A`: zero, or the monomial with
coefficient 1 whose exponents are the minima over the terms of `A`.
:::

::: function
int fmpz_mod_mpoly_content_vars(fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_t A, slong * vars, slong vars_length, const fmpz_mod_mpoly_ctx_t ctx)

Sets `g` to the GCD of the coefficients of `A` viewed as a polynomial in
the variables `vars[i]`, and returns 1 for success.
:::

::: function
int fmpz_mod_mpoly_gcd(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_gcd_cofactors(fmpz_mod_mpoly_t G, fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

Set `G` to the monic GCD of `A` and `B` (zero if both are zero), and
`fmpz_mod_mpoly_gcd_cofactors` also $\bar{A} = A/G$ and $\bar{B} = B/G$.
For a word-size modulus the computation is done by `nmod_mpoly`.
:::

::: function
int fmpz_mod_mpoly_gcd_brown(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_gcd_hensel(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_gcd_zippel(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_gcd_zippel2(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_gcd_subresultant(fmpz_mod_mpoly_t G, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

Set `G` to the monic GCD by one algorithm: Brown's dense interpolation,
Hensel lifting, Zippel's sparse interpolation, sparse interpolation with
two main variables, or the subresultant pseudo-remainder sequence in one
variable. Return 0 if the algorithm cannot be applied.
:::

::: function
void fmpz_mod_mpoly_deflation(fmpz * shift, fmpz * stride, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_deflate(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_inflate(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mod_mpoly_ctx_t ctx)

As `nmod_mpoly_deflation`, `nmod_mpoly_deflate` and `nmod_mpoly_inflate`:
the shifts and strides of the exponents of each variable, and the
substitutions $e \mapsto (e - \mathit{shift})/\mathit{stride}$ and $e
\mapsto \mathit{shift} + \mathit{stride} \cdot e$.
:::

## Resultants and discriminants {#resultant}

::: function
int fmpz_mod_mpoly_resultant(fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_discriminant(fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_t A, slong var, const fmpz_mod_mpoly_ctx_t ctx)

Set `R` to the resultant of `A` and `B`, or the discriminant of `A`, with
respect to the variable of index `var`, by the subresultant algorithm of
Ducos [@Duc2000]. Return 1 for success. The modulus must be prime.
:::

## Univariate representation {#univar}

An `fmpz_mod_mpoly_univar_t` holds a polynomial in a main variable with
`fmpz_mod_mpoly_t` coefficients in the other variables, as in
[](nmod_mpoly#univar).

::: type
fmpz_mod_mpoly_univar_struct
fmpz_mod_mpoly_univar_t

The fields are `coeffs`, `exps` (`fmpz` exponents), `alloc` and
`length`.
:::

::: function
void fmpz_mod_mpoly_univar_init(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_univar_clear(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_univar_fit_length(fmpz_mod_mpoly_univar_t A, slong length, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_univar_swap(fmpz_mod_mpoly_univar_t A, fmpz_mod_mpoly_univar_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_univar_zero(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_univar_print_pretty(const fmpz_mod_mpoly_univar_t A, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_univar_assert_canonical(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)

Memory management, printing and consistency check.
:::

::: function
void fmpz_mod_mpoly_to_univar(fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_from_univar(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_univar_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)

Convert `B` to a polynomial in the variable of index `var`, and back
(the coefficients must not involve that variable).
:::

::: function
int fmpz_mod_mpoly_univar_degree_fits_si(const fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)
slong fmpz_mod_mpoly_univar_length(const fmpz_mod_mpoly_univar_t A, const fmpz_mod_mpoly_ctx_t ctx)
slong fmpz_mod_mpoly_univar_get_term_exp_si(fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_univar_get_term_coeff(fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_univar_swap_term_coeff(fmpz_mod_mpoly_t c, fmpz_mod_mpoly_univar_t A, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_univar_set_coeff_ui(fmpz_mod_mpoly_univar_t A, ulong e, const fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx)

Access the terms: whether the degree fits an `slong`, the number of terms,
the exponent and the coefficient of the term of index $i$, and setting the
coefficient of $x^e$.
:::

::: function
int fmpz_mod_mpoly_univar_pseudo_gcd(fmpz_mod_mpoly_univar_t Gx, const fmpz_mod_mpoly_univar_t Ax, const fmpz_mod_mpoly_univar_t Bx, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_univar_resultant(fmpz_mod_mpoly_t R, const fmpz_mod_mpoly_univar_t Ax, const fmpz_mod_mpoly_univar_t Bx, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_univar_discriminant(fmpz_mod_mpoly_t D, const fmpz_mod_mpoly_univar_t Fx, const fmpz_mod_mpoly_ctx_t ctx)

A GCD up to factors free of the main variable, the resultant and the
discriminant with respect to the main variable, by subresultants. Return
1 for success.
:::

## Vectors {#vectors}

::: type
fmpz_mod_mpoly_vec_struct
fmpz_mod_mpoly_vec_t

A vector of polynomials, with the fields `p` (the entries), `alloc` and
`length`.
:::

::: macro
fmpz_mod_mpoly_vec_entry(vec, i)

A pointer to the entry of index $i$ of `vec`.
:::

::: function
void fmpz_mod_mpoly_vec_init(fmpz_mod_mpoly_vec_t vec, slong len, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_vec_clear(fmpz_mod_mpoly_vec_t vec, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_vec_swap(fmpz_mod_mpoly_vec_t x, fmpz_mod_mpoly_vec_t y, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_vec_fit_length(fmpz_mod_mpoly_vec_t vec, slong len, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_vec_set_length(fmpz_mod_mpoly_vec_t vec, slong len, const fmpz_mod_mpoly_ctx_t ctx)

Initialise `vec` to `len` zero polynomials, free it, swap two vectors,
ensure room for `len` entries, and set the length (truncating or
appending zeros).
:::

::: function
void fmpz_mod_mpoly_vec_set(fmpz_mod_mpoly_vec_t dest, const fmpz_mod_mpoly_vec_t src, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_vec_append(fmpz_mod_mpoly_vec_t vec, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_ctx_t ctx)
slong fmpz_mod_mpoly_vec_insert_unique(fmpz_mod_mpoly_vec_t vec, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_ctx_t ctx)

Copy a vector, append a copy of `f`, and insert `f` unless it is already
present, returning its index.
:::

::: function
void fmpz_mod_mpoly_vec_print(const fmpz_mod_mpoly_vec_t F, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_vec_randtest_not_zero(fmpz_mod_mpoly_vec_t vec, flint_rand_t state, slong len, slong poly_len, ulong exp_bound, fmpz_mod_mpoly_ctx_t ctx)

Print a vector, and set it to `len` random nonzero polynomials with at
most `poly_len` terms and exponents below `exp_bound`.
:::

## Ideals and Gröbner bases {#groebner}

These functions work with ideals of $(\Z/n\Z)[x_0, \ldots, x_{m-1}]$ for
prime $n$ (leading coefficients are inverted), using monic polynomials as
normalised generators. They follow [@BW1993] and are meant for small
examples.

::: function
void fmpz_mod_mpoly_spoly(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_ctx_t ctx)

Sets `res` to the S-polynomial of `f` and `g`, computed after making `f`
and `g` monic: $\mathrm{lcm}/\mathrm{lm}(f) \cdot f/\mathrm{lc}(f) -
\mathrm{lcm}/\mathrm{lm}(g) \cdot g/\mathrm{lc}(g)$, where lcm is the
least common multiple of the leading monomials.
:::

::: function
void fmpz_mod_mpoly_reduction_monic_part(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_vec_t Iv, const fmpz_mod_mpoly_ctx_t ctx)

Sets `res` to the remainder of `f` on division by the polynomials in
`Iv` (with `fmpz_mod_mpoly_divrem_ideal`), made monic unless it is zero.
:::

::: function
int fmpz_mod_mpoly_vec_is_groebner(const fmpz_mod_mpoly_vec_t G, const fmpz_mod_mpoly_vec_t F, const fmpz_mod_mpoly_ctx_t ctx)

Returns whether `G` is a Gröbner basis of the ideal it generates (all
S-polynomials reduce to zero) and, if `F` is not `NULL`, every element of
`F` reduces to zero modulo `G`, so that the ideal of `F` is contained in
that of `G`. The converse inclusion is not checked.
:::

::: function
int fmpz_mod_mpoly_vec_is_autoreduced(const fmpz_mod_mpoly_vec_t G, const fmpz_mod_mpoly_ctx_t ctx)

Returns whether `G` is autoreduced (inter-reduced).
:::

::: function
void fmpz_mod_mpoly_vec_autoreduction_groebner(fmpz_mod_mpoly_vec_t H, const fmpz_mod_mpoly_vec_t G, const fmpz_mod_mpoly_ctx_t ctx)

Sets `H` to the autoreduction of the Gröbner basis `G`, which is the
reduced Gröbner basis of the ideal, unique up to the order of the
entries.
:::

::: function
void fmpz_mod_mpoly_vec_set_monic_unique(fmpz_mod_mpoly_vec_t G, const fmpz_mod_mpoly_vec_t F, const fmpz_mod_mpoly_ctx_t ctx)

Sets `G` to the nonzero entries of `F` made monic, without duplicates.
:::

::: function
void fmpz_mod_mpoly_buchberger_naive(fmpz_mod_mpoly_vec_t G, const fmpz_mod_mpoly_vec_t F, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_buchberger_naive_with_limits(fmpz_mod_mpoly_vec_t G, const fmpz_mod_mpoly_vec_t F, slong ideal_len_limit, slong poly_len_limit, const fmpz_mod_mpoly_ctx_t ctx)

Set `G` to a Gröbner basis of the ideal generated by `F`, by Buchberger's
algorithm in the form GROEBNERNEW2 of [@BW1993], which discards useless
pairs with the criteria of Gebauer and Möller [@GM1988]. The basis is not
reduced in general; `fmpz_mod_mpoly_vec_autoreduction_groebner` reduces
it. The second version stops if the basis gets more than
`ideal_len_limit` elements or a polynomial gets more than
`poly_len_limit` terms, and returns 0 in that
case (`G` then generates the ideal but need not be a Gröbner basis);
otherwise it returns 1.
:::

## Conversions {#conversions}

::: function
void fmpz_mod_mpoly_set_fmpz_mpoly(fmpz_mod_mpoly_t res, const fmpz_mpoly_t f, fmpz_mod_mpoly_ctx_t ctxm, fmpz_mpoly_ctx_t ctx)

Sets `res` to the reduction modulo $n$ of the integer polynomial `f` of
the context `ctx`, which must have the same number of variables and
ordering as `ctxm`.
:::

::: function
void fmpz_mod_mpoly_get_fmpz_mpoly(fmpz_mpoly_t res, const fmpz_mod_mpoly_t f, fmpz_mpoly_ctx_t ctx)

Sets `res` to `f` with its coefficients lifted to integers in $[0, n)$.
:::

## Internal functions {#internal}

::: function
void _fmpz_mod_mpoly_fit_length(fmpz ** coeffs, slong * coeffs_alloc, ulong ** exps, slong * exps_alloc, slong N, slong length)
void _fmpz_mod_mpoly_set_length(fmpz_mod_mpoly_t A, slong newlen, const fmpz_mod_mpoly_ctx_t ctx)

Ensure room for `length` coefficients and exponent vectors of `N` words in
the given arrays, and set the length field.
:::

::: function
void _fmpz_mod_mpoly_radix_sort1(fmpz *, ulong *, slong left, slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask)
void _fmpz_mod_mpoly_radix_sort(fmpz *, ulong *, slong left, slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask)
void _fmpz_mod_mpoly_push_exp_ffmpz(fmpz_mod_mpoly_t A, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_push_exp_pfmpz(fmpz_mod_mpoly_t A, fmpz * const * exp, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_push_exp_ui(fmpz_mod_mpoly_t A, const ulong * exp, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_set_coeff_fmpz_fmpz(fmpz_mod_mpoly_t A, const fmpz_t c, const fmpz * exp, const fmpz_mod_mpoly_ctx_t ctx)

Radix sort of the terms `left` to `right - 1` (coefficients and exponent
arrays given explicitly), appending a term with undefined coefficient, and
setting a coefficient given an array of `fmpz` exponents.
:::

::: function
void _fmpz_mod_mpoly_eval_all_fmpz_mod(fmpz_t eval, const fmpz * Acoeffs, const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, const fmpz * alphas, const mpoly_ctx_t mctx, const fmpz_mod_ctx_t fctx)
void _fmpz_mod_mpoly_compose_mat(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mat_t M, const fmpz_mod_mpoly_ctx_t ctxB, const fmpz_mod_mpoly_ctx_t ctxAC)

Evaluation of term arrays at a point, and substitution of monomials by a
matrix acting on the exponent vectors.
:::

::: function
void _fmpz_mod_mpoly_mul_johnson_maxfields(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, fmpz * maxBfields, const fmpz_mod_mpoly_t C, fmpz * maxCfields, const fmpz_mod_mpoly_ctx_t ctx)
int _fmpz_mod_mpoly_mul_dense_maxfields(fmpz_mod_mpoly_t P, const fmpz_mod_mpoly_t A, fmpz * maxAfields, const fmpz_mod_mpoly_t B, fmpz * maxBfields, const fmpz_mod_mpoly_ctx_t ctx)
int _fmpz_mod_mpoly_divides_dense_maxfields(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, fmpz * maxAfields, const fmpz_mod_mpoly_t B, fmpz * maxBfields, const fmpz_mod_mpoly_ctx_t ctx)
int _fmpz_mod_mpoly_divides_monagan_pearce_maxfields(fmpz_mod_mpoly_t Q, const fmpz_mod_mpoly_t A, fmpz * maxAfields, const fmpz_mod_mpoly_t B, fmpz * maxBfields, const fmpz_mod_mpoly_ctx_t ctx)

Multiplication and division methods given the maximum exponent fields of
the inputs.
:::

::: function
void _fmpz_mod_mpoly_init_dense_mock(fmpz_mod_poly_t D, const fmpz_mod_mpoly_t A, const slong * Adeg_bounds, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_clear_dense_mock(fmpz_mod_poly_t D)

Set up `D` as a dense univariate polynomial holding the Kronecker
substitution of `A` with the given degree bounds, and free it.
:::

::: function
void _fmpz_mod_mpoly_from_univar(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mod_mpoly_univar_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_set_fmpz_mod_poly(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz * Bcoeffs, slong Blen, slong var, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_to_fmpz_mod_poly_deflate(fmpz_mod_poly_t A, const fmpz_mod_mpoly_t B, slong var, const ulong * Bshift, const ulong * Bstride, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_from_fmpz_mod_poly_inflate(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mod_poly_t B, slong var, const ulong * Ashift, const ulong * Astride, const fmpz_mod_mpoly_ctx_t ctx)

Conversions to and from univariate representations with given exponent
bits, shifts and strides.
:::

::: function
void _fmpz_mod_mpoly_get_nmod_mpoly(nmod_mpoly_t nA, const nmod_mpoly_ctx_t nctx, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_set_nmod_mpoly(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, const nmod_mpoly_t nA, const nmod_mpoly_ctx_t nctx)

Convert between an `fmpz_mod_mpoly` with a word-size modulus and an
`nmod_mpoly` in a context with the same variables, ordering and modulus.
:::

::: function
void fmpz_mod_mpolyl_lead_coeff(fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_t A, slong num_vars, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyl_content(fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_t A, slong num_vars, const fmpz_mod_mpoly_ctx_t ctx)

For `A` viewed as a polynomial in its first `num_vars` variables (in a
lexicographic context): the leading coefficient, and the GCD of the
coefficients (returning 1 for success).
:::

::: function
void fmpz_mod_mpoly_to_mpolyl_perm_deflate(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t lctx, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride)
void fmpz_mod_mpoly_from_mpolyl_perm_inflate(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t ctx, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t lctx, const slong * perm, const ulong * shift, const ulong * stride)

Permute, deflate and convert to a lexicographic context `lctx`, and undo
this with `Abits`-bit exponents.
:::

::: function
void fmpz_mod_mpoly_monic_part(fmpz_mod_mpoly_t res, const fmpz_mod_mpoly_t f, const fmpz_mod_mpoly_ctx_t ctx)

Sets `res` to `f` divided by its leading coefficient, or to zero if `f`
is zero. Throws if the leading coefficient is not invertible. Not used
in FLINT.
:::

::: function
void mpoly_void_ring_init_fmpz_mod_mpoly_ctx(mpoly_void_ring_t R, const fmpz_mod_mpoly_ctx_t ctx)

Sets up the generic ring interface `R` used by the parser and the
subresultant code.
:::

::: function
void fmpz_mod_mpoly_remainder_strongtest(const fmpz_mod_mpoly_t r, const fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_ctx_t ctx)

Throws if some term of `r` is divisible by the leading monomial of `g`.
:::

::: type
fmpz_mod_mpoly_geobucket_struct
fmpz_mod_mpoly_geobucket_t

A geobucket [@Yan1998], as for `nmod_mpoly_geobucket_t`.
:::

::: function
void fmpz_mod_mpoly_geobucket_init(fmpz_mod_mpoly_geobucket_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_geobucket_clear(fmpz_mod_mpoly_geobucket_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_geobucket_fit_length(fmpz_mod_mpoly_geobucket_t B, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_geobucket_set(fmpz_mod_mpoly_geobucket_t B, fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_geobucket_add(fmpz_mod_mpoly_geobucket_t B, fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_geobucket_sub(fmpz_mod_mpoly_geobucket_t B, fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_geobucket_empty(fmpz_mod_mpoly_t p, fmpz_mod_mpoly_geobucket_t B, const fmpz_mod_mpoly_ctx_t ctx)

As the `nmod_mpoly_geobucket` functions: `set` swaps `p` into the
geobucket, `add` and `sub` add or subtract `p`, `empty` returns the sum in
`p`.
:::
