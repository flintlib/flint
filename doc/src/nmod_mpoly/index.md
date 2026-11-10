---
title: nmod_mpoly.h
header: nmod_mpoly.h
description: Sparse multivariate polynomials over $\Z/n\Z$ for a word-size modulus $n$, with heap, array and Kronecker multiplication, heap division, GCDs with extension-field fallbacks for small fields, square roots and resultants.
category: Multivariate polynomials
order: 7
---

This module implements multivariate polynomials over $\Z/n\Z$ for a
word-size modulus $n \ge 1$. A polynomial is an `nmod_mpoly_t`, stored in
sparse distributed form: an array of terms, each a coefficient (a `ulong`
reduced modulo $n$) and a packed exponent vector, sorted by a monomial
ordering ([](algorithms/multivariate/representation)). The number of
variables, the ordering and the modulus are held in a separate context
object `nmod_mpoly_ctx_t` that every function takes as
its last argument. The exponent machinery (orderings, packing of exponent
vectors into fields of `bits` bits, parsing) is shared with the other
multivariate modules and described in [](mpoly); the interface mirrors
`fmpz_mpoly`, and `fmpz_mod_mpoly`, `fq_nmod_mpoly` and `fq_zech_mpoly`
offer the same functions over other coefficient rings.

Multiplication chooses between Johnson's heap algorithm [@Joh1974], dense
arrays indexed by exponents, and Kronecker substitution into an
`nmod_poly`, and uses several threads for large inputs
([](algorithms/multivariate/multiplication)). Exact division
uses the heap algorithm of Monagan and Pearce [@MonPea2011]
([](algorithms/multivariate/division)). GCDs choose among
Brown's dense interpolation [@Bro1971], Zippel's sparse interpolation
[@Zip1979] and Hensel lifting ([](algorithms/multivariate/gcd),
[](algorithms/multivariate/hensel-lifting)), and fall back to images over extension
fields $\mathbb{F}_{p^k}$ when the prime field is too small for random
evaluation points; see [](algorithms/multivariate/modp-gcd). Resultants
are computed by evaluation and interpolation after a Kronecker substitution
([](algorithms/multivariate/modp-bivariate-resultant)). Factorisation is in
`nmod_mpoly_factor.h`.

The guide [](nmod_mpoly/getting-started) shows how to set up a context,
build, print and inspect polynomials.

::: summary algorithms/multivariate/overview :::

## Types and conventions {#conventions}

::: type
nmod_mpoly_struct
nmod_mpoly_t

A multivariate polynomial over $\Z/n\Z$. The structure has the fields
`coeffs` (the coefficients), `exps` (the packed exponent vectors, $N$
words per term, where $N$ depends on `bits` and the context), `length`
(the number of terms), `bits` (the number of bits per exponent field) and
`coeffs_alloc`, `exps_alloc` (the allocated sizes in words). An
`nmod_mpoly_t` is an array of length one of `nmod_mpoly_struct`, so it is
passed by reference.
:::

::: type
nmod_mpoly_ctx_struct
nmod_mpoly_ctx_t

The context of a polynomial ring $(\Z/n\Z)[x_0, \ldots, x_{m-1}]$: the
field `minfo`, an `mpoly_ctx_t` holding the number $m$ of variables and the
monomial ordering, and the field `mod`, an `nmod_t` for the modulus $n$
(see [](nmod)). An `nmod_mpoly_ctx_t` is an array of length one of
`nmod_mpoly_ctx_struct`.
:::

- Every polynomial belongs to the context it was initialised with, and all
  polynomials passed to one function (except where two contexts are
  explicitly taken) must belong to the same context. This is not checked.
- The variables are indexed $0, \ldots, m-1$; variable 0 is the most
  significant one in every ordering. Without names, they print as `x1`,
  `x2`, and so on.
- A polynomial is in *canonical form* if its terms are sorted in strictly
  decreasing order of their monomials, no coefficient is zero, and every
  coefficient is reduced, that is, lies in $[0, n)$. All functions require
  canonical inputs and produce canonical outputs, except the container
  operations, which can break canonical form and have functions to repair
  it ([](#container)).
- Exponents are unbounded. They are packed into fields of `bits` bits,
  and the field width grows automatically when a result needs it; the
  functions with exponents passed as `fmpz` work for any size, those with
  `ulong` or `slong` exponents require them to fit.
- Unless stated otherwise, the output may be the same object as any input.
  Several outputs of one function (quotient and remainder, or GCD and
  cofactors) must be distinct.
- Ring operations (addition, multiplication, powering, evaluation,
  composition) work for any modulus. Division needs an invertible leading
  coefficient of the divisor and throws `FLINT_IMPINV` otherwise. GCDs,
  resultants, discriminants, square roots and factorisation assume that $n$
  is prime, which is not checked.
- Functions that return an `int` success flag fail (return 0) only when an
  intermediate result would need unreasonably large exponents or memory, or,
  for GCDs, when the inputs cannot be repacked into one word per exponent
  field even after deflation. For inputs of reasonable size they succeed.
- Multiplication, exact division and GCDs use the threads set with
  `flint_set_num_threads` when the inputs are large enough.

## Context objects {#context}

::: function
void nmod_mpoly_ctx_init(nmod_mpoly_ctx_t ctx, slong nvars, const ordering_t ord, ulong modulus)

Initialises `ctx` for polynomials in `nvars` variables with the monomial
ordering `ord`, one of `ORD_LEX`, `ORD_DEGLEX` and `ORD_DEGREVLEX` (see
[](mpoly)), and coefficients modulo $n = $ `modulus`. The modulus must be
at least 1.
:::

::: function
void nmod_mpoly_ctx_init_rand(nmod_mpoly_ctx_t ctx, flint_rand_t state, slong max_nvars, ulong modulus)

Initialises `ctx` with a random number of variables in $[0,$
`max_nvars`$]$ and a random ordering, and the given modulus. Used for
testing.
:::

::: function
void nmod_mpoly_ctx_set_modulus(nmod_mpoly_ctx_t ctx, ulong modulus)

Replaces the modulus of `ctx` by `modulus`, keeping the variables and the
ordering. Polynomials of the context must be reduced modulo the new
modulus before they are used with it.
:::

::: function
slong nmod_mpoly_ctx_nvars(const nmod_mpoly_ctx_t ctx)
ordering_t nmod_mpoly_ctx_ord(const nmod_mpoly_ctx_t ctx)
ulong nmod_mpoly_ctx_modulus(const nmod_mpoly_ctx_t ctx)

Return the number of variables, the monomial ordering and the modulus of
`ctx`.
:::

::: function
void nmod_mpoly_ctx_clear(nmod_mpoly_ctx_t ctx)

Frees the memory used by `ctx`. Polynomials of the context must be cleared
first.
:::

## Memory management {#memory}

::: function
void nmod_mpoly_init(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_init2(nmod_mpoly_t A, slong alloc, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_init3(nmod_mpoly_t A, slong alloc, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)

Initialise `A` as the zero polynomial. `nmod_mpoly_init` allocates
nothing; `nmod_mpoly_init2` allocates space for `alloc` terms with the
minimum exponent field width `MPOLY_MIN_BITS`, and `nmod_mpoly_init3`
space for `alloc` terms with `bits` bits per exponent field.
:::

::: function
void nmod_mpoly_clear(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Frees the memory used by `A`.
:::

::: function
void nmod_mpoly_fit_length(nmod_mpoly_t A, slong length, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_realloc(nmod_mpoly_t A, slong alloc, const nmod_mpoly_ctx_t ctx)

`nmod_mpoly_fit_length` ensures that `A` has space for at least `len`
terms, growing the allocation at least geometrically.
`nmod_mpoly_realloc` changes the allocation to exactly `alloc` terms,
which must be at least the current length; with `alloc` = 0 it frees the
memory and sets `A` to zero.
:::

::: function
void nmod_mpoly_fit_length_fit_bits(nmod_mpoly_t A, slong len, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_fit_length_reset_bits(nmod_mpoly_t A, slong len, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)

Ensure space for `len` terms with exponent fields of `bits` bits.
`nmod_mpoly_fit_length_fit_bits` repacks the existing exponents if `bits`
exceeds `A->bits`; `nmod_mpoly_fit_length_reset_bits` sets `A->bits` to
`bits` and leaves the exponent data undefined, for use before overwriting
all terms.
:::

::: function
void nmod_mpoly_truncate(nmod_mpoly_t A, slong newlen, const nmod_mpoly_ctx_t ctx)

Reduces the length of `A` to `newlen` if it is longer, keeping the
`newlen` leading terms. The result is canonical.
:::

::: function
int nmod_mpoly_repack_bits(nmod_mpoly_t A, const nmod_mpoly_t B, flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_repack_bits_inplace(nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx)

Set `A` to `B` (resp. modify `A`) with the exponents repacked into fields
of `Abits` bits. Return 1 for success and 0 if some exponent does not fit
into `Abits` bits, in which case `A` is set to zero.
:::

## Input and output {#io}

The array `x` gives the names of the variables, starting with variable 0.
If `x` is `NULL`, the variables are named `x1`, `x2`, and so on.

::: function
char * nmod_mpoly_get_str_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)

Returns a string representing `A`, such as `x^2*y+4*x+3`, which the caller
must free with `flint_free`. Terms appear in the order of the context and
coefficients as integers in $[0, n)$.
:::

::: function
int nmod_mpoly_fprint_pretty(FILE * file, const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_print_pretty(const nmod_mpoly_t A, const char ** x, const nmod_mpoly_ctx_t ctx)

Print `A` in the format of `nmod_mpoly_get_str_pretty` to `file` or to
`stdout`. The return value is positive on success and negative on failure.
:::

::: function
int nmod_mpoly_set_str_pretty(nmod_mpoly_t A, const char * str, const char ** x, const nmod_mpoly_ctx_t ctx)

Sets `A` to the polynomial in the null-terminated string `str`, written
with integers, the variable names `x`, parentheses and the operators `+`,
`-`, `*`, `/` and `^` (followed by a nonnegative integer exponent).
Integers are reduced modulo $n$. A division must be exact, which allows
division by invertible constants. Returns 0 on success; returns $-1$ and
sets `A` to zero if the string cannot be parsed or a division is not
exact.
:::

## Basic manipulation {#basic}

::: function
void nmod_mpoly_gen(nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)

Sets `A` to the variable of index `var`.
:::

::: function
int nmod_mpoly_is_gen(const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)

If `var` $\ge 0$, returns whether `A` is the variable of index `var`; if
`var` $< 0$, returns whether `A` is any variable. For $n = 1$ every
polynomial is zero, and so equal to every variable.
:::

::: function
void nmod_mpoly_set(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_swap(nmod_mpoly_t A, nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Set `A` to `B`, and swap `A` and `B` efficiently.
:::

::: function
int nmod_mpoly_equal(const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Returns whether `A` and `B` are equal.
:::

## Constants {#constants}

::: function
int nmod_mpoly_is_ui(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
ulong nmod_mpoly_get_ui(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

`nmod_mpoly_is_ui` returns whether `A` is a constant (including zero).
`nmod_mpoly_get_ui` returns the value of the constant `A`, and throws if
`A` is not constant.
:::

::: function
void nmod_mpoly_set_ui(nmod_mpoly_t A, ulong c, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_set_fmpz(nmod_mpoly_t A, const fmpz_t c, const nmod_mpoly_ctx_t ctx)

Set `A` to the constant $c \bmod n$.
:::

::: function
void nmod_mpoly_zero(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_one(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Set `A` to 0 or to 1.
:::

::: function
int nmod_mpoly_equal_ui(const nmod_mpoly_t A, ulong c, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_is_zero(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_is_one(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Return whether `A` equals the constant $c \bmod n$, zero or one.
:::

## Degrees {#degrees}

The degree of the zero polynomial is $-1$ in every variable, and so is its
total degree.

::: function
int nmod_mpoly_degrees_fit_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_degrees_si(slong * degs, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_degrees_fmpz(fmpz ** degs, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Set the array `degs` to the degrees of `A` in each variable;
`nmod_mpoly_degrees_fit_si` returns whether all of them fit an `slong`,
which `nmod_mpoly_degrees_si` requires.
:::

::: function
slong nmod_mpoly_degree_si(const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_degree_fmpz(fmpz_t deg, const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)

Return or set `deg` to the degree of `A` in the variable of index `var`.
The `slong` version requires the degree to fit.
:::

::: function
int nmod_mpoly_total_degree_fits_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_total_degree_si(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_total_degree_fmpz(fmpz_t td, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Return or set `td` to the total degree of `A`;
`nmod_mpoly_total_degree_fits_si` returns whether it fits an `slong`.
:::

::: function
void nmod_mpoly_used_vars(int * used, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Sets `used[i]` to a nonzero value if the variable of index $i$ appears in
`A` and to zero otherwise, for each variable.
:::

## Coefficients {#coefficients}

::: function
ulong nmod_mpoly_get_coeff_ui_monomial(const nmod_mpoly_t A, const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_set_coeff_ui_monomial(nmod_mpoly_t A, ulong c, const nmod_mpoly_t M, const nmod_mpoly_ctx_t ctx)

Get or set the coefficient of `A` at the monomial `M`, which must be a
monomial with any nonzero coefficient (otherwise the function throws).
Setting a coefficient to zero removes the term.
:::

::: function
ulong nmod_mpoly_get_coeff_ui_ui(const nmod_mpoly_t A, const ulong * exp, const nmod_mpoly_ctx_t ctx)
ulong nmod_mpoly_get_coeff_ui_fmpz(const nmod_mpoly_t A, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_set_coeff_ui_ui(nmod_mpoly_t A, ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_set_coeff_ui_fmpz(nmod_mpoly_t A, ulong c, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)

Get or set the coefficient of the monomial with exponent vector `exp`
(one entry per variable). Setting to zero removes the term, setting a new
monomial inserts it in its place. Each call costs time linear in the length
of `A`; use `nmod_mpoly_push_term_ui_ui` to build polynomials term by
term.
:::

::: function
ulong nmod_mpoly_get_term_ui_ui(const nmod_mpoly_t poly, const ulong * exp, const nmod_mpoly_ctx_t ctx)
ulong nmod_mpoly_get_term_ui_fmpz(const nmod_mpoly_t poly, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)

Return the coefficient of the monomial with exponent vector `exp`, as
`nmod_mpoly_get_coeff_ui_ui` and `nmod_mpoly_get_coeff_ui_fmpz` (older
names).
:::

::: function
void nmod_mpoly_get_coeff_vars_ui(nmod_mpoly_t C, const nmod_mpoly_t A, const slong * vars, const ulong * exps, slong length, const nmod_mpoly_ctx_t ctx)

Sets `C` to the coefficient of $\prod_i x_{v_i}^{e_i}$ in `A` viewed as a
polynomial in the variables $x_{v_i}$ with coefficients in the other
variables, where $v_i$ = `vars[i]` and $e_i$ = `exps[i]` for $0 \le i <$
`length`. The variables must be distinct, and $0 <$ `length` $\le m$. The
result does not involve the variables in `vars`.
:::

::: function
ulong nmod_mpoly_leadcoeff(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Returns the leading coefficient of `A`, which must be nonzero.
:::

## Conversion to univariate polynomials {#univariate-poly}

::: function
int nmod_mpoly_is_nmod_poly(const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)

Returns whether `A` involves no variable other than the one of index
`var`.
:::

::: function
int nmod_mpoly_get_nmod_poly(nmod_poly_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_get_n_poly(n_poly_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)

If `B` is a polynomial in the variable of index `var` only, set `A` to it
as a univariate polynomial and return 1; otherwise return 0. `A` must have
the same modulus as the context.
:::

::: function
void nmod_mpoly_set_nmod_poly(nmod_mpoly_t A, const nmod_poly_t B, slong var, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_set_n_poly_mod(nmod_mpoly_t A, const n_poly_t B, slong var, const nmod_mpoly_ctx_t ctx)

Set `A` to the univariate polynomial `B` in the variable of index `var`.
:::

## Comparison {#comparison}

::: function
int nmod_mpoly_cmp(const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Returns 1, $-1$ or 0 according as `A` is greater than, less than or equal
to `B` in some total order on polynomials that is fixed for the context.
For monomials it agrees with the monomial ordering of the context.
:::

## Container operations {#container}

These functions access the terms of `A` by their index $i$ in the sorted
array, $0 \le i <$ `nmod_mpoly_length(A)`, the leading term having index
0. They can produce non-canonical polynomials (zero coefficients, unsorted
or repeated monomials), which `nmod_mpoly_sort_terms` and
`nmod_mpoly_combine_like_terms` repair. An index $i \ge$ length throws an
exception; negative indices are not allowed.

::: function
ulong * nmod_mpoly_term_coeff_ref(nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

Returns a pointer to the coefficient of the term of index $i$.
:::

::: function
int nmod_mpoly_is_canonical(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_assert_canonical(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

`nmod_mpoly_is_canonical` returns whether `A` is in canonical form
(sorted, no zero coefficients, reduced coefficients, valid exponent
packing). `nmod_mpoly_assert_canonical` throws if it is not.
:::

::: function
slong nmod_mpoly_length(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Returns the number of terms of `A`, which for canonical `A` is the number
of nonzero coefficients.
:::

::: function
void nmod_mpoly_resize(nmod_mpoly_t A, slong new_length, const nmod_mpoly_ctx_t ctx)

Sets the length of `A` to `new_length`, deleting terms from the end or
appending terms with zero coefficient and zero exponent.
:::

::: function
ulong nmod_mpoly_get_term_coeff_ui(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_set_term_coeff_ui(nmod_mpoly_t A, slong i, ulong c, const nmod_mpoly_ctx_t ctx)

Get or set the coefficient of the term of index $i$. Setting a zero
coefficient makes `A` non-canonical.
:::

::: function
int nmod_mpoly_term_exp_fits_ui(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_term_exp_fits_si(const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

Return whether all exponents of the term of index $i$ fit a `ulong`
(resp. an `slong`).
:::

::: function
void nmod_mpoly_get_term_exp_fmpz(fmpz ** exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_get_term_exp_ui(ulong * exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_get_term_exp_si(slong * exp, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

Set `exp` to the exponent vector of the term of index $i$. The `ulong`
and `slong` versions throw if an exponent does not fit.
:::

::: function
ulong nmod_mpoly_get_term_var_exp_ui(const nmod_mpoly_t A, slong i, slong var, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_get_term_var_exp_si(const nmod_mpoly_t A, slong i, slong var, const nmod_mpoly_ctx_t ctx)

Return the exponent of the variable of index `var` in the term of index
$i$, throwing if it does not fit.
:::

::: function
void nmod_mpoly_set_term_exp_fmpz(nmod_mpoly_t A, slong i, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_set_term_exp_ui(nmod_mpoly_t A, slong i, const ulong * exp, const nmod_mpoly_ctx_t ctx)

Set the exponent vector of the term of index $i$, which usually makes `A`
unsorted.
:::

::: function
void nmod_mpoly_get_term(nmod_mpoly_t M, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_get_term_monomial(nmod_mpoly_t M, const nmod_mpoly_t A, slong i, const nmod_mpoly_ctx_t ctx)

Set `M` to the term of index $i$, or to its monomial with coefficient 1.
:::

::: function
void nmod_mpoly_push_term_ui_ui(nmod_mpoly_t A, ulong c, const ulong * exp, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_push_term_ui_fmpz(nmod_mpoly_t A, ulong c, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_push_term_ui_ffmpz(nmod_mpoly_t A, ulong c, const fmpz * exp, const nmod_mpoly_ctx_t ctx)

Append a term with coefficient $c \bmod n$ and exponent vector `exp` to
`A`, in constant amortised time. The `_ffmpz` version takes the exponents
as an array of `fmpz`, the `_fmpz` version as an array of pointers to
`fmpz`. After all terms are pushed, call `nmod_mpoly_sort_terms` and
`nmod_mpoly_combine_like_terms` unless they were pushed in decreasing
order with distinct monomials and nonzero coefficients.
:::

::: function
void nmod_mpoly_sort_terms(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Sorts the terms of `A` into the order of the context, without combining
equal monomials or removing zero coefficients, by a radix sort in time
linear in the size of `A`.
:::

::: function
void nmod_mpoly_combine_like_terms(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Combines adjacent terms with equal monomials and removes zero
coefficients. If `A` is sorted, the result is canonical.
:::

::: function
void nmod_mpoly_reverse(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Sets `A` to `B` with the order of its terms reversed, which is not
canonical unless `B` has at most one term.
:::

## Random generation {#random}

::: function
void nmod_mpoly_randtest_bound(nmod_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_randtest_bounds(nmod_mpoly_t A, flint_rand_t state, slong length, ulong * exp_bounds, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_randtest_bits(nmod_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t exp_bits, const nmod_mpoly_ctx_t ctx)

Set `A` to a random polynomial with at most `length` terms and random
coefficients, with the exponents of each variable in $[0,$
`exp_bound`$)$, of the variable of index $i$ in $[0,$
`exp_bounds[i]`$)$, or of at most `exp_bits` bits.
:::

## Addition and subtraction {#addition}

::: function
void nmod_mpoly_add_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_sub_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)

Set `A` to $B + c$ or $B - c$, with $c$ reduced modulo $n$.
:::

::: function
void nmod_mpoly_add(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_sub(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

Set `A` to $B + C$ or $B - C$.
:::

## Scalar operations {#scalar}

::: function
void nmod_mpoly_neg(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Sets `A` to $-B$.
:::

::: function
void nmod_mpoly_scalar_mul_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_scalar_mul_nmod_invertible(nmod_mpoly_t A, const nmod_mpoly_t B, ulong c, const nmod_mpoly_ctx_t ctx)

Set `A` to $c B$. `nmod_mpoly_scalar_mul_nmod_invertible` requires $c$ to
be reduced and invertible modulo $n$, so that no coefficient can vanish.
:::

::: function
void nmod_mpoly_scalar_addmul_ui(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, ulong d, const nmod_mpoly_ctx_t ctx)

Sets `A` to $B + d C$.
:::

::: function
void nmod_mpoly_make_monic(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Sets `A` to `B` divided by its leading coefficient. Throws if `B` is zero
or its leading coefficient is not invertible.
:::

## Differentiation {#derivative}

::: function
void nmod_mpoly_derivative(nmod_mpoly_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)

Sets `A` to the derivative of `B` with respect to the variable of index
`var`. In characteristic $p$ the terms whose exponent in that variable is
divisible by $p$ disappear.
:::

## Evaluation and composition {#evaluation}

::: function
ulong nmod_mpoly_evaluate_all_ui(const nmod_mpoly_t A, const ulong * vals, const nmod_mpoly_ctx_t ctx)

Returns the value of `A` at the point whose coordinates are the entries of
`vals`, one per variable, which must be reduced modulo $n$.
:::

::: function
void nmod_mpoly_evaluate_one_ui(nmod_mpoly_t A, const nmod_mpoly_t B, slong var, ulong val, const nmod_mpoly_ctx_t ctx)

Sets `A` to `B` with the variable of index `var` replaced by the constant
`val`.
:::

::: function
int nmod_mpoly_compose_nmod_poly(nmod_poly_t A, const nmod_mpoly_t B, nmod_poly_struct * const * C, const nmod_mpoly_ctx_t ctx)

Sets `A` to `B` evaluated at the univariate polynomials `C[0]`, `C[1]`,
..., one per variable of the context. Returns 1 for success and 0 if the
result would be too large to compute.
:::

::: function
int nmod_mpoly_compose_nmod_mpoly(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
int nmod_mpoly_compose_nmod_mpoly_horner(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)
int nmod_mpoly_compose_nmod_mpoly_geobucket(nmod_mpoly_t A, const nmod_mpoly_t B, nmod_mpoly_struct * const * C, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)

Set `A` to `B` evaluated at the polynomials `C[0]`, `C[1]`, ..., one per
variable of `ctxB`, where `A` and the `C[i]` belong to `ctxAC`. `A` must
not alias `B` or any `C[i]`. Return 1 for success and 0 if the result
would be too large to compute. If every `C[i]` is zero or a monomial with
coefficient 1, `nmod_mpoly_compose_nmod_mpoly` only transforms the
exponent vectors by a matrix; otherwise it uses a Horner scheme (if some
`C[i]` has more than one term) or a geobucket sum of the powered terms.
:::

::: function
void nmod_mpoly_compose_nmod_mpoly_gen(nmod_mpoly_t A, const nmod_mpoly_t B, const slong * c, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)

Sets `A` to `B` with the variable of index $i$ in `ctxB` replaced by the
variable of index `c[i]` in `ctxAC`, or by zero if `c[i]` is negative.
The array `c` has one entry per variable of `ctxB`, and nonnegative
entries must be valid variable indices of `ctxAC`.
:::

## Multiplication {#multiplication}

::: function
void nmod_mpoly_mul(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

Sets `A` to $B C$. For short inputs it uses `nmod_mpoly_mul_johnson`. For
larger inputs with one word per exponent it tries Kronecker substitution
(`nmod_mpoly_mul_dense`) if the product is dense enough, then the array
method (`nmod_mpoly_mul_array`) for 2 to 7 variables if the array is not
too large and sparse, and otherwise the heap method; the array and heap
methods use several threads when the inputs have more than about 512 terms
per thread.
:::

::: function
void nmod_mpoly_mul_johnson(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_mul_heap_threaded(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

Set `A` to $B C$ by Johnson's heap algorithm [@Joh1974], which merges the
products $b_i C$ with a heap of size $\min(\mathrm{len}(B),
\mathrm{len}(C))$ and costs $O(\mathrm{len}(B) \mathrm{len}(C) \log
\min(\mathrm{len}(B), \mathrm{len}(C)))$ monomial comparisons. The first
version uses one thread, the second the available threads.
:::

::: function
int nmod_mpoly_mul_array(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_mul_array_threaded(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

Try to set `A` to $B C$ by accumulating the products in a dense array
indexed by the exponents (by the exponents bounded by the degrees for
`ORD_LEX`, by the total degree for the degree orderings). Return 0 without
changing `A` if the array would be too large or the exponents need more
than one word. The first version uses one thread.
:::

::: function
int nmod_mpoly_mul_dense(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t C, const nmod_mpoly_ctx_t ctx)

Try to set `A` to $B C$ by Kronecker substitution into univariate
polynomials and `nmod_poly_mul`. Returns 0 if the dense representation
would be too large.
:::

## Powering {#powering}

::: function
int nmod_mpoly_pow_ui(nmod_mpoly_t A, const nmod_mpoly_t B, ulong k, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_pow_fmpz(nmod_mpoly_t A, const nmod_mpoly_t B, const fmpz_t k, const nmod_mpoly_ctx_t ctx)

Set `A` to $B^k$, with $B^0 = 1$ (zero if $n = 1$). `nmod_mpoly_pow_fmpz`
throws if $k < 0$, and returns 0 (failure) if $k$ does not fit an `slong`
and `B` has more than one term; otherwise both return 1. For a prime
modulus $p \le 99999$, `nmod_mpoly_pow_ui` writes $k = \sum_j k_j p^j$ and
uses $B^{p^j} = B(x_0^{p^j}, \ldots)$, so that only the powers $B^{k_j}$
with $k_j < p$ are computed by repeated multiplication.
:::

## Division {#division}

The division functions require the leading coefficient of the divisor to be
invertible modulo $n$; they need not assume $n$ prime.

::: function
int nmod_mpoly_divides(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

If `B` divides `A`, sets `Q` to the exact quotient $A/B$ and returns 1;
otherwise sets `Q` to zero and returns 0. Throws `FLINT_DIVZERO` if `B` is
zero and `A` is not (unless $n = 1$), and `FLINT_IMPINV` if the leading
coefficient of `B` is not invertible. For long `A` the heap division runs
in several threads.
:::

::: function
void nmod_mpoly_divexact(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Sets `Q` to $A/B$, assuming that the division is exact; throws if it is
not.
:::

::: function
int nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_divides_dense(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Variants of `nmod_mpoly_divides`. The first uses the heap algorithm of
Monagan and Pearce [@MonPea2011] in one thread, which stops as soon as a
remainder term cannot be divided. `nmod_mpoly_divides_dense` uses
Kronecker substitution and univariate division; it returns $-1$ if the
dense representation would be too large, and otherwise 0 or 1 as
`nmod_mpoly_divides`.
:::

::: function
int nmod_mpoly_divides_heap_threaded(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
int _nmod_mpoly_divides_heap_threaded_pool(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles)

As `nmod_mpoly_divides`, using the heap and the available threads; the
underscore version takes the thread handles. Only
declared when the configuration knows that the target is strongly ordered
(the macro `nmod_mpoly_divides_heap_threaded` is then defined); in that
case `nmod_mpoly_divides` already uses it for large inputs.
:::

::: function
void nmod_mpoly_div(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_divrem(nmod_mpoly_t Q, nmod_mpoly_t R, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Set `Q` and `R` to the quotient and the remainder of the division of `A`
by `B` with respect to the monomial ordering: $A = QB + R$ where no term of
$R$ is divisible by the leading monomial of $B$. `nmod_mpoly_div` only
computes $Q$; when the division is known to be exact it is faster than
`nmod_mpoly_divides`.
:::

::: function
void nmod_mpoly_divrem_ideal(nmod_mpoly_struct ** Q, nmod_mpoly_t R, const nmod_mpoly_t A, nmod_mpoly_struct * const * B, slong len, const nmod_mpoly_ctx_t ctx)

Divides `A` by the `len` polynomials `B[0]`, ..., `B[len-1]`: sets the
quotients `Q[i]` and the remainder `R` with $A = R + \sum_i Q_i B_i$, where
no term of $R$ is divisible by the leading monomial of any $B_i$. At each
step the leading remainder term is divided by the first $B_i$ whose leading
monomial divides it. The result depends on the order of the divisors
unless they form a Gröbner basis.
:::

::: function
void nmod_mpoly_div_monagan_pearce(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_divrem_monagan_pearce(nmod_mpoly_t q, nmod_mpoly_t r, const nmod_mpoly_t poly2, const nmod_mpoly_t poly3, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_divrem_ideal_monagan_pearce(nmod_mpoly_struct ** Q, nmod_mpoly_t R, const nmod_mpoly_t A, nmod_mpoly_struct * const * B, slong len, const nmod_mpoly_ctx_t ctx)

The heap implementations [@MonPea2011] of `nmod_mpoly_div`,
`nmod_mpoly_divrem` and `nmod_mpoly_divrem_ideal`, which those functions
call.
:::

## Square roots {#sqrt}

::: function
int nmod_mpoly_sqrt(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_sqrt_heap(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_is_square(const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

If $A = Q^2$ for some $Q$, set `Q` to such a square root and return 1;
otherwise return 0 (and set `Q` to zero). `nmod_mpoly_is_square` only
returns whether `A` is a square. The modulus must be prime. For odd $p$ the
root is computed term by term with a heap, as for division; for $p = 2$,
$A$ is a square exactly when all its exponents are even, and the root is
obtained by halving them. Even moduli other than 2 throw `FLINT_IMPINV`.
:::

::: function
int nmod_mpoly_quadratic_root(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

If $Q^2 + AQ = B$ has a solution, sets `Q` to a solution and returns 1;
otherwise returns 0. The modulus must be prime. For odd $p$ this completes
the square; in characteristic 2, where that is impossible, the equation is
solved term by term ([](algorithms/multivariate/modp-factoring#quadratic)).
:::

## Greatest common divisors {#gcd}

The GCD functions assume that the modulus is prime, which is not checked.
They return 1 for success and 0 for failure, which only happens when the
exponents of the inputs need more than one word per field even after
removing common powers and strides.

::: function
void nmod_mpoly_term_content(nmod_mpoly_t M, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Sets `M` to the GCD of the terms of `A`: zero if `A` is zero, and
otherwise the monomial with coefficient 1 whose exponent in each variable is
the minimum over the terms of `A`.
:::

::: function
int nmod_mpoly_content_vars(nmod_mpoly_t g, const nmod_mpoly_t A, slong * vars, slong vars_length, const nmod_mpoly_ctx_t ctx)

Sets `g` to the GCD of the coefficients of `A` viewed as a polynomial in
the `vars_length` variables `vars[i]` with coefficients in the other
variables, and returns 1 for success. The result does not involve the
variables in `vars`.
:::

::: function
int nmod_mpoly_gcd(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_gcd_cofactors(nmod_mpoly_t G, nmod_mpoly_t Abar, nmod_mpoly_t Bbar, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Set `G` to the monic GCD of `A` and `B`; the GCD of zero and zero is zero.
`nmod_mpoly_gcd_cofactors` also sets $\bar{A} = A/G$ and $\bar{B} = B/G$.
Return 1 for success and 0 for failure, in which case the outputs are
undefined. The algorithm is chosen from the shapes and estimated densities
of the inputs; over small fields the evaluation-based algorithms work with
images over extension fields.
:::

::: function
int nmod_mpoly_gcd_brown(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_gcd_hensel(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_gcd_zippel(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_gcd_zippel2(nmod_mpoly_t G, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)

Set `G` to the monic GCD of `A` and `B` by one algorithm: Brown's dense
interpolation, Hensel lifting, Zippel's sparse interpolation in one main
variable, or sparse interpolation with two main variables and the leading
coefficient known in advance. Return 0 if the algorithm cannot be applied
or fails. The trivial cases (monomials, missing variables, univariate
inputs) are handled first, as in `nmod_mpoly_gcd`.
:::

::: function
void nmod_mpoly_deflation(fmpz * shift, fmpz * stride, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

For each variable $v$, let $S_v$ be the set of its exponents in `A`. Sets
`shift[v]` to $\min S_v$ and `stride[v]` to $\gcd \{ e - \min S_v : e \in
S_v \}$. If `A` is zero, all shifts and strides are zero.
:::

::: function
void nmod_mpoly_deflate(nmod_mpoly_t A, const nmod_mpoly_t B, const fmpz * shift, const fmpz * stride, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_inflate(nmod_mpoly_t A, const nmod_mpoly_t B, const fmpz * shift, const fmpz * stride, const nmod_mpoly_ctx_t ctx)

`nmod_mpoly_deflate` replaces each exponent $e$ of the variable $v$ by $(e
-$ `shift[v]`$)/$`stride[v]` (which must be exact), with the convention
that a zero stride gives 0. `nmod_mpoly_inflate` applies $e \mapsto$
`shift[v]` $+$ `stride[v]` $e$. Shifts and strides must be nonnegative.
Deflating by the values of `nmod_mpoly_deflation` and inflating gives back
the polynomial.
:::

## Resultants and discriminants {#resultant}

::: function
int nmod_mpoly_resultant(nmod_mpoly_t R, const nmod_mpoly_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)

Sets `R` to the resultant of `A` and `B` with respect to the variable of
index `var`, that is, the determinant of their Sylvester matrix as
polynomials in that variable. Returns 1 for success. The modulus must be
prime. When the inputs are of sufficient size, the resultant is computed
by `nmod_mpoly_resultant_kronecker`; otherwise by the subresultant
algorithm of Ducos on the univariate forms.
:::

::: function
int nmod_mpoly_resultant_kronecker(nmod_mpoly_t R, const nmod_mpoly_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)

Tries to set `R` to the resultant of `A` and `B` with respect to the
variable of index `var` by packing the other variables into one by
Kronecker substitution, after removing monomial factors and deflating,
and evaluating that variable at many points and interpolating
([](algorithms/multivariate/modp-bivariate-resultant)). Returns 1 for
success, and 0, leaving `R` unchanged, if an input has exponents of more
than one word or degree 0 in the variable of index `var`, if the dense box
of the resultant has more than $2^{32}$ coefficients, or if the field has
too few elements for that many evaluation points. The modulus must be
prime.
:::

::: function
int nmod_mpoly_discriminant(nmod_mpoly_t R, const nmod_mpoly_t A, slong var, const nmod_mpoly_ctx_t ctx)

Sets `R` to the discriminant of `A` with respect to the variable of index
`var`, and returns 1 for success. The modulus must be prime.
:::

## Univariate representation {#univar}

An `nmod_mpoly_univar_t` holds a polynomial in one main variable $x_v$
with `nmod_mpoly_t` coefficients that do not involve $x_v$, as a sparse
list of (exponent, coefficient) pairs sorted by decreasing exponent. It
rewrites an element of $(\Z/n\Z)[x_0, \ldots, x_{m-1}]$ as an element of
$(\Z/n\Z)[x_0, \ldots, \hat{x}_v, \ldots, x_{m-1}][x_v]$ and back; the
coefficients still belong to the original context.

::: type
nmod_mpoly_univar_struct
nmod_mpoly_univar_t

The fields are `coeffs` (an array of `nmod_mpoly_struct`), `exps` (an
array of `fmpz` exponents), `alloc` and `length`.
:::

::: function
void nmod_mpoly_univar_init(nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_univar_clear(nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_univar_fit_length(nmod_mpoly_univar_t A, slong length, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_univar_swap(nmod_mpoly_univar_t A, nmod_mpoly_univar_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_univar_zero(nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)

Initialise, clear, allocate space for `length` terms, swap, and set to
zero.
:::

::: function
void nmod_mpoly_to_univar(nmod_mpoly_univar_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_from_univar(nmod_mpoly_t A, const nmod_mpoly_univar_t B, slong var, const nmod_mpoly_ctx_t ctx)

`nmod_mpoly_to_univar` sets `A` to `B` viewed as a polynomial in the
variable of index `var`. `nmod_mpoly_from_univar` converts back; its
result is undefined if the coefficients of `B` involve that variable.
:::

::: function
int nmod_mpoly_univar_degree_fits_si(const nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_univar_length(const nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_univar_get_term_exp_si(nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)

Return whether the degree in the main variable fits an `slong`, the
number of terms, and the exponent of the term of index $i$.
:::

::: function
void nmod_mpoly_univar_get_term_coeff(nmod_mpoly_t c, const nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_univar_swap_term_coeff(nmod_mpoly_t c, nmod_mpoly_univar_t A, slong i, const nmod_mpoly_ctx_t ctx)

Set `c` to (resp. swap `c` with) the coefficient of the term of index $i$.
:::

::: function
void nmod_mpoly_univar_set_coeff_ui(nmod_mpoly_univar_t A, ulong e, const nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)

Sets the coefficient of $x_v^e$ to `c`, inserting or removing a term as
needed.
:::

::: function
int nmod_mpoly_univar_pseudo_gcd(nmod_mpoly_univar_t Gx, const nmod_mpoly_univar_t Ax, const nmod_mpoly_univar_t Bx, const nmod_mpoly_ctx_t ctx)

Sets `Gx` to a GCD of `Ax` and `Bx` up to a factor that does not involve
the main variable, computed by the subresultant pseudo-remainder sequence
of Ducos. Returns 1 for success. The modulus must be prime.
:::

::: function
int nmod_mpoly_univar_resultant(nmod_mpoly_t R, const nmod_mpoly_univar_t Ax, const nmod_mpoly_univar_t Bx, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_univar_discriminant(nmod_mpoly_t D, const nmod_mpoly_univar_t Fx, const nmod_mpoly_ctx_t ctx)

Set `R` to the resultant of `Ax` and `Bx`, and `D` to the discriminant of
`Fx`, with respect to the main variable, by subresultants. Return 1 for
success. The modulus must be prime.
:::

::: function
void nmod_mpoly_univar_print_pretty(const nmod_mpoly_univar_t A, const char ** x, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_univar_assert_canonical(nmod_mpoly_univar_t A, const nmod_mpoly_ctx_t ctx)

Print `A` as a sum of terms `(c)*X^e`, the coefficients `c` with the
variable names `x` and the main variable as `X`; and throw if `A` is not
canonical.
:::

## Internal functions {#internal}

These functions are used by the multiplication, GCD and factorisation
code. Their interfaces may change without notice.

### Low-level arithmetic and sorting

::: function
void _nmod_mpoly_fit_length(ulong ** coeffs, slong * coeffs_alloc, ulong ** exps, slong * exps_alloc, slong N, slong length)

Ensures that the coefficient array `coeffs` (with allocation
`coeffs_alloc`) has room for `length` coefficients and the exponent array
`exps` (with allocation `exps_alloc`) room for `length` exponent vectors of
`N` words, growing them at least geometrically.
:::

::: function
void _nmod_mpoly_set_length(nmod_mpoly_t A, slong newlen, const nmod_mpoly_ctx_t ctx)

Sets the length field of `A` to `newlen`, which must not exceed the
allocation.
:::

::: function
slong _nmod_mpoly_add(ulong * coeff1, ulong * exp1, const ulong * coeff2, const ulong * exp2, slong len2, const ulong * coeff3, const ulong * exp3, slong len3, slong N, const ulong * cmpmask, nmod_t fctx)
slong _nmod_mpoly_sub(ulong * coeff1, ulong * exp1, const ulong * coeff2, const ulong * exp2, slong len2, const ulong * coeff3, const ulong * exp3, slong len3, slong N, const ulong * cmpmask, nmod_t fctx)

Merge the sorted term arrays of lengths `len2` and `len3` (with `N` words
per exponent and comparison mask `cmpmask`) into `coeff1`, `exp1`, adding
or subtracting coefficients modulo `fctx`, and return the length of the
result. The output must have room for `len2 + len3` terms.
:::

::: function
void _nmod_mpoly_radix_sort1(nmod_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask)
void _nmod_mpoly_radix_sort(nmod_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask)

Sort the terms `left` to `right - 1` of `A` by a most-significant-bit
radix sort starting at bit `pos`, for one-word exponents (with masks
`cmpmask` and `totalmask`) or exponents of `N` words.
:::

::: function
void _nmod_mpoly_push_exp_ui(nmod_mpoly_t A, const ulong * exp, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_push_exp_ffmpz(nmod_mpoly_t A, const fmpz * exp, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_push_exp_pfmpz(nmod_mpoly_t A, fmpz * const * exp, const nmod_mpoly_ctx_t ctx)

Append a term with exponent vector `exp` and undefined coefficient to `A`,
repacking `A` if the exponent needs more bits.
:::

::: function
void _nmod_mpoly_set_coeff_ui_fmpz(nmod_mpoly_t A, ulong c, const fmpz * exp, const nmod_mpoly_ctx_t ctx)
ulong _nmod_mpoly_get_term_ui_fmpz(const nmod_mpoly_t poly, const fmpz * exp, const nmod_mpoly_ctx_t ctx)

Versions of `nmod_mpoly_set_coeff_ui_fmpz` and `nmod_mpoly_get_term_ui_fmpz`
that take the exponents as an array of `fmpz`.
:::

::: function
slong _nmod_mpoly_mul_johnson(nmod_mpoly_t A, const ulong * coeff2, const ulong * exp2, slong len2, const ulong * coeff3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N, const ulong * cmpmask, nmod_t fctx)

Sets `A` to the product of the term arrays (`coeff2`, `exp2`, `len2`) and
(`coeff3`, `exp3`, `len3`), with `N` words of `bits`-bit fields per
exponent, by Johnson's heap algorithm, and returns the length of the
product.
:::

::: function
void _nmod_mpoly_mul_johnson_maxfields(nmod_mpoly_t A, const nmod_mpoly_t B, fmpz * maxBfields, const nmod_mpoly_t C, fmpz * maxCfields, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_mul_heap_threaded_pool_maxfields(nmod_mpoly_t A, const nmod_mpoly_t B, fmpz * maxBfields, const nmod_mpoly_t C, fmpz * maxCfields, const nmod_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles)
int _nmod_mpoly_mul_array_DEG(nmod_mpoly_t A, const nmod_mpoly_t B, fmpz * maxBfields, const nmod_mpoly_t C, fmpz * maxCfields, const nmod_mpoly_ctx_t ctx)
int _nmod_mpoly_mul_array_LEX(nmod_mpoly_t A, const nmod_mpoly_t B, fmpz * maxBfields, const nmod_mpoly_t C, fmpz * maxCfields, const nmod_mpoly_ctx_t ctx)
int _nmod_mpoly_mul_array_threaded_pool_DEG(nmod_mpoly_t A, const nmod_mpoly_t B, fmpz * maxBfields, const nmod_mpoly_t C, fmpz * maxCfields, const nmod_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles)
int _nmod_mpoly_mul_array_threaded_pool_LEX(nmod_mpoly_t A, const nmod_mpoly_t B, fmpz * maxBfields, const nmod_mpoly_t C, fmpz * maxCfields, const nmod_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles)
int _nmod_mpoly_mul_dense(nmod_mpoly_t P, const nmod_mpoly_t A, fmpz * maxAfields, const nmod_mpoly_t B, fmpz * maxBfields, const nmod_mpoly_ctx_t ctx)

The multiplication methods behind `nmod_mpoly_mul`, given the maximum
exponent fields of the inputs (computed with `mpoly_max_fields_fmpz`) and,
for the threaded versions, the thread handles. The array and dense
versions return 0 if they cannot be applied.
:::

::: function
void _nmod_mpoly_mul_array_chunked_LEX(nmod_mpoly_t P, const nmod_mpoly_t A, const nmod_mpoly_t B, const ulong * mults, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_mul_array_chunked_DEG(nmod_mpoly_t P, const nmod_mpoly_t A, const nmod_mpoly_t B, ulong degb, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_addmul_array1_ulong1(ulong * poly1, const ulong * poly2, const ulong * exp2, slong len2, const ulong * poly3, const ulong * exp3, slong len3)
void _nmod_mpoly_addmul_array1_ulong2(ulong * poly1, const ulong * poly2, const ulong * exp2, slong len2, const ulong * poly3, const ulong * exp3, slong len3)
void _nmod_mpoly_addmul_array1_ulong3(ulong * poly1, const ulong * poly2, const ulong * exp2, slong len2, const ulong * poly3, const ulong * exp3, slong len3)

Kernels of the array multiplication: the products are accumulated chunk by
chunk in arrays of one, two or three words per coefficient, which are
reduced modulo $n$ when a chunk is complete.
:::

::: function
void _nmod_mpoly_pow_rmul(nmod_mpoly_t A, const ulong * Bcoeffs, const ulong * Bexps, slong Blen, ulong k, slong N, const ulong * cmpmask, nmod_t mod, nmod_mpoly_t T)
void nmod_mpoly_pow_rmul(nmod_mpoly_t A, const nmod_mpoly_t B, ulong k, const nmod_mpoly_ctx_t ctx)

Set `A` to $B^k$ by repeated multiplication; the underscore version works
on a term array and uses `T` as scratch space.
:::

::: function
int _nmod_mpoly_divides_threaded_pool(nmod_mpoly_t Q, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles)
int _nmod_mpoly_divides_monagan_pearce(nmod_mpoly_t Q, const ulong * coeff2, const ulong * exp2, slong len2, const ulong * coeff3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N, const ulong * cmpmask, nmod_t fctx)

Exact division with the given thread handles, and the heap exact
division of term arrays used by `nmod_mpoly_divides_monagan_pearce`, which
returns 0 if the division is not exact.
:::

::: function
ulong _nmod_mpoly_eval_all_ui(const ulong * Acoeffs, const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, const ulong * alphas, const mpoly_ctx_t mctx, nmod_t mod)

Evaluates the polynomial given by the term arrays at the point `alphas`.
:::

::: function
int _ff_poly_pow_ui_is_not_feasible(slong length, ulong e)
int _ff_poly_pow_fmpz_is_not_feasible(slong length, const fmpz_t e)

Return whether raising a univariate polynomial of the given length to the
power $e$ would exceed the addressable memory; used by the composition
functions.
:::

::: function
void _nmod_mpoly_compose_mat(nmod_mpoly_t A, const nmod_mpoly_t B, const fmpz_mat_t M, const nmod_mpoly_ctx_t ctxB, const nmod_mpoly_ctx_t ctxAC)

Sets `A` to `B` with each exponent vector multiplied by the matrix `M`
(with one more row than the number of fields of `ctxAC`), which
implements substitutions of monomials for variables.
:::

::: function
void _nmod_mpoly_set_nmod_poly(nmod_mpoly_t A, flint_bitcnt_t Abits, const ulong * Bcoeffs, slong Blen, slong var, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_to_nmod_poly_deflate(nmod_poly_t A, const nmod_mpoly_t B, slong var, const ulong * Bshift, const ulong * Bstride, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_from_nmod_poly_inflate(nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_poly_t B, slong var, const ulong * Ashift, const ulong * Astride, const nmod_mpoly_ctx_t ctx)

Conversions between `nmod_mpoly` and univariate polynomials in the
variable of index `var`, with given exponent bits, shifts and strides.
:::

::: function
void _nmod_mpoly_from_univar(nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_mpoly_univar_t B, slong var, const nmod_mpoly_ctx_t ctx)

Version of `nmod_mpoly_from_univar` producing exponents of `Abits` bits.
:::

::: function
void nmod_mpoly_remainder_strongtest(const nmod_mpoly_t r, const nmod_mpoly_t g, const nmod_mpoly_ctx_t ctx)

Throws if some term of `r` is divisible by the leading monomial of `g`;
used in tests of `nmod_mpoly_divrem`.
:::

::: function
void mpoly_void_ring_init_nmod_mpoly_ctx(mpoly_void_ring_t R, const nmod_mpoly_ctx_t ctx)

Sets up the generic ring interface `R` (used by the parser and the
subresultant code in `mpoly`) for polynomials of the context.
:::

### Geobuckets

A geobucket [@Yan1998] holds a sum of polynomials as polynomials of
geometrically increasing lengths (powers of 4), so that adding many short
polynomials costs a logarithmic factor instead of a quadratic one.

::: type
nmod_mpoly_geobucket_struct
nmod_mpoly_geobucket_t

The fields are the arrays `polys` and `temps` of `FLINT_BITS/2`
polynomials and the number `length` of buckets in use.
:::

::: function
void nmod_mpoly_geobucket_init(nmod_mpoly_geobucket_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_geobucket_clear(nmod_mpoly_geobucket_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_geobucket_fit_length(nmod_mpoly_geobucket_t B, slong i, const nmod_mpoly_ctx_t ctx)

Initialise an empty geobucket, free it, and make room for `i` buckets.
:::

::: function
void nmod_mpoly_geobucket_set(nmod_mpoly_geobucket_t B, nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_geobucket_add(nmod_mpoly_geobucket_t B, nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_geobucket_sub(nmod_mpoly_geobucket_t B, nmod_mpoly_t p, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_geobucket_empty(nmod_mpoly_t p, nmod_mpoly_geobucket_t B, const nmod_mpoly_ctx_t ctx)

`nmod_mpoly_geobucket_set` sets `B` to `p`, swapping `p` into a bucket
(so `p` is changed). `nmod_mpoly_geobucket_add` and
`nmod_mpoly_geobucket_sub` add or subtract `p`, and
`nmod_mpoly_geobucket_empty` sets `p` to the sum and empties `B`.
:::

### Recursive and dense representations

The GCD algorithms work with polynomials in a few main variables whose
coefficients are polynomials in the others. An `nmod_mpolyu_t` is a sparse
univariate polynomial with `nmod_mpoly_t` coefficients, all with the same
exponent bits and lexicographic ordering. An `nmod_mpolyn_t` (defined in
`mpoly_types.h`) is a multivariate polynomial whose coefficients are
univariate `n_poly_t` polynomials in the last variable, and an
`nmod_mpolyun_t` a sparse univariate polynomial with `nmod_mpolyn_t`
coefficients. An `nmod_mpolyd_t` is a dense polynomial: with degree bounds
$(r_0, \ldots, r_{k-1})$, the coefficient of $x_0^{e_0} \cdots
x_{k-1}^{e_{k-1}}$ is stored at index $e_{k-1} + r_{k-1}(e_{k-2} +
r_{k-2}(\cdots))$.

::: type
nmod_mpolyu_struct
nmod_mpolyu_t
nmod_mpolyd_struct
nmod_mpolyd_t
nmod_mpolyd_ctx_struct
nmod_mpolyd_ctx_t

Sparse univariate polynomials with `nmod_mpoly_t` coefficients (fields
`coeffs`, `exps`, `alloc`, `length`, `bits`), dense polynomials (fields
`nvars`, `degb_alloc`, `deg_bounds`, `coeff_alloc`, `coeffs`), and the
variable permutation used with dense polynomials.
:::

::: function
void nmod_mpolyd_ctx_init(nmod_mpolyd_ctx_t dctx, slong nvars)
void nmod_mpolyd_ctx_clear(nmod_mpolyd_ctx_t dctx)
void nmod_mpolyd_init(nmod_mpolyd_t poly, slong nvars)
void nmod_mpolyd_clear(nmod_mpolyd_t poly)
void nmod_mpolyd_swap(nmod_mpolyd_t poly1, nmod_mpolyd_t poly2)
void nmod_mpolyd_fit_length(nmod_mpolyd_t poly, slong len)
void nmod_mpolyd_zero(nmod_mpolyd_t poly)
void nmod_mpolyd_set_nvars(nmod_mpolyd_t poly, slong nvars)
int nmod_mpolyd_set_degbounds(nmod_mpolyd_t A, slong * bounds)
int nmod_mpolyd_set_degbounds_perm(nmod_mpolyd_t A, const nmod_mpolyd_ctx_t dctx, slong * bounds)
slong nmod_mpolyd_length(const nmod_mpolyd_t A)
void nmod_mpolyd_print(nmod_mpolyd_t poly)

Memory management and basic operations of dense polynomials. The
`set_degbounds` functions set the degree bounds and allocate the
coefficient array, returning 0 if its size would overflow.
:::

::: function
void nmod_mpoly_convert_to_nmod_mpolyd(nmod_mpolyd_t A, const nmod_mpolyd_ctx_t dctx, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_convert_to_nmod_mpolyd_degbound(nmod_mpolyd_t A, const nmod_mpolyd_ctx_t dctx, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_convert_from_nmod_mpolyd(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, const nmod_mpolyd_t B, const nmod_mpolyd_ctx_t dctx)

Conversions between sparse and dense polynomials, used by
`nmod_mpoly_divides_dense`.
:::

::: function
void nmod_mpolyu_init(nmod_mpolyu_t A, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_clear(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
void nmod_mpolyu_swap(nmod_mpolyu_t A, nmod_mpolyu_t B, const nmod_mpoly_ctx_t uctx)
void nmod_mpolyu_zero(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
void nmod_mpolyu_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
int nmod_mpolyu_is_one(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx)
void nmod_mpolyu_fit_length(nmod_mpolyu_t A, slong length, const nmod_mpoly_ctx_t uctx)
void nmod_mpolyu_set(nmod_mpolyu_t A, const nmod_mpolyu_t B, const nmod_mpoly_ctx_t uctx)
void nmod_mpolyu_print_pretty(const nmod_mpolyu_t poly, const char ** x, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_degrees_si(slong * degs, const nmod_mpolyu_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_repack_bits_inplace(nmod_mpolyu_t A, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
nmod_mpoly_struct * _nmod_mpolyu_get_coeff(nmod_mpolyu_t A, ulong pow, const nmod_mpoly_ctx_t uctx)
ulong nmod_mpolyu_leadcoeff(nmod_mpolyu_t A, const nmod_mpoly_ctx_t ctx)

Memory management, assignment, printing and access for `nmod_mpolyu_t`.
`_nmod_mpolyu_get_coeff` returns the coefficient of $x^{\mathit{pow}}$,
inserting a zero term if necessary, and `nmod_mpolyu_leadcoeff` the
leading coefficient of the leading coefficient.
:::

::: function
void nmod_mpolyu_shift_right(nmod_mpolyu_t A, ulong s)
void nmod_mpolyu_shift_left(nmod_mpolyu_t A, ulong s)
void nmod_mpolyu_scalar_mul_nmod(nmod_mpolyu_t A, ulong c, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyu_content_mpoly(nmod_mpoly_t g, const nmod_mpolyu_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_divexact_mpoly_inplace(nmod_mpolyu_t A, nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_mul_mpoly(nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_mul_mpoly_inplace(nmod_mpolyu_t A, nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyuu_divides(nmod_mpolyu_t Q, const nmod_mpolyu_t A, const nmod_mpolyu_t B, slong nmainvars, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_setform(nmod_mpolyu_t A, nmod_mpolyu_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_setform_mpolyun(nmod_mpolyu_t A, nmod_mpolyun_t B, const nmod_mpoly_ctx_t ctx)

Arithmetic on `nmod_mpolyu_t`: shifting the main exponents by $s$,
scaling, the GCD of the coefficients (returning 1 for success), exact
division and multiplication by a polynomial in the other variables, exact
division with `nmainvars` main variables, and copying the monomial
structure ("form") of `B` into `A` with undefined coefficients, as needed
by sparse interpolation.
:::

::: function
void nmod_mpolyu_cvtto_poly(nmod_poly_t a, nmod_mpolyu_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_cvtfrom_poly(nmod_mpolyu_t A, nmod_poly_t a, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_cvtfrom_poly_notmain(nmod_mpolyu_t A, nmod_poly_t a, slong var, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_cvtto_mpolyun(nmod_mpolyun_t A, const nmod_mpolyu_t B, slong k, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu_cvtfrom_mpolyun(nmod_mpolyu_t A, const nmod_mpolyun_t B, slong var, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_to_mpolyu_perm_deflate_threaded_pool(nmod_mpolyu_t A, const nmod_mpoly_ctx_t uctx, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride, const thread_pool_handle * handles, slong num_handles)
void nmod_mpoly_from_mpolyu_perm_inflate(nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx, const nmod_mpolyu_t B, const nmod_mpoly_ctx_t uctx, const slong * perm, const ulong * shift, const ulong * stride)
void nmod_mpoly_to_mpolyl_perm_deflate(nmod_mpoly_t A, const nmod_mpoly_ctx_t lctx, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride)
void nmod_mpoly_from_mpolyl_perm_inflate(nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx, const nmod_mpoly_t B, const nmod_mpoly_ctx_t lctx, const slong * perm, const ulong * shift, const ulong * stride)
void nmod_mpoly_cvtto_mpolyn(nmod_mpolyn_t A, const nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_cvtfrom_mpolyn(nmod_mpoly_t A, const nmod_mpolyn_t B, slong var, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_to_mpolyn_perm_deflate_threaded_pool(nmod_mpolyn_t A, const nmod_mpoly_ctx_t nctx, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride, const thread_pool_handle * handles, slong num_handles)
void nmod_mpoly_from_mpolyn_perm_inflate(nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx, const nmod_mpolyn_t B, const nmod_mpoly_ctx_t nctx, const slong * perm, const ulong * shift, const ulong * stride)
void nmod_mpoly_from_mpolyun_perm_inflate(nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_mpoly_ctx_t ctx, const nmod_mpolyun_t B, const nmod_mpoly_ctx_t uctx, const slong * perm, const ulong * shift, const ulong * stride)

Conversions between the representations. The `perm_deflate` functions
permute the variables by `perm`, deflate their exponents by `shift` and
`stride` and convert into a context `uctx`, `lctx` or `nctx` with
lexicographic ordering; the `perm_inflate` functions undo this with
`Abits`-bit exponents.
:::

::: function
void nmod_mpolyn_init(nmod_mpolyn_t A, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_clear(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_swap(nmod_mpolyn_t A, nmod_mpolyn_t B)
void nmod_mpolyn_zero(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_one(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_is_zero(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_fit_length(nmod_mpolyn_t A, slong length, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_set_length(nmod_mpolyn_t A, slong newlen, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_fit_bits(nmod_mpolyn_t A, slong bits, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_set(nmod_mpolyn_t A, const nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_set_mod(nmod_mpolyn_t A, const nmod_t mod)
int nmod_mpolyn_is_canonical(const nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_print_pretty(const nmod_mpolyn_t A, const char ** x_in, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_equal(const nmod_mpolyn_t A, const nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_is_nonzero_nmod(const nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
ulong nmod_mpolyn_leadcoeff(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
n_poly_struct * nmod_mpolyn_leadcoeff_poly(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)
slong nmod_mpolyn_lastdeg(nmod_mpolyn_t A, const nmod_mpoly_ctx_t ctx)

Memory management, assignment and inspection of `nmod_mpolyn_t`.
`nmod_mpolyn_lastdeg` returns the degree in the last variable, which is
the variable of the `n_poly_t` coefficients.
:::

::: function
void nmod_mpolyn_scalar_mul_nmod(nmod_mpolyn_t A, ulong c, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_mul_poly(nmod_mpolyn_t A, const nmod_mpolyn_t B, const n_poly_t c, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_mul_last(nmod_mpolyn_t A, n_poly_t b, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_content_last(n_poly_t a, nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_divexact_last(nmod_mpolyn_t A, n_poly_t b, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_divides(nmod_mpolyn_t Q, const nmod_mpolyn_t A, const nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_divides_threaded_pool(nmod_mpolyn_t Q, const nmod_mpolyn_t A, const nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles)

Arithmetic on `nmod_mpolyn_t`: scaling, multiplication and exact division
by univariate polynomials in the last variable, the GCD of the `n_poly_t`
coefficients, and exact division (returning 1 if it is exact).
:::

::: function
void nmod_mpolyun_init(nmod_mpolyun_t A, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_clear(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_swap(nmod_mpolyun_t A, nmod_mpolyun_t B)
void nmod_mpolyun_zero(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_one(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_fit_length(nmod_mpolyun_t A, slong length, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_set(nmod_mpolyun_t A, const nmod_mpolyun_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_set_mod(nmod_mpolyun_t A, const nmod_t mod)
int nmod_mpolyun_is_canonical(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_print_pretty(const nmod_mpolyun_t poly, const char ** x, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyun_equal(const nmod_mpolyun_t A, const nmod_mpolyun_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyun_is_nonzero_nmod(const nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
ulong nmod_mpolyun_leadcoeff(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
n_poly_struct * nmod_mpolyun_leadcoeff_poly(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
slong nmod_mpolyun_lastdeg(nmod_mpolyun_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_shift_right(nmod_mpolyun_t A, ulong s)
void nmod_mpolyun_shift_left(nmod_mpolyun_t A, ulong s)
void nmod_mpolyun_scalar_mul_nmod(nmod_mpolyun_t A, ulong c, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_mul_poly(nmod_mpolyun_t A, const nmod_mpolyun_t B, const n_poly_t c, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_mul_last(nmod_mpolyun_t A, n_poly_t b, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_content_last(n_poly_t a, nmod_mpolyun_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_divexact_last(nmod_mpolyun_t A, n_poly_t b, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyun_divides(nmod_mpolyun_t Q, const nmod_mpolyun_t A, const nmod_mpolyun_t B, const nmod_mpoly_ctx_t ctx)

The same operations for `nmod_mpolyun_t`.
:::

::: function
void nmod_mpolyl_lead_coeff(nmod_mpoly_t c, const nmod_mpoly_t A, slong num_vars, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyl_content(nmod_mpoly_t g, const nmod_mpoly_t A, slong num_vars, const nmod_mpoly_ctx_t ctx)

For `A` in a lexicographic context, viewed as a polynomial in the first
`num_vars` variables: set `c` to its leading coefficient, and `g` to the
GCD of its coefficients (returning 1 for success).
:::

### Stacks of temporaries

::: type
nmod_poly_stack_struct
nmod_poly_stack_t

A stack of preallocated `n_poly_t`, `nmod_mpolyn_t` and `nmod_mpolyun_t`
objects for a fixed context and exponent bits, used to avoid repeated
allocation in the recursive GCD algorithms.
:::

::: function
void nmod_poly_stack_init(nmod_poly_stack_t S, flint_bitcnt_t bits, const nmod_mpoly_ctx_t ctx)
void nmod_poly_stack_clear(nmod_poly_stack_t S)
void nmod_poly_stack_set_ctx(nmod_poly_stack_t S, const nmod_mpoly_ctx_t ctx)
n_poly_struct ** nmod_poly_stack_fit_request_poly(nmod_poly_stack_t S, slong k)
nmod_mpolyn_struct ** nmod_poly_stack_fit_request_mpolyn(nmod_poly_stack_t S, slong k)
nmod_mpolyun_struct ** nmod_poly_stack_fit_request_mpolyun(nmod_poly_stack_t S, slong k)
n_poly_struct ** nmod_poly_stack_request_poly(nmod_poly_stack_t S, slong k)
nmod_mpolyn_struct ** nmod_poly_stack_request_mpolyn(nmod_poly_stack_t S, slong k)
nmod_mpolyun_struct ** nmod_poly_stack_request_mpolyun(nmod_poly_stack_t S, slong k)
n_poly_struct * nmod_poly_stack_take_top_poly(nmod_poly_stack_t S)
nmod_mpolyn_struct * nmod_poly_stack_take_top_mpolyn(nmod_poly_stack_t S)
nmod_mpolyun_struct * nmod_poly_stack_take_top_mpolyun(nmod_poly_stack_t S)
void nmod_poly_stack_give_back_poly(nmod_poly_stack_t S, slong k)
void nmod_poly_stack_give_back_mpolyn(nmod_poly_stack_t S, slong k)
void nmod_poly_stack_give_back_mpolyun(nmod_poly_stack_t S, slong k)
slong nmod_poly_stack_size_poly(const nmod_poly_stack_t S)
slong nmod_poly_stack_size_mpolyn(const nmod_poly_stack_t S)
slong nmod_poly_stack_size_mpolyun(const nmod_poly_stack_t S)

Initialise and clear a stack, and change its context. The `fit_request`
functions ensure that `k` more objects of a kind are available and return
a pointer to them, the `request` functions also reserve them,
`take_top` reserves one previously fitted object, `give_back` releases the
`k` most recent objects, and `size` returns the number of reserved
objects.
:::

::: function
n_poly_struct * evil_cast_nmod_poly_to_n_poly(nmod_poly_struct * a)
const n_poly_struct * evil_const_cast_nmod_poly_to_n_poly(const nmod_poly_struct * a)

Reinterpret an `nmod_poly_struct` as an `n_poly_struct`, whose first
three fields have the same layout.
:::

### GCD algorithms

These implement the methods selected by `nmod_mpoly_gcd`; see
[](algorithms/multivariate/modp-gcd). The suffix `smprime` marks versions
working over $\mathbb{F}_p$ itself, `medprime` versions that work over an
extension $\mathbb{F}_{p^k}$ represented with Zech logarithms, and
`lgprime` versions that reduce modulo irreducible polynomials in one
variable, that is, work over $\mathbb{F}_p[x]/(f)$, and combine the
images by Chinese remaindering in $\mathbb{F}_p[x]$.

::: function
int _nmod_mpoly_gcd_algo(nmod_mpoly_t G, nmod_mpoly_t Abar, nmod_mpoly_t Bbar, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx, unsigned int algo)
int _nmod_mpoly_gcd_algo_small(nmod_mpoly_t G, nmod_mpoly_t Abar, nmod_mpoly_t Bbar, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx, unsigned int algo)

Compute the GCD `G` and, when `Abar` and `Bbar` are not `NULL`, the
cofactors of the nonzero `A` and `B`, using the algorithms in the bit mask
`algo` (`MPOLY_GCD_USE_BROWN`, `MPOLY_GCD_USE_HENSEL`,
`MPOLY_GCD_USE_ZIPPEL`, `MPOLY_GCD_USE_ZIPPEL2`, or `MPOLY_GCD_USE_ALL`).
The `_small` version requires exponents of at most `FLINT_BITS` bits.
Return 1 for success.
:::

::: function
int nmod_mpolyn_gcd_brown_smprime(nmod_mpolyn_t G, nmod_mpolyn_t Abar, nmod_mpolyn_t Bbar, nmod_mpolyn_t A, nmod_mpolyn_t B, slong var, const nmod_mpoly_ctx_t ctx, const mpoly_gcd_info_t Iv, nmod_poly_stack_t Sp)
int nmod_mpolyn_gcd_brown_smprime_threaded_pool(nmod_mpolyn_t G, nmod_mpolyn_t Abar, nmod_mpolyn_t Bbar, nmod_mpolyn_t A, nmod_mpolyn_t B, slong var, const nmod_mpoly_ctx_t ctx, const mpoly_gcd_info_t Iv, const thread_pool_handle * handles, slong num_workers)
int nmod_mpolyn_gcd_brown_smprime_bivar(nmod_mpolyn_t G, nmod_mpolyn_t Abar, nmod_mpolyn_t Bbar, nmod_mpolyn_t A, nmod_mpolyn_t B, const nmod_mpoly_ctx_t ctx, nmod_poly_stack_t Sp)
int nmod_mpolyn_gcd_brown_lgprime(nmod_mpolyn_t G, nmod_mpolyn_t Abar, nmod_mpolyn_t Bbar, nmod_mpolyn_t A, nmod_mpolyn_t B, slong var, const nmod_mpoly_ctx_t ctx)

Brown's dense algorithm on `nmod_mpolyn_t` inputs: evaluation of the
variable of index `var` at points of $\mathbb{F}_p$ (`smprime`), or
reduction modulo irreducible polynomials of increasing degree in it
(`lgprime`), recursion, and interpolation with a degree-based termination
test. Return 1 for success, also setting the cofactors.
:::

::: function
int nmod_mpolyu_gcdm_zippel(nmod_mpolyu_t G, nmod_mpolyu_t Abar, nmod_mpolyu_t Bbar, nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpoly_ctx_t ctx, flint_rand_t randstate)
int nmod_mpolyu_gcdp_zippel(nmod_mpolyu_t G, nmod_mpolyu_t Abar, nmod_mpolyu_t Bbar, nmod_mpolyu_t A, nmod_mpolyu_t B, slong var, const nmod_mpoly_ctx_t ctx, flint_rand_t randstate)
nmod_gcds_ret_t nmod_mpolyu_gcds_zippel(nmod_mpolyu_t G, nmod_mpolyu_t A, nmod_mpolyu_t B, nmod_mpolyu_t f, slong var, const nmod_mpoly_ctx_t ctx, flint_rand_t randstate, slong * degbound)

Zippel's sparse interpolation on `nmod_mpolyu_t` inputs. `gcdp` works
with evaluation points in $\mathbb{F}_p$, `gcds` solves for the
coefficients of an assumed monomial structure `f` by evaluating
`var` variables and returns a status code (`nmod_gcds_ret_t`) that
distinguishes success from bad evaluation points and a wrong form, and
`gcdm` first tries `gcdp` and then falls back to images over
$\mathbb{F}_p[x_{m-1}]/(f)$ for irreducible $f$ of increasing degree.
:::

::: function
int nmod_mpolyl_gcd_zippel_smprime(nmod_mpoly_t rG, const slong * rGdegs, nmod_mpoly_t rAbar, nmod_mpoly_t rBbar, const nmod_mpoly_t A, const slong * Adegs, const nmod_mpoly_t B, const slong * Bdegs, const nmod_mpoly_t gamma, const slong * gammadegs, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyl_gcd_zippel_lgprime(nmod_mpoly_t rG, const slong * rGdegs, nmod_mpoly_t rAbar, nmod_mpoly_t rBbar, const nmod_mpoly_t A, const slong * Adegs, const nmod_mpoly_t B, const slong * Bdegs, const nmod_mpoly_t gamma, const slong * gammadegs, const nmod_mpoly_ctx_t ctx)

Sparse interpolation with two main variables (the algorithm of
`nmod_mpoly_gcd_zippel2`), given degree bounds and the leading coefficient
`gamma` of the GCD, over $\mathbb{F}_p$ or over extension fields.
:::

::: function
int nmod_mpolyl_gcd_hensel_smprime(nmod_mpoly_t G, slong Gdeg, nmod_mpoly_t Abar, nmod_mpoly_t Bbar, const nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyl_gcd_hensel_medprime(nmod_mpoly_t G, slong Gdeg, nmod_mpoly_t Abar, nmod_mpoly_t Bbar, const nmod_mpoly_t smA, const nmod_mpoly_t smB, const nmod_mpoly_ctx_t smctx)

GCD by Hensel lifting of a bivariate image factorisation $A = G \bar{A}$,
$B = G \bar{B}$ (with $G$ of degree `Gdeg` in the main variable), with
evaluation points in $\mathbb{F}_p$ or in an extension $\mathbb{F}_{p^k}$
with $p^k < 10^6$.
:::

::: function
int mpoly_gcd_get_use_first(slong rGdeg, slong Adeg, slong Bdeg, slong gammadeg)
int nmod_mpoly_gcd_get_use_new(slong rGdeg, slong Adeg, slong Bdeg, slong gammadeg, slong degxAB, slong degyAB, slong numABgamma, const n_polyun_t G, const n_polyun_t Abar, const n_polyun_t Bbar)

Heuristics of the two-variable sparse interpolation: decide which of the
GCD and the cofactors to interpolate, from their estimated degrees and
sizes.
:::

::: function
void _nmod_mpoly_monomial_evals_cache(n_poly_t E, const ulong * Aexps, flint_bitcnt_t Abits, slong Alen, const ulong * betas, slong start, slong stop, const mpoly_ctx_t mctx, nmod_t mod)
void _nmod_mpoly_monomial_evals2_cache(n_polyun_t E, const ulong * Aexps, flint_bitcnt_t Abits, slong Alen, const ulong * betas, slong m, const mpoly_ctx_t ctx, nmod_t mod)
void _nmod_poly_eval2_pow(ulong * vp, ulong * vm, n_poly_t P, n_poly_t alphapow, nmod_t fctx)

Evaluate all monomials of a polynomial at given values of some variables,
caching the powers; and evaluate a univariate polynomial at $\alpha$ and
$-\alpha$ simultaneously, given the powers of $\alpha$.
:::

::: function
void nmod_mpolyn_interp_reduce_2sm_poly(n_poly_t E, n_poly_t F, const nmod_mpolyn_t A, n_poly_t alphapow, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_interp_lift_2sm_poly(slong * lastdeg_, nmod_mpolyn_t F, const n_poly_t A, const n_poly_t B, ulong alpha, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_interp_crt_2sm_poly(slong * lastdeg_, nmod_mpolyn_t F, nmod_mpolyn_t T, const n_poly_t A, const n_poly_t B, const n_poly_t modulus, n_poly_t alphapow, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_interp_lift_sm_bpoly(nmod_mpolyn_t F, n_bpoly_t A, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_interp_crt_sm_bpoly(slong * lastdeg, nmod_mpolyn_t F, nmod_mpolyn_t T, n_bpoly_t A, n_poly_t modulus, n_poly_t alphapow, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_interp_reduce_2sm_mpolyn(nmod_mpolyn_t E, nmod_mpolyn_t F, nmod_mpolyn_t A, slong var, n_poly_t alphapow, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_interp_lift_2sm_mpolyn(slong * lastdeg, nmod_mpolyn_t T, nmod_mpolyn_t A, nmod_mpolyn_t B, slong var, ulong alpha, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_interp_crt_2sm_mpolyn(slong * lastdeg, nmod_mpolyn_t F, nmod_mpolyn_t T, nmod_mpolyn_t A, nmod_mpolyn_t B, slong var, n_poly_t modulus, n_poly_t alphapow, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_interp_reduce_sm_mpolyu(nmod_mpolyu_t B, nmod_mpolyun_t A, ulong alpha, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyn_interp_lift_sm_mpoly(nmod_mpolyn_t A, const nmod_mpoly_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyun_interp_lift_sm_mpolyu(nmod_mpolyun_t A, const nmod_mpolyu_t B, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_interp_crt_sm_mpoly(slong * lastdeg, nmod_mpolyn_t F, nmod_mpolyn_t T, nmod_mpoly_t A, n_poly_t modulus, ulong alpha, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyun_interp_crt_sm_mpolyu(slong * lastdeg, nmod_mpolyun_t F, nmod_mpolyun_t T, nmod_mpolyu_t A, n_poly_t modulus, ulong alpha, const nmod_mpoly_ctx_t ctx)
int nmod_mpolyn_interp_mcrt_sm_mpoly(slong * lastdeg_, nmod_mpolyn_t F, const nmod_mpoly_t A, const n_poly_t modulus, n_poly_t alphapow, const nmod_mpoly_ctx_t ctx)

Interpolation steps of Brown's algorithm: `reduce` evaluates the last
variable at a point $\alpha$ (the `2sm` versions at $\pm\alpha$ at once),
`lift` turns a first image into a polynomial in the last variable, and
`crt` combines the current interpolant `F` (modulo `modulus`) with a new
image by the Chinese remainder theorem, returning whether `F` changed and
updating `lastdeg`, the degree in the last variable.
:::

### Resultants

::: function
int n_bpoly_mod_resultant(n_poly_t res, const n_bpoly_t A, const n_bpoly_t B, nmod_t mod)
int _n_bpoly_mod_resultant(n_poly_t res, const n_poly_struct * A, slong lenA, const n_poly_struct * B, slong lenB, nmod_t mod)
int _n_bpoly_mod_resultant_len(n_poly_t res, const n_poly_struct * A, slong lenA, const n_poly_struct * B, slong lenB, slong rlen, nmod_t mod)
int _n_bpoly_mod_resultant_multipoint(n_poly_t res, const n_poly_struct * A, slong lenA, const n_poly_struct * B, slong lenB, slong rlen, nmod_t mod)
int n_bpoly_mod_resultant_cutoff(slong lenA, slong lenB, slong npoints)

`n_bpoly_mod_resultant` sets `res` to the resultant with respect to the
main variable $y$ of the bivariate polynomials `A` and `B` over
$\mathbb{F}_p$ (stored as polynomials in $y$ whose coefficients are
polynomials in $x$). It returns 0, leaving `res` unchanged, when no
algorithm applies, which happens when $\mathbb{F}_p$ does not contain
enough suitable points. `_n_bpoly_mod_resultant` is the same for
coefficient arrays with $\mathit{lenA} \ge \mathit{lenB} \ge 2$ and
nonzero leading coefficients, and `_n_bpoly_mod_resultant_len` is the same
again given that the resultant has length at most `rlen` in $x$ (or
`WORD_MAX`), which saves work when that is less than the bound read off the
degrees.

These choose the algorithm, by the sizes, and every dense bivariate
resultant over $\mathbb{F}_p$ in FLINT goes through them: those of
`nmod_mpoly_resultant` after Kronecker substitution, the images modulo
primes of `fmpz_bpoly_resultant` (and so of `fmpz_mpoly_resultant` and
`fmpq_mpoly_resultant`), and the resultants of `gr_poly` over polynomial
rings over `nmod`. An algorithm added to `_n_bpoly_mod_resultant_len` is
therefore used by all of these in its range. For now the only one is
`_n_bpoly_mod_resultant_multipoint`, evaluation at many values of $x$ and
interpolation (see [](algorithms/multivariate/modp-bivariate-resultant)),
which returns 0 when $\mathbb{F}_p$ does not contain enough suitable
points. The algorithms are threaded. `n_bpoly_mod_resultant_cutoff`
returns whether the best of them is expected to beat the subresultant
algorithm for the given lengths in $y$ and number of points.
:::

### Dense coefficient arrays

::: function
slong nmod_mpoly_append_array_sm1_LEX(nmod_mpoly_t P, slong Plen, ulong * coeff_array, const ulong * mults, slong num, slong array_size, slong top, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_append_array_sm2_LEX(nmod_mpoly_t P, slong Plen, ulong * coeff_array, const ulong * mults, slong num, slong array_size, slong top, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_append_array_sm3_LEX(nmod_mpoly_t P, slong Plen, ulong * coeff_array, const ulong * mults, slong num, slong array_size, slong top, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_append_array_sm1_DEGLEX(nmod_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_append_array_sm2_DEGLEX(nmod_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_append_array_sm3_DEGLEX(nmod_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_append_array_sm1_DEGREVLEX(nmod_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_append_array_sm2_DEGREVLEX(nmod_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb, const nmod_mpoly_ctx_t ctx)
slong nmod_mpoly_append_array_sm3_DEGREVLEX(nmod_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb, const nmod_mpoly_ctx_t ctx)

Append the entries of a dense array of unreduced coefficients (one chunk
of an array product: main variable or total degree `top`), reduced modulo
$n$ and omitting those that reduce to zero, to `P` in decreasing monomial
order, starting at term
index `Plen`, clear the array, and return the new length. `sm1`, `sm2` and
`sm3` arrays hold one, two or three words per coefficient.
:::
