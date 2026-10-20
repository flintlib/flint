---
title: fmpz_mpoly.h
header: fmpz_mpoly.h
description: Sparse multivariate polynomials over the integers, with heap, array and dense multiplication and division, square roots, GCDs, resultants and Gröbner bases.
category: Multivariate polynomials
order: 2
---

This module implements sparse multivariate polynomials in $\Z[x_0,
\ldots, x_{n-1}]$. A polynomial is an `fmpz_mpoly_t`, a list of terms
sorted by a monomial ordering fixed in a context object `fmpz_mpoly_ctx_t`,
which is passed to every function. Coefficients are `fmpz`, and exponents
are packed into machine words as described in the `mpoly` module; they
grow automatically and can be arbitrarily large. Multiplication chooses
between Johnson's heap method, dense array methods and Kronecker
substitution, and uses threads for large inputs. Exact and inexact division
use the heap algorithms of Monagan and Pearce. GCDs are computed by
modular algorithms (Brown's dense algorithm, Zippel's and the
Berlekamp–Massey based sparse algorithms) chosen by a cost model; Hensel
lifting and subresultants are available on request.
Factorisation is in `fmpz_mpoly_factor`, rational functions in
`fmpz_mpoly_q`, and polynomials over $\Q$ in `fmpq_mpoly`.

Start with the guide [](fmpz_mpoly/getting-started), which covers
contexts, orderings, building, parsing, printing and iterating over terms.

::: summary algorithms/multivariate/overview :::

## Types and conventions {#conventions}

::: type
fmpz_mpoly_struct
fmpz_mpoly_t

A multivariate polynomial over $\Z$. The fields are `coeffs` (an array of
`alloc` coefficients), `exps` (an array of `alloc` packed monomials of $N$
words each, where `N = mpoly_words_per_exp(bits, ctx->minfo)`), `alloc`,
`length` (the number of terms) and `bits` (the number of bits per packed
exponent field). Term $i$ is `coeffs[i]` times the monomial at `exps +
N*i`. An `fmpz_mpoly_t` is an array of one structure, so it is passed by
reference.
:::

::: type
fmpz_mpoly_ctx_struct
fmpz_mpoly_ctx_t

The context of a polynomial ring $\Z[x_0, \ldots, x_{n-1}]$: the number of
variables and the monomial ordering, stored in an `mpoly_ctx_t` called
`minfo`.
:::

- Every function takes the context as its last argument, and all
  polynomials passed to one function must belong to the same context
  (except for the composition functions, which take two contexts).
- The variables are indexed from 0, and the variable of index 0 is the
  most significant in the ordering. Exponent vectors are arrays indexed in
  the same way (see [](mpoly#conventions)).
- A polynomial is *canonical* if its coefficients are nonzero, its
  monomials are strictly decreasing in the ordering of the context, and
  its exponents are packed correctly with no field larger than allowed by
  `bits`. Every function expects canonical inputs and produces canonical
  outputs, except the container functions in [](#container), which can
  create noncanonical polynomials on purpose; `fmpz_mpoly_sort_terms` and
  `fmpz_mpoly_combine_like_terms` restore canonical form.
- The packing width `bits` of an output is chosen by the function and may
  differ from that of the inputs; it only grows when necessary.
- Unless stated otherwise, outputs may be the same objects as inputs.
  Several outputs of one function (quotient and remainder, GCD and
  cofactors) must be distinct.
- Functions that return an `int` and are documented as possibly failing
  return 1 on success and 0 on failure. Failure means that the result would
  be unreasonably large (an evaluation or power with a huge exponent) or
  that an algorithm with restricted applicability did not apply; the
  output is then unspecified unless stated otherwise.
- Many functions are thin dispatchers to one of several algorithms that
  are also exported (for example `fmpz_mpoly_mul_johnson`); the
  dispatchers are the ones to use.

## Context object {#context}

::: function
void fmpz_mpoly_ctx_init(fmpz_mpoly_ctx_t ctx, slong nvars, const ordering_t ord)

Initialises `ctx` for the ring $\Z[x_0, \ldots, x_{\mathit{nvars}-1}]$
with the monomial ordering `ord`, which is one of `ORD_LEX`, `ORD_DEGLEX`
and `ORD_DEGREVLEX`.
:::

::: function
void fmpz_mpoly_ctx_init_rand(fmpz_mpoly_ctx_t mctx, flint_rand_t state, slong max_nvars)

Initialises `mctx` with a random ordering and a random number of
variables between 0 and `max_nvars`, for testing.
:::

::: function
slong fmpz_mpoly_ctx_nvars(const fmpz_mpoly_ctx_t ctx)
ordering_t fmpz_mpoly_ctx_ord(const fmpz_mpoly_ctx_t ctx)

Return the number of variables and the ordering of `ctx`.
:::

::: function
void fmpz_mpoly_ctx_clear(fmpz_mpoly_ctx_t ctx)

Releases any memory used by `ctx`. Polynomials must be cleared before
their context.
:::

## Memory management {#memory}

::: function
void fmpz_mpoly_init(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Initialises `A` as the zero polynomial, without allocating memory, with
exponents packed into `MPOLY_MIN_BITS` bits.
:::

::: function
void fmpz_mpoly_init2(fmpz_mpoly_t A, slong alloc, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_init3(fmpz_mpoly_t A, slong alloc, flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)

Initialise `A` as zero with room for `alloc` terms, with exponents packed
into `MPOLY_MIN_BITS` bits, respectively into `bits` bits. `bits` must be
a valid packing width (see `mpoly_fix_bits`).
:::

::: function
void fmpz_mpoly_realloc(fmpz_mpoly_t A, slong alloc, const fmpz_mpoly_ctx_t ctx)

Changes the allocation of `A` to `alloc` terms. Terms beyond `alloc` are
discarded, which can leave a polynomial that is not the truncation of
interest but is still canonical. With `alloc` equal to zero, `A` is freed
and reinitialised.
:::

::: function
void fmpz_mpoly_fit_length(fmpz_mpoly_t A, slong len, const fmpz_mpoly_ctx_t ctx)

Ensures that `A` has room for at least `len` terms, at least doubling the
allocation when it grows, so that appending terms one by one costs
amortised constant time per term.
:::

::: function
void fmpz_mpoly_fit_bits(fmpz_mpoly_t A, flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)

Ensures that the exponents of `A` are packed with at least `bits` bits,
repacking the existing terms if necessary. `bits` must be a valid packing
width.
:::

::: function
void fmpz_mpoly_fit_length_reset_bits(fmpz_mpoly_t A, slong len, flint_bitcnt_t bits, const fmpz_mpoly_ctx_t ctx)

Ensures room for `len` terms with exponents packed into `bits` bits and
sets `A->bits` to `bits`. The value of `A` is destroyed (its monomials are
not repacked); the caller writes new terms and sets the length.
:::

::: function
void fmpz_mpoly_clear(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Frees the memory used by `A`.
:::

::: function
void _fmpz_mpoly_set_length(fmpz_mpoly_t A, slong newlen, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_truncate(fmpz_mpoly_t A, slong newlen, const fmpz_mpoly_ctx_t ctx)

Set the length of `A` to `newlen`, demoting the coefficients that are
dropped so that no memory is held by them. `_fmpz_mpoly_set_length`
requires `newlen` to be at most the allocation and may also increase the
length; `fmpz_mpoly_truncate` only decreases it, keeping the first
`newlen` terms.
:::

::: function
void _fmpz_mpoly_realloc(fmpz ** Acoeff, ulong ** Aexp, slong * Aalloc, slong len, slong N)
void _fmpz_mpoly_fit_length(fmpz ** Acoeff, ulong ** Aexp, slong * Aalloc, slong len, slong N)

Versions of `fmpz_mpoly_realloc` and `fmpz_mpoly_fit_length` for raw
arrays of coefficients and of monomials of $N$ words, with allocation
`*Aalloc`. New coefficients are zero.
:::

## Input and output {#io}

A variable name array `x` has one string per variable, the most
significant variable first. If `x` is `NULL`, the variables are called
`x1`, `x2`, and so on.

::: function
char * fmpz_mpoly_get_str_pretty(const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx)
char * _fmpz_mpoly_get_str_pretty(const fmpz * poly, const ulong * exps, slong len, const char ** x, slong bits, const mpoly_ctx_t mctx)

Return a string such as `3*x^2*y-y+1` representing `A` (respectively the
terms `(poly, exps, len)` packed with `bits` bits), terms in the order in
which they are stored. The string must be freed with `flint_free`.
:::

::: function
int fmpz_mpoly_print_pretty(const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_fprint_pretty(FILE * file, const fmpz_mpoly_t A, const char ** x, const fmpz_mpoly_ctx_t ctx)
int _fmpz_mpoly_print_pretty(const fmpz * poly, const ulong * exps, slong len, const char ** x, slong bits, const mpoly_ctx_t mctx)
int _fmpz_mpoly_fprint_pretty(FILE * file, const fmpz * poly, const ulong * exps, slong len, const char ** x_in, flint_bitcnt_t bits, const mpoly_ctx_t mctx)

Print `A` in the format of `fmpz_mpoly_get_str_pretty` to standard output
or to `file`. The return value is that of the underlying output functions
(positive on success).
:::

::: function
int fmpz_mpoly_set_str_pretty(fmpz_mpoly_t A, const char * str, const char ** x, const fmpz_mpoly_ctx_t ctx)

Sets `A` to the polynomial given by the expression `str` in the variables
named by `x`, and returns 0, or returns $-1$ if the string cannot be
parsed. The expression may contain nonnegative integer constants, the
variable names, `+`, `-` (binary or unary), `*`, `/`, `^`, parentheses and
spaces. Division must be exact, and `^` must be followed immediately by a
nonnegative integer exponent; powers associate to the left. Implicit
multiplication such as `2x` is not supported. Variable names are matched
greedily, longest first. On failure the value of `A` is unspecified.
:::

## Basic manipulation {#basic}

::: function
void fmpz_mpoly_gen(fmpz_mpoly_t poly, slong i, const fmpz_mpoly_ctx_t ctx)

Sets `poly` to the variable of index $i$.
:::

::: function
int fmpz_mpoly_is_gen(const fmpz_mpoly_t poly, slong k, const fmpz_mpoly_ctx_t ctx)

Returns whether `poly` is the variable of index $k$, or, if $k$ is
negative, whether it is any variable.
:::

::: function
void fmpz_mpoly_set(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_set(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong n, slong N)

Set `A` to `B`. The underscore version copies $n$ terms with monomials of
$N$ words between raw arrays.
:::

::: function
int fmpz_mpoly_equal(const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
int _fmpz_mpoly_equal(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong n, slong N)

Return whether `A` and `B` are equal. The polynomials may be packed with
different widths. The underscore version compares $n$ terms with
monomials of $N$ words.
:::

::: function
void fmpz_mpoly_swap(fmpz_mpoly_t A, fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Swaps `A` and `B` in constant time.
:::

::: function
slong fmpz_mpoly_max_bits(const fmpz_mpoly_t A)

Returns the maximum bit length $b$ of the absolute values of the
coefficients of `A`, as $b$ if all coefficients are nonnegative and as
$-b$ otherwise.
:::

::: function
int _fmpz_mpoly_fits_small(const fmpz * poly, slong len)

Returns whether every coefficient of the array is a small `fmpz` (not an
`mpz`), that is, has at most `SMALL_FMPZ_BITCOUNT_MAX` bits. Products of
two such coefficients fit in two words, and sums of up to $2^{\bits/4}$
or so of them in three, which the heap algorithms exploit.
:::

## Constants {#constants}

::: function
int fmpz_mpoly_is_fmpz(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Returns whether `A` is a constant (including zero).
:::

::: function
void fmpz_mpoly_get_fmpz(fmpz_t c, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Sets $c$ to the constant `A`. Throws `FLINT_ERROR` if `A` is not
constant.
:::

::: function
void fmpz_mpoly_set_fmpz(fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_ui(fmpz_mpoly_t A, ulong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_si(fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx)

Set `A` to the constant $c$.
:::

::: function
void fmpz_mpoly_zero(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_one(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Set `A` to 0, respectively 1.
:::

::: function
int fmpz_mpoly_equal_fmpz(const fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_equal_ui(const fmpz_mpoly_t A, ulong c, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_equal_si(const fmpz_mpoly_t A, slong c, const fmpz_mpoly_ctx_t ctx)

Return whether `A` equals the constant $c$.
:::

::: function
int fmpz_mpoly_is_zero(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_is_one(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Return whether `A` is 0, respectively 1.
:::

## Degrees {#degrees}

The degree of the zero polynomial is $-1$ in every variable, and so is its
total degree.

::: function
int fmpz_mpoly_degrees_fit_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Returns whether the degree of `A` in each variable fits an `slong`.
:::

::: function
void fmpz_mpoly_degrees_fmpz(fmpz ** degs, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_degrees_si(slong * degs, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Set `degs[i]` to the degree of `A` in the variable of index $i$. The
`_si` version requires the degrees to fit an `slong`.
:::

::: function
void fmpz_mpoly_degree_fmpz(fmpz_t deg, const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)
slong fmpz_mpoly_degree_si(const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)

Set or return the degree of `A` in the variable of index `var`.
:::

::: function
int fmpz_mpoly_total_degree_fits_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_total_degree_fmpz(fmpz_t td, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
slong fmpz_mpoly_total_degree_si(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Test whether the total degree of `A` fits an `slong`, and set or return
it. For degree orderings it is read from the first term; for lex all
terms are scanned.
:::

::: function
void fmpz_mpoly_used_vars(int * used, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Sets `used[i]` to a nonzero value if the variable of index $i$ occurs in
`A` and to zero otherwise.
:::

## Coefficients {#coefficients}

These functions access a coefficient by its monomial, given as an exponent
vector (`_ui`: `const ulong *`, `_fmpz`: `fmpz * const *`) or as a
monomial polynomial. A lookup is a binary search, $O(\log \mathit{len})$
monomial comparisons; inserting a new term or removing one moves the
following terms, $O(\mathit{len})$.

::: function
void fmpz_mpoly_get_coeff_fmpz_fmpz(fmpz_t c, const fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_get_coeff_fmpz_ui(fmpz_t c, const fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
ulong fmpz_mpoly_get_coeff_ui_fmpz(const fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
ulong fmpz_mpoly_get_coeff_ui_ui(const fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
slong fmpz_mpoly_get_coeff_si_fmpz(const fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
slong fmpz_mpoly_get_coeff_si_ui(const fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

Set or return the coefficient of the monomial with exponent vector `exp`
(zero if it does not occur). The `_ui` and `_si` return values are only
meaningful if the coefficient fits.
:::

::: function
void fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_coeff_fmpz_ui(fmpz_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_coeff_ui_fmpz(fmpz_mpoly_t A, const ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_coeff_ui_ui(fmpz_mpoly_t A, const ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_coeff_si_fmpz(fmpz_mpoly_t A, const slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_coeff_si_ui(fmpz_mpoly_t A, const slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_set_coeff_fmpz_fmpz(fmpz_mpoly_t A, const fmpz_t c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx)

Set the coefficient of the monomial with exponent vector `exp` to $c$,
inserting a term, changing one, or removing one if $c$ is zero, so that
`A` stays canonical. The packing width of `A` grows if needed. The
underscore version takes the exponents as a flat `fmpz` array.
:::

::: function
void fmpz_mpoly_get_coeff_fmpz_monomial(fmpz_t c, const fmpz_mpoly_t A, const fmpz_mpoly_t M, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_coeff_fmpz_monomial(fmpz_mpoly_t A, const fmpz_t c, const fmpz_mpoly_t M, const fmpz_mpoly_ctx_t ctx)

Get or set the coefficient of the monomial of the one-term polynomial
`M` (its coefficient is ignored). Throw `FLINT_ERROR` if `M` does not
have exactly one term.
:::

::: function
void fmpz_mpoly_get_coeff_vars_ui(fmpz_mpoly_t C, const fmpz_mpoly_t A, const slong * vars, const ulong * exps, slong length, const fmpz_mpoly_ctx_t ctx)

Sets $C$ to the coefficient of $\prod_k x_{\mathit{vars}[k]}^{\mathit{exps}[k]}$
in `A` viewed as a polynomial in the variables `vars` with coefficients
in the other variables: the sum of the terms of `A` whose exponents in
these variables are exactly `exps`, divided by that monomial. The arrays
have `length` entries, with $0 < \mathit{length} \le$ `nvars` and
distinct variables.
:::

## Comparison {#comparison}

::: function
int fmpz_mpoly_cmp(const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Returns $1$, $0$ or $-1$ according to an arbitrary but fixed total order
on polynomials: first by length, then by the monomials of the terms in
order, then by the coefficients. For two monomials with the same
coefficient this agrees with the monomial ordering.
:::

## Conversion to and from univariate polynomials {#conversion}

::: function
int fmpz_mpoly_is_fmpz_poly(const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)

Returns whether `A` involves no variable other than the one of index
`var`.
:::

::: function
int fmpz_mpoly_get_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)

Sets `A` to `B` as a univariate polynomial in the variable of index
`var`, and returns 1. Returns 0 if `B` involves another variable or a
degree does not fit an `slong`.
:::

::: function
void fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_t A, const fmpz_poly_t B, slong v, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_gen_fmpz_poly(fmpz_mpoly_t res, slong var, const fmpz_poly_t pol, const fmpz_mpoly_ctx_t ctx)

Set the multivariate polynomial to the univariate polynomial `B` (or
`pol`) in the variable of index `v` (or `var`). The second version throws
if the context has no variables.
:::

::: function
void _fmpz_mpoly_set_fmpz_poly(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz * Bcoeffs, slong Blen, slong var, const fmpz_mpoly_ctx_t ctx)

Sets `A` to $\sum_i \mathit{Bcoeffs}[i] x_{\mathit{var}}^i$ with
exponents packed into `Abits` bits, which must be large enough for the
degree $\mathit{Blen} - 1$.
:::

::: function
void _fmpz_mpoly_set_fmpz_poly_one_var(fmpz_mpoly_t A, flint_bitcnt_t Aminbits, fmpz * Acoeffs, slong Adeg, const fmpz_mpoly_ctx_t ctx)

For a context with exactly one variable, sets `A` to $\sum_{i \le
\mathit{Adeg}} \mathit{Acoeffs}[i] x^i$ with at least `Aminbits` bits per
field. The coefficients `Acoeffs` are moved into `A` and cleared, so the
array must not be used afterwards.
:::

## Container operations {#container}

These functions access terms by their index in the term list and can
produce noncanonical polynomials. A term index must lie in $[0,
\mathit{length})$; the functions that check it throw `FLINT_ERROR`
otherwise.

::: function
slong fmpz_mpoly_length(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Returns the number of terms of `A`, which for a canonical polynomial is
the number of nonzero coefficients.
:::

::: function
int fmpz_mpoly_is_canonical(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_assert_canonical(const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Test whether `A` is canonical: no coefficient is zero, the monomials are
strictly decreasing, and the packed exponents are valid (no guard bit
set, and for degree orderings the total degree field is correct). The
`assert` version throws `FLINT_ERROR` with a message describing the
problem instead of returning 0.
:::

::: function
fmpz * fmpz_mpoly_term_coeff_ref(fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)
fmpz * fmpz_mpoly_leadcoeff(const fmpz_mpoly_t A)

Return a pointer to the coefficient of term $i$, respectively of the
leading term (term 0) of the nonzero polynomial `A`. The index is only
checked in builds with assertions.
:::

::: function
void fmpz_mpoly_resize(fmpz_mpoly_t A, slong new_length, const fmpz_mpoly_ctx_t ctx)

Sets the length of `A` to `new_length`: terms are removed from the end,
or terms with zero coefficient and monomial 1 are appended (which makes
the polynomial noncanonical until they are filled in).
:::

::: function
void fmpz_mpoly_get_term_coeff_fmpz(fmpz_t c, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)
ulong fmpz_mpoly_get_term_coeff_ui(const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)
slong fmpz_mpoly_get_term_coeff_si(const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

Set or return the coefficient of term $i$. The `_ui` and `_si` results
are only meaningful if the coefficient fits.
:::

::: function
void fmpz_mpoly_set_term_coeff_fmpz(fmpz_mpoly_t A, slong i, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_term_coeff_ui(fmpz_mpoly_t A, slong i, ulong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_term_coeff_si(fmpz_mpoly_t A, slong i, slong c, const fmpz_mpoly_ctx_t ctx)

Set the coefficient of term $i$ to $c$. A zero coefficient is stored as
such.
:::

::: function
int fmpz_mpoly_term_exp_fits_si(const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_term_exp_fits_ui(const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

Return whether every exponent of term $i$ fits an `slong`, respectively a
`ulong`.
:::

::: function
void fmpz_mpoly_get_term_exp_fmpz(fmpz ** exp, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_get_term_exp_ui(ulong * exp, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_get_term_exp_si(slong * exp, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

Set `exp` to the exponent vector of term $i$. With exponents packed into
more than one word per field, the `_ui` and `_si` versions throw if an
exponent does not fit.
:::

::: function
ulong fmpz_mpoly_get_term_var_exp_ui(const fmpz_mpoly_t A, slong i, slong var, const fmpz_mpoly_ctx_t ctx)
slong fmpz_mpoly_get_term_var_exp_si(const fmpz_mpoly_t A, slong i, slong var, const fmpz_mpoly_ctx_t ctx)

Return the exponent of the variable of index `var` in term $i$. They
throw if the exponent does not fit the return type.
:::

::: function
void fmpz_mpoly_set_term_exp_fmpz(fmpz_mpoly_t A, slong i, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_term_exp_ui(fmpz_mpoly_t A, slong i, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

Set the monomial of term $i$ to the exponent vector `exp`, repacking `A`
if it does not fit. The terms may no longer be in order afterwards.
:::

::: function
void fmpz_mpoly_get_term(fmpz_mpoly_t M, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_get_term_monomial(fmpz_mpoly_t M, const fmpz_mpoly_t A, slong i, const fmpz_mpoly_ctx_t ctx)

Set `M` to term $i$ of `A`, respectively to its monomial with
coefficient 1.
:::

::: function
void fmpz_mpoly_push_term_fmpz_fmpz(fmpz_mpoly_t A, const fmpz_t c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_push_term_fmpz_ffmpz(fmpz_mpoly_t A, const fmpz_t c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_push_term_fmpz_ui(fmpz_mpoly_t A, const fmpz_t c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_push_term_ui_fmpz(fmpz_mpoly_t A, ulong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_push_term_ui_ffmpz(fmpz_mpoly_t A, ulong c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_push_term_ui_ui(fmpz_mpoly_t A, ulong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_push_term_si_fmpz(fmpz_mpoly_t A, slong c, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_push_term_si_ffmpz(fmpz_mpoly_t A, slong c, const fmpz * exp, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_push_term_si_ui(fmpz_mpoly_t A, slong c, const ulong * exp, const fmpz_mpoly_ctx_t ctx)

Append the term $c x^{\mathit{exp}}$ to `A` in amortised constant time,
repacking `A` if the exponents do not fit. Nothing else is done: the
result is canonical only if the terms are appended in decreasing order
with distinct monomials and nonzero coefficients. Otherwise call
`fmpz_mpoly_sort_terms` and `fmpz_mpoly_combine_like_terms` after the last
term.
:::

::: function
void _fmpz_mpoly_push_exp_ui(fmpz_mpoly_t A, const ulong * exp, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_push_exp_ffmpz(fmpz_mpoly_t A, const fmpz * exp, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_push_exp_pfmpz(fmpz_mpoly_t A, fmpz * const * exp, const fmpz_mpoly_ctx_t ctx)

Append a term with monomial `exp` and an unspecified coefficient, which
the caller then sets through `A->coeffs + A->length - 1`.
:::

::: function
void fmpz_mpoly_sort_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Sorts the terms of `A` into decreasing order of their monomials. Like
terms are not combined and zero coefficients are not removed. The sort is
a most significant bit first radix sort on the packed monomials, which
costs $O(\mathit{len} \cdot b)$ for $b$ significant bits of the monomials.
:::

::: function
void _fmpz_mpoly_radix_sort1(fmpz_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, ulong cmpmask, ulong totalmask)
void _fmpz_mpoly_radix_sort(fmpz_mpoly_t A, slong left, slong right, flint_bitcnt_t pos, slong N, ulong * cmpmask)

The recursive steps of `fmpz_mpoly_sort_terms`: sort the terms of index in
$[\mathit{left}, \mathit{right})$ by the bits of their monomials at
positions $\mathit{pos}, \mathit{pos} - 1, \ldots, 0$, assuming the higher
bits are equal. The first version is for one-word monomials and skips the
bit positions not set in `totalmask`.
:::

::: function
void fmpz_mpoly_combine_like_terms(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Adds up adjacent terms with equal monomials and removes terms with zero
coefficient, in linear time. If the terms are sorted, the result is
canonical.
:::

::: function
void fmpz_mpoly_reverse(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Sets `A` to `B` with the order of the terms reversed (which is not
canonical unless `B` has at most one term).
:::

## Random generation {#random}

::: function
void fmpz_mpoly_randtest_bound(fmpz_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t coeff_bits, ulong exp_bound, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_randtest_bounds(fmpz_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t coeff_bits, ulong * exp_bounds, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_randtest_bits(fmpz_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t coeff_bits, flint_bitcnt_t exp_bits, const fmpz_mpoly_ctx_t ctx)

Set `A` to a random polynomial with at most `length` terms (fewer after
like terms are combined), coefficients of about `coeff_bits` bits (as
from `fmpz_randtest`), and exponents in $[0, \mathit{exp\_bound})$,
in $[0, \mathit{exp\_bounds}[i])$ for the variable of index $i$, or such
that the monomials can be packed into `exp_bits` bits. Intended for
testing.
:::

## Addition and subtraction {#addition}

::: function
void fmpz_mpoly_add_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_add_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_add_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_sub_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_sub_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_sub_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)

Set `A` to $B + c$, respectively $B - c$. Only the constant term changes,
so with `A` aliased to `B` this costs constant time.
:::

::: function
void fmpz_mpoly_add(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_sub(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

Set `A` to $B + C$, respectively $B - C$, by merging the sorted term lists
in time linear in the total length.
:::

::: function
slong _fmpz_mpoly_add(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len2, const fmpz * poly3, const ulong * exps3, slong len3, slong N, const ulong * cmpmask)
slong _fmpz_mpoly_sub(fmpz * poly1, ulong * exps1, const fmpz * poly2, const ulong * exps2, slong len2, const fmpz * poly3, const ulong * exps3, slong len3, slong N, const ulong * cmpmask)

Merge two sorted term lists with monomials of $N$ words packed with the
same width into `(poly1, exps1)`, which must have room for $\mathit{len2}
+ \mathit{len3}$ terms and must not overlap the inputs, and return the
length of the result.
:::

## Scalar operations {#scalar}

::: function
void fmpz_mpoly_neg(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Sets `A` to $-B$.
:::

::: function
void fmpz_mpoly_scalar_mul_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_scalar_mul_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_scalar_mul_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)

Set `A` to $c B$.
:::

::: function
void fmpz_mpoly_scalar_fmma(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_t D, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)

Sets `A` to $c B + e D$ in one pass over both polynomials.
:::

::: function
void fmpz_mpoly_scalar_divexact_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_scalar_divexact_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_scalar_divexact_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)

Set `A` to $B/c$ for a nonzero $c$ that divides every coefficient of $B$.
This is not checked; if the division is not exact the result is
meaningless.
:::

::: function
int fmpz_mpoly_scalar_divides_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_scalar_divides_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong c, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_scalar_divides_si(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong c, const fmpz_mpoly_ctx_t ctx)

If the nonzero integer $c$ divides every coefficient of $B$, set `A` to
$B/c$ and return 1; otherwise set `A` to zero and return 0.
:::

## Differentiation and integration {#calculus}

::: function
void fmpz_mpoly_derivative(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)

Sets `A` to the partial derivative of `B` with respect to the variable of
index `var`.
:::

::: function
void fmpz_mpoly_integral(fmpz_mpoly_t A, fmpz_t scale, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)

Sets `A` and `scale` such that $A = \mathit{scale} \int B \, dx_v$ with
respect to the variable $x_v$ of index `var` (with zero constant of
integration), where $\mathit{scale} > 0$ is the least positive integer for
which this has integer coefficients.
:::

## Evaluation {#evaluation}

The functions that return an `int` return 0 instead of computing a result
that would need an unreasonable amount of memory, that is, when a value
of $b$ bits would have to be raised to a power $e$ with $b \cdot e$
larger than roughly `WORD_MAX`$/16$ (see `_fmpz_pow_ui_is_not_feasible`).

::: function
int fmpz_mpoly_evaluate_all_fmpz(fmpz_t ev, const fmpz_mpoly_t A, fmpz * const * vals, const fmpz_mpoly_ctx_t ctx)

Sets `ev` to $A(\mathit{vals}[0], \ldots, \mathit{vals}[n-1])$. The powers
of each value are computed by repeated squaring once, and each term is
evaluated with them. Returns 1 for success and 0 for failure.
:::

::: function
ulong fmpz_mpoly_evaluate_all_nmod(const fmpz_mpoly_t A, const ulong * alphas, const fmpz_mpoly_ctx_t ctx, nmod_t fpctx)
void fmpz_mpoly_evaluate_all_fmpz_mod(fmpz_t ev, const fmpz_mpoly_t A, const fmpz * alphas, const fmpz_mpoly_ctx_t ctx, const fmpz_mod_ctx_t fpctx)

Return or set the value of $A \bmod n$ at the point `alphas`, whose
entries must be reduced modulo the modulus $n$ of `fpctx`.
:::

::: function
void fmpz_mpoly_evaluate_acb(acb_t res, const fmpz_mpoly_t pol, acb_srcptr x, slong prec, const fmpz_mpoly_ctx_t ctx)

Sets `res` to the value of `pol` at the complex ball vector `x` (one
entry per variable) at precision `prec`, by direct term evaluation for at
most 6 terms and a Horner scheme otherwise; the result is indeterminate if
the evaluation fails.
:::

::: function
int fmpz_mpoly_evaluate_one_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_t val, const fmpz_mpoly_ctx_t ctx)

Sets `A` to `B` with the variable of index `var` replaced by `val`.
Returns 1 for success and 0 for failure.
:::

::: function
int fmpz_mpoly_compose_fmpz_poly(fmpz_poly_t A, const fmpz_mpoly_t B, fmpz_poly_struct * const * C, const fmpz_mpoly_ctx_t ctxB)

Sets `A` to $B(C_0, \ldots, C_{n-1})$ for univariate polynomials $C_i$,
one per variable of `ctxB`. Returns 1 for success and 0 for failure.
:::

::: function
int fmpz_mpoly_compose_fmpz_mpoly(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
int fmpz_mpoly_compose_fmpz_mpoly_geobucket(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)
int fmpz_mpoly_compose_fmpz_mpoly_horner(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz_mpoly_struct * const * C, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)

Set `A` to $B(C_0, \ldots, C_{n-1})$, where `B` belongs to `ctxB` and `A`
and the $n$ polynomials $C_i$ (one per variable of `ctxB`) belong to
`ctxAC`. `A` must not alias `B` or any $C_i$. Return 1 for success and 0
for failure. If every $C_i$ is zero or a monomial with coefficient 1,
`fmpz_mpoly_compose_fmpz_mpoly` maps the exponents by a matrix
(`_fmpz_mpoly_compose_mat`) and cannot fail; otherwise it uses the Horner
method if some $C_i$ has more than one term and the geobucket method
otherwise. The geobucket version computes powers of the $C_i$ and sums the
terms in a geobucket; the Horner version evaluates a recursive Horner
scheme.
:::

::: function
void fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_t A, const fmpz_mpoly_t B, const slong * c, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)

Sets `A` to `B` with the variable of index $i$ of `ctxB` replaced by the
variable of index `c[i]` of `ctxAC`, or by zero if `c[i]` is negative.
Each nonnegative `c[i]` must be less than the number of variables of
`ctxAC`. This renames, permutes, merges or drops variables at the cost of
mapping and sorting the monomials. `A` may alias `B` if the contexts are
the same.
:::

::: function
void _fmpz_mpoly_compose_mat(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mat_t M, const fmpz_mpoly_ctx_t ctxB, const fmpz_mpoly_ctx_t ctxAC)

Sets `A` to the image of `B` under the monomial map given by the matrix
`M`, as filled by `mpoly_compose_mat_gen` or
`mpoly_compose_mat_fill_column`: each unpacked monomial of `B` is
multiplied by `M`, terms whose image has a nonzero last entry are
dropped, and the result is sorted and combined. `A` must not alias `B`.
:::

::: function
int _fmpz_pow_ui_is_not_feasible(flint_bitcnt_t bbits, ulong e)
int _fmpz_pow_fmpz_is_not_feasible(flint_bitcnt_t bbits, const fmpz_t e)

Return whether raising an integer of `bbits` bits to the power $e \ge 0$
should be refused because the result would have about `WORD_MAX`$/16$ bits
or more.
:::

## Multiplication {#multiplication}

::: function
void fmpz_mpoly_mul(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

Sets `A` to $B C$. The algorithm (Johnson's heap method, an array method
or a dense Kronecker substitution, possibly with threads) is chosen from
the lengths, degrees and ordering of the inputs. See [](fmpz_mpoly/mul).
:::

::: function
void fmpz_mpoly_mul_johnson(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_mul_heap_threaded(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

Set `A` to $B C$ by Johnson's heap method [@Joh1974], which works for any
inputs and uses memory proportional to the input and output. The first
version uses one thread; the second divides the output among the
available threads [@MonPea2009].
:::

::: function
int fmpz_mpoly_mul_array(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_mul_array_threaded(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

Try to set `A` to $B C$ by accumulating the products in a dense array,
one slice of the main variable (lex) or of the total degree (degree
orderings) at a time. Return 0, leaving `A` unchanged, if the inputs do
not have one-word monomials or the array would be too large. The first
version uses one thread.
:::

::: function
int fmpz_mpoly_mul_dense(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

Tries to set `A` to $B C$ by Kronecker substitution into a univariate
product computed by `fmpz_poly_mul`. Returns 0 if the exponents of an
input do not fit a word or the dense size of the product overflows.
:::

::: function
void fmpz_mpoly_mul_monomial(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t C, const fmpz_mpoly_ctx_t ctx)

Sets `A` to $B C$ where `C` has exactly one term, in linear time.
:::

::: function
slong _fmpz_mpoly_mul_johnson(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N, const ulong * cmpmask)
void _fmpz_mpoly_mul_johnson_maxfields(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx)

The heap multiplication on raw term lists: multiply `(poly2, exp2, len2)`
by `(poly3, exp3, len3)` (both nonzero, monomials of $N$ words packed with
`bits` bits, which must suffice for the product), writing the product to
`(*poly1, *exp1)` with allocation `*alloc`, which are reallocated as
needed, and return its length. The heap has at most `len2` entries, so the
shorter input should come first. No aliasing is allowed. The `maxfields`
version takes precomputed field-wise maxima of the inputs (and clobbers
`maxBfields`), chooses the packing width and handles aliasing.
:::

::: function
void _fmpz_mpoly_mul_heap_threaded_pool_maxfields(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles)

Threaded heap multiplication with the given thread pool handles (and the
calling thread) and precomputed field maxima; `maxBfields` is clobbered.
:::

::: function
int _fmpz_mpoly_mul_array_LEX(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx)
int _fmpz_mpoly_mul_array_DEG(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx)
int _fmpz_mpoly_mul_array_threaded_pool_LEX(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles)
int _fmpz_mpoly_mul_array_threaded_pool_DEG(fmpz_mpoly_t A, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_t C, fmpz * maxCfields, const fmpz_mpoly_ctx_t ctx, const thread_pool_handle * handles, slong num_handles)

Array multiplication for lex and for degree orderings, single threaded or
with a thread pool, given the field-wise maxima of the nonzero inputs,
which must have one-word monomials. Return 0 without changing `A` if the
array would be too large: for lex the main variable must have degree less
than 300 in the product and the other variables must span at most
300000 array entries; for degree orderings $d^{n-1} \le 300000$ where $d$
is one more than the total degree of the product.
:::

::: function
int _fmpz_mpoly_mul_dense(fmpz_mpoly_t P, const fmpz_mpoly_t A, fmpz * maxAfields, const fmpz_mpoly_t B, fmpz * maxBfields, const fmpz_mpoly_ctx_t ctx)

Dense multiplication given the field-wise maxima of the nonzero inputs,
which must have single precision exponents. Returns 0 if the dense
representation does not fit.
:::

## Powering {#powering}

::: function
int fmpz_mpoly_pow_ui(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong k, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_pow_fmpz(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_t k, const fmpz_mpoly_ctx_t ctx)

Set `A` to $B^k$ and return 1, or return 0 if the result would be
unreasonably large. For $k \ge 3$ the power is computed by the
Monagan–Pearce FPS algorithm (`fmpz_mpoly_pow_fps`); squares are computed
by multiplication. `fmpz_mpoly_pow_ui` fails when $k (\mathit{len}(B) -
1)$ exceeds about `WORD_MAX`$/16$. `fmpz_mpoly_pow_fmpz` throws for
negative $k$, and for $k$ that does not fit an `slong` succeeds only if
$B$ is zero or a monomial with coefficient $\pm 1$.
:::

::: function
void fmpz_mpoly_pow_fps(fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong k, const fmpz_mpoly_ctx_t ctx)

Sets `A` to $B^k$ for nonzero $B$ and $k \ge 2$ by the FPS algorithm of
Monagan and Pearce [@MonPea2012], a heap-based sparse version of the
recurrence obtained from $B (B^k)' = k B' B^k$, which computes the terms
of $B^k$ and $B^{k-1}$ in order using $O(k \cdot \mathit{len}(B) \cdot
\mathit{len}(B^k))$ monomial operations in the worst case but far fewer
coefficient operations than repeated multiplication. See
[](algorithms/multivariate/multiplication#powering).
:::

## Division {#division}

Over $\Z$, division with remainder depends on the leading coefficient of
the divisor. The *divrem* functions perform the division algorithm with
floor division of coefficients, and the *quasidivrem* functions multiply
the dividend by a power of the leading coefficient first, so that their
remainder is fully reduced. See [](fmpz_mpoly/division).

::: function
int fmpz_mpoly_divides(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

If $B$ divides $A$ exactly, sets $Q$ to $A/B$ and returns 1; otherwise sets
$Q$ to zero and returns 0. Throws `FLINT_DIVZERO` if $B$ is zero. Uses
heap division, threaded for large dividends.
:::

::: function
void fmpz_mpoly_divexact(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Sets $Q$ to $A/B$, throwing `FLINT_ERROR` if the division is not exact.
:::

::: function
void fmpz_mpoly_divrem(fmpz_mpoly_t Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Sets $Q$ and $R$ such that $A = QB + R$, by the division algorithm in
which each quotient coefficient is the floor quotient of the current
coefficient by $\operatorname{lc}(B)$. A term of $R$ whose monomial is
divisible by $\operatorname{lm}(B)$ has a coefficient with the sign of
$\operatorname{lc}(B)$ and absolute value less than
$|\operatorname{lc}(B)|$. If $\operatorname{lc}(B) = \pm 1$ this is the
usual multivariate division with remainder.
:::

::: function
void fmpz_mpoly_div(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Sets $Q$ to the quotient of `fmpz_mpoly_divrem`, without computing the
remainder. It is much faster than `fmpz_mpoly_divides` when the division
is known to be exact.
:::

::: function
void fmpz_mpoly_quasidivrem(fmpz_t scale, fmpz_mpoly_t Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_quasidiv(fmpz_t scale, fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Set `scale`, $Q$ and $R$ such that $\mathit{scale} \cdot A = QB + R$ and no
monomial of $R$ is divisible by $\operatorname{lm}(B)$, where $\mathit{scale}
> 0$ is the least positive integer for which this division over $\Z$ is
possible (a divisor of a power of $\operatorname{lc}(B)$). The second
version does not compute $R$.
:::

::: function
void fmpz_mpoly_divrem_ideal(fmpz_mpoly_struct ** Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B, slong len, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_quasidivrem_ideal(fmpz_t scale, fmpz_mpoly_struct ** Q, fmpz_mpoly_t R, const fmpz_mpoly_t A, fmpz_mpoly_struct * const * B, slong len, const fmpz_mpoly_ctx_t ctx)

Versions of `fmpz_mpoly_divrem` and `fmpz_mpoly_quasidivrem` for a list of
`len` divisors $B_i$: set the quotients $Q_i$ and the remainder $R$ with
$(\mathit{scale} \cdot) A = \sum_i Q_i B_i + R$, by the multivariate
division algorithm that reduces each leading term by the first divisor
whose leading monomial divides it.
:::

::: function
int fmpz_mpoly_divides_monagan_pearce(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_divides_heap_threaded(fmpz_mpoly_t Q, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Exact division testing by the heap algorithm of Monagan and Pearce
[@MonPea2007], single threaded, respectively with threads
[@MonPea2010]. The threaded version, and the threaded path of
`fmpz_mpoly_divides`, only exist if FLINT was configured for a machine
with a strongly ordered memory model (`FLINT_KNOW_STRONG_ORDER`); test
with `#ifdef fmpz_mpoly_divides_heap_threaded`.
:::

::: function
int fmpz_mpoly_divides_array(fmpz_mpoly_t poly1, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

Exact division testing by a dense array method: returns 1 and sets
`poly1` to the quotient if it is exact, 0 if it is not, and $-1$ without
a result if the array method does not apply (the exponents need more than
one word, the ordering is `ORD_DEGREVLEX`, or the array would have more
than 300000 entries). It is efficient for dense inputs.
:::

::: function
void fmpz_mpoly_divrem_monagan_pearce(fmpz_mpoly_t q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_div_monagan_pearce(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

The heap algorithms behind `fmpz_mpoly_divrem` and `fmpz_mpoly_div`.
:::

::: function
int fmpz_mpoly_divrem_array(fmpz_mpoly_t q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

Division with remainder as in `fmpz_mpoly_divrem` by a dense array
method. Returns 0 without a result if the array would be too large or the
exponents do not fit one word, and 1 otherwise.
:::

::: function
void fmpz_mpoly_quasidivrem_heap(fmpz_t scale, fmpz_mpoly_t q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_quasidiv_heap(fmpz_t scale, fmpz_mpoly_t q, const fmpz_mpoly_t poly2, const fmpz_mpoly_t poly3, const fmpz_mpoly_ctx_t ctx)

The heap algorithms behind `fmpz_mpoly_quasidivrem` and
`fmpz_mpoly_quasidiv`. They throw `FLINT_DIVZERO` if the divisor is zero.
:::

::: function
void fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3, slong len, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_quasidivrem_ideal_heap(fmpz_t scale, fmpz_mpoly_struct ** q, fmpz_mpoly_t r, const fmpz_mpoly_t poly2, fmpz_mpoly_struct * const * poly3, slong len, const fmpz_mpoly_ctx_t ctx)

The heap algorithms behind `fmpz_mpoly_divrem_ideal` and
`fmpz_mpoly_quasidivrem_ideal`, with one heap for all divisors.
:::

::: function
slong _fmpz_mpoly_divides_monagan_pearce(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, flint_bitcnt_t bits, slong N, const ulong * cmpmask)

Exact division of nonzero raw term lists with monomials of $N$ words
packed with `bits` bits: writes the quotient to `(*poly1, *exp1)` with
allocation `*alloc`, reallocating as needed, and returns its length if
the division is exact and 0 otherwise. No aliasing is allowed.
:::

::: function
slong _fmpz_mpoly_div_monagan_pearce(fmpz ** polyq, ulong ** expq, slong * allocq, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong bits, slong N, const ulong * cmpmask)
slong _fmpz_mpoly_divrem_monagan_pearce(slong * lenr, fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong bits, slong N, const ulong * cmpmask)

Division as in `fmpz_mpoly_div` and `fmpz_mpoly_divrem` on raw nonzero term
lists. They return the length of the quotient, and the second sets
`*lenr` to the length of the remainder. No aliasing is allowed.
:::

::: function
slong _fmpz_mpoly_divrem_ideal_monagan_pearce(fmpz_mpoly_struct ** polyq, fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2, slong len2, fmpz_mpoly_struct * const * poly3, ulong * const * exp3, slong len, slong N, slong bits, const fmpz_mpoly_ctx_t ctx, const ulong * cmpmask)

Division by `len` divisors on raw term lists; `exp3[i]` holds the
monomials of `poly3[i]` repacked to `bits` bits. Returns the length of
the remainder.
:::

::: function
slong _fmpz_mpoly_divides_array(fmpz ** poly1, ulong ** exp1, slong * alloc, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, slong num, slong bits)
slong _fmpz_mpoly_divrem_array(slong * lenr, fmpz ** polyq, ulong ** expq, slong * allocq, fmpz ** polyr, ulong ** expr, slong * allocr, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3, slong * mults, slong num, slong bits)

Array division on raw term lists with one-word monomials in `num` fields
of `bits` bits, using the bases `mults` for the tight packing of the
array index. The first returns the length of the quotient if the division
is exact and 0 otherwise; the second returns the length of the quotient
and sets `*lenr`. No aliasing is allowed.
:::

::: function
int mpoly_divides_select_exps(fmpz_mpoly_t S, fmpz_mpoly_ctx_t zctx, slong nworkers, ulong * Aexp, slong Alen, ulong * Bexp, slong Blen, flint_bitcnt_t bits)

Chooses the monomial boundaries at which the dividend is split among
`nworkers` threads in the threaded exact division, storing them as the
monomials of `S`. Returns nonzero if it detects on the way that the
division cannot be exact. It is also used by `nmod_mpoly`.
:::

## Square roots {#sqrt}

::: function
int fmpz_mpoly_sqrt_heap(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx, int check)

If `poly2` is the square of a polynomial, sets `q` to the square root with
positive leading coefficient and returns 1; otherwise sets `q` to zero
and returns 0. The terms of the root are computed in order with a heap,
as in division [@MonPea2011]. With `check` equal to 0 the input is
assumed to be a square and the remaining terms are not verified, which
saves time but gives a meaningless result for non-squares.
:::

::: function
int fmpz_mpoly_sqrt(fmpz_mpoly_t q, const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_is_square(const fmpz_mpoly_t poly2, const fmpz_mpoly_ctx_t ctx)

`fmpz_mpoly_sqrt_heap` with checking, and a test whether `poly2` is a
square (which computes the root).
:::

::: function
slong _fmpz_mpoly_sqrt_heap(fmpz ** polyq, ulong ** expq, slong * allocq, const fmpz * poly2, const ulong * exp2, slong len2, flint_bitcnt_t bits, const mpoly_ctx_t mctx, int check)

The square root on raw term lists with monomials packed with `bits` bits;
returns the length of the root, or 0 if the input is not a square.
:::

## Content and GCD {#gcd}

::: function
void fmpz_mpoly_term_content(fmpz_mpoly_t M, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Sets $M$ to the GCD of the terms of $A$: the monomial whose exponents are
the minimal exponents of $A$, with coefficient the content (GCD of the
coefficients) of $A$, which is positive. If $A$ is zero, so is $M$.
:::

::: function
void fmpz_mpoly_primitive_part(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)

Sets `res` to $f$ divided by the content of its coefficients, with the sign
chosen so that the leading coefficient is positive. Zero stays zero.
:::

::: function
int fmpz_mpoly_content_vars(fmpz_mpoly_t g, const fmpz_mpoly_t A, slong * vars, slong vars_length, const fmpz_mpoly_ctx_t ctx)

Sets $g$ to the GCD of the coefficients of $A$ viewed as a polynomial in
the variables listed in `vars` (with coefficients in $\Z$ and the other
variables), and returns 1, or returns 0 if a GCD computation failed. On
success $g$ does not involve the variables `vars`. Throws if a variable
index is out of range. With no variables, $g = A$.
:::

::: function
int fmpz_mpoly_gcd(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Tries to set $G$ to the GCD of $A$ and $B$, normalised to have a positive
leading coefficient, and returns 1, or returns 0 if no algorithm
succeeded (which does not happen in practice). The GCD of zero and zero is
zero. See [](fmpz_mpoly/gcd).
:::

::: function
int fmpz_mpoly_gcd_cofactors(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

As `fmpz_mpoly_gcd`, and also sets $\bar A = A/G$ and $\bar B = B/G$.
If $A$ and $B$ are zero, all three outputs are zero; if only $A$ is zero,
$G = \pm B$, $\bar A = 0$ and $\bar B = \pm 1$. The three outputs must be
distinct.
:::

::: function
int fmpz_mpoly_gcd_brown(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_gcd_hensel(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_gcd_subresultant(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_gcd_zippel(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_gcd_zippel2(fmpz_mpoly_t G, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

Try to compute the GCD with one specific algorithm, after the same
preprocessing as `fmpz_mpoly_gcd` (which may solve easy cases directly):
Brown's dense modular algorithm, Hensel lifting (EEZ-GCD), the
subresultant pseudo-remainder sequence in a well chosen main variable,
Zippel's sparse modular algorithm, and the sparse algorithm based on
Berlekamp–Massey interpolation of Hu and Monagan. They return 0 if the
algorithm is not applicable or fails, for example `fmpz_mpoly_gcd_zippel2`
with fewer than three variables. Intended for testing and benchmarking.
:::

::: function
void fmpz_mpoly_deflation(fmpz * shift, fmpz * stride, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

For each variable $v$, with $S_v$ the set of exponents of $v$ in $A$, sets
`shift[v]` to $\min S_v$ and `stride[v]` to $\gcd\{e - \min S_v : e \in
S_v\}$ (zero if $S_v$ has one element). If $A$ is zero, all shifts and
strides are zero.
:::

::: function
void fmpz_mpoly_deflate(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_inflate(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz * shift, const fmpz * stride, const fmpz_mpoly_ctx_t ctx)

Apply $e \mapsto (e - \mathit{shift}[v])/\mathit{stride}[v]$ (deflate) or
$e \mapsto \mathit{shift}[v] + \mathit{stride}[v] \, e$ (inflate) to the
exponent of each variable $v$. Shifts and strides must be nonnegative.
When deflating, every $e - \mathit{shift}[v]$ must be a nonnegative
multiple of the stride, and a zero stride requires $e =
\mathit{shift}[v]$, which is mapped to zero; deflating with the output of
`fmpz_mpoly_deflation` therefore always works, and inflating undoes it.
:::

## Univariate polynomials with multivariate coefficients {#univar}

An `fmpz_mpoly_univar_t` holds a polynomial in one main variable $x_v$
whose coefficients are polynomials in the other variables, as a sparse
list of terms with `fmpz` exponents in decreasing order. It is the form
used by the resultant, discriminant and pseudo-remainder functions, and
the coefficients belong to the same context as the original polynomial
(they just do not involve $x_v$).

::: type
fmpz_mpoly_univar_struct
fmpz_mpoly_univar_t

A univariate polynomial with `length` terms, coefficients `coeffs` (of type
`fmpz_mpoly_struct`) and exponents `exps` (of type `fmpz`), with room for
`alloc` terms.
:::

::: function
void fmpz_mpoly_univar_init(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_univar_clear(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_univar_fit_length(fmpz_mpoly_univar_t A, slong length, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_univar_swap(fmpz_mpoly_univar_t A, fmpz_mpoly_univar_t B, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_univar_zero(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)

Initialise as zero, free, ensure room for `length` terms, swap in constant
time, and set to zero.
:::

::: function
void fmpz_mpoly_to_univar(fmpz_mpoly_univar_t A, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)

Sets `A` to `B` written as a polynomial in the variable of index `var`.
:::

::: function
void fmpz_mpoly_from_univar(fmpz_mpoly_t A, const fmpz_mpoly_univar_t B, slong var, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_from_univar(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_univar_t B, slong var, const fmpz_mpoly_ctx_t ctx)

Set `A` to the polynomial $\sum_i B_i x_{\mathit{var}}^{e_i}$. The
coefficients of `B` must not involve the variable `var`. The underscore
version packs the result into `Abits` bits, which must suffice.
:::

::: function
void fmpz_mpoly_univar_set_coeff_ui(fmpz_mpoly_univar_t A, ulong e, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

Sets the coefficient of $x^e$ in `A` to $c$, inserting or removing a term
as needed.
:::

::: function
int fmpz_mpoly_univar_degree_fits_si(const fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)
slong fmpz_mpoly_univar_length(const fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)
slong fmpz_mpoly_univar_get_term_exp_si(fmpz_mpoly_univar_t A, slong i, const fmpz_mpoly_ctx_t ctx)

Return whether the degree fits an `slong`, the number of terms, and the
exponent of term $i$ (which must fit).
:::

::: function
void fmpz_mpoly_univar_get_term_coeff(fmpz_mpoly_t c, const fmpz_mpoly_univar_t A, slong i, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_univar_swap_term_coeff(fmpz_mpoly_t c, fmpz_mpoly_univar_t A, slong i, const fmpz_mpoly_ctx_t ctx)

Set $c$ to, respectively swap $c$ with, the coefficient of term $i$.
:::

::: function
void fmpz_mpoly_univar_print_pretty(const fmpz_mpoly_univar_t A, const char ** x, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_univar_assert_canonical(fmpz_mpoly_univar_t A, const fmpz_mpoly_ctx_t ctx)

Print `A` with the main variable written `X` and the coefficients in the
variables `x`, and check that the exponents are strictly decreasing and
the coefficients nonzero and canonical (throwing otherwise).
:::

::: function
int fmpz_mpoly_univar_pseudo_gcd(fmpz_mpoly_univar_t gx, const fmpz_mpoly_univar_t ax, const fmpz_mpoly_univar_t bx, const fmpz_mpoly_ctx_t ctx)

Sets `gx` to the last nonzero subresultant of `ax` and `bx` (Ducos'
algorithm), a pseudo-GCD: a multiple of their GCD by an element of the
coefficient ring. Requires $\deg \mathit{ax} \ge \deg \mathit{bx}$ and
nonzero inputs. Returns 1 for success and 0 for failure.
:::

::: function
int fmpz_mpoly_univar_resultant(fmpz_mpoly_t d, const fmpz_mpoly_univar_t ax, const fmpz_mpoly_univar_t bx, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_univar_discriminant(fmpz_mpoly_t d, const fmpz_mpoly_univar_t fx, const fmpz_mpoly_ctx_t ctx)

Set $d$ to the resultant of `ax` and `bx`, respectively the discriminant
of `fx`, in the main variable, by the subresultant algorithm. Return 1 for
success and 0 for failure.
:::

## Resultants and discriminants {#resultants}

::: function
int fmpz_mpoly_resultant(fmpz_mpoly_t R, const fmpz_mpoly_t A, const fmpz_mpoly_t B, slong var, const fmpz_mpoly_ctx_t ctx)

Sets $R$ to the resultant of $A$ and $B$ with respect to the variable of
index `var`, and returns 1, or returns 0 on failure. If the context has
two variables and the inputs are dense enough and of large enough degree,
the resultant is computed by a dense multimodular algorithm
(`fmpz_bpoly_resultant`, which uses the evaluation and interpolation
algorithm of [](algorithms/multivariate/modp-bivariate-resultant) modulo
each prime); otherwise by Ducos' subresultant algorithm on
the univariate forms [@Duc2000].
:::

::: function
int fmpz_mpoly_discriminant(fmpz_mpoly_t R, const fmpz_mpoly_t A, slong var, const fmpz_mpoly_ctx_t ctx)

Sets $R$ to the discriminant of $A$ with respect to the variable of index
`var`, $(-1)^{d(d-1)/2} \operatorname{res}(A, A')/\operatorname{lc}(A)$
up to the convention for low degrees used by
`mpoly_univar_discriminant`, and returns 1, or returns 0 on failure.
:::

## Special polynomials {#special}

::: function
void fmpz_mpoly_symmetric_gens(fmpz_mpoly_t res, ulong k, slong * vars, slong n, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_symmetric(fmpz_mpoly_t res, ulong k, const fmpz_mpoly_ctx_t ctx)

Set `res` to the elementary symmetric polynomial $e_k$ in all variables,
respectively in the $n$ distinct variables of indices `vars`: the sum of
all products of $k$ distinct variables. $e_0 = 1$, and $e_k = 0$ for $k >
n$.
:::

## Vectors of polynomials {#vectors}

::: type
fmpz_mpoly_vec_struct
fmpz_mpoly_vec_t

A vector of polynomials in one context: an array `p` of `alloc`
initialised polynomials, of which the first `length` are in use.
:::

::: macro
fmpz_mpoly_vec_entry(vec, i)

A pointer to entry $i$ of `vec`.
:::

::: function
void fmpz_mpoly_vec_init(fmpz_mpoly_vec_t vec, slong len, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_vec_clear(fmpz_mpoly_vec_t vec, const fmpz_mpoly_ctx_t ctx)

Initialise `vec` as a vector of `len` zero polynomials, and free it.
:::

::: function
void fmpz_mpoly_vec_print(const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)

Prints the vector as `[f0, f1, ...]` with the default variable names.
:::

::: function
void fmpz_mpoly_vec_swap(fmpz_mpoly_vec_t x, fmpz_mpoly_vec_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_vec_fit_length(fmpz_mpoly_vec_t vec, slong len, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_vec_set(fmpz_mpoly_vec_t dest, const fmpz_mpoly_vec_t src, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_vec_set_length(fmpz_mpoly_vec_t vec, slong len, const fmpz_mpoly_ctx_t ctx)

Swap two vectors, ensure room for `len` entries, copy, and change the
length (truncating, or appending zero polynomials).
:::

::: function
void fmpz_mpoly_vec_append(fmpz_mpoly_vec_t vec, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)
slong fmpz_mpoly_vec_insert_unique(fmpz_mpoly_vec_t vec, const fmpz_mpoly_t f, const fmpz_mpoly_ctx_t ctx)

Append a copy of $f$. `fmpz_mpoly_vec_insert_unique` appends it only if
it is not already present (by linear search) and returns its index.
:::

::: function
void fmpz_mpoly_vec_randtest_not_zero(fmpz_mpoly_vec_t vec, flint_rand_t state, slong len, slong poly_len, slong bits, ulong exp_bound, fmpz_mpoly_ctx_t ctx)

Sets `vec` to `len` random nonzero polynomials generated by
`fmpz_mpoly_randtest_bound` with the given parameters.
:::

::: function
void fmpz_mpoly_vec_set_primitive_unique(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)

Sets `G` to the primitive parts of the nonzero entries of `F`, without
duplicates, in an unspecified order.
:::

## Ideals and Gröbner bases {#groebner}

These functions treat ideals of $\Q[x_0, \ldots, x_{n-1}]$, using
primitive integer polynomials with positive leading coefficient as
normalised generators instead of monic rational polynomials. The
implementation of Buchberger's algorithm is simple and meant for small
problems. See [@BW1993] and [@CLO2015].

::: function
void fmpz_mpoly_spoly(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_t g, const fmpz_mpoly_ctx_t ctx)

Sets `res` to the S-polynomial of $f$ and $g$, scaled to integer
coefficients: with $L = \operatorname{lcm}(\operatorname{lm} f,
\operatorname{lm} g)$ and $c = \operatorname{lcm}(\operatorname{lc} f,
\operatorname{lc} g)$, $\frac{cL}{\operatorname{lt} f} f -
\frac{cL}{\operatorname{lt} g} g$.
:::

::: function
void fmpz_mpoly_reduction_primitive_part(fmpz_mpoly_t res, const fmpz_mpoly_t f, const fmpz_mpoly_vec_t Iv, const fmpz_mpoly_ctx_t ctx)

Sets `res` to the primitive part of the remainder of $f$ on division by
the polynomials `Iv` (`fmpz_mpoly_quasidivrem_ideal`).
:::

::: function
int fmpz_mpoly_vec_is_groebner(const fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)

Returns whether `G` is a Gröbner basis (all S-polynomials reduce to zero)
and, if `F` is not `NULL`, whether every element of `F` reduces to zero
modulo `G`, so that `G` is a Gröbner basis of the ideal generated by `F`
if `G` is contained in it.
:::

::: function
int fmpz_mpoly_vec_is_autoreduced(const fmpz_mpoly_vec_t G, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_vec_autoreduction(fmpz_mpoly_vec_t H, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_vec_autoreduction_groebner(fmpz_mpoly_vec_t H, const fmpz_mpoly_vec_t G, const fmpz_mpoly_ctx_t ctx)

Test whether the vector is autoreduced (inter-reduced): no term of an
element is divisible by the leading monomial of another element; and set
`H` to the autoreduction of `F`. The `groebner` version assumes that `G`
is a Gröbner basis and produces the reduced Gröbner basis, which is unique
up to the order of the entries.
:::

::: function
void fmpz_mpoly_buchberger_naive(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_buchberger_naive_with_limits(fmpz_mpoly_vec_t G, const fmpz_mpoly_vec_t F, slong ideal_len_limit, slong poly_len_limit, slong poly_bits_limit, const fmpz_mpoly_ctx_t ctx)

Set `G` to a Gröbner basis of the ideal generated by `F`, by Buchberger's
algorithm [@Buc1965] with the usual criteria to skip pairs, following
[@BW1993]. The second version gives up and returns 0 if the basis gets
longer than `ideal_len_limit`, a polynomial longer than `poly_len_limit`,
or a coefficient larger than `poly_bits_limit` bits; `G` then generates
the ideal but need not be a Gröbner basis. It returns 1 on success.
:::

## Internal functions {#internal}

These functions are used by other modules and may change without notice.

::: function
void mpoly_void_ring_init_fmpz_mpoly_ctx(mpoly_void_ring_t R, const fmpz_mpoly_ctx_t ctx)

Fills `R` with the function pointers of the ring $\Z[x_0, \ldots,
x_{n-1}]$, for the generic code in `mpoly`.
:::

::: function
void fmpz_mpolyl_lead_coeff(fmpz_mpoly_t c, const fmpz_mpoly_t A, slong num_vars, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpolyl_content(fmpz_mpoly_t g, const fmpz_mpoly_t A, slong num_vars, const fmpz_mpoly_ctx_t ctx)

For a lex context and $0 < \mathit{num\_vars} < n$, view $A \ne 0$ as a
polynomial in $x_0, \ldots, x_{\mathit{num\_vars}-1}$ with coefficients in
the other variables, and set $c$ to its leading coefficient, respectively
$g$ to the GCD of its coefficients (returning 0 if a GCD failed). The
outputs must not alias $A$ and have the packing width of $A$.
:::

::: function
void fmpz_mpoly_to_fmpz_poly(fmpz_poly_t poly1, slong * poly1_shift, const fmpz_mpoly_t poly2, slong var, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_from_fmpz_poly(fmpz_mpoly_t poly1, const fmpz_poly_t poly2, slong shift2, slong var, const fmpz_mpoly_ctx_t ctx)

Convert between a polynomial in the single variable `var` and a
univariate polynomial with a shift: the first sets `poly1` and
`*poly1_shift` so that `poly2` equals $\mathit{poly1}(x) x^{\mathit{shift}}$
(it throws for multiprecision exponents), and the second sets `poly1` to
$\mathit{poly2}(x_{\mathit{var}}) x_{\mathit{var}}^{\mathit{shift2}}$.
:::

::: function
void _fmpz_mpoly_to_fmpz_poly_deflate(fmpz_poly_t A, const fmpz_mpoly_t B, slong var, const ulong * Bshift, const ulong * Bstride, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_from_fmpz_poly_inflate(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_poly_t B, slong var, const ulong * Ashift, const ulong * Astride, const fmpz_mpoly_ctx_t ctx)

Convert between a polynomial whose exponents are $\mathit{shift} + k \cdot
\mathit{stride}$ in each variable, with only `var` varying, and a
univariate polynomial in $k$. Used by the GCD code.
:::

::: function
int fmpz_mpoly_repack_bits(fmpz_mpoly_t A, const fmpz_mpoly_t B, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_repack_bits_inplace(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx)

Set `A` to `B` (or change `A` in place) with exponents packed into
`Abits` bits, and return 1, or return 0 if the exponents do not fit,
leaving `A` unchanged. The first rounds `Abits` with `mpoly_fix_bits`;
the second requires a valid width.
:::

::: function
void fmpz_mpoly_height(fmpz_t max, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_heights(fmpz_t max, fmpz_t sum, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Set `max` to the largest absolute value of a coefficient of `A`, and
`sum` to the sum of the absolute values.
:::

::: type
fmpz_pow_cache_t

A cache of the powers $1, v, v^2, \ldots$ of an integer $v$, used when
evaluating polynomials at $v$.
:::

::: function
void fmpz_pow_cache_init(fmpz_pow_cache_t T, const fmpz_t val)
void fmpz_pow_cache_clear(fmpz_pow_cache_t T)

Initialise a cache for the powers of `val`, and free it.
:::

::: function
int fmpz_pow_cache_mulpow_ui(fmpz_t a, const fmpz_t b, ulong k, fmpz_pow_cache_t T)
int fmpz_pow_cache_mulpow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t k, fmpz_pow_cache_t T)

Set $a$ to $b v^k$, extending the cache for $k \le 100$ and computing the
power directly for larger $k$, and return 1; the `fmpz` version returns 0
if the power cannot be computed.
:::

::: type
fmpz_mpoly_geobucket_struct
fmpz_mpoly_geobucket_t

A geobucket [@Yan1998]: an array of polynomials of geometrically
increasing lengths (powers of 4) whose sum is the value. Adding a short
polynomial to a long sum costs time proportional to the short one,
amortised.
:::

::: function
void fmpz_mpoly_geobucket_init(fmpz_mpoly_geobucket_t B, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_geobucket_clear(fmpz_mpoly_geobucket_t B, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_geobucket_fit_length(fmpz_mpoly_geobucket_t B, slong i, const fmpz_mpoly_ctx_t ctx)

Initialise a geobucket with value zero, free it, and make sure it has at
least `i` (zeroed) levels.
:::

::: function
void fmpz_mpoly_geobucket_set(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_geobucket_add(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_geobucket_sub(fmpz_mpoly_geobucket_t B, fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_geobucket_empty(fmpz_mpoly_t p, fmpz_mpoly_geobucket_t B, const fmpz_mpoly_ctx_t ctx)

Set the value of `B` to $p$, add $p$ to it, subtract $p$ from it, and move
the value to $p$, leaving `B` zero. The first three functions take $p$
over (its value is destroyed).
:::

::: type
fmpz_mpolyd_struct
fmpz_mpolyd_t
fmpz_mpolyd_ctx_struct
fmpz_mpolyd_ctx_t

A dense polynomial in `nvars` variables with degree bounds
`deg_bounds[i]` (exclusive): the coefficient of $x_0^{e_0} \cdots
x_{n-1}^{e_{n-1}}$ is stored at index $e_{n-1} + r_{n-1}(e_{n-2} +
r_{n-2}(\cdots))$. Used by the dense multiplication.
:::

::: function
void fmpz_mpolyd_init(fmpz_mpolyd_t poly, slong nvars)
void fmpz_mpolyd_fit_length(fmpz_mpolyd_t poly, slong len)
void fmpz_mpolyd_clear(fmpz_mpolyd_t poly)

Initialise, ensure room for `len` coefficients, and free.
:::

::: type
fmpz_mpoly_stripe_struct
fmpz_mpoly_stripe_t

Working memory and parameters (exponent ranges, packing, coefficient
size) of one thread of the threaded multiplication and division.
:::

::: function
void _fmpz_mpoly_mul_array_chunked_LEX(fmpz_mpoly_t P, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const ulong * mults, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_mul_array_chunked_DEG(fmpz_mpoly_t P, const fmpz_mpoly_t A, const fmpz_mpoly_t B, ulong degb, const fmpz_mpoly_ctx_t ctx)

The array multiplication proper: split both inputs into chunks by the
main variable (lex, with array bases `mults`) or by total degree (degree
orderings, with base `degb`), and compute the chunks of the product one
after another, each in a dense array, with coefficients accumulated in
one, two or three words when a bound on their size allows it.
:::

::: function
void _fmpz_mpoly_addmul_array1_slong1(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3)
void _fmpz_mpoly_addmul_array1_slong2(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3)
void _fmpz_mpoly_addmul_array1_slong(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3)
void _fmpz_mpoly_addmul_array1_fmpz(fmpz * poly1, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3)

Add the products of all pairs of terms of the two term lists, whose
monomials are tightly packed array indices, to the dense array `poly1`,
with coefficients of one, two or three words (two's complement), or of
type `fmpz`. The inputs are processed in blocks for cache locality.
:::

::: function
void _fmpz_mpoly_submul_array1_slong1(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3)
void _fmpz_mpoly_submul_array1_slong2(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3)
void _fmpz_mpoly_submul_array1_slong(ulong * poly1, const slong * poly2, const ulong * exp2, slong len2, const slong * poly3, const ulong * exp3, slong len3)
void _fmpz_mpoly_submul_array1_fmpz(fmpz * poly1, const fmpz * poly2, const ulong * exp2, slong len2, const fmpz * poly3, const ulong * exp3, slong len3)
void _fmpz_mpoly_submul_array1_slong2_1(ulong * poly1, slong d, const ulong exp2, const slong * poly3, const ulong * exp3, slong len3)
void _fmpz_mpoly_submul_array1_slong_1(ulong * poly1, slong d, const ulong exp2, const slong * poly3, const ulong * exp3, slong len3)
void _fmpz_mpoly_submul_array1_fmpz_1(fmpz * poly1, const fmpz_t d, ulong exp2, const fmpz * poly3, const ulong * exp3, slong len3)

As the `addmul` functions, subtracting; the `_1` versions subtract the
product of a single term $d \cdot \mathit{exp2}$ with a term list. Used by
the array division.
:::

::: function
slong fmpz_mpoly_append_array_sm1_LEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, const ulong * mults, slong num, slong array_size, slong top)
slong fmpz_mpoly_append_array_sm2_LEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, const ulong * mults, slong num, slong array_size, slong top)
slong fmpz_mpoly_append_array_sm3_LEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, const ulong * mults, slong num, slong array_size, slong top)
slong fmpz_mpoly_append_array_fmpz_LEX(fmpz_mpoly_t P, slong Plen, fmpz * coeff_array, const ulong * mults, slong num, slong array_size, slong top)
slong fmpz_mpoly_append_array_sm1_DEGLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb)
slong fmpz_mpoly_append_array_sm2_DEGLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb)
slong fmpz_mpoly_append_array_sm3_DEGLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb)
slong fmpz_mpoly_append_array_fmpz_DEGLEX(fmpz_mpoly_t P, slong Plen, fmpz * coeff_array, slong top, slong nvars, slong degb)
slong fmpz_mpoly_append_array_sm1_DEGREVLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb)
slong fmpz_mpoly_append_array_sm2_DEGREVLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb)
slong fmpz_mpoly_append_array_sm3_DEGREVLEX(fmpz_mpoly_t P, slong Plen, ulong * coeff_array, slong top, slong nvars, slong degb)
slong fmpz_mpoly_append_array_fmpz_DEGREVLEX(fmpz_mpoly_t P, slong Plen, fmpz * coeff_array, slong top, slong nvars, slong degb)

Append the nonzero entries of a dense coefficient array (one chunk of a
product: main variable or total degree `top`) to `P` in decreasing
monomial order, starting at term index `Plen`, clear the array, and return
the new length. `sm1`, `sm2` and `sm3` arrays hold one, two or three words
per coefficient.
:::

::: function
slong _fmpz_mpoly_from_ulong_array1(fmpz ** poly1, ulong ** exp1, slong * alloc, ulong * poly2, const slong * mults, slong num, slong bits, slong k)
slong _fmpz_mpoly_from_ulong_array2(fmpz ** poly1, ulong ** exp1, slong * alloc, ulong * poly2, const slong * mults, slong num, slong bits, slong k)
slong _fmpz_mpoly_from_ulong_array(fmpz ** poly1, ulong ** exp1, slong * alloc, ulong * poly2, const slong * mults, slong num, slong bits, slong k)
slong _fmpz_mpoly_from_fmpz_array(fmpz ** poly1, ulong ** exp1, slong * alloc, fmpz * poly2, const slong * mults, slong num, slong bits, slong k)

Convert a dense coefficient array indexed by tightly packed monomials with
the bases `mults` (one, two or three words per coefficient, or `fmpz`)
into a term list with exponents packed into `bits` bits, starting at term
index `k`, and return the new length. Used by the array division.
:::

::: function
void _fmpz_mpoly_to_ulong_array1(ulong * p, const fmpz * coeffs, const ulong * exps, slong len)
void _fmpz_mpoly_to_ulong_array2(ulong * p, const fmpz * coeffs, const ulong * exps, slong len)
void _fmpz_mpoly_to_ulong_array(ulong * p, const fmpz * coeffs, const ulong * exps, slong len)
void _fmpz_mpoly_to_fmpz_array(fmpz * p, const fmpz * coeffs, const ulong * exps, slong len)

Scatter a term list whose monomials are tightly packed array indices into
a dense array with one, two or three words per coefficient, or with `fmpz`
coefficients.
:::

::: function
void _fmpz_mpoly_add_uiuiui_fmpz(ulong * c, const fmpz_t d)
void _fmpz_mpoly_sub_uiuiui_fmpz(ulong * c, const fmpz_t d)
void _fmpz_mpoly_addmul_uiuiui_fmpz(ulong * c, slong d1, slong d2)
void _fmpz_mpoly_submul_uiuiui_fmpz(ulong * c, slong d1, slong d2)

Add or subtract $d$, or $d_1 d_2$, to or from the three-word two's
complement integer `(c[0], c[1], c[2])`. In the first two, $d$ must fit
in three words.
:::

::: function
void fmpz_mpoly_remainder_test(const fmpz_mpoly_t r, const fmpz_mpoly_t g, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_remainder_strongtest(const fmpz_mpoly_t r, const fmpz_mpoly_t g, const fmpz_mpoly_ctx_t ctx)

Consistency checks for the test suite: throw `FLINT_ERROR` if $r$ is not
a valid remainder for division by the nonzero $g$, that is, if a term of
$r$ divisible by $\operatorname{lm}(g)$ has a coefficient of absolute
value at least $|\operatorname{lc}(g)|$ (first version), or if any term of
$r$ is divisible by $\operatorname{lm}(g)$ at all (second version, for
quasi-division).
:::
