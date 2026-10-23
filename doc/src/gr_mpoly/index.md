---
title: gr_mpoly.h
header: gr_mpoly.h
description: Sparse multivariate polynomials over any generic ring, with heap-based multiplication and division (also multithreaded), division with remainder by ideals and pseudo-division, as a generic ring in its own right.
category: Generic rings
order: 9
---

This module implements sparse multivariate polynomials $f \in R[X_1,
\ldots, X_n]$ over a ring $R$ given by a generic context. A polynomial
is stored as an array of coefficients in $R$ and an array of packed
exponent vectors, sorted by a monomial ordering, in the same format as
`fmpz_mpoly` and the other multivariate modules ([](mpoly)). Unlike
`gr_poly`, the context of this module is the context of the polynomial
ring itself: a `gr_mpoly_ctx_t` is a `gr_ctx_t` whose elements are
`gr_mpoly_t`, so the polynomials can be used with the generic functions
(`gr_add`, `gr_mul`, `gr_div`, `gr_set_str`, ...), in matrices
(`gr_mat`), as coefficients of univariate polynomials and in any other
generic algorithm.

Products are computed with Johnson's heap algorithm and exact quotients
and divisions with remainder with the heap algorithm of Monagan and
Pearce [@MonPea2007; @MonPea2009], in parallel for large inputs over
thread-safe coefficient rings ([](algorithms/multivariate/multiplication),
[](algorithms/multivariate/division)); see
[](algorithms/multivariate/overview) for the specialised multivariate
modules
(`fmpz_mpoly`, `nmod_mpoly`, `fq_nmod_mpoly`, ...), which are much faster
over their rings and also provide GCDs and factorisation.

## Types and conventions {#conventions}

::: type
gr_mpoly_struct
gr_mpoly_t

A multivariate polynomial: the fields are `coeffs` (an array of
coefficients in $R$), `exps` (the packed exponent vectors, `bits` bits
per exponent field), `length` (the number of terms), `bits`, and the
allocated sizes `coeffs_alloc` and `exps_alloc`. The terms are sorted
in decreasing order of their monomials with respect to the ordering of
the context, and no two terms have the same monomial. A `gr_mpoly_t` is
an array of one `gr_mpoly_struct`.
:::

A polynomial is *normalised* by removing terms whose coefficient is
provably zero. Over rings with undecidable zero tests (balls, exact
real numbers) a coefficient that may be zero is kept, as in $3xy^2 +
[\pm 0.01] xy$, so that the length and the degrees are upper bounds;
functions that need the exact term structure, such as division, return
`GR_UNABLE` for such inputs. Functions returning `int` return status
codes as described in [](gr), and functions without a leading
underscore allow aliasing of inputs and outputs.

::: type
gr_mpoly_ctx_struct
gr_mpoly_ctx_t

The context of a polynomial ring $R[X_1, \ldots, X_n]$. It is a
`gr_ctx_t`, so the generic `gr` functions can be used with it: `gr_add`
on `gr_mpoly_t` elements of this context is `gr_mpoly_add`. The
context data (of type `_gr_mpoly_ctx_struct`) consists of the
coefficient context `cctx`, which is not owned by the polynomial
context and must stay alive as long as it is used; an `mpoly_ctx_t`
`mctx` with the number of variables and the monomial ordering, which is
owned; and an optional array `vars` of names of the generators, set
with `gr_mpoly_ctx_set_gen_names` (by default `NULL`, meaning default
names). The names are used for printing, for parsing with `gr_set_str`
and for conversions between rings.
:::

::: macro
GR_MPOLY_CCTX(ctx)
GR_MPOLY_MCTX(ctx)
GR_MPOLY_VARS(ctx)
GR_MPOLY_NVARS(ctx)

The coefficient context, the `mpoly_ctx_struct` of the exponent format
and ordering, the array of generator names, and the number of variables
of a `gr_mpoly` context.
:::

## Contexts {#contexts}

::: function
void gr_mpoly_ctx_init(gr_mpoly_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord)

Initialises `ctx` as the ring of polynomials in `nvars` variables over
`base_ring`, with the monomial ordering `ord` (`ORD_LEX`, `ORD_DEGLEX`
or `ORD_DEGREVLEX`).
:::

::: function
void gr_mpoly_ctx_init_rand(gr_mpoly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring, slong max_nvars)

Initialises `ctx` with a random number of variables between 0 and
`max_nvars` and a random ordering, for testing.
:::

::: function
void gr_mpoly_ctx_clear(gr_mpoly_ctx_t ctx)

Clears the context (but not the coefficient ring).
:::

::: function
int gr_mpoly_ctx_set_gen_names(gr_mpoly_ctx_t ctx, const char ** s)

Sets the names of the $n$ generators to copies of the strings `s[0]`,
…, `s[n-1]`.
:::

::: function
int gr_mpoly_ctx_write(gr_stream_t out, gr_mpoly_ctx_t ctx)
truth_t gr_mpoly_ctx_is_ring(gr_mpoly_ctx_t ctx)
truth_t gr_mpoly_ctx_is_zero_ring(gr_mpoly_ctx_t ctx)
truth_t gr_mpoly_ctx_is_commutative_ring(gr_mpoly_ctx_t ctx)
truth_t gr_mpoly_ctx_is_integral_domain(gr_mpoly_ctx_t ctx)
truth_t gr_mpoly_ctx_is_field(gr_mpoly_ctx_t ctx)
truth_t gr_mpoly_ctx_is_threadsafe(gr_mpoly_ctx_t ctx)

Context methods: a description of the ring, and structural predicates
inherited from the coefficient ring. The polynomial ring is a field
only if it has no variables and the coefficient ring is a field.
:::

## Memory management {#memory}

::: function
void gr_mpoly_init(gr_mpoly_t A, gr_mpoly_ctx_t ctx)
void gr_mpoly_init2(gr_mpoly_t A, slong alloc, gr_mpoly_ctx_t ctx)
void gr_mpoly_init3(gr_mpoly_t A, slong alloc, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx)

Initialise `A` as the zero polynomial; `init2` and `init3` allocate
space for `alloc` terms, with the default number of bits per exponent
or with `bits` bits.
:::

::: function
void gr_mpoly_clear(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

Clears `A`.
:::

::: function
void gr_mpoly_fit_length(gr_mpoly_t A, slong len, gr_mpoly_ctx_t ctx)
void gr_mpoly_fit_bits(gr_mpoly_t A, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx)
void gr_mpoly_fit_length_fit_bits(gr_mpoly_t A, slong len, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx)
void gr_mpoly_fit_length_reset_bits(gr_mpoly_t A, slong len, flint_bitcnt_t bits, gr_mpoly_ctx_t ctx)

Ensure space for `len` terms, or exponent fields of at least `bits` bits
(repacking the existing exponents), or both. The *reset_bits* version
sets the number of bits to `bits` without preserving the exponents,
for use before writing a new polynomial.
:::

::: function
void _gr_mpoly_fit_length(gr_ptr * coeffs, slong * coeffs_alloc, ulong ** exps, slong * exps_alloc, slong N, slong length, gr_ctx_t cctx)

Grows the arrays `*coeffs` and `*exps` (with `N` words per exponent
vector) to hold at least `length` terms, updating the allocation
counters. Here `cctx` is the coefficient ring.
:::

::: function
void _gr_mpoly_set_length(gr_mpoly_t A, slong newlen, gr_mpoly_ctx_t ctx)
slong gr_mpoly_length(const gr_mpoly_t x, gr_mpoly_ctx_t ctx)

Set the number of terms of `A` (which must be allocated), without
normalising, and return the number of terms.
:::

::: function
void _gr_mpoly_normalise(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

Removes the terms whose coefficient is provably zero. All functions
return normalised polynomials; this is for code that writes terms
directly.
:::

::: function
void gr_mpoly_swap(gr_mpoly_t A, gr_mpoly_t B, gr_mpoly_ctx_t ctx)
void gr_mpoly_set_shallow(gr_mpoly_t res, const gr_mpoly_t poly, gr_mpoly_ctx_t ctx)

Swap two polynomials in constant time, and make `res` a shallow copy of
`poly` sharing its arrays (which must then not be cleared twice).
:::

## Vectors of polynomials {#vectors}

::: type
gr_mpoly_vec_struct
gr_mpoly_vec_t

A resizable vector of polynomials over one context, with the fields
`entries`, `alloc` and `length`. It has the layout of a `gr_vec_t`
whose elements are polynomials, and can be cast to and from it.
:::

::: function
void gr_mpoly_vec_init(gr_mpoly_vec_t vec, slong len, gr_mpoly_ctx_t ctx)
void gr_mpoly_vec_clear(gr_mpoly_vec_t vec, gr_mpoly_ctx_t ctx)

Initialise `vec` as `len` zero polynomials, and clear it.
:::

::: function
gr_mpoly_struct * gr_mpoly_vec_entry_ptr(gr_mpoly_vec_t vec, slong i, gr_mpoly_ctx_t ctx)
const gr_mpoly_struct * gr_mpoly_vec_entry_srcptr(const gr_mpoly_vec_t vec, slong i, gr_mpoly_ctx_t ctx)
slong gr_mpoly_vec_length(const gr_mpoly_vec_t vec, gr_mpoly_ctx_t ctx)

Return a pointer to the entry with index $0 \le i <$ `length`, and the
length.
:::

::: function
void gr_mpoly_vec_fit_length(gr_mpoly_vec_t vec, slong len, gr_mpoly_ctx_t ctx)
void gr_mpoly_vec_set_length(gr_mpoly_vec_t vec, slong len, gr_mpoly_ctx_t ctx)
int gr_mpoly_vec_set(gr_mpoly_vec_t res, const gr_mpoly_vec_t src, gr_mpoly_ctx_t ctx)
int gr_mpoly_vec_append(gr_mpoly_vec_t vec, const gr_mpoly_t f, gr_mpoly_ctx_t ctx)

Allocate space for `len` entries, change the length (appending zero
polynomials), copy a vector, and append a copy of `f`.
:::

## Assignment, generators and comparisons {#basic}

::: function
int gr_mpoly_set(gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_zero(gr_mpoly_t A, gr_mpoly_ctx_t ctx)
int gr_mpoly_one(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

Set `A` to a copy of `B`, to zero or to one.
:::

::: function
int gr_mpoly_set_scalar(gr_mpoly_t A, gr_srcptr c, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_ui(gr_mpoly_t A, ulong c, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_si(gr_mpoly_t A, slong c, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_fmpz(gr_mpoly_t A, const fmpz_t c, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_fmpq(gr_mpoly_t A, const fmpq_t c, gr_mpoly_ctx_t ctx)

Set `A` to the constant polynomial $c$.
:::

::: function
int gr_mpoly_set_other(gr_mpoly_t res, gr_srcptr A, gr_ctx_t A_ctx, gr_mpoly_ctx_t ctx)

Converts the element `A` of the ring `A_ctx` into the polynomial ring.

::: details
If `A_ctx` is a multivariate polynomial ring, the coefficients are
converted and the generators matched: by name if both rings have named
generators (every generator used in `A` must have a namesake in
`ctx`), and directly if both rings have the same number of unnamed
generators and the same ordering. Other cases are not supported.
Otherwise `A` is converted as a constant into the coefficient ring.
Absorbing the generators of nested rings, as from $R[x, y][s, t]$ to
$R[x, y, s, t]$, is not supported.
:::
:::

::: function
int gr_mpoly_gen(gr_mpoly_t A, slong var, gr_mpoly_ctx_t ctx)
int gr_mpoly_gens(gr_vec_t res, gr_mpoly_ctx_t ctx)
int gr_mpoly_gens_recursive(gr_vec_t vec, gr_mpoly_ctx_t ctx)

Set `A` to the generator with index `var` (from 0), `res` to the vector
of the generators $X_1, \ldots, X_n$, or `vec` to the generators of the
coefficient ring (recursively, as constants) followed by $X_1, \ldots,
X_n$.
:::

::: function
truth_t gr_mpoly_is_zero(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)
truth_t gr_mpoly_is_one(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)
truth_t gr_mpoly_is_gen(const gr_mpoly_t A, slong var, gr_mpoly_ctx_t ctx)
truth_t gr_mpoly_equal(const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

Return whether `A` is zero, one, the generator with index `var`, or
equal to `B`, comparing terms and deciding coefficients in the
coefficient ring.
:::

## Random polynomials {#random}

::: function
int gr_mpoly_randtest_bits(gr_mpoly_t A, flint_rand_t state, slong length, flint_bitcnt_t exp_bits, gr_mpoly_ctx_t ctx)
int gr_mpoly_randtest_bound(gr_mpoly_t A, flint_rand_t state, slong length, ulong exp_bound, gr_mpoly_ctx_t ctx)
int _gr_mpoly_randtest_default(gr_mpoly_t res, flint_rand_t state, gr_mpoly_ctx_t ctx)

Set `A` to a random polynomial with up to `length` terms, with
exponents of up to `exp_bits` bits or less than `exp_bound`, and random
coefficients (`gr_randtest`). The *default* version, used for
`gr_randtest`, chooses small random parameters.
:::

## Input and output {#io}

::: function
int gr_mpoly_write_pretty(gr_stream_t out, const gr_mpoly_t A, gr_mpoly_ctx_t ctx)
int gr_mpoly_write(gr_stream_t out, gr_mpoly_t poly, gr_mpoly_ctx_t ctx)
int gr_mpoly_print_pretty(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)

Write the polynomial as an expression in the generator names of the
context, with the terms in decreasing order; `gr_mpoly_write` is the
same as `gr_mpoly_write_pretty`. Polynomials can be parsed with
`gr_set_str`.
:::

## Coefficients and terms {#terms}

::: function
int gr_mpoly_get_coeff_scalar_fmpz(gr_ptr c, const gr_mpoly_t A, const fmpz * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_get_coeff_scalar_ui(gr_ptr c, const gr_mpoly_t A, const ulong * exp, gr_mpoly_ctx_t ctx)

Set `c` to the coefficient of the monomial with exponent vector `exp`
(zero if there is no such term).
:::

::: function
int gr_mpoly_set_coeff_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_coeff_ui_fmpz(gr_mpoly_t A, ulong c, const fmpz * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_coeff_si_fmpz(gr_mpoly_t A, slong c, const fmpz * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_coeff_fmpz_fmpz(gr_mpoly_t A, const fmpz_t c, const fmpz * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_coeff_fmpq_fmpz(gr_mpoly_t A, const fmpq_t c, const fmpz * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_coeff_scalar_ui(gr_mpoly_t poly, gr_srcptr c, const ulong * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_coeff_ui_ui(gr_mpoly_t A, ulong c, const ulong * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_coeff_si_ui(gr_mpoly_t A, slong c, const ulong * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_coeff_fmpz_ui(gr_mpoly_t A, const fmpz_t c, const ulong * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_set_coeff_fmpq_ui(gr_mpoly_t A, const fmpq_t c, const ulong * exp, gr_mpoly_ctx_t ctx)

Set the coefficient of the monomial with exponent vector `exp` (given
as `fmpz` or `ulong` entries) to $c$, inserting or removing the term as
needed. Each call costs time linear in the number of terms; to build a
polynomial term by term, use the push functions below.
:::

::: function
void _gr_mpoly_push_exp_ui(gr_mpoly_t A, const ulong * exp, gr_mpoly_ctx_t ctx)
void _gr_mpoly_push_exp_fmpz(gr_mpoly_t A, const fmpz * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_push_term_scalar_ui(gr_mpoly_t A, gr_srcptr c, const ulong * exp, gr_mpoly_ctx_t ctx)
int gr_mpoly_push_term_scalar_fmpz(gr_mpoly_t A, gr_srcptr c, const fmpz * exp, gr_mpoly_ctx_t ctx)

Append a term with exponent vector `exp` at the end of `A`, increasing
the number of bits per exponent if necessary; the underscore versions
append only the monomial and leave the new coefficient for the caller
to set. The result is in general not sorted and may contain repeated
monomials or zero coefficients: call `gr_mpoly_sort_terms` and
`gr_mpoly_combine_like_terms` afterwards.
:::

::: function
void gr_mpoly_sort_terms(gr_mpoly_t A, gr_mpoly_ctx_t ctx)
int gr_mpoly_combine_like_terms(gr_mpoly_t A, gr_mpoly_ctx_t ctx)

Sort the terms of `A` by the monomial ordering, and add the
coefficients of terms with equal monomials, removing provably zero
terms; the second assumes sorted terms.
:::

::: function
truth_t gr_mpoly_is_canonical(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)
void gr_mpoly_assert_canonical(const gr_mpoly_t A, gr_mpoly_ctx_t ctx)

Check that `A` is in canonical form: terms strictly sorted, exponents
fitting in the bit size, and no provably zero coefficients.
`gr_mpoly_assert_canonical` aborts with an error message otherwise.
:::

## Arithmetic {#arithmetic}

::: function
int gr_mpoly_neg(gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_add(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)
int gr_mpoly_sub(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)

Set `A` to $-B$, $B + C$ or $B - C$, by merging the sorted term lists.
:::

::: function
int gr_mpoly_mul_scalar(gr_mpoly_t A, const gr_mpoly_t B, gr_srcptr c, gr_mpoly_ctx_t ctx)
int gr_mpoly_mul_si(gr_mpoly_t A, const gr_mpoly_t B, slong c, gr_mpoly_ctx_t ctx)
int gr_mpoly_mul_ui(gr_mpoly_t A, const gr_mpoly_t B, ulong c, gr_mpoly_ctx_t ctx)
int gr_mpoly_mul_fmpz(gr_mpoly_t A, const gr_mpoly_t B, const fmpz_t c, gr_mpoly_ctx_t ctx)
int gr_mpoly_mul_fmpq(gr_mpoly_t A, const gr_mpoly_t B, const fmpq_t c, gr_mpoly_ctx_t ctx)

Set `A` to $B c$ for a scalar $c$ (removing terms that become zero).
:::

::: function
int gr_mpoly_mul(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)
int gr_mpoly_mul_heap(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)
int gr_mpoly_mul_heap_threaded(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)
int gr_mpoly_mul_monomial(gr_mpoly_t A, const gr_mpoly_t B, const gr_mpoly_t C, gr_mpoly_ctx_t ctx)

Set `A` to the product $BC$. The *heap* version uses Johnson's
algorithm, which merges the $|B|$ sorted rows $b_i C$ with a heap of
size $|B|$, with $O(|B| |C| \log |B|)$ monomial comparisons and memory
proportional to the output; the *heap_threaded* version divides the
output among threads. The *monomial* version requires $C$ to have one
term.

::: details
The default `gr_mpoly_mul` squares when $B$ and $C$ are the same object
(see `gr_mpoly_sqr`), uses the monomial algorithm if one factor has a
single term, and otherwise the heap algorithm, threaded when $|B| |C| >
10000$, both have more than two terms and threads are available.
Products with more than the size limit of the context in terms return
`GR_UNABLE`.
:::
:::

::: function
int gr_mpoly_sqr(gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
int gr_mpoly_sqr_monomial(gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
int gr_mpoly_sqr_commutative_heap(gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)
int gr_mpoly_sqr_commutative_heap_threaded(gr_mpoly_t poly1, const gr_mpoly_t poly2, gr_mpoly_ctx_t ctx)

Set `poly1` to the square of `poly2`. The *commutative* versions
require a commutative (or approximately commutative) coefficient ring
and compute only the products $b_i b_j$ with $i \le j$, doubling the
off-diagonal ones; the default uses them (threaded above 100 terms) and
the monomial version for a single term.
:::

## Division {#division}

::: function
int gr_mpoly_divides(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_divides_heap(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_divides_heap_threaded(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

Checked exact division. If $B$ divides $A$, set $Q = A / B$ and return
`GR_SUCCESS`; if not, set $Q$ to zero and return `GR_DOMAIN`; if
divisibility cannot be decided, return `GR_UNABLE`. Division by zero
gives `GR_DOMAIN` (except in the zero ring).

::: details
The *heap* version is the heap division algorithm of Monagan and
Pearce, which stops as soon as a remainder term appears that the
leading monomial of $B$ does not divide. The *heap_threaded* version
splits the work among threads for large inputs over thread-safe rings
and otherwise calls the serial algorithm. The default handles a
single-term divisor directly, requires an integral domain for the
general case (`GR_UNABLE` otherwise, since this is also the method
behind `gr_div`), and uses the threaded version when $|A| > 500$, $|B| >
2$, the ring is thread-safe and threads are available.
:::
:::

::: function
int gr_mpoly_divexact(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_divexact_heap_threaded(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

Set $Q = A / B$, assuming without checking that $B$ divides $A$, which
saves the divisibility checks. The default is threaded under the same
conditions as `gr_mpoly_divides`.
:::

::: function
int gr_mpoly_divrem(gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_divrem_heap(gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_divrem_heap_threaded(gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_div(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_div_heap_threaded(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

Multivariate division with remainder: $A = QB + R$ where no monomial of
$R$ is divisible by the leading monomial of $B$. The *div* versions
return only $Q$. Over a field this always works; over other rings it
works whenever each division by the leading coefficient of $B$ that
occurs is exact, and returns `GR_DOMAIN` (with $Q$ and $R$ set to zero)
otherwise. The defaults use the threaded algorithm for $|A| > 500$,
$|B| > 2$ and a thread-safe ring when threads are available.
:::

::: function
int gr_mpoly_divrem_weak(gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_divrem_weak_heap_threaded(gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_div_weak(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_div_weak_heap_threaded(gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)

Weak division with remainder over rings that are not fields: when a
coefficient division by the leading coefficient of $B$ is not exact,
the coefficients are divided with remainder (`gr_euclidean_divrem` in
the coefficient ring, the analogue of `fmpz_fdiv_qr` over $\Z$) and the
remainder is kept, instead of failing. The result satisfies $A = QB + R$, but monomials of $R$ may
still be divisible by the leading monomial of $B$.
:::

::: function
int gr_mpoly_divrem_ideal(gr_mpoly_vec_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_vec_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_divrem_ideal_heap_threaded(gr_mpoly_vec_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_vec_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_divrem_ideal_weak(gr_mpoly_vec_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_vec_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_divrem_ideal_weak_heap_threaded(gr_mpoly_vec_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_vec_t B, gr_mpoly_ctx_t ctx)

Division by a vector of divisors $B_1, \ldots, B_k$: $A = \sum_i Q_i B_i
+ R$, where no monomial of $R$ is divisible by any leading monomial
$\operatorname{lm}(B_i)$ (the multivariate division algorithm,
[@CLO2015, §2.3]); the result depends on the order of the divisors
unless they form a Gröbner basis. The *weak* versions accumulate inexact
coefficient remainders as for `gr_mpoly_divrem_weak`. The output vector
`Q` is resized to $k$ entries.
:::

::: function
int gr_mpoly_quasidivrem(gr_ptr scale, gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_quasidivrem_heap(gr_ptr scale, gr_mpoly_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_quasidiv(gr_ptr scale, gr_mpoly_t Q, const gr_mpoly_t A, const gr_mpoly_t B, gr_mpoly_ctx_t ctx)
int gr_mpoly_quasidivrem_ideal(gr_ptr scale, gr_mpoly_vec_t Q, gr_mpoly_t R, const gr_mpoly_t A, const gr_mpoly_vec_t B, gr_mpoly_ctx_t ctx)

Pseudo-division: compute a scalar `scale` (an element of the coefficient
ring, a product of powers of leading coefficients of the divisors) with
$\mathit{scale} \cdot A = QB + R$, respectively $\mathit{scale} \cdot A
= \sum_i Q_i B_i + R$, and no monomial of $R$ divisible by a leading
monomial of the divisors. This divides over the fraction field of an
integral domain without leaving the ring.
:::

::: function
int gr_mpoly_inv(gr_mpoly_t res, const gr_mpoly_t poly, gr_mpoly_ctx_t ctx)

Sets `res` to the inverse of `poly` if it is a constant with an
invertible coefficient. Returns `GR_DOMAIN` if `poly` is zero (outside
the zero ring) or, over an integral domain, nonconstant; returns
`GR_UNABLE` if invertibility cannot be decided (over rings with
nilpotents a nonconstant polynomial can be a unit).
:::

::: function
int gr_mpoly_canonical_associate(gr_mpoly_t res, gr_mpoly_t u, const gr_mpoly_t poly, gr_mpoly_ctx_t ctx)

Sets `res` to the canonical associate of `poly`: `poly` multiplied by
the unit $c$ that makes its leading coefficient canonical in the
coefficient ring (`gr_canonical_associate`: monic over a field,
positive over $\Z$). If `u` is not `NULL` it is set to $c$, so that
$\mathit{res} = u \cdot \mathit{poly}$.
:::

## Derivatives and integrals {#calculus}

::: function
int gr_mpoly_derivative(gr_mpoly_t A, const gr_mpoly_t B, slong var, gr_mpoly_ctx_t ctx)

Sets `A` to the partial derivative of `B` with respect to the variable
with index `var`.
:::

::: function
int gr_mpoly_integral(gr_mpoly_t A, const gr_mpoly_t B, slong var, gr_mpoly_ctx_t ctx)

Sets `A` to the integral of `B` with respect to the variable with index
`var`, with zero constant of integration; this divides coefficients by
the new exponents (`gr_div_ui`), which must be possible in the
coefficient ring.
:::
