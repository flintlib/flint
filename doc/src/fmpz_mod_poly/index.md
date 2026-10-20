---
title: fmpz_mod_poly.h
header: fmpz_mod_poly.h
description: Dense univariate polynomials over $\Z/n\Z$ for a modulus of any size, with fast multiplication, division, GCDs, modular composition, evaluation and power series.
category: Univariate polynomials
order: 8
---

This module implements dense univariate polynomials over $\Z/n\Z$ for an
arbitrary modulus $n \ge 1$, typically a prime of more than one word.
A polynomial is an `fmpz_mod_poly_t`, an array of `fmpz` coefficients
reduced to $[0, n)$; the modulus is not stored in the polynomial but
passed to every function as an `fmpz_mod_ctx_t` (see `fmpz_mod`). For
moduli that fit in a word, `nmod_poly` is faster, and its documentation
describes the same operations in more detail. For moduli of 2 to 16
limbs, the generic polynomials `gr_poly` over the ring `mpn_mod`
(`gr_ctx_init_mpn_mod`) are usually faster still; see
[](algorithms/arithmetic/multiword-overview).

Most of the arithmetic is delegated: products are computed over $\Z$ by
`fmpz_poly` (Kronecker substitution, Schönhage–Strassen or multimodular
FFTs, depending on the sizes) and reduced afterwards, while division,
power series, GCDs, resultants, powering modulo a polynomial,
composition and modular composition are thin wrappers around the generic
algorithms of [`gr_poly`](gr_poly) over the `fmpz_mod` context. These
algorithms are quasi-linear:
Newton iteration for division and series ([](algorithms/polynomials/newton)),
the half-GCD algorithm for GCDs and resultants
([](algorithms/polynomials/half-gcd)) and the Brent–Kung algorithm for
modular composition ([](algorithms/polynomials/composition)).
Factorisation, irreducibility tests and root finding are in
`fmpz_mod_poly_factor.h`.

The guide [](fmpz_mod_poly/getting-started) shows how to create and use
polynomials and how to work with a modulus that may be composite.

## Types and conventions {#conventions}

::: type
fmpz_mod_poly_struct
fmpz_mod_poly_t

A polynomial over $\Z/n\Z$. The structure has the fields `coeffs` (an
array of `alloc` coefficients of type `fmpz`), `alloc` and `length`, in
the same layout as `gr_poly_struct`. An `fmpz_mod_poly_t` is an array of
one `fmpz_mod_poly_struct`, so it is passed by reference. The coefficient
of $x^i$ is `coeffs[i]`.
:::

- All polynomials passed to one function must belong to the context
  `ctx`. Some functions do not use `ctx` but take it for uniformity.
- A polynomial is *normalised* if its length is zero or its leading
  coefficient is nonzero, and its coefficients are *reduced* if they lie
  in $[0, n)$. All functions require normalised inputs with reduced
  coefficients and produce such outputs, unless stated otherwise.
- Functions without a leading underscore manage memory and allow every
  output to be the same object as any input, unless stated otherwise.
  Functions with a leading underscore work on coefficient arrays, assume
  nonzero lengths and enough space in the output, and do not allow
  aliasing unless stated otherwise.
- Division by a polynomial needs an invertible leading coefficient, and
  GCDs, modular inverses, minimal polynomials and square roots need a
  prime modulus. This is not checked; a non-invertible element usually
  raises an exception.
- The functions with suffix `_f` support a composite modulus $n$: they
  take an `fmpz_t f` as first argument and either set $f = 1$ and
  compute the result, or, if they meet a non-invertible leading
  coefficient, set $f$ to a nontrivial factor of $n$ (and leave the other
  outputs undefined).

## Memory management {#memory}

::: function
void fmpz_mod_poly_init(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_init2(fmpz_mod_poly_t poly, slong alloc, const fmpz_mod_ctx_t ctx)

Initialise `poly` as the zero polynomial. `fmpz_mod_poly_init` allocates
nothing, `fmpz_mod_poly_init2` allocates space for `alloc` coefficients.
:::

::: function
void fmpz_mod_poly_clear(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Frees the memory used by `poly`, which must be initialised again before
it is used.
:::

::: function
void fmpz_mod_poly_realloc(fmpz_mod_poly_t poly, slong alloc, const fmpz_mod_ctx_t ctx)

Changes the allocated space of `poly` to `alloc` coefficients. If the
length exceeds `alloc`, the polynomial is truncated to length `alloc`
and normalised. For `alloc` = 0 the polynomial is cleared and
reinitialised.
:::

::: function
void fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, slong len, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_fit_length(fmpz_mod_poly_t poly, slong len)

Ensure that `poly` has space for at least `len` coefficients. If it has to
grow, it grows to at least twice the old allocation, so that repeated
calls cost amortised linear time. The length and coefficients are not
changed.
:::

::: function
void _fmpz_mod_poly_normalise(fmpz_mod_poly_t poly)

Decreases the length of `poly` until the leading coefficient is nonzero or
the length is zero.
:::

::: function
void _fmpz_mod_poly_set_length(fmpz_mod_poly_t poly, slong len)

Sets the length of `poly` to `len`. If `len` is smaller than the current
length, the discarded coefficients are set to zero (freeing their
memory); otherwise the coefficients between the old and the new length
must already be valid. The result is not normalised.
:::

::: function
void fmpz_mod_poly_truncate(fmpz_mod_poly_t poly, slong len, const fmpz_mod_ctx_t ctx)

Truncates `poly` to length `len` (if it is longer) and normalises it.
:::

::: function
void fmpz_mod_poly_set_trunc(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, slong n, const fmpz_mod_ctx_t ctx)

Sets `res` to `poly` truncated to length $n$, that is, to $\mathit{poly}
\bmod x^n$.
:::

## Length, degree and basic properties {#properties}

::: function
slong fmpz_mod_poly_length(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
slong fmpz_mod_poly_degree(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Return the length and the degree of `poly`. The zero polynomial has
length 0 and degree $-1$.
:::

::: function
fmpz * fmpz_mod_poly_lead(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Returns a pointer to the leading coefficient of `poly`, or `NULL` if
`poly` is zero.
:::

::: function
int fmpz_mod_poly_is_zero(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_is_one(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_is_gen(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Return whether `poly` is 0, 1, or the generator $x$. Modulo $n = 1$ every
polynomial is zero; `fmpz_mod_poly_is_gen` then returns 1, while
`fmpz_mod_poly_is_one` returns 0.
:::

::: function
int fmpz_mod_poly_is_monic(const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Returns whether `f` is nonzero with leading coefficient 1.
:::

::: function
int fmpz_mod_poly_is_unit(const fmpz_mod_poly_t op, const fmpz_mod_ctx_t ctx)

Returns whether `op` is a constant that is invertible modulo $n$ (or
whether $n = 1$). For prime $n$ these are exactly the units of
$(\Z/n\Z)[x]$; for composite $n$ there are also nonconstant units, such
as $1 + 2x$ modulo 4, which this function does not recognise.
:::

::: function
int fmpz_mod_poly_is_canonical(const fmpz_mod_poly_t A, const fmpz_mod_ctx_t ctx)

Returns whether `A` is normalised and has all its coefficients reduced
modulo $n$. Useful for checking the result of low-level manipulations.
:::

::: function
slong fmpz_mod_poly_hamming_weight(const fmpz_mod_poly_t A, const fmpz_mod_ctx_t ctx)

Returns the number of nonzero coefficients of `A`.
:::

::: function
int fmpz_mod_poly_equal(const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_equal_trunc(const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx)

Return whether `poly1` and `poly2` are equal, respectively equal modulo
$x^n$.
:::

## Randomisation {#randomisation}

::: function
void fmpz_mod_poly_rand(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_rand_monic(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)

Set `f` to a random polynomial of length at most `len` with uniformly
distributed coefficients; `fmpz_mod_poly_rand_monic` makes the
coefficient of $x^{\mathit{len} - 1}$ equal to 1 ($\mathit{len} > 0$).
:::

::: function
void fmpz_mod_poly_rand_irreducible(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)

Sets `f` to a random irreducible polynomial of length at most `len` with
uniformly distributed coefficients, by testing random polynomials with
`fmpz_mod_poly_is_irreducible`. Requires $\mathit{len} \ge 1$ and a prime
modulus.
:::

::: function
void fmpz_mod_poly_randtest(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_randtest_not_zero(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_randtest_monic(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)

Set `f` to a random polynomial of length at most `len`, nonzero
($\mathit{len} \ge 1$), or monic of length exactly `len`. The
coefficients are drawn with `fmpz_randtest_mod`, which favours special
values (small values, values close to $n$), as appropriate for testing.
:::

::: function
void fmpz_mod_poly_randtest_irreducible(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_randtest_monic_irreducible(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)

Set `f` to a random irreducible polynomial of length at most `len`
(respectively a monic one of length exactly `len`), for testing.
Requires $\mathit{len} \ge 1$ and a prime modulus.
:::

::: function
void fmpz_mod_poly_randtest_monic_primitive(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)

Sets `f` to a random monic irreducible polynomial of length `len` whose
roots generate the multiplicative group of $\mathbb{F}_{p^{\mathit{len} -
1}}$ (a primitive polynomial), tested with `fq_is_primitive`.
:::

::: function
void fmpz_mod_poly_randtest_trinomial(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_randtest_pentomial(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)

Set `f` to a random monic polynomial of length `len` of the form $x^{d} +
a x^k + b$ with $0 < k < d$ (trinomial, $\mathit{len} \ge 3$), or
$x^{d} + a_3 x^3 + a_2 x^2 + a_1 x + a_0$ (pentanomial, $\mathit{len} \ge
5$), where $d = \mathit{len} - 1$ and some of the random coefficients may
be zero.
:::

::: function
int fmpz_mod_poly_randtest_trinomial_irreducible(fmpz_mod_poly_t f, flint_rand_t state, slong len, slong max_attempts, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_randtest_pentomial_irreducible(fmpz_mod_poly_t f, flint_rand_t state, slong len, slong max_attempts, const fmpz_mod_ctx_t ctx)

Try up to `max_attempts` random trinomials or pentanomials of length
`len` (without limit if `max_attempts` is 0) and return 1 as soon as one
is irreducible, with `f` set to it; return 0 if none was found.
:::

::: function
void fmpz_mod_poly_randtest_sparse_irreducible(fmpz_mod_poly_t f, flint_rand_t state, slong len, const fmpz_mod_ctx_t ctx)

Sets `f` to a sparse monic irreducible polynomial of length `len`: an
irreducible trinomial if one is found in $2 \mathit{len}$ attempts,
otherwise a pentanomial (likewise), otherwise a random monic irreducible
polynomial with few terms.
:::

## Assignment and basic manipulation {#assignment}

::: function
void fmpz_mod_poly_set(fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

Sets `poly1` to `poly2`.
:::

::: function
void fmpz_mod_poly_swap(fmpz_mod_poly_t poly1, fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

Swaps `poly1` and `poly2` efficiently, by swapping the structures.
:::

::: function
void fmpz_mod_poly_zero(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_one(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_gen(fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Set `poly` to 0, 1 or the generator $x$ (all three are zero for $n =
1$).
:::

::: function
void fmpz_mod_poly_zero_coeffs(fmpz_mod_poly_t poly, slong i, slong j, const fmpz_mod_ctx_t ctx)

Sets the coefficients of $x^k$ for $i \le k < j$ to zero, normalising the
result if the top coefficients were cleared. Indices outside $[0,
\mathit{len})$ are ignored.
:::

::: function
void _fmpz_mod_poly_reverse(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_mod_poly_reverse(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, slong n, const fmpz_mod_ctx_t ctx)

Set `res` to the reversal $x^{n-1} \mathit{poly}(1/x)$ of `poly` viewed
as a polynomial of length $n$ (truncated or zero-padded). The
underscore version writes exactly $n$ coefficients, requires $\mathit{len}
\le n$ and allows `res` and `poly` to be the same array; the other
version normalises the result, which can therefore be shorter than $n$.
:::

## Conversion {#conversion}

::: function
void fmpz_mod_poly_set_ui(fmpz_mod_poly_t f, ulong x, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_set_fmpz(fmpz_mod_poly_t poly, const fmpz_t c, const fmpz_mod_ctx_t ctx)

Set the polynomial to the constant `x` or `c`, reduced modulo $n$.
:::

::: function
void fmpz_mod_poly_set_fmpz_poly(fmpz_mod_poly_t f, const fmpz_poly_t g, const fmpz_mod_ctx_t ctx)

Sets `f` to the integer polynomial `g` reduced modulo $n$.
:::

::: function
void fmpz_mod_poly_get_fmpz_poly(fmpz_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx)

Sets `f` to the integer polynomial whose coefficients are the
representatives in $[0, n)$ of the coefficients of `g`.
:::

::: function
void fmpz_mod_poly_set_nmod_poly(fmpz_mod_poly_t f, const nmod_poly_t g)
void fmpz_mod_poly_get_nmod_poly(nmod_poly_t f, const fmpz_mod_poly_t g)

Copy the coefficients between an `fmpz_mod_poly_t` and an `nmod_poly_t`,
assuming that both have the same modulus, which is not checked (no
context is passed). For `fmpz_mod_poly_get_nmod_poly`, the modulus of `f`
must already be set.
:::

## Getting and setting coefficients {#coefficients}

::: function
void fmpz_mod_poly_set_coeff_fmpz(fmpz_mod_poly_t poly, slong n, const fmpz_t x, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_set_coeff_ui(fmpz_mod_poly_t poly, slong n, ulong x, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_set_coeff_si(fmpz_mod_poly_t poly, slong n, slong x, const fmpz_mod_ctx_t ctx)

Set the coefficient of $x^n$, $n \ge 0$, to `x` reduced modulo the
modulus, extending or normalising the polynomial as needed.
:::

::: function
void fmpz_mod_poly_get_coeff_fmpz(fmpz_t x, const fmpz_mod_poly_t poly, slong n, const fmpz_mod_ctx_t ctx)

Sets `x` to the coefficient of $x^n$ in `poly` (zero if $n$ is at least
the length).
:::

## Shifting {#shifting}

::: function
void _fmpz_mod_poly_shift_left(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_mod_poly_shift_left(fmpz_mod_poly_t f, const fmpz_mod_poly_t g, slong n, const fmpz_mod_ctx_t ctx)

Set the output to the input multiplied by $x^n$. The underscore version
writes `len + n` coefficients and allows `res` and `poly` to be the same
array.
:::

::: function
void _fmpz_mod_poly_shift_right(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_mod_poly_shift_right(fmpz_mod_poly_t f, const fmpz_mod_poly_t g, slong n, const fmpz_mod_ctx_t ctx)

Set the output to the input divided by $x^n$, discarding the low
coefficients; the result is zero if $n$ is at least the length. The
underscore version writes `len - n` coefficients, requires $\mathit{len}
> n$ and allows `res` and `poly` to be the same array.
:::

## Inflation and deflation {#deflation}

::: function
void fmpz_mod_poly_inflate(fmpz_mod_poly_t result, const fmpz_mod_poly_t input, ulong inflation, const fmpz_mod_ctx_t ctx)

Sets `result` to $p(x^m)$ where $p$ is `input` and $m$ = `inflation`.
:::

::: function
void fmpz_mod_poly_deflate(fmpz_mod_poly_t result, const fmpz_mod_poly_t input, ulong deflation, const fmpz_mod_ctx_t ctx)

Sets `result` to $p(x^{1/m})$ where $p$ is `input` and $m \ge 1$ =
`deflation`, which must divide the exponents of all nonzero terms of
`input` (see `fmpz_mod_poly_deflation`). $m = 0$ raises an exception.
:::

::: function
ulong fmpz_mod_poly_deflation(const fmpz_mod_poly_t input, const fmpz_mod_ctx_t ctx)

Returns the largest $m$ such that `input` is a polynomial in $x^m$, that
is, the gcd of the exponents of its nonzero terms; 0 for the zero
polynomial and 1 for a nonzero constant.
:::

## Addition and subtraction {#addition}

::: function
void _fmpz_mod_poly_add(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_sub(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz_mod_ctx_t ctx)

Set `res` to the sum or difference of `(poly1, len1)` and `(poly2,
len2)`, writing $\max(\mathit{len1}, \mathit{len2})$ coefficients without
normalising. `res` may be one of the inputs.
:::

::: function
void fmpz_mod_poly_add(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_sub(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

Set `res` to $\mathit{poly1} + \mathit{poly2}$, respectively
$\mathit{poly1} - \mathit{poly2}$.
:::

::: function
void fmpz_mod_poly_add_series(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_sub_series(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx)

Set `res` to the sum or difference of `poly1` and `poly2` truncated to
length $n$.
:::

::: function
void fmpz_mod_poly_add_si(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, slong c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_add_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_sub_si(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, slong c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_sub_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_si_sub(fmpz_mod_poly_t res, slong c, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_fmpz_sub(fmpz_mod_poly_t res, const fmpz_t c, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Set `res` to $\mathit{poly} + c$, $\mathit{poly} - c$ or $c -
\mathit{poly}$ for an integer constant $c$, which is reduced modulo $n$.
:::

::: function
void _fmpz_mod_poly_neg(fmpz * res, const fmpz * poly, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_neg(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Set `res` to $-\mathit{poly}$.
:::

## Scalar multiplication and division {#scalar}

::: function
void _fmpz_mod_poly_scalar_mul_fmpz(fmpz * res, const fmpz * poly, slong len, const fmpz_t x, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_scalar_mul_ui(fmpz * res, const fmpz * poly, slong len, ulong x, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_scalar_mul_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_t x, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_scalar_mul_ui(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, ulong x, const fmpz_mod_ctx_t ctx)

Set `res` to $x \cdot \mathit{poly}$, with $x$ reduced modulo $n$ (it
need not be reduced on input).
:::

::: function
void fmpz_mod_poly_scalar_addmul_fmpz(fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_t x, const fmpz_mod_ctx_t ctx)

Sets $A$ to $A + x B$.
:::

::: function
void _fmpz_mod_poly_scalar_div_fmpz(fmpz * res, const fmpz * poly, slong len, const fmpz_t x, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_scalar_div_fmpz(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_t x, const fmpz_mod_ctx_t ctx)

Set `res` to $x^{-1} \cdot \mathit{poly}$. If $x$ is not invertible
modulo $n$, an exception is raised (`FLINT_DIVZERO` for $x \equiv 0$,
`FLINT_ERROR` otherwise).
:::

## Multiplication {#multiplication}

::: function
void _fmpz_mod_poly_mul(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_mul(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

Set `res` to the product of `poly1` and `poly2`. The underscore version
writes $\mathit{len1} + \mathit{len2} - 1$ coefficients, requires positive
lengths (in either order) and allows zero-padded inputs.

::: details
The product is computed over $\Z$ with `_fmpz_poly_mul` and then
reduced; for large inputs this uses Kronecker substitution or FFT
multiplication of integer polynomials, so its cost is roughly that of
one integer product of about $(\mathit{len1} + \mathit{len2})(2 \log_2 n +
\log_2 \mathit{len})$ bits.
:::
:::

::: function
void _fmpz_mod_poly_mullow(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_mullow(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong n, const fmpz_mod_ctx_t ctx)

Set `res` to the product of `poly1` and `poly2` truncated to length $n$.
The underscore version requires positive lengths and $0 < n \le
\mathit{len1} + \mathit{len2} - 1$, and does not allow aliasing.
:::

::: function
void _fmpz_mod_poly_mulmid(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_mulmid(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong nlo, slong nhi, const fmpz_mod_ctx_t ctx)

Set `res` to the coefficients of index $\mathit{nlo} \le i < \mathit{nhi}$
of the product, that is, to $((\mathit{poly1} \cdot \mathit{poly2}) \bmod
x^{\mathit{nhi}}) / x^{\mathit{nlo}}$. The underscore version requires
positive lengths and $0 \le \mathit{nlo} < \mathit{nhi} \le \mathit{len1}
+ \mathit{len2} - 1$, writes $\mathit{nhi} - \mathit{nlo}$ coefficients
and does not allow aliasing.
:::

::: function
void fmpz_mod_poly_mulhigh(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, slong start, const fmpz_mod_ctx_t ctx)

Sets the coefficients of index at least `start` of `res` to those of the
product of `poly1` and `poly2`. The coefficients of index less than
`start` are arbitrary (they are not set to zero).
:::

::: function
void _fmpz_mod_poly_sqr(fmpz * res, const fmpz * poly, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_sqr(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Set `res` to the square of `poly`; the underscore version writes $2
\mathit{len} - 1$ coefficients.
:::

::: function
void _fmpz_mod_poly_mulmod(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz * f, slong lenf, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_mulmod(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Set `res` to $\mathit{poly1} \cdot \mathit{poly2} \bmod f$. The
underscore version requires $\mathit{len1} + \mathit{len2} - \mathit{lenf}
> 0$ (so that a reduction takes place), writes $\mathit{lenf} - 1$
coefficients and does not allow `res` to alias `f`. The leading
coefficient of $f$ must be invertible (`FLINT_IMPINV` otherwise).
:::

::: function
void _fmpz_mod_poly_mulmod_preinv(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz * f, slong lenf, const fmpz * finv, slong lenfinv, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_mulmod_preinv(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, const fmpz_mod_ctx_t ctx)

Like `fmpz_mod_poly_mulmod`, given `finv`, the inverse of the reversal of
$f$ modulo $x^{\mathit{lenf}}$ (computed with `fmpz_mod_poly_reverse`
and `fmpz_mod_poly_inv_series`). Both factors must be reduced modulo $f$
(`len1`, `len2` $< \mathit{lenf}$); otherwise an exception is raised.
:::

## Powering {#powering}

::: function
void _fmpz_mod_poly_pow(fmpz * rop, const fmpz * op, slong len, ulong e, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_pow(fmpz_mod_poly_t rop, const fmpz_mod_poly_t op, ulong e, const fmpz_mod_ctx_t ctx)

Set `rop` to $\mathit{op}^e$, with $0^0 = 1$. The underscore version
requires $e > 1$ and $\mathit{len} > 0$, writes $e(\mathit{len} - 1) + 1$
coefficients and does not allow aliasing.
:::

::: function
void _fmpz_mod_poly_pow_trunc(fmpz * res, const fmpz * poly, ulong e, slong trunc, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_pow_trunc(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, ulong e, slong trunc, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_pow_trunc_binexp(fmpz * res, const fmpz * poly, ulong e, slong trunc, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_pow_trunc_binexp(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, ulong e, slong trunc, const fmpz_mod_ctx_t ctx)

Set `res` to $\mathit{poly}^e \bmod x^{\mathit{trunc}}$ by binary
exponentiation with truncated products. The underscore versions require
`poly` to be zero-padded to length `trunc`, $\mathit{trunc} > 0$ and $e >
1$, and do not allow aliasing.
:::

::: function
void _fmpz_mod_poly_powmod_ui_binexp(fmpz * res, const fmpz * poly, ulong e, const fmpz * f, slong lenf, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_powmod_ui_binexp(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, ulong e, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_powmod_fmpz_binexp(fmpz * res, const fmpz * poly, const fmpz_t e, const fmpz * f, slong lenf, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_powmod_fmpz_binexp(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_t e, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Set `res` to $\mathit{poly}^e \bmod f$ by binary exponentiation, for $e
\ge 0$. The underscore versions require $e > 0$, $\mathit{lenf} > 1$ and
`poly` reduced modulo $f$ and zero-padded to length $\mathit{lenf} - 1$;
they write $\mathit{lenf} - 1$ coefficients.
:::

::: function
void _fmpz_mod_poly_powmod_ui_binexp_preinv(fmpz * res, const fmpz * poly, ulong e, const fmpz * f, slong lenf, const fmpz * finv, slong lenfinv, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_powmod_ui_binexp_preinv(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, ulong e, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_powmod_fmpz_binexp_preinv(fmpz * res, const fmpz * poly, const fmpz_t e, const fmpz * f, slong lenf, const fmpz * finv, slong lenfinv, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_powmod_fmpz_binexp_preinv(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_t e, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, const fmpz_mod_ctx_t ctx)

Like the previous functions, given `finv`, the inverse of the reversal
of $f$ modulo $x^{\mathit{lenf}}$, so that each reduction modulo $f$
costs two multiplications.
:::

::: function
void _fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz * res, const fmpz_t e, const fmpz * f, slong lenf, const fmpz * finv, slong lenfinv, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_powmod_x_fmpz_preinv(fmpz_mod_poly_t res, const fmpz_t e, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, const fmpz_mod_ctx_t ctx)

Set `res` to $x^e \bmod f$ for $e \ge 0$ by sliding-window
exponentiation, given `finv` as above. Multiplications by $x$ are shifts,
which makes this cheaper than powering a general polynomial. The
underscore version requires $e > 0$ and $\mathit{lenf} > 2$ and writes
$\mathit{lenf} - 1$ coefficients.
:::

::: function
void fmpz_mod_poly_powmod_linear_fmpz_preinv(fmpz_mod_poly_t res, const fmpz_t a, const fmpz_t e, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, const fmpz_mod_ctx_t ctx)

Sets `res` to $(x + a)^e \bmod f$ for $e \ge 0$, given `finv` as above,
by binary exponentiation in which the multiplications by $x + a$ are
done in linear time. Used in root finding. A negative $e$ raises an
exception.
:::

::: function
void _fmpz_mod_poly_powers_mod_preinv_naive(fmpz ** res, const fmpz * f, slong flen, slong n, const fmpz * g, slong glen, const fmpz * ginv, slong ginvlen, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_powers_mod_naive(fmpz_mod_poly_struct * res, const fmpz_mod_poly_t f, slong n, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx)

Set the $n$ entries of `res` to $f^0, f^1, \ldots, f^{n-1} \bmod g$ by
repeated multiplication. In the underscore version $f$ must be reduced
modulo $g$ (with length $\mathit{flen} > 0$, possibly zero-padded), each
`res[i]` must have space for $\mathit{glen} - 1$ coefficients, and `ginv`
is the inverse of the reversal of $g$. The entries of `res` must not
alias the inputs.
:::

::: function
void _fmpz_mod_poly_powers_mod_preinv_threaded_pool(fmpz ** res, const fmpz * f, slong flen, slong n, const fmpz * g, slong glen, const fmpz * ginv, slong ginvlen, const fmpz_mod_ctx_t ctx, thread_pool_handle * threads, slong num_threads)
void fmpz_mod_poly_powers_mod_bsgs(fmpz_mod_poly_struct * res, const fmpz_mod_poly_t f, slong n, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx)

Like the previous functions, computing the powers with a baby-step
giant-step scheme in which blocks of consecutive powers are computed in
parallel using the threads of the given thread pool (respectively the
threads available through `flint_get_num_threads`).
:::

::: type
fmpz_mod_poly_frobenius_powers_2exp_struct
fmpz_mod_poly_frobenius_powers_2exp_t
fmpz_mod_poly_frobenius_powers_struct
fmpz_mod_poly_frobenius_powers_t

Tables of Frobenius powers $x^{p^k} \bmod f$, with fields `pow` (an array
of polynomials) and `len`.
:::

::: function
void fmpz_mod_poly_frobenius_powers_2exp_precomp(fmpz_mod_poly_frobenius_powers_2exp_t pow, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, ulong m, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_frobenius_powers_2exp_clear(fmpz_mod_poly_frobenius_powers_2exp_t pow, const fmpz_mod_ctx_t ctx)

Compute (and free) the table of $x^{p^{2^i}} \bmod f$ for $2^i \le m$,
where $p$ is the modulus, given `finv`, the inverse of the reversal of
$f$. The first entry is computed by exponentiation and each following one
by modular composition of the previous one with itself.
:::

::: function
void fmpz_mod_poly_frobenius_power(fmpz_mod_poly_t res, fmpz_mod_poly_frobenius_powers_2exp_t pow, const fmpz_mod_poly_t f, ulong m, const fmpz_mod_ctx_t ctx)

Sets `res` to $x^{p^m} \bmod f$ using the table `pow` from
`fmpz_mod_poly_frobenius_powers_2exp_precomp` (computed for a bound at
least $m$), by modular compositions of the table entries along the binary
expansion of $m$ (see
[](algorithms/polynomials/factoring-finite-fields#frobenius)).
:::

::: function
void fmpz_mod_poly_frobenius_powers_precomp(fmpz_mod_poly_frobenius_powers_t pow, const fmpz_mod_poly_t f, const fmpz_mod_poly_t finv, ulong m, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_frobenius_powers_clear(fmpz_mod_poly_frobenius_powers_t pow, const fmpz_mod_ctx_t ctx)

Compute (and free) the table of all $x^{p^k} \bmod f$ for $0 \le k \le m$,
given `finv`.
:::

## Division {#division}

::: function
void _fmpz_mod_poly_divrem(fmpz * Q, fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_divrem(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Compute $Q$ and $R$ with $A = BQ + R$ and $\deg R < \deg B$. The leading
coefficient of $B$ must be invertible; `invB` is its inverse. The
underscore version requires $\mathit{lenA} \ge \mathit{lenB} > 0$, writes
$\mathit{lenA} - \mathit{lenB} + 1$ coefficients of $Q$ and $\mathit{lenB}
- 1$ of $R$, allows zero-padding of $A$ and no aliasing. A zero $B$ raises
`FLINT_DIVZERO` (except modulo 1).

::: details
Small cases (`lenB` at most 30 or quotients of length at most 11) use
the classical algorithm; larger ones use Newton iteration
(`_gr_poly_divrem_newton`, see
[](algorithms/polynomials/newton#divrem)), in $O(M(\mathit{lenA}))$
operations.
:::
:::

::: function
void _fmpz_mod_poly_div(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_div(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Set $Q$ to the quotient of the division of $A$ by $B$, with the same
conventions as the division with remainder.
:::

::: function
void _fmpz_mod_poly_rem(fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_rem(fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Set $R$ to $A \bmod B$, with the same conventions as the division with
remainder.
:::

::: function
void _fmpz_mod_poly_divrem_basecase(fmpz * Q, fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_divrem_basecase(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_rem_basecase(fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_rem_basecase(fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Division with remainder, respectively remainder, by the classical
algorithm in $O(\mathit{lenB} (\mathit{lenA} - \mathit{lenB} + 1))$
operations. In the underscore versions $R$ (of length $\mathit{lenA}$ for
`_fmpz_mod_poly_divrem_basecase`, $\mathit{lenB} - 1$ for
`_fmpz_mod_poly_rem_basecase`) may alias $A$, but no other aliasing is
allowed.
:::

::: function
void _fmpz_mod_poly_divrem_newton_n_preinv(fmpz * Q, fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz * Binv, slong lenBinv, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_divrem_newton_n_preinv(fmpz_mod_poly_t Q, fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_poly_t Binv, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_div_newton_n_preinv(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz * Binv, slong lenBinv, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_div_newton_n_preinv(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_poly_t Binv, const fmpz_mod_ctx_t ctx)

Division with remainder, respectively quotient only, given `Binv`, the
inverse of the reversal of $B$ modulo $x^{\mathit{lenB}}$: the quotient is
a truncated product of the reversed $A$ with `Binv`, and the remainder
costs one more product. Requires $\mathit{lenA} \le 2 \mathit{lenB} - 2$;
the non-underscore versions raise an exception otherwise.
:::

::: function
void _fmpz_mod_poly_divrem_f(fmpz_t f, fmpz * Q, fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_divrem_f(fmpz_t f, fmpz_mod_poly_t Q, fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_rem_f(fmpz_t f, fmpz_mod_poly_t R, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Division with remainder (respectively remainder) for a possibly
composite modulus: if the leading coefficient of $B$ is invertible, sets
$f = 1$ and computes the result; otherwise sets $f$ to the nontrivial
factor $\gcd(\mathrm{lc}(B), n)$ of $n$ and leaves $Q$ and $R$
unchanged. A zero $B$ raises a division by zero exception.
:::

::: function
ulong fmpz_mod_poly_remove(fmpz_mod_poly_t f, const fmpz_mod_poly_t p, const fmpz_mod_ctx_t ctx)

Divides `f` by the largest power of `p` that divides it and returns the
exponent. `p` must have an invertible leading coefficient. An exception
is raised if `p` is zero or constant.
:::

## Divisibility testing {#divisibility}

::: function
int _fmpz_mod_poly_divides(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_divides(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
int _fmpz_mod_poly_divides_classical(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_divides_classical(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Return 1 and set $Q = A/B$ if $B$ divides $A$, and otherwise return 0 and
set $Q$ to zero. The leading coefficient of $B$ must be invertible. The
underscore versions require $\mathit{lenA} \ge \mathit{lenB} > 0$ and
write $\mathit{lenA} - \mathit{lenB} + 1$ coefficients.

::: details
The quotient is computed first, and the product $BQ$ is compared with $A$
coefficient by coefficient from the bottom, so that a non-divisor is
usually detected early; the classical version does the comparison with
dot products and is used by `fmpz_mod_poly_divides` for
$\mathit{lenA} < 40$ and $\mathit{lenB} < 20$.
:::
:::

## Power series {#series}

::: function
void _fmpz_mod_poly_inv_series(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_inv_series(fmpz_mod_poly_t Qinv, const fmpz_mod_poly_t Q, slong n, const fmpz_mod_ctx_t ctx)

Set `Qinv` to the inverse of the power series $Q$ modulo $x^n$, $n \ge 1$,
by Newton iteration (`_gr_poly_inv_series`). The constant term of $Q$
must be invertible: zero raises `FLINT_DIVZERO`, a non-invertible nonzero
value (composite modulus) another exception. The underscore version
writes $n$ coefficients and does not allow aliasing.
:::

::: function
void fmpz_mod_poly_inv_series_f(fmpz_t f, fmpz_mod_poly_t Qinv, const fmpz_mod_poly_t Q, slong n, const fmpz_mod_ctx_t ctx)

Like `fmpz_mod_poly_inv_series` for a possibly composite modulus: sets $f
= 1$ and `Qinv` to the inverse, or sets $f$ to a nontrivial factor of $n$
(when the constant term is not invertible) with `Qinv` undefined.
:::

::: function
void _fmpz_mod_poly_div_series(fmpz * Q, const fmpz * A, slong Alen, const fmpz * B, slong Blen, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_div_series(fmpz_mod_poly_t Q, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, slong n, const fmpz_mod_ctx_t ctx)

Set $Q$ to $A/B \bmod x^n$. The constant term of $B$ must be invertible
(an exception is raised otherwise). The underscore version requires
$\mathit{Alen}, \mathit{Blen} \le n$.
:::

::: function
void _fmpz_mod_poly_sqrt_series(fmpz * g, const fmpz * h, slong hlen, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_sqrt_series(fmpz_mod_poly_t g, const fmpz_mod_poly_t h, slong n, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_invsqrt_series(fmpz * g, const fmpz * h, slong hlen, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_invsqrt_series(fmpz_mod_poly_t g, const fmpz_mod_poly_t h, slong n, const fmpz_mod_ctx_t ctx)

Set $g$ to $\sqrt{h}$, respectively $1/\sqrt{h}$, modulo $x^n$, the series
being defined by the binomial expansion of $(1 + y)^{\pm 1/2}$. The
constant term of $h$ must be 1 (an exception is raised otherwise), and 2
must be invertible modulo $n$. The underscore versions require $n > 0$
and $\mathit{hlen} > 0$ and do not allow aliasing.
:::

## Greatest common divisors and modular inverses {#gcd}

::: function
void fmpz_mod_poly_make_monic(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Sets `res` to `poly` divided by its leading coefficient (zero if `poly`
is zero). A non-invertible leading coefficient raises `FLINT_IMPINV`.
:::

::: function
void fmpz_mod_poly_make_monic_f(fmpz_t f, fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Either sets $f = 1$ and `res` to `poly` made monic, or sets $f$ to a
nontrivial factor of $n$ (the gcd of $n$ with the leading coefficient)
and leaves `res` undefined.
:::

::: function
slong _fmpz_mod_poly_gcd(fmpz * G, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_gcd(fmpz_mod_poly_t G, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Set $G$ to the greatest common divisor of $A$ and $B$, which is made
monic by the non-underscore version (the gcd of two zero polynomials is
zero). The underscore version requires $\mathit{lenA} \ge \mathit{lenB} >
0$ and space for $\mathit{lenB}$ coefficients in $G$, returns the length
of $G$ and does not normalise its leading coefficient. The modulus must
be prime.

::: details
For $\min(\mathit{lenA}, \mathit{lenB}) < 256$ the Euclidean algorithm is
used (`_gr_poly_gcd_euclidean`), and above that the half-GCD algorithm
(`_gr_poly_gcd_hgcd` with base case cutoff 128), in $O(M(n) \log n)$
operations; see [](algorithms/polynomials/half-gcd).
:::
:::

::: function
slong _fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz * G, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_gcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
slong _fmpz_mod_poly_gcd_f(fmpz_t f, fmpz * G, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_gcd_f(fmpz_t f, fmpz_mod_poly_t G, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

GCD for a possibly composite modulus by the Euclidean algorithm: either
sets $f = 1$ and computes $G$ as in `fmpz_mod_poly_gcd` (the underscore
versions returning its length), or sets $f$ to a nontrivial factor of $n$
met as a non-invertible leading coefficient of a remainder, with $G$
undefined (the underscore versions then return 0). The `_gcd_f` names
are aliases of the `_gcd_euclidean_f` functions. The underscore versions
require $\mathit{lenA} \ge \mathit{lenB} > 0$ and do not allow aliasing.
:::

::: function
slong _fmpz_mod_poly_hgcd(fmpz ** M, slong * lenM, fmpz * A, slong * lenA, fmpz * B, slong * lenB, const fmpz * a, slong lena, const fmpz * b, slong lenb, const fmpz_mod_ctx_t ctx)

One step of the half-GCD algorithm (`_gr_poly_hgcd`): for $\mathit{lena} >
\mathit{lenb} > 0$, computes the $2 \times 2$ matrix $M$ (the entries
`M[0]`, ..., `M[3]` with lengths `lenM[0]`, ..., `lenM[3]`, each with space
for `lena` coefficients) and the two consecutive remainders $A$, $B$ of
the Euclidean remainder sequence of $(a, b)$ that straddle half the
degree of $a$, such that $(A, B)^T = \sigma M^{-1} (a, b)^T$, and returns
the sign $\sigma = \pm 1$. $A$ and $B$ need space for `lena` and `lenb`
coefficients; their lengths are returned in `*lenA`, `*lenB`. See
[](algorithms/polynomials/half-gcd).
:::

::: function
slong _fmpz_mod_poly_xgcd(fmpz * G, fmpz * S, fmpz * T, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_xgcd(fmpz_mod_poly_t G, fmpz_mod_poly_t S, fmpz_mod_poly_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Compute the GCD $G$ of $A$ and $B$ and cofactors with $SA + TB = G$. The
non-underscore version makes $G$ monic, with $\deg S < \deg B$ and $\deg T
< \deg A$ (if $A$ or $B$ is nonzero; if both are zero, all outputs are
zero). The underscore version requires $\mathit{lenA} \ge \mathit{lenB}
\ge 1$ with $(\mathit{lenA}, \mathit{lenB}) \ne (1, 1)$, does not make $G$
monic, needs space for $\mathit{lenB}$, $\mathit{lenB} - 1$ and
$\mathit{lenA} - 1$ coefficients in $G$, $S$, $T$, returns the length of
$G$, and does not allow aliasing (`invB` is not used). The modulus must
be prime.

::: details
The Euclidean algorithm is used for $\min(\mathit{lenA}, \mathit{lenB}) <
256$, and the half-GCD algorithm above that (`_gr_poly_xgcd_euclidean`,
`_gr_poly_xgcd_hgcd`).
:::
:::

::: function
slong _fmpz_mod_poly_xgcd_euclidean_f(fmpz_t f, fmpz * G, fmpz * S, fmpz * T, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_xgcd_euclidean_f(fmpz_t f, fmpz_mod_poly_t G, fmpz_mod_poly_t S, fmpz_mod_poly_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
slong _fmpz_mod_poly_xgcd_f(fmpz_t f, fmpz * G, fmpz * S, fmpz * T, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_xgcd_f(fmpz_t f, fmpz_mod_poly_t G, fmpz_mod_poly_t S, fmpz_mod_poly_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Extended GCD for a possibly composite modulus by the Euclidean
algorithm: either set $f = 1$ and compute the outputs as in
`fmpz_mod_poly_xgcd` (with the same conventions for the underscore
versions, where `invB` is the inverse of the leading coefficient of $B$),
or set $f$ to a nontrivial factor of $n$. The `_xgcd_f` names are aliases
of the `_xgcd_euclidean_f` functions.
:::

::: function
slong _fmpz_mod_poly_gcdinv(fmpz * G, fmpz * S, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_gcdinv(fmpz_mod_poly_t G, fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
slong _fmpz_mod_poly_gcdinv_euclidean(fmpz * G, fmpz * S, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invA, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_gcdinv_euclidean(fmpz_mod_poly_t G, fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Compute $G = \gcd(A, B)$ and $S$, both reduced modulo $B$, with $G \equiv
SA \pmod B$; when $G = 1$, $S$ is the inverse of $A$ modulo $B$. Requires
$\mathit{len}(B) \ge 2$ (an exception is raised otherwise); if $A \equiv 0
\pmod B$, sets $G = S = 0$. The underscore versions require $0 <
\mathit{lenA} < \mathit{lenB}$, write `lenA` coefficients of $G$ and
$\mathit{lenB} - 1$ of $S$, and return the length of $G$, which is not
made monic. `fmpz_mod_poly_gcdinv` uses the Euclidean algorithm for
$\mathit{lenB} < 16$ and the extended half-GCD algorithm otherwise. The
modulus must be prime.
:::

::: function
slong _fmpz_mod_poly_gcdinv_f(fmpz_t f, fmpz * G, fmpz * S, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_gcdinv_f(fmpz_t f, fmpz_mod_poly_t G, fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
slong _fmpz_mod_poly_gcdinv_euclidean_f(fmpz_t f, fmpz * G, fmpz * S, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_t invA, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_gcdinv_euclidean_f(fmpz_t f, fmpz_mod_poly_t G, fmpz_mod_poly_t S, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)

Versions of the previous functions for a possibly composite modulus:
either set $f = 1$ and compute $G$ and $S$, or set $f$ to a nontrivial
factor of $n$.
:::

::: function
int _fmpz_mod_poly_invmod(fmpz * A, const fmpz * B, slong lenB, const fmpz * P, slong lenP, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_invmod(fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_poly_t P, const fmpz_mod_ctx_t ctx)

If $B$ is invertible modulo $P$, set $A$ to the inverse (reduced modulo
$P$) and return 1; otherwise return 0, with $A$ undefined. The modulus
must be prime. `fmpz_mod_poly_invmod` raises an exception if $\deg P < 1$
and reduces $B$ modulo $P$ first. The underscore version requires $0 <
\mathit{lenB} < \mathit{lenP}$ (allowing zero-padding of $B$), writes
$\mathit{lenP} - 1$ coefficients and does not allow aliasing.
:::

::: function
int _fmpz_mod_poly_invmod_f(fmpz_t f, fmpz * A, const fmpz * B, slong lenB, const fmpz * P, slong lenP, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_invmod_f(fmpz_t f, fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_poly_t P, const fmpz_mod_ctx_t ctx)

Versions of the previous functions for a possibly composite modulus:
either set $f = 1$ and behave as `fmpz_mod_poly_invmod`, or set $f$ to a
nontrivial factor of $n$.
:::

## Resultant and discriminant {#resultant}

::: function
void _fmpz_mod_poly_resultant(fmpz_t res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_resultant(fmpz_t res, const fmpz_mod_poly_t f, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx)

Set `res` to the resultant of the two polynomials (zero if one of them is
zero), computed with `_gr_poly_resultant`, which uses the Euclidean or
the half-GCD algorithm. The underscore version requires $\mathit{len1}
\ge \mathit{len2} > 0$. The modulus should be prime; for a composite
modulus the function raises an exception if it meets a non-invertible
leading coefficient.
:::

::: function
void _fmpz_mod_poly_discriminant(fmpz_t d, const fmpz * poly, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_discriminant(fmpz_t d, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Set $d$ to the discriminant of $f$, normalised as
$$
\operatorname{disc}(f) = (-1)^{n(n-1)/2}\, \operatorname{res}(f, f') \,
\operatorname{lc}(f)^{-(n - m)}
$$
with $n = \deg f$ and $m = \deg f'$ (so $\operatorname{lc}(f)^{-1}$ when
$p \nmid n$); then $\operatorname{disc}(f) = \operatorname{lc}(f)^{2n-2}
\prod_{i < j} (r_i - r_j)^2$ over the roots $r_i$ of $f$. The
discriminant of a polynomial of degree less than 1, or with $f' = 0$, is
zero. The underscore version requires $\mathit{len} > 1$.
:::

## Derivative {#derivative}

::: function
void _fmpz_mod_poly_derivative(fmpz * res, const fmpz * poly, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_derivative(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Set `res` to the derivative of `poly`. The underscore version writes
$\mathit{len} - 1$ coefficients and allows aliasing.
:::

## Evaluation {#evaluation}

::: function
void _fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz * poly, slong len, const fmpz_t a, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_evaluate_fmpz(fmpz_t res, const fmpz_mod_poly_t poly, const fmpz_t a, const fmpz_mod_ctx_t ctx)

Set `res` to the value of the polynomial at $a$ (reduced modulo $n$ first
if necessary), by Horner's rule. `res` may alias $a$ in the
non-underscore version, but not a coefficient of `poly`.
:::

::: function
void _fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys, const fmpz * coeffs, slong len, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_evaluate_fmpz_vec(fmpz * ys, const fmpz_mod_poly_t poly, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)

Evaluate the polynomial at the $n$ reduced points `xs`, writing the values
to `ys`. Uses Horner's rule for each point if the length of the
polynomial is less than 32, and fast multipoint evaluation otherwise.
:::

::: function
void _fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys, const fmpz * coeffs, slong len, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_evaluate_fmpz_vec_iter(fmpz * ys, const fmpz_mod_poly_t poly, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)

Evaluate the polynomial at the $n$ reduced points `xs` by Horner's rule
for each point, in $O(n \cdot \mathit{len})$ operations.
:::

::: function
void _fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys, const fmpz * poly, slong plen, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_evaluate_fmpz_vec_fast(fmpz * ys, const fmpz_mod_poly_t poly, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_evaluate_fmpz_vec_fast_precomp(fmpz * vs, const fmpz * poly, slong plen, fmpz_poly_struct * const * tree, slong len, const fmpz_mod_ctx_t ctx)

Evaluate the polynomial at the reduced points by fast multipoint
evaluation: the remainders modulo the products in a subproduct tree of
the points are computed from the root to the leaves, in $O(M(n) \log n)$
operations ([](algorithms/polynomials/subproduct-tree)). The `_precomp`
version takes the tree built by `_fmpz_mod_poly_tree_build` for `len`
points.
:::

::: function
fmpz_poly_struct ** _fmpz_mod_poly_tree_alloc(slong len)
void _fmpz_mod_poly_tree_free(fmpz_poly_struct ** tree, slong len)
void _fmpz_mod_poly_tree_build(fmpz_poly_struct ** tree, const fmpz * roots, slong len, const fmpz_mod_ctx_t ctx)

Allocate, free and fill a subproduct tree for `len` points: level 0 holds
the linear factors $x - r_i$ and each following level the products of
pairs from the level below; the root (the product of all factors) is not
computed.
:::

## Products of linear factors and roots {#products}

::: function
void _fmpz_mod_poly_product_roots_fmpz_vec(fmpz * poly, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_product_roots_fmpz_vec(fmpz_mod_poly_t poly, const fmpz * xs, slong n, const fmpz_mod_ctx_t ctx)

Set `poly` to the monic polynomial $(x - x_0)(x - x_1) \cdots (x -
x_{n-1})$ for reduced $x_i$, by a quadratic algorithm for $n < 20$ and by
splitting the roots into halves otherwise. The underscore version writes
$n + 1$ coefficients and does not allow aliasing.
:::

::: function
int fmpz_mod_poly_find_distinct_nonzero_roots(fmpz * roots, const fmpz_mod_poly_t P, const fmpz_mod_ctx_t ctx)

If the nonzero polynomial $P$ has $\deg P$ distinct nonzero roots in
$\mathbb{F}_p$, writes them to `roots[0]`, ..., `roots[deg(P) - 1]` and
returns 1; otherwise returns 0. The modulus must be prime.

::: details
The function first splits $P$ by $\gcd(P, x^{(p-1)/2} \mp 1)$, which
decides whether $P$ has the required form, and then splits the factors
with `_fmpz_mod_poly_split_rabin`. See [](nmod_poly/roots) for the
word-size version and [](algorithms/polynomials/root-finding).
:::
:::

::: function
void _fmpz_mod_poly_split_rabin(fmpz_mod_poly_t a, fmpz_mod_poly_t b, const fmpz_mod_poly_t f, const fmpz_t halfp, fmpz_mod_poly_t t, fmpz_mod_poly_t t2, flint_rand_t randstate, const fmpz_mod_ctx_t ctx)

Splits $f$ of degree at least 2, assumed to have $\deg f$ distinct
nonzero roots in $\mathbb{F}_p$, into two nontrivial factors $a$ and $b$
with $\deg a \ge \deg b$, as $a = \gcd(f, (x + \delta)^{(p-1)/2} - 1)$ for
random $\delta$ (degree 2 is solved directly with a square root).
`halfp` must be $(p - 1)/2$, and $t$, $t_2$ are scratch polynomials.
:::

## Composition {#composition}

::: function
void _fmpz_mod_poly_compose(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_compose(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_ctx_t ctx)

Set `res` to the composition $\mathit{poly1}(\mathit{poly2})$, by
`_gr_poly_compose`. The underscore version requires nonzero inputs,
writes $(\mathit{len1} - 1)(\mathit{len2} - 1) + 1$ coefficients and does
not allow aliasing.
:::

## Modular composition {#compose-mod}

All functions in this section compute $f(g) \bmod h$, where $h$ is
nonzero with invertible leading coefficient and $g$ is reduced modulo
$h$; the underscore versions require $g$ zero-padded to length
$\mathit{lenh} - 1$, write $\mathit{lenh} - 1$ coefficients and do not
allow aliasing. They are wrappers of the corresponding `gr_poly`
functions; see [](algorithms/polynomials/composition#modular).

::: function
void _fmpz_mod_poly_compose_mod(fmpz * res, const fmpz * f, slong lenf, const fmpz * g, const fmpz * h, slong lenh, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_compose_mod(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx)

Set `res` to $f(g) \bmod h$ (in the non-underscore version $f, g, h$ are
`poly1`, `poly2`, `poly3`), choosing between Horner's rule and the
Brent–Kung algorithm according to the lengths.
:::

::: function
void _fmpz_mod_poly_compose_mod_horner(fmpz * res, const fmpz * f, slong lenf, const fmpz * g, const fmpz * h, slong lenh, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_compose_mod_horner(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx)

Set `res` to $f(g) \bmod h$ by Horner's rule, with $\mathit{lenf} - 1$
multiplications modulo $h$.
:::

::: function
void _fmpz_mod_poly_compose_mod_brent_kung(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, const fmpz * poly3, slong len3, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_compose_mod_brent_kung(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3, const fmpz_mod_ctx_t ctx)

Set `res` to $f(g) \bmod h$ by the Brent–Kung baby-step giant-step
algorithm [@BrentKung1978], with about $2\sqrt{n}$ multiplications
modulo $h$ and one product of a $\sqrt{n} \times \sqrt{n}$ by a $\sqrt{n}
\times n$ matrix, $n = \deg h$. Requires $\deg f < \deg h$.
:::

::: function
void _fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, const fmpz * poly3, slong len3, const fmpz * poly3inv, slong len3inv, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_compose_mod_brent_kung_preinv(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv, const fmpz_mod_ctx_t ctx)

Like `fmpz_mod_poly_compose_mod_brent_kung`, given `poly3inv`, the
inverse of the reversal of $h$ = `poly3` modulo $x^{\mathit{len3}}$.
:::

::: function
void _fmpz_mod_poly_precompute_matrix(fmpz_mat_t A, const fmpz * poly1, const fmpz * poly2, slong len2, const fmpz * poly2inv, slong len2inv, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_precompute_matrix(fmpz_mat_t A, const fmpz_mod_poly_t poly1, const fmpz_mod_poly_t poly2, const fmpz_mod_poly_t poly2inv, const fmpz_mod_ctx_t ctx)

Set the rows of $A$ to the coefficients of $g^i \bmod h$ for $0 \le i <
m$, where $g$ = `poly1` is reduced modulo $h$ = `poly2`, `poly2inv` is
the inverse of the reversal of $h$, and $A$ must have been initialised
as an $m \times \deg h$ matrix with $m = \lfloor \sqrt{\deg h} \rfloor +
1$. This is the precomputation for repeated compositions with the same
$g$ and $h$.
:::

::: function
void _fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz * res, const fmpz * poly1, slong len1, const fmpz_mat_t A, const fmpz * poly3, slong len3, const fmpz * poly3inv, slong len3inv, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_compose_mod_brent_kung_precomp_preinv(fmpz_mod_poly_t res, const fmpz_mod_poly_t poly1, const fmpz_mat_t A, const fmpz_mod_poly_t poly3, const fmpz_mod_poly_t poly3inv, const fmpz_mod_ctx_t ctx)

Set `res` to $f(g) \bmod h$ for $f$ = `poly1` with $\deg f < \deg h$,
$h$ = `poly3`, given the matrix $A$ of powers of $g$ computed by
`fmpz_mod_poly_precompute_matrix` and the inverse `poly3inv` of the
reversal of $h$. Only the matrix product and the giant steps remain,
which makes repeated compositions with the same $g$, such as the
Frobenius iterations in factorisation, cheaper.
:::

::: function
void _fmpz_mod_poly_reduce_matrix_mod_poly(fmpz_mat_t A, const fmpz_mat_t B, const fmpz_mod_poly_t f, const fmpz_mod_ctx_t ctx)

Initialises $A$ and sets its rows to the rows of $B$ reduced modulo $f$:
$B$ is a matrix of powers $g^i \bmod h$ as above and $f$ a divisor of $h$,
so that $A$ is the corresponding matrix modulo $f$, with $\lfloor
\sqrt{\deg f} \rfloor + 1$ rows and $\deg f$ columns.
:::

## Square roots {#sqrt}

::: function
int _fmpz_mod_poly_sqrt(fmpz * s, const fmpz * p, slong len, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_sqrt(fmpz_mod_poly_t b, const fmpz_mod_poly_t a, const fmpz_mod_ctx_t ctx)

If the polynomial is a perfect square, set the output to a square root
and return 1; otherwise return 0 (and set `b` to zero). The modulus must
be prime. The underscore version writes $\mathit{len}/2 + 1$
coefficients.

::: details
After removing an even power of $x$, the constant coefficient must be a
square $c^2$ in $\mathbb{F}_p$; the square root is $c$ times the power
series square root of $a/a_0$, computed to half the length, and the
result is verified by squaring. For $p = 2$, $a$ is a square if and only
if it is a polynomial in $x^2$.
:::
:::

## Minimal polynomials of sequences {#minpoly}

A *minimal generating polynomial* of a sequence $s_0, \ldots,
s_{\ell-1}$ is a monic $f = x^d + c_{d-1} x^{d-1} + \cdots + c_0$ of least
degree such that $s_{i+d} + \sum_{j<d} c_j s_{i+j} = 0$ for all $0 \le i <
\ell - d$. The modulus must be prime.

::: function
slong _fmpz_mod_poly_minpoly(fmpz * poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_minpoly(fmpz_mod_poly_t poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx)

Set `poly` to a minimal generating polynomial of `(seq, len)`, choosing
the Berlekamp–Massey algorithm for $\mathit{len} < \max(200, 530 - 22s)$,
$s$ being the number of limbs of the modulus, and the half-GCD algorithm
otherwise. The underscore version needs space for $\mathit{len} + 1$
coefficients, returns the length of `poly` and does not allow aliasing.
:::

::: function
slong _fmpz_mod_poly_minpoly_bm(fmpz * poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_minpoly_bm(fmpz_mod_poly_t poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx)

Like `fmpz_mod_poly_minpoly`, using the Berlekamp–Massey algorithm
[@Ber1968] [@Mas1969], in $O(\mathit{len} \cdot d)$ operations for a
minimal polynomial of degree $d$.
:::

::: function
slong _fmpz_mod_poly_minpoly_hgcd(fmpz * poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_minpoly_hgcd(fmpz_mod_poly_t poly, const fmpz * seq, slong len, const fmpz_mod_ctx_t ctx)

Like `fmpz_mod_poly_minpoly`, using a half-GCD of $x^{\mathit{len}}$ and
the reversed sequence, in $O(M(\mathit{len}) \log \mathit{len})$
operations regardless of the degree.
:::

### Incremental Berlekamp–Massey {#berlekamp-massey}

An `fmpz_mod_berlekamp_massey_t` processes an unbounded stream of points
$a_1, a_2, \ldots$. After $n$ points have been added,
`fmpz_mod_berlekamp_massey_reduce` computes polynomials $U$, $V$, $R$ in
the extended Euclidean remainder sequence of $x^n$ and $a_1 x^{n-1} +
\cdots + a_n$ with
$$
U x^n + V (a_1 x^{n-1} + \cdots + a_{n-1} x + a_n) = R, \qquad \deg U <
\deg V \le n/2, \qquad \deg R < n/2.
$$
Writing $v_i$ for the coefficients of $V$, this means that $\sum_i v_i
a_{j+i} = 0$ for $1 \le j < n - \deg R$: $V$ annihilates the sequence if
$\deg R < \deg V$. Unlike `fmpz_mod_poly_minpoly`, which can return a
polynomial of degree up to $n - 1$ that annihilates all the points, the
result is unique and the points can be added incrementally, as in sparse
interpolation where the number of points needed is not known in advance.

::: type
fmpz_mod_berlekamp_massey_struct
fmpz_mod_berlekamp_massey_t

The state of an incremental Berlekamp–Massey computation: the points
added so far and the polynomials of the remainder sequence.
:::

::: function
void fmpz_mod_berlekamp_massey_init(fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_berlekamp_massey_clear(fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx)

Initialise `B` with an empty stream of points, respectively free it.
:::

::: function
void fmpz_mod_berlekamp_massey_start_over(fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx)

Empties the stream of points of `B`.
:::

::: function
void fmpz_mod_berlekamp_massey_add_points(fmpz_mod_berlekamp_massey_t B, const fmpz * a, slong count, const fmpz_mod_ctx_t ctx)
void fmpz_mod_berlekamp_massey_add_zeros(fmpz_mod_berlekamp_massey_t B, slong count, const fmpz_mod_ctx_t ctx)
void fmpz_mod_berlekamp_massey_add_point(fmpz_mod_berlekamp_massey_t B, const fmpz_t a, const fmpz_mod_ctx_t ctx)
void fmpz_mod_berlekamp_massey_add_point_ui(fmpz_mod_berlekamp_massey_t B, ulong a, const fmpz_mod_ctx_t ctx)

Append points (reduced modulo $n$) to the stream of `B`: `count` points
from `a`, `count` zeros, or a single point. $V$ and $R$ are not updated
until the next call to `fmpz_mod_berlekamp_massey_reduce`.
:::

::: function
int fmpz_mod_berlekamp_massey_reduce(fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx)

Brings $V$ and $R$ up to date with the points added so far and returns 1
if $V$ changed, 0 otherwise. For example, after the points $1, 1, 2, 3$
the call returns 1 with $\deg R < \deg V = 2$, since the Fibonacci
recurrence has been recognised; after adding $5, 8$ another call returns
0.
:::

::: function
slong fmpz_mod_berlekamp_massey_point_count(const fmpz_mod_berlekamp_massey_t B)
const fmpz * fmpz_mod_berlekamp_massey_points(const fmpz_mod_berlekamp_massey_t B)

Return the number of points stored in `B` and a pointer to them (which
may be `NULL` if there are none).
:::

::: function
const fmpz_mod_poly_struct * fmpz_mod_berlekamp_massey_V_poly(const fmpz_mod_berlekamp_massey_t B)
const fmpz_mod_poly_struct * fmpz_mod_berlekamp_massey_R_poly(const fmpz_mod_berlekamp_massey_t B)

Return pointers to the polynomials $V$ and $R$ of `B`, as of the last
call to `fmpz_mod_berlekamp_massey_reduce`.
:::

::: function
void fmpz_mod_berlekamp_massey_print(const fmpz_mod_berlekamp_massey_t B, const fmpz_mod_ctx_t ctx)

Prints $V$ and the stored points of `B` to standard output, for
debugging.
:::

## Radix conversion {#radix}

Given a nonconstant polynomial $R$ with invertible leading coefficient,
every $F$ has a unique expansion
$$
F = \sum_{i=0}^{N} B_i R^i, \qquad \deg B_i < \deg R, \qquad N =
\lfloor \deg F / \deg R \rfloor.
$$
The functions below compute it by divide and conquer with divisions by
$R^{2^i}$, in $O(M(\deg F) \log \deg F)$ operations, reusing the powers of
$R$ and the power series inverses of their reversals, which are stored in
an `fmpz_mod_poly_radix_t`. The modulus need not be prime; this is used,
for example, with prime power moduli in $p$-adic computations.

::: type
fmpz_mod_poly_radix_struct
fmpz_mod_poly_radix_t

Precomputed data for radix conversion with respect to a fixed $R$: the
powers $R^{2^i}$, the inverses of their reversals, the degree of $R$, the
number $k$ of powers and the inverse of the leading coefficient of $R$.
:::

::: function
void fmpz_mod_poly_radix_init(fmpz_mod_poly_radix_t D, const fmpz_mod_poly_t R, slong degF, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_radix_clear(fmpz_mod_poly_radix_t D)

Initialise `D` with the data for radix conversion with respect to $R$ of
polynomials of degree at most `degF`, respectively free it. $R$ must have
degree at least 1 and an invertible leading coefficient.
:::

::: function
void fmpz_mod_poly_radix(fmpz_mod_poly_struct ** B, const fmpz_mod_poly_t F, const fmpz_mod_poly_radix_t D, const fmpz_mod_ctx_t ctx)

Sets `B[0]`, ..., `B[N]` to the digits $B_i$ of $F$ with respect to the
radix $R$ of `D`, where $N = \lfloor \deg F / \deg R \rfloor$ and $\deg F$
is at most the bound given to `fmpz_mod_poly_radix_init`.
:::

::: function
void _fmpz_mod_poly_radix_init(fmpz ** Rpow, fmpz ** Rinv, const fmpz * R, slong lenR, slong k, const fmpz_t invL, const fmpz_mod_ctx_t ctx)

Computes the powers $R^{2^i}$ and the inverses of their reversals modulo
$x^{2^i \deg R}$ for $0 \le i < k$, where `invL` is the inverse of the
leading coefficient of $R$. `Rpow[i]` and `Rinv[i]` need space for $2^i
\deg R + 1$ and $2^i \deg R$ coefficients. The data can be used for all
$F$ with $\mathrm{len}(F) \le 2^k \deg R$.
:::

::: function
void _fmpz_mod_poly_radix(fmpz ** B, const fmpz * F, fmpz ** Rpow, fmpz ** Rinv, slong degR, slong k, slong i, fmpz * W, const fmpz_mod_ctx_t ctx)

The recursive step of `fmpz_mod_poly_radix`: $F$ has length $2^{i+1}
\deg R$ (zero-padded), $k$ is the number of digits already accounted for
(0 at the top level), every `B[j]` has space for $\deg R$ coefficients,
and $W$ is scratch space of length $2^{i+1} \deg R$.
:::

## Input and output {#io}

A polynomial is written as its length, the modulus, and the coefficients,
for example `"4 6  1 2 0 5"` for $5x^3 + 2x + 1$ modulo 6 (two spaces
before the first coefficient), and in the pretty format as `5*x^3+2*x+1`.
The printing functions return a positive value on success and a
non-positive value on failure.

::: function
char * fmpz_mod_poly_get_str(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
char * fmpz_mod_poly_get_str_pretty(const fmpz_mod_poly_t poly, const char * x, const fmpz_mod_ctx_t ctx)

Return a string for `poly` in the format of `fmpz_poly_get_str` (the
length and the coefficients, without the modulus), respectively the
pretty format with variable name `x`. The string must be freed with
`flint_free`.
:::

::: function
int _fmpz_mod_poly_fprint(FILE * file, const fmpz * poly, slong len, const fmpz_t p)
int fmpz_mod_poly_fprint(FILE * file, const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)
int _fmpz_mod_poly_print(const fmpz * poly, slong len, const fmpz_t p)
int fmpz_mod_poly_print(const fmpz_mod_poly_t poly, const fmpz_mod_ctx_t ctx)

Print the polynomial, with its length and modulus, to `file` or to
standard output.
:::

::: function
int fmpz_mod_poly_fprint_pretty(FILE * file, const fmpz_mod_poly_t poly, const char * x, const fmpz_mod_ctx_t ctx)
int fmpz_mod_poly_print_pretty(const fmpz_mod_poly_t poly, const char * x, const fmpz_mod_ctx_t ctx)

Print the polynomial in the pretty format with variable name `x`.
:::

::: function
int fmpz_mod_poly_fread(FILE * file, fmpz_mod_poly_t poly, fmpz_mod_ctx_t ctx)

Reads a polynomial in the format of `fmpz_mod_poly_fprint` from `file`,
sets the modulus of `ctx` to the modulus read (which must be at least 2)
and `poly` to the polynomial, and returns 1 on success and 0 on failure.
:::
