---
title: fmpq_poly.h
header: fmpq_poly.h
description: Dense univariate polynomials over $\mathbb{Q}$, stored as an integer polynomial with a common denominator, with arithmetic, division, GCDs, power series and transcendental functions, composition, reversion and interpolation.
category: Univariate polynomials
order: 3
---

This module implements dense univariate polynomials over $\mathbb{Q}$. A
polynomial $f \in \mathbb{Q}[x]$ is stored as an integer polynomial (the
*numerator*, an `fmpz` coefficient array) together with a single positive
integer *denominator*, $f = (a_0 + a_1 x + \cdots + a_{n-1} x^{n-1}) / d$.
A common denominator is usually much more efficient than a vector of
rational coefficients: most operations reduce to `fmpz_poly` arithmetic
on the numerators, with a few integer gcds to keep the representation
in lowest terms.

Multiplication uses Kronecker substitution and FFTs over $\Z$, division
[pseudo-division](fmpz_poly/divrem#pseudo), GCDs the modular and
heuristic algorithms of `fmpz_poly`
([](algorithms/polynomials/gcd-integers)), and power series Newton
iteration with fraction-free basecases. The guide
[](fmpq_poly/getting-started) explains the representation, canonical
forms and the rules for the underscore functions.

## Types and conventions {#conventions}

::: type
fmpq_poly_struct
fmpq_poly_t

A polynomial over $\mathbb{Q}$. The structure has the fields `coeffs` (an array of
`alloc` integer coefficients of the numerator), `alloc`, `length` and
`den` (an `fmpz_t` holding the denominator). An `fmpq_poly_t` is an array
of one `fmpq_poly_struct`, so it is passed by reference.
:::

- A polynomial is *normalised* if its length is zero or its leading
  numerator coefficient is nonzero. It is *canonical* if moreover the
  denominator is positive and coprime to the content of the numerator;
  the zero polynomial is canonical with denominator 1. The canonical
  representation is unique.
- Functions without a leading underscore take and return canonical
  polynomials and, unless stated otherwise, allow aliasing between all
  inputs and outputs.
- Functions with a leading underscore take a polynomial as a triple
  `(poly, den, len)` of a numerator array, a denominator and a length.
  They generally expect inputs in lowest terms without zero padding and
  do not allow aliasing; each function states its conditions and whether
  the output is in canonical form.
- Some underscore functions and the functions with suffix `_can` take a
  flag `can`; with `can = 0` they skip the final reduction to lowest terms,
  which saves gcds in long chains of operations (*weak
  canonicalisation*). The result is then mathematically correct but not
  canonical; call `fmpq_poly_canonicalise` before using it with other
  functions.
## Memory management

::: function
void fmpq_poly_init(fmpq_poly_t poly)

Initialises the polynomial for use.  The length is set to zero.

:::

::: function
void fmpq_poly_init2(fmpq_poly_t poly, slong alloc)

Initialises the polynomial with space for at least `alloc`
coefficients and sets the length to zero. The `alloc` coefficients
are all set to zero.

:::

::: function
void fmpq_poly_realloc(fmpq_poly_t poly, slong alloc)

Reallocates the given polynomial to have space for `alloc`
coefficients. If `alloc` is zero then the polynomial is cleared
and then reinitialised.  If the current length is greater than
`alloc` then `poly` is first truncated to length
`alloc`. Note that this might leave the rational polynomial in
non-canonical form.

:::

::: function
void fmpq_poly_fit_length(fmpq_poly_t poly, slong len)

If `len` is greater than the number of coefficients currently
allocated, then the polynomial is reallocated to have space for at
least `len` coefficients, at least doubling the allocation so that
many calls in small increments are efficient. No data is lost.

:::

::: function
void _fmpq_poly_set_length(fmpq_poly_t poly, slong len)

Sets the length of the numerator polynomial to `len`, demoting
coefficients beyond the new length.  Note that this method does
not guarantee that the rational polynomial is in canonical form.

:::

::: function
void fmpq_poly_clear(fmpq_poly_t poly)

Clears the given polynomial, releasing any memory used. The polynomial
must be reinitialised in order to be used again.

:::

::: function
void _fmpq_poly_normalise(fmpq_poly_t poly)

Sets the length of `poly` so that the top coefficient is
non-zero. If all coefficients are zero, the length is set to zero.
Note that this function does not guarantee the coprimality of the
numerator polynomial and the integer denominator.

:::

::: function
void _fmpq_poly_canonicalise(fmpz * poly, fmpz_t den, slong len)

Puts `(poly, den)` of length `len` into canonical form.

::: details
It is assumed that the array `poly` contains a non-zero entry in
position `len - 1` whenever `len > 0`.  Assumes that `den`
is non-zero.

:::
:::

::: function
void fmpq_poly_canonicalise(fmpq_poly_t poly)

Puts the polynomial `poly` into canonical form.  Firstly, the length
is set to the actual length of the numerator polynomial.  For non-zero
polynomials, it is then ensured that the numerator and denominator are
coprime and that the denominator is positive.  The canonical form of the
zero polynomial is a zero numerator polynomial and a one denominator.

:::

::: function
int _fmpq_poly_is_canonical(const fmpz * poly, const fmpz_t den, slong len)

Returns whether the polynomial is in canonical form.

:::

::: function
int fmpq_poly_is_canonical(const fmpq_poly_t poly)

Returns whether the polynomial is in canonical form.

:::

## Polynomial parameters

::: function
slong fmpq_poly_degree(const fmpq_poly_t poly)

Returns the degree of `poly`, which is one less than its length, as
a `slong`.

:::

::: function
slong fmpq_poly_length(const fmpq_poly_t poly)

Returns the length of `poly`.

:::

## Accessing the numerator and denominator

::: macro
fmpq_poly_numref(poly)
fmpq_poly_denref(poly)

The numerator coefficient array (`fmpz *`, which may be `NULL` for a zero
polynomial with no allocated space) and the denominator (`fmpz_t`) of
`poly`, as macros returning `(poly)->coeffs` and `(poly)->den`. After
modifying them, restore the canonical form with `fmpq_poly_canonicalise`.
:::

::: function
void fmpq_poly_get_numerator(fmpz_poly_t res, const fmpq_poly_t poly)

Sets `res` to the numerator of `poly`, the integer polynomial $d \cdot
\mathit{poly}$ where $d$ is the denominator.

:::

::: function
void fmpq_poly_get_denominator(fmpz_t den, const fmpq_poly_t poly)

Sets `den` to the denominator of `poly`.

:::

## Random testing

The functions `fmpq_poly_randtest_foo` provide random
polynomials suitable for testing.  On an integer level, this
means that long strings of zeros and ones in the binary
representation are favoured as well as the special absolute
values $0$, $1$, `COEFF_MAX`, and `WORD_MAX`.  On a
polynomial level, the integer numerator has a reasonable chance
to have a non-trivial content.

::: function
void fmpq_poly_randtest(fmpq_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)

Sets $f$ to a random polynomial with coefficients up to the given
length and where each coefficient has up to the given number of bits.
The coefficients are signed randomly.  One must call
`flint_rand_init` before calling this function.

:::

::: function
void fmpq_poly_randtest_unsigned(fmpq_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)

Sets $f$ to a random polynomial with coefficients up to the given length
and where each coefficient has up to the given number of bits.  One must
call `flint_rand_init` before calling this function.

:::

::: function
void fmpq_poly_randtest_not_zero(fmpq_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)

As for `fmpq_poly_randtest` except that `len` and `bits`
may not be zero and the polynomial generated is guaranteed not to be the
zero polynomial.  One must call `flint_rand_init` before calling
this function.

:::

## Assignment, swap, negation

::: function
void fmpq_poly_set(fmpq_poly_t poly1, const fmpq_poly_t poly2)

Sets `poly1` to equal `poly2`.

:::

::: function
void fmpq_poly_set_si(fmpq_poly_t poly, slong x)

Sets `poly` to the integer $x$.

:::

::: function
void fmpq_poly_set_ui(fmpq_poly_t poly, ulong x)

Sets `poly` to the integer $x$.

:::

::: function
void fmpq_poly_set_fmpz(fmpq_poly_t poly, const fmpz_t x)

Sets `poly` to the integer $x$.

:::

::: function
void fmpq_poly_set_fmpq(fmpq_poly_t poly, const fmpq_t x)

Sets `poly` to the rational $x$, which is assumed to be
given in lowest terms, which is not checked.

:::

::: function
void fmpq_poly_set_fmpz_poly(fmpq_poly_t rop, const fmpz_poly_t op)

Sets the rational polynomial `rop` to the same value
as the integer polynomial `op`.

:::

::: function
void fmpq_poly_set_nmod_poly(fmpq_poly_t rop, const nmod_poly_t op)

Sets the coefficients of `rop` to the residues in `op`,
normalised to the interval $-m/2 \le r < m/2$ where $m$ is the modulus.

:::

::: function
void fmpq_poly_get_nmod_poly(nmod_poly_t rop, const fmpq_poly_t op)
void fmpq_poly_get_nmod_poly_den(nmod_poly_t rop, const fmpq_poly_t op, int den)

Set `rop` to the reduction of `op` modulo the modulus $m$ of `rop`:
the numerator coefficients are reduced modulo $m$ and, by
`fmpq_poly_get_nmod_poly` and by `fmpq_poly_get_nmod_poly_den` with `den`
nonzero, multiplied by the inverse of the denominator modulo $m$, which
must exist (not checked). With `den = 0` only the numerator is reduced.
:::

::: function
int _fmpq_poly_set_str(fmpz * poly, fmpz_t den, const char * str, slong len)

Sets `(poly, den)` to the polynomial specified by the
null-terminated string `str` of `len` coefficients. The input
format is a sequence of coefficients separated by one space.

::: details
The result is only guaranteed to be in lowest terms if all
coefficients in the input string are in lowest terms.

Returns $0$ if no error occurred. Otherwise, returns -1
in which case the resulting value of `(poly, den)` is undefined.
If `str` is not null-terminated, calling this method might result
in a segmentation fault.

:::
:::

::: function
int fmpq_poly_set_str(fmpq_poly_t poly, const char * str)

Sets `poly` to the polynomial specified by the null-terminated
string `str`. The input format is the same as the output format
of `fmpq_poly_get_str`: the length given as a decimal integer,
then two spaces, then the list of coefficients separated by one space.

::: details
The result is only guaranteed to be in canonical form if all
coefficients in the input string are in lowest terms.

Returns $0$ if no error occurred.  Otherwise, returns -1 in which case
the resulting value of `poly` is set to zero. If `str` is not
null-terminated, calling this method might result in a segmentation fault.

:::
:::

::: function
char * fmpq_poly_get_str(const fmpq_poly_t poly)

Returns the string representation of `poly`.

:::

::: function
char * fmpq_poly_get_str_pretty(const fmpq_poly_t poly, const char * var)

Returns the pretty representation of `poly`, using the
null-terminated string `var` not equal to `"\0"` as
the variable name.

:::

::: function
void fmpq_poly_zero(fmpq_poly_t poly)

Sets `poly` to zero.

:::

::: function
void fmpq_poly_one(fmpq_poly_t poly)

Sets `poly` to the constant polynomial $1$.

:::

::: function
void fmpq_poly_neg(fmpq_poly_t poly1, const fmpq_poly_t poly2)

Sets `poly1` to the additive inverse of `poly2`.

:::

::: function
void fmpq_poly_inv(fmpq_poly_t poly1, const fmpq_poly_t poly2)

Sets `poly1` to $1/\mathit{poly2}$, which requires `poly2` to be a
nonzero constant (the only units of $\mathbb{Q}[x]$); otherwise an exception
(`FLINT_ERROR`) is raised.
:::

::: function
void fmpq_poly_swap(fmpq_poly_t poly1, fmpq_poly_t poly2)

Efficiently swaps the polynomials `poly1` and `poly2`.

:::

::: function
void fmpq_poly_truncate(fmpq_poly_t poly, slong n)

Truncates `poly` to length $n$ (if it is longer) and brings the result
into canonical form.

:::

::: function
void fmpq_poly_set_trunc(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Sets `res` to a copy of `poly`, truncated to length `n`.

:::

::: function
void fmpq_poly_get_slice(fmpq_poly_t rop, const fmpq_poly_t op, slong i, slong j)

Sets `rop` to the sum of the terms of `op` of degree $k$ with $i \le k <
j$ (indices are clamped to $[0, \mathrm{len})$), in canonical form.
:::

::: function
void fmpq_poly_reverse(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

This function considers the polynomial `poly` to be of length $n$,
notionally truncating and zero padding if required, and reverses
the result.  Since the function normalises its result `res` may be
of length less than $n$.

:::

## Getting and setting coefficients {#coefficients}

::: function
void fmpq_poly_get_coeff_fmpz(fmpz_t x, const fmpq_poly_t poly, slong n)

Sets $x$ to the coefficient of $x^n$ in the *numerator* of `poly` (zero
if $n$ is at least the length). The coefficient of `poly` itself is $x$
divided by the denominator.
:::

::: function
void fmpq_poly_get_coeff_fmpq(fmpq_t x, const fmpq_poly_t poly, slong n)

Sets $x$ to the coefficient of $x^n$ in `poly`, in lowest terms.
:::

::: function
void fmpq_poly_set_coeff_si(fmpq_poly_t poly, slong n, slong x)
void fmpq_poly_set_coeff_ui(fmpq_poly_t poly, slong n, ulong x)
void fmpq_poly_set_coeff_fmpz(fmpq_poly_t poly, slong n, const fmpz_t x)
void fmpq_poly_set_coeff_fmpq(fmpq_poly_t poly, slong n, const fmpq_t x)

Set the coefficient of $x^n$, $n \ge 0$, in `poly` to $x$, keeping the
polynomial canonical (which may change the denominator and the other
numerator coefficients).
:::

## Comparison

::: function
int fmpq_poly_equal(const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Returns $1$ if `poly1` is equal to `poly2`,
otherwise returns $0$.

:::

::: function
int _fmpq_poly_equal_trunc(const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)

Returns $1$ if `poly1` and `poly2` notionally truncated to length
$n$ are equal, otherwise returns $0$.

:::

::: function
int fmpq_poly_equal_trunc(const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)

Returns $1$ if `poly1` and `poly2` notionally truncated to length
$n$ are equal, otherwise returns $0$.

:::

::: function
int _fmpq_poly_cmp(const fmpz * lpoly, const fmpz_t lden, const fmpz * rpoly, const fmpz_t rden, slong len)

Compares two non-zero polynomials, assuming they have the same length
`len > 0`.

::: details
The polynomials are expected to be provided in canonical form.

:::
:::

::: function
int fmpq_poly_cmp(const fmpq_poly_t left, const fmpq_poly_t right)

Compares the two polynomials `left` and `right`.

::: details
Returns $-1$, $0$, or $1$ as `left` is less than, equal to, or greater
than `right`.  The comparison is first done by the degree, and
then, in case of a tie, by the individual coefficients from highest
to lowest.

:::
:::

::: function
int fmpq_poly_is_one(const fmpq_poly_t poly)

Returns $1$ if `poly` is the constant polynomial $1$, otherwise
returns $0$.

:::

::: function
int fmpq_poly_is_zero(const fmpq_poly_t poly)

Returns $1$ if `poly` is the zero polynomial, otherwise returns $0$.

:::

::: function
int fmpq_poly_is_gen(const fmpq_poly_t poly)

Returns $1$ if `poly` is the degree $1$ polynomial $x$, otherwise returns
$0$.

:::

## Addition and subtraction

::: function
void _fmpq_poly_add(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)

Forms the sum `(rpoly, rden)` of `(poly1, den1, len1)` and
`(poly2, den2, len2)`, placing the result into canonical form.

::: details
Assumes that `rpoly` is an array of length the maximum of
`len1` and `len2`.  The input operands are assumed to
be in canonical form and are also allowed to be of length $0$.

`(rpoly, rden)` and `(poly1, den1)` may be aliased,
but `(rpoly, rden)` and `(poly2, den2)` may *not*
be aliased.

:::
:::

::: function
void _fmpq_poly_add_can(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, int can)

As per `_fmpq_poly_add` except that one can specify whether to
canonicalise the output or not, for use with weak canonicalisation to
prevent explosion in memory usage.

:::

::: function
void fmpq_poly_add(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Sets `res` to the sum of `poly1` and `poly2`, using
Henrici's algorithm.

:::

::: function
void fmpq_poly_add_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, int can)

As per `fmpq_poly_add` except that one can specify whether to
canonicalise the output or not, for use with weak canonicalisation to
prevent explosion in memory usage.

:::

::: function
void _fmpq_poly_add_series(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)

As per `_fmpq_poly_add` but the inputs are first notionally truncated
to length $n$. If $n$ is less than `len1` or `len2` then the
output only needs space for $n$ coefficients. We require $n \geq 0$.

:::

::: function
void _fmpq_poly_add_series_can(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n, int can)

As per `_fmpq_poly_add_can` but the inputs are first notionally
truncated to length $n$. If $n$ is less than `len1` or `len2`
then the output only needs space for $n$ coefficients. We require
$n \geq 0$.

:::

::: function
void fmpq_poly_add_series(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)

As per `fmpq_poly_add` but the inputs are first notionally
truncated to length $n$.

:::

::: function
void fmpq_poly_add_series_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n, int can)

As per `fmpq_poly_add_can` but the inputs are first notionally
truncated to length $n$.

:::

::: function
void _fmpq_poly_sub(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)

Forms the difference `(rpoly, rden)` of `(poly1, den1, len1)`
and `(poly2, den2, len2)`, placing the result into canonical form.

::: details
Assumes that `rpoly` is an array of length the maximum of
`len1` and `len2`.  The input operands are assumed to be in
canonical form and are also allowed to be of length $0$.

`(rpoly, rden)` and `(poly1, den1, len1)` may be aliased,
but `(rpoly, rden)` and `(poly2, den2, len2)` may *not* be
aliased.

:::
:::

::: function
void _fmpq_poly_sub_can(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, int can)

As per `_fmpq_poly_sub` except that one can specify whether to
canonicalise the output or not, for use with weak canonicalisation to
prevent explosion in memory usage.

:::

::: function
void fmpq_poly_sub(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Sets `res` to the difference of `poly1` and `poly2`,
using Henrici's algorithm.

:::

::: function
void fmpq_poly_sub_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, int can)

As per `_fmpq_poly_sub` except that one can specify whether to
canonicalise the output or not, for use with weak canonicalisation to
prevent explosion in memory usage.

:::

::: function
void _fmpq_poly_sub_series(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)

As per `_fmpq_poly_sub` but the inputs are first notionally truncated
to length $n$. If $n$ is less than `len1` or `len2` then the
output only needs space for $n$ coefficients. We require $n \geq 0$.

:::

::: function
void _fmpq_poly_sub_series_can(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n, int can)

As per `_fmpq_poly_sub_can` but the inputs are first notionally
truncated to length $n$. If $n$ is less than `len1` or `len2`
then the output only needs space for $n$ coefficients. We require
$n \geq 0$.

:::

::: function
void fmpq_poly_sub_series(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)

As per `fmpq_poly_sub` but the inputs are first notionally
truncated to length $n$.

:::

::: function
void fmpq_poly_sub_series_can(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n, int can)

As per `fmpq_poly_sub_can` but the inputs are first notionally
truncated to length $n$.

:::

::: function
void fmpq_poly_add_si(fmpq_poly_t res, const fmpq_poly_t poly, slong c)
void fmpq_poly_add_fmpz(fmpq_poly_t res, const fmpq_poly_t poly, const fmpz_t c)
void fmpq_poly_add_fmpq(fmpq_poly_t res, const fmpq_poly_t poly, const fmpq_t c)
void fmpq_poly_sub_si(fmpq_poly_t res, const fmpq_poly_t poly, slong c)
void fmpq_poly_sub_fmpz(fmpq_poly_t res, const fmpq_poly_t poly, const fmpz_t c)
void fmpq_poly_sub_fmpq(fmpq_poly_t res, const fmpq_poly_t poly, const fmpq_t c)
void fmpq_poly_si_sub(fmpq_poly_t res, slong c, const fmpq_poly_t poly)
void fmpq_poly_fmpz_sub(fmpq_poly_t res, const fmpz_t c, const fmpq_poly_t poly)
void fmpq_poly_fmpq_sub(fmpq_poly_t res, const fmpq_t c, const fmpq_poly_t poly)

Set `res` to $\mathit{poly} + c$, $\mathit{poly} - c$ or $c -
\mathit{poly}$ for a constant $c$.
:::

## Scalar multiplication and division

::: function
void _fmpq_poly_scalar_mul_si(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, slong c)

Sets `(rpoly, rden, len)` to the product of $c$ of
`(poly, den, len)`.

::: details
If the input is normalised, then so is the output, provided it is
non-zero.  If the input is in lowest terms, then so is the output.
However, even if neither of these conditions are met, the result
will be (mathematically) correct.

Supports exact aliasing between `(rpoly, den)`
and `(poly, den)`.

:::
:::

::: function
void _fmpq_poly_scalar_mul_ui(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, ulong c)

Sets `(rpoly, rden, len)` to the product of $c$ of
`(poly, den, len)`.

::: details
If the input is normalised, then so is the output, provided it is
non-zero.  If the input is in lowest terms, then so is the output.
However, even if neither of these conditions are met, the result
will be (mathematically) correct.

Supports exact aliasing between `(rpoly, den)`
and `(poly, den)`.

:::
:::

::: function
void _fmpq_poly_scalar_mul_fmpz(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t c)

Sets `(rpoly, rden, len)` to the product of $c$ of
`(poly, den, len)`.

::: details
If the input is normalised, then so is the output, provided it is
non-zero.  If the input is in lowest terms, then so is the output.
However, even if neither of these conditions are met, the result
will be (mathematically) correct.

Supports exact aliasing between `(rpoly, den)`
and `(poly, den)`.

:::
:::

::: function
void _fmpq_poly_scalar_mul_fmpq(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t r, const fmpz_t s)

Sets `(rpoly, rden)` to the product of $r/s$ and
`(poly, den, len)`, in lowest terms.

::: details
Assumes that `(poly, den, len)` and $r/s$ are provided in lowest
terms.  Assumes that `rpoly` is an array of length `len`.
Supports aliasing of `(rpoly, den)` and `(poly, den)`.
The `fmpz_t`'s $r$ and $s$ may not be part of `(rpoly, rden)`.

:::
:::

::: function
void fmpq_poly_scalar_mul_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)
void fmpq_poly_scalar_mul_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c)
void fmpq_poly_scalar_mul_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)

Sets `rop` to $c$ times `op`.

:::

::: function
void fmpq_poly_scalar_mul_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)

Sets `rop` to $c$ times `op`.  Assumes that the `fmpz_t c`
is not part of `rop`.

:::

::: function
void _fmpq_poly_scalar_div_fmpz(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t c)

Sets `(rpoly, rden, len)` to `(poly, den, len)` divided by $c$,
in lowest terms.

::: details
Assumes that `len` is positive.  Assumes that $c$ is non-zero.
Supports aliasing between `(rpoly, rden)` and `(poly, den)`.
Assumes that $c$ is not part of `(rpoly, rden)`.

:::
:::

::: function
void _fmpq_poly_scalar_div_si(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, slong c)

Sets `(rpoly, rden, len)` to `(poly, den, len)` divided by $c$,
in lowest terms.

::: details
Assumes that `len` is positive.  Assumes that $c$ is non-zero.
Supports aliasing between `(rpoly, rden)` and `(poly, den)`.

:::
:::

::: function
void _fmpq_poly_scalar_div_ui(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, ulong c)

Sets `(rpoly, rden, len)` to `(poly, den, len)` divided by $c$,
in lowest terms.

::: details
Assumes that `len` is positive.  Assumes that $c$ is non-zero.
Supports aliasing between `(rpoly, rden)` and `(poly, den)`.

:::
:::

::: function
void _fmpq_poly_scalar_div_fmpq(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t r, const fmpz_t s)

Sets `(rpoly, rden, len)` to `(poly, den, len)` divided by $r/s$,
in lowest terms.

::: details
Assumes that `len` is positive.  Assumes that $r/s$ is non-zero and
in lowest terms.  Supports aliasing between `(rpoly, rden)` and
`(poly, den)`. The `fmpz_t`'s $r$ and $s$ may not be part of
`(rpoly, poly)`.

:::
:::

::: function
void fmpq_poly_scalar_div_si(fmpq_poly_t rop, const fmpq_poly_t op, slong c)
void fmpq_poly_scalar_div_ui(fmpq_poly_t rop, const fmpq_poly_t op, ulong c)
void fmpq_poly_scalar_div_fmpz(fmpq_poly_t rop, const fmpq_poly_t op, const fmpz_t c)
void fmpq_poly_scalar_div_fmpq(fmpq_poly_t rop, const fmpq_poly_t op, const fmpq_t c)

Sets `rop` to `op` divided by the scalar `c`.

:::

## Multiplication

::: function
void _fmpq_poly_mul(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)

Sets `(rpoly, rden, len1 + len2 - 1)` to the product of
`(poly1, den1, len1)` and `(poly2, den2, len2)`. If the
input is provided in canonical form, then so is the output.

::: details
Assumes `len1, len2 > 0`; the operands may be given in either
order.  Allows zero-padding in the input.
Does not allow aliasing between the inputs and outputs.

:::
:::

::: function
void fmpq_poly_mul(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Sets `res` to the product of `poly1` and `poly2`.

:::

::: function
void _fmpq_poly_mullow(fmpz * rpoly, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)

Sets `(rpoly, rden, n)` to the low $n$ coefficients of
`(poly1, den1)` and `(poly2, den2)`.  The output is
not guaranteed to be in canonical form.

::: details
Assumes `len1, len2 > 0` (the operands may be given in either
order) and `0 < n <= len1 + len2 - 1`.
Allows for zero-padding in the inputs.  Does not allow aliasing between
the inputs and outputs.

:::
:::

::: function
void fmpq_poly_mullow(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)

Sets `res` to the product of `poly1` and `poly2`,
truncated to length $n$.

:::

::: function
void fmpq_poly_addmul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2)

Adds the product of `op1` and `op2` to `rop`.

:::

::: function
void fmpq_poly_submul(fmpq_poly_t rop, const fmpq_poly_t op1, const fmpq_poly_t op2)

Subtracts the product of `op1` and `op2` from `rop`.

:::

## Powering

::: function
void _fmpq_poly_pow(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, ulong e)

Sets `(rpoly, rden)` to `(poly, den)^e`, assuming
`e, len > 0`.  Assumes that `rpoly` is an array of
length at least `e * (len - 1) + 1`.  Supports aliasing
of `(rpoly, den)` and `(poly, den)`.

:::

::: function
void fmpq_poly_pow(fmpq_poly_t res, const fmpq_poly_t poly, ulong e)

Sets `res` to `poly^e`, where the only special case $0^0$ is
defined as $1$.

:::

::: function
void _fmpq_poly_pow_trunc(fmpz * res, fmpz_t rden, const fmpz * f, const fmpz_t fden, slong flen, ulong exp, slong len)

Sets `(rpoly, rden, len)` to `(poly, den)^e` truncated to length `len`,
where `len` is at most `e * (flen - 1) + 1`.

:::

::: function
void fmpq_poly_pow_trunc(fmpq_poly_t res, const fmpq_poly_t poly, ulong e, slong n)

Sets `res` to `poly^e` truncated to length `n`.

:::

## Shifting

::: function
void fmpq_poly_shift_left(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Set `res` to `poly` shifted left by $n$ coefficients. Zero
coefficients are inserted.

:::

::: function
void fmpq_poly_shift_right(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Set `res` to `poly` shifted right by $n$ coefficients.
If $n$ is equal to or greater than the current length of `poly`,
`res` is set to the zero polynomial.

:::

## Euclidean division

::: function
void _fmpq_poly_divrem(fmpz * Q, fmpz_t q, fmpz * R, fmpz_t r, const fmpz * A, const fmpz_t a, slong lenA, const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)

Finds the quotient `(Q, q)` and remainder `(R, r)` of the
Euclidean division of `(A, a)` by `(B, b)`.

::: details
Assumes that `lenA >= lenB > 0`.  Assumes that $R$ has space for
`lenA` coefficients, although only the bottom `lenB - 1` will
carry meaningful data on exit.  Supports no aliasing between the two
outputs, or between the inputs and the outputs.

An optional precomputed inverse of the leading coefficient of $B$ from
`fmpz_preinvn_init` can be supplied. Otherwise `inv` should be
`NULL`.

Note: `fmpz.h` has to be included before `fmpq_poly.h` in order for the
latter to declare this function.

:::
:::

::: function
void fmpq_poly_divrem(fmpq_poly_t Q, fmpq_poly_t R, const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Finds the quotient $Q$ and remainder $R$ of the Euclidean division of
`poly1` by `poly2`, so that $\mathit{poly1} = Q \cdot \mathit{poly2} + R$
with $\deg R < \deg \mathit{poly2}$. $Q$ and $R$ must be distinct objects
(an exception is raised otherwise), and a zero `poly2` raises
`FLINT_DIVZERO`.

::: details
The division is done over $\Z$ by pseudo-division
(`_fmpz_poly_pseudo_divrem`, see
[division and pseudo-division](fmpz_poly/divrem#pseudo)),
$\mathrm{lc}(B)^k A = QB + R$, after which the denominators are
adjusted; this avoids rational arithmetic in the inner loop.
:::
:::

::: function
void _fmpq_poly_div(fmpz * Q, fmpz_t q, const fmpz * A, const fmpz_t a, slong lenA, const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)

Finds the quotient `(Q, q)` of the Euclidean division
of `(A, a)` by `(B, b)`.

::: details
Assumes that `lenA >= lenB > 0`.  Supports no aliasing
between the inputs and the outputs.

An optional precomputed inverse of the leading coefficient of $B$ from
`fmpz_preinvn_init` can be supplied. Otherwise `inv` should be
`NULL`.

Note: `fmpz.h` has to be included before `fmpq_poly.h` in order for the
latter to declare this function.

:::
:::

::: function
void fmpq_poly_div(fmpq_poly_t Q, const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Sets $Q$ to the quotient of the Euclidean division of `poly1` by
`poly2`. A zero `poly2` raises `FLINT_DIVZERO`.

:::

::: function
void _fmpq_poly_rem(fmpz * R, fmpz_t r, const fmpz * A, const fmpz_t a, slong lenA, const fmpz * B, const fmpz_t b, slong lenB, const fmpz_preinvn_t inv)

Finds the remainder `(R, r)` of the Euclidean division
of `(A, a)` by `(B, b)`.

::: details
Assumes that `lenA >= lenB > 0`.  Supports no aliasing between
the inputs and the outputs.

An optional precomputed inverse of the leading coefficient of $B$ from
`fmpz_preinvn_init` can be supplied. Otherwise `inv` should be
`NULL`.

Note: `fmpz.h` has to be included before `fmpq_poly.h` in order for the
latter to declare this function.

:::
:::

::: function
void fmpq_poly_rem(fmpq_poly_t R, const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Finds the remainder $R$ of the Euclidean division of `poly1` by
`poly2`. A zero `poly2` raises `FLINT_DIVZERO`.

:::

## Division with precomputed powers {#powers-precomp}

These functions reduce polynomials of length up to $2n - 1$ modulo a
fixed $B$ of length $n$ by a linear combination of the precomputed
remainders $x^i \bmod B$, which is faster than repeated division when
many polynomials are reduced modulo the same $B$, for example in number
field arithmetic.

::: type
fmpq_poly_powers_precomp_struct
fmpq_poly_powers_precomp_t

The precomputed powers $x^i \bmod B$, $0 \le i < 2n - 1$: fields `powers`
(an array of `fmpq_poly_struct`) and `len` ($n$).
:::

::: function
fmpq_poly_struct * _fmpq_poly_powers_precompute(const fmpz * B, const fmpz_t denB, slong len)

Computes `2*len - 1` powers of $x$ modulo the polynomial $B$ of
the given length. This is used as a kind of precomputed inverse in
the remainder routine below.

:::

::: function
void fmpq_poly_powers_precompute(fmpq_poly_powers_precomp_t pinv, fmpq_poly_t poly)

Computes the powers $x^i \bmod B$, $0 \le i < 2n - 1$, for the polynomial
$B$ = `poly` of length $n \ge 1$ (a zero `poly` raises `FLINT_DIVZERO`).

:::

::: function
void _fmpq_poly_powers_clear(fmpq_poly_struct * powers, slong len)

Clean up resources used by precomputed powers which have been computed
by `_fmpq_poly_powers_precompute`.

:::

::: function
void fmpq_poly_powers_clear(fmpq_poly_powers_precomp_t pinv)

Clean up resources used by precomputed powers which have been computed
by `fmpq_poly_powers_precompute`.

:::

::: function
void _fmpq_poly_rem_powers_precomp(fmpz * A, fmpz_t denA, slong m, const fmpz * B, const fmpz_t denB, slong n, fmpq_poly_struct * const powers)

Set $A$ to the remainder of $A$ divide $B$ given precomputed powers mod $B$
provided by `_fmpq_poly_powers_precompute`. No aliasing is allowed.

::: details
This function is only faster if $m \leq 2\cdot n - 1$.

The output of this function is *not* canonicalised.

:::
:::

::: function
void fmpq_poly_rem_powers_precomp(fmpq_poly_t R, const fmpq_poly_t A, const fmpq_poly_t B, const fmpq_poly_powers_precomp_t B_inv)

Set $R$ to the remainder of $A$ divide $B$ given precomputed powers mod $B$
provided by `fmpq_poly_powers_precompute`.

::: details
This function is only faster if `A->length <= 2*B->length - 1`.

The output of this function is *not* canonicalised.

:::
:::

## Divisibility testing

::: function
int _fmpq_poly_divides(fmpz * qpoly, fmpz_t qden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)

Return $1$ if `(poly2, den2, len2)` divides `(poly1, den1, len1)` and
set `(qpoly, qden, len1 - len2 + 1)` to the quotient. Otherwise return
$0$. Requires that `qpoly` has space for `len1 - len2 + 1`
coefficients and that `len1 >= len2 > 0`.

:::

::: function
int fmpq_poly_divides(fmpq_poly_t q, const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Return $1$ if `poly2` divides `poly1` and set `q` to the quotient.
Otherwise return $0$.

:::

::: function
slong fmpq_poly_remove(fmpq_poly_t q, const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Sets `q` to the quotient of `poly1` by the highest power of `poly2`
which divides it, and returns the power. The divisor `poly2` must not be
constant or an exception is raised.

:::

## Power series division

::: function
void _fmpq_poly_inv_series_newton(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, slong n)

Computes the first $n$ terms of the inverse power series of
`(poly, den, len)` using Newton iteration.

::: details
The result is produced in canonical form.

Assumes that $n \geq 1$ and that `poly` has non-zero constant term;
this is not checked.
Does not support aliasing.

:::
:::

::: function
void fmpq_poly_inv_series_newton(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Computes the first $n$ terms of the inverse power series of `poly`
using Newton iteration, for $n \ge 1$. A zero constant term raises
`FLINT_DIVZERO`.

::: details
If the input has integer coefficients and constant term $\pm 1$, the
inverse is computed over $\Z$ with `_fmpz_poly_inv_series`; otherwise a
fraction-free basecase is lifted by Newton iteration
([](algorithms/polynomials/newton#inverse)).
:::

:::

::: function
void _fmpq_poly_inv_series(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong den_len, slong n)

Computes the first $n$ terms of the inverse power series of
`(poly, den, len)`.

::: details
The result is produced in canonical form.

Assumes that $n \geq 1$ and that `poly` has non-zero constant term;
this is not checked.
Does not support aliasing.

:::
:::

::: function
void fmpq_poly_inv_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Computes the first $n$ terms of the inverse power series of `poly`, for
$n \ge 1$. This is currently the same as `fmpq_poly_inv_series_newton`;
a zero constant term raises `FLINT_DIVZERO`.

:::

::: function
void _fmpq_poly_div_series(fmpz * Q, fmpz_t denQ, const fmpz * A, const fmpz_t denA, slong lenA, const fmpz * B, const fmpz_t denB, slong lenB, slong n)

Divides `(A, denA, lenA)` by `(B, denB, lenB)` as power series
over $\mathbb{Q}$, assuming $B$ has non-zero constant term and that
all lengths are positive. The constant term is not checked.

::: details
Aliasing is not supported.

This function ensures that the numerator and denominator
are coprime on exit.

:::
:::

::: function
void fmpq_poly_div_series(fmpq_poly_t Q, const fmpq_poly_t A, const fmpq_poly_t B, slong n)

Performs power series division in $\mathbb{Q}[[x]] / (x^n)$.  The function
considers the polynomials $A$ and $B$ as power series of length $n$
starting with the constant terms. Requires $n \ge 1$; a zero $B$ raises
`FLINT_DIVZERO`, and $B$ must have a nonzero constant term (not
checked).

:::

## Greatest common divisor

::: function
void _fmpq_poly_gcd(fmpz * G, fmpz_t denG, const fmpz * A, slong lenA, const fmpz * B, slong lenB)

Computes the monic greatest common divisor $G$ of $A$ and $B$.

::: details
Assumes that $G$ has space for $\operatorname{len}(B)$ coefficients,
where $\operatorname{len}(A) \geq \operatorname{len}(B) > 0$.

Aliasing between the output and input arguments is not supported.

Does not support zero-padding.

:::
:::

::: function
void fmpq_poly_gcd(fmpq_poly_t G, const fmpq_poly_t A, const fmpq_poly_t B)

Computes the monic greatest common divisor $G$ of $A$ and $B$.

::: details
In the special case when $A = B = 0$, sets $G = 0$.

:::
:::

::: function
void _fmpq_poly_xgcd(fmpz * G, fmpz_t denG, fmpz * S, fmpz_t denS, fmpz * T, fmpz_t denT, const fmpz * A, const fmpz_t denA, slong lenA, const fmpz * B, const fmpz_t denB, slong lenB)

Computes polynomials $G$, $S$, and $T$ such that
$G = \gcd(A, B) = S A + T B$, where $G$ is the monic
greatest common divisor of $A$ and $B$.

::: details
Assumes that $G$, $S$, and $T$ have space for $\operatorname{len}(B)$,
$\operatorname{len}(B)$, and $\operatorname{len}(A)$ coefficients, respectively,
where it is also assumed that $\operatorname{len}(A) \geq \operatorname{len}(B) > 0$.

Does not support zero padding of the input arguments.

:::
:::

::: function
void fmpq_poly_xgcd(fmpq_poly_t G, fmpq_poly_t S, fmpq_poly_t T, const fmpq_poly_t A, const fmpq_poly_t B)

Computes polynomials $G$, $S$, and $T$ such that
$G = \gcd(A, B) = S A + T B$, where $G$ is the monic
greatest common divisor of $A$ and $B$.

::: details
The outputs $G$, $S$, $T$ must be distinct objects (an exception is raised
otherwise). Corner cases are handled as follows. If $A = B = 0$, returns
$G = S = T = 0$.  If $A \neq 0$, $B = 0$, returns the suitable
scalar multiple of $G = A$, $S = 1$, and $T = 0$.  The case
when $A = 0$, $B \neq 0$ is handled similarly.

:::
:::

::: function
void _fmpq_poly_lcm(fmpz * L, fmpz_t denL, const fmpz * A, slong lenA, const fmpz * B, slong lenB)

Computes the monic least common multiple $L$ of $A$ and $B$.

::: details
Assumes that $L$ has space for $\operatorname{len}(A) + \operatorname{len}(B) - 1$ coefficients,
where $\operatorname{len}(A) \geq \operatorname{len}(B) > 0$.

Aliasing between the output and input arguments is not supported.

Does not support zero-padding.

:::
:::

::: function
void fmpq_poly_lcm(fmpq_poly_t L, const fmpq_poly_t A, const fmpq_poly_t B)

Computes the monic least common multiple $L$ of $A$ and $B$.

::: details
In the special case when $A = B = 0$, sets $L = 0$.

:::
:::

::: function
void _fmpq_poly_resultant(fmpz_t rnum, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)

Sets `(rnum, rden)` to the resultant of the two input
polynomials.

::: details
Assumes that `len1 >= len2 > 0`.  Does not support zero-padding
of the input polynomials.  Does not support aliasing of the input and
output arguments.

:::
:::

::: function
void fmpq_poly_resultant(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g)

Returns the resultant of $f$ and $g$.

::: details
Enumerating the roots of $f$ and $g$ over $\bar{\mathbf{Q}}$ as
$r_1, \dotsc, r_m$ and $s_1, \dotsc, s_n$, respectively, and
letting $x$ and $y$ denote the leading coefficients, the resultant
is defined as

$$
x^{n} y^{m} \prod_{i=1}^{m} \prod_{j=1}^{n} (r_i - s_j).
$$

We handle special cases as follows:  if one of the polynomials is zero,
the resultant is zero.  Otherwise, if one of the polynomials is
constant, the last term in the above expression is the empty product.

:::
:::

::: function
void _fmpq_poly_resultant_div(fmpz_t rnum, fmpz_t rden, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, const fmpz_t divisor, slong nbits)
void fmpq_poly_resultant_div(fmpq_t r, const fmpq_poly_t f, const fmpq_poly_t g, const fmpz_t divisor, slong nbits)

Set the output to the resultant of $f$ and $g$ divided by `divisor`,
assuming that this quotient is an integer of at most `nbits` bits. The
resultant of the numerators is then computed by a multimodular algorithm
that needs only enough primes for `nbits` bits, which is much faster
than `fmpq_poly_resultant` when a good bound is known (for example when
computing norms in number fields). A zero `divisor` gives zero. The
underscore version requires $\mathit{len1} \ge \mathit{len2} > 0$.
:::

## Discriminant

::: function
void _fmpq_poly_discriminant(fmpq_t res, const fmpz * poly, const fmpz_t den, slong len)

Set `res` to the discriminant of `(poly, den, len)`. Assumes
`len > 1`.

:::

::: function
void fmpq_poly_discriminant(fmpq_t res, const fmpq_poly_t poly)

Set `res` to the discriminant of `poly`, normalised so that
$\operatorname{disc}(f) = (-1)^{n(n-1)/2} \operatorname{res}(f,
f')/\operatorname{lc}(f)$, and thus $\operatorname{disc}(f) =
\operatorname{lc}(f)^{2n - 2} \prod_{i < j} (r_i - r_j)^2$, where $n$ is the
degree of $f$ and the $r_i$ are its roots. It is computed as the
discriminant of the numerator (`_fmpz_poly_discriminant`) divided by
$d^{2n-2}$, $d$ being the denominator. Polynomials of degree less than 1
have discriminant zero.

:::

## Derivative and integral

::: function
void _fmpq_poly_derivative(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len)

Sets `(rpoly, rden, len - 1)` to the derivative of
`(poly, den, len)`.  Does nothing if `len <= 1`.
Supports aliasing between the two polynomials.

:::

::: function
void fmpq_poly_derivative(fmpq_poly_t res, const fmpq_poly_t poly)

Sets `res` to the derivative of `poly`.

:::

::: function
void _fmpq_poly_nth_derivative(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, ulong n, slong len)

Sets `(rpoly, rden, len - n)` to the nth derivative of
`(poly, den, len)`.  Does nothing if `len <= n`.
Supports aliasing between the two polynomials.

:::

::: function
void fmpq_poly_nth_derivative(fmpq_poly_t res, const fmpq_poly_t poly, ulong n)

Sets `res` to the nth derivative of `poly`.

:::

::: function
void _fmpq_poly_integral(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len)

Sets `(rpoly, rden, len)` to the integral of
`(poly, den, len - 1)`.  Assumes `len >= 0`.
Supports aliasing between the two polynomials.
The output will be in canonical form if the input is
in canonical form.

:::

::: function
void fmpq_poly_integral(fmpq_poly_t res, const fmpq_poly_t poly)

Sets `res` to the integral of `poly`. The constant
term is set to zero.

:::

::: function
void _fmpq_poly_integral_offset(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, slong m)

Sets `(rpoly, rden, len)` to the offset integral of
`(poly, den, len1)` mapping $c_k$ to $c_k / (k + m)$.
Assumes `len >= 0`.
Supports aliasing between the two polynomials.
The output will be in canonical form if the input is
in canonical form.

:::

## Square roots

::: function
void _fmpq_poly_sqrt_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the
square root of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 1.
Does not support aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_sqrt_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the square root of `f`
to order `n > 1`. Requires `f` to have constant term 1.

:::

::: function
void _fmpq_poly_invsqrt_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the inverse
square root of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 1.
Does not support aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_invsqrt_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the inverse square root of
`f` to order `n > 0`. Requires `f` to have constant term 1.

:::

## Power sums

::: function
void _fmpq_poly_power_sums(fmpz * res, fmpz_t rden, const fmpz * poly, slong len, slong n)

Compute the (truncated) power sums series of the polynomial
`(poly,len)` up to length $n$ using Newton identities.

:::

::: function
void fmpq_poly_power_sums(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Compute the (truncated) power sum series of the monic polynomial
`poly` up to length $n$ using Newton identities. That is the power
series whose coefficient of degree $i$ is the sum of the $i$-th power of
all (complex) roots of the polynomial `poly`.

:::

::: function
void _fmpq_poly_power_sums_to_poly(fmpz * res, const fmpz * poly, const fmpz_t den, slong len)

Compute an integer polynomial given by its power sums series `(poly,den,len)`.

:::

::: function
void fmpq_poly_power_sums_to_fmpz_poly(fmpz_poly_t res, const fmpq_poly_t Q)

Compute the integer polynomial with content one and positive leading
coefficient given by its power sums series `Q`.

:::

::: function
void fmpq_poly_power_sums_to_poly(fmpq_poly_t res, const fmpq_poly_t Q)

Compute the monic polynomial from its power sums series `Q`.

:::

## Transcendental functions

::: function
void _fmpq_poly_log_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the
logarithm of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 1.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_log_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the logarithm of `f`
to order `n > 0`. Requires `f` to have constant term 1.

:::

::: function
void _fmpq_poly_exp_series(fmpz * g, fmpz_t gden, const fmpz * h, const fmpz_t hden, slong hlen, slong n)

Sets `(g, gden, n)` to the series expansion of the
exponential function of `(h, hden, hlen)`.  Assumes
`n > 0, hlen > 0` and
that `(h, hden, hlen)` has constant term 0.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_exp_series(fmpq_poly_t res, const fmpq_poly_t h, slong n)

Sets `res` to the series expansion of the exponential function
of `h` to order `n > 0`. Requires `h` to have
constant term 0.

:::

::: function
void _fmpq_poly_exp_expinv_series(fmpz * res1, fmpz_t res1den, fmpz * res2, fmpz_t res2den, const fmpz * h, const fmpz_t hden, slong hlen, slong n)

The same as `fmpq_poly_exp_series`, but simultaneously computes
the exponential (in `res1`, `res1den`) and its multiplicative inverse
(in `res2`, `res2den`).
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_exp_expinv_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t h, slong n)

The same as `fmpq_poly_exp_series`, but simultaneously computes
the exponential (in `res1`) and its multiplicative inverse
(in `res2`).

:::

::: function
void _fmpq_poly_atan_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the
inverse tangent of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_atan_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the inverse tangent of `f`
to order `n > 0`. Requires `f` to have constant term 0.

:::

::: function
void _fmpq_poly_atanh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the inverse
hyperbolic tangent of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_atanh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the inverse hyperbolic
tangent of `f` to order `n > 0`. Requires `f` to have
constant term 0.

:::

::: function
void _fmpq_poly_asin_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the
inverse sine of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_asin_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the inverse sine of `f`
to order `n > 0`. Requires `f` to have constant term 0.

:::

::: function
void _fmpq_poly_asinh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the inverse
hyperbolic sine of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_asinh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the inverse hyperbolic
sine of `f` to order `n > 0`. Requires `f` to have
constant term 0.

:::

::: function
void _fmpq_poly_tan_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the
tangent function of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Does not support aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_tan_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the tangent function
of `f` to order `n > 0`. Requires `f` to have
constant term 0.

:::

::: function
void _fmpq_poly_sin_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the
sine of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_sin_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the sine of `f`
to order `n > 0`. Requires `f` to have constant term 0.

:::

::: function
void _fmpq_poly_cos_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the
cosine of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_cos_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the cosine of `f`
to order `n > 0`. Requires `f` to have constant term 0.

:::

::: function
void _fmpq_poly_sin_cos_series(fmpz * s, fmpz_t sden, fmpz * c, fmpz_t cden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(s, sden, n)` to the series expansion of the
sine of `(f, fden, flen)`, and `(c, cden, n)` to the series
expansion of the cosine.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_sin_cos_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t f, slong n)

Sets `res1` to the series expansion of the sine of `f`
to order `n > 0`, and `res2` to the series expansion
of the cosine. Requires `f` to have constant term 0.

:::

::: function
void _fmpq_poly_sinh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the
hyperbolic sine of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Does not support aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_sinh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the hyperbolic sine of `f`
to order `n > 0`. Requires `f` to have constant term 0.

:::

::: function
void _fmpq_poly_cosh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the hyperbolic
cosine of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Does not support aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_cosh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the hyperbolic cosine of
`f` to order `n > 0`. Requires `f` to have constant term 0.

:::

::: function
void _fmpq_poly_sinh_cosh_series(fmpz * s, fmpz_t sden, fmpz * c, fmpz_t cden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(s, sden, n)` to the series expansion of the hyperbolic
sine of `(f, fden, flen)`, and `(c, cden, n)` to the series
expansion of the hyperbolic cosine.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Supports aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_sinh_cosh_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t f, slong n)

Sets `res1` to the series expansion of the hyperbolic sine of `f`
to order `n > 0`, and `res2` to the series expansion
of the hyperbolic cosine. Requires `f` to have constant term 0.

:::

::: function
void _fmpq_poly_tanh_series(fmpz * g, fmpz_t gden, const fmpz * f, const fmpz_t fden, slong flen, slong n)

Sets `(g, gden, n)` to the series expansion of the
hyperbolic tangent of `(f, fden, flen)`.  Assumes `n > 0` and
that `(f, fden, flen)` has constant term 0.
Does not support aliasing between the input and output polynomials.

:::

::: function
void fmpq_poly_tanh_series(fmpq_poly_t res, const fmpq_poly_t f, slong n)

Sets `res` to the series expansion of the hyperbolic tangent of
`f` to order `n > 0`. Requires `f` to have constant term 0.

:::

::: function
void _fmpq_poly_sin_cos_series_basecase(fmpz * S, fmpz_t Sden, fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden, slong Alen, slong n)
void _fmpq_poly_sin_cos_series_newton(fmpz * S, fmpz_t Sden, fmpz * C, fmpz_t Cden, const fmpz * h, const fmpz_t hden, slong hlen, slong cutoff, slong n)

Alternative algorithms implementing `_fmpq_poly_sin_cos_series`.

:::

## Orthogonal polynomials

::: function
void _fmpq_poly_legendre_p(fmpz * coeffs, fmpz_t den, ulong n)

Sets `coeffs` to the coefficient array of the Legendre polynomial
$P_n(x)$, defined by $(n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)$,
for $n\ge0$. Sets `den` to the overall denominator.
The coefficients are calculated using a hypergeometric recurrence
in integer arithmetic.
The length of the array will be `n+1`. The
denominator is given by
$\gcd(n!,2^n) = 2^{\lfloor n/2 \rfloor + \lfloor n/4 \rfloor + \ldots}.$
See `fmpz_poly` for the shifted Legendre polynomials.

:::

::: function
void fmpq_poly_legendre_p(fmpq_poly_t poly, ulong n)

Sets `poly` to the Legendre polynomial $P_n(x)$, defined
by $(n+1) P_{n+1}(x) = (2n+1) x P_n(x) - n P_{n-1}(x)$, for $n\ge0$.
The coefficients are calculated using a hypergeometric recurrence
in integer arithmetic. The denominator is given by
$\gcd(n!,2^n) = 2^{\lfloor n/2 \rfloor + \lfloor n/4 \rfloor + \ldots}.$
See `fmpz_poly` for the shifted Legendre polynomials.

:::

::: function
void _fmpq_poly_laguerre_l(fmpz * coeffs, fmpz_t den, ulong n)

Sets `coeffs` to the coefficient array of the Laguerre polynomial
$L_n(x)$, defined by $(n+1) L_{n+1}(x) = (2n+1-x) L_n(x) - n L_{n-1}(x)$,
for $n\ge0$. Sets `den` to the overall denominator.
The coefficients are calculated using a hypergeometric recurrence.
The length of the array will be `n+1`.

:::

::: function
void fmpq_poly_laguerre_l(fmpq_poly_t poly, ulong n)

Sets `poly` to the Laguerre polynomial $L_n(x)$, defined by
$(n+1) L_{n+1}(x) = (2n+1-x) L_n(x) - n L_{n-1}(x)$, for $n\ge0$.
The coefficients are calculated using a hypergeometric recurrence.

:::

::: function
void _fmpq_poly_gegenbauer_c(fmpz * coeffs, fmpz_t den, ulong n, const fmpq_t a)

Sets `coeffs` to the coefficient array of the Gegenbauer
(ultraspherical) polynomial
$C^{(\alpha)}_n(x) = \frac{(2\alpha)_n}{n!}\,{}_2F_1\left(-n,2\alpha+n;
\alpha+\frac12;\frac{1-x}{2}\right)$, for integer $n\ge0$ and rational
$\alpha>0$. Sets `den` to the overall denominator.
The coefficients are calculated using a hypergeometric recurrence.

:::

::: function
void fmpq_poly_gegenbauer_c(fmpq_poly_t poly, ulong n, const fmpq_t a)

Sets `poly` to the Gegenbauer (ultraspherical) polynomial
$C^{(\alpha)}_n(x) = \frac{(2\alpha)_n}{n!}\,{}_2F_1\left(-n,2\alpha+n;
\alpha+\frac12;\frac{1-x}{2}\right)$, for integer $n\ge0$ and rational
$\alpha>0$.
The coefficients are calculated using a hypergeometric recurrence.

:::

## Evaluation

::: function
void _fmpq_poly_evaluate_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t a)

Evaluates the polynomial `(poly, den, len)` at the integer $a$ and
sets `(rnum, rden)` to the result in lowest terms.

:::

::: function
void fmpq_poly_evaluate_fmpz(fmpq_t res, const fmpq_poly_t poly, const fmpz_t a)

Evaluates the polynomial `poly` at the integer $a$ and sets
`res` to the result.

:::

::: function
void _fmpq_poly_evaluate_fmpq(fmpz_t rnum, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t anum, const fmpz_t aden)

Evaluates the polynomial `(poly, den, len)` at the rational
`(anum, aden)` and sets `(rnum, rden)` to the result in
lowest terms.  Aliasing between `(rnum, rden)` and
`(anum, aden)` is not supported.

:::

::: function
void fmpq_poly_evaluate_fmpq(fmpq_t res, const fmpq_poly_t poly, const fmpq_t a)

Evaluates the polynomial `poly` at the rational $a$ and
sets `res` to the result.

:::

## Interpolation

::: function
void _fmpq_poly_interpolation_weights(fmpz * w, fmpz_t wden, const fmpq * xs, slong len)
void _fmpq_poly_interpolate_fast_precomp(fmpz * poly, fmpz_t den, const fmpq * ys, fmpz * const * tree, const fmpz * weights, slong len)
void _fmpq_poly_interpolate_fast(fmpz * poly, fmpz_t den, const fmpq * xs, const fmpq * ys, slong len)
void fmpq_poly_interpolate_fast(fmpq_poly_t poly, const fmpq * xs, const fmpq * ys, slong n)

Fast polynomial interpolation using a subproduct tree assuming the evaluation points $x_i$ are distinct.
The *precomp* version requires a precomputed subproduct tree generated using
`_fmpz_poly_tree_build_fmpq_vec` and precomputed interpolation weights
generated using `_fmpq_poly_interpolation_weights`.

:::

::: function
void _fmpq_poly_interpolate_barycentric(fmpz * poly, fmpz_t den, const fmpq * xs, const fmpq * ys, slong n)
void fmpq_poly_interpolate_barycentric(fmpq_poly_t poly, const fmpq * xs, const fmpq * ys, slong n)
void _fmpq_poly_interpolate_multi_mod(fmpz * poly, fmpz_t den, const fmpq * xs, const fmpq * ys, slong n)
void fmpq_poly_interpolate_multi_mod(fmpq_poly_t poly, const fmpq * xs, const fmpq * ys, slong n)
int _fmpq_poly_interpolate_fmpq_vec(fmpz * poly, fmpz_t den, const fmpq * xs, const fmpq * ys, slong n)
int fmpq_poly_interpolate_fmpz_vec(fmpq_poly_t poly, const fmpz * xs, const fmpz * ys, slong n)
int fmpq_poly_interpolate_fmpz_fmpq_vec(fmpq_poly_t poly, const fmpz * xs, const fmpq * ys, slong n)
int fmpq_poly_interpolate_fmpq_vec(fmpq_poly_t poly, const fmpq * xs, const fmpq * ys, slong n)

Sets `poly` (or `poly` / `den`) to the unique interpolating polynomial of
degree at most $n - 1$ satisfying $f(x_i) = y_i$ for every pair $x_i, y_i$
in `xs` and `ys`.

::: details
The *barycentric* uses a naive implementation of Lagrange interpolation in its barycentric form,
clearing denominators to avoid working with fractions.
It assumes that the evaluation points $x_i$ are distinct.

The *multi_mod* algorithm performs interpolation modulo one or more prime numbers, combines
the results using the Chinese Remainder Theorem (CRT) with adaptive termination, and recovers
the rational coefficients by [rational reconstruction](algorithms/number-theory/rational-reconstruction).
This function is particularly efficient when the interpolated polynomial is expected
to have small-bitsize coefficients.
It assumes that the evaluation points $x_i$ are distinct.

The default implementation first checks that the $x_i$ are distinct and, if so, automatically chooses
between the *barycentric*, *multi_mod* and *fast* algorithms.
It returns 1 if the $x_i$ are distinct, otherwise 0.

:::
:::

## Composition

::: function
void _fmpq_poly_compose(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2)

Sets `(res, den)` to the composition of `(poly1, den1, len1)`
and `(poly2, den2, len2)`, assuming `len1, len2 > 0`.

::: details
Assumes that `res` has space for `(len1 - 1) * (len2 - 1) + 1`
coefficients.  Does not support aliasing.

:::
:::

::: function
void fmpq_poly_compose(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)

Sets `res` to the composition of `poly1` and `poly2`.

:::

::: function
void _fmpq_poly_rescale(fmpz * res, fmpz_t denr, const fmpz * poly, const fmpz_t den, slong len, const fmpz_t anum, const fmpz_t aden)

Sets `(res, denr, len)` to `(poly, den, len)` with the
indeterminate rescaled by `(anum, aden)`.

::: details
Assumes that `len > 0` and that `(anum, aden)` is non-zero and
in lowest terms.  Supports aliasing between `(res, denr, len)` and
`(poly, den, len)`.

:::
:::

::: function
void fmpq_poly_rescale(fmpq_poly_t res, const fmpq_poly_t poly, const fmpq_t a)

Sets `res` to `poly` with the indeterminate rescaled by $a$.

:::

## Power series composition {#compose-series}

The functions in this section set the output to the composition
$\mathit{poly1}(\mathit{poly2}) \bmod x^n$, where the constant term of
`poly2` must be zero (the non-underscore versions raise an exception
otherwise). The underscore versions require $\mathit{len1},
\mathit{len2}, n > 0$ and $\mathit{len1}, \mathit{len2} \le n$, write $n$
coefficients and do not allow aliasing; the Horner and Brent–Kung
versions also require $(\mathit{len1} - 1)(\mathit{len2} - 1) + 1 \le n$.
See [](algorithms/polynomials/composition#series) for the algorithms.

::: function
void _fmpq_poly_compose_series(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
void fmpq_poly_compose_series(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)

Compose with an automatic choice of algorithm: if `poly2` has integer
coefficients, the composition is computed over $\Z$ with
`_fmpz_poly_compose_series`; otherwise Horner's rule is used for
$\min(\mathit{len1}, n) < 8$, the Brent–Kung algorithm for $n < 250$ and
the Kinoshita–Li algorithm above.
:::

::: function
void _fmpq_poly_compose_series_horner(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
void fmpq_poly_compose_series_horner(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)

Compose by Horner's rule, with $\mathit{len1} - 1$ truncated
multiplications.
:::

::: function
void _fmpq_poly_compose_series_brent_kung(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
void fmpq_poly_compose_series_brent_kung(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)

Compose by the baby-step giant-step algorithm 2.1 of Brent and Kung
[@BrentKung1978], with about $2\sqrt{n}$ truncated multiplications and a
matrix product over $\Z$.
:::

::: function
void _fmpq_poly_compose_series_kinoshita_li(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, const fmpz * poly2, const fmpz_t den2, slong len2, slong n)
void fmpq_poly_compose_series_kinoshita_li(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2, slong n)

Compose by the algorithm of Kinoshita and Li [@KL2024], which uses
$O(M(n) \log n)$ operations on bivariate polynomials. The underscore
version does not need the condition on $(\mathit{len1} - 1)(\mathit{len2}
- 1) + 1$.
:::

## Power series reversion {#revert-series}

The functions in this section set the output to the compositional
inverse $g$ of `poly` modulo $x^n$, that is, $\mathit{poly}(g(x)) \equiv
g(\mathit{poly}(x)) \equiv x \pmod{x^n}$. `poly` must have constant term
zero and a nonzero coefficient of $x$ (the non-underscore versions raise
an exception otherwise). The underscore versions require $n > 0$ and do
not allow aliasing. See [](algorithms/polynomials/composition) and
[](algorithms/polynomials/newton).

::: function
void _fmpq_poly_revert_series(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, slong n)
void fmpq_poly_revert_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Revert with an automatic choice of algorithm: if `poly` has integer
coefficients and the coefficient of $x$ is $\pm 1$, the reversion is
computed over $\Z$ with `_fmpz_poly_revert_series`; otherwise the fast
Lagrange inversion is used for $n < 20$ and Newton iteration above.
:::

::: function
void _fmpq_poly_revert_series_lagrange(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, slong n)
void fmpq_poly_revert_series_lagrange(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Revert by the Lagrange inversion formula, $[x^k] g = \frac1k [x^{k-1}]
(x/\mathit{poly})^k$, computing the powers one after the other ($O(n M(n))$
operations).
:::

::: function
void _fmpq_poly_revert_series_lagrange_fast(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, slong n)
void fmpq_poly_revert_series_lagrange_fast(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Revert by the Lagrange inversion formula with a baby-step giant-step
evaluation of the needed coefficients of the powers, in $O(\sqrt{n} M(n))$
operations.
:::

::: function
void _fmpq_poly_revert_series_newton(fmpz * res, fmpz_t den, const fmpz * poly1, const fmpz_t den1, slong len1, slong n)
void fmpq_poly_revert_series_newton(fmpq_poly_t res, const fmpq_poly_t poly, slong n)

Revert by Newton iteration on $\mathit{poly}(g) - x = 0$, doubling the
precision in each step, with one power series composition and one
inversion per step.
:::

## Gaussian content

::: function
void _fmpq_poly_content(fmpq_t res, const fmpz * poly, const fmpz_t den, slong len)

Sets `res` to the content of `(poly, den, len)`.
If `len == 0`, sets `res` to zero.

:::

::: function
void fmpq_poly_content(fmpq_t res, const fmpq_poly_t poly)

Sets `res` to the content of `poly`.  The content of the zero
polynomial is defined to be zero.

:::

::: function
void _fmpq_poly_primitive_part(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len)

Sets `(rpoly, rden, len)` to the primitive part, with non-negative
leading coefficient, of `(poly, den, len)`.  Assumes that
`len > 0`.  Supports aliasing between the two polynomials.

:::

::: function
void fmpq_poly_primitive_part(fmpq_poly_t res, const fmpq_poly_t poly)

Sets `res` to the primitive part, with non-negative leading
coefficient, of `poly`.

:::

::: function
int _fmpq_poly_is_monic(const fmpz * poly, const fmpz_t den, slong len)

Returns whether the polynomial `(poly, den, len)` is monic.
The zero polynomial is not monic by definition.

:::

::: function
int fmpq_poly_is_monic(const fmpq_poly_t poly)

Returns whether the polynomial `poly` is monic. The zero
polynomial is not monic by definition.

:::

::: function
void _fmpq_poly_make_monic(fmpz * rpoly, fmpz_t rden, const fmpz * poly, const fmpz_t den, slong len)

Sets `(rpoly, rden, len)` to the monic scalar multiple of
`(poly, den, len)`.  Assumes that `len > 0`.  Supports
aliasing between the two polynomials.

:::

::: function
void fmpq_poly_make_monic(fmpq_poly_t res, const fmpq_poly_t poly)

Sets `res` to the monic scalar multiple of `poly` whenever
`poly` is non-zero.  If `poly` is the zero polynomial, sets
`res` to zero.

:::

## Square-free

::: function
int fmpq_poly_is_squarefree(const fmpq_poly_t poly)

Returns whether the polynomial `poly` is square-free.  A non-zero
polynomial is defined to be square-free if it has no non-unit square
factors.  We also define the zero polynomial to be square-free.

:::

## Input and output

::: function
int _fmpq_poly_print(const fmpz * poly, const fmpz_t den, slong len)

Prints the polynomial `(poly, den, len)` to `stdout`.

::: details
In case of success, returns a positive value.  In case of failure,
returns a non-positive value.

:::
:::

::: function
int fmpq_poly_print(const fmpq_poly_t poly)

Prints the polynomial to `stdout`.

::: details
In case of success, returns a positive value.  In case of failure,
returns a non-positive value.

:::
:::

::: function
int _fmpq_poly_print_pretty(const fmpz * poly, const fmpz_t den, slong len, const char * x)

Prints the pretty representation of `(poly, den, len)` to standard
output, with variable name `x`.
:::

::: function
int fmpq_poly_print_pretty(const fmpq_poly_t poly, const char * var)

Prints the pretty representation of `poly` to `stdout`, using
the null-terminated string `var` not equal to `"\0"` as the
variable name.

::: details
In the current implementation always returns $1$.

:::
:::

::: function
int _fmpq_poly_fprint(FILE * file, const fmpz * poly, const fmpz_t den, slong len)

Prints the polynomial `(poly, den, len)` to the stream `file`.

::: details
In case of success, returns a positive value.  In case of failure,
returns a non-positive value.

:::
:::

::: function
int fmpq_poly_fprint(FILE * file, const fmpq_poly_t poly)

Prints the polynomial to the stream `file`.

::: details
In case of success, returns a positive value.  In case of failure,
returns a non-positive value.

:::
:::

::: function
int _fmpq_poly_fprint_pretty(FILE * file, const fmpz * poly, const fmpz_t den, slong len, const char * x)

Prints the pretty representation of `(poly, den, len)` to `file`, with
variable name `x`.
:::

::: function
int fmpq_poly_fprint_pretty(FILE * file, const fmpq_poly_t poly, const char * var)

Prints the pretty representation of `poly` to `file`, using
the null-terminated string `var` not equal to `"\0"` as the
variable name.

::: details
In the current implementation, always returns $1$.

:::
:::

::: function
int fmpq_poly_read(fmpq_poly_t poly)

Reads a polynomial from `stdin`, storing the result
in `poly`.

::: details
In case of success, returns a positive number.  In case of failure,
returns a non-positive value.

:::
:::

::: function
int fmpq_poly_fread(FILE * file, fmpq_poly_t poly)

Reads a polynomial from the stream `file`, storing the result
in `poly`.

::: details
In case of success, returns a positive number.  In case of failure,
returns a non-positive value.

:::
:::

::: function
int fmpq_poly_debug(const fmpq_poly_t poly)

Prints the internal representation of `poly` (allocation, length, all
allocated numerator coefficients and the denominator) to standard output,
for debugging. Returns 1.
:::
