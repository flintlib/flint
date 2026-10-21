---
title: fmpz_poly.h
header: fmpz_poly.h
description: Dense univariate polynomials over the integers, from basic arithmetic to fast multiplication, division, GCDs, resultants, evaluation, composition, power series, real root isolation and Hensel lifting.
category: Univariate polynomials
order: 1
---

This module implements dense univariate polynomials over $\Z$. An
`fmpz_poly_t` is an array of `fmpz` coefficients together with its
length; each coefficient is an integer of arbitrary size and is managed
separately. Multiplication chooses
among the classical algorithm, word-size and packed-limb kernels,
Karatsuba and Toom–Cook interpolation in the polynomial variable,
[Kronecker substitution](algorithms/polynomials/kronecker) on top of
the FFTs of `fft_small`, and Schönhage–Strassen convolutions, according
to the lengths and the coefficient sizes. GCDs and resultants are
computed by multimodular algorithms that reduce to polynomial
arithmetic modulo word-size primes, division by divide-and-conquer
algorithms, and power series by Newton iteration. Factorisation is in
`fmpz_poly_factor.h`; polynomials with rational coefficients are in
`fmpq_poly.h`, matrices of polynomials in `fmpz_poly_mat.h`, and
rigorous numerical root finding in `arb_fmpz_poly.h` ([](arb_fmpz_poly/roots),
[](algorithms/ball-arithmetic/polynomial-roots)).

Start with the guide [](fmpz_poly/getting-started), which explains the
conventions below with examples.

::: summary algorithms/polynomials/multiplication-integers :::

::: summary algorithms/polynomials/gcd-integers :::

## Types and conventions {#conventions}

::: type
fmpz_poly_struct
fmpz_poly_t

A polynomial over $\Z$. The structure has the fields `coeffs` (an
array of `alloc` coefficients of type `fmpz`), `alloc` and `length`.
An `fmpz_poly_t` is an array of one `fmpz_poly_struct`, so it is passed
by reference. The coefficient of $x^i$ is `coeffs[i]`.
:::

- A polynomial is *normalised* if its length is zero or its leading
  coefficient `coeffs[length - 1]` is nonzero. The zero polynomial has
  length 0 and degree $-1$. All functions require normalised inputs and
  produce normalised outputs, unless stated otherwise.
- The allocated coefficients beyond the length are always zero: every
  function that shortens a polynomial clears the coefficients it drops,
  and newly allocated coefficients are zero. Code that writes
  coefficients directly must keep this invariant.
- Functions without a leading underscore manage memory and allow every
  output to be the same object as any input, unless stated otherwise.
  Several outputs of one function (such as a quotient and a remainder)
  must be distinct objects.
- Functions with a leading underscore work on coefficient arrays
  `(poly, len)` and leave memory management to the caller. They usually
  require positive lengths and an ordering of the lengths, write a
  fixed number of coefficients (including leading zeros), and do not
  allow the output to overlap the inputs. Each function states its
  conditions. *Zero padding* of an input `(poly, len)` means that the
  top coefficients may be zero.
- Division by a polynomial whose leading coefficient is not $\pm 1$
  cannot in general be carried out over $\Z$. The division functions
  either assume that the division is exact, or define the quotient by
  rounding each coefficient (see [](fmpz_poly/divrem)), or work with
  pseudo-division.

## Memory management {#memory}

::: function
void fmpz_poly_init(fmpz_poly_t poly)
void fmpz_poly_init2(fmpz_poly_t poly, slong alloc)

Initialise `poly` as the zero polynomial. `fmpz_poly_init` allocates
nothing; `fmpz_poly_init2` allocates `alloc` coefficients, all zero.
A polynomial must be cleared with `fmpz_poly_clear` after use.
:::

::: function
void fmpz_poly_realloc(fmpz_poly_t poly, slong alloc)

Changes the allocated space of `poly` to exactly `alloc`
coefficients. If the length exceeds `alloc`, the polynomial is first
truncated to length `alloc` and normalised. For `alloc` = 0 the
coefficients are freed and `poly` becomes zero.
:::

::: function
void fmpz_poly_fit_length(fmpz_poly_t poly, slong len)

Ensures that `poly` has space for at least `len` coefficients. If it
has to grow, it grows to at least twice the old allocation, so that
repeated calls with increasing lengths cost amortised linear time. The
length and the coefficients are not changed.
:::

::: function
void fmpz_poly_clear(fmpz_poly_t poly)

Frees the memory used by `poly`, including that of its coefficients.
The polynomial must be initialised again before it is used.
:::

::: function
void _fmpz_poly_set_length(fmpz_poly_t poly, slong newlen)

Sets the length of `poly` to `newlen`. If the length decreases, the
dropped coefficients are set to zero (freeing their memory), which
keeps the invariant that coefficients beyond the length are zero. No
normalisation is done; `poly` must have space for `newlen`
coefficients.
:::

::: function
void _fmpz_poly_normalise(fmpz_poly_t poly)

Decreases the length of `poly` until the leading coefficient is
nonzero or the length is zero. This is the usual last step after
writing coefficients with an underscore function.
:::

::: function
void fmpz_poly_attach_truncate(fmpz_poly_t trunc, const fmpz_poly_t poly, slong n)
void fmpz_poly_attach_shift(fmpz_poly_t trunc, const fmpz_poly_t poly, slong n)

Make the uninitialised `trunc` a read-only view of part of `poly`
without copying: `fmpz_poly_attach_truncate` of the low $\min(n,
\operatorname{len})$ coefficients, `fmpz_poly_attach_shift` of the
coefficients from index $n$ on, that is of $\lfloor \mathit{poly} /
x^n \rfloor$ (zero if $n \ge \operatorname{len}$). The view shares
memory with `poly`; it must not be cleared, modified or used as an
output, and becomes invalid when `poly` changes. A truncated view
may have a zero leading coefficient, so it is not necessarily
normalised.
:::

## Length and degree {#properties}

::: function
slong fmpz_poly_length(const fmpz_poly_t poly)
slong fmpz_poly_degree(const fmpz_poly_t poly)

Return the length and the degree of `poly`. The zero polynomial has
length 0 and degree $-1$.
:::

## Coefficients {#coefficients}

::: function
void fmpz_poly_get_coeff_fmpz(fmpz_t x, const fmpz_poly_t poly, slong n)
slong fmpz_poly_get_coeff_si(const fmpz_poly_t poly, slong n)
ulong fmpz_poly_get_coeff_ui(const fmpz_poly_t poly, slong n)

Return (or set `x` to) the coefficient of $x^n$ in `poly`, which is
zero if $n$ is at least the length. The `si` and `ui` versions return
the result of `fmpz_get_si` or `fmpz_get_ui` and are only meaningful
when the coefficient fits.
:::

::: function
void fmpz_poly_set_coeff_fmpz(fmpz_poly_t poly, slong n, const fmpz_t x)
void fmpz_poly_set_coeff_si(fmpz_poly_t poly, slong n, slong x)
void fmpz_poly_set_coeff_ui(fmpz_poly_t poly, slong n, ulong x)

Set the coefficient of $x^n$ in `poly` to $x$, extending the
polynomial with zero coefficients if $n$ is beyond the length, and
normalising if the leading coefficient is set to zero.
:::

::: macro
fmpz * fmpz_poly_get_coeff_ptr(const fmpz_poly_t poly, slong n)

Returns a pointer to the coefficient of $x^n$, or `NULL` if $n$ is at
least the length. The coefficient can be read or changed with the
functions of [](fmpz); changing the leading coefficient to zero
requires a call to `_fmpz_poly_normalise`.
:::

::: macro
fmpz * fmpz_poly_lead(const fmpz_poly_t poly)

Returns a pointer to the leading coefficient of `poly`, or `NULL` if
`poly` is zero.
:::

## Assignment and basic manipulation {#assignment}

::: function
void fmpz_poly_set(fmpz_poly_t poly1, const fmpz_poly_t poly2)

Sets `poly1` to a copy of `poly2`.
:::

::: function
void fmpz_poly_set_ui(fmpz_poly_t poly, ulong c)
void fmpz_poly_set_si(fmpz_poly_t poly, slong c)
void fmpz_poly_set_fmpz(fmpz_poly_t poly, const fmpz_t c)

Set `poly` to the constant polynomial $c$.
:::

::: function
void fmpz_poly_zero(fmpz_poly_t poly)
void fmpz_poly_one(fmpz_poly_t poly)

Set `poly` to the constant polynomial 0, respectively 1.
:::

::: function
void fmpz_poly_swap(fmpz_poly_t poly1, fmpz_poly_t poly2)

Swaps `poly1` and `poly2` in constant time by swapping their
structures.
:::

::: function
void fmpz_poly_zero_coeffs(fmpz_poly_t poly, slong i, slong j)

Sets the coefficients of $x^i, \ldots, x^{j-1}$ to zero, clamping the
range to the length of `poly`, and normalises.
:::

::: function
void fmpz_poly_truncate(fmpz_poly_t poly, slong newlen)
void fmpz_poly_set_trunc(fmpz_poly_t res, const fmpz_poly_t poly, slong n)

Reduce modulo $x^n$: `fmpz_poly_truncate` truncates `poly` in place
to length at most `newlen`, and `fmpz_poly_set_trunc` sets `res` to
`poly` truncated to length $n$. Both normalise the result.
:::

::: function
void _fmpz_poly_reverse(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_poly_reverse(fmpz_poly_t res, const fmpz_poly_t poly, slong n)

Set `res` to the reversal $x^{n-1} f(1/x)$ of $f$ = `poly` considered
as a polynomial of length $n$, truncating or zero-padding it first.
The underscore version writes $n$ coefficients and requires $0 <
\mathit{len} \le n$; `res` may be equal to `poly` but must not
partially overlap it. The non-underscore version normalises, so the
result can be shorter than $n$.
:::

## Comparison {#comparison}

::: function
int fmpz_poly_equal(const fmpz_poly_t poly1, const fmpz_poly_t poly2)
int fmpz_poly_equal_trunc(const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)

Return whether `poly1` and `poly2` are equal, respectively equal
modulo $x^n$.
:::

::: function
int fmpz_poly_equal_fmpz(const fmpz_poly_t poly, const fmpz_t c)

Returns whether `poly` is the constant polynomial $c$ (the zero
polynomial if $c = 0$).
:::

::: macro
int fmpz_poly_is_zero(const fmpz_poly_t poly)

Returns whether `poly` is zero.
:::

::: function
int _fmpz_poly_is_one(const fmpz * poly, slong len)
int fmpz_poly_is_one(const fmpz_poly_t poly)

Return whether the polynomial is the constant 1. The underscore
version allows zero padding and returns 0 for `len` = 0.
:::

::: function
int fmpz_poly_is_unit(const fmpz_poly_t poly)

Returns whether `poly` is a unit of $\Z[x]$, that is, the constant
$1$ or $-1$.
:::

::: function
int fmpz_poly_is_gen(const fmpz_poly_t poly)

Returns whether `poly` is the generator $x$.
:::

## Randomisation {#random}

::: function
void fmpz_poly_randtest(fmpz_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)
void fmpz_poly_randtest_unsigned(fmpz_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)
void fmpz_poly_randtest_not_zero(fmpz_poly_t f, flint_rand_t state, slong len, flint_bitcnt_t bits)

Set $f$ to a random polynomial of length at most `len` whose
coefficients have at most `bits` bits, with random signs (`randtest`),
nonnegative (`randtest_unsigned`), or, for `randtest_not_zero`, not
the zero polynomial, which requires `len` and `bits` to be positive.
The coefficients are generated by `fmpz_randtest`, which favours
special values (zero, small values, values close to powers of two),
so these functions are meant for testing.
:::

::: function
void fmpz_poly_randtest_no_real_root(fmpz_poly_t p, flint_rand_t state, slong len, flint_bitcnt_t bits)

Sets `p` to a random polynomial of length at most `len` with
coefficients of at most `bits` bits that has no real root.
:::

::: function
void fmpz_poly_randtest_irreducible1(fmpz_poly_t pol, flint_rand_t state, slong len, flint_bitcnt_t bits)
void fmpz_poly_randtest_irreducible2(fmpz_poly_t pol, flint_rand_t state, slong len, flint_bitcnt_t bits)
void fmpz_poly_randtest_irreducible(fmpz_poly_t pol, flint_rand_t state, slong len, flint_bitcnt_t bits)

Set `pol` to a random polynomial of length at most `len` with
coefficients of at most `bits` bits that is irreducible over $\Z$.
`irreducible1` lifts a random irreducible polynomial modulo a random
prime to $\Z$ (an irreducible image modulo a prime that does not
divide the leading coefficient proves irreducibility); `irreducible2`
factors a random polynomial and returns a random factor; the third
function chooses one of the two at random.
:::

## Input and output {#io}

The functions in this section are meant for debugging and for exchanging
small amounts of data, not for speed. There are two formats.

- The *plain* format lists the length followed by the coefficients from
  the constant term up, in decimal and separated by spaces; the length
  is followed by two spaces. For example $5x^3 - x + 1$ is
  `"4  1 -1 0 5"` and the zero polynomial is `"0"`.
- The *pretty* format writes the nonzero terms from the highest degree
  down, as in `5*x^3-x+1`, with a variable name chosen by the caller.
  Coefficients $\pm 1$ are written as a sign only, and the constant term
  without the variable.

::: function
int _fmpz_poly_set_str(fmpz * poly, const char * str)
int fmpz_poly_set_str(fmpz_poly_t poly, const char * str)

Set `poly` to the polynomial given in plain format by the
null-terminated string `str`. Return 0 on success and a nonzero value
if the string is malformed, in which case `fmpz_poly_set_str` sets
`poly` to zero and the underscore version leaves the array in an
undefined state. The underscore version requires that `poly` has room
for the number of coefficients given in the string.
:::

::: function
char * _fmpz_poly_get_str(const fmpz * poly, slong len)
char * fmpz_poly_get_str(const fmpz_poly_t poly)
char * _fmpz_poly_get_str_pretty(const fmpz * poly, slong len, const char * x)
char * fmpz_poly_get_str_pretty(const fmpz_poly_t poly, const char * x)

Return a newly allocated string with the polynomial in plain format,
respectively in pretty format with the variable name `x`. The string
must be freed with `flint_free`.
:::

::: function
int _fmpz_poly_print(const fmpz * poly, slong len)
int fmpz_poly_print(const fmpz_poly_t poly)
int _fmpz_poly_fprint(FILE * file, const fmpz * poly, slong len)
int fmpz_poly_fprint(FILE * file, const fmpz_poly_t poly)

Print the polynomial in plain format to `stdout`, respectively to
`file`. Return a positive value on success and a nonpositive value on
failure.
:::

::: function
int _fmpz_poly_print_pretty(const fmpz * poly, slong len, const char * x)
int fmpz_poly_print_pretty(const fmpz_poly_t poly, const char * x)
int _fmpz_poly_fprint_pretty(FILE * file, const fmpz * poly, slong len, const char * x)
int fmpz_poly_fprint_pretty(FILE * file, const fmpz_poly_t poly, const char * x)

Print the polynomial in pretty format with the variable name `x` to
`stdout`, respectively to `file`. Return a positive value on success
and a nonpositive value on failure.
:::

::: function
int fmpz_poly_read(fmpz_poly_t poly)
int fmpz_poly_fread(FILE * file, fmpz_poly_t poly)

Read a polynomial in plain format from `stdin`, respectively from
`file`, into `poly`. Return a positive value on success and a
nonpositive value on failure.
:::

::: function
int fmpz_poly_read_pretty(fmpz_poly_t poly, char ** x)
int fmpz_poly_fread_pretty(FILE * file, fmpz_poly_t poly, char ** x)

Read a polynomial in pretty format from `stdin`, respectively from
`file`, into `poly`, and set `*x` to a newly allocated string with
the variable name found in the input, to be freed with `flint_free`.
Return the number of characters read on success, and a nonpositive
value on a read error or malformed input.
:::

::: function
void fmpz_poly_debug(const fmpz_poly_t poly)

Prints the fields `alloc` and `length` and all `alloc` allocated
coefficients of `poly` to `stdout`, for debugging.
:::

## Addition and subtraction {#addition}

::: function
void _fmpz_poly_add(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void _fmpz_poly_sub(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)

Set `(res, max(len1, len2))` to the sum, respectively difference, of
`(poly1, len1)` and `(poly2, len2)`. The output may be equal to either
input.
:::

::: function
void fmpz_poly_add(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
void fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
void fmpz_poly_neg(fmpz_poly_t res, const fmpz_poly_t poly)

Set `res` to `poly1 + poly2`, `poly1 - poly2`, respectively `-poly`.
:::

::: function
void fmpz_poly_add_series(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)
void fmpz_poly_sub_series(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)

Set `res` to the sum, respectively difference, of `poly1` and `poly2`
truncated to length $n$.
:::

::: function
void fmpz_poly_add_si(fmpz_poly_t res, const fmpz_poly_t poly, slong c)
void fmpz_poly_add_fmpz(fmpz_poly_t res, const fmpz_poly_t poly, fmpz_t c)
void fmpz_poly_sub_si(fmpz_poly_t res, const fmpz_poly_t poly, slong c)
void fmpz_poly_sub_fmpz(fmpz_poly_t res, const fmpz_poly_t poly, fmpz_t c)
void fmpz_poly_si_sub(fmpz_poly_t res, slong c, const fmpz_poly_t poly)
void fmpz_poly_fmpz_sub(fmpz_poly_t res, fmpz_t c, const fmpz_poly_t poly)

Set `res` to `poly + c`, `poly - c`, respectively `c - poly`, for an
integer constant $c$.
:::

## Scalar multiplication and division {#scalar}

::: function
void fmpz_poly_scalar_abs(fmpz_poly_t res, const fmpz_poly_t poly)

Sets `res` to the polynomial whose coefficients are the absolute values
of those of `poly`.
:::

::: function
void fmpz_poly_scalar_mul_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)
void fmpz_poly_scalar_mul_si(fmpz_poly_t poly1, const fmpz_poly_t poly2, slong x)
void fmpz_poly_scalar_mul_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong x)
void fmpz_poly_scalar_mul_2exp(fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong exp)

Set `poly1` to `poly2` multiplied by $x$, respectively by $2^{exp}$.
:::

::: function
void fmpz_poly_scalar_addmul_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)
void fmpz_poly_scalar_addmul_si(fmpz_poly_t poly1, const fmpz_poly_t poly2, slong x)
void fmpz_poly_scalar_addmul_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong x)
void fmpz_poly_scalar_submul_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)

Set `poly1` to `poly1 + x * poly2`, respectively `poly1 - x * poly2`.
:::

::: function
void fmpz_poly_scalar_fdiv_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)
void fmpz_poly_scalar_fdiv_si(fmpz_poly_t poly1, const fmpz_poly_t poly2, slong x)
void fmpz_poly_scalar_fdiv_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong x)
void fmpz_poly_scalar_tdiv_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)
void fmpz_poly_scalar_tdiv_si(fmpz_poly_t poly1, const fmpz_poly_t poly2, slong x)
void fmpz_poly_scalar_tdiv_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong x)

Set `poly1` to `poly2` divided by $x$, rounding each coefficient
towards $-\infty$ (`fdiv`) or towards zero (`tdiv`). Division by zero
raises an exception.
:::

::: function
void fmpz_poly_scalar_fdiv_2exp(fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong exp)
void fmpz_poly_scalar_tdiv_2exp(fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong exp)

Set `poly1` to `poly2` divided by $2^{exp}$, rounding each coefficient
towards $-\infty$, respectively towards zero.
:::

::: function
void fmpz_poly_scalar_divexact_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)
void fmpz_poly_scalar_divexact_si(fmpz_poly_t poly1, const fmpz_poly_t poly2, slong x)
void fmpz_poly_scalar_divexact_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong x)

Set `poly1` to `poly2` divided by $x$, assuming that $x$ divides every
coefficient. This is not checked; the result is undefined otherwise.
Division by zero raises an exception.
:::

::: function
void fmpz_poly_scalar_mod_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)
void fmpz_poly_scalar_smod_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)

Set `poly1` to `poly2` with each coefficient reduced modulo $x > 0$,
into $[0, x)$, respectively into the symmetric range $(-x/2, x/2]$.
:::

## Shifting, inflation and deflation {#shifting}

::: function
void _fmpz_poly_shift_left(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_poly_shift_left(fmpz_poly_t res, const fmpz_poly_t poly, slong n)

Set `res` to `poly` multiplied by $x^n$. The underscore version writes
`len + n` coefficients; `res` may be equal to `poly`.
:::

::: function
void _fmpz_poly_shift_right(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_poly_shift_right(fmpz_poly_t res, const fmpz_poly_t poly, slong n)

Set `res` to $\lfloor \mathit{poly} / x^n \rfloor$, dropping the low
$n$ coefficients; the result is zero if $n$ is at least the length. The
underscore version requires $\mathit{len} > n$ and writes `len - n`
coefficients; `res` may be equal to `poly`, in which case the top $n$
entries of the array are left unchanged.
:::

::: function
void fmpz_poly_inflate(fmpz_poly_t result, const fmpz_poly_t input, ulong inflation)

Sets `result` to $p(x^n)$, where $p$ is `input` and $n$ = `inflation`.
For $n = 0$ the result is the constant $p(1)$.
:::

::: function
void fmpz_poly_deflate(fmpz_poly_t result, const fmpz_poly_t input, ulong deflation)

Sets `result` to $p(x^{1/n})$, where $p$ is `input` and $n$ =
`deflation`, that is, keeps the coefficients of $x^0, x^n, x^{2n},
\ldots$. Requires that $n \ge 1$ divides `fmpz_poly_deflation(input)`;
other coefficients are dropped and the result is not necessarily
normalised. Division by $n = 0$ raises an exception.
:::

::: function
ulong _fmpz_poly_deflation(const fmpz * a, slong len)
ulong fmpz_poly_deflation(const fmpz_poly_t input)

Return the largest $n$ such that the polynomial is a polynomial in
$x^n$, that is, the GCD of the exponents of its nonzero terms. Return
0 for the zero polynomial and 1 for a nonzero constant.
:::

## Bit sizes and norms {#norms}

::: function
ulong fmpz_poly_max_limbs(const fmpz_poly_t poly)

Returns the largest number of limbs needed to store the absolute value
of a coefficient of `poly`, or 0 if `poly` is zero.
:::

::: function
slong fmpz_poly_max_bits(const fmpz_poly_t poly)

Returns the largest bit length $b$ of the absolute value of a
coefficient of `poly` if all coefficients are nonnegative, and $-b$
otherwise.
:::

::: function
void fmpz_poly_height(fmpz_t res, const fmpz_poly_t poly)

Sets `res` to the height of `poly`, the largest absolute value of a
coefficient (the $\infty$-norm), which is 0 for the zero polynomial.
:::

::: function
void _fmpz_poly_2norm(fmpz_t res, const fmpz * poly, slong len)
void fmpz_poly_2norm(fmpz_t res, const fmpz_poly_t poly)

Set `res` to $\lfloor \|f\|_2 \rfloor$, the integer square root of the
sum of the squares of the coefficients.
:::

::: function
flint_bitcnt_t _fmpz_poly_2norm_normalised_bits(const fmpz * poly, slong len)

Returns an upper bound, at most one too large, for the number of bits
of $\|f\|_2 / |\operatorname{lc}(f)|$, the Euclidean norm of
`(poly, len)` divided by the absolute value of its leading
coefficient. Requires $\mathit{len} > 0$ and a nonzero leading
coefficient. This is used for Landau–Mignotte bounds.
:::

::: function
slong _fmpz_poly_hamming_weight(const fmpz * a, slong len)

Returns the number of nonzero coefficients of `(a, len)`.
:::

## Multiplication {#mul}

All products, truncated products, middle products and squares are
computed by one dispatcher, `_fmpz_poly_mulmid`, which chooses the
algorithm from the lengths, the requested range of coefficients and
the bit sizes of the coefficients; see [](fmpz_poly/mul) and
[](algorithms/polynomials/multiplication-integers). Squaring is
detected when both operands are the same array with the same length.

::: function
void _fmpz_poly_mul(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_mul(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Set `res` to the product of `poly1` and `poly2`. The underscore
version writes `len1 + len2 - 1` coefficients, requires positive
lengths (in either order), allows zero padding and does not allow
`res` to overlap the inputs.
:::

::: function
void _fmpz_poly_mullow(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong n)
void fmpz_poly_mullow(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)

Set `res` to the product of `poly1` and `poly2` truncated to length
$n$, that is, reduced modulo $x^n$. The underscore version writes $n$
coefficients and requires positive lengths (in either order), $0 < n
\le \mathit{len1} + \mathit{len2} - 1$ and no overlap of `res` with
the inputs.
:::

::: function
void _fmpz_poly_mulmid(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi)
void fmpz_poly_mulmid(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong nlo, slong nhi)

Set `res` to the coefficients of $x^{nlo}, \ldots, x^{nhi - 1}$ of
the product of `poly1` and `poly2`, shifted down to start at index 0,
that is, to $\lfloor (\mathit{poly1} \cdot \mathit{poly2} \bmod
x^{nhi}) / x^{nlo} \rfloor$. The non-underscore version clamps
$\mathit{nhi}$ to the length of the product and returns zero for an
empty range. The underscore version writes $\mathit{nhi} -
\mathit{nlo}$ coefficients and requires positive lengths (in either
order), $0 \le \mathit{nlo} < \mathit{nhi} \le \mathit{len1} +
\mathit{len2} - 1$ and no overlap of `res` with the inputs.
:::

::: function
void fmpz_poly_mulhigh_n(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)

Sets `res` to a polynomial of length $\mathit{len1} + \mathit{len2} -
1$ whose top $n$ coefficients are those of the product of `poly1` and
`poly2`; the coefficients below are arbitrary. Both operands are used
at their full lengths.
:::

::: function
void _fmpz_poly_mulhigh(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong start)

Sets the coefficients `res[start]`, ..., `res[len1 + len2 - 2]` to
those of the product of `(poly1, len1)` and `(poly2, len2)`, leaving
`res[0]`, ..., `res[start - 1]` unchanged. Requires positive lengths
(in either order), $0 \le \mathit{start} < \mathit{len1} +
\mathit{len2} - 1$ and no overlap of `res` with the inputs.
:::

::: function
void _fmpz_poly_sqr(fmpz * rop, const fmpz * op, slong len)
void fmpz_poly_sqr(fmpz_poly_t rop, const fmpz_poly_t op)

Set `rop` to the square of `op`. The underscore version writes `2*len
- 1` coefficients, requires $\mathit{len} > 0$, allows zero padding and
does not allow `rop` to overlap `op`.
:::

::: function
void _fmpz_poly_sqrlow(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_poly_sqrlow(fmpz_poly_t res, const fmpz_poly_t poly, slong n)

Set `res` to the square of `poly` truncated to length $n$. The
underscore version writes $n$ coefficients and requires $\mathit{len}
> 0$, $0 < n \le 2 \mathit{len} - 1$ and no overlap of `res` with
`poly`.
:::

## Multiplication algorithms {#mul-algorithms}

The functions in this section implement one algorithm each and are
called by `_fmpz_poly_mulmid`; they are public for testing, tuning and
special applications. Unless stated otherwise, each family has full
(`mul`), low (`mullow`), high (`mulhigh`) and middle (`mulmid`)
products and squares with the same conventions as the functions of the
previous section, and the non-underscore versions allow aliasing. See
[](algorithms/polynomials/multiplication-integers) for the algorithms
and the selection.

::: function
void _fmpz_poly_mul_classical(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_mul_classical(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
void _fmpz_poly_mullow_classical(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong n)
void fmpz_poly_mullow_classical(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)
void _fmpz_poly_mulmid_classical(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi)
void fmpz_poly_mulmid_classical(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong nlo, slong nhi)
void _fmpz_poly_sqr_classical(fmpz * rop, const fmpz * op, slong len)
void fmpz_poly_sqr_classical(fmpz_poly_t rop, const fmpz_poly_t op)
void _fmpz_poly_sqrlow_classical(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_poly_sqrlow_classical(fmpz_poly_t res, const fmpz_poly_t poly, slong n)

Classical (schoolbook) multiplication: each output coefficient is
computed as one dot product of `fmpz` values (`_fmpz_vec_dot_general`),
and squares use the symmetry of the products. The cost is proportional
to the number of coefficient products in the requested range. The
underscore versions require positive lengths and do not allow `res` to
overlap the inputs.
:::

::: function
void _fmpz_poly_mulhigh_classical(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong start)
void fmpz_poly_mulhigh_classical(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong start)

Set the coefficients of `res` from index `start` on to those of the
product of `poly1` and `poly2` by the classical algorithm, and the
coefficients below `start` to zero; the result has length
$\mathit{len1} + \mathit{len2} - 1$, or is zero if `start` is at least
that. The underscore version requires positive lengths and does not
allow `res` to overlap the inputs.
:::

::: function
void _fmpz_poly_mulmid_mpn(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi)
void _fmpz_poly_mul_mpn(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_mulmid_mpn(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong nlo, slong nhi)
void fmpz_poly_mul_mpn(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Middle and full products computed on packed fixed-size limb arrays
instead of `fmpz` values, with classical, Karatsuba and Karatsuba
middle-product algorithms (`_flint_mpn_poly_mulmid`). Intended for
lengths up to a few dozen and coefficients up to about a thousand
bits.

::: details
Coefficients that cannot contribute to the requested window are
trimmed first. `_flint_mpn_poly_mulmid_params` then chooses from the
bit bounds and the lengths how the inputs are packed (as unsigned,
two's complement, sign-magnitude or biased $n$-limb integers; when all
coefficients are small `fmpz` values, the `fmpz` arrays are used
directly as one-limb two's complement integers), and the number of
limbs $s$ of the output; each output coefficient is computed exactly
modulo $2^{\bits s}$ and converted back, with a correction for the
bias when needed. Squaring is detected when `poly1` and `poly2` are the
same array with the same length.
:::
:::

::: function
void _fmpz_poly_mulmid_mpn_bits(fmpz * res, const fmpz * poly1, slong len1, slong bits1, const fmpz * poly2, slong len2, slong bits2, slong nlo, slong nhi)

As `_fmpz_poly_mulmid_mpn`, taking the signed bit bounds `bits1` and
`bits2` of the coefficients, as returned by `_fmpz_vec_max_bits`, from
the caller. The absolute values of the bounds need not be tight, but a
bound must be negative if the polynomial has a negative coefficient,
and a bound 0 means that the polynomial is zero.
:::

::: function
int _fmpz_poly_mul_toom_scalar(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
int _fmpz_poly_mul_toom_scalar_nmax(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nmax)
int fmpz_poly_mul_toom_scalar(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Multiply by Toom–Cook evaluation and interpolation in the polynomial
variable, at the points $0, \pm 1, \ldots, \pm p$, one more point
when the output length is odd, and infinity, using exactly
$\mathit{len1} + \mathit{len2} - 1$ coefficient multiplications.
Return 1 on success and 0, leaving `res` unchanged, if the output
length $\mathit{len1} + \mathit{len2} - 1$ exceeds
`FMPZ_POLY_TOOM_SCALAR_N_MAX` (20 on 64-bit and 13 on 32-bit
machines), respectively $\min(\mathit{nmax},$
`FMPZ_POLY_TOOM_SCALAR_N_MAX`$)$.

::: details
The interpolation is table driven, runs in the output array without
temporary allocation, and divides each row exactly by one word. The
length bound is where the precomputed interpolation constants stop
fitting in a word. Since the number of multiplications is minimal
while the linear work grows quadratically with the length, the method
pays off for short polynomials with large coefficients; the dispatcher
uses it through `_fmpz_poly_mul_toom_karatsuba`, which caps the number
of points according to the coefficient size.

The underscore versions require $\mathit{len1} \ge \mathit{len2} > 0$,
allow zero padding and do not allow `res` to overlap the inputs. When
`poly1` and `poly2` are the same array with the same length, the
evaluations are shared and the pointwise products are squares. The
non-underscore version allows aliasing.
:::
:::

::: function
int _fmpz_poly_mulmid_toom_scalar(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi)
int _fmpz_poly_mulmid_toom_scalar_nmax(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi, slong nmax)
int fmpz_poly_mulmid_toom_scalar(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong nlo, slong nhi)

Compute the middle product window $[\mathit{nlo}, \mathit{nhi})$ as
`_fmpz_poly_mulmid`, by Toom–Cook interpolation. With $m$, $n$ the
lengths of the operands after removing coefficients that cannot
contribute and $w$ the length of the window, the number of coefficient
multiplications is the sum of the two smallest of $m, n, w$, minus
one: either a full product of the trimmed operands, or, transposed, a
middle product of the window against the shorter operand. Return 0,
leaving `res` unchanged, when both routes exceed the length bound of
`_fmpz_poly_mul_toom_scalar_nmax`.

::: details
The count is symmetric in the three lengths, as the transposition
principle requires, and is optimal when the window touches an end of
the product. A narrow window can be supported when the full product
is not: four coefficients of a $16 \times 16$ product need 19
multiplications, the full product 31.

The underscore versions take the lengths in either order and require
$0 \le \mathit{nlo} < \mathit{nhi} \le \mathit{len1} + \mathit{len2}
- 1$ and no overlap of `res` with the inputs.
:::
:::

::: function
void _fmpz_poly_mul_toom_karatsuba(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_mul_toom_karatsuba(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
void _fmpz_poly_mulmid_toom_karatsuba(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi)
void fmpz_poly_mulmid_toom_karatsuba(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong nlo, slong nhi)

Karatsuba multiplication with the Toom–Cook functions
`_fmpz_poly_mul_toom_scalar_nmax` and
`_fmpz_poly_mulmid_toom_scalar_nmax` as base cases, for short
polynomials with large coefficients. Full products use Karatsuba with
the longer operand cut into pieces of the length of the shorter one;
low products use Mulders' short product [@Mul2000], high products a
low product of the reversed operands, and other windows the Karatsuba
middle product of Hanrot, Quercia and Zimmermann [@HQZ2004], so that a
middle-product-shaped window costs about half a full product.

::: details
The base case is used up to a product length that grows with the
coefficient size: 5 below 2500 bits, 11 below 4000 bits, 15 below 6000
bits and `FMPZ_POLY_TOOM_SCALAR_N_MAX` above. The underscore versions
take the lengths in either order, require positive lengths and
$0 \le \mathit{nlo} < \mathit{nhi} \le \mathit{len1} + \mathit{len2}
- 1$, and do not allow `res` to overlap the inputs. The non-underscore
versions clamp the window to the product and allow aliasing.
:::
:::

::: function
void _fmpz_poly_mul_KS(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_mul_KS(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
void _fmpz_poly_mullow_KS(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong n)
void fmpz_poly_mullow_KS(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)
void _fmpz_poly_mulmid_KS(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi)
void fmpz_poly_mulmid_KS(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong nlo, slong nhi)
void _fmpz_poly_sqr_KS(fmpz * rop, const fmpz * op, slong len)
void fmpz_poly_sqr_KS(fmpz_poly_t rop, const fmpz_poly_t op)
void _fmpz_poly_sqrlow_KS(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_poly_sqrlow_KS(fmpz_poly_t res, const fmpz_poly_t poly, slong n)

Multiply by [Kronecker substitution](algorithms/polynomials/kronecker):
the coefficients are packed into bit fields of $b_1 + b_2 + \lceil
\log_2 (\min(\mathit{len1}, \mathit{len2}) + 1) \rceil$ bits, plus one
for signs (where $b_i$ bounds the bits of the coefficients), the two
integers are multiplied with `flint_mpn_mul` (or `flint_mpn_mulmid`
for windows, which computes only the needed limbs), and the product is
unpacked. The underscore versions require positive lengths (in either
order), allow zero padding and zero polynomials, and allow `res` to
be equal to an input, since the inputs are packed before anything is
written.
:::

::: function
void _fmpz_poly_mul_SS(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_mul_SS(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
void _fmpz_poly_mullow_SS(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong n)
void fmpz_poly_mullow_SS(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)
void _fmpz_poly_mulmid_SS(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi)
void fmpz_poly_mulmid_SS(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong nlo, slong nhi)

Multiply by a Schönhage–Strassen convolution
([](algorithms/fft/schonhage-strassen)): each coefficient becomes one
element of $\Z/(2^N + 1)\Z$ large enough for the coefficients of the
product, and the product is a convolution of power-of-two length
computed with `fft_convolution` (or `fft_convolution_sd_fft`, whose
pointwise products use `fft_small`). For a window that does not reach
the top of the product, a cyclic convolution of a power-of-two length
$\ge \max(\mathit{nhi}, \mathit{len1} + \mathit{len2} - 1 -
\mathit{nlo})$ is used when it is shorter than the truncated
transform, so that the wrapped-around terms land below
$\mathit{nlo}$ (`fft_convolution_choose_mid`). The underscore versions
require positive lengths (in either order), allow zero padding and
allow `res` to be equal to an input. The convolutions use
`flint_get_num_threads` threads.
:::

::: function
int _fmpz_poly_mulmid_classical_fft_small(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong nlo, slong nhi)
int fmpz_poly_mulmid_classical_fft_small(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong nlo, slong nhi)

Classical multiplication in which the coefficient products are
computed in the transform space of `fft_small`
([](fft_small/transform-space)): the coefficients of the shorter
operand are transformed once, the longer operand runs through a window
of transformed coefficients, and each output coefficient is a sum of
pointwise products transformed back once. Squares share transforms and
pair symmetric terms. This is efficient for short polynomials with
very large coefficients. Return 1 on success and 0, writing nothing,
when the transform representation is not available (small
coefficients, or a build without `fft_small`). The underscore version
takes the lengths in either order and does not allow `res` to overlap
the inputs.
:::

## Precached multiplication {#precache}

When many polynomials are multiplied by the same polynomial with large
coefficients, its Schönhage–Strassen transform can be computed once.
See also [](fft/precache).

::: type
fmpz_poly_mul_precache_struct
fmpz_poly_mul_precache_t

Stores the transform of a fixed polynomial `poly2` (a copy of which is
kept in the field `poly2`), with the transform length and ring chosen
for products with polynomials of bounded length and coefficient size.
:::

::: function
void fmpz_poly_mul_SS_precache_init(fmpz_poly_mul_precache_t pre, slong len1, slong bits1, const fmpz_poly_t poly2)

Precomputes the transform of `poly2` for multiplications by
polynomials of length at most `len1` whose coefficients have at most
`bits1` bits. The transform length is the power of two $2^k \ge
\mathit{len1} + \mathit{len2} - 1$; products whose length is at most
$2^{k-1}$ are computed with the same, too long, transform. `bits1` may
be negative, as returned by `fmpz_poly_max_bits` for polynomials with
negative coefficients; only its absolute value is used.
:::

::: function
void fmpz_poly_mul_precache_clear(fmpz_poly_mul_precache_t pre)

Frees the memory used by `pre`.
:::

::: function
void _fmpz_poly_mullow_SS_precache(fmpz * output, const fmpz * input1, slong len1, fmpz_poly_mul_precache_t pre, slong trunc)
void fmpz_poly_mullow_SS_precache(fmpz_poly_t res, const fmpz_poly_t poly1, fmpz_poly_mul_precache_t pre, slong n)
void fmpz_poly_mul_SS_precache(fmpz_poly_t res, const fmpz_poly_t poly1, fmpz_poly_mul_precache_t pre)

Set `res` to the product of `poly1` and the polynomial stored in
`pre`, truncated to length $n$ (`trunc`) for the `mullow` versions.
The length and bit size of `poly1` must not exceed those given to
`fmpz_poly_mul_SS_precache_init`. The non-underscore versions use
classical multiplication when one of the lengths or $n$ is at most 2.
The underscore version requires `trunc` to be at most the length of
the product and does not allow `output` to overlap `input1`.
:::

## Bit packing {#bit-packing}

These functions convert between polynomials and integers by evaluation
at $2^b$, the basis of Kronecker substitution.

::: function
void _fmpz_poly_bit_pack(nn_ptr arr, const fmpz * poly, slong len, flint_bitcnt_t bit_size, int negate)

Adds the coefficients of `(poly, len)`, negated if `negate` is $-1$,
into consecutive fields of `bit_size` bits of the zero-initialised
limb array `arr`, representing negative coefficients by borrows from
the next field, so that the array holds $\pm f(2^b)$ as a nonnegative
integer when the leading coefficient times the sign is positive.
`negate` must be 0 or $-1$, and `bit_size` must exceed the bit length
of every coefficient.
:::

::: function
int _fmpz_poly_bit_unpack(fmpz * poly, slong nlo, slong nhi, nn_srcptr arr, flint_bitcnt_t bit_size, int negate)
void _fmpz_poly_bit_unpack_unsigned(fmpz * poly, slong nlo, slong nhi, nn_srcptr arr, flint_bitcnt_t bit_size)

Unpack the signed (respectively unsigned) coefficients in the fields
$\mathit{nlo}, \ldots, \mathit{nhi} - 1$ of `bit_size` bits of
`arr`, written by `_fmpz_poly_bit_pack`, to `poly[0]`, ...,
`poly[nhi - nlo - 1]`, negating them if `negate` is $-1$. The signed
version returns the final borrow: if it is nonzero, a coefficient
$\pm 1$ (with the sign given by `negate`) has to be added at index
`nhi`.
:::

::: function
void fmpz_poly_bit_pack(fmpz_t f, const fmpz_poly_t poly, flint_bitcnt_t bit_size)

Sets $f$ to $\mathit{poly}(2^b)$ for $b$ = `bit_size`, which must
exceed the bit length of every coefficient plus one. The sign of $f$
is that of the leading coefficient.
:::

::: function
void fmpz_poly_bit_unpack(fmpz_poly_t poly, const fmpz_t f, flint_bitcnt_t bit_size)
void fmpz_poly_bit_unpack_unsigned(fmpz_poly_t poly, const fmpz_t f, flint_bitcnt_t bit_size)

Set `poly` to the polynomial with signed (respectively nonnegative)
coefficients of absolute value less than $2^{b-1}$ (respectively
$2^b$) such that $\mathit{poly}(2^b) = f$, for $b$ = `bit_size`.
The unsigned version requires $f \ge 0$.
:::

## Powering {#pow}

::: function
void _fmpz_poly_pow(fmpz * res, const fmpz * poly, slong len, ulong e)
void fmpz_poly_pow(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)

Set `res` to `poly` raised to the power $e$, with $0^0 = 1$. The
function uses repeated squaring and multiplication for $e \le 4$, the
binomial theorem for binomials, the J. C. P. Miller recurrence
(`_fmpz_poly_pow_multinomial`) for short polynomials with small
coefficients compared to $e$, and binary exponentiation otherwise. The
underscore version requires $e > 0$ and $\mathit{len} > 0$, writes
$e(\mathit{len} - 1) + 1$ coefficients and does not allow aliasing.
:::

::: function
void _fmpz_poly_pow_small(fmpz * res, const fmpz * poly, slong len, ulong e)

Sets `res` to `poly` raised to the power $e \le 4$ by squarings and
one multiplication. Requires $\mathit{len} > 0$ and room for
$e(\mathit{len} - 1) + 1$ coefficients; no aliasing.
:::

::: function
void _fmpz_poly_pow_binexp(fmpz * res, const fmpz * poly, slong len, ulong e)
void fmpz_poly_pow_binexp(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)

Set `res` to `poly` raised to the power $e$ by left-to-right binary
exponentiation [@Knu1997, §4.6.3]. The underscore version requires
$\mathit{len} > 0$, $e > 1$
and room for $e(\mathit{len} - 1) + 1$ coefficients, and does not
allow aliasing.
:::

::: function
void _fmpz_poly_pow_addchains(fmpz * res, const fmpz * poly, slong len, const int * a, int n)
void fmpz_poly_pow_addchains(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)

Set `res` to `poly` raised to the power $e$ along an addition chain.
The underscore version takes a star chain $1 = a_0 < a_1 < \cdots <
a_n = e$ (each $a_i$ is $a_{i-1} + a_j$ for some $j < i$) with $n > 1$,
requires $\mathit{len} > 0$ and does not allow aliasing. The
non-underscore version uses a table of chains for $e \le 148$ and
raises an exception for larger $e$.
:::

::: function
void _fmpz_poly_pow_binomial(fmpz * res, const fmpz * poly, ulong e)
void fmpz_poly_pow_binomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)

Set `res` to $(a + bx)^e$ for `poly` $= a + bx$ of length exactly 2,
by the binomial theorem. The non-underscore version raises an
exception for any other length; the underscore version requires $e >
0$ and does not allow aliasing.
:::

::: function
void _fmpz_poly_pow_multinomial(fmpz * res, const fmpz * poly, slong len, ulong e)
void fmpz_poly_pow_multinomial(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)

Set `res` to `poly` raised to the power $e$ using the J. C. P. Miller
recurrence [@Knu1997, §4.6.1, exercise 16], [@Zei1995]: if $P(x) =
\sum_{i=0}^m p_i x^i$ with $p_0 \ne 0$ and $P(x)^e = \sum_k a_k x^k$,
then $a_0 = p_0^e$ and

$$
a_k = \frac{1}{k p_0} \sum_{i=1}^{\min(k, m)} p_i \bigl((e + 1) i - k\bigr) a_{k-i},
\qquad k \ge 1.
$$

Each coefficient costs $O(m)$ multiplications and one exact division.
Low zero coefficients are removed first. The underscore version
requires $\mathit{len} > 0$ and $e > 0$ and does not allow aliasing.
:::

::: function
void _fmpz_poly_pow_trunc(fmpz * res, const fmpz * poly, ulong e, slong n)
void fmpz_poly_pow_trunc(fmpz_poly_t res, const fmpz_poly_t poly, ulong e, slong n)

Set `res` to `poly` raised to the power $e$ and truncated to length
$n$, computed with truncated products throughout, which is much faster
than powering and truncating; with $0^0 = 1$ for $n > 0$. The
underscore version reads $n$ coefficients of `poly` (so it must be zero padded to length
$n$), requires $e, n > 0$ and does not allow aliasing.
:::

## Division {#division}

Over $\Z$, the quotient of $A$ by $B$ is defined coefficient by
coefficient from the top: at each step the leading remaining
coefficient $r$ is divided by $\ell = \operatorname{lc}(B)$ with
floor rounding if $|r| \ge |\ell|$, and the quotient coefficient is
zero otherwise. This gives $A = BQ + R$ in which every coefficient of
$R$ of degree at least $\deg B$ has absolute value less than $|\ell|$.
If $\ell = \pm 1$, or if $B$ divides $A$, this is division in
$\mathbb{Q}[x]$ and $\deg R < \deg B$. See [](fmpz_poly/divrem) for
details, pseudo-division and examples.

::: function
int _fmpz_poly_divrem(fmpz * Q, fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, int exact)
void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)

Set $Q$ and $R$ to the quotient and remainder of $A$ by $B$ as defined
above. Division by zero raises an exception. The underscore version
writes `(Q, lenA - lenB + 1)` and `(R, lenA)`, requires $\mathit{lenA}
\ge \mathit{lenB} > 0$, allows zero padding of $A$ and does not allow
aliasing. It uses the basecase algorithm for $\mathit{lenB} < 6$ and
divide and conquer otherwise.

::: details
If `exact` is 1, the underscore functions of this and the next section
stop and return 0 as soon as a quotient coefficient is not an exact
quotient (the outputs are then undefined), and return 1 otherwise; a
return value 1 does not mean that the remainder is zero, only that the
quotient is that of $\mathbb{Q}[x]$. This is used for divisibility
testing and exact series division. With `exact` = 0 they always return
1.
:::
:::

::: function
int _fmpz_poly_div(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB, int exact)
void fmpz_poly_div(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)

Set $Q$ to the quotient of $A$ by $B$ as defined above, without
computing the remainder. Division by zero raises an exception. The
underscore version writes `lenA - lenB + 1` coefficients, requires
$\mathit{lenA} \ge \mathit{lenB} > 0$, allows zero padding of $A$ and
does not allow aliasing; `exact` is as for `_fmpz_poly_divrem`.
:::

::: function
void _fmpz_poly_rem(fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB)
void fmpz_poly_rem(fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)

Set $R$ to the remainder of $A$ by $B$ as defined above. Division by
zero raises an exception. The underscore version writes `lenA`
coefficients, requires $\mathit{lenA} \ge \mathit{lenB} > 0$ and does
not allow aliasing.
:::

::: function
void _fmpz_poly_divexact(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB)
void fmpz_poly_divexact(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)

Set $Q$ to $A / B$, assuming that $B$ divides $A$; this is not
checked, and the result is undefined otherwise. Division by zero
raises an exception. Short quotients and short divisors are computed
from both ends at once, with coefficients no larger than those of $Q$;
long ones use `_fmpz_poly_div`. The underscore version requires
$\mathit{lenA} \ge \mathit{lenB} > 0$ and does not allow aliasing.
:::

::: function
int _fmpz_poly_divides(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB)
int fmpz_poly_divides(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)

Return whether $B$ divides $A$ in $\Z[x]$, and if so set $Q$ to the
quotient. Inexact divisions are usually detected quickly: by the
constant coefficients, by the values at $x = 1$, and by the first
inexact quotient coefficient. If $B$ does not divide $A$, the value of
$Q$ is undefined. Division by zero raises an exception. The underscore
version requires $\mathit{lenA} \ge \mathit{lenB} > 0$, room for
$\mathit{lenA} - \mathit{lenB} + 1$ coefficients in $Q$ and no
aliasing.
:::

::: function
slong fmpz_poly_remove(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Returns the largest $k$ such that `poly2`$^k$ divides `poly1` and sets
`res` to `poly1` divided by `poly2`$^k$. The value of $k$ is first
bounded using the values at $x = 1$. Raises an exception if `poly2` is
zero or $\pm 1$.
:::

::: function
void _fmpz_poly_div_root_fmpz(fmpz * Q, const fmpz * A, slong len, const fmpz_t c)
void fmpz_poly_div_root_fmpz(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_t c)

Set $Q$ to the quotient of $A$ by $x - c$ by synthetic division
(Horner's rule); the remainder $A(c)$ is discarded. The underscore
version writes `len - 1` coefficients; $Q$ may be equal to $A$ but must
not partially overlap it.
:::

::: function
void _fmpz_poly_divexact_root_fmpq(fmpz * Q, const fmpz * A, slong len, const fmpq_t c)
void fmpz_poly_divexact_root_fmpq(fmpz_poly_t Q, const fmpz_poly_t A, const fmpq_t c)

Set $Q$ to $A / (qx - p)$ where $c = p/q$ in lowest terms, assuming
that $c$ is a root of $A$; this is not checked. The underscore version
writes `len - 1` coefficients; $Q$ may be equal to $A$ but must not
partially overlap it.
:::

::: function
void _fmpz_poly_pseudo_divrem(fmpz * Q, fmpz * R, ulong * d, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_preinvn_t inv)
void fmpz_poly_pseudo_divrem(fmpz_poly_t Q, fmpz_poly_t R, ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)

Set $Q$, $R$ and $d$ such that $\ell^d A = BQ + R$ with $\deg R <
\deg B$, where $\ell = \operatorname{lc}(B)$ and $0 \le d \le
\mathit{lenA} - \mathit{lenB} + 1$. Factors $\ell$ are introduced
only where a quotient coefficient would otherwise not be an integer,
so $d$ is small when the division is nearly exact. These are the
divide-and-conquer functions `fmpz_poly_pseudo_divrem_divconquer`
under another name. `Q` and `R` must be distinct, and division by zero
raises an exception.

::: details
The underscore version requires $\mathit{lenA} \ge \mathit{lenB} > 0$,
room for $\mathit{lenA} - \mathit{lenB} + 1$ coefficients in $Q$ and
$\mathit{lenA}$ in $R$, of which only the low $\mathit{lenB} - 1$ are
meaningful on exit, and does not allow aliasing. `inv` is `NULL` or a
precomputed inverse of $\ell$ from `fmpz_preinvn_init`. The
underscore pseudo-division functions are declared only if `fmpz.h` is
included before `fmpz_poly.h`.
:::
:::

::: function
void _fmpz_poly_pseudo_div(fmpz * Q, ulong * d, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_preinvn_t inv)
void fmpz_poly_pseudo_div(fmpz_poly_t Q, ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)
void _fmpz_poly_pseudo_rem(fmpz * R, ulong * d, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_preinvn_t inv)
void fmpz_poly_pseudo_rem(fmpz_poly_t R, ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)

As `fmpz_poly_pseudo_divrem`, returning only the quotient,
respectively only the remainder, together with $d$.
:::

## Division algorithms {#division-algorithms}

::: function
int _fmpz_poly_divrem_basecase(fmpz * Q, fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, int exact)
void fmpz_poly_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)

Division with remainder by the schoolbook algorithm, with the
semantics of `fmpz_poly_divrem`. The underscore version requires
positive lengths, allows zero padding of $A$, and allows $R$ to be
equal to $A$; no other aliasing is allowed.
:::

::: function
int _fmpz_poly_div_basecase(fmpz * Q, fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, int exact)
void fmpz_poly_div_basecase(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)

Quotient by the schoolbook algorithm, with the semantics of
`fmpz_poly_div`; only the coefficients of the remainder that are
needed for the quotient are updated. The underscore version takes a
temporary array $R$ of length `lenA`, or `NULL` to allocate one; $R$
may be equal to $A$; no other aliasing is allowed.
:::

::: function
int _fmpz_poly_divrem_divconquer(fmpz * Q, fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB, int exact)
void fmpz_poly_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
int _fmpz_poly_div_divconquer(fmpz * Q, const fmpz * A, slong lenA, const fmpz * B, slong lenB, int exact)
void fmpz_poly_div_divconquer(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)

Division with remainder, respectively quotient only, by recursive
division of the top halves, so that the cost is $O(M(n) \log n)$
coefficient operations for a product cost $M(n)$, with the semantics
of `fmpz_poly_divrem` and `fmpz_poly_div`. The underscore versions
require $\mathit{lenA} \ge \mathit{lenB} > 0$, allow zero padding of
$A$ and do not allow aliasing.
:::

::: function
int _fmpz_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, fmpz * W, const fmpz * A, const fmpz * B, slong lenB, int exact)
int _fmpz_poly_divremlow_divconquer_recursive(fmpz * Q, fmpz * QB, const fmpz * A, const fmpz * B, slong lenB, int exact)
int _fmpz_poly_div_divconquer_recursive(fmpz * Q, fmpz * temp, const fmpz * A, const fmpz * B, slong lenB, int exact)

The recursive steps of divide-and-conquer division for the balanced
case $\operatorname{len}(A) = 2 \mathit{lenB} - 1$: they compute the
quotient `(Q, lenB)`, and `divrem` also the product `(BQ, 2 lenB - 1)`
using the temporary `(W, 2 lenB - 1)`; `divremlow` computes only the
low $\mathit{lenB} - 1$ coefficients of $BQ$ correctly (the array
must still have room for $2 \mathit{lenB} - 1$); `div` uses a
temporary `(temp, 2 lenB - 1)` and does not compute $BQ$ [@Mul2000].
The low $\mathit{lenB} - 1$ coefficients of $A$ are not read. They
require $\mathit{lenB} > 0$ and do not allow aliasing; `exact` is as
for `_fmpz_poly_divrem`.
:::

::: function
void _fmpz_poly_rem_basecase(fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB)
void fmpz_poly_rem_basecase(fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)

Remainder by the schoolbook algorithm, with the semantics of
`fmpz_poly_rem`. The underscore version writes `lenA` coefficients and
allows $R$ to be equal to $A$; no other aliasing is allowed.
:::

::: function
void _fmpz_poly_preinvert(fmpz * B_inv, const fmpz * B, slong n)
void fmpz_poly_preinvert(fmpz_poly_t B_inv, const fmpz_poly_t B)

Set `B_inv` to the inverse of the reversal of $B$ as a power series
of length $n = \operatorname{len}(B)$, by Newton iteration, for use
with `fmpz_poly_div_preinv` and `fmpz_poly_divrem_preinv`. $B$ must
have leading coefficient $\pm 1$ (for example, be monic). Division by
zero raises an exception. The underscore version requires $n > 0$ and
no aliasing.
:::

::: function
void _fmpz_poly_div_preinv(fmpz * Q, const fmpz * A, slong len1, const fmpz * B, const fmpz * B_inv, slong len2)
void fmpz_poly_div_preinv(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B, const fmpz_poly_t B_inv)

Set $Q$ to the quotient of $A$ by $B$, given the precomputed inverse
`B_inv` of $B$, using a truncated product with the inverse (and
splitting $A$ into blocks if its quotient is longer than $B$). $B$ must
have leading coefficient $\pm 1$. The underscore version requires
$\mathit{len1} \ge \mathit{len2}$, room for $\mathit{len1} -
\mathit{len2} + 1$ coefficients in $Q$, and no aliasing.
:::

::: function
void _fmpz_poly_divrem_preinv(fmpz * Q, fmpz * A, slong len1, const fmpz * B, const fmpz * B_inv, slong len2)
void fmpz_poly_divrem_preinv(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B, const fmpz_poly_t B_inv)

As `fmpz_poly_div_preinv`, also computing the remainder. The
underscore version replaces the low $\mathit{len2} - 1$ coefficients
of $A$ by the remainder.
:::

::: type
fmpz_poly_powers_precomp_struct
fmpz_poly_powers_precomp_t

Stores the residues $x^i \bmod B$ for $0 \le i < 2 \operatorname{len}(B)
- 1$, for repeated reduction modulo a fixed monic polynomial $B$.
:::

::: function
fmpz ** _fmpz_poly_powers_precompute(const fmpz * B, slong len)
void fmpz_poly_powers_precompute(fmpz_poly_powers_precomp_t pinv, fmpz_poly_t poly)

Compute the residues of $x^i$ modulo the monic polynomial $B$ of
length `len` for $0 \le i < 2 \mathit{len} - 1$, each as a vector of
$\mathit{len} - 1$ coefficients. Division by zero raises an
exception.
:::

::: function
void _fmpz_poly_powers_clear(fmpz ** powers, slong len)
void fmpz_poly_powers_clear(fmpz_poly_powers_precomp_t pinv)

Free the residues computed by `_fmpz_poly_powers_precompute`,
respectively `fmpz_poly_powers_precompute`.
:::

::: function
void _fmpz_poly_rem_powers_precomp(fmpz * A, slong m, const fmpz * B, slong n, fmpz ** const powers)
void fmpz_poly_rem_powers_precomp(fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B, const fmpz_poly_powers_precomp_t B_inv)

Set $R$ to $A \bmod B$ for the monic $B$ of length $n$, as a linear
combination of the precomputed residues of $x^i$; this costs
$O(\operatorname{len}(A) \cdot n)$ coefficient operations and is fast
for $\operatorname{len}(A) < 2n$. The underscore version reduces
`(A, m)` in place, leaving the remainder in its low $n - 1$
coefficients, and does not allow aliasing.
:::

::: function
void fmpz_poly_divlow_smodp(fmpz * res, const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_t p, slong n)

Sets `(res, n)` to the low $n$ coefficients of $f / g$ modulo $p$,
reduced into the symmetric range, assuming that the division is exact
modulo $p$. Common low zero coefficients are removed first; the lowest
nonzero coefficient of $g$ must be invertible modulo $p$, which is not
checked. Used by `_fmpz_poly_factor_CLD_mat`.
:::

::: function
void fmpz_poly_divhigh_smodp(fmpz * res, const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_t p, slong n)

Sets `(res, n)` to the high $n$ coefficients of $f / g$ modulo $p$,
reduced into the symmetric range, where $f$ is the output of
`fmpz_poly_mulhigh_n` applied to $g$ and a polynomial of length $n$
(only the top $n$ coefficients of $f$ are read). The leading
coefficient of $g$ must be invertible modulo $p$, which is not
checked. Used by `_fmpz_poly_factor_CLD_mat`.
:::

::: function
void _fmpz_poly_pseudo_divrem_basecase(fmpz * Q, fmpz * R, ulong * d, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_preinvn_t inv)
void fmpz_poly_pseudo_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R, ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)
void _fmpz_poly_pseudo_divrem_divconquer(fmpz * Q, fmpz * R, ulong * d, const fmpz * A, slong lenA, const fmpz * B, slong lenB, const fmpz_preinvn_t inv)
void fmpz_poly_pseudo_divrem_divconquer(fmpz_poly_t Q, fmpz_poly_t R, ulong * d, const fmpz_poly_t A, const fmpz_poly_t B)

Pseudo-division with the semantics of `fmpz_poly_pseudo_divrem`, by
the schoolbook algorithm, which multiplies the partial remainder by
$\ell$ only when a quotient coefficient would be inexact, respectively
by divide and conquer. The underscore basecase version allows $R$ to
be equal to $A$; the divide-and-conquer version only sets the low
$\mathit{lenB} - 1$ coefficients of `(R, lenA)` correctly and allows
no aliasing. `inv` is as for `_fmpz_poly_pseudo_divrem`.
:::

::: function
void _fmpz_poly_pseudo_divrem_cohen(fmpz * Q, fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB)
void fmpz_poly_pseudo_divrem_cohen(fmpz_poly_t Q, fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
void _fmpz_poly_pseudo_rem_cohen(fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB)
void fmpz_poly_pseudo_rem_cohen(fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)

Pseudo-division with the exponent fixed at $d = \max(0,
\operatorname{len}(A) - \operatorname{len}(B) + 1)$, by Algorithm
3.1.2 of [@Coh1996]: $\ell^d A = BQ + R$ with $\deg R < \deg B$. This
is the pseudo-remainder of the subresultant algorithm; it is quadratic
and meant for short polynomials. The underscore versions require
$\mathit{lenA} \ge \mathit{lenB} > 0$, room for $\mathit{lenA}$
coefficients in $R$, and allow $R$ to be equal to $A$; no other
aliasing.
:::

## Content and primitive part {#content}

The *content* of a nonzero $f \in \Z[x]$ is the nonnegative GCD of
its coefficients, and the *primitive part* is $f$ divided by its
content, with the sign chosen so that the leading coefficient is
positive. By Gauss' lemma, the content of a product is the product of
the contents.

::: function
void _fmpz_poly_content(fmpz_t res, const fmpz * poly, slong len)
void fmpz_poly_content(fmpz_t res, const fmpz_poly_t poly)

Set `res` to the nonnegative content of the polynomial; the content
of the zero polynomial is zero. `res` may be one of the coefficients
of `poly` in the non-underscore version, but not in the underscore
version.
:::

::: function
void _fmpz_poly_primitive_part(fmpz * res, const fmpz * poly, slong len)
void fmpz_poly_primitive_part(fmpz_poly_t res, const fmpz_poly_t poly)

Set `res` to the primitive part of `poly`, with positive leading
coefficient; the primitive part of zero is zero. The underscore version
requires a nonzero polynomial and allows `res` to be equal to `poly`.
:::

## Greatest common divisors {#gcd}

The GCD of two polynomials in $\Z[x]$ is defined up to sign; these
functions return the one with positive leading coefficient, so that
$\gcd(f, g) = \gcd(\operatorname{cont} f, \operatorname{cont} g) \cdot
\gcd(\operatorname{pp} f, \operatorname{pp} g)$. The GCD of $f$ and
$0$ is $\pm f$, and that of $0$ and $0$ is $0$. See [](fmpz_poly/gcd)
and [](algorithms/polynomials/gcd-integers).

::: function
void _fmpz_poly_gcd(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_gcd(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Set `res` to the GCD of `poly1` and `poly2` with positive leading
coefficient. Common powers of $x$ are removed first. Polynomials of
length less than 6 use the subresultant algorithm; otherwise, when the
coefficients have at most 40 bits in total, the heuristic GCD is tried,
and the multimodular algorithm is used if it fails or the
coefficients are larger. The underscore version requires
$\mathit{len1} \ge \mathit{len2} > 0$ and writes `len2` coefficients
(zero padded); `res` may be equal to an input.
:::

::: function
void _fmpz_poly_gcd_modular(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_gcd_modular(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

GCD by the multimodular algorithm: the primitive parts are reduced
modulo batches of word-size primes, their monic GCDs modulo each prime
are scaled to a common leading (or constant) coefficient and combined
by the Chinese remainder theorem, and the candidate is certified by
checking that it divides both inputs. The result is always correct.
The underscore version has the conventions of `_fmpz_poly_gcd`.

::: details
The batches grow geometrically and the modular GCDs of a batch are
computed in parallel when threads are available. A full
reconstruction is attempted only when a random signed sum of the
images has stabilised, so the number of primes adapts to the size of
the GCD rather than to an a priori bound. The divisibility checks are
also multimodular: the quotients are computed modulo the primes
already used, and the product of candidate and quotient is checked by
coefficient bounds or by one multiplication. The first prime is a
20-bit prime unless the coefficients are small, since it usually
suffices to show that the GCD is 1; further primes are 62-bit primes,
or 50-bit primes that suit `fft_small` for polynomials of length at
least 2500. See [](algorithms/polynomials/gcd-integers#modular).
:::
:::

::: function
int _fmpz_poly_gcd_heuristic(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
int fmpz_poly_gcd_heuristic(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Attempt to compute the GCD by the heuristic GCD algorithm GCDHEU
[@CGG1989], evaluating the primitive parts at a power of two, taking
the GCD of the two integers and reconstructing a polynomial from it.
Return 1 on success, when the reconstructed polynomial is verified to
divide both inputs, and 0 on failure, in which case `res` is
undefined. The underscore version has the conventions of
`_fmpz_poly_gcd`.
:::

::: function
void _fmpz_poly_gcd_subresultant(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_gcd_subresultant(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

GCD by the subresultant pseudo-remainder sequence, Algorithm 3.3.1
of [@Coh1996], which keeps the intermediate coefficients small by
exact divisions. It is quadratic in the degree and fast only for short
polynomials. The underscore version has the conventions of
`_fmpz_poly_gcd`.
:::

::: function
void _fmpz_poly_xgcd(fmpz_t r, fmpz * s, fmpz * t, const fmpz * f, slong len1, const fmpz * g, slong len2)
void fmpz_poly_xgcd(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t, const fmpz_poly_t f, const fmpz_poly_t g)
void _fmpz_poly_xgcd_modular(fmpz_t r, fmpz * s, fmpz * t, const fmpz * f, slong len1, const fmpz * g, slong len2)
void fmpz_poly_xgcd_modular(fmpz_t r, fmpz_poly_t s, fmpz_poly_t t, const fmpz_poly_t f, const fmpz_poly_t g)

Set $r$ to the resultant of $f$ and $g$ and, if $r \ne 0$, find $s$,
$t$ with $\deg s < \deg g$, $\deg t < \deg f$ and $sf + tg = r$. If
$r = 0$ the function returns after setting $r$; $s$ and $t$ are then
not computed and their values are undefined. The functions without
`_modular` are inline wrappers of the multimodular algorithm, which
computes the extended GCD modulo word-size primes and combines the
cofactors by the Chinese remainder theorem until the identity is
proved by a size bound.

::: details
If both $f$ and $g$ are nonzero, $s$ and $t$ are unique and their
coefficients are minors of the Sylvester matrix, so they are integers.
$f$ and $g$ need not be primitive.

The underscore versions require $\mathit{len1} \ge \mathit{len2} > 0$,
write $\mathit{len2}$ coefficients to $s$ and $\mathit{len1}$ to $t$
(zero padded), and do not allow aliasing. When
$\operatorname{len}(f) < \operatorname{len}(g)$, the non-underscore
functions swap the inputs and correct the sign, so that $r$ is
$\operatorname{res}(f, g)$ in all cases.
:::
:::

::: function
void _fmpz_poly_lcm(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_lcm(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Set `res` to the least common multiple $fg / \gcd(f, g)$ of the two
polynomials with nonnegative leading coefficient, so that $fg = \pm
\gcd(f, g) \operatorname{lcm}(f, g)$. The LCM is zero if either input
is zero. The underscore version requires $\mathit{len1} \ge
\mathit{len2} > 0$, writes $\mathit{len1} + \mathit{len2} - 1$
coefficients (zero padded) and does not allow aliasing.
:::

## Resultants and discriminants {#resultant}

For nonzero $f = a \prod_{i=1}^m (x - \alpha_i)$ and $g = b
\prod_{j=1}^n (x - \beta_j)$, the resultant is

$$
\operatorname{res}(f, g) = a^n b^m \prod_{i, j} (\alpha_i - \beta_j)
= a^n \prod_i g(\alpha_i),
$$

the determinant of the Sylvester matrix; it is zero if either
polynomial is zero. See [](fmpz_poly/resultant) and
[](algorithms/polynomials/resultant-integers).

::: function
void _fmpz_poly_resultant(fmpz_t res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_resultant(fmpz_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Set `res` to the resultant of `poly1` and `poly2`. The function uses
the subresultant algorithm (`fmpz_poly_resultant_euclidean`) for short
polynomials with small coefficients and the multimodular algorithm
otherwise: when $\mathit{len2} > 144$ or $\mathit{len2}^3
(b_1 + b_2) > 6 \cdot 10^6$, where $b_1$, $b_2$ are the bit sizes of
the coefficients. The underscore version requires $\mathit{len1} \ge
\mathit{len2} > 0$.
:::

::: function
void _fmpz_poly_resultant_modular(fmpz_t res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_resultant_modular(fmpz_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Resultant by the multimodular algorithm [@Col1971]: the contents are
removed, the resultant of the primitive parts is computed modulo
enough 64-bit primes to exceed twice the Hadamard bound
$\|f\|_2^{\deg g} \|g\|_2^{\deg f}$, and the images are combined by
the Chinese remainder theorem. The underscore version requires
$\mathit{len1} \ge \mathit{len2} > 0$.
:::

::: function
void _fmpz_poly_resultant_modular_div(fmpz_t res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, const fmpz_t divisor, slong nbits)
void fmpz_poly_resultant_modular_div(fmpz_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t divisor, slong nbits)

Set `res` to the resultant of `poly1` and `poly2` divided by
`divisor`, assuming that the division is exact and that the quotient
has at most `nbits` bits in absolute value. Only enough primes for
`nbits` bits are used, instead of the general bound; this is used
when the size of the result is known in advance, as for discriminants
of number fields. If `divisor` is zero, `res` is set to zero. The
underscore version requires $\mathit{len1} \ge \mathit{len2} > 0$.
The number of primes is chosen so that their product exceeds
$2^{\mathit{nbits} + 1}$, as the signed reconstruction requires.
:::

::: function
void _fmpz_poly_resultant_euclidean(fmpz_t res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_resultant_euclidean(fmpz_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Resultant by the subresultant algorithm, Algorithm 3.3.7 of
[@Coh1996], with exact divisions that keep the coefficients of the
remainders of the size of minors of the Sylvester matrix. The
underscore version requires $\mathit{len1} \ge \mathit{len2} > 0$.
:::

::: function
void _fmpz_poly_discriminant(fmpz_t res, const fmpz * poly, slong len)
void fmpz_poly_discriminant(fmpz_t res, const fmpz_poly_t poly)

Set `res` to the discriminant of $f$ = `poly` of degree $n$,

$$
\operatorname{disc}(f) = \frac{(-1)^{n(n-1)/2}}{\operatorname{lc}(f)} \operatorname{res}(f, f')
= \operatorname{lc}(f)^{2n-2} \prod_{i < j} (\alpha_i - \alpha_j)^2,
$$

where $\alpha_1, \ldots, \alpha_n$ are the roots of $f$. The
discriminant of a polynomial of degree 1 is 1, and that of a constant
or zero polynomial is defined to be 0. The underscore version requires
$\mathit{len} > 1$.
:::

## Squarefree polynomials {#squarefree}

::: function
int _fmpz_poly_is_squarefree(const fmpz * poly, slong len)
int fmpz_poly_is_squarefree(const fmpz_poly_t poly)

Return whether the polynomial has no repeated factor of positive
degree, that is, whether it is squarefree in $\mathbb{Q}[x]$; the
content is ignored, so $4x + 4$ counts as squarefree. Polynomials of
length at most 2 are squarefree, including zero; quadratics are
tested with the discriminant, and longer polynomials by checking that
$\gcd(f, f')$ is constant.
:::

::: function
void fmpz_poly_squarefree_part(fmpz_poly_t res, const fmpz_poly_t poly)

Sets `res` to $f / \gcd(f, f')$ for $f$ = `poly`, with positive
leading coefficient: the product of the distinct irreducible factors
of positive degree of $f$, times a constant. The squarefree part of a
nonzero constant is 1 and that of zero is zero.
:::

## Derivatives {#derivative}

::: function
void _fmpz_poly_derivative(fmpz * rpoly, const fmpz * poly, slong len)
void fmpz_poly_derivative(fmpz_poly_t res, const fmpz_poly_t poly)

Set `res` to the derivative of `poly`. The underscore version writes
`len - 1` coefficients, handles $\mathit{len} \le 1$, and allows
`rpoly` to be equal to `poly`.
:::

::: function
void _fmpz_poly_nth_derivative(fmpz * rpoly, const fmpz * poly, ulong n, slong len)
void fmpz_poly_nth_derivative(fmpz_poly_t res, const fmpz_poly_t poly, ulong n)

Set `res` to the $n$-th derivative of `poly`. The underscore version
writes `len - n` coefficients, handles $\mathit{len} \le n$, and
allows `rpoly` to be equal to `poly`.
:::

## Evaluation {#evaluation}

::: function
void _fmpz_poly_evaluate_fmpz(fmpz_t res, const fmpz * f, slong len, const fmpz_t a)
void fmpz_poly_evaluate_fmpz(fmpz_t res, const fmpz_poly_t f, const fmpz_t a)

Set `res` to $f(a)$ for an integer $a$, by Horner's rule for length at
most 50 and by the divide-and-conquer algorithm otherwise. In the
non-underscore version `res` may be equal to $a$; in both versions
`res` must not be a coefficient of $f$, and in the underscore version
it must not be $a$. See [](fmpz_poly/evaluate).
:::

::: function
void _fmpz_poly_evaluate_horner_fmpz(fmpz_t res, const fmpz * f, slong len, const fmpz_t a)
void fmpz_poly_evaluate_horner_fmpz(fmpz_t res, const fmpz_poly_t f, const fmpz_t a)
void _fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz * poly, slong len, const fmpz_t a)
void fmpz_poly_evaluate_divconquer_fmpz(fmpz_t res, const fmpz_poly_t poly, const fmpz_t a)

Evaluate at an integer by Horner's rule, respectively by splitting
$f = f_0 + a^m f_1$ recursively, which balances the sizes of the
integer products so that they can use fast multiplication. The
aliasing rules are those of `fmpz_poly_evaluate_fmpz`; the
divide-and-conquer underscore version requires $\mathit{len} \ge 1$
and allows zero padding.
:::

::: function
void _fmpz_poly_evaluate_fmpq(fmpz_t rnum, fmpz_t rden, const fmpz * f, slong len, const fmpz_t anum, const fmpz_t aden)
void fmpz_poly_evaluate_fmpq(fmpq_t res, const fmpz_poly_t f, const fmpq_t a)

Set `res` to $f(a)$ for a rational $a$, in lowest terms, by Horner's
rule for short polynomials or large denominators, and by divide and
conquer otherwise. The underscore version takes $a$ =
`anum/aden` in lowest terms with positive denominator and sets
`(rnum, rden)`; the outputs must not be aliased with the inputs or the
coefficients. The non-underscore version allows `res` to be equal to
$a$.
:::

::: function
void _fmpz_poly_evaluate_horner_fmpq(fmpz_t rnum, fmpz_t rden, const fmpz * f, slong len, const fmpz_t anum, const fmpz_t aden)
void fmpz_poly_evaluate_horner_fmpq(fmpq_t res, const fmpz_poly_t f, const fmpq_t a)
void _fmpz_poly_evaluate_divconquer_fmpq(fmpz_t rnum, fmpz_t rden, const fmpz * f, slong len, const fmpz_t anum, const fmpz_t aden)
void fmpz_poly_evaluate_divconquer_fmpq(fmpq_t res, const fmpz_poly_t f, const fmpq_t a)

Evaluate at a rational number by Horner's rule, respectively by
divide and conquer, with the conventions of `fmpz_poly_evaluate_fmpq`.
The divide-and-conquer underscore version requires $\mathit{len} \ge
1$.
:::

::: function
ulong _fmpz_poly_evaluate_mod(const fmpz * poly, slong len, ulong a, ulong n, ulong ninv)
ulong fmpz_poly_evaluate_mod(const fmpz_poly_t poly, ulong a, ulong n)

Return $f(a) \bmod n$ for a word-size modulus $n \ge 1$, by Horner's
rule with each coefficient reduced modulo $n$. The underscore version
takes the precomputed inverse `ninv = n_preinvert_limb(n)`.
:::

::: function
void fmpz_poly_evaluate_fmpz_vec(fmpz * res, const fmpz_poly_t f, const fmpz * a, slong n)

Sets `res[i]` to $f(a_i)$ for the $n$ integers in `a`, one at a time
with `fmpz_poly_evaluate_fmpz`.
:::

::: function
double _fmpz_poly_evaluate_horner_d(const fmpz * poly, slong n, double d)
double fmpz_poly_evaluate_horner_d(const fmpz_poly_t poly, double d)

Return an approximation of $f(d)$ computed by Horner's rule in double
precision, with no control of the error and no protection against
overflow. Meant for quick estimates, for example for polynomials with
nonnegative coefficients evaluated at $d > 0$.
:::

::: function
double _fmpz_poly_evaluate_horner_d_2exp(slong * exp, const fmpz * poly, slong n, double d)
double fmpz_poly_evaluate_horner_d_2exp(slong * exp, const fmpz_poly_t poly, double d)
double _fmpz_poly_evaluate_horner_d_2exp2(slong * exp, const fmpz * poly, slong n, double d, slong dexp)
double fmpz_poly_evaluate_horner_d_2exp2(slong * exp, const fmpz_poly_t poly, double d, slong dexp)

Return $m$ and set `*exp` to $e$ such that $m 2^e$ approximates $f(d)$,
respectively $f(d \cdot 2^{\mathit{dexp}})$, by Horner's rule in a
floating-point format with a double mantissa and a separate exponent,
which cannot overflow. As for `fmpz_poly_evaluate_horner_d`, the error
is not controlled; cancellation can make the result meaningless for
polynomials with coefficients of both signs.
:::

## Multipoint evaluation and interpolation {#interpolation}

::: function
int _fmpz_poly_interpolate(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
int fmpz_poly_interpolate(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n)
int _fmpz_poly_interpolate_newton(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
int fmpz_poly_interpolate_newton(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n)
int _fmpz_poly_interpolate_multi_mod(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
int fmpz_poly_interpolate_multi_mod(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n)

Set `poly` to the polynomial $f$ of degree less than $n$ with $f(x_i)
= y_i$ for $0 \le i < n$ and return 1, if the $x_i$ are distinct and
$f$ has integer coefficients; otherwise return 0. The *newton* version
computes divided differences, failing at the first inexact division,
and converts from the Newton basis; the *multi_mod* version
interpolates modulo word-size primes (with a fast algorithm for many
points) and combines the results
by the Chinese remainder theorem until they stabilise and are checked.
The default uses Newton interpolation for $n \le 25$ and chooses by
size otherwise. The underscore versions write $n$ coefficients.
:::

::: function
void _fmpz_poly_interpolate_exact(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
void fmpz_poly_interpolate_exact(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n)
void _fmpz_poly_interpolate_exact_newton(fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
void fmpz_poly_interpolate_exact_newton(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n)

As `fmpz_poly_interpolate` and `fmpz_poly_interpolate_newton`, but
assuming that the $x_i$ are distinct and that the interpolating
polynomial has integer coefficients; this is not checked, and the
result is undefined otherwise. The default uses Newton interpolation
for $n \le 70$ and the multimodular algorithm otherwise.
:::

::: function
void fmpz_poly_interpolate_fmpz_vec(fmpz_poly_t poly, const fmpz * xs, const fmpz * ys, slong n)

As `fmpz_poly_interpolate`, but raises a `FLINT_INEXACT` exception
if the $x_i$ are not distinct or the interpolating polynomial does not
have integer coefficients.
:::

::: function
void _fmpz_poly_monomial_to_newton(fmpz * poly, const fmpz * roots, slong n)
void _fmpz_poly_newton_to_monomial(fmpz * poly, const fmpz * roots, slong n)

Convert `(poly, n)` in place between the monomial basis and the Newton
basis $1, (x - r_0), (x - r_0)(x - r_1), \ldots, (x - r_0) \cdots
(x - r_{n-2})$ for the integers $r_i$ in `roots`: `monomial_to_newton`
finds the $c_i$ with $f = \sum_i c_i \prod_{j < i} (x - r_j)$ by
repeated synthetic division, and `newton_to_monomial` evaluates such a
sum by Horner's rule.
:::

::: function
void _fmpz_poly_product_roots_fmpz_vec(fmpz * poly, const fmpz * xs, slong n)
void fmpz_poly_product_roots_fmpz_vec(fmpz_poly_t poly, const fmpz * xs, slong n)

Set `poly` to $(x - x_0)(x - x_1) \cdots (x - x_{n-1})$. The underscore
version writes $n + 1$ coefficients and does not allow aliasing.
:::

::: function
void _fmpz_poly_product_roots_fmpq_vec(fmpz * poly, const fmpq * xs, slong n)
void fmpz_poly_product_roots_fmpq_vec(fmpz_poly_t poly, const fmpq * xs, slong n)

Set `poly` to $(q_0 x - p_0) \cdots (q_{n-1} x - p_{n-1})$ for the
rational numbers $x_i = p_i / q_i$ in `xs`, the primitive polynomial
in $\Z[x]$ with roots $x_i$ when they are given in lowest terms. The
underscore version writes $n + 1$ coefficients.
:::

::: function
fmpz ** _fmpz_poly_tree_alloc(slong len)
void _fmpz_poly_tree_free(fmpz ** tree, slong len)

Allocate and free a subproduct tree for `len` linear factors. Level
$i$ of the tree is one array that holds the $\lfloor
\mathit{len}/2^i \rfloor$ products of $2^i$ linear factors, each
stored with $2^i + 1$ coefficients, followed by the product of the
remaining factors if $2^i$ does not divide `len`.
:::

::: function
void _fmpz_poly_tree_build_fmpq_vec(fmpz ** tree, const fmpq * roots, slong len)

Fills the tree allocated by `_fmpz_poly_tree_alloc` with the products
of the linear factors $q_i x - p_i$ for the roots $p_i / q_i$; the
product at the top level is not computed.
:::

## Composition and Taylor shift {#composition}

::: function
void _fmpz_poly_compose(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_compose(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Set `res` to the composition $g(h(x))$ of $g$ = `poly1` and $h$ =
`poly2`. Inner polynomials of the form $ax^n + b$ use a Taylor shift
and a scaling, short outer polynomials Horner's rule, and others the
divide-and-conquer algorithm. The underscore version requires
positive lengths, writes $(\mathit{len1} - 1)(\mathit{len2} - 1) + 1$
coefficients and does not allow aliasing. See [](fmpz_poly/compose).
:::

::: function
void _fmpz_poly_compose_horner(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_compose_horner(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)
void _fmpz_poly_compose_divconquer(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2)
void fmpz_poly_compose_divconquer(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2)

Composition by Horner's rule, respectively by the divide-and-conquer
algorithm that splits $g$ into blocks evaluated at $h$ and combines
them with the powers $h^{2^k}$. The underscore versions have the
conventions of `_fmpz_poly_compose`.
:::

::: function
void _fmpz_poly_taylor_shift(fmpz * poly, const fmpz_t c, slong n)
void fmpz_poly_taylor_shift(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c)

Set $g$ to $f(x + c)$; the underscore version shifts `(poly, n)` in
place. It uses Horner's scheme for length less than 64 and the
divide-and-conquer algorithm otherwise.
:::

::: function
void _fmpz_poly_taylor_shift_horner(fmpz * poly, const fmpz_t c, slong n)
void fmpz_poly_taylor_shift_horner(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c)

Taylor shift by Horner's scheme, $O(n^2)$ additions and
multiplications by $c$ (only additions for $c = \pm 1$).
:::

::: function
void _fmpz_poly_taylor_shift_divconquer(fmpz * poly, const fmpz_t c, slong n)
void fmpz_poly_taylor_shift_divconquer(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c)

Taylor shift by divide and conquer, $f_0(x + c) + (x + c)^m f_1(x +
c)$ with the powers of $x + c$ computed by the binomial theorem; the
two halves are shifted in parallel when threads are available. Short
inputs use Horner's scheme, with a cutoff that grows with the size of
the coefficients.
:::

::: function
void _fmpz_poly_taylor_shift_multi_mod(fmpz * poly, const fmpz_t c, slong n)
void fmpz_poly_taylor_shift_multi_mod(fmpz_poly_t g, const fmpz_poly_t f, const fmpz_t c)

Taylor shift by a multimodular algorithm: the shift is computed modulo
enough word-size primes for a bound on the result and combined by the
Chinese remainder theorem, using `flint_get_num_threads` threads.
:::

## Power series {#series}

These functions treat polynomials as power series truncated to length
$n$, that is, as elements of $\Z[[x]]/(x^n)$. See [](fmpz_poly/series)
and [](algorithms/polynomials/newton).

::: function
void _fmpz_poly_inv_series(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
void fmpz_poly_inv_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)

Set `Qinv` to the inverse of the power series $Q$ modulo $x^n$, which
requires the constant term of $Q$ to be $\pm 1$; this is not checked.
The basecase recurrence is used when $Q$ or $n$ is shorter than 64,
Newton iteration otherwise. A zero $Q$ or $n = 0$ raises an exception.
The underscore version requires $n \ge 1$, writes $n$ coefficients
and does not allow aliasing.
:::

::: function
void _fmpz_poly_inv_series_basecase(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
void fmpz_poly_inv_series_basecase(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)
void _fmpz_poly_inv_series_newton(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
void fmpz_poly_inv_series_newton(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)

Power series inversion by the recurrence $b_k = -b_0 \sum_{j=1}^{k}
q_j b_{k-j}$, with word-size arithmetic while the sizes allow and a
special case for binomials, respectively by Newton iteration. The
conventions are those of
`fmpz_poly_inv_series`.
:::

::: function
void _fmpz_poly_div_series(fmpz * Q, const fmpz * A, slong Alen, const fmpz * B, slong Blen, slong n)
void fmpz_poly_div_series(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B, slong n)

Set $Q$ to the quotient $A/B$ of power series modulo $x^n$. The
constant term of $B$ must be nonzero and the quotient must have
integer coefficients, which holds in particular if $B(0) = \pm 1$; an
inexact division raises an exception. For $B(0) = \pm 1$ and lengths
at least 72 the function multiplies by the inverse of $B$; otherwise
it uses the basecase recurrence or, for long inputs, divide and
conquer. A zero $B$ raises an exception. The underscore version
requires $n \ge 1$, writes $n$ coefficients and does not allow
aliasing.
:::

::: function
void _fmpz_poly_div_series_basecase(fmpz * Q, const fmpz * A, slong Alen, const fmpz * B, slong Blen, slong n)
void fmpz_poly_div_series_basecase(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B, slong n)
void _fmpz_poly_div_series_divconquer(fmpz * Q, const fmpz * A, slong Alen, const fmpz * B, slong Blen, slong n)
void fmpz_poly_div_series_divconquer(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B, slong n)

Power series division by the recurrence, respectively by
divide-and-conquer division of the reversed polynomials, with the
conventions of `fmpz_poly_div_series`.
:::

::: function
void _fmpz_poly_compose_series(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong n)
void fmpz_poly_compose_series(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)

Set `res` to $g(h(x)) \bmod x^n$ for $g$ = `poly1` and $h$ = `poly2`,
where $h(0)$ must be zero; otherwise an exception is raised. The
function uses Horner's scheme for $\mathit{len1} \le 10$, the
Brent–Kung baby-step giant-step algorithm [@BrentKung1978] for $n <
500$ and the algorithm of Kinoshita and Li [@KL2024] otherwise. The
underscore version requires $0 < \mathit{len1}, \mathit{len2} \le n$,
writes $n$ coefficients and does not allow aliasing; the
non-underscore version first reduces $n$ to at most
$(\mathit{len1} - 1)(\mathit{len2} - 1) + 1$.
:::

::: function
void _fmpz_poly_compose_series_horner(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong n)
void fmpz_poly_compose_series_horner(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)
void _fmpz_poly_compose_series_brent_kung(fmpz * res, const fmpz * poly1, slong len1, const fmpz * poly2, slong len2, slong n)
void fmpz_poly_compose_series_brent_kung(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)

Power series composition by Horner's scheme, respectively by
Algorithm 2.1 of [@BrentKung1978], with the conventions of
`fmpz_poly_compose_series`.
:::

::: function
void _fmpz_poly_revert_series(fmpz * Qinv, const fmpz * Q, slong Qlen, slong n)
void fmpz_poly_revert_series(fmpz_poly_t Qinv, const fmpz_poly_t Q, slong n)

Set `Qinv` to the compositional inverse of the power series $Q$
modulo $x^n$, so that $Q(Q^{-1}(x)) = Q^{-1}(Q(x)) = x \bmod x^n$.
Requires $Q(0) = 0$ and $Q'(0) = \pm 1$, which is not checked.
:::

::: function
int _fmpz_poly_sqrt_series(fmpz * res, const fmpz * poly, slong len, slong n)
int fmpz_poly_sqrt_series(fmpz_poly_t b, const fmpz_poly_t a, slong n)

If the power series $a$ of precision $n$ has a square root in
$\Z[[x]]$, set $b$ to it (the one with positive lowest coefficient)
and return 1; otherwise return 0. If $a$ has valuation $2v > 0$, the
square root is only determined to precision $n - v$, and the
underscore version pads with zeros. The underscore version writes $n$
coefficients.
:::

## Square roots {#sqrt}

::: function
int _fmpz_poly_sqrt(fmpz * res, const fmpz * poly, slong len)
int fmpz_poly_sqrt(fmpz_poly_t b, const fmpz_poly_t a)

If $a$ is the square of a polynomial in $\Z[x]$, set $b$ to the square
root with positive leading coefficient and return 1; otherwise return
0. The heuristic `_fmpz_poly_sqrt_KS` is tried first, and the
divide-and-conquer algorithm is used if it fails. The underscore
version writes `len / 2 + 1` coefficients.
:::

::: function
int _fmpz_poly_sqrt_classical(fmpz * res, const fmpz * poly, slong len, int exact)
int fmpz_poly_sqrt_classical(fmpz_poly_t b, const fmpz_poly_t a)
int _fmpz_poly_sqrt_divconquer(fmpz * res, const fmpz * poly, slong len, int exact)
int fmpz_poly_sqrt_divconquer(fmpz_poly_t b, const fmpz_poly_t a)

Square root with the conventions of `fmpz_poly_sqrt`, computed from
the top coefficient down, by the classical $O(n^2)$ algorithm,
respectively by divide and conquer. Most non-squares are detected
early by cheap tests: the length must be odd, the valuation even, the
leading and constant coefficients squares, and certain coefficients
even (the square root must exist modulo 2). If `exact` is
0, the underscore versions only compute the top half of the root and
allow a remainder, which is not computed.
:::

::: function
int _fmpz_poly_sqrt_KS(fmpz * rop, const fmpz * op, slong len)
int fmpz_poly_sqrt_KS(fmpz_poly_t b, const fmpz_poly_t a)

Heuristic square root by Kronecker substitution: after the same
quick tests as the classical algorithm, the integer square root of
$a(2^k)$, for $k$ a few bits more than the coefficients of $a$, is
unpacked and checked. Return $-1$ if the heuristic is inconclusive,
and otherwise 1 or 0 as `fmpz_poly_sqrt`.
:::

::: function
int _fmpz_poly_sqrtrem_classical(fmpz * res, fmpz * r, const fmpz * poly, slong len)
int fmpz_poly_sqrtrem_classical(fmpz_poly_t b, fmpz_poly_t r, const fmpz_poly_t a)
int _fmpz_poly_sqrtrem_divconquer(fmpz * res, fmpz * r, const fmpz * poly, slong len, fmpz * temp)
int fmpz_poly_sqrtrem_divconquer(fmpz_poly_t b, fmpz_poly_t r, const fmpz_poly_t a)

If $a$ can be written as $b^2 + r$ with $b \in \Z[x]$ and $\deg r <
\deg a / 2$ (equivalently, $\deg r < \deg b$), set $b$ and $r$ and return
1; otherwise return 0. The underscore versions write $\mathit{len}/2
+ 1$ coefficients of $b$ and require room for `len` coefficients in
$r$, which may be equal to `poly`; the divide-and-conquer version also
takes a temporary array of `len` coefficients.
:::

## Power sums {#power-sums}

For a monic $f = \prod_{i=1}^d (x - \alpha_i)$, the power sums are
$p_k = \sum_i \alpha_i^k$, and the power sum series is $\sum_{k \ge 0}
p_k x^k = d - x f_{\mathrm{rev}}'(x) / f_{\mathrm{rev}}(x)$, where
$f_{\mathrm{rev}} = x^d f(1/x)$. These sums determine $f$ by Newton's
identities.

::: function
void fmpz_poly_power_sums(fmpz_poly_t res, const fmpz_poly_t poly, slong n)
void _fmpz_poly_power_sums_naive(fmpz * res, const fmpz * poly, slong len, slong n)
void fmpz_poly_power_sums_naive(fmpz_poly_t res, const fmpz_poly_t poly, slong n)

Set `res` to the power sums $p_0, \ldots, p_{n-1}$ of the roots of the
monic polynomial `poly`, computed with Newton's identities; the roots
at zero are counted in $p_0$ only. A zero polynomial raises an
exception in `fmpz_poly_power_sums`. The underscore version writes $n$
coefficients.
:::

::: function
void _fmpz_poly_power_sums_to_poly(fmpz * res, const fmpz * poly, slong len)
void fmpz_poly_power_sums_to_poly(fmpz_poly_t res, const fmpz_poly_t Q)

Set `res` to the monic polynomial whose roots have the power sums
given by $Q$, with $p_0 = Q(0)$ the degree, using Newton's
identities; the power sums must be those of a polynomial in $\Z[x]$.
:::

## Integer and rational roots {#roots}

::: function
slong fmpz_poly_roots_fmpz(fmpz * res, slong * exp, const fmpz_poly_t poly)
slong fmpz_poly_roots_fmpq(fmpq * res, slong * exp, const fmpz_poly_t poly)

Compute the distinct integer (respectively rational) roots of the
nonzero polynomial `poly`, write them to `res` in increasing order and
return their number. If `exp` is not `NULL`, the multiplicities are
written to `exp`. The arrays must have room for $\deg(\mathit{poly})$
entries. A zero polynomial raises an exception.

::: details
After removing roots at zero, if the constant and leading coefficients
are small with few divisors, the candidates $a/b$ with $a$ dividing the
constant and $b$ the leading coefficient are tested directly.
Otherwise the roots are found modulo a suitable prime $p$, lifted
$p$-adically by Newton iteration to the precision needed to recover
rational roots of bounded height, and verified. This is much faster
than factoring the polynomial. Multiplicities are found with a
squarefree factorisation, computed only if the polynomial is not
squarefree.
:::
:::

::: function
slong _fmpz_poly_roots_fmpz_squarefree(fmpz * res, const fmpz * poly, slong len)
slong _fmpz_poly_roots_fmpq_squarefree(fmpq * res, const fmpz * poly, slong len)

Write the integer (respectively rational) roots of `(poly, len)` to
`res` in no particular order and return their number. The polynomial
must be squarefree, with $\mathit{len} \ge 1$ and a nonzero constant
coefficient; `res` must have room for $\mathit{len} - 1$ entries.
:::

::: function
void _fmpz_poly_bound_roots(fmpz_t bound, const fmpz * poly, slong len)
void fmpz_poly_bound_roots(fmpz_t bound, const fmpz_poly_t poly)

Set `bound` to a nonnegative integer that is intended to bound the
absolute values of all complex roots, computed from Fujiwara's bound
[@Fuj1916]

$$
2 \max \left( \left|\frac{a_{n-1}}{a_n}\right|,
\left|\frac{a_{n-2}}{a_n}\right|^{1/2}, \ldots,
\left|\frac{a_1}{a_n}\right|^{1/(n-1)},
\left|\frac{a_0}{2 a_n}\right|^{1/n} \right)
$$

for $f = \sum_{i=0}^n a_i x^i$, with each term rounded up to an integer
(for $n \ge 2$, as $\lfloor \cdot \rfloor + 1$). The bound is zero for
constant polynomials.
:::

## Real roots {#real-roots}

Unless stated otherwise, the functions in this section assume that the
polynomial is squarefree; the results are undefined otherwise, and
some functions raise an exception. Roots at zero are counted where
noted. See [](algorithms/polynomials/real-roots-integers). For
numerical enclosures of the real or complex roots, see
[](arb_fmpz_poly/roots).

::: function
slong _fmpz_poly_num_real_roots(const fmpz * pol, slong len)
slong fmpz_poly_num_real_roots(const fmpz_poly_t poly)

Return the number of real roots of the squarefree polynomial `poly`,
including a root at zero. A zero polynomial raises an exception.

::: details
Degrees up to 4 use closed formulas based on the discriminant.
Otherwise Descartes' rule of signs is tried first; it decides the
answer when $f(x)$ and $f(-x)$ each have at most one sign variation.
Short polynomials with small coefficients use the subresultant Sturm
sequence; larger ones first try a primitive Sturm sequence
(`_fmpz_poly_num_real_roots_sturm_bounded`) with a bound that detects
coefficient growth, and fall back on the Vincent–Collins–Akritas
algorithm.
:::
:::

::: function
void _fmpz_poly_num_real_roots_sturm(slong * n_neg, slong * n_pos, const fmpz * pol, slong len)
slong fmpz_poly_num_real_roots_sturm(const fmpz_poly_t pol)

Count real roots with a Sturm sequence computed as the subresultant
pseudo-remainder sequence (with `_fmpz_poly_pseudo_rem_cohen`). The
underscore version sets the numbers of negative and positive roots and
requires a squarefree polynomial of degree at least 2 with nonzero
constant coefficient.
:::

::: function
slong _fmpz_poly_num_real_roots_sturm_bounded(const fmpz * pol, slong len, int on_0_1, slong max_size)

Counts the real roots (all of them if `on_0_1` is zero, those in $(0,
1)$ otherwise) with a Sturm sequence computed as a primitive
pseudo-remainder sequence, removing the content after each step.
Returns $-1$ as soon as a remainder has length times coefficient bit
size larger than `max_size`; pass `WORD_MAX` to always finish.

::: details
The primitive sequence is much cheaper than the subresultant sequence
for polynomials with special structure, such as classical orthogonal
polynomials or sparse polynomials whose degrees drop quickly, and
usually more expensive for generic polynomials, whose coefficients grow
along the sequence; `max_size` detects the latter case cheaply.
:::
:::

::: function
slong _fmpz_poly_num_real_roots_vca(const fmpz * pol, slong len)
slong fmpz_poly_num_real_roots_vca(const fmpz_poly_t pol)

Count real roots with the Vincent–Collins–Akritas algorithm, by
isolating them with `_fmpz_poly_isolate_real_roots_0_1_vca` after
scaling the positive and negative roots into $(0, 1)$.
:::

::: function
slong fmpz_poly_num_real_roots_0_1(const fmpz_poly_t pol)
slong fmpz_poly_num_real_roots_0_1_sturm(const fmpz_poly_t pol)
slong fmpz_poly_num_real_roots_0_1_vca(const fmpz_poly_t pol)

Return the number of roots of the squarefree polynomial `pol` in the
open interval $(0, 1)$. The default first applies Descartes' rule of
signs on $(0, 1)$, which settles most generic polynomials, then tries
the primitive Sturm sequence with a growth bound and falls back on the
Vincent–Collins–Akritas algorithm. A zero polynomial raises an
exception.
:::

::: function
slong fmpz_poly_num_real_roots_upper_bound(const fmpz_poly_t pol)

Returns an upper bound for the number of real roots of `pol`, counted
with multiplicity: the multiplicity of the root at zero plus the
Descartes bounds for the positive and the negative roots, each
computed after scaling the roots into $(0, 1)$. The polynomial need
not be squarefree. A zero polynomial raises an exception.
:::

::: function
slong _fmpz_poly_descartes_bound_0_1(const fmpz * p, slong len, slong bound)

Returns the number of sign variations of $(1 + x)^{n} p(1/(1 + x))$,
which bounds the number of roots of `(p, len)` in $(0, 1)$ and has the
same parity (Descartes' rule of signs), or `WORD_MAX` as soon as the
count exceeds `bound`.
:::

::: function
int _fmpz_poly_has_real_root(const fmpz * p, slong len)
int fmpz_poly_has_real_root(const fmpz_poly_t pol)

Return whether the squarefree polynomial has a real root. Odd degree, a
root at zero, a sign change between $0$ and $\infty$, or an odd number
of sign variations of $f(x)$ or $f(-x)$ answer immediately; otherwise
the real roots are counted.
:::

::: function
slong _fmpz_poly_positive_root_upper_bound_2exp_local_max(const fmpz * pol, slong len)
slong _fmpz_poly_positive_root_upper_bound_2exp(const fmpz * pol, slong len)
slong fmpz_poly_positive_root_upper_bound_2exp(const fmpz_poly_t pol)

Return $e$ such that $2^e$ is an upper bound (not necessarily strict)
for the positive real roots of `pol`, computed with the local-max
bound of Akritas, Strzeboński and Vigklas [@ASV2008] as implemented in
SLV by Tsigaridas; `WORD_MIN` is returned if the bound shows that
there are no positive roots. The non-underscore version returns 0 for
the zero polynomial and ignores roots at zero.
:::

::: function
void _fmpz_poly_isolate_real_roots_0_1_vca(fmpq * exact_roots, slong * n_exact_roots, fmpz * c_array, slong * k_array, slong * n_intervals, const fmpz * pol, slong len)

Isolates the roots of `(pol, len)` in the open interval $(0, 1)$ with
the Vincent–Collins–Akritas bisection algorithm. Roots found exactly
(dyadic numbers) are appended to `exact_roots`, incrementing
`*n_exact_roots`; every other root is enclosed in an interval
$(c\,2^k, (c + 1) 2^k)$ containing exactly one root, appended as $c$ =
`c_array[i]` and $k$ = `k_array[i]`, incrementing `*n_intervals`. If
`exact_roots`, respectively `c_array` or `k_array`, is `NULL`, only the
counters are updated.
:::

::: function
void fmpz_poly_isolate_real_roots(fmpq * exact_roots, slong * n_exact, fmpz * c_array, slong * k_array, slong * n_interval, const fmpz_poly_t pol)
void fmpz_poly_isolate_positive_roots(fmpq * exact_roots, slong * n_exact, fmpz * c_array, slong * k_array, slong * n_interval, const fmpz_poly_t pol)

Isolate the real (respectively strictly positive) roots of the
squarefree polynomial `pol`. Exact dyadic roots are written to
`exact_roots` and counted in `*n_exact`; the other roots are enclosed
in disjoint open intervals $(c\,2^k, (c + 1) 2^k)$ with $c$ =
`c_array[i]` and $k$ = `k_array[i]`, counted in `*n_interval`. The
arrays must have room for $\deg(\mathit{pol})$ entries.

::: details
Three algorithms are used. For polynomials of large degree, bisection
with a Sturm sequence is tried first; the sequence is computed as a
primitive remainder sequence and evaluated at each bisection point in
$O(n)$ operations using the pseudo-quotients. This is much faster than
VCA when the Sturm sequence has small coefficients and there are many
real roots (classical orthogonal polynomials) or clusters of close
roots (sparse polynomials). Otherwise, if the polynomial appears to
have only real roots (Newton's inequalities hold and the numbers of
sign variations of $f(x)$ and $f(-x)$ add up to the degree) spread
over many orders of magnitude, as for Eulerian polynomials, sign
changes are located by evaluating at dyadic points in ball arithmetic,
and the result is certified by checking that the number of sign
changes equals the number of sign variations. This avoids the huge
coefficients that VCA produces for such polynomials. If these methods
do not apply or fail, the Vincent–Collins–Akritas algorithm is used.
:::
:::

::: function
void _fmpz_poly_signature(slong * r1, slong * r2, const fmpz * poly, slong len)
void fmpz_poly_signature(slong * r1, slong * r2, const fmpz_poly_t poly)

Set $(r_1, r_2)$ to the signature of the squarefree polynomial `poly`:
it has $r_1$ real roots and $r_2$ pairs of complex conjugate roots.
Uses a Sturm sequence, Algorithm 4.1.11 of [@Coh1996]. The zero
polynomial gives $(0, 0)$; for a polynomial that is not squarefree the
result is undefined and an exception may be raised.
:::

## Modular reduction and reconstruction {#crt}

::: function
void fmpz_poly_get_nmod_poly(nmod_poly_t res, const fmpz_poly_t poly)

Sets `res` to `poly` reduced modulo the modulus of `res`.
:::

::: function
void fmpz_poly_set_nmod_poly(fmpz_poly_t res, const nmod_poly_t poly)
void fmpz_poly_set_nmod_poly_unsigned(fmpz_poly_t res, const nmod_poly_t poly)

Set `res` to the lift of `poly` with coefficients in the symmetric
range $[-m/2, m/2)$, respectively in $[0, m)$, where $m$ is the
modulus of `poly`.
:::

::: function
void _fmpz_poly_CRT_ui_precomp(fmpz * res, const fmpz * poly1, slong len1, const fmpz_t m1, nn_srcptr poly2, slong len2, ulong m2, ulong m2inv, fmpz_t m1m2, ulong c, int sign)
void _fmpz_poly_CRT_ui(fmpz * res, const fmpz * poly1, slong len1, const fmpz_t m1, nn_srcptr poly2, slong len2, ulong m2, ulong m2inv, int sign)

Set `res` to the polynomial congruent to `(poly1, len1)` modulo $m_1$
and to `(poly2, len2)` modulo the word-size $m_2$, with coefficients in
$[0, m_1 m_2)$ if `sign` = 0 and in $[-m_1 m_2/2, m_1 m_2/2)$ if
`sign` = 1, writing $\max(\mathit{len1}, \mathit{len2})$ coefficients.
The moduli must be coprime. `m2inv` is `n_preinvert_limb(m2)`; the
`precomp` version also takes $m_1 m_2$ and $c = m_1^{-1} \bmod m_2$,
which the other version computes, raising an exception if $m_1$ is not
invertible modulo $m_2$. `res` may be equal to `poly1`.
:::

::: function
void fmpz_poly_CRT_ui(fmpz_poly_t res, const fmpz_poly_t poly1, const fmpz_t m1, const nmod_poly_t poly2, int sign)

Sets `res` to the CRT combination of `poly1` modulo $m_1$ and `poly2`
modulo its modulus $m_2$, coprime to $m_1$, with coefficients in $[0,
m_1 m_2)$ if `sign` = 0 and in $[-m_1 m_2/2, m_1 m_2/2)$ if `sign` = 1.
:::

## Hensel lifting {#hensel}

These functions lift a factorisation $f \equiv \ell\, g_1 \cdots g_r
\pmod p$ into pairwise coprime monic factors to a factorisation modulo
$p^N$, as needed for factoring over $\Z$; see
[](algorithms/polynomials/hensel-lifting) for the algorithm and
`fmpz_poly_hensel_lift_once` for the simplest entry point. The factors
are lifted along a binary tree whose nodes hold products of factors
and Bézout cofactors, with quadratic lifting $p^e \to p^{2e}$.

::: function
void fmpz_poly_hensel_lift_once(fmpz_poly_factor_t lifted_fac, const fmpz_poly_t f, const nmod_poly_factor_t local_fac, slong N)

Lifts the $r \ge 2$ factors in `local_fac` of $f$ modulo a prime $p$
(the modulus of the factors) to factors modulo $p^N$, stored in the
same order in `lifted_fac`. The factors must be monic, pairwise
coprime, of product $f / \operatorname{lc}(f)$ modulo $p$, and
$\operatorname{lc}(f)$ must not be divisible by $p$; the lifted
factors are monic with product $f / \operatorname{lc}(f)$ modulo
$p^N$. The lift cannot be continued afterwards.
:::

::: function
slong _fmpz_poly_hensel_start_lift(fmpz_poly_factor_t lifted_fac, slong * link, fmpz_poly_t * v, fmpz_poly_t * w, const fmpz_poly_t f, const nmod_poly_factor_t local_fac, slong target_exp)

Builds the Hensel tree for the $r \ge 2$ local factors in `local_fac`
of $f$ (with the conditions of `fmpz_poly_hensel_lift_once`) and lifts
them to precision $p^N$ with $N$ = `target_exp` $\ge 1$, storing the
lifted factors in `lifted_fac` in the original order. `link`, `v` and
`w` must be arrays of $2r - 2$ entries (initialised polynomials for
`v`, `w`), which hold the tree on return so that the lift can be
continued. Returns the exponent $e$ to which the cofactors are known,
to be passed as `prev` to `_fmpz_poly_hensel_continue_lift`.
:::

::: function
slong _fmpz_poly_hensel_continue_lift(fmpz_poly_factor_t lifted_fac, slong * link, fmpz_poly_t * v, fmpz_poly_t * w, const fmpz_poly_t f, slong prev, slong curr, slong N, const fmpz_t p)

Continues a lift started by `_fmpz_poly_hensel_start_lift`: the
factors in the tree `(link, v, w)` are known modulo $p^{\mathit{curr}}$
and the cofactors modulo $p^{\mathit{prev}}$, as returned by the
previous call; the factors are lifted to $p^N$ for $N >
\mathit{curr} \ge 1$ and written to `lifted_fac`. Returns the new
value of `prev`.
:::

::: function
void fmpz_poly_hensel_build_tree(slong * link, fmpz_poly_t * v, fmpz_poly_t * w, const nmod_poly_factor_t fac)

Builds a Hensel tree from the $r$ local factors in `fac`. The arrays
`link`, `v` and `w` have $2r - 2$ entries. The tree is built by
repeatedly multiplying the two factors (or products) of smallest
degree, so that it is balanced by degree; `v` holds the factors and
products, with `v[2r - 4]` and `v[2r - 3]` the two children of the
root, and `w` holds the Bézout cofactors, $v_{2i} w_{2i} + v_{2i+1}
w_{2i+1} \equiv 1 \pmod p$. A negative `link[j]` means that `v[j]` is
the original factor number `-link[j] - 1`; a nonnegative `link[j]` means
that `v[j]` is the product of `v[link[j]]` and `v[link[j] + 1]`.
:::

::: function
void fmpz_poly_hensel_lift(fmpz_poly_t Gout, fmpz_poly_t Hout, fmpz_poly_t Aout, fmpz_poly_t Bout, const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_poly_t h, const fmpz_poly_t a, const fmpz_poly_t b, const fmpz_t p, const fmpz_t p1)

One Hensel step for a factorisation into two factors: given $f \equiv
gh$ and $ag + bh \equiv 1$ modulo $p$, with $g$ and $h$ monic, computes
$G \equiv g$, $H \equiv h$, $A \equiv a$, $B \equiv b$ modulo $p$ with
$f \equiv GH$ and $AG + BH \equiv 1$ modulo $P = p p_1$:

$$
\begin{aligned}
G &= g + p \left( \tfrac{f - gh}{p}\, b \bmod g \right), &
H &= h + p \left( \tfrac{f - gh}{p}\, a \bmod h \right), \\
B &= b + p \left( \tfrac{1 - aG - bH}{p}\, b \bmod g \right), &
A &= a + p \left( \tfrac{1 - aG - bH}{p}\, a \bmod h \right),
\end{aligned}
$$

with the products reduced modulo $p_1$. Requires $1 < p_1 \le p$
(so $p_1 \mid p$ for prime powers of the same prime), $f, g, h$ of
degree at least 1, and $a, b$ nonzero. The outputs $G, H, A, B$ may be
aliased only with $g, h, a, b$ respectively.
:::

::: function
void _fmpz_poly_hensel_lift_without_inverse(fmpz * G, fmpz * H, const fmpz * f, slong lenF, const fmpz * g, slong lenG, const fmpz * h, slong lenH, const fmpz * a, slong lenA, const fmpz * b, slong lenB, const fmpz_t p, const fmpz_t p1)
void fmpz_poly_hensel_lift_without_inverse(fmpz_poly_t Gout, fmpz_poly_t Hout, const fmpz_poly_t f, const fmpz_poly_t g, const fmpz_poly_t h, const fmpz_poly_t a, const fmpz_poly_t b, const fmpz_t p, const fmpz_t p1)

As `fmpz_poly_hensel_lift`, lifting only the factors $g$ and $h$.
:::

::: function
void _fmpz_poly_hensel_lift_only_inverse(fmpz * A, fmpz * B, const fmpz * G, slong lenG, const fmpz * H, slong lenH, const fmpz * a, slong lenA, const fmpz * b, slong lenB, const fmpz_t p, const fmpz_t p1)
void fmpz_poly_hensel_lift_only_inverse(fmpz_poly_t Aout, fmpz_poly_t Bout, const fmpz_poly_t G, const fmpz_poly_t H, const fmpz_poly_t a, const fmpz_poly_t b, const fmpz_t p, const fmpz_t p1)

As `fmpz_poly_hensel_lift`, lifting only the cofactors, given the
already lifted factors $G$ and $H$.
:::

::: function
void fmpz_poly_hensel_lift_tree_recursive(slong * link, fmpz_poly_t * v, fmpz_poly_t * w, fmpz_poly_t f, slong j, slong inv, const fmpz_t p0, const fmpz_t p1)

Lifts the subtree of the Hensel tree `(link, v, w)` rooted at the pair
`v[j]`, `v[j + 1]`, whose product is $f$, from modulus $p_0$ to $P =
p_0 p_1$ with $1 < p_1 \le p_0$. The cofactors are lifted if `inv` is
1, and not lifted if it is 0 (for the last step); with $-1$ the
cofactors are lifted from the previous precision (when a lift is
restarted). $f$ may be aliased with an entry of `v` or `w`, but the
entries of `v` and `w` must be distinct.

::: details
At each node the two halves of the lift are fused: since $p_1 \mid
p_0$, the moduli $g \bmod p_1$ and $h \bmod p_1$ and their Newton
preinverses serve all four divisions, and when $p \nmid f(0)$ the new
second factor is computed as an exact quotient $H = f/G$ modulo $P$
instead of a second lift. The arithmetic runs in `gr` rings chosen by
the size of the modulus (`nmod`, `mpn_mod` or `fmpz_mod`). Small nodes
use `fmpz_poly_hensel_lift`.
:::
:::

::: function
void fmpz_poly_hensel_lift_tree(slong * link, fmpz_poly_t * v, fmpz_poly_t * w, fmpz_poly_t f, slong r, const fmpz_t p, slong e0, slong e1, slong inv)

Lifts the whole Hensel tree for the $r$ factors of the monic $f$ from
$p^{e_0}$ to $p^{e_1}$, where $e_0 < e_1 \le 2 e_0$, by calling
`fmpz_poly_hensel_lift_tree_recursive` at the root with $p_0 =
p^{e_0}$ and $p_1 = p^{e_1 - e_0}$. `inv` is as for that function.
:::

## Special polynomials {#special}

::: function
void _fmpz_poly_cyclotomic(fmpz * a, ulong n, nn_ptr factors, slong num_factors, ulong phi)
void fmpz_poly_cyclotomic(fmpz_poly_t poly, ulong n)

Set `poly` to the cyclotomic polynomial $\Phi_n(x) = \prod_\omega (x -
\omega)$, where $\omega$ runs over the primitive $n$-th roots of unity.
The non-underscore version writes $n = qs$ with $q$ squarefree and
uses $\Phi_n(x) = \Phi_q(x^s)$.

::: details
The underscore version requires $n \ge 3$ squarefree, the array
`factors` of its `num_factors` prime factors (with 2 first if $n$ is
even) and $\varphi(n)$ as `phi`, and writes only the low $\lfloor
(\varphi(n) + 1)/2 \rfloor$ coefficients; the others follow by
symmetry, since $\Phi_n$ is palindromic for $n \ne 1$. It uses the
sparse power series algorithm (Algorithm 4) of [@ArnoldMonagan2011],
based on $\Phi_n(x) = \prod_{d \mid n} (x^d - 1)^{\mu(n/d)}$, in which
multiplications and divisions by $x^d - 1$ cost $O(\varphi(n))$
additions, for $O(2^k \varphi(n))$ operations in total for $k$ prime
factors. Machine integers are used when no overflow is possible: by
Table 6 of [@ArnoldMonagan2011], all coefficients of $\Phi_n$ have at
most 27 bits for $n < 10163195$ and at most 60 bits for $n <
169828113$, and they are $0, \pm 1$ when $n$ has two prime factors.
The cases $n = p$ and $n = 2m$ with $m$ odd ($\Phi_n(x) = \Phi_m(-x)$)
are handled directly.
:::
:::

::: function
ulong _fmpz_poly_is_cyclotomic(const fmpz * poly, slong len)
ulong fmpz_poly_is_cyclotomic(const fmpz_poly_t poly)

If `poly` is a cyclotomic polynomial $\Phi_n$, return $n$; otherwise
return 0.
:::

::: function
void _fmpz_poly_cos_minpoly(fmpz * f, ulong n)
void fmpz_poly_cos_minpoly(fmpz_poly_t f, ulong n)

Set $f$ to the minimal polynomial of $2 \cos(2\pi/n)$, which is monic
with integer coefficients [@WaktinsZeitlin1993]; substituting $x \to
2x$ gives the minimal polynomial of $\cos(2\pi/n)$ up to a constant.
Its degree is $\varphi(n)/2$ for $n > 2$ and 1 for $n = 1, 2$; for $n
= 0$ the result is defined to be 1. Suitable $n$ give the minimal
polynomials of $2 \cos(a\pi)$ and $2 \sin(a\pi)$ for all rational $a$.
:::

::: function
void _fmpz_poly_swinnerton_dyer(fmpz * T, ulong n)
void fmpz_poly_swinnerton_dyer(fmpz_poly_t poly, ulong n)

Set `poly` to the Swinnerton-Dyer polynomial $S_n = \prod (x \pm
\sqrt{2} \pm \sqrt{3} \pm \sqrt{5} \pm \cdots \pm \sqrt{p_n})$, where
$p_n$ is the $n$-th prime and all sign combinations are taken. It has
degree $2^n$ and is irreducible over $\Z$, but splits into factors of
degree at most 2 modulo every prime, which makes it the standard worst
case for the Zassenhaus algorithm.
:::

::: function
void _fmpz_poly_chebyshev_t(fmpz * coeffs, ulong n)
void fmpz_poly_chebyshev_t(fmpz_poly_t poly, ulong n)
void _fmpz_poly_chebyshev_u(fmpz * coeffs, ulong n)
void fmpz_poly_chebyshev_u(fmpz_poly_t poly, ulong n)

Set `poly` to the Chebyshev polynomial of the first kind $T_n$, with
$T_n(\cos \theta) = \cos(n\theta)$, respectively of the second kind
$U_n = T'_{n+1}/(n + 1)$, for $n \ge 0$. The coefficients are computed
by their hypergeometric recurrence. The underscore versions write
$n + 1$ coefficients.
:::

::: function
void _fmpz_poly_legendre_pt(fmpz * coeffs, ulong n)
void fmpz_poly_legendre_pt(fmpz_poly_t poly, ulong n)

Set `poly` to the shifted Legendre polynomial $\tilde{P}_n(x) =
P_n(2x - 1)$, which has integer coefficients, computed by a
hypergeometric recurrence. The Legendre polynomials themselves are in
`fmpq_poly`. The underscore version writes $n + 1$ coefficients.
:::

::: function
void _fmpz_poly_hermite_h(fmpz * coeffs, ulong n)
void fmpz_poly_hermite_h(fmpz_poly_t poly, ulong n)
void _fmpz_poly_hermite_he(fmpz * coeffs, ulong n)
void fmpz_poly_hermite_he(fmpz_poly_t poly, ulong n)

Set `poly` to the physicists' Hermite polynomial $H_n$ (with $H_n' =
2n H_{n-1}$), respectively the probabilists' Hermite polynomial
$\mathit{He}_n(x) = 2^{-n/2} H_n(x / \sqrt{2})$, computed by a
hypergeometric recurrence. The underscore versions write $n + 1$
coefficients.
:::

::: function
void _fmpz_poly_fibonacci(fmpz * coeffs, ulong n)
void fmpz_poly_fibonacci(fmpz_poly_t poly, ulong n)

Set `poly` to the $n$-th Fibonacci polynomial $F_n$, with $F_0 = 0$,
$F_1 = 1$ and $F_{n+1} = x F_n + F_{n-1}$, computed by a hypergeometric
recurrence.
:::

::: function
void fmpz_poly_eulerian_polynomial(fmpz_poly_t poly, ulong n)

Sets `poly` to the Eulerian polynomial $A_n(x) = \sum_{m=0}^n A(n, m)
x^m$, with $A_0 = 1$, computed by a recurrence. The coefficients are
the Eulerian numbers

$$
A(n, m) = \sum_{k=0}^{m+1} (-1)^k \binom{n+1}{k} (m + 1 - k)^n,
$$

and the polynomials have the exponential generating function $(x -
1)/(x - e^{(x - 1)t}) = \sum_n A_n(x) t^n / n!$. They are not the
Euler polynomials.
:::

::: function
void _fmpz_poly_eta_qexp(fmpz * f, slong e, slong n)
void fmpz_poly_eta_qexp(fmpz_poly_t f, slong e, slong n)

Set $f$ to the $q$-expansion to length $n$ of $(q^{-1/24} \eta(q))^e =
\prod_{k \ge 1} (1 - q^k)^e$, the Dedekind eta function without its
leading factor, raised to the power $e$. For $e = -1$ this is the
generating function of the partition function, and for $e = 24$,
after multiplication by $q$, the modular discriminant $\Delta$, which
generates the Ramanujan $\tau$ function. Sparse formulas are used for
$e = 1, 2, 3, 4, 6$, and power series arithmetic reduces other $e$ to
those.
:::

::: function
void _fmpz_poly_theta_qexp(fmpz * f, slong e, slong n)
void fmpz_poly_theta_qexp(fmpz_poly_t f, slong e, slong n)

Set $f$ to the $q$-expansion to length $n$ of $\vartheta(q)^e$, where
$\vartheta(q) = 1 + 2 \sum_{k \ge 1} q^{k^2}$ is the Jacobi theta
function. Sparse formulas are used for $e = 1, 2$, and power series
arithmetic reduces other $e$ to those.
:::

## Internal functions {#internal}

::: function
void fmpz_poly_CLD_bound(fmpz_t res, const fmpz_poly_t f, slong n)

Sets `res` to a bound for the absolute value of the coefficient of
$x^n$ in $f g'/g$ for every factor $g$ of $f$ in $\Z[x]$ (the
coefficients of the logarithmic derivative used by van Hoeij's
algorithm, [](algorithms/polynomials/factoring-integers#cld)).
Requires $0 \le n < \deg f$.

::: details
With $f = \sum_{i=0}^N a_i x^i$, $B_1(r) = r^{-(n+1)} \sum_{i \le n}
|a_i| r^i$ and $B_2(r) = r^{-(n+1)} \sum_{i > n} |a_i| r^i$, the bound
is $N \max(B_1(r), B_2(r))$ for an $r > 0$ that approximately minimises
the maximum, found by a search on $\log_2 r$ in double precision with
an extra safety margin [@HvHN2011].
:::
:::

::: function
slong _fmpz_poly_remove_content_2exp(fmpz * pol, slong len)

Divides `(pol, len)` by $2^k$, the largest power of two dividing all
coefficients, and returns $k$; returns 0 for the zero polynomial.
:::

::: function
void _fmpz_poly_scale_2exp(fmpz * pol, slong len, slong k)

Replaces `(pol, len)` by $f(2^k x)$ for $k \ge 0$, or by $2^{-k
(\mathit{len} - 1)} f(2^k x)$ for $k < 0$ (so that the result has
integer coefficients), and then divides by the largest power of two
dividing all coefficients.
:::
