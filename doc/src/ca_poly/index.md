---
title: ca_poly.h
header: ca_poly.h
description: Dense univariate polynomials with exact real and complex coefficients: arithmetic, division, GCDs, power series, squarefree factorisation and roots.
category: Exact real and complex numbers
order: 4
---

This module implements dense univariate polynomials whose coefficients
are exact real or complex numbers of type `ca_t` (see [](ca)), that is,
elements of $\mathbb{C}[x]$ with coefficients in the subfield of
$\mathbb{C}$ that `ca_t` can represent. Besides arithmetic and power
series operations, it computes GCDs, squarefree factorisations and,
where the roots can be expressed exactly, the roots of a polynomial
with multiplicities.

Many operations are thin wrappers around the generic polynomial
functions of [](gr_poly) with a `ca` context; others have dedicated
code that recognises coefficients lying in $\mathbb{Q}$ or in a common
number field and then uses fast arithmetic over $\Z$ ([](fmpz_poly),
[](fmpq_poly)). Because coefficients are tested for zero with
`ca_check_is_zero`, which may fail to decide, every operation that
depends on a degree or on the invertibility of a leading coefficient
can fail; such functions return an `int` success flag.

## Types and conventions {#conventions}

::: type
ca_poly_struct
ca_poly_t

A polynomial with `ca_t` coefficients. The structure has the fields
`coeffs` (an array of `alloc` initialised coefficients), `length` and
`alloc`; the coefficient of $x^i$ is `coeffs[i]`. A `ca_poly_t` is an
array of one structure, so it is passed by reference.
:::

::: function
ca_ptr ca_poly_coeff_ptr(ca_poly_t poly, slong i)

Returns a pointer to the coefficient of $x^i$. The index must be less
than the length (or the allocation, when writing coefficients before
setting the length); this is not checked.
:::

- *Normalisation.* A polynomial is normalised by removing top
  coefficients that are proved to be zero. A top coefficient that
  cannot be proved zero is kept, so the length can exceed the degree of
  the mathematical polynomial plus one. When the exact degree matters,
  check the leading coefficient, for example with `ca_poly_is_proper`.
- *Special values.* *Undefined* and infinities are not meaningful as
  coefficients. The functions in this module do not produce them from
  numerical input unless stated otherwise, and assume that the user
  does not store them. *Unknown* can occur as a coefficient to signal a
  number that could not be computed.
- A polynomial is *proper* if all its coefficients are numbers and its
  leading coefficient is proved nonzero (or the polynomial is zero).
- *Aliasing.* Functions without a leading underscore manage memory and
  allow outputs to be the same objects as inputs. Functions with a
  leading underscore work on coefficient arrays `(poly, len)` and leave
  memory management to the caller; unless stated otherwise they require
  positive lengths and do not allow the output to overlap the inputs.

## Memory management {#memory}

::: function
void ca_poly_init(ca_poly_t poly, ca_ctx_t ctx)
void ca_poly_init2(ca_poly_t poly, slong len, ca_ctx_t ctx)

Initialises `poly` as the zero polynomial, with space for `len`
coefficients for `ca_poly_init2`.
:::

::: function
void ca_poly_clear(ca_poly_t poly, ca_ctx_t ctx)

Clears `poly`, freeing the coefficients and the coefficient array.
:::

::: function
void ca_poly_fit_length(ca_poly_t poly, slong len, ca_ctx_t ctx)

Ensures that `poly` has at least `len` initialised coefficients
allocated, growing the array to at least twice its size if it must
grow. The length is not changed.
:::

::: function
void _ca_poly_set_length(ca_poly_t poly, slong len, ca_ctx_t ctx)

Sets the length of `poly` to `len`, which must not exceed the
allocation. Coefficients beyond the new length are set to zero; nothing
is normalised.
:::

::: function
void _ca_poly_normalise(ca_poly_t poly, ca_ctx_t ctx)

Decreases the length of `poly` while the top coefficient is proved to be
zero by `ca_check_is_zero`.
:::

::: function
void ca_poly_swap(ca_poly_t poly1, ca_poly_t poly2, ca_ctx_t ctx)

Swaps `poly1` and `poly2` in constant time.
:::

## Assignment and simple values {#assignment}

::: function
void ca_poly_zero(ca_poly_t poly, ca_ctx_t ctx)
void ca_poly_one(ca_poly_t poly, ca_ctx_t ctx)
void ca_poly_x(ca_poly_t poly, ca_ctx_t ctx)

Sets `poly` to the zero polynomial, the constant 1, or the monomial $x$.
:::

::: function
void ca_poly_set_ca(ca_poly_t poly, const ca_t c, ca_ctx_t ctx)
void ca_poly_set_si(ca_poly_t poly, slong c, ca_ctx_t ctx)

Sets `poly` to the constant polynomial $c$ (normalised, so the zero
polynomial if $c$ is proved zero).
:::

::: function
void ca_poly_set(ca_poly_t res, const ca_poly_t src, ca_ctx_t ctx)
void ca_poly_set_fmpz_poly(ca_poly_t res, const fmpz_poly_t src, ca_ctx_t ctx)
void ca_poly_set_fmpq_poly(ca_poly_t res, const fmpq_poly_t src, ca_ctx_t ctx)

Sets `res` to the polynomial `src`.
:::

::: function
void ca_poly_set_coeff_ca(ca_poly_t poly, slong n, const ca_t x, ca_ctx_t ctx)

Sets the coefficient of $x^n$ to `x`, extending `poly` with zeros if
needed, and normalises.
:::

::: function
void ca_poly_transfer(ca_poly_t res, ca_ctx_t res_ctx, const ca_poly_t src, ca_ctx_t src_ctx)

Sets `res`, which belongs to the context `res_ctx`, to the polynomial
`src` belonging to `src_ctx`, transferring each coefficient with
`ca_transfer`.
:::

## Conversions {#conversions}

::: function
void ca_poly_get_fexpr(fexpr_t res, const ca_poly_t A, ulong flags, ca_ctx_t ctx)

Sets `res` to a symbolic expression `Polynomial(List(c_0, ..., c_n))`
for `A`, wrapped in `Where(..., Def(a_1, ...), ...)` that defines the
extension numbers occurring in the coefficients, as for `ca_get_fexpr`,
which also explains `flags`.
:::

## Random generation {#random}

::: function
void ca_poly_randtest(ca_poly_t poly, flint_rand_t state, slong len, slong depth, slong bits, ca_ctx_t ctx)

Sets `poly` to a random polynomial of length up to `len` whose
coefficients are generated by `ca_randtest` with the parameters `depth`
and `bits`.
:::

::: function
void ca_poly_randtest_rational(ca_poly_t poly, flint_rand_t state, slong len, slong bits, ca_ctx_t ctx)

Sets `poly` to a random polynomial of length up to `len` with rational
coefficients of up to `bits` bits.
:::

## Input and output {#io}

::: function
void ca_poly_print(const ca_poly_t poly, ca_ctx_t ctx)

Prints the line `ca_poly of length n:` followed by the coefficients,
from the constant term upwards, one per line with `ca_print` and
indented by four spaces, and a blank line.
:::

::: function
void ca_poly_printn(const ca_poly_t poly, slong digits, ca_ctx_t ctx)

Prints decimal approximations of the coefficients with `digits`
significant digits, comma-separated in square brackets, followed by a
newline.
:::

## Degree and leading coefficient {#degree}

::: function
int ca_poly_is_proper(const ca_poly_t poly, ca_ctx_t ctx)

Returns 1 if all coefficients of `poly` are numbers (not special values)
and the leading coefficient is proved nonzero, so that `poly` is an
element of $\mathbb{C}[x]$ of known degree; also returns 1 for the zero
polynomial. Returns 0 otherwise.
:::

::: function
int ca_poly_make_monic(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)

Sets `res` to `poly` divided by its leading coefficient and returns 1,
or returns 0 if `poly` is zero or its leading coefficient is not proved
nonzero.
:::

::: function
void _ca_poly_reverse(ca_ptr res, ca_srcptr poly, slong len, slong n, ca_ctx_t ctx)
void ca_poly_reverse(ca_poly_t res, const ca_poly_t poly, slong n, ca_ctx_t ctx)

Sets `res` to the reversal $x^{n-1} f(1/x)$ of `poly` considered as a
polynomial of length `n`, truncating or zero-padding as needed. The
underscore version requires $0 < \mathit{len} \le n$. Wrappers of
`_gr_poly_reverse`.
:::

## Comparisons {#comparisons}

::: function
truth_t _ca_poly_check_equal(ca_srcptr poly1, slong len1, ca_srcptr poly2, slong len2, ca_ctx_t ctx)
truth_t ca_poly_check_equal(const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)

Checks whether `poly1` and `poly2` are equal as polynomials, comparing
coefficients with `ca_check_equal` (and the excess coefficients of the
longer one with zero). The underscore version requires $\mathit{len1}
\ge \mathit{len2}$. Wrappers of `_gr_poly_equal`.
:::

::: function
truth_t ca_poly_check_is_zero(const ca_poly_t poly, ca_ctx_t ctx)
truth_t ca_poly_check_is_one(const ca_poly_t poly, ca_ctx_t ctx)

Checks whether `poly` is the zero polynomial or the constant 1.
:::

## Arithmetic {#arithmetic}

::: function
void _ca_poly_shift_left(ca_ptr res, ca_srcptr poly, slong len, slong n, ca_ctx_t ctx)
void ca_poly_shift_left(ca_poly_t res, const ca_poly_t poly, slong n, ca_ctx_t ctx)
void _ca_poly_shift_right(ca_ptr res, ca_srcptr poly, slong len, slong n, ca_ctx_t ctx)
void ca_poly_shift_right(ca_poly_t res, const ca_poly_t poly, slong n, ca_ctx_t ctx)

Sets `res` to `poly` multiplied by $x^n$, or divided by $x^n$ with the
$n$ lowest coefficients discarded. Wrappers of `_gr_poly_shift_left`
and `_gr_poly_shift_right`.
:::

::: function
void ca_poly_neg(ca_poly_t res, const ca_poly_t src, ca_ctx_t ctx)

Sets `res` to $-\mathit{src}$.
:::

::: function
void _ca_poly_add(ca_ptr res, ca_srcptr poly1, slong len1, ca_srcptr poly2, slong len2, ca_ctx_t ctx)
void ca_poly_add(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
void _ca_poly_sub(ca_ptr res, ca_srcptr poly1, slong len1, ca_srcptr poly2, slong len2, ca_ctx_t ctx)
void ca_poly_sub(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)

Sets `res` to the sum or difference of `poly1` and `poly2`. The
underscore versions write $\max(\mathit{len1}, \mathit{len2})$
coefficients and allow aliasing. Wrappers of `_gr_poly_add` and
`_gr_poly_sub`.
:::

::: function
void _ca_poly_mul(ca_ptr C, ca_srcptr A, slong lenA, ca_srcptr B, slong lenB, ca_ctx_t ctx)
void ca_poly_mul(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)
void _ca_poly_mullow(ca_ptr C, ca_srcptr A, slong lenA, ca_srcptr B, slong lenB, slong n, ca_ctx_t ctx)
void ca_poly_mullow(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, slong n, ca_ctx_t ctx)

Sets `res` to the product of `poly1` and `poly2`, or to the product
truncated to length `n`. The underscore versions require $\mathit{lenA}
\ge \mathit{lenB} \ge 1$ and $0 < n \le \mathit{lenA} + \mathit{lenB} -
1$, and write `lenA + lenB - 1` or `n` coefficients.

::: details
If both factors have rational coefficients (and length at least 4), the
product is computed over $\Z$ after clearing denominators. If all
coefficients lie in one number field $\mathbb{Q}(a)$ and the lengths are
large enough, `_ca_poly_mullow_same_nf` is used. Otherwise the
classical algorithm is used, with a dedicated squaring routine; each
coefficient is a dot product computed with `ca_dot`.
:::
:::

::: function
void _ca_poly_mullow_same_nf(ca_ptr C, ca_srcptr A, slong Alen, ca_srcptr B, slong Blen, slong len, ca_field_t K, ca_ctx_t ctx)

Multiplies polynomials all of whose coefficients are rational numbers or
elements of the number field `K` = $\mathbb{Q}(a)$ of degree $m$, by
Kronecker substitution: the coefficients are written as polynomials in
$a$ over a common denominator, packed into a single polynomial over $\Z$
of length about $\mathit{Alen} \cdot (2m - 1)$, multiplied with
`_fmpz_poly_mullow` and unpacked, followed by reduction modulo the
minimal polynomial of $a$. Requirements as for `_ca_poly_mullow`.
:::

::: function
void ca_poly_mul_ca(ca_poly_t res, const ca_poly_t poly, const ca_t c, ca_ctx_t ctx)
void ca_poly_div_ca(ca_poly_t res, const ca_poly_t poly, const ca_t c, ca_ctx_t ctx)
void ca_poly_div_fmpz(ca_poly_t res, const ca_poly_t poly, const fmpz_t c, ca_ctx_t ctx)

Sets `res` to `poly` multiplied or divided by the scalar `c`. Division
by zero gives special values as coefficients (see
`_ca_vec_scalar_div_ca`).
:::

::: function
void _ca_poly_divrem_basecase(ca_ptr Q, ca_ptr R, ca_srcptr A, slong lenA, ca_srcptr B, slong lenB, const ca_t invB, ca_ctx_t ctx)
int ca_poly_divrem_basecase(ca_poly_t Q, ca_poly_t R, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
void _ca_poly_divrem(ca_ptr Q, ca_ptr R, ca_srcptr A, slong lenA, ca_srcptr B, slong lenB, const ca_t invB, ca_ctx_t ctx)
int ca_poly_divrem(ca_poly_t Q, ca_poly_t R, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
int ca_poly_div(ca_poly_t Q, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
int ca_poly_rem(ca_poly_t R, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)

If the leading coefficient of `B` is proved nonzero, sets `Q` and `R` to
the quotient and remainder of the division of `A` by `B`, with $\deg R <
\deg B$, and returns 1. Otherwise, in particular if `B` is zero, returns
0 and leaves the outputs unspecified.

::: details
The underscore versions take the precomputed inverse `invB` of the
leading coefficient of `B`, require $\mathit{lenA} \ge \mathit{lenB} \ge
1$, and write $\mathit{lenA} - \mathit{lenB} + 1$ coefficients to `Q`
and $\mathit{lenB} - 1$ to `R`. The algorithm is classical long
division (the default versions call the basecase versions).
:::
:::

::: function
void _ca_poly_pow_ui_trunc(ca_ptr res, ca_srcptr f, slong flen, ulong exp, slong len, ca_ctx_t ctx)
void ca_poly_pow_ui_trunc(ca_poly_t res, const ca_poly_t poly, ulong exp, slong len, ca_ctx_t ctx)
void _ca_poly_pow_ui(ca_ptr res, ca_srcptr f, slong flen, ulong exp, ca_ctx_t ctx)
void ca_poly_pow_ui(ca_poly_t res, const ca_poly_t poly, ulong exp, ca_ctx_t ctx)

Sets `res` to `poly` raised to the power `exp`, truncated to length
`len` for the *trunc* versions. Exponents larger than 2 use
`_gr_poly_pow_series_ui_binexp`. The underscore *trunc* version writes
`len` coefficients and the other one $\mathit{exp} (\mathit{flen} - 1) +
1$.
:::

## Evaluation and composition {#evaluation}

::: function
void _ca_poly_evaluate_horner(ca_t res, ca_srcptr f, slong len, const ca_t x, ca_ctx_t ctx)
void ca_poly_evaluate_horner(ca_t res, const ca_poly_t f, const ca_t a, ca_ctx_t ctx)
void _ca_poly_evaluate(ca_t res, ca_srcptr f, slong len, const ca_t x, ca_ctx_t ctx)
void ca_poly_evaluate(ca_t res, const ca_poly_t f, const ca_t a, ca_ctx_t ctx)

Sets `res` to $f(a)$, by Horner's rule or by the default algorithm of
`_gr_poly_evaluate`, of which these are wrappers.
:::

::: function
void _ca_poly_compose(ca_ptr res, ca_srcptr poly1, slong len1, ca_srcptr poly2, slong len2, ca_ctx_t ctx)
void ca_poly_compose(ca_poly_t res, const ca_poly_t poly1, const ca_poly_t poly2, ca_ctx_t ctx)

Sets `res` to the composition $\mathit{poly1}(\mathit{poly2}(x))$.
Wrappers of `_gr_poly_compose`; the underscore version writes
$(\mathit{len1} - 1)(\mathit{len2} - 1) + 1$ coefficients.
:::

## Derivative and integral {#calculus}

::: function
void _ca_poly_derivative(ca_ptr res, ca_srcptr poly, slong len, ca_ctx_t ctx)
void ca_poly_derivative(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)

Sets `res` to the derivative of `poly`. The underscore version writes
$\mathit{len} - 1$ coefficients and is a wrapper of
`_gr_poly_derivative`.
:::

::: function
void _ca_poly_integral(ca_ptr res, ca_srcptr poly, slong len, ca_ctx_t ctx)
void ca_poly_integral(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)

Sets `res` to the integral of `poly` with constant term zero. The
underscore version writes $\mathit{len} + 1$ coefficients.
:::

## Greatest common divisor {#gcd}

::: function
slong _ca_poly_gcd_euclidean(ca_ptr G, ca_srcptr A, slong lenA, ca_srcptr B, slong lenB, ca_ctx_t ctx)
int ca_poly_gcd_euclidean(ca_poly_t G, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)
slong _ca_poly_gcd(ca_ptr G, ca_srcptr A, slong lenA, ca_srcptr B, slong lenB, ca_ctx_t ctx)
int ca_poly_gcd(ca_poly_t G, const ca_poly_t A, const ca_poly_t B, ca_ctx_t ctx)

Sets `G` to the monic greatest common divisor of `A` and `B` (zero if
both are zero) and returns 1. Returns 0, leaving `G` unspecified, if a
zero test of a leading coefficient fails during the computation (in
particular if the leading coefficient of an input is not proved
nonzero).

::: details
The underscore versions require $\mathit{lenA} \ge \mathit{lenB} \ge 1$
and nonzero leading coefficients, write at most $\mathit{lenB}$
coefficients of a GCD that is not necessarily monic, and return its
length, or 0 on failure.

The *euclidean* versions run the Euclidean algorithm, testing each
remainder's leading coefficient for zero. The default versions first
handle rational polynomials by an integer GCD computation over $\Z$
(`_fmpz_poly_gcd`), then try to prove numerically that $A$ and $B$ are
coprime by evaluating their Sylvester resultant in ball arithmetic at
precision `CA_OPT_LOW_PREC` (a resultant ball excluding zero proves
$\gcd = 1$), and otherwise use the Euclidean algorithm.
:::
:::

## Power series {#series}

::: function
void _ca_poly_inv_series(ca_ptr res, ca_srcptr f, slong flen, slong len, ca_ctx_t ctx)
void ca_poly_inv_series(ca_poly_t res, const ca_poly_t f, slong len, ca_ctx_t ctx)

Sets `res` to the power series inverse $1/f \bmod x^{\mathit{len}}$.
The constant term must be invertible: if it is zero the result consists
of special values (unsigned infinity or *Undefined*), and *Unknown*
values if it cannot be decided whether it is zero. Rational series are
inverted over $\mathbb{Q}$; otherwise the classical recurrence or, for
long series with fast multiplication, Newton iteration is used.
:::

::: function
void _ca_poly_div_series(ca_ptr res, ca_srcptr f, slong flen, ca_srcptr g, slong glen, slong len, ca_ctx_t ctx)
void ca_poly_div_series(ca_poly_t res, const ca_poly_t f, const ca_poly_t g, slong len, ca_ctx_t ctx)

Sets `res` to the power series quotient $f/g \bmod x^{\mathit{len}}$.
As for `ca_poly_inv_series`, the constant term of $g$ must be nonzero;
remove common factors $x^k$ first when an exact cancellation is
expected.
:::

::: function
void _ca_poly_exp_series(ca_ptr res, ca_srcptr f, slong flen, slong len, ca_ctx_t ctx)
void ca_poly_exp_series(ca_poly_t res, const ca_poly_t f, slong len, ca_ctx_t ctx)

Sets `res` to the power series $\exp(f) \bmod x^{\mathit{len}}$. The
constant term is $\exp(f_0)$ computed with `ca_exp`. Rational series
use `_fmpq_poly_exp_series`; otherwise the basecase recurrence or
Newton iteration is used.
:::

::: function
void _ca_poly_log_series(ca_ptr res, ca_srcptr f, slong flen, slong len, ca_ctx_t ctx)
void ca_poly_log_series(ca_poly_t res, const ca_poly_t f, slong len, ca_ctx_t ctx)

Sets `res` to the power series $\log(f) \bmod x^{\mathit{len}}$,
computed as $\log(f_0) + \int f'/f$; the constant term must be nonzero.
Binomials $a + b x^d$ are handled directly.
:::

## Roots and factorisation {#roots}

::: function
int ca_poly_factor_squarefree(ca_t c, ca_poly_vec_t fac, ulong * exp, const ca_poly_t F, ca_ctx_t ctx)

Computes the squarefree factorisation $F = c \, f_1^{e_1} \cdots
f_k^{e_k}$, where the $f_i$ are monic, squarefree, pairwise coprime and
nonconstant and $e_1 < \cdots < e_k$. Sets `c` to the leading
coefficient, `fac` to the vector of the $f_i$ and the first $k$ entries
of `exp` to the $e_i$, and returns 1. Returns 0 if `F` is not proper or
if a GCD computation fails. The array `exp` must have room for $\deg F$
entries.
:::

::: function
int ca_poly_squarefree_part(ca_poly_t res, const ca_poly_t poly, ca_ctx_t ctx)

Sets `res` to the monic squarefree part $f / \gcd(f, f')$ of `poly` and
returns 1, or returns 0 if a GCD computation fails.
:::

::: function
void _ca_poly_set_roots(ca_ptr poly, ca_srcptr roots, const ulong * exp, slong len, ca_ctx_t ctx)
void ca_poly_set_roots(ca_poly_t poly, ca_vec_t roots, const ulong * exp, ca_ctx_t ctx)

Sets `poly` to the monic polynomial $(x - r_0)^{e_0} \cdots (x -
r_{n-1})^{e_{n-1}}$ with the roots $r_i$ in `roots` (of length `len`, or
the length of the vector) and the multiplicities $e_i$ in `exp`. The
product is computed by binary splitting. The underscore version writes
$e_0 + \cdots + e_{n-1} + 1$ coefficients.
:::

::: function
int _ca_poly_roots(ca_ptr roots, ca_srcptr poly, slong len, ca_ctx_t ctx)
int ca_poly_roots(ca_vec_t roots, ulong * exp, const ca_poly_t poly, ca_ctx_t ctx)

Attempts to compute all complex roots of `poly`. On success, returns 1,
sets `roots` to the distinct roots in no particular order and `exp` to
their multiplicities. Returns 0 on failure, leaving the outputs
unspecified. See [](ca_poly/roots).

::: details
The non-underscore version computes the squarefree factorisation and
then the roots of each factor; `exp` must have room for $\deg f$
entries. The underscore version assumes a squarefree polynomial, which
is not checked, and writes $\mathit{len} - 1$ roots. Roots of factors
with rational coefficients are computed as algebraic numbers with
`qqbar_roots_fmpz_poly`; other factors are solved by the formulas for
degree 1, 2 and 3. The computation fails if the leading coefficient is
not proved nonzero, if the squarefree factorisation fails, or if a
factor of degree 4 or higher has non-rational coefficients.
:::
:::

## Vectors of polynomials {#vectors}

::: type
ca_poly_vec_struct
ca_poly_vec_t

A vector of polynomials, with the fields `entries` (an array of `alloc`
initialised polynomials), `length` and `alloc`. A `ca_poly_vec_t` is an
array of one structure.
:::

::: function
ca_poly_struct * _ca_poly_vec_init(slong len, ca_ctx_t ctx)
void ca_poly_vec_init(ca_poly_vec_t res, slong len, ca_ctx_t ctx)

Allocates an array of `len` zero polynomials, or initialises `res` as a
vector of `len` zero polynomials.
:::

::: function
void _ca_poly_vec_fit_length(ca_poly_vec_t vec, slong len, ca_ctx_t ctx)

Ensures that `vec` has space for `len` polynomials, without changing
its length.
:::

::: function
void ca_poly_vec_set_length(ca_poly_vec_t vec, slong len, ca_ctx_t ctx)

Sets the length of `vec` to `len`, appending zero polynomials or
truncating as needed.
:::

::: function
void _ca_poly_vec_clear(ca_poly_struct * v, slong len, ca_ctx_t ctx)
void ca_poly_vec_clear(ca_poly_vec_t vec, ca_ctx_t ctx)

Clears the polynomials and frees the array.
:::

::: function
void ca_poly_vec_append(ca_poly_vec_t vec, const ca_poly_t f, ca_ctx_t ctx)

Appends a copy of `f` to `vec`.
:::
