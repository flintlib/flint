---
title: fq_poly.h
header: fq_poly.h
template-headers: fq_poly_templates
description: Dense univariate polynomials over finite fields $\mathbb{F}_q$ of arbitrary characteristic, with Kronecker substitution multiplication and generic fast division, GCD, composition and power series algorithms.
category: Finite fields
order: 8
---

This module implements dense univariate polynomials over a finite field
$\mathbb{F}_q$ in the representation of [](fq). Multiplication maps the
problem to a product of integer polynomials, by Kronecker substitution
(packing each coefficient into an integer) or by substituting
$y \mapsto x^{2d-1}$ for the generator of $\mathbb{F}_q$, so that it
benefits from FLINT's fast multiplication in $\Z[x]$. Most other
algorithms (division, GCD, composition, power series, evaluation) are thin
wrappers around the generic implementations of [](gr_poly), which call the
multiplication of this module. Factorisation is in `fq_poly_factor.h`.

The interface is generated from `fq_poly_templates.h` and is the same for
`fq_nmod_poly` and `fq_zech_poly`, apart from the multiplication
algorithms. `fq_default_poly` provides the same functions for any finite
field, and is the recommended interface.

::: summary algorithms/polynomials/multiplication :::

## Types and conventions {#types}

::: type
fq_poly_struct
fq_poly_t

A polynomial over $\mathbb{F}_q$: an array `coeffs` of `alloc`
initialised elements (`fq_struct`), of which the first `length` are the
coefficients, the coefficient of $x^i$ at index $i$. `fq_poly_t` is an
array of one structure. The context is not stored and must be passed to
every function.
:::

- A polynomial is *normalised* if its length is 0 or its leading
  coefficient is nonzero. All functions require normalised inputs and
  produce normalised outputs. The zero polynomial has length 0 and
  degree $-1$.
- Functions without a leading underscore allow outputs to be aliased with
  inputs, unless stated otherwise; several outputs of one function (such
  as a quotient and a remainder) must be distinct.
- Functions with a leading underscore work on arrays `(poly, len)` of
  initialised elements, require the output to have room for the stated
  number of coefficients, write that number of coefficients including
  leading zeros, and in general do not allow aliasing; each one states
  its conditions.
- The context must describe a field; division by a polynomial with
  a non-invertible leading coefficient throws `FLINT_ERROR` (division by
  zero). The `_f` variants of the division and GCD functions also work
  for a reducible defining polynomial and return a factor of it instead
  of failing.

## Memory management {#memory}

::: function
void fq_poly_init(fq_poly_t poly, const fq_ctx_t ctx)
void fq_poly_init2(fq_poly_t poly, slong alloc, const fq_ctx_t ctx)

Initialise `poly` as the zero polynomial. `fq_poly_init2` allocates space
for `alloc` coefficients.
:::

::: function
void fq_poly_realloc(fq_poly_t poly, slong alloc, const fq_ctx_t ctx)

Changes the allocation of `poly` to `alloc` coefficients, truncating and
normalising the polynomial if it is longer. For `alloc` = 0 the
coefficients are freed.
:::

::: function
void fq_poly_fit_length(fq_poly_t poly, slong len, const fq_ctx_t ctx)

Ensures that `poly` has space for at least `len` coefficients, growing
the allocation to at least twice its old size if it must grow. The
polynomial is not changed.
:::

::: function
void fq_poly_clear(fq_poly_t poly, const fq_ctx_t ctx)

Frees the memory used by `poly`.
:::

::: function
void _fq_poly_set_length(fq_poly_t poly, slong len, const fq_ctx_t ctx)

Sets the length of `poly` to `len`, setting the coefficients from index
`len` up to the old length to zero. Does not normalise. `poly` must have
space for `len` coefficients.
:::

::: function
void _fq_poly_normalise(fq_poly_t poly, const fq_ctx_t ctx)
void _fq_poly_normalise2(const fq_struct * poly, slong * length, const fq_ctx_t ctx)

Decrease the length of `poly`, respectively the variable `*length` of
the array `(poly, *length)`, until the leading coefficient is nonzero or
the length is zero.
:::

::: function
void fq_poly_truncate(fq_poly_t poly, slong len, const fq_ctx_t ctx)

Truncates `poly` to length at most `len` (that is, reduces it modulo
$x^{\mathit{len}}$) and normalises it.
:::

::: function
void fq_poly_set_trunc(fq_poly_t poly1, fq_poly_t poly2, slong len, const fq_ctx_t ctx)

Sets `poly1` to `poly2` truncated to length `len`. `poly2` is not modified
unless it is the same object as `poly1`.
:::

## Length, degree and comparison {#properties}

::: function
slong fq_poly_length(const fq_poly_t poly, const fq_ctx_t ctx)
slong fq_poly_degree(const fq_poly_t poly, const fq_ctx_t ctx)

Return the length, respectively the degree, of `poly`; the zero
polynomial has length 0 and degree $-1$.
:::

::: function
fq_struct * fq_poly_lead(const fq_poly_t poly, const fq_ctx_t ctx)

Returns a pointer to the leading coefficient of `poly`, or `NULL` for the
zero polynomial.
:::

::: function
slong _fq_poly_hamming_weight(const fq_struct * op, slong len, const fq_ctx_t ctx)
slong fq_poly_hamming_weight(const fq_poly_t op, const fq_ctx_t ctx)

Return the number of nonzero coefficients.
:::

::: function
int fq_poly_equal(const fq_poly_t poly1, const fq_poly_t poly2, const fq_ctx_t ctx)

Returns whether the two polynomials are equal.
:::

::: function
int fq_poly_equal_trunc(const fq_poly_t poly1, const fq_poly_t poly2, slong n, const fq_ctx_t ctx)

Returns whether the two polynomials are equal modulo $x^n$.
:::

::: function
int fq_poly_equal_fq(const fq_poly_t poly, const fq_t c, const fq_ctx_t ctx)

Returns whether `poly` is the constant polynomial $c$.
:::

::: function
int fq_poly_is_zero(const fq_poly_t poly, const fq_ctx_t ctx)
int fq_poly_is_one(const fq_poly_t op, const fq_ctx_t ctx)
int fq_poly_is_gen(const fq_poly_t poly, const fq_ctx_t ctx)

Return whether the polynomial is 0, 1, respectively $x$.
:::

::: function
int fq_poly_is_unit(const fq_poly_t op, const fq_ctx_t ctx)

Returns whether `op` is a unit of $\mathbb{F}_q[x]$, that is, a nonzero
constant.
:::

## Assignment and conversion {#assignment}

::: function
void _fq_poly_set(fq_struct * rop, const fq_struct * op, slong len, const fq_ctx_t ctx)
void fq_poly_set(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

Copy `op` into `rop`.
:::

::: function
void fq_poly_set_fq(fq_poly_t poly, const fq_t c, const fq_ctx_t ctx)

Sets `poly` to the constant polynomial $c$.
:::

::: function
void fq_poly_set_fmpz_mod_poly(fq_poly_t rop, const fmpz_mod_poly_t op, const fq_ctx_t ctx)
void fq_poly_set_nmod_poly(fq_poly_t rop, const nmod_poly_t op, const fq_ctx_t ctx)

Set `rop` to the polynomial `op` over $\mathbb{F}_p$, with coefficients
mapped to $\mathbb{F}_p \subseteq \mathbb{F}_q$. The modulus of `op` must
be $p$; it is not checked.
:::

::: function
void fq_poly_swap(fq_poly_t op1, fq_poly_t op2, const fq_ctx_t ctx)

Swaps `op1` and `op2` efficiently.
:::

::: function
void _fq_poly_zero(fq_struct * rop, slong len, const fq_ctx_t ctx)
void fq_poly_zero(fq_poly_t poly, const fq_ctx_t ctx)

Set the `len` coefficients of `rop`, respectively the polynomial `poly`,
to zero.
:::

::: function
void fq_poly_one(fq_poly_t poly, const fq_ctx_t ctx)
void fq_poly_gen(fq_poly_t f, const fq_ctx_t ctx)

Set the polynomial to 1, respectively to $x$.
:::

::: function
void _fq_poly_make_monic(fq_struct * rop, const fq_struct * op, slong length, const fq_ctx_t ctx)
void fq_poly_make_monic(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

Set `rop` to `op` divided by its leading coefficient. The zero polynomial
is mapped to zero; the underscore version requires a nonzero leading
coefficient `op[length - 1]` and allows aliasing.
:::

::: function
void _fq_poly_reverse(fq_struct * res, const fq_struct * poly, slong len, slong n, const fq_ctx_t ctx)
void fq_poly_reverse(fq_poly_t res, const fq_poly_t poly, slong n, const fq_ctx_t ctx)

Set `res` to the reversal $x^{n-1} f(1/x)$ of the polynomial $f$ =
`poly` truncated to length $n$, considered as a polynomial of length $n$.
The underscore version requires $\mathit{len} \le n$, writes $n$
coefficients and allows `res` to be `poly`.
:::

::: function
ulong fq_poly_deflation(const fq_poly_t input, const fq_ctx_t ctx)

Returns the largest $k$ such that `input` is a polynomial in $x^k$. For
a polynomial of length at most 1 the result is its length.
:::

::: function
void fq_poly_deflate(fq_poly_t result, const fq_poly_t input, ulong deflation, const fq_ctx_t ctx)
void fq_poly_inflate(fq_poly_t result, const fq_poly_t input, ulong inflation, const fq_ctx_t ctx)

`fq_poly_deflate` sets `result` to $g$ with $g(x^k) = f$ for the
polynomial $f$ = `input` and $k$ = `deflation`, which must divide the
exponents of all nonzero terms of $f$; it throws `FLINT_DIVZERO` for
$k = 0$. `fq_poly_inflate` sets `result` to $f(x^k)$ for
$k$ = `inflation`.
:::

## Coefficients {#coefficients}

::: function
void fq_poly_get_coeff(fq_t x, const fq_poly_t poly, slong n, const fq_ctx_t ctx)

Sets $x$ to the coefficient of $x^n$ in `poly`, which is zero for $n$ at
least the length.
:::

::: function
void fq_poly_set_coeff(fq_poly_t poly, slong n, const fq_t x, const fq_ctx_t ctx)
void fq_poly_set_coeff_fmpz(fq_poly_t poly, slong n, const fmpz_t x, const fq_ctx_t ctx)

Set the coefficient of $x^n$, $n \ge 0$, to the element $x$, respectively
to the image of the integer $x$ in $\mathbb{F}_p$, extending and
normalising the polynomial as needed.
:::

## Random polynomials {#random}

::: function
void fq_poly_randtest(fq_poly_t f, flint_rand_t state, slong len, const fq_ctx_t ctx)
void fq_poly_randtest_not_zero(fq_poly_t f, flint_rand_t state, slong len, const fq_ctx_t ctx)

Set $f$ to a random polynomial of length at most `len`, with coefficients
from `fq_randtest`; the second version returns a nonzero polynomial and
requires `len` $> 0$ (it throws otherwise).
:::

::: function
void fq_poly_randtest_monic(fq_poly_t f, flint_rand_t state, slong len, const fq_ctx_t ctx)

Sets $f$ to a random monic polynomial of length `len` $\ge 1$.
:::

::: function
void fq_poly_randtest_irreducible(fq_poly_t f, flint_rand_t state, slong len, const fq_ctx_t ctx)

Sets $f$ to a random monic irreducible polynomial of length `len`
$\ge 2$. Random monic polynomials are tested by computing
$\gcd(f, x^{q^i} - x)$ for $i \le (\mathit{len} - 1)/2$ until one is
irreducible, which takes about `len` attempts on average.
:::

## Addition and subtraction {#addition}

::: function
void _fq_poly_add(fq_struct * res, const fq_struct * poly1, slong len1, const fq_struct * poly2, slong len2, const fq_ctx_t ctx)
void _fq_poly_sub(fq_struct * res, const fq_struct * poly1, slong len1, const fq_struct * poly2, slong len2, const fq_ctx_t ctx)

Set `(res, max(len1, len2))` to the sum, respectively the difference, of
the two arrays. Aliasing of `res` with the inputs is allowed.
:::

::: function
void fq_poly_add(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)
void fq_poly_sub(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

Set `rop` to $\mathit{op1} + \mathit{op2}$, respectively
$\mathit{op1} - \mathit{op2}$.
:::

::: function
void fq_poly_add_si(fq_poly_t rop, const fq_poly_t op1, slong c, const fq_ctx_t ctx)

Sets `rop` to $\mathit{op1} + c$, where $c$ is mapped to $\mathbb{F}_p$.
:::

::: function
void fq_poly_add_series(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)
void fq_poly_sub_series(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)

Set `rop` to the sum, respectively the difference, of the two
polynomials truncated to length $n$.
:::

::: function
void _fq_poly_neg(fq_struct * rop, const fq_struct * op, slong len, const fq_ctx_t ctx)
void fq_poly_neg(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

Set `rop` to $-\mathit{op}$.
:::

## Scalar multiplication and division {#scalar}

::: function
void _fq_poly_scalar_mul_fq(fq_struct * rop, const fq_struct * op, slong len, const fq_t x, const fq_ctx_t ctx)
void fq_poly_scalar_mul_fq(fq_poly_t rop, const fq_poly_t op, const fq_t x, const fq_ctx_t ctx)

Set `rop` to $x \cdot \mathit{op}$.
:::

::: function
void _fq_poly_scalar_div_fq(fq_struct * rop, const fq_struct * op, slong len, const fq_t x, const fq_ctx_t ctx)
void fq_poly_scalar_div_fq(fq_poly_t rop, const fq_poly_t op, const fq_t x, const fq_ctx_t ctx)

Set `rop` to $\mathit{op} / x$. Throws `FLINT_ERROR` if $x$ is zero.
:::

::: function
void _fq_poly_scalar_addmul_fq(fq_struct * rop, const fq_struct * op, slong len, const fq_t x, const fq_ctx_t ctx)
void fq_poly_scalar_addmul_fq(fq_poly_t rop, const fq_poly_t op, const fq_t x, const fq_ctx_t ctx)
void _fq_poly_scalar_submul_fq(fq_struct * rop, const fq_struct * op, slong len, const fq_t x, const fq_ctx_t ctx)
void fq_poly_scalar_submul_fq(fq_poly_t rop, const fq_poly_t op, const fq_t x, const fq_ctx_t ctx)

Add $x \cdot \mathit{op}$ to `rop`, respectively subtract it. The
underscore versions work on the first `len` coefficients of both arrays.
:::

## Multiplication {#multiplication}

[](fq_poly/mul) explains the multiplication algorithms and when each is
used.

::: function
void _fq_poly_mul(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, const fq_ctx_t ctx)
void fq_poly_mul(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

Set `rop` to $\mathit{op1} \cdot \mathit{op2}$, choosing the algorithm.
The underscore version writes `len1 + len2 - 1` coefficients, requires
`len1`, `len2` $\ge 1$ and does not allow aliasing.

::: details
For `fq` the classical algorithm is used if both lengths are below
`FQ_MUL_CLASSICAL_CUTOFF` (6), the reordering algorithm
`_fq_poly_mul_reorder` if the field degree is below 4, and Kronecker
substitution `_fq_poly_mul_KS` otherwise.
:::
:::

::: function
void _fq_poly_mul_classical(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, const fq_ctx_t ctx)
void fq_poly_mul_classical(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

Set `rop` to the product, computed by the classical algorithm with one
dot product (reduced once) per coefficient (a wrapper of
`_gr_poly_mullow_classical`). The underscore version requires
`len1`, `len2` $\ge 1$ and does not allow aliasing.
:::

::: function
void _fq_poly_mul_KS(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, const fq_ctx_t ctx)
void fq_poly_mul_KS(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

Set `rop` to the product, computed by Kronecker substitution: every
coefficient is packed into one integer, the integer polynomials are
multiplied with `_fmpz_poly_mul`, and the coefficients are unpacked and
reduced. The underscore version allows zero padding of the inputs (it
writes `len1 + len2 - 1` coefficients) and aliasing.
:::

::: function
void _fq_poly_mul_univariate(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, const fq_ctx_t ctx)
void fq_poly_mul_univariate(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

Set `rop` to the product, computed by viewing the polynomials as
bivariate polynomials in $x$ and the generator $y$ of $\mathbb{F}_q$ and
substituting $y = x^{2d-1}$, which gives one product of integer
polynomials of length about $(2d - 1)\,\mathit{len}$. The underscore
version requires `len1`, `len2` $\ge 1$, a nonzero leading coefficient
of both inputs and writes `len1 + len2 - 1` coefficients.
:::

::: function
void _fq_poly_mul_reorder(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, const fq_ctx_t ctx)
void fq_poly_mul_reorder(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

Set `rop` to the product, computed by exchanging the roles of $x$ and
the generator $y$: the inputs are written as $\sum_{j < d} y^j a_j(x)$
with $a_j \in \Z[x]$, the $d^2$ products $a_j b_k$ are computed with
fast integer polynomial multiplication, and the result is reduced
modulo the defining polynomial. This is efficient for small $d$ and long
polynomials. The underscore version requires `len1`, `len2` $\ge 1$ and
allows aliasing.
:::

::: function
void _fq_poly_mullow(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, slong n, const fq_ctx_t ctx)
void fq_poly_mullow(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)

Set `rop` to the product truncated to length $n$, that is, the low $n$
coefficients of $\mathit{op1} \cdot \mathit{op2}$. The underscore version
requires $0 < n \le \mathit{len1} + \mathit{len2} - 1$ and
`len1`, `len2` $\ge 1$, writes $n$ coefficients and does not allow
aliasing. The classical algorithm is used for $n$ below
`FQ_MULLOW_CLASSICAL_CUTOFF` (6) or short inputs, and Kronecker
substitution otherwise.
:::

::: function
void _fq_poly_mullow_classical(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, slong n, const fq_ctx_t ctx)
void fq_poly_mullow_classical(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)
void _fq_poly_mullow_KS(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, slong n, const fq_ctx_t ctx)
void fq_poly_mullow_KS(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)
void _fq_poly_mullow_univariate(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, slong n, const fq_ctx_t ctx)
void fq_poly_mullow_univariate(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong n, const fq_ctx_t ctx)

Truncated products with a fixed algorithm: classical, Kronecker
substitution, or the substitution $y = x^{2d-1}$. The underscore versions
require $0 < n \le \mathit{len1} + \mathit{len2} - 1$ and write $n$
coefficients; `_fq_poly_mullow_KS` allows zero padding and zero inputs.
:::

::: function
void _fq_poly_mulhigh(fq_struct * res, const fq_struct * poly1, slong len1, const fq_struct * poly2, slong len2, slong n, fq_ctx_t ctx)
void fq_poly_mulhigh(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong start, const fq_ctx_t ctx)

Compute the coefficients of index at least $n$ (`start`) of the product.
`fq_poly_mulhigh` sets `rop` to a polynomial of length
`len1 + len2 - 1` whose coefficients from `start` on are those of the
product and whose lower coefficients are zero; it uses the classical
algorithm. `_fq_poly_mulhigh` requires `len1` $\ge$ `len2` $\ge 1$ and
space for `len1 + len2 - 1` coefficients in `res`; its coefficients below
$n$ are arbitrary (it computes the full product by Kronecker substitution
for inputs of length at least 6).
:::

::: function
void _fq_poly_mulhigh_classical(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, slong start, const fq_ctx_t ctx)
void fq_poly_mulhigh_classical(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, slong start, const fq_ctx_t ctx)

As the previous functions, with the classical algorithm, which computes
only the coefficients from `start` on and sets the lower ones to zero.
The underscore version requires `len1` $\ge$ `len2` $\ge 1$ and does not
allow aliasing.
:::

::: function
void _fq_poly_mulmod(fq_struct * res, const fq_struct * poly1, slong len1, const fq_struct * poly2, slong len2, const fq_struct * f, slong lenf, const fq_ctx_t ctx)
void fq_poly_mulmod(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t f, const fq_ctx_t ctx)

Set `res` to $\mathit{poly1} \cdot \mathit{poly2} \bmod f$. The underscore
version requires $\mathit{len1} + \mathit{len2} > \mathit{lenf}$ (the
product must need reduction), writes $\mathit{lenf} - 1$ coefficients
and needs space for `len1 + len2 - 1` coefficients in `res`, which must
not be aliased with $f$. `fq_poly_mulmod` throws `FLINT_DIVZERO` if $f = 0$.
:::

::: function
void _fq_poly_mulmod_preinv(fq_struct * res, const fq_struct * poly1, slong len1, const fq_struct * poly2, slong len2, const fq_struct * f, slong lenf, const fq_struct * finv, slong lenfinv, const fq_ctx_t ctx)
void fq_poly_mulmod_preinv(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)

As `fq_poly_mulmod`, with the precomputed inverse `finv` of the reversal
of $f$ modulo $x^{\deg f + 1}$ (`fq_poly_reverse` followed by
`fq_poly_inv_series`). The product must satisfy
$\mathit{len1} + \mathit{len2} - 1 \le 2\,\mathit{lenf} - 2$; the inputs
need not be reduced. The underscore version does not allow aliasing.
:::

## Squaring {#squaring}

::: function
void _fq_poly_sqr(fq_struct * rop, const fq_struct * op, slong len, const fq_ctx_t ctx)
void fq_poly_sqr(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

Set `rop` to $\mathit{op}^2$, choosing the algorithm as for
multiplication. The underscore version writes `2 * len - 1`
coefficients, requires `len` $\ge 1$ and does not allow aliasing.
:::

::: function
void _fq_poly_sqr_classical(fq_struct * rop, const fq_struct * op, slong len, const fq_ctx_t ctx)
void fq_poly_sqr_classical(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)
void _fq_poly_sqr_KS(fq_struct * rop, const fq_struct * op, slong len, const fq_ctx_t ctx)
void fq_poly_sqr_KS(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)
void _fq_poly_sqr_reorder(fq_struct * rop, const fq_struct * op, slong len, const fq_ctx_t ctx)
void fq_poly_sqr_reorder(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

Squaring with a fixed algorithm: classical, Kronecker substitution, or
reordering of the variables (see `fq_poly_mul_reorder`). The underscore
versions write `2 * len - 1` coefficients; the classical version does not
allow aliasing.
:::

## Powering {#powering}

::: function
void _fq_poly_pow(fq_struct * rop, const fq_struct * op, slong len, ulong e, const fq_ctx_t ctx)
void fq_poly_pow(fq_poly_t rop, const fq_poly_t op, ulong e, const fq_ctx_t ctx)

Set `rop` to $\mathit{op}^e$, with $\mathit{op}^0 = 1$. The underscore
version requires $e > 0$, `len` $> 0$ and space for
$e(\mathit{len} - 1) + 1$ coefficients, and does not allow aliasing. It
uses binary exponentiation (`_gr_poly_pow_ui_binexp`).
:::

::: function
void _fq_poly_pow_trunc(fq_struct * res, const fq_struct * poly, ulong e, slong trunc, const fq_ctx_t ctx)
void fq_poly_pow_trunc(fq_poly_t res, const fq_poly_t poly, ulong e, slong trunc, const fq_ctx_t ctx)
void _fq_poly_pow_trunc_binexp(fq_struct * res, const fq_struct * poly, ulong e, slong trunc, const fq_ctx_t ctx)
void fq_poly_pow_trunc_binexp(fq_poly_t res, const fq_poly_t poly, ulong e, slong trunc, const fq_ctx_t ctx)

Set `res` to $\mathit{poly}^e \bmod x^{\mathit{trunc}}$, by binary
exponentiation with truncated multiplications. The underscore versions
require the input to be zero-padded to length `trunc`, `trunc` $> 0$,
$e > 1$ and space for `trunc` coefficients, and do not allow aliasing.
:::

::: function
void _fq_poly_powmod_ui_binexp(fq_struct * res, const fq_struct * poly, ulong e, const fq_struct * f, slong lenf, const fq_ctx_t ctx)
void fq_poly_powmod_ui_binexp(fq_poly_t res, const fq_poly_t poly, ulong e, const fq_poly_t f, const fq_ctx_t ctx)
void _fq_poly_powmod_fmpz_binexp(fq_struct * res, const fq_struct * poly, const fmpz_t e, const fq_struct * f, slong lenf, const fq_ctx_t ctx)
void fq_poly_powmod_fmpz_binexp(fq_poly_t res, const fq_poly_t poly, const fmpz_t e, const fq_poly_t f, const fq_ctx_t ctx)

Set `res` to $\mathit{poly}^e \bmod f$ for $e \ge 0$, by binary
exponentiation. The underscore versions require $e > 0$,
`lenf` $> 1$, `poly` reduced modulo $f$ and zero-padded to length
`lenf - 1`, and space for `lenf - 1` coefficients in `res`.
:::

::: function
void _fq_poly_powmod_ui_binexp_preinv(fq_struct * res, const fq_struct * poly, ulong e, const fq_struct * f, slong lenf, const fq_struct * finv, slong lenfinv, const fq_ctx_t ctx)
void fq_poly_powmod_ui_binexp_preinv(fq_poly_t res, const fq_poly_t poly, ulong e, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)
void _fq_poly_powmod_fmpz_binexp_preinv(fq_struct * res, const fq_struct * poly, const fmpz_t e, const fq_struct * f, slong lenf, const fq_struct * finv, slong lenfinv, const fq_ctx_t ctx)
void fq_poly_powmod_fmpz_binexp_preinv(fq_poly_t res, const fq_poly_t poly, const fmpz_t e, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)

As the previous functions, with the precomputed inverse `finv` of the
reversal of $f$, so that every reduction is a Newton division with two
multiplications.
:::

::: function
void _fq_poly_powmod_fmpz_sliding_preinv(fq_struct * res, const fq_struct * poly, const fmpz_t e, ulong k, const fq_struct * f, slong lenf, const fq_struct * finv, slong lenfinv, const fq_ctx_t ctx)
void fq_poly_powmod_fmpz_sliding_preinv(fq_poly_t res, const fq_poly_t poly, const fmpz_t e, ulong k, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)

As `fq_poly_powmod_fmpz_binexp_preinv`, using sliding window
exponentiation with windows of $k$ bits; $k = 0$ chooses $k$ from the
size of $e$.
:::

::: function
void _fq_poly_powmod_x_fmpz_preinv(fq_struct * res, const fmpz_t e, const fq_struct * f, slong lenf, const fq_struct * finv, slong lenfinv, const fq_ctx_t ctx)
void fq_poly_powmod_x_fmpz_preinv(fq_poly_t res, const fmpz_t e, const fq_poly_t f, const fq_poly_t finv, const fq_ctx_t ctx)

Set `res` to $x^e \bmod f$, using sliding window exponentiation, in
which multiplications by $x$ are shifts. The underscore version requires
$e > 0$, `lenf` $> 2$ and space for `lenf - 1` coefficients. This is the
first step of distinct-degree factorisation and irreducibility testing,
with $e = q$.
:::

## Shifting {#shifting}

::: function
void _fq_poly_shift_left(fq_struct * rop, const fq_struct * op, slong len, slong n, const fq_ctx_t ctx)
void fq_poly_shift_left(fq_poly_t rop, const fq_poly_t op, slong n, const fq_ctx_t ctx)

Set `rop` to $x^n \cdot \mathit{op}$. The underscore version writes
`len + n` coefficients and allows aliasing.
:::

::: function
void _fq_poly_shift_right(fq_struct * rop, const fq_struct * op, slong len, slong n, const fq_ctx_t ctx)
void fq_poly_shift_right(fq_poly_t rop, const fq_poly_t op, slong n, const fq_ctx_t ctx)

Set `rop` to $\lfloor \mathit{op} / x^n \rfloor$, dropping the $n$
lowest coefficients; the result is zero if $n$ is at least the length.
The underscore version requires `len` $> n$, writes `len - n`
coefficients and allows aliasing (the top $n$ coefficients of `op` are
then permuted, not zeroed).
:::

## Euclidean division {#division}

Division uses the generic implementations of `gr_poly`: the classical
algorithm with a precomputed inverse of the leading coefficient for short
quotients, and Newton iteration ([](algorithms/polynomials/newton)) for
long ones.

::: summary algorithms/polynomials/newton :::

::: function
void _fq_poly_divrem(fq_struct * Q, fq_struct * R, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_t invB, const fq_ctx_t ctx)

Computes `(Q, lenA - lenB + 1)` and `(R, lenB - 1)` with $A = BQ + R$ and
$\deg R < \deg B$, given the inverse `invB` of the leading coefficient of
$B$. Requires `lenA` $\ge$ `lenB` $> 0$. Allows aliasing of $R$ and $A$
but no other aliasing.
:::

::: function
void fq_poly_divrem(fq_poly_t Q, fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

Computes $Q$ and $R$ with $A = BQ + R$ and $\deg R < \deg B$. $B$ must
be nonzero. If $\deg A < \deg B$, then $Q = 0$ and $R = A$.
:::

::: function
void _fq_poly_div(fq_struct * Q, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_t invB, const fq_ctx_t ctx)
void fq_poly_div(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

Compute only the quotient $Q$ of the division of $A$ by $B$. The
underscore version requires `lenA` $\ge$ `lenB` $> 0$ and writes
`lenA - lenB + 1` coefficients; `fq_poly_div` requires $B \ne 0$.
:::

::: function
void fq_poly_div_basecase(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

The same as `fq_poly_div`; kept for compatibility with FLINT 2.
:::

::: function
void _fq_poly_rem(fq_struct * R, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_t invB, const fq_ctx_t ctx)
void fq_poly_rem(fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

Compute only the remainder $R$ of the division of $A$ by $B \ne 0$. The
underscore version writes `lenB - 1` coefficients and also accepts
`lenA` $<$ `lenB`; `invB` is not used.
:::

::: function
void _fq_poly_divrem_f(fq_t f, fq_struct * Q, fq_struct * R, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_ctx_t ctx)
void fq_poly_divrem_f(fq_t f, fq_poly_t Q, fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

As `fq_poly_divrem`, for a possibly reducible defining polynomial: if the
leading coefficient of $B$ is invertible, sets $f = 1$ and computes $Q$
and $R$; otherwise sets $f$ to a nontrivial factor of the defining
polynomial (see `fq_gcdinv`) and leaves $Q$ and $R$ unchanged. $B$ must
be nonzero.
:::

::: function
void _fq_poly_div_newton_n_preinv(fq_struct * Q, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_struct * Binv, slong lenBinv, const fq_ctx_t ctx)
void fq_poly_div_newton_n_preinv(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B, const fq_poly_t Binv, const fq_ctx_t ctx)

Compute the quotient $Q$ of $A$ by $B$ with two multiplications, given
the inverse `Binv` of the reversal of $B$ modulo $x^{\mathit{lenB}}$:
the reversed quotient is the product of the reversal of $A$ and `Binv`
modulo $x^{\mathit{lenA} - \mathit{lenB} + 1}$. Requires
$\mathit{lenA} \le 2\,\mathit{lenB} - 1$; the underscore version writes
`lenA - lenB + 1` coefficients and does not read $B$. The non-underscore
version throws `FLINT_DIVZERO` for $B = 0$.
:::

::: function
void _fq_poly_divrem_newton_n_preinv(fq_struct * Q, fq_struct * R, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_struct * Binv, slong lenBinv, const fq_ctx_t ctx)
void fq_poly_divrem_newton_n_preinv(fq_poly_t Q, fq_poly_t R, const fq_poly_t A, const fq_poly_t B, const fq_poly_t Binv, const fq_ctx_t ctx)

As the previous functions, also computing the remainder
$R = A - BQ$, with one more multiplication. Requires
$\mathit{lenA} \le 2\,\mathit{lenB} - 1$.
:::

::: function
int _fq_poly_divides(fq_struct * Q, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_t invB, const fq_ctx_t ctx)
int fq_poly_divides(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

Return 1 and set $Q = A / B$ if $B$ divides $A$, and return 0 otherwise
(then $Q$ is undefined). The underscore version requires
`lenA` $\ge$ `lenB` $> 0$, space for `lenA - lenB + 1` coefficients in
$Q$ and no aliasing. `fq_poly_divides` throws `FLINT_DIVZERO` if $B = 0$. The test performs a full division with remainder.
:::

::: function
ulong fq_poly_remove(fq_poly_t f, const fq_poly_t g, const fq_ctx_t ctx)

Removes the highest power of $g$ that divides $f$ from $f$ and returns
its exponent. $g$ must have degree at least 1.
:::

## Power series {#series}

::: function
void _fq_poly_inv_series_newton(fq_struct * Qinv, const fq_struct * Q, slong n, const fq_t cinv, const fq_ctx_t ctx)
void _fq_poly_inv_series(fq_struct * Qinv, const fq_struct * Q, slong n, const fq_t cinv, const fq_ctx_t ctx)

Set `(Qinv, n)` to the inverse of the power series `(Q, n)` modulo
$x^n$, given the inverse `cinv` of the constant term. Requires $n > 0$
and no aliasing. The basecase algorithm is used for small $n$ and Newton
iteration (`_gr_poly_inv_series_newton`) otherwise; the two functions
are the same.
:::

::: function
void fq_poly_inv_series_newton(fq_poly_t Qinv, const fq_poly_t Q, slong n, const fq_ctx_t ctx)
void fq_poly_inv_series(fq_poly_t Qinv, const fq_poly_t Q, slong n, const fq_ctx_t ctx)

Set `Qinv` to the inverse of the power series $Q$ modulo $x^n$. The
constant term of $Q$ must be nonzero (otherwise `FLINT_ERROR` is thrown)
and $n \ge 1$.
:::

::: function
void _fq_poly_div_series(fq_struct * Q, const fq_struct * A, slong Alen, const fq_struct * B, slong Blen, slong n, const fq_ctx_t ctx)
void fq_poly_div_series(fq_poly_t Q, const fq_poly_t A, const fq_poly_t B, slong n, const fq_ctx_t ctx)

Set $Q$ to $A / B \bmod x^n$ for power series with $B(0) \ne 0$. The
underscore version requires $0 < \mathit{Alen}, \mathit{Blen} \le n$,
writes $n$ coefficients and does not allow aliasing.
`fq_poly_div_series` throws `FLINT_DIVZERO` if $B$ is zero modulo $x^n$
and `FLINT_ERROR` if its constant term is zero.
:::

::: function
void _fq_poly_sqrt_series(fq_struct * g, const fq_struct * h, slong n, fq_ctx_t ctx)
void fq_poly_sqrt_series(fq_poly_t g, const fq_poly_t h, slong n, fq_ctx_t ctx)
void _fq_poly_invsqrt_series(fq_struct * g, const fq_struct * h, slong n, fq_ctx_t ctx)
void fq_poly_invsqrt_series(fq_poly_t g, const fq_poly_t h, slong n, fq_ctx_t ctx)

Set $g$ to $\sqrt{h}$, respectively $1/\sqrt{h}$, modulo $x^n$, for a
power series $h$ with constant term 1 (the square root with constant term
1), by Newton iteration. The underscore versions require $n > 0$ and $h$
zero-padded to length $n$, and do not allow aliasing. All four throw
`FLINT_ERROR` in characteristic 2 for $n > 1$; the others throw if the
constant term of $h$ is not 1.
:::

::: function
int _fq_poly_sqrt(fq_struct * s, const fq_struct * p, slong len, fq_ctx_t ctx)
int fq_poly_sqrt(fq_poly_t b, const fq_poly_t a, fq_ctx_t ctx)

If the polynomial is the square of a polynomial, set the output to a
square root and return 1; otherwise return 0. The underscore version
writes `len / 2 + 1` coefficients.

::: details
In characteristic 2 a square is a polynomial in $x^2$, whose square root
is computed coefficientwise with `fq_sqrt`. Otherwise the valuation must
be even, the constant term must be a square in $\mathbb{F}_q$, and the
candidate square root is computed as a power series square root of
length `len / 2 + 1` and checked by squaring its high part.
:::
:::

## Derivative {#derivative}

::: function
void _fq_poly_derivative(fq_struct * rop, const fq_struct * op, slong len, const fq_ctx_t ctx)
void fq_poly_derivative(fq_poly_t rop, const fq_poly_t op, const fq_ctx_t ctx)

Set `rop` to the formal derivative of `op`. The underscore version
writes `len - 1` coefficients and allows aliasing.
:::

## Greatest common divisor {#gcd}

GCDs use the generic Euclidean algorithm below a cutoff and the
[half-GCD algorithm](algorithms/polynomials/half-gcd) above it. The
cutoffs depend on the size of $p$: for `fq` the Euclidean algorithm is
used if the shorter input has length below `FQ_POLY_SMALL_GCD_CUTOFF`
(80) for $p < 2^8$ and below `FQ_POLY_GCD_CUTOFF` (90) otherwise, and
the half-GCD recursion switches to the basecase below
`FQ_POLY_HGCD_CUTOFF` (30).

::: summary algorithms/polynomials/half-gcd :::

::: function
slong _fq_poly_gcd(fq_struct * G, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_ctx_t ctx)

Sets $G$ to a GCD of `(A, lenA)` and `(B, lenB)` and returns its length.
Requires `lenA` $\ge$ `lenB` $> 0$ and space for `lenB` coefficients in
$G$. The GCD is not made monic. This is a wrapper of `_gr_poly_gcd_euclidean`
and `_gr_poly_gcd_hgcd`.
:::

::: function
void fq_poly_gcd(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

Sets `rop` to the monic GCD of `op1` and `op2`. The GCD of two zero
polynomials is zero.
:::

::: function
slong _fq_poly_gcd_euclidean_f(fq_t f, fq_struct * G, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_ctx_t ctx)
void fq_poly_gcd_euclidean_f(fq_t f, fq_poly_t G, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

GCD by the Euclidean algorithm for a possibly reducible defining
polynomial: either set $f = 1$ and $G$ to the GCD (monic in the
non-underscore version) and return its length, or set $f$ to a nontrivial
factor of the defining polynomial when a leading coefficient is not
invertible, leaving $G$ undefined. The underscore version requires
`lenA` $\ge$ `lenB` $> 0$ and space for `lenB` coefficients.
:::

::: function
slong _fq_poly_xgcd(fq_struct * G, fq_struct * S, fq_struct * fq, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_ctx_t ctx)

Computes a GCD $G$ of $A$ and $B$ and cofactors $S$, $T$ (the argument
named `fq`) with $SA + TB = G$, and returns the length of $G$. Requires
`lenA` $\ge$ `lenB` $\ge 1$ and $(\mathit{lenA}, \mathit{lenB}) \ne
(1, 1)$, space for `lenB` coefficients in $G$, `lenB - 1` in $S$ and
`lenA - 1` in $T$, and no aliasing. $G$ is not made monic. This is a
wrapper of `_gr_poly_xgcd_euclidean` and `_gr_poly_xgcd_hgcd`.
:::

::: function
void fq_poly_xgcd(fq_poly_t G, fq_poly_t S, fq_poly_t fq, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

Computes the monic GCD $G$ of $A$ and $B$ and cofactors $S$, $T$ (the
argument named `fq`) with $SA + TB = G$, where $\deg S < \deg B -
\deg G$ and $\deg T < \deg A - \deg G$ when both are nonzero. The GCD of
two zero polynomials is zero, with zero cofactors. $G$, $S$ and $T$ must
be distinct objects.
:::

::: function
slong _fq_poly_xgcd_euclidean_f(fq_t f, fq_struct * G, fq_struct * S, fq_struct * fq, const fq_struct * A, slong lenA, const fq_struct * B, slong lenB, const fq_ctx_t ctx)
void fq_poly_xgcd_euclidean_f(fq_t f, fq_poly_t G, fq_poly_t S, fq_poly_t fq, const fq_poly_t A, const fq_poly_t B, const fq_ctx_t ctx)

Extended Euclidean algorithm for a possibly reducible defining
polynomial: either set $f = 1$ and compute $G$, $S$, $T$ as
`fq_poly_xgcd` (the underscore version as `_fq_poly_xgcd`), or set $f$
to a nontrivial factor of the defining polynomial, leaving the outputs
undefined.
:::

## Evaluation {#evaluation}

::: function
void _fq_poly_evaluate_fq(fq_t rop, const fq_struct * op, slong len, const fq_t a, const fq_ctx_t ctx)
void fq_poly_evaluate_fq(fq_t res, const fq_poly_t f, const fq_t a, const fq_ctx_t ctx)

Set the result to $f(a)$, by Horner's rule. Any length is allowed.
:::

::: function
void _fq_poly_evaluate_fq_vec(fq_struct * ys, const fq_struct * coeffs, slong len, const fq_struct * xs, slong n, const fq_ctx_t ctx)
void fq_poly_evaluate_fq_vec(fq_struct * ys, const fq_poly_t poly, const fq_struct * xs, slong n, const fq_ctx_t ctx)

Set `ys[i]` to $f(x_i)$ for the $n$ points `xs`, using `_vec_iter`
for polynomials of length below 32 and `_vec_fast` otherwise.
:::

::: function
void _fq_poly_evaluate_fq_vec_iter(fq_struct * ys, const fq_struct * coeffs, slong len, const fq_struct * xs, slong n, const fq_ctx_t ctx)
void fq_poly_evaluate_fq_vec_iter(fq_struct * ys, const fq_poly_t poly, const fq_struct * xs, slong n, const fq_ctx_t ctx)

Evaluate at the $n$ points one by one with Horner's rule, at the cost of
$O(n \cdot \mathit{len})$ operations.
:::

::: function
void _fq_poly_evaluate_fq_vec_fast(fq_struct * ys, const fq_struct * poly, slong plen, const fq_struct * xs, slong n, const fq_ctx_t ctx)
void fq_poly_evaluate_fq_vec_fast(fq_struct * ys, const fq_poly_t poly, const fq_struct * xs, slong n, const fq_ctx_t ctx)

Evaluate at the $n$ points by building the subproduct tree of the
points and reducing the polynomial down the tree
([](algorithms/polynomials/subproduct-tree)), at the cost of
$O(M(n) \log n)$ operations plus one initial division.
:::

::: function
void _fq_poly_evaluate_fq_vec_fast_precomp(fq_struct * vs, const fq_struct * poly, slong plen, fq_poly_struct * const * tree, slong len, const fq_ctx_t ctx)

Evaluates `(poly, plen)` at the `len` points whose subproduct tree
`tree` was built by `_fq_poly_tree_build`, writing the values to `vs`.
:::

::: function
fq_poly_struct ** _fq_poly_tree_alloc(slong len, const fq_ctx_t ctx)
void _fq_poly_tree_free(fq_poly_struct ** tree, slong len, const fq_ctx_t ctx)

Allocate, respectively free, a subproduct tree for `len` points, with
$\lceil \log_2 \mathit{len} \rceil + 1$ levels of initialised polynomials.
:::

::: function
void _fq_poly_tree_build(fq_poly_struct ** tree, const fq_struct * roots, slong len, const fq_ctx_t ctx)

Fills the subproduct tree for the `len` points `roots`: level 0 holds the
linear factors $x - a_i$, and each higher level the products of pairs of
the level below.
:::

## Composition {#composition}

Composition and modular composition are wrappers of the generic
implementations of `gr_poly`, which use Horner's rule or the
baby-step giant-step algorithm of Brent and Kung [@BrentKung1978] with
a matrix multiplication (`fq_mat_mul`); see
[](algorithms/polynomials/composition). Modular composition is the main
step of the Frobenius computations in factorisation.

::: summary algorithms/polynomials/composition :::

::: function
void _fq_poly_compose(fq_struct * rop, const fq_struct * op1, slong len1, const fq_struct * op2, slong len2, const fq_ctx_t ctx)
void fq_poly_compose(fq_poly_t rop, const fq_poly_t op1, const fq_poly_t op2, const fq_ctx_t ctx)

Set `rop` to the composition $g(h(x))$ of $g$ = `op1` and $h$ = `op2`.
The underscore version requires nonzero inputs, space for
$(\mathit{len1} - 1)(\mathit{len2} - 1) + 1$ coefficients and no
aliasing.
:::

::: function
void _fq_poly_compose_mod(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_ctx_t ctx)
void fq_poly_compose_mod(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t poly3, const fq_ctx_t ctx)

Set `res` to $f(g) \bmod h$ for $f$ = `poly1`, $g$ = `poly2` and
$h$ = `poly3` $\ne 0$, choosing between Horner's rule and the
Brent–Kung algorithm. The underscore version requires $g$ reduced modulo
$h$ and zero-padded to length `lenh - 1`, writes `lenh - 1`
coefficients and does not allow aliasing.
:::

::: function
void _fq_poly_compose_mod_preinv(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_struct * hinv, slong lenhinv, const fq_ctx_t ctx)
void fq_poly_compose_mod_preinv(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t poly3, const fq_poly_t poly3inv, const fq_ctx_t ctx)

As `fq_poly_compose_mod`, given the inverse of the reversal of $h$.
:::

::: function
void _fq_poly_compose_mod_horner(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_ctx_t ctx)
void fq_poly_compose_mod_horner(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t poly3, const fq_ctx_t ctx)
void _fq_poly_compose_mod_horner_preinv(fq_struct * res, const fq_struct * f, slong lenf, const fq_struct * g, const fq_struct * h, slong lenh, const fq_struct * hinv, slong lenhinv, const fq_ctx_t ctx)
void fq_poly_compose_mod_horner_preinv(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t poly3, const fq_poly_t poly3inv, const fq_ctx_t ctx)

Modular composition by Horner's rule, with $\deg f$ multiplications
modulo $h$. Requirements as for `fq_poly_compose_mod`; the `_preinv`
versions also require $\mathit{lenf} < \mathit{lenh}$ and the inverse of
the reversal of $h$.
:::

::: function
void _fq_poly_compose_mod_brent_kung(fq_struct * res, const fq_struct * poly1, slong len1, const fq_struct * poly2, const fq_struct * poly3, slong len3, const fq_ctx_t ctx)
void fq_poly_compose_mod_brent_kung(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t poly3, const fq_ctx_t ctx)
void _fq_poly_compose_mod_brent_kung_preinv(fq_struct * res, const fq_struct * poly1, slong len1, const fq_struct * poly2, const fq_struct * poly3, slong len3, const fq_struct * poly3inv, slong len3inv, const fq_ctx_t ctx)
void fq_poly_compose_mod_brent_kung_preinv(fq_poly_t res, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t poly3, const fq_poly_t poly3inv, const fq_ctx_t ctx)

Modular composition $f(g) \bmod h$ by the Brent–Kung algorithm: the
powers $g^i \bmod h$ for $i < m = \lceil \sqrt{\mathit{lenf}} \rceil$
are computed, $f$ is split into $m$ blocks, the block compositions are
obtained by one $m \times m$ by $m \times \deg h$ matrix product, and
the results are combined by Horner's rule in $g^m$. Requires
$\deg f < \deg h$ and $g$ reduced modulo $h$; the `_preinv` versions also
require the inverse of the reversal of $h$. No aliasing in the
underscore versions.
:::

::: function
void _fq_poly_precompute_matrix(fq_mat_t A, const fq_struct * poly1, const fq_struct * poly2, slong len2, const fq_struct * poly2inv, slong len2inv, const fq_ctx_t ctx)
void fq_poly_precompute_matrix(fq_mat_t A, const fq_poly_t poly1, const fq_poly_t poly2, const fq_poly_t poly2inv, const fq_ctx_t ctx)

Set the rows of the matrix $A$ to the coefficients of the powers
$g^i \bmod h$, $0 \le i < m$, for $g$ = `poly1` and $h$ = `poly2`, as
needed by `fq_poly_compose_mod_brent_kung_precomp_preinv`. $A$ must have
exactly $m = \lfloor \sqrt{\deg h} \rfloor + 1$ rows and $\deg h$
columns (otherwise `fq_poly_precompute_matrix` throws `FLINT_ERROR`), and
`poly2inv` must be the inverse of the reversal of $h$. In the underscore
version $g$ must be reduced and zero-padded to length `len2 - 1`.
:::

::: function
void _fq_poly_compose_mod_brent_kung_precomp_preinv(fq_struct * res, const fq_struct * poly1, slong len1, const fq_mat_t A, const fq_struct * poly3, slong len3, const fq_struct * poly3inv, slong len3inv, const fq_ctx_t ctx)
void fq_poly_compose_mod_brent_kung_precomp_preinv(fq_poly_t res, const fq_poly_t poly1, const fq_mat_t A, const fq_poly_t poly3, const fq_poly_t poly3inv, const fq_ctx_t ctx)

Set `res` to $f(g) \bmod h$ for $f$ = `poly1` with $\deg f < \deg h$,
using the matrix $A$ of powers of $g$ computed by
`fq_poly_precompute_matrix`. This amortises the precomputation over many
compositions with the same $g$, as in the iterated Frobenius
computations of factorisation.
:::

::: function
void _fq_poly_reduce_matrix_mod_poly(fq_mat_t A, const fq_mat_t B, const fq_poly_t f, const fq_ctx_t ctx)

Sets the rows of $A$ to the rows of $B$, considered as polynomials,
reduced modulo $f$; $A$ must have $\deg f$ columns. Used for the
precomputation of the Brent–Kung matrix.
:::

## Characteristic and minimal polynomial of a matrix {#charpoly}

::: function
void fq_mat_charpoly(fq_poly_t p, const fq_mat_t M, const fq_ctx_t ctx)

Sets $p$ to the characteristic polynomial $\det(xI - M)$ of the square
matrix $M$; throws `FLINT_ERROR` if $M$ is not square. This is a wrapper
of `gr_mat_charpoly`, which uses the Berkowitz algorithm for
$n \le 12$ and Danilevsky's algorithm otherwise (for `fq_zech`, from
$n = 5$); see [](algorithms/linear-algebra/charpoly).
:::

::: function
void fq_mat_minpoly(fq_poly_t p, const fq_mat_t X, const fq_ctx_t ctx)

Sets $p$ to the minimal polynomial of the square matrix $X$, the monic
polynomial of least degree with $p(X) = 0$; throws `FLINT_ERROR` if $X$
is not square. This is a wrapper of `gr_mat_minpoly_field`.
:::

## Input and output {#io}

::: function
int _fq_poly_fprint(FILE * file, const fq_struct * poly, slong len, const fq_ctx_t ctx)
int fq_poly_fprint(FILE * file, const fq_poly_t poly, const fq_ctx_t ctx)
int _fq_poly_print(const fq_struct * poly, slong len, const fq_ctx_t ctx)
int fq_poly_print(const fq_poly_t poly, const fq_ctx_t ctx)

Print the polynomial to `file`, respectively `stdout`: the length, two
spaces, and the coefficients in the format of `fq_fprint` separated by
single spaces. Return a positive value on success and a non-positive value
on failure.
:::

::: function
int _fq_poly_fprint_pretty(FILE * file, const fq_struct * poly, slong len, const char * x, const fq_ctx_t ctx)
int fq_poly_fprint_pretty(FILE * file, const fq_poly_t poly, const char * x, const fq_ctx_t ctx)
int _fq_poly_print_pretty(const fq_struct * poly, slong len, const char * x, const fq_ctx_t ctx)
int fq_poly_print_pretty(const fq_poly_t poly, const char * x, const fq_ctx_t ctx)

Print the polynomial in the variable `x`, with coefficients other than 1
in parentheses in the format of `fq_fprint_pretty`, for example
`x^2+(a+1)*x+(a)`. Return a positive value on success.
:::

::: function
char * _fq_poly_get_str(const fq_struct * poly, slong len, const fq_ctx_t ctx)
char * fq_poly_get_str(const fq_poly_t poly, const fq_ctx_t ctx)
char * _fq_poly_get_str_pretty(const fq_struct * poly, slong len, const char * x, const fq_ctx_t ctx)
char * fq_poly_get_str_pretty(const fq_poly_t poly, const char * x, const fq_ctx_t ctx)

Return the plain, respectively pretty, representation as a newly
allocated string, which the caller must free with `flint_free`.
:::
