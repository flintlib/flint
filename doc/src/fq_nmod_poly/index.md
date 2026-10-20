---
title: fq_nmod_poly.h
header: fq_nmod_poly.h
template-headers: fq_poly_templates
description: Dense univariate polynomials over finite fields $\mathbb{F}_q$ in the word-size characteristic.
category: Finite fields
order: 14
---

This module implements dense univariate polynomials over a finite field $\mathbb{F}_q$ in the representation of [](fq_nmod) (a prime $p < 2^{\bits}$, elements stored as `nmod_poly`). It is generated from the same template as [](fq_poly) and has the same functions, with `fq_nmod` in place of `fq`; each function links to the `fq_poly` version, which has the details. As there, division, GCD, composition and power series functions are wrappers of the generic `gr_poly` algorithms.

The differences are in multiplication: `fq_nmod_poly_mul` uses the substitution $y = x^{2d-1}$ of the generator, which turns the product into one product in $(\Z/p\Z)[x]$ computed by `nmod_poly_mul` (with `fft_small` for long products), except for very short products (see [](fq_poly/mul)). The GCD functions switch from the Euclidean algorithm to the half-GCD algorithm at length 110 for $p < 2^8$ (`FQ_NMOD_POLY_SMALL_GCD_CUTOFF`) and 120 otherwise (`FQ_NMOD_POLY_GCD_CUTOFF`), with the half-GCD basecase below 25 (`FQ_NMOD_POLY_HGCD_CUTOFF`).

::: summary algorithms/polynomials/multiplication :::

## Types and conventions {#types}

::: type
fq_nmod_poly_struct
fq_nmod_poly_t

A polynomial over $\mathbb{F}_q$: an array `coeffs` of `alloc`
initialised elements (`fq_nmod_struct`), of which the first `length` are the
coefficients, the coefficient of $x^i$ at index $i$. `fq_nmod_poly_t` is an
array of one structure. The context is not stored and must be passed to
every function.
:::

## Memory management {#memory}

::: function
void fq_nmod_poly_init(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_init2(fq_nmod_poly_t poly, slong alloc, const fq_nmod_ctx_t ctx)

Initialise `poly` as the zero polynomial. `fq_nmod_poly_init2` allocates space
for `alloc` coefficients.

::: details
Generated from the same template as `fq_poly_init`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_realloc(fq_nmod_poly_t poly, slong alloc, const fq_nmod_ctx_t ctx)

Changes the allocation of `poly` to `alloc` coefficients, truncating and
normalising the polynomial if it is longer. For `alloc` = 0 the
coefficients are freed.

::: details
Generated from the same template as `fq_poly_realloc`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_fit_length(fq_nmod_poly_t poly, slong len, const fq_nmod_ctx_t ctx)

Ensures that `poly` has space for at least `len` coefficients, growing
the allocation to at least twice its old size if it must grow. The
polynomial is not changed.

::: details
Generated from the same template as `fq_poly_fit_length`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_clear(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)

Frees the memory used by `poly`.

::: details
Generated from the same template as `fq_poly_clear`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_set_length(fq_nmod_poly_t poly, slong len, const fq_nmod_ctx_t ctx)

Sets the length of `poly` to `len`, setting the coefficients from index
`len` up to the old length to zero. Does not normalise. `poly` must have
space for `len` coefficients.

::: details
Generated from the same template as `_fq_poly_set_length`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_normalise(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_normalise2(const fq_nmod_struct * poly, slong * length, const fq_nmod_ctx_t ctx)

Decrease the length of `poly`, respectively the variable `*length` of
the array `(poly, *length)`, until the leading coefficient is nonzero or
the length is zero.

::: details
Generated from the same template as `_fq_poly_normalise`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_truncate(fq_nmod_poly_t poly, slong len, const fq_nmod_ctx_t ctx)

Truncates `poly` to length at most `len` (that is, reduces it modulo
$x^{\mathit{len}}$) and normalises it.

::: details
Generated from the same template as `fq_poly_truncate`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_set_trunc(fq_nmod_poly_t poly1, fq_nmod_poly_t poly2, slong len, const fq_nmod_ctx_t ctx)

Sets `poly1` to `poly2` truncated to length `len`. `poly2` is not modified
unless it is the same object as `poly1`.

::: details
Generated from the same template as `fq_poly_set_trunc`, whose documentation has the details.
:::
:::

## Length, degree and comparison {#properties}

::: function
slong fq_nmod_poly_length(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
slong fq_nmod_poly_degree(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)

Return the length, respectively the degree, of `poly`; the zero
polynomial has length 0 and degree $-1$.

::: details
Generated from the same template as `fq_poly_length`, whose documentation has the details.
:::
:::

::: function
fq_nmod_struct * fq_nmod_poly_lead(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)

Returns a pointer to the leading coefficient of `poly`, or `NULL` for the
zero polynomial.

::: details
Generated from the same template as `fq_poly_lead`, whose documentation has the details.
:::
:::

::: function
slong _fq_nmod_poly_hamming_weight(const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
slong fq_nmod_poly_hamming_weight(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)

Return the number of nonzero coefficients.

::: details
Generated from the same template as `_fq_poly_hamming_weight`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_poly_equal(const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_ctx_t ctx)

Returns whether the two polynomials are equal.

::: details
Generated from the same template as `fq_poly_equal`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_poly_equal_trunc(const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, slong n, const fq_nmod_ctx_t ctx)

Returns whether the two polynomials are equal modulo $x^n$.

::: details
Generated from the same template as `fq_poly_equal_trunc`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_poly_equal_fq_nmod(const fq_nmod_poly_t poly, const fq_nmod_t c, const fq_nmod_ctx_t ctx)

Returns whether `poly` is the constant polynomial $c$.

::: details
Generated from the same template as `fq_poly_equal_fq`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_poly_is_zero(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_is_one(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_is_gen(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)

Return whether the polynomial is 0, 1, respectively $x$.

::: details
Generated from the same template as `fq_poly_is_zero`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_poly_is_unit(const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)

Returns whether `op` is a unit of $\mathbb{F}_q[x]$, that is, a nonzero
constant.

::: details
Generated from the same template as `fq_poly_is_unit`, whose documentation has the details.
:::
:::

## Assignment and conversion {#assignment}

::: function
void _fq_nmod_poly_set(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_set(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)

Copy `op` into `rop`.

::: details
Generated from the same template as `_fq_poly_set`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_set_fq_nmod(fq_nmod_poly_t poly, const fq_nmod_t c, const fq_nmod_ctx_t ctx)

Sets `poly` to the constant polynomial $c$.

::: details
Generated from the same template as `fq_poly_set_fq`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_set_fmpz_mod_poly(fq_nmod_poly_t rop, const fmpz_mod_poly_t op, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_set_nmod_poly(fq_nmod_poly_t rop, const nmod_poly_t op, const fq_nmod_ctx_t ctx)

Set `rop` to the polynomial `op` over $\mathbb{F}_p$, with coefficients
mapped to $\mathbb{F}_p \subseteq \mathbb{F}_q$. The modulus of `op` must
be $p$; it is not checked.

::: details
Generated from the same template as `fq_poly_set_fmpz_mod_poly`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_swap(fq_nmod_poly_t op1, fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)

Swaps `op1` and `op2` efficiently.

::: details
Generated from the same template as `fq_poly_swap`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_zero(fq_nmod_struct * rop, slong len, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_zero(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)

Set the `len` coefficients of `rop`, respectively the polynomial `poly`,
to zero.

::: details
Generated from the same template as `_fq_poly_zero`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_one(fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_gen(fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)

Set the polynomial to 1, respectively to $x$.

::: details
Generated from the same template as `fq_poly_one`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_make_monic(fq_nmod_struct * rop, const fq_nmod_struct * op, slong length, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_make_monic(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)

Set `rop` to `op` divided by its leading coefficient. The zero polynomial
is mapped to zero; the underscore version requires a nonzero leading
coefficient `op[length - 1]` and allows aliasing.

::: details
Generated from the same template as `_fq_poly_make_monic`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_reverse(fq_nmod_struct * res, const fq_nmod_struct * poly, slong len, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_reverse(fq_nmod_poly_t res, const fq_nmod_poly_t poly, slong n, const fq_nmod_ctx_t ctx)

Set `res` to the reversal $x^{n-1} f(1/x)$ of the polynomial $f$ =
`poly` truncated to length $n$, considered as a polynomial of length $n$.
The underscore version requires $\mathit{len} \le n$, writes $n$
coefficients and allows `res` to be `poly`.

::: details
Generated from the same template as `_fq_poly_reverse`, whose documentation has the details.
:::
:::

::: function
ulong fq_nmod_poly_deflation(const fq_nmod_poly_t input, const fq_nmod_ctx_t ctx)

Returns the largest $k$ such that `input` is a polynomial in $x^k$. For
a polynomial of length at most 1 the result is its length.

::: details
Generated from the same template as `fq_poly_deflation`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_deflate(fq_nmod_poly_t result, const fq_nmod_poly_t input, ulong deflation, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_inflate(fq_nmod_poly_t result, const fq_nmod_poly_t input, ulong inflation, const fq_nmod_ctx_t ctx)

`fq_nmod_poly_deflate` sets `result` to $g$ with $g(x^k) = f$ for the
polynomial $f$ = `input` and $k$ = `deflation`, which must divide the
exponents of all nonzero terms of $f$; it throws `FLINT_DIVZERO` for
$k = 0$. `fq_nmod_poly_inflate` sets `result` to $f(x^k)$ for
$k$ = `inflation`.

::: details
Generated from the same template as `fq_poly_deflate`, whose documentation has the details.
:::
:::

## Coefficients {#coefficients}

::: function
void fq_nmod_poly_get_coeff(fq_nmod_t x, const fq_nmod_poly_t poly, slong n, const fq_nmod_ctx_t ctx)

Sets $x$ to the coefficient of $x^n$ in `poly`, which is zero for $n$ at
least the length.

::: details
Generated from the same template as `fq_poly_get_coeff`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_set_coeff(fq_nmod_poly_t poly, slong n, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_set_coeff_fmpz(fq_nmod_poly_t poly, slong n, const fmpz_t x, const fq_nmod_ctx_t ctx)

Set the coefficient of $x^n$, $n \ge 0$, to the element $x$, respectively
to the image of the integer $x$ in $\mathbb{F}_p$, extending and
normalising the polynomial as needed.

::: details
Generated from the same template as `fq_poly_set_coeff`, whose documentation has the details.
:::
:::

## Random polynomials {#random}

::: function
void fq_nmod_poly_randtest(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_randtest_not_zero(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)

Set $f$ to a random polynomial of length at most `len`, with coefficients
from `fq_nmod_randtest`; the second version returns a nonzero polynomial and
requires `len` $> 0$ (it throws otherwise).

::: details
Generated from the same template as `fq_poly_randtest`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_randtest_monic(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)

Sets $f$ to a random monic polynomial of length `len` $\ge 1$.

::: details
Generated from the same template as `fq_poly_randtest_monic`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_randtest_irreducible(fq_nmod_poly_t f, flint_rand_t state, slong len, const fq_nmod_ctx_t ctx)

Sets $f$ to a random monic irreducible polynomial of length `len`
$\ge 2$. Random monic polynomials are tested by computing
$\gcd(f, x^{q^i} - x)$ for $i \le (\mathit{len} - 1)/2$ until one is
irreducible, which takes about `len` attempts on average.

::: details
Generated from the same template as `fq_poly_randtest_irreducible`, whose documentation has the details.
:::
:::

## Addition and subtraction {#addition}

::: function
void _fq_nmod_poly_add(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_sub(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, const fq_nmod_ctx_t ctx)

Set `(res, max(len1, len2))` to the sum, respectively the difference, of
the two arrays. Aliasing of `res` with the inputs is allowed.

::: details
Generated from the same template as `_fq_poly_add`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_add(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_sub(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op1} + \mathit{op2}$, respectively
$\mathit{op1} - \mathit{op2}$.

::: details
Generated from the same template as `fq_poly_add`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_add_si(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, slong c, const fq_nmod_ctx_t ctx)

Sets `rop` to $\mathit{op1} + c$, where $c$ is mapped to $\mathbb{F}_p$.

::: details
Generated from the same template as `fq_poly_add_si`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_add_series(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_sub_series(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)

Set `rop` to the sum, respectively the difference, of the two
polynomials truncated to length $n$.

::: details
Generated from the same template as `fq_poly_add_series`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_neg(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_neg(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)

Set `rop` to $-\mathit{op}$.

::: details
Generated from the same template as `_fq_poly_neg`, whose documentation has the details.
:::
:::

## Scalar multiplication and division {#scalar}

::: function
void _fq_nmod_poly_scalar_mul_fq_nmod(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_scalar_mul_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)

Set `rop` to $x \cdot \mathit{op}$.

::: details
Generated from the same template as `_fq_poly_scalar_mul_fq`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_scalar_div_fq_nmod(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_scalar_div_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op} / x$. Throws `FLINT_ERROR` if $x$ is zero.

::: details
Generated from the same template as `_fq_poly_scalar_div_fq`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_scalar_addmul_fq_nmod(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_scalar_addmul_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_scalar_submul_fq_nmod(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_t x, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_scalar_submul_fq_nmod(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_t x, const fq_nmod_ctx_t ctx)

Add $x \cdot \mathit{op}$ to `rop`, respectively subtract it. The
underscore versions work on the first `len` coefficients of both arrays.

::: details
Generated from the same template as `_fq_poly_scalar_addmul_fq`, whose documentation has the details.
:::
:::

## Multiplication {#multiplication}

::: function
void _fq_nmod_poly_mul(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mul(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op1} \cdot \mathit{op2}$, choosing the algorithm.
The underscore version writes `len1 + len2 - 1` coefficients, requires
`len1`, `len2` $\ge 1$ and does not allow aliasing.

::: details
Generated from the same template as `_fq_poly_mul`, whose documentation has the details. The algorithm selection of each representation is described in [](fq_poly/mul).
:::
:::

::: function
void _fq_nmod_poly_mul_classical(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mul_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)

Set `rop` to the product, computed by the classical algorithm with one
dot product (reduced once) per coefficient (a wrapper of
`_gr_poly_mullow_classical`). The underscore version requires
`len1`, `len2` $\ge 1$ and does not allow aliasing.

::: details
Generated from the same template as `_fq_poly_mul_classical`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_mul_KS(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mul_KS(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)

Set `rop` to the product, computed by Kronecker substitution: every
coefficient is packed into one integer, the integer polynomials are
multiplied with `_fmpz_poly_mul`, and the coefficients are unpacked and
reduced. The underscore version allows zero padding of the inputs (it
writes `len1 + len2 - 1` coefficients) and aliasing.

::: details
Generated from the same template as `_fq_poly_mul_KS`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_mul_univariate(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mul_univariate(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)

Set `rop` to the product, computed by substituting $y = x^{2d-1}$ for the
generator of $\mathbb{F}_q$, which gives one product in $(\Z/p\Z)[x]$ of
length about $(2d - 1)\,\mathit{len}$ computed with `_nmod_poly_mul`, and
reducing each block of $2d - 1$ coefficients modulo the defining
polynomial. The underscore version requires `len1`, `len2` $\ge 1$ and
writes `len1 + len2 - 1` coefficients.
:::

::: function
void _fq_nmod_poly_mullow(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mullow(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)

Set `rop` to the product truncated to length $n$, that is, the low $n$
coefficients of $\mathit{op1} \cdot \mathit{op2}$. The underscore version
requires $0 < n \le \mathit{len1} + \mathit{len2} - 1$ and
`len1`, `len2` $\ge 1$, writes $n$ coefficients and does not allow
aliasing. The classical algorithm is used for $n$ below
`FQ_NMOD_MULLOW_CLASSICAL_CUTOFF` (6) or short inputs, and Kronecker
substitution otherwise.

::: details
Generated from the same template as `_fq_poly_mullow`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_mullow_classical(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mullow_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_mullow_KS(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mullow_KS(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_mullow_univariate(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mullow_univariate(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong n, const fq_nmod_ctx_t ctx)

Truncated products with a fixed algorithm: classical, Kronecker
substitution, or the substitution $y = x^{2d-1}$. The underscore versions
require $0 < n \le \mathit{len1} + \mathit{len2} - 1$ and write $n$
coefficients; `_fq_nmod_poly_mullow_KS` allows zero padding and zero inputs.

::: details
Generated from the same template as `_fq_poly_mullow_classical`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_mulhigh(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, slong n, fq_nmod_ctx_t ctx)
void fq_nmod_poly_mulhigh(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong start, const fq_nmod_ctx_t ctx)

Compute the coefficients of index at least $n$ (`start`) of the product.
`fq_nmod_poly_mulhigh` sets `rop` to a polynomial of length
`len1 + len2 - 1` whose coefficients from `start` on are those of the
product and whose lower coefficients are zero; it uses the classical
algorithm. `_fq_nmod_poly_mulhigh` requires `len1` $\ge$ `len2` $\ge 1$ and
space for `len1 + len2 - 1` coefficients in `res`; its coefficients below
$n$ are arbitrary (it computes the full product by Kronecker substitution
for inputs of length at least 6).

::: details
Generated from the same template as `_fq_poly_mulhigh`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_mulhigh_classical(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, slong start, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mulhigh_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, slong start, const fq_nmod_ctx_t ctx)

As the previous functions, with the classical algorithm, which computes
only the coefficients from `start` on and sets the lower ones to zero.
The underscore version requires `len1` $\ge$ `len2` $\ge 1$ and does not
allow aliasing.

::: details
Generated from the same template as `_fq_poly_mulhigh_classical`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_mulmod(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, const fq_nmod_struct * f, slong lenf, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mulmod(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)

Set `res` to $\mathit{poly1} \cdot \mathit{poly2} \bmod f$. The underscore
version requires $\mathit{len1} + \mathit{len2} > \mathit{lenf}$ (the
product must need reduction), writes $\mathit{lenf} - 1$ coefficients
and needs space for `len1 + len2 - 1` coefficients in `res`, which must
not be aliased with $f$. `fq_nmod_poly_mulmod` throws `FLINT_DIVZERO` if $f = 0$.

::: details
Generated from the same template as `_fq_poly_mulmod`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_mulmod_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, slong len2, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_mulmod_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)

As `fq_nmod_poly_mulmod`, with the precomputed inverse `finv` of the reversal
of $f$ modulo $x^{\deg f + 1}$ (`fq_nmod_poly_reverse` followed by
`fq_nmod_poly_inv_series`). The product must satisfy
$\mathit{len1} + \mathit{len2} - 1 \le 2\,\mathit{lenf} - 2$; the inputs
need not be reduced. The underscore version does not allow aliasing.

::: details
Generated from the same template as `_fq_poly_mulmod_preinv`, whose documentation has the details.
:::
:::

## Squaring {#squaring}

::: function
void _fq_nmod_poly_sqr(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_sqr(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op}^2$, choosing the algorithm as for
multiplication. The underscore version writes `2 * len - 1`
coefficients, requires `len` $\ge 1$ and does not allow aliasing.

::: details
Generated from the same template as `_fq_poly_sqr`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_sqr_classical(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_sqr_classical(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_sqr_KS(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_sqr_KS(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)

Squaring with a fixed algorithm: classical or Kronecker substitution. The
underscore versions write `2 * len - 1` coefficients; the classical
version does not allow aliasing.
:::

## Powering {#powering}

::: function
void _fq_nmod_poly_pow(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, ulong e, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_pow(fq_nmod_poly_t rop, const fq_nmod_poly_t op, ulong e, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op}^e$, with $\mathit{op}^0 = 1$. The underscore
version requires $e > 0$, `len` $> 0$ and space for
$e(\mathit{len} - 1) + 1$ coefficients, and does not allow aliasing. It
uses binary exponentiation (`_gr_poly_pow_ui_binexp`).

::: details
Generated from the same template as `_fq_poly_pow`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_pow_trunc(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_pow_trunc(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_pow_trunc_binexp(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_pow_trunc_binexp(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, slong trunc, const fq_nmod_ctx_t ctx)

Set `res` to $\mathit{poly}^e \bmod x^{\mathit{trunc}}$, by binary
exponentiation with truncated multiplications. The underscore versions
require the input to be zero-padded to length `trunc`, `trunc` $> 0$,
$e > 1$ and space for `trunc` coefficients, and do not allow aliasing.

::: details
Generated from the same template as `_fq_poly_pow_trunc`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_powmod_ui_binexp(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, const fq_nmod_struct * f, slong lenf, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_powmod_ui_binexp(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_powmod_fmpz_binexp(fq_nmod_struct * res, const fq_nmod_struct * poly, const fmpz_t e, const fq_nmod_struct * f, slong lenf, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_powmod_fmpz_binexp(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fmpz_t e, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)

Set `res` to $\mathit{poly}^e \bmod f$ for $e \ge 0$, by binary
exponentiation. The underscore versions require $e > 0$,
`lenf` $> 1$, `poly` reduced modulo $f$ and zero-padded to length
`lenf - 1`, and space for `lenf - 1` coefficients in `res`.

::: details
Generated from the same template as `_fq_poly_powmod_ui_binexp`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_powmod_ui_binexp_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly, ulong e, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_powmod_ui_binexp_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly, ulong e, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_powmod_fmpz_binexp_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly, const fmpz_t e, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_powmod_fmpz_binexp_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fmpz_t e, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)

As the previous functions, with the precomputed inverse `finv` of the
reversal of $f$, so that every reduction is a Newton division with two
multiplications.

::: details
Generated from the same template as `_fq_poly_powmod_ui_binexp_preinv`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_powmod_fmpz_sliding_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly, const fmpz_t e, ulong k, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_powmod_fmpz_sliding_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly, const fmpz_t e, ulong k, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)

As `fq_nmod_poly_powmod_fmpz_binexp_preinv`, using sliding window
exponentiation with windows of $k$ bits; $k = 0$ chooses $k$ from the
size of $e$.

::: details
Generated from the same template as `_fq_poly_powmod_fmpz_sliding_preinv`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_powmod_x_fmpz_preinv(fq_nmod_struct * res, const fmpz_t e, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * finv, slong lenfinv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_powmod_x_fmpz_preinv(fq_nmod_poly_t res, const fmpz_t e, const fq_nmod_poly_t f, const fq_nmod_poly_t finv, const fq_nmod_ctx_t ctx)

Set `res` to $x^e \bmod f$, using sliding window exponentiation, in
which multiplications by $x$ are shifts. The underscore version requires
$e > 0$, `lenf` $> 2$ and space for `lenf - 1` coefficients. This is the
first step of distinct-degree factorisation and irreducibility testing,
with $e = q$.

::: details
Generated from the same template as `_fq_poly_powmod_x_fmpz_preinv`, whose documentation has the details.
:::
:::

## Shifting {#shifting}

::: function
void _fq_nmod_poly_shift_left(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_shift_left(fq_nmod_poly_t rop, const fq_nmod_poly_t op, slong n, const fq_nmod_ctx_t ctx)

Set `rop` to $x^n \cdot \mathit{op}$. The underscore version writes
`len + n` coefficients and allows aliasing.

::: details
Generated from the same template as `_fq_poly_shift_left`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_shift_right(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_shift_right(fq_nmod_poly_t rop, const fq_nmod_poly_t op, slong n, const fq_nmod_ctx_t ctx)

Set `rop` to $\lfloor \mathit{op} / x^n \rfloor$, dropping the $n$
lowest coefficients; the result is zero if $n$ is at least the length.
The underscore version requires `len` $> n$, writes `len - n`
coefficients and allows aliasing (the top $n$ coefficients of `op` are
then permuted, not zeroed).

::: details
Generated from the same template as `_fq_poly_shift_right`, whose documentation has the details.
:::
:::

## Euclidean division {#division}

::: summary algorithms/polynomials/newton :::

::: function
void _fq_nmod_poly_divrem(fq_nmod_struct * Q, fq_nmod_struct * R, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)

Computes `(Q, lenA - lenB + 1)` and `(R, lenB - 1)` with $A = BQ + R$ and
$\deg R < \deg B$, given the inverse `invB` of the leading coefficient of
$B$. Requires `lenA` $\ge$ `lenB` $> 0$. Allows aliasing of $R$ and $A$
but no other aliasing.

::: details
Generated from the same template as `_fq_poly_divrem`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_divrem(fq_nmod_poly_t Q, fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

Computes $Q$ and $R$ with $A = BQ + R$ and $\deg R < \deg B$. $B$ must
be nonzero. If $\deg A < \deg B$, then $Q = 0$ and $R = A$.

::: details
Generated from the same template as `fq_poly_divrem`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_div(fq_nmod_struct * Q, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_div(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

Compute only the quotient $Q$ of the division of $A$ by $B$. The
underscore version requires `lenA` $\ge$ `lenB` $> 0$ and writes
`lenA - lenB + 1` coefficients; `fq_nmod_poly_div` requires $B \ne 0$.

::: details
Generated from the same template as `_fq_poly_div`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_div_basecase(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

The same as `fq_nmod_poly_div`; kept for compatibility with FLINT 2.

::: details
Generated from the same template as `fq_poly_div_basecase`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_rem(fq_nmod_struct * R, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_rem(fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

Compute only the remainder $R$ of the division of $A$ by $B \ne 0$. The
underscore version writes `lenB - 1` coefficients and also accepts
`lenA` $<$ `lenB`; `invB` is not used.

::: details
Generated from the same template as `_fq_poly_rem`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_divrem_f(fq_nmod_t f, fq_nmod_struct * Q, fq_nmod_struct * R, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_divrem_f(fq_nmod_t f, fq_nmod_poly_t Q, fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

As `fq_nmod_poly_divrem`, for a possibly reducible defining polynomial: if the
leading coefficient of $B$ is invertible, sets $f = 1$ and computes $Q$
and $R$; otherwise sets $f$ to a nontrivial factor of the defining
polynomial (see `fq_nmod_gcdinv`) and leaves $Q$ and $R$ unchanged. $B$ must
be nonzero.

::: details
Generated from the same template as `_fq_poly_divrem_f`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_div_newton_n_preinv(fq_nmod_struct * Q, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_struct * Binv, slong lenBinv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_div_newton_n_preinv(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_poly_t Binv, const fq_nmod_ctx_t ctx)

Compute the quotient $Q$ of $A$ by $B$ with two multiplications, given
the inverse `Binv` of the reversal of $B$ modulo $x^{\mathit{lenB}}$:
the reversed quotient is the product of the reversal of $A$ and `Binv`
modulo $x^{\mathit{lenA} - \mathit{lenB} + 1}$. Requires
$\mathit{lenA} \le 2\,\mathit{lenB} - 1$; the underscore version writes
`lenA - lenB + 1` coefficients and does not read $B$. The non-underscore
version throws `FLINT_DIVZERO` for $B = 0$.

::: details
Generated from the same template as `_fq_poly_div_newton_n_preinv`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_divrem_newton_n_preinv(fq_nmod_struct * Q, fq_nmod_struct * R, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_struct * Binv, slong lenBinv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_divrem_newton_n_preinv(fq_nmod_poly_t Q, fq_nmod_poly_t R, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_poly_t Binv, const fq_nmod_ctx_t ctx)

As the previous functions, also computing the remainder
$R = A - BQ$, with one more multiplication. Requires
$\mathit{lenA} \le 2\,\mathit{lenB} - 1$.

::: details
Generated from the same template as `_fq_poly_divrem_newton_n_preinv`, whose documentation has the details.
:::
:::

::: function
int _fq_nmod_poly_divides(fq_nmod_struct * Q, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_t invB, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_divides(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

Return 1 and set $Q = A / B$ if $B$ divides $A$, and return 0 otherwise
(then $Q$ is undefined). The underscore version requires
`lenA` $\ge$ `lenB` $> 0$, space for `lenA - lenB + 1` coefficients in
$Q$ and no aliasing. `fq_nmod_poly_divides` throws `FLINT_DIVZERO` if $B = 0$. The test performs a full division with remainder.

::: details
Generated from the same template as `_fq_poly_divides`, whose documentation has the details.
:::
:::

::: function
ulong fq_nmod_poly_remove(fq_nmod_poly_t f, const fq_nmod_poly_t g, const fq_nmod_ctx_t ctx)

Removes the highest power of $g$ that divides $f$ from $f$ and returns
its exponent. $g$ must have degree at least 1.

::: details
Generated from the same template as `fq_poly_remove`, whose documentation has the details.
:::
:::

## Power series {#series}

::: function
void _fq_nmod_poly_inv_series_newton(fq_nmod_struct * Qinv, const fq_nmod_struct * Q, slong n, const fq_nmod_t cinv, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_inv_series(fq_nmod_struct * Qinv, const fq_nmod_struct * Q, slong n, const fq_nmod_t cinv, const fq_nmod_ctx_t ctx)

Set `(Qinv, n)` to the inverse of the power series `(Q, n)` modulo
$x^n$, given the inverse `cinv` of the constant term. Requires $n > 0$
and no aliasing. The basecase algorithm is used for small $n$ and Newton
iteration (`_gr_poly_inv_series_newton`) otherwise; the two functions
are the same.

::: details
Generated from the same template as `_fq_poly_inv_series_newton`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_inv_series_newton(fq_nmod_poly_t Qinv, const fq_nmod_poly_t Q, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_inv_series(fq_nmod_poly_t Qinv, const fq_nmod_poly_t Q, slong n, const fq_nmod_ctx_t ctx)

Set `Qinv` to the inverse of the power series $Q$ modulo $x^n$. The
constant term of $Q$ must be nonzero (otherwise `FLINT_ERROR` is thrown)
and $n \ge 1$.

::: details
Generated from the same template as `fq_poly_inv_series_newton`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_div_series(fq_nmod_struct * Q, const fq_nmod_struct * A, slong Alen, const fq_nmod_struct * B, slong Blen, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_div_series(fq_nmod_poly_t Q, const fq_nmod_poly_t A, const fq_nmod_poly_t B, slong n, const fq_nmod_ctx_t ctx)

Set $Q$ to $A / B \bmod x^n$ for power series with $B(0) \ne 0$. The
underscore version requires $0 < \mathit{Alen}, \mathit{Blen} \le n$,
writes $n$ coefficients and does not allow aliasing.
`fq_nmod_poly_div_series` throws `FLINT_DIVZERO` if $B$ is zero modulo $x^n$
and `FLINT_ERROR` if its constant term is zero.

::: details
Generated from the same template as `_fq_poly_div_series`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_sqrt_series(fq_nmod_struct * g, const fq_nmod_struct * h, slong n, fq_nmod_ctx_t ctx)
void fq_nmod_poly_sqrt_series(fq_nmod_poly_t g, const fq_nmod_poly_t h, slong n, fq_nmod_ctx_t ctx)
void _fq_nmod_poly_invsqrt_series(fq_nmod_struct * g, const fq_nmod_struct * h, slong n, fq_nmod_ctx_t ctx)
void fq_nmod_poly_invsqrt_series(fq_nmod_poly_t g, const fq_nmod_poly_t h, slong n, fq_nmod_ctx_t ctx)

Set $g$ to $\sqrt{h}$, respectively $1/\sqrt{h}$, modulo $x^n$, for a
power series $h$ with constant term 1 (the square root with constant term
1), by Newton iteration. The underscore versions require $n > 0$ and $h$
zero-padded to length $n$, and do not allow aliasing. All four throw
`FLINT_ERROR` in characteristic 2 for $n > 1$; the others throw if the
constant term of $h$ is not 1.

::: details
Generated from the same template as `_fq_poly_sqrt_series`, whose documentation has the details.
:::
:::

::: function
int _fq_nmod_poly_sqrt(fq_nmod_struct * s, const fq_nmod_struct * p, slong len, fq_nmod_ctx_t ctx)
int fq_nmod_poly_sqrt(fq_nmod_poly_t b, const fq_nmod_poly_t a, fq_nmod_ctx_t ctx)

If the polynomial is the square of a polynomial, set the output to a
square root and return 1; otherwise return 0. The underscore version
writes `len / 2 + 1` coefficients.

::: details
Generated from the same template as `_fq_poly_sqrt`, whose documentation has the details.
:::
:::

## Derivative {#derivative}

::: function
void _fq_nmod_poly_derivative(fq_nmod_struct * rop, const fq_nmod_struct * op, slong len, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_derivative(fq_nmod_poly_t rop, const fq_nmod_poly_t op, const fq_nmod_ctx_t ctx)

Set `rop` to the formal derivative of `op`. The underscore version
writes `len - 1` coefficients and allows aliasing.

::: details
Generated from the same template as `_fq_poly_derivative`, whose documentation has the details.
:::
:::

## Greatest common divisor {#gcd}

::: summary algorithms/polynomials/half-gcd :::

::: function
slong _fq_nmod_poly_gcd(fq_nmod_struct * G, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_ctx_t ctx)

Sets $G$ to a GCD of `(A, lenA)` and `(B, lenB)` and returns its length.
Requires `lenA` $\ge$ `lenB` $> 0$ and space for `lenB` coefficients in
$G$. The GCD is not made monic. This is a wrapper of `_gr_poly_gcd_euclidean`
and `_gr_poly_gcd_hgcd`.

::: details
Generated from the same template as `_fq_poly_gcd`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_gcd(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)

Sets `rop` to the monic GCD of `op1` and `op2`. The GCD of two zero
polynomials is zero.

::: details
Generated from the same template as `fq_poly_gcd`, whose documentation has the details.
:::
:::

::: function
slong _fq_nmod_poly_gcd_euclidean_f(fq_nmod_t f, fq_nmod_struct * G, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_gcd_euclidean_f(fq_nmod_t f, fq_nmod_poly_t G, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

GCD by the Euclidean algorithm for a possibly reducible defining
polynomial: either set $f = 1$ and $G$ to the GCD (monic in the
non-underscore version) and return its length, or set $f$ to a nontrivial
factor of the defining polynomial when a leading coefficient is not
invertible, leaving $G$ undefined. The underscore version requires
`lenA` $\ge$ `lenB` $> 0$ and space for `lenB` coefficients.

::: details
Generated from the same template as `_fq_poly_gcd_euclidean_f`, whose documentation has the details.
:::
:::

::: function
slong _fq_nmod_poly_xgcd(fq_nmod_struct * G, fq_nmod_struct * S, fq_nmod_struct * fq_nmod, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_ctx_t ctx)

Computes a GCD $G$ of $A$ and $B$ and cofactors $S$, $T$ (the argument
named `fq_nmod`) with $SA + TB = G$, and returns the length of $G$. Requires
`lenA` $\ge$ `lenB` $\ge 1$ and $(\mathit{lenA}, \mathit{lenB}) \ne
(1, 1)$, space for `lenB` coefficients in $G$, `lenB - 1` in $S$ and
`lenA - 1` in $T$, and no aliasing. $G$ is not made monic. This is a
wrapper of `_gr_poly_xgcd_euclidean` and `_gr_poly_xgcd_hgcd`.

::: details
Generated from the same template as `_fq_poly_xgcd`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_xgcd(fq_nmod_poly_t G, fq_nmod_poly_t S, fq_nmod_poly_t fq_nmod, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

Computes the monic GCD $G$ of $A$ and $B$ and cofactors $S$, $T$ (the
argument named `fq_nmod`) with $SA + TB = G$, where $\deg S < \deg B -
\deg G$ and $\deg T < \deg A - \deg G$ when both are nonzero. The GCD of
two zero polynomials is zero, with zero cofactors. $G$, $S$ and $T$ must
be distinct objects.

::: details
Generated from the same template as `fq_poly_xgcd`, whose documentation has the details.
:::
:::

::: function
slong _fq_nmod_poly_xgcd_euclidean_f(fq_nmod_t f, fq_nmod_struct * G, fq_nmod_struct * S, fq_nmod_struct * fq_nmod, const fq_nmod_struct * A, slong lenA, const fq_nmod_struct * B, slong lenB, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_xgcd_euclidean_f(fq_nmod_t f, fq_nmod_poly_t G, fq_nmod_poly_t S, fq_nmod_poly_t fq_nmod, const fq_nmod_poly_t A, const fq_nmod_poly_t B, const fq_nmod_ctx_t ctx)

Extended Euclidean algorithm for a possibly reducible defining
polynomial: either set $f = 1$ and compute $G$, $S$, $T$ as
`fq_nmod_poly_xgcd` (the underscore version as `_fq_nmod_poly_xgcd`), or set $f$
to a nontrivial factor of the defining polynomial, leaving the outputs
undefined.

::: details
Generated from the same template as `_fq_poly_xgcd_euclidean_f`, whose documentation has the details.
:::
:::

## Evaluation {#evaluation}

::: function
void _fq_nmod_poly_evaluate_fq_nmod(fq_nmod_t rop, const fq_nmod_struct * op, slong len, const fq_nmod_t a, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_evaluate_fq_nmod(fq_nmod_t res, const fq_nmod_poly_t f, const fq_nmod_t a, const fq_nmod_ctx_t ctx)

Set the result to $f(a)$, by Horner's rule. Any length is allowed.

::: details
Generated from the same template as `_fq_poly_evaluate_fq`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_evaluate_fq_nmod_vec(fq_nmod_struct * ys, const fq_nmod_struct * coeffs, slong len, const fq_nmod_struct * xs, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_evaluate_fq_nmod_vec(fq_nmod_struct * ys, const fq_nmod_poly_t poly, const fq_nmod_struct * xs, slong n, const fq_nmod_ctx_t ctx)

Set `ys[i]` to $f(x_i)$ for the $n$ points `xs`, using `_vec_iter`
for polynomials of length below 32 and `_vec_fast` otherwise.

::: details
Generated from the same template as `_fq_poly_evaluate_fq_vec`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_evaluate_fq_nmod_vec_iter(fq_nmod_struct * ys, const fq_nmod_struct * coeffs, slong len, const fq_nmod_struct * xs, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_evaluate_fq_nmod_vec_iter(fq_nmod_struct * ys, const fq_nmod_poly_t poly, const fq_nmod_struct * xs, slong n, const fq_nmod_ctx_t ctx)

Evaluate at the $n$ points one by one with Horner's rule, at the cost of
$O(n \cdot \mathit{len})$ operations.

::: details
Generated from the same template as `_fq_poly_evaluate_fq_vec_iter`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_evaluate_fq_nmod_vec_fast(fq_nmod_struct * ys, const fq_nmod_struct * poly, slong plen, const fq_nmod_struct * xs, slong n, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_evaluate_fq_nmod_vec_fast(fq_nmod_struct * ys, const fq_nmod_poly_t poly, const fq_nmod_struct * xs, slong n, const fq_nmod_ctx_t ctx)

Evaluate at the $n$ points by building the subproduct tree of the
points and reducing the polynomial down the tree
([](algorithms/polynomials/subproduct-tree)), at the cost of
$O(M(n) \log n)$ operations plus one initial division.

::: details
Generated from the same template as `_fq_poly_evaluate_fq_vec_fast`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_evaluate_fq_nmod_vec_fast_precomp(fq_nmod_struct * vs, const fq_nmod_struct * poly, slong plen, fq_nmod_poly_struct * const * tree, slong len, const fq_nmod_ctx_t ctx)

Evaluates `(poly, plen)` at the `len` points whose subproduct tree
`tree` was built by `_fq_nmod_poly_tree_build`, writing the values to `vs`.

::: details
Generated from the same template as `_fq_poly_evaluate_fq_vec_fast_precomp`, whose documentation has the details.
:::
:::

::: function
fq_nmod_poly_struct ** _fq_nmod_poly_tree_alloc(slong len, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_tree_free(fq_nmod_poly_struct ** tree, slong len, const fq_nmod_ctx_t ctx)

Allocate, respectively free, a subproduct tree for `len` points, with
$\lceil \log_2 \mathit{len} \rceil + 1$ levels of initialised polynomials.

::: details
Generated from the same template as `_fq_poly_tree_alloc`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_tree_build(fq_nmod_poly_struct ** tree, const fq_nmod_struct * roots, slong len, const fq_nmod_ctx_t ctx)

Fills the subproduct tree for the `len` points `roots`: level 0 holds the
linear factors $x - a_i$, and each higher level the products of pairs of
the level below.

::: details
Generated from the same template as `_fq_poly_tree_build`, whose documentation has the details.
:::
:::

## Composition {#composition}

::: summary algorithms/polynomials/composition :::

::: function
void _fq_nmod_poly_compose(fq_nmod_struct * rop, const fq_nmod_struct * op1, slong len1, const fq_nmod_struct * op2, slong len2, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_compose(fq_nmod_poly_t rop, const fq_nmod_poly_t op1, const fq_nmod_poly_t op2, const fq_nmod_ctx_t ctx)

Set `rop` to the composition $g(h(x))$ of $g$ = `op1` and $h$ = `op2`.
The underscore version requires nonzero inputs, space for
$(\mathit{len1} - 1)(\mathit{len2} - 1) + 1$ coefficients and no
aliasing.

::: details
Generated from the same template as `_fq_poly_compose`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_compose_mod(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_compose_mod(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t poly3, const fq_nmod_ctx_t ctx)

Set `res` to $f(g) \bmod h$ for $f$ = `poly1`, $g$ = `poly2` and
$h$ = `poly3` $\ne 0$, choosing between Horner's rule and the
Brent–Kung algorithm. The underscore version requires $g$ reduced modulo
$h$ and zero-padded to length `lenh - 1`, writes `lenh - 1`
coefficients and does not allow aliasing.

::: details
Generated from the same template as `_fq_poly_compose_mod`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_compose_mod_preinv(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_struct * hinv, slong lenhinv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_compose_mod_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t poly3, const fq_nmod_poly_t poly3inv, const fq_nmod_ctx_t ctx)

As `fq_nmod_poly_compose_mod`, given the inverse of the reversal of $h$.

::: details
Generated from the same template as `_fq_poly_compose_mod_preinv`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_compose_mod_horner(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_compose_mod_horner(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t poly3, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_compose_mod_horner_preinv(fq_nmod_struct * res, const fq_nmod_struct * f, slong lenf, const fq_nmod_struct * g, const fq_nmod_struct * h, slong lenh, const fq_nmod_struct * hinv, slong lenhinv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_compose_mod_horner_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t poly3, const fq_nmod_poly_t poly3inv, const fq_nmod_ctx_t ctx)

Modular composition by Horner's rule, with $\deg f$ multiplications
modulo $h$. Requirements as for `fq_nmod_poly_compose_mod`; the `_preinv`
versions also require $\mathit{lenf} < \mathit{lenh}$ and the inverse of
the reversal of $h$.

::: details
Generated from the same template as `_fq_poly_compose_mod_horner`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_compose_mod_brent_kung(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, const fq_nmod_struct * poly3, slong len3, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_compose_mod_brent_kung(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t poly3, const fq_nmod_ctx_t ctx)
void _fq_nmod_poly_compose_mod_brent_kung_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_struct * poly2, const fq_nmod_struct * poly3, slong len3, const fq_nmod_struct * poly3inv, slong len3inv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_compose_mod_brent_kung_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t poly3, const fq_nmod_poly_t poly3inv, const fq_nmod_ctx_t ctx)

Modular composition $f(g) \bmod h$ by the Brent–Kung algorithm: the
powers $g^i \bmod h$ for $i < m = \lceil \sqrt{\mathit{lenf}} \rceil$
are computed, $f$ is split into $m$ blocks, the block compositions are
obtained by one $m \times m$ by $m \times \deg h$ matrix product, and
the results are combined by Horner's rule in $g^m$. Requires
$\deg f < \deg h$ and $g$ reduced modulo $h$; the `_preinv` versions also
require the inverse of the reversal of $h$. No aliasing in the
underscore versions.

::: details
Generated from the same template as `_fq_poly_compose_mod_brent_kung`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_precompute_matrix(fq_nmod_mat_t A, const fq_nmod_struct * poly1, const fq_nmod_struct * poly2, slong len2, const fq_nmod_struct * poly2inv, slong len2inv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_precompute_matrix(fq_nmod_mat_t A, const fq_nmod_poly_t poly1, const fq_nmod_poly_t poly2, const fq_nmod_poly_t poly2inv, const fq_nmod_ctx_t ctx)

Set the rows of the matrix $A$ to the coefficients of the powers
$g^i \bmod h$, $0 \le i < m$, for $g$ = `poly1` and $h$ = `poly2`, as
needed by `fq_nmod_poly_compose_mod_brent_kung_precomp_preinv`. $A$ must have
exactly $m = \lfloor \sqrt{\deg h} \rfloor + 1$ rows and $\deg h$
columns (otherwise `fq_nmod_poly_precompute_matrix` throws `FLINT_ERROR`), and
`poly2inv` must be the inverse of the reversal of $h$. In the underscore
version $g$ must be reduced and zero-padded to length `len2 - 1`.

::: details
Generated from the same template as `_fq_poly_precompute_matrix`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_compose_mod_brent_kung_precomp_preinv(fq_nmod_struct * res, const fq_nmod_struct * poly1, slong len1, const fq_nmod_mat_t A, const fq_nmod_struct * poly3, slong len3, const fq_nmod_struct * poly3inv, slong len3inv, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_compose_mod_brent_kung_precomp_preinv(fq_nmod_poly_t res, const fq_nmod_poly_t poly1, const fq_nmod_mat_t A, const fq_nmod_poly_t poly3, const fq_nmod_poly_t poly3inv, const fq_nmod_ctx_t ctx)

Set `res` to $f(g) \bmod h$ for $f$ = `poly1` with $\deg f < \deg h$,
using the matrix $A$ of powers of $g$ computed by
`fq_nmod_poly_precompute_matrix`. This amortises the precomputation over many
compositions with the same $g$, as in the iterated Frobenius
computations of factorisation.

::: details
Generated from the same template as `_fq_poly_compose_mod_brent_kung_precomp_preinv`, whose documentation has the details.
:::
:::

::: function
void _fq_nmod_poly_reduce_matrix_mod_poly(fq_nmod_mat_t A, const fq_nmod_mat_t B, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)

Sets the rows of $A$ to the rows of $B$, considered as polynomials,
reduced modulo $f$; $A$ must have $\deg f$ columns. Used for the
precomputation of the Brent–Kung matrix.

::: details
Generated from the same template as `_fq_poly_reduce_matrix_mod_poly`, whose documentation has the details.
:::
:::

## Characteristic and minimal polynomial of a matrix {#charpoly}

::: function
void fq_nmod_mat_charpoly(fq_nmod_poly_t p, const fq_nmod_mat_t M, const fq_nmod_ctx_t ctx)

Sets $p$ to the characteristic polynomial $\det(xI - M)$ of the square
matrix $M$; throws `FLINT_ERROR` if $M$ is not square. This is a wrapper
of `gr_mat_charpoly`, which uses the Berkowitz algorithm for
$n \le 12$ and Danilevsky's algorithm otherwise (for `fq_zech`, from
$n = 5$); see [](algorithms/linear-algebra/charpoly).

::: details
Generated from the same template as `fq_mat_charpoly`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_mat_minpoly(fq_nmod_poly_t p, const fq_nmod_mat_t X, const fq_nmod_ctx_t ctx)

Sets $p$ to the minimal polynomial of the square matrix $X$, the monic
polynomial of least degree with $p(X) = 0$; throws `FLINT_ERROR` if $X$
is not square. This is a wrapper of `gr_mat_minpoly_field`.

::: details
Generated from the same template as `fq_mat_minpoly`, whose documentation has the details.
:::
:::

## Input and output {#io}

::: function
int _fq_nmod_poly_fprint(FILE * file, const fq_nmod_struct * poly, slong len, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_fprint(FILE * file, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
int _fq_nmod_poly_print(const fq_nmod_struct * poly, slong len, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_print(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)

Print the polynomial to `file`, respectively `stdout`: the length, two
spaces, and the coefficients in the format of `fq_nmod_fprint` separated by
single spaces. Return a positive value on success and a non-positive value
on failure.

::: details
Generated from the same template as `_fq_poly_fprint`, whose documentation has the details.
:::
:::

::: function
int _fq_nmod_poly_fprint_pretty(FILE * file, const fq_nmod_struct * poly, slong len, const char * x, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_fprint_pretty(FILE * file, const fq_nmod_poly_t poly, const char * x, const fq_nmod_ctx_t ctx)
int _fq_nmod_poly_print_pretty(const fq_nmod_struct * poly, slong len, const char * x, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_print_pretty(const fq_nmod_poly_t poly, const char * x, const fq_nmod_ctx_t ctx)

Print the polynomial in the variable `x`, with coefficients other than 1
in parentheses in the format of `fq_nmod_fprint_pretty`, for example
`x^2+(a+1)*x+(a)`. Return a positive value on success.

::: details
Generated from the same template as `_fq_poly_fprint_pretty`, whose documentation has the details.
:::
:::

::: function
char * _fq_nmod_poly_get_str(const fq_nmod_struct * poly, slong len, const fq_nmod_ctx_t ctx)
char * fq_nmod_poly_get_str(const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)
char * _fq_nmod_poly_get_str_pretty(const fq_nmod_struct * poly, slong len, const char * x, const fq_nmod_ctx_t ctx)
char * fq_nmod_poly_get_str_pretty(const fq_nmod_poly_t poly, const char * x, const fq_nmod_ctx_t ctx)

Return the plain, respectively pretty, representation as a newly
allocated string, which the caller must free with `flint_free`.

::: details
Generated from the same template as `_fq_poly_get_str`, whose documentation has the details.
:::
:::
