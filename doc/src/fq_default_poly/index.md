---
title: fq_default_poly.h
header: fq_default_poly.h
description: Dense univariate polynomials over any finite field $\mathbb{F}_q$, in the representation chosen by `fq_default`.
category: Finite fields
order: 3
---

This module implements dense univariate polynomials over a finite field
given by an `fq_default_ctx_t` (see [](fq_default)). A polynomial is a
union of the polynomial types of the five representations
(`fq_zech_poly`, `fq_nmod_poly`, `fq_poly`, `nmod_poly` and
`fmpz_mod_poly`), and every function calls the function of the
representation of the context. The semantics are those of [](fq_poly),
where the algorithms are described: polynomials are normalised, outputs
may be aliased with inputs, and the fast algorithms (Kronecker
substitution multiplication, Newton division, half-GCD, Brent–Kung
modular composition) are used automatically. Over prime fields the
functions of `nmod_poly` and `fmpz_mod_poly` are used. Factorisation and
root finding are in `fq_default_poly_factor.h`.

::: summary algorithms/polynomials/multiplication :::

## Types {#types}

::: type
fq_default_poly_struct
fq_default_poly_t

A polynomial over an `fq_default` field: a union of an `fq_poly_t`, an
`fq_nmod_poly_t`, an `fq_zech_poly_t`, an `nmod_poly_t` and an
`fmpz_mod_poly_t`, of which the member matching the context is used.
`fq_default_poly_t` is an array of one union.
:::

## Memory management {#memory}

::: function
void fq_default_poly_init(fq_default_poly_t poly, const fq_default_ctx_t ctx)
void fq_default_poly_init2(fq_default_poly_t poly, slong alloc, const fq_default_ctx_t ctx)

Initialise `poly` as the zero polynomial; `fq_default_poly_init2`
allocates space for `alloc` coefficients.
:::

::: function
void fq_default_poly_realloc(fq_default_poly_t poly, slong alloc, const fq_default_ctx_t ctx)
void fq_default_poly_fit_length(fq_default_poly_t poly, slong len, const fq_default_ctx_t ctx)

Change the allocation to `alloc` coefficients (truncating if needed),
respectively ensure room for at least `len` coefficients.
:::

::: function
void fq_default_poly_clear(fq_default_poly_t poly, const fq_default_ctx_t ctx)

Frees the memory used by `poly`.
:::

::: function
void _fq_default_poly_set_length(fq_default_poly_t poly, slong len, const fq_default_ctx_t ctx)

Sets the length of `poly` to `len` without normalising, as
`_fq_poly_set_length`; for internal use.
:::

::: function
void fq_default_poly_truncate(fq_default_poly_t poly, slong len, const fq_default_ctx_t ctx)
void fq_default_poly_set_trunc(fq_default_poly_t poly1, fq_default_poly_t poly2, slong len, const fq_default_ctx_t ctx)

Truncate `poly` to length `len`, respectively set `poly1` to `poly2`
truncated to length `len`.
:::

## Length, degree and comparison {#properties}

::: function
slong fq_default_poly_length(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
slong fq_default_poly_degree(const fq_default_poly_t poly, const fq_default_ctx_t ctx)

Return the length, respectively the degree ($-1$ for zero).
:::

::: function
slong fq_default_poly_hamming_weight(const fq_default_poly_t op, const fq_default_ctx_t ctx)

Returns the number of nonzero coefficients.
:::

::: function
int fq_default_poly_equal(const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_ctx_t ctx)
int fq_default_poly_equal_trunc(const fq_default_poly_t poly1, const fq_default_poly_t poly2, slong n, const fq_default_ctx_t ctx)
int fq_default_poly_equal_fq_default(const fq_default_poly_t poly, const fq_default_t c, const fq_default_ctx_t ctx)

Return whether the polynomials are equal, equal modulo $x^n$, respectively
whether `poly` is the constant $c$.
:::

::: function
int fq_default_poly_is_zero(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
int fq_default_poly_is_one(const fq_default_poly_t op, const fq_default_ctx_t ctx)
int fq_default_poly_is_gen(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
int fq_default_poly_is_unit(const fq_default_poly_t op, const fq_default_ctx_t ctx)

Return whether the polynomial is 0, 1, $x$, respectively a nonzero
constant.
:::

## Assignment and conversion {#assignment}

::: function
void fq_default_poly_set(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)
void fq_default_poly_swap(fq_default_poly_t op1, fq_default_poly_t op2, const fq_default_ctx_t ctx)

Copy, respectively swap, polynomials.
:::

::: function
void fq_default_poly_zero(fq_default_poly_t poly, const fq_default_ctx_t ctx)
void fq_default_poly_one(fq_default_poly_t poly, const fq_default_ctx_t ctx)
void fq_default_poly_gen(fq_default_poly_t f, const fq_default_ctx_t ctx)
void fq_default_poly_set_fq_default(fq_default_poly_t poly, const fq_default_t c, const fq_default_ctx_t ctx)

Set the polynomial to 0, 1, $x$, respectively the constant $c$.
:::

::: function
void fq_default_poly_set_nmod_poly(fq_default_poly_t rop, const nmod_poly_t op, const fq_default_ctx_t ctx)
void fq_default_poly_set_fmpz_mod_poly(fq_default_poly_t rop, const fmpz_mod_poly_t op, const fq_default_ctx_t ctx)
void fq_default_poly_set_fmpz_poly(fq_default_poly_t rop, const fmpz_poly_t op, const fq_default_ctx_t ctx)

Set `rop` to the image of a polynomial over $\mathbb{F}_p$ (whose modulus
must be $p$), respectively over $\Z$ (with coefficients reduced modulo
$p$), in $\mathbb{F}_q[x]$.
:::

::: function
void fq_default_poly_make_monic(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)

Sets `rop` to `op` divided by its leading coefficient (zero stays zero).
:::

::: function
void fq_default_poly_reverse(fq_default_poly_t res, const fq_default_poly_t poly, slong n, const fq_default_ctx_t ctx)

Sets `res` to the reversal of `poly` as a polynomial of length $n$,
$x^{n-1} f(1/x)$ for $f$ = `poly` truncated to length $n$.
:::

::: function
ulong fq_default_poly_deflation(const fq_default_poly_t input, const fq_default_ctx_t ctx)
void fq_default_poly_deflate(fq_default_poly_t result, const fq_default_poly_t input, ulong deflation, const fq_default_ctx_t ctx)
void fq_default_poly_inflate(fq_default_poly_t result, const fq_default_poly_t input, ulong inflation, const fq_default_ctx_t ctx)

Return the largest $k$ such that `input` is a polynomial in $x^k$, write
$f(x^k)$ as $g(x)$, respectively compute $f(x^k)$, as for `fq_poly`.
:::

## Coefficients {#coefficients}

::: function
void fq_default_poly_get_coeff(fq_default_t x, const fq_default_poly_t poly, slong n, const fq_default_ctx_t ctx)

Sets $x$ to the coefficient of $x^n$ (zero beyond the length).
:::

::: function
void fq_default_poly_set_coeff(fq_default_poly_t poly, slong n, const fq_default_t x, const fq_default_ctx_t ctx)
void fq_default_poly_set_coeff_fmpz(fq_default_poly_t poly, slong n, const fmpz_t x, const fq_default_ctx_t ctx)

Set the coefficient of $x^n$ to the element $x$, respectively to the image
of the integer $x$ in $\mathbb{F}_p$.
:::

## Random polynomials {#random}

::: function
void fq_default_poly_randtest(fq_default_poly_t f, flint_rand_t state, slong len, const fq_default_ctx_t ctx)
void fq_default_poly_randtest_not_zero(fq_default_poly_t f, flint_rand_t state, slong len, const fq_default_ctx_t ctx)
void fq_default_poly_randtest_monic(fq_default_poly_t f, flint_rand_t state, slong len, const fq_default_ctx_t ctx)
void fq_default_poly_randtest_irreducible(fq_default_poly_t f, flint_rand_t state, slong len, const fq_default_ctx_t ctx)

Random polynomials of length at most `len` (nonzero for the second
function, which requires `len` $> 0$), random monic polynomials of length
`len`, and random monic irreducible polynomials of length `len`, for
testing.
:::

## Arithmetic {#arithmetic}

::: function
void fq_default_poly_add(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)
void fq_default_poly_sub(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)
void fq_default_poly_neg(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)

Set `rop` to the sum, the difference, respectively the negation.
:::

::: function
void fq_default_poly_add_si(fq_default_poly_t rop, const fq_default_poly_t op1, slong c, const fq_default_ctx_t ctx)

Sets `rop` to $\mathit{op1} + c$, with $c$ mapped to $\mathbb{F}_p$.
:::

::: function
void fq_default_poly_add_series(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, slong n, const fq_default_ctx_t ctx)
void fq_default_poly_sub_series(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, slong n, const fq_default_ctx_t ctx)

Set `rop` to the sum, respectively difference, truncated to length $n$.
:::

::: function
void fq_default_poly_scalar_mul_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
void fq_default_poly_scalar_div_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
void fq_default_poly_scalar_addmul_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)
void fq_default_poly_scalar_submul_fq_default(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_t x, const fq_default_ctx_t ctx)

Set `rop` to $x \cdot \mathit{op}$, $\mathit{op} / x$ (for $x \ne 0$),
respectively add $x \cdot \mathit{op}$ to or subtract it from `rop`.
:::

::: function
void fq_default_poly_mul(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)
void fq_default_poly_sqr(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)

Set `rop` to $\mathit{op1} \cdot \mathit{op2}$, respectively
$\mathit{op}^2$; see [](fq_poly/mul) for the algorithms over extension
fields.
:::

::: function
void fq_default_poly_mullow(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, slong n, const fq_default_ctx_t ctx)

Sets `rop` to the product truncated to length $n$.
:::

::: function
void fq_default_poly_mulhigh(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, slong start, const fq_default_ctx_t ctx)

Sets `rop` to a polynomial whose coefficients of index at least `start`
are those of the product $\mathit{op1} \cdot \mathit{op2}$; the lower
coefficients are zero for the extension field representations and
unspecified (zero or those of the product) for prime fields.
:::

::: function
void fq_default_poly_mulmod(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_poly_t f, const fq_default_ctx_t ctx)

Sets `res` to $\mathit{poly1} \cdot \mathit{poly2} \bmod f$ for $f \ne 0$.
:::

::: function
void fq_default_poly_pow(fq_default_poly_t rop, const fq_default_poly_t op, ulong e, const fq_default_ctx_t ctx)
void fq_default_poly_pow_trunc(fq_default_poly_t res, const fq_default_poly_t poly, ulong e, slong trunc, const fq_default_ctx_t ctx)

Set the output to $\mathit{op}^e$, respectively
$\mathit{poly}^e \bmod x^{\mathit{trunc}}$, with $\mathit{op}^0 = 1$.
:::

::: function
void fq_default_poly_powmod_ui_binexp(fq_default_poly_t res, const fq_default_poly_t poly, ulong e, const fq_default_poly_t f, const fq_default_ctx_t ctx)
void fq_default_poly_powmod_fmpz_binexp(fq_default_poly_t res, const fq_default_poly_t poly, const fmpz_t e, const fq_default_poly_t f, const fq_default_ctx_t ctx)

Set `res` to $\mathit{poly}^e \bmod f$ for $e \ge 0$ and $f \ne 0$, by
binary exponentiation.
:::

::: function
void fq_default_poly_shift_left(fq_default_poly_t rop, const fq_default_poly_t op, slong n, const fq_default_ctx_t ctx)
void fq_default_poly_shift_right(fq_default_poly_t rop, const fq_default_poly_t op, slong n, const fq_default_ctx_t ctx)

Set `rop` to $x^n \cdot \mathit{op}$, respectively
$\lfloor \mathit{op} / x^n \rfloor$.
:::

::: function
void fq_default_poly_derivative(fq_default_poly_t rop, const fq_default_poly_t op, const fq_default_ctx_t ctx)

Sets `rop` to the formal derivative of `op`.
:::

## Division and GCD {#division}

::: function
void fq_default_poly_divrem(fq_default_poly_t Q, fq_default_poly_t R, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)
void fq_default_poly_rem(fq_default_poly_t R, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)

Compute the quotient $Q$ and the remainder $R$ of the division of $A$ by
$B \ne 0$, with $A = BQ + R$ and $\deg R < \deg B$, respectively only the
remainder. $Q$ and $R$ must be distinct.
:::

::: function
int fq_default_poly_divides(fq_default_poly_t Q, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)

Returns 1 and sets $Q = A/B$ if $B$ divides $A$, and returns 0 otherwise.
$B$ must be nonzero.
:::

::: function
ulong fq_default_poly_remove(fq_default_poly_t f, const fq_default_poly_t g, const fq_default_ctx_t ctx)

Removes the highest power of $g$ (of degree at least 1) dividing $f$ and
returns its exponent.
:::

::: function
void fq_default_poly_gcd(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)

Sets `rop` to the monic GCD of the two polynomials (zero if both are
zero).
:::

::: function
void fq_default_poly_xgcd(fq_default_poly_t G, fq_default_poly_t S, fq_default_poly_t T, const fq_default_poly_t A, const fq_default_poly_t B, const fq_default_ctx_t ctx)

Computes the monic GCD $G$ of $A$ and $B$ and cofactors with
$SA + TB = G$, as `fq_poly_xgcd`. $G$, $S$ and $T$ must be distinct.
:::

## Power series {#series}

::: function
void fq_default_poly_inv_series(fq_default_poly_t Qinv, const fq_default_poly_t Q, slong n, const fq_default_ctx_t ctx)
void fq_default_poly_div_series(fq_default_poly_t Q, const fq_default_poly_t A, const fq_default_poly_t B, slong n, const fq_default_ctx_t ctx)

Set the output to $1/Q \bmod x^n$, respectively $A/B \bmod x^n$, for a
power series with nonzero constant term.
:::

::: function
void fq_default_poly_sqrt_series(fq_default_poly_t rop, const fq_default_poly_t op, slong n, fq_default_ctx_t ctx)
void fq_default_poly_invsqrt_series(fq_default_poly_t rop, const fq_default_poly_t op, slong n, fq_default_ctx_t ctx)

Set `rop` to $\sqrt{\mathit{op}}$, respectively $1/\sqrt{\mathit{op}}$,
modulo $x^n$, for a power series with constant term 1, in odd
characteristic.
:::

::: function
int fq_default_poly_sqrt(fq_default_poly_t rop, const fq_default_poly_t op, fq_default_ctx_t ctx)

If `op` is the square of a polynomial, sets `rop` to a square root and
returns 1; otherwise returns 0.
:::

## Evaluation and composition {#evaluation}

::: function
void fq_default_poly_evaluate_fq_default(fq_default_t res, const fq_default_poly_t f, const fq_default_t a, const fq_default_ctx_t ctx)

Sets `res` to $f(a)$.
:::

::: function
void fq_default_poly_compose(fq_default_poly_t rop, const fq_default_poly_t op1, const fq_default_poly_t op2, const fq_default_ctx_t ctx)

Sets `rop` to $g(h(x))$ for $g$ = `op1` and $h$ = `op2`.
:::

::: function
void fq_default_poly_compose_mod(fq_default_poly_t res, const fq_default_poly_t poly1, const fq_default_poly_t poly2, const fq_default_poly_t poly3, const fq_default_ctx_t ctx)

Sets `res` to $f(g) \bmod h$ for $f$ = `poly1`, $g$ = `poly2` and $h$ =
`poly3` $\ne 0$ ([](algorithms/polynomials/composition)).
:::

## Characteristic and minimal polynomial of a matrix {#charpoly}

::: function
void fq_default_mat_charpoly(fq_default_poly_t p, const fq_default_mat_t M, const fq_default_ctx_t ctx)
void fq_default_mat_minpoly(fq_default_poly_t p, const fq_default_mat_t X, const fq_default_ctx_t ctx)

Set $p$ to the characteristic polynomial $\det(xI - M)$, respectively the
minimal polynomial, of a square matrix over the field (see
`fq_mat_charpoly` and `fq_mat_minpoly`).
:::

## Input and output {#io}

::: function
int fq_default_poly_fprint(FILE * file, const fq_default_poly_t poly, const fq_default_ctx_t ctx)
int fq_default_poly_print(const fq_default_poly_t poly, const fq_default_ctx_t ctx)
char * fq_default_poly_get_str(const fq_default_poly_t poly, const fq_default_ctx_t ctx)

Print the polynomial, or return it as a newly allocated string, in the
plain format of the representation.
:::

::: function
int fq_default_poly_fprint_pretty(FILE * file, const fq_default_poly_t poly, const char * x, const fq_default_ctx_t ctx)
int fq_default_poly_print_pretty(const fq_default_poly_t poly, const char * x, const fq_default_ctx_t ctx)
char * fq_default_poly_get_str_pretty(const fq_default_poly_t poly, const char * x, const fq_default_ctx_t ctx)

Print the polynomial in the variable `x`, or return it as a newly
allocated string, with coefficients in parentheses in the readable
format of the field (for example `x^2+(a+1)*x+(a)`).
:::
