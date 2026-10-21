---
title: padic_poly.h
header: padic_poly.h
description: Polynomials over the $p$-adic numbers, stored as $p^v f(x)$ with $f \in \Z[x]$ reduced modulo $p^{N - v}$.
category: p-adic numbers
order: 2
---

This module implements dense polynomials over $\mathbf{Q}_p$ in the
fixed-precision model of `padic.h`. A polynomial is stored as
$p^v f(x)$ with $v \in \Z$ and $f \in \Z[x]$, together with a precision
$N$ that applies to every coefficient: the polynomial represents the
class of $p^v f$ modulo $p^N \mathbf{Z}_p[x]$. As in `padic`, every
function computes its result exactly from the stored inputs and reduces
it to the precision of the output. The module provides arithmetic,
powering, power series inversion, derivatives, evaluation and
composition; multiplication is done in $\Z[x]$ by `fmpz_poly`. The type
is also the representation of elements of unramified extensions in
`qadic.h`. Polynomials over $\mathbf{Q}_p$ with error tracking, division
and GCDs are available through the generic rings, as `gr_poly` over
`padic_radix`.

::: summary algorithms/p-adic/overview :::

## Representation {#representation}

::: type
padic_poly_struct
padic_poly_t

A polynomial $p^v f(x)$ over $\mathbf{Q}_p$. The fields are `coeffs`
(an array of `alloc` coefficients of $f$), `alloc`, `length`, `val`
($v$) and `N` (the precision). A `padic_poly_t` is an array of one
`padic_poly_struct`.
:::

- The polynomial is *normalised* if its length is zero or its leading
  coefficient is nonzero.
- It is *canonical* if it is normalised and either some coefficient of
  $f$ is not divisible by $p$, or $f = 0$ and $v = 0$. Then $v$ is the
  minimum of the valuations of the coefficients.
- It is *reduced* if it is canonical and all coefficients of $f$ lie in
  $[0, p^{N - v})$; a reduced polynomial with $v \ge N$ is zero.

Functions require canonical inputs and produce reduced outputs unless
stated otherwise, and allow aliasing between inputs and outputs. A
single precision for all coefficients means that a coefficient of
valuation $w > v$ is stored to relative precision $N - w$ only
([](algorithms/p-adic/precision#fixed)).

The functions with a leading underscore take the polynomial as
`(fmpz * poly, slong val, slong len)` and the output as
`(fmpz * rop, slong * rval)` with the precision $N$ of the output and a
`padic_ctx_t` for $p$; they require sufficient space in the output
array and do not normalise it.

## Memory management {#memory}

::: function
void padic_poly_init(padic_poly_t poly)
void padic_poly_init2(padic_poly_t poly, slong alloc, slong prec)

Initialise `poly` to zero, with precision `PADIC_DEFAULT_PREC` (that
is, 20) and no allocated coefficients, respectively with precision
`prec` and space for `alloc` coefficients.
:::

::: function
void padic_poly_clear(padic_poly_t poly)

Clears `poly`.
:::

::: function
void padic_poly_realloc(padic_poly_t poly, slong alloc, const fmpz_t p)

Reallocates `poly` to have space for exactly `alloc` coefficients. If
the length is larger, the polynomial is first truncated and
canonicalised with respect to $p$. With `alloc` $= 0$ the polynomial is
cleared and reinitialised, which also resets its precision to
`PADIC_DEFAULT_PREC`.
:::

::: function
void padic_poly_fit_length(padic_poly_t poly, slong len)

Ensures that `poly` has space for at least `len` coefficients, at least
doubling the allocation when it grows. The value is not changed.
:::

::: function
void _padic_poly_set_length(padic_poly_t poly, slong len)

Sets the length to `len`, first zeroing (and freeing) the coefficients
beyond `len`. The result may not be normalised or canonical.
:::

::: function
void _padic_poly_normalise(padic_poly_t poly)

Decreases the length until the leading coefficient is nonzero.
:::

::: function
void _padic_poly_canonicalise(fmpz * poly, slong * v, slong len, const fmpz_t p)
void padic_poly_canonicalise(padic_poly_t poly, const fmpz_t p)

Bring a normalised polynomial into canonical form by dividing $f$ by
the largest power of $p$ dividing all its coefficients and adding the
exponent to $v$; set $v = 0$ for the zero polynomial. No reduction is
done.
:::

::: function
void padic_poly_reduce(padic_poly_t poly, const padic_ctx_t ctx)

Reduces the canonical `poly` modulo $p^N$: reduces the coefficients
modulo $p^{N - v}$ and normalises, or sets `poly` to zero if $v \ge N$.
:::

::: function
void padic_poly_truncate(padic_poly_t poly, slong n, const fmpz_t p)

Truncates `poly` to length at most $n$ and brings it into canonical
form. A reduced input gives a reduced output.
:::

## Parameters {#parameters}

::: function
slong padic_poly_degree(const padic_poly_t poly)
slong padic_poly_length(const padic_poly_t poly)

Return the degree ($-1$ for zero) and the length of `poly`.
:::

::: function
slong padic_poly_val(const padic_poly_t poly)

Returns the valuation $v$, which for a canonical polynomial is the
minimum valuation of its coefficients (0 for the zero polynomial). The
name is also defined as a macro, so `padic_poly_val(poly)` can be
assigned to.
:::

::: macro
padic_poly_prec(poly)

The precision $N$ of `poly`, as an lvalue. After lowering it, call
`padic_poly_reduce`.
:::

## Random generation {#random}

::: function
void padic_poly_randtest(padic_poly_t f, flint_rand_t state, slong len, const padic_ctx_t ctx)

Sets $f$ to a random reduced polynomial of length at most `len`, with
valuation drawn as for `padic_randtest`. With `len` $= 0$, $f$ is not
changed.
:::

::: function
void padic_poly_randtest_not_zero(padic_poly_t f, flint_rand_t state, slong len, const padic_ctx_t ctx)

Sets $f$ to a random nonzero reduced polynomial of length at most
`len`. Throws if `len` is zero.
:::

::: function
void padic_poly_randtest_val(padic_poly_t f, flint_rand_t state, slong val, slong len, const padic_ctx_t ctx)

Sets $f$ to a random reduced polynomial of length at most `len` with
valuation `val`, or to zero if `val` $\ge N$. With `len` $= 0$, $f$ is
not changed.
:::

## Assignment and conversion {#conversion}

The functions that set a `padic_poly_t` reduce the value to the
precision of the output.

::: function
void padic_poly_set(padic_poly_t poly1, const padic_poly_t poly2, const padic_ctx_t ctx)

Sets `poly1` to `poly2`.
:::

::: function
void padic_poly_set_padic(padic_poly_t poly, const padic_t x, const padic_ctx_t ctx)

Sets `poly` to the constant polynomial $x$.
:::

::: function
void padic_poly_set_si(padic_poly_t poly, slong x, const padic_ctx_t ctx)
void padic_poly_set_ui(padic_poly_t poly, ulong x, const padic_ctx_t ctx)
void padic_poly_set_fmpz(padic_poly_t poly, const fmpz_t x, const padic_ctx_t ctx)
void padic_poly_set_fmpq(padic_poly_t poly, const fmpq_t x, const padic_ctx_t ctx)

Set `poly` to the constant polynomial $x$.
:::

::: function
void padic_poly_set_fmpz_poly(padic_poly_t rop, const fmpz_poly_t op, const padic_ctx_t ctx)
void padic_poly_set_fmpq_poly(padic_poly_t rop, const fmpq_poly_t op, const padic_ctx_t ctx)

Set `rop` to the integer or rational polynomial `op`. For a rational
polynomial the power of $p$ in the common denominator becomes a negative
valuation and the rest of the denominator is inverted modulo
$p^{N - v}$.
:::

::: function
int padic_poly_get_fmpz_poly(fmpz_poly_t rop, const padic_poly_t op, const padic_ctx_t ctx)

If $v \ge 0$, sets `rop` to the integer polynomial $p^v f$ and returns
1. Otherwise returns 0 and leaves `rop` unchanged.
:::

::: function
void padic_poly_get_fmpq_poly(fmpq_poly_t rop, const padic_poly_t op, const padic_ctx_t ctx)

Sets `rop` to the rational polynomial $p^v f$.
:::

::: function
void padic_poly_zero(padic_poly_t poly)

Sets `poly` to zero.
:::

::: function
void padic_poly_one(padic_poly_t poly)

Sets `poly` to 1, or to zero if its precision is not positive.
:::

::: function
void padic_poly_swap(padic_poly_t poly1, padic_poly_t poly2)

Swaps `poly1` and `poly2`, including their precisions, in constant
time.
:::

## Coefficients {#coefficients}

::: function
void padic_poly_get_coeff_padic(padic_t c, const padic_poly_t poly, slong n, const padic_ctx_t ctx)

Sets $c$ to the coefficient of $x^n$, which is zero for $n$ beyond the
length, reduced to the precision of $c$ (which is not changed).
:::

::: function
void padic_poly_set_coeff_padic(padic_poly_t poly, slong n, const padic_t c, const padic_ctx_t ctx)

Sets the coefficient of $x^n$ to $c$, reduced to the precision of
`poly`. If $c$ has smaller valuation than `poly`, all coefficients are
rescaled, so the cost can be linear in the length. The result is
normalised: setting the leading coefficient to zero (or
to an element of valuation at least the precision of `poly`) reduces
the length.
:::

## Comparison {#comparison}

::: function
int padic_poly_equal(const padic_poly_t poly1, const padic_poly_t poly2)

Returns whether the stored values are equal, that is, whether the
valuations, lengths and coefficients agree. The precisions are not
compared.
:::

::: function
int padic_poly_is_zero(const padic_poly_t poly)
int padic_poly_is_one(const padic_poly_t poly)

Return whether `poly` is zero, respectively the constant 1.
:::

## Addition and subtraction {#addition}

::: function
void _padic_poly_add(fmpz * rop, slong * rval, slong N, const fmpz * op1, slong val1, slong len1, slong N1, const fmpz * op2, slong val2, slong len2, slong N2, const padic_ctx_t ctx)
void _padic_poly_sub(fmpz * rop, slong * rval, slong N, const fmpz * op1, slong val1, slong len1, slong N1, const fmpz * op2, slong val2, slong len2, slong N2, const padic_ctx_t ctx)

Set `(rop, *rval)` of length $\max(len_1, len_2)$ to the sum,
respectively difference, of the two polynomials, reduced modulo $p^N$
but not normalised. The input precisions $N_1, N_2$ are ignored. Allow
`rop` to be equal to `op1` or `op2`, provided it has space for
$\max(len_1, len_2)$ coefficients.
:::

::: function
void padic_poly_add(padic_poly_t f, const padic_poly_t g, const padic_poly_t h, const padic_ctx_t ctx)
void padic_poly_sub(padic_poly_t f, const padic_poly_t g, const padic_poly_t h, const padic_ctx_t ctx)

Set $f$ to $g + h$, respectively $g - h$.
:::

::: function
void padic_poly_neg(padic_poly_t f, const padic_poly_t g, const padic_ctx_t ctx)

Sets $f$ to $-g$.
:::

## Scalar multiplication {#scalar}

::: function
void _padic_poly_scalar_mul_padic(fmpz * rop, slong * rval, slong N, const fmpz * op, slong val, slong len, const padic_t c, const padic_ctx_t ctx)

Sets `(rop, *rval, len)` to `(op, val, len)` times $c$, reduced modulo
$p^N$. If $c = 0$ or the valuation of the product is at least $N$, the
array is zeroed and `*rval` set to 0. Requires a nonzero polynomial for
the result to be canonical.
:::

::: function
void padic_poly_scalar_mul_padic(padic_poly_t rop, const padic_poly_t op, const padic_t c, const padic_ctx_t ctx)

Sets `rop` to $c \cdot op$.
:::

## Multiplication and powering {#multiplication}

::: function
void _padic_poly_mul(fmpz * rop, slong * rval, slong N, const fmpz * op1, slong val1, slong len1, const fmpz * op2, slong val2, slong len2, const padic_ctx_t ctx)

Sets `(rop, *rval, len1 + len2 - 1)` to the product, reduced modulo
$p^N$. Requires $len_1 \ge len_2 \ge 1$ and $val_1 + val_2 < N$. Does
not allow aliasing.
:::

::: function
void padic_poly_mul(padic_poly_t res, const padic_poly_t poly1, const padic_poly_t poly2, const padic_ctx_t ctx)

Sets `res` to the product of `poly1` and `poly2`.

::: details
The integer parts are multiplied in $\Z[x]$ with `_fmpz_poly_mul`
([](algorithms/polynomials/multiplication-integers)) and reduced
modulo $p^{N - v_1 - v_2}$. No canonicalisation is needed: by Gauss's
lemma the product of two polynomials that are nonzero modulo $p$ is
nonzero modulo $p$.
:::
:::

::: function
void _padic_poly_pow(fmpz * rop, slong * rval, slong N, const fmpz * op, slong val, slong len, ulong e, const padic_ctx_t ctx)

Sets `(rop, *rval, e(len - 1) + 1)` to `(op, val, len)` raised to the
power $e$, reduced modulo $p^N$. Requires $e \ge 2$, $len \ge 1$ and
$e \cdot val < N$. Does not allow aliasing.
:::

::: function
void padic_poly_pow(padic_poly_t rop, const padic_poly_t op, ulong e, const padic_ctx_t ctx)

Sets `rop` to $op^e$, with $op^0 = 1$. The power is computed in
$(\Z/p^{N - ev}\Z)[x]$ by `_fmpz_mod_poly_pow`.

::: details
If `op` has valuation $v < 0$ and is only known modulo $p^N$, its power
is only known modulo $p^{N + (e - 1)v}$
([](algorithms/p-adic/precision#prop-ring)); the function does not
account for this.
:::
:::

## Power series inversion {#inv-series}

::: function
void padic_poly_inv_series(padic_poly_t g, const padic_poly_t f, slong n, const padic_ctx_t ctx)

Sets $g$ to the power series inverse of $f$ modulo $x^n$, $n \ge 1$,
reduced to the precision of $g$. The constant coefficient of $f$ must
have the minimal valuation of all coefficients, that is, the constant
coefficient of the stored integer polynomial must be a $p$-adic unit;
otherwise, and if it is zero, the function throws. If $-v \ge N$ the
result is zero.

::: details
The result has valuation $-v$. Its integer part is computed modulo
$p^{N + v}$ by Newton iteration ([](algorithms/polynomials/newton)),
starting from the inverse of the constant coefficient computed by
`_padic_inv` and using middle products over $\Z$ followed by reduction.
:::
:::

## Derivative {#derivative}

::: function
void _padic_poly_derivative(fmpz * rop, slong * rval, slong N, const fmpz * op, slong val, slong len, const padic_ctx_t ctx)

Sets `(rop, *rval, len - 1)` to the derivative of `(op, val, len)`,
reduced modulo $p^N$ and canonicalised. Allows aliasing.
:::

::: function
void padic_poly_derivative(padic_poly_t rop, const padic_poly_t op, const padic_ctx_t ctx)

Sets `rop` to the derivative of `op`.
:::

## Shifting {#shifting}

::: function
void padic_poly_shift_left(padic_poly_t rop, const padic_poly_t op, slong n, const padic_ctx_t ctx)

Sets `rop` to $x^n \cdot op$, $n \ge 0$. Throws if the precision of
`rop` is smaller than that of `op`. No reduction is done for $n > 0$.
:::

::: function
void padic_poly_shift_right(padic_poly_t rop, const padic_poly_t op, slong n, const padic_ctx_t ctx)

Sets `rop` to the quotient of `op` by $x^n$, $n \ge 0$, discarding the
coefficients of $x^0, \ldots, x^{n-1}$, and canonicalises the result.
For $n > 0$ no reduction is done, so the result is reduced only if the
precision of `rop` is at least that of `op`.
:::

## Evaluation and composition {#evaluation}

::: function
void _padic_poly_evaluate_padic(fmpz_t u, slong * v, slong N, const fmpz * poly, slong val, slong len, const fmpz_t a, slong b, const padic_ctx_t ctx)
void padic_poly_evaluate_padic(padic_t y, const padic_poly_t poly, const padic_t a, const padic_ctx_t ctx)

Set $y = p^v u$ to the value of `poly` at $a = p^b a'$, reduced modulo
$p^N$ (the precision of $y$). `padic_poly_evaluate_padic` allows $y$ to
be the same object as $a$.

::: details
For $b \ge 0$ the integer polynomial is evaluated at $p^b a'$ by Horner's
rule modulo $p^{N - val}$. For $b < 0$ the polynomial
$p^{-nb} f(x)$ with $n = len - 1$ is integral at $x = a'$ after the
substitution, so $p^{-nb} f(p^b a')$ is evaluated modulo
$p^{N - val - nb}$ and the valuation adjusted by $nb$. If `poly` and
$a$ are only known to precision $N$ and $b < 0$, the value is only
known to precision $N + nb$; this is not taken into account.
:::
:::

::: function
void _padic_poly_compose(fmpz * rop, slong * rval, slong N, const fmpz * op1, slong val1, slong len1, const fmpz * op2, slong val2, slong len2, const padic_ctx_t ctx)

Sets `(rop, *rval, (len1 - 1)(len2 - 1) + 1)` to the composition
$f(g(x))$ of $f$ = `op1` and $g$ = `op2`, reduced modulo $p^N$ and
canonicalised. Requires $len_1 \ge 1$. Does not allow aliasing.
:::

::: function
void padic_poly_compose(padic_poly_t rop, const padic_poly_t op1, const padic_poly_t op2, const padic_ctx_t ctx)

Sets `rop` to $op_1(op_2(x))$. The composition is computed in
$(\Z/p^M\Z)[x]$ by `_fmpz_mod_poly_compose`, after scaling the
coefficients so that they are integral, as for evaluation.
:::

::: function
void _padic_poly_compose_pow(fmpz * rop, slong * rval, slong N, const fmpz * op, slong val, slong len, slong k, const padic_ctx_t ctx)

Sets `(rop, *rval, (len - 1)k + 1)` to $f(x^k)$ for $f$ = `(op, val,
len)`, $k \ge 1$, $len \ge 1$. Only the length-one case is reduced
modulo $p^N$. Allows aliasing.
:::

::: function
void padic_poly_compose_pow(padic_poly_t rop, const padic_poly_t op, slong k, const padic_ctx_t ctx)

Sets `rop` to $op(x^k)$, $k \ge 1$. The coefficients are moved, not
reduced, so the result is reduced only if the precision of `rop` is at
least that of `op`.
:::

## Input and output {#io}

::: function
int _padic_poly_fprint(FILE * file, const fmpz * poly, slong val, slong len, const padic_ctx_t ctx)
int padic_poly_fprint(FILE * file, const padic_poly_t poly, const padic_ctx_t ctx)
int _padic_poly_print(const fmpz * poly, slong val, slong len, const padic_ctx_t ctx)
int padic_poly_print(const padic_poly_t poly, const padic_ctx_t ctx)

Print the polynomial as its length, two spaces and the coefficients
separated by spaces, each in the printing mode of the context (all three
modes are supported); the zero polynomial prints as `0`. Return 1.
:::

::: function
int _padic_poly_fprint_pretty(FILE * file, const fmpz * poly, slong len, slong val, const char * var, const padic_ctx_t ctx)
int padic_poly_fprint_pretty(FILE * file, const padic_poly_t poly, const char * var, const padic_ctx_t ctx)
int _padic_poly_print_pretty(const fmpz * poly, slong len, slong val, const char * var, const padic_ctx_t ctx)
int padic_poly_print_pretty(const padic_poly_t poly, const char * var, const padic_ctx_t ctx)

Print the polynomial as a sum of terms in the variable `var`, highest
degree first, each nonconstant coefficient in parentheses. The
underscore versions take `len, val` in the opposite order to
`_padic_poly_fprint`. Return 1.
:::

::: function
int padic_poly_debug(const padic_poly_t poly)

Prints the raw fields `alloc`, `length`, `val`, `N` and all allocated
coefficients. Returns 1.
:::

## Testing {#testing}

::: function
int _padic_poly_is_canonical(const fmpz * op, slong val, slong len, const padic_ctx_t ctx)
int padic_poly_is_canonical(const padic_poly_t op, const padic_ctx_t ctx)

Return whether the polynomial is in canonical form.
:::

::: function
int _padic_poly_is_reduced(const fmpz * op, slong val, slong len, slong N, const padic_ctx_t ctx)
int padic_poly_is_reduced(const padic_poly_t op, const padic_ctx_t ctx)

Return whether the polynomial is reduced modulo $p^N$: canonical, with
$v < N$ (or zero) and all coefficients in $[0, p^{N - v})$.
:::

## Internal functions {#internal}

::: function
slong _fmpz_vec_ord_p(const fmpz * vec, slong len, const fmpz_t p)

Returns the minimum $p$-adic valuation of the nonzero entries of
`(vec, len)`, or 0 if all entries are zero or `len` is zero.
:::
