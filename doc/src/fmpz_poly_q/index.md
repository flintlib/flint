---
title: fmpz_poly_q.h
header: fmpz_poly_q.h
description: Rational functions in $\mathbb{Q}(x)$, stored as canonical quotients of integer polynomials.
category: Univariate polynomials
order: 4
---

This module implements arithmetic in the field $\mathbb{Q}(x)$ of rational
functions. An element is an `fmpz_poly_q_t`, stored as a quotient of two
integer polynomials of type `fmpz_poly_t` in a canonical form: numerator
and denominator are coprime in $\Z[x]$ and the denominator has positive
leading coefficient. Every function except `fmpz_poly_q_canonicalise`
assumes that its inputs are canonical and returns a canonical result.
Sums, products and quotients are computed with Henrici's algorithms
[@Hen1956], which keep the intermediate polynomials small by taking gcds
of the factors before multiplying, and derivatives with a variant of the
algorithm of Horowitz [@Hor1972]; the guide
[](fmpz_poly_q/representation) explains these algorithms and shows how
to work with the numerator and denominator directly.

For multivariate rational functions see `fmpz_mpoly_q`; for polynomials
with rational coefficients, `fmpq_poly`.

## Types and canonical form {#types}

::: type
fmpz_poly_q_struct
fmpz_poly_q_t

A rational function over $\mathbb{Q}$. The structure has two fields `num` and
`den`, pointers to separately allocated `fmpz_poly_struct`s holding the
numerator and the denominator. An `fmpz_poly_q_t` is an array of one
`fmpz_poly_q_struct`, so it is passed by reference.
:::

- A rational function is *canonical* if its denominator is nonzero with
  positive leading coefficient and the numerator and denominator are
  coprime in $\Z[x]$, which includes that their contents are coprime. The
  zero function is $0/1$. This representation is unique.
- All functions allow aliasing between their inputs and outputs. Two
  distinct `fmpz_poly_q_t` must not share their numerator or denominator
  objects.

::: macro
fmpz_poly_q_numref(op)
fmpz_poly_q_denref(op)

Pointers (`fmpz_poly_struct *`) to the numerator and the denominator of
`op`, which can be passed to any `fmpz_poly` function. If they are
modified, the rational function must be brought back into canonical form
with `fmpz_poly_q_canonicalise`.
:::

::: function
void fmpz_poly_q_canonicalise(fmpz_poly_q_t rop)

Brings `rop` into canonical form by dividing the numerator and the
denominator by their gcd in $\Z[x]$ and making the leading coefficient of
the denominator positive. The denominator must be nonzero; a zero
denominator causes an exception (`FLINT_ERROR`).
:::

::: function
int fmpz_poly_q_is_canonical(const fmpz_poly_q_t op)

Returns whether `op` is in canonical form.
:::

## Memory management {#memory}

::: function
void fmpz_poly_q_init(fmpz_poly_q_t rop)

Initialises `rop` as the zero function $0/1$, allocating its numerator and
denominator.
:::

::: function
void fmpz_poly_q_clear(fmpz_poly_q_t rop)

Frees the memory used by `rop`.
:::

## Assignment and basic operations {#assignment}

::: function
void fmpz_poly_q_set(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

Sets `rop` to `op`.
:::

::: function
void fmpz_poly_q_set_si(fmpz_poly_q_t rop, slong op)

Sets `rop` to the constant function `op`.
:::

::: function
void fmpz_poly_q_swap(fmpz_poly_q_t op1, fmpz_poly_q_t op2)

Swaps `op1` and `op2` efficiently, by exchanging the pointers to the
numerators and the denominators.
:::

::: function
void fmpz_poly_q_zero(fmpz_poly_q_t rop)
void fmpz_poly_q_one(fmpz_poly_q_t rop)

Set `rop` to 0 or 1.
:::

::: function
void fmpz_poly_q_neg(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

Sets `rop` to $-\mathit{op}$.
:::

::: function
void fmpz_poly_q_inv(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

Sets `rop` to $1/\mathit{op}$, by exchanging numerator and denominator and
fixing the sign. Raises an exception (`FLINT_ERROR`) if `op` is zero.
:::

## Randomisation {#randomisation}

::: function
void fmpz_poly_q_randtest(fmpz_poly_q_t poly, flint_rand_t state, slong len1, flint_bitcnt_t bits1, slong len2, flint_bitcnt_t bits2)
void fmpz_poly_q_randtest_not_zero(fmpz_poly_q_t poly, flint_rand_t state, slong len1, flint_bitcnt_t bits1, slong len2, flint_bitcnt_t bits2)

Set `poly` to a random rational function, the canonical form of $a/b$
with $a$ a random polynomial of length at most `len1` and coefficients
of at most `bits1` bits (`fmpz_poly_randtest`) and $b$ a random nonzero
polynomial of length at most `len2` and at most `bits2` bits. In the
`_not_zero` version $a$ is nonzero as well. Lengths and bit sizes less
than 1 are replaced by 1 where a nonzero polynomial is needed.
:::

## Comparison {#comparison}

::: function
int fmpz_poly_q_is_zero(const fmpz_poly_q_t op)
int fmpz_poly_q_is_one(const fmpz_poly_q_t op)

Return whether `op` is 0, respectively 1.
:::

::: function
int fmpz_poly_q_equal(const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

Returns whether `op1` and `op2` are equal. Since the canonical form is
unique, this compares numerators and denominators.
:::

## Addition and subtraction {#addition}

::: function
void fmpz_poly_q_add(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
void fmpz_poly_q_sub(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

Set `rop` to $\mathit{op1} + \mathit{op2}$, respectively $\mathit{op1} -
\mathit{op2}$, by Henrici's algorithm: with $g$ the gcd of the
denominators, only $\gcd(\text{numerator}, g)$ has to be cancelled.
:::

::: function
void fmpz_poly_q_add_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op)
void fmpz_poly_q_sub_in_place(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

Set `rop` to $\mathit{rop} + \mathit{op}$, respectively $\mathit{rop} -
\mathit{op}$. `fmpz_poly_q_add` and `fmpz_poly_q_sub` call these when the
output is aliased with the first input.
:::

::: function
void fmpz_poly_q_addmul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)
void fmpz_poly_q_submul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

Set `rop` to $\mathit{rop} + \mathit{op1} \cdot \mathit{op2}$,
respectively $\mathit{rop} - \mathit{op1} \cdot \mathit{op2}$, using a
temporary for the product.
:::

## Scalar multiplication and division {#scalar}

::: function
void fmpz_poly_q_scalar_mul_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x)
void fmpz_poly_q_scalar_mul_fmpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpz_t x)
void fmpz_poly_q_scalar_mul_fmpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpq_t x)

Set `rop` to $x \cdot \mathit{op}$.
:::

::: function
void fmpz_poly_q_scalar_div_si(fmpz_poly_q_t rop, const fmpz_poly_q_t op, slong x)
void fmpz_poly_q_scalar_div_fmpz(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpz_t x)
void fmpz_poly_q_scalar_div_fmpq(fmpz_poly_q_t rop, const fmpz_poly_q_t op, const fmpq_t x)

Set `rop` to $\mathit{op} / x$. A zero $x$ raises an exception
(`FLINT_DIVZERO`).
:::

## Multiplication and division {#multiplication}

::: function
void fmpz_poly_q_mul(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

Sets `rop` to $\mathit{op1} \cdot \mathit{op2}$. Writing $\mathit{op1} =
a/b$ and $\mathit{op2} = c/d$, Henrici's algorithm cancels $\gcd(a, d)$
and $\gcd(c, b)$ before multiplying, which gives the canonical result
directly.
:::

::: function
void fmpz_poly_q_div(fmpz_poly_q_t rop, const fmpz_poly_q_t op1, const fmpz_poly_q_t op2)

Sets `rop` to $\mathit{op1} / \mathit{op2}$, computed like the product
with $\mathit{op2}$ inverted. A zero `op2` raises an exception
(`FLINT_DIVZERO`).
:::

## Powering {#powering}

::: function
void fmpz_poly_q_pow(fmpz_poly_q_t rop, const fmpz_poly_q_t op, ulong exp)

Sets `rop` to $\mathit{op}^{\mathit{exp}}$, by raising the numerator and
the denominator to the power `exp`. For `exp` = 0 the result is 1, also
when `op` is zero.
:::

## Derivative {#derivative}

::: function
void fmpz_poly_q_derivative(fmpz_poly_q_t rop, const fmpz_poly_q_t op)

Sets `rop` to the derivative of `op`.

::: details
For $\mathit{op} = a/b$ with $b$ nonconstant, let $d = \gcd(b, b')$. Then
$$
\left(\frac{a}{b}\right)' = \frac{a' (b/d) - a (b'/d)}{b \cdot (b/d)},
$$
and the numerator and denominator of this expression can only have a
constant common factor, which is removed by dividing by the gcd of their
contents [@Hor1972].
:::
:::

## Evaluation {#evaluation}

::: function
int fmpz_poly_q_evaluate_fmpq(fmpq_t rop, const fmpz_poly_q_t f, const fmpq_t a)

If the denominator of `f` does not vanish at the rational number $a$,
sets `rop` to $f(a)$ and returns 0. Otherwise returns a nonzero value and
leaves `rop` unchanged.
:::

## Input and output {#io}

The string format is based on that of `fmpz_poly_get_str`: a polynomial
$a_0 + a_1 x + \cdots + a_n x^n$ is written as `"n+1  a_0 a_1 ... a_n"`,
with two spaces after the length, and the zero polynomial as `"0"`. A
rational function with numerator and denominator strings `num` and `den`
is written as `num/den`, or as `num` alone if the denominator is 1. The
pretty format uses `fmpz_poly_get_str_pretty` for both parts and puts
parentheses around nonconstant ones, for example `(3*t+1)/2`. These
functions are intended for debugging and small inputs, not as fast
parsers.

::: function
int fmpz_poly_q_set_str(fmpz_poly_q_t rop, const char * s)

Sets `rop` to the rational function described by `s`, in the format
above, and canonicalises it. Returns 0 on success and a nonzero value if
`s` is not valid, in which case the value of `rop` is unspecified. A
denominator string that describes the zero polynomial raises an
exception.
:::

::: function
char * fmpz_poly_q_get_str(const fmpz_poly_q_t op)
char * fmpz_poly_q_get_str_pretty(const fmpz_poly_q_t op, const char * x)

Return the string representation of `op`, in the plain format or in the
pretty format with the variable name `x`. The string is allocated with
`flint_malloc` and must be freed with `flint_free`.
:::

::: function
int fmpz_poly_q_print(const fmpz_poly_q_t op)
int fmpz_poly_q_print_pretty(const fmpz_poly_q_t op, const char * x)

Print the plain or the pretty representation of `op` to standard output,
and return 1.
:::
