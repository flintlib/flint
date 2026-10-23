---
title: fmpz_mpoly_q.h
header: fmpz_mpoly_q.h
description: Multivariate rational functions over the rationals, as reduced fractions of integer polynomials, optionally with factored denominators.
category: Multivariate polynomials
order: 4
---

This module implements the field $\Q(x_0, \ldots, x_{n-1})$ of rational
functions. An `fmpz_mpoly_q_t` is a fraction $N/D$ of two polynomials in
$\Z[x_0, \ldots, x_{n-1}]$ (see [](fmpz_mpoly)), kept in canonical form:
$N$ and $D$ are coprime in $\Z[x_0, \ldots, x_{n-1}]$ (no common
polynomial factor and no common integer content), and the leading
coefficient of $D$ is positive. Zero is $0/1$. Canonical form makes the
representation unique, so equality is structural.

All functions assume that their inputs are canonical and produce canonical
output; this is not checked. Values built directly through
`fmpz_mpoly_q_numref` and `fmpz_mpoly_q_denref` must be passed through
`fmpz_mpoly_q_canonicalise`. Arithmetic needs polynomial GCDs
(`fmpz_mpoly_gcd`, see [](fmpz_mpoly/gcd)), which dominate the cost.

The second half of the module, `fmpz_mpoly_q_factored_t`, keeps the
denominator as a product of irreducible factors, which replaces most GCD
computations by trial divisions. Both types are available as generic rings
through `gr_ctx_init_fmpz_mpoly_q` and `gr_ctx_init_fmpz_mpoly_q_factored`.

There is no separate context type: the functions take the
`fmpz_mpoly_ctx_t` of the numerator and denominator, which fixes the
number of variables and the monomial ordering.

```c example name=fmpz_mpoly_q_factored_sum
#include <flint/flint.h>
#include <flint/fmpz_mpoly_q.h>

int main(void)
{
    fmpz_mpoly_ctx_t ctx;
    fmpz_mpoly_q_t a, b, s;
    fmpz_mpoly_q_factored_t fa, fb, fs;
    const char * vars[] = {"x", "y"};

    fmpz_mpoly_ctx_init(ctx, 2, ORD_LEX);
    fmpz_mpoly_q_init(a, ctx);
    fmpz_mpoly_q_init(b, ctx);
    fmpz_mpoly_q_init(s, ctx);
    fmpz_mpoly_q_factored_init(fa, ctx);
    fmpz_mpoly_q_factored_init(fb, ctx);
    fmpz_mpoly_q_factored_init(fs, ctx);

    fmpz_mpoly_q_set_str_pretty(a, "1/(x^2 - y^2)", vars, ctx);
    fmpz_mpoly_q_set_str_pretty(b, "x/(2*(x + y)^2)", vars, ctx);
    fmpz_mpoly_q_add(s, a, b, ctx);
    flint_printf("plain:    ");
    fmpz_mpoly_q_print_pretty(s, vars, ctx);
    flint_printf("\n");

    fmpz_mpoly_q_factored_set_fmpz_mpoly_q(fa, a, ctx);
    fmpz_mpoly_q_factored_set_fmpz_mpoly_q(fb, b, ctx);
    fmpz_mpoly_q_factored_add(fs, fa, fb, ctx);
    flint_printf("factored: ");
    fmpz_mpoly_q_factored_print_pretty(fs, vars, ctx);
    flint_printf("\n");

    fmpz_mpoly_q_factored_get_fmpz_mpoly_q(a, fs, ctx);
    flint_printf("equal: %d\n", fmpz_mpoly_q_equal(a, s, ctx));

    fmpz_mpoly_q_factored_clear(fa, ctx);
    fmpz_mpoly_q_factored_clear(fb, ctx);
    fmpz_mpoly_q_factored_clear(fs, ctx);
    fmpz_mpoly_q_clear(a, ctx);
    fmpz_mpoly_q_clear(b, ctx);
    fmpz_mpoly_q_clear(s, ctx);
    fmpz_mpoly_ctx_clear(ctx);
    return 0;
}
```

```output
plain:    (x^2-x*y+2*x+2*y)/(2*x^3+2*x^2*y-2*x*y^2-2*y^3)
factored: (x^2-x*y+2*x+2*y)/(2*(x-y)*(x+y)^2)
equal: 1
```

## Types and macros {#types}

::: type
fmpz_mpoly_q_struct
fmpz_mpoly_q_t

A rational function, with fields `num` and `den` of type
`fmpz_mpoly_struct`.
:::

::: macro
fmpz_mpoly_q_numref(x)
fmpz_mpoly_q_denref(x)

Pointers to the numerator and denominator of `x`, usable as
`fmpz_mpoly_t`.
:::

## Memory management {#memory}

::: function
void fmpz_mpoly_q_init(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)

Initialises `res` and sets it to zero.
:::

::: function
void fmpz_mpoly_q_clear(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)

Frees the memory used by `res`.
:::

## Assignment {#assignment}

::: function
void fmpz_mpoly_q_swap(fmpz_mpoly_q_t x, fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)

Swaps $x$ and $y$ in constant time.
:::

::: function
void fmpz_mpoly_q_set(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_set_fmpq(fmpz_mpoly_q_t res, const fmpq_t x, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_set_fmpz(fmpz_mpoly_q_t res, const fmpz_t x, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_set_si(fmpz_mpoly_q_t res, slong x, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x$.
:::

## Canonicalisation {#canonical}

::: function
void fmpz_mpoly_q_canonicalise(fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

Puts $x$ in canonical form by dividing the numerator and the denominator
by their GCD (as polynomials over $\Z$, including the integer content)
and making the leading coefficient of the denominator positive. The terms
of both must be sorted and the denominator nonzero.
:::

::: function
int fmpz_mpoly_q_is_canonical(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

Returns whether $x$ is canonical. Besides coprimality and the sign of the
leading coefficient of the denominator, this checks that the denominator
is nonzero and that both polynomials are canonical `fmpz_mpoly` values,
which makes it useful as a consistency check in tests.
:::

## Properties {#properties}

::: function
int fmpz_mpoly_q_is_zero(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_q_is_one(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

Return whether $x$ is 0, respectively 1.
:::

::: function
int fmpz_mpoly_q_is_fmpz(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_q_is_fmpq(const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

Return whether $x$ is an integer constant, respectively a rational
constant.
:::

::: function
void fmpz_mpoly_q_used_vars(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_used_vars_num(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_used_vars_den(int * used, const fmpz_mpoly_q_t f, const fmpz_mpoly_ctx_t ctx)

Set `used[i]` to a nonzero value if the variable of index $i$ occurs in
$f$ (respectively in its numerator or its denominator), and to zero
otherwise.
:::

## Special values {#special}

::: function
void fmpz_mpoly_q_zero(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_one(fmpz_mpoly_q_t res, const fmpz_mpoly_ctx_t ctx)

Set `res` to 0, respectively 1.
:::

::: function
void fmpz_mpoly_q_gen(fmpz_mpoly_q_t res, slong i, const fmpz_mpoly_ctx_t ctx)

Sets `res` to the variable of index $i$, where $0 \le i < n$.
:::

## Input and output {#io}

The variable names `x` are indexed by variable number; if `x` is `NULL`,
the names `x1`, `x2`, ... are used.

::: function
void fmpz_mpoly_q_print_pretty(const fmpz_mpoly_q_t f, const char ** x, const fmpz_mpoly_ctx_t ctx)
char * fmpz_mpoly_q_get_str_pretty(const fmpz_mpoly_q_t f, const char ** x, const fmpz_mpoly_ctx_t ctx)

Print $f$ to standard output, or return a string to be freed with
`flint_free`, in the form `(num)/(den)`, `(num)/c` for a constant
denominator $c$, or just the numerator if the denominator is 1.
:::

::: function
int fmpz_mpoly_q_set_str_pretty(fmpz_mpoly_q_t res, const char * s, const char ** x, fmpz_mpoly_ctx_t ctx)

Sets `res` to the rational function given by the null-terminated string
`s`, which may use integers, the variable names, parentheses, `+`, `-`,
`*`, `/` and `^` followed by an integer exponent. Returns 0 on success and $-1$ if parsing fails or a division by
zero occurs, in which case `res` is unspecified.
:::

## Random generation {#random}

::: function
void fmpz_mpoly_q_randtest(fmpz_mpoly_q_t res, flint_rand_t state, slong length, ulong coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx)

Sets `res` to a random rational function whose numerator and denominator
(before canonicalisation) have at most `length` terms, coefficients of at
most `coeff_bits` bits and exponents less than `exp_bound`.
:::

## Comparison {#comparison}

::: function
int fmpz_mpoly_q_equal(const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)

Returns whether $x = y$.
:::

## Arithmetic {#arithmetic}

Division by zero throws `FLINT_DIVZERO`.

::: function
void fmpz_mpoly_q_neg(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

Sets `res` to $-x$.
:::

::: function
void fmpz_mpoly_q_add(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_add_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_add_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_add_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_sub(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_sub_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_sub_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x + y$, respectively $x - y$. The sum $a/b + c/d$ is
formed as in the classical algorithm for fractions: with $g = \gcd(b,
d)$, the numerator $a (d/g) + c (b/g)$ only needs to be reduced by its GCD
with $g$.
:::

::: function
void fmpz_mpoly_q_mul(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_mul_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_mul_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x y$, cancelling the GCDs of each numerator with the other
denominator before multiplying.
:::

::: function
void fmpz_mpoly_q_div(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_div_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpq_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_div_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_div_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, slong y, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x / y$ for nonzero $y$.
:::

::: function
void fmpz_mpoly_q_inv(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)

Sets `res` to $1/x$ for nonzero $x$.
:::

## Content {#content}

::: function
void fmpz_mpoly_q_content(fmpq_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_q_content(fmpz_t num, fmpz_t den, const fmpz_mpoly_t xnum, const fmpz_mpoly_t xden, const fmpz_mpoly_ctx_t ctx)

Set `res` (or `num`/`den`) to the content of the numerator divided by the
content of the denominator, the positive rational $c$ such that $x/c$ has
primitive numerator and denominator. For zero the result is 1.
:::

## Evaluation {#evaluation}

::: function
void fmpz_mpoly_q_evaluate_acb(acb_t res, const fmpz_mpoly_q_t f, acb_srcptr x, slong prec, const fmpz_mpoly_ctx_t ctx)

Sets `res` to an enclosure of $f$ evaluated at the complex balls $x_i$
(`x` has one entry per variable), at working precision `prec`. If the
enclosure of the denominator contains zero, `res` is set to an
indeterminate value.
:::

## Rational functions with factored denominators {#factored}

An `fmpz_mpoly_q_factored_t` represents an element of $\Q(x_0, \ldots,
x_{n-1})$ as $N/D$ with an expanded numerator $N$ (an `fmpz_mpoly_t`) and
a denominator $D = c \, p_1^{e_1} \cdots p_k^{e_k}$ stored in factored
form as an `fmpz_mpoly_factor_t` (see [](fmpz_mpoly_factor)).

This representation pays off when denominators are products of a modest
number of factors that occur repeatedly. Since the
irreducible factors of the denominators are known, addition, subtraction
and multiplication only need trial divisions by these factors instead of
GCDs, and denominators are never expanded. In exchange, division and
inversion need the factorisation of the numerator of the divisor.

The canonical form is:

- the constant $c$ is a positive integer, and the field `constant_den` of
  the denominator is 1;
- each $p_i$ is a nonconstant irreducible polynomial over $\Z$ (hence
  primitive) with positive leading coefficient, and each $e_i$ is
  positive;
- the $p_i$ are distinct and sorted in increasing order by
  `fmpz_mpoly_cmp`;
- $N$ is not divisible by any $p_i$ and its content is coprime to $c$, so
  that $N$ and the expanded $D$ are coprime;
- zero is represented with $N = 0$ and $D = 1$.

The canonical form is unique, so equality is structural. All functions
assume canonical inputs and produce canonical output, without checking;
values built directly must be passed through
`fmpz_mpoly_q_factored_canonicalise`. Unless stated otherwise, outputs
may alias inputs.

::: type
fmpz_mpoly_q_factored_struct
fmpz_mpoly_q_factored_t

A rational function with fields `num` (an `fmpz_mpoly_struct`) and `den`
(an `fmpz_mpoly_factor_struct`).
:::

::: macro
fmpz_mpoly_q_factored_numref(x)
fmpz_mpoly_q_factored_denref(x)

Pointers to the numerator of `x`, usable as an `fmpz_mpoly_t`, and to its
denominator, usable as an `fmpz_mpoly_factor_t`.
:::

::: function
void fmpz_mpoly_q_factored_init(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_clear(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)

Initialise `res` to zero, and free its memory.
:::

::: function
void fmpz_mpoly_q_factored_swap(fmpz_mpoly_q_factored_t x, fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

Swaps $x$ and $y$ in constant time.
:::

::: function
void fmpz_mpoly_q_factored_set(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_set_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_t x, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_set_fmpq(fmpz_mpoly_q_factored_t res, const fmpq_t x, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_set_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_t x, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_set_si(fmpz_mpoly_q_factored_t res, slong x, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x$.
:::

::: function
void fmpz_mpoly_q_factored_zero(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_one(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_ctx_t ctx)

Set `res` to 0, respectively 1.
:::

::: function
void fmpz_mpoly_q_factored_gen(fmpz_mpoly_q_factored_t res, slong i, const fmpz_mpoly_ctx_t ctx)

Sets `res` to the variable of index $i$, where $0 \le i < n$.
:::

::: function
void fmpz_mpoly_q_factored_set_fmpz_mpoly_q(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_get_fmpz_mpoly_q(fmpz_mpoly_q_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

Convert from an `fmpz_mpoly_q_t`, which factors its denominator, and to
one, which expands the denominator.
:::

::: function
void fmpz_mpoly_q_factored_canonicalise(fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

Puts $x$ in canonical form. The input may have any nonzero integer
constant and `constant_den` in the denominator, and bases that are any
nonzero polynomials (not necessarily irreducible, distinct or sorted)
with any integer exponents; bases with negative exponents are moved to
the numerator. Each base with positive exponent is factored into
irreducibles, and common factors of numerator and denominator are then
removed. Throws if the denominator is zero.
:::

::: function
int fmpz_mpoly_q_factored_is_canonical(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

Returns whether $x$ is canonical. This factors every base of the
denominator to check irreducibility and is meant for tests.
:::

::: function
int fmpz_mpoly_q_factored_is_zero(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_q_factored_is_one(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

Return whether $x$ is 0, respectively 1.
:::

::: function
int fmpz_mpoly_q_factored_equal(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

Returns whether $x = y$.
:::

::: function
void fmpz_mpoly_q_factored_print_pretty(const fmpz_mpoly_q_factored_t x, const char ** vars, const fmpz_mpoly_ctx_t ctx)
char * fmpz_mpoly_q_factored_get_str_pretty(const fmpz_mpoly_q_factored_t x, const char ** vars, const fmpz_mpoly_ctx_t ctx)

Print $x$ to standard output, or return a string to be freed with
`flint_free`, with the denominator in factored form, for example
`(x1+1)/(2*(x2)*(x1-x2)^3)`. If `vars` is `NULL`, the names `x1`, `x2`,
... are used.
:::

::: function
int fmpz_mpoly_q_factored_set_str_pretty(fmpz_mpoly_q_factored_t res, const char * s, const char ** vars, const fmpz_mpoly_ctx_t ctx)

Sets `res` to the rational function in the string `s`, with the syntax of
`fmpz_mpoly_q_set_str_pretty`; this parses the output of both
`get_str_pretty` functions. The expression is evaluated with the
arithmetic of this type, so each division by a nonconstant expression
factors it. Returns 0 on success and $-1$ if parsing fails or a division by
zero occurs, in which case `res` is unspecified.
:::

::: function
void fmpz_mpoly_q_factored_randtest(fmpz_mpoly_q_factored_t res, flint_rand_t state, slong length, ulong coeff_bits, slong exp_bound, const fmpz_mpoly_ctx_t ctx)

Sets `res` to a random rational function whose numerator and the random
polynomials multiplied to form the denominator have at most `length`
terms, coefficients of at most `coeff_bits` bits and exponents less than
`exp_bound`.
:::

::: function
void fmpz_mpoly_q_factored_neg(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

Sets `res` to $-x$.
:::

::: function
void fmpz_mpoly_q_factored_add(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_sub(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x + y$, respectively $x - y$. The common denominator takes
the larger exponent of each irreducible factor. Only factors with the same
exponent in both denominators can divide the new numerator; they are
removed by trial division.
:::

::: function
void fmpz_mpoly_q_factored_mul(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

Sets `res` to $x y$. The factors of the denominator of $x$ that are not in
the denominator of $y$ are removed from the numerator of $y$ by trial
division, and vice versa.
:::

::: function
void fmpz_mpoly_q_factored_div(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, const fmpz_mpoly_ctx_t ctx)

Sets `res` to $x / y$, which factors the numerator of $y$. Throws if $y$
is zero.
:::

::: function
void fmpz_mpoly_q_factored_inv(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_ctx_t ctx)

Sets `res` to $1/x$, which factors the numerator of $x$ and expands its
denominator. Throws if $x$ is zero.
:::

::: function
void fmpz_mpoly_q_factored_pow_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_pow_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong e, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x^e$ for any integer $e$, with $0^0 = 1$. A negative
exponent factors the numerator of $x$. Throw if $x$ is zero and $e$ is
negative.
:::

::: function
void fmpz_mpoly_q_factored_add_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_add_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_add_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_add_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_sub_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_sub_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_sub_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_sub_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x + c$, respectively $x - c$.
:::

::: function
void fmpz_mpoly_q_factored_mul_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_mul_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_mul_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_mul_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x c$. A scalar only changes the numerator and the constant
of the denominator; for a polynomial, the factors of the denominator of
$x$ are removed from $c$ by trial division.
:::

::: function
void fmpz_mpoly_q_factored_div_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_div_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_div_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_q_factored_div_fmpz_mpoly(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

Set `res` to $x / c$; division by a polynomial factors $c$. Throw if $c$
is zero.
:::

## Internal functions {#internal}

The following functions work on numerators and denominators passed
separately; the inputs are canonical fractions and the outputs must not
alias the inputs unless the implementation allows it.

::: function
void _fmpz_mpoly_q_add(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den, const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den, const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_q_sub(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den, const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den, const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_q_mul(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den, const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den, const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_q_div(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den, const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den, const fmpz_mpoly_t y_num, const fmpz_mpoly_t y_den, const fmpz_mpoly_ctx_t ctx)

Set the canonical fraction `res_num/res_den` to the sum, difference,
product or quotient of `x_num/x_den` and `y_num/y_den`.
:::

::: function
void _fmpz_mpoly_q_add_fmpq(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den, const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den, const fmpz_t y_num, const fmpz_t y_den, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den, const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den, const fmpz_t y_num, const fmpz_t y_den, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_t res_num, fmpz_mpoly_t res_den, const fmpz_mpoly_t x_num, const fmpz_mpoly_t x_den, const fmpz_t y_num, const fmpz_t y_den, const fmpz_mpoly_ctx_t ctx)

The same with the canonical rational number `y_num/y_den` as second
operand.
:::

::: function
void _fmpz_vec_content2(fmpz_t res, const fmpz * vec, slong len, const fmpz_t inp)

Sets `res` to the GCD of $|\mathit{inp}|$ and the `len` entries of
`vec`, stopping early when it reaches 1.
:::

::: function
void fmpz_mpoly_gcd_assert_successful(fmpz_mpoly_t res, const fmpz_mpoly_t x, const fmpz_mpoly_t y, const fmpz_mpoly_ctx_t ctx)

Sets `res` to $\gcd(x, y)$ with `fmpz_mpoly_gcd`, throwing `FLINT_ERROR`
if the GCD computation fails.
:::

::: function
void _fmpz_mpoly_q_mpoly_divexact(fmpz_mpoly_t res, const fmpz_mpoly_t x, const fmpz_mpoly_t y, const fmpz_mpoly_ctx_t ctx)

Sets `res` to $x/y$, assuming that the division is exact, with a fast
path for constant $y$.
:::
