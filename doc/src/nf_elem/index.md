---
title: nf_elem.h
header: nf_elem.h
description: Elements of number fields $\mathbf{Q}[x]/(f)$, with special representations for linear and quadratic fields, arithmetic with delayed reduction, norms, traces, representation matrices and reduction modulo integers.
category: Algebraic numbers and number fields
order: 3
---

This module implements arithmetic with elements of a number field
$K = \mathbf{Q}(\theta) = \mathbf{Q}[x]/(f)$ given by an `nf_t`
(`nf.h`). An element is a polynomial in $\theta$ of degree less than
$d = \deg f$ with rational coefficients, stored as an integer
numerator polynomial and a common positive denominator. Fields of
degree 1 and 2 use a fixed-size representation, and the field
$\mathbf{Q}(i)$ has its own multiplication. Products are reduced modulo
$f$ using precomputed powers of $\theta$; inverses are computed by the
extended Euclidean algorithm, norms as resultants, traces from
precomputed power sums. Every function takes the `nf_t` as its last
argument, and all elements passed to a function must belong to the same
field.

The guide [](nf_elem/getting-started) shows the basic operations. How a
fixed number field compares with the other exact
representations (`qqbar`, `ca`) is discussed in
[](algorithms/exact/overview).

::: summary algorithms/number-theory/number-fields :::

## Types and representation {#types}

::: type
nf_elem_struct
nf_elem_t

An element of a number field: a union of an `fmpq_poly_t` (general
fields), an `lnf_elem_t` (degree 1) and a `qnf_elem_t` (degree 2). Which
member is used is determined by the `flag` of the `nf_t`.
:::

::: type
lnf_elem_struct
lnf_elem_t

An element of a field of degree 1, as a numerator `num` and a
denominator `den`.
:::

::: type
qnf_elem_struct
qnf_elem_t

An element $(a_0 + a_1\theta + a_2\theta^2)/q$ of a field of degree 2,
as an array `num` of three `fmpz` and a denominator `den`. The third
coefficient allows products to be formed without immediate reduction;
in a reduced element it is zero.
:::

::: macro
NF_ELEM(x)
LNF_ELEM(x)
QNF_ELEM(x)

The `fmpq_poly_t`, `lnf_elem_t` and `qnf_elem_t` members of the union.
:::

::: macro
NF_ELEM_NUMREF(x)
NF_ELEM_DENREF(x)
LNF_ELEM_NUMREF(x)
LNF_ELEM_DENREF(x)
QNF_ELEM_NUMREF(x)
QNF_ELEM_DENREF(x)

Pointers to the numerator coefficients and the denominator in each
representation.
:::

An element is *reduced* if its numerator has degree less than $d$ and
*canonical* if moreover the numerator and denominator are coprime (their
content is 1) and the denominator is positive, with denominator 1 for
zero. Functions without a leading underscore require canonical inputs
and produce canonical outputs unless stated otherwise, and allow
aliasing.

## Memory management {#memory}

::: function
void nf_elem_init(nf_elem_t a, const nf_t nf)
void nf_elem_clear(nf_elem_t a, const nf_t nf)

Initialise $a$ to zero in `nf`, respectively clear it.
:::

::: function
void nf_elem_swap(nf_elem_t a, nf_elem_t b, const nf_t nf)

Swaps $a$ and $b$ in constant time.
:::

::: function
void nf_elem_randtest(nf_elem_t a, flint_rand_t state, flint_bitcnt_t bits, const nf_t nf)
void nf_elem_randtest_not_zero(nf_elem_t a, flint_rand_t state, flint_bitcnt_t bits, const nf_t nf)

Set $a$ to a random canonical element, respectively a random nonzero
one, whose numerator coefficients and denominator have at most `bits`
bits.
:::

## Canonicalisation and reduction {#reduction}

::: function
void nf_elem_canonicalise(nf_elem_t a, const nf_t nf)

Removes the common content of numerator and denominator and makes the
denominator positive (1 for zero). Does not reduce modulo $f$.
:::

::: function
void _nf_elem_reduce(nf_elem_t a, const nf_t nf)

Reduces the numerator of $a$ modulo the defining polynomial, without
canonicalising. Used after products with delayed reduction
(`nf_elem_mul_red` with `red` $= 0$).
:::

::: function
void nf_elem_reduce(nf_elem_t a, const nf_t nf)

Reduces $a$ modulo the defining polynomial and canonicalises it.
:::

::: function
int _nf_elem_invertible_check(nf_elem_t a, const nf_t nf)

Returns whether $a$ is invertible modulo the defining polynomial, that
is, whether $\gcd(a, f) = 1$ in $\mathbf{Q}[x]$. This matters only when
$f$ is reducible, which is not excluded by `nf_init`; the function is
meant for test code.
:::

## Assignment and conversion {#conversion}

::: function
void nf_elem_zero(nf_elem_t a, const nf_t nf)
void nf_elem_one(nf_elem_t a, const nf_t nf)

Set $a$ to 0, respectively 1.
:::

::: function
void nf_elem_gen(nf_elem_t a, const nf_t nf)

Sets $a$ to the generator $\theta$, the class of $x$. In a field of
degree 1 this is the root $-f_0/f_1$ of $f = f_1 x + f_0$.
:::

::: function
void nf_elem_set(nf_elem_t a, const nf_elem_t b, const nf_t nf)

Sets $a$ to $b$.
:::

::: function
void nf_elem_set_si(nf_elem_t a, slong c, const nf_t nf)
void nf_elem_set_ui(nf_elem_t a, ulong c, const nf_t nf)
void nf_elem_set_fmpz(nf_elem_t a, const fmpz_t c, const nf_t nf)
void nf_elem_set_fmpq(nf_elem_t a, const fmpq_t c, const nf_t nf)

Set $a$ to the rational number $c$.
:::

::: function
void nf_elem_set_fmpq_poly(nf_elem_t a, const fmpq_poly_t pol, const nf_t nf)

Sets $a$ to the image of the canonical polynomial `pol` in $K$,
reducing it modulo $f$ if its degree is at least $d$.
:::

::: function
void nf_elem_get_fmpq_poly(fmpq_poly_t pol, const nf_elem_t a, const nf_t nf)

Sets `pol` to the polynomial in $x$ representing the reduced element $a$.
:::

::: function
void nf_elem_set_fmpz_mat_row(nf_elem_t b, const fmpz_mat_t M, const slong i, fmpz_t den, const nf_t nf)

Sets $b$ to the element whose numerator coefficients are the entries of
row $i$ of $M$ (column 0 holding the constant coefficient) and whose
denominator is `den`, and canonicalises it.
:::

::: function
void nf_elem_get_fmpz_mat_row(fmpz_mat_t M, const slong i, fmpz_t den, const nf_elem_t b, const nf_t nf)

Sets row $i$ of $M$, which must have $d$ columns, to the numerator
coefficients of $b$ (column 0 holding the constant coefficient) and
`den` to its denominator.
:::

::: function
void _nf_elem_get_nmod_poly(nmod_poly_t pol, const nf_elem_t a, const nf_t nf)

Sets `pol` to the numerator of $a$ reduced modulo the modulus of `pol`,
ignoring the denominator.
:::

::: function
void nf_elem_get_nmod_poly_den(nmod_poly_t pol, const nf_elem_t a, const nf_t nf, int den)
void nf_elem_get_nmod_poly(nmod_poly_t pol, const nf_elem_t a, const nf_t nf)

Set `pol` to the reduction of the numerator of $a$ modulo the modulus
$n$ of `pol`; if `den` is nonzero (always for `nf_elem_get_nmod_poly`),
multiply it by the inverse of the denominator modulo $n$, which must
exist.
:::

::: function
void _nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf, const fmpz_mod_ctx_t ctx)
void nf_elem_get_fmpz_mod_poly_den(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf, int den, const fmpz_mod_ctx_t ctx)
void nf_elem_get_fmpz_mod_poly(fmpz_mod_poly_t pol, const nf_elem_t a, const nf_t nf, const fmpz_mod_ctx_t ctx)

Versions of the three functions above for a multiword modulus given by
`ctx`.
:::

## Coefficients and denominator {#coefficients}

::: function
void nf_elem_get_den(fmpz_t d, const nf_elem_t b, const nf_t nf)
void nf_elem_set_den(nf_elem_t b, fmpz_t d, const nf_t nf)

Get and set the denominator of $b$. `nf_elem_set_den` requires $d > 0$
and does not canonicalise.
:::

::: function
int nf_elem_den_is_one(const nf_elem_t a, const nf_t nf)

Returns whether the denominator of $a$ is 1, that is, whether the
canonical $a$ lies in $\Z[\theta]$.
:::

::: function
void nf_elem_get_coeff_fmpq(fmpq_t a, const nf_elem_t b, slong i, const nf_t nf)

Sets $a$ to the coefficient of $\theta^i$ in $b$, a rational number.
:::

::: function
void nf_elem_get_coeff_fmpz(fmpz_t a, const nf_elem_t b, slong i, const nf_t nf)

Sets $a$ to the coefficient of $\theta^i$ in the numerator of $b$,
ignoring the denominator.
:::

::: function
void _nf_elem_set_coeff_num_fmpz(nf_elem_t a, slong i, const fmpz_t b, const nf_t nf)

Sets the coefficient of $\theta^i$ in the numerator of $a$ to $b$,
keeping the denominator, and canonicalises. Allows $i$ up to
$2(d - 1)$, so that the result may be unreduced; throws for larger $i$.
:::

## Comparison {#comparison}

::: function
int _nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf)

Returns whether $a = b$ for reduced but not necessarily canonical $a$
and $b$, by cross-multiplication, with a cheap test of the bit sizes
first.
:::

::: function
int nf_elem_equal(const nf_elem_t a, const nf_elem_t b, const nf_t nf)

Returns whether $a = b$ for canonical $a$ and $b$, by comparing the
representations.
:::

::: function
int nf_elem_equal_si(const nf_elem_t a, const slong b, const nf_t nf)
int nf_elem_equal_ui(const nf_elem_t a, const ulong b, const nf_t nf)
int nf_elem_equal_fmpz(const nf_elem_t a, const fmpz_t b, const nf_t nf)
int nf_elem_equal_fmpq(const nf_elem_t a, const fmpq_t b, const nf_t nf)

Return whether the canonical element $a$ equals the rational number $b$.
:::

::: function
int nf_elem_is_zero(const nf_elem_t a, const nf_t nf)
int nf_elem_is_one(const nf_elem_t a, const nf_t nf)

Return whether $a$ is 0, respectively 1.
:::

::: function
int nf_elem_is_gen(const nf_elem_t a, const nf_t nf)

Returns whether $a$ equals the generator $\theta$. In a field of
degree 1 with defining polynomial $f = f_1 x + f_0$, the generator is
the rational number $-f_0/f_1$.
:::

::: function
int nf_elem_is_integer(const nf_elem_t a, const nf_t nf)
int nf_elem_is_rational(const nf_elem_t a, const nf_t nf)

Return whether the canonical element $a$ lies in $\Z$, respectively in
$\mathbf{Q}$.
:::

## Input and output {#io}

::: function
void nf_elem_print_pretty(const nf_elem_t a, const nf_t nf, const char * var)

Prints $a$ to `stdout` as a polynomial in the variable `var`, a
nonempty string.
:::

::: function
char * nf_elem_get_str_pretty(const nf_elem_t a, const char * var, const nf_t nf)

Returns the string printed by `nf_elem_print_pretty`, allocated with
`flint_malloc`.
:::

## Arithmetic {#arithmetic}

::: function
void nf_elem_neg(nf_elem_t a, const nf_elem_t b, const nf_t nf)

Sets $a$ to $-b$.
:::

::: function
void _nf_elem_add(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
void _nf_elem_sub(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)

Set $a$ to $b + c$, respectively $b - c$, without canonicalising.
:::

::: function
void nf_elem_add(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
void nf_elem_sub(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)

Set $a$ to $b + c$, respectively $b - c$.
:::

::: function
void nf_elem_add_si(nf_elem_t a, const nf_elem_t b, slong c, const nf_t nf)
void nf_elem_add_fmpz(nf_elem_t a, const nf_elem_t b, const fmpz_t c, const nf_t nf)
void nf_elem_add_fmpq(nf_elem_t a, const nf_elem_t b, const fmpq_t c, const nf_t nf)
void nf_elem_sub_si(nf_elem_t a, const nf_elem_t b, slong c, const nf_t nf)
void nf_elem_sub_fmpz(nf_elem_t a, const nf_elem_t b, const fmpz_t c, const nf_t nf)
void nf_elem_sub_fmpq(nf_elem_t a, const nf_elem_t b, const fmpq_t c, const nf_t nf)
void nf_elem_si_sub(nf_elem_t a, slong c, const nf_elem_t b, const nf_t nf)
void nf_elem_fmpz_sub(nf_elem_t a, const fmpz_t c, const nf_elem_t b, const nf_t nf)
void nf_elem_fmpq_sub(nf_elem_t a, const fmpq_t c, const nf_elem_t b, const nf_t nf)

Set $a$ to $b + c$, $b - c$, respectively $c - b$, for a rational
number $c$.
:::

::: function
void nf_elem_scalar_mul_si(nf_elem_t a, const nf_elem_t b, slong c, const nf_t nf)
void nf_elem_scalar_mul_fmpz(nf_elem_t a, const nf_elem_t b, const fmpz_t c, const nf_t nf)
void nf_elem_scalar_mul_fmpq(nf_elem_t a, const nf_elem_t b, const fmpq_t c, const nf_t nf)

Set $a$ to $cb$.
:::

::: function
void nf_elem_scalar_div_si(nf_elem_t a, const nf_elem_t b, slong c, const nf_t nf)
void nf_elem_scalar_div_fmpz(nf_elem_t a, const nf_elem_t b, const fmpz_t c, const nf_t nf)
void nf_elem_scalar_div_fmpq(nf_elem_t a, const nf_elem_t b, const fmpq_t c, const nf_t nf)

Set $a$ to $b/c$ for $c \ne 0$.
:::

::: function
void nf_elem_mul_gen(nf_elem_t a, const nf_elem_t b, const nf_t nf)

Sets $a$ to $\theta b$, by a shift of the coefficients followed by one
reduction step.
:::

::: function
void _nf_elem_mul(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
void _nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf, int red)

Set $a$ to $bc$ without canonicalising, reducing modulo $f$ only if
`red` is nonzero (always for `_nf_elem_mul`). The inputs must be
reduced, and $a$ must not be aliased with $b$ or $c$.
:::

::: function
void nf_elem_mul(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
void nf_elem_mul_red(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf, int red)

Set $a$ to $bc$, canonicalised. With `red` $= 0$ the product is not
reduced modulo $f$, which allows sums of products to be reduced once
with `nf_elem_reduce`. The inputs must be reduced.

::: details
The numerators are multiplied in $\Z[x]$ and the denominators in $\Z$.
For $d \ge 3$ the product of length up to $2d - 1$ is reduced with the
precomputed powers $\theta^d, \ldots, \theta^{2d-2}$ if
$d < 30$, and by division with remainder otherwise; for $d = 2$ one
step $\theta^2 = -(f_1\theta + f_0)/f_2$ suffices, and for
$\mathbf{Q}(i)$ a dedicated formula with a gcd for the denominator is
used ([](algorithms/number-theory/number-fields#mul)).
:::
:::

::: function
void _nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf)
void nf_elem_inv(nf_elem_t a, const nf_elem_t b, const nf_t nf)

Set $a$ to $b^{-1}$ for $b \ne 0$, by the extended Euclidean algorithm
for $b$ and $f$ over $\mathbf{Q}$. `_nf_elem_inv` does not allow
aliasing. The result is undefined for $b = 0$.
:::

::: function
void _nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
void nf_elem_div(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)

Set $a$ to $b/c$ for $c \ne 0$, as $b \cdot c^{-1}$. `_nf_elem_div`
does not canonicalise and does not allow $a$ to be aliased with $b$.
:::

::: function
void _nf_elem_pow(nf_elem_t res, const nf_elem_t a, ulong e, const nf_t nf)

Sets `res` to $a^e$ by left-to-right binary exponentiation
[@Knu1997, §4.6.3], for $a \ne 0$ and $e \ge 2$. Does not allow
aliasing.
:::

::: function
void nf_elem_pow(nf_elem_t res, const nf_elem_t a, ulong e, const nf_t nf)

Sets `res` to $a^e$, with $a^0 = 1$ for all $a$ including zero.
:::

## Norm, trace and representation matrix {#norm-trace}

::: function
void _nf_elem_norm(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf)
void nf_elem_norm(fmpq_t res, const nf_elem_t a, const nf_t nf)

Set the rational number `rnum/rden`, respectively `res`, to the absolute
norm $N_{K/\mathbf{Q}}(a)$, the product of the conjugates of $a$. It is
computed as $\operatorname{Res}(f, a)$ divided by the appropriate power
of the leading coefficient of $f$; the resultant of the primitive parts
is computed by `_fmpz_poly_resultant`, which uses the multimodular
algorithm for large inputs and a pseudo-remainder (subresultant-type) Euclidean algorithm otherwise
([](algorithms/number-theory/number-fields#norm)).
:::

::: function
void _nf_elem_norm_div(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf, const fmpz_t divisor, slong nbits)
void nf_elem_norm_div(fmpq_t res, const nf_elem_t a, const nf_t nf, const fmpz_t divisor, slong nbits)

Set the result to $N_{K/\mathbf{Q}}(a)/\mathtt{divisor}$, assuming that
the quotient of the resultant by `divisor` is an integer of at most
`nbits` bits; the multimodular resultant
(`_fmpz_poly_resultant_modular_div`) then needs only enough primes for
the quotient.
For `divisor` $= 0$ the result is 0.
:::

::: function
void _nf_elem_trace(fmpz_t rnum, fmpz_t rden, const nf_elem_t a, const nf_t nf)
void nf_elem_trace(fmpq_t res, const nf_elem_t a, const nf_t nf)

Set the result to the absolute trace $\operatorname{Tr}_{K/\mathbf{Q}}(a)$,
computed as $\sum_i a_i S_i$ from the coefficients $a_i$ of $a$ and the
power sums $S_i$ of the roots of $f$ stored in the `nf_t`.
:::

::: function
void nf_elem_rep_mat(fmpq_mat_t res, const nf_elem_t a, const nf_t nf)

Sets the $d \times d$ matrix `res` to the matrix of multiplication by
$a$ with respect to the basis $1, \theta, \ldots, \theta^{d-1}$, in the
row convention: row $j$ holds the coordinates of $\theta^j a$.
:::

::: function
void nf_elem_rep_mat_fmpz_mat_den(fmpz_mat_t res, fmpz_t den, const nf_elem_t a, const nf_t nf)

Sets `res` and `den` so that `res`/`den` is the matrix of
`nf_elem_rep_mat`, with `res` primitive, that is, with content 1.
:::

## Reduction modulo integers {#modular}

These functions replace an element by a representative with reduced
coefficients modulo an integer $m$, as needed by multimodular and
$p$-adic algorithms over number fields. Let $q$ be the denominator of
$a$.

::: function
void _nf_elem_mod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf, int sign)

Sets `res` to the integral element whose coefficients are those of the
numerator of $a$ reduced modulo $m$ = `mod`, ignoring the denominator,
in $[0, m)$ if `sign` is zero and in the symmetric range otherwise.
:::

::: function
void nf_elem_mod_fmpz_den(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf, int den)
void nf_elem_smod_fmpz_den(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf, int den)

If `den` is zero, set `res` to the integral element with the
coefficients of $qa$ reduced modulo $m$. If `den` is nonzero, set `res`
to $z/q$, where $z$ is the numerator of $a$ reduced modulo $mq$, so
that $a - \mathtt{res}$ is integral with coefficients divisible by $m$.
The `mod` versions use the residues in $[0, m)$ (or $[0, mq)$), the
`smod` versions the symmetric residues.
:::

::: function
void nf_elem_mod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
void nf_elem_smod_fmpz(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)

The case `den` $= 1$ of the functions above.
:::

::: function
void nf_elem_coprime_den(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)
void nf_elem_coprime_den_signed(nf_elem_t res, const nf_elem_t a, const fmpz_t mod, const nf_t nf)

Set `res` to an element congruent to $a$ modulo $m$ whose denominator
divides $m^\infty$: writing $q = q_1 q_2$ with $q_1$ composed of the
primes dividing $m$ and $q_2$ coprime to $m$, the numerator of $a$ is
multiplied by $q_2^{-1}$ modulo $m q_1$, reduced into $[0, m q_1)$
(or the symmetric range for the signed version) and divided by $q_1$.
Then $a - \mathtt{res}$ has denominator coprime to $m$ and numerator
coefficients divisible by $m$.
:::

## Internal functions {#internal}

::: function
void _nf_elem_add_lf(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf, int can)
void _nf_elem_sub_lf(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf, int can)

Addition and subtraction in a field of degree 1, canonicalising the
result if `can` is nonzero.
:::

::: function
void _nf_elem_add_qf(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf, int can)
void _nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf, int can)
void nf_elem_add_qf(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)
void nf_elem_sub_qf(nf_elem_t a, const nf_elem_t b, const nf_elem_t c, const nf_t nf)

Addition and subtraction in a field of degree 2, on all three numerator
coefficients (so that unreduced elements can be added), canonicalising
the result if `can` is nonzero (always for the versions without
underscore). The underscore versions do not allow $a$ to be aliased
with $c$.
:::
