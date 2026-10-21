---
title: fq.h
header: fq.h
template-headers: fq_templates
description: Finite fields $\mathbb{F}_{p^d}$ of arbitrary characteristic, with elements stored as polynomials over $\Z/p\Z$ with multiprecision coefficients.
category: Finite fields
order: 5
---

This module implements arithmetic in a finite field
$\mathbb{F}_q = \mathbb{F}_p[x]/(f)$, $q = p^d$, for a prime $p$ of any
size. An element is stored as its unique representative of degree less
than $d$, an `fmpz_poly_t` whose coefficients lie in $[0, p)$. A product
is computed as an integer polynomial product followed by reduction
modulo $f$, which is fast for the sparse defining polynomials that FLINT
chooses by default ([](algorithms/finite-fields/reduction)). The modules
`fq_vec`, `fq_mat`, `fq_poly`, `fq_poly_factor` and `fq_embed` provide
vectors, matrices, polynomials, factorisation and embeddings over these
fields.

The interface is generated, like that of `fq_nmod` (word-size $p$) and
`fq_zech` (Zech logarithms), from the template header `fq_templates.h`;
the three modules have the same functions apart from the conversions and
a few context functions. The common interface is documented here in
detail. Most users should use `fq_default` (see
[](fq_default/getting-started)), which chooses the best of the three
representations, and use this module directly only for fields with a
characteristic that does not fit in a word, or for its low-level
functions.

::: summary algorithms/finite-fields/overview :::

## Types and conventions {#types}

::: type
fq_ctx_struct
fq_ctx_t

A finite field context. It holds an `fmpz_mod_ctx_t` for $p$ (`ctxp`),
the defining polynomial $f$ as an `fmpz_mod_poly_t` (`modulus`), the
power series inverse of its reversal (`inv`, used for reduction by Newton
division), the nonzero coefficients of the monic multiple of $f$ and
their exponents (`a`, `j`, `len`, used for sparse reduction), a flag
`sparse_modulus` set if $f$ has at most five nonzero terms, a flag
`is_conway` set if $f$ is a Conway polynomial from FLINT's table, and the
name of the generator (`var`). `fq_ctx_t` is an array of one structure.
:::

::: type
fq_struct
fq_t

An element of $\mathbb{F}_q$: an `fmpz_poly_struct` holding the
representative of degree less than $d$, with coefficients in $[0, p)$,
normalised (the leading coefficient is nonzero, and zero has length 0).
`fq_t` is an array of one structure.
:::

- All functions take the context as their last argument, even those that
  do not use it.
- Inputs must be reduced representatives and outputs are reduced, except
  for the functions with a leading underscore, which work on coefficient
  arrays and state their conditions.
- Outputs may be the same objects as inputs (aliasing), unless stated
  otherwise.
- The defining polynomial must be irreducible and $p$ must be prime; this
  is not checked. Functions that divide throw `FLINT_ERROR` when asked to
  invert zero.

## Context management {#context}

::: function
void fq_ctx_init(fq_ctx_t ctx, const fmpz_t p, slong d, const char * var)

Initialises `ctx` as a context for $\mathbb{F}_{p^d}$, $d \ge 1$, with the
generator printed as `var`, a nonempty null-terminated string. The
defining polynomial is the Conway polynomial for $(p, d)$ if it is in
FLINT's table, and otherwise a random sparse irreducible polynomial
(a trinomial or pentanomial if one is found). $p$ must be prime; this is
not checked.

::: details
The table covers all primes below 110000 for some degrees: degrees up to
409 for $p = 2$ and up to about 200 for $p < 20$, fewer for larger $p$, and
only $d \le 4$ above 11000 (see
[](algorithms/finite-fields/defining-polynomials)). The random choice
uses a fixed random state, so it is reproducible, but it depends on the
FLINT version.
:::
:::

::: function
void fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char * var)
void fq_ctx_init_conway_ui(fq_ctx_t ctx, ulong p, slong d, const char * var)
int _fq_ctx_init_conway(fq_ctx_t ctx, const fmpz_t p, slong d, const char * var)
int _fq_ctx_init_conway_ui(fq_ctx_t ctx, ulong p, slong d, const char * var)

Initialise `ctx` as a context for $\mathbb{F}_{p^d}$ defined by the Conway
polynomial for $(p, d)$. The functions without an underscore throw
`FLINT_ERROR` if the polynomial is not in FLINT's table. The functions
with an underscore return 1 on success and 0, leaving `ctx`
uninitialised, if the polynomial is not in the table.

::: details
A Conway polynomial is primitive, so its root, the generator returned by
`fq_gen`, generates the multiplicative group. Conway polynomials for
$\mathbb{F}_{p^m}$ and $\mathbb{F}_{p^n}$ with $m \mid n$ are compatible:
$x^{(p^n - 1)/(p^m - 1)}$ in the larger field is a root of the
polynomial for the smaller one. The two functions with an underscore do
not set the flag `is_conway` of the context (it only affects how
`fq_sqrt` finds a non-square); the other two do.
:::
:::

::: function
void fq_ctx_init_modulus(fq_ctx_t ctx, const fmpz_mod_poly_t modulus, const fmpz_mod_ctx_t ctxp, const char * var)

Initialises `ctx` as a context for $\mathbb{F}_p[x]/(f)$, where $f$ is
`modulus`, an irreducible polynomial of degree $d \ge 1$ over
$\mathbb{F}_p$, and $p$ is the modulus of `ctxp`. $f$ need not be monic.
The polynomial and `ctxp` are copied. Irreducibility is not checked.
:::

::: function
void fq_ctx_init_randtest(fq_ctx_t ctx, flint_rand_t state, int type)

Initialises `ctx` as a random field for testing: it creates a random
`fq_nmod` context with `fq_nmod_ctx_init_randtest` and copies it, so the
prime fits in a word and the defining polynomial is a Conway polynomial,
possibly multiplied by a random constant (to test non-monic moduli). The
`type` selects the ranges as in `_nmod_poly_conway_rand`: primes of up to
17 bits for `type` 0 and 1, of up to 10 bits for 2 and 3, and degrees
below 15 for odd `type`.
:::

::: function
void fq_ctx_init_randtest_reducible(fq_ctx_t ctx, flint_rand_t state, int type)

Initialises `ctx` with a random monic modulus that is usually reducible,
for testing the functions that also work over $\mathbb{F}_p[x]/(f)$ for
reducible $f$. For `type` 0 the prime has up to 30 bits and the degree
is at most 30, for 1 up to 30 bits and degree 15, for 2 up to 10 bits and
degree 30, for 3 up to 10 bits and degree 15.
:::

::: function
void fq_ctx_clear(fq_ctx_t ctx)

Frees the memory used by `ctx`.
:::

::: function
const fmpz_mod_poly_struct * fq_ctx_modulus(const fq_ctx_t ctx)

Returns a pointer to the defining polynomial $f$, an `fmpz_mod_poly` over
the context `ctx->ctxp`.
:::

::: function
slong fq_ctx_degree(const fq_ctx_t ctx)

Returns the degree $d = [\mathbb{F}_q : \mathbb{F}_p]$.
:::

::: function
const fmpz * fq_ctx_prime(const fq_ctx_t ctx)

Returns a pointer to the characteristic $p$.
:::

::: function
void fq_ctx_order(fmpz_t f, const fq_ctx_t ctx)

Sets `f` to the number of elements $q = p^d$.
:::

::: function
int fq_ctx_fprint(FILE * file, const fq_ctx_t ctx)
void fq_ctx_print(const fq_ctx_t ctx)

Print $p$, $d$ and the defining polynomial in the variable `X` to `file`,
respectively `stdout`, on three lines. `fq_ctx_fprint` returns a positive
value on success.
:::

## Memory management {#memory}

::: function
void fq_init(fq_t rop, const fq_ctx_t ctx)
void fq_init2(fq_t rop, const fq_ctx_t ctx)

Initialise `rop` and set it to zero. `fq_init2` also allocates space for
$d$ coefficients.
:::

::: function
void fq_clear(fq_t rop, const fq_ctx_t ctx)

Frees the memory used by `rop`.
:::

## Assignment and conversion {#assignment}

::: function
void fq_set(fq_t rop, const fq_t op, const fq_ctx_t ctx)

Sets `rop` to `op`.
:::

::: function
void fq_swap(fq_t op1, fq_t op2, const fq_ctx_t ctx)

Swaps `op1` and `op2` efficiently.
:::

::: function
void fq_zero(fq_t rop, const fq_ctx_t ctx)
void fq_one(fq_t rop, const fq_ctx_t ctx)

Set `rop` to 0, respectively 1.
:::

::: function
void fq_gen(fq_t rop, const fq_ctx_t ctx)

Sets `rop` to the generator of the field, the class of $x$, a root of
the defining polynomial. For $d = 1$ this is the root $-c_0/c_1$ of
$f = c_1 x + c_0$. The generator generates the multiplicative group if
and only if $f$ is primitive, which holds for Conway polynomials.
:::

::: function
void fq_set_ui(fq_t rop, const ulong x, const fq_ctx_t ctx)
void fq_set_si(fq_t rop, const slong x, const fq_ctx_t ctx)
void fq_set_fmpz(fq_t rop, const fmpz_t x, const fq_ctx_t ctx)

Set `rop` to the image of the integer $x$ in $\mathbb{F}_p \subseteq
\mathbb{F}_q$.
:::

::: function
int fq_get_fmpz(fmpz_t a, const fq_t b, const fq_ctx_t ctx)

If `b` lies in $\mathbb{F}_p$, sets `a` to its representative in
$[0, p)$ and returns 1; otherwise returns 0 and leaves `a` unchanged.
:::

::: function
void fq_get_fmpz_poly(fmpz_poly_t a, const fq_t b, const fq_ctx_t ctx)
void fq_get_fmpz_mod_poly(fmpz_mod_poly_t a, const fq_t b, const fq_ctx_t ctx)

Set `a` to the representative of `b` of degree less than $d$, as an
integer polynomial with coefficients in $[0, p)$, respectively as a
polynomial over $\mathbb{F}_p$ (the `fmpz_mod_poly` must have been
initialised with modulus $p$).
:::

::: function
void fq_set_fmpz_poly(fq_t a, const fmpz_poly_t b, const fq_ctx_t ctx)
void fq_set_fmpz_mod_poly(fq_t a, const fmpz_mod_poly_t b, const fq_ctx_t ctx)

Set `a` to the class of the polynomial `b` in $\mathbb{F}_p[x]/(f)$,
reducing the coefficients modulo $p$ and the polynomial modulo $f$. `b`
may have any length.
:::

::: function
void fq_get_fmpz_mod_mat(fmpz_mod_mat_t col, const fq_t a, const fq_ctx_t ctx)
void fq_set_fmpz_mod_mat(fq_t a, const fmpz_mod_mat_t col, const fq_ctx_t ctx)

Convert between `a` and the column vector of its coefficients
$(a_0, \ldots, a_{d-1})^T$ with respect to the basis
$1, x, \ldots, x^{d-1}$. `col` must be a $d \times 1$ matrix over
$\mathbb{F}_p$ with reduced entries.
:::

## Comparison {#comparison}

::: function
int fq_equal(const fq_t op1, const fq_t op2, const fq_ctx_t ctx)

Returns whether `op1` and `op2` are equal.
:::

::: function
int fq_is_zero(const fq_t op, const fq_ctx_t ctx)
int fq_is_one(const fq_t op, const fq_ctx_t ctx)

Return whether `op` is 0, respectively 1.
:::

::: function
int fq_is_invertible(const fq_t op, const fq_ctx_t ctx)

Returns whether `op` is invertible, that is, nonzero.
:::

::: function
int fq_is_invertible_f(fq_t f, const fq_t op, const fq_ctx_t ctx)

Returns whether `op` is invertible modulo the defining polynomial, using
`fq_gcdinv`. If it is not, `f` is set to a nontrivial factor of the
modulus (the monic GCD of the representative of `op` and the modulus);
this can only happen if the modulus is reducible. Otherwise `f` is set to
one. If `op` is zero, the function returns 0 and sets `f` to zero.
:::

## Arithmetic {#arithmetic}

::: function
void fq_add(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
void fq_sub(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)

Set `rop` to $\mathit{op1} + \mathit{op2}$, respectively
$\mathit{op1} - \mathit{op2}$.
:::

::: function
void fq_sub_one(fq_t rop, const fq_t op1, const fq_ctx_t ctx)

Sets `rop` to $\mathit{op1} - 1$.
:::

::: function
void fq_neg(fq_t rop, const fq_t op1, const fq_ctx_t ctx)

Sets `rop` to $-\mathit{op1}$.
:::

::: function
void fq_mul(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)
void fq_sqr(fq_t rop, const fq_t op, const fq_ctx_t ctx)

Set `rop` to $\mathit{op1} \cdot \mathit{op2}$, respectively
$\mathit{op}^2$. The product of the representatives is computed with
`fmpz_poly_mul` (or `fmpz_poly_sqr`) and reduced with `fq_reduce`.
:::

::: function
void fq_mul_fmpz(fq_t rop, const fq_t op, const fmpz_t x, const fq_ctx_t ctx)
void fq_mul_si(fq_t rop, const fq_t op, slong x, const fq_ctx_t ctx)
void fq_mul_ui(fq_t rop, const fq_t op, ulong x, const fq_ctx_t ctx)

Set `rop` to `op` times the image of the integer $x$ in $\mathbb{F}_p$.
:::

::: function
void fq_inv(fq_t rop, const fq_t op1, const fq_ctx_t ctx)

Sets `rop` to $\mathit{op1}^{-1}$, computed with the extended Euclidean
algorithm in $\mathbb{F}_p[x]$. Throws
`FLINT_ERROR` if `op1` is zero.
:::

::: function
void fq_gcdinv(fq_t rop, fq_t inv, const fq_t op, const fq_ctx_t ctx)

Computes the monic GCD $g$ of the representative of `op` and the
defining polynomial $f$, and sets `rop` to $g$ and `inv` to an $s$ with
$s \cdot \mathit{op} \equiv g \pmod f$. In a field $g = 1$ for nonzero
`op`, and `inv` is the inverse; a nontrivial $g$ is a factor of a
reducible modulus. If `op` is zero, `rop` is set to zero and `inv` is not
changed.
:::

::: function
void fq_div(fq_t rop, const fq_t op1, const fq_t op2, const fq_ctx_t ctx)

Sets `rop` to $\mathit{op1} / \mathit{op2}$, computed as
$\mathit{op1} \cdot \mathit{op2}^{-1}$. Throws `FLINT_ERROR` if `op2` is
zero.
:::

::: function
void fq_pow(fq_t rop, const fq_t op1, const fmpz_t e, const fq_ctx_t ctx)
void fq_pow_ui(fq_t rop, const fq_t op, const ulong e, const fq_ctx_t ctx)

Set `rop` to $\mathit{op}^e$ for $e \ge 0$, with $\mathit{op}^0 = 1$ for
every `op`. `fq_pow` throws `FLINT_ERROR` if $e < 0$. If $e \ge p$, the
exponent is first reduced modulo $q - 1$.

::: details
The power is computed by left-to-right binary exponentiation with
`_fq_pow`, which multiplies by the original (short) representative of
`op` and reduces after every step, so the cost is about
$\log_2 \min(e, q)$ squarings and as many multiplications in
$\mathbb{F}_q$.
:::
:::

::: function
void _fq_pow(fmpz * rop, const fmpz * op, slong len, const fmpz_t e, const fq_ctx_t ctx)

Sets `(rop, 2d - 1)` to `(op, len)` raised to the power $e \ge 0$,
reduced modulo $f$: the first $d$ coefficients hold the reduced result and
the others are zero. Requires $1 \le \mathit{len} \le d$ and space for
$2d - 1$ coefficients in `rop`. Does not allow aliasing. For $e = 0$ the
result is 1.
:::

## Roots {#roots}

::: function
int fq_sqrt(fq_t rop, const fq_t op, const fq_ctx_t ctx)

If `op` is a square, sets `rop` to a square root of it and returns 1;
otherwise returns 0, and `rop` is unspecified.

::: details
For $p = 2$ the square root is the inverse of the Frobenius map,
computed by `fq_pth_root`. For odd $p$ the function uses the
Tonelli–Shanks algorithm [@Ton1891] [@Sha1973] in the cyclic group
$\mathbb{F}_q^\times$ of order $q - 1 = 2^s t$, $t$ odd: it needs a
quadratic non-residue $z$, which is the generator if the context was
created from a Conway polynomial (a primitive element is never a square)
and is found by random sampling otherwise, and then $O(s^2)$
multiplications after three exponentiations.
:::
:::

::: function
int fq_is_square(const fq_t op, const fq_ctx_t ctx)

Returns whether `op` is a square, by Euler's criterion
$\mathit{op}^{(q-1)/2} = 1$. If $p = 2$, every element is a square.
:::

::: function
void fq_pth_root(fq_t rop, const fq_t op1, const fq_ctx_t ctx)

Sets `rop` to the unique $p$-th root of `op1`, which is
$\mathit{op1}^{p^{d-1}}$, computed by $d - 1$ successive $p$-th powers.
:::

## Frobenius map, norm, trace and order {#frobenius}

The Galois group of $\mathbb{F}_q/\mathbb{F}_p$ is cyclic of order $d$,
generated by the Frobenius automorphism $\sigma(a) = a^p$. The norm and
the trace of $a \in \mathbb{F}_q$ are
$N(a) = \prod_{i=0}^{d-1} \sigma^i(a)$ and
$\operatorname{Tr}(a) = \sum_{i=0}^{d-1} \sigma^i(a)$; they are the
determinant and the trace of the $\mathbb{F}_p$-linear map
$x \mapsto ax$ and lie in $\mathbb{F}_p$.

::: summary algorithms/finite-fields/frobenius-norm-trace :::

::: function
void fq_frobenius(fq_t rop, const fq_t op, slong e, const fq_ctx_t ctx)

Sets `rop` to $\sigma^e(\mathit{op}) = \mathit{op}^{p^e}$. Any $e$ is
allowed, including negative values; it is reduced modulo $d$. The power
is computed by binary exponentiation, at the cost of about
$e \log_2 p$ squarings for the reduced $e$.
:::

::: function
void _fq_frobenius(fmpz * rop, const fmpz * op, slong len, slong e, const fq_ctx_t ctx)

Sets `(rop, 2d - 1)` to the image of `(op, len)` under $\sigma^e$, with
the reduced result in the first $d$ coefficients. Requires
$1 \le \mathit{len} \le d$, $1 \le e < d$, space for $2d - 1$
coefficients in `rop`, and no aliasing.
:::

::: function
void fq_norm(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
void _fq_norm(fmpz_t rop, const fmpz * op, slong len, const fq_ctx_t ctx)

Set `rop` to the norm of `op`, respectively of the nonzero element
`(op, len)` with $1 \le \mathit{len} \le d$, as an integer in $[0, p)$.
The norm of an element with representative $a(x)$ is computed as the
resultant $\operatorname{Res}(f, a) / \operatorname{lc}(f)^{\deg a}$,
or as $a_0^d$ for an element of $\mathbb{F}_p$.
:::

::: function
void fq_trace(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
void _fq_trace(fmpz_t rop, const fmpz * op, slong len, const fq_ctx_t ctx)

Set `rop` to the trace of `op`, respectively of the nonzero element
`(op, len)` with $1 \le \mathit{len} \le d$, as an integer in $[0, p)$.
The traces $\operatorname{Tr}(x^i)$ of the basis elements are computed
from the coefficients of $f$ by Newton's identities, and the trace of
$a = \sum a_i x^i$ is $\sum a_i \operatorname{Tr}(x^i)$.
:::

::: function
int fq_multiplicative_order(fmpz * ord, const fq_t op, const fq_ctx_t ctx)

Computes the order of `op` in the multiplicative group, which divides
$q - 1$, and stores it in `ord` unless `ord` is `NULL`. Returns 0 if `op`
is zero (and sets the order to 0), 1 if `op` is a primitive element
(its order is $q - 1$), and $-1$ otherwise.

::: details
The function factors $q - 1$ with `fmpz_factor`, which can take a long
time for large $q$ when $q - 1$ has large prime factors, and then
removes prime factors from $q - 1$ as long as the power of `op` stays 1.
This also tests whether a defining polynomial is primitive: its root
`fq_gen` is primitive if and only if this function returns 1.
:::
:::

::: function
int fq_is_primitive(const fq_t op, const fq_ctx_t ctx)

Returns whether `op` generates the multiplicative group
$\mathbb{F}_q^\times$, using `fq_multiplicative_order`.
:::

## Random elements {#random}

::: function
void fq_rand(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
void fq_rand_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)

Set `rop` to a uniformly random element, respectively a uniformly random
nonzero element.
:::

::: function
void fq_randtest(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)
void fq_randtest_not_zero(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)

Set `rop` to a random element, respectively a random nonzero element, for
testing: with probability 1/2 the coefficients are uniform, otherwise
most of them are zero.
:::

::: function
void fq_randtest_dense(fq_t rop, flint_rand_t state, const fq_ctx_t ctx)

Sets `rop` to a random element whose representative has degree
exactly $d - 1$ and leading coefficient 1, for testing.
:::

## Input and output {#io}

::: function
int fq_fprint(FILE * file, const fq_t op, const fq_ctx_t ctx)
void fq_print(const fq_t op, const fq_ctx_t ctx)

Print the representative of `op` to `file`, respectively `stdout`, in the
format of `fmpz_poly_print`: the length, two spaces, and the
coefficients from the constant term up, separated by spaces.
`fq_fprint` returns a positive value on success.
:::

::: function
int fq_fprint_pretty(FILE * file, const fq_t op, const fq_ctx_t ctx)
int fq_print_pretty(const fq_t op, const fq_ctx_t ctx)

Print the representative of `op` as a polynomial in the variable name of
the context, such as `3*a^2+1`, to `file`, respectively `stdout`. Return a
positive value on success.
:::

::: function
char * fq_get_str(const fq_t op, const fq_ctx_t ctx)
char * fq_get_str_pretty(const fq_t op, const fq_ctx_t ctx)

Return the plain, respectively readable, representation of `op` as a
newly allocated string, which the caller must free with `flint_free`.
:::

## Bit packing {#bit-packing}

These functions are used by the Kronecker substitution algorithms of
`fq_poly` and `fq_mat`.

::: function
void fq_bit_pack(fmpz_t f, const fq_t op, flint_bitcnt_t bit_size, const fq_ctx_t ctx)

Sets `f` to $\sum_i a_i 2^{i \cdot \mathit{bit\_size}}$, where the $a_i$
are the coefficients of `op`. Each coefficient must fit in `bit_size`
bits, that is, $p \le 2^{\mathit{bit\_size}}$.
:::

::: function
void fq_bit_unpack(fq_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_ctx_t ctx)

Sets `rop` to the class of the polynomial whose coefficients are the
consecutive `bit_size`-bit fields of the nonnegative integer `f`, reducing
the coefficients modulo $p$ and the polynomial modulo $f$. This is the
inverse of `fq_bit_pack` and recovers a product from the packed product
of two packed operands if the fields are wide enough.
:::

## Reduction modulo the defining polynomial {#reduction}

::: summary algorithms/finite-fields/reduction :::

::: function
void fq_reduce(fq_t rop, const fq_ctx_t ctx)

Reduces the integer polynomial `rop`, with arbitrary integer coefficients,
modulo $p$ and modulo the defining polynomial, in place, so that it
becomes a valid element. The length of `rop` must be at most $2d$ unless
the defining polynomial is sparse (see below); products of two elements
satisfy this.

::: warning
For a defining polynomial with more than five nonzero terms, the
reduction is a Newton division that is only valid for inputs of length
at most $2d$; longer inputs give wrong results without any error. Use
`fq_set_fmpz_poly`, which handles any length, to reduce long
polynomials.
:::
:::

::: function
void _fq_reduce(fmpz * R, slong lenR, const fq_ctx_t ctx)
void _fq_sparse_reduce(fmpz * R, slong lenR, const fq_ctx_t ctx)
void _fq_dense_reduce(fmpz * R, slong lenR, const fq_ctx_t ctx)

Reduce the polynomial `(R, lenR)`, with arbitrary integer coefficients and
$\mathit{lenR} \le 2d$ (any length for `_fq_sparse_reduce`),
modulo $p$ and modulo the defining polynomial $f$ in place: the first
$\min(\mathit{lenR}, d)$ coefficients receive the reduced value (and are
not normalised). `_fq_sparse_reduce` eliminates the coefficients of
degree at least $d$ from the top, using the nonzero terms of the monic
version of $f$, and reduces modulo $p$ at the end; it costs about
$(\mathit{lenR} - d) \cdot t$ multiplications for $f$ with $t$ nonzero
terms. `_fq_dense_reduce` reduces modulo $p$ and then divides by $f$
with the precomputed inverse (Newton division,
`_fmpz_mod_poly_divrem_newton_n_preinv`). `_fq_reduce` calls the first if
$f$ has at most five nonzero terms and the second otherwise.
:::

## Internal functions {#internal}

::: function
void __fq_ctx_prime(fmpz_t p, fq_ctx_t ctx)

Sets `p` to the characteristic. Equivalent to `fmpz_set(p,
fq_ctx_prime(ctx))`; kept for compatibility.
:::
