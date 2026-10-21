---
title: fq_nmod.h
header: fq_nmod.h
template-headers: fq_templates
description: Finite fields $\mathbb{F}_{p^d}$ of word-size characteristic, with elements stored as polynomials over $\Z/p\Z$ (`nmod_poly`).
category: Finite fields
order: 11
---

This module implements arithmetic in a finite field
$\mathbb{F}_q = \mathbb{F}_p[x]/(f)$, $q = p^d$, for a prime
$p < 2^{\bits}$. An element is stored as its unique representative of
degree less than $d$, an `nmod_poly_t` with coefficients in $[0, p)$, and
a product is an `nmod_poly` product followed by reduction modulo $f$
([](algorithms/finite-fields/reduction)). It is the representation chosen
by `fq_default` for extension fields of word-size characteristic that are
too large for Zech logarithms. The modules `fq_nmod_vec`, `fq_nmod_mat`,
`fq_nmod_poly`, `fq_nmod_poly_factor` and `fq_nmod_embed` provide
vectors, matrices, polynomials, factorisation and embeddings over these
fields.

The interface is the one documented in detail for [](fq), with `fq_nmod`
in place of `fq`, word-size integers (`ulong`, `nmod_poly`, `nmod_mat`) in
place of multiprecision ones (`fmpz`, `fmpz_mod_poly`, `fmpz_mod_mat`) for
the characteristic and the conversions, and a few additional context
functions.

::: summary algorithms/finite-fields/overview :::

## Types and conventions {#types}

::: type
fq_nmod_ctx_struct
fq_nmod_ctx_t

A finite field context: the `nmod_t` of $p$ (`mod`), the defining
polynomial (`modulus`) and the power series inverse of its reversal
(`inv`), both `nmod_poly_t`, the nonzero coefficients and exponents of
the monic multiple of $f$ (`a`, `j`, `len`), the flags `sparse_modulus`
(at most five nonzero terms) and `is_conway`, and the name of the
generator (`var`).
:::

::: type
fq_nmod_struct
fq_nmod_t

An element of $\mathbb{F}_q$: an `nmod_poly_struct` holding the reduced
representative, normalised, with the modulus $p$ stored in it.
:::

The conventions are those of [](fq): elements must be reduced, outputs
may alias inputs, the context is the last argument of every function,
and the defining polynomial must be irreducible and $p$ prime (not
checked).

## Context management {#context}

::: function
void fq_nmod_ctx_init_ui(fq_nmod_ctx_t ctx, ulong prime, slong deg, const char * var)

Initialises `ctx` as a context for $\mathbb{F}_{p^d}$, $p$ = `prime`,
$d$ = `deg` $\ge 1$, with the generator printed as `var`. The defining
polynomial is the Conway polynomial if $(p, d)$ is in FLINT's table, and
otherwise a polynomial with the least possible number of nonzero terms
(`nmod_poly_minimal_irreducible`).
:::

::: function
void fq_nmod_ctx_init_conway_ui(fq_nmod_ctx_t ctx, ulong prime, slong deg, const char * var)
int _fq_nmod_ctx_init_conway_ui(fq_nmod_ctx_t ctx, ulong prime, slong deg, const char * var)

Initialise `ctx` with the Conway polynomial for $(p, d)$. The first
function throws `FLINT_ERROR` if the polynomial is not in FLINT's table;
the second returns 0 in that case (leaving `ctx` uninitialised) and 1 on
success, and does not set the flag `is_conway`.
:::

::: function
void fq_nmod_ctx_init_minimal_weight_ui(fq_nmod_ctx_t ctx, ulong prime, slong deg, const char * var)

Initialises `ctx` with a monic irreducible polynomial of degree $d$ with
the least possible number of nonzero terms, as computed by
`nmod_poly_minimal_irreducible`: a binomial $x^d + a$ if one exists,
otherwise a trinomial, tetranomial or pentanomial, the smallest in a
fixed order. The polynomial is deterministic but in general not
primitive.
:::

::: function
void fq_nmod_ctx_init_modulus(fq_nmod_ctx_t ctx, const nmod_poly_t modulus, const char * var)

Initialises `ctx` as a context for $\mathbb{F}_p[x]/(f)$, where $f$ is
`modulus`, an irreducible polynomial of degree $d \ge 1$ over
$\mathbb{F}_p$ and $p$ is its modulus. $f$ need not be monic and is
copied. Irreducibility is not checked.
:::

::: function
void fq_nmod_ctx_init(fq_nmod_ctx_t ctx, fmpz_t p, slong d, const char * var)
void fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, fmpz_t p, slong d, const char * var)
int _fq_nmod_ctx_init_conway(fq_nmod_ctx_t ctx, fmpz_t p, slong d, const char * var)

Deprecated versions of `fq_nmod_ctx_init_ui`, `fq_nmod_ctx_init_conway_ui`
and `_fq_nmod_ctx_init_conway_ui` taking the prime as an `fmpz_t`, which
must fit in a word.
:::

::: function
void fq_nmod_ctx_init_randtest(fq_nmod_ctx_t ctx, flint_rand_t state, int type)

Initialises `ctx` as a random field for testing, defined by a Conway
polynomial, multiplied by a random nonzero constant with probability
1/2 (to test non-monic moduli). The prime and degree are chosen by
`_nmod_poly_conway_rand` according to `type` as for `fq_ctx_init_randtest`.
:::

::: function
void fq_nmod_ctx_init_randtest_reducible(fq_nmod_ctx_t ctx, flint_rand_t state, int type)

Initialises `ctx` with a random monic modulus that is usually reducible,
for testing: primes of up to 20 bits for `type` 0 and 1, of up to 10 bits
for 2 and 3, degrees up to 30 for even and up to 15 for odd `type`.
:::

::: function
void fq_nmod_ctx_clear(fq_nmod_ctx_t ctx)

Frees the memory used by `ctx`.
:::

::: function
const nmod_poly_struct * fq_nmod_ctx_modulus(const fq_nmod_ctx_t ctx)

Returns a pointer to the defining polynomial.
:::

::: function
slong fq_nmod_ctx_degree(const fq_nmod_ctx_t ctx)

Returns the degree $d$.
:::

::: function
ulong fq_nmod_ctx_prime(const fq_nmod_ctx_t ctx)

Returns the characteristic $p$.
:::

::: function
void fq_nmod_ctx_order(fmpz_t f, const fq_nmod_ctx_t ctx)

Sets `f` to $q = p^d$.
:::

::: function
int fq_nmod_ctx_fprint(FILE * file, const fq_nmod_ctx_t ctx)
void fq_nmod_ctx_print(const fq_nmod_ctx_t ctx)

Print $p$, $d$ and the defining polynomial to `file`, respectively
`stdout`, as `fq_ctx_fprint`.
:::

## Memory management {#memory}

::: function
void fq_nmod_init(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
void fq_nmod_init2(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

Initialise `rop` as zero; `fq_nmod_init2` allocates space for $d$
coefficients. The context is needed for the modulus $p$.
:::

::: function
void fq_nmod_clear(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

Frees the memory used by `rop`.
:::

## Assignment and conversion {#assignment}

::: function
void fq_nmod_set(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
void fq_nmod_swap(fq_nmod_t op1, fq_nmod_t op2, const fq_nmod_ctx_t ctx)
void fq_nmod_zero(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
void fq_nmod_one(fq_nmod_t rop, const fq_nmod_ctx_t ctx)
void fq_nmod_gen(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

Copy, swap, and set to 0, 1 or the generator (the class of $x$, or the
root of the defining polynomial if $d = 1$), as `fq_set`, `fq_swap`,
`fq_zero`, `fq_one` and `fq_gen`.
:::

::: function
void fq_nmod_set_ui(fq_nmod_t rop, const ulong x, const fq_nmod_ctx_t ctx)
void fq_nmod_set_si(fq_nmod_t rop, const slong x, const fq_nmod_ctx_t ctx)
void fq_nmod_set_fmpz(fq_nmod_t rop, const fmpz_t x, const fq_nmod_ctx_t ctx)

Set `rop` to the image of the integer $x$ in $\mathbb{F}_p$.
:::

::: function
int fq_nmod_get_fmpz(fmpz_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx)

If `b` lies in $\mathbb{F}_p$, sets `a` to its representative in
$[0, p)$ and returns 1; otherwise returns 0.
:::

::: function
void fq_nmod_get_nmod_poly(nmod_poly_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx)
void fq_nmod_set_nmod_poly(fq_nmod_t a, const nmod_poly_t b, const fq_nmod_ctx_t ctx)

Convert between an element and its representative in $\mathbb{F}_p[x]$.
`fq_nmod_get_nmod_poly` also sets the modulus of `a` to $p$.
`fq_nmod_set_nmod_poly` reduces a polynomial of any length modulo the
defining polynomial; its modulus must be $p$.
:::

::: function
void fq_nmod_get_nmod_mat(nmod_mat_t col, const fq_nmod_t a, const fq_nmod_ctx_t ctx)
void fq_nmod_set_nmod_mat(fq_nmod_t a, const nmod_mat_t col, const fq_nmod_ctx_t ctx)

Convert between `a` and the $d \times 1$ column vector of its
coefficients, as `fq_get_fmpz_mod_mat` and `fq_set_fmpz_mod_mat`.
:::

## Comparison {#comparison}

::: function
int fq_nmod_equal(const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
int fq_nmod_is_zero(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
int fq_nmod_is_one(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

Return whether `op1` equals `op2`, whether `op` is 0, respectively 1.
:::

::: function
int fq_nmod_cmp(const fq_nmod_t a, const fq_nmod_t b, const fq_nmod_ctx_t ctx)

Compares the representatives of $a$ and $b$ and returns $-1$, 0 or 1: a
shorter representative is smaller, and representatives of the same
length are compared coefficientwise from the constant term up. This is
a total order for sorting, unrelated to the field structure.
:::

::: function
int fq_nmod_is_invertible(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
int fq_nmod_is_invertible_f(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

Return whether `op` is invertible; the second version also sets `rop` to
a nontrivial factor of a reducible modulus, as `fq_is_invertible_f`.
:::

## Arithmetic {#arithmetic}

::: function
void fq_nmod_add(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
void fq_nmod_sub(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
void fq_nmod_sub_one(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
void fq_nmod_neg(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op1} + \mathit{op2}$, $\mathit{op1} - \mathit{op2}$,
$\mathit{op1} - 1$, respectively $-\mathit{op1}$.
:::

::: function
void fq_nmod_mul(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)
void fq_nmod_sqr(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op1} \cdot \mathit{op2}$, respectively
$\mathit{op}^2$, with `nmod_poly_mul` followed by `fq_nmod_reduce`.
:::

::: function
void fq_nmod_mul_fmpz(fq_nmod_t rop, const fq_nmod_t op, const fmpz_t x, const fq_nmod_ctx_t ctx)
void fq_nmod_mul_si(fq_nmod_t rop, const fq_nmod_t op, slong x, const fq_nmod_ctx_t ctx)
void fq_nmod_mul_ui(fq_nmod_t rop, const fq_nmod_t op, ulong x, const fq_nmod_ctx_t ctx)

Set `rop` to `op` times the image of the integer $x$ in $\mathbb{F}_p$.
:::

::: function
void fq_nmod_inv(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)
void fq_nmod_div(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_t op2, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op1}^{-1}$, respectively $\mathit{op1} /
\mathit{op2}$. Throw `FLINT_ERROR` on division by zero.
:::

::: function
void fq_nmod_gcdinv(fq_nmod_t rop, fq_nmod_t inv, const fq_nmod_t op, const fq_nmod_ctx_t ctx)

Sets `rop` to the monic GCD of the representative of `op` and the
modulus and `inv` to the corresponding inverse, as `fq_gcdinv`.
:::

::: function
void fq_nmod_pow(fq_nmod_t rop, const fq_nmod_t op1, const fmpz_t e, const fq_nmod_ctx_t ctx)
void fq_nmod_pow_ui(fq_nmod_t rop, const fq_nmod_t op1, const ulong e, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op1}^e$ for $e \ge 0$ (with $\mathit{op1}^0 = 1$);
`fq_nmod_pow` throws `FLINT_ERROR` for $e < 0$ and reduces $e$ modulo
$q - 1$ when $e \ge p$.
:::

::: function
void _fq_nmod_pow(ulong * rop, const ulong * op, slong len, const fmpz_t e, const fq_nmod_ctx_t ctx)

Sets `(rop, 2d - 1)` to `(op, len)` raised to the power $e \ge 0$, reduced
(the result in the first $d$ coefficients), as `_fq_pow`. Requires
$1 \le \mathit{len} \le d$ and no aliasing.
:::

## Roots {#roots}

::: function
int fq_nmod_sqrt(fq_nmod_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
int fq_nmod_is_square(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
void fq_nmod_pth_root(fq_nmod_t rop, const fq_nmod_t op1, const fq_nmod_ctx_t ctx)

Square root (returning 0 if `op` is not a square), squareness test and
$p$-th root, as `fq_sqrt`, `fq_is_square` and `fq_pth_root`.
:::

## Frobenius map, norm, trace and order {#frobenius}

::: function
void fq_nmod_frobenius(fq_nmod_t rop, const fq_nmod_t op, slong e, const fq_nmod_ctx_t ctx)
void _fq_nmod_frobenius(ulong * rop, const ulong * op, slong len, slong e, const fq_nmod_ctx_t ctx)

Set `rop` to $\mathit{op}^{p^e}$, as `fq_frobenius` and `_fq_frobenius`
(the underscore version requires $1 \le \mathit{len} \le d$,
$1 \le e < d$ and space for $2d - 1$ coefficients).
:::

::: function
void fq_nmod_norm(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
void _fq_nmod_norm(fmpz_t rop, const ulong * op, slong len, const fq_nmod_ctx_t ctx)
void fq_nmod_trace(fmpz_t rop, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
void _fq_nmod_trace(fmpz_t rop, const ulong * op, slong len, const fq_nmod_ctx_t ctx)

Set `rop` to the norm, respectively the trace, of the element, an
integer in $[0, p)$, computed by a resultant, respectively from the
traces of the powers of the generator, as `fq_norm` and `fq_trace`. The
underscore versions require a nonzero element with
$1 \le \mathit{len} \le d$.
:::

::: function
int fq_nmod_multiplicative_order(fmpz * ord, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
int fq_nmod_is_primitive(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

Multiplicative order (returning 0 for zero, 1 for a primitive element
and $-1$ otherwise) and primitivity test, as `fq_multiplicative_order`
and `fq_is_primitive`.
:::

## Random elements {#random}

::: function
void fq_nmod_rand(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)
void fq_nmod_rand_not_zero(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)
void fq_nmod_randtest(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)
void fq_nmod_randtest_not_zero(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)
void fq_nmod_randtest_dense(fq_nmod_t rop, flint_rand_t state, const fq_nmod_ctx_t ctx)

Uniformly random elements (`rand`), random elements for testing
(`randtest`, often sparse) and random elements with a monic
representative of degree $d - 1$ (`randtest_dense`), as for `fq`.
:::

## Input and output {#io}

::: function
int fq_nmod_fprint(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
void fq_nmod_print(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

Print the representative of `op` in the format of `nmod_poly_fprint`:
the length, the modulus $p$, and the coefficients from the constant term
up.
:::

::: function
int fq_nmod_fprint_pretty(FILE * file, const fq_nmod_t op, const fq_nmod_ctx_t ctx)
void fq_nmod_print_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

Print the representative as a polynomial in the variable name of the
context, such as `3*a^2+1`.
:::

::: function
char * fq_nmod_get_str(const fq_nmod_t op, const fq_nmod_ctx_t ctx)
char * fq_nmod_get_str_pretty(const fq_nmod_t op, const fq_nmod_ctx_t ctx)

Return the plain, respectively readable, representation as a newly
allocated string, to be freed with `flint_free`.
:::

## Bit packing {#bit-packing}

::: function
void fq_nmod_bit_pack(fmpz_t f, const fq_nmod_t op, flint_bitcnt_t bit_size, const fq_nmod_ctx_t ctx)
void fq_nmod_bit_unpack(fq_nmod_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_nmod_ctx_t ctx)

Pack the coefficients of `op` into `bit_size`-bit fields of an integer,
and unpack and reduce, as `fq_bit_pack` and `fq_bit_unpack`.
:::

## Reduction modulo the defining polynomial {#reduction}

::: function
void fq_nmod_reduce(fq_nmod_t rop, const fq_nmod_ctx_t ctx)

Reduces the polynomial `rop`, of length at most $2d$ with reduced
coefficients, modulo the defining polynomial, in place. Longer inputs
are only allowed for a sparse defining polynomial (they trigger an
assertion failure in debug builds and give wrong results otherwise);
`fq_nmod_set_nmod_poly` reduces polynomials of any length.
:::

::: function
void _fq_nmod_reduce(ulong * R, slong lenR, const fq_nmod_ctx_t ctx)
void _fq_nmod_sparse_reduce(ulong * R, slong lenR, const fq_nmod_ctx_t ctx)
void _fq_nmod_dense_reduce(ulong * R, slong lenR, const fq_nmod_ctx_t ctx)

Reduce `(R, lenR)`, with coefficients in $[0, p)$ and $\mathit{lenR} \le
2d$ (any length for the sparse version), modulo the defining polynomial in
place, leaving the result in the first $d$ coefficients.
`_fq_nmod_sparse_reduce` eliminates the high coefficients with the
nonzero terms of $f$, with special code for $p = 2$ (exclusive or), for
moduli whose nonzero coefficients are all 1 (subtractions only), and for
$p < 2^{\bits/2}$ (Barrett reduction of a single-word product);
`_fq_nmod_dense_reduce` uses Newton division with the precomputed
inverse. `_fq_nmod_reduce` chooses the sparse version if $f$ has at most
five nonzero terms. See [](algorithms/finite-fields/reduction).
:::
