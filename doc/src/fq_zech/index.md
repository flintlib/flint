---
title: fq_zech.h
header: fq_zech.h
template-headers: fq_templates
description: Small finite fields $\mathbb{F}_q$, $q < 2^{64}$, with elements stored as discrete logarithms (Zech logarithm representation).
category: Finite fields
order: 17
---

This module implements arithmetic in a finite field $\mathbb{F}_q$,
$q = p^d < 2^{\bits}$, in the Zech logarithm representation: a nonzero
element is stored as the integer $k \in [0, q - 2]$ with $a = g^k$, where
$g$ is the generator of the field, a root of a *primitive* defining
polynomial, and zero is stored as $q - 1$. Multiplication, division and
powering are additions, subtractions and multiplications of logarithms
modulo $q - 1$, and addition uses a table of Zech logarithms
$Z(k)$, defined by $1 + g^k = g^{Z(k)}$, with $q$ entries. The context
therefore holds tables of $2q + p$ words, built in $O(q)$ operations when
it is created, which limits the representation to small fields; `fq_default`
uses it for $q < 2^{16}$.

The interface is the one documented in detail for [](fq), with `fq_zech`
in place of `fq` and word-size integers (`ulong`, `nmod_poly`, `nmod_mat`)
in place of multiprecision ones, plus conversions to and from the
`fq_nmod` representation, on which a Zech context is built. Multivariate polynomials
with coefficients in this representation are in [](fq_zech_mpoly).

::: summary algorithms/finite-fields/zech-logarithms :::

## Types {#types}

::: type
fq_zech_ctx_struct
fq_zech_ctx_t

A Zech logarithm context: $q - 1$ (`qm1`), $(q - 1)/2$ (`qm1o2`, 0 if
$p = 2$; the logarithm of $-1$), $(q - 1)/(p - 1)$ (`qm1opm1`), $p$
with a floating-point inverse, the norm of the generator (`prime_root`,
a primitive root modulo $p$), the tables `zech_log_table` (the Zech
logarithms $Z(k)$), `prime_field_table` (the logarithms of
$1, \ldots, p - 1$) and `eval_table` (for each $k$, the integer
$\sum_i a_i p^i$ encoding the coefficients $a_i$ of $g^k$), a pointer to
the underlying `fq_nmod` context, a flag `owns_fq_nmod_ctx` and the flag
`is_conway`.
:::

::: type
fq_zech_struct
fq_zech_t

An element: a structure with one field `value`, the discrete logarithm
of the element, or $q - 1$ for zero. Elements need no memory management,
but `fq_zech_init` and `fq_zech_clear` exist for uniformity.
:::

## Context management {#context}

::: function
void fq_zech_ctx_init_ui(fq_zech_ctx_t ctx, ulong p, slong d, const char * var)

Initialises `ctx` as a Zech context for $\mathbb{F}_{p^d}$, with the
generator printed as `var`, defined by the Conway polynomial if it is in
FLINT's table and by a random primitive polynomial otherwise. Requires
$p^d < 2^{\bits}$ (otherwise `FLINT_ERROR` is thrown); the tables need
$2q + p$ words.
:::

::: function
void fq_zech_ctx_init_conway_ui(fq_zech_ctx_t ctx, ulong p, slong d, const char * var)
int _fq_zech_ctx_init_conway_ui(fq_zech_ctx_t ctx, ulong p, slong d, const char * var)

Initialise `ctx` with the Conway polynomial for $(p, d)$, which is
primitive. The first function throws `FLINT_ERROR` if the polynomial is
not in FLINT's table; the second returns 0 in that case and 1 on success.
:::

::: function
void fq_zech_ctx_init_random_ui(fq_zech_ctx_t ctx, ulong p, slong d, const char * var)

Initialises `ctx` with a random monic primitive polynomial of degree $d$
(from `nmod_poly_randtest_monic_primitive`, with an internal random
state).
:::

::: function
void fq_zech_ctx_init_modulus(fq_zech_ctx_t ctx, const nmod_poly_t modulus, const char * var)
int fq_zech_ctx_init_modulus_check(fq_zech_ctx_t ctx, const nmod_poly_t modulus, const char * var)

Initialise `ctx` with the given defining polynomial, which must be
primitive; it is copied into a new `fq_nmod` context owned by `ctx`.
While the tables are built, primitivity is checked: the powers
$g^0, \ldots, g^{q-2}$ must be distinct. `fq_zech_ctx_init_modulus`
throws `FLINT_ERROR` for a non-primitive polynomial;
`fq_zech_ctx_init_modulus_check` returns 1 on success and 0 if the
polynomial is not primitive. After a return value of 0, all memory has
already been freed: `ctx` must not be used, and need not be cleared
(clearing it is allowed and has no effect).
`fq_default_ctx_init_modulus_nmod` falls back to `fq_nmod` for a
non-primitive modulus.

::: details
In FLINT 3.6 and earlier,
`fq_zech_ctx_init_modulus` did not throw but returned an unusable
context, and clearing a context after a failed check freed the tables
twice.
:::
:::

::: function
void fq_zech_ctx_init_fq_nmod_ctx(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctxn)
int fq_zech_ctx_init_fq_nmod_ctx_check(fq_zech_ctx_t ctx, fq_nmod_ctx_t ctxn)

Initialise `ctx` as the Zech representation of the field of the
`fq_nmod` context `ctxn`, which is referenced, not copied: it must stay
alive while `ctx` is used and is not cleared with `ctx`. The defining
polynomial must be primitive; `fq_zech_ctx_init_fq_nmod_ctx` throws
`FLINT_ERROR` otherwise, and the `_check` version returns 0 instead of 1
(after which the memory allocated for `ctx` has been freed: `ctx` must
not be used, clearing it is allowed and has no effect, and `ctxn`
remains owned by the caller).
:::

::: function
void fq_zech_ctx_init(fq_zech_ctx_t ctx, fmpz_t p, slong d, const char * var)
void fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, fmpz_t p, slong d, const char * var)
int _fq_zech_ctx_init_conway(fq_zech_ctx_t ctx, fmpz_t p, slong d, const char * var)
void fq_zech_ctx_init_random(fq_zech_ctx_t ctx, fmpz_t p, slong d, const char * var)

Deprecated versions of the `_ui` functions, taking the prime as an
`fmpz_t` that must fit in a word.
:::

::: function
void fq_zech_ctx_init_randtest(fq_zech_ctx_t ctx, flint_rand_t state, int type)
void fq_zech_ctx_init_randtest_reducible(fq_zech_ctx_t ctx, flint_rand_t state, int type)

Initialise `ctx` as a random field for testing, with a random primitive
modulus: primes of up to 5 bits for `type` 0 and 1, of up to 4 bits for
2 and 3, and degree at most 4 for even and 3 for odd `type`. The
reducible version is the same, since a Zech context needs an irreducible
(indeed primitive) modulus.
:::

::: function
void fq_zech_ctx_clear(fq_zech_ctx_t ctx)

Frees the tables, and the `fq_nmod` context if `ctx` owns it.
:::

::: function
const nmod_poly_struct * fq_zech_ctx_modulus(const fq_zech_ctx_t ctx)
slong fq_zech_ctx_degree(const fq_zech_ctx_t ctx)
ulong fq_zech_ctx_prime(const fq_zech_ctx_t ctx)

Return the defining polynomial, the degree $d$ and the characteristic
$p$.
:::

::: function
ulong fq_zech_ctx_order_ui(const fq_zech_ctx_t ctx)
void fq_zech_ctx_order(fmpz_t f, const fq_zech_ctx_t ctx)

Return, respectively set `f` to, the number of elements $q$.
`fq_zech_ctx_order` is deprecated.
:::

::: function
int fq_zech_ctx_fprint(FILE * file, const fq_zech_ctx_t ctx)
void fq_zech_ctx_print(const fq_zech_ctx_t ctx)

Print `Zech Representation:` followed by the description of the
underlying `fq_nmod` context to `file`, respectively `stdout`.
:::

## Memory management {#memory}

::: function
void fq_zech_init(fq_zech_t rop, const fq_zech_ctx_t ctx)
void fq_zech_init2(fq_zech_t rop, const fq_zech_ctx_t ctx)
void fq_zech_clear(fq_zech_t rop, const fq_zech_ctx_t ctx)

`fq_zech_init` and `fq_zech_init2` set `rop` to zero; `fq_zech_clear` does
nothing.
:::

::: function
void fq_zech_reduce(fq_zech_t rop, const fq_zech_ctx_t ctx)

Subtracts $q$ from the value of `rop` if it is at least $q$. Values in
$[0, q - 1]$ are valid elements and are not changed.
:::

## Assignment and conversion {#assignment}

::: function
void fq_zech_set(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
void fq_zech_swap(fq_zech_t op1, fq_zech_t op2, const fq_zech_ctx_t ctx)
void fq_zech_zero(fq_zech_t rop, const fq_zech_ctx_t ctx)
void fq_zech_one(fq_zech_t rop, const fq_zech_ctx_t ctx)
void fq_zech_gen(fq_zech_t rop, const fq_zech_ctx_t ctx)

Copy, swap, and set to 0 (logarithm $q - 1$), 1 (logarithm 0) or the
generator (logarithm 1).
:::

::: function
void fq_zech_set_ui(fq_zech_t rop, const ulong x, const fq_zech_ctx_t ctx)
void fq_zech_set_si(fq_zech_t rop, const slong x, const fq_zech_ctx_t ctx)
void fq_zech_set_fmpz(fq_zech_t rop, const fmpz_t x, const fq_zech_ctx_t ctx)

Set `rop` to the image of the integer $x$ in $\mathbb{F}_p$, using the
table of logarithms of the prime field.
:::

::: function
int fq_zech_get_fmpz(fmpz_t a, const fq_zech_t op, const fq_zech_ctx_t ctx)

If `op` lies in $\mathbb{F}_p$, sets `a` to its representative in
$[0, p)$ and returns 1; otherwise returns 0.
:::

::: function
void fq_zech_get_fq_nmod(fq_nmod_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
void fq_zech_set_fq_nmod(fq_zech_t rop, const fq_nmod_t op, const fq_zech_ctx_t ctx)

Convert between the Zech representation and the polynomial
representation of the underlying `fq_nmod` context `ctx->fq_nmod_ctx`.
`fq_zech_get_fq_nmod` reads the coefficients from `eval_table`;
`fq_zech_set_fq_nmod` sums $c_i g^i$ with $d$ table lookups.
:::

::: function
void fq_zech_get_nmod_poly(nmod_poly_t a, const fq_zech_t b, const fq_zech_ctx_t ctx)
void fq_zech_set_nmod_poly(fq_zech_t a, const nmod_poly_t b, const fq_zech_ctx_t ctx)

Convert between an element and its representative in $\mathbb{F}_p[x]$.
`fq_zech_get_nmod_poly` also sets the modulus of `a` to $p$.
`fq_zech_set_nmod_poly` accepts a polynomial of any length, using
$g^{q-1} = 1$ to reduce exponents.
:::

::: function
void fq_zech_get_nmod_mat(nmod_mat_t col, const fq_zech_t a, const fq_zech_ctx_t ctx)
void fq_zech_set_nmod_mat(fq_zech_t a, const nmod_mat_t col, const fq_zech_ctx_t ctx)

Convert between `a` and the $d \times 1$ column vector of the
coefficients of its polynomial representative.
:::

## Comparison {#comparison}

::: function
int fq_zech_equal(const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
int fq_zech_is_zero(const fq_zech_t op, const fq_zech_ctx_t ctx)
int fq_zech_is_one(const fq_zech_t op, const fq_zech_ctx_t ctx)

Return whether the logarithms are equal, whether `op` is 0, respectively
1.
:::

::: function
int fq_zech_is_invertible(const fq_zech_t op, const fq_zech_ctx_t ctx)
int fq_zech_is_invertible_f(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)

Return whether `op` is nonzero. The second version sets `rop` to one if
it is (a Zech context always describes a field).
:::

## Arithmetic {#arithmetic}

::: function
void fq_zech_add(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
void fq_zech_sub(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)

Set `rop` to $\mathit{op1} \pm \mathit{op2}$, with one lookup in the Zech
table: $g^i + g^j = g^{j + Z(i - j)}$, and $g^i - g^j = g^i + g^{j +
(q-1)/2}$ for odd $p$.
:::

::: function
void fq_zech_sub_one(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
void fq_zech_neg(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)

Set `rop` to $\mathit{op1} - 1$, respectively $-\mathit{op1}$ (adding
$(q - 1)/2$ to the logarithm, or nothing if $p = 2$).
:::

::: function
void fq_zech_mul(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)
void fq_zech_sqr(fq_zech_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)

Set `rop` to $\mathit{op1} \cdot \mathit{op2}$, respectively
$\mathit{op}^2$, by adding logarithms modulo $q - 1$.
:::

::: function
void fq_zech_mul_fmpz(fq_zech_t rop, const fq_zech_t op, const fmpz_t x, const fq_zech_ctx_t ctx)
void fq_zech_mul_si(fq_zech_t rop, const fq_zech_t op, slong x, const fq_zech_ctx_t ctx)
void fq_zech_mul_ui(fq_zech_t rop, const fq_zech_t op, ulong x, const fq_zech_ctx_t ctx)

Set `rop` to `op` times the image of the integer $x$ in $\mathbb{F}_p$.
:::

::: function
void fq_zech_inv(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
void fq_zech_div(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2, const fq_zech_ctx_t ctx)

Set `rop` to $\mathit{op1}^{-1}$ (logarithm $q - 1 - k$), respectively
$\mathit{op1} / \mathit{op2}$. Throw `FLINT_ERROR` on division by zero.
:::

::: function
void fq_zech_gcdinv(fq_zech_t rop, fq_zech_t inv, const fq_zech_t op, const fq_zech_ctx_t ctx)

For nonzero `op`, sets `rop` to one and `inv` to the inverse; for zero
`op`, sets `rop` to zero.
:::

::: function
void fq_zech_pow(fq_zech_t rop, const fq_zech_t op1, const fmpz_t e, const fq_zech_ctx_t ctx)
void fq_zech_pow_ui(fq_zech_t rop, const fq_zech_t op1, const ulong e, const fq_zech_ctx_t ctx)

Set `rop` to $\mathit{op1}^e$ for $e \ge 0$ (with $\mathit{op1}^0 = 1$),
by multiplying the logarithm by $e$ modulo $q - 1$. `fq_zech_pow` throws
`FLINT_ERROR` for $e < 0$.
:::

## Roots {#roots}

::: function
int fq_zech_sqrt(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
int fq_zech_is_square(const fq_zech_t op1, const fq_zech_ctx_t ctx)

For odd $p$, an element is a square if and only if its logarithm is
even, and then half of it is the logarithm of a square root. For $p = 2$
every element is a square, with square root of logarithm $k/2$ or
$(k + q - 1)/2$. `fq_zech_sqrt` returns 1 and sets `rop` to a square root
if `op1` is a square, and returns 0 otherwise.
:::

::: function
void fq_zech_pth_root(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)

Sets `rop` to the $p$-th root of `op1`, by multiplying the logarithm by
$p^{d-1}$ modulo $q - 1$.
:::

## Frobenius map, norm, trace and order {#frobenius}

::: function
void fq_zech_frobenius(fq_zech_t rop, const fq_zech_t op, slong e, const fq_zech_ctx_t ctx)

Sets `rop` to $\mathit{op}^{p^e}$ (any $e$, taken modulo $d$), by
multiplying the logarithm by $p^e \bmod (q - 1)$.
:::

::: function
void fq_zech_norm(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)

Sets `rop` to the norm of $\mathit{op} = g^k$, which is $N(g)^k$, where
$N(g) = (-1)^d f(0)$ for the monic defining polynomial $f$ is stored in
the context (`prime_root`); one modular exponentiation in
$\mathbb{F}_p$.
:::

::: function
void fq_zech_trace(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)

Sets `rop` to the trace of `op`, computed as the sum of the $d$
conjugates $\mathit{op}^{p^i}$ with Zech additions, and converted from a
logarithm $l = m (q - 1)/(p - 1)$ in the prime field to the integer
$N(g)^m \bmod p$.
:::

::: function
int fq_zech_multiplicative_order(fmpz * ord, const fq_zech_t op, const fq_zech_ctx_t ctx)
int fq_zech_is_primitive(const fq_zech_t op, const fq_zech_ctx_t ctx)

Multiplicative order (returning 0 for zero, 1 for a primitive element and
$-1$ otherwise) and primitivity test, as `fq_multiplicative_order` and
`fq_is_primitive`; the order of $g^k$ is $(q - 1)/\gcd(k, q - 1)$, but
the generic algorithm, which factors $q - 1$, is used.
:::

## Random elements {#random}

::: function
void fq_zech_rand(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
void fq_zech_rand_not_zero(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
void fq_zech_randtest(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)
void fq_zech_randtest_not_zero(fq_zech_t rop, flint_rand_t state, const fq_zech_ctx_t ctx)

Set `rop` to a uniformly random element, respectively nonzero element;
the `randtest` versions use a different random generator.
:::

## Input and output {#io}

::: function
int fq_zech_fprint(FILE * file, const fq_zech_t op, const fq_zech_ctx_t ctx)
void fq_zech_print(const fq_zech_t op, const fq_zech_ctx_t ctx)
char * fq_zech_get_str(const fq_zech_t op, const fq_zech_ctx_t ctx)

Print, respectively return as a newly allocated string, the logarithm of
`op` as a decimal integer ($q - 1$ for zero).
:::

::: function
int fq_zech_fprint_pretty(FILE * file, const fq_zech_t op, const fq_zech_ctx_t ctx)
void fq_zech_print_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx)
char * fq_zech_get_str_pretty(const fq_zech_t op, const fq_zech_ctx_t ctx)

Print, respectively return as a newly allocated string, `op` as a power
of the generator, such as `a^17`, and zero as `0`.
:::

## Bit packing {#bit-packing}

::: function
void fq_zech_bit_pack(fmpz_t f, const fq_zech_t op, flint_bitcnt_t bit_size, const fq_zech_ctx_t ctx)
void fq_zech_bit_unpack(fq_zech_t rop, const fmpz_t f, flint_bitcnt_t bit_size, const fq_zech_ctx_t ctx)

Pack the coefficients of the polynomial representative of `op` into
`bit_size`-bit fields of an integer, and unpack and convert back, via the
`fq_nmod` representation, as `fq_bit_pack` and `fq_bit_unpack`.
:::
