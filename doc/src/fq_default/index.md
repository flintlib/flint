---
title: fq_default.h
header: fq_default.h
description: Finite fields $\mathbb{F}_q$ of any order, with the representation (Zech logarithms, word-size or multiprecision characteristic, prime field) chosen automatically.
category: Finite fields
order: 1
---

This module is the recommended interface to finite fields $\mathbb{F}_q$,
$q = p^d$. A context `fq_default_ctx_t` describes the field, and elements
are `fq_default_t` values whose meaning depends on the context. When the
context is created, the module chooses one of five representations: Zech
logarithms (`fq_zech`) for small fields, polynomials over a word-size
prime field (`fq_nmod`), polynomials over a multiprecision prime field
(`fq`), and, for prime fields ($d = 1$), plain residues modulo $p$
(`nmod` or `fmpz_mod`). Every function dispatches on the representation,
at the cost of a few comparisons or an indirect call per operation.
Polynomials, matrices and factorisation over the same fields are provided
by `fq_default_poly.h`, `fq_default_mat.h` and
`fq_default_poly_factor.h`.

The guide [](fq_default/getting-started) shows how to create fields and
compute with their elements. The representations and the algorithms
they use are compared on [](algorithms/finite-fields/overview).

::: summary algorithms/finite-fields/overview :::

## Types and representations {#types}

::: type
fq_default_ctx_struct
fq_default_ctx_t

A finite field context. The structure is a generic ring context
`gr_ctx_struct` (see [](gr)) for one of the rings `fq_zech`, `fq_nmod`,
`fq`, `nmod` or `fmpz_mod`; `fq_default_ctx_t` is an array of one
structure. A context must be initialised with one of the
`fq_default_ctx_init` functions and cleared with `fq_default_ctx_clear`.
Because it is a `gr` context, `FQ_DEFAULT_GR_CTX(ctx)` can be passed to
any generic function, and an `fq_default_t` is a valid element for it.
:::

::: type
fq_default_struct
fq_default_t

An element of a finite field: a union of an `fq_t`, an `fq_nmod_t`, an
`fq_zech_t`, a `ulong` (for `nmod`) and an `fmpz_t` (for `fmpz_mod`).
Which member is valid is determined by the context, which must be passed
to every function. `fq_default_t` is an array of one union, so it is
passed by reference.
:::

::: macro
FQ_DEFAULT_FQ_ZECH
FQ_DEFAULT_FQ_NMOD
FQ_DEFAULT_FQ
FQ_DEFAULT_NMOD
FQ_DEFAULT_FMPZ_MOD

The type codes 1, 2, 3, 4 and 5 returned by `fq_default_ctx_type` and
accepted by the `_type` initialisation functions. They select the
representations of `fq_zech.h` (Zech logarithms), `fq_nmod.h`
(polynomials over $\mathbb{F}_p$, $p$ word size), `fq.h` (polynomials
over $\mathbb{F}_p$, $p$ multiprecision), `nmod.h` ($\mathbb{F}_p$,
$p$ word size) and `fmpz_mod.h` ($\mathbb{F}_p$, any $p$).
:::

Unless a type is requested, a context for $\mathbb{F}_{p^d}$ uses

- `FQ_DEFAULT_FQ_ZECH` if $d > 1$ and $\operatorname{bits}(p) \cdot d
  \le 16$, so that $q < 2^{16}$ and the tables of the Zech
  representation take at most about a megabyte;
- `FQ_DEFAULT_FQ_NMOD` if $d > 1$ and $p < 2^{\bits}$;
- `FQ_DEFAULT_NMOD` if $d = 1$ and $p < 2^{\bits}$;
- `FQ_DEFAULT_FMPZ_MOD` if $d = 1$ and $p \ge 2^{\bits}$;
- `FQ_DEFAULT_FQ` otherwise.

For example, $\mathbb{F}_{2^8}$, $\mathbb{F}_{3^8}$ and $\mathbb{F}_{7^5}$
use Zech logarithms, while $\mathbb{F}_{2^{16}}$ and $\mathbb{F}_{65537^2}$
use `fq_nmod`.

The rules for elements are those of the underlying modules: every
element must be initialised for the context in which it is used, outputs
may be the same objects as inputs, and the result of mixing elements of
different contexts is undefined. Nothing is checked.

## Context management {#context}

::: function
void fq_default_ctx_init(fq_default_ctx_t ctx, const fmpz_t p, slong d, const char * var)
void fq_default_ctx_init_type(fq_default_ctx_t ctx, const fmpz_t p, slong d, const char * var, int type)

Initialise `ctx` as a context for $\mathbb{F}_{p^d}$, for a prime $p$
and $d \ge 1$, with the name `var` for the generator (used for
printing). The defining polynomial is a Conway polynomial if $(p, d)$ is
in FLINT's table; otherwise it is chosen as in the underlying module.
`fq_default_ctx_init` chooses the representation by the rules
[above](#types); `fq_default_ctx_init_type` uses the given `type`, one
of the five type codes, or 0 for the automatic choice.

::: details
Neither the primality of $p$ nor the compatibility of `type` with $p$ and
$d$ is checked. `FQ_DEFAULT_NMOD` and `FQ_DEFAULT_FMPZ_MOD` create the
prime field $\mathbb{F}_p$ whatever $d$ is; `FQ_DEFAULT_FQ_ZECH`,
`FQ_DEFAULT_FQ_NMOD` and `FQ_DEFAULT_NMOD` need $p < 2^{\bits}$, and
`FQ_DEFAULT_FQ_ZECH` also needs $q = p^d < 2^{\bits}$ and builds tables
of $2q + p$ words.

If no Conway polynomial is available (see
[](algorithms/finite-fields/defining-polynomials)), `fq_zech` uses a
random primitive polynomial, `fq_nmod` a polynomial with the fewest
possible nonzero terms (`fq_nmod_ctx_init_ui`) and `fq` a random sparse
irreducible polynomial (`fq_ctx_init`). A prime field created by these
functions has the defining polynomial $x$, so its generator
`fq_default_gen` is 0.
:::
:::

::: function
void fq_default_ctx_init_modulus(fq_default_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var)
void fq_default_ctx_init_modulus_type(fq_default_ctx_t ctx, const fmpz_mod_poly_t modulus, fmpz_mod_ctx_t mod_ctx, const char * var, int type)
void fq_default_ctx_init_modulus_nmod(fq_default_ctx_t ctx, const nmod_poly_t modulus, const char * var)
void fq_default_ctx_init_modulus_nmod_type(fq_default_ctx_t ctx, const nmod_poly_t modulus, const char * var, int type)

Initialise `ctx` as a context for $\mathbb{F}_p[x]/(f)$, where $f$ is
`modulus`, an irreducible polynomial of degree $d \ge 1$ over
$\mathbb{F}_p$, and $p$ is the modulus of `mod_ctx`, respectively of the
`nmod_poly`. The generator, printed as `var`, is the class of $x$. The
representation is chosen as by `fq_default_ctx_init`, or given by
`type` (0 for the automatic choice).

::: details
The polynomial need not be monic. Its irreducibility and the primality
of $p$ are not checked.

The Zech representation needs a *primitive* polynomial, one whose roots
generate the multiplicative group. If the Zech representation is chosen
or requested and $f$ is not primitive, these functions fall back to
`FQ_DEFAULT_FQ_NMOD`; the `fmpz_mod_poly` versions do the same when $p$
does not fit in a word, which `fq_nmod` cannot handle, so a Zech
representation must not be requested for such $p$. For $d = 1$ with the
type `FQ_DEFAULT_NMOD` or `FQ_DEFAULT_FMPZ_MOD`, the field is
$\mathbb{F}_p$ and the generator is the root $a = -c_0/c_1$ of
$f = c_1 x + c_0$.
:::
:::

::: function
void fq_default_ctx_clear(fq_default_ctx_t ctx)

Frees the memory used by `ctx`. Elements, polynomials and matrices over
the field must be cleared first.
:::

::: function
int fq_default_ctx_type(const fq_default_ctx_t ctx)

Returns the representation used by `ctx` as one of the type codes
`FQ_DEFAULT_FQ_ZECH` (1), `FQ_DEFAULT_FQ_NMOD` (2), `FQ_DEFAULT_FQ` (3),
`FQ_DEFAULT_NMOD` (4) and `FQ_DEFAULT_FMPZ_MOD` (5).
:::

::: function
void * fq_default_ctx_inner(const fq_default_ctx_t ctx)

Returns a pointer to the context of the underlying representation: an
`fq_zech_ctx_struct *`, `fq_nmod_ctx_struct *`, `fq_ctx_struct *`,
`nmod_t *` or `fmpz_mod_ctx_struct *`, according to
`fq_default_ctx_type`. Together with the corresponding member of
`fq_default_struct`, this gives direct access to the functions of the
underlying module.
:::

::: function
slong fq_default_ctx_degree(const fq_default_ctx_t ctx)

Returns the degree $d = [\mathbb{F}_q : \mathbb{F}_p]$.
:::

::: function
void fq_default_ctx_prime(fmpz_t prime, const fq_default_ctx_t ctx)

Sets `prime` to the characteristic $p$.
:::

::: function
void fq_default_ctx_order(fmpz_t f, const fq_default_ctx_t ctx)

Sets `f` to the number of elements $q = p^d$.
:::

::: function
void fq_default_ctx_modulus(fmpz_mod_poly_t p, const fq_default_ctx_t ctx)

Sets `p` to the defining polynomial of the field. The `fmpz_mod_poly`
must have been initialised with a context whose modulus is the
characteristic. For a prime field the result is $x - a$, where $a$ is the
generator.
:::

::: function
int fq_default_ctx_fprint(FILE * file, const fq_default_ctx_t ctx)
void fq_default_ctx_print(const fq_default_ctx_t ctx)

Print a description of the field to `file`, respectively to `stdout`: the
characteristic, and for extension fields the degree and the defining
polynomial, preceded by `Zech Representation:` for the Zech
representation. `fq_default_ctx_fprint` returns a positive value on
success.
:::

::: function
void fq_default_ctx_init_randtest(fq_default_ctx_t ctx, flint_rand_t state)

Initialises `ctx` as a random finite field in a random one of the five
representations, for testing: up to 4-bit primes and degree 3 for Zech
logarithms, up to 30-bit primes and degree 15 for `fq_nmod` and `nmod`,
and up to 63-bit primes for `fq` (degree up to 7) and `fmpz_mod`.
:::

## Memory management {#memory}

::: function
void fq_default_init(fq_default_t rop, const fq_default_ctx_t ctx)
void fq_default_init2(fq_default_t rop, const fq_default_ctx_t ctx)

Initialise `rop` for the field `ctx` and set it to zero.
`fq_default_init2` also allocates space for $d$ coefficients in the
polynomial representations, so that later operations need not allocate.
:::

::: function
void fq_default_clear(fq_default_t rop, const fq_default_ctx_t ctx)

Frees the memory used by `rop`.
:::

## Assignment and conversion {#assignment}

::: function
void fq_default_set(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)

Sets `rop` to `op`.
:::

::: function
void fq_default_swap(fq_default_t op1, fq_default_t op2, const fq_default_ctx_t ctx)

Swaps `op1` and `op2` efficiently.
:::

::: function
void fq_default_zero(fq_default_t rop, const fq_default_ctx_t ctx)
void fq_default_one(fq_default_t rop, const fq_default_ctx_t ctx)

Set `rop` to 0, respectively 1.
:::

::: function
void fq_default_gen(fq_default_t rop, const fq_default_ctx_t ctx)

Sets `rop` to the generator of the field: the class of $x$ in
$\mathbb{F}_p[x]/(f)$, a root of the defining polynomial $f$. It
generates $\mathbb{F}_q$ as an algebra over $\mathbb{F}_p$; it generates
the multiplicative group if $f$ is primitive, which is the case for
Conway polynomials and always in the Zech representation.
:::

::: function
void fq_default_set_ui(fq_default_t rop, const ulong x, const fq_default_ctx_t ctx)
void fq_default_set_si(fq_default_t rop, const slong x, const fq_default_ctx_t ctx)
void fq_default_set_fmpz(fq_default_t rop, const fmpz_t x, const fq_default_ctx_t ctx)

Set `rop` to the image of the integer $x$ in the prime field
$\mathbb{F}_p \subseteq \mathbb{F}_q$.
:::

::: function
int fq_default_get_fmpz(fmpz_t z, const fq_default_t op, const fq_default_ctx_t ctx)

If `op` lies in the prime field $\mathbb{F}_p$, sets `z` to its
representative in $[0, p)$ and returns 1; otherwise returns 0 and leaves
`z` unchanged.
:::

::: function
void fq_default_get_coeff_fmpz(fmpz_t c, fq_default_t op, slong n, const fq_default_ctx_t ctx)

Sets `c` to the coefficient of $x^n$ in the representative of `op` as a
polynomial of degree less than $d$ in the generator, or to zero if $n
\ge d$. For a prime field, the representative is a constant.
:::

::: function
void fq_default_get_nmod_poly(nmod_poly_t poly, const fq_default_t op, const fq_default_ctx_t ctx)
void fq_default_set_nmod_poly(fq_default_t op, const nmod_poly_t poly, const fq_default_ctx_t ctx)

Convert between `op` and its representative of degree less than $d$ in
$\mathbb{F}_p[x]$. `fq_default_set_nmod_poly` accepts a polynomial of
any degree and reduces it modulo the defining polynomial. The modulus of
`poly` must be $p$; it is not checked. Only for the types
`FQ_DEFAULT_FQ_ZECH`, `FQ_DEFAULT_FQ_NMOD` and `FQ_DEFAULT_NMOD`; the
other types throw `FLINT_ERROR`.
:::

::: function
void fq_default_get_fmpz_mod_poly(fmpz_mod_poly_t poly, const fq_default_t op, const fq_default_ctx_t ctx)
void fq_default_set_fmpz_mod_poly(fq_default_t op, const fmpz_mod_poly_t poly, const fq_default_ctx_t ctx)

Convert between `op` and its representative in $\mathbb{F}_p[x]$, as
the `nmod_poly` versions but for all representations. The modulus of the
`fmpz_mod_poly` context must be $p$; it is not checked.
:::

::: function
void fq_default_get_fmpz_poly(fmpz_poly_t a, const fq_default_t b, const fq_default_ctx_t ctx)
void fq_default_set_fmpz_poly(fq_default_t a, const fmpz_poly_t b, const fq_default_ctx_t ctx)

Convert between an element and an integer polynomial.
`fq_default_get_fmpz_poly` sets `a` to the representative of `b` of
degree less than $d$, with coefficients in $[0, p)$.
`fq_default_set_fmpz_poly` reduces the coefficients of `b` modulo $p$ and
the polynomial modulo the defining polynomial.
:::

## Comparison {#comparison}

::: function
int fq_default_is_zero(const fq_default_t op, const fq_default_ctx_t ctx)
int fq_default_is_one(const fq_default_t op, const fq_default_ctx_t ctx)

Return whether `op` is 0, respectively 1.
:::

::: function
int fq_default_equal(const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)

Returns whether `op1` and `op2` are equal.
:::

::: function
int fq_default_is_invertible(const fq_default_t op, const fq_default_ctx_t ctx)

Returns whether `op` is invertible, that is, nonzero.
:::

## Arithmetic {#arithmetic}

::: function
void fq_default_add(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)
void fq_default_sub(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)

Set `rop` to $\mathit{op1} + \mathit{op2}$, respectively
$\mathit{op1} - \mathit{op2}$.
:::

::: function
void fq_default_sub_one(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)

Sets `rop` to $\mathit{op1} - 1$.
:::

::: function
void fq_default_neg(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)

Sets `rop` to $-\mathit{op1}$.
:::

::: function
void fq_default_mul(fq_default_t rop, const fq_default_t op1, const fq_default_t op2, const fq_default_ctx_t ctx)
void fq_default_sqr(fq_default_t rop, const fq_default_t op, const fq_default_ctx_t ctx)

Set `rop` to $\mathit{op1} \cdot \mathit{op2}$, respectively
$\mathit{op}^2$.
:::

::: function
void fq_default_mul_fmpz(fq_default_t rop, const fq_default_t op, const fmpz_t x, const fq_default_ctx_t ctx)
void fq_default_mul_si(fq_default_t rop, const fq_default_t op, slong x, const fq_default_ctx_t ctx)
void fq_default_mul_ui(fq_default_t rop, const fq_default_t op, ulong x, const fq_default_ctx_t ctx)

Set `rop` to `op` times the image of the integer $x$ in $\mathbb{F}_p$.
:::

::: function
void fq_default_inv(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)

Sets `rop` to $\mathit{op1}^{-1}$. Throws `FLINT_ERROR` if `op1` is
zero.
:::

::: function
void fq_default_div(fq_default_t rop, fq_default_t op1, fq_default_t op2, const fq_default_ctx_t ctx)

Sets `rop` to $\mathit{op1} / \mathit{op2}$. Throws `FLINT_ERROR` if
`op2` is zero. The inputs are not modified, although they are not
declared `const`.
:::

::: function
void fq_default_pow(fq_default_t rop, const fq_default_t op1, const fmpz_t e, const fq_default_ctx_t ctx)
void fq_default_pow_ui(fq_default_t rop, const fq_default_t op, const ulong e, const fq_default_ctx_t ctx)

Set `rop` to $\mathit{op}^e$, with $0^0 = 1$. In `fq_default_pow` the
exponent may be negative if `op` is nonzero; $0^e$ with $e < 0$ throws
`FLINT_ERROR`.

::: details
The power is computed by `gr_pow_fmpz` or `gr_pow_ui` for the
representation. The polynomial representations reduce $e$ modulo $q - 1$
when $e \ge p$ and use binary exponentiation, which costs
$O(\log \min(e, q))$ multiplications in $\mathbb{F}_q$; the Zech
representation needs one multiplication of the discrete logarithm by $e$
modulo $q - 1$.
:::
:::

## Roots {#roots}

::: function
int fq_default_sqrt(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)

If `op1` is a square, sets `rop` to a square root of it and returns 1.
Otherwise returns 0, and the value of `rop` is unspecified.
:::

::: function
int fq_default_is_square(const fq_default_t op, const fq_default_ctx_t ctx)

Returns whether `op` is a square. If $p = 2$, every element is a square.
:::

::: function
void fq_default_pth_root(fq_default_t rop, const fq_default_t op1, const fq_default_ctx_t ctx)

Sets `rop` to the unique $p$-th root of `op1`, which is
$\mathit{op1}^{p^{d-1}}$, the image under the inverse of the Frobenius
map.
:::

## Frobenius map, norm and trace {#frobenius}

The Galois group of $\mathbb{F}_q/\mathbb{F}_p$ is cyclic of order $d$,
generated by the Frobenius automorphism $\sigma(a) = a^p$. The norm and
the trace of $a \in \mathbb{F}_q$ are
$N(a) = \prod_{i=0}^{d-1} \sigma^i(a)$ and
$\operatorname{Tr}(a) = \sum_{i=0}^{d-1} \sigma^i(a)$, which are the
determinant and the trace of the $\mathbb{F}_p$-linear map
$x \mapsto ax$. The algorithms are described in
[](algorithms/finite-fields/frobenius-norm-trace).

::: function
void fq_default_frobenius(fq_default_t rop, const fq_default_t op, slong e, const fq_default_ctx_t ctx)

Sets `rop` to $\sigma^e(\mathit{op}) = \mathit{op}^{p^e}$. Any $e$ is
allowed, including negative values; it is taken modulo $d$.
:::

::: function
void fq_default_norm(fmpz_t rop, const fq_default_t op, const fq_default_ctx_t ctx)
void fq_default_trace(fmpz_t rop, const fq_default_t op, const fq_default_ctx_t ctx)

Set `rop` to the norm $N(\mathit{op})$, respectively the trace
$\operatorname{Tr}(\mathit{op})$, as an integer in $[0, p)$.
:::

## Random elements {#random}

::: function
void fq_default_rand(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)
void fq_default_rand_not_zero(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)

Set `rop` to a uniformly random element, respectively a uniformly random
nonzero element, of the field.
:::

::: function
void fq_default_randtest(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)
void fq_default_randtest_not_zero(fq_default_t rop, flint_rand_t state, const fq_default_ctx_t ctx)

Set `rop` to a random element, respectively a random nonzero element, for
testing. In the polynomial representations, representatives with many
zero coefficients are more likely than with uniform sampling.
:::

## Input and output {#io}

::: function
int fq_default_fprint(FILE * file, const fq_default_t op, const fq_default_ctx_t ctx)
void fq_default_print(const fq_default_t op, const fq_default_ctx_t ctx)

Print `op` to `file`, respectively `stdout`, in the plain format of the
representation: the format of `nmod_poly_print` (length, modulus,
coefficients from the constant term up) for `fq_nmod`, that of
`fmpz_poly_print` (length, coefficients) for `fq`, the discrete
logarithm $k$ of $\mathit{op} = g^k$ for the Zech representation (with
$q - 1$ for zero), and the residue for prime fields.
`fq_default_fprint` returns a positive value on success.
:::

::: function
int fq_default_fprint_pretty(FILE * file, const fq_default_t op, const fq_default_ctx_t ctx)
void fq_default_print_pretty(const fq_default_t op, const fq_default_ctx_t ctx)

Print `op` to `file`, respectively `stdout`, in a readable form using the
variable name of the context: a polynomial such as `a+3` in the
polynomial representations, a power such as `a^36` in the Zech
representation (`0` for zero), and the residue in prime fields.
`fq_default_fprint_pretty` returns a positive value on success.
:::

::: function
char * fq_default_get_str(const fq_default_t op, const fq_default_ctx_t ctx)
char * fq_default_get_str_pretty(const fq_default_t op, const fq_default_ctx_t ctx)

Return the plain, respectively readable, representation of `op` as a
newly allocated string, which the caller must free with `flint_free`.
:::

## Internal macros {#internal}

::: macro
FQ_DEFAULT_GR_CTX(ctx)
_FQ_DEFAULT_TYPE(ctx)

The context as a `gr_ctx_struct *`, and its field `which_ring`, the
internal type code (`GR_CTX_FQ_ZECH`, `GR_CTX_FQ_NMOD`, `GR_CTX_FQ`,
`GR_CTX_NMOD` or `GR_CTX_FMPZ_MOD`, also available as
`_FQ_DEFAULT_FQ_ZECH` and so on).
:::

::: macro
FQ_DEFAULT_CTX_FQ_ZECH(ctx)
FQ_DEFAULT_CTX_FQ_NMOD(ctx)
FQ_DEFAULT_CTX_FQ(ctx)
FQ_DEFAULT_CTX_NMOD(ctx)
FQ_DEFAULT_CTX_FMPZ_MOD(ctx)

The context of the underlying representation: a pointer to the
`fq_zech`, `fq_nmod`, `fq` or `fmpz_mod` context, or the `nmod_t` itself
for `FQ_DEFAULT_CTX_NMOD`. Valid only for the matching type.
:::

::: macro
FQ_DEFAULT_CTX_NMOD_A(ctx)
FQ_DEFAULT_CTX_FMPZ_MOD_A(ctx)

A pointer to the generator $a$ of a prime field context of type
`FQ_DEFAULT_NMOD` (a `ulong`), respectively `FQ_DEFAULT_FMPZ_MOD` (an
`fmpz`); the defining polynomial is $x - a$.
:::
