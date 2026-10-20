---
title: fq_nmod_poly_factor.h
header: fq_nmod_poly_factor.h
template-headers: fq_poly_factor_templates
description: Factorisation of polynomials over finite fields $\mathbb{F}_q$ in the word-size characteristic.
category: Finite fields
order: 15
---

This module factors polynomials over a finite field $\mathbb{F}_q$ in the representation of [](fq_nmod) (a prime $p < 2^{\bits}$, elements stored as `nmod_poly`), tests irreducibility and squarefreeness, and finds roots. It is generated from the same template as [](fq_poly_factor) and has the same functions, with `fq_nmod` in place of `fq`; each function links to the `fq_poly_factor` version, which has the details. All algorithms are those of the generic module [](gr_poly), see [](algorithms/polynomials/factoring-finite-fields).

::: summary algorithms/polynomials/factoring-finite-fields :::

## Factorisation structure {#type}

::: type
fq_nmod_poly_factor_struct
fq_nmod_poly_factor_t

A factorisation: an array `poly` of `alloc` initialised polynomials and an
array `exp` of exponents, of which the first `num` pairs are used. The
factored polynomial is the product of `poly[i]^exp[i]` (times a leading
coefficient that is returned separately). `fq_nmod_poly_factor_t` is an array
of one structure.
:::

::: function
void fq_nmod_poly_factor_init(fq_nmod_poly_factor_t fac, const fq_nmod_ctx_t ctx)

Initialises `fac` as an empty factorisation, with room for 5 factors.

::: details
Generated from the same template as `fq_poly_factor_init`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_clear(fq_nmod_poly_factor_t fac, const fq_nmod_ctx_t ctx)

Frees the memory used by `fac`.

::: details
Generated from the same template as `fq_poly_factor_clear`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_realloc(fq_nmod_poly_factor_t fac, slong alloc, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_factor_fit_length(fq_nmod_poly_factor_t fac, slong len, const fq_nmod_ctx_t ctx)

`fq_nmod_poly_factor_realloc` changes the allocation to `alloc` factors
(clearing factors beyond it; `alloc` = 0 reinitialises `fac`).
`fq_nmod_poly_factor_fit_length` ensures room for `len` factors, at least
doubling the allocation when it grows.

::: details
Generated from the same template as `fq_poly_factor_realloc`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_set(fq_nmod_poly_factor_t res, const fq_nmod_poly_factor_t fac, const fq_nmod_ctx_t ctx)

Sets `res` to a copy of `fac`.

::: details
Generated from the same template as `fq_poly_factor_set`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_insert(fq_nmod_poly_factor_t fac, const fq_nmod_poly_t poly, slong exp, const fq_nmod_ctx_t ctx)

Multiplies the factorisation by `poly^exp`: if `poly` is already one of
the factors, `exp` is added to its exponent, otherwise `poly` is appended
with exponent `exp`. Constant polynomials are ignored.

::: details
Generated from the same template as `fq_poly_factor_insert`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_concat(fq_nmod_poly_factor_t res, const fq_nmod_poly_factor_t fac, const fq_nmod_ctx_t ctx)

Inserts all factors of `fac` into `res` with `fq_nmod_poly_factor_insert`.

::: details
Generated from the same template as `fq_poly_factor_concat`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_pow(fq_nmod_poly_factor_t fac, slong exp, const fq_nmod_ctx_t ctx)

Multiplies all exponents by `exp`, raising the factored polynomial to the
power `exp`.

::: details
Generated from the same template as `fq_poly_factor_pow`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_get_poly(fq_nmod_poly_t z, const fq_nmod_poly_factor_t fac, slong i, const fq_nmod_ctx_t ctx)

Sets `z` to the factor of index $i$, $0 \le i < \mathit{num}$.

::: details
Generated from the same template as `fq_poly_factor_get_poly`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_print(const fq_nmod_poly_factor_t fac, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_factor_print_pretty(const fq_nmod_poly_factor_t fac, const char * var, const fq_nmod_ctx_t ctx)

Print the factors to `stdout`, one per line, each followed by ` ^ ` and
its exponent, in the format of `fq_nmod_poly_print`, respectively of
`fq_nmod_poly_print_pretty` with the variable name `var`.

::: details
Generated from the same template as `fq_poly_factor_print`, whose documentation has the details.
:::
:::

## Squarefree and irreducible polynomials {#tests}

::: function
int _fq_nmod_poly_is_squarefree(const fq_nmod_struct * f, slong len, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_is_squarefree(const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)

Return whether $f$ is squarefree, that is, not divisible by the square of
a polynomial of positive degree. Polynomials of length at most 2 other
than 0 are squarefree, and 0 is not. This is a wrapper of
`gr_poly_is_squarefree`, which computes $\gcd(f, f')$ and handles
$p$-th powers.

::: details
Generated from the same template as `_fq_poly_is_squarefree`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_squarefree(fq_nmod_poly_factor_t res, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)

Sets `res` to a squarefree factorisation of $f$: pairwise coprime monic
squarefree polynomials $g_i$ with distinct exponents $e_i$ such that
$f = c \prod g_i^{e_i}$, where $c$ is the leading coefficient of $f$.
Previous contents of `res` are discarded; for $f$ of length at most 1 the
factorisation is empty.

::: details
Generated from the same template as `fq_poly_factor_squarefree`, whose documentation has the details.
:::
:::

::: function
int fq_nmod_poly_is_irreducible(const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_is_irreducible_ddf(const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
int fq_nmod_poly_is_irreducible_ben_or(const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)

Return whether $f$ is irreducible. Unlike `gr_poly_is_irreducible`, these
functions return 1 for all polynomials of length at most 2, including
constants and zero.

::: details
Generated from the same template as `fq_poly_is_irreducible`, whose documentation has the details.
:::
:::

## Factorisation {#factorisation}

::: function
void fq_nmod_poly_factor(fq_nmod_poly_factor_t result, fq_nmod_t leading_coeff, const fq_nmod_poly_t input, const fq_nmod_ctx_t ctx)

Factors `input` as $c \prod g_i^{e_i}$ with distinct monic irreducible
$g_i$: sets `result` to the factors and exponents (discarding previous
contents) and `leading_coeff` to $c$. For the zero polynomial $c = 0$ and
the factorisation is empty; for a constant, $c$ is the constant.

::: details
Generated from the same template as `fq_poly_factor`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_with_cantor_zassenhaus(fq_nmod_poly_factor_t result, fq_nmod_t leading_coeff, const fq_nmod_poly_t input, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_factor_with_berlekamp(fq_nmod_poly_factor_t result, fq_nmod_t leading_coeff, const fq_nmod_poly_t input, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_factor_with_kaltofen_shoup(fq_nmod_poly_factor_t result, fq_nmod_t leading_coeff, const fq_nmod_poly_t input, const fq_nmod_ctx_t ctx)

As `fq_nmod_poly_factor`, but factoring the squarefree parts with the
Cantor–Zassenhaus algorithm (distinct-degree factorisation followed by
equal-degree splitting), Berlekamp's algorithm (a basis of the kernel of
the Frobenius map minus the identity, followed by random splitting), or
the Kaltofen–Shoup algorithm (distinct-degree factorisation with
baby-step giant-step Frobenius computations by modular composition).

::: details
Generated from the same template as `fq_poly_factor_with_cantor_zassenhaus`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_cantor_zassenhaus(fq_nmod_poly_factor_t res, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_factor_berlekamp(fq_nmod_poly_factor_t factors, const fq_nmod_poly_t f, const fq_nmod_ctx_t ctx)
void fq_nmod_poly_factor_kaltofen_shoup(fq_nmod_poly_factor_t res, const fq_nmod_poly_t poly, const fq_nmod_ctx_t ctx)

As the `_with_` functions, without returning the leading coefficient:
the output is the factorisation of the monic polynomial $f / c$. Previous
contents of the output are discarded.

::: details
Generated from the same template as `fq_poly_factor_cantor_zassenhaus`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_distinct_deg(fq_nmod_poly_factor_t res, const fq_nmod_poly_t poly, slong * const * degs, const fq_nmod_ctx_t ctx)

Distinct-degree factorisation of the squarefree polynomial `poly`:
appends to `res`, with exponent 1, polynomials $g_1, g_2, \ldots$ such
that $\operatorname{monic}(\mathit{poly}) = \prod g_i$ and each $g_i$ is
the product of all irreducible factors of degree $d_i$, and writes the
degrees $d_i$ to `(*degs)[0]`, `(*degs)[1]`, .... The caller must
provide room for $\deg(\mathit{poly})$ entries in `*degs`. The result is
undefined if `poly` is not squarefree. This is a wrapper of
`gr_poly_factor_distinct_deg`.

::: details
Generated from the same template as `fq_poly_factor_distinct_deg`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_equal_deg(fq_nmod_poly_factor_t factors, const fq_nmod_poly_t pol, slong d, const fq_nmod_ctx_t ctx)

Equal-degree factorisation: appends to `factors`, with exponent 1, the
irreducible factors of `pol`, which must be a monic squarefree product of
irreducible polynomials of degree $d$. The factors are found by random
splitting (Cantor–Zassenhaus), using an internal random state.

::: details
Generated from the same template as `fq_poly_factor_equal_deg`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_factor_split_single(fq_nmod_poly_t linfactor, const fq_nmod_poly_t input, const fq_nmod_ctx_t ctx)

Sets `linfactor` to one monic linear factor $x - a$ of `input`, which must
split into linear factors (not necessarily distinct). If `input` has
degree 1 it is copied. Throws `FLINT_ERROR` if no root is found. This is
the step used by `fq_nmod_embed_gens` to find a root of a defining
polynomial in an extension field.

::: details
Generated from the same template as `fq_poly_factor_split_single`, whose documentation has the details.
:::
:::

::: function
void fq_nmod_poly_iterated_frobenius_preinv(fq_nmod_poly_t * rop, slong n, const fq_nmod_poly_t v, const fq_nmod_poly_t vinv, const fq_nmod_ctx_t ctx)

Sets `rop[i]` to $x^{q^i} \bmod v$ for $0 \le i < n$, where `vinv` is the
inverse of the reversal of $v$ (see `fq_nmod_poly_reverse` and
`fq_nmod_poly_inv_series`). The array `rop` must contain $n$ initialised
polynomials. If `FQ_NMOD_POLY_ITERATED_FROBENIUS_CUTOFF` returns 1, the
powers after the first are computed by modular composition with a
precomputed Brent–Kung matrix ($x^{q^{i}} = x^{q^{i-1}}(x^q)$), otherwise
by exponentiation.

::: details
Generated from the same template as `fq_poly_iterated_frobenius_preinv`, whose documentation has the details.
:::
:::

::: function
int FQ_NMOD_POLY_ITERATED_FROBENIUS_CUTOFF(const fq_nmod_ctx_t ctx, slong length)

Returns 1 if `fq_nmod_poly_iterated_frobenius_preinv` should use modular
composition for a modulus of length `length`, which is the case if
$2 \log_2 q < 3(\lfloor \sqrt{\mathit{length}} \rfloor + 1)$.
:::

## Roots {#roots}

::: function
void fq_nmod_poly_roots(fq_nmod_poly_factor_t r, const fq_nmod_poly_t f, int with_multiplicity, const fq_nmod_ctx_t ctx)

Sets `r` to the roots of the nonzero polynomial $f$ in $\mathbb{F}_q$,
stored as the monic linear factors $x - a$. If `with_multiplicity` is
nonzero, the exponents are the multiplicities of the roots; otherwise
they are 1. Previous contents of `r` are discarded. This is a wrapper of
`gr_poly_roots_finite_field`, which splits $\gcd(f, x^q - x)$ by
equal-degree factorisation; see
[](algorithms/polynomials/root-finding).

::: details
Generated from the same template as `fq_poly_roots`, whose documentation has the details.
:::
:::

## Example {#example}
