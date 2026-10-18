---
title: fmpz_mod_mpoly_factor.h
header: fmpz_mod_mpoly_factor.h
description: Factorisation of multivariate polynomials over $\mathbb{F}_p$ for a prime $p$ of any size, and the GCD and lifting machinery behind `fmpz_mod_mpoly`.
category: Multivariate polynomials
order: 10
---

This module factors polynomials of `fmpz_mod_mpoly.h` over a prime field
$\mathbb{F}_p$, where $p$ is given as an `fmpz`. Factorisations are stored
in an `fmpz_mod_mpoly_factor_t`, a constant times a product of powers of
polynomials, with the same interface as `nmod_mpoly_factor` (see
[](nmod_mpoly_factor)).

If $p$ fits in a word, the factorisation functions convert to `nmod_mpoly`
and use `nmod_mpoly_factor_algo`, which handles small fields through
extension fields. Otherwise $p > 2^{\bits}$ and random evaluation points
in $\mathbb{F}_p$ are good with overwhelming probability, so only the
algorithms over $\mathbb{F}_p$ itself are implemented here: content
removal, separable factorisation, bivariate factorisation by Hensel lifting
and recombination, and lifting by Wang's method, sparse Hensel lifting or
Zassenhaus recombination [@Wan1978] [@MonTun2020] [@Zas1969]. See
[](algorithms/multivariate/modp-factoring) for the methods and
[](algorithms/multivariate/hensel-lifting) for the lifting steps.

The header also declares the dense bivariate types and the GCD algorithms
used by `fmpz_mod_mpoly_gcd` ([](algorithms/multivariate/modp-gcd)),
documented in the internal section.

## Types and memory management {#types}

::: type
fmpz_mod_mpoly_factor_struct
fmpz_mod_mpoly_factor_t

A factorisation $c \prod_i f_i^{e_i}$, with the fields `constant` (an
`fmpz_t`), `poly` (the bases), `exp` (the `fmpz` exponents), `num` and
`alloc`.
:::

::: function
void fmpz_mod_mpoly_factor_init(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)

Initialises `f` as the empty product with constant 1.
:::

::: function
void fmpz_mod_mpoly_factor_clear(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_factor_realloc(fmpz_mod_mpoly_factor_t f, slong alloc, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_factor_fit_length(fmpz_mod_mpoly_factor_t f, slong len, const fmpz_mod_mpoly_ctx_t ctx)

Free `f`, change its allocation to `alloc` factors, or ensure space for
`len` factors.
:::

## Basic manipulation {#basic}

::: function
void fmpz_mod_mpoly_factor_swap(fmpz_mod_mpoly_factor_t f, fmpz_mod_mpoly_factor_t g, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_factor_one(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)

Swap two factorisations, and set `f` to the empty product with constant
1.
:::

::: function
slong fmpz_mod_mpoly_factor_length(const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_factor_get_constant_fmpz(fmpz_t c, const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)

Return the number of factors, and set `c` to the constant.
:::

::: function
void fmpz_mod_mpoly_factor_get_base(fmpz_mod_mpoly_t p, const fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_factor_swap_base(fmpz_mod_mpoly_t p, fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)
slong fmpz_mod_mpoly_factor_get_exp_si(fmpz_mod_mpoly_factor_t f, slong i, const fmpz_mod_mpoly_ctx_t ctx)

Set `p` to the base of index $i$ or swap it with `p`, and return the
exponent of index $i$, which must fit an `slong`.
:::

::: function
void fmpz_mod_mpoly_factor_append_fmpz_swap(fmpz_mod_mpoly_factor_t f, fmpz_mod_mpoly_t A, const fmpz_t e, const fmpz_mod_mpoly_ctx_t ctx)

Appends the factor $A^e$ to `f`, swapping `A` into it (so `A` is
changed).
:::

::: function
int fmpz_mod_mpoly_factor_expand(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_factor_matches(const fmpz_mod_mpoly_t a, const fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)

Set `A` to the product represented by `f` (returning 0 for negative or
too large exponents), and test whether `f` represents `a`.
:::

::: function
void fmpz_mod_mpoly_factor_sort(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_factor_cmp(const fmpz_mod_mpoly_factor_t A, const fmpz_mod_mpoly_factor_t B, const fmpz_mod_mpoly_ctx_t ctx)

Sort the factors by exponent and then by base, and compare two
factorisations (constant, number of factors, then factor by factor).
:::

::: function
void fmpz_mod_mpoly_factor_print_pretty(const fmpz_mod_mpoly_factor_t f, const char ** vars, const fmpz_mod_mpoly_ctx_t ctx)

Prints the constant followed by one line `*(base)^exp` per factor.
:::

## Factorisation {#factorisation}

These functions assume that the modulus is prime, which is not checked.
They set `f` to a factorisation of `A` whose constant is the leading
coefficient of `A` and whose bases are monic, and return 1 for success
(0 only when exponents do not fit one word per field after removing the
content).

::: function
int fmpz_mod_mpoly_factor_content(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Sets `f` to a factorisation into variables and polynomials primitive with
respect to every variable.
:::

::: function
int fmpz_mod_mpoly_factor_squarefree(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_factor_separable(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, int sep)

Set `f` to a factorisation with squarefree, pairwise coprime bases. With
`sep` = 1, `fmpz_mod_mpoly_factor_separable` also makes each base $g$
satisfy $\partial g/\partial x_i = 0$ or $\gcd(g, \partial g/\partial
x_i) = 1$ for every variable, as `nmod_mpoly_factor_separable`.
:::

::: function
int fmpz_mod_mpoly_factor(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)

Sets `f` to the factorisation of `A` into monic irreducible factors.
:::

::: function
int fmpz_mod_mpoly_factor_zassenhaus(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_factor_wang(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_factor_zippel(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_factor_algo(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, unsigned int algo)

As `fmpz_mod_mpoly_factor` with only one lifting method, or the methods
in the bit mask `algo` (`MPOLY_FACTOR_USE_ZAS`, `MPOLY_FACTOR_USE_WANG`,
`MPOLY_FACTOR_USE_ZIP`, `MPOLY_FACTOR_USE_ALL`).
:::

## Internal functions {#internal}

These functions support the factorisation and the GCD algorithms of
`fmpz_mod_mpoly`. Their interfaces may change without notice.

### Dense and recursive types

::: type
fmpz_mod_bpoly_struct
fmpz_mod_bpoly_t
fmpz_mod_tpoly_struct
fmpz_mod_tpoly_t

A dense bivariate polynomial: an array (`coeffs`, `alloc`, `length`) of
`fmpz_mod_poly_struct`, the coefficient of index $i$ being the coefficient
of $x^i$ as a polynomial in $y$; and an array of bivariate polynomials
(used for lists of factors).
:::

::: type
fmpz_mod_polyu_struct
fmpz_mod_polyu_t
fmpz_mod_polyun_struct
fmpz_mod_polyun_t

Sparse polynomials in up to three variables with packed exponents
(`exps`) and `fmpz` coefficients (`polyu`), or with `fmpz_mod_poly`
coefficients in one more variable (`polyun`).
:::

::: type
fmpz_mod_mpolyu_struct
fmpz_mod_mpolyu_t
fmpz_mod_mpolyn_struct
fmpz_mod_mpolyn_t

Sparse univariate polynomials with `fmpz_mod_mpoly` coefficients, and
multivariate polynomials with `fmpz_mod_poly` coefficients in the last
variable, as the corresponding `nmod_mpoly` types.
:::

::: type
fmpz_mod_mpolyv_struct
fmpz_mod_mpolyv_t
fmpz_mod_mpoly_pfrac_struct
fmpz_mod_mpoly_pfrac_t

Vectors of polynomials (`coeffs`, `alloc`, `length`), and precomputed
data for multivariate partial fractions.
:::

::: type
fmpz_mod_poly_stack_struct
fmpz_mod_poly_stack_t
fmpz_mod_bpoly_stack_struct
fmpz_mod_bpoly_stack_t
fmpz_mod_polyun_stack_struct
fmpz_mod_polyun_stack_t
fmpz_mod_mpolyn_stack_struct
fmpz_mod_mpolyn_stack_t
fmpz_mod_poly_bpoly_stack_struct
fmpz_mod_poly_bpoly_stack_t
fmpz_mod_poly_polyun_stack_struct
fmpz_mod_poly_polyun_stack_t
fmpz_mod_poly_polyun_mpolyn_stack_struct
fmpz_mod_poly_polyun_mpolyn_stack_t

Stacks of preallocated temporaries, and bundles of such stacks.
:::

### Stacks

::: function
void fmpz_mod_poly_stack_init(fmpz_mod_poly_stack_t S)
void fmpz_mod_poly_stack_clear(fmpz_mod_poly_stack_t S)
fmpz_mod_poly_struct ** fmpz_mod_poly_stack_fit_request(fmpz_mod_poly_stack_t S, slong k)
fmpz_mod_poly_struct ** fmpz_mod_poly_stack_request(fmpz_mod_poly_stack_t S, slong k)
fmpz_mod_poly_struct * fmpz_mod_poly_stack_take_top(fmpz_mod_poly_stack_t S)
void fmpz_mod_poly_stack_give_back(fmpz_mod_poly_stack_t S, slong k)
slong fmpz_mod_poly_stack_size(const fmpz_mod_poly_stack_t S)
void fmpz_mod_bpoly_stack_init(fmpz_mod_bpoly_stack_t S)
void fmpz_mod_bpoly_stack_clear(fmpz_mod_bpoly_stack_t S)
fmpz_mod_bpoly_struct ** fmpz_mod_bpoly_stack_fit_request(fmpz_mod_bpoly_stack_t S, slong k)
fmpz_mod_bpoly_struct ** fmpz_mod_bpoly_stack_request(fmpz_mod_bpoly_stack_t S, slong k)
fmpz_mod_bpoly_struct * fmpz_mod_bpoly_stack_take_top(fmpz_mod_bpoly_stack_t S)
void fmpz_mod_bpoly_stack_give_back(fmpz_mod_bpoly_stack_t S, slong k)
slong fmpz_mod_bpoly_stack_size(const fmpz_mod_bpoly_stack_t S)
void fmpz_mod_polyun_stack_init(fmpz_mod_polyun_stack_t S)
void fmpz_mod_polyun_stack_clear(fmpz_mod_polyun_stack_t S)
fmpz_mod_polyun_struct ** fmpz_mod_polyun_stack_fit_request(fmpz_mod_polyun_stack_t S, slong k)
fmpz_mod_polyun_struct ** fmpz_mod_polyun_stack_request(fmpz_mod_polyun_stack_t S, slong k)
fmpz_mod_polyun_struct * fmpz_mod_polyun_stack_take_top(fmpz_mod_polyun_stack_t S)
void fmpz_mod_polyun_stack_give_back(fmpz_mod_polyun_stack_t S, slong k)
slong fmpz_mod_polyun_stack_size(const fmpz_mod_polyun_stack_t S)
void fmpz_mod_mpolyn_stack_init(fmpz_mod_mpolyn_stack_t S, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_stack_clear(fmpz_mod_mpolyn_stack_t S, const fmpz_mod_mpoly_ctx_t ctx)
fmpz_mod_mpolyn_struct ** fmpz_mod_mpolyn_stack_fit_request(fmpz_mod_mpolyn_stack_t S, slong k, const fmpz_mod_mpoly_ctx_t ctx)
fmpz_mod_mpolyn_struct ** fmpz_mod_mpolyn_stack_request(fmpz_mod_mpolyn_stack_t S, slong k, const fmpz_mod_mpoly_ctx_t ctx)
fmpz_mod_mpolyn_struct * fmpz_mod_mpolyn_stack_take_top(fmpz_mod_mpolyn_stack_t S)
void fmpz_mod_mpolyn_stack_give_back(fmpz_mod_mpolyn_stack_t S, slong k)
slong fmpz_mod_mpolyn_stack_size(const fmpz_mod_mpolyn_stack_t S)

Initialise and clear a stack; `fit_request` ensures that `k` more objects
are available and returns a pointer to them, `request` also reserves
them, `take_top` reserves one fitted object, `give_back` releases the `k`
most recently reserved objects, and `size` returns the number of reserved
objects.
:::

### Univariate helpers

::: function
slong _fmpz_mod_poly_degree(const fmpz_mod_poly_t a)
void fmpz_mod_poly_scalar_addmul_fmpz_mod(fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_poly_t C, const fmpz_t d0, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_addmul_linear(fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_poly_t C, const fmpz_t d1, const fmpz_t d0, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(fmpz_mod_poly_t A, slong k, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_eval_pow(fmpz_t eval, const fmpz_mod_poly_t P, fmpz_mod_poly_t alphapow, const fmpz_mod_ctx_t ctx)
void fmpz_mod_poly_eval2_pow(fmpz_t evalp, fmpz_t evalm, const fmpz_mod_poly_t P, fmpz_mod_poly_t alphapow, const fmpz_mod_ctx_t ctx)

The degree of `a`; $A = B + d_0 C$; $A = B + (d_1 x + d_0) C$; $A = x^k A
+ c$; evaluation of `P` at $\alpha$, and at $\pm\alpha$ simultaneously,
using the cached powers `alphapow` of $\alpha$.
:::

::: function
slong _fmpz_mod_poly_vec_max_degree(const fmpz_mod_poly_struct * A, slong Alen, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_vec_content(fmpz_mod_poly_t g, const fmpz_mod_poly_struct * A, slong Alen, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_vec_remove_content(fmpz_mod_poly_t g, fmpz_mod_poly_struct * A, slong Alen, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_vec_mul_poly(fmpz_mod_poly_struct * A, slong Alen, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_vec_divexact_poly(fmpz_mod_poly_struct * A, slong Alen, const fmpz_mod_poly_t g, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_poly_vec_mul_fmpz_mod(fmpz_mod_poly_struct * A, slong Alen, const fmpz_t g, const fmpz_mod_ctx_t ctx)

Operations on vectors of univariate polynomials: maximum degree, monic
GCD of the entries (and division by it), multiplication and exact division
by `g`, and scaling.
:::

::: function
void fmpz_mod_pow_cache_start(const fmpz_t b, fmpz_mod_poly_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_pow_cache_mulpow_ui(fmpz_t a, const fmpz_t b, ulong e, fmpz_mod_poly_t c, const fmpz_mod_ctx_t ctx)

Start a cache `c` of powers of `b`, and set $a = a \cdot b^e$ using and
extending the cache.
:::

::: function
int fmpz_mod_mat_is_reduced(const fmpz_mod_mat_t N, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mat_init_nullspace_tr(fmpz_mod_mat_t X, fmpz_mod_mat_t tmp, const fmpz_mod_ctx_t ctx)

As `nmod_mat_is_reduced` and `nmod_mat_init_nullspace_tr`: test whether
the rows of `N` describe a partition of the columns, and compute a basis
of the right nullspace of `tmp` as the rows of `X`.
:::

### Bivariate and sparse polynomials

::: function
void fmpz_mod_bpoly_init(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_clear(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_swap(fmpz_mod_bpoly_t A, fmpz_mod_bpoly_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_fit_length(fmpz_mod_bpoly_t A, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_normalise(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
int fmpz_mod_bpoly_is_canonical(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_zero(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_one(fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
int fmpz_mod_bpoly_is_one(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_set(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_ctx_t ctx)
int fmpz_mod_bpoly_equal(const fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_get_coeff(fmpz_t c, const fmpz_mod_bpoly_t A, slong e0, slong e1, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_set_coeff(fmpz_mod_bpoly_t A, slong xi, slong yi, const fmpz_t c, const fmpz_mod_ctx_t ctx)
slong fmpz_mod_bpoly_degree0(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
slong fmpz_mod_bpoly_degree1(const fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_set_poly_gen0(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_set_poly_gen1(fmpz_mod_bpoly_t A, const fmpz_mod_poly_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_print_pretty(const fmpz_mod_bpoly_t A, const char * xvar, const char * yvar, const fmpz_mod_ctx_t ctx)

Memory management, assignment, coefficient access, degrees in the main
(`gen0`) and the inner (`gen1`) variable, conversion from univariate
polynomials in either variable, and printing of dense bivariate
polynomials.
:::

::: function
void fmpz_mod_bpoly_add(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_sub(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_mul(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_mul_series(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_bpoly_t C, slong order, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_divrem_series(fmpz_mod_bpoly_t Q, fmpz_mod_bpoly_t R, const fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, slong order, const fmpz_mod_ctx_t ctx)
int fmpz_mod_bpoly_divides(fmpz_mod_bpoly_t Q, const fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_make_monic_series(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, slong order, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_make_primitive(fmpz_mod_poly_t g, fmpz_mod_bpoly_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_derivative_gen0(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_taylor_shift_gen0(fmpz_mod_bpoly_t A, const fmpz_t alpha, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_taylor_shift_gen1(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_t c, const fmpz_mod_ctx_t ctx)
void fmpz_mod_bpoly_reverse_vars(fmpz_mod_bpoly_t A, const fmpz_mod_bpoly_t B, const fmpz_mod_ctx_t ctx)

Arithmetic of bivariate polynomials: products, and products and division
with the inner variable truncated to `order` terms; exact division
(returning 1 if exact); making the polynomial monic as a series in the
inner variable, or removing its content `g` (a polynomial in the inner
variable); the derivative with respect to the main variable; the
substitutions $x \mapsto x + \alpha$ and $y \mapsto y + c$; and swapping
the variables.
:::

::: function
void fmpz_mod_tpoly_init(fmpz_mod_tpoly_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_tpoly_clear(fmpz_mod_tpoly_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_tpoly_swap(fmpz_mod_tpoly_t A, fmpz_mod_tpoly_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_tpoly_fit_length(fmpz_mod_tpoly_t A, slong len, const fmpz_mod_ctx_t ctx)

Memory management of arrays of bivariate polynomials.
:::

::: function
void fmpz_mod_mpoly_get_fmpz_mod_bpoly(fmpz_mod_bpoly_t A, const fmpz_mod_mpoly_t B, slong var0, slong var1, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_fmpz_mod_bpoly(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mod_bpoly_t B, slong var0, slong var1, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_set_fmpz_mod_bpoly_var1_zero(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mod_bpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_get_polyu1n(fmpz_mod_polyun_t A, const fmpz_mod_mpoly_t B, slong varx, slong vary, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_set_polyu1n(fmpz_mod_mpoly_t B, const fmpz_mod_polyun_t A, slong varx, slong vary, const fmpz_mod_mpoly_ctx_t ctx)

Conversions between polynomials in two variables of the context and dense
or sparse bivariate representations.
:::

::: function
void fmpz_mod_polyu_init(fmpz_mod_polyu_t A)
void fmpz_mod_polyu_clear(fmpz_mod_polyu_t A)
void fmpz_mod_polyu_swap(fmpz_mod_polyu_t A, fmpz_mod_polyu_t B)
void fmpz_mod_polyu_realloc(fmpz_mod_polyu_t A, slong len)
void fmpz_mod_polyu_fit_length(fmpz_mod_polyu_t a, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyu3_degrees(slong * deg0, slong * deg1, slong * deg2, const fmpz_mod_polyu_t A)
void fmpz_mod_polyu3_print_pretty(const fmpz_mod_polyu_t A, const char * var0, const char * var1, const char * var2, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyun_init(fmpz_mod_polyun_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyun_clear(fmpz_mod_polyun_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyun_swap(fmpz_mod_polyun_t A, fmpz_mod_polyun_t B)
void fmpz_mod_polyun_realloc(fmpz_mod_polyun_t A, slong len, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyun_fit_length(fmpz_mod_polyun_t A, slong len, const fmpz_mod_ctx_t ctx)
int fmpz_mod_polyun_is_canonical(const fmpz_mod_polyun_t A, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyun_set(fmpz_mod_polyun_t A, const fmpz_mod_polyun_t B, const fmpz_mod_ctx_t ctx)
int fmpz_mod_polyun_equal(fmpz_mod_polyun_t A, const fmpz_mod_polyun_t B, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyun_one(fmpz_mod_polyun_t A, const fmpz_mod_ctx_t ctx)
const fmpz * fmpz_mod_polyun_leadcoeff(const fmpz_mod_polyun_t A)
ulong fmpz_mod_polyu1n_bidegree(const fmpz_mod_polyun_t A)
void fmpz_mod_polyu1n_print_pretty(const fmpz_mod_polyun_t A, const char * var0, const char * varlast, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyu2n_print_pretty(const fmpz_mod_polyun_t A, const char * var0, const char * var1, const char * varlast, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyu3n_print_pretty(const fmpz_mod_polyun_t A, const char * var0, const char * var1, const char * var2, const char * varlast, const fmpz_mod_ctx_t ctx)

Memory management, assignment, degrees and printing of the sparse types
with one, two or three packed main variables.
:::

::: function
void fmpz_mod_mpolyn_init(fmpz_mod_mpolyn_t A, flint_bitcnt_t bits, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_clear(fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_swap(fmpz_mod_mpolyn_t A, fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_fit_length(fmpz_mod_mpolyn_t A, slong length, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_one(fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_set(fmpz_mod_mpolyn_t A, const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyn_equal(const fmpz_mod_mpolyn_t A, const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyn_is_canonical(const fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyn_is_nonzero_fmpz(const fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_ctx_t ctx)
const fmpz * fmpz_mod_mpolyn_leadcoeff(const fmpz_mod_mpolyn_t A)
slong fmpz_mod_mpolyn_lastdeg(const fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_scalar_mul_fmpz_mod(fmpz_mod_mpolyn_t A, const fmpz_t c, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyn_divides(fmpz_mod_mpolyn_t Q, const fmpz_mod_mpolyn_t A, const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_print_pretty(const fmpz_mod_mpolyn_t poly, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_cvtto_mpolyn(fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_cvtfrom_mpolyn(fmpz_mod_mpoly_t A, const fmpz_mod_mpolyn_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_to_mpolyn_perm_deflate(fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_ctx_t nctx, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride)
void fmpz_mod_mpoly_from_mpolyn_perm_inflate(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t ctx, const fmpz_mod_mpolyn_t B, const fmpz_mod_mpoly_ctx_t nctx, const slong * perm, const ulong * shift, const ulong * stride)

Operations on `fmpz_mod_mpolyn_t` and conversions from and to
`fmpz_mod_mpoly_t` (with the variable `var` moved into the univariate
coefficients, or with a permutation and deflation of the variables).
:::

### Interpolation

::: function
void fmpz_mod_polyu1n_interp_reduce_2sm_poly(fmpz_mod_poly_t E, fmpz_mod_poly_t F, const fmpz_mod_polyun_t A, fmpz_mod_poly_t alphapow, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyu1n_interp_lift_2sm_poly(slong * lastdeg, fmpz_mod_polyun_t F, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_t alpha, const fmpz_mod_ctx_t ctx)
int fmpz_mod_polyu1n_interp_crt_2sm_poly(slong * lastdeg, fmpz_mod_polyun_t F, fmpz_mod_polyun_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t B, const fmpz_mod_poly_t modulus, fmpz_mod_poly_t alphapow, const fmpz_mod_ctx_t ctx)
void fmpz_mod_mpolyn_interp_reduce_sm_poly(fmpz_mod_poly_t E, const fmpz_mod_mpolyn_t A, const fmpz_t alpha, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_interp_lift_sm_poly(fmpz_mod_mpolyn_t A, const fmpz_mod_poly_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyn_interp_crt_sm_poly(slong * lastdeg_, fmpz_mod_mpolyn_t F, fmpz_mod_mpolyn_t T, const fmpz_mod_poly_t A, const fmpz_mod_poly_t modulus, const fmpz_t alpha, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_interp_reduce_2sm_mpolyn(fmpz_mod_mpolyn_t E, fmpz_mod_mpolyn_t F, fmpz_mod_mpolyn_t A, slong var, fmpz_mod_poly_t alphapow, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_interp_lift_2sm_mpolyn(slong * lastdeg, fmpz_mod_mpolyn_t T, fmpz_mod_mpolyn_t A, fmpz_mod_mpolyn_t B, slong var, const fmpz_t alpha, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyn_interp_crt_2sm_mpolyn(slong * lastdeg, fmpz_mod_mpolyn_t F, fmpz_mod_mpolyn_t T, fmpz_mod_mpolyn_t A, fmpz_mod_mpolyn_t B, slong var, fmpz_mod_poly_t modulus, fmpz_mod_poly_t alphapow, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyn_interp_lift_sm_mpoly(fmpz_mod_mpolyn_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyn_interp_crt_sm_mpoly(slong * lastdeg, fmpz_mod_mpolyn_t F, fmpz_mod_mpolyn_t T, fmpz_mod_mpoly_t A, fmpz_mod_poly_t modulus, fmpz_mod_poly_t alphapow, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(slong * lastdeg, fmpz_mod_mpolyn_t F, const fmpz_mod_mpoly_t A, const fmpz_mod_poly_t modulus, fmpz_mod_poly_t alphapow, const fmpz_mod_mpoly_ctx_t ctx)

The evaluation (`reduce`, at $\alpha$ or at $\pm\alpha$), first-image
(`lift`) and Chinese remaindering (`crt`) steps of Brown's algorithm, as
for `nmod_mpoly`; the `crt` functions return whether the interpolant
changed and update `lastdeg`.
:::

::: function
int _fmpz_mod_zip_vand_solve(fmpz * coeffs, const fmpz * monomials, slong mlength, const fmpz * evals, slong elength, const fmpz * master, fmpz * scratch, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_zip_eval_step(fmpz_t ev, fmpz * cur, const fmpz * inc, const fmpz * coeffs, slong length, const fmpz_mod_ctx_t ctx)
slong fmpz_mod_polyun_product_roots(fmpz_mod_polyun_t M, const fmpz_mod_polyun_t H, const fmpz_mod_ctx_t ctx)
void fmpz_mod_polyun_zip_start(fmpz_mod_polyun_t Z, fmpz_mod_polyun_t H, slong req_images, const fmpz_mod_ctx_t fctx)
void fmpz_mod_polyu2n_zip_eval_cur_inc_coeff(fmpz_mod_polyun_t E, fmpz_mod_polyun_t Acur, const fmpz_mod_polyun_t Ainc, const fmpz_mod_polyun_t Acoeff, const fmpz_mod_ctx_t ctx)
int fmpz_mod_polyun_zip_solve(fmpz_mod_mpoly_t A, fmpz_mod_polyun_t Z, fmpz_mod_polyun_t H, fmpz_mod_polyun_t M, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_mock_eval_coeff(fmpz_mod_polyun_t mock, const fmpz_mod_mpoly_t A, const fmpz_mod_polyun_t Aeh_inc, const fmpz_mod_mpoly_ctx_t ctx)
void mpoly_monomial_evals_fmpz_mod(fmpz_mod_poly_t EH, const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, fmpz_mod_poly_struct * alpha_caches, slong start, slong stop, const mpoly_ctx_t mctx, const fmpz_mod_ctx_t fpctx)
void mpoly1_monomial_evals_fmpz_mod(fmpz_mod_polyun_t EH, const ulong * Aexps, flint_bitcnt_t Abits, const ulong * Amarks, slong Amarkslen, fmpz_mod_poly_struct * alpha_caches, slong m, const mpoly_ctx_t mctx, const fmpz_mod_ctx_t fpctx)
void mpoly2_monomial_evals_fmpz_mod(fmpz_mod_polyun_t EH, const ulong * Aexps, flint_bitcnt_t Abits, ulong * Amarks, slong Amarkslen, fmpz_mod_poly_struct * alpha_caches, slong m, const mpoly_ctx_t mctx, const fmpz_mod_ctx_t fpctx)

Sparse interpolation of polynomials with known support: the monomial
evaluations at powers of the evaluation points, the master polynomials
$\prod (z - m_i)$ of the monomial values (`product_roots`, returning the
maximum number of terms), the incremental evaluation of the next image,
and the solution of the transposed Vandermonde systems (returning 1 for
success and 0 or $-1$ if the images are inconsistent).
:::

### GCD algorithms

::: function
int _fmpz_mod_mpoly_gcd_algo(fmpz_mod_mpoly_t G, fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx, unsigned int algo)

Computes the GCD `G` and, if `Abar` and `Bbar` are not `NULL`, the
cofactors of the nonzero `A` and `B` with the algorithms in the bit mask
`algo` (`MPOLY_GCD_USE_BROWN`, `MPOLY_GCD_USE_HENSEL`,
`MPOLY_GCD_USE_ZIPPEL`, `MPOLY_GCD_USE_ZIPPEL2`, `MPOLY_GCD_USE_PRS` or
`MPOLY_GCD_USE_ALL`). For a word-size modulus it calls
`_nmod_mpoly_gcd_algo_small`. Returns 1 for success.
:::

::: function
int fmpz_mod_polyu1n_gcd_brown_smprime(fmpz_mod_polyun_t G, fmpz_mod_polyun_t Abar, fmpz_mod_polyun_t Bbar, fmpz_mod_polyun_t A, fmpz_mod_polyun_t B, const fmpz_mod_ctx_t ctx, fmpz_mod_poly_stack_t St_poly, fmpz_mod_polyun_stack_t St_polyun)
int fmpz_mod_mpolyn_gcd_brown_smprime(fmpz_mod_mpolyn_t G, fmpz_mod_mpolyn_t Abar, fmpz_mod_mpolyn_t Bbar, fmpz_mod_mpolyn_t A, fmpz_mod_mpolyn_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx, const mpoly_gcd_info_t Iv, fmpz_mod_poly_polyun_mpolyn_stack_t St)
int fmpz_mod_mpolyl_gcdp_zippel(fmpz_mod_mpoly_t G, fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar, fmpz_mod_mpoly_t A, fmpz_mod_mpoly_t B, slong var, const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state)
int fmpz_mod_mpolyl_gcd_zippel2_smprime(fmpz_mod_mpoly_t rG, const slong * rGdegs, fmpz_mod_mpoly_t rAbar, fmpz_mod_mpoly_t rBbar, const fmpz_mod_mpoly_t A, const slong * Adegs, const fmpz_mod_mpoly_t B, const slong * Bdegs, const fmpz_mod_mpoly_t gamma, const slong * gammadegs, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpolyl_gcd_hensel_smprime(fmpz_mod_mpoly_t G, slong Gdeg, fmpz_mod_mpoly_t Abar, fmpz_mod_mpoly_t Bbar, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_ctx_t ctx)

The GCD algorithms with evaluation points in $\mathbb{F}_p$: Brown's
dense algorithm (bivariate and multivariate), Zippel's algorithm in the
variables up to `var`, sparse interpolation with two main variables given
the leading coefficient of the GCD, and Hensel lifting. Return 1 for
success, also setting the cofactors.
:::

### Factoring algorithms

::: function
void _fmpz_mod_mpoly_get_lead0(fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_set_lead0(fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx)

Get, or replace by `c`, the leading coefficient of a polynomial viewed as
a polynomial in variable 0.
:::

::: function
void _fmpz_mod_mpoly_factor_set_nmod_mpoly_factor(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_ctx_t ctx, const nmod_mpoly_factor_t nf, const nmod_mpoly_ctx_t nctx)

Sets `f` to the factorisation `nf` computed in the `nmod_mpoly` context
`nctx` with the same modulus.
:::

::: function
void fmpz_mod_mpolyv_init(fmpz_mod_mpolyv_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyv_clear(fmpz_mod_mpolyv_t A, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyv_swap(fmpz_mod_mpolyv_t A, fmpz_mod_mpolyv_t B, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyv_fit_length(fmpz_mod_mpolyv_t A, slong length, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyv_set_coeff(fmpz_mod_mpolyv_t A, slong i, fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpolyv_print_pretty(const fmpz_mod_mpolyv_t poly, const char ** x, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_to_mpolyv(fmpz_mod_mpolyv_t A, const fmpz_mod_mpoly_t B, const fmpz_mod_mpoly_t xalpha, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_from_mpolyv(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mod_mpolyv_t B, const fmpz_mod_mpoly_t xalpha, const fmpz_mod_mpoly_ctx_t ctx)

Vectors of polynomials, and conversion between `B` and the coefficients
of its expansion in powers of `xalpha` $= x_j - \alpha$.
:::

::: function
int _fmpz_mod_mpoly_vec_content_mpoly(fmpz_mod_mpoly_t g, const fmpz_mod_mpoly_struct * A, slong Alen, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_vec_divexact_mpoly(fmpz_mod_mpoly_struct * A, slong Alen, const fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx)
void _fmpz_mod_mpoly_vec_mul_mpoly(fmpz_mod_mpoly_struct * A, slong Alen, const fmpz_mod_mpoly_t c, const fmpz_mod_mpoly_ctx_t ctx)

GCD of `Alen` polynomials (returning 1 for success), and exact division
or multiplication of all of them by `c`.
:::

::: function
void fmpz_mod_mpoly_compression_do(fmpz_mod_mpoly_t L, const fmpz_mod_mpoly_ctx_t Lctx, fmpz * Acoeffs, slong Alen, mpoly_compression_t M)
void fmpz_mod_mpoly_compression_undo(fmpz_mod_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mod_mpoly_ctx_t Actx, fmpz_mod_mpoly_t L, const fmpz_mod_mpoly_ctx_t Lctx, mpoly_compression_t M)

Apply and undo a monomial compression computed by
`mpoly_compression_set`, as `nmod_mpoly_compression_do` and
`nmod_mpoly_compression_undo`.
:::

::: function
int _fmpz_mod_mpoly_factor_separable(fmpz_mod_mpoly_factor_t f, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, int sep)

As `fmpz_mod_mpoly_factor_separable` for a monic primitive `A`.
:::

::: function
int fmpz_mod_bpoly_factor_smprime(fmpz_mod_poly_t c, fmpz_mod_tpoly_t F, fmpz_mod_bpoly_t B, int allow_shift, const fmpz_mod_ctx_t ctx)

Factors the bivariate `B`, separable in the main variable: sets `c` to its
content and `F` to its irreducible factors, by evaluation of the inner
variable at a point of $\mathbb{F}_p$ (only 0 unless `allow_shift`),
Hensel lifting and recombination by linear algebra [@Lec2006]. Returns 0
if no suitable point is found.
:::

::: function
int fmpz_mod_bpoly_hlift2(fmpz_mod_bpoly_t A, fmpz_mod_bpoly_t B0, fmpz_mod_bpoly_t B1, const fmpz_t alpha, slong degree_inner, const fmpz_mod_ctx_t ctx, fmpz_mod_poly_bpoly_stack_t St)
int fmpz_mod_bpoly_hlift(slong r, fmpz_mod_bpoly_t A, fmpz_mod_bpoly_struct * B, const fmpz_t alpha, slong degree_inner, const fmpz_mod_ctx_t ctx, fmpz_mod_poly_bpoly_stack_t St)
int fmpz_mod_polyu3_hlift(slong r, fmpz_mod_polyun_struct * BB, fmpz_mod_polyu_t A, fmpz_mod_polyu_struct * B, const fmpz_t beta, slong degree_inner, const fmpz_mod_ctx_t ctx)

Bivariate (and trivariate) Hensel lifting of two or
`r` factors modulo powers of $y - \alpha$, as the `n_bpoly_mod` versions
in `nmod_mpoly_factor`.
:::

::: function
int fmpz_mod_mpoly_pfrac_init(fmpz_mod_mpoly_pfrac_t Iv, flint_bitcnt_t bits, slong l, slong r, const fmpz_mod_mpoly_struct * betas, const fmpz * alpha, const fmpz_mod_mpoly_ctx_t ctx)
void fmpz_mod_mpoly_pfrac_clear(fmpz_mod_mpoly_pfrac_t Iv, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_pfrac(slong r, fmpz_mod_mpoly_t t, const slong * deg, fmpz_mod_mpoly_pfrac_t Iv, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_hlift(slong m, fmpz_mod_mpoly_struct * f, slong r, const fmpz * alpha, const fmpz_mod_mpoly_t A, const slong * degs, const fmpz_mod_mpoly_ctx_t ctx)
int fmpz_mod_mpoly_hlift_zippel(slong m, fmpz_mod_mpoly_struct * B, slong r, const fmpz * alpha, const fmpz_mod_mpoly_t A, const slong * degs, const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state)

Multivariate partial fractions and Hensel lifting, dense or by sparse
interpolation, as the corresponding `nmod_mpoly` functions.
:::

::: function
int fmpz_mod_mpoly_factor_lcc_wang(fmpz_mod_mpoly_struct * lc_divs, const fmpz_mod_mpoly_factor_t lcAfac, const fmpz_mod_poly_t Auc, const fmpz_mod_bpoly_struct * Auf, slong r, const fmpz_mod_poly_struct * alpha, const fmpz_mod_mpoly_ctx_t ctx)

Wang's leading coefficient distribution [@Wan1978], as
`nmod_mpoly_factor_lcc_wang`.
:::

::: function
int fmpz_mod_mpoly_factor_irred_smprime_zassenhaus(fmpz_mod_mpolyv_t fac, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state)
int fmpz_mod_mpoly_factor_irred_smprime_wang(fmpz_mod_mpolyv_t fac, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_factor_t lcAfac, const fmpz_mod_mpoly_t lcA, const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state)
int fmpz_mod_mpoly_factor_irred_smprime_zippel(fmpz_mod_mpolyv_t fac, const fmpz_mod_mpoly_t A, const fmpz_mod_mpoly_factor_t lcAfac, const fmpz_mod_mpoly_t lcA, const fmpz_mod_mpoly_ctx_t ctx, flint_rand_t state)

Irreducible factorisation of a monic squarefree polynomial in at least
three variables with evaluation points in $\mathbb{F}_p$. Return 1 for
success, 0 if they ran out of evaluation points and $-1$ for a failure
that other points would not fix.
:::

::: function
int _fmpz_mod_mpoly_evaluate_rest_fmpz_mod_poly(fmpz_mod_poly_struct * E, slong * starts, slong * ends, slong * stops, ulong * es, const fmpz * Acoeffs, const ulong * Aexps, slong Alen, slong var, const fmpz_mod_poly_struct * alphas, const slong * offsets, const slong * shifts, slong N, ulong mask, slong nvars, const fmpz_mod_ctx_t ctx)
void _fmpz_mod_mpoly_eval_rest_to_fmpz_mod_bpoly(fmpz_mod_bpoly_t E, const fmpz_mod_mpoly_t A, const fmpz_mod_poly_struct * alphabetas, const fmpz_mod_mpoly_ctx_t ctx)

Evaluate all but the first one or two variables, giving univariate or
bivariate images.
:::
