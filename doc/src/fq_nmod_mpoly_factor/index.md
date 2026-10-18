---
title: fq_nmod_mpoly_factor.h
header: fq_nmod_mpoly_factor.h
description: Factorisation of multivariate polynomials over finite fields $\mathbb{F}_q$ of word-size characteristic, with extension-field methods for small fields.
category: Multivariate polynomials
order: 13
---

This module factors polynomials of `fq_nmod_mpoly.h` over $\mathbb{F}_q$,
$q = p^d$. A factorisation is stored in an `fq_nmod_mpoly_factor_t`, a
constant in $\mathbb{F}_q$ times a product of powers of polynomials, with
the same interface as `nmod_mpoly_factor` (see [](nmod_mpoly_factor)).

The algorithms are those of `nmod_mpoly_factor`, adapted to
$\mathbb{F}_q$. In the squarefree decomposition, a polynomial all of
whose exponents are divisible by $p^k$ is a $p^k$-th power; its root is
obtained by dividing the exponents by $p^k$ and applying the inverse of
the Frobenius automorphism $c \mapsto c^{p^k}$ to the coefficients, which
is $c \mapsto c^{p^{j}}$ with $j \equiv -k \pmod d$. The irreducible
factorisation evaluates at points of $\mathbb{F}_q$ when there are enough
of them (`smprime`), and otherwise over extensions $\mathbb{F}_{q^n}$
presented through the embeddings of `fq_nmod_mpoly.h` (`lgprime`),
combining the factors found there into factors over $\mathbb{F}_q$. See
[](algorithms/multivariate/modp-factoring) and
[](algorithms/multivariate/hensel-lifting).

## Types and memory management {#types}

::: type
fq_nmod_mpoly_factor_struct
fq_nmod_mpoly_factor_t

A factorisation $c \prod_i f_i^{e_i}$, with the fields `constant` (an
`fq_nmod_t`), `poly` (the bases), `exp` (the `fmpz` exponents), `num` and
`alloc`.
:::

::: function
void fq_nmod_mpoly_factor_init(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_factor_clear(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)

Initialise `f` as the empty product with constant 1, and free it.
:::

::: function
void fq_nmod_mpoly_factor_realloc(fq_nmod_mpoly_factor_t f, slong alloc, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_factor_fit_length(fq_nmod_mpoly_factor_t f, slong len, const fq_nmod_mpoly_ctx_t ctx)

Change the allocation of `f` to `alloc` factors, or ensure space for `len`
factors.
:::

## Basic manipulation {#basic}

::: function
void fq_nmod_mpoly_factor_swap(fq_nmod_mpoly_factor_t A, fq_nmod_mpoly_factor_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_factor_set(fq_nmod_mpoly_factor_t a, const fq_nmod_mpoly_factor_t b, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_factor_one(fq_nmod_mpoly_factor_t a, const fq_nmod_mpoly_ctx_t ctx)

Swap two factorisations, copy `b` into `a`, and set `a` to the empty
product with constant 1.
:::

::: function
slong fq_nmod_mpoly_factor_length(const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_factor_get_constant_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)

Return the number of factors, and set `c` to the constant.
:::

::: function
void fq_nmod_mpoly_factor_get_base(fq_nmod_mpoly_t p, const fq_nmod_mpoly_factor_t f, slong i, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_factor_swap_base(fq_nmod_mpoly_t p, const fq_nmod_mpoly_factor_t f, slong i, const fq_nmod_mpoly_ctx_t ctx)
slong fq_nmod_mpoly_factor_get_exp_si(fq_nmod_mpoly_factor_t f, slong i, const fq_nmod_mpoly_ctx_t ctx)

Set `p` to the base of index $i$ or swap it with `p` (despite the `const`
in the prototype, `fq_nmod_mpoly_factor_swap_base` modifies `f`), and
return the exponent of index $i$, which must fit an `slong`.
:::

::: function
void fq_nmod_mpoly_factor_append_ui(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, ulong e, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_factor_append_fmpz(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fmpz_t e, const fq_nmod_mpoly_ctx_t ctx)

Append the factor $A^e$ to `f`, copying `A`.
:::

::: function
int fq_nmod_mpoly_factor_expand(fq_nmod_mpoly_t A, const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_factor_matches(const fq_nmod_mpoly_t a, const fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)

Set `A` to the product represented by `f` (returning 0 for negative or
too large exponents), and test whether `f` represents `a`.
:::

::: function
void fq_nmod_mpoly_factor_sort(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_factor_cmp(const fq_nmod_mpoly_factor_t A, const fq_nmod_mpoly_factor_t B, const fq_nmod_mpoly_ctx_t ctx)

Sort the factors by exponent and then by base, and compare two
factorisations (constant, number of factors, then factor by factor).
:::

::: function
void fq_nmod_mpoly_factor_print_pretty(const fq_nmod_mpoly_factor_t f, const char ** vars, const fq_nmod_mpoly_ctx_t ctx)

Prints the constant followed by one line `*(base)^exp` per factor.
:::

## Factorisation {#factorisation}

These functions set `f` to a factorisation of `A` whose constant is the
leading coefficient of `A` and whose bases are monic, and return 1 for
success (0 only when exponents do not fit one word per field after
removing the content, or an algorithm runs out of extension fields).

::: function
int fq_nmod_mpoly_factor_content(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Sets `f` to a factorisation into variables and polynomials primitive with
respect to every variable.
:::

::: function
int fq_nmod_mpoly_factor_squarefree(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_factor_separable(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx, int sep)

Set `f` to a factorisation with squarefree, pairwise coprime bases. With
`sep` = 1, `fq_nmod_mpoly_factor_separable` also makes each base $g$
satisfy $\partial g/\partial x_i = 0$ or $\gcd(g, \partial g/\partial
x_i) = 1$ for every variable, as `nmod_mpoly_factor_separable`.
:::

::: function
int fq_nmod_mpoly_factor(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)

Sets `f` to the factorisation of `A` into monic irreducible factors over
$\mathbb{F}_q$.
:::

::: function
int fq_nmod_mpoly_factor_zassenhaus(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_factor_wang(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_factor_zippel(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_factor_algo(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx, unsigned int algo)

As `fq_nmod_mpoly_factor` with only one lifting method, or the methods in
the bit mask `algo` (`MPOLY_FACTOR_USE_ZAS`, `MPOLY_FACTOR_USE_WANG`,
`MPOLY_FACTOR_USE_ZIP`, `MPOLY_FACTOR_USE_ALL`).
:::

## Internal functions {#internal}

These functions are used by the factoring algorithms. Their interfaces
may change without notice.

### Bivariate polynomials over $\mathbb{F}_q$

The `n_fq_bpoly` functions work on `n_bpoly_t` (see [](n_poly)) whose
coefficients are univariate polynomials with $d$ words per coefficient,
the `n_fq` format.

::: function
void n_fq_bpoly_add(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_sub(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_mul(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_mul_series(n_bpoly_t A, const n_bpoly_t B, const n_bpoly_t C, slong order, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_divrem_series(n_bpoly_t Q, n_bpoly_t R, const n_bpoly_t A, const n_bpoly_t B, slong order, const fq_nmod_ctx_t ctx)
int n_fq_bpoly_divides(n_bpoly_t Q, const n_bpoly_t A, const n_bpoly_t B, const fq_nmod_ctx_t ctx)
void n_fq_bpoly_make_primitive(n_poly_t g, n_bpoly_t A, const fq_nmod_ctx_t ctx)

Sum, difference and product; product and division with remainder with the
inner variable truncated to `order` terms; exact division (returning 1 if
exact); and removal of the content `g` (a polynomial in the inner
variable).
:::

::: function
void fq_nmod_mpoly_get_n_fq_bpoly(n_bpoly_t A, const fq_nmod_mpoly_t B, slong varx, slong vary, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_set_n_fq_bpoly(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const n_bpoly_t B, slong varx, slong vary, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_set_n_fq_bpoly_gen1_zero(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const n_bpoly_t B, slong var, const fq_nmod_mpoly_ctx_t ctx)

Convert between a polynomial in the variables `varx` (main) and `vary`
only and a dense bivariate polynomial.
:::

::: function
int n_fq_bpoly_factor_smprime(n_poly_t c, n_tpoly_t F, n_bpoly_t B, int allow_shift, const fq_nmod_ctx_t ctx)
int n_fq_bpoly_factor_lgprime(n_poly_t c, n_tpoly_t F, n_bpoly_t B, const fq_nmod_ctx_t ctx, flint_rand_t state)

Factor the bivariate `B`, separable in the main variable: set `c` to its
content and `F` to its irreducible factors. The `smprime` version
evaluates the inner variable at a point of $\mathbb{F}_q$ (only 0 unless
`allow_shift`), lifts by Hensel lifting and recombines by linear algebra
[@Lec2006], returning 0 if no point works; the `lgprime` version reduces
the inner variable modulo irreducible polynomials over $\mathbb{F}_q$
(working over $\mathbb{F}_{q^n}$) and always succeeds.
:::

::: function
int n_fq_bpoly_hlift2(n_bpoly_t A, n_bpoly_t B0, n_bpoly_t B1, const fq_nmod_t alpha, slong degree_inner, const fq_nmod_ctx_t ctx, n_poly_bpoly_stack_t St)
int n_fq_bpoly_hlift2_cubic(n_fq_bpoly_t A, n_fq_bpoly_t B0, n_fq_bpoly_t B1, const fq_nmod_t alpha_, slong degree_inner, const fq_nmod_ctx_t ctx, nmod_eval_interp_t E, n_poly_bpoly_stack_t St)
int n_fq_bpoly_hlift(slong r, n_bpoly_t A, n_bpoly_struct * B, const fq_nmod_t alpha, slong degree_inner, const fq_nmod_ctx_t ctx, n_poly_bpoly_stack_t St)
int n_fq_polyu3_hlift(slong r, n_polyun_struct * BB, n_polyu_t A, n_polyu_struct * B, const fq_nmod_t beta, slong degree_inner, const fq_nmod_ctx_t ctx, n_poly_bpoly_stack_t St)

Hensel lifting of a factorisation of a bivariate (trivariate for
`n_fq_polyu3_hlift`) polynomial `A` modulo $y - \alpha$ into two factors
`B0`, `B1` or `r` factors `B`, up to degree `degree_inner` in $y$; the
`cubic` version evaluates and interpolates in the other variable. Return 1
for success.
:::

::: function
int n_polyu_fq_is_canonical(const n_polyu_t A, const fq_nmod_ctx_t ctx)
int n_polyun_fq_is_canonical(const n_polyun_t A, const fq_nmod_ctx_t ctx)
void n_polyu3_fq_print_pretty(const n_polyu_t A, const char * var0, const char * var1, const char * var2, const fq_nmod_ctx_t ctx)
void n_polyu2n_fq_print_pretty(const n_polyun_t A, const char * var0, const char * var1, const char * varlast, const fq_nmod_ctx_t ctx)
void n_polyu3n_fq_print_pretty(const n_polyun_t A, const char * var0, const char * var1, const char * var2, const char * varlast, const fq_nmod_ctx_t ctx)

Consistency checks and printing of sparse polynomials with `n_fq`
coefficients.
:::

### Multivariate lifting

::: type
fq_nmod_mpolyv_struct
fq_nmod_mpolyv_t
fq_nmod_mpoly_pfrac_struct
fq_nmod_mpoly_pfrac_t

Vectors of polynomials (`coeffs`, `alloc`, `length`), and precomputed
data for multivariate partial fractions.
:::

::: function
void fq_nmod_mpolyv_init(fq_nmod_mpolyv_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyv_clear(fq_nmod_mpolyv_t A, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyv_swap(fq_nmod_mpolyv_t A, fq_nmod_mpolyv_t B, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyv_fit_length(fq_nmod_mpolyv_t A, slong length, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyv_set_coeff(fq_nmod_mpolyv_t A, slong i, fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpolyv_print_pretty(const fq_nmod_mpolyv_t poly, const char ** x, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_to_mpolyv(fq_nmod_mpolyv_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t xalpha, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_from_mpolyv(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpolyv_t B, const fq_nmod_mpoly_t xalpha, const fq_nmod_mpoly_ctx_t ctx)

Vectors of polynomials, and conversion between `B` and the coefficients
of its expansion in powers of `xalpha` $= x_j - \alpha$.
:::

::: function
int _fq_nmod_mpoly_vec_content_mpoly(fq_nmod_mpoly_t g, const fq_nmod_mpoly_struct * A, slong Alen, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_vec_divexact_mpoly(fq_nmod_mpoly_struct * A, slong Alen, const fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_vec_mul_mpoly(fq_nmod_mpoly_struct * A, slong Alen, const fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)

GCD of `Alen` polynomials (returning 1 for success), and exact division
or multiplication of all of them by `c`.
:::

::: function
void _fq_nmod_mpoly_get_lead0(fq_nmod_mpoly_t c, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
void _fq_nmod_mpoly_set_lead0(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B, const fq_nmod_mpoly_t c, const fq_nmod_mpoly_ctx_t ctx)

Get, or replace by `c`, the leading coefficient of a polynomial viewed as
a polynomial in variable 0.
:::

::: function
void fq_nmod_mpoly_compression_do(fq_nmod_mpoly_t L, const fq_nmod_mpoly_ctx_t Lctx, ulong * Acoeffs, slong Alen, mpoly_compression_t M)
void fq_nmod_mpoly_compression_undo(fq_nmod_mpoly_t A, flint_bitcnt_t Abits, const fq_nmod_mpoly_ctx_t Actx, fq_nmod_mpoly_t L, const fq_nmod_mpoly_ctx_t Lctx, mpoly_compression_t M)

Apply and undo a monomial compression computed by
`mpoly_compression_set`.
:::

::: function
int fq_nmod_mpoly_pfrac_init(fq_nmod_mpoly_pfrac_t Iv, flint_bitcnt_t bits, slong l, slong r, const fq_nmod_mpoly_struct * betas, const fq_nmod_struct * alpha, const fq_nmod_mpoly_ctx_t ctx)
void fq_nmod_mpoly_pfrac_clear(fq_nmod_mpoly_pfrac_t Iv, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_pfrac(slong r, fq_nmod_mpoly_t t, const slong * deg, fq_nmod_mpoly_pfrac_t Iv, const fq_nmod_mpoly_ctx_t ctx)
int fq_nmod_mpoly_hlift(slong m, fq_nmod_mpoly_struct * f, slong r, const fq_nmod_struct * alpha, const fq_nmod_mpoly_t A, const slong * degs, const fq_nmod_mpoly_ctx_t ctx)

Multivariate partial fractions and Hensel lifting, as the corresponding
`nmod_mpoly` functions.
:::

::: function
int fq_nmod_mpoly_factor_lcc_wang(fq_nmod_mpoly_struct * lc_divs, const fq_nmod_mpoly_factor_t lcAfac, const n_poly_t Auc, const n_bpoly_struct * Auf, slong r, const n_poly_struct * alpha, const fq_nmod_mpoly_ctx_t ctx)

Wang's leading coefficient distribution [@Wan1978], as
`nmod_mpoly_factor_lcc_wang`.
:::

::: function
int _fq_nmod_mpoly_factor_separable(fq_nmod_mpoly_factor_t f, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx, int sep)

As `fq_nmod_mpoly_factor_separable` for a monic primitive `A`.
:::

::: function
int fq_nmod_mpoly_factor_irred_smprime_zassenhaus(fq_nmod_mpolyv_t fac, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t state)
int fq_nmod_mpoly_factor_irred_smprime_wang(fq_nmod_mpolyv_t fac, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_factor_t lcAfac, const fq_nmod_mpoly_t lcA, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t state)
int fq_nmod_mpoly_factor_irred_smprime_zippel(fq_nmod_mpolyv_t fac, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_factor_t lcAfac, const fq_nmod_mpoly_t lcA, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t state)
int fq_nmod_mpoly_factor_irred_lgprime_zassenhaus(fq_nmod_mpolyv_t fac, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t state)
int fq_nmod_mpoly_factor_irred_lgprime_wang(fq_nmod_mpolyv_t Af, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_factor_t lcAfac, const fq_nmod_mpoly_t lcA, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t state)
int fq_nmod_mpoly_factor_irred_lgprime_zippel(fq_nmod_mpolyv_t Af, const fq_nmod_mpoly_t A, const fq_nmod_mpoly_factor_t lcAfac, const fq_nmod_mpoly_t lcA, const fq_nmod_mpoly_ctx_t ctx, flint_rand_t state)

Irreducible factorisation of a monic squarefree polynomial in at least
three variables, given for the Wang and Zippel methods the factorisation
`lcAfac` of the leading coefficient `lcA` in variable 0. The `smprime`
versions use evaluation points in $\mathbb{F}_q$; the `lgprime` versions
map `A` into $\mathbb{F}_{q^n}$ with the embeddings of `fq_nmod_mpoly.h`,
factor there and multiply together the factors in each orbit of the
Frobenius automorphism $c \mapsto c^q$ of $\mathbb{F}_{q^n}$ over
$\mathbb{F}_q$. Return 1 for
success, 0 if they ran out of evaluation points or extensions, and $-1$
for a failure that other points would not fix.
:::

::: function
int _fq_nmod_mpoly_eval_rest_n_fq_poly(n_poly_struct * E, slong * starts, slong * ends, slong * stops, ulong * es, const ulong * Acoeffs, const ulong * Aexps, slong Alen, slong var, const n_fq_poly_struct * alphas, const slong * offsets, const slong * shifts, slong N, ulong mask, slong nvars, const fq_nmod_ctx_t ctx)
void _fq_nmod_mpoly_eval_rest_to_n_fq_bpoly(n_bpoly_t E, const fq_nmod_mpoly_t A, const n_poly_struct * alphabetas, const fq_nmod_mpoly_ctx_t ctx)

Evaluate all but the first one or two variables, giving univariate or
bivariate images.
:::
