---
title: fq_zech_mpoly_factor.h
header: fq_zech_mpoly_factor.h
description: Factorisation of multivariate polynomials over small finite fields in Zech logarithm representation, and the bivariate and lifting machinery used by the `medprime` factoring algorithms of `nmod_mpoly`.
category: Multivariate polynomials
order: 15
---

This module factors polynomials of `fq_zech_mpoly.h`. The public
factorisation functions convert to `fq_nmod_mpoly` and call the
corresponding `fq_nmod_mpoly_factor` function, so their semantics are
those of [](fq_nmod_mpoly_factor). The rest of the module is a port of
the `smprime` factoring algorithms (bivariate factorisation, Wang's
leading coefficients, Hensel lifting, sparse lifting and Zassenhaus
recombination) to Zech logarithm arithmetic, which
`nmod_mpoly_factor_irred_medprime_*` uses to factor a polynomial over
$\mathbb{F}_p$ in an extension $\mathbb{F}_{p^k}$ with $p^k \le 10^6$
(see [](algorithms/multivariate/modp-factoring) and
[](algorithms/multivariate/hensel-lifting)).

## Types {#types}

::: type
fq_zech_mpoly_factor_struct
fq_zech_mpoly_factor_t

A factorisation $c \prod_i f_i^{e_i}$, with the fields `constant` (an
`fq_zech_t`), `poly`, `exp` (`fmpz` exponents), `num` and `alloc`.
:::

## Memory management and basic manipulation {#basic}

::: function
void fq_zech_mpoly_factor_init(fq_zech_mpoly_factor_t f, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_factor_clear(fq_zech_mpoly_factor_t f, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_factor_realloc(fq_zech_mpoly_factor_t f, slong alloc, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_factor_fit_length(fq_zech_mpoly_factor_t f, slong len, const fq_zech_mpoly_ctx_t ctx)

Initialise `f` as the empty product with constant 1, free it, change its
allocation to `alloc` factors, or ensure space for `len` factors.
:::

::: function
void fq_zech_mpoly_factor_swap(fq_zech_mpoly_factor_t A, fq_zech_mpoly_factor_t B, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_factor_one(fq_zech_mpoly_factor_t a, const fq_zech_mpoly_ctx_t ctx)

Swap two factorisations, and set `a` to the empty product with constant
1.
:::

::: function
void fq_zech_mpoly_factor_print_pretty(const fq_zech_mpoly_factor_t f, const char ** vars, const fq_zech_mpoly_ctx_t ctx)

Prints the constant followed by one line `*(base)^exp` per factor.
:::

## Factorisation {#factorisation}

::: function
int fq_zech_mpoly_factor(fq_zech_mpoly_factor_t f, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_factor_zassenhaus(fq_zech_mpoly_factor_t f, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_factor_wang(fq_zech_mpoly_factor_t f, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_factor_zippel(fq_zech_mpoly_factor_t f, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_factor_algo(fq_zech_mpoly_factor_t f, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx, unsigned int algo)

Set `f` to the factorisation of `A` into monic irreducible factors, with
the constant equal to the leading coefficient of `A`, using all methods,
one lifting method, or the methods in the bit mask `algo` (see
`fq_nmod_mpoly_factor_algo`). The computation is done by
`fq_nmod_mpoly_factor_algo` after conversion. Return 1 for success.
:::

## Internal functions {#internal}

### Univariate and bivariate polynomials

::: type
fq_zech_bpoly_struct
fq_zech_bpoly_t
fq_zech_tpoly_struct
fq_zech_tpoly_t
fq_zech_polyu_struct
fq_zech_polyu_t
fq_zech_polyun_struct
fq_zech_polyun_t

Dense bivariate polynomials (arrays of `fq_zech_poly_struct`, the entry
of index $i$ being the coefficient of $x^i$ as a polynomial in $y$),
arrays of bivariate polynomials, and sparse polynomials in up to three
packed variables with `fq_zech` coefficients (`polyu`) or with
`fq_zech_poly` coefficients in one more variable (`polyun`).
:::

::: function
void fq_zech_bpoly_init(fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_clear(fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_swap(fq_zech_bpoly_t A, fq_zech_bpoly_t B, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_realloc(fq_zech_bpoly_t A, slong len, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_fit_length(fq_zech_bpoly_t A, slong len, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_normalise(fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
int fq_zech_bpoly_is_canonical(const fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_zero(fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_one(fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
int fq_zech_bpoly_is_zero(const fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_set(fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_ctx_t ctx)
int fq_zech_bpoly_equal(const fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_get_coeff(fq_zech_t c, const fq_zech_bpoly_t A, slong e0, slong e1, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_set_coeff_fq_zech(fq_zech_bpoly_t A, slong xi, slong yi, const fq_zech_t c, const fq_zech_ctx_t ctx)
slong fq_zech_bpoly_degree0(const fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
slong fq_zech_bpoly_degree1(const fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_set_fq_zech_poly_var0(fq_zech_bpoly_t A, const fq_zech_poly_t B, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_set_fq_zech_poly_var1(fq_zech_bpoly_t A, const fq_zech_poly_t B, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_print_pretty(const fq_zech_bpoly_t A, const char * var0, const char * var1, const fq_zech_ctx_t ctx)

Memory management, assignment, coefficient access, the degrees in the
main (`0`) and inner (`1`) variable, and conversion from univariate
polynomials in either variable.
:::

::: function
void fq_zech_bpoly_add(fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_bpoly_t C, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_sub(fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_bpoly_t C, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_mul(fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_bpoly_t C, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_mul_series(fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_bpoly_t C, slong order, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_divrem_series(fq_zech_bpoly_t Q, fq_zech_bpoly_t R, const fq_zech_bpoly_t A, const fq_zech_bpoly_t B, slong order, const fq_zech_ctx_t ctx)
int fq_zech_bpoly_divides(fq_zech_bpoly_t Q, const fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_make_primitive(fq_zech_poly_t g, fq_zech_bpoly_t A, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_derivative(fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_taylor_shift_var0(fq_zech_bpoly_t A, const fq_zech_t alpha, const fq_zech_ctx_t ctx)
void fq_zech_bpoly_taylor_shift_var1(fq_zech_bpoly_t A, const fq_zech_bpoly_t B, const fq_zech_t c_, const fq_zech_ctx_t ctx)

Arithmetic of bivariate polynomials: products and divisions, possibly
truncated to `order` terms in the inner variable; exact division
(returning 1 if exact); removal of the content `g`; the derivative with
respect to the main variable; and the shifts $x \mapsto x + \alpha$ and $y
\mapsto y + c$.
:::

::: function
void fq_zech_tpoly_init(fq_zech_tpoly_t A, const fq_zech_ctx_t ctx)
void fq_zech_tpoly_clear(fq_zech_tpoly_t A, const fq_zech_ctx_t ctx)
void fq_zech_tpoly_swap(fq_zech_tpoly_t A, fq_zech_tpoly_t B, const fq_zech_ctx_t ctx)
void fq_zech_tpoly_fit_length(fq_zech_tpoly_t A, slong len, const fq_zech_ctx_t ctx)
void fq_zech_polyu_init(fq_zech_polyu_t A, const fq_zech_ctx_t ctx)
void fq_zech_polyu_clear(fq_zech_polyu_t A, const fq_zech_ctx_t ctx)
void fq_zech_polyu_swap(fq_zech_polyu_t A, fq_zech_polyu_t B, const fq_zech_ctx_t ctx)
void fq_zech_polyu_realloc(fq_zech_polyu_t A, slong len, const fq_zech_ctx_t ctx)
void fq_zech_polyu_fit_length(fq_zech_polyu_t A, slong len, const fq_zech_ctx_t ctx)
int fq_zech_polyu_is_canonical(const fq_zech_polyu_t A, const fq_zech_ctx_t ctx)
void fq_zech_polyu3_degrees(slong * deg0, slong * deg1, slong * deg2, const fq_zech_polyu_t A)
void fq_zech_polyu3_print_pretty(const fq_zech_polyu_t A, const char * var0, const char * var1, const char * var2, const fq_zech_ctx_t ctx)
void fq_zech_polyun_init(fq_zech_polyun_t A, const fq_zech_ctx_t ctx)
void fq_zech_polyun_clear(fq_zech_polyun_t A, const fq_zech_ctx_t ctx)
void fq_zech_polyun_swap(fq_zech_polyun_t A, fq_zech_polyun_t B, const fq_zech_ctx_t ctx)
void fq_zech_polyun_realloc(fq_zech_polyun_t A, slong len, const fq_zech_ctx_t ctx)
void fq_zech_polyun_fit_length(fq_zech_polyun_t A, slong len, const fq_zech_ctx_t ctx)
int fq_zech_polyun_is_canonical(const fq_zech_polyun_t A, const fq_zech_ctx_t ctx)
void fq_zech_polyu2n_print_pretty(const fq_zech_polyun_t A, const char * var0, const char * var1, const char * varlast, const fq_zech_ctx_t ctx)
void fq_zech_polyu3n_print_pretty(const fq_zech_polyun_t A, const char * var0, const char * var1, const char * var2, const char * varlast, const fq_zech_ctx_t ctx)

Memory management, consistency checks, degrees and printing of the other
types.
:::

::: function
void fq_zech_mpoly_get_fq_zech_bpoly(fq_zech_bpoly_t A, const fq_zech_mpoly_t B, slong varx, slong vary, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_set_fq_zech_bpoly(fq_zech_mpoly_t A, flint_bitcnt_t Abits, const fq_zech_bpoly_t B, slong varx, slong vary, const fq_zech_mpoly_ctx_t ctx)
void _fq_zech_mpoly_set_fq_zech_bpoly_var1_zero(fq_zech_mpoly_t A, flint_bitcnt_t Abits, const fq_zech_bpoly_t B, slong var, const fq_zech_mpoly_ctx_t ctx)

Convert between a polynomial in the variables `varx` (main) and `vary`
only and a dense bivariate polynomial.
:::

::: function
int fq_zech_mpoly_is_fq_zech_poly(const fq_zech_mpoly_t A, slong var, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_get_fq_zech_poly(fq_zech_poly_t A, const fq_zech_mpoly_t B, slong var, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_set_fq_zech_poly(fq_zech_mpoly_t A, const fq_zech_poly_t B, slong var, const fq_zech_mpoly_ctx_t ctx)
void _fq_zech_mpoly_set_fq_zech_poly(fq_zech_mpoly_t A, flint_bitcnt_t Abits, const fq_zech_struct * Bcoeffs, slong Blen, slong var, const fq_zech_mpoly_ctx_t ctx)

Test whether `A` involves only the variable of index `var`, convert such
a polynomial to an `fq_zech_poly` (returning 0 otherwise), and set `A` to
a univariate polynomial in that variable.
:::

::: function
int fq_zech_bpoly_factor_smprime(fq_zech_poly_t c, fq_zech_tpoly_t F, fq_zech_bpoly_t B, int allow_shift, const fq_zech_ctx_t ctx)

Factors the bivariate `B`, separable in the main variable: sets `c` to its
content and `F` to its irreducible factors, by evaluation of the inner
variable at a point of $\mathbb{F}_q$ (only 0 unless `allow_shift`),
Hensel lifting and recombination by linear algebra [@Lec2006]. Returns 0
if no suitable point exists.
:::

::: function
int fq_zech_bpoly_hlift2(fq_zech_bpoly_t A, fq_zech_bpoly_t B0, fq_zech_bpoly_t B1, const fq_zech_t alpha, slong degree_inner, const fq_zech_ctx_t ctx)
int fq_zech_bpoly_hlift(slong r, fq_zech_bpoly_t A, fq_zech_bpoly_struct * B, const fq_zech_t alpha, slong degree_inner, const fq_zech_ctx_t ctx)
int fq_zech_polyu3_hlift(slong r, fq_zech_polyun_struct * BB, fq_zech_polyu_t A, fq_zech_polyu_struct * B, const fq_zech_t beta, slong degree_inner, const fq_zech_ctx_t ctx)

Hensel lifting of a factorisation modulo $y - \alpha$ into two or `r`
factors, up to degree `degree_inner` in $y$. Return 1 for success.
:::

### Multivariate lifting

::: type
fq_zech_mpolyv_struct
fq_zech_mpolyv_t
fq_zech_mpoly_pfrac_struct
fq_zech_mpoly_pfrac_t

Vectors of polynomials, and precomputed data for multivariate partial
fractions.
:::

::: function
void fq_zech_mpolyv_init(fq_zech_mpolyv_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpolyv_clear(fq_zech_mpolyv_t A, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpolyv_swap(fq_zech_mpolyv_t A, fq_zech_mpolyv_t B, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpolyv_fit_length(fq_zech_mpolyv_t A, slong length, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpolyv_set_coeff(fq_zech_mpolyv_t A, slong i, fq_zech_mpoly_t c, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpolyv_print_pretty(const fq_zech_mpolyv_t poly, const char ** x, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_to_mpolyv(fq_zech_mpolyv_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_t xalpha, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_from_mpolyv(fq_zech_mpoly_t A, const fq_zech_mpolyv_t B, const fq_zech_mpoly_t xalpha, const fq_zech_mpoly_ctx_t ctx)

Vectors of polynomials, and conversion between `B` and the coefficients
of its expansion in powers of `xalpha` $= x_j - \alpha$.
:::

::: function
int fq_zech_mpoly_univar_content_mpoly(fq_zech_mpoly_t g, const fq_zech_mpoly_univar_t A, const fq_zech_mpoly_ctx_t ctx)

Sets `g` to the GCD of the coefficients of `A` and returns 1 for success.
:::

::: function
void _fq_zech_mpoly_get_lead0(fq_zech_mpoly_t c, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
void _fq_zech_mpoly_set_lead0(fq_zech_mpoly_t A, const fq_zech_mpoly_t B, const fq_zech_mpoly_t c, const fq_zech_mpoly_ctx_t ctx)

Get, or replace by `c`, the leading coefficient of a polynomial viewed as
a polynomial in variable 0.
:::

::: function
int fq_zech_mpoly_pfrac_init(fq_zech_mpoly_pfrac_t Iv, flint_bitcnt_t bits, slong l, slong r, const fq_zech_mpoly_struct * betas, const fq_zech_struct * alpha, const fq_zech_mpoly_ctx_t ctx)
void fq_zech_mpoly_pfrac_clear(fq_zech_mpoly_pfrac_t Iv, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_pfrac(slong r, fq_zech_mpoly_t t, const slong * deg, fq_zech_mpoly_pfrac_t Iv, const fq_zech_mpoly_ctx_t ctx)
int fq_zech_mpoly_hlift(slong m, fq_zech_mpoly_struct * f, slong r, const fq_zech_struct * alpha, const fq_zech_mpoly_t A, const slong * degs, const fq_zech_mpoly_ctx_t ctx)

Multivariate partial fractions and Hensel lifting, as the corresponding
`nmod_mpoly` functions.
:::

::: function
int fq_zech_mpoly_factor_lcc_wang(fq_zech_mpoly_struct * lc_divs, const fq_zech_mpoly_factor_t lcAfac, const fq_zech_poly_t Auc, const fq_zech_bpoly_struct * Auf, slong r, const fq_zech_poly_struct * alpha, const fq_zech_mpoly_ctx_t ctx)

Wang's leading coefficient distribution [@Wan1978].
:::

::: function
int fq_zech_mpoly_factor_irred_smprime_zassenhaus(fq_zech_mpolyv_t fac, const fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx, flint_rand_t state)
int fq_zech_mpoly_factor_irred_smprime_wang(fq_zech_mpolyv_t fac, const fq_zech_mpoly_t A, const fq_zech_mpoly_factor_t lcAfac, const fq_zech_mpoly_t lcA, const fq_zech_mpoly_ctx_t ctx, flint_rand_t state)
int fq_zech_mpoly_factor_irred_smprime_zippel(fq_zech_mpolyv_t fac, const fq_zech_mpoly_t A, const fq_zech_mpoly_factor_t lcAfac, const fq_zech_mpoly_t lcA, const fq_zech_mpoly_ctx_t ctx, flint_rand_t state)

Irreducible factorisation of a monic squarefree polynomial in at least
three variables with evaluation points in $\mathbb{F}_q$, as the
`fq_nmod_mpoly` versions. Return 1 for success, 0 if they ran out of
evaluation points, and $-1$ for a failure that other points would not
fix.
:::

::: function
void _fq_zech_mpoly_monomial_evals(fq_zech_struct * E, const ulong * Aexps, flint_bitcnt_t Abits, slong Alen, const fq_zech_struct * alpha, slong vstart, const fq_zech_mpoly_ctx_t ctx)
int _fq_zech_mpoly_eval_rest_fq_zech_poly(fq_zech_poly_struct * E, slong * starts, slong * ends, slong * stops, ulong * es, const fq_zech_struct * Acoeffs, const ulong * Aexps, slong Alen, slong var, const fq_zech_poly_struct * alphas, const slong * offsets, const slong * shifts, slong N, ulong mask, slong nvars, const fq_zech_ctx_t ctx)
void _fq_zech_mpoly_eval_to_bpoly(fq_zech_bpoly_t E, const fq_zech_mpoly_t A, const fq_zech_poly_struct * alphabetas, const fq_zech_mpoly_ctx_t ctx)

Evaluate the monomials of a polynomial at given values of some variables,
and evaluate all but the first one or two variables to get univariate or
bivariate images.
:::

::: function
int fq_zech_next(fq_zech_t x, const fq_zech_ctx_t ctx)

Steps `x` to the next element of $\mathbb{F}_q$ in a fixed enumeration;
returns 0 when it wraps around.
:::
