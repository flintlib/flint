---
title: nmod_mpoly_factor.h
header: nmod_mpoly_factor.h
description: Factorisation of multivariate polynomials over $\mathbb{F}_p$ for a word-size prime $p$, into content, squarefree, separable and irreducible factors, with extension-field methods for small fields.
category: Multivariate polynomials
order: 8
---

This module factors polynomials of `nmod_mpoly.h` over a prime field
$\mathbb{F}_p$ with word-size $p$. A factorisation is stored in an
`nmod_mpoly_factor_t` as a constant times a product of powers of
polynomials. `nmod_mpoly_factor` computes the irreducible factorisation,
`nmod_mpoly_factor_squarefree` a squarefree one; the content of a
polynomial with respect to each variable is split off first, and in
characteristic $p$ the squarefree decomposition also extracts $p$-th
powers, which over a prime field are obtained by dividing exponents by $p$.

The irreducible factors of a squarefree polynomial are found by evaluating
all variables but two at random points, factoring the bivariate image (by
Hensel lifting and linear-algebra recombination of a univariate
factorisation), and lifting the bivariate factors back one variable at a
time, either by Wang's leading coefficient method with Hensel lifting
[@Wan1978], by sparse Hensel lifting in the style of Zippel and of Monagan
and Tuncer [@MonTun2020], or by Zassenhaus-style recombination
[@Zas1969]. When $\mathbb{F}_p$ is too small to supply good evaluation
points, the same algorithms run over extension fields $\mathbb{F}_{p^k}$,
and the factors found there are combined into factors over
$\mathbb{F}_p$ along their orbits under the Frobenius automorphism; see
[](algorithms/multivariate/modp-factoring), and
[](algorithms/multivariate/hensel-lifting) for the lifting steps.

::: summary algorithms/multivariate/modp-factoring :::

The factorisation functions are described in detail in
[](nmod_mpoly_factor/factor).

## Types and memory management {#types}

::: type
nmod_mpoly_factor_struct
nmod_mpoly_factor_t

A factorisation $c \prod_{i=0}^{k-1} f_i^{e_i}$ of a polynomial over
$\mathbb{F}_p$, with the fields `constant` ($c$, a `ulong`), `poly` (the
bases $f_i$, an array of `nmod_mpoly_struct`), `exp` (the exponents
$e_i$, an array of `fmpz`), `num` ($k$) and `alloc`. All bases belong to
the context passed with the factorisation.
:::

::: function
void nmod_mpoly_factor_init(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_factor_init2(nmod_mpoly_factor_t f, slong alloc, const nmod_mpoly_ctx_t ctx)

Initialise `f` as the empty product with constant 1;
`nmod_mpoly_factor_init2` also allocates space for `alloc` factors.
:::

::: function
void nmod_mpoly_factor_clear(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

Frees the memory used by `f`.
:::

::: function
void nmod_mpoly_factor_realloc(nmod_mpoly_factor_t f, slong alloc, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_factor_fit_length(nmod_mpoly_factor_t f, slong len, const nmod_mpoly_ctx_t ctx)

Change the allocation of `f` to `alloc` factors (truncating the product if
it is longer), or ensure space for at least `len` factors.
:::

## Basic manipulation {#basic}

::: function
void nmod_mpoly_factor_swap(nmod_mpoly_factor_t f, nmod_mpoly_factor_t g, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_factor_set(nmod_mpoly_factor_t f, const nmod_mpoly_factor_t g, const nmod_mpoly_ctx_t ctx)

Swap `f` and `g` efficiently, and set `f` to a copy of `g`.
:::

::: function
void nmod_mpoly_factor_one(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

Sets `f` to the empty product with constant 1.
:::

::: function
slong nmod_mpoly_factor_length(const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
ulong nmod_mpoly_factor_get_constant_ui(const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

Return the number $k$ of factors and the constant $c$ of `f`.
:::

::: function
void nmod_mpoly_factor_get_base(nmod_mpoly_t p, const nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_factor_swap_base(nmod_mpoly_t p, nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t ctx)

Set `p` to the base $f_i$ of index $i$, or swap `p` with it, where $0 \le
i < k$.
:::

::: function
slong nmod_mpoly_factor_get_exp_si(nmod_mpoly_factor_t f, slong i, const nmod_mpoly_ctx_t ctx)

Returns the exponent $e_i$ of index $i$, which must fit an `slong`.
:::

::: function
void nmod_mpoly_factor_append_ui(nmod_mpoly_factor_t f, const nmod_mpoly_t A, ulong e, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_factor_append_fmpz(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const fmpz_t e, const nmod_mpoly_ctx_t ctx)

Append the factor $A^e$ to `f`, copying `A`.
:::

::: function
int nmod_mpoly_factor_expand(nmod_mpoly_t A, const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_factor_matches(const nmod_mpoly_t a, const nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

`nmod_mpoly_factor_expand` sets `A` to the product $c \prod_i f_i^{e_i}$
and returns 1, or returns 0 if an exponent is negative or a power is too
large to compute. `nmod_mpoly_factor_matches` returns whether `f` expands
to `a`.
:::

::: function
void nmod_mpoly_factor_sort(nmod_mpoly_factor_t f, const nmod_mpoly_ctx_t ctx)

Sorts the factors of `f` by increasing exponent, and factors with equal
exponents by `nmod_mpoly_cmp`.
:::

::: function
int nmod_mpoly_factor_cmp(const nmod_mpoly_factor_t A, const nmod_mpoly_factor_t B, const nmod_mpoly_ctx_t ctx)

Compares two factorisations, first by constant, then by number of
factors, then factor by factor (exponent, then base by `nmod_mpoly_cmp`).
Returns $-1$, 0 or 1. Factorisations of the same polynomial compare equal
after `nmod_mpoly_factor_sort`.
:::

::: function
void nmod_mpoly_factor_print_pretty(const nmod_mpoly_factor_t f, const char ** vars, const nmod_mpoly_ctx_t ctx)

Prints `f` as the constant followed by one line `*(base)^exp` per factor,
with the variable names `vars`.
:::

## Factorisation {#factorisation}

These functions assume that the modulus is prime, which is not checked.
They set `f` to a factorisation of `A`, independently of the initial value
of `f`, and return 1 for success. They return 0 (leaving `f` undefined)
only if the exponents of `A` need more than one word per field after the
content has been removed, or an algorithm runs out of evaluation points or
extension fields, which does not happen for inputs of practical size. In
the result, the constant is the leading coefficient of `A` and the bases
are monic. A zero polynomial gives the constant 0 and no factors.

::: function
int nmod_mpoly_factor_content(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Sets `f` to a factorisation of `A` whose bases are variables and
polynomials that are primitive with respect to every variable, that is,
have no factor independent of some variable that they involve. This is
the first step of the other functions.
:::

::: function
int nmod_mpoly_factor_squarefree(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Sets `f` to a factorisation of `A` whose bases are squarefree and pairwise
coprime. Bases with equal exponents are not multiplied together; to get
the product of all irreducible factors of a given multiplicity, sort with
`nmod_mpoly_factor_sort` and multiply the bases with that exponent.
:::

::: function
int nmod_mpoly_factor_separable(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, int sep)

With `sep` = 0 this is `nmod_mpoly_factor_squarefree`. With `sep` = 1,
each base $g$ in addition satisfies, for every variable $x_i$, either
$\partial g/\partial x_i = 0$ or $\gcd(g, \partial g/\partial x_i) = 1$,
with $\partial g/\partial x_i \ne 0$ for at least one $i$; this is the
form required by the irreducible factorisation.
:::

::: function
int nmod_mpoly_factor(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)

Sets `f` to the factorisation of `A` into irreducible factors over
$\mathbb{F}_p$: the bases are monic, irreducible and pairwise distinct.
:::

::: function
int nmod_mpoly_factor_zassenhaus(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_factor_wang(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_factor_zippel(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_factor_algo(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, unsigned int algo)

As `nmod_mpoly_factor`, but lifting the bivariate factorisation with only
one method: Zassenhaus-style recombination after lifting a univariate
factorisation, Wang's leading coefficient method with dense Hensel lifting,
or sparse Hensel lifting. `nmod_mpoly_factor_algo` takes a bit mask
`algo` of `MPOLY_FACTOR_USE_ZAS`, `MPOLY_FACTOR_USE_WANG`,
`MPOLY_FACTOR_USE_ZIP` (or `MPOLY_FACTOR_USE_ALL`). Polynomials in one or
two variables are always factored by the univariate and bivariate
algorithms.
:::

## Internal functions {#internal}

These functions are used by the factoring algorithms and by the GCD
algorithms of `nmod_mpoly`. Their interfaces may change without notice.

### Conversions and vectors

::: function
void nmod_mpoly_get_bpoly(n_bpoly_t A, const nmod_mpoly_t B, slong var0, slong var1, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_set_bpoly(nmod_mpoly_t A, flint_bitcnt_t Abits, const n_bpoly_t B, slong var0, slong var1, const nmod_mpoly_ctx_t ctx)

Convert between `B`, a polynomial in the variables of indices `var0` and
`var1` only, and the dense bivariate `n_bpoly_t` `A` in which `var0` is
the main (outer) variable; `nmod_mpoly_set_bpoly` produces exponents of
`Abits` bits.
:::

::: function
void nmod_mpoly_get_polyu1n(n_polyun_t A, const nmod_mpoly_t B, slong varx, slong vary, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_set_polyu1n(nmod_mpoly_t B, const n_polyun_t A, slong varx, slong vary, const nmod_mpoly_ctx_t ctx)

Convert between a polynomial in the variables `varx` and `vary` and an
`n_polyun_t` with exponents in `varx` and `n_poly_t` coefficients in
`vary`.
:::

::: function
void _nmod_mpoly_get_lead0(nmod_mpoly_t c, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_set_lead0(nmod_mpoly_t A, const nmod_mpoly_t B, const nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)

`_nmod_mpoly_get_lead0` sets `c` to the leading coefficient of `A` viewed
as a polynomial in variable 0 (a polynomial in the other variables).
`_nmod_mpoly_set_lead0` sets `A` to `B` with that leading coefficient
replaced by `c`.
:::

::: type
nmod_mpolyv_struct
nmod_mpolyv_t

A vector of `nmod_mpoly_t`, with the fields `coeffs`, `alloc` and
`length`.
:::

::: function
void nmod_mpolyv_init(nmod_mpolyv_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyv_clear(nmod_mpolyv_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyv_swap(nmod_mpolyv_t A, nmod_mpolyv_t B, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyv_fit_length(nmod_mpolyv_t A, slong length, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyv_set_coeff(nmod_mpolyv_t A, slong i, nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyv_print_pretty(const nmod_mpolyv_t poly, const char ** x, const nmod_mpoly_ctx_t ctx)

Memory management of vectors; `nmod_mpolyv_set_coeff` swaps `c` into
position $i$, extending the vector with zeros.
:::

::: function
void nmod_mpoly_to_mpolyv(nmod_mpolyv_t A, const nmod_mpoly_t B, const nmod_mpoly_t xalpha, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_from_mpolyv(nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_mpolyv_t B, const nmod_mpoly_t xalpha, const nmod_mpoly_ctx_t ctx)

Convert between `B` and the vector of its coefficients in the expansion
$B = \sum_i A_i\, \mathit{xalpha}^i$ in powers of `xalpha`, which is of
the form $x_j - \alpha$; the coefficients do not involve $x_j$.
`nmod_mpoly_from_mpolyv` produces exponents of `Abits` bits.
:::

::: function
int _nmod_mpoly_vec_content_mpoly(nmod_mpoly_t g, const nmod_mpoly_struct * A, slong Alen, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_vec_divexact_mpoly(nmod_mpoly_struct * A, slong Alen, const nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_vec_mul_mpoly(nmod_mpoly_struct * A, slong Alen, const nmod_mpoly_t c, const nmod_mpoly_ctx_t ctx)

The GCD of the `Alen` polynomials `A` (returning 1 for success), and
exact division or multiplication of all of them by `c`.
:::

::: function
slong _n_poly_vec_max_degree(const n_poly_struct * A, slong Alen)
void _n_poly_vec_mul_nmod_intertible(n_poly_struct * A, slong Alen, ulong c, nmod_t ctx)
void _n_poly_vec_mod_mul_poly(n_poly_struct * A, slong Alen, const n_poly_t g, const nmod_t ctx)
void _n_poly_vec_mod_divexact_poly(n_poly_struct * A, slong Alen, const n_poly_t g, nmod_t ctx)
void _n_poly_vec_mod_content(n_poly_t g, const n_poly_struct * A, slong Alen, nmod_t ctx)
void _n_poly_vec_mod_remove_content(n_poly_t g, n_poly_struct * A, slong Alen, nmod_t ctx)

Operations on vectors of univariate polynomials: maximum degree,
multiplication by an invertible scalar or by `g`, exact division by `g`,
the monic GCD `g` of the entries, and division by it.
:::

::: function
void nmod_mpoly_compression_do(nmod_mpoly_t L, const nmod_mpoly_ctx_t Lctx, ulong * Acoeffs, slong Alen, mpoly_compression_t M)
void nmod_mpoly_compression_undo(nmod_mpoly_t A, flint_bitcnt_t Abits, const nmod_mpoly_ctx_t Actx, nmod_mpoly_t L, const nmod_mpoly_ctx_t Lctx, mpoly_compression_t M)

Apply the compression `M` (computed by `mpoly_compression_set`, see
[](mpoly)) to the terms `Acoeffs` of length `Alen`, giving `L` in the
context `Lctx` with fewer variables, and undo it, giving `A` with `Abits`
exponent bits. A compression is an invertible monomial change of variables
that removes unused variables and redundancy in the exponents.
:::

::: function
int nmod_mpolyu_is_canonical(const nmod_mpolyu_t A, const nmod_mpoly_ctx_t ctx)
void nmod_mpolyu3_print_pretty(const nmod_mpolyu_t A, const char * var0, const char * var1, const char * var2, const char ** vars, const nmod_mpoly_ctx_t ctx)

Check an `nmod_mpolyu_t`, and print one whose exponents pack three main
variables with the names `var0`, `var1`, `var2`.
:::

### Bivariate factorisation

::: function
int n_bpoly_mod_factor_smprime(n_poly_t c, n_tpoly_t F, n_bpoly_t B, int allow_shift, nmod_t ctx)

Factors the bivariate polynomial `B` in $x$ (main variable) and $y$,
which must be separable in $x$: sets `c` to its content with respect to
$x$ and `F` to its irreducible factors. Evaluates $y = \alpha$ for
$\alpha \in \mathbb{F}_p$ (only $\alpha = 0$ if `allow_shift` is zero),
factors the univariate image, lifts the factors $y$-adically by Hensel
lifting and recombines them with the linear algebra on logarithmic
derivatives of Lecerf [@Lec2006]
([](algorithms/multivariate/modp-factoring#bivariate)). Returns 0 if no suitable $\alpha$
exists in $\mathbb{F}_p$. `B` is destroyed.
:::

::: function
void n_bpoly_mod_factor_lgprime(n_poly_t c, n_tpoly_t F, n_bpoly_t B, nmod_t ctx)

As `n_bpoly_mod_factor_smprime`, but reduces $y$ modulo irreducible
polynomials $m(y)$ of increasing degree, factors the image over
$\mathbb{F}_p[y]/(m) \cong \mathbb{F}_{p^d}$, lifts $m$-adically and
recombines. Always succeeds.
:::

::: function
int nmod_mat_is_reduced(const nmod_mat_t N)
void nmod_mat_init_nullspace_tr(nmod_mat_t X, nmod_mat_t tmp)

Helpers of the recombination: `nmod_mat_is_reduced` returns whether every
column of `N` contains exactly one nonzero entry, equal to 1, so that the
rows describe a partition of the local factors.
`nmod_mat_init_nullspace_tr` initialises `X` to a matrix whose rows are a
basis of the right nullspace of `tmp` (which is destroyed).
:::

::: function
int n_bpoly_mod_pfrac(slong r, n_bpoly_struct * C, slong * C_deg1_bound, n_bpoly_t A, n_bpoly_struct * B, nmod_t mod)
int n_bpoly_mod_hlift2(n_bpoly_t A, n_bpoly_t B0, n_bpoly_t B1, ulong alpha, slong degree_inner, nmod_t mod, n_poly_bpoly_stack_t St)
int n_bpoly_mod_hlift2_cubic(n_bpoly_t A, n_bpoly_t B0, n_bpoly_t B1, ulong alpha, slong degree_inner, nmod_t ctx, nmod_eval_interp_t E, n_poly_bpoly_stack_t St)
int n_bpoly_mod_hlift(slong r, n_bpoly_t A, n_bpoly_struct * B, ulong alpha, slong degree_inner, nmod_t mod, n_poly_bpoly_stack_t St)
int n_bpoly_mod_hlift_cubic(slong r, n_bpoly_t A, n_bpoly_struct * B, ulong alpha, slong degree_inner, nmod_t mod, nmod_eval_interp_t E, n_poly_bpoly_stack_t St)
int n_polyu3_mod_hlift(slong r, n_polyun_struct * BB, n_polyu_t A, n_polyu_struct * B, ulong beta, slong degree_inner, nmod_t ctx)

Hensel lifting of a factorisation of a bivariate (or, for
`n_polyu3_mod_hlift`, trivariate) polynomial `A` modulo $y - \alpha$ to a
factorisation modulo $(y - \alpha)^{d}$ with $d$ = `degree_inner` + 1,
for two factors `B0`, `B1` or `r` factors `B`; the `cubic` versions
evaluate and interpolate in the other variable. `n_bpoly_mod_pfrac`
solves the corresponding partial fraction equation. They return 1 for
success and 0 if the lifted factors do not multiply to `A`.
:::

### Multivariate lifting

::: type
nmod_mpoly_pfrac_struct
nmod_mpoly_pfrac_t

Precomputed data for solving the multivariate partial fraction equations
$\sum_i \delta_i \prod_{j \ne i} \beta_j = t$ during Hensel lifting.
:::

::: function
int nmod_mpoly_pfrac_init(nmod_mpoly_pfrac_t Iv, flint_bitcnt_t bits, slong l, slong r, const nmod_mpoly_struct * betas, const ulong * alpha, const nmod_mpoly_ctx_t ctx)
void nmod_mpoly_pfrac_clear(nmod_mpoly_pfrac_t Iv, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_pfrac(slong r, nmod_mpoly_t t, const slong * deg, nmod_mpoly_pfrac_t Iv, const nmod_mpoly_ctx_t ctx)

`nmod_mpoly_pfrac_init` prepares `Iv` for the `r` factors `betas` in
`l` + 1 variables with the evaluation points `alpha`, and returns 0 if the
factors are not pairwise coprime at the evaluation points.
`nmod_mpoly_pfrac` solves the equation for the right-hand side `t` in the
first `r` + 1 variables with degree bounds `deg`, returning 1 for success.
:::

::: function
int nmod_mpoly_hlift(slong m, nmod_mpoly_struct * f, slong r, const ulong * alpha, const nmod_mpoly_t A, const slong * degs, const nmod_mpoly_ctx_t ctx)
int nmod_mpoly_hlift_zippel(slong m, nmod_mpoly_struct * B, slong r, const ulong * alpha, const nmod_mpoly_t A, const slong * degs, const nmod_mpoly_ctx_t ctx, flint_rand_t state)

Given factors `f[0]`, ..., `f[r-1]` with $A = \prod_i f_i$ modulo $x_m -
\alpha_{m-1}$ and degree bounds `degs`, lift them to a factorisation of
`A` in the variables $x_0, \ldots, x_m$: by dense Hensel lifting, or by
sparse interpolation of the lifted coefficients. Return 1 for success.
:::

::: function
int nmod_mpoly_factor_lcc_wang(nmod_mpoly_struct * lc_divs, const nmod_mpoly_factor_t lcAfac, const n_poly_t Auc, const n_bpoly_struct * Auf, slong r, const n_poly_struct * alpha, const nmod_mpoly_ctx_t ctx)

Wang's leading coefficient distribution [@Wan1978] along a curve: given
the factorisation `lcAfac` of the leading coefficient of $A$, and the
content `Auc` $\in \mathbb{F}_p[Y]$ and the `r` factors `Auf` of the
bivariate image of $A$ under the substitution of the polynomials
`alpha[i]` $\in \mathbb{F}_p[Y]$ for the variables other than the main
one, sets `lc_divs` to the true leading coefficients of the factors (up
to constants) ([](algorithms/multivariate/modp-factoring#lcc)). Returns 1
for success and 0 if the curve is unlucky.
:::

::: function
int nmod_mpoly_factor_irred_smprime_zassenhaus(nmod_mpolyv_t fac, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, flint_rand_t state)
int nmod_mpoly_factor_irred_smprime_wang(nmod_mpolyv_t fac, const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac, const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state)
int nmod_mpoly_factor_irred_smprime_zippel(nmod_mpolyv_t fac, const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac, const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state)
int nmod_mpoly_factor_irred_medprime_zassenhaus(nmod_mpolyv_t fac, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, flint_rand_t state)
int nmod_mpoly_factor_irred_medprime_wang(nmod_mpolyv_t Af, const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac, const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state)
int nmod_mpoly_factor_irred_medprime_zippel(nmod_mpolyv_t Af, const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac, const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state)
int nmod_mpoly_factor_irred_lgprime_zassenhaus(nmod_mpolyv_t fac, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, flint_rand_t state)
int nmod_mpoly_factor_irred_lgprime_wang(nmod_mpolyv_t Af, const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac, const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state)
int nmod_mpoly_factor_irred_lgprime_zippel(nmod_mpolyv_t Af, const nmod_mpoly_t A, const nmod_mpoly_factor_t lcAfac, const nmod_mpoly_t lcA, const nmod_mpoly_ctx_t ctx, flint_rand_t state)

Set `fac` (or `Af`) to the irreducible factors of the monic squarefree
`A` in lexicographic order with at least three variables, given for the
Wang and Zippel methods the factorisation `lcAfac` of the leading
coefficient `lcA` in variable 0. The `smprime` versions use evaluation
points in $\mathbb{F}_p$; the `medprime` versions map `A` to
$\mathbb{F}_{p^k}$ represented by `fq_zech` (with $p^k \le 10^6$), factor
there and combine Frobenius-conjugate factors; the `lgprime` versions do
the same with `fq_nmod` and no bound on $k$. Return 1 for success, 0 if
they ran out of evaluation points or extension degrees, and $-1$ for a
failure that other evaluation points would not fix (for example exponents
that do not fit one word).
:::

::: function
int _nmod_mpoly_factor_separable(nmod_mpoly_factor_t f, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx, int sep)

As `nmod_mpoly_factor_separable` for a monic `A` whose content has been
removed.
:::

::: function
int _nmod_mpoly_evaluate_rest_n_poly(n_poly_struct * E, slong * starts, slong * ends, slong * stops, ulong * es, const ulong * Acoeffs, const ulong * Aexps, slong Alen, slong var, const n_poly_struct * alphas, const slong * offsets, const slong * shifts, slong N, ulong mask, slong nvars, nmod_t ctx)
void _nmod_mpoly_eval_rest_to_n_bpoly(n_bpoly_t E, const nmod_mpoly_t A, const n_poly_struct * alphabetas, const nmod_mpoly_ctx_t ctx)
void _nmod_mpoly_set_n_bpoly_var1_zero(nmod_mpoly_t A, flint_bitcnt_t Abits, const n_bpoly_t B, slong var, const nmod_mpoly_ctx_t ctx)

Evaluate all variables after the first one or two at given points (or at
given univariate polynomials), giving univariate or bivariate images; and
set `A` from a bivariate polynomial in variable 0 and `var`.
:::

### Sparse interpolation

::: function
int nmod_mpolyl_gcdp_zippel_smprime(nmod_mpoly_t G, nmod_mpoly_t Abar, nmod_mpoly_t Bbar, nmod_mpoly_t A, nmod_mpoly_t B, slong var, const nmod_mpoly_ctx_t ctx, flint_rand_t state)
int nmod_mpolyl_gcds_zippel(nmod_mpoly_t G, const ulong * Gmarks, slong Gmarkslen, nmod_mpoly_t A, nmod_mpoly_t B, slong *perm, slong l, slong var, const nmod_mpoly_ctx_t ctx, flint_rand_t state, slong * Gdegbound, n_poly_t Amarks, n_poly_t Bmarks)

Zippel's GCD algorithm in a lexicographic context, as used by the
factoring code: `gcdp` interpolates the variables up to `var` with points
in $\mathbb{F}_p$, `gcds` solves for the coefficients of the assumed
monomial structure `G` (marked by `Gmarks`).
:::

::: function
void mpoly_monomial_evals_nmod(n_poly_t EH, const ulong * Aexps, slong Alen, flint_bitcnt_t Abits, n_poly_struct * alpha_caches, slong start, slong stop, const mpoly_ctx_t mctx, const nmod_t fpctx)
void mpoly1_monomial_evals_nmod(n_polyun_t EH, const ulong * Aexps, flint_bitcnt_t Abits, const ulong * Amarks, slong Amarkslen, n_poly_struct * alpha_caches, slong m, const mpoly_ctx_t mctx, const nmod_t fpctx)
void mpoly2_monomial_evals_nmod(n_polyun_t EH, const ulong * Aexps, flint_bitcnt_t Abits, ulong * Amarks, slong Amarkslen, n_poly_struct * alpha_caches, slong m, const mpoly_ctx_t mctx, const nmod_t fpctx)

Evaluate the monomials of a polynomial at powers of the evaluation points
cached in `alpha_caches`, keeping zero, one or two main variables
symbolic; used to set up the Vandermonde systems of sparse interpolation.
:::

::: function
void n_polyun_zip_start(n_polyun_t Z, n_polyun_t H, slong req_images)
int n_polyu2n_add_zip_must_match(n_polyun_t Z, const n_bpoly_t A, slong cur_length)
int n_polyun_zip_solve(nmod_mpoly_t A, n_polyun_t Z, n_polyun_t H, n_polyun_t M, const nmod_mpoly_ctx_t ctx)

Sparse interpolation of the coefficients of a polynomial with known
monomial structure `H`: start collecting `req_images` images in `Z`, add
the image `A` (returning 0 if its support does not match), and solve the
transposed Vandermonde systems, returning 1 for success and 0 or $-1$ if
the images are inconsistent with the assumed structure.
:::
