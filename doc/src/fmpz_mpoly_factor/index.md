---
title: fmpz_mpoly_factor.h
header: fmpz_mpoly_factor.h
description: Factorisation of multivariate polynomials over the integers into content, squarefree and irreducible factors, with the multimodular GCD, Hensel lifting and dense bivariate machinery behind it.
category: Multivariate polynomials
order: 3
---

This module factors polynomials in $\Z[x_0, \ldots, x_{n-1}]$ (see
[](fmpz_mpoly)). A factorisation is stored as an `fmpz_mpoly_factor_t`:
an integer constant times a product of polynomial bases raised to positive
exponents. `fmpz_mpoly_factor_content` splits off the content and the
powers of variables, `fmpz_mpoly_factor_squarefree` computes a squarefree
factorisation, and `fmpz_mpoly_factor` the factorisation into irreducible
polynomials. The factoring strategy and an example are on the function page
[](fmpz_mpoly_factor/factor), and the algorithms are explained in
[](algorithms/multivariate/factorisation) and
[](algorithms/multivariate/hensel-lifting).

The header also declares the internal machinery shared by the factoring
code and by the GCD code of `fmpz_mpoly`: GCD algorithms on lex contexts,
modular reduction and reconstruction, multivariate partial fractions and
Hensel lifting, and dense bivariate polynomials with their factorisation
and resultant. These functions are documented below for developers; their
interfaces may change between releases.

## Types {#types}

::: type
fmpz_mpoly_factor_struct
fmpz_mpoly_factor_t

A factored polynomial $c \prod_{i < \mathit{num}} p_i^{e_i}$. The fields
are `constant` (the integer $c$), `constant_den` (a denominator for $c$,
which stays 1 in this module and is used by the factored denominators of
`fmpz_mpoly_q`), `poly` (the bases $p_i$, an array of
`fmpz_mpoly_struct`), `exp` (the exponents $e_i$, an array of `fmpz`),
`num` and `alloc`. The bases of a factorisation returned by the factoring
functions are nonconstant with positive leading coefficient, and the
exponents are positive.
:::

## Memory management {#memory}

::: function
void fmpz_mpoly_factor_init(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_factor_init2(fmpz_mpoly_factor_t f, slong alloc, const fmpz_mpoly_ctx_t ctx)

Initialise `f` as the empty product with constant 1, with room for
`alloc` factors in the second version.
:::

::: function
void fmpz_mpoly_factor_realloc(fmpz_mpoly_factor_t f, slong alloc, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_factor_fit_length(fmpz_mpoly_factor_t f, slong len, const fmpz_mpoly_ctx_t ctx)

Change the allocation to exactly `alloc` factors (freeing the factors
beyond it, and everything if $\mathit{alloc} \le 0$), respectively ensure
room for at least `len` factors.
:::

::: function
void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

Frees the memory used by `f`.
:::

## Basic manipulation {#basic}

::: function
void fmpz_mpoly_factor_swap(fmpz_mpoly_factor_t f, fmpz_mpoly_factor_t g, const fmpz_mpoly_ctx_t ctx)

Swaps `f` and `g` in constant time.
:::

::: function
void fmpz_mpoly_factor_set(fmpz_mpoly_factor_t f, const fmpz_mpoly_factor_t g, const fmpz_mpoly_ctx_t ctx)

Sets `f` to a copy of `g` (the field `constant_den` is not copied).
:::

::: function
slong fmpz_mpoly_factor_length(const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

Returns the number of bases of `f`.
:::

::: function
void fmpz_mpoly_factor_get_constant_fmpz(fmpz_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

Set $c$ to the constant of `f`. The `fmpq` version uses `constant_den` as
the denominator.
:::

::: function
void fmpz_mpoly_factor_get_base(fmpz_mpoly_t p, const fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_factor_swap_base(fmpz_mpoly_t p, fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t ctx)

Set $p$ to, respectively swap $p$ with, the base of index $i$.
:::

::: function
slong fmpz_mpoly_factor_get_exp_si(fmpz_mpoly_factor_t f, slong i, const fmpz_mpoly_ctx_t ctx)

Returns the exponent of the base of index $i$, which must fit an `slong`.
:::

::: function
void fmpz_mpoly_factor_set_fmpz(fmpz_mpoly_factor_t f, const fmpz_t a, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_factor_zero(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_factor_one(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

Set `f` to the constant $a$, 0 or 1, with no bases.
:::

::: function
void fmpz_mpoly_factor_append_ui(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, ulong e, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_factor_append_fmpz_swap(fmpz_mpoly_factor_t f, fmpz_mpoly_t A, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)

Append the factor $A^e$ to `f`, copying $A$, respectively swapping it in
(so that $A$ is left with an unspecified value). No normalisation is done.
:::

::: function
void fmpz_mpoly_factor_sort(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

Sorts the factors by increasing exponent, and factors with equal exponents
by `fmpz_mpoly_cmp`. Since the irreducible factorisation returned by
`fmpz_mpoly_factor` has distinct bases with positive leading coefficients,
two such factorisations of the same polynomial are identical after
sorting.
:::

::: function
int fmpz_mpoly_factor_cmp(const fmpz_mpoly_factor_t f, const fmpz_mpoly_factor_t g, const fmpz_mpoly_ctx_t ctx)

Compares `f` and `g` by their constants, then their numbers of factors,
then factor by factor (exponent, then base by `fmpz_mpoly_cmp`). Returns
a negative, zero or positive value. This is a total order on the stored
data, not a comparison of the products.
:::

::: function
int fmpz_mpoly_factor_expand(fmpz_mpoly_t A, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

Sets $A$ to the product represented by `f` (ignoring `constant_den`).
Returns 0 if an exponent is negative or a power cannot be computed, and 1
otherwise.
:::

::: function
int fmpz_mpoly_factor_matches(const fmpz_mpoly_t A, const fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx)

Returns whether `f` expands to $A$.
:::

::: function
void fmpz_mpoly_factor_print_pretty(const fmpz_mpoly_factor_t f, const char ** vars, const fmpz_mpoly_ctx_t ctx)

Prints `f` in the form `c*(p1)^e1*(p2)^e2*...`, using the variable names
`vars` (or `x1`, `x2`, ... if `vars` is `NULL`).
:::

## Factorisation {#factorisation}

The factoring functions return 1 for success and 0 for failure, in which
case `f` is unspecified. They fail only when an intermediate result is
unreasonably large, for example when a polynomial with exponents beyond
one word cannot be repacked, and not for ordinary inputs. They overwrite
`f`: its initial value does not matter.

::: function
int fmpz_mpoly_factor_content(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Sets `f` to a factorisation of $A$ into an integer constant (the content
of $A$ with the sign of its leading coefficient), powers of the
variables, and polynomials with exponent 1 that are primitive with
respect to each variable they involve. The latter are found by repeatedly
splitting off the content with respect to one variable; they are neither
squarefree nor irreducible in general. For constant $A$, `f` is that
constant with no bases.
:::

::: function
int fmpz_mpoly_factor_squarefree(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Sets `f` to a factorisation of $A$ whose bases are primitive, squarefree
and pairwise coprime. A base with exponent $k$ is a product of irreducible
factors of $A$ of multiplicity exactly $k$, but the factors of a given
multiplicity may be spread over several bases (one for each piece found by
`fmpz_mpoly_factor_content`). To get a single base for each multiplicity,
sort with `fmpz_mpoly_factor_sort` and multiply the bases with equal
exponents.
:::

::: function
int fmpz_mpoly_factor(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Sets `f` to the factorisation of $A$ into an integer constant and
distinct irreducible polynomials with positive leading coefficients. The
constant is the signed content of $A$ and is not factored into primes.
See [](fmpz_mpoly_factor/factor).
:::

::: function
int fmpz_mpoly_factor_wang(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_factor_zippel(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
int fmpz_mpoly_factor_zassenhaus(fmpz_mpoly_factor_t f, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Versions of `fmpz_mpoly_factor` restricted to one method for the hard
case (a squarefree primitive polynomial that is not handled by the
special cases for few variables or low degree): Wang's algorithm with
dense Hensel lifting, sparse Hensel lifting in Zippel's style, or lifting
without leading coefficient information followed by recombination. The
result is the same; these functions exist for testing and benchmarking.
:::

::: function
int fmpz_mpoly_factor_irred(fmpz_mpoly_factor_t f, const fmpz_mpoly_ctx_t ctx, unsigned int algo)

Replaces each base of `f` by its irreducible factors, with the exponent
of the base, using the methods in the bit set `algo` of
`MPOLY_FACTOR_USE_ZAS`, `MPOLY_FACTOR_USE_WANG` and `MPOLY_FACTOR_USE_ZIP`
(or `MPOLY_FACTOR_USE_ALL`). Each base must be squarefree, have positive
leading coefficient and be primitive with respect to each variable it
involves.
:::

::: function
int fmpz_mpoly_factor_bound_si(fmpz_t B, const fmpz_t A, const slong * degs, slong nvars)

Sets $B$ to a bound on the coefficients of the divisors of a
polynomial $f$ whose coefficients are bounded by $A$ in absolute value
and whose degrees in its $n$ = `nvars` variables are $d_0, \ldots,
d_{n-1}$ (`degs`): $B = |A| \cdot (\lfloor \sqrt{\lceil \prod_i
(d_i + 1) / 2^{n} \rceil} \rfloor + 1) \cdot 2^{\sum_i d_i}$,
and $B = 0$ if some $d_i < 0$. Returns 0 if the degree sum overflows and
1 otherwise. `fmpz_mpolyl_gcd_zippel` uses this bound to decide when to
give up on a modular reconstruction and start over.
:::

## GCD algorithms on lex contexts {#gcd}

The GCD driver works in a lex context in the variables that occur in the
inputs, in an order chosen for the algorithm, with exponents shifted by
the minimal exponents and divided by common strides. The functions named
`fmpz_mpolyl_*` require such a lex context and single-word exponents, with
all inputs and outputs packed with the same number of bits. See
[](fmpz_mpoly/gcd) and [](algorithms/multivariate/gcd).

::: function
int _fmpz_mpoly_gcd_algo(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, unsigned int algo)

The GCD driver behind `fmpz_mpoly_gcd` and its variants, for nonzero $A$
and $B$: handles the easy cases, estimates degrees and term counts of the
GCD, and dispatches to the algorithms allowed by the bit set `algo` of
`MPOLY_GCD_USE_*` flags. `Abar` and `Bbar` may be `NULL` if the cofactors
are not wanted. Returns 1 for success.
:::

::: function
void fmpz_mpoly_to_mpolyl_perm_deflate(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t lctx, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const slong * perm, const ulong * shift, const ulong * stride)
void fmpz_mpoly_from_mpolyl_perm_inflate(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t ctx, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t lctx, const slong * perm, const ulong * shift, const ulong * stride)

Convert from `ctx` to the context `lctx`, whose variable $k$ is variable
`perm[k]` of `ctx` with exponents mapped by $e \mapsto (e - s)/t$ for the
`shift` $s$ and `stride` $t$ of that variable, and back, packing the
result into `Abits` bits. The results are sorted. See
`mpoly_to_mpolyl_perm_deflate`.
:::

::: function
int fmpz_mpolyl_gcd_brown(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, fmpz_mpoly_t A, fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const mpoly_gcd_info_t Iv)
int fmpz_mpolyl_gcd_brown_threaded_pool(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, fmpz_mpoly_t A, fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, const mpoly_gcd_info_t Iv, const thread_pool_handle * handles, slong num_handles)

Brown's dense multimodular algorithm [@Bro1971] for the GCD and cofactors
of nonzero $A$ and $B$. The integer contents are removed from $A$ and $B$,
which are modified. Images are computed with
`nmod_mpolyn_gcd_brown_smprime` modulo primes $p > 2^{\mathrm{FLINT\_BITS}-1}$
not dividing the GCD $\gamma$ of the leading coefficients, scaled to have
leading coefficient $\gamma \bmod p$ and combined by Chinese remaindering,
discarding images of too large degree. The loop stops when the product $m$
of the primes exceeds twice $\gamma$ times the larger height of $A$ and
$B$ and the reconstructed $H$, $\bar A_H$, $\bar B_H$ satisfy
$2 \min(\|H\|_\infty \|\bar A_H\|_1, \|H\|_1 \|\bar A_H\|_\infty) < m$ and
likewise for $B$; this proves $H \bar A_H = \gamma A$ and $H \bar B_H =
\gamma B$ without a division. `Iv` (possibly `NULL`) supplies degree
bounds for the images. The threaded version computes images for
different primes in parallel using the given thread pool handles; it is
not used by `fmpz_mpoly_gcd`. Return 0 if the primes run out.
:::

::: function
int fmpz_mpolyl_gcd_zippel(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx, flint_rand_t randstate)

Zippel's sparse modular GCD [@Zip1979] for $A$ and $B$ in at least two
variables that are primitive as polynomials in $x_0$ over $\Z[x_1, \ldots,
x_{n-1}]$: the GCD modulo one prime fixes the monomials (the shape), and the
images modulo further primes are computed by solving for the coefficients
of that shape. The coefficients are reconstructed by Chinese remaindering,
and the result is verified by division. The bound of
`fmpz_mpoly_factor_bound_si` on the inputs decides when a reconstruction is
abandoned.
:::

::: function
int fmpz_mpolyl_gcd_zippel2(fmpz_mpoly_t G, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t Gamma, const fmpz_mpoly_ctx_t ctx)

The sparse GCD of Hu and Monagan [@HuMon2016] for $A$ and $B$ in at least
three variables that are primitive as polynomials in $x_0, x_1$ over
$\Z[x_2, \ldots, x_{n-1}]$: bivariate GCDs in $x_0, x_1$ modulo primes, with
the remaining variables evaluated at powers of random points, are
interpolated with the Berlekamp–Massey algorithm, and the result is verified
by division. `Gamma`, a polynomial in $x_2, \ldots$ that is a multiple of
the leading coefficient of the GCD in $x_0, x_1$ (in practice the GCD of
the leading coefficients of $A$ and $B$), is imposed on the images.
:::

::: function
int fmpz_mpolyl_gcd_hensel(fmpz_mpoly_t G, slong Gdeg, fmpz_mpoly_t Abar, fmpz_mpoly_t Bbar, const fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)

The EEZ-GCD algorithm [@Wan1980] for $A$ and $B$ in at least two variables
that are primitive as polynomials in $x_0$, given an upper bound `Gdeg` on
the degree in $x_0$ of the GCD: the variables $x_1, \ldots$ are evaluated
at small integers, the univariate GCD is computed, and the factorisation
$G \cdot \bar A$ of $A$, of $B$, or of a combination $\mu_1 A + \mu_2 B$ in
which the two factors are coprime, is lifted with
`fmpz_mpoly_hlift`.
:::

## Modular reduction and reconstruction {#interp}

::: function
void fmpz_mpoly_interp_reduce_p(nmod_mpoly_t Ap, const nmod_mpoly_ctx_t ctxp, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_interp_lift_p(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx, nmod_mpoly_t Ap, const nmod_mpoly_ctx_t ctxp)

Set $A_p = A \bmod p$, and $A$ to the lift of $A_p$ with coefficients in
the symmetric range $(-p/2, p/2]$. Both contexts have the same variables
and the polynomials are packed with the same number of bits.
:::

::: function
int fmpz_mpoly_interp_mcrt_p(flint_bitcnt_t * coeffbits, fmpz_mpoly_t H, const fmpz_mpoly_ctx_t ctx, const fmpz_t m, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctxp)

Replaces $H$, known modulo $m$, by the polynomial with coefficients in
the symmetric range modulo $mp$ that is congruent to $H$ modulo $m$ and to
$A$ modulo $p$. $H$ and $A$ must have the same monomials. Returns whether
$H$ changed, and sets `*coeffbits` to the largest bit length of its new
coefficients.
:::

::: function
void fmpz_mpoly_interp_reduce_p_mpolyn(nmod_mpolyn_t E, const nmod_mpoly_ctx_t pctx, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_interp_lift_p_mpolyn(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx, const nmod_mpolyn_t B, const nmod_mpoly_ctx_t pctx)
int fmpz_mpoly_interp_crt_p_mpolyn(fmpz_mpoly_t F, fmpz_mpoly_t T, const fmpz_mpoly_ctx_t ctx, fmpz_t modulus, const nmod_mpolyn_t A, const nmod_mpoly_ctx_t pctx)

The same operations with the image in recursive form
$\F_p[x_0, \ldots, x_{n-2}][x_{n-1}]$ (an `nmod_mpolyn_t`, whose
coefficients are univariate in the last variable), as used by Brown's
algorithm. The CRT version updates $F$ so that it is unchanged modulo
`modulus` and congruent to $A$ modulo $p$, without assuming that the
monomials match; it uses $T$ as scratch space, returns whether $F$
changed, and does not update `modulus`.
:::

## Vectors of polynomials {#mpolyv}

::: type
fmpz_mpolyv_struct
fmpz_mpolyv_t

A growable vector of polynomials (fields `coeffs`, `alloc`, `length`),
used for lists of factors and for expansions in powers of $x_i - \alpha$.
:::

::: function
void fmpz_mpolyv_init(fmpz_mpolyv_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpolyv_clear(fmpz_mpolyv_t A, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpolyv_swap(fmpz_mpolyv_t A, fmpz_mpolyv_t B, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpolyv_fit_length(fmpz_mpolyv_t A, slong length, const fmpz_mpoly_ctx_t ctx)

Initialise an empty vector, free it, swap two vectors in constant time,
and ensure room for `length` entries.
:::

::: function
void fmpz_mpolyv_print_pretty(const fmpz_mpolyv_t poly, const char ** x, const fmpz_mpoly_ctx_t ctx)

Prints the entries, one per line, as `coeff[i]: ...`.
:::

::: function
void fmpz_mpolyv_set_coeff(fmpz_mpolyv_t A, slong i, fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

Sets entry $i$ to the nonzero polynomial $c$ by swapping (so $c$ is left
with an unspecified value), extending the vector with zero entries if
needed.
:::

::: function
void fmpz_mpoly_to_mpolyv(fmpz_mpolyv_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t xalpha, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_from_mpolyv(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpolyv_t B, const fmpz_mpoly_t xalpha, const fmpz_mpoly_ctx_t ctx)

The first writes $B = \sum_i A_i \cdot \mathit{xalpha}^i$ by repeated
division with remainder by `xalpha`; for $\mathit{xalpha} = x_j -
\alpha$ the $A_i$ are the Taylor coefficients of $B$ at $x_j = \alpha$.
The second evaluates such an expansion by Horner's rule and packs the
result into `Abits` bits.
:::

::: function
int _fmpz_mpoly_vec_content_mpoly(fmpz_mpoly_t g, const fmpz_mpoly_struct * A, slong Alen, const fmpz_mpoly_ctx_t ctx)

Sets $g$ to the GCD of the `Alen` polynomials in $A$, with positive
leading coefficient, starting with the two shortest ones. Returns 0 if a
GCD computation failed.
:::

::: function
void _fmpz_mpoly_vec_divexact_mpoly(fmpz_mpoly_struct * A, slong Alen, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_vec_mul_mpoly(fmpz_mpoly_struct * A, slong Alen, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

Divide each polynomial of $A$ exactly by $c$, respectively multiply each
by $c$.
:::

::: function
void fmpz_mpoly_unit_normalize(fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)

Negates the nonzero polynomial $A$ if its leading coefficient is
negative.
:::

## Partial fractions and Hensel lifting {#hensel}

These functions implement the multivariate Hensel lifting described in
[](algorithms/multivariate/hensel-lifting). Given $A \equiv \prod_i f_i
\bmod (x_m - \alpha)$ with the correct leading coefficients, the factors
are lifted one power of $x_m - \alpha$ at a time, each step solving a
multivariate partial fraction problem $t = \sum_i \delta_i \prod_{j \ne i}
f_j$ in the variables $x_0, \ldots, x_{m-1}$, recursively down to a
univariate problem over $\Z$.

::: type
fmpz_poly_pfrac_struct
fmpz_poly_pfrac_t

Precomputed data for univariate partial fraction problems over $\Z$ with
fixed denominators $b_0, \ldots, b_{r-1}$, solved modulo a power of a
prime and checked by division.
:::

::: function
void fmpz_poly_pfrac_init(fmpz_poly_pfrac_t Iv)
void fmpz_poly_pfrac_clear(fmpz_poly_pfrac_t Iv)

Initialise and free.
:::

::: function
int fmpz_poly_pfrac_precompute(fmpz_poly_pfrac_t Iv, const fmpz_poly_struct * b, slong r)

Prepares `Iv` for the $r \ge 2$ nonconstant denominators $b_i \in \Z[x]$:
computes the products $b_{i+1} \cdots b_{r-1}$, the resultants that bound
the size of solutions, a prime $p$ for which the $b_i$ stay coprime
modulo $p$, and the corresponding inverses. Returns 0 if $r < 2$, some
$b_i$ is constant, or two of them are not coprime.
:::

::: function
int fmpz_poly_pfrac_precomp(fmpz_poly_struct * c, const fmpz_poly_t A, fmpz_poly_pfrac_t Iv)

Solves $A / \prod_i b_i = \sum_i c_i / b_i$ for $c_i \in \Z[x]$ with
$\deg c_i < \deg b_i$ for $i < r - 1$. The solution is computed modulo a
power $p^k$ of the prime and checked by exact division; on failure the
precision is raised until a Cramer–Hadamard bound on the solution shows
that no integer solution exists. Returns 1 if a solution exists and 0
otherwise, so a return value of 0 is a proof.
:::

::: type
fmpz_mpoly_pfrac_struct
fmpz_mpoly_pfrac_t

Precomputed data for multivariate partial fraction problems with $r$
fixed denominators $\beta_0, \ldots, \beta_{r-1}$ in the variables $x_0,
\ldots, x_w$: their images $\beta_j^{(l)}$ under $x_k = \alpha_{k-1}$ for
$k > l$ (`mbetas`), the products $\prod_{j \ne i} \beta_j^{(l)}$ and their
Taylor expansions in $x_l - \alpha_{l-1}$, and the solutions `deltas` at
each level $l$.
:::

::: function
int fmpz_mpoly_pfrac_init(fmpz_mpoly_pfrac_t Iv, flint_bitcnt_t bits, slong r, slong w, const fmpz_mpoly_struct * betas, const fmpz * alpha, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_pfrac_clear(fmpz_mpoly_pfrac_t Iv, const fmpz_mpoly_ctx_t ctx)

Prepare for the $r$ denominators `betas` in the variables $x_0, \ldots,
x_w$ and the evaluation point $(\alpha_0, \ldots, \alpha_{w-1})$ for
$x_1, \ldots, x_w$, in a lex context with `bits` bits. The univariate
images in $x_0$ must keep their degrees and be pairwise coprime; if not,
the initialisation throws an error. Free the data.
:::

::: function
int fmpz_mpoly_pfrac(slong l, fmpz_mpoly_t t, const slong * degs, fmpz_mpoly_pfrac_t Iv, const fmpz_mpoly_ctx_t ctx)

Solves $t = \sum_i \delta_i \prod_{j \ne i} \beta_j^{(l)}$ in the variables
$x_0, \ldots, x_l$ with $\deg_{x_0} \delta_i < \deg_{x_0} \beta_i$ and
$\deg_{x_k} \delta_i$ at most what the bound `degs[k]` allows, storing
the solution in `Iv->deltas + l*r`. For $l = 0$ this is a univariate
problem; for $l > 0$ the Taylor coefficients of $t$ in $x_l - \alpha_{l-1}$
are processed in turn, each by a call at level $l - 1$. $t$ is
destroyed. Returns 1 if a solution exists, 0 if no solution satisfies the
degree bounds, and $-1$ if $t$ cannot be repacked to the working bits.
:::

::: function
int fmpz_mpoly_hlift(slong m, fmpz_mpoly_struct * f, slong r, const fmpz * alpha, const fmpz_mpoly_t A, const slong * degs, const fmpz_mpoly_ctx_t ctx)

Lifts the $r \ge 2$ factors $f_i$ of $A \equiv \prod_i f_i \bmod (x_m -
\alpha_{m-1})$, in a lex context with single-word exponents, to factors in
the variables $x_0, \ldots, x_m$. The leading coefficients of the $f_i$ in
$x_0$ must already be the correct ones (polynomials in all of $x_1,
\ldots, x_m$), and `degs` are the degrees of $A$. The lifting runs up to
degree `degs[m]` in $x_m - \alpha_{m-1}$; it uses the error update of the
quartic method, specialised for $r = 2$, when $r < 20$, and recomputes the
product of the factors (quintic method) otherwise. Returns 0 if a partial
fraction problem has no solution or the degrees exceed `degs[m]`, in which
case there is no factorisation of $A$ with these images and leading
coefficients. Only the quintic version checks at the end that $A =
\prod_i f_i$; the callers verify the factors by division in any case. See
[](algorithms/multivariate/hensel-lifting#hlift).
:::

::: function
void _fmpz_mpoly_get_lead0(fmpz_mpoly_t c, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx)
void _fmpz_mpoly_set_lead0(fmpz_mpoly_t A, const fmpz_mpoly_t B, const fmpz_mpoly_t c, const fmpz_mpoly_ctx_t ctx)

In a lex context, set $c$ to the leading coefficient of $A$ as a
polynomial in $x_0$ (a polynomial in the other variables), and set $A$
to $B$ with that leading coefficient replaced by $c$.
:::

## Dense bivariate polynomials {#bpoly}

::: type
fmpz_bpoly_struct
fmpz_bpoly_t

A dense bivariate polynomial $\sum_i c_i(y) x^i$ over $\Z$, stored as an
array `coeffs` of `length` univariate polynomials (`fmpz_poly_struct`) in
$y$, with room for `alloc`. The leading coefficient
`coeffs[length - 1]` is nonzero.
:::

::: function
void fmpz_bpoly_init(fmpz_bpoly_t A)
void fmpz_bpoly_clear(fmpz_bpoly_t A)
void fmpz_bpoly_swap(fmpz_bpoly_t A, fmpz_bpoly_t B)
void fmpz_bpoly_realloc(fmpz_bpoly_t A, slong len)
void fmpz_bpoly_fit_length(fmpz_bpoly_t A, slong len)
void fmpz_bpoly_zero(fmpz_bpoly_t A)

Memory management: initialise as zero, free, swap in constant time,
ensure room for at least `len` coefficients in $x$ (`realloc` never
shrinks), and set to zero.
:::

::: function
fmpz_poly_struct * fmpz_bpoly_lead(fmpz_bpoly_t A)
slong fmpz_bpoly_degree0(const fmpz_bpoly_t A)
slong fmpz_bpoly_degree1(const fmpz_bpoly_t A)

Return a pointer to the leading coefficient in $x$, the degree in $x$,
and the degree in $y$ (both $-1$ for zero).
:::

::: function
void fmpz_bpoly_set_coeff(fmpz_bpoly_t A, slong exp0, slong exp1, const fmpz_t c)

Sets the coefficient of $x^{\mathit{exp0}} y^{\mathit{exp1}}$ to the
nonzero integer $c$.
:::

::: function
void fmpz_bpoly_print_pretty(fmpz_bpoly_t A, const char * var0, const char * var1)

Prints $A$ as a sum of terms `(c(y))*x^i`, with the given names for $x$
and $y$.
:::

::: function
void fmpz_mpoly_get_bpoly(fmpz_bpoly_t A, const fmpz_mpoly_t B, slong var0, slong var1, const fmpz_mpoly_ctx_t ctx)
void fmpz_mpoly_set_fmpz_bpoly(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_bpoly_t B, slong var0, slong var1, const fmpz_mpoly_ctx_t ctx)

Convert between a multivariate polynomial in the variables of index
`var0` (for $x$) and `var1` (for $y$) and a dense bivariate polynomial.
The first requires single-word exponents and ignores the exponents of the
other variables; the second packs the result into `Abits` bits.
:::

::: function
int fmpz_bpoly_resultant(fmpz_poly_t res, const fmpz_bpoly_t A, const fmpz_bpoly_t B, int proved)
int _fmpz_bpoly_resultant(fmpz_poly_t res, fmpz_poly_struct * A, slong lenA, fmpz_poly_struct * B, slong lenB, int proved)
int _fmpz_bpoly_resultant_len(fmpz_poly_t res, fmpz_poly_struct * A, slong lenA, fmpz_poly_struct * B, slong lenB, slong rlen, int proved)

Set `res` to the resultant of $A$ and $B$ with respect to $x$, a
polynomial in $y$, by reduction modulo word-size primes, resultants of
the images (by `_n_bpoly_mod_resultant_len`, with whichever algorithm it
picks) and Chinese remaindering. If `proved` is nonzero, enough primes
are used for their product to exceed twice a Hadamard-type bound on the
coefficients of the resultant, which makes the result exact; otherwise the
reconstruction stops once it has been left unchanged by primes whose
product exceeds $2^{100}$, which is a heuristic. Return 0 if no image could
be computed. The underscore version works in place on the coefficient
arrays, whose contents it removes, and requires $\mathit{lenA} \ge
\mathit{lenB} \ge 2$; `_fmpz_bpoly_resultant_len` is the same given that
the resultant has length at most `rlen`, which saves evaluation points
when that is less than the bound read off the degrees. This is the dense
case of `fmpz_mpoly_resultant`
and the modular resultant of the generic polynomial module; see
[](algorithms/multivariate/modp-bivariate-resultant).
:::

::: type
fmpz_tpoly_struct
fmpz_tpoly_t

A growable vector of `fmpz_bpoly_struct` (fields `coeffs`, `alloc`,
`length`), used for lists of bivariate factors.
:::

::: function
void fmpz_tpoly_init(fmpz_tpoly_t A)
void fmpz_tpoly_clear(fmpz_tpoly_t A)
void fmpz_tpoly_swap(fmpz_tpoly_t A, fmpz_tpoly_t B)
void fmpz_tpoly_fit_length(fmpz_tpoly_t A, slong len)

Initialise an empty vector, free it, swap in constant time, and ensure
room for `len` entries.
:::

::: function
void fmpz_bpoly_factor(fmpz_poly_t c, fmpz_tpoly_t F, fmpz_bpoly_t B)

Factors the squarefree bivariate polynomial $B$ of positive degree in $x$:
sets $c$ to its content as a polynomial in $x$ over $\Z[y]$ and $F$ to the
irreducible factors of the primitive part. The point $\alpha$ is the first
of $0, 1, -1, 2, -2, \ldots$ for which $B(x, \alpha)$ keeps its degree and
is squarefree; $B(x, \alpha)$ is factored with `fmpz_poly_factor`, the
factors are lifted in $(\Z/p^k\Z)[[y - \alpha]][x]$ to a precision given by
a bound on the factors, and recombined by trying subsets. $B$ is
destroyed. If $B$ is not squarefree, the search for $\alpha$ does not
terminate.
:::

::: function
int fmpz_bpoly_factor_ordered(fmpz_poly_t c, fmpz_tpoly_t F, fmpz_bpoly_t B, const fmpz_t alpha, const fmpz_poly_factor_t Bevalf)

As `fmpz_bpoly_factor`, with the point $\alpha$ and the squarefree
factorisation `Bevalf` of $B(x, \alpha)$ given (with at least two
factors). On success there is one factor in $F$ for each factor of
`Bevalf`, in the same order. Returns 1 on success, 0 if the lifted
factors do not correspond one to one with the factors of `Bevalf`, and
$-1$ if $B$ is not primitive (with its content in $c$). Used by Kaltofen's
leading coefficient method.
:::

## Irreducible factorisation internals {#internal}

::: function
int _fmpz_mpoly_factor_squarefree(fmpz_mpoly_factor_t f, fmpz_mpoly_t A, const fmpz_t e, const fmpz_mpoly_ctx_t ctx)

Appends a squarefree factorisation of $A^e$ to `f`, where $A$
is primitive with respect to each variable. Constant and monomial $A$ are
handled directly; otherwise Yun's algorithm [@Yun1976] is run with respect
to the first variable in which $A$ has a nonzero derivative. $A$ is
destroyed. Returns 0 if a GCD failed for every variable.
:::

::: function
int fmpz_mpoly_factor_irred_zassenhaus(fmpz_mpolyv_t fac, const fmpz_mpoly_t A, const fmpz_mpoly_ctx_t ctx, zassenhaus_prune_t Z)
int fmpz_mpoly_factor_irred_wang(fmpz_mpolyv_t fac, const fmpz_mpoly_t A, const fmpz_mpoly_factor_t lcAfac, int lcAfac_irred, const fmpz_mpoly_t lcA, const fmpz_mpoly_ctx_t ctx, flint_rand_t state, zassenhaus_prune_t Z, int allow_shift)
int fmpz_mpoly_factor_irred_zippel(fmpz_mpolyv_t fac, const fmpz_mpoly_t A, const fmpz_mpoly_factor_t lcAfac, int lcAfac_irred, const fmpz_mpoly_t lcA, const fmpz_mpoly_ctx_t ctx, flint_rand_t state, zassenhaus_prune_t Z)

Set `fac` to the irreducible factors of $A$, which must be squarefree,
primitive with respect to each variable and have positive leading
coefficient, in a lex context with at least two variables. `lcA` is the
leading coefficient of $A$ in $x_0$ and `lcAfac` a squarefree
factorisation of it, irreducible if `lcAfac_irred` is nonzero. `Z`
collects the degrees of the factors of univariate images and can prove
irreducibility early.

The Wang version evaluates $x_1, \ldots$ at $0$ and, if `allow_shift` is
nonzero, then at random points of growing size, distributes the leading
coefficient with `fmpz_mpoly_factor_lcc_wang` or
`fmpz_mpoly_factor_lcc_kaltofen`, and lifts the factors one variable at a
time with `fmpz_mpoly_hlift`. The Zippel version lifts modulo a prime with
sparse Hensel lifting and then $p$-adically. The Zassenhaus version lifts
without leading coefficient information and recombines subsets of the
lifted factors. Return 1 on success, 0 if the method failed (another
method or point may succeed), and $-1$ on an exceptional failure such as
exponents that no longer fit one word.
:::

::: function
int fmpz_mpoly_factor_lcc_wang(fmpz_mpoly_struct * lc_divs, const fmpz_mpoly_factor_t lcAfac, const fmpz_t Auc, const fmpz_poly_struct * Auf, slong r, const fmpz * alpha, const fmpz_mpoly_ctx_t ctx)

Wang's leading coefficient determination [@Wan1978]. Given the
irreducible factorisation `lcAfac` of the leading coefficient of $A$ in
$x_0$ and the factorisation $\mathit{Auc} \prod_{i < r} \mathit{Auf}_i$ of
$A(x_0, \alpha)$, sets `lc_divs[i]` to the predicted leading coefficient
of the $i$-th factor of $A$: the product of the powers of the factors of
the leading coefficient whose values at $\alpha$ divide the leading
coefficient of $\mathit{Auf}_i$, times an integer. Returns 0 if the
evaluation point fails Wang's condition, that each factor of the leading
coefficient evaluates to an integer with a prime divisor not dividing
$\mathit{Auc}$ times the constant and the values of the other factors.
:::

::: function
int fmpz_mpoly_factor_lcc_kaltofen(fmpz_mpoly_struct * divs, const fmpz_mpoly_factor_t lcAf_, const fmpz_mpoly_t A, slong r, const fmpz * alpha, slong * degs, const fmpz_poly_factor_t uf, const fmpz_mpoly_ctx_t ctx)

Kaltofen's leading coefficient determination [@Kal1985]. For each
variable $x_v$, $v \ge 1$, the bivariate image of $A$ in $x_0, x_v$ (all
other variables evaluated at $\alpha$) is factored with
`fmpz_bpoly_factor_ordered` in the order of the univariate factors `uf`,
and the leading coefficients of its factors show how the squarefree
factors `lcAf_` of the leading coefficient of $A$ involving $x_v$ are
distributed; these are multiplied into `divs`. Returns 1 if the whole
leading coefficient has been distributed, 0 if only part of it was, and
$-1$ if the evaluation point is unusable.
:::

::: function
int fmpz_mpoly_factor_lcc_kaltofen_step(fmpz_mpoly_struct * divs, slong r, fmpz_mpoly_factor_t Af, const fmpz_poly_struct * Au, slong v, const fmpz * alphas, const fmpz_mpoly_ctx_t ctx)

One step of `fmpz_mpoly_factor_lcc_kaltofen` for the variable $x_v$:
given the squarefree factorisation `Af` of the remaining part of the
leading coefficient and the univariate leading coefficients `Au[i]` in
$x_v$ of the bivariate factors, lifts the matching factors of the image of
`Af` to all variables, multiplies them into `divs[i]` and removes them from
`Af`. Returns 1 on success and 0 if the images do not match.
:::

::: function
int fmpz_mpoly_evaluate_rest_except_one(fmpz_poly_t e, const fmpz_mpoly_t A, const fmpz * alphas, slong v, const fmpz_mpoly_ctx_t ctx)

Sets $e$ to $A$ with $x_i = \alpha_{i-1}$ for every $i \ge 1$ other than
$v$, as a polynomial in $x_v$. Returns 0 if an evaluation failed or the
result involves $x_0$.
:::

::: function
void _fmpz_mpoly_eval_rest_to_poly(fmpz_poly_t E, const fmpz_mpoly_t A, const fmpz * alphas, const fmpz_mpoly_ctx_t ctx)
int _fmpz_mpoly_evaluate_rest_fmpz(fmpz * E, slong * starts, slong * ends, slong * stops, ulong * es, const fmpz * Acoeffs, const ulong * Aexps, slong Alen, slong var, const fmpz * alphas, const slong * offsets, const slong * shifts, slong N, ulong mask, slong nvars)

The first sets $E = A(x_0, \alpha_0, \ldots, \alpha_{n-2})$ for $A$ in a
lex context with single-word exponents. The second is its kernel: for a
block of `Alen` terms of a lex polynomial whose exponents agree in the
variables before `var`, it evaluates the variables from `var` on at
`alphas` (indexed from `var`) by a recursive Horner scheme and puts the
value in `E[0]`; `starts`, `ends`, `stops`, `es` and the further entries
of $E$ are work arrays of length `nvars`, and `offsets` and `shifts`
locate the exponent fields.
:::

::: function
void fmpz_mpoly_compression_do(fmpz_mpoly_t L, const fmpz_mpoly_ctx_t Lctx, fmpz * Acoeffs, slong Alen, mpoly_compression_t M)
void fmpz_mpoly_compression_undo(fmpz_mpoly_t A, flint_bitcnt_t Abits, const fmpz_mpoly_ctx_t Actx, fmpz_mpoly_t L, const fmpz_mpoly_ctx_t Lctx, mpoly_compression_t M)

Apply the change of variables found by `mpoly_compression_set`. The first
builds $L$ in the compressed lex context `Lctx` from the coefficients
`Acoeffs`, which are swapped in, and the compressed exponents stored in
`M`, then sorts it and makes its leading coefficient positive. The second
maps $L$ (typically a factor of the compressed polynomial, whose
coefficients are swapped out) back to `Actx` with `Abits` bits, divides by
the largest monomial dividing it, sorts, and makes the leading coefficient
positive.
:::

::: function
void tuple_print(fmpz * alpha, slong n)
void tuple_next(fmpz * alpha, slong n)
void tuple_saturate(fmpz * alpha, slong n, slong m)

Helpers for enumerating evaluation points in $\N^n$: print a tuple, step
to the next tuple (all tuples with a given sum are visited before the sum
increases), and move the entries after index $m$ into entry $m$ (making it
1 if the first $m + 1$ entries would all be zero) so that one of the first
$m + 1$ entries changes at the next step.
:::
