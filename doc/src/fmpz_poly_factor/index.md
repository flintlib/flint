---
title: fmpz_poly_factor.h
header: fmpz_poly_factor.h
description: Factorisation of polynomials over the integers into irreducible factors, by squarefree factorisation, Hensel lifting and the Zassenhaus and van Hoeij recombination algorithms.
category: Univariate polynomials
order: 2
---

This module factors polynomials in $\Z[x]$ into irreducible factors.
A factorisation is stored as an `fmpz_poly_factor_t`: a signed
content $c$ together with distinct primitive irreducible factors $f_i$
of positive degree and positive leading coefficient, and multiplicities
$e_i$, so that $f = c \prod_i f_i^{e_i}$. The main function is
`fmpz_poly_factor`, which lifts a factorisation modulo a small prime
$p$-adically and recombines the modular factors by the Zassenhaus
algorithm when there are few and by van Hoeij's lattice reduction
algorithm when there are many. The guide
[](fmpz_poly_factor/factor) shows how to use it, and
[](algorithms/polynomials/factoring-integers) describes the
algorithms.

::: summary algorithms/polynomials/factoring-integers :::

## Types and memory management {#memory}

::: type
fmpz_poly_factor_struct
fmpz_poly_factor_t

A factorisation $c \prod_{i < \mathit{num}} p_i^{e_i}$, with fields `c`
(an `fmpz`, the constant factor), `p` (an array of `alloc`
polynomials, of which the first `num` are used), `exp` (the
multiplicities) and `num`, `alloc`. An `fmpz_poly_factor_t` is an
array of one structure. A freshly initialised structure represents the
empty product, with $c = 1$.
:::

::: function
void fmpz_poly_factor_init(fmpz_poly_factor_t fac)
void fmpz_poly_factor_init2(fmpz_poly_factor_t fac, slong alloc)

Initialise `fac` as the empty factorisation with $c = 1$, allocating
nothing, respectively room for `alloc` factors.
:::

::: function
void fmpz_poly_factor_realloc(fmpz_poly_factor_t fac, slong alloc)

Changes the allocated space to exactly `alloc` factors; for `alloc` =
0 the structure is cleared and reinitialised. The number of factors is
not adjusted, so `alloc` should be at least `fac->num`.
:::

::: function
void fmpz_poly_factor_fit_length(fmpz_poly_factor_t fac, slong len)

Ensures that `fac` has room for at least `len` factors, growing the
allocation to at least twice its old size if it has to grow.
:::

::: function
void fmpz_poly_factor_clear(fmpz_poly_factor_t fac)

Frees the memory used by `fac`.
:::

## Manipulating factorisations {#manipulation}

::: function
void fmpz_poly_factor_set(fmpz_poly_factor_t res, const fmpz_poly_factor_t fac)

Sets `res` to a copy of `fac`.
:::

::: function
void fmpz_poly_factor_insert(fmpz_poly_factor_t fac, const fmpz_poly_t p, slong e)

Multiplies the factorisation by $p^e$: if $p$ is equal to one of the
stored factors, its multiplicity is increased by $e$, otherwise $p$ is
appended with multiplicity $e$. The content is not changed. Meant for
primitive $p$ of positive degree and $e \ne 0$, which is not checked.
:::

::: function
void fmpz_poly_factor_concat(fmpz_poly_factor_t res, const fmpz_poly_factor_t fac)

Multiplies `res` by the factorisation `fac`: the contents are
multiplied and the factors of `fac` are inserted with
`fmpz_poly_factor_insert`. `res` and `fac` must be distinct.
:::

::: function
void fmpz_poly_factor_get_fmpz(fmpz_t z, const fmpz_poly_factor_t F)
void fmpz_poly_factor_get_fmpz_poly(fmpz_poly_t z, const fmpz_poly_factor_t F, slong i)

Set $z$ to the constant factor $c$ of `F`, respectively to the factor
with index $i$, which must satisfy $0 \le i <$ `F->num`.
:::

::: function
void fmpz_poly_factor_print(const fmpz_poly_factor_t fac)

Prints the content on one line, then each factor in the plain format
of `fmpz_poly_print` followed by `^` and its multiplicity, one per
line.
:::

## Factoring {#factoring}

::: function
void fmpz_poly_factor(fmpz_poly_factor_t fac, const fmpz_poly_t G)

Sets `fac` to the factorisation of $G$ into its signed content and
irreducible primitive factors with positive leading coefficient and
their multiplicities. The zero polynomial gives content 0 and no
factors, a constant $c$ gives content $c$ and no factors. `fac` must
be initialised; its previous factors are discarded. The factors are
not sorted.

::: details
Factors $x^k$ are split off, and if $G = g(x^d)$ with $d > 1$, the
deflated $g$ is factored first and each factor $h$ of $g$ is inflated
again: $h(x^p)$ for a prime $p \mid d$ is either certified irreducible
by `_fmpz_poly_factor_inflation_is_irreducible_capelli` or factored
again at degree $p \deg h$, which is much smaller than refactoring
$G$ itself. Otherwise the content is removed and
`fmpz_poly_factor_squarefree` splits the polynomial into squarefree
parts, each of which is factored by `_fmpz_poly_factor_zassenhaus`
with van Hoeij's algorithm enabled when there are more than 8 modular
factors. Quadratics and cubics are factored by special methods.
:::
:::

::: function
void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t fac, const fmpz_poly_t G)

Factors $G$ like `fmpz_poly_factor`, but without deflation and using
only the Zassenhaus recombination. Its running time is exponential in
the number of modular factors, so it is meant for comparison and
testing. `fac` should be freshly initialised, since the factors are
appended to those already present.
:::

::: function
void fmpz_poly_factor_squarefree(fmpz_poly_factor_t fac, const fmpz_poly_t F)

Sets the content of `fac` to the signed content $c$ of $F$ (with the
sign of the leading coefficient) and appends squarefree, pairwise
coprime, primitive factors $g_i$ of positive degree with distinct
multiplicities $e_i$, so that $F = c \prod_i g_i^{e_i}$ with
$\gcd(g_i, g_i') = 1$. Each $g_i$ is the product of all irreducible
factors of $F$ of multiplicity $e_i$. Uses Yun's algorithm
[@Yun1976] ([](algorithms/polynomials/squarefree)). `fac` should be
freshly initialised.
:::

::: function
void _fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, slong exp, const fmpz_poly_t f, slong cutoff, int use_van_hoeij)

Appends the irreducible factors of the squarefree, primitive
polynomial $f$ of length at least 2 with nonzero constant coefficient
and positive leading coefficient, each with multiplicity `exp`, to
`final_fac`; none of the conditions is checked.

::: details
Three primes $p$ for which $f$ is squarefree modulo $p$ with the same
degree and nonzero constant coefficient are tried, starting from 2,
and the one with fewest modular factors is chosen; the factor degrees
seen modulo all three primes are combined to rule out impossible
degrees of true factors. With one modular factor, $f$ is irreducible.
If the number $r$ of modular factors exceeds `cutoff` and
`use_van_hoeij` is nonzero, `fmpz_poly_factor_van_hoeij` is called;
otherwise the factors are lifted to the precision given by the
Mignotte bound and recombined by
`fmpz_poly_factor_zassenhaus_recombination_with_prune`. Quadratics and
cubics are passed to the special functions.
:::
:::

::: function
void fmpz_poly_factor_van_hoeij(fmpz_poly_factor_t final_fac, const nmod_poly_factor_t fac, const fmpz_poly_t f, slong exp, ulong p)

Appends the irreducible factors of the squarefree primitive $f$, each
with multiplicity `exp`, to `final_fac`, given its factorisation `fac`
into $r \ge 2$ monic irreducible factors modulo the prime $p$, by van
Hoeij's knapsack algorithm with gradual feeding of CLD data
[@vHo2002], [@HvHN2011]. The modular factors are lifted to an initial
precision, a lattice is reduced with LLL after each new column of CLD
data, and the precision is doubled until a 0-1 basis of the
recombination lattice is found and verified.
:::

::: function
void _fmpz_poly_factor_quadratic(fmpz_poly_factor_t fac, const fmpz_poly_t f, slong exp)
void _fmpz_poly_factor_cubic(fmpz_poly_factor_t fac, const fmpz_poly_t f, slong exp)

Insert the factorisation of the primitive quadratic (respectively
cubic) polynomial $f$ into `fac`, with multiplicity `exp`, without
changing the content of `fac`. A quadratic splits if and only if its
discriminant is a square; a cubic is reducible if and only if it has a
rational root, which is located from real and 2-adic approximations of
the roots and tested over $\Z$
([](algorithms/polynomials/factoring-integers#cubic)). The content of $f$ must be 1, which is
not checked.
:::

::: function
void fmpz_poly_factor_mignotte(fmpz_t B, const fmpz_poly_t f)

Sets $B$ to a bound for the absolute values of the coefficients of all
factors in $\Z[x]$ of $f$ of degree $m \ge 2$: the maximum of $|f(0)|$,
$|\operatorname{lc}(f)|$ and of $\binom{m-1}{j} N + \binom{m-1}{j-1}
|\operatorname{lc}(f)|$ for $1 \le j < m$, where $N \ge \|f\|_2$ is
an integer [@Coh1996, Theorem 3.5.1]. The Zassenhaus algorithm lifts to
$p^a > 2 |\operatorname{lc}(f)| B$.
:::

## Irreducibility of inflations {#capelli}

::: function
int _fmpz_poly_factor_inflation_is_irreducible_capelli(const fmpz_poly_t T, ulong p)

Given $T$ irreducible over $\mathbb{Q}$ of degree at least 1 and a prime
$p$, attempts to prove that $T(x^p)$ is irreducible over $\mathbb{Q}$.
Returns 1 if a proof was found and 0 if the search was inconclusive,
in which case $T(x^p)$ may still be irreducible. The behaviour is
undefined if $T$ is reducible. $T$ need not be primitive and may have
negative leading coefficient.

::: details
By Capelli's theorem, for prime $p$ the polynomial $T(x^p)$ is
irreducible over $\mathbb{Q}$ if and only if $\theta$, the image of $x$
in $K = \mathbb{Q}[x]/(T)$, is not a $p$-th power in $K$. For $\deg T =
1$, and for $\deg T = 2$ with $p = 2$, this is decided exactly (the
return value is 1 precisely when $T(x^p)$ is irreducible). Otherwise the
function searches for a rational prime $q$, with $T$ squarefree modulo
$q$ and $q$ dividing neither the leading nor the constant coefficient
of $T$, such that $\theta$ reduces to a non-$p$-th power in some
residue field of $K$ above $q$, which is detected with Euler's
criterion $\theta^{(q^f - 1)/p} \ne 1$ in a residue field of degree
$f$; any such witness proves $\theta \notin K^p$. The search budget is
linear in $\deg T$ (for inputs with abelian splitting fields of
exponent $p$, such as deflations of Swinnerton-Dyer polynomials,
witnesses have density only $1/(p \deg T)$ among the primes and occur
only in degree-1 residue fields), and the search stops early when the
accumulated local evidence makes irreducibility unlikely.
:::
:::

::: function
int fmpz_poly_factor_inflation_is_irreducible_capelli(const fmpz_poly_t T, ulong d)

Given $T$ irreducible over $\mathbb{Q}$ of degree at least 1 and $d \ge
1$, returns 1 if $T(x^d)$ is proved irreducible over $\mathbb{Q}$ and
0 if the attempt was inconclusive. The prime-step certificate is
chained over the prime factorisation of $d$, using $T(x^d) =
U(x^{d/p})$ with $U = T(x^p)$. The behaviour is undefined if $T$ is
reducible. `fmpz_poly_factor` uses this to avoid refactoring
inflations of irreducible factors of deflated polynomials.
:::

## Recombination {#recombination}

::: function
void fmpz_poly_factor_zassenhaus_recombination(fmpz_poly_factor_t final_fac, const fmpz_poly_factor_t lifted_fac, const fmpz_poly_t F, const fmpz_t P, slong exp)
void fmpz_poly_factor_zassenhaus_recombination_with_prune(fmpz_poly_factor_t final_fac, const fmpz_poly_factor_t lifted_fac, const fmpz_poly_t F, const fmpz_t P, slong exp, const zassenhaus_prune_t Z)

Find the irreducible factors of the squarefree primitive $F$ from its
monic factors `lifted_fac` modulo $P = p^a$, where $P > 2
|\operatorname{lc}(F)| B$ for the Mignotte bound $B$, by trying subsets
of the modular factors in order of increasing size: for each subset
$S$, the product $\operatorname{lc}(F) \prod_{i \in S} g_i$ is reduced
into the symmetric range modulo $P$, its primitive part is tested as a
divisor of $F$, and true factors are removed together with their
modular factors. Each factor is appended to `final_fac` with
multiplicity `exp`. The `with_prune` version skips subsets whose total
degree is impossible according to `Z`.
:::

::: function
slong _fmpz_poly_factor_CLD_mat(fmpz_mat_t res, const fmpz_poly_t f, fmpz_poly_factor_t lifted_fac, fmpz_t P, ulong k)

Fills the $(r + 1) \times 2k$ matrix `res` with the CLD data of van
Hoeij's algorithm for the $r$ lifted factors $g_i$ of $f$ modulo $P$:
row $i$ holds low and high coefficients of $f g_i' / g_i$ modulo $P$
in the symmetric range, and row $r$ the corresponding bounds from
`fmpz_poly_CLD_bound`. Of the low $k$ and high $k$ coefficients, only
those whose bound times $\sqrt{\deg f}$ is small enough compared to
$P$ are used; the function returns the number of columns filled, low
ones first.
:::

::: function
int fmpz_poly_factor_van_hoeij_check_if_solved(fmpz_mat_t M, fmpz_poly_factor_t final_fac, fmpz_poly_factor_t lifted_fac, const fmpz_poly_t f, fmpz_t P, slong exp, fmpz_t lc)

Checks whether the reduced lattice basis $M$ of van Hoeij's algorithm
determines the factorisation: if the first $r$ columns of $M$ define a
partition of the $r$ lifted factors (a 0-1 basis), the products of
the factors in each part (times `lc`, reduced modulo $P$, made
primitive) are tested as divisors of $f$, first modulo 2 and then over
$\Z$. If all of them divide, the factors are appended to `final_fac`
with multiplicity `exp` and 1 is returned; otherwise 0 is returned and
`final_fac` is unchanged.
:::

## Internal functions {#internal}

These functions support the Zassenhaus recombination. A subset of the
$r$ modular factors is encoded in an array `s` of length $r$ whose
entries are factor indices $i \ge 0$ for the factors in the subset and
$-i - 1$ for the others.

::: function
void zassenhaus_subset_first(slong * s, slong r, slong m)

Sets `s` to the first subset of size $m$ of its $r$ entries, the first
$m$ entries.
:::

::: function
int zassenhaus_subset_next(slong * s, slong r)

Advances `s` to the next subset of the same size in a fixed order and
returns 1, or returns 0 if it was the last.
:::

::: function
slong zassenhaus_subset_next_disjoint(slong * s, slong r)

Removes the entries of the current subset from `s`, moving the
remaining entries to the front, and sets it to the next subset of the
same size $m$ of the remaining entries, returning 1, or returns 0 if
there is none.
:::

::: type
zassenhaus_prune_struct
zassenhaus_prune_t

Records the degrees $d$ for which a factor of degree $d$ is still
possible, given the degrees of the modular factors modulo several
primes: a true factor of degree $d$ is a product of modular factors of
total degree $d$ modulo every prime.
:::

::: function
void zassenhaus_prune_init(zassenhaus_prune_t Z)
void zassenhaus_prune_clear(zassenhaus_prune_t Z)

Initialise and free `Z`.
:::

::: function
void zassenhaus_prune_set_degree(zassenhaus_prune_t Z, slong d)

Sets the degree of the polynomial to be factored to $d \ge 1$ and marks
all degrees $0, \ldots, d$ as possible.
:::

::: function
void zassenhaus_prune_start_add_factors(zassenhaus_prune_t Z)
void zassenhaus_prune_add_factor(zassenhaus_prune_t Z, slong deg, slong exp)
void zassenhaus_prune_end_add_factors(zassenhaus_prune_t Z)

Add the degrees of the modular factors for one prime: after
`zassenhaus_prune_start_add_factors`, each call to
`zassenhaus_prune_add_factor` records `exp` factors of degree `deg`,
and `zassenhaus_prune_end_add_factors` intersects the set of possible
degrees with the set of subset sums of the recorded degrees. The
degrees must add up to the degree given to
`zassenhaus_prune_set_degree`; otherwise an exception is raised.
:::

::: function
int zassenhaus_prune_degree_is_possible(const zassenhaus_prune_t Z, slong d)

Returns whether a factor of degree $d$ is still possible.
:::

::: function
int zassenhaus_prune_must_be_irreducible(const zassenhaus_prune_t Z)

Returns whether no degree strictly between 0 and the degree of the
polynomial is possible, which proves the polynomial irreducible.
:::
