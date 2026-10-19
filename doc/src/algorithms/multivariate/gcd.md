---
title: Greatest common divisors
summary: How FLINT computes GCDs of multivariate polynomials over the integers: reductions that make most inputs easy (monomial content, deflation, missing variables, univariate cases), degree bounds from random images, and the four modular algorithms it chooses between (Brown's dense algorithm, Zippel's sparse algorithm, the Berlekamp–Massey based algorithm of Hu and Monagan, and EEZ-GCD Hensel lifting), with their correctness proofs and the cost model behind the choice.
complexity: Brown: polynomial in the dense sizes $\prod_i (d_i + 1)$ of the inputs and the coefficient size. Zippel and Hu–Monagan: polynomial in the number of terms of the GCD and cofactors, the number of variables, the degrees and the coefficient size, for inputs without unlucky evaluations. Hensel lifting: like dense multivariate Hensel lifting of two factors.
order: 40
---

$\Z[x_0, \ldots, x_{n-1}]$ is a unique factorisation domain, so two
polynomials $A$, $B$ have a greatest common divisor $G$, unique up to
sign; FLINT returns the one with positive leading coefficient, and also
the cofactors $\bar A = A/G$ and $\bar B = B/G$ if requested. Euclid's
algorithm with pseudo-remainders works in principle, but the
coefficients of the intermediate remainders grow and their supports fill
up even when $G$ and the cofactors are small and sparse. All competitive
algorithms are therefore *modular*: they compute images of $G$ modulo
primes and at evaluation points, where the problem is univariate or
bivariate over a finite field, and reconstruct $G$ from these images
[@Bro1971; @GCL1992, Ch. 7; @GG2003, §6]. The difficulty is to
know the shape of $G$ in advance or to detect the images that do not
reflect it.

The algorithms below are those used for `fmpz_mpoly` (and through it
`fmpq_mpoly`); the function-level description, with an example, is
[](fmpz_mpoly/gcd). The algorithms over finite fields are the same
inside, with the additional problem that a small field may not have
enough evaluation points, which is treated in
[](algorithms/multivariate/modp-gcd). The reconstruction of sparse
polynomials from values is described in [](algorithms/multivariate/sparse-interpolation)
and the lifting of factorisations on
[](algorithms/multivariate/hensel-lifting). There is no generic
multivariate GCD in `gr_mpoly`; the `gr` contexts of these modules call
the specialised functions.

## Contents and normalisation {#prelim}

Write $R = \Z[x_1, \ldots, x_{n-1}]$, so that $A \in R[x_0]$. The
*content* $\operatorname{cont}_{x_0}(A) \in R$ is the GCD of the
coefficients of $A$ as a polynomial in $x_0$, and $A$ is *primitive* if
this content is 1. By Gauss's lemma products of primitive polynomials are
primitive, so
$$
\gcd(A, B) = \gcd(\operatorname{cont}_{x_0} A, \operatorname{cont}_{x_0}
B) \cdot \gcd(\operatorname{pp}_{x_0} A, \operatorname{pp}_{x_0} B),
$$
and the algorithms below assume primitive inputs, computing contents by
recursive GCDs in fewer variables (`fmpz_mpolyl_content`). For an integer
content this is the usual $\gcd(\operatorname{cont}(A),
\operatorname{cont}(B))$ in $\Z$.

::: lemma Leading coefficient {#lem-gamma}
Let $A$, $B$ be primitive with respect to $x_0$ and $G = \gcd(A, B)$. Then
$\operatorname{lc}_{x_0}(G)$ divides $\gamma = \gcd(\operatorname{lc}_{x_0}(A),
\operatorname{lc}_{x_0}(B)) \in R$, and $\tilde G = (\gamma /
\operatorname{lc}_{x_0}(G)) \, G$ is the unique multiple of $G$ by an
element of $R$ with leading coefficient $\gamma$. It satisfies $\tilde G
\bar A = \gamma A$ and $\tilde G \bar B = \gamma B$ for $\bar A = A/G$,
$\bar B = B/G$.
:::

::: proof
$G \mid A$ gives $\operatorname{lc}(G) \mid \operatorname{lc}(A)$, and
similarly for $B$, so $\operatorname{lc}(G) \mid \gamma$. Multiplying $G
\bar A = A$ by $\gamma / \operatorname{lc}(G) \in R$ gives $\tilde G \bar A
= \gamma A$, and similarly for $B$. Two multiples of $G$ by elements of $R$
with the same leading coefficient are equal, since $R$ is a domain.
:::

The images of $\tilde G$ modulo a prime or at an evaluation point are
determined up to the normalisation of the leading coefficient, which is
known: this is how images computed independently are made consistent
before they are combined. At the end, $G$ is the primitive part of
$\tilde G$.

The same reasoning applies with several main variables, which gives a
common frame for the algorithms below [@Sch-mpoly]. View $A$ and $B$ in
$R'[x_{m+1}, \ldots, x_{n-1}]$ with $R' = D[x_0, \ldots, x_m]$, make them
primitive over $D[x_{m+1}, \ldots]$, compute GCDs in $m + 1$ variables of
images at evaluation points of the other variables, impose the image of
$\gamma = \gcd(\operatorname{lc}(A), \operatorname{lc}(B))$ (leading
coefficients with respect to $x_0, \ldots, x_m$) as their leading
coefficient, and interpolate. If the primitive part of the interpolated
polynomial divides $A$ and $B$, it is the GCD. Brown's algorithm takes
$m = n - 2$, with dense interpolation in the last variable and recursion
down to the Euclidean algorithm; Zippel's algorithm takes $m = 0$, and
the algorithm of Hu and Monagan $m = 1$, both with sparse
interpolation.

## Easy cases {#easy}

Before any modular computation, `fmpz_mpoly_gcd` removes as much
structure as possible; the list of steps is described in [](fmpz_mpoly/gcd#easy).
They rest on three facts.

- *Monomials*: the GCD of the monomial contents (term-wise minima of the
  exponents) can be split off, and $\gcd(c x^e, B)$ is $\gcd(c,
  \operatorname{cont}(B)) \cdot \min(x^e, \operatorname{mcont}(B))$.
- *Deflation*: by [](algorithms/multivariate/representation#lem-inflate),
  exponents can be divided by common strides without changing the GCD,
  which reduces degrees, and inputs can then be repacked into one word.
- *Missing variables*: if $x_v$ occurs in $A$ but not in $B$, then
  $\gcd(A, B) = \gcd(\operatorname{cont}_{x_v} A, B)$, because every
  divisor of $B$ is free of $x_v$ and divides $A$ if and only if it
  divides all coefficients of $A$ with respect to $x_v$.

With one common variable the GCD is univariate (`fmpz_poly_gcd`, see
[](algorithms/polynomials/gcd-integers)), and when one input is a
monomial times a univariate polynomial of degree at most 100 in some
variable, the GCD is computed with the coefficients of the other input
with respect to the remaining variables.

## Degree bounds {#bounds}

::: lemma Degrees of images {#lem-image-deg}
Let $\pi$ be a ring homomorphism from $\Z[x_0, \ldots, x_{n-1}]$ to
$\F_p[x_v]$ that reduces modulo $p$ and evaluates all variables except
$x_v$, and assume $\deg_{x_v} \pi(A) = \deg_{x_v} A$. Then $\pi(G)$ divides
$\gcd(\pi(A), \pi(B))$ and $\deg_{x_v} G \le \deg \gcd(\pi(A), \pi(B))$.
Equality holds unless $\pi$ maps a certain nonzero subresultant of $A/G$
and $B/G$ (with respect to $x_v$) to zero.
:::

::: proof
$\pi(A) = \pi(G) \pi(\bar A)$ and similarly for $B$, so $\pi(G)$ divides
both images. Since $\operatorname{lc}_{x_v}(G)$ divides
$\operatorname{lc}_{x_v}(A)$, whose image is nonzero, $\deg \pi(G) =
\deg_{x_v} G$. The images $\pi(\bar A)$ and $\pi(\bar B)$ are coprime
unless their resultant vanishes, and the resultant of the images is the
image of the resultant when the leading coefficients survive; the general
statement with subresultants is in [@GG2003, §6.7].
:::

`mpoly_gcd_info_set_estimates` uses this with a random 63-bit prime and
random evaluation points, for each variable in turn: the degrees of the
univariate GCDs are upper bounds for the degrees of $G$ (exactly the
degrees with high probability), and the numbers of their nonzero terms
estimate the term counts of $G$. If all bounds are zero, $G$ is the GCD
of the contents; if they equal the degrees of $A$ (or $B$), the code
tests whether $A \mid B$ (or $B \mid A$) by exact division. The estimates
then feed the cost model of [](#choice).

## Brown's algorithm {#brown}

Brown's algorithm [@Bro1971] computes images of $\tilde G$, $\bar A$ and
$\bar B$ modulo primes $p$ and reconstructs them by the Chinese remainder
theorem. Modulo $p$ it recurses on the variables: the last variable is
evaluated at points of $\F_p$, the GCDs of the images in one variable
fewer are computed recursively, normalised to the image of the leading
coefficient $\gamma$, and interpolated
(`nmod_mpolyn_gcd_brown_smprime`, described in
[](algorithms/multivariate/modp-gcd)). Over $\Z$:

::: algorithm Brown's algorithm over the integers {#alg-brown}
Input: nonzero $A, B \in \Z[x_0, \ldots, x_{n-1}]$ with their integer
contents removed (their GCD is multiplied back at the end); $\gamma =
\gcd(\operatorname{lc}(A), \operatorname{lc}(B)) \in \Z$ for the leading
coefficients in lex order, which are integers.
Output: $G$, $\bar A$, $\bar B$.

1. Set $\beta = 2 \gamma \max(\|A\|_\infty, \|B\|_\infty)$ and $m = 1$.
2. Take the next prime $p > 2^{63}$ not dividing $\gamma$. Compute
   $g_p = \gcd(A \bmod p, B \bmod p)$ and the cofactors over $\F_p$, scaled
   so that $\operatorname{lc}(g_p) = \gamma \bmod p$. If $g_p = 1$, return
   $G = 1$.
3. Compare the leading monomial of $g_p$ with that of the current $H$: if
   larger, discard $p$ (unlucky) and go to 2; if smaller, discard all
   previous primes and set $m = 1$.
4. Combine $H \equiv g_p$, $\bar A_H \equiv \bar a_p$, $\bar B_H \equiv
   \bar b_p$ modulo $mp$ with symmetric remainders, and set $m \gets mp$.
5. If $m \le \beta$, go to 2.
6. If $2 \min(\|H\|_\infty \|\bar A_H\|_1, \|H\|_1 \|\bar A_H\|_\infty) < m$
   and the same holds for $\bar B_H$, return $G = \operatorname{pp}(H)$,
   $\bar A = \bar A_H \cdot \operatorname{cont}(H)/\gamma$, and similarly
   $\bar B$. Otherwise go to 2.
:::

::: theorem Correctness of Brown's algorithm {#thm-brown}
If [](#alg-brown) returns, its output is correct, and it returns after
finitely many primes.
:::

::: proof
Each image satisfies $g_p \bar a_p = \gamma A \bmod p$, and so does every
combination: $H \bar A_H \equiv \gamma A \pmod m$. The coefficients of the
product $H \bar A_H$ are bounded in absolute value by $\min(\|H\|_\infty
\|\bar A_H\|_1, \|H\|_1 \|\bar A_H\|_\infty)$, since each coefficient is a
sum of products $h_i a_j$ in which every $h_i$, and every $a_j$, occurs at
most once. By the test of step 6 these are less than $m/2$, and those of
$\gamma A$ are less than $m/2$ by step 5. Two polynomials with
coefficients in $(-m/2, m/2)$ that are congruent modulo $m$ are equal, so
$H \bar A_H = \gamma A$, and similarly $H \bar B_H = \gamma B$. Hence $H$
divides $\gamma A$ and $\gamma B$, its primitive part divides $A$ and $B$
(Gauss's lemma), and so divides $G$. Conversely, all images used have
leading monomial at most that of $G$ by [](#lem-image-deg), so the leading
monomial of $H$ is that of $G$, and a divisor of $G$ with the same leading
monomial and positive leading coefficient is $G$ up to a unit, since $G$
is primitive. For termination: only finitely many primes divide $\gamma$
or the nonzero integer that detects unluckiness ([](#lem-image-deg)
applied to reduction modulo $p$), and once the combination uses only good
primes it equals $\tilde G$, $\bar A$, $\bar B$ times $\gamma$-scalings,
whose norms are fixed, so step 6 eventually succeeds.
:::

The test of step 6 replaces trial division, which is the expensive
final step of most modular GCD algorithms. Brown's algorithm costs about
as much as a few dense multiplications of the inputs: it is the method of
choice for dense inputs and for two variables, and hopeless for sparse
inputs of high degree in many variables, where the dense box has
$\prod_i (d_i + 1)$ entries.

## Zippel's algorithm {#zippel}

Zippel's algorithm [@Zip1979] exploits sparsity through the
*shape* assumption: the image of $G$ modulo a random prime has the same
monomials as $G$. The first prime gives the monomials; for the others,
only the coefficients are unknown, and they are found from few
evaluations by solving linear systems (*sparse interpolation*,
[](algorithms/multivariate/sparse-interpolation#zippel)).

::: algorithm Zippel's GCD over the integers {#alg-zippel}
Input: $A$, $B$ in at least two variables, primitive with respect to
$x_0$ over $\Z[x_1, \ldots]$; $\gamma = \gcd(\operatorname{lc}(A),
\operatorname{lc}(B)) \in \Z$ for the lex leading coefficients.
Output: $G$, $\bar A$, $\bar B$, or failure.

1. Choose a prime $p$ with $\gamma \not\equiv 0$ and $A, B \not\equiv 0$.
   Compute $g_p = \gcd(A \bmod p, B \bmod p)$ by recursive sparse
   interpolation over $\F_p$ (`nmod_mpolyl_gcdp_zippel_smprime`); if $1$,
   return $G = 1$; if its degree in $x_0$ exceeds the current bound, try
   another prime; normalise to $\operatorname{lc}(g_p) = \gamma$ and lift
   it to $H$ with symmetric coefficients.
2. For further primes $p'$: compute the image $g_{p'}$ with the monomials
   of $H$ by solving for its coefficients
   (`nmod_mpolyl_gcds_zippel`); on inconsistency go back to 1 with a new
   first prime. Normalise and combine with $H$ by the Chinese remainder
   theorem.
3. When $H$ no longer changes, divide $H$ by its content and test whether
   it divides $A$ and $B$ exactly. If so, return. If $H$ changes but its
   coefficients exceed a bound, go back to 1.
:::

The result is correct because of the final divisions and the degree
bound: a primitive common divisor of $A$ and $B$ whose degree in $x_0$ is
the degree of the GCD of images, which by [](#lem-image-deg) is at least
$\deg_{x_0} G$, is $G$ itself. The leading coefficient of $G$ in $x_0$ is
in general a polynomial in $x_1, \ldots$, not a monomial, so the
univariate images computed during sparse interpolation are only known up
to scalar multiples; FLINT treats these scalars as additional unknowns
of the linear systems, as proposed by de Kleine, Monagan and Wittkopf
[@dKMW2005].

The bound of step 3, in bits, is the bit size of the smaller of the
bounds of `fmpz_mpoly_factor_bound_si` for the divisors of $A$ and of
$B$, which use the degrees in all variables, plus the bit size of
$\gamma$.

## Sparse interpolation with Berlekamp–Massey {#zippel2}

Zippel's algorithm needs one recursive sparse interpolation per variable.
Hu and Monagan [@HuMon2016] instead treat two *main variables* $x_0,
x_1$ densely and all others at once: the coefficients of $G$ as a
polynomial in $x_0, x_1$ are sparse polynomials in $x_2, \ldots$, which
are recovered by the algorithm of Ben-Or and Tiwari [@BenOrTiw1988] from
their values at the powers of one point
([](algorithms/multivariate/sparse-interpolation#bma)).

::: algorithm GCD with Berlekamp–Massey interpolation {#alg-zippel2}
Input: $A$, $B$ in at least three variables, primitive with respect to
$x_0, x_1$ over $\Z[x_2, \ldots]$; $\Gamma = \gcd(\operatorname{lc}_{x_0,
x_1}(A), \operatorname{lc}_{x_0, x_1}(B)) \in \Z[x_2, \ldots]$.
Output: $G$, $\bar A$, $\bar B$, or failure.

1. Bound the bidegree of $G$ in $(x_0, x_1)$ from bivariate images at
   random points modulo a word-size prime.
2. Choose Kronecker substitution degrees for $x_2, \ldots, x_{n-1}$ (degree
   bounds plus one), so that the coefficients become univariate in one
   variable $z$ with exponents below $D = \prod_{i \ge 2} (d_i + 1)$, and a
   prime $p > D$ for which $p - 1$ is smooth (`fmpz_next_smooth_prime`),
   so that discrete logarithms modulo $p$ are cheap (Pohlig–Hellman
   [@PohHel1978]).
3. For $k = 1, 2, \ldots$: evaluate $z = \omega^{k}$ for a generator
   $\omega$, compute the bivariate GCD of the images and its cofactors
   (`n_polyu1n_mod_gcd_brown_smprime`), scaled to have leading
   coefficient the image of $\Gamma$, and feed each coefficient to a
   Berlekamp–Massey instance. Stop when the instances for $\tilde G$, or
   for one of the cofactors, have stabilised and their roots give
   monomials $\omega^{e}$ with $e < D$ (discrete logarithms); check the
   candidate at a random point.
4. With the monomials known, compute images modulo further primes by
   solving transposed Vandermonde systems (`zip_solve`), and combine them
   by the Chinese remainder theorem until stable.
5. Make the result primitive and verify by exact division.
:::

The number of evaluations is about twice the number of terms of the
largest coefficient of $\tilde G$ (or of the cofactor used), against
roughly the number of terms times the number of variables for Zippel's
algorithm. Using $\Gamma$ to normalise the images solves the leading
coefficient problem at the price of possibly introducing a factor
$\Gamma / \operatorname{lc}(G)$, removed at the end.

## Hensel lifting {#hensel}

The EEZ-GCD algorithm [@MosYun1973; @Wan1980] views the GCD problem as a
factorisation problem: if $g = \gcd(A(x_0, \alpha), B(x_0, \alpha))$ for an
evaluation point $\alpha \in \Z^{n-1}$ and $g$ is coprime to the
cofactor $\bar a = A(x_0, \alpha)/g$, then $A = G \bar A$ is the unique
lifting of $A(x_0, \alpha) = g \cdot \bar a$ with the correct leading
coefficients ([](algorithms/multivariate/hensel-lifting#lem-unique)).

::: algorithm EEZ-GCD {#alg-eez}
Input: $A$, $B$ primitive with respect to $x_0$; an upper bound $d$ for
$\deg_{x_0} G$.
Output: $G$, $\bar A$, $\bar B$, or failure.

1. Choose $\alpha$ (first $0$, then random values of increasing size) with
   $\deg_{x_0} A(x_0, \alpha) = \deg_{x_0} A$ and likewise for $B$. Compute
   $g$ and the cofactors $\bar a$, $\bar b$ of the univariate images. If
   $\deg g = 0$, return $G = 1$; if $\deg g > d$, try another point; if
   $\deg g < d$, lower $d$ and try another point (the bound now has two
   witnesses only when two points agree).
2. If $\deg g$ equals $\deg_{x_0} A$, test whether $A \mid B$, and
   symmetrically.
3. Choose $H = A$ if $\gcd(g, \bar a) = 1$, else $H = B$ if $\gcd(g, \bar
   b) = 1$, else $H = \mu_1 A + \mu_2 B$ for small rationals $\mu_1/\mu_2$
   (enumerated in Calkin–Wilf order) such that $H$ keeps its degree and
   its image cofactor is coprime to $g$.
4. Let $\gamma = \gcd(\operatorname{lc}_{x_0} A, \operatorname{lc}_{x_0}
   B)$. Lift the factorisation $\gamma H(x_0, \alpha) = g' \cdot h'$, where
   $g'$, $h'$ are $g$ and $H(x_0,\alpha)/g$ scaled to have leading
   coefficients $\gamma(\alpha)$ and $\operatorname{lc}(H)(\alpha)$, to a
   factorisation $\gamma H = G' H'$ with leading coefficients $\gamma$ and
   $\operatorname{lc}_{x_0} H$, one variable at a time with
   `fmpz_mpoly_hlift`. On failure, try another point.
5. Let $G = \operatorname{pp}_{x_0}(G')$ and verify that $G$ divides $A$
   and $B$.
:::

::: theorem Correctness of EEZ-GCD {#thm-eez}
If [](#alg-eez) returns $G$, then $G = \gcd(A, B)$.
:::

::: proof
$G$ divides $A$ and $B$ by step 5, hence divides $\gcd(A, B)$. Its degree in
$x_0$ is $\deg g$, which by [](#lem-image-deg) is at least $\deg_{x_0}
\gcd(A, B)$, because the evaluation preserves the degrees of $A$ and $B$.
So $\gcd(A, B)/G$ has degree 0 in $x_0$, and as a divisor of the primitive
polynomial $\gcd(A, B)$ it is a unit.
:::

When the GCD has few terms but a dense structure in few variables, Hensel
lifting is efficient, and it needs only one good evaluation point, which
matters over small finite fields ([](algorithms/multivariate/modp-gcd)).
Its weakness is the cost of dense lifting for sparse polynomials in many
variables. `fmpz_mpoly_gcd` does not use it; it is available as
`fmpz_mpoly_gcd_hensel`. A last method, `fmpz_mpoly_gcd_subresultant`,
computes the last nonzero subresultant in a well chosen main variable
with Ducos' algorithm [@Duc2000] and is useful only for very small
inputs.

## Choosing an algorithm {#choice}

After the easy cases, the inputs are converted to a lex context in the
variables that occur, ordered differently for each algorithm (the
variables of largest degree become the dense variables of Brown's
algorithm and the main variables of the others), and the estimates of
[](#bounds) give costs:

- for Brown's algorithm, roughly the size of the dense boxes of $A$, $B$,
  $G$ and the cofactors, adjusted by the densities $\ell_A/\prod_i
  (\deg_{x_i} A + 1)$ and the expected number of points before
  stabilisation (`mpoly_gcd_info_measure_brown`);
- for the Berlekamp–Massey algorithm, the number of evaluations (the
  maximal number of terms of a coefficient of $G$ in the two main
  variables) times the cost of a bivariate GCD and of evaluating $A$ and
  $B$ (`mpoly_gcd_info_measure_bma`).

With at most two variables Brown's algorithm is used. With three or more,
the Berlekamp–Massey algorithm is tried first if its estimated time is
less than Brown's and either the inputs are sparse
($n(\delta_A + \delta_B) < 1$ for the densities $\delta$) or its estimate
is below 1% of Brown's; otherwise Brown's algorithm comes first. If both
fail, Zippel's algorithm is tried. Each algorithm can fail (for instance
when it runs out of good points, or when the exponents do not fit), never
returning a wrong result, and the next one is used.

## Use in FLINT {#use}

- Driver: `fmpz_mpoly_gcd`, `fmpz_mpoly_gcd_cofactors` and the forced
  variants through `_fmpz_mpoly_gcd_algo`, with the bookkeeping of
  `mpoly_gcd_info_t` (`mpoly_gcd_info_limits`, `mpoly_gcd_info_stride`,
  `mpoly_gcd_info_set_perm`, `mpoly_gcd_info_measure_brown`,
  `mpoly_gcd_info_measure_bma`, `mpoly_gcd_info_measure_zippel`,
  `mpoly_gcd_info_measure_hensel`).
- Brown: `fmpz_mpolyl_gcd_brown` ([](#alg-brown), [](#thm-brown)), with
  `nmod_mpolyn_gcd_brown_smprime` for the images and the interpolation
  helpers `fmpz_mpoly_interp_reduce_p_mpolyn`,
  `fmpz_mpoly_interp_crt_p_mpolyn`; a threaded version over several primes
  is `fmpz_mpolyl_gcd_brown_threaded_pool`.
- Zippel: `fmpz_mpolyl_gcd_zippel` ([](#alg-zippel)), with
  `nmod_mpolyl_gcdp_zippel_smprime` and `nmod_mpolyl_gcds_zippel`.
- Berlekamp–Massey: `fmpz_mpolyl_gcd_zippel2` ([](#alg-zippel2)).
- Hensel: `fmpz_mpolyl_gcd_hensel` ([](#alg-eez)), with
  `fmpz_mpoly_hlift`.
- Subresultants: `fmpz_mpoly_gcd_subresultant`.
- `fmpq_mpoly_gcd` and friends apply these to the integer parts; the
  finite field versions are described in [](algorithms/multivariate/modp-gcd).

::: functions fmpz_mpoly_gcd fmpz_mpoly_gcd_cofactors fmpz_mpoly_gcd_brown fmpz_mpoly_gcd_zippel fmpz_mpoly_gcd_zippel2 fmpz_mpoly_gcd_hensel fmpz_mpoly_gcd_subresultant fmpz_mpolyl_gcd_brown fmpz_mpolyl_gcd_zippel fmpz_mpolyl_gcd_zippel2 fmpz_mpolyl_gcd_hensel :::
