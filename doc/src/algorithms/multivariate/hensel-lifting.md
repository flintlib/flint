---
title: Multivariate Hensel lifting
summary: Lifting a factorisation of $A(x_0, \alpha_1, \ldots, \alpha_m)$ to a factorisation of $A$ one variable at a time: uniqueness given the leading coefficients, the multivariate partial fraction (diophantine) equations solved by Taylor expansion, the quartic and quintic error updates of FLINT's `fmpz_mpoly_hlift`, bivariate lifting modulo prime powers, sparse lifting, and the leading coefficient methods of Wang and Kaltofen.
complexity: Dense lifting of $r$ factors of $A$ with degrees $d_0, \ldots, d_m$ costs $O(r \prod_i (d_i + 1)^2)$ coefficient operations up to logarithmic factors in the worst case; sparse lifting is polynomial in the number of terms of the factors.
order: 60
---

Hensel lifting is the step that turns a factorisation of a univariate (or
bivariate) image of $A \in \Z[x_0, \ldots, x_m]$ into a factorisation of
$A$ itself. FLINT uses it in the factoring algorithms
([](algorithms/multivariate/factorisation)), in the EEZ-GCD algorithm
([](algorithms/multivariate/gcd#hensel)) and, over finite fields, in the
corresponding `nmod_mpoly` and `fq_nmod_mpoly` code. The setting is that of
Wang [@WanRot1975; @Wan1978] with the notation of Schultz's notes
[@Sch-mpoly]; see also [@GCL1992, Ch. 6] and [@GG2003, §15].

Throughout, $x_0$ is the *main variable*, $\alpha = (\alpha_1, \ldots,
\alpha_m)$ is the evaluation point for the other variables, and for a
polynomial $F$ we write $F^{(j)} = F(x_0, \ldots, x_j, \alpha_{j+1}, \ldots,
\alpha_m)$. Lifting proceeds from $F^{(0)}$ to $F^{(1)}$, then $F^{(2)}$,
and so on: at stage $j$ the factorisation $A^{(j-1)} = \prod_i F_i^{(j-1)}$
is lifted to $A^{(j)} = \prod_i F_i^{(j)}$ by computing the Taylor
coefficients of the $F_i^{(j)}$ in $x_j - \alpha_j$ one at a time.

## Uniqueness {#uniqueness}

A factorisation can only be lifted if it is determined by its image, and
the leading coefficients in $x_0$ are the obstruction: if $A = FG$, then
also $A = (uF)(u^{-1}G)$ for any unit $u$ of the coefficient ring, and
over $\Z[x_1, \ldots, x_m]$ a factor of the leading coefficient can move
between $F$ and $G$ without changing the images up to scalars. The
lifting therefore assumes that the leading coefficients are known
([](#lcc) explains how they are found).

::: lemma Uniqueness of the lifted factors {#lem-unique}
Let $D$ be an integral domain, $A, F_1, \ldots, F_r \in D[x_0, \ldots,
x_m]$, $\alpha \in D$, and suppose that $A = \prod_i F_i$, that
$\operatorname{lc}_{x_0}(F_i)(x_m = \alpha) \ne 0$, and that the images
$f_i = F_i(x_m = \alpha)$ are pairwise coprime in $K[x_0]$, where $K$ is the
fraction field of $D[x_1, \ldots, x_{m-1}]$. If $G_1, \ldots, G_r$ satisfy
$A = \prod_i G_i$, $G_i(x_m = \alpha) = f_i$ and $\operatorname{lc}_{x_0}(G_i)
= \operatorname{lc}_{x_0}(F_i)$, then $G_i = F_i$ for all $i$.
:::

::: proof
Suppose not, and let $k \ge 1$ be the largest integer such that $G_i
\equiv F_i \pmod{(x_m - \alpha)^k}$ for all $i$. Write $G_i = F_i + (x_m -
\alpha)^k H_i$, with $h_i = H_i(x_m = \alpha)$ not all zero. Since the leading
coefficients agree, $\deg_{x_0} h_i < \deg_{x_0} f_i$. Expanding $\prod_i G_i
- \prod_i F_i = 0$ modulo $(x_m - \alpha)^{k+1}$ gives $\sum_i h_i \prod_{j
\ne i} f_j = 0$. Then $f_i$ divides $h_i \prod_{j \ne i} f_j$, hence divides
$h_i$ by coprimality, and $h_i = 0$ by the degree bound, for every $i$; a
contradiction.
:::

## Multivariate partial fractions {#pfrac}

Each lifting step needs the solution of a *multivariate diophantine
equation*: given pairwise coprime $\beta_1, \ldots, \beta_r$ and $t$, find
$\delta_i$ with
$$
t = \sum_{i=1}^r \delta_i \prod_{j \ne i} \beta_j, \qquad
\deg_{x_0} \delta_i < \deg_{x_0} \beta_i. \tag{$\ast$}
$$
Dividing by $\prod_j \beta_j$, this is the partial fraction decomposition
$t / \prod_j \beta_j = \sum_i \delta_i / \beta_i$.

::: lemma Univariate partial fractions {#lem-pfrac-uni}
Over a field $K$, if $\beta_1, \ldots, \beta_r \in K[x_0]$ are pairwise
coprime and $\deg t < \sum_i \deg \beta_i$, then $(\ast)$ has a unique
solution, given by $\delta_i = t \cdot (\prod_{j \ne i} \beta_j)^{-1} \bmod
\beta_i$.
:::

::: proof
Reducing $(\ast)$ modulo $\beta_i$ gives $\delta_i \prod_{j \ne i} \beta_j
\equiv t$, and $\prod_{j \ne i} \beta_j$ is invertible modulo $\beta_i$;
together with $\deg \delta_i < \deg \beta_i$ this determines $\delta_i$. With
these $\delta_i$, the difference of the two sides of $(\ast)$ is divisible
by every $\beta_i$, hence by their product, and has smaller degree, so it
is zero.
:::

Over $\Z$ the solution has rational coefficients in general; the lifting
needs to know whether an integral solution exists. `fmpz_poly_pfrac_precomp`
solves $(\ast)$ modulo a power $p^k$ of a prime for which the $\beta_i$ stay
coprime, checks the candidate by exact division, and increases $k$ until
a bound on the solutions obtained from Cramer's rule and Hadamard's
inequality (in terms of the norms of the $\beta_i$ and their resultants)
shows that no integral solution exists if the check still fails.

With more variables, the solution is computed variable by variable, as a
dense interpolation in Taylor coefficients.

::: algorithm Multivariate partial fractions {#alg-pfrac}
Input: $l \ge 0$; $t \in D[x_0, \ldots, x_l]$; $\beta_i^{(l)} \in D[x_0,
\ldots, x_l]$ whose images $\beta_i^{(0)}$ are pairwise coprime; degree
bounds $d_1, \ldots, d_l$.
Output: $\delta_i$ satisfying $(\ast)$ for the $\beta_i^{(l)}$, or failure.

1. If $l = 0$, solve with [](#lem-pfrac-uni) (over $\Z$: as described
   above) and return.
2. Otherwise, for $k = 0, 1, \ldots, d_l$: let $\tilde t$ be the $k$-th
   Taylor coefficient in $x_l - \alpha_l$ of $t - \sum_i \delta_i \prod_{j \ne
   i} \beta_j^{(l)}$, solve $(\ast)$ for $\tilde t$ and the
   $\beta_i^{(l-1)}$ recursively, giving $\tilde \delta_i$, and add
   $\tilde\delta_i (x_l - \alpha_l)^k$ to $\delta_i$. Fail if a recursive call
   fails or a $\delta_i$ exceeds the degree bound in $x_l$.
:::

::: theorem {#thm-pfrac}
If $(\ast)$ has a solution with $\deg_{x_l} \delta_i \le d_l$ (and
similarly for the other variables), [](#alg-pfrac) returns it, and the
solution is unique.
:::

::: proof
Uniqueness follows as in [](#lem-unique): two solutions differ by a
solution of the homogeneous equation, whose lowest nonzero Taylor
coefficients in $x_l - \alpha_l$ would solve the homogeneous equation for
the $\beta_i^{(l-1)}$, which by induction (down to [](#lem-pfrac-uni)) only
has the zero solution. The Taylor coefficient of order $k$ of $(\ast)$,
given the coefficients of order less than $k$ of the $\delta_i$, is an
equation of the same form for the coefficients of order $k$ with the
$\beta_i^{(l-1)}$, so the recursion computes the solution coefficient by
coefficient.
:::

`fmpz_mpoly_pfrac` implements [](#alg-pfrac), keeping the products
$\prod_{j \ne i} \beta_j^{(l)}$ and their Taylor expansions for every level
in an `fmpz_mpoly_pfrac_t` prepared by `fmpz_mpoly_pfrac_init`. When
$\alpha_l = 0$ the Taylor coefficients are just the coefficients in
$x_l$, and the expansion is replaced by a conversion to a univariate
polynomial in $x_l$. Wang [@Wan1978] recommends lifting modulo a power
$p^k$ of a prime; FLINT's dense lifting instead works with exact integer
coefficients in all variables, and only the univariate base case
computes modulo prime powers [@Sch-mpoly].

## The lifting step {#hlift}

::: algorithm Hensel lifting in one variable {#alg-hlift}
Input: $A^{(m)}$; factors $F_1, \ldots, F_r$ of $A^{(m-1)}$ whose leading
coefficients in $x_0$ have already been replaced by the correct ones in
$D[x_1, \ldots, x_m]$; the degree $d$ of $A$ in $x_m$.
Output: factors $F_i$ with $A^{(m)} = \prod_i F_i$, or failure.

1. Write $F_i = \sum_k F_{i,k} (x_m - \alpha_m)^k$, with only $F_{i,0}$
   known (apart from the leading coefficients).
2. For $k = 1, \ldots, d$: let $t$ be the $k$-th Taylor coefficient of
   $A^{(m)} - \prod_i F_i$; solve $t = \sum_i \delta_i \prod_{j \ne i} F_{j,0}$
   with [](#alg-pfrac) at level $m - 1$; add $\delta_i$ to $F_{i,k}$. Fail
   if there is no solution or if $\sum_i \deg_{x_m} F_i > d$.
:::

The Taylor coefficient $t$ only depends on the coefficients $F_{i,l}$ with
$l < k$, and FLINT has three ways to compute it:

- the *quintic* version (`_hlift_quintic`, for many factors) keeps the error
  $e = A - \prod_i F_i$ and recomputes it after each step, dividing by
  $(x_m - \alpha_m)^k$ to extract $t$; it finishes by checking $e = 0$;
- the *quartic* version (`_hlift_quartic`, for $3 \le r < 20$) computes
  $t = [(x_m - \alpha_m)^k] A - \sum_{k_1 + \cdots + k_r = k} \prod_i F_{i,k_i}$
  incrementally from the Taylor coefficients, using partial products
  $F_1 F_2$, $F_1 F_2 F_3$, ...;
- `_hlift_quartic2` specialises this to $r = 2$: $t = [(x_m -
  \alpha_m)^k] A - \sum_{l=0}^{k} F_{1,l} F_{2,k-l}$.

The names refer to the exponent of the degree in the dense cost, as in
[@Sch-mpoly]. The quartic versions do not verify the final product; the
callers check the factors by exact division.

::: theorem Correctness of the lifting step {#thm-hlift}
If $A^{(m)}$ has a factorisation $\prod_i G_i$ with $G_i^{(m-1)} = F_i$ (up to
the known leading coefficients) and the conditions of [](#lem-unique)
hold, then [](#alg-hlift) computes it. If it returns failure, no such
factorisation exists.
:::

::: proof
By induction on $k$, after step $k$ the $F_i$ agree with the $G_i$ modulo
$(x_m - \alpha_m)^{k+1}$: the Taylor coefficient of order $k$ of $A -
\prod_i G_i = 0$ is an equation $(\ast)$ for the coefficients $G_{i,k}$,
which by [](#thm-pfrac) has a unique solution, found by the algorithm.
Since $\deg_{x_m} G_i \le d$, the lifting is complete after $d$ steps. If
the factorisation does not exist, either some equation has no solution,
the degrees exceed $d$, or the final product differs from $A$.
:::

### Bivariate lifting modulo prime powers {#bivariate}

For the bivariate factorisation of `fmpz_bpoly_factor` the factors of
$B(x, \alpha) \in \Z[x]$ have no known leading coefficients. They are
made monic modulo $p^k$ and lifted in $(\Z/p^k\Z)[[y - \alpha]][x]$ with the
monic normalisation (`_bivar_lift_quartic2`, `_bivar_lift_quartic` and
`_bivar_lift_quintic` in `fmpz_mpoly_factor/bpoly_factor.c`), to a
precision in $y - \alpha$ given by the degree of $B$ and its leading
coefficient; the true factors are then found by recombination
([](algorithms/multivariate/factorisation#bivariate)). The prime power
$p^k$ is chosen above a bound on the coefficients of the factors.

## Sparse Hensel lifting {#sparse}

Dense lifting costs a polynomial in the dense size $\prod_i (d_i + 1)$,
which is prohibitive for sparse polynomials in many variables. Zippel's
sparse Hensel lifting [@Zip1981] and Kaltofen's variant [@Kal1985]
assume, as in sparse interpolation, that each Taylor coefficient of a
factor in the new variable $x_j$ has the same support in the previous
variables as the factor itself; the coefficients are then found from
bivariate (in $x_0, x_1$) factorisations of $A$ at the points $x_i =
\beta_i^k$ ($2 \le i < j$) by solving transposed Vandermonde systems
([](algorithms/multivariate/sparse-interpolation#hensel)). Monagan and
Tuncer [@MonTun2020] analyse the cost and the probability of failure.

FLINT's `fmpz_mpoly_factor_irred_zippel` reduces modulo a word-size prime
$p$, lifts the first two variables densely and the others with
`nmod_mpoly_hlift_zippel` over $\F_p$, and then lifts the factors from
$\F_p$ to $\Z$ by $p$-adic Newton iteration on the same supports
(the static function `fmpz_mfactor_lift_prime_power_zippel`). The result
is verified by division.

## Leading coefficients {#lcc}

Before lifting, the leading coefficient $\ell = \operatorname{lc}_{x_0}(A)
\in \Z[x_1, \ldots, x_m]$ must be distributed among the factors. If $\ell$ is
an integer, it can be imposed on every factor, at the price of a power
$\ell^{r-1}$ multiplied into $A$ and removed at the end (the integer
content of each lifted factor). Otherwise the factors of $\ell$ must be
assigned to the right factors of $A$. Two methods are implemented, both
needing the factorisation of $\ell$ (FLINT factors it recursively).

### Wang's method {#wang}

::: lemma Wang's leading coefficient condition {#lem-wang}
Let $\ell = c \prod_{j=1}^{s} L_j^{e_j}$ with $c \in \Z$ and distinct
irreducible nonconstant $L_j \in \Z[x_1, \ldots, x_m]$. Let $A = \prod_{i=1}^r
A_i$ be the factorisation of the primitive squarefree $A$ into irreducible
factors, and suppose that $A(x_0, \alpha) = \delta \prod_{i=1}^r u_i$ with
$\delta \in \Z$ and $u_i$ primitive irreducible, $A_i(x_0, \alpha) = \kappa_i
u_i$ with $\kappa_i \in \Z$. If $q$ is a prime that divides $L_j(\alpha)$ but
divides neither $c\delta$ nor any $L_k(\alpha)$, $k \ne j$, then the exponent
$s_{ij}$ of $L_j$ in $\operatorname{lc}_{x_0}(A_i)$ is
$$
s_{ij} = \frac{v_q(\delta \operatorname{lc}(u_i))}{v_q(L_j(\alpha))}.
$$
:::

::: proof
Write $\operatorname{lc}(A_i) = c_i \prod_j L_j^{s_{ij}}$ with $c_i \in \Z$.
Then $\prod_i c_i = c$ and $\prod_i \kappa_i = \delta$, and evaluating at
$\alpha$ gives $\kappa_i \operatorname{lc}(u_i) = c_i \prod_k
L_k(\alpha)^{s_{ik}}$. Multiplying by $\delta/\kappa_i \in \Z$, $\delta
\operatorname{lc}(u_i) = c_i (\delta / \kappa_i) \prod_k L_k(\alpha)^{s_{ik}}$,
where $c_i (\delta/\kappa_i)$ divides $c\delta$ and is therefore not
divisible by $q$. Taking $q$-adic valuations gives $v_q(\delta
\operatorname{lc}(u_i)) = s_{ij} v_q(L_j(\alpha))$.
:::

If every $L_j(\alpha)$ has such a prime, all exponents $s_{ij}$ are
determined [@Wan1978]; the integer parts $c_i$ are not, and are handled by
multiplying $A$ by a power of the remaining integer factor, which is
removed again as the content of the lifted factors. `fmpz_mpoly_factor_lcc_wang` checks a sequential form of the condition
(each $L_j(\alpha)$ must keep a nontrivial factor after removing all
primes of $c\delta$ and of $L_1(\alpha), \ldots, L_{j-1}(\alpha)$; it
returns 0 otherwise, and another point is tried) and computes the
exponents by removing the values $|L_j(\alpha)|$ from $\delta
\operatorname{lc}(u_i)$ in the order $j = s, \ldots, 1$ with `fmpz_remove`. The condition fails for
evaluation points with small values, in particular for $\alpha = 0$ when
$\ell$ has factors without constant term, which is why FLINT first tries
the point $0$ only when the input is dense enough, then random points.

For example, let $A = A_1 A_2$ with $A_1 = x_1(x_1 + x_2) x_0 + 1$ and $A_2
= x_2 x_0 + x_1 - 1$, so that $\ell = x_1 x_2 (x_1 + x_2)$, $c = 1$, and $L_1
= x_1$, $L_2 = x_2$, $L_3 = x_1 + x_2$. At $\alpha = (5, 7)$ the values
$L_j(\alpha) = 5, 7, 12$ are pairwise coprime, and $A(x_0, \alpha) = (60
x_0 + 1)(7 x_0 + 4)$ with $\delta = 1$. Removing $12$, $7$ and $5$ from
$\operatorname{lc}(u_1) = 60$ leaves $1$ after one division by $12$ and
one by $5$, so $\operatorname{lc}(A_1) = L_1 L_3$; from
$\operatorname{lc}(u_2) = 7$ only $7$ is removed, so
$\operatorname{lc}(A_2) = L_2$. At $\alpha = (2, 4)$, on the other hand,
the values $2, 4, 6$ have no private primes and the point is rejected.

### Kaltofen's method {#kaltofen}

Kaltofen [@Kal1985] determines the distribution of the factors of $\ell$
one variable at a time from *bivariate* factorisations: for each $v \ge
1$, the bivariate image $A(x_0, \alpha_1, \ldots, x_v, \ldots, \alpha_m)$ is
factored with factors matched to the univariate factors $u_i$ (by
`fmpz_bpoly_factor_ordered`, which lifts the $u_i$ in $x_v$); the leading
coefficients of these bivariate factors are univariate in $x_v$ and show
which factors of the image of $\ell$ go where; these are lifted back to
the factors of $\ell$ involving $x_v$. This works without Wang's
condition, and `fmpz_mpoly_factor_irred_wang` falls back to it when the
condition fails (`fmpz_mpoly_factor_lcc_kaltofen`,
`fmpz_mpoly_factor_lcc_kaltofen_step`).

Since only coprime factors can be lifted, the step for
$x_v$ works with squarefree parts. The leading coefficients of the
bivariate factors (polynomials in $x_v$, primitive) are replaced by their
squarefree decompositions, which are then refined by GCDs into a single
set of pairwise coprime squarefree bases. Let $P$ be the product of the
distinct irreducible factors of $\ell$ that involve $x_v$, a squarefree
polynomial. If the product of the bases equals
the primitive part of $P(\alpha_1, \ldots, x_v, \ldots, \alpha_m)$, the
bases are lifted by [](#alg-hlift) to a factorisation of $P$ into
polynomials in all variables (the images of the bases are coprime, so
[](#lem-unique) applies). Each lifted base, raised to its multiplicity, is
assigned to the factors whose bivariate leading coefficient contained
it, and divided out of the factorisation of $\ell$. When all of $\ell$
has been distributed, the integer content is assigned by comparing the
contents of the bivariate leading coefficients with those of the
assigned parts.

If neither method distributes the whole leading coefficient, the
remaining part is imposed on every factor as for an integer leading
coefficient, or the Zassenhaus method of
[](algorithms/multivariate/factorisation#zassenhaus), which lifts without
leading coefficient information, is used.

## Use in FLINT {#use}

- Partial fractions: `fmpz_poly_pfrac_precompute`,
  `fmpz_poly_pfrac_precomp` ([](#lem-pfrac-uni)), `fmpz_mpoly_pfrac_init`,
  `fmpz_mpoly_pfrac` ([](#alg-pfrac)).
- Lifting: `fmpz_mpoly_hlift` ([](#alg-hlift)), dispatching to
  `_hlift_quartic2`, `_hlift_quartic` and `_hlift_quintic`;
  `_fmpz_mpoly_get_lead0` and `_fmpz_mpoly_set_lead0` to impose leading
  coefficients; the analogues `nmod_mpoly_hlift` and
  `fq_nmod_mpoly_hlift` over finite fields.
- Bivariate lifting: `fmpz_bpoly_factor`, `fmpz_bpoly_factor_ordered`.
- Sparse lifting: `nmod_mpoly_hlift_zippel`,
  `fmpz_mpoly_factor_irred_zippel`.
- Leading coefficients: `fmpz_mpoly_factor_lcc_wang` ([](#lem-wang)),
  `fmpz_mpoly_factor_lcc_kaltofen`.
- Users: `fmpz_mpoly_factor_irred_wang`, `fmpz_mpoly_factor_irred_zippel`,
  `fmpz_mpoly_factor_irred_zassenhaus`, `fmpz_mpolyl_gcd_hensel`, and the
  finite field factoring and GCD code
  ([](algorithms/multivariate/modp-factoring),
  [](algorithms/multivariate/modp-gcd)).

::: functions fmpz_mpoly_hlift fmpz_mpoly_pfrac fmpz_poly_pfrac_precomp fmpz_mpoly_factor_lcc_wang fmpz_mpoly_factor_lcc_kaltofen fmpz_bpoly_factor nmod_mpoly_hlift_zippel :::
