---
title: Factorisation
summary: How FLINT factors multivariate polynomials over the integers: removal of contents and powers of variables, squarefree decomposition by Yun's algorithm, compression of the exponents by unimodular changes of variables, irreducibility certificates from Newton polytopes, direct formulas for low degree, dense bivariate factorisation with Hensel lifting and recombination, degree pruning with univariate images, and the three lifting strategies of Wang, Zippel and Zassenhaus.
complexity: Polynomial in the dense size of the input for the dense methods and in the number of terms for sparse Hensel lifting, in practice; the recombination of local factors is exponential in their number in the worst case.
order: 70
---

The factorisation of a polynomial $A \in \Z[x_0, \ldots, x_{n-1}]$ into
irreducibles, `fmpz_mpoly_factor`, proceeds in stages that each reduce
the problem: contents and monomial factors are removed, then repeated
factors, then the exponents are compressed so that as few variables and
degrees as possible remain, and finally each squarefree primitive part is
split into irreducible factors by a method chosen from its size and
density. Most of these stages are cheap and settle many inputs alone;
the expensive last stage is the classical reduction to a univariate (or
bivariate) problem by evaluation, followed by Hensel lifting
[@Mus1975; @WanRot1975; @Wan1978; @Kal1985]. The function-level
description and examples are in [](fmpz_mpoly_factor/factor), the lifting
algorithms in [](algorithms/multivariate/hensel-lifting), and the
specific problems of finite fields in
[](algorithms/multivariate/modp-factoring). The code is by Schultz, whose
notes [@Sch-mpoly] contain further details.

## The pipeline {#pipeline}

1. *Content* (`fmpz_mpoly_factor_content`): the signed integer content
   and the largest monomial dividing $A$ are removed, and the rest is
   split by contents with respect to single variables (below).
2. *Squarefree decomposition* (`_fmpz_mpoly_factor_squarefree`) of each
   piece.
3. *Compression* (`mpoly_compression_set`): a unimodular change of the
   exponents reduces the number of variables and the degrees, and the
   Newton polytope may prove irreducibility.
4. *Irreducible factorisation* of each squarefree piece in the compressed
   variables: special cases by degree and number of variables, otherwise
   evaluation, pruning and lifting.
5. The factors are mapped back by the inverse change of variables, made
   primitive with positive leading coefficient, and collected with their
   multiplicities.

Each stage only needs GCDs, exact divisions and recursive factorisations
of smaller polynomials, and every factor returned is either proved
irreducible or obtained by a computation that is verified by exact
division.

## Contents {#content}

For a variable $x_v$, write $A = \sum_k a_k x_v^k$ with $a_k$ free of $x_v$.
The content $\operatorname{cont}_{x_v}(A) = \gcd_k(a_k)$ divides $A$, and
$A / \operatorname{cont}_{x_v}(A)$ is primitive with respect to $x_v$.
`fmpz_mpoly_factor_content` repeatedly picks a variable, computes the
content with `_fmpz_mpoly_vec_content_mpoly`, and splits $A$ into the
content and the cofactor, until every piece is primitive with respect to
every variable it involves. The pieces are pairwise coprime: a piece
obtained as a content with respect to $x_v$ is free of $x_v$, while the
cofactor is primitive with respect to $x_v$ and has no factor free of
$x_v$.

## Squarefree decomposition {#squarefree}

::: lemma Yun's algorithm {#lem-yun}
Let $K$ be a field of characteristic zero, $A \in K[x_0, \ldots, x_{n-1}]$
primitive with respect to $x_v$ with $\partial A / \partial x_v \ne 0$, and
$A = \prod_{k \ge 1} S_k^k$ with $S_k$ squarefree and pairwise coprime (the
squarefree decomposition with respect to factors involving $x_v$). Let $G
= \gcd(A, A')$ with $A' = \partial A/\partial x_v$, $C_1 = A/G$ and $D_1 =
A'/G - C_1'$. Then $S_1 = \gcd(C_1, D_1)$, and with $C_{k+1} = C_k / S_k$,
$D_{k+1} = D_k / S_k - C_{k+1}'$, also $S_{k+1} = \gcd(C_{k+1},
D_{k+1})$ [@Yun1976].
:::

::: proof
Every irreducible factor $P$ of $A$ involves $x_v$ (otherwise it would
divide $\operatorname{cont}_{x_v} A = 1$), so $\partial P / \partial x_v
\ne 0$ and $P$ divides $A'$ exactly $k - 1$ times if $P \mid S_k$; hence
$G = \prod_k S_k^{k-1}$ and $A'/G = C_1 \sum_j j S_j'/S_j$. One checks by
induction that $C_k = \prod_{j \ge k} S_j$ and $D_k = C_k \sum_{j \ge k} (j
- k) S_j'/S_j$. An irreducible factor of $S_k$ divides every term of this
sum (the term $j = k$ vanishes, the others contain $S_k$), so it divides
$D_k$. For $P \mid S_j$ with $j > k$, $D_k \equiv (C_k/S_j)(j - k) S_j'
\not\equiv 0 \pmod P$ in characteristic zero, because $S_j$ is squarefree.
So $\gcd(C_k, D_k) = S_k$.
:::

Over $\Z$ the same computation with primitive GCDs gives the squarefree
decomposition, since the integer content has been removed and the pieces
from the content stage are primitive with respect to every variable they
involve, which is the hypothesis of the lemma for any such variable. FLINT
uses the first variable
with nonzero derivative (`_fmpz_mpoly_factor_squarefree`). In
characteristic $p$ the derivative can vanish identically and $p$-th roots
are needed ([](algorithms/multivariate/modp-factoring#separable)).

## Compression {#compression}

Many polynomials met in practice depend on fewer "true" variables than
they seem to, for instance $x^2 y^4 + x y^2 z + z^2$ is a polynomial in
$x y^2$ and $z$. FLINT detects this with integer linear algebra on the
exponent vectors.

::: lemma Unimodular changes of variables {#lem-compress}
Let $A \in K[x_1, \ldots, x_n]$ (with $K$ a field or $\Z$) be divisible by no
variable. Let $V \in \mathrm{GL}_n(\Z)$ and $D \in \Z^n$ be such that every
exponent vector of $A$ is $e = V e' + D$ with $e' \in \N^n$, and let $L \in
K[y_1, \ldots, y_n]$ have the terms $c\, y^{e'}$ for the terms $c\, x^e$ of
$A$. Then the factorisations of $A$ and of $L$ (with the largest monomial
factor of $L$ removed) correspond: each irreducible factor $P$ of $L$
not divisible by a variable gives the irreducible factor of $A$ obtained by
mapping the exponents of $P$ with $e' \mapsto V e'$ and dividing by the
largest monomial factor, and conversely.
:::

::: proof
In the Laurent polynomial ring $K[x_1^{\pm 1}, \ldots, x_n^{\pm 1}]$ the
units are $c\, x^e$ with $c$ a unit of $K$, and $x^{e'} \mapsto x^{V e'}$
defines a ring automorphism $\varphi$, because $V$ is invertible over $\Z$.
Thus $\varphi$ maps irreducible Laurent polynomials to irreducible ones,
and $A = x^D \varphi(L)$. A polynomial divisible by no variable is
irreducible in $K[x]$ if and only if it is irreducible in the Laurent ring
(a factorisation in the Laurent ring becomes a polynomial factorisation
after multiplying the factors by suitable monomials, and the monomial
factors cancel since the variables are prime and do not divide $A$). Each
Laurent irreducible is associated to a unique polynomial divisible by no
variable, which gives the correspondence.
:::

`mpoly_compression_set` builds the exponent matrix, subtracts the
minimal exponents, and applies greedy unimodular column operations to
reduce the number of nonzero columns and the product of the degrees
(`_mpoly_compress_exps`); a column that becomes zero is a variable that disappears. If the result
is the identity, the polynomial is factored as it is; otherwise it is
factored in the new variables (`fmpz_mpoly_compression_do`, lex order) and
the factors are mapped back (`fmpz_mpoly_compression_undo`). A change of
variables may make a polynomial non-primitive or non-squarefree in the new
variables, so contents and the squarefree decomposition are recomputed
there.

### Irreducibility from the Newton polytope {#polytope}

The Newton polytope $\operatorname{Newt}(A)$ is the convex hull of the
exponent vectors of $A$. Ostrowski's theorem [@Ost1921] states that
$\operatorname{Newt}(FG) = \operatorname{Newt}(F) + \operatorname{Newt}(G)$
(Minkowski sum) over any field, so if $\operatorname{Newt}(A)$ cannot be
written as a Minkowski sum of two lattice polytopes with more than one
point each (*integrally indecomposable*), then $A$ has no nontrivial
factorisation over any field, unless one factor is a monomial; for $A$
divisible by no variable, $A$ is then absolutely irreducible. Gao and
Lauder [@Gao2001; @GaoLau2001] give a practical test: project the
exponents to the plane with a random integer matrix and decide whether the
projected polygon is integrally indecomposable, which reduces to a subset
sum problem on the edge vectors of the polygon. FLINT runs a few such
projections (`_mpoly_test_irreducible`, called from
`mpoly_compression_set`), with bounds on the memory used; success proves
irreducibility, failure says nothing. Polynomials with a few terms and
sparse "random" polynomials are typically proved irreducible this way
without any factorisation.

FLINT does not compute absolute factorisations (over an algebraic
closure of the coefficient field); the polytope test is only used as a
certificate of irreducibility.

## Small cases {#small}

After compression and squarefree decomposition, a primitive squarefree
$A$ in lex order with main variable $x_0$ is handled directly when:

- $\deg_{x_0} A = 1$: $A$ is irreducible, since a factor of degree 0 in
  $x_0$ would divide the content;
- $\deg_{x_0} A = 2$: $A$ factors by [](#lem-quadratic), with the square
  root computed by `fmpz_mpoly_sqrt` (`_apply_quadratic`);
- one variable: `fmpz_poly_factor`
  ([](algorithms/polynomials/factoring-integers));
- two variables and a dense polynomial ($(d_0 + 1)(d_1 + 1) \le 400$, or at
  least a fifth of the dense box occupied): dense bivariate factorisation
  ([](#bivariate)).

::: lemma Quadratics {#lem-quadratic}
Let $R$ be a unique factorisation domain of characteristic not 2 and $A =
a x^2 + b x + c \in R[x]$ primitive with $a \ne 0$. Then $A$ is reducible if
and only if $b^2 - 4ac$ is a square $d^2$ in $R$, and then $A$ is, up to a
unit, the product of the primitive parts of $2ax + b - d$ and $2ax + b + d$.
:::

::: proof
If $A = (px + q)(rx + s)$, then $b^2 - 4ac = (ps - qr)^2$. Conversely, $4aA =
(2ax + b)^2 - d^2 = (2ax + b - d)(2ax + b + d)$; by Gauss's lemma the
primitive parts of the two factors multiply to the primitive part of
$4aA$, which is $A$ up to a unit.
:::

## Dense bivariate factorisation {#bivariate}

`fmpz_bpoly_factor` factors $B \in \Z[y][x]$, squarefree and primitive in
$x$, in the classical way [@Mus1975; @GG2003, §15.6]:

1. Find $\alpha \in \{0, 1, -1, 2, -2, \ldots\}$ such that $B(x, \alpha)$
   has the same degree in $x$ and is squarefree; factor it with
   `fmpz_poly_factor`. If it is irreducible, so is $B$.
2. Choose a prime $p$ not dividing the leading coefficient and a precision
   $k$ such that $p^k$ exceeds a bound on the coefficients of the factors
   of $B$; make the local factors monic modulo $p^k$ and lift them in
   $(\Z/p^k\Z)[[y - \alpha]][x]$ to a precision in $y - \alpha$ determined by
   the degrees of $B$ and of its leading coefficient
   ([](algorithms/multivariate/hensel-lifting#bivariate)).
3. *Recombination*: for subsets $S$ of the lifted factors, of increasing
   size, form $\operatorname{lc}(B) \prod_{i \in S} f_i$, truncate, reduce
   symmetrically modulo $p^k$, take the primitive part and test whether it
   divides $B$. Each successful test gives an irreducible factor and
   removes $S$.

::: lemma Recombination {#lem-recomb}
Every irreducible factor $F$ of $B$ is, up to an integer scalar, the
primitive part of $\operatorname{lc}(F) \prod_{i \in S} f_i$ truncated at
the lifting precision and reduced modulo $p^k$, for a unique subset $S$,
provided the precision and $p^k$ exceed the bounds on the degrees and
coefficients of the factors of $B$.
:::

::: proof
$F(x, \alpha)$ divides $B(x, \alpha)$, which is squarefree, so it is a
product of some of the local factors; by the uniqueness of Hensel lifting
for coprime monic factors, $F/\operatorname{lc}(F)$ is the product of the
corresponding lifted factors in $(\Z/p^k\Z)[[y-\alpha]][x]$. Multiplying by
$\operatorname{lc}(F)$, which divides $\operatorname{lc}(B)$, gives a
polynomial with integer coefficients smaller than $p^k/2$ and degree in
$y$ smaller than the precision, which is $F$; multiplying by
$\operatorname{lc}(B)$ instead gives a multiple of $F$ by a polynomial in
$y$, whose primitive part is $F$.
:::

The recombination costs $O(2^r)$ trial divisions for $r$ local factors in
the worst case; the polynomial-time recombination of van Hoeij
[@vHo2002] and Lecerf's bounds [@Lec2006] are not used. In practice $r$
is small for bivariate images, because a random univariate image of an
irreducible bivariate polynomial rarely has many factors.

## Evaluation and degree pruning {#pruning}

For more variables, or sparse bivariate inputs, the variables $x_1,
\ldots, x_{n-1}$ are evaluated at integer points, first at $0$ and then
at random points with coordinates in ranges $[-m/2, m/2)$ for $m = 6,
\ldots, 10$, until three images that preserve the degree in $x_0$ and are
squarefree have been factored (or the points run out).

::: lemma Degree pruning {#lem-prune}
Let $A$ be primitive and squarefree with $\deg_{x_0} A = d$, and let
$A(x_0, \alpha) = c \prod_i u_i^{m_i}$ be a univariate image of the same
degree. If $A$ has a factor of degree $k$ in $x_0$, then $k$ is a sum of
the degrees of a sub-multiset of the $u_i$ (with multiplicities).
:::

::: proof
If $A = FG$ then $F(x_0, \alpha)$ divides $A(x_0, \alpha)$ and has the same
degree as $F$ because the leading coefficient of $A$, and hence of $F$,
does not vanish at $\alpha$. By unique factorisation in $\Z[x_0]$, it is a
product of some of the $u_i$.
:::

`zassenhaus_prune_t` keeps, for each $0 < k < d$, whether $k$ is a subset
sum for every image seen so far (`zassenhaus_prune_end_add_factors`
intersects the sets). If no $k$ survives, $A$ is irreducible
(`zassenhaus_prune_must_be_irreducible`), and no lifting is needed; this
settles most irreducible inputs after three images. The same structure is
used inside the lifting methods to stop early.

## Lifting methods {#methods}

What remains is a squarefree primitive $A$ that may be reducible. FLINT
factors the leading coefficient $\operatorname{lc}_{x_0}(A)$ (recursively,
into irreducibles when its factors are short) and then tries, in this
order, the methods allowed by the `algo` flags:

- *Wang* (`fmpz_mpoly_factor_irred_wang`), first at the point $\alpha = 0$
  if $A$ is dense enough (more than $0.2\%$ of its dense box) and $0$ was
  a good evaluation point, then at random points if $A$ is denser than
  $4\%$: factor $A(x_0, \alpha)$, distribute the leading coefficient with
  Wang's or Kaltofen's method
  ([](algorithms/multivariate/hensel-lifting#lcc)), multiply $A$ by a power
  of the unassigned part, and lift one variable at a time with dense Hensel
  lifting ([](algorithms/multivariate/hensel-lifting#hlift)).
- *Zippel* (`fmpz_mpoly_factor_irred_zippel`): the same, but the lifting
  is done modulo a word-size prime, densely for $x_1$ and with sparse
  Hensel lifting for the other variables, followed by $p$-adic lifting
  of the coefficients
  ([](algorithms/multivariate/hensel-lifting#sparse)); this is the method
  for sparse inputs in many variables [@Zip1981; @MonTun2020].
- *Wang* again at random points, if Zippel's method failed.
- *Zassenhaus* (`fmpz_mpoly_factor_irred_zassenhaus`) as the last resort,
  described below.

Every method returns failure rather than a wrong answer: the lifted
factors are verified by exact division, and a method that cannot proceed
(bad evaluation points after several tries, exponents too large, a
leading coefficient that cannot be distributed) gives way to the next.

### Zassenhaus-style lifting {#zassenhaus}

Without leading coefficient information, the factors can still be lifted
if the leading coefficient is imposed on *every* factor: if $A = \prod_{i
\le r} A_i$ and $\ell = \operatorname{lc}_{x_0}(A)$, then $\ell^{r-1} A =
\prod_i (\ell / \operatorname{lc}(A_i)) A_i$ is a factorisation with all
leading coefficients equal to $\ell$, which is unique by
[](algorithms/multivariate/hensel-lifting#lem-unique). The method:

1. Evaluate $x_2, \ldots, x_{n-1}$ at a point and factor the bivariate
   image in $x_0, x_1$ with `fmpz_bpoly_factor`.
2. For $m = 2, \ldots, n - 1$, try to lift all current factors to the
   variable $x_m$ with the leading coefficient imposed on each. If this
   fails, the image has more factors than $A^{(m)}$: try to lift
   splittings into two factors, a product of a subset $S$ of the current
   factors and the product of the rest, for subsets of increasing size, as
   in [](#lem-recomb); each success splits off a factor.
3. The factors found at the last variable are made primitive.

With two local factors that cannot be lifted, or one, $A$ is irreducible.
This method is exponential in the number of bivariate factors in the worst
case, but it needs no factorisation of the leading coefficient and no
condition on the evaluation point.

## Use in FLINT {#use}

- Driver: `fmpz_mpoly_factor`, `fmpz_mpoly_factor_wang`,
  `fmpz_mpoly_factor_zippel`, `fmpz_mpoly_factor_zassenhaus`,
  `fmpz_mpoly_factor_irred`, and the static functions `_factor`,
  `_factor_irred`, `_factor_irred_compressed` in
  `src/fmpz_mpoly_factor/factor.c`.
- Contents and squarefree parts: `fmpz_mpoly_factor_content`,
  `fmpz_mpoly_factor_squarefree`, `_fmpz_mpoly_factor_squarefree`
  ([](#lem-yun)).
- Compression and polytopes: `mpoly_compression_set`,
  `fmpz_mpoly_compression_do`, `fmpz_mpoly_compression_undo`
  ([](#lem-compress)), `mpoly_test_irreducible`.
- Bivariate: `fmpz_bpoly_factor`, `fmpz_bpoly_factor_ordered`
  ([](#lem-recomb)).
- Pruning: `zassenhaus_prune_t` and its functions ([](#lem-prune)).
- Lifting: `fmpz_mpoly_factor_irred_wang`, `fmpz_mpoly_factor_irred_zippel`,
  `fmpz_mpoly_factor_irred_zassenhaus`, `fmpz_mpoly_hlift`,
  `fmpz_mpoly_factor_lcc_wang`, `fmpz_mpoly_factor_lcc_kaltofen`.
- `fmpq_mpoly_factor` reduces to `fmpz_mpoly_factor`; the factorisation
  over finite fields (`nmod_mpoly_factor`, `fq_nmod_mpoly_factor`,
  `fmpz_mod_mpoly_factor`) follows the same pipeline with the changes
  described in [](algorithms/multivariate/modp-factoring).

::: functions fmpz_mpoly_factor fmpz_mpoly_factor_squarefree fmpz_mpoly_factor_content fmpz_mpoly_factor_wang fmpz_mpoly_factor_zippel fmpz_mpoly_factor_zassenhaus fmpz_bpoly_factor mpoly_compression_set :::
