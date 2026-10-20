---
title: Factoring over finite fields
summary: The multivariate factorisation pipeline of `nmod_mpoly_factor` and `fq_nmod_mpoly_factor`, with what is specific to characteristic $p$: separable factorisation with $p$-th roots, deflation by powers of $p$, factoring over extension fields with recombination of Frobenius orbits when the field is too small for good evaluation points, quadratic equations in characteristic 2, bivariate factorisation with recombination by linear algebra, and Wang's leading coefficients along random curves.
complexity: Polynomial in the degrees and the number of terms in practice (randomised); the extension degree is $O(\log_p)$ of the number of terms, and the Frobenius recombination costs $O(k)$ coefficient powerings per factor for an extension of degree $k$.
order: 75
---

FLINT factors a polynomial $A \in \mathbb{F}_q[x_1, \ldots, x_n]$, $q =
p^d$ with a word-size prime $p$, in the same stages as over the integers
(see [](algorithms/multivariate/factorisation)): content removal,
squarefree (separable) decomposition, reduction to a bivariate problem by
evaluation, bivariate factorisation, and lifting of the bivariate factors
to all variables by Hensel lifting [@Wan1978] [@MonTun2020] or recombination
[@Zas1969]. Two things differ in positive characteristic. A polynomial can
have zero derivative without being constant, which the squarefree
decomposition must handle; and the field can be too small to provide the
evaluation points that the reduction to two variables needs. The second
problem is solved by factoring over an extension field and recombining
the factors found there along the orbits of the Frobenius automorphism,
an idea that goes back to von zur Gathen and Kaltofen's analysis of
factoring over finite fields. The implementation follows the notes of its
author [@Sch-mpoly].

## The pipeline {#pipeline}

`nmod_mpoly_factor_algo` (and `fq_nmod_mpoly_factor_algo`) proceeds as
follows. Throughout, factors are made monic and the leading coefficient
of $A$ is kept as the constant of the factorisation.

1. *Content.* `nmod_mpoly_factor_content` removes the monomial content
   (powers of variables dividing $A$) and splits $A$ into factors that are
   primitive with respect to every variable, by repeatedly computing the
   content with respect to one variable.
2. *Compression.* Each factor is compressed by an invertible affine change
   of the exponent lattice (`mpoly_compression_set`) so that it uses as
   few variables and as small exponents as possible. A randomised test on
   projections of its Newton polytope may prove the factor irreducible at
   this point (an indecomposable Newton polytope implies irreducibility
   [@Gao2001]).
3. *Separable decomposition* (`nmod_mpoly_factor_separable` with `sep`
   = 1): see [](#separable). Each resulting factor $g$ is squarefree and
   satisfies, for every variable, $\partial g/\partial x_i = 0$ or
   $\gcd(g, \partial g / \partial x_i) = 1$, with at least one nonzero
   derivative.
4. *Deflation by powers of $p$* ([](#deflation)) and a permutation that
   puts a variable with $\gcd(g, \partial g/\partial x_0) = 1$ first.
5. *Irreducible factorisation* of each factor: univariate factorisation
   (`nmod_poly_factor`) for one variable; bivariate factorisation
   (`n_bpoly_mod_factor_smprime`, else `n_bpoly_mod_factor_lgprime`) for
   two; for more variables, polynomials of degree 1 in $x_0$ are
   irreducible, those of degree 2 are split by solving a quadratic
   equation with `nmod_mpoly_quadratic_root` ([](#quadratic)), and the
   others are factored by the lifting algorithms of [](#irreducible).
6. The factors are inflated, uncompressed, and multiplied by the
   multiplicities of the separable decomposition.

## Separable factorisation in characteristic $p$ {#separable}

In characteristic zero, Yun's algorithm [@Yun1976] computes the squarefree
decomposition from GCDs of $A$ with its derivatives. In characteristic $p$
the derivative $\partial A/\partial x_i$ vanishes on the terms whose
exponent in $x_i$ is a multiple of $p$, so factors of multiplicity
divisible by $p$, and factors that are polynomials in $x_i^p$, are
invisible to it [@GT1996] [@Lec2008]. The following lemma deals with what
remains.

::: lemma $p$-th roots {#lem-proot}
Let $K$ be a perfect field of characteristic $p$, and $F \in K[x_1,
\ldots, x_n]$ with $\partial F/\partial x_i = 0$ for all $i$. Then $F =
H^p$, where $H$ is obtained from $F$ by dividing every exponent by $p$
and replacing every coefficient $c$ by $c^{1/p}$. For $K = \mathbb{F}_q$,
$c^{1/p} = c^{p^{d-1}}$, and over $\mathbb{F}_p$, $c^{1/p} = c$.
:::

::: proof
The coefficient of $x^e / x_i$ in $\partial F/\partial x_i$ is $e_i c_e$,
so $\partial F/\partial x_i = 0$ means that $p \mid e_i$ whenever $c_e \ne
0$. Thus $F = \sum_e c_e x^{e}$ with all exponents divisible by $p$. As
the Frobenius map $u \mapsto u^p$ is a ring endomorphism in characteristic
$p$, $H^p = \sum_e (c_e^{1/p})^p x^{e} = F$. In $\mathbb{F}_q$,
$c^{p^d} = c$, so $c^{1/p} = c^{p^{d-1}}$.
:::

`_nmod_mpoly_factor_separable` works as follows on a monic primitive $A$.
For each variable in turn (choosing first the variable whose derivative
has the fewest terms), it runs Yun's algorithm with respect to that
variable on the current cofactor $C$: with $G = \partial C/\partial x$,
$C \gets \gcd(C, G)$, $W = C_{\mathrm{old}}/C$, $V = G/C$, and while $V -
\partial W/\partial x \ne 0$ and the multiplicity $k$ is less than $p -
1$, the factor of multiplicity $k$ is $\gcd(W, V - \partial W/\partial
x)$. The factors found are split further against the derivatives with
respect to the remaining variables, to obtain the property of step 3. When
all variables have been processed, the remaining $C$ has zero derivative
in every variable; with $p^k$ the largest power of $p$ dividing all its
exponents, $C = D^{p^k}$ by [](#lem-proot), $D$ is computed by deflation
(and, over $\mathbb{F}_q$, by the coefficient map $c \mapsto c^{p^{j}}$
with $j \equiv -k \pmod d$), its separable decomposition is computed
recursively and raised to the power $p^k$, and the two decompositions are
merged into one with pairwise coprime bases by GCDs.

## Deflation by powers of $p$ {#deflation}

After the separable decomposition, a factor $g$ may still have all its
exponents in some variables divisible by a power of $p$, for example $g =
x^p - y$. FLINT divides the exponents of each variable $x_i$ by the
largest power $p^{s_i}$ of $p$ dividing all of them, factors the deflated
polynomial $h$ (so that $g = h(x_1^{p^{s_1}}, \ldots, x_n^{p^{s_n}})$)
and inflates the factors. This is correct because of the following lemma.

::: lemma {#lem-pdefl}
Let $K$ be a perfect field of characteristic $p$, and $f \in K[x, y_1,
\ldots, y_m]$ irreducible with $\partial f/\partial y_j \ne 0$ for some
$j$. Then $f(x^p, y_1, \ldots, y_m)$ is irreducible.
:::

::: proof
Write $F = f(x^p, y)$ and let $\sigma$ be the automorphism of $K[x, z_1,
\ldots, z_m]$ that maps each coefficient $c$ to $c^{1/p}$ and fixes the
variables. In $K[x, z]$ with $z_i^p = y_i$, the Frobenius gives $F(x, z^p)
= \sigma(f)(x, z)^p$, and $\sigma(f)(x, z)$ is irreducible. Suppose $F =
GH$ in $K[x, y]$. Then $G(x, z^p) H(x, z^p) = \sigma(f)(x, z)^p$, so by
unique factorisation $G(x, z^p) = c\, \sigma(f)(x, z)^a$ for a unit $c$
and some $0 \le a \le p$. The left side has zero derivative with respect
to $z_j$, while the derivative of the right side is $a c\,
\sigma(f)^{a-1} \partial \sigma(f)/\partial z_j$, where $\partial
\sigma(f) / \partial z_j = \sigma(\partial f/\partial y_j)(x, z) \ne 0$.
Hence $p \mid a$, that is $a \in \{0, p\}$, and $G$ or $H$ is a unit.
:::

::: corollary {#cor-pdefl}
Let $g$ be squarefree, and $h$ obtained by dividing the exponents of
some variables by powers of $p$ as above, with at least one variable
undeflated. If $h = \prod_i f_i$ is the irreducible factorisation, so is
$g = \prod_i f_i(x_1^{p^{s_1}}, \ldots, x_n^{p^{s_n}})$.
:::

::: proof
Each $f_i$ has nonzero derivative with respect to some undeflated
variable: otherwise $f_i$ is a polynomial in the $p$-th powers of the
undeflated variables, and its inflation is a polynomial in $p$-th powers
of all variables, hence a $p$-th power by [](#lem-proot), contradicting
squarefreeness of $g$. Deflating one variable and one power of $p$ at a
time, [](#lem-pdefl) shows that each inflated $f_i$ stays irreducible.
:::

## Quadratics {#quadratic}

Let $A = a X^2 + b X + c$ with $a, b, c \in R = \mathbb{F}_q[x_1, \ldots,
x_n]$ be squarefree and primitive in $X = x_0$. It is reducible if and
only if it has a factor of degree 1 in $X$. Since $a A = Y^2 + b Y + ac$
with $Y = aX$, and the roots $Y = a r$ of this monic polynomial, $r$ a
root of $A$, are integral over the unique factorisation domain $R$, the
polynomial $A$ is reducible if and only if $Y^2 + bY + ac$ has a root $y \in
R$; then $A$ is the product of $(a/g) X - y/g$ and $g X - cg/y$ with $g
= \gcd(a, y)$. In odd characteristic $y = (-b + \sqrt{b^2 - 4ac})/2$, and
the square root is computed by the heap algorithm of
[](algorithms/multivariate/division#sqrt). In characteristic 2 there is no
such formula, and FLINT solves the equation $Y^2 + AY = B$ (signs do
not matter) directly, term by term, by an algorithm from Schultz's notes
[@Sch-mpoly] that simplifies §§4–5 of Cherly, Gallardo, Vaserstein and
Wheland [@CGVW1998], whose method divides by polynomials and does not
preserve sparsity.

Let $A, B \in R = \mathbb{F}_{2^k}[x_1, \ldots, x_n]$ with $A, B \ne 0$
(if $A = 0$ the root is $\sqrt{B}$, if $B = 0$ it is 0), and write
$\mathrm{lt}(A) = a M$ with $a \in \mathbb{F}_{2^k}$ and $M$ a monomial. If
$Y_0$ is a root, so is $Y_0 + A$, and the coefficients of $M$ in $Y_0$ and
$Y_0 + A$ differ by $a$, so at most one of them equals $a$: we may look
for a root in which the coefficient of $M$ is not $a$.

::: lemma Leading term of a root {#lem-quad-lt}
Let $Y_0$ be a root of $Y^2 + AY = B$ in which the coefficient of $M$ is
not $a$, and let $m = \gamma\, x^{\mu}$ be its leading term. Then exactly
one of the following holds:

1. $\mathrm{lm}(B) > M^2$: then $x^\mu > M$ and $m = \sqrt{\mathrm{lt}(B)}$,
   which requires $\mathrm{lm}(B)$ to be a square;
2. $\mathrm{lm}(B) = M^2$: then $x^\mu = M$ and $m = a \zeta M$ with $\zeta^2 +
   \zeta = \mathrm{lc}(B)/a^2$;
3. $\mathrm{lm}(B) < M^2$: then $x^\mu < M$ and $m = \mathrm{lt}(B) /
   \mathrm{lt}(A)$, which requires this quotient to be a monomial
   smaller than $M$.

Moreover $Y_0 + m$ is a root of $Y^2 + AY = B + m^2 + Am$, the leading
monomial of $B + m^2 + Am$ is smaller than that of $B$, and the
coefficient of $M$ in $Y_0 + m$ is not $a$ (in case 2 it is zero).
:::

::: proof
In characteristic 2, $Y_0^2$ is the sum of the squares of the terms of
$Y_0$, so its leading term is $m^2 = \gamma^2 x^{2\mu}$, while the leading
term of $A Y_0$ is $a \gamma M x^\mu$. If $x^\mu > M$, the first dominates
and $\mathrm{lt}(B) = m^2$; if $x^\mu < M$, the second does and
$\mathrm{lt}(B) = \mathrm{lt}(A) m$. If $x^\mu = M$, the two have the same
monomial and their sum is $(\gamma^2 + a\gamma) M^2$, which is nonzero
because $\gamma \ne 0, a$; it equals $\mathrm{lt}(B)$, and $\zeta = \gamma/a$
satisfies $\zeta^2 + \zeta = \mathrm{lc}(B)/a^2$. In the three cases
$\mathrm{lm}(B)$ is $x^{2\mu} > M x^{\mu} > M^2$, $M^2$, or $M x^\mu < M^2$,
so the case is determined by $\mathrm{lm}(B)$, and $m$ is determined by
$\mathrm{lt}(B)$ (in case 2 up to the two roots $\zeta$, $\zeta + 1$).
Finally $(Y_0 + m)^2 + A(Y_0 + m) = B + m^2 + Am$, and $m^2 + Am$ has
leading term $\mathrm{lt}(B)$ by the above.
:::

::: algorithm Quadratic equations in characteristic 2 {#alg-quad2}
Input: nonzero $A, B \in \mathbb{F}_{2^k}[x_1, \ldots, x_n]$.
Output: $Y$ with $Y^2 + AY = B$, or "no root".

1. Set $Y = 0$ and $E = B$.
2. While $E \ne 0$: determine $m$ from $\mathrm{lt}(E)$ by the case of
   [](#lem-quad-lt) that applies, using case 1 only while all terms of $Y$
   are larger than $M$, and case 3 once a term $\le M$ has been chosen; if
   the required square root, quadratic root or monomial quotient does not
   exist, return "no root". Set $Y \gets Y + m$ and $E \gets E + m^2 + Am$.
3. Return $Y$.
:::

By [](#lem-quad-lt) and induction, if a root exists the algorithm finds
one, term by term; if it returns $Y$, then $E = 0$ and $Y$ is a root.
The leading monomial of $E$ decreases strictly, so the terms of $Y$ form
a strictly decreasing sequence of monomials, which is finite because a
monomial ordering is a well-ordering; so the algorithm terminates also
when there is no root. The terms of $E$ are produced in decreasing order
by a heap, as in [](algorithms/multivariate/division#heap), merging $B$,
the products $A m$ and the squares $m^2$ (no cross terms appear in
characteristic 2), so the heap has at most $\ell_A + 2$ entries and the
cost is that of a heap division with $\ell_A \ell_Y$ products.

Case 2 needs a root of $\zeta^2 + \zeta = c$ in $\mathbb{F}_{2^k}$.

::: lemma Artin–Schreier equations in characteristic 2 {#lem-trace}
Let $\mathrm{Tr}(c) = \sum_{i=0}^{k-1} c^{2^i}$ be the absolute trace of
$\mathbb{F}_{2^k}$. The equation $\zeta^2 + \zeta = c$ has a root if and
only if $\mathrm{Tr}(c) = 0$, and then, for any $u$ with $\mathrm{Tr}(u) =
1$,
$$
\zeta = \sum_{i=1}^{k-1} c^{2^i} \sum_{j=0}^{i-1} u^{2^j}
$$
is a root. If $\mathbb{F}_{2^k} = \mathbb{F}_2[\theta]/(P)$ with $k \ge 2$,
then $u = 1/(\theta P'(\theta))$ has trace 1.
:::

::: proof
The map $\zeta \mapsto \zeta^2 + \zeta$ is $\mathbb{F}_2$-linear with kernel
$\mathbb{F}_2$, so its image has $2^{k-1}$ elements; it lies in the kernel
of the trace, because $\mathrm{Tr}(\zeta^2) = \mathrm{Tr}(\zeta)$, and the
trace is onto, so the image is the kernel of the trace [@LN1997, Ch. 2].
For the formula, let $U_i = \sum_{j < i} u^{2^j}$, so that $U_i^2 = U_{i+1}
+ u$, $U_1 = u$ and $U_k = \mathrm{Tr}(u) = 1$. Then $\zeta^2 = \sum_{i=2}^{k}
c^{2^i}(U_i + u)$, and in $\zeta^2 + \zeta$ the terms with $2 \le i \le k-1$
cancel, leaving $c^{2^k} U_k + c^2 U_1 + u \sum_{i=2}^{k} c^{2^i} = c + u
\,\mathrm{Tr}(c) = c$. For $u$: the partial fraction decomposition
$1/(xP(x)) = P(0)^{-1}/x + \sum_\rho (\rho P'(\rho))^{-1}/(x - \rho)$ over
the roots $\rho$ of $P$, multiplied by $x$ and compared at $x \to \infty$,
gives $\sum_\rho 1/(\rho P'(\rho)) = -1/P(0) = 1$, and the left side is the
trace of $1/(\theta P'(\theta))$.
:::

`fq_nmod_mpoly_quadratic_root` and `nmod_mpoly_quadratic_root` (to which
`fmpz_mod_mpoly_quadratic_root` reduces for a word-size modulus)
implement this: in odd characteristic
by a square root, in characteristic 2 by [](#alg-quad2) with the static
heap function `_fq_nmod_mpoly_quadratic_root_heap` and the trace formula
of [](#lem-trace). Over $\mathbb{F}_2$ itself case 2 never applies, since
$\mathrm{Tr}(1) = 1$. The factoring code uses them for squarefree
factors of degree 2 in the main variable (`_apply_quadratic`).

## Evaluation points and extension fields {#irreducible}

To factor a squarefree $A$ in $n \ge 3$ variables, the lifting algorithms
choose a point $\alpha = (\alpha_2, \ldots, \alpha_{n-1})$ and factor the
bivariate image $A(x_0, x_1, \alpha)$. The point must keep the leading
coefficient in $x_0$ nonzero and the image separable in $x_0$, and Wang's
method needs in addition that the factors of the leading coefficient stay
distinguishable. These conditions fail on the zero set of nonzero
polynomials, so over a large field random points work with high
probability ([](algorithms/multivariate/modp-gcd#lem-sz)). Over a small
field they may all fail, and the bivariate factorisation itself needs a
value $y = \alpha$ for which the univariate image is separable.

FLINT therefore tries, in order: points in $\mathbb{F}_p$ (`smprime`),
points in an extension $\mathbb{F}_{p^k}$ with $p^k \le 10^6$, where field
arithmetic uses Zech logarithms (`medprime`, via `fq_zech_mpoly`), and
points in extensions of any degree (`lgprime`, via `fq_nmod_mpoly`). The
extension degree starts at about $\max(2, 1 + \log_p \ell)$ for $\ell$
terms (half that for Wang's and Zassenhaus's methods) and is increased
each time the algorithm runs out of points. In the bivariate algorithm,
`n_bpoly_mod_factor_lgprime` reduces $y$ modulo irreducible polynomials
$m(y)$ of increasing degree and lifts $m$-adically instead of
$(y-\alpha)$-adically.

A factorisation over $\mathbb{F}_{q^k}$ is turned into one over
$\mathbb{F}_q$ as follows.

::: theorem Frobenius orbits {#thm-frob}
Let $f \in \mathbb{F}_q[x_1, \ldots, x_n]$ be irreducible and separable,
that is, $\gcd(f, \partial f/\partial x_i) = 1$ for some $i$ with
$\partial f/\partial x_i \ne 0$, and let $\phi$ be the automorphism $c
\mapsto c^q$ of $L = \mathbb{F}_{q^k}$ over $\mathbb{F}_q$, acting on
$L[x_1, \ldots, x_n]$ through the coefficients. Then the monic
irreducible factors of $f$ in $L[x_1, \ldots, x_n]$ are pairwise
distinct and form a single orbit $\{g, \phi(g), \ldots, \phi^{r-1}(g)\}$
under $\phi$, with $r \mid k$, and $f = \mathrm{lc}(f) \prod_{j=0}^{r-1}
\phi^j(g)$.
:::

::: proof
$f$ remains separable over $L$, because $\gcd(f, \partial f/\partial x_i)
= 1$ over $L$ by [](algorithms/multivariate/modp-gcd#lem-ext); so $f$ is
squarefree over $L$ and its monic irreducible factors $g_1, \ldots, g_s$
are distinct. Since $\phi$ fixes $f$ and is a ring automorphism preserving
monicity, it permutes the $g_i$. Let $O$ be the orbit of $g_1$ and $P$ the
product of the polynomials in $O$. Then $\phi(P) = P$, so $P$ has
coefficients in $\mathbb{F}_q$ (the fixed field of $\phi$), and $P$
divides $f$ in $\mathbb{F}_q[x_1, \ldots, x_n]$ (division by a monic
polynomial is the same over $\mathbb{F}_q$ and over $L$). As $f$ is
irreducible over $\mathbb{F}_q$, $P = f/\mathrm{lc}(f)$, so all factors
lie in $O$. The orbit length divides the order $k$ of $\phi$.
:::

Hence, for a squarefree $A$ factored over $\mathbb{F}_q$ via its factors
$g_1, \ldots, g_s$ over $\mathbb{F}_{q^k}$, the irreducible factors over
$\mathbb{F}_q$ are the products over the orbits. FLINT computes them by
taking a factor $t$, applying $\phi$ to its coefficients $k - 1$ times,
removing each conjugate $\phi^j(t)$ from the list of remaining factors,
and multiplying the orbit together; the product has coefficients in
$\mathbb{F}_q$ and is mapped back.

## Bivariate factorisation {#bivariate}

Over $\mathbb{Q}$, an irreducible bivariate polynomial stays irreducible
under a generic specialisation $y = y_0$, so a univariate image usually has
as many factors as the polynomial itself and recombination is cheap
([](algorithms/multivariate/factorisation#bivariate)). Over a finite
field this fails: a univariate image of an irreducible polynomial
typically splits into several factors, and recombination needs a better
method than trying subsets [@BHKS2009]. `n_bpoly_mod_factor_smprime`
factors $A \in \mathbb{F}_p[y][x]$, separable in $x$, as follows:

1. Remove the content in $\mathbb{F}_p[y]$.
2. For $\alpha = 0, 1, 2, \ldots$ with $\deg_x A(x, \alpha) = \deg_x A$,
   factor $A(x, \alpha)$ with `nmod_poly_factor`, record the degrees of
   the factors for the degree pruning of
   [](algorithms/multivariate/factorisation#lem-prune) (which may prove
   $A$ irreducible), and skip $\alpha$ if the image is not squarefree. Of
   the first two good points, keep the one with fewer factors.
3. Shift $y \gets y + \alpha$, make $A$ monic in $x$ as a power series in
   $y$, and lift the $r$ monic local factors to $g_1, \ldots, g_r \in
   \mathbb{F}_p[[y]][x]$ modulo $y^\lambda$.
4. Find the true factors as products of the $g_i$ by linear algebra over
   $\mathbb{F}_p$, increasing $\lambda$ until this succeeds.

The linear algebra is the function field version of van Hoeij's method
[@vHo2002; @BHKS2009], with the logarithmic derivative bounds used by
Lecerf [@Lec2006]. For $A$ of degree $d$ in $x$, let $H$ be the smallest
concave function on $[0, d]$ with $H(i) \ge \deg_y [x^i] A$ whenever $[x^i]
A \ne 0$, the upper boundary of the Newton polygon of $A$ with respect to
$\deg_y$ (`mpoly_bivar_cld_bounds`).

::: lemma Logarithmic derivative bounds {#lem-cld}
If $F \in K[y][x]$ divides $A$, then $\deg_y [x^i] \bigl(A \,
\partial_x F / F\bigr) \le H(i + 1)$ for $0 \le i < d$.
:::

::: proof
Write $A = FG$, so that $A\, \partial_x F/F = G\, \partial_x F$. Newton
polygons with respect to the valuation $-\deg_y$ are additive under
multiplication (Dumas' theorem, the one-variable case of Ostrowski's
theorem of [](algorithms/multivariate/factorisation#polytope)), so $H$ is
the sup-convolution $H(n) = \max_{a + b = n} (H_G(a) + H_F(b))$ of the
corresponding functions for $G$ and $F$. Since $[x^b] \partial_x F = (b+1)
[x^{b+1}] F$, the coefficient of $x^i$ in $G\, \partial_x F$ has degree in
$y$ at most $\max_{a + b = i} (H_G(a) + H_F(b + 1)) \le H(i + 1)$.
:::

Let $\ell_i = (A / g_i) \, \partial_x g_i \bmod y^\lambda$. If $F$ is a
true factor with $F \equiv \operatorname{lc}_x(F) \prod_{i \in S} g_i$,
the logarithmic derivative of the product gives $\sum_{i \in S} \ell_i
\equiv A\, \partial_x F / F$, so by [](#lem-cld) the coefficients of $y^j
x^i$ in $\sum_i e_i \ell_i$ vanish for $H(i + 1) < j < \lambda$ when $e \in
\{0, 1\}^r$ is the indicator vector of $S$. These coefficients are linear
forms in $e$ over $\mathbb{F}_p$, and the indicator vectors of the
irreducible factors of $A$ lie in their common kernel. FLINT keeps a basis
of this kernel in reduced row echelon form, starting from the identity
matrix and intersecting with the kernel of the new equations
(`nmod_mat_init_nullspace_tr`, tested with `nmod_mat_is_reduced`) as the
precision grows:

- if the kernel has dimension 1, $A$ is irreducible, since the indicator
  vectors of its irreducible factors are linearly independent and lie in
  the kernel;
- if the reduced basis consists of 0-1 vectors with exactly one 1 in each
  column, its rows partition the local factors, and the candidate factors
  are the products over the rows, multiplied by $\operatorname{lc}_x(A)$,
  reduced modulo $y^{\deg_y A + 1}$ and made primitive; they are tested
  by division, and if a test fails, products of two candidates are tried
  as well;
- otherwise the precision is increased.

When groups of local factors are known to belong together, they are
multiplied together and the lifting continues with fewer factors. Each
returned factor is verified by division, so the result is correct
whatever the precision. `n_bpoly_mod_factor_lgprime` does the same with
$y - \alpha$ replaced by an irreducible polynomial $m(y)$ of increasing
degree, for fields without good points: $A$ is factored over
$\mathbb{F}_p[y]/(m)$ and the factors are lifted $m$-adically. The
`fq_nmod` versions `n_fq_bpoly_factor_smprime` and
`n_fq_bpoly_factor_lgprime` are analogous.

## Leading coefficients over a finite field {#lcc}

Wang's method ([](algorithms/multivariate/hensel-lifting#wang)) assigns
the factors of $\ell = \operatorname{lc}_{x_0}(A)$ to the factors of $A$
by looking at prime numbers that divide the values $L_j(\alpha)$ of the
factors of $\ell$. Over a finite field the values are field elements and
there are no primes; moreover a univariate image of an irreducible
polynomial is usually reducible, so lifting from a univariate
factorisation would require recombination. FLINT therefore replaces the
evaluation point by a *curve*: each variable $x_i$, $i \ge 1$, is
replaced by a random polynomial $x_i = \alpha_i + \beta_i Y + \gamma_i Y^2 +
\cdots$ in a new variable $Y$, and $A(x_0, \alpha(Y))$ is factored as a
bivariate polynomial in $\mathbb{F}_p[Y][x_0]$ [@Sch-mpoly]. Irreducible
polynomials in $\mathbb{F}_p[Y]$ play the role of the primes, and the
proof of [](algorithms/multivariate/hensel-lifting#lem-wang) applies
with $\Z$ replaced by $\mathbb{F}_p[Y]$: if an irreducible $q(Y)$ divides
$L_j(\alpha(Y))$ but neither the content of the bivariate image nor the
other $L_k(\alpha(Y))$, the exponent of $L_j$ in the leading coefficient
of the $i$-th factor is the multiplicity of $q$ in the leading
coefficient of the $i$-th bivariate factor times the content, divided by
the multiplicity of $q$ in $L_j(\alpha(Y))$.

`nmod_mpoly_factor_lcc_wang` checks the sequential form of the condition
(each monic $L_j(\alpha(Y))$ must keep a factor of positive degree after
removing its common factors with the content and with $L_1(\alpha(Y)),
\ldots, L_{j-1}(\alpha(Y))$) and computes the exponents with
`n_poly_mod_remove`. The lifting itself starts from the factors at $Y = 0$,
which are the images of the bivariate factors at the point $(\alpha_1,
\ldots, \alpha_n)$; since they come from a bivariate factorisation, they
correspond to the factors of $A$ much more often than the factors of an
arbitrary univariate image would. If the bivariate image is irreducible,
so is $A$, because the substitution preserves the degree in $x_0$.
`nmod_mpoly_factor_irred_smprime_wang` tries a few random curves through
each random point $\alpha$ (one at which the degrees are preserved and
the univariate image is squarefree), raising the degree of the curves in
$Y$ from 1 up to 9 as attempts fail, before giving up.

## Choosing among the methods {#choice}

For a factor with at least three variables after deflation,
`nmod_mpoly_factor_algo` with all methods enabled first factors the
leading coefficient in $x_0$ (needed by Wang's and Zippel's methods). It
then computes the density $\ell / \binom{D + n}{n}$ of the factor ($\ell$
terms, total degree $D$, $n$ variables) and tries in order:

- density above 0.005: Wang over $\mathbb{F}_p$, Wang over
  $\mathbb{F}_{p^k}$ (`medprime`), sparse lifting over $\mathbb{F}_p$, sparse
  lifting `medprime`; otherwise the same with sparse lifting first;
- density above 0.001: Wang `lgprime`, then sparse lifting `lgprime`;
  otherwise the other order;
- Zassenhaus's method over $\mathbb{F}_p$, `medprime` and `lgprime`, as a
  last resort (and first if the leading coefficient cannot be factored).

The first method that succeeds determines the factorisation. Over
$\mathbb{F}_q$, `fq_nmod_mpoly_factor_algo` has the same structure with
only the `smprime` and `lgprime` levels, the extensions being presented
through the embeddings of `fq_nmod_mpoly.h` and the orbits taken under $c
\mapsto c^q$. `fmpz_mod_mpoly_factor_algo` calls `nmod_mpoly_factor_algo`
for a word-size modulus and uses only the `smprime` level otherwise.

## Use in FLINT {#use}

- Public functions: `nmod_mpoly_factor`, `nmod_mpoly_factor_squarefree`,
  `nmod_mpoly_factor_separable`, `nmod_mpoly_factor_algo` and its
  variants; the same for `fq_nmod_mpoly`, `fmpz_mod_mpoly` and
  `fq_zech_mpoly` (which converts to `fq_nmod_mpoly`).
- Separable decomposition: `_nmod_mpoly_factor_separable`,
  `_fq_nmod_mpoly_factor_separable` ([](#separable)).
- Bivariate: `n_bpoly_mod_factor_smprime` and `n_bpoly_mod_factor_lgprime`,
  `n_fq_bpoly_factor_smprime` and `n_fq_bpoly_factor_lgprime`, with
  Hensel lifting and linear-algebra recombination ([](#lem-cld)), using
  `mpoly_bivar_cld_bounds`.
- Quadratics: `nmod_mpoly_quadratic_root`, `fq_nmod_mpoly_quadratic_root`
  ([](#alg-quad2), [](#lem-trace)).
- Leading coefficients: `nmod_mpoly_factor_lcc_wang` ([](#lcc)), used by
  `nmod_mpoly_factor_irred_smprime_wang`.
- Lifting: `nmod_mpoly_factor_irred_{sm,med,lg}prime_{zassenhaus,wang,zippel}`
  and the `fq_nmod_mpoly` and `fq_zech_mpoly` analogues; Frobenius
  recombination ([](#thm-frob)) in the `medprime` and `lgprime` versions.

::: functions nmod_mpoly_factor nmod_mpoly_factor_squarefree nmod_mpoly_factor_separable fq_nmod_mpoly_factor fmpz_mod_mpoly_factor :::
