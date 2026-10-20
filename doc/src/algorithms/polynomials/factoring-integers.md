---
title: Factoring polynomials over the integers
summary: The Zassenhaus algorithm (factor modulo a prime, Hensel lift, recombine subsets of the modular factors) and van Hoeij's knapsack method with CLD data and gradual feeding to LLL, as implemented in fmpz_poly_factor, with the reductions that come first (content, powers of x, deflation with Capelli certificates, squarefree factorisation), the choice of the prime and the bounds that make the result provable.
complexity: Zassenhaus recombination tries up to $2^{r-1}$ subsets of the $r$ modular factors; van Hoeij's algorithm with gradual feeding runs in time polynomial in the degree and the coefficient size, and in practice its cost is dominated by the modular factorisation and the Hensel lifting.
order: 70
---

Factoring a polynomial $f \in \Z[x]$ means writing it as $f = c
\prod_i f_i^{e_i}$ with $c \in \Z$ and distinct primitive irreducible
$f_i$ of positive degree. By Gauss' lemma this is the same as
factoring over $\mathbb{Q}$, up to the content. All practical methods
follow the plan of Zassenhaus [@Zas1969]: factor $f$ modulo a small
prime $p$, which is easy ([](algorithms/polynomials/factoring-finite-fields)),
lift the modular factorisation to one modulo a large power $p^a$
([](algorithms/polynomials/hensel-lifting)), and find out which
products of modular factors are true factors. The last step, the
*recombination*, is the difficult one: a polynomial with $r$ modular
factors has $2^r$ candidate subsets. Lenstra, Lenstra and Lovász showed
that factoring is possible in polynomial time [@LLL1982], but their
algorithm was not competitive. Van Hoeij [@vHo2002] turned the
recombination into a knapsack problem solved by lattice reduction,
which is fast in practice and, with the refinements of Belabas, van
Hoeij, Klüners and Steel [@BHKS2009], van Hoeij and Novocin
[@vHN2010], and Hart, van Hoeij and Novocin [@HvHN2011], provably
polynomial. The implementation in FLINT is due to Novocin and Hart and
follows [@HvHN2011] and [@Nov2008].

Throughout, $f$ has degree $n$ and leading coefficient $\ell =
\operatorname{lc}(f) > 0$, and $\|f\|_2$ is the Euclidean norm of its
coefficient vector.

## Reductions {#reductions}

`fmpz_poly_factor` first reduces to the factorisation of squarefree
primitive polynomials with nonzero constant term:

1. The content is removed, with the sign of $\ell$. Polynomials of
   degree at most 3 are handled directly: a quadratic splits if and
   only if its discriminant is a square, and a cubic is reducible if
   and only if it has a rational root, which `_fmpz_poly_factor_cubic`
   locates from real and 2-adic approximations of the roots
   ([](#cubic)).
2. The factor $x^k$ is split off.
3. If $f = g(x^d)$ with $d > 1$, the polynomial $g$ of degree $n/d$ is
   factored, and the inflation $h(x^d)$ of each irreducible factor $h$
   of $g$ is handled separately ([](#deflation)).
4. Otherwise `fmpz_poly_factor_squarefree` computes the squarefree
   factorisation $f = \prod_i g_i^i$ by Yun's algorithm [@Yun1976]
   ([](algorithms/polynomials/squarefree)), with the GCDs of
   [](algorithms/polynomials/gcd-integers), and each $g_i$ is factored
   by `_fmpz_poly_factor_zassenhaus`.

### Deflation and Capelli's theorem {#deflation}

If $h$ is irreducible, $h(x^d)$ need not be: $x - 4$ is irreducible
but $x^2 - 4$ is not. For prime exponents, Capelli's theorem decides
the question [@Sch2000].

::: theorem Capelli {#thm-capelli}
Let $h \in \mathbb{Q}[x]$ be irreducible, $\theta$ a root of $h$, $K =
\mathbb{Q}(\theta)$ and $p$ a prime. Then $h(x^p)$ is irreducible over
$\mathbb{Q}$ if and only if $\theta$ is not a $p$-th power in $K$.
:::

FLINT applies the prime case along a chain $h(x^d) = U(x^{d/p})$ with
$U = h(x^p)$. A *witness* that $\theta$ is not a $p$-th power is a
prime $q$ and a residue field $\mathbb{F}_{q^f}$ of $K$ above $q$ in
which the image of $\theta$ is not a $p$-th power, which by Euler's
criterion means $\theta^{(q^f - 1)/p} \ne 1$ when $p \mid q^f - 1$.
`_fmpz_poly_factor_inflation_is_irreducible_capelli` searches for a
witness with a budget linear in $\deg h$; if it finds one, $h(x^p)$ is
irreducible, and otherwise $h(x^p)$ is factored by the general
algorithm at degree $p \deg h$, which is still smaller than factoring
$f$ from scratch.

### Cubics {#cubic}

`_fmpz_poly_factor_cubic` looks for rational roots instead of factoring
modulo a prime, by a method from Schultz's notes [@Sch-mpoly]. Let $f =
a x^3 + b x^2 + c x + d$ with $a > 0$. With $X = 3a x + b$,
$$
27 a^2 f = X^3 - 3 A X - B, \qquad A = b^2 - 3ac, \qquad B = (9ac - 2b^2) b -
27 a^2 d,
$$
and a rational root of the monic $g = X^3 - 3AX - B \in \Z[X]$ is an
integer. Each integer root $r$ of $g$ gives the factor
$\operatorname{pp}(3a x + b - r)$ of $f$. The integer roots of $g$ are
found as follows (the cases $B = 0$ and $B = \pm 1$ are immediate).

1. *Small primes.* If $A$ and $B$ are odd, $g \equiv X^3 + X + 1 \pmod 2$
   has no root, so $f$ is irreducible. If $g$ has no root modulo one of
   the primes $5, 7, \ldots, 29$, $f$ is irreducible; this is tested at
   once by computing $\prod_{x=0}^{28} g(x)$ modulo the product of these
   primes, which fits in a word.
2. *Discriminant.* Let $D = B^2 - 4A^3$, so that the discriminant of $g$
   is $-27D$. If $D = 0$, the roots are $B/A$ and the double root
   $-B/(2A)$. If $D > 0$, $g$ has one real root, $\bigl(\sqrt[3]{4B + 4\sqrt
   D} + \sqrt[3]{4B - 4\sqrt D}\bigr)/2$ by Cardano's formula; it is
   computed with integer square and cube roots, and the three integers
   around the approximation are tested.
3. *Three real roots* ($D < 0$, hence $A > 0$). With $X = 2\sqrt{A} \cos
   \theta$ the equation becomes $2 A^{3/2} \cos 3\theta = B$, so all roots
   satisfy $|r| \le 2\sqrt A$. After a test modulo the primes $31, \ldots,
   47$, the roots are computed in the 2-adic integers $\Z_2$ to a
   precision $2^k$ with $2^{k - 1} > 2\sqrt A$, reduced to the symmetric
   range, and tested exactly.

For the 2-adic roots write $-3A = 2^{\alpha} a'$ and $-B = 2^{\beta} b'$ with
$a', b'$ odd, so that $g = X^3 + 2^\alpha a' X + 2^\beta b'$. The Newton
polygon of $g$ at 2 gives three cases.

- $2\beta = 3\alpha$: with $X = 2^{\beta/3} Y$, $g = 2^\beta (Y^3 + a' Y +
  b')$, which has no root modulo 2; $g$ has no root in $\Z_2$, and $f$ is
  irreducible.
- $2\beta < 3\alpha$: all roots have valuation $\beta/3$, so if $3 \nmid
  \beta$ there is none in $\Z_2$. Otherwise $X = 2^{\beta/3} Y$ gives
  $2^\beta (Y^3 + 2^{\alpha - 2\beta/3} a' Y + b')$, which is $(Y + 1)(Y^2 + Y
  + 1)$ modulo 2 with an irreducible quadratic factor, so there is
  exactly one root in $\Z_2$, found by Hensel lifting.
- $2\beta > 3\alpha$: with $X = 2^{\beta - \alpha} Y$ and $e = 2\beta - 3\alpha >
  0$, $g = 2^\beta (2^e Y^3 + a' Y + b')$, which has the simple root $Y \equiv
  1$ modulo 2. Newton iteration, doubling the precision at each step,
  lifts the factorisation
  $$
  2^e Y^3 + a' Y + b' = (Y + r)(2^e Y^2 - 2^e r Y + s), \qquad r, s \in
  \Z_2 \text{ odd},
  $$
  that is, $g = (X + 2^{\beta - \alpha} r)(X^2 - 2^{\beta - \alpha} r X +
  2^\alpha s)$. The roots of the quadratic have valuation $\alpha/2$. If
  $\alpha$ is odd, $-2^{\beta - \alpha} r$ is the only root in $\Z_2$. If
  $\alpha$ is even, the other roots are
  $$
  2^{\beta - \alpha - 1} r \pm 2^{\alpha/2} \sqrt{2^{2\beta - 3\alpha - 2} r^2 -
  s},
  $$
  with a 2-adic square root (an odd 2-adic integer is a square if and
  only if it is $1$ modulo 8). The square root loses precision, so $r$
  and $s$ are lifted further when the result is not precise enough.

Once an integer root $r_0$ has been found, the other roots are those of
$X^2 + r_0 X + r_0^2 - 3A$, rational if and only if $3(4A - r_0^2)$ is a
square. The method costs a bounded number of integer square and cube
roots and 2-adic Newton iterations on numbers of the size of the
coefficients, so it is much faster than the general algorithm for
cubics with large coefficients.

## The modular factorisation {#prime}

Let $f$ be squarefree and primitive with $f(0) \ne 0$. A prime $p$ is
*admissible* if $p \nmid \ell$, $p \nmid f(0)$ and $f \bmod p$ is
squarefree; all but finitely many primes are, since $f \bmod p$ is
squarefree when $p$ divides neither $\ell$ nor the discriminant. For an
admissible prime,

$$
f \equiv \ell \, g_1 g_2 \cdots g_r \pmod p
$$

with distinct monic irreducible $g_i \in \mathbb{F}_p[x]$, and every
irreducible factor $g$ of $f$ in $\Z[x]$ reduces to $\operatorname{lc}(g)
\prod_{i \in S} g_i$ for a subset $S$, the subsets of distinct factors
being disjoint. In particular $f$ has at most $r$ irreducible factors,
and it is irreducible if $r = 1$.

`_fmpz_poly_factor_zassenhaus` factors $f$ modulo the first three
admissible primes, starting from 2, and keeps the one with the fewest
factors, since the cost of recombination grows with $r$. The factor
degrees modulo all three primes are combined: a true factor of degree
$k$ is a product of modular factors of total degree $k$ modulo each
prime, so $k$ must be a subset sum of the modular degrees for all
three. The `zassenhaus_prune` functions maintain the set of such
degrees; if it contains no degree strictly between 0 and $n$, the
polynomial is irreducible, and otherwise subsets of an impossible
degree are skipped during recombination. By Chebotarev's density
theorem, the factor degrees modulo a random prime are the cycle lengths
of a random element of the Galois group, so for a polynomial with
Galois group $S_n$, the typical case, $r$ is about $\log n$ on
average.

## Bounds and precision {#bounds}

Recombination recovers a factor from its image modulo $p^a$, which
requires a bound for its coefficients.

::: theorem Mignotte {#thm-mignotte}
If $g = \sum_{j=0}^m b_j x^j \in \Z[x]$ divides $f$ in $\Z[x]$, then
for all $j$

$$
|b_j| \le \binom{m-1}{j} \|f\|_2 + \binom{m-1}{j-1} |\operatorname{lc}(f)|,
$$

and in particular $\|g\|_\infty \le 2^m \|f\|_2$ [@Mig1974],
[@Coh1996, Theorem 3.5.1].
:::

`fmpz_poly_factor_mignotte` computes $B$, the maximum of this bound
over $1 \le j < n$ with $\|f\|_2$ rounded up, of $|f(0)|$ and of
$\ell$, so that every proper factor of $f$ has coefficients bounded by
$B$ (the extreme coefficients of a factor divide those of $f$).

::: lemma Recovering a factor {#lem-recover}
Let $g \in \Z[x]$ be a primitive divisor of $f$ with positive leading
coefficient, let $f \equiv \ell \prod_{i=1}^r g_i \pmod{p^a}$ with
monic $g_i$, and suppose that $g \equiv \operatorname{lc}(g) \prod_{i
\in S} g_i \pmod{p^a}$. If $p^a > 2 \ell B$, then $g$ is the primitive
part of the polynomial $\tilde g$ with coefficients in $(-p^a/2,
p^a/2]$ that is congruent to $\ell \prod_{i \in S} g_i$ modulo $p^a$.
:::

::: proof
The polynomial $(\ell / \operatorname{lc}(g))\, g$ has integer
coefficients, because $\operatorname{lc}(g) \mid \ell$, and it is
congruent to $\ell \prod_{i \in S} g_i$ modulo $p^a$, since
$\operatorname{lc}(g)$ is invertible modulo $p$. Its coefficients are
bounded by $\ell B < p^a/2$ in absolute value, so it is equal to
$\tilde g$, and its primitive part is $g$.
:::

The Zassenhaus algorithm therefore lifts to the least $a$ with $p^a
\ge 2 \ell B + 1$.

## Zassenhaus recombination {#zassenhaus}

::: algorithm Zassenhaus recombination {#alg-zassenhaus}
**Input:** squarefree primitive $f$, monic $g_1, \ldots, g_r$ with $f
\equiv \ell \prod g_i \pmod{p^a}$, where $p^a > 2 \ell B$.

**Output:** the irreducible factors of $f$.

1. $T \gets \{1, \ldots, r\}$, $F \gets f$, $k \gets 1$.
2. While $2k \le |T|$: for each $S \subseteq T$ with $|S| = k$ and a
   possible total degree, let $h$ be the primitive part of the
   symmetric residue of $\operatorname{lc}(F) \prod_{i \in S} g_i$
   modulo $p^a$; if $h \mid F$, output $h$ and set $F \gets F/h$ and $T
   \gets T \setminus S$. Then $k \gets k + 1$.
3. Output $F$ if $\deg F > 0$.
:::

::: theorem {#thm-zassenhaus}
[](#alg-zassenhaus) outputs the irreducible factors of $f$.
:::

::: proof
Every irreducible factor of $F$ corresponds to a subset of $T$, and the
subsets of distinct factors partition $T$. A factor $h$ found for a
subset $S$ divides $F$; if it were reducible, its irreducible factors
would correspond to proper subsets of $S$, which were tried earlier
and would have been removed from $F$. So every output is irreducible.
Since $F \mid f$, the bound $B$ applies to the factors of $F$, and
$\operatorname{lc}(F) \mid \ell$; by [](#lem-recover) (with $F$ and
$\operatorname{lc}(F)$ in place of $f$ and $\ell$) the factor of $F$
with subset $S$ is found when $S$ is tried. When all subsets of size at
most $|T|/2$ have been tried, $F$ is irreducible, since a reducible $F$
would have a factor with at most half of the remaining modular
factors.
:::

The number of subsets tried is at most $\sum_{k \le r/2}
\binom{r}{k} \le 2^{r-1}$, and each test is a trial division, made
cheap by `fmpz_poly_divides`, which first compares the constant
coefficients and the values at 1. This is fast when $r$ is small, the
typical case, and hopeless when $r$ is large. FLINT uses Zassenhaus
recombination for $r \le 8$.

## Van Hoeij's algorithm {#van-hoeij}

Van Hoeij's idea [@vHo2002] is to represent subsets by indicator
vectors $v \in \{0, 1\}^r$ and to find the vectors of the true factors
as short vectors of a lattice. Let $g_1, \ldots, g_r$ be monic with
$f \equiv \ell \prod g_i \pmod{p^a}$, and for $S \subseteq \{1,
\ldots, r\}$ let $g_S = \prod_{i \in S} g_i$.

### Logarithmic derivatives {#cld}

The map $g \mapsto f g'/g$ turns products into sums:

$$
\frac{f g_S'}{g_S} = \sum_{i \in S} \frac{f g_i'}{g_i}.
$$

The *coefficients of the logarithmic derivative* (CLD) are the
coefficients of the polynomials $f g_i' / g_i$, which have degree less
than $n$ and are well defined modulo $p^a$, since $g_i$ divides $f$
there. If $g$ is a true factor of $f$, then $f g'/g = (f/g) g' \in
\Z[x]$, and its coefficients are much smaller than $p^a$: on $|x| =
\rho$, $f g'/g = f \sum_\alpha 1/(x - \alpha)$ over the roots $\alpha$
of $g$, and splitting the coefficients of $f$ at index $j$ gives the
bound computed by `fmpz_poly_CLD_bound` [@HvHN2011]: with $f = \sum
a_i x^i$ and any $\rho > 0$, the coefficient of $x^j$ in $f g'/g$ is at
most

$$
n \max\left( \rho^{-(j+1)} \sum_{i \le j} |a_i| \rho^i,\;
\rho^{-(j+1)} \sum_{i > j} |a_i| \rho^i \right),
$$

which is minimised numerically over $\rho$. For the $k$ lowest and
$k$ highest coefficients $j$, let $c_{ij}$ be the coefficient of $x^j$
in $f g_i'/g_i$ reduced into $(-p^a/2, p^a/2]$
(`_fmpz_poly_factor_CLD_mat`). If $v$ is the indicator vector of a true
factor $g$, then $\sum_i v_i c_{ij} \equiv [x^j]\, f g'/g \pmod{p^a}$
is small. The low coefficients are computed with `fmpz_poly_mullow`
and `fmpz_poly_divlow_smodp`, the high ones with
`fmpz_poly_mulhigh_n` and `fmpz_poly_divhigh_smodp`, so that only
$O(k)$ coefficients of each product are needed.

### The knapsack lattice {#lattice}

For a set $J$ of CLD columns, let $L$ be the lattice spanned by the
rows of

$$
\begin{pmatrix} 2^U I_r & C \\ 0 & p^a I_{|J|} \end{pmatrix},
\qquad C = (c_{ij})_{i \le r,\, j \in J},
$$

with each column of $C$ scaled down so that the values produced by
true factors become small. The vector $(2^U v, \ldots)$ of a true
factor then has length at most about $2^U \sqrt{r + 1}$, while a
random 0-1 combination produces entries of the size of $p^a$. The
vectors of the true factors span a sublattice $W$ whose rank is the
number of irreducible factors, and the reduced echelon basis of $W$
consists of the indicator vectors of the factors, a *0-1 basis*.
Lattice reduction finds a basis of a lattice $L' \supseteq W$ of
decreasing rank, by the following lemma.

::: lemma Short vectors after reduction {#lem-lll-removal}
Let $b_1, \ldots, b_m$ be a basis of a lattice with Gram–Schmidt
vectors $b_1^*, \ldots, b_m^*$, and suppose $\|b_i^*\| > \beta$ for all
$i > k$. Then every lattice vector of norm at most $\beta$ lies in the
span of $b_1, \ldots, b_k$.
:::

::: proof
Write $v = \sum_{i \le j} c_i b_i$ with integers $c_i$ and $c_j \ne 0$.
The component of $v$ orthogonal to $b_1, \ldots, b_{j-1}$ is $c_j
b_j^*$, so $\|v\| \ge |c_j| \|b_j^*\| \ge \|b_j^*\|$. If $\|v\| \le
\beta$, then $\|b_j^*\| \le \beta$, so $j \le k$.
:::

After each LLL reduction, `fmpz_lll_with_removal_ulll` drops the final
basis vectors whose Gram–Schmidt norm exceeds the bound $\beta^2 = (r
+ 1) 2^{2U}$ ([LLL with removals](algorithms/linear-algebra/lll#removals), with the
truncation strategy [](algorithms/linear-algebra/lll#ulll));
by the lemma, $W$ survives. Each new column of CLD data
can only lower the rank of $L'$, and enough columns at enough precision
make $L' = W$ [@vHo2002], [@HvHN2011].

::: algorithm Van Hoeij's algorithm with gradual feeding {#alg-van-hoeij}
**Input:** squarefree primitive $f$ of degree $n$ and its factorisation
modulo an admissible prime $p$ into $r$ monic factors.

**Output:** the irreducible factors of $f$.

1. Choose an initial precision $a$ (see below) and lift the factors to
   precision $p^a$. Set $M \gets 2^U I_r$ with $U = \lceil \log_2
   (\max(r, 20) + 1) \rceil$.
2. If the first $r$ columns of $M$ define a partition of $\{1, \ldots,
   r\}$, form the corresponding products; if they all divide $f$,
   output them and stop.
3. Compute the CLD data for $k$ low and $k$ high coefficients, leaving
   out coefficients whose CLD bound is too large compared to $p^a$.
4. For each new column, taken alternately from the low and the high
   end: append it to $M$, keeping only the bits of the entries above
   the CLD bound (`fmpz_mat_next_col_van_hoeij`); if this makes the
   matrix worth reducing, run LLL with removal, and go to step 2.
5. When the columns are used up, double $k$ (up to about $n/2$); when
   all coefficients have been used, double the precision $a$, continue
   the Hensel lift, and go to step 2.
:::

The test in step 2 is `fmpz_poly_factor_van_hoeij_check_if_solved`:
it reads a partition from the echelon form of the first $r$ columns
(`fmpz_mat_col_partition`), forms the products of the lifted factors
in each part, multiplied by $\ell$ and reduced modulo $p^a$, and
checks that their primitive parts divide $f$, first modulo 2 and then
over $\Z$. Since the answer is verified by division, the output is
always correct; the bounds only ensure that the algorithm makes
progress and terminates.

*Gradual feeding*, adding one column at a time, keeps each LLL call on
a matrix whose new column has only a few more bits than the reduced
basis, which is what makes the total work polynomial [@vHN2010],
[@HvHN2011]. The starting precision is about the least $a$ with

$$
p^a \ge 2^{2.5 r + b} \sqrt{n + 1},
$$

where $b$ is the bit size of the smaller of the CLD bounds of the
lowest and the highest coefficient, but not more than the Zassenhaus
precision. This is usually far below the Mignotte bound and enough for
the first few columns to determine the factorisation; otherwise the
precision is doubled and the Hensel lift continued with
`_fmpz_poly_hensel_continue_lift`. The first rounds use $k = 30$ (or
50 for $r > 200$) coefficients from each end instead of 10 when $3r >
n + 1$.

## Complexity {#complexity}

For a polynomial of degree $n$ with $b$-bit coefficients, the
Mignotte precision is $a \log p = O(n + b)$ bits, and the Hensel lift
to that precision costs $O(M(n) \log r)$ operations on integers of
that size ([](algorithms/polynomials/hensel-lifting)). Zassenhaus
recombination adds up to $2^{r-1}$ trial divisions. Van Hoeij's
algorithm with gradual feeding runs in time polynomial in $n$ and $b$
[@HvHN2011], and usually needs much less than the Mignotte precision;
the LLL work is small unless $r$ is large and there are many factors
of similar degree. The classical worst cases of the Zassenhaus
algorithm, such as the Swinnerton-Dyer polynomials of degree $n = 2^k$,
which have $n/2$ factors modulo every prime, are easy: proving
$S_7(x + 1)$, of degree 128, irreducible takes about 0.1 seconds on
the x86-64 machine used for this documentation.

## Use in FLINT {#use}

- `fmpz_poly_factor` performs the reductions of [](#reductions) and
  calls `_fmpz_poly_factor_zassenhaus` with cutoff 8 and van Hoeij's
  algorithm enabled; `fmpz_poly_factor_zassenhaus` uses only the
  Zassenhaus recombination.
- `_fmpz_poly_factor_zassenhaus` chooses the prime, lifts with
  `fmpz_poly_hensel_lift_once` to the Mignotte precision and calls
  `fmpz_poly_factor_zassenhaus_recombination_with_prune`, or calls
  `fmpz_poly_factor_van_hoeij` when $r > 8$.
- `fmpz_poly_factor_van_hoeij` uses `_fmpz_poly_hensel_start_lift`,
  `_fmpz_poly_hensel_continue_lift`, `_fmpz_poly_factor_CLD_mat`,
  `fmpz_mat_next_col_van_hoeij`, `fmpz_lll_with_removal_ulll` (with
  $\delta = 0.75$ and $\eta = 0.81$; see [](algorithms/linear-algebra/lll))
  and `fmpz_poly_factor_van_hoeij_check_if_solved`.
- `fmpz_poly_factor_squarefree`, `fmpz_poly_factor_mignotte`,
  `fmpz_poly_CLD_bound`, `_fmpz_poly_factor_inflation_is_irreducible_capelli`,
  `_fmpz_poly_factor_quadratic` and `_fmpz_poly_factor_cubic` provide
  the remaining steps.
- Multivariate polynomials over $\Z$ are factored by `fmpz_mpoly_factor`,
  which reduces to univariate factorisations of this kind
  ([](algorithms/multivariate/overview)); polynomials over number fields
  and finite fields have their own modules.

::: functions fmpz_poly_factor fmpz_poly_factor_zassenhaus _fmpz_poly_factor_zassenhaus fmpz_poly_factor_van_hoeij fmpz_poly_factor_squarefree fmpz_poly_factor_zassenhaus_recombination_with_prune _fmpz_poly_factor_inflation_is_irreducible_capelli :::
