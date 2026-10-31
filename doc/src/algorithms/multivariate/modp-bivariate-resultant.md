---
title: Bivariate resultants by evaluation and interpolation
summary: The resultant $\mathrm{Res}_y(A, B) \in \mathbb{F}_p[x]$ of two bivariate polynomials, computed from univariate resultants at many values of $x$ (a geometric progression or roots of unity) and interpolation, with the degree bound, the correction for specialisations that lower the degree, the reduction of more variables to two by Kronecker substitution, and FLINT's cutoffs.
complexity: $O((n_A + n_B) M(N) + N M(n_A) \log n_A)$ operations in $\mathbb{F}_p$ for $N = n_B d_A + n_A d_B + 1$ points, where $n_A, n_B$ are the degrees in the eliminated variable and $d_A, d_B$ those in the other one.
order: 80
---

Let $A, B \in \mathbb{F}_p[x][y]$ have degrees $n_A = \deg_y A \ge 1$ and
$n_B = \deg_y B \ge 1$ in $y$ and degrees $d_A$, $d_B$ in $x$. Their
resultant with respect to $y$ is the determinant of the Sylvester matrix,
a polynomial $R = \mathrm{Res}_y(A, B) \in \mathbb{F}_p[x]$. The
subresultant algorithm computes it by pseudo-remainder sequences in
$\mathbb{F}_p[x][y]$, whose intermediate coefficients have growing degree
in $x$. When $\mathbb{F}_p$ is large enough, it is faster to evaluate $x$
at many points, compute univariate resultants over $\mathbb{F}_p$, and
interpolate, as in the classical modular resultant algorithm [@GG2003,
§6.9]. FLINT implements this method for `nmod_mpoly` with fast evaluation
and interpolation at structured points: a geometric progression
[@BS2005], or roots of unity when $p$ is an FFT prime.

## Degree bound and specialisation {#theory}

::: lemma Degree bound {#lem-deg}
$\deg_x R \le n_B d_A + n_A d_B$.
:::

::: proof
The Sylvester matrix has $n_B$ rows of coefficients of $A$ (polynomials
in $x$ of degree at most $d_A$) and $n_A$ rows of coefficients of $B$
(degree at most $d_B$). Every term of the expansion of the determinant
is a product of one entry from each row.
:::

Hence $R$ is determined by its values at $N = n_B d_A + n_A d_B + 1$
distinct points, the number Collins [@Col1971] takes for each variable in
turn. The value at $x = a$ must be computed from the
specialisations $f = A(a, y)$ and $g = B(a, y)$, whose degrees may be
lower than $n_A$ and $n_B$ when $a$ is a root of a leading coefficient.

::: lemma Specialisation {#lem-spec}
Let $a \in \mathbb{F}_p$, $f = A(a, y)$, $g = B(a, y)$, and $\ell_A$,
$\ell_B$ the leading coefficients of $A$ and $B$ in $y$.

1. If $\deg f = n_A$, then $R(a) = \ell_A(a)^{n_B - \deg g}
   \mathrm{Res}(f, g)$.
2. If $\deg g = n_B$ and $\deg f < n_A$, then $R(a) =
   (-1)^{n_B (n_A - \deg f)} \ell_B(a)^{n_A - \deg f} \mathrm{Res}(f, g)$.
3. If both degrees drop, $R(a) = 0$.

Here $\mathrm{Res}(f, g)$ is the resultant of $f$ and $g$ with respect to
their actual degrees, and is 0 if $f$ or $g$ is zero.
:::

::: proof
Evaluation commutes with the determinant, so $R(a)$ is the determinant of
the Sylvester matrix of $f$ and $g$ built with the formal degrees $n_A$
and $n_B$. In case 1, the leading $n_B - \deg g$ coefficients of $g$ in
that matrix are zero, and expanding along the first columns, each of
which then has the single nonzero entry $\ell_A(a)$, gives the formula.
Case 2 follows from case 1 applied to $\mathrm{Res}_y(B, A) =
(-1)^{n_A n_B} R$, with the sign of the exchange for the actual degrees.
In case 3, the first column of the matrix contains only $\ell_A(a) =
\ell_B(a) = 0$.
:::

## Algorithm {#algorithm}

::: algorithm Bivariate resultant by multipoint evaluation {#alg-res}
Input: $A, B$ with $n_A \ge n_B \ge 1$ and a prime $p$.
Output: $R = \mathrm{Res}_y(A, B)$, or failure if $\mathbb{F}_p$ is too
small.

1. Set $N = n_B d_A + n_A d_B + 1$.
2. Choose $N$ distinct points: if FLINT has `fft_small` and $p - 1$ is
   divisible by $2^k \ge \max(N, 16)$ (an FFT prime), the powers of a primitive $2^k$-th root of unity $w$;
   otherwise a geometric progression $1, q, q^2, \ldots, q^{N-1}$ with
   $q = r^2$ of order at least $N$, where $r$ is the first of $2, 3, 4,
   \ldots$ that works. If $p < 2N + 1$, no such $q$ exists: fail.
3. Evaluate each of the $n_A + n_B + 2$ coefficients of $A$ and $B$ (as
   polynomials in $x$) at the points, by FFTs or by the geometric
   progression algorithm [@BS2005], in blocks of points sized to bound
   the memory.
4. At each point, compute the univariate resultant of the specialisations
   with `_nmod_poly_resultant` and apply [](#lem-spec).
5. Interpolate $R$ from its $N$ values, by an inverse FFT or geometric
   interpolation.
:::

The evaluation of one coefficient of length at most $\max(d_A, d_B) + 1$
at $N$ points of a geometric progression costs $O(M(N))$ operations, so
step 3 costs $O((n_A + n_B) M(N))$, step 5 costs $O(M(N))$, and step 4
costs $N$ univariate resultants of degree $n_A$, each $O(M(n_A) \log
n_A)$ with the half-gcd algorithm ([](algorithms/polynomials/half-gcd))
or $O(n_A n_B)$ with Euclid's algorithm for small degrees. In the
FFT case the points are split into $2^{\mathrm{depth} - m}$ cosets of the
$2^m$-th roots of unity, each evaluated by a substitution $x \mapsto w^t
x$ and one transform of length $2^m$; the values come out in the order
that the truncated inverse transform needs, so interpolation is a single
inverse transform (or a geometric interpolation with $r = \sqrt{w}$ when
that is cheaper). Both loops over the points run in parallel on the
threads set with `flint_set_num_threads` when the work is large enough.

## More variables: Kronecker substitution {#kronecker}

Let now $A, B \in \mathbb{F}_p[x_1, \ldots, x_k][y]$, of degrees $n_A$,
$n_B \ge 1$ in $y$ and $d_{A,j}$, $d_{B,j}$ in $x_j$. A ring homomorphism
$\phi$ of the coefficient ring that keeps the leading coefficients in $y$
nonzero commutes with the resultant, since the resultant is a polynomial
in the coefficients: $\phi(\mathrm{Res}_y(A, B)) = \mathrm{Res}_y(\phi(A),
\phi(B))$. Three such maps reduce the problem to [](#alg-res).

1. *Monomial factors.* If $A = x^a A'$, then $\mathrm{Res}_y(A, B) =
   x^{a n_B} \mathrm{Res}_y(A', B)$, and similarly for $B$; both are
   shifted down to their minimal exponent in each $x_j$.
2. *Deflation.* If every exponent of $x_j$ left in $A$ and $B$ is a
   multiple of $g_j$, the map $x_j \mapsto x_j^{g_j}$ is injective and its
   inverse is applied to both.
3. *Kronecker substitution* [@Kro1882]. With $D_j$ one more than a
   bound on the degree in $x_j$ of the resultant of the shifted and
   deflated polynomials, the map $x_j \mapsto t^{S_j}$ with $S_j = D_1
   \cdots D_{j-1}$ is injective on the polynomials of degree less than
   $D_j$ in each $x_j$.

The resultant is then recovered from its image in $\mathbb{F}_p[t]$ by
reading off the digits of each exponent of $t$ in the mixed radix $(D_1,
\ldots, D_k)$, and [](#alg-res) needs $D_1 \cdots D_k$ points: exactly one
per coefficient of the dense box of the resultant. The packed
coefficients of $A$ and $B$ hold at most as many entries, so the memory
used is of the order of the size of that box. The bounds come from the
following refinement of [](#lem-deg).

::: lemma Weighted degree bound {#lem-wdeg}
Let $w \ge 0$, and let $a(w)$ be the largest value of $e + w m$ over the
terms $x_j^e y^m$ of $A$ (the other variables being ignored), and $b(w)$
that over the terms of $B$. Then
$$
\deg_{x_j} \mathrm{Res}_y(A, B) \le n_B\, a(w) + n_A\, b(w) - w\, n_A n_B.
$$
:::

::: proof
Row $i < n_B$ of the Sylvester matrix holds the coefficient of $y^{c-i}$
in $A$ in column $c$, a polynomial of degree at most $a(w) + w i - w c$ in
$x_j$, and row $n_B + i$ holds that of $B$, of degree at most $b(w) + w i
- w c$. A term of the determinant takes one entry in each row and each
column, so its degree is at most the sum of these over the rows less $w$
times the sum of the column indices, which is the right-hand side.
:::

The weight $w = 0$ gives $n_B d_{A,j} + n_A d_{B,j}$ as in
[](#lem-deg), which is reached by polynomials dense in a box. The weight
$w = 1$ gives the Bézout bound $d_A d_B$ when $A$ and $B$ have total
degrees $d_A$ and $d_B$ in $(x_j, y)$ and the same degrees in $y$: half as
much, and the box of the resultant of two such polynomials in three
variables has four times fewer coefficients. The right-hand side is a
convex, piecewise linear function of $w$, and $a(w)$ only depends on the
upper convex hull of the points $(m, e)$, so its minimum is at $w = 0$ or
at a slope $w$ of that hull for $A$ or for $B$; FLINT takes the minimum
over these. Since $b(w) \ge w n_B$, the bound is at least $n_B d_{A,j}$
and $n_A d_{B,j}$, so the substitution is also injective on $A$ and $B$,
and keeps their leading coefficients in $y$ nonzero.

The images of $A$ and $B$ have degree at most $\sum_j d_{A,j} S_j$ and
$\sum_j d_{B,j} S_j$ in $t$, so that the bound of [](#lem-deg) on the
degree of their resultant is
$$
\sum_j (n_B d_{A,j} + n_A d_{B,j}) S_j + 1 \ge \sum_j (D_j - 1) S_j + 1
  = D_1 \cdots D_k,
$$
with equality when every $D_j$ comes from the weight 0, and
[](#alg-res) is told the number of points $D_1 \cdots D_k$. Applied to
the substitution $x_j \mapsto s x_j$ for all $j$ at once,
[](#lem-wdeg) also bounds the total degree of the resultant in the $x_j$,
which tells how much of the box the resultant can fill: about half of it
for inputs dense in total degree with two variables packed, a sixth with
three.

Over $\mathbb{Z}$ the substitution changes neither the coefficients nor
their norms, so the multimodular algorithm and its coefficient bound
apply to the images unchanged.

Collins [@Col1971] eliminates the variables one at a time instead: each
$x_j$ in turn is evaluated at $n_B d_{A,j} + n_A d_{B,j} + 1$ points,
down to univariate resultants, which are interpolated back one variable
at a time. Both take one univariate resultant per point of the box. The
substitution turns the nested evaluations and interpolations into a
single one in $t$, which [](#alg-res) does with fast algorithms on
structured points. Hu and Monagan [@HuMon2016] pack the variables of a
sparse GCD problem the same way, leaving the main variable alone. For
GCDs some substitutions are unlucky and have to be detected. There are
none for the resultant: it commutes with every substitution that keeps
the leading coefficients in $y$ nonzero, so being injective on the
resultant is enough.

## Example {#example}

Two random polynomials of degree 8 in each variable modulo the prime
$10^6 + 3$; `nmod_mpoly_resultant` takes the evaluation and interpolation
path, and the result agrees with the subresultant algorithm applied to the
univariate forms. The resultant reaches the degree bound of
[](#lem-deg):

```c example name=nmod_mpoly_bivariate_resultant
#include <flint/flint.h>
#include <flint/nmod_mpoly.h>

int main(void)
{
    nmod_mpoly_ctx_t ctx;
    nmod_mpoly_t a, b, r1, r2;
    nmod_mpoly_univar_t ua, ub;
    ulong bounds[2] = {9, 9};     /* degrees < 9 in x and in y */
    flint_rand_t state;

    flint_rand_init(state);
    nmod_mpoly_ctx_init(ctx, 2, ORD_LEX, UWORD(1000003));
    nmod_mpoly_init(a, ctx); nmod_mpoly_init(b, ctx);
    nmod_mpoly_init(r1, ctx); nmod_mpoly_init(r2, ctx);
    nmod_mpoly_univar_init(ua, ctx); nmod_mpoly_univar_init(ub, ctx);

    nmod_mpoly_randtest_bounds(a, state, 81, bounds, ctx);
    nmod_mpoly_randtest_bounds(b, state, 81, bounds, ctx);

    /* eliminate y (variable 1): evaluation/interpolation in x */
    nmod_mpoly_resultant(r1, a, b, 1, ctx);

    /* the same resultant by subresultants over F_p[x] */
    nmod_mpoly_to_univar(ua, a, 1, ctx);
    nmod_mpoly_to_univar(ub, b, 1, ctx);
    nmod_mpoly_univar_resultant(r2, ua, ub, ctx);

    flint_printf("deg_y a = %wd, deg_x a = %wd\n", nmod_mpoly_degree_si(a, 1, ctx),
                 nmod_mpoly_degree_si(a, 0, ctx));
    flint_printf("deg_x res = %wd (bound %wd), equal: %d\n",
                 nmod_mpoly_degree_si(r1, 0, ctx),
                 nmod_mpoly_degree_si(b, 1, ctx) * nmod_mpoly_degree_si(a, 0, ctx)
                 + nmod_mpoly_degree_si(a, 1, ctx) * nmod_mpoly_degree_si(b, 0, ctx),
                 nmod_mpoly_equal(r1, r2, ctx));

    nmod_mpoly_univar_clear(ua, ctx); nmod_mpoly_univar_clear(ub, ctx);
    nmod_mpoly_clear(a, ctx); nmod_mpoly_clear(b, ctx);
    nmod_mpoly_clear(r1, ctx); nmod_mpoly_clear(r2, ctx);
    nmod_mpoly_ctx_clear(ctx);
    flint_rand_clear(state);
    return 0;
}
```
```output
deg_y a = 8, deg_x a = 8
deg_x res = 128 (bound 128), equal: 1
```

## Use in FLINT {#use}

- `nmod_mpoly_resultant_kronecker` applies [](#kronecker) and then
  [](#alg-res) to polynomials in two or more variables with single-word
  exponents and positive degrees in the eliminated variable, packing the
  terms directly into the dense bivariate form; it declines when the dense
  form would have more than $2^{32}$ entries or points, or when
  $\mathbb{F}_p$ is too small ([](#alg-res) fails).
- `nmod_mpoly_resultant` uses it when the dense form does not exceed the
  input by more than 64 times beyond $2^{26}$ entries and a gain is
  predicted. With a single variable packed into $t$ (the bivariate case),
  the inputs must not be tiny and `n_bpoly_mod_resultant_cutoff` must
  predict a gain: $(n_A + 1)(n_B + 1) \ge 36$, or $\ge 16$ with at least
  4096 points. With two or more, both inputs must have at least 5 terms
  and $(T^2/6)^{L - 2} (N'/N)^{16} \ge e^{5.5} N^{1/2}$, where $T$ is the
  smaller number of terms, $L$ the smaller length in $y$, $N$ the number
  of points and $N' \le N$ the number of monomials of the box within the
  bound on the total degree: the subresultant algorithm pays for the terms
  of the resultant, and the dense algorithm for its whole box. Otherwise it
  uses the subresultant algorithm of Ducos [@Duc2000] on the univariate
  forms (`nmod_mpoly_univar_resultant`).
- `_n_bpoly_mod_resultant_multipoint` implements [](#alg-res) on
  `n_bpoly_t`, taking the number of points from the caller. It is reached
  through `_n_bpoly_mod_resultant_len` (and `n_bpoly_mod_resultant`,
  `_n_bpoly_mod_resultant`), where the bivariate algorithm is chosen by the
  sizes for every dense bivariate resultant over $\mathbb{F}_p$: those of
  the Kronecker forms above, of the images modulo primes below and of the
  generic polynomial rings. [](#alg-res) is the only algorithm there so
  far; one added there, an algorithm of lower complexity at large degrees
  for instance [@Vil2018], is used by all of these in its range.
- `fmpz_mpoly_resultant_kronecker` applies [](#kronecker) and then
  `_fmpz_bpoly_resultant_len`, which uses `_n_bpoly_mod_resultant_len`
  modulo many primes and combines the results by the Chinese remainder
  theorem (see [](fmpz_mpoly)), up to $2^{30}$ entries and $2^{28}$
  points.
  `fmpz_mpoly_resultant` uses it under the same conditions as over
  $\mathbb{F}_p$, except that with two or more variables packed it needs
  at least 4 terms and $(T^2/6)^{L - 2} (N'/N)^{16} \ge e^{8}$, the
  number of points mattering little against the coefficient growth of the
  subresultant algorithm. `fmpq_mpoly_resultant` goes through `fmpz_mpoly_resultant`.
- `gr_poly_resultant` over a polynomial ring over `nmod` passes its
  arguments to `_n_bpoly_mod_resultant` when `n_bpoly_mod_resultant_cutoff`
  predicts a gain, and `gr_poly_resultant_multipoint` to
  `_n_bpoly_mod_resultant_multipoint` (the memory layouts of `gr_poly`
  over `nmod` and of `n_bpoly` coincide).

::: functions nmod_mpoly_resultant nmod_mpoly_resultant_kronecker n_bpoly_mod_resultant n_bpoly_mod_resultant_cutoff :::
