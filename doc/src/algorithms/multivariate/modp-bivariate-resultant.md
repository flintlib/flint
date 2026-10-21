---
title: Bivariate resultants by evaluation and interpolation
summary: The resultant $\mathrm{Res}_y(A, B) \in \mathbb{F}_p[x]$ of two bivariate polynomials, computed from univariate resultants at many values of $x$ (a geometric progression or roots of unity) and interpolation, with the degree bound, the correction for specialisations that lower the degree, and FLINT's cutoffs.
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
distinct points. The value at $x = a$ must be computed from the
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

- `nmod_mpoly_resultant` uses the algorithm for polynomials in exactly two
  variables with positive degrees in the eliminated variable, after
  converting to dense bivariate form (`nmod_mpoly_get_bpoly`), when the
  dense boxes are not too large (at most $2^{32}$ entries and points, and
  not more than 64 times the number of terms beyond $2^{26}$ entries), the
  inputs are not tiny, and `n_bpoly_mod_resultant_cutoff` predicts a gain:
  $(n_A + 1)(n_B + 1) \ge 36$, or $\ge 16$ with at least 4096 points.
  Otherwise, or if $\mathbb{F}_p$ is too small ([](#alg-res) fails), it
  uses the subresultant algorithm of Ducos [@Duc2000] on the univariate
  forms (`nmod_mpoly_univar_resultant`).
- `n_bpoly_mod_resultant` and `_n_bpoly_mod_resultant` implement
  [](#alg-res) on `n_bpoly_t`.
- `fmpz_mpoly_resultant` uses `_n_bpoly_mod_resultant` modulo many primes
  and combines the results by the Chinese remainder theorem (see
  [](fmpz_mpoly)).
- `gr_poly_resultant_multipoint`, for polynomials over $\mathbb{F}_p[x]$
  in the generic rings framework, passes its arguments to
  `_n_bpoly_mod_resultant` (the memory layouts of `gr_poly` over `nmod`
  and of `n_bpoly` coincide).

::: functions nmod_mpoly_resultant n_bpoly_mod_resultant n_bpoly_mod_resultant_cutoff :::
