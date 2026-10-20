---
title: Sparse interpolation
summary: Recovering a sparse polynomial from its values at a few points: Zippel's method for a known support (transposed Vandermonde systems solved in quadratic time, with extra unknowns for scaled images), and the method of Ben-Or and Tiwari for an unknown support (Berlekamp–Massey, roots of the recurrence polynomial, discrete logarithms modulo a smooth prime after Kronecker substitution) with the Berlekamp–Massey algorithm computed incrementally by half-gcd steps, as used by FLINT's sparse GCD and factoring algorithms.
complexity: With a known support of $t$ monomials, $t$ evaluations and $O(t^2)$ field operations to solve for the coefficients; with an unknown support, $2t$ evaluations, $O(t^2)$ operations for Berlekamp–Massey, root finding of a degree-$t$ polynomial and $t$ discrete logarithms.
order: 50
---

The sparse GCD and factoring algorithms of FLINT reduce a problem in many
variables to many problems in one or two variables, whose solutions are
the values of the unknown polynomial (a GCD, a cofactor, a lifted factor)
at chosen points. *Sparse interpolation* reconstructs a polynomial with
$t$ terms from $O(t)$ such values, independently of the number of
variables and of the degrees, which would determine the cost of dense
interpolation. The two classical approaches are Zippel's, which needs the
support (set of monomials) of the result and finds it variable by variable
[@Zip1979; @Zip1981], and Ben-Or and Tiwari's, which finds support and
coefficients together from twice as many values [@BenOrTiw1988;
@KalLak1988]. Both work over a finite field $\F_p$ with $p$ large enough
for the points chosen to be good with high probability; over $\Z$ the
coefficients are then lifted or combined over several primes.

All methods need strict degree bounds $\deg_{x_j} f < r_j$. The number of
values of $f$ (*probes*) they use, for $f$ with $t$ terms in $m$
variables, compares as follows; the presentation below, in particular
of the Berlekamp–Massey algorithm as a half-gcd computation, follows
Schultz's notes [@Sch-mpoly].

| Method | Probes | Main difficulties |
|---|---|---|
| dense, variable by variable (Newton) | $\prod_j r_j$ | enough evaluation points |
| Zippel ([](#zippel-shape)) | about $t \sum_j r_j$ | enough points; inconsistent or underdetermined systems when the shape assumption fails; the linear algebra |
| Ben-Or–Tiwari ([](#bma)) | $2t$ | an element of order at least $\prod_j r_j$; the linear algebra, root finding and discrete logarithms |

For Zippel's method the count comes from adding one variable $x_j$ at a
time: each of the $r_j$ images needed for the new variable is found from
as many probes as the current support has terms, at most $t$.

Throughout, $f = \sum_{i=1}^{t} c_i x^{e_i} \in \F_p[x_1, \ldots, x_m]$
with distinct $e_i$ and nonzero $c_i$. For a point $\alpha \in (\F_p^*)^m$
write $\mu_i = \alpha^{e_i}$ for the value of the $i$-th monomial.

## Known support: transposed Vandermonde systems {#zippel}

If the $e_i$ are known, evaluating $f$ at the powers $\alpha, \alpha^2,
\ldots, \alpha^t$ (componentwise) gives
$$
v_k = f(\alpha^k) = \sum_{i=1}^{t} c_i \mu_i^k, \qquad k = 1, \ldots, t,
$$
a linear system for the $c_i$ whose matrix $(\mu_i^k)$ is a transposed
Vandermonde matrix scaled by the $\mu_i$. It is nonsingular if and only if
the $\mu_i$ are distinct, and it can be solved in quadratic time.

::: lemma Transposed Vandermonde solution {#lem-vand}
Let $M(z) = \prod_{i=1}^{t} (z - \mu_i)$ with distinct nonzero $\mu_i$, and
let $q_i(z) = M(z)/(z - \mu_i) = \sum_{k=0}^{t-1} q_{i,k} z^k$. Then
$$
c_i = \frac{\sum_{k=0}^{t-1} q_{i,k} \, v_{k+1}}{\mu_i \, q_i(\mu_i)} .
$$
All $c_i$ can be computed with $O(t^2)$ operations from the $v_k$ and the
coefficients of $M$.
:::

::: proof
$\sum_k q_{i,k} v_{k+1} = \sum_j c_j \mu_j \sum_k q_{i,k} \mu_j^k = \sum_j c_j
\mu_j q_i(\mu_j)$, and $q_i(\mu_j) = 0$ for $j \ne i$, while $q_i(\mu_i) =
M'(\mu_i) \ne 0$. Given $M$, the coefficients of $q_i$ follow by synthetic
division in $O(t)$ operations, and the numerator and denominator are
accumulated during that division.
:::

`_nmod_zip_vand_solve` implements this. With more than $t$ values it also
checks the remaining equations, which detects a wrong support: it
returns $-1$ if two $\mu_i$ coincide (choose another point) and $0$ if
the system is inconsistent (the assumed support was wrong).

### Zippel's variable-by-variable method {#zippel-shape}

Zippel's algorithm applies this inside a recursion on the variables.
Suppose $f$ has been interpolated in $x_1, \ldots, x_{j}$ with the
remaining variables at random values; write it as $\sum_\ell
f_\ell(x_1, \ldots, x_j) \, x_{j+1}^{\ell}$, the *shape*. To add $x_{j+1}$,
Zippel assumes that each coefficient of $f$ as a polynomial in $x_{j+1},
\ldots$ has the same monomials in $x_1, \ldots, x_j$ as in the image; new
images at other values of $x_{j+1}$ are then found by
[](#lem-vand) instead of a full recursive computation, and the
dependence on $x_{j+1}$ is recovered by dense univariate interpolation.

::: lemma Probability of a wrong shape {#lem-shape}
Let $f \in \F_p[x_1, \ldots, x_m]$ and let $f = \sum_{u} f_u(x_{j+1},
\ldots, x_m) \, x^{u}$ over monomials $x^u$ in $x_1, \ldots, x_j$. If the
values $\beta$ of $x_{j+1}, \ldots, x_m$ are chosen uniformly at random from
$S \subseteq \F_p$, the image $f(x_1, \ldots, x_j, \beta)$ has fewer
monomials than $f$ has distinct $x^u$ with probability at most $t \cdot
\deg f / |S|$.
:::

::: proof
The image lacks $x^u$ exactly when $f_u(\beta) = 0$; each $f_u$ is a
nonzero polynomial of degree at most $\deg f$, so this has probability at
most $\deg f / |S|$ by the Schwartz–Zippel lemma
([](algorithms/multivariate/modp-gcd#lem-sz)), and there are at most $t$
of them.
:::

The method needs $O(t d)$ values per variable for degree $d$, against
$d^m$ for dense interpolation. A wrong shape is detected later, because
the linear systems become inconsistent or the final verification (a
division) fails, and the computation is restarted with new points.

### Scaled images {#scaling}

In GCD computations the values are not values of $f$ itself: an image
$\gcd(A(x_0, \alpha^k), B(x_0, \alpha^k))$ is only determined up to a
nonzero scalar $s_k$ unless the leading coefficient of the GCD is known.
If the GCD is $\sum_\ell x_0^\ell \sum_i c_{\ell i} x^{e_{\ell i}}$ and the
images are $\sum_\ell x_0^\ell g_{\ell k}$, then
$$
\sum_i c_{\ell i} \mu_{\ell i}^{k} = s_k \, g_{\ell k} \quad \text{for all } \ell, k,
$$
a linear system in the $c_{\ell i}$ and the $s_k$, homogeneous, and
normalised by fixing one coefficient. De Kleine, Monagan and Wittkopf
[@dKMW2005] solve it by eliminating the $c_{\ell i}$ with the structure of
[](#lem-vand) for each $\ell$, which leaves a small system for the
scalars. `nmod_mpolyl_gcds_zippel` does this (with a shortcut when the
leading coefficient is a monomial, which makes the images monic).

## Unknown support: Ben-Or–Tiwari {#bma}

::: theorem Ben-Or–Tiwari {#thm-bot}
With $v_k = f(\alpha^k) = \sum_{i=1}^{t} c_i \mu_i^k$ for $k \ge 0$ and
distinct $\mu_i$, the sequence $(v_k)$ satisfies the linear recurrence
with characteristic polynomial $\Lambda(z) = \prod_{i=1}^t (z - \mu_i) =
z^t + \lambda_{t-1} z^{t-1} + \cdots + \lambda_0$, that is,
$v_{k+t} + \sum_{j<t} \lambda_j v_{k+j} = 0$ for all $k \ge 0$, and no
recurrence of smaller order. $\Lambda$ is determined by $v_0, \ldots,
v_{2t-1}$ and can be computed from them by the Berlekamp–Massey algorithm
in $O(t^2)$ operations.
:::

::: proof
$v_{k+t} + \sum_j \lambda_j v_{k+j} = \sum_i c_i \mu_i^k \Lambda(\mu_i) =
0$. A recurrence of order $s < t$ would give $\sum_i c_i \mu_i^k P(\mu_i) = 0$
for a nonzero $P$ of degree $s$ and all $k$, that is, a vanishing
combination of the geometric sequences $(\mu_i^k)$ with coefficients
$c_i P(\mu_i)$, not all zero since $P$ has fewer than $t$ roots; this
contradicts the independence of geometric sequences with distinct ratios
(a Vandermonde determinant). The minimal polynomial of a linearly
recurrent sequence of order $t$ is determined by its first $2t$ terms, and
the Berlekamp–Massey algorithm [@Mas1969; @ADL2006] finds it with $O(t^2)$
operations.
:::

The roots of $\Lambda$ are the monomial values $\mu_i$, from which the
exponents must be recovered, and then the coefficients follow by
[](#lem-vand). Two devices make the recovery cheap:

- *Kronecker substitution*: with degree bounds $d_j > \deg_{x_j} f$, the
  substitution $x_j \mapsto z^{D_j}$, $D_j = \prod_{l > j} d_l$, maps the
  monomials of $f$ injectively to powers of one variable $z$, with
  exponents $E_i < D = \prod_j d_j$ ([](algorithms/multivariate/multiplication#lem-kronecker)).
- *Discrete logarithms*: choose a prime $p > D$ for which $p - 1$ has only
  small prime factors, and a generator $\omega$ of $\F_p^*$. Evaluating at
  $z = \omega^k$ gives $\mu_i = \omega^{E_i}$, the $\mu_i$ are distinct
  because $E_i < D < p - 1$, and $E_i$ is the discrete logarithm of
  $\mu_i$, computed by the Pohlig–Hellman algorithm [@PohHel1978] in time
  polynomial in the largest prime factor of $p - 1$.

In practice the number of terms $t$ is unknown: the Berlekamp–Massey
algorithm is run incrementally as values arrive, and the candidate
$\Lambda$ is accepted when it stops changing, its roots are all in
$\F_p$ and distinct, the discrete logarithms are below $D$, and the
resulting polynomial passes a check at a random point. Compared with
Zippel's method this needs about $2t$ values in total instead of $O(t)$
per variable, but the prime must be larger than $D$, which grows with the
degrees and the number of variables.

`fmpz_mpolyl_gcd_zippel2` uses this for the coefficients of a GCD with
respect to its two main variables ([](algorithms/multivariate/gcd#zippel2)):
`fmpz_next_smooth_prime` chooses the prime, `nmod_berlekamp_massey_t`
(or `fmpz_mod_berlekamp_massey_t` for primes beyond a word) runs the
recurrence, the roots are found with `nmod_poly` root finding, and
`nmod_discrete_log_pohlig_hellman_run` gives the exponents. Images modulo
further primes reuse the support and only solve transposed Vandermonde
systems.

## Berlekamp–Massey as a partial Euclidean algorithm {#bm-hgcd}

FLINT computes the recurrence of [](#thm-bot) by rational reconstruction
with the extended Euclidean algorithm, in a form that accepts the values
in batches and reuses the previous work, with the half-gcd for large
steps [@ThullYap1990]. For values $a_1, a_2, \ldots \in K$ ($K$ a field)
write
$$
S_n = a_1 z^{n-1} + a_2 z^{n-2} + \cdots + a_n \in K[z],
\qquad
\sum_{i \ge 1} a_i z^{-i} \in K[[z^{-1}]],
$$
so that $S_n / z^n$ is the truncation of the series. If $a_i = f(\alpha^{s
+ i - 1})$ with $f$, $\alpha$, $\mu_k$ as above and a shift $s \in \Z$,
then
$$
\sum_{i \ge 1} a_i z^{-i} = \sum_{k=1}^{t} \frac{c_k \mu_k^{s}}{z - \mu_k},
$$
a rational function whose reduced denominator is $\Lambda = \prod_k (z -
\mu_k)$ when the $\mu_k$ are distinct, since every numerator $c_k
\mu_k^s$ is nonzero. Finding $\Lambda$ is therefore a rational
reconstruction problem.

::: theorem Rational reconstruction {#thm-ratrec}
Suppose that $\sum_{i \ge 1} a_i z^{-i} = -\bar u / \bar v$ with $\deg \bar
u < \deg \bar v \le n/2$, and that $u, v, r \in K[z]$ satisfy
$$
u z^n + v S_n = r, \qquad \deg u < \deg v \le n/2, \qquad
\deg r < n/2, \qquad \deg r < \deg v .
$$
Then $u / v = \bar u / \bar v$.
:::

::: proof
Dividing by $v z^n$ gives $u/v + S_n/z^n = r/(v z^n)$. Here $S_n / z^n =
-\bar u/\bar v + O(z^{-n-1})$, and $r/(v z^n) = O(z^{-n-1})$ because $\deg
r < \deg v$. Hence $u/v - \bar u/\bar v = (u \bar v - \bar u v)/(v \bar v)$
is $O(z^{-n-1})$. A nonzero numerator would make it at least of order
$z^{-\deg v - \deg \bar v}$, and $\deg v + \deg \bar v \le n$, so the
numerator is zero.
:::

The Euclidean remainder sequence $r_0 = z^n$, $r_1 = S_n$, $r_2, \ldots$
of [](algorithms/polynomials/half-gcd) comes with Bézout coefficients
$u_j z^n + v_j S_n = r_j$, $\deg v_j = n - \deg r_{j-1}$. At the index
$j$ of the half-gcd, $2 \deg r_{j-1} \ge n > 2 \deg r_j$, so $\deg v_j \le
n/2$ and $\deg r_j < n/2$: if also $\deg r_j < \deg v_j$, the theorem
applies to $(u_j, v_j, r_j)$. Conversely, if the series has a denominator
of degree $t \le n/2$, then $\bar u z^n + \bar v S_n$ has degree less than
$t$, and by the uniqueness of rational reconstruction [@GG2003, §5.7]
$(\bar u, \bar v)$ is a multiple of the row $(u_j, v_j)$, a scalar
multiple because $\bar u$ and $\bar v$ are coprime. So $2t$ values
determine $\Lambda$, as the half-gcd row $v_j$ made monic.

Since $t$ is not known, the computation is repeated as values are added,
and $v_j$ is accepted once further values no longer change it. The next
theorem shows how the half-gcd for $n$ values is extended to $n + m$
values. We use the notation of [](algorithms/polynomials/half-gcd): $(a,
b)^T = M (A, B)^T$ with $M$ a product of quotient matrices, and $\deg a$
the degree of the first entry.

::: theorem Extending the half-gcd {#thm-bm-extend}
Let $M$ be the matrix of the half-gcd of $(z^n, S_n)$, with remainders
$(R_0, R_1)$, $2 \deg R_0 \ge n > 2 \deg R_1$, and let $S_{n+m} = z^m S_n +
s$ with $\deg s < m$. Define
$$
\begin{pmatrix} A \\ B \end{pmatrix} = M^{-1} \begin{pmatrix} z^{n+m} \\
S_{n+m} \end{pmatrix} = z^m \begin{pmatrix} R_0 \\ R_1 \end{pmatrix} +
M^{-1} \begin{pmatrix} 0 \\ s \end{pmatrix}.
$$
Then $A$, $B$ are consecutive remainders of $(z^{n+m}, S_{n+m})$ with $2
\deg A \ge n + m$ and $\deg B < m + n/2$. If $2 \deg B < n + m$, then $M$
is the matrix of the half-gcd of $(z^{n+m}, S_{n+m})$. Otherwise let $C = B$,
$D = A \bmod B$, $q = A \operatorname{quo} B$, $l = \deg C$ and $k = n + m
- l$. Then $0 < k \le l$, and if $H$ is the matrix of the half-gcd of
$(C \operatorname{quo} z^k, D \operatorname{quo} z^k)$, the matrix of
the half-gcd of $(z^{n+m}, S_{n+m})$ is $M Q H$ with $Q = \bigl(
\begin{smallmatrix} q & 1 \\ 1 & 0 \end{smallmatrix} \bigr)$.
:::

::: proof
The entries of $M$, and hence of $M^{-1} = \pm \operatorname{adj} M$, have
degree at most $n - \deg R_0 \le n/2$. So the second term has degree less
than $m + n/2$, while $z^m R_0$ has degree $m + \deg R_0 \ge m + n/2$:
thus $\deg A = m + \deg R_0 \ge (n + m)/2$, and $\deg B < m + n/2$
because $z^m R_1$ also has degree less than $m + n/2$. As $\deg A > \deg
B$, [](algorithms/polynomials/half-gcd#lem-unique) shows that $A, B$
are consecutive remainders of $(z^{n+m}, S_{n+m})$ (this is the lifting
lemma [](algorithms/polynomials/half-gcd#lem-lift) with $a_0 = 0$). If $2
\deg B < n + m$ they are the half-gcd remainders. Otherwise $2l \ge n + m$
gives $k \le l$, and $l < m + n/2 \le n + m$ gives $k > 0$. The pair $(C,
D)$ consists of the next consecutive remainders. By the lifting lemma
applied to $(C, D)$ split at $z^k$, the remainders $(C', D')^T = H^{-1}(C,
D)^T$ are consecutive with $\deg C' = k + \deg C^* \ge k + (l - k)/2 = (n
+ m)/2$, where $C^*$ is the first half-gcd remainder of $(C
\operatorname{quo} z^k, D \operatorname{quo} z^k)$, whose first entry has
degree $l - k$. The bound on the degrees of the entries of $H^{-1}$ used
above gives, in the same way, $\deg D' < k + (l - k)/2 = (n + m)/2$. So
$(C', D')$ are the half-gcd remainders of $(z^{n+m}, S_{n+m})$.
:::

The algorithm keeps the remainders $R_0, R_1$ and the second Bézout
coefficients $V_0, V_1$ (the first ones are not needed), with $U_i z^n +
V_i S_n = R_i$; multiplying $(R_0, R_1)$ and $(V_0, V_1)$ by the same
matrix on the left preserves these relations.

::: algorithm Incremental Berlekamp–Massey {#alg-bm-incr}
State: $n$; $R_0, R_1, V_0, V_1$ as above, initially $n = 0$, $R_0 = 1$,
$R_1 = 0$, $V_0 = 0$, $V_1 = 1$.
Input: new values $a_{n+1}, \ldots, a_{n+m}$.
Output: the updated state; $V_1$ is the candidate denominator.

1. Let $s = \sum_{i=1}^{m} a_{n+i} z^{m-i}$. Set $R_i \gets z^m R_i + V_i s$
   for $i = 0, 1$ and $n \gets n + m$.
2. If $2 \deg R_1 < n$, return: $V_1$ has not changed.
3. Do one Euclidean step: $q = R_0 \operatorname{quo} R_1$, $(R_0, R_1)
   \gets (R_1, R_0 - q R_1)$, $(V_0, V_1) \gets (V_1, V_0 - q V_1)$.
4. Let $l = \deg R_0$ and $k = n - l$. If $l - k$ is small, do Euclidean
   steps until $2 \deg R_1 < n$. Otherwise compute the half-gcd matrix
   $H$ of $(R_0 \operatorname{quo} z^k, R_1 \operatorname{quo} z^k)$ and
   replace $(R_0, R_1)^T$ and $(V_0, V_1)^T$ by $H^{-1} (R_0, R_1)^T$ and
   $H^{-1}(V_0, V_1)^T$.
:::

By [](#thm-bm-extend) the state after step 4 is the half-gcd of $(z^n,
S_n)$ with its Bézout coefficients, and by [](#thm-ratrec) $V_1$ is the
denominator of the series as soon as $n \ge 2t$. When the new values do
not change the reconstruction, only step 1 is executed, at the cost of two
products of polynomials of degrees at most $n/2$ and less than $m$. When
they do, the half-gcd works on polynomials of degree $l - k = 2l - n$,
about twice the number of Euclidean steps still needed, rather than on
polynomials of degree $n$.

`nmod_berlekamp_massey_reduce` implements [](#alg-bm-incr), using
Euclidean steps when $l - k < 10$ and `nmod_poly_hgcd` otherwise; it
returns whether $V_1$ changed. `fmpz_mod_berlekamp_massey_reduce` has the
same structure but always uses Euclidean steps. In
`fmpz_mpolyl_gcd_zippel2` each coefficient of the bivariate GCD has its
own instance; the values $a_i$ are the images at the successive powers
$\omega^{s}, \omega^{s+1}, \ldots$ of the evaluation point, the reduction
is run after every eighth value, and a reduction that leaves a set of
candidates unchanged triggers the root finding of [](#bma) and the
verification. The coefficients $c_k$ then follow from the roots and the
first $t$ values by [](#lem-vand), with the shift $\mu_k^{s}$ divided out.

## Sparse interpolation in Hensel lifting {#hensel}

In sparse Hensel lifting [@Zip1981; @Kal1985; @MonTun2020] the unknowns
are the coefficients of the lifted factors as polynomials in the new
variable $x_j$: each coefficient of $x_j^\ell$ is assumed to have the same
support in the previous variables as the factor itself, and its values are
obtained from bivariate factorisations at the points $(\beta_2^k, \ldots,
\beta_{j-1}^k)$ for $k = 1, 2, \ldots$. The values of all factors at the same
points are only determined up to scalars whose product is known, and the
leading coefficients are imposed as in the GCD case.
`nmod_mpoly_hlift_zippel` (over $\F_p$) and the $p$-adic lifting inside
`fmpz_mpoly_factor_irred_zippel` work this way; see
[](algorithms/multivariate/hensel-lifting#sparse).

## Use in FLINT {#use}

- Transposed Vandermonde systems ([](#lem-vand)): `_nmod_zip_vand_solve`,
  `_n_fq_zip_vand_solve`, `_n_fqp_zip_vand_solve`, with the evaluation
  helpers `_nmod_zip_eval_step` and `n_polyun_zip_start`,
  `n_polyun_zip_solve`.
- Zippel's GCD with scaled images ([](#scaling)): `nmod_mpolyl_gcds_zippel`,
  `nmod_mpolyl_gcdp_zippel_smprime`, used by `fmpz_mpolyl_gcd_zippel` and
  `nmod_mpoly_gcd_zippel`.
- Ben-Or–Tiwari ([](#thm-bot)): `fmpz_mpolyl_gcd_zippel2`, with
  `nmod_berlekamp_massey_t`, `fmpz_mod_berlekamp_massey_t`,
  `nmod_discrete_log_pohlig_hellman_t` and
  `fmpz_mod_discrete_log_pohlig_hellman_t`.
- Incremental Berlekamp–Massey ([](#alg-bm-incr)):
  `nmod_berlekamp_massey_reduce` (with `nmod_poly_hgcd`),
  `fmpz_mod_berlekamp_massey_reduce`.
- Sparse Hensel lifting: `nmod_mpoly_hlift_zippel`,
  `fmpz_mpoly_factor_irred_zippel`, `nmod_mpoly_factor_irred_smprime_zippel`.

::: functions fmpz_mpolyl_gcd_zippel fmpz_mpolyl_gcd_zippel2 nmod_mpolyl_gcds_zippel nmod_mpoly_hlift_zippel _nmod_zip_vand_solve :::
