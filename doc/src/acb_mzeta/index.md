---
title: acb_mzeta.h
header: acb_mzeta.h
description: Multivariate zeta functions; multiple zeta values, multiple polylogarithms, Epstein zeta functions of lattices in any dimension.
category: Special functions
order: 6
---

This module collects zeta functions in several variables or attached to
higher-dimensional objects, including the classical multiple zeta values,
multiple polylogarithms, and the Epstein zeta function of a lattice.
The algorithms, with their error bounds, are described in
[](algorithms/special-functions/multiple-zeta) and
[](algorithms/special-functions/epstein).

::: summary algorithms/special-functions/multiple-zeta :::

::: summary algorithms/special-functions/epstein :::

## Multiple zeta values

::: function
void acb_mzeta_mzv(arb_t res, const ulong * s, slong r, slong prec)

Sets `res` to the multiple zeta value

$$
\zeta(s_1, \ldots, s_r) = \sum_{n_1 > n_2 > \cdots > n_r \ge 1} \frac{1}{n_1^{s_1} \cdots n_r^{s_r}},
$$

with the convention of PARI/GP's `zetamult` (the first index belongs to
the largest summation variable; $s_1 \ge 2$ for convergence, otherwise
the result is indeterminate). The empty index gives 1, depth 1 uses
`arb_zeta_ui`.

For depth $r \ge 2$ the value is computed from the iterated integral
representation by Akhilesh's double tails [@Akh2017], a recursion over
the subwords of the index in fixed-point arithmetic, with cost
$O(k^2 P^2)$ for weight $k$ and $P$ bits; from 25000 bits it is instead
reduced by Hölder convolution [@BBBL2001] to multiple polylogarithms at
$\tfrac12$, evaluated by binary splitting in quasi-linear time. See
[](algorithms/special-functions/multiple-zeta#double-tails) and
[](algorithms/special-functions/multiple-zeta#holder).
:::

::: function
void acb_mzeta_mzv_star(arb_t res, const ulong * s, slong r, slong prec)
void acb_mzeta_mzv_t(acb_t res, const ulong * s, slong r, const acb_t t, slong prec)

The multiple zeta star value
$\zeta^\star(s_1, \ldots, s_r) = \sum_{n_1 \ge \cdots \ge n_r \ge 1} n_1^{-s_1} \cdots n_r^{-s_r}$,
respectively Yamamoto's interpolation [@Yam2013]

$$
\zeta^t(s) = \sum_{p} t^{\,r - \operatorname{depth}(p)} \zeta(p),
$$

the sum over the $2^{r-1}$ indices $p$ obtained from $s$ by replacing
some of the commas by plus signs ($t = 0$ gives $\zeta(s)$, $t = 1$ gives
$\zeta^\star(s)$), as with PARI/GP's `zetamult(s, t)`. The values
$\zeta(p)$ are computed together, sharing their subwords, and the
result is a polynomial in $t$ evaluated at `t`; for $t$ far from
$[0, 1]$ the absolute error is relative to the largest term. Depth
at most 24.
:::

::: macro
ACB_MZETA_MZV_STAR
ACB_MZETA_MZV_WEIGHT

Flags for `acb_mzeta_mzv_all`: star values instead of multiple zeta
values; only the weight $k$ instead of all weights up to $k$. The
values are those of PARI/GP's `zetamultall` with the bits 0 and 2 of its
flag.
:::

::: function
void acb_mzeta_mzv_all(arb_ptr res, ulong k, int flags, slong prec)

Sets `res` to all multiple zeta values of weight $2, \ldots, k$
($2^{k-1} - 1$ values) or of weight $k$ only ($2^{k-2}$ values), ordered
by their index $m$ (see `acb_mzeta_mzv_index`), i.e. starting at
$m = 1$, respectively $m = 2^{k-2}$. The words are evaluated together
by the double tail recursion, so that shared subwords are computed once.
Requires $2 \le k \le$ `FLINT_BITS` $- 2$.
:::

::: function
ulong acb_mzeta_mzv_index(const ulong * s, slong r)
slong acb_mzeta_mzv_index_composition(ulong * s, ulong m)

Conversions between an admissible index $s$ and the integer
$m = (1 y)_2$ where $x_0 y x_1$ is the word of $s$ (written with the
digits $0 = x_0$, $1 = x_1$), as in PARI/GP's `zetamultconvert`: the
indices of weight $k$ are $2^{k-2} \le m < 2^{k-1}$. The first function
returns 0 if $s$ is not admissible or its weight exceeds `FLINT_BITS`;
the second writes the index to `s` (room for $k - 1$ entries) and
returns its depth.
:::

## Multiple polylogarithms and colored multiple zeta values

::: function
void acb_mzeta_polylog(acb_t res, const ulong * s, acb_srcptr x, slong r, slong prec)

Sets `res` to the multiple polylogarithm

$$
\operatorname{Li}_{s_1, \ldots, s_r}(x_1, \ldots, x_r) = \sum_{n_1 > \cdots > n_r \ge 1} \frac{x_1^{n_1} \cdots x_r^{n_r}}{n_1^{s_1} \cdots n_r^{s_r}}
$$

(the convention of PARI/GP's `polylogmult`). The value is computed from
Goncharov's iterated integrals $G(a_1, \ldots, a_L; 1)$ with letters
$a_i \in \{0, 1/y_1, \ldots, 1/y_r\}$, $y_i = x_1 \cdots x_i$, as power
series at 0 by binary splitting, either directly (when all nonzero
letters have absolute value at least 2) or after splitting the path at
$\tfrac12$ (when the nonzero letters $a_i$ and $1 - a_i$ have absolute
value at least $\tfrac34$); see
[](algorithms/special-functions/multiple-zeta#polylog). Other arguments,
for example $x$ outside the unit disk or letters close to $\tfrac12$,
are not supported yet, and the result is then indeterminate, as it is
for divergent sums ($s_1 = 1$, $x_1 = 1$).
:::

::: function
void acb_mzeta_mzv_colored(acb_t res, const ulong * s, const ulong * e, slong r, ulong N, slong prec)

The colored multiple zeta value
$\operatorname{Li}_{s_1, \ldots, s_r}(\zeta^{e_1}, \ldots, \zeta^{e_r})$,
$\zeta = e^{2 \pi i/N}$; for $N = 2$ these are the alternating sums (Euler
sums). Supported for $N \le 8$, where the algorithm of
`acb_mzeta_polylog` applies to all indices; the result is indeterminate
for larger $N$.
:::

## Cache

::: type
acb_mzeta_cache_struct
acb_mzeta_cache_t

An opaque store of data that algorithms precompute and may reuse between
calls with the same lattice or parameters: the reduced lattice
(LLL-reduced basis, Cholesky factors, tail-bound constants), and the
tables of the double-precision algorithm (one per value of $\nu$). A
cache may be shared by different algorithms; the most recent entries of
each kind are kept. A cache is not thread-safe: use one per thread.
:::

::: function
void acb_mzeta_cache_init(acb_mzeta_cache_t cache)
void acb_mzeta_cache_clear(acb_mzeta_cache_t cache)

Initializes or frees a cache.
:::

## Epstein zeta function

The Epstein zeta function of a lattice in $\R^d$ with shifts in the
lattice and in its dual is defined, following the conventions of
Buchheit, Busse and Gutendorf [@BBG2024], by

$$
Z_{\Lambda,\nu}\!\begin{vmatrix} x \\ y \end{vmatrix}
= \sum_{z \in \Lambda}{}' \frac{e^{-2 \pi i\, y \cdot z}}{|x - z|^{\nu}},
\qquad \Lambda = A \Z^d,
$$

where the columns of $A$ are a basis of the lattice and the prime omits
$z = x$. The sum converges for $\operatorname{Re}(\nu) > d$ and continues
to a meromorphic function of $\nu$ with at most a simple pole at $\nu = d$,
which occurs when $y$ is in the dual lattice $\Lambda^* = A^{-T} \Z^d$.
Lattice sums of physics are special cases: the Madelung constant of NaCl
is $-Z_{\Z^3,1}(0; (\tfrac12, \tfrac12, \tfrac12))$, and $Z_{\Z,\nu}(0; 0) = 2 \zeta(\nu)$.

All the variants below share one interface. The lattice is given by a
$d \times d$ basis matrix `A` (columns), the shifts `x` and `y` by vectors
of length $d$, and $\nu$ is any complex ball. With the flag
`ACB_MZETA_EPSTEIN_REGULARIZED` they compute instead the regularised
function

$$
Z^{\mathrm{reg}}_{\Lambda,\nu}\!\begin{vmatrix} x \\ y \end{vmatrix}
= e^{2 \pi i\, x \cdot y} Z_{\Lambda,\nu}\!\begin{vmatrix} x \\ y \end{vmatrix}
- \frac{\hat s_\nu(y)}{|\det A|},
\qquad
\hat s_\nu(y) = \frac{\pi^{\nu/2} \Gamma((d - \nu)/2)}{\Gamma(\nu/2)} (\pi y^2)^{(\nu - d)/2},
$$

which is holomorphic in $y$ near $0$ (and equals $Z$ at $y = 0$). The
case $\nu \in d + 2 \Z_{\ge 0}$, where $\hat s_\nu$ needs a logarithmic
definition, is not supported (the output is then indeterminate).

The functions return 1 when `res` has been computed, and 0 (with `res`
indeterminate) when the algorithm does not apply to the input. A
computed result may be indeterminate, for example at the pole, or when
the ball `x` (resp. `y`) contains a point of $\Lambda$ (resp.
$\Lambda^*$) without being exactly that point. The `cache` may be `NULL`.

::: macro
ACB_MZETA_EPSTEIN_REGULARIZED

Flag selecting the regularised function.
:::

::: function
int acb_mzeta_epstein(acb_t res, const acb_t nu, const arb_mat_t A, arb_srcptr x, arb_srcptr y, int flags, slong prec)
int acb_mzeta_epstein_generic(acb_t res, const acb_t nu, const arb_mat_t A, arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)

Sets `res` to the Epstein zeta function, choosing an algorithm:
`acb_mzeta_epstein_shell` when the input is exact and rational, above 95
bits; `acb_mzeta_epstein_crandall_dfloat` up to 95 bits for real $\nu$,
unless its result has lost more than 12 bits of accuracy;
`acb_mzeta_epstein_theta` in dimension 3 from 800 bits, 4 from 400 bits
and 5 to 8 from 250 bits; and otherwise
`acb_mzeta_epstein_crandall`. Always returns 1.

The first form keeps a thread-local cache with the data of the most
recently used lattices and $\nu$ (the typical use being many evaluations
with a fixed lattice), released by `flint_cleanup`; the second form,
like the individual algorithms below, takes an explicit cache (or
`NULL` for none).
:::

::: function
int acb_mzeta_epstein_crandall(acb_t res, const acb_t nu, const arb_mat_t A, arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)

Evaluates Crandall's representation [@BBG2024, Thm. 2.11], which writes
$Z$ as two exponentially convergent sums of incomplete gamma functions
over $\Lambda$ and $\Lambda^*$, term by term in ball arithmetic, for any
$d$, $\nu$ and precision. The sums are truncated rigorously with the
tail bound of theta series; their number of terms grows like $P^{d/2}$
for $P$ bits. With an exact Gram matrix, the incomplete gamma function
is evaluated once per value of the quadratic form. See
[](algorithms/special-functions/epstein#crandall) and
[](algorithms/special-functions/epstein#arb).
:::

::: function
int acb_mzeta_epstein_crandall_dfloat(acb_t res, const acb_t nu, const arb_mat_t A, arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)

The same sums in double-precision ball arithmetic ([`dfloat`](dfloat)
balls `d1b` for `prec` at most 45, `d2b` for `prec` at most 95), with the
lattice points processed in batches by vector operations and the
incomplete gamma function taken from rigorous tables kept in the cache
(one per value of $\nu$); see
[](algorithms/special-functions/epstein#dfloat). Applies (returns 1)
when $\nu$ is real and exact as a double, $|\nu| \le 64$, $d \le 8$,
`prec` is at most 95 and [](dfloat_is_supported) holds.
:::

::: function
int acb_mzeta_epstein_theta(acb_t res, const acb_t nu, const arb_mat_t A, arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)

Writes each of the two Crandall sums as a Mellin integral
$\int_1^\infty w^{a-1} \Theta(w)\, dw$ of a Riemann theta function and
computes it by Gauss–Legendre quadrature with rigorous error bounds,
using `acb_theta_one` for small $w$ and direct sums over a few lattice
points for large $w$; see [](algorithms/special-functions/epstein#theta).
The cost is about $P$ integrand values per sum at $P$ bits instead of
$O(P^{d/2})$ terms, so this method is preferable to
`acb_mzeta_epstein_crandall` at high precision. Applies for $d \le 8$.
The computation is parallel with up to `flint_get_num_threads()`
threads, and the result does not depend on the number of threads.
:::

::: function
int acb_mzeta_epstein_shell(acb_t res, const acb_t nu, const arb_mat_t A, arb_srcptr x, arb_srcptr y, int flags, acb_mzeta_cache_t cache, slong prec)
int acb_mzeta_epstein_gram_shell(acb_t res, const acb_t nu, const fmpq_mat_t G, const fmpq * u, const fmpq * v, int flags, slong prec)

The algorithm by shells, for exact rational input: the first form
requires $A$, $x$, $y$ to be exact (dyadic, so that the Gram matrix
$A^T A$ and the coordinates are rational); the second takes the Gram
matrix $G$ of the lattice and the coordinates $u$, $v$ with $x = A u$,
$y = A^{-T} v$ (then $|x - z|^2 = (u - n)^T G (u - n)$ and $y \cdot z = v \cdot n$
for $z = A n$), so that lattices such as the hexagonal one can be given
exactly. $\nu$ must be real and exact with a short dyadic value. Returns
0 (not applicable) for the regularised function, for $d > 8$, for
denominators of $u$, $v$ above 64 and when the estimated cost of the
shell counts is too large.

The lattice points are grouped into shells of equal norm, counted
exactly (by products of theta series of low-rank blocks, by enumeration,
or from sums of three squares), so that each side of the Crandall
splitting needs one incomplete gamma value per shell; the splitting
parameter is chosen so that one side has exact rational arguments
(binary splitting) and, for integer indices such as the Coulomb case
$\nu = d - 2$, the other side is a power series with rational
coefficients. See [](algorithms/special-functions/epstein#shells). The
computation is parallel with up to `flint_get_num_threads()` threads,
and the result does not depend on the number of threads.
:::

## Madelung constant

::: function
void acb_mzeta_madelung_nacl(arb_t res, slong prec)

Sets `res` to the Madelung constant of sodium chloride,

$$
M = -\sum_{(i,j,k) \ne 0} \frac{(-1)^{i+j+k}}{\sqrt{i^2+j^2+k^2}} = 1.7475645946331821906\ldots,
$$

computed as $-Z_{\Z^3,1}\bigl(0; (\tfrac12,\tfrac12,\tfrac12)\bigr)$ by
`acb_mzeta_epstein_gram_shell`, which here reduces to Crandall's
expansion [@Cra1999, eq. (1.2)] as an alternating sum of $r_3(N)$ times
complementary error functions and a power series in
$e^{-2\pi^2/\lambda}$ with coefficients $r_3(8k+3)/(8k+3)$; see
[](algorithms/special-functions/epstein#madelung). The computation is
parallel with up to `flint_get_num_threads()` threads, and the result
does not depend on the number of threads.
:::

## Internal functions

::: function
void _acb_mzeta_gamma_upper_scaled(acb_t res, const acb_t a, const arb_t t, const acb_t ga, slong prec)

Sets `res` to $t^{-a} \Gamma(a, t) = \int_1^\infty u^{a-1} e^{-tu}\, du$
for real $t > 0$ to a relative accuracy of about `prec` bits; `ga` is
$\Gamma(a)$ (at a sufficient precision) or `NULL`.

The method (asymptotic series, convergent series with the working
precision raised to cover its cancellation, or an expression in $E_1$
for nonpositive integers $a$) and the working precision are chosen here
for the lattice sums, which need many values with one $a$. A ball $t$
is evaluated at its midpoint, with the radius propagated through a
bound for the derivative; see
[](algorithms/special-functions/epstein#gamma-upper).
:::
