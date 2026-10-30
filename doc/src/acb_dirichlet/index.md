---
title: acb_dirichlet.h
header: acb_dirichlet.h
description: Dirichlet L-functions, the Riemann and Hurwitz zeta functions, character sums, Hardy Z-functions and the zeros of the Riemann zeta function.
category: Special functions
order: 6
---

This module evaluates Dirichlet L-functions

$$
L(s,\chi) = \sum_{k=1}^\infty \frac{\chi(k)}{k^s}, \qquad \operatorname{Re}(s) > 1,
$$

and their analytic continuations, where $\chi$ is a Dirichlet character
of the [`dirichlet`](dirichlet) module, together with the special cases
and relatives of $L(s, \chi)$: the Riemann zeta function (the trivial
character), the Hurwitz zeta function $\zeta(s, a)$, the Lerch
transcendent and the Stieltjes constants. It also computes Gauss and
Jacobi sums and theta series of characters, discrete Fourier transforms
on groups of characters, the Hardy Z-function, and the nontrivial zeros
of $\zeta(s)$: counting them with Turing's method, isolating and
refining them by the Riemann–Siegel formula, or many at a time at large
height by Platt's method. All results are rigorous enclosures.

`acb_zeta` and `acb_hurwitz_zeta` of the [`acb`](acb) module call
`acb_dirichlet_zeta` and `acb_dirichlet_hurwitz`. The Euler–Maclaurin
summation behind the Hurwitz zeta function is `_acb_poly_zeta_cpx_series`
of [`acb_poly`](acb_poly), and polylogarithms are `acb_polylog`. The
methods are described in [](algorithms/special-functions/zeta) and
[](algorithms/special-functions/zeta-zeros), and the tutorial
[](tutorials/numerics/zeta) computes zeros and L-function values in
complete programs.

::: summary algorithms/special-functions/zeta :::

::: summary algorithms/special-functions/zeta-zeros :::

## Conventions {#conventions}

- A character is given by a group `G` (the characters modulo $q$) and an
  element `chi` of `dirichlet_char_t`, or by its Conrey number; see
  [`dirichlet`](dirichlet) and
  [](algorithms/number-theory/dirichlet-characters). Functions that accept `G` = `chi` = `NULL`
  then use the trivial character modulo 1, that is, the Riemann zeta
  function. This is the case for `acb_dirichlet_l`, `acb_dirichlet_l_jet`,
  `acb_dirichlet_l_fmpq` and the Hardy Z-functions.
- The parity of $\chi$ is $\delta = 0$ if $\chi(-1) = 1$ and $\delta = 1$
  if $\chi(-1) = -1$, and $\overline{\chi}$ is the conjugate character.
- Power series in $s$ ("jets") are written as vectors of Taylor
  coefficients: entry $k$ of the output of `acb_dirichlet_zeta_jet` is
  $\zeta^{(k)}(s)/k!$.
- Outputs may alias inputs unless stated otherwise.

## Roots of unity {#roots}

::: type
acb_dirichlet_roots_struct
acb_dirichlet_roots_t

Precomputed data for evaluating the roots of unity $e^{2\pi i k/n}$ of a
fixed order $n$, as used for character values.
:::

::: function
void acb_dirichlet_roots_init(acb_dirichlet_roots_t t, ulong order, slong num, slong prec)
void acb_dirichlet_roots_clear(acb_dirichlet_roots_t t)

Initialise `t` for the roots of unity of order $n$ = `order`, optimised
for `num` evaluations at precision `prec`, and free it.

::: details
For very small `num`, only $e^{2\pi i/n}$ is precomputed and then raised
to powers; for small `prec` and large $n$ even this may be skipped if
evaluating each root from scratch is estimated to be faster. For large
`num`, all roots in the first quadrant are precomputed, unless this would
use too much memory, and for intermediate `num` a baby-step giant-step
table is built.
:::
:::

::: function
void acb_dirichlet_root(acb_t z, const acb_dirichlet_roots_t t, ulong n, slong prec)

Sets `z` to $e^{2\pi i n / \text{order}}$.
:::

## Character values {#chi}

::: function
void acb_dirichlet_chi(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, ulong n, slong prec)

Sets `res` to $\chi(n)$, which is 0 if $\gcd(n, q) > 1$ and a root of
unity otherwise.
:::

::: function
void acb_dirichlet_chi_vec(acb_ptr v, const dirichlet_group_t G, const dirichlet_char_t chi, slong nv, slong prec)

Sets `v` to $(\chi(0), \chi(1), \ldots, \chi(\text{nv} - 1))$, sharing
the powers of one root of unity.
:::

::: function
void acb_dirichlet_pairing(acb_t res, const dirichlet_group_t G, ulong m, ulong n, slong prec)
void acb_dirichlet_pairing_char(acb_t res, const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b, slong prec)

Set `res` to the value $\chi_q(m, n)$ of the Conrey pairing, that is, the
value at $n$ of the character with Conrey number $m$; the second form
takes the two characters as input.
:::

## Gauss and Jacobi sums {#sums}

::: function
void acb_dirichlet_gauss_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
void acb_dirichlet_gauss_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
void acb_dirichlet_gauss_sum_factor(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
void acb_dirichlet_gauss_sum_order2(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
void acb_dirichlet_gauss_sum_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

Set `res` to the Gauss sum $G_q(\chi) = \sum_{x \bmod q} \chi(x) e^{2\pi i x/q}$.
The `naive` version sums the definition; `factor` writes the sum as a
product of local Gauss sums by the Chinese remainder theorem; `order2`
assumes that $\chi$ is real and primitive and returns $i^\delta \sqrt{q}$;
`theta` assumes that $\chi$ is primitive and obtains the sum from the
functional equation of the theta series at $t = 1$. The default version
reduces an imprimitive character to the primitive one inducing it and
then chooses among the methods.

::: details
For an imprimitive $\chi$ modulo $q$ induced by $\chi_0$ of conductor
$q_0$, $G_q(\chi) = \mu(q/q_0) \chi_0(q/q_0) G_{q_0}(\chi_0)$.

The `theta` method fails when the theta series vanishes at $t = 1$.
Only four characters are known with this property, $\chi_{300}(71, \cdot)$,
$\chi_{300}(131, \cdot)$, $\chi_{600}(11, \cdot)$ and
$\chi_{600}(491, \cdot)$, and none with prime power modulus;
`acb_dirichlet_gauss_sum_theta` throws an exception for them and for
imprimitive characters.
:::
:::

::: function
void acb_dirichlet_jacobi_sum(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)
void acb_dirichlet_jacobi_sum_naive(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)
void acb_dirichlet_jacobi_sum_factor(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)
void acb_dirichlet_jacobi_sum_gauss(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)
void acb_dirichlet_jacobi_sum_ui(acb_t res, const dirichlet_group_t G, ulong a, ulong b, slong prec)

Set `res` to the Jacobi sum $J_q(\chi_1, \chi_2) = \sum_{x \bmod q} \chi_1(x) \chi_2(1-x)$.
The `naive` version sums the definition; `factor` writes it as a
product of local Jacobi sums; `gauss` assumes that $\chi_1 \chi_2$ is
primitive and uses $J_q(\chi_1, \chi_2) G_q(\chi_1\chi_2) = G_q(\chi_1) G_q(\chi_2)$.
The default version handles trivial and inverse pairs by closed formulas
(and returns 0 for even $q$) and otherwise chooses among the methods. The
`ui` version takes the Conrey numbers $a$ and $b$.
:::

## Theta series {#theta}

The theta series of a character $\chi$ modulo $q$ with parity $\delta$ is

$$
\Theta_q(\chi, t) = \sum_{n \in \Z} \chi(n) n^\delta e^{-\pi n^2 t^2 / q},
\qquad t > 0.
$$

For primitive $\chi$ it satisfies
$t^{1/2+\delta} \, \Theta_q(\chi, t) = \epsilon(\chi) \, \Theta_q(\overline{\chi}, 1/t)$,
where $\epsilon(\chi)$ is the root number, which gives the Gauss sum and
the root number from $\Theta_q(\chi, 1)$.

::: function
void acb_dirichlet_theta_arb(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t t, slong prec)
void acb_dirichlet_ui_theta_arb(acb_t res, const dirichlet_group_t G, ulong a, const arb_t t, slong prec)

Set `res` to $\Theta_q(\chi, t)$ for real $t > 0$; the `ui` version
takes the Conrey number $a$. The series converges quickly only for $t$
not too small; for $t < 1$ the caller should use the functional equation
above, which is not applied automatically (to avoid recomputing the
Gauss sum).

::: details
The number of terms is `acb_dirichlet_theta_length`. For characters of
order less than 30 the series is evaluated with
`acb_dirichlet_qseries_arb_powers_smallorder`, otherwise with
`acb_dirichlet_qseries_arb_powers_naive`, and the tail bound of
`mag_tail_kexpk2_arb` is added.
:::
:::

::: function
ulong acb_dirichlet_theta_length(ulong q, const arb_t t, slong prec)
ulong acb_dirichlet_theta_length_d(ulong q, double t, slong prec)

Return a number of terms $N$ such that the tail of the theta series of
modulus $q$ at $t$ is less than $2^{-\text{prec}}$. The first version
uses a lower bound for $t$.
:::

::: function
void acb_dirichlet_qseries_arb_powers_naive(acb_t res, const arb_t x, int parity, const ulong * a, const acb_dirichlet_roots_t z, slong len, slong prec)
void acb_dirichlet_qseries_arb_powers_smallorder(acb_t res, const arb_t x, int parity, const ulong * a, const acb_dirichlet_roots_t z, slong len, slong prec)

Set `res` to $\sum_{n<\text{len}} n^p z^{a_n} x^{n^2}$, where $p$ is
`parity` (0 or 1), $z$ is the root of unity of `z`, and the exponents
$a_n$ are given (with a special value for $\chi(n) = 0$, as produced by
`dirichlet_chi_vec_order`). The `naive` version sums the series as
written; the `smallorder` version sums it in the quotient ring by the
cyclotomic polynomial and evaluates at the root of unity only at the end,
ignoring `z` except for its order.
:::

::: function
void acb_dirichlet_qseries_arb(acb_t res, acb_srcptr a, const arb_t x, slong len, slong prec)

Sets `res` to $\sum_{k=1}^{\text{len}-1} a_k x^{k^2}$ (the entry $a_0$
is ignored), computing $x^{k^2}$ by two multiplications per term.
:::

::: function
void acb_dirichlet_arb_quadratic_powers(arb_ptr v, slong nv, const arb_t x, slong prec)

Sets `v` to $(x^{k^2})_{0 \le k < \text{nv}}$.
:::

::: function
void acb_dirichlet_vec_mellin_arb(acb_ptr res, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, const arb_t t, slong n, slong prec)

Sets `res[k]`, for $0 \le k < n$, to
$t^{k/2} \sum_{m=1}^{\text{len}-1} \chi(m) m^\delta e^{-\pi m^2 t^{2k}/q}$
(half the theta series at the geometric sequence $t^k$, scaled), with
the tail bound for the truncation after `len` terms. These are samples
for computing the Mellin transform of the theta series.
:::

## Power sums {#powsum}

::: function
void acb_dirichlet_powsum_term(acb_ptr res, arb_t log_prev, ulong * prev, const acb_t s, ulong k, int integer, int critical_line, slong len, slong prec)

Sets `res` to $k^{-(s+x)}$ as a power series in $x$ of length `len`. The
flags `integer` and `critical_line` select optimisations for integer $s$
and for $\operatorname{Re}(s) = 1/2$.

::: details
On input, `log_prev` holds $\log(\text{prev})$ for the integer `*prev`;
if `*prev` is close to $k$, $\log k$ is computed from it. When
$\log k$ is computed, `log_prev` and `*prev` are overwritten with
$\log k$ and $k$, so that consecutive calls in a power sum reuse it.
:::
:::

::: function
void acb_dirichlet_powsum_sieved(acb_ptr z, const acb_t s, ulong n, slong len, slong prec)

Sets `z` to $\sum_{k=1}^n k^{-(s+x)}$ as a power series in $x$ of length
`len`. Powers of composite $k = ij$ are computed as $i^{-s} j^{-s}$ from
a table, and the even terms are grouped as a polynomial in $2^{-(s+x)}$.
This needs about $n / \log n$ powers, $n/2$ multiplications and storage
for $n/6$ power series, so it is only faster than the naive sum for
small `len`.
:::

::: function
void acb_dirichlet_powsum_smooth(acb_ptr z, const acb_t s, ulong n, slong len, slong prec)

Sets `z` to $\sum_{k=1}^n k^{-(s+x)}$ as a power series in $x$ of length
`len`, by partial sieving that adds the multiples of 5-smooth numbers
into separate buckets. It computes about $4/15$ of the powers, more than
`acb_dirichlet_powsum_sieved`, but needs only logarithmic extra space and
is faster for large `len`, since power series multiplications are
replaced by additions.
:::

## The Riemann zeta function {#zeta}

::: function
void acb_dirichlet_zeta(acb_t res, const acb_t s, slong prec)

Sets `res` to $\zeta(s)$, choosing the algorithm from $s$ and `prec`:
exact values at integers (Bernoulli numbers for negative integers,
`arb_zeta_ui` for positive ones), the Riemann–Siegel formula for
$|\operatorname{Im}(s)| \ge 24 \, p^{3/2}$ and
$|\operatorname{Re}(s)| \le 10 + p/10$ at precision $p$, the approximate
functional equation (`acb_dirichlet_l_fmpq_afe_jet`) when $s$ is an exact
real number which is a rational number of small height (in practice a
dyadic number such as $1/2$ or $-3/4$) and the precision is high enough
for this to be faster, and
Euler–Maclaurin summation otherwise, after the reflection formula for
$\operatorname{Re}(s) < 0$. The result is indeterminate if $s$ contains
the pole $s = 1$.
:::

::: function
void acb_dirichlet_zeta_jet(acb_t res, const acb_t s, int deflate, slong len, slong prec)

Sets `res` (a vector of length `len`) to the Taylor coefficients of
$\zeta(s + x)$, or of the deflated function $\zeta(s + x) - 1/(s + x - 1)$
if `deflate` is nonzero. Euler–Maclaurin summation is used (after the
reflection formula when $\operatorname{Re}(s) \le -4$, or $\le -2.5$ at
precision below 128 bits; for larger negative $\operatorname{Re}(s)$ it is
cheaper and loses fewer bits to sum directly at $s$, with guard bits
proportional to $-\operatorname{Re}(s) \log_2 p$), the Riemann–Siegel
formula for `len` $\le 2$ high on the critical line, or, for an exact
rational number $s$ of small height (including integers) at very high
precision, the approximate functional equation
(`acb_dirichlet_l_fmpq_afe_jet`), whose cost is $\tilde O(p^{1.5})$
instead of $\tilde O(p^2)$ and which needs far less memory.
:::

::: function
void acb_dirichlet_zeta_bound(mag_t res, const acb_t s)

Sets `res` to an upper bound for $|\zeta(s)|$, computed quickly: in and
slightly around the critical strip by formula (43.3) of [@Rad1973], to the
right by the value at $\operatorname{Re}(s)$, and to the left by the
functional equation.
:::

::: function
void acb_dirichlet_zeta_deriv_bound(mag_t der1, mag_t der2, const acb_t s)

Sets `der1` and `der2` to upper bounds for $|\zeta'(s)|$ and
$|\zeta''(s)|$, meant for the critical strip. The bounds are not tight.
:::

::: function
void acb_dirichlet_eta(acb_t res, const acb_t s, slong prec)

Sets `res` to the Dirichlet eta function
$\eta(s) = \sum_{k=1}^{\infty} (-1)^{k+1} k^{-s} = (1-2^{1-s}) \zeta(s)$
(the alternating character is not a Dirichlet character). The removable
singularity at $s = 1$, where $\eta(1) = \log 2$, is handled.
:::

::: function
void acb_dirichlet_xi(acb_t res, const acb_t s, slong prec)

Sets `res` to the Riemann xi function
$\xi(s) = \frac{1}{2} s (s-1) \pi^{-s/2} \Gamma(\tfrac{s}{2}) \zeta(s)$,
which satisfies $\xi(1-s) = \xi(s)$.
:::

## The Riemann–Siegel formula {#rs}

The Riemann–Siegel formula is implemented following Arias de Reyna
[@Ari2011]. For $s = \sigma + it$ with $t > 0$,

$$
\zeta(s) = \mathcal{R}(s) + X(s) \overline{\mathcal{R}}(1-s), \qquad
X(s) = \pi^{s-1/2} \frac{\Gamma((1-s)/2)}{\Gamma(s/2)},
$$

$$
\mathcal{R}(s) = \sum_{k=1}^N \frac{1}{k^s} + (-1)^{N-1} U a^{-\sigma} \left[ \sum_{k=0}^K \frac{C_k(p)}{a^k} + RS_K \right],
$$

where $U = \exp(-i[\frac{t}{2} \log(\frac{t}{2\pi})-\frac{t}{2}-\frac{\pi}{8}])$,
$a = \sqrt{t/(2\pi)}$, $N = \lfloor a \rfloor$ and $p = 1-2(a-N)$. The
coefficients $C_k(p)$ are expressed through auxiliary coefficients
$d_j^{(k)}$ and the derivatives $F^{(j)}(p)$. The cost is $O(\sqrt{t})$
terms, against $O(t)$ for Euler–Maclaurin summation; see
[](algorithms/special-functions/zeta-zeros#riemann-siegel). Because
of artificial discontinuities, $s$ should be exact inside the
evaluation; inexact input is handled by bounding the derivative.

::: variable
int acb_dirichlet_zeta_rs_use_dfloat

If nonzero (the default), the main sum of the Riemann–Siegel formula is
computed with `dfloat_powsum_sieved` (double-word arithmetic, see
[`dfloat`](dfloat)) when that gives enough accuracy. Set it to 0 to
always use `acb_dirichlet_powsum_sieved`. This is a global variable,
which should not be changed while other threads evaluate zeta.
:::

::: function
void acb_dirichlet_zeta_rs_f_coeffs(acb_ptr c, const arb_t p, slong N, slong prec)

Sets `c` to $F^{(j)}(p)$ for $0 \le j < N$, by power series division.
This breaks down at $p = \pm 1/2$, which does not occur for exact $s$.
:::

::: function
void acb_dirichlet_zeta_rs_d_coeffs(arb_ptr d, const arb_t sigma, slong k, slong prec)

Updates `d` in place from $d_j^{(k-1)}$ to $d_j^{(k)}$ for
$0 \le j \le \lfloor 3k/2 \rfloor + 1$ (for $k = 0$, initialises it).
:::

::: function
void acb_dirichlet_zeta_rs_bound(mag_t err, const acb_t s, slong K)

Sets `err` to a bound for the error term $RS_K$, by Theorem 4.2 of
[@Ari2011].
:::

::: function
void acb_dirichlet_zeta_rs_r(acb_t res, const acb_t s, slong K, slong prec)

Sets `res` to $\mathcal{R}(s)$ for $s$ in the upper half-plane, with
exactly $K$ asymptotic terms if $K > 0$ and an automatic choice
otherwise.

::: details
The main sum is computed with `dfloat_powsum_sieved` when it can deliver
about $\text{prec} + 10 - \log_2 |t|$ bits after the binary point and
has at least 16 terms, and with
`acb_dirichlet_powsum_sieved` otherwise, or always if
`acb_dirichlet_zeta_rs_use_dfloat` is 0. The double-word power sums are
described in [](dfloat/powsum).
:::
:::

::: function
void acb_dirichlet_zeta_rs(acb_t res, const acb_t s, slong K, slong prec)

Sets `res` to $\zeta(s)$ by the Riemann–Siegel formula, with exactly $K$
asymptotic terms if $K > 0$ and an automatic choice otherwise.
:::

::: function
void acb_dirichlet_zeta_jet_rs(acb_ptr res, const acb_t s, slong len, slong prec)

Sets `res` to the first `len` Taylor coefficients of $\zeta$ at $s$ by
the Riemann–Siegel formula. Only `len` $\le 2$ is supported; the
derivative is computed with a finite difference.
:::

## The Hurwitz zeta function {#hurwitz}

::: function
void acb_dirichlet_hurwitz(acb_t res, const acb_t s, const acb_t a, slong prec)

Sets `res` to the Hurwitz zeta function
$\zeta(s, a) = \sum_{k=0}^{\infty} (k + a)^{-s}$, analytically continued.
The case $a = 1$ calls `acb_dirichlet_zeta`, and $s = 0$ and small
negative integers $s$ use Bernoulli polynomials; otherwise
Euler–Maclaurin summation is used (`_acb_poly_zeta_cpx_series`). The
result is indeterminate if $s$ contains 1.
:::

::: type
acb_dirichlet_hurwitz_precomp_struct
acb_dirichlet_hurwitz_precomp_t

Precomputed data for evaluating $\zeta(s, a)$ at many $a \in (0, 1]$ for
a fixed $s$.
:::

::: function
void acb_dirichlet_hurwitz_precomp_init(acb_dirichlet_hurwitz_precomp_t pre, const acb_t s, int deflate, slong A, slong K, slong N, slong prec)

Precomputes a grid of $N$ Taylor polynomials of length $K$ for
$\zeta(s, a + A)$, for fast evaluation of $\zeta(s, a)$ at $a \in (0,1]$
with fixed $s$. This costs $NK$ Hurwitz zeta values; each later
evaluation costs $2K$ arithmetic operations plus $A$ powers, with a
truncation error $O(1/(2AN)^K)$. With `deflate` nonzero, the deflated
function (without the pole at $s = 1$) is used.

::: details
With $A = 0$, no precomputation is done and evaluation is the same as
calling `acb_dirichlet_hurwitz`. Otherwise $A$, $K$ and $N$ must be
positive (an exception is thrown otherwise), and the error bound is finite only if
$K + \operatorname{Re}(s) > 1$. To avoid an initial bump that loses
precision and slows convergence, $AN$ should be at least about as large
as $|s|$, for instance $AN > 0.5 |s|$. The error bound is that of
[](algorithms/special-functions/zeta#hurwitz-taylor).
:::
:::

::: function
void acb_dirichlet_hurwitz_precomp_init_num(acb_dirichlet_hurwitz_precomp_t pre, const acb_t s, int deflate, double num_eval, slong prec)
void acb_dirichlet_hurwitz_precomp_choose_param(ulong * A, ulong * K, ulong * N, const acb_t s, double num_eval, slong prec)

Choose $A$, $K$ and $N$ to minimise the cost of `num_eval` evaluations at
$s$ and precision `prec`, and (for the first function) initialise `pre`
with them. If evaluating each value from scratch is estimated to be
cheaper, $A = K = N = 0$.
:::

::: function
void acb_dirichlet_hurwitz_precomp_clear(acb_dirichlet_hurwitz_precomp_t pre)

Frees the precomputed data.
:::

::: function
void acb_dirichlet_hurwitz_precomp_bound(mag_t res, const acb_t s, slong A, slong K, slong N)

Sets `res` to a bound for the truncation error (not including rounding
errors) of the evaluation of $\zeta(s, a)$ with parameters $A$, $K$,
$N$, for $0 < a \le 1$.
:::

::: function
void acb_dirichlet_hurwitz_precomp_eval(acb_t res, const acb_dirichlet_hurwitz_precomp_t pre, ulong p, ulong q, slong prec)

Sets `res` to $\zeta(s, p/q)$ using the precomputed data, for
$0 < p/q \le 1$. Throws an exception if $p > q$.
:::

## The Lerch transcendent {#lerch}

::: function
void acb_dirichlet_lerch_phi(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)
void acb_dirichlet_lerch_phi_direct(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)
void acb_dirichlet_lerch_phi_integral(acb_t res, const acb_t z, const acb_t s, const acb_t a, slong prec)

Set `res` to the Lerch transcendent
$\Phi(z,s,a) = \sum_{k=0}^{\infty} z^k (k+a)^{-s}$, analytically
continued to $|z| \ge 1$. The `direct` version sums the truncated series;
the `integral` version uses a Hankel contour integral; the default
version checks special cases that reduce to the Hurwitz zeta function or
polylogarithms and otherwise chooses between the two.

::: details
The integral representation is

$$
\Phi(z,s,a) = -\frac{\Gamma(1-s)}{2 \pi i} \int_C \frac{(-t)^{s-1} e^{-a t}}{1 - z e^{-t}} dt,
$$

with the path deformed as needed to avoid the poles and branch cuts of
the integrand.
:::
:::

## Stieltjes constants {#stieltjes}

::: function
void acb_dirichlet_stieltjes(acb_t res, const fmpz_t n, const acb_t a, slong prec)

Sets `res` to the generalised Stieltjes constant $\gamma_n(a)$ for an
integer $n \ge 0$ (an exception is thrown for $n < 0$), the coefficient in the Laurent series
$\zeta(s,a) = \frac{1}{s-1} + \sum_{n=0}^\infty \frac{(-1)^n}{n!} \gamma_n(a) (s-1)^n$.
With $a = 1$ these are the ordinary Stieltjes constants.

::: details
The function uses an integral representation that allows fast computation
for extremely large $n$ [@JB2018]. If $n$ is moderate and the precision
high enough, it instead expands the Hurwitz zeta function in a power
series at $s = 1$ and reads off the last coefficient. For $a = 1$, the
power series of the deflated zeta function
$\zeta(1+x) - 1/x = \sum_n (-1)^n \gamma_n x^n / n!$ comes from
`acb_dirichlet_zeta_jet`, which at very high precision uses the
approximate functional equation (`acb_dirichlet_l_fmpq_afe_jet`) in
time $\tilde O(p^{1.5})$. To compute all of
$\gamma_0(a), \ldots, \gamma_n(a)$, expanding $\zeta(s, a)$ once (for
instance with `acb_poly_zeta_series`) is faster than calling this function
repeatedly, unless $n$ is at least several hundred.
:::
:::

## Dirichlet L-functions {#l}

For a primitive character $\chi$ of conductor $q$ and parity $\delta$,
$L(s, \chi)$ satisfies the functional equation

$$
\left(\frac{q}{\pi}\right)^{\frac{s+\delta}2}\Gamma\left(\frac{s+\delta}2\right) L(s, \chi) = \epsilon(\chi) \left(\frac{q}{\pi}\right)^{\frac{1-s+\delta}2}\Gamma\left(\frac{1-s+\delta}2\right) L(1 - s, \overline\chi)
$$

with the root number $|\epsilon(\chi)| = 1$.

::: function
void acb_dirichlet_root_number(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
void acb_dirichlet_root_number_theta(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

Set `res` to the root number $\epsilon(\chi)$ of a primitive character.
The `theta` version uses $\epsilon(\chi) = \Theta_q(\chi, 1) / \overline{\Theta_q(\chi, 1)}$;
the default version uses the Gauss sum, $\epsilon(\chi) = G_q(\chi)/(i^\delta \sqrt q)$,
when $q$ has several prime factors, and the theta series otherwise. The
default version throws an exception if $\chi$ is not primitive.
:::

::: function
void acb_dirichlet_l(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
void _acb_dirichlet_l_generic(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

Sets `res` to $L(s,\chi)$, choosing the algorithm. If $s$ is an exact
integer, `acb_dirichlet_l_si` is used. If $s$ is an exact real number which
is a rational number of small height (see `_acb_dirichlet_get_small_fmpq`),
the approximate functional equation `acb_dirichlet_l_fmpq_afe_jet` is used
when `_acb_dirichlet_l_fmpq_use_afe` estimates it to be faster (at high
precision or for large conductor). Otherwise (and always in the
underscore version): the Riemann zeta
function for $q = 1$, the functional equation for primitive $\chi$ and
$\operatorname{Re}(s) < -1/2$ (or near $s = 0$ for even $\chi$), the Euler
product for large $\operatorname{Re}(s)$, and otherwise the decomposition
into Hurwitz zeta functions with a precomputation shared by the $\phi(q)$
values.
:::

::: function
void acb_dirichlet_l_si(acb_t res, slong s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

Sets `res` to $L(s,\chi)$ for an integer $s$, using closed forms where
they exist. An imprimitive character is reduced to the primitive
character $\chi_0$ modulo $q_0$ inducing it, using
$L(s,\chi) = L(s,\chi_0) \prod_{p \mid q} (1 - \chi_0(p) p^{-s})$, and
$q_0 = 1$ gives the Riemann zeta function. For primitive $\chi$ modulo
$q > 1$ with parity $\delta$ ($\chi(-1) = (-1)^{\delta}$):

* If $s = 1 - n \le 0$, then $L(s,\chi) = -B_{n,\chi}/n$ in terms of a
  generalized Bernoulli number, which is zero unless $n \equiv \delta \pmod 2$.
  (The functional equation and the Euler product are used instead when $n$
  is large compared to the precision.)

* If $s \ge 1$ and $s \equiv \delta \pmod 2$, then

  $$L(s,\chi) = (-1)^{1 + (s - \delta)/2} \frac{\tau(\chi)}{2 i^{\delta}} \left(\frac{2 \pi}{q}\right)^s \frac{B_{s,\bar \chi}}{s!}$$

  where $\tau(\chi)$ is the Gauss sum (the Euler product is used instead
  for $s$ large compared to the precision). The generalized Bernoulli number
  $B_{s,\bar\chi} = q^{s-1} \sum_{a=1}^{q} \bar\chi(a) B_s(a/q)$ is
  computed with exact integer arithmetic grouped by the classes of residues
  with the same character value, so that only $\operatorname{ord}(\chi)$
  roots of unity are needed; at high precision the cost is dominated by the
  Gauss sum, i.e. by about $q$ multiplications.

* If $s = 1$ and $\chi$ is even, then
  $L(1,\chi) = -(\tau(\chi)/q) \sum_{a=1}^{q-1} \bar\chi(a) \log|2 \sin(\pi a/q)|$,
  computed with one logarithm per class of residues with the same character
  value (a single logarithm for real $\chi$).

* Otherwise ($s \ge 2$ and $s \not\equiv \delta \pmod 2$, where no
  closed form is known), the algorithm is chosen among the Euler product (for
  $s$ large compared to the precision), Catalan's constant
  (`arb_const_catalan`, for $q = 4$, $s = 2$),
  `acb_dirichlet_l_ui_borwein_bsplit` (at high precision, for moderate $q$),
  the approximate functional equation (`acb_dirichlet_l_fmpq_afe`, for
  larger $q$) and the Hurwitz zeta function.
:::

::: function
void acb_dirichlet_l_ui_borwein_bsplit(acb_t res, ulong s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

Sets `res` to $L(s,\chi)$ for an integer $s \ge 1$ and a primitive
character $\chi$ modulo $q$ (with `G` = `NULL` for $q = 1$), using a
generalization of Borwein's algorithm for the Riemann zeta function
(`arb_zeta_ui_borwein_bsplit`) with binary splitting, which costs
$\tilde O(\phi(q) p)$ at precision $p$.

If $4 \mid q$, the 2-component of $\chi$ is primitive, which implies that
$\chi(n + Q) = -\chi(n)$ with $Q = q/2$; if $q$ is odd, we instead use
$\sum_n (-1)^{n-1} \chi(n) n^{-s} = (1 - \chi(2) 2^{1-s}) L(s,\chi)$ with
$Q = q$. In both cases, $L(s,\chi)$ is a linear combination over the residues
$a$ modulo $Q$ of alternating series
$\sum_{j \ge 0} (-1)^j (Q j + a)^{-s}$, whose terms are moments of a positive
measure on $[0,1]$ with total mass $a^{-s}$. Borwein's acceleration with $n$ terms
therefore has error at most $3 a^{-s} / (3+\sqrt{8})^n$ for each series. The
series are evaluated together by one binary splitting over vectors
(several in parallel when threads are available), and the character values only
enter in the final linear combination.

The output is indeterminate if the character is not primitive, if
$q \equiv 2 \pmod 4$ (there are then no primitive characters), and if
$s = 1$, $q$ is odd and $\chi(2) = 1$ (when the factor $1 - \chi(2) 2^{1-s}$
vanishes).
:::

::: function
void acb_dirichlet_l_hurwitz(acb_t res, const acb_t s, const acb_dirichlet_hurwitz_precomp_t precomp, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

Sets `res` to $L(s,\chi) = q^{-s}\sum_{k=1}^q \chi(k) \zeta(s, k/q)$. For
$s = 1$ and non-principal $\chi$, the deflated Hurwitz zeta function is
used to avoid the poles. If `precomp` is `NULL`, each Hurwitz zeta value
is computed directly; otherwise the precomputation is used.
:::

::: function
void acb_dirichlet_l_euler_product(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s, const signed char * chi, int mod, int reciprocal, slong prec)

Set `res` to $L(s,\chi)$ by the Euler product, which is efficient for
large $\operatorname{Re}(s)$; the result is finite only if
$\operatorname{Re}(s) \ge 2$. The underscore version takes an integer
$s$ and a real character given by the list `chi` of its values at
$0, \ldots, \text{mod}-1$, and computes $1/L(s, \chi)$ if `reciprocal`
is nonzero.

::: details
The truncation error is bounded with `mag_hurwitz_zeta_uiui`: since
$1/L(s,\chi) = \prod_{p} (1 - \chi(p) p^{-s}) = \sum_{k \ge 1} \mu(k)\chi(k) k^{-s}$
and the truncated product contains all terms with smooth index,

$$
\left|\prod_{p < N} \left(1 - \frac{\chi(p)}{p^s}\right) - \frac{1}{L(s,\chi)}\right|
\le \sum_{k=N}^{\infty} \frac{1}{k^{\operatorname{Re}(s)}} = \zeta(\operatorname{Re}(s),N).
$$
:::
:::

::: function
void acb_dirichlet_l_fmpq(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)
void acb_dirichlet_l_fmpq_afe(acb_t res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

Set `res` to $L(s,\chi)$ for rational $s$. The `afe` version uses the
approximate functional equation, a rapidly convergent series of
incomplete gamma functions (evaluated by binary splitting and the
[bit-burst algorithm](algorithms/special-functions/error-functions#bit-burst)),
which is fastest at very high precision or large modulus; the default version uses
`acb_dirichlet_l_si` when $s$ is an integer, the approximate functional
equation when this is estimated to be faster, and otherwise calls
`_acb_dirichlet_l_generic`. Imprimitive characters are reduced to the
primitive characters inducing them, as in `acb_dirichlet_l_fmpq_afe_jet`.
:::

::: function
int _acb_dirichlet_l_fmpq_use_afe(ulong q, const fmpq_t s, slong prec)
int _acb_dirichlet_l_jet_use_afe(ulong q, const fmpq_t s, int deflate, slong len, slong prec)

Return whether the approximate functional equation is estimated to be
faster than the generic algorithms for $L(s,\chi)$, respectively for its
expansion to length `len`, where $q$ is the conductor of $\chi$
($q = 1$ for the Riemann zeta function). These are the heuristics used by
`acb_dirichlet_l`, `acb_dirichlet_l_fmpq`, `acb_dirichlet_l_jet`,
`acb_dirichlet_zeta` and `acb_dirichlet_zeta_jet` (and therefore by
`acb_dirichlet_stieltjes` and `arb_const_glaisher`, which call
`acb_dirichlet_zeta_jet`). The cutoffs come from timings (wall time)
with 1 and 8 threads. They decrease with $q$ (for $q > 50$, values at
half-integers, and for $q > 200$, jets, always use the approximate
functional equation), and with the number of threads $T$
(`flint_get_num_threads`), since the terms of the approximate functional
equation are computed in parallel with nearly perfect speedup while
Euler–Maclaurin summation and the Hurwitz zeta functions gain less from
threads. For $q = 1$ they are
about $4 \cdot 10^4 T^{-1/2}$ bits for $\zeta(1/2)$,
$10^5 T^{-0.9}$ bits for other values, $1.5 \cdot 10^5 T^{-1}$ bits for
jets at $s = 1/2$ and $4.4 \cdot 10^5 T^{-0.8}$ bits for other jets of
length 2 (somewhat higher for longer jets), the latter capped at
$3.5 \cdot 10^5$ bits because of the $O(p^2)$ memory of Euler–Maclaurin
summation (4 GB at $3.3 \cdot 10^5$ bits). The choice of algorithm, and
hence the radius of the output, can therefore depend on the number of
threads.
:::

::: function
int _acb_dirichlet_get_small_fmpq(fmpq_t x, const acb_t s)

If $s$ is an exact real number equal to a rational number with numerator
of at most 40 bits and denominator of at most 20 bits (i.e. a dyadic
number $m 2^{-e}$ with $e \le 20$), sets $x$ to it and returns 1;
otherwise returns 0.
:::

::: function
void acb_dirichlet_l_fmpq_afe_jet(acb_ptr res, const fmpq_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec)

Sets `res` to the Taylor coefficients $L(s), L'(s), \ldots, L^{(\text{len}-1)}(s)/(\text{len}-1)!$
for rational $s$, with `deflate` as in `acb_dirichlet_l_jet`.
This uses the same approximate functional equation and precision
management as `acb_dirichlet_l_fmpq_afe`, with every quantity replaced by a
power series in $x$ where $s \to s + x$. The required expansions
of the incomplete gamma functions $\Gamma(a + x/2, \pi n^2 / q)$ with respect to the
parameter are computed by jet versions of the binary splitting and
bit-burst algorithms (see `_arb_gamma_upper_fmpq_step_bsplit_jet`), with
truncation errors bounded separately for each coefficient.
For fixed `len`, the cost is $\tilde O(p^{1.5})$ at precision $p$,
roughly `len` times the cost of `acb_dirichlet_l_fmpq_afe`; this
makes it possible to compute, for example, Glaisher's constant (via
$\zeta'(-1)$), Stieltjes constants and their generalizations to Dirichlet
$L$-functions, or the second-order Landau-Ramanujan constant (via
$\zeta'(s)$ and $L'(s, \chi_{-4})$ at $s = 2^k$) in quasi-$p^{1.5}$ time.
An imprimitive character $\chi$ modulo $q$ induced by $\chi_0$ is handled by
evaluating $L(s + x, \chi_0)$ and multiplying by the Euler factors
$\prod_{p \mid q} (1 - \chi_0(p) p^{-s-x})$.
The deflation only has an effect for $q = 1$ (as in `acb_dirichlet_zeta_jet`);
for the Riemann zeta function at $s = 1$, `deflate` must be set, and for
a principal character modulo $q > 1$, `deflate` must not be set;
otherwise the output is indeterminate.

The terms of the approximate functional equation are independent, and are
computed in parallel when FLINT is allowed to use several threads
(see `flint_set_num_threads`). The terms are added up in a fixed order,
so the output does not depend on the number of threads. Each term is
added as soon as all earlier terms are done; a term that finishes early
waits in a small buffer (a few terms per thread), so the memory usage is
$O(t \cdot \text{len} \cdot p)$ bits with $t$ threads, independent of
the number of terms $O(\sqrt{q p})$. The same applies to
`acb_dirichlet_l_fmpq_afe`, which is the case `len = 1` of this function.
:::

::: function
void acb_dirichlet_l_vec_hurwitz(acb_ptr res, const acb_t s, const acb_dirichlet_hurwitz_precomp_t precomp, const dirichlet_group_t G, slong prec)

Sets `res` (of length $\phi(q)$) to $L(s,\chi)$ for all characters
$\chi$ modulo $q$, from the $\phi(q)$ values $\zeta(s, k/q)$ and one
discrete Fourier transform. The values are ordered by lexicographically
ordered Conrey logarithms, as in `acb_dirichlet_dft_index`. If `precomp`
is `NULL`, each Hurwitz zeta value is computed directly.
:::

::: function
void acb_dirichlet_l_jet(acb_ptr res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec)

Sets `res` to the Taylor coefficients $L(s), L'(s), \ldots, L^{(\text{len}-1)}(s)/(\text{len}-1)!$.
With `deflate` nonzero and principal $\chi$, the expansion of
$L(s,\chi) - \frac{\phi(q)}{q (s-1)}$ is computed instead, which at
$s = 1$ is the regular part of the Laurent expansion; for non-principal
$\chi$, `deflate` has no effect. For $q = 1$ this calls
`acb_dirichlet_zeta_jet`. Otherwise, if $s$ is an exact rational number of small height
(see `_acb_dirichlet_get_small_fmpq`), the approximate functional
equation `acb_dirichlet_l_fmpq_afe_jet` is used when
`_acb_dirichlet_l_jet_use_afe` estimates it to be faster (this is not
done for principal characters with deflation); otherwise the Hurwitz
zeta expansions are summed.
:::

::: function
void _acb_dirichlet_l_series(acb_ptr res, acb_srcptr s, slong slen, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec)
void acb_dirichlet_l_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, slong len, slong prec)

Set `res` to $L(s(x),\chi)$ truncated to length `len`, for a power series
$s(x)$, with `deflate` as in `acb_dirichlet_l_jet`.
:::

## Discrete Fourier transforms on characters {#dft}

The discrete Fourier transform of $f : \Z/q\Z \to \mathbb{C}$ is the function
$\hat f(\chi) = \sum_{x \bmod q} \overline{\chi(x)} f(x)$ on the
characters modulo $q$. Using the Conrey isomorphism between
$(\Z/q\Z)^\times$ and its dual, it is computed as
$g(a) = \sum_{b \bmod q} \overline{\chi_q(a,b)} f(b)$. Both functions
are thin wrappers around `gr_dft_acb_dirichlet_index` and
`gr_dft_acb_dirichlet` of [`gr_dft`](gr_dft), which compute a product DFT
over the cyclic components of the group in fixed-point arithmetic with
rigorous error bounds.

::: function
void acb_dirichlet_dft_index(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)

Sets `w` to the DFT of `v`, both vectors of length $\phi(q)$ indexed by
the lexicographically ordered Conrey logarithms (the order of
`dirichlet_char_next` and `dirichlet_char_index`).

::: details
For example, for $q = 15$ the elements are stored in the following
order:

| index | log $= [e, f]$ | number $= 7^e 11^f$ |
|---|---|---|
| 0 | $[0, 0]$ | 1 |
| 1 | $[0, 1]$ | 7 |
| 2 | $[0, 2]$ | 4 |
| 3 | $[0, 3]$ | 13 |
| 4 | $[1, 0]$ | 11 |
| 5 | $[1, 1]$ | 2 |
| 6 | $[1, 2]$ | 14 |
| 7 | $[1, 3]$ | 8 |
:::
:::

::: function
void acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const dirichlet_group_t G, slong prec)

Sets `w` to the DFT of `v`, both vectors of length $q$ indexed by Conrey
numbers. The entries at indices not coprime to $q$ are ignored.
:::

## Hardy Z-functions {#hardy}

For a primitive character $\chi$ (or for $\zeta$, with `G` and `chi`
`NULL`), the Hardy Z-function is $Z(t) = e^{i\theta(t)} L(1/2+it, \chi)$,
where the phase is

$$
\theta(t) = -\frac{t}{2} \log(\pi/q) - \frac{i \log(\epsilon)}{2}
+ \frac{\log \Gamma((s+\delta)/2) - \log \Gamma((1-s+\delta)/2)}{2i},
\qquad s = \tfrac12 + it.
$$

$Z(t)$ is real for real $t$, and its real zeros are the zeros of
$L(s, \chi)$ on the critical line. The functions throw an exception if
$\chi$ is not primitive.

::: function
void acb_dirichlet_hardy_theta(acb_ptr res, const acb_t t, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)

Sets `res` to the first `len` Taylor coefficients of $\theta(t)$, using
the log-gamma function.
:::

::: function
void acb_dirichlet_hardy_z(acb_ptr res, const acb_t t, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)

Sets `res` to the first `len` Taylor coefficients of $Z(t)$, computed
as $e^{i\theta(t)} L(1/2 + it, \chi)$ with `acb_dirichlet_l_jet`. For the
Riemann zeta function, this uses the Riemann–Siegel formula when $t$ is
large and `len` $\le 2$. If $t$ is real, the imaginary parts are set to
exactly zero.
:::

::: function
void _acb_dirichlet_hardy_theta_series(acb_ptr res, acb_srcptr s, slong slen, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)
void acb_dirichlet_hardy_theta_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)
void _acb_dirichlet_hardy_z_series(acb_ptr res, acb_srcptr s, slong slen, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)
void acb_dirichlet_hardy_z_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, slong len, slong prec)

Set `res` to $\theta(t(x))$, respectively $Z(t(x))$, truncated to length
`len`, for a power series $t(x)$ (given as `s`).
:::

## Gram points and counting zeros {#counting}

These functions concern the Riemann zeta function only. Let $N(t)$ be the
number of zeros of $\zeta(s)$, counted with multiplicity, with
$0 < \operatorname{Im}(s) \le t$, and $S(t) = N(t) - \theta(t)/\pi - 1$.
They implement Turing's method as described in
[](algorithms/special-functions/zeta-zeros#turing), following the
implementation in mpmath [@Joh2018b] by Arias de Reyna [@Ari2012].

::: function
void acb_dirichlet_gram_point(arb_t res, const fmpz_t n, const dirichlet_group_t G, const dirichlet_char_t chi, slong prec)

Sets `res` to the $n$-th Gram point $g_n$, the unique solution in
$[7, \infty)$ of $\theta(g_n) = \pi n$, for $n \ge -1$. Only the Riemann
zeta function is supported: `G` and `chi` must be `NULL`; otherwise, or
for $n < -1$, the result is indeterminate.
:::

::: function
ulong acb_dirichlet_turing_method_bound(const fmpz_t p)

Returns a number $B$ of consecutive good Gram blocks that suffices for
Turing's method [@Tur1953] with the bounds of [@Leh1970], [@Bre1979] and
[@Tru2011]: if $B$ consecutive Gram blocks with union $[g_n, g_p)$
satisfy Rosser's rule, then $N(g_n) \le n + 1$ and $N(g_p) \ge p + 1$.
:::

::: function
void acb_dirichlet_zeta_nzeros(arb_t res, const arb_t t, slong prec)
void _acb_dirichlet_exact_zeta_nzeros(fmpz_t res, const arf_t t)

Set `res` to $N(t)$. For exact $t$ the result is an exact integer; for
an inexact ball, it is the ball containing $N(t_1), \ldots, N(t_2)$ where
$[t_1, t_2]$ encloses $t$. The underscore version takes an exact $t$ and
sets an integer.
:::

::: function
void acb_dirichlet_zeta_nzeros_gram(fmpz_t res, const fmpz_t n)

Sets `res` to $N(g_n)$, for $n \ge -1$.
:::

::: function
void acb_dirichlet_backlund_s(arb_t res, const arb_t t, slong prec)

Sets `res` to $S(t) = \frac{1}{\pi}\operatorname{arg}\zeta(\frac{1}{2} + it)$,
where the argument is defined by continuous variation from $s = 2$ to
$2 + it$ and then to $\frac12 + it$, not the principal value. It is
computed as $N(t) - \theta(t)/\pi - 1$.
:::

::: function
void acb_dirichlet_backlund_s_bound(mag_t res, const arb_t t)

Sets `res` to an upper bound for $|S(t)|$, by Theorem 1 and (1.2) of
[@Tru2014].
:::

::: function
slong acb_dirichlet_backlund_s_gram(const fmpz_t n)

Returns $S(g_n)$, for $n \ge -1$.
:::

## Zeros of the Riemann zeta function {#zeros}

The $n$-th zero of the Hardy Z-function is $\gamma_n$, where
$\rho_n = 1/2 + i\gamma_n$ is the $n$-th nontrivial zero of $\zeta(s)$
with positive imaginary part (assuming, as verified for every $n$ handled
by these functions, that the zeros are on the critical line and simple).
[](acb_dirichlet/zeta-zeros) explains the methods and gives an
example.

::: function
void acb_dirichlet_hardy_z_zero(arb_t res, const fmpz_t n, slong prec)
void acb_dirichlet_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, slong prec)

Set `res` to $\gamma_n$, respectively to the `len` consecutive values
$\gamma_n, \ldots, \gamma_{n+\text{len}-1}$, accurate to about `prec`
bits. Throws an exception if $n < 1$. The zeros are isolated rigorously and computed
by the Riemann–Siegel method, or by Platt's method at large height when
that is expected to be faster.
:::

::: function
void acb_dirichlet_zeta_zero(acb_t res, const fmpz_t n, slong prec)
void acb_dirichlet_zeta_zeros(acb_ptr res, const fmpz_t n, slong len, slong prec)

Set `res` to $\rho_n = 1/2 + i\gamma_n$, respectively to `len`
consecutive zeros starting with $\rho_n$. Requires $n \ge 1$.
:::

::: function
void _acb_dirichlet_hardy_z_zeros_rs(arb_ptr res, const fmpz_t n, slong len, slong prec)

The same as `acb_dirichlet_hardy_z_zeros`, always with the Riemann–Siegel
method: isolation with Gram points and Turing's method, then refinement
with `_acb_dirichlet_refine_hardy_z_zero`.
:::

::: function
slong _acb_dirichlet_hardy_z_zeros_use_platt(const fmpz_t n, slong len, slong prec)

Returns 0 if `acb_dirichlet_hardy_z_zeros` uses the Riemann–Siegel method
for `len` zeros from the $n$-th at target precision `prec` (or at the
default precision $64 + \lceil \log_2 n \rceil$ if `prec` is not
positive), and otherwise the working precision of the large height
method. The choice minimises a cost model measured with 2 threads.

::: details
The model compares three options.

- The Riemann–Siegel method costs about the same per zero; beyond the
  default precision, the last stage of the refinement takes four
  evaluations of $Z$ at the target precision (about 6 times the cost of a
  zero at the default precision at 192 bits, 12 times at 320 bits).
- The large height method with the double-word sums over $j$, at a
  working precision of at most $\log_2 n + 104$ bits (about
  $\log_2 n + 68$ bits of accuracy), has a fixed cost of a few seconds up
  to $n = 10^{15}$ for one multi-evaluation, which isolates thousands of
  zeros, then about 10 ms per zero, plus the last stage of the refinement
  for each zero if a higher precision is requested.
- The large height method with ball arithmetic sums at working precision
  `prec` + 48 costs 3 to 12 times as much per call and is limited in
  accuracy (its parameters depend only on the height) to about 253, 234,
  219, 201, 187, 170 and 155 bits at $n = 10^{11}, \ldots, 10^{17}$.

At the default precision, the large height method is used from 158 zeros
at $n = 10^{11}$, 51 at $10^{12}$, 14 at $10^{13}$, 10 at $10^{14}$, 3 at
$10^{15}$, 2 from $10^{16}$, and (extrapolated) for a single zero from
$10^{20}$. Without double-word arithmetic (x87 builds), the
large height method is used from $10^{15}$, or from $10^{11}$ for more
than 100 zeros. It is only used for $10^{11} \le n < 10^{23}$.
:::
:::

::: function
void acb_dirichlet_isolate_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n)
void _acb_dirichlet_isolate_gram_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n)
void _acb_dirichlet_isolate_rosser_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n)
void _acb_dirichlet_isolate_turing_hardy_z_zero(arf_t a, arf_t b, const fmpz_t n)

Set $(a, b)$ to an interval that contains $\gamma_n$ and no other zero of
$Z$. The `gram` version uses Gram's law and requires $1 \le n \le 126$;
the `rosser` version uses Rosser's rule and requires
$1 \le n \le 13999526$; the `turing` version uses Turing's method and
requires $n \ge 2$; the default version chooses among them and requires
$n \ge 1$.
:::

::: function
int _acb_dirichlet_definite_hardy_z(arb_t res, const arf_t t, slong * pprec)

Sets `res` to $Z(t)$, increasing the precision `*pprec` as needed until
the sign of $Z(t)$ is determined, and returns the sign.
:::

::: function
void _acb_dirichlet_refine_hardy_z_zero(arb_t res, const arf_t a, const arf_t b, slong prec)

Sets `res` to the unique zero of $Z$ in $(a, b)$, to `prec` bits.

::: details
Below $4 \log_2 t + 40$ bits, the Illinois method is used, with iterates
at $\text{prec} + \log_2 t + 8$ bits but $Z$ evaluated at only
$\text{prec} + 12$ bits (its value near the zero is needed to about
$\text{prec} - \log_2 t$ bits after the binary point), so that the
Riemann–Siegel main sum stays within reach of `dfloat_powsum_sieved`;
when an iterate is so close to the zero that the sign of $Z$ is
undetermined, the zero is bracketed by the points at a quarter of the
tolerance on either side instead of raising the precision. From
$\log_2 t + 96$ bits, the zero is first refined to
$p_1 = (\text{prec} + \log_2 t)/2 + 16$ bits, where $Z$ is much cheaper,
and finished by one secant step through the endpoints of that bracket and
the signs of $Z$ at a quarter of the tolerance on either side of the
result: four evaluations at the target precision rather than a dozen for
the Illinois method alone. Above $4 \log_2 t + 40$ bits, Newton's method is used.
:::
:::

::: function
void _acb_dirichlet_refine_hardy_z_zero_ball(arb_t res, const arb_t z, slong prec)

Sets `res` to the zero of $Z$ in the ball `z`, which must contain exactly
one zero (for example a zero computed by Platt's method), to `prec` bits:
`z` rounded if it is accurate enough, otherwise the last stage of the
refinement above if `z` is accurate to about half the precision, and the
full refinement otherwise.
:::

## Platt's method for zeros at large height {#platt}

Platt's method [@Pla2017] evaluates the real function
$\Lambda(t) e^{\pi t/4}$, which has the same zeros as $Z(t)$, on a grid
of $AB$ points with spacing $1/A$ by a convolution computed with fast
Fourier transforms, then interpolates between the grid points. One
multi-evaluation isolates thousands of zeros; see
[](algorithms/special-functions/zeta-zeros#platt). The functions without
underscore choose all parameters from the height and expect
$10^4 \le n < 10^{37}$; for smaller $n$ they return 0. They use
`flint_get_num_threads()` threads.

::: function
void acb_dirichlet_platt_scaled_lambda(arb_t res, const arb_t t, slong prec)

Sets `res` to $\Lambda(t) e^{\pi t/4}$, where
$\Lambda(t) = \pi^{-it/2} \Gamma(\frac{1/2+it}{2}) \zeta(\frac12 + it)$
is defined at the beginning of §3 of [@Pla2017]. By the functional
equation it is real for real $t$, it has the same zeros as
$\zeta(1/2 + it)$ [@Pla2011], and the exponential factor cancels the
decay of the gamma factor.
:::

::: function
void acb_dirichlet_platt_scaled_lambda_vec(arb_ptr res, const fmpz_t T, slong A, slong B, slong prec)

Sets `res` to the values of `acb_dirichlet_platt_scaled_lambda` at the
$N = AB$ grid points $T - B/2 + k/A$, $0 \le k < N$, by evaluating each
point separately. $N$ must be even.
:::

::: function
void acb_dirichlet_platt_multieval(arb_ptr out, const fmpz_t T, slong A, slong B, const arb_t h, const fmpz_t J, slong K, slong sigma, slong prec)
void acb_dirichlet_platt_multieval_threaded(arb_ptr out, const fmpz_t T, slong A, slong B, const arb_t h, const fmpz_t J, slong K, slong sigma, slong prec)

Set `out` to the values of `acb_dirichlet_platt_scaled_lambda` on the
same grid as `acb_dirichlet_platt_scaled_lambda_vec`, all at once with
discrete Fourier transforms, in the notation of [@Pla2017]: $h$ is the
width of the Gaussian, $J$ the number of terms of the Dirichlet series,
$K$ the number of Taylor terms and `sigma` an odd integer used in the
error bounds. The `threaded` version uses `flint_get_num_threads()`
threads; the default version decides whether to use threads.

::: details
The default version computes the sums over $j \le J$ (the bulk of the
work) with `_dfloat_platt_smk_dd` ([](dfloat#powsum)) when `prec` is at
most $\log_2 T + 110$ bits, which gives a grid about as accurate as the ball arithmetic version
at $\log_2 T + 106$ bits; setting `acb_dirichlet_platt_use_dfloat` to 0
forces ball arithmetic. In both versions, the part after these sums works
at $\text{prec} - \log_2 T + 16$ bits (its inputs, of size about 1,
carry that accuracy), generates the rows of the gamma factor table for
each group of four values of $k$ from $N$ stored factors and bases
instead of storing $KN$ entries, and sums the $K$ convolutions in the
frequency domain, in parallel over the groups with one group per thread
at a time. The result does not depend on the number of threads.
:::
:::

::: variable
int acb_dirichlet_platt_use_dfloat

If nonzero (the default), `acb_dirichlet_platt_multieval` computes the
sums over $j$ with double-word arithmetic when the precision allows it.
Set it to 0 to always use ball arithmetic.
:::

::: function
void acb_dirichlet_platt_ws_interpolation(arb_t res, arf_t deriv, const arb_t t0, arb_srcptr p, const fmpz_t T, slong A, slong B, slong Ns_max, const arb_t H, slong sigma, slong prec)

Sets `res` to `acb_dirichlet_platt_scaled_lambda` at $t_0$, by
Gaussian-windowed Whittaker–Shannon interpolation of the grid values `p`.
If `deriv` is not `NULL`, it is set to an approximation of the
derivative. `Ns_max` is the maximum number of grid points used on either
side of $t_0$, $H$ is the standard deviation of the Gaussian window, and
`sigma` is an odd positive integer used in the error bounds.
:::

::: function
slong acb_dirichlet_platt_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, slong prec)
slong acb_dirichlet_platt_local_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, slong prec)
slong _acb_dirichlet_platt_local_hardy_z_zeros(arb_ptr res, const fmpz_t n, slong len, const fmpz_t T, slong A, slong B, const arb_t h, const fmpz_t J, slong K, slong sigma_grid, slong Ns_max, const arb_t H, slong sigma_interp, slong prec)

Set at most `len` entries of `res` to consecutive zeros $\gamma_n, \ldots$
of the Hardy Z-function and return the number found. The `local`
versions make one grid evaluation; `acb_dirichlet_platt_hardy_z_zeros`
repeats them until `len` zeros are found. The parameters of the
underscore version have the meanings of `acb_dirichlet_platt_multieval`
and `acb_dirichlet_platt_ws_interpolation`.
:::

::: function
slong acb_dirichlet_platt_isolate_local_hardy_z_zeros(arf_interval_ptr res, const fmpz_t n, slong len, slong prec)
slong _acb_dirichlet_platt_isolate_local_hardy_z_zeros(arf_interval_ptr res, const fmpz_t n, slong len, const fmpz_t T, slong A, slong B, const arb_t h, const fmpz_t J, slong K, slong sigma_grid, slong Ns_max, const arb_t H, slong sigma_interp, slong prec)

Like the `local` functions above, but set `res` to isolating intervals of
the zeros instead of refining them.
:::

::: function
slong acb_dirichlet_platt_zeta_zeros(acb_ptr res, const fmpz_t n, slong len, slong prec)

Sets at most `len` entries of `res` to consecutive zeros $\rho_n, \ldots$
of $\zeta(s)$ by Platt's method, and returns the number found.
:::

## The secondary zeta function {#secondary}

::: function
void acb_dirichlet_secondary_zeta(acb_t res, const acb_t s, slong prec)

Sets `res` to the secondary zeta function $Z(s) = \sum_{n=1}^{\infty} \gamma_n^{-s}$,
where $\rho_n = 1/2 + i\gamma_n$ runs over the zeros of $\zeta(s)$ with
positive imaginary part. The Riemann hypothesis is not assumed. The
algorithm follows [@Ari2020], with the decomposition
$Z(s) = A(s) - P(s) + E(s) - S(s)$ depending on a parameter $a > 0$, an
improved tail bound for $P$, and rigorous bounds for the tails of $A$ and
$S$ that replace the heuristic bounds of [@Ari2020].
:::

## Internal functions {#internal}

::: function
void _acb_dirichlet_theta_argument_at_arb(arb_t xt, ulong q, const arb_t t, slong prec)

Sets `xt` to $\pi t^2 / q$, the exponent of the theta series.
:::

::: function
void _acb_dirichlet_theta_arb_series(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, const arb_t xt, slong len, slong prec)

Sets `res` to $\sum_{n=1}^{\text{len}-1} \chi(n) n^\delta x^{n^2}$ for
$x$ = `xt`, using `acb_dirichlet_qseries_arb`.
:::

::: function
void mag_tail_kexpk2_arb(mag_t res, const arb_t a, ulong n)

Sets `res` to an upper bound for $\sum_{k > n} k e^{-a k^2}$.
:::

::: function
void acb_dirichlet_si_poly_evaluate(acb_t res, slong * v, slong len, const acb_t z, slong prec)

Sets `res` to the polynomial with the `len` coefficients `v` evaluated at
$z$, by the baby-step giant-step method.
:::

::: function
ulong jacobi_one_prime(ulong p, ulong e, ulong pe, ulong cond)

Returns the local factor at $p$ (with $p^e \| q$ and `pe` $= p^e$) of
the Jacobi sum $J_q(1, \chi)$ for a character of conductor `cond`. The
value is a signed integer cast to `ulong`.
:::

::: function
void acb_dirichlet_platt_beta(arb_t res, const arb_t t, slong prec)

Sets `res` to $\beta(t) = 1/6 + \log\log t / \log t$, the exponent in the
bound $|\zeta(1/2 + it)| = O(t^{\beta(t)})$ used in [@Pla2017].
:::

::: function
void acb_dirichlet_platt_lemma_32(arb_t out, const arb_t h, const arb_t t0, const arb_t x, slong prec)
void acb_dirichlet_platt_lemma_A5(arb_t out, slong B, const arb_t h, slong k, slong prec)
void acb_dirichlet_platt_lemma_A7(arb_t out, slong sigma, const arb_t t0, const arb_t h, slong k, slong A, slong prec)
void acb_dirichlet_platt_lemma_A9(arb_t out, slong sigma, const arb_t t0, const arb_t h, slong A, slong prec)
void acb_dirichlet_platt_lemma_A11(arb_t out, const arb_t t0, const arb_t h, slong B, slong prec)
void acb_dirichlet_platt_lemma_B1(arb_t out, slong sigma, const arb_t t0, const arb_t h, const fmpz_t J, slong prec)
void acb_dirichlet_platt_lemma_B2(arb_t out, slong K, const arb_t h, const arb_t xi, slong prec)

Set `out` to the error bounds of Lemma 3.2 and Lemmas A.5, A.7, A.9,
A.11, B.1 and B.2 of [@Pla2017], which bound the errors of the grid
evaluation in `acb_dirichlet_platt_multieval` (truncation of the
Dirichlet series, of the Taylor expansion, aliasing of the discrete
Fourier transform, and so on). Lemma A.5 requires $B > h\sqrt{k}$; when
a lemma's hypotheses fail, the bound is infinite.
:::

::: type
acb_dirichlet_platt_c_precomp_struct
acb_dirichlet_platt_c_precomp_t
acb_dirichlet_platt_i_precomp_struct
acb_dirichlet_platt_i_precomp_t
acb_dirichlet_platt_ws_precomp_struct
acb_dirichlet_platt_ws_precomp_t

Precomputed constants of the interpolation error bounds.
:::

::: function
void acb_dirichlet_platt_c_precomp_init(acb_dirichlet_platt_c_precomp_t pre, slong sigma, const arb_t h, ulong k, slong prec)
void acb_dirichlet_platt_c_precomp_clear(acb_dirichlet_platt_c_precomp_t pre)
void acb_dirichlet_platt_c_bound_precomp(arb_t res, const acb_dirichlet_platt_c_precomp_t pre, slong sigma, const arb_t t0, const arb_t h, slong k, slong prec)
void acb_dirichlet_platt_c_bound(arb_t res, slong sigma, const arb_t t0, const arb_t h, slong k, slong prec)

Set `res` to the bound $C(\sigma, t_0, h, k)$ of [@Pla2017] for the
$k$-th derivative of the scaled function in the strip used for the error
analysis, with or without precomputed constants that depend only on
$\sigma$, $h$ and $k$.
:::

::: function
void acb_dirichlet_platt_i_precomp_init(acb_dirichlet_platt_i_precomp_t pre, slong A, const arb_t H, slong sigma, slong prec)
void acb_dirichlet_platt_i_precomp_clear(acb_dirichlet_platt_i_precomp_t pre)
void acb_dirichlet_platt_i_bound_precomp(arb_t res, const acb_dirichlet_platt_i_precomp_t pre_i, const acb_dirichlet_platt_c_precomp_t pre_c, const arb_t t0, slong A, const arb_t H, slong sigma, slong prec)
void acb_dirichlet_platt_i_bound(arb_t res, const arb_t t0, slong A, const arb_t H, slong sigma, slong prec)

Set `res` to the bound $I$ of [@Pla2017] for the error of the
Whittaker–Shannon interpolation caused by the Gaussian window and the
sampling rate $A$, with or without precomputed constants.
:::

::: function
void acb_dirichlet_platt_ws_precomp_init(acb_dirichlet_platt_ws_precomp_t pre, slong A, const arb_t H, slong sigma, slong prec)
void acb_dirichlet_platt_ws_precomp_clear(acb_dirichlet_platt_ws_precomp_t pre)
void acb_dirichlet_platt_ws_interpolation_precomp(arb_t res, arf_t deriv, const acb_dirichlet_platt_ws_precomp_t pre, const arb_t t0, arb_srcptr p, const fmpz_t T, slong A, slong B, slong Ns_max, const arb_t H, slong sigma, slong prec)

The interpolation `acb_dirichlet_platt_ws_interpolation` with the
constants of both bounds precomputed, for many evaluations with the same
$A$, $H$ and $\sigma$.
:::

::: function
void acb_dirichlet_platt_bound_C3(arb_t res, const arb_t t0, slong A, const arb_t H, slong Ns, slong prec)

Sets `res` to the bound $C_3$ of [@Pla2017] for the error of truncating
the interpolation sum to $N_s$ points on either side.
:::
