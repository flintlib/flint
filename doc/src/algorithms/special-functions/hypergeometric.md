---
title: Hypergeometric functions
summary: The general machinery for hypergeometric functions in ball arithmetic — convergent series with rigorous tail bounds (also for power series parameters), asymptotic expansions with Olver's error bounds, the choice between them by estimated cancellation, connection formulas and regularization, and the summation algorithms (forward recurrence, rectangular splitting, binary splitting, fast multipoint evaluation, bit-burst continuation).
complexity: N terms of a series at p bits — O(N) full multiplications forward, O(N^{1/2}) full and O(N) scalar multiplications by rectangular splitting, O(M(p) log² p) by binary splitting for arguments with O(1) bits, and O(M(p) log² p) for any argument with the bit-burst algorithm (holonomic functions).
order: 40
---

A hypergeometric series is a series $\sum_k T(k)$ whose term ratio
$T(k+1)/T(k)$ is a rational function of $k$. In standard form,
$$
{}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z) = \sum_{k=0}^\infty
\frac{(a_1)_k\cdots(a_p)_k}{(b_1)_k\cdots(b_q)_k} \frac{z^k}{k!}.
$$
Most classical special functions are hypergeometric series, or finite
combinations of such series with exponential, power and gamma function
prefactors. Evaluating them rigorously to arbitrary precision for all
complex parameters and arguments requires three ingredients: the
convergent series with a bound for the remainder, asymptotic
expansions with error bounds for large arguments, and transformations
and connection formulas that map the remaining cases to one of the
first two. The implementation in [`acb_hypgeom`](acb_hypgeom) is
described in Johansson [@Joh2016], from which the error bounds below
are taken; the summation algorithms are surveyed in [@Joh2014a].

Throughout this section,
$$
T(k) = \frac{\prod_{i=0}^{p-1} (a_i)_k}{\prod_{i=0}^{q-1} (b_i)_k} z^k,
\qquad {}_pf_q(z) = \sum_{k=0}^\infty T(k) = {}_{p+1}F_q(a_0, \ldots,
a_{p-1}, 1; b_0, \ldots, b_{q-1}; z),
$$
so that ${}_pF_q$ is ${}_pf_{q+1}$ with the extra lower parameter
$b_q = 1$. This is the convention of the low-level functions
`acb_hypgeom_pfq_sum` and `acb_hypgeom_pfq_direct`.

## Convergent series {#tail-bound}

If $p \le q$, the series ${}_pf_q(z)$ converges for all $z$, and the
terms eventually decrease faster than geometrically; if $p = q + 1$, it
converges for $|z| < 1$. A bound for the tail $\sum_{k \ge n} T(k)$ in
terms of the first omitted term $T(n)$ follows from a bound for the
term ratio.

::: theorem Tail bound {#thm-tail}
Let $p \le q$ and $n \ge 0$ with $\operatorname{Re}(b_i + n) > 0$ for
all $i$, and let
$$
D = |z| \prod_{i=0}^{p-1} \left(1 + \frac{|a_i - b_i|}{|b_i + n|}\right)
\prod_{i=p}^{q-1} \frac{1}{|b_i + n|}.
$$
If $D < 1$, then $\left|\sum_{k=n}^{\infty} T(k)\right| \le
\frac{|T(n)|}{1 - D}$.
:::

::: proof
For $k \ge n$, $T(k+1)/T(k) = z \prod_{i<p} \frac{a_i + k}{b_i + k}
\prod_{p \le i < q} \frac{1}{b_i+k}$. Since $\operatorname{Re}(b_i + n)
> 0$, $|b_i + k| \ge |b_i + n|$ for $k \ge n$ (adding a positive real
number to a number with positive real part increases its modulus).
Hence
$$
\left|\frac{a_i + k}{b_i + k}\right| = \left|1 + \frac{a_i - b_i}{b_i + k}\right|
\le 1 + \frac{|a_i - b_i|}{|b_i + n|}, \qquad \frac{1}{|b_i + k|} \le
\frac{1}{|b_i + n|},
$$
so $|T(k+1)| \le D |T(k)|$ for all $k \ge n$, and the tail is bounded
by $|T(n)| \sum_{j \ge 0} D^j$.
:::

`acb_hypgeom_pfq_bound_factor` computes $C = (1-D)^{-1}$ with upper
bounds in `mag_t` arithmetic, pairing the upper parameter $a_i$ with
the lower parameter $b_i$; $C$ is infinite if $p > q$, if some
$\operatorname{Re}(b_i + n) \le 0$, or if $D \ge 1$. The bound can
therefore be infinite for a small $n$ even if the series converges;
the number of terms must be large enough that the terms have started
to decrease. `acb_hypgeom_pfq_direct` sums $n$ terms with
`acb_hypgeom_pfq_sum`, which also returns $T(n)$, and adds $C|T(n)|$
to the radius. If $T(n)$ is exactly zero, the series has terminated
and no error is added; this is how polynomials such as ${}_pF_q$ with a
nonpositive integer upper parameter are computed for any $p$, $q$ and
$z$.

### Power series parameters {#series-tail}

The same principle gives bounds for each coefficient of $\sum_{k \ge
N} T(k) \in \mathbb{C}[[x]]/(x^n)$ when the parameters $a_i, b_i$ and the
argument $z$ are themselves power series, which is how derivatives
with respect to parameters and limits at singular parameter values are
computed. Write $A_{[k]}$ for the coefficient of $x^k$ in $A$,
$A_{[m:n]} = \sum_{k=m}^{n-1} A_{[k]} x^k$, $|A| = \sum_k |A_{[k]}|
x^k \in \mathbb{R}_{\ge 0}[[x]]$, $A \le B$ if $|A|_{[k]} \le |B|_{[k]}$ for all
$k$, and $\mathcal{R}(B) = |B_{[0]}| - |B_{[1:\infty]}|$.

::: lemma {#lem-majorants}
For power series $A$ and $B$ with $B_{[0]} \ne 0$,
$$
|A + B| \le |A| + |B|, \qquad |AB| \le |A| |B|, \qquad
|A / B| \le |A| / \mathcal{R}(B).
$$
:::

::: proof
The first two follow from $(AB)_{[k]} = \sum_{j \le k} A_{[j]}
B_{[k-j]}$ and the triangle inequality. For the third, write $B =
B_{[0]}(1 - E)$ with $E = -B_{[1:\infty]}/B_{[0]}$; then $1/B =
B_{[0]}^{-1} \sum_j E^j$ and $|1/B| \le |B_{[0]}|^{-1} \sum_j |E|^j =
1/(|B_{[0]}| - |B_{[1:\infty]}|)$, where the right side is expanded as
a power series in $x$.
:::

::: theorem {#thm-series-tail}
If $p \le q$ and $\operatorname{Re}(b_{i[0]} + N) > 0$ for all $i$, let
$$
D = |z| \prod_{i=0}^{p-1} \left(1 + \frac{|a_i-b_i|}{\mathcal{R}(b_i+N)}\right)
\prod_{i=p}^{q-1} \frac{1}{\mathcal{R}(b_i + N)} \in \mathbb{R}_{\ge 0}[[x]].
$$
If $D_{[0]} < 1$, then $\sum_{k=N}^\infty T(k) \le (1 - D)^{-1} |T(N)|$
coefficientwise.
:::

The proof is that of [](#thm-tail), with [](#lem-majorants) in place of
the scalar inequalities. In the implementation
(`acb_hypgeom_pfq_series_direct`), sums and products of the majorant
series are computed with point-valued upper bounds for the absolute
values, so no interval arithmetic is needed; for $\mathcal{R}(B)$ a
lower bound for $|B_{[0]}|$ and upper bounds for the other
coefficients are used.

### Number of terms and working precision {#choose-n}

`acb_hypgeom_pfq_choose_n` chooses $n$ by simulating the term
magnitudes in double precision: it stops when the current term is
`prec` bits smaller than the largest term seen so far, since the
cancellation between large terms has to be compensated anyway. This
also picks a sensible truncation point for divergent asymptotic series
(where the terms decrease and then grow), stops at the end of a
terminating series, and steps over parameters close to nonpositive
integers. The functions built on the series also estimate the
cancellation and raise the working precision: the sum of a series
whose largest term is $2^m$ times its value loses about $m$ bits. For
example ${}_1F_1(a; b; z)$ with $z = -x < 0$ has terms of size up to
about $e^{x}$ but a value of size $x^{-a}$, so Kummer's transformation
$M(a,b,z) = e^z M(b-a, b, -z)$ is applied first, and ${}_0F_1(a; z)$
loses about $2\sqrt{|z|}/\log 2$ bits for negative $z$.

## Summation algorithms {#summation}

Let $T(k+1) = T(k)\, P(k)/Q(k)$ with polynomials $P$ and $Q$ (absorbing
$z$ into $P$). `acb_hypgeom_pfq_sum` and the variants for rational
parameters in [](algorithms/special-functions/hypergeometric-sums)
choose among four ways to compute $S = \sum_{k<n} T(k)$.

**Forward recurrence.** Compute $T(k)$ and the partial sums term by
term: $n$ multiplications and $n$ divisions at full precision.
Divisions can be delayed by accumulating numerators and denominators
separately (as `arb_hypgeom_sum_fmpq_arb_forward` does). This is best
for few terms or when the parameters are full-precision numbers.

**Rectangular splitting.** When the parameters are short (integers or
rationals with small numerators and denominators) but $z$ is a
full-precision number, the expensive operations are the
multiplications by $z$. Writing the sum in Horner form backwards,
$$
S = 1 + \frac{P(0)}{Q(0)}\left(1 + \frac{P(1)}{Q(1)}\left(1 + \cdots\right)\right),
$$
and grouping the powers of $z$ into blocks of length $m$, the sum
becomes a polynomial in $z^m$ whose coefficients are polynomials in
$z$ of degree less than $m$ with short coefficients. With $z, \ldots,
z^m$ precomputed, the full multiplications number about $m + n/m$, and
the remaining $O(n)$ operations are multiplications and divisions by
short numbers. This is Smith's method of concurrent series
[@Smi1989], a variant of the Paterson–Stockmeyer algorithm for
polynomial evaluation [@PS1973]; Johansson [@Joh2014b] generalizes it to
parameters that are themselves full-precision numbers by expanding the
products of consecutive factors as polynomials in the parameter.
`acb_hypgeom_pfq_sum_rs` applies it with $m \approx \sqrt n$ (at most
150), and `arb_hypgeom_sum_fmpq_arb_rs` combines up to 8 consecutive
term ratios over a common integer denominator.

**Binary splitting.** The recurrence
$$
\begin{pmatrix} T(k+1) \\ S(k+1) \end{pmatrix} =
\frac{1}{Q(k)} \begin{pmatrix} P(k) & 0 \\ Q(k) & Q(k) \end{pmatrix}
\begin{pmatrix} T(k) \\ S(k) \end{pmatrix}
$$
gives $S$ from the product of $n$ matrices, which binary splitting
computes as a balanced product tree
([](algorithms/integers/binary-splitting)) [@Bre1976; @CC1988;
@BZ2011]. If $z$ and the parameters have $O(\log n)$ bits, the entries
of a subproduct of $\ell$ matrices have $O(\ell \log n)$ bits, and the
whole product costs $O(M(n \log n) \log n)$; for $p$-bit precision
with $n = O(p/\log p)$ terms, this is $O(M(p) \log^2 p)$, quasi-linear
in $p$. If $z$ has $p$ bits, the entries grow to $p$ bits after few
levels and the method degrades to $O(n)$ full multiplications. In
complex ball arithmetic, binary splitting also limits the
overestimation of the radii (the wrapping effect), since each number
goes through only $O(\log n)$ operations, so
`acb_hypgeom_pfq_sum` prefers it to the forward recurrence for $n > 8$
even for long parameters. `acb_hypgeom_pfq_sum_bs_invz` puts $1/z$
into the denominators, which helps when $1/z$ is short, as in
asymptotic expansions at integer arguments.

**Fast multipoint evaluation.** If $P$ and $Q$ have degree $d$, the
product of the first $m \approx \sqrt n$ matrices can be computed as a
matrix of polynomials of degree $O(dm)$ in the shift $k$; evaluating it
at the $w \approx n/m$ points $0, m, 2m, \ldots$ with fast multipoint
evaluation and multiplying the $w$ resulting matrices gives the whole
product with $\tilde O(\sqrt n)$ full-precision operations, instead of
$O(\sqrt n)$ full and $O(n)$ scalar operations [@CC1988; @BGS2007]. `acb_hypgeom_pfq_sum_fme` implements this
baby-step giant-step method for $p, q \le 3$, using product trees for
the evaluation (`_acb_poly_evaluate_vec_fast_precomp`). The numerical
stability of the polynomial arithmetic in ball arithmetic is poorer,
so it pays off only at high precision (from 1500 bits, for many
terms).

`acb_hypgeom_pfq_sum` uses rectangular splitting when the parameters
are short compared to the precision (binary splitting if moreover $z$
is short), fast multipoint evaluation for many terms at high precision,
and otherwise binary splitting or, for at most 8 terms, the forward
recurrence. The power series versions
(`acb_hypgeom_pfq_series_sum_forward`, `_bs`, `_rs`) perform the same
algorithms with truncated power series arithmetic.

## The bit-burst algorithm {#bit-burst}

Binary splitting is quasi-linear only for short arguments. For a
holonomic function $f$ — one that satisfies a linear differential
equation with polynomial coefficients, like every ${}_pF_q$ — the
bit-burst algorithm [@CC1990; @vdH1999; @Mez2010] removes this
restriction. Write the $p$-bit argument as $z = z_0 + (z_1 - z_0) +
(z_2 - z_1) + \cdots$ where $z_j$ is $z$ truncated to $2^{j} c$ bits.
The Taylor coefficients of $f$ at $z_j$ satisfy a recurrence with
polynomial coefficients in $k$ and $z_j$, and the step from $z_j$ to
$z_{j+1}$ needs only $O(p/2^j)$ terms of a series whose argument
$z_{j+1} - z_j$ has $2^{j+1}c$ bits and magnitude at most $2^{-2^jc}$.
The total size of the numbers in the binary splitting tree of each
step is therefore $O(p)$, and each step costs $O(M(p) \log p)$
(evaluating the matrix of a fundamental system of solutions); with
$O(\log p)$ steps this gives $O(M(p) \log^2 p)$ overall. The
error bounds for the local Taylor series come from Cauchy majorants of
the differential equation [@Hoe2001; @vdH1999].

FLINT uses the bit-burst algorithm for the dilogarithm at very high
precision (`acb_hypgeom_dilog_bitburst`, see
[](algorithms/special-functions/dilogarithm)), for the error function
and the incomplete gamma function $\Gamma(a, z)$ with rational $a$
above 30000 bits (`arb_hypgeom_erf_bb`, see
[](algorithms/special-functions/error-functions)), and in a simple form
(one or two continuation steps with Taylor series, without the
doubling of bits) for the corner case of ${}_2F_1$
([](algorithms/special-functions/gauss-2f1)). The elementary functions
and constants, which are the most important holonomic functions, are
covered in [](algorithms/special-functions/elementary) and
[](algorithms/special-functions/constants).

## Asymptotic expansions {#asymptotic}

For large $|z|$ the convergent series needs about $|z|$ terms and
suffers from cancellation of order $e^{|z|}$. Asymptotic expansions
then converge quickly in the sense that the first terms decrease
rapidly, but they diverge, so the error has to be bounded explicitly
for each truncation. All confluent hypergeometric functions reduce to
the function $U(a,b,z)$, the confluent hypergeometric function of the
second kind with its principal branch cut on $(-\infty, 0]$. With $U^*(a,b,z)
= z^a U(a,b,z)$,
$$
U^{*}(a,b,z) \sim {}_2F_0(a, a-b+1; ; -1/z) = \sum_{k=0}^{n-1}
\frac{(a)_k (a-b+1)_k}{k! (-z)^k} + \varepsilon_n(z).
$$
The error term is bounded following Olver [@NIST2012, §13.7], with a
case distinction according to the region $R$ in which $z$ lies; the
formula for the bound increases with $R$, so if $z$ is in two regions
the larger index can be used. Let $r = |b - 2a|$.

- If $\operatorname{Re}(z) \ge r$, then $R = 1$.
- Otherwise, if $|\operatorname{Im}(z)| \ge r$, or $\operatorname{Re}(z)
  \ge 0$ and $|z| \ge r$, then $R = 2$.
- Otherwise, if $|z| \ge 2r$, then $R = 3$.
- Otherwise the bound is infinite.

With the auxiliary quantities
$$
\sigma = |(b-2a)/z|, \qquad
\nu = \left(\tfrac{1}{2} + \tfrac{1}{2}\sqrt{1-4\sigma^2}\right)^{-1/2}
\le 1 + 2 \sigma^2, \qquad
\chi(n) = \sqrt{\pi}\, \frac{\Gamma(\tfrac{1}{2}n+1)}{\Gamma(\tfrac{1}{2} n
+ \tfrac{1}{2})},
$$
$$
C_n = \begin{cases}
1 & R = 1 \\
\chi(n) & R = 2 \\
(\chi(n) + \sigma \nu^2 n) \nu^n & R = 3,
\end{cases} \qquad
\sigma' = \begin{cases}
\sigma & R \ne 3 \\
\nu \sigma & R = 3,
\end{cases}
$$
$$
\alpha = (1 - \sigma')^{-1}, \qquad
\rho = \tfrac{1}{2} |2a^2-2ab+b| + \sigma' (1+ \tfrac{1}{4} \sigma')
(1-\sigma')^{-2},
$$
the error is bounded by
$$
|\varepsilon_n(z)| \le 2 \alpha C_n \left|\frac{(a)_n (a-b+1)_n}{n!
z^n} \right| \exp\left(\frac{2 \alpha \rho C_1}{|z|}\right).
$$
`acb_hypgeom_u_asymp` evaluates this bound with `mag_t` arithmetic,
chooses the number of terms with `acb_hypgeom_pfq_choose_n` (raising
the target by the size of the exponential factor), and sums the series
in powers of $w = -1/z$ with `acb_hypgeom_pfq_sum_invz`. If $a$ or $a -
b + 1$ is a nonpositive integer, the series terminates and is exact,
for every $z \ne 0$.

The Bessel, Airy, error, incomplete gamma, exponential and
trigonometric integral and Coulomb functions all use this expansion
through their representations by $U$. The Airy functions have
dedicated bounds, given in
[](algorithms/special-functions/bessel-airy#airy-asymp).

## Confluent hypergeometric functions {#confluent}

The convergent and the asymptotic expansions are tied together by
connection formulas. For $z \ne 0$ and $b \notin \Z$ [@NIST2012,
13.2.42],
$$
U(a,b,z) = \frac{\Gamma(1-b)}{\Gamma(a-b+1)} M(a,b,z)
         + \frac{\Gamma(b-1)}{\Gamma(a)} z^{1-b} M(a-b+1,2-b,z),
$$
and for all $z \ne 0$
$$
\frac{{}_1F_1(a,b,z)}{\Gamma(b)} = \frac{(-z)^{-a}}{\Gamma(b-a)}
U^{*}(a,b,z) + \frac{z^{a-b} e^z}{\Gamma(a)} U^{*}(b-a,b,-z),
$$
which is equivalent to [@NIST2012, 13.2.41] but simpler in form. The
first gives $U$ for small $|z|$ from two convergent series
(`acb_hypgeom_u_1f1`); the second gives $M$ for large $|z|$ from two
asymptotic series (`acb_hypgeom_m_asymp`).

**Choosing between them.** For $|z|$ below 4 the convergent series is
used, and for $|z| > 2^{64}$ the asymptotic one. In between,
`acb_hypgeom_m` and `acb_hypgeom_u` compare the accuracy that the
asymptotic expansion can reach, about $|z|/\log 2$ bits (the smallest
term of ${}_2F_0(\cdot;-1/z)$ is about $e^{-|z|}$ for small
parameters), with the accuracy of the convergent series after the
estimated cancellation, and choose the method and the working
precision accordingly. If the asymptotic expansion reaches the target
precision it is preferred, since it is cheaper. For large parameters
the asymptotic expansion is useless, and the convergent series needs
a very high working precision; for real arguments the functions then
try rigorous numerical integration ([](#integration)).

**Limits at integer parameters.** The formula for $U$ has removable
singularities at integer $b$. `acb_hypgeom_u_1f1` handles this case by
evaluating the formula with $b$ replaced by the power series $b + x$
and dividing out the zero of the denominator $\sin(\pi b)$, using the
power series algorithms and tail bounds above; the regularized
functions ${}_1\mathbf{F}_1 = {}_1F_1/\Gamma(b)$ are summed so that
the factor $1/\Gamma(b + k)$ vanishes for the terms that would divide
by zero. The same technique gives $K_\nu$ and $Y_\nu$ at integer
$\nu$, and the limits of the ${}_2F_1$ transformations.

## Numerical integration {#integration}

For large real parameters, where both expansions are slow, the
integral representations
$$
{}_1F_1(a,b,z) = \frac{\Gamma(b)}{\Gamma(a)\Gamma(b-a)} \int_0^1
e^{zt} t^{a-1}(1-t)^{b-a-1} dt, \qquad U(a,b,z) =
\frac{1}{\Gamma(a)}\int_0^\infty e^{-zt} t^{a-1} (1+t)^{b-a-1} dt,
$$
and the analogous ones for ${}_2F_1$, $\Gamma(s,z)$, $I_\nu$ and
$K_\nu$ are evaluated with the rigorous integration of
`acb_calc_integrate` [@Joh2018a]. The integrands are smooth and
positive for the parameter ranges in which the functions are
implemented ($a \ge 1$, $b - a \ge 1$ and so on), and the cost of the
integration grows only slowly with the parameters. `acb_hypgeom_m`,
`acb_hypgeom_u`, `acb_hypgeom_2f1`, `acb_hypgeom_gamma_upper`,
`acb_hypgeom_bessel_i` and `acb_hypgeom_bessel_k` call the integration
functions of [`arb_hypgeom`](arb_hypgeom) when their other methods give
less than half of the requested accuracy and all arguments are real
and below $2^{60}$, and keep the better enclosure.

## Use in FLINT {#use}

- Series summation: `acb_hypgeom_pfq_sum` and its variants
  `acb_hypgeom_pfq_sum_forward`, `acb_hypgeom_pfq_sum_rs`,
  `acb_hypgeom_pfq_sum_bs`, `acb_hypgeom_pfq_sum_fme`,
  `acb_hypgeom_pfq_sum_invz`; for power series
  `acb_hypgeom_pfq_series_sum` and its variants; for rational
  parameters `arb_hypgeom_sum_fmpq_arb` and
  `arb_hypgeom_sum_fmpq_imag_arb`.
- Tail bounds: `acb_hypgeom_pfq_bound_factor`,
  `acb_hypgeom_pfq_direct`, `acb_hypgeom_pfq_series_direct`; choice of
  the number of terms: `acb_hypgeom_pfq_choose_n`,
  `acb_hypgeom_pfq_series_choose_n`.
- Asymptotic expansions: `acb_hypgeom_u_asymp` and everything built on
  it (`acb_hypgeom_m_asymp`, `acb_hypgeom_0f1_asymp`,
  `acb_hypgeom_bessel_j_asymp`, `acb_hypgeom_bessel_k_asymp`,
  `acb_hypgeom_erf_asymp`, `acb_hypgeom_gamma_upper_asymp`,
  `acb_hypgeom_ei_asymp`, `acb_hypgeom_si_asymp`, the Coulomb wave
  functions).
- Connection formulas and the choice of method: `acb_hypgeom_u`,
  `acb_hypgeom_u_1f1`, `acb_hypgeom_m`, `acb_hypgeom_0f1`,
  `acb_hypgeom_pfq`.
- The top-level page [](algorithms/special-functions/overview) maps
  every function family to the methods it uses; the generic
  multiple-precision floating-point type [`mp_real`](mp_real) has its
  own hypergeometric binary splitting for constants, see
  [](algorithms/floating-point/overview).
