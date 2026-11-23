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
the lower parameter $b_i$; $C$ is infinite if $p > q$. The bound can
be infinite for a small $n$ even if the series converges; the number
of terms must be large enough that the terms have started to
decrease.

The bound $1 + |a_i - b_i| / |b_i + n|$ for $|a_i + k| / |b_i + k|$ is
far from the truth if $a_i - b_i$ is large while $a_i + k$ and $b_i +
k$ are of similar size, or $a_i + k$ is small. Two sharper bounds for
the term ratio are used where $D \ge 1$.

::: lemma Ratio bounds {#lem-ratio}
Let $a, b$ be complex and $k_0$ an integer with $x_0 =
\operatorname{Re}(b + k_0) > 0$, and let $d = \operatorname{Re}(a -
b)$. For all $k \ge k_0$,
$$
\left|\frac{a + k}{b + k}\right|^2 \le
\max\left(\left(1 + \frac{d}{x_0}\right)^2, 1\right) +
\left(\frac{\operatorname{Im} a}{x_0}\right)^2.
$$
For integers $k_0 \le k < k_1$, $|a + k| \le \max(|a + k_0|, |a + k_1 -
1|)$.
:::

::: proof
With $x = \operatorname{Re}(b + k) \ge x_0$, $|a + k|^2 = (x + d)^2 +
(\operatorname{Im} a)^2$ and $|b + k|^2 \ge x^2$, so the ratio is at
most $(1 + d/x)^2 + (\operatorname{Im} a / x)^2$. The second term
decreases with $x$. If $d \ge 0$, $1 + d/x$ lies between 1 and $1 +
d/x_0$; if $d < 0$, it lies between $1 + d/x_0$ and 1, so its square is
at most the larger of $(1 + d/x_0)^2$ and 1. The second statement holds
because $k \mapsto |a + k|$ is convex.
:::

The first bound replaces the factors of $D$ in the theorem (with
$a_i$ paired with $b_i$ in the given order, or each $a_i$ with the
remaining $b_j$ nearest in real part; any pairing is valid, and the
smallest $D$ is used). It is close to 1 when $a_i + k$ and $b_i + k$
are both large, which is the case for the ${}_2F_1$ series with large
$b$ and $c$. The second bound is used on finite ranges of indices,
described next.

If $D \ge 1$ for the bounds on $[n, \infty)$, or if
$\operatorname{Re}(b_i + n) \le 0$ for some $i$ (then $b_i$ is close to
or beyond a negative integer $-m_i$ with $m_i \ge n$, and the terms can
grow again near $k = m_i$), the range $[n, \infty)$ is split into
subranges $[k_j, k_{j+1})$, the ratio is bounded by $D_j$ on the $j$-th
subrange, and
$$
\sum_{k \ge n} \frac{|T(k)|}{|T(n)|} \le \sum_j P_j S_j, \qquad
P_j = \prod_{l < j} D_l^{k_{l+1} - k_l}, \qquad
S_j = \sum_{t=0}^{k_{j+1} - k_j - 1} D_j^{t},
$$
with the last subrange infinite, taken as soon as all
$\operatorname{Re}(b_i + k_j) > 0$ and the bound on $[k_j, \infty)$ is
less than 1. On a finite subrange, $D_j$ is the smaller of the bound
of the theorem and $|z| \prod_i \sup |a_i + k| / \prod_i \inf |b_i +
k|$, where the supremum is attained at an end point and $|b_i + k|$ is
monotone on the subrange, or is at least the distance from $b_i$ to
the integers if the subrange contains $m_i$. The subranges end at
points refined geometrically on both sides of each $m_i$, and have
length about $\frac{1}{2}\sqrt{k}$ while $D_j > \frac{1}{2}$ (the bounds
for the ratio on a subrange then differ from its values by a factor $1
+ O(k^{-1/2})$, and the product $P_j$ by a bounded factor), doubling
where the ratio is smaller. The computation stops with an infinite
bound when the sum exceeds $2^{64}$.

This allows the series to be truncated before the indices where a
lower parameter passes a nonpositive integer, when the terms are
already small there, and where an upper parameter has a large negative
real part, so that the terms decrease long before the bound of the
theorem does. `acb_hypgeom_pfq_choose_n` takes such a truncation point
only if the terms stay below the tolerance beyond it: near a
nonpositive integer $-b_i$, and after $k = -\operatorname{Re} a_i$,
where the terms can grow to a second maximum. Before a nonpositive
integer $-b_i$ it also requires that the bound $C|T(n)|$ at the
truncation point, computed with `acb_hypgeom_pfq_bound_factor`, is
below the tolerance: $C$ covers the range up to the pole with suprema
of the ratios and can exceed the sum of the terms by many orders of
magnitude, or be infinite. Up to four truncation points are tried ($n$
and 16, 48 and 112 terms more); if the bound is too large at all of
them, the series is summed past the pole. (For ${}_2F_1(2.03125,
-1.203125; -81.8984375; 0.42)$ the terms fall below $2^{-64}$ after 30
terms and stay there, but $C = 2^{52}$ at $n = 30$.) The bound costs
as much as tens of terms when it covers the range up to a pole, and
`acb_hypgeom_pfq_direct` uses the one that was computed for the choice
of $n$. `acb_hypgeom_pfq_direct`
sums $n$ terms with `acb_hypgeom_pfq_sum`, which also returns $T(n)$,
and adds $C|T(n)|$ to the radius. If $T(n)$ is exactly zero, the
series has terminated and no error is added; this is how polynomials
such as ${}_pF_q$ with a nonpositive integer upper parameter are
computed for any $p$, $q$ and $z$.

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
terminating series, steps over parameters close to nonpositive
integers, and looks ahead for a second maximum of the terms when an
upper parameter has negative real part.
`acb_hypgeom_pfq_series_choose_n` does the same from the constant
terms of power series parameters. If an upper parameter is $-m + x +
\ldots$ with a nonnegative integer $m$ (the limit of a parameter, as in
the connection formulas of ${}_2F_1$), the terms beyond $k = m$ have a
multiple of $x$ in place of the factor zero and are of the size of the
other factors; the factor is left out of the simulation, so that the
number of terms is that of the terms which matter. The functions built on the series also estimate the
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
recurrence. With complex parameters, the products of the factors $a_i +
k$ in rectangular splitting are balls that are multiplied by complex
numbers at every step, and their radii grow by a factor of about $1 +
|\operatorname{Im} a_i| / |\operatorname{Re} a_i + k|$ each time, in
total by about $\exp\bigl(\sum_i |\operatorname{Im} a_i| \log((|
\operatorname{Re} a_i| + n)/(|\operatorname{Re} a_i| + 1))\bigr)$. If
this is estimated to exceed about $p/10$ bits at precision $p$, binary splitting
or fast multipoint evaluation is used instead. In an example with $8
\cdot 10^4$ terms and parameters of size $10^5 (1 + 0.3 i)$, the
radius from rectangular splitting is larger by a factor $e^{24000}$. The power series versions
(`acb_hypgeom_pfq_series_sum_forward`, `_bs`, `_rs`) perform the same
algorithms with truncated power series arithmetic.
`acb_hypgeom_pfq_series_sum` uses rectangular splitting above 900 bits
for short parameters and a long argument, binary splitting above 600
bits if both are short (and for series of length above 20), and the
forward recurrence otherwise. The forward recurrence multiplies a ball
by $z$ and by the ratio of the terms at every step; with complex
numbers the radius grows by about 0.16 bits per term (the wrapping
effect), which is all of the precision if there are several times more
terms than bits: the limit cases of the connection formulas of
${}_2F_1$ at $|w| = 0.85$ returned no result below 900 bits. If the
series has more terms than bits (and at least 64), short parameters and
complex data, rectangular splitting is therefore used at any precision
(it multiplies by $z$ only about $\sqrt{n}$ times, at 1.3 to 1.7 times
the cost), or binary splitting if $z$ is short as well. With real
data the forward recurrence does not lose accuracy in this way and is
kept.

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
the target by the size of the exponential factor, up to 10 times the
precision, and of $C_n$), and sums the series
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

**Choosing between them.** `acb_hypgeom_m` and `acb_hypgeom_u` choose
between the convergent series, the asymptotic expansion and numerical
integration ([](#integration)) by estimating what each would cost; see
[](#choice).

**Integers $0 < a < b$.** Both series of the second formula then
terminate, after $b - a$ and $a$ terms, and
$$
{}_1F_1(a,b,z) = \frac{\Gamma(b)}{\Gamma(b-a)} (-z)^{-a}
\sum_{k<b-a} \frac{(a)_k (a-b+1)_k}{k!\, (-z)^k}
+ \frac{\Gamma(b)}{\Gamma(a)}\, e^z z^{a-b}
\sum_{k<a} \frac{(b-a)_k (1-a)_k}{k!\, z^k}
$$
exactly, for any $z \ne 0$. `acb_hypgeom_m` evaluates this finite sum
where it is shorter than the series, which has at least $|z| + p/8$
terms at $p$ bits. For small $|z|$ the terms grow and the two parts
cancel; the working precision is raised by the logarithm of the
largest term (with its factor) minus that of the function, which is
about the larger part for large $|z|$ and at most $e^{\max(0,
\operatorname{Re} z)}$ in size, and the sum is not used if this
exceeds $2p + 100$ bits. ${}_1F_1(10, 20, 1000 (1 + i))$ takes 0.05
milliseconds at 3333 bits instead of 3.2, and ${}_1F_1(2, 5, -10)$ at
64 bits 6 instead of 17 microseconds.

**Limits at integer parameters.** The formula for $U$ has removable
singularities at integer $b$. `acb_hypgeom_u_1f1` handles this case by
evaluating the formula with $b$ replaced by the power series $b + x$
and dividing out the zero of the denominator $\sin(\pi b)$, using the
power series algorithms and tail bounds above; the regularized
functions ${}_1\mathbf{F}_1 = {}_1F_1/\Gamma(b)$ are summed so that
the factor $1/\Gamma(b + k)$ vanishes for the terms that would divide
by zero. The same technique gives $K_\nu$ and $Y_\nu$ at integer
$\nu$, and the limits of the ${}_2F_1$ transformations.

**Parameters close to an integer.** For $b = n + \varepsilon$ with a
small $\varepsilon \ne 0$ the two terms of the formula for $U$ are of
the size $1/\varepsilon$ and cancel. With an exact $b$ this costs
working precision; on a ball for $b$ of radius $r$ the result has at
most $\log_2 (|\varepsilon| / r)$ accurate bits whatever the precision
is, $\log_2 (1/|\varepsilon|)$ fewer than the ball determines, and a
ball which contains $n$ gives nothing. But $U$ is an entire
function of $b$, and
$$
U(a, n + d, z) = (-1)^n \pi\, \frac{A(d) - B(d)}{\sin \pi d}, \qquad
A = \frac{{}_1\mathbf{F}_1(a, b, z)}{\Gamma(a - b + 1)}, \quad
B = \frac{z^{1-b}\, {}_1\mathbf{F}_1(a - b + 1, 2 - b, z)}{\Gamma(a)},
$$
where $A$ and $B$ are entire functions of $d$ which agree at $d = 0$.
The Taylor coefficients $u_j$ of $U$ at $d = 0$ are the coefficients
of the power series $U(a, n + x, z)$ of the previous paragraph. On the
circle $|d| = R = 1/4$ one has $|\sin \pi d| \ge \sin(\pi/4)$, so that
by Cauchy's estimate
$$
|u_j| \le M R^{-j}, \qquad M = \pi \sqrt{2}\, \sup_{|d| \le R} (|A| + |B|),
\qquad
\Bigl| U(a, n + \varepsilon, z) - \sum_{j < J} u_j \varepsilon^j \Bigr|
\le \frac{M q^J}{1 - q}, \quad q = \frac{|\varepsilon|}{R} .
$$
The supremum is bounded from the absolute values of the terms of the
two series ([](#inexact)): with $K$ the number of terms up to the
nonpositive integers in the disc for the lower parameter $\beta$, the
first $K$ terms are bounded one by one with $1/\Gamma(\beta + k) =
(\beta + k) \cdots (\beta + K - 1) / \Gamma(\beta + K)$, and the
others by the bound for ${}_1F_1(\alpha + K, \beta + K, z)$, whose
lower parameter has a real part of at least $2 - R$. `acb_hypgeom_u_1f1`
uses this expansion for $|\varepsilon| < 2^{-12}$, with $J \approx (p +
16) / \log_2 (1/q)$ terms at $p$ bits (and more if $M$, which is
larger than $U$ by the cancellation of $A$ and $B$, turns out to be
what limits the result). The expansion of order $J$ costs about $0.6 +
0.35 J$ times the limit $J = 1$. It is used for a ball which contains
$n$ if $J \le 60$; for a point if $J \le 2$ below 200 bits, $J \le 3$
below 1500 bits and $J \le 6$ beyond (also for a ball which is so
small that the formula keeps the accuracy, $\log_2 (|\varepsilon| /
r) \ge p - 4$); and for another ball which does not contain $n$ if $J
\le 12$, $J \le 4$ and $J \le 6$. For larger orders the formula at the
midpoint with more precision costs less, with the bound for the
variation over the ball of [](#inexact): these are the orders where
the timings cross at 64, 333 and 3333 bits (at 3333 bits, 8.6 against
12 milliseconds for $J = 6$ and 15 against 10 for $J = 8$). A ball
below 200 bits is the exception, because there `acb_hypgeom_u`
accepts the result of the formula on the ball if it has lost less
than 30 bits, which the expansion does not lose (0.35 milliseconds
for $J = 8$ at 64 bits).
For $U(1/3, b, 1)$ at 333 bits, $b = 1 + 10^{-100}$ takes 0.6
milliseconds ($J = 3$), and the ball $b = 1 \pm 10^{-10}$ gives the 34
bits that it determines.

## Numerical integration {#integration}

For large parameters, where both expansions are slow or cancel, the
integral representations
$$
{}_1F_1(a,b,z) = \frac{\Gamma(b)}{\Gamma(a)\Gamma(b-a)} \int_0^1
e^{zt} t^{a-1}(1-t)^{b-a-1} dt, \qquad U(a,b,z) =
\frac{1}{\Gamma(a)}\int_0^\infty e^{-zt} t^{a-1} (1+t)^{b-a-1} dt,
$$
and the analogous one for ${}_2F_1$ are evaluated by rigorous
quadrature. The integrands are $\exp(\varphi(t))$ with $\varphi(t) =
\lambda t + \sum_j \mu_j \log(\alpha_j + \beta_j t)$, for which bounds
on Bernstein ellipses follow from an expansion of $\varphi$ to the
second order. For real parameters up to $2^{40}$ the integral is taken
on the real line, where the integrand is positive
(`arb_hypgeom_1f1_integration`, `arb_hypgeom_u_integration`,
`arb_hypgeom_2f1_integration`); for complex or larger parameters, along
paths of steepest descent through the saddle points of $\varphi$
(`acb_hypgeom_1f1_integration`, `acb_hypgeom_u_integration`,
`acb_hypgeom_2f1_integration`). The cost depends little on the size of
the parameters, but grows quickly with the precision: an integral at
$p$ bits costs about as much as a few thousand terms of a series at
64 bits and a few hundred thousand at 4096 bits, and the quadrature
nodes are computed on the first call at each precision.

The first integral requires $\operatorname{Re} a > 0$ and
$\operatorname{Re}(b - a) > 0$ (the code uses it for real parts of at
least 1). For $\operatorname{Re} a < 1$ and $\operatorname{Re}(b - a)
\ge 1$, which includes the Laguerre-type functions ${}_1F_1(-\nu, b,
z)$ and $U(-\nu, b, z)$ that oscillate for $z > 0$, both functions are
obtained from the integral for $U$ with the first parameter $b - a$ on
two sheets,
$$
I_{\pm} = \Gamma(b - a)\, U(b-a, b, z e^{\pm \pi i}) = \int_0^{\infty e^{i
\theta_{\pm}}} e^{zt} t^{b-a-1} (1+t)^{a-1} \, dt,
$$
over rays that pass below and above the singular point $-1$
respectively (with $|\theta_{\pm} + \arg z \pm \pi| < \pi/2$, which
is possible for $z$ in the right half-plane). The connection formula
[DLMF, 13.2.41](https://dlmf.nist.gov/13.2.E41),
$$
\frac{{}_1F_1(a,b,z)}{\Gamma(b)} = \frac{e^{\mp \pi i a}}{\Gamma(b-a)} U(a,b,z)
+ \frac{e^{\pm \pi i (b-a)}}{\Gamma(a)} e^z U(b-a, b, z e^{\pm \pi i}),
$$
holds with both signs; eliminating $U(a,b,z)$, and the left-hand side
respectively, gives
$$
\frac{{}_1F_1(a,b,z)}{\Gamma(b)} = \frac{e^z\, \Gamma(1-a)}{2 \pi i\, \Gamma(b-a)}
\left(e^{\pi i b} I_+ - e^{-\pi i b} I_-\right), \qquad
U(a,b,z) = \frac{e^z\, \Gamma(1-a)}{2 \pi i}
\left(e^{\pi i (b-a)} I_+ - e^{-\pi i (b-a)} I_-\right)
$$
(using $1/(\Gamma(a) \sin \pi a) = \Gamma(1-a)/\pi$). The two integrals
are complex conjugates for real $a$, $b$ and $z > 0$, and each has the
size of the oscillating function; where the function does not
oscillate (the two saddle points of the integrand are real), they
cancel, and the formulas are not used. Kummer's transformation gives
${}_1F_1$ for $\operatorname{Re} a \ge 1$, $\operatorname{Re}(b - a) <
1$ and $z$ in the left half-plane in the same way, and $U(a,b,z) =
z^{1-b} U(a-b+1, 2-b, z)$ gives $U$ for $\operatorname{Re} a \le 0$.

For ${}_2F_1$, Euler's integral
$$
{}_2F_1(a,b;c;z) = \frac{\Gamma(c)}{\Gamma(b)\Gamma(c-b)} \int_0^1
t^{b-1} (1-t)^{c-b-1} (1 - zt)^{-a} dt
$$
requires $\operatorname{Re} c > \operatorname{Re} b > 0$ (used for
$\operatorname{Re} b \ge 1$, $\operatorname{Re}(c - b) \ge 1$, or the
same for $a$). There is a counterpart of the formulas above where the
integrand is not integrable at $1$: for $a, b$ larger than $c$, as in
${}_2F_1(a, b; c; z)$ with $c - a$ and $c - b$ negative, a function of
Jacobi type which oscillates for $z < 0$. The integrand has zeros at $0$
and at infinity (if $\operatorname{Re}(a - c + 1) > 0$), and the
integrals $E_\pm$ from $0$ to infinity which pass above and below the
point $1$ agree up to $1$; beyond it, their integrands differ by the
factors $e^{\mp \pi i (c - b - 1)}$, so that $e^{\pi i (c-b)} E_+ -
e^{-\pi i (c-b)} E_-$ does not depend on the part beyond $1$, and is
$2 i \sin(\pi (c - b))$ times the integral over $[0, 1]$ where that
exists. With $t \to t/(t-1)$ the two integrals become $e^{\pm \pi i b}
J_\mp$,
$$
J = \int_0^1 t^{b-1} (1 - t)^{a-c} (1 - x t)^{-a} dt, \qquad x = 1 - z,
$$
taken from $0$ to $1$ past the point $1/x$ on the lower side ($J_-$)
and on the upper side ($J_+$), and
$$
{}_2F_1(a,b;c;z) = \frac{\Gamma(c)\, \Gamma(b - c + 1)}{2 \pi i\, \Gamma(b)}
\left(e^{\pi i c} J_- - e^{-\pi i c} J_+\right)
$$
(with $1/(\Gamma(c-b) \sin \pi(c-b)) = \Gamma(b-c+1)/\pi$), an identity
between analytic functions of the parameters which holds for
$\operatorname{Re} b > 0$ and $\operatorname{Re}(a - c) > -1$. For $z <
0$, $J_\mp$ are the Euler integrals of ${}_2F_1(a, b; a+b-c+1; x \mp
i0)$, the values of that function on the two sides of its branch cut
$x > 1$; for other $z$, one of them is its integral over $[0, 1]$ and
the other the continuation across the cut. They are computed along
paths of steepest descent which are homotopic to the polygons from $0$
to $1$ by way of $\tfrac12 \mp i Y$, with $Y$ so large that the
triangle contains $1/x$ where it has to; the homotopy is verified
through the branch cuts of the logarithms as for the other integrals.
For real parameters and $z < 0$ the two integrals are complex
conjugates. They have the size of the function where the saddle points
of the integrand are complex, which is where the function oscillates;
`acb_hypgeom_2f1_integration` uses them for $\operatorname{Re} z < 0$,
$\operatorname{Re} b \ge 1$ and $\operatorname{Re}(a - c) \ge 0$ (or the
same with $a$ and $b$ exchanged) where Euler's integral does not
apply. With parameters $a, b, c = 423, 412, 153$ and $z = -1$, where
the terminating series in $z/(z-1)$ loses about 600 bits, the
integrals give 61 bits at 64 bits of precision in about a millisecond,
and 37 bits if the four arguments are only known to a relative accuracy
of $2^{-53}$, which determines 43 bits of the function.

The same two integrals give the functions of the type of the Jacobi
polynomials, $a < 0 < c < b$, which oscillate for $0 < z < 1$ and
whose Euler integral does not exist at $1$ either. The integral over
the loop from $0$ around $1$ in the positive sense
[DLMF, 15.6.3](https://dlmf.nist.gov/15.6.E3),
$$
{}_2F_1(a,b;c;z) = \frac{\Gamma(c)\, \Gamma(b - c + 1)}{2 \pi i\, \Gamma(b)}
\int_0^{(1+)} t^{b-1} (t - 1)^{c-b-1} (1 - z t)^{-a} dt
$$
($\operatorname{Re} b > 0$, $c - b \ne 1, 2, \ldots$), has an integrand
which vanishes at $1/z$ if $\operatorname{Re} a < 0$. For
$\operatorname{Re}(1/z) > 1$, that is $|z - \tfrac12| < \tfrac12$,
the loop is the path from $0$ to $1/z$ below the point $1$ minus the
path above it, with $(t - 1)^{c-b-1} = e^{\mp \pi i (c - b - 1)} (1 -
t)^{c-b-1}$ on them, and $t = s/z$ turns the two paths into the
integrals $J_\mp$ above with the parameters $A = b - c + 1$, $B = b$,
$C = a + b - c + 1$ and the argument $1 - 1/z$ in the left half-plane:
$$
{}_2F_1(a,b;c;z) = \frac{\Gamma(c)\, \Gamma(A)}{2 \pi i\, \Gamma(b)}\, z^{-b}
\left(e^{\pi i A} J_- - e^{-\pi i A} J_+\right), \qquad
J = \int_0^1 s^{b-1} (1 - s)^{-a} (1 - s/z)^{c-b-1} ds.
$$
This differs from the formula for ${}_2F_1(A, B; C; 1 - 1/z)$ in the
phases and in the factors only, and is used for $\operatorname{Re} a
\le 0$, $\operatorname{Re} b \ge 1$ (or the same with $a$ and $b$
exchanged) where the saddle points are complex. The argument $1 -
1/z$ is computed with as many extra bits as the largest parameter
has, since it is raised to that power. The Legendre polynomial
$P_n(0) = {}_2F_1(-n, n+1; 1; \tfrac12)$ takes 1.3 milliseconds at 64
bits for $n = 10^4$ and 2.1 for $n = 10^6$ (62 bits), 10 and 37
milliseconds at 333 bits; the terms of the terminating series add up
to $P_n(2) \approx (2 + \sqrt 3)^n$ in absolute value, so that it
loses about $1.9 n$ bits.

On the branch cut $z > 1$ the point $1/z$ lies on the segment of
Euler's integral. The value of the function there is the limit from
the lower half-plane, where $1/z$ is above the segment: it is the
integral over a path from $0$ to $1$ which passes below $1/z$,
computed as above for the polygon by way of $\tfrac12 - i Y$ with
the term $1 - z t$ given as $-i (1 - z t)$, whose principal logarithm
is continuous along the polygon (and by way of $\tfrac12 + i Y$ for
$\operatorname{Im} z > 0$, $\operatorname{Re} z > 1$, where the
segment itself passes $1/z$ too closely). With real parameters the
saddle points can both be real, on the two sides of $1/z$, and the
line of steepest descent from the first runs along the real axis into
the second instead of leading to the end point. The path is then
planned for a kernel in which the singular point (if the integrand is
large there, $\operatorname{Re} a > 0$) is moved off the
axis, away from the polygon, by $\min(0.05, 2/|a|)$ times its distance
to the nearer end point: this changes the exponent by a phase of the
order of one, the line passes the second saddle point on the side of
the polygon and turns to the end point. The plan is heuristic in any
case; the path is integrated with the kernel as it is, and the
homotopy is verified for it. This plan is tried first if $1/z$ is on
the segment or within $10^{-6}$ of it (relative to its distance from
the end points), where the other one follows rounding errors, and
second otherwise. ${}_2F_1(1000, 2000; 5000; z)$ at 64 bits
takes about a millisecond both at $z = 5/2$ (complex saddle points)
and at $z = 77/64$ (real ones), with 62 bits. The polygon itself is not
integrated, so that there is no result if the parameters are small
(the series and the connection formulas apply there).

For $\operatorname{Re} z > 1$ there is a second pair of integrals:
Euler's transformation ${}_2F_1(a,b;c;z) = (1 - z)^{c-a-b}\,
{}_2F_1(c-a, c-b; c; z)$ exchanges $b$ and $c - b$, so that its
integral exists for the same parameters, with the exponent $a - c$ at
$1/z$ in place of $-a$. For $a > c$ the first integrand has a pole at
$1/z$ and can have both saddle points on the real axis outside of $[0,
1]$; its path then leads from $0$ back to the saddle point on the
left, around the pole through the lower half-plane and by way of the
saddle point on the right to $1$, and was not found for ${}_2F_1(38191
+ 37/125, 22309; 27150 + 209/250; z)$ at $z = 3373/1000$ and close to
it. The second integrand has a zero at $1/z$, and its integral is the
sum of two integrals over the real segments on the two sides of it.
An integral with a zero at $1/z$ is therefore tried first; that
function takes 1.2 milliseconds at 64 bits (62 bits). On the cut the
factor $(1 - z)^{c-a-b}$ has the argument $\pi (c - a - b)$, the limit
from below.

The Bessel and incomplete gamma functions are integrated along
steepest-descent contours of their own integral representations, and
choose between the series, the asymptotic expansions and the
integration by estimates of their costs
([](algorithms/special-functions/steepest-descent)).

## Choice of algorithm for large parameters {#choice}

`acb_hypgeom_m`, `acb_hypgeom_u` and `acb_hypgeom_2f1` aim at a result
that is accurate to nearly the full precision when the input is exact.
For parameters of moderate size the methods of ${}_1F_1$ and $U$ are
chosen from $z$ alone, as described above (not if that choice is the
asymptotic expansion, the other parameter is large and the terms of
the expansion do not get small enough by the estimate below, as for
${}_1F_1(1, 1000 + 1/2, 1000)$, where the attempt cost eight times as
much as the series at 333 bits); ${}_2F_1$ does so for
parameters up to 30 where a series at $z$ or $z/(z-1)$ converges with
ratio $3/4$ (also with larger parameters if the terms of
the series at $z$ decrease by
a factor of four or more each until they reach the precision, so that
their sum is within a third of the first one: ${}_2F_1(1000 i, 1/3;
1/2 + 2000 i; 1/2)$ takes 20 microseconds at 64 bits instead of 150
with the estimates), and
for parameters up to 12 where a connection formula
has an argument of at most $1/2$ and is not a limit case
([](algorithms/special-functions/gauss-2f1)). Otherwise the series
lose accuracy to cancellation that depends on all the arguments, their
numbers of terms differ by large factors, and the functions estimate,
in double precision from the midpoints of the input:

* *The largest term* of each series, from $\log |(a)_k| =
  \operatorname{Re} \log \Gamma(a + k) - \operatorname{Re} \log
  \Gamma(a)$ (with the reflection formula for negative real parts,
  where the terms almost vanish near $k = -a$ and can be large near $k
  = -b$), and the index where the terms have decreased by the working
  precision. This gives the number of terms $T$.
* *The size of the function*, from the integrand $e^{\varphi}$ of its
  integral representation at a saddle point $t_s$: $\log |f| \approx
  \operatorname{Re} \varphi(t_s) + \frac{1}{2}\log(2\pi) - \frac{1}{2}
  \log |\varphi''(t_s)|$ plus the logarithm of the gamma factors. The
  saddle points are the roots of a quadratic equation. For real
  arguments in the domain of the integral, the function is given by
  the saddle point where the integrand has its maximum on the interval
  of integration; for ${}_1F_1$ with $a < 0$ or $b - a < 0$, by a
  saddle point on the real axis where $\varphi$ has a minimum along
  the axis (the integral over a loop around 0 or 1 crosses the axis
  there). With complex arguments, the saddle point that satisfies the
  same conditions on the real parts is tried first, then the other
  one. Where two saddle points almost coincide (a turning point: $z
  \approx -4a$ for $U$ with a large negative $a$) $\varphi''$ is
  small, and the width of the saddle is that of the cubic term,
  $\frac{2}{3} \Gamma(\frac{1}{3})\, (6 / |\varphi'''|)^{1/3}$, if this is
  less than $(2\pi / |\varphi''|)^{1/2}$.
* *The cancellation* in a series: the logarithm of its largest term
  (with the factors in front of the series) minus that of the
  function. The series is summed with this many additional bits, and
  with a limit on the number of terms that depends on $T$ instead of
  the precision.
* *The asymptotic expansion*: whether the terms of ${}_2F_0$ reach the
  target with the factor $\exp(2 \alpha \rho C_1 / |z|)$ of the error
  bound ([](#asymptotic)), which is large when $|a (a - b)|$ is not
  small compared to $|z|$, and the cancellation if the terms grow
  before they decrease.
* *Costs*: the number of terms times the cost of a term at the
  working precision, and the cost of an integral at the target
  precision, both from timings.

The cheapest method is tried first. If the result is less accurate
than the target, its accuracy shows how much the cancellation was
underestimated, and the same method is repeated with the missing bits
if it is still the cheapest; if the result carries no information, the
estimate from the other saddle point is used, or the precision is
doubled. The attempts stop when their total cost exceeds a few times
that of an integral, so that the functions do not sum millions of
terms for huge parameters where no method applies.

The asymptotic series of $U$ is an exception to the rule that the
accuracy of a result shows what is missing: close to a turning point
the factor of its error bound exceeds the smallest term, whatever the
cancellation is, and the result carries no information. It is then
repeated with twice the precision (up to twice, within a third of the
remaining budget), which reaches smaller terms: $U(-1000.5, 3, 4008)$
takes about a millisecond in this way at any precision up to 3333
bits, where the error bound of the series at 64 bits is $2^{801}$
times the value and the integral takes 20 milliseconds at 128 bits.

For ${}_1F_1$ the candidates are the series and its Kummer
transformation $e^z {}_1F_1(b - a, b, -z)$, which have different
largest terms; for $U$, the two series of the connection formula with
their factors $\Gamma(1 - b)/\Gamma(a - b + 1)$ and $\Gamma(b - 1)
z^{1-b} / \Gamma(a)$ (for an integer $b$, the regularized series: the
one whose lower parameter $b$ or $2 - b$ is a nonpositive integer $-m$
begins at the term $m + 1$, and both are summed beyond it; the limit
costs 3.5 times as much per term for real arguments and 1.4 times for
complex ones, and the expansion of order $J$ at a nearby integer
([](#confluent)) $0.6 + 0.35 J$ times that); for ${}_2F_1$, the four
series
$$
F(a,b,c,z), \quad (1-z)^{c-a-b} F(c-a,c-b,c,z), \quad
(1-z)^{-a} F\!\left(a,c-b,c,\tfrac{z}{z-1}\right), \quad
(1-z)^{-b} F\!\left(c-a,b,c,\tfrac{z}{z-1}\right)
$$
for $|z| < 1$ or $|z/(z-1)| < 1$, the connection formulas and the
corner algorithm, described next.

### Connection formulas of ${}_2F_1$ {#choice-2f1}

The connection formulas at $w = 1/z$, $1/(1-z)$, $1 - z$ and $1 - 1/z$
[@NIST2012, 15.8.2--15.8.5] are each a sum of two series in $w$ with
factors of gamma functions and powers, for example
$$
\frac{{}_2F_1(a,b,c,z)}{\Gamma(c)} =
\frac{\Gamma(b-a)\,(-z)^{-a}}{\Gamma(b)\Gamma(c-a)}\,
\frac{{}_2F_1(a, a-c+1, a-b+1, \tfrac1z)}{\Gamma(a-b+1)}
+ \frac{\Gamma(a-b)\,(-z)^{-b}}{\Gamma(a)\Gamma(c-b)}\,
\frac{{}_2F_1(b, b-c+1, b-a+1, \tfrac1z)}{\Gamma(b-a+1)} .
$$
The choice by the smallest $|w|$ (below $3/4$; the corner algorithm
otherwise) ignores that the parameters of the two series change the
number of terms by large factors, that the two terms can cancel, that
the limit cases cost three times as much per term, and that the series
at $z$ may converge faster than its $|z|$ suggests (${}_2F_1(1, 1000;
2000; 0.99)$ has terms that decrease like $2^{-k}$, and the formula at
$1 - z$ needs 1000 terms for its limit). In an experiment with 750
random inputs outside of $|z|, |z/(z-1)| \le 3/4$ (parameters up to 10
or 50, real, integer or complex; $z$ near the unit circle, near 1,
near $e^{\pm \pi i/3}$ or up to 30), the algorithm chosen by $z$ was
within a factor 1.25 of the fastest formula in half of the cases, and
with the retries that an accurate result required it took 3.6 times as
long in total as the fastest formula at 64 bits, and 1.9 times as long
at 333 bits; the corner algorithm was slower than the fastest formula by
a factor 5.7 in the median.

`acb_hypgeom_2f1` therefore treats as candidates, with estimated
costs: the four series above (for $\log |z|$ or $\log |z/(z-1)|$ below
$-0.01$), the four connection formulas (for $\log |w| < -0.01$), the
formula at $1 - 1/z$ with $a$ and $b$ exchanged, whose series are
different, and the corner algorithm where it would be chosen by $z$.

* *The largest term of a connection formula* is the larger one of
  the two products of the factor and the largest term of the series,
  $\Lambda = \max_j (\log |g_j| + \operatorname{Re}(e_j \log w') +
  \max_k \log |t_{j,k}|)$, with $\log |\Gamma|$ from the Stirling series
  or the reflection formula. The cancellation is $(\Lambda - \log
  |f|) / \log 2$ bits, as for the series. In a limit case ($b - a$ an
  integer for $1/z$ and $1/(1-z)$, $c - a - b$ for $1 - z$ and $1 -
  1/z$) both terms have poles which cancel; what remains has the size
  of the terms with the gamma function at its pole replaced by its
  size away from the poles, and with the Pochhammer symbol of a lower
  parameter $-m$ beyond its pole replaced by $m!\,(k - m - 1)!$. The
  series with the lower parameter $1 - m$ is summed past $k = m$, so
  it has at least $m + 1$ terms. A limit case is a candidate only if
  the difference is an exact integer or the caller's flags say so.
* *The terms* of a series with parameters up to 256 are computed in
  double precision (the squares of their absolute values, with the
  exponent kept apart, a few nanoseconds per term), which gives the
  largest term and the number of terms down to a given size where the
  terms are not monotone. With real parameters, the ratio of the terms
  beyond the negative parameters is $|z| (a + k)(b + k) / ((c + k)(k +
  1))$; this exceeds 1 between the roots of a quadratic polynomial, so
  the largest term is among the first ones or at the larger root,
  without a scan of the $(a + b - c - 1)/(1 - |z|)$ terms up to it.
  The terms of a candidate are counted in the order of a first guess
  of the costs and only as far as the cost of the best candidate so
  far, so that the estimates for a series with $10^5$ terms cost no
  more than those for the series which is used.
* *The size of the function* is at most the sum of the absolute
  values of the terms of any candidate: $\log |f| \le \Lambda + \log
  N$ with $N$ the number of terms of the size of the largest one,
  taken as the index of the largest term plus $1 / (-\log |w|)$ (plus
  a constant). The saddle point estimate is capped by the smallest of
  these bounds. If no saddle point lies between the end points of a
  convergent Euler integral, the integrand can have its maximum at an
  end point; then $f \approx \Gamma(c) \lambda^{-b} / \Gamma(c - b)$ with
  $\lambda = c - b - 1 - a z$ if $\operatorname{Re} \lambda > 0$ (and
  the corresponding expression at $t = 1$) is used if it is larger.
  For ${}_2F_1(1000, 1; 2000; 0.99) = 1.98$ the saddle points are
  beyond $t = 1$, where the integrand is $2^{-4677}$.
* *Costs* in microseconds, from timings of 2500 evaluations each at
  64 and 333 bits and 700 at 3333 bits. Per term: 0.38, 0.38 and 0.86 for real data, 0.55,
  0.55 and 1.45 for complex $z$, and 0.8, 1.0 and 2.65 for complex
  parameters, the same for the series and the connection formulas;
  the limit cases cost 3.2 times as much per term; parameters with
  many bits ($1/3$ rounded at 3333 bits) add a multiplication of that
  size per term (at most three multiplications at the working
  precision, so that the cost of a series grows like $p^2 \log p$ for
  bounded parameters at $p$ bits and remains below that of an
  integral, about $p^{2.7}$, at a high enough precision). A long $z$ with short parameters does not: in the
  rectangular splitting $z$ enters by its powers, and a term costs a
  multiplication by a small number (with $a = 1000$, $b = 10$, $c =
  3000$ at 3333 bits, 0.57 microseconds per term for a short $z$, 1.1
  for a real $z$ of full length and 1.7 for a complex one; counting a
  full multiplication made the function choose a connection formula
  for $z = 9/10$ which takes 22 milliseconds instead of 3.3). Fixed
  costs: 1 to 5, 8 to 14 and 50 to 500 for a
  series; 13 to 16, 45 to 50 and 700 to 1100 for the gamma functions
  and powers of a connection formula (43, 170 and 4100 with complex
  parameters). The corner algorithm costs 2.2 to 2.8 per term of a
  series at $|z| = 0.73$ for parameters below 5, and 5 to 20 for
  parameters between 20 and 60 (the Taylor series of its two
  continuation steps get longer).
* *Guard bits*, from the measured loss beyond the estimated
  cancellation: 12 bits (4 to 13 measured) plus one per 64 bits; 18
  for the limit cases, plus 0.6 bits per unit of the largest parameter
  above 12 (up to 60), or 2 bits per unit if $a$, $b$ and $c$ are
  integers (up to 150); 0.03 bits per term for connection formulas
  with complex parameters (up to 60); $2.5$ bits per unit of the
  largest parameter for the corner algorithm.

The cheapest candidate is evaluated, and the others follow as above if
its result is not accurate. If the algorithm chosen by $z$ was
evaluated first (small parameters) and has lost bits, their number
replaces the estimate of $\log |f|$. In the experiment above, the
function takes 0.73 times as long in total as the fastest single
formula at the precision that formula needs (which is not known in
advance) at 64 bits, and 0.82 times at 333 bits; it is slower than the
fastest formula by a factor 1.1 in the median at 64 bits, 2.5 at the
90th percentile and 4.3 at the 99th (1.0, 1.7 and 2.8 at 333 bits), and
accurate in all cases. The estimates themselves cost about 10
microseconds.

Polynomials (${}_2F_1$ with $a$, $b$, $c - a$ or $c - b$ a nonpositive
integer; $U$ with $a$ or $a - b + 1$ a nonpositive integer) are summed
by formulas of their own, and are not modelled. If the result is not
accurate and the input is exact, the same computation is repeated with
the bits which are found missing, or with twice the precision if the
result carries no information, as long as the estimated cost of the
terms stays within half that of the integral which follows if the
polynomial fails (at least 50 milliseconds; where no integral applies,
16 times the cost of the polynomial at the precision of the caller)
and the precision below $16 p + 8 n$
(plus a constant) for $n$ terms at $p$ bits; the limit is there because
a polynomial which vanishes at $z$ is not accurate at any precision.
(At high precision an integral costs seconds to minutes, and the
polynomial with all the bits it loses is usually much cheaper:
${}_2F_1(-10^4, 10^4 + 1; 1; 1/2)$ at 33333 bits takes 11 milliseconds
by the series at 53000 bits.) An
integral is used first where it applies and costs less than four
evaluations of the polynomial at twice the precision, and in the end
if the polynomial has not become accurate. With inexact input, the
integral on the input is used if it is that cheap and accurate (within
18 bits of what the input allows: the alternatives cost two or more
further integrals); if
only the argument is inexact, the Taylor expansion of the polynomial at
the midpoint of $z$; and otherwise the polynomial at the midpoints,
with a bound for the variation ([](#inexact)). A terminating series whose tail bound is infinite where it was
truncated (a polynomial at $|z| \ge 1$, for example) is summed to its
last term.

## Inexact input {#inexact}

The error bounds of a series for inexact parameters are inflated by
the cancellation: with terms $2^K$ times larger than the sum, $K$ of
the bits of the input are lost, and more precision does not help. The
three functions instead evaluate at the midpoints $m$ of the input
with the algorithm for exact input, and bound the variation of the
function over the input.

*Cauchy bounds.* Let $f$ be analytic on the closed polydisc with the
centre $m$ and radii $R_i$, and $|f| \le M$ there. The Taylor
coefficients at $m$ satisfy $|c_k| \le M / R^k$, so that for $|d_i| \le
r_i$ with $\rho_i = r_i / R_i$ and $s = \sum_i \rho_i < 1$,
$$
|f(m + d) - f(m)| \le M \Bigl(\prod_i \frac{1}{1 - \rho_i} - 1\Bigr) \le \frac{M s}{1 - s},
\qquad
\Bigl|f(m + d) - f(m) - \sum_i \partial_i f(m)\, d_i\Bigr| \le \frac{M s^2}{1 - s} .
$$
The radii are chosen as the scales on which the function varies by a
bounded factor: of the order $1 / \log(|a| + |b| + |z|)$ for the
parameters, and $1 / (1 + |a/b|)$, $|z| / (|z| + |a| + |b|)$ or $|1 - z|
/ (1 + \min(|a|, |b|, |ab/c|))$ for $z$, with discs for $z$ that keep
away from the branch cuts. If $M$ is within a small factor of
$|f(m)|$, the first bound is within a small factor of the true
variation.

*Divided differences.* If $M$ is too large by a factor $2^K$ (a bound
obtained from a series in ball arithmetic, which is too large by the
cancellation), the second bound is used, in which the remainder has
the additional factor $s$: it is negligible if the input is accurate
to more than $K$ bits. The derivatives at $m$ are bounded by divided
differences,
$$
|\partial_i f(m)| \le \frac{|f(m + h e_i) - f(m)|}{h} + \frac{M h}{R_i^2 (1 - h / R_i)},
\qquad 0 < h < R_i ,
$$
(from the Taylor series in the variable $i$), with $h = 2^{-K-4} R_i /
n$. This costs one value of $f$ per inexact variable at a precision of
about $K$ bits, since the difference has to be resolved and not $f$
itself: at high precision, where $K$ is small compared to the
precision, the function on inexact input costs little more than on
exact input, and the result is accurate to within a few bits of what
the input allows.

*Taylor expansion in the argument.* For a polynomial of the kinds of
[](#choice) (also ${}_1F_1$ with $a$ or $b - a$ a nonpositive integer)
with exact parameters, the derivatives with respect to $z$ are
polynomials of the same kind, with $c_j = f^{(j)}/j!$ equal to
$$
\frac{(a)_j (b)_j}{(c)_j\, j!}\, {}_2F_1(a+j, b+j, c+j, z), \qquad
\frac{(a)_j}{(b)_j\, j!}\, {}_1F_1(a+j, b+j, z), \qquad
\frac{(-1)^j (a)_j}{j!}\, U(a+j, b+j, z),
$$
and for $z = m + d$ in a box on which $f$ is analytic,
$$
f(m + d) = \sum_{j < J} c_j(m)\, d^j + R_J, \qquad
|R_J| \le \sup\nolimits_{\text{box}} |c_J| \, |d|^J .
$$
The supremum is bounded by $c_J$ on the box in ball arithmetic at 64
bits, which is too large by the cancellation $2^K$ of the terms; each
order gains about the accuracy $p$ of the input, so that $J \approx 2
+ K/p$ orders bring the remainder below an eighth of the first-order
term $|c_1(m)| r$, which is the variation of $f$ over a box of radius
$r$ and cannot be avoided. $c_0(m)$ and $c_1(m)$ are computed with
exact arguments by the function itself (which repeats a polynomial
with more precision, [](#choice)), and the coefficients from the
second on from the recurrences which the differential equations give,
$$
m (1 - m) (j+2)(j+1)\, c_{j+2} = (j + a)(j + b)\, c_j - (j + 1) \bigl((1 - 2m) j + c - (a + b + 1) m\bigr) c_{j+1}
$$
for the hypergeometric equation and
$$
m (j+2)(j+1)\, c_{j+2} = (j + a)\, c_j - (j + 1)(j + b - m)\, c_{j+1}
$$
for Kummer's, in ball arithmetic (and from a value of the derivative
where a ball has lost its accuracy; only the magnitudes of these
coefficients are used). The cost is that of two values of the
function whatever $K$ is, with orders up to 60. The crude values are
computed first, and nothing else if they show that no order will do
(a box which is not small compared to the scale of the function). The
box must not meet the cut of the factor by which the function differs
from a polynomial in $z$: $(1 - z)^{c-a-b}$ for ${}_2F_1$ with $c - a$
or $c - b$ a nonpositive integer, and $z^{1-b}$ for $U$ with $a - b +
1$ one, unless the exponent is an integer. The result is within a few
bits of what the input allows: the Legendre polynomial
$P_{1000}(1/3)$ as ${}_2F_1(-1000, 1001, 1, z)$ with $z$ rounded at 53
bits, whose terms cancel by 1650 bits, has 42 accurate bits, which is
what the rounded argument determines.

*Bounds from the terms.* The sum of the absolute values of the terms
of a series bounds the function,
$$
|{}_pF_q(a, b, z)| \le \sum_k t_k, \qquad
t_{k+1} = t_k\, \frac{\sup |z| \prod_i \sup |a_i + k|}{(k+1) \prod_i \inf |b_i + k|},
$$
with the suprema over the balls. For $j \ge k$ the ratios are at most
$q_k = |z| \prod_i \max(1, (U_i + k)/(L_i + k))$ with $U_i = |\operatorname{Re}
a_i| + |\operatorname{Im} a_i| + \operatorname{rad}(a_i)$ and $L_i =
\operatorname{Re} b_i - \operatorname{rad}(b_i)$ (and $1/(L_i + k)$ for
a lower parameter without a partner, counting $k!$ as the lower
parameter 1), so that the tail after $t_k$ is at most $t_k q_k / (1 -
q_k)$. The upper parameters can be paired with the lower ones in
either order, and the smaller $q_k$ is used: for ${}_2F_1(1, 1000;
2000; z)$ the pairs $(1000, 2000)$ and $(1, k + 1)$ give $q_k = |z|$
from the first term on, the other pairs only after $10^5$ terms. The
sum stops when $q_k \le 1/2$ or when the tail bound is at most $1/8$ of
the sum (for ${}_2F_1$, where $q_k$ stays above $|z|$), so that the
bound is within a factor $9/8$ of the sum of the $t_k$. The sum for
a further series stops as soon as it exceeds the smallest bound found
so far divided by the factor in front of the series, or after four
times as many terms as the shortest sum so far (for ${}_2F_1(1000,
10; 3000; z)$ on a disc around $9/10$, the series with the parameters
$c - a$ and $c - b$ would take 24000 terms, 1.6 milliseconds, to
arrive at a bound which is of no use, where the first series takes
100). The factors are computed in double precision with outward
rounding and the terms as magnitudes, at 50 to 100 nanoseconds per
term, much less than the value of the function. For ${}_1F_1$ this
is done for the series and for $e^z {}_1F_1(b - a, b, -z)$, and for
${}_2F_1$ for the four series of [](#choice) with their factors. The
bound is close to the function where the terms of one of the series
have equal signs, which is the case for real arguments in the domain
of the Euler integral ($0 < a < b$ for ${}_1F_1$; $a > 0$, $0 < b < c$
and $z < 1$ for ${}_2F_1$), and too large by the cancellation
otherwise: it then takes the place of the series in ball arithmetic.

*Bounds from integrals with positive integrands.* In the domains of
the Euler integrals there are bounds for which $K$ is small. For
$\operatorname{Re} a > 0$ and $\operatorname{Re}(b - a) > 0$, with
$\alpha, \beta, x$ the real parts of $a, b, z$,
$$
|{}_1F_1(a,b,z)| \le G\, {}_1F_1(\alpha, \beta, x), \qquad
G = \frac{|\Gamma(b)|\, \Gamma(\alpha)\, \Gamma(\beta - \alpha)}{|\Gamma(a)|\, |\Gamma(b - a)|\, \Gamma(\beta)},
$$
because $|e^{zt} t^{a-1} (1-t)^{b-a-1}|$ is the integrand at the real
parts. Here $|\Gamma(b)| \le \Gamma(\beta)$, and
$$
\frac{\Gamma(x)}{|\Gamma(x + iy)|} = \prod_{k \ge 0} \Bigl(1 +
\frac{y^2}{(x+k)^2}\Bigr)^{1/2} \le \exp\Bigl(\frac{y^2 \psi'(x)}{2}\Bigr)
\le \exp\Bigl(\frac{y^2}{2} \Bigl(\frac{1}{x} + \frac{1}{x^2}\Bigr)\Bigr).
$$
${}_1F_1(\alpha, \beta, x)$ is the mean of $e^{xT}$ for $T$ with the
beta distribution with parameters $(\alpha, \beta - \alpha)$, which
increases stochastically with $\alpha$ and decreases with $\beta$: the
function increases with $x$, and for $x \ge 0$ it increases with
$\alpha$ and decreases with $\beta$ (the opposite for $x \le 0$). Its
supremum over a box is therefore the value at one corner, computed at
low precision with exact arguments. In the same way, for
$\operatorname{Re} a > 0$ and $\operatorname{Re} z > 0$,
$$
|U(a,b,z)| \le \frac{\Gamma(\alpha)}{|\Gamma(a)|}\, U(\alpha, \beta, x),
$$
where the integral $\Gamma(\alpha) U(\alpha, \beta, x)$ decreases with
$\alpha$ and $x$ and increases with $\beta$ (also for $z^{1-b} U(a - b +
1, 2 - b, z)$ if $\operatorname{Re}(a - b + 1) > 0$), and for
$\operatorname{Re} b > 0$, $\operatorname{Re}(c - b) > 0$ and
$\operatorname{Re} z < 1$, with $\gamma = \operatorname{Re} c$ and $T$
distributed as $\mathrm{Beta}(\beta, \gamma - \beta)$,
$$
|{}_2F_1(a,b,c,z)| \le G\, e^{|\operatorname{Im} a|\, |\arg(1 - z)|}\,
\mathbf{E}\bigl[\,|1 - zT|^{-\alpha}\bigr], \qquad
|1 - zt|^{-\alpha} \le \begin{cases} (1 - xt)^{-\alpha}, & \alpha \ge 0, \\
(1 + |z| t)^{-\alpha}, & \alpha \le 0, \end{cases}
$$
so that the expectation is at most ${}_2F_1(\alpha, \beta, \gamma, x)$
or ${}_2F_1(\alpha, \beta, \gamma, -|z|)$, again monotone in each
parameter, with $G$ as above for $b$ and $c - b$ (and the same with
$a$ and $b$ exchanged). These bounds are close to the function unless
the arguments are far from real.

*Bounds from the integrals on balls.* The integrals along paths of
steepest descent ([](#integration)) evaluated at 30 bits on the
polydisc as balls give enclosures whose size is that of the integral
of the absolute value of the integrand along the path, which is close
to $|f|$ on such a path. This is the bound for complex arguments,
for $\operatorname{Re} a < 1$ through the two integrals for $U$, and
for ${}_2F_1$ with $a, b > c$ and $\operatorname{Re} z < 0$ through
its two integrals. It
costs as much as the integrals do (one or two milliseconds), which is
much more than a series for parameters of moderate size.

These bounds are evaluated in the order of their costs: from the
terms, from the positive integrands (a value of the function at 40
bits), from the asymptotic series or the algorithm chosen by $z$ on
balls, and from the integrals on balls. If the bounds from the series
do not exist or are far too large on the polydisc (small parameters),
they are computed on a smaller one. For $U$ with $b$ at an integer
$n$ or close to one, where the series of the connection formula on a
ball for $b$ that contains $n$ has poles, the disc for $b$ is kept
within $|b - n| \le 1/4$ and the bound is $M$ from the terms of the
limit ([](#confluent)): $U(-1000.5, 3, 500)$ with all three arguments
rounded at 3325 bits has 3310 accurate bits in 19 milliseconds.

*What is attempted.* The costs are estimated from the model of
[](#choice): the cost of a value of $f$ at a low precision $p$ by a
series that loses $c$ bits, compared to the value at the full
precision, is that of a term at $p + c$ bits compared to one at
$\mathit{prec} + c$ bits (the number of terms depends little on the
precision; a value at 40 bits costs a quarter to a half of one at a
few hundred bits, and a tenth of one at a few thousand). With a bound
from a positive integrand, the first-order bound is used, improved by
the divided differences if they add at most half to the cost, which
is the case from about 2000 bits (they remove the factors by
which the $R_i$ are too small, and the loss is then one to three bits
instead of five to ten). Otherwise, the cheaper one of the
divided differences and the bound from the integrals on balls is used.

If the input is as accurate as the working precision, so that a
caller who needs more accuracy can increase both, this is done only if
it costs no more than twice the value of the function itself (less
than the function at twice the precision); otherwise the
result is the first-order bound with the bound from the series, which
loses $K$ bits as the series in ball arithmetic would, and nothing is
computed if the estimated cancellation exceeds the accuracy of the
input and no integral applies at a reasonable cost: at most four
times the estimated cost of the series at the precision which their
cancellation calls for, which is what a caller would pay in the end.
${}_1F_1(-1000 - 1/3, 2000 + 1/2, 1000)$ with the parameters rounded
at 64 bits, whose series cancel by 1644 bits, has 48 of the 55 bits
that the input determines after 0.7 milliseconds from the two
integrals for $U$ (a loop that doubles the precision takes 3.3
milliseconds); at 333 bits the integrals would cost 4.2 times the
series, the result is indeterminate after 10 microseconds, and the
loop takes 2.5 milliseconds. A loop that
doubles the precision until the result is accurate then spends its
time on the series at the final precision and not on integrals at low
precision, and at high precision the divided differences give the
full accuracy at once. If the input is less accurate than the working
precision (more precision cannot help), the bounds are computed
whatever they cost.

Polynomials with exact parameters and an inexact argument are computed
by the Taylor expansion above, whether or not the input is as accurate
as the precision.

Four shortcuts keep the cost at that of the function on exact input.
Where the value itself would be computed by an integral (the Euler
integral, or for $U$ the two integrals of [](#integration) for
$\operatorname{Re} a < 1$), the integral is evaluated on the input
balls (and accepted within 30 bits of what the input allows): $U(-10^4
- 1/3, 1, 10^4)$ with the arguments rounded at 64 bits, whose series
cancel by more bits than the input has, is accurate to 43 bits in 0.6
milliseconds. The same is done with an integral of up to four times
the cost of the series where the series on the input give nothing:
for input as accurate as the precision where nothing would be
computed, as described above, and for input which is less accurate
than the precision when the cancellation exceeds its accuracy, for
which the bounds would be computed whatever they cost. At low precision (up to 100 bits) a series
that loses at most 8 bits, or the asymptotic series, is evaluated on
the input balls (for ${}_2F_1$ also a connection formula or the corner
algorithm, [](#choice-2f1)). For ${}_2F_1$ this is also done up to
1000 bits with a candidate which loses at most 3 bits, such as a
series with positive terms: the value at the midpoints with a bound
costs twice as much or more, and is not more accurate
(${}_2F_1(1000, 10; 3000; 9/10)$ rounded at 333 bits takes 0.18
milliseconds instead of 0.6, and ${}_2F_1(10^6, 1/3; 2 \cdot 10^6 +
1/2; -1)$ 0.26 instead of 2.7). Beyond 1000 bits a result on the
input is a few bits short of what is accepted there, and the work
would be done twice; for ${}_1F_1$ and $U$ nothing is gained. For parameters of moderate size, the
algorithm chosen by $z$ is first applied to the input balls as for
exact input, and the result is accepted within 12 bits of what the
input allows (within 30 bits up to a precision of 1000 bits, where the
bounds cost several times as much as the function). What the input
allows is taken from the relative accuracy of the arguments, and for
${}_2F_1$ with a ball $z$ close to 1 also from that of $1 - z$, on
which the function depends there: with $s = c - a - b$ it behaves like
$(1 - z)^s$ for $\operatorname{Re} s < 0$ and its derivative like
$(1 - z)^{s-1}$ for $\operatorname{Re} s < 1$, so that of the bits by
which $1 - z$ is less accurate than $z$ the part $\min(1, \max(0, 1 -
\operatorname{Re} s))$ is counted as lost (less $\log_2 (1/|s|)$ for a
small $|s|$). ${}_2F_1(1/3, 1/4; 1/5; 1 - 10^{-10})$ rounded at 64
bits has 31 accurate bits, which the connection formula at $1 - z$
gives in 10 microseconds. A candidate which
cancels is evaluated on inexact input with no more precision than the
input has, and once it has been limited by the input the candidates
which cancel as much are not tried. If the ball for $z$ meets the cut
$[1, \infty)$ from both sides, only the algorithm chosen by $z$ is
used. Input which is more accurate than the working
precision by at least 20 bits more than a series is estimated to lose
is treated as exact input: the series is evaluated on the input balls
with the precision that makes up for the cancellation, and is repeated
with the bits which are found missing as long as the input provides
them (likewise the algorithm chosen by $z$, and the polynomials of
[](#choice) up to the accuracy of the input). Input which is more
accurate than the working precision by 30 bits or more is not taken
for the result of rounding at this precision, and the bounds are
computed whatever they cost.

Not covered are ${}_2F_1$ outside of the domains of its integrals
when the cancellation exceeds the accuracy of the input (the two
integrals of [](#integration) cover $a, b > c$ for $\operatorname{Re} z
< 0$ and the Jacobi-type functions with $a < 0 < c < b$ for $|z -
\tfrac12| < \tfrac12$ where the function oscillates, not the
neighbourhoods of the turning points on the other side), the
confluent functions close to their
turning points outside of the oscillatory region, and $U$ with an
inexact $z$ on the branch cut; the result is then a wide ball.

## Profiling {#profiling}

The program `src/acb_hypgeom/profile/p-catalog.c` (built by
`make profile MOD=acb_hypgeom`) times the Bessel functions, the
incomplete gamma functions, ${}_1F_1$, $U$ and ${}_2F_1$ on a catalog
of several hundred inputs and prints the accuracy of each result. The
inputs are written as expressions which are evaluated at the working
precision, and are sorted into classes: exact, rounded and inexact
input (with a fixed radius, or one which follows the precision); real
and complex; integer, half-integer and nearly integer parameters; tiny,
normal, large and huge parameters and arguments; and by hand the
oscillatory and the monotonic regions, turning points, zeros, branch
cuts, polynomials and the limit cases of the connection formulas.
Functions, precisions, classes and single cases are selected on the
command line (`-help` lists the options, `-list` the catalog), and
`-input "2f1(5, 2, 400, 11/10*i)"` runs an input which is not in the
catalog.

With `-alg` the default function is compared with its separate
algorithms (the convergent series, the asymptotic expansion, the
functions without and with integration only, the transformations of
${}_2F_1$, the functions of `arb_hypgeom`), and with loops which call
an algorithm at the precisions $p, 2p, 4p, \ldots$ until the result is
accurate; the summary lists the cases in which an algorithm is as
accurate as the default function and faster, or more accurate. With
`-ref` the results are compared with the accuracy which the inputs
allow, estimated from values of the function at the midpoints of the
input balls and at points on their boundaries, and are checked to
contain these values. All results of a case, from all algorithms and at
all precisions, are checked to overlap. Each measurement runs in a
child process, with a limit on the time and the memory, so that an
algorithm can be forced on input for which it is unsuitable. A time is
that of a call with the caches of FLINT filled, except for a call
which takes more than half a second: this is timed once, as the first
call in a new process, which includes the quadrature nodes of an
integral, unless `-warm` asks for a second call.

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
