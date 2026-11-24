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
bounds in `mag_t` arithmetic; $C$ is infinite if $p > q$, and can be
infinite for a small $n$ even if the series converges, since the terms
must have started to decrease. The factor $1 + |a_i - b_i| / |b_i +
n|$ is far from the truth when $a_i - b_i$ is large while $a_i + k$
and $b_i + k$ are of similar size, and the theorem does not apply at
all while a lower parameter lies beyond a negative integer. Two
sharper bounds for the term ratio are used where $D \ge 1$.

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

The first bound replaces the factors of $D$ (any pairing of upper
with lower parameters is valid, and the best one is used); it is close
to 1 when $a_i + k$ and $b_i + k$ are both large, as for ${}_2F_1$ with
large $b$ and $c$. The second is used piecewise: $[n, \infty)$ is split
into subranges $[k_j, k_{j+1})$ on which the ratio is at most $D_j$,
and
$$
\sum_{k \ge n} \frac{|T(k)|}{|T(n)|} \le \sum_j P_j S_j, \qquad
P_j = \prod_{l < j} D_l^{k_{l+1} - k_l}, \qquad
S_j = \sum_{t=0}^{k_{j+1} - k_j - 1} D_j^{t}.
$$
On a finite subrange $D_j$ is $|z| \prod_i \sup |a_i + k| / \prod_i
\inf |b_i + k|$, where $\inf |b_i + k|$ is at least the distance from
$b_i$ to the integers if the subrange contains a pole $-m_i$; the last
subrange is infinite, with the bound of the theorem. The subranges are
short (of length about $\sqrt k$) where the ratio is close to 1, so
that the bound stays within a bounded factor of the sum of the terms.
This allows a series to be truncated before a lower parameter passes a
nonpositive integer when the terms are already small there, and where
an upper parameter with a large negative real part makes the terms
decrease long before the bound of the theorem does.

`acb_hypgeom_pfq_direct` sums $n$ terms with `acb_hypgeom_pfq_sum`,
which also returns $T(n)$, and adds $C|T(n)|$ to the radius. If $T(n)$
is exactly zero, the series has terminated and no error is added; this
is how polynomials such as ${}_pF_q$ with a nonpositive integer upper
parameter are computed for any $p$, $q$ and $z$, and a terminating
series whose bound is infinite where it was truncated is summed to its
last term. A series which does not terminate and whose bound is
infinite gives an indeterminate result.

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
upper parameter has negative real part. Before a nonpositive integer
$-b_i$ it truncates only where the bound $C|T(n)|$ of
`acb_hypgeom_pfq_bound_factor` is below the tolerance as well: $C$
covers the range up to the pole with suprema of the ratios and can be
many orders of magnitude larger than the terms, so a few truncation
points are tried before the series is summed past the pole.
`acb_hypgeom_pfq_series_choose_n` does the same from the constant terms
of power series parameters; an upper parameter $-m + x + \ldots$ does
not end the series, and the terms beyond $k = m$ are counted without
the factor which vanishes at $x = 0$. The functions built on the series also estimate the
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
recurrence. Complex balls need care because of the wrapping effect. In
rectangular splitting the products of the factors $a_i + k$ are
multiplied by complex numbers at every step, and their radii grow by a
factor of about $1 + |\operatorname{Im} a_i| / |\operatorname{Re} a_i +
k|$ each time; where this is estimated to cost more than about $p/10$
bits at precision $p$, binary splitting or fast multipoint evaluation
is used instead (with parameters of size $10^5 (1 + 0.3i)$ and $8 \cdot
10^4$ terms the difference is a factor $e^{24000}$ in the radius). The
power series versions (`acb_hypgeom_pfq_series_sum_forward`, `_bs`,
`_rs`) perform the same algorithms with truncated power series
arithmetic. Their forward recurrence loses about 0.16 bits per term in
the same way with complex data, which is all of the precision for a
series with several times more terms than bits, so
`acb_hypgeom_pfq_series_sum` uses rectangular splitting (or binary
splitting for a short argument) for such series, and otherwise
chooses by the precision and the lengths of the parameters and the
argument.

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
where it is shorter than the convergent series, with the working
precision raised by the cancellation between the two parts, which
occurs for small $|z|$.

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
working precision, and a ball for $b$ of radius $r$ gives at most
$\log_2 (|\varepsilon| / r)$ accurate bits whatever the precision, and
nothing if it contains $n$. But $U$ is an entire function of $b$:
$$
U(a, n + d, z) = (-1)^n \pi\, \frac{A(d) - B(d)}{\sin \pi d}, \qquad
A = \frac{{}_1\mathbf{F}_1(a, b, z)}{\Gamma(a - b + 1)}, \quad
B = \frac{z^{1-b}\, {}_1\mathbf{F}_1(a - b + 1, 2 - b, z)}{\Gamma(a)},
$$
with $A$ and $B$ entire in $d$ and equal at $d = 0$. The Taylor
coefficients $u_j$ of $U$ at $d = 0$ are the coefficients of the power
series of the previous paragraph, and since $|\sin \pi d| \ge
\sin(\pi/4)$ on $|d| = R = 1/4$, Cauchy's estimate gives
$$
\Bigl| U(a, n + \varepsilon, z) - \sum_{j < J} u_j \varepsilon^j \Bigr|
\le \frac{M q^J}{1 - q}, \qquad q = \frac{|\varepsilon|}{R}, \qquad
M = \pi \sqrt{2}\, \sup_{|d| \le R} (|A| + |B|),
$$
with the supremum bounded from the absolute values of the terms of the
two series ([](#inexact)). `acb_hypgeom_u_1f1` uses this expansion for
$|\varepsilon| < 2^{-12}$ where it costs less than the formula with
more precision: for a ball which contains $n$ always (up to the order
60), and otherwise for the first few orders. A ball then gives the
bits that it determines: $U(1/3, 1 \pm 10^{-10}, 1)$ has 34 bits at any
precision.

## Numerical integration {#integration}

For large parameters, where both expansions are slow or cancel, the
integral representations
$$
{}_1F_1(a,b,z) = \frac{\Gamma(b)}{\Gamma(a)\Gamma(b-a)} \int_0^1
e^{zt} t^{a-1}(1-t)^{b-a-1} dt,
$$
$$
U(a,b,z) = \frac{1}{\Gamma(a)}\int_0^\infty e^{-zt} t^{a-1} (1+t)^{b-a-1} dt,
$$
$$
{}_2F_1(a,b;c;z) = \frac{\Gamma(c)}{\Gamma(b)\Gamma(c-b)} \int_0^1
t^{b-1} (1-t)^{c-b-1} (1 - zt)^{-a} dt
$$
are evaluated by rigorous quadrature. They require $\operatorname{Re} a
> 0$ and $\operatorname{Re}(b - a) > 0$, $\operatorname{Re} a > 0$, and
$\operatorname{Re} c > \operatorname{Re} b > 0$ (or the same with $a$
and $b$ exchanged); the code uses them for real parts of at least 1.
Their cost depends little on the size of the parameters but grows
quickly with the precision, and the quadrature nodes are computed (and
cached) on the first call at each precision; the choice of method
([](#choice)) weighs this. All three integrands are $e^{\varphi(t)}$
with
$$
\varphi(t) = \lambda t + \sum_j \mu_j \log(\alpha_j + \beta_j t),
$$
and one implementation serves them all.

**On the real line.** For real parameters up to $2^{40}$ the integrand
is positive on the interval of integration, and the integral is
computed there (`arb_hypgeom_1f1_integration`,
`arb_hypgeom_u_integration`, `arb_hypgeom_2f1_integration`). The
interval is divided into cells, starting from the part around the peak
of the integrand outside which it is negligible. On a cell with
midpoint $m$,
$$
\varphi(m + u) - \varphi(m) = \varphi'(m)\, u + \sum_j \mu_j
\Bigl(\log\bigl(1 + \tfrac{u}{d_j}\bigr) - \tfrac{u}{d_j}\Bigr),
\qquad d_j = m + \frac{\alpha_j}{\beta_j},
$$
in which the large terms of $\varphi$ have cancelled. This bounds the
integrand on Bernstein ellipses around the cell, which gives the number
of nodes of a Gauss–Legendre rule and the bound for its error
([](algorithms/ball-arithmetic/integration#thm-gl)); a cell is bisected
if its halves need fewer nodes together, and the cells where the bound
is below the tolerance are bounded instead of integrated. At an end
point where the integrand is $t^\mu$ times an analytic function, the
first cell is integrated term by term with the power series of that
function, so that the singularity costs no subdivision.

**Along paths of steepest descent.** For complex or larger parameters
the integral is taken along a path through the saddle points of
$\varphi$, on which the integrand does not oscillate
(`acb_hypgeom_1f1_integration`, `acb_hypgeom_u_integration`,
`acb_hypgeom_2f1_integration`). The path is planned in double
precision: the saddle points are the roots of a quadratic; the lines of
steepest descent from each are followed until they reach an end point,
a zero of the integrand or a valley at infinity; and the path is a
chain of such lines from one end of the integral to the other, which
leaves a line for the straight chord to the zero it leads to as soon as
the exponent descends along that chord. Correctness does not depend on
this plan. Each logarithm in $\varphi$ is given a branch cut, $\log w =
\operatorname{Log}(c w) + \gamma$ with a Gaussian integer $c$, which
avoids both the path and the original contour; the integrand is then
analytic in a simply connected domain which contains both, and the two
integrals agree. That the cuts avoid the path is verified on every
cell of the quadrature, and that they avoid the contour separately,
both in ball arithmetic. A poor path costs accuracy or time, not
correctness. The quadrature is the cell scheme above on the chords of
the path, with the exponent expanded to the second order on each chord
to bound it on the ellipse.

With parameters of the size $N$ a saddle point has the width
$N^{-1/2}$, which double precision does not resolve beyond about
$2^{40}$, and $\varphi$, a sum of terms of the size $N$, cannot be
evaluated there to an absolute accuracy of 1. Each saddle point then
gets a *frame*: the kernel is translated to an exact base point $t_0$
close to the root and written as
$$
\varphi(t_0 + \delta) - \varphi(t_0) = \varphi'(t_0)\, \delta + \sum_j \mu_j
\bigl(\log(1 + \zeta_j) - \zeta_j\bigr), \qquad \zeta_j = \frac{\beta_j \delta}{\alpha_j + \beta_j t_0},
$$
with $\varphi'(t_0)$ and the translated coefficients computed once in
ball arithmetic. The path is planned and integrated in these
coordinates, so that the working precision of the quadrature does not
grow with the parameters (which have been tested up to $10^{100}$).
In both cases the integrand is evaluated in vectorized
[`dfloat`](dfloat) arithmetic up to about 200 bits and in `arb`
arithmetic beyond, with the precision lowered at the nodes where the
integrand is small.

**Two integrals past a singular point.** Where the integrals above do
not exist, combinations of two integrals of the same form give the
functions in the regions where they oscillate. For $\operatorname{Re}
a < 1 \le \operatorname{Re}(b - a)$, which includes the Laguerre-type
functions ${}_1F_1(-\nu, b, z)$ and $U(-\nu, b, z)$ that oscillate for
$z > 0$, let
$$
I_{\pm} = \Gamma(b - a)\, U(b-a, b, z e^{\pm \pi i}) = \int_0^{\infty e^{i
\theta_{\pm}}} e^{zt} t^{b-a-1} (1+t)^{a-1} \, dt
$$
over rays which pass below and above the point $-1$ (for $z$ in the
right half-plane). The connection formula
[DLMF, 13.2.41](https://dlmf.nist.gov/13.2.E41) with both signs gives
$$
\frac{{}_1F_1(a,b,z)}{\Gamma(b)} = \frac{e^z\, \Gamma(1-a)}{2 \pi i\, \Gamma(b-a)}
\left(e^{\pi i b} I_+ - e^{-\pi i b} I_-\right), \qquad
U(a,b,z) = \frac{e^z\, \Gamma(1-a)}{2 \pi i}
\left(e^{\pi i (b-a)} I_+ - e^{-\pi i (b-a)} I_-\right).
$$
The two integrals have the size of the function where its saddle points
are complex, which is where it oscillates; where they are real, the
integrals cancel and the formulas are not used. Kummer's
transformation gives ${}_1F_1$ for $\operatorname{Re}(b - a) < 1 \le
\operatorname{Re} a$ in the left half-plane in the same way, and $U(a,
b, z) = z^{1-b} U(a - b + 1, 2 - b, z)$ gives $U$ for $\operatorname{Re}
a \le 0$.

For ${}_2F_1$ with $a$, $b$ larger than $c$, a function of Jacobi type
which oscillates for $z < 0$, Euler's integrand is not integrable at
1. The integrals from $0$ to infinity which pass above and below the
point 1 differ beyond 1 by the factors $e^{\mp \pi i (c - b - 1)}$, so
that a combination of them is $2i \sin(\pi(c - b))$ times the integral
over $[0, 1]$ where that exists; with $t \to t/(t - 1)$ they become
$$
J_\mp = \int_0^1 t^{b-1} (1 - t)^{a-c} (1 - x t)^{-a} dt, \qquad x = 1 - z,
$$
taken past the point $1/x$ on the lower and the upper side, and
$$
{}_2F_1(a,b;c;z) = \frac{\Gamma(c)\, \Gamma(b - c + 1)}{2 \pi i\, \Gamma(b)}
\left(e^{\pi i c} J_- - e^{-\pi i c} J_+\right),
$$
an identity between analytic functions of the parameters for
$\operatorname{Re} b > 0$ and $\operatorname{Re}(a - c) > -1$. The two
paths are planned to be homotopic to the polygons from 0 to 1 by way of
$\tfrac12 \mp iY$, which is verified through the branch cuts as above.
They are used for $\operatorname{Re} z < 0$, $\operatorname{Re} b \ge
1$, $\operatorname{Re}(a - c) \ge 0$ (or with $a$ and $b$ exchanged):
for ${}_2F_1(423, 412; 153; -1)$, where the terminating series in
$z/(z-1)$ loses about 600 bits, they give 61 bits at 64 bits in about
a millisecond.

The functions of the type of the Jacobi polynomials, $a < 0 < c < b$,
oscillate for $0 < z < 1$, and their Euler integral does not exist at
1 either. The integral over the loop from 0 around 1
[DLMF, 15.6.3](https://dlmf.nist.gov/15.6.E3),
$$
{}_2F_1(a,b;c;z) = \frac{\Gamma(c)\, \Gamma(b - c + 1)}{2 \pi i\, \Gamma(b)}
\int_0^{(1+)} t^{b-1} (t - 1)^{c-b-1} (1 - z t)^{-a} dt,
$$
has an integrand which vanishes at $1/z$ if $\operatorname{Re} a < 0$.
For $|z - \tfrac12| < \tfrac12$ the loop is the path from 0 to $1/z$
below the point 1 minus the path above it, and $t = s/z$ turns these
into the integrals $J_\mp$ with the parameters $b - c + 1$, $b$, $a + b
- c + 1$ at $1 - 1/z$, which differs from the formula for that function
only in its phases and factors. This is used for $\operatorname{Re} a
\le 0$, $\operatorname{Re} b \ge 1$ where the saddle points are complex.
The Legendre polynomial $P_n(0) = {}_2F_1(-n, n+1; 1; \tfrac12)$ with
$n = 10^6$ takes 2 milliseconds at 64 bits, where the terminating
series loses about $1.9 n$ bits.

On the branch cut $z > 1$ the value of ${}_2F_1$ is the limit from the
lower half-plane, the integral over a path from 0 to 1 which passes
below $1/z$ (by way of $\tfrac12 - iY$). With two real saddle points on
the two sides of $1/z$, the line of descent from the first runs into
the second instead of the end point; the path is then planned for a
kernel in which the singular point is moved slightly off the axis, and
integrated with the kernel as it is. For $\operatorname{Re} z > 1$ the
integrals of $(1 - z)^{c-a-b}\, {}_2F_1(c-a, c-b; c; z)$, which exist
for the same parameters, are also candidates, and are tried first when
their integrand has a zero at $1/z$ rather than a pole. The polygon
itself is not integrated, so that there is no result for some small
parameters, where the series and connection formulas apply.

The Bessel and incomplete gamma functions are integrated along
steepest-descent contours of their own integral representations
([](algorithms/special-functions/steepest-descent)).

## Choice of algorithm for large parameters {#choice}

`acb_hypgeom_m`, `acb_hypgeom_u` and `acb_hypgeom_2f1` aim at a result
which is accurate to nearly the full precision when the input is exact.
For parameters of moderate size the classical choice from $z$ alone
(above, and [](algorithms/special-functions/gauss-2f1)) is made first,
since it is cheap and usually right. Otherwise the cancellation of the
series depends on all the arguments and their numbers of terms differ
by large factors, and the functions estimate, in double precision from
the midpoints of the input:

* *the largest term* of each candidate series and its number of terms,
  from $\log |(a)_k| = \operatorname{Re} \log \Gamma(a + k) -
  \operatorname{Re} \log \Gamma(a)$, with the reflection formula for
  negative real parts;
* *the size of the function*, from the integrand of its integral
  representation at a saddle point $t_s$, $\log |f| \approx
  \operatorname{Re} \varphi(t_s) + \frac12 \log(2\pi /
  |\varphi''(t_s)|)$ plus the logarithms of the gamma factors (near a
  turning point, where two saddle points merge, with the width of the
  cubic term);
* *the cancellation* of each series, its largest term (with its
  factor) minus the size of the function, which is the number of guard
  bits it needs;
* whether *the asymptotic expansion* reaches the target, from the
  factor $\exp(2 \alpha \rho C_1 / |z|)$ of its error bound
  ([](#asymptotic));
* *the costs*: the number of terms times the cost of a term at the
  working precision, and the cost of an integral, both fitted to
  timings.

The cheapest method is tried first, at the precision which its
estimate calls for. Every result is an enclosure, so the estimates
affect only the time: a result which falls short shows by how much
the cancellation was underestimated, and the method is repeated with
the missing bits if it is still the cheapest; a result without any
information leads to the estimate from another saddle point, or to a
doubled precision. The attempts are limited by a budget of a few times
the cost of an integral, so that the functions do not sum millions of
terms where no method applies. (The asymptotic series of $U$ near a
turning point is limited by its error bound rather than by
cancellation, and is repeated with twice the precision instead.)

The cost of a term with long parameters is bounded by a few
multiplications at the working precision, so that for bounded
parameters the cost of a series grows like $p^2 \log p$ at $p$ bits,
against about $p^{2.7}$ for an integral: at a high enough precision the
series are chosen.

The candidates are, for ${}_1F_1$, the series and its Kummer
transformation; for $U$, the two series of the connection formula (the
regularized ones for an integer $b$, and the Taylor expansion of
[](#confluent) near an integer); and for ${}_2F_1$, the four series
$$
F(a,b,c,z), \quad (1-z)^{c-a-b} F(c-a,c-b,c,z), \quad
(1-z)^{-a} F\!\left(a,c-b,c,\tfrac{z}{z-1}\right), \quad
(1-z)^{-b} F\!\left(c-a,b,c,\tfrac{z}{z-1}\right)
$$
for $|z| < 1$ or $|z/(z-1)| < 1$, the connection formulas and the
corner algorithm.

Polynomials (${}_2F_1$ with $a$, $b$, $c - a$ or $c - b$ a nonpositive
integer, $U$ with $a$ or $a - b + 1$ one) have formulas of their own
and are not modelled. With exact input, a result which is not accurate
is repeated with the missing bits, within a budget relative to the
cost of the integral which would follow (a polynomial which vanishes
at $z$ is not accurate at any precision); the integral is used first
where it is much cheaper, and last otherwise. Inexact input is treated
in [](#inexact).

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
The choice by the smallest $|w|$ ignores that the parameters change the
numbers of terms by large factors, that the two terms can cancel, that
the limit cases cost several times as much per term, and that the
series at $z$ can converge much faster than $|z|$ suggests
(${}_2F_1(1, 1000; 2000; 0.99)$ has terms which decrease like $2^{-k}$,
while the formula at $1 - z$ needs 1000 terms for its limit).
`acb_hypgeom_2f1` therefore estimates for each formula:

* *its largest term*, the larger of the products of the factor and the
  largest term of each of its two series. In a limit case the poles of
  the two terms cancel, and what remains has the size of the terms with
  the gamma function at its pole replaced by its size away from it;
* *the terms themselves*, in double precision for parameters up to
  256, in the order of a first guess of the costs and only as far as
  the best candidate so far, so that the estimate for a series with
  $10^5$ terms costs no more than the series which is used;
* *the size of the function*, capped by the sum of the absolute values
  of the terms of any candidate, and raised to the end-point value
  where the integrand of Euler's integral has its maximum at an end
  point (for ${}_2F_1(1000, 1; 2000; 0.99)$ the saddle points are
  beyond $t = 1$);
* *costs and guard bits* fitted to timings, with extra guard bits for
  the limit cases which grow with the parameters.

The cheapest candidate is evaluated, and the others follow as above if
its result is not accurate. On 750 random inputs outside $|z|,
|z/(z-1)| \le 3/4$, the choice by $z$ alone took 3.6 times as long as
the fastest formula at 64 bits, including the retries which an accurate
result needed, and this choice 0.73 times: less than any single
formula, which has to be repeated with the bits it lacks. The estimates
cost about 10 microseconds.

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
|f(m + d) - f(m)| \le \frac{M s}{1 - s},
\qquad
\Bigl|f(m + d) - f(m) - \sum_i \partial_i f(m)\, d_i\Bigr| \le \frac{M s^2}{1 - s} .
$$
The radii $R_i$ are the scales on which the function varies by a
bounded factor (of the order $1 / \log(|a| + |b| + |z|)$ for the
parameters), with discs for $z$ which keep away from the branch cuts.
If $M$ is within a small factor of $|f(m)|$, the first bound is within
a small factor of the true variation. If $M$ is too large by a factor
$2^K$, the second bound is used, whose remainder is negligible if the
input is accurate to more than $K$ bits, with the derivatives bounded
by divided differences,
$$
|\partial_i f(m)| \le \frac{|f(m + h e_i) - f(m)|}{h} + \frac{M h}{R_i^2 (1 - h / R_i)},
\qquad 0 < h < R_i .
$$
This costs one value of $f$ per inexact variable at a precision of
about $K$ bits, since only the difference has to be resolved: at high
precision the function on inexact input costs little more than on
exact input.

*The bounds $M$*, in the order of their costs:

- the sum of the absolute values of the terms of a series, computed in
  double and `mag` arithmetic with the suprema of the parameters over
  the balls and a bound for the tail (50 to 100 nanoseconds per term).
  It is close to the function where the terms have equal signs, which
  is the case for real arguments in the domain of the Euler integral,
  and too large by the cancellation otherwise;
- in the domains of the Euler integrals, the integrals at the real
  parts $\alpha, \beta, \gamma, x$ of $a, b, c, z$, whose integrands
  dominate the absolute values of the complex ones:
  $$
  |{}_1F_1(a,b,z)| \le G\, {}_1F_1(\alpha, \beta, x), \qquad
  |U(a,b,z)| \le \frac{\Gamma(\alpha)}{|\Gamma(a)|}\, U(\alpha, \beta, x),
  $$
  with $G = |\Gamma(b)|\, \Gamma(\alpha)\, \Gamma(\beta - \alpha) /
  (|\Gamma(a)|\, |\Gamma(b - a)|\, \Gamma(\beta))$, and for ${}_2F_1$ the
  analogous bound by ${}_2F_1(\alpha, \beta, \gamma, x)$ or
  ${}_2F_1(\alpha, \beta, \gamma, -|z|)$ with a factor $e^{|\operatorname{Im}
  a|\, |\arg(1 - z)|}$. These real functions are means of $e^{xT}$ or $(1 -
  xT)^{-\alpha}$ over beta distributions, hence monotone in each
  parameter, so their suprema over a box are values at one corner,
  computed at low precision;
- the asymptotic series or the algorithm chosen by $z$ on balls;
- the integrals along paths of steepest descent ([](#integration))
  evaluated on the balls at 30 bits, whose radius is the integral of the
  absolute value of the integrand along the path, close to $|f|$. This
  is the bound for complex arguments and for the oscillating functions
  of the two-integral formulas.

*Taylor expansion in the argument.* For a polynomial with exact
parameters and an inexact $z$, the derivatives with respect to $z$ are
polynomials of the same kind, $c_j = f^{(j)}/j!$ being for example
$(a)_j (b)_j / ((c)_j\, j!)\; {}_2F_1(a+j, b+j, c+j, z)$, and
$$
f(m + d) = \sum_{j < J} c_j(m)\, d^j + R_J, \qquad
|R_J| \le \sup\nolimits_{\text{box}} |c_J| \, |d|^J .
$$
The supremum, computed in ball arithmetic at 64 bits, is too large by
the cancellation $2^K$, but each order gains about the accuracy $p$ of
the input, so that $J \approx 2 + K/p$ orders bring the remainder below
the first-order term, which is the unavoidable variation over the box.
$c_0(m)$ and $c_1(m)$ are computed by the function itself with exact
arguments, and the higher coefficients from the recurrences which the
hypergeometric and Kummer equations give (only their magnitudes are
needed). The cost is that of two values of the function whatever $K$
is. The Legendre polynomial $P_{1000}(1/3)$ as ${}_2F_1(-1000, 1001, 1,
z)$ with $z$ rounded at 53 bits, whose terms cancel by 1650 bits, has 42
accurate bits, which is what the rounded argument determines.

*Policy.* If the input is less accurate than the working precision,
more precision cannot help, and the bounds are computed whatever they
cost. If it is as accurate as the working precision (as a rule the
result of rounding), a caller who needs more accuracy can increase
both, and the bounds are computed only if they cost no more than about
twice the function itself; otherwise the result loses the cancellation
as the series on balls would, or is indeterminate if the cancellation
exceeds the accuracy of the input and no integral applies at a cost in
proportion. A loop which doubles the precision then spends its time on
the series at the final precision and not on integrals at low
precision. Input which is more accurate than the working precision by
more bits than a series loses is treated as exact input. Cheap
shortcuts come first: where the value would be computed by an integral
anyway, the integral is evaluated on the input balls; at low precision
a series which loses only a few bits is evaluated on the balls; and for
${}_2F_1$ with a ball $z$ close to 1 the accuracy which the input allows
is taken from that of $1 - z$, on which the function then depends.

*Not covered* are ${}_2F_1$ outside of the domains of its integrals
when the cancellation exceeds the accuracy of the input, the confluent
functions close to their turning points outside of the oscillatory
region, and $U$ with an inexact $z$ on the branch cut; the result is
then a wide ball.

## Profiling {#profiling}

The program `src/acb_hypgeom/profile/p-catalog.c` (built by `make
profile MOD=acb_hypgeom`) times the Bessel functions, the incomplete
gamma functions, ${}_1F_1$, $U$ and ${}_2F_1$ on a catalog of several
hundred inputs, sorted into classes (exact, rounded and inexact input;
real and complex; small to huge parameters; turning points, zeros,
branch cuts, polynomials, limit cases), and prints the accuracy of each
result. With `-alg` the default functions are compared with their
separate algorithms and with loops which double the precision until the
result is accurate, and with `-ref` with the accuracy which the input
allows; `-help` lists the options. Each measurement runs in a child
process with limits on time and memory.

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
