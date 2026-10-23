---
title: The gamma function
summary: Rigorous evaluation of Γ, 1/Γ and log Γ in ball arithmetic — the Stirling series with Olver's error bound and argument reduction, rising factorials by rectangular and binary splitting, a tabulated Taylor series of 1/Γ for low and medium precision, fast formulas for rational arguments, and the principal branch of log Γ.
complexity: Stirling series at p bits — O(p) terms and an argument shift by O(p), with O(p^{1/2}) full multiplications for each part using rectangular splitting; rational arguments with small denominators — O(M(p) log² p) by binary splitting.
order: 30
---

FLINT computes $\Gamma(z)$, $1/\Gamma(z)$ and
$\log\Gamma(z)$ for real and complex balls with a combination of four
methods:

| Method | Arguments | Implemented by |
|---|---|---|
| [Stirling series](#stirling) with argument reduction | general; large arguments; high precision | `arb_hypgeom_gamma_stirling`, `acb_hypgeom_gamma_stirling` |
| [Taylor series of $1/\Gamma$](#taylor) with tabulated coefficients | $\lvert z\rvert \lesssim 160$ near the real axis, up to about 3400 bits | `arb_hypgeom_gamma_taylor`, `acb_hypgeom_gamma_taylor` |
| [Rational arguments](#rational) | exact $p/q$ | `arb_hypgeom_gamma_fmpq` |
| [Rising factorials](#rising) | argument reduction $\Gamma(z) = \Gamma(z+r)/(z)_r$ | `arb_hypgeom_rising_ui_rec`, `acb_hypgeom_rising_ui_rec` |

and computes the [principal branch of $\log\Gamma$](#lgamma) with a
branch correction. The methods are described by Johansson [@Joh2023g];
the Stirling series code goes back to Arb [@Joh2017; @Joh2014a].
[Choosing a method](#choosing) summarizes which one is used when.

## The Stirling series {#stirling}

For $z \notin (-\infty, 0]$ and $n \ge 1$,
$$
\log \Gamma(z) = \left(z-\frac{1}{2}\right)\log z - z +
\frac{\log 2 \pi}{2} + \sum_{k=1}^{n-1}  \frac{B_{2k}}{2k(2k-1)z^{2k-1}}
+ R_n(z),
$$
where $B_{2k}$ are Bernoulli numbers and the remainder is exactly
[@Olv1997, pp. 293–295]
$$
R_n(z) = \int_0^{\infty} \frac{B_{2n} - \tilde B_{2n}(x)}{2n(x+z)^{2n}} dx,
$$
with $\tilde B_{2n}(x)$ the periodic Bernoulli function. To evaluate
the gamma function of a power series (as `arb_poly_gamma_series` and
`acb_poly_gamma_series` do), $z$ is replaced by $z + t \in \mathbb{C}[[t]]$.

::: lemma {#lem-olver}
Let $z = |z| e^{i\theta}$ with $|\theta| < \pi$ and $x \ge 0$. Then
$|x + z| \ge (x + |z|) \cos(\theta/2)$.
:::

::: proof
$|x+z|^2 - (x+|z|)^2 \cos^2(\theta/2) = x^2 + 2x|z|\cos\theta + |z|^2 -
(x+|z|)^2 (1 + \cos\theta)/2 = (1 - \cos\theta)(x - |z|)^2/2 \ge 0$.
:::

::: theorem Error bound for the Stirling series {#thm-stirling}
Let $b = 1/\cos(\arg(z)/2)$. For $n \ge 1$ and $k \ge 0$, the
coefficients of the remainder of the Stirling series at $z + t$
satisfy
$$
|[t^k] R_n(z+t)| \le 2 |B_{2n}|
\frac{\Gamma(2n+k-1)}{\Gamma(k+1) \Gamma(2n+1)}
\; |z| \; \left(\frac{b}{|z|}\right)^{2n+k}.
$$
:::

::: proof
The numerator of the integrand satisfies $|B_{2n} - \tilde B_{2n}(x)|
\le 2|B_{2n}|$, since $|\tilde B_{2n}(x)| \le |B_{2n}|$. Expanding
$(x+z+t)^{-2n}$ in powers of $t$ gives $|[t^k](x+z+t)^{-2n}| =
\frac{(2n)_k}{k!} |x+z|^{-2n-k}$, and by [](#lem-olver) $|x+z|^{-1} \le
b/(x+|z|)$. Hence
$$
|[t^k] R_n(z+t)| \le \frac{2|B_{2n}|}{2n} \frac{(2n)_k}{k!} b^{2n+k}
\int_0^\infty \frac{dx}{(x+|z|)^{2n+k}} = \frac{2|B_{2n}|}{2n}
\frac{(2n)_k}{k!} \frac{b^{2n+k}}{(2n+k-1)|z|^{2n+k-1}},
$$
and $\frac{(2n)_k}{2n (2n+k-1) k!} = \frac{\Gamma(2n+k-1)}{\Gamma(2n+1)
\Gamma(k+1)}$.
:::

For $z = x + yi$ in the right half-plane, $b = \sqrt{1 + u^2}$ with $u =
y/(\sqrt{x^2+y^2} + x)$, and in the left half-plane $u = (\sqrt{x^2+y^2}
- x)/y$; these are the two forms of $\tan(\arg(z)/2)$ that avoid
cancellation. The bound is computed by `acb_gamma_stirling_bound` in
[`acb`](acb), with upper and lower bounds for $|z|$ and $|\arg z|$ over
the ball. On the positive real line $b = 1$. The bound is infinite if
the ball contains points of $(-\infty, 0]$.

### Argument reduction {#reduction}

Since $|B_{2n}| \approx 2 (2n)!/(2\pi)^{2n}$, the terms of the series
decrease only while $2n \lesssim 2\pi |z|$, and the smallest term is
about $e^{-2\pi|z|}$. An error below $2^{-p}$ therefore requires $|z|
\ge \beta p$ with $\beta > \log 2/(2\pi) \approx 0.11$. If $z$ is
smaller, it is shifted with
$$
\Gamma(z) = \frac{\Gamma(z+r)}{(z)_r}, \qquad (z)_r = z(z+1)\cdots(z+r-1),
$$
until $\operatorname{Re}(z) + r \ge 1$ and $|z + r| \ge \max(1, \beta
p)$. A larger $\beta$ means fewer terms but a longer rising factorial;
FLINT uses $\beta = 0.17$ up to 32768 bits, $0.20$ up to 131072 bits
and $0.24$ above, which balances the two costs. The number of terms
$n$ is the smallest one for which the bound of [](#thm-stirling) with
$k = 0$, estimated in double precision with `bernoulli_bound_2exp_si`,
is below $2^{-p}$. For $\operatorname{Re}(z) < -5$ the reflection
formula
$$
\Gamma(z) = \frac{\pi}{\sin(\pi z)\, \Gamma(1-z)}
$$
is applied first. Arguments with a huge real or imaginary part need no
shift and few terms; this case is detected from the exponents without
converting to doubles.

The working precision is the target precision plus a few guard bits,
plus the number of bits in the exponent of $z$ for large $|z|$ (to
compensate for the cancellation in $(z - \tfrac12)\log z - z$), but at
most the relative accuracy of the input plus 20 bits.

The Bernoulli numbers $B_2, \ldots, B_{2n}$ come from the global cache
of the [`bernoulli`](bernoulli) module
([](algorithms/number-theory/bernoulli)). Since $n = O(p)$ and $B_{2n}$
has $O(n \log n)$ bits, filling the cache costs more than one
evaluation at the same precision, but it is done only once.

### Evaluating the sum {#stirling-sum}

The sum $\sum_{k=1}^{n-1} B_{2k}/(2k(2k-1)z^{2k-1})$ is computed by
`arb_hypgeom_gamma_stirling_sum_horner` (and its complex counterpart)
up to 128 bits, and for moderately many terms up to about 2048 bits:
Horner's scheme in $w = 1/z^2$, where the small terms are added at a
reduced precision.

At higher precision `arb_hypgeom_gamma_stirling_sum_improved` uses the
following reformulation [@Joh2023g]. With $B_{2k} = (-1)^{k+1} 2
(2k)!\, \zeta(2k)/(2\pi)^{2k}$ and $\zeta(2k) = \sum_{j \ge 1}
j^{-2k}$, each term of the sum splits into the contributions
$$
\frac{(-1)^{k+1} \, 2 (2k-2)!\, z}{(2\pi j z)^{2k}}, \qquad j = 1, 2,
\ldots
$$
The terms with $k < M$ are computed with the exact Bernoulli numbers,
by rectangular splitting in $w$ with coefficients combined over common
denominators in groups of up to four. For $k \ge M$, the contributions
of $j = 1, \ldots, K$ are summed separately: for fixed $j$ they form a
polynomial in $u_j = -1/(2\pi j z)^2$ whose consecutive coefficients
differ by the small integer factors $(2k)(2k-1)$, so that the sum needs
one full multiplication only every $\sqrt{n}$ terms. The contributions
of $j > K$ are bounded with the Hurwitz zeta function $\zeta(2M,
K+1)$, and the cutoffs $M_j \ge M$ for each $j$ are chosen so that all
omitted parts are below $2^{-p}$. The splitting parameter $K$ (at most
100) grows like $\sqrt{p}$.

## Rising factorials {#rising}

FLINT has three algorithms for the rising factorial $(x)_n$.

- **Forward recurrence.** $n-1$ multiplications at a precision of $p +
  \log_2 n$ bits, enough to absorb the rounding errors. The increments
  $x + k$ are formed in place, without a full addition, when the
  exponent permits.
- **Binary splitting.** The product $\prod_{k<n}(x+k)$ as a balanced
  tree ([](algorithms/integers/binary-splitting)). If $x$ has $b$
  bits, the subproducts of length $\ell$ have about $\ell b$ bits until
  they reach the precision $p$, so the method costs $O(M(p) \log n)$
  when $nb \le p$, much less than $n$ full multiplications.
- **Rectangular splitting** [@Smi2001; @Joh2014b]. Split the product
  into $n/m$ blocks of $m$ factors and expand each block as a
  polynomial in $x$:
  $$
  (x + k)_m = x^m + \sum_{i=0}^{m-1} c_{k,i}\, x^i,
  $$
  with integer coefficients $c_{k,i}$ computed by the recurrence of
  `_arb_hypgeom_rising_coeffs_1` (in one word if they fit, in two
  words, or as `fmpz`). With the powers $x^2, \ldots, x^m$ computed
  once, each block is one dot product (`arb_dot_ui`, `arb_dot_uiui` or
  `arb_dot_fmpz`) of $m$ balls with small integers, and the blocks are
  multiplied together. The cost is about $m + n/m$ multiplications of
  full-precision numbers and $n$ multiplications of balls by small
  integers (plus $O(nm)$ word operations for the coefficients), instead
  of $n$ full multiplications; with $m \approx \sqrt n$ this is
  $O(\sqrt n)$ full multiplications. In practice $m$ is limited to 60
  so that the coefficients stay small.

`arb_hypgeom_rising_ui_rec` uses the forward recurrence for small $n$,
binary splitting when $x$ has at most $p/8$ bits, and rectangular
splitting otherwise; the complex version is analogous. For $n \ge
\max(p, 100)$, `arb_hypgeom_rising_ui` computes $\Gamma(x+n)/\Gamma(x)$
instead.

The jets $(x+t)_n \in \mathbb{R}[t]/(t^\ell)$ (`arb_hypgeom_rising_ui_jet`)
use the same methods with truncated polynomial arithmetic, plus a
direct method for small $n$: the coefficients of $(x+t)_n$ are integer
combinations of the powers of $x$ (with Stirling numbers of the first
kind as coefficients). They are used for the power series gamma
function and for the log-rising factorial.

## Taylor series of the reciprocal gamma function {#taylor}

The reciprocal gamma function is entire, and
$$
\frac{1}{\Gamma(1+x)} = \sum_{k=0}^{\infty} c_k x^k, \qquad c_0 = 1,
\; c_1 = \gamma, \; c_2 = \frac{\gamma^2}{2} - \frac{\pi^2}{12}, \ldots
$$
with coefficients that decrease faster than geometrically. FLINT
stores the first 536 coefficients to 3456 bits as static data
(`arb_hypgeom_gamma_coeffs`) and evaluates $1/\Gamma(1+u)$ for $|u|
\le 1/2$ with Horner's scheme (`arb_hypgeom_gamma_taylor`) or with
rectangular splitting and dot products (`acb_hypgeom_gamma_taylor`),
at a precision that decreases with the size of the terms; when the
target precision is low, the first terms are summed in double
precision. The number of terms $N$ is the first one for which the
tabulated exponent of $c_N$ and $|u|^N$ show that the tail is below
$2^{-p}$. Since the coefficients decay super-exponentially over the
range of the table, the first omitted term bounds the tail up to a
small factor, which is included; the complex version uses an explicit
bound for the tail that is valid for $|u| \le 16$ and the truncation
orders in the table.

For $x = u + r$ with an integer $r \ge 1$,
$$
\Gamma(x) = \frac{(1+u)_{r-1}}{1/\Gamma(1+u)},
$$
and for negative $r$ the formula is combined with $\Gamma(x) =
\Gamma(x + 1 - r)/(x)_{1-r}$. The method is fast because it needs no
logarithm, exponential or Bernoulli numbers. It is used for $|x| <
160$, precisions up to 3456 bits, and shifts that are not long
compared to the precision ($-(40 + (p - 40)/4) \le r \le 70 + (p -
40)/8$). For complex $z$ the imaginary part must also be small (at most
4 to 10, depending on the precision).

For inexact real input, the result is computed at the midpoint and the
propagated error is bounded with $|\Gamma'(x)| = |\Gamma(x)\psi(x)|$,
using a piecewise polynomial upper bound for $|\psi|$. Wide real balls
(radius above $2^{-16}$ or precision below 40 bits) are enclosed
directly from double-precision upper and lower bounds for $\Gamma$ at
the endpoints, using that $\Gamma$ decreases on $(0, x_0)$ and
increases on $(x_0, \infty)$ with $x_0 \approx 1.4616$.

## Rational arguments {#rational}

For $x = p/q$ with small $q$, faster methods exist.

**Integers.** $\Gamma(n) = (n-1)!$ is computed exactly with
`fmpz_fac_ui` when $n < 1200$ or $n\log n$ is small compared to the
precision.

**Denominators 2, 3, 4 and 6.** With $x = a + m$, $0 < a \le 1$, the
rising factorial $(a)_m$ (computed with `arb_rising_fmpq_ui`) reduces
the problem to $\Gamma(a)$, and
$$
\Gamma(\tfrac12) = \sqrt\pi, \qquad
\Gamma(\tfrac23) = \frac{2 \pi}{\sqrt 3\, \Gamma(1/3)}, \qquad
\Gamma(\tfrac34) = \frac{\sqrt 2\, \pi}{\Gamma(1/4)},
$$
$$
\Gamma(\tfrac16) = \frac{\Gamma(1/3)^2}{(\pi/3)^{1/2}\, 2^{1/3}},
\qquad
\Gamma(\tfrac56) = \frac{2 \pi (\pi/3)^{1/2}\, 2^{1/3}}{\Gamma(1/3)^2}
$$
reduce it to the two constants $\Gamma(1/3)$ and $\Gamma(1/4)$. These
are computed by the [`mp_real`](mp_real) module and cached.
$\Gamma(1/3) = (810^{1/4} \pi / X)^{1/3}$, where $X$ is a rapidly
convergent hypergeometric series found by Guillera (2023), summed by
binary splitting in parallel with $\pi$. $\Gamma(1/4)$ comes from a
series of Ebisu for the lemniscate constant and, from about four
million bits, from
$$
\Gamma(\tfrac14) = \sqrt{\frac{(2\pi)^{3/2}}{\operatorname{agm}(1, \sqrt 2)}},
$$
whose arithmetic-geometric mean costs $O(M(p) \log p)$. Formulas for
other small fractions through elliptic integrals are given by Borwein
and Zucker [@BZ1992].

**General rational numbers at high precision.** For other denominators,
above about $7000 + 300\log_2 q$ bits and for $|x|$ not too large,
`arb_hypgeom_gamma_fmpq` writes $x = a + m$ with $0 < a \le 1$,
computes $\Gamma(1 + a)$ and splits it at $N \approx \alpha p\log 2$
(with $\alpha = 0.52$, and $N$ rounded to a number with trailing zero
bits) as $\Gamma(s) = \gamma(s, N) + \Gamma(s, N)$, $s = 1 + a$:
$$
\Gamma(s) = N^s e^{-N} \left(\frac{1}{s}\,{}_1F_1(1; 1+s; N) +
\frac{1}{N}\,{}_2F_0(1, 1-s; ; -1/N)\right).
$$
The convergent series is truncated after $n$ terms and bounded by its
first omitted term times a geometric series; the asymptotic series is
truncated after $n_2$ terms, with error at most $n_2!/N^{n_2}$, and is
computed at a lower precision since its contribution is multiplied by
$1/N$. The numbers of terms are determined with the Lambert $W$
function. Both series have rational term ratios in $N$ and the numerator
and denominator of $s$, so they are summed by binary splitting in
$O(M(p) \log^2 p)$ [@Bor1987; @Bre1976]. Below this precision, or for
large denominators, the Taylor or the Stirling series is applied to
the rational number rounded to a ball, with the working precision
raised by the size of $x$.

## The logarithm of the gamma function {#lgamma}

$\log\Gamma(z)$ denotes the analytic continuation of the real function
from $(0, \infty)$ to $\mathbb{C} \setminus (-\infty, 0]$, which differs from
$\log(\Gamma(z))$ by multiples of $2\pi i$. The Stirling series gives
this branch directly. The difficulty lies in the argument reduction,
which needs the log-rising factorial
$$
\log\,(z)_r = \sum_{k=0}^{r-1} \log(z + k)
$$
with the sum of the principal logarithms, not the principal logarithm
of the product. `acb_hypgeom_log_rising_ui` computes $(z)_r$, takes
the principal logarithm and adds $2\pi i m$, where $m$ is found with
Hare's algorithm [@Hare1997]: the partial products $(z)_k$ are formed
in double precision, rescaled to avoid overflow, and $m$ counts the
crossings of the negative real axis. This is valid when the input is
accurate enough and not too close to the real axis, so that each
partial product is determined to within a quadrant; otherwise the
imaginary part is computed as $\sum_k \arg(z + k)$, first at low
precision to identify $m$. The reflection formula becomes
$$
\log\Gamma(z) = \log \pi - \log\sin(\pi z) - \log\Gamma(1-z),
$$
where `acb_log_sin_pi` computes the branch of $\log \sin(\pi z)$ that
is analytic in the upper and lower half-planes.

At low and medium precision near the real axis, `acb_hypgeom_lgamma`
uses the Taylor series instead: it reduces the real part to about
$[1/2, 3/2]$, computes $-\log(1/\Gamma(u))$, adds $k\pi i$ according to
a table of linear approximations of the curves where
$\operatorname{Im}\log\Gamma(u) = (k + \tfrac12)\pi$, and undoes the
shift with the log-rising factorial (or the reflection formula for
real part below $-3$).

For real $x > 0$, `arb_hypgeom_lgamma` applies the Stirling series to
$\log\Gamma$ directly; for exact small fractions and for arguments in
the range of the Taylor series it takes the logarithm of $\Gamma(x)$.

## Choosing a method {#choosing}

| Argument | Method |
|---|---|
| exact integer $n$ | exact factorial for moderate $n$; otherwise as a real number |
| exact $p/q$, $q \in \{2,3,4,6\}$, moderate $p$ | rising factorial and the constants $\sqrt\pi$, $\Gamma(1/3)$, $\Gamma(1/4)$ |
| exact $p/q$, other $q$, very high precision | hypergeometric series by binary splitting |
| real or complex, $\lvert z\rvert < 160$, precision $\le 3456$ bits, near the real axis | Taylor series of $1/\Gamma$ |
| everything else | Stirling series with reflection and argument reduction |

`arb_hypgeom_gamma` tries the methods in the order of the table (exact
multiples of $1/4$ go to the rational code), and `acb_hypgeom_gamma`
uses the real code for real input. The reciprocal gamma function uses
the same methods with the reciprocal formulas; it is zero at the poles,
and for real balls containing 0 it is computed as $x \cdot
(1/\Gamma(x+1))$. For complex balls on which the Stirling series gives
no finite result, `acb_hypgeom_rgamma` returns the bound $|1/\Gamma(z)|
\le R^{R+1/2} e^{\pi R/2}$ for $|z| \le R$.

## Use in FLINT {#use}

- `arb_hypgeom_gamma`, `arb_hypgeom_rgamma`, `arb_hypgeom_lgamma`,
  `arb_hypgeom_gamma_fmpq`, `arb_hypgeom_gamma_fmpz` and the complex
  `acb_hypgeom_gamma`, `acb_hypgeom_rgamma`, `acb_hypgeom_lgamma` are
  the implementations; `arb_gamma`, `arb_rgamma`, `arb_lgamma`,
  `arb_gamma_fmpq`, `arb_gamma_fmpz`, `acb_gamma`, `acb_rgamma`,
  `acb_lgamma` and the generic `gr_gamma`, `gr_rgamma` and `gr_lgamma`
  for ball contexts call them.
- `arb_hypgeom_gamma_stirling`, `arb_hypgeom_gamma_taylor`,
  `arb_hypgeom_gamma_stirling_sum_horner`,
  `arb_hypgeom_gamma_stirling_sum_improved` and their complex
  counterparts expose the individual methods.
- The digamma functions `arb_digamma` and `acb_digamma` and the power
  series functions `arb_poly_gamma_series`, `arb_poly_lgamma_series`,
  `arb_poly_rgamma_series`, `arb_poly_digamma_series` and their complex
  versions use the same parameter choice and the bound of
  [](#thm-stirling) for higher derivatives (`acb_gamma_stirling_bound`).
- The rising factorials `arb_hypgeom_rising_ui` and
  `acb_hypgeom_rising_ui` and their jets implement `arb_rising_ui`,
  `acb_rising_ui` and the related functions of [`arb`](arb) and
  [`acb`](acb).
- `arb_hypgeom_central_bin_ui` uses the analogous asymptotic expansion
  of the central binomial coefficient with Brent's error bound
  [@Bre2016].
