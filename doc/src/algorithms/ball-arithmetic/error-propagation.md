---
title: Error propagation
summary: Bounds for the propagated error of arithmetic operations, roots and elementary functions on balls, with proofs, together with the tail bounds for series that are used throughout the special function code.
complexity: Each bound costs O(1) operations on mag_t numbers.
order: 20
---

A ball function evaluates its result at the midpoint of the input and
adds to the radius a bound for the *propagated error*: how much the exact
function can differ between the midpoint and any other point of the input
ball. For a function $f$ and a ball $[x \pm r]$, the task is to bound

$$
|f(x + a) - f(x)| \quad \text{for all } |a| \le r,
$$

in terms of $|x|$ and $r$. Most bounds below increase with $|a|$, so they
are stated in terms of $|a|$ and applied with $|a| = r$. The rounding
error of the midpoint computation is added separately (see
[](algorithms/ball-arithmetic/representation)). All bounds are evaluated
with `mag_t` arithmetic, rounding upward, from upper bounds for the
quantities that occur in numerators and lower bounds for those in
denominators.

## Arithmetic {#arithmetic}

Addition is immediate: $|(x + a) + (y + b) - (x + y)| \le |a| + |b|$, so
`arb_add` adds the radii.

::: theorem Product {#thm-mul}
For $x, y, a, b \in \mathbb{C}$,
$$
|(x + a)(y + b) - xy| \le |x||b| + |y||a| + |a||b|.
$$
:::

::: proof
Expand $(x + a)(y + b) - xy = xb + ya + ab$ and apply the triangle
inequality.
:::

The bound is attained for real numbers when all of $x, y, a, b$ have the
same sign, so for real balls it is optimal up to rounding. `arb_mul`
computes the radius $|m_x| r_y + |m_y| r_x + r_x r_y$ with
`mag_fast_mul` and `mag_fast_addmul`, using upper bounds of $|m_x|$ and
$|m_y|$ obtained from the top 30 bits of the midpoints.

::: theorem Quotient {#thm-div}
For $x, y, a, b \in \mathbb{C}$ with $|b| < |y|$,
$$
\left| \frac{x}{y} - \frac{x + a}{y + b} \right| \le
\frac{|x||b| + |y||a|}{|y| (|y| - |b|)}.
$$
:::

::: proof
We have $\dfrac{x}{y} - \dfrac{x + a}{y + b} = \dfrac{x(y + b) - y(x +
a)}{y(y + b)} = \dfrac{xb - ya}{y(y + b)}$. The triangle inequality
bounds the numerator by $|x||b| + |y||a|$, and the reverse triangle
inequality gives $|y + b| \ge |y| - |b| > 0$.
:::

When $r_y$ is much smaller than $|m_y|$, `arb_div` uses the simplified
form

$$
\frac{|x| r_y + |y| r_x}{|y|(|y| - r_y)} \le (1 + \varepsilon)
\frac{|x/y| r_y + r_x}{|y|}, \quad \frac{|y|}{|y| - r_y} \le 1 +
\varepsilon,
$$

with the computed quotient in place of $x/y$ and a fixed factor $1 +
2^{-16}$ for $1 + \varepsilon$, which is valid when the relative accuracy
of $y$ is more than 20 bits. For less accurate input it evaluates the
quotient at the endpoints instead, which gives a tighter enclosure.

## Roots {#roots}

::: theorem Square root {#thm-sqrt}
For real $x$, $a$ with $0 \le |a| \le x$,
$$
|\sqrt{x + a} - \sqrt{x}| \le \sqrt{x} \left(1 - \sqrt{1 - \frac{|a|}{x}}\right)
\le \frac{\sqrt{x}}{2} \left( \frac{|a|}{x} + \frac{|a|^2}{x^2} \right),
$$
where the first inequality is an equality for $a \le 0$ (and the limit
value is 0 when $x = a = 0$).
:::

::: proof
Since the square root is concave, a perturbation of size $|a|$ changes it
most when it is negative, so the error is at most $\sqrt{x} - \sqrt{x -
|a|} = \sqrt{x}(1 - \sqrt{1 - t})$ with $t = |a|/x \in [0, 1]$. For the
second inequality we need $\sqrt{1 - t} \ge 1 - (t + t^2)/2$. The right
side is nonnegative on $[0, 1]$, and squaring turns the inequality into
$1 - t \ge 1 - t - t^2 + (t + t^2)^2/4$, that is, $t^2 (t + 3)(t - 1) / 4
\le 0$, which holds for $0 \le t \le 1$.
:::

::: theorem Reciprocal square root {#thm-rsqrt}
For real $x, a$ with $0 \le |a| < x$,
$$
\left| \frac{1}{\sqrt{x + a}} - \frac{1}{\sqrt{x}} \right| \le
\frac{|a|}{2 (x - |a|)^{3/2}}.
$$
:::

::: proof
By the mean value theorem the difference is $|a| \cdot |f'(\xi)|$ for
some $\xi$ between $x$ and $x + a$, where $f(t) = t^{-1/2}$ and $|f'(t)| =
t^{-3/2}/2$ is decreasing, so $|f'(\xi)| \le (x - |a|)^{-3/2}/2$.
:::

::: theorem $k$-th root {#thm-root}
For $k > 1$ and real $x, a$ with $0 \le |a| \le x$,
$$
|(x + a)^{1/k} - x^{1/k}| \le x^{1/k} \min\left(1, \frac{1}{k}
\log\left(1 + \frac{|a|}{x - |a|}\right)\right).
$$
:::

::: proof
By concavity the error is largest for $a = -r$, $r = |a|$, and

$$
x^{1/k} - (x - r)^{1/k} = x^{1/k} \left[1 - \exp\left(\tfrac{1}{k} \log(1
- r/x)\right)\right] \le x^{1/k} \min\left(1, -\tfrac{1}{k} \log(1 -
r/x)\right),
$$

using $1 - e^{-u} \le \min(1, u)$ for $u \ge 0$. Finally $-\log(1 - r/x) =
\log(1 + r/(x - r))$.
:::

## Elementary functions {#elementary}

The following bound is the general principle; the theorems after it are
its sharp special cases.

::: theorem Generic first-order bound {#thm-generic}
Let $f$ be differentiable on a neighbourhood of the segment from $x$ to
$x + a$ (real or complex). Then
$$
|f(x + a) - f(x)| \le \min(2 C_0, C_1 |a|), \quad C_0 = \sup_{|t| \le
|a|} |f(x + t)|, \quad C_1 = \sup_{|t| \le |a|} |f'(x + t)|.
$$
:::

::: proof
The first bound is the triangle inequality, the second the mean value
inequality $|f(x + a) - f(x)| = |\int_0^1 f'(x + sa)\, a\, ds| \le C_1
|a|$.
:::

::: theorem Exponential {#thm-exp}
For real $x, a$: $|e^{x + a} - e^x| = e^x |e^a - 1| \le e^x (e^{|a|} -
1)$, with equality if $a \ge 0$.
:::

::: proof
For $a \ge 0$ this is an identity. For $a < 0$, $|e^a - 1| = 1 - e^{-|a|}
\le e^{|a|} - 1$, since $e^u + e^{-u} \ge 2$.
:::

`arb_exp` therefore adds $e^m (e^r - 1)$, computed with `mag_expm1` and
an upper bound for the computed $e^m$, to the radius.

::: theorem Logarithm {#thm-log}
For real $x, a$ with $0 \le |a| < x$,
$$
|\log(x + a) - \log(x)| \le \log\left(1 + \frac{|a|}{x - |a|}\right),
$$
with equality if $a \le 0$.
:::

::: proof
The logarithm is concave and increasing, so the largest change is for $a
= -|a|$, where $\log x - \log(x - |a|) = \log(x/(x - |a|)) = \log(1 +
|a|/(x - |a|))$.
:::

`arb_log` uses this bound with `mag_log1p` when the input has more than
about 20 bits of relative accuracy, and otherwise evaluates $\log$ at both
endpoints with the `mag_t` logarithm bounds.

::: theorem Sine and cosine {#thm-sin}
For real $x, a$: $|\sin(x + a) - \sin x| \le \min(2, |a|, |\cos x| |a| +
|a|^2/2)$, and the same holds for $\cos$ with $|\sin x|$ in place of
$|\cos x|$.
:::

::: proof
The bounds 2 and $|a|$ follow from [](#thm-generic). For the third,
Taylor's theorem with the Lagrange remainder gives $\sin(x + a) = \sin x
+ a \cos x - \tfrac{a^2}{2} \sin \xi$, and $|\sin \xi| \le 1$.
:::

`arb_sin_cos` uses the third form, with $|\cos x|$ bounded from an 8-bit
table entry, which is much tighter than $|a|$ near the extrema of the
sine. Wide balls (radius above about $2^{-24}$) are instead handled by
evaluating at the endpoints and locating the extrema in between
(`arb_sin_cos_wide`).

::: theorem Arctangent {#thm-atan}
For real $x, a$,
$$
|\operatorname{atan}(x + a) - \operatorname{atan}(x)| \le \min(\pi, C_1
|a|), \quad C_1 = \sup_{|t| \le |a|} \frac{1}{1 + (x + t)^2}.
$$
If $|a| < |x|$, one may take $C_1 = (1 + (|x| - |a|)^2)^{-1}$. Exactly,
if $|a| < |x|$ or $|x(x + a)| < 1$, the difference equals
$\operatorname{atan}(|a| / |1 + x(x + a)|)$ up to the choice of sign.
:::

::: proof
The first statement is [](#thm-generic) with $f' = 1/(1 + t^2)$ and the
range of $\operatorname{atan}$ of length $\pi$. If $|a| < |x|$, every
point between $x$ and $x + a$ has absolute value at least $|x| - |a|$,
which gives the monotone choice of $C_1$. The exact form is the
subtraction formula $\operatorname{atan} u - \operatorname{atan} v =
\operatorname{atan}((u - v)/(1 + uv))$, valid when $uv > -1$.
:::

`arb_atan` uses $|a| / (1 + (|x| - |a|)^2)$, capped by $\pi$, for accurate
input and the endpoints otherwise.

## Sums and series {#series}

Many special functions are evaluated by truncating a series. The
truncated sum is computed in ball arithmetic, and a bound for the tail is
added to the radius. The following bounds are the most common.

::: theorem Geometric bound {#thm-geometric}
If $|c_k| \le C$ for all $k \ge N$ and $|z| \le D < 1$, then
$$
\left| \sum_{k=N}^{\infty} c_k z^k \right| \le \frac{C D^N}{1 - D}.
$$
:::

::: proof
Compare termwise with the geometric series $\sum_{k \ge N} C D^k$.
:::

::: theorem Integral bound {#thm-integral}
If $f$ is nonnegative and decreasing on $[N, \infty)$, then
$$
\int_N^{\infty} f(t)\, dt \le \sum_{k=N}^{\infty} f(k) \le f(N) +
\int_N^{\infty} f(t)\, dt.
$$
:::

::: proof
For $k \ge N$, $f(k + 1) \le \int_k^{k+1} f \le f(k)$; sum over $k$.
:::

For example, `mag_hurwitz_zeta_uiui` bounds $\zeta(s, a) \le a^{-s} +
a^{1-s}/(s - 1)$ in this way, and `mag_exp_tail` uses the geometric bound
for the exponential series: if $N \ge 2|x|$, consecutive terms of $\sum_{k
\ge N} x^k/k!$ decrease by a factor at most $1/2$, so the tail is at most
twice its first term.

For analytic functions, the coefficients can be bounded by Cauchy's
integral formula. If $f(z) = \sum_{k=0}^{\infty} c_k z^k$ is analytic on
a neighbourhood of the disk $|z| \le R$, then

$$
c_k = \frac{1}{2\pi i} \int_{|z| = R} \frac{f(z)}{z^{k+1}}\, dz,
$$

which gives the following two bounds.

::: corollary Cauchy bounds {#cor-cauchy}
With $C = \max_{|z| = R} |f(z)|$, the coefficients satisfy $|c_k| \le C /
R^k$, and for $0 \le r < R$ and $|z| \le r$,
$$
\left| \sum_{k=N}^{\infty} c_k z^k \right| \le \frac{C D^N}{1 - D},
\quad D = \frac{r}{R}.
$$
:::

::: proof
Bound the integrand by $C / R^{k+1}$ on a contour of length $2\pi R$; the
tail bound is then [](#thm-geometric) applied to $c_k R^k \cdot (z/R)^k$.
:::

The Euler-Maclaurin formula and its remainder bounds, used for the zeta
function, are discussed in [](algorithms/special-functions/zeta).

## Use in FLINT {#use}

- `arb_add`, `arb_sub`, `arb_mul`, `arb_addmul` and their variants use
  the sum and [](#thm-mul); `acb_mul` applies [](#thm-mul) to complex
  balls in rectangular form componentwise.
- `arb_div` and `arb_inv` use [](#thm-div) or the endpoints.
- `arb_sqrt`, `arb_rsqrt` and `arb_root_ui` use [](#thm-sqrt),
  [](#thm-rsqrt) and [](#thm-root).
- `arb_exp`, `arb_log`, `arb_sin_cos` and `arb_atan` use the
  corresponding theorems, with endpoint evaluation for inaccurate input.
- The tail bounds are implemented by `mag_geom_series`, `mag_exp_tail`,
  `mag_hurwitz_zeta_uiui` and `mag_polylog_tail`, and are used throughout
  [arb_hypgeom](arb_hypgeom), [acb_hypgeom](acb_hypgeom) and the zeta
  function code.

::: functions arb_mul arb_div arb_sqrt arb_exp arb_log arb_sin_cos arb_atan :::
