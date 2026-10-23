---
title: acb_hypgeom.h
header: acb_hypgeom.h
description: Rigorous evaluation of generalized hypergeometric functions and their special cases (gamma, error, exponential-integral, Bessel, Airy, Coulomb, Legendre, Gauss hypergeometric, orthogonal polynomials, dilogarithm) for complex ball and power series arguments.
category: Special functions
order: 2
---

This module evaluates hypergeometric functions of complex arguments
given as balls (`acb_t`) or as complex power series (`acb_poly_t`),
with rigorous error bounds. The generalized hypergeometric function is
defined by the series
$$
{}_pF_q(a_1,\ldots,a_p;b_1,\ldots,b_q;z) =
\sum_{k=0}^\infty \frac{(a_1)_k\cdots(a_p)_k}{(b_1)_k\cdots(b_q)_k} \frac{z^k}{k!},
$$
interpreted by analytic continuation (or as a formal series together
with an asymptotic expansion) where it does not converge. In a wider
sense, "hypergeometric functions" are linear combinations of such
functions with prefactors that are products of exponentials, powers and
gamma functions. This covers most classical special functions, and
this module implements the gamma function and rising factorials, the
confluent functions ${}_0F_1$, ${}_1F_1$ and $U$, error functions and
Fresnel integrals, incomplete gamma and beta functions, exponential,
sine and cosine integrals, Bessel, Airy and Coulomb wave functions, the
Gauss function ${}_2F_1$, Legendre functions, orthogonal polynomials,
spherical harmonics and the dilogarithm. The real versions are in
[`arb_hypgeom`](arb_hypgeom), which calls this module for most
functions, and the generic interface [`gr_special`](gr_special) calls
both for its real and complex ball contexts.

Each function chooses between convergent series, asymptotic
expansions, transformation formulas, analytic continuation and, for
real arguments, numerical integration, and sets its working precision
from an estimate of the cancellation; the algorithms and error bounds
are described in [](algorithms/special-functions/hypergeometric), and
[](algorithms/special-functions/overview) maps each function
to its methods. The guide [](acb_hypgeom/using) explains the
conventions with examples.

::: summary algorithms/special-functions/hypergeometric :::

## Conventions {#conventions}

- **Enclosures.** Every function returns a ball containing the value of
  the function at every point of the input balls. If the function is
  undefined somewhere in the input (a pole, or an input that straddles
  a branch cut in a way the function cannot handle), or if no finite
  enclosure could be computed, the output is indeterminate (see
  `acb_indeterminate`). A wide or indeterminate output for a valid
  input usually means that the precision is too low for the size of
  the parameters; increasing `prec` helps.
- **Precision.** `prec` is the target precision. The functions estimate
  the cancellation in the series they use and increase the working
  precision internally, and reduce it to what the accuracy of the input
  justifies. For inexact input they usually evaluate at the midpoint
  and bound the propagated error with derivative bounds.
- **Branch cuts.** Multivalued functions use the principal branch, with
  the branch cut on $(-\infty, 0]$ for functions built from $\log z$
  and $z^a$ ($U$, $\operatorname{Ei}$, $\operatorname{Ci}$, Coulomb
  functions, the log-gamma function) and on $[1, +\infty)$ for
  functions built from $\log(1-z)$ (${}_2F_1$, the dilogarithm). On a
  cut the value follows `acb_log`, which is continuous from the upper
  half-plane on $(-\infty, 0)$; for a cut on $(1, +\infty)$ this means
  the limit from the lower half-plane. The conventions for the
  exponential and trigonometric integrals match *Mathematica*.
- **Regularization.** An `int` argument `regularized` (or `flags` for
  ${}_2F_1$) selects the regularized function, divided by the gamma
  functions of the lower parameters so that it is entire in them; the
  regularized function is defined at the nonpositive integers where the
  ordinary one has poles.
- **Aliasing and optional outputs.** Inputs and outputs may be the same
  object unless stated otherwise. Functions with several outputs accept
  `NULL` for unwanted outputs. The underscore power series functions
  require positive lengths, write exactly `len` coefficients and do not
  allow the output to overlap the input.
- **Wrappers.** For real input, several functions call the faster real
  implementation in [`arb_hypgeom`](arb_hypgeom) (gamma, rising
  factorials, error function, Legendre polynomials).

## Rising factorials {#rising}

::: function
void acb_hypgeom_rising_ui_forward(acb_t res, const acb_t x, ulong n, slong prec)
void acb_hypgeom_rising_ui_bs(acb_t res, const acb_t x, ulong n, slong prec)
void acb_hypgeom_rising_ui_rs(acb_t res, const acb_t x, ulong n, ulong m, slong prec)
void acb_hypgeom_rising_ui_rec(acb_t res, const acb_t x, ulong n, slong prec)
void acb_hypgeom_rising_ui(acb_t res, const acb_t x, ulong n, slong prec)
void acb_hypgeom_rising(acb_t res, const acb_t x, const acb_t n, slong prec)

Sets `res` to the rising factorial $(x)_n = x(x+1)\cdots(x+n-1)$;
`acb_hypgeom_rising` allows any complex $n$ and computes
$\Gamma(x+n)/\Gamma(x)$.

::: details
The `forward`, `bs` and `rs` versions use the forward recurrence,
binary splitting and rectangular splitting (with block size $m$, or 0
to choose it), as for the real functions `arb_hypgeom_rising_ui_forward`
and its relatives. `acb_hypgeom_rising_ui_rec` calls the real version
for real $x$, and otherwise uses the forward recurrence for $n \le 5$
up to 512 bits, binary splitting for $n \ge 20$ if $x$ has fewer than
$\mathit{prec}/8$ bits, and rectangular splitting otherwise; it never
calls the gamma function. `acb_hypgeom_rising_ui` and
`acb_hypgeom_rising` use it for integers $0 \le n < \max(\mathit{prec},
100)$ and the gamma function otherwise. These are the functions behind
`acb_rising_ui` and `acb_rising`. See
[](algorithms/special-functions/gamma#rising).
:::
:::

::: function
void acb_hypgeom_rising_ui_jet_powsum(acb_ptr res, const acb_t x, ulong n, slong len, slong prec)
void acb_hypgeom_rising_ui_jet_bs(acb_ptr res, const acb_t x, ulong n, slong len, slong prec)
void acb_hypgeom_rising_ui_jet_rs(acb_ptr res, const acb_t x, ulong n, ulong m, slong len, slong prec)
void acb_hypgeom_rising_ui_jet(acb_ptr res, const acb_t x, ulong n, slong len, slong prec)

Sets the vector `res` of length `len` to the coefficients of the
polynomial $(x + t)_n \in \mathbb{C}[t]$, truncated if $\mathit{len} < n+1$ and
padded with zeros if $\mathit{len} > n+1$.

::: details
The algorithms and the automatic choice are those of the real versions
`arb_hypgeom_rising_ui_jet_powsum`, `arb_hypgeom_rising_ui_jet_bs`,
`arb_hypgeom_rising_ui_jet_rs` and `arb_hypgeom_rising_ui_jet`. The
output must not overlap $x$.
:::
:::

::: function
void acb_hypgeom_log_rising_ui(acb_ptr res, const acb_t x, ulong n, slong prec)

Sets `res` to the log-rising factorial $\log\,(x)_n = \sum_{k=0}^{n-1}
\log(x+k)$, the sum of principal logarithms, which differs from the
principal logarithm of $(x)_n$ by a multiple of $2\pi i$.

::: details
The function computes $(x)_n$, takes its principal logarithm and
determines the branch correction $2\pi i m$. The correction is found
with Hare's algorithm [@Hare1997] in double-precision arithmetic, which
counts how often the partial products cross the negative real axis,
when this is safe ($n$ and $|x|$ at most $10^6$, $|\operatorname{Im}
x| > 10^{-6}$ and an accurate input); otherwise the imaginary part is
computed directly as $\sum_{k} \arg(x+k)$, first at low precision to
determine $m$ and if necessary at full precision. Real $x$ needs no
correction for $x > 0$ and gets the imaginary part $\pi \min(n,
\lceil -x \rceil)$ for negative non-integer $x$; a real $x$ containing
a nonpositive integer gives an indeterminate result.
:::
:::

::: function
void acb_hypgeom_log_rising_ui_jet(acb_ptr res, const acb_t x, ulong n, slong len, slong prec)

Sets the vector `res` to the coefficients of $\log\,(x+t)_n$ as a power
series in $t$, truncated to length `len`, with the constant term as in
`acb_hypgeom_log_rising_ui`.
:::

## Gamma function {#gamma}

The algorithms for the gamma function are described in
[](algorithms/special-functions/gamma).

::: function
void acb_hypgeom_gamma(acb_t res, const acb_t x, slong prec)

Sets `res` to the gamma function $\Gamma(x)$. The result is
indeterminate if $x$ contains a pole.

::: details
Real input is passed to `arb_hypgeom_gamma`. Otherwise the function
uses the tabulated Taylor series of `acb_hypgeom_gamma_taylor` near the
real axis at low and moderate precision and the Stirling series of
`acb_hypgeom_gamma_stirling` elsewhere. See [](acb_hypgeom/gamma).
`acb_gamma` calls this function.
:::
:::

::: function
void acb_hypgeom_rgamma(acb_t res, const acb_t x, slong prec)

Sets `res` to the reciprocal gamma function $1/\Gamma(x)$, which is
entire.

::: details
Uses the same methods as `acb_hypgeom_gamma`. If the Stirling series
gives no finite enclosure (a wide ball), it returns the bound
$|1/\Gamma(z)| \le R^{R + 1/2} e^{\pi R/2}$ for $|z| \le R$ instead.
`acb_rgamma` calls this function.
:::
:::

::: function
void acb_hypgeom_lgamma(acb_t res, const acb_t x, slong prec)

Sets `res` to the principal branch of the log-gamma function $\log
\Gamma(x)$, which is analytic except on the branch cut $(-\infty, 0]$
and differs from $\log(\Gamma(x))$ by a multiple of $2\pi i$.

::: details
Positive real input is passed to `arb_hypgeom_lgamma`. Near the real
axis at moderate precision, the function computes $\log(1/\Gamma(u))$
with the Taylor series for $u$ with real part in $[1/2, 3/2]$ (where
the branch is determined by tabulated bounds for the imaginary part),
and corrects the shift and the reflection with
`acb_hypgeom_log_rising_ui` and `acb_log_sin_pi`. Elsewhere it uses
the Stirling series, which gives $\log\Gamma$ directly, with argument
reduction by the log-rising factorial and, for negative real parts,
the reflection formula $\log\Gamma(z) = \log\pi - \log\sin(\pi z) -
\log\Gamma(1-z)$ with the branch of `acb_log_sin_pi`. See
[](algorithms/special-functions/gamma#lgamma). `acb_lgamma` calls this
function.
:::
:::

::: function
void acb_hypgeom_gamma_stirling_sum_horner(acb_t s, const acb_t z, slong N, slong prec)
void acb_hypgeom_gamma_stirling_sum_improved(acb_t s, const acb_t z, slong N, slong K, slong prec)

Sets `s` to $\sum_{n=1}^{N-1} B_{2n} / (2n(2n-1) z^{2n-1})$, the sum
in the Stirling series truncated before the term of index $N$, using
Horner's scheme or the improved algorithm with splitting parameter $K$
(0 to choose it), as `arb_hypgeom_gamma_stirling_sum_horner` and
`arb_hypgeom_gamma_stirling_sum_improved` do for real $z$. No bound for
the truncation error is added.
:::

::: function
void acb_hypgeom_gamma_stirling(acb_t res, const acb_t x, int reciprocal, slong prec)

Sets `res` to $\Gamma(x)$, or $1/\Gamma(x)$ if `reciprocal` is
nonzero, computed with the Stirling series after argument reduction by
a rising factorial and, for real part below $-5$, the reflection
formula. See `arb_hypgeom_gamma_stirling` and
[](algorithms/special-functions/gamma#stirling).
:::

::: function
int acb_hypgeom_gamma_taylor(acb_t res, const acb_t x, int reciprocal, slong prec)

Attempts to compute $\Gamma(x)$, or $1/\Gamma(x)$ if `reciprocal` is
nonzero, with the tabulated Taylor series of $1/\Gamma(1+x)$ and
argument reduction. Returns 1 on success and 0, without changing
`res`, if $x$ or `prec` is outside the supported range.

::: details
The method applies when $|\operatorname{Im} x|$ is small (at most 4 to
10, depending on the precision), $|\operatorname{Re} x|$ is not too
large compared to the precision, and the working precision is at most
`ARB_HYPGEOM_GAMMA_TAB_PREC` bits. The Taylor polynomial is evaluated
by rectangular splitting with a rigorous bound for the truncation
error. See [](algorithms/special-functions/gamma#taylor).
:::
:::

## Generalized hypergeometric function {#pfq}

::: function
void acb_hypgeom_pfq(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, int regularized, slong prec)

Sets `res` to the generalized hypergeometric function ${}_pF_q(a_1,
\ldots, a_p; b_1, \ldots, b_q; z)$, or the regularized function
${}_pF_q / (\Gamma(b_1) \cdots \Gamma(b_q))$ if `regularized` is
nonzero.

::: details
The orders $(p, q) = (0,0), (1,0), (0,1), (1,1), (2,1)$ are passed to
`acb_exp`, the power $(1-z)^{-a}$, `acb_hypgeom_0f1`, `acb_hypgeom_m`
and `acb_hypgeom_2f1`, which implement analytic continuation and
asymptotic expansions. All other orders use direct summation of the
series with `acb_hypgeom_pfq_direct` (with a power series variant for
regularized functions whose lower parameters may be nonpositive
integers), which converges for $p \le q$ and for $p = q+1$, $|z| < 1$,
and gives a finite result for $p > q+1$ only when the series
terminates. Parameters that occur in both $a$ and $b$ are not removed
(except a single upper parameter equal to 1, which cancels the $k!$);
removing them first saves time. See [](acb_hypgeom/pfq).
:::
:::

## Confluent hypergeometric functions {#confluent}

In this section $U(a,b,z)$ is the confluent hypergeometric function of
the second kind with the principal branch cut, and $U^{*}(a,b,z) = z^a
U(a,b,z)$, which tends to 1 as $z \to \infty$ in the sector
$|\arg z| < 3\pi/2$.

::: function
void acb_hypgeom_m_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
void acb_hypgeom_m_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
void acb_hypgeom_m(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)
void acb_hypgeom_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)

Sets `res` to Kummer's confluent hypergeometric function $M(a,b,z) =
{}_1F_1(a,b,z)$, or $\mathbf{M}(a,b,z) = {}_1F_1(a,b,z)/\Gamma(b)$ if
`regularized` is nonzero. `acb_hypgeom_1f1` is the same function as
`acb_hypgeom_m`.

::: details
`acb_hypgeom_m_asymp` uses the asymptotic expansions of $U^{*}$ through
the connection formula
$$
\mathbf{M}(a,b,z) = \frac{(-z)^{-a}}{\Gamma(b-a)} U^{*}(a,b,z) +
\frac{z^{a-b} e^z}{\Gamma(a)} U^{*}(b-a,b,-z).
$$
`acb_hypgeom_m_1f1` sums the convergent series, after Kummer's
transformation $M(a,b,z) = e^z M(b-a,b,-z)$ if $\operatorname{Re} z <
0$ (unless the series terminates), and handles $\mathbf{M}$ at
nonpositive integer $b$ by the limit formula. `acb_hypgeom_m` chooses:
terminating series are summed directly; otherwise it compares the
accuracy of the asymptotic expansion (about $|z|/\log 2$ bits) with
that of the convergent series after cancellation, uses Kummer's
transformation for negative real part, and sets the working precision
accordingly. If the result has less than half the requested accuracy
and all arguments are real and moderate, it also tries
`arb_hypgeom_1f1_integration` and keeps the better enclosure. See
[](algorithms/special-functions/hypergeometric#confluent).
:::
:::

::: function
void acb_hypgeom_u_asymp(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong n, slong prec)

Sets `res` to $U^{*}(a,b,z)$ computed with $n$ terms of the asymptotic
series $U^{*}(a,b,z) \sim {}_2F_0(a, a-b+1; ; -1/z)$, with a rigorous
bound for the truncation error. Requires $n \ge 0$; $n < 0$ chooses the
number of terms automatically.

::: details
The error bound is that of Olver [@NIST2012, §13.7], with the case
distinction described in
[](algorithms/special-functions/hypergeometric#asymptotic). If $a$ or
$a - b + 1$ is a nonpositive integer the series terminates and the
result is exact up to rounding for every $z \ne 0$.
:::
:::

::: function
int acb_hypgeom_u_use_asymp(const acb_t z, slong prec)

Returns nonzero if the asymptotic series is expected to give
$U(a,b,z)$ to `prec` bits, namely if $|z| > \mathit{prec} \log 2$ or
$|z| > 2^{64}$ (assuming that $a$ and $b$ are small).
:::

::: function
void acb_hypgeom_u_1f1_series(acb_poly_t res, const acb_poly_t a, const acb_poly_t b, const acb_poly_t z, slong len, slong prec)

Sets `res` to $U(a,b,z)$ as a power series truncated to length `len`,
for power series $a, b, z \in \mathbb{C}[[x]]$, using the connection formula
with two ${}_1F_1$ series. If the constant term of $b$ is an integer,
one extra derivative is computed and the removable singularity is
divided out; it is then assumed that the coefficient of $x$ in $b$ is
nonzero. The output is indeterminate if the constant term of $b$ is
inexact and contains an integer.
:::

::: function
void acb_hypgeom_u_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)

Sets `res` to $U(a,b,z)$ computed from the convergent series by the
connection formula
$$
U(a,b,z) = \frac{\Gamma(1-b)}{\Gamma(a-b+1)} M(a,b,z) +
\frac{\Gamma(b-1)}{\Gamma(a)} z^{1-b} M(a-b+1,2-b,z).
$$
If $b$ is an integer the limit is computed with
`acb_hypgeom_u_1f1_series`; the output is indeterminate if $b$ is
inexact and contains an integer.
:::

::: function
void acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t z, slong prec)

Sets `res` to $U(a,b,z)$.

::: details
For small nonpositive integer $a$, $U$ is a polynomial in $z$ computed
by its three-term recurrence (with binary splitting for $-a \ge 20$).
If $a$ or $a - b + 1$ is a nonpositive integer, or `acb_hypgeom_u_use_asymp`
returns true, the (terminating or convergent enough) asymptotic series
is used. Otherwise the function compares the expected accuracy of the
asymptotic series with that of `acb_hypgeom_u_1f1`, whose cancellation
it estimates as $|z|/\log 2$ bits, and chooses the working precision.
As with `acb_hypgeom_m`, real arguments with an inaccurate result are
retried with `arb_hypgeom_u_integration`.
:::
:::

::: function
void acb_hypgeom_0f1_asymp(acb_t res, const acb_t a, const acb_t z, int regularized, slong prec)
void acb_hypgeom_0f1_direct(acb_t res, const acb_t a, const acb_t z, int regularized, slong prec)
void acb_hypgeom_0f1(acb_t res, const acb_t a, const acb_t z, int regularized, slong prec)

Sets `res` to the confluent hypergeometric limit function ${}_0F_1(a,
z)$, or ${}_0F_1(a,z)/\Gamma(a)$ if `regularized` is nonzero, using
asymptotic expansions, direct summation, or an automatic choice.

::: details
The asymptotic version uses the connection formulas with Bessel
functions
$$
\frac{{}_0F_1(a,z)}{\Gamma(a)} = (-z)^{(1-a)/2} J_{a-1}(2 \sqrt{-z}) =
z^{(1-a)/2} I_{a-1}(2 \sqrt{z}),
$$
using the $J$ form in the left half-plane and the $I$ form in the right
half-plane, so that the square root is never evaluated on its branch
cut, and the asymptotic expansions of the Bessel functions (through
`acb_hypgeom_u_asymp`). The automatic version uses the asymptotic
expansion when $|z|$ is large enough for the target precision.
:::
:::

## Error functions and Fresnel integrals {#erf}

The representations and the choice of method for this section and the
next two are described in [](algorithms/special-functions/error-functions).

::: function
void acb_hypgeom_erf_1f1a(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_erf_1f1b(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_erf_asymp(acb_t res, const acb_t z, int complementary, slong prec, slong prec2)

Sets `res` to $\operatorname{erf}(z)$ computed respectively from
$$
\operatorname{erf}(z) = \frac{2z}{\sqrt{\pi}}
{}_1F_1(\tfrac{1}{2}, \tfrac{3}{2}, -z^2), \qquad
\operatorname{erf}(z) = \frac{2z e^{-z^2}}{\sqrt{\pi}}
{}_1F_1(1, \tfrac{3}{2}, z^2),
$$
$$
\operatorname{erf}(z) = \frac{z}{\sqrt{z^2}}
\left(1 - \frac{e^{-z^2}}{\sqrt{\pi}} U(\tfrac{1}{2}, \tfrac{1}{2},
z^2)\right) = \frac{z}{\sqrt{z^2}} - \frac{e^{-z^2}}{z \sqrt{\pi}}
U^{*}(\tfrac{1}{2}, \tfrac{1}{2}, z^2).
$$

::: details
The first series has terms of alternating sign when $z$ is near the
real axis, the second has positive terms there and is the one to use
for real $z$; near the imaginary axis the roles are exchanged. The
asymptotic version computes the $U$ term at the second precision
`prec2`, which can be much lower than `prec` when the term is small,
and computes $\operatorname{erfc}(z)$ instead if `complementary` is
nonzero.
:::
:::

::: function
void acb_hypgeom_erf(acb_t res, const acb_t z, slong prec)

Sets `res` to the error function $\operatorname{erf}(z) =
\frac{2}{\sqrt\pi} \int_0^z e^{-t^2} dt$.

::: details
Real input is passed to `arb_hypgeom_erf`. Otherwise the function uses
the asymptotic expansion when its error term is below $2^{-\mathit{prec}}$
(with a reduced precision for the $U$ term), and otherwise the
${}_1F_1$ series better suited to the direction of $z$, at a working
precision increased by the $2 \min(x^2, y^2)/\log 2$ bits (for $z = x +
yi$) lost to cancellation. For inexact $z$ the series is evaluated at
the midpoint and the propagated error is bounded with
`acb_hypgeom_erf_propagated_error`.
:::
:::

::: function
void acb_hypgeom_erf_propagated_error(mag_t re, mag_t im, const acb_t z)

Sets `re` and `im` to upper bounds for the errors in the real and
imaginary parts that result from replacing $\operatorname{erf}(z)$ by
$\operatorname{erf}$ evaluated at the midpoint of $z$, using a bound for
the first derivative $\frac{2}{\sqrt\pi} e^{-z^2}$ on the ball.
:::

::: function
void _acb_hypgeom_erf_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_erf_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

Sets `res` to the error function of the power series $z$, truncated to
length `len`.
:::

::: function
void acb_hypgeom_erfc(acb_t res, const acb_t z, slong prec)

Sets `res` to the complementary error function
$\operatorname{erfc}(z) = 1 - \operatorname{erf}(z)$, avoiding the
catastrophic cancellation of the subtraction when $\operatorname{Re} z$
is large and positive.
:::

::: function
void _acb_hypgeom_erfc_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_erfc_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

Sets `res` to the complementary error function of the power series $z$,
truncated to length `len`.
:::

::: function
void acb_hypgeom_erfi(acb_t res, const acb_t z, slong prec)

Sets `res` to the imaginary error function $\operatorname{erfi}(z) =
-i\operatorname{erf}(iz)$, computed with `acb_hypgeom_erf`.
:::

::: function
void _acb_hypgeom_erfi_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_erfi_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

Sets `res` to the imaginary error function of the power series $z$,
truncated to length `len`.
:::

::: function
void acb_hypgeom_fresnel(acb_t res1, acb_t res2, const acb_t z, int normalized, slong prec)

Sets `res1` to the Fresnel sine integral $S(z)$ and `res2` to the
Fresnel cosine integral $C(z)$; either output may be `NULL`. If
`normalized` is 0, $S(z) = \int_0^z \sin(t^2) dt$ and $C(z) = \int_0^z
\cos(t^2) dt$; if `normalized` is 1, the integrands are $\sin(\tfrac12
\pi t^2)$ and $\cos(\tfrac12 \pi t^2)$, the convention of Abramowitz and
Stegun and of the DLMF. The values are computed from
$\operatorname{erf}((1 \pm i) w)$, where $w$ is $z$ times a real
constant, evaluated at the midpoint of $z$ with a bound for the
propagated error.
:::

::: function
void _acb_hypgeom_fresnel_series(acb_ptr res1, acb_ptr res2, acb_srcptr z, slong zlen, int normalized, slong len, slong prec)
void acb_hypgeom_fresnel_series(acb_poly_t res1, acb_poly_t res2, const acb_poly_t z, int normalized, slong len, slong prec)

Sets `res1` and `res2` to the Fresnel sine and cosine integrals of the
power series $z$, truncated to length `len`. Either output may be
`NULL`.
:::

## Incomplete gamma and beta functions {#incomplete}

::: function
void acb_hypgeom_gamma_upper_asymp(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)
void acb_hypgeom_gamma_upper_1f1a(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)
void acb_hypgeom_gamma_upper_1f1b(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)
void acb_hypgeom_gamma_upper_singular(acb_t res, slong s, const acb_t z, int regularized, slong prec)
void acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)

Sets `res` to the upper incomplete gamma function $\Gamma(s,z) =
\int_z^\infty t^{s-1} e^{-t} dt$ if `regularized` is 0, to $Q(s,z) =
\Gamma(s,z) / \Gamma(s)$ if `regularized` is 1, and to the generalized
exponential integral $z^{-s} \Gamma(s,z) = E_{1-s}(z)$ if `regularized`
is 2 (mainly for internal use; `acb_hypgeom_expint` is the interface
for $E_s$).

::: details
The methods implement, in this order,
$$
\Gamma(s,z) = e^{-z} U(1-s,1-s,z), \qquad
\Gamma(s,z) = \Gamma(s) - \frac{z^s}{s} {}_1F_1(s, s+1, -z),
$$
$$
\Gamma(s,z) = \Gamma(s) - \frac{z^s e^{-z}}{s} {}_1F_1(1, s+1, z),
$$
$$
\Gamma(-n,z) = \frac{(-1)^n}{n!} (\psi(n+1) - \log(z))
            + \frac{(-1)^n}{(n+1)!} z \, {}_2F_2(1,1;2,2+n;-z)
            - z^{-n} \sum_{k=0}^{n-1} \frac{(-z)^k}{(k-n) k!}
$$
for an integer $n = -s \ge 0$ (the singular case, where the second and
third formulas have removable singularities; it sums the finite sum
directly and assumes that $n$ is not too large). The automatic version
handles $z = 0$ and $s = 1, 2, 3$ by closed forms, returns 0 for $Q$ at
nonpositive integer $s$, uses the asymptotic expansion for large $|z|$
(with a region test for complex $s$ and $z$), the singular formula for
integer $s \le 0$, and otherwise the first ${}_1F_1$ form for
$\operatorname{Re} z \le 0$ and the second for $\operatorname{Re} z >
0$. Real arguments with an inaccurate result are retried with
`arb_hypgeom_gamma_upper_integration`.
:::
:::

::: function
void _acb_hypgeom_gamma_upper_series(acb_ptr res, const acb_t s, acb_srcptr z, slong zlen, int regularized, slong n, slong prec)
void acb_hypgeom_gamma_upper_series(acb_poly_t res, const acb_t s, const acb_poly_t z, int regularized, slong n, slong prec)

Sets `res` to the upper incomplete gamma function of the constant $s$
and the power series $z$, truncated to length $n$, with `regularized`
as in `acb_hypgeom_gamma_upper`.
:::

::: function
void acb_hypgeom_gamma_lower(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)

Sets `res` to the lower incomplete gamma function $\gamma(s,z) =
\int_0^z t^{s-1} e^{-t} dt = \frac{z^s}{s} {}_1F_1(s, s+1, -z)$ if
`regularized` is 0, to $P(s,z) = \gamma(s,z) / \Gamma(s)$ if
`regularized` is 1, and to $\gamma^{*}(s,z) = z^{-s} P(s,z)$, which is
entire in $s$ and $z$, if `regularized` is 2.
:::

::: function
void _acb_hypgeom_gamma_lower_series(acb_ptr res, const acb_t s, acb_srcptr z, slong zlen, int regularized, slong n, slong prec)
void acb_hypgeom_gamma_lower_series(acb_poly_t res, const acb_t s, const acb_poly_t z, int regularized, slong n, slong prec)

Sets `res` to the lower incomplete gamma function of the constant $s$
and the power series $z$, truncated to length $n$, with `regularized`
as in `acb_hypgeom_gamma_lower`.
:::

::: function
void acb_hypgeom_beta_lower(acb_t res, const acb_t a, const acb_t b, const acb_t z, int regularized, slong prec)

Sets `res` to the incomplete beta function $B(a,b;z) = \int_0^z
t^{a-1} (1-t)^{b-1} dt$, or the regularized incomplete beta function
$I(a,b;z) = B(a,b;z) / B(a,b;1)$ if `regularized` is nonzero.

::: details
In general the integral is interpreted by analytic continuation. The
precise definitions for all parameters are
$$
B(a,b;z) = \frac{z^a}{a} {}_2F_1(a, 1-b, a+1, z), \qquad
I(a,b;z) = \frac{\Gamma(a+b)}{\Gamma(b)} z^a\,
{}_2\mathbf{F}_1(a, 1-b, a+1, z).
$$
With these definitions both functions are undefined for nonpositive
integer $a$, and $I$ is undefined for nonpositive integer $a + b$.
:::
:::

::: function
void _acb_hypgeom_beta_lower_series(acb_ptr res, const acb_t a, const acb_t b, acb_srcptr z, slong zlen, int regularized, slong n, slong prec)
void acb_hypgeom_beta_lower_series(acb_poly_t res, const acb_t a, const acb_t b, const acb_poly_t z, int regularized, slong n, slong prec)

Sets `res` to $B(a,b;z)$ (or $I(a,b;z)$ if `regularized` is nonzero)
for constants $a$ and $b$ and the power series $z$, truncated to
length $n$.
:::

## Exponential and trigonometric integrals {#expint}

::: function
void acb_hypgeom_expint(acb_t res, const acb_t s, const acb_t z, slong prec)

Sets `res` to the generalized exponential integral $E_s(z) = z^{s-1}
\Gamma(1-s, z)$, computed with `acb_hypgeom_gamma_upper`.
:::

::: function
void acb_hypgeom_ei_asymp(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_ei_2f2(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_ei(acb_t res, const acb_t z, slong prec)

Sets `res` to the exponential integral $\operatorname{Ei}(z)$, computed
respectively from
$$
\operatorname{Ei}(z) = -e^z U(1,1,-z) - \log(-z)
    + \frac{1}{2} \left(\log(z) - \log\left(\frac{1}{z}\right) \right),
$$
$$
\operatorname{Ei}(z) = z\, {}_2F_2(1, 1; 2, 2; z) + \gamma
    + \frac{1}{2} \left(\log(z) - \log\left(\frac{1}{z}\right) \right),
$$
or by an automatic choice based on the expected accuracy of the two
methods. The result is indeterminate if $z$ contains zero.
:::

::: function
void _acb_hypgeom_ei_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_ei_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

Sets `res` to the exponential integral of the power series $z$,
truncated to length `len`.
:::

::: function
void acb_hypgeom_si_asymp(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_si_1f2(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_si(acb_t res, const acb_t z, slong prec)

Sets `res` to the sine integral $\operatorname{Si}(z)$, computed
respectively from
$$
\operatorname{Si}(z) = \frac{i}{2} \left[
    e^{iz} U(1,1,-iz) - e^{-iz} U(1,1,iz) +
    \log(-iz) - \log(iz) \right],
$$
$$
\operatorname{Si}(z) = z\, {}_1F_2(\tfrac{1}{2}; \tfrac{3}{2}, \tfrac{3}{2}; -\tfrac{z^2}{4}),
$$
or by an automatic choice. Real input is passed to `arb_hypgeom_si`.
:::

::: function
void _acb_hypgeom_si_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_si_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

Sets `res` to the sine integral of the power series $z$, truncated to
length `len`.
:::

::: function
void acb_hypgeom_ci_asymp(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_ci_2f3(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_ci(acb_t res, const acb_t z, slong prec)

Sets `res` to the cosine integral $\operatorname{Ci}(z)$, computed
respectively from
$$
\operatorname{Ci}(z) = \log(z) - \frac{1}{2} \left[
    e^{iz} U(1,1,-iz) + e^{-iz} U(1,1,iz) +
    \log(-iz) + \log(iz) \right],
$$
$$
\operatorname{Ci}(z) = -\tfrac{z^2}{4}\,
    {}_2F_3(1, 1; 2, 2, \tfrac{3}{2}; -\tfrac{z^2}{4})
    + \log(z) + \gamma,
$$
or by an automatic choice. Positive real input is passed to
`arb_hypgeom_ci`.
:::

::: function
void _acb_hypgeom_ci_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_ci_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

Sets `res` to the cosine integral of the power series $z$, truncated to
length `len`.
:::

::: function
void acb_hypgeom_shi(acb_t res, const acb_t z, slong prec)

Sets `res` to the hyperbolic sine integral $\operatorname{Shi}(z) = -i
\operatorname{Si}(iz)$, computed with `acb_hypgeom_si`.
:::

::: function
void _acb_hypgeom_shi_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_shi_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

Sets `res` to the hyperbolic sine integral of the power series $z$,
truncated to length `len`.
:::

::: function
void acb_hypgeom_chi_asymp(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_chi_2f3(acb_t res, const acb_t z, slong prec)
void acb_hypgeom_chi(acb_t res, const acb_t z, slong prec)

Sets `res` to the hyperbolic cosine integral $\operatorname{Chi}(z)$,
computed respectively from
$$
\operatorname{Chi}(z) = -\frac{1}{2} \left[
    e^{z} U(1,1,-z) + e^{-z} U(1,1,z) +
    \log(-z) - \log(z) \right],
$$
$$
\operatorname{Chi}(z) = \tfrac{z^2}{4}\,
    {}_2F_3(1, 1; 2, 2, \tfrac{3}{2}; \tfrac{z^2}{4})
    + \log(z) + \gamma,
$$
or by an automatic choice.
:::

::: function
void _acb_hypgeom_chi_series(acb_ptr res, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_chi_series(acb_poly_t res, const acb_poly_t z, slong len, slong prec)

Sets `res` to the hyperbolic cosine integral of the power series $z$,
truncated to length `len`.
:::

::: function
void acb_hypgeom_li(acb_t res, const acb_t z, int offset, slong prec)

Sets `res` to the logarithmic integral $\operatorname{li}(z) =
\operatorname{Ei}(\log(z))$ if `offset` is zero, or the offset
logarithmic integral $\operatorname{Li}(z) = \operatorname{li}(z) -
\operatorname{li}(2)$ if `offset` is nonzero.
:::

::: function
void _acb_hypgeom_li_series(acb_ptr res, acb_srcptr z, slong zlen, int offset, slong len, slong prec)
void acb_hypgeom_li_series(acb_poly_t res, const acb_poly_t z, int offset, slong len, slong prec)

Sets `res` to the logarithmic integral (the offset version if `offset`
is nonzero) of the power series $z$, truncated to length `len`.
:::

## Bessel functions {#bessel}

The methods for the Bessel, Airy and Coulomb wave functions are
described in [](algorithms/special-functions/bessel-airy).

::: function
void acb_hypgeom_bessel_j_asymp(acb_t res, const acb_t nu, const acb_t z, slong prec)

Sets `res` to the Bessel function of the first kind $J_\nu(z)$
computed with the asymptotic expansion of `acb_hypgeom_u_asymp`.

::: details
For all complex $\nu, z$,
$$
J_{\nu}(z) = \frac{z^{\nu}}{2^{\nu} e^{iz} \Gamma(\nu+1)}
    {}_1F_1(\nu+\tfrac{1}{2}, 2\nu+1, 2iz) = A_{+} B_{+} + A_{-} B_{-}
$$
where
$$
A_{\pm} = z^{\nu} (z^2)^{-\tfrac{1}{2}-\nu} (\mp i z)^{\tfrac{1}{2}+\nu}
(2 \pi)^{-1/2} = (\pm iz)^{-1/2-\nu} z^{\nu} (2 \pi)^{-1/2}, \qquad
B_{\pm} = e^{\pm i z} U^{*}(\nu+\tfrac{1}{2}, 2\nu+1, \mp 2iz).
$$
Simpler forms of $A_{\pm}$ hold in special cases and are used when they
apply: if $\nu + \tfrac{1}{2} = n \in \Z$, then $A_{\pm} = (\pm i)^{n}
(2 \pi z)^{-1/2}$, and if $\operatorname{Re}(z) > 0$, then $A_{\pm} =
\exp(\mp i [(2\nu+1)/4] \pi) (2 \pi z)^{-1/2}$.
:::
:::

::: function
void acb_hypgeom_bessel_j_0f1(acb_t res, const acb_t nu, const acb_t z, slong prec)

Sets `res` to $J_\nu(z)$ computed from the convergent series
$$
J_{\nu}(z) = \frac{1}{\Gamma(\nu+1)} \left(\frac{z}{2}\right)^{\nu}
    {}_0F_1\left(\nu+1, -\frac{z^2}{4}\right),
$$
using $J_{-n}(z) = (-1)^n J_n(z)$ for negative integer $\nu$.
:::

::: function
void acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, slong prec)

Sets `res` to the Bessel function of the first kind $J_\nu(z)$.

::: details
For $|z| < 8$ the function sums the convergent series. Otherwise it
uses the asymptotic expansion if $|z| \ge 2^{64}$ or $2|z| \ge
\mathit{prec}$, and the convergent series otherwise, with the working
precision increased by an estimate $(|z| - |\operatorname{Im} z|)/\log 2$
of the cancellation (no increase when the terms have equal signs or
$\nu > |z|^2/4$). For inexact input with integer $\nu$ the series is
evaluated at the midpoint and the propagated error is bounded with
`acb_hypgeom_bessel_j_deriv_bound`. If the enclosure has less than one
bit of accuracy, it is intersected with a direct bound for
$|J_\nu(z)|$.
:::
:::

::: function
void acb_hypgeom_bessel_j_deriv_bound(mag_t res, const acb_t nu, const acb_t z, ulong d)

Sets `res` to an upper bound, possibly crude, for $|J_\nu^{(d)}(z)|$
($|J_\nu(z)|$ for $d = 0$). Implemented only for integer $\nu$, using
for example $|J_n^{(d)}(z)| \le e^{|\operatorname{Im} z|}$; for
non-integer $\nu$ the bound is infinite.
:::

::: function
void _acb_hypgeom_bessel_j_series(acb_ptr res, const acb_t nu, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_bessel_j_series(acb_poly_t res, const acb_t nu, const acb_poly_t z, slong len, slong prec)

Sets `res` to $J_\nu(z)$ for the constant $\nu$ and the power series
$z$, truncated to length `len`. Both functions call the generic
`_gr_poly_bessel_j_series` with a complex ball context; as described
for `arb_hypgeom_bessel_j_series`, the result is indeterminate if the
constant term of $z$ is an inexact ball containing zero.
:::

::: function
void acb_hypgeom_bessel_y(acb_t res, const acb_t nu, const acb_t z, slong prec)

Sets `res` to the Bessel function of the second kind
$$
Y_{\nu}(z) = \frac{\cos(\nu \pi) J_{\nu}(z) - J_{-\nu}(z)}{\sin(\nu \pi)},
$$
or for $\nu = n$ an integer the limit
$$
Y_n(z) = -\frac{2}{\pi} \left( i^n K_n(iz) +
    \left[\log(iz)-\log(z)\right] J_n(z) \right).
$$
The output is indeterminate if $\nu$ is inexact and contains an
integer.
:::

::: function
void acb_hypgeom_bessel_jy(acb_t res1, acb_t res2, const acb_t nu, const acb_t z, slong prec)

Sets `res1` to $J_\nu(z)$ and `res2` to $Y_\nu(z)$, computed together.
The Hankel functions $H_\nu^{(1)}(z), H_\nu^{(2)}(z) = J_\nu(z) \pm i
Y_\nu(z)$ can be formed from them.
:::

## Modified Bessel functions {#modified-bessel}

::: function
void acb_hypgeom_bessel_i_asymp(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
void acb_hypgeom_bessel_i_0f1(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)
void acb_hypgeom_bessel_i(acb_t res, const acb_t nu, const acb_t z, slong prec)
void acb_hypgeom_bessel_i_scaled(acb_t res, const acb_t nu, const acb_t z, slong prec)

Sets `res` to the modified Bessel function of the first kind $I_\nu(z)
= z^\nu (iz)^{-\nu} J_\nu(iz)$, computed respectively with the
asymptotic expansion (see `acb_hypgeom_bessel_j_asymp`), the convergent
series
$$
I_{\nu}(z) = \frac{1}{\Gamma(\nu+1)} \left(\frac{z}{2}\right)^{\nu}
             {}_0F_1\left(\nu+1, \frac{z^2}{4}\right),
$$
or an automatic choice. `acb_hypgeom_bessel_i_scaled` computes
$e^{-z} I_\nu(z)$; the `asymp` and `0f1` versions compute the scaled
function if `scaled` is nonzero. For real arguments with an inaccurate
result, the automatic versions also try
`arb_hypgeom_bessel_i_integration`.
:::

::: function
void acb_hypgeom_bessel_k_asymp(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)

Sets `res` to the modified Bessel function of the second kind
$$
K_{\nu}(z) = \left(\frac{2z}{\pi}\right)^{-1/2} e^{-z}
    U^{*}(\nu+\tfrac{1}{2}, 2\nu+1, 2z),
$$
valid for all $\nu$ and $z \ne 0$, using `acb_hypgeom_u_asymp`; or
$e^z K_\nu(z)$ if `scaled` is nonzero.
:::

::: function
void acb_hypgeom_bessel_k_0f1_series(acb_poly_t res, const acb_poly_t nu, const acb_poly_t z, int scaled, slong len, slong prec)

Sets `res` to $K_\nu(z)$ (or $e^z K_\nu(z)$ if `scaled` is nonzero) as
a power series truncated to length `len`, for power series $\nu, z \in
\mathbb{C}[[x]]$, using
$$
K_{\nu}(z) = \frac{1}{2} \frac{\pi}{\sin(\pi \nu)} \left[
    \left(\frac{z}{2}\right)^{-\nu}
        {}_0{\widetilde F}_1\left(1-\nu, \frac{z^2}{4}\right) -
    \left(\frac{z}{2}\right)^{\nu}
        {}_0{\widetilde F}_1\left(1+\nu, \frac{z^2}{4}\right) \right].
$$
If the constant term of $\nu$ is an integer, one extra derivative is
computed and the singularity is removed (the coefficient of $x$ in
$\nu$ must then be nonzero). The output is indeterminate if the
constant term of $\nu$ is inexact and contains an integer.
:::

::: function
void acb_hypgeom_bessel_k_0f1(acb_t res, const acb_t nu, const acb_t z, int scaled, slong prec)

Sets `res` to $K_\nu(z)$ (or $e^z K_\nu(z)$ if `scaled` is nonzero)
computed from
$$
K_{\nu}(z) = \frac{1}{2} \left[
    \left(\frac{z}{2}\right)^{-\nu} \Gamma(\nu)
        {}_0F_1\left(1-\nu, \frac{z^2}{4}\right) -
    \left(\frac{z}{2}\right)^{\nu}
        \frac{\pi}{\nu \sin(\pi \nu) \Gamma(\nu)}
        {}_0F_1\left(\nu+1, \frac{z^2}{4}\right) \right]
$$
for $\nu \notin \Z$, and the limit computed with
`acb_hypgeom_bessel_k_0f1_series` for $\nu \in \Z$. The output is
indeterminate if $\nu$ is inexact and contains an integer.
:::

::: function
void acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, slong prec)
void acb_hypgeom_bessel_k_scaled(acb_t res, const acb_t nu, const acb_t z, slong prec)

Sets `res` to the modified Bessel function of the second kind
$K_\nu(z)$, or to $e^z K_\nu(z)$, choosing between the asymptotic
expansion and the convergent series. For real arguments with an
inaccurate result they also try `arb_hypgeom_bessel_k_integration`.
:::

## Airy functions {#airy}

The Airy functions $\operatorname{Ai}(z)$ and $\operatorname{Bi}(z)$
are linearly independent solutions of $y'' - zy = 0$; both are entire.
For negative $z$ both oscillate; for positive $z$ the first decreases
and the second increases exponentially. They can be expressed with
Bessel functions of order $\pm 1/3$, but these formulas hold only
piecewise (because of the Stokes phenomenon), so the module provides a
dedicated implementation. The functions below compute
$(\operatorname{Ai}(z), \operatorname{Ai}'(z), \operatorname{Bi}(z),
\operatorname{Bi}'(z))$ together; any of the four outputs may be `NULL`,
which saves work.

::: function
void acb_hypgeom_airy_direct(acb_t ai, acb_t ai_prime, acb_t bi, acb_t bi_prime, const acb_t z, slong n, slong prec)

Computes the Airy functions with the convergent series at the origin
truncated after $n$ terms, including a bound for the truncation error.
:::

::: function
void acb_hypgeom_airy_asymp(acb_t ai, acb_t ai_prime, acb_t bi, acb_t bi_prime, const acb_t z, slong n, slong prec)

Computes the Airy functions with the asymptotic expansions truncated
after $n$ terms, including a bound for the truncation error. The
bounds and connection formulas are given in
[](algorithms/special-functions/bessel-airy#airy-asymp).
:::

::: function
void acb_hypgeom_airy_bound(mag_t ai, mag_t ai_prime, mag_t bi, mag_t bi_prime, const acb_t z)

Sets the outputs to upper bounds for $|\operatorname{Ai}(z)|$,
$|\operatorname{Ai}'(z)|$, $|\operatorname{Bi}(z)|$,
$|\operatorname{Bi}'(z)|$ on the ball $z$, from first-order asymptotic
expansions with error bounds. This is faster than calling
`acb_hypgeom_airy_asymp` with $n = 1$.
:::

::: function
void acb_hypgeom_airy(acb_t ai, acb_t ai_prime, acb_t bi, acb_t bi_prime, const acb_t z, slong prec)

Computes the Airy functions with an automatic choice of method.

::: details
The asymptotic expansion is used whenever it reaches full accuracy
($|z| \ge 4$ and enough terms decrease), and the convergent series
otherwise. In the latter case the working precision is estimated in
double precision by comparing the leading asymptotic estimate of the
function values with the size of the largest term of the series (the
estimate ignores the cancellation near the zeros of the functions).
The series is evaluated at the midpoint of $z$, and the propagated
error is bounded with the derivative bounds of `acb_hypgeom_airy_bound`
together with $\operatorname{Ai}''(z) = z \operatorname{Ai}(z)$.
:::
:::

::: function
void acb_hypgeom_airy_jet(acb_ptr ai, acb_ptr bi, const acb_t z, slong len, slong prec)

Sets the vectors `ai` and `bi` to the Taylor coefficients of
$\operatorname{Ai}$ and $\operatorname{Bi}$ at the point $z$, truncated
to length `len`. Either output may be `NULL`; $z$ must not alias the
outputs. The expansions of the derivatives are obtained by computing
one more coefficient and applying `_acb_poly_derivative`.
:::

::: function
void _acb_hypgeom_airy_series(acb_ptr ai, acb_ptr ai_prime, acb_ptr bi, acb_ptr bi_prime, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_airy_series(acb_poly_t ai, acb_poly_t ai_prime, acb_poly_t bi, acb_poly_t bi_prime, const acb_poly_t z, slong len, slong prec)

Sets the outputs to the Airy functions and their derivatives evaluated
at the power series $z$, truncated to length `len`. Any of the outputs
may be `NULL`.
:::

## Coulomb wave functions {#coulomb}

Coulomb wave functions are solutions of the Coulomb wave equation
$$
y'' + \left(1 - \frac{2 \eta}{z} - \frac{\ell(\ell+1)}{z^2}\right) y = 0,
$$
the radial Schrödinger equation for a charged particle in a Coulomb
potential $1/z$, where $\ell$ is the orbital angular momentum and
$\eta$ the Sommerfeld parameter. The standard solutions are
$F_\ell(\eta,z)$, regular at the origin, and $G_\ell(\eta,z)$,
irregular at the origin; the irregular solutions $H^{\pm}_\ell(\eta,z)
= G_\ell(\eta,z) \pm i F_\ell(\eta,z)$ are also used. They are special
cases of confluent hypergeometric functions; the normalizations and
connection formulas are discussed in [@DYF1999], [@Gas2018], [@Mic2007]
and [@NIST2012, chapter 33]. Here the analytic continuations of all
functions are defined so that the branch cut with respect to $z$ is
the negative real axis; precise definitions are given in
[@Fungrim] (topic "Coulomb wave functions"), whose formulas the
implementation follows.

The functions below compute $F_\ell(\eta,z)$, $G_\ell(\eta,z)$,
$H^{+}_\ell(\eta,z)$ and $H^{-}_\ell(\eta,z)$ together; any output may
be `NULL`. The redundant functions $H^{\pm}$ are provided because the
linear combination of $F$ and $G$ suffers from cancellation in parts of
the complex plane.

::: function
void acb_hypgeom_coulomb(acb_t F, acb_t G, acb_t Hpos, acb_t Hneg, const acb_t l, const acb_t eta, const acb_t z, slong prec)

Sets `F`, `G`, `Hpos`, `Hneg` to the respective Coulomb wave functions.
Any of the outputs may be `NULL`.

::: details
All four are computed from $U(1 + \ell \pm i\eta, 2\ell + 2, \mp 2iz)$
(by the asymptotic expansion when $|z|$ is large, by `acb_hypgeom_u`
otherwise) and the log-gamma values $\log\Gamma(1 + \ell \pm i\eta)$,
with formulas chosen according to the sign of $\operatorname{Re} z$ and
$\operatorname{Im} z$ so that each output is continuous across the
real axis except on the branch cut.
:::
:::

::: function
void acb_hypgeom_coulomb_jet(acb_ptr F, acb_ptr G, acb_ptr Hpos, acb_ptr Hneg, const acb_t l, const acb_t eta, const acb_t z, slong len, slong prec)

Sets `F`, `G`, `Hpos`, `Hneg` to the Taylor coefficients of the Coulomb
wave functions at the point $z$, truncated to length `len`, using the
differential equation to obtain the higher coefficients. Any of the
outputs may be `NULL`.
:::

::: function
void _acb_hypgeom_coulomb_series(acb_ptr F, acb_ptr G, acb_ptr Hpos, acb_ptr Hneg, const acb_t l, const acb_t eta, acb_srcptr z, slong zlen, slong len, slong prec)
void acb_hypgeom_coulomb_series(acb_poly_t F, acb_poly_t G, acb_poly_t Hpos, acb_poly_t Hneg, const acb_t l, const acb_t eta, const acb_poly_t z, slong len, slong prec)

Sets the outputs to the Coulomb wave functions evaluated at the power
series $z$, truncated to length `len`. Any of the outputs may be
`NULL`.
:::

## Gauss hypergeometric function {#gauss}

The functions in this section compute the Gauss hypergeometric function
$$
F(z) = {}_2F_1(a,b,c,z) = \sum_{k=0}^{\infty} \frac{(a)_k (b)_k}{(c)_k} \frac{z^k}{k!}
$$
or, if the flag `regularized` (`ACB_HYPGEOM_2F1_REGULARIZED`) is set,
the regularized function $\mathbf{F}(z) = \mathbf{F}(a,b,c,z) =
{}_2F_1(a,b,c,z) / \Gamma(c)$, which is entire in $a$, $b$ and $c$. The
branch cut is $[1, +\infty)$; on the cut, the value is the limit from
the lower half-plane. The algorithms are described in
[](algorithms/special-functions/gauss-2f1).

::: macro
ACB_HYPGEOM_2F1_REGULARIZED
ACB_HYPGEOM_2F1_AB
ACB_HYPGEOM_2F1_AC
ACB_HYPGEOM_2F1_BC
ACB_HYPGEOM_2F1_ABC

Bit flags for the `flags` argument of `acb_hypgeom_2f1` and
`acb_hypgeom_2f1_transform`, with the values 1, 2, 4, 8 and 16:
compute the regularized function, and declare that $a - b$, $a - c$,
$b - c$ or $a + b - c$ respectively is an exact integer.
:::

::: function
void acb_hypgeom_2f1(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int flags, slong prec)

Sets `res` to $F(z)$, or $\mathbf{F}(z)$ if `flags` contains
`ACB_HYPGEOM_2F1_REGULARIZED`, with an automatic choice of algorithm.
The other flags declare that a parameter difference is an exact integer
even if the input balls are inexact, so that the correct limit is
computed where a transformation formula is singular.

::: details
The flags `ACB_HYPGEOM_2F1_AB` ($a - b \in \Z$) and
`ACB_HYPGEOM_2F1_ABC` ($a + b - c \in \Z$) are used in this way; `AC`
and `BC` are accepted but currently unused. For example, to evaluate
${}_2F_1(\sqrt{2}, 1/2, \sqrt{2}+3/2, 9/10)$, set `ACB_HYPGEOM_2F1_ABC`;
otherwise the transformation $z \to 1 - z$ divides by
$\sin(\pi(c-a-b))$ computed on a ball containing zero and the result is
indeterminate. If the parameters are exact (including exact integers
and half-integers), the limits are detected automatically and the flags
are unnecessary.

The function handles special parameters first (terminating series,
$z = 0$ and $z = 1$, $a = c$ or $b = c$ for $\mathbf{F}$, reduction to a
polynomial by the Pfaff transformation when $c - a$ or $c - b$ is a
nonpositive integer), then chooses with `acb_hypgeom_2f1_choose` between
direct summation, one of five transformations and the corner case, and
retries real arguments with `arb_hypgeom_2f1_integration` if the result
is inaccurate. See [](acb_hypgeom/2f1) and
[](algorithms/special-functions/gauss-2f1).
:::
:::

::: function
void acb_hypgeom_2f1_direct(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, slong prec)

Sets `res` to $F(z)$ (or $\mathbf{F}(z)$) by direct summation of the
hypergeometric series, which converges for $|z| < 1$.
:::

::: function
void acb_hypgeom_2f1_series_direct(acb_poly_t res, const acb_poly_t a, const acb_poly_t b, const acb_poly_t c, const acb_poly_t z, int regularized, slong len, slong prec)

Sets `res` to $F(z)$ (or $\mathbf{F}(z)$) for power series parameters
and argument, truncated to length `len`, by direct summation of the
hypergeometric series.
:::

::: function
void acb_hypgeom_2f1_transform(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int flags, int which, slong prec)
void acb_hypgeom_2f1_transform_limit(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, int which, slong prec)

Sets `res` to $F(z)$ (or $\mathbf{F}(z)$) computed with the linear
fractional transformation selected by `which`: 1 for $z/(z-1)$, 2 for
$1/z$, 3 for $1/(1-z)$, 4 for $1-z$ and 5 for $1-1/z$.

::: details
Transformation 1 (Pfaff) maps the function to one ${}_2F_1$ series;
the others express it as a combination of two series with gamma
function factors [@NIST2012, §15.8], which is singular when $b - a$
(for 2 and 3) or $c - a - b$ (for 4 and 5) is an integer.
`acb_hypgeom_2f1_transform` detects this case from the exactness of the
parameters or from `flags` (with the meaning of `acb_hypgeom_2f1`) and
then calls `acb_hypgeom_2f1_transform_limit`, which assumes that
`which` is not 1, that $b - a$ is an exact integer if `which` is 2 or 3
and that $c - a - b$ is an exact integer if `which` is 4 or 5, and
computes the limit by differentiating with respect to a perturbation of
$a$ (power series of length 2). The result is indeterminate if $z$
contains 0 or 1.
:::
:::

::: function
void acb_hypgeom_2f1_corner(acb_t res, const acb_t a, const acb_t b, const acb_t c, const acb_t z, int regularized, slong prec)

Sets `res` to $F(z)$ (or $\mathbf{F}(z)$) for $z$ near the points
$\exp(\pm \pi i/3)$, where none of the transformations gives a
convergent series, by analytic continuation along a path from a point
where the direct series converges.
:::

::: function
void acb_hypgeom_2f1_continuation(acb_t res0, acb_t res1, const acb_t a, const acb_t b, const acb_t c, const acb_t z0, const acb_t z1, const acb_t f0, const acb_t f1, slong prec)

Given $F(z_0)$ and $F'(z_0)$ in `f0` and `f1`, sets `res0` and `res1`
to $F(z_1)$ and $F'(z_1)$, computed by integrating the hypergeometric
differential equation along the straight segment from $z_0$ to $z_1$
with Taylor series and rigorous bounds. The points should be well
isolated from the singular points 0 and 1. If $F$ denotes the
regularized function, $f_0, f_1$ must be values of the regularized
function.
:::

::: function
int acb_hypgeom_2f1_choose(const acb_t z)

Chooses the method for computing ${}_2F_1$ from the location of $z$:
returns 0 for direct summation, 1 to 5 for the transformation with that
index in `acb_hypgeom_2f1_transform`, or 6 for the corner case.

::: details
Direct summation is chosen if $|z| \le 1/4$, the Pfaff transformation
if $|z/(z-1)| \le 1/4$, the smaller of the two if either is at most
$3/4$, and otherwise the transformation that minimizes the modulus of
the new argument, if that modulus is at most $3/4$; the remaining
region, around $\exp(\pm\pi i/3)$, is the corner case.
:::
:::

## Orthogonal polynomials and functions {#orthogonal}

The algorithms for this section are described in
[](algorithms/special-functions/legendre).

::: function
void acb_hypgeom_chebyshev_t(acb_t res, const acb_t n, const acb_t z, slong prec)
void acb_hypgeom_chebyshev_u(acb_t res, const acb_t n, const acb_t z, slong prec)

Sets `res` to the Chebyshev polynomial or Chebyshev function of the
first or second kind
$$
T_n(z) = {}_2F_1\left(-n,n,\frac{1}{2},\frac{1-z}{2}\right), \qquad
U_n(z) = (n+1)\, {}_2F_1\left(-n,n+2,\frac{3}{2},\frac{1-z}{2}\right).
$$

::: details
For word-size integer $n$ the functions call `acb_chebyshev_t_ui` and
`acb_chebyshev_u_ui`. Otherwise the hypergeometric series is used near
$z = 1$, and the trigonometric representations $T_n(\cos \theta) =
\cos(n\theta)$ and $U_n(\cos\theta) = \sin((n+1)\theta)/\sin\theta$
elsewhere.
:::
:::

::: function
void acb_hypgeom_jacobi_p(acb_t res, const acb_t n, const acb_t a, const acb_t b, const acb_t z, slong prec)

Sets `res` to the Jacobi polynomial (or Jacobi function)
$$
P_n^{(a,b)}(z) = \frac{(a+1)_n}{\Gamma(n+1)}\,
{}_2F_1\left(-n,n+a+b+1,a+1,\frac{1-z}{2}\right).
$$
For nonnegative integer $n$ this is a polynomial in $a$, $b$ and $z$,
also for parameters where the hypergeometric series is undefined; it is
then evaluated by direct methods.
:::

::: function
void acb_hypgeom_gegenbauer_c(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec)

Sets `res` to the Gegenbauer polynomial (or Gegenbauer function)
$$
C_n^{m}(z) = \frac{(2m)_n}{\Gamma(n+1)}\,
{}_2F_1\left(-n,2m+n,m+\frac{1}{2},\frac{1-z}{2}\right).
$$
For nonnegative integer $n$ this is a polynomial in $m$ and $z$, also
where the hypergeometric series is undefined; it is then evaluated by
direct methods.
:::

::: function
void acb_hypgeom_laguerre_l(acb_t res, const acb_t n, const acb_t m, const acb_t z, slong prec)

Sets `res` to the Laguerre polynomial (or Laguerre function)
$$
L_n^{m}(z) = \frac{(m+1)_n}{\Gamma(n+1)}\, {}_1F_1\left(-n,m+1,z\right).
$$
For nonnegative integer $n$ this is a polynomial in $m$ and $z$, also
where the hypergeometric series is undefined; it is then evaluated by
direct methods.

::: details
There are at least two incompatible extensions to negative integer
$n$: for $m = 0$ one can define $L_{-n}^0(z) = e^z L_{n-1}^0(-z)$, or
one can use the recurrence $L_{n-1}^m(z) + L_n^{m-1}(z) = L_n^m(z)$.
This case is left undefined: the result is indeterminate.
:::
:::

::: function
void acb_hypgeom_hermite_h(acb_t res, const acb_t n, const acb_t z, slong prec)

Sets `res` to the Hermite polynomial (or Hermite function)
$$
H_n(z) = 2^n \sqrt{\pi} \left(
    \frac{1}{\Gamma((1-n)/2)}\, {}_1F_1\left(-\frac{n}{2},\frac{1}{2},z^2\right)
    - \frac{2z}{\Gamma(-n/2)}\, {}_1F_1\left(\frac{1-n}{2},\frac{3}{2},z^2\right)\right).
$$
:::

::: function
void acb_hypgeom_legendre_p(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, slong prec)

Sets `res` to the associated Legendre function of the first kind of
degree $n$, order $m$ and argument $z$; for $m = 0$ this is the
Legendre polynomial $P_n(z)$. Several branch cut conventions are in
use. If `type` is 0, the function is
$$
P_n^m(z) = \frac{(1+z)^{m/2}}{(1-z)^{m/2}}\,
    \mathbf{F}\left(-n, n+1, 1-m, \frac{1-z}{2}\right),
$$
with cuts $(-\infty, -1]$ and $[1, \infty)$, and if `type` is 1, the
function is
$$
{\mathcal P}_n^m(z) = \frac{(z+1)^{m/2}}{(z-1)^{m/2}}\,
    \mathbf{F}\left(-n, n+1, 1-m, \frac{1-z}{2}\right),
$$
with the cut $(-\infty, 1]$. Type 0 and type 1 correspond to type 2 and
type 3 in *Mathematica* and *mpmath*; other values of `type` raise an
exception.

::: details
For real $z$, $m = 0$ and a nonnegative word-size integer $n$, the
function calls `arb_hypgeom_legendre_p_ui`. For integers $0 \le m \le
n$ and type 0 it evaluates $(1-z^2)^{m/2}$ times the polynomial of
`acb_hypgeom_legendre_p_uiui_rec` at low precision, near $z = 1$, or
for small $n - m$; if $n - m > 2\,\mathit{prec}$ the result is
indeterminate, since the recurrence would lose all accuracy. Otherwise
the hypergeometric representation is evaluated with `acb_hypgeom_2f1`.
:::
:::

::: function
void acb_hypgeom_legendre_q(acb_t res, const acb_t n, const acb_t m, const acb_t z, int type, slong prec)

Sets `res` to the associated Legendre function of the second kind of
degree $n$, order $m$ and argument $z$; for $m = 0$ this is the
Legendre function $Q_n(z)$. If `type` is 0, the function is
$$
Q_n^m(z) = \frac{\pi}{2 \sin(\pi m)}
    \left( \cos(\pi m) P_n^m(z) -
    \frac{\Gamma(1+m+n)}{\Gamma(1-m+n)} P_n^{-m}(z)\right),
$$
and if `type` is 1,
$$
\mathcal{Q}_n^m(z) = \frac{\pi}{2 \sin(\pi m)} e^{\pi i m}
    \left( \mathcal{P}_n^m(z) -
    \frac{\Gamma(1+m+n)}{\Gamma(1-m+n)} \mathcal{P}_n^{-m}(z)\right).
$$
Type 0 and type 1 correspond to type 2 and type 3 in *Mathematica* and
*mpmath*; other values of `type` raise an exception.

::: details
For integer $m$ either expression is interpreted as a limit. The
connection formulas [@WQ3a], [@WQ3b] and [@WQ3c] make it possible to
compute the function also in the limiting case. (The formula [@WQ3d]
would be useful too, but it is incorrect in the lower half-plane.)
:::
:::

::: function
void acb_hypgeom_legendre_p_uiui_rec(acb_t res, ulong n, ulong m, const acb_t z, slong prec)

For nonnegative integers $n$ and $m$, sets `res` to $(1-z^2)^{-m/2}
P_n^m(z)$, which is a polynomial in $z$, computed with the three-term
recurrence in the degree starting from $P_m^m$. The result is zero for
$m > n$, and indeterminate if $(n - m)/4 > \mathit{prec}$, where the
recurrence would lose all accuracy.
:::

::: function
void acb_hypgeom_spherical_y(acb_t res, slong n, slong m, const acb_t theta, const acb_t phi, slong prec)

Sets `res` to the spherical harmonic of degree $n$ and order $m$ at
latitude angle $\theta$ and longitude angle $\phi$, normalized so that
$$
Y_n^m(\theta, \phi) = \sqrt{\frac{2n+1}{4\pi} \frac{(n-m)!}{(n+m)!}}\,
e^{im\phi} P_n^m(\cos(\theta)).
$$
The definition is extended to negative $m$ and $n$ by symmetry. The
function is a polynomial in $\cos\theta$ and $\sin\theta$ and is
evaluated with `acb_hypgeom_legendre_p_uiui_rec`.
:::

## Dilogarithm {#dilog}

The dilogarithm is $\operatorname{Li}_2(z) = -\int_0^z
\frac{\log(1-t)}{t} dt = z\, {}_3F_2(1,1,1;2,2;z)$, with the branch cut
$[1, +\infty)$ (the value on the cut is the limit from the lower
half-plane). The algorithms are described in
[](algorithms/special-functions/dilogarithm).

::: function
void acb_hypgeom_dilog_bernoulli(acb_t res, const acb_t z, slong prec)

Sets `res` to $\operatorname{Li}_2(z)$ computed with the series
expansion in $w = \log(z)$, which converges like $|w/(2\pi)|^n$. This
gives good convergence near $z = e^{\pm i \pi / 3}$, where the
hypergeometric series fail. Since the coefficients involve Bernoulli
numbers, the method is meant for moderate precision.
:::

::: function
void acb_hypgeom_dilog_zero_taylor(acb_t res, const acb_t z, slong prec)

Sets `res` to $\operatorname{Li}_2(z)$ computed with the hypergeometric
series at 0, which is efficient only for $|z| \ll 1$.
:::

::: function
void acb_hypgeom_dilog_zero(acb_t res, const acb_t z, slong prec)

Sets `res` to $\operatorname{Li}_2(z)$ for $z$ close to 0, using the
Taylor series below 40000 bits (or for tiny $|z|$) and the bit-burst
algorithm otherwise.
:::

::: function
void acb_hypgeom_dilog_transform(acb_t res, const acb_t z, int algorithm, slong prec)

Sets `res` to $\operatorname{Li}_2(z)$ computed by one of the
transformations $1/z$, $1-z$, $z/(z-1)$, $1/(1-z)$, selected by
`algorithm` from 1 to 4, followed by `acb_hypgeom_dilog_zero` at the
transformed argument. For `algorithm` between 5 and 7, it instead
starts from the respective point $\pm i$, $(1\pm i)/2$ or $1 \pm i$
(with the sign of the imaginary part of $z$), where $\operatorname{Li}_2$
has a closed form in terms of $\pi$, $\log 2$ and Catalan's constant,
and continues to $z$ with the bit-burst method; these require that the
imaginary part of $z$ does not contain 0. The result is indeterminate
if the transformation is not valid on the whole input ball.
:::

::: function
void acb_hypgeom_dilog_continuation(acb_t res, const acb_t a, const acb_t z, slong prec)

Sets `res` to $\operatorname{Li}_2(z) - \operatorname{Li}_2(a)$
computed with the Taylor expansion at $a$, summed by binary splitting.
Both $a$ and $z$ should be well isolated from 0 and 1, except that $a$
may be exactly 0. If the segment from $a$ to $z$ crosses the branch
cut, the result is the analytic continuation along the segment rather
than the principal branch.
:::

::: function
void acb_hypgeom_dilog_bitburst(acb_t res, acb_t z0, const acb_t z, slong prec)

Sets `z0` to a point with a short binary expansion close to $z$ and
`res` to $\operatorname{Li}_2(z) - \operatorname{Li}_2(z_0)$, computed
with the bit-burst algorithm: $z_0$ is $z$ truncated to 16 fractional
bits, and continuation steps lead from it through the truncations of
$z$ to 32, 64, 128, ... bits to $z$. If $z$ is close to the
branch cut, it sets `z0` to $z$ and `res` to 0.
:::

::: function
void acb_hypgeom_dilog(acb_t res, const acb_t z, slong prec)

Sets `res` to the dilogarithm $\operatorname{Li}_2(z)$ with an
automatic choice of method: the series at 0 after the transformation
that minimizes $|z|$, the Bernoulli series near $e^{\pm i\pi/3}$ below
10000 bits, and expansions at the points of
`acb_hypgeom_dilog_transform` otherwise.
:::

## Summation of hypergeometric series {#series}

The functions in this section are the building blocks of the functions
above. With
$$
T(k) = \frac{\prod_{i=0}^{p-1} (a_i)_k}{\prod_{i=0}^{q-1} (b_i)_k} z^k,
$$
they compute partial sums of
$$
{}_pf_{q}(a_0,\ldots,a_{p-1}; b_0 \ldots b_{q-1}; z) =
{}_{p+1}F_{q}(a_0,\ldots,a_{p-1},1; b_0 \ldots b_{q-1}; z) =
\sum_{k=0}^{\infty} T(k).
$$
There is no implicit $k!$: to compute the usual ${}_pF_q$, call the
${}_pf_{q+1}$ functions with the extra lower parameter $b_q = 1$, or
remove an upper parameter equal to 1 if there is one. The algorithms
and the error bounds are described in
[](algorithms/special-functions/hypergeometric).

::: function
void acb_hypgeom_pfq_bound_factor(mag_t C, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, ulong n)

Sets `C` to a factor such that $\left|\sum_{k=n}^{\infty} T(k)\right|
\le C |T(n)|$. The bound is infinite if $p > q$, if $\operatorname{Re}(b_i
+ n) > 0$ fails for some $i$, or if $n$ is too small for the term ratio
to be bounded below 1, even when the series converges. See
[](algorithms/special-functions/hypergeometric#tail-bound).
:::

::: function
slong acb_hypgeom_pfq_choose_n(acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong prec)

Heuristically chooses the number of terms $n$ for summing the series
at `prec` bits of working precision.

::: details
The function estimates the terms in double precision and stops when
the current term is `prec` bits smaller than the largest term
encountered, since cancellation is assumed to be significant. It also
picks a reasonable truncation point for divergent (asymptotic) series,
where the terms first decrease and then increase, stops at the last
nonzero term of a terminating series, skips past parameters close to
nonpositive integers, and returns at most $50 + 10\,\mathit{prec}$. It
can give poor results if the parameters are extremely large or
extremely close to nonpositive integers.
:::
:::

::: function
slong acb_hypgeom_pfq_series_choose_n(const acb_poly_struct * a, slong p, const acb_poly_struct * b, slong q, const acb_poly_t z, slong len, slong prec)

Like `acb_hypgeom_pfq_choose_n` for power series parameters and
argument, using their constant terms; the result is at least `len`.
:::

::: function
void acb_hypgeom_pfq_sum_forward(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
void acb_hypgeom_pfq_sum_rs(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
void acb_hypgeom_pfq_sum_bs(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
void acb_hypgeom_pfq_sum_fme(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)
void acb_hypgeom_pfq_sum(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)

Sets $s = \sum_{k=0}^{n-1} T(k)$ and $t = T(n)$. Requires $n \ge 0$.
The outputs must not alias the inputs. No truncation error is added.

::: details
- `acb_hypgeom_pfq_sum_forward` uses the forward recurrence $T(k+1) =
  T(k) \prod (a_i + k) / \prod (b_i + k) \cdot z$, and stops early if a
  term is exactly zero.
- `acb_hypgeom_pfq_sum_rs` sums in reverse order with rectangular
  splitting over the powers of $z$ (Horner's scheme in $z^m$ with $m
  \approx \sqrt n$, at most 150), so that only $O(\sqrt n)$
  multiplications by powers of $z$ are of full cost. It computes only
  a magnitude bound for $t$ (the midpoint is zero).
- `acb_hypgeom_pfq_sum_bs` uses binary splitting on the $2 \times 2$
  matrix recurrence for $(T(k), S(k))$, which reduces the bit
  complexity when $z$ and the parameters have few bits and also limits
  the growth of the error bounds of complex balls.
- `acb_hypgeom_pfq_sum_fme` uses baby-step giant-step summation with
  fast multipoint evaluation of the matrix polynomials, for $p, q \le
  3$ (otherwise it falls back to the forward recurrence).
- `acb_hypgeom_pfq_sum` uses rectangular splitting if $n > 4$, $\mathit{prec}
  \ge 128$ and the parameters have few bits (and binary splitting if
  in addition $\mathit{prec} \ge 256$ and $z$ has fewer than
  $\mathit{prec}/100$ bits), fast multipoint evaluation from 1500
  bits when $n \ge 30 + 100000/(\mathit{prec} - 1000)$, and otherwise
  binary splitting for $n > 8$ and the forward recurrence for $n \le
  8$.
:::
:::

::: function
void acb_hypgeom_pfq_sum_bs_invz(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t w, slong n, slong prec)
void acb_hypgeom_pfq_sum_invz(acb_t s, acb_t t, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, const acb_t w, slong n, slong prec)

Like `acb_hypgeom_pfq_sum_bs` and `acb_hypgeom_pfq_sum`, but given $w
= 1/z$ (and also $z$ for the second function), which is advantageous
when $w$ has few bits, as in asymptotic expansions at an integer or
short argument: binary splitting then divides by $w$ instead of
multiplying by the long $z$.
:::

::: function
void acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, slong p, acb_srcptr b, slong q, const acb_t z, slong n, slong prec)

Sets `res` to
$$
{}_pf_{q}(z) = \sum_{k=0}^{\infty} T(k) = \sum_{k=0}^{n-1} T(k) +
\varepsilon,
$$
computed directly from the series with `acb_hypgeom_pfq_sum` and a
rigorous bound for the truncation error $\varepsilon$ from
`acb_hypgeom_pfq_bound_factor`. If $n < 0$, the number of terms is
chosen with `acb_hypgeom_pfq_choose_n`. If $T(n)$ is exactly zero (a
terminating series), no error is added. For real input the error is
added to the real part only.
:::

::: function
void acb_hypgeom_pfq_series_sum_forward(acb_poly_t s, acb_poly_t t, const acb_poly_struct * a, slong p, const acb_poly_struct * b, slong q, const acb_poly_t z, int regularized, slong n, slong len, slong prec)
void acb_hypgeom_pfq_series_sum_bs(acb_poly_t s, acb_poly_t t, const acb_poly_struct * a, slong p, const acb_poly_struct * b, slong q, const acb_poly_t z, int regularized, slong n, slong len, slong prec)
void acb_hypgeom_pfq_series_sum_rs(acb_poly_t s, acb_poly_t t, const acb_poly_struct * a, slong p, const acb_poly_struct * b, slong q, const acb_poly_t z, int regularized, slong n, slong len, slong prec)
void acb_hypgeom_pfq_series_sum(acb_poly_t s, acb_poly_t t, const acb_poly_struct * a, slong p, const acb_poly_struct * b, slong q, const acb_poly_t z, int regularized, slong n, slong len, slong prec)

Sets $s = \sum_{k=0}^{n-1} T(k)$ and $t = T(n)$ for parameters and
argument that are power series, with all arithmetic truncated to length
`len`. Requires $n \ge 0$ and $\mathit{len} \ge 1$. The outputs must
not alias the inputs.

::: details
If `regularized` is set, the terms are those of the regularized series,
divided by $\prod_i \Gamma(b_i)$, computed so that a lower parameter
whose constant term is a nonpositive integer causes no division by
zero. The `forward`, `bs` and `rs` versions use the forward recurrence,
binary splitting and rectangular splitting, and the default version
chooses automatically.
:::
:::

::: function
void acb_hypgeom_pfq_series_direct(acb_poly_t res, const acb_poly_struct * a, slong p, const acb_poly_struct * b, slong q, const acb_poly_t z, int regularized, slong n, slong len, slong prec)

Sets `res` to ${}_pf_{q}(z)$ (regularized if `regularized` is set) for
parameters and argument that are power series, truncated to length
`len` $\ge 1$, summing $n$ terms of the series and adding a rigorous
bound for the truncation error in each coefficient (see
[](algorithms/special-functions/hypergeometric#series-tail)). If $n <
0$, the number of terms is chosen automatically.
:::
