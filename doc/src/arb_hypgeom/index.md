---
title: arb_hypgeom.h
header: arb_hypgeom.h
description: Gamma, error, exponential-integral, Bessel, Airy, Coulomb, Legendre and other hypergeometric functions of real ball arguments, with real-specific algorithms for the gamma and error functions, Legendre polynomials and their roots, and sums of hypergeometric series with rational parameters.
category: Special functions
order: 1
---

This module evaluates hypergeometric and related special functions of
real arguments given as balls (`arb_t`) or real power series
(`arb_poly_t`), with rigorous error bounds. It provides the functions
of [`acb_hypgeom`](acb_hypgeom) for real variables. Most of them
convert their arguments to complex balls, call the complex function,
and return the real part. A smaller set is implemented directly for
real arguments because a dedicated algorithm is much faster there: the
gamma function and rising factorials, the error functions and their
inverses, the sine and cosine integrals, central binomial
coefficients, Legendre polynomials $P_n(x)$ of integer degree and their
roots, the real zeros of the Airy functions, the incomplete gamma
function at rational parameters, and finite hypergeometric sums with
rational parameters. Several functions also have variants based on
rigorous numerical integration, which work for large parameters where
series expansions fail.

[](algorithms/special-functions/overview) explains how FLINT
evaluates special functions and which methods each function uses;
[](algorithms/special-functions/gamma) and
[](algorithms/special-functions/hypergeometric) describe the two main
families of algorithms in detail. The same functions are available
through the generic interface of [`gr_special`](gr_special): for the
context of real balls (`gr_ctx_init_real_arb`), `gr_gamma`, `gr_erf`,
`gr_hypgeom_2f1`, `gr_bessel_j` and their relatives call the functions
documented here.

::: summary algorithms/special-functions/overview :::

## Conventions {#conventions}

- **Precision and accuracy.** Every function returns a ball that
  contains the exact value of the function at every point of the input
  balls. The working precision `prec` bounds the relative accuracy of
  the result, which is close to `prec` bits for exact inputs away from
  zeros and singularities of the function; the functions add guard bits
  internally to compensate for the cancellation they can predict. For
  inexact inputs, most functions evaluate at the midpoint and add a
  bound for the propagated error, and lower the working precision to
  what the accuracy of the input justifies.
- **Undefined values.** If the function is undefined or not real at
  some point of the input (a pole, a branch cut, a real argument for
  which the value is complex), or if no finite enclosure could be
  computed, the output is indeterminate (NaN midpoint, infinite radius;
  see `arb_indeterminate`). The functions that wrap a complex function
  return an indeterminate result whenever the complex result is not
  finite or has an imaginary part that is not exactly zero. An
  indeterminate or very wide result for a valid input means that the
  input is too wide or the precision too low, not that the function is
  undefined; increasing `prec` usually helps.
- **Aliasing.** Inputs and outputs of the functions taking `arb_t`
  arguments may be the same object, except where stated (the jet
  functions). The underscore power series functions require positive
  lengths, write exactly `len` coefficients and must not have the
  output overlap the input; the non-underscore versions accept any
  lengths and aliasing.
- **Optional outputs.** Functions that compute several related values
  (Fresnel integrals, Airy functions, Coulomb functions, Legendre
  polynomials with their derivative) accept `NULL` for the outputs that
  are not wanted, and then skip their computation when this saves time.
- **Regularized and normalized variants.** An `int` argument
  `regularized` selects the regularized function (divided by gamma
  functions, so that it is entire in the parameters) instead of the
  ordinary one; the exact meaning is given for each function.

The guide [](acb_hypgeom/using) illustrates these conventions with
examples.

## Rising factorials {#rising}

::: function
void arb_hypgeom_rising_ui_forward(arb_t res, const arb_t x, ulong n, slong prec)
void arb_hypgeom_rising_ui_bs(arb_t res, const arb_t x, ulong n, slong prec)
void arb_hypgeom_rising_ui_rs(arb_t res, const arb_t x, ulong n, ulong m, slong prec)
void arb_hypgeom_rising_ui_rec(arb_t res, const arb_t x, ulong n, slong prec)
void arb_hypgeom_rising_ui(arb_t res, const arb_t x, ulong n, slong prec)
void arb_hypgeom_rising(arb_t res, const arb_t x, const arb_t n, slong prec)

Sets `res` to the rising factorial $(x)_n = x (x+1) \cdots (x+n-1)$,
with $(x)_0 = 1$; `arb_hypgeom_rising` allows any real $n$ and computes
$\Gamma(x+n)/\Gamma(x)$.

::: details
The variants differ in the algorithm; see
[](algorithms/special-functions/gamma#rising).

- `arb_hypgeom_rising_ui_forward` multiplies the factors one by one at
  a working precision of `prec` plus $\log_2 n$ bits, incrementing the
  factor in place.
- `arb_hypgeom_rising_ui_bs` multiplies the factors in a balanced
  binary tree, which is efficient when $x$ has few bits compared to
  `prec` (see [](algorithms/integers/binary-splitting)).
- `arb_hypgeom_rising_ui_rs` uses rectangular splitting: it groups the
  factors in blocks of $m$, expands each block
  $(x+k)(x+k+1)\cdots(x+k+m-1)$ as a polynomial in $x$ with integer
  coefficients (one or two words each when they fit, see
  `_arb_hypgeom_rising_coeffs_1`) and evaluates it as a dot product
  with the precomputed powers $x, \ldots, x^m$. Most multiplications of
  two full-precision numbers become multiplications by small integers
  [@Smi2001; @Joh2014b]. The tuning parameter $m$ can be given as 0 to
  choose it automatically (for large $n$, $m = \min(\sqrt n, 8 + 0.2
  \max(0, \mathit{prec} - 4096)^{0.4}, 60)$).
- `arb_hypgeom_rising_ui_rec` chooses between the three by $n$,
  `prec` and the bit length of $x$: the forward recurrence for small
  $n$ (up to 20 below 512 bits, up to $\min(80, 6000/\mathit{prec})$
  otherwise), binary splitting for $n \ge 20$ if $x$ has fewer than
  $\mathit{prec}/8$ bits, and rectangular splitting otherwise. It
  never calls the gamma function, which is why the gamma function
  itself uses it for argument reduction.
- `arb_hypgeom_rising_ui` uses `arb_hypgeom_rising_ui_rec` if $n <
  \max(\mathit{prec}, 100)$, and computes $\Gamma(x+n) \cdot
  (1/\Gamma(x))$ otherwise. `arb_hypgeom_rising` does the same for an
  exact nonnegative integer $n$ in that range and uses the gamma
  function otherwise.

These are the functions behind `arb_rising_ui` and `arb_rising` of
[`arb`](arb).
:::
:::

::: function
void arb_hypgeom_rising_ui_jet_powsum(arb_ptr res, const arb_t x, ulong n, slong len, slong prec)
void arb_hypgeom_rising_ui_jet_bs(arb_ptr res, const arb_t x, ulong n, slong len, slong prec)
void arb_hypgeom_rising_ui_jet_rs(arb_ptr res, const arb_t x, ulong n, ulong m, slong len, slong prec)
void arb_hypgeom_rising_ui_jet(arb_ptr res, const arb_t x, ulong n, slong len, slong prec)

Sets the vector `res` of length `len` to the coefficients of the
polynomial $(x + t)_n = (x+t)(x+t+1)\cdots(x+t+n-1) \in \mathbb{R}[t]$,
truncated to length `len` if $\mathit{len} < n+1$ and padded with zeros
if $\mathit{len} > n+1$. The coefficient of $t^k$ is the $k$-th
derivative of $(x)_n$ with respect to $x$ divided by $k!$.

::: details
`arb_hypgeom_rising_ui_jet_powsum` computes the powers of $x$ and forms
the coefficients as integer linear combinations of them;
`arb_hypgeom_rising_ui_jet_bs` uses binary splitting of the product of
the linear polynomials $x + k + t$; `arb_hypgeom_rising_ui_jet_rs` uses
rectangular splitting with block size $m$ (0 chooses it), multiplying
the block polynomials as truncated power series.
`arb_hypgeom_rising_ui_jet` chooses: `len` = 1 is
`arb_hypgeom_rising_ui_rec`, $n \le 7$ uses the power sums, and
otherwise the choice depends on $n$, `len`, `prec` and the bit length
of $x$ (binary splitting for $\mathit{len} \ge 64$, or for $n \ge 32$
and $x$ with fewer than $\mathit{prec}/1024$ bits). The output must not
overlap $x$.
:::
:::

## Gamma function {#gamma}

::: function
void arb_hypgeom_gamma(arb_t res, const arb_t x, slong prec)
void arb_hypgeom_gamma_fmpq(arb_t res, const fmpq_t x, slong prec)
void arb_hypgeom_gamma_fmpz(arb_t res, const fmpz_t x, slong prec)

Sets `res` to the gamma function $\Gamma(x)$. The result is
indeterminate if $x$ contains a pole $0, -1, -2, \ldots$.

::: details
`arb_hypgeom_gamma` passes exact arguments that are multiples of $1/4$
(and above 1000 bits, of $2^{-\mathit{prec}/50}$) to
`arb_hypgeom_gamma_fmpq`; otherwise it uses the tabulated Taylor
series of `arb_hypgeom_gamma_taylor` if $|x|$ and `prec` are small
enough, and the Stirling series of `arb_hypgeom_gamma_stirling`
otherwise. `arb_hypgeom_gamma_fmpq` computes $\Gamma(n)$ as an exact
factorial for moderate integers, reduces arguments with denominator 2,
3, 4 or 6 to the constants $\sqrt\pi$, $\Gamma(1/3)$ and $\Gamma(1/4)$,
uses a hypergeometric series summed by binary splitting at very high
precision for other denominators, and otherwise the Taylor or Stirling
series. The selection is described in [](arb_hypgeom/gamma).
`arb_gamma`, `arb_gamma_fmpq` and `arb_gamma_fmpz` of [`arb`](arb) call
these functions.
:::
:::

::: function
void arb_hypgeom_rgamma(arb_t res, const arb_t x, slong prec)

Sets `res` to the reciprocal gamma function $1/\Gamma(x)$, which is
entire and zero at the nonpositive integers.

::: details
Uses the same methods as `arb_hypgeom_gamma`, applied to the
reciprocal, so that $x$ may contain a pole of $\Gamma$. This is the
function behind `arb_rgamma`.
:::
:::

::: function
void arb_hypgeom_lgamma(arb_t res, const arb_t x, slong prec)

Sets `res` to $\log \Gamma(x)$ for $x > 0$. The result is indeterminate
unless $x$ is positive; `acb_hypgeom_lgamma` computes the log-gamma
function of negative and complex arguments.

::: details
Computes $\log \Gamma(x)$ directly from the Stirling series, with
argument reduction by the logarithm of a rising factorial, so that
huge $x$ cause no overflow; exact rational arguments and arguments for
which the Taylor series applies go through the gamma function followed
by a logarithm. This is the function behind `arb_lgamma`.
:::
:::

::: function
void arb_hypgeom_gamma_stirling_sum_horner(arb_t res, const arb_t z, slong N, slong prec)
void arb_hypgeom_gamma_stirling_sum_improved(arb_t res, const arb_t z, slong N, slong K, slong prec)

Sets `res` to the sum $\sum_{n=1}^{N-1} B_{2n} / (2n(2n-1) z^{2n-1})$
of the Stirling series for $\log \Gamma(z)$, truncated before the term
of index $N$. No bound for the truncation error is added.

::: details
`arb_hypgeom_gamma_stirling_sum_horner` uses Horner's scheme in $1/z^2$,
lowering the precision for the smaller terms.
`arb_hypgeom_gamma_stirling_sum_improved` evaluates the first terms
with exact Bernoulli numbers by rectangular splitting and rewrites the
remaining terms, using $B_{2n} = (-1)^{n+1} 2 (2n)! \zeta(2n) /
(2\pi)^{2n}$ and $\zeta(2n) = \sum_k k^{-2n}$, as at most $K$
polynomials in $-1/(2\pi k z)^2$ with integer coefficients plus a
bounded remainder; $K = 0$ chooses the splitting parameter
automatically. See [](algorithms/special-functions/gamma#stirling-sum).
:::
:::

::: function
void arb_hypgeom_gamma_stirling(arb_t res, const arb_t x, int reciprocal, slong prec)

Sets `res` to $\Gamma(x)$, or $1/\Gamma(x)$ if `reciprocal` is
nonzero, computed with the Stirling series after argument reduction.

::: details
For $x < -5$ the reflection formula $\Gamma(x) = \pi / (\sin(\pi x)
\Gamma(1-x))$ is applied first. The argument is then shifted by a
rising factorial to $x + r \ge \beta \cdot \mathit{prec}$, where
$\beta$ is between 0.17 and 0.24 depending on the precision, and the
number of terms is chosen so that the bound for the remainder of the
Stirling series is below $2^{-\mathit{prec}}$. The working precision is
increased by the number of exponent bits of large $x$ and limited by
the accuracy of the input. Arguments with more than about
$10\,\mathit{prec} + 4096$ exponent bits give an indeterminate result.
See [](algorithms/special-functions/gamma#stirling).
:::
:::

::: function
int arb_hypgeom_gamma_taylor(arb_t res, const arb_t x, int reciprocal, slong prec)

Attempts to compute $\Gamma(x)$, or $1/\Gamma(x)$ if `reciprocal` is
nonzero, using a precomputed Taylor series. Returns 1 and sets `res` on
success; returns 0 without changing `res` if $x$ or `prec` is too large
for the method.

::: details
The method reduces $x$ to $u = x - r \in [-1/2, 1/2]$ with $r$ an
integer, evaluates the Taylor series of $1/\Gamma(1+u)$ with the
coefficients stored in `arb_hypgeom_gamma_coeffs`, and multiplies or
divides by a rising factorial. It applies for $|x| < 160$, for $r$
between about $-(40 + (\mathit{prec} - 40)/4)$ and $70 +
(\mathit{prec}-40)/8$, and up to a working precision of
`ARB_HYPGEOM_GAMMA_TAB_PREC` = 3456 bits. Wide inputs (below 40 bits
of precision, or a radius above $2^{-16}$) get an enclosure from
double-precision bounds for $\Gamma$ at the endpoints, taking the
minimum of $\Gamma$ on $(0, \infty)$ into account. Otherwise the
result has a radius of a few ulp, plus a bound for the error
propagated from the radius of $x$ computed with an upper bound for the
digamma function. See [](algorithms/special-functions/gamma#taylor).
:::
:::

## Binomial coefficients {#binomial}

::: function
void arb_hypgeom_central_bin_ui(arb_t res, ulong n, slong prec)

Sets `res` to the central binomial coefficient $\binom{2n}{n}$.

::: details
For $n \le 17$ the value comes from a table and for $n < 6\,\mathit{prec}
+ 200$ from the exact integer `fmpz_bin_uiui`. For larger $n$ the
function sums the asymptotic expansion of $\binom{2n}{n} \sqrt{\pi n} /
4^n$ in powers of $1/n^2$, whose coefficients involve Bernoulli
numbers, with the rigorous error bound of Brent [@Bre2016].
:::
:::

## Generalized hypergeometric function {#pfq}

::: function
void arb_hypgeom_pfq(arb_t res, arb_srcptr a, slong p, arb_srcptr b, slong q, const arb_t z, int regularized, slong prec)

Sets `res` to the generalized hypergeometric function
${}_pF_q(a_1, \ldots, a_p; b_1, \ldots, b_q; z)$, or the regularized
function ${}_pF_q / (\Gamma(b_1) \cdots \Gamma(b_q))$ if `regularized`
is nonzero. Wrapper of `acb_hypgeom_pfq`.
:::

## Confluent hypergeometric functions {#confluent}

::: function
void arb_hypgeom_0f1(arb_t res, const arb_t a, const arb_t z, int regularized, slong prec)

Sets `res` to the confluent hypergeometric limit function
${}_0F_1(a, z)$, or ${}_0F_1(a,z)/\Gamma(a)$ if `regularized` is
nonzero. Wrapper of `acb_hypgeom_0f1`.
:::

::: function
void arb_hypgeom_m(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)
void arb_hypgeom_1f1(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)

Sets `res` to Kummer's confluent hypergeometric function $M(a,b,z) =
{}_1F_1(a,b,z)$, or $\mathbf{M}(a,b,z) = {}_1F_1(a,b,z)/\Gamma(b)$ if
`regularized` is nonzero. The two names denote the same function, a
wrapper of `acb_hypgeom_m`.
:::

::: function
void arb_hypgeom_u(arb_t res, const arb_t a, const arb_t b, const arb_t z, slong prec)

Sets `res` to Tricomi's confluent hypergeometric function $U(a,b,z)$,
which is real for $z > 0$. Wrapper of `acb_hypgeom_u`.
:::

::: function
void arb_hypgeom_1f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)

Computes ${}_1F_1(a,b,z)$ (regularized if `regularized` is nonzero) by
rigorous numerical integration of
$$
{}_1F_1(a,b,z) = \frac{\Gamma(b)}{\Gamma(a) \Gamma(b-a)} \int_0^1
e^{zt} t^{a-1} (1-t)^{b-a-1} \, dt.
$$
The result is finite only if $a \ge 1$ and $b - a \ge 1$.

::: details
The integral is computed with `acb_calc_integrate`, which bounds the
error of Gauss–Legendre quadrature rigorously [@Joh2018a]; the bounds
for the integrand on the edges of the integration region are computed
with intervals of doubles
([interval arithmetic with doubles](algorithms/floating-point/double-intervals#mean-value)). The method
is useful when the parameters are large, where the series expansions
need many terms and lose accuracy to cancellation. `acb_hypgeom_m`
falls back to it for real arguments when its other methods give less
than half the requested accuracy.
:::
:::

::: function
void arb_hypgeom_u_integration(arb_t res, const arb_t a, const arb_t b, const arb_t z, slong prec)

Computes $U(a,b,z)$ by rigorous numerical integration of
$$
U(a,b,z) = \frac{1}{\Gamma(a)} \int_0^{\infty} e^{-zt} t^{a-1}
(1+t)^{b-a-1} \, dt.
$$
The result is finite only if $a \ge 1$ and $z > 0$. `acb_hypgeom_u`
uses it as a fallback for real arguments in the same way as
`arb_hypgeom_1f1_integration`.
:::

## Gauss hypergeometric function {#gauss}

::: function
void arb_hypgeom_2f1(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, slong prec)

Sets `res` to the Gauss hypergeometric function ${}_2F_1(a,b,c,z)$,
or $\mathbf{F}(a,b,c,z) = {}_2F_1(a,b,c,z)/\Gamma(c)$ if the bit
`ACB_HYPGEOM_2F1_REGULARIZED` (= 1) is set in `regularized`. Wrapper of
`acb_hypgeom_2f1`; the flags of that function, which declare that a
parameter difference is an exact integer, can be set in `regularized`
too. The value is real for $z < 1$; on the branch cut $z > 1$ the
result is indeterminate unless the function is a polynomial.
:::

::: function
void arb_hypgeom_2f1_integration(arb_t res, const arb_t a, const arb_t b, const arb_t c, const arb_t z, int regularized, slong prec)

Computes ${}_2F_1(a,b,c,z)$ (regularized if `regularized` is nonzero)
by rigorous numerical integration of
$$
{}_2F_1(a,b,c,z) = \frac{\Gamma(c)}{\Gamma(b) \Gamma(c-b)} \int_0^1
t^{b-1} (1-t)^{c-b-1} (1-zt)^{-a} \, dt.
$$
The result is finite only if $b \ge 1$, $c - b \ge 1$ and $z < 1$,
possibly after exchanging $a$ and $b$. `acb_hypgeom_2f1` uses it as a
fallback for real arguments.
:::

## Error functions and Fresnel integrals {#erf}

::: function
void arb_hypgeom_erf(arb_t res, const arb_t z, slong prec)

Sets `res` to the error function $\operatorname{erf}(z) =
\frac{2}{\sqrt\pi} \int_0^z e^{-t^2} dt$.

::: details
Implemented for real arguments with the series $\operatorname{erf}(z) =
\frac{2z}{\sqrt\pi} e^{-z^2} \sum_k \frac{(z^2)^k}{(3/2)_k}$, which has
positive terms and is summed by rectangular splitting; the asymptotic
expansion of $\operatorname{erfc}$ for large $|z|$; the bound
$|\operatorname{erf}(z) - \operatorname{sgn}(z)| \le
e^{-z^2}/(|z|\sqrt\pi)$ when it is below $2^{-\mathit{prec}}$; and above
30000 bits, for $|z|$ up to about $0.6 \sqrt{\mathit{prec}}$, a
bit-burst algorithm for the incomplete gamma function $\Gamma(1/2,
z^2)$. See [](arb_hypgeom/erf).
:::
:::

::: function
void arb_hypgeom_erfc(arb_t res, const arb_t z, slong prec)

Sets `res` to the complementary error function
$\operatorname{erfc}(z) = 1 - \operatorname{erf}(z)$, accurately also
for large positive $z$ where $\operatorname{erfc}(z)$ is tiny.

::: details
For $z \le 1$ this is $1 - \operatorname{erf}(z)$. For larger $z$ the
function uses the asymptotic expansion when it reaches the target
accuracy, the bit-burst algorithm above 30000 bits, and otherwise the
convergent series at a working precision increased by the
$(z^2 + \log z)/\log 2$ bits lost to cancellation. See
[](arb_hypgeom/erf).
:::
:::

::: function
void arb_hypgeom_erfi(arb_t res, const arb_t z, slong prec)

Sets `res` to the imaginary error function $\operatorname{erfi}(z) =
-i \operatorname{erf}(iz) = \frac{2}{\sqrt\pi} \int_0^z e^{t^2} dt$.
Wrapper of `acb_hypgeom_erfi`.
:::

::: function
void arb_hypgeom_erfinv(arb_t res, const arb_t x, slong prec)
void arb_hypgeom_erfcinv(arb_t res, const arb_t x, slong prec)

Sets `res` to the inverse error function $\operatorname{erf}^{-1}(x)$
for $-1 \le x \le 1$, or the inverse complementary error function
$\operatorname{erfc}^{-1}(x) = \operatorname{erf}^{-1}(1-x)$ for $0 \le
x \le 2$. The values at the endpoints are $\pm\infty$; arguments
outside the domain give an indeterminate result.

::: details
A starting approximation (in double precision [@Bla1976], or from the
first terms of the asymptotic expansion near $\pm 1$) is refined by
Newton iterations with doubling precision, and the final enclosure is
proved by an interval Newton step. `arb_hypgeom_erfcinv` works with $x$
itself rather than $1 - x$ near 0, so that tiny arguments such as
$10^{-1000}$ lose no accuracy. The error propagated from an inexact
argument is bounded with $|\frac{d}{dx}\operatorname{erf}^{-1}(x)| \le
\frac{\sqrt\pi}{2} (1-|x|)^{-1}$. See [](arb_hypgeom/erf).
:::
:::

::: function
void _arb_hypgeom_erf_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_erf_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)
void _arb_hypgeom_erfc_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_erfc_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)
void _arb_hypgeom_erfi_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_erfi_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

Sets `res` to $\operatorname{erf}$, $\operatorname{erfc}$ or
$\operatorname{erfi}$ of the power series $z$, truncated to length
`len`.

::: details
The constant term is computed with the scalar function; the other
coefficients come from integrating $\pm\frac{2}{\sqrt\pi} z'(x) e^{\mp
z(x)^2}$ as a power series.
:::
:::

::: function
void arb_hypgeom_fresnel(arb_t res1, arb_t res2, const arb_t z, int normalized, slong prec)

Sets `res1` to the Fresnel sine integral $S(z)$ and `res2` to the
Fresnel cosine integral $C(z)$; either output may be `NULL`. If
`normalized` is 0, $S(z) = \int_0^z \sin(t^2) dt$ and $C(z) = \int_0^z
\cos(t^2) dt$; if `normalized` is 1, $S(z) = \int_0^z \sin(\tfrac12
\pi t^2) dt$ and $C(z) = \int_0^z \cos(\tfrac12 \pi t^2) dt$, the
convention of Abramowitz and Stegun and of the DLMF [@NIST2012].
Wrapper of `acb_hypgeom_fresnel`.
:::

::: function
void _arb_hypgeom_fresnel_series(arb_ptr res1, arb_ptr res2, arb_srcptr z, slong zlen, int normalized, slong len, slong prec)
void arb_hypgeom_fresnel_series(arb_poly_t res1, arb_poly_t res2, const arb_poly_t z, int normalized, slong len, slong prec)

Sets `res1` and `res2` to the Fresnel sine and cosine integrals of the
power series $z$, truncated to length `len`. Either output may be
`NULL`.
:::

## Incomplete gamma and beta functions {#incomplete}

::: function
void arb_hypgeom_gamma_upper(arb_t res, const arb_t s, const arb_t z, int regularized, slong prec)

Sets `res` to the upper incomplete gamma function $\Gamma(s,z) =
\int_z^\infty t^{s-1} e^{-t} dt$ if `regularized` is 0, to $Q(s,z) =
\Gamma(s,z)/\Gamma(s)$ if `regularized` is 1, and to $z^{-s}
\Gamma(s,z) = E_{1-s}(z)$ if `regularized` is 2. Wrapper of
`acb_hypgeom_gamma_upper`.

::: details
The value is real for $z > 0$, and for $z < 0$ only in special cases
(such as $s$ a positive integer). Option 2 is mainly for internal use;
`arb_hypgeom_expint` is the interface for the exponential integral.
:::
:::

::: function
void arb_hypgeom_gamma_upper_integration(arb_t res, const arb_t s, const arb_t z, int regularized, slong prec)

Computes the upper incomplete gamma function (with the same meaning of
`regularized` as `arb_hypgeom_gamma_upper`) from $\Gamma(s,z) = z^s
e^{-z} U(1, 1+s, z)$, using `arb_hypgeom_u_integration`. The result is
finite only for $z > 0$. `acb_hypgeom_gamma_upper` uses it as a
fallback for real arguments.
:::

::: function
void _arb_hypgeom_gamma_upper_series(arb_ptr res, const arb_t s, arb_srcptr z, slong zlen, int regularized, slong n, slong prec)
void arb_hypgeom_gamma_upper_series(arb_poly_t res, const arb_t s, const arb_poly_t z, int regularized, slong n, slong prec)

Sets `res` to the upper incomplete gamma function $\Gamma(s, z)$ of the
constant $s$ and the power series $z$, truncated to length $n$.
`regularized` has the same meaning as in `arb_hypgeom_gamma_upper`.
:::

::: function
void arb_hypgeom_gamma_lower(arb_t res, const arb_t s, const arb_t z, int regularized, slong prec)

Sets `res` to the lower incomplete gamma function $\gamma(s,z) =
\int_0^z t^{s-1} e^{-t} dt = \frac{z^s}{s} {}_1F_1(s, s+1, -z)$ if
`regularized` is 0, to $P(s,z) = \gamma(s,z)/\Gamma(s)$ if
`regularized` is 1, and to $\gamma^{*}(s,z) = z^{-s} P(s,z)$, which is
entire in both variables, if `regularized` is 2. Wrapper of
`acb_hypgeom_gamma_lower`.
:::

::: function
void _arb_hypgeom_gamma_lower_series(arb_ptr res, const arb_t s, arb_srcptr z, slong zlen, int regularized, slong n, slong prec)
void arb_hypgeom_gamma_lower_series(arb_poly_t res, const arb_t s, const arb_poly_t z, int regularized, slong n, slong prec)

Sets `res` to the lower incomplete gamma function of the constant $s$
and the power series $z$, truncated to length $n$, with `regularized`
as in `arb_hypgeom_gamma_lower`.
:::

::: function
void arb_hypgeom_beta_lower(arb_t res, const arb_t a, const arb_t b, const arb_t z, int regularized, slong prec)

Sets `res` to the incomplete beta function $B(a,b;z) = \int_0^z
t^{a-1} (1-t)^{b-1} dt$, or the regularized function $I(a,b;z) =
B(a,b;z)/B(a,b;1)$ if `regularized` is nonzero. The precise definition
by analytic continuation is given for `acb_hypgeom_beta_lower`, of
which this is a wrapper.
:::

::: function
void _arb_hypgeom_beta_lower_series(arb_ptr res, const arb_t a, const arb_t b, arb_srcptr z, slong zlen, int regularized, slong n, slong prec)
void arb_hypgeom_beta_lower_series(arb_poly_t res, const arb_t a, const arb_t b, const arb_poly_t z, int regularized, slong n, slong prec)

Sets `res` to $B(a,b;z)$ (or $I(a,b;z)$ if `regularized` is nonzero)
for constants $a, b$ and the power series $z$, truncated to length
$n$.
:::

## Exponential and trigonometric integrals {#expint}

::: function
void arb_hypgeom_expint(arb_t res, const arb_t s, const arb_t z, slong prec)

Sets `res` to the generalized exponential integral $E_s(z) = z^{s-1}
\Gamma(1-s, z)$. Wrapper of `acb_hypgeom_expint`.
:::

::: function
void arb_hypgeom_ei(arb_t res, const arb_t z, slong prec)

Sets `res` to the exponential integral $\operatorname{Ei}(z) =
-\int_{-z}^\infty e^{-t}/t \, dt$ (a principal value for $z > 0$),
defined for real $z \ne 0$. Wrapper of `acb_hypgeom_ei`.
:::

::: function
void _arb_hypgeom_ei_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_ei_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

Sets `res` to the exponential integral of the power series $z$,
truncated to length `len`.
:::

::: function
void arb_hypgeom_si(arb_t res, const arb_t z, slong prec)

Sets `res` to the sine integral $\operatorname{Si}(z) = \int_0^z \sin(t)/t
\, dt$.

::: details
Implemented for real arguments: the asymptotic expansion
$\operatorname{Si}(z) = \frac{\pi}{2} - \frac{\cos z}{z} f(z) - \frac{\sin
z}{z} g(z)$ when it reaches the target accuracy, and otherwise the
series $\operatorname{Si}(z) = z\, {}_1F_2(\tfrac12; \tfrac32, \tfrac32;
-z^2/4)$ at a working precision increased by the $|z|/\log 2$ bits lost
to cancellation. Both are summed with `arb_hypgeom_sum_fmpq_arb` and
`arb_hypgeom_sum_fmpq_imag_arb`. For inexact $z$ the series is
evaluated at the midpoint and the propagated error is bounded with
$|\operatorname{Si}'(z)| \le \min(1, 1/|z|)$.
:::
:::

::: function
void _arb_hypgeom_si_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_si_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

Sets `res` to the sine integral of the power series $z$, truncated to
length `len`.
:::

::: function
void arb_hypgeom_ci(arb_t res, const arb_t z, slong prec)

Sets `res` to the cosine integral $\operatorname{Ci}(z) = -\int_z^\infty
\cos(t)/t \, dt$ for $z > 0$. The result is indeterminate unless $z$ is
positive, since the value is complex for $z < 0$.

::: details
Implemented for real arguments with the same two methods as
`arb_hypgeom_si`, using $\operatorname{Ci}(z) = \gamma + \log z -
\frac{z^2}{4}\, {}_2F_3(1, 1; 2, 2, \tfrac32; -z^2/4)$ for the
convergent series.
:::
:::

::: function
void _arb_hypgeom_ci_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_ci_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

Sets `res` to the cosine integral of the power series $z$, truncated to
length `len`.
:::

::: function
void arb_hypgeom_shi(arb_t res, const arb_t z, slong prec)

Sets `res` to the hyperbolic sine integral $\operatorname{Shi}(z) = -i
\operatorname{Si}(iz) = \int_0^z \sinh(t)/t \, dt$. Wrapper of
`acb_hypgeom_shi`.
:::

::: function
void _arb_hypgeom_shi_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_shi_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

Sets `res` to the hyperbolic sine integral of the power series $z$,
truncated to length `len`.
:::

::: function
void arb_hypgeom_chi(arb_t res, const arb_t z, slong prec)

Sets `res` to the hyperbolic cosine integral $\operatorname{Chi}(z) =
\gamma + \log z + \int_0^z (\cosh(t) - 1)/t \, dt$ for $z > 0$. The
result is indeterminate unless $z$ is positive. Wrapper of
`acb_hypgeom_chi`.
:::

::: function
void _arb_hypgeom_chi_series(arb_ptr res, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_chi_series(arb_poly_t res, const arb_poly_t z, slong len, slong prec)

Sets `res` to the hyperbolic cosine integral of the power series $z$,
truncated to length `len`.
:::

::: function
void arb_hypgeom_li(arb_t res, const arb_t z, int offset, slong prec)

Sets `res` to the logarithmic integral $\operatorname{li}(z) =
\operatorname{Ei}(\log z)$ if `offset` is zero, or the offset
logarithmic integral $\operatorname{Li}(z) = \operatorname{li}(z) -
\operatorname{li}(2)$ if `offset` is nonzero, for $z \ge 0$. The result
is indeterminate unless $z$ is nonnegative. Wrapper of
`acb_hypgeom_li`.
:::

::: function
void _arb_hypgeom_li_series(arb_ptr res, arb_srcptr z, slong zlen, int offset, slong len, slong prec)
void arb_hypgeom_li_series(arb_poly_t res, const arb_poly_t z, int offset, slong len, slong prec)

Sets `res` to the logarithmic integral (the offset version if `offset`
is nonzero) of the power series $z$, truncated to length `len`.
:::

## Bessel functions {#bessel}

::: function
void arb_hypgeom_bessel_j(arb_t res, const arb_t nu, const arb_t z, slong prec)
void arb_hypgeom_bessel_y(arb_t res, const arb_t nu, const arb_t z, slong prec)

Sets `res` to the Bessel function of the first kind $J_\nu(z)$ or of
the second kind $Y_\nu(z)$. Wrappers of `acb_hypgeom_bessel_j` and
`acb_hypgeom_bessel_y`; the result is indeterminate where the value is
not real, for example $Y_\nu(z)$ with $z < 0$.
:::

::: function
void arb_hypgeom_bessel_jy(arb_t res1, arb_t res2, const arb_t nu, const arb_t z, slong prec)

Sets `res1` to $J_\nu(z)$ and `res2` to $Y_\nu(z)$, computed together,
which is faster than two separate calls. Wrapper of
`acb_hypgeom_bessel_jy`.
:::

::: function
void arb_hypgeom_bessel_i(arb_t res, const arb_t nu, const arb_t z, slong prec)
void arb_hypgeom_bessel_i_scaled(arb_t res, const arb_t nu, const arb_t z, slong prec)

Sets `res` to the modified Bessel function of the first kind $I_\nu(z)
= z^\nu (iz)^{-\nu} J_\nu(iz)$, or to the scaled function $e^{-z}
I_\nu(z)$, which does not overflow for large $z$. Wrappers of
`acb_hypgeom_bessel_i` and `acb_hypgeom_bessel_i_scaled`.
:::

::: function
void arb_hypgeom_bessel_k(arb_t res, const arb_t nu, const arb_t z, slong prec)
void arb_hypgeom_bessel_k_scaled(arb_t res, const arb_t nu, const arb_t z, slong prec)

Sets `res` to the modified Bessel function of the second kind
$K_\nu(z)$, or to the scaled function $e^z K_\nu(z)$, which does not
underflow for large $z$. Wrappers of `acb_hypgeom_bessel_k` and
`acb_hypgeom_bessel_k_scaled`.
:::

::: function
void arb_hypgeom_bessel_i_integration(arb_t res, const arb_t nu, const arb_t z, int scaled, slong prec)
void arb_hypgeom_bessel_k_integration(arb_t res, const arb_t nu, const arb_t z, int scaled, slong prec)

Computes $I_\nu(z)$ or $K_\nu(z)$ (multiplied by $e^{-z}$ or $e^{z}$ if
`scaled` is nonzero) by rigorous numerical integration of an integral
representation. These functions are intended for large $\nu$, where the
series expansions are slow; `acb_hypgeom_bessel_i` and
`acb_hypgeom_bessel_k` use them as fallbacks for real arguments.
:::

::: function
void _arb_hypgeom_bessel_j_series(arb_ptr res, const arb_t nu, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_bessel_j_series(arb_poly_t res, const arb_t nu, const arb_poly_t z, slong len, slong prec)

Sets `res` to $J_\nu(z)$ for the constant $\nu$ and the power series
$z$, truncated to length `len`.

::: details
Both functions call the generic `_gr_poly_bessel_j_series` with a real
ball context. It computes the Taylor coefficients of $J_\nu$ at the
constant term $z_0$ of $z$ with `gr_bessel_j_jet`, from $J_\nu(z_0)$,
$J_{\nu \pm 1}(z_0)$ and a recurrence derived from the Bessel
differential equation, and composes them with $z - z_0$. The
recurrence divides by $z_0^2$, so the result is indeterminate if $z_0$
is an inexact ball containing zero; for $z_0 = 0$ exactly and an
integer $\nu$ the coefficients are computed directly.
:::
:::

## Airy functions {#airy}

::: function
void arb_hypgeom_airy(arb_t ai, arb_t ai_prime, arb_t bi, arb_t bi_prime, const arb_t z, slong prec)

Sets `ai`, `ai_prime`, `bi`, `bi_prime` to the Airy functions
$\operatorname{Ai}(z)$, $\operatorname{Ai}'(z)$, $\operatorname{Bi}(z)$,
$\operatorname{Bi}'(z)$. Any of the outputs may be `NULL`, and omitting
outputs speeds up the evaluation. Wrapper of `acb_hypgeom_airy`.
:::

::: function
void arb_hypgeom_airy_jet(arb_ptr ai, arb_ptr bi, const arb_t z, slong len, slong prec)

Sets the vectors `ai` and `bi` to the Taylor coefficients of
$\operatorname{Ai}$ and $\operatorname{Bi}$ at the point $z$, truncated
to length `len`. Either output may be `NULL`; $z$ must not alias the
outputs.

::: details
The coefficients follow from the function values and their first
derivatives at $z$ and the recurrence implied by $y'' = zy$. The
expansions of $\operatorname{Ai}'$ and $\operatorname{Bi}'$ are not
returned separately: compute one more coefficient and differentiate
with `_arb_poly_derivative`.
:::
:::

::: function
void _arb_hypgeom_airy_series(arb_ptr ai, arb_ptr ai_prime, arb_ptr bi, arb_ptr bi_prime, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_airy_series(arb_poly_t ai, arb_poly_t ai_prime, arb_poly_t bi, arb_poly_t bi_prime, const arb_poly_t z, slong len, slong prec)

Sets the outputs to the Airy functions and their derivatives evaluated
at the power series $z$, truncated to length `len`. Any of the outputs
may be `NULL`.
:::

::: function
void arb_hypgeom_airy_zero(arb_t a, arb_t a_prime, arb_t b, arb_t b_prime, const fmpz_t n, slong prec)

Sets `a`, `a_prime`, `b`, `b_prime` to the $n$-th real zero $a_n$,
$a'_n$, $b_n$, $b'_n$ of $\operatorname{Ai}$, $\operatorname{Ai}'$,
$\operatorname{Bi}$, $\operatorname{Bi}'$ respectively; any of them may
be `NULL`. The zeros are negative and indexed by increasing magnitude
starting with $n = 1$, as in the literature. An index $n \le 0$ raises
an exception.

::: details
A starting enclosure comes from a table for $n \le 10$ and otherwise
from the asymptotic expansions $a_n \sim -T(\tfrac38 \pi (4n-1))$,
$b_n \sim -T(\tfrac38\pi(4n-3))$, $a'_n \sim -U(\tfrac38\pi(4n-3))$,
$b'_n \sim -U(\tfrac38\pi(4n-1))$ [@NIST2012, §9.9], whose truncations
give proven bounds for $a_n$ and $b_n$ [@PS1991]. The enclosure is
refined by interval Newton iteration with doubling precision; for the
zeros of the derivatives, for which the asymptotic bounds are not
proved, an interval Newton step verifies that the enclosure contains
the zero.
:::
:::

## Coulomb wave functions {#coulomb}

::: function
void arb_hypgeom_coulomb(arb_t F, arb_t G, const arb_t l, const arb_t eta, const arb_t z, slong prec)

Sets `F` and `G` to the regular and irregular Coulomb wave functions
$F_\ell(\eta, z)$ and $G_\ell(\eta, z)$. Either output may be `NULL`.
An output is indeterminate if its complex value is not real. Wrapper of
`acb_hypgeom_coulomb`, which gives the definitions.
:::

::: function
void arb_hypgeom_coulomb_jet(arb_ptr F, arb_ptr G, const arb_t l, const arb_t eta, const arb_t z, slong len, slong prec)

Sets the vectors `F` and `G` to the Taylor coefficients of the Coulomb
wave functions at the point $z$, truncated to length `len`. Either
output may be `NULL`.
:::

::: function
void _arb_hypgeom_coulomb_series(arb_ptr F, arb_ptr G, const arb_t l, const arb_t eta, arb_srcptr z, slong zlen, slong len, slong prec)
void arb_hypgeom_coulomb_series(arb_poly_t F, arb_poly_t G, const arb_t l, const arb_t eta, const arb_poly_t z, slong len, slong prec)

Sets `F` and `G` to the Coulomb wave functions evaluated at the power
series $z$, truncated to length `len`. Either output may be `NULL`.
:::

## Orthogonal polynomials and functions {#orthogonal}

::: function
void arb_hypgeom_chebyshev_t(arb_t res, const arb_t nu, const arb_t z, slong prec)
void arb_hypgeom_chebyshev_u(arb_t res, const arb_t nu, const arb_t z, slong prec)
void arb_hypgeom_jacobi_p(arb_t res, const arb_t n, const arb_t a, const arb_t b, const arb_t z, slong prec)
void arb_hypgeom_gegenbauer_c(arb_t res, const arb_t n, const arb_t m, const arb_t z, slong prec)
void arb_hypgeom_laguerre_l(arb_t res, const arb_t n, const arb_t m, const arb_t z, slong prec)
void arb_hypgeom_hermite_h(arb_t res, const arb_t nu, const arb_t z, slong prec)

Sets `res` to the Chebyshev polynomial $T_\nu(z)$ or $U_\nu(z)$, the
Jacobi polynomial $P_n^{(a,b)}(z)$, the Gegenbauer polynomial
$C_n^{m}(z)$, the Laguerre polynomial $L_n^{m}(z)$ or the Hermite
polynomial $H_\nu(z)$, or their extensions to non-integer degree.
Wrappers of the functions of the same name in
[`acb_hypgeom`](acb_hypgeom), which give the definitions.
:::

::: function
void arb_hypgeom_legendre_p(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, slong prec)
void arb_hypgeom_legendre_q(arb_t res, const arb_t n, const arb_t m, const arb_t z, int type, slong prec)

Sets `res` to the associated Legendre function of the first kind
$P_n^m(z)$ or of the second kind $Q_n^m(z)$, with the branch cut
convention selected by `type` (0 or 1) as described for
`acb_hypgeom_legendre_p` and `acb_hypgeom_legendre_q`. For $m = 0$ and
a nonnegative integer $n < 2^{\bits - 1}$, `arb_hypgeom_legendre_p`
calls `arb_hypgeom_legendre_p_ui`; otherwise both functions wrap the
complex versions.
:::

::: function
void arb_hypgeom_legendre_p_ui_zero(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong K, slong prec)
void arb_hypgeom_legendre_p_ui_one(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong K, slong prec)
void arb_hypgeom_legendre_p_ui_asymp(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong K, slong prec)
void arb_hypgeom_legendre_p_ui_rec(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong prec)
void arb_hypgeom_legendre_p_ui(arb_t res, arb_t res_prime, ulong n, const arb_t x, slong prec)

Sets `res` to the Legendre polynomial $P_n(x)$ and `res_prime` to its
derivative $P'_n(x)$. Either output may be `NULL`.

::: details
The algorithm is that of Johansson and Mezzarobba [@JM2018]; see
[](algorithms/special-functions/legendre).

- `arb_hypgeom_legendre_p_ui_zero` and `arb_hypgeom_legendre_p_ui_one`
  sum the hypergeometric expansions at $x = 0$ and $x = 1$, and
  `arb_hypgeom_legendre_p_ui_asymp` the asymptotic expansion in $n$
  valid on $(-1, 1)$ with the error bound of Bogaert, Michiels and
  Fostier [@Bog2012]. The parameter $K$ is the exact number of terms;
  if the truncated expansion is not the exact polynomial, a bound for
  the remainder is added.
- `arb_hypgeom_legendre_p_ui_rec` uses the three-term recurrence in
  fixed-point arithmetic; it is meant for $x \in (-1,1)$, moderate $n$
  and modest precision.
- `arb_hypgeom_legendre_p_ui` reduces to $x \ge 0$ by symmetry, picks
  the cheapest method from estimates of the number of terms and of the
  cancellation in each series, raises the working precision to
  compensate for the cancellation, and bounds the error propagated
  from an inexact $x$ with `arb_hypgeom_legendre_p_ui_deriv_bound`.
:::
:::

::: function
void arb_hypgeom_legendre_p_ui_root(arb_t res, arb_t weight, ulong n, ulong k, slong prec)

Sets `res` to the $k$-th root $x_k$ of the Legendre polynomial $P_n$,
with the roots in decreasing order $1 > x_0 > x_1 > \cdots > x_{n-1} >
-1$, and, if `weight` is not `NULL`, sets `weight` to the corresponding
weight $w_k = 2 / ((1 - x_k^2) P'_n(x_k)^2)$ of Gauss–Legendre
quadrature on $[-1, 1]$. Requires $0 \le k < n$. See
[](arb_hypgeom/legendre-roots).
:::

## Dilogarithm {#dilog}

::: function
void arb_hypgeom_dilog(arb_t res, const arb_t z, slong prec)

Sets `res` to the dilogarithm $\operatorname{Li}_2(z) = -\int_0^z
\log(1-t)/t \, dt$ for $z \le 1$; the result is indeterminate for $z >
1$, where the value is complex. Wrapper of `acb_hypgeom_dilog`.
:::

## Hypergeometric sums {#sums}

These functions sum a finite number of terms of a hypergeometric series
with rational parameters and a real ball argument. They are the
building blocks of several real functions (sine and cosine integrals,
error functions, incomplete gamma functions) and are useful on their
own for constants and functions given by such series. No truncation
error is added: the caller chooses $N$ and bounds the tail. See
[](arb_hypgeom/sums) and
[](algorithms/special-functions/hypergeometric-sums).

::: function
void arb_hypgeom_sum_fmpq_arb_forward(arb_t res, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
void arb_hypgeom_sum_fmpq_arb_rs(arb_t res, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
void arb_hypgeom_sum_fmpq_arb_bs(arb_t res, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
void arb_hypgeom_sum_fmpq_arb(arb_t res, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)

Sets `res` to the finite sum $\sum_{n=0}^{N-1} \frac{(a_1)_n \cdots
(a_{\mathit{alen}})_n}{(b_1)_n \cdots (b_{\mathit{blen}})_n} w^n$ with
rational parameters $a_i$ and $b_j$, where $w = z$, or $w = 1/z$ if
`reciprocal` is nonzero. There is no implicit $n!$: include $b_j = 1$
for the usual ${}_pF_q$ normalisation.

::: details
`arb_hypgeom_sum_fmpq_arb_forward` uses the forward recurrence for the
terms, delaying the divisions; `arb_hypgeom_sum_fmpq_arb_rs` uses
rectangular splitting, so that most operations are multiplications of
the full-precision powers of $w$ by integers; `arb_hypgeom_sum_fmpq_arb_bs`
uses binary splitting, which is fastest when $z$ has few bits.
`arb_hypgeom_sum_fmpq_arb` uses the forward recurrence for few terms
($N \le 2$, $N \le 8$ up to 1024 bits, $N \le 4$ up to 4096 bits),
binary splitting from 8192 bits if $z$ has at most $\mathit{prec}/1000$
bits, and rectangular splitting otherwise.
:::
:::

::: function
void arb_hypgeom_sum_fmpq_imag_arb_forward(arb_t res1, arb_t res2, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
void arb_hypgeom_sum_fmpq_imag_arb_rs(arb_t res1, arb_t res2, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
void arb_hypgeom_sum_fmpq_imag_arb_bs(arb_t res1, arb_t res2, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)
void arb_hypgeom_sum_fmpq_imag_arb(arb_t res1, arb_t res2, const fmpq * a, slong alen, const fmpq * b, slong blen, const arb_t z, int reciprocal, slong N, slong prec)

Sets `res1` and `res2` to the real and imaginary parts of the finite
sum $\sum_{n=0}^{N-1} \frac{(a_1)_n \cdots}{(b_1)_n \cdots} (iw)^n$ with
$w = z$, or $w = 1/z$ if `reciprocal` is nonzero, using only real
arithmetic.

::: details
The real and imaginary parts are the sums over the even and the odd
indices, with alternating signs. The algorithm variants are as for
`arb_hypgeom_sum_fmpq_arb`; the default uses the forward recurrence for
$N \le 10$ (or $N \le 16$ up to 1024 bits). The sine and cosine
integrals use these sums for their asymptotic expansions.
:::
:::

## Internal functions {#internal}

::: function
void _arb_hypgeom_rising_coeffs_1(ulong * c, ulong k, slong n)
void _arb_hypgeom_rising_coeffs_2(ulong * c, ulong k, slong n)
void _arb_hypgeom_rising_coeffs_fmpz(fmpz * c, ulong k, slong n)

Sets $c_0, \ldots, c_{n-1}$ to the coefficients of the monic polynomial
$(X+k)_n = X^n + c_{n-1} X^{n-1} + \cdots + c_0$ (the leading 1 is not
stored). Requires $n \ge 2$ and raises an exception otherwise.

::: details
The `1` version stores each coefficient in one word and the `2`
version in two words (low word first, so $c$ must have room for $2n$
words); neither checks for overflow, and the caller must ensure that
the coefficients fit, for example by checking that $n
\lceil\log_2(k+n)\rceil$ is at most one or two word sizes. The `fmpz`
version has no restriction; it uses the quadratic recurrence $(X+k)_{i+1} =
(X+k)_i (X + k + i)$, which is not asymptotically fast. These are the
block polynomials of rectangular splitting in
`arb_hypgeom_rising_ui_rs` and `arb_hypgeom_rising_ui_jet_rs`.
:::
:::

::: function
void _arb_hypgeom_gamma_stirling_term_bounds(slong * bound, const mag_t zinv, slong N)

For $1 \le n < N$, sets `bound[n]` to an exponent $e$ such that the
$n$-th term $|B_{2n}/(2n(2n-1) z^{2n-1})|$ of the Stirling series is
at most $2^e$, given an upper bound `zinv` for $|1/z|$. Does not check
for overflow or underflow of the exponents.
:::

::: macro
ARB_HYPGEOM_GAMMA_TAB_NUM
ARB_HYPGEOM_GAMMA_TAB_PREC

The number of entries (536) of the table of Taylor coefficients of
$1/\Gamma(1+x)$ and their precision in bits (3456).
:::

::: type
arb_hypgeom_gamma_coeff_t

Describes one tabulated coefficient: its exponent (`exp`), its
position in the limb table (`tab_pos`), its number of limbs
(`nlimbs`) and its sign (`negative`).
:::

::: variable
arb_hypgeom_gamma_coeffs

The table `arb_hypgeom_gamma_coeff_t
arb_hypgeom_gamma_coeffs[ARB_HYPGEOM_GAMMA_TAB_NUM]` describing the
Taylor coefficients $c_k$ of $1/\Gamma(1+x) = \sum_k c_k x^k$, used by
`arb_hypgeom_gamma_taylor` and `acb_hypgeom_gamma_taylor`.
:::

::: function
int _arb_hypgeom_gamma_coeff_shallow(arf_t c, mag_t err, slong i, slong prec)

Sets `c` to a shallow copy of the tabulated coefficient $c_i$ of
$1/\Gamma(1+x)$, $1 \le i <$ `ARB_HYPGEOM_GAMMA_TAB_NUM`, truncated to
`prec` bits (rounded up to whole limbs), and `err` (if not `NULL`) to a
bound for the truncation error. Returns 0, without setting anything, if
`prec` exceeds the precision of the table and 1 otherwise. `c` points
into static data and must not be modified or cleared; $i$ is not
checked.
:::

::: function
void _arb_hypgeom_gamma_lower_sum_rs_1(arb_t res, ulong p, ulong q, const arb_t z, slong N, slong prec)

Sets `res` to $\sum_{k=0}^{N-1} z^k / (a)_k$ with $a = p/q$, computed
by rectangular splitting. Requires that $p + qN$ fits in a word. Used
for the error function.
:::

::: function
void _arb_hypgeom_gamma_upper_sum_rs_1(arb_t res, ulong p, ulong q, const arb_t z, slong N, slong prec)

Sets `res` to $\sum_{k=0}^{N-1} (a)_k / z^k$ with $a = p/q$, computed
by rectangular splitting. Requires that $p + qN$ fits in a word. Used
for the asymptotic expansion of the error function.
:::

::: function
slong _arb_hypgeom_gamma_upper_fmpq_inf_choose_N(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol)

Returns a number of terms $N$ of the asymptotic series of $\Gamma(a,
z)$ at infinity that gives an absolute error at most `abs_tol`, and
sets `err` to a bound for the truncation error. Returns $-1$ if the
tolerance cannot be reached. Assumes $z > 0$.
:::

::: function
void _arb_hypgeom_gamma_upper_fmpq_inf_bsplit(arb_t res, const fmpq_t a, const arb_t z, slong N, slong prec)

Sets `res` to the asymptotic series of $\Gamma(a,z)$ at infinity,
including the prefactor $z^{a-1} e^{-z}$, truncated before the term of
index $N$ and computed by binary splitting. The truncation error is
not included.
:::

::: function
slong _arb_hypgeom_gamma_lower_fmpq_0_choose_N(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol)

Returns a number of terms $N$ of the Taylor series of $\gamma(a,z)$ at
zero that gives an absolute error at most `abs_tol`, and sets `err` to
a bound for the truncation error. Assumes $z > 0$.
:::

::: function
void _arb_hypgeom_gamma_lower_fmpq_0_bsplit(arb_t res, const fmpq_t a, const arb_t z, slong N, slong prec)

Sets `res` to the series $\gamma(a,z) = \frac{z^a e^{-z}}{a} \sum_n
\frac{z^n}{(a+1)_n}$ truncated before the term of index $N$ and
computed by binary splitting. The truncation error is not included.
:::

::: function
slong _arb_hypgeom_gamma_upper_singular_si_choose_N(mag_t err, slong n, const arb_t z, const mag_t abs_tol)

Returns a number of terms $N$ of the series at zero of $\Gamma(-n,z)$,
for an integer $n \ge 0$, that gives an absolute error at most
`abs_tol`, and sets `err` to a bound for the truncation error.
:::

::: function
void _arb_hypgeom_gamma_upper_singular_si_bsplit(arb_t res, slong n, const arb_t z, slong N, slong prec)

Sets `res` to the series at zero of $\Gamma(-n,z)$, for an integer $n
\ge 0$, truncated before the term of index $N$ and computed by binary
splitting. The truncation error is not included.
:::

The following functions compute jets $\Gamma(a + \varepsilon, z)$ (truncated
power series in $\varepsilon$ of length `len`) for rational $a$, and are
used by `acb_dirichlet_l_fmpq_afe_jet`. The `choose_N` functions set `err`
to a bound for the truncation error which is valid simultaneously for every
coefficient of the jet; this bound is computed from majorant series
for the differentiated terms of the corresponding scalar series.
With `len = 1`, they reduce to the corresponding scalar functions.
Where a parameter `logz` (or `logz0`) appears, it may be passed an
enclosure of $\log(z)$ to avoid recomputing it, or `NULL`.

::: function
void _arb_hypgeom_gamma_fmpq_jet(arb_ptr res, const fmpq_t a, slong len, slong prec)

Sets `res` to the jet of $\Gamma(a + \varepsilon)$ for rational $a$ which is not a
nonpositive integer, computed as $\gamma(a + \varepsilon, Z) + \Gamma(a + \varepsilon, Z)$ with
an integer $Z \approx 0.55 \log(2) p$ at precision $p$, as in
`arb_hypgeom_gamma_fmpq`: the lower series at $Z$ and the asymptotic
series at $Z$ (at reduced precision, since $\Gamma(a+\varepsilon, Z)$ is
small) are both evaluated by binary splitting. If the asymptotic
series does not converge to the required accuracy (at low precision or
for large $a$), $Z$ is increased until it does, or until the second term is
negligible. With `len` = 1, this calls `arb_gamma_fmpq`.
:::

::: function
slong _arb_hypgeom_gamma_upper_fmpq_inf_choose_N_jet(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol, slong len)
void _arb_hypgeom_gamma_upper_fmpq_inf_bsplit_jet(arb_ptr res, const fmpq_t a, const arb_t z, const arb_t logz, slong N, slong len, slong prec)

Jet versions of `_arb_hypgeom_gamma_upper_fmpq_inf_choose_N` and
`_arb_hypgeom_gamma_upper_fmpq_inf_bsplit`. The truncation error
$(a-1)\cdots(a-N) \Gamma(a-N,z)$ is bounded using majorants for the
polynomial factor and for $\int_z^{\infty} t^{a-N-1} e^{-t} \log^j(t) / j! \, dt$.
The `choose_N` function requires $z \ge 1$ (otherwise it returns $-1$).
:::

::: function
slong _arb_hypgeom_gamma_lower_fmpq_0_choose_N_jet(mag_t err, const fmpq_t a, const arb_t z, const mag_t abs_tol, slong len)
void _arb_hypgeom_gamma_lower_fmpq_0_bsplit_jet(arb_ptr res, const fmpq_t a, const arb_t z, const arb_t logz, slong N, slong len, slong prec)

Jet versions of `_arb_hypgeom_gamma_lower_fmpq_0_choose_N` and
`_arb_hypgeom_gamma_lower_fmpq_0_bsplit`.
:::

::: function
slong _arb_hypgeom_gamma_upper_singular_si_choose_N_jet(mag_t err, slong n, const arb_t z, const mag_t abs_tol, slong len)
void _arb_hypgeom_gamma_upper_singular_si_bsplit_jet(arb_ptr res, slong n, const arb_t z, const arb_t logz, slong N, arb_srcptr gamma1, slong len, slong prec)

Jet versions of `_arb_hypgeom_gamma_upper_singular_si_choose_N` and
`_arb_hypgeom_gamma_upper_singular_si_bsplit`, computing the jet of
$\Gamma(-n + \varepsilon, z)$ for an integer $n \ge 0$. This uses

$$\Gamma(-n + \varepsilon, z) = \frac{(-1)^n}{n!} \frac{1}{\varepsilon}\left[\frac{\Gamma(1+\varepsilon)}{\prod_{i=1}^n (1 - \varepsilon/i)} - z^{\varepsilon}\right] - z^{-n+\varepsilon} \sum_{k \ne n} \frac{(-z)^k}{k! (k - n + \varepsilon)}.$$

If `gamma1` is not `NULL`, it should contain the jet of
$\Gamma(1 + \varepsilon)$ of length `len + 1`; otherwise this is
computed internally.
:::

::: function
void _arb_gamma_upper_fmpq_step_bsplit_jet(arb_ptr Gz1, const fmpq_t a, const arb_t z0, const arb_t logz0, const arb_t z1, arb_srcptr Gz0, const arb_t expmz0, const mag_t abs_tol, slong len, slong prec)

Jet version of `_arb_gamma_upper_fmpq_step_bsplit`: given the jet `Gz0`
of $\Gamma(a + \varepsilon, z_0)$ and `expmz0` containing $e^{-z_0}$, sets
`Gz1` to the jet of $\Gamma(a + \varepsilon, z_1)$. The matrix binary splitting
is done with power series entries. The truncation error is bounded
using Cauchy estimates on a circle of radius $R$ around $z_0$, with an
extra factor $\max_{j < \text{len}} \Lambda^j / j!$ where $\Lambda$ bounds
$|\log(t)|$ on the circle; the radius $R$ is chosen with this factor
taken into account.
:::

::: function
void _arb_gamma_upper_fmpq_step_bsplit(arb_t Gz1, const fmpq_t a, const arb_t z0, const arb_t z1, const arb_t Gz0, const arb_t expmz0, const mag_t abs_tol, slong prec)

Given `Gz0` and `expmz0` containing $\Gamma(a, z_0)$ and $e^{-z_0}$,
sets `Gz1` to $\Gamma(a, z_1)$, computed from the Taylor series at
$z_0$ by binary splitting with an absolute truncation error at most
`abs_tol` (included in the output). Assumes $z_0, z_1 > 0$. This is one
step of the bit-burst algorithm.
:::

::: function
int arb_hypgeom_erf_bb(arb_t res, const arb_t z, int complementary, slong prec)

Attempts to compute $\operatorname{erf}(z)$, or
$\operatorname{erfc}(z)$ if `complementary` is nonzero, with the
bit-burst algorithm applied to $\Gamma(1/2, z^2)$. Returns 1 on
success. Returns 0 without changing `res` if $z$ is not finite, if
$|z| > \mathit{prec}$ or if $|z| < 2^{-\mathit{prec}/16}$.
:::

::: function
void _arb_hypgeom_si_asymp(arb_t res, const arb_t z, slong N, slong prec)
void _arb_hypgeom_si_1f2(arb_t res, const arb_t z, slong N, slong wp, slong prec)

Compute $\operatorname{Si}(z)$ with $N$ terms of the asymptotic
expansion, or of the ${}_1F_2$ series summed at working precision `wp`
and rounded to `prec`, including a bound for the truncation error.
These are the two methods of `arb_hypgeom_si`.
:::

::: function
void _arb_hypgeom_ci_asymp(arb_t res, const arb_t z, slong N, slong prec)
void _arb_hypgeom_ci_2f3(arb_t res, const arb_t z, slong N, slong wp, slong prec)

Compute $\operatorname{Ci}(z)$ with $N$ terms of the asymptotic
expansion or of the ${}_2F_3$ series, including a bound for the
truncation error. These are the two methods of `arb_hypgeom_ci`.
:::

::: function
void arb_hypgeom_legendre_p_ui_deriv_bound(mag_t dp, mag_t dp2, ulong n, const arb_t x, const arb_t x2sub1)

Sets `dp` to an upper bound for $|P'_n(x)|$ and `dp2` to an upper bound
for $|P''_n(x)|$, for $x$ representing real numbers with $|x| \le 1$,
given the precomputed `x2sub1` $= 1 - x^2$ (or $x^2 - 1$). Used to
bound the propagated error in `arb_hypgeom_legendre_p_ui`.
:::
