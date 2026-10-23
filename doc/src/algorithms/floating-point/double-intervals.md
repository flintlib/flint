---
title: Interval arithmetic with doubles
summary: Fast rigorous enclosures with double endpoints computed in round-to-nearest mode and widened by a relative perturbation, with a proof that the widening is sufficient, and the mean value form used to bound integrands.
complexity: A constant number of double operations per interval operation, about 3 to 10 times the cost of the plain operation.
order: 90
---

Interval arithmetic [@Moo1966] represents a real number by an interval
$[a, b]$ that contains it, and replaces each operation by one that maps
intervals to an interval containing all possible results. With
floating-point endpoints, the lower endpoint must be rounded down and
the upper endpoint up. The IEEE 1788 standard [@IEEE1788] and most
libraries do this by switching the rounding mode of the processor, which
is slow on many processors, interferes with other code, and is not
portable in C. The `double_interval` module instead computes the
endpoints in the default round-to-nearest mode and then moves them
outward by a few units in the last place. This gives slightly wider
intervals, but each operation is branch-free apart from sign tests, and
it is valid whether or not the rounded operation was exact. A similar
approach, with tighter predecessor and successor functions, is analysed
in [@RZBM2009].

FLINT has three kinds of enclosures with machine precision:

| Type | Representation | Use |
|---|---|---|
| `di_t` | interval $[a, b]$ of doubles | crude bounds, e.g. on integrands |
| `d1b_t` (`dfloat`) | ball: double midpoint and double radius | general rigorous arithmetic at 53 bits, elementary functions |
| `mag_t` | upper bound with 30-bit mantissa and unbounded exponent | radii of `arb_t` balls |

Midpoint-radius (ball) arithmetic gives tighter results for narrow
enclosures, since the radius is computed separately from the midpoint
[@Rum2010]; endpoint arithmetic is simpler and tighter for wide
intervals and monotone operations.

## Outward rounding by perturbation {#perturbation}

Let $u = 2^{-53}$. The functions `_di_below` and `_di_above` take a
double $x$ that is the round-to-nearest value of some real $t$ and
return a lower and an upper bound for $t$:
$$\mathtt{below}(x) = \operatorname{RN}\bigl(x - \operatorname{RN}(\operatorname{RN}(|x| + \tau) \cdot 4u)\bigr),
\qquad \tau = 10^{-300},$$
and symmetrically for $\mathtt{above}$, for $|x| \le 10^{300}$.

::: theorem {#thm-below}
Let $t$ be a real number such that $x = \operatorname{RN}(t)$ is finite
with $|x| \le 10^{300}$. Then $\mathtt{below}(x) \le t \le \mathtt{above}(x)$.
:::

::: proof
If $x$ is a normal double, $|t - x| \le \operatorname{ulp}(x)/2 \le u|x|$;
if $x$ is subnormal or zero, $|t - x| \le 2^{-1075}$. The perturbation
$\delta = \operatorname{RN}(\operatorname{RN}(|x| + \tau) \cdot 4u)$
satisfies $\delta \ge 4u(1 - u)^2(|x| + \tau)$, since both roundings
are of normal numbers. The final subtraction has an error of at most
$u|x - \delta| + 2^{-1075}$, so
$$\mathtt{below}(x) \le x - \delta + u(|x| + \delta) + 2^{-1075}
\le x - (3u - 12u^2)|x| - 3u\tau + 2^{-1075},$$
which is at most $x - u|x| - 2^{-1075} \le t$ because
$3u\tau \approx 2^{-1049} > 2 \cdot 2^{-1075}$. The case of $\mathtt{above}$
is symmetric.
:::

For $|x| > 10^{300}$, including infinities produced by overflow,
`_di_below` returns $10^{300}$ for positive $x$ (a valid lower bound,
since $t > 10^{300}$ when $\operatorname{RN}(t) > 10^{300}$) and the
perturbed value otherwise; NaN, which arises from $\infty - \infty$ in
endpoint arithmetic, gives $-\infty$. The widening is
$2^{-51}$ relative plus about $2^{-1047}$ absolute, four times the
maximal rounding error, so a sum of two exact operands gets width about
$2^{-49}$ relative to the result.

## Operations {#operations}

Since round-to-nearest is monotone, the endpoint formulas of exact
interval arithmetic remain valid when each endpoint expression is
computed with a single rounding and then perturbed with
[](#thm-below):

- $[a, b] + [c, d] = [\mathtt{below}(a + c), \mathtt{above}(b + d)]$,
  and similarly for subtraction;
- for products, the endpoints are selected by the signs of the
  operands (four cases with one product per endpoint), and in the
  remaining case (an operand containing zero) as the minimum and maximum
  of the four endpoint products; a product $0 \cdot \infty$ makes the
  result $[-\infty, +\infty]$;
- for quotients, the same with the divisor of fixed sign; a divisor
  containing zero gives $[-\infty, +\infty]$;
- squares are computed by sign cases, with lower endpoint exactly 0 when
  the interval contains zero;
- the logarithm of a nonnegative interval uses the bounds
  `mag_d_log_lower_bound` and `mag_d_log_upper_bound` from the `mag`
  module, which bound the error of the C library logarithm.

## Bounding functions on boxes {#mean-value}

A typical use is to bound $\sup |f|$ on a box in the complex plane, for
example to choose the path or the number of nodes of a numerical
integration. Evaluating $f$ in interval arithmetic on the whole box
overestimates because of the dependency problem; the *mean value form*
$$f(X) \subseteq f(m) + f'(X)(X - m), \qquad m = \operatorname{mid}(X),$$
gives a much tighter enclosure for narrow $X$, since the overestimation
of $f'(X)$ is multiplied by the small width of $X - m$. The integration
code of `arb_hypgeom` (`arb_hypgeom_1f1_integration`,
`arb_hypgeom_2f1_integration`, `arb_hypgeom_u_integration`) uses
`di_t` intervals to bound the logarithm of the integrand
$t^{a-1}(1 - t)^{b - a - 1}e^{zt}$ on subintervals of the edges of a
rectangle, in exactly this form: the value at the midpoint
(`di_fast_mid`) plus the derivative on the subinterval times
$[0, r]$ with $r$ from `di_fast_ubound_radius`.

## Use in FLINT {#use}

::: functions di_fast_add di_fast_mul di_fast_div di_fast_log_nonnegative di_fast_mid _di_below :::
