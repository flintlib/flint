---
title: The gamma function of a real ball
functions: arb_hypgeom_gamma arb_hypgeom_rgamma arb_hypgeom_lgamma
---

`arb_hypgeom_gamma`, `arb_hypgeom_rgamma` and `arb_hypgeom_lgamma`
compute $\Gamma(x)$, $1/\Gamma(x)$ and $\log\Gamma(x)$ for a real ball
$x$; `arb_hypgeom_gamma_fmpq` and `arb_hypgeom_gamma_fmpz` take an exact
rational or integer argument. They are the implementations behind
`arb_gamma`, `arb_rgamma`, `arb_lgamma`, `arb_gamma_fmpq` and
`arb_gamma_fmpz`. Each combines three methods: exact or
hypergeometric formulas for rational arguments with small denominators,
a tabulated Taylor series of $1/\Gamma$ at low and medium precision,
and the Stirling series with argument reduction for large arguments
and high precision. The algorithms are described with their error
bounds in [](algorithms/special-functions/gamma) and in
[@Joh2023g].

::: summary algorithms/special-functions/gamma :::

## Values and domains {#domains}

- $\Gamma$ has poles at $0, -1, -2, \ldots$; if $x$ contains one,
  `arb_hypgeom_gamma` returns an indeterminate result.
- $1/\Gamma$ is entire; `arb_hypgeom_rgamma` returns exactly zero at
  the poles of $\Gamma$ and a finite enclosure for balls containing
  them (if necessary by computing $x \cdot (1/\Gamma(x+1))$).
- `arb_hypgeom_lgamma` is defined only for $x > 0$, where $\Gamma(x) >
  0$; for other $x$ the result is indeterminate. The log-gamma function
  of negative and complex arguments, which is the principal branch
  $\log \Gamma$ (not $\log |\Gamma|$), is `acb_hypgeom_lgamma`.
- Ball exponents are arbitrary-precision integers, so the functions do
  not overflow: $\Gamma(10^{30})$ is a ball with a binary exponent of
  about $10^{32}$, and `arb_hypgeom_lgamma` computes $\log\Gamma(x)$
  without forming $\Gamma(x)$ unless $x$ is an exact small fraction.
  Arguments larger than about $2^{10\,\mathit{prec}}$ give an
  indeterminate result.

## Choice of method {#selection}

`arb_hypgeom_gamma` and `arb_hypgeom_rgamma` proceed as follows.

1. **Exact small fractions.** If $x$ is exact, $|x| < 2^{\mathit{prec}}$
   and $4x$ is an integer (above 1000 bits also if $x$ has at most
   $\mathit{prec}/50$ fractional bits), $x$ is converted to a rational
   number and passed to `arb_hypgeom_gamma_fmpq`. Infinite and NaN
   midpoints, and the poles, are handled at this point too.
2. **Taylor series.** `arb_hypgeom_gamma_taylor` is tried next. It
   applies to $|x| < 160$, up to a working precision of 3456 bits, as
   long as the shift $r$ of the argument reduction (the nearest integer
   to $x$) satisfies $-(40 + (\mathit{prec}-40)/4) \le r \le 70 +
   (\mathit{prec}-40)/8$. In this range it is the fastest method by
   far: the work is one polynomial evaluation with tabulated
   coefficients and at most a short rising factorial.
3. **Stirling series.** Otherwise `arb_hypgeom_gamma_stirling` applies
   the reflection formula for $x < -5$, shifts the argument up to $x +
   r \approx \beta\,\mathit{prec}$ with $\beta \approx 0.2$, and sums
   the Stirling series.

For inexact $x$, the Taylor method evaluates at the midpoint and adds
a bound for the propagated error, while the Stirling method lowers the
working precision to the accuracy of the input plus 20 bits; very wide
balls are bounded directly using the monotonicity of $\Gamma$ on each
side of its minimum at $x \approx 1.4616$.

`arb_hypgeom_gamma_fmpq` uses the following methods for $x = p/q$:

- $\Gamma(n) = (n-1)!$ computed exactly for a positive integer $n$ when
  $n < 1200$ or $n \log n$ is small compared to the precision;
- for $q \in \{2, 3, 4, 6\}$ and $|p|$ not too large, reduction to
  $\Gamma(a)$ with $0 < a \le 1$ by a rising factorial, and
  $\Gamma(1/2) = \sqrt\pi$, $\Gamma(1/3)$, $\Gamma(1/4)$ and the
  reflection and multiplication formulas for the others; the two
  constants are computed once per precision and cached;
- for other denominators above about $7000 + 300 \log_2 q$ bits, a
  hypergeometric series for $\Gamma(a)$ summed by binary splitting,
  which costs $O(M(n) \log^2 n)$ at $n$ bits;
- otherwise the Taylor or Stirling series at the rational argument
  rounded to a ball.

`arb_hypgeom_lgamma` uses the same exact and Taylor methods followed by
a logarithm, and otherwise a Stirling series for $\log\Gamma$, so that
no exponential of a huge number is formed.

## Example {#example}

The program computes $\Gamma(1/3)$, $\Gamma$ of a ball, the values at a
pole and $\log\Gamma(10^{30})$:

```c example name=arb_hypgeom_gamma_values
#include <flint/fmpq.h>
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>

int main(void)
{
    arb_t x, y;
    fmpq_t q;

    arb_init(x);
    arb_init(y);
    fmpq_init(q);

    /* Gamma(1/3) from the cached constant */
    fmpq_set_si(q, 1, 3);
    arb_hypgeom_gamma_fmpq(y, q, 128);
    flint_printf("Gamma(1/3)      = "); arb_printn(y, 35, 0); flint_printf("\n");

    /* a ball argument: the result contains Gamma(t) for all t in x */
    arb_set_str(x, "3.25 +/- 1e-10", 128);
    arb_hypgeom_gamma(y, x, 128);
    flint_printf("Gamma(3.25+-e)  = "); arb_printn(y, 20, 0); flint_printf("\n");

    /* the reciprocal gamma function is zero at the poles */
    arb_set_si(x, -3);
    arb_hypgeom_rgamma(y, x, 128);
    flint_printf("1/Gamma(-3)     = "); arb_printn(y, 20, 0); flint_printf("\n");
    arb_hypgeom_gamma(y, x, 128);
    flint_printf("Gamma(-3)       = "); arb_printn(y, 20, 0); flint_printf("\n");

    /* log Gamma(10^30) does not overflow */
    arb_set_str(x, "1e30", 128);
    arb_hypgeom_lgamma(y, x, 128);
    flint_printf("lgamma(10^30)   = "); arb_printn(y, 38, 0); flint_printf("\n");

    arb_clear(x);
    arb_clear(y);
    fmpq_clear(q);
    flint_cleanup();
    return 0;
}
```
```output
Gamma(1/3)      = [2.6789385347077476336556929409746776 +/- 4.42e-35]
Gamma(3.25+-e)  = [2.549256967 +/- 5.45e-10]
1/Gamma(-3)     = 0
Gamma(-3)       = nan
lgamma(10^30)   = [68077552789821370520539743640497.306390 +/- 3.10e-7]
```

For the ball $3.25 \pm 10^{-10}$, the radius of the result is of the
order of $|\Gamma'(3.25)| \cdot 10^{-10}$: the error propagated from
the input dominates, and the midpoint is printed only to the digits
that the radius justifies. The value $\log\Gamma(10^{30})
\approx 6.8 \cdot 10^{31}$ has 126 correct bits, the accuracy
requested.
