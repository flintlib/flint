---
title: Elementary functions
summary: What the exponential, logarithmic, trigonometric and hyperbolic functions of arb compute, how accurate they are, what happens outside their domains and for wide or huge arguments, and which variants to use near cancellation.
order: 2
---

The elementary functions of `arb` (`arb_exp`, `arb_log`, `arb_sin`,
`arb_atan`, `arb_sinh` and so on) take a ball and a precision and return
an enclosure of the range of the function on the ball. For exact or
accurate input, the result is accurate to within a few units in the last
place at the given precision, so they behave like correctly rounded
floating-point functions with an attached error bound. The algorithms
are described in [](algorithms/special-functions/elementary).

## Accuracy {#accuracy}

For an input ball $[m \pm r]$, the midpoint of the result is $f(m)$
computed to `prec` bits of relative accuracy, and the radius is the
error of that computation plus a bound for the propagated error
$\sup_{|a| \le r} |f(m + a) - f(m)|$. The propagated error bound is
nearly optimal for narrow balls (it uses the derivative of $f$ near $m$);
for balls with less than about 20 bits of relative accuracy, the
functions compute the range from the endpoints instead, so that even very
wide balls get good enclosures: $\sin([0 \pm 10])$ is $[\pm 1]$, not
$[\pm 10]$.

Relative accuracy is preserved where it makes sense: `arb_exp` of a tiny
argument, `arb_log` near 1, `arb_sin` near a multiple of $\pi$ and
`arb_atan` near 0 all give results with full relative accuracy, as long
as the input is exact or accurate. But the functions cannot recover
accuracy that the input has lost. An expression such as $e^x - 1$ for
small $x$ or $\log(1 + x)$ computed from $1 + x$ suffers from
cancellation in the subtraction or in forming $1 + x$; use the dedicated
functions `arb_expm1`, `arb_log1p`, `arb_sqrt1pm1` and `arb_log_hypot`
instead. In the same way, `arb_sin_pi(x)` computes $\sin(\pi x)$
accurately near the zeros and exactly at integers, which multiplying an
enclosure of $\pi$ by $x$ cannot do.

## Domains {#domains}

Where a function is not real-valued on (part of) the input ball, the
result is indeterminate, $[\operatorname{NaN} \pm \infty]$:

- `arb_log`, `arb_sqrt` and `arb_rsqrt` when the ball contains
  nonpositive (for the square root: negative) numbers;
- `arb_asin` and `arb_acos` outside $[-1, 1]$, `arb_atanh` outside
  $(-1, 1)$ and `arb_acosh` below 1;
- `arb_pow` with a base containing negative numbers and a non-integer
  exponent.

This also happens when only a tiny part of the ball lies outside the
domain, for example for $\sqrt{x}$ with $x = [0 \pm 10^{-20}]$, the
result of a computation that should give 0. When the negative part is
known to be spurious, use `arb_sqrtpos`, which computes the square root
of the nonnegative part of the ball. The complex functions of `acb` are
defined everywhere except at singularities, with branch cuts on the
standard rays.

## Huge arguments {#huge}

The cost of evaluating a function is bounded by a polynomial in the
precision, whatever the input (see
[](guide/types/ball-arithmetic#polynomial-time)). For arguments so
large that an accurate result would be too expensive, crude bounds are
returned: $\sin$ and $\cos$ give $[\pm 1]$ above $2^{\max(65536,
4\,\mathit{prec})}$, and $\exp$ gives $[\pm \infty]$ or a tiny ball near
0 above $2^{\max(128, 2\,\mathit{prec})}$. Below these limits, large
arguments are reduced exactly: $\sin(10^{30})$ is computed with full
accuracy at 64 bits, because the reduction modulo $\pi/4$ uses as many
bits of $\pi$ as needed.

An inexact large argument, however, has a large absolute error. $10^{30}$
rounded to 64 bits has an absolute error of about $2^{36}$, and the sine of
that ball is honestly $[\pm 1]$.

```c example name=arb_elementary
#include <flint/flint.h>
#include <flint/arb.h>

#define SHOW(label, x) do { flint_printf("%-22s", label); \
    arb_printn(x, 20, 0); flint_printf("\n"); } while (0)

int main(void)
{
    arb_t x, y, z;
    slong prec = 64;
    arb_init(x); arb_init(y); arb_init(z);

    /* accurate results for exact input, including large arguments */
    arb_set_str(x, "1e30", 128);      /* exact: 10^30 has 70 bits */
    arb_sin(y, x, prec);              SHOW("sin(10^30)", y);
    arb_set_str(x, "1e30", prec);     /* rounded to 64 bits */
    arb_sin(y, x, prec);              SHOW("sin([10^30 +/- r])", y);
    arb_set_ui(x, 1000);
    arb_exp(y, x, prec);              SHOW("exp(1000)", y);
    arb_set_str(x, "1e-30", prec);
    arb_expm1(y, x, prec);            SHOW("expm1(10^-30)", y);
    arb_exp(y, x, prec);
    arb_sub_ui(y, y, 1, prec);        SHOW("exp(10^-30) - 1", y);

    /* propagated error and wide input */
    arb_set_str(x, "[2 +/- 1e-10]", prec);
    arb_log(y, x, prec);              SHOW("log([2 +/- 1e-10])", y);
    arb_set_str(x, "[0 +/- 10]", prec);
    arb_sin(y, x, prec);              SHOW("sin([0 +/- 10])", y);
    arb_set_str(x, "[0.5 +/- 1]", prec);
    arb_log(y, x, prec);              SHOW("log([0.5 +/- 1])", y);
    arb_sqrt(y, x, prec);             SHOW("sqrt([0.5 +/- 1])", y);
    arb_sqrtpos(y, x, prec);          SHOW("sqrtpos([0.5 +/- 1])", y);

    /* functions of pi x are exact at simple rational points */
    arb_set_str(x, "1000000.5", prec);
    arb_cos_pi(y, x, prec);           SHOW("cos(pi * 1000000.5)", y);
    arb_const_pi(z, prec);
    arb_mul(z, z, x, prec);
    arb_cos(y, z, prec);              SHOW("cos(1000000.5 pi)", y);

    arb_clear(x); arb_clear(y); arb_clear(z);
    flint_cleanup();
    return 0;
}
```
```output
sin(10^30)            [-0.0901169019121380580 +/- 3.16e-20]
sin([10^30 +/- r])    [+/- 1.00]
exp(1000)             [1.970071114017046994e+434 +/- 2.66e+415]
expm1(10^-30)         [1.0000000000000000000e-30 +/- 6.36e-50]
exp(10^-30) - 1       [+/- 1.01e-30]
log([2 +/- 1e-10])    [0.6931471806 +/- 9.01e-11]
sin([0 +/- 10])       [+/- 1.00]
log([0.5 +/- 1])      nan
sqrt([0.5 +/- 1])     nan
sqrtpos([0.5 +/- 1])  [+/- 1.23]
cos(pi * 1000000.5)   0
cos(1000000.5 pi)     [+/- 3.64e-13]
```

## The functions {#list}

| Function | Notes |
|---|---|
| `arb_exp`, `arb_expm1`, `arb_exp_invexp` | `arb_expm1` for $e^x - 1$ near 0 |
| `arb_log`, `arb_log1p`, `arb_log_ui`, `arb_log_base_ui`, `arb_log_hypot` | `arb_log_base_ui` is exact for exact powers |
| `arb_pow`, `arb_pow_ui`, `arb_pow_fmpq`, `arb_root_ui`, `arb_sqrt`, `arb_rsqrt` | integer and half-integer exponents by repeated multiplication |
| `arb_sin`, `arb_cos`, `arb_sin_cos`, `arb_tan`, `arb_cot`, `arb_sec`, `arb_csc`, `arb_sinc` | `arb_sin_cos` computes both for the price of one |
| `arb_sin_pi`, `arb_cos_pi`, `arb_tan_pi`, `arb_cot_pi`, `arb_csc_pi`, `arb_sinc_pi`, `arb_sin_cos_pi_fmpq` | exact at integers and half-integers |
| `arb_atan`, `arb_atan2`, `arb_asin`, `arb_acos` | `arb_atan2` is the argument of $a + bi$ |
| `arb_sinh`, `arb_cosh`, `arb_sinh_cosh`, `arb_tanh`, `arb_coth`, `arb_sech`, `arb_csch` | accurate near 0 |
| `arb_asinh`, `arb_acosh`, `arb_atanh` | |

The constants $\pi$, $e$ and $\log 2$ are `arb_const_pi`, `arb_const_e`
and `arb_const_log2`; they are cached and cost nothing after the first
call at a given precision. Power series of the elementary functions
(Taylor expansions, for automatic differentiation) are in
[arb_poly](arb_poly), complex versions in [acb](acb).
