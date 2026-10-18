---
title: Correctly rounded elementary and special functions
functions: decfloat_pi decfloat_exp decfloat_via_arb
---

The elementary functions of `decfloat_t` (and, through the generic
methods, the special functions such as `gr_gamma`, `gr_zeta`,
`gr_erf` and `gr_bessel_j` with a decimal context) are correctly
rounded: the result is the exact function value rounded once to the
context precision in the context rounding mode. The functions are
evaluated with `arb` and rounded by Ziv's strategy [@Ziv1991], with the
cases where that strategy cannot succeed handled separately.

## Ziv's strategy {#ziv}

The argument is converted to an `arb` ball exactly, the function is
evaluated at a working precision somewhat above the target, and the
resulting ball is rounded: if both of its endpoints round to the same
`prec`-digit number, that number is the correctly rounded result. If
not, the working precision is increased and the evaluation repeated.
This terminates unless the function value is a representable number or
lies extremely close to one; the loop gives up with `GR_UNABLE` at 64
times the target precision (for the elementary functions, plus twice
the number of digits of the argument, since an argument with more
digits than the target can still matter near a cancellation, as in
$\log_2(1 + 10^{-10^6})$).

## Cases decided without Ziv's strategy {#special}

Three kinds of arguments would defeat Ziv's strategy systematically, and
are detected before any evaluation, from the sign, the exponent and a
few digits of the argument, in time independent of its length:

- **Exact values**, such as $\Gamma(n) = (n - 1)!$, $\log_2(2^n)$,
  $\sin(\pi k/4)$ for the functions of $\pi x$, $\zeta(-2k) = 0$,
  $\operatorname{acos}(1) = 0$, $\operatorname{asin}(1/2)/\pi = 1/6$,
  and powers $x^{p/q}$ with exact results.
- **Tiny arguments**, where $f(x) = S + t$ with an exact decimal leading
  term $S$ ($x$, $1$, $\pm 1/x$, $x - 1$, ...) and a tail $t$ of known
  sign far below the rounding position: the result is $S$ rounded with
  $t$ as a sticky bit (`_decfloat_round_with_tail`). For instance
  $\sin(10^{-10^9})$, $\Gamma(10^{-10^9})$ and $\zeta(1 + 10^{-10^9})$
  are correctly rounded at negligible cost.
- **Large arguments** of functions converging exponentially to a
  representable limit ($\tanh$, $\operatorname{erf}$,
  $\operatorname{expm1}$, $\zeta$) or behaving like $1/x$
  ($\operatorname{acot}$, ...), handled in the same way.

Points outside the real domain (such as $\log(-1)$ or
$\operatorname{acosh}(1/2)$) give `GR_DOMAIN` before any evaluation, and
arguments whose evaluation would be unreasonably expensive (non-elementary
functions at arguments of magnitude $10^{10^6}$ and beyond, for example)
give `GR_UNABLE`. The special functions other than those listed rely on
Ziv's strategy alone and could in principle return `GR_UNABLE` at an
argument where the value is an undetected representable number.

The ball versions (`decball_exp` and so on) instead evaluate the
function once on the `arb` ball enclosing the argument and convert the
result back, and the complex versions use `acb`
([](decimal/index#complex-functions)).

## Example {#example}

```c example name=decimal_functions
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/gr_special.h>
#include <flint/decimal.h>

int main(void)
{
    gr_ctx_t ctx;
    decfloat_t x, y;
    int status = GR_SUCCESS;

    gr_ctx_init_decfloat(ctx, 25, 0);
    decfloat_init(x, ctx);
    decfloat_init(y, ctx);

    status |= decfloat_set_ui(x, 1, ctx);
    status |= decfloat_exp(y, x, ctx);
    flint_printf("exp(1)         = "); gr_println(y, ctx);
    status |= decfloat_atan(y, x, ctx);
    flint_printf("atan(1)        = "); gr_println(y, ctx);

    /* exact values are recognised */
    status |= decfloat_set_ui(x, 5, ctx);
    status |= gr_gamma(y, x, ctx);
    flint_printf("gamma(5)       = "); gr_println(y, ctx);
    status |= decfloat_set_str(x, "1024", ctx);
    status |= gr_log2(y, x, ctx);
    flint_printf("log2(1024)     = "); gr_println(y, ctx);
    status |= decfloat_set_str(x, "0.25", ctx);
    status |= gr_sin_pi(y, x, ctx);
    flint_printf("sin(pi/4)      = "); gr_println(y, ctx);
    status |= decfloat_set_str(x, "0.5", ctx);
    status |= gr_sin_pi(y, x, ctx);
    flint_printf("sin(pi/2)      = "); gr_println(y, ctx);

    /* a tiny argument in directed rounding */
    decimal_ctx_set_rnd(ctx, DECIMAL_RND_DOWN);
    status |= decfloat_set_str(x, "1e-1000000000", ctx);
    status |= decfloat_sin(y, x, ctx);
    flint_printf("sin(1e-10^9)   = "); gr_println(y, ctx);
    status |= decfloat_cos(y, x, ctx);
    flint_printf("cos(1e-10^9)   = "); gr_println(y, ctx);

    /* outside the real domain */
    status |= decfloat_set_si(x, -1, ctx);
    flint_printf("log(-1) domain error: %d\n", decfloat_log(y, x, ctx) == GR_DOMAIN);

    flint_printf("status %d\n", status);
    decfloat_clear(x, ctx);
    decfloat_clear(y, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```

```output
exp(1)         = 2.718281828459045235360287
atan(1)        = 0.7853981633974483096156608
gamma(5)       = 24
log2(1024)     = 10
sin(pi/4)      = 0.7071067811865475244008444
sin(pi/2)      = 1
sin(1e-10^9)   = 9.999999999999999999999999e-1000000001
cos(1e-10^9)   = 0.9999999999999999999999999
log(-1) domain error: 1
status 0
```

In the directed rounding mode toward zero, $\sin(10^{-10^9})$ is the
largest 25-digit number below $10^{-10^9}$ and $\cos(10^{-10^9})$ the
largest below 1.
