---
title: Getting started with decimal arithmetic
summary: The four decimal rings and their contexts, precision in digits, rounding modes, the exact ring of decimal fractions, and how the decimal types work with the generic rings interface.
order: 1
---

The `decimal` module provides four rings, each described by a context:

| Context | Elements | Results |
|---|---|---|
| `gr_ctx_init_decfloat` | `decfloat_t`, real numbers $\pm M \cdot 10^k$ | correctly rounded |
| `gr_ctx_init_deccfloat` | `deccfloat_t`, complex numbers | correctly rounded per component |
| `gr_ctx_init_decball` | `decball_t`, real balls $[m \pm r]$ | rigorous enclosures |
| `gr_ctx_init_deccball` | `deccball_t`, complex rectangular balls | rigorous enclosures |

A context is created with a precision in decimal digits and flags, used
for every operation on its elements, and cleared with `gr_ctx_clear`.
Its settings can be changed at any time with `decimal_ctx_set_prec`,
`decimal_ctx_set_rnd`, `decimal_ctx_set_rad_prec`,
`decimal_ctx_set_exp_limits` and `decimal_ctx_set_flags`; elements are
not tied to the settings with which they were computed, and an
operation only rounds its result.

## Precision and rounding {#rounding}

The precision is a number of significant decimal digits. Every
floating-point operation (arithmetic, conversions, and the elementary
and special functions) returns the exact result rounded once to that
many digits in the context's rounding mode: toward zero, away from
zero, toward $\pm\infty$, or to nearest with ties to even, away from
zero or toward zero (`DECIMAL_RND_DOWN`, ..., `DECIMAL_RND_NEAR_ZERO`).
Precision and exponents are unbounded. Decimal constants
such as $0.1$ are represented exactly, so decimal inputs and outputs
need no binary conversion: in the example below, $0.1 + 0.2 = 0.3$
holds exactly. The `_round` versions of the arithmetic functions take
the precision and rounding mode as arguments instead.

Exponents are unbounded by default. Limits on the scientific exponent
can be set with `decimal_ctx_set_exp_limits`; the flags
`DECIMAL_ALLOW_INF`, `DECIMAL_ALLOW_NAN` and `DECIMAL_ALLOW_UNDERFLOW`
then decide whether an overflow gives an infinity, an invalid operation
a NaN, and an underflow zero, instead of returning an error status.

The precision `DECIMAL_PREC_EXACT` gives the exact ring
$\Z[1/10]$ of decimal fractions: nothing is rounded, and an operation
whose result is not a decimal fraction, such as $1/3$, returns
`GR_UNABLE`.

## Status codes {#status}

Every function returns `GR_SUCCESS`, `GR_DOMAIN` (the result is
undefined, such as $\log(-1)$ in the real rings or $1/0$ without
infinities) or `GR_UNABLE` (the result could not be computed, such as an
inexact result in the exact ring, a NaN operand without
`DECIMAL_ALLOW_NAN`, or a ball whose sign cannot be decided). Results
are only meaningful after `GR_SUCCESS`; the usual idiom is to
accumulate `status |= ...` and check once.

## The generic interface {#generic}

The decimal rings are `gr` rings, so the generic functions and
structures work with them: `gr_println` prints an element, `gr_gamma`,
`gr_zeta`, `gr_bessel_j` and the other special functions dispatch to
correctly rounded decimal implementations, and polynomials (`gr_poly`),
matrices (`gr_mat`) and power series over the decimal rings are
available. Generic algorithms compose several rounded operations, so
their results (a matrix product, a polynomial evaluation) are not
correctly rounded as a whole; see [](decimal/arithmetic#generic).

## Example {#example}

```c example name=decimal_getting_started
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/decimal.h>

int main(void)
{
    gr_ctx_t ctx;
    decfloat_t x, y, z;
    int status = GR_SUCCESS;

    gr_ctx_init_decfloat(ctx, 30, 0);
    decfloat_init(x, ctx);
    decfloat_init(y, ctx);
    decfloat_init(z, ctx);

    gr_ctx_println(ctx);

    /* 0.1 + 0.2 = 0.3 exactly in decimal arithmetic */
    status |= decfloat_set_str(x, "0.1", ctx);
    status |= decfloat_set_str(y, "0.2", ctx);
    status |= decfloat_add(z, x, y, ctx);
    status |= decfloat_set_str(x, "0.3", ctx);
    flint_printf("0.1 + 0.2 == 0.3: %d\n", decfloat_equal(z, x, ctx) == T_TRUE);

    /* 2/3 in three rounding modes */
    status |= decfloat_set_ui(x, 2, ctx);
    status |= decfloat_div_ui(y, x, 3, ctx);
    flint_printf("near:  "); gr_println(y, ctx);
    decimal_ctx_set_rnd(ctx, DECIMAL_RND_DOWN);
    status |= decfloat_div_ui(y, x, 3, ctx);
    flint_printf("down:  "); gr_println(y, ctx);
    decimal_ctx_set_rnd(ctx, DECIMAL_RND_UP);
    status |= decfloat_div_ui(y, x, 3, ctx);
    flint_printf("up:    "); gr_println(y, ctx);
    decimal_ctx_set_rnd(ctx, DECIMAL_RND_NEAR);

    /* elementary functions are correctly rounded */
    status |= decfloat_sqrt(y, x, ctx);
    flint_printf("sqrt(2) = "); gr_println(y, ctx);
    status |= decfloat_set_str(x, "1e-1000", ctx);
    status |= decfloat_exp(y, x, ctx);
    flint_printf("exp(1e-1000), near = "); gr_println(y, ctx);
    decimal_ctx_set_rnd(ctx, DECIMAL_RND_UP);
    status |= decfloat_exp(y, x, ctx);
    flint_printf("exp(1e-1000), up   = "); gr_println(y, ctx);

    /* the exact ring of decimal fractions */
    decimal_ctx_set_prec(ctx, DECIMAL_PREC_EXACT);
    status |= decfloat_set_ui(x, 1, ctx);
    flint_printf("1/8 exact: %d, ", decfloat_div_ui(y, x, 8, ctx) == GR_SUCCESS);
    gr_println(y, ctx);
    flint_printf("1/3 exact: %d\n", decfloat_div_ui(y, x, 3, ctx) == GR_SUCCESS);

    flint_printf("status %d\n", status);
    decfloat_clear(x, ctx);
    decfloat_clear(y, ctx);
    decfloat_clear(z, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```

```output
Decimal floating-point numbers (prec 30, rnd near)
0.1 + 0.2 == 0.3: 1
near:  0.666666666666666666666666666667
down:  0.666666666666666666666666666666
up:    0.666666666666666666666666666667
sqrt(2) = 1.41421356237309504880168872421
exp(1e-1000), near = 1
exp(1e-1000), up   = 1.00000000000000000000000000001
1/8 exact: 1, 0.125
1/3 exact: 0
status 0
```

The exponential of $10^{-1000}$ is $1 + 10^{-1000} + \cdots$; rounded to
nearest at 30 digits it is 1, and rounded up it is the next number
above 1. Ziv's strategy alone could not decide this rounding, since the
value lies $10^{-1000}$ above a representable number; the function
detects the tiny argument instead ([](decimal/functions)).
