---
title: Getting started with padic_radix
summary: A first program with radix p-adic numbers: creating a context with relative or absolute precision, exact and inexact elements, how error terms appear and propagate, status codes, and when to use padic_radix instead of padic.
order: 1
---

A `padic_radix` context is a generic ring (`gr_ctx_t`) for
$\mathbf{Q}_p$. Its elements are balls $p^v u + O(p^N)$ or exact
elements, and every operation computes an error term that is guaranteed
to contain the error of the result. The context chooses how much
precision results keep: a relative precision $r$ (digits after the
leading one), an absolute precision $a$ (digits below $p^a$), both, or
neither.

## A first program {#first}

The program works in $\mathbf{Q}_7$ with relative precision 10 and no
absolute bound.

```c example name=padic_radix_first
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/padic_radix.h>

int main(void)
{
    gr_ctx_t ctx;
    padic_radix_t x, y, z;
    int status = GR_SUCCESS;

    /* Q_7, relative precision 10 digits, no absolute precision cap */
    gr_ctx_init_padic_radix(ctx, 7, 10, PADIC_RADIX_PREC_INF, 0);
    gr_ctx_println(ctx);
    padic_radix_init(x, ctx);
    padic_radix_init(y, ctx);
    padic_radix_init(z, ctx);

    status |= padic_radix_set_ui(x, 50, ctx);
    status |= padic_radix_set_ui(y, 1, ctx);
    status |= padic_radix_sub(z, x, y, ctx);       /* exact: 49 */
    flint_printf("50 - 1  = "); gr_println(z, ctx);
    status |= padic_radix_div(z, y, z, ctx);       /* exact: 1/49 */
    flint_printf("1/49    = "); gr_println(z, ctx);

    status |= padic_radix_set_ui(x, 3, ctx);
    status |= padic_radix_inv(x, x, ctx);          /* 1/3 to 10 digits */
    flint_printf("1/3     = "); gr_println(x, ctx);

    status |= padic_radix_set_ui(z, 3, ctx);
    status |= padic_radix_mul(z, z, x, ctx);       /* 3 * (1/3) */
    status |= padic_radix_sub(z, z, y, ctx);       /* ... - 1 */
    flint_printf("3*(1/3) - 1 = "); gr_println(z, ctx);
    flint_printf("is zero: %s\n",
        padic_radix_is_zero(z, ctx) == T_UNKNOWN ? "unknown" : "decided");

    status |= padic_radix_set_ui(x, 2, ctx);
    status |= padic_radix_sqrt(z, x, ctx);         /* 2 = 3^2 mod 7 */
    flint_printf("sqrt(2) = "); gr_println(z, ctx);
    status |= padic_radix_set_ui(x, 7, ctx);
    status |= padic_radix_exp(z, x, ctx);
    status |= padic_radix_log(y, z, ctx);
    flint_printf("log(exp(7)) = "); gr_println(y, ctx);
    flint_printf("status = %d\n", status);

    padic_radix_clear(x, ctx);
    padic_radix_clear(y, ctx);
    padic_radix_clear(z, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```
```output
Radix 7-adic numbers (rel prec 10, abs prec inf)
50 - 1  = (1) * 7^2
1/49    = (1) * 7^-2
1/3     = 188316833 + O(7^10)
3*(1/3) - 1 = 0 + O(7^10)
is zero: unknown
sqrt(2) = 266983762 + O(7^10)
log(exp(7)) = (1) * 7^1 + O(7^10)
status = 0
```

The integers 50 and 1 fit in 10 digits, so they and their difference
$49 = 7^2$ are exact, and so is $1/49$, whose expansion terminates. The
inverse of 3 has an infinite expansion and is truncated to 10 digits,
which the error term $O(7^{10})$ records. Multiplying it by 3 and
subtracting 1 cancels every known digit: the result is $0 + O(7^{10})$,
a ball that contains zero and nonzero numbers, so `padic_radix_is_zero`
answers `T_UNKNOWN`, and dividing by it would return `GR_UNABLE`. In
`padic` the same computation would return an exact-looking zero. The
square root of 2 and the exponential and logarithm are computed to the
relative precision of the context.

## Choosing the precision {#precision}

| `prec_rel` | `prec_abs` | Behaviour |
|---|---|---|
| $r$ | `PADIC_RADIX_PREC_INF` | floating point: every result keeps $r$ digits after its leading digit |
| `PADIC_RADIX_PREC_INF` | $a$ | fixed point: results are computed modulo $p^a$, like `padic` with precision $a$ |
| $r$ | $a$ | both caps apply |
| `PADIC_RADIX_PREC_INF` | `PADIC_RADIX_PREC_INF` | exact arithmetic in $\Z[1/p]$; operations with an infinite result return `GR_UNABLE` |

In every configuration, inexact inputs limit the precision of the
result further, as in the example. The precision of a result can be
read with `padic_radix_get_error`, and `padic_radix_is_exact` tells
whether it is exact.

## Status codes {#status}

The element functions return `GR_SUCCESS`, or `GR_DOMAIN` when the
operation is undefined (division by an exact zero, the square root of a
non-square, the logarithm of an element not congruent to 1 modulo $p$),
or `GR_UNABLE` when the answer is not determined by the inputs (division
by a ball containing zero) or cannot be represented (an infinite
expansion in an exact context, valuations beyond
`PADIC_RADIX_ERR_MAX`). Generic code should check the status as for any
other ring; `GR_MUST_SUCCEED` aborts on failure.

## Generic structures {#generic}

Since the context is a `gr_ctx_t`, the generic modules work over it:
`gr_poly` for polynomials and power series arithmetic, `gr_mat` for
matrices (products use the dot products of this module, determinants the
generic algorithms for fields), and `gr_series` for power series with
error tracking in both the $p$-adic and the $x$-adic direction. Compared
with `padic_poly` and `padic_mat`, every coefficient has its own
valuation and precision instead of sharing one, and precision loss in
division and elimination is detected.
