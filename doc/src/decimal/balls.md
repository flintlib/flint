---
title: Decimal ball arithmetic
summary: Decimal balls with a correctly rounded midpoint and a short decimal radius, how radii are computed and printed, and how balls relate to arb.
order: 3
---

A `decball_t` is a ball $[m \pm r]$ with a `decfloat_t` midpoint $m$ and
a `decmag_t` radius $r$, and represents every real number within $r$ of
$m$. Ball operations return balls that contain the exact result for
every choice of points in the input balls. `deccball_t` is the complex
version, a pair of real balls (a rectangle, as in `acb`).

## Midpoints and radii {#mid-rad}

The midpoint of a result is computed by the corresponding
floating-point operation, correctly rounded to the context precision,
and the radius is the propagated input radii plus the rounding error of
the midpoint. The rounding error is the exact one, rounded up to the
radius precision, rather than a whole ulp (unless the flag
`DECIMAL_SLOPPY_RADIUS` is set), so that an exact input often gives a
radius well below one ulp of the midpoint.

A radius is a `decmag_t`: a mantissa of `rad_prec` digits in a single
word (by default 4, at most `DECMAG_MAX_PREC` = 9 on 64-bit machines)
with an `fmpz` exponent, or zero or $+\infty$. All radius operations
round in the direction that preserves the bound, and reduce to a few
word operations. An infinite radius denotes the whole real line; the
midpoint of a ball is always finite.

The functions `decball_rel_accuracy_digits` and `decball_trim` measure
and remove the digits of the midpoint that the radius makes
meaningless; `decball_mid`, `decball_rad`, `decball_lower` and
`decball_upper` extract parts and endpoints as exact balls.

## Balls and arb {#arb}

The basic arithmetic is computed in decimal. The elementary and special
functions are computed by a rigorous roundtrip through `arb`: the input
balls are converted to `arb` balls, the function is evaluated, and the
result is converted back with the conversion errors added to the radius
(`decball_get_arb`, `decball_set_arb`). The conversion from `arb` is
tight: the midpoint is correctly rounded and the exact rounding error is
added to the radius.

## Strings {#strings}

Balls print as `[mid +/- rad]`, or as `mid` when exact, with all digits
of the midpoint and the radius. `decball_set_str` parses this format and
general expressions with the `+/-` operator, also inside strings for
polynomials and matrices over balls.

## Example {#example}

```c example name=decimal_balls
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/decimal.h>

int main(void)
{
    gr_ctx_t ctx;
    decball_t x, y;
    int status = GR_SUCCESS;

    gr_ctx_init_decball(ctx, 20, 0);
    decball_init(x, ctx);
    decball_init(y, ctx);

    /* a measured quantity */
    status |= decball_set_str(x, "[0.785398 +/- 1e-6]", ctx);
    status |= decball_mul(y, x, x, ctx);
    flint_printf("x^2      = "); gr_println(y, ctx);
    status |= decball_sin(y, x, ctx);
    flint_printf("sin(x)   = "); gr_println(y, ctx);
    flint_printf("digits   = %wd\n", decball_rel_accuracy_digits(y, ctx));
    status |= decball_trim(y, y, ctx);
    flint_printf("trimmed  = "); gr_println(y, ctx);

    /* exact input: only rounding errors */
    status |= decball_set_ui(x, 2, ctx);
    status |= decball_sqrt(y, x, ctx);
    flint_printf("sqrt(2)  = "); gr_println(y, ctx);
    status |= decball_pi(y, ctx);
    flint_printf("pi       = "); gr_println(y, ctx);

    /* sqrt of a ball containing negative numbers */
    status |= decball_set_str(x, "[0 +/- 1e-10]", ctx);
    flint_printf("sqrt([0 +/- 1e-10]) unable: %d\n",
        decball_sqrt(y, x, ctx) == GR_UNABLE);

    flint_printf("status %d\n", status);
    decball_clear(x, ctx);
    decball_clear(y, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```

```output
x^2      = [0.616850018404 +/- 1.572e-6]
sin(x)   = [0.70710666564674062862 +/- 7.082e-7]
digits   = 6
trimmed  = [0.70710666565 +/- 7.083e-7]
sqrt(2)  = [1.4142135623730950488 +/- 1.689e-21]
pi       = [3.1415926535897932385 +/- 3.737e-20]
sqrt([0 +/- 1e-10]) unable: 1
status 0
```

The square root of a ball containing negative numbers returns
`GR_UNABLE`, since the result is not a real ball; a ball lying entirely
below zero returns `GR_DOMAIN`.
