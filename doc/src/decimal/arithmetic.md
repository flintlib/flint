---
title: Correctly rounded arithmetic
functions: decfloat_neg decfloat_neg_round _decfloat_add _decfloat_set_round_limbs decfloat_vec_dot decimal_rounding_info
---

The arithmetic operations on `decfloat_t` return the exact result of
the operation on the exact operands, rounded once to the precision in
the rounding mode; the result is exact whenever it is representable.
The plain functions use the context's precision and rounding mode, the
`_round` functions take them as arguments, and the underscore functions
additionally report what the rounding did.

## Rounding {#rounding}

A result is rounded to `prec` significant digits, counted from its
leading nonzero digit. Rounding never depends on how the operands were
stored: the trailing zero digits of a mantissa's lowest limb, which are
an artefact of the limb alignment, do not count as significant. The
seven rounding modes are toward zero, away from zero, toward
$-\infty$, toward $+\infty$, and to nearest with ties to even, away
from zero or toward zero. `decfloat_set_round_info` and
the underscore functions report in a `decimal_rounding_info` whether the
result is inexact, whether its magnitude was rounded up, and whether
the exponent limits caused an underflow or overflow, and optionally
write an upper bound for the rounding error to a `decmag_t`; the ball
arithmetic adds this bound to the radius, which is why ball radii carry
the exact rounding error rather than a whole ulp.

## Algorithms {#algorithms}

The mantissas are integers in the radix $B = 10^{19}$ (on 64-bit
machines) and the exponents count limbs, so an addition aligns its
operands by whole limbs and a rounding drops whole limbs and clears
digits in one boundary limb. An addition whose smaller operand lies
entirely below the rounding position of the larger is computed as the
larger plus a sticky tail, in $O(\mathit{prec})$ time whatever the
difference of exponents. A product is first computed as a truncated
high product with an error bound; if both ends of the resulting interval
round to the same value, that is the correctly rounded result, and
otherwise (rarely) the full product is formed. Division and square root
compute a truncated quotient or root with at least `prec + 1` digits and
use the remainder as a sticky bit. See
[](algorithms/floating-point/decimal-arithmetic).

## Composite operations {#generic}

Operations reached through the generic fallbacks of the `gr` interface
are not correctly rounded, since they compose several rounded
operations. This includes `decfloat_vec_dot` and hence matrix and
polynomial arithmetic, `gr_pow_fmpq`, `gr_hypot`, `gr_fac`, `gr_bin`,
power series functions, and `gr_set_str` applied to an expression such
as `"1/3 + 2/7"`. The complex arithmetic of `deccfloat_t`, on the other
hand, rounds each component of a product or quotient once
(`_decfloat_dot2`).

## Example {#example}

The example rounds a tie in three modes, adds numbers whose exponents
differ by $10^{12}$, and asks for the rounding error of $1/7$ at three
digits.

```c example name=decimal_arithmetic
#include <flint/flint.h>
#include <flint/gr.h>
#include <flint/decimal.h>

int main(void)
{
    gr_ctx_t ctx;
    decfloat_t x, y, z;
    decmag_t err;
    decimal_rounding_info info;
    int status = GR_SUCCESS;

    gr_ctx_init_decfloat(ctx, 10, 0);
    decfloat_init(x, ctx);
    decfloat_init(y, ctx);
    decfloat_init(z, ctx);
    _decmag_init(err, ctx);

    /* ties: 0.125 to two digits */
    status |= decfloat_set_str(x, "0.125", ctx);
    status |= decfloat_set_round(y, x, 2, DECIMAL_RND_NEAR, ctx);
    flint_printf("near:      "); gr_println(y, ctx);
    status |= decfloat_set_round(y, x, 2, DECIMAL_RND_NEAR_AWAY, ctx);
    flint_printf("near_away: "); gr_println(y, ctx);
    status |= decfloat_set_round(y, x, 2, DECIMAL_RND_FLOOR, ctx);
    flint_printf("floor:     "); gr_println(y, ctx);

    /* a sum of operands with very different exponents costs O(prec) */
    status |= decfloat_set_ui(x, 1, ctx);
    status |= decfloat_set_str(y, "1e-1000000000000", ctx);
    status |= decfloat_add_round(z, x, y, 10, DECIMAL_RND_UP, ctx);
    flint_printf("1 + 1e-10^12, up:  "); gr_println(z, ctx);
    status |= decfloat_add_round(z, x, y, 10, DECIMAL_RND_NEAR, ctx);
    flint_printf("1 + 1e-10^12, near: "); gr_println(z, ctx);

    /* the rounding error of 1/7 at 10 digits */
    status |= decfloat_set_ui(x, 1, ctx);
    status |= decfloat_div_ui(y, x, 7, ctx);
    status |= decfloat_set_round_info(z, y, 3, DECIMAL_RND_NEAR, &info, err, ctx);
    flint_printf("1/7 = "); gr_println(y, ctx);
    flint_printf("to 3 digits: "); gr_println(z, ctx);
    flint_printf("inexact %d, increased %d, error <= ", info.inexact, info.increased);
    {
        char * t = _decmag_get_str(err, ctx);
        flint_printf("%s\n", t);
        flint_free(t);
    }

    flint_printf("status %d\n", status);
    decfloat_clear(x, ctx);
    decfloat_clear(y, ctx);
    decfloat_clear(z, ctx);
    _decmag_clear(err, ctx);
    gr_ctx_clear(ctx);
    return 0;
}
```

```output
near:      0.12
near_away: 0.13
floor:     0.12
1 + 1e-10^12, up:  1.000000001
1 + 1e-10^12, near: 1
1/7 = 0.1428571429
to 3 digits: 0.143
inexact 1, increased 1, error <= 1.429e-4
status 0
```
