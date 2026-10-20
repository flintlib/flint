---
title: Rounding, precision and the inexact flag
summary: How arf_t functions round their results, what the precision argument and the return value mean, and how to compute exactly or with a single rounding.
order: 1
---

Every `arf_t` function that can produce a result not representable with
`prec` bits takes a precision `prec` and a rounding mode `rnd`, and
returns an `int`. The result is the exact mathematical result rounded
once to a mantissa of at most `prec` bits, and the return value says
whether that rounding changed anything.

## Rounding modes {#modes}

The program below computes $-1/3$ with an 8-bit mantissa in each strict
rounding mode, and then a quotient that is exactly representable.

```c example name=arf_rounding_modes
#include <flint/flint.h>
#include <flint/arf.h>

int main(void)
{
    arf_t x, y, z;
    arf_rnd_t modes[5] = { ARF_RND_DOWN, ARF_RND_UP, ARF_RND_FLOOR,
                           ARF_RND_CEIL, ARF_RND_NEAR };
    const char * names[5] = { "DOWN ", "UP   ", "FLOOR", "CEIL ", "NEAR " };
    int i, inexact;

    arf_init(x); arf_init(y); arf_init(z);

    /* -1/3 rounded to 8 bits in each direction */
    arf_set_si(x, -1);
    arf_set_ui(y, 3);
    for (i = 0; i < 5; i++)
    {
        inexact = arf_div(z, x, y, 8, modes[i]);
        flint_printf("%s ", names[i]); arf_print(z);
        flint_printf("  inexact = %d\n", inexact);
    }

    /* an exact result: 3/4 at 8 bits */
    arf_set_ui(x, 3);
    inexact = arf_div_ui(z, x, 4, 8, ARF_RND_NEAR);
    flint_printf("3/4: "); arf_print(z); flint_printf("  inexact = %d\n", inexact);

    arf_clear(x); arf_clear(y); arf_clear(z);
    return 0;
}
```
```output
DOWN  (-85 * 2^-8)  inexact = 1
UP    (-171 * 2^-9)  inexact = 1
FLOOR (-171 * 2^-9)  inexact = 1
CEIL  (-85 * 2^-8)  inexact = 1
NEAR  (-171 * 2^-9)  inexact = 1
3/4: (3 * 2^-2)  inexact = 0
```

`arf_print` shows the exact value with an odd mantissa, so a rounded
result can print with fewer than `prec` bits: $-85/256$ has the 8-bit
mantissa $-170 \cdot 2^{-9}$ with a trailing zero. `ARF_RND_DOWN` and
`ARF_RND_UP` round towards and away from zero; `ARF_RND_FLOOR` and
`ARF_RND_CEIL` towards $-\infty$ and $+\infty$, so for negative numbers
FLOOR agrees with UP. `ARF_RND_NEAR` rounds to the nearest representable
number and breaks ties by choosing an even mantissa.

Ball arithmetic uses `ARF_RND_DOWN` for the midpoint and adds the
rounding error to the radius, because truncation is the cheapest mode:
the discarded limbs need not be inspected beyond deciding whether they
are all zero.

## The precision is a bound on the mantissa {#precision}

An `arf_t` has no precision of its own. The result of an operation has
*at most* `prec` bits: an exact result with fewer bits is stored with
fewer bits, which keeps integers and simple fractions small and makes
operations on them fast. The precision passed to one operation does not
affect the next one; `arf_set` copies all bits of its input, and
`arf_set_round` rounds to a given precision.

For the same reason a result can need far more bits than its inputs.
`arf_add` of $2^{100}$ and $1$ is exact at 101 bits and inexact at 100
bits. The special precision `ARF_PREC_EXACT` asks for the exact result.
It is safe only when that result is known to be small: the exact sum of
$1$ and $2^{-10^{12}}$ has $10^{12}$ bits. To sum many terms with a single
rounding, use `arf_sum` (or `arf_dot` for products), which never forms
such exact intermediate results.

## Exactness and the return value {#inexact}

The return value is 0 exactly when the result is the exact value, and 1
otherwise. Ball arithmetic relies on it: `arb_add` adds a rounding error
to the radius only when `arf_add` reports one, so that exact inputs with
an exact result give an exact ball. The flag is computed from the
discarded bits, or, in division, by multiplying the quotient back when
its low bits are zero.

The relaxed modes `ARF_RND_FAST` and `ARF_RND_ACCURATE` of division and
square roots trade correct rounding for speed at high precision: they
allow an error of up to one unit (respectively 0.51 units) in the last
place. Their return value is still the inexact flag: it is 0 exactly
when the result equals the exact value.

## Special values {#special}

Infinities and NaN follow IEEE 754 where it applies: $\infty - \infty$,
$0 \cdot \infty$ and $0/0$ give NaN. Since there is no signed zero,
$1/0$ is NaN rather than an infinity. Comparisons with NaN return 0, as
if the numbers were equal, and `arf_equal` considers NaN equal to itself;
test with `arf_is_nan` when NaN is possible. Results that are special
values are exact and return 0.
