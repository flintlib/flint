---
title: Computing bounds with mag
summary: How to combine upper and lower bounds correctly with mag_t, with an example, and when to use arb_t instead.
order: 1
---

A `mag_t` is a one-sided quantity: each function either rounds up (the
plain name) or rounds down (suffix `_lower`), and there is no
round-to-nearest. To compute a rigorous upper bound for an expression,
every intermediate quantity must be rounded in the direction that makes
the final result larger. For an increasing operation that is the same
direction as the result; for a decreasing one, such as a denominator or
the argument of $e^{-x}$, it is the opposite direction.

## Monotonicity decides the direction {#direction}

To bound $(x + 1)/(y + 1)$ from above, for nonnegative $x$ and $y$, the
numerator must be rounded up and the denominator rounded down, starting
from an upper bound for $x$ and a lower bound for $y$. For a lower bound
everything is reversed. Forgetting this is the most common mistake when
writing bounds with `mag_t`: the results look plausible but are not
bounds, and the error only shows up for rare inputs.

```c example name=mag_bounds_quotient
#include <flint/flint.h>
#include <flint/mag.h>

/* Bounds for (x + 1) / (y + 1), given bounds xlo <= x <= xhi
   and ylo <= y <= yhi for nonnegative x and y. */
int main(void)
{
    mag_t xlo, xhi, ylo, yhi, num, den, upper, lower;

    mag_init(xlo); mag_init(xhi); mag_init(ylo); mag_init(yhi);
    mag_init(num); mag_init(den); mag_init(upper); mag_init(lower);

    mag_set_d_lower(xlo, 2.5);  mag_set_d(xhi, 2.75);   /* x in [2.5, 2.75] */
    mag_set_ui_lower(ylo, 9);   mag_set_ui(yhi, 10);    /* y in [9, 10] */

    /* upper bound: numerator up, denominator down */
    mag_add_ui(num, xhi, 1);
    mag_add_ui_lower(den, ylo, 1);
    mag_div(upper, num, den);

    /* lower bound: numerator down, denominator up */
    mag_add_ui_lower(num, xlo, 1);
    mag_add_ui(den, yhi, 1);
    mag_div_lower(lower, num, den);

    flint_printf("lower = "); mag_printd(lower, 10); flint_printf("\n");
    flint_printf("upper = "); mag_printd(upper, 10); flint_printf("\n");

    /* tail of the series for exp(1/2) after 10 terms */
    mag_set_ui_2exp_si(num, 1, -1);
    mag_exp_tail(upper, num, 10);
    flint_printf("exp tail = "); mag_printd(upper, 5); flint_printf("\n");

    mag_clear(xlo); mag_clear(xhi); mag_clear(ylo); mag_clear(yhi);
    mag_clear(num); mag_clear(den); mag_clear(upper); mag_clear(lower);
    return 0;
}
```
```output
lower = 0.3181818170
upper = 0.3750000014
exp tail = 5.3823e-10
```

The exact range is $[3.5/11, 3.75/10] = [0.31818\ldots, 0.375]$. Even
inputs such as 2.75 and 10, which are exactly representable, come out
slightly enlarged: `mag_set_d` and `mag_add_ui` always add a unit in the
last place rather than checking whether rounding was needed. This is
deliberate; with only 30 bits it costs nothing in practice and keeps the
code branch-free. The last line uses the tail bound
$\sum_{k \ge N} x^k/k! \le 2 x^N / N!$ for $N \ge 2x$, about twice the
true value $2.8 \cdot 10^{-10}$.

## Special values {#special}

The functions treat $+\infty$ as an ordinary value, so a bound that
overflows, or a division by a lower bound that is zero, simply propagates
$+\infty$ and the result is still valid. Upper-bound functions return
$+\infty$ for meaningless operations such as $0 \cdot \infty$, and
lower-bound functions return 0 for them. Code that computes a bound
therefore rarely needs special cases, but it should check
`mag_is_finite` on the result when an infinite bound means that the
computation must be redone with other parameters.

## When to use arb_t instead {#arb}

`mag_t` is the right type for error bounds inside an algorithm: the
radius of a ball, the tail of a series, a bound on a derivative. Its
arithmetic is fast because it never needs more than one word of
mantissa. It is inconvenient when the expression involves subtraction of
quantities that are not ordered, or long chains where the rounding
direction alternates. In that case, compute an enclosure with `arb_t`
balls at a low precision such as 30 bits and convert it with
`arb_get_mag` (upper bound for the absolute value) or `arb_get_mag_lower`
(lower bound). The radius of the ball then takes care of the rounding
directions automatically.

The representation of `mag_t` and the proofs that its operations round in
the right direction are described in [](algorithms/ball-arithmetic/representation).
