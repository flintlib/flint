---
title: Complex arithmetic
functions: acb_add acb_mul acb_sqr acb_inv acb_div
---

Complex arithmetic on rectangles reduces to real ball arithmetic on the
real and imaginary parts. Addition and subtraction act on the parts
separately and are as tight as the real operations. Multiplication and
division mix the parts, and since a rectangle is not closed under these
operations (the product of two rectangles is not a rectangle), the result
is a rectangle enclosing the true set, with some unavoidable
overestimation.

## Multiplication {#mul}

For $x = a + bi$ and $y = c + di$, `acb_mul` computes

$$
xy = (ac - bd) + (ad + bc) i.
$$

The midpoints of both parts are computed with `arf_complex_mul`, each
correctly rounded from the exact products (using three real
multiplications instead of four for long mantissas), and the radius of
each part is bounded from the radii of the four real products with the
product bound of [](algorithms/ball-arithmetic/error-propagation#thm-mul).
Several cases are recognised first:

- if one of $a, b, c, d$ is exactly zero, the product needs only two real
  multiplications and no addition, so real and purely imaginary factors
  are handled exactly as in real arithmetic (`acb_mul_arb` is the same
  for a real `arb_t` factor);
- if `x` and `y` are the same object, the squaring formula $(a^2 - b^2) +
  2abi$ is used, which is allowed because identical pointers mean the same
  number (see [](guide/types/ball-arithmetic#aliasing)); `acb_sqr`
  calls it;
- `acb_cube` uses two real squarings and two real multiplications.

`acb_mul_naive` computes the four real products with `arb_mul` and adds
them, which gives slightly wider results; it is kept for testing.

## Division {#div}

`acb_div` divides the parts directly by a real or purely imaginary
divisor. Otherwise it computes $x \cdot (1/y)$ with `acb_inv`, which
computes the midpoint of $1/(a + bi) = (a - bi)/(a^2 + b^2)$ with $a^2 +
b^2$ rounded once, and bounds the propagated error of each part directly
from the exact formula for the perturbation (given with `acb_inv` on the
module page). If the divisor contains zero, or one of its parts is
infinite or NaN, the result is indeterminate.

```c example name=acb_arithmetic
#include <flint/flint.h>
#include <flint/acb.h>

#define SHOW(label, x) do { flint_printf("%-15s", label); \
    acb_printd(x, 6); flint_printf("\n"); } while (0)

int main(void)
{
    acb_t x, y, z;
    mag_t e;
    slong prec = 64;

    acb_init(x); acb_init(y); acb_init(z); mag_init(e);

    /* x = [1 +/- 2^-7] + [2 +/- 2^-7] i, and an independent copy y */
    acb_set_si_si(x, 1, 2);
    mag_set_ui_2exp_si(e, 1, -7);
    acb_add_error_mag(x, e);
    acb_set(y, x);

    acb_mul(z, x, y, prec);          SHOW("x * y", z);
    acb_mul_arb(z, x, acb_realref(x), prec);  SHOW("x * re(x)", z);

    acb_div(z, x, y, prec);          SHOW("x / y", z);
    acb_inv(z, x, prec);             SHOW("1 / x", z);

    acb_set_si_si(y, 0, 0);
    mag_set_ui_2exp_si(e, 1, -3);
    acb_add_error_mag(y, e);         /* y contains 0 */
    acb_inv(z, y, prec);             SHOW("1 / [+/- 1/8]", z);

    acb_clear(x); acb_clear(y); acb_clear(z); mag_clear(e);
    flint_cleanup();
    return 0;
}
```
```output
x * y          (-3.00000 + 4.00000j)  +/-  (0.0470, 0.0470j)
x * re(x)      (1.00000 + 2.00000j)  +/-  (0.0157, 0.0235j)
x / y          (1.00000 + 0j)  +/-  (0.0123, 0.0123j)
1 / x          (0.200000 - 0.400000j)  +/-  (0.00253, 0.00253j)
1 / [+/- 1/8]  (nan + nanj)  +/-  (+inf, +infj)
```

Because `x` and `y` are different objects, $x / y$ is computed as the
quotient of two independent numbers in the same rectangle, which is a
rectangle around 1 rather than exactly 1. The product with the real ball
$\operatorname{Re} x$ uses only two real multiplications.

## Precision {#precision}

As in `arb`, `prec` bounds the length of the midpoints of the result, and
an exact result is returned exactly: products and sums of Gaussian
integers are exact if they fit in `prec` bits. For long sums of products,
`acb_dot` is faster and more accurate than repeated `acb_addmul`.
