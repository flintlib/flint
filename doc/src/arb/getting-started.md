---
title: Getting started with arb
summary: A first program with real balls, how to compile it, initialisation and the precision argument, how to read the results, and how to get numbers out as integers, doubles or strings.
order: 1
---

An `arb_t` is a real number represented as a ball: a midpoint and an
error bound. The user guide [](guide/types/ball-arithmetic) explains the
concepts in more depth.

## A first program {#first}

The following program computes $e^{\pi \sqrt{163}}$, which is famously
close to an integer, at three precisions, and then finds the nearest
integer and the distance to it.

```c example name=arb_getting_started
#include <flint/flint.h>
#include <flint/arb.h>

int main(void)
{
    arb_t x, y;
    fmpz_t n;
    slong prec;

    arb_init(x);
    arb_init(y);
    fmpz_init(n);

    /* exp(pi sqrt(163)) at increasing precision */
    for (prec = 32; prec <= 128; prec *= 2)
    {
        arb_sqrt_ui(x, 163, prec);
        arb_const_pi(y, prec);
        arb_mul(x, x, y, prec);
        arb_exp(x, x, prec);
        flint_printf("prec %3wd: ", prec);
        arb_printn(x, 35, 0);
        flint_printf("\n");
    }

    /* x is not an integer, but we can find the nearest one */
    if (!arb_get_unique_fmpz(n, x))
        flint_printf("x does not contain a unique integer\n");
    arf_get_fmpz(n, arb_midref(x), ARF_RND_NEAR);
    arb_sub_fmpz(y, x, n, 128);
    flint_printf("n = "); fmpz_print(n); flint_printf("\n");
    flint_printf("x - n = "); arb_printn(y, 10, 0); flint_printf("\n");
    flint_printf("as a double: %.17g\n", arf_get_d(arb_midref(y), ARF_RND_NEAR));

    arb_clear(x);
    arb_clear(y);
    fmpz_clear(n);
    flint_cleanup();
    return 0;
}
```
```output
prec  32: [2.625374e+17 +/- 1.70e+10]
prec  64: [2.6253741264076874e+17 +/- 5.40]
prec 128: [262537412640768743.99999999999925007 +/- 2.65e-18]
x does not contain a unique integer
n = 262537412640768744
x - n = [-7.49927e-13 +/- 5.70e-19]
as a double: -7.4992746048967683e-13
```

With FLINT installed under a prefix such as `/usr/local`, the program is
compiled with

```
cc first.c -o first -lflint -lmpfr -lgmp
```

(add `-I<prefix>/include -L<prefix>/lib` for another prefix; see
[](guide/getting-started/building)). The ball arithmetic modules are part
of the FLINT library, so no other library is needed; their headers are
`<flint/arb.h>`, `<flint/acb.h>` and so on.

## Reading the program {#reading}

The interface follows the conventions of FLINT and GMP.

- Every `arb_t` is initialised with `arb_init`, which sets it to 0, and
  freed with `arb_clear`. Variables are passed by reference; there is no
  need to take addresses.
- Output arguments come first, then inputs, then the precision `prec` in
  bits. Any output may be the same variable as an input, as in
  `arb_mul(x, x, y, prec)`.
- The precision is the precision of the *midpoint*: the result is rounded
  to `prec` bits and the rounding error is added to the radius. It is not a
  promise about the accuracy of the result. Here the argument of the
  exponential is about 40, and its rounding error, a few units in the
  last place, becomes a relative error of about $40 \cdot 2^{-64}$ in the
  result, which is about $2^{58}$. At 64 bits the result is therefore
  only known to within $\pm 5.4$.
- `arb_printn(x, 35, 0)` prints at most 35 significant digits, but only
  digits that are certain (up to one unit in the last place), followed by
  an error bound such that the printed decimal ball contains the binary
  one. At 32 bits only 7 digits are known.
- `flint_cleanup` frees caches, here the cached value of $\pi$. It is
  optional.

Reading the output: at 128 bits the ball is
$[262537412640768743.99999999999925007 \pm 2.65 \cdot 10^{-18}]$, which
proves that $e^{\pi \sqrt{163}}$ is not an integer and that it lies
within $7.5 \cdot 10^{-13}$ of one. `arb_get_unique_fmpz` returns an
integer only if the ball contains exactly one, which is the standard way
to recover an exact integer that a computation is known to produce; here
it correctly fails.

## Getting numbers out {#output}

- To a string: `arb_get_str`, or `arb_printn` to print. See
  [](arb/printing).
- To a `double`: `arf_get_d(arb_midref(x), ARF_RND_NEAR)` returns the
  midpoint, and `mag_get_d(arb_radref(x))` an upper bound for the radius.
  `arb_can_round_arf` tells whether the midpoint rounds to the correctly
  rounded value of the number.
- To an integer: `arb_get_unique_fmpz` if the ball contains exactly one
  integer, or `arf_get_fmpz` of the midpoint for the nearest integer.
- To bounds: `arb_get_lbound_arf` and `arb_get_ubound_arf` give the
  endpoints, `arb_get_mag` an upper bound for the absolute value.

To get numbers in, use the exact conversions `arb_set_ui`, `arb_set_si`,
`arb_set_fmpz` and `arb_set_d`, or the rounding conversions
`arb_set_fmpq` and `arb_set_str`. Remember that a `double` constant such
as `0.1` is not the decimal number $1/10$.

## Where to go next {#next}

- [](arb/arithmetic): what the arithmetic operations compute exactly,
  and how precision and exactness interact.
- [](arb/predicates): comparisons in ball arithmetic.
- [](arb/elementary): the elementary functions, their domains and their
  accuracy.
- [](guide/types/ball-arithmetic) and [](guide/types/ball-arithmetic#precision):
  the guess-and-verify strategy and how to choose the precision.
- [acb](acb) for complex numbers, [arb_poly](arb_poly) and
  [arb_mat](arb_mat) for polynomials and matrices, and
  [arb_hypgeom](arb_hypgeom) for special functions.
- The tutorials in [](tutorials/numerics) for complete programs.
