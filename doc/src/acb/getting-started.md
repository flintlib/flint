---
title: Getting started with acb
summary: Complex numbers as rectangles of real balls - a first program, working with the real and imaginary parts, printing, and how exactness and the precision behave.
order: 1
---

An `acb_t` is a complex number $a + bi$ whose real and imaginary parts
are real balls (`arb_t`). Everything said about precision, predicates and
printing in [](arb/getting-started) applies to each part.

## A first program {#first}

```c example name=acb_getting_started
#include <flint/flint.h>
#include <flint/acb.h>

int main(void)
{
    acb_t z, w;
    arb_t r;
    slong prec = 64;

    acb_init(z); acb_init(w); arb_init(r);

    acb_set_si_si(z, 1, 2);                   /* z = 1 + 2i */
    flint_printf("z        = "); acb_printn(z, 15, 0); flint_printf("\n");

    acb_exp(w, z, prec);
    flint_printf("exp(z)   = "); acb_printn(w, 15, 0); flint_printf("\n");
    acb_log(w, z, prec);
    flint_printf("log(z)   = "); acb_printn(w, 15, 0); flint_printf("\n");
    acb_sqrt(w, z, prec);
    flint_printf("sqrt(z)  = "); acb_printn(w, 15, 0); flint_printf("\n");
    acb_gamma(w, z, prec);
    flint_printf("gamma(z) = "); acb_printn(w, 15, 0); flint_printf("\n");

    acb_abs(r, z, prec);
    flint_printf("|z|      = "); arb_printn(r, 15, 0); flint_printf("\n");
    acb_arg(r, z, prec);
    flint_printf("arg(z)   = "); arb_printn(r, 15, 0); flint_printf("\n");

    /* exp(pi i) + 1 = 0, up to the error of the enclosure of pi */
    acb_const_pi(w, prec);
    acb_mul_onei(w, w);
    acb_exp(w, w, prec);
    acb_add_ui(w, w, 1, prec);
    flint_printf("exp(pi i) + 1 = "); acb_printn(w, 15, 0);
    flint_printf(", contains 0: %d\n", acb_contains_zero(w));

    /* the same with an exact function of pi z */
    acb_one(z);
    acb_exp_pi_i(w, z, prec);
    acb_add_ui(w, w, 1, prec);
    flint_printf("exp_pi_i(1) + 1 = "); acb_printn(w, 15, 0); flint_printf("\n");

    acb_clear(z); acb_clear(w); arb_clear(r);
    flint_cleanup();
    return 0;
}
```
```output
z        = 1.00000000000000 + 2.00000000000000*I
exp(z)   = [-1.13120438375681 +/- 3.64e-15] + [2.47172667200482 +/- 1.08e-15]*I
log(z)   = [0.804718956217050 +/- 1.88e-16] + [1.10714871779409 +/- 5.04e-16]*I
sqrt(z)  = [1.27201964951407 +/- 1.04e-15] + [0.786151377757423 +/- 2.87e-16]*I
gamma(z) = [0.151904002670036 +/- 1.38e-16] + [0.0198048801618550 +/- 1.81e-17]*I
|z|      = [2.23606797749979 +/- 3.05e-16]
arg(z)   = [1.10714871779409 +/- 5.04e-16]
exp(pi i) + 1 = [+/- 1.10e-19] + [+/- 1.59e-19]*I, contains 0: 1
exp_pi_i(1) + 1 = 0
```

## Real and imaginary parts {#parts}

`acb_realref(z)` and `acb_imagref(z)` give the two parts as `arb_t`
pointers that can be passed to any `arb` function, for reading or for
writing. Many operations that the `acb` module does not provide are easily
written this way. The parts are independent balls: the rectangle may have
an exact real part and an inexact imaginary part, and a
real number stored in an `acb_t` has an imaginary part that is exactly
zero. Functions use this: `acb_mul` of a real and a complex number does two
real multiplications, and `acb_is_real` tests whether the imaginary part
is exactly zero. Conversely, an imaginary part such as $[\pm 10^{-30}]$
coming from rounding errors is not exactly zero, so a computation that is
known to give a real number should end with `arb_zero(acb_imagref(z))` or
use `acb_get_real` if the imaginary part is not needed.

The functions that return real values, such as `acb_abs`, `acb_arg` and
`acb_csgn`, take an `arb_t` output.

## Exactness and the precision {#precision}

As for real balls, the precision bounds the length of the midpoints of
the result, exact results are exact (`1 + 2i` above is printed without a
radius), and the radii are rigorous bounds, not estimates. The last two
lines of the program illustrate the difference between computing with an
enclosure of $\pi$, which gives a small rectangle around 0, and using a
function of $\pi z$ such as `acb_exp_pi_i`, `acb_sin_pi` or
`acb_cos_pi`, which can return exact results at rational points.

## Branch cuts {#branch-cuts}

The logarithm, the square root, powers and the inverse trigonometric and
hyperbolic functions are multivalued, and `acb` returns their principal
values with the standard branch cuts. A rectangle that crosses a branch
cut gives an enclosure of the values on both sides, which can be large;
the guide [](acb/branch-cuts) explains the conventions, the `analytic`
variants and the functions that avoid the problem.

## Where to go next {#next}

- [](acb/arithmetic) for the error bounds of complex arithmetic.
- [acb_hypgeom](acb_hypgeom), [acb_elliptic](acb_elliptic),
  [acb_modular](acb_modular) and [acb_dirichlet](acb_dirichlet) for
  special functions.
- [acb_poly](acb_poly), [acb_mat](acb_mat) and [acb_calc](acb_calc) for
  polynomials, matrices and rigorous integration.
