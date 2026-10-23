---
title: Special functions in ball arithmetic
summary: Calling special functions on balls, how input uncertainty and cancellation show up in the output, when to raise the precision, and accurate double-precision values with arb_fpwrap; the example programs fpwrap, complex_plot and functions_benchmark.
order: 4
---

FLINT evaluates a large collection of special functions in ball
arithmetic: elementary functions in [](arb) and [](acb), the gamma
function and its relatives, the error function, Bessel, Airy and
hypergeometric functions in [](arb_hypgeom) and [](acb_hypgeom), elliptic
and modular functions in [](acb_elliptic) and [](acb_modular), theta
functions in [](acb_theta), and zeta and $L$-functions in
[](acb_dirichlet). All of them follow the same contract: given balls as
input and a working precision, they return a ball that contains the
function value for every point of the input balls. This tutorial looks
at what the output tells about the input and the formula, and at the
floating-point interface [](arb_fpwrap) for programs that only want
accurate `double` values.

## Exact input, inexact input, and cancellation {#uncertainty}

The program below evaluates the gamma function at $1/3$, first rounded
to 64 bits (a ball of radius about $2^{-66}$) and then with a radius of
$10^{-10}$; then $\sin(10^{30})$ for an exact argument; and finally a
formula with catastrophic cancellation,
$$
f(x) = \frac{e^x - 1 - x}{x^2}, \qquad x = 2^{-40},
$$
at increasing precision until the result has 48 accurate bits.

```c example name=numerics_special_functions
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>

/* f(x) = (exp(x) - 1 - x) / x^2, evaluated naively */
void f_naive(arb_t res, const arb_t x, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_exp(t, x, prec);
    arb_sub_ui(t, t, 1, prec);
    arb_sub(t, t, x, prec);
    arb_sqr(res, x, prec);
    arb_div(res, t, res, prec);
    arb_clear(t);
}

int main(void)
{
    arb_t x, y;
    slong prec;

    arb_init(x);
    arb_init(y);

    /* exact and inexact input */
    arb_set_ui(x, 1);
    arb_div_ui(x, x, 3, 64);                /* 1/3, rounded to 64 bits */
    arb_gamma(y, x, 64);
    flint_printf("gamma(1/3)            = "); arb_printn(y, 20, 0); flint_printf("\n");
    mag_set_d(arb_radref(x), 1e-10);        /* x = [1/3 +/- 1e-10] */
    arb_gamma(y, x, 64);
    flint_printf("gamma(1/3 +/- 1e-10)  = "); arb_printn(y, 20, 0); flint_printf("\n");

    /* a huge exact argument */
    arb_set_ui(x, 10);
    arb_pow_ui(x, x, 30, 128);
    arb_sin(y, x, 64);
    flint_printf("sin(10^30)            = "); arb_printn(y, 15, 0); flint_printf("\n");

    /* catastrophic cancellation: increase the precision until accurate */
    arb_one(x);
    arb_mul_2exp_si(x, x, -40);             /* x = 2^-40, exact */
    for (prec = 32; ; prec *= 2)
    {
        f_naive(y, x, prec);
        flint_printf("prec = %3wd: f(2^-40) = ", prec);
        arb_printn(y, 15, 0);
        flint_printf("\n");
        if (arb_rel_accuracy_bits(y) >= 48)
            break;
    }

    arb_clear(x);
    arb_clear(y);
    flint_cleanup();
    return 0;
}
```

```output
gamma(1/3)            = [2.678938534707747634 +/- 9.48e-19]
gamma(1/3 +/- 1e-10)  = [2.67893853 +/- 5.58e-9]
sin(10^30)            = [-0.0901169019121381 +/- 4.20e-17]
prec =  32: f(2^-40) = [+/- 5.50e+12]
prec =  64: f(2^-40) = [+/- 4.01]
prec = 128: f(2^-40) = [0.50000000000015 +/- 7.91e-15]
prec = 256: f(2^-40) = [0.500000000000152 +/- 4.18e-16]
```

The three parts show three different sources of error:

- *Input uncertainty.* With the input $[1/3 \pm 10^{-10}]$, the radius
  of the output is $8.7 \cdot 10^{-10}$, close to
  $|\Gamma'(1/3)| \cdot 10^{-10} \approx 8.4 \cdot 10^{-10}$ (the
  printed radius is larger because it includes the rounding of the
  midpoint to 8 digits). No working precision can do better: the output describes
  the gamma function on the whole input ball. When the inputs are
  measured or rounded data, this is exactly the information wanted;
  when they are exact numbers, they should be represented exactly or
  recomputed at each precision (see
  [](tutorials/numerics/adaptive-precision#inexact-input)).
- *A large exact argument.* The argument $10^{30}$ is exact, and
  `arb_sin` reduces it modulo $2\pi$ using enough bits of $\pi$
  internally, so 64 bits of working precision give a result with 62
  accurate bits. In floating point, $\sin(10^{30})$ is also
  well-defined for the exactly representable argument, but an argument
  that has itself been rounded to 53 bits would carry an absolute error
  of about $10^{14}$, and its sine would be meaningless. In ball
  arithmetic that situation would show up as a radius larger than 1.
- *Cancellation.* $e^x - 1 - x \approx x^2 / 2 = 2^{-81}$ is computed as
  a difference of numbers of size 1, so about 81 bits cancel. At 32 and
  64 bits nothing is left, and the balls are honestly wide; at 128 bits
  about 47 bits remain. Raising the precision is a correct but expensive
  fix. The better fix is a formula without cancellation, here
  `arb_expm1` (which computes $e^x - 1$ accurately for small $x$)
  followed by one subtraction, or the Taylor series
  $\tfrac12 + \tfrac{x}{6} + \tfrac{x^2}{24} + \cdots$ with a tail bound.

The precision of library functions is handled internally: a function
such as `arb_gamma` or `arb_hypgeom_bessel_j` increases its internal
precision where its algorithm needs it (for example in argument
reduction or in series with cancellation), and aims to return a result
whose accuracy is close to the working precision when the input is
exact. It does not guarantee this in all cases, so a program that needs
a given accuracy should check `arb_rel_accuracy_bits` of the result and
retry at a higher precision if needed, as in the loop above. Near
singularities or for huge parameters, the result may be a wide or
infinite ball; it is still correct, just not informative. The details of
the algorithms are described in [](algorithms/special-functions/overview).

## Accurate floating-point values with arb_fpwrap {#fpwrap}

Programs that work in `double` precision can use the wrappers in
[](arb_fpwrap). Each wrapper converts its floating-point arguments to
exact balls, evaluates the function in ball arithmetic, and increases
the working precision until the result rounds to a `double` with the
requested accuracy: by default, a real result is accurate up to the
rounding of the last bit (which is not guaranteed to be rounded
optimally), and a complex result is accurate to 53 bits relative to its
absolute value; with `FPWRAP_CORRECT_ROUNDING`, the result is correctly
rounded, and with `FPWRAP_ACCURATE_PARTS`, both the real and the
imaginary part of a complex result have full relative accuracy. The return value is `FPWRAP_SUCCESS`
(0), or `FPWRAP_UNABLE` (1) if the accuracy could not be reached within
the work limit, for instance at a pole; the result is then NaN.

```c example name=numerics_fpwrap
#include <math.h>
#include <flint/arb_fpwrap.h>

int main(void)
{
    double x, y;
    complex_double z, w;
    int status;

    /* gamma close to a pole */
    x = -3.0 + 1e-8;
    status = arb_fpwrap_double_gamma(&y, x, 0);
    printf("gamma(%.17g) = %.16g (status %d)\n", x, y, status);

    status = arb_fpwrap_double_bessel_j(&y, 0.0, 1000.0, FPWRAP_CORRECT_ROUNDING);
    printf("J_0(1000) = %.17g (status %d)\n", y, status);

    z.real = 0.5;
    z.imag = 1000.0;
    status = arb_fpwrap_cdouble_zeta(&w, z, FPWRAP_ACCURATE_PARTS);
    printf("zeta(0.5 + 1000i) = %.16g + %.16gi (status %d)\n", w.real, w.imag, status);

    status = arb_fpwrap_double_zeta(&y, 1.0, 0);
    printf("zeta(1): isnan = %d (status %d)\n", isnan(y) != 0, status);

    flint_cleanup();
    return 0;
}
```

```output
gamma(-2.9999999900000001) = -16666666.9773108 (status 0)
J_0(1000) = 0.024786686152420176 (status 0)
zeta(0.5 + 1000i) = 0.356334367194396 + 0.9319978312329936i (status 0)
zeta(1): isnan = 1 (status 1)
```

The argument $-3 + 10^{-8}$ is not exactly representable; the wrapper
computes $\Gamma$ at the `double` nearest to it, printed in the first
line with 17 digits, and that value is what the result approximates.
The wrappers are convenient for testing other implementations and for
values that are hard to get right in floating point (near poles and
zeros, for large arguments, for functions with many parameters). They
are slower than a dedicated `double` implementation, since every call
goes through arbitrary-precision ball arithmetic and possibly several
attempts at increasing precision.

## The example programs {#programs}

`examples/fpwrap.c` is the minimal version of the program above:

```
> build/examples/fpwrap
zeta(2) = 1.644934066848226
zeta(0.5 + 123i) = 0.006252861175594465 + 0.08206030514520983i
```

`examples/complex_plot.c` plots a function of a complex variable on a
rectangle by domain colouring and writes the image to `arbplot.ppm`
(convert it with, for example, `convert arbplot.ppm arbplot.png`):

```
complex_plot [-range xa xb ya yb] [-size xn yn] [-color n] [-threads n] <func>
```

The default range is $[-10, 10] + [-10, 10] i$ at $512 \times 512$
pixels. The functions include `sin`, `gamma`, `digamma`, `lgamma`,
`zeta`, `erf`, the Airy functions `ai` and `bi`, the Bessel functions
`besselj`, `bessely`, `besseli`, `besselk` (of order 0), modular
functions (`modj`, `modeta`, `modlambda`, and the variants `modjq`,
`modetaq`, `modlambdaq` as functions of the nome $q$), the Weierstrass
functions `ellipp`, `ellipzeta`, `ellipsigma` for a square lattice,
`barnesg`, `agm`, and the Fresnel integrals `fresnels` and `fresnelc`.
`-color` selects one of seven colour schemes (0: phase as hue and
magnitude as brightness; 1: phase only; 2 to 6: combinations of
these). Each pixel is computed by evaluating the function at one point,
at a precision doubled from 30 bits until the value has a few accurate
bits (pixels where this fails below 500 bits are drawn grey), so the
colours are correct even where `double` arithmetic would fail, for
example for $\zeta(s)$ high in the critical strip; features smaller than
a pixel can be missed.
For example, the following plots the zeta function in a part of the
critical strip:

```
> build/examples/complex_plot zeta -range -10 10 100 140 -size 256 512
```

`examples/functions_benchmark.c` measures the time to compute $\pi$,
$\gamma$, $\exp$, $\log$, $\sin$, $\operatorname{atan}$,
$\operatorname{erf}$, $\Gamma$, $\zeta$ and the Dedekind eta function to
$n$ digits for $n$ from 10 up to $10^8$ (depending on the function), as
well as Bernoulli numbers and partition numbers, with 1 and 8 threads
and for a first and a repeated call (the repeated call profits from
cached constants and precomputed data), with the argument
$x = \sqrt{2} - 1$. The largest cases take a long time and much memory.
