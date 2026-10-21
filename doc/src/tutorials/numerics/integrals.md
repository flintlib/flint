---
title: Rigorous numerical integration
summary: Integrals with guaranteed error bounds using acb_calc_integrate: writing the integrand, holomorphy and branch cuts, piecewise functions, infinite intervals by truncation with a tail bound, and accuracy goals; the example programs integrals, taylor_integrals and integrals_double_exp.
order: 8
---

Numerical integration is usually done with error *estimates*: a
quadrature rule is applied with two different numbers of nodes, and the
difference is taken as the error. Such estimates can be fooled by
oscillating or nearly singular integrands. `acb_calc_integrate` in
[](acb_calc) computes integrals with *rigorous* error bounds, and it is
fast for the typical case of a piecewise analytic integrand: the number
of function evaluations grows only linearly with the number of digits.
This tutorial shows how to write integrands for it and how to handle
non-analytic points and infinite intervals.

## How it works {#method}

`acb_calc_integrate(res, f, param, a, b, goal, tol, options, prec)`
computes $\int_a^b f(x)\,dx$ along the straight segment from $a$ to $b$
in the complex plane [@Joh2018a]. On a segment, it applies
Gauss–Legendre quadrature with $n$ nodes, whose error is bounded using
the following classical estimate: if $f$ is analytic and bounded by $M$
inside the ellipse $E_\rho$ with foci at the endpoints and sum of
semi-axes $\rho$ times half the length of the segment, the error is at
most about $M \rho^{-2n}$ (up to a factor depending on $\rho$ and the
length). The algorithm evaluates $f$ on a ball containing such an
ellipse to get $M$, and chooses $\rho$ and $n$ to meet the accuracy
goal. If no admissible $n$ exists, for example because $f$ is not
analytic near the segment, the segment is bisected. On small segments
where nothing else works, a direct enclosure $(b - a) f([a, b])$ is
used. The result is the sum of all segment contributions with their
error bounds, so it is rigorous whatever subdivision is chosen; the
subdivision only affects the cost. The nodes and weights are computed
rigorously and cached [@JM2018].

## Writing the integrand {#integrand}

The integrand is a callback of type `acb_calc_func_t`:
```
int f(acb_ptr out, const acb_t z, void * param, slong order, slong prec)
```
`acb_calc_integrate` calls it with `order` 0 or 1:

- with `order` = 0, it must write $f(z)$, enclosing the values of $f$ on
  the ball $z$; $f$ may be any function, even discontinuous;
- with `order` = 1, it must also check that $f$ is *holomorphic* on the
  ball $z$, and write a non-finite value (for instance with
  `acb_indeterminate`) if it is not.

For functions built from arithmetic operations and meromorphic functions
(exp, sin, gamma, zeta, ...), nothing special is needed: at a pole the
ball arithmetic produces an infinite result automatically. Functions
with branch cuts or that are not analytic must perform the check. FLINT
provides versions of the usual functions with a flag `analytic` for this
purpose: `acb_sqrt_analytic`, `acb_rsqrt_analytic`, `acb_log_analytic`,
`acb_pow_analytic`, and for piecewise definitions on the real line
`acb_real_abs`, `acb_real_sgn`, `acb_real_heaviside`, `acb_real_floor`,
`acb_real_ceil`, `acb_real_max`, `acb_real_min` and `acb_real_sqrtpos`.
With `analytic` nonzero, they return an indeterminate result when the
input ball touches a point where the function is not holomorphic (a
branch cut, a sign change, a jump), which makes the integration
algorithm bisect there. The callback passes `order != 0` as the flag.

## Three integrals {#program}

The program computes three integrals at 128 bits:
$$
\int_0^8 \sin(x + e^x)\,dx, \qquad
\int_0^1 |x^4 + 10x^3 + 19x^2 - 6x - 6|\, e^x\,dx, \qquad
\int_0^\infty e^{-x^2}\,dx.
$$
The first, due to Rump, oscillates faster and faster: near $x = 8$ the
period is about $2\pi e^{-8} \approx 0.002$. Many quadrature programs
give wrong results for it without warning. The second, from a question
by Helfgott, has an absolute value whose argument changes sign inside
$[0, 1]$, so the integrand is not analytic there; `acb_real_abs` with
the `analytic` flag lets the algorithm find the break point by
bisection. The third has an infinite interval.

```c example name=numerics_integrals
#include <flint/acb.h>
#include <flint/acb_calc.h>

/* sin(x + exp(x)) */
int f_rump(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();          /* derivatives are not needed */
    acb_exp(res, z, prec);
    acb_add(res, res, z, prec);
    acb_sin(res, res, prec);
    return 0;
}

/* |x^4 + 10x^3 + 19x^2 - 6x - 6| exp(x) */
int f_helfgott(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;
    if (order > 1)
        flint_abort();
    acb_init(t);
    acb_add_si(res, z, 10, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, 19, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, -6, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, -6, prec);
    acb_real_abs(res, res, order != 0, prec);   /* not holomorphic at sign changes */
    acb_exp(t, z, prec);
    acb_mul(res, res, t, prec);
    acb_clear(t);
    return 0;
}

/* exp(-x^2) */
int f_gauss(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();
    acb_mul(res, z, z, prec);
    acb_neg(res, res);
    acb_exp(res, res, prec);
    return 0;
}

int main(void)
{
    acb_t a, b, s, t;
    mag_t tol;
    slong prec = 128;

    acb_init(a); acb_init(b); acb_init(s); acb_init(t);
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -prec);

    acb_set_si(a, 0);
    acb_set_si(b, 8);
    acb_calc_integrate(s, f_rump, NULL, a, b, prec, tol, NULL, prec);
    flint_printf("int_0^8 sin(x + exp(x)) dx = "); acb_printn(s, 30, 0); flint_printf("\n");

    acb_set_si(a, 0);
    acb_set_si(b, 1);
    acb_calc_integrate(s, f_helfgott, NULL, a, b, prec, tol, NULL, prec);
    flint_printf("int_0^1 |...| exp(x) dx    = "); acb_printn(s, 30, 0); flint_printf("\n");

    /* int_0^inf exp(-x^2) dx = int_0^N + tail, tail <= exp(-N^2) */
    acb_set_si(a, 0);
    acb_set_si(b, 10);              /* exp(-100) < 2^-144 */
    acb_calc_integrate(s, f_gauss, NULL, a, b, prec, tol, NULL, prec);
    acb_mul(t, b, b, prec);
    acb_neg(t, t);
    acb_exp(t, t, prec);
    arb_add_error(acb_realref(s), acb_realref(t));
    flint_printf("int_0^inf exp(-x^2) dx     = "); acb_printn(s, 30, 0); flint_printf("\n");

    acb_const_pi(t, prec);
    acb_sqrt(t, t, prec);
    acb_mul_2exp_si(t, t, -1);
    flint_printf("sqrt(pi)/2                 = "); acb_printn(t, 30, 0); flint_printf("\n");
    flint_printf("overlap: %d\n", acb_overlaps(s, t));

    acb_clear(a); acb_clear(b); acb_clear(s); acb_clear(t);
    mag_clear(tol);
    flint_cleanup();
    return 0;
}
```

```output
int_0^8 sin(x + exp(x)) dx = [0.347400172657247807879512159120 +/- 1.07e-31]
int_0^1 |...| exp(x) dx    = [11.1473105500571397339159020843 +/- 4.47e-29]
int_0^inf exp(-x^2) dx     = [0.886226925452758013649083741671 +/- 4.28e-31]
sqrt(pi)/2                 = [0.886226925452758013649083741671 +/- 4.28e-31]
overlap: 1
```

The results are real balls (the imaginary parts are exactly zero and are
not printed) because the path and the integrands are real.

## Accuracy goals and precision {#goals}

The arguments `goal` and `tol` set the target: the algorithm tries to
reach a relative error of $2^{-\text{goal}}$ or an absolute error of
`tol`, whichever is larger, on each segment. The usual choice, used
above, is `goal` = `prec` and `tol` = $2^{-\text{prec}}$. The working
precision `prec` is used for all evaluations of $f$; the algorithm does
not increase it, so if the integrand loses accuracy (for example through
cancellation, or because $|f|$ is huge compared with the integral), the
caller must raise `prec` above the goal. For integrals that are very
small or very large, `tol` should be set to about the magnitude of the
integral times $2^{-\text{goal}}$ (or to 0, for a purely relative goal),
since an absolute tolerance of $2^{-\text{prec}}$ is meaningless for an
integral of size $10^{-400}$. The goals are not guarantees: the returned
ball is always correct, but it can be wider than requested, for example
when a limit on the number of evaluations is reached. The return value
is `ARB_CALC_SUCCESS` if all segments converged and
`ARB_CALC_NO_CONVERGENCE` otherwise.

The options (an `acb_calc_integrate_opt_t`, initialised with
`acb_calc_integrate_opt_init`; passing `NULL` uses the defaults) limit
the quadrature degree (`deg_limit`), the number of function evaluations
(`eval_limit`) and the number of queued segments (`depth_limit`), select
a heap instead of a stack for the queue of segments (`use_heap`, which
processes the segments with the largest error first and helps for
integrands with singularities), and print diagnostic information
(`verbose`).

## Infinite intervals and singularities {#improper}

`acb_calc_integrate` requires a finite path on which $f$ is bounded.
Improper integrals are reduced to proper ones by hand, as in the program:
integrate over $[0, N]$ and add a rigorous bound for the rest, here
$$
\int_N^\infty e^{-x^2}\,dx \le \int_N^\infty \frac{x}{N} e^{-x^2}\,dx
= \frac{e^{-N^2}}{2N} \le e^{-N^2}.
$$
`arb_add_error` adds this bound to the radius. As with the truncation of
series ([](tutorials/numerics/series)), the bound must be valid and
should be smaller than the target accuracy; here $N = 10$ gives
$e^{-100} < 2^{-144}$. For integrable singularities at an endpoint, the
same technique applies (integrate from $\varepsilon$ and bound
$\int_0^\varepsilon$), or a change of variables removes the singularity.
Integrands with a singularity inside the interval are handled
automatically if they are bounded (the segments around a discontinuity
get small, and the direct enclosure takes over), but unbounded
integrands require splitting the interval.

## Complex paths {#complex}

The endpoints $a$ and $b$ are complex numbers, and $f$ can be any complex
function, so contour integrals are computed by summing integrals over
segments. This gives rigorous access to quantities defined by Cauchy's
integral formula: Taylor and Laurent coefficients, residues, and, by the
argument principle, the number of zeros of an analytic function in a
region. The example program below counts the zeros of $\zeta(s)$ in a
region this way.

## The example programs {#programs}

`examples/integrals.c` computes any of 39 test integrals, selected with
`-i n` (or `-i all`):

```
integrals -i n [-prec p] [-goal p] [-tol eps] [-twice] [-heap] [-verbose]
          [-verbose2] [-deg n] [-eval n] [-depth n] [-threads n]
```

`-prec` sets the working precision (default 64 bits), `-goal` and `-tol`
the relative and absolute goals (default: the precision and
$2^{-\text{prec}}$), and the remaining options correspond to the options
structure above; `-twice` runs the computation twice to show the cost
of computing quadrature nodes the first time. The integrals include
the three above (`-i 4`, `-i 6`, `-i 14`), $\int_0^1 \sqrt{1 - x^2}\,dx$
(`-i 3`, square root branch point at the endpoint), $\int_1^{101}
\lfloor x \rfloor\,dx$ (`-i 5`, many discontinuities),
$\int_0^1 \sin(1/x)\,dx$ and $\int_0^1 x \sin(1/x)\,dx$ (`-i 8`, `-i 9`,
infinitely many oscillations; use `-heap`), a residue of $\zeta(s)$ by a
closed contour around $s = 1$ (`-i 7`), the 10th Laurent coefficient of
the Weierstrass $\wp$-function (`-i 21`), the zero counting function
$N(1000)$ of $\zeta(s)$ by the argument principle (`-i 22`), integrals
of $\Gamma$, Lambert $W$, Airy, error, Bessel and sine integral
functions, and several infinite intervals with explicit tail bounds:

```
> build/examples/integrals -i 3 -prec 333 -tol 1e-80
I3 = 4 int_0^1 sqrt(1-x^2) dx ...
cpu/wall(s): 0.013 0.014
I3 = [3.141592653589793238462643383279502884197169399375105820974944592307816406286209 +/- 4.66e-79]

> build/examples/integrals -i 9 -heap
I9 = int_0^1 x sin(1/x) dx  (slow convergence, use -heap and/or -tol) ...
cpu/wall(s): 0.021 0.022
I9 = [0.378530017 +/- 6.12e-10]

> build/examples/integrals -i 22
I22 = N(1000) = count zeros with 0 < t <= 1000 of zeta(s) using argument principle ...
cpu/wall(s): 2.412 2.515
I22 = [649.000000000000000 +/- 3.12e-16]
```

Since $N(1000)$ is an integer, the last result proves that $\zeta(s)$ has
exactly 649 zeros with $0 < \operatorname{Im}(s) \le 1000$ in the
critical strip.

`examples/taylor_integrals.c` uses the older function
`acb_calc_integrate_taylor`, which integrates by Taylor expansion and
needs the integrand's Taylor coefficients (`order` > 1): it computes
$\int_0^{100} \sin t\,dt$, the elliptic integral
$F(\phi, m) = \int_0^\phi (1 - m \sin^2 t)^{-1/2}\,dt$ with
$\phi = 6 + 6i$, $m = 1/2$ along the path $0 \to 6 \to 6 + 6i$, and the
Bessel function $J_{10}(20 + 10i)$ from its integral representation,
to $d$ digits (`taylor_integrals d`). This method is usually slower than
`acb_calc_integrate`, which supersedes it.

`examples/integrals_double_exp.c`
(`integrals_double_exp [-p digits] [-t threads] [-v]`, by default at 1024
bits) demonstrates the double exponential (tanh-sinh) quadrature on
integrands
with logarithmic and square-root endpoint singularities and on infinite
intervals, where it converges quickly without any special treatment of
the singularities:

```
> build/examples/integrals_double_exp
example: integral_0^oo 1/(1+x^3) dx = 2*pi/(3*sqrt(3))

	result = 1.209199576156145 +/- 4.6453e-305, |result - exact| = 0e-304
...
```

::: warning
The radius printed by `integrals_double_exp` only accounts for the
rounding errors in the evaluation of the quadrature sum. The
discretisation and truncation errors of the double exponential rule are
estimated heuristically (by comparing successive step sizes and
stopping when terms no longer change the sum) and are not included in
the ball, so this program, unlike `acb_calc_integrate`, does not produce
rigorous enclosures. It illustrates the general rule that ball
arithmetic bounds only the errors of the operations that are actually
performed.
:::
