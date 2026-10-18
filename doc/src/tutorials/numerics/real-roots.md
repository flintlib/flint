---
title: Isolating the real roots of a function
summary: Certified isolation of all roots of a real function on an interval by bisection and sign tests, refinement by bisection and by Newton's method with a rigorous error bound, with arb_calc; the example program real_roots.
order: 5
---

Finding the roots of a function numerically is easy; proving that the
computed list is complete, and that each reported root really is a root
known to a stated accuracy, is not. The module [](arb_calc) does both for
a real function $f$ given by a program that evaluates $f$ and a few
derivatives in ball arithmetic. This tutorial uses it to find all zeros
of the Bessel function $J_0$ on $[0, 20]$ to 30 digits, and explains the
mathematics that makes the result rigorous.

## Describing the function {#callback}

The function is passed as a callback of type `arb_calc_func_t`:
```
int f(arb_ptr out, const arb_t x, void * param, slong order, slong prec)
```
It must write to `out` the first `order` Taylor coefficients of $f$ at
the ball $x$, that is $f(x), f'(x), f''(x)/2, \ldots$, each a ball that
contains the corresponding value for every point of $x$. The root
isolation needs `order` up to 2 and the Newton refinement up to 3. The
`param` pointer passes extra data through; below it points to a counter
of evaluations. For $J_0$ we use $J_0' = -J_1$ and the differential
equation $x^2 y'' + x y' + x^2 y = 0$, which gives
$J_0''(x)/2 = (J_1(x)/x - J_0(x))/2$. When the ball $x$ contains 0,
the division gives an infinite ball; that is correct (the formula is not
valid there) and only means that this ball cannot be used to prove
anything.

For a function built from elementary and special functions, the Taylor
coefficients can be computed mechanically with power series arithmetic:
`examples/real_roots.c` evaluates $\sin(x^2)$ by composing
`_arb_poly_mullow` and `_arb_poly_sin_series` on the series $x + t$, and
the Riemann–Siegel $Z$-function with `_arb_poly_riemann_siegel_z_series`.

## Isolation {#isolation}

`arb_calc_isolate_roots` subdivides the interval $[a, b]$ recursively. On
each subinterval $I$ it evaluates $f$ on the ball containing $I$:

- if the result does not contain zero, $I$ contains no root and is
  discarded;
- if $f$ has opposite signs at the endpoints of $I$ and the ball for
  $f'(I)$ does not contain zero, then $I$ contains exactly one root and is
  reported as isolated;
- otherwise $I$ is bisected, up to a maximum depth and a maximum number
  of evaluations.

::: proposition {#prop-isolated}
If $f$ is continuously differentiable on $I = [u, v]$,
$f(u) f(v) < 0$ and $f'$ has no zero on $I$, then $f$ has exactly one
root in $I$, and it is simple.
:::

::: proof
By the intermediate value theorem there is a root. Since $f'$ is
continuous and nonzero on $I$, it has constant sign, so $f$ is strictly
monotonic and the root is unique; it is simple because $f'$ does not
vanish there.
:::

The output is a list of intervals with flags: flag 1 means "exactly one
root", flag 2 means "could not be decided" (a root of even multiplicity,
a cluster of roots, a root at an endpoint, or simply a limit on the
depth or the number of evaluations reached). The function has no roots
in $[a, b]$ outside the reported intervals, so if all flags are 1, the
list is complete. The interval endpoints are `arf_t` floating-point numbers, so
the isolation itself is exact; the working precision only needs to be
large enough to decide the signs, and 30 bits are usually plenty.

## Refinement {#refinement}

An isolating interval is refined in two stages. A few bisection steps
with `arb_calc_refine_root_bisect` (keeping the half where the sign
changes) give a smaller interval at low precision. Then Newton's method
gives quadratic convergence, with an error bound based on the following
estimate.

::: proposition Error of a Newton step {#prop-newton}
Let $f$ be twice continuously differentiable on an interval $X$ containing
a simple root $\xi$ and a point $m$, with $f' \ne 0$ on $X$, and let
$$
C = \frac{\max_{X} |f''|}{2 \min_{X} |f'|}.
$$
Then $\left| m - f(m)/f'(m) - \xi \right| \le C |m - \xi|^2$.
:::

::: proof
By Taylor's theorem, $0 = f(\xi) = f(m) + f'(m)(\xi - m) + \tfrac12
f''(\eta)(\xi - m)^2$ for some $\eta$ between $m$ and $\xi$. Dividing by
$f'(m)$ gives $m - f(m)/f'(m) - \xi = \tfrac12 f''(\eta) (\xi - m)^2 /
f'(m)$.
:::

`arb_calc_newton_conv_factor` computes an upper bound for $C$ over a
region $X$ from one evaluation with `order` 3 on the ball $X$. Given a
ball $x = [m \pm r]$ containing the root, `arb_calc_newton_step`
evaluates $f(m)$ and $f'(m)$ at the exact point $m$, computes
$m - f(m)/f'(m)$ in ball arithmetic and adds $C r^2$ to the radius; the
result contains $\xi$ by the proposition. The step is accepted only if
the new ball is contained in $X$ and is smaller than $x$.
`arb_calc_refine_root_newton` repeats this with a precision that
roughly doubles in each step, starting from the accuracy of the initial
ball, so that the final steps, which dominate the cost, are done at the
target precision only once.

## The program {#program}

```c example name=numerics_real_roots_j0
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>
#include <flint/arb_calc.h>

/* Taylor coefficients of J_0 at x: J_0, J_0' = -J_1, J_0''/2 = (J_1/x - J_0)/2 */
int bessel_j0(arb_ptr out, const arb_t x, void * param, slong order, slong prec)
{
    arb_t nu, j1;
    arb_init(nu);
    arb_init(j1);

    arb_hypgeom_bessel_j(out, nu, x, prec);     /* nu = 0 */
    if (order >= 2)
    {
        arb_one(nu);
        arb_hypgeom_bessel_j(j1, nu, x, prec);
        arb_neg(out + 1, j1);
    }
    if (order >= 3)
    {
        arb_div(out + 2, j1, x, prec);
        arb_sub(out + 2, out + 2, out, prec);
        arb_mul_2exp_si(out + 2, out + 2, -1);
    }

    (*(slong *) param)++;      /* count evaluations */
    arb_clear(nu);
    arb_clear(j1);
    return 0;
}

int main(void)
{
    arf_interval_t interval, block;
    arf_interval_ptr blocks;
    int * flags;
    slong i, num, evals = 0;
    slong prec = 30, digits = 30, high_prec = digits * 3.33 + 10;
    arb_t v, w, root;
    arf_t C;

    arf_interval_init(interval);
    arf_interval_init(block);
    arb_init(v); arb_init(w); arb_init(root);
    arf_init(C);

    arf_set_si(&interval->a, 0);
    arf_set_si(&interval->b, 20);

    /* step 1: isolation */
    num = arb_calc_isolate_roots(&blocks, &flags, bessel_j0, &evals,
        interval, 30, 100000, 100000, prec);

    for (i = 0; i < num; i++)
    {
        if (flags[i] != 1)
        {
            flint_printf("could not decide [%.6f, %.6f]\n",
                arf_get_d(&blocks[i].a, ARF_RND_DOWN), arf_get_d(&blocks[i].b, ARF_RND_UP));
            continue;
        }

        /* step 2: a few bisection steps, at low precision */
        arb_calc_refine_root_bisect(block, bessel_j0, &evals, blocks + i, 10, prec);

        /* step 3: Newton iteration with error bounds, at high precision */
        arf_interval_get_arb(v, blocks + i, high_prec);       /* region of convergence */
        arb_calc_newton_conv_factor(C, bessel_j0, &evals, v, prec);
        arf_interval_get_arb(w, block, high_prec);            /* starting ball */
        if (arb_calc_refine_root_newton(root, bessel_j0, &evals, w, v, C, 10, high_prec)
                != ARB_CALC_SUCCESS)
            flint_printf("Newton iteration failed\n");

        arb_printn(root, digits, 0);
        flint_printf("\n");
    }
    flint_printf("%wd blocks, %wd function evaluations\n", num, evals);

    for (i = 0; i < num; i++)
        arf_interval_clear(blocks + i);
    flint_free(blocks);
    flint_free(flags);
    arf_interval_clear(interval);
    arf_interval_clear(block);
    arb_clear(v); arb_clear(w); arb_clear(root);
    arf_clear(C);
    flint_cleanup();
    return 0;
}
```

```output
[2.40482555769577276862163187933 +/- 3.55e-30]
[5.52007811028631064959660411281 +/- 3.03e-30]
[8.65372791291101221695419871266 +/- 9.47e-31]
[11.7915344390142816137430449119 +/- 2.55e-29]
[14.9309177084877859477625939974 +/- 1.14e-29]
[18.0710639679109225431478829756 +/- 1.82e-29]
6 blocks, 1776 function evaluations
```

All six intervals have flag 1, so $J_0$ has exactly six zeros in
$[0, 20]$, and each is known to 30 digits. Notes:

- `arb_calc_isolate_roots` allocates the arrays `blocks` and `flags`;
  the caller clears the intervals and frees both arrays with
  `flint_free`.
- The arguments `30, 100000, 100000` are the maximum recursion depth,
  the approximate maximum number of evaluations, and the maximum number
  of isolated roots to find before stopping.
- The extra precision `10` passed to `arb_calc_refine_root_newton` is
  added to the working precision for evaluating $f$, to compensate for
  rounding errors in the evaluation.
- Most of the evaluations are spent in the isolation and bisection at 30
  bits; the Newton steps are few but at higher precision.

For roots of polynomials with exact coefficients, specialised and much
faster methods exist; see [](tutorials/numerics/polynomial-roots). For
an alternative certification of a single root, the interval Newton
method, see [](tutorials/numerics/interval-newton).

## The example program {#programs}

`examples/real_roots.c` isolates the roots of one of eight built-in
functions on an interval $(a, b)$ given as `double` values:

```
real_roots function a b [-refine d] [-verbose] [-maxdepth n] [-maxeval n] [-maxfound n] [-prec n]
```

The functions are 0: the Hardy $Z$-function $Z(x)$ (with
`-character q n`, of the Dirichlet $L$-function of the character with
Conrey label $(q, n)$ instead of the Riemann zeta function); 1:
$\sin(x)$; 2: $\sin(x^2)$; 3: $\sin(1/x)$; 4 to 7: the Airy functions
$\operatorname{Ai}(x)$, $\operatorname{Ai}'(x)$, $\operatorname{Bi}(x)$,
$\operatorname{Bi}'(x)$. With `-refine d`, each isolated root is
refined to about $d$ digits as in the program above (two rounds of 5
bisection steps, then Newton) and printed; `-maxdepth`, `-maxeval` and
`-maxfound` set the limits of the isolation (defaults 30, 100000 and
100000) and `-prec` its working precision (default 30). The program
reports the number of roots found, the number of undecided intervals,
and the number of function evaluations.

The zeros of $Z(t)$ are the imaginary parts of the zeros of $\zeta(s)$
on the critical line, so this finds and certifies them, missing none:

```
> build/examples/real_roots 0 0 50 -verbose
interval: [0, 50.0000000000000]
maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
found isolated root in: [14.0625000000000, 14.1601562500000]
found isolated root in: [20.9960937500000, 21.0937500000000]
found isolated root in: [25.0000000000000, 25.0488281250000]
found isolated root in: [30.4199218750000, 30.4687500000000]
found isolated root in: [32.9101562500000, 32.9589843750000]
found isolated root in: [37.5488281250000, 37.5976562500000]
found isolated root in: [40.9179687500000, 40.9667968750000]
found isolated root in: [43.3105468750000, 43.3593750000000]
found isolated root in: [47.9980468750000, 48.0468750000000]
found isolated root in: [49.7558593750000, 49.8046875000000]
---------------------------------------------------------------
Found roots: 10
Subintervals possibly containing undetected roots: 0
Function evaluations: 1783
```

(`examples/zeta_zeros.c` computes these zeros far more efficiently; see
[](tutorials/numerics/zeta).) The zeros of the Airy function, refined to
30 digits:

```
> build/examples/real_roots 4 -10 0 -refine 30
interval: [-10.0000000000000, 0]
maxdepth = 30, maxeval = 100000, maxfound = 100000, low_prec = 30
refined root (0/6):
[-9.0226508533409803801581908398801 +/- 1.08e-32]

refined root (1/6):
[-7.9441335871208531231382805557983 +/- 3.15e-32]

...

refined root (5/6):
[-2.3381074104597670384891972524467 +/- 3.55e-32]

---------------------------------------------------------------
Found roots: 6
Subintervals possibly containing undetected roots: 0
Function evaluations: 184
```

The method fails gracefully where it must. On $(0, 100)$, $\sin(x^2)$
has a double root at the endpoint 0, which cannot be isolated, and the
program reports one undecided interval (`real_roots 2 0 100` finds 3183
roots and 1 undecided interval; `real_roots 2 1 100` finds the same 3183
roots and none undecided). The function $\sin(1/x)$ has infinitely many
roots accumulating at 0; `real_roots 3 0 1` stops at the evaluation
limit with about 14000 isolated roots and about 22000
undecided intervals near 0. In all cases the isolated roots are correct and every
root lies in a reported interval.
