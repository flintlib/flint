---
title: Isolating and refining real roots
functions: arb_calc_isolate_roots arb_calc_refine_root_bisect arb_calc_partition
---

`arb_calc_isolate_roots` finds the simple real roots of an analytic
function $f$ on an interval $[a, b]$ by recursive bisection. The
function is given by an evaluator of type `arb_calc_func_t` that
encloses $f$ and $f'$ on a ball. The result is a list of subintervals
with flags: the subintervals with flag 1 contain exactly one root each,
and no root of $f$ in $[a, b]$ lies outside the returned subintervals.
Isolated roots can then be refined to high precision by bisection
(`arb_calc_refine_root_bisect`, linear convergence, one bit per
evaluation) or by Newton's method (`arb_calc_refine_root_newton`,
quadratic convergence).

::: summary algorithms/ball-arithmetic/real-root-isolation :::

## The bisection algorithm {#algorithm}

The function first determines the signs of $f$ at the endpoints $a$ and
$b$ (0 if the value contains zero). Each subinterval $I = [u, v]$ with
known endpoint signs is then processed as follows:

1. If the enclosure of $f(I)$ excludes zero, $I$ contains no root and is
   discarded.
2. If $f(u)$ and $f(v)$ have opposite nonzero signs and the enclosure of
   $f'(I)$ excludes zero, then $f$ is strictly monotonic on $I$ and has
   exactly one root there, which is simple: $I$ is returned with flag 1.
3. Otherwise $I$ is split at its exact midpoint $m$ by
   `arb_calc_partition`, which also computes the sign of $f(m)$, and the
   two halves are processed recursively, left before right.

The recursion stops after `maxdepth` levels (returning the subinterval
with flag 2, meaning "unknown"), after `maxeval` subintervals have been
examined (returning the unexamined ones with flag 2), or when `maxfound`
roots have been isolated. Each subinterval costs one evaluation of $f$
on $I$ and, in case 2, one of $f$ and $f'$ on $I$, plus one evaluation
at the midpoint when it is split.

Correctness only depends on the enclosures computed by the evaluator,
which must contain $f(I)$ and $f'(I)$ for the *whole* interval $I$.
Overestimation (for example from the dependency problem when $f$ is
evaluated naively on a wide ball) only causes more bisections. A root of
multiplicity greater than 1 is never isolated, because $f'$ vanishes
there; nor is a root that happens to be
exactly a bisection point, since the sign of $f$ there is zero. Such
roots end up in small subintervals with flag 2.

## Writing the evaluator {#evaluator}

The evaluator receives a ball `inp`, an `order` (1 or 2 for root
isolation) and a precision, and writes the Taylor coefficients $f(x)$
and $f'(x)$ to `out`. The simplest way is to compute $f(x + t)$ as a
power series in $t$ with the functions of [](arb_poly), starting from
the series $x + t$ (see [](arb_poly/power-series)); for elementary
functions, `_arb_poly_sin_series` and its relatives accept the
coefficient array directly. If the precision passed is too low for the
function, the enclosures are wide and the search stalls with many
flag-2 intervals.

## Refinement {#refinement}

`arb_calc_refine_root_bisect(r, f, param, start, iter, prec)` halves an
interval with a sign change `iter` times. It needs the signs at the
endpoints of `start` to be determined and different, and stops early
if the sign at a midpoint cannot be determined.

`arb_calc_refine_root_newton` performs Newton steps $m \leftarrow m -
f(m)/f'(m)$ with the error bound $C r^2$, where $C$ bounds $\tfrac12
|f''|/|f'|$ on a region $I$ around the root, at precisions that double
up to the target. It converges quadratically once the starting ball is
small enough that $C r < 1$, so a few bisection steps are usually
needed first. The factor $C$ can be computed by
`arb_calc_newton_conv_factor`.

## Example {#example}

The program isolates the roots of $f(x) = \sin x + x/10$ on $[-20, 20]$
at 64 bits, then refines the largest root to 200 bits. The function has
five roots in this range. The root $x = 0$ is a bisection point, so it
is not isolated: it ends up in two tiny intervals with flag 2, after 30
bisection levels.

```c example name=arb_calc_isolate_roots
#include <flint/arb.h>
#include <flint/arb_poly.h>
#include <flint/arb_calc.h>

/* f(x) = sin(x) + x/10; out = Taylor coefficients of f at inp */
static int
f(arb_ptr out, const arb_t inp, void * param, slong order, slong prec)
{
    arb_struct x[2];

    arb_init(x + 0);
    arb_init(x + 1);
    arb_set(x + 0, inp);
    arb_one(x + 1);                           /* the series inp + t */
    _arb_poly_sin_series(out, x, FLINT_MIN(2, order), order, prec);
    arb_div_ui(x + 0, inp, 10, prec);
    arb_add(out + 0, out + 0, x + 0, prec);   /* + x/10 */
    if (order > 1)
    {
        arb_set_ui(x + 1, 10);
        arb_inv(x + 1, x + 1, prec);
        arb_add(out + 1, out + 1, x + 1, prec);
    }
    arb_clear(x + 0);
    arb_clear(x + 1);
    return 0;
}

int main(void)
{
    arf_interval_t interval, J;
    arf_interval_ptr blocks;
    int * flags;
    slong i, n;
    arb_t r, I;
    arf_t C;

    arf_interval_init(interval);
    arf_interval_init(J);
    arb_init(r);
    arb_init(I);
    arf_init(C);

    arf_set_si(&interval->a, -20);
    arf_set_si(&interval->b, 20);
    n = arb_calc_isolate_roots(&blocks, &flags, f, NULL, interval,
            30, 10000, WORD_MAX, 64);

    for (i = 0; i < n; i++)
    {
        flint_printf("flag %d: ", flags[i]);
        arf_interval_printd(blocks + i, 8);
        flint_printf("\n");
    }

    /* refine the last root: 20 bisections, then Newton to 200 bits */
    arf_interval_get_arb(I, blocks + n - 1, 64);
    arb_calc_newton_conv_factor(C, f, NULL, I, 64);
    arb_calc_refine_root_bisect(J, f, NULL, blocks + n - 1, 20, 64);
    arf_interval_get_arb(r, J, 64);
    if (arb_calc_refine_root_newton(r, f, NULL, r, I, C, 0, 200) == ARB_CALC_SUCCESS)
    {
        arb_printn(r, 50, 0);
        flint_printf("\n");
    }

    _arf_interval_vec_clear(blocks, n);
    flint_free(flags);
    arf_interval_clear(interval);
    arf_interval_clear(J);
    arb_clear(r);
    arb_clear(I);
    arf_clear(C);
    return 0;
}
```
```output
flag 1: [-7.5000000, -5.0000000]
flag 1: [-3.7500000, -2.5000000]
flag 2: [-3.7252903e-8, 0]
flag 2: [0, 3.7252903e-8]
flag 1: [2.5000000, 3.7500000]
flag 1: [5.0000000, 7.5000000]
[5.6792077963144036704201845877840800977761539706112 +/- 3.68e-50]
```

The 20 bisection steps reduce the interval to a width of about
$2.4 \cdot 10^{-6}$, from which Newton iteration reaches 200 bits in
a few steps.
