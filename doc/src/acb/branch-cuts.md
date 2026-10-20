---
title: Branch cuts
summary: The principal branches used by the multivalued complex functions, what happens when a rectangle touches or crosses a branch cut, the analytic variants that return indeterminate results there, and functions that avoid branch cuts altogether.
order: 2
---

The complex logarithm, square root, powers, inverse trigonometric and
inverse hyperbolic functions, the polylogarithm, the log-gamma function
and others are multivalued. The functions of `acb` return their principal
values, defined by branch cuts that follow the conventions of the DLMF
[@NIST2012] and of most computer algebra systems.

## The principal branches {#conventions}

| Function | Branch cut | Value on the cut |
|---|---|---|
| `acb_log`, `acb_arg` | $(-\infty, 0]$ | continuous from above: $\arg x = \pi$ for $x < 0$ |
| `acb_sqrt`, `acb_rsqrt`, `acb_root_ui`, `acb_pow` | $(-\infty, 0)$ | continuous from above: $\sqrt{-1} = i$ |
| `acb_asin`, `acb_acos` | $(-\infty, -1) \cup (1, +\infty)$ | |
| `acb_atan` | $(-i\infty, -i) \cup (i, +i\infty)$ | |
| `acb_asinh` | $(-i\infty, -i) \cup (i, +i\infty)$ | |
| `acb_acosh` | $(-\infty, 1)$ | |
| `acb_atanh` | $(-\infty, -1) \cup (1, +\infty)$ | |
| `acb_lgamma`, `acb_log_barnes_g` | $(-\infty, 0]$ | |
| `acb_polylog` | $(1, +\infty)$ | |
| `acb_agm1` | $(-\infty, 0)$ | continuous from above |
| `acb_lambertw` | see [@CGHJK1996] | |

The inverse functions are defined by formulas in terms of the logarithm
and the square root, so their branch cuts follow from those of `acb_log`
and `acb_sqrt`; the definitions are given with each function on the
module page.

## Rectangles on and across a cut {#crossing}

A function evaluated at a point exactly on a branch cut returns the value
given by the convention above. A rectangle that only touches or crosses
the cut contains points on both sides, where the function values differ
by a jump, so the result must enclose both values: the imaginary part of
`acb_log` of a rectangle that crosses the negative real axis encloses
$[-\pi, \pi]$, and the square root of a rectangle around $-4$ encloses
both $2i$ and $-2i$.

```c example name=acb_branch_cuts
#include <flint/flint.h>
#include <flint/acb.h>

#define SHOW(label, x) do { flint_printf("%-26s", label); \
    acb_printn(x, 8, 0); flint_printf("\n"); } while (0)

int main(void)
{
    acb_t z, w, w2;
    slong prec = 64;

    acb_init(z); acb_init(w); acb_init(w2);

    acb_set_si(z, -4);                                /* on the cut */
    acb_sqrt(w, z, prec);                SHOW("sqrt(-4)", w);
    acb_log(w, z, prec);                 SHOW("log(-4)", w);

    arb_set_str(acb_imagref(z), "-1e-10", prec);     /* just below it */
    acb_sqrt(w, z, prec);                SHOW("sqrt(-4 - 1e-10 i)", w);

    arb_set_str(acb_imagref(z), "[+/- 1e-10]", prec); /* across it */
    acb_sqrt(w, z, prec);                SHOW("sqrt(-4 + [+/- 1e-10] i)", w);
    acb_sqrt_analytic(w, z, 1, prec);    SHOW("analytic version", w);
    acb_sqrts(w, w2, z, prec);           SHOW("sqrts, first", w);
                                         SHOW("sqrts, second", w2);

    acb_clear(z); acb_clear(w); acb_clear(w2);
    flint_cleanup();
    return 0;
}
```
```output
sqrt(-4)                  2.0000000*I
log(-4)                   [1.3862944 +/- 3.89e-8] + [3.1415927 +/- 4.65e-8]*I
sqrt(-4 - 1e-10 i)        [2.5000000e-11 +/- 9.66e-31] + [-2.0000000 +/- 3.50e-21]*I
sqrt(-4 + [+/- 1e-10] i)  [+/- 2.51e-11] + [+/- 2.01]*I
analytic version          nan + nan*I
sqrts, first              [+/- 2.51e-11] + [2.0000000 +/- 1.57e-22]*I
sqrts, second             [+/- 2.51e-11] + [-2.0000000 +/- 1.57e-22]*I
```

This is a frequent source of unexpectedly wide results: a quantity that
is mathematically a negative real number but computed with a tiny
imaginary error term $[\pm \varepsilon] i$ straddles the cut. When the
imaginary part is known to be zero, set it to zero exactly
(`arb_zero(acb_imagref(z))`) before calling the function.

## Analytic variants {#analytic}

Rigorous integration (`acb_calc_integrate`) and other algorithms that rely
on analyticity must know when the integrand is not analytic on the input.
The functions `acb_sqrt_analytic`, `acb_rsqrt_analytic`,
`acb_log_analytic` and `acb_pow_analytic` take a flag `analytic`; when it
is set and the input rectangle touches the branch cut (including a point
exactly on it), they return an indeterminate result instead of an
enclosure of both branches. The piecewise real functions `acb_real_abs`,
`acb_real_sgn`, `acb_real_floor`, `acb_real_max` and others behave in the
same way at their discontinuities. Most special functions in
[acb_hypgeom](acb_hypgeom) and elsewhere do not detect branch cuts
automatically; their documentation states the cuts, and an integrand that
uses them near a cut must check this itself.

## Avoiding branch cuts {#avoiding}

- `acb_sqrts` returns both square roots, with no loss of accuracy near
  the cut, when the choice of branch does not matter or is made later.
- `acb_log_sin_pi` and `acb_lgamma` implement $\log \sin(\pi z)$ and
  $\log \Gamma(z)$ with the branch structure that makes functional
  equations hold everywhere, rather than as $\log$ of the function.
- `acb_agm1` and `acb_agm` choose the optimal square root at each step of
  the AGM iteration and are continuous except on a single cut.
- `acb_lambertw` provides nonstandard branches (`ACB_LAMBERTW_LEFT` and
  `ACB_LAMBERTW_MIDDLE`) that are continuous across the standard cuts.
